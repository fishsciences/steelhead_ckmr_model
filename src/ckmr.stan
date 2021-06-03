data {
    /*
     * Observations are carried out in T consecutive
     * years
     */
    int<lower=1> T;

    /*
     * number of fish observed in each state, in each 
     * observation year
     */
    int<lower=0> observed_count[T, 5];

    /*
     * number of recaptures observed in each year, in
     * each state
     */
    int<lower=0> recaptured_count[T, 5];

    /*
     * number of offsprings observed each year, in
     * each state
     */
    int<lower=0> offspring_count[T, 5]; 

}


transformed data {
    /*
     * states of the Markov chain that models
     * the fish life cycle
     */
    int S = 5;
    int<lower=1, upper=S> STATE_JR = 1; // juvenile fish in river
    int<lower=1, upper=S> STATE_JO = 2; // juvenile fish in ocean
    int<lower=1, upper=S> STATE_AO = 3; // adult fish in ocean
    int<lower=1, upper=S> STATE_AR = 4; // adult fish in river
    int<lower=1, upper=S> STATE_D  = 5; // dead fish
    int<lower=1, upper=S> next_state[S - 1] = {STATE_JO, STATE_AO, STATE_AR, STATE_AO};

    /*
     * the three populations of fish we track
     */
    int P = 3;
    int<lower=1, upper=P> POP_OBS;  // directly observed fish
    int<lower=1, upper=P> POP_MARK; // offsprings of directly observed fish
    int<lower=1, upper=P> POP_UNKN; // every other fish

    /*
     * (uniform) initial distribution
     */
    simplex[S] init_dist = [0.25, 0.25, 0.25, 0.25, 0.0]'; 
     
    /*
     * pure distributions
     */
    simplex[S] pure_dist[S] = {
	    [1.0, 0.0, 0.0, 0.0, 0.0]' ,
	    [0.0, 1.0, 0.0, 0.0, 0.0]' ,
	    [0.0, 0.0, 1.0, 0.0, 0.0]' ,
	    [0.0, 0.0, 0.0, 1.0, 0.0]' ,
	    [0.0, 0.0, 0.0, 0.0, 1.0]' 
    };

    /*
     * newly observed fish quantity and distribution
     */
    vector[T] new_obs_pop;
    simplex[S] new_obs_dist[T];

    /*
     * recaptured fish quantity and distribution
     */
    vector[T] reobs_pop;
    simplex[S] reobs_dist[T];

    
    init_dist = rep_vector(1.0, S);
    init_dist = init_dist / S;

    for (t in 1:T) {
        vector[S] tmp = to_vector(observed_count[t, :]) - to_vector(recaptured_count[t, :]);
        new_obs_pop[t] = sum(tmp);
	new_obs_dist[t] = tmp / new_obs_pop[t];
	reobs_pop[t] = sum(to_vector(recaptured_count[t, :]));
	reobs_dist[t] = to_vector(recaptured_count[t, :]) / reobs_pop[t];
    }

     
}


parameters {
    /* 
     * fish mortality in each of the alive states
     */
    vector<lower=0, upper=1>[S - 1] mortality;

    /*
     * probability that a fish remains in its current
     * state for another time step
     */
     vector<lower=0, upper=1>[S - 1] prob_stasis;

    /*
     * the initial size of each tracked population
     */
    real<lower=0> init_pop;


    /*
     * the probability that an adult fish survives
     * the migration long enough to reproduce
     */
    real<lower=0, upper=1> prob_rep;

    /* 
     * The distribution of the number number of offsprings 
     * an individual will spawn in a particular year, provided
     * the individual survives long enough to reproduce 
     */
    real<lower=0> mu_fert;
    real<lower=0> sigma_fert;

    /* reproductive population sizes */
    real marked_rep_pop[T];
    real unknown_rep_pop[T];


    /* spawn size (adds to population each year) */ 
    real marked_spawn_size[T];
    real unknown_spawn_size[T];
}


transformed parameters {
    matrix<lower=0, upper=1>[S, S] markov_trans;
    real<lower=0> pop_size[T, P];
    simplex[S] pop_dist_in[T, P];
    simplex[S] pop_dist_out[T, P];
    real<lower=0, upper=1> recapture_prob[T, S];
    real<lower=0, upper=1> offspring_prob[T, S];

    /*
     * populate the non-zero entries of the Markov 
     * transition matrix with the right values
     */
    markov_trans = rep_matrix(0.0, S, S); 
    for (s in 1:S) {
        if (s == STATE_D) {
	    /* no undead fish, please */
	    markov_trans[s, s] = 1.0;
	}
	else {
            markov_trans[s, STATE_D] = mortality[s];
	    markov_trans[s, s] = prob_stasis[s];
	    markov_trans[s, next_state[s]] = 1.0 - mortality[s] - prob_stasis[s];
	}
    }



    /*
     * time evolution for populations and their 
     * distribution over states
     */
    for (t in 1:T) {

	/* 
	 * Markov evolution for one year to the next
	 */
	for (p in 1:P) {
	    pop_dist_in[t, p] =  markov_trans * (t > 1 ? pop_dist_out[t - 1, p] : init_dist);
	}

        /*
	 * intra-year evolution for each population
	 */
        {
	    real beta;
            vector[S] alpha;
	    real old_pop = t > 1 ? pop_size[t - 1, POP_OBS] : 0;
	    real new_pop = 0;

	    for (s in 1:S) {
	        alpha[s] = observed_count[t, s];
		old_pop -= recaptured_count[t, s];
	        new_pop += observed_count[t, s];
	    }

	    pop_size[t, POP_OBS] = old_pop + new_pop;
	    beta = old_pop / pop_size[t, POP_OBS];
	    alpha /= new_pop;
	    pop_dist_out[t, POP_OBS] = beta * pop_dist_in[t, POP_OBS] + (1 - beta) * alpha;
        }

	{
	    real beta;
	    real old_pop = t > 1 ? pop_size[t - 1, POP_MARK] : 0;
	    real new_pop = marked_spawn_size[t];
	    vector[S] alpha = old_pop * pop_dist_in[t, POP_MARK];

	    for (s in 1:S) {
		old_pop -= offspring_count[t, s];
		alpha[s] -= offspring_count[t, s];
            }

	    pop_size[t, POP_MARK] = old_pop + new_pop;
	    alpha /= t > 1 ? old_pop : 1.0;
	    beta = old_pop / pop_size[t, POP_MARK];
	    pop_dist_out[t, POP_MARK] = beta * alpha + (1 - beta) * pure_dist[STATE_JR];
	}

        {
            real beta;
	    real old_pop = t > 1 ? pop_size[t - 1, POP_UNKN] : init_pop;
	    real new_pop = unknown_spawn_size[t];
	    vector[S] alpha = old_pop * pop_dist_in[t, POP_UNKN];

	    for (s in 1:S) {
                real delta = observed_count[t, s] - recaptured_count[t, s] - offspring_count[t, s];
		old_pop -= delta;
		alpha[s] -= delta;
            }

	    pop_size[t, POP_UNKN] = old_pop + new_pop;
	    alpha /= old_pop;
	    beta = old_pop / pop_size[t, POP_UNKN];
	    pop_dist_out[t, POP_UNKN] = beta * pop_dist_in[t, POP_UNKN] + (1 - beta) * pure_dist[STATE_JR];
        }

        /* 
	 * probabilities of various kinds of observations 
	 */
	{
            real total_pop = pop_size[t, POP_OBS] + pop_size[t, POP_MARK] + pop_size[t, POP_UNKN];
	    real prob_pop[P];

	    for (p in 1:P) {
	        prob_pop[p] = pop_size[t, p] / total_pop;
	    }

	    for (s in 1:S) {
		real denom = 0;

		for (p in 1:P) {
		    denom += prob_pop[p] * pop_dist_in[t, p][s];
		}
		
                recapture_prob[t, s] = prob_pop[POP_OBS] * pop_dist_in[t, POP_OBS][s] / denom;
	        offspring_prob[t, s] = prob_pop[POP_MARK] * pop_dist_in[t, POP_MARK][s] / denom;
            }
        }
    }
}


model {
    /* FIXME: add hyperpriors */

    for (t in 1:T) {
	{
            real pop = t > 1 ? pop_size[t - 1, POP_OBS] : 0;
            real mu = prob_rep * pop;	    
            real sigma = sqrt(prob_rep * (1 - prob_rep) * pop);
            marked_rep_pop[t] ~ normal(mu, sigma);
        }
        {
            real pop = t > 1 ? pop_size[t - 1, POP_MARK] + pop_size[t - 1, POP_UNKN] : init_pop;
	    real mu = prob_rep * pop;
	    real sigma = sqrt(prob_rep * (1 - prob_rep) * pop);
	    unknown_rep_pop[t] ~ normal(mu, sigma);
        }
	{
            real mu = marked_rep_pop[t] * mu_fert;
            real sigma = sqrt(marked_rep_pop[t] * sigma_fert);
	    marked_spawn_size[t] ~ normal(mu, sigma);	
	}
        {
            real mu = unknown_rep_pop[t] * mu_fert;
            real sigma = sqrt(unknown_rep_pop[t] * sigma_fert);
	    unknown_spawn_size[t] ~ normal(mu, sigma);	
        }
        for (s in 1:S) {
	    {
	        int n = observed_count[t, s];
                recaptured_count[t, s] ~ binomial(n, recapture_prob[t, s]);
            }
	    {
	        int n = observed_count[t, s] - recaptured_count[t, s];
                offspring_count[t, s] ~ binomial(n, offspring_prob[t, s]);
	    }
	}
        
    }
}


generated quantities {
  /* FIXME: what are we really trying to compute? */
}
