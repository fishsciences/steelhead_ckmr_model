functions {
  real spawn_prob(real age, vector as, vector ps, real radius) {
    vector[size(as)] tmp;
    if (size(ps) != size(as)) {
      reject("log_spawn_prob: dimension mismatch: AS and PS must have the same dimension");
    }
    tmp = age - as; 
    tmp /= radius;
    tmp = -tmp .* tmp;
    tmp = softmax(tmp);
    return dot_product(ps, tmp);
  }

  real log_surv_prob(real t, real lambda, real kappa_0, real kappa_1) {
    real res = (t / lambda)^(kappa_1 * t + kappa_0);
    return -res;
  }

  real switch_fn(real x, real x0, real mu ) {
    real z;
    real f1;
    real f2;
    if (mu < 0) {
      reject("switch_fn: the scale parameter mu must be positive");
    }
    /* shift and scale */
    z = (x - x0) / mu;
    f1 = z > 0 ? exp(-1/z) : 0;
    f2 = z < 1 ? exp(1/(z - 1)) : 0;
    return f1 / (f1 + f2);
  }

  /* 
   * this function computes the probability that
   * a fish which hatched at time t_hatch and was 
   * last observed to be alive at time t_alive would
   * be alive at time t.
   *
   * Because this function would have discontinuous 
   * gradient, we used a mollified version that has 
   * smooth gradient and is equal to the exact function
   * outside a small neighborhood of a discontinuity.
   * The `mollify` parameter is the redius of this 
   * neighborhood
   */ 
  real alive_prob(real t, real t_hatch, real t_alive, real lambda, real kappa_0, real kappa_1, real mollify) {
    real srv, mu1, mu2;
    mu1 = switch_fn(t, t_hatch, mollify);
    mu2 = switch_fn(t,  t_alive, mollify);
    srv  = log_surv_prob(t - t_hatch, lambda, kappa_0, kappa_1);
    srv -= log_surv_prob(t_alive - t_hatch, lambda, kappa_0, kappa_1);
    srv  = exp(srv); 
    return mu1 * (1 - mu2) + mu2 * srv;
  }
  
}


data {
  /* number of (adult) age groups */
  int<lower=1> G;

  /* number of years */
  int<lower=1> Y;

  /* known fish */
  int<lower=1> N_known;                              // total number of known fish
  array[N_known] int<lower=1, upper=3> sex;          // fish sex assignments (1 = male, 2 = female, 3 = unknown)
  array[N_known] int<lower=1, upper=Y> first_year;   // first year the fish was observed
  array[N_known, 2] real<lower=0, upper=Y+1> t_span; // times of first and last observation of this fish;
  array[N_known] int<lower=1, upper=2> life_stage;   // lifes stage at first observation (1 = juvenile, 2 = adult)
  array[N_known] int<lower=1, upper=3> life_style;   // lifestyle 1 = resident, 2 = migratory, 3 = unknown 
  array[N_known] int<lower=0, upper=1> genotyped;    // was this fish genotyped
  matrix[N_known, 2] hatch;                          // hatch time low and high bounds 

  /* juvenile observations */
  int<lower=1> N_juv;                                // number of juveniles observations 
  array[N_juv] int<lower=1> juv_count;               // number of times a juvenile was observed
  array[N_juv] int<lower=1, upper=N_known> juv_fish; // which fish is this?
  array[N_juv] int<lower=1, upper=Y> juv_year;       // year of observation (data from one year is contiguous)
  array[N_juv] int<lower=1, upper=4> juv_ckmr;       // ckmr outcome: 1 = no parent, 2 = only mother, 3 = only father, 4 = both mother and father  

  /* adult observations */
  int<lower=1> N_adult;                                // number of adult observations 
  array[N_adult] int<lower=0> adult_count;             // number of times an adult was observed
  array[N_juv] int<lower=1, upper=N_known> adult_fish; // which fish is this?
  array[N_adult] int<lower=1, upper=Y> adult_year;     // year of observation (data from one year is contiguous)

  /* CKMR observations */
  array[Y, 4] int<lower=1, upper=N_known> N_ckmr_juv;    // number of juvenile CKMR probes in each year, with outcome:
                                                         //     1 = no parent found
							 //     2 = only mother
							 //     3 = only father
							 //     4 = both mother and father


  /* von Bertalanffy observations */
  int<lower=1> N_vB;                                 // number of size observations
  array[N_vB] int<lower=1, upper=N_known> vB_fish;   // observed fish (observation for same fish are consecutive) 
  vector<lower=0, upper=Y+1>[N_vB] vB_time;          // observation time
  vector<lower=0>[N_vB] vB_fl;                       // fork length 
  /** 
   ** priors for von Bertalanffy model indexed by:
   **   sex (1 = male, 2 = female)
   **   migratory status (1 = migratory, 2 = resident)
   **/
  array[2, 2] real<lower=0> vB_K_mu;
  array[2, 2] real<lower=0> vB_K_sigma;
  array[2, 2] real<lower=0> vB_L_inf_mu;
  array[2, 2] real<lower=0> vB_L_inf_sigma;
  array[2, 2] real<upper=0> vB_t0_mu;
  array[2, 2] real<lower=0> vB_t0_sigma;


  /* survival observations */
  int<lower=1> N_surv;                                 // number of survival observations 
  array[N_surv] int<lower=1, upper=N_known> surv_fish; // observed fish
  vector<lower=0, upper=Y>[N_surv] surv_time;          // last time observed to be alive 

  /* spawn event observations */
  int<lower=1> N_spawn;                                  // number of spawn event observations
  array[N_spawn] int<lower=1, upper=N_known> spawn_fish; // fish the observation refers to
  array[N_spawn, 2] real spawn_time;                     // lower and upper bounds on spawn event time
  array[N_spawn] int<lower=0, upper=1> spawn_abs;        // is the observation time absolute
  vector<lower=0, upper=1>[G] spawn_mu;                  // mean for spawn probability priors for age groups 
  vector<lower=2>[G] spawn_kappa;                        // evidence level for age group spawn probability prior
  real<lower=0> spawn_radius_epsilon;                    // prior parameter for the smoothing radius of the  spawn probability function       

  /* mollifier parameter */
  real<lower=0> mollify;
}

transformed data {
  /* constants */
  int MALE = 1;
  int FEMALE = 2;
  int MIGRATORY = 1;
  int RESIDENT = 2;
  int JUVENILE = 1;
  int ADULT = 2;
  int UNKNOWN = 3;
  int FIRST = 1;
  int LAST = 2;
  int LOW = 1;
  int HIGH = 2;

  /* number of genotyped fish */
  int<lower=1, upper=N_known> N_geno = sum(genotyped); 

  /* female sample distribution */
  array[Y] int N_sex_female;
  array[Y] int N_sex_total;

  /* migratory fish sample distribution */
  array[Y] int N_mig_mark;
  array[Y] int N_mig_total;

  /* number of distinct juveniles observed each year */
  array[Y] int<lower=1, upper=N_known> N_juv_known;
  array[Y, 2] int<lower=1, upper=N_juv> juv_span;  

  /* number of distinct adults observed each year */
  array[Y] int<lower=1, upper=N_known> N_adult_known;
  array[Y, 2] int<lower=1, upper=N_adult> adult_span;  

  vector<lower=0>[G+1] spawn_age_points; 
  
  /* sex observations */
  for (i in 1:N_known) {
    if (life_stage[i] == JUVENILE) continue;
    if (sex[i] != UNKNOWN) {
      N_sex_total[first_year[i]] += 1;
    }
    if (sex[i] == FEMALE) {
      N_sex_female[first_year[i]] += 1;
    }
  }

  /* migratory check observations */
  for (i in 1:N_known) {
    if (life_stage[i] == JUVENILE) continue;
    if (life_style[i] != UNKNOWN) {
      N_mig_total[first_year[i]] += 1;
    }
    if (life_style[i] == MIGRATORY) {
      N_mig_mark[first_year[i]] += 1;
    }
  }

  /* juvenile observations */
  {
    int count = 0;
    int year = 0;
    for (i in 1:N_juv) {
      if (juv_year[i] != year) {
	/* new year */
	if (year > 0) {
          juv_span[year, LAST] = i - 1;
	  N_juv_known[year] = count;
	}
	year = juv_year[i];
	juv_span[year, FIRST] = i;
      }
      count += 1;
    }
    juv_span[year, LAST] = N_juv;
  }

  /* adult observations */
  {
    int count = 0;
    int year = 0;
    for (i in 1:N_adult) {
      if (adult_year[i] != year) {
	/* new year */
	if (year > 0) {
          adult_span[year, LAST] = i - 1;
	  N_adult_known[year] = count;
	}
	year = adult_year[i];
	adult_span[year, FIRST] = i;
      }
      count += 1;
    }
    adult_span[year, LAST] = N_juv;
  }

  /* spawn probability hint locations*/ 
  for (i in 1:(G+1)) {
    spawn_age_points[i] = i;
  }
}

parameters {
  /* female ratio each year */
  vector<lower=0, upper=1>[Y] female_ratio;

  /* migratory fish ratio */
  real<lower=0, upper=1> mig_ratio;

  /* unknown fish in each year */
  array[Y] real<lower=0> N_juv_unknown;
  array[Y] real<lower=0> N_adult_unknown;

  /* number of active CKMR markers by sex (1 = male, 2 = female) */
  array[Y, 2] real<lower=0> N_ckmr_active;

  /* number of active CKMR duds by sex (1 = male, 2 = female) */
  array[Y, 2] real<lower=0> N_ckmr_dud;

  /* hatch time parameter for each known fish */
  vector<lower=0, upper=1>[N_known] hatch_theta;

  /* survival parameters, by sex (1 = male, 2 = female) 
   * and migratory status (1 = migratory, 2 = resident) 
   */
  array[2, 2] real<lower=0> surv_lambda;
  array[2, 2] real<lower=0> surv_kappa_1;
  array[2, 2] real<lower=0> surv_kappa_0;

  /* von Bertalanffy parameters (1 = male, 2 = female) */
  array[2, 2] real<lower=0> vB_L_inf;
  array[2, 2] real<lower=0> vB_K;
  array[2, 2] real<lower=0, upper=1> vB_theta;
  array[2, 2] real<lower=0> vB_sigma;

  /* spawn parameters */
  array[2, 2] vector<lower=0, upper=1>[G] spawn_p;  // spawn probability hints by sex and migratory status
  array[2, 2] real<lower=0, upper=1> spawn_radius;  // ditto for spawn sprobability smoothing radius
  vector<lower=0, upper=1>[N_spawn] spawn_theta;    // spawn event temporal location within bounds
}


transformed parameters {

  /* hatch times for known fish */
  vector[N_known] hatch_time = hatch[:, LOW] + (hatch[:, HIGH] - hatch[:, LOW]) .* hatch_theta;

  /* log probabilities of survival observations */
  real surv_lp;

  /* log probabilities of spawn event observations */
  real spawn_lp;

  /* log probability for vonBeratlanffy observations */
  real vB_lp;

  /* probabilistic sex assignments */
  array[N_known, 2] real sex_prob;

  /* migratory assesment */
  array[N_known, 2] real mig_prob; 

  /* overall fish type */
  array[N_known, 2, 2] real<lower=0, upper=1> fish_type_prob;

  /* number of active CKMR markers by year and sex (1 = male, 2 = female)*/
  array[Y, 2] real ckmr_active_mu;
  array[Y, 2] real ckmr_active_sigma2;

  /* number of active CKMR duds by year and sex (1 = male, 2 = female)*/
  array[Y, 2] real ckmr_dud_mu;
  array[Y, 2] real ckmr_dud_sigma2;

  /* 
   * fish type calculations 
   */

  /** sex assignment **/
  for (i in 1:N_known) {
    if (sex[i] == MALE) {
      sex_prob[i, MALE] = 1.0;
      sex_prob[i, FEMALE] = 0.0;
    }
    else if (sex[i] == FEMALE) {
      sex_prob[i, FEMALE] = 1.0;
      sex_prob[i, MALE] = 0.0;
    }
    else {
      sex_prob[i, MALE] = 1 - female_ratio[first_year[i]];
      sex_prob[i, FEMALE] = female_ratio[first_year[i]];
    }
  }

  /** migratory assesment **/
  for (i in 1:N_known) {
    if (life_style[i] == MIGRATORY) {
      mig_prob[i, MIGRATORY] = 1.0;
      mig_prob[i, RESIDENT] = 0.0;
    }
    else if (life_style[i] == RESIDENT) {
      mig_prob[i, RESIDENT] = 1.0;
      mig_prob[i, MIGRATORY] = 0.0;
    }
    else {
      mig_prob[i, RESIDENT] = 1 - mig_ratio;
      mig_prob[i, MIGRATORY] = mig_ratio;
    }
  }

  /** fish type **/
  for (i in 1:N_known) {
    for (j in 1:2) {
      for (k in 1:2) {
	fish_type_prob[i, j, k] = sex_prob[i, j] * mig_prob[i, j];
      }
    }
  }


  /* von Bertalanffy calculations */
  {
    vector[4] theta = to_vector(to_array_1d(vB_theta));
    vector[4] L_inf = to_vector(to_array_1d(vB_L_inf));
    vector[4] K = to_vector(to_array_1d(vB_K));
    vector[4] sigma0 = to_vector(to_array_1d(vB_sigma));
    vector[4] mix;
    vector[N_vB] lps;
    int fish = -1; // out of band value
    for (i in 1:N_vB) {
      real time;
      vector[4] x;
      vector[4] mu;
      vector[4] sigma;
      vector[4] jacobi;
      vector[4] lp;
      if (vB_fish[i] != fish) {
	/* first observation for this fish */
	time = vB_time[i] - hatch_time[i];
        x = (L_inf - vB_fl[i]) ./ (L_inf .* (1 - theta));
	jacobi = log(L_inf - vB_fl[i]) - 2 * (log(L_inf) - log(1 - theta));
	fish = vB_fish[i];
      }
      else {
	time = vB_time[i] - vB_time[i-1]; 
	x = (L_inf - vB_fl[i]) ./ (L_inf - vB_fl[i - 1]);
	jacobi = log(vB_fl[i] - vB_fl[i-1]) - 2 * log(L_inf - vB_fl[i-1]);
      }
      mu = -time * K;
      sigma = sqrt(time) * sigma0;
      mix = log(to_vector(to_array_1d(fish_type_prob[i])));
      mix += lognormal_lpdf(x | mu, sigma);
      mix -= jacobi;
      lps[i] = log_sum_exp(mix);
    }
    vB_lp = sum(lps);
  }

  /* survival calculations */
  {
    vector[N_surv] tmp;
    for (i in 1:N_surv) {
      real age = surv_time[i] - hatch_time[surv_fish[i]];
      real lambda;
      real kappa_0;
      real kappa_1;
      array[2, 2] real lsp;
      for (j in 1:2) {
        for (k in 1:2) {
          lambda = surv_lambda[j, k];
          kappa_0 = surv_kappa_0[j, k];
          kappa_1 = surv_kappa_1[j, k];
	  lsp[j, k] = log(fish_type_prob[i, j, k]);
          lsp[j, k] += log_surv_prob(age, lambda, kappa_0, kappa_1);
	}
      }
      tmp[i] = log_sum_exp(to_vector(to_array_1d((lsp))));
    }
    surv_lp = sum(tmp);
  }

  /* spawn model calculations */
  { 
    vector[N_spawn] tmp;
    for (i in 1:N_spawn) {
      real t = spawn_time[i, 1] + spawn_theta[spawn_fish[i]] .* spawn_time[i, 2];
      real age = spawn_abs[i] ? t - hatch_time[spawn_fish[i]] : t;
      matrix[2, 2] sp;
      vector[4] sp_vec;
      vector[4] ft_vec;
      for (j in 1:2) {
        for (k in 1:2) {
          vector[G+1] ps;
          ps[1] = 0;
          ps[2:] = spawn_p[j, k];
          sp[j, k] = spawn_prob(age, spawn_age_points, ps, spawn_radius[j, k]);
        }
      }
      sp_vec = to_vector(to_array_1d(sp));
      ft_vec = to_vector(to_array_1d(fish_type_prob[i, :, :]));
      tmp[i] = dot_product(sp_vec, ft_vec);
    }
    spawn_lp = sum(log(tmp));
  }

  /* CKMR calculations */
  for (i in 1:N_known) {
    /* CKMR marker */      
    for (year in 1:Y) {
      real t = Y - 0.5; 
      array[2, 2] real p_active;
      for (i_sex in 1:2) {
        for (i_mig in 1:2) {
          real lambda = surv_lambda[i_sex, i_mig];
          real kappa_0 = surv_kappa_0[i_sex, i_mig];
          real kappa_1 = surv_kappa_1[i_sex, i_mig];
	  real radius = spawn_radius[i_sex, i_mig];
          vector[G+1] ps;
          ps[1] = 0;
          ps[2:] = spawn_p[i_sex, i_mig];
          p_active[i_sex, i_mig] = fish_type_prob[i, i_sex, i_mig];
          p_active[i_sex, i_mig] *= alive_prob(t, hatch_time[i], t_span[i, LAST], lambda, kappa_0, kappa_1, mollify);
          p_active[i_sex, i_mig] *= spawn_prob(t - hatch_time[i], spawn_age_points, ps, radius); 
        }
      }
      for (i_sex in 1:2) {
        if (genotyped[i]) {
          /* CKMR marker */
          real p = p_active[i_sex, MIGRATORY] + p_active[i_sex, RESIDENT];
          ckmr_active_mu[year, i_sex] += p;
          ckmr_active_sigma2[year, i_sex] += p * (1 - p);
        }
        else {
	  /* CKMR dud */
          real p = p_active[i_sex, MIGRATORY] + p_active[i_sex, RESIDENT];
          ckmr_dud_mu[year, i_sex] += p;
          ckmr_dud_sigma2[year, i_sex] += p * (1 - p);
        }
      }
    }
  }
}


model {   
  /* survival */
  to_vector(to_array_1d(surv_lambda)) ~ cauchy(0, 10);
  to_vector(to_array_1d(surv_kappa_0)) ~ cauchy(1, 1);
  to_vector(to_array_1d(surv_kappa_1)) ~ exponential(1);
  target += surv_lp;


  /* spawn events */
  for (i_sex in 1:2) {
    for (i_mig in 1:2) {
      spawn_radius[i_sex, i_mig] ~ beta_proportion(0.5, 2 + spawn_radius_epsilon);
      spawn_p[i_sex, i_mig] ~ beta_proportion(spawn_mu, spawn_kappa);
    }
  }
  target += spawn_lp;


  /* von Bertalanffy growth model */
  {
    vector[4] mu = to_vector(to_array_1d(vB_K_mu));
    vector[4] sigma = to_vector(to_array_1d(vB_K_sigma));
    to_vector(to_array_1d(vB_K)) ~ normal(mu, sigma);
  }
  {
    vector[4] mu = to_vector(to_array_1d(vB_L_inf_mu));
    vector[4] sigma = to_vector(to_array_1d(vB_L_inf_sigma));
    to_vector(to_array_1d(vB_L_inf)) ~ normal(mu, sigma);
  }
  {
    
    vector[4] mu = to_vector(to_array_1d(vB_t0_mu));
    vector[4] sigma = to_vector(to_array_1d(vB_t0_sigma));
    vector[4] theta = to_vector(to_array_1d(vB_theta));
    vector[4] K = to_vector(to_array_1d(vB_K));
    vector[4] t0 = log(1 - theta) ./ K;
    t0 ~ normal(mu, sigma);
    target += -sum(log(K) + log(1 - theta)); // because we sample t0, a transformed variable 
  }
  target += vB_lp;

  /* juvenile within year recapture observations */
  for (y in 1:Y) {
    real N = N_juv_known[y] + N_juv_unknown[y];
    real theta_1 = 1.0 / N;
    real theta_2 = N_juv_unknown[y] / N;
    vector[N_juv_known[y] + 1] theta;
    array[N_juv_known[y] + 1] int draw;

    theta[1:N_juv_known[y]] = rep_vector(theta_1, N_juv_known[y]);
    theta[N_juv_known[y]+1] = theta_2;
    draw[1:N_juv_known[y]] = juv_count[juv_span[y, FIRST] : juv_span[y, LAST]];
    draw[N_juv_known[y]+1] = 0;

    draw ~ multinomial(theta);
  }

  /* adult within year recapture observations */
  for (y in 1:Y) {
    real N = N_adult_known[y] + N_adult_unknown[y];
    real theta_1 = 1.0 / N;
    real theta_2 = N_adult_unknown[y] / N;
    vector[N_adult_known[y] + 1] theta;
    array[N_adult_known[y] + 1] int draw;

    theta[1:N_adult_known[y]] = rep_vector(theta_1, N_adult_known[y]);
    theta[N_adult_known[y]+1] = theta_2;
    draw[1:N_adult_known[y]] = adult_count[adult_span[y, FIRST] : adult_span[y, LAST]];
    draw[N_adult_known[y]+1] = 0;

    draw ~ multinomial(theta);
  }

  /* sex ratio */
  N_sex_female ~ binomial(N_sex_total, female_ratio);

  /* migratory ratio */
  N_mig_mark ~ binomial(N_mig_total, mig_ratio);

  /* CKMR model */
  {
    vector[2*Y] mu = to_vector(to_array_1d(ckmr_active_mu));
    vector[2*Y] sigma = sqrt(to_vector(to_array_1d(ckmr_active_sigma2)));

    to_vector(to_array_1d(N_ckmr_active)) ~ normal(mu, sigma);
  }
  {
    vector[2*Y] mu = to_vector(to_array_1d(ckmr_dud_mu));
    vector[2*Y] sigma = sqrt(to_vector(to_array_1d(ckmr_dud_sigma2)));

    to_vector(to_array_1d(N_ckmr_dud)) ~ normal(mu, sigma);
  }
  for (y in 1:Y) {
    real F0 = N_adult_unknown[y] * female_ratio[y] + N_ckmr_dud[y, FEMALE];
    real M0 = N_adult_unknown[y] * (1 - female_ratio[y]) + N_ckmr_dud[y, MALE];
    real F = N_ckmr_active[y, FEMALE];
    real M = N_ckmr_active[y, MALE];
    vector[4] theta;

    theta[1] = log(M0) - log(M + M0) + log(F0) - log(F + F0);   // neither parent observed
    theta[2] = log(F) - log(F + F0) + log(M0) - log(M + M0);    // only mother observed
    theta[3] = log(M) - log(M + M0) + log(F0) - log(F + F0);    // only father observed
    theta[4] = log(M) + log(F) - log(M + M0) - log(F + F0);     // both parrents observed

    N_ckmr_juv[y, :] ~ multinomial(softmax(theta));
  }
}
