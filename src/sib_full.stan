data {
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> F;
    int<lower=max(M, F), upper=M*F> G;
    int<lower=1, upper=N> grp_size[G];

    int<lower=1, upper=2> method;

    real osr_mu;
    real<lower=0> osr_sigma;
}


transformed data {
    int<lower=0, upper=N> n[G+1];

    n[1:G] = grp_size;
    n[G+1] = 0;
    
}


parameters {
    real osr_logit;
    real<lower=F+M> N0;
    real<lower=M, upper=N0-F> M0;
    real<lower=G, upper=M0*(N0-M0)> G0;
}


transformed parameters {
    real<lower=0, upper=1> osr = inv_logit(osr_logit);
    real<lower=F, upper=N0-M> F0 = N0 - M0;
    simplex[G+1] theta;
    
    theta[1:G] = rep_vector(1/G0, G);
    theta[G+1] = 1 - G/G0;
}


model {
    osr_logit ~ normal(osr_mu, osr_sigma);

    if (method == 1) { 
	/* use normal approximation for binomial */    
	real mu = N0 * osr;
	real sigma = sqrt(mu * (1 - osr));

        M0 ~ normal(mu, sigma);
    }
    else {
        /* use relationship between binomial and beta distribution */
        target += -log(N0 + 1) + beta_lpdf(osr | M0 + 1, N0 - M0 + 1);
    }

    //target += lchoose(M0*F0 - G, G0 - G) - lchoose(M0*F0, G0);
    G ~ binomial(M*F, G0 / (M0 * F0));

    n ~ multinomial(theta);
    target += lgamma(G0 + 1) - lgamma(G0 - G + 1);
}
