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
    int<lower=0, upper=(N-1)> k[G];
    int<lower=0, upper=(N-1)> n[G];
    
    {
        int total = N;
        for (i in 1:G) {
            k[i] = grp_size[i] - 1;
            n[i] = total - 1;
            total -= grp_size[i];
        }
    }
}


parameters {
    real osr_logit;
    real<lower=F+M> N0;
    real<lower=M, upper=N0-F> M0;
    real<lower=G, upper=M0*(N0-M0)> G0;
}


transformed parameters {
    real<lower=0, upper=1> osr = inv_logit(osr_logit);
    real<lower=F, upper=N0> F0 = N0 - M0;
    real<lower=0, upper=1> phi = G0 / (M0 * F0);
    vector<lower=0, upper=1>[G] theta;
    
    
    {
        vector[G] tmp = cumulative_sum(rep_vector(1.0, G)) - 1.0;
        tmp = G0 - tmp;
        theta = 1.0 ./ tmp;
    }
}


model {
    osr_logit ~ normal(osr_mu, osr_sigma);

    if (method == 1) { 
	/* use normal approximation for binomial */    
	real mu = N0 * osr;
	real sigma = sqrt(mu * (1 - osr));

        M0 ~ normal(N0 * osr, sqrt(N0 * osr *(1 - osr)));
    }
    else {
        /* use relationship between binomial and beta distribution */
        target += -log(N0 + 1) + beta_lpdf(osr | M0 + 1, N0 - M0 + 1);
    }

    if (method == 1) {
	real mu = phi * M * F;
	real sigma = sqrt(mu * (1 - phi));
	G ~ normal(mu, sigma);
    }
    else {
        target += -log(M) -log(F) + beta_lpdf(phi | G + 1, M*F - G + 1);
    }

    k ~ binomial(n, theta);
}
