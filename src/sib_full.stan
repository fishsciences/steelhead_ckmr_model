data {
    int<lower=1> N;
    int<lower=1, upper=N> F;
    int<lower=1, upper=N> size[F];

    int<lower=1, upper=2> method;

    real osr_mu;
    real<lower=0> osr_sigma;
    real prom_mu;
    real<lower=0> prom_sigma;

}


transformed data {
    int<lower=0, upper=(N-1)> k[F];
    int<lower=0, upper=(N-1)> n[F];
    
    {
        int total = N;
        for (i in 1:F) {
            k[i] = size[i] - 1;
            n[i] = total - 1;
            total -= size[i];
        }
    }
}


parameters {
    real osr_logit;
    real prom_logit;
    real<lower=2*sqrt(F)> N0;
    real<lower=1, upper=N0-1> M;
    real<lower=F, upper=N0*N0/4> F0;
}


transformed parameters {
    real<lower=0, upper=1> osr = inv_logit(osr_logit);
    real<lower=0, upper=1> prom = inv_logit(prom_logit);
    real<lower=F> F_max = M * (N0 - M);
    vector<lower=0, upper=1>[F] theta;
    
    
    {
        vector[F] tmp = cumulative_sum(rep_vector(1.0, F)) - 1.0;
        tmp = F0 - tmp;
        theta = 1.0 ./ tmp;
    }
}


model {
    osr_logit ~ normal(osr_mu, osr_sigma);
    prom_logit ~ normal(prom_mu, prom_sigma);

    if (method == 1) { 
	/* use normal approximation for binomial */    
        M ~ normal(N0 * osr, sqrt(N0 * osr *(1 - osr)));
        F0 ~ normal(F_max * prom, sqrt(F_max * prom * (1 - prom)));
    }
    else {
        /* use relationship between binomial and beta distribution */
        target += -log(N0 + 1) + beta_lpdf(osr | M + 1, N0 - M + 1);
        target += -log(F_max + 1) + beta_lpdf(prom | F0 + 1, F_max - F0 + 1);
    }

    k ~ binomial(n, theta);
}
