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
    real<lower=0> osr;
    real<lower=fmax(F, M/osr)> F0;
    real<lower=G, upper=osr*F0*F0> G0;
}


transformed parameters {
    real<lower=M> M0 = osr * F0;
    simplex[G+1] theta;
    
    theta[1:G] = rep_vector(1/G0, G);
    theta[G+1] = 1 - G/G0;
}


model {
    osr ~ lognormal(osr_mu, osr_sigma);

    target += lchoose(M0, M) + lchoose(F0, F) - lchoose(M0*F0, G);

    n ~ multinomial(theta);
    target += lgamma(G0 + 1) - lgamma(G0 - G + 1);
}

generated quantities {
    real<lower=M+F> N0 = F0 + M0;
}
