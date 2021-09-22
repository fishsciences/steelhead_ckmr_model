data {
    int<lower=1> N;
    int<lower=1, upper=N> F;
    int<lower=1, upper=N> size[F];
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
    real<lower=F> F0;
}

transformed parameters {
    vector<lower=0, upper=1>[F] theta;
    
    
    {
        vector[F] tmp = cumulative_sum(rep_vector(1.0, F)) - 1.0;
        tmp = F0 - tmp;
        theta = 1.0 ./ tmp;
    }
}

model {
    k ~ binomial(n, theta);
}
