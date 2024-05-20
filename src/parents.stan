data {
    int<lower=0> M; // number of adult males collected
    int<lower=0> F; // number of adult females collected
    array[4] int<lower=0> by_parents;
    real<lower=0, upper=1> mu;
    real<lower=0> kappa;
}

transformed data {
  int<lower=0> A = M + F; //number of adults collected
  real<lower=0> alpha = mu * kappa;
  real<lower=0> beta = kappa - alpha;
}

parameters { 
  real<lower=1> F_unk;
  real<lower=1> M_unk;
}

transformed parameters {
  real F0 = F + F_unk;
  real M0 = M + M_unk;
  vector[4] theta;

  theta[1] = log(M_unk) - log(M0) + log(F_unk) - log(F0);   // neither parent observed
  theta[2] = log(F) - log(F0) + log(M_unk) - log(M0);       // only mother observed
  theta[3] = log(M) - log(M0) + log(F_unk) - log(F0);       // only father observed
  theta[4] = log(M) + log(F) - log(M0) - log(F0);           // both parrents observed

}

model {
  by_parents ~ multinomial(softmax(theta));
}

