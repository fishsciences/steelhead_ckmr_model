data {
    int<lower=0> M; // number of adult males collected
    int<lower=0> F; // number of adult females collected
    array[4] int<lower=0> by_parents;
    real<lower=0, upper=1> mu;
    real<lower=0> kappa;
    /*
    int<lower=M+F> A0_gt;
    int<lower=F> F0_gt;
    */
}

transformed data {
  int<lower=0> A = M + F; //number of adults collected
  array[2] int<lower=0> fm = {F, M};
  real<lower=0> alpha = mu * kappa;
  real<lower=0> beta = kappa - alpha;
}

parameters { 
  real<lower=A> A0;
  real<lower=F/A0, upper=1-M/A0> phi; // female ratio of adult population
}

transformed parameters {
  real F0 = phi * A0;
  real M0 = A0 - F0;
  vector[4] theta;

  theta[1] = log(M0 - M) - log(M0) + log(F0 - F) - log(F0);   // neither parent observed
  theta[2] = log(F) - log(F0) + log(M0 - M) - log(M0);        // only mother observed
  theta[3] = log(M) - log(M0) + log(F0 - F) - log(F0);        // only father observed
  theta[4] = log(M) + log(F) - log(M0) - log(F0);             // both parrents observed

}

model {
  phi ~ beta(alpha,beta) T [F/A0, 1-M/A0];
  //omega ~ lognormal(mu, sigma);
  
  by_parents ~ multinomial(softmax(theta));
}

/*
generated quantities {
  real A0_err = (A0 - A0_gt) / A0_gt;
  real F0_err = (F0 - F0_gt) / F0_gt;
  real A0_ppq = (A0 < A0_gt);
  real F0_ppq = (F0 < F0_gt);
}
*/
