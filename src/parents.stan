functions {
    real wallenius(data array[] int fm, array[] real FM, real omg);

    real wallenius(data array[] int fm, array[] real FM, real omg) {
	if (fm[1] == 0 && fm[2] == 0) {
	    return 0;
	}
	else if (fm[1] == 0 && fm[2] > 0) {
	    real w = wallenius({0, fm[2] - 1}, {FM[1], FM[2] - 1}, omg);
	    real p = log(FM[2]) - log(omg*FM[1] + FM[2]);
	    return p + w;
	}
	else if (fm[2] == 0 && fm[1] > 0) {
	    real w = wallenius({fm[1] - 1, 0}, {FM[1] - 1, FM[2]}, omg);
	    real p = log(omg) + log(FM[1]) - log(omg*FM[1] + FM[2]);
	    return p + w;
	}
	else {
	    row_vector[2] w = [
	        wallenius({fm[1] - 1, fm[2]}, {FM[1] - 1, FM[2]}, omg),
		wallenius({fm[1], fm[2] - 1}, {FM[1], FM[2] - 1}, omg)
	    ];
	    row_vector[2] p = [
	        log(omg) + log(FM[1]) - log(omg*FM[1] + FM[2]),
	        log(FM[2]) - log(omg*FM[1] + FM[2])
	    ];
	    return log_sum_exp(p + w);
	}
    }
}

data {
    int<lower=0> M; // number of adult males collected
    int<lower=0> F; // number of adult females collected
    array[4] int<lower=0> by_parents;
    real<lower=0> alpha;
    real<lower=0> beta;
    real mu;
    real<lower=0> sigma;
}

transformed data {
  int<lower=0> A = M + F; //number of adults collected
  array[2] int<lower=0> fm = {F, M};
}

parameters { 
  real<lower=A> A0;
  real<lower=F/A0, upper=1-M/A0> phi; // female ratio of adult population
  //real<lower=0> omega;                // Wallenius bias for sampling
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

  /* 
   * observed M, F follow a Wallenius non-central 
   * hypergeometric distribution
   */
  //target += wallenius(fm, {F0, M0}, omega);
  
  by_parents ~ multinomial(softmax(theta));
}
