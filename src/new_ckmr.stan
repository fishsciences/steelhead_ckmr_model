data {
  /* number of years in study */
  int<lower=2> Y;

  /* number of obsevation methods */
  int<lower=2> M;

  /* number of age groups */
  int<lower=2> A;

  /* number of observed fish */
  int<lower=2> F;

  /* number of observations */
  int<lower=F> Obs;

  /* number of parent-offspring pairs */
  int<lower=1> P; 

  /* sex data */
  array[F] int<lower=1, upper=2> fish_sex;

  /* observation data */
  array[Obs] int<lower=1, upper=M> obs_method;
  array[Obs] int<lower=1, upper=F> obs_fish;
  array[Obs] int<lower=1, upper=Y> obs_year;
  array[Obs] int<lower=0, upper=A> obs_age;    // zero means no age estimate

  /* parentage data */
  array[P] int<lower=1, upper=F> parent;
  array[P] int<lower=1, upper=F> offspring;

  /* applicability map */
  array[M, 2, A] int<lower=0, upper=1> app;
}

transformed data {
  /* Sex indices */
  int MALE = 1;
  int FEM  = 2;


  /* age of each fish in each study year */
  array[F, Y] int age;

  /* offspring detected from each fish in each year */
  array[F, Y] int<lower=0> Bk = rep_array(0, F, Y);

  /* hatch year for each fish */
  array[F] int hatch_year;

  /* oldest age the fish is known to have reached */
  array[F] int<lower=1, upper=A> last_age;

  /* 
   * age of fish when it entered the study 
   * (this is 1 for fish hatched during the study)
   */
  array[F] int<lower=1, upper=A> first_age;

  /* first and last observation year for each fish */
  array[F] int<lower=1, upper=Y> first_year = rep_array(1, F);
  array[F] int<lower=1, upper=Y> last_year  = rep_array(Y, F);

  /* effective juvenile sample in each year */
  array[Y] int Jsmp = rep_array(0, Y);

  /* known spawnners in each year */
  array[2, Y] int<lower=0> Sk = rep_array(0, 2, Y); 

  /* fish detections */
  array[F, M, Y] int<lower=0, upper=1> fish_det = rep_array(0, F, M, Y);


  /* 
   * calculate the age of each fish
   * in each of the study years
   */
  {
    array[F] int n_obs = rep_array(0, F);
    array[F, Y] real age_acc = rep_array(0, F, Y);
    for (i in 1:Obs) {
      if (obs_age[i] <= 0) continue;
      n_obs[obs_fish[i]] += 1;
      for (y in 1:Y) {
        age_acc[obs_fish[i], y] += obs_age[i] + y - obs_year[i];
      }
    }
    for (f in 1:F) {
      for (y in 1:Y) {
        age_acc[f, y] /= n_obs[f];
      }
    }
    age = to_int(round(age_acc));
  }

  /* compute the birth year for every fish */
  for (f in 1:F) {
    hatch_year[f] = 2 - age[f, 1];
  }

  /* compute effective juvenile sample in each year */
  for (f in 1:F) {
    for (y in 1:Y) {
      Jsmp[y] += (age[f, y] == 1);
    }
  }

  /* 
   * compute the number of known spawnners in each year,
   * as well as the number of offspring detected for each
   * known fish, in each year
   */
  {
    array[F, Y] int unseen = rep_array(1, F, Y);
    for (i in 1:P) {
      int sex  = fish_sex[parent[i]];
      int year = hatch_year[offspring[i]];
      int fish = parent[i];
      if (year < 1 || year > Y) continue;
      Sk[sex, year] += unseen[fish, year];
      Bk[fish, year] += 1;
      unseen[fish, year] = 0;
    }
  }

  /* 
   * compute the first and last observation
   * for each fish
   */
  for (i in 1:Obs) {
    int fish = obs_fish[i];
    int year = obs_year[i];
    first_year[fish] = year < first_year[fish] ? year : first_year[fish];
    last_year[fish]  = year > last_year[fish]  ? year : last_year[fish];
  }
  for (i in 1:P) {
    int y = hatch_year[offspring[i]];
    int f = parent[i];
    first_year[f] = y < first_year[f] ? max(y, 1) : first_year[f];
    last_year[f]  = y > last_year[f]  ? min(y, Y) : last_year[f];
  }

  /* 
   * compute the age of each fish 
   * at the last observation 
   */
  for (f in 1:F) {
    last_age[f] = age[f, last_year[f]];
  }

  /* compute the age of each fish when it entered the study */
  for (f in 1:F) {
    first_age[f] = hatch_year[f] < 1 ? age[f, 1] : 1; 
  }

  /* compute detection matrix */
  for (i in 1:Obs) {
    int f = obs_fish[i];
    int y = obs_year[i];
    int m = obs_method[i];
    fish_det[f, m, y] = 1;
  }
  
}



parameters {
  /* logarithmic force of mortality */
  matrix<upper=0>[A, 2] log_fom; 

  /* detection probability for each method */
  vector<lower=0, upper=1>[M] met_prob;

  /* unknown spawners in each year */
  array[2, Y] real<lower=0> Su;
}



transformed parameters {
  /* Total spawners in each year */
  array[2, Y] real<lower=0> S;

  /* detection probability  */
  array[M, 2, A] real<lower=0, upper=1> det_prob;

  /* compute total number of spawners */
  for (s in 1:2) {
    for (y in 1:Y) {
      S[s,y] = Su[s, y] + Sk[s, y];
    }
  }

  /* compute detection probability  */
  for (m in 1:M) {
    for (s in 1:2) {
      for (a in 1:A) {
        det_prob[m, s, a] = met_prob[m] * app[m, s, a];
      }
    }
  }
}

model {
  array[2, Y] int J = {Jsmp, Jsmp};

  /* survival */
  {
    to_vector(log_fom) ~ cauchy(0, 5);
    vector[F] surv;
    for (f in 1:F) {
      surv[f] = sum(log_fom[1:last_age[f], fish_sex[f]]);
    }
    target += sum(surv);
  }


  /* detection */
  {
    for (f in 1:F) {
      int fa = max(1, age[f, 1]);
      int la = last_age[f];
      int na = la - fa + 1;
      int s = fish_sex[f];
      int fy = max(1, hatch_year[f]);
      int ly = last_year[f];
      /* direct detections */
      for (m in 1:M) {
        vector[na] p = to_vector(det_prob[m, s, fa:la]) ;
        fish_det[f, m, fy:ly] ~ bernoulli(p);
      }
      /* CKMR detections */
      for (y in fy:ly) {
        Bk[f, y] ~ binomial(J[s, y], 1/S[s, y]);
        J[s, y] -= Bk[f, y];
      }
    }
  }

  /* end of life history */
  for (f in 1:F) {
    int ny = Y - last_year[f];
    int s  = fish_sex[f];
    int fy = last_year[f] + 1;
    array[ny] int as = age[f, fy:Y]; 
    vector[ny] pnd = rep_vector(0, ny);
    vector[ny] chi; 
    if (ny <= 0) continue;
    /* non-detection by direct methods */
    for (m in 1:M) {
      pnd += log(1 - to_vector(det_prob[m, s, as]));
    }
    /* non-detection by CKMR */
    for (i in 1:ny) {
      int y = fy + i - 1;
      pnd[i] += binomial_lpmf(0 | J[s, y], 1/S[s, y]);  
    }
    for (i in 0:(ny-1)) {
      int y = Y - i;
      int j = ny - i;
      if (i == 0) {
        chi[j] = pnd[j];
      }
      else {
        real mu = exp(log_fom[age[f, y], s]);
	chi[j] = log_mix(mu, 0, pnd[j] + chi[j+1]);
      }
    }
    target += sum(chi);
  }
}
