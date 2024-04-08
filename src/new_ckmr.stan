data {
  /* number of years in study */
  int<lower=2> Y;

  /* number of obsevation methods */
  int<lower=3> M;

  /* number of age groups */
  int<lower=2> A;

  /* number of observed fish */
  int<lower=T> F;

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
}

transformed data {
  /* Sex indices */
  int MALE = 1;
  int FEM  = 2;


  /* age of each fish in each study year */
  array[F, Y] int age;

  /* offspring detected from each fish in each year */
  array[F, Y] int<lower=0> Bk;

  /* hatch year for each fish */
  array[F] int hatch_year;

  /* first and last observation year for each fish */
  array[F] int<lower=1, upper=Y> first = rep_array(1, F);
  array[F] int<lower=1, upper=Y> last  = rep_array(Y, F);

  /* effective juvenile sample in each year */
  array[Y] int Jsmp = rep_array(0, Y);

  /* known spawnners in each year */
  array[2, Y] int<lower=0> Sk = rep_array(0, 2, Y); 

  /* fish detections */
  array[F, Y, M] int<lower=0, upper=1> fish_det = rep_array(0, F, Y, M);


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
    first[fish] = year < first[fish] ? year : first[fish];
    last[fish]  = year > last[fish]  ? year : last[fish];
  }

  /* compute detection matrix */
  for (i in 1:Obs) {
    int f = obs_fish[i];
    int y = obs_year[i];
    int m = obs_method[i];
    fish_det[f, y, m] = 1;
  }
  
}



parameters {
  /* logarithmic force of mortality */
  array[2, A] real<upper=0> fom; 

  /* spawn probability */
  array[2, A] real<lower=0, upper=1> sp;

  /* detection probability for each method */
  vector<lower=0, upper=1>[M] det_prob;

  /* unknown spawners in each year */
  array[2, Y] real<lower=0> Su;
}



transformed parameters {
  /* logarithmic survival to age probabilities */
  array[2, A] real<upper=0> surv;

  /* Total spawners in each year */
  array[2, Y] real<lower=0> S;

  /* CKMR detection probability */
  array[F, Y] real<lower=0, upper=1> ckmr_prob;

  /* compute survival to age probabilities */
  surv[MALE, :] = cumulative_sum(fom[MALE, :]); 
  surv[FEM,  :] = cumulative_sum(fom[FEM,  :]); 

  /* compute total number of spawners */
  for (s in 1:2) {
    for (y in 1:Y) {
      S[s,y] = Su[s, y] + Sk[s, y];
    }
  }
}

model {

  /* survival */
  for (f in 1:F) {
    /* survival probability to last observation for each fish */
    target += surv[fish_sex[f], age[f, last[f]]];
  }

  /* detection */
  for (f in 1:F) {
    for (y in hatch_year[f]:last[f]) {
      fish_det[f, y, :] ~ bernoulli(det_prob);
    }
  }

  /* drop-off */
  rep_array(1, F) ~ bernoulli(off_prob);

  /* CKMR */
  { 
    for (f in 1:F) {
      array[Y] int J = Jsmp;
      for (y in 1:Y) {
        if (Bk[f, y] < 1) continue;
	Bk[f, y] ~ binomial(J[y], ckmr_prob[y]);
	J[y] -= Bk[f, y];
    }
  }

}
