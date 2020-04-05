
data {
  int<lower=1> N;               // Number of subjects
  int<lower=1> J;               // the number of groups (Number of studies from which data comes)
  int<lower=1,upper=J> id[N];   // vector of group ids (indices)
  // matrix[N,A] X;                // fixed effects model integer array
  row_vector<lower=0>[J] id_wt;
  int<lower=0, upper=1> Y[N];    // response variable - Outcome for each subect (died=1, survived=0)
  vector<lower=0, upper=99>[N] age_min;
  vector<lower=0, upper=99>[N] age_max;
  int<lower=0> preds; // number of preds
  vector<lower=0, upper=99>[preds] pred_min;
  vector<lower=0, upper=99>[preds] pred_max;
}

transformed data {
  vector<lower=0>[N] age;
  // vector<lower=0>[N] age2;
  // vector<lower=0>[N] age3;

  for(i in 1:N){
     age[i] = uniform_rng(age_min[i], age_max[i]);
     // age2[i] = age[i]^2;
     // age3[i] = age[i]^3;
  }
}

parameters {
  real beta0;
  real beta1;               // Fixed effects
  // real beta2;
  // real beta3;
  real<lower=0> sigma;    // SDs for random effect
  vector[J] beta_j;                  // Study Intercept (prob coefficient in each study)
}

model {
  // Priors on group parameters
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  // beta2 ~ normal(0, 100);
  // beta3 ~ normal(0, 100);

  beta_j ~ normal(beta0, sigma);

  // The likelihood
    Y ~ bernoulli_logit(beta_j[id] + beta1*age);
    // Y ~ bernoulli_logit(beta_j[id] + beta1*age + beta2*age2 + beta3*age3);
}

generated quantities {         // simulate quantities of interest
  vector[preds] age_grp;
  // vector[preds] age_grp2;
  // vector[preds] age_grp3;
  vector[preds] prob_age;
  // real log_lik[N];

  for(i in 1:preds){
    age_grp[i] = uniform_rng(pred_min[i], pred_max[i]);
    // age_grp2[i] = age_grp[i]^2;
    // age_grp3[i] = age_grp[i]^3;
  }
  prob_age = inv_logit(dot_product(id_wt, beta_j) + beta1*age_grp);
  // prob_age = inv_logit(dot_product(id_wt, beta_j) + beta1*age_grp + beta2*age_grp2 + beta3*age_grp3);

  // for(i in 1:N){
  //   // log_lik[i] = bernoulli_logit_lpmf(Y[i] | dot_product(id_wt, beta_j) + beta1*age[i]);
  //   log_lik[i] = bernoulli_logit_lpmf(Y[i] | dot_product(id_wt, beta_j) + beta1*age[i] + beta2*age2[i] + beta3*age3[i]);
  // }

}

