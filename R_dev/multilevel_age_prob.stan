
data {
  int<lower=1> S;               // Number of subjects
  int<lower=1> J;               // the number of groups (Number of studies from which data comes)
  int<lower=1> A;               // max number of age groups (Number of studies from which data comes)
  int<lower=1,upper=A> agegrp[S];   // vector of age group ids (indices)
  int<lower=1,upper=J> id[S];   // vector of group ids (indices)
  // matrix[S,A] X;                // fixed effects model integer array
  int<lower=1> N[S];
  int<lower=0> Y[S];    // response variable - Outcome for each subect (died=1, survived=0)
}

// transformed data {
//   int<lower=0> K_real;
// }

parameters {
  real beta0;
  vector[A] beta;               // Fixed effects
  real<lower=0> sigma;    // SDs for random effect
  vector[J] beta_j;                  // Study Intercept (prob coefficient in each study)
}

model {
  // Priors on group parameters
  // beta0 ~ normal(0, 100);
  beta ~ normal(0, 100);

  beta_j ~ normal(beta0, sigma);

  // The likelihood
  Y ~ binomial_logit(N, beta_j[id] + beta[agegrp]);

}

generated quantities {         // simulate quantities of interest
  vector[A] prob_age;

  prob_age = inv_logit(beta0 + beta);

}

