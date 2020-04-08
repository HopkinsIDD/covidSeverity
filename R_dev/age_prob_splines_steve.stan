// generated with brms 2.11.1
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  // data for spline s(ageL,bs="cs")
  // int nb_1;  // number of bases
  int knots_1;  // number of knots
  // basis function matrices
  matrix[N, knots_1] Bmat;
  // matrix[N, knots_1[1]] Bmat;
  // data for group-level effects of ID 1
  int<lower=1> J;  // number of grouping levels
  int<lower=1> id[N];  // grouping indicator per observation
  row_vector<lower=0>[J] id_wt; // weights for predicting out of sample
  int<lower=0> preds; // number of groups to predict
  // vector<lower=0, upper=99>[preds] pred_min; // lower bound of groups to predict
  // vector<lower=0, upper=99>[preds] pred_max; // upper bound of groups to predict
  matrix[preds, knots_1] Bmat_preds;
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept;
  // parameters for spline s(ageL,bs="cs")
  // standarized spline coefficients
  vector[knots_1] spline_coef;
  // standard deviations of the coefficients
  real coef_mu;
  real<lower=0> coef_sd;
  real group_mu; // overall effect
  real<lower=0> group_sd;  // group-level standard deviations
  // standardized group-level effects
  // vector[J] z_1;
  vector[J] beta_j;
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1] s_1_1 = coef_sd * spline_coef;
  // actual group-level effects
  vector[J] r_1_1 = (group_sd * (beta_j)); //should that [1] be there?
}
model {
  // initialize linear predictor term
  // vector[N] mu = Intercept + rep_vector(0, N) + Bmat * spline_coef;
  vector[N] mu = Intercept + rep_vector(0, N) + Bmat * s_1_1;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[id[n]];
  }
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 0, 10);
  target += normal_lpdf(beta_j | group_mu, 1);
  target += normal_lpdf(spline_coef | coef_mu, coef_sd);
  target += normal_lpdf(spline_coef | 0, 1);
  target += student_t_lpdf(coef_sd | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(group_sd | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  // target += normal_lpdf(z_1[1] | 0, 1);
  // likelihood including all constants
  target += bernoulli_logit_lpmf(Y | mu);
  // Intercept ~ student_t(3, 0, 10);
  // beta_j ~ normal(group_mu, group_sd);
  // spline_coef ~ normal(coef_mu, coef_sd);
  // Y ~ bernoulli_logit(mu);
}
generated quantities {
  vector<lower=0,upper=1>[preds] age_preds = inv_logit(Intercept + rep_vector(0, preds) + Bmat_preds * spline_coef + group_mu);
  // // vector[preds] age_grp;
  // //
  // // for(i in 1:preds){
  // //   age_grp[i] = (pred_min[i], pred_max[i]);
  // // }
  // age_preds = ;
}
