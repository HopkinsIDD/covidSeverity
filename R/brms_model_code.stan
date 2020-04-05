// generated with brms 2.11.1
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  // data for spline s(ageL,bs="cs")
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept;
  // parameters for spline s(ageL,bs="cs")
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  // standard deviations of the coefficients
  real<lower=0> sds_1_1;
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  // standardized group-level effects
  vector[N_1] z_1[M_1];
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1 = sds_1_1 * zs_1_1;
  // actual group-level effects
  vector[N_1] r_1_1 = (sd_1[1] * (z_1[1]));
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0, N) + Zs_1_1 * s_1_1;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
  }
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 0, 10);
  target += normal_lpdf(zs_1_1 | 0, 1);
  target += student_t_lpdf(sds_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_1[1] | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    target += bernoulli_logit_lpmf(Y | mu);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
