// monotone_spline_epsilon_v2.stan
data {
  int<lower=1> N;                 // observed points
  vector[N] y;                    // observed outputs (Epsilon)
  int<lower=1> N_grid;            // grid size
  vector[N_grid] x_grid;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1> obs_idx;

  int<lower=1> K;                 // number of spline basis functions
  matrix[N_grid, K] I_grid;       // I-spline evaluated on grid
  matrix[N_obs, K] I_obs;         // I-spline at observed locations

  real<lower=0> range_r;
  real<lower=0> sd_y_prior;
}

parameters {
  real Logit_L;
  real<lower=0> kappa;            // scale mapping f -> mu
  real<lower=0> sigma_y;
  real f0;

  // new: log-alpha parameterization (unconstrained)
  vector[K] log_alpha;            // alpha = exp(log_alpha) >= 0
  real<lower=0> sigma_logalpha;   // hierarchical scale for log_alpha
  real mu_logalpha;               // prior mean for log_alpha, lets us center alpha
}

transformed parameters {
  real<lower=0, upper=1> L_asym = inv_logit(Logit_L);
  vector[N_grid] f_grid;
  vector[N_obs] f_obs;
  vector[K] alpha;
  alpha = exp(log_alpha);         // positive derivative coefficients

  // f_grid = f0 + I_grid * alpha
  f_grid = rep_vector(f0, N_grid) + I_grid * alpha;
  for (n in 1:N_obs) f_obs[n] = f_grid[ obs_idx[n] ];
}

model {
  // PRIORS (more informative)
  Logit_L ~ normal(-5.5, 1.0);    // prefer small L but allow some mass
  // center kappa around 1/range_r (scale to convert f to mu)
  kappa ~ lognormal(log(1.0 / (range_r + 1e-12)), 0.40);

  sigma_y ~ normal(sd_y_prior, sd_y_prior * 0.6);

  // hierarchical log-alpha prior; allows multipliers > 1 if data demands
  mu_logalpha ~ normal(0.0, 1.0);               // center log-alpha around 0 by default
  sigma_logalpha ~ normal(0.5, 0.5);            // half-normal-ish; adjust if needed
  log_alpha ~ normal(mu_logalpha, sigma_logalpha);

  // prior for f0: center near 0 so left mu ≈ 0.5 by default but not stuck
  f0 ~ normal(0.0, 1.0);

  // LIKELIHOOD
  for (n in 1:N_obs) {
    real mu_n = L_asym + (0.5 - L_asym) * exp(- kappa * f_obs[n]);
    y[n] ~ normal(mu_n, sigma_y);
  }
}

generated quantities {
  vector[N_grid] mu_grid;
  vector[N_grid] y_pred_grid;
  for (m in 1:N_grid) {
    mu_grid[m] = L_asym + (0.5 - L_asym) * exp(- kappa * f_grid[m]);
    y_pred_grid[m] = normal_rng(mu_grid[m], sigma_y);
  }
}
