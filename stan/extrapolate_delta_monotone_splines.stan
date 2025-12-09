// monotone_spline_delta.stan
// Monotone integrated spline for delta curve (sigmoid-like drop from 1 -> L).
data {
  int<lower=1> N;                 // observed points
  vector[N] y;                    // observed outputs (Delta)
  int<lower=1> N_grid;            // grid size for predictions / integration
  vector[N_grid] x_grid;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1> obs_idx;
  real<lower=0> kappa0;  // fixed kappa provided from R
  int<lower=1> K;                 // number of spline basis functions
  matrix[N_grid, K] I_grid;       // I-spline basis evaluated on grid
  matrix[N_obs, K] I_obs;         // I-spline basis at observed x positions

  real<lower=0> range_r;
  real<lower=0> sd_y_prior;
}

parameters {
  real Logit_L;                   // lower asymptote L via inv_logit
  real log_kappa;            // scale mapping f -> inv_logit
  real<lower=0> sigma_y;
  real f0;

  // log-alpha parameterization to improve sampling
  vector[K] log_alpha;
  real mu_logalpha;
  real<lower=0> sigma_logalpha;
}

transformed parameters {
  real<lower=0, upper=1> L_asym = inv_logit(Logit_L);
  real<lower=0> kappa = exp(log_kappa);
  vector[N_grid] f_grid;
  vector[N_obs] f_obs;
  vector[K] alpha;
  alpha = exp(log_alpha);

  f_grid = rep_vector(f0, N_grid) + I_grid * alpha;
  for (n in 1:N_obs) f_obs[n] = f_grid[ obs_idx[n] ];
}

model {
  // PRIORS
  Logit_L ~ normal(-6, 1);   // favors small lower asymptote near 0
  // kappa prior centered at 1/range to scale f->inv_logit
  log_kappa ~ normal( log( kappa0 ), 0.45 );
  sigma_y ~ normal(sd_y_prior, sd_y_prior * 0.6);

  mu_logalpha ~ normal(log(0.08), 0.8);
  sigma_logalpha ~ normal(0.5, 0.4);
  log_alpha ~ normal(mu_logalpha, sigma_logalpha);

  // f0 strongly negative to place left mu near 1 by default
  f0 ~ normal(-6.0, 0.6);

  // LIKELIHOOD: mu = 1 - (1-L) * inv_logit(kappa * f)
  for (n in 1:N_obs) {
    real f_scaled_n = f_obs[n] / range_r;  // range_r supplied in data
    real mu_n = 1.0 - (1.0 - L_asym) * inv_logit( kappa * f_scaled_n );
    y[n] ~ normal(mu_n, sigma_y);
  }
}

generated quantities {
  vector[N_grid] mu_grid;
  vector[N_grid] y_pred_grid;
  for (m in 1:N_grid) {
    mu_grid[m] = 1.0 - (1.0 - inv_logit(Logit_L)) * inv_logit( kappa * (f_grid[m] / range_r) );
    y_pred_grid[m] = normal_rng(mu_grid[m], sigma_y);
  }
}
