// extrapolate_epsilon_monotone.stan
// Monotone GP via integrated nonnegative process
// g(x) = softplus(h(x)) >= 0, h ~ GP(0, k)
// f(x) = f0 + integral_0^x g(t) dt  (approx via trapezoid on grid)
// mu(x) = L + (0.5 - L) * inv_logit(-f(x))
// y ~ Normal(mu, sigma_y)

functions {
  real k_ff(real x, real xp, real sf, real ell) {
    real r2 = (x - xp) * (x - xp);
    return sf * sf * exp(-0.5 * r2 / (ell * ell));
  }
}

data {
  int<lower=1> N;                // # observed points
  vector[N] x_obs;               // observed input locations
  vector[N] y;                   // observed outputs

  int<lower=1> M;                // size of computation grid (prediction + integration)
  vector[M] x_grid;              // strictly increasing grid covering observed + prediction domain
  array[N] int<lower=1, upper=M> obs_idx; // mapping from observed rows to grid indices

  // prior helpers
  real<lower=0> range_r;         // domain range = max(x_grid)-min(x_grid)
  real<lower=0> sd_y_prior;      // typical observation noise scale
}

transformed data {
  // grid increments for trapezoid rule
  vector[M-1] dx;
  for (m in 1:(M-1)) dx[m] = x_grid[m+1] - x_grid[m];
}

parameters {
  real Logit_L;
  real<lower=0> sf_h;
  real<lower=1e-3> ell_h;
  real<lower=0> sigma_y;
  real f0;
  real<lower=0> kappa;      // NEW: scale controlling how f affects mu
  vector[M] z_h;
}

transformed parameters {
  real<lower=0, upper=1> L_asym = inv_logit(Logit_L);
  matrix[M, M] K_h;              // covariance for h on grid
  matrix[M, M] Lchol_h;
  vector[M] h;                   // latent field
  vector[M] g;                   // nonnegative derivative approx via softplus
  vector[M] f_grid;              // integrated f on grid (f' = g >= 0)
  real jitter;

  // build covariance K_h (symmetric)
  jitter = 1e-9 + 1e-8 * sf_h * sf_h;
  for (i in 1:M) {
    for (j in i:M) {
      real v = k_ff(x_grid[i], x_grid[j], sf_h, ell_h);
      if (i == j) v = v + jitter;
      K_h[i,j] = v;
      K_h[j,i] = v;
    }
  }

  // cholesky and non-centered transform
  Lchol_h = cholesky_decompose(K_h);
  h = Lchol_h * z_h;

  // softplus: numerically stable softplus = log1p(exp(h))
  for (m in 1:M) {
    g[m] = log1p(exp(h[m]));   // strictly > 0
  }

  // integrate g to get f_grid
  // set f_grid[1] = f0 (user-provided prior handles its plausible range)
  f_grid[1] = f0;
  for (m in 2:M) {
    f_grid[m] = f_grid[m-1] + 0.5 * (g[m-1] + g[m]) * dx[m-1];
  }
}

model {
  Logit_L ~ normal(-6, 1.3);
  sf_h ~ normal(0, 0.8);             // tighter amplitude prior
  ell_h ~ lognormal(log(range_r * 0.20), 0.45); // make lengthscale a bit larger on average
  sigma_y ~ normal(sd_y_prior, sd_y_prior * 0.75);
  kappa ~ lognormal(log(1.0 / range_r + 1e-6), 0.8); // center ~ 1/range (scales f->mu)
  z_h ~ normal(0, 1);
  f0 ~ normal(-1.0, 1.0);            // start f small negative (avoid huge immediate drop)

  // LIKELIHOOD
  for (n in 1:N) {
  real f_n = f_grid[ obs_idx[n] ];
  real mu = L_asym + (0.5 - L_asym) * exp( - kappa * f_n ); // changed mapping
  y[n] ~ normal(mu, sigma_y);
}
}

generated quantities {
  vector[M] mu_grid;
  vector[M] y_pred_grid;
  for (m in 1:M) mu_grid[m] = L_asym + (0.5 - L_asym) * exp( - kappa * f_grid[m] );
  for (m in 1:M) y_pred_grid[m] = normal_rng(mu_grid[m], sigma_y);

}
