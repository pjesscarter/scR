// extrapolate_epsilon_checked.stan
// GP model for a decreasing curve that starts near 0.5 at the left edge
// and flattens to a lower asymptote L (≈0 but possibly >0).
// Non-centered parameterization for the GP latent.

functions {
  real k_ff(real x, real xp, real sf, real ell) {
    real r2 = (x - xp) * (x - xp);
    return sf * sf * exp(-0.5 * r2 / (ell * ell));
  }
  real k_fd(real x, real xp, real sf, real ell) {
    real r = x - xp;
    real k = k_ff(x, xp, sf, ell);
    return k * ( r / (ell * ell) );
  }
  real k_dd(real x, real xp, real sf, real ell) {
    real r = x - xp;
    real k = k_ff(x, xp, sf, ell);
    real ell2 = ell * ell;
    return k * ( (1.0 / ell2) - (r * r) / (ell2 * ell2) );
  }
}

data {
  int<lower=1> N;              // observed points
  vector[N] x;                 // observed inputs
  vector[N] y;                 // observed outputs (in [0,1] typically)

  int<lower=0> Nd;             // number of derivative pseudo-observations
  vector[Nd] x_deriv;          // derivative locations
  vector[Nd] deriv_target;     // target values for f'(x_deriv)  (should be > 0 to make mu decrease)
  vector<lower=1e-9>[Nd] sigma_deriv;

  // optional soft anchor: index of observed point to anchor near 0.5 (left edge)
  int<lower=0> idx_anchor;     // 0 = disabled, otherwise index in 1..N of observed x to anchor
  real anchor_f_target;        // e.g., -6 or -8 to push inv_logit(-f) ~ 1 => mu ~ 0.5
  real<lower=0> sigma_anchor;  // small e.g., 0.5

  // helpers for priors
  real<lower=0> range_r;       // domain range (used to set ell prior)
  real<lower=0> sd_y_prior;    // typical noise scale for y
}

transformed data {
  int L = N + Nd;              // joint latent dimension: [f_obs; f_deriv]
}

parameters {
  real Logit_L;                // inverse-logit gives lower asymptote in [0,1]
  real<lower=0> sf;            // GP amplitude
  real<lower=1e-3> ell;        // lengthscale
  real<lower=0> sigma_y;       // obs noise

  vector[L] z_u;               // non-centered latent (standard normal)
}

transformed parameters {
  real<lower=0, upper=1> L_asym = inv_logit(Logit_L);
  vector[N] f_obs;
  vector[Nd] f_deriv;

  {
    // build K for [f_obs; f_deriv]
    matrix[L, L] K;
    array[L] int types;
    array[L] real xs;

    for (i in 1:N) {
      types[i] = 0;
      xs[i] = x[i];
    }
    for (i in 1:Nd) {
      types[N + i] = 1;
      xs[N + i] = x_deriv[i];
    }

    for (i in 1:L) {
      for (j in i:L) {
        real v;
        if (types[i] == 0 && types[j] == 0) {
          v = k_ff(xs[i], xs[j], sf, ell);
        } else if (types[i] == 0 && types[j] == 1) {
          v = k_fd(xs[i], xs[j], sf, ell);
        } else if (types[i] == 1 && types[j] == 0) {
          v = k_fd(xs[j], xs[i], sf, ell);
        } else {
          v = k_dd(xs[i], xs[j], sf, ell);
        }
        if (i == j) {
          // scaled jitter: small absolute floor plus a fraction of sf^2
          real jitter = 1e-9 + 1e-8 * sf * sf;
          v = v + jitter;
        }
        K[i, j] = v;
        K[j, i] = v;
      }
    }

    // cholesky on jittered K
    matrix[L, L] Lchol = cholesky_decompose(K);
    vector[L] u = Lchol * z_u;   // non-centered transform

    for (i in 1:N) f_obs[i] = u[i];
    for (i in 1:Nd) f_deriv[i] = u[N + i];
  }
}

model {
  // PRIORS (data-aware)
  Logit_L ~ normal(-6, 1.3);       // favor small L but allow moderate values

  // amplitude prior: prefer modest amplitudes (reduce extreme-K risk)
  sf ~ normal(0, 1.5);

  // lengthscale: data-aware prior (centered on a fraction of the range)
  // keep moderate spread but avoid vanishing ell
  ell ~ lognormal(log(range_r * 0.15), 0.5);

  // observation noise: center on observed sd to avoid mass at zero
  sigma_y ~ normal(sd_y_prior, sd_y_prior);

  // non-centered z
  z_u ~ normal(0, 1);

  // Likelihood: mu(x) = L + (0.5 - L) * inv_logit(-f(x))
  for (n in 1:N) {
    real mu = L_asym + (0.5 - L_asym) * inv_logit(-f_obs[n]);
    y[n] ~ normal(mu, sigma_y);
  }

  // derivative pseudo-obs (encourage f' > 0 so mu decreases)
  for (j in 1:Nd) {
    f_deriv[j] ~ normal(deriv_target[j], sigma_deriv[j]);
  }

  // optional soft anchor
  if (idx_anchor > 0) {
    f_obs[idx_anchor] ~ normal(anchor_f_target, sigma_anchor);
  }
}
