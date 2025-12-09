// gp_decreasing_asymptote_nc.stan
// GP on latent f(x) with fixed upper asymptote 1 and unknown lower asymptote L in [0,1].
// Non-centered parameterization for GP latent (much better HMC geometry).
//
// Model:
//   y ~ Normal( 1 - (1 - L) * inv_logit(f(x)), sigma_y )
// derivative pseudo-observations: f'(x_deriv) ~ Normal( deriv_target, sigma_deriv )

functions {
  // Squared-exponential (RBF) kernel between scalars
  real k_ff(real x, real xp, real sf, real ell) {
    real r2 = (x - xp) * (x - xp);
    return sf * sf * exp(-0.5 * r2 / (ell * ell));
  }

  // cov( f(x), f'(xp) ) = d/dxp k(x,xp) = k(x,xp) * (x - xp) / ell^2
  real k_fd(real x, real xp, real sf, real ell) {
    real r = x - xp;
    real k = k_ff(x, xp, sf, ell);
    return k * ( r / (ell * ell) );
  }

  // cov( f'(x), f'(xp) ) = d/dx d/dxp k(x,xp)
  real k_dd(real x, real xp, real sf, real ell) {
    real r = x - xp;
    real k = k_ff(x, xp, sf, ell);
    real ell2 = ell * ell;
    return k * ( (1.0 / ell2) - (r * r) / (ell2 * ell2) );
  }
}

data {
  int<lower=1> N;            // data points
  vector[N] x;               // observed inputs (original n; can be scaled in R if you prefer)
  vector[N] y;               // observed outputs (delta in [0,1])

  int<lower=1> Nd;           // number of derivative pseudo-observations
  vector[Nd] x_deriv;        // derivative locations
  vector[Nd] deriv_target;   // target values for f'(x_deriv)
  vector<lower=1e-9>[Nd] sigma_deriv; // noise SDs for derivative pseudo-obs

  // helpers from R for scale-aware priors
  real<lower=0> range_r;     // = max(x) - min(x)
  real<lower=0> sd_y_prior;  // e.g. max(sd(y), 0.01)
}

transformed data {
  int L = N + Nd;            // dimension of joint latent [f_obs; f_deriv]
}

parameters {
  real Logit_L;                      // unconstrained; L_asym = inv_logit(Logit_L)
  real<lower=0> sf;                  // GP amplitude
  real<lower=1e-3> ell;              // lengthscale > 1e-3
  real<lower=0> sigma_y;             // observation noise

  vector[L] z_u;                     // non-centered latent (standard normal)
}

transformed parameters {
  real<lower=0, upper=1> L_asym = inv_logit(Logit_L);
  vector[N] f_obs;
  vector[Nd] f_deriv;

  // Build covariance K for the L = N + Nd joint vector ([f_obs; f_deriv])
  {
    matrix[L, L] K;
    array[L] int types;
    array[L] real xs;

    // fill coordinates and types
    for (i in 1:N) {
      types[i] = 0;
      xs[i] = x[i];
    }
    for (i in 1:Nd) {
      types[N + i] = 1;
      xs[N + i] = x_deriv[i];
    }

    // fill symmetric covariance matrix K using analytic expressions
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
        if (i == j) v = v + 1e-6;  // jitter for numerical stability
        K[i,j] = v;
        K[j,i] = v;
      }
    }

    // Cholesky for non-centered transform (add small jitter to diagonal already)
    matrix[L, L] Lchol = cholesky_decompose(K);

    // Non-centered parameterization: u = Lchol * z_u
    vector[L] u = Lchol * z_u;

    // split u into observed f and derivative latents
    for (i in 1:N) f_obs[i] = u[i];
    for (j in 1:Nd) f_deriv[j] = u[N + j];
  }
}

model {
  // Priors (data-aware)
  Logit_L ~ normal(-6, 1.3);                      // favors L near 0
  sf ~ normal(0, 3);                              // half-normal on positive sf
  ell ~ lognormal(log(0.5 * range_r), 0.3);
     // lengthscale prior centered on 0.5*range
  sigma_y ~ normal(0, sd_y_prior);                // half-normal for observation noise

  // Prior for non-centered latent
  z_u ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    real mu = 1.0 - (1.0 - L_asym) * inv_logit(f_obs[n]);
    y[n] ~ normal(mu, sigma_y);
  }

  // derivative pseudo-observations (soft constraints)
  for (j in 1:Nd) {
    f_deriv[j] ~ normal(deriv_target[j], sigma_deriv[j]);
  }
}


