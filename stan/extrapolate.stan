
// Gaussian Process model on latent f(x) with lower asymptote L:
// y ~ Normal( 1 - (1 - L) * inv_logit(f(x)), sigma_y )
// Include derivative pseudo-observations at x_deriv with targets deriv_target
// (use positive targets to encourage f' > 0, and near-zero targets on the tail to enforce flattening).

// Analytic checks implemented:
//  k(x,x')   = sf^2 * exp(- (x-x')^2 / (2 * ell^2))
//  d/dx' k   = k * ( (x - x') / ell^2 )
//  d/dx d/dx' k = k * ( 1/ell^2 - (x-x')^2 / ell^4 )

functions {
  // Radial kernel between scalars x and xp k(x,xp)
  real k_ff(real x, real xp, real sf, real ell) {
    real r2 = (x - xp) * (x - xp);
    return sf * sf * exp(-0.5 * r2 / (ell * ell));
  }

  // Covariance cov( f(x), f'(xp) ) = d/dxp k(x, xp)
  // Using analytic derivative:= k(x,xp) * ( (x - xp) / ell^2 )
  real k_fd(real x, real xp, real sf, real ell) {
    real r = x - xp;
    real k = k_ff(x, xp, sf, ell);
    return k * ( r / (ell * ell) );
  }

  // Covariance cov( f'(x), f'(xp) ) = d/dx d/dxp k(x,xp)
  // Using analytic derivative:= k(x,xp) * ( 1/ell^2 - (x - xp)^2 / ell^4 )
  real k_dd(real x, real xp, real sf, real ell) {
    real r = x - xp;
    real k = k_ff(x, xp, sf, ell);
    real ell2 = ell * ell;
    return k * ( (1.0 / ell2) - (r * r) / (ell2 * ell2) );
  }
}

data {
  int<lower=1> N;            // observed data points
  vector[N] x;               // observed pilot observation points
  vector[N] y;               // observed pilot delta values

  int<lower=1> Np;           // prediction locations
  vector[Np] x_pred;         // prediction inputs (grid for extrapolation)

  // Derivative pseudo-observations for shape constraint
  int<lower=1> Nd;           // number of derivative pseudo-observations (>=1)
  vector[Nd] x_deriv;        // derivative locations
  vector[Nd] deriv_target;   // target values for f'(x_deriv) 
  vector<lower=0>[Nd] sigma_deriv; // noise SDs for derivative pseudo-observations
}
transformed data {
  int L = N + Np + Nd;       // total latent vector length (f_obs, f_pred, f_deriv)
}
parameters {
  real<lower=0,upper=1> L_asym; // Lower asymptote
  real<lower=0> sf;          // GP marginal std-dev (sigma_f)
  real<lower=1e-6> ell;      // GP lengthscale
  real<lower=0> sigma_y;     // observation noise for y
  vector[L] u;               // joint latent vector: [f(x_obs); f(x_pred); f'(x_deriv)]
}
transformed parameters {
  vector[N] f_obs = u[1:N];                 // latent f at observed inputs
  vector[Np] f_pred = u[N + 1 : N + Np];    // latent f at prediction grid
  vector[Nd] f_deriv = u[N + Np + 1 : L];   // latent derivatives at x_deriv
}
model {
  // Priors
  // ----------------------
  // Weakly informative prior for asymptote L: assume strongly more likely to be close to 0
  L_asym ~ normal(-6, 2);

  // GP hyperpriors (weakly informative)
  sf ~ normal(0, 1.0);              // constrained positive 
  ell ~ lognormal(0, 1.0);          // positive
  sigma_y ~ normal(0, 0.5);         // constrained positive

  // ----------------------
  
  // Build GP covariance matrix K for joint vector u
  // block indices: 1..N => f(x), N+1..N+Np => f(x_pred), rest => derivatives
  // types: 0 => f, 1 => derivative
  {
    matrix[L, L] K;
    int types[L];
    real xs[L];

    // fill coordinates and types
    for (i in 1:N) {
      types[i] = 0;
      xs[i] = x[i];
    }
    for (i in 1:Np) {
      types[N + i] = 0;
      xs[N + i] = x_pred[i];
    }
    for (i in 1:Nd) {
      types[N + Np + i] = 1;
      xs[N + Np + i] = x_deriv[i];
    }

    // fill symmetric covariance matrix K using analytic expressions
    for (i in 1:L) {
      for (j in i:L) {
        real v;
        if (types[i] == 0 && types[j] == 0) {
          // cov(f(x_i), f(x_j))
          v = k_ff(xs[i], xs[j], sf, ell);
        } else if (types[i] == 0 && types[j] == 1) {
          // cov(f(x_i), f'(x_j)) = d/dx' k(x_i, x_j)
          v = k_fd(xs[i], xs[j], sf, ell);
        } else if (types[i] == 1 && types[j] == 0) {
          v = k_fd(xs[j], xs[i], sf, ell);
        } else {
          // cov(f'(x_i), f'(x_j)) = d/dx d/dx' k(x_i, x_j)
          v = k_dd(xs[i], xs[j], sf, ell);
        }
        if (i == j) v = v + 1e-9; // tiny jitter for numerical stability
        K[i, j] = v;
        K[j, i] = v;
      }
    }

    // Prior on joint latent u is zero-mean multivariate normal with covariance K
    u ~ multi_normal_cholesky(rep_vector(0, L), cholesky_decompose(K));
  }

  // Likelihood: mu = 1 - (1 - L_asym) * inv_logit(f_obs[n])
  for (n in 1:N) {
    real mu = 1.0 - (1.0 - L_asym) * inv_logit(f_obs[n]);
    y[n] ~ normal(mu, sigma_y);
  }


  // ----------------------
  // Derivative pseudo-observations
  // f_deriv[j] should be close to deriv_target[j] with sd sigma_deriv[j]
  // ----------------------
  for (j in 1:Nd) {
    f_deriv[j] ~ normal(deriv_target[j], sigma_deriv[j]);
  }
}
generated quantities {
  vector[Np] y_pred;
  for (i in 1:Np) {
    real mu_pred = 1.0 - (1.0 - L_asym) * inv_logit(f_pred[i]);
    y_pred[i] = mu_pred + normal_rng(0, sigma_y);
  }
}


