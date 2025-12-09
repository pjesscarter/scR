// sample_u_given_theta.stan
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
  int<lower=1> N;            // observed points
  vector[N] x;
  vector[N] y;

  int<lower=1> Np;
  vector[Np] x_pred;

  int<lower=1> Nd;
  vector[Nd] x_deriv;
  vector[Nd] deriv_target;
  vector<lower=0>[Nd] sigma_deriv;

  // fixed hyperparameters from posterior draw
  real<lower=0> sf;
  real<lower=1e-9> ell;
  real<lower=0> sigma_y;
  real logit_L;
}
transformed data {
  int Ltot = N + Np + Nd;
}
parameters {
  vector[Ltot] u;   // joint latent: f_obs, f_pred, f_deriv
}
transformed parameters {
  real<lower=0, upper=1> L_asym = inv_logit(logit_L);
  vector[N] f_obs = u[1:N];
  vector[Np] f_pred = u[N + 1 : N + Np];
  vector[Nd] f_deriv = u[N + Np + 1 : N + Np + Nd];
}
model {
  // prior on u is MVN(0, K) built from sf and ell
  {
    int types[Ltot];
    real xs[Ltot];
    matrix[Ltot, Ltot] K;

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

    for (i in 1:Ltot) {
      for (j in i:Ltot) {
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
        if (i == j) v = v + 1e-9;
        K[i, j] = v;
        K[j, i] = v;
      }
    }
    u ~ multi_normal_cholesky(rep_vector(0, Ltot), cholesky_decompose(K));
  }

  // likelihood from observed data (maps f_obs -> mu via L_asym and inv_logit)
  for (n in 1:N) {
    real mu = 1.0 - (1.0 - L_asym) * inv_logit(f_obs[n]);
    y[n] ~ normal(mu, sigma_y);
  }

  // derivative pseudo-observations
  for (j in 1:Nd) {
    f_deriv[j] ~ normal(deriv_target[j], sigma_deriv[j]);
  }
}
generated quantities {
  vector[Np] y_pred;
  for (i in 1:Np) {
    real mu_pred = 1.0 - (1.0 - inv_logit(logit_L)) * inv_logit(u[N + i]);
    y_pred[i] = mu_pred + normal_rng(0, sigma_y);
  }
}