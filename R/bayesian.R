library(readr); library(ggplot2); library(viridis)
library(cmdstanr)
library(parallel)
library(dplyr)
library(posterior)
library(zoo)

options(mc.cores = detectCores())

mit_boot <- read_rds("../../../Dropbox/SampleComplexity/AJPS/Data/Mitts/full_boot/mit_boot_99_full.RDS")$Summary
pilot_n <- 2000
pilot <- mit_boot %>% filter(n <= pilot_n)

ggplot(mit_boot,aes(x=n,y=Delta)) + geom_point()
ggplot(pilot,aes(x=n,y=Delta)) + geom_point()

x <- as.numeric(pilot$n)
y <- as.numeric(pilot$Delta)
N <- length(x)
x_min <- min(x)
x_max <- max(x)

x_pred <- seq(x_min, max(mit_boot$n), length.out = 300)
Np <- length(x_pred)

Nd_inside <- 4
Nd_tail <- 4

estimate_empirical_fprime <- function(df, x_col = "n", y_col = "delta",
                                      smooth_span = NULL,
                                      S_lower = 0.05, S_upper = 0.95,
                                      loess_control = NULL){
  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  S_obs <- pmin(pmax(1 - y, 0), 1)
  if (is.null(smooth_span)) {
    n <- length(x)
    smooth_span <- min(0.9, max(0.2, 3 / sqrt(n)))
  }
  if(!is.null(loess_control)){
    loess_fit <- do.call(stats::loess, c(list(formula = S_obs ~ x, span = smooth_span, control = loess_control), list()))
  } else{
    loess_fit <- do.call(stats::loess, c(list(formula = S_obs ~ x, span = smooth_span), list()))
  }
  x_grid <- seq(min(x), max(x), length.out = max(200, length(x)*5))
  S_smooth <- predict(loess_fit, newdata = data.frame(x = x_grid))
  if (any(is.na(S_smooth))) {
    S_smooth <- zoo::na.approx(S_smooth, x = x_grid, na.rm = FALSE)
    S_smooth[is.na(S_smooth)] <- approx(x, S_obs, xout = x_grid, rule = 2)$y[is.na(S_smooth)]
  }
  dS_dx <- c(diff(S_smooth) / diff(x_grid), NA)
  stable_idx <- which(S_smooth >= S_lower & S_smooth <= S_upper)
  stable_idx <- stable_idx[!is.na(dS_dx[stable_idx])]
  if (length(stable_idx) == 0) {
    warning("No stable region for S found (S in [0.05,0.95]) — estimator may be unstable.")
  }
  denom <- S_smooth * (1 - S_smooth)
  denom[denom < 1e-6] <- NA
  fprime_est <- rep(NA_real_, length(x_grid))
  good <- !is.na(dS_dx) & !is.na(denom)
  fprime_est[good] <- dS_dx[good] / denom[good]
  stable_vals <- fprime_est[stable_idx]
  stable_vals <- stable_vals[is.finite(stable_vals)]
  fprime_peak_est <- if (length(stable_vals) > 0) max(stable_vals, na.rm = TRUE) else NA_real_
  list(x_grid = x_grid,
       S_smooth = S_smooth,
       S_deriv = dS_dx,
       fprime_est = fprime_est,
       fprime_peak_est = fprime_peak_est,
       loess_fit = loess_fit)
}


make_pseudo_deriv_from_empirical <- function(df, x_col = "n", y_col = "Delta",
                                             xp_extra_factor = 1.5,
                                             Nd_inside = 4, Nd_tail = 4,
                                             safety_mult = 1.2,
                                             min_fprime = 1e-3,
                                             interior_fraction = c(0.5, 0.8),
                                             weak_factor = 0.5) {
  est <- estimate_empirical_fprime(df, x_col = x_col, y_col = y_col)
  fpeak_obs <- est$fprime_peak_est
  if (is.na(fpeak_obs) || !is.finite(fpeak_obs) || fpeak_obs <= 0) {
    fpeak_obs <- min_fprime
  }
  fprime_peak <- max(min_fprime, safety_mult * fpeak_obs)
  xr <- range(df[[x_col]])
  x_min <- xr[1]; x_max <- xr[2]; range_r <- x_max - x_min
  q_low <- interior_fraction[1]; q_high <- interior_fraction[2]
  interior_x <- seq(quantile(df[[x_col]], q_low), quantile(df[[x_col]], q_high), length.out = Nd_inside)
  tail_x <- seq(x_max, x_max * xp_extra_factor, length.out = Nd_tail)
  x_deriv <- c(interior_x, tail_x)
  if (!is.na(est$fprime_peak_est)) {
    idx_peak <- which.max(est$fprime_est)
    peak_x <- est$x_grid[idx_peak]
  } else {
    peak_x <- median(df[[x_col]])
  }
  sigma_bump <- max(range_r * 0.15, 0.1 * (x_max - x_min))
  deriv_target_interior <- fprime_peak * exp(- (interior_x - peak_x)^2 / (2 * sigma_bump^2))
  deriv_target_tail <- rep(0, length(tail_x))
  deriv_target <- c(deriv_target_interior, deriv_target_tail)
  sigma_deriv_interior <- pmax(weak_factor * deriv_target_interior, 1e-3)
  sigma_deriv_tail <- rep(max(0.02, 0.5 * median(sigma_deriv_interior, na.rm = TRUE)), length(tail_x))
  sigma_deriv <- c(sigma_deriv_interior, sigma_deriv_tail)
  list(x_deriv = x_deriv,
       deriv_target = deriv_target,
       sigma_deriv = sigma_deriv,
       fprime_peak = fprime_peak,
       peak_x = peak_x,
       empirical = est)
}


Nd_inside <- 6        # interior pseudo-deriv points (in observed region)
Nd_tail <- 12         # how many pseudo-deriv points to place across the extrapolation region
safety_mult <- 1.2    # allow a bit more slope than observed peak
weak_factor <- 0.05   # interior sigma = weak_factor * target  (smaller => stronger)
sigma_deriv_tail_val <- 0.005  # tight noise for tail points to enforce flattening
xp_extra_max <- max(x_pred)    # maximum x for extrapolation domain (we already computed x_pred)

# Use your existing empirical-based function to get an interior design
pseudo_info <- make_pseudo_deriv_from_empirical(
  pilot,
  x_col = "n", y_col = "Delta",
  xp_extra_factor = 1.0,      # interior function will place its tail points at x_max by default
  Nd_inside = Nd_inside,
  Nd_tail = 0,                # ask the helper for only interior points (we'll build tail grid separately)
  safety_mult = safety_mult,
  min_fprime = 1e-4,
  interior_fraction = c(0.45, 0.95),
  weak_factor = weak_factor
)

# interior results
x_deriv_interior <- pseudo_info$x_deriv
deriv_target_interior <- pseudo_info$deriv_target
sigma_deriv_interior <- pseudo_info$sigma_deriv

# Now create a future grid of derivative pseudo-observations that spans from just beyond observed max
x_obs_max <- max(x)
# place tail points starting slightly beyond observed region, up to xp_extra_max (prediction max)
tail_start <- x_obs_max + 0.01 * (x_obs_max - min(x))   # small buffer beyond last observed
if (tail_start >= xp_extra_max) {
  # If prediction max is same as observed max (unlikely), just put a few points beyond it
  tail_start <- x_obs_max + 0.1 * (x_obs_max + 1)
}
x_deriv_tail <- seq(tail_start, xp_extra_max, length.out = Nd_tail)

# tail derivative targets: strongly encourage flattening -> near 0 slope
deriv_target_tail <- rep(0, length(x_deriv_tail))
# tail sigma: small to enforce flattening but not absolutely hard
sigma_deriv_tail <- rep(sigma_deriv_tail_val, length(x_deriv_tail))

# Optionally: add a few soft derivative constraints near the expected peak region slightly beyond the observed region
# (This nudges the model to place the peak within the extrapolation band rather than before it.)
# Locate empirical peak location (if available from pseudo_info)
peak_x <- pseudo_info$peak_x
# One or two extra points near peak (optional)
extra_peak_x <- c(peak_x + 0.5 * (x_obs_max - peak_x)) # halfway toward observed max
extra_peak_target <- pseudo_info$fprime_peak * 0.6      # somewhat smaller than peak
extra_peak_sigma <- max(1e-3, 0.08 * extra_peak_target) # moderate precision

# assemble full vectors (order: interior then peak extras then tail)
x_deriv_full <- c(x_deriv_interior,
                  extra_peak_x,
                  x_deriv_tail)
deriv_target_full <- c(deriv_target_interior,
                       extra_peak_target,
                       deriv_target_tail)
sigma_deriv_full <- c(sigma_deriv_interior,
                      rep(extra_peak_sigma, length(extra_peak_x)),
                      sigma_deriv_tail)
# debug print
cat("Interior derivative points:", length(x_deriv_interior), "\n")
cat("Extra peak points:", length(extra_peak_x), "\n")
cat("Tail derivative points:", length(x_deriv_tail), "\n")
cat("Total pseudo-deriv points:", length(x_deriv_full), "\n")
print(data.frame(x = x_deriv_full, deriv_target = deriv_target_full, sigma = sigma_deriv_full))

# Build stan_data (note Nd is the total number)
stan_data <- list(
  N = N,
  x = x,
  y = y,
  Nd = length(x_deriv_full),
  x_deriv = as.numeric(x_deriv_full),
  deriv_target = as.numeric(deriv_target_full),
  sigma_deriv = as.numeric(sigma_deriv_full)
)

# Use the full-domain range for the prior (observed + deriv + predictions)
range_r <- max(c(x, x_deriv_full, x_pred)) - min(c(x, x_deriv_full, x_pred))
sd_y_prior <- max(sd(y, na.rm = TRUE), 0.01)
stan_data$range_r <- range_r
stan_data$sd_y_prior <- sd_y_prior

# Updated init function (z_u length uses new Nd)
init_fun <- function(chain_id = 1) {
  L <- N + stan_data$Nd
  base_range <- range_r
  base_ell <- max(base_range / 5, 1e-3)
  base_sf  <- 0.5
  base_sigma_y <- sd(y, na.rm = TRUE) / 2
  base_Logit_L <- qlogis(0.005)
  
  list(
    Logit_L = base_Logit_L + rnorm(1, 0, 0.1),
    sf = abs(base_sf + rnorm(1, 0, 0.1)),
    ell = abs(base_ell * exp(rnorm(1, 0, 0.1))),
    sigma_y = abs(base_sigma_y + rnorm(1, 0, base_sigma_y * 0.1)),
    z_u = rnorm(L, 0, 0.01)
  )
}


stan_file <- "../Stan/extrapolate.stan"
mod <- cmdstan_model(stan_file) 

fit_cmd <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.98,   # tighten a bit
  max_treedepth = 15,
  refresh = 250,
  init = init_fun
)

print(fit_cmd)
# helper kernels (same as Stan)
k_ff <- function(x, xp, sf, ell) {
  r2 <- (x - xp)^2
  sf^2 * exp(-0.5 * r2 / (ell^2))
}
k_fd <- function(x, xp, sf, ell) {
  r <- x - xp
  k <- k_ff(x, xp, sf, ell)
  k * (r / (ell^2))
}
k_dd <- function(x, xp, sf, ell) {
  r <- x - xp
  k <- k_ff(x, xp, sf, ell)
  ell2 <- ell^2
  k * ((1/ell2) - (r^2) / (ell2^2))
}
robust_mvnorm_sample <- function(mu, Sigma, jitter_start = 1e-12, max_jitter = 1e-2, warn = TRUE) {
  # ensure numeric matrix and symmetry
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  n <- length(mu)
  if (n == 0) return(numeric(0))
  if (!all(dim(Sigma) == c(n, n))) stop("robust_mvnorm_sample: dimension mismatch between mu and Sigma")
  
  # Try chol with increasing jitter
  eps <- jitter_start
  for (attempt in seq_len(10)) {
    # if eps exceeds max_jitter stop trying further
    if (eps > max_jitter) break
    ok <- tryCatch({
      R <- chol(Sigma + diag(eps, n))
      TRUE
    }, error = function(e) FALSE)
    if (ok) {
      # perform draw via chol factor R (upper triangular)
      R <- chol(Sigma + diag(eps, n))
      z <- rnorm(n)
      draw <- as.numeric(mu + t(R) %*% z)  # note: if R is upper chol, t(R) is lower
      return(draw)
    }
    eps <- eps * 10
  }
  
  # If chol attempts failed, fallback to eigen-decomposition with eigenvalue flooring
  eig <- eigen(Sigma, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  
  # Floor small negative values to zero, and very small positives to a small floor
  vals[!is.finite(vals)] <- 0
  vals_floor <- pmax(vals, 0)
  tol <- 1e-16
  # if all eigenvalues are effectively zero, return the mean
  if (all(vals_floor <= tol)) {
    if (warn) warning("robust_mvnorm_sample: covariance numerically zero; returning mean")
    return(as.numeric(mu))
  }
  
  # make numeric robust: floor very small eigenvalues to a tiny positive for sampling stability
  vals_floor <- pmax(vals_floor, 1e-12)
  sqrt_vals <- sqrt(vals_floor)
  z <- rnorm(length(sqrt_vals))
  draw <- as.numeric(mu + vecs %*% (sqrt_vals * z))
  return(draw)
}
draws_df <- as_draws_df(fit_cmd$draws(variables = c("Logit_L","sf","ell","sigma_y", "z_u")))
n_draws_total <- nrow(draws_df)
set.seed(8536)
n_use <- min(n_draws_total, 400)
idx <- sample(seq_len(n_draws_total), n_use)

z_names <- grep("^z_u\\[", names(draws_df), value = TRUE)
Luz <- length(z_names)

for_draws_y <- matrix(NA, nrow = n_use, ncol = length(x_pred))
x_obs <- x
Llen_obs <- length(c(x_obs, x_deriv_full))

for (k in seq_along(idx)) {
  cat(sprintf("Processing draw %i of %i\n", k, n_use))
  r <- idx[k]
  Logit_L <- draws_df$Logit_L[r]
  L_asym <- 1 / (1 + exp(-Logit_L))
  sf <- draws_df$sf[r]
  ell <- draws_df$ell[r]
  sigma_y <- draws_df$sigma_y[r]
  
  z_vec <- as.numeric(draws_df[r, z_names])
  if (length(z_vec) != (N + length(x_deriv_full))) stop("z_u length mismatch")
  
  Loo_x <- c(x_obs, x_deriv_full)
  types <- c(rep(0, length(x_obs)), rep(1, length(x_deriv_full)))
  K_oo <- matrix(0, Llen_obs, Llen_obs)
  for (i in 1:Llen_obs) {
    for (j in i:Llen_obs) {
      xi <- Loo_x[i]; xj <- Loo_x[j]
      if (types[i]==0 && types[j]==0) {
        v <- k_ff(xi, xj, sf, ell)
      } else if (types[i]==0 && types[j]==1) {
        v <- k_fd(xi, xj, sf, ell)
      } else if (types[i]==1 && types[j]==0) {
        v <- k_fd(xj, xi, sf, ell)
      } else {
        v <- k_dd(xi, xj, sf, ell)
      }
      if (i == j) v <- v + 1e-6
      K_oo[i,j] <- v
      K_oo[j,i] <- v
    }
  }
  
  chol_ok <- TRUE
  cholKoo <- tryCatch(chol(K_oo), error = function(e) { chol_ok <<- FALSE; NULL })
  if (!chol_ok) {
    eig <- eigen(K_oo, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-8)
    K_oo_pd <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
    cholKoo <- chol(K_oo_pd)
  }
  
  u_vec <- as.numeric(cholKoo %*% z_vec)
  
  u_obs_draw <- u_vec[1:length(x_obs)]
  u_deriv_draw <- u_vec[(length(x_obs)+1):length(u_vec)]
  
  K_po <- matrix(0, length(x_pred), Llen_obs)
  for (i in seq_along(x_pred)) {
    for (j in 1:Llen_obs) {
      xi <- x_pred[i]; xj <- Loo_x[j]
      if (types[j]==0) K_po[i,j] <- k_ff(xi, xj, sf, ell) else K_po[i,j] <- k_fd(xi, xj, sf, ell)
    }
  }
  K_pp <- outer(x_pred, x_pred, Vectorize(function(a,b) k_ff(a,b,sf,ell)))
  
  #solve K_oo^{-1} * u_o
  Koo_inv_u <- backsolve(cholKoo, forwardsolve(t(cholKoo), c(u_obs_draw, u_deriv_draw)))
  mu_fpred <- K_po %*% Koo_inv_u
  
  #compute K_po %*% K_oo^{-1} 
  solveKoo_tKpo <- backsolve(cholKoo, forwardsolve(t(cholKoo), t(K_po)))
  cov_fpred <- K_pp - K_po %*% solveKoo_tKpo
  cov_fpred <- (cov_fpred + t(cov_fpred)) / 2
  
  #draw f_pred (ensure positive definite)
  f_draw <- robust_mvnorm_sample(mu = as.numeric(mu_fpred), Sigma = cov_fpred)
  
  #transform to y
  Sdraw <- 1 / (1 + exp(-f_draw))
  mu_y <- 1 - (1 - L_asym) * Sdraw
  y_draw <- mu_y + rnorm(length(mu_y), 0, sigma_y)
  
  for_draws_y[k, ] <- y_draw
}

y_med <- apply(for_draws_y, 2, median)
y_lo <- apply(for_draws_y, 2, quantile, 0.35)
y_hi <- apply(for_draws_y, 2, quantile, 0.65)
# assemble plot data
# pilot rows (yellow)
pilot_df <- tibble(n = pilot$n, Delta = pilot$Delta, type = "Pilot")
# benchmark rows (purple) are all mit_boot rows that are not in pilot
bench_df <- mit_boot %>% filter(n > pilot_n) %>% transmute(n = n, Delta = Delta, type = "Benchmark")
# extrapolated (x_grid) rows
pred_df <- tibble(n = x_pred, Delta = y_med, Delta_lo = y_lo, Delta_hi = y_hi, type = "Extrapolated")

plot_df <- bind_rows(
  pilot_df,
  bench_df,
  pred_df %>% transmute(n, Delta, type)
)

# basic plot: lines for each type + ribbon for extrapolation
p <- ggplot() +
  geom_line(data = pred_df, aes(x = n, y = Delta, colour = "Extrapolated"), size = 0.8) +
  geom_ribbon(data = pred_df, aes(x = n, ymin = Delta_lo, ymax = Delta_hi), fill = "grey70", alpha = 0.25) +
  geom_line(data = pilot_df, aes(x = n, y = Delta, colour = "Pilot"), size = 1.1) +
  geom_line(data = bench_df, aes(x = n, y = Delta, colour = "Benchmark"), size = 0.8) +
  geom_vline(xintercept = max(pilot$n), linetype = "solid", colour = "black") +
  scale_colour_manual(name = "type", values = c("Benchmark" = "#440154FF", "Extrapolated" = "#3B9983", "Pilot" = "#FDE725FF")) +
  theme_bw() + labs(x = "n", y = "Delta") +
  theme(legend.position = "right")

print(p)

save.image("/Users/pjc504/Library/CloudStorage/Dropbox/SampleComplexity/AJPS/Data/Mitts/GPs/delta.RData")

#Epsilon
load("/Users/pjc504/Library/CloudStorage/Dropbox/SampleComplexity/AJPS/Data/Mitts/GPs/delta.RData")

