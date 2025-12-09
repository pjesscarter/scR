x_obs <- as.numeric(pilot$n)
y_obs <- as.numeric(pilot$Epsilon)
N <- length(x_obs)

# build grid (make sure grid contains the observed x exactly)
# we choose 300 points across the whole domain (observed + extrapolation)
x_pred_full <- seq(min(x_obs), max(mit_boot$n), length.out = 300)
x_grid <- sort(unique(c(x_obs, x_pred_full)))
M_try <- 120
x_grid <- sort(unique(c(x_obs, seq(min(x_obs), max(mit_boot$n), length.out = M_try))))
M <- length(x_grid)
obs_idx <- sapply(x_obs, function(xx) which.min(abs(x_grid - xx)))


init_fun <- function() {
  list(
    Logit_L = qlogis(0.001),
    sf_h = 0.6,
    ell_h = max(range_r / 6, 1e-3),
    sigma_y = sd(y_obs)/3,
    f0 = -1.0,
    kappa = 1.0 / max(1, range_r),
    z_h = rnorm(M, 0, 0.03)
  )
}
# helpers for priors
range_r <- max(x_grid) - min(x_grid)
sd_y_prior <- max(sd(y_obs, na.rm = TRUE), 0.01)

# build stan data
stan_data <- list(
  N = N,
  x_obs = x_obs,
  y = y_obs,
  M = M,
  x_grid = x_grid,
  obs_idx = as.integer(obs_idx),
  range_r = range_r,
  sd_y_prior = sd_y_prior
)

# compile model (path relative to this script)
stan_file <- "../Stan/extrapolate_epsilon_monotone.stan"
mod <- cmdstan_model(stan_file)


# fit: moderate length; increase iter if you want tighter posteriors
fit <- mod$sample(
  data = stan_data,
  seed = 2026,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 600,
  iter_sampling = 600,
  adapt_delta = 0.95,
  max_treedepth = 12,
  init = init_fun
)

print(fit)    # quick diagnostics (rhat, ess, etc.)

# Extract posterior mu_grid draws (generated quantities mu_grid)
draws_df <- as_draws_df(fit$draws(variables = "mu_grid"))
mu_cols <- grep("^mu_grid\\[", names(draws_df), value = TRUE)
mu_mat <- as.matrix(draws_df[, mu_cols])   # rows = draws, cols = grid positions (1..M)

# posterior summaries on grid
mu_med <- apply(mu_mat, 2, median)
mu_lo  <- apply(mu_mat, 2, quantile, probs = 0.05)
mu_hi  <- apply(mu_mat, 2, quantile, probs = 0.95)

# assemble plot data
# pilot rows (yellow)
pilot_df <- tibble(n = pilot$n, Delta = pilot$Epsilon, type = "Pilot")
# benchmark rows (purple) are all mit_boot rows that are not in pilot
bench_df <- mit_boot %>% filter(n > pilot_n) %>% transmute(n = n, Delta = Epsilon, type = "Benchmark")
# extrapolated (x_grid) rows
pred_df <- tibble(n = x_grid, Delta = mu_med, Delta_lo = mu_lo, Delta_hi = mu_hi, type = "Extrapolated")

plot_df <- bind_rows(
  pilot_df,
  bench_df,
  pred_df %>% transmute(n, Delta, type)
)

# basic plot: lines for each type + ribbon for extrapolation
p <- ggplot() +
  geom_line(data = pred_df, aes(x = n, y = Delta, colour = "Extrapolated"), size = 0.8,alpha=0.7) +
  geom_ribbon(data = pred_df, aes(x = n, ymin = Delta_lo, ymax = Delta_hi), fill = "grey70", alpha = 0.25) +
  geom_line(data = pilot_df, aes(x = n, y = Delta, colour = "Pilot"), size = 1.1,alpha=0.7) +
  geom_line(data = bench_df, aes(x = n, y = Delta, colour = "Benchmark"), size = 0.8,alpha=0.7) +
  geom_vline(xintercept = max(pilot$n), linetype = "solid", colour = "black") +
  scale_colour_manual(name = "type", values = c("Benchmark" = "#440154FF", "Extrapolated" = "#3B9983", "Pilot" = "#FDE725FF")) +
  theme_bw() + labs(x = "n", y = "Epsilon") +
  theme(legend.position = "right")

print(p)

save.image("/Users/pjc504/Library/CloudStorage/Dropbox/SampleComplexity/AJPS/Data/Mitts/GPs/epsilon_monotone.RData")
