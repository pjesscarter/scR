# fit_monotone_spline_epsilon.R
library(cmdstanr); library(posterior)
library(ggplot2); library(dplyr); library(tidyr)
library(pacman)
library(splines2)
p_load(stringr,
       scR,
       ggplot2,
       nlraa,
       nlme,
       dplyr,
       minpack.lm,
       parallel,
       pbapply
)
extract_num <- function(x,sample_size) {
  as.numeric(sub(paste(".*_(\\d+)_",sample_size,"\\.RDS",sep=""), "\\1", x))
}
set.seed(123)
datadir <- "../../../Dropbox/SampleComplexity/AJPS/Data/Mitts/boot"
benchmarkdir <- "../../../Dropbox/SampleComplexity/AJPS/Data/Mitts/full_boot"
files <- sort(list.files(datadir))
benchmark <- sort(list.files(benchmarkdir))
initial_sample_size <- 2000
# Sort based on extracted number
files <- files[str_detect(files,paste("_",initial_sample_size,".RDS",sep=""))]
files <- files[order(sapply(files, extract_num,sample_size=initial_sample_size))]
benchmark <- sapply(benchmark,function(x){paste(benchmarkdir,x,sep="/")})
benchmark <- lapply(benchmark,readRDS)
paths <- sapply(files,function(x){paste(datadir,x,sep="/")})
dat <- lapply(paths,readRDS)
outdir <- "../../Plots"
epsilon <- 0.05
delta <- 0.05
alpha <- 0.05
truedata <- bind_rows(lapply(benchmark, function(x) getpac(x, epsilon, delta)$Summary)) %>%
  group_by(n) %>%
  summarise(
    CI_Upper_Delta = quantile(Delta, 1 - alpha / 2, na.rm = TRUE),
    CI_Upper_Epsilon = quantile(Epsilon, 1 - alpha / 2, na.rm = TRUE),
    CI_Lower_Delta = quantile(Delta, alpha / 2, na.rm = TRUE),
    CI_Lower_Epsilon = quantile(Epsilon, alpha / 2, na.rm = TRUE),
    Delta = mean(Delta, na.rm = TRUE),
    Epsilon = mean(Epsilon, na.rm = TRUE)
  ) %>%
  rename(N = n)

upper_epsilon <-c(Asym = (epsilon-0.001),b2=Inf,b3=Inf)
# load data (same as before)
mit_boot <- readr::read_rds("../../../Dropbox/SampleComplexity/AJPS/Data/Mitts/full_boot/mit_boot_99_full.RDS")$Summary
pilot_n <- 2000
pilot <- mit_boot %>% filter(n <= pilot_n)

# Observations
x_obs <- as.numeric(pilot$n)
y_obs <- as.numeric(pilot$Epsilon)
N <- length(x_obs)

# set up prediction grid (you can control resolution)
N_grid <- 200                              # smaller => faster; increase for final
x_pred_full <- seq(min(x_obs), max(mit_boot$n), length.out = N_grid)

# choose knots for M-splines / I-splines
# we'll make interior knots across the observed range, plus boundary knots extending beyond max to allow slow decline
K_interior <- 8      # number of interior basis segments (adjust; 6-12 typical)
x_min <- min(x_obs)
x_max <- max(mit_boot$n)
# place interior knots between lower quantile and upper quantile of observed (to avoid too many left tail knots)
interior_knots <- seq(quantile(x_obs, 0.05), quantile(x_obs, 0.95), length.out = K_interior)
# choose boundary knots a bit beyond prediction max to allow slow decline (Option B)
right_extend <- 1.25   # extend 25% beyond maximum prediction
boundary_min <- x_min
boundary_max <- x_max * right_extend

# Build I-spline basis (I-splines are integrals of M-splines). Use degree = 3 (cubic)
degree <- 3
# Create the I-spline evaluation matrix at the grid and at observed points
I_grid <- iSpline(x_pred_full, knots = interior_knots, degree = degree,
                  intercept = TRUE, Boundary.knots = c(boundary_min, boundary_max))
I_obs  <- iSpline(x_obs,      knots = interior_knots, degree = degree,
                  intercept = TRUE, Boundary.knots = c(boundary_min, boundary_max))

# Convert to plain matrices (I_grid: N_grid x K)
I_grid_mat <- as.matrix(I_grid)
I_obs_mat  <- as.matrix(I_obs)
K <- ncol(I_grid_mat)   # number of basis functions


# write Stan file path and compile
stan_file <- "../Stan/extrapolate_epsilon_monotone_splines.stan"
mod <- cmdstan_model(stan_file)

stan_data <- list(
  N = length(y_obs),
  y = y_obs,
  N_grid = length(x_pred_full),
  x_grid = x_pred_full,
  N_obs = length(y_obs),
  obs_idx = sapply(x_obs, function(xx) which.min(abs(x_pred_full - xx))),
  K = ncol(I_grid_mat),
  I_grid = I_grid_mat,
  I_obs = I_obs_mat,
  range_r = max(x_pred_full) - min(x_pred_full),
  sd_y_prior = max(sd(y_obs, na.rm = TRUE), 0.01)
)

init_fun <- function() {
  K <- stan_data$K
  list(
    Logit_L = qlogis(0.001),
    kappa = 1.0 / max(1, stan_data$range_r),   # initial sensible scale
    sigma_y = sd(y_obs) / 3,
    f0 = 0.0,                                  # start near zero so mu ~ 0.5 but moveable
    mu_logalpha = 0.0,
    sigma_logalpha = 0.5,
    log_alpha = rep(log(0.05), K)             # small positive initial alpha: exp(log_alpha)=0.05
  )
}
# sampling settings (fast-ish for exploratory; increase iter for production)
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 800,
  iter_sampling = 800,
  adapt_delta = 0.99,    # tighter to reduce divergences
  max_treedepth = 15,
  init = init_fun
)

print(fit)

# Extract posterior predictive mu_grid median + 90% CI
draws <- as_draws_df(fit$draws(variables = c("mu_grid")))
mu_cols <- grep("^mu_grid\\[", names(draws), value = TRUE)
mu_mat <- as.matrix(draws[, mu_cols])
mu_med <- apply(mu_mat, 2, median)
mu_lo  <- apply(mu_mat, 2, quantile, probs = 0.05)
mu_hi  <- apply(mu_mat, 2, quantile, probs = 0.95)

# Assemble plotting data (pilot, benchmark, extrapolated)
pilot_df <- tibble(n = pilot$n, Epsilon = pilot$Epsilon, type = "Pilot")
bench_df <- mit_boot %>% filter(n > pilot_n) %>% transmute(n = n, Epsilon = Epsilon, type = "Benchmark")
pred_df  <- tibble(n = x_pred_full, Epsilon = mu_med, Epsilon_lo = mu_lo, Epsilon_hi = mu_hi, type = "Extrapolated")

# Plot similar to before
p <- ggplot() +
  geom_line(data = pred_df, aes(x = n, y = Epsilon, colour = "Extrapolated"), size = 0.9) +
  geom_ribbon(data = pred_df, aes(x = n, ymin = Epsilon_lo, ymax = Epsilon_hi),
              fill = "grey70", alpha = 0.25) +
  geom_point(data = pilot_df, aes(x = n, y = Epsilon, colour = "Pilot"), size = 1.2) +
  geom_point(data = bench_df, aes(x = n, y = Epsilon, colour = "Benchmark"), size = 1.0) +
  geom_vline(xintercept = max(pilot$n), linetype = "solid", colour = "black") +
  scale_colour_manual(name = "Type", values = c("Benchmark" = "#440154FF", "Extrapolated" = "#3B9983", "Pilot" = "#FDE725FF")) +
  theme_bw() +
  labs(x = "n", y = "Epsilon") +
  theme(legend.position = "right")

print(p)

# diagnostic: look at posterior median of alpha (shape of derivative)
alpha_draws <- as_draws_df(fit$draws(variables = "alpha"))
alpha_med <- apply(as.matrix(alpha_draws), 2, median)[1:K]
alpha_tbl <- tibble(basis = 1:K, alpha = alpha_med)
print(alpha_tbl)

# save fit if desired
# fit$save_object(file = "fit_monotone_spline_epsilon.rds")
# ggsave("monotone_spline_epsilon_fit.png", p, width = 10, height = 5, dpi = 200)
save.image("/Users/pjc504/Library/CloudStorage/Dropbox/SampleComplexity/AJPS/Data/Mitts/GPs/epsilon_monotone_spline.RData")
