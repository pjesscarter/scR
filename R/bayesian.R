library(readr); library(ggplot2) #Not needed
library(rstan)
library(parallel)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = detectCores())
mit_boot <- read_rds("../../../Dropbox/SampleComplexity/AJPS/Data/Mitts/full_boot/mit_boot_99_full.RDS")$Summary
pilot_n <- 2000
pilot <- mit_boot %>% filter(n <= pilot_n)
x <- as.numeric(pilot$n)
y <- as.numeric(pilot$Delta)
N <- length(x)
x_min <- min(x)
x_max <- max(x)
# values to predict for
x_pred <- seq(x_min, max(mit_boot$n), length.out = 300)
Np <- length(x_pred)

Nd_inside <- 4
Nd_tail <- 4
#Provide pseudo-derivative observations for monotonicity and lower asymptotic behavior
x_deriv_inside <- seq(x_min + 0.05*(x_max - x_min), x_max - 0.05*(x_max - x_min), length.out = Nd_inside)
x_deriv_tail <- seq(x_max, x_max * 2.5, length.out = Nd_tail)
x_deriv <- c(x_deriv_inside, x_deriv_tail)
Nd <- length(x_deriv)

eps_inside <- 0.02   # positive small derivative target inside data
# Force near-zero derivatives in the tail (flattening)
eps_tail <- 1e-10

deriv_target <- c(rep(eps_inside, Nd_inside), rep(eps_tail, Nd_tail))
# Noise for derivative pseudo-obs (smaller noise means stronger enforcement)
sigma_deriv <- c(rep(0.02, Nd_inside), rep(0.01, Nd_tail))

stan_data <- list(
  N = N,
  x = x,
  y = y,
  Np = Np,
  x_pred = x_pred,
  Nd = Nd,
  x_deriv = x_deriv,
  deriv_target = deriv_target,
  sigma_deriv = sigma_deriv
)

sm <- stan_model("../Stan/extrapolate.stan")

fit <- sampling(sm,
                data = stan_data,
                iter = 2000,
                warmup = 1000,
                chains = 4,
                control = list(adapt_delta = 0.95, max_treedepth = 12),
                pars = c("L_asym", "sf", "ell", "sigma_y"))

print(fit, probs = c(0.05, 0.5, 0.95))

ggplot(mit_boot,aes(x=n,y=Delta)) + geom_point()
ggplot(pilot,aes(x=n,y=Delta)) + geom_point()