#' Utility function to fit extrapolation model, for use with [conduct_interpolation()]
#'
#' @param formula A formula object giving the model to be fit.
#' @param data A data frame giving the data the model is to be fit on.
#' @param N_grid An integer vector of N values to conduct interpolation and extrapolation on.
#' @param maxN_obs A positive integer giving value of the largest N in the observed data.
#' @param start Optional named vector of starting values for the model parameters.
#' @param lower Optional named vector of lower bounds for the model parameters.
#' @param upper Optional named vector of upper bounds for the model parameters.
#' @return A fitted model object of the chosen type.
#' @importFrom minpack.lm nlsLM
#' @importFrom stats predict approxfun

fit_and_predict <- function(formula, data, N_grid, maxN_obs, start = NULL, lower = NULL, upper = NULL) {
  tryCatch({
    args <- list(formula = formula, data = data)

    if (!is.null(start)) args$start <- start
    if (!is.null(lower)) args$lower <- lower
    if (!is.null(upper)) args$upper <- upper

    fit <- do.call(nlsLM, args)

    pred_extrap <- predict(fit, newdata = data.frame(n = N_grid[N_grid > maxN_obs]))
    full_fit <- c(
      approxfun(data$n, fitted(fit), method = "linear", rule = 2)(N_grid[N_grid <= maxN_obs]),
      pred_extrap
    )
    return(full_fit)
  }, error = function(e) {
    return(rep(NA, length(N_grid)))
  })
}


#' Wrapper function for fitting extrapolation model to a single object of class "scb_data". Allows fitting of custom functions, but for general use interpolate_scb should be used instead.
#'
#' @param scbobject An object of class "scb_data" for interpolation to be conducted on.
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param maxN A positive integer giving value of the largest N in the observed data.
#' @param delta_formula Formula of the form Delta ~ model(n,...) giving the NLS model to be applied to the delta curve.
#' @param epsilon_formula Formula of the form Epsilon ~ model(n,...) giving the NLS model to be applied to the epsilon curve.
#' @param delta_lower_bounds Optional named vector of upper bounds for the model parameters.
#' @param epsilon_lower_bounds A formula object giving the model to be fit.
#' @param delta_upper_bounds A data frame giving the data the model is to be fit on.
#' @param epsilon_upper_bounds An integer vector of N values to conduct interpolation and extrapolation on.
#' @param delta_start A positive integer giving value of the largest N in the observed data.
#' @param epsilon_start Optional named vector of starting values for the model parameters.
#' @return A named `list` containing the interpolated dataframe, the original input data frame, and the given values of epsilon, delta, and maxN.
#' @seealso [interpolate_scb()] is the main wrapper for interpolation on a list.
#' @export
conduct_interpolation <- function(
    scbobject,
    epsilon,
    delta,
    maxN,
    delta_formula,
    epsilon_formula,
    delta_lower_bounds=NULL,
    epsilon_lower_bounds=NULL,
    delta_upper_bounds=NULL,
    epsilon_upper_bounds=NULL,
    delta_start = NULL,
    epsilon_start = NULL){
  dat <- getpac(scbobject,epsilon,delta)$Summary
  train <- dat
  nvals = unique(dat$n)
  maxN_obs <- max(nvals)
  N_grid = seq(min(nvals),maxN,by=1)
  full_fit_delta <- fit_and_predict(
    formula = delta_formula,
    data = train,
    start = delta_start,
    lower = delta_lower_bounds,
    upper = delta_upper_bounds,
    N_grid = N_grid,
    maxN_obs = maxN_obs
  )
  full_fit_epsilon <- fit_and_predict(
    formula = epsilon_formula,
    data = train,
    start = epsilon_start,
    lower = epsilon_lower_bounds,
    upper = epsilon_upper_bounds,
    N_grid = N_grid,
    maxN_obs = maxN_obs
  )
  interp_plot_df <- data.frame(
    N = N_grid,
    Fit_Delta = full_fit_delta,
    Fit_Epsilon = full_fit_epsilon,
    Type = ifelse(N_grid > maxN_obs, "Extrapolated", "Observed")
  )
  out <- list(interpolated = interp_plot_df,original = dat,epsilon = epsilon,delta = delta,maxN = maxN_obs)
  class(out) <- "empirical_scb"
  return(out)
}

#' Wrapper function for fitting extrapolation model to a single object of class "scb_data". Allows fitting of custom functions, but for general use interpolate_scb should be used instead.
#'
#' @param scbobject An object of class "scb_data" for interpolation to be conducted on.
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param maxN A positive integer giving value of the largest N in the observed data.
#' @param delta_formula Formula of the form Delta ~ model(n,...) giving the NLS model to be applied to the delta curve.
#' @param epsilon_formula Formula of the form Epsilon ~ model(n,...) giving the NLS model to be applied to the epsilon curve.
#' @param delta_lower_bounds Optional named vector of upper bounds for the model parameters.
#' @param epsilon_lower_bounds A formula object giving the model to be fit.
#' @param delta_upper_bounds A data frame giving the data the model is to be fit on.
#' @param epsilon_upper_bounds An integer vector of N values to conduct interpolation and extrapolation on.
#' @param delta_start A positive integer giving value of the largest N in the observed data.
#' @param epsilon_start Optional named vector of starting values for the model parameters.
#' @return A named `list` containing the interpolated dataframe, the original input data frame, and the given values of epsilon, delta, and maxN.
#' @seealso [interpolate_scb()] is the main wrapper for interpolation on a list.
#' @export
interpolate_scb <- function(data_list,
                             delta_interp_fun = c("logis","logis5","logis4","declin"),
                             epsilon_interp_fun = c("gompertz","exp_plateau","weibull","quad_plateau"),
                             epsilon,
                             delta,
                             maxN,
                             delta_lower_bounds=NULL,
                             epsilon_lower_bounds=NULL,
                             delta_upper_bounds=NULL,
                             epsilon_upper_bounds=NULL,
                             delta_start = NULL,
                             epsilon_start = NULL
){
  delta_formula <- switch(match.arg(delta_interp_fun),
                          logis = Delta ~ SSlogis(n, Asym, xmid, scal),
                          logis5 = Delta ~ SSlogis5(n, asym1, asym2, xmid, iscal, theta),
                          logis4 = Delta ~ SSfpl(n, A, B, xmid, scal),
                          declin = Delta ~ SSdlf(n,asym, a2, xmid, scal)
                          )
  epsilon_formula <- switch(match.arg(epsilon_interp_fun),
                          gompertz = Epsilon ~ SSgompertz(n, Asym, b2, b3),
                          exp_plateau = Epsilon ~ SSexpfp(n, a, c, xs),
                          weibull = Epsilon ~ SSweibull(n, Asym, Drop, lrc, pwr),
                          quad_plateau = Epsilon ~ SSquadp3(n, a, b, c)
  )
  out <- pblapply(data_list,conduct_interpolation,
                  epsilon=epsilon,
                  delta=delta,
                  maxN=maxN,
                  delta_formula = delta_formula,
                  epsilon_formula = epsilon_formula,
                  delta_lower_bounds=delta_lower_bounds,
                  epsilon_lower_bounds=epsilon_lower_bounds,
                  delta_upper_bounds=delta_upper_bounds,
                  epsilon_upper_bounds=epsilon_upper_bounds,
                  delta_start = delta_start,
                  epsilon_start = epsilon_start)
  class(out) <- "empirical_scb_list"
  return(out)
}

#' Plot method for an `empirical_scb_list` object
#'
#' Visualizes bootstrap-estimated empirical sample complexity bounds (SCB) for either delta or epsilon.
#'
#' @param x An object of class `"empirical_scb_list"` containing extrapolated SCB values.
#' @param truedata A bootstrapped list of benchmark simulations, each of class `"scb_data"`.
#' @param alpha Numeric between 0 and 1. Significance level used to compute bootstrap confidence intervals.
#' @param initial Integer. The size of the initial training subsample used in estimation.
#' @param plot_type Character string. Determines which SCB to plot: `"Delta"` (default) or `"Epsilon"`.
#' @param include_legend Logical. Whether to display a legend (default: `TRUE`).
#' @param include_title Logical. Whether to include a title (default: `TRUE`).
#' @return A \link[ggplot2]{ggplot} object displaying either the SCB-Delta or SCB-Epsilon curve with bootstrap confidence bands.
#' @seealso [interpolate_scb()] in order to prepare input data..
#' @export
#' @import dplyr
#' @import ggplot2
#' @method plot empirical_scb_list

plot.empirical_scb_list <- function(x,
                               truedata,
                               alpha,
                               initial,
                               plot_type = c("Delta","Epsilon"),
                               include_legend=TRUE,
                               include_title=TRUE){
  plot_type <- match.arg(plot_type)
  epsilon <- x[[1]]$epsilon
  delta <- x[[1]]$delta
  maxN <- x[[1]]$maxN

  interp <- bind_rows(lapply(x, function(y) y$interpolated)) %>%
    group_by(N) %>%
    summarise(
      CI_Upper_Delta = quantile(Fit_Delta, 1 - alpha / 2, na.rm = TRUE),
      CI_Upper_Epsilon = quantile(Fit_Epsilon, 1 - alpha / 2, na.rm = TRUE),
      CI_Lower_Delta = quantile(Fit_Delta, alpha / 2, na.rm = TRUE),
      CI_Lower_Epsilon = quantile(Fit_Epsilon, alpha / 2, na.rm = TRUE),
      Fit_Delta = mean(Fit_Delta, na.rm = TRUE),
      Fit_Epsilon = mean(Fit_Epsilon, na.rm = TRUE),
      Type = first(Type)
    )

  truedata <- bind_rows(lapply(truedata, function(x) getpac(x, epsilon, delta)$Summary)) %>%
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

  # Determine SCB values
  scb_value_delta <- min(interp$N[interp$Fit_Delta <= delta], na.rm = TRUE)
  scb_status_delta <- ifelse(scb_value_delta <= maxN, "Observed", "Extrapolated")

  scb_value_epsilon <- min(interp$N[interp$Fit_Epsilon <= epsilon], na.rm = TRUE)
  scb_status_epsilon <- ifelse(scb_value_epsilon <= maxN, "Observed", "Extrapolated")

  base_theme <- theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

  add_labels <- function(p, title_text, y_label) {
    if (include_legend & include_title) {
      p + labs(
        title = title_text, x = "Training Sample Size (N)", y = y_label,
        color = "Segment Type", linetype = "Segment Type"
      ) +
        base_theme +
        theme(
          legend.box.background = element_rect(color = "black", size = 0.4),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5)
        )
    } else if (include_legend) {
      p + labs(
        x = "Training Sample Size (N)", y = y_label,
        color = "Segment Type", linetype = "Segment Type"
      ) +
        base_theme +
        theme(legend.box.background = element_rect(color = "black", size = 0.4),
              legend.position = "bottom")
    } else if (include_title) {
      p + labs(
        title = title_text, x = "Training Sample Size (N)", y = y_label
      ) +
        base_theme +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", hjust = 0.5))
    } else {
      p + labs(x = "Training Sample Size (N)", y = y_label) +
        base_theme +
        theme(legend.position = "none")
    }
  }

  if (plot_type == "Delta") {
    p <- ggplot(interp, aes(x = N)) +
      geom_ribbon(aes(ymin = CI_Lower_Delta, ymax = CI_Upper_Delta), fill = "#0072B2", alpha = 0.5) +
      geom_line(aes(y = Fit_Delta, color = Type), size = 1.2) +
      geom_point(data = truedata, aes(x = N, y = Delta), size = 2, color = "black", shape = 3) +
      geom_ribbon(data = truedata, aes(x = N, ymin = CI_Lower_Delta, ymax = CI_Upper_Delta), fill = "#D55E00", alpha = 0.3) +
      geom_hline(yintercept = delta, linetype = "dashed", color = "red") +
      geom_vline(xintercept = scb_value_delta, linetype = "dotted", color = "darkgreen") +
      geom_vline(xintercept = initial, linetype = "solid", color = "black") +
      annotate("text", x = scb_value_delta - 100, y = 0.1,
               label = paste("SCB ≈", scb_value_delta, "\n(", scb_status_delta, ")"),
               color = "darkgreen", size = 4.5, hjust = 0) +
      annotate("text", x = initial + 100, y = 0.5,
               label = paste("Initial Subsample =", initial),
               color = "black", size = 4.5, hjust = 0) +
      scale_color_manual(values = c("Observed" = "blue", "Extrapolated" = "darkorange"))

    return(add_labels(p, "Interpolated SCB (δ) Curve with Bootstrap Confidence Bands", "P(error ≤ ε)"))

  } else {  # Epsilon
    p <- ggplot(interp, aes(x = N)) +
      geom_ribbon(aes(ymin = CI_Lower_Epsilon, ymax = CI_Upper_Epsilon), fill = "#0072B2", alpha = 0.5) +
      geom_line(aes(y = Fit_Epsilon, color = Type, linetype = Type), size = 1.2) +
      geom_point(data = truedata, aes(x = N, y = Epsilon), size = 2, color = "black", shape = 3) +
      geom_ribbon(data = truedata, aes(x = N, ymin = CI_Lower_Epsilon, ymax = CI_Upper_Epsilon), fill = "#D55E00", alpha = 0.3) +
      geom_hline(yintercept = epsilon, linetype = "dashed", color = "red") +
      geom_vline(xintercept = scb_value_epsilon, linetype = "dotted", color = "darkgreen") +
      geom_vline(xintercept = initial, linetype = "solid", color = "black") +
      annotate("text", x = scb_value_epsilon - 100, y = 0.1,
               label = paste("SCB ≈", scb_value_epsilon, "\n(", scb_status_epsilon, ")"),
               color = "darkgreen", size = 4.5, hjust = 0) +
      annotate("text", x = initial + 100, y = 0.5,
               label = paste("Initial Subsample =", initial),
               color = "black", size = 4.5, hjust = 0) +
      scale_color_manual(values = c("Observed" = "blue", "Extrapolated" = "darkorange"))

    return(add_labels(p, "Interpolated SCB (ε) Curve with Bootstrap Confidence Band", "ε | δ"))
  }
}


#' Summary of empirical sample complexity bound results
#'
#' For an \code{empirical_scb_list} object, finds
#' 1. the “mean‐fit” SCB crossing \(N\) (where the average bootstrap curve first drops below the target),
#' 2. a lower bound on that crossing (the smallest \(N\) at which the *lower* CI envelope crosses), and
#' 3. an upper bound on the crossing (the smallest \(N\) at which the *upper* CI envelope crosses).
#'
#' @param object An object of class \code{"empirical_scb_list"} containing bootstrap replicates.
#' @param truedata (ignored) Kept for signature consistency.
#' @param alpha Numeric in \((0,1)\).  Two‑sided CI level for the envelope (e.g.\ 0.05 for 95%).
#' @param initial Integer.  The size of the initial subsample used.
#' @param ... Additional args (ignored).
#'
#' @return Invisibly, a list of components
#' \describe{
#'   \item{\code{delta}}{List with
#'     \code{SCB_N} (mean‐curve crossing),
#'     \code{status} (“Observed”/“Extrapolated”),
#'     \code{CI_lower_N}, \code{CI_upper_N} (bounds on the crossing).
#'   }
#'   \item{\code{epsilon}}{Same four elements for the \(\epsilon\) target.}
#'   \item{\code{alpha}}{The CI level.}
#'   \item{\code{initial}}{The initial subsample size.}
#' }
#' @seealso \code{\link{plot.empirical_scb_list}}, \code{\link{getpac}}
#' @export
#' @method summary empirical_scb_list
summary.empirical_scb_list <- function(object, truedata, alpha, initial, ...) {
  # pull targets and maxN
  delta   <- object[[1]]$delta
  epsilon <- object[[1]]$epsilon
  maxN    <- object[[1]]$maxN

  # build envelope + mean curve
  env <- dplyr::bind_rows(lapply(object, `[[`, "interpolated")) %>%
    dplyr::group_by(N) %>%
    dplyr::summarise(
      CI_Lower_Delta   = quantile(Fit_Delta,   alpha/2,    na.rm=TRUE),
      CI_Upper_Delta   = quantile(Fit_Delta,   1-alpha/2,  na.rm=TRUE),
      CI_Lower_Epsilon = quantile(Fit_Epsilon, alpha/2,    na.rm=TRUE),
      CI_Upper_Epsilon = quantile(Fit_Epsilon, 1-alpha/2,  na.rm=TRUE),
      Fit_Delta        = mean(Fit_Delta,      na.rm=TRUE),
      Fit_Epsilon      = mean(Fit_Epsilon,    na.rm=TRUE)
    )

  scb_mean_delta   <- min(env$N[env$Fit_Delta   <= delta],   na.rm=TRUE)
  scb_mean_epsilon <- min(env$N[env$Fit_Epsilon <= epsilon], na.rm=TRUE)

  delta_lowN  <- min(env$N[env$CI_Lower_Delta   <= delta],   na.rm=TRUE)
  delta_highN <- min(env$N[env$CI_Upper_Delta   <= delta],   na.rm=TRUE)
  eps_lowN    <- min(env$N[env$CI_Lower_Epsilon <= epsilon], na.rm=TRUE)
  eps_highN   <- min(env$N[env$CI_Upper_Epsilon <= epsilon], na.rm=TRUE)

  status_delta   <- if (scb_mean_delta   <= maxN) "Observed" else "Extrapolated"
  status_epsilon <- if (scb_mean_epsilon <= maxN) "Observed" else "Extrapolated"

  out <- list(
    delta   = list(SCB_N = scb_mean_delta,
                   status = status_delta,
                   CI_lower_N = delta_lowN,
                   CI_upper_N = delta_highN),
    epsilon = list(SCB_N = scb_mean_epsilon,
                   status = status_epsilon,
                   CI_lower_N = eps_lowN,
                   CI_upper_N = eps_highN),
    alpha   = alpha,
    initial = initial
  )
  class(out) <- "summary_empirical_scb_list"
  invisible(out)
}

#' @export
#' @method print summary_empirical_scb_list
print.summary_empirical_scb_list <- function(x, ...) {
  cat(sprintf(
    "Empirical SCB Summary  (alpha = %g,  initial = %d)\n\n",
    x$alpha, x$initial
  ))

  cat("Delta SCB:\n")
  cat(sprintf("  SCB N (mean‐curve): %d  (%s)\n",
              x$delta$SCB_N, x$delta$status))
  cat(sprintf("  %g%% CI on N crossing: [%d, %d]\n\n",
              (1 - x$alpha) * 100,
              x$delta$CI_lower_N, x$delta$CI_upper_N))

  cat("Epsilon SCB:\n")
  cat(sprintf("  SCB N (mean‐curve): %d  (%s)\n",
              x$epsilon$SCB_N, x$epsilon$status))
  cat(sprintf("  %g%% CI on N crossing: [%d, %d]\n",
              (1 - x$alpha) * 100,
              x$epsilon$CI_lower_N, x$epsilon$CI_upper_N))

  invisible(x)
}



#' Utility function for creating custom classification function for use in SCB calculations.
#'
#' @param model_fun A binary classification model supplied by the user, e.g. glm
#' @param extra_args A list of additional default arguments to be passed to the classification model, e.g. family = binomial(link = "logit")
#' @param split Logical indicating whether the model expects a single data argument or separate x/y values.
#' @param arg_map Named list giving mappings for names of formula, data, x, and y arguments expected in `model_fun`. Must include formula and data if split is FALSE or x and y otherwise.
#' @return A function taking supplied arguments that can be passed to other package functions such as [estimate_accuracy()]
#' @export

create_scb_model <- function(model_fun,
                             extra_args = list(),
                             split = FALSE,
                             arg_map = list(formula = "formula",
                                            data    = "data",
                                            x       = "x",
                                            y       = "y")) {
  if(split){
    if (is.null(arg_map$x) || is.null(arg_map$y)) {
      stop("When split = TRUE, 'arg_map' must include mappings for 'x' and 'y'.")
    }
    function(formula, data,...) {
      if (missing(formula) || missing(data)) {
        stop("Both 'formula' and 'data' must be supplied to the wrapper.")
      }
      x <- data[,-ncol(data)]
      y <- as.vector(data[,ncol(data)])

      wrapper_args <- list(...)
      mapped <- list()
      mapped[[ arg_map$x ]] <- X
      mapped[[ arg_map$y ]] <- y
      all_args <- c(mapped,
                    extra_args,
                    wrapper_args[!names(wrapper_args) %in% names(arg_map)])
      fit <- do.call(model_fun, all_args)
      class(fit) <- c("svrclass", class(fit))
      return(fit)
    }
  } else{
    if (is.null(arg_map$formula) || is.null(arg_map$data)) {
      stop("When split = FALSE, 'arg_map' must include mappings for 'formula' and 'data'.")
    }
    function(formula, data,...) {
      if (missing(formula) || missing(data)) {
        stop("Both 'formula' and 'data' must be supplied to the wrapper.")
      }
      wrapper_args <- list(...)
      mapped <- list()
      mapped[[ arg_map$formula ]] <- formula
      mapped[[ arg_map$data    ]] <- data
      all_args <- c(mapped,
                    extra_args,
                    wrapper_args[!names(wrapper_args) %in% names(arg_map)])
      fit <- do.call(model_fun, all_args)
      class(fit) <- c("svrclass", class(fit))
      return(fit)
    }
  }
}


#' Utility function for creating custom prediction function for use in SCB calculations.
#'
#' @param predict_fun A binary prediction model supplied by the user, e.g. predict.glm
#' @param extra_args A list of additional default arguments to be passed to the classification model, e.g. type="response"
#' @param transform_fn Function giving the transformation, if any, to be applied to the prediction output to generate binary predictions. Defaults to identity()
#' @param arg_map Named list giving mappings for names of `m` and `newdata` arguments expected in `predict_fun`.
#' @return A prediction function taking supplied arguments that can be passed to other package functions such as [estimate_accuracy()]
#' @export

create_scb_prediction <- function(predict_fun,
                             extra_args = list(),
                             transform_fn = identity,
                             arg_map = list(m = "m",
                                            newdata = "newdata")) {
    if (is.null(arg_map$m) || is.null(arg_map$newdata)) {
      stop("When split = TRUE, 'arg_map' must include mappings for 'm' and 'newdata'.")
    }
    function(m, newdata,...) {
      if (missing(m) || missing(newdata)) {
        stop("Both 'm' and 'newdata' must be supplied to the wrapper.")
      }
      wrapper_args <- list(...)
      mapped <- list()
      mapped[[ arg_map$m ]] <- m
      mapped[[ arg_map$newdata ]] <- newdata
      all_args <- c(mapped,
                    extra_args,
                    wrapper_args[!names(wrapper_args) %in% names(arg_map)])
      out <- do.call(predict_fun,all_args)
      out <- factor(transform_fn(out), levels = c("0", "1")) #Important - must specify levels to account for possibility of all observations being classified into the same class in smaller samples
      return(out)
    }
  }
}


#' Utility function to generate data points for estimation of the VC Dimension of a user-specified binary classification algorithm given a specified sample size.
#'
#' @param x An integer giving the desired sample size for which the target function is to be approximated.
#' @param l A positive integer giving dimension (number of input features) of the model.
#' @param m A positive integer giving the number of simulations to be performed at each design point (sample size value). Higher values give more accurate results but increase computation time.
#' @param model A binary classification model supplied by the user. Must take arguments `formula` and `data`
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param sparse Logical indicating whether sparse matrix generation should be used to save on memory. Defaults to false for better accuracy.
#' @param density Real number between 0 and 1 giving the proportion of non 0 entries in the sparse matrix. Used only if sparse is TRUE.
#' @param ... Additional model parameters to be specified by the user.
#' @importFrom stats rnorm predict
#' @importFrom Matrix rsparsematrix
#' @return A real number giving the estimated value of Xi(n), the bounding function


risk_bounds <- function(x, l, m, model,predictfn, sparse,density=NULL,...) {
  n <- x
  if(!is.null(predictfn)){
    predict.svrclass <- predictfn
  }
  xis <- numeric(m)  # Preallocate memory
  for (j in 1:m) {
    for (attempt in 1:100) {  # Limit attempts to avoid infinite loop
      if(sparse){

        x <- rsparsematrix(nrow = 2 * n, ncol = l,density=density)  # Use sparse matrix to save memory
      } else{

        x <- matrix(rnorm(2 * n * l), nrow = 2 * n, ncol = l)
      }
      coeff <- rnorm(l)
      y <- as.numeric(((x %*% coeff) > 0))

      W_idx <- sample(1:(2 * n), n)
      W <- x[W_idx, ]
      Wprime <- x[-W_idx, ]
      Wprime_y <- 1 - y[-W_idx]  # Flip labels

      y_train <- c(y[W_idx], Wprime_y)
      x_train <- rbind(W, Wprime)

      # Train model
      fhat <- tryCatch(
        model(formula = y ~ ., data = cbind(x_train,y_train), ...
              ),
        error = function(e) NULL
      )
      if (!is.null(fhat)) break  # Exit loop if model trains successfully

      gc()
      if (attempt == 100) stop("Failed to train model after 100 attempts.")
    }

    # Predict and compute RW and RWprime
    RW <- suppressWarnings(mean(predict(fhat, W) != factor(y[W_idx])))
    RWprime <- suppressWarnings(mean(predict(fhat, Wprime) != factor(1 - Wprime_y)))
    xis[j] <- abs(RW - RWprime)
    #gc()
  }
  return(mean(xis))
}




#' Utility function to define the least-squares loss function to be optimized for [simvcd()]
#'
#' @param h A positive real number giving the current guess at VC dimension
#' @param ngrid Vector of sample sizes for which the bounding function is estimated.
#' @param xi Vector of estimated values of the bounding function, usually obtained from [risk_bounds()]
#' @param a Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param a1 Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param a11 Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @return A real number giving the estimated value of the MSE given the current guess.
#' @seealso [simvcd()], the user-facing function for simulating VC dimension and [risk_bounds()] to generate estimates for xi.
loss <- function(h,ngrid,xi,a=0.16,a1=1.2,a11=0.14927){
  #These constants are calculated in Vapnik, Levin and Le Cun 1994
  #based on the known VC dimension of linear discriminant functions
  #and rely on the assumption that they are universal for all classifiers
  ratio <- ngrid/h
  phi <- a*((log((2*ratio))+1)/(ratio-a11))*(1+sqrt((1+((a1*(ratio - a11))/(1+log(2*ratio))))))
  test <- ngrid < (h/2)
  phihats <- ifelse(test,1,phi)
  devs <- (xi - phihats)^2
  out <- sum(devs)
  return(out)
  #TODO - add gradient (very messy function but straightforward to derive)
}

#' Estimate the Vapnik-Chervonenkis (VC) dimension of an arbitrary binary classification algorithm.
#'
#' @param model A binary classification model supplied by the user. Must take arguments `formula` and `data`
#' @param dim A positive integer giving dimension (number of input features) of the model.
#' @param maxn Gives the vertical dimension of the data (number of observations) to be generated.
#' @param m A positive integer giving the number of simulations to be performed at each design point (sample size value). Higher values give more accurate results but increase computation time.
#' @param k A positive integer giving the number of design points (sample size values) for which the bounding function is to be estimated. Higher values give more accurate results but increase computation time.
#' @param parallel Boolean indicating whether or not to use parallel processing.
#' @param coreoffset If `parallel` is true, a positive integer indicating the number of free threads to be kept unused. Should not be larger than the number of CPU cores.
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param a Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param a1 Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param a11 Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param minn Optional argument to set a different minimum n than the dimension of the algorithm. Useful with e.g. regularized regression models such as elastic net.
#' @param sparse Logical indicating whether sparse matrix generation should be used to save on memory. Defaults to false for better accuracy.
#' @param density Real number between 0 and 1 giving the proportion of non 0 entries in the sparse matrix. Used only if sparse is TRUE.
#' @param backend One of the parallel backends used by [future::plan()]. See function documentation for more details.
#' @param packages A `list` of strings giving the names of packages to be loaded in order to estimate the model.
#' @param ... Additional arguments that need to be passed to `model`
#' @return A real number giving the estimated value of the VC dimension of the supplied model.
#' @seealso [scb()], to calculate sample complexity bounds given estimated VCD.
#' @importFrom parallel detectCores
#' @importFrom future plan cluster
#' @importFrom furrr future_map_dbl furrr_options
#' @import dplyr
#' @importFrom stats optim
#' @importFrom progressr progressor
#' @examples
#' mylogit <- function(formula, data){
#' m <- structure(
#'   glm(formula=formula,data=data,family=binomial(link="logit")),
#'   class=c("svrclass","glm")  #IMPORTANT - must use the class svrclass to work correctly
#' )
#' return(m)
#' }
#' mypred <- function(m,newdata){
#' out <- predict.glm(m,newdata,type="response")
#' out <- factor(ifelse(out>0.5,1,0),levels=c("0","1"))
#' #Important - must specify levels to account for possibility of all
#' #observations being classified into the same class in smaller samples
#' return(out)
#' }
#' library(parallel)
#' vcd <- simvcd(model=mylogit,dim=7,m=10,k=10,maxn=50,predictfn = mypred,
#'     coreoffset = (detectCores() -2))
#' @export
simvcd <- function(model,dim,m=1000,k=1000,maxn=5000,parallel = TRUE,coreoffset=0, predictfn = NULL, a=0.16,a1=1.2,a11=0.14927,minn = (dim+1),sparse=FALSE,density=NULL,backend = c("multisession","multicore","cluster","sequential"),packages=list(), ...){
  force(minn)
  ngrid <- round(seq(minn,maxn,(maxn/k)),0)
  backend <- match.arg(backend)
  if(parallel){
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      cl <- 2L
    } else {
      # use all cores in devtools::test()
      cl <- detectCores() -coreoffset
    }
  } else{
    cl <- 1
  }

  l<-dim
  plan(get(backend), workers = cl)  # Use chosen backend
  p <- progressor(steps = length(ngrid))
  temp <- function(x,l,m,model,packages,predictfn,sparse,density,...){
    p()
    #set.seed(as.numeric(Sys.time()))
    lapply(packages, library, character.only = TRUE)
    r <- risk_bounds(x=x,l=l,m=m,model=model,predictfn=predictfn,sparse=sparse,density=density,...)
    return(r)
  }
  xihats <- future_map_dbl(ngrid, temp,l=l,m=m,model=model,packages=packages,predictfn=predictfn,sparse=sparse,density=density,...,.options = furrr_options(seed = TRUE))
  vcd <- optim((l+1),loss,ngrid=ngrid,xi=xihats, a=a,a1=a1,a11=a11, method="Brent",lower=1,upper = 2*(max(ngrid)),...)
  return(vcd$par)
}

#' Calculate sample complexity bounds for a classifier given target accuracy
#'
#' @param vcd The Vapnik-Chervonenkis dimension (VCD) of the chosen classifier. If `theor` is `FALSE`, this can be left unspecified and [simvcd()] will be called to estimate the VCD
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param eta A real number between 0 and 1 giving the probability of misclassification error in the training data.
#' @param theor A Boolean indicating whether the theoretical VCD is to be used. If `FALSE`, it will instead be estimated using [simvcd()]
#' @param ... Arguments to be passed to [simvcd()]
#' @return A real number giving the sample complexity bound for the specified parameters.
#' @seealso [simvcd()], to calculate VCD for a chosen model
#' @examples
#' mylogit <- function(formula, data){
#' m <- structure(
#'   glm(formula=formula,data=data,family=binomial(link="logit")),
#'   class=c("svrclass","glm")  #IMPORTANT - must use the class svrclass to work correctly
#' )
#' return(m)
#' }
#' mypred <- function(m,newdata){
#' out <- predict.glm(m,newdata,type="response")
#' out <- factor(ifelse(out>0.5,1,0),levels=c("0","1"))
#' #Important - must specify levels to account for possibility of all
#' #observations being classified into the same class in smaller samples
#' return(out)
#' }
#' library(parallel)
#' scb(epsilon=0.05,delta=0.05,eta=0.05,theor=FALSE,
#' model=mylogit,dim=7,m=10,k=10,maxn=50,predictfn = mypred,
#'     coreoffset = (detectCores() -2))
#' vcd <- 7
#' scb(vcd,epsilon=0.05,delta=0.05,eta=0.05)
#' @export
scb <- function(vcd=NULL,epsilon=NULL,delta=NULL,eta=NULL,theor=TRUE,...){
  if(!theor){
    vc <- simvcd(...) #Pass extra args to the simvcd function if theoretical value is unknown
  } else if(is.null(vcd)){
    simpleError("No VCD specified")
  } else {
    vc <- vcd
  }
  m <- (1/((epsilon)*((1-2*eta)^2)))*(vc + log(1/delta))
  return(ceiling(m))
}

