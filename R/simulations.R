#' Utility function to generate accuracy metrics, for use with [estimate_accuracy()]
#'
#' @param n An integer giving the desired sample size for which the target function is to be calculated.
#' @param method An optional string stating the distribution from which data is to be generated. Default is i.i.d. uniform sampling. Currently also supports "Class Imbalance". Can also take a function outputting a vector of probabilities if the user wishes to specify a custom distribution.
#' @param p If method is 'Class Imbalance', gives the degree of weight placed on the positive class.
#' @param dat A rectangular `data.frame` or matrix-like object giving the full data from which samples are to be drawn. If left unspecified, [gendata()] is called to produce synthetic data with an appropriate structure.
#' @param model A function giving the model to be estimated
#' @param eta A real number between 0 and 1 giving the probability of misclassification error in the training data.
#' @param nsample A positive integer giving the number of samples to be generated for each value of $n$. Larger values give more accurate results.
#' @param outcome A string giving the name of the outcome variable.
#' @param power A logical indicating whether experimental power based on the predictions should also be reported
#' @param effect_size If `power` is `TRUE`, a real number indicating the scaled effect size the user would like to be able to detect.
#' @param powersims If `power` is `TRUE`, an integer indicating the number of simulations to be conducted at each step to calculate power.
#' @param alpha If `power` is `TRUE`, a real number between 0 and 1 indicating the probability of Type I error to be used for hypothesis testing. Default is 0.05.
#' @param split A logical indicating whether the data was passed as a single data frame or separately.
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param replacement A logical flag indicating whether sampling should be performed with replacement.
#' @param ... Additional model parameters to be specified by the user.
#' @return A data frame giving performance metrics for the specified sample size.
#' @importFrom caret precision recall F_meas
#' @importFrom stats rbinom predict rnorm lm
#' @import dplyr
#' @export
acc_sim <- function(n, method, p, dat, model, eta, nsample, outcome, power, effect_size, powersims, alpha, split, predictfn, replacement,...) {
  if (!is.null(predictfn)) {
    predict.svrclass <- predictfn
  }

  accuracy <- numeric(nsample)
  prec <- numeric(nsample)
  rec <- numeric(nsample)
  fscore <- numeric(nsample)
  pwr <- if (power) numeric(nsample) else NULL

  oc <- if (split) as.vector(dat[, ncol(dat)]) else dat[[outcome]]

  for (j in seq_len(nsample)) {
    for (attempt in 1:100) {
      if (method == "Uniform") {
        indices <- sample(seq_len(nrow(dat)), n,replace = replacement)
        samp <- dat[indices, , drop = FALSE]
        error <- rbinom(n, 1, eta)
      } else if (method == "Class Imbalance") {
        if (is.null(p)) stop("Class Imbalance method selected. Please supply a class imbalance parameter.")
        probs <- ifelse(oc == 1, p, 1 - p)
        indices <- sample(seq_len(nrow(dat)), n, prob = probs,replace = replacement)
        samp <- dat[indices, , drop = FALSE]
        error <- rbinom(n, 1, eta)
      } else {
        probs <- method(dat)
        indices <- sample(seq_len(nrow(dat)), n, prob = probs,replace = replacement)
        samp <- dat[indices, , drop = FALSE]
        error <- rbinom(n, 1, eta)
      }

      if (!split) {
        soc <- samp[[outcome]]
        if(is.factor(soc)){
          samp$outobs <- factor(ifelse(error, !as.numeric(as.character(soc)), as.numeric(as.character(soc))), levels = c("0", "1"))
        } else{
          samp$outobs <- factor(ifelse(error, (1-soc), soc), levels = c("0", "1"))
        }
        samp <- samp[, !(names(samp) %in% outcome), drop = FALSE]
      } else {
        soc <- as.vector(samp[, ncol(samp)])
        samp <- cbind(samp[, -ncol(samp), drop = FALSE], as.numeric(ifelse(error,!soc,soc)))
      }

      m <- tryCatch(model(outobs ~ ., data = samp, ...),
                    error = function(e) {
                      warning("Failed to train model, regenerating data, see below error:")
                      warning(e)
                      return(NULL)
                      }
      )

      if (!is.null(m)) break
      gc()
      if(attempt == 100) stop("Failed to train model after 100 attempts.")
    }

    foc <- factor(oc, levels = c("0", "1"))
    if(split){
      pred <- suppressWarnings(predict(m, dat[, -ncol(dat), drop = FALSE]))
    } else{
      pred <- suppressWarnings(predict(m, dat[,-which(names(dat)==outcome), drop = FALSE]))
    }

    accuracy[j] <- mean(pred == foc)
    prec[j] <- tryCatch(precision(table(pred, foc), relevant = 1), error = function(e) NA)
    rec[j] <- tryCatch(recall(table(pred, foc), relevant = 1), error = function(e) NA)
    fscore[j] <- tryCatch(F_meas(table(pred, foc), relevant = 1), error = function(e) NA)

    if (power) {
      Dobs <- as.numeric(levels(pred)[pred])
      Dtrue <- if(is.factor(oc)) as.numeric(levels(oc)[oc]) else oc
      reject <- replicate(powersims, {
        Y <- effect_size * Dtrue + rnorm(length(Dobs))
        mdl <- lm(Y ~ Dobs)
        tryCatch(summary(mdl)$coefficients[2, 4] < alpha,
                 error = function(e) {
                   warning("Failed to generate experimental data, model might have predicted all observations of a single class")
                   return(NA)
                 }
        )
      })
      pwr[j] <- mean(reject, na.rm = TRUE)
    }
  }

  data.frame(accuracy, prec, rec, fscore, n, pwr = if (power) pwr else NA)
}

#' Estimate sample complexity bounds for a binary classification algorithm using either simulated or user-supplied data.
#'
#' @param formula A `formula` that can be passed to the `model` argument to define the classification algorithm
#' @param model A binary classification model supplied by the user. Must take arguments `formula` and `data`
#' @param data Optional. A rectangular `data.frame` object giving the full data from which samples are to be drawn. If left unspecified, [gendata()] is called to produce synthetic data with an appropriate structure.
#' @param dim Required if `data` is unspecified. Gives the horizontal dimension of the data (number of predictor variables) to be generated.
#' @param maxn Required if `data` is unspecified. Gives the vertical dimension of the data (number of observations) to be generated.
#' @param sparse Optional. A logical giving whether to generate sparse data, if data was not given.
#' @param density Real number between 0 and 1 giving the proportion of non 0 entries in the sparse matrix. Used only if sparse is TRUE.
#' @param upperlimit Optional. A positive integer giving the maximum sample size to be simulated, if data was supplied.
#' @param nsample A positive integer giving the number of samples to be generated for each value of $n$. Larger values give more accurate results.
#' @param steps A positive integer giving the interval of values of $n$ for which simulations should be conducted. Larger values give more accurate results.
#' @param eta A real number between 0 and 1 giving the probability of misclassification error in the training data.
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param subsample_size An integer giving the size of the initial 'pilot' sample to be simulated. If left as NULL, the input data size will be used (benchmark SCB).
#' @param nboot An integer giving the number of SCB bootstraps to be performed.
#' @param power A logical indicating whether experimental power based on the predictions should also be reported
#' @param effect_size If `power` is `TRUE`, a real number indicating the scaled effect size the user would like to be able to detect.
#' @param powersims If `power` is `TRUE`, an integer indicating the number of simulations to be conducted at each step to calculate power.
#' @param alpha If `power` is `TRUE`, a real number between 0 and 1 indicating the probability of Type I error to be used for hypothesis testing. Default is 0.05.
#' @param parallel Boolean indicating whether or not to use parallel processing.
#' @param coreoffset If `parallel` is true, a positive integer indicating the number of free threads to be kept unused. Should not be larger than the number of CPU cores.
#' @param packages A list of packages that need to be loaded in order to run `model`.
#' @param method An optional string stating the distribution from which data is to be generated. Default is i.i.d. uniform sampling. Can also take a function outputting a vector of probabilities if the user wishes to specify a custom distribution.
#' @param p If method is 'Class Imbalance', gives the degree of weight placed on the positive class.
#' @param minn Optional argument to set a different minimum n than the dimension of the algorithm. Useful with e.g. regularized regression models such as elastic net.
#' @param x Optional argument for methods that take separate predictor and outcome data. Specifies a matrix-like object containing predictors. Note that if used, the x and y objects are bound together columnwise; this must be handled in the user-supplied helper function.
#' @param y Optional argument for methods that take separate predictor and outcome data. Specifies a vector-like object containing outcome values. Note that if used, the x and y objects are bound together columnwise; this must be handled in the user-supplied helper function.
#' @param backend One of the parallel backends used by [future::plan()]. See function documentation for more details.
#' @param replacement A logical flag indicating whether sampling should be performed with replacement.function.
#' @param ... Additional arguments that need to be passed to `model`
#' @return A `list` containing two named elements. `Raw` gives the exact output of the simulations, while `Summary` gives a table of accuracy metrics, including the achieved levels of \eqn{\epsilon} and \eqn{\delta} given the specified values. Alternative values can be calculated using [getpac()]
#' @seealso [plot_accuracy()], to represent simulations visually, [getpac()], to calculate summaries for alternate values of \eqn{\epsilon} and \eqn{\delta} without conducting a new simulation, and [gendata()], to generated synthetic datasets.
#' @importFrom parallel detectCores
#' @importFrom future plan cluster
#' @importFrom furrr future_map furrr_options
#' @import dplyr
#' @importFrom stats model.frame quantile
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
#' \donttest{
#' }
#' @export


estimate_accuracy <- function(formula,
                              model,
                              data=NULL,
                              dim=NULL,
                              maxn=NULL,
                              sparse=FALSE,
                              density=NULL,
                              upperlimit=NULL,
                              nsample= 30,
                              steps= 50,
                              eta=0.05,
                              delta=0.05,
                              epsilon=0.05,
                              predictfn = NULL,
                              subsample_size = NULL,
                              nboot = 1L,
                              power = FALSE,
                              effect_size=NULL,
                              powersims=NULL,
                              alpha=0.05,
                              parallel = TRUE,
                              coreoffset=0,
                              packages=list(),
                              method = c("Uniform","Class Imbalance"),
                              p=NULL,
                              minn = ifelse(is.null(data),ifelse(is.null(x),(dim+1),(ncol(x)+1)),(ncol(data)+1)),
                              x=NULL,
                              y=NULL,
                              backend = c("multisession","multicore","cluster","sequential"),
                              replacement=TRUE,
                              ...){
  force(minn)
  split <- FALSE
  backend <- match.arg(backend)
  if(is.null(data) & is.null(x)){
    if(is.null(dim)){
      simpleError("No data provided but no dimension given for data simulation")
    }
    names <- all.vars(formula)
    data <- gendata(model,dim,maxn,predictfn,names,sparse,density,...)
    if(sparse){
      split <- TRUE
      dat <- data
    } else{
      dat <- model.frame(formula,data)
    }
  } else if(!is.null(x)){
    if(is.null(y)){
      simpleError("Predictor matrix specified but no outcome given")
    } else if(is.factor(y)){
      if(length(levels(y))!=2){
        simpleError(paste("Expected binary outcome vector but got number of levels:", length(levels(y))))
      } else{
        y <- ifelse(as.numeric(y)==1,0,1)
      }
    }
    data <- cbind(x,y)
    split <- TRUE
    dat <- data
  } else{
    dat <- model.frame(formula,data)
  }
  if(!is.null(subsample_size)){
    ind <- sample(seq_len(nrow(data)), subsample_size, replace = F)
  } else{
    ind <- sample(seq_len(nrow(data)), nrow(data), replace = F)
  }
  data <- data[ind, , drop = FALSE]
  method <- match.arg(method)
  results <- list()
  outcome <- all.vars(formula)[1]
  nvalues <- round(seq(minn,ifelse(is.null(upperlimit),nrow(dat),upperlimit),steps),0)

  # if(!is.null(predictfn)){
  #   predict.svrclass <- predictfn
  # }
  ####function printing out the minimum sample size that achieves the highest accuracy
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
    cl <- 1L
    if(backend != "sequential"){
      stop("Parallel set to false but parallel backend has been given. Please set backend to 'sequential'")
    }
  }
  if(backend == "sequential"){
    plan(get(backend))
  } else if(backend== "cluster"){
    cl <- makeClusterPSOCK(availableWorkers(), revtunnel = FALSE)
    plan(get(backend), workers = cl)
  } else{
    plan(get(backend), workers = cl)  # Use chosen backend
  }

  progressor_p  <- progressor(steps = length(nvalues))
  temp <- function(x,method,p,dat,model,eta,nsample,outcome,power,effect_size,powersims,alpha,split,predictfn,replacement,...){
    progressor_p()
    #set.seed(as.numeric(Sys.time()))
    lapply(packages, library, character.only = TRUE)
    r <- acc_sim(n=x,method=method,p=p,dat=dat,model=model,eta=eta,nsample=nsample,outcome=outcome,power=power,effect_size=effect_size,powersims=powersims,alpha=alpha,split=split,predictfn=predictfn,replacement=replacement,...)
    return(r)
  }
  out <- vector("list",nboot)
  for(i in seq_len(nboot)){
    results <- future_map(nvalues, temp,method=method,p=p,dat=dat,model=model,eta=eta,nsample=nsample,outcome=outcome,power=power,effect_size=effect_size,powersims=powersims,alpha=alpha,split=split,predictfn=predictfn,replacement=replacement,...,.options = furrr_options(seed = TRUE))
    results <- bind_rows(results)
    summtable <- results %>% group_by(n) %>% summarise(Accuracy = mean(accuracy,na.rm=T),
                                                       Precision = mean(prec,na.rm=T),
                                                       Recall = mean(rec,na.rm=T),
                                                       Fscore = mean(fscore,na.rm=T),
                                                       Delta = mean((1-accuracy) > epsilon,na.rm=T),
                                                       Epsilon = quantile((1-accuracy),(1-delta),na.rm=T),
                                                       Power = mean(pwr,na.rm=T))
    out[[i]] <- list("Raw"=results,"Summary"=summtable)
    class(out[[i]]) <- "scb_data"
  }
  if(length(out)>1) {class(out) <- "scb_list"}
  if(nboot==1L) return(out[[1]]) else return(out)
}
#' Recalculate achieved sample complexity bounds given different parameter inputs
#'
#' @param table A list containing an element named `Raw`. Should always be used with the output of [estimate_accuracy()]
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @return A `list` containing two named elements. `Raw` gives the exact output of the simulations, while `Summary` gives a table of accuracy metrics, including the achieved levels of \eqn{\epsilon} and \eqn{\delta} given the specified values. Alternative values can be calculated using [getpac()] again.
#' @seealso [plot_accuracy()], to represent simulations visually, [getpac()], to calculate summaries for alternate values of \eqn{\epsilon} and \eqn{\delta} without conducting a new simulation, and [gendata()], to generated synthetic datasets.
#' @import dplyr
#' @importFrom stats quantile
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
#' \donttest{
#' }
#' @export
getpac <- function(table,epsilon=0.05,delta=0.05){
  results <- table$Raw
  summtable <- results %>% group_by(n) %>% summarise(Accuracy = mean(accuracy,na.rm=T),
                                                     Precision = mean(prec,na.rm=T),
                                                     Recall = mean(rec,na.rm=T),
                                                     Fscore = mean(fscore,na.rm=T),
                                                     Delta = mean((1-accuracy) > (epsilon),na.rm=T),
                                                     Epsilon = quantile((1-accuracy),(1-delta),na.rm=T),
                                                     Power = mean(pwr,na.rm=T))
  out <- list("Raw"=results,"Summary"=summtable)
  class(out) <- "scb_data"
  return(out)
}
#' Plot method for simulated sample complexity bounds (`scb_data` object)
#'
#' This method plots performance metrics estimated using [estimate_accuracy()] for objects of class `"scb_data"`.
#'
#' @param x An object of class `"scb_data"`, typically the output of [estimate_accuracy()]
#' @param metrics A character vector containing the metrics to display in the plot. Can include any of "Accuracy", "Precision", "Recall", "Fscore", "Delta", "Epsilon", "Power".
#' @param plottype A string indicating the graphics system to use. Must be either `"ggplot"` or `"plotly"`.
#' @param letters A string specifying whether `"Delta"` and `"Epsilon"` should appear as Greek letters (`"greek"`) or Latin (`"latin"`) in the legend. Defaults to `"greek"`.
#'
#' @return A \link[ggplot2]{ggplot} or \link[plotly]{plot_ly} object, depending on the value of `plottype`.
#' @seealso [estimate_accuracy()], which generates the data used for plotting.
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @importFrom ggplot2 ggplot aes geom_line
#' @importFrom plotly plot_ly layout
#' @examples
#' \donttest{
#' }
#' @export
#' @method plot scb_data


plot.scb_data <- function(x,
                          metrics = c("Accuracy", "Precision", "Recall", "Fscore", "Delta", "Epsilon", "Power"),
                          plottype = c("ggplot", "plotly"),
                          letters = c("greek", "latin")){
  letters <- match.arg(letters)
  toplot <- x$Summary %>%
   select(n,all_of(metrics)) %>%
      pivot_longer(cols=Accuracy:Power,names_to = "Metric",values_to = "Value")
  if(("Delta" %in% metrics) & (letters == "greek")){
   toplot$Metric[which(toplot$Metric == "Delta")] <-'\u03B4'
  }
  if(("Epsilon" %in% metrics) & (letters== "greek")){
    toplot$Metric[which(toplot$Metric == "Epsilon")] <-'\u03B5'
  }
  plottype <- match.arg(plottype)
  if(plottype=="ggplot"){
    plot<- ggplot(toplot,aes(x=n,y=Value,col=Metric,linetype=Metric)) + geom_line()
  } else if(plottype=="plotly"){
    plot <- plot_ly(toplot, x= ~n, y = ~Value, color = ~Metric, mode = 'lines') %>%
      layout(yaxis = list(title= ""), legend = list(orientation = ''))  #TODO - add ability to change plotly options
  }else{
    simpleError("Invalid Plot Type")
  }
  return(plot)
}


