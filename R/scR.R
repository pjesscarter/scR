

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
#' @importFrom dplyr bind_rows
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
      y <- as.numeric(rowSums(x * coeff) > 0)  # Vectorized operation

      W_idx <- sample(1:(2 * n), n)
      W <- x[W_idx, ]
      Wprime <- x[-W_idx, ]
      Wprime_y <- 1 - y[-W_idx]  # Flip labels

      y_train <- c(y[W_idx], Wprime_y)
      x_train <- rbind(W, Wprime)

      # Train model
      fhat <- tryCatch(
        model(formula = y ~ ., data = cbind(y_train,x_train), ...
              ),
        error = function(e) NULL
      )
      if (!is.null(fhat)) break  # Exit loop if model trains successfully

      gc()
      if (attempt == 10) stop("Failed to train model after 100 attempts.")
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
#' @importFrom pbapply pbsapply
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
  temp <- function(x,l,m,model,packages,predictfn,density, ...){
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

