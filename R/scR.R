#Does not need to be user-facing
#Works for any function that uses formula+data interface
#TODO - add support for custom helper function to prepare data for use


#' Utility function to generate data points for estimation of the VC Dimension of a user-specified binary classification algorithm given a specified sample size.
#'
#' @param x An integer giving the desired sample size for which the target function is to be approximated.
#' @param ... Additional model parameters to be specified by the user.
#' @return A real number giving the estimated value of Xi(n), the bounding function
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
#' out <- factor(ifelse(out>0.5,1,0),levels=c("0","1")) #Important - must specify levels to account for possibility of all observations being classified into the same class in smaller samples
#' return(out)
#' }
#' predict.svrclass <- mypred
#' model <- mylogit
#' risk_bounds(5)

risk_bounds <- function(x,...){
  n <- x
  xis <- vector()
  for(j in 1:m){
    skip <- T
    while(skip){
      skip <- F
      x <- replicate(l,rnorm(2*n))
      #In theory the complexity of the target concept is irrelevant - is this true?
      coeff <- rnorm(l)
      y <- as.numeric(apply(x,1,FUN=function(r){r %*% coeff}) > 0)
      dat <- data.frame(x,y)
      indices <- sample(1:n,n)
      W <- dat[indices,]
      Wprime <- dat[-indices,]
      #Flip labels
      Wprime$y <- 1- Wprime$y
      dat <- bind_rows(W, Wprime)
      #Flip labels again to recover correct values
      Wprime$y <- 1- Wprime$y
      traindata <- dat
      traindata$y <- factor(traindata$y,levels=c("0","1"))
      #Models tend to use different syntax so we might need to provide a list of supported functions
      fhat <- tryCatch({model(formula = y ~ .,data = traindata)},
                       error=function(e){
                         warning("Model failed to compute, regenerating training data")
                         skip <<- T}
      )
    }
    
    #Different models also use different predict methods
    RW <- mean(predict(fhat,W[,-(l+1)]) != factor(W$y,levels=c("0","1")))
    RWprime <- mean(predict(fhat,Wprime[,-(l+1)]) != factor(Wprime$y,levels=c("0","1")))
    xis[j] <- abs(RW - RWprime)
  }
  xihat <- mean(xis)
  return(xihat)
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
  #Todo - add gradient (very messy function but straightforward to derive)
}

#' Estimate the Vapnik-Chervonenkis (VC) dimension of an arbitrary binary classification algorithm.
#'
#' @param model A binary classification model supplied by the user. Must take arguments `formula` and `data`
#' @param dim A positive integer giving dimension (number of input features) of the model.
#' @param packages A `list` of strings giving the names of packages to be loaded in order to estimate the model.
#' @param m A positive integer giving the number of simulations to be performed at each design point (sample size value). Higher values give more accurate results but increase computation time.
#' @param k A positive integer giving the number of design points (sample size values) for which the bounding function is to be estimated. Higher values give more accurate results but increase computation time.
#' @param parallel Boolean indicating whether or not to use parallel processing. 
#' @param coreoffset If `parallel` is true, a positive integer indicating the number of free threads to be kept unused. Should not be larger than the number of CPU cores.
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param a Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param a1 Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param a11 Scaling coefficient for the bounding function. Defaults to the value given by Vapnik, Levin and Le Cun 1994.
#' @param ... Additional arguments that need to be passed to `model`
#' @return A real number giving the estimated value of the VC dimension of the supplied model.
#' @seealso [scb()], to calculate sample complexity bounds given estimated VCD.
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
#' out <- factor(ifelse(out>0.5,1,0),levels=c("0","1")) #Important - must specify levels to account for possibility of all observations being classified into the same class in smaller samples
#' return(out)
#' }
#' vcd <- simvcd(model=mylogit,dim=7,m=100,k=100,maxn=500,predictfn = mypred) #Takes about 20 minutes to run on mid-range test machine (Intel i5-11600K CPU)
#' @export
simvcd <- function(model,dim,packages=list(),m=1000,k=1000,maxn=5000,parallel = T,coreoffset=0, predictfn = NULL, a=0.16,a1=1.2,a11=0.14927, ...){
  ngrid <- seq(dim,maxn,(maxn/k))
  if(parallel){
    cl <- detectCores() -coreoffset
    cl <- makeCluster(cl)
  } else{
    cl <- 1
    cl <- makeCluster(cl)
  }

  l<-dim
  clusterExport(cl,varlist = c("l","k","m","model","packages","predictfn"),envir = environment())
  clusterEvalQ(cl=cl,expr={
    library(dplyr)
    lapply(packages, library, character.only = TRUE)
    if(!is.null(predictfn)){
      predict.svrclass <- predictfn
    }
  })
  #Need to work on passing function arguments with pbapply - some issues around parallelization
  xihats <- suppressWarnings({pbsapply(ngrid,risk_bounds,...,cl=cl)})
  stopCluster(cl)
  vcd <- optim((l+1),loss,ngrid=ngrid,xi=xihats, , a=a,a1=a1,a11=a11, method="Brent",lower=1,upper = 2*(max(ngrid)), ...)
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
#' out <- factor(ifelse(out>0.5,1,0),levels=c("0","1")) #Important - must specify levels to account for possibility of all observations being classified into the same class in smaller samples
#' return(out)
#' }
#' scb(epsilon=0.05,delta=0.05,eta=0.05,theor=F,model=mylogit,dim=7,m=100,k=100,maxn=500,predictfn = mypred)
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

