#Does not need to be user-facing
#Works for any function that uses formula+data interface
#TODO - add support for custom helper function to prepare data for use
simvcd <- function(x,...){
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
loss <- function(h,ngrid,xi){
  #These constants are calculated in Vapnik, Levin and Le Cun 1994
  #based on the known VC dimension of linear discriminant functions 
  #and rely on the assumption that they are universal for all classifiers
  a <- 0.16
  a1 <- 1.2
  a11 <- 0.14927
  ratio <- ngrid/h
  phi <- a*((log((2*ratio))+1)/(ratio-a11))*(1+sqrt((1+((a1*(ratio - a11))/(1+log(2*ratio))))))
  test <- ngrid < (h/2)
  phihats <- ifelse(test,1,phi)
  devs <- (xi - phihats)^2
  out <- sum(devs)
  return(out)
  #Todo - add gradient (very messy function but straightforward to derive)
}

#Note - make parallelization optionally removable
risk_bounds <- function(model,dim,packages,m=1000,k=1000,maxn=5000,coreoffset=0,...){
  ngrid <- seq(dim,maxn,(maxn/k))
  cl <- detectCores() -coreoffset
  cl <- makeCluster(cl)
  l<-dim
  clusterExport(cl,varlist = c("l","k","m","model","packages"),envir = environment())
  clusterEvalQ(cl=cl,expr={
    library(dplyr)
    lapply(packages, library, character.only = TRUE)
  })
  #Need to work on passing function arguments with pbapply - some issues around parallelization
  xihats <- pbsapply(ngrid,simvcd,...,cl=cl)
  stopCluster(cl)
  vcd <- optim((l+1),loss,ngrid=ngrid,xi=xihats,method="Brent",lower=1,upper = 2*(max(ngrid)))
  return(vcd$par)
}

