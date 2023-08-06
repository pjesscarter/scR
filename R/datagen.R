gendata <- function(model,dim,maxn,predictfn=NULL,varnames=NULL, ...){
  if(!is.null(predictfn)){
    predict.svrclass <- predictfn
  }
  dat <- matrix(runif((dim*maxn),-100,100),nrow=maxn)
  coefs <- runif(dim,0,1) 
  y <- dat %*% coefs #TODO - support nonlinear boundaries
  dat <- data.frame(y=as.numeric(y>0),dat)
  m <- model(y~.,dat,...)
  outpred <- suppressWarnings(predict(m,dat[,2:(dim+1)]))
  dat$y <- outpred
  colnames(dat) <- varnames
  return(dat)
}
