gendata <- function(model,dim,maxn,predictfn=NULL,...){
  if(!is.null(predictfn)){
    predict.svrclass <- predictfn
  }
  dat <- data.frame(matrix(runif((dim*maxn),-100,100),nrow=maxn),y=rbinom(maxn,1,0.5))
  m <- model(y~.,dat,...)
  outpred <- suppressWarnings(predict(m,dat[,1:dim]))
  dat$y <- outpred
}
