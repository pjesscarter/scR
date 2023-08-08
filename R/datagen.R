#' Simulate data with appropriate structure to be used in estimating sample complexity bounds
#'
#' @param model A binary classification model supplied by the user. Must take arguments `formula` and `data`
#' @param dim Gives the horizontal dimension of the data (number of predictor variables) to be generated.
#' @param maxn Gives the vertical dimension of the data (number of observations) to be generated.
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param varnames An optional character vector giving the names of variables to be used for the generated data
#' @param ... Additional arguments that need to be passed to `model`
#' @return A `data.frame` containing the simulated data. 
#' @seealso [estimate_accuracy()], to estimate sample complexity bounds given the generated data
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
#' formula <- two_year_recid ~ race + sex + age + juv_fel_count + juv_misd_count + priors_count + charge_degree..misd.fel.,mylogit
#' br <- gendata(mylogit,7,7214,mypred,all.vars(formula))
#' results <- estimate_accuracy(formula,br,predictfn = mypred)
#' @export

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
