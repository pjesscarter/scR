#' Utility function to generate accuracy metrics, for use with [estimate_accuracy()]
#'
#' @param n An integer giving the desired sample size for which the target function is to be calculated.
#' @param method An optional string stating the distribution from which data is to be generated. Default is i.i.d. uniform sampling. Currently also supports "Class Imbalance". Can also take a function outputting a vector of probabilities if the user wishes to specify a custom distribution.
#' @param p If method is 'Class Imbalance', gives the degree of weight placed on the positive class.
#' @param ntest An integer giving the size of the test set to be drawn.
#' @param rtest An integer giving the number of the test sets to be drawn.
#' @param ... Additional model parameters to be specified by the user.
#' @return A data frame giving performance metrics for the specified sample size.
#' @export
acc_sim <- function(n,method = "Uniform",p=NULL,ntest=1000,rtest=100,...){
  accuracy <- vector()
  prec <- vector()
  rec <- vector()
  fscore <- vector()
  pwr <- vector()
  for(j in seq_len(nsample)){
    skip <- T
    while(skip){
      skip <- F
      if(method=="Uniform"){
        indices <- sample(seq_len(nrow(dat)),n)
        samp <- dat[indices,]
        error <- rbinom(nrow(samp),1,eta)
      } else if(method=="Class Imbalance"){
        if(is.null(p)){
          stop("Class Imbalance method selected. Please supply a class imbalance parameter.")
        }
        probs <- ifelse(dat[[outcome]]==1,p,(1-p))
        indices <- sample(seq_len(nrow(dat)),n,replace=T,prob=probs)
        samp <- dat[indices,]
        error <- rbinom(nrow(samp),1,eta)
      } else{
        probs <- method(dat)
        indices <- sample(seq_len(nrow(dat)),n,replace=T,prob=probs)
        samp <- dat[indices,]
        error <- rbinom(nrow(samp),1,eta)
      }
      
      if(is.factor(samp[[outcome]])){
        samp$outobs <- factor(ifelse(error,!as.numeric(as.character(samp[[outcome]])),as.numeric(as.character(samp[[outcome]]))),levels=c("0","1"))
      } else{
        samp$outobs <- factor(ifelse(error,!samp[[outcome]],samp[[outcome]]),levels=c("0","1"))
      }
      samp <- samp %>% select(!all_of(outcome))
      m <- tryCatch({model(update(formula, outobs~.),data=samp#,...
      )},
      error=function(e){
        warning("Model failed to compute, regenerating training data")
        skip <<- T} #TODO - provide useful error message to help diagnose misuse
      )
    }
    accuracyk <- vector()
    preck <- vector()
    reck <- vector()
    fscorek <- vector()
    pwrk <- vector()
    for(k in seq_len(rtest)){

      if(method=="Uniform"){
        indices <- sample(seq_len(nrow(dat)),ntest)
        test <- dat[indices,] %>% select(!all_of(outcome))
        eval <- dat[indices,]
      } else if(method=="Class Imbalance"){
        probs <- ifelse(dat[[outcome]]==1,p,(1-p))
        indices <- sample(seq_len(nrow(dat)),ntest,replace=T,prob=probs)
        test <- dat[indices,] %>% select(!all_of(outcome))
        eval <- dat[indices,]
      } else{
        probs <- method(dat)
        indices <- sample(seq_len(nrow(dat)),ntest,replace=T,prob=probs)
        test <- dat[indices,] %>% select(!all_of(outcome))
        eval <- dat[indices,]
      }
      pred <- suppressWarnings({predict(m,test)})
      accuracyk[k] <- mean(as.numeric(levels(pred)[pred])== factor(eval[[outcome]],levels=c("0","1")))
      preck[k] <- tryCatch({precision(table(levels(pred)[pred],factor(eval[[outcome]],levels=c("0","1"))), relevant = 1)},
                          error = function(e){return(NA)})
      reck[k] <- tryCatch({recall(table(levels(pred)[pred],factor(eval[[outcome]],levels=c("0","1"))), relevant = 1)},
                         error = function(e){return(NA)})
      fscorek[k] <- tryCatch({F_meas(table(levels(pred)[pred],factor(eval[[outcome]],levels=c("0","1"))), relevant = 1)},
                            error = function(e){return(NA)})
      
      if(power){
        reject <- vector()
        Dobs <- as.numeric(levels(pred)[pred])
        Dtrue <- if(is.factor(eval[[outcome]])){as.numeric(as.character(eval[[outcome]]))} else{eval[[outcome]]}
        for(r in 1:powersims){
          Y <- effect_size * Dtrue + rnorm(length(Dobs))
          X <- data.frame(D = Dobs, Y = Y)
          mdl <- lm(Y ~ D, data=X)
          reject[r] <- tryCatch({summary(mdl)$coefficients[2,4] < alpha},
                                error = function(e){return(NA)})
        }
        pwrk[k] <- mean(reject,na.rm=T)
        
      } else{pwrk[k] <- NA}
    }
    accuracy[j] <- mean(accuracyk, na.rm=T)
    prec[j] <- mean(preck, na.rm=T)
    rec[j] <- mean(reck, na.rm=T)
    fscore[j] <- mean(fscorek, na.rm=T)
    pwr[j] <- mean(pwrk, na.rm=T)
    }      
  out <- tryCatch({data.frame(accuracy,prec,rec,fscore,n,pwr)},
                  error = function(e){return(NA)})
  return(out)
}
#' Estimate sample complexity bounds for a binary classification algorithm using either simulated or user-supplied data.
#'
#' @param formula A `formula` that can be passed to the `model` argument to define the classification algorithm
#' @param model A binary classification model supplied by the user. Must take arguments `formula` and `data`
#' @param data Optional. A rectangular `data.frame` object giving the full data from which samples are to be drawn. If left unspecified, [gendata()] is called to produce synthetic data with an appropriate structure.
#' @param dim Required if `data` is unspecified. Gives the horizontal dimension of the data (number of predictor variables) to be generated.
#' @param maxn Required if `data` is unspecified. Gives the vertical dimension of the data (number of observations) to be generated.
#' @param upperlimit Optional. A positive integer giving the maximum sample size to be simulated, if data was supplied.
#' @param nsample A positive integer giving the number of samples to be generated for each value of $n$. Larger values give more accurate results.
#' @param steps A positive integer giving the number of values of $n$ for which simulations should be conducted. Larger values give more accurate results.
#' @param eta A real number between 0 and 1 giving the probability of misclassification error in the training data.
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @param predictfn An optional user-defined function giving a custom predict method. If also using a user-defined model, the `model` should output an object of class `"svrclass"` to avoid errors.
#' @param power A logical indicating whether experimental power based on the predictions should also be reported
#' @param effect_size If `power` is `TRUE`, a real number indicating the scaled effect size the user would like to be able to detect.
#' @param powersims If `power` is `TRUE`, an integer indicating the number of simulations to be conducted at each step to calculate power.
#' @param alpha If `power` is `TRUE`, a real number between 0 and 1 indicating the probability of Type I error to be used for hypothesis testing. Default is 0.05.
#' @param parallel Boolean indicating whether or not to use parallel processing. 
#' @param coreoffset If `parallel` is true, a positive integer indicating the number of free threads to be kept unused. Should not be larger than the number of CPU cores.
#' @param packages A list of packages that need to be loaded in order to run `model`.
#' @param method An optional string stating the distribution from which data is to be generated. Default is i.i.d. uniform sampling. Can also take a function outputting a vector of probabilities if the user wishes to specify a custom distribution.
#' @param p If method is 'Class Imbalance', gives the degree of weight placed on the positive class.
#' @param ntest An integer giving the size of the test set to be drawn.
#' @param rtest An integer giving the number of the test sets to be drawn.
#' @param ... Additional arguments that need to be passed to `model`
#' @return A `list` containing two named elements. `Raw` gives the exact output of the simulations, while `Summary` gives a table of accuracy metrics, including the achieved levels of $\epsilon$ and $\delta$ given the specified values. Alternative values can be calculated using [getpac()]
#' @seealso [plot_accuracy()], to represent simulations visually, [getpac()], to calculate summaries for alternate values of $\epsilon$ and $\delta$ without conducting a new simulation, and [gendata()], to generated synthetic datasets.
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
#' br <- scR::br
#' results <- estimate_accuracy(two_year_recid ~ race + sex + age + juv_fel_count + juv_misd_count + priors_count + charge_degree..misd.fel.,mylogit,br,predictfn = mypred)
#' @export


estimate_accuracy <- function(formula, model,data=NULL, dim=NULL,maxn=NULL,upperlimit=NULL,nsample= 30, steps= 50,eta=0.05,delta=0.05,epsilon=0.05,predictfn = NULL,power = F,effect_size=NULL,powersims=NULL,alpha=0.05,parallel = T,coreoffset=0,packages=list(),method = "Uniform",p=NULL,ntest=1000,rtest=100,...){
  if(is.null(data)){
    names <- all.vars(formula)
    data <- gendata(model,dim,maxn,predictfn,names,...)
  }
  results <- list()
  outcome <- all.vars(formula)[1]
  dat <- model.frame(formula,data)
  #nvalues <- seq(4,300,15)
  nvalues <- seq((ncol(dat)+1),ifelse(is.null(upperlimit),nrow(dat),upperlimit),steps)
  # if(!is.null(predictfn)){
  #   predict.svrclass <- predictfn
  # }
  ####function printing out the minimum sample size that achieves the highest accuracy
  if(parallel){
    cl <- detectCores() -coreoffset
    cl <- makeCluster(cl)
  } else{
    cl <- 1
    cl <- makeCluster(cl)
  }
  clusterExport(cl,varlist = c("formula","dat","model","eta","packages","predictfn","nsample","outcome","power","effect_size","powersims","alpha"),envir = environment())
  clusterEvalQ(cl=cl,expr={
    library(dplyr)
    library(caret)
    lapply(packages, library, character.only = TRUE)
    if(!is.null(predictfn)){
      predict.svrclass <- predictfn
    }
  })
  results <-   suppressWarnings({pblapply(nvalues,acc_sim,method,p,ntest,rtest,cl=cl)})
  stopCluster(cl)
  results <- bind_rows(results)
  summtable <- results %>% group_by(n) %>% summarise(Accuracy = mean(accuracy,na.rm=T), 
                                                     Precision = mean(prec,na.rm=T), 
                                                     Recall = mean(rec,na.rm=T),
                                                     Fscore = mean(fscore,na.rm=T),
                                                     Delta = mean((1-accuracy) > epsilon,na.rm=T),
                                                     Epsilon = quantile((1-accuracy),(1-delta),na.rm=T),
                                                     Power = mean(pwr,na.rm=T))
  return(list("Raw"=results,"Summary"=summtable))
}
#' Recalculate achieved sample complexity bounds given different parameter inputs
#'
#' @param table A list containing an element named `Raw`. Should always be used with the output of [estimate_accuracy()]
#' @param delta A real number between 0 and 1 giving the targeted maximum probability of observing an OOS error rate higher than `epsilon`
#' @param epsilon A real number between 0 and 1 giving the targeted maximum out-of-sample (OOS) error rate
#' @return A `list` containing two named elements. `Raw` gives the exact output of the simulations, while `Summary` gives a table of accuracy metrics, including the achieved levels of $\epsilon$ and $\delta$ given the specified values. Alternative values can be calculated using [getpac()] again.
#' @seealso [plot_accuracy()], to represent simulations visually, [getpac()], to calculate summaries for alternate values of $\epsilon$ and $\delta$ without conducting a new simulation, and [gendata()], to generated synthetic datasets.
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
#' br <- scR::br
#' results <- estimate_accuracy(two_year_recid ~ race + sex + age + juv_fel_count + juv_misd_count + priors_count + charge_degree..misd.fel.,mylogit,br,predictfn = mypred)
#' resultsalt <- getpac(results,epsilon=0.5,delta=0.3)
#' print(resultsalt$Summary)
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
  return(list("Raw"=results,"Summary"=summtable))
}
#' Represent simulated sample complexity bounds graphically
#'
#' @param table A list containing an element named `Raw`. Should always be used with the output of [estimate_accuracy()]
#' @param metrics A character vector containing the metrics to display in the plot. Can be any of "Accuracy", "Precision", "Recall", "Fscore", "delta", "epsilon"
#' @param plottype A string giving the graphics package to be used to generate the plot. Can be one of "ggplot" or "plotly"
#' @return Either a [ggplot2] or [plotly] plot object, depending on the chosen option of `plottype`.
#' @seealso [estimate_accuracy()], to generate estimated sample complexity bounds.
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
#' br <- scR::br
#' results <- estimate_accuracy(two_year_recid ~ race + sex + age + juv_fel_count + juv_misd_count + priors_count + charge_degree..misd.fel.,mylogit,br,predictfn = mypred)
#' fig <- plot_accuracy(results)
#' fig
#' @export
plot_accuracy <- function(table,metrics=c("Accuracy","Precision","Recall","Fscore","Delta","Epsilon","Power"),plottype = c("ggplot","plotly")){
  
  toplot <- table$Summary %>% 
   select(n,all_of(metrics)) %>%
      pivot_longer(cols=Accuracy:Power,names_to = "Metric",values_to = "Value")
  if("Delta" %in% metrics){
   toplot$Metric[which(toplot$Metric == "Delta")] <-'\u03B4'
  }
  if("Epsilon" %in% metrics){
    toplot$Metric[which(toplot$Metric == "Epsilon")] <-'\u03B5'
  }
  plottype <- plottype[1]
  if(plottype=="ggplot"){
    plot<- toplot %>% 
      ggplot(.,aes(x=n,y=Value,col=Metric,linetype=Metric)) + geom_line() 
  } else if(plottype=="plotly"){
    plot <- plot_ly(toplot, x= ~n, y = ~Value, color = ~Metric, mode = 'lines') %>%
      layout(yaxis = list(title= ""), legend = list(orientation = ''))  #TODO - add ability to change plotly options
  }else{
    simpleError("Invalid Plot Type")
  }
  return(plot)
    
}


