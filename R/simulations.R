estimate_accuracy <- function(formula, model, data=NULL, dim=NULL,maxn=NULL,nsample= 30, steps= 50,eta=0.05,delta=0.05,epsilon=0.05, predictfn = NULL,...){
  if(is.null(data)){
    data <- gendata(model,dim,maxn,predictfn,...)
  }
  results <- list()
  outcome <- all.vars(formula)[1]
  dat <- model.frame(formula,data)
  #nvalues <- seq(4,300,15)
  nvalues <- seq((ncol(dat)),nrow(dat),steps)
  if(!is.null(predictfn)){
    predict.svrclass <- predictfn
  }
  ####function printing out the minimum sample size that achieves the highest accuracy
  for(i in seq_along(nvalues)){
    n <- nvalues[i]
    accuracy <- vector()
    prec <- vector()
    rec <- vector()
    fscore <- vector()
    for(j in seq_len(nsample)){
      skip <- T
      while(skip){
        skip <- F
        indices <- sample(seq_len(nrow(dat)),n)
        samp <- dat[indices,]
        error <- rbinom(nrow(samp),1,eta)
        samp$outobs <- factor(ifelse(error,!samp[[outcome]],samp[[outcome]]),levels=c("0","1"))
        samp <- samp %>% select(!all_of(outcome))
        m <- tryCatch({model(outobs ~.,data=samp,...)},
                       error=function(e){
                         warning("Model failed to compute, regenerating training data")
                         skip <<- T} #TODO - provide useful error message to help diagnose misuse
        )
      }
      pred <- suppressWarnings({predict(m,dat %>% select(!all_of(outcome)))})
      accuracy[j] <- mean(as.numeric(levels(pred)[pred])== factor(dat[[outcome]],levels=c("0","1")))
      prec[j] <- tryCatch({precision(table(levels(pred)[pred],factor(dat[[outcome]],levels=c("0","1"))), relevant = 1)},
                          error = function(e){return(NA)})
      rec[j] <- tryCatch({recall(table(levels(pred)[pred],factor(dat[[outcome]],levels=c("0","1"))), relevant = 1)},
                         error = function(e){return(NA)})
      fscore[j] <- tryCatch({F_meas(table(levels(pred)[pred],factor(dat[[outcome]],levels=c("0","1"))), relevant = 1)},
                            error = function(e){return(NA)})
      
    }
    results[[i]] <- tryCatch({data.frame(accuracy,prec,rec,fscore,n)},
                             error = function(e){return(NA)})
  }
  results <- bind_rows(results)
  summtable <- results %>% group_by(n) %>% summarise(Accuracy = mean(accuracy,na.rm=T), 
                                                     Precision = mean(prec,na.rm=T), 
                                                     Recall = mean(rec,na.rm=T),
                                                     Fscore = mean(fscore,na.rm=T),
                                                     delta = mean((1-accuracy) > epsilon,na.rm=T),
                                                     epsilon = quantile(1-accuracy,1-delta,na.rm=T))
  return(summtable)
}
plot_accuracy <- function(table,metrics=c("Accuracy","Precision","Recall","Fscore","delta","epsilon"),plottype = c("ggplot","plotly")){
  
  toplot <- table %>% 
   select(n,all_of(metrics)) %>%
      pivot_longer(cols=Accuracy:epsilon,names_to = "Metric",values_to = "Value")
  if("delta" %in% metrics){
   toplot$Metric[which(toplot$Metric == "delta")] <-'\u03B4'
  }
  if("epsilon" %in% metrics){
    toplot$Metric[which(toplot$Metric == "epsilon")] <-'\u03B5'
  }
  plottype <- plottype[1]
  if(plottype=="ggplot"){
    plot<- toplot %>% 
      ggplot(.,aes(x=n,y=Value,col=Metric)) + geom_line() 
  } else if(plottype=="plotly"){
    plot <- plot_ly(toplot, x= ~n, y = ~Value, color = ~Metric, mode = 'lines') %>%
      layout(yaxis = list(title= ""), legend = list(orientation = ''))  #TODO - add ability to change plotly options
  }else{
    simpleError("Invalid Plot Type")
  }
  return(plot)
    
}


