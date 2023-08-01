estimate_accuracy <- function(data, nsample= 30, steps= 50){
  results <- list()
  #nvalues <- seq(4,300,15)
  nvalues <- seq((ncol(data)+1),nrow(data),steps)
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
        indices <- sample(seq_len(nrow(br)),n)
        sample <- br[indices,c(2:7,9,11)] %>% ungroup() %>% #Note that errors should be drawn at sampling stage,
          #otherwise there is a 0 prob of observing some values ex ante
          mutate(error = rbinom(nrow(.),1,eta),
                 outobs = as.factor(ifelse(!error, two_year_recid,!two_year_recid)))
        m <- tryCatch({train(sample[,1:7],as.factor(sample$outobs),
                             method = "rf",
                             #family = "binomial",
                             trControl= trainControl("cv"))},
                      error=function(e){
                        skip <<- T}
        ) #Think carefully about how training choices matter
        # form <- as.formula(paste("outobs ~ ",paste(colnames(sample)[1:7],sep=" + ",collapse=" + ")))
        # m <- tryCatch({glm(form,data=sample,family=binomial())},
        #               error=function(e){
        #                 skip <<- T}
        #) #Think carefully about how training choices matter
      }
      pred <- predict(m,br[,c(2:7,9)])
      accuracy[j] <- mean(as.numeric(levels(pred)[pred])== br$two_year_recid)
      prec[j] <- tryCatch({precision(table(levels(pred)[pred],br$two_year_recid), relevant = 1)},
                          error = function(e){return(NA)})
      rec[j] <- tryCatch({recall(table(levels(pred)[pred],br$two_year_recid), relevant = 1)},
                         error = function(e){return(NA)})
      fscore[j] <- tryCatch({F_meas(table(levels(pred)[pred],br$two_year_recid), relevant = 1)},
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

#print(summtable[which.max(summtable$Accuracy), ]$n)
#Fix delta = 0.1, find smallest n that guarantees epsilon = 0.05 - no such n exists
#print(summtable[cumsum(summtable$epsilon < epsilon)==1, ]$n)
#summtable2 <- summtable[,c(1:5,7)]
#summtable %>% 
#  pivot_longer(cols=Accuracy:epsilon,names_to = "Metric",values_to = "Value") %>%
#  ggplot(.,aes(x=n,y=Value,col=Metric)) + geom_line() 

#epsilon <- '\u03B5'
#colnames(summtable2)[6] <- epsilon

#data<- summtable2 %>% 
  #pivot_longer(cols=Accuracy:epsilon,names_to = "Metric",values_to = "Value") 

#fig <- plot_ly(data, x= ~n, y = ~Value, color = ~Metric, mode = 'lines') %>%
  #layout(yaxis = list(title= ""), legend = list(orientation = ''))

