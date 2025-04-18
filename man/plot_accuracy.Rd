% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{plot_accuracy}
\alias{plot_accuracy}
\title{Represent simulated sample complexity bounds graphically}
\usage{
plot_accuracy(
  table,
  metrics = c("Accuracy", "Precision", "Recall", "Fscore", "Delta", "Epsilon", "Power"),
  plottype = c("ggplot", "plotly"),
  letters = c("greek", "latin")
)
}
\arguments{
\item{table}{A list containing an element named \code{Raw}. Should always be used with the output of \code{\link[=estimate_accuracy]{estimate_accuracy()}}}

\item{metrics}{A character vector containing the metrics to display in the plot. Can be any of "Accuracy", "Precision", "Recall", "Fscore", "delta", "epsilon"}

\item{plottype}{A string giving the graphics package to be used to generate the plot. Can be one of "ggplot" or "plotly"}

\item{letters}{A string determining whether delta and epsilon should be given as greek letters in the plot legend. Defaults to Greek lettering but available in case of rendering issues.}
}
\value{
Either a \link[ggplot2]{ggplot} or \link[plotly]{plot_ly} plot object, depending on the chosen option of \code{plottype}.
}
\description{
Represent simulated sample complexity bounds graphically
}
\examples{
mylogit <- function(formula, data){
m <- structure(
  glm(formula=formula,data=data,family=binomial(link="logit")),
  class=c("svrclass","glm")  #IMPORTANT - must use the class svrclass to work correctly
)
return(m)
}
mypred <- function(m,newdata){
out <- predict.glm(m,newdata,type="response")
out <- factor(ifelse(out>0.5,1,0),levels=c("0","1"))
#Important - must specify levels to account for possibility of all
#observations being classified into the same class in smaller samples
return(out)
}
\donttest{
library(parallel)
results <- estimate_accuracy(two_year_recid ~ race + sex + age +
      juv_fel_count + juv_misd_count + priors_count +
      charge_degree..misd.fel.,mylogit,br,predictfn = mypred,
    nsample=10,
    steps=1000,
    coreoffset = (detectCores() -2))

fig <- plot_accuracy(results,letters="latin")
fig
}
}
\seealso{
\code{\link[=estimate_accuracy]{estimate_accuracy()}}, to generate estimated sample complexity bounds.
}
