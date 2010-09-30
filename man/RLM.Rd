\name{RLM}
\alias{RLM}
\title{ Robust Fitting of Linear Models }
\description{
  Fit a linear model by robust regression using the Huber estimator. 
}
\usage{
RLM(formula, maxit=20, k=1.345, data, model=TRUE,na.action,
method=c("joint","rlm"), x=TRUE, y=TRUE,
offset,cov.formula=c("weighted","asymptotic"), start=NULL,...)
}    
\arguments{
  \item{formula}{ a formula of the form y ~ x1 + x2 + ... }
  \item{maxit}{ the limit on the number of IWLS iterations. }
  \item{k}{ tuning constant used for Huber proposal 2 scale estimation. }
  \item{data}{ data frame from which variables specified in formula are preferentially to be taken. }
  \item{model}{ should the model frame be returned in the object? }
  \item{na.action}{ A function to specify the action to be taken if NAs
    are found. The 'factory-fresh' default action in R is
    \code{\link{na.omit}}, and can be changed by
    \code{\link{options}}. }
  \item{method}{ currently, method="rlm" and "joint" are supported. }
  \item{x}{ should the model frame be returned in the object? }
  \item{y}{ should the model matrix be returned in the object? }
  \item{offset}{ numeric of length n. This can be used to specify an a
    priori known component to be included in the linear predictor during
    fitting. }
  \item{cov.formula}{ are the methods to compute covariance matrix, currently either weighted or asymptotic. }
  \item{start}{ vector containing starting values for the parameters in the predictor. }
  \item{\dots}{ \code{\dots} }
}                                 
\details{
  Fitting is done by iterated re-weighted least squares (IWLS). This
  customized version of robust linear model deal with wild ouliers using
  log link in joint modelling heterogeneous variance of covariates.
}

\value{
An object of class "RLM" inheriting from "lm".
}
\references{ Pawitan, Y. 'In All Likelihood: Statistical Modeling and
  Inference Using Likelihood', (2001, Oxford University Press); Huber,
  P. J. , Robust Statistics, (1981. Wiley).  }
\author{ Stefano Calza <stefano.calza@biostatistics.it>, Suo Chen and Yudi Pawitan.}
\seealso{ \code{RLM} is modified from \code{"\link[=rlm]{rlm}"} in the \code{MASS}, \code{"\link[=rlmFit]{rlmFit}"}}
\examples{

set.seed(133)
n <- 9 
p <- 3
X <- matrix(rnorm(n * p), n,p) 
y <- rnorm(n)

fit <- RLM(y~X-1) #no intercept

}
