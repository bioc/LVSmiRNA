\name{rlmFit}
\alias{rlmFit}
\title{ Fitter Functions for Robust Linear Models }
\description{
  These are the basic computing engines called by \code{\link{RLM}} used to fit robust linear models. These should not be used
directly unless by experienced users.
}
\usage{
rlmFit(x, y, maxit=20L, k=1.345, offset=NULL,method=c("joint","rlm"),
cov.formula=c("weighted","asymptotic"),start=NULL, error.limit=0.01)
}
\arguments{
  \item{x}{ design matrix of dimension n * p. }
  \item{y}{ vector of observations of length n, or a matrix with n rows. }
  \item{maxit}{ the limit on the number of IWLS iterations. }
  \item{k}{ tuning constant used for Huber proposal 2 scale estimation. }
  \item{offset}{ numeric of length n. This can be used to specify an a
    priori known component to be included in the linear predictor during
    fitting. }
    \item{method}{ currently, only method="rlm.fit" is supported. }
  \item{cov.formula}{ are the methods to compute covariance matrix, currently either weighted or asymptotic. }
  \item{start}{ vector containing starting values for the paramter estimates. }
  \item{error.limit}{ the convergence criteria during iterative estimation. }
}
\value{
a list with components
  \item{coeffecients }{p vector}
  \item{Std.Error }{p vector}
  \item{t.value }{p vector}
  \item{cov.matrix }{matrix of dimension p*p}
  \item{res.SD }{value of residual standard deviation}
  ...
}
\references{ Yudi Pawitan: In All Likelihood: Statistical modeling and
  inference using likelihood. Oxford University  Press. 2001. }

\author{ Stefano Calza <stefano.calza@biostatistics.it>, Suo Chen and Yudi Pawitan. }

\seealso{ \code{\link{RLM}} which you should use for robust linear regression usually. }

\examples{

set.seed(133)
n <- 9 
p <- 3
X <- matrix(rnorm(n * p), n,p) #no intercept
y <- rnorm(n)

RLM.fit <- rlmFit (x=X, y=y)
}
