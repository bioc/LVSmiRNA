\name{plot-method}
\docType{methods}
\alias{plot-method}
\alias{plot,RA-method}
\title{Plot of Residual	Variance and Array Effect}
\description{Plots results from \code{estVC}}
\usage{
\S4method{plot}{RA}(x,Atransf=c("both","sqrt","log"), abline=c("none","rq"),df=3,proportion=.7,
                    col="black",col.rq="red")
}

\arguments{
  \item{x}{An object of class \code{RA} resulting from \code{estVC}.}
  \item{Atransf}{Transformation to apply at Array Effect}
  \item{abline}{Add a line to the plot representing a quantile fit}
  \item{df}{Degrees of freedom of the quantile regression}
  \item{proportion}{Quantile to fit}
  \item{col}{Color for plotting points}
  \item{col.rq}{Color for plotting quantile line}
}

\references{ Calza et al., 'Normalization of oligonucleotide arrays based on the least variant set of genes', (2008, BMCBioinformatics); Pawitan, Y. 'In All Likelihood: Statistical Modeling and Inference Using Likelihood', (2001, Oxford University Press); Huber, P. J., 'Robust estimation of a location parameter', (1964, Annuas of Mathematical Statistics).
}
\author{ Stefano Calza <stefano.calza@biostatistics.it>, Suo Chen and Yudi Pawitan.}

\seealso{\code{\link{estVC}},\code{\link{rq}}}
\examples{
\dontrun{

# Starting from an EList object called MIR
data("MIR-spike-in")
AA <- estVC(MIR,method="joint")
plot(AA)
 }}
\keyword{ normalization }
\keyword{ miRNA }
