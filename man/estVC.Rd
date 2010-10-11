\name{estVC}
\alias{estVC}
\alias{estVC.EList}
\alias{estVC.RGList}
\title{ Robust Linear Model to Estimate Residual Variance and Array Effect}
\description{
  Given intensities from microRNA data, fits a robust linear model at probe level and return the residual standard deviations and the array effects.
}
\usage{

estVC(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),clName,verbose=FALSE)
\S3method{estVC}{RGList}(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),clName,verbose=FALSE)
\S3method{estVC}{EList}(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),clName,verbose=FALSE)

}

\arguments{
  \item{object}{ an object of class \code{EList} or \code{RGList}. }
  \item{method}{ character string specifying the estimating algorithm to be used. Choices are "joint" and "rlm".}
  \item{cov.formula}{ character string specifying the covariance formula to be used. Choices are "weighted" and "asymptotic".}
  \item{clName}{Cluster object produced by \code{makeCluster}
    function. Used only if \code{snow} is loaded.}

  \item{verbose}{Print some debug messages.}
}

\details{
\code{estVC} is the first step in LVS normalization. It fits a robust linear model at the probe-level data in order to estimate the variability of probe intensities due to array-to-array variability. Depending on whether probes show considerable differences in within-probe variance, user can choose the more complex \code{joint} model to accommodate the potential heteroscedasticity or standard robust linear model if within-probe variance can be ignored.

The array effects are then captured by the chi-square statistic. The covariance matrix can be estimated based either on the sandwich form of weighted covariance matrix or an asymptotic form.
}

\value{
An object of class \code{RA}  containing three components as follows:

 \item{ArrayEffects}{a matrix containing the array effect with samples as columns and miRNAs as rows.}
 \item{ArrayChi2}{vector giving chi-square statisitcs of the miRNAs as a measure of array-to-array variability.}
 \item{logStdDev}{vector giving standard deviations of the genes on log scale.}
}
\references{ Calza et al., 'Normalization of oligonucleotide arrays based on the least variant set of genes', (2008, BMCBioinformatics); Pawitan, Y. 'In All Likelihood: Statistical Modeling and Inference Using Likelihood', (2001, Oxford University Press); Huber, P. J., 'Robust estimation of a location parameter', (1964, Annuas of Mathematical Statistics).
}
\author{ Stefano Calza <stefano.calza@biostatistics.it>, Suo Chen and Yudi Pawitan.}

\seealso{\code{\link{read.mir}}, \code{\link{lvs}}}
\examples{
\dontrun{

# Starting from an EList object called MIR
data("MIR-spike-in")
AA <- estVC(MIR,method="joint")

# Parellel execution using multicore

library(multicore)

# use this to set the desided number of
#cores. Otherwise multicore would use all the available
options(cores=8)

AA <- estVC(MIR,method="joint")

detach('package:multicore')

# Parellel execution using snow


library(snow)

cl <- makeCluster(8,type="SOCK")

# Or also...see ?makeCluster
# cl <- makeCluster(8,type="MPI")

AA <- estVC(MIR,method="joint",clName=cl)

}}
\keyword{ normalization }
\keyword{ miRNA }
