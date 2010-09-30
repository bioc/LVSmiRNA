\name{estVC}
\alias{estVC}
\alias{estVC.EList}
\alias{estVC.RGList}
\title{ Robust Linear Model to Estimate Residual Variance and Array Effect}
\description{
  Given intensities from microRNA data, fits a robust linear model at probe level and return the residual standard deviations and the array effects.
}
\usage{
estVC(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),
                      do.parallel=FALSE,cores=NULL,use.snow=FALSE,type.cluster=c("MPI","SOCK"),
                      name.cluster = "cl")
\S3method{estVC}{RGList}(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),
                      do.parallel=FALSE,cores=NULL,use.snow=FALSE,type.cluster=c("MPI","SOCK"),
                      name.cluster = "cl")
\S3method{estVC}{EList}(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),
                      do.parallel=FALSE,cores=NULL,use.snow=FALSE,type.cluster=c("MPI","SOCK"),
                      name.cluster = "cl")
}

\arguments{
  \item{object}{ an object of class \code{EList} or \code{RGList}. }
  \item{method}{ character string specifying the estimating algorithm to be used. Choices are "joint" and "rlm".}
  \item{cov.formula}{ character string specifying the covariance formula to be used. Choices are "weighted" and "asymptotic".}

  \item{do.parallel}{Perform computation in parallel. Default to FALSE}
  \item{cores}{Number of cores/clusters to use if do.parallel=TRUE}
  \item{use.snow}{Parallel computing based on package
    \code{snow}. Default to FALSE which means using \code{multicore}}
  \item{type.cluster}{Argument for \code{snow}}
  \item{name.cluster}{Character. Cluster name. Make sure no other object with that
    name exists.}
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
}}
\keyword{ normalization }
\keyword{ miRNA }
