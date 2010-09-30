\name{summarize}
\alias{summarize}
\alias{summarize.EList}
\alias{summarize.RGList}
\Rdversion{1.1}

\docType{package}
\title{ LVSmiRNA Summarization Function(s) for microRNA Microarray}

\description{
  Summarize microRNA microarray data objects.
}

\author{ Stefano Calza <stefano.calza@biostatistics.it>, Suo Chen and Yudi Pawitan. }

\usage{
summarize(object, ...)
\method{summarize}{EList}(object,RA,remove.ctrl=FALSE,is.log=!is.null(object$preprocessing$Normalization),
method=c("rlm","medianpolish","mean"),verbose=FALSE,make.exprs=FALSE,...)
\method{summarize}{RGList}(object,RA,remove.ctrl=FALSE,is.log=!is.null(object$preprocessing$Normalization),
method=c("rlm","medianpolish","mean"),verbose=FALSE,make.exprs=FALSE,...)
}

\arguments{
  \item{object}{ an object for which a summary is desired.}
  \item{RA}{ an object from estVC.}
  \item{remove.ctrl}{ logical, indicating whether to remove control
    probes.}
  \item{is.log}{Are data already logged?}
  \item{method}{ currently, method "medianpolish","mean" and "rlm" are
    supported.}
  \item{verbose}{More output}
  \item{make.exprs}{Should the output be and \code{exprSet} object?}
  \item{\ldots}{\ldots}
}
\details{
For multi-probe, multi-replicate microarray, intensities need to be summarized into a single expression value for each miRNA. The data objects are summarized as if they were lists. 
}
\value{  
  An Elist object containing components as follows:
 
 \item{G}{ matrix containing the summarized intensities for each array with miRNAs as rows and arrays as columns.} 
 \item{Gb}{ matrix containing the background intensities for each array with probes as rows and arrays as columns.}
 \item{targets}{ data frame with column \code{FileName} giving the names of the files read, with column \code{Sample} giving the names of the samplse.} 
 \item{genes}{ data frame containing annotation information about the probes, for examples gene names and IDs and positions on the array.}
 \item{source}{ character string giving the image analysis program name.}
 \item{preprocessing}{ list with components \code{Background}, \code{Normalization}, \code{is.log}, \code{Summarization} indicate which pre-processing step has been done.} 
}

\references{
Irizarry et al., 'Exploration, normalization, and summaries of high density oligonucleotide array probe level data', (2003a, Biostatistics); Huber, P. J., 'Robust estimation of a location parameter', (1964, Annuas of Mathematical Statistics)
}

\seealso{
\code{\link{lvs}}, \code{\link{estVC}}
}

\examples{
\dontrun{

data("MIR-spike-in")
AA <- estVC(MIR,method="joint")
dd <- summarize(MIR,RA=AA,method="rlm")

##summarization methods other than rlm, object RA is not required
dd1 <- summarize(MIR,method="medianpolish")
dd2 <- summarize(MIR,method="mean")
}
}
