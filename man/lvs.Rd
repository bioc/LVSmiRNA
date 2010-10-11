\name{lvs}
\alias{lvs}
\alias{lvs.EList}
\alias{lvs.RGList}
\alias{normalize.lvs}
\title{ Least Variant Set selection and Normalization Function(s) }
\description{
  Selects the Least Variant Set of mircoRNAs, according to the chosen proportion of miRNAs expected not to vary between arrays. Then performs normalization.
}
\usage{
lvs(RG,RA,ref,proportion=0.7,df=3,method=c("joint","rlm"),
                cov.formula=c("weighted","asymptotic"),
                spar=NULL,normalize.method=c("vsn","smooth.spline","mixed"),
                summarize.args=NULL,stratify=TRUE,n.strata=3,
                level=c("mir","probe"),Atransf=c("sqrt","log"),keep.iset=FALSE,clName,
                verbose=FALSE,...)
\S3method{lvs}{RGList}(RG,RA,ref,proportion=0.7,df=3,method=c("joint","rlm"),
                cov.formula=c("weighted","asymptotic"),
                spar=NULL,normalize.method=c("vsn","smooth.spline","mixed"),
                summarize.args=NULL,stratify=TRUE,n.strata=3,
                level=c("mir","probe"),Atransf=c("sqrt","log"),
                keep.iset=FALSE,clName,verbose=FALSE,...)
\S3method{lvs}{EList}(RG,RA,ref,proportion=0.7,df=3,method=c("joint","rlm"),
                cov.formula=c("weighted","asymptotic"),
                spar=NULL,normalize.method=c("vsn","smooth.spline","mixed"),
                summarize.args=NULL,stratify=TRUE,n.strata=3,
                level=c("mir","probe"),Atransf=c("sqrt","log"),keep.iset=FALSE,clName,
                verbose=FALSE,...)
}
\arguments{
  \item{RG}{ an object of class \code{EList} or \code{RGList}}
  \item{RA}{ a list contaning components residual standard deviations,
    chi-square statistics and array effects. It can be computed by
    \code{estVC}. If not provided it will
    computed (slower),}
  \item{proportion}{ the proportion below which miRNAs are expected not to vary between arrays. Default is set to 0.7.}
  \item{ref}{ reference array to be used for normalization. Default is set to mean of array effects across samples. }
  \item{df}{ the desired equivalent number of degrees of freedom(trace of the smooth matrix) in smoothing spline. }
  \item{method}{ character string specifying the estimating algorithm to be used. Choices are "joint" and "rlm".}
  \item{cov.formula}{ character string specifying the covariance formula to be used. Choices are "weighted" and "asymptotic".}
  \item{spar}{ smoothing parameter, typicallly in (0,1].}
  \item{normalize.method}{ character string specifying the normalization method to be used. Choices are "smooth.spline" and "vsn".}
  \item{summarize.args}{ a named list containnig components from argument of \code{summarize}.}
  \item{stratify}{ logical, if TRUE selection of least variant set will be stratified by expression level.}
  \item{n.strata}{ integer giving the number of strata. }
  \item{level}{ character string specifying the normalization performed
    at miRNA level or probe-level.}

  \item{Atransf}{Which transformation to use for Array Effect}
  \item{keep.iset}{return the LVS ids}
  \item{clName}{Cluster object. See \code{estVC}.}
  \item{verbose}{Verbose computation}
  \item{\dots}{ \code{\dots} }
}
\details{
\code{lvs} works by first identifying least variant set (LVS) with the smallest array-to-array variation. The total information extracted from probe-level intensity data of all samples is modeled as a function of array and probe effect in order to select the reference set for normalization. If the residual variances and array effects are available, \code{lvs} runs faster because the step of robust linear modeling has already been done. 

Once the LVS miRNAs are identified, the normalization is performed using \code{VSN} or \code{smooth.spline}.
}

\value{
An object of the same class as RG.

 \item{G}{ matrix containing the normalized intensities for each array with miRNAs as rows and arrays as columns.} 
 \item{Gb}{ matrix containing the background intensities for each array with probes as rows and arrays as columns.}
 \item{targets}{ data frame with column \code{FileName} giving the names of the files read, with column \code{Sample} giving the names of the samplse.} 
 \item{genes}{ data frame containing annotation information about the probes, for examples miRNA names and IDs and positions on the array.}
 \item{source}{ character string giving the image analysis program name.}
 \item{preprocessing}{list with components \code{Background}, \code{Normalization}, \code{is.log}, \code{Summarization} indicate which pre-processing step has been done.}
}
 
\references{ Calza et al., 'Normalization of oligonucleotide arrays based on the least variant set of genes' (2008, BMCBioinformatics).}
\author{ Stefano Calza <stefano.calza@biostatistics.it>, Suo Chen and Yudi Pawitan. }
\seealso{ \code{\link{estVC}}, \code{\link{summarize}}}
\examples{
\dontrun{

# Starting from an Elist object called MIR
data("MIR-spike-in")
AA <- estVC(MIR,method="joint")
bb <- lvs(MIR,RA=AA,level="probe")

##It can also run with object RA missing, but taking longer time
cc <- lvs(MIR)
 }}
\keyword{ normalization }
\keyword{ LVS }
