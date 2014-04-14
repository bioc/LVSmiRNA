\name{read.mir}
\alias{read.mir}

\title{ Read in miRNA Data from Agilent Feature Extractiion Output Files  }
\description{
  Reads intensitity data from a set of one-color microarray image analysis output files.
}
\usage{                                                                                 
read.mir(files=NULL, source="agilent.median", path=NULL, ext=NULL, names=NULL,
              columns=NULL, other.columns=NULL, annotation=NULL, green.only=TRUE,
              wt.fun=NULL, verbose=TRUE, sep="\t", quote=NULL, remove.ctrl=TRUE,\dots)
}
\arguments{
  \item{files}{character vector giving the names of the files containing image analysis output or, for Imagene data, a character matrix of names of files.
  Alternatively, it can be a data.frame containing a column called \code{FileName}.
  If omitted, then all files with extension \code{ext} in the specified directory will be read in alphabetical order.}
  \item{source}{character string specifying the image analysis program which produced the output files.  Choices are \code{"agilent.median"}, \code{"agilent.mean"}.}
  \item{path}{character string giving the directory containing the files.
  The default is the current working directory.}
  \item{ext}{character string giving optional extension to be added to each file name}
  \item{names}{character vector of names to be associated with each array as column name.
  Defaults to \code{removeExt(files)}.}
  \item{columns}{list, or named character vector.
  For two color data, this should have fields \code{R}, \code{G}, \code{Rb} and \code{Gb} giving the column names to be used for red and green foreground and background or, in the case of Imagene data, a list with fields \code{f} and \code{b}.
  For single channel data, the fields are usually \code{E} and \code{Eb}.
  This argument is optional if \code{source} is specified, otherwise it is required.}
  \item{other.columns}{character vector of names of other columns to be read containing spot-specific information}
  \item{annotation}{character vector of names of columns containing annotation information about the probes}
  \item{green.only}{logical, for use with \code{source}, should the
    green (Cy3) channel only be read, or are both red and green
    required?. Standard Agilent MIR data have only one channel so
    defaults to \code{TRUE}.}
  \item{wt.fun}{function to calculate spot quality weights}
  \item{verbose}{logical, \code{TRUE} to report each time a file is read}
  \item{sep}{the field separator character}
  \item{quote}{character string of characters to be treated as quote marks}
  \item{remove.ctrl}{ logical, if \code{TRUE} control probes will not be read}
  \item{\dots}{any other arguments are passed to \code{read.table}}
}

\details{
  This is the main data input function for the LVSmiRNA package for
  one-color microRNA data. It was originally designed to extract the green channel intensities from
  a series of files, produced by Agilent Feature Extractiion software, and
  assembles them into the components of one list. Data from some other
  image analysis programs can be read if the appropriate column names
  containing the intensities are specified using the \code{columns}
  argument (This will work if the column names are unique and if there
  are no incomplete rows in the file after the last line of data. Header
  lines are ok, if appropriately skipped).
  The function is a simple wrapper for
  \code{"\link[=read.maimages]{read.maimages}"} in \code{\link{limma}}
  package so it shares all its features (though right now the input
  source is restricted to agilent type file).
  
  The argument \code{files} should be a matrix with two columns at
  least. One column should contain the names of the samples and the other
  column should contain names of files containing intensity data.
  
  The argument \code{other.columns} allows arbitrary columns of the image
  analysis output files to be reserved in the data object. These become
  matrices in the 'other' component.
}

\value{
  An Elist object.
  
  \item{G}{matrix containing the intensities for each array with probes as rows and arrays as columns.} 
  \item{Gb}{matrix containing the background intensities for each array with probes as rows and arrays as columns.}
  \item{targets}{data frame with column \code{FileName} giving the names of the files read, with column \code{Sample} giving the names of the samplse.} 
  \item{genes}{data frame containing annotation information about the probes, for examples miRNA names and IDs and positions on the array.}
  \item{source}{character string giving the image analysis program name.}
  \item{preprocessing}{list with components \code{Background},
    \code{Normalization}, \code{is.log}, \code{Summarization} indicate
    which pre-processing step has beendone.}  
  
}

\note{ All image analysis files being read are assumed to contain data for the same genelist in the same order. No checking is done to confirm that this is true. Probe annotation information is read from one of the files only.}
\author{ Stefano Calza <stefano.calza@unibs.it>, Suo Chen and Yudi Pawitan. }

\seealso{\code{read.mir} is based on \code{"\link[=read.table]{read.table}"} in the base
  package and modified from \code{"\link[=read.maimages]{read.maimages}"} in the \code{\link{limma}}
  package.}

\examples{
#  Read all intensity files from current working directory
\dontrun{
dir.files <- system.file("extdata", package="LVSmiRNA")
taqman.data <- read.table(file.path(dir.files,"Comparison_Array.txt"),header=TRUE,as.is=TRUE)
MIR <- read.mir(taqman.data)
}}

