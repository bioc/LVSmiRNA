\name{MIR-spike-in}
\alias{MIR-spike-in}
\alias{MIR_RA}
\alias{MIR.RA}
\alias{MIR}
\docType{data}
\title{Data example}
\description{
  Data from a micro-RNA spike-in experiment, extracted from
  scanned images using Agilent Feature Extraction Software.
}
\usage{
data("MIR-spike-in")
data("MIR_RA")
}
\details{
  This dataset is derived from a library of synthetic RNA sequences, cor-
  responding to human mature miRNAs as well as in-house miRNAs with
  particularly similar sequences hybridized on an Agilent Human miRNA Mi-
  croarray 2.0. Data consist of a total of 799 miRNA species
  (excluding control features) for 4 samples organized in two groups A and
  B.
  
  Data, colected with the Agilent Feature Extraction Software, are stored
  in a RGList object with the following components: 
  
   - MIR\$G:          		         'gMeanSignal'  
   - MIR\$Gb:          		         'gProcessedSignal'
   - MIR\$gBGMedianSignal:	         'gBGMedianSignal'
   - MIR\$targets      		         'targets'
   - MIR\$Row       		         'Row'
   - MIR\$Col       		         'Column'
   - MIR\$ProbeUID       		 'Probe ID'
   - MIR\$genes\$ControlType		 'FLAG to specify the sort of feature'
   - MIR\$genes\$ProbeName		  'Probe Name'
   - MIR\$genes\$GeneName		  'microRNA Name'
   - MIR\$genes\$SystematicName		  'microRNA Name'
   - MIR\$genes\$Description		  'Description (not used)'

   \code{MIR.RA} holds an object of class \code{RA} obtained from using
   \code{estVC} on the example data.
 }
 \references{
   Willenbrock H, Salomon J, Barken KIMB, Nielsen FC, Litman T. 2009. Quantitative
   miRNA expression analysis: Comparing microarrays with next-generation sequencing.
   RNA 15: 2028-2034.

   Data available from Genome Expression Omnibus (GEO) database under the series accession number
   GSE14511
 }
 \author{Stefano Calza}
 \seealso{
   \code{read.mir}, \code{estVC}
 }
 \keyword{datasets}
