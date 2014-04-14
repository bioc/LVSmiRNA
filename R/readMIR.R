#############################################################
##                                                         ##
## readMIR.R                                               ## 
## Author: SC                                              ##
## First version: 11 Aug 2009                              ##
##                                                         ##
## 12 April 2014: wrap read.maimages                       ##
##                                                         ##
## Reads in miRNA data                                     ##
#############################################################



read.mir <- function (files = NULL, source = "agilent.median", path = NULL, ext = NULL, 
                      names = NULL, columns = NULL, other.columns = NULL, annotation = NULL, 
                      green.only = TRUE, wt.fun = NULL, verbose = TRUE, sep = "\t", 
                      quote = NULL, remove.ctrl=TRUE, ...)
{
  
  source <- match.arg(source, c("agilent", "agilent.mean","agilent.median"))
#   , "arrayvision", "arrayvision.ARM", "arrayvision.MTM", 
#                                 "bluefuse", "genepix", "genepix.mean", "genepix.median", 
#                                 "genepix.custom", "imagene", "imagene9", "quantarray", 
#                                 "scanarrayexpress", "smd.old", "smd", "spot", "spot.close.open"))
  
  RG <- read.maimages(files=files, source=source, path=path, ext=ext, names=names, columns=columns,
                      other.columns=other.columns, 
                      annotation=annotation, green.only=green.only, wt.fun=wt.fun, verbose=verbose, sep=sep,
                      quote=quote,...)
  
  
    RG$preprocessing <- list(Background=NULL,Normalization=NULL,
                             is.log=FALSE,Summarization=NULL)
    if(remove.ctrl)
      RG <- RG[RG$genes$ControlType == 0,]

  return(RG)
  
}
