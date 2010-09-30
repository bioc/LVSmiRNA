#############################################################
##                                                         ##
## summarize.R                                             ## 
## Author: SC                                              ##
## First version: 11 Aug 2009                              ##
##                                                         ##
## Given a EList object this code summarize probe level    ##
## values using a RMA-like  method                         ##
#############################################################


summarize <- function(object,...)
  UseMethod("summarize")

## Returns data on log2 scale
summarize.EList <- function(object,RA,remove.ctrl=FALSE,
                            is.log=!is.null(object$preprocessing$Normalization),
                            method=c("rlm","medianpolish","mean"),
                            verbose=FALSE,make.exprs=FALSE,...)
  {

        
    method <- match.arg(method)


    if(!is.null(object$preprocessing$Summarization))
      stop(class(object)[1]," object already summarized (method:",
           object$preprocessing$Summarization,")")
    
    list.names <- unique(featureNames(object))        
    
    
    if(method == "rlm")
      {
        if(missing(RA))
          stop("RA argument required for method \"rlm\"")
        out <- RA[["ArrayEffects"]]

        if(remove.ctrl)
          out <- out[object$genes$ControlType == 0,]
        
      }
    else
      {
        
        if(preproc(object)$is.log)
          is.log <- preproc(object)$is.log

        ## control sequence are removed while reading
        if(remove.ctrl)
          object <- object[object$genes$ControlType == 0,]
        
        out <- NULL
        
        if (verbose)
          {
            ## This requires affy to be loaded
            pbt <- new("ProgressBarText", length(list.names), barsteps = as.integer(20))
            open(pbt)
          }
        
        if(!is.log)
          exprs(object) <- log2(exprs(object))
        
        for(i in 1:length(list.names))
          {
            
            if(verbose)
              updateMe(pbt)
            
            ## SC: bug fixed 11 Sept 09

            PROBE <- object$genes[object$genes$SystematicName == list.names[i],"ProbeName"]
            xdati <- exprs(object)[object$genes$SystematicName == list.names[i],]
            rownames(xdati) <- make.names(PROBE,unique=TRUE)
            
            if(method == "medianpolish")
              {
                
                xdati <- reShapeMed(xdati,PROBE)
                fit <-medpolish2(xdati,trace.iter=FALSE,na.rm=TRUE)
                out <- rbind(out,fit$overall+fit$col) ## overall+column=array effect
              }
            if(method=="mean")
              {
                fit <- colMeans(xdati)
                out <- rbind(out,fit)
              }
          }
        
        if(verbose)
          close(pbt)
      } ## END else
    
    exprs(object) <- out
    rownames(exprs(object)) <- list.names

    preproc(object) <- list(Background=object$preprocessing$Background,
                            Normalization=object$preprocessing$Normalization,
                            is.log=TRUE,Summarization=method)
    if(make.exprs)
      {
        phdata <- as(object$target,"AnnotatedDataFrame")
        object <- new("ExpressionSet",exprs=exprs(object),pheno=phdata)
        preproc(object) <- list(Background=object$preprocessing$Background,
                                Normalization=object$preprocessing$Normalization,
                                is.log=TRUE,Summarization=method)
        
      }
    else
      object$genes <- object$genes[!duplicated(object$genes$SystematicName),
                                   -which(names(object$genes)=="ProbeName"),
                                   drop=FALSE]

    
    
    return(object)
    
  }

summarize.RGList <- function(object,RA,remove.ctrl=FALSE,
                            is.log=!is.null(object$preprocessing$Normalization),
                            method=c("rlm","medianpolish","mean"),
                            verbose=FALSE,make.exprs=FALSE,...)
  {
    summarize.EList(object=object,RA=RA,remove.ctrl=remove.ctrl,
                    is.log=is.log,
                    method=method,verbose=verbose,make.exprs=make.exprs,...)
    
  }

