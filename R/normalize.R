## Time-stamp: <11-10-2010 13:42:16 on Goliath>
#############################################################
##                                                         ##
## normalize.R                                             ## 
## Author: SC                                              ##
## First version: 11 Aug 2009                              ##
##                                                         ##
## Given a EList object this code get LVS miRNAs and       ##
## normalize the arrays using a smooth spline              ##
#############################################################

lvs <- function(RG,RA,ref,proportion=0.7,df=3,method=c("joint","rlm"),
                cov.formula=c("weighted","asymptotic"),
                spar=NULL,normalize.method=c("vsn","smooth.spline","mixed"),
                summarize.args=NULL,stratify=TRUE,n.strata=3,
                #subtract.bg=FALSE,## used only with vsn
                level=c("mir","probe"),Atransf=c("sqrt","log"),
                keep.iset=FALSE,clName,verbose=FALSE,...)
  UseMethod("lvs")

lvs.EList <- function(RG,RA,ref,proportion=0.7,df=3,method=c("joint","rlm"),
                      cov.formula=c("weighted","asymptotic"),
                      spar=NULL,normalize.method=c("vsn","smooth.spline","mixed"),
                      summarize.args=NULL,stratify=TRUE,n.strata=3,
                      ## subtract.bg=FALSE,## used only with vsn
                      level=c("mir","probe"),Atransf=c("sqrt","log"),
                      keep.iset=FALSE, clName, verbose=FALSE,...)
  {
    
    method <- match.arg(method)
    cov.formula <- match.arg(cov.formula)
    normalize.method <- match.arg(normalize.method)
    level <- match.arg(level)
    Atransf <- match.arg(Atransf)
    i.vsn=NULL
    
    if(!is.null(preproc(RG)) && !is.null(preproc(RG)$Summarization))
      stop("Object must not be summarized")

   
    
    if(missing(RA))
      RA <- estVC(RG,method=method,cov.formula=cov.formula,clName=clName,verbose=verbose)

    RA$Y <- eval(call(Atransf,RA$ArrayChi2))


    if(stratify)
      {

        ref <- rowMeans(RA$ArrayEffect)
        
        q.ref <- factor(cut(ref,quantile(ref,seq(0,1,1/n.strata)),right=TRUE,include=TRUE))

        out.list <- split(as.data.frame(do.call("cbind",RA[c("logStdDev","Y")])),
                          q.ref)
        names(out.list) <- NULL
        i.set <- lapply(out.list,function(x)
                        {
                          fit.rq <- rq(x$Y~bs(x$logStdDev,df),tau=proportion)
                          i.tmp <- x$Y < fitted(fit.rq)
                          names(i.tmp) <- rownames(x)
                          return(i.tmp)
                        })

        i.set <- unlist(i.set)
        i.set <- i.set[match(names(RA$logStdDev),names(i.set))]

        if(normalize.method=="mixed")
          {
            i.vsn <- lapply(out.list,function(x)
                            {
                              fit.rq <- rq(x$Y~bs(x$logStdDev,df),tau=0.9)
                              i.tmp <- x$Y < fitted(fit.rq)
                              names(i.tmp) <- rownames(x)
                              return(i.tmp)
                            })
            
            i.vsn <- unlist(i.vsn)
            i.vsn <- i.vsn[match(names(RA$logStdDev),names(i.vsn))]

          }
      }
    else
      {
        fit.rq <- with(RA,rq(Y~bs(logStdDev,df),tau=proportion))
        i.set <- RA$Y < fitted(fit.rq)
        names(i.set) <- names(RA$ArrayChi2)

        if(normalize.method=="mixed")
          {
            fit.rq <- with(RA,rq(Y~bs(logStdDev,df),tau=0.9))
            i.vsn <- RA$Y < fitted(fit.rq)
            names(i.vsn) <- names(RA$ArrayChi2)
          }
          
      }
    
    if(!is.null(summarize.args) && !is.list(summarize.args))
      stop("summarize.args must be a named list. See args(summarize.RGList) or args(summarize.EList).")


    if(level == "probe")
      {
        if(is.null(summarize.args))
          summarize.args <- list(method="medianpolish")
        else
          {
            if(any(summarize.args == "rlm"))
              warning("\rlm\" method is implemented only at \"mir\" level. Switching to default method (medianpolish)")
            summarize.args[["method"]] <- "medianpolish"
          }
        
      }
    else
      {
        if(is.null(summarize.args))
          summarize.args <- list(RA=RA)
        else
          summarize.args[["RA"]] <- RA
      }
    
    switch(level,
           ## summarize + normalize mirs
           "mir"=
           {
             if(!is.null(summarize.args))
               {
                 summarize.args$object = RG
                 RG <- do.call("summarize",summarize.args)
                 obj <- exprs(RG)
                                  
               }
             else
               {
                 RG <- summarize(RG)
                 obj <- exprs(RG)
               }
             
             ## Remember VSN needs original scale!
             if(normalize.method!="vsn")
               {
                 if(missing(ref))
                   ref <- rowMeans(obj)
               }
             else
               {
                 obj <- 2^obj
               }
             
             out <- normalize.lvs(obj,ref=ref,i.set=i.set,spar=spar,
                                  normalize.method=normalize.method,i.vsn=i.vsn,verbose=verbose)

             exprs(RG) <- out
           },
           ## Normalize probes + summarize
           "probe"=
           {
             i.set <- i.set[match(RG$genes$SystematicName,names(i.set))]
             if(normalize.method=="mixed")
               i.vsn <- i.vsn[match(RG$genes$SystematicName,names(i.vsn))]
             
             if(normalize.method!="vsn")
               {
                 if(preproc(RG)$is.log)
                   obj=exprs(RG)
                 else
                   {
                     obj=log2(exprs(RG))
                     preproc(RG)$is.log <- TRUE
                   }
                 
                 if(missing(ref))
                   ref <- rowMeans(obj)
               }
             else
               {
                 if(preproc(RG)$is.log)
                   exprs(RG) <- 2^exprs(RG)
                 
                 obj <- exprs(RG)

                 preproc(RG)$is.log <- TRUE ## it will transformed by vsn
               }
             
             out <- normalize.lvs(obj,ref=ref,i.set=i.set,spar=spar,
                                  normalize.method=normalize.method,i.vsn=i.vsn,verbose=verbose)
             
             exprs(RG) <- out
             
             if(!is.null(summarize.args))
               {
                 if(is.list(summarize.args))
                   {
                     summarize.args$object = RG
                     RG <- do.call("summarize",summarize.args)
                     obj <- exprs(RG)
                   }
                 else
                   stop("summarize.args must be a names list. See args(summarize.RGList).")
               }
             else
               {
                 RG <- summarize(RG)
                 obj <- exprs(RG)
               }
           },
           ## Using ArrayEffects out of RA object
           {
             obj <- RA[["ArrayEffects"]]

             if(normalize.method!="vsn")
               {
                 if(missing(ref))
                   ref <- rowMeans(obj)
               }
             else
               {
                 obj <- 2^obj
                 preproc(RG)$is.log <- TRUE ## it will be transformed after vsn
               }
             
             out <- normalize.lvs(obj,ref=ref,i.set=i.set,spar=spar,
                                  normalize.method=normalize.method,i.vsn=i.vsn,
                                  verbose=verbose)

             exprs(RG) <- out
             
           })

    preproc(RG)$Normalization <- paste("LVS",ifelse(normalize.method=="vsn","-VSN","-smooth"),paste="")
    preproc(RG)$is.log=TRUE

    if(keep.iset)
      attr(RG,"i.set") <- names(i.set)[i.set]

    return(RG)

  }


lvs.RGList <- lvs.EList



normalize.lvs <- function(object,ref,i.set,spar=NULL,
                          normalize.method=c("mixed","vsn","smooth.spline"),
                          i.vsn=NULL,verbose=FALSE)
  {
    
    
    normalize.method <- match.arg(normalize.method)
    
    
    switch(normalize.method,
           smooth.spline=
           {
                          
             for(i in 1:ncol(object))
               {
                 n.curve <- smooth.spline(ref[i.set], object[i.set,i],spar=spar)
                 tmp <- as.numeric(approx(n.curve$y, n.curve$x, 
                                          xout = object[, i], rule = 2)$y)
                 
                 object[,i] <- tmp
               }
           },
           vsn=
           {

             require(vsn)
             
             fit <- vsn2(object[i.set,],verbose=FALSE)
             object <- predict(fit,newdata=object)
           },
           mixed =
           {
             require(vsn)
             
             ## smooth-spline normalization
             for (i in 1:ncol(object))
               {
                 n.curve <- smooth.spline(ref[i.set], object[i.set,i], spar = spar)
                 tmp <- as.numeric(approx(n.curve$y, n.curve$x, xout = object[,i], rule = 2)$y)
                 object[, i] <- tmp
                 
               }

             ## VSN variance stabilization but no calibration
             object <- 2^object ## go back to raw scale

             if(!is.null(i.vsn))
               {
                 fit <- vsnMatrix(object[i.vsn, ], verbose = FALSE, calib="none",
                                  returnData=FALSE)
                 object <- vsn:::vsn2trsf(object,fit@coefficients,strata=rep(1L, nrow(object)),
                                          hoffset=fit@hoffset,calib="none")
               }
             else
               {
                 fit <- vsnMatrix(object, verbose = FALSE,calib="none",
                                  returnData=TRUE)
                 object <- exprs(fit) ## extract fit@hx
               }
           }
           )
    
    object
    
    
  }
