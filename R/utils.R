setMethod("boxplot",signature="EList",function(x,...)
          {
            boxplot(as.data.frame(x$E),outline=FALSE)
          })


setMethod("exprs",signature="EList",function(object)
          {
            object$E
          })

setMethod("exprs<-",signature="EList",function(object,value)
          {
            object$E <- value
            object
          })



setMethod("preproc",signature="EList",function(object)
          {
            object$preprocessing
          })


setMethod("preproc<-",signature="EList",function(object,value)
          {
            object$preprocessing <- value
            object
            
          })


setMethod("sampleNames",signature="EList",function(object)
          {
            colnames(exprs(object))
            
          })


setMethod("featureNames",signature="EList",function(object)
          {
            unique(object$genes$SystematicName)
            
          })


setMethod("probeNames",signature="EList",function(object)
          {
            object$genes$ProbeName
            
          })

## RGList

setMethod("boxplot",signature="RGList",function(x,...)
          {
            boxplot(as.data.frame(x$G),outline=FALSE)
          })


setMethod("exprs",signature="RGList",function(object)
          {
            object$G
          })

setMethod("exprs<-",signature="RGList",function(object,value)
          {
            object$G <- value
            object
          })



setMethod("preproc",signature="RGList",function(object)
          {
            object$preprocessing
          })


setMethod("preproc<-",signature="RGList",function(object,value)
          {
            object$preprocessing <- value
            object
            
          })


setMethod("sampleNames",signature="RGList",function(object)
          {
            colnames(exprs(object))
            
          })

setMethod("featureNames",signature="RGList",function(object)
          {
            unique(object$genes$SystematicName)
            
          })


setMethod("probeNames",signature="RGList",function(object)
          {
            object$genes$ProbeName
            
          })


################################################################################
################################################################################

setClass("RA",
         representation=representation("list"))

setMethod("plot","RA",
          function(x,Atransf=c("both","sqrt","log"),
                   abline=c("none","rq"),df=3,proportion=.7,
                   col="black",col.rq="red")
          {
            
            Atransf <- match.arg(Atransf)
            
            abline <- match.arg(abline)
            
            switch(Atransf,
                   "sqrt"=
                   {
                     ylab = expression(sqrt("Array Effect"))
                     y = sqrt(x$ArrayChi2)
                     plot(y=y,x=x$logStdDev,pch=16,cex=.3,xlab="log Std Dev",
                          ylab=ylab,col=col)
                     switch(abline,
                            rq = 
                            {
                              fit.rq <- rq(y~bs(x$logStdDev,df=df),tau=proportion)
                              lines(sort(x$logStdDev),fitted(fit.rq)[order(x$logStdDev)],col=col.rq,lwd=1.5)
                            },
                            invisible(NULL))
                     
                   },
                   "log"=
                   {
                     ylab = "log(Array Effect)"
                     y = log(x$ArrayChi2)
                     plot(y=y,x=x$logStdDev,pch=16,cex=.3,xlab="log Std Dev",
                          ylab=ylab,col=col)
                     switch(abline,
                            rq = 
                            {
                              fit.rq <- rq(y~bs(x$logStdDev,df=df),tau=proportion)
                              lines(sort(x$logStdDev),fitted(fit.rq)[order(x$logStdDev)],col=col.rq,lwd=1.5)
                            },
                            invisible(NULL))
                   },
                   {
                     par(mfrow=c(1,2))
                     ylab = expression(sqrt("Array Effect"))
                     y = sqrt(x$ArrayChi2)
                     plot(y=y,x=x$logStdDev,pch=16,cex=.3,xlab="log Std Dev",
                          ylab=ylab,col=col)
                     switch(abline,
                            rq = 
                            {
                              fit.rq <- rq(y~bs(x$logStdDev,df=df),tau=proportion)
                              lines(sort(x$logStdDev),fitted(fit.rq)[order(x$logStdDev)],col=col.rq,lwd=1.5)
                            },
                            invisible(NULL))
                     ylab = "log(Array Effect)"
                     y = log(x$ArrayChi2)
                     plot(y=y,x=x$logStdDev,pch=16,cex=.3,xlab="log Std Dev",
                          ylab=ylab,col=col)
                     switch(abline,
                            rq = 
                            {
                              fit.rq <- rq(y~bs(x$logStdDev,df=df),tau=proportion)
                              lines(sort(x$logStdDev),fitted(fit.rq)[order(x$logStdDev)],col=col.rq,lwd=1.5)
                            },
                            invisible(NULL))
                     
                   },invisible(NULL))
            
          }
          )



## SC 18 Sept 09 (Ryanair flight!): allow for replications in median polish
reShapeMed <- function(mat,probe)
  {
    tab <- table(probe)
    nrep <- max(tab)            
    nr <- length(unique(probe))
    nc <- ncol(mat)
    if(is.null(rownames(mat)))
      rownames(mat) <- make.names(probe,unique=TRUE)

    cnames <- colnames(mat)
    ## SC 19 Sept 09: allows for unbalanced design
    if(sum(abs(diff(tab))))
      {
        
        mat.2 <- matrix(NA,ncol=ncol(mat),nrow=nrep*nr)
        rownames(mat.2) <- make.names(rep(unique(probe),each=nrep),unique=TRUE)
        mat.2[rownames(mat),] <- mat
        mat <- mat.2
      }

    ## This will order rows according to probe names.
    mat <- mat[order(rownames(mat)),]
    
    oo <- rep(1:nrep,nr*nc)
    
    mat <- array(as.vector(mat)[order(oo)],c(nr,nc,nrep))

    if(is.factor(probe))
      dimnames(mat)[1:2] <- list(levels(probe),cnames)
    else
      dimnames(mat)[1:2] <- list(sort(unique(probe)),cnames)
    
    return(mat)
  }

## 1=rows
## 2=cols
## 3=replication

medpolish2 <- 
  function (x,eps = 0.01, maxiter = 10, trace.iter = FALSE, na.rm = TRUE) 
{
  z <- as.array(x)

  converged <- 0
  
  nr <- dim(z)[1L]
  nc <- dim(z)[2L]
  t <- 0
  r <- numeric(nr)
  c <- numeric(nc)
  oldsum <- 0
  for (iter in 1L:maxiter) {
    rdelta <- apply(z, 1L, median, na.rm = na.rm)
    z <- z - array(rdelta, dim=dim(x))
    r <- r + rdelta
    delta <- median(c, na.rm = na.rm)
    c <- c - delta
    t <- t + delta
    cdelta <- apply(z, 2L, median, na.rm = na.rm)
    z <- z - array(rep(cdelta,each=dim(x)[1]), dim=dim(x))
    
    c <- c + cdelta
    delta <- median(r, na.rm = na.rm)
    r <- r - delta
    t <- t + delta
    newsum <- sum(abs(z), na.rm = na.rm)
    converged <- newsum == 0 || abs(newsum - oldsum) < eps * newsum
    if (converged) 
      break
    oldsum <- newsum
    if (trace.iter) 
      cat(iter, ":", newsum, "\n")
  }
  if (converged) {
    if (trace.iter) 
      cat("Final:", newsum, "\n")
  }
  else warning(gettextf("medianpolish() did not converge in %d iterations", 
                        maxiter), domain = NA)
  names(r) <- dimnames(z)[[1]]
  names(c) <- dimnames(z)[[2]]
  ans <- list(overall = t, row = r, col = c, residuals = z, 
              name = deparse(substitute(x)),iteration=iter)
##  ans <- list(constant = t, probe = r, array = c)
  
  class(ans) <- "medpolish"
  ans
}





## Modified from ShortRead

.fapply <- function(X,FUN,...,clName,verbose)
{
  if (is.loaded("mc_fork", PACKAGE="multicore"))
    {
      mcLapply <- get('mclapply', envir=getNamespace('multicore'))
      if (verbose)
        message("using 'mclapply'")
      results <- mcLapply(X, FUN, ...)
    }
  else if(!is.na(match("package:snow", search())))
    {
      if(missing(clName))
        stop("clName must be supplied. See ?estVC")
      
      mcLapply <- get('parLapply', 'package:snow')
      clEvalQ <- get('clusterEvalQ', 'package:snow')
      
      clEvalQ(clName, library(LVSmiRNA))
      
      if (verbose)
        message("using snow 'parLapply'")
      results <- mcLapply(clName, X, FUN, ...)
    }
  else {
    if (verbose)
      message("using lapply")
    results <- lapply(X, FUN, ...)
  }

  return(results)
}
