#############################################################
##                                                         ##
## estVC.R                                                 ## 
## Author: SC                                              ##
## First version: 11 Aug 2009                              ##
## Update
## YP: 4 Oct 2009
## 	a big change on estVC.std: 
##	- option for method='rlm' or 'joint'
##	  cov.formula
##	- output list includes ArrayEffects
##                                                         ##
## Given a EList object estimates Array and probe effects  ##
#############################################################

estVC <- function(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),clName,
                  verbose=FALSE)
  UseMethod("estVC")


estVC.EList <- function(object,method=c("joint","rlm"),cov.formula=c("weighted","asymptotic"),clName,
                        verbose=FALSE)
  {


    method <- match.arg(method)
    cov.formula <- match.arg(cov.formula)
    list.names <- featureNames(object)

    ord <- order(sampleNames(object))
    object <- object[,ord]
    
    ## Parallel version

    ## let's create a list and apply a function in a parallel fashion

    dati.full <- data.frame(exprs(object),probe=factor(probeNames(object)))
    in.list <- split(dati.full,object$genes$SystematicName)

    do.fun <- function(x,is.log)
      {

        probe <- factor(x$probe)
        x <- as.matrix(x[,-which(names(x)=="probe")])
        rownames(x) <- as.character(probe)
        
        if(is.log)
          y <- as.vector(x)
        else
          y <- log2(as.vector(x))


        if(nlevels(probe) >= 2)
          {
            dMat <- data.frame(sample = factor(rep(colnames(x),each=nrow(x))),
                               probe = factor(rep(probe,ncol(x))))
            ## starting value
            ##             cdat = cbind(y,dMat$sample,dMat$probe)
            ##             medpol= medpolish0(cdat)
            
            cdat <- reShapeMed(x,probe=probe)
            medpol <- medpolish2(cdat)
            ##             arr.effect = medpol$array[-1]- medpol$array[1]
            ##             prb.effect = medpol$probe[-1]- medpol$probe[1]
            ##             start = c(medpol$const, arr.effect, prb.effect)
            arr.effect = medpol$col[-1]- medpol$col[1]
            prb.effect = medpol$row[-1]- medpol$row[1]
            start = c(medpol$overall, arr.effect, prb.effect)
            
            fit <- RLM(y~sample+probe,data=dMat, start=start, 
                        cov.formula=cov.formula, method=method)
          }
        else
          {
            dMat <- data.frame(sample = factor(rep(colnames(x),each=nrow(x))))
            tmp = tapply(y, dMat$sample, median)
            start = c(tmp[1], tmp[-1]- tmp[1])
            fit <- RLM(y~sample, data=dMat, start=start, 
                       cov.formula=cov.formula, method=method)
          }
        wf <- grep("sample",names(coef(fit)))
        arr.effect = fit$coef[wf]
        covmat <- fit$cov.matrix

        A <- covmat[wf,wf]
        B <- as.matrix(arr.effect)
        storage.mode(A) <- "double"
        storage.mode(B) <- "double"
        C <- .Call("La_dgesv", A, B, .Machine$double.eps , PACKAGE = "base")
        
        chi2 =  sum(arr.effect * C)
        ## chi2 =  sum(arr.effect * solve(covmat[wf,wf],arr.effect))
        
        ArrayEffects = c(fit$coef[1], fit$coef[1]+arr.effect)

        ans <- c(ArrayEffects,ArrayChi2=chi2,StdDev=fit$res.SD)

        return(ans)
        
      }
    

    out.lst <- .fapply(X=in.list,FUN=do.fun,is.log=preproc(object)$is.log,clName=clName,verbose=verbose)

    out <- do.call("rbind",out.lst)

    out <- out[list.names,] ## Put the features in the original order. Recall: split would use alphabetical order
    
    ArrayEffects <- out[,-which(colnames(out) %in% c("ArrayChi2","StdDev"))]
    ArrayChi2 <- out[,"ArrayChi2"]
    StdDev <- out[,"StdDev"]
    rm(out)


    colnames(ArrayEffects) = sampleNames(object) ## Recall, object has now a different order compared to input!
    ## reset original ordering in ArrayEffects.
    ArrayEffects <- ArrayEffects[,order(ord)]
    names(ArrayChi2) = names(StdDev) = list.names


        
    out <- list(ArrayEffects=ArrayEffects,
                ArrayChi2=ArrayChi2,
                logStdDev=log(StdDev))

    out <- new("RA",out)

    return(out)
  }


estVC.RGList <- estVC.EList
