####################################################
## YP: 4 October 2009
## Customized version on the rlm function to allow
## 	- flexible weights
## 	- joint mean-variance modelling
##
####################################################

RLM <- function(formula, maxit=20, k=1.345, data, model=TRUE,
                na.action, method=c("joint","rlm"), x=TRUE, y=TRUE, offset,
                cov.formula=c("weighted","asymptotic"), start=NULL,...)
{

  method <- match.arg(method)
  cov.formula <- match.arg(cov.formula)
  
  
  call <- match.call()
  if(missing(data)) data <- environment(formula)
  
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula","data"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf,parent.frame())
  mt <- attr(mf,"terms")
  
  Y <- model.response(mf,"numeric")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)){
    if(length(offset)==1)
      offset <- rep(offset,NROW(Y))
    else if(length(offset)!=NROW(Y))
      stop(gettextf("number of offsets is %d, should equal %d(number of observations)"),
           length(offset),NROW(y),domain=NA)
  }
  
  if(is.empty.model(mt)){
    X <- NULL
    fit <- list(coefficients=if(is.matrix(Y)) matrix(,0,3) else numeric(0),
                residuals=if(is.matrix(Y)) nrow(Y) else length(Y))
    if(!is.null(offset)){
      fit$fitted.values <- offset
      fit$residuals <- Y-offset
    }
  }
  else{
    X <- model.matrix(mt,mf)
    X <- as.matrix(X)
  }
  
  fit <- rlmFit(x=X, y=Y, maxit=maxit, k=k, offset=offset,
                cov.formula=cov.formula, method=method, start=start,...)
  if(model)
    fit$model <- mf
  fit$na.action <- attr(mf,"na.action")
  fit$x <- X
  fit$y <- Y
  fit <- c(fit,list(call = call, formula=formula, terms = mt,
                    offset=offset, xlevels = .getXlevels(mt,mf)))
  class(fit) <- c("rlm","lm")
  fit
}

## RLM function
rlmFit <- function(x, y, maxit=20L, k=1.345, offset=NULL, 
                   method=c("joint","rlm"), cov.formula=c("weighted","asymptotic"),
                   start=NULL, error.limit=0.01)
{



  method <- match.arg(method)
  cov.formula <- match.arg(cov.formula)
  
  if(is.null(n <- nrow(x)))
    stop("'x'must be a matrix")
  if(n==0)
    stop("0 cases")
  n_obs <- NROW(y)
  n_vars <- ncol(x)
  
  ifelse(is.null(offset),offset <- rep.int(0,n_obs),y <- y-offset)
  if(n_obs!=n)
    stop("incompatible dimensions")

  
  ## if NULL start starting values
  initialized = as.numeric(!is.null(start))
  storage.mode(initialized) <- "integer"
  storage.mode(maxit) <- "integer"
  
  x <- as.matrix(x)
  
  ## mean-model only
  switch(method,
         rlm={fit <- .Call("rlm_fit",x,y,start,maxit,initialized,.Machine$double.eps,PACKAGE="LVSmiRNA")},
         joint={fit <- .Call("joint_fit",x,y,start,maxit,initialized,.Machine$double.eps,PACKAGE="LVSmiRNA")},
         stop("Error in method argument"))

  res <- fit$resids
  new <- fit$Beta
  sigma <- fit$sigma
  wt <- fit$weights
  pred <- fit$predicted
  iterations <- fit$iter
  converged <- fit$converged


  ## covariance matrix

  switch(cov.formula,
         asymptotic={

           psi <- res
           psi[abs(res) > k*sigma] <- k*sigma
           
           psi_d <- rep(0,length(res))
           psi_d[abs(res) <= k*sigma] <- 1
           
           A <- t(x)%*%x
           B <- diag(1, nrow(A))
           storage.mode(B) <- "double"

           C <- solve(A, B)
           
           cov_matrix <- (mean(psi^2)/(mean(psi_d))^2 * C)
           std_error =  sqrt(diag(cov_matrix))
         },
         weighted={
           ## use sandwich formula: See Pawitan IOL, page 394
           wt2 = wt^2 * res^2
           J = crossprod(x*wt2,x)
           
           A <- crossprod(x*wt,x)
           B <- diag(1, nrow(A))
           storage.mode(B) <- "double"
           Iinv <- solve(A,B)
           
           cov_matrix = Iinv %*% J %*% Iinv
           std_error <- sqrt(diag(cov_matrix))
         })
  
  res_sd <- sqrt(sum(res^2)/(n_obs-n_vars))
  names(new) = colnames(x)

  return(list(x=x, wt=wt, sigma=sigma, coefficients=new,
              Std.Error=std_error, t.value=new/std_error,
              cov.matrix=cov_matrix, res.SD=res_sd, residuals=res,
              fitted.value=pred,
              converged=converged,iteration=iterations))
  

}
