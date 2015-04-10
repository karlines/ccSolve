## ==========================================================================
##
##  Adaptation to compiled code of the R-function uniroot 
##
##  Copyright (C) for uniroot 2000-12 The R Core Team
##
##  note: based on the simpler uniroot version pre R 3.1 
##  ignoring flower and fupper input
##  adaptation by Karline Soetaert
##
## ==========================================================================

ccuniroot <- function (f, interval, ..., lower = min(interval), upper = max(interval), 
    tol = .Machine$double.eps^0.25, maxiter = 1000, 
    dllname = NULL, rpar = NULL, ipar = NULL) 
{

    cl <- attributes(f)$call
    if (! is.null(cl))
      if (!cl %in% c("compile.uniroot", "compile.optimize"))
        stop ("problem is not compiled with 'compile.uniroot'")

  if (is.list(f)) {       
      if (!is.null(dllname) & "dllname" %in% names(f))
         stop("If 'f' is a list that contains dllname, argument 'dllname' should be NULL")

     dllname <- f$dllname
     f <- f$func
  }
 
    if (! isvalid(f))
      stop("'f' should be either a compiled function or a character vector")
    if (is.character(f) & is.null(dllname))  
      stop("'dllname' should have a value if 'fn' is a character")
      if (class (f) == "CFunc")
        f <- body(f)[[2]]
      else if (is.loaded(f, PACKAGE = dllname, type = "") ||
        is.loaded(f, PACKAGE = dllname, type = "Fortran"))  {
        f <- getNativeSymbolInfo(f, PACKAGE = dllname)$address
      } else 
        stop(paste("'f' not loaded ", f))

    if (!missing(interval) && length(interval) != 2L) 
        stop("'interval' must be a vector of length 2")
    if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
        upper) 
        stop("lower < upper  is not fulfilled")
    if (is.null(rpar))
       rpar <- 0.
    if (is.null(ipar))
       ipar <- 0  

    val <- .Call("cc_zeroin", f, as.double(lower), 
        as.double(upper), as.double(tol), as.integer(maxiter),
        as.double(rpar), as.integer(ipar), PACKAGE = "ccSolve")
    iter <- as.integer(val[2L])
    if (iter < 0) {
        warning("_NOT_ converged in ", maxiter, " iterations")
        iter <- maxiter
    }
    list(root = val[1L], f.root = val[4L], iter = iter, 
        estim.prec = val[3L])
}

# Adaptation of function uniroot.all from the R-package rootSolve.

ccuniroot.all <- function (f, interval, lower = min(interval), upper = max(interval), 
    tol = .Machine$double.eps^0.25, maxiter = 1000, n = 100, ..., dllname = NULL, 
    rpar = NULL, ipar = NULL) 
{

  if (is.list(f)) {            ### IF a list
      if (!is.null(dllname) & "dllname" %in% names(f))
         stop("If 'f' is a list that contains dllname, argument 'dllname' should be NULL")

     dllname <- f$dllname
     f <- f$func
  }
    if (is.null(rpar))
       rpar <- 0.
    if (is.null(ipar))
       ipar <- 0  

    fn <- f
    if (! isvalid(f))
      stop("'f' should be either a compiled function or a character vector")
    if (is.character(f) & is.null(dllname))  
      stop("'dllname' should have a value if 'fn' is a character")
      if (class (f) == "CFunc") {
        f <- body(f)[[2]]
        func <- function(x) fn(x, f = 1, rpar, ipar)$f
      } else if (is.loaded(f, PACKAGE = dllname, type = "") ||
        is.loaded(f, PACKAGE = dllname, type = "Fortran"))  {
        f <- getNativeSymbolInfo(f, PACKAGE = dllname)$address
        func <- function(x) .Call(fn, as.double(x), f = 1, as.double(rpar), as.integer(ipar))$f
      } else 
        stop(paste("'f' not loaded ", f))

    if (!missing(interval) && length(interval) != 2L) 
        stop("'interval' must be a vector of length 2")
    if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
        upper) 
        stop("lower < upper  is not fulfilled")

    xseq <- seq (lower, upper, length.out = n+1)
    mod <- mapply(xseq, FUN = function(x) func(x))
    Equi <- xseq[which(mod == 0)]
    ss <- mod[1:n] * mod[2:(n + 1)]
    ii <- which(ss < 0)
    for (i in ii)
      Equi <- c(Equi, .Call("cc_zeroin", f, as.double(xseq[i]), 
        as.double(xseq[i+1]), as.double(tol), as.integer(maxiter),
        as.double(rpar), as.integer(ipar), PACKAGE = "ccSolve")[1])
    return(Equi)
}