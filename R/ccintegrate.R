## ==========================================================================
##
##  adapted from File src/library/stats/R/integrate.R
##  by Karline Soetaert
##  Original Copyright (C) 1995-2012 The R Core Team
##
## ==========================================================================

ccintegrate <- function(f, lower, upper, ..., subdivisions = 100L,
                      rel.tol = .Machine$double.eps^.25,
                      abs.tol = rel.tol, stop.on.error = TRUE,
                      keep.xy = FALSE, aux = NULL, dllname = NULL,
                      rpar = NULL, ipar = NULL)
{
    cl <- attributes(f)$call
    if (! is.null(cl))
      if (cl != "compile.integrate")
        stop ("problem is not compiled with 'compile.integrate'")

    if (is.list(f)) {       
      if (!is.null(dllname) & "dllname" %in% names(f))
         stop("If 'f' is a list that contains dllname, argument 'dllname' should be NULL")

     dllname <- f$dllname
     f <- f$func
    }
    ff <- NULL
    if (! isvalid(f))
      stop("'f' should be either a compiled function or a character vector")
    if (is.character(f) & is.null(dllname))  
      stop("'dllname' should have a value if 'f' is a character")
    if (class (f) == "CFunc") {
        f <- body(f)[[2]]
        ff <- function(x) f(n = length(x), x, f = rep(1, n), rpar, ipar)$f

      } else if (is.loaded(f, PACKAGE = dllname, type = "") ||
        is.loaded(f, PACKAGE = dllname, type = "Fortran"))  {
        f <- getNativeSymbolInfo(f, PACKAGE = dllname)$address
        n <- length(x)
        ff <- function(x) .Call("f", n = length(x), as.double(x), f = rep(1, n), as.double(rpar), as.integer(ipar))$f

      } else 
        stop(paste("'f' not loaded ", f))

    if (is.null(rpar))
       rpar <- 0.
    if (is.null(ipar))
       ipar <- 0  


    limit <- as.integer(subdivisions)
    if (limit < 1L || (abs.tol <= 0 &&
 	     rel.tol < max(50*.Machine$double.eps, 0.5e-28)))
   	  stop("invalid parameter values")
    if(is.finite(lower) && is.finite(upper)) {
	    wk <- .Call("cc_call_dqags", f, as.double(lower), as.double(upper),
	  		as.double(abs.tol), as.double(rel.tol), limit = limit, 
        as.double(rpar), as.integer(ipar))
    } else { # indefinite integral
	if(is.na(lower) || is.na(upper)) stop("a limit is missing")
	if (is.finite(lower)) {
	    inf <- 1
	    bound <- lower
	} else if (is.finite(upper)) {
	    inf <- -1
	    bound <- upper
	} else {
	    inf <- 2
	    bound <- 0.0
	}
	wk <- .Call("cc_call_dqagi", f, as.double(bound), as.integer(inf),
			as.double(abs.tol), as.double(rel.tol),	limit = limit,
      as.double(rpar), as.integer(ipar))
    }
    res <- wk[c("value", "abs.error", "subdivisions")]
    res$message <-
	switch(wk$ierr + 1L,
	       "OK",
	       "maximum number of subdivisions reached",
	       "roundoff error was detected",
	       "extremely bad integrand behaviour",
	       "roundoff error is detected in the extrapolation table",
	       "the integral is probably divergent",
	       "the input is invalid")
    if(wk$ierr == 6L || (wk$ierr > 0L && stop.on.error)) stop(res$message)
    res$call <- match.call()
    class(res) <- "integrate"
    res
}