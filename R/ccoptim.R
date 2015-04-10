## ==========================================================================
##
##  Adaptation to compiled code of the R-functions optim and optimize
##
##  Copyright (C) for optim and optimize 2000-12 The R Core Team
##
##  Adaptation for use with compiled code by Karline Soetaert
##
## ==========================================================================

ccoptim <-
    function(par, fn, gr = NULL, ...,
             method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),   #,Brent not allowed: use optimize
             lower = -Inf, upper = Inf, 
             control = list(), hessian = FALSE, 
             data = NULL, dllname = NULL, initfunc = NULL,
             rpar = NULL, ipar = NULL)
{

  cl <- attributes(fn)$call
  if (! is.null(cl))
    if (cl != "compile.optim")
      stop ("problem is not compiled with 'compile.optim'")
      
  if (is.list(fn)) {            
      
      if (! is.null(data) & ! is.null(Dnames <- attributes(fn)$datanames))
       if (getnames(data) != Dnames )
         stop ("'data' not compatible with data used to compile 'fn'")
      if (!is.null(gr) & "jacfunc" %in% names(fn))
         stop("If 'fn' is a list that contains 'jacfunc', argument 'gr' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(fn))
         stop("If 'fn' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(fn))
         stop("If 'fn' is a list that contains initfunc, argument 'initfunc' should be NULL")

     dllname <- fn$dllname
     gr <- fn$jacfunc  
     initfunc <- fn$initfunc
     fn <- fn$func
  }

    if (! isvalid(fn))
      stop("'fn' should be either a compiled function or a character vector")

    if (! is.null(gr))
      if (! isvalid(gr))
        stop("'gr' should be either a compiled function or a character vector")

    if (is.character(fn) & is.null(dllname))  
      stop("'dllname' should have a value if 'fn' is a character")

    if (sum(duplicated (c(fn, gr, initfunc))) > 0)
      stop("fn, gr, initfunc cannot be the same")

    if (class (fn) == "CFunc")
       fn1 <- body(fn)[[2]]
     else if (is.loaded(fn, PACKAGE = dllname, type = "") ||
       is.loaded(fn, PACKAGE = dllname, type = "Fortran"))  {
       fn1 <- getNativeSymbolInfo(fn, PACKAGE = dllname)$address
     } else 
       stop(paste("'fn' not loaded ", fn))

    gr1 <- NULL
    if (!is.null(gr)) {
      if (class (gr) == "CFunc")
        gr1 <- body(gr)[[2]]
      else if (is.loaded(gr, PACKAGE = dllname, type = "") ||
        is.loaded(gr, PACKAGE = dllname, type = "Fortran"))  {
        gr1 <- getNativeSymbolInfo(gr, PACKAGE = dllname)$address
      } else 
        stop(paste("'gr' not loaded ", gr))
    }

    initfunc1 <- NULL
    ndat <- ncoldat <- 0
    if (!is.null(initfunc)) { # initfunc is associated with initialisation of data
      if (is.null(data))
        stop("if 'initfunc' has a value, 'data' should be passed")
      data <- as.matrix(data)  
      ndat <- nrow(data)
      ncoldat <- ncol(data)    
      storage.mode(data) <- "double"  
    
      if (class (initfunc) == "CFunc")
        initfunc1 <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        initfunc1 <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else 
        stop(paste("'initfunc' not loaded ", initfunc))
    } else {
      if (!is.null(data))
        stop("if 'data' is passed, 'initfunc' should be present to initialise them")
    }
    
    method <- match.arg(method)
    if((length(lower) > 1L || length(upper) > 1L ||
       lower[1L] != -Inf || upper[1L] != Inf)
       && !any(method == c("L-BFGS-B"))) {
	      warning("bounds can only be used with method L-BFGS-B ")
      	method <- "L-BFGS-B"  
    }
    npar <- length(par)
    
    if (is.null(rpar))
      rpar <- 0.
    if (is.null(ipar))
      ipar <- 0  
    
    ## Defaults :
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, npar),
		ndeps = rep.int(1e-3, npar),
		maxit = 100L, abstol = -Inf, reltol = sqrt(.Machine$double.eps),
		alpha = 1.0, beta = 0.5, gamma = 2.0,
		REPORT = 10,
		type = 1,
		lmm = 5, factr = 1e7, pgtol = 0,
		tmax = 10, temp = 10.0)
    nmsC <- names(con)
    if (method == "Nelder-Mead") con$maxit <- 500
    if (method == "SANN") {
   	 con$maxit <- 10000
	   con$REPORT <- 100
    }
    con[(namc <- names(control))] <- control
    if(length(noNms <- namc[!namc %in% nmsC]))
	    warning("unknown names in control: ", paste(noNms,collapse=", "))
    if(con$trace < 0)
	    warning("read the documentation for 'trace' more carefully")
    else if (method == "SANN" && con$trace && as.integer(con$REPORT) == 0)
	    stop("'trace != 0' needs 'REPORT >= 1'")
    if (method == "L-BFGS-B" && any(!is.na(match(c("reltol","abstol"), namc))))
   	  warning("method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")
    if(npar == 1 && method == "Nelder-Mead")
      warning("one-dimensional optimization by Nelder-Mead is unreliable:\nuse optimize()")
    lower <- as.double(rep_len(lower, npar))
    upper <- as.double(rep_len(upper, npar))
    res <- .Call("call_optim", par, fn1, gr1, initfunc1, method, con, lower, 
      upper, as.integer(ndat), as.integer(ncoldat), data, as.double(rpar), as.integer(ipar), package = "ccSolve")
    if (hessian)
        res$hessian <- .Call("call_optimhess", res$par, fn1, gr1, con, 
          as.double(rpar), as.integer(ipar), PACKAGE = "ccSolve")
    res
}

# ------------------------------------------------------------------------------

ccoptimize <- function(f, interval, ...,
		     lower=min(interval), upper=max(interval),
		     maximum=FALSE, tol=.Machine$double.eps^0.25, 
         dllname = NULL, rpar = NULL, ipar = NULL)
{
  tt <- attributes(f)$call
  if (! is.null(tt))
    if (!tt %in% c("compile.optimize", "compile.uniroot"))
      stop ("problem is not compiled with 'compile.optimize'")

  if (is.list(f)) { 
    if (!is.null(dllname) & "dllname" %in% names(f))
      stop("If 'f' is a list that contains dllname, argument 'dllname' should be NULL")

    dllname <- f$dllname
    f <- f$func
  }

  if (! isvalid(f))
    stop("'f' should be either a compiled function or a character vector")
  
  if (class (f) == "CFunc")
    f1 <- body(f)[[2]]
  else if (is.loaded(f, PACKAGE = dllname, type = "") ||
    is.loaded(f, PACKAGE = dllname, type = "Fortran"))  {
      f1 <- getNativeSymbolInfo(f, PACKAGE = dllname)$address
  } else 
    stop(paste("'f' not loaded ", f))

  if (is.null(rpar))
    rpar <- 0.
  if (is.null(ipar))
    ipar <- 0  

  val <- .Call("cc_do_fmin", f1, as.double(lower), as.double(upper), as.double(tol), 
     as.integer(maximum), as.double(rpar), as.integer(ipar), PACKAGE = "ccSolve")

  if (maximum)
   	list(maximum = val[1], objective = val[2])
  else
  	list(minimum = val[1], objective = val[2])
}

## 
ccoptimise <- ccoptimize

ccfunc.optimize <- function(f, x, rpar = NULL, ipar = NULL)
   f$func(x, f = 1, rpar, ipar)$f


ccfunc.optim <-
    function(fn, par, data = NULL, dllname = NULL, initfunc = NULL,
             rpar = NULL, ipar = NULL)
{

  cl <- attributes(fn)$call
  if (! is.null(cl))
    if (cl != "compile.optim")
      stop ("problem is not compiled with 'compile.optim'")
      
  if (is.list(fn)) {            
      
      if (! is.null(data) & ! is.null(Dnames <- attributes(fn)$datanames))
       if (getnames(data) != Dnames )
         stop ("'data' not compatible with data used to compile 'fn'")
      if (!is.null(dllname) & "dllname" %in% names(fn))
         stop("If 'fn' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(fn))
         stop("If 'fn' is a list that contains initfunc, argument 'initfunc' should be NULL")

     dllname <- fn$dllname
     initfunc <- fn$initfunc
     fn <- fn$func
  }

    if (! isvalid(fn))
      stop("'fn' should be either a compiled function or a character vector")

    if (is.character(fn) & is.null(dllname))  
      stop("'dllname' should have a value if 'fn' is a character")

    if (sum(duplicated (c(fn, initfunc))) > 0)
      stop("fn, initfunc cannot be the same")

    if (class (fn) == "CFunc")
       fn1 <- body(fn)[[2]]
     else if (is.loaded(fn, PACKAGE = dllname, type = "") ||
       is.loaded(fn, PACKAGE = dllname, type = "Fortran"))  {
       fn1 <- getNativeSymbolInfo(fn, PACKAGE = dllname)$address
     } else 
       stop(paste("'fn' not loaded ", fn))

    initfunc1 <- NULL

    ndat <- ncoldat <- 0
    if (!is.null(initfunc)) { # initfunc is associated with initialisation of data
      if (is.null(data))
        stop("if 'initfunc' has a value, 'data' should be passed")
      data <- as.matrix(data)  
      ndat <- nrow(data)
      ncoldat <- ncol(data)    
      storage.mode(data) <- "double"  
    
      if (class (initfunc) == "CFunc")
        initfunc1 <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        initfunc1 <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else 
        stop(paste("'initfunc' not loaded ", initfunc))
    } else {
      if (!is.null(data))
        stop("if 'data' is passed, 'initfunc' should be present to initialise them")
    }
    
    npar <- length(par)
    
    if (is.null(rpar))
      rpar <- 0.
    if (is.null(ipar))
      ipar <- 0  
    

    out <- .Call("call_optim_func", par, fn1, initfunc1,  
      as.integer(ndat), as.integer(ncoldat), data, as.double(rpar), 
      as.integer(ipar), package = "ccSolve")
    out
}
