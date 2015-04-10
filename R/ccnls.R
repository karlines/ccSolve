## ==========================================================================
##
## Compiled version of nls.lm - nonlinear least squares using a modified
## Levenberg-Marquardt algorithm - based on R-package minpack.lm
## adapted from minpack.lm function 'nls.lm'
## by Karline Soetaert
##
## ==========================================================================

ccnls <- function (fn, data, par, lower = NULL, upper = NULL, jac = NULL, 
  control = nls.lm.control(), dllname = NULL, initfunc = NULL,
    rpar = NULL, ipar = NULL) {

    cl <- attributes(fn)$call
    if (! is.null(cl))
      if (cl != "compile.nls")
        stop ("problem is not compiled with 'compile.nls'")

    if (is.list(fn)) {       
      if (! is.null(Dnames <- attributes(fn)$datanames))
         if (getnames(data) != Dnames )
        stop ("'data' not compatible with data used to compile 'fn'")
      if (!is.null(dllname) & "dllname" %in% names(fn))
         stop("If 'fn' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(jac) & "jac" %in% names(fn))
         stop("If 'fn' is a list that contains jac, argument 'jac' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(fn))
         stop("If 'fn' is a list that contains initfunc, argument 'initfunc' should be NULL")

      dllname <- fn$dllname
      jac <- fn$jacfunc
      initfunc <- fn$initfunc
      fn <- fn$func
   }
 
    if (! isvalid(fn))
      stop("'fn' should be either a compiled function or a character vector")
    if (is.character(fn) & is.null(dllname))  
      stop("'dllname' should have a value if 'fn' is a character")

      if (class (fn) == "CFunc")
        fn1 <- body(fn)[[2]]
      else if (is.loaded(fn, PACKAGE = dllname, type = "") ||
        is.loaded(fn, PACKAGE = dllname, type = "Fortran"))  {
        fn1 <- getNativeSymbolInfo(fn, PACKAGE = dllname)$address
      } else 
        stop(paste("'fn' not loaded ", fn))

      jac1 <- jac
      if (! is.null(jac)) {
      if (class (jac) == "CFunc")
        jac1 <- body(jac)[[2]]
      else if (is.loaded(jac, PACKAGE = dllname, type = "") ||
        is.loaded(jac, PACKAGE = dllname, type = "Fortran"))  {
        jac1 <- getNativeSymbolInfo(jac, PACKAGE = dllname)$address
      } else 
        stop(paste("'jac' not loaded ", jac))
      }
      
      if (class (initfunc) == "CFunc")
        initfunc1 <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        initfunc1 <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else 
        stop(paste("'initfunc' not loaded ", initfunc))

    if (is.null(rpar))
       rpar <- 0.
    if (is.null(ipar))
       ipar <- 0  

    if (is.null(lower)) 
        lower <- rep(-Inf, length(par))
    if (is.null(upper)) 
        upper <- rep(Inf, length(par))
    if (length(lower) != length(par)) 
        stop("length(lower) must be equal to length(par)")
    if (length(upper) != length(par)) 
        stop("length(upper) must be equal to length(par)")

    if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper) 
        stop("lower < upper  is not fulfilled")

    ctrl <- nls.lm.control()
    if (!missing(control)) {
        control <- as.list(control)
        ctrl[names(control)] <- control
    }
    if (length(ctrl[["maxfev"]]) == 0) 
        ctrl[["maxfev"]] <- 100 * (length(unlist(par)) + 1)

    data <- as.matrix(data)  
    ndat <- nrow(data)
    ncoldat <- ncol(data)    
    storage.mode(data) <- "double"  
    out <- .Call("cc_nls_lm", as.double(unlist(par)), fn1, jac1, initfunc1, 
      ctrl, as.double(lower), as.double(upper), as.integer(ndat),
      as.integer(ncoldat), data, as.double(rpar),
      as.integer(ipar), PACKAGE = "ccSolve")
      
    out$hessian <- matrix(out$hessian, nrow = length(unlist(par)))
    names(out$par) <- rownames(out$hessian) <- colnames(out$hessian) <- names(out$diag) <- names(par)
    out$deviance <- sum(out$fvec^2)
    class(out) <- c("ccnls","nls.lm")
    out
}

## to calculate values 
ccfunc.ccnls <- function (fn, data, par, dllname = NULL, initfunc = NULL,
    rpar = NULL, ipar = NULL) {

    if (!is.numeric(par)) 
        stop("`par' must be numeric")

    cl <- attributes(fn)$call
    if (! is.null(cl))
      if (cl != "compile.nls")
        stop ("problem is not compiled with 'compile.nls'")

    if (is.list(fn)) {       
      if (! is.null(Dnames <- attributes(fn)$datanames))
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

      if (class (fn) == "CFunc")
        fn1 <- body(fn)[[2]]
      else if (is.loaded(fn, PACKAGE = dllname, type = "") ||
        is.loaded(fn, PACKAGE = dllname, type = "Fortran"))  {
        fn1 <- getNativeSymbolInfo(fn, PACKAGE = dllname)$address
      } else 
        stop(paste("'fn' not loaded ", fn))

      if (class (initfunc) == "CFunc")
        initfunc1 <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        initfunc1 <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else 
        stop(paste("'initfunc' not loaded ", initfunc))

    if (is.null(rpar))
       rpar <- 0.
    if (is.null(ipar))
       ipar <- 0  

    data <- as.matrix(data)  
    ndat <- nrow(data)
    ncoldat <- ncol(data)    
    storage.mode(data) <- "double"  

    out <- .Call("cc_nls_func", as.double(unlist(par)), fn1, initfunc1, 
      as.integer(ndat), as.integer(ncoldat), data, as.double(rpar),
      as.integer(ipar), PACKAGE = "ccSolve")
      
    out
}
