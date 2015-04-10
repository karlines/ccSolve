## ============================================================================
## function evaluation of compiled code objects
## ============================================================================

ccfunc <- function(fn, ...) {
  Call <- attributes (fn)$call

  if (Call == "compile.nls")
    out <- ccfunc.ccnls(fn, ...)
  else if (Call %in% c("compile.ode", "compile.steady"))
    out <- ccfunc.de(fn, ...)
  else if (Call %in% c("compile.bvp"))
    out <- ccfunc.bvp(fn, ...)
  else if (Call == "compile.multiroot")
    out <- ccfunc.multiroot(fn, ...)
  else if (Call == "compile.optim")
    out <- ccfunc.optim(fn, ...)
  else if (Call %in% c("compile.optimize", "compile.uniroot"))
    out <- ccfunc.optimize(fn, ...)
  else if (Call == "compile.dae")
    out <- ccfunc.dae(fn, ...)
  else  
    stop("No function evaluation possible for fn compiled with ", Call)

  out    
}

ccfunc.multiroot <- function(fn, start, parms = NULL, ...) {
  ccfunc.de(fn, times = 0, y = start, parms = parms, ...)$dy  
}
ccfunc.bvp <- function(fn, yini, x, ...) {
  ccfunc.de(fn, times = x, y = yini, ...)$dy  
}

## ==========================================================================

ccfunc.de <- function (func, times, y, parms, dllname = NULL, initfunc = dllname, 
    rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
    initforc = NULL, fcontrol = NULL) 
{
    cl <- attributes(func)$call
    if (! is.null(cl))
      if (!cl %in% c("compile.ode", "compile.steady", "compile.bvp", "compile.multiroot"))
        stop ("problem is not compiled with 'compile.ode' or 'compile.steady'")

    if (is.list(func))
      func <- func[[1]]
    out <- DLLfunc(func, times, y, parms, dllname, initfunc, 
      rpar, ipar, nout, outnames, forcings, initforc, fcontrol)
    return(out)
}

ccfunc.dae <- function (res, times, y, dy, parms, dllname = NULL, initfunc = dllname, 
    rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
    initforc = NULL, fcontrol = NULL) 
{
    cl <- attributes(res)$call
    if (! is.null(cl))
      if (!cl == "compile.dae")
        stop ("problem is not compiled with 'compile.dae'")

    out <- DLLres(res, times, y, dy, parms, dllname, initfunc, 
      rpar, ipar, nout, outnames, forcings, initforc, fcontrol)
    return(out)
}
  
