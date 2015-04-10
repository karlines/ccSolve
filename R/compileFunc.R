# =============================================================================
#
# Creates the initialisers for parameter and forcings 
#
# =============================================================================

create.initfunc <- function (par, language = "F95", pname = "xcbpar") 
  create.init (par, pname, language, subname = "initpar", parms = "parms")
  
create.initforc <- function (forc, language = "F95", fname = "xcbforc") 
  create.init (forc, fname, language, subname = "initforc", parms = "forcs")

# -----------------------------------------------------------------------------
# main initialiser function
# -----------------------------------------------------------------------------

create.init <- function (pars, pname, language, subname, parms, ...) {

  if (! class(pars) %in% c("numeric", "list", "character"))
    stop ("'pars' should be a vector or a list ")

  if (is.list(pars))  # parameters with a dimension
    npar <- sum(unlist(lapply(pars, FUN = length)))
  else
    npar <- length(pars)

  if (is.character(pars))
    pnames <- pars
  else
    pnames <- names(pars)
  
  if (language %in% c("F95", "Fortran")) {  # common blocks in subroutines
     lead <- ifelse (language == "Fortran", "       ","") 
     code <- paste(lead, "SUBROUTINE ", subname, 
    "(de",parms,")\n", lead, "EXTERNAL de",parms,
                 "\n", lead, "DOUBLE PRECISION ",parms,"(",npar,")\n", sep="")
 
    code <- paste(code, 
                  lead, "COMMON / ", pname, " / ",parms,
                 "\n", lead, "CALL de",parms,"(", npar,", ",parms,
             ")\n \n", lead, "END\n", sep = "")

    np <- paste(pnames, collapse = ", ")
    if (is.list(pars)) 
      dimpnames <- paste(getdims(pars), collapse = ",")
    else
      dimpnames <- np

    
    if (language == "Fortran") {
      head <- paste("\n", lead, "DOUBLE PRECISION ", dimpnames , "\n")
      head <- toFortran(head)
      h2 <- paste( lead, "COMMON / ", pname, " / ", np)
      h2 <- toFortran(h2)
      head <- paste(head, h2)
    } else {
    head <- paste("\n", lead, "double precision ", dimpnames , "\n")
    head <- toF95(head)
    h2 <- paste( lead, "common / ", pname, " / ", np)
    h2 <- toF95(h2)
    head <- paste(head, h2)

    }
       
  } else if (length(pnames) == npar) {    
  
  # C-code : main variables with #define statements in front of subroutines
    
    code <- paste("static double ", parms, "[",npar, "];\n", sep = "")
    defines <- paste("#define ", pnames, " ", 
      parms,"[", 0:(npar-1), "]\n", sep = "", collapse = "")

    code <- paste(code, defines, sep = "")
    if (language == "C++")
      code <- paste(code, "extern \"C\" void", subname, "(void(*));\n") 
    
    code <- paste(code, "\nvoid ", subname, 
      "(void (* ode",parms,")(int *, double *)){\n    int N =", npar, 
              ";\n    ode",parms,"(&N, ",parms,");\n}", sep = "")
    head <- ""          
    
  } else { # multidimensional parameters - par initialisation in odeparms
  
    code <- paste("static double ", parms, "[",npar, "];\n", sep = "")
    if (is.list(pars)) 
      dimpnames <- paste(getdims(pars, language = "C"), collapse = ",")
    else
      dimpnames <- np

    defines <- paste("\n", "double ", dimpnames , ";\n", sep = "", collapse = ", ")

    code <- paste(code, defines, "\nvoid ", subname, 
      "(void (* ode",parms,")(int *, double *)){\n int i;\n    int N =", npar, 
              ";\n    ode",parms,"(&N, ",parms,");", sep = "")
    cdim <- unlist(lapply(pars, FUN = length))
    cl <- c(0, cumsum(cdim))

    for (i  in 1:length(pnames))  {
      if (cdim[i] > 1)
        if (i > 1)
          code <- paste(code, "\nfor (i = ", cl[i], "; i < ", cl[i+1], "; i++)\n ",
          pnames[i], "[i-",cl[i],"] = parms[i];",sep = "")
        else
         code <- paste(code, "\nfor (i = ", cl[i], "; i < ", cl[i+1], "; i++)\n ",
          pnames[i], "[i] = parms[i];",sep = "")
           
      else
        code <- paste(code, "\n", pnames[i], " = parms[",cl[i],"];",sep = "")
          
    }
    code <- paste(code, ";\n}")          
    head <- ""          
  
  }
  return(list(code = code, head = head))
} 

# =============================================================================
#
# Compiler functions for bvp, ode, pde, root
#
# =============================================================================

compile.bvp <- function(func, jacfunc = NULL, bound = NULL, jacbound = NULL, 
  parms = NULL, yini = NULL, forcings = NULL, outnames = NULL, 
  declaration = character(), includes = character(), language = "F95", ...) {

  DD <- declare.ynames (yini, language, declaration)
       
  out <- compileDE(func, jacfunc = jacfunc, bound = bound, jacbound = jacbound, 
     parms = parms, forcings = forcings, outnames = outnames, 
     header = DD$header, header2 = DD$header2, tail = DD$tail,
     includes = includes, language = language, type = 2, ...)
  attr(out, "call") <- "compile.bvp"
  out   
}

# =============================================================================

compile.ode <- function(func, jacfunc = NULL, rootfunc = NULL, eventfunc = NULL, 
  parms = NULL, y = NULL, forcings = NULL, outnames = NULL, 
  declaration = character(), includes = character(), language = "F95", ...) {

  DD <- declare.ynames (y, language, declaration)

  out <- compileDE(func = func, jacfunc = jacfunc, 
     bound = NULL, jacbound = NULL, rootfunc = rootfunc, eventfunc = eventfunc,
     parms = parms, forcings = forcings, outnames = outnames, 
     header = DD$header, header2 = DD$header2, tail = DD$tail, type = 1,
     includes = includes, language = language, ...)

  attr(out, "call") <- "compile.ode"

  out   
}

compile.steady <- function(func, jacfunc = NULL, rootfunc = NULL, eventfunc = NULL,
      parms = NULL, y = NULL, forcings = NULL, outnames = NULL, 
      declaration = character(), includes = character(), language = "F95", ...) {
  out <- compile.ode(func, jacfunc, rootfunc, eventfunc,
      parms, y, forcings, outnames, declaration, includes, language,...)
  attr(out, "call") <- "compile.steady"

  out   
}  

# =============================================================================

compile.multiroot <- function(func, jacfunc = NULL, parms = NULL, x = NULL,
  declaration = character(), includes = character(), language = "F95", ...) {

  DD <- declare.ynames (x, language, declaration, "x", FALSE)

  out <- compileDE(func = func, jacfunc = jacfunc, parms = parms, 
     header = DD$header, header2 = DD$header2, includes = includes, 
     language = language, type = 4, ...)
  attr(out, "call") <- "compile.multiroot"
  out   
}

# =============================================================================

compile.integrate <- function(func, declaration = character(), includes = character(),      
    language = "F95", ...) {

   out <- compileDE(func, header = declaration, 
     includes = includes, language = language, type = 7, ...)

  attr(out, "call") <- "compile.integrate"
  out   

}    

# =============================================================================

compile.optim <- function(func, jacfunc = NULL, data = NULL, par = NULL, 
  declaration = character(), includes = character(), language = "F95", ...) {
  if (is.null(data)) {
    DD <- declare.ynames(par, language, declaration, "x", FALSE)
    out <- compileDE(func = func, jacfunc = jacfunc, parms = NULL, 
          header = DD$header, header2 = DD$header2, 
          includes = includes, language = language, 
          type = 3, ...)
  } else {
   DD <- declare.ynames(par, language, declaration, "x", FALSE)
   if(any(mapply(data[1,], FUN = is.factor))) stop ("data should not contain factors")

   dn <- declare.dnames (data, language, includes, FALSE)
       
   out <- compileDE(func, Data = data, jacfunc = jacfunc,
     header = DD$header, header2 = DD$header2, module = dn$module, 
     includes = dn$includes, language = language, type = 3, ...)
  }
  attr(out, "call") <- "compile.optim"
  out   
}

# =============================================================================

compile.uniroot <- function(func,  
  declaration = character(), includes = character(), language = "F95", ...) {
  out <- compileDE(func = func, jacfunc = NULL, parms = NULL, 
     header = declaration, includes = includes, language = language, 
     type = 6, ...)
  attr(out, "call") <- "compile.uniroot"
  out   
}

compile.optimize <- compile.optimise <- function(func,  
  declaration = character(), includes = character(), language = "F95", ...) {
  out <- compileDE(func = func, jacfunc = NULL, parms = NULL, 
     header = declaration, includes = includes, language = language, 
     type = 6, ...)
  attr(out, "call") <- "compile.optimize"
  out   
}
# =============================================================================

compile.dae <- function(res, rootfunc = NULL, eventfunc = NULL,  
  parms = NULL,  y = NULL, forcings = NULL, outnames = NULL,
  declaration = character(), includes = character(), language = "F95", ...) {
  DD <- declare.ynames (y, language, declaration)

  out <- compileDE(res = res, 
     bound = NULL, jacbound = NULL, rootfunc = rootfunc, eventfunc = eventfunc,
     parms = parms, forcings = forcings, outnames = outnames, 
     header = DD$header, header2 = DD$header2, tail = DD$tail,
     includes = includes, language = language, ...)
  attr(out, "call") <- "compile.dae"
  out   
}


# =============================================================================

compileDE <- function(func = NULL, jacfunc = NULL, bound = NULL, jacbound = NULL, 
  rootfunc = NULL, res = NULL, eventfunc = NULL,
  parms = NULL, forcings = NULL, Data = NULL, outnames = NULL, module = NULL,
  header = NULL, header2 = NULL, tail = NULL, language = "F95", includes = "", 
  type = 1, ...) {

  if (language == "Fortran") lead <- "      "
  else lead <- ""

  if (language %in% c("Fortran", "F95")) {
    convention <- ".Fortran"  
  } else if (language == "C") {
    convention <- ".C"  
  } else
    stop ("language", language, " not (yet) supported")

  if (! is.null(parms)) {
    F <- create.initfunc(parms, language )
    includes <- paste(includes, "\n", F$code)
    header <- paste(F$head, "\n", header, "\n")
  } 
  
  if (! is.null(forcings)) {
    F <- create.initforc(forcings, language )
    includes <- paste(includes, "\n",F$code)
    header <- paste(F$head, "\n", header, "\n")
  } 
  
  if (is.null (func) & is.null (res))
    stop ("'func' and 'res' cannot both be NULL")
  if (!is.null (func) & !is.null (res))
    stop ("'func' and 'res' cannot both be given")
  
  outhead <- ""
  outtail <- ""
  if (! is.null(outnames)) {
    nout <- length(outnames)
    allnames <- paste(outnames, collapse = ",")
    if (language %in% c("F95", "Fortran")) {
      outhead <- paste(lead, "double precision", allnames)
      if (language == "Fortran")
        outhead <- toFortran(outhead)
      else
        outhead <- toF95(outhead)
      
      header <- paste(header, outhead, "\n")  
      outhead <- paste("\n",lead, "if (ipar(1) < ", nout, ") call rexit('nout should be >= ",nout," ')\n")
      outtail <- paste("rpar(", 1:nout, ") = ", outnames, "\n", sep = "", collapse = "")
    } else {
      outhead <- paste("double ", allnames, ";\n")
      header <- paste(header, outhead, "\n")  

      outhead <- paste("char str[5]; str[0]='n'; str[1]='o'; str[2]='u'; str[3]='t'; str[4]='!';")
      outhead <- paste(outhead,"if (ipar[0] < ", nout, ")\n", "error(str);\n")
      outtail <- paste("rpar[", 0:(nout-1), "] = ", outnames, ";\n", sep = "", collapse = "")
 
    }
  }
  
  if (! is.null(func)) {
    F <- create.func(paste(outhead,header2, func,tail,outtail), header, language, 
      convention, type = type)
    fnames <- "func"
  } else {
    if (! is.null(jacfunc) )
      stop ("cannot combine 'res' with 'jacfunc'")  
    F <- create.res(paste(outhead,res,outtail), header, language, convention)
    fnames <- "res"
  }  
  ii <- 1
  if (! is.null(jacfunc)) {
    F2 <- create.jacfunc(jacfunc, header, language, convention, type = type)
      
    fnames <- c(fnames, "jacfunc")
    ii <- c(ii, 2)
  } else
    F2 <- NULL
  
  if (! is.null(bound)) {
    F3 <- create.bound(bound, header, language, convention)
    fnames <- c(fnames, "bound")
    ii <- c(ii, 3)
  } else
    F3 <- NULL

  if (! is.null(jacbound)) {
    F4 <- create.jacbound(jacbound, header, language, convention)
    fnames <- c(fnames, "jacbound")
    ii <- c(ii, 4)
  } else 
    F4 <- NULL
    
  if (! is.null(rootfunc)) {
    F5 <- create.rootfunc(rootfunc, header, language, convention)
    fnames <- c(fnames, "rootfunc")
    ii <- c(ii, 5)
  } else 
    F5 <- NULL

  if (! is.null(eventfunc)) {
    F6 <- create.event(eventfunc, header, language, convention)
    fnames <- c(fnames, "eventfunc")
    ii <- c(ii, 6)
  } else 
    F6 <- NULL
  
  if (! is.null(F2) | ! is.null(F3) | ! is.null(F4) | ! is.null(F5) | ! is.null(F6)) {
    F$sig  <- list(F$sig, F2$sig, F3$sig, F4$sig, F5$sig, F6$sig)       [ii]
    F$body <- list(F$body, F2$body, F3$body, F4$body, F5$body, F6$body) [ii]
    F$dim  <- list(F$dim, F2$dim, F3$dim, F4$dim, F5$dim, F6$dim)       [ii]
  } else 
    F$sig <- list(F$sig)  # this to make sure that the subroutine has correct name!
  names(F$sig) <- fnames   
 
  CF <- cfunction(sig = F$sig, body = F$body,  dim = F$dim, implicit = "none", 
    includes = includes, convention = convention, language = language, module = module, ...) 
  if (! is.null(parms) | ! is.null(forcings) | ! is.null(Data)) {
    DLL <- getDynLib(CF)   # THIS IF cfunction from inline can be used
#    DLL <- dyn.load(attributes(CF)$DLL)
    
    if (! is.null(parms)) {
      res <- new("CFunc", code = includes)
      fn <- function() NULL
      body <- list()
      body[[1]] <- ""
      body[[2]] <- getNativeSymbolInfo("initpar", DLL)$address
      body(fn) <- body 
      res@.Data <- fn
      remove(list = c("body", "fn"))
      fnames <- c(fnames, "initfunc")
      if (class(CF) == "CFunc") {
        CFL <- CF
        CF <- new("CFuncList", list())
        CF[[1]] <- CFL
      }  
      CF$initfunc <- res
    }
    if (! is.null(forcings)) {
      res <- new("CFunc", code = includes)
      fn <- function() NULL
      body <- list()
      body[[1]] <- ""
      body[[2]] <- getNativeSymbolInfo("initforc", DLL)$address
      body(fn) <- body 
      res@.Data <- fn
      remove(list = c("body", "fn"))
      fnames <- c(fnames, "initforc")
      if (class(CF) == "CFunc") {
        CFL <- CF
        CF <- new("CFuncList", list())
        CF[[1]] <- CFL
      }  
      CF$initforc <- res
    } 
    if (! is.null(Data)) {
      res <- new("CFunc", code = includes)
      fn <- function() NULL
      body <- list()
      body[[1]] <- ""
      body[[2]] <- getNativeSymbolInfo("initdat", DLL)$address
      body(fn) <- body 
      res@.Data <- fn
      remove(list = c("body", "fn"))
      if (is.null(CF$initfunc))
        fnames <- c(fnames, "initfunc")           # this should occur at initialisation; there will be no initpar
      if (class(CF) == "CFunc") {
        CFL <- CF
        CF <- new("CFuncList", list())
        CF[[1]] <- CFL
      }  
      CF$initfunc <- res
    }

  }
  if (length(CF) > 1) {
     names(CF) <- fnames 
  }

  return(CF)
}

# =============================================================================
# Suitable derivative function subroutine
#
# func(n, x, y, f, rpar, ipar) for bvp
# func(n, t, y, f, rpar, ipar) for ivp
# func(n, t, x, f, rpar, ipar) for root
#
# integer n, ipar(*)
# double precision x, t, y(*), x(*), F(*), rpar(*)
# =============================================================================

create.func <- function (body, header, language, convention, 
  type = 1) {

  checkinput(body, header, "func")

  body <- paste("\n",header, "\n" , body, "\n")
   dim <- c("", "", rep("(*)", 4))

   if (type == 2)                       # bvp
     sig <- signature(n = "integer", x = "numeric", 
       y = "numeric", f = "numeric", rpar = "numeric", ipar = "integer")
   else if (type == 1)                  # ode
     sig <- signature(n = "integer", t = "numeric", 
       y = "numeric", f = "numeric", rpar = "numeric", ipar = "integer")
   else if (type == 4)                  # root
     sig <- signature(n = "integer", t = "numeric", 
       x = "numeric", f = "numeric", rpar = "numeric", ipar = "integer")
   else if (type == 3) {                # optim
     sig <- signature(n = "integer", x = "numeric", 
       f = "numeric", rpar = "numeric", ipar = "integer")
    dim <- c("", "(n)", "", rep("(*)", 2))
   } else if (type == 5) {                # nls
     sig <- signature(n = "integer", ndat = "integer", x = "numeric", 
       f = "numeric", rpar = "numeric", ipar = "integer")
    dim <- c("", "", "(n)", "(ndat)", rep("(*)", 2))
   } else if (type == 6) {                # uniroot
     sig <- signature(x = "numeric", 
       f = "numeric", rpar = "numeric", ipar = "integer")
     dim <- c("", "", rep("(*)", 2))
   }    else if (type == 7) {                # integrate
     sig <- signature(n = "integer", x = "numeric", 
       f = "numeric", rpar = "numeric", ipar = "integer")
    dim <- c("", rep("(n)", 2), rep("(*)", 2))
   } else
     stop ("'type' not known")
  
  list(body = body, sig = sig, dim = dim)    
}

# -----------------------------------------------------------------------------

create.res <- function (body, header, language, convention) {

  checkinput(body, header, "func")

  body <- paste("\n",header, "\n" , body, "\n")

  sig <- signature(t = "numeric", 
      y = "numeric", dy = "numeric", cj = "numeric", r = "numeric", ires = "integer",
       rpar = "numeric", ipar = "integer")
  dim <- c("", rep("(*)", 2), "", "(*)", "", rep("(*)", 2))

  list(body = body, sig = sig, dim = dim)    
}

# =============================================================================
# Suitable jacobian function subroutine
# jacfunc(n, x, y, df, rpar, ipar) 
# =============================================================================

create.jacfunc <- function(body, header, language, convention, type = 1) {

  checkinput(body, header, "jacfunc")
  dim <- c("", "", "(*)", "", "", "(nrowpd,*)","", "(*)","(*)")

  if (type == 2) {     #bvp
    sig <- signature(n = "integer", x = "numeric", 
       y = "numeric", df = "numeric", 
       rpar = "numeric", ipar = "integer")
    dim <- c("", "", "(*)", "(n,*)","(*)","(*)")
  } else if (type == 1) #ode
    sig <- signature(n = "integer", t = "numeric", 
      y = "numeric", ml = "integer", mu = "integer", 
      df = "numeric", nrowpd = "integer", 
      rpar = "numeric", ipar = "integer")
  else if (type == 4) #root
    sig <- signature(n = "integer", t = "numeric", 
      x = "numeric", ml = "integer", mu = "integer", 
      df = "numeric", nrowpd = "integer", 
      rpar = "numeric", ipar = "integer")
  else if (type == 5) {                # nls   n, ndat, ldfjac, x, df, rpar, ipar
     sig <- signature(n = "integer", ndat = "integer", ldfjac = "integer",
       x = "numeric", df = "numeric", rpar = "numeric", ipar = "integer")
    dim <- c("", "", "", "(n)", "(n,ldfjac)", rep("(*)", 2))
  } else if (type == 3) { #min
    sig <- signature(n = "integer", x = "numeric", df = "numeric", 
      rpar = "numeric", ipar = "integer")
    dim <- c("", "(*)", "(*)", "(*)", "(*)")
  }
  if (language == "Fortran") {
    head <- "\n        integer ix, jx\n"
    if (type == 5) {
      head <- "\n        integer ix, jx\n"
      head <- paste("\n", header, head,
            "\n        DO ix = 1, n\n          DO jx = 1, ldfjac\n            DF(ix,jx) = 0.D0\n          ENDDO\n        ENDDO\n\n")
    } else if (type != 3){
      head <- "\n        integer ix, jx\n"
      head <- paste("\n", header, head,
            "\n        DO ix = 1, n\n          DO jx = 1, n\n            DF(ix,jx) = 0.D0\n          ENDDO\n        ENDDO\n\n")
    } else {
      head <- "\n        integer ix\n"
      head <- paste("\n", header, head,
            "\n        DO ix = 1, N\n          DF(ix) = 0.D0\n        ENDDO\n\n")
    }        
    body <- paste(head, body, "\n")
   } else if (language == "F95")  {   
    if (type == 5) {
      head <- "\n integer ix, jx\n do ix = 1, n\n"
      head <- paste("\n", header, head," do jx = 1, ldfjac\n   df(ix,jx) = 0.d0\n  enddo\n enddo\n")
    } else if (type != 3) {
      head <- "\n integer ix, jx\n do ix = 1, n\n"
      head <- paste("\n", header, head," do jx = 1, n\n   df(ix,jx) = 0.d0\n  enddo\n enddo\n")
    } else {
      head <- "\n integer ix\n do ix = 1, n\n"
      head <- paste("\n", header, head," df(ix) = 0.d0\n enddo\n")
    }
    body <- paste(head, body, "\n")
  } else {
    if (type == 5) 
      head <- "\n int ix;\n for(ix = 0; ix < *n * *ldfjac; ix++)\n"
    else if (type != 3)
      head <- "\n int ix;\n for(ix = 0; ix < *n * *n; ix++)\n"
    else
      head <- "\n int ix;\n for(ix = 0; ix < *n; ix++)\n"
    
    head <- paste("\n", header, head," df[ix] = 0.;\n")
    body <- paste(head, body, "\n")
    body <- paste(header,"\n", body, "\n")
  }
  list(body = body, sig = sig, dim = dim)    

}

# =============================================================================
# Suitable boundary function subroutine
# gsub(i, n, x, y, g, rpar, ipar) 
# =============================================================================

create.bound <- function(body, header = NULL, language = "F95", convention) {

  checkinput(body, header, "bound")

  sig <- signature(i = "integer", n = "integer", 
     y = "numeric", g = "numeric", rpar = "numeric", ipar = "integer")

  if (language %in% c("Fortran", "F95")) 
    dim <- c("", "", "(*)", "", rep("(*)", 2))
  else
    dim <- ""
  body <- paste("\n", header, "\n" , body, "\n")

  list(body = body, sig = sig, dim = dim)    
}

# =============================================================================
# Suitable jacobian of boundary function subroutine
# dgsub(i, n, x, y, dg, rpar, ipar) 
# =============================================================================

create.jacbound <- function(body, header = NULL, language = "F95", convention) {

  checkinput(body, header, "jacbound")

  sig <- signature(i = "integer", n = "integer", 
     y = "numeric", dg = "numeric", rpar = "numeric", ipar = "integer")

  dim <- c("", "", rep("(*)", 4))

  if (language == "Fortran") { 
    body <- paste(header, "\n        INTEGER jx\n        DO jx = 1, N\n          DG(jx) = 0.D0\n        ENDDO\n", body , "\n") 
  } else if (language == "F95") {  # \n implicit none\n
    body <- paste(header, "\n integer jx\n do jx = 1, n\n dg(jx) = 0.d0\n enddo\n", body, "\n") 

  } else {
    body <- paste(header, "\n  int ix;\n for(ix = 0; ix < *n; ix++)\n dg[ix] = 0.;\n", body, "\n")
    dim <- ""
  }

  list(body = body, sig = sig, dim = dim)    

}

# =============================================================================
# rootfunction
# rootfunc(n, t, y, nroot, root, rpar, ipar) 
# =============================================================================

create.rootfunc <- function(body, header = NULL, language = "F95", convention) {

  checkinput(body, header, "rootfunc")

  sig <- signature(n = "integer", t = "numeric", 
     y = "numeric", nroot = "integer", root = "numeric", rpar = "numeric", ipar = "integer")

  if (language %in% c("Fortran", "F95")) 
    dim <- c("", "", "(*)", "", "(*)", rep("(*)", 2))
  else
    dim <- ""
  body <- paste("\n", header, "\n" , body, "\n")

  list(body = body, sig = sig, dim = dim)    
}

# =============================================================================
# eventfunction
# eventfunc(n, x, y) 
# =============================================================================

create.event <- function(body, header = NULL, language = "F95", convention) {

  checkinput(body, header, "eventfunc")

  sig <- signature(n = "integer", t = "numeric", 
     y = "numeric")

  if (language %in% c("Fortran", "F95")) 
    dim <- c("", "", "(*)")
  else
    dim <- ""
  body <- paste("\n", header, "\n" , body, "\n")

  list(body = body, sig = sig, dim = dim)    
}


