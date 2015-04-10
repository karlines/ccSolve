declare.pnames <- function(parms, language, header) {
  header2 <- character()
  if (language == "Fortran")
    lead <- "        "
  else
    lead <- ""  
  if (! is.null(parms))
    if (is.character(parms))
      pnames <- parms
    else
      pnames <- names(parms)
  else
    pnames <- NULL  
  if (! is.null(pnames) & language %in% c("F95", "Fortran")){
     header <- paste(header, 
         "\n", lead,"double precision ", paste(pnames, collapse = ", "), sep ="")
     dol <- "\n"
     for (i in 1:length(pnames))
       dol <- paste(dol,    
           lead, pnames[i]," = x(", i,")\n", collapse = "", sep = " ")
     header2 <- dol
  

  } else if (! is.null(pnames)){ # C
     header <- paste(header, "\n double ", paste(pnames, collapse = ", "),";\n")
     dol <- "\n"
     for (i in 1:length(pnames))
       dol <- paste(dol, " ",
            pnames[i]," = x[", i-1,"];\n", collapse = "", sep = "")
     header2 <- dol
  }
  list(header = header, header2 = header2)
}

# the data - they stay the same during the fitting - put in a module (F95;
# or as global value (C) - refuse to use Fortran 

declare.dnames <- function(d, language, includes, isnls = TRUE) {

  if (language == "Fortran")
    lead <- "        "
  else
    lead <- ""  

  if (! is.null(d))
    if (is.character(d))
      dnames <- d
    else if (is.matrix(d))
      dnames <- colnames(d)
    else  
      dnames <- names(d)
  else
    dnames <- NULL

  nvar <- ncol(d)
  if (is.null(nvar))
    nvar <- length(d)
    
  if (is.null(dnames) & isnls)
    stop("code that uses data without varnames is not supported")
  
  module <- NULL    
  if (! is.null(dnames) & language == "Fortran")
    stop ("cannot use named variables in Fortran - switch to F95 or C")
  
  if (language == "F95"){
     dol <- paste(lead,  "module modnlsdata \n implicit none \n", sep ="")
     dol <- paste(dol, lead, " integer, parameter :: nvar = ", nvar, "\n", sep ="")
     if (! isnls)
       dol <- paste(dol, lead, " integer :: ndata \n", sep ="")
     dol <- paste(dol, 
         lead," double precision, dimension (:), allocatable ::", paste(dnames, collapse = ", "), sep ="")
     dol <- paste(dol, "\n", lead,"end module modnlsdata\n\n", sep ="")
     module <- "modnlsdata"

     dol <- paste(dol, lead, "subroutine initdat(nlsdat, m)\n use modnlsdata\n external nlsdat \n integer m\n\n", sep ="")
     if (! isnls) {
       dol <- paste(dol, lead, " if (m >= 0) ndata = m\n", sep ="")
     }
     for (i in 1:length(dnames)) 
       dol <- paste(dol,lead, " if (allocated(",dnames[i],")) deallocate(",dnames[i],")\n", sep ="")

     dol <- paste(dol,lead, " if (m <= 0) return\n", sep ="")
     for (i in 1:length(dnames)) 
       dol <- paste(dol,lead, " allocate(",dnames[i],"( m))\n", sep ="")
     for (i in 1:length(dnames)) 
       dol <- paste(dol,lead, " call nlsdat(", i, ",",dnames[i]," )\n", sep ="")
     dol <- paste(dol,lead,"\nreturn\nend\n", sep ="")   

  } else  { # C
     dol <- paste("int nvar = ", nvar,";\n", sep = "")
     if (! isnls)
       dol <- paste(dol, " int ndata;\n", sep ="")
     
     dol <- paste(dol, "\ndouble *", paste(dnames, collapse = ", *"),";\n", sep = "")
     dol <- paste(dol,"\n")
     
#     dol <- paste(dol, "extern \"C\" void initdat(void(*));\n") 
    
     dol <- paste(dol, "void initdat(void (* nlsdat)(int *, double *), int *m){\n int i;\n", sep = "")    
       dol <- paste(dol, " if (*m <= 0) return;\n", sep ="")

     if (! isnls)
       dol <- paste(dol, " ndata = *m;\n", sep ="")
     for (i in 1:length(dnames)) {
       dol <- paste(dol, dnames[i]," =(double *)R_alloc(*m, sizeof(double));\n");     
       dol <- paste(dol, "i = ", i,"; nlsdat(&i,", dnames[i] , ");\n", sep = " ")
     }  
     dol <- paste(dol, "}\n", sep = "")
       
  }
  list(module  = module, includes = paste(includes, "\n", dol, sep = ""))
}

# =============================================================================
#
# Compiler functions for nonlinear least squares
#
# =============================================================================
#ccDNase <- compile.nls(func = "y = 1/(1 + dexp((xmid - dlog(conc))/scal))",
#  parms = c(xmid = 0, scal = 1), data = DNase1)  

compile.nls <- function(func, jacfunc = NULL, data, par = NULL, 
  declaration = character(), includes = character(), language = "F95", ...) {

   DD <- declare.ynames(par, language, declaration, "x", FALSE)
   if(any(mapply(data[1,], FUN = is.factor))) stop ("data should not contain factors")

  dn <- declare.dnames (data, language, includes)
       
  out <- compileDE(func, Data = data, jacfunc = jacfunc,
     header = DD$header, header2 = DD$header2, module = dn$module, 
     includes = dn$includes, language = language, type = 5, ...)

# make sure that a direct call to out$func does not break R: (will not work anyway)
  fn <- out$func
  body(fn)$ndat <- 0
  out$func@.Data <- fn
  
  attr(out, "call") <- "compile.nls"
  attr(out, "Dnames") <- getnames(data)
  out   
}

