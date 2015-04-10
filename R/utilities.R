## ============================================================================
## Small utilities - not exported
## ============================================================================

getnames <- function (x) 
  if(is.matrix(x)) return(colnames(x)) else return(names(x))

isvalid <- function(f) {        # check if a function points to compiled code
  is.character(f) | class(f) == "CFunc"
}

## ============================================================================
## String functions
## ============================================================================

# ---------------------------------------------------------------------------
# if FORTRAN text longer than 72 characters - split it
# ---------------------------------------------------------------------------

toFortran <- function (text) {
   if ((cl <- nchar(text)) >= 72) {
     fstring <- substr(text, 72, cl)
     text <- substr(text, 1, 71)
     while ((cf <- nchar(fstring)) > 66) {
        text <- paste(text, "\n     &", substr(fstring, 1, 66), sep = "")
        fstring <- substr(fstring, 67, cf)
     } 
     if (cf > 0)   
        text <- paste(text, "\n     &", fstring, sep = "")
        text <- paste(text, "\n")   
     } 
   return (text)
} 
toF95 <- function (text) {
   if ((cl <- nchar(text)) >= 100) {
     fstring <- text
     text <- ""
     while ((cf <- nchar(fstring)) > 99) {
       poscomma <- gregexpr(',',fstring)[[1]]
       ii <- poscomma[poscomma < 100]
       il <- ii[length(ii)]
       text <- paste(text, substr(fstring, 1, il), "&\n",sep = "")
       fstring <- substr(fstring, il+1, cl)
     } 
     if (cf > 0)   
        text <- paste(text, fstring, sep = "")
        text <- paste(text, "\n")   
     } 
   return (text)
} 

# ---------------------------------------------------------------------------
# get and print dimensionality of object
# A <- list(B = 1, C= 1:2, D = diag(2))  => "B,C(2),D(2,2)"
# ---------------------------------------------------------------------------

getdims <- function (pars, language = "F95") {
  dimfun.f <- function(x) {
    nn <- ""
    if (!is.vector(x)) 
      paste(nn, "(",paste(dim(x), collapse = ","),")", sep = "") 
    else if (length(x) == 1) 
      nn 
    else 
      paste(nn,"(",length(x),")",sep="")
  }                             

  dimfun.C <- function(x) {
    nn <- ""
    if (!is.vector(x)) 
      paste(nn, "[",paste(dim(x), collapse = ","),"]", sep = "") 
    else if (length(x) == 1) 
      nn 
    else 
      paste(nn,"[",length(x),"]",sep="")
  }                             
  if (language == "C")
    return(paste(names(pars),lapply(pars, FUN = dimfun.C),sep=""))
  else
    return(paste(names(pars),lapply(pars, FUN = dimfun.f),sep=""))
}

# ---------------------------------------------------------------------------

create.ynamesc <- function (ynames) {
    npar <- length(ynames)
    head <- paste("#define ", ynames, " y[", 0:(npar-1), "]\n", sep = "", collapse = "")
    head
}

# ---------------------------------------------------------------------------

declare.ynames <- function(y, language, header, label = "y", dy = TRUE) {

  if (language == "Fortran")
    lead <- "        "
  else
    lead <- ""  

  tail <- character() 
  header2 <- character()
  if (! is.null(y))
    if (is.character(y))
      ynames <- y
    else
      ynames <- names(y)
  else
    ynames <- NULL  

  if (! is.null(ynames) & language %in% c("F95", "Fortran")){
     dynames <- paste("d", ynames, sep="")
     header <- paste(header, 
         "\n", lead, " double precision ", paste(ynames, collapse = ", "))
     if (dy) 
       header <- paste(header, 
         "\n", lead, " double precision ", paste(dynames, collapse = ", "),"\n")

     dol <- "\n"
     for (i in 1:length(ynames))
       dol <- paste(dol,    
           lead , " ", ynames[i]," = ", label, "(", i,")\n", collapse = "", sep = "")
     header2 <- dol
  
     if (dy) {
      tail <- "\n"
      for (i in 1:length(dynames))
        tail <- paste(tail,    
            lead, " f(" , i,") = ",dynames[i],"\n", collapse = "", sep = "")
     }       
  } else if (! is.null(ynames)){ # C
     dynames <- paste("d", ynames, sep="")
     header <- paste(header, "\n double ", paste(ynames, collapse = ", "),";\n")
     if (dy)
       header <- paste(header, "\n double ", paste(dynames, collapse = ", "),";\n")

     dol <- "\n"
     for (i in 1:length(ynames))
       dol <- paste(dol, 
           lead, " " , ynames[i]," = ",label,"[", i-1,"];\n", collapse = "", sep = "")
     header2 <- dol
  
     if (dy) {
      tail <- "\n"
      for (i in 1:length(dynames))
        tail <- paste(tail,    
           lead, " f[" , i-1,"] = ",dynames[i],";\n", collapse = "", sep = "")
     }
  }
  list(header = header, header2 = header2, tail = tail)
}

# ---------------------------------------------------------------------------
# Check input 
# ---------------------------------------------------------------------------

checkinput <- function (body, header, fun) {
  if (! is.character(body))
    stop ("'body' should be a character string in ", fun)
    
  if (! is.null(header))
    if (! is.character(header))
    stop ("'header' should be a character string or NULL, in ", fun)
} 
  
