\name{ccfunc}
\alias{ccfunc}
\title{
  Utilities for ccSolve.
}
\description{
  Retrieves function values of compiled code problems.
}

\usage{
ccfunc(fn, ...)
}
\arguments{
  \item{fn }{The compiled object whose function value should be retrieved.
  } 
  \item{... }{Arguments required to estimate the function value, as passed to the 
    solution method. 
  } 

}

\note{
  As the compilation of a problem is the time-limiting factor, a compiled object
  can be saved with \link{writeDynLib} and then later loaded with \link{readDynLib}
  from the package \code{inline}.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
# generate a code and print it
String <- " 
 f(1) = 3*x(1) + x(2)*x(3) +2
 f(2) = 2*x(2) + 3
 f(3) = x(3)*(x(1)+1)
"
rfun <- compile.multiroot(String)
code(rfun)

(res <- multiroot(rfun, 1:3))

# function value
ccfunc(rfun, 1:3)
ccfunc(rfun, res$root)
\dontrun{

# save the compiled code for later use - will create two files
save.cc(rfun, file = "test.cc")

# load it, call it rr
rr <- load.cc("test.cc")
multiroot(rr, 1:3)
}
}

\references{
A lot of the code used to create the compiled codes is modified from the R-package
inline.
The modifications relate to the part of function \code{cfunction} 
to create inline code for fortran and fortran 95. 

Oleg Sklyar, Duncan Murdoch, Mike Smith, Dirk Eddelbuettel and Romain Francois (2013). 
inline: Inline C, C++, Fortran function calls from R. R package version
  0.3.13. http://CRAN.R-project.org/package=inline 
}



