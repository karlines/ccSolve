\name{ccSolve-package}
\alias{ccSolve-package}
\alias{ccSolve}
\docType{package}
\title{
  creating compiled code problems for solver packages deSolve, bvpSolve, rootSolve, 
  and for nls, optim, uniroot.
}
\description{
  Functions to solve numerical problems using inline 
   compiled code, written as the body of Fortran, F95, C (and C++) code.
   
  Can be used with solvers from packages deSolve, bvpSolve, deTestSet, 
  and rootSolve.
}
\details{
  \tabular{ll}{
    Package: \tab ccSolve\cr
    Type:    \tab Package\cr
    License: \tab GNU Public License 2 or above\cr
  }

  The numerical problem (e.g. system of differential equations,
    function whose root is to be found),
  is written as a string that contains the body of Fortran, F95 or C code.    
  
  The functions from package \code{ccSolve} then complete this code body,
  by adding headers and footers, and sometimes initialisation functions.
  The completed code is then compiled and the resulting address stored, with 
  the code in a list of class \code{CFunc}, as from the R-package \code{inline}.   
  
  The solvers from the packages \code{deSolve}, \code{rootSolve}, \code{bvpSolve}
  and \code{deTestSet} have been adapted to accept problems specified by \code{CFunc}
  objects.
  
  Package \code{ccSolve} contains slightly modified versions of the R-functions
  \code{optim}, \code{optimize} and \link{uniroot} to work with compiled code
  (\link{ccoptim}, \link{ccoptimize}, \link{ccuniroot}).
  
  Finally, a modified version of \link{nls.lm}, (from R-package \code{minpack.lm})
  is called \link{ccnls}
}

\author{
  Karline Soetaert (Maintainer)
}
\references{
R Core Team (2014). R: A language and environment for statistical computing. R
  Foundation for Statistical Computing, Vienna, Austria. URL
  http://www.R-project.org/

Oleg Sklyar, Duncan Murdoch, Mike Smith, Dirk Eddelbuettel and Romain Francois (2013). inline: Inline C, C++, Fortran function calls from R. R package version 0.3.13. http://CRAN.R-project.org/package=inline 

Timur V. Elzhov, Katharine M. Mullen, Andrej-Nikolai Spiess and Ben Bolker (2013). minpack.lm: R interface to the
  Levenberg-Marquardt nonlinear least-squares algorithm found in MINPACK, plus support for bounds. R package version 1.1-8.
  http://CRAN.R-project.org/package=minpack.lm

Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010).  Solving Differential Equations in R: Package deSolve Journal of
  Statistical Software, 33(9), 1--25.  URL http://www.jstatsoft.org/v33/i09/.

Francesca Mazzia, Jeff R. Cash, Karline Soetaert (2014).  Solving boundary value problems in the open source software R: Package bvpSolve.  Opuscula mathematica, 34(2), 387--403.  URL http://dx.doi.org/10.7494/OpMath.2014.34.2.387

Soetaert K. (2009).  rootSolve: Nonlinear root finding, equilibrium and steady-state analysis of ordinary differential equations.  R-package version 1.6

Karline Soetaert, Jeff Cash and Francesca Mazzia (2014). deTestSet: Testset for differential equations. R package version 1.1.2
}
\seealso{
  \code{\link{ode}} solver from the deSolve package.

  \code{\link{bvptwp}},  \code{\link{bvpcol}}, solvers from the bvpSolve package.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(compile.bvp)

example(compile.ode)

example(compile.dae)

example(compile.multiroot)

example(utilities)
}
}

\keyword{ package }
