## =============================================================================
## Example 1.
## Find the 4th eigenvalue of Mathieu's equation:
## y''+(lam-10cos2t)y=0   on the interval [0,pi]
## y(0)=1, y'(0)=0  and y'(pi)=0
##
## 2nd order problem is rewritten as:
## dy=y2
## dy2= -(y(3)-10cos(2t))*y
## dy3 = 0 
## =============================================================================

require(bvpSolve)

## ------------------------
## Problem specifications
## ------------------------

x      <- seq(0, pi, by = 0.01)
init   <- c(y = 1,dy = 0, lambda = NA)
xguess <- c(0, 1, 2*pi)
yguess <- matrix(nr = 3, rep(1, 9))
rownames(yguess) <- c("y", "dy", "lambda")

## ------------------------
## Inline fortran derivative
## ------------------------

Mathieu <- "
       f[0] = y[1] ;
       f[1] = -(y[2]-10*cos(2* *x))*y[0] ;
       f[2] = 0. ;
"
mathieu2 <- compile.bvp(Mathieu, language = "C")   # compile it

## ------------------------
## Solution
## ------------------------

print(system.time(                  # run it
sol <- bvptwp(yini = init, yend = c(NA, 0, NA), x = x,
        func = mathieu2, xguess = xguess, yguess = yguess)
))
plot(sol)

## =============================================================================
## PROBLEM measels
## Models the spread of measels in three equations
## U. M. Ascher, R. M. R. Mattheij, and R. D. Russell. Numerical Solution of
## Boundary Value Problems for Ordinary Differential Equations. Prentice{Hall,
## Englewood Cliffs, NJ, USA, 1988.
## =============================================================================

## ------------------------
## Problem specifications
## ------------------------

x <- seq (0, 1, by = 0.01)
yguess <- matrix(ncol = length(x), nrow = 6, data = 1)
rownames(yguess) <- paste("y", 1:6, sep = "")

## ------------------------
## Inline C
## ------------------------

parms <- c(mu = 0.02, lam = 0.0279, vv = 0.1)

Measel <- "
   double pi = 3.141593  ;
   bet = 1575.0*(1.+cos(2*pi*x[0]));
   f[0] = mu - bet*y[0]*y[2] ;
   f[1] = bet*y[0]*y[2] - y[1]/lam  ;
   f[2] = y[1]/lam-y[2]/vv         ;
   f[3] = 0.0 ;
   f[4] = 0.0 ;
   f[5] = 0.0 ;
"

bound <- "
  if ( (*i == 1) | (*i == 4)) g[0] = (y[0] - y[3]) ;
  if ( (*i == 2) | (*i == 5)) g[0] = (y[1] - y[4]) ;
  if ( (*i == 3) | (*i == 6)) g[0] = (y[2] - y[5]) ; 
"                             
cBound <- compile.bvp(func = Measel, bound = bound, parms = parms, 
  language = "C", declaration = "double bet;")

## ------------------------
## Solution
## ------------------------

print(system.time(
  sol <- bvptwp(func = cBound, xguess = x, yguess = yguess,
    x=x, leftbc = 3, parms = parms, ncomp = 6, 
    nmax = 100000, atol = 1e-8)                     
))
print(system.time(
  sol2 <- bvptwp(func = cBound, xguess = x, yguess = yguess,
    x=x, leftbc = 3, parms = parms*5, ncomp = 6, 
    nmax = 100000, atol = 1e-8)                     
))
                                                     
plot(sol, sol2)



## ==================  tubular reactor with axial dispersion ===================
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## y'(0) = Pe (y(0)-1),
## y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
##
## The initial condition y'(0) is a function of y(0)
## =============================================================================

parms <- c(Pe = 1, R = 2)

## ------------------------
## Inline C++
## ------------------------

# 2nd order implementation
freac2 <- "
  f[0] = Pe * (y[1] + R*(y[0]*y[0]));
"

jreac2 <- "
  df[0] = Pe * R*2 *y[0];
  df[1] = Pe * y[1];
"

greac <- "
  if (*i == 1) g[0] = (y[1]-Pe*(y[0]-1.0));
  if (*i == 2) g[0] = (y[1]);
"

cfreac2 <- compile.bvp(func = freac2, jacfunc = jreac2, bound = greac,
  parms = parms, language = "C")

## ------------------------
## Solution
## ------------------------

sol <- bvpcol(func = cfreac2, x = seq(0, 1, by = 0.01), order = 2, 
              leftbc = 1, ynames = c("y", "dy"), parms = parms)
sol2 <- bvpcol(func = cfreac2, x = seq(0, 1, by = 0.01), order = 2, nmax = 1e5,
              leftbc = 1, ynames = c("y", "dy"), parms = parms*0.5)
plot(sol, sol2)


## =============================================================================
## "Swirling Flow III", a test problem in Ascher, Mattheij,
## and Russell, Numerical Solution of Boundary Value Problems
## for Ordinary Differential Equations", Classics in Applied
## Mathematics Series, SIAM, Philadelphia, 1995].
## g'' = (g f' - f g')/eps
## f'''' = (-ff'''-gg')/eps
## g(0)=-1,f(0)=0,f'(0)=0, g(1)=1,f(1)=0,f'(1)=0
##
## 1st order system (y1=g, y3=f)
## y1' = y2
## y2' = (y1*y4 -y3*y2)/eps
## y3'=y4
## y4'=y5
## y5'=y6
## y6'=(-y3y6-y1y2)/eps
## y1(0)=-1,y3(0)=0,y4(0)=0, y1(1)=1,y3(1)=0,y4(1)=0
## =============================================================================

x       <- seq(0, 1, 0.01)
yini <- c(-1, NA, 0, 0, NA, NA)
yend <- c(1 , NA, 0, 0, NA, NA)

## ------------------------
## Inline C
## ------------------------

fswirl <- '
  f[0] = y[1];
  f[1] = (y[0]*y[3] - y[2]*y[1])/eps ;
  f[2] = y[3];
 	f[3] = y[4];
 	f[4] = y[5];
  f[5] = (-y[2]*y[5] - y[0]*y[1])/eps ;
'
parms <- c(eps = 0.01)
cswirl <- compile.bvp(func = fswirl, parms = parms, language = "C")


## ------------------------
## Solution
## ------------------------

print(system.time(sol <- bvpcol(x = x, func = cswirl, 
                  yini = yini, yend = yend, #eps = 0.0001, 
                  parms = c(eps = 0.1))))

# For successively smaller values of eps:
print(system.time(sol2 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol[,1], yguess = t(sol[,-1]),
                  yini = yini, yend = yend, eps=0.001, epsini = 0.01,
                  parms=0.001)))

print(system.time(sol3 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol2[,1], yguess = t(sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))

