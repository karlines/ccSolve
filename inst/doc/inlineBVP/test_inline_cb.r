## ================================================================     comments
## fluid injection problem
## Original definition:
## f'''- R[(f')^2 -f*f''] + A = 0
## h'' + R*f*h' + 1 = 0
## O'' + P*f*O' = 0
## A is unknown, P = 0.7*R
##
## rewritten as:
## df=f1                   #f
## df1=f2                  #f1=f'
## df2=R(f1^2-f*f2)-A      #f2=f''
## dh=h1
## dh1= -Rfh1-1
## dO=O1
## dO1 = O2
## dO2 = -P*f*O1
## dA = 0                  # the constant to be estimated
## Version with parameter common block
## ================================================================ end comments
# load the package with the solvers
require(ccSolve)

pars <- c(R = 200, P = 0.7*200)

times  <- seq(0, 1, by = 0.1)      
yini   <- c(f=0, f1=0, f2=NA, h=0, h1=NA, O=0, O1=NA, A=NA)
yend   <- c(f=1, f1=0, f2=NA, h=0, h1=NA, O=1, O1=NA, A=NA)

# the derivative function; R and P are parameters
ffluid <- "
 f(1) = y(2)
 f(2) = y(3)
 f(3) = R*(y(2)**2-y(1)*y(3))-y(8)
 f(4) = y(5)
 f(5) = -R*y(1)*y(5)-1.d0
 f(6) = y(7)
 f(7) = -P*y(1)*y(7)
 f(8) = 0.d0
"
pars <- c(R = 200, P = 0.7*200)
cfluid <- compile.bvp(func = ffluid, parms = pars)

soltwp <- bvptwp(func = cfluid$func, x = times, parms = pars, initfunc = cfluid$initfunc,
                nmax = 195, cond = TRUE, yini = yini, yend = yend)
pars2 <- c(R = 5000, P = 0.7*5000)

soltwp2 <- bvptwp(func = cfluid$func, x = times, parms = pars2, initfunc = cfluid$initfunc,
                cond = TRUE, yini = yini, yend = yend, nmax = 5000)

plot(soltwp, soltwp2)







## =============================================================================
##   This is the example for MUSN in U. Ascher, R. Mattheij, and R. Russell,
##   Numerical Solution of Boundary Value Problems for Ordinary Differential
##   Equations, SIAM, Philadelphia, PA, 1995.  MUSN is a multiple shooting
##   code for nonlinear BVPs.  The problem is
##
##      u' =  0.5*u*(w - u)/v
##      v' = -0.5*(w - u)
##      w' = (0.9 - 1000*(w - y) - 0.5*w*(w - u))/z
##      z' =  0.5*(w - u)
##      y' = -100*(y - w)
##
##   The interval is [0 1] and the boundary conditions are
##
##      u(0) = v(0) = w(0) = 1,  z(0) = -10,  w(1) = y(1)
##
## note: there are two solutions...
## =============================================================================

## ------------------------
## Problem specifications
## ------------------------

x <- seq(0, 1, by = 0.05)
xguess <- seq(0, 1, len = 5)
yguess <- matrix(nc = 5,data = rep(c(1,1,1,-10,0.91),5))
rownames(yguess) <- c("u", "v", "w", "z", "y")

## ------------------------
## Inline Fortran
## ------------------------

fmusn = "
   F(1) = 0.5*y(1)*(y(3)-y(1))/y(2)
   F(2) = -0.5*(y(3)-y(1))
   F(3) = (0.9-1000*(y(3)-y(5))-0.5*y(3)*(y(3)-y(1)))/y(4)
   F(4) = 0.5*(y(3)-y(1))
   F(5) = -100*(y(5)-y(3))
"

# boundaries 
fbound <- "
    if (i == 1) g = (y(1)-1)
    if (i == 2) g = (y(2)-1)
    if (i == 3) g = (y(3)-1)
    if (i == 4) g = (y(4)+10)
    if (i == 5) g = (y(3)-y(5))
"

musn <- compile.bvp(func = fmusn, bound = fbound)

## ------------------------
## Solution
## ------------------------

print(system.time(
  sol <- bvptwp(yini = NULL, x = x, func = musn$func, bound = musn$bound,
                xguess = xguess, yguess = yguess, leftbc = 4, atol = 1e-10)
))
plot(sol)



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

# Runtime with R: 0.88, 8.11, 13.20, not solved, 26.21
require(bvpSolve)

# Declare global problem dependent parameters.
x       <- seq(0, 1, 0.01)
yini <- c(-1, NA, 0, 0, NA, NA)
yend <- c(1 , NA, 0, 0, NA, NA)

## ============================================================================
## using parameter common block
## ============================================================================

fswirl3 <- "
  f(1) = Y(2)
  f(2) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(3) = Y(4)
 	f(4) = Y(5)
 	f(5) = Y(6)
  f(6) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"

cswirl <- compile.bvp(fswirl3, parms = c(eps = 0.1))

print(system.time(solb <- bvptwp(x = x, func = cswirl, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))

# For successively smaller values of eps:
print(system.time(sol2 <- bvptwp(x = x, func = cswirl, 
                  xguess = solb[,1], yguess = t(solb[,-1]),
                  yini = yini, yend = yend, eps=0.001, epsini = 0.01,
                  parms=0.001)))

print(system.time(sol3 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol2[,1], yguess = t(sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))

print(system.time(sol3b <- bvptwp(x = x, func = cswirl3, 
                  xguess = sol2[,1], yguess = t(sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))

#print(system.time(sol4 <- bvptwp(x = x, func = cswirl, 
#                  xguess = Sol3[,1], yguess = t(Sol3[,-1]),
#                  yini = yini, yend = yend, eps = 5e-5, 
#                  epsini = 1e-4, parms = 5e-5, nmax =1e5)))

# with conditioning 
print(system.time(sol4 <- bvptwp(atol = 1e-5, x = x, func = cswirl, cond = TRUE, 
                  xguess = sol3[,1], yguess = t(sol3[,-1]),
                  yini = yini, yend = yend, eps = 5e-5 , parms=5e-5)))
plot(solb, sol2, sol3, sol4)

#