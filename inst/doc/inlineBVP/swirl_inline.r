## test with eps, as rpar AND common block

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
## Parameter is in rpar
## ============================================================================

fswirl <- "
  eps = rpar(1)
  f(1) = Y(2)
  f(2) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(3) = Y(4)
 	f(4) = Y(5)
 	f(5) = Y(6)
  f(6) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl <- compile.bvp(func = fswirl, declaration = "double precision:: eps")
print(system.time(sol <- bvptwp(x = x, func = cswirl, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))

fswirl2 <- "
  eps = rpar(1)
  f(1) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(2) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl2 <- compile.bvp(func = fswirl2, declaration = "double precision:: eps")
print(system.time(sola <- bvpcol(x = x, func = cswirl2, order = c(2, 4), 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))

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

cswirl3 <- compile.bvp(fswirl3, parms = c(eps = 0.1))

print(system.time(solb <- bvptwp(x = x, func = cswirl3, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))

# For successively smaller values of eps:
print(system.time(sol2 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol[,1], yguess = t(sol[,-1]),
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
plot(sol, sol2, sol3, sol4)

# DOES NOT SOLVE:
#print(system.time(sol4b <- bvpcol(x = x, func = cswirl, 
#                  xguess = sol3[,1], yguess = t(sol3[,-1]),
#                  yini = yini, yend = yend, eps = 5e-5 , parms=5e-5)))

print(system.time(sol4c <- bvpcol(x = x, func = cswirl2, order = c(2, 4), 
                  xguess = sol3[,1], yguess = t(sol3[,-1]),
                  yini = yini, yend = yend, eps = 5e-5 , epsini = 1e-4, parms=5e-5)))
