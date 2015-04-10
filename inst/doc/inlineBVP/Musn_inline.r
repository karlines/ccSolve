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

require(bvpSolve)

## =============================================================================
## R-implementation
## =============================================================================

# Derivatives

musn <- function(t,Y,pars)  {
  with (as.list(Y),   {
   du <- 0.5*u*(w-u)/v
   dv <- -0.5*(w-u)
   dw <- (0.9-1000*(w-y)-0.5*w*(w-u))/z
   dz <- 0.5*(w-u)
   dy <- -100*(y-w)
   return(list(c(du, dv, dw, dz, dy)))
  })
}

x <- seq(0,1,by=0.05)


# 
bound <- function(i, y, pars) {
  with (as.list(y), {
    if (i == 1) return (u-1)
    if (i == 2) return (v-1)
    if (i == 3) return (w-1)
    if (i == 4) return (z+10)
    if (i == 5) return (w-y)
 })
}

xguess <- seq(0, 1, len = 5)
yguess <- matrix(nc = 5,data = rep(c(1,1,1,-10,0.91),5))
rownames(yguess) <- c("u", "v", "w", "z", "y")

print(system.time(
Sol <- bvptwp(yini = NULL, x = x, func = musn, bound = bound,
              xguess = xguess, yguess = yguess, leftbc = 4,
              atol = 1e-10)
))
plot(Sol)



## =============================================================================
## Inline implementation
## =============================================================================

# Derivatives

fmusn = "
   F(1) = 0.5*y(1)*(y(3)-y(1))/y(2)
   F(2) = -0.5*(y(3)-y(1))
   F(3) = (0.9-1000*(y(3)-y(5))-0.5*y(3)*(y(3)-y(1)))/y(4)
   F(4) = 0.5*(y(3)-y(1))
   F(5) = -100*(y(5)-y(3))
" 



# 
fbound <- "
    if (i == 1) g = (y(1)-1)
    if (i == 2) g = (y(2)-1)
    if (i == 3) g = (y(3)-1)
    if (i == 4) g = (y(4)+10)
    if (i == 5) g = (y(3)-y(5))
"
cmusn <- compile.bvp(func = fmusn, bound = fbound)

x <- seq(0, 1, by = 0.05)
xguess <- seq(0, 1, len = 5)
yguess <- matrix(nc = 5,data = rep(c(1,1,1,-10,0.91),5))
rownames(yguess) <- c("u", "v", "w", "z", "y")

print(system.time(
sol <- bvptwp(yini = NULL, x = x, func = cmusn, 
              xguess = xguess, yguess = yguess, leftbc = 4,
              atol = 1e-10)
))
plot(sol)

