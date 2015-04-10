## =============================================================================
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

mathieu <- function(t,y,parms) {
 list(c(y[2],
        -(y[3]-10*cos(2*t))*y[1],
        0 ))
}
jac <- function(x, y ,p) {
  df <- matrix(nr = 3, nc = 3, 0)
  df[1,2] <- 1
  df[2,1] <- -(y[3]-10*cos(2*x))
  df[2,3] <- -y[1]
  df
}

bound <- function(i,y,parms){
  if (i ==1) return(y[1]-1)
  if (i ==2) return(y[2])
  if (i ==3) return(y[2])
}

jacbound <- function(i,y,parms){
  if (i ==1) return(c(1,0,0))
  else return(c(0,1,0))
}

## =============================================================================
## Problem specifications
## =============================================================================

x      <- seq(0, pi, by = 0.01)
init   <- c(y = 1,dy = 0, lambda = NA)
xguess <-  c(0, 1, 2*pi)
yguess <- matrix(nr = 3, rep(1, 9))
rownames(yguess) <- c("y", "dy", "lambda")


## =============================================================================
## Inline fortran derivative
## =============================================================================
Mathieu <- "
       F(1) = y(2)
       F(2) = -(y(3)-10*cos(2*x))*y(1)
       F(3) = 0
"

## =============================================================================
## Inline fortran jacobian
## =============================================================================
Jac <- "
  df(1,2) = 1
  df(2,1) = -(y(3)-10*cos(2*x))
  df(2,3) = -y(1)
"

## =============================================================================
## Inline fortran boundary
## =============================================================================
Bound  <- "
  if (i ==1) g = (y(1)-1.d0)
  if (i ==2) g = (y(2))
  if (i ==3) g = (y(2))
"

## =============================================================================
## Inline fortran boundary jacobian
## =============================================================================

Jacbound <- "
  if (i ==1) THEN
    DG(1) = 1
  else
    DG(2) = 1  
  endif
"

cmathieu <- compile.bvp(func = Mathieu, jacfunc = Jac, bound = Bound, jacbound = Jacbound)
sol5  <- bvptwp(func = cmathieu, leftbc = 2, x = x,
        xguess = xguess, yguess = yguess)

plot(sol5)