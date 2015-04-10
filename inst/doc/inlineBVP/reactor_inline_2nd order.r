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
require(bvpSolve)

# 2nd order equation solved  

Reactor2 <- function(x, y, parms)  {
  list(Pe*(y[2]+R*(y[1]^n)))
}
Pe <- 1
R  <- 2
n  <- 2

# boundary function...

bound <- function(i,y,p) {
  if (i == 1) return(y[2]-Pe*(y[1]-1))
  if (i == 2) return(y[2])
}

Sol2 <- bvpcol(func = Reactor2, x = seq(0, 1, by = 0.01), order = 2,
             leftbc = 1, ynames = c("y", "dy"), bound = bound)

freac2 <- "
  F(1) = Pe * (y(2) + R*(Y(1)**2.))
"
greac <- "
  if (i == 1) g = (y(2)-Pe*(y(1)-1.d0))
  if (i == 2) g = (y(2))
"

parms <- c(Pe = 1, R = 2)

cfreac2 <- compile.bvp(func = freac2, bound = greac, parms = parms)

## ------------------------
## Solution
## ------------------------

sol <- bvpcol(func = cfreac2, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), parms = parms)

plot(sol)

