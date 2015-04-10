## ==================  tubular reactor with axial dispersion ===================
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## y'(0) = Pe (y(0)-1),
## y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
##
## The initial condition y'(0) is a function of y(0)
## Can be solved with bvpshoot
## =============================================================================
require(bvpSolve)

# Not so simple to solve with bvpshoot
Reactor <- function(x, y, parms)  {
  list(c(y[2], Pe*(y[2]+R*(y[1]^n))))
}


bound <- function(i,y,p) {
  if (i == 1) return(y[2]-Pe*(y[1]-1))
  if (i == 2) return(y[2])
}

Pe <- 1
R  <- 2
n  <- 2
Sol<-bvptwp(func = Reactor, x = seq(0, 1, by = 0.01),
            leftbc = 1, ynames = c("y", "dy"), bound = bound)
plot(Sol)

freac <- "
  Pe = rpar(1)
  R = rpar(2)

  F(1) = y(2)
  F(2) = Pe * (y(2) + R*(Y(1)**2.))
"

greac <- "
  Pe = rpar(1)
  if (i == 1) g = (y(2)-Pe*(y(1)-1.d0))
  if (i == 2) g = (y(2))
"
cfreac <- compile.bvp(func = freac, bound = greac, declaration = " DOUBLE PRECISION :: Pe, R")

sol <- bvptwp(func = cfreac, x = seq(0, 1, by = 0.01),
              leftbc = 1, ynames = c("y", "dy"), rpar = c(Pe, R))
plot(sol)
