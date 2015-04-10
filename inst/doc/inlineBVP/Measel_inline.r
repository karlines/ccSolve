## =============================================================================
## PROBLEM measels
## Models the spread of measels in three equations
## U. M. Ascher, R. M. R. Mattheij, and R. D. Russell. Numerical Solution of
## Boundary Value Problems for Ordinary Differential Equations. Prentice{Hall,
## Englewood Cliffs, NJ, USA, 1988.
## =============================================================================

require(bvpSolve)

## =============================================================================
## R implementation
## =============================================================================

measel <- function(t, y, pars, vv)  {
  bet <- 1575*(1+cos(2*pi*t))
  dy1 <- mu-bet*y[1]*y[3]
  dy2 <- bet*y[1]*y[3]-y[2]/lam
  dy3 <-y[2]/lam-y[3]/vv
  dy4 <- 0
  dy5 <-0
  dy6 <-0
  
  list(c(dy1, dy2, dy3, dy4, dy5, dy6))
}

dmeasel <- function(t, y, pars, vv) {
  df <- matrix (data = 0, nrow = 6, ncol = 6)
  bet <- 1575*(1+cos(2*pi*t))
  df[1,1] <-  -bet*y[3]
  df[1,3] <-  -bet*y[1]

  df[2,1] <-  bet*y[3]
  df[2,2] <-  -1/lam
  df[2,3] <-  bet*y[1]

  df[3,2] <- 1/lam 
  df[3,3] <- -1/vv
  
  return(df)
}

bound <- function(i, y, pars,vv) {
  if ( i == 1 | i == 4) return(y[1]-y[4])
  if ( i == 2 | i == 5) return(y[2]-y[5])
  if ( i == 3 | i == 6) return(y[3]-y[6])  
}

dbound <- function(i, y, pars,vv) {
  if ( i == 1 | i == 4) return(c(1,0,0,-1,0,0))
  if ( i == 2 | i == 5) return(c(0,1,0,0,-1,0))
  if ( i == 3 | i == 6) return(c(0,0,1,0,0,-1))
}

mu  <- 0.02
lam <- 0.0279
v   <- 0.1

yguess <- matrix(ncol = length(x), nrow = 6, data = 1)
rownames(yguess) <- paste("y", 1:6, sep = "")


print(system.time(
  solR <- bvptwp(fun = measel, bound = bound, 
    xguess = x, yguess = yguess,
    x=x, leftbc = 3, vv = v, ncomp = 6, 
    nmax = 100000, atol = 1e-4)
))

plot(sol)

## =============================================================================
## Inline implementation
## =============================================================================

measel <- " 
   bet = 1575d0*(1.+cos(2*pi*x))
   F(1) = mu - bet*y(1)*y(3)
   F(2) = bet*y(1)*y(3) - y(2)/lam
   F(3) = y(2)/lam-y(3)/vv
   F(4) = 0.d0
   F(5) = 0.d0
   F(6) = 0.d0
"

dmeasel <- "
  bet = 1575d0*(1+cos(2*pi*x))
  df(1,1) =  -bet*y(3)
  df(1,3) =  -bet*y(1)

  df(2,1) =  bet*y(3)
  df(2,2) =  -1.d0/lam
  df(2,3) =  bet*y(1)

  df(3,2) = 1.d0/lam 
  df(3,3) = -1.d0/vv
"

bound <- "
  if ( i == 1 .OR. i == 4) g = (y(1) - y(4))
  if ( i == 2 .OR. i == 5) g = (y(2) - y(5))
  if ( i == 3 .OR. i == 6) g = (y(3) - y(6))  
"

dbound <- "
  if ( i == 1 .OR. i == 4) THEN
    dg(1) = 1.
    dg(4) = -1.
  else if ( i == 2 .OR. i == 5) then
    dg(2) = 1.
    dg(5) = -1.
  else
    dg(3) = 1.
    dg(6) = -1.
  end if
"  

parms <- c(vv = 0.1, mu = 0.02, lam = 0.0279)
cMeasel <- compile.bvp(func = measel, jacfunc = dmeasel, 
  bound = bound, jacbound = dbound, parms  = parms, 
  declaration = "double precision, parameter :: pi = 3.141592653589793116d0\n double precision :: bet")

x <- seq (0, 1, by = 0.01)
yguess <- matrix(ncol = length(x), nrow = 6, data = 1)
rownames(yguess) <- paste("y", 1:6, sep = "")

print(system.time(
  sol1 <- bvptwp(func = cMeasel, 
    xguess = x, yguess = yguess,
    x = x, leftbc = 3, parms = parms, ncomp = 6, 
    nmax = 100000, atol = 1e-8)
))


print(system.time(
  sol2 <- bvptwp(func = cMeasel, 
    xguess = x, yguess = yguess,
    x=x, leftbc = 3, parms = parms * c(1, 2, 2) , ncomp = 6, 
    nmax = 100000, atol = 1e-8)
))

plot(sol1, sol2)
