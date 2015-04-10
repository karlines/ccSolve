## =============================================================================
##   This is a nerve impulse model considered in Example 7.1 of
##   R. Seydel, From Equilibrium to Chaos, Elsevier, 1988.  The
##   differential equations
##
##      y1' = 3*t*(y1 + y2 - 1/3*y1^3 + lambda)
##      y2' = -t*(y1 - 0.7 + 0.8*y2)/3
##
##   are to be solved subject to periodic boundary conditions
##
##      y1(0) = y1(1)
##      y2(0) = y2(1)
##
##
##   The parameter lambda has the value -1.3.  The
##   period T is unknown, i.e the length of the interval is not known.
##
##   An extra equation to estimate T is: -T(y1(0)-0.7+0.8*y2(0))/3=1
##
## =============================================================================


require(bvpSolve)

## =============================================================================
## use augmented equation with bvptwp... 
## solvable only when good initial conditions are used...
## =============================================================================

nerve3 <- function (t, y ,p)
  list(c( 3*y[3]*(y[1] + y[2] - 1/3*(y[1]^3) - 1.3),
        (-1/3)*y[3]*(y[1] - 0.7 + 0.8*y[2]) ,
        0,
        0,
        0)
  )

# derivate function
dnerve3 <- function (t, y, p)  {
  df <- matrix(nr = 5, nc = 5, data = 0)
  
  df[1,1] <- 3.0*y[3] -3*y[3]*y[1]^2.
  df[1,2] <- 3.0*y[3]
	df[1,3] <- 3.0*(y[1] + y[2] - 1/3*(y[1]^3) - 1.3)
  df[2,1] <- (-1.0/3)*y[3]
	df[2,2] <- 0.80* (-1.0/3)*y[3]
  df[2,3] <- (-1.0/3)*(y[1] - 0.7 + 0.8*y[2]) 
  return(df)
}
# boundary function
bound <- function(i, y, p) {
 if (i == 1) return(y[3]*(-1/3)*(y[1] - 0.7 + 0.8*y[2]) - 1)
 if (i == 2) return(y[1]-y[4])
 if (i == 3) return(y[2]-y[5])
 if (i == 4) return(y[1]-y[4])
 if (i == 5) return(y[2]-y[5])
}

# needs good initial guesses of the solution.
xguess <- seq(0, 1, by = 0.1)
yguess <- matrix(nr = 5, nc = length(xguess), 5.)
yguess[1,] <- sin(2*pi*xguess)
yguess[2,] <- cos(2*pi*xguess)

print (system.time(
Sol  <- bvpcol(func = nerve3, bound = bound, jacfunc = dnerve3, 
        x = seq(0, 1, by = 0.01), 
        ynames = c("y","dy","T","yi","yj"),
        leftbc = 3, xguess = xguess, yguess = yguess)
))
plot(Sol)

print (system.time(
Sol2  <- bvptwp(func = nerve3, bound = bound, jacfunc = dnerve3, 
        x = seq(0, 1, by = 0.01), 
        ynames = c("y","dy","T","yi","yj"),
        leftbc = 3, xguess = xguess, yguess = yguess)
))        
plot(Sol2)

### inline fortran:
# derivative function
fnerve = "
 F(1) = 3.d0*y(3)*(y(1) + y(2) - 1.d0/3*(y(1)**3) - 1.3)
 F(2) = (-1.d0/3)*y(3)*(y(1) - 0.7d0 + 0.8*y(2)) 
 F(3) = 0.d0
 F(4) = 0.d0
 F(5) = 0.d0"

# jacobian
dfnerve  = "
  
  df(1,1) = 3.0d0*y(3) -3.d0*y(3)*y(1)**2.
  df(1,2) = 3.0d0*y(3)
	df(1,3) = 3.0d0*(y(1) + y(2) - 1.d0/3*(y(1)**3) - 1.3d0)
  df(2,1) = (-1.0d0/3)*y(3)
	df(2,2) = 0.80d0* (-1.0d0/3)*y(3)
  df(2,3) = (-1.0d0/3)*(y(1) - 0.7d0 + 0.8d0*y(2)) 
"

# boundary function
boundnerve <- "
 if (i == 1) g = (y(3)*(-1.d0/3)*(y(1) - 0.7d0 + 0.8d0*y(2)) - 1.d0)
 if (i == 2) g = (y(1)-y(4))
 if (i == 3) g = (y(2)-y(5))
 if (i == 4) g = (y(1)-y(4))
 if (i == 5) g = (y(2)-y(5))
"

cfnerve <- compile.bvp(func = fnerve, jacfunc = dfnerve, bound = boundnerve)

print (system.time(
Sol  <- bvpcol(func = cfnerve, 
        x = seq(0, 1, by = 0.01), 
        ynames = c("y", "dy", "T", "yi", "yj"),
        leftbc = 3, xguess = xguess, yguess = yguess)
))
plot(Sol)
