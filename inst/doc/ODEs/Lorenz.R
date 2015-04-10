## =========================================================================
## The Lorenz equation 
## =========================================================================
# --------------------------------------------------------------------------
# implemented in R
# --------------------------------------------------------------------------

require(deSolve)

chaos <- function(t, state, parameters) {
  with(as.list(c(state)), {

    dxx     <- -8/3 * xx + yy * zz
    dyy     <- -10 * (yy - zz)
    dzz     <- -xx * yy + 28 * yy - zz

    list(c(dxx, dyy, dzz))
  })
}
state <- c(xx = 1, yy = 1, zz = 1)
times <- seq(0, 100, 0.01)
print(system.time(
  out   <- vode(state, times, chaos, 0)
))
# --------------------------------------------------------------------------
# In  Fortran 95, passing state variable names:
# --------------------------------------------------------------------------
chaos.f95 <- " 
    dxx     = -8.d0/3 * xx + yy * zz
    dyy     = -10.d0 * (yy - zz)
    dzz     = -xx * yy + 28d0 * yy - zz
"
cChaos <- compile.ode(chaos.f95, y = state) 
print(system.time(
  cout   <- vode(state, times, func = cChaos, parms = 0)
))

# --------------------------------------------------------------------------
# in C, passing parameter values, but not the state variable names:
# --------------------------------------------------------------------------
parms <- c(a = -8.0/3, b = -10.0, c = 28.0)

chaos.C <- " 
   f[0] = a * y[0] + y[1]*y[2];
   f[1] = b*(y[1]-y[2]);
   f[2] = -y[0]*y[1] + c * y[1] - y[2];
"
parms <- c(a = -8.0/3, b = -10.0, c = 28.0)
cChaos2 <- compile.ode(chaos.C, language = "C", parms = parms) 
print(system.time(
  cout2  <- vode(state, times, func = cChaos2, parms = parms)
))

# output is the same
plot(out, cout, cout2)    
plot(out[,"xx"], out[,"yy"], type = "l", main = "Lorenz butterfly",
  xlab = "x", ylab = "y")
