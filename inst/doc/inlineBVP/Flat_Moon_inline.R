## =============================================================================
##
## flat moon problem
## finding the optimal ascent trajectory for launch from a flat moon to 
## a 100 nautical mile circular orbat
##
## assumptions: constant thrust, constant mass, no drag
##  uniform flat-moon gravity
##
## implementation: Karline Soetaert
## =============================================================================

require(bvpSolve)

x    <- seq(0, 1, length.out = 100)   # non-dimensional time vector

h    <- 185.2e3 # metres, final altitude at 100nmi
VC   <- 1627.0  # m/s circular spped at 100nmi
g    <- 1.62    # m/s2 gravitational acceleration
t2w  <- 3       # thrust to weight ratio for ascent vehicle, in lunar g's
A    <- g* 3    # accelleration

# boundary conditions
yini <- c(xx  = 0, yy = 0, Vx = 0,  Vy = 0, l2 = NA, l4 = NA, tau =NA)
yend <- c(xx = NA, yy = h, Vx = VC, Vy = 0, l2 = NA, l4 = NA, tau = NA)


flatmoon <- function(t, y, p) {
   res <- c(y[3],                                        #dx /dtau
            y[4],                                        #dy /dtau
            A * (1/sqrt(1+y[6]^2)),                      #dVx/dtau
            A * (y[6]/sqrt(1+y[6]^2)) - g,               #dVy/dtau
            0.,                                          #dl2/dtau
            -y[5],                                       #dl4/dtau
            0.                                           #dt /dtau
           ) 
    list (res * y[7])
}

print (system.time(
sol <-  bvptwp(func = flatmoon, x = x, yini = yini, yend = yend, xguess = x,
  yguess = matrix (nrow = 7, ncol = length(x), data = c(0, 0, 0, 0, 0, 0, 700)))
))

sol.unscaled <- sol 
sol.unscaled[,1] <- sol[,1]*sol[,ncol(sol)]

plot(sol.unscaled)

fflat <- "

      dxx = Vx * tau                           
      dyy = Vy * tau                             
      dvx = A*(1.d0/sqrt(1+y(6)**2))* tau         
      dvy = (A*(y(6)/sqrt(1+y(6)**2)) - g)* tau   
      dl2 = 0.d0                                  
      dl4 = -y(5) * tau                           
      dtau = 0.d0                                  
"
x    <- seq(0, 1, length.out = 100)   # non-dimensional time vector

h    <- 185.2e3 # metres, final altitude at 100nmi
VC   <- 1627.0  # m/s circular spped at 100nmi
g    <- 1.62    # m/s2 gravitational acceleration
t2w  <- 3       # thrust to weight ratio for ascent vehicle, in lunar g's
A    <- g* 3    # accelleration

# boundary conditions
yini <- c(xx  = 0, yy = 0, Vx = 0,  Vy = 0, l2 = NA, l4 = NA, tau =NA)
yend <- c(xx = NA, yy = h, Vx = VC, Vy = 0, l2 = NA, l4 = NA, tau = NA)

parms <- c(g = g, A = 3*g)

cflat <- compile.bvp(fflat, yini = yini, parms = parms)

# needs good initial conditions
print (system.time(
sol2 <- bvptwp(func = cflat, x = x, yini = yini, yend = yend, xguess = x, 
  yguess = matrix (nrow = 7, ncol = length(x), data = c(0, 0, 0, 0, 0, 0, 700)),
  parms = parms,  nmax = 1e5)
))

print (system.time(
sol2 <- bvptwp(func = cflat, x = x, yini = yini, yend = yend, xguess = x, 
  yguess = matrix (nrow = 7, ncol = length(x), data = c(0, 0, 0, 0, 0, 0, 700)),
  parms = parms*2,  nmax = 1e5)
))


plot(sol, sol2)
  