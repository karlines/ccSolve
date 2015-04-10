## test with DAE
## =============================================================================
## Example 5 - a bvp DAE
## =============================================================================
require(ccSolve)

# The R-functions
bvpdae <- function(t, x, ks, ...) {
   p1  <- p2 <- sin(t)
   dp1 <- dp2 <- cos(t)
   
   dx1 <- (ks + x[2] - p2)*x[4] + dp1
   dx2 <- dp2
   dx3 <- x[4]
   res <- (x[1] - p1)*(x[4] - exp(t))
 
   list(c(dx1, dx2, dx3, res), res = res)
 }


boundfun <- function(i,  x, par, ...) {
   if (i == 1) return(x[1] - sin(0))
   if (i == 2) return(x[3] - 1)
   if (i == 3) return(x[2] - sin(1))
   if (i == 4) return((x[1] - sin(1))*(x[4] - exp(1)))  # Not used here..
 }

x <- seq(0, 1, by = 0.01)

mass <- diag(nrow = 4)  ; mass[4, 4] <- 0

# solved using boundfun
out <- bvpcol (func = bvpdae, bound = boundfun, x = x, 
                parms = 1e-4, ncomp = 4, leftbc = 2,
                dae = list(index = 2,  nalg = 1)) 

# solved using yini, yend
out1 <- bvpcol (func = bvpdae, x = x, parms = 1e-4, 
                 yini = c(sin(0), NA, 1, NA), 
                 yend = c(NA, sin(1), NA, NA),
                 dae = list(index = 2,  nalg = 1)) 

# the analytic solution
ana <- cbind(x, "1" = sin(x), "2" = sin(x), "3" = 1, "4" = 0, res = 0)

plot(out, out1, obs = ana)

# FORTRAN 95 IMPLEMENTATION

fbvpdae <- "
   f(1) = (ks + y(2) - sin(x))*y(4) + cos(x)
   f(2) = cos(x)
   f(3) = y(4)
   f(4) = (y(1) - sin(x))*(y(4) - exp(x))
"

fboundfun <- "
   if (i == 1) then 
     g = y(1) - sin(0.d0)
   else if (i == 2) then 
     g = y(3) - 1.d0
   else if (i == 3) then
     g = y(2) - sin(1.d0)
   else if (i == 4) then
     g = (y(1) - sin(1.d0))*(y(4) - exp(1.d0))  
   endif
"
parms <- c(ks = 1e-4)
cdae <- compile.bvp(func = fbvpdae, bound = fboundfun, parms = parms)
cdae2 <- compile.bvp(func = fbvpdae, parms = parms)
outc <- bvpcol (func = cdae, x = x, parms = 1e-4, 
                ncomp = 4, leftbc = 2,
                 dae = list(index = 2,  nalg = 1)) 
# solved using yini, yend
out1 <- bvpcol (func = cdae2, x = x, parms = 1e-4, 
                 yini = c(sin(0), NA, 1, NA), 
                 yend = c(NA, sin(1), NA, NA),
                 dae = list(index = 2,  nalg = 1)) 
