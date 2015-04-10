
## =============================================================================
## Standard linear problem with boundary layer at the origin
##
##
## d2y/dt^2=-3py/(p+t^2)^2
## y(t= -0.1)=-0.1/sqrt(p+0.01)
## y(t=  0.1)= 0.1/sqrt(p+0.01)
## where p = 1e-5
##
## analytical solution y(t) = t/sqrt(p + t^2).
##
## The problem is rewritten as a system of 2 ODEs:
## dy=y2
## dy2=-3p*y/(p+t^2)^2
##
## Solved using shooting and bvptwp
## =============================================================================

require(bvpSolve)

#--------------------------------
# Functions
#--------------------------------
fbnd <- "
  F(1) = y(2)
  F(2) = - 3.d0 * rpar(1) * y(1)/(rpar(1) + x*x)**2.
"
fjac <- "
  df(1, 2) = 1.d0
  df(2, 1) = -3.d0*rpar(1)/(rpar(1) +x*x)**2.
"
fbound <- "
  if (i == 1) then
     g = (y(1) + 0.1/sqrt(rpar(1) + 0.01d0))
  else  
     g = (y(1) - 0.1/sqrt(rpar(1) + 0.01d0))
  end if
"
fjacbound <- "
  dg(1) = 1.d0
"
# parameter value


p    <-1e-5
x <- seq(-1, 1, length.out = 100)
print(system.time(Sol <- bvptwp(x = x, leftbc = 1, func = cbnd, 
        ncomp = 2, verbose = FALSE, rpar = p)))

# all compiled in one go
cbnd <- compile.bvp (func = fbnd, jacfunc = fjac, bound = fbound, jacbound = fjacbound)
print(system.time(sol <- bvptwp(x = x, leftbc = 1, func = cbnd, 
        ncomp = 2, verbose = FALSE, rpar = p)))
plot(sol, which = 1)




for (pp in 0:9){
  Soln <- bvptwp(x = x, leftbc = 1, func = cbnd, nmax = 1e5,
        ncomp = 2, rpar = 10^(-pp))
  lines(Soln[,1], Soln[,2])      
}  
        
