## ==========================================================================
## A Simple BVP Example implemented in Fortran and in C
## ==========================================================================

#  problem 7 from the test problems available
#  from \url{http://www.ma.ic.ac.uk/~jcash/BVP_software/readme.php} ):
#
#  \begin{eqnarray*}
#    \xi y'' + x y' - y &=&  -(1 + \xi \pi ^2) \cos(\pi x) -\pi x \sin(\pi x)\\
#    y(-1) &=& -1 \\
#    y(1) &=& 1
#  \end{eqnarray*}

fun.f95 <- "
 f(1) = 1/ks * (-x * y(2) + y(1)-(1 + ks*3.14159**2)*cos(3.14159*x)-
    3.14159*x*sin(3.14159*x))
"

fun.C <- "
 f[0] = 1/ks * (-1 * *x *y[1]+y[0]-(1+ks*3.14159*3.14159)*cos(3.14159* *x)-
    3.14159* *x*sin(3.14159* *x));
"

ks <- 0.1
x  <- seq(-1, 1, by = 0.01)
cfun <- compile.bvp(fun.C, parms = c(ks = 0.1), language = "C")
sol1  <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = x, 
  parms = 0.1, func = cfun, order = 2)
sol2  <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = x, 
  parms = 0.01, func = cfun, order = 2)
sol3  <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = x, 
  parms = 0.001, func = cfun, order = 2)
sol4  <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = x, 
  parms = 0.0001, func = cfun, order = 2)
plot(sol1, sol2, sol3, sol4, type = "l", main = "test problem 7, ksi=0.1..0.0001",
     lwd = 2, lty = 1)
