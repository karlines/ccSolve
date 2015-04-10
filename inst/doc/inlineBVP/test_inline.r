## uses rpar....

## ================================================================     comments
## fluid injection problem
## Original definition:
## f'''- R[(f')^2 -f*f''] + A = 0
## h'' + R*f*h' + 1 = 0
## O'' + P*f*O' = 0
## A is unknown, P = 0.7*R
##
## rewritten as:
## df=f1                   #f
## df1=f2                  #f1=f'
## df2=R(f1^2-f*f2)-A      #f2=f''
## dh=h1
## dh1= -Rfh1-1
## dO=O1
## dO1 = O2
## dO2 = -P*f*O1
## dA = 0                  # the constant to be estimated
## ================================================================ end comments


# load the package with the solver
require(bvpSolve)

times  <- seq(0, 1, by = 0.1)      
yini   <- c(f=0, f1=0, f2=NA, h=0, h1=NA, O=0, O1=NA, A=NA)
yend   <- c(f=1, f1=0, f2=NA, h=0, h1=NA, O=1, O1=NA, A=NA)

# the derivative function
ffluid <- "
 R    = rpar(1)
 P    = 0.7d0*R
 f(1) = y(2)
 f(2) = y(3)
 f(3) = R*(y(2)**2-y(1)*y(3))-y(8)
 f(4) = y(5)
 f(5) = -R*y(1)*y(5)-1.d0
 f(6) = y(7)
 f(7) = -P*y(1)*y(7)
 f(8) = 0.d0
"
cfluid <- compile.bvp(ffluid, declaration = "Double precision :: R, P")

ffluid2 <- "
       R    = rpar(1)
       P    = 0.7d0*R
       f(1) = y(2)
       f(2) = y(3)
       f(3) = R*(y(2)**2-y(1)*y(3))-y(8)
       f(4) = y(5)
       f(5) = -R*y(1)*y(5)-1.d0
       f(6) = y(7)
       f(7) = -P*y(1)*y(7)
       f(8) = 0.d0
"
cfluid2 <- compile.bvp(ffluid2, declaration = "       Double precision :: R, P", language = "Fortran")


ffluid3 <- "
       R    = rpar[0];
       P    = 0.7*R;
       f[0] = y[1];
       f[1] = y[2];
       f[2] = R*(y[1]*y[1]-y[0]*y[2])-y[7];
       f[3] = y[4];
       f[4] = -R*y[0]*y[4]-1.0;
       f[5] = y[6];
       f[6] = -P*y[0]*y[6];
       f[7] = 0.0;
"
cfluid3 <- compile.bvp(ffluid3, declaration = " double R, P;", language = "C")

soltwp <- bvptwp(func = cfluid, x = times, parms = NULL, rpar = 200,
                nmax = 195, cond = TRUE, yini = yini, yend = yend)
diagnostics(soltwp)

# Solving the model, using twpbvpC + very large value of R + conditioning
system.time(
soltwp2 <- bvptwp(func = cfluid, x = times, parms = NULL, rpar = 5000,
                nmax = 1100, cond = TRUE, #verbose = TRUE,    #
                yini = yini, yend = yend, xguess = soltwp[,1],
                yguess = t(soltwp[,-1]))
)

system.time(
soltwp2 <- bvptwp(func = cfluid2, x = times, parms = NULL, rpar = 5000,
                nmax = 1100, cond = TRUE, #verbose = TRUE,    #
                yini = yini, yend = yend, xguess = soltwp[,1],
                yguess = t(soltwp[,-1]))
)

system.time(
soltwp2 <- bvptwp(func = cfluid3, x = times, parms = NULL, rpar = 5000,
                nmax = 1100, cond = TRUE, #verbose = TRUE,    #
                yini = yini, yend = yend, xguess = soltwp[,1],
                yguess = t(soltwp[,-1]))
)
diagnostics(soltwp2)

# plot the results
plot(soltwp, soltwp2, type = "l", lwd = 2)

