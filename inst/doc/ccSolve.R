### R code from vignette source 'ccSolve.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: preliminaries
###################################################
library("ccSolve")
options(prompt = " ")
options(continue = "  ")
options(width=70)


###################################################
### code chunk number 2: ccSolve.Rnw:158-160
###################################################
args(compile.ode)
args(compile.bvp)


###################################################
### code chunk number 3: ccSolve.Rnw:164-165
###################################################
args(compile.optim)


###################################################
### code chunk number 4: ccSolve.Rnw:217-218
###################################################
args(compile.uniroot)


###################################################
### code chunk number 5: ccSolve.Rnw:228-233
###################################################
rootR <- function (x) 1/cos(1+x^2)

print(system.time(
 for (i in 1:100) AA <- uniroot.all(rootR, c(-10, 10))
))


###################################################
### code chunk number 6: ccSolve.Rnw:243-253
###################################################
croot.f95 <- compile.uniroot("f = 1.d0 / cos(1.d0+x*x)")
croot.C   <- compile.uniroot("f[0] = 1.0/cos(1.0+x[0]*x[0]);", language = "C")

print(system.time(
 for (i in 1:100) A2 <- ccuniroot.all(croot.f95, c(-10, 10))
))
print(system.time(
 for (i in 1:100) A3 <- ccuniroot.all(croot.C, c(-10, 10))
))
max(abs(AA-A2))


###################################################
### code chunk number 7: ccSolve.Rnw:260-262
###################################################
code(croot.f95)
code(croot.C)


###################################################
### code chunk number 8: ccSolve.Rnw:265-266
###################################################
ccfunc(croot.f95, x = 1)


###################################################
### code chunk number 9: ccSolve.Rnw:276-277
###################################################
args(compile.multiroot)


###################################################
### code chunk number 10: ccSolve.Rnw:285-286
###################################################
args(multiroot)


###################################################
### code chunk number 11: ccSolve.Rnw:293-299
###################################################
fun.R <- function(x){
  c(x[1] - 4*x[1]^2 - x[1]*x[2],
    2*x[2] - x[2]^2 + 3*x[1]*x[2] )
}
sol <- multiroot(f = fun.R, start = c(1, 1))
sol


###################################################
### code chunk number 12: ccSolve.Rnw:305-311
###################################################
fun.f95 <- "
 f(1) = x(1) - 4.d0*x(1)**2. - x(1) *x(2)
 f(2) = 2.d0*x(2) - x(2)**2 + 3.d0*x(1)*x(2)
"
cfun.f95 <- compile.multiroot(fun.f95)
multiroot(f = cfun.f95, start = c(1, 1))


###################################################
### code chunk number 13: ccSolve.Rnw:320-328
###################################################
jacfun.f95 <- "
  df (1, 1) = 1.d0 - 8.d0*x(1) - x(2)
  df (1, 2) = -x(1)
  df (2, 1) = 3.d0*x(2)
  df (2, 2) = 2.d0 - 2.d0*x(2) + 3.d0*x(1)
"
cfunjac.f95 <- compile.multiroot(func = fun.f95, jacfunc = jacfun.f95)
multiroot(f = cfunjac.f95, start = c(1, 1), jactype = "fullusr")


###################################################
### code chunk number 14: ccSolve.Rnw:331-332
###################################################
code(cfunjac.f95)


###################################################
### code chunk number 15: ccSolve.Rnw:350-359
###################################################
sixeq.f95 <- "
 f(1) = x(1) + x(2)/x(6) + x(3) + a*x(4) - b
 f(2) = a*x(3) + c*x(4) + x(5) + x(6) - d
 f(3) = x(1) + b*x(2) + exp(x(4)) + x(5) + x(6) + e
 f(4) = a*x(3) + x(3)*x(5) - x(2)*x(3) - ff*(x(5)**2)
 f(5) = g*(x(3)**2) - x(4)*x(6)
 f(6) = h*x(1)*x(6) - x(2)*x(5)
"
parms <- c(a = 2, b = 2, c = 3, d = 4, e = 8, ff = 0.1, g = 8, h = 50)


###################################################
### code chunk number 16: ccSolve.Rnw:362-364
###################################################
csixeq <- compile.multiroot(sixeq.f95, parms = parms)
multiroot(f = csixeq, start = rep(1, 6), parms = parms)


###################################################
### code chunk number 17: ccSolve.Rnw:375-376
###################################################
code(csixeq)


###################################################
### code chunk number 18: ccSolve.Rnw:382-391
###################################################
rosenbrock.R <- function(x) {
  f[i.uneven] <- 1 - x[i.uneven]
  f[i.even]   <- 10 *(x[i.even] - x[i.uneven]^2)
  f
}
n <- 100000
i.uneven <- seq(1, n-1, by = 2)
i.even <- i.uneven + 1
f <- vector(length = n)


###################################################
### code chunk number 19: ccSolve.Rnw:404-407
###################################################
print(system.time(
AR <- multiroot.1D(f = rosenbrock.R, start = runif(n), nspec = 1))
)


###################################################
### code chunk number 20: ccSolve.Rnw:411-422
###################################################
rosenbrock.f95 <- "
 integer i
 do i  = 1, n-1, 2
  f(i) = 1 - x(i)
 enddo
 do i  = 2, n, 2
  f(i) = 10 *(x(i) - x(i-1)**2)
 enddo
"

cRosenbrock.f95 <- compile.multiroot(rosenbrock.f95)


###################################################
### code chunk number 21: ccSolve.Rnw:425-433
###################################################
rosenbrock.C <- "
 int i;
 for(i = 0; i < *n-1; i = i+2)
  f[i] = 1 - x[i];
 for(i = 1; i < *n; i = i+2)
  f[i] = 10 *(x[i] - x[i-1]*x[i-1]);
"
cRosenbrock.C <- compile.multiroot(rosenbrock.C, language = "C")


###################################################
### code chunk number 22: ccSolve.Rnw:438-441
###################################################
print(system.time(
  A <- multiroot.1D(f = cRosenbrock.f95, start = runif(100000),  nspec = 1))
)


###################################################
### code chunk number 23: ccSolve.Rnw:444-445
###################################################
range(A$root)


###################################################
### code chunk number 24: ccSolve.Rnw:471-481
###################################################
DNase1 <- subset(DNase, Run == 1)
head (DNase1)  

print(system.time(
for (i in 1:100)
fm <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(xmid = 0, scal = 1))
))
summary(fm)


###################################################
### code chunk number 25: ccSolve.Rnw:498-503
###################################################
fm.f95 = "f = density - 1.d0/(1.d0 + dexp((xmid - dlog(conc))/scal))"

cfm.f95 <- compile.nls(func = fm.f95, par = c(xmid = 0, scal = 1), 
  data = DNase1[,-1])  
code(cfm.f95)


###################################################
### code chunk number 26: ccSolve.Rnw:506-512
###################################################
cfm.C <- compile.nls(func = '
 int i; 
 for (i = 0; i < *ndat; i++)
   f[i] = density[i] - 1.0/(1.0 + exp((xmid - log(conc[i]))/scal));',
 parms = c(xmid = 0, scal = 1), 
 data = DNase1[,-1], language = "C")  


###################################################
### code chunk number 27: ccSolve.Rnw:514-520
###################################################
print(system.time(
for (i in 1:100)
 fm2 <- ccnls(fn = cfm.f95, data = DNase1[,-1], 
   par = c(xmid = 0, scal = 1))
))
summary(fm2)


###################################################
### code chunk number 28: ccSolve.Rnw:539-558
###################################################
Opt.R <- function(p) {
  with (DNase1[,-1],
   sum((density - 1/(1 + exp((p[1] - log(conc))/p[2])))^2)  
  ) 
}
print(system.time(for (i in 1:100)
A <- optim (fn = Opt.R, par = c(xmid = 0, scal = 1), 
  method = "CG")))

Opt.f95 = "f = sum((density - 1.d0/(1.d0 + dexp((xmid - dlog(conc))/scal)))**2)"

cOpt.f95 <- compile.optim(func = Opt.f95, par = c(xmid = 0, scal = 1), 
  data = DNase1[,-1])  
code(cOpt.f95)

print(system.time(for (i in 1:100)
AA <- ccoptim (fn = cOpt.f95, par = c(xmid = 0, scal = 1), 
  method = "CG", data = DNase1[,-1])))
AA


###################################################
### code chunk number 29: ccSolve.Rnw:564-574
###################################################
brown.R <- function(p) {
 sum((p[odd]^2)^(p[even]^2 + 1) + (p[even]^2)^(p[odd]^2 + 1))
}
npar <- 100  
p0 <- rnorm(npar, sd = 2)
n <- npar
odd <- seq(1, n, by = 2)
even <- seq(2, n, by = 2)
print(system.time(
  ans.opt <- optim(par = p0, fn = brown.R, method = "BFGS")))


###################################################
### code chunk number 30: ccSolve.Rnw:577-586
###################################################
brown.f95 <- "integer i
  f = 0.d0
  do i = 1, n-1, 2
     f = f + (x(i)**2)**(x(i+1)**2 + 1.d0) + (x(i+1)**2)**(x(i)**2 + 1.d0)
  enddo
"
ccbrown <- compile.optim(brown.f95)
print(system.time(
  ans.cc <- ccoptim(par = p0, fn = ccbrown, method = "BFGS")))


###################################################
### code chunk number 31: ccSolve.Rnw:595-596
###################################################
args(compile.ode)


###################################################
### code chunk number 32: ccSolve.Rnw:607-619
###################################################
parms <- c(rH = 2.82, A = 100, ks = 1)

parasite.R <- function (t, y, parms) {
  with (as.list(parms), {
   P <- y[1]
   H <- y[2]
   f <- A * P / (ks +H)
   Pnew <- H* (1-exp(-f))
   Hnew <- H * exp(rH*(1.-H) - f)
   list (c(Pnew, Hnew))   
  })
} 


###################################################
### code chunk number 33: ccSolve.Rnw:623-633
###################################################
declaration <- "        double precision ff"
parasite.f90 <- "
        ff = A * P / (ks + H)
        dP = H * (1.d0 - exp(-ff))
        dH = H * exp (rH * (1.d0 - H) - ff)
"
parms <- c(rH = 2.82, A = 100, ks = 15)
yini <- c(P = 0.5, H = 0.5)
cParasite <- compile.ode(func = parasite.f90, parms = parms, 
  y = yini, declaration = declaration, language = "Fortran")


###################################################
### code chunk number 34: ccSolve.Rnw:635-639
###################################################
system.time(out <- ode (func = parasite.R, y = yini, parms = parms, times = 0:1000,
    method = "iteration"))
system.time(outc <- ode (func = cParasite, y = yini, parms = parms, times = 0:1000,
     method = "iteration"))


###################################################
### code chunk number 35: it
###################################################
plot(out, outc, xlim = c(800, 1000), which = "P")


###################################################
### code chunk number 36: it
###################################################
plot(out, outc, xlim = c(800, 1000), which = "P")


###################################################
### code chunk number 37: ccSolve.Rnw:662-683
###################################################
declaration <- "double precision :: Ll, Lr, xb, yb" 
caraxis.f95 <- "
    yb = r * sin(w * t)
    xb = sqrt(L * L - yb * yb)
    Ll = sqrt(xl**2 + yl**2)
    Lr = sqrt((xr - xb)**2 + (yr - yb)**2)
        
    dxl = ul
    dyl = vl
    dxr = ur 
    dyr = vr
        
    dul  = (L0-Ll) * xl/Ll      + 2.0 * lam2 * (xl-xr) + lam1*xb
    dvl  = (L0-Ll) * yl/Ll      + 2.0 * lam2 * (yl-yr) + lam1*yb - k * g
               
    dur  = (L0-Lr) * (xr-xb)/Lr - 2.0 * lam2 * (xl-xr)
    dvr  = (L0-Lr) * (yr-yb)/Lr - 2.0 * lam2 * (yl-yr) - k * g
        
    dlam1 = xb * xl + yb * yl
    dlam2 = (xl - xr)**2 + (yl - yr)**2. - L * L
"


###################################################
### code chunk number 38: ccSolve.Rnw:687-699
###################################################
eps <- 0.01; M <- 10; k <- M * eps^2/2; 
L <- 1; L0 <- 0.5; r <- 0.1; w <- 10; g <- 1

parameter <- c(eps = eps, M = M, k = k, L = L, L0 = L0, 
               r = r, w = w, g = g)

yini <- c(xl = 0, yl = L0, xr = L, yr = L0,
          ul = -L0/L, vl = 0,
          ur = -L0/L, vr = 0,
          lam1 = 0, lam2 = 0)
ccaraxis <- compile.ode(caraxis.f95, parms = parameter, y = yini, 
  declaration = declaration)


###################################################
### code chunk number 39: ccSolve.Rnw:702-703
###################################################
index <- c(4, 4, 2)


###################################################
### code chunk number 40: ccSolve.Rnw:707-721
###################################################
Mass      <- diag(nrow = 10, 1)
Mass[5,5] <- Mass[6,6] <- Mass[7,7] <- Mass[8,8] <- M * eps * eps/2
Mass[9,9] <- Mass[10,10] <- 0
Mass

times <- seq(0, 3, by = 0.01)
outDLL <- daspk(y = yini, mass = Mass, times = times, func = ccaraxis,
                parms = parameter, nind = index)
p2 <- parameter; p2["r"] <- 0.2
outDLL2 <- daspk(y = yini, mass = Mass, times = times, func = ccaraxis,
                parms = p2, nind = index)
p2["r"] <- 0.05
outDLL3 <- daspk(y = yini, mass = Mass, times = times, func = ccaraxis,
                parms = p2, nind = index)


###################################################
### code chunk number 41: dllimp
###################################################
plot(outDLL, outDLL2, outDLL3, which = 1, type = "l", lwd = 2)


###################################################
### code chunk number 42: dllimp
###################################################
plot(outDLL, outDLL2, outDLL3, which = 1, type = "l", lwd = 2)


###################################################
### code chunk number 43: ccSolve.Rnw:758-773
###################################################
declaration <- "  double precision :: Min, oxicmin, anoxicmin
"  

cBiogeo.f95 <- "
  Min       = r*OM
  oxicmin   = Min*(O2/(O2+ks))
  anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2)

  dOM  = Flux - oxicmin - anoxicmin
  dO2  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)
  dSO4 = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4)
  dHS  = 0.5*anoxicmin   -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)

  SumS = SO4 + HS
"


###################################################
### code chunk number 44: ccSolve.Rnw:780-782
###################################################
pars <- c(D = 1, Flux = 100, r = 0.1, rox = 1,
          ks = 1, ks2 = 1, BO2 = 100, BSO4 = 10000, BHS = 0)


###################################################
### code chunk number 45: ccSolve.Rnw:784-788
###################################################
y <- c(OM = 1, O2 = 1, SO4 = 1, HS = 1)

cBiogeo <- compile.ode(func = cBiogeo.f95, parms = pars, y = y,
  outnames = "SumS", declaration = declaration)


###################################################
### code chunk number 46: ccSolve.Rnw:799-800
###################################################
code(cBiogeo)


###################################################
### code chunk number 47: ccSolve.Rnw:806-813
###################################################
ST <- stode (y = y, func = cBiogeo, parms = pars, 
  pos = TRUE, outnames = "SumS", nout = 1)
ST
pars["Flux"] <- 200
ST2 <- stode (y = y, func = cBiogeo, parms = pars, 
  pos = TRUE, outnames = "SumS", nout = 1)
ST2


###################################################
### code chunk number 48: ccSolve.Rnw:816-819
###################################################
out <- ode(y = y, func = cBiogeo, times = 0:50, parms = pars, 
  outnames = "sumS", nout = 1)
tail(out, n = 2)  


###################################################
### code chunk number 49: ccSolve.Rnw:844-845
###################################################
args(compile.bvp)


###################################################
### code chunk number 50: ccSolve.Rnw:879-884
###################################################
require(bvpSolve)

x       <- seq(0, 1, 0.01)
yini <- c(-1, NA, 0, 0, NA, NA)
yend <- c(1 , NA, 0, 0, NA, NA)


###################################################
### code chunk number 51: ccSolve.Rnw:891-903
###################################################
fswirl <- "
  eps = rpar(1)
  f(1) = Y(2)
  f(2) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(3) = Y(4)
 	f(4) = Y(5)
 	f(5) = Y(6)
  f(6) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl <- compile.bvp(func = fswirl, declaration = "double precision:: eps")
print(system.time(sol <- bvptwp(x = x, func = cswirl, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))


###################################################
### code chunk number 52: ccSolve.Rnw:907-915
###################################################
fswirl2 <- "
  eps = rpar(1)
  f(1) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(2) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl2 <- compile.bvp(func = fswirl2, declaration = "double precision:: eps")
print(system.time(sola <- bvpcol(x = x, func = cswirl2, order = c(2, 4), 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))


###################################################
### code chunk number 53: ccSolve.Rnw:921-933
###################################################
fswirl3 <- "
  f(1) = Y(2)
  f(2) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(3) = Y(4)
 	f(4) = Y(5)
 	f(5) = Y(6)
  f(6) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl3 <- compile.bvp(fswirl3, parms = c(eps = 0.1))

print(system.time(solb <- bvptwp(x = x, func = cswirl3, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))


###################################################
### code chunk number 54: ccSolve.Rnw:936-950
###################################################
print(system.time(sol2 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol[,1], yguess = t(sol[,-1]),
                  yini = yini, yend = yend, eps=0.001, epsini = 0.01,
                  parms=0.001)))

print(system.time(sol3 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol2[,1], yguess = t(sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))

print(system.time(sol3b <- bvptwp(x = x, func = cswirl3, 
                  xguess = sol2[,1], yguess = t(sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))


###################################################
### code chunk number 55: ccSolve.Rnw:953-956
###################################################
print(system.time(sol4 <- bvptwp(atol = 1e-5, x = x, func = cswirl, cond = TRUE, 
                  xguess = sol3[,1], yguess = t(sol3[,-1]),
                  yini = yini, yend = yend, eps = 5e-5 , parms=5e-5)))


###################################################
### code chunk number 56: swirl
###################################################
plot(sol, sol2, sol3, sol4)


###################################################
### code chunk number 57: swirl
###################################################
plot(sol, sol2, sol3, sol4)


###################################################
### code chunk number 58: ccSolve.Rnw:979-1019
###################################################
require(bvpSolve)
measel.R <- function(t, y, pars)  {
  bet <- 1575*(1+cos(2*pi*t))
  dy1 <- mu-bet*y[1]*y[3]
  dy2 <- bet*y[1]*y[3]-y[2]/lam
  dy3 <- y[2]/lam-y[3]/vv
  dy4 <- 0
  dy5 <- 0
  dy6 <-0
  
  list(c(dy1, dy2, dy3, dy4, dy5, dy6))
}

dmeasel.R <- function(t, y, pars) {
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

bound.R <- function(i, y, pars) {
  if ( i == 1 | i == 4) return(y[1] - y[4])
  if ( i == 2 | i == 5) return(y[2] - y[5])
  if ( i == 3 | i == 6) return(y[3] - y[6])  
}

dbound.R <- function(i, y, pars,vv) {
  if ( i == 1 | i == 4) return(c(1, 0, 0, -1 ,0, 0))
  if ( i == 2 | i == 5) return(c(0, 1, 0, 0, -1, 0))
  if ( i == 3 | i == 6) return(c(0, 0, 1, 0, 0, -1))
}


###################################################
### code chunk number 59: ccSolve.Rnw:1024-1039
###################################################
mu  <- 0.02
lam <- 0.0279
vv  <- 0.1

x <- seq (0, 1, by = 0.01)
yguess <- matrix(ncol = length(x), nrow = 6, data = 1)
rownames(yguess) <- paste("y", 1:6, sep = "")

print(system.time(
  solR <- bvptwp(func = measel.R, jacfunc = dmeasel.R, 
    bound = bound.R, jacbound = dbound.R, 
    xguess = x, yguess = yguess,
    x=x, leftbc = 3, ncomp = 6, 
    nmax = 100000, atol = 1e-4)
))


###################################################
### code chunk number 60: ccSolve.Rnw:1042-1106
###################################################
measel.f95 <- " 
   bet = 1575d0*(1.+cos(2*pi*x))
   f(1) = mu - bet*y(1)*y(3)
   f(2) = bet*y(1)*y(3) - y(2)/lam
   f(3) = y(2)/lam-y(3)/vv
   f(4) = 0.d0
   f(5) = 0.d0
   f(6) = 0.d0
"

dmeasel.f95 <- "
  bet = 1575d0*(1+cos(2*pi*x))
  df(1,1) =  -bet*y(3)
  df(1,3) =  -bet*y(1)

  df(2,1) =  bet*y(3)
  df(2,2) =  -1.d0/lam
  df(2,3) =  bet*y(1)

  df(3,2) = 1.d0/lam 
  df(3,3) = -1.d0/vv
"

bound.f95 <- "
  if ( i == 1 .OR. i == 4) g = (y(1) - y(4))
  if ( i == 2 .OR. i == 5) g = (y(2) - y(5))
  if ( i == 3 .OR. i == 6) g = (y(3) - y(6))  
"

dbound.f95 <- "
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
cMeasel <- compile.bvp(func = measel.f95, jacfunc = dmeasel.f95, 
  bound = bound.f95, jacbound = dbound.f95, parms  = parms, 
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


###################################################
### code chunk number 61: mes
###################################################
plot(sol1, sol2)


###################################################
### code chunk number 62: mes
###################################################
plot(sol1, sol2)


###################################################
### code chunk number 63: ccSolve.Rnw:1134-1178
###################################################
require(deSolve)

chaos.R <- function(t, state, parameters) {
    list(
    c(-8/3 * state[1] + state[2] * state[3],
      -10 * (state[2] - state[3]),
      -state[1] * state[2] + 28 * state[2] - state[3]))
}

state <- c(xx = 1, yy = 1, zz = 1)
times <- seq(0, 200, 0.01)
print(system.time(
  out   <- vode(state, times, chaos.R, 0)
))

# --------------------------------full compiled code -------------------------
chaos.f95 <- " 
    f(1)     = -8.d0/3 * y(1) + y(2) * y(3)
    f(2)     = -10.d0 * (y(2) - y(3))
    f(3)     = -y(1) * y(2) + 28d0 * y(2) - y(3)
"
cChaos <- compile.ode(chaos.f95) 
print(system.time(
  cout   <- vode(state, times, func = cChaos, parms = 0)
))

# ----------------------- calling compiled code in R -------------------------

rchaos <- function(t, state, parameters) {
   list(cChaos$func(3, t, state, f = 1:3, 1, 1)$f)
}

print(system.time(
  cout2  <- vode(state, times, func = rchaos, parms = 0)
))

# ----------------------- bitwise compilation in R -------------------------
require(compiler)
bchaos <- cmpfun(chaos.R)

print(system.time(
  cout3  <- vode(state, times, func = bchaos, parms = 0)
))



