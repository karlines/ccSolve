## =========================================================================
## A boundary value differential algebraic equation
## =========================================================================

# last equation is algebraic equation:
daebvp <- "
  f(1) = (ks + y(2) - sin(x))*y(4) + cos(x)
  f(2) = cos(x)
  f(3) = y(4)
  f(4) = (y(1) - sin(x))*(y(4) - exp(x))
"
bounddae <- "
  if (i == 1) then
    g = (y(1) - sin(0.d0))
  else if (i == 2) then
    g = y(3) - 1
  else if (i == 3) then
    g = y(2) - sin(1.d0)
  else
    g = (y(1) - sin(1.d0))*(y(4) - exp(1.d0))
  endif  
"
cdaebvp <- compile.bvp(func = daebvp, bound = bounddae, parms = c(ks = 1e-4))

x <- seq(0, 1, by = 0.01)
mass <- diag(nrow = 4)  ; mass[4, 4] <- 0

out <- bvpcol (func = cdaebvp, x = x, atol = 1e-10, rtol = 1e-10,
               parms = 1e-4, ncomp = 4, leftbc = 2,
               dae = list(index = 2,  nalg = 1)) 

# the analytic solution
ana <- cbind(x, "1" = sin(x), "2" = sin(x), "3" = 1, "4" = 0, res = 0)
plot(out, obs = ana)
