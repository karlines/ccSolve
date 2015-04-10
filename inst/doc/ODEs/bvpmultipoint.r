## ==========================================================================
## a multipoint boundary value problem
## ==========================================================================
#
#  \begin{eqnarray*}
#     y_1' = (y_2 - 1)/2\\\
#     y_2' = (y_1*y_2 - x) / \mu\\
#     y_1(1) &=& 0\\
#     y_2(0.5) &=& 1
#  \end{eqnarray*}

multip <- "
f(1) = (y(2) - 1)/2
f(2) = (y(1)*y(2) - x)/mu
"

bound <- "
  if (i == 1) then
    g = y(2) -1    
  else 
    g = y(1)              
  endif  
"  

cmultip <- compile.bvp(func = multip, bound = bound, parms = c(mu = 0.1))

mu  <- 0.1
sol <- bvpcol(func = cmultip, x = seq(0, 1, 0.01), posbound = c(0.5, 1), parms = mu)

# check boundary value
sol[sol[,1] == 0.5,]
plot(sol)
