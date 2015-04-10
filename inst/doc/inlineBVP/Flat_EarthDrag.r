## ---------------------------------------------------------------------------
## the flat earth with drag boundary value problem,
## implemented with inline F95 code AND in R
## ---------------------------------------------------------------------------

require(bvpSolve)

x <- seq (0, 1, length.out = 100)
yg <- c(0, 0, 0, 0, 0, -1, 0, 100)   # guess
  
f       <- 2.1*10^6
h       <- 180000
m       <- 60880
g_accel <- 9.80665
vc      <- 1000*sqrt((3.986004*10^5)/(6378.14+(h/1000)))
beta    <- 180000/840
eta     <- 1.225*0.5*7.069/2

earthdrag <- function(tau, X, parms) { 
    
  dx <- X[3] * (vc/h)
  dy <- X[4] * (vc/h)
  dVx <- (f/vc * (-X[6]/sqrt(X[6]^2+X[7]^2))
    - eta*exp(-X[2]*beta)*X[3]*sqrt(X[3]^2+X[4]^2)*vc)/m
  dVy <- (f/vc * (-X[7]/sqrt(X[6]^2+X[7]^2)) 
    - eta*exp(-X[2]*beta)*X[4]*sqrt(X[3]^2+X[4]^2)*vc)/m - g_accel/vc

  if (sqrt(X[3]^2 + X[4]^2) == 0){
    dlambda_2 <- 0
    dlambda_3 <- 0
    dlambda_4 <- -X[5]*(vc/h) 
  } else {
    dlambda_2 <- -(X[6]*X[3]+X[7]*X[4]) *eta*beta*sqrt(X[3]^2 + X[4]^2)*exp(-X[2]*beta)*vc/m
    dlambda_3 <- eta*exp(-X[2]*beta)*vc*(X[6]*(2*X[3]^2+X[4]^2)+X[7]*X[3]*X[4])/sqrt(X[3]^2+X[4]^2)/m
    dlambda_4 <- -X[5]*vc/h+eta*exp(-X[2]*beta)*vc*(X[7]*(X[3]^2)+2*X[4]^2+X[6]*X[3]*X[4])/sqrt(X[3]^2+X[4]^2)/m
  }
    
  f <- (X[8] * c(dx, dy, dVx, dVy, dlambda_2, dlambda_3, dlambda_4, 0))
  return(list(f))
}

bound <- function  (i, y, parms) {
          
  if (i == 1) return (y[1])
  if (i == 2) return (y[2])
  if (i == 3) return (y[3])
  if (i == 4) return (y[4])
  if (i == 5) return (y[2]-1.0e0)
  if (i == 6) return (y[3]-1.0e0)
  if (i == 7) return (y[4])
  if (i == 8) return ((-sqrt(y[6]^2+y[7]^2)*f/m/vc-
        (y[6]*y[3])*eta*exp(-beta)*sqrt(y[3]^2)*vc/m
        - y[7]*g_accel/vc)*y[8]+ 1)
    
}
  
print (system.time(
sol <-  bvptwp(func = earthdrag, x = x, xguess = x, leftbc = 4,
  yguess = matrix (nrow = 8, ncol = length(x), data = yg), bound = bound)
))



# all this very sensitive to double precision...
parms <- c(
  fr      = 2.1*10^6,
  h       = 180000,
  m       = 60880,
  g_accel = 9.80665,
  vc      = 1000*sqrt((3.986004*10^5)/(6378.14+(180000/1000))),
  beta    = 180000/840,
  eta     = 1.225*0.5*7.069/2)

fflatearth <- "
      dx  = Y(3)*(vc/h)
      dy  = Y(4)*(vc/h)
      dVx = (fr/vc*(-Y(6)/sqrt(Y(6)**2.0d0+Y(7)**2.0d0))                      &
             -eta*exp(-Y(2)*beta)*Y(3)*sqrt(Y(3)**2.0d0+Y(4)**2.0d0)          &
              *vc)/m
      dVy = (fr/vc*(-Y(7)/sqrt(Y(6)**2.0d0+Y(7)**2.0d0))-                     &
         eta*exp(-Y(2)*beta)*Y(4)*sqrt(Y(3)**2.0d0+Y(4)**2.0d0)*vc)           &
          /m-g_accel/vc
      if (sqrt(Y(3)**2.0d0 + Y(4)**2.0d0) == 0) then
         dlambda_2 = 0.0d0
         dlambda_3 = 0.0d0
         dlambda_4 = -Y(5)*(vc/h)
      else
         dlambda_2 = -(Y(6)*Y(3)+Y(7)*Y(4))*eta*beta*                         &
              sqrt(Y(3)**2.0d0+Y(4)**2.0d0)*exp(-Y(2)*beta)*vc/m
         dlambda_3 = eta*exp(-Y(2)*beta)*vc*(Y(6)*                            &
              (2*Y(3)**2.0d0+Y(4)**2.0d0)+Y(7)*Y(3)*Y(4))/                    &
               sqrt(Y(3)**2.0d0+Y(4)**2.0d0)/m
         dlambda_4 =-Y(5)*vc/h+eta*exp(-Y(2)*beta)                            &
                *vc*(Y(7)*(Y(3)**2.0d0)+2*Y(4)**2.0d0+                        &
                 Y(6)*Y(3)*Y(4))/sqrt(Y(3)**2.0d0+Y(4)**2.0d0)/m 
      end if
      
      tf = Y(8)


      F(1) = tf*dx
      F(2) = tf*dy
      F(3) = tf*dVx
      F(4) = tf*dVy
      F(5) = tf*dlambda_2
      F(6) = tf*dlambda_3
      F(7) = tf*dlambda_4
      F(8) = 0.d0
"
fblatearth <- "

      if (i .eq. 1) g=Y(1)
      if (i .eq. 2) g=Y(2)
      if (i .eq. 3) g=Y(3)
      if (i .eq. 4) g=Y(4)
        
      if (i .eq. 5) g=Y(2)-1.0d0
      if (i .eq. 6) g=Y(3)-1.0d0
      if (i .eq. 7) g=Y(4)
      if (i .eq. 8) g=(-sqrt(Y(6)**2.0d0+y(7)**2.0d0)*fr/m/vc                &
                   -(y(6)*y(3))*eta*exp(-beta)*                              &
                    sqrt(y(3)**2.0d0)*vc/m-y(7)*g_accel/vc)*y(8)+ 1.0d0

"

cflatearth <- compile.bvp(func = fflatearth, bound = fblatearth, parms  = parms,
   declaration = "double precision :: tf,dx, dy, dVx, dVy, dlambda_2, dlambda_3, dlambda_4")



print (system.time(
csol <-  bvptwp(func = cflatearth, x = x, xguess = x, leftbc = 4, parms = parms,
  yguess = matrix (nrow = 8, ncol = length(x), data = yg))
))


