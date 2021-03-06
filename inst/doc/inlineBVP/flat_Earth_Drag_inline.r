## ---------------------------------------------------------------------------
## the flat earth with drag boundary value problem,
## implemented with inline F95 code
## ---------------------------------------------------------------------------

# -------------------------
# parameter declarations
# -------------------------

h <- 180000

parms <- c(fr = 2.1e6, h = h, m = 60880, g_accel = 9.80665,
      vc = 1000*sqrt((3.986004e5)/(6378.14 + (h/1000))), beta =180000/840,
      eta = 1.225*0.5*7.069/2)

x <- seq (0, 1, length.out = 100)
yg <- c(0, 0, 0, 0, 0, -1, 0, 100)   # guess


# -------------------------
# derivative function
# -------------------------

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
# -------------------------
# boundary function
# -------------------------

bflatearth <- "

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

# -------------------------
# compile it- 
# passing parameters
# -------------------------

cflatearth <- compile.bvp(func = fflatearth, bound = bflatearth, parms = parms, 
  declaration = "double precision :: tf, dx, dy, dVx, dVy, dlambda_2, dlambda_3, dlambda_4")



x <- seq (0, 1, length.out = 100)
yg <- c(0, 0, 0, 0, 0, -1, 0, 100)   # guess

print (system.time(
csol <-  bvptwp(func = cflatearth, x = x, xguess = x, leftbc = 4,
  yguess = matrix (nrow = 8, ncol = length(x), data = yg), 
  parms = parms)
))




print (system.time(
csol2 <-  bvpcol(func = cflatearth, x = x, xguess = x, leftbc = 4, nmax = 5e5,
  yguess = matrix (nrow = 8, ncol = length(x), data = yg), parms = parms)
))

plot(csol, csol2)
