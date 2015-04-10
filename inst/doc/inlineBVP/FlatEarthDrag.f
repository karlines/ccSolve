c
c     This file is part of the Test Set for BVP solvers
c     http://www.dm.uniba.it/~bvpsolvers/
c
c        Problem Flat Earth Drag
c        ODE of dimension 8
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~bvpsolvers/testsetbvpsolvers/

c-----------------------------------------------------------------------
c
      integer function pidate()
      pidate = 20130531
      return 
      end
c-----------------------------------------------------------------------
      subroutine prob(fullnm,problm,type,
     +                neqn,aleft,aright,nlbc,ms,
     +                numjac,numbcjac,linear,rpar,ipar)
      implicit none
      character*(20), intent(out) :: fullnm
      character*(9), intent(out) ::problm
      character*(5), intent(out) :: type
      integer, intent(out) ::  neqn,ipar(1), nlbc, ms(8)
      double precision, intent(out) :: aleft, aright, rpar(2)
      logical, intent(out) :: numjac, numbcjac,linear



      fullnm = ' '
      problm = ' '
      type   = ' '

      fullnm = 'Problem Flat Earth Drag'
      problm = 'FlatEarthDrag'
      type   = 'BVP'

      neqn = 8
      nlbc = 4
      ms(1:8)  = 1

      aleft    = 0.0d0
      aright   = 1.0d0
      numjac   = .true.
      numbcjac = .true.
      linear   = .false.

c for this problem rpar is a dummy variable
      rpar(1)  =  0.0d0
c for this problem ipar is a dummy variable
      ipar(1)  = 0
      return
      end
c
c-----------------------------------------------------------------
c
      subroutine settolerances(neqn,rtol,atol,tolvec,ntol,ltol)
      integer, intent(in) :: neqn
      integer, intent(out) :: ntol, ltol(*)
      logical, intent(out) :: tolvec
      double precision, intent(in out) :: rtol(neqn), atol(neqn)

      tolvec  = .true.
      ntol = neqn
      ltol(1)=1

      DO i=2,ntol
         ltol(i) = i
         rtol(i) = rtol(1)
         atol(i) = atol(1)
      ENDDO
      

      return
      end
c---------------------------------------------------------------------
      subroutine setoutput(neqn,solref,printsolout,
     +                    nindsol,indsol)

      logical, intent(out) :: solref, printsolout
      integer, intent(in) :: neqn
      integer, intent(out) ::  nindsol
      integer, intent(out) :: indsol(*)

c the reference solution is available
      solref = .false.

c output file is required
      printsolout = .true.

c default values if printsolout is .true.
      nindsol = neqn
c only nindsol component of indsol are referenced
      do i=1,nindsol
          indsol(i) = i
      end do  


      return
      end
c-----------------------------------------------------
c Initialiser for parameter common block
c used only for the R interface
c
c
      SUBROUTINE initrparms(bvpparms)
      EXTERNAL bvpparms
c
c For Singularly Perturbed BVP (SPBVP) eps should be the
c first elements of parms and the first element of rpar
c
      DOUBLE PRECISION parms(1)
      COMMON / pars / parms

      CALL bvpparms(1, parms)

      RETURN
      END



c-----------------------------------------------------------------------
      subroutine init(neqn,ms,yval0, nmsh,givmsh,givey,
     +       xguess,yguess,dmguess,nudim)
      integer, intent(in) :: neqn,nudim,ms(neqn)
      integer, intent(out) :: nmsh
      double precision, intent(out) :: xguess(*),yguess(nudim,*)
      double precision, intent(out) :: dmguess(nudim,*)
      double precision, intent(out) :: yval0(1:neqn)
      logical, intent(out) :: givmsh, givey
      double precision :: h, tf_guess


      givmsh = .true.
      givey = .true.
      tf_guess = 100d0
     
      yval0(1:neqn) = 0.0d0
      yval0(6) = -1.0d0
      yval0(neqn)=tf_guess
      nmsh = 11

c---- dummy values, not used if givmsh=FALSE, givey=false
      h = 1.0d0/(nmsh-1)
      xguess(1) = 0
      do i=2,nmsh
        xguess(i)=xguess(i-1)+h
      end do
      xguess(nmsh)= 1

      DO i=1,neqn
        yguess(i,1:nmsh) = yval0(i)
      end do

      DO i=1,neqn
           dmguess(i,1:nmsh) = 0.0d0
      end do

      return
      end


c-----------------------------------------------------------------------
      subroutine feval(neqn,x,y,f,rpar,ipar)
      implicit none
      integer, intent(in) :: neqn,ipar(*)
      double precision, intent(in):: x,y(1:neqn)
      double precision, intent(in):: rpar(*)
      double precision, intent(out):: f(1:neqn)

      double precision :: fr,h,m,g_accel,vc,beta,eta,
     +      xbardot,ybardot,Vxbardot,Vybardot
     +      lambda_2_bar,lambda_3_bar,lambda_4_bar
    
      
      fr = 2.1*1d6
      h = 180000d0
      m = 60880d0
      g_accel = 9.80665d0
      vc = 1000d0*sqrt((3.986004*1d5)/(6378.14d0+(h/1000)))
      beta =180000d0/840d0
      eta = 1.225d0*0.5d0*7.069d0/2d0

      xbardot = Y(3)*(vc/h)
      ybardot = Y(4)*(vc/h)
      Vxbardot = (fr/vc*(-Y(6)/sqrt(Y(6)**2.0d0+Y(7)**2.0d0))
     +         -eta*exp(-Y(2)*beta)*Y(3)*sqrt(Y(3)**2.0d0+Y(4)**2.0d0)
     +          *vc)/m
      Vybardot = (fr/vc*(-Y(7)/sqrt(Y(6)**2.0d0+Y(7)**2.0d0))-
     +   eta*exp(-Y(2)*beta)*Y(4)*sqrt(Y(3)**2.0d0+Y(4)**2.0d0)*vc)
     +    /m-g_accel/vc
      if (sqrt(Y(3)**2.0d0 + Y(4)**2.0d0) == 0) then
         lambda_2_bar = 0.0d0
         lambda_3_bar = 0.0d0
         lambda_4_bar = -Y(5)*(vc/h)
      else
         lambda_2_bar = -(Y(6)*Y(3)+Y(7)*Y(4))*eta*beta*
     +        sqrt(Y(3)**2.0d0+Y(4)**2.0d0)*exp(-Y(2)*beta)*vc/m
         lambda_3_bar = eta*exp(-Y(2)*beta)*vc*(Y(6)*
     +        (2*Y(3)**2.0d0+Y(4)**2.0d0)+Y(7)*Y(3)*Y(4))/
     +         sqrt(Y(3)**2.0d0+Y(4)**2.0d0)/m
         lambda_4_bar =-Y(5)*vc/h+eta*exp(-Y(2)*beta)
     +          *vc*(Y(7)*(Y(3)**2.0d0)+2*Y(4)**2.0d0+
     +           Y(6)*Y(3)*Y(4))/sqrt(Y(3)**2.0d0+Y(4)**2.0d0)/m}
      end if
      
      tf = Y(8)


      F(1) = tf*xbardot
      F(2) = tf*ybardot
      F(3) = tf*Vxbardot
      F(4) = tf*Vybardot
      F(5) = tf*lambda_2_bar
      F(6) = tf*lambda_3_bar
      F(7) = tf*lambda_4_bar
      F(8) = 0

      

      return
      end

c-----------------------------------------------------------------------
      subroutine bceval(i,neqn,y,bc,rpar,ipar)
      integer, INTENT(IN) :: i, neqn, ipar(*)
      double precision, INTENT(IN) :: y(neqn), rpar(*)
      double precision, INTENT(OUT) :: bc
      

      double precision :: fr,h,m,g_accel,vc,beta,eta

      fr = 2.1*1d6
      h = 180000d0
      m = 60880d0
      g_accel = 9.80665d0
      vc = 1000d0*sqrt((3.986004*1d5)/(6378.14d0+(h/1000)))
      beta =180000d0/840d0
      eta = 1.225d0*0.5d0*7.069d0/2d0


      if (i .eq. 1) bc=Y(1)
      if (i .eq. 2) bc=Y(2)
      if (i .eq. 3) bc=Y(3)
      if (i .eq. 4) bc=Y(4)
        
      if (i .eq. 5) bc=Y(2)-1.0d0
      if (i .eq. 6) bc=Y(3)-1.0d0
      if (i .eq. 7) bc=Y(4)
      if (i .eq. 8) bc=(-sqrt(Y(6)**2.0d0+y(7)**2.0d0)*fr/m/vc
     +             -(y(6)*y(3))*eta*exp(-beta)*
     +              sqrt(y(3)**2.0d0)*vc/m-y(7)*g_accel/vc)*y(8)+ 1.0d0


      return
      end

c
c----------------------------------------------------------------------
c  Numerical Jacobian
c
      SUBROUTINE  jeval(R, X0, Y0, JF0, RPAR,IPAR)
c      SUBROUTINE  fnumjac(R, X0, Y0, JF0, RPAR,IPAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: R, IPAR(*)
      DOUBLE PRECISION, INTENT(IN) :: X0, Y0(R), RPAR(*)
      DOUBLE PRECISION,  INTENT(OUT) :: JF0(1:R,1:R)

      EXTERNAL feval
      DOUBLE PRECISION  :: DN(R), Y0N(R),YSAFE, DELT, X,FP(R),EPSPREC
      INTEGER :: I
                EPSPREC= epsilon(1.0d0)
                CALL feval(R,X0,Y0,FP,RPAR,IPAR)
                Y0N(1:R) = Y0(1:R)
                DO I=1,R
                  YSAFE=Y0(I)
                  DELT=SQRT(1.D-1*EPSPREC*MAX(1.D-5,ABS(YSAFE)))
                  Y0N(I)=YSAFE+DELT
                  X = X0
                  CALL feval(R,X,Y0N,DN,RPAR,IPAR)
                  JF0(1:R,I)=(DN(1:R)-FP(1:R))/DELT
                  Y0N(I)=YSAFE
               END DO

      RETURN
      END
c
c----------------------------------------------------------------------
c  computation of the numerical jacobian for the boundary conditions
c
      SUBROUTINE  dbceval(I,R,y,dbc,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  I, R
      DOUBLE PRECISION, INTENT(OUT)::dbc(R)
      DOUBLE PRECISION, INTENT(IN ):: y(R)
      INTEGER, INTENT(IN):: ipar(*)
      DOUBLE PRECISION, INTENT(IN) :: rpar(*)

      DOUBLE PRECISION :: DN, Y0N(R),YSAFE, DELT, BC, EPSPREC
      INTEGER :: J
      EXTERNAL bceval
      EPSPREC= epsilon(1.0d0)


                Y0N(1:R) = Y(1:R)
                call bceval(i,R,y,bc,rpar,ipar)
                DO J=1,R
                  YSAFE=Y(J)
                  DELT=SQRT(1.D-1*EPSPREC*MAX(1.D-5,ABS(YSAFE)))
                  Y0N(J)=YSAFE+DELT
                  call bceval(i,R,y0N,DN,rpar,ipar)
                  DBC(J)=(DN-BC)/DELT
                  Y0N(J)=YSAFE
                END DO


      RETURN
      END


c-----------------------------------------------------------------------
      subroutine solut(neqn,x,y,nmsh,rpar,ipar)
      integer, intent(in) ::  neqn, ipar(*)
      double precision, intent(in) :: x(nmsh),rpar(*)
      double precision, intent(out) :: y(1:neqn,1:nmsh)

      double precision :: lambda, lambda2
      common /pars/ lambda


      return
      end
