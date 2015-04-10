c
c     This file is part of the Test Set for BVP solvers
c     http://www.dm.uniba.it/~bvpsolvers/
c
c        Problem Flat Earth
c        ODE of dimension 7
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
     +                neqn,aleft,aright,nlbc,
     +                numjac,numbcjac,linear,rpar,ipar)
      implicit none
      character*(20), intent(out) :: fullnm
      character*(9), intent(out) ::problm
      character*(5), intent(out) :: type
      integer, intent(out) ::  neqn,ipar(1), nlbc
      double precision, intent(out) :: aleft, aright, rpar(2)
      logical, intent(out) :: numjac, numbcjac,linear

      double precision :: EPS
      common /pars/ EPS
c
c  For Singularly perturbed BVP (SPBVP) eps should be the
c first elements of parms and the first element of rpar
c

      fullnm = 'Problem Flat Earth'
      problm = 'FlatEarth'
      type   = 'BVP'

      neqn = 7
      nlbc = 4
      aleft    = 0.0d0
      aright   = 1.0d0
      numjac   = .true.
      numbcjac = .true.
      linear   = .false.

c
c      EPS is an input parameter in report.f
c
c rpar(1) may contains epsin, the starting value for the ontinuation codes
c rpar(1) = 0 means that epsin is the default value
      rpar(1)  =  0.5d0
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
      subroutine init(neqn,yval0, nmsh,givmsh,givey,
     +       xguess,yguess,nudim)
      integer, intent(in) :: neqn,nudim
      integer, intent(out) :: nmsh
      double precision, intent(out) :: xguess(*),yguess(nudim,*)
      double precision, intent(out) :: yval0(1:neqn)
      logical, intent(out) :: givmsh, givey
      double precision :: h, tf_guess


      givmsh = .true.
      givey = .true.
      tf_guess = 700d0
     
      yval0(1:neqn) = 0.0d0
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

      return
      end


c-----------------------------------------------------------------------
      subroutine feval(neqn,x,y,f,rpar,ipar)
      implicit none
      integer, intent(in) :: neqn,ipar(*)
      double precision, intent(in):: x,y(1:neqn)
      double precision, intent(in):: rpar(*)
      double precision, intent(out):: f(1:neqn)

      double precision :: TF,g,vc,h,acc,  thrust2Weight
    
      VC=sqrt(3.986004d5/(6378.14d0+300d0))*1000d0 
	  h=300000d0
      g = 9.80665;
      thrust2Weight = 3
      Acc= Thrust2Weight*g
	  
      
     
      TF=Y(7)
      

      F(1)=(Y(3)*(vc/h))*TF
      F(2)=(Y(4)*(vc/h))*TF
      F(3)=(acc*(1/(ABS(vc)*sqrt(1+Y(6)**2))))*TF
      F(4)=(acc*(Y(6)/(ABS(vc)*sqrt(1+Y(6)**2)))- (g/vc))*TF
      F(5)=0
      F(6)=(-Y(5)*(vc/h))*TF
      F(7)=0
      

      return
      end
c-----------------------------------------------------------------------
c      subroutine jeval(neqn,x,y,dfdy,rpar,ipar)
c      integer, intent(in) ::  neqn,ipar(*)
c      double precision, intent(in):: x,y(1:neqn),rpar(*)
c      double precision, intent(out):: dfdy(1:neqn,1:neqn)

c      double precision :: EPS
c      common /pars/ EPS


c      dfdy(1,1) = 0d0
c      dfdy(1,2) = 1d0
c      dfdy(2,1) = 1d0/EPS
c      dfdy(2,2) = 0d0

c      return
c      end

c-----------------------------------------------------------------------
      subroutine bceval(i,neqn,y,bc,rpar,ipar)
      integer, INTENT(IN) :: i, neqn, ipar(*)
      double precision, INTENT(IN) :: y(neqn), rpar(*)
      double precision, INTENT(OUT) :: bc
      double precision :: h, VC
      VC=sqrt(3.986004d5/(6378.14d0+300d0))*1000d0 
	  h=300000d0
      
      if (i .eq. 1) bc=Y(1)
      if (i .eq. 2) bc=Y(2)
      if (i .eq. 3) bc=Y(3)
      if (i .eq. 4) bc=Y(4)
        
      if (i .eq. 5) bc=Y(2)-h/h
      if (i .eq. 6) bc=Y(3)-VC/VC
      if (i .eq. 7) bc=Y(4)
                


      return
      end

c      subroutine dbceval(i, ncomp,y, dbc, rpar,ipar)
c      integer, intent(in) ::  i,ncomp
c      double precision, INTENT(IN) ::  y(ncomp)
c      double precision, INTENT(OUT) :: dbc(ncomp)

c      double precision :: EPS
c      common /pars/ EPS

c      dbc(1)=1.D0
c      dbc(2)=0.D0

c      return
c      end

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
c     SUBROUTINE  BCNUMJAC(I,R,y,dbc,RPAR,IPAR)

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
