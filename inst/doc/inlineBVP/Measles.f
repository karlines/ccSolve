c
c     This file is part of the Test Set for BVP solvers
c     http://www.dm.uniba.it/~bvpsolvers/
c
c        Problem MEASLES
c        ODE of dimension 2
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=34
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~bvpsolvers/src/problems/Kidney_Model.f
c
c
c
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
      character*(7), intent(out) ::problm
      character*(3), intent(out) :: type
      integer, intent(out) ::  neqn,ipar(1), nlbc
      double precision, intent(out) :: aleft, aright, rpar(4)
      logical, intent(out) :: numjac, numbcjac,linear
      double precision :: pi

      fullnm = 'Problem Measles'
      problm = 'Measles'
      type   = 'BVP'
      neqn   = 6
      nlbc = 3
      aleft   = 0.0d0
      aright  = 1.0d0
      numjac = .true.
      numbcjac = .true.
      linear = .false.

      rpar(1)  =  0.5
      PI=4.0D0*DATAN(1.0D0)
      rpar(2) = pi
      ipar(1)  = 0

      return
      end

c
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

c-----------------------------------------------------------------------
      subroutine init(neqn,yval0, nmsh, givmsh,givey,
     +        xguess,yguess,nudim)
      integer, intent(in) :: neqn
      integer, intent(out) :: nmsh
      double precision, intent(out) :: xguess(*),yguess(nudim,*)
      double precision, intent(out) :: yval0(1:neqn)
      logical, intent(out) :: givmsh, givey

      givmsh = .true.
      givey = .true.
      yval0(1:neqn) = 0.0d0
      nmsh = 11

c---- dummy values, not used if givmsh=FALSE, givey=false
      xguess(1) = 0
      do i=2,nmsh
        xguess(i)=xguess(i-1)+0.1d0
      end do
      xguess(nmsh)= 1

      DO i=1,neqn
        yguess(i,1:nmsh) = 0.5d0
      end do

      return
      end


c-----------------------------------------------------------------------
      subroutine feval(neqn,t,y,f,rpar,ipar)
      integer, intent(in) :: neqn,ipar(*)
      double precision, intent(in):: t,y(1:neqn),rpar(*)
      double precision, intent(out):: f(1:neqn)
      double precision :: S,N,L,I, mu, lambda, eta, B0, bt,pigreco
	  
c      Y(1)=S/N
c      y(2)=L/N
c      Y(3)=I/N
	  mu=0.02d0
      lambda=0.0279d0
      eta=0.01d0
      B0=1575d0
      pigreco = rpar(2)
      bt=B0*(1.0d0+(COS(pigreco))*t) 


c       call dblepr('eps', -1, rpar, 1)

      F(1)=mu-(bt*Y(1)*Y(2))
      F(2)=(bt*Y(1)*Y(3))-(Y(2)/lambda)
      F(3)=(Y(2)/lambda)-(Y(3)/eta)
      F(4)=0
      F(5)=0
      F(6)=0
      
      return
      end


c-----------------------------------------------------------------------
      subroutine bceval(i,neqn,y,bc,rpar,ipar)
      integer, INTENT(IN) :: i, neqn, ipar(*)
      double precision, INTENT(IN) :: y(neqn), rpar(*)
      double precision, INTENT(OUT) :: bc


      
        if (i .eq. 1) bc=Y(1)-Y(4)
        if (i .eq. 2) bc=Y(2)-Y(5)
        if (i .eq. 3) bc=Y(3)-Y(6)
        if (i .eq. 4) bc=Y(1)-Y(4)
        if (i .eq. 5) bc=Y(2)-Y(5)
        if (i .eq. 6) bc=Y(3)-Y(6)



      
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
      double precision lambda, lambda2


      

      return
      end
