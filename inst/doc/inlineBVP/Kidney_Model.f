c
c     This file is part of the Test Set for BVP solvers
c     http://www.dm.uniba.it/~bvpsolvers/
c
c        Problem VAN DER POL
c        ODE of dimension 2
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=34
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~bvpsolvers/src/problems/bvpT1.f
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
      character*(8), intent(out) ::problm
      character*(3), intent(out) :: type
      integer, intent(out) ::  neqn,ipar(1), nlbc
      double precision, intent(out) :: aleft, aright, rpar(4)
      logical, intent(out) :: numjac, numbcjac,linear

      fullnm = 'Problem Kidney Model'
      fullnm = fullnm(1:13)
      problm = 'Kidney_Model'
      problm = problm(1:5)
      type   = 'BVP'
      neqn   = 14
      nlbc = 8
      aleft   = 0.0d0
      aright  = 1.0d0
      numjac = .false.
      numbcjac = .false.
      linear = .true.

      rpar(1)  =  1d-4
      ipar(1)  = 0
      rpar(3) = 10
      rpar(4)= 0
      return
      end

c
c
      subroutine settolerances(neqn,rtol,atol,tolvec,ntol,ltol)
      integer, intent(in) :: neqn
      integer, intent(out) :: ntol, ltol(2)
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
      solref = .true.  

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
      subroutine init(neqn,yval0, nmsh, givmsh,givey,xguess,yguess)
      integer, intent(in) :: neqn
      integer, intent(out) :: nmsh
      double precision, intent(out) :: xguess(*),yguess(neqn,*)
      double precision, intent(out) :: yval0(1:neqn)
      logical, intent(out) :: givmsh, givey

       write(*,*) 'input inside'
 

      givmsh = .false.
      givey = .false.
      yval0(1) = 0.0d0
      yval0(2) = 0.0d0

c---- dummy values, not used if givmsh=FALSE, givey=false
      xguess(1) = 0
      yguess(1,1) = 0
      nmsh = 10

      return
      end


c-----------------------------------------------------------------------
      subroutine feval(neqn,x,y,f,rpar,ipar)
      integer, intent(in) :: neqn,ipar(*)
      double precision, intent(in):: x,y(1:neqn),rpar(*)
      double precision, intent(out):: f(1:neqn)
      double precision :: eps, h1v, h3v, h2v, h4v, h5v,c11, F1v, F2v,
     +		J21, k1, c12, c21, c22, c31, c32, F4v, c41, c42, c51, c52, 
     +		c61,c62, J3v, J4v, J5v, J2v, J1v, f3v, f5v, f6v, 
     +		J41, J42,J12, J22, J32, J11
	  
   	  c12=Y(1)
      c21=Y(2)
      c22=Y(3)
      c31=Y(4)
      c32=Y(5)
      k1=Y(6)
      F4v=Y(7)
      c41=Y(8)
      c42=Y(9)
      c51=Y(10)
      c52=Y(11)
      c61=Y(12)
      c22=Y(13)
      
      
	  h1v= RPAR(3)
      h2v= RPAR(4)      
      h3v=RPAR(3)
      h2v=h4v
      h2v=h5v

      c11=20*Y(1)
      F1v=0.05/Y(1)
      F2v= -0.05/Y(3)


      if ((x.ge.0) .AND. (x.le.0.4)) then 
        J21=1.8
        J32=0
      else if ((x.gt.0.4) .AND. (x.lt.0.5)) then 
        J21=1.8+(-18+100*(c12-c41)*(x-0.4))
        J32=0.1*(c32-c42)*(x-0.4)
      else if ((x.ge.0.5) .and. (x.le.1)) then 
        J21=10*(c21-c41)
        J32=0.01*(c32-c42)
      end if
      
      
	  J3v= h3v*((c41-c31)+(c42-c32))
      J4v=-(J1v+J2v+J3v+J5v)
      J1v=h1v*(c42-c11)+(c42-c12)
      J2v=h2v*(c41-c21)+(c41-c22)     
      f3v=k1/c31
      f5v=5
      f6v=0.05/c62
      J11=0
      J12=0
      J22=0
      J31=0
      J51=1000*(c51-c41)
      J52=1000*(c52-c42)
      J61=(0.75*c61)/(1+c61)
      J62=0
      J41=-(J11+J21+J31+J51)
      J42=- (J12+J22+J32+J52)
c       call dblepr('eps', -1, rpar, 1)

      F(1)=(20*h1v*(c21)**2)*(c41+c42-21*c12)
      F(2)=20*Y(3)*J21
      F(3)=0
      F(4)=((h3v/k1)*(c31)**2)*(c41+c42-c31-c32)
      F(5)=(c31/k1)*((J3v*c32)-J32)
      F(6)=0
      F(7)=-J4v
      F(8)=(1/F4v)*(J4v*c41-J41)
      F(9)=(1/F4v)*(J4v*c42-J42)
      F(10)= -200*(c51-c41)
      F(11)= -200*(c52-c42)
      F(12)= (20*c62)*(J62*c61-J61)
	  F(13)= 20*(c62**2)
      F(14)=0
      return
      end


c-----------------------------------------------------------------------
      subroutine bceval(i,neqn,y,bc,rpar,ipar)
      integer, INTENT(IN) :: i, neqn, ipar(*)
      double precision, INTENT(IN) :: y(neqn), rpar(*)
      double precision, INTENT(OUT) :: bc
      c12=Y(1)
      c21=Y(2)
      c22=Y(3)
      c31=Y(4)
      c32=Y(5)
      k1=Y(6)
      F4v=Y(7)
      c41=Y(8)
      c42=Y(9)
      c51=Y(10)
      c52=Y(11)
      c61=Y(12)
      c22=Y(13)

      
        if (i .eq. 1) bc=c12-0.05
        if (i .eq. 2) bc=c51-1
        if (i .eq. 3) bc=c52-0.05
        if (i .eq. 4) bc=fv4+5          
        if (i .eq. 5) bc=(c31-20*k1*c31)-0
        if (i .eq. 6) bc=c22-c62
        if (i .eq. 7) bc=c61-c21
        if (i .eq. 8) bc=c31-F(14)
c con 1
        if (i .eq. 9) bc=c12-c22
        if (i .eq. 10) bc=c21-20*c12
        if (i .eq. 11) bc=c41-c51
        if (i .eq. 12) bc=c42+c52          
        if (i .eq. 13) bc=(c61-20*k1*c62)-0
        if (i .eq. 14) bc=c31-F(14)




      
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
c
c
c
      lambda = rpar(1)
      lambda2 = lambda**(1.0d0/2.0d0)
      DO i=1,nmsh

      y(1,i) = (exp(-X(i)/sqrt(lambda))-exp(-(2d0-X(i))/sqrt(lambda)))/
     +           (1.d0-exp(-2.d0/sqrt(lambda)))
      y(2,i) = (1./(lambda2*exp(X(i)/lambda2))+ exp((X(i) - 2.0d0)/
     +  lambda2)/lambda2)/(1/exp(2.0d0/lambda2) - 1)


      END DO

      return
      end
