!-----------------------------------------------------------------------------------------------------------------
! some subroutines related to combination C(n,m)
!
!-----------------------------------------------------------------------------------------------------------------


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  this module is reconstruct by wenliang li, but the origin code is $ 
!  get from the internet. the origin code is .f which was convert to $ 
!   .f90 by using to_f90 by alan miller.                             $
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! Code converted using TO_F90 by Alan Miller
SUBROUTINE allnr(n, r, numcombin,jout, ifault)
implicit none

!        Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3

!        When called once, generates all possible combinations
!        from a group of N items.  Each combination (represented in j as
!        r ordered integers between 1 and n) is processed within allnr.

!        Parameters:-

!        n        integer             input:  The size of the group from which
!                                             the combinations are selected.

!        r        integer             input:  The size of each comination.

!        j        integer array(r)  workspace: Used by allnr to store
!                                              combinations.

!        ifault   integer            output:  Fault indicator, equal to:
!                                             0 if 1 le R le N;
!                                             1 otherwise.



INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: r
INTEGER, INTENT(IN OUT)                  :: numcombin
INTEGER, INTENT(OUT)                     :: jout(numcombin,r)
INTEGER, INTENT(OUT)                     :: ifault
INTEGER :: j(r)
integer :: kount,nmr,i,ip1,l,ipan

ifault = 1
IF (r < 1 .OR. r > n) RETURN
ifault = 0
kount = 0
nmr = n - r

!        Initialize J(1) to lower limit separately, since lower limit for
!        each index depends on lower limit for previous index

i = 1
j(1) = 1

!        Initialize indices for loops i=1,...,r to lower limits

1 IF (i == r) GO TO 3
ip1 = i + 1
DO  l = ip1, r
  j(l) = j(l - 1) + 1
END DO

!        Update the count (kount) of combinations and process the current
!        combination.  The call to Subroutine job may be replaced by
!        statements to process the current combination.

3 kount = kount + 1
!CALL job(n, r, j, kount)

!        Increment the first possible index (of loop i) among indices of
!        loops R, R-1,...,1

! PanShi:      Write(*,*) 'PanShi-CeShi-09-06-01',J,Kount
! PanShi:==============================================================
DO ipan=1,r
  jout(kount,ipan)=j(ipan)
END DO
! PanShi:==============================================================
i = r
4 IF (j(i) < nmr + i) GO TO 5
i = i - 1

!        Return after all indices have achieved their upper limits

IF (i <= 0) RETURN
GO TO 4

5 j(i) = j(i) + 1
GO TO 1
END SUBROUTINE allnr





SUBROUTINE ncr(n, r, ncomb, ier)
 
! Code converted using TO_F90 by Alan Miller

!     Calculate the number of different combinations of r objects out of n.

!     ier = 0 if no error is detected
!         = 1 if n < 1
!         = 2 if r < 0
!         = 3 if r > n
!         = 4 if nCr > 1.e+308, i.e. if it overflows.  In this case, the
!                natural log of nCr is returned.

!     Programmer: Alan.Miller @ cmis.csiro.au
!     Latest revision - 28 July 1988 (Fortran 77 version)

IMPLICIT NONE
INTEGER, PARAMETER      :: dp =kind(0.0d0)! SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: r
REAL (dp), INTENT(OUT)  :: ncomb
INTEGER, INTENT(OUT)    :: ier

INTEGER :: rr, i, nn

INTERFACE
  FUNCTION lngamma(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp =kind(0.0d0) ! SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION lngamma
END INTERFACE

IF (n < 1) THEN
  ier = 1
ELSE IF (r < 0) THEN
  ier = 2
ELSE IF (r > n) THEN
  ier = 3
ELSE
  ier = 0
END IF
IF (ier /= 0) RETURN

IF (r <= n-r) THEN
  rr = r
ELSE
  rr = n - r
END IF

IF (rr == 0) THEN
  ncomb = 1.0_dp
  RETURN
END IF

IF (rr > 25) THEN
  ncomb = lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1))
  IF (ncomb > 709._dp) THEN
    ier = 4
  ELSE
    ncomb = EXP(ncomb)
  END IF
  RETURN
END IF

ncomb = n
i = 1
nn = n
DO
  IF (i == rr) RETURN
  nn = nn - 1
  i = i + 1
  ncomb = (ncomb * nn) / REAL(i)
END DO

RETURN
END SUBROUTINE nCr

!PROGRAM test_nCr
Subroutine FactorialCnm(N,MR,ResultPan)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp =kind(0.0d0)! SELECTED_REAL_KIND(12, 60)

INTERFACE
  SUBROUTINE ncr(n, r, ncomb, ier)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp =kind(0.0d0)! SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: r
    REAL (dp), INTENT(OUT)  :: ncomb
    INTEGER, INTENT(OUT)    :: ier
  END SUBROUTINE ncr
END INTERFACE

INTEGER    :: n, r, ier,Mr
REAL (dp)  :: resultPan
r=Mr

!DO
!  WRITE(*, '(a)', ADVANCE='NO') ' Enter n, r : '
!  READ(*, *) n, r
  CALL nCr(n, r, resultPan, ier)
  IF (ier /= 0) THEN
    WRITE(*, *) ' Error, IER = ', ier
    IF (ier == 4) WRITE(*, '(a, f12.5)') ' Ln(nCr) = ', resultPan
  ELSE
!    WRITE(*, '(a, g16.8)') ' nCr = ', resultPan
  END IF
!END DO
Return
!STOP
!END PROGRAM test_nCr
End Subroutine FactorialCnm


FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

IMPLICIT NONE
INTEGER, PARAMETER    :: dp =kind(0.0d0)! SELECTED_REAL_KIND(12, 60)
REAL(dp), INTENT(IN)  :: z
REAL(dp)              :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
                              -1259.139216722289D0, 771.3234287757674D0, &
                              -176.6150291498386D0, 12.50734324009056D0, &
                              -0.1385710331296526D0, 0.9934937113930748D-05, &
                               0.1659470187408462D-06 /), zero = 0.D0,   &
                               one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
                               half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER          :: j

IF (z <= zero) THEN
  WRITE(*, *) 'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma


      subroutine getsign(iae,ibe,ndim,signdicp)
      implicit none
      integer,intent(in) :: ndim
      integer,intent(out) :: signdicp
      integer,dimension(ndim) :: iae,ibe,ibep
      integer idicp,isturn,jdicp,is,IsIncrease,ISmediaTemp
 
      ibep = ibe

    Is=0                                
    IsTurn=0                                          
    do idicp=1,ndim                                    
        IsIncrease=0                                   
        do jdicp=1,ndim                             
            If(iae(idicp)==ibep(jdicp)) Then
             is=is+1
             isincrease=isincrease+1          
               
            If(idicp/=jdicp) Then
              isturn=isturn+1
              ISmediaTemp=ibep(idicp)
              ibep(idicp)=ibep(jdicp) 
              ibep(jdicp)=ISmediaTemp 
            endIf               
           endIf 
           if(isincrease.eq.1) exit    
          enddo
                         
          If(isincrease>1) then 
                 write(*,*) 'error happen in subroutine get'
          endIf
        enddo
        signdicp=(-1)**isturn   
        if((ndim-is)/= 0) then 
         write(*,*) 'error happen in subroutine get'
         endIf  
      return
    end subroutine getsign


  subroutine get(iae,ibe,ndim,numdiffer)
    implicit none
    integer,intent(in) :: ndim
    integer,intent(in),dimension(ndim) :: iae,ibe
    integer,intent(out) :: numdiffer
    integer :: ie,is,isincrease,idicp,jdicp
    is=0                                                                       
    do idicp=1,ndim                                    
     isincrease=0                                   
     do jdicp=1,ndim                 
         if(iae(idicp)==ibe(jdicp)) Then
          is=is+1
          isincrease=isincrease+1 
         endif 
          if(isincrease.eq.1) exit 
         enddo
      enddo
      numdiffer=ndim-is
    return
   end subroutine get
!==================================================================
subroutine reg_mat(ca_daummy,ca_daummy_reg,nn,eps)
!==================================================================

  implicit none

  integer(4) :: lwork,nn,info,ip,iq,ir
  real(8)    :: eps,rwork(3*nn-2),ee(nn)
  complex(8) :: cs,ca_daummy(nn,nn),ca_daummy_reg(nn,nn),cwork(8*nn)

!---- regulalization of a matrix ----
  lwork=6*nn
  call zheev('v','u',nn,ca_daummy,nn,ee,cwork,lwork,rwork,info)

  if(info.ne.0) stop '***** ERROR IN REG_MAT *****'
!      write(*,*) (ee(ip),ip=1,nn)

!  eps=1.d-7
  do ip=1,nn
     ee(ip)=ee(ip)+eps*dexp(-dabs(ee(ip))/eps)
  enddo
!  write(*,*) (ee(ip),ip=1,nn),"<------"
     
  do ip=1,nn
  do iq=1,nn
     cs=dcmplx(0.d0,0.d0)
     do ir=1,nn
        cs=cs+ca_daummy(ip,ir)*dcmplx(ee(ir),0.d0)*dconjg(ca_daummy(iq,ir))
     enddo
     ca_daummy_reg(ip,iq)=cs
  enddo
  enddo
  
return
end subroutine reg_mat

!! DONE BY JUAN

!! This function sign_m1(m) gives the result of (-1)**m as a real(kind=k1)
Function sign_m1(m)
implicit none

!! INPUT
integer, intent(in) :: m !! exponent of (-1)

!! OUTPUT
real*8:: sign_m1 

If (mod(abs(m),2).eq.0) sign_m1=1.0d0
If (mod(abs(m),2).eq.1) sign_m1=-1.0d0

End Function sign_m1

!-------------------------------------------------------------------------
! This function calculates the wigner3j symbols using Racah formula
! 
!   _________________
!  /                  \
! |   a     b     c    |
! | alfa  beta  gamma  |
!  \ _________________/
!  
!-------------------------------------------------------------------------

  function wigner3j(a,b,c,alfa,beta,gamma)
    Implicit none
    Integer :: a, b, c
    Integer :: alfa, beta, gamma
    Integer :: t, i, tini, tfin
    Real*8 :: x, wigner3j

    wigner3j=0.0d0

      If (alfa+beta+gamma.ne.0) return
      If (c.lt.abs(a-b)) return
      If (c.gt.(a+b)) return
      If (a.lt.0.or.b.lt.0.or.c.lt.0) return
            
      tini=max(0,max(-c+a+beta,-c+b-alfa))
      tfin=min(a+b-c,min(a-alfa,b+beta))
     
      Do t=tini,tfin
         x=1.0d0
         
         Do i=1,t
            x=x*dsqrt(dble(i))/dble(i)
         EndDo

         Do i=t+1,a+b-c
            x=x*dsqrt(dble(i))/dble(i-t)
         EndDo

         Do i=1,c-a-beta+t
            x=x*dsqrt(dble(i))/dble(i)
         EndDo

         Do i=c-a-beta+t+1,c+b-a
            x=x*dsqrt(dble(i))/dble(i-c+a+beta-t)
         Enddo

         Do i=1,c-b+alfa+t
            x=x*dsqrt(dble(i))/dble(i)
         EndDo

         Do i=c-b+alfa+t+1,a-b+c
            x=x*dsqrt(dble(i))/dble(i-c+b-alfa-t)
         EndDo

         x=(-1.0d0)**(dble(t))*x
        
         wigner3j=wigner3j+x
      EndDo

      Do i=1,a-abs(alfa)
         wigner3j=wigner3j*dble(i)/dsqrt(dble(i))
      EndDo

      Do i=a-abs(alfa)+1, a
         wigner3j=wigner3j*dsqrt(dble(i))*dsqrt(dble(abs(alfa)+i))/dsqrt(dble(i))     
      EndDo

      Do i=1, b-abs(beta)
         wigner3j=wigner3j*dble(i)/dsqrt(dble(a+i))
      EndDo

      Do i=b-abs(beta)+1,b
         wigner3j=wigner3j*dsqrt(dble(i))*dsqrt(dble(abs(beta)+i))/dsqrt(dble(a+i))
      EndDo

      Do i=1, c-abs(gamma)
         wigner3j=wigner3j*dble(i)/dsqrt(dble(a+b+i))
      EndDo

      Do i=c-abs(gamma)+1,c
         wigner3j=wigner3j*dsqrt(dble(i))*dsqrt(dble(abs(gamma)+i))/dsqrt(dble(a+b+i))
      EndDo

      wigner3j=wigner3j/dsqrt(dble(a+b+c+1))

      wigner3j=wigner3j*(-1.0d0)**dble(a-b-gamma)

      return
    end Function wigner3j

!!========================================================================
!! Modified from NUMERICAL RECIPES
!! 
!! gauleg calculates the nodes and the weights of a Gauss-Legendre 
!! quadrature in a finite interval
!!
!!========================================================================

      SUBROUTINE gauleg(x1,x2,x,w,n)
        implicit none
!! INPUT
        INTEGER, intent(in):: n !! number of points of the quadrature
      REAL*8, intent(in):: x1,x2 !! lower and upper limit of the interval

!! OUTPUT    
      real*8, allocatable, intent(out) :: x(:),w(:) !! nodes and weights of the quadrature
      

!! Auxiliar variables
      
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      Real*8::  p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do i=1,m
        z=cos(acos(-1.0d0)*(i-.25d0)/(n+.5d0))
        Do
           p1=1.d0
           p2=0.d0
           do  j=1,n
              p3=p2
              p2=p1
              p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
           End do
           pp=n*(z*p1-p2)/(z*z-1.d0)
           z1=z
           z=z1-p1/pp
           if(abs(z-z1).lt.EPS) exit
        End Do
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
     End do
      return
    END SUBROUTINE gauleg
!!C  (C) Copr. 1986-92 Numerical Recipes Software ]k1">"@w.

!!==============================================================
!! Calculation of the associated Legendre polynomials with m>=0
!!==============================================================

      FUNCTION plgndr(l,m,x)
      INTEGER l,m
      REAL*8:: plgndr,x
      INTEGER i,ll
      REAL*8 fact,pll,pmm,pmmp1,somx2

      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) write(*,*) 'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do  i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
       End do
    endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
         End do
          plgndr=pll
        endif
      endif
      return
    END FUNCTION plgndr
!!C  (C) Copr. 1986-92 Numerical Recipes Software ]k1">"@w.

!!==============================================================
!! Calculation of the Spherical Harmonics in radians
!!==============================================================

    Function y_lm(l,m,theta,phi)
      !! Spherical Harmonic Y_{lm}(theta,phi)
      complex*16 y_lm
      integer l, m
      real*8 theta, phi
      !! AUXILIAR VARIABLES
      real*8 aux, iaux !! stores the factorial
      integer i, j, k !! indexes for the loops
      
      Interface
         FUNCTION plgndr(l,m,x)
           INTEGER l,m
           REAL*8:: plgndr,x
         End FUNCTION plgndr
         Function sign_m1(m)
           integer, intent(in) :: m 
           real*8:: sign_m1 
         End Function sign_m1
      End Interface

      iaux=1

      Do i=l-abs(m)+1,l+abs(m)
         iaux=iaux*i
      End Do

      If (m.ne.0) aux=1.0d0/dble(iaux)
      If (m.eq.0) aux=1.0d0

      aux=(2.0d0*dble(l)+1.0d0)/(4.0d0*acos(-1.0d0))*aux

      aux=sqrt(aux)

      y_lm=aux*plgndr(l,abs(m),cos(theta))*(cos(dble(abs(m))*phi)+cmplx(0.0d0,1.0d0)*sin(dble(abs(m))*phi))
      
      If (m.lt.0) y_lm=cmplx(sign_m1(m))*conjg(y_lm)
      
      return

    End Function y_lm
