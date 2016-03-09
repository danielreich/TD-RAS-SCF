
C***********************************************************************
C 
C                             LA-LIB                                 
C
C  Library module containing general routines for matrix and vector   
C  operations i.e. those not classifiable as matrix*matrix (mm), 
C  matrix*vector (mv) etc.
C
C Contents:
C   mattens:   multiplies a matrix with a tensor to give a matrix
C              b(i,j)=t(i,j,k)*a(i,k)
C   mktens:    makes a tensor from the product of a complex matrix with
C              a complex vector
C              t(i,j,k)=v(i)*m(j,k)
C   mktens_zd: makes a tensor from the product of a complex matrix with
C              a real vector
C              t(i,j,k)=v(i)*m(j,k)
C   normspf:   Normalizes a number of complex vectors (e.g. all 
C              single-particle functions of a specified degree of
C              freedom) simultaneously
C   normspfr:  Normalizes a number of real vectors simultaneously
c   vrminmax:  finds min and max of a real*4 vector simultaneously
c   vdminmax:  finds min and max of a real*8 vector simultaneously
C
C***********************************************************************


c-----------------------------------------------------------------------
c Library subroutine mattens
c
c Multiplication of a complex matrix with a complex tensor of third
c order:
c           B(i,j) = Sum_k[ T(i,j,k)*A(i,k) ]
c           i=1,...,p;  j=1,...,q;  k=1,...,r
c
c Input-variables:  A     : complex matrix
c                   T     : complex tensor of third order
c                   p,q,r : dimensions of the different arrays
c                   add   : if true, product is added to B
c Output-variables: B     : resulting complex matrix
c-----------------------------------------------------------------------

      subroutine mattens (B,A,T,p,q,r,add)

      implicit none

      integer    i,j,k,p,q,r
      complex*16 A(p,r),B(p,q),T(p,q,r)
      logical    add

        
      if (.not. add) then
         do j=1,q
            do i=1,p
               B(i,j)=(0.0d0,0.0d0)
            enddo
         enddo
      endif

      do j=1,q
         do k=1,r
            do i=1,p
               B(i,j)=B(i,j)+T(i,j,k)*A(i,k)
            enddo
         enddo
      enddo

      return
      end
      
c-----------------------------------------------------------------------
c Library subroutine mktens
c
c Makes a tensor of third order by elementwise multiplying a vector with
c a matrix:
c           T(i,j,k) = V(i)*M(j,k)   (no summation)
c           i=1,...,a;  j=1,...,b;  k=1,...,c
c
c Input-variables:  V     : complex vector
c                   M     : complex matrix
c                   a,b,c : dimensions of the different arrays
c                   add   : if true, result is added to T
c Output-variables: T     : resulting complex tensor of third order
c-----------------------------------------------------------------------

      subroutine mktens (T,V,M,a,b,c,add)

      implicit none

      integer    i,j,k,a,b,c
      complex*16 T(a,b,c),V(a),M(b,c)
      logical    add

      if (.not. add) then
         do k=1,c
            do j=1,b
               do i=1,a
                  T(i,j,k)=(0.0d0,0.0d0)
               enddo
            enddo
         enddo
      endif
      do k=1,c
         do j=1,b
            do i=1,a
               T(i,j,k)=T(i,j,k)+V(i)*M(j,k)
            enddo
         enddo
      enddo

      return
      end
      
c-----------------------------------------------------------------------
c Library subroutine mktens_zd
c
c Makes a tensor of third order by elementwise multiplying a real vector by
c a complex matrix:
c           T(i,j,k) = V(i)*M(j,k)   (no summation)
c           i=1,...,a;  j=1,...,b;  k=1,...,c
c
c Input-variables:  V     : complex vector
c                   M     : complex matrix
c                   a,b,c : dimensions of the different arrays
c                   add   : if true, result is added to T
c Output-variables: T     : resulting complex tensor of third order
c-----------------------------------------------------------------------

      subroutine mktens_zd (T,V,M,a,b,c,add)

      implicit none

      integer    i,j,k,a,b,c
      real*8     V(a)
      complex*16 T(a,b,c),M(b,c)
      logical    add

      if (.not. add) then
         do k=1,c
            do j=1,b
               do i=1,a
                  T(i,j,k)=(0.0d0,0.0d0)
               enddo
            enddo
         enddo
      endif
      do k=1,c
         do j=1,b
            do i=1,a
               T(i,j,k)=T(i,j,k)+V(i)*M(j,k)
            enddo
         enddo
      enddo

      return
      end

C***********************************************************************
C
C           Subroutine       NORMSPF
C
C   Normalizes a number of complex vectors (e.g. all single-particle 
C   functions of a specified degree of freedom) simultaneously
C
C   may/96  AJ
C
C***********************************************************************

      subroutine normspf(psi,subdim,dim)

      implicit none
      
      integer    e,g,subdim,dim
      real*8     norm
      complex*16 psi(subdim,dim)

      do e=1,dim
         norm=0d0
         do g=1,subdim
            norm=norm+dconjg(psi(g,e))*psi(g,e)
         enddo
         norm=1d0/sqrt(norm)
         do g=1,subdim
            psi(g,e)=norm*psi(g,e)
         enddo
      enddo
 
      return
      end

C***********************************************************************
C
C           Subroutine       NORMSPFR
C
C   Normalizes a number of real vectors simultaneously
C
C   may/96  AJ
C
C***********************************************************************

      subroutine normspfr(rvec,subdim,dim)

      implicit none
      
      integer    e,g,subdim,dim
      real*8     rvec(subdim,dim),norm

      do e=1,dim
         norm=0d0
         do g=1,subdim
            norm=norm+rvec(g,e)*rvec(g,e)
         enddo
         norm=1d0/sqrt(norm)
         do g=1,subdim
            rvec(g,e)=norm*rvec(g,e)
         enddo
      enddo
 
      return
      end


      subroutine vrminmax(v,n,vmin,vmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     finds the minimum and the maximum of a real*4 vector v
c
c     v (i)
c     n (i)
c     vmin, vmax (o)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer n,i,imax
      real*4 v(n),vmin,vmax,x

      imax=1
      vmin=v(1)
      vmax=vmin

      do i=2,n
         x=v(i)
         if (x.lt.vmin) then
            vmin=x
         elseif (x.gt.vmax) then
            imax=i
            vmax=x
         endif
      enddo

      return
      end


      subroutine vdminmax(v,n,vmin,vmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     finds the minimum and the maximum of a real*8 vector v
c
c     v (i)
c     n (i)
c     vmin, vmax (o)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer n,i
      real*8 v(n),vmin,vmax,x

      vmin=v(1)
      vmax=vmin

      do i=2,n
         x=v(i)
         if (x.lt.vmin) then
            vmin=x
         elseif (x.gt.vmax) then
            vmax=x
         endif
      enddo

      return
      end
      
