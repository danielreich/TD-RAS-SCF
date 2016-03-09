C **********************************************************************
C                                                                     
C                              MTLIB                                  
C                                                                     
C  Library module containing linear algebra routines that involve the 
C  multiplication of a tensor of third order with a matrix.           
C                                                                     
C  Nomenclature:
C    Each name has 6 basic characters:
C    First 2 characters denote the objects being multiplied:
C       q: quadratic matrix
C       m: general (rectangular) matrix
C       h: hermitian matrix
C       p: positive definite matrix
C       s: symmetric matrix
C       u: unitary matrix
C       d: diagonal matrix (only diagonal elements are supplied as a
C          vector)
C       t: 2nd index of tensor of third order
C       t1: 1st index of tensor of third order
C       v: vector
C       x: scalar
C       e.g. 'qt' denotes the operation (quadratic matrix * tensor 2nd
C            index)
C       e.g. 'qt1' denotes the operation (quadratic matrix * tensor 1st
C            index)
C    Character 3 denotes how first object is used:
C       x: unchanged from input
C       t: transpose of input
C       a: adjoint of input
C       c: complex conjugate of input
C    Character 4 denotes how second object is used:
C       see Character 3 above.
C    Character 5, 6 denote data types of first, second object
C       respectively:
C       z: complex double precision (complex*16)
C       c: complex single precision (complex*8)
C       d: real double precision (real*8)
C       r: real single precision (real*4)
C       y: complex matrix stored as two double precision (real*8)
C          matrices
C     Further characters, if present, give more informaion:
C       a: the result is added to a further object
C       r: the result is subtracted (removed) from a further object
C       c: the input matrices commute
C       h: the resulting matrix is hermitian
C       s: the resulting matrix is symmetric
C       1: the physical dimensions of the matrices differs from those
C          used.
C
C  The suffix _s after the name means that the routine works with a
C  "selected" vector", i.e. not all elements are present.
C   
C  Contents:                                                          
C  In the following list of available subroutines, matrices/tensors on
C  the LHS of the definition are input, that on the RHS output. The
C  usual summation convention is used i.e. a sum is made over repeated
C  indices on the LHS.
C
C   addtxxzz (a,b,dim1,dim2,dim3)
C       Definition: a(j) + b(i,j,k) = b(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3)
C
C   qtxxzz (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxzd (a,b,c,dim1,dim2,dim3)
C       Definition: dble(a(l,j))*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxzza (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) + c(i,l,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxzda (a,b,c,dim1,dim2,dim3)
C       Definition: dble(a(l,j))*b(i,j,k) + c(i,l,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qttxzza (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,l)*b(i,j,k) + c(i,l,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qttxzzs (a,b,c,dim1,dim2,dim3)
C       Definition: c(i,l,k) - a(j,l)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxzzr (a,b,c,dim1,dim2,dim3)
C       Definition: c(i,l,k) - a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   tqxxzza (a,b,c,dim1,dim2,dim3)
C       Definition: a(i,j,k)*b(j,l) + c(i,l,k) = c(i,l,k) .
C       Dimensions: a(dim1,dim2,dim3),b(dim2,dim2),c(dim1,dim2,dim3)
C
C   tqxazza (a,b,c,dim1,dim2,dim3)
C       Definition: a(i,j,k)*dconjg(b(l,j)) + c(i,l,k) = c(i,l,k) .
C       Dimensions: a(dim1,dim2,dim3),b(dim2,dim2),c(dim1,dim2,dim3)
C
C   tqxazz (a,b,c,dim1,dim2,dim3)
C       Definition: a(i,j,k)*dconjg(b(l,j)) = c(i,l,k) .
C       Dimensions: a(dim1,dim2,dim3),b(dim2,dim2),c(dim1,dim2,dim3)
C
C   mtxxzz (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mtcxzz (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: conjg(a(l,j))*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mtxxdz (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mttxzz (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: a(j,l)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mttxdd (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: a(j,l)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mttxdr (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: a(j,l)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mtxxdd (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
C   mtaxzz (a,b,c,dim1,dim2,dim3,dim4)
C       Definition: dconjg(a(j,l))*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)
C
CC  qtxxdd (a,b,c,dim1,dim2,dim3)
CC      Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
CC      Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtcxzz (a,b,c,dim1,dim2,dim3)
C       Definition: dconjg(a(l,j))*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qttxdd (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   dtxxzz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j)*b(i,j,k) = c(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   dtxxzzo (a,b,dim1,dim2,dim3)
C       Definition: a(j)*b(i,j,k) = b(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3)
C
C   dtxxdd (a,b,c,dim1,dim2,dim3)
C       Definition: a(j)*b(i,j,k) = c(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   dtxxddo (a,b,dim1,dim2,dim3)
C       Definition: a(j)*b(i,j,k) = b(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3)
C
C   dtxxdz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j)*b(i,j,k) = c(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   dtxxdzo (a,b,dim1,dim2,dim3)
C       Definition: a(j)*b(i,j,k) = b(i,j,k) .
C       Dimensions: a(dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtaxzz (a,b,c,dim1,dim2,dim3)
C       Definition: dconjg(a(j,l))*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtatzz (a,b,c,dim1,dim2,dim3)
C       Definition: dconjg(a(j,l))*b(i,k,j) = c(i,l,k) .
C       Dimension:  a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxdz (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxdza (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) + c(i,l,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qt1xxdz (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,i)*b(i,j,k) = c(l,j,k) .
C       Dimensions: a(dim1,dim1),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   dt1xxdz (a,b,c,dim1,dim2,dim3)
C       Definition: a(i)*b(i,j,k) = c(i,j,k) .
C       Dimensions: a(dim1),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qt1txdz (a,b,c,dim1,dim2,dim3)
C       Definition: a(i,l)*b(i,j,k) = c(l,j,k) .
C       Dimensions: a(dim1,dim1),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qttxdz (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qttxzz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,l)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxyz (a,b,c,dim1,dim2,dim3)
C       Definition: a(l,j)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qttxyz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,l)*b(i,j,k) = c(i,l,k) .
C       Dimensions: a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   dtxxyz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,1)*b(i,j,k) + (0,1)*a(j,2)*b(i,j,k) = c(i,j,k) .
C       Dimensions: a(dim2,2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)
C
C   qtxxzz_s (a,b,c,dim1,dim2,dim3,index,index1)
C       Definition: index(i,j,k)=x
C                   index1(i,l,k)=x1
C                   a(l,j)*b(x) = c(x1) .
C       Dimensions: a(dim2,dim2),index(dim1,dim2,dim3),
C                   index1(dim1,dim2,dim3),b(max(x)),c(max(x1))
C
C   mtxxzz_s (a,b,c,dim1,dim2,dim3,dim4,dim5,index,index1)
C       Definition: index(i,j,k)=x
C                   index1(i,l,k)=x1
C                   a(l,j)*b(x) = c(x1) .
C       Dimensions: a(dim4,dim5),index(dim1,dim2,dim3),
C                   index1(dim1,dim2,dim3),b(max(x)),c(max(x1))
C
C   qttxzz_s (a,b,c,dim1,dim2,dim3,index,index1)
C       Definition: index(i,j,k)=x
C                   index1(i,l,k)=x1
C                   a(j,l)*b(x) = c(x1) .
C       Dimensions: a(dim2,dim2),index(dim1,dim2,dim3),
C                   index1(dim1,dim2,dim3),b(max(x)),c(max(x1))
C
C   qtxxzza_s (a,b,c,dim1,dim2,dim3,index,index1)
C       Definition: index(i,j,k)=x
C                   index1(i,l,k)=x1
C                   a(l,j)*b(x) + c(x1) = c(x1) .
C       Dimensions: a(dim2,dim2),index(dim1,dim2,dim3),
C                   index1(dim1,dim2,dim3),b(max(x)),c(max(x1))
C
C   qtxxzzr_s (a,b,c,dim1,dim2,dim3,index,index1)
C       Definition: index(i,j,k)=x
C                   index1(i,l,k)=x1
C                   c(x1) - a(l,j)*b(x) = c(x1) .
C       Dimensions: a(dim2,dim2),index(dim1,dim2,dim3),
C                   index1(dim1,dim2,dim3),b(max(x)),c(max(x1))
C
C **********************************************************************

c-----------------------------------------------------------------------
c Library subroutine addtxxzz
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order:
c           a(j) + b(i,j,k) =  b(i,j,k).
c
c Input-variables:  a      - complex vector 
c                   b      - complex tensor of third order
c Output-variables: b      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine addtxxzz (a,b,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      complex*16    a(dim2),b(dim1,dim2,dim3)

      do k=1,dim3
         do j=1,dim2
            do i=1,dim1
               b(i,j,k)=b(i,j,k)+a(j)
            enddo
         enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
c Library subroutine qtxxzz
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order:
c           a(l,j)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtxxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=a(l,1)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=a(l,1)*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=a(l,1)*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qtxxzd
c
c Multiplication of a complex quadratic matrix with a real tensor
c of third order:
c           dble(a(l,j))*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix (only real part used)
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine qtxxzd (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        b(dim1,dim2,dim3),c(dim1,dim2,dim3)
      complex*16    a(dim2,dim2)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=dble(a(l,1))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+dble(a(l,j))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=dble(a(l,1))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+dble(a(l,j))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=dble(a(l,1))*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dble(a(l,j))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=dble(a(l,1))*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+dble(a(l,j))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end


c-----------------------------------------------------------------------
c Library subroutine qtxxzza
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and added to a different tensor:
c           a(l,j)*b(i,j,k) + c(i,l,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c                   c      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
C
C NB input c matrix is overwritten on output.
c-----------------------------------------------------------------------

      subroutine qtxxzza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qtxxzda
c
c Multiplication of a complex quadratic matrix with a complex tensor
c of third order and added to a different tensor:
c           dble(a(l,j))*b(i,j,k) + c(i,l,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix
c                   b      - real tensor of third order
c                   c      - real tensor of third order
c Output-variables: c      - resulting complex tensor
C
C NB input c matrix is overwritten on output.
c-----------------------------------------------------------------------

      subroutine qtxxzda (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        b(dim1,dim2,dim3),c(dim1,dim2,dim3)
      complex*16    a(dim2,dim2)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+dble(a(l,j))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+dble(a(l,j))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dble(a(l,j))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+dble(a(l,j))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qttxzza
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and added to a different tensor:
c           a(j,l)*b(i,j,k) + c(i,l,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c                   c      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
C
C NB input c matrix is overwritten on output.
c-----------------------------------------------------------------------

      subroutine qttxzza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 


c-----------------------------------------------------------------------
c Library subroutine qttxzzs
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and substracted from a different tensor:
c            c(i,l,k) - a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c                   c      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
C
C NB input c matrix is overwritten on output.
c-----------------------------------------------------------------------

      subroutine qttxzzs (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)-a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)-a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)-a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)-a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 


c-----------------------------------------------------------------------
c Library subroutine qtxxzzr
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and subtracted from a different tensor:
c           c(i,l,k) - a(l,j)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c                   c      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtxxzzr (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)-a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)-a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)-a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)-a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 


c-----------------------------------------------------------------------
c Library subroutine tqxxzza
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and added to a different tensor:
c           a(i,j,k)*b(j,l) + c(i,l,k) = c(i,l,k).
c
c Input-variables:  a      - complex tensor of third order 
c                   b      - complex matrix
c                   c      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
C
C NB input c matrix is overwritten on output.
c-----------------------------------------------------------------------

      subroutine tqxxzza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim1,dim2,dim3),b(dim2,dim2),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(i,j,k)*b(j,l)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  c(1,l,k)=c(1,l,k)+a(1,j,k)*b(j,l)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(i,j,1)*b(j,l)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               c(1,l,1)=c(1,l,1)+a(1,j,1)*b(j,l)
            enddo
         enddo
      endif

      return
      end ! of tqxxzza



c-----------------------------------------------------------------------
c Library subroutine tqxazza
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and added to a different tensor:
c           a(i,j,k)*dconjg(b(l,j)) + c(i,l,k) = c(i,l,k).
c
c Input-variables:  a      - complex tensor of third order 
c                   b      - complex matrix
c                   c      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
C
C NB input c matrix is overwritten on output.
c-----------------------------------------------------------------------

      subroutine tqxazza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim1,dim2,dim3),b(dim2,dim2),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(i,j,k)*dconjg(b(l,j))
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  c(1,l,k)=c(1,l,k)+a(1,j,k)*dconjg(b(l,j))
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(i,j,1)*dconjg(b(l,j))
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               c(1,l,1)=c(1,l,1)+a(1,j,1)*dconjg(b(l,j))
            enddo
         enddo
      endif

      return
      end ! of tqxazza



c-----------------------------------------------------------------------
c Library subroutine tqxazz
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order:
c           a(i,j,k)*dconjg(b(l,j)) = c(i,l,k).
c
c Input-variables:  a      - complex tensor of third order 
c                   b      - complex matrix
c
c Output-variables: c      - resulting complex tensor
C
c-----------------------------------------------------------------------

      subroutine tqxazz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim1,dim2,dim3),b(dim2,dim2),c(dim1,dim2,dim3)

       If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=a(i,1,k)*dconjg(b(l,1))
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(i,j,k)*dconjg(b(l,j))
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=a(1,1,k)*dconjg(b(l,1))
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(1,j,k)*dconjg(b(l,j))
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=a(i,1,1)*dconjg(b(l,1))
            enddo
         enddo
         do l=1,dim2
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(i,j,1)*dconjg(b(l,j))
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=a(1,1,1)*dconjg(b(l,1))
         enddo
         do l=1,dim2
           do j=2,dim2
              c(1,l,1)=c(1,l,1)+a(1,j,1)*dconjg(b(l,j))
            enddo
         enddo
      endif

      return
      end ! of tqxazz


c-----------------------------------------------------------------------
c Library subroutine mtxxzz
c
c Multiplication of a complex rectangular matrix with a complex tensor 
c of third order:
c           a(l,j)*b(i,j,k) = c(i,l,k) .
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine mtxxzz (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      complex*16    a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=a(l,1)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=a(l,1)*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim4
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=a(l,1)*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim4
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine mtcxzz
c
c Multiplication of a complex rectangular matrix with a complex tensor 
c of third order:
c           conjg(a(l,j))*b(i,j,k) = c(i,l,k) .
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine mtcxzz (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      complex*16    a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=conjg(a(l,1))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+conjg(a(l,j))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=conjg(a(l,1))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  c(1,l,k)=c(1,l,k)+conjg(a(l,j))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=conjg(a(l,1))*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim4
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+conjg(a(l,j))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=conjg(a(l,1))*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim4
               c(1,l,1)=c(1,l,1)+conjg(a(l,j))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 


c-----------------------------------------------------------------------
c Library subroutine mtxxdz
c
c Multiplication of a real rectangular matrix with a complex tensor 
c of third order:
c           a(l,j)*b(i,j,k) = c(i,l,k) .
c
c Input-variables:  a      - real rectangular matrix
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine mtxxdz (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      real*8        a(dim4,dim2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=a(l,1)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=a(l,1)*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim4
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=a(l,1)*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim4
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 


c-----------------------------------------------------------------------
c Library subroutine mttxzz
c
c Multiplication of a transposed complex rectangular matrix with a 
c complex tensor of third order:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine mttxzz (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      complex*16    a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=a(1,l)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=a(1,l)*b(i,1,1)
            enddo
         enddo
         do l=1,dim4
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim4
            do j=2,dim2
               c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine mttxdd
c
c Multiplication of a transposed real rectangular matrix with a 
c real tensor of third order:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine mttxdd (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      real*8        a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=a(1,l)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=a(1,l)*b(i,1,1)
            enddo
         enddo
         do l=1,dim4
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim4
            do j=2,dim2
               c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine mttxdr
c
c Multiplication of a transposed real rectangular matrix with a 
c real tensor of third order:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine mttxdr (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      real*8        a(dim2,dim4),c(dim1,dim4,dim3)
      real*4        b(dim1,dim2,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=a(1,l)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=a(1,l)*b(i,1,1)
            enddo
         enddo
         do l=1,dim4
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim4
            do j=2,dim2
               c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine mtxxdd
c
c Multiplication of a real rectangular matrix with a complex tensor 
c of third order:
c           a(l,j)*b(i,j,k) = c(i,l,k) .
c
c Input-variables:  a      - real matrix 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine mtxxdd (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      real*8        a(dim4,dim2),b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=a(l,1)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim4
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
                c(i,l,1)=a(l,1)*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim4
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
        do l=1,dim4
            c(1,l,1)=a(l,1)*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim4
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine mtaxzz
c
c Multiplication of the adjoint of a complex rectangular matrix with a 
c complex tensor of third order:
c           dconj(a(j,l))*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine mtaxzz (a,b,c,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      complex*16    a(dim2,dim4),b(dim1,dim2,dim3),c(dim1,dim4,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do i=1,dim1
                  c(i,l,k)=dconjg(a(1,l))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+dconjg(a(j,l))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               c(1,l,k)=dconjg(a(1,l))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim4
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+dconjg(a(j,l))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do i=1,dim1
               c(i,l,1)=dconjg(a(1,l))*b(i,1,1)
            enddo
         enddo
         do l=1,dim4
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dconjg(a(j,l))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            c(1,l,1)=dconjg(a(1,l))*b(1,1,1)
         enddo
         do l=1,dim4
            do j=2,dim2
               c(1,l,1)=c(1,l,1)+dconjg(a(j,l))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qtxxdd
c
c Multiplication of a real quadratic matrix with a real tensor 
c of third order:
c           a(l,j)*b(i,j,k) = c(i,l,k) .
c
c Input-variables:  a      - real matrix 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

C     subroutine qtxxdd (a,b,c,dim1,dim2,dim3)

C     implicit none

C     integer       dim1,dim2,dim3, i, j, k, l
C     real*8        a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

C     if (dim1.ne.1 .and. dim3.ne.1) then
C        do k=1,dim3
C           do l=1,dim2
C              do i=1,dim1
C                 c(i,l,k)=a(l,1)*b(i,1,k)
C              enddo
C           enddo
C        enddo
C        do k=1,dim3
C           do j=2,dim2
C              do l=1,dim2
C                 do i=1,dim1
C                    c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
C                 enddo
C              enddo
C           enddo
C        enddo
C     else if (dim1.eq.1 .and. dim3.ne.1) then
C        do k=1,dim3
C           do l=1,dim2
C              c(1,l,k)=a(l,1)*b(1,1,k)
C           enddo
C        enddo
C        do k=1,dim3
C           do j=2,dim2
C              do l=1,dim2
C                 c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
C              enddo
C           enddo
C        enddo
C     else if (dim1.ne.1 .and. dim3.eq.1) then
C        do l=1,dim2
C           do i=1,dim1
C              c(i,l,1)=a(l,1)*b(i,1,1)
C           enddo
C        enddo
C        do j=2,dim2
C           do l=1,dim2
C              do i=1,dim1
C                 c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
C              enddo
C           enddo
C        enddo
C     else if (dim1.eq.1 .and. dim3.eq.1) then
C        do l=1,dim2
C           c(1,l,1)=a(l,1)*b(1,1,1)
C        enddo
C        do j=2,dim2
C           do l=1,dim2
C              c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
C           enddo
C        enddo
C     endif

C     return
C     end 

c-----------------------------------------------------------------------
c Library subroutine qtcxzz
c
c Multiplication of the complex conjugate of a complex quadratic matrix 
c with a complex tensor of third order:
c           dconjg(a(l,j))*b(i,j,k) = c(i,l,k) .
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtcxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=dconjg(a(l,1))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+dconjg(a(l,j))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=dconjg(a(l,1))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+dconjg(a(l,j))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=dconjg(a(l,1))*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dconjg(a(l,j))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=dconjg(a(l,1))*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+dconjg(a(l,j))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qttxdd
c
c Multiplication of a transposed real quadratic matrix with a 
c real tensor of third order:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine qttxdd (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=a(1,l)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=a(1,l)*b(i,1,1)
            enddo
         enddo
         do l=1,dim2
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim2
           do j=2,dim2
              c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine dtxxzz
c
c Multiplication of a complex diagonal matrix with a complex tensor 
c of third order. Only diagonal matrix elements are given as a vector:
c           a(j)*b(i,j,k) = c(i,j,k) .
c
c Input-variables:  a      - complex vector with diagonal matrix elements 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine dtxxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      complex*16    a(dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  c(i,j,k)=a(j)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               c(1,j,k)=a(j)*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               c(i,j,1)=a(j)*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            c(1,j,1)=a(j)*b(1,j,1)
         enddo
      endif

      return
      end  

c-----------------------------------------------------------------------
c Library subroutine dtxxzzo
c
c Multiplication of a complex diagonal matrix with a complex tensor
c of third order. Only diagonal matrix elements are given as a vector,
c and the result is sored in the input tensor:
c           a(j)*b(i,j,k) = b(i,j,k) .
c
c Input-variables:  a      - complex vector with diagonal matrix elements
c                   b      - complex tensor of third order
c Output-variables: b      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine dtxxzzo (a,b,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      complex*16    a(dim2),b(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  b(i,j,k)=a(j)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               b(1,j,k)=a(j)*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               b(i,j,1)=a(j)*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            b(1,j,1)=a(j)*b(1,j,1)
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine dtxxdd
c
c Multiplication of a real diagonal matrix with a real tensor 
c of third order. Only diagonal matrix elements are given as a vector:
c           a(j)*b(i,j,k) = c(i,j,k) .
c
c Input-variables:  a      - real vector with diagonal matrix elements 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine dtxxdd (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      real*8        a(dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  c(i,j,k)=a(j)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               c(1,j,k)=a(j)*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               c(i,j,1)=a(j)*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            c(1,j,1)=a(j)*b(1,j,1)
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine dtxxddo
c
c Multiplication of a real diagonal matrix with a real tensor
c of third order. Only diagonal matrix elements are given as a vector,
c and the result is sored in the input tensor:
c           a(j)*b(i,j,k) = b(i,j,k) .
c
c Input-variables:  a      - real vector with diagonal matrix elements
c                   b      - real tensor of third order
c Output-variables: b      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine dtxxddo (a,b,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      real*8        a(dim2),b(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  b(i,j,k)=a(j)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               b(1,j,k)=a(j)*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               b(i,j,1)=a(j)*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            b(1,j,1)=a(j)*b(1,j,1)
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine dtxxdz
c
c Multiplication of a real diagonal matrix with a complex tensor 
c of third order. Only diagonal matrix elements are given as a vector:
c           a(j)*b(i,j,k) = c(i,j,k) .
c
c Input-variables:  a      - real vector with diagonal matrix elements 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine dtxxdz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      real*8        a(dim2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  c(i,j,k)=a(j)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               c(1,j,k)=a(j)*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               c(i,j,1)=a(j)*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            c(1,j,1)=a(j)*b(1,j,1)
         enddo
      endif

      return
      end  

c-----------------------------------------------------------------------
c Library subroutine dtxxdzo
c
c Multiplication of a real diagonal matrix with a complex tensor
c of third order. Only diagonal matrix elements are given as a vector,
c and the result is sored in the input tensor:
c           a(j)*b(i,j,k) = b(i,j,k) .
c
c Input-variables:  a      - complex vector with diagonal matrix elements
c                   b      - complex tensor of third order
c Output-variables: b      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine dtxxdzo (a,b,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      real*8        a(dim2)
      complex*16    b(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  b(i,j,k)=a(j)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               b(1,j,k)=a(j)*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               b(i,j,1)=a(j)*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            b(1,j,1)=a(j)*b(1,j,1)
         enddo
      endif

      return
      end


c-----------------------------------------------------------------------
c Library subroutine qtaxzz
c
c Multiplication of the adjoint of a complex quadratic matrix with a 
c complex tensor of third order:
c           dconj(a(j,l))*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtaxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=dconjg(a(1,l))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+dconjg(a(j,l))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=dconjg(a(1,l))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+dconjg(a(j,l))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=dconjg(a(1,l))*b(i,1,1)
            enddo
         enddo
         do l=1,dim2
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dconjg(a(j,l))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=dconjg(a(1,l))*b(1,1,1)
         enddo
         do l=1,dim2
            do j=2,dim2
               c(1,l,1)=c(1,l,1)+dconjg(a(j,l))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qtxxdz
c
c Multiplication of a real quadratic matrix with a complex tensor 
c of third order:
c           a(l,j)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtxxdz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim2,dim2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=a(l,1)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=a(l,1)*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=a(l,1)*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qtxxdza
c
c Multiplication of a real quadratic matrix with a complex tensor 
c of third order:
c           a(l,j)*b(i,j,k) + c(i,l,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtxxdza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim2,dim2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(l,j)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+a(l,j)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(l,j)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+a(l,j)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qt1xxdz
c
c Multiplication of a real quadratic matrix with the first index of a 
c  complex tensor of third order:
c           a(l,i)*b(i,j,k) = c(l,j,k) .
c
c Input-variables:  a      - real matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qt1xxdz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim1,dim1)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim2.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim1
                  c(l,j,k)=a(l,1)*b(1,j,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=1,dim2
               do i=2,dim1
                  do l=1,dim1
                     c(l,j,k)=c(l,j,k)+a(l,i)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim2.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim1
               c(l,1,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do i=2,dim1
               do l=1,dim1
                  c(l,1,k)=c(l,1,k)+a(l,i)*b(i,1,k)
               enddo
            enddo
         enddo
      else if (dim2.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim1
               c(l,j,1)=a(l,1)*b(1,j,1)
            enddo
         enddo
         do j=1,dim2
            do i=2,dim1
               do l=1,dim1
                  c(l,j,1)=c(l,j,1)+a(l,i)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim1
            c(l,1,1)=a(l,1)*b(1,1,1)
         enddo
         do i=2,dim1
            do l=1,dim1
               c(l,1,1)=c(l,1,1)+a(l,i)*b(i,1,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine dt1xxzz
c
c Multiplication of a complex diagonal matrix with the first index of a 
c  complex tensor of third order. Only diagonal matrix elements are
c  supplied as  vector:
c           a(i)*b(i,j,k) = c(i,j,k) .
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine dt1xxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      complex*16    a(dim1),b(dim1,dim2,dim3),c(dim1,dim2,dim3)


      if (dim2.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do i=1,dim1
                  c(i,j,k)=a(i)*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim2.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do i=1,dim1
               c(i,1,k)=a(i)*b(i,1,k)
            enddo
         enddo
      else if (dim2.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do i=1,dim1
               c(i,j,1)=a(i)*b(i,j,1)
            enddo
         enddo
      else if (dim2.eq.1 .and. dim3.eq.1) then
         do i=1,dim1
            c(i,1,1)=a(i)*b(i,1,1)
         enddo
      endif

      return
      end  

c-----------------------------------------------------------------------
c Library subroutine qt1txdz
c
c Multiplication of a transposed real quadratic matrix with the 
c  first index of a complex tensor of third order:
c           a(i,l)*b(i,j,k) = c(l,j,k) .
c
c Input-variables:  a      - real matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qt1txdz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim1,dim1)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim2.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim1
                  c(l,j,k)=a(1,l)*b(1,j,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=1,dim2
               do l=1,dim1
                  do i=2,dim1
                     c(l,j,k)=c(l,j,k)+a(i,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim2.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim1
               c(l,1,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim1
               do i=2,dim1
                  c(l,1,k)=c(l,1,k)+a(i,l)*b(i,1,k)
               enddo
            enddo
         enddo
      else if (dim2.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim1
               c(l,j,1)=a(1,l)*b(1,j,1)
            enddo
         enddo
         do j=1,dim2
            do l=1,dim1
               do i=2,dim1
                  c(l,j,1)=c(l,j,1)+a(i,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim1
            c(l,1,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim1
            do i=2,dim1
               c(l,1,1)=c(l,1,1)+a(i,l)*b(i,1,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qttxdz
c
c Multiplication of a transposed real quadratic matrix with a 
c complex tensor of third order:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qttxdz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim2,dim2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=a(1,l)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=a(1,l)*b(i,1,1)
            enddo
         enddo
         do l=1,dim2
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim2
           do j=2,dim2
              c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qt1xxzz
c
c Multiplication of a real quadratic matrix with the first index of a 
c  complex tensor of third order:
c           a(l,i)*b(i,j,k) = c(l,j,k) .
c
c Input-variables:  a      - real matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qt1xxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim1,dim1),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim2.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim1
                  c(l,j,k)=a(l,1)*b(1,j,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=1,dim2
               do i=2,dim1
                  do l=1,dim1
                     c(l,j,k)=c(l,j,k)+a(l,i)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim2.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim1
               c(l,1,k)=a(l,1)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do i=2,dim1
               do l=1,dim1
                  c(l,1,k)=c(l,1,k)+a(l,i)*b(i,1,k)
               enddo
            enddo
         enddo
      else if (dim2.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim1
               c(l,j,1)=a(l,1)*b(1,j,1)
            enddo
         enddo
         do j=1,dim2
            do i=2,dim1
               do l=1,dim1
                  c(l,j,1)=c(l,j,1)+a(l,i)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim1
            c(l,1,1)=a(l,1)*b(1,1,1)
         enddo
         do i=2,dim1
            do l=1,dim1
               c(l,1,1)=c(l,1,1)+a(l,i)*b(i,1,1)
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qttxzz
c
c Multiplication of a transposed complex quadratic matrix with a 
c complex tensor of third order:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - real matrix 
c                   b      - real tensor of third order
c Output-variables: c      - resulting real tensor
c-----------------------------------------------------------------------

      subroutine qttxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    a(dim2,dim2),b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      If (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=a(1,l)*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  do i=1,dim1
                     c(i,l,k)=c(i,l,k)+a(j,l)*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=a(1,l)*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do l=1,dim2
               do j=2,dim2
                  c(1,l,k)=c(1,l,k)+a(j,l)*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=a(1,l)*b(i,1,1)
            enddo
         enddo
         do l=1,dim2
            do j=2,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+a(j,l)*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=a(1,l)*b(1,1,1)
         enddo
         do l=1,dim2
           do j=2,dim2
              c(1,l,1)=c(1,l,1)+a(j,l)*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qtxxyz
c
c Multiplication of a complex quadratic matrix with a complex tensor
c of third order, where the complex matrix is stored as two real
c matrices:
c           a(l,j)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtxxyz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim2,dim2,2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=dcmplx(a(l,1,1),a(l,1,2))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  do i=1,dim1
                    c(i,l,k)=c(i,l,k)+dcmplx(a(l,j,1),a(l,j,2))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=dcmplx(a(l,1,1),a(l,1,2))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+dcmplx(a(l,j,1),a(l,j,2))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=dcmplx(a(l,1,1),a(l,1,2))*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dcmplx(a(l,j,1),a(l,j,2))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=dcmplx(a(l,1,1),a(l,1,2))*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+dcmplx(a(l,j,1),a(l,j,2))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qttxyz
c
c Multiplication of the trabspose of a complex quadratic matrix with 
c a complex tensor of third order, where the complex matrix is stored 
c as two real matrices:
c           a(j,l)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qttxyz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        a(dim2,dim2,2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do i=1,dim1
                  c(i,l,k)=dcmplx(a(1,l,1),a(1,l,2))*b(i,1,k)
               enddo
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  do i=1,dim1
                    c(i,l,k)=c(i,l,k)+dcmplx(a(j,l,1),a(j,l,2))*b(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               c(1,l,k)=dcmplx(a(1,l,1),a(1,l,2))*b(1,1,k)
            enddo
         enddo
         do k=1,dim3
            do j=2,dim2
               do l=1,dim2
                  c(1,l,k)=c(1,l,k)+dcmplx(a(j,l,1),a(j,l,2))*b(1,j,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do i=1,dim1
               c(i,l,1)=dcmplx(a(1,l,1),a(1,l,2))*b(i,1,1)
            enddo
         enddo
         do j=2,dim2
            do l=1,dim2
               do i=1,dim1
                  c(i,l,1)=c(i,l,1)+dcmplx(a(j,l,1),a(j,l,2))*b(i,j,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            c(1,l,1)=dcmplx(a(1,l,1),a(1,l,2))*b(1,1,1)
         enddo
         do j=2,dim2
            do l=1,dim2
               c(1,l,1)=c(1,l,1)+dcmplx(a(j,l,1),a(j,l,2))*b(1,j,1)
            enddo
         enddo
      endif

      return
      end


c-----------------------------------------------------------------------
c Library subroutine dtxxyz
c
c Multiplication of a complex quadratic diagonal matrix with a complex
c tensor of third order, where the complex matrix is stored as two real
c matrices:
c           a(j)*b(i,j,k) = c(i,l,k).
c
c Input-variables:  a      - complex matrix
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine dtxxyz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k
      real*8        a(dim2,2)
      complex*16    b(dim1,dim2,dim3),c(dim1,dim2,dim3)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=2,dim2
               do i=1,dim1
                  c(i,j,k)=dcmplx(a(j,1),a(j,2))*b(i,j,k)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=2,dim2
               c(1,j,k)=dcmplx(a(j,1),a(j,2))*b(1,j,k)
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=2,dim2
            do i=1,dim1
               c(i,j,1)=dcmplx(a(j,1),a(j,2))*b(i,j,1)
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=2,dim2
            c(1,j,1)=dcmplx(a(j,1),a(j,2))*b(1,j,1)
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine qtxxzz_s
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order. The tensor is not completely stored, and the indices
c are managed by the index tensors:
c           index(i,j,k)=x
c           index1(i,l,k)=x1
c           a(l,j)*b(x) = c(x1).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qtxxzz_s (a,b,c,dim1,dim2,dim3,index,index1)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l,index(dim1,dim2,dim3),
     +              index1(dim1,dim2,dim3),x,x1
      complex*16    a(dim2,dim2),b(*),c(*)

      do k=1,dim3
         do l=1,dim2
            do i=1,dim1
               x1=index1(i,l,k)
               if (x1 .ne. 0) c(x1)=0.0d0
            enddo
         enddo
      enddo

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        c(x1)=c(x1)+a(l,j)*b(x)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  c(x1)=c(x1)+a(l,j)*b(x)
               endif
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine mtxxzz_s
c
c Multiplication of a complex matrix with a complex tensor 
c of third order. The tensor is not completely stored, and the indices
c are managed by the index tensors:
c           index(i,j,k)=x
c           index1(i,l,k)=x1
c           a(l,j)*b(x) = c(x1).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine mtxxzz_s (a,b,c,dim1,dim2,dim3,dim4,dim5,index,index1)

      implicit none

      integer       dim1,dim2,dim3,dim4,dim5,i,j,k,l,
     +              index(dim1,dim2,dim3),index1(dim1,dim2,dim3),x,x1
      complex*16    a(dim5,dim4),b(*),c(*)

      do k=1,dim3
         do l=1,dim2
            do i=1,dim1
               x1=index1(i,l,k)
               if (x1 .ne. 0) c(x1)=0.0d0
            enddo
         enddo
      enddo
 
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim4
               do l=1,dim5
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        c(x1)=c(x1)+a(l,j)*b(x)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim4
               do l=1,dim5
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim4
            do l=1,dim5
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim4
            do l=1,dim5
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  c(x1)=c(x1)+a(l,j)*b(x)
               endif
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qttxzz_s
c
c Multiplication of a transposed complex quadratic matrix with a 
c complex tensor c of third order. The tensor is not completely stored, 
c and the indices are managed by the index tensors:
c           index(i,j,k)=x
c           index1(i,l,k)=x1
c           a(j,l)*b(x) = c(x1).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c-----------------------------------------------------------------------

      subroutine qttxzz_s (a,b,c,dim1,dim2,dim3,index,index1)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l,index(dim1,dim2,dim3),
     +              index1(dim1,dim2,dim3),x,x1
      complex*16    a(dim2,dim2),b(*),c(*)

      do k=1,dim3
         do l=1,dim2
            do i=1,dim1
               x1=index1(i,l,k)
               if (x1 .ne. 0) c(x1)=0.0d0
            enddo
         enddo
      enddo

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        c(x1)=c(x1)+a(j,l)*b(x)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(j,l)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(j,l)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  c(x1)=c(x1)+a(j,l)*b(x)
               endif
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qtxxzza_s
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and added to a different tensor. The tensors are not 
c completely stored, and the indices are managed by the index tensors:
c           index(i,j,k)=x
c           index1(i,l,k)=x1
c           a(l,j)*b(x) + c(x1) = c(x1).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c
c NB tensor c is overwritten on output
c-----------------------------------------------------------------------

      subroutine qtxxzza_s (a,b,c,dim1,dim2,dim3,index,index1)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l,index(dim1,dim2,dim3),
     +              index1(dim1,dim2,dim3),x,x1
      complex*16    a(dim2,dim2),b(*),c(*)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        c(x1)=c(x1)+a(l,j)*b(x)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)+a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  c(x1)=c(x1)+a(l,j)*b(x)
               endif
            enddo
         enddo
      endif

      return
      end 

c-----------------------------------------------------------------------
c Library subroutine qtxxzzr_s
c
c Multiplication of a complex quadratic matrix with a complex tensor 
c of third order and subtracted from a different tensor. The tensors are 
c not completely stored, and the indices are managed by the index tensors:
c           index(i,j,k)=x
c           index1(i,l,k)=x1
c           c(x1) - a(l,j)*b(x) = c(x1).
c
c Input-variables:  a      - complex matrix 
c                   b      - complex tensor of third order
c Output-variables: c      - resulting complex tensor
c
c NB tensor c is overwritten on output
c-----------------------------------------------------------------------

      subroutine qtxxzzr_s (a,b,c,dim1,dim2,dim3,index,index1)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l,index(dim1,dim2,dim3),
     +              index1(dim1,dim2,dim3),x,x1
      complex*16    a(dim2,dim2),b(*),c(*)

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        c(x1)=c(x1)-a(l,j)*b(x)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim2
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)-a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     c(x1)=c(x1)-a(l,j)*b(x)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim2
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  c(x1)=c(x1)-a(l,j)*b(x)
               endif
            enddo
         enddo
      endif

      return
      end 


      

