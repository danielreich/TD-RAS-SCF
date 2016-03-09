C **********************************************************************
C                                                                     
C                              RMLIB                                  
C                                                                     
C  Library module containing linear algebra routines that involve the 
C  formation of a "density" matrix
C              
C  Nomenclature:
C    Each name is 8 or 9 characters long.
C    First 2 Characters:
C       'rm' denote the operation 'formation of "density" matrix'
C       'rm1' denote the operation 'formation of "density" matrix',
C             where a tensor is multiplied with itself. 
C    Character 3 denotes type of matrix formed:
C       q: quadratic
C       m: general (rectangular) matrix
C       h: hermitian
C       a: anti-hermitian
C       s: symmetric
C       u: unitary
C       d: diagonal (only diagonal elements are saved as a vector)
C    Character 4 denotes how matrix is saved:
C       x: unchanged from calculated
C       t: transpose of calculated
C       a: adjoint of calculated
C       c: complex conjugate of calculated
C    Character 5 denotes how 'bra' tensor is used: 
C       x: complex conjugate of input (standard definition)
C       c: unchanged from input
C    Character 6 denotes how 'ket' tensor is used:
C       x: unchanged from input
C       c: complex conjugate of input
C    Character 7 and 8 denote data types of tensors:
C       z: complex double precision (complex*16)
C       c: complex single precision (complex*8)
C       d: real double precision (real*8)
C       r: real single precision (real*4)
C    Character 9 (optional) denotes data type of matrix:
C       z: complex double precision (complex*16)
C       c: complex single precision (complex*8)
C       d: real double precision (real*8)
C       r: real single precision (real*4)
C   
C  The suffix _s after the name means that the routine works with a
C  "selected" vector", i.e. not all elements are present.
C
C  Contents:                                                          
C  In the following list of available subroutines, matrices/tensors on the 
C  LHS of the definition are input, that on the RHS output. The usual summation
C  convention is used i.e. a sum is made over repeated indices on the LHS.
C
C    rmqxxxzz (bra,ket,mat,dim1,dim2,dim3)
C        Definition: dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
C    rmqxcxzz (bra,ket,mat,dim1,dim2,dim3)
C        Definition: bra(i,j,k)*ket(i,l,k) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
C    rmmxxxzz (bra,ket,mat,dim1,dim2,dim3,dim4)
C        Definition: dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l) .
C        Dimensions: mat(dim2,dim4),bra(dim1,dim2,dim3),ket(dim1,dim4,dim3)
C
C    rmmtxxzz (bra,ket,mat,dim1,dim2,dim3,dim4)
C        Definition: dconjg(bra(i,j,k))*ket(i,l,k) = mat(l,j) .
C        Dimensions: mat(dim4,dim2),bra(dim1,dim2,dim3),ket(dim1,dim4,dim3)
C
C    rmhxxxzz (bra,ket,mat,dim1,dim2,dim3)
C        Definition: dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
C    rmaxxxzz (bra,ket,mat,dim1,dim2,dim3)
C        Definition: dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
C    rmsxcxdd (bra,ket,mat,dim1,dim2,dim3)
C        Definition: bra(i,j,k)*ket(i,l,k) = mat(j,l) .
C     Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
C    rmsxcxrrd (bra,ket,mat,dim1,dim2,dim3)
C        Definition: bra(i,j,k)*ket(i,l,k) = mat(j,l) .
C     Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
CC   rmqxxxdd (bra,ket,mat,dim1,dim2,dim3)
CC       Definition: bra(i,j,k)*ket(i,l,k) = mat(j,l) .
CC       Dimensions: mat(dim2,dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C
CC   rmdxxxdd (bra,ket,mat,dim1,dim2,dim3)
CC       Definition: bra(i,j,k)*ket(i,j,k) = mat(j) .
CC       Dimensions: mat(dim2),bra(dim1,dim2,dim3),ket(dim1,dim2,dim3)
C                                                                    
C    rm1hxxxzz (ten,mat,dim1,dim2,dim3)
C        Definition: dconjg(ten(i,j,k))*ten(i,l,k) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),ten(dim1,dim2,dim3)
C
C    rm1hxxxzz_s (ten,mat,dim1,dim2,dim3,index)
C        Definition: index(i,j,k)=x
C                    index(i,l,k)=x1
C                    dconjg(ten(x))*ten(x1) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),index(dim1,dim2,dim3),ten(max(x))
C
C    rmhxxxzz_s (bra,ket,mat,dim1,dim2,dim3,index,index1)
C        Definition: index(i,j,k)=x
C                    index1(i,l,k)=x1
C                    dconjg(bra(x))*ket(x1) = mat(j,l) .
C        Dimensions: mat(dim2,dim2),index(dim1,dim2,dim3),
C                    index1(dim1,dim2,dim3),bra(max(x)),ket(max(x1))
C
C    rmmxxxzz_s (bra,ket,mat,dim1,dim2,dim3,dim4,dim5,index,index1)
C        Definition: index(i,j,k)=x
C                    index1(i,l,k)=x1
C                    dconjg(bra(x))*ket(x1) = mat(j,l) .
C        Dimensions: mat(dim4,dim5),index(dim1,dim2,dim3),
C                    index1(dim1,dim2,dim3),bra(max(x)),ket(max(x1))
C
C **********************************************************************


c-----------------------------------------------------------------------
c Library subroutine rmqxxxzz
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order:
c
c            dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmqxxxzz (bra,ket,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    bra(dim1,dim2,dim3), ket(dim1,dim2,dim3), 
     +              mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+dconjg(bra(i,j,k))*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  mat(j,l)=mat(j,l)+dconjg(bra(1,j,k))*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+dconjg(bra(i,j,1))*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               mat(j,l)=mat(j,l)+dconjg(bra(1,j,1))*ket(1,l,1)
            enddo
         enddo
      endif

      return
      end

C-----------------------------------------------------------------------
C Library subroutine rmqxcxzz
C
C Formation of density type matrix by multiplication of two complex 
C tensors of third order:
C
C            bra(i,j,k)*ket(i,l,k) = mat(j,l).
C
C NB bra is complex conjugate of standard definition i.e. the c.c. is
C    not taken
C
C-----------------------------------------------------------------------

      subroutine rmqxcxzz (bra,ket,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    bra(dim1,dim2,dim3), ket(dim1,dim2,dim3), 
     +              mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+bra(i,j,k)*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  mat(j,l)=mat(j,l)+bra(1,j,k)*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+bra(i,j,1)*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=2,dim2
               mat(j,l)=mat(j,l)+bra(1,j,1)*ket(1,l,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmmxxxzz
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order:
c
c            dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmmxxxzz (bra,ket,mat,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      complex*16    bra(dim1,dim2,dim3), ket(dim1,dim4,dim3), 
     +              mat(dim2,dim4)

      do l=1,dim4
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do j=1,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+dconjg(bra(i,j,k))*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim4
               do j=1,dim2
                  mat(j,l)=mat(j,l)+dconjg(bra(1,j,k))*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim4
            do j=1,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+dconjg(bra(i,j,1))*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim4
            do j=1,dim2
               mat(j,l)=mat(j,l)+dconjg(bra(1,j,1))*ket(1,l,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmmtxxzz
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order, where the transpose of the matrix is stored:
c
c            dconjg(bra(i,j,k))*ket(i,l,k) = mat(l,j).
c
c-----------------------------------------------------------------------

      subroutine rmmtxxzz (bra,ket,mat,dim1,dim2,dim3,dim4)

      implicit none

      integer       dim1,dim2,dim3,dim4, i, j, k, l
      complex*16    bra(dim1,dim2,dim3), ket(dim1,dim4,dim3), 
     +              mat(dim4,dim2)

      do j=1,dim2
         do l=1,dim4
            mat(l,j)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim4
                  do i=1,dim1
                     mat(l,j)=mat(l,j)+dconjg(bra(i,j,k))*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do j=1,dim2
               do l=1,dim4
                  mat(l,j)=mat(l,j)+dconjg(bra(1,j,k))*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim4
               do i=1,dim1
                  mat(l,j)=mat(l,j)+dconjg(bra(i,j,1))*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do j=1,dim2
            do l=1,dim4
               mat(l,j)=mat(l,j)+dconjg(bra(1,j,1))*ket(1,l,1)
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmhxxxzz
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order, where the matrix formed is hermitian:
c
c            dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmhxxxzz (bra,ket,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    bra(dim1,dim2,dim3), ket(dim1,dim2,dim3), 
     +              mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+dconjg(bra(i,j,k))*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  mat(j,l)=mat(j,l)+dconjg(bra(1,j,k))*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+dconjg(bra(i,j,1))*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               mat(j,l)=mat(j,l)+dconjg(bra(1,j,1))*ket(1,l,1)
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         mat(l,l)=dble(mat(l,l))
      enddo
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=dconjg(mat(l,j))
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmaxxxzz
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order, where the matrix formed is anti-hermitian:
c
c            dconjg(bra(i,j,k))*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmaxxxzz (bra,ket,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        x
      complex*16    bra(dim1,dim2,dim3), ket(dim1,dim2,dim3), 
     +              mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+dconjg(bra(i,j,k))*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  mat(j,l)=mat(j,l)+dconjg(bra(1,j,k))*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+dconjg(bra(i,j,1))*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               mat(j,l)=mat(j,l)+dconjg(bra(1,j,1))*ket(1,l,1)
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         x=dimag(mat(l,l))
         mat(l,l)=dcmplx(0.0d0,x)
      enddo
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=-dconjg(mat(l,j))
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmsxcxdd
c
c Formation of density type matrix by multiplication of two real
c tensors of third order, where the matrix formed is symmetric:
c
c            bra(i,j,k)*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmsxcxdd (bra,ket,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*8        bra(dim1,dim2,dim3), ket(dim1,dim2,dim3),
     +              mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+bra(i,j,k)*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  mat(j,l)=mat(j,l)+bra(1,j,k)*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+bra(i,j,1)*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               mat(j,l)=mat(j,l)+bra(1,j,1)*ket(1,l,1)
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=mat(l,j)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmsxcxrrd
c
c Formation of density type matrix by multiplication of two real
c tensors of third order, where the matrix formed is symmetric:
c
c            bra(i,j,k)*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmsxcxrrd (bra,ket,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      real*4        bra(dim1,dim2,dim3), ket(dim1,dim2,dim3)
      real*8        mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+bra(i,j,k)*ket(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  mat(j,l)=mat(j,l)+bra(1,j,k)*ket(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+bra(i,j,1)*ket(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               mat(j,l)=mat(j,l)+bra(1,j,1)*ket(1,l,1)
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=mat(l,j)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmqxxxdd
c
c Formation of density type matrix by multiplication of two real
c tensors of third order:
c
c            bra(i,j,k)*ket(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

C     subroutine rmqxxxdd (bra,ket,mat,dim1,dim2,dim3)

C     implicit none

C     integer       dim1,dim2,dim3, i, j, k, l
C     real*8        bra(dim1,dim2,dim3), ket(dim1,dim2,dim3), 
C    +              mat(dim2,dim2)

C     do l=1,dim2
C        do j=1,dim2
C           mat(j,l)=0.0d0
C        enddo
C     enddo
C     if (dim1.ne.1 .and. dim3.ne.1) then
C        do k=1,dim3
C           do l=1,dim2
C              do j=1,dim2
C                 do i=1,dim1
C                    mat(j,l)=mat(j,l)+bra(i,j,k)*ket(i,l,k)
C                 enddo
C              enddo
C           enddo
C        enddo
C     else if (dim1.eq.1 .and. dim3.ne.1) then
C        do k=1,dim3
C           do l=1,dim2
C              do j=1,dim2
C                 mat(j,l)=mat(j,l)+bra(1,j,k)*ket(1,l,k)
C              enddo
C           enddo
C        enddo
C     else if (dim1.ne.1 .and. dim3.eq.1) then
C        do l=1,dim2
C           do j=1,dim2
C              do i=1,dim1
C                 mat(j,l)=mat(j,l)+bra(i,j,1)*ket(i,l,1)
C              enddo
C           enddo
C        enddo
C     else if (dim1.eq.1 .and. dim3.eq.1) then
C        do l=1,dim2
C           do j=1,dim2
C              mat(j,l)=mat(j,l)+bra(1,j,1)*ket(1,l,1)
C           enddo
C        enddo
C     endif

C     return
C     end

C-----------------------------------------------------------------------
C Library subroutine rmdxxxdd
C
C Formation of density type matrix by multiplication of two complex 
C tensors of third order:
C
C            bra(i,j,k)*ket(i,j,k) = mat(j).
C
C-----------------------------------------------------------------------

C     subroutine rmdxxxdd (bra,ket,mat,dim1,dim2,dim3)

C     implicit none

C     integer       dim1,dim2,dim3, i, j, k
C     real*8        bra(dim1,dim2,dim3), ket(dim1,dim2,dim3), 
C    +              mat(dim2)

C     do j=1,dim2
C        mat(j)=0.0d0
C     enddo
C     if (dim1.ne.1 .and. dim3.ne.1) then
C        do k=1,dim3
C           do j=1,dim2
C              do i=1,dim1
C                 mat(j)=mat(j)+bra(i,j,k)*ket(i,j,k)
C              enddo
C           enddo
C        enddo
C     else if (dim1.eq.1 .and. dim3.ne.1) then
C        do k=1,dim3
C           do j=1,dim2
C              mat(j)=mat(j)+bra(1,j,k)*ket(1,j,k)
C           enddo
C        enddo
C     else if (dim1.ne.1 .and. dim3.eq.1) then
C        do j=1,dim2
C           do i=1,dim1
C              mat(j)=mat(j)+bra(i,j,1)*ket(i,j,1)
C           enddo
C        enddo
C     else if (dim1.eq.1 .and. dim3.eq.1) then
C        do j=1,dim2
C           mat(j)=mat(j)+dconjg(bra(1,j,1))*ket(1,j,1)
C        enddo
C     endif

C     return
C     end




c-----------------------------------------------------------------------
c Library subroutine rm1hxxxzz
c
c Formation of density type matrix by multiplication of a complex 
c tensor of third order with itself. The matrix formed is hermitian:
c
c            dconjg(ten(i,j,k))*ten(i,l,k) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rm1hxxxzz (ten,mat,dim1,dim2,dim3)

      implicit none

      integer       dim1,dim2,dim3, i, j, k, l
      complex*16    ten(dim1,dim2,dim3),mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo
      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     mat(j,l)=mat(j,l)+dconjg(ten(i,j,k))*ten(i,l,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  mat(j,l)=mat(j,l)+dconjg(ten(1,j,k))*ten(1,l,k)
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  mat(j,l)=mat(j,l)+dconjg(ten(i,j,1))*ten(i,l,1)
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               mat(j,l)=mat(j,l)+dconjg(ten(1,j,1))*ten(1,l,1)
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         mat(l,l)=dble(mat(l,l))
      enddo
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=dconjg(mat(l,j))
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rm1hxxxzz_s
c
c Formation of density type matrix by multiplication of a complex 
c tensor of third order with itself. The tensor is not completely
c stored, and the indices are managed by the index tensor. The matrix 
c formed is hermitian:
c
c            index(i,j,k)=x 
c            index(i,l,k)=x1 
c            dconjg(ten(x))*ten(x1) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rm1hxxxzz_s (ten,mat,dim1,dim2,dim3,index)

      implicit none

      integer       dim1,dim2,dim3,i,j,k,l,index(dim1,dim2,dim3),x,x1
      complex*16    ten(*),mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        mat(j,l)=mat(j,l)+dconjg(ten(x))*ten(x1)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  x=index(1,j,k)
                  x1=index(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(ten(x))*ten(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(ten(x))*ten(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               x=index(1,j,1)
               x1=index(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  mat(j,l)=mat(j,l)+dconjg(ten(x))*ten(x1)
               endif
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         mat(l,l)=dble(mat(l,l))
      enddo
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=dconjg(mat(l,j))
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmhxxxzz_s
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order. The tensors are not completely stored, and the 
c indices are managed by the index tensors. The matrix formed is hermitian:
c
c            index(i,j,k)=x 
c            index1(i,l,k)=x1 
c            dconjg(bra(x))*ket(x1) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmhxxxzz_s (bra,ket,mat,dim1,dim2,dim3,index,index1)

      implicit none

      integer       dim1,dim2,dim3,i,j,k,l,index(dim1,dim2,dim3),x,x1,
     +              index1(dim1,dim2,dim3)
      complex*16    bra(*),ket(*),mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=l,dim2
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=l,dim2
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
               endif
            enddo
         enddo
      endif
C
C now form other half of matrix
C
      do l=1,dim2
         mat(l,l)=dble(mat(l,l))
      enddo
      do l=1,dim2
         do j=1,l-1
            mat(j,l)=dconjg(mat(l,j))
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmmxxxzz_s
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order. The tensors are not completely stored, and the 
c indices are managed by the index tensors:
c
c            index(i,j,k)=x 
c            index1(i,l,k)=x1 
c            dconjg(bra(x))*ket(x1) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmmxxxzz_s (bra,ket,mat,dim1,dim2,dim3,dim4,dim5,index,
     +           index1)

      implicit none

      integer       dim1,dim2,dim3,dim4,dim5,i,j,k,l,
     +              index(dim1,dim2,dim3),x,x1,index1(dim1,dim2,dim3)
      complex*16    bra(*),ket(*),mat(dim4,dim5)

      do l=1,dim4
         do j=1,dim5
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim5
               do j=1,dim4
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim5
               do j=1,dim4
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim5
            do j=1,dim4
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim5
            do j=1,dim4
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
               endif
            enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c Library subroutine rmqxxxzz_s
c
c Formation of density type matrix by multiplication of two complex 
c tensors of third order. The tensors are not completely stored, and the 
c indices are managed by the index tensors:
c
c            index(i,j,k)=x 
c            index1(i,l,k)=x1 
c            dconjg(bra(x))*ket(x1) = mat(j,l).
c
c-----------------------------------------------------------------------

      subroutine rmqxxxzz_s (bra,ket,mat,dim1,dim2,dim3,index,
     +           index1)

      implicit none

      integer       dim1,dim2,dim3,i,j,k,l,
     +              index(dim1,dim2,dim3),x,x1,index1(dim1,dim2,dim3)
      complex*16    bra(*),ket(*),mat(dim2,dim2)

      do l=1,dim2
         do j=1,dim2
            mat(j,l)=(0.0d0,0.0d0)
         enddo
      enddo

      if (dim1.ne.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  do i=1,dim1
                     x=index(i,j,k)
                     x1=index1(i,l,k)
                     if (x .ne. 0 .and. x1 .ne. 0) then
                        mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.ne.1) then
         do k=1,dim3
            do l=1,dim2
               do j=1,dim2
                  x=index(1,j,k)
                  x1=index1(1,l,k)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.ne.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               do i=1,dim1
                  x=index(i,j,1)
                  x1=index1(i,l,1)
                  if (x .ne. 0 .and. x1 .ne. 0) then
                     mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
                  endif
               enddo
            enddo
         enddo
      else if (dim1.eq.1 .and. dim3.eq.1) then
         do l=1,dim2
            do j=1,dim2
               x=index(1,j,1)
               x1=index1(1,l,1)
               if (x .ne. 0 .and. x1 .ne. 0) then
                  mat(j,l)=mat(j,l)+dconjg(bra(x))*ket(x1)
               endif
            enddo
         enddo
      endif

      return
      end
