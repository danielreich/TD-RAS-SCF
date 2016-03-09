C-----------------------------------------------------------------------
C                      OP2LIB
C  SUBROUTINES ACTING ON TWO OBJECTS OF THE SAME TYPE (LEVEL 2 ROUTINES)
C
C NOMENCLATURE:
C    Each of the following routines starts with:
C       add:  the objects are added
C       sub:  the second object is subtracted from the first
C       sum2: the products of the elements are summed
C       swap: the objects are swapped
C    This is followed by three chracters:
C    Character 1 denotes the objects on which the action is performed:
C       s: scalar
C       v: vector
C       q: quadratic matrix
C       m: general (rectangular) matrix
C       h: hermitian matrix
C       x: this is a blank
C    Character 2 specifies information about how first object is used:
C       x: this is a blank
C       a: first entry is taken as complex-conjugate
C    Character 3 specifies information about how second objects is used:
C       x: this is a blank
C       a: second entry is taken as complex-conjugate
C    Character 4 defines the data type of the objects:
C       i: integer
C       s: real single precision (real*4)
C       d: real double precision (real*8)
C       c: complex single precision (complex*8)
C       z: complex double precision (complex*16)
C    Further characters if present give more information:
C       o: the result is stored (overwritten) in second input object
C       o1: the result is stored (overwritten) in first input object
C
C Contents:
C In the following list of available subroutines, objects on the LHS of
C the definition are input, that on the RHS output. The usual summation
C convention is used i.e. a sum is made over repeated indices on the LHS
C
C    addmxxz (a,b,c,dim1,dim2)
C       Definition: a(i,j) + b(i,j) = c(i,j)
C       Dimensions: a(dim1,dim2),b(dim1,dim2),c(dim1,dim2)
C
C    addmxazo1(a,b,dim1,dim2)
C       Definition: a(i,j) = a(i,j) + conjg(b(j,i))
C       Dimensions: a(dim2,dim1), b(dim1,dim2)
C
C    submxxz (a,b,c,dim1,dim2)
C       Definition: a(j,i) - b(j,i) = c(j,i)
C       Dimensions: a(dim1,dim2),b(dim1,dim2),c(dim1,dim2)
C
C    submxazo1(a,b,dim1,dim2)
C       Definition: a(i,j) = a(i,j) - conjg(b(j,i))
C       Dimensions: a(dim2,dim1), b(dim1,dim2)
C
C    sum2qxxz (a,b,x,dim)
C       Definition: a(j,i)*b(j,i) = x
C       Dimensions: a(dim,dim),b(dim,dim),x
C
C    addmxxdo (a,b,dim1,dim2)
C       Definition: a(j,i) + b(j,i) = b(j,i)
C       Dimensions: a(dim1,dim2),b(dim1,dim2)
C
C    addmxxzo (a,b,dim1,dim2)
C       Definition: a(j,i) + b(j,i) = b(j,i)
C       Dimensions: a(dim1,dim2),b(dim1,dim2)
C
C    addvxxdo (a,b,dim)
C       Definition: a(i) + b(i) = b(i)
C       Dimensions: a(dim),b(dim)
C
C    addvxxzo (a,b,dim)
C       Definition: a(i) + b(i) = b(i)
C       Dimensions: a(dim),b(dim)
C
C    addvcxzo (a,b,dim)
C       Definition: conjg(a(i)) + b(i) = b(i)
C       Dimensions: a(dim), b(dim)
C
C    submxxzo1 (a,b,dim1,dim2)
C       Definition: a(j,i) - b(j,i) = a(j,i)
C       Dimensions: a(dim1,dim2),b(dim1,dim2)
C
C    subvxxz (a,b,c,dim)
C       Definition: a(i) - b(i) = c(i)
C       Dimensions: a(dim),b(dim),c(dim)
C
C    subvxxzo1 (a,b,dim)
C       Definition: a(i) - b(i) = a(i)
C       Dimensions: a(dim),b(dim)
C
C    subvcxzo (a,b,dim)
C       Definition: b(i) - conjg(a(i)) = b(i)
C       Dimensions: a(dim), b(dim)
C
C    swapsxxi (i,j)
C       Definition: i <-->j 
C
C    swapsxxd (a,b)
C       Definition: a <--> b
C
C
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C Library subroutine addmxxz
C
C adds two complex rectangular matrices together
C     a(i,j) + b(i,j) = c(i,j)
C-----------------------------------------------------------------------

      subroutine addmxxz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,i) + b(j,i)
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine addmxazo1
C
C     Definition: a(i,j) = a(i,j) + conjg(b(j,i))
C     Dimensions: a(dim2,dim1), b(dim1,dim2)
C-----------------------------------------------------------------------

      subroutine  addmxazo1(a,b,dim1,dim2)

      implicit none

      integer dim1, dim2, i, j
      complex*16 a(dim2,dim1), b(dim1,dim2)

      do i=1,dim1
         do j=1,dim2
            a(j,i) = a(j,i) + conjg(b(i,j))
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine submxxz
C
C subtracts one complex rectangular matrix from another
C     a(i,j) - b(i,j) = c(i,j)
C-----------------------------------------------------------------------

      subroutine submxxz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,i) - b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine submxazo1
C
C     Definition: a(j,i) - b(j,i) = a(j,i)
C     Dimensions: a(dim1,dim2),b(dim1,dim2)
C-----------------------------------------------------------------------

      subroutine  submxazo1(a,b,dim1,dim2)

      implicit none

      integer dim1, dim2, i, j
      complex*16 a(dim2,dim1), b(dim1,dim2)

      do i=1,dim1
         do j=1,dim2
            a(j,i) = a(j,i) - conjg(b(i,j))
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine sum2qxxz
C
C sums the product of elements of two complex quadratic matrices
C     a(j,i)*b(j,i) = x
C-----------------------------------------------------------------------

      subroutine sum2qxxz (a,b,x,dim)

      implicit none

      integer     dim,i,j
      complex*16  a(dim,dim),b(dim,dim),x

      x=(0.0d0,0.0d0)
      do i = 1,dim
         do j = 1,dim
            x = x + a(j,i)*b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine addmxxdo
C
C adds two real rectangular matrices together
C     a(i,j) + b(i,j) = b(i,j)
C
C NB input matrix b is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine addmxxdo (a,b,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      real*8  a(dim1,dim2),b(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            b(j,i) = a(j,i) + b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine addvxxdo
C
C adds two real vectors together
C     a(i) + b(i) = b(i)
C
C NB input vector b is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine addvxxdo (a,b,dim)

      implicit none

      integer     dim,i
      real*8  a(dim),b(dim)

      do i = 1,dim
         b(i) = a(i) + b(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine addvxxzo
C
C adds two complex vectors together
C     a(i) + b(i) = b(i)
C
C NB input vector b is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine addvxxzo (a,b,dim)

      implicit none

      integer     dim,i
      complex*16  a(dim),b(dim)

      do i = 1,dim
         b(i) = a(i) + b(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine addvcxzo
C
C adds two complex vectors together
C     conjg(a(i)) + b(i) = b(i)
C
C NB input vector b is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine addvcxzo (a,b,dim)

      implicit none

      integer     dim,i
      complex*16  a(dim),b(dim)

      do i = 1,dim
         b(i) = conjg(a(i)) + b(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine addmxxzo
C
C adds two complex rectangular matrices together
C     a(i,j) + b(i,j) = b(i,j)
C
C NB input matrix b is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine addmxxzo (a,b,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),b(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            b(j,i) = a(j,i) + b(j,i)
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine submxxzo1
C
C subtracts a complex rectangular matrix from another
C     a(i,j) - b(i,j) = a(i,j)
C
C NB input matrix a is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine submxxzo1 (a,b,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),b(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            a(j,i) = a(j,i) - b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine subvxxz
C 
C subtracts a complex vector from another
C     a(i) - b(i) = c(i)
C 
C-----------------------------------------------------------------------

      subroutine subvxxz (a,b,c,dim)

      implicit none

      integer     dim,i
      complex*16  a(dim),b(dim),c(dim)

      do i = 1,dim
         c(i) = a(i) - b(i)
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine subvxxzo1
C
C subtracts a complex vector from another
C     a(i) - b(i) = a(i)
C
C NB input vector a is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine subvxxzo1 (a,b,dim)

      implicit none

      integer     dim,i
      complex*16  a(dim),b(dim)

      do i = 1,dim
         a(i) = a(i) - b(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine subvcxzo
C
C adds two complex vectors together
C     b(i) - conjg(a(i)) = b(i)
C
C NB input vector b is overwritten and holds output at end of routine
C-----------------------------------------------------------------------

      subroutine subvcxzo (a,b,dim)

      implicit none

      integer     dim,i
      complex*16  a(dim),b(dim)

      do i = 1,dim
         b(i) = b(i) - conjg(a(i)) 
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine swapsxxi
C
C swappes two integers: a <--> b
C
C-----------------------------------------------------------------------

      subroutine swapsxxi (int1,int2)

      implicit none

      integer int1,int2,aux

      aux  = int1
      int1 = int2
      int2 = aux

      return
      end

C-----------------------------------------------------------------------
C Library subroutine swapsxxd
C
C swappes two real variables: a <--> b
C
C-----------------------------------------------------------------------

      subroutine swapsxxd (x1,x2)

      implicit none

      real*8 x1,x2,aux

      aux  = x1
      x1 = x2
      x2 = aux

      return
      end

