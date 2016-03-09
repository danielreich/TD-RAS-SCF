
C **********************************************************************
C *                                                                    *
C *                        LINEAR EQUATION (lineq.f)                   *
C *                                                                    *
C * Library module storing routines that solve linear equations.       *
C *                                                                    *
C * Contains:                                                          *
C *   lineqz:      Solves the complex linear equation Ax = b employing *
C *                the Gauss algorithm.                                *
C *   choleskyzhp: Decomposition of a hermitian positive definite      *
C *                matrix.                                             *
C *   matinvzhp:   Inverts a hermitian positive definite matrix.       *
C *                                                                    *
C *                                                                    *
C **********************************************************************

C ----------------------------------------------------------------------
C                                 LINEQD
C
C Solves the linear equation Ax = b employing the Gauss algorithm.
C All floating point parameters are double precision
C
C   Input parameters:
C     dim:       Leading dimension of A.
C     n:         True dimension of A, x and b.
C     a:         Coefficient matrix A (destroyed on output).
C     b:         Constant vector b (destroyed on output).
C
C   Output parameters:
C     x:         Solution vector.
C     err:       Error flag having the following meaning:
C                0: everything was o. k.,
C                1: the coefficient matrix is (numerically) singular.
C ----------------------------------------------------------------------

      subroutine lineqd (dim,n,a,b,x,err)

      implicit none

      integer    dim,n,err
      real*8     a(dim,n),b(n),x(n)

      integer    pivrow,i,j,k
      real*8     pivabs,pivabs2,pivot,swap,tmp

c ----------------------------------------------------------------------
c Initialise variables
c ----------------------------------------------------------------------
      err = 0
      k   = 1

c ----------------------------------------------------------------------
c Loop over all columns
c ----------------------------------------------------------------------
 100  continue

c ----------------------------------------------------------------------
c Find pivot element in current column
c ----------------------------------------------------------------------
      pivot  = a(k,k)
      pivrow = k
      pivabs = dble(abs(pivot))
      do i = k+1,n
         pivabs2 = dble(abs(a(i,k)))
         if (pivabs2 .gt. pivabs) then
            pivot  = a(i,k)
            pivrow = i
            pivabs = pivabs2
         endif
      enddo

c ----------------------------------------------------------------------
c Abort if matrix is singular
c ----------------------------------------------------------------------
      if (pivabs .eq. 0.0d0) then
         err = 1
         return
      endif

c ----------------------------------------------------------------------
c Swap current and pivot row
c ----------------------------------------------------------------------
      if (k .ne. pivrow) then
         do j = k,n
            swap        = a(k,j)
            a(k,j)      = a(pivrow,j)
            a(pivrow,j) = swap
         enddo
         swap      = b(k)
         b(k)      = b(pivrow)
         b(pivrow) = swap
      endif

c ----------------------------------------------------------------------
c Eliminate subdiagonal elements in current column
c ----------------------------------------------------------------------
      do j = k+1,n
         tmp = a(k,j)/a(k,k)
         do i = k+1,n
            a(i,j) = a(i,j)-tmp*a(i,k)
         enddo
      enddo
      tmp = b(k)/a(k,k)
      do i = k+1,n
         b(i) = b(i)-tmp*a(i,k)
      enddo

c ----------------------------------------------------------------------
c Continue with next column
c ----------------------------------------------------------------------
      k = k+1
      if (k .le. n) goto 100

c ----------------------------------------------------------------------
c Determine solution vector
c ----------------------------------------------------------------------
      do i = n,1,-1
         x(i) = b(i)
         do j = i+1,n
            x(i) = x(i)-a(i,j)*x(j)
         enddo
         x(i) = x(i)/a(i,i)
      enddo

      return
      end

c #######################################################################

C ----------------------------------------------------------------------
C                                 LINEQZ
C
C Solves the complex linear equation Ax = b employing the Gauss
C algorithm. All floating point parameters are double precision
C (i. e. complex*16).
C
C   Input parameters:
C     dim:       Leading dimension of A.
C     n:         True dimension of A.
C     a:         Coefficient matrix A (destroyed on output).
C     b:         Constant vector b (destroyed on output).
C
C   Output parameters:
C     x:         Solution vector.
C     err:       Error flag having the following meaning:
C                0: everything was o. k.,
C                1: the coefficient matrix is (numerically) singular.
C ----------------------------------------------------------------------

      subroutine lineqz (dim,n,a,b,x,err)

      implicit none

      integer    dim,n,err
      complex*16 a(dim,n),b(n),x(n)

      integer    pivrow,i,j,k
      real*8     pivabs,pivabs2
      complex*16 pivot,swap,tmp

c ----------------------------------------------------------------------
c Initialise variables
c ----------------------------------------------------------------------
      err = 0
      k   = 1

c ----------------------------------------------------------------------
c Loop over all columns
c ----------------------------------------------------------------------
 100  continue

c ----------------------------------------------------------------------
c Find pivot element in current column
c ----------------------------------------------------------------------
      pivot  = a(k,k)
      pivrow = k
      pivabs = dble(abs(pivot))
      do i = k+1,n
         pivabs2 = dble(abs(a(i,k)))
         if (pivabs2 .gt. pivabs) then
            pivot  = a(i,k)
            pivrow = i
            pivabs = pivabs2
         endif
      enddo

c ----------------------------------------------------------------------
c Abort if matrix is singular
c ----------------------------------------------------------------------
      if (pivabs .eq. 0.0d0) then
         err = 1
         return
      endif

c ----------------------------------------------------------------------
c Swap current and pivot row
c ----------------------------------------------------------------------
      if (k .ne. pivrow) then
         do j = k,n
            swap        = a(k,j)
            a(k,j)      = a(pivrow,j)
            a(pivrow,j) = swap
         enddo
         swap      = b(k)
         b(k)      = b(pivrow)
         b(pivrow) = swap
      endif

c ----------------------------------------------------------------------
c Eliminate subdiagonal elements in current column
c ----------------------------------------------------------------------
      do j = k+1,n
         tmp = a(k,j)/a(k,k)
         do i = k+1,n
            a(i,j) = a(i,j)-tmp*a(i,k)
         enddo
      enddo
      tmp = b(k)/a(k,k)
      do i = k+1,n
         b(i) = b(i)-tmp*a(i,k)
      enddo

c ----------------------------------------------------------------------
c Continue with next column
c ----------------------------------------------------------------------
      k = k+1
      if (k .le. n) goto 100

c ----------------------------------------------------------------------
c Determine solution vector
c ----------------------------------------------------------------------
      do i = n,1,-1
         x(i) = b(i)
         do j = i+1,n
            x(i) = x(i)-a(i,j)*x(j)
         enddo
         x(i) = x(i)/a(i,i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C                               CHOLESKYZHP
C
C Decomposites a complex hermitian (H), positive definite (P) matrix
C A into a product of a left lower triangular matrix L and its adjoint,
C i. e. A = L L^dagger where L_ij = 0 for j > i and L_ii is real and
C positive. All calculations are performed with double precision (Z).
C
C Input parameters:
C   a:     The matrix to be decomposited (a(i,j) is needed for j <= i
C          only).
C   dim:   The number of allocated rows (and columns) of "a".
C   n:     The number of used rows (and columns) of "a".
C   tol:   A (relative) tolerance for the check wether "a" is hermitian
C          and positive definite.
C
C Output parameters:
C   a:     The left lower triangular matrix L. (The part of "a" right of
C          the diagonal still contains the corresponding part of A.)
C   error: An error flag having the following meaning:
C          0: Everything was o. k.
C          1: The matrix is not hermitian and positive definite.
C-----------------------------------------------------------------------

      subroutine choleskyzhp (a,dim,n,tol,error)

      implicit none

      integer    dim,n,error,i,j,k
      real*8     tol,y
      complex*16 a(dim,dim),x

      error = 0
      do j = 1,n
         do i = j,n
            x = a(i,j)
            do k = 1,j-1
               x = x-a(i,k)*dconjg(a(j,k))
            enddo
            if (i .eq. j) then
               if (dabs(dimag(x)) .gt. tol*dble(x)) then
                  error = 1
                  return
               else
                  a(i,j) = dcmplx(dsqrt(dble(x)))
                  y = 1.0d0/dble(a(i,j))
               endif
            else
               a(i,j) = x*y
            endif
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C                               CHOLESKYDSP
C
C Decomposites a real symmetric (S), positive definite (P) matrix
C A into a product of a left lower triangular matrix L and its adjoint,
C i. e. A = L L^dagger where L_ij = 0 for j > i and L_ii is real and
C positive. All calculations are performed with double precision (D).
C
C Input parameters:
C   a:     The matrix to be decomposited (a(i,j) is needed for j <= i
C          only).
C   dim:   The number of allocated rows (and columns) of "a".
C   n:     The number of used rows (and columns) of "a".
C   tol:   A (relative) tolerance for the check wether "a" is symmetric
C          and positive definite.
C
C Output parameters:
C   a:     The left lower triangular matrix L. (The part of "a" right of
C          the diagonal still contains the corresponding part of A.)
C   error: An error flag having the following meaning:
C          0: Everything was o. k.
C          1: The matrix is not symmetric and positive definite.
C-----------------------------------------------------------------------

      subroutine choleskydsp (a,dim,n,error)

      implicit none

      integer    dim,n,error,i,j,k
      real*8     y,a(dim,dim),x

      error = 0
      do j = 1,n
         do i = j,n
            x = a(i,j)
            do k = 1,j-1
               x = x-a(i,k)*a(j,k)
            enddo
            if (i .eq. j) then
               if (x.le.0.0) then
                  error=1
                  return
               else
                  a(i,j) = sqrt(x)
                  y = 1.0d0/a(i,j)
               endif
            else
               a(i,j) = x*y
            endif
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C                               MATINVZHP
C
C Inverts a double precision complex (Z), Hermitian (H), and positive
C definite (P) matrix A employing a Cholesky-decomposition of A.
C
C Input parameters:
C   a:     The matrix to be inverted (a(i,j) is needed for j <= i only).
C   dim:   The number of allocated rows (and columns) of "a".
C   n:     The number of used rows (and columns) of "a".
C   tol:   A (relative) tolerance for the check whether "a" is Hermitian
C          and positive definite.
C Output parameters:
C   b:     The inverse of "a".
C   a:     The matrix "a" is overwritten.
C   error: An error flag having the following meaning:
C          0: Everything was o. k.
C          1: The matrix "a" is not hermitian and positive definite.
C Remarks:
C   Library routine "choleskyzhp" is used.
C-----------------------------------------------------------------------

      subroutine matinvzhp (b,a,dim,n,tol,error)

      implicit none

      complex*16 zeroz,onez
      parameter (zeroz = (0.0d0,0.0d0),onez = (1.0d0,0.0d0))

      integer    dim,n,error,i,j,k
      real*8     tol
      complex*16 b(dim,dim),a(dim,dim),x

C 1. Construct the A = L L^dagger Cholesky-decomposition.
C 2. For every column j of the identity matrix do:
C    a. Solve L b = 1 (first do-loop over i).
C    b. Solve L_dagger u = b (u is stored directly in the j-th column of
C       the matrix b) (second do-loop over i).

      call choleskyzhp(a,dim,n,tol,error)

      if (error .ne. 0) return
      do j = 1,n
         do i = 1,n
            if (i .lt. j) then
               b(i,j) = zeroz
            elseif (i .eq. j) then
               b(i,j) = onez/a(i,i)
            else
               x = zeroz
               do k = j,i-1
                  x = x-a(i,k)*b(k,j)
               enddo
               b(i,j) = x/a(i,i)
            endif
         enddo
         do i = n,1,-1
            x = b(i,j)
            do k = i+1,n
               x = x-dconjg(a(k,i))*b(k,j)
            enddo
            b(i,j) = x/a(i,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C                                 LINEQCHD
C
C Solves the linear equation Ax = b employing the Cholesky
C decomposition. All floating point parameters are double precision
C
C   Input parameters:
C     dim:       Leading dimension of A.
C     n:         True dimension of A, x and b.
C     a:         Coefficient matrix A (destroyed on output).
C     b:         Constant vector b (destroyed on output), holds solution
C
C   Output parameters:
C     x:         auxiliary vector.
C     err:       Error flag having the following meaning:
C                0: everything was o. k.,
C                1: the coefficient matrix is (numerically) singular.
C ----------------------------------------------------------------------

      subroutine lineqchd (dim,n,a,c,b,err)

      implicit none

      integer    dim,n,err
      real*8     a(dim,n),b(n),c(n),x

      integer    i,j

c ----------------------------------------------------------------------
c Initialise variables
c ----------------------------------------------------------------------
      err = 0

C-----------------------------------------------------------------------
C Do the cholesky decomposition
C-----------------------------------------------------------------------
      call choleskydsp(a,dim,n,err)
      if (err .ne. 0) return

C-----------------------------------------------------------------------
C Solve system G*x`=b
C-----------------------------------------------------------------------
      do i=1,n
         x=c(i)
         do j=1,i-1
            x=x-a(i,j)*b(j)
         enddo
         b(i)=x/a(i,i)
      enddo

      do i=1,n
         x=0.0
         do j=1,i
            x=x+b(j)*a(i,j)
         enddo
         c(i)=x
      enddo

c ----------------------------------------------------------------------
c Determine solution vector: Solve system G^T*x=x`
c ----------------------------------------------------------------------
      do i=0,n-1
         x=b(n-i)
         do j=n-i+1,n
            x=x-a(j,n-i)*c(j)
         enddo
         c(n-i)=x/a(n-i,n-i)
      enddo

      do i=1,n
         x=0.0
         do j=i,n
            x=x+a(j,i)*c(j)
         enddo
         b(i)=x
      enddo

      return
      end
