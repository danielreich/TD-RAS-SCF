************************************************************************
*
* Library for arithmetic with side-diagonal matrices.
* (Currently only for quadratic matrices.)
* FO 08/2006
*
* The side-diagonal matrices are passed as a pair of arguments:
*   integer k        : specifies which side-diagonal is populated.
*                      k>0 is below the diagonal,
*                      k<0 is above the diagonal.
*   <numeric_type> a : array which stores the side-diagonal elements.
*
* It's possible to use diagonal matrices (stored as vector) as input,
* simply pass k=0.
*
* Naming Convention (same as in xvlib):
*   First two characters specify the objects being multiplied
*     o : side-diagonal matrix (stored as vector)
*     t : tensor of rank 3
*     v : vector
*   Next two characters specify how the objetcs are used:
*     x : unchanged
*   Next two characters specify numeric type:
*     d : real*8
*     y : complex*16, stored as two real*8 arrays
*     z : complex*16
*   Further characters give more options:
*     t : compute the tensor product, not the matrix product
*     o : result overwrites second input object
*   
************************************************************************



* ----------------------------------------------------------------------
      subroutine ooxxdd(a,ka,b,kb,c,k,n)
* ----------------------------------------------------------------------
* Multiplication of two side-diagonal matrices a and b with shifts
* ka and kb, respectively. Result is another side-diagonal matrix c
* with shift k. All matrices are nxn.
* ----------------------------------------------------------------------
      implicit none
      integer  ka,kb,k,n,i
      real*8   a(n),b(n),c(n)

      k=ka+kb
      if (ka.ge.0) then
         if (kb.ge.0) then
            do i=1,n-k
               c(i)=a(i+kb)*b(i)
            enddo
         else
            if (k.ge.0) then
               do i=1,-kb
                  c(i)=0.d0
               enddo
               do i=-kb+1,n-k
                  c(i)=a(i+kb)*b(i+kb)
               enddo
            else
               do i=1,ka
                  c(i)=0.d0
               enddo
               do i=ka+1,n+k
                  c(i)=a(i-ka)*b(i-ka)
               enddo
            endif
         endif
      else 
         if (kb.ge.0) then
            if (k.ge.0) then
               do i=1,n-kb
                  c(i)=a(i+k)*b(i)
               enddo
               do i=n-kb+1,n-k
                  c(i)=0.d0
               enddo
            else
               do i=1,n+ka
                  c(i)=a(i)*b(i-k)
               enddo
               do i=n+ka+1,n+k
                  c(i)=0.d0
               enddo
            endif
         else
            do i=1,n+k
               c(i)=a(i)*b(i-ka)
            enddo
         endif
      endif

      return
      end ! of ooxxdd



* ----------------------------------------------------------------------
      subroutine ooxxddo(a,ka,b,kb,n)
* ----------------------------------------------------------------------
* Multiplication of two side-diagonal nxn matrices a and b with shifts
* ka and kb, respectively. Result is another side-diagonal matrix, it is
* returned in b/kb.
* This is the same as ooxxdd, but some loops must be done backwards.
* ----------------------------------------------------------------------
      implicit none
      integer  ka,kb,n,k,i
      real*8   a(n),b(n)

      k=ka+kb
      if (ka.ge.0) then
         if (kb.ge.0) then
            do i=1,n-k
               b(i)=a(i+kb)*b(i)
            enddo
         else
            if (k.ge.0) then
               do i=n-k,-kb+1,-1
                  b(i)=a(i+kb)*b(i+kb)
               enddo
               do i=-kb,1,-1
                  b(i)=0.d0
               enddo
            else
               do i=n+k,ka+1,-1
                  b(i)=a(i-ka)*b(i-ka)
               enddo
               do i=ka,1,-1
                  b(i)=0.d0
               enddo
            endif
         endif
      else 
         if (kb.ge.0) then
            if (k.ge.0) then
               do i=1,n-kb
                  b(i)=a(i+k)*b(i)
               enddo
               do i=n-kb+1,n-k
                  b(i)=0.d0
               enddo
            else
               do i=1,n+ka
                  b(i)=a(i)*b(i-k)
               enddo
               do i=n+ka+1,n+k
                  b(i)=0.d0
               enddo
            endif
         else
            do i=1,n+k
               b(i)=a(i)*b(i-ka)
            enddo
         endif
      endif
      kb=k

      return
      end ! of ooxxddo



* ----------------------------------------------------------------------
      subroutine ooxxddt(a,ka,b,kb,c,k,n,m)
* ----------------------------------------------------------------------
* Tensor product of an nxn side-diagonal (shift ka) matrix a and
* an mxm side-diagonal (shift kb) matrix b. Result is an (nm)x(nm)
* side-diagonal (shift k) matrix c.
* ----------------------------------------------------------------------
      implicit none
      integer ka,kb,k,n,m,i,ia,ib
      real*8  a(n),b(m),c(n*m)

      k=n*kb+ka
      if (k.ge.0) then
         do i=1,n*m-k
            ia=mod(i-1,n)+1
            ib=(i-1)/n+1
            if (ka.ge.0) then
               if ((ia.le.n-ka).and.(ib.le.m-kb)) then
                  c(i)=a(ia)*b(ib)
               else
                  c(i)=0.d0
               endif
            else
               if ((ia+ka.gt.0).and.(ib.le.m-kb)) then
                  c(i)=a(ia+ka)*b(ib)
               else
                  c(i)=0.d0
               endif
            endif
         enddo
      else
         do i=1,n*m+k
            ia=mod(i-1,n)+1
            ib=(i-1)/n+1
            if (ka.le.0) then
               if ((ia.le.n+ka).and.(ib.le.m+kb)) then
                  c(i)=a(ia)*b(ib)
               else
                  c(i)=0.d0
               endif
            else
               if ((ia-ka.gt.0).and.(ib.le.m+kb)) then
                  c(i)=a(ia-ka)*b(ib)
               else
                  c(i)=0.d0
               endif
            endif
         enddo
      endif

      return
      end ! of ooxxddt



* ----------------------------------------------------------------------
      subroutine ooxxzzt(a,ka,b,kb,c,k,n,m)
* ----------------------------------------------------------------------
* Tensor product of an nxn side-diagonal (shift ka) matrix a and
* an mxm side-diagonal (shift kb) matrix b. Result is an (nm)x(nm)
* side-diagonal (shift k) matrix c.
* ----------------------------------------------------------------------
      implicit none
      integer ka,kb,k,n,m,i,ia,ib
      complex*16 a(n),b(m),c(n*m)

      k=n*kb+ka
      if (k.ge.0) then
         do i=1,n*m-k
            ia=mod(i-1,n)+1
            ib=(i-1)/n+1
            if (ka.ge.0) then
               if ((ia.le.n-ka).and.(ib.le.m-kb)) then
                  c(i)=a(ia)*b(ib)
               else
                  c(i)=(0.d0,0.d0)
               endif
            else
               if ((ia+ka.gt.0).and.(ib.le.m-kb)) then
                  c(i)=a(ia+ka)*b(ib)
               else
                  c(i)=(0.d0,0.d0)
               endif
            endif
         enddo
      else
         do i=1,n*m+k
            ia=mod(i-1,n)+1
            ib=(i-1)/n+1
            if (ka.le.0) then
               if ((ia.le.n+ka).and.(ib.le.m+kb)) then
                  c(i)=a(ia)*b(ib)
               else
                  c(i)=(0.d0,0.d0)
               endif
            else
               if ((ia-ka.gt.0).and.(ib.le.m+kb)) then
                  c(i)=a(ia-ka)*b(ib)
               else
                  c(i)=(0.d0,0.d0)
               endif
            endif
         enddo
      endif

      return
      end ! of ooxxzzt



* ----------------------------------------------------------------------
      subroutine ovxxdd(a,k,v,w,n)
* ----------------------------------------------------------------------
* Multiplication of an nxn side-diagonal (shift k) matrix a with a
* vector v. Result is vector w.
* ----------------------------------------------------------------------
      implicit none
      integer  k,n,i
      real*8   a(n),v(n),w(n)

      if (k.ge.0) then
         do i=1,k
            w(i)=0.d0
         enddo
         do i=k+1,n
            w(i)=a(i-k)*v(i-k)
         enddo
      else
         do i=1,n+k
            w(i)=a(i)*v(i-k)
         enddo
         do i=n+k+1,n
            w(i)=0.d0
         enddo
      endif

      return
      end ! of ovxxdd



* ----------------------------------------------------------------------
      subroutine ovxxddo(a,k,v,n)
* ----------------------------------------------------------------------
* Multiplication of an nxn side-diagonal (shift k) matrix a with a
* vector v. Result is returned in v.
* ----------------------------------------------------------------------
      implicit none
      integer  k,n,i
      real*8   a(n),v(n)

      if (k.ge.0) then
         do i=n,k+1,-1
            v(i)=a(i-k)*v(i-k)
         enddo
         do i=k,1,-1
            v(i)=0.d0
         enddo
      else
         do i=1,n+k
            v(i)=a(i)*v(i-k)
         enddo
         do i=n+k+1,n
            v(i)=0.d0
         enddo
      endif

      return
      end ! of ovxxddo



* ----------------------------------------------------------------------
      subroutine ovxxdz(a,k,v,w,n)
* ----------------------------------------------------------------------
* Multiplication of an nxn side-diagonal (shift k) matrix a with a
* complex vector v. Result is vector w.
* ----------------------------------------------------------------------
      implicit none
      integer  k,n,i
      real*8   a(n)
      complex*16 v(n),w(n)

      if (k.ge.0) then
         do i=1,k
            w(i)=(0.d0,0.d0)
         enddo
         do i=k+1,n
            w(i)=a(i-k)*v(i-k)
         enddo
      else
         do i=1,n+k
            w(i)=a(i)*v(i-k)
         enddo
         do i=n+k+1,n
            w(i)=(0.d0,0.d0)
         enddo
      endif

      return
      end ! of ovxxdz



* ----------------------------------------------------------------------
      subroutine otxxdz(a,k,v,w,vdim,n,ndim)
* ----------------------------------------------------------------------
* Multiplication of an nxn side-diagonal (shift k) matrix a with a
* complex tensor of 3rd order. Result is the tensor w.
* ----------------------------------------------------------------------
      implicit none
      integer  k,n,vdim,ndim,i,vi,ni
      real*8   a(n)
      complex*16 v(vdim,n,ndim),w(vdim,n,ndim)
        
      if (k.ge.0) then
         do ni=1,ndim
            do i=1,k
               do vi=1,vdim
                  w(vi,i,ni)=(0.d0,0.d0)
               enddo
            enddo
            do i=k+1,n
               do vi=1,vdim
                  w(vi,i,ni)=a(i-k)*v(vi,i-k,ni)
               enddo
            enddo
         enddo
      else
         do ni=1,ndim
            do i=1,n+k
               do vi=1,vdim
                  w(vi,i,ni)=a(i)*v(vi,i-k,ni)
               enddo
            enddo
            do i=n+k+1,n
               do vi=1,vdim
                  w(vi,i,ni)=(0.d0,0.d0)
               enddo
            enddo
         enddo
      endif

      return
      end ! of otxxdz



* ----------------------------------------------------------------------
      subroutine ovxxyz(a,k,v,w,n)
* ----------------------------------------------------------------------
* Multiplication of an nxn side-diagonal (shift k) matrix a with a
* complex vector v, where a is complex but stored as two real arrays.
* Result is vector w.
* ----------------------------------------------------------------------
      implicit none
      integer  k,n,i
      real*8   a(n,2)
      complex*16 v(n),w(n)

      if (k.ge.0) then
         do i=1,k
            w(i)=(0.d0,0.d0)
         enddo
         do i=k+1,n
            w(i)=dcmplx(a(i-k,1),a(i-k,2))*v(i-k)
         enddo
      else
         do i=1,n+k
            w(i)=dcmplx(a(i,1),a(i,2))*v(i-k)
         enddo
         do i=n+k+1,n
            w(i)=(0.d0,0.d0)
         enddo
      endif

      return
      end ! of ovxxyz



* ----------------------------------------------------------------------
      subroutine otxxyz(a,k,v,w,vdim,n,ndim)
* ----------------------------------------------------------------------
* Multiplication of an nxn side-diagonal (shift k) matrix a with a
* complex tensor of 3rd order, where a is complex but stored as two
* real arrays. Result is tensor w.
* ----------------------------------------------------------------------
      implicit none
      integer  k,n,vdim,ndim,i,vi,ni
      real*8   a(n,2)
      complex*16 v(vdim,n,ndim),w(vdim,n,ndim)

      if (k.ge.0) then
         do ni=1,ndim
            do i=1,k
               do vi=1,vdim
                  w(vi,i,ni)=(0.d0,0.d0)
               enddo
            enddo
            do i=k+1,n
               do vi=1,vdim
                  w(vi,i,ni)=dcmplx(a(i-k,1),a(i-k,2))*v(vi,i-k,ni)
               enddo
            enddo
         enddo
      else
         do ni=1,ndim
            do i=1,n+k
               do vi=1,vdim
                  w(vi,i,ni)=dcmplx(a(i,1),a(i,2))*v(vi,i-k,ni)
               enddo
            enddo
            do i=n+k+1,n
               do vi=1,vdim
                  w(vi,i,ni)=(0.d0,0.d0)
               enddo
            enddo
         enddo
      endif

      return
      end ! of otxxyz



* ----------------------------------------------------------------------
      subroutine omattens(B,A,T,k,n,p,q,add)
* ----------------------------------------------------------------------
* Multiplication of a complex matrix with a complex tensor of third
* order, similar to "mattens". But while in mattens T corresponds to
* a matrix that is diagonal in its first index, here T corresponds to
* a side-diagonal matrix with shift k.
* ----------------------------------------------------------------------
      implicit none
      integer    k,n,p,q,g,i,j
      complex*16 A(n,q),B(n,p),T(n,p,q)
      logical    add

      if (.not.add) then
         do i=1,p
            do g=1,n
               B(g,i)=(0.0d0,0.0d0)
            enddo
         enddo
      endif

      if (k.ge.0) then
         do j=1,q
            do i=1,p
               do g=k+1,n
                  B(g,i)=B(g,i)+T(g-k,i,j)*A(g-k,j)
               enddo
            enddo
         enddo
      else
         do j=1,q
            do i=1,p
               do g=1,n+k
                  B(g,i)=B(g,i)+T(g,i,j)*A(g-k,j)
               enddo
            enddo
         enddo
      endif

      return
      end ! of omattens

* ----------------------------------------------------------------------
