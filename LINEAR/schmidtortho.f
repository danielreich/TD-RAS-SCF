C########################################################################
C                                                                       #
C                        SCHMIDTORTHO                                   #
C                                                                       #
C Library module making Schmidt-orthogonalization.                      #
C                                                                       #
C Contains:                                                             #
C       schmidtortho: Schmidt-orthogonalizes a set of complex vectors   #
C                     in column form.                                   #
C                                                                       #
C V7.0 MB                                                               #
C########################################################################


C-----------------------------------------------------------------------
C                         SCHMIDTORTHO
C
C Schmidt-orthogonalizes a set of complex vectors,
C (e.g. single-particle functions).
C The vectors are the column-vectors of psi(gdim,dim)
C-----------------------------------------------------------------------

      subroutine schmidtortho (psi,gdim,dim,ierr)

      implicit none

      integer    gdim,dim,e,e1,ierr
      complex*16 psi(gdim,dim),overlap
      real*8     norm

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the set of functions.
C The orthonormalization is made twice to remove numerical inaccuracies.
C-----------------------------------------------------------------------
      ierr = 0
      do e=1,dim
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxzzr(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         if (norm .le. 1.0d-99 ) norm = 1.0d99
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      do e=1,dim
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxzzr(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         if( norm .le. 0.8d0 ) then
            ierr = e 
            return
         end if
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      return
      end


C-----------------------------------------------------------------------
C                            SCHMIDTORTHOD
C
C Schmidt-orthogonalizes (row-orthogonalization) a real quadratic
C matrix psi(dim,dim).
C
C SS 11/98
C-----------------------------------------------------------------------

      subroutine schmidtorthod (psi,dim)

      implicit none

      integer    dim,e,e1,i
      real*8 psi(dim,dim),overlap,norm

      call tranqxd(psi,dim)

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the matrix.
C The orthonormalization is made twice to remove numerical inaccuracies.
C-----------------------------------------------------------------------
      do i=1,2
         do e=1,dim
            do e1=1,e-1
               call vvtxdd(psi(1,e1),psi(1,e),overlap,dim)
               call xvxxddr(overlap,psi(1,e1),psi(1,e),dim)
            enddo
            call normvxd(psi(1,e),norm,dim)
c            if(abs(norm-1.d0) .gt. 0.01d0) norm=1.d99
            call xvixddo(norm,psi(1,e),dim)
         enddo
      enddo

      call tranqxd(psi,dim)

      return
      end


C-----------------------------------------------------------------------
C                         SCHMIDTORTHO2
C
C Schmidt-orthogonalizes a set of complex vectors,
C (e.g. single-particle functions).
C The vectors are the column-vectors of psi(gdim,dim)
C This routine is identical to schmidtortho on lib/linear/schmidtortho.f
C except for the additional parameter mnorm (minmal norm).
C HDM 10/05
C-----------------------------------------------------------------------

      subroutine schmidtortho2(psi,gdim,dim,ierr,mnorm)

      implicit none

      integer    gdim,dim,e,e1,ierr
      complex*16 psi(gdim,dim),overlap
      real*8     norm, mnorm

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the set of functions.
C The orthonormalization is made twice to remove numerical inaccuracies.
C-----------------------------------------------------------------------
      mnorm = 1.d9
      ierr = 0
      do e=1,dim
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxzzr(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         mnorm = min(norm,mnorm)
         if (norm .le. 1.0d-99 ) norm = 1.0d99
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      do e=1,dim
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxzzr(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         if( norm .le. 0.8d0 ) then
            ierr = e
            return
         end if
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      return
      end


C-----------------------------------------------------------------------
C                         SCHMIDTORTHO3
C
C Schmidt-orthogonalizes a single vector phi against a set of 
C complex vectors, (e.g. single-particle functions).
C The vectors are the column-vectors of psi(gdim,dim)
C This routine is similar to schmidtortho on lib/linear/schmidtortho.f
C HDM 10/06
C-----------------------------------------------------------------------

      subroutine schmidtortho3(psi,phi,gdim,dim,norm)

      implicit none

      integer    gdim,dim,e
      real*8     norm
      complex*16 psi(gdim,dim),phi(gdim),overlap

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the set of functions.
C The orthonormalization is made only once.
C-----------------------------------------------------------------------
      do e=1,dim
         call vvaxzz(psi(1,e),phi,overlap,gdim)
         call xvxxzzr(overlap,psi(1,e),phi,gdim)
      enddo
      call normvxz(phi,norm,gdim)
      if (norm .ge. 1.d-50 )
     +  call xvixdzo(norm,phi,gdim)

      return
      end

C-----------------------------------------------------------------------
C                         SCHMIDTORTHOR
C
C Schmidt-orthogonalizes a set of real vectors,
C (e.g. natpots).
C The vectors are the column-vectors of psi(gdim,dim)
C This routine is identical to schmidtortho2 on lib/linear/schmidtortho.f
C except for real*8 arithmetic.
C HDM 10/06
C-----------------------------------------------------------------------

      subroutine schmidtorthor_lwl(psi,gdim,dim,ierr,mnorm)

      implicit none

      integer    gdim,dim,e,e1,ierr
      real*8     psi(gdim,dim),overlap
      real*8     norm, mnorm

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the set of functions.
C The orthonormalization is made twice to remove numerical inaccuracies.
C-----------------------------------------------------------------------
      mnorm = 1.d9
      ierr = 0
      do e=1,dim
         do e1=1,e-1
            call vvtxdd(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxddr(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxd(psi(1,e),norm,gdim)
         mnorm = min(norm,mnorm)
         if (norm .le. 1.0d-99 ) norm = 1.0d99
         call xvixddo(norm,psi(1,e),gdim)
      enddo

      do e=1,dim
         do e1=1,e-1
            call vvtxdd(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxddr(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxd(psi(1,e),norm,gdim)
         if( norm .le. 0.8d0 ) then
            ierr = e
            return
         end if
         call xvixddo(norm,psi(1,e),gdim)
      enddo

      return
      end


C#######################################################################
