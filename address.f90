!=======================================================================
! find the index of the path by using alphastring and betastring method
! the details can be found in the paper, j.chem.phys. 89(4),1988
! Jeppe Olsen etc.
!-----------------------------------------------------------------------
! input, ibox(nel), nel is the number of the electron,
!        ibox store the index of the spin orbitals
! output: ii_add is the index of the path
!    0     1    2
! 0  (1)
!     |  \
!     |   \           
! 1  (1)  (1)
!     | \  | \ 
!     |  \ |  \ 
! 2  (1)  (2)  (1)
!       \  | \  |
!        \ |  \ |
! 3       (3)  (3)   
!           \   |
!            \  |
! 4           (6)  
!
! the w_vertex is the weight of the vertex, please read the details in 
! the paper, the above example: ne = 2, nptot = 4
!=======================================================================
subroutine address(ii_add,ibox)
!=======================================================================
  use global
  implicit none
  integer :: ibox(system%numelectron),ibox_2nd(system%nptot*2)
  integer :: ii,mm,km,ii_a,ii_b,ii_add,ie
!
! ibox_2nd is the second quantum representation of the slater ibox
! for example, ibox(4) = /1,2,4,5/
! then, ibox_2nd = /1,1,0,1,1/
!
  do mm=1,system%nptot*2
     ii=0
     do ie=1,system%numelectron
        if (mm.eq.ibox(ie)) ii=1
     enddo
     ibox_2nd(mm)=ii
  enddo
!
! find the index of alphastring path
!
  ii_a=0
  do mm=1,system%nptot
     km=0
     do ii=1,mm
        km=km+ibox_2nd(ii)
     enddo
     ii_a=ii_a+ibox_2nd(mm)*W_VERTEX(mm-1,km)
  enddo
!
! find the index of betastring path
!
  ii_b=0
  do mm=1,system%nptot
     km=0
     do ii=1,mm
        km=km+ibox_2nd(ii+system%nptot)
     enddo
     ii_b=ii_b+ibox_2nd(mm+system%nptot)*W_VERTEX(mm-1,km)
  enddo
!
! find the total slater index
!
  ii_add=1+ii_a+ii_b*N_STRING_A_FCI
  return
end subroutine address
