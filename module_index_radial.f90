 
!
! the spare operator are using in the calcuations
!
 module not_full_operator_index
  use global
  use fedvr3d_basis_set
  implicit none
  integer :: n_total_kinetic  ! the number of kinetic element with non-zero 
  integer,allocatable,save :: index_kinetic_basis(:,:)
  integer,allocatable,save :: index_ham_act(:,:)
  integer,allocatable,save :: index_kinetic_operator(:,:)
  contains

 
  subroutine index_of_kinetic
   implicit none
   integer :: idicp,isum,ibasis,jbasis,e,f,itemp,jtemp
   integer :: i,j
!---------------------------------------------------------------
!                radial part
! figure out how many element in kinetic operator are non-zero
!---------------------------------------------------------------
   isum = 0
   isum = (fedvr3d%fedvr_nb(1) -1)**2
   
   do idicp =2, fedvr3d%number_of_element -1 
     isum = isum + (fedvr3d%fedvr_nb(idicp))**2 
   enddo

    isum = isum + (fedvr3d%fedvr_nb(fedvr3d%number_of_element) -1)**2
    
    n_total_kinetic = isum - (fedvr3d%number_of_element - 1)

    n_total_kinetic = (n_total_kinetic - fedvr3d%nb_r)/2 + fedvr3d%nb_r 
!
! there are total n_total_kinetic non-zero elements in radial part kinetic operator matrix element
!
   allocate( index_kinetic_basis(n_total_kinetic,2)  )
   allocate( index_ham_act(fedvr3d%nb_r,2))
   allocate(index_kinetic_operator(fedvr3d%nb_r,fedvr3d%nb_r))

    itemp = 0
    do ibasis = 1, fedvr3d%nb_r 
       
     ! jtemp = 0
      do jbasis =1,ibasis
 
       e = which_element(ibasis) ; i = which_basis(ibasis)
       f = which_element(jbasis) ; j = which_basis(jbasis)


      if(e==f .or. ((e==(f-1)).and. (i==fedvr3d%fedvr_nb(e))) .or.((e==(f+1)).and.(j==fedvr3d%fedvr_nb(f)))) then
        itemp = itemp +1
        index_kinetic_basis(itemp,1) = ibasis
        index_kinetic_basis(itemp,2) = jbasis 

        index_kinetic_operator(ibasis,jbasis) = itemp


      !  jtemp = jtemp +1
      !  if(jtemp ==1) then
      !      index_ham_act(ibasis,1) = jbasis
      !  endif 
      !  index_ham_act(ibasis,2) = jbasis
       endif
      enddo
   enddo

!
!=================================================================================
!
    do ibasis = 1, fedvr3d%nb_r 
       
      jtemp = 0
      do jbasis =1,fedvr3d%nb_r
 
       e = which_element(ibasis) ; i = which_basis(ibasis)
       f = which_element(jbasis) ; j = which_basis(jbasis)
      if(e==f .or. ((e==(f-1)).and. (i==fedvr3d%fedvr_nb(e))) .or.((e==(f+1)).and.(j==fedvr3d%fedvr_nb(f)))) then
        jtemp = jtemp +1
if(jtemp ==1) then
  index_ham_act(ibasis,1) = jbasis
endif 
          index_ham_act(ibasis,2) = jbasis
      endif
      enddo
   enddo
 




    if(itemp /= n_total_kinetic) then
       write(*,*) 'error happen in set_kinetic_energy_operator'
       stop 
    endif

   return
  end subroutine index_of_kinetic








 end module not_full_operator_index

!
! the full
!
subroutine drive_not_full_operator_index
 use not_full_operator_index
 implicit none

 call index_of_kinetic

return
end subroutine drive_not_full_operator_index











