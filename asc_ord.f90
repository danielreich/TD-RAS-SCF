!============================================================================================================
! input: ibox(1:numelectron) is array, store the spin orbital, (A slater Derter.)
! ouput: kk: if kk ==0, this slater are not exist, or it is not a real slater(not satifsy the pauli principle)
!            if kk /=0, kk is the index of the slater, 
!            if kk /=0, ipm store the sign, + / -  
!============================================================================================================
subroutine asc_order(ibox,i_sign,i_index)

!=================================================================
  use global
  
  implicit none
  integer,intent(out) :: i_sign,i_index
  integer  :: ibox(system%numelectron)
  integer :: judge0,ie,iep,ii_add,nn
!
! if two electron occupied the same orbital, not real slater, return,kk=0 
!    
  judge0=0
  do ie=2,system%numelectron
    do iep=1,ie-1
      if (ibox(ie)==ibox(iep)) judge0=1
    enddo
  enddo

!----------------------------------------------------------------------------
! if it is a real slater, then we have to try find its position (index)
!
  if (judge0.eq.0) then !-----------------------------------------------------
!----------------------------------------------------------------------------
!
  i_sign=1
!
!exchange the sequence of the orbital, we need ordered the orbitals 
!from small one to large one, the ipm remember the exchange time
! 
  do ie=1,system%numelectron-1
   do iep=ie+1,system%numelectron
     if (ibox(ie).gt.ibox(iep)) then
        nn=ibox(ie)
        ibox(ie)=ibox(iep)
        ibox(iep)=nn
        i_sign=i_sign*(-1)
     endif
   enddo
  enddo
!
!find the index of the slater derter.
!
  call address(ii_add,ibox)  
  i_index = n_string_ab_inv(ii_add)  

!----------------------------------------------------------------------------
  else !---------------------------------------------------------------------
!----------------------------------------------------------------------------
  i_index = 0
!----------------------------------------------------------------------------
  endif !--------------------------------------------------------------------
!----------------------------------------------------------------------------
  return
end subroutine asc_order
 
