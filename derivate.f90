! ========================================================================
! find the derivate of A and spatial orbital
! ========================================================================
subroutine derivate(time,i)
  use global
  use auxiliary
  use density
  use operator_spatial_space
  use solver
  implicit none
  real(kind=k1),intent(in) :: time
  integer, intent(in) :: i
  integer :: idicp,jdicp,kdicp,ldicp,method

  
!--------------------------------------------------------------
! calc. the rho1 and rho2
!--------------------------------------------------------------
     call rho1()
     call rho2()

!    if(system%io_fcore/=1) then

      call inverse_rho1()

      call rho3()
!    else
!      call invers_rho1_fc()
!    endif


!--------------------------------------------------------------
! calc. the matrix expression of ham. and v(r1,r2) in time-depe
! ndent spatial orbitals
! one electron operator and two electron operator
!--------------------------------------------------------------   


   if(laser%tdornot) then !! Real time propagation
      if (laser%gauge.eq.'l') then !! Length gauge
         call laserpara(time)
      Else if (laser%gauge.eq.'v') then !! Velocity gauge
         call laserpara_a(time)
         Else 
            write(*,*) ''
            write(*,*) '================================'
            write(*,*) 'Select l or v for length or'
            write(*,*) 'velocity gauge, respectively'
            write(*,*) '================================'
            write(*,*) ''
            STOP
      End if
   endif

    call update_ham_onebody()
    
!! If we store the matrix elements of the two body interactions 
!! fedvr3d%store=1-> Using uncoupled basis
!! fedvr3d%store=2-> Using coupled representation


    If (fedvr3d%store.eq.1.or.fedvr3d%store.eq.2) call update_vv_twobody() 

!! commented by Wenliang
!
!        if(laser%tdornot) then
!         call update_zmat_spatial_space()
!        endif
!


    do idicp =1,system%nptot
       do jdicp =1,idicp          
          ch_dummy(idicp,jdicp) = ham_onebody(ia(idicp)+jdicp)
          if(idicp/=jdicp) then  
             ch_dummy(jdicp,idicp) = dconjg(ch_dummy(idicp,jdicp))
          endif
       enddo
    enddo


!--------------------------------------------------------------
! two different way to solve Q_space equations
! Q_space equations
!--------------------------------------------------------------
!    if(system%io_fcore/=1) then

    
    select case (fedvr3d%store)
    case(0)
!!$       call Qspace_opt()    ! Juan: I did this to optimize
!!$       call Qspace_omp()    ! Juan: I did this to parallelize
       call Qspace_omp_select()    ! Juan: I did this to optimize and parallelize
    case(1)
       call Qspace()    ! usually
    case(2)
       !!  
       !! This subroutine is using the coupled representation technique.
       call Qspace_coupled()
    end select

!    else
!      call Qspace_fc()   ! hf fcore
!    endif


   if((system%io_orb_eqs_p/=0).and.(system%io_fcore/=1)) then



!-----------------------------------------
       if(system%io_sdtq/=0) then
          call zet_mat1
          call zet_mat2
          call orb_eqs_p_zet
       endif
!-----------------------------------------
     ! call Pspace_eqs(system%NP1,system%NP2,system%NP0,system%NP0+system%NP1) ! ch_dummy

       call orb_eqs_p()

       call acoff_eqs()  

       call combin_p_q()

    else

!--------------------------------------------------------------
! find derivate of coff. A
!--------------------------------------------------------------

    call acoff_eqs()   

 endif

return


end subroutine derivate

   
