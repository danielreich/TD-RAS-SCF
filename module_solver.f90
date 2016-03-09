!
! module solver, to solver the equation of A and phi
!
 module solver
   use global
   use auxiliary
   use wfunction

   implicit none
!   public 
   complex(kind=k2),allocatable,save :: ka(:),kphi(:,:)
   complex(kind=k2),allocatable,save :: suma0(:)
   complex(kind=k2),allocatable,save :: sumphi0(:,:)
   contains
!=========================runge-kutta 4 order ======================================
   subroutine solver0(t,i)
    implicit none
    real(kind=k1),intent(in) :: t
    integer, intent(in) :: i
    integer(4) :: idicp,jdicp,itime
! 
!==========================k1========================================================
!   

   call derivate(t,i)
! acoff    

   do idicp =1,n_string_ab
     acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(0.5_k1)*prop%chstep  
     suma0(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(1.0_k1/6.0_k1)*prop%chstep
   enddo

! phi

   do jdicp =1,system%nptot
    do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle

     phi(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(0.5_k1)*prop%chstep
     sumphi0(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(1.0_k1/6.0_k1)*prop%chstep
    enddo
  enddo

! 
!==========================k2========================================================
! 

  call  derivate(t+0.5_k1*prop%hstep,i)

! acoff
  do idicp =1,n_string_ab
     acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(0.5_k1*prop%chstep)  
     suma0(idicp) = suma0(idicp) + ka(idicp)*dcmplx(prop%chstep/3.0_k1)
  enddo




! phi
  do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
    do jdicp =1,system%nptot

     phi(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(0.5_k1*prop%chstep)
     sumphi0(idicp,jdicp) = sumphi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(prop%chstep/3.0_k1)
   enddo
  enddo
! 
!==========================k3========================================================
!
!!$  print*, ka
!!$     print*, acoff
!!$     print*, 'module_solver.f90'
!!$     stop

  call   derivate(t+0.5_k1*prop%hstep,i)

! acoff
  do idicp =1,n_string_ab
     acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(prop%chstep)  
     suma0(idicp) = suma0(idicp) + ka(idicp)*dcmplx(prop%chstep/3.0_k1)
  enddo
! phi
  do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
    do jdicp =1,system%nptot

     phi(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(prop%chstep)
     sumphi0(idicp,jdicp) = sumphi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(prop%chstep/3.0_k1)
    enddo
  enddo

! 
!==========================k4========================================================
!
  call  derivate(t+prop%hstep,i)
! acoff
  do idicp =1,n_string_ab  
     acoff(idicp) = suma0(idicp) + ka(idicp)*dcmplx(prop%chstep/6.0_k1)
  enddo
! phi
  do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
    do jdicp =1,system%nptot
     phi(idicp,jdicp) = sumphi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(prop%chstep/6.0_k1)
    enddo
  enddo

 return
end subroutine solver0

end module solver


subroutine configure_solver
  use global
  use auxiliary
  use solver
  implicit none

   allocate(ka(n_string_ab)) 
   allocate(kphi(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot))
   allocate(suma0(n_string_ab))
   allocate(sumphi0(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot))
  
  return
end subroutine configure_solver




       
     












 
