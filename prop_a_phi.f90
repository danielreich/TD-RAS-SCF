!
! solve the equations of acofficient and spatial orbitals
!
 subroutine prop_a_phi()
  use global
  use wfunction
  use analysis_part
  use solver
  use fedvr3d_basis_set
  use operator_3d
  use density

  implicit none
  integer :: iab,iorb,idvr,ierr,icount,element_basis
  real(kind=k1) :: time,mnorm
  complex(kind=k2) :: norm
  real*8 :: xd,yd,zd !! components of the dipole
!! storing of the components of the dipole
  complex*16, allocatable :: xd_store(:), yd_store(:), zd_store(:)
  real*8 :: xdtt,ydtt,zdtt !! components of the dipole acceleratio
!! storing of the components of the dipole acceleration
  complex*16, allocatable :: xdtt_store(:), ydtt_store(:), zdtt_store(:)
  integer :: time_step_n,ijk
  integer :: idicp, jdicp
 integer, allocatable :: global_to_local_rho(:,:)
  complex*16, allocatable :: somega(:)
  real*8, allocatable :: omega(:)
  real*8 :: tstep_fourier
 real(kind=k1), allocatable :: fedvr_k_global(:) !! fedvr node in the momentum space.
 real(kind=k1), allocatable :: fedvr_k_weight(:,:), fedvr_k_node(:,:) !! weights and nodes in the element in the first argument and for the point in the element, given by the second argument.
 complex(kind=k2), allocatable :: phik(:,:)

complex(kind=k2), allocatable :: rhor_phi(:,:) !! monoparticular 3D density for each orbitals
!! Density of a single orbital.
 complex(kind=k2), allocatable :: rhor_all(:) !! 3D density for all the electrons
!! Density the whole wavefunction.
Character*6:: rhorall
Character*17:: rhorphi
 complex(kind=k2), allocatable :: rho_phi(:,:) !! monoparticular 3D density for each orbitals
!! Density of a single orbital.
 complex(kind=k2), allocatable :: rho_all(:) !! 3D density for all the electrons
!! Density the whole wavefunction.

!! Auxiliar variables
Character*6:: rhoall
Character*17:: rhophi


!
! configure the array in solver0
!
  call configure_solver
!
! count the time of cycle
!
!  icount = 0 ! the time of cycle, the time of imaginary propgation


  icount = 0
  time = prop%t0
!!$  do while(icount<=prop%icount_max) !! Wenliang
!!$    icount = icount +1 !! Wenliang
  Do icount=1,prop%icount_max

    if(laser%tdornot) then
       if(time>prop%t_end) then
          print*, 'Exit of the time propagtion, because the final time is reached'
        exit
      endif 
    endif
 !
 !  runge-kutta method, the begining wave function is acoff0 ,phi0 
 !

    if(prop%solver_which== 0 ) then    

       call solver0(time,icount)

    else
       write(*,*) ' ----new solver method---not yet implement-- '
       stop 
    endif    

    time = time + prop%hstep

!
! now we get a new wave function phi and acoff
!

    if(.not.laser%tdornot) then
!===========================================================================
! normalize  acoff
!
     norm = zzero

     do iab =1,n_string_ab
      norm = norm + real(conjg(acoff(iab))*acoff(iab))
     enddo
      
     do iab =1,n_string_ab
      acoff(iab) = acoff(iab)/sqrt(real(norm))
     enddo 
     
!==========================================================================
! JUAN's CHANGE: THE ORBITALS IN NP0 CLOSE SHELL

     call np0_lm(phi,system%np0,fedvr3d%l_max,fedvr3d%m_max)
     print*, 'setting core to mean_field_approximation'
     

 !==========================================================================
 !orth-normalize the spatial orbital
 !
    call schmidtortho2(phi,fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot,ierr,mnorm)
    if(ierr/=0) then
      write(*,*) '-----orth-normalize wrong in tdrasscf-----'
      stop 
   endif
endif

!
! prepare acoff0 and phi0, for next cycle   
! 
     do iab=1,n_string_ab 
      acoff0(iab) = acoff(iab)
     enddo
!
! prepare phi0, for the next cycle
!
    do iorb = 1, system%nptot
       do idvr= 1,fedvr3d%nb_r*fedvr3d%nb_angle
        phi0(idvr,iorb) = phi(idvr,iorb) 
       enddo
    enddo

!==============================================================================================
! analysis part
!==============================================================================================

!!$    If (icount.le.prop%icount_max_pulse) then


!! Block IF done by Juan: Ensures that only prop%icount_store points are saved
!!$    If (mod(icount,prop%icount_max_pulse/prop%icount_store).eq.0.or.icount.eq.prop%icount_max_pulse) then 

    If (laser%tdornot) then !! Analysis for real time propagation
       
       If (mod(icount,prop%icount_max_pulse/prop%icount_store).eq.0.and.icount.le.prop%icount_max_pulse) then !! Storage during the pulse
          
          call drive_analysis_part(time,xd,yd,zd,xdtt,ydtt,zdtt,icount)
          
       End If       
       
       If (mod(icount,prop%icount_max/prop%icount_ionization).eq.0.and.icount.ge.prop%icount_max_pulse) then !! Storage after the pulse
          
          call drive_analysis_part(time,xd,yd,zd,xdtt,ydtt,zdtt,icount)
          
       End If
       
    End If


    !! Store the dipole moment to perform the calculation of the HHG

    if (analysis%i_hhg==1.and.laser%tdornot.and.icount.le.prop%icount_max_pulse.and.mod(icount*prop%icount_store,prop%icount_max_pulse).eq.0) then
          
          !initialize the dipole store
       If (.not.allocated(xd_store)) then
          allocate(xd_store(1:prop%icount_store))
          allocate(yd_store(1:prop%icount_store))
          allocate(zd_store(1:prop%icount_store))
          xd_store=zzero
          yd_store=zzero
          zd_store=zzero
       End If
       
       !initialize the dipole acceleration store
       
       If (.not.allocated(xdtt_store)) then
          allocate(xdtt_store(1:prop%icount_store))
          allocate(ydtt_store(1:prop%icount_store))
          allocate(zdtt_store(1:prop%icount_store))
          xdtt_store=zzero
          ydtt_store=zzero
          zdtt_store=zzero
       End If
       
       time_step_n=icount*prop%icount_store/prop%icount_max_pulse

       !! store the dipole 
       print*, 'hola', time_step_n, prop%icount_max_pulse
       xd_store(time_step_n)=cmplx(xd,0.d0)
       yd_store(time_step_n)=cmplx(yd,0.d0)
       zd_store(time_step_n)=cmplx(zd,0.d0)
       
       !! store the dipole acceleration
       
       xdtt_store(time_step_n)=cmplx(xdtt,0.d0)
       ydtt_store(time_step_n)=cmplx(ydtt,0.d0)
       zdtt_store(time_step_n)=cmplx(zdtt,0.d0)
       
    End if
    
!!$    End If
!!$ End If

!! Store during imaginary time propagation

    if(.not.laser%tdornot.and.mod(icount,prop%icount_max/prop%icount_store).eq.0.and.icount.le.prop%icount_max) then
!
! store the energy, dipole moment and dipole acceleration
!
       call drive_analysis_part(time,xd,yd,zd,xdtt,ydtt,zdtt,icount)

       write(*,*)
       write(*,*) 'Storing the orbitals and coefficients'
       write(*,*)

!store the coefficients of each configuration to read by user

       Open(217,file='coefficients.dat')

       do iab = 1,n_string_ab
          write(217,*) 'Coefficient', acoff(iab)
          write(217,*) 'Orbitals occupied:', i_string_ab(iab,:)
          write(217,*)
       enddo
       
       Close(217)

!! save the output used as input for the real time propagation.       

       Open(217,file='wf0.txt') 

       write(217,*) '#INPUT FILE FOR THE REAL PROPAGATION'
       write(217,*) '#==================================='
       write(217,*) '#'
       write(217,*) '# Nuclear charge (Z), Number of electrons'
       write(217,*) '#'
       write(217,'(F7.3,I3.2)') system%zz, system%numelectron
       write(217,*) '#'
       write(217,*) '# RAS SCHEME: NP2, NP1, NP0 AND SDT'
       write(217,*) '#'
       write(217,'(3I3.2,I2.1)') system%np2, system%np1, system%np0, system%io_sdtq
       write(217,*) '#'
       write(217,*) '# r_begin, r_end inner region, number of inner element, r_end, number of elements, nodes per element'
       write(217,*) '#'
       write(217,'(F6.3,F7.3,I4.3,F7.3,I4.3,I3.2)')  fedvr3d%r_begin,fedvr3d%r_inner,fedvr3d%number_of_element_inner,fedvr3d%r_end, fedvr3d%number_of_element,fedvr3d%fedvr_nb(1)
       write(217,*) '#'
       write(217,*) '# l_max, m_max'
       write(217,*) '#'
       write(217,'(I3.2,I3.2)') fedvr3d%l_max, fedvr3d%m_max
       write(217,*) '#'
       write(217,*) '# COEFFICIENTS'
       write(217,*) '#'
       Do iab=1,n_string_ab !! Running in the Slater Determinants   
          write(217,*) acoff(iab) 
       End Do
       write(217,*) '#'
       write(217,*) '# ORBITALS'
       write(217,*) '#'

       Do iorb = 1, system%nptot  !! Running in the orbitals

          write(217,*) 'Orbital ',iorb
          write(217,*) '---------'
          write(217,*) '#'

          do idvr= 1,fedvr3d%nb_r*fedvr3d%nb_angle !! Running in the spatial basis
             write(217,*) phi(idvr,iorb)

          enddo
          
          write(217,*) '#'
          
       End Do

       write(217,*) '#'
       write(217,*) '#'
       write(217,*) '#'

       Close(217)

!! store the coefficients together with the orbitals, to read by the 
!! code.
       
        open(ifile,file='imag.txt') 
        open(217, file='orbitals.txt')

       do iab = 1,n_string_ab !! Running in the Slater Determinants
          write(ifile,*) acoff(iab)
       enddo
       
       write(ifile,*)
       write(ifile,*)             

       do iorb = 1, system%nptot  !! Runing in the orbitals
          write(217,*)
          write(217,'(a4,i2.2)') 'phi_',iorb
          write(217,'(a6)') '======'
          write(217,*) ''
          write(217,*) 'Element           Basis           l          m    phi'
          write(217,*) ''
          do idvr= 1,fedvr3d%nb_r*fedvr3d%nb_angle !! Runing in the spatial basis
             element_basis=global_to_local(idvr,1)
             write(ifile,*) phi(idvr,iorb)

             write(217,*) which_element(element_basis),which_basis(element_basis),global_to_local(idvr,2), global_to_local(idvr,3), phi(idvr,iorb)

          enddo
          write(ifile,*)
          write(ifile,*)
       enddo
 
       close(ifile)
       close(217)

       !! Finish imaginary time propagation if dipole moment and the dipole acceleration are zero.

       If (abs(xd).lt.1.d-5.and.abs(yd).lt.1.d-5.and.abs(zd).lt.1.d-5.and.abs(xdtt).lt.1.d-5.and.abs(ydtt).lt.1.d-5.and.abs(zdtt).lt.1.d-5) then
          close(101) !!close dipole file
          close(108) !! close dipole acceleration file
         write(*,*) 
         write(*,*) 'Dipole zero=> STOP'
         write(*,*) 
         exit !! Finish the propagation
       End If

    End if


 enddo


 if(.not.laser%tdornot) then
!
! store the final ground energy and wf
!

!coefficients of each configuration

Open(217,file='coefficients.dat')

  do iab = 1,n_string_ab
   write(217,*) 'Coefficient', acoff(iab)
   write(217,*) 'Orbitals occupied:', i_string_ab(iab,:)
   write(217,*)
  enddo

Close(217)

 open(ifile,file='imag.txt') !!! Changed by Juan

!storing coefficient and wf to be input of the next cycle
   do iab = 1,n_string_ab
      write(ifile,*) acoff(iab)
   enddo

   write(ifile,*)
   write(ifile,*)

  do iorb = 1, system%nptot
    do idvr= 1,fedvr3d%nb_r*fedvr3d%nb_angle 
      write(ifile,*) phi(idvr,iorb)
    enddo
    write(ifile,*)
    write(ifile,*)
  enddo

  close(ifile)
 
Else !! If there is a time propagation

!coefficients of each configuration

Open(217,file='coefficients_real.dat')

  do iab = 1,n_string_ab
   write(217,*) 'Coefficient', acoff(iab)
   write(217,*) 'Orbitals occupied:', i_string_ab(iab,:)
   write(217,*)
  enddo

Close(217)

Open(1120,file='real.txt')

!storing coefficient and wf to be input of the next cycle
   do iab = 1,n_string_ab
     write(1120,*) acoff(iab)
   enddo

   write(1120,*)
   write(1120,*)

  do iorb = 1, system%nptot
    do idvr= 1,fedvr3d%nb_r*fedvr3d%nb_angle 
      write(1120,*) phi(idvr,iorb)
    enddo
    write(1120,*)
    write(1120,*)
  enddo

  close(1120)   

if(analysis%i_hhg==1) then !! calculation of HHG
             tstep_fourier=prop%icount_store*prop%hstep

!! Calculation of the FFT to obtain the HHG
!! <Z>
             Call fastft(zd_store,tstep_fourier,somega,omega) !!calculate Fast Fourier Transform

             Do iab=1,size(somega)
                write(file4,*) omega(iab), real(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0)),aimag(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0))
             End Do

close(file4)

!!<Y>

open(file4,file='hhg_y.txt')

             Call fastft(yd_store,tstep_fourier,somega,omega) !!calculate Fast Fourier Transform

             Do iab=1,size(somega)
                write(file4,*) omega(iab), real(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0)), aimag(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0))
             End Do

close(file4)

!!<X>

open(file4,file='hhg_x.txt')

             Call fastft(xd_store,tstep_fourier,somega,omega) !!calculate Fast Fourier Transform

             Do iab=1,size(somega)
                write(file4,*) omega(iab), real(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0)), aimag(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0))
             End Do

close(file4)

!!<\partial_{tt} Z>

open(file4,file='hhg_z_tt.txt')

             Call fastft(zdtt_store,tstep_fourier,somega,omega) !!calculate Fast Fourier Transform

             Do iab=1,size(somega)
                write(file4,*) omega(iab), real(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0)),aimag(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0))
             End Do

close(file4)

!!<\partial_{tt} Y>

open(file4,file='hhg_y_tt.txt')

             Call fastft(ydtt_store,tstep_fourier,somega,omega) !!calculate Fast Fourier Transform

             Do iab=1,size(somega)
                write(file4,*) omega(iab), real(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0)),aimag(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0))
             End Do

close(file4)

!!<\partial_{tt} X>

open(file4,file='hhg_x_tt.txt')

             Call fastft(xdtt_store,tstep_fourier,somega,omega) !!calculate Fast Fourier Transform

             Do iab=1,size(somega)
                write(file4,*) omega(iab), real(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0)),aimag(dconjg(somega(iab))*somega(iab)*(omega(iab)**2.0d0))
             End Do


          End if

       endif

!!PLOT DENSITY (for the last time step) for REAL AND IMAGINARY TIME PROPAGATION

   print*, 
   print*, 'P(r) at the end of the propagation'
   print*, 

!! Update the rho1 and rho2 to calculate the one and two body density functions

   call rho1()

!! Radial one body density function.

   Call plot_phi_fedvr3d(phi,fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_w, cden1)

   Call rho2()
   
!! Radial two body density function.

   print*, 
   print*, 'P(r_1,r_2) at the end of the propagation'
   print*, 'Using the coupled basis'
   print*, 

   call plot_phi2_fedvr3d_coupled(phi,fedvrx_global,which_element,which_basis, fedvr_w,cden2)

!!$print*, 'prop_a_phi.f90'
!!$stop
!!$
!!$   Call plot_phi2_fedvr3d(phi,fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_w, cden2)

!! Calculation of the spatial basis set for the density in 3D.

print*, 'prop_a_phi.f90'
print*, 'call basis_set_rho3d'

   call basis_set_rho3d(which_element,which_basis,fedvr3d%l_max,fedvr3d%m_max,global_to_local_rho)

!! Calculation of the 3D monoparticular density.
!! rhox_phi is the density of each orbital.
!! rhox_all is the total density.

   print*, 
   print*, 'Calculation of \rho(\vec{r})'
   print*, 'Using the coupled basis'
   print*, 

   call rho3d_coupled(phi, cden1,rhor_phi, rhor_all)   

!!$   print*, 
!!$   print*, 'Calculation of \rho(\vec{r})'
!!$   print*, 'Not using the coupled basis'
!!$   print*, 
!!$
!!$call rho3d(phi,fedvrx_global,global_to_local, global_to_local_rho,cden1,rhor_phi, rhor_all) !! CONSTRUCTION OF THE TOTAL DENSITY

!! Write the results in files

If (allocated(rhor_phi)) then !! DEBUG 04.03.2015

Do jdicp=1,size(rhor_phi(1,:))
   write(rhorphi,'(a15,i2.2)') 'rhor3d_orbital_',jdicp
   Open(197+jdicp,file=rhorphi)
      write(197+jdicp,*) '# Time= ', time
      write(197+jdicp,*) ''
   Do idicp=1,size(rhor_phi(:,1))
      write(197+jdicp,*) idicp,rhor_phi(idicp,jdicp)
   End Do
   close(197+jdicp)
End Do

write(rhorall,'(a6)') 'rhor3d'

Open(198,file=rhorall)

   Do idicp=1,size(rhor_phi(:,1))
      write(198,*) idicp,rhor_all(idicp)
   End Do

close(198)

End If !! DEBUG 04.03.2015
stop

   xd=0.d0
   yd=0.d0
   zd=0.d0

 return
end subroutine prop_a_phi
