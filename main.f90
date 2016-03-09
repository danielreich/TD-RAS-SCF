!===============================================================================
! main drive of tdrasscf calcuations of atomic and molecular
!===============================================================================

program tdrasscf_3d
 use operator_3d
 implicit none

!
! open the input file, print the log information
!
 call log_information
!
! read in the argument 
! 
call readin
!
! information about the basis set
!
 call drive_fedvr3d_basis_set

!
! information about the matrix element in basis set
!
 call drive_not_full_operator_index

!
! set the operator elements in the fedvr
!
 call drive_operator_radial
!
!calculate the kinetic energy operator, nuclear attractive potential operator
!

 call drive_operator_3d

!
!calculate two-electron integrals in spherical harmonic coordinate 
!
 call drive_twoe_basis_set

!
! auxiliary array 
!

 call drive_auxiliary

!
! initialization of wf
!

 call drive_wfunction

!
! allocate the necessary address for the work array
!
 call configure_operator_spatial_space

 call configure_density

!
! propagate the orbitals and the cofficients
!

 call prop_a_phi()

stop 'God bless you'
end program tdrasscf_3d









