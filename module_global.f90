!-----------------------------------------------------------
! module global parameter used in the calculations
!-----------------------------------------------------------
 module global
      implicit none      
      integer,parameter:: in_channel =10,ifile=1010
      integer :: file1,file2,file3,file4,file5,file6,file7
      integer,parameter :: k1 =KIND(0.0d0)
      integer,parameter :: k2 =KIND((0.0d0,0.0d0))
      integer,parameter :: ia_max = 6000,max_element = 1000,max_basis = 10000
      real(kind=K1),parameter :: PI=3.1415926535897932d0,&
      zero=0.0d0,eps_Q = 1.0d-10,speedc = 2.99792458d0
      integer,parameter :: izero = 0 
      complex(kind=K2),parameter :: zzero =(0.0d0,0.0d0),zone=(1.0d0,0.0d0),ci=(0.0d0,1.0d0) 
      integer :: ia(100000)
! system 
      type system_prime
        character(len=70) :: title ! the calculation title
        integer :: numelectron
        real(kind=k1) :: zz  ! the number of the electron, the charge of the nulcei
        integer :: nptot           ! the total number spatial number used in the calcuations
        integer :: np2,np1,np0     ! the number of the spatial orbitals in different space np2,np1,np0
        integer :: io_sdtq         ! =0, D, =1 S, =2 SD
        integer :: io_fcore        ! =0, no core orbital, /=0, core orbital are fixed
        integer :: nptot_fc,n_offset,io_ir
        integer :: io_orb_eqs_p
     end type system_prime

!fedvr3d
      type fedvr3d_prime
       real(kind=k1) :: r_begin,r_end   ! the boundary of the radial part {r_begin         r_end}
       real(kind=k1) :: r_inner ! boundary of the inner region in the fedvr. In this region there are number_of_element_inner elements
       integer :: number_of_element     ! the number of the element
       integer :: number_of_element_inner ! the number of elements in the inner region
       integer :: fedvr_nb(max_element) ! the nodes in every element 
       integer :: l_max,m_max           ! the maximum l,m in the calculation
       integer :: l_max_core            ! the maximum l for the core
       integer :: m_max_core ! the maximum m for the core
       integer :: l_max_coupled !! maximum L allowed in the coupled basis
       integer :: nb_r,nb_angle         ! the number of the basis function in radial part and angle part respectively 
       integer :: store !! if 0, cww is not stored
                        !! if 1, cww is stored
       real(kind=k1) :: r_out !! r_out defines the outer region.
       real(kind=k1) :: delta !! delta is the parameter of the hamming window, used to obtain the photoelectron spectrum
       real(kind=k1) :: k_max,k_min!! maximum and minimum linear momentum
    end type fedvr3d_prime

!! Nonzero matrix elements
    
    type h_nonzero_elements
       real(kind=k1), allocatable :: value(:) !! value of the Hamiltonian one body operator. value(:)=<bra(:)|h|ket(:)>
       integer, allocatable :: bra(:), ket(:),l(:) !! ket and bra, as it is said in the previous line.
    End type h_nonzero_elements

! condition on the angular momentum to calculate the mean field operator
! It is the states which fulfill the condition m4-m3=m1-m2

    type m4m3
       integer, allocatable :: ang3(:,:), ang4(:,:)
       integer, allocatable :: total(:) !! total number of elements for each m4-m3
    end type m4m3



!! matrix to transform from uncoupled to the coupled basis in angular momentum.

    type llmm_lm_type
       !! alpha and beta correspond to the uncoupled basis
       integer,allocatable :: alpha(:) !! is the spherical harmonic conjg(alpha)
       integer, allocatable :: beta(:) !! is the spherical harmonic beta
       !! gamma corresponds to the coupled basis
       integer, allocatable :: gamma(:) !! is the spherical harmonic gamma
       !! is the coefficient which links them.
       real (kind=k1), allocatable :: value(:)
    end type llmm_lm_type

    complex(kind=k2), allocatable :: phi_coupled(:,:,:)

!laser field
      type laser_prime
	logical :: tdornot                              ! realtime or not
        integer :: ncycle                               ! how many cycle are used
	real(kind=k1):: ex_max,ey_max,ez_max,e_max,ellipticity !! field indensity and ellipicity in the calcuations 
        real(kind=k1):: omega,lamda_nm                  ! frequency and wavelength in the calcuations
        real(kind=k1):: ex_t,ey_t,ez_t                  ! time dependent laser field
        real(kind=k1):: ax_t,ay_t,az_t                  ! time dependent laser field
        real(kind=k1) :: cep !! Carrier phase envelope
        character*1 :: gauge !! Gauge: 'l-> length' and 'v-> velocity' gauge
     end type laser_prime

! propagator information
      type prop_prime               
	real(kind=k1) :: t0,t_end,hstep  ! the beginning time,end time,and the propagation step
        real(kind=k1) :: ncycle_extra     ! cycles after to continue the propagation
        integer :: solver_which          ! which solver method are used in the calcuations
        integer :: icount_max            ! the maximum loop in the imaginary time propgations
        integer :: icount_max_pulse      ! the maximum number of time steps during the pulse
        integer :: icount_store          ! number of points to store for HHG
        integer :: icount_ionization     ! number of points to store the ionization
        complex(kind=k2):: chstep        ! propagation step
        integer :: stiffness             ! =1, consider the stiffness effect, =0, not considered
     end type prop_prime

     type absorb_prime
      real(kind=k1) :: r_left, r_right  ! the absorb potential, the aboundary 
      real(kind=k1) :: ramp             ! the argument of absorb potential
      real(kind=k1) :: strength         ! the argument of absorb potential
   end type absorb_prime

     type analysis_prime
      integer :: i_energy              ! energy are calculated 
      integer :: i_auto                ! auto function are calcualted 
      integer :: i_ionization          ! ionization probability are caculated
      integer :: i_hhg                 ! hhg are calcuated in the calculation
      integer :: i_twoe_density        ! calc. the two-e density
      integer :: i_momentum            ! momentum spectrum
      integer :: i_atas                ! attosecond transit absorb spectrum  
   end type analysis_prime
 
      type(analysis_prime)  :: analysis
      type(absorb_prime) :: absorb
      type(system_prime) :: system
      type(laser_prime)  :: laser
      type(prop_prime)   :: prop
      type(fedvr3d_prime) :: fedvr3d


      type(llmm_lm_type) :: llmm_lm !! transform from uncoupled to couple, that is
!! Integrals of Y_{alpha}* Y_{beta} Y_{gamma}*, where alpha and beta are in the uncoupled basis and gamma is the coupled basis.
      type(llmm_lm_type) :: integral_llmm_lm
      !! Integrals of Y_{alpha}* Y_{beta} Y_{gamma}, where alpha and beta are in the uncoupled basis and gamma is the coupled basis.

      !! nonzero one body Hamiltonian matrix elements

      type(h_nonzero_elements), save :: h_nonzero
    
end module global
