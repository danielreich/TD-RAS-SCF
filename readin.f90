!
! read the necessary parameter from the input files
! input the necessary data, need in the calculation
! the details of the input parameter can be found in the 
! code. 
subroutine readin
   use global 
   implicit none
   character(len=70) :: title
   logical :: tdornot
   integer :: nb_sub
   integer :: numelectron,np2,np1,np0,io_sdtq,io_fcore
   integer :: number_of_element,number_of_element_inner,l_max,m_max, l_max_core
   integer :: l_max_coupled
   integer :: ncycle
   integer :: solver_which,icount_max,stiffness, icount_store,icount_ionization
   integer :: i_energy,i_auto,i_ionization,i_hhg,i_twoe_density,i_momentum,i_atas
   !!Name of the files
   character*89 :: energy_file
!========== real type data =======================================================================
   real(kind=k1) :: tcycle_fs, tcycle_au, energy_ev,zz
   real(kind=k1) :: r_begin,r_end 
   real (kind=k1) :: r_inner
   real(kind=k1) :: ex_max,ey_max,ez_max, e_max,ellipticity,omega,lamda_nm,cep
   character*1 :: gauge !! 'l-> length' and 'v-> velocity' gauge
   real(kind=k1) :: t0,t_end,hstep,r_left,r_right,ramp,strength
   real(kind=k1) :: t_end_hhg,ncycle_extra
   real(kind=k1) :: r_out, delta
   real(kind=k1) :: k_max,k_min
   integer :: store, exponent_2

!==================================================================================================   
   namelist /sys/ title,numelectron,zz,np2,np1,np0,io_sdtq,io_fcore
   namelist /basis/r_begin,r_end,number_of_element,number_of_element_inner,r_inner,l_max,m_max,l_max_core,l_max_coupled,nb_sub,store,r_out,k_max,k_min,delta
   namelist /laserfield/ tdornot,ncycle,e_max,ellipticity,omega,lamda_nm,cep,gauge
   namelist /prop_plus_absorb/t0,t_end,ncycle_extra,hstep,solver_which,icount_max,icount_store,icount_ionization,r_left,r_right,&
   ramp,strength,stiffness    
   namelist/ana/ i_energy,i_auto,i_ionization,i_hhg,i_twoe_density,i_momentum,i_atas  

!==================================================================================================
! read the data from the input file
!==================================================================================================     
    read(in_channel,sys,end=100)
100 continue
    read(in_channel,basis,end=101)
101 continue   
    read(in_channel,laserfield,end=102)
102 continue
    read(in_channel,prop_plus_absorb,end=103)
103 continue   
    read(in_channel,ana,end=104)
104 continue

!--------------------------------------------------------------------------------------------------
! input the title and the information of the system
!--------------------------------------------------------------------------------------------------
   system%title = title   ! the titel of the calculation, such as, which system, which method
   system%numelectron = numelectron
   system%zz = zz ! the number of the electron, the charge of the nulcear

   write(*,*) '================================calculating system==================================' 
   write(*,*) system%title
   write(*,*)
   write(*,*) 'the number of electrons is ',system%numelectron
   write(*,*)
   write(*,'(a,1f12.3)') ' the number of nuclei charge is ',system%zz
   write(*,*)
!---------------------------------------------------------------------------------------------------
!  the restricted active space information(ras) 
! the detials can be found in haru's papers.
! the number of the spatial orbial in np0, np1,np2,different sub-space.
!---------------------------------------------------------------------------------------------------
   system%np2 = np2
   system%np1 = np1
   system%np0 = np0 
 
! io_sdtq: S/D/SD, which are used in the calculations,   io_fcore: if the core orbital are used in the calculations
! if there are some core orbital are used, (np0)the number of the spatial orbital must not equal to zero.   

   system%io_sdtq = io_sdtq
   system%io_fcore = io_fcore
   write(*,*)
   write(*,*) '================================excitation and core =================================='
   write(*,*) 'the io_sdtq is equal to ',system%io_sdtq
   write(*,*)
   if(system%io_sdtq==0) then
        write(*,*) ' Only double excitation are considered in the calculations'
   elseif(system%io_sdtq==1) then
        write(*,*) 'Only the single excitation are considered in the calculations'
   elseif(system%io_sdtq==2) then
        write(*,*) 'single and double excitation are considered in the calculations'
   endif
   write(*,*)
   write(*,*) 'the io_fcore is equal to ',system%io_fcore

   if(system%io_fcore==0) then
    write(*,*)
    write(*,*) ' the calculation of tdrasscf, all the orbital are active, no core ' 
    write(*,*)
   elseif(system%io_fcore==1) then
    write(*,*)
    write(*,*) ' the calculation of tdrasscf, with hf-frozen-core app. '
    write(*,*)
   elseif(system%io_fcore==2) then 
    write(*,*)
    write(*,*) ' the calculation of tdrasscf, with casscf-frozen-core app. '
    write(*,*)
   endif
!
! total spatial orbital used in the calculations
!
   system%nptot = system%np2 + system%np1 + system%np0    ! the total number of spatial orbital used in the calculations



   system%nptot_fc = system%nptot
   system%n_offset = 0
   
   if(system%io_fcore/=0) then
     system%nptot_fc = system%nptot - system%np0
     system%n_offset = system%np0
   endif 

!
! if the number of spatial orbital are not enought to accomodate all the elelctrons
!
   if(2*system%nptot.lt.system%numelectron) then 
   
     write(*,*) 'the spatials cannot accommodate all the electrons'
     stop 
   endif     
!
! all the electron are put in the np0 space
!   
   if(2*system%np0.gt.system%numelectron) then
   
     write(*,*) 'all electrons cannot be put in the np0 space'
     stop 
   endif   
!
! the np1 space are empty, are not allowed
!
   if(system%np1.eq.0) then
   
     write(*,*) 'the np1 space do not allow to be empty'
     stop 
   endif   
!
! 
!
   if((system%np2.eq.0).and.(system%io_sdtq.ne.0)) then
   
     write(*,*) 'if io_sdtq/=0, then there are some spatial orbital in &
	 np2 space avaiable '

     stop 
   endif 


   if((2*system%np0.eq.system%numelectron) .and. (system%np1.ne.0)) then
   
     write(*,*) 'all electron are fill in np0 space'
     stop 
   endif 

   if((system%np0.eq.0) .and. (system%io_fcore.ne.0)) then
   
     write(*,*) 'all electron are fill in np0 space'
     stop 
   endif
   write(*,*) '================================spaces================================================'
   write(*,*)
   write(*,*) 'there are total spatial orbitals ======= ', system%nptot
   write(*,*)
   write(*,*) 'there are =====', system%np2, '=====spatial orbital in np2 space'
   write(*,*)
   write(*,*) 'there are =====', system%np1, '=====spatial orbital in np1 space'
   write(*,*)
   write(*,*) 'there are =====', system%np0, '=====spatial orbital in np0 space'
   write(*,*)
!------------------------------------------------------------------------------------------------------
! read the information of the basis set
!------------------------------------------------------------------------------------------------------  
!  two type, Gauss-lobatto FEDVR, Gauss-radau + Gauss+ lobatto FEDVR
   write(*,*) '================================fedvr information======================================' 

     fedvr3d%r_begin = r_begin
     fedvr3d%r_end = r_end     
     fedvr3d%number_of_element = number_of_element 
     fedvr3d%number_of_element_inner = number_of_element_inner !! elements in the inner region      

     If (number_of_element.lt.number_of_element_inner) then
        write(*,*) ''
        write(*,*) 'ERROR in readin.f90'
        write(*,*) 'The number of elements in the inner region MUST be smaller or equal'
        write(*,*) 'than the total number of elements'
        write(*,*) ''
        stop
     End If
     
!!$     If (abs(r_inner-r_out).gt.1d-10) then !! They should be the same
!!$        write(*,*) ''
!!$        write(*,*) 'ERROR in reading.f90'            
!!$        write(*,*) 'r_inner should be the same than r_out'
!!$        write(*,*) 'This should be include in the input file'
!!$        write(*,*) ''
!!$        stop
!!$     End If

     fedvr3d%r_inner=r_inner !! radius of the inner region. This inner region is used to define the FE-DVR.

     If (r_end.lt.r_inner) then
        write(*,*) ''
        write(*,*) 'ERROR in readin.f90'
        write(*,*) 'The radius of the inner region MUST be smaller or equal'
        write(*,*) 'than the maximum value of the radius'
        write(*,*) ''
        stop
     End If

     fedvr3d%store=store
     fedvr3d%r_out=r_out !! r_out defines the outer region. Used to perform the Fourier transform     
     fedvr3d%delta=delta !! r_out defines the outer region. Used to perform the Fourier transform     

    
     write(*,*) 
     write(*,'(a,1f12.3,a,1f12.3)')  ' from ', fedvr3d%r_begin , ' to ',fedvr3d%r_end 
     write(*,*)
     write(*,*)  'number of element is   ', fedvr3d%number_of_element
     write(*,*)
     write(*,*)  'number of sub in every element ', nb_sub

    if(fedvr3d%number_of_element < 1 ) then
      write(*,*) 'the number of element is not correct !!!'
      stop '--- error happen in initial_fedvr3d_radial' 
    endif

    if(fedvr3d%number_of_element > max_element) then
      write(*,*) '--- you must be kidding, the number of element is larger'
      write(*,*) '--- please increase the max_element'
      stop '--- error happen in initial_fedvr3d_radial'
    endif 

 ! spherical coordinate    
    if(fedvr3d%r_begin<=-0.00000001d0 .or. fedvr3d%r_begin >= fedvr3d%r_end - 0.000001d0) then
       write(*,*) 'the range r_begin-----r_end is not correct!!!'
       stop '--- error happen in initial_fedvr3d_radial' 
    endif

    if(nb_sub<=0) then
      write(*,*) 'input file nb_sub is not correct !!!!'
      stop 
    endif

     fedvr3d%fedvr_nb(1:fedvr3d%number_of_element) = nb_sub 

!
! in spherical coordinate, for angular part, l and m 
!
   write(*,*) '================================angular information==================================' 
    fedvr3d%l_max = l_max
    fedvr3d%m_max = m_max
    fedvr3d%l_max_core = l_max_core
    fedvr3d%l_max_coupled = l_max_coupled
    fedvr3d%m_max_core = min(fedvr3d%l_max_core,fedvr3d%m_max)

    If (l_max_coupled.lt.0) then
       write(*,*) 'ERROR IN readin.f90'
       write(*,*)
       write(*,*) 'l_max_coupled is smaller than 0'
       stop
    End If

    If (l_max_coupled.gt.2*l_max) then
       write(*,*) 'ERROR IN readin.f90'
       write(*,*)
       write(*,*) 'l_max_coupled is larger than 2*l_max'
       stop
    End If   

     if(fedvr3d%l_max<0 .or. fedvr3d%m_max>fedvr3d%l_max) then  
         write(*,*) 'you must be kidding, the input l_max,m_max is not correct !!!'
         stop 
     endif 

     if(fedvr3d%l_max_core.gt.fedvr3d%l_max) then  
         write(*,*) 'l_max is larger than the orbital quantum number of the core l_max_core'
         stop 
     endif 

     write(*,*) 'in the calculation the maximum l is :', l_max
     write(*,*)
     write(*,*) 'in the calculation the maximum m is :', m_max
     write(*,*)
     write(*,*) 'in the calculation the maximum l for the core is :', l_max_core

   write(*,*) '===============================momentum space information==============================='

   fedvr3d%k_max=k_max
   fedvr3d%k_min=k_min

   write(*,*) 'Maximum linear momentum, k_max= ', fedvr3d%k_max
   write(*,*) 'Minimum linear momentum, k_min= ', fedvr3d%k_min
   
   write(*,*) '===============================laser field information==============================='
!
! input laser parameter if time-dependent in real propagation
! 
    laser%tdornot = tdornot

    if(laser%tdornot) then
      write(*,*) 'the strong laser field is present here '

     laser%cep = cep
     laser%gauge = gauge
     laser%ellipticity = ellipticity
     laser%ex_max = ex_max
     laser%ey_max = ey_max
     laser%ez_max = ez_max
     laser%e_max = e_max
!
! laser1%ex_max,ey_max,ez_max : input, unit : w/cm2
!
    !! transform from w/cm2 to atomic unit                    
    laser%e_max = sqrt(laser%e_max)*5.338d-9 
    laser%ez_max = laser%e_max*sqrt(zone-ellipticity**2.0d0)
    laser%ex_max = laser%e_max*(ellipticity)
    laser%ey_max = sqrt(laser%ey_max)*5.338d-9

     

!  laser1%omega: unit, eV, 
!  laser1%lamda_nm, unit, nm
!  
    laser%omega = omega
    laser%lamda_nm = lamda_nm
    laser%ncycle = ncycle
    prop%ncycle_extra = ncycle_extra

    laser%omega = laser%omega/27.211396d0  ! transfer from eV to a.u.
 

    if(laser%omega<=0.0d0) then
      tcycle_fs = laser%lamda_nm/speedc*1.0d-2 
      tcycle_au = tcycle_fs/(2.4188d-2) 
      laser%omega = 2.0d0*pi/(tcycle_au)   
    endif

    if(laser%lamda_nm<=0.0d0) then
      tcycle_au = 2.0d0*pi/laser%omega 
      tcycle_fs = tcycle_au*2.4188d-2
      laser%lamda_nm = tcycle_fs*speedc*1.0d2
    endif
 
    energy_ev = 1.24d0/(laser%lamda_nm*1.0d-3)

     write(*,*) 'the central wavelength of laser used in the calcuations is, unit(nm) : ',&
        laser%lamda_nm
     write(*,*) 'the freqency is unit(a.u.) : ',laser%omega
     write(*,*) 'the energy of the photon,unit(eV): ', energy_ev 
     write(*,*) 'there are ', laser%ncycle, ' optical cycle in the calculations'
     write(*,*) 'Pulse duration (a.u.) ', dble(laser%ncycle)*2.0d0*pi/laser%omega
     write(*,*) 'Carrier Envelope Phase (radians): ', laser%cep

     write(*,*) '---------laser field parameter'
     write(*,*) 'laser%ex_max,  laser%ey_max,  laser%ez_max (atomic units)'
     write(*,*)  laser%ex_max, laser%ey_max, laser%ez_max
     write(*,*)
     write(*,*) 'Ellipticity: ',ellipticity
     write(*,*)
     write(*,*) 'Duration of the propagation (a.u.) ', (dble(laser%ncycle)+dble(prop%ncycle_extra))*2.0d0*pi/laser%omega
     write(*,*)

   else

     write(*,*) 'the strong laser field is absent, the imag. time propagation '

   endif

  write(*,*) '================================propagation==================================' 

!=======================================================================================
! prop t0, t_end, hstep
!=======================================================================================
   prop%hstep = hstep 
   prop%ncycle_extra = ncycle_extra
   prop%solver_which  = solver_which 
   prop%icount_max = icount_max
   prop%icount_store = icount_store
   prop%icount_ionization = icount_ionization
   prop%stiffness = stiffness
   
   If (icount_max.lt.icount_store) then
      print*, 
      print*, 'ERROR IN reading.f90'
      print*, 'icount_max must be larger than icount_store'
      print*, 
      stop
   End If

   If (icount_max.lt.icount_ionization) then
      print*, 
      print*, 'ERROR IN reading.f90'
      print*, 'icount_max must be larger than icount_ionization'
      print*, 
      stop
   End If

   If (icount_ionization.gt.9999.or.icount_ionization.lt.1) then
      print*, 
      print*, 'ERROR IN reading.f90'
      print*, 'icount_ionization must be in the interval [1:9999]'
      print*, 
      stop
   End If

   If (abs(dble(icount_max/icount_store)-dble(icount_max)/dble(icount_store)).gt.1d-10) then
      print*, 
      print*, 'ERROR IN reading.f90'
      print*, 'icount_max must be a multiple of icount_store'
      print*, 
      stop
   End If

         write(*,*) 
   select case (fedvr3d%store)
      case(0)
         write(*,*) 'The Mean field operator is not stored'
      case(1)
         write(*,*) 'The Mean field operator is stored'
      case(2)
         write(*,*) 'The Mean field operator is calculated using the coupled representation'
      case default
         write(*,*) 'ERROR in the value of store in the input file'
         stop
      end select
         write(*,*) 

   if(prop%stiffness==1.and.prop%stiffness==2) then
    write(*,*)
    write(*,*) 'the stiffness effect are considered in the propagation '
    write(*,*)
  else
   write(*,*) 
   write(*,*) '----------Be careful about the step ----------------------'
   write(*,*) 'the stiffness effect are not considered in the propagation'
   write(*,*) 
endif

!! Files to store the dipole and the dipole acceleration

             Open(101,file='dipole')
write(101,*) '# Time                        <x> (a. u.)               <y> (a. u.)              <z> (a. u.)'
             write(101,*) ''

             Open(108,file='dipole_acceleration')
             write(108,*) '# Time                       <\partial_tt x> (a. u.)    <\partial_tt y> (a. u.)   <\partial_tt z> (a. u.)'
             write(108,*) ''


 
   if(laser%tdornot) then
      prop%t0 =  -laser%ncycle*2.0d0*pi/laser%omega*.5d0!! Initial time: -Duration of the pulse/2
!!$     prop%t_end = t_end !! Final time taken from the input
      t_end_hhg = laser%ncycle*2.0d0*pi/laser%omega*.5d0 !! Final time: Duration of the pulse/2
      prop%t_end = (dble(laser%ncycle)/2.0d0+dble(prop%ncycle_extra))*2.0d0*pi/laser%omega !! Final time taken as the time of ncycle+1 cycles.      
      t_end= prop%t_end 

!! To perform the HHG we change the number of points of icount_store, to have a power of two

          if(i_hhg==1) then
             exponent_2=int(log(dble(prop%icount_store))/log(dble(2.0d0)))+1
             prop%icount_store=2**exponent_2

             !! Initialize the files to store the dipole and the dipole acceleration
          End if


!! TIME STEPS IN THE PULSE TO CALCULATE THE HHG

!! Number of time steps associated to the time step described
          prop%icount_max = int(( t_end_hhg - prop%t0 )/ prop%hstep) 

!! Change to a number of time steps multiple of the stored points
          prop%icount_max = prop%icount_max-mod(prop%icount_max,prop%icount_store)+prop%icount_store

!! New time step (associated to the number of maximum points during the pulse)
         
          prop%hstep=(t_end_hhg - prop%t0)/dble(prop%icount_max)
          
          prop%icount_max_pulse=prop%icount_max !! Number of loops in the pulse

!! TIME STEPS DURING THE PROPAGATION TO STORE THE IONIZATION
         
          prop%icount_max = int(dble(prop%icount_max_pulse)/dble(laser%ncycle)*dble(laser%ncycle+prop%ncycle_extra))
          !! prop%icount_max is proportional to i_count_max_pulse.

          prop%icount_max = prop%icount_max-mod(prop%icount_max,prop%icount_ionization)
          !! prop%icount_max has to be scaled to be a multiple of icount_ionization

print*, ''
print*, '====================================================================='
print*, 'Parameters in the time propagation'
print*, '----------------------------------'
print*, 
print*, 'Time step (prop%hstep):', prop%hstep
print*, 'Number of time steps in the propagation (icount_max):',prop%icount_max
print*, 'Number of time steps during the pulse (icount_max_pulse):',prop%icount_max_pulse
print*, '====================================================================='
print*, ''

     write(*,'(a,1f12.3,1f12.3)') 'the begin and end time is (laser duration time + field free region)', prop%t0,prop%t_end !! note that we add the time of #ncycle_extra cycle

  endif !! End if it is time propagation
   
   write(*,'(a,1f12.6)') 'the initial integral step is', prop%hstep
   write(*,*) '==========================absorb potential=================================='
!======================================================================================
! absorb part
!======================================================================================
absorb%r_right = r_right
absorb%r_left = r_left
absorb%ramp = ramp
absorb%strength = strength
 write(*,'(a,f12.3,a,f12.3)') 'absorb right position ', absorb%r_right, ' absorb left position ',absorb%r_left
 write(*,*)
 write(*,'(a,f12.3,a,f12.3)') 'ramp parameter ', ramp, ' absorb strength ',absorb%strength
 write(*,*) '================================analysis part=================================='
!======================================================================================
! analysis part
!====================================================================================== 
   analysis%i_energy = i_energy 
   if(i_energy==1) then
     file1 = 71
!!$     open(file1,file='energy.txt',status='unknown')

     write(energy_file,'(a12,es8.2,a4,i3.3,a4,i2.2,a6,i2.2,a6,i2.2,a3,i2.2,a4,i2.2,a5,i2.2,a5,i2.2,a5,i2.2,a6,i2.2)') 'energy_rmax_',fedvr3d%r_end,'_el_',fedvr3d%number_of_element,'_nb_',fedvr3d%fedvr_nb(1),'_lmax_',fedvr3d%l_max,'_mmax_',fedvr3d%m_max,'_Z_',int(system%zz),'_ne_',int(system%numelectron),'_np2_',system%np2,'_np1_',system%np1,'_np0_',system%np0,'_sdtq_',system%io_sdtq

     open(file1,file=energy_file,status='unknown')

     write(file1,*) '# Time (atomic units)       Energy (atomic units)'
     write(file1,*) ''

     write(*,*)
     write(*,*) ' the energy of system will be calculated in file energy.txt '
     write(*,*)
   endif
   
   write(*,*)
   analysis%i_auto = i_auto

!! Calculation of the autoprojection

   if(i_auto==1) then
     file2 = 72
     open(file2,file='auto.txt',status='unknown')
     write(*,*) ' the auto-function will be calculated in file auto.txt '
   endif
   write(*,*)


   analysis%i_ionization = i_ionization
   if(i_ionization==1.and.laser%tdornot) then
      file3 = 73
      Open(file3,file='ionization.txt')
     write(*,*) ' the total ionization probabitity will be calculated in file ionization.txt '
     write(73,*) '# 1st: Time (atomic units)'
     write(73,*) '# 2nd: Norm (normalize to the number of electrons)'
     write(73,*) '# 3rd: number of electrons - electrons from 0 to r_out=', fedvr3d%r_out,'to infinity'
     write(73,*) '# 4th: Norm of the Two-Body density (normalize to n_e(n_e-1)/2)'
     write(73,*) '# 5th: n_e(n_e-1)/2 - yield in the inner region for the Two-Body'
     write(73,*) '# '
     write(73,*) ''

   endif
   write(*,*)
   analysis%i_hhg = i_hhg
   if(i_hhg==1) then
     file4 = 74
     open(file4,file='hhg.txt',status='unknown')
     write(*,*) ' the hhg will be calculated in file hhg.txt '
   endif
   write(*,*)
   analysis%i_twoe_density = i_twoe_density
   if(i_twoe_density==1) then
     file5 = 75
     open(file5,file='i_twoe_density.txt',status='unknown')
     write(*,*) ' the i_twoe_density of system will be calculated in file i_twoe_density.txt '
   endif
   write(*,*)
   analysis%i_momentum  =i_momentum
   if(i_momentum==1) then
     file6 = 76
     open(file6,file='i_momentum.txt',status='unknown')
     write(*,*) ' the momentum spectrum of system will be calculated in file i_momentum.txt '
   endif   
   write(*,*)
   analysis%i_atas = i_atas
   if(i_atas==1) then
     file7 = 77
     open(file7,file='atas.txt',status='unknown')
     write(*,*) ' the atas will be calculated in file atas.txt '
   endif
  return
end subroutine readin

