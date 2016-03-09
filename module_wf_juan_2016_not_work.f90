
module wfunction
  use global
  use auxiliary
  use operator_3d

  implicit none
  complex(kind=k2),allocatable,save :: acoff(:),phi(:,:)
!check to debug CMFA
  complex(kind=k2),allocatable,save :: phi_cmfa(:,:) !! To check Constant Mean Field approximation
!check to debug
  complex(kind=k2),allocatable,save :: acoff0(:),phi0(:,:)
  complex(kind=k2),allocatable,save :: acoff_old(:),phi_old(:,:)
  integer, allocatable,save :: phi_rang(:)!! This array gives the basis set functions which contribute for any of the orbitals
!#1 is the basis set wave function included.
  integer :: phi_rang_total !! This integer gives the maximun number of basis set functions which have a contribution.
  integer, allocatable, save :: phi_rangang_1(:),phi_rangang_2(:)
  integer,  save :: phi_rangang_total
!! phi_rangang1(#1) and phi_rangang_2(#1) gives the global function basis set of the fedvr+angle+angle contribution to the two body integrals.
  character*1 :: auxiliar_char !! to read the groundstate wf
  integer :: auxiliar_np2, auxiliar_np1, auxiliar_np0, auxiliar_io_sdtq, auxiliar_ne, auxiliar_l_max, auxiliar_m_max
  real (kind=k1) :: auxiliar_zz
  real (kind=k1) :: auxiliar_r_begin, auxiliar_r_end, auxiliar_inner_elements, auxiliar_nb_sub
  
  contains
!
! initial the wf:  acoff and phi
!
     subroutine wf_beginning
     implicit none
     integer :: idicp,jdicp,kdicp,ldicp,ierr,ix,iy,ip,info
     real(kind=k1) :: mnorm 
     logical :: lexist
     real(kind=k1),allocatable,save :: w_temp(:),fv1(:),fv2(:)
     real(kind=k1),allocatable,save :: z_temp(:,:),h_temp(:,:)
     real(kind=k1), allocatable :: h_band(:,:)
     complex(kind=k2), allocatable :: phi_state(:)
     real*8 :: z
     integer :: n, lh, mh
     integer :: i_angle,i_r,i_row_r,i_column_r,i_row,i_column
     integer :: counter, ii,kd
     real(kind=k1) :: E_biggest,step_biggest,sumtemp
     real(kind=k1) :: e_min
     real(kind=k1) ::  rstep_smallest,  rstep_aux
     !! DEFINITIONS FOR THE BANDED STORAGE
     integer :: m_band!! number of eigenvalues found
     
!!$     real(kind=k1),allocatable :: w_band(:) !! eigenvalues found
!!$     real(kind=k1), allocatable:: z_band(:,:) !! eigenfunctions found
     !! the first argument corresponds to the component of the eigenvector
     !! and the second one is the order of the eigenvector.
     real(kind=k1), allocatable:: work_band(:) !!Work arrays
     integer, allocatable:: iwork_band(:) !!Work arrays
     integer, allocatable:: ifail_band(:)
     integer :: ldq
     real(kind=k1), allocatable:: q_band(:,:) !!matrix used in the reduction
     real :: seed
     
     allocate(acoff(n_string_ab), phi(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) )

     allocate(acoff0(n_string_ab), phi0(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) )
     
     if(analysis%i_auto==1) then
       allocate(acoff_old(n_string_ab), phi_old(fedvr3d%nb_r*fedvr3d%nb_angle,&
       system%nptot) )
     endif

 if(prop%stiffness==1) then
!!$  allocate(w_temp(fedvr3d%nb_r*fedvr3d%nb_angle),fv1(fedvr3d%nb_r*fedvr3d%nb_angle),&
!!$  fv2(fedvr3d%nb_r*fedvr3d%nb_angle))

!!$  allocate(z_temp(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle))
  allocate(h_temp(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle))

  allocate(h_stiffness(fedvr3d%nb_r*(fedvr3d%nb_r+1)/2,fedvr3d%nb_angle))
endif




rstep_smallest=fedvrx_global(size(fedvrx_global))

Do ii=1,size(fedvrx_global)-1
    rstep_aux=fedvrx_global(ii+1)-fedvrx_global(ii) 
    !! set the step in the radial part for the node in ii and ii-1
    If (rstep_aux.lt.rstep_smallest) rstep_smallest=rstep_aux
End Do

!! For the stiffness

step_biggest=(rstep_smallest/3.0d0)**2.d0 !! Biggest step

select case (prop%stiffness)

   case(0)
    
!!     step_biggest =  1.0d0/w_temp(fedvr3d%nb_r*fedvr3d%nb_angle)*3.1415926d0*0.5d0


!!$step_biggest=(.25d0+dble(fedvr3d%l_max*(fedvr3d%l_max+1)))/(rstep_smallest)**2.0d0
!!$
!!$step_biggest=1.0d0/step_biggest !!Juan: I think this is a bad choise

       write(*,*) 
       write(*,*)  'the advice step is : ', step_biggest 
       write(*,*)  'Current step is: ', prop%hstep
     if(prop%hstep > step_biggest) then
       write(*,*)  'Warning: The step should be decrease '  
       write(*,*)
!       stop
     endif
     write(*,*)

case(1) !!stiffness

   write(*,*) 
   write(*,*)  'the advice step is : ', step_biggest 
   write(*,*)  'Current step is: ', prop%hstep
   write(*,*) 

   write(*,*)
   write(*,*) 'Stiffness is taken into account'
   write(*,*)

    h_temp = zzero

!!=====================================================================
!! Building the matrix to diagonalize.
!! Wenliang: Build a full matrix to diagonalize

    do i_angle =1, fedvr3d%nb_angle        
       do i_r = 1,n_total_kinetic
          
          i_row_r = index_kinetic_basis(i_r,1)
          i_column_r = index_kinetic_basis(i_r,2)
          i_row = (i_angle-1)*fedvr3d%nb_r + i_row_r
          i_column = (i_angle-1)*fedvr3d%nb_r + i_column_r
          
          if(i_row == i_column) then 
             h_temp(i_row,i_column)=(tmat_3d(i_row_r,i_angle) + vmat_radial(i_row_r) + tmat_radial(i_r))       
          endif
          
          if(i_row /= i_column) then
             h_temp(i_row,i_column) =  tmat_radial(i_r)
             h_temp(i_column,i_row) =  tmat_radial(i_r) 
          endif
          
       enddo
    enddo

!! Juan: Storing the matrix as banded.
!! First, we use h_temp as an entry to rebuild it in another storage.
!! *** This must be changed at some point.

!! The way of storing the matrix may be check in Lapack help for subroutine dsbevx

!! The number of superdiagonals above the diagonal is, as maximum, the number of nodes in each element, fedvr3d%fedvr_nb.

    kd=fedvr3d%fedvr_nb(1)

    allocate(h_band(1:kd+1,1:fedvr3d%nb_r*fedvr3d%nb_angle))

    h_band=zero

    Do i_column=1,fedvr3d%nb_r*fedvr3d%nb_angle !! Run for the columns
       Do i_row=max(1,i_column-kd),i_column !! Run for the arrows
          h_band(kd+1+i_row-i_column,i_column)=h_temp(i_row,i_column)
       End Do
    End Do

!! DIAGONALIZATION OF THE BAND MATRIX
    
    allocate(q_band(1:fedvr3d%nb_r*fedvr3d%nb_angle,1:fedvr3d%nb_r*fedvr3d%nb_angle))
    allocate(w_band(1:fedvr3d%nb_r*fedvr3d%nb_angle))
    allocate(z_band_temp(1:fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle))
    allocate(work_band(1:fedvr3d%nb_r*fedvr3d%nb_angle*7))
    allocate(iwork_band(1:fedvr3d%nb_r*fedvr3d%nb_angle*5))
    allocate(ifail_band(1:fedvr3d%nb_r*fedvr3d%nb_angle))

    write(*,*) 'Lowest bound of the energy',-((system%zz+.1d0)**2.0d0)/2.0d0!Juan
    E_biggest = 1.0d0/prop%hstep! Juan

    write(*,*) 'Largest energy allowed',E_biggest 
  
    call DSBEVX( 'V', 'V', 'U', fedvr3d%nb_r*fedvr3d%nb_angle, kd, h_band,&
         kd+1, q_band, fedvr3d%nb_r*fedvr3d%nb_angle, -((system%zz+.1d0)**2.0d0)/2.0d0, E_biggest, 1,&
         2, 1d-20, m_band, w_band, z_band_temp, fedvr3d%nb_r*fedvr3d%nb_angle,&
         work_band, iwork_band, ifail_band, info )
      if(info/=0) then
        write(*,*) 'error in subroutine DSBEVX'
        write(*,*) 'info=', info
        stop 
      endif

    deallocate(q_band,iwork_band,ifail_band)


!!=======================================================



!! END of Building the matrix to diagonalize
!!=====================================================================

!!=====================================================================
!! diag. the full matrix matrix. Juan: I comment this

!!$      call rs(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle,h_temp,& 
!!$      w_temp,1,z_temp,fv1,fv2,ierr)
!!$
!!$      print*, w_temp(1)
!!$
!!$
!!$      if(ierr/=0) then
!!$        write(*,*) 'error in subroutine initial_wp'
!!$        stop 
!!$      endif
!!$
!!$!
!!$! stiffness one body operatror
!!$!     
!!$!     E_biggest = 1.0d0/prop%hstep  ! Juan
!!$     E_biggest = pi*0.5d0/prop%hstep! Juan
!!$
!!$     do idicp =1,fedvr3d%nb_r*fedvr3d%nb_angle   
!!$       If( w_temp(idicp).gt.E_biggest ) then
!!$          w_temp(idicp) = zero
!!$       End If
!!$     enddo 
!!$
!!$     do idicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! number of the state in the row
!!$      do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! number of the state in the column
!!$       sumtemp = zero
!!$        do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! Number of eigenvalues
!!$          sumtemp = sumtemp + z_temp(idicp,kdicp)*w_temp(kdicp)*z_temp(jdicp,kdicp)
!!$        enddo
!!$       h_temp(idicp,jdicp) = sumtemp
!!$      enddo
!!$     enddo
!! End of diag. the full matrix matrix. Juan: I comment this
!!=====================================================================
!!$

!! Build the matrix after the stiffness

    h_temp=zzero

     do idicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! number of the state in the row
      do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! number of the state in the column
       sumtemp = zero
        do kdicp =1,m_band !! Number of eigenvalues
          sumtemp = sumtemp + z_band_temp(idicp,kdicp)*w_band(kdicp)*z_band_temp(jdicp,kdicp)
        enddo
       h_temp(idicp,jdicp) = sumtemp
      enddo
     enddo

!
! store the stiffness matrix in the h_stiffness matrix
!
 
   do idicp =1,fedvr3d%nb_angle !! Running in the angular part
      ldicp =0 
      do jdicp=(idicp-1)*fedvr3d%nb_r+1,(idicp-1)*fedvr3d%nb_r + fedvr3d%nb_r  !! Running in the radial part corresponding to the angular part idicp.
       do kdicp =(idicp-1)*fedvr3d%nb_r+1,jdicp !! 
       ldicp = ldicp +1
        h_stiffness(ldicp,idicp) = h_temp(jdicp,kdicp)
      enddo
   enddo
 enddo

!! ! diag. the matrix. Juan: Wenliang
!!$     call rs(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle,h_temp,&
!!$      w_temp,1,z_temp,fv1,fv2,ierr)
!!$
!!$     print*, w_temp(6)
!! ! diag. the matrix. Juan: Wenliang

 write(*,*) 'Number of eigenfunctions found', m_band

 If (m_band.lt.fedvr3d%nb_r*fedvr3d%nb_angle/2) then
    write(*,*)
    write(*,*) '=================================================='
    write(*,*) 'Decrease the time step'
    write(*,*) 'The number of the eigenvalues found is very small'
    write(*,*) '=================================================='
    stop
 End If

 allocate(z_band(1:fedvr3d%nb_r*fedvr3d%nb_angle,1:m_band))

 z_band=zzero

 Do idicp=1, size(z_band(:,1))
    Do jdicp=1, size(z_band(1,:))
       z_band(idicp,jdicp)=z_band_temp(idicp,jdicp)
    End Do
 End Do


deallocate(z_band_temp)

!!$     print*, 'module_wf'
!!$     stop

!!$
!
! diag. hcore get spatial orbital
!     

     deallocate(h_temp)

case(2) !!stiffness with an optimized procedure

   write(*,*) 
   write(*,*)  'the advice step is : ', step_biggest 
   write(*,*)  'Current step is: ', prop%hstep
   write(*,*) 
  
   E_biggest=1.0d0/(1.0*prop%hstep) !! Larger energy allowed   

!! In the case of Imaginary time propagation, we turn off the stiffness, because all the states are included in the random wave function
   if (.not.laser%tdornot) then

      write(*,*) 'Imaginary time propagation: STIFFNESS IS OFF'
      E_biggest=1e15
   End if

   write(*,*)
   write(*,*) 'Largest energy allowed',E_biggest 
   write(*,*)
   write(*,*)
   write(*,*) 'Stiffness is taken into account!!'

   e_min=-((system%zz+.1d0)**2.0d0)/2.0d0 !! Lower bound of the one body
   !! hamiltonian eigenenergies.

!! now we obtain the stiffness Hamiltonian and the initial states.

   call hl_stif(e_biggest,e_min, system%zz,fedvrx_global,index_kinetic_basis,tmat_radial,hl_stiffness,phi) 

End select

!
! for real time propagation, we read from wf0.txt
!

if(laser%tdornot) then

    prop%chstep = dcmplx(prop%hstep,0.0d0)

    inquire(file='wf0.txt',exist=lexist)
    if (.not.lexist) then
       write(*,*) '---The ground state (wf0.txt) does not exist---'
       stop 
    endif

    !!Read the parameters fo the ground state to check that they are compatible

    open(ifile,file='wf0.txt')
    
    Do idicp=1,5  !! read the heading
       read(ifile,*) auxiliar_char
    End Do
    
    read(ifile,*) auxiliar_zz, auxiliar_ne !! Nuclear charge and number of electrons

    If (abs(auxiliar_zz-system%zz).gt.1d-4) then
       write(*,*) 'THE NUCLEAR CHARGE IS NOT THE SAME FOR THE SYSTEMS'
       stop
    End If
    If (auxiliar_ne.ne.system%numelectron) then
       write(*,*) 'THE NUMBER OF ELECTRONS IS NOT THE SAME FOR THE SYSTEMS'
       stop
    End If

    Do idicp=1,3  !! read the empty spaces and titles
       read(ifile,*) auxiliar_char
    End Do

    read(ifile,*) auxiliar_np2, auxiliar_np1, auxiliar_np0, auxiliar_io_sdtq

    If (auxiliar_np2.ne.system%np2) then
       write(*,*) 'RAS SCHEME DIFFERS ON NP2'
       stop
    End If

    If (auxiliar_np1.ne.system%np1) then
       write(*,*) 'RAS SCHEME DIFFERS ON NP1'
       stop
    End If

    If (auxiliar_np0.ne.system%np0) then
       write(*,*) 'RAS SCHEME DIFFERS ON NP0'
       stop
    End If

    If (auxiliar_io_sdtq.ne.system%io_sdtq) then
       write(*,*) 'RAS SCHEME DIFFERS ON THE EXCITATION SCHEME'
       stop
    End If

    Do idicp=1,3  !! read the empty spaces and titles
       read(ifile,*) auxiliar_char
    End Do

    read(ifile,*) auxiliar_r_begin, auxiliar_r_end, auxiliar_inner_elements, auxiliar_nb_sub

    if (abs(auxiliar_r_begin-fedvr3d%r_begin).gt.1.0d-3) then
       write(*,*) 'NOT THE SAME r_begin'
    End if

    if (abs(auxiliar_r_end-fedvr3d%r_end).gt.1.0d-3) then
       write(*,*) 'NOT THE SAME r_end'
    End if

    if (auxiliar_inner_elements.ne.fedvr3d%number_of_element_inner) then
       write(*,*) 'FEDVR IS NOT THE SAME'
       stop
    End if

    if (auxiliar_nb_sub.ne.fedvr3d%fedvr_nb(1)) then
       write(*,*) 'FEDVR IS NOT THE SAME'
       stop
    End if

    Do idicp=1,3  !! read the empty spaces and titles
       read(ifile,*) auxiliar_char
    End Do

    read(217,*) auxiliar_l_max,auxiliar_m_max

    If (auxiliar_l_max.gt.fedvr3d%l_max) then
       write(*,*) 'l_max of the groundstate is larger than l_max for the propagation'
    End If

    If (auxiliar_m_max.gt.fedvr3d%m_max) then
       write(*,*) 'm_max of the groundstate is larger than m_max for the propagation'
    End If

    print*, 'module_wf_juan.f90'
    stop

!
! reading the acoff from the imag.txt file    
! 

        open(ifile,file='imag.txt') !!! Changed by Juan

      do idicp =1, n_string_ab
       read(ifile,*) acoff0(idicp)
      enddo
!
! readin the phi orbitals from the imag.txt file
! 
     do idicp =1, system%nptot
       do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
        read(ifile,*) phi0(jdicp,idicp)
      enddo 
    enddo

  acoff = acoff0
  phi = phi0

  if(analysis%i_auto==1) then
   acoff_old = acoff0
   phi_old = phi0
  endif

close(ifile)

!===============================================================================
! imaginary time propgation 
!===============================================================================

else

  prop%chstep = dcmplx(0.0d0, -1.0_k1*prop%hstep) !! set the imaginary step

  open(ifile,file='imag.txt',status='unknown')

!
! imaginary propagation
! prepare the initial A cofficient
!

!! Initialize the seedn to choose the random initial guess function
         call cpu_time(seed)
         call srand(int(10000.0d0*seed))
!! End of initialize the seed


   do idicp =1, n_string_ab
!     acoff(1) = dcmplx(1.0d0)
   enddo

   Select case(prop%stiffness)
      case(20)
         continue !! calculated in hl_stif
         !! The orbitals are the eigenvectors of the one body Hamiltonian
      case default

!
! initial the spatial orbitals using radial hydrogen functions
!

!!$   allocate(phi_state(1:fedvr3d%nb_r*fedvr3d%nb_angle))
!!$     counter=0

!!$Do n=1,system%nptot !! run in principal quantum number
!!$   Do lh=0,min(n-1,fedvr3d%l_max),1 !! run in the total angular momentum
!!$      Do mh=-min(fedvr3d%m_max,lh),min(fedvr3d%m_max,lh),1 !! run in the magnetic quantum number
!!$         counter=counter+1
!!$         call radial_hydrogen_dvr(fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m,fedvr_w,n,lh,mh,system%zz,&
!!$phi_state)
!!$         phi(:,counter)=phi_state         
!!$
!!$         If (counter.eq.system%nptot) exit
!!$      End Do
!!$      If (counter.eq.system%nptot) exit
!!$   End Do
!!$   If (counter.eq.system%nptot) exit
!!$End Do
!!$
!!$deallocate(phi_state)

         Do counter=1,size(phi(1,:))               
            Do jdicp= 1,size(phi(:,1)) !! Random choice
                  phi(jdicp,counter)=(rand(0)-0.50d0)*2.0d0
            End Do
         End Do
!!$         
   End Select
!!$
!
! you can select any one to use
! We construct substract some noise of the ortogonality

     call schmidtortho2(phi,fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot,ierr,mnorm)

      if(ierr /= 0) then
        write(*,*) 'error in schmidtortho2'
        stop 
      endif

     phi0 = phi
     acoff0 = acoff

  endif
  return
 end subroutine wf_beginning

!!========================================================================
!! Functions

!! The function laguerrel calculates the value of the Associated Laguerre Polynomials of degree 'nu' and index 'alpha' at the real value 'x'.

Function laguerre(n,alpha, x)
implicit none
!!INPUT
integer, intent(in) :: n !! degree of the Associated Laguerre Polynomials.
integer, intent(in) :: alpha !! index of the Associated Laguerre Polynomials.
real*8, intent(in) :: x !! argument of the Associated Laguerre Polynomials.
!!OUTPUT
real*8 :: laguerre !! Associated Laguerre Polynomial.
!!Auxiliar variables
integer :: ii, jj, kk,nu
real*8 :: lnu_minus_1, lnu

lnu_minus_1=1.0d0
lnu=-x+dble(alpha+1)

select case (n)
case(0)
   laguerre=lnu_minus_1 !! Associated Laguerre Polynomials of degree 0
   
case(1)
   laguerre=lnu !! Associated Laguerre Polynomials of degree 1
   
case default

   Do nu=1,n-1 !! Loop for the recursive relation
      laguerre=0.0d0       !! initialize the value of the polynomial
      laguerre=(dble(2*nu+alpha+1)-x)*lnu-dble(nu+alpha)*lnu_minus_1
      laguerre=laguerre/(dble(nu+1))
      
      !! update the auxiliar terms
      
      lnu_minus_1=lnu
      lnu=laguerre
   End Do

   continue

end select

return
    End Function laguerre

!! The function radial_hydrogen gives the radial function for a hydrogen-like atom with a nuclear charge z WITHOUT NORMALIZED and in atomic units. 

Function radial_hydrogen(n,l,z,r)
implicit none
!!INPUT
integer, intent(in) :: n, l !! principal and orbital angular momentum quantum numbers
real*8, intent(in) :: z !! nuclear charge
real*8, intent(in) :: r !! radial coordinate
!!OUTPUT
real*8 :: radial_hydrogen

If (n-l-1.lt.0) then
   print*, 'ERROR IN radial_hydrogen'
   print*, 'n must be greater than l+1'
end If

radial_hydrogen=0.0d0 !!initialize the function

radial_hydrogen=laguerre(n-l-1,2*l+1,2.0d0*z*r/dble(n))*exp(-z*r/dble(n))
radial_hydrogen=(2.0d0*z*r/dble(n))**(dble(l+1))*radial_hydrogen

End Function radial_hydrogen

!! Subroutine to write the radial hydrogen-like function in the DVR basis together with the angular part.

 subroutine radial_hydrogen_dvr(xglobal,nb_angle,element,basis,l,m, weights,n,lh,mh,z,phiuni)
  implicit none

!!INPUT
!! phi in terms of the fedvr global basis.
 real*8, intent(in), allocatable :: xglobal(:)
!! xglobal is the radial basis fedvr
 real*8, intent(in), allocatable :: weights(:,:)
!! weights are the weights for the element in 1st argument and basis
!! in the element for the 2nd argument
 integer, intent(in) :: element(:), basis(:) 
!! element and basis in this element of the function chosen
 integer, intent(in) :: nb_angle
!! number of angular functions
 integer, intent(in), allocatable :: l(:),m(:)
!! L and M for the angular function in the argument.
integer, intent(in) :: n
!! n is the principal quantum number of the hydrogen-like function
integer, intent(in) :: lh,mh
!! lh and mh are the angular momentum and the magnetic quantum number of the hydrogen-like function.
real*8,intent(in) :: z
!! z is the nuclear charge
!!OUTPUT
 complex(kind=k2), intent(out) :: phiuni(:)
!! this is the radial wavefunction squared
!!Auxiliary aspects
integer :: i,j,k,ko, counter, jangle

phiuni=zzero

Do j=1,nb_angle !! to set the angular momentum basis
   If (l(j).eq.lh.and.m(j).eq.mh) then
      jangle= j 
      exit
   End If
End Do

Do i=1,size(xglobal) !! run for the radial coordinate
   counter=(jangle-1)*size(xglobal)+i
   phiuni(counter)=radial_hydrogen(n,lh,z,xglobal(i))
   If (i.lt.size(xglobal)) then
      If (basis(i).lt.basis(i+1)) then !! basis(i) is not the last
           !! node of the element.
         phiuni(counter)=phiuni(counter)*sqrt(weights(element(i),basis(i)))
      Else !! basis(i) is the last node of the element.
         phiuni(counter)=phiuni(counter)*sqrt(weights(element(i),basis(i))+weights(element(i+1),1))
      Endif
   End If
End Do

!! normalize the wavefunction

phiuni=phiuni/sqrt(dot_product(phiuni,phiuni))

return

End subroutine radial_hydrogen_dvr

!!====================================================================
!! Subroutines to prepare matrix to stiffness.
!!====================================================================

!!--------------------------------------------------------------------
!! hl_stif builds a matrix getting rid of the eigenvectors of the
!! one body hamiltonian with the higher eigenvalues.
!! hl_stif also stores the lowest states.
!!--------------------------------------------------------------------

subroutine hl_stif(energy_cutoff,energy_min, z,xglobal,index_t,t_r,hl,orb)
  implicit none
  
  !! INPUT
  
  real(kind=k1),intent(in) :: energy_cutoff, energy_min
  !! Energy cutoff used in the stiffness and a lower bound of the 
  !! energy.
  
  real*8, intent(in), allocatable :: xglobal(:)
  !! xglobal is the radial basis fedvr
  
  real(kind=k1), intent(in) :: z 
  !! nuclear charge
  
  integer, allocatable, intent(in) :: index_t(:,:) 
  !! index for the radial kinetic energy.
  !! Arguments
  !!
  !! #1: Number of the matrix elements of the radial kinetic 
  !!     energy
  !! #2: If #2=1-> bra function
  !!     If #2=2-> ket function
  
  real(kind=k1), allocatable, intent(in) :: t_r(:)
  !! matrix element of the
  
  !! OUTPUT
  
  real(kind=k1), allocatable, intent(out) :: hl(:,:,:)
  !! Hamiltonian after the stiffness procedure.
  !! The hamiltonian matrix is the same for each value of the 
  !! angular momentum l. 
  !! Arguments
  !!
  !! #1: FE-DVR function of the bra.
  !! #2: FE-DVR function of the ket.
  !! #3: Angular momentum of the Hamiltonian, l.
  !!
  
  complex(kind=k2), allocatable,intent(inout) :: orb(:,:)
  !! orb(,) stores the eigenvectors with the lowest energy.
  !! orb stores np0+np1+np2 eigenvectors.
  !! Arguments
  !!
  !! #1: FE-DVR+Angular momentum
  !! #2: Number of orbital
  
  !! AUXILIAR VARIABLES
  
  integer :: kd !! distance to the superdiagonal.
  integer :: il, it !! indexes for the loops
  integer :: i_row, i_column
  real(kind=k1),allocatable :: htemp(:,:)
  !! #1: Distance to the super diagonal+1
  !! #2: FE-DVR function of the ket.
  
  integer :: i_eigen!! run over states
  
  !! Auxiliar for the diagonalization

  integer :: m_band!! number of eigenvalues found    
  real(kind=k1), allocatable:: work_band(:) !!Work arrays
  integer, allocatable:: iwork_band(:) !!Work arrays
  integer, allocatable:: ifail_band(:)
  integer :: ldq, info
  real(kind=k1), allocatable:: q_band(:,:) !!matrix used in the reduction

  !! Auxiliar to calculate the initial orbitals
  complex(kind=k2), allocatable:: orb_aux(:,:) !! store as the orbitals.
  real(kind=k1), allocatable:: orb_ener(:) !! energies for each orbital.   
  integer :: norbital, nbasis
  integer :: ij, ik, jk, jk1
  integer :: multiplicity !! number of states with the same energy.
  integer :: im, ilprime !! loop m and auxiliar in ilprime
  integer :: orb_ener_order, counter
  

  kd=fedvr3d%fedvr_nb(1) !! the distance to the superdiagonal.
  
  allocate(hl(1:fedvr3d%nb_r,1:fedvr3d%nb_r,0:fedvr3d%l_max))
  
  hl=zero
  
  !! allocate the auxiliar matrix to diagonalize
  
  allocate(htemp(1:kd+1,1:fedvr3d%nb_r))
  
  htemp=zero

  !! Store auxiliar variables for the diagonalization
  
  allocate(q_band(1:fedvr3d%nb_r,1:fedvr3d%nb_r))
  allocate(w_band(1:fedvr3d%nb_r))
  allocate(z_band(1:fedvr3d%nb_r,1:fedvr3d%nb_r))
  allocate(work_band(1:fedvr3d%nb_r*7))
  allocate(iwork_band(1:fedvr3d%nb_r*5))
  allocate(ifail_band(1:fedvr3d%nb_r))
  
  !! Initialize the auxiliar orbitals and energies

  norbital=size(orb(1,:)) !! number of orbitals
  nbasis=size(orb(:,1)) !! number of basis functions

  allocate(orb_aux(1:nbasis,1:norbital)) !! store as the orbitals.
  allocate(orb_ener(1:norbital)) !! energies for each orbital.   

  orb_aux=zzero
 
  write(*,*)   
  write(*,*) 'Construct the Hamiltonian applying stiffness' 
  write(*,*) '--------------------------------------------'
  write(*,*)   
  
  Do il=0,fedvr3d%l_max
       htemp=zero
     write(*,*) 'Building for L=',il
     write(*,*) ''
       
     !! initialize the auxiliar variables of the diagonalization to 0.
     q_band=zero
     w_band=zero
     work_band=zero
     iwork_band=0
     ifail_band=0
     Do it=1,size(t_r(:))
!! First, we store the band matrix 'htemp' to diagonalize.
!! The kinetic terms
        i_column=index_t(it,1) !! is taken in this way because DSBEVX diagonalize is written to diagonalize upper matrices
        i_row=index_t(it,2)          
        If (i_row.le.i_column) then 
           htemp(kd+1+i_row-i_column,i_column)=t_r(it)
        End If
!! The terms in the diagonal include the effective potential
        If (i_row.eq.i_column) then
           htemp(kd+1+i_row-i_column,i_column)=htemp(kd+1+i_row-i_column,i_column)-z/xglobal(i_column)+dble(il*(il+1))/(2.0d0*xglobal(i_column)**2.0d0)
        End If
        !! we store it in the hl matrix
        !! for the stiffness we substract the contribution of the
        !! eigenvectors with E> E_cutoff
        hl(i_row,i_column,il)=htemp(kd+1+i_row-i_column,i_column)
        hl(i_column,i_row,il)=htemp(kd+1+i_row-i_column,i_column)
     End Do
       
     !! The band matrix 'h_temp' is stored. Now, we diagonalize it.
     
     !! Diagonalize
     
     z_band=zzero
     w_band=zzero

     call DSBEVX( 'V', 'V', 'U', fedvr3d%nb_r, kd, htemp,&
          kd+1, q_band, fedvr3d%nb_r, energy_min, 100*energy_cutoff, 1,&
          2, 1d-20, m_band, w_band, z_band, fedvr3d%nb_r,&
          work_band, iwork_band, ifail_band, info)
     
     !! Now we construct substract the eigenenergies which are not included
     !! we substract from hl -\sum_j E_j |j><j|, with E_j>E_cutoff

     Do i_eigen=fedvr3d%nb_r,1,-1 
        if (w_band(i_eigen).lt.energy_cutoff) then
           print*, fedvr3d%nb_r-i_eigen, 'states removed for L=', il
           exit !! E_{i_eigen}<E_cutoff, exit
        end if
        Do i_column=1,fedvr3d%nb_r
           Do i_row=1,fedvr3d%nb_r
              hl(i_row,i_column,il)=hl(i_row,i_column,il)-w_band(i_eigen)*z_band(i_row,i_eigen)*z_band(i_column,i_eigen)
           End Do
        End Do
     End Do

     !! We obtain here the initial orbitals. They are the ones with
     !! lower eigenenergies.

     !! First, we store the ones with L=0
     If (il.eq.0) then !! IF L=0.
        Do ij=1, norbital !! here we initialize the energies, and store the contribution of L=0
           orb_ener(ij)=w_band(ij) !! storing the energies
           Do jk=1,fedvr3d%nb_r !! storing the eigenvalues
              orb_aux(jk,ij)=cmplx(z_band(jk,ij))
           End Do
        End Do
     End If

     If (il.gt.0) then !! IF L> 0
        !! In this case, we find degeneration in M, that we have to take
        !! into account.

        orb_ener_order=0 !! orb_ener_order is the position of the last energy substitued with a fixed 'il'

        Do ij=1,norbital-1 !! running in orb_ener energies !! CHECK 16.04.2015 (before ij=1,norbital)             
           If (ij.le.orb_ener_order) cycle !! If this energy has already been stored-> cycle
              Do ik=1,norbital !! running in the w_band energies
                 If (abs(w_band(ik)-orb_ener(ij)).ge.1.0d-10.and.w_band(ik).lt.orb_ener(ij+1)) then  !! Energy between two stored energies
                    !! We calculate the multiplicity of this factor
                    multiplicity=0
                    Do im=-min(il,fedvr3d%m_max),min(il,fedvr3d%m_max)
                       !! sum in allowed m's, which are degenerated.
                       multiplicity=multiplicity+1
                    End Do
                    
                    orb_ener_order=min(ij+1+multiplicity,norbital) !! last energy we are storing according to the multiplicity and the states which are left.
                    

                    Do jk=ij+1,orb_ener_order !! We leave free the positions to store the energies and eigenfunctions
                       If (jk+orb_ener_order.gt.norbital) exit 
                       orb_ener(orb_ener_order+jk)=orb_ener(jk)
                       orb_aux(:,orb_ener_order+jk)=orb_aux(:,jk)
                    End Do

                    !! To store the orbitals, we look for the first global spatial functions with angular momentum 'il'
                    counter=0
                    Do ilprime=0,il-1 !! Run in all the functions with L<= il
                       Do im=-min(fedvr3d%m_max,ilprime),min(fedvr3d%m_max,ilprime)
                          counter=counter+1
                       End Do
                    End Do


                       !! store energies                    
                       Do jk=ij+1,orb_ener_order
                          If (jk+orb_ener_order.gt.norbital) exit 
                          orb_ener(jk)=w_band(ik)                          
                          
                          !! Now, we run for these functions with angular momentum 'il'
                          !! fedvr3d%nb_r*(counter-1)+1 is the first function with the 
                          orb_aux(:,jk)=zzero !!initialize

                          Do jk1=1,fedvr3d%nb_r                      
                             orb_aux(fedvr3d%nb_r*(counter-1+(jk-ij-1))+jk1,jk)=cmplx(z_band(jk1,jk))
                             !! we add (jk-ij-1) to run in the angular functions
                          End Do
                       End Do
                       
                       cycle !! cycle to the following energy in orb_ener
                    End If
                 End Do !! loop in ik
              End Do !! loop in ij
           End If

        End Do !! loop in il

  write(*,*) 'Finished the construction of the Hamiltonian'
  write(*,*) 'applying stiffness'
  write(*,*) '============================================'
  write(*,*)
 
!! store the orbitals in the output
orb=zzero 
orb=orb_aux

  !! deallocate auxiliar variables for diagonalization

  deallocate(work_band,q_band,iwork_band,ifail_band)
  deallocate(w_band,z_band)
  deallocate(htemp,orb_aux)

!! Normalize the result

  Do ij=1,size(orb(1,:))
     orb(:,ij)=orb(:,ij)/sqrt(dot_product(orb(:,ij),orb(:,ij)))
  End Do

  return
End subroutine hl_stif


end module wfunction

!
! drive module wfunction
!
subroutine drive_wfunction
 use global
 use wfunction
 implicit none

  call wf_beginning

!
! if hf core, load the hf orbital as the core orbitals
!

  if(system%io_fcore==1) then

  endif

 return
end subroutine drive_wfunction
 
