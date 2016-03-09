
module wfunction
  use global
  use auxiliary
  use operator_3d

  implicit none
  complex(kind=k2),allocatable,save :: acoff(:),phi(:,:)
  complex(kind=k2),allocatable,save :: acoff0(:),phi0(:,:)
  complex(kind=k2),allocatable,save :: acoff_old(:),phi_old(:,:)
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
     real(kind=k1) ::  rstep_smallest,  rstep_aux
     

     allocate(acoff(n_string_ab), phi(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) )

     allocate(acoff0(n_string_ab), phi0(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) )
     
     if(analysis%i_auto==1) then
       allocate(acoff_old(n_string_ab), phi_old(fedvr3d%nb_r*fedvr3d%nb_angle,&
       system%nptot) )
     endif
!
! if real time propagation
!

if(laser%tdornot) then

    prop%chstep = dcmplx(prop%hstep,0.0d0)

        inquire(file='imag.txt',exist=lexist)
        if (.not.lexist) then
          write(*,*) '---real time propagation, but the ground state wf is not exist---'
          stop 
       endif
!
! readin the acoff from the imag.txt file    
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

!===============================================================================
! imaginary time propgation 
!===============================================================================

else

  prop%chstep = dcmplx(0.0d0, -1.0_k1*prop%hstep) !! set the imaginary step

  open(ifile,file='imag.txt',status='unknown')

 if(prop%stiffness==1) then
  allocate(w_temp(fedvr3d%nb_r*fedvr3d%nb_angle),fv1(fedvr3d%nb_r*fedvr3d%nb_angle),&
  fv2(fedvr3d%nb_r*fedvr3d%nb_angle))

  allocate(z_temp(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle),&
  h_temp(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle)) 

  allocate(h_stiffness(fedvr3d%nb_r*(fedvr3d%nb_r+1)/2,fedvr3d%nb_angle))
 endif


!
! imaginary propagation
! prepare the initial A cofficient
!

   do idicp =1, n_string_ab
     acoff(idicp) = dcmplx(sqrt( 1.0d0/dble(n_string_ab) ),0.0d0) !! Commented by Juan
!!$     acoff(idicp) = zzero !! Added by Juan
   enddo
!!$   acoff(1)=dcmplx(1.0d0,0.0d0) !! Added by Juan

!
! initial the spatial orbitals
!

select case (prop%stiffness)

   case(0)

     write(*,*)
     write(*,*) 'the stiffness effect are not considered in the running '

     allocate(phi_state(1:fedvr3d%nb_r*fedvr3d%nb_angle))
     counter=0

Do n=1,system%nptot !! run in principal quantum number
   Do lh=0,min(n-1,fedvr3d%l_max),1 !! run in the total angular momentum
      Do mh=-min(fedvr3d%m_max,lh),min(fedvr3d%m_max,lh),1 !! run in the magnetic quantum number
         counter=counter+1
         call radial_hydrogen_dvr(fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m,fedvr_w,n,lh,mh,system%zz,&
phi_state)
         phi(:,counter)=phi_state
         If (counter.eq.system%nptot) exit
      End Do
      If (counter.eq.system%nptot) exit
   End Do
   If (counter.eq.system%nptot) exit
End Do

deallocate(phi_state)

    
!!     step_biggest =  1.0d0/w_temp(fedvr3d%nb_r*fedvr3d%nb_angle)*3.1415926d0*0.5d0

rstep_smallest=fedvrx_global(size(fedvrx_global))

Do ii=1,size(fedvrx_global)-1
    rstep_aux=fedvrx_global(ii+1)-fedvrx_global(ii) 
    !! set the step in the radial part for the node in ii and ii-1
    If (rstep_aux.lt.rstep_smallest) rstep_smallest=rstep_aux
End Do

!!$step_biggest=(.25d0+dble(fedvr3d%l_max*(fedvr3d%l_max+1)))/(rstep_smallest)**2.0d0
!!$
!!$step_biggest=1.0d0/step_biggest !!Juan: I think this is a bad choise

step_biggest=(rstep_smallest/3.0d0)**2.d0

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

   write(*,*) 'Stiffness is taken into account'

    h_temp = 0.0d0


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

!! DIAGONALIZATION

    call DSBEVX( 'V', 'V', 'U', fedvr3d%nb_r*fedvr3d%nb_angle, kd, h_band,&
         kd+1, q_matrix, fedvr3d%nb_r*fedvr3d%nb_angle, -200, 100, 1,&
         2, 1d-12, , w_band, z_band, fedvr3d%nb_r*fedvr3d%nb_angle,&
         work_band, iwork_band, ifail_band, info )

    print*, 'module_wf'
    stop

!! END of Building the matrix to diagonalize
!!=====================================================================

! diag. the full matrix matrix

      call rs(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle,h_temp,& 
      w_temp,1,z_temp,fv1,fv2,ierr)

      if(ierr/=0) then
        write(*,*) 'error in subroutine initial_wp'
        stop 
      endif

!
! stiffness one body operatror
!     
!     E_biggest = 1.0d0/prop%hstep  ! Juan
     E_biggest = pi*0.5d0/prop%hstep! Juan

     do idicp =1,fedvr3d%nb_r*fedvr3d%nb_angle   
       If( w_temp(idicp).gt.E_biggest ) then
          w_temp(idicp) = zero
       End If
     enddo 

     do idicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! number of the state in the row
      do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! number of the state in the column
       sumtemp = zero
        do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle !! Number of eigenvalues
          sumtemp = sumtemp + z_temp(idicp,kdicp)*w_temp(kdicp)*z_temp(jdicp,kdicp)
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

 ! diag. the matrix
     call rs(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle,h_temp,&
      w_temp,1,z_temp,fv1,fv2,ierr)

!!$
!
! diag. hcore get spatial orbital
!     

      do idicp =1,system%nptot
       do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
        phi(jdicp,idicp) = dcmplx(z_temp(jdicp,idicp),0.d0) 
       enddo
     enddo  

     deallocate(h_temp,z_temp,w_temp,fv1,fv2)

  End select   


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

!! Function

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
 
