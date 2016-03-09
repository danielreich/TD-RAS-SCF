!===============================================================================
! calc. the matrix element of kinetic energy and nuclei attractive
! potential, dipole matrix in fedvr_3d, spherial coordiante 
!===============================================================================
 module operator_spatial_space
  use global
  use operator_3d
  use wfunction
  use twoe_basis_set
  use fedvr3d_basis_set

  implicit none  
  complex(kind=k2),allocatable,save :: ham_onebody(:),xmat_spatial_space(:),&
  ymat_spatial_space(:),zmat_spatial_space(:),dvdxmat_spatial_space(:),&
  dvdymat_spatial_space(:),dvdzmat_spatial_space(:)
  complex(kind=k2),allocatable,save :: ch_dummy(:,:)
  complex(kind=k2),allocatable,save :: tei_spatial(:,:,:,:),cww(:,:,:,:)
  complex(kind=k2), allocatable, save :: mean_field_lm(:,:,:)
!! mean_field_lm is the Mean Field Operator written in the coupled basis.
!! Arguments
!!
!! #1 Coupled basis function
!! #2 Orbital which is complex conjugated (coming from the bra).
!! #3 Orbital which is not complex conjujated (coming from the ket).
!!
  
  contains
!
! calc. the onebody operator in the spatial space
!
   subroutine update_ham_onebody_spatial_space()
    implicit none
    integer :: i_spatial,j_spatial,k_temp
    integer :: i_angle,i_r
    integer :: i_row_r,i_column_r,i_row,i_column,index_here
    integer :: idicp,jdicp,kdicp,idicp_prime,jdicp_prime
    real(kind=k1) ::h_onebody
    complex(kind=k2):: sumtemp,sumtemp_laser_term 
    integer:: iangle,ibra,iket,iorb !! indexes for the loops.
    integer :: iangle_bra, iangle_ket !! indexes for the loops in angular functions
    integer:: norbital, nbasis !! number of orbitals and basis fuctions.
    integer:: nbasis_r !! number of radial fucntions
    integer:: l_basis !!angular momentum
    real(kind=k1) :: ang_x, ang_y, ang_z !! \cos\theta,\sin\theta\cos\phi, and \sin\theta\sin\phi contributions
    !! for the differents functions.
    complex(kind=k2),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: korb

interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

nbasis_r=fedvr3d%nb_r !! number of radial functions

!    allocate(ham_onebody(system%nptot*(system%nptot+1)/2))

    k_temp = 0
    do i_spatial =1, system%nptot
     do j_spatial =1,i_spatial

       k_temp = k_temp +1

!
! the total number of non-zero matrix element of ham_onebody
!        
       sumtemp = zzero

    do i_angle =1, fedvr3d%nb_angle        
     do i_r = 1,n_total_kinetic
          
      i_row_r = index_kinetic_basis(i_r,1)
      i_column_r = index_kinetic_basis(i_r,2)


      i_row = (i_angle-1)*fedvr3d%nb_r + i_row_r
      i_column = (i_angle-1)*fedvr3d%nb_r + i_column_r

!      index_here = ia(max(i_row_r,i_column_r) + min(i_row_r,i_column_r))


!      if(i_row == i_column) then
!        h_onebody = tmat_3d(i_row_r,i_angle) + vmat_radial(i_row_r) + tmat_radial(index_here)
!      endif

!      if(i_row /= i_column) then
!        h_onebody = tmat_radial(index_here)
!      endif
!      sumtemp = sumtemp + dconjg(phi(i_row,i_spatial))*phi(i_column,j_spatial)*h_onebody   
  
      if(i_row==i_column) then

       sumtemp = sumtemp + dconjg(phi(i_row,i_spatial))*phi(i_column,j_spatial)*&
       (tmat_3d(i_row_r,i_angle) +vmat_radial(i_row_r) + tmat_radial(i_r) )

      endif   
  
      if(i_row /= i_column) then
        sumtemp = sumtemp + dconjg(phi(i_row,i_spatial))*phi(i_column,j_spatial)* tmat_radial(i_r) + &
        dconjg(phi(i_column,i_spatial))*phi(i_row,j_spatial)* tmat_radial(i_r) 
      endif     
   

   enddo
enddo


!=========================================================================================================
!     laser field term, add here
!========================================================================================================= 
 
     sumtemp_laser_term = zzero

If (laser%tdornot) then

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(iangle_bra,ang_x,ang_y,ang_z,iket)
!$OMP DO SCHEDULE(dynamic)
   
      Do iangle_ket=1,size(lm_l(:))!! run in the angular functions of the ket
         Do iangle_bra=1,size(lm_l(:)) !! run in the angular functions of the bra
            ang_x=xmat_3d(iangle_bra,iangle_ket) !! angular contribution of x direction
            ang_y=ymat_3d(iangle_bra,iangle_ket) !! angular contribution of y direction
            ang_z=zmat_3d(iangle_bra,iangle_ket) !! angular contribution of z direction
            Do iket=1,nbasis_r !! run in the radial basis               
               korb(nbasis_r*(iangle_bra-1)+iket,j_spatial)=korb(nbasis_r*(iangle_bra-1)+iket,j_spatial)+phi(nbasis_r*(iangle_ket-1)+iket,j_spatial)*fedvrx_global(iket)*(laser%ex_t*ang_x+laser%ey_t*ang_y*ci+laser%ez_t*ang_z)               
               !! multiply by the imaginary unity the contribution in the Y direction.
               !! Add the absorbing potential
               
               If (iangle_ket.eq.iangle_bra) then
                  korb(nbasis_r*(iangle_ket-1)+iket,j_spatial)=korb(nbasis_r*(iangle_ket-1)+iket,j_spatial)+phi(nbasis_r*(iangle_ket-1)+iket,j_spatial)*vabsorb_pot_fedvr3d(iket)*ci !! we add the imaginary unit here, which is part of the absorbing potential.
               End If

            End Do !! iket
            !! initialize the variables
            ang_x=zero
            ang_y=zero
            ang_z=zero
         End Do !! iangle_bra
      End Do !! iangle_ket


!$OMP END DO
!$OMP END PARALLEL

sumtemp_laser_term=dot_product(phi(:,i_spatial),korb(:,j_spatial)) !! laser and absorbing potential contribution

End If

ham_onebody(k_temp) = sumtemp +  sumtemp_laser_term

enddo
enddo

return
end subroutine update_ham_onebody_spatial_space

!
! stiffness onebody operator
!

   subroutine update_ham_onebody_spatial_space2()
    implicit none
    integer :: i_spatial,j_spatial
    integer :: idicp,jdicp,kdicp,ldicp,idicp_prime,jdicp_prime   
    complex(kind=k2) :: sumtemp,sumtemp_laser_term    

   do i_spatial =1,system%nptot
     do j_spatial =1,i_spatial

    sumtemp = zzero
    do idicp =1,fedvr3d%nb_angle
      ldicp = 0
      do jdicp =(idicp-1)*fedvr3d%nb_r+1,(idicp-1)*fedvr3d%nb_r+fedvr3d%nb_r
        do kdicp =(idicp-1)*fedvr3d%nb_r+1,jdicp
        
         ldicp = ldicp +1
         if(jdicp/=kdicp) then
          sumtemp = sumtemp + (dconjg(phi(jdicp,i_spatial))*phi(kdicp,j_spatial) + &
          dconjg(phi(kdicp,i_spatial))*phi(jdicp,j_spatial))*h_stiffness(ldicp,idicp) 
         else
          sumtemp = sumtemp + dconjg(phi(jdicp,i_spatial))*phi(kdicp,j_spatial)*h_stiffness(ldicp,idicp)
         endif
      enddo
     enddo
    enddo
!=========================================================================================================
!     laser field term, add here
!========================================================================================================= 

     sumtemp_laser_term = zzero
     if(laser%tdornot) then

        do idicp =1,fedvr3d%nb_angle
          do jdicp =1, fedvr3d%nb_angle

            do kdicp =1,fedvr3d%nb_r

              idicp_prime  =  (idicp-1)*fedvr3d%nb_r + kdicp
              jdicp_prime  =  (jdicp-1)*fedvr3d%nb_r + kdicp
              sumtemp_laser_term = sumtemp_laser_term + dconjg(phi(idicp_prime,i_spatial ))* phi(jdicp_prime,j_spatial)*&
              (zmat_3d(idicp,jdicp)*fedvrx_global(kdicp)*laser%ez_t)

              If (idicp_prime.eq.jdicp_prime) then
                 sumtemp_laser_term=sumtemp_laser_term+ci*vabsorb_pot_fedvr3d(kdicp)*dconjg(phi(idicp_prime,i_spatial ))* phi(jdicp_prime,j_spatial)
              End If


!!$              print*, 'vabsorb_pot_fedvr3d is not well used here -stif-'
!!$              print*, 'module_operator_spatial_space.f90'
!!$              stop

             ! sumtemp_laser_term = sumtemp_laser_term + dconjg(phi(idicp_prime,i_spatial ))* phi(jdicp_prime,j_spatial)*&
             ! xmat_3d(idicp,jdicp)*fedvrx_global(kdicp)*laser%ex_t


             ! sumtemp_laser_term = sumtemp_laser_term + dconjg(phi(idicp_prime,i_spatial ))* phi(jdicp_prime,j_spatial)*&
             ! ymat_3d(idicp,jdicp)*fedvrx_global(kdicp)*ci*laser%ey_t

            enddo
          enddo
        enddo
     endif

 ham_onebody(ia(i_spatial)+j_spatial ) = sumtemp + sumtemp_laser_term
enddo
enddo
 

return
end subroutine update_ham_onebody_spatial_space2

!
! stiffness onebody operator !! Done by Juan
!

   subroutine update_ham_onebody_spatial_space_stif_2(orb,hl,laserp,ham1body)
    implicit none
    integer:: i_spatial, j_spatial,counter

    !! INPUT
    
    complex(kind=k2), allocatable, intent(in) :: orb(:,:)
    !! orb stores the orbitals
    !! Arguments
    !! 
    !! #1: component of the global basis (FE-DVR+Angular).
    !! #2: number of the orbital.

    real(kind=k1), allocatable, intent(in) :: hl(:,:,:)
    !! Hamiltonian after the stiffness procedure.
    !! The hamiltonian matrix is the same for each value of the 
    !! angular momentum l. 
    !! Arguments
    !!
    !! #1: FE-DVR function of the bra.
    !! #2: FE-DVR function of the ket.
    !! #3: Angular momentum of the Hamiltonian, l.
    !!

    type(laser_prime), intent(in) :: laserp
    !! laserp stores the properties of the laser

    
    !! OUTPUT
    
    complex(kind=k2), allocatable, intent(inout) :: ham1body(:)
    !! ham1body stores the one body hamiltonian <phi_i|H|phi_j>. 
    !! <phi_i|H|phi_j>=ham1body(i(i-1)/2+j), with j<=i

    !! AUXILIAR
    
    complex(kind=k2), dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: korb_aux
    !! korb_aux stores the orbitals after applying the one body 
    !! Hamiltonian.
    !! Arguments
    !! 
    !! #1: component of the global basis (FE-DVR+Angular).
    !! #2: number of the orbital.

    !! First, apply the one body Hamiltonian on the orbitals
    
    call act_stif_2(orb,hl, laserp,korb_aux)    

    ham1body=zzero !! initialize

    do i_spatial =1,system%nptot
       do j_spatial =1,i_spatial
          counter=i_spatial*(i_spatial-1)/2+j_spatial
          !! counter stores half of the matrix          
          ham1body(counter)=ham1body(counter)+dot_product(orb(:,i_spatial),korb_aux(:,j_spatial))
       enddo
    enddo

!!$    if (allocated(korb_aux)) deallocate(korb_aux) !! DEBUG 11/5/2015

return
end subroutine update_ham_onebody_spatial_space_stif_2

!
! twobody operator in spatial space
!
!
!
! calc. the CWW matrix, mean field operator
!
!
!
! calc. the mean field operator matrix in fedvr-3d
!
    subroutine update_twobody_spatial_space
     implicit none
     integer :: imo,jmo,kmo,lmo,i_fedvr3d,j_fedvr3d,k_fedvr3d,l_fedvr3d
     integer :: n1,l1,m1,n2,l2,m2,n3,l3,m3,n4,l4,m4
     complex(kind=k2) :: sumtemp
     real(kind=k1) :: rtemp
     integer :: idicp,jdicp,itwo,kdicp,ldicp,mdicp,ndicp,i_global,j_global      
     real(kind=k1) :: start,finish     


     select case (fedvr3d%store)

     case (0)  !! do not store cww

        call  two_body_orbitals(tei_spatial)

     case (1)  !! stores cww

     cww = zzero

     do jdicp =1, fedvr3d%nb_angle
        do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
           do itwo = 1, num_two(idicp,jdicp)
              
              kdicp = index_two_storage(1,itwo,idicp,jdicp)
              ldicp = index_two_storage(2,itwo,idicp,jdicp)
              
              do ndicp =1, system%nptot
                 do mdicp = 1, system%nptot
                    cww(idicp,jdicp,mdicp,ndicp) = cww(idicp,jdicp,mdicp,ndicp) +  &
                         dconjg(phi(kdicp,mdicp)) * phi(ldicp,ndicp) *two_storage(itwo,idicp,jdicp) 
                 enddo
              enddo
           enddo
        enddo
     enddo

! cal. two-electron repulsive in spatial space
!
!       *            *
!  fi(1)  fj(1) fk(2)  fl(2)
!
     
     do jmo =1, system%nptot
        do lmo = 1,system%nptot
           do kmo =1, system%nptot
              do imo =1, system%nptot
                 
                 sumtemp = zzero
                 do j_fedvr3d =1,fedvr3d%nb_angle
                    do i_fedvr3d =1,fedvr3d%nb_r*fedvr3d%nb_angle
                       
                       
                       i_global = index_two(1,i_fedvr3d,j_fedvr3d)  
                       j_global = index_two(2,i_fedvr3d,j_fedvr3d)
                       
                       sumtemp = sumtemp + dconjg(phi(i_global,imo))*phi(j_global,jmo)*cww(i_fedvr3d,j_fedvr3d,kmo,lmo)           
                    enddo
                 enddo
                 tei_spatial(imo,kmo,lmo,jmo) = sumtemp  
              enddo
           enddo
        enddo
     enddo

      continue 

      case(2) !! coupled case

         !! Transforms the product of two orbitals to the coupled representation
         call mean_field_and_phikl(twoe_radial_store, phi, llmm_lm, mean_field_lm, phi_coupled)
         !! calculation of the two-body element
         call two_body_coupled(mean_field_lm,phi_coupled,tei_spatial)         
         
         continue

      End select
            
      return
    end subroutine update_twobody_spatial_space

!
! update xyzmat_spatial_space
!
subroutine update_xyzmat_spatial_space(xyzmat_3d,ndim,xyzmat_spatial,mdim)
  implicit none
  integer :: i_spatial,j_spatial,ktemp,i_row,i_column
  integer,intent(in) :: ndim,mdim
  real(kind=k1),intent(in) :: xyzmat_3d(ndim,ndim)
  complex(kind=k2),intent(out) :: xyzmat_spatial(mdim*(mdim+1)/2)
  complex(kind=k2) :: sumtemp
  integer :: i_angle,j_angle,k_r
  ktemp =0
  do i_spatial =1, mdim  ! nptot
    do j_spatial =1,i_spatial
      ktemp = ktemp +1
      sumtemp =zzero

      do i_angle =1, ndim ! fedvr3d%nb_angle
        do j_angle =1, ndim ! fedvr3d%nb_angle

         if(xyzmat_3d(i_angle,j_angle)/=0.0d0) then !!??? !! This line is not correct!!! : JUAN found mistake

          do k_r =1,fedvr3d%nb_r

            i_row = (i_angle-1)*ndim + k_r
            i_column = (j_angle-1)*ndim + k_r 
  
            sumtemp = sumtemp + dconjg(phi(i_row,i_spatial))*phi(i_column,j_spatial)*fedvrx_global(k_r)*&
            xyzmat_3d(i_angle,j_angle)   
          enddo 
         endif
        enddo
       enddo
     xyzmat_spatial(ktemp) =  sumtemp  
    enddo
 enddo


 return
end subroutine update_xyzmat_spatial_space


!
! update dvdxyz_spatial_space for hhg
!
 subroutine update_dvdxyzmat_spatial_space(xyzmat_3d,ndim,dvxyzmat_spatial,mdim)
  implicit none
  integer :: i_spatial,j_spatial,ktemp,i_row,i_column
  complex(kind=k2) :: sumtemp
  integer,intent(in) :: ndim,mdim
  real(kind=k1),intent(in) :: xyzmat_3d(ndim,ndim)
  complex(kind=k2),intent(out) :: dvxyzmat_spatial(mdim*(mdim+1)/2)
  integer :: i_angle,j_angle,k_r


  ktemp =0
  do i_spatial =1, mdim !system%nptot
    do j_spatial =1,i_spatial
      ktemp = ktemp +1
      sumtemp =zzero

      do i_angle =1,ndim !fedvr3d%nb_angle
       do j_angle =1,ndim !fedvr3d%nb_angle

         if(xmat_3d(i_angle,j_angle)/=0.0d0) then 

          do k_r =1,fedvr3d%nb_r

            i_row = (i_angle-1)*ndim + k_r
            i_column = (j_angle-1)*ndim + k_r            
            sumtemp = sumtemp + dconjg(phi(i_row,i_spatial))*phi(i_column,j_spatial)*&
            xyzmat_3d(i_angle,j_angle)   
          enddo 
         endif
        enddo
       enddo
     dvxyzmat_spatial(ktemp) =  sumtemp  
    enddo
   enddo
 
  return
 end subroutine update_dvdxyzmat_spatial_space

!! This subroutines are done by Juan to avoid the storing of the two body operators

!! Calculate the average the due to the two electrons interaction integrating in one coordinate for two orbitals ph1 and phi2

Function cww_calc(chi1,chi2_ang,phi1,phi2)    
implicit none
integer :: chi1, chi2, chi2_ang !! FEDVR functions
!! chi1 is the global FEDVR of the first DVR.
!! chi2 is the angular part of the second DVR function.
integer :: phi1, phi2 !! orbitals
complex(kind=k2) :: cww_calc
!!Auxiliar variables
integer :: kdicp, ldicp, ldicp_ang,kdicp_ang, kdicp_radial
integer :: kdicp_end,ldicp_ang_end
integer :: n1, l1,m1,n2,l2,m2,n3,l3,m3,l4,m4
integer :: k, kp,ll, ang1,ang2, lmin,lmax,ang3,ang4, ii
real(kind=k1) :: twoe_fedvr_aux, twoe_fedvr_angle,twoe_fedvr_radial

!!$allocate(a(1:fedvr3d%nb_r,1:fedvr3d%nb_r,1:2*fedvr3d%nb_angle))
!!$print*, size(a)
!!$stop

If (chi1.gt.fedvr3d%nb_r*fedvr3d%nb_angle) then
   print*, 'ERROR IN cww_calc'
   stop
end If

If (chi2_ang.gt.fedvr3d%nb_angle) then
   print*, 'ERROR IN cww_calc'
   stop
end If

      !! For the FEDVR function
      n1=global_to_local(chi1,1)
      l1=global_to_local(chi1,2)
      m1=global_to_local(chi1,3)

      chi2=(chi2_ang-1)*fedvr3d%nb_r+n1

      n2=n1
      l2=global_to_local(chi2,2)
      m2=global_to_local(chi2,3)

cww_calc=cmplx(0.0d0,0.0d0)                                            

If (phi1.gt.system%np0+system%np1) then !! Conditions for the core
   kdicp_end=fedvr3d%nb_r*fedvr3d%nb_angle
   Else
   kdicp_end=fedvr3d%nb_r*ang_max_core
End If

If (phi2.gt.system%np0+system%np1) then !! Conditions for the core
   ldicp_ang_end=fedvr3d%nb_angle
   Else
   ldicp_ang_end=ang_max_core
End If

ang1=(chi1-n1)/fedvr3d%nb_r+1 !angular for the first function
ang2=chi2_ang         ! angular for the second function

Do kdicp=1, kdicp_end!fedvr3d%nb_r*fedvr3d%nb_angle ! Juan: 25 November

   If (abs(phi(kdicp,phi1)).lt.1d-10) cycle           

   n3=global_to_local(kdicp,1) !!Juan: 25 november  
   kdicp_radial=n3
   kdicp_ang=(kdicp-n3)/fedvr3d%nb_r+1
      
      l3=lm_l(kdicp_ang) 
      m3=lm_m(kdicp_ang)

   Do  ldicp_ang=1,ldicp_ang_end!fedvr3d%nb_angle !! Run in the angular FEDVR basis for the angular part (since the radial is the same).      

      ldicp=(ldicp_ang-1)*fedvr3d%nb_r+n3

      If (abs(phi(ldicp,phi2)).lt.1d-10) cycle   
      If (abs(conjg(phi(kdicp,phi1))*phi(ldicp,phi2)).lt.1d-15) cycle   

      l4=lm_l(ldicp_ang)
      m4=lm_m(ldicp_ang)

      If ((m1-m2).ne.(m4-m3)) cycle 
      If (max(abs(l1-l2), abs(l3-l4)).gt.min((l1+l2), (l3+l4))) cycle


      Do ll=max(abs(l1-l2), abs(l3-l4)), min((l1+l2), (l3+l4))
         
         twoe_fedvr_angle=fedvr3dbase_angpart(ll,l1,m1,l2,m2,l3,m3,l4,m4)!!twoe_angle_store(ang1,chi2_ang,kdicp_ang,ldicp_ang,ll)
         
!!$         If (abs(twoe_fedvr_angle).lt.1d-15) cycle
         
!!$         twoe_fedvr_aux=zero
         
         twoe_fedvr_aux=twoe_radial_store(n1,n3,ll)*twoe_fedvr_angle !!twoe_angle_store(ang1,ang2,ang3,ang4,ll) !! Juan: 25 November
         cww_calc = cww_calc + conjg(phi(kdicp,phi1))*phi(ldicp,phi2)*twoe_fedvr_aux !! Juan: I have included this (November 24th)     
         
      End Do
      
   Enddo
   
End Do


return
End Function cww_calc


Subroutine two_body_orbitals(tei_spatial)
implicit none
complex(kind=k2),allocatable, intent(inout) :: tei_spatial(:,:,:,:) !! orbitals
!!Auxiliary variable
 complex(kind=k2) :: sumtemp
integer :: imo,jmo,kmo,lmo,i_fedvr3d,j_fedvr3d
integer :: i_global, j_global
complex(kind=k2),allocatable,save :: tei_spatial_proc(:,:,:,:,:)
integer :: num_proc,proc !! number of processes
interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

If (allocated(tei_spatial)) then
   continue
   else
 !! If tei_spatial is not allocated, do it!
   allocate(tei_spatial(1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))
End If

tei_spatial=zzero !! initialize the two body integrals

!! AUXILIARY VARIABLES FOR PARALLELIZATION

!! SET MAX NUM OF THREADS
!!$      call OMP_set_num_threads(4)

      proc=OMP_get_num_procs()

if (allocated(tei_spatial_proc)) then
   continue
else
   allocate(tei_spatial_proc(0:proc-1,1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))   !! to parallelize
   !! this is the contribution of each processor
End if

tei_spatial_proc=zzero    

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(j_fedvr3d,i_global,j_global,imo,jmo,kmo,lmo) 
!$OMP DO SCHEDULE(dynamic)!, 2200)

do i_fedvr3d =1,fedvr3d%nb_r*fedvr3d%nb_angle
   do j_fedvr3d =1,fedvr3d%nb_angle
      
      i_global = i_fedvr3d
      j_global = (j_fedvr3d-1)*fedvr3d%nb_r+global_to_local(i_fedvr3d,1)
            
      do imo =1, system%nptot !! Run in orbitals

         If (abs(phi(i_global,imo)).lt.1d-10) cycle !! Juan: change 14 november AQUI

         do jmo =1, system%nptot !! Run in orbitals
            If (abs(conjg(phi(i_global,imo))*phi(j_global,jmo)).lt.1d-10) cycle !! Juan: change 13 november AQUI

            do kmo =1, system%nptot !! Run in orbitals      
               do lmo = 1,system%nptot !! Run in orbitals
                  
                   tei_spatial_proc(OMP_get_thread_num(),imo,kmo,lmo,jmo)=&
                                tei_spatial_proc(OMP_get_thread_num(),imo,kmo,lmo,jmo)+&
                                dconjg(phi(i_global,imo))*cww_calc(i_fedvr3d,j_fedvr3d,kmo,lmo)*phi(j_global,jmo)     
                   
               enddo
            enddo
            
         enddo
      enddo
   enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

            Do num_proc=0,proc-1
               !! sum all the contributions
               tei_spatial(:,:,:,:)=tei_spatial(:,:,:,:)+tei_spatial_proc(num_proc,:,:,:,:)
            End Do
    
           deallocate(tei_spatial_proc)

return

End Subroutine two_body_orbitals

!! End of the subroutines to avoid the storing of the two body matrix

!! SUBROUTINES TO BUILD THE MEAN FIELD OPERATOR USING THE FOURIER TRANSFORM.

Subroutine al_make(k_max,fedvr_r,l_max,al,fl)  !!CHECK 
  implicit none

!! INPUT
!! Maximum linear momentum considered
  real (kind=k1) :: k_max
!! Nodes of the FE-DVR for the position space.
  real (kind=k1), intent(in), allocatable :: fedvr_r(:)
!! Maximum value of the orbital angular momentum.
  integer, intent(in) :: l_max

!! OUTPUT
!! Operator A_L, which performs the Fourier Transform for fixed angular momentum
!! L of reduced wavefunctions (wavefunctions multiplied by the radial coordinate)
!! The Fourier Transform goes as F\psi= \sqrt(2/\pi)\times (-i)**L A_L\psi

!! Arguments
!! #1 is the row of the matrix. It changes with the momentum
!! #2 is the row of the matrix. It changes with the position
!! #3 labels the angular momentum of each matrix.

  real (kind=k1), intent(out), allocatable :: al(:,:,:)
  
!! fl is the matrix which performs the operation to obtain the  mean field operator. It collects the operations F(-1)(W(k)F(phi*phi)). In matricial form, this operation can be approximated as f_l=a_l**(inverse) W(k) a_l = a_l**(transpose)*a_l.
  real (kind=k1), intent(out), allocatable :: fl(:,:,:)
!! The transpose of al is the Inverse Fourier Transform applied to the Fourier transform of the potential.

!! AUXILIAR VARIABLES
  integer :: nk, nr !! number of nodes on each FE-DVR, momentum and position.
  integer :: i, j,k,l !! indexes to run the loop

!! needed for initial_fedvr3d_radial_momentum
  integer, allocatable :: element(:), basis(:)
  real(kind=k1), allocatable :: fedvr_k_node(:,:), fedvr_k_weight(:,:)
  real(kind=k1), allocatable :: fedvr_k_total(:), fedvr_k_weight_total(:)

!! SPHERICAL BESSEL FUNCTIONS

INTERFACE
   Function besselj(n,x)
     integer :: n
     real*8 :: x
     real*8 :: besselj
   End Function besselj
END INTERFACE
  
nr=size(fedvr_r)  !! size of the position FE-DVR

!! Construct the FE-DVR in the momentum space to build the matrix al and fl.

call initial_fedvr3d_radial_momentum(k_max, 0.0d0,which_element,which_basis,fedvr_k_weight,fedvr_k_node, fedvr_k_total)

!! Allocate the FE-DVR in momentum space with the same size than in spatial space. This can be modified.

nk=size(fedvr_k_total)  !! size of the momentum FE-DVR.
allocate(fedvr_k_weight_total(1:nk)) !! allocate the global weights.


!! Rearrange the weights

Do j=1,size(which_element)-1
   If (which_basis(j).gt.which_basis(j+1)) then !! if the functions is a border of the element.
      fedvr_k_weight_total(j)=fedvr_k_weight(which_element(j),which_basis(j))+fedvr_k_weight(which_element(j)+1,1)
   Else  !! if it is not a border of the element
      fedvr_k_weight_total(j)=fedvr_k_weight(which_element(j),which_basis(j))
   End If
End Do

fedvr_k_weight_total(nk)= fedvr_k_weight(which_element(nk),which_basis(nk))!! Last weight function


!! Deallocate the variables that we do not need.

deallocate(fedvr_k_weight,fedvr_k_node)

!! We construct the matrices for each angular momentum

allocate(al(1:nk,1:nr,0:2*l_max)) !! allocate the memory for the matrices.
allocate(fl(1:nr,1:nr,0:2*l_max)) !! allocate the memory for the matrices.

al=zero !! Initialize the matrices.
fl=zero

Do l=0,2*l_max !! Run for the angular momentum
   
   Do j=1,nr !! Run for the position grid
      
      Do i=1,nk !! Run for the momentum grid
         
         al(i,j,l)=besselj(l,fedvr_k_total(i)*fedvr_r(j))         

      End Do

   End Do

!! calculation of the the matrix fl.

   Do j=1,nr !! Running in the indexes of fl. These indexes correspond to the
      !! position space.
      Do i=1,nr 
         Do k=1,nk !! Runing in the internal indexes fo al. This index corresponds to the momentum
            
            fl(i,j,l)=fl(i,j,l)+al(k,i,l)*al(k,j,l)*fedvr_k_weight_total(k)
            
         End Do
      End Do
   End Do
   
End Do

deallocate(fedvr_k_weight_total,fedvr_k_total)

fl=4.0d0*pi*fl  !! we include the factor 4\pi, which comes from the transform of the coulomb potential

return

End Subroutine al_make

!!=====================================================================
!!
!! phikphil_to_phikl transforms two orbitals from uncoupled basis to
!! the coupled basis. Note that the operation is 
!!
!! phi_k* phi_l -> phi_{kl}
!!
!!=====================================================================

Subroutine phikphil_to_phikl(phi, c_llmm_lm, phi_coupled)
implicit none

!! INPUT

complex (kind=k2), allocatable, intent(in):: phi(:,:) !! orbitals in the uncoupled basis

type (llmm_lm_type), intent(in) :: c_llmm_lm !! coefficients which link the uncouple basis with the coupled basis for Spherical Harmonics. 

!! OUTPUT

complex (kind=k2), allocatable, intent(out) :: phi_coupled(:,:,:) !! orbitals in the coupled basis

!! AUXILIAR VARIABLES

integer :: i, j, k, l
integer :: alpha, beta, gamma !! number of angular functions.
real(kind=k1) :: aux_coupling
integer :: alpha_max, beta_max, gamma_max !! total number of angular functions.
integer :: nr !! number of radial points.
integer :: num_couplings !! length of the c_llmm_lm matrix


interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

!! Check which is the largest value of the gamma.

alpha_max=maxval(c_llmm_lm%alpha(:))
beta_max=maxval(c_llmm_lm%beta(:))
gamma_max=maxval(c_llmm_lm%gamma(:))

nr=size(phi(:,1))/alpha_max !! since phi(1:nr*alpha_max)

!! allocate the coupling functions.

If (allocated(phi_coupled)) then
   continue
else
   allocate(phi_coupled(1:nr*gamma_max,size(phi(1,:)),size(phi(1,:))))
End If

!! initialize the coupling functions.
   phi_coupled=zzero

!! number of couplings among angular momentum

num_couplings=size(c_llmm_lm%alpha(:)) 

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(l,i,alpha,beta,gamma,aux_coupling,j) 
!$OMP DO SCHEDULE(dynamic)!, 2200)

Do k=1,size(phi(1,:)) !! Run in the orbitals for the complex conjugate
   Do l=1,size(phi(1,:)) !! Run in the orbitals
      Do i=1, num_couplings !! Run in the non-vanishing couplings of the angular functions.
         alpha=c_llmm_lm%alpha(i) !! Y_alpha*
         beta=c_llmm_lm%beta(i) !! Y_beta
         gamma=c_llmm_lm%gamma(i) !! Y_gamma
         aux_coupling=c_llmm_lm%value(i) !! Coupling(alpha, beta, gamma)
         Do j=1,nr !! Run in the FE-DVR
            phi_coupled(nr*(gamma-1)+j,k,l)= phi_coupled(nr*(gamma-1)+j,k,l)+aux_coupling*conjg(phi(nr*(alpha-1)+j,k))*phi(nr*(beta-1)+j,l)
         End Do
      End Do
   End Do
End Do

!$OMP END DO
!$OMP END PARALLEL

return

End Subroutine phikphil_to_phikl

!!=============================================================
!!
!! Calculation of the Mean Field Operator using the angular
!! coupled basis.
!!
!! The output is an array of vectors written in the coupled
!! basis. Each of them corresponds to the vector phi_k* phi_l.
!!
!! Subroutine: mean_field_coupled
!!
!!=============================================================

Subroutine mean_field_coupled(fl,phi_coupled,lmax,mmax,lmax_coupled, mean_field)
implicit none

!! INPUT

real(kind=k1), allocatable, intent(in) :: fl(:,:,:) 
!! fl performs the operation to obtain the mean field operator using
!! the coupled basis
!!
!! Arguments
!! #1 is the row of the L block.
!! #2 is the column of the L block.
!! #3 is the L block matrix.

complex(kind=k2), allocatable, intent(in) :: phi_coupled(:,:,:)
!! phi_coupled(:,k,l) is the product of two orbital functions 
!! phi(:,k)* phi(:,l) in the coupled basis.
!!
!! Arguments
!! #1 is the coefficient in the coupled basis.
!! #2 and #3 just explained above.

integer, intent(in) :: lmax,mmax ,lmax_coupled
!! lmax is the maximum value of the orbital angular momentum.
!! mmax is the maximum value of the orbital magnetic quantum number.
!! lmax_coupled is the maximum component in the coupled representation

!! OUTPUT

complex(kind=k2), allocatable, intent(inout) :: mean_field(:,:,:)
!! Mean field is the the mean field operator.
!! Arguments
!!
!! #1 is the coefficient in the coupled basis.
!! #2 orbital complex conjugated.
!! #3 orbital without complex conjugate.

!! AUXILIAR

integer :: l_i,m_i,k, counter
integer :: l !! auxiliar variable to store temporarily the 
             !! total angular momentum.
integer :: nr, nangle, nrangle, norbital
!! nr is the number of radial poits
!! nangle is the number of angular functions (in the coupled basis)
!! nrangle is the number of functions of the radial and angular functions
!!         in the coupled basis
!! norbital is the number of orbitals
integer :: orb1, orb2 !! indeces in the loops for the orbitals.
integer :: r1, r2 !! indeces for the loops in the radial part.
integer :: i_angle !! index for the loop in the angular functions.
integer, allocatable :: lcoupled(:), mcoupled(:) !! l and m in the coupled basis


interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

!! First we allocate the mean_field operator array.

nrangle=size(phi_coupled(:,1,1))
nr=fedvr3d%nb_r
norbital=size(phi_coupled(1,1,:))

!! Check that it is not allocated

If (.not.allocated(mean_field)) allocate(mean_field(1:nrangle,1:norbital,1:norbital))

!! Initialize it to zero

mean_field=zzero

!! Calculate the lcoupled and mcoupled indexes
!! 1. Number of angular functions

counter=0
Do l_i=0,lmax_coupled
   Do m_i=-min(2*mmax,l_i),min(2*mmax,l_i)
      counter=counter+1
   End Do
End Do

If (counter.ne.nrangle/nr) then
   print*, counter, nrangle/nr
   write(*,*) 'ERROR IN SUBROUTINE mean_field_coupled'
   write(*,*)
   write(*,*) 'The number of coefficients in the coupled basis'
   write(*,*) 'does not correspond to l_max_coupled'
   stop
End If

!! 2. storing the indexes of the angular momentum and magnetic quantum number.

nangle=counter

allocate(lcoupled(nangle)) 
allocate(mcoupled(nangle))

!! 3. initialize the indexes

lcoupled=0
mcoupled=0

!! 4. angular momentum and magnetic quantum numbers for the coupled basis.

counter=0
Do l_i=0,lmax_coupled
   Do m_i=-min(2*mmax,l_i),min(2*mmax,l_i)
      counter=counter+1
      lcoupled(counter)=l_i
      mcoupled(counter)=m_i
   End Do
End Do

!! ^                                            ^
!! |                                            |
!!
!! This allocation of the lcoupled is always the same.
!! It can be included as a input in the subroutine.

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(orb1,i_angle,l,r2,r1)
!$OMP DO SCHEDULE(dynamic)!, 2200)

Do orb2=1,norbital !! Run in the orbital without complex conjugate.
   Do orb1=1,norbital !! Run in the orbital complex conjugate.
      Do i_angle=1,nangle !! Run in the angular functions of the coupled basis.
         l=lcoupled(i_angle)
         Do r2=1,nr !! Index of the radial functions we use to calculate
                    !! Mean field operator. They are under the integral.
            Do r1=1,nr !! Radial dependence of the Mean Field Operator
               
               mean_field((i_angle-1)*nr+r1,orb1,orb2)=mean_field((i_angle-1)*nr+r1,orb1,orb2)+fl(r1,r2,l)*phi_coupled((i_angle-1)*nr+r2,orb1,orb2)

            End Do
         End Do
      End Do
   End Do
End Do

!$OMP END DO
!$OMP END PARALLEL

deallocate(lcoupled,mcoupled)

return

End Subroutine mean_field_coupled

!!=======================================================================
!! SUBROUTINES TO SOLVE THE Q-SPACE EQUATIONS USING THE COUPLED 
!! REPRESENTATION.
!!=======================================================================

!!-----------------------------------------------------------------------
!! mean-field_phikl calculates the mean field operator in the 
!! coupled representation and the products \phi_p\dagger \phi_q
!!-----------------------------------------------------------------------

Subroutine mean_field_and_phikl(rl, orb, c_llmm_lm, mf, orb_coupled)
use global
Implicit none

!! INPUT

real(kind=k1), allocatable, intent(in) :: rl(:,:,:) 
!! Array to obtain the Mean-Field-Operator (W_{pq}(r_1)), R_L(r1,r2,L)
!! Arguments:
!! 
!! #1: Radial coordinate of the mean field operator.
!! #2: Radial coordinate of the integration.
!! #3: Angular momentum of the coupled representation. It goes from 0 to 2*l_max

complex(kind=k2), allocatable, intent(in) :: orb(:,:)
!! orb stores the orbitals
!! Arguments
!! 
!! #1: component of the global basis (FE-DVR+Angular).
!! #2: number of the orbital.

type (llmm_lm_type), intent(in) :: c_llmm_lm !! coefficients which link the uncouple basis with the coupled basis for Spherical Harmonics. 

!! OUTPUT
!!
!!
complex (kind=k2), allocatable, intent(out) :: orb_coupled(:,:,:)
!! orb_coupled(:,k,l) is the product of two orbital functions 
!! orb(:,k)* phi(:,l) in the coupled basis.
!!
!! Arguments
!! #1 is the coefficient in the coupled basis.
!! #2 Orbital which is complex conjugated (coming from the bra).
!! #3 Orbital which is not complex conjujated (coming from the ket).

!!
complex(kind=k2), allocatable, intent(out) :: mf(:,:,:)
!! mf is the Mean Field Operator written in the coupled basis.
!! Arguments
!!
!! #1 Coupled basis function
!! #2 Orbital which is complex conjugated (coming from the bra).
!! #3 Orbital which is not complex conjujated (coming from the ket).
!!

!! AUXILIAR VARIABLES

!! First, we calculate the product of the two orbitals.

call phikphil_to_phikl(orb,c_llmm_lm, orb_coupled)

!! Second, we calculate the mean field operator

call mean_field_coupled(rl,orb_coupled,fedvr3d%l_max,fedvr3d%m_max,fedvr3d%l_max_coupled, mf)  !! calculate the mean field operator in the coupled representation

return
End Subroutine mean_field_and_phikl

!!--------------------------------------------------------------------------------------------------
!! mean_field_phi calculates the application of the density operators and the mean field operator 
!! on the orbitals. 
!!--------------------------------------------------------------------------------------------------

Subroutine mean_field_phi(mf,orb,int_llmm_lm, c3,korb)  !! 24.02.2015 CHECK ***!!!
Implicit none

!! INPUT
!!
complex(kind=k2), allocatable, intent(in) :: mf(:,:,:)
!! mf is the Mean Field Operator written in the coupled basis.
!! Arguments
!!
!! #1 Coupled basis function
!! #2 Orbital which is complex conjugated (coming from the bra).
!! #3 Orbital which is not complex conjujated (coming from the ket).
!!
complex(kind=k2), allocatable,intent(in) :: orb(:,:)
!! orb stores the orbitals
!! Arguments
!! 
!! #1: component of the global basis (FE-DVR+Angular).
!! #2: number of the orbital.
!!
type (llmm_lm_type), intent(in) :: int_llmm_lm 
!! Integrals of Y_{alpha}* Y_{beta} Y_{gamma}, where alpha and beta are in the uncoupled basis and gamma is in the coupled basis.
!!
!! %alpha : angular function of Y_{alpha}
!! %beta : angular function of Y_{beta}
!! %gamma : angular function of Y_{gamma}
!! %value : angular function of Y_{alpha} 
!!
complex(kind=k2), allocatable, intent(in) :: c3(:,:,:,:)
!!
!! c3 is the result of the inverse one body density matrix times the two body density matrix
!! The arguments are the orbitals.
!!
!! OUTPUT
!!
complex(kind=k2), intent(out) :: korb(:,:)
!! orb stores the orbitals
!! Arguments
!! 
!! #1: component of the global basis (FE-DVR+Angular).
!! #2: number of the orbital.
!!

!! AUXILIAR VARIABLES
integer :: i,j,k,nglobal, norb,ix,iy_here,iz
integer :: ij,il,ik, i_spatial !! are the indexes for the loops in the orbitals
integer :: nangle, nr !! total number of angular

integer :: orb_index(1:system%nptot**4,1:4) !! running in the useful orbitals
!! Arguments: 
!!
!! #1 is the counter. It runs over the possible number of orbitals.
!! #2 labels the orbital. It goes from 1 to 4. 1-> ij, 2-> il, 3-> ik, 4-> i_spatial
integer:: counter,total
real*8:: finish, start !! To measure CPU time. This is used to check the time for debugging.

complex(kind=k2), allocatable :: korb_omp(:,:,:) !! OMP auxiliar variable. It stores the information for the different threads.
!! Arguments:
!!
!! #1: component of the global basis (FE-DVR+Angular).
!! #2: number of the orbital.
!! #3: thread.

integer :: proc !! number of processors
integer :: number_proc !! threads

interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

!! initialize the output

nglobal=size(orb(:,1)) !! global basis
norb=size(orb(1,:)) !! number of orbitals
nr=fedvr3d%nb_r !! number of radial points
nangle=nglobal/nr !! number of angular functions

!! Initialize the index for the combination of orbitals which take part
!! in the calculation orbitals.
 
orb_index=0
total=0
Do ij=1,norb
   do il=1,norb
      do ik=1,norb
         do i_spatial =1,norb
          If (abs(c3(i_spatial,ik,il,ij)).gt.1d-15) then !! contribution to the projection
             total=total+1
             orb_index(total,1)=ij
             orb_index(total,2)=ik
             orb_index(total,3)=il
             orb_index(total,4)=i_spatial
          End If
         End do
      End do
   End do
End Do

!! Now we store the 

!! Prepare the output !! JUAN: taken out because it affects to a variable defined elsewhere.
!! in the main code it is not allocatable

!!$If (allocated(korb)) then 
!!$   continue
!!$Else
!!$   allocate(korb(1:nglobal,1:norb)) !! allocate memory
   korb=zzero !! initialize to zero
!!$End If

!!$call cpu_time(start)

!! Set the number of processors

proc=OMP_get_num_procs()
!!      call OMP_set_num_threads(8)
!! initialize the auxiliar variable for parallelization

allocate(korb_omp(1:nglobal,1:norb,0:proc-1))

korb_omp=zzero

!! Loop to solve the equations


!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(ij,ik,il,i_spatial,j,i,ix,iy_here,iz) 
!$OMP DO SCHEDULE(dynamic)

Do counter=1, total
   ij=orb_index(counter,1)
   ik=orb_index(counter,2) 
   il=orb_index(counter,3)
   i_spatial=orb_index(counter,4)
   Do j=1,size(int_llmm_lm%alpha(:)) !! Run in the angular couplings.
      Do i=1,nr !! Run in the radial grid
         ix=(int_llmm_lm%alpha(j)-1)*nr+i !! global function of the result
         iy_here=(int_llmm_lm%beta(j)-1)*nr+i !! global function of the ket
         iz=(int_llmm_lm%gamma(j)-1)*nr+i !! global function of the coupled representation
!         korb(ix,i_spatial)=korb(ix,i_spatial)+c3(i_spatial,ik,il,ij)*mf(iz,ik,il)*phi(iy_here,ij)*int_llmm_lm%value(j)          
         korb_omp(ix,i_spatial,OMP_get_thread_num())=korb_omp(ix,i_spatial,OMP_get_thread_num())+c3(i_spatial,ik,il,ij)*mf(iz,ik,il)*phi(iy_here,ij)*int_llmm_lm%value(j)          
      End Do !! loop in i
   End Do !! loop in j
End Do !! loop in counter

!$OMP END DO
!$OMP END PARALLEL


!!$call cpu_time(finish)
!!$print*, finish-start
!!$print*, 'module_operator_spatial_space_omp.f90'

!! sum all the contributions

Do number_proc=0,proc-1
   korb(:,:)=korb(:,:)+korb_omp(:,:,number_proc)
End Do

!! deallocate auxiliar variables for parallelization.

deallocate(korb_omp)

return
End Subroutine mean_field_phi 

!!--------------------------------------------------------------------------------------------------
!! two_body_coupled calculates the two body matrix elements
!!--------------------------------------------------------------------------------------------------

Subroutine two_body_coupled(mf,orb_coupled,tei_coupled) !! 24.02.2015 ***CHECK***!!
implicit none

!! INPUT
!!
complex (kind=k2), allocatable, intent(in) :: mf(:,:,:) 
!! mf is the Mean Field Operator written in the coupled basis.
!! Arguments
!!
!! #1 Coupled basis function
!! #2 Orbital which is complex conjugated (coming from the bra).
!! #3 Orbital which is not complex conjujated (coming from the ket).
!!
complex (kind=k2), allocatable, intent(in) :: orb_coupled(:,:,:)
!! orb_coupled(:,k,l) is the product of two orbital functions 
!! orb(:,k)* phi(:,l) in the coupled basis.
!!
!! Arguments
!! #1 is the coefficient in the coupled basis.
!! #2 Orbital which is complex conjugated (coming from the bra).
!! #3 Orbital which is not complex conjujated (coming from the ket).

!! OUTPUT
!!
complex (kind=k2), allocatable, intent(out) :: tei_coupled(:,:,:,:)
!!
!! two body matrix element
!!
!! Arguments
!! #1-#4 are the orbitals.

!! AUXILIAR VARIABLES
!!
integer :: ij, il, ik, i_spatial !! loops in the orbitals
integer :: ix !! global basis in the coupled basis.
integer :: norb !! number of orbitals
integer :: nglobal !! number of basis functions


interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

norb=size(orb_coupled(1,1,:)) !! number of orbitals
nglobal=size(orb_coupled(:,1,1)) !! number of basis functions

!! Prepare the output

if (allocated(tei_coupled)) then
   continue
Else
   allocate(tei_coupled(1:norb,1:norb,1:norb,1:norb)) !! allocate memory
   tei_coupled=zzero !! initialize to cmplx(0.0d0, 0.0d0)
End if

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(ix,il,ik,i_spatial)
!$OMP DO SCHEDULE(dynamic)

Do ij=1,norb
   Do ix=1,nglobal
      Do il=1,norb
         Do ik=1,norb
            Do i_spatial=1,norb
!!$               tei_coupled(i_spatial,ik,il,ij)=tei_coupled(i_spatial,ik,il,ij)+orb_coupled(ix,i_spatial,ij)*mf(ix,ik,il)
               tei_coupled(i_spatial,ik,il,ij)=tei_coupled(i_spatial,ik,il,ij)+conjg(orb_coupled(ix,ij,i_spatial))*mf(ix,ik,il) !!27th January
            End Do !! Loop in i_spatial
         End Do !! Loop in ik
      End Do !! Loop in il
   End Do !! Loop in ix
End Do !! Loop in ij

!$OMP END DO
!$OMP END PARALLEL

return
End Subroutine two_body_coupled

!!===================================================================
!! Subroutine to apply the Hamiltonian on a orbital
!!===================================================================

!! subroutine act_stif_2 performs the result of applying the 
!! one body operator on the orbitals.

Subroutine act_stif_2(orb,hl, laserp,korb)

!! INPUT 
complex(kind=k2), allocatable, intent(in) :: orb(:,:)
!! orb stores the orbitals
!! Arguments
!! 
!! #1: component of the global basis (FE-DVR+Angular).
!! #2: number of the orbital.

real(kind=k1), allocatable, intent(in) :: hl(:,:,:)
  !! Hamiltonian after the stiffness procedure.
  !! The hamiltonian matrix is the same for each value of the 
  !! angular momentum l. 
  !! Arguments
  !!
  !! #1: FE-DVR function of the bra.
  !! #2: FE-DVR function of the ket.
  !! #3: Angular momentum of the Hamiltonian, l.
  !!

type(laser_prime), intent(in) :: laserp
!! laserp stores the properties of the laser

!! OUTPUT

complex(kind=k2),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot), intent(out) :: korb
!! orb stores the orbitals
!! Arguments
!! 
!! #1: component of the global basis (FE-DVR+Angular).
!! #2: number of the orbital.

!! AUXILIAR VARIABLES
complex(kind=k2), allocatable :: korb_aux(:,:)
integer:: iangle,ibra,iket,iorb !! indexes for the loops.
integer :: iangle_bra, iangle_ket !! indexes for the loops in angular functions
integer:: norbital, nbasis !! number of orbitals and basis fuctions.
integer:: nbasis_r !! number of radial fucntions
integer:: l_basis !!angular momentum
real(kind=k1) :: ang_x, ang_y, ang_z !! \cos\theta,\sin\theta\cos\phi, and \sin\theta\sin\phi contributions
!! for the differents functions.

interface
   function OMP_get_thread_num()
     integer :: OMP_get_thread_num
   end function OMP_get_thread_num
   function OMP_get_num_procs()
     integer :: OMP_get_num_procs
   end function OMP_get_num_procs
   function OMP_get_num_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_max_threads()
     integer :: OMP_get_num_threads
   end function OMP_get_max_threads
   subroutine omp_set_num_threads(num_threads)
     integer, intent(in) :: num_threads
   end subroutine omp_set_num_threads
!!$        function  omp_set_num_threads()
!!$          integer :: OMP_set_num_threads
!!$        End function omp_set_num_threads
   function OMP_get_nested()
     logical :: OMP_get_nested
   end function OMP_get_nested
   subroutine OMP_set_nested(enable)
     logical, intent(in) :: enable
   end subroutine OMP_set_nested
end interface

norbital=size(orb(1,:)) !! number of orbitals
nbasis=size(orb(:,1)) !! number of basis functions
nbasis_r=size(hl(:,1,1)) !! number of radial functions

!! Initialize the output, korb.

!!$If (allocated(korb)) then !! DEBUG 11/5/2015
!!$   continue
!!$else
!!$   allocate(korb(1:nbasis,1:norbital))
!!$End If

korb=zzero

!! First we evaluate in the field free case.
!! To run in the angular functions we use lm_l(#1) which returns
!! the angular momentum of the basis function #1.

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(iangle,l_basis,ibra,iket)
!$OMP DO SCHEDULE(dynamic)

Do iorb=1,norbital !! Run in the orbitals
   Do iangle=1, size(lm_l(:)) !! run in the angular functions
      l_basis=lm_l(iangle) !! angular momentum
      Do ibra=1,nbasis_r !! run in the radial basis
         Do iket=1,nbasis_r !! run in the radial basis
            korb(nbasis_r*(iangle-1)+ibra,iorb)=korb(nbasis_r*(iangle-1)+ibra,iorb)+orb(nbasis_r*(iangle-1)+iket,iorb)*hl(ibra,iket,l_basis)
         End Do !! iket
      End Do !! ibra
   End Do !! iangle
End Do !! iorb

!$OMP END DO
!$OMP END PARALLEL

If (laserp%tdornot) then

   Select case(laserp%gauge)
   case('l')
   
!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(iangle_ket,iangle_bra,ang_x,ang_y,ang_z,iket)
!$OMP DO SCHEDULE(dynamic)
   
   Do iorb=1,norbital !! run in the orbitals
      Do iangle_ket=1,size(lm_l(:))!! run in the angular functions of the ket
         Do iangle_bra=1,size(lm_l(:)) !! run in the angular functions of the bra
            ang_x=xmat_3d(iangle_bra,iangle_ket) !! angular contribution of x direction
            ang_y=ymat_3d(iangle_bra,iangle_ket) !! angular contribution of y direction
            ang_z=zmat_3d(iangle_bra,iangle_ket) !! angular contribution of z direction
            Do iket=1,nbasis_r !! run in the radial basis               
               korb(nbasis_r*(iangle_bra-1)+iket,iorb)=korb(nbasis_r*(iangle_bra-1)+iket,iorb)+orb(nbasis_r*(iangle_ket-1)+iket,iorb)*fedvrx_global(iket)*(laserp%ex_t*ang_x+laserp%ey_t*ang_y*ci+laserp%ez_t*ang_z)               
               !! multiply by the imaginary unity the contribution in the Y direction.
               !! Add the absorbing potential
               
               If (iangle_ket.eq.iangle_bra) then
                  korb(nbasis_r*(iangle_ket-1)+iket,iorb)=korb(nbasis_r*(iangle_ket-1)+iket,iorb)+orb(nbasis_r*(iangle_ket-1)+iket,iorb)*vabsorb_pot_fedvr3d(iket)*ci !! we add the imaginary unit here, which is part of the absorbing potential.
               End If

            End Do !! iket
            !! initialize the variables
            ang_x=zero
            ang_y=zero
            ang_z=zero
         End Do !! iangle_bra
      End Do !! iangle_ket
   End Do !! iorb

!$OMP END DO
!$OMP END PARALLEL

case('v')

!! Absorbing potential

!!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(iangle_ket,iket)
!$OMP DO SCHEDULE(dynamic)
   
   Do iorb=1,norbital !! run in the orbitals
      Do iangle_ket=1,size(lm_l(:))!! run in the angular functions of the ket
         Do iket=1,nbasis_r !! run in the radial basis               
            
              korb(nbasis_r*(iangle_ket-1)+iket,iorb)=korb(nbasis_r*(iangle_ket-1)+iket,iorb)+orb(nbasis_r*(iangle_ket-1)+iket,iorb)*vabsorb_pot_fedvr3d(iket)*ci !! we add the imaginary unit here, which is part of the absorbing potential.
            
         End Do !! iangle_ket
      End Do !! iket
   End Do !! iorb
   
!$OMP END DO
!$OMP END PARALLEL

   !! DEBUG: JUAN
   !! Velocity gauge   

   call act_index(orb,pz_element,index_pz_element,korb_aux)  !! Apply pz
   korb=korb+korb_aux*laserp%az_t !! times the vector potential

   deallocate(korb_aux)

!!$   call act_index(orb,py_element,index_py_element,korb_aux)  !! Apply py
!!$
!!$   korb=korb+korb_aux*laserp%ay_t !! times the vector potential
!!$
!!$   deallocate(korb_aux)
!!$
!!$   call act_index(orb,px_element,index_px_element,korb_aux) !! Apply px
!!$
!!$   korb=korb+korb_aux*laserp%ax_t !! times the vector potential
!!$   
!!$   deallocate(korb_aux)
   
case default

   write(*,*) ''
   write(*,*) 'ERROR'
   write(*,*) ''
   write(*,*) 'We should use length (l) or velocity (v) gauge'
   write(*,*) ''
   
   stop

End Select

End If

return
End Subroutine act_stif_2

end module operator_spatial_space

!
! configure the address
!
 subroutine configure_operator_spatial_space
  use operator_spatial_space
  implicit none

   allocate(ch_dummy(system%nptot,system%nptot))

   allocate(ham_onebody(system%nptot*(system%nptot+1)/2))
   allocate(xmat_spatial_space(system%nptot*(system%nptot+1)/2))
   allocate(dvdxmat_spatial_space(system%nptot*(system%nptot+1)/2))
   allocate(ymat_spatial_space(system%nptot*(system%nptot+1)/2))
   allocate(dvdymat_spatial_space(system%nptot*(system%nptot+1)/2))
   allocate(zmat_spatial_space(system%nptot*(system%nptot+1)/2))
   allocate(dvdzmat_spatial_space(system%nptot*(system%nptot+1)/2))
  
   allocate(tei_spatial(system%nptot,system%nptot,system%nptot,system%nptot))
   allocate(cww(fedvr3d%nb_angle*fedvr3d%nb_r,fedvr3d%nb_angle,system%nptot,&
   system%nptot))

  return
 end subroutine configure_operator_spatial_space

!
! calc. the onebody operator in the spatial space
!
   subroutine update_ham_onebody()
    use operator_spatial_space
    implicit none

!=================================================================================
! the stiffness are not considered
!================================================================================= 
    select case (prop%stiffness)
       case(0)
       call update_ham_onebody_spatial_space()
!=================================================================================
! the stiffness are considered
!=================================================================================
    case(1) !! Wenliang's case
       call update_ham_onebody_spatial_space2()         
    case(2) !! Juan's case        
        call update_ham_onebody_spatial_space_stif_2(phi,hl_stiffness,laser,ham_onebody)
     case default !! Error in the input.txt
        write(*,*) 'ERROR in update_ham_onebody in'
        write(*,*) 'module_operator_spatial_space_omp.f90'
        stop       
     End select

    return
   end subroutine update_ham_onebody

!
! calc. the twobody operator in the spatial space
!
   subroutine update_vv_twobody()
    use operator_spatial_space
    implicit none
     call update_twobody_spatial_space()
     
    return
   end subroutine update_vv_twobody

!
! calc. the x dipole operator in the spatial space
!
   subroutine update_xmat_spatial_space()
    use global
    use operator_3d
    use operator_spatial_space
    implicit none
     call update_xyzmat_spatial_space(xmat_3d,fedvr3d%nb_angle,xmat_spatial_space,system%nptot)
    return
   end subroutine update_xmat_spatial_space



!
! calc. the dvdx  operator in the spatial space
!
   subroutine update_dvdxmat_spatial_space()
    use global
    use operator_3d
    use operator_spatial_space
    implicit none
     call update_dvdxyzmat_spatial_space(xmat_3d,fedvr3d%nb_angle,dvdxmat_spatial_space,system%nptot)
    return
   end subroutine update_dvdxmat_spatial_space

 
!================================ Y ==========================================================
!
! calc. the y dipole operator in the spatial space
!
   subroutine update_ymat_spatial_space()
    use global
    use operator_3d
    use operator_spatial_space
    implicit none
     call update_xyzmat_spatial_space(ymat_3d,fedvr3d%nb_angle,ymat_spatial_space,system%nptot)
     return
   end subroutine update_ymat_spatial_space

!
! calc. the dvdy  operator in the spatial space
!
   subroutine update_dvdymat_spatial_space()
    use global
    use operator_3d
    use operator_spatial_space
    implicit none
     call update_dvdxyzmat_spatial_space(ymat_3d,fedvr3d%nb_angle,dvdymat_spatial_space,system%nptot)
    return
   end subroutine update_dvdymat_spatial_space


!====================================== Z ===================================================

!
! calc. the z dipole operator in the spatial space
!
   subroutine update_zmat_spatial_space()
    use global
    use operator_3d
    use operator_spatial_space
    implicit none
     call update_xyzmat_spatial_space(zmat_3d,fedvr3d%nb_angle,zmat_spatial_space,system%nptot)
     return
   end subroutine update_zmat_spatial_space


!
! calc. the dvdz  operator in the spatial space
!
   subroutine update_dvdzmat_spatial_space()
    use global
    use operator_3d
    use operator_spatial_space 
    implicit none

     call update_dvdxyzmat_spatial_space(zmat_3d,fedvr3d%nb_angle,dvdzmat_spatial_space,system%nptot)
    return
   end subroutine update_dvdzmat_spatial_space







