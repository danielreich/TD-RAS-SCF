!===============================================================================
! calc. the matrix element of kinetic energy and nuclei attractive
! potential, dipole matrix in fedvr_3d, spherial coordiante 
!===============================================================================
 module operator_spatial_space
  use global
  use operator_3d
  use wfunction
  use twoe_basis_set

  implicit none 
  complex(kind=k2),allocatable,save :: ham_onebody(:),xmat_spatial_space(:),&
  ymat_spatial_space(:),zmat_spatial_space(:),dvdxmat_spatial_space(:),&
  dvdymat_spatial_space(:),dvdzmat_spatial_space(:)
  complex(kind=k2),allocatable,save :: ch_dummy(:,:)
  complex(kind=k2),allocatable,save :: tei_spatial(:,:,:,:),cww(:,:,:,:)

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

!!$              print*, 'vabsorb_pot_fedvr3d is not well used here'
!!$              print*, 'module_operator_spatial_space.f90'


             ! sumtemp_laser_term = sumtemp_laser_term + dconjg(phi(idicp_prime,i_spatial ))* phi(jdicp_prime,j_spatial)*&
             ! xmat_3d(idicp,jdicp)*fedvrx_global(kdicp)*laser%ex_t


             ! sumtemp_laser_term = sumtemp_laser_term + dconjg(phi(idicp_prime,i_spatial ))* phi(jdicp_prime,j_spatial)*&
             ! ymat_3d(idicp,jdicp)*fedvrx_global(kdicp)*ci*laser%ey_t

           enddo
        enddo
     enddo
     endif 

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
              (zmat_3d(idicp,jdicp)*fedvrx_global(kdicp)*laser%ez_t + vabsorb_pot_fedvr3d(kdicp))

              print*, 'vabsorb_pot_fedvr3d is not well used here -stif-'
              print*, 'module_operator_spatial_space.f90'
              stop

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
     
     select case (fedvr3d%store)

     case (0)  !! do not store cww

        call  two_body_orbitals(tei_spatial)
                
     case (1)  !! stores cww

     cww = zzero

     do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
       do jdicp =1, fedvr3d%nb_angle
        do itwo = 1, num_two(idicp,jdicp)

            kdicp = index_two_storage(1,itwo,idicp,jdicp)
            ldicp = index_two_storage(2,itwo,idicp,jdicp)

            do mdicp = 1, system%nptot
              do ndicp =1, system%nptot
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

    do imo =1, system%nptot
      do jmo =1, system%nptot
       do kmo =1, system%nptot
        do lmo = 1,system%nptot
        
          sumtemp = zzero
          do i_fedvr3d =1,fedvr3d%nb_r*fedvr3d%nb_angle
            do j_fedvr3d =1,fedvr3d%nb_angle

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

         if(xmat_3d(i_angle,j_angle)/=0.0d0) then 

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
integer :: kdicp, ldicp, ldicp_ang,kdicp_ang
integer :: n1, l1,m1,n2,l2,m2,n3,l3,m3,l4,m4
double precision:: twoe_aux

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

If (n1.ne.n2) return !! this contribution is zero

cww_calc=cmplx(0.0d0,0.0d0)

!$ omp do

Do kdicp=1,fedvr3d%nb_angle*fedvr3d%nb_r !! Run in the FEDVR basis
   If (abs(phi(kdicp,phi1)).lt.1d-15) cycle
   Do  ldicp_ang=1,fedvr3d%nb_angle !! Run in the angular FEDVR basis for the angular part (since the radial is the same).

      n3=global_to_local(kdicp,1)
      l3=global_to_local(kdicp,2)
      m3=global_to_local(kdicp,3)

      ldicp=(ldicp_ang-1)*fedvr3d%nb_r+n3

      l4=global_to_local(ldicp,2)
      m4=global_to_local(ldicp,3)

      If (abs(dconjg(phi(kdicp,phi1))*phi(ldicp,phi2)).lt.1d-20) then
!         write(123,*) kdicp,ldicp, abs(dconjg(phi(kdicp,phi1))*phi(ldicp,phi2))
!         stop
         cycle !! change 12 november AQUI
      end If

      if ((m1-m2).ne.(m4-m3)) cycle !! This contribution is zero

      cww_calc = cww_calc + dconjg(phi(kdicp,phi1))*phi(ldicp,phi2)*fedvr3d_twoe(n1,l1,m1,l2,m2,n3,l3,m3,l4,m4)! 

   Enddo
Enddo

!$ omp end do

return
End Function cww_calc

!! Calculation of the two body integral in the orbitals. This a subroutine which stores in tei_spatial(imo,kmo,lmo,jmo)

!!tei_spatial(imo,kmo,lmo,jmo) is the integral phi_imo*(r1)phi_kmo*(r2)(1/r12)phi_jmo(r1)phi_lmo(r2)

Subroutine two_body_orbitals(tei_spatial)
implicit none
complex(kind=k2),allocatable, intent(inout) :: tei_spatial(:,:,:,:) !! orbitals
!!Auxiliary variable
 complex(kind=k2) :: sumtemp
integer :: imo,jmo,kmo,lmo,i_fedvr3d,j_fedvr3d
integer :: i_global, j_global

If (allocated(tei_spatial)) then
   tei_spatial=zzero
   else
 !! If tei_spatial is not allocated, do it!
   allocate(tei_spatial(1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))
   tei_spatial=zzero
End If
!!$
!!$do imo =1, system%nptot !! Run in orbitals
!!$   do kmo =imo, system%nptot !! Run in orbitals
!!$      do jmo =1, system%nptot !! Run in orbitals
!!$         do lmo = 1,system%nptot !! Run in orbitals
!!$            
!!$            sumtemp = zzero
!!$            do i_fedvr3d =1,fedvr3d%nb_r*fedvr3d%nb_angle
!!$
!!$                  i_global = i_fedvr3d
!!$
!!$               If (abs(dconjg(phi(i_global,imo))).lt.1d-15) cycle
!!$
!!$               do j_fedvr3d =1,fedvr3d%nb_angle
!!$                  
!!$                  j_global = (j_fedvr3d-1)*fedvr3d%nb_r+global_to_local(i_fedvr3d,1)
!!$                  
!!$                  If (abs(dconjg(phi(i_global,imo))*phi(j_global,jmo)).lt.1d-15) cycle !! Juan: change 13 november AQUI
!!$
!!$                  sumtemp = sumtemp + dconjg(phi(i_global,imo))*phi(j_global,jmo)*cww_calc(i_fedvr3d,j_fedvr3d,kmo,lmo)           
!!$
!!$               enddo
!!$            enddo
!!$            tei_spatial(imo,kmo,lmo,jmo) = sumtemp  
!!$
!!$!! calculating the elements related            
!!$
!!$            tei_spatial(kmo,imo,jmo,lmo)= tei_spatial(imo,kmo,lmo,jmo)
!!$            tei_spatial(lmo,jmo,imo,kmo) = conjg(tei_spatial(imo,kmo,lmo,jmo))
!!$            tei_spatial(jmo,lmo,kmo,imo)= conjg(tei_spatial(imo,kmo,lmo,jmo))
!!$
!!$         enddo
!!$      enddo
!!$   enddo
!!$enddo

!! Juan: I have put this. New version done in 14 november

do i_fedvr3d =1,fedvr3d%nb_r*fedvr3d%nb_angle
   do j_fedvr3d =1,fedvr3d%nb_angle
      
      i_global = i_fedvr3d
      j_global = (j_fedvr3d-1)*fedvr3d%nb_r+global_to_local(i_fedvr3d,1)
            
      do imo =1, system%nptot !! Run in orbitals

         If (abs(phi(i_global,imo)).lt.1d-15) cycle !! Juan: change 14 november AQUI

         do jmo =1, system%nptot !! Run in orbitals
            If (abs(dconjg(phi(i_global,imo))*phi(j_global,jmo)).lt.1d-15) cycle !! Juan: change 13 november AQUI
            do kmo =imo, system%nptot !! Run in orbitals      
               do lmo = 1,system%nptot !! Run in orbitals
                  
                  sumtemp = zzero                 
                  
                  sumtemp = dconjg(phi(i_global,imo))*phi(j_global,jmo)*cww_calc(i_fedvr3d,j_fedvr3d,kmo,lmo)           
                 
                  tei_spatial(imo,kmo,lmo,jmo) = tei_spatial(imo,kmo,lmo,jmo)+sumtemp  
 
               enddo
            enddo
            
         enddo
      enddo
   enddo
enddo

do imo =1, system%nptot !! Run in orbitals
   do kmo =imo, system%nptot !! Run in orbitals
      do jmo =1, system%nptot !! Run in orbitals
         do lmo = 1,system%nptot !! Run in orbitals

!! calculating the elements related            
            tei_spatial(kmo,imo,jmo,lmo)= tei_spatial(imo,kmo,lmo,jmo)
            tei_spatial(lmo,jmo,imo,kmo) = conjg(tei_spatial(imo,kmo,lmo,jmo))
            tei_spatial(jmo,lmo,kmo,imo)= conjg(tei_spatial(imo,kmo,lmo,jmo))

         End do
      End do
   End do
End do

End Subroutine two_body_orbitals

!! End of the subroutines to avoid the storing of the two body matrix

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
    if(prop%stiffness==0) then
     call update_ham_onebody_spatial_space()
!=================================================================================
! the stiffness are considered
!=================================================================================
    else
     call update_ham_onebody_spatial_space2()  
    endif

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







