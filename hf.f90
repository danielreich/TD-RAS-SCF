
  subroutine hf_calc()
  use global
  use operator_3d
  implicit none
  integer :: idicp,jdicp,kdicp,ldicp,ierr,ilamda,isigma
  real(kind=k1) :: rsum,global_sum,two1,two2,energy 
  integer :: imiu,iv,k_imiu,l_imiu,m_imiu,k_ilamda,l_ilamda,m_ilamda,k_isigma,l_isigma,m_isigma
  integer :: k_v,l_v,m_v,i_cycle,ii1
  integer :: isum,i_k,k,l,m,nb_angle,imin,il,num_lm,im

  real(kind=k1),allocatable,save :: density_alpha(:,:)
  real(kind=k1),allocatable,save :: fock(:,:),h_core(:,:)
  real(kind=k1),allocatable,private :: cwave(:,:),cwave_old(:,:),w_temp(:),z_temp(:,:),fv1p(:),fv2p(:)






  allocate( h_core(numbasis,numbasis),density_alpha(numbasis,numbasis),fock(numbasis,numbasis),&
  cwave(numbasis,system1%numelectron),cwave_old(numbasis,system1%numelectron),w_temp(numbasis),z_temp(numbasis,numbasis),&
  fv1p(numbasis),fv2p(numbasis))


  !allocate(h_core(fedvr3d%nb_r*fedvr3d%nb_angle, fedvr3d%nb_r*fedvr3d%nb_angle))
  allocate(h_small(fedvr3d%nb_r,fedvr3d%nb_r))  
  allocate(w_small(fedvr3d%nb_r,fedvr3d%nb_angle))
  allocate(z_small(fedvr3d%nb_r,fedvr3d%nb_r,fedvr3d%nb_angle))
  allocate(fv1_small(fedvr3d%nb_r),fv2_small(fedvr3d%nb_r))

  allocate(energy_small(fedvr3d%nb_r*fedvr3d%nb_angle))
  allocate(index_b(fedvr3d%nb_r*fedvr3d%nb_angle,2)) 
  allocate(index_small(fedvr3d%nb_r*fedvr3d%nb_angle))
  allocate()




   kdicp = 0
   do idicp =1,fedvr3d%nb_angle
     do jdicp =1,fedvr3d%nb_r
      kdicp = kdicp+1

      index_small(kdicp) = kdicp 
      index_b(kdicp,1) = idicp
      index_b(kdicp,2) = jdicp
     enddo
  enddo
    



!================================================================================================================
! diag. the l-dependent onebody hamitoinan
!================================================================================================================
    do i_angle =1, fedvr3d%nb_angle
      h_small = zero
     do i_r = 1,n_total_kinetic

      i_row_r = index_kinetic_basis(i_r,1)
      i_column_r = index_kinetic_basis(i_r,2)

      if(i_row_r == i_column_r) then
         h_small(i_row_r,i_column_r)=(tmat_3d(i_row_r,i_angle) + vmat_radial(i_row_r) + tmat_radial(i_r))
      endif

      if(i_row_r /= i_column_r) then
        h_small(i_row_r,i_column_r) =  tmat_radial(i_r)
        h_small(i_column_r,i_row_r) =  tmat_radial(i_r)
      endif

      enddo
      call rs(fedvr3d%nb_r,fedr3d%nb_r,h_small,energy_small(  (i_angle-1)*fedvr3d%nb_r + 1  ),1,z_small(:,:,i_angle),&
           fv1_small,fv2_small,ierr)   
     enddo
!
! after diag. the ham., then sort the energy value in ascending order
!
  
    call sort2( fedvr3d%nb_r*fedvr3d%nb_angle,  energy_small, index_small)
!================================================================================================================
! initial the orbital 
!================================================================================================================     
    cwave_old = zero

    do idicp =1,system%numelectron/2

      jdicp = index_small(idicp)
       
      which_angle = index_b(jdicp,1)
      which_r = index_b(jdicp,2)
      ii1 = 0 
     do i1 =(which_angle - 1)*fedvr3d%nb_r+1, (which_angle - 1)*fedvr3d%nb_r + fedvr3d%nb_r 
       ii1 = ii1 + 1
       cwave_old(i1,idicp) = z_small(ii1,which_r,which_angle)
     enddo
   enddo
    

!
! begin the hf-cycle calculations
!
   i_cycle = 0 
   global_sum = 100.0d0
   do while(global_sum>0.0000001d0) 

   i_cycle = i_cycle + 1

   do idicp =1,numbasis
     do jdicp =1,numbasis

       rsum = 0.0d0
       do kdicp =1, system1%numelectron/2
  
         rsum = rsum + cwave_old(idicp,kdicp)*cwave_old(jdicp,kdicp)
        enddo
       density_alpha(idicp,jdicp) = 2.0d0*rsum
      enddo
  enddo


!
! two electron part
!
 do imiu =1, fedvr3d%nb_r*fedvr3d%nb_angle
    do iv =1, fedvr3d%nb_r*fedvr3d%nb_angle
     k_imiu = global_to_local(imiu,1)
     l_imiu = global_to_local(imiu,2)
     m_imiu = global_to_local(imiu,3)

     k_v = global_to_local(iv,1)
     l_v = global_to_local(iv,2)
     m_v = global_to_local(iv,3)

     rsum = 0.0d0

      do ilamda=1, fedvr3d%nb_r*fedvr3d%nb_angle
        do isigma =1,fedvr3d%nb_r*fedvr3d%nb_angle


         k_ilamda = global_to_local(ilamda,1)
         l_ilamda = global_to_local(ilamda,2)
         m_ilamda = global_to_local(ilamda,3)

         k_isigma = global_to_local(isigma,1)
         l_isigma = global_to_local(isigma,2)
         m_isigma = global_to_local(isigma,3)
        
       two1 = 0.0d0
       if(k_imiu == k_v .and. k_ilamda==k_isigma) then
          two1 =  fedvr3d_twoe(k_imiu,l_imiu,m_imiu,l_v,m_v,k_ilamda,l_ilamda,m_ilamda,l_isigma,m_isigma)
       endif

       two2 = 0.0d0
       if(k_imiu == k_ilamda .and. k_v==k_isigma) then
          two2 = fedvr3d_twoe(k_imiu,l_imiu,m_imiu,l_ilamda,m_ilamda,k_v,l_isigma,m_isigma,l_v,m_v)
       endif 

       rsum = rsum + density_alpha(ilamda,isigma) * (two1 - two2*0.5d0)

      enddo
    enddo
    fock(imiu,iv ) = rsum + dreal(h_basis(imiu,iv))
  enddo
enddo

 rsum = 0.0d0
 do idicp =1,numbasis !system1%numelectron/2
   do jdicp=1,numbasis !system1%numelectron/2
    rsum = rsum + density_alpha(jdicp,idicp)*(dreal(h_basis(idicp,jdicp)) + fock(idicp,jdicp)     )
   enddo
enddo

  energy = rsum*0.5d0

 
 write(*,*) i_cycle,energy
 

call rs(numbasis,numbasis,fock,w_temp,1,z_temp,fv1p,fv2p,ierr)

!write(*,*) w_temp(1:10)





cwave(:,1:system1%numelectron) = z_temp(:,1:system1%numelectron) 

rsum = 0.0d0
do idicp =1,system1%numelectron 
 do jdicp =1,numbasis
  
  rsum = rsum + abs(cwave_old(jdicp,idicp) - cwave(jdicp,idicp))
  enddo
 enddo

!   do idicp =1,numbasis
 !    do jdicp =1,numbasis

  !     rsum = 0.0d0
   !        do kdicp =1, system1%numelectron/2
  
    !     rsum = rsum + cwave(idicp,kdicp)*cwave(jdicp,kdicp)
     !      enddo
      ! density_alpha(idicp,jdicp) = 2.0d0*rsum
      !   enddo
  !enddo






! rsum = 0.0d0
! do idicp =1,system1%numelectron/2
!   do jdicp=1,system1%numelectron/2
!    rsum = rsum + density_alpha(jdicp,idicp)*(dreal(h_basis(idicp,jdicp)) + fock(idicp,jdicp)     )
!   enddo
!enddo

!  energy = rsum*0.5d0


 global_sum = rsum
 cwave_old = cwave

enddo


 
 do idicp =1,numbasis

   if(w_temp(idicp)<0.0d0) then
     m_auto = idicp

   else
     exit
   endif
 enddo


 m_auto =10

 allocate(phi_ionization(numbasis,m_auto))


 phi_ionization(:,1:m_auto) = z_temp(:,1:m_auto)


end subroutine hf_calc


