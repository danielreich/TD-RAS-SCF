
   subroutine Qspace_fc()
     use global
     use operator_radial
     use density
     use wfunction
     use operator_spatial_space
     use operator_3d
     use twoe_basis_set
     implicit none
      integer :: ij,il,ik,ii,ip,ix,iy,iip,icoloum,i_spatial
      complex(kind=k2) :: cs,cff,csum
      complex(kind=k2),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1,cbox2
      integer :: idicp,jdicp
      integer :: indexp,k1p,k2p,index_here,iy_here 
   


     do i_spatial =1,system%nptot

!
! not consider the stiffness problem in the calculations
!
     if(prop%stiffness==0) then

      call act(phi(:,i_spatial),cbox1(:,i_spatial))

     else
      call act2(phi(:,i_spatial),cbox1(:,i_spatial))

     endif

    enddo




    do ix =1,fedvr3d%nb_r*fedvr3d%nb_angle
     do ii =1,np_tot_fc

       cs = zzero
       do ij =1,system%nptot
          cs = cs + cden1(ii+n_offset,ij)*cbox1(ix,ij)
       enddo
      cbox2(ix,ii) = cs
     enddo
    enddo

     
     do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle

      do ii =1,nptot_fc
       cs=zzero
        do il=1,system%nptot
         do ik=1,system%nptot
          do ij=1,system%nptot
           do iy = 1,fedvr3d%nb_angle
            iy_here = index_two(2,ix,iy)
            cs=cs+cden2(ii,ik,il,ij)*cww(ix,iy,ik,il)*phi(iy_here,ij)
           enddo
          enddo
         enddo
        enddo
        cbox3(ix,i_spatial)= cs + cbox2(ix,i_spatial)
      enddo
    enddo

!

  do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle
    do ii=1,nptot_fc
      cs = zzero
      do ij =1,nptot_fc
       cs = cs + cden1_reg_inv_fc(ii,ij)*cbox3(ix,ij)
      enddo
     cbox1(ix,ii) = cs

   enddo
   enddo


!
! the last step, projcetor acting
!
   call projector_Q_hf(cbox1)

   return
  end subroutine Qspace_fc
  !#######################################################################################
! projector operator in Q space, i is spatial orbital, n is the num of spatial orbital  
!                  _N___                 __N____
!                  \                     \ _N___
!				    \                     \\     
!             1 -   / | i > <i|  == 1 -   //___  |i > [O]_ij <j| 
!                  /___                  /____
!				   i=1                   i,j=1
! O matrix is the inverse of  <i|j> matrix, the reference can be found in 
! M.H. Beck, Physics Report, 324, 2000, Page 25                                         
!                                                                                       
!########################################################################################

   subroutine projector_Q_hf(cbox1)
   use global
   use wfunction
   use solver
   implicit none
   integer :: idicp,jdicp,kdicp,ldicp
   integer,dimension(system%nptot) :: ipiv2
   complex(8),dimension(system%nptot,system%nptot) :: phij,inv_phij
   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle) :: proj_Q
   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1,cbox3
   complex(8) :: sumtemp
   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle) :: temp
   integer :: info
    
 !
 ! construct the matrix  phij = < i | j >, i,j spatial orbital
 ! 
   do idicp =1,system%nptot
     do jdicp =1, system%nptot
         sumtemp = zzero

        do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle 
           sumtemp = sumtemp + dconjg(phi(kdicp,idicp))*phi(kdicp,jdicp)
        enddo
        phij(idicp,jdicp) = sumtemp
        inv_phij(idicp,jdicp) = zzero
        enddo
          inv_phij(idicp,idicp) = zone
        enddo  
!
! find the inverse matrix of phij
! 

    call zgesv(system%nptot,system%nptot,phij,system%nptot,ipiv2,inv_phij,system%nptot,info)
    if(info.ne.0)  then
        write(*,*) '____________________________________________________________'
        write(*,*) 'fail to find the inverse of phij <phi | phj> in projector_Q'
        write(*,*) '------------------------------------------------------------'
        stop 
    endif

!
! construct the Q.space Projector operator
!

    proj_Q = zzero
 
    do idicp =1,system%nptot
       do jdicp =1,system%nptot
         
           do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
            do ldicp =1,fedvr3d%nb_r*fedvr3d%nb_angle

            proj_Q(kdicp,ldicp) =  proj_Q(kdicp,ldicp) + phi(kdicp,idicp)*inv_phij(idicp,jdicp)*dconjg(phi(ldicp,jdicp)) 

           enddo
         enddo
         
        enddo
      enddo
!
! projector operator act on the intermediate wavefunction
!
      do idicp =1,nptot_fc
          do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
           sumtemp = zzero
                   do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
                     sumtemp = sumtemp + proj_Q(jdicp,kdicp)*cbox1(kdicp,idicp)
                   enddo
          
           temp(jdicp) =  sumtemp

         enddo
   
!
!  1  -  Projector
!

         do ldicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
           cbox3(ldicp,idicp) =  (-ci*(cbox1(ldicp,idicp) - temp(ldicp)))
         enddo

      enddo

    
    do ip=1,nptot_fc
      do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle

       kphi(ix,ip+n_offset) = cbox3(ix,ip)
      enddo
    enddo

    return
 end subroutine projector_Q_hf









!===============================================================================================
! not considered the stiffness effect,  the matrix is hightly sparse
!===============================================================================================
subroutine act_fc(phi_in,phi_out)
use global
use operator_radial
use operator_3d
implicit none
integer :: idicp,jdicp,kdicp,ldicp
complex(kind=k2) :: phi_in(fedvr3d%nb_r,fedvr3d%nb_angle),phi_out(fedvr3d%nb_r,fedvr3d%nb_angle)
complex(kind=k2) :: temp1,temp2
integer :: mdicp,ndicp



do idicp =1,fedvr3d%nb_angle
    do jdicp =1,fedvr3d%nb_r

      temp1 = phi_in(jdicp,idicp)*(vmat_radial(jdicp)  +  tmat_3d(jdicp,idicp)) 
      temp2 = zzero
       do kdicp =index_ham_act(jdicp,1),index_ham_act(jdicp,2)
          
         ldicp = index_kinetic_operator(max(jdicp,kdicp), min(jdicp,kdicp))
         temp2 = temp2 + tmat_radial(ldicp)*phi_in(kdicp,idicp) 
       enddo
     
    phi_out(jdicp,idicp)  = temp1 + temp2

   enddo
enddo



!=================================================================================================
! the laser interaction term
!=================================================================================================
if(laser%tdornot) then

 do idicp=1,fedvr3d%nb_angle
  do jdicp =1,fedvr3d%nb_r
   
    temp1 = zzero
    do ldicp =1,fedvr3d%nb_angle  
 
     temp1 = temp1 + laser%ez_t*zmat_3d(idicp,ldicp)*fedvrx_global(jdicp)*phi_in(jdicp,ldicp)
    enddo

  phi_out(jdicp,idicp) = phi_out(jdicp,idicp) + temp1  

 enddo
enddo

endif

return
end subroutine act_fc


!===============================================================================================
! considered the stiffness effect,  the matrix is full with element
!===============================================================================================
subroutine act2_fc(phi_in,phi_out)
use global
use operator_radial
use operator_3d
implicit none
integer :: idicp,jdicp,kdicp,ldicp
complex(kind=k2) :: phi_in(fedvr3d%nb_r,fedvr3d%nb_angle),phi_out(fedvr3d%nb_r,fedvr3d%nb_angle)
complex(kind=k2) :: temp1
 
 do idicp =1,fedvr3d%nb_angle
   do jdicp =1,fedvr3d%nb_r
     temp1 = zzero
     do kdicp =1,fedvr3d%nb_r


       temp1 = temp1 + h_stiffness(ia(max(jdicp,kdicp))+min(jdicp,kdicp),idicp)*phi_in(kdicp,idicp)
     enddo
  phi_out(jdicp,idicp) = temp1
 enddo

enddo


!=================================================================================================
! the laser interaction term
!=================================================================================================
if(laser%tdornot) then

 do idicp=1,fedvr3d%nb_angle
  do jdicp =1,fedvr3d%nb_r    
    temp1 = zzero
    do ldicp =1,fedvr3d%nb_angle

     temp1 = temp1 + laser%ez_t*zmat_3d(idicp,ldicp)*fedvrx_global(jdicp)*phi_in(jdicp,ldicp)
    enddo


  phi_out(jdicp,idicp) = phi_out(jdicp,idicp) + temp1

 enddo
enddo

endif

return
end subroutine act2_fc
















