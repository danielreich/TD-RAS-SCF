
   subroutine Qspace()
     use global
     use operator_radial
     use density
     use wfunction
     use operator_spatial_space
     use operator_3d
     use twoe_basis_set
     implicit none
      integer :: ij,il,ik,ii,ip,ix,iy,iip,icoloum,i_spatial
      complex(kind=k2) :: cs,cff,csum, cww_auxiliar
      complex(kind=k2),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1,cbox2
      integer :: idicp,jdicp
      integer :: indexp,k1p,k2p,index_here
      integer :: iy_here !! global basis (FEDVR+ANGULAR)
      integer :: iy_ang !! staring angular for the 
   
!!===========================================================
!! This is done by Juan

      cbox1=zzero

      do i_spatial =1,system%nptot
         !
         ! not consider the stiffness problem in the calculations
         !
         if(prop%stiffness==0) then       
            
            call act(phi(:,i_spatial),cbox2(:,i_spatial))
            
         else
            call act2(phi(:,i_spatial),cbox2(:,i_spatial))
            
         endif
      End do

      do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle
       
         do iy = 1,fedvr3d%nb_angle!! 18/nov/2014 1,fedvr3d%nb_angle !! Juan: I have put it here!! CHECK HERE
            
            select case (fedvr3d%store)  !! Juan: case to find the FEDVR+angular of the second spatial function!! Juan: I have put it here
               
            case(0)  !! Juan: I have put it here
               iy_here = (iy-1)*fedvr3d%nb_r+global_to_local(ix,1) !! iy_here is the index in the FEDVR global basis             !! Juan: I have put it here
                 
            case(1)  !! Juan: I have put it here
               iy_here = index_two(2,ix,iy)!! Juan: I have put it here
               
              end select!! Juan: I have put it here
              
              do ij=1,system%nptot  !! Juan: I have put it here
                 
                 If (abs(phi(iy_here,ij)).lt.1d-20) cycle !! Juan: I have put it here
                 
                 do il=1,system%nptot 
                    do ik=1,system%nptot !! Juan: it can be reduced AQUI?
                       cww_auxiliar=cww_calc(ix,iy,ik,il) !! here we calculate the cww_calc which are common for the Q-space equations
                       If (abs(cww_auxiliar).lt.1d-20) cycle !! Juan: I have put it here
                       do i_spatial =1,system%nptot                        
                          
                          select case (fedvr3d%store)
                          case(0) 

!!$                             If (abs(cden3(i_spatial,ik,il,ij)*phi(iy_here,ij)).lt.1d-15) cycle!! change 12 november AQUI
                             
!!$                          cs=cs+cden3(i_spatial,ik,il,ij)*cww_calc(ix,iy,ik,il)*phi(iy_here,ij)
                             cbox1(ix,i_spatial)=cbox1(ix,i_spatial)+cden3(i_spatial,ik,il,ij)*cww_auxiliar*phi(iy_here,ij) !! Juan: I put this here. 17 november
               
                                                         
                          case(1)                             !! if cww is store
                             
                             iy_here = index_two(2,ix,iy)
                             cbox1(ix,i_spatial)=cbox1(ix,i_spatial)+cden3(i_spatial,ik,il,ij)*cww(ix,iy,ik,il)*phi(iy_here,ij)
                             
                          end select
                       enddo
                    enddo
                 enddo
              enddo
              
           enddo
        enddo
                
cbox1=cbox1+ cbox2

    
!! End of This is done by Juan
!!===========================================================

!!$     do i_spatial =1,system%nptot
!!$
!!$!
!!$! not consider the stiffness problem in the calculations
!!$!
!!$     if(prop%stiffness==0) then       
!!$
!!$        call act(phi(:,i_spatial),cbox2(:,i_spatial))
!!$
!!$     else
!!$        call act2(phi(:,i_spatial),cbox2(:,i_spatial))
!!$
!!$     endif
!!$     
!!$     do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle
!!$        cs=zzero        
!!$           do iy = 1,fedvr3d%nb_angle !! Juan: I have put it here
!!$              select case (fedvr3d%store)  !! Juan: case to find the FEDVR+angular of the second function!! Juan: I have put it here
!!$              case(0)!! Juan: I have put it here
!!$                 iy_here = (iy-1)*fedvr3d%nb_r+global_to_local(ix,1) !! iy_here is the index in the FEDVR global basis             !! Juan: I have put it here
!!$              case(1)!! Juan: I have put it here
!!$                 iy_here = index_two(2,ix,iy)!! Juan: I have put it here
!!$              end select!! Juan: I have put it here
!!$
!!$              do ij=1,system%nptot  !! Juan: I have put it here
!!$
!!$                 If (abs(phi(iy_here,ij)).lt.1d-15) cycle !! Juan: I have put it here
!!$
!!$                 do il=1,system%nptot 
!!$                    do ik=1,system%nptot  !! Wenliang has put it here
!!$              do ij=1,system%nptot  !! Wenliang has put it here
!!$                 do iy = 1,fedvr3d%nb_angle !! Wenliang has put it here
!!$                       select case (fedvr3d%store)
!!$                       case(0) 
!!$                          iy_here = (iy-1)*fedvr3d%nb_r+global_to_local(ix,1) !! iy_here is the index in the FEDVR global basis   
!!$                          If (abs(cden3(i_spatial,ik,il,ij)*phi(iy_here,ij)).lt.1d-15) cycle!! change 12 november AQUI
!!$                          
!!$                          cs=cs+cden3(i_spatial,ik,il,ij)*cww_calc(ix,iy,ik,il)*phi(iy_here,ij)
!!$                          cs=cs+cden3(i_spatial,ik,il,ij)*cww_auxiliar*phi(iy_here,ij) !! Juan: I put this here. 17 november
!!$                          
!!$                       case(1) !! iif cww is store
!!$                          iy_here = index_two(2,ix,iy)
!!$                          cs=cs+cden3(i_spatial,ik,il,ij)*cww(ix,iy,ik,il)*phi(iy_here,ij)
!!$                       end select
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$           enddo
!!$           cbox1(ix,i_spatial)= cs + cbox2(ix,i_spatial)
!!$        enddo
!!$     enddo
!
     ! the last step, projcetor acting
!
   call projector_Q(cbox1)

   return
  end subroutine Qspace
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

   subroutine projector_Q(cbox1)
   use global
   use wfunction
   use solver
   implicit none
   integer :: idicp,jdicp,kdicp,ldicp
   integer,dimension(system%nptot) :: ipiv2
   complex(8),dimension(system%nptot,system%nptot) :: phij,inv_phij
   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle) :: proj_Q
   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1
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
         If (abs(inv_phij(idicp,jdicp)).lt.1d-20) cycle !! Juan: I have put this
           do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
              If (abs(phi(kdicp,idicp)).lt.1d-15) cycle !! Juan: I have put this
            do ldicp =1,fedvr3d%nb_r*fedvr3d%nb_angle

            proj_Q(kdicp,ldicp) =  proj_Q(kdicp,ldicp) + phi(kdicp,idicp)*inv_phij(idicp,jdicp)*dconjg(phi(ldicp,jdicp)) 

           enddo
         enddo
         
        enddo
      enddo
!
! projector operator act on the intermediate wavefunction
!
      do idicp =1,system%nptot
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
            kphi(ldicp,idicp) =  (-ci*(cbox1(ldicp,idicp) - temp(ldicp)))
         enddo
         
      enddo

    return
 end subroutine projector_Q









!===============================================================================================
! not considered the stiffness effect,  the matrix is hightly sparse
!===============================================================================================
subroutine act(phi_in,phi_out)
use global
use operator_radial
use operator_3d
implicit none
integer :: idicp,jdicp,kdicp,ldicp
complex(kind=k2) :: phi_in(fedvr3d%nb_r,fedvr3d%nb_angle),phi_out(fedvr3d%nb_r,fedvr3d%nb_angle)
complex(kind=k2) :: temp1,temp2
integer :: mdicp,ndicp

!!AQUI REVISAR QUE ESTA BIEN
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
end subroutine act


!===============================================================================================
! considered the stiffness effect,  the matrix is full with element
!===============================================================================================
subroutine act2(phi_in,phi_out)
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
end subroutine act2
















