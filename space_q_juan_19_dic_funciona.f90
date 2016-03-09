
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
!!$

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

      select case (fedvr3d%store)
      case(0)
         If (allocated(tei_spatial)) then
            tei_spatial=zzero
         else
            !! If tei_spatial is not allocated, do it!
            allocate(tei_spatial(1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))
            tei_spatial=zzero
         End If
      case default
         continue
      end select

      do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle
!!$         print*, ix !! Juan: Used to debug
         do iy = 1,fedvr3d%nb_angle!! 18/nov/2014 1,fedvr3d%nb_angle !! Juan: I have put it here!! CHECK HERE

            select case (fedvr3d%store)  !! Juan: case to find the FEDVR+angular of the second spatial function!! Juan: I have put it here
               
            case(0)  !! Juan: I have put it here
               iy_here = (iy-1)*fedvr3d%nb_r+global_to_local(ix,1) !! iy_here is the index in the FEDVR global basis             !! Juan: I have put it here
                 
            case(1)  !! Juan: I have put it here

               iy_here = index_two(2,ix,iy)!! Juan: I have put it here
               
              end select!! Juan: I have put it here
              
              do ij=1,system%nptot  !! Juan: I have put it here

                 If (abs(phi(iy_here,ij)).lt.1d-15) cycle !! Juan: I have put it here
                 
                 do il=1,system%nptot 

                    do ik=1,system%nptot 
                                         
                       if (fedvr3d%store.eq.0) then
                          cww_auxiliar=cww_calc(ix,iy,ik,il)*phi(iy_here,ij) !! here we calculate the cww_calc which are common for the Q-space equations
                          If (abs(cww_auxiliar).lt.1d-20) cycle !! Juan: I have put it here
                       end if

                       do i_spatial =1,system%nptot                        
!!$                           If (abs(cden3(i_spatial,ik,il,ij)).lt.1d-15) cycle!! change 12 november AQUI            

                          select case (fedvr3d%store)
                          case(0) 
                                                         
!!$                          cs=cs+cden3(i_spatial,ik,il,ij)*cww_calc(ix,iy,ik,il)*phi(iy_here,ij)
                             cbox1(ix,i_spatial)=cbox1(ix,i_spatial)+cden3(i_spatial,ik,il,ij)*cww_auxiliar!! Juan: I put this here. 17 november               

!! JUAN: HERE WE CALCULATE THE TWO BODY INTEGRALS AT THE SAME TIME WE SOLVE THE Q EQUATIONS
                             tei_spatial(i_spatial,ik,il,ij)=tei_spatial(i_spatial,ik,il,ij)+dconjg(phi(ix,i_spatial))*cww_auxiliar
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
     ! the last step, projector acting
!
   call projector_Q(cbox1)

   return
  end subroutine Qspace


!!=============================================================================================
!! Subroutine Qspace_opt is a subroutine to solve the Q space equations and calculate the two
!! body integrals without calling the subroutines which calculates the mean field operators
!!=============================================================================================

   subroutine Qspace_opt()
     use global
    use operator_radial
     use density
     use wfunction
     use operator_spatial_space
     use operator_3d
     use twoe_basis_set
     implicit none

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

      integer :: ij,il,ik,ii,ip,ix,iy,iip,icoloum,i_spatial
      complex(kind=k2) :: cs,cff,csum, cww_auxiliar
      complex(kind=k2),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1,cbox2
      complex(kind=k2), allocatable :: cbox1_proc(:,:,:) !! used to paralellize
      complex(kind=k2),allocatable,save :: tei_spatial_proc(:,:,:,:,:)
      integer :: proc !! number of processors
      integer :: num_proc !! number of process
      integer :: idicp,jdicp
      integer :: indexp,k1p,k2p,index_here, kdicp, kdicp_radial, kdicp_ang,ldicp_ang
      integer :: n1, n2, n3, ll, l1,l2,l3,l4,m1,m2,m3,m4, ldicp
      integer :: iy_here !! global basis (FEDVR+ANGULAR)
      integer :: iy_ang !! staring angular for the 
      integer :: ang1,ang2,ang3,ang4
      integer :: n1_kdicp_radial !! This variable runs for n1*kdicp
      integer :: orb_tot
      real (kind=k1) :: twoe_fedvr_angle,twoe_fedvr_aux,twoe_fedvr_total
   
!!===========================================================
!! This is done by Juan

!!$      call omp_set_num_threads(1)

      cbox1=zzero

      do i_spatial =1,system%nptot
         !
         ! not consider the stiffness problem in the calculations
         !
         if(prop%stiffness==0) then       
            call act(phi(:,i_spatial),cbox2(:,i_spatial))
         else
            call act2(phi(:,i_spatial),cbox2(:,i_spatial))
!!$            call act3(phi(:,i_spatial),cbox2(:,i_spatial)). made by Juan: This subroutine uses the expansion
!            in the eigenvectors of the field-free one body hamitonian.
         endif
      End do     

!! Calculate the two body contributions

      If (allocated(tei_spatial)) then
         tei_spatial=zzero
      else
         !! If tei_spatial is not allocated, do it!
         allocate(tei_spatial(1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))
      End If

!! SET MAX NUM OF THREADS

      proc=OMP_get_num_procs()

!!$      call OMP_set_num_threads(proc)

         if (allocated(tei_spatial_proc).and.allocated(cbox1_proc)) then
            continue
         else
            allocate(tei_spatial_proc(0:proc-1,1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))   !! to parallelize
            allocate(cbox1_proc(0:proc-1,fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot))!! for the paralellization            
         !! this is the contribution of each processor
         End if
         cbox1_proc=zzero
         tei_spatial=zzero
         tei_spatial_proc=zzero    

            !!PARALLELIZE HERE
!$OMP PARALLEL PRIVATE(kdicp_radial,ij,il,ik,i_spatial,ang1,iy,l1,m1,l2,m2,kdicp_ang,l3,m3,ldicp_ang,l4,m4,ll,twoe_fedvr_total,n1,ii)!n1_kdicp_radial)! Before ii (runing in angular)
!!PRIVATE(kdicp_radial,ij,il,ik,i_spatial,ang1,iy,l1,m1,l2,m2,kdicp_ang,l3,m3,ldicp_ang,l4,m4,ll,twoe_fedvr_angle,n1_kdicp_radial) Before n1
!!PRIVATE(kdicp_radial,ij,il,ik,i_spatial) !! if it is in the angular part            
            !$OMP DO           

            Do n1_kdicp_radial= 1, fedvr3d%nb_r*(fedvr3d%nb_r+1)/2 !! running in the radial part without nested loops            

               Do ii=1,size(twoe_ang1(:)) !! Running for the combinations of maximum angular momentum

               n1=twoe_radial_r1(n1_kdicp_radial) !! radial part of r1
               kdicp_radial=twoe_radial_r2(n1_kdicp_radial) !! radial part of r2
               
               ang1=twoe_ang1(ii)             !! first angular function
               
               iy=twoe_ang2(ii)               !! second angular function
               
               kdicp_ang=twoe_ang3(ii)        !! third angular function
               
               ldicp_ang=twoe_ang4(ii)        !! fourth angular function            
               
               l1=lm_l(ang1)  !! total angular momentum of the first basis function
               l2=lm_l(iy)    !! total angular momentum of the second basis function
               l3=lm_l(kdicp_ang)   !! total angular momentum of the third basis function
               l4=lm_l(ldicp_ang)   !! total angular momentum of the fourth basis function
               
               twoe_fedvr_total=zero  !initialize the value of the integral
               
!!$               twoe_fedvr_total=dot_product(twoe_radial_store(n1,kdicp_radial,:),twoe_integral_ll(ii,:))
!!$
               Do ll=max(abs(l1-l2),abs(l3-l4)),min(l1+l2,l3+l4)
                  twoe_fedvr_total=twoe_fedvr_total+twoe_radial_store(n1,kdicp_radial,ll)*twoe_integral_ll(ii,ll) !! contribution for each element
               End Do
               
               if (abs(twoe_fedvr_total).lt.1d-15) cycle 
               
               do ij=1,system%nptot  !! Run in the orbitals
                  If (abs(phi((iy-1)*fedvr3d%nb_r+n1,ij)).lt.1d-15.and.abs(phi((iy-1)*fedvr3d%nb_r+kdicp_radial,ij)).lt.1d-15) cycle !! Juan: I have put it here
                  
                  do il=1,system%nptot !! Run in the orbitals
!!$                                 If (abs(phi(ldicp,il)).lt.1d-15) cycle
                     If (abs(phi((ldicp_ang-1)*fedvr3d%nb_r+kdicp_radial,il)).lt.1d-15.and.abs(phi((ldicp_ang-1)*fedvr3d%nb_r+n1,il)).lt.1d-15) cycle
                     
                     do ik=1,system%nptot !! Run in the orbitals
!!$                                    If (abs(phi(kdicp,ik)).lt.1d-15) cycle
                        If (abs(phi((kdicp_ang-1)*fedvr3d%nb_r+kdicp_radial,ik)).lt.1d-15.and.abs(phi((kdicp_ang-1)*fedvr3d%nb_r+n1,ik)).lt.1d-15) cycle
                        
                        do i_spatial =1,system%nptot                    
                           
                           !!  HERE WE CALCULATE THE TWO BODY INTEGRALS AT THE SAME TIME WE SOLVE THE Q EQUATIONS
                           
                           cbox1_proc(OMP_get_thread_num(),(ang1-1)*fedvr3d%nb_r+n1,i_spatial)=&
                                cbox1_proc(OMP_get_thread_num(),(ang1-1)*fedvr3d%nb_r+n1,i_spatial)+&
                                cden3(i_spatial,ik,il,ij)*&
                                dconjg(phi((kdicp_ang-1)*fedvr3d%nb_r+kdicp_radial,ik))*&
                                phi((ldicp_ang-1)*fedvr3d%nb_r+kdicp_radial,il)*&
                                twoe_fedvr_total*&!twoe_radial_store(n1,kdicp_radial,ll)*twoe_fedvr_angle*&
                                phi((iy-1)*fedvr3d%nb_r+n1,ij) !! Term corresponding to the mean field operator in Q-space equations 
                           
          
                        
                           tei_spatial_proc(OMP_get_thread_num(),i_spatial,ik,il,ij)=&
                                tei_spatial_proc(OMP_get_thread_num(),i_spatial,ik,il,ij)+&
                                dconjg(phi((ang1-1)*fedvr3d%nb_r+n1,i_spatial))&
                                *dconjg(phi((kdicp_ang-1)*fedvr3d%nb_r+kdicp_radial,ik))&
                                *phi((ldicp_ang-1)*fedvr3d%nb_r+kdicp_radial,il)*&
                                twoe_fedvr_total*&!*twoe_radial_store(n1,kdicp_radial,ll)*twoe_fedvr_angle*&
                                phi((iy-1)*fedvr3d%nb_r+n1,ij) !! Calculation of the two body operator, which will be used in 
                           
                           If (n1.ne.kdicp_radial) then
                              
                              cbox1_proc(OMP_get_thread_num(),(ang1-1)*fedvr3d%nb_r+kdicp_radial,i_spatial)=&
                                   cbox1_proc(OMP_get_thread_num(),(ang1-1)*fedvr3d%nb_r+kdicp_radial,i_spatial)+&
                                   cden3(i_spatial,ik,il,ij)*&
                                   dconjg(phi((kdicp_ang-1)*fedvr3d%nb_r+n1,ik))*&
                                   phi((ldicp_ang-1)*fedvr3d%nb_r+n1,il)*&
                                   twoe_fedvr_total*&!twoe_radial_store(n1,kdicp_radial,ll)*twoe_fedvr_angle*&
                                   phi((iy-1)*fedvr3d%nb_r+kdicp_radial,ij) !! Term corresponding to the mean field operator in Q-space equations 
                              
                              tei_spatial_proc(OMP_get_thread_num(),i_spatial,ik,il,ij)=&
                                   tei_spatial_proc(OMP_get_thread_num(),i_spatial,ik,il,ij)+&
                                   dconjg(phi((ang1-1)*fedvr3d%nb_r+kdicp_radial,i_spatial))&
                                   *dconjg(phi((kdicp_ang-1)*fedvr3d%nb_r+n1,ik))&
                                   *phi((ldicp_ang-1)*fedvr3d%nb_r+n1,il)*&
                                   twoe_fedvr_total*&!*twoe_radial_store(n1,kdicp_radial,ll)*twoe_fedvr_angle*&
                                   phi((iy-1)*fedvr3d%nb_r+kdicp_radial,ij) !! Calculation of the two body operator, which will be used in 
                              
                           End If
                                                                        
!!$                                       tei_spatial(i_spatial,ik,il,ij)=tei_spatial(i_spatial,ik,il,ij)+&
!!$                                            dconjg(phi((ang1-1)*fedvr3d%nb_r+n1,i_spatial))&
!!$                                            *dconjg(phi((kdicp_ang-1)*fedvr3d%nb_r+kdicp_radial,ik))&
!!$                                            *phi((ldicp_ang-1)*fedvr3d%nb_r+kdicp_radial,il)&
!!$                                            *twoe_radial_store(n1,kdicp_radial,ll)*twoe_fedvr_angle*&
!!$                                            phi((iy-1)*fedvr3d%nb_r+n1,ij) !! Calculation of the two body operator, which will be used in 
                        end do
                     enddo
                  enddo
               enddo
            enddo
         End Do
            
            !$OMP END DO
      !$OMP END PARALLEL


            Do num_proc=0,OMP_get_num_procs()-1
               !! sum all the contributions
               tei_spatial(:,:,:,:)=tei_spatial(:,:,:,:)+tei_spatial_proc(num_proc,:,:,:,:)
               cbox1(:,:)=cbox1(:,:)+cbox1_proc(num_proc,:,:)
            End Do

            tei_spatial_proc=zzero
            cbox1_proc=zzero

!!$                     enddo
!!$                  enddo
!!$               end do
!!$            end do                     
    
            deallocate(tei_spatial_proc,cbox1_proc)

            cbox1=cbox1+ cbox2

!!===========================================================================================    
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
     ! the last step, projector acting
!

!!$   call projector_Q_opt(cbox1)
            print*, 'hola'
            call projector_Q(cbox1)

   return
 end subroutine Qspace_opt

!!=============================================================================================
!! Subroutine Qspace_omp is the parallelized version of the subroutine Qspace
!!=============================================================================================

   subroutine Qspace_omp()
     use global
     use operator_radial
     use density
     use wfunction
     use operator_spatial_space
     use operator_3d
     use twoe_basis_set
     implicit none
      integer :: ij,il,ik,ii,ip,ix,iy,iip,icoloum,i_spatial
      integer :: ix_y
      integer :: i_spatial_init, ij_init
      complex(kind=k2) :: cs,cff,csum, cww_auxiliar
      complex(kind=k2),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1,cbox2
      integer :: idicp,jdicp
      integer :: indexp,k1p,k2p,index_here
      integer :: iy_here !! global basis (FEDVR+ANGULAR)
      integer :: iy_ang !! staring angular for the 
      integer :: proc,num_proc !! parameter of the parallelization
      complex(kind=k2), allocatable :: cbox1_proc(:,:,:) !! used to paralellize
      complex(kind=k2),allocatable,save :: tei_spatial_proc(:,:,:,:,:)
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
   
!!===========================================================
!! This is done by Juan

!!$      call OMP_set_nested(.TRUE.)

      cbox1=zzero
!!$

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

      select case (fedvr3d%store)
      case(0)
         If (allocated(tei_spatial)) then
            tei_spatial=zzero
         else
            !! If tei_spatial is not allocated, do it!
            allocate(tei_spatial(1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))
            tei_spatial=zzero
         End If
      case default
         continue
      end select


!! SET MAX NUM OF THREADS
!!$      call OMP_set_num_threads(8)

      proc=OMP_get_num_procs()

         if (allocated(tei_spatial_proc).and.allocated(cbox1_proc)) then
            continue
         else
            allocate(tei_spatial_proc(0:proc-1,1:system%nptot,1:system%nptot,1:system%nptot,1:system%nptot))   !! to parallelize
            allocate(cbox1_proc(0:proc-1,fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot))!! to paralellize           
         !! this is the contribution of each processor
         End if

         cbox1_proc=zzero
         tei_spatial=zzero
         tei_spatial_proc=zzero    

         !!PARALLELIZE HERE
         !$OMP PARALLEL PRIVATE(iy,iy_here,ij,il,ik,cww_auxiliar,i_spatial,i_spatial_init,ij_init) 
         !$OMP DO SCHEDULE(dynamic)!, 2200)
            
!!$         do ix_y=0,fedvr3d%nb_r*fedvr3d%nb_angle*fedvr3d%nb_angle-1
!!$
!!$            ix=mod(ix_y,fedvr3d%nb_r*fedvr3d%nb_angle)+1
!!$            iy=(ix_y-(ix-1))/(fedvr3d%nb_r*fedvr3d%nb_angle)+1


         do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle

            If (global_to_local(ix,2).gt.fedvr3d%l_max_core) then
               i_spatial_init=system%np0+system%np1+1
            Else
               i_spatial_init=1
            End If
            
            do iy = 1,fedvr3d%nb_angle!! 18/nov/2014 1,fedvr3d%nb_angle !! Juan: I have put it here!! CHECK HERE

               iy_here = (iy-1)*fedvr3d%nb_r+global_to_local(ix,1) !! iy_here is the index in the FEDVR global basis             !! Juan: I have put it here
                 
!! This if block takes into account that the core has maximum L, different to the other one
               If (global_to_local(iy_here,2).gt.fedvr3d%l_max_core) then
                  ij_init=system%np0+system%np1+1
               Else
                  ij_init=1
               End If

              do ij=ij_init,system%nptot  !! Juan: I have put it here

                 If (abs(phi(iy_here,ij)).lt.1d-10) cycle !! Juan: I have put it here
                
                 do il=1,system%nptot 

                    do ik=1,system%nptot                                

                       cww_auxiliar=cww_calc(ix,iy,ik,il)*phi(iy_here,ij) !! here we calculate the cww_calc which are common for the Q-space equations
                      
                       If (abs(cww_auxiliar).lt.1d-15) cycle !! Juan: I have put it here                    

                       do i_spatial =i_spatial_init,system%nptot                

                           cbox1_proc(OMP_get_thread_num(),ix,i_spatial)= cbox1_proc(OMP_get_thread_num(),ix,i_spatial)+&
                                cden3(i_spatial,ik,il,ij)*cww_auxiliar !! Term corresponding to the mean field operator in Q-space equations 

!! JUAN: HERE WE CALCULATE THE TWO BODY INTEGRALS AT THE SAME TIME WE SOLVE THE Q EQUATIONS
                                                  
                           tei_spatial_proc(OMP_get_thread_num(),i_spatial,ik,il,ij)=&
                                tei_spatial_proc(OMP_get_thread_num(),i_spatial,ik,il,ij)+&
                                dconjg(phi(ix,i_spatial))*cww_auxiliar

                        enddo
                    enddo
                 enddo
              enddo
              
           enddo
!!$           enddo
           enddo
        
            !$OMP END DO
      !$OMP END PARALLEL

            Do num_proc=0,proc-1
               !! sum all the contributions
               cbox1(:,:)=cbox1(:,:)+cbox1_proc(num_proc,:,:)
               tei_spatial(:,:,:,:)=tei_spatial(:,:,:,:)+tei_spatial_proc(num_proc,:,:,:,:)
            End Do          
    
           deallocate(tei_spatial_proc,cbox1_proc)

           cbox1=cbox1+ cbox2


   call projector_Q_omp(cbox1)


   return
 end subroutine Qspace_omp

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

!! PARALLELIZED VERSION OF PROJECTOR_Q

   subroutine projector_Q_omp(cbox1)
   use global
   use wfunction
   use solver
   implicit none
   integer :: idicp,jdicp,kdicp,ldicp,kdicp_end
   integer :: idicp_init, jdicp_init,ldicp_end
   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,system%nptot) :: cbox1 !! it is an input
   integer, allocatable :: ipiv2(:)
   complex(8),allocatable :: phij(:,:),inv_phij(:,:)
   complex(8),allocatable :: proj_Q(:,:)
!!$   integer,dimension(system%nptot) :: ipiv2
!!$   complex(8),dimension(system%nptot,system%nptot) :: phij,inv_phij
!!$   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle) :: proj_Q
   complex(8) :: sumtemp
!!$   complex(8),dimension(fedvr3d%nb_r*fedvr3d%nb_angle) :: temp
   integer :: info
!! To parallelize using OPENMP
   integer :: proc, num_proc !! processes and the number we are following
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

 !
 ! construct the matrix  phij = < i | j >, i,j spatial orbital
 ! 

     allocate(phij(system%nptot,system%nptot))
     allocate(inv_phij(system%nptot,system%nptot))

   do idicp =1,system%nptot
     do jdicp =1, system%nptot

!! To change the range for the maximum angular momentum for the core
         If (idicp.le.system%np0+system%np1) then
            kdicp_end=fedvr3d%nb_r*ang_max_core
            else
            kdicp_end=fedvr3d%nb_r*fedvr3d%nb_angle 
         End If

         If (jdicp.le.system%np0+system%np1) then
            kdicp_end=fedvr3d%nb_r*ang_max_core
            else
            kdicp_end=fedvr3d%nb_r*fedvr3d%nb_angle 
         End If

         sumtemp = zzero

        do kdicp =1,kdicp_end!fedvr3d%nb_r*fedvr3d%nb_angle 
           sumtemp = sumtemp + conjg(phi(kdicp,idicp))*phi(kdicp,jdicp)
        enddo

        phij(idicp,jdicp) = sumtemp
        inv_phij(idicp,jdicp) = zzero
        enddo
          inv_phij(idicp,idicp) = zone
       enddo

!
! find the inverse matrix of phij
! 

       allocate(ipiv2(system%nptot)) !! auxiliar variable
       ipiv2=0

    call zgesv(system%nptot,system%nptot,phij,system%nptot,ipiv2,inv_phij,system%nptot,info)

       deallocate(ipiv2)    

    if(info.ne.0)  then
        write(*,*) '____________________________________________________________'
        write(*,*) 'fail to find the inverse of phij <phi | phj> in projector_Q'
        write(*,*) '------------------------------------------------------------'
        stop 
    endif

!
! construct the Q.space Projector operator
!

    allocate(proj_Q(fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_r*fedvr3d%nb_angle))

    proj_Q = zzero

    do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
       If (global_to_local(kdicp,2).gt.fedvr3d%l_max_core) then
          idicp_init=system%np0+system%np1+1
       Else
          idicp_init=1
       End If
       do idicp =idicp_init,system%nptot
          If (abs(phi(kdicp,idicp)).lt.1d-15) cycle !! Juan: I have put this
          do jdicp =1,system%nptot
             If (abs(inv_phij(idicp,jdicp)).lt.1d-15) cycle !! Juan: I have put this
             If (jdicp.le.system%np0+system%np1+1) then !! this condition to set the maximun value of the core.
                ldicp_end=fedvr3d%nb_r*ang_max_core
             else
                ldicp_end=fedvr3d%nb_r*fedvr3d%nb_angle 
             End If
             do ldicp =1,ldicp_end!fedvr3d%nb_r*fedvr3d%nb_angle
                
                proj_Q(kdicp,ldicp) =  proj_Q(kdicp,ldicp) + phi(kdicp,idicp)*inv_phij(idicp,jdicp)*dconjg(phi(ldicp,jdicp)) 

             enddo
          enddo
       enddo
    enddo

!
! projector operator act on the intermediate wavefunction
!

kphi=cbox1 !! the vector which will store the projection is equal to the vector before project.


do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
 If (global_to_local(kdicp,2).gt.fedvr3d%l_max_core) then
          idicp_init=system%np0+system%np1+1
       Else
          idicp_init=1
       End If
   do idicp =idicp_init,system%nptot
      If (abs(cbox1(kdicp,idicp)).lt.1d-15) cycle
      do jdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle                  
         !
         !  1  -  Projector
         !         
         !! here we substract the projection         
         !! the second term is the projection, which is substracted
         kphi(jdicp,idicp) =  kphi(jdicp,idicp) - proj_Q(jdicp,kdicp)*cbox1(kdicp,idicp)
         
      enddo      
   enddo         
enddo

      kphi=-ci*kphi

      deallocate(proj_Q)

    return
  end subroutine projector_Q_omp


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


!===============================================================================================
! considered the stiffness effect, using Juan's method.
! We take the stiffness hamiltonian as the sum of the vectors times the
! corresponding 
!===============================================================================================
subroutine act3(phi_in,phi_out)
use global
use operator_radial
use operator_3d
implicit none
integer :: idicp,jdicp,kdicp,ldicp
complex(kind=k2),intent(in) :: phi_in(1:fedvr3d%nb_r*fedvr3d%nb_angle)
complex(kind=k2), intent(out) :: phi_out(1:fedvr3d%nb_r*fedvr3d%nb_angle)
complex(kind=k2) :: temp1

phi_out=zzero

Do idicp=1,size(z_band(1,:)) !! Run for all the eigenvectors

phi_out=phi_out+dot_product(z_band(:,idicp),phi_in)*w_band(idicp)*z_band(:,idicp)

End Do

!!$
!!$
!!$
!!$do idicp =1,fedvr3d%nb_angle
!!$   do jdicp =1,fedvr3d%nb_r
!!$      temp1 = zzero
!!$      do kdicp =1,fedvr3d%nb_r
!!$         
!!$         temp1 = temp1 + h_stiffness(ia(max(jdicp,kdicp))+min(jdicp,kdicp),idicp)*phi_in(kdicp,idicp)
!!$         
!!$      enddo
!!$      phi_out(jdicp,idicp) = temp1
!!$   enddo
!!$   
!!$enddo
!!$

!!$!=================================================================================================
!!$! the laser interaction term
!=================================================================================================
if(laser%tdornot) then

 do idicp=1,fedvr3d%nb_angle
  do jdicp =1,fedvr3d%nb_r    
    temp1 = zzero
    do ldicp =1,fedvr3d%nb_angle

     temp1 = temp1 + laser%ez_t*zmat_3d(idicp,ldicp)*fedvrx_global(jdicp)*phi_out((ldicp-1)*fedvr3d%nb_r+jdicp)
    enddo

  phi_out((idicp-1)*fedvr3d%nb_r+jdicp) = phi_out((idicp-1)*fedvr3d%nb_r+jdicp) + temp1


 enddo
enddo

endif
 
return
end subroutine act3
















