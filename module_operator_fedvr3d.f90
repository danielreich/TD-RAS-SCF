!===============================================================================
! calc. the matrix element of kinetic energy and nuclei attractive
! potential in fedvr_3d, spherical coordinate 
!===============================================================================
 module operator_3d
  use global
  use fedvr3d_basis_set
  use operator_radial
  implicit none 
  ! the matrix elements in the fedvr_3d of kinetic energy operator
  real(kind=k1),allocatable,save :: tmat_3d(:,:) 
  real(kind=k1),allocatable,save :: tmat_inv(:,:,:)
  ! the matrxi elements in the fedvr_3d of x,y,z
  real(kind=k1),allocatable,save :: xmat_3d(:,:),ymat_3d(:,:),zmat_3d(:,:)
  complex(kind=k2), allocatable,save :: px_element(:),py_element(:),pz_element(:) 
  !! index_pz(#1,1:2). The matrix element #1 = <index_pz(#1,1)|pz|index_pz(#1,2)>
  integer, allocatable,save :: index_px_element(:,:),index_py_element(:,:),index_pz_element(:,:)
  real(kind=k1),allocatable,save :: h_stiffness(:,:), vabsorb_pot_fedvr3d(:)  
  real(kind=k1),allocatable,save :: hl_stiffness(:,:,:)
  real(kind=k1),allocatable, save :: w_band(:), z_band(:,:),z_band_temp(:,:)
 
  contains

   subroutine calc_tmat_3d
    implicit none
    integer :: i_angle,l_in,e,i,idicp
    allocate(tmat_3d(fedvr3d%nb_r,fedvr3d%nb_angle))

!
! calc. the angular part of kinetic operator 
!	
    do i_angle =1,fedvr3d%nb_angle
      l_in=lm_l(i_angle)
      do idicp =1,fedvr3d%nb_r
  
      e = which_element(idicp)
      i = which_basis(idicp)

      tmat_3d(idicp,i_angle) = 1.0_k1*l_in*(l_in+1)/(2.0_k1*fedvr_x(e,i )**2) 
     enddo 
   enddo 
   return
   end subroutine calc_tmat_3d
!calc. the dipole moment matrix Z in spherical coordinate
!================================================================================================
!
! <n1,l1,m1  |Z| n2,l2,m2> = <n1,l1,m1  | r*cos(theta) | n2,l2,m2>  
!  
!
! because: cos(theta)  = Y  (theta, fphi) 
!                         10   
!
!       X  (r)                                                    X  (r) 
!        n1                                                        n2 
! = < ~~~~~~~~~~ *  Y   (theta, fphi)  | r * Y   (theta, fphi)| ~~~~~~~~~~~ *  Y  (theta, fphi) >  r^2 * sin(theta)* d(theta)*d(fphi)
!         r          l1m1                      10                    r          l2m2
!
!         
! = delta(n1,n2) * r(n1) * integrate { Y     * Y   *  Y     * sin(theta)}
!                                       l1m1    10     l2m2
!====================================================================================================

!=========================================================================================================================================
!
! <n1,l1,m1  |x| n2,l2,m2> = <n1,l1,m1  | r*sin(theta)*cos(fphi) | n2,l2,m2>  
!  
!                                     /~~~~~
! because: sin(theta)*cos(fphi)  =   / 2*pi    [   Y  (theta, fphi)  -  Y  (theta, fphi) ]  
!                                  |/ ~~~~~~        1-1                  11 
!                                       3 
!                           /~~~~ 
! = delta(n1,n2) * r(n1) * / 2*pi    [ integrate { Y     * Y   *  Y      * sin(theta)} - integrate { Y     * Y   *  Y      * sin(theta)}]
!                        \/ ~~~~~                   l1m1    1-1     l2m2                              l1m1    11     l2m2  
!                             3
!==========================================================================================================================================

 subroutine calc_xmat_3d
  implicit none
  integer :: l1,m1,l2,m2
  integer :: i_angle,j_angle
  real(kind=k1) :: getgaunt
  allocate(xmat_3d(fedvr3d%nb_angle,fedvr3d%nb_angle))

  do i_angle =1,fedvr3d%nb_angle
    l1 = lm_l(i_angle)
    m1 = lm_m(i_angle)  

    do j_angle =1,fedvr3d%nb_angle
       l2 = lm_l(j_angle)
       m2 = lm_m(j_angle)
     
       xmat_3d(i_angle,j_angle) = dsqrt(2.0d0*pi/3.0d0)*(getgaunt(l1,1,l2,m1,-1,m2) - getgaunt(l1,1,l2,m1,1,m2) )

    enddo
   enddo
  return
 end subroutine calc_xmat_3d

!=========================================================================================================================================
!
! <n1,l1,m1  |y| n2,l2,m2> = <n1,l1,m1  | r*sin(theta)*sin(fphi) | n2,l2,m2>  
!  
!                                     /~~~~~
! because: sin(theta)*cos(fphi)  =   / 2*pi  i  [   Y  (theta, fphi)  +  Y  (theta, fphi) ]  
!                                  |/ ~~~~~~         1-1                  11 
!                                       3 
!                           /~~~~ 
! = delta(n1,n2) * r(n1) * / 2*pi  i  [ integrate { Y     * Y   *  Y      * sin(theta)} + integrate { Y     * Y   *  Y      * sin(theta)}]
!                        \/ ~~~~~                    l1m1    1-1     l2m2                              l1m1    11     l2m2  
!                             3
!==========================================================================================================================================
 subroutine calc_ymat_3d
  implicit none
  integer :: l1,m1,l2,m2
  integer :: i_angle,j_angle
  real(kind=k1) :: getgaunt

  allocate(ymat_3d(1:fedvr3d%nb_angle,1:fedvr3d%nb_angle))
!
! can be improved here !!!!
!
  do i_angle =1,fedvr3d%nb_angle
    l1 = lm_l(i_angle)
    m1 = lm_m(i_angle)  

    do j_angle =1,fedvr3d%nb_angle
       l2 = lm_l(j_angle)
       m2 = lm_m(j_angle)
    
       ymat_3d(i_angle,j_angle) = sqrt(2.0d0*pi/3.0d0)* (getgaunt(l1,1,l2,m1,-1,m2) + getgaunt(l1,1,l2,m1,1,m2) )
     enddo
   enddo

  return
 end subroutine calc_ymat_3d


 subroutine calc_zmat_3d
  implicit none
  integer :: l1,m1,l2,m2
  integer :: i_angle,j_angle
  real(kind=k1) :: getgaunt

  allocate(zmat_3d(fedvr3d%nb_angle,fedvr3d%nb_angle))
!
! can be improved here !!!!
!

  do i_angle =1,fedvr3d%nb_angle
    l1 = lm_l(i_angle)
    m1 = lm_m(i_angle)  

    do j_angle =1,fedvr3d%nb_angle
       l2 = lm_l(j_angle)
       m2 = lm_m(j_angle)
     
       zmat_3d(i_angle,j_angle) = dsqrt(4.0d0*pi/3.0d0)*getgaunt(l1,1,l2,m1,0,m2)  

    enddo
  enddo

  return
 end subroutine calc_zmat_3d

!!==========================================================================
!! CALCULATION OF P_Z. Done by Juan Omiste
!!
!! In general it is complex.
!! It is necessary to calculate the interaction with the laser in the
!! velocity gauge.
!!==========================================================================

Subroutine pz_fedvr(pz,index_pz)
implicit none

!!INPUT
!! The input needed is defined as global variables.

!!OUTPUT
!! pz(#1) is the matrix element #1 of the linear momentum in z.
complex(kind=k2), allocatable,intent(out):: pz(:) 
!! index_pz(#1,1:2). The matrix element #1 = <index_pz(#1,1)|pz|index_pz(#1,2)>
integer, allocatable, intent(out):: index_pz(:,:)

!!AUXILIAR VARIABLES
integer :: i, j, ij,ji, counter
real(kind=k1) :: t_aux

!! First, we construct the index.

!! We calculate the elements which are not zero

counter=0

Do i=1,size(index_kinetic_basis(:,1)) !! Run in the FE-DVR functions (given by the second argument 1-> bra, and 2-> ket).
   If (index_kinetic_basis(i,1).eq.index_kinetic_basis(i,2)) cycle !! If r is the same, pz=0

   Do ij=1,fedvr3d%nb_angle
      Do ji=1,fedvr3d%nb_angle

         If (abs(lm_l(ij)-lm_l(ji)).eq.1.and.lm_m(ij).eq.lm_m(ji)) then 
!! First condition |Delta L|=+/- 1
!! Second condition |Delta M|=0
            counter=counter+1
         End If

      End Do
   End Do

End Do

!! Allocate the memory for the index_pz and the pz

allocate(index_pz(1:counter,1:2)) 
allocate(pz(1:counter))

!! Initialize

index_pz=0
pz=zzero

!! Store the value of the index

counter=0

Do i=1,size(index_kinetic_basis(:,1)) !! Run in the FE-DVR functions (given by the second argument 1-> bra, and 2-> ket).
   If (index_kinetic_basis(i,1).eq.index_kinetic_basis(i,2)) cycle !! If r is the same, pz=0
   Do ij=1,fedvr3d%nb_angle 
      Do ji=1,fedvr3d%nb_angle

         If (abs(lm_l(ij)-lm_l(ji)).eq.1.and.lm_m(ij).eq.lm_m(ji)) then 
!! First condition Delta L=+/- 1
!! Second condition Delta M=0
            counter=counter+1

!!INDEX storing

            index_pz(counter,1)=(ij-1)*fedvr3d%nb_r+index_kinetic_basis(i,1)
            index_pz(counter,2)=(ji-1)*fedvr3d%nb_r+index_kinetic_basis(i,2)

!!Pz

            pz(counter)=-ci*cmplx(tmat_radial(i)*zmat_3d(ij,ji)*(fedvrx_global(index_kinetic_basis(i,1))-fedvrx_global(index_kinetic_basis(i,2))))            
            
         End If

      End Do
   End Do
End Do

End Subroutine pz_fedvr

!!==========================================================================
!! CALCULATION OF P_X. Done by Juan Omiste
!!
!! In general it is complex.
!! It is necessary to calculate the interaction with the laser in the
!! velocity gauge.
!!==========================================================================

Subroutine px_fedvr(px,index_px)
implicit none

!!INPUT
!! The input needed is defined as global variables.

!!OUTPUT
!! px(#1) is the matrix element #1 of the linear momentum in z.
complex(kind=k2), allocatable,intent(out):: px(:) 
!! index_px(#1,1:2). The matrix element #1 = <index_px(#1,1)|px|index_px(#1,2)>
integer, allocatable, intent(out):: index_px(:,:)

!!AUXILIAR VARIABLES
integer :: i, j, ij,ji, counter
real(kind=k1) :: t_aux

!! First, we construct the index.

!! We calculate the elements which are not zero

counter=0

Do i=1,size(index_kinetic_basis(:,1)) !! Run in the FE-DVR functions (given by the second argument 1-> bra, and 2-> ket).
   If (index_kinetic_basis(i,1).eq.index_kinetic_basis(i,2)) cycle !! If r is the same, px=0
   Do ij=1,fedvr3d%nb_angle
      Do ji=1,fedvr3d%nb_angle

         If (abs(lm_l(ij)-lm_l(ji)).eq.1.and.abs(lm_m(ij)-lm_m(ji)).eq.1) then 
!! First condition |Delta L|=+/- 1
!! Second condition |Delta M|=+/- 1
            counter=counter+1
         End If

      End Do
   End Do
End Do

!! Allocate the memory for the index_px and the px

allocate(index_px(1:counter,1:2)) 
allocate(px(1:counter))

!! Initialize

index_px=0
px=zzero

!! Store the value of the index

counter=0

Do i=1,size(index_kinetic_basis(:,1)) !! Run in the FE-DVR functions (given by the second argument 1-> bra, and 2-> ket).
   If (index_kinetic_basis(i,1).eq.index_kinetic_basis(i,2)) cycle !! If r is the same, px=0
   Do ij=1,fedvr3d%nb_angle
      Do ji=1,fedvr3d%nb_angle

         If (abs(lm_l(ij)-lm_l(ji)).eq.1.and.abs(lm_m(ij)-lm_m(ji)).eq.1) then 
!! First condition Delta L=+/- 1
!! Second condition Delta M=+/- 1
            counter=counter+1

!!INDEX storing

            index_px(counter,1)=((ij-1)*fedvr3d%nb_r+index_kinetic_basis(i,1))
            index_px(counter,2)=((ji-1)*fedvr3d%nb_r+index_kinetic_basis(i,2))

!!Px

            px(counter)=-ci*cmplx(tmat_radial(i)*xmat_3d(ij,ji)*(fedvrx_global(index_kinetic_basis(i,1))-fedvrx_global(index_kinetic_basis(i,2))))

            
         End If

      End Do
   End Do
End Do


End Subroutine px_fedvr

!!==========================================================================
!! CALCULATION OF P_Y. Done by Juan Omiste
!!
!! In general it is complex.
!! It is necessary to calculate the interaction with the laser in the
!! velocity gauge.
!!==========================================================================

Subroutine py_fedvr(py,index_py)
implicit none

!!INPUT
!! The input needed is defined as global variables.

!!OUTPUT
!! py(#1) is the matrix element #1 of the linear momentum in y.
complex(kind=k2), allocatable,intent(out):: py(:) 
!! index_py(#1,1:2). The matrix element #1 = <index_py(#1,1)|py|index_py(#1,2)>
integer, allocatable, intent(out):: index_py(:,:)

!!AUXILIAR VARIABLES
integer :: i, j, ij,ji, counter
real(kind=k1) :: t_aux

write(*,*) 'py_fedvr in module_operator_fedvr3d.f90'
write(*,*) 'use that the kinetic energy index is not a full index'
return

!! First, we construct the index.

!! We calculate the elements which are not zero

counter=0

Do i=1,size(index_kinetic_basis(:,1)) !! Run in the FE-DVR functions (given by the second argument 1-> bra, and 2-> ket).
   Do ij=1,fedvr3d%nb_angle
      Do ji=1,fedvr3d%nb_angle

         If (abs(lm_l(ij)-lm_l(ji)).eq.1.and.abs(lm_m(ij)-lm_m(ji)).eq.1) then 
!! First condition Delta L=+/- 1
!! Second condition Delta M=+/- 1
            counter=counter+1
         End If

      End Do
   End Do
End Do

!! Allocate the memory for the index_py and the py

allocate(index_py(1:counter,1:2)) 
allocate(py(1:counter))

!! Initialize

index_py=0
py=zzero

!! Store the value of the index

counter=0

Do i=1,size(index_kinetic_basis(:,1)) !! Run in the FE-DVR functions (given by the second argument 1-> bra, and 2-> ket).
   Do ij=1,fedvr3d%nb_angle
      Do ji=1,fedvr3d%nb_angle

         If (abs(lm_l(ij)-lm_l(ji)).eq.1.and.abs(lm_m(ij)-lm_m(ji)).eq.1) then 
!! First condition Delta L=+/- 1
!! Second condition Delta M=+/- 1
            counter=counter+1

!!INDEX storing

            index_py(counter,1)=((ij-1)*fedvr3d%nb_r+index_kinetic_basis(i,1))
            index_py(counter,2)=((ji-1)*fedvr3d%nb_r+index_kinetic_basis(i,2))

!!Py

!!Calculation of <a|Ty|b>
!!Contribution of the contribution of T

            If (index_kinetic_basis(i,1).eq.index_kinetic_basis(i,2)) then
               t_aux=tmat_radial(i)+tmat_3d(index_kinetic_basis(i,1),ij)
            Else
               t_aux=tmat_radial(i)
            End If

!!Including z

            py(counter)=t_aux*ymat_3d(ij,ji)*fedvrx_global(index_kinetic_basis(i,2))*cmplx(0.0d0, 1.0d0)

!!Calculation of <a|yT|b>
!!Contribution of the contribution of T

            If (index_kinetic_basis(i,1).eq.index_kinetic_basis(i,2)) then
               t_aux=tmat_radial(i)+tmat_3d(index_kinetic_basis(i,1),ji)
            Else
               t_aux=tmat_radial(i)
            End If

!!Including x and added to the one done before

            py(counter)=py(counter)-t_aux*ymat_3d(ij,ji)*fedvrx_global(index_kinetic_basis(i,1))*cmplx(0.0d0, 1.0d0)    

         End If

      End Do
   End Do
End Do

py=py*ci

End Subroutine py_fedvr

!
! calc. the tkk -1
!
!
   subroutine calc_tmat_inv
    implicit none
	real(kind=k1),allocatable :: rtmp(:,:),work(:)
    integer :: nb_r,ibasis,jbasis
    integer,allocatable :: ipiv(:)
    integer :: idicp,jdicp,info,l,itemp,kdicp

        nb_r = fedvr3d%nb_r

       allocate(tmat_inv(nb_r,nb_r,0:2*fedvr3d%l_max ))
       allocate(rtmp(nb_r,nb_r ),work(nb_r),ipiv(nb_r))
       tmat_inv = 0.0_k1
    
       do l=0,2*fedvr3d%l_max

       rtmp = 0.0_k1

!       do idicp=1,nb_r
!        do jdicp =1,nb_r

!          itemp = index_kinetic_operator(idicp,jdicp)
!          if(itemp /=0 ) then
!              rtmp(idicp,jdicp)  = 2.0_k1*tmat_fedvr3d_radial(itemp)
!          else
!              rtmp(idicp,jdicp)  = 0.0_k1 
!          endif
!        if(idicp == jdicp) then
!          rtmp(idicp,jdicp) =rtmp(idicp,jdicp) + 1.0_k1*l*(l+1)/(fedvrx_global(idicp)**2) 
!        endif
!       enddo
!      enddo
        rtmp = 0.0d0
        do idicp =1, n_total_kinetic
          ibasis = index_kinetic_basis(idicp,1)
          jbasis = index_kinetic_basis(idicp,2)

           if(ibasis/=jbasis) then 
             rtmp(ibasis,jbasis) = 2.0d0*tmat_radial(idicp)
             rtmp(jbasis,ibasis) = 2.0d0*tmat_radial(idicp) 
           else
             rtmp(ibasis,jbasis) = 2.0d0*tmat_radial(idicp) + 1.0_k1*l*(l+1)/(fedvrx_global(ibasis)**2) 
           endif
         enddo  


         ! computing matrix inverse
         ! first step: calculate LU factorization
         call dgetrf(Nb_r, Nb_r, rtmp, Nb_r, ipiv, info)
         ! second step: invert LU matrix
         call dgetri(Nb_r, rtmp, Nb_r, ipiv, work, Nb_r,info )

         do idicp =1,fedvr3d%nb_r
           do jdicp =1,fedvr3d%nb_r
            tmat_inv(idicp,jdicp,l) = rtmp(idicp,jdicp)
           enddo
          enddo
       enddo
    return
   end subroutine calc_tmat_inv

!
! absorb potential in radial part
!

  subroutine calc_absorb3d()
   implicit none
   integer(4) :: idicp
   real(kind=k1) :: fun_absorb_rightside,xgrid
   
   allocate(vabsorb_pot_fedvr3d(fedvr3d%nb_r))
   vabsorb_pot_fedvr3d = 0.0d0
  
   do idicp =1, fedvr3d%nb_r
    xgrid = fedvrx_global(idicp)  
    
    vabsorb_pot_fedvr3d(idicp) = fun_absorb_rightside(xgrid)
 
  enddo
  return
 end subroutine calc_absorb3d

!
! in the spherical coordiate ,calc. the matrix element of kinetic energy operator(angluar part)
!
 subroutine drive_operator_3d
  implicit none
  integer :: idicp,jdicp
!calc. the kinetic energy element matrix in fedvr3d
   call  calc_tmat_3d
 
   call calc_tmat_inv

!!$  to calculate the dipole moment

   call  calc_xmat_3d
   call  calc_ymat_3d
   call  calc_zmat_3d


  if(laser%tdornot)  then
!calc. the elements in fedvr3d of x,y,z dipole matrix

     if (laser%gauge.eq.'v') then
        
        !We need them to calculate the dipole and the interaction with the
        !laser in the length gauge.
        
        !Calculation in the fedvr3d of p_x, p_y, p_z for the velocity gauge
        
        write(*,*) ''
        write(*,*) 'CALCULATION OF THE MOMENTUM OPERATOR TO USE THE VELOCITY GAUGE'
        write(*,*) ''
        
        call px_fedvr(px_element,index_px_element)
        call py_fedvr(py_element,index_py_element)
        call pz_fedvr(pz_element,index_pz_element)
        
        write(*,*) ''
        write(*,*) 'DONE'
        write(*,*) ''
        
     End if
!Calculation of the absorbing potential

     write(*,*) ''
     write(*,*) 'CALCULATION OF THE ABSORBING POTENTIAL'
     write(*,*) ''

   call  calc_absorb3d()

     write(*,*) ''
     write(*,*) 'DONE'
     write(*,*) ''

 endif


  return
 end subroutine drive_operator_3d


!! The subroutine act_index applies an operator to the orbitals, where only the nonzero elements are stored and the upper diagonal are stored. The elements are complex.

 Subroutine act_index(phi,op,index_op,op_phi)
   implicit none
   !! INPUT
   !! phi(#1,#2) is the coefficient corresponding to the #1 spatial basis function and #2 is the orbital.
   complex (kind=k2), allocatable, intent(in):: phi(:,:) 
   !! op(#1) is the matrix element #1 of the linear momentum in z.
   complex(kind=k2), allocatable,intent(in):: op(:) 
   !! index_op(#1,1:2). The matrix element #1 = <index_op(#1,1)|op|index_op(#1,2)>
   integer, allocatable, intent(in):: index_op(:,:)
   !! OUTPUT
   !! op_phi(#1,#2) is the coefficient corresponding to the #1 spatial basis function and #2 is the orbital after applying the operator 'op'.
   complex (kind=k2), allocatable,intent(out):: op_phi(:,:) 
   !! AUXILIAR
   integer :: j, k  !! loops
   integer :: ntotal, norbital !! size of the arrays
   integer :: nelements !! number of elements of the operator
   integer :: nbra, nket !! index in bra and ket
   
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


   !! allocate and initialize the output
   
   ntotal = size(phi(:,1)) !! number of spatial basis functions
   norbital = size(phi(1,:)) !! number of orbitals

   if (.not.allocated(op_phi)) allocate(op_phi(1:ntotal, 1:norbital))

   op_phi = zzero

   !! # elements of the operator

   nelements = size(op(:)) !! number of elements

   !!PARALLELIZE HERE
   !$OMP PARALLEL PRIVATE(j,nbra,nket)
   !$OMP DO SCHEDULE(dynamic) 

   Do k=1,norbital
      Do j=1,nelements
         nbra=index_op(j,1) !! function in the bra
         nket=index_op(j,2) !! function in the ket
         op_phi(nbra,k)=op_phi(nbra,k) + op(j)*phi(nket,k) !! apply the operator
         If (nbra.ne.nket) op_phi(nket,k)=op_phi(nket,k)+conjg(op(j))*phi(nbra,k) !! if (nbra,nket) is not in the diagonal
      End Do
   End Do
   
   !$OMP END DO
   !$OMP END PARALLEL

   return
 End Subroutine act_index


 end module operator_3d

