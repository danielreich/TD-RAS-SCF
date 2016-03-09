!
! calc. the two-e integrals in basis set
!
module twoe_basis_set
  use global
  use auxiliary
  use fedvr3d_basis_set
  use operator_3d
  implicit none
  integer, allocatable,save :: index_two(:,:,:)
  integer, allocatable,save ::  num_two(:,:)
  real(kind=k1),allocatable,save :: two_storage(:,:,:)
  integer,allocatable,save :: index_two_storage(:,:,:,:) 
  !Included by Juan
  Double precision, allocatable, save :: twoe_radial_store(:,:,:) !! r1,r2 and l
  Integer, allocatable, save :: twoe_radial_r1(:),twoe_radial_r2(:) !! r1 and r2 corresponding to a certain index
  Double precision, allocatable, save :: twoe_angle_store(:,:,:,:,:)
  !End of included by Juan

!  type(twoe_angle_index_type) :: twoe_angle_index
  integer,allocatable :: twoe_ang1(:),twoe_ang2(:),twoe_ang3(:), twoe_ang4(:)
  Double precision,allocatable :: twoe_integral_ll(:,:)
  Double precision, allocatable :: twoe_angle_index_integral(:)

  contains

subroutine twoe_radial(rl)
implicit none

!! INPUT
!!
!! The input is stored in module_global. That is:
!!
!! Maximum orbital angular momentum: fedvr3d%l_max
!! Number of nodes in the radial grid: 1:fedvr3d%nb_r

!! OUTPUT
!!
Real(kind=k1), allocatable, intent(out) :: rl(:,:,:)
!! rl is the radial contribution to the two electron integrals.
!! Arguments
!! #1 is the radial part of r_1
!! #2 is the radial part of r_2. The integration to obtain the 
!!    mean field operator is over r_2.
!! #3 is L of the operator. L goes from 0 to 2*fevr3d%l_max
!! 

!! AUXILIAR VARIABLES

integer :: r1,r2, l !! indeces for the loops.

!! First we allocate and initialize rl

allocate(rl(1:fedvr3d%nb_r,1:fedvr3d%nb_r,0:2*fedvr3d%l_max))
rl=zero

Do l=0,2*fedvr3d%l_max !! Run in l
   Do r2=1,fedvr3d%nb_r !! Run in r2 radial grid
      Do r1=1,fedvr3d%nb_r !! Run in r1 radial grid
     
               rl(r1,r2,l)=fedvr3d_twoe_radial(r1,r2,l)

!! We include the factor 4\pi/(2*l+1)
               
               rl(r1,r2,l)=rl(r1,r2,l)*4.0d0*acos(-1.0d0)/dble(2*l+1)

       End Do
    End Do
 End Do

!! This subroutine can be improved by using that each block with the same l is real and symmetric

return
End subroutine twoe_radial

!! This subroutines prepares the two electron integrals to calculate the mean field operator.
!! If fedvr3d%store is equal to 0, the contributions to the radial part and to the angular part are stored separately.
!! If fedvr3d%store is equal to 1, the contributions are stored together.

  subroutine satis_two
   implicit none
   integer(4) :: idicp,jdicp,kdicp,ldicp,isum,n1,l1,m1,n2,l2,m2,n3,l3,m3,n4,l4,m4
   integer(4) :: isum_total
   double precision :: twotwo
   integer(4) :: i_global, j_global,k_global,l_global
   integer(4) :: kk, kkp,ll!! added by Juan
   integer(4) :: ang1,ang2,ang3,ang4,angaux!! added by Juan

   select case (fedvr3d%store)
   case (0)

!! Radial part
      
      allocate(twoe_radial_store(1:fedvr3d%nb_r,1:fedvr3d%nb_r,0:2*fedvr3d%l_max))
      twoe_radial_store=zero
      

      Do ll=0,2*fedvr3d%l_max
         Do kk=1, fedvr3d%nb_r
            twoe_radial_store(kk,kk,ll)=fedvr3d_twoe_radial(kk,kk,ll)
            Do kkp=kk+1, fedvr3d%nb_r
               twoe_radial_store(kk,kkp,ll)=fedvr3d_twoe_radial(kk,kkp,ll)
               twoe_radial_store(kkp,kk,ll)=twoe_radial_store(kk,kkp,ll)
            End Do
         End Do
      End Do
 
!! Angular part
      
      allocate(twoe_angle_store(1:fedvr3d%nb_angle,1:fedvr3d%nb_angle,1:fedvr3d%nb_angle,1:fedvr3d%nb_angle,0:2*fedvr3d%l_max))

      print*, 'Size of the radial part of the two electron integral', size(twoe_radial_store)
      print*, 'Size of the angular part of the two electron integral',size(twoe_angle_store)

      Do ang1=1,fedvr3d%nb_angle
         Do ang2=1,fedvr3d%nb_angle
            Do ang3=1,fedvr3d%nb_angle
               Do ang4=1,fedvr3d%nb_angle                     
                  Do ll=0,2*fedvr3d%l_max
                                          
                     l1=  lm_l(ang1)
                     m1=  lm_m(ang1)
                     l2=  lm_l(ang2)
                     m2=  lm_m(ang2)
                     l3=  lm_l(ang3)
                     m3=  lm_m(ang3)
                     l4=  lm_l(ang4)
                     m4=  lm_m(ang4)
                     
                     twoe_angle_store(ang1,ang2,ang3,ang4,ll)=fedvr3dbase_angpart(ll,l1,m1,l2,m2,l3,m3,l4,m4)

                  End Do
               End Do
            End Do
         End Do
      End Do

!!$      End Do
      
      return !! no store cww
   case (1) 
      continue !! stores cww
   End select
   
   allocate(index_two(2,fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_angle)) 
   allocate(num_two( fedvr3d%nb_r*fedvr3d%nb_angle, fedvr3d%nb_angle ))
   allocate(index_two_storage(4,fedvr3d%nb_angle**2*fedvr3d%nb_r,fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_angle))
   allocate(two_storage(fedvr3d%nb_angle**2*fedvr3d%nb_r,fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_angle ))
!
! initilization the index_two
!
   do idicp= 1, fedvr3d%nb_r
     do jdicp =1, fedvr3d%nb_angle
      index_two(1,idicp,jdicp) = idicp 
      index_two(2,idicp,jdicp) = (jdicp -1 ) *fedvr3d%nb_r + idicp 
     enddo
   enddo


  do idicp = fedvr3d%nb_r + 1 , fedvr3d%nb_r*fedvr3d%nb_angle

      index_two(1,idicp,:) =  index_two(1,idicp-fedvr3d%nb_r,:)  +  fedvr3d%nb_r 
      index_two(2,idicp,:) =  index_two(2,idicp-fedvr3d%nb_r,:) 

  enddo
!
!
! < n1,l1,m1(1) n2,l2,m2(1)  n3l3m3(2) ,n4l4m4(2) >
!
   
   isum_total = 0 
   do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
    do jdicp = 1,fedvr3d%nb_angle

     i_global = index_two(1,idicp,jdicp)
     j_global = index_two(2,idicp,jdicp)

       

    n1 =  global_to_local(i_global,1)
    l1 =  global_to_local(i_global,2)
    m1 =  global_to_local(i_global,3)
   
   
    n2 =  global_to_local(j_global,1)
    l2 =  global_to_local(j_global,2)
    m2 =  global_to_local(j_global,3)

   isum = 0 
   do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
      do ldicp =1,fedvr3d%nb_angle
              
      k_global = index_two(1,kdicp,ldicp)
      l_global = index_two(2,kdicp,ldicp)

      n3 =  global_to_local(k_global,1)
      l3 =  global_to_local(k_global,2)
      m3 =  global_to_local(k_global,3)

      n4 =  global_to_local(l_global,1)
      l4 =  global_to_local(l_global,2)
      m4 =  global_to_local(l_global,3)

     if((n1==n2).and.(n3==n4) .and. ((m1-m2)==(m4-m3))) then   
      twotwo  =  fedvr3d_twoe(n1,l1,m1,l2,m2,n3,l3,m3,l4,m4) 
      if(abs(twotwo)>0.000000001d0) then

        isum = isum  + 1 
        isum_total = isum_total + 1

        if(isum> fedvr3d%nb_angle**2*fedvr3d%nb_r) then
           write(*,*) 'increase the dim', isum, fedvr3d%nb_r*fedvr3d%nb_angle
           stop ' satis_two_integrals.f90'
        endif
       
       index_two_storage(1,isum,idicp,jdicp) = k_global
       index_two_storage(2,isum,idicp,jdicp) = l_global

       two_storage(isum,idicp,jdicp) = twotwo
      endif 
   endif
enddo
enddo
  num_two(idicp,jdicp) = isum
 enddo
enddo

  write(*,*) 'Only ', isum_total, ' two electron integrals is useding '
  
  return
end subroutine satis_two

!! This subroutines prepares the two electron integrals to calculate the mean field operator.
!! If fedvr3d%store is equal to 0, the contributions to the radial part and to the angular part are stored separately. satis_two_opt is optimized to calculate only the non-zero contributions.
!! If fedvr3d%store is equal to 1, the contributions are stored together.

  subroutine satis_two_opt
   implicit none
   integer(4) :: idicp,jdicp,kdicp,ldicp,isum,n1,l1,m1,n2,l2,m2,n3,l3,m3,n4,l4,m4
   integer(4) :: isum_total
   double precision :: twotwo
   integer(4) :: i_global, j_global,k_global,l_global
   integer(4) :: kk, kkp, kk_kkp,ll!! added by Juan
   integer(4) :: ang1,ang2,ang3,ang4,angaux!! added by Juan
   integer(4) :: counter
   

   select case (fedvr3d%store)
   case (0)

!! Radial part
!! Storing the index that we use in the optimize region

      allocate(twoe_radial_r1(1:fedvr3d%nb_r*(fedvr3d%nb_r+1)/2))
      allocate(twoe_radial_r2(1:fedvr3d%nb_r*(fedvr3d%nb_r+1)/2))
      
      twoe_radial_r1=0
      twoe_radial_r2=0
      
      counter=0

      Do kk=1,fedvr3d%nb_r
         Do kkp=kk,fedvr3d%nb_r
            counter=counter+1
            twoe_radial_r1(counter)=kk
            twoe_radial_r2(counter)=kkp
         End Do
      End Do

      counter=0

!! Storing the integrals
      
      allocate(twoe_radial_store(1:fedvr3d%nb_r,1:fedvr3d%nb_r,0:2*fedvr3d%l_max))
      twoe_radial_store=zero      

      Do ll=0,2*fedvr3d%l_max
         Do kk=1, fedvr3d%nb_r
            twoe_radial_store(kk,kk,ll)=fedvr3d_twoe_radial(kk,kk,ll)
            Do kkp=kk+1, fedvr3d%nb_r
               twoe_radial_store(kk,kkp,ll)=fedvr3d_twoe_radial(kk,kkp,ll)
               twoe_radial_store(kkp,kk,ll)=twoe_radial_store(kk,kkp,ll)
            End Do
         End Do
      End Do

       print*, 'Size of the radial part of the two electron integral', size(twoe_radial_store)

!! Angular part

!!First, we check the index in the integral

      counter=0
      
      Do ang1=1,fedvr3d%nb_angle
         Do ang2=1,fedvr3d%nb_angle
            Do ang3=1,fedvr3d%nb_angle
               Do ang4=1,fedvr3d%nb_angle                     

                     l1=  lm_l(ang1)
                     m1=  lm_m(ang1)
                     l2=  lm_l(ang2)
                     m2=  lm_m(ang2)
                     l3=  lm_l(ang3)
                     m3=  lm_m(ang3)
                     l4=  lm_l(ang4)
                     m4=  lm_m(ang4)

            If ((m1-m2).ne.(m4-m3)) cycle !! These integrals are zero
            If (max(abs(l1-l2),abs(l3-l4)).gt.min((l1+l2),(l3+l4))) cycle !! These integrals are zero

               counter=counter+1                        
            
         End Do
      End Do
   End Do
End Do

      print*, 'Number of angular integrals of the two electron integral', counter

!! Second, we allocate the memory to store twoe_integral

      allocate(twoe_ang1(1:counter))
      allocate(twoe_ang2(1:counter))
      allocate(twoe_ang3(1:counter))
      allocate(twoe_ang4(1:counter))
      allocate(twoe_integral_ll(1:counter,0:2*fedvr3d%l_max))
      
!! Finally, we run again to store the index
      
      counter=0
      twoe_integral_ll=zzero

      Do ang1=1,fedvr3d%nb_angle
         Do ang2=1,fedvr3d%nb_angle
            Do ang3=1,fedvr3d%nb_angle
               Do ang4=1,fedvr3d%nb_angle 
                  
                  l1=  lm_l(ang1)
                  m1=  lm_m(ang1)
                  l2=  lm_l(ang2)
                  m2=  lm_m(ang2)
                  l3=  lm_l(ang3)
                  m3=  lm_m(ang3)
                  l4=  lm_l(ang4)
                  m4=  lm_m(ang4)
                  
            If ((m1-m2).ne.(m4-m3)) cycle !! This integrals are zero
            
            If (max(abs(l1-l2),abs(l3-l4)).gt.min((l1+l2),(l3+l4))) cycle !! These integrals are zero
            
               counter=counter+1

               twoe_ang1(counter)=ang1
               twoe_ang2(counter)=ang2
               twoe_ang3(counter)=ang3
               twoe_ang4(counter)=ang4                    

            Do ll=max(abs(l1-l2),abs(l3-l4)), min(l1+l2,l3+l4)              
               
               twoe_integral_ll(counter,ll)=&
                    fedvr3dbase_angpart(ll,l1,m1,l2,m2,l3,m3,l4,m4)
               
            End Do
         End Do
      End Do
   End Do
End Do

      return !! no store cww
   case (1) 
      !! stores cww

   allocate(index_two(2,fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_angle)) 
   allocate(num_two( fedvr3d%nb_r*fedvr3d%nb_angle, fedvr3d%nb_angle ))
   allocate(index_two_storage(4,fedvr3d%nb_angle**2*fedvr3d%nb_r,fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_angle))
   allocate(two_storage(fedvr3d%nb_angle**2*fedvr3d%nb_r,fedvr3d%nb_r*fedvr3d%nb_angle,fedvr3d%nb_angle ))
!
! initilization the index_two
!
   do idicp= 1, fedvr3d%nb_r
     do jdicp =1, fedvr3d%nb_angle
      index_two(1,idicp,jdicp) = idicp 
      index_two(2,idicp,jdicp) = (jdicp -1 ) *fedvr3d%nb_r + idicp 
     enddo
   enddo

  do idicp = fedvr3d%nb_r + 1 , fedvr3d%nb_r*fedvr3d%nb_angle

      index_two(1,idicp,:) =  index_two(1,idicp-fedvr3d%nb_r,:)  +  fedvr3d%nb_r 
      index_two(2,idicp,:) =  index_two(2,idicp-fedvr3d%nb_r,:) 

  enddo



!
!
! < n1,l1,m1(1) n2,l2,m2(1)  n3l3m3(2) ,n4l4m4(2) >
!
   
   isum_total = 0 
   do idicp =1, fedvr3d%nb_r*fedvr3d%nb_angle
    do jdicp = 1,fedvr3d%nb_angle

     i_global = index_two(1,idicp,jdicp)
     j_global = index_two(2,idicp,jdicp)

       

    n1 =  global_to_local(i_global,1)
    l1 =  global_to_local(i_global,2)
    m1 =  global_to_local(i_global,3)
   
   
    n2 =  global_to_local(j_global,1)
    l2 =  global_to_local(j_global,2)
    m2 =  global_to_local(j_global,3)

   isum = 0 
   do kdicp =1,fedvr3d%nb_r*fedvr3d%nb_angle
      do ldicp =1,fedvr3d%nb_angle
              
      k_global = index_two(1,kdicp,ldicp)
      l_global = index_two(2,kdicp,ldicp)

      n3 =  global_to_local(k_global,1)
      l3 =  global_to_local(k_global,2)
      m3 =  global_to_local(k_global,3)

      n4 =  global_to_local(l_global,1)
      l4 =  global_to_local(l_global,2)
      m4 =  global_to_local(l_global,3)

     if((n1==n2).and.(n3==n4) .and. ((m1-m2)==(m4-m3))) then   
      twotwo  =  fedvr3d_twoe(n1,l1,m1,l2,m2,n3,l3,m3,l4,m4) 
      if(abs(twotwo)>0.000000001d0) then

        isum = isum  + 1 
        isum_total = isum_total + 1

        if(isum> fedvr3d%nb_angle**2*fedvr3d%nb_r) then
           write(*,*) 'increase the dim', isum, fedvr3d%nb_r*fedvr3d%nb_angle
           stop ' satis_two_integrals.f90'
        endif
       
       index_two_storage(1,isum,idicp,jdicp) = k_global
       index_two_storage(2,isum,idicp,jdicp) = l_global

       two_storage(isum,idicp,jdicp) = twotwo
      endif 
   endif
enddo
enddo
  num_two(idicp,jdicp) = isum
 enddo
enddo

  write(*,*) 'Only ', isum_total, ' two electron integrals is useding '

case(2) !! stores only the radial part of the 

write(*,*)
write(*,*) 'CONSTRUCTION OF THE RADIAL PART'
write(*,*)

     call twoe_radial(twoe_radial_store)

  End select

  return
end subroutine satis_two_opt


!
! calc. the two-e integrals
!
!============================================================================================
! two electron part
!
! < f_k,l1,m1, f_k,l2,m2  ||f_kp,l3,m3, f_kp,l4,m4 >
!
!============================================================================================
double precision function fedvr3d_twoe(k,l1,m1,l2,m2,kp,l3,m3,l4,m4)
  implicit none
  integer(4),intent(in) :: k,l1,m1,l2,m2,kp,l3,m3,l4,m4 
  integer(4) :: lmin,lmax,l,idicp,jdicp,ldicp,kdicp
  real(kind=k1) :: rsum 

  if((m1 - m2) /= (m4 - m3)) then
     fedvr3d_twoe = 0.0_k1
     return
  endif

  lmin = max(abs(l1-l2), abs(l3-l4))
  
  lmax = min((l1+l2), (l3+l4))
  
  rsum = 0.0_k1
  do l= lmin,lmax
     select case (fedvr3d%store)
     case(0)
        if (abs(twoe_radial_store(k,kp,l)).lt.1d-15) cycle        
        rsum = rsum + fedvr3dbase_angpart(L,l1,m1,l2,m2,l3,m3,l4,m4)*&
             twoe_radial_store(k,kp,l)        
     case(1)
        rsum = rsum + fedvr3dbase_angpart(L,l1,m1,l2,m2,l3,m3,l4,m4)*&
             fedvr3d_twoe_radial(k,kp,l)
     end select
  enddo
  fedvr3d_twoe = rsum
  return
end function fedvr3d_twoe

!
!two electrons integral, radial part
!
  double precision function fedvr3d_twoe_radial(k,kp,l)
   implicit none
   integer(4) :: k,kp,l,e,i,f,j,nb_e,nb_f
   
   e = which_element(k);  i = which_basis(k)
   f = which_element(kp); j = which_basis(kp)
   nb_e = fedvr3d%fedvr_nb(e)
   nb_f = fedvr3d%fedvr_nb(f) 

  if(i<nb_e) then
  
     if(j<nb_f) then


      fedvr3d_twoe_radial = (2.0_k1*l + 1.0_k1) *Tmat_inv(k,kp,l)/ &
      (fedvr_x(e,i)*fedvr_x(f,j)*dsqrt( fedvr_w(e,i) * fedvr_w(f,j)))  +   &
      (fedvr_x(e,i)**l*fedvr_x(f,j)**l/(fedvr3d%r_end**(2.0_k1*l+1.0_k1)))

      elseif(j==nb_f .and. f/= fedvr3d%number_of_element ) then

      fedvr3d_twoe_radial = (2.0_k1*l + 1.0_k1) *Tmat_inv(k,kp,l)/ &
      (fedvr_x(e,i)*fedvr_x(f,j)*dsqrt( fedvr_w(e,i) * (fedvr_w(f,j) + fedvr_w(f+1,1)) ))  +   &
      (fedvr_x(e,i)**l*fedvr_x(f,j)**l/(fedvr3d%r_end**(2.0_k1*l+1.0_k1)))


     else
       fedvr3d_twoe_radial = 0.0_k1
     endif
  elseif((i==nb_e) .and. (e/= fedvr3d%number_of_element)) then
  
    if(j<nb_f) then
      fedvr3d_twoe_radial = (2.0_k1*l + 1.0_k1) *Tmat_inv(k,kp,l)/ &
      (fedvr_x(e,i)*fedvr_x(f,j)*dsqrt( (fedvr_w(e,i)+fedvr_w(e+1,1)) * fedvr_w(f,j)))  +   &
      (fedvr_x(e,i)**l*fedvr_x(f,j)**l/(fedvr3d%r_end**(2.0_k1*l+1.0_k1)))
     elseif(j==nb_f.and. f/=fedvr3d%number_of_element) then

      fedvr3d_twoe_radial = (2.0_k1*l + 1.0_k1) *Tmat_inv(k,kp,l)/ &
      (fedvr_x(e,i)*fedvr_x(f,j)*dsqrt( (fedvr_w(e,i)+fedvr_w(e+1,1)) * (fedvr_w(f,j)+fedvr_w(f+1,1)  )))  +   &
      (fedvr_x(e,i)**l*fedvr_x(f,j)**l/(fedvr3d%r_end**(2.0_k1*l+1.0_k1)))

     else
      fedvr3d_twoe_radial = 0.0_k1
     endif 
  else 
    fedvr3d_twoe_radial = 0.0_k1
  endif
end  function fedvr3d_twoe_radial

!
! angular part
!
    double precision function fedvr3dbase_angpart(L,l1,m1,l2,m2,l3,m3,l4,m4)
      implicit none
      integer :: L
      integer :: l1,l2,l3,l4
      integer :: m1,m2,m3,m4
	  real(kind=k1) :: pi, getgaunt

	  pi = acos(-1.0d0)!! Juan: I included this !3.1415926d0 !Wenliang

      fedvr3dbase_angpart=(-1.0d0)**(m4-m3)*4.0d0*pi/(2.0d0*L+1.0d0)*&
           getgaunt(l1,l2,l,m1,m2,m1-m2)*getgaunt(l3,l,l4,m3,m3-m4,m4)
       end function fedvr3dbase_angpart

end module twoe_basis_set

subroutine drive_twoe_basis_set
 use twoe_basis_set
 implicit none

 If (fedvr3d%store.eq.1.or.fedvr3d%store.eq.2) then
!! Storing the integrals in the global basis.
!!$   call satis_two
    call satis_two_opt !! Juan: added on 8th December
 Else 
     
    continue
 End If

 return
end subroutine drive_twoe_basis_set
