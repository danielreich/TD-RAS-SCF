
!======================================================================================
! calculate the matrix elment of kinetic and nuclei attractive potential operator
! in radial part
!======================================================================================
 module operator_radial
  use global
  use fedvr3d_basis_set
  use not_full_operator_index
  implicit none
  ! the matrix element of kinetic operator in radial part
  real(kind=k1),allocatable,save :: tmat_radial(:)
  ! the matrix element of v potential in radial part 
  real(kind=k1),allocatable,save :: vmat_radial(:)

  contains

!==============================================================================
! calc. the kinetic energy operator in fedvr_radial part space
!calc. the kinetic operator in fe-dvr.                                                                e        f
!                                                                                                 d  x (x)  d x (x)
!      e           1   d     f          1                                                             i        j 
!   < x  (x)  | - ~~~ ~~~ | x (x) > = ~~~~~ (delta(e,f) + delta(e,f-1) + delta(e,f+1))* Int { dx * ~~~~~~~ * ~~~~~~~} 
!      i           2   dx    j          2                                                           dx         dx
!
!==============================================================================
    subroutine calc_tmat_radial
     implicit none
     integer :: idicp,e,i,f,j,ee,ii
     integer :: ibasis,jbasis
     real(kind=k1) :: sumtemp
  
    allocate(tmat_radial(n_total_kinetic))
 
     tmat_radial = 0.0_k1

     do idicp =1, n_total_kinetic

        ibasis = index_kinetic_basis(idicp,1); jbasis = index_kinetic_basis(idicp,2) 

! from global to local e,i,f,j
        e = which_element(ibasis)
        i = which_basis(ibasis)
        f = which_element(jbasis)
        j = which_basis(jbasis)
!===================================================================================================
     if( e==f )  then

            if((i/= fedvr3d%fedvr_nb(e)).and.(j/=fedvr3d%fedvr_nb(f))) then        
 
               sumtemp = 0.0_k1
               do ii =1,fedvr3d%fedvr_nb(e)

                 sumtemp = sumtemp + fedvr_w (e,ii) * derivate_x(e,i,e,ii) * derivate_x(f,j,e,ii)

               enddo 


            elseif(((i==fedvr3d%fedvr_nb(e)).and.(j/=fedvr3d%fedvr_nb(f))) &
           .or.((i/=fedvr3d%fedvr_nb(e)).and.(j==fedvr3d%fedvr_nb(f))) ) then

               sumtemp = 0.0_k1
               do ii =1,fedvr3d%fedvr_nb(e)

                 sumtemp = sumtemp + fedvr_w (e,ii) * derivate_x(e,i,e,ii) * derivate_x(f,j,e,ii)

               enddo


            
            elseif((i== fedvr3d%fedvr_nb(e)).and.(j==fedvr3d%fedvr_nb(f))) then            


             sumtemp = 0.0_k1
             do ee =e, e+1

               do ii =1, fedvr3d%fedvr_nb(ee)

                sumtemp = sumtemp + fedvr_w (ee,ii) * derivate_x(e,i,ee,ii) * derivate_x(f,j,ee,ii) 
               enddo 
             enddo           
           endif
!==================================================================================================      
     elseif(e==(f-1)) then


            sumtemp = 0.0_k1
           do ii =1,fedvr3d%fedvr_nb(f)

             sumtemp = sumtemp + fedvr_w (f,ii) * derivate_x(e,i,f,ii) * derivate_x(f,j,f,ii)

           enddo
!===================================================================================================
     elseif(e==(f+1)) then

           sumtemp = 0.0_k1
           do ii =1,fedvr3d%fedvr_nb(e)

                 sumtemp = sumtemp + fedvr_w (e,ii) * derivate_x(e,i,e,ii) * derivate_x(f,j,e,ii)

           enddo
!====================================================================================================
     endif

     tmat_radial(idicp) = sumtemp*0.5d0

     enddo
 return
 end subroutine calc_tmat_radial

!==============================================================================
! calc. the kinetic energy operator in fedvr_radial part space
!==============================================================================
!====================================================================================
! the fedvr1d , fedvr3d is the same one
!    e        f
! < x  | v | x  >  == v ( x(e,i) ) * delta(e,f)*delta(i,j)
!    i        j
!
!                 e     1      f                    e       1         f       e                                   1                                        f
! for fedvr-1d,< x  | ~~~~~ | x   >, for,fedvr-3d <x   |~~~~~~~~~~~| x  > = <x  |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|x  > 
!                 i     r      j                    i      -->        j       i  sqrt((rsin[theta]cos[phi])^2+(rsin[theta]sin[phi])^2) + (r*cos[theta])^2  j
!                                                           r 
!     e     1      f
! = <x  |~~~~~~ | x >
!     i     r      j
!====================================================================================  
 subroutine calc_vmat_radial
  implicit none
    integer :: idicp, e,i

    allocate(vmat_radial(fedvr3d%nb_r))  

    do idicp =1,fedvr3d%nb_r
        
       e = which_element(idicp)
       i = which_basis(idicp)      
       vmat_radial(idicp) = - system%zz/( fedvr_x(e,i ) ) 
     enddo
  return
 end subroutine calc_vmat_radial


!==============================================================================================
! assume : num_of_element, the number of the element 
!          n(e), the number of the basis function in e-th element 
! 
!                   e
! find derivate of x 
!                   i
!
!
! for element 1,2,3....., num_of_element -1, there are i =2,3,4,...,ne
!
! 
!               e                     e        e+1
!              f  (x)                f (x)  + f  (x)  
!    e          i              e      ne       1
!   x  (x) =~~~~~~~~~~~~~ ,   x   = ~~~~~~~~~~~~~~~~~~~  , {i=2, n(e) -1 ,  ne  }             
!    i               e         ne           e     e+1
!             sqrt( w  )             sqrt( w   + w    )
!                    i                      ne    1
!
!       e                 e
!    d x  (x)          d f (x) 
!       i                 i            1
!  ~~~~~~~~~~~  ==  ~~~~~~~~~~~~ * ~~~~~~~~~~~      ,  for i = 2, ne -1    
!     d x                dx               e
!                                  sqrt( w  )
!                                         i
!
!       e                 e              e+1
!    d x  (x)          d f  (x)       d f   (x)
!       ne                ne             1                1
!  ~~~~~~~~~~~  == [ ~~~~~~~~~~~~~ + ~~~~~~~~~~ ] * ~~~~~~~~~~~~~    ,     for i = ne    
!     d x                dx             dx              e     e+1
!                                                 sqrt(w   + w   )   
!                                                       ne    1 
!
!  For num_of_element-th,  the last element have no brige function 
!  
!
!
!               e                    
!              f  (x)                 
!    e          i              
!   x  (x) =~~~~~~~~~~~~~    {i =2, n(num_of_element) - 1 }               
!    i               e         
!             sqrt( w  )            
!                    i                  
!
!
!       e                 e
!    d x  (x)          d f (x) 
!       i                 i            1
!  ~~~~~~~~~~~  ==  ~~~~~~~~~~~~ * ~~~~~~~~~~~      ,  for i = 2, ne -1    
!     d x                dx               e
!                                  sqrt( w  )
!                                         i
! there are no ne  !!!!! Attention here 
!=================================================================================================
!
!
!    e
! d x  (x)   |
!    i       |         ===  derivate_x(e,i,f,j)
!~~~~~~~~~~~ |      f
!   dx       | x = x  
!                   j
!
!=================================================================================================
    double precision function derivate_x(e,i,f,j)
      implicit none
      integer :: e,i,f,j,ne_here
      
      if( e<1 .or. e>fedvr3d%number_of_element .or. i<2 .or. i > fedvr3d%fedvr_nb(e)) then
 
        write(6,*) 'error happen in function deriate_x'
        stop 
      endif
      if((e==fedvr3d%number_of_element) .and. (i > fedvr3d%fedvr_nb(e)-1)) then
 
        write(6,*) 'error happen in function deriate_x'
        stop
      endif
      if( f<1 .or. f>fedvr3d%number_of_element .or. i<1 .or. i > fedvr3d%fedvr_nb(f)) then
 
        write(6,*) 'error happen in function deriate_x'
        stop 
      endif
!
!  1< e < num_of_element - 1 
!  
!
     ne_here = fedvr3d%fedvr_nb(e) 
     if((e>=1) .and. (e<= fedvr3d%number_of_element -1)) then
         if((i>=2) .and. (i<= fedvr3d%fedvr_nb(e) -1 )) then
           derivate_x = derivate_lobatto(e,i,f,j)/(sqrt(fedvr_w(e,i)))
         elseif(i ==fedvr3d%fedvr_nb(e) ) then
          derivate_x = (derivate_lobatto(e,ne_here,f,j) + derivate_lobatto(e+1,1,f,j) ) &
          /(sqrt( fedvr_w(e,ne_here) + fedvr_w(e+1,1) ))
         else
           write(6,*) 'error happen in function derivate_x'
          endif
         elseif(e == fedvr3d%number_of_element) then
         if((i>=2) .and. (i<= fedvr3d%fedvr_nb(e) -1 )) then
           derivate_x = derivate_lobatto(e,i,f,j)/(sqrt(fedvr_w(e,i)))
          else
           write(6,*) 'error happen in function derivate_x'
          endif
     else
       write(6,*) 'error element in function derviate_x'
      endif
  return
 end function derivate_x

! derivative of Lobatto shape function
    ! a speed up is possible by circumventing somehow the calculation of indices
    ! df_{e,i}(x_{f,j})
!=============================================================================================
!calc. the derivative of Lobatto shape function
!
!    e
! d f (x)  |
!    i     |                 =  0, e/=f     
!~~~~~~~~~ |       f
!  d x     | x = x    
!                  j         
!                                 1
!                           = ~~~~~~~~~~*( delta(i, n_e) - delta(i, 1)  ) , i==j , e==f                    
!                                    e
!                                2* w      
!                                    i 
!                                                  e     e
!                                1        ____    x  -  x 
!                         = ~~~~~~~~~~~ *  ||      j     k  
!                             e     e      ||   ~~~~~~~~~~~~ , i /= j , e==f 
!                            x     x     k/=i,j    e      e 
!                             i     j             x   -  x  
!                                                  i      k
!===============================================================================================
    double precision function derivate_lobatto(e,i,f,j)
      use global
      implicit none
      integer :: e,i,f,j,k
      real(kind=k1) :: sum_temp,product_temp 
      if (e/=f) then
         derivate_lobatto=0.0_k1
      else
         if (i==j) then
         
             sum_temp = 0.0_k1
             if(i == fedvr3d%fedvr_nb(e)) then
                 sum_temp = sum_temp + 1.0_k1
             endif
             if(i == 1) then
                 sum_temp = sum_temp - 1.0_k1
             endif
 
             derivate_lobatto=1.0d0/(2.0d0*fedvr_w(e,i))*(sum_temp)
 
           else
            product_temp=1.0_k1
            do k=1,fedvr3d%fedvr_nb(e)
               if ((k/=i).and.(k/=j)) then
                  product_temp=product_temp*(fedvr_x(e,j)-fedvr_x(e,k))/(fedvr_x(e,i)-fedvr_x(e,k))
                endif
            enddo
            derivate_lobatto=product_temp/(fedvr_x(e,i)-fedvr_x(e,j))
         endif
      endif
    end function derivate_lobatto


 end module operator_radial



 
!
! drive subroutine, find the information about the radial part and angular part
!
 subroutine drive_operator_radial
  use operator_radial
  implicit none
!
! radial part, kinetic energy operator and nuclei attractive potential
!
  call calc_tmat_radial
  call calc_vmat_radial

  return
 end subroutine drive_operator_radial
