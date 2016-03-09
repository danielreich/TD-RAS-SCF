!
! contains the information of basis set (fedvr3d) in spherical coordinates
!
 module fedvr3d_basis_set
  use global
  implicit none
  real(kind=k1),allocatable,save :: fedvr_w(:,:),fedvr_x(:,:) 
  real(kind=k1),allocatable,save :: fedvrx_global(:)
  integer,allocatable,save :: which_element(:),which_basis(:) 
  integer,allocatable,save :: lm_l(:),lm_m(:) 
  integer, allocatable, save :: ang_max_core !! Maximum angular function for the core
  integer,allocatable,save :: global_to_local(:,:)

  contains


     
! fedvr_r(1)      fedvr_r(2)          fedvr_r(3)          fedvr_r(4)
!     |________________|__________________|___________________|
!
!
   subroutine initial_fedvr3d_radial
    implicit none
    integer(4) :: ir,i_temp,e,nb_e,i,ie,ibasis,i_all
    real(kind=k1) :: dr_step
    real(kind=k1),allocatable :: a(:),b(:),lxi(:),lwi(:),fedvr_r(:)
    integer(4) :: max_nb,idicp,is,nbt

! num_of_element : number of finite element
!
! find the maximum nodes in the num_of_element elements
!
    max_nb = 0
    do idicp=1, fedvr3d%number_of_element
   
      if(fedvr3d%fedvr_nb(idicp) <= 0 ) then
          write(*,*) 'the basis in element ', idicp, ' is not correct !!!'
          stop '--- error happen in initial_fedvr3d_radial'                                 
      endif  
      if(fedvr3d%fedvr_nb(idicp) > max_nb ) then
        max_nb = fedvr3d%fedvr_nb(idicp)
      endif
    enddo
!
! dr_step: is the length of every element, equally distribution in the radial part
!
     dr_step  = (fedvr3d%r_end - fedvr3d%r_begin)/fedvr3d%number_of_element
! please read the beginning
     allocate ( fedvr_r(fedvr3d%number_of_element+1))
     fedvr_r(1)  = fedvr3d%r_begin
     do ir = 2,fedvr3d%number_of_element +1
       fedvr_r(ir) = fedvr_r(ir-1) +dr_step 
     enddo
!
! figure out the number of the basis function in the radial part
!
     i_temp = 0 
     do e =1,fedvr3d%number_of_element - 1
       do i=2, fedvr3d%fedvr_nb(e)
         i_temp = i_temp + 1
       enddo
     enddo

    fedvr3d%nb_r = i_temp + fedvr3d%fedvr_nb(fedvr3d%number_of_element) - 2

    write(*,*) 
    write(*,*) 'the basis function in the radial part is here ', fedvr3d%nb_r
    write(*,*)

!
! find the weights and nodes in every element
!
    allocate (a(max_nb),b(max_nb),lxi(max_nb),lwi(max_nb))
    allocate ( fedvr_w(fedvr3d%number_of_element,max_nb),fedvr_x(fedvr3d%number_of_element,max_nb))
    allocate(fedvrx_global(fedvr3d%nb_r))

     nbt = 1
     do e=1,fedvr3d%number_of_element
         nb_e =fedvr3d%fedvr_nb(e) 
  
         do i =1,nb_e
            a(i)=0.0d0
            b(i)=(i-1)*(i-1)/((2.0d0*(i-1)-1)*(2.0d0*(i-1)+1.0d0))
         enddo
!
! from sebastain
! calc. the weights and nodes in [-1,1] 
!
         call lobatto(a,b,2.0d0,-1.0d0,1.0d0,nb_e,lxi,lwi)
!
! transform from [-1,1] to [fedvr_r(i)    fedvr_r(i+1)]
!       
         do i=1,nb_e
            fedvr_w(e,i)=0.50d0*(fedvr_r(e+1)-fedvr_r(e))*lwi(i)
            fedvr_x(e,i)=0.50d0*((fedvr_r(e+1)-fedvr_r(e))*lxi(i)+(fedvr_r(e+1)+fedvr_r(e)))
         enddo
!
! remember the total global nodes
!
         is =1
         if(e==1) is =2
         ie = fedvr3d%fedvr_nb(e)-1
         do i=is,ie
            fedvrx_global(nbt)=0.50d0*((fedvr_r(e+1)-fedvr_r(e))*lxi(i)+(fedvr_r(e+1)+fedvr_r(e)))
            nbt=nbt+1
         enddo 
      enddo

      if(nbt- 1 /= fedvr3d%nb_r) then
        write(*,*) 'error happen in fedvr3d_radial'
        stop 
       endif
!
! find the which_element, which_basis
!
!=======================================================================================
! for example, number_of_element = 3, fenb() =4  
!
!
!        1      2       3      4      5      6      7     8
! |______*______*_______|______&______&______|______$_____$_____| 
!        1      1       1      2      2      2      3     3   
!       x      x       x      x      x      x      x     x 
!        2      3       4      2      3      4      2     3
!
! which_element( 5 ) = 2 , the 5-th global radial basis is in 2-the element
! which_basis (5) = 3 ,                                    and 3-the basis function
!======================================================================================
    allocate(which_element(fedvr3d%nb_r),which_basis(fedvr3d%nb_r))

    i_all = 0
    do ie=1,fedvr3d%number_of_element-1
      do ibasis=2, fedvr3d%fedvr_nb(ie)
  
        i_all = i_all +1
      
        which_element(i_all) = ie 
        which_basis(i_all) = ibasis
      enddo
    enddo

    do ibasis = 2,fedvr3d%fedvr_nb(fedvr3d%number_of_element)-1
    
      i_all = i_all +1
      which_element(i_all) = fedvr3d%number_of_element
      which_basis(i_all) = ibasis
    enddo

   deallocate(a,b,lxi,lwi,fedvr_r)
  end subroutine initial_fedvr3d_radial

!
! the basis in the angular part
!
!==========================================================================================
!
! fe-dvr, angular part, basis set { l ,m  } 
! spherical harmonics function
! for example, l_max = 2, m_max = 2
!
! basis set :  l,  m
!              0   0
!              1  -1
!              1   0   
!              1   1
!              2  -2
!              2  -1
!              2   0
!              2   1
!              2   2  
!
! another axample as following, l_max =2, m_max =1
!
! basis set :  l,  m
!              0   0
!              1  -1
!              1   0   
!              1   1
!              2  -1
!              2   0
!              2   1
!==========================================================================================
! calc. the angular basis set, { lm }
!
   subroutine inital_fedvr3d_angle
    implicit none
    integer :: num_lm,isum
    integer :: il,imin,im,idicp
    integer,parameter :: max_angle = 10000 
    integer,allocatable,save :: lm_lp(:),lm_mp(:)
!
! find the combination of (l,m) which satisfied the conditions
!
      allocate(lm_lp(max_angle),lm_mp(max_angle))
      
      num_lm=0
      do il=0,fedvr3d%l_max
         if (il .lt. fedvr3d%m_max) then
            imin=il
         else
            imin=fedvr3d%m_max
         endif
         do im=-imin,imin
            num_lm= num_lm+1

            if(num_lm>max_angle) then
              write(*,*) 'increase the number of angular part basis set'
              stop 
            endif
            lm_lp(num_lm) = il
            lm_mp(num_lm) = im
           enddo
      enddo
           
      fedvr3d%nb_angle = num_lm
      write(*,*)
      write(*,*) 'the number of basis function in the the angular part ', num_lm
      write(*,*) 
      write(*,*) 'the global basis function is ', fedvr3d%nb_r*fedvr3d%nb_angle
      write(*,*)
      allocate(lm_l(fedvr3d%nb_angle),lm_m(fedvr3d%nb_angle))
 !
 ! store the (l,m) combination in the arrays lm_l and lm_m 
 !
      do idicp =1,fedvr3d%nb_angle
        lm_l(idicp) = lm_lp(idicp)
        lm_m(idicp) = lm_mp(idicp)
      enddo    

     deallocate(lm_lp,lm_mp)    

  !! Calculate the maximum angular momentum for the core

     Do idicp=1,,fedvr3d%nb_angle
        If (lm_l(idicp).gt.fedvr3d%l_max_core) then
           ang_max_core=idicip
        End If
     End Do

     print*, ang_max_core
     stop

   end subroutine inital_fedvr3d_angle



!
! from global basis to local basis {radial,angluar}
!  
 subroutine figure_out_global_to_local()
  implicit none
  integer :: isum,idicp,i_k,k,l,m

  allocate(global_to_local(fedvr3d%nb_r*fedvr3d%nb_angle,4))

      isum = 0 
      do idicp =1,  fedvr3d%nb_angle
        do i_k = 1, fedvr3d%nb_r
          isum = isum +1
          l = lm_l(idicp)
          m = lm_m(idicp)
          k = i_k
          global_to_local(isum,1) = k
          global_to_local(isum,2) = l
          global_to_local(isum,3) = m  
          global_to_local(isum,4) = idicp
        enddo
      enddo
  return
 end subroutine figure_out_global_to_local

!! To plot the monoparticular orbitals.

 subroutine plot_phi_fedvr3d(phiuni,xglobal,nb_angle,element,basis,l,m, weights,rho)
  implicit none

! This subroutine produces the density for the monoparticular
! functions.
!!INPUT
 complex(kind=k2), intent(in) :: phiuni(:,:)
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
 complex(8), allocatable, intent(in) :: rho(:,:)
!! First arguments is the annihilation and the second creation
!! wavefunction
real*8, allocatable :: density(:)
!! this is the radial wavefunction squared
!!Auxiliary aspects
integer :: i,j,k,ko, counter
character*6 :: file_orbital
real*8 :: exp_r2,exp_r,exp_rm1 !! <r**2>, <r>, <r**-1>
real*8 :: exp_rm2
real*8 :: exp_l2, exp_m2, exp_m, norm
real*8 :: ionization_single

allocate(density(size(xglobal)))

Do i=1,size(phiuni(1,:)) !! Run in the orbitals

   !Initialize the varibles
   density=0.0d0 !!Density
   exp_r2=0.0d0 !! <r2>
   exp_r=0.0d0 !! <r>
   exp_rm1=0.0d0 !! <r**-1>
   exp_rm2=0.0d0 !! <r**-2>
   exp_l2=0.0d0 !! <L**2>
   exp_m2=0.0d0 !! <M**2>
   exp_m=0.0d0 !! <M>

   counter=0
   write(file_orbital,'(a4,i2.2)') 'phi_', i

   Open(122+i,file=file_orbital)

   Do j=1,nb_angle !! Run for the spherical harmonics      
      Do k=1,size(density) !! Run for the radial points
         counter=counter+1 !! Run in the global basis (radial+angular)
      !! Density of the wavefunction
         density(k)=density(k)+real(conjg(phiuni(counter,i))*phiuni(counter,i))
         exp_l2=exp_l2+real(conjg(phiuni(counter,i))*phiuni(counter,i))*dble(l(j)*(l(j)+1))
         exp_m2=exp_m2+real(conjg(phiuni(counter,i))*phiuni(counter,i))*dble(m(j)*m(j))
         exp_m=exp_m+real(conjg(phiuni(counter,i))*phiuni(counter,i))*dble(m(j))
      End Do
   End Do

!! Calculation of the expectation values

      Do k=1,size(density)
         exp_r2=exp_r2+density(k)*xglobal(k)**2.0d0
         exp_r=exp_r+density(k)*xglobal(k)
         exp_rm1=exp_rm1+density(k)*xglobal(k)**(-1.0d0)
         exp_rm2=exp_rm2+density(k)*xglobal(k)**(-2.0d0)
      End Do
      
      write(122+i,'(a28,i2.2,a3)') ' # This is the density |phi_',i,'|^2'
      write(122+i,*) '# <r**2>=', exp_r2
      write(122+i,*) '# <r>=', exp_r
      write(122+i,*) '# <r**-1>=', exp_rm1
      write(122+i,*) '# <r**-2>=', exp_rm2
      write(122+i,*) '# <L**2>=', exp_l2
      write(122+i,*) '# <M>=', exp_m
      write(122+i,*) '# <M**2>=', exp_m2
      write(122+i,*) ''
      write(122+i,*) ''
      
      write(122+i,*) 0.0d0,0.0d0
   Do k=1,size(density)
      !! Plot the wave function.
      !! Take into account the weights on the points.
      If (k.lt.size(density)) then
         If (basis(k).lt.basis(k+1)) then !! basis(k) is not the last
                                      !! node of the element.
            write(122+i,*) xglobal(k), density(k)/weights(element(k),basis(k))
         Else !! basis(k) is the last node of the element.
            write(122+i,*) xglobal(k), density(k)/(weights(element(k),basis(k))+weights(element(k+1),1))
         Endif
      End If
   End Do
   Close(122+i)
End Do

!! initialize the density

density=0.0d0

!! Calculation of the one particle density.
!! It is normalized to the number of electrons.

Do i=1,size(rho(1,:)) !! loop of the orbitals
   Do j=i+1,size(rho(:,1)) !! loop of the orbitals
      counter=0
      Do ko=1,nb_angle
         Do k=1,size(density) !! run the number of points
            counter=counter+1
            !! Average density
            density(k)=density(k)+2.0d0*real(conjg(phiuni(counter,i))*phiuni(counter,j)*rho(i,j))
         End Do
      End Do
   End Do
   !!For j=i
   counter=0
   Do ko=1,nb_angle
      Do k=1,size(density) !! run the number of points
         counter=counter+1
         !! Average density
         density(k)=density(k)+real(conjg(phiuni(counter,i))*phiuni(counter,i)*rho(i,i))
      End Do
   End Do
End Do

!! Expectation value for the density
!! Initialize the expectation values

exp_r2=0.0d0
exp_r=0.0d0
exp_rm1=0.0d0
exp_rm2=0.0d0
norm=0.0d0
!! Calculation of the expectation values

Do k=1,size(density)
   norm=norm+density(k) !! Norm of the density
   exp_r2=exp_r2+density(k)*xglobal(k)**2.0d0
   exp_r=exp_r+density(k)*xglobal(k)
   exp_rm1=exp_rm1+density(k)*xglobal(k)**(-1.0d0)
   exp_rm2=exp_rm2+density(k)*xglobal(k)**(-2.0d0)
End Do

exp_r2=exp_r2/norm
exp_r=exp_r/norm
exp_rm1=exp_rm1/norm
exp_rm2=exp_rm2/norm

!! calculating ionization for single ionization

ionization_single=0.0d0

Open(510,file='ionization_single')

Do k=size(density),1,-1
   ionization_single=ionization_single+density(k)
   write(510,*) xglobal(k),ionization_single
End Do
close(510)

!! Opening the density file

Open(122,file='rho_r')

!! Writting the heading of the file

write(122,*) '# <r**2>=', exp_r2
write(122,*) '# <r>=', exp_r
write(122,*) '# <r**-1>=', exp_rm1
write(122,*) '# <r**-2>=', exp_rm2
write(122,*) '# Norm=', norm
write(122,*) ''
write(122,*) ''

!! Prepare the density to plot it by dividing by the weights.

   write(122,*) 0.0d0,0.0d0
   Do k=1,size(density)
      !! Plot the wave function.
      !! Take into account the weights on the points.
      If (k.lt.size(density)) then
         If (basis(k).lt.basis(k+1)) then !! basis(k) is not the last
                                      !! node of the element.
            write(122,*) xglobal(k), density(k)/weights(element(k),basis(k))/norm
         Else !! basis(k) is the last node of the element.
            write(122,*) xglobal(k), density(k)/(weights(element(k),basis(k))+weights(element(k+1),1))/norm
         Endif
      End If
   End Do
   Close(122)

Deallocate(density)

End subroutine plot_phi_fedvr3d


 subroutine plot_phi2_fedvr3d(phiuni,xglobal,nb_angle,element,basis,l,m, weights,rho2)
  implicit none

! This subroutine produces the two particle density distribution.
! rho(r_1,r_2)
!!INPUT
 complex(kind=k2), intent(in) :: phiuni(:,:)
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
 complex(8), allocatable, intent(in) :: rho2(:,:,:,:)
!! First arguments is the annihilation and the second creation
!! wavefunction
real*8, allocatable :: density2(:,:),density1(:)
real*8 :: ionization_double
!! density2(r1:r2) = rho(r1, r2)
!!Auxiliary aspects
real*8 :: norm,aux
complex(8) :: auxcomplex
integer :: i,j,k,ko, counter
integer :: kk1, kk2
integer :: i2, j2, k2, ko2, counter2
integer :: p, q, r, s !! index for orbitals
character*6 :: file_orbital

allocate(density2(1:size(xglobal),1:size(xglobal)))

!! initialize the density

density2=0.0d0
norm=0.0d0

!! Obtaining rho(r1,r2). At this point it is written in the DVR

Do p=1, size(rho2(:,1,1,1)) !! Loop in the orbitals a_p+
   Do r=1,size(rho2(1,:,1,1)) !! Loop in the orbitals a_r+
      Do s=1,size(rho2(1,1,:,1)) !! Loop in the orbitals a_s
         Do q=1, size(rho2(1,1,1,:)) !! Loop in the orbitals a_q
            
            counter=0
            Do ko=1,nb_angle !!running in the basis to separate the
                             !! angle and the DVR functions
               Do k=1,size(density2(:,1)) !! run the number of points of r1
                  counter=counter+1
                  
                  counter2=0
                  Do ko2=1,nb_angle !!running in the basis to separate the
                                    !! angle and the DVR functions
                     Do k2=1,size(density2(1,:)) !! run the number of points of r2
                        counter2=counter2+1

       auxcomplex=conjg(phiuni(counter2,p))*conjg(phiuni(counter,r))
       auxcomplex=auxcomplex*phiuni(counter,s)*phiuni(counter2,q)
       aux=real(auxcomplex*rho2(p,r,s,q))
       aux=aux/2.0d0 !! due to the factor out of the sum
                        
       density2(k,k2)=density2(k,k2)+aux
    End Do
                  End Do
               End Do
            End Do
            
            
         End Do
      End Do
   End Do
End Do

!! Calculate the norm

Do k=1,size(density2(:,1))
   Do k2=1,size(density2(1,:))
      norm=norm+density2(k,k2)
   End Do
End Do

!! calculating ionization for double ionization

ionization_double=0.0d0

Open(510,file='ionization_double')

Do kk1=size(density2(:,1)),1,-1
   ionization_double=ionization_double+density2(kk1,kk1)
      write(510,*) xglobal(kk1),xglobal(kk1),ionization_double
   Do kk2=kk1-1,1,-1
      ionization_double=ionization_double+2.0d0*density2(kk1,kk2)
      write(510,*) xglobal(kk1),xglobal(kk2),ionization_double
   End Do
   write(510,*)
End Do
close(510)


!! Opening the density file

      Open(122,file='rho_r1_r2')

!! Writting the heading of the file

      write(122,*) '#norm =', norm
write(122,*) ''
   write(122,*) '# r1   r2     rho(r1,r2)'



!! Prepare the density to plot it by dividing by the weights.
   write(122,*) 0.0d0,0.0d0,0.0d0
   Do k2=1,size(density2(1,:))-1
      write(122,*) 0.0d0,xglobal(k2),0.0d0
   End Do
   
   write(122,*) ''

   Do k=1,size(density2(1,:))-1 !! it is not plotting the last point
      write(122,*) xglobal(k),0.0d0, 0.0d0
      Do k2=1,size(density2(:,1))-1 !! it is not plotting the last point
         !! Plot rho2
         !! Take into account the weights on the points.
         If (basis(k).lt.basis(k+1).and.basis(k2).lt.basis(k2+1)) then 
            !! basis(k) and basis(k2) are not the last nodes of the element.
            write(122,*) xglobal(k),xglobal(k2), density2(k,k2)/(weights(element(k),basis(k))*weights(element(k2),basis(k2)))/norm
         Elseif (basis(k).gt.basis(k+1).and.basis(k2).lt.basis(k2+1)) then 
            !! basis(k) is the last node of the element, but not basis(k2).
            aux=((weights(element(k),basis(k))+weights(element(k+1),1))*weights(element(k2),basis(k2)))
            write(122,*) xglobal(k),xglobal(k2), density2(k,k2)/(aux*norm)
         Elseif (basis(k2).gt.basis(k2+1).and.basis(k).lt.basis(k+1)) then 
            !! basis(k2) is the last node of the element, but not basis(k).
            aux=((weights(element(k2),basis(k2))+weights(element(k2+1),1))*weights(element(k),basis(k)))
            write(122,*) xglobal(k),xglobal(k2), density2(k,k2)/(aux*norm)
         Elseif (basis(k2).gt.basis(k2+1).and.basis(k).gt.basis(k+1)) then 
            aux=((weights(element(k2),basis(k2))+weights(element(k2+1),1))*(weights(element(k),basis(k))+weights(element(k+1),1)))
            write(122,*) xglobal(k),xglobal(k2), density2(k,k2)/(aux*norm)
         Endif
      End Do
      write(122,*) ''
   End Do
   Close(122)

Deallocate(density2)

End subroutine plot_phi2_fedvr3d

!! This subroutine calculates <Psi|dipole|Psi>

 subroutine dipole_expectation(phiuni,xglobal,l,m,rho,xmat_3d,ymat_3d,zmat_3d,nb_angle,xdipole,ydipole,zdipole)
  implicit none

! This subroutine produces the density for the monoparticular
! functions.
!!INPUT
 complex(kind=k2), intent(in) :: phiuni(:,:)
!! phi in terms of the fedvr global basis.
 real*8, intent(in), allocatable :: xglobal(:)
!! xglobal is the radial basis fedvr
 integer, intent(in), allocatable :: l(:),m(:)
!! L and M for the angular function in the argument.
 complex(8), allocatable, intent(in) :: rho(:,:)
!! First arguments is the creation and the second is the annihilation operator
!! wavefunction <Psi|a_1+a_2|Psi>
real*8, allocatable, intent(in) :: xmat_3d(:,:), ymat_3d(:,:), zmat_3d(:,:)
!! dipole matrix elements
 integer, intent(in) :: nb_angle
!! number of angular functions
!!OUTPUT 
real*8 :: xdipole, ydipole, zdipole
!! x, y and z
!!Auxiliary aspects
integer :: i,j,k,ko
integer :: i1,i2, counter1, counter2

!! AQUI: IT IS ONLY DONE FOR THE z COMPONENT

!! initialize the value for t=0

xdipole=0.0d0
ydipole=0.0d0
zdipole=0.0d0

Do i1=1,size(phiuni(1,:)) !! Run in the orbitals for bra

   Do i2=1,size(phiuni(1,:)) !! Run in the orbitals for ket  !! AQUI !! This loop can be improved by starting with i2=i1,...

      Do j=1,nb_angle !! Run for the spherical harmonics for the bra
         
         !! Loop for the same angular function
!           |                                |
!           V                                V

         Do ko=1,size(xglobal) !! Run for the radial points            
 
            counter1=(j-1)*size(xglobal)+ko
           
            zdipole=zdipole+zmat_3d(j,j)*(conjg(phiuni(counter1,i1))*phiuni(counter1,i2)*rho(i1,i2))*xglobal(ko)
         
!!$            write(111,*) j,j,ko,zdipole
   
         End Do
         
         !! Loop for the rest of angular functions
                  
         Do k=j+1,nb_angle !! Run for the spherical harmonics for the ket
            
            If (abs(l(j)-l(k)).ge.2) exit           
            If (abs(m(j)-m(k)).ge.2) exit
            
            Do ko=1,size(xglobal) !! Run for the radial points

            counter1=(j-1)*size(xglobal)+ko
            counter2=(k-1)*size(xglobal)+ko
                             
            zdipole=zdipole+zmat_3d(j,k)*2.d0*real(conjg(phiuni(counter1,i1))*phiuni(counter2,i2)*rho(i1,i2))*xglobal(ko)

!!$            write(111,*) j,k,ko,zdipole,xglobal(ko)

            End Do
            
         End Do

      End Do
   End Do
End Do

End subroutine dipole_expectation


end module fedvr3d_basis_set
 
!
! drive subroutine, find the information about the radial part and angular part
!
 subroutine drive_fedvr3d_basis_set
  use fedvr3d_basis_set
  implicit none

  write(*,*) '======begin initialization of the basis set, radial part======='
! radial part
  call initial_fedvr3d_radial
  write(*,*) '======begin initialization of the basis set , angular part====='
! angular part
  call inital_fedvr3d_angle

! global_to_local array
  call figure_out_global_to_local 
 
  return
 end subroutine drive_fedvr3d_basis_set






!
! taken from sebastain's code
!

module anglib
! Library of angular momentum coupling coefficient routines in fortran 90
! Paul Stevenson, Oxford University/Oak Ridge National Laboratory.
! spaul@mail.phy.ornl.gov
  integer, parameter :: rk = selected_real_kind(p=15)
contains
  function cleb(j1,m1,j2,m2,j,m)
    implicit none
    ! calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
    ! arguments are integer and twice the true value. 
    double precision    :: cleb,factor,sum
    integer :: j1,m1,j2,m2,j,m,par,z,zmin,zmax
    ! some checks for validity (let's just return zero for bogus arguments)
    if (2*(j1/2)-int(2*(j1/2.0)) /= 2*(abs(m1)/2)-int(2*(abs(m1)/2.0)) .or. &
         2*(j2/2)-int(2*(j2/2.0)) /= 2*(abs(m2)/2)-int(2*(abs(m2)/2.0)) .or. &
         2*(j/2)-int(2*(j/2.0)) /= 2*(abs(m)/2)-int(2*(abs(m)/2.0)) .or. &
         j1<0 .or. j2<0 .or. j<0 .or. abs(m1)>j1 .or. abs(m2)>j2 .or.&
         abs(m)>j .or. j1+j2<j .or. abs(j1-j2)>j .or. m1+m2/=m) then
       cleb= 0.0d0
    else
    
       factor = 0.0d0
       factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2)
       factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
       factor = sqrt(factor)
       
       zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
       
       sum=0.0d0
       do z = zmin,zmax
          par=1
          if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
          sum=sum+par*binom((j1+j2-j)/2,z)*binom((j1-j2+j)/2,(j1-m1)/2-z)*&
               binom((-j1+j2+j)/2,(j2+m2)/2-z)
       end do
       
       cleb = factor*sum
    end if
  end function cleb
  function sixj(a,b,c,d,e,f)
    implicit none
    integer, intent(in) :: a,b,c,d,e,f
    double precision :: sixj
    integer :: phase, nlo, nhi, n
    double precision :: outfactors, sum, sumterm
    ! calculates a Wigner 6-j symbol. Argument a-f are integer and are
    ! twice the true value of the 6-j's arguments, in the form
    ! { a b c }
    ! { d e f }
    ! Calculated using binomial coefficients to allow for (reasonably) high
    ! arguments.
    ! First check for consistency of arguments:
    sixj=0.0d0
    if(mod(a+b,2)/=mod(c,2)) return
    if(mod(c+d,2)/=mod(e,2)) return
    if(mod(a+e,2)/=mod(f,2)) return
    if(mod(b+d,2)/=mod(f,2)) return
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(c-d)>e .or. c+d<e) return
    if(abs(a-e)>f .or. a+e<f) return
    if(abs(b-d)>f .or. b+d<f) return
    phase=(-1)**((a+c+d+f)/2)
    outfactors = angdelta(a,e,f)/angdelta(a,b,c)
    outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e)
!    write(6,*) outfactors
    nlo = max( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
    nhi = min( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)
    sum=0.0
    do n=nlo,nhi
       sumterm = (-1)**n
       sumterm = sumterm * binom(n+1,n-(a+b+c)/2)
       sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2)
       sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2)
       sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2)
!       write(6,*) ',sumterm: ',sumterm
       sum=sum+sumterm
    end do
    sixj = phase * sum * outfactors
  end function sixj
  function angdelta(a,b,c)
    implicit none
    integer :: a,b,c
    double precision    :: angdelta, scr1
    ! calculate the function delta as defined in varshalovich et al. for
    ! use in 6-j symbol:
    scr1= factorial((a+b-c)/2)
    scr1=scr1/factorial((a+b+c)/2+1)
    scr1=scr1*factorial((a-b+c)/2)
    scr1=scr1*factorial((-a+b+c)/2)
    angdelta=sqrt(scr1)
  end function angdelta
  function ninej(a,b,c,d,e,f,g,h,i)
    implicit none
    integer :: a,b,c,d,e,f,g,h,i
    double precision    :: ninej, sum
    integer :: xlo, xhi
    integer :: x
    ! calculate a 9-j symbol. The arguments are given as integers twice the
    ! value of the true arguments in the form
    ! { a b c }
    ! { d e f }
    ! { g h i }
    ninej=0.0
    ! first check for bogus arguments (and return zero if so)
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(d-e)>f .or. d+e<f) return
    if(abs(g-h)>i .or. g+h<i) return
    if(abs(a-d)>g .or. a+d<g) return
    if(abs(b-e)>h .or. b+e<h) return
    if(abs(c-f)>i .or. c+f<i) return
    
    xlo = max(abs(b-f),abs(a-i),abs(h-d))
    xhi = min(b+f,a+i,h+d)
    
    sum=0.0d0
    do x=xlo,xhi,2
       sum=sum+(-1)**x*(x+1)*sixj(a,b,c,f,i,x)*sixj(d,e,f,b,x,h)*&
            sixj(g,h,i,x,a,d)
    end do
    ninej=sum
  end function ninej
  recursive function factorial(n) result(res)
    implicit none
    integer :: n
    double precision :: res
    if (n==0 .or. n==1) then
       res=1.0
    else
       res=n*factorial(n-1)
    end if
  end function factorial
  recursive function binom(n,r) result(res)
    implicit none
    integer :: n,r
    double precision :: res
    if(n==r .or. r==0) then
       res = 1.0d0
    else if (r==1) then
       res = real(n,rk)
    else
       res = real(n,rk)/real(n-r,rk)*binom(n-1,r)
    end if
  end function binom
end module anglib
    double precision function fedvr3dbase_angpart(L,l1,m1,l2,m2,l3,m3,l4,m4)
      implicit none
      integer :: L
      integer :: l1,l2,l3,l4
      integer :: m1,m2,m3,m4
  real(8) :: pi, getgaunt
  pi = 3.1415926d0
      fedvr3dbase_angpart=(-1.0d0)**(m4-m3)*4.0d0*pi/(2.0d0*L+1.0d0)*&
           getgaunt(l1,l2,l,m1,m2,m1-m2)*getgaunt(l3,l,l4,m3,m3-m4,m4)
       end function fedvr3dbase_angpart

!
! taken from sebastain 's code
!
    double precision function getgaunt(l1,l2,l3,m1,m2,m3)
      use anglib
      implicit none
      integer::l1,l2,l3,m1,m2,m3
      real(8) :: pi
      pi = 3.1415926d0
      
      if (l2.lt.abs(l1-l3)) then
         getgaunt=0.0d0
         return
      endif
      if (l2.gt.l1+l3) then
         getgaunt=0.0d0
         return
      endif
      if (m1.ne.(m2+m3)) then
         getgaunt=0.0d0
         return
      endif
      if ((abs(m1).gt.l1) .or. (abs(m2).gt.l2) .or. (abs(m3).gt.l3)) then
         getgaunt=0.0d0
         return
      endif
      if (mod(l1+l2+l3,2).ne.0) then
         getgaunt=0.0d0
         return
      endif
      getgaunt= cleb(2*l2,2*m2,2*l3,2*m3,2*l1,2*m1)*cleb(2*l2,0,2*l3,0,2*l1,0)*&
           sqrt((2.0d0*l2+1.0d0)*(2.0d0*l3+1.0d0)/(4.0d0*pi*(2.0d0*l1+1.0d0)))
    end function getgaunt


    ! subroutine calculates lobatto points and weigths within interval a,b
    ! input:
    ! nbe: number of basis functions in element under consideration
    ! adapted from Numerical Recipes
    ! gets lxi,lwi arrays, on return contains Lobatto points and weights
    
    subroutine lobatto(a,b,amu0,x1,xn,n,lxi,lwi)
      implicit none
      double precision :: amu0,x1,xn
      integer :: n,i
      double precision :: a(n), b(n)
      double precision :: pl, pr,pm1l,pm1r,p1l,p1r
      double precision :: det
      double precision :: lxi(n), lwi(n)
      
      pl=x1-a(1)
      pr=xn-a(1)
      pm1l=1.0d0
      pm1r=1.0d0
      p1l=pl
      p1r=pr

      do i=2,n-1
         pl=(x1-a(i))*p1l-b(i)*pm1l
         pr=(xn-a(i))*p1r-b(i)*pm1r
         pm1l=p1l
         pm1r=p1r
         p1l=pl
         p1r=pr   
      enddo
      
      det=pl*pm1r-pr*pm1l
      a(n)=(x1*pl*pm1r-xn*pr*pm1l)/det
      b(n)=(xn-x1)*pl*pr/det
      
      call gaucof(a,b,amu0,n,lxi,lwi);
      
    end subroutine lobatto
    
    ! subroutine calculates Gauss coefficients for integration
    ! adapted from numerical recipes
    subroutine gaucof(a,b,amu0,n,lxi,lwi)
      implicit none
      integer n
      double precision :: lxi(n), lwi(n)
      double precision :: a(n),b(n)
      double precision :: amu0
      double precision :: M(n,n)
      
      integer :: i
      integer :: info
      integer :: lwork
      double precision :: work(3*n-1)
      double precision :: w(n)
      
      ! initialize array to zero
      M=0.0d0

      do i=1,n
         if (i .gt. 1)  b(i)=sqrt(b(i));
      enddo
      ! building matrix for diagonalization
      do i=1,n
         M(i,i)=a(i)
         if (i .gt. 1) M(i,i-1)=b(i)
         if (i .lt. n) M(i,i+1)=b(i+1)
      enddo
      
      lwork=3*n-1
      ! find Eigenvalues and Eigenvectors
      call dsyev('V','U',n,M,n,W,work,lwork,info )
      
      do i =1,n
         lxi(i)=W(i)
         lwi(i)=amu0*M(1,i)*M(1,i)
      enddo
    end subroutine gaucof

 




