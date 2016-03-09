!
! module solver, to solver the equation of A and phi
!
 module solver
   use global
   use coffa
   use public_data
   use density
   implicit none

   complex(kind=k2),allocatable,save :: suma0(:)
   complex(kind=k2),allocatable,save :: sumphi0(:,:)
   complex(kind=k2),allocatable,save :: suma1(:,:)
   complex(kind=k2),allocatable,save :: sumphi1(:,:,:)
   complex(kind=k2),allocatable,save :: suma3(:,:)
   complex(kind=k2),allocatable,save :: sumphi3(:,:,:)
     
!  solver1 
   real(kind=k1),dimension(5),save :: s1a1
   real(kind=k1),dimension(4),save :: s1a2
   real(kind=k1),dimension(3),save :: s1a3
   real(kind=k1),dimension(2),save :: s1a4
   real(kind=k1),dimension(1),save :: s1a5
   real(kind=k1),dimension(4),save :: s1c
   real(kind=k1),dimension(5),save :: s1b
 ! solver3
   real(kind=k1),dimension(12),save :: s3a1
   real(kind=k1),dimension(11),save :: s3a2
   real(kind=k1),dimension(10),save :: s3a3
   real(kind=k1),dimension(9), save :: s3a4
   real(kind=k1),dimension(8), save :: s3a5
   real(kind=k1),dimension(7), save :: s3a6
   real(kind=k1),dimension(6), save :: s3a7
   real(kind=k1),dimension(5), save :: s3a8
   real(kind=k1),dimension(4), save :: s3a9 
   real(kind=k1),dimension(3), save :: s3a10   
   real(kind=k1),dimension(2), save :: s3a11
   real(kind=k1),dimension(1), save :: s3a12
   real(kind=k1),dimension(11),save :: s3c
   real(kind=k1),dimension(7),save :: s3b
   complex(kind=k2),save :: chstep

   contains
!=========================runge-kutta 4 order ======================================
   subroutine solver0(t)
    implicit none
    real(kind=k1),intent(in) :: t
    integer(4) :: idicp,jdicp,itime

! 
!==========================k1========================================================
!

   
   call derivate(t)



! acoff    

   do idicp =1,n_string_ab
	 acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(0.5_k1)*chstep  
     suma0(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(1.0_k1/6.0_k1)*chstep
   enddo


! phi
   do jdicp =1,ras1%nptot
   do idicp =1, numbasis

     phi(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(0.5_k1)*chstep
     sumphi0(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(1.0_k1/6.0_k1)*chstep
	enddo
  enddo


! 
!==========================k2========================================================
!
  
  call  derivate(t+0.5_k1*prop1%hstep)
 

! acoff
  do idicp =1,n_string_ab
	 acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(0.5_k1*chstep)  
     suma0(idicp) = suma0(idicp) + ka(idicp)*dcmplx(chstep/3.0_k1)
  enddo
! phi
  do idicp =1, numbasis
    do jdicp =1,ras1%nptot

     phi(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(0.5_k1*chstep)
     sumphi0(idicp,jdicp) = sumphi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(chstep/3.0_k1)
	enddo
  enddo
! 
!==========================k3========================================================
!
  call   derivate(t+0.5_k1*prop1%hstep)

! acoff
  do idicp =1,n_string_ab
	 acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(chstep)  
     suma0(idicp) = suma0(idicp) + ka(idicp)*dcmplx(chstep/3.0_k1)
  enddo
! phi
  do idicp =1, numbasis
    do jdicp =1,ras1%nptot

     phi(idicp,jdicp) = phi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(chstep)
     sumphi0(idicp,jdicp) = sumphi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(chstep/3.0_k1)
	enddo
  enddo

! 
!==========================k4========================================================
!
  call  derivate(t+prop1%hstep)
! acoff
  do idicp =1,n_string_ab  
     acoff(idicp) = suma0(idicp) + ka(idicp)*dcmplx(chstep/6.0_k1)
  enddo
! phi
  do idicp =1, numbasis
    do jdicp =1,ras1%nptot
     phi(idicp,jdicp) = sumphi0(idicp,jdicp) + kphi(idicp,jdicp)*dcmplx(chstep/6.0_k1)
    enddo
  enddo
 return
end subroutine solver0

!
!========================solver_a ========================================================
!
!=========================runge-kutta 4 order ======================================
   subroutine solver_a(t)
    implicit none
    real(kind=k1),intent(in) :: t
    integer(4) :: idicp,jdicp,itime

! 
!==========================k1========================================================
!
   


   call derivate_a(t)

! acoff    

   do idicp =1,n_string_ab
	 acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(0.5_k1)*chstep  
     suma0(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(1.0_k1/6.0_k1)*chstep
   enddo
! 
!==========================k2========================================================
!
  
  call  derivate_a(t+0.5_k1*prop1%hstep)
 

! acoff
  do idicp =1,n_string_ab
	 acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(0.5_k1*chstep)  
     suma0(idicp) = suma0(idicp) + ka(idicp)*dcmplx(chstep/3.0_k1)
  enddo



! 
!==========================k3========================================================
!
  call   derivate_a(t+0.5_k1*prop1%hstep)

! acoff
  do idicp =1,n_string_ab
	 acoff(idicp) = acoff0(idicp) + ka(idicp)*dcmplx(chstep)  
     suma0(idicp) = suma0(idicp) + ka(idicp)*dcmplx(chstep/3.0_k1)
  enddo
! 
!==========================k4========================================================
!
 
  call  derivate_a(t+prop1%hstep)
! acoff
  do idicp =1,n_string_ab  
     acoff(idicp) = suma0(idicp) + ka(idicp)*dcmplx(chstep/6.0_k1)
  enddo

 
 return

end subroutine solver_a





!
!========================runge-kutta 5(6) order ===========================================
!
subroutine solver1(t)
  implicit none
    integer(4) :: idicp,jdicp,kdicp
    real(kind=k1),intent(in) :: t


! initilzation suma,sumphi
    do idicp =1,6
     suma1(:,idicp) = acoff0(:)
     sumphi1(:,:,idicp) = phi0(:,:)
    enddo

!
!===============k1====================================================== 
!  
    call derivate(t)
    call subsolver1(1,99,s1a1,s1b(1))
!
!==============k2=======================================================
!
   call derivate(t+s1c(1)*prop1%hstep)
   call subsolver1(2,0,s1a2,0.0_k1)
!
!==============k3=======================================================
!
  call derivate(t+s1c(2)*prop1%hstep )
  call subsolver1(3,99,s1a3,s1b(2))
!
!==============k4=======================================================
!
   call derivate(t+s1c(3)*prop1%hstep)       		
   call subsolver1(4,99,s1a4,s1b(3))
!
!==============k5=======================================================
!
   call derivate(t+s1c(4)*prop1%hstep)  
   call subsolver1(5,99,s1a5,s1b(4))

   call derivate(t+prop1%hstep)

!acoff ===============================================================================================
   do idicp =1,n_string_ab
    acoff(idicp) = suma1(idicp,6) + ka(idicp) * dcmplx(s1b(5)*chstep)
   enddo
! phi
   do idicp =1,numbasis
       do jdicp =1,ras1%nptot
         phi(idicp,jdicp) = sumphi1(idicp,jdicp,6) + kphi(idicp,jdicp)*dcmplx(s1b(5)*chstep)
      enddo
 enddo

 return
 end subroutine solver1


!
!========================runge-kutta 7(8) order ===========================================
!
 subroutine solver3(t)
    implicit none
    real(kind=k1),intent(in) :: t
    integer(4) :: idicp,jdicp
    
!initial suma3,sumphi3
    do idicp =1,14
      suma3(:,idicp) = acoff0(:)
      sumphi3(:,:,idicp) = phi0(:,:)
    enddo

!========================k1===============================================================	 
    call derivate(t) 
    call subsolver3(1,99,s3a1,s3b(1))

!========================k2===============================================================  
  	call derivate(t+s3c(1)*prop1%hstep) 
    call subsolver3(2,0,s3a2,0.0d0)
  
!========================k3===============================================================   
    call derivate(t+s3c(2)*prop1%hstep)
    call subsolver3(3,0,s3a3,0.0d0)

!========================k4===============================================================
    call derivate(t+s3c(3)*prop1%hstep)      
    call subsolver3(4,0,s3a4,0.0d0)

!========================k5===============================================================
    call derivate(t+s3c(4)*prop1%hstep)
    call subsolver3(5,0,s3a5,0.0d0)

!========================k6===============================================================
    call derivate(t+s3c(5)*prop1%hstep)
    call subsolver3(6,99,s3a6,s3b(2))
	
!========================k7===============================================================
    call derivate(t+s3c(6)*prop1%hstep)
    call subsolver3(7,99,s3a7,s3b(3))

!========================k8===============================================================
    call derivate(t+s3c(7)*prop1%hstep)
    call subsolver3(8,99,s3a8,s3b(4))

!========================k9===============================================================
    call derivate(t+s3c(8)*prop1%hstep)
    call subsolver3(9,99,s3a9,s3b(5))

!========================k10===============================================================
    call derivate(t+s3c(9)*prop1%hstep)
    call subsolver3(10,99,s3a10,s3b(6))

!========================k11===============================================================   
    call derivate(t+s3c(10)*prop1%hstep)
    call subsolver3(11,99,s3a11,s3b(7))

!========================k12===============================================================
    !call derivate(t+s3c(11)*prop1%hstep)
    !call subsolver3(12,0,s3a12,0.0d0)

!========================k13===============================================================
    !call derivate(t+prop1%hstep)


!=========================================================================================
!acoff
   do idicp =1,n_string_ab
    acoff(idicp) = suma3(idicp,14) 
   enddo
! phi
   do idicp =1,numbasis
       do jdicp =1,ras1%nptot
         phi(idicp,jdicp) = sumphi3(idicp,jdicp,14) 
	   enddo
	enddo

  return
 end subroutine solver3


!
! runge-kutta 5(6) order
!
 subroutine subsolver1(k11,bb,s1a,sbb)
    implicit none
    integer(4),intent(in) :: bb,k11
    integer(4) :: mdicp,idicp,kdicp,jdicp
    real(kind=k1),intent(in),dimension(6-k11) :: s1a
    real(kind=k1),intent(in) :: sbb
!=======================================================================================
! acoff
    do idicp =1,n_string_ab   
     acoff(idicp) = suma1(idicp,k11) + ka(idicp)*dcmplx(s1a(1)*chstep)
     if(bb/=0) then
	   suma1(idicp,6) = suma1(idicp,6) + ka(idicp)*dcmplx(sbb*chstep)
	 endif 
     mdicp = 1
     do jdicp =k11+1,5  
        mdicp = mdicp +1
        suma1(idicp,jdicp) = suma1(idicp,jdicp) + ka(idicp)*dcmplx(s1a(mdicp)*chstep)
     enddo
   enddo
! phi
   do idicp =1,numbasis
       do jdicp =1,ras1%nptot
	    phi(idicp,jdicp) = sumphi1(idicp,jdicp,k11) + kphi(idicp,jdicp)*dcmplx(s1a(1)*chstep)
        if(bb/=0) then
	      sumphi1(idicp,jdicp,6) = sumphi1(idicp,jdicp,6) + kphi(idicp,jdicp)*dcmplx(sbb*chstep)
	    endif
		mdicp =1
        do kdicp = k11+1,5
		 mdicp = mdicp +1
         sumphi1(idicp,jdicp,kdicp) = sumphi1(idicp,jdicp,kdicp) + kphi(idicp,jdicp)*dcmplx(s1a(mdicp)*chstep)
		enddo 
	   enddo
	 enddo
  return
 end subroutine subsolver1


  subroutine subsolver3(k11,bb,s1a,sbb)
    implicit none
    integer(4),intent(in) :: bb,k11
    integer(4) :: mdicp,idicp,kdicp,jdicp
    real(kind=k1),intent(in),dimension(14-k11) :: s1a
    real(kind=k1),intent(in) :: sbb
!=======================================================================================
! acoff

    loop1:do idicp =1,n_string_ab
      acoff(idicp) = suma3(idicp,k11) + ka(idicp)*dcmplx(s1a(1)*chstep)

      if(bb/=0) then
	    suma3(idicp,14) = suma3(idicp,14) + ka(idicp)*dcmplx(sbb*chstep)
	  endif 

      mdicp = 1
      do jdicp =k11+1,13  
        mdicp = mdicp +1
        suma3(idicp,jdicp) = suma3(idicp,jdicp) + ka(idicp)*dcmplx(s1a(mdicp)*chstep)
      enddo

   enddo loop1
! phi
   loop2:do idicp =1,numbasis
     
       do jdicp =1,ras1%nptot
	 
	    phi(idicp,jdicp) = sumphi3(idicp,jdicp,k11) + kphi(idicp,jdicp)*dcmplx(s1a(1)*chstep)
	    
        if(bb/=0) then
	      sumphi3(idicp,jdicp,14) = sumphi3(idicp,jdicp,14) + kphi(idicp,jdicp)*dcmplx(sbb*chstep)
	    endif


		mdicp =1
        do kdicp = k11+1,13
		 mdicp = mdicp +1
         sumphi3(idicp,jdicp,kdicp) = sumphi3(idicp,jdicp,kdicp) + kphi(idicp,jdicp)*dcmplx(s1a(mdicp)*chstep)

		enddo 

	   
	   enddo
	 enddo loop2

  return
 end subroutine subsolver3

end module solver







		    
 subroutine initial_solver
  use global
  use solver
  implicit none

! solver1 parameter
!================================================================
 !a21-a61
  s1a1(1) = 0.2_k1 
  s1a1(2) = 3.0_k1/40.0_k1
  s1a1(3) = 44.0_k1/45.0_k1
  s1a1(4) = 19372.0_k1/6561.0_k1
  s1a1(5) = 9017.0_k1/3168.0_k1
!a32 - a62
  s1a2(1) = 9.0_k1/40.0_k1
  s1a2(2) = -56.0_k1/15.0_k1        
  s1a2(3) = -25360.0_k1/2187.0_k1
  s1a2(4) = -355.0_k1/33.0_k1

!a43 - a63							
  s1a3(1) = 32.0_k1/9.0_k1 
  s1a3(2) = 64448.0_k1/6561.0_k1
  s1a3(3) = 46732.0_k1/5247.0_k1

!a54 -a64
  s1a4(1) = -212.0_k1/729.0_k1
  s1a4(2) =  49.0_k1/176.0_k1                             

! a65
  s1a5(1) = -5103.0_k1/18656.0_k1

! c2-c5
 s1c(1) = 0.2_k1
 s1c(2) = 0.3_k1
 s1c(3) = 0.8_k1
 s1c(4) = 8.0_k1/9.0_k1

! b1,b3,b4,b5,b6

 s1b(1) = 35.0_k1/384.0_k1
 s1b(2) = 500.0_k1/1113.0_k1
 s1b(3) = 125.0_k1/192.0_k1
 s1b(4) = -2187.0_k1/6784.0_k1
 s1b(5) = 11.0_k1/84.0_k1

!================================================================

  
   s3a1(1) = 2.0_k1/27.0_k1  ! a21
   s3a1(2) = 1.0_k1/36.0_k1  ! a31
   s3a1(3) = 1.0_k1/24.0_k1  ! a41
   s3a1(4) = 5.0_k1/12.0_k1  ! a51 
   s3a1(5) = 1.0_k1/20.0_k1  ! a61
   s3a1(6)  = -25.0_k1/108.0_k1    ! a71
   s3a1(7)  = 31.0_k1/300.0_k1     ! a81
   s3a1(8)  = 2.0_k1               ! a91
   s3a1(9)  = -91.0_k1/108.0_k1    ! a101 
   s3a1(10) = 2383.0_k1/4100_k1    ! a111
   s3a1(11) = 3.0_k1/205.0_k1      ! a121 
   s3a1(12) = -1777.0_k1/4100.0_k1 ! a131 


   s3a2(1) = 1.0_k1/12.0_k1  ! a32
   s3a2(2) = 0.0_k1     
   s3a2(3) = 0.0_k1
   s3a2(4) = 0.0_k1
   s3a2(5) = 0.0_k1
   s3a2(6) = 0.0_k1
   s3a2(7) = 0.0_k1
   s3a2(8) = 0.0_k1
   s3a2(9) = 0.0_k1
   s3a2(10) = 0.0_k1
   s3a2(11) = 0.0_k1        !a132
   

! a43----a133
   s3a3(1) =1.0_k1/8.0_k1
   s3a3(2) = -25.0_k1/16.0_k1
   s3a3(3) = 0.0_k1
   s3a3(4) = 0.0_k1
   s3a3(5) = 0.0_k1
   s3a3(6) = 0.0_k1
   s3a3(7) = 0.0_k1
   s3a3(8) = 0.0_k1
   s3a3(9) = 0.0_k1
   s3a3(10) = 0.0_k1

! a54 - a134
   s3a4(1) = 25.0_k1/16.0_k1
   s3a4(2) = 1.0_k1/4.0_k1
   s3a4(3) = 125.0_k1/108.0_k1
   s3a4(4) = 0.0_k1
   s3a4(5) = -53.0_k1/6.0_k1
   s3a4(6) = 23.0_k1/108.0_k1
   s3a4(7) = -341.0_k1/164.0_k1 
   s3a4(8) = 0.0_k1
   s3a4(9) = -341.0_k1/164.0_k1

!a65-- a135
   s3a5(1) = 1.0_k1/5.0_k1
   s3a5(2) = -65.0_k1/27.0_k1
   s3a5(3) = 61.0_k1/225.0_k1
   s3a5(4) = 704.0_k1/45.0_k1
   s3a5(5) = -976.0_k1/135.0_k1
   s3a5(6)= 4496.0_k1/1025.0_k1
   s3a5(7) = 0.0_k1
   s3a5(8) = 4496.0_k1/1025.0_k1

!a76-a136
   s3a6(1) = 125.0_k1/54.0_k1
   s3a6(2) = -2.0_k1/9.0_k1
   s3a6(3) = -107.0_k1/9.0_k1
   s3a6(4) = 311.0_k1/54.0_k1
   s3a6(5) = -301.0_k1/82.0_k1
   s3a6(6) = -6.0_k1/41.0_k1
   s3a6(7) = -289.0_k1/82.0_k1

!a87--a137
   s3a7(1) = 13.0_k1/900.0_k1
   s3a7(2) = 67.0_k1/90.0_k1
   s3a7(3) = -19.0_k1/60.0_k1
   s3a7(4) = 2133.0_k1/4100.0_k1
   s3a7(5) = -3.0_k1/205.0_k1
   s3a7(6) = 2193.0_k1/4100.0_k1

!a98--a138

  s3a8(1) = 3.0_k1
  s3a8(2) = 17.0_k1/6.0_k1
  s3a8(3) = 45.0_k1/82.0_k1
  s3a8(4) = -3.0_k1/41.0_k1
  s3a8(5) = 51.0_k1/82.0_k1
!a109-a139 
  s3a9(1) = -1.0_k1/12.0_k1  
  s3a9(2) = 45.0_k1/164.0_k1
  s3a9(3) = 3.0_k1/41.0_k1
  s3a9(4) = 33.0_k1/164.0_k1
!a1110-a1310
  s3a10(1) = 18.0_k1/41.0_k1
  s3a10(2) = 6.0_k1/41.0_k1
  s3a10(3) = 19.0_k1/41.0_k1
!a1211-a1311
  s3a11(1) = 0.0_k1
  s3a11(2) = 0.0_k1
! a1312 
  s3a12(1) = 1.0_k1



!b1,b6,b7,b8,b9,b10,b11
  s3b(1) = 41.0_k1/840_k1
  s3b(2) = 34.0_k1/105.0_k1
  s3b(3) = 9.0_k1/35.0_k1
  s3b(4) = 9.0_k1/35.0_k1 
  s3b(5) = 9.0_k1/280.0_k1
  s3b(6) = 9.0_k1/280.0_k1
  s3b(7) = 41.0_k1/840.0_k1 

! c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12
 
  s3c(1) = 2.0_k1/27.0_k1
  s3c(2) =1.0_k1/9.0_k1
  s3c(3) = 1.0_k1/6.0_k1
  s3c(4) = 5.0_k1/12.0_k1
  s3c(5) =0.5_k1
  s3c(6) =5.0_k1/6.0_k1
  s3c(7) =1.0_k1/6.0_k1
  s3c(8) = 2.0_k1/3.0_k1
  s3c(9) = 1.0_k1/3.0_k1
  s3c(10) =1.0_k1
  s3c(11) =0.0_k1




  return
 end subroutine initial_solver
	 


subroutine configure_solver
  use global
  use coffa
  use solver
  implicit none

   
   allocate(suma0(n_string_ab))
   allocate(sumphi0(numbasis,ras1%nptot))
   allocate(suma1(n_string_ab,6))
   allocate(sumphi1(numbasis,ras1%nptot,6))

   allocate(suma3(n_string_ab,14))
   allocate(sumphi3(numbasis,ras1%nptot,14))

  
  return
end subroutine configure_solver




       
     












 
