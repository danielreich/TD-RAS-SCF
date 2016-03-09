
module analysis_part
  use global
  use density
  use fedvr3d_basis_set
  use operator_spatial_space
  implicit none
contains

!
! calc. the expectation of the onebody operator
!==================================================================================
!
! calc. onebody energy part   < WF | onebody  | WF >
!
!             ______                     p
! = < WF |    \        { onebody }_pq * E       |WF>
!             /  p q                     q
!             ------ 
!
!=   _______
!    \
!     \      { rho 1 }_qp * { onebody }_pq  
!     / p,q  
!    /______
!
!                                         *
!  {onebody } _pq =  Integrate (   {fhi_p}   onebody {fhi_q}  )
!
!  fhi_p is the p-th spatial orbital
!=================================================================================

subroutine onebody_expectation(onebody_spatial,rdensity1,ndim,value_return)
implicit none
integer :: ip,iq
integer,intent(in) :: ndim
complex(kind=k2),intent(in) :: onebody_spatial(ndim*(ndim+1)/2)
complex(kind=k2),intent(in) :: rdensity1(ndim,ndim)
complex(kind=k2),intent(out) :: value_return
complex(kind=k2) :: ctemp


value_return = zzero

do ip = 1, ndim
 do iq =1, ndim

  if(ip>=iq) then
    ctemp = onebody_spatial( ia(ip) + iq )
  else
    ctemp = conjg(onebody_spatial( ia(iq) + ip))
  endif

  value_return = value_return + rdensity1(iq,ip)*ctemp

enddo
enddo

return
end subroutine onebody_expectation

!
! calc. twobody operator expectation 
!
! it can be improved by symmetry
!
subroutine twobody_expectation(twobody_spatial,rdensity2,ndim,value_return)
implicit none
integer :: ip,iq,is,ir
integer,intent(in) :: ndim
complex(kind=k2),intent(in) :: twobody_spatial(ndim,ndim,ndim,ndim),&
rdensity2(ndim,ndim,ndim,ndim)
complex(kind=k2),intent(out) :: value_return


value_return = zzero
do ip =1, ndim
   do iq = 1,ndim
      do ir=1,ndim
         do is=1,ndim
            
            value_return = value_return + twobody_spatial(ip,ir,is,iq)*rdensity2(ip,ir,is,iq)
            
         enddo
      enddo
   enddo
enddo

value_return= value_return*0.5d0


return
end subroutine twobody_expectation



!! Subroutines to perform the Fast Fourier Transform

SUBROUTINE four1(data,nn,isign)
!!Subroutine taken from Numerical recipes. It is a bit changed to adapt the code. Further information in pages 501-502 from Numerical recipes in Fortran 77 (Second edition)
  INTEGER isign,nn !! isign is the sign of the exponential
  REAL*8 data(1:2*nn) !! real array of length 2*nn. nn MUST be an integer power of 2
  INTEGER i,istep,j,m,mmax,n
  REAL*8 tempi,tempr
  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
  n=2*nn
  j=1
  do  i=1,n,2
     if(j.gt.i)then
        tempr=data(j)
        tempi=data(j+1)
        data(j)=data(i)
        data(j+1)=data(i+1)
        data(i)=tempr
        data(i+1)=tempi
     endif
     m=n/2
1    if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        goto 1
     endif
     j=j+m
  End do
     mmax=2
2    if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do  m=1,mmax,2
           do  i=m,n,istep
              j=i+mmax
              tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
              tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
              data(j)=data(i)-tempr
              data(j+1)=data(i+1)-tempi
              data(i)=data(i)+tempr
              data(i+1)=data(i+1)+tempi
              continue
           End do
           wtemp=wr
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
           
        End do
        mmax=istep
        goto 2
     endif
     return
   END SUBROUTINE four1


!! Subroutine fastft() gives the Fast Fourier Transform for a complex function on an equidistance grid starting from t_0=0, and tstep=t_i+1-t_i.

Subroutine fastft(ft,tstep,ff,fgrid)
implicit none
!!INPUT
complex*16, intent(in) :: ft(:) !! function in time
real*8, intent(in) :: tstep !! time step
!!OUTPUT
complex*16, allocatable, intent(out) :: ff(:) !! fast fourier transform
real*8, allocatable, intent(out) :: fgrid(:) !! grid of the fast fourier transform
!!Auxiliar variables
integer :: ii, jj, kk
real*8 :: dist, fstep
real*8, allocatable :: ftaux(:) !! Auxiliar array

If (allocated(ff)) deallocate(ff)

!! To check that the number of points in grid ft is a power of two.

dist=log(dble(size(ft)))/log(2.0d0)

If (abs(dble(int(dist))-dist).gt.1d-7) then
   write(*,*) 'ERROR IN fastft'
   write(*,*) 'The number of temporal time step ',size(ft), ' is not a power of 2'
stop
End If

allocate(ff(1:size(ft)+1)) !! Allocate the ouput
allocate(fgrid(1:size(ft)+1)) !! Allocate the ouput

ff=cmplx(0.0d0,0.0d0) !! initialize the output

!!initialize the workarray

allocate(ftaux(1:2*size(ft)))
ftaux=0.0d0

!!rearrange the input in the workarray

Do ii=1,size(ft)
   ftaux(2*ii-1)=real(ft(ii))
   ftaux(2*ii)=aimag(ft(ii))
End Do

!!Fast Fourier Transform

call four1(ftaux,size(ft),-1)

!!frequency step

fstep=1.0d0/(size(ft)*tstep)
fgrid=0.0d0

!!Rearrange to the output (see pages 501-502 in Numerical recipes in Fortran)

Do ii=size(ft)+2,2*size(ft),2
   jj=(ii-size(ft))/2
   ff(jj)=cmplx(ftaux(ii-1),ftaux(ii))
   fgrid(jj)=-1.0d0/(2.0d0*tstep)+dble(jj-1)*fstep
End Do

Do ii=2,size(ft),2
   jj=(ii+size(ft))/2
   ff(jj)=cmplx(ftaux(ii-1),ftaux(ii))
   fgrid(jj)=-1.0d0/(2.0d0*tstep)+dble(jj-1)*fstep
End Do

ff(size(ft)+1)=cmplx(ftaux(size(ft)+1),ftaux(size(ft)+2))
fgrid(size(ft)+1)=1.0d0/(2.0d0*tstep)

!!include normalization factor

ff=ff*tstep/sqrt(2.0d0*acos(-1.0d0)) !! normalize the wavefunction
fgrid=fgrid*2.0d0*acos(-1.0d0) !! change to angular frequency

!! Delete auxiliary arrays

deallocate(ftaux)

End Subroutine fastft

!!==========================================================================
!! SUBROUTINES TO CALCULATE THE WAVEFUNCTION IN THE MOMENTUM SPACE
!!==========================================================================

!! To calculate the wavefunction in momentum space. phi is the "wavefunction times r"
!! It calculates the transformation for a single point in k radial, for all the Spherical Harmonics.
!! It is written in FE-DVR basis+Spherical Harmonics

 subroutine phi_momentum(phix,rout,l,m,xglobal, weights,element,basis,fedvr_k_global,weights_k,phik)
  implicit none

!! This functions calculates the function in momentum space for a radial linear momentum. 
!! IMPORTANT: The input function is a radial function times a Spherical Harmonics
!! INPUT
 complex(kind=k2),intent(in) :: phix(:,:)
!! phix is the radial function times the spherical harmonic
real(kind=k1), intent(in) :: rout 
!! rout is the outer region 
 integer, intent(in), allocatable :: l(:),m(:)
!! L and M for the angular function in the argument.
 real*8, intent(in), allocatable :: xglobal(:)
!! xglobal is the radial basis fedvr
 real*8, intent(in), allocatable :: weights(:,:)
!! weights are the weights for the element in 1st argument and basis
!! in the element for the 2nd argument
 integer, intent(in) :: element(:), basis(:) 
!! element and basis in this element of the function chosen
  real(kind=k1),allocatable, intent(in) :: fedvr_k_global(:) 
!! all the nodes of the FE-DVR
  real(kind=k1),allocatable, intent(in) :: weights_k(:,:)
!! OUTPUT
complex(kind=k2), allocatable, intent(inout) :: phik(:,:)

!! AUXILIAR VARIABLES
integer :: i, j, ij ! auxiliar variable for loops
integer :: iorbital ! auxiliar variable running for orbitals
integer :: nr !! number of radial points in position space.
integer :: nk !! number of radial points in momentum space.
integer :: nangle !! number of angular functions
integer :: norbitals !! number of orbitals
integer :: laux !! auxiliar variable for the total angular momentum

INTERFACE

!! To obtain the spherical bessel

Function besselj(n,x)
  implicit none
  integer :: n !! order of the Spherical Bessel Function
  real*8 :: x !! argument of the Spherical Bessel Function
  real*8 :: besselj !! Spherical Bessel Function  
End Function besselj

END INTERFACE

nr=size(xglobal) !! number of radial points in position space.
nk=size(fedvr_k_global) !! number of radial points in momentum space.
nangle=size(phix(:,1))/size(xglobal) !! number of angular functions.
norbitals=size(phix(1,:)) !! number of orbitals

!! initialize the wavefunction in momentum space

If (allocated(phik)) then
   deallocate(phik)
Else
   continue
End If

allocate(phik(1:nangle*nk,1:size(phi(1,:)))) !! allocate memory for the function.
!#1st argument is the coordinate and the #2nd argument is the orbital.

phik=zzero !! initialize the momentum space function

Do iorbital=1,norbitals !! running in the orbitals
   
   Do i=1,nangle !! Run in the spherical harmonics
      
      laux=l(i) !! set the variable for the angular momentum.
      
      Do j=1,nk !! Run in the radial coordinates of the momentum space.
         
         Do ij=nr,1,-1 !! run in the radial part (a FE-DVR). 
            
            If (xglobal(ij).lt.rout) exit !! If it is in the inner region, do not take into account

            !! Take into account the weights on the points. 
            
            If (ij.eq.nr) then !! If this is the last point of the grid.

               phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)+phix((i-1)*nr+ij,iorbital)*besselj(laux, fedvr_k_global(j)*xglobal(ij))*sqrt(weights(element(ij),basis(ij)))*xglobal(ij)

            else !! this is not the last point of the grid.
               
               If (basis(ij).lt.basis(ij+1)) then !! basis(ij) is not the last node in an element
                  phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)+phix((i-1)*nr+ij,iorbital)*besselj(laux, fedvr_k_global(j)*xglobal(ij))*sqrt(weights(element(ij),basis(ij)))*xglobal(ij) !! We include the factor xglobal(ij) to have the reduced radial wavefuction
                  
               Else !! basis(ij) is the last node of the element.
                  phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)+phix((i-1)*nr+ij,iorbital)*besselj(laux, fedvr_k_global(j)*xglobal(ij))*sqrt((weights(element(ij),basis(ij))+weights(element(ij+1),1)))*xglobal(ij) !! We include the factor xglobal(ij) to have the reduced radial wavefuction
               Endif

            End If
         End Do
         
         phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt(2.0/acos(-1.0d0))*(-ci)**dble(laux)*fedvr_k_global(j) !! including the factor multiplying each element
         
         !! Multiply by the weights to obtain the FE-DVR in the momentum space
         !! Remember, the number of weights and nodes is the same than in the case of position space.

         If (j.eq.nk) then !! last point of the grid
            phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt(weights_k(element(j),basis(j))) !! not DVR functions
         else !! not the last point of the grid
            If (basis(j).lt.basis(j+1)) then !! basis(j) is not the last node in an element
               phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt(weights_k(element(j),basis(j))) !! not DVR functions
            else
               phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt((weights_k(element(j),basis(j))+weights_k(element(j+1),1))) !! bridge function
            End If
         End If
      End Do
   End Do
End Do

!! in phik(#1,#2) we have stored the radial reduced function in momentum space. That is, "the wavefunction times k". #1 is the FE-DVR+spherical harmonics and #2 is the orbital number.

return

End subroutine phi_momentum

!! To calculate the wavefunction in momentum space. phi is the "wavefunction times r"
!! It calculates the transformation for a single point in k radial, for all the Spherical Harmonics.
!! It is written in FE-DVR basis+Spherical Harmonics

!! Here we include a Hamming Window h(r), which multiplies the wavefunction to get rid of the tail of the wavefunction after r_out. It is defined as

!!       {   1-cos(pi/2*(r-r_out)/Delta), r_out <= r <= r_out+Delta
!! h(r)= {
!!       {  1, r_out+Delta< r

!! where Delta is the parameter of the Hamming window

 subroutine phi_momentum_hamming(phix,rout,rend,delta,l,m,xglobal, weights,element,basis,fedvr_k_global,weights_k,phik)
  implicit none

!! This functions calculates the function in momentum space for a radial linear momentum. 
!! IMPORTANT: The input function is a radial function times a Spherical Harmonics
!!INPUT
 complex(kind=k2),intent(in) :: phix(:,:)
!! phix is the radial function times the spherical harmonic
real(kind=k1), intent(in) :: rout 
!! rout is the outer region 
real(kind=k1), intent(in) :: rend
!! rend is the end of the box
real(kind=k1), intent(in) :: delta 
!! delta is the parameter of the Hamming window
 integer, intent(in), allocatable :: l(:),m(:)
!! L and M for the angular function in the argument.
 real*8, intent(in), allocatable :: xglobal(:)
!! xglobal is the radial basis fedvr
 real*8, intent(in), allocatable :: weights(:,:)
!! weights are the weights for the element in 1st argument and basis
!! in the element for the 2nd argument
 integer, intent(in) :: element(:), basis(:) 
!! element and basis in this element of the function chosen
  real(kind=k1),allocatable, intent(in) :: fedvr_k_global(:) 
!! all the nodes of the FE-DVR
  real(kind=k1),allocatable, intent(in) :: weights_k(:,:)
!! OUTPUT
complex(kind=k2), allocatable, intent(inout) :: phik(:,:)

!! AUXILIAR VARIABLES
integer :: i, j, ij ! auxiliar variable for loops
integer :: iorbital ! auxiliar variable running for orbitals
integer :: nr !! number of radial points in position space.
integer :: nk !! number of radial points in momentum space.
integer :: nangle !! number of angular functions
integer :: norbitals !! number of orbitals
integer :: laux !! auxiliar variable for the total angular momentum
 complex(kind=k2), allocatable :: phix_hamming(:,:)
!! phix_hamming is the radial function times the spherical harmonic times
!! the Hamming Window, h(r)

INTERFACE

!! To obtain the spherical bessel

Function besselj(n,x)
  implicit none
  integer :: n !! order of the Spherical Bessel Function
  real*8 :: x !! argument of the Spherical Bessel Function
  real*8 :: besselj !! Spherical Bessel Function  
End Function besselj

END INTERFACE

nr=size(xglobal) !! number of radial points in position space.
nk=size(fedvr_k_global) !! number of radial points in momentum space.
nangle=size(phix(:,1))/size(xglobal) !! number of angular functions.
norbitals=size(phix(1,:)) !! number of orbitals

!! initialize the wavefunction in position space

allocate(phix_hamming(1:nr*nangle,1:norbitals))
phix_hamming=zzero

!! Now, we multiply the orbitals times the Hamming Window, h(r)

Do i=1,norbitals !! Running in the orbitals
   Do j=1, nangle !! Running in the angular part
      Do ij=1,nr !! Running in the radial part
         
         If (xglobal(ij).gt.rout.and.xglobal(ij).lt.rout+delta) then
            phix_hamming((j-1)*nr+ij,i)=(1.0d0-cos(acos(-1.0d0)*(xglobal(ij)-rout)/(2.0d0*delta)))*phix((j-1)*nr+ij,i)
         Else If (xglobal(ij).gt.rend-delta.and.xglobal(ij).lt.rend) then
            phix_hamming((j-1)*nr+ij,i)=cos(acos(-1.0d0)*(xglobal(ij)-rend+delta)/(2.0d0*delta))*phix((j-1)*nr+ij,i)
         Else
            phix_hamming((j-1)*nr+ij,i)=phix((j-1)*nr+ij,i)
         End If
         
      End Do
   End Do
End Do
!! initialize the wavefunction in momentum space

If (allocated(phik)) then
   deallocate(phik)
Else
   continue
End If

allocate(phik(1:nangle*nk,1:size(phi(1,:)))) !! allocate memory for the function.
!#1st argument is the coordinate and the #2nd argument is the orbital.

phik=zzero !! initialize the momentum space function

Do iorbital=1,norbitals !! running in the orbitals
   
   Do i=1,nangle !! Run in the spherical harmonics
      
      laux=l(i) !! set the variable for the angular momentum.
      
      Do j=1,nk !! Run in the radial coordinates of the momentum space.
         
         Do ij=nr,1,-1 !! run in the radial part (a FE-DVR). 
            
            If (xglobal(ij).lt.rout) exit !! If it is in the inner region, do not take into account

            !! Take into account the weights on the points. 
            
            If (ij.eq.nr) then !! If this is the last point of the grid.

               phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)+phix_hamming((i-1)*nr+ij,iorbital)*besselj(laux, fedvr_k_global(j)*xglobal(ij))*sqrt(weights(element(ij),basis(ij)))*xglobal(ij)

            else !! this is not the last point of the grid.
               
               If (basis(ij).lt.basis(ij+1)) then !! basis(ij) is not the last node in an element
                  phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)+phix_hamming((i-1)*nr+ij,iorbital)*besselj(laux, fedvr_k_global(j)*xglobal(ij))*sqrt(weights(element(ij),basis(ij)))*xglobal(ij) !! We include the factor xglobal(ij) to have the reduced radial wavefuction
                  
               Else !! basis(ij) is the last node of the element.
                  phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)+phix_hamming((i-1)*nr+ij,iorbital)*besselj(laux, fedvr_k_global(j)*xglobal(ij))*sqrt((weights(element(ij),basis(ij))+weights(element(ij+1),1)))*xglobal(ij) !! We include the factor xglobal(ij) to have the reduced radial wavefuction
               Endif

            End If
         End Do
         
         phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt(2.0/acos(-1.0d0))*(-ci)**dble(laux)*fedvr_k_global(j) !! including the factor multiplying each element
         
         !! Multiply by the weights to obtain the FE-DVR in the momentum space
         !! Remember, the number of weights and nodes is the same than in the case of position space.

         If (j.eq.nk) then !! last point of the grid
            phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt(weights_k(element(j),basis(j))) !! not DVR functions
         else !! not the last point of the grid
            If (basis(j).lt.basis(j+1)) then !! basis(j) is not the last node in an element
               phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt(weights_k(element(j),basis(j))) !! not DVR functions
            else
               phik((i-1)*nk+j,iorbital)=phik((i-1)*nk+j,iorbital)*sqrt((weights_k(element(j),basis(j))+weights_k(element(j+1),1))) !! bridge function
            End If
         End If
      End Do
   End Do
End Do

!! in phik(#1,#2) we have stored the radial reduced function in momentum space. That is, "the wavefunction times p". #1 is the FE-DVR+spherical harmonics and #2 is the orbital number.

deallocate(phix_hamming)

return

End subroutine phi_momentum_hamming



!! To plot the monoparticular orbitals and the radial density in momentum space in momentum space.

 subroutine plot_phik_fedvr3d(phiuni,xglobal,nb_angle,element,basis,l,m, weights,rho)
  implicit none

! This subroutine produces the density for the monoparticular
! functions in momentum space.
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
character*7 :: file_orbital
real*8 :: exp_k2,exp_k,exp_km1 !! <k**2>, <k>, <k**-1>
real*8 :: exp_km2 
real*8 :: exp_l2, exp_m2, exp_m, norm

allocate(density(size(xglobal)))

Do i=1,size(phiuni(1,:)) !! Run in the orbitals

   !Initialize the varibles
   density=0.0d0 !!Density
   exp_k2=0.0d0 !! <k2>
   exp_k=0.0d0 !! <k>
   exp_km1=0.0d0 !! <k**-1>
   exp_km2=0.0d0 !! <k**-2>
   exp_l2=0.0d0 !! <L**2>
   exp_m2=0.0d0 !! <M**2>
   exp_m=0.0d0 !! <M>

   counter=0
   write(file_orbital,'(a5,i2.2)') 'phik_', i

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
         exp_k2=exp_k2+density(k)*xglobal(k)**2.0d0
         exp_k=exp_k+density(k)*xglobal(k)
         exp_km1=exp_km1+density(k)*xglobal(k)**(-1.0d0)
         exp_km2=exp_km2+density(k)*xglobal(k)**(-2.0d0)
      End Do
      
      write(122+i,'(a29,i2.2,a3)') ' # This is the density |phik_',i,'|^2'
      write(122+i,*) '# <k**2>=', exp_k2
      write(122+i,*) '# <k>=', exp_k
      write(122+i,*) '# <k**-1>=', exp_km1
      write(122+i,*) '# <k**-2>=', exp_km2
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

exp_k2=0.0d0
exp_k=0.0d0
exp_km1=0.0d0
exp_km2=0.0d0
norm=0.0d0
!! Calculation of the expectation values

Do k=1,size(density)
   norm=norm+density(k) !! Norm of the density
   exp_k2=exp_k2+density(k)*xglobal(k)**2.0d0
   exp_k=exp_k+density(k)*xglobal(k)
   exp_km1=exp_km1+density(k)*xglobal(k)**(-1.0d0)
   exp_km2=exp_km2+density(k)*xglobal(k)**(-2.0d0)
End Do

exp_k2=exp_k2/norm
exp_k=exp_k/norm
exp_km1=exp_km1/norm
exp_km2=exp_km2/norm

!! Opening the density file in momentum space

Open(122,file='rho_k')

!! Writting the heading of the file

write(122,*) '# <p**2>=', exp_k2
write(122,*) '# <p>=', exp_k
write(122,*) '# <p**-1>=', exp_km1
write(122,*) '# <p**-2>=', exp_km2
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

End subroutine plot_phik_fedvr3d

!!====================================================================
!!
!! CALCULATION OF THE DENSITY FUNCTION INCLUDING THE ANGLES.
!!
!! The density function is written in FE-DVR + SPHERICAL HARMONICS,
!! as the orbitals, but including more spherical harmonics.
!! The coefficients include the weights of the FE-DVR, that is, to 
!! perform the integrals we do not need to include them again.
!! This subroutines can be used in both, position and momentum space.
!!
!!====================================================================


!!--------------------------------------------------------------------
!! The subroutine basis_set_rho produces the appropiate basis set in 
!! both FE-DVR + SPHERICAL COORDINATES
!!--------------------------------------------------------------------

Subroutine basis_set_rho3d(element,basis,lmax,mmax,label)
  implicit none
!! INPUT
 integer, allocatable,intent(in) :: element(:), basis(:) 
!! element and basis in this element of the function chosen
 integer, intent(in) :: lmax, mmax
!! maximum value of single particle angular momentum and magnetic 
!! quantum number
!! OUTPUT
integer, allocatable,intent(out) :: label(:,:)
!! label(#1,1)-> FE-DVR function of #1
!! label(#1,2)-> angular momentum of #1
!! label(#1,3)-> magnetic quantum number of #1

!! AUXILIAR
integer :: i, j, ij, l, m
integer :: counter

!! First, how many elements are there in the new basis set?

counter=0

Do l=0,2*lmax !! Running in the total angular momentum
   Do m=max(-l,-2*mmax),min(l,2*mmax)
      Do i=1,size(element) !! FE-DVR basis set
         counter=counter+1
      End Do
   End Do
End Do

!! Now, we can allocate the label(:,:)

allocate(label(1:counter,1:3))

label=0 !! Initialize label

!! Now, we assign the labels

counter=0

Do l=0,2*lmax !! Running in the total angular momentum
   Do m=max(-l,-2*mmax),min(l,2*mmax) !! Running magnetic quantum number
      Do i=1,size(element) !! FE-DVR basis set
         counter=counter+1
         label(counter,1)=i
         label(counter,2)=l
         label(counter,3)=m
      End Do
   End Do
End Do

!! label assigned 
return

End Subroutine basis_set_rho3d

!!---------------------------------------------------------------------------
!! This subroutine calculates the 3D one density function. ***CHECK***
!!---------------------------------------------------------------------------
!! 
!! Find details of the structure of the outcoming result just before 
!! the subroutine basis_set_rho3d.
!!
!!---------------------------------------------------------------------------

Subroutine rho3d(phiuni,xglobal,label_phi, label_rho3d,rho,n_phi, n_all)
implicit none
!! This subroutine calculates the 3D monoparticular density and the 3D density
!! for all the electrons
!!INPUT
 complex(kind=k2),allocatable, intent(in) :: phiuni(:,:)
!! phi in terms of the fedvr global basis.
real (kind=k1), allocatable,intent(in) :: xglobal(:)
!! fedvr in the radial basis
 integer, allocatable, intent(in) :: label_phi(:,:) !! refered to the wavefunction
!! label_phi(#1,1)-> FE-DVR function of #1
!! label_phi(#1,2)-> angular momentum of #1
!! label_phi(#1,3)-> magnetic quantum number of #1
 integer, allocatable, intent(in) :: label_rho3d(:,:) !! refered to the density
!! label_rho3d(#1,1)-> FE-DVR function of #1
!! label_rho3d(#1,2)-> angular momentum of #1
!! label_rho3d(#1,3)-> magnetic quantum number of #1
 complex(8), allocatable, intent(in) :: rho(:,:)
!! First arguments is the annihilation and the second creation
!! wavefunction

!!OUTPUT
 complex(kind=k2), allocatable, intent(out) :: n_phi(:,:) !! monoparticular 3D density for each orbitals
 complex(kind=k2), allocatable, intent(out) :: n_all(:) !! 3D density for all the electrons

!!Auxiliary quantities 
integer :: norbital !! number of orbitals.
integer :: nbasis_phi !! number of basis set functions for phi.
integer :: nbasis_rho3d !! number of basis set functions for the 3D density.
integer :: iorbital, jorbital !! Running in the orbitals
integer :: ibasis_phi, jbasis_phi, nbradial !! running in the basis of the wavefunction.
integer :: ibasis_angle, jbasis_angle, nbangle
integer :: ibasis_rho3d, jbasis_rho3d !! counters for spatial basis for the density.
integer :: lrho3d, l,lp
integer :: mrho3d, m,mp

!!Functions

INTERFACE

!! FUNCTION TO OBTAIN THE 3J SYMBOLS

!!--------------------------------------------------------
!! The function wigner3j calculates the 3J wigner symbols 
!!   _________________
!!  /                  \
!! |   a     b     c    |
!! | alfa  beta  gamma  |
!!  \ _________________/
!!  
!!--------------------------------------------------------

  FUNCTION wigner3j(a,b,c,alfa,beta,gamma)
    IMPLICIT NONE
    Integer :: a, b, c
    Integer :: alfa, beta, gamma
    Real*8 :: wigner3j
  END FUNCTION wigner3j

!! FUNCTION TO OBTAIN THE sign_m1=(-1)**m

  FUNCTION sign_m1(m) 
    IMPLICIT NONE
    Integer, INTENT(IN)  :: m
    REAL*8              :: sign_m1
  END FUNCTION sign_m1

END INTERFACE

print*,
print*, 'Skip subroutine rho3d'
print*, 'in module_analysis.f90'
print*, 'It takes too long!!'
return
print*,
!return

norbital= size(phiuni(1,:)) !! number of orbitals
nbasis_phi=size(label_phi(:,1)) !! number of basis set functions for wavefuctions
nbasis_rho3d=size(label_rho3d(:,1)) !! number of basis set functions for density

nbangle=nbasis_phi/size(xglobal) !! number of spherical harmonics for the wavefunctions
nbradial=size(xglobal) !! number of the radial part of the FE-DVR

!! First, we allocate the array for the density of the orbitals

allocate(n_phi(1:nbasis_rho3d,1:norbital))
allocate(n_all(1:nbasis_rho3d))

!! We initialize the output

n_phi=zzero
n_all=zzero

!! Now we obtain the density functions for the orbitals.
      Do ibasis_rho3d=1,nbasis_rho3d !! Running in the spatial basis for the density.
         
         lrho3d=label_rho3d(ibasis_rho3d,2)
         mrho3d=label_rho3d(ibasis_rho3d,3)
         
         Do ibasis_angle=1,nbangle !! Run in angle for the wavefunction, since 
            !! the FE-DVR is the same
            
            ibasis_phi=(ibasis_angle-1)*nbradial+label_rho3d(ibasis_rho3d,1)
            !! number of basis set for the wave function of the bra.
            
            l=label_phi(ibasis_phi,2)
            m=label_phi(ibasis_phi,3)
            
            Do jbasis_angle=1,nbangle
               
               jbasis_phi=(jbasis_angle-1)*nbradial+label_rho3d(ibasis_rho3d,1)
               !! number of basis set for the wave function of the ket.         
               
               lp=label_phi(jbasis_phi,2)
               mp=label_phi(jbasis_phi,3)
               
               !! Check the allowed couplings
               
               If (mrho3d.ne.mp-m) cycle
               !! If M\ne m-m' then cycle
               If (abs(lrho3d-lp).gt.l) cycle
               !! If |L-l'|>l cycle
               If (lrho3d+lp.lt.l) cycle
               !! If L+l'<l cycle
                             
               !! Contribution for each 

               Do iorbital=1, size(phiuni(1,:)) !! run or the orbitals
   
                  Do jorbital=iorbital, size(phiuni(1,:)) !! run or the orbitals          
               
                     If (iorbital.eq.jorbital) then !! for the contributions coming from a_iorbital^+a_iorbital

               !! Expression of the density for each state.
               
                        n_phi(ibasis_rho3d,iorbital)=n_phi(ibasis_rho3d,iorbital)+ conjg(phiuni(ibasis_phi,iorbital))*phiuni(jbasis_phi,iorbital)*sign_m1(m+mrho3d)*wigner3j(lrho3d,l,lp,0,0,0)*wigner3j(lrho3d,l,lp,-mrho3d,-m,mp)*sqrt(dble((2*lrho3d+1)*(2*l+1)*(2*lp+1))/(4.0d0*acos(-1.0d0)))


                        n_all(ibasis_rho3d)=n_all(ibasis_rho3d) + conjg(phiuni(ibasis_phi,iorbital))*phiuni(jbasis_phi,jorbital)*sign_m1(m+mrho3d)*wigner3j(lrho3d,l,lp,0,0,0)*wigner3j(lrho3d,l,lp,-mrho3d,-m,mp)*sqrt(dble((2*lrho3d+1)*(2*l+1)*(2*lp+1))/(4.0d0*acos(-1.0d0)))*rho(iorbital,jorbital)
                  
                     Else !! for contributions coming from a_iorbital^+a_jorbital

                        n_all(ibasis_rho3d)=n_all(ibasis_rho3d) + 2.0d0*real(conjg(phiuni(ibasis_phi,iorbital))*phiuni(jbasis_phi,jorbital)*rho(iorbital,jorbital))*sign_m1(m+mrho3d)*wigner3j(lrho3d,l,lp,0,0,0)*wigner3j(lrho3d,l,lp,-mrho3d,-m,mp)*sqrt(dble((2*lrho3d+1)*(2*l+1)*(2*lp+1))/(4.0d0*acos(-1.0d0)))
                        
                     End If

                  End Do
            
               End Do
               
            End Do
            
         End Do
         
      End Do

return

End Subroutine rho3d


!! Using the coupled representations subroutines

Subroutine rho3d_coupled(orb, rho,n_phi, n_all)
implicit none

!!INPUT

complex (kind=k2), allocatable, intent(in) :: orb(:,:) !! orbitals in the uncoupled basis

 complex(8), allocatable, intent(in) :: rho(:,:)
!! First arguments is the annihilation and the second creation
!! wavefunction

!!OUTPUT

 complex(kind=k2), allocatable, intent(out) :: n_phi(:,:) !! monoparticular 3D density for each orbitals

 complex(kind=k2), allocatable, intent(out) :: n_all(:) !! 3D density for all the electrons

!!AUXILIAR VARIBALES
complex(kind=k2), allocatable :: orb_coupled(:,:,:)
!! conjg(orb(rl'm',#1))orb(rl'm',#2) -> orb_coupled(rlm,#1,#2)
integer :: i_orb, j_orb
integer :: nglobal, norbital
type (llmm_lm_type) :: llmm_lm_aux !! coefficients which link the uncouple basis with the coupled basis for Spherical Harmonics. 

!! Build the coefficients to link the uncoupled basis with the 
!! coupled basis with maximum accesible angular momentum.

call uncoupled_to_coupled(fedvr3d%l_max,fedvr3d%m_max,2*fedvr3d%l_max,llmm_lm_aux)

!! Calculate conjg(orb(rl'm',#1))orb(rl'm',#2)

call phikphil_to_phikl(orb, llmm_lm_aux, orb_coupled) !! we obtain the 
!! orbitals in the coupled basis 
!! conjg(orb(rl'm',#1))orb(rl'm',#1)->orb_coupled(rlm,#1,#2)

nglobal=size(orb_coupled(:,1,1)) !! global basis 
norbital=size(orb(1,:)) !! number of orbitals

!! allocate the density associated to the orbitals

allocate(n_phi(1:nglobal,1:norbital)) 
n_phi=zzero !! initialize the density associated to the orbitals

!! allocate the total density

allocate(n_all(1:nglobal))
n_all=zzero !! initialize the total density

Do i_orb=1, norbital !! loop in the orbitals
   n_phi(:,i_orb)=orb_coupled(:,i_orb,i_orb) !! store in n_phi
   Do j_orb=1, norbital !! loop in the orbitals
      n_all(:)=n_all(:)+rho(j_orb,i_orb)*orb_coupled(:,j_orb,i_orb)
   End Do
End Do

deallocate(orb_coupled)

End Subroutine rho3d_coupled

!!====================================================================
!!
!! END OF CALCULATION OF THE DENSITY FUNCTION INCLUDING THE ANGLES.
!!
!!====================================================================

!!====================================================================
!!
!! EVALUATION OF THE DENSITY FUNCTION INCLUDING THE ANGLES.
!!
!!====================================================================

!! Subroutine rho3d_evaluation gives the density for a fixed angle
!! \theta_k, fixing phi=0

subroutine rho3d_evaluation(rho3d,rho3d_hamming,global_to_local,fedvr_k,diff,diff_hamming)
implicit none
!! INPUT
complex*16, allocatable, intent(in) :: rho3d(:),rho3d_hamming(:)
!! 3D density without and with Hamming Window
integer, allocatable, intent(in) :: global_to_local(:,:)
!! global_to_local(#1,1)-> FE-DVR function of #1
!! global_to_local(#1,2)-> angular momentum of #1
!! global_to_local(#1,3)-> magnetic quantum number of #1 
real(kind=k1),allocatable, intent(in) :: fedvr_k(:) 
!! all the nodes of the FE-DVR in radial momentum
!! OUTPUT
complex(kind=k2), allocatable, intent(out) :: diff(:,:), diff_hamming(:,:)
!! Differential probability distribution without and with 
!! Hamming Window
!! 
!! Arguments
!!
!! #1: FEDVR node
!! #2: angle. Cases: 1 => 0, 2 => pi/4, 3 => pi/2, 4 => 3pi/2, 5 => pi
!! 
!! AUXILIAR VARIABLES
integer :: ntotal !! number of FE-DVR + Spherical Harmonics in coupled basis
integer :: nr !! number of FE-DVR
integer :: i, j, k !! index in the loops
integer :: lprime, mprime !! orbital angular momentum and magnetic quantum number in the loop
integer :: kprime !! node in radial momentum
real(kind=k2) :: pi
real(kind=k2) :: angle !! angle in the loop

Interface
    Function y_lm(l,m,theta,phi)
      !! Spherical Harmonic Y_{lm}(theta,phi)
      complex*16 y_lm
      integer l, m
      real*8 theta, phi
    End Function y_lm
End Interface

ntotal=size(rho3d) !! number of FE-DVR + Spherical Harmonics in coupled basis
nr=size(fedvr_k)
pi=acos(-1.0d0)

!! initialize output

allocate(diff(1:nr,1:5),diff_hamming(1:nr,1:5))

diff=zzero
diff_hamming=zzero

Do i=1, ntotal
   kprime=global_to_local(i,1) !! running in the radial momentum
   lprime=global_to_local(i,2) !! running in angular momentum
   mprime=global_to_local(i,3) !! running in magnetic quantum number
   Do j=1,5 !! Running in the angles for the evaluation
      angle=dble((j-1))/4.0d0*pi
      diff(kprime,j)=diff(kprime,j)+ y_lm(lprime,mprime,angle,0.0d0)*rho3d(i)
      diff_hamming(kprime,j)=diff_hamming(kprime,j)+ y_lm(lprime,mprime,angle,0.0d0)*rho3d_hamming(i)
   End Do  
EndDo

return
End subroutine rho3d_evaluation

!! Subroutine pxpz plane gives the density in the px-pz plane, assuming
!! py = 0.
!! The density is given in polar coordinates.

subroutine pxpz_density(rho3d,global_to_local,fedvr_k,iangle,pxpz)
implicit none
!! INPUT
complex*16, allocatable, intent(in) :: rho3d(:)
!! 3D density 
integer, allocatable, intent(in) :: global_to_local(:,:)
!! global_to_local(#1,1)-> FE-DVR function of #1
!! global_to_local(#1,2)-> angular momentum of #1
!! global_to_local(#1,3)-> magnetic quantum number of #1 
real(kind=k1),allocatable, intent(in) :: fedvr_k(:) 
!! all the nodes of the FE-DVR in radial momentum
integer, intent(in) :: iangle
!! Number of the point with the same radial momentum k_r

!! OUTPUT
complex(kind=k2), allocatable, intent(out) :: pxpz(:,:)
!! Density on the pxpz plane, with py=0
!! 
!! Arguments
!!
!! #1: k_r
!! #2: theta_2D !! If theta_2D=0 then kx=0, kz>0 
                !! for theta_2D=\pi/2 then kz=0, kx>0
                !! for theta_2D=\pi then kx=0, kz<0
                !! for theta_2D=3\pi/2 then kz=0, kx<0
!! 

!! AUXILIAR VARIABLES
integer :: i, j, k, jk !! index for the loops
integer :: nr, ntotal,n_pxpz
real(kind=k1) :: pi
integer :: kprime, lprime, mprime !! FE-DVR, L, and M the 3D monoparticular momentum density.
real(kind=k1) :: angle !! angle \theta_k 

Interface
    Function y_lm(l,m,theta,phi)
      !! Spherical Harmonic Y_{lm}(theta,phi)
      complex*16 y_lm
      integer l, m
      real*8 theta, phi
    End Function y_lm
End Interface

ntotal=size(rho3d) !! number of FE-DVR + Spherical Harmonics in coupled basis
nr=size(fedvr_k) !! number of FE-DVR 
n_pxpz=nr*iangle !! number of nodes per each side
pi=acos(-1.0d0)

!! initialize output

allocate(pxpz(1:nr,1:iangle)) !! allocate memory

pxpz=zzero !! initialize the density

Do i=1, ntotal
   kprime=global_to_local(i,1) !! running in the radial momentum
   lprime=global_to_local(i,2) !! running in angular momentum
   mprime=global_to_local(i,3) !! running in magnetic quantum number
   Do j=1,iangle !! Running in the angles for the evaluation

      !! the relation between angles.
      !! j runs in the polar angle (theta_polar), which is related to theta_k (angle) as
      
      !! theta_polar=arctan(kz/kx)=arctan(cos(theta_k)/sin(theta_k))=
      !! = arctan(sin(pi/2-theta_k)/cos(pi/2-theta_k))=pi/2-theta_k
      
      angle=pi*.5d0 - dble(j-1)/dble(iangle)*pi*2.0d0     
      pxpz(kprime,j)=pxpz(kprime,j)+ y_lm(lprime,mprime,angle,0.0d0)*rho3d(i)
   End Do  
EndDo

return
End subroutine pxpz_density



!!====================================================================
!!
!! END EVALUATION OF THE DENSITY FUNCTION INCLUDING THE ANGLES.
!!
!!====================================================================



!! subroutine plot_phi2_fedvr3d_coupled calculates and plot the two body radial density using the !! coupled basis.

subroutine plot_phi2_fedvr3d_coupled(orb,xglobal,element,basis, weights,rho2)
  implicit none
  ! This subroutine produces the two particle density distribution.
  ! rho(r_1,r_2)
  !!INPUT
  complex(kind=k2), allocatable,intent(in) :: orb(:,:)
  !! phi in terms of the fedvr global basis.
  real*8, intent(in), allocatable :: xglobal(:)
  !! xglobal is the radial basis fedvr
  real*8, intent(in), allocatable :: weights(:,:)
  !! weights are the weights for the element in 1st argument and basis
  !! in the element for the 2nd argument
  integer, intent(in) :: element(:), basis(:) 
  !! OUTPUT
  complex(kind=k2), allocatable:: density2(:,:)
  !! element and basis in this element of the function chosen
  complex(8), allocatable, intent(in) :: rho2(:,:,:,:)
  
  !! OUTPUT
  !! The output is written in the subroutine.

  !! AUXILIAR VARIABLES
  integer :: nglobal, norbital
  type (llmm_lm_type) :: llmm_lm_aux !! coefficients which link the uncouple basis with the coupled basis for Spherical Harmonics. 
  complex(kind=k2), allocatable:: orb_coupled(:,:,:)
  real(kind=k1) :: norm
  integer :: k, k2, p, q, r, s
  integer :: kk1, kk2 !! indexes in the loops to calculate ionization
  real(kind=k1) :: ionization_double, aux !! ionization and 
  real(kind=k1),allocatable :: ionization_double_aux(:)
  !! auxiliar variable to calculate it.

!! initialize the size of the arrays

nglobal=size(xglobal) !! global basis 
norbital=size(orb(1,:)) !! number of orbitals

!! Build the coefficients to link the uncoupled basis with the 
!! coupled basis with maximum accesible angular momentum.

call uncoupled_to_coupled(fedvr3d%l_max,fedvr3d%m_max,0,llmm_lm_aux) 
!! The only we only calculate coupling to l=zero. The other contributions are zero.

!! Calculate conjg(orb(rl'm',#1))orb(rl''m'',#2)

call phikphil_to_phikl(orb, llmm_lm_aux, orb_coupled) !! we obtain the 
!! orbitals in the coupled basis 
!! conjg(orb(rl'm',#1))orb(rl''m'',#2)->orb_coupled(rlm,#1,#2)
 
!! allocate the density2
allocate(density2(1:size(xglobal),1:size(xglobal)))
!! set to zzero
density2=zzero


!! Obtaining rho(r1,r2). At this point it is written in the DVR

Do p=1, size(rho2(:,1,1,1)) !! Loop in the orbitals a_p+
   Do r=1,size(rho2(1,:,1,1)) !! Loop in the orbitals a_r+
      Do s=1,size(rho2(1,1,:,1)) !! Loop in the orbitals a_s
         Do q=1, size(rho2(1,1,1,:)) !! Loop in the orbitals a_q          
            Do k2=1,size(density2(:,1)) !! run the number of points of r1, only with Y_00
               Do k=1,size(density2(1,:)) !! run the number of points of r2, only with Y_00
                  
                  density2(k,k2)=density2(k,k2)+orb_coupled(k2,r,s)*orb_coupled(k,p,q)*rho2(p,r,s,q) !!CHECK 02.06.2015

                  
               End Do !! k2
            End Do !! k
         End Do !! q
      End Do !! s
   End Do !! r
End Do !! p

!! Now, density2(k,k2) is a written as f(r_1,r_2)*Y_00(\Omega_1)*Y_00(\Omega_2)

!! We include a factor 4\pi, to take into account another Spherical Harmonic
!! for each component.

!! We divide by two because of the combinatorials of the electrons 

density2=density2*4.0d0*acos(-1.0d0)/2

!! Calculate the norm
norm=0.0d0 !! norm of the 2D radial density
Do k2=1,size(density2(1,:))
   Do k=1,size(density2(:,1))
      norm=norm+real(density2(k,k2))
   End Do
End Do

!! calculating ionization for double ionization P(r_1,r_2), that is, 
!! the probability of finding an electron further than r_1 and the 
!! other further thatn r_2.

!! P(r_1,r_2)=\int r_1^\infty (\int_r_2^\infty \rho(r_1,r_2) d r_2) dr_1=
!! = \int r_1 \infty \rho_{r_2}(r_1) dr_1.

!! First, we prepare the array to the variable to store \rho_{r_2}(r_1),
!! that is, ionization_double_aux

allocate(ionization_double_aux(nglobal))

ionization_double=zero !!
ionization_double_aux=zero !!

Open(510,file='ionization_double')

Do kk2=nglobal,1,-1
   ionization_double=zero
   Do kk1=nglobal,1,-1
   ionization_double_aux(kk1)=ionization_double_aux(kk1)+real(density2(kk1,kk2))
   ionization_double=ionization_double+ionization_double_aux(kk1)
   write(510,*) xglobal(kk1), xglobal(kk2),ionization_double
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
            write(122,*) xglobal(k),xglobal(k2), real(density2(k,k2))/(weights(element(k),basis(k))*weights(element(k2),basis(k2)))/norm
         Elseif (basis(k).gt.basis(k+1).and.basis(k2).lt.basis(k2+1)) then 
            !! basis(k) is the last node of the element, but not basis(k2).
            aux=((weights(element(k),basis(k))+weights(element(k+1),1))*weights(element(k2),basis(k2)))
            write(122,*) xglobal(k),xglobal(k2), real(density2(k,k2))/(aux*norm)
         Elseif (basis(k2).gt.basis(k2+1).and.basis(k).lt.basis(k+1)) then 
            !! basis(k2) is the last node of the element, but not basis(k).
            aux=((weights(element(k2),basis(k2))+weights(element(k2+1),1))*weights(element(k),basis(k)))
            write(122,*) xglobal(k),xglobal(k2), real(density2(k,k2))/(aux*norm)
         Elseif (basis(k2).gt.basis(k2+1).and.basis(k).gt.basis(k+1)) then 
            aux=((weights(element(k2),basis(k2))+weights(element(k2+1),1))*(weights(element(k),basis(k))+weights(element(k+1),1)))
            write(122,*) xglobal(k),xglobal(k2), real(density2(k,k2))/(aux*norm)
         Endif
      End Do
      write(122,*) ''
   End Do
   Close(122)


deallocate(ionization_double_aux,density2)
  return
End subroutine plot_phi2_fedvr3d_coupled

!! subroutine phi2_fedvr3d_coupled calculates the two body radial density using the coupled basis.

subroutine phi2_fedvr3d_coupled(orb,xglobal,rho2,density2)
  implicit none
  ! This subroutine produces the two particle density distribution.
  ! rho(r_1,r_2)
  !!INPUT
  complex(kind=k2), allocatable,intent(in) :: orb(:,:)
  !! phi in terms of the fedvr global basis.
  real*8, intent(in), allocatable :: xglobal(:)
  !! xglobal is the radial basis fedvr
  complex(8), allocatable, intent(in) :: rho2(:,:,:,:)
  !! density a_p^+ a_r^+ a_s a_q
  !! OUTPUT
  complex(kind=k2), allocatable,intent(out):: density2(:,:)
  !! 2 body density written in the FE-DVR
  
  !! OUTPUT
  !! The output is written in the subroutine.

  !! AUXILIAR VARIABLES
  integer :: nglobal, norbital
  type (llmm_lm_type) :: llmm_lm_aux !! coefficients which link the uncouple basis with the coupled basis for Spherical Harmonics. 
  complex(kind=k2), allocatable:: orb_coupled(:,:,:)
  real(kind=k1) :: norm
  integer :: k, k2, p, q, r, s
  integer :: kk1, kk2 !! indexes in the loops to calculate ionization
  real(kind=k1) :: ionization_double, aux !! ionization and 
  real(kind=k1),allocatable :: ionization_double_aux(:)
  !! auxiliar variable to calculate it.

!! initialize the size of the arrays

nglobal=size(xglobal) !! global basis 
norbital=size(orb(1,:)) !! number of orbitals

!! Build the coefficients to link the uncoupled basis with the 
!! coupled basis with maximum accesible angular momentum.

call uncoupled_to_coupled(fedvr3d%l_max,fedvr3d%m_max,0,llmm_lm_aux) 
!! The only we only calculate coupling to l=zero. The other contributions are zero.

!! Calculate conjg(orb(rl'm',#1))orb(rl''m'',#2)

call phikphil_to_phikl(orb, llmm_lm_aux, orb_coupled) !! we obtain the 
!! orbitals in the coupled basis 
!! conjg(orb(rl'm',#1))orb(rl''m'',#2)->orb_coupled(rlm,#1,#2)
 
!! allocate the density2
allocate(density2(1:size(xglobal),1:size(xglobal)))
!! set to zzero
density2=zzero


!! Obtaining rho(r1,r2). At this point it is written in the DVR

Do p=1, size(rho2(:,1,1,1)) !! Loop in the orbitals a_p+
   Do r=1,size(rho2(1,:,1,1)) !! Loop in the orbitals a_r+
      Do s=1,size(rho2(1,1,:,1)) !! Loop in the orbitals a_s
         Do q=1, size(rho2(1,1,1,:)) !! Loop in the orbitals a_q          
            Do k2=1,size(density2(:,1)) !! run the number of points of r1, only with Y_00
               Do k=1,size(density2(1,:)) !! run the number of points of r2, only with Y_00
                  
                  density2(k,k2)=density2(k,k2)+orb_coupled(k2,r,s)*orb_coupled(k,p,q)*rho2(p,r,s,q) !!CHECK 02.06.2015

                  
               End Do !! k2
            End Do !! k
         End Do !! q
      End Do !! s
   End Do !! r
End Do !! p

!! Now, density2(k,k2) is a written as f(r_1,r_2)*Y_00(\Omega_1)*Y_00(\Omega_2)

!! We include a factor 4\pi, to take into account another Spherical Harmonic
!! for each component.

!! We divide by two because of the combinatorials of the electrons 

density2=density2*4.0d0*acos(-1.0d0)/2

  return
End subroutine phi2_fedvr3d_coupled

!! subroutine plot_density2 calculates two-body radial density just to plot normalized to (number_e-1)x number_e/2

subroutine plot_density2(element,basis, weights,density2,density2_plot)
  implicit none
  ! This subroutine produces the two particle density distribution.
  ! rho(r_1,r_2)
  !!INPUT
  real*8, intent(in), allocatable :: weights(:,:)
  !! weights are the weights for the element in 1st argument and basis
  !! in the element for the 2nd argument
  integer, intent(in) :: element(:), basis(:) 
  !! element and basis of the argument
  complex(kind=k2), allocatable,intent(in) :: density2(:,:)
  !! density written in the FE-DVR

  !! OUTPUT
  real(kind=k1), allocatable,intent(out) :: density2_plot(:,:)
  !! value of the two-body density to plot normalized to (number_e-1)x number_e/2 to plot

  !! AUXILIAR VARIABLES
  integer :: nradial
  integer :: k, k2 !! indexes in the loops to calculate ionization
  real(kind=k1) :: aux
  !! auxiliar variable to calculate it.

  !! Initialize the output
  
  nradial=size(density2(:,1))

  allocate(density2_plot(1:nradial,1:nradial))
  density2_plot=zero

   Do k=1,size(density2(1,:))-1 !! it is not plotting the last point
      Do k2=1,size(density2(:,1))-1 !! it is not plotting the last point
         !! Plot rho2
         !! Take into account the weights on the points.
         If (basis(k).lt.basis(k+1).and.basis(k2).lt.basis(k2+1)) then 
            !! basis(k) and basis(k2) are not the last nodes of the element.
            density2_plot(k,k2)= real(density2(k,k2))/(weights(element(k),basis(k))*weights(element(k2),basis(k2)))
         Elseif (basis(k).gt.basis(k+1).and.basis(k2).lt.basis(k2+1)) then 
            !! basis(k) is the last node of the element, but not basis(k2).
            aux=((weights(element(k),basis(k))+weights(element(k+1),1))*weights(element(k2),basis(k2)))
            density2_plot(k,k2)= real(density2(k,k2))/(aux)
         Elseif (basis(k2).gt.basis(k2+1).and.basis(k).lt.basis(k+1)) then 
            !! basis(k2) is the last node of the element, but not basis(k).
            aux=((weights(element(k2),basis(k2))+weights(element(k2+1),1))*weights(element(k),basis(k)))
            density2_plot(k,k2)= real(density2(k,k2))/(aux)
         Elseif (basis(k2).gt.basis(k2+1).and.basis(k).gt.basis(k+1)) then 
            aux=((weights(element(k2),basis(k2))+weights(element(k2+1),1))*(weights(element(k),basis(k))+weights(element(k+1),1)))
            density2_plot(k,k2)= real(density2(k,k2))/(aux)
         Endif
      End Do
   End Do

  return
End subroutine plot_density2

!! The subroutine double_ionization calculates the probability of the electrons to be localized
!! further than a radius rout

subroutine double_ionization(rout,xglobal,density2,d_ionization)
Implicit none

!! INPUT
real(kind=k1), intent(in) :: rout
!! rout defines the outer region
real(kind=k1), intent(in), allocatable :: xglobal(:)
!! nodes of the radial FE-DVR
complex(kind=k1), allocatable,intent(in) :: density2(:,:)
!! two body density in the FE-DVR basis

!! OUTPUT
real(kind=k1), intent(out) :: d_ionization
!! number_electrons x (number_electrons-1)/2 - P(r1<r_out,r2<rout)

!! AUXILIAR VARIABLES
integer :: k, k2
!! variables in loops

!! initialize the output

d_ionization = zero

Do k=1, size(density2(:,1))
   Do k2=1, size(density2(:,1))
         If (xglobal(k).gt.rout.and.xglobal(k2).gt.rout) exit !! we only sum the contribution of both electrons in the inner region
         !! or the one in the inner region and another in the outer region.
      d_ionization=d_ionization+real(density2(k,k2)) !! yield in the inner region
   End Do
End Do


!!$Do k=1, size(density2(:,1))
!!$   If (xglobal(k).gt.rout) exit !! we only sum the contribution of the inner region
!!$   Do k2=1, size(density2(:,1))
!!$      If (xglobal(k2).gt.rout) exit !! we only sum the contribution of the inner region
!!$      d_ionization=d_ionization+real(density2(k,k2)) !! yield in the inner region
!!$   End Do
!!$End Do

d_ionization=dble(system%numelectron*(system%numelectron-1))/2.0d0-d_ionization !! yield in the outer region

return
End subroutine double_ionization

end module analysis_part


subroutine drive_analysis_part(time,xd,yd,zd,xdtt,ydtt,zdtt,icount)
 use analysis_part
 use operator_spatial_space
 use fedvr3d_basis_set 

 implicit none
 complex(kind=k2) :: onebody_ener,twobody_ener,total_energy
 integer :: icount
 real(kind=k1) :: time
 real*8 :: xd,yd,zd !! expectation value of the dipole
 real*8 :: xdtt,ydtt,zdtt !! expectation value of the dipole
 complex(kind=k2), allocatable :: phik(:,:) !! Orbitals in momentum space
 complex(kind=k2), allocatable :: phik_hamming(:,:) !! Orbitals in momentum space using Hamming Window
 real(kind=k1), allocatable :: fedvr_k_weight(:,:), fedvr_k_node(:,:) !! weights and nodes in the element in the first argument and for the point in the element, given by the second argument.
 real(kind=k1), allocatable :: fedvr_k_global(:) !! fedvr node in the momentum space.
 integer :: idicp,jdicp,iangle
 real(kind=k1) :: angle_polar !! Polar angle in the pxpz plane
!! Variables to define the 3D density
 integer, allocatable :: global_to_local_rho(:,:)
!! global_to_local_rho(#1,1)-> FE-DVR function of #1
!! global_to_local_rho(#1,2)-> angular momentum of #1
!! global_to_local_rho(#1,3)-> magnetic quantum number of #1 

 complex(kind=k2), allocatable :: rho_phi(:,:) !! monoparticular 3D density for each orbitals
!! Density of a single orbital.
 complex(kind=k2), allocatable :: rho_all(:) !! 3D density for all the electrons
!! Density the whole wavefunction.
 complex(kind=k2), allocatable :: rho_all_hamming(:) !! 3D density for all the electrons using the hamming window for the photoelectron spectrum
!! Density the whole wavefunction.
complex(kind=k2), allocatable :: diff(:,:), diff_hamming(:,:) !! differential probability
complex(kind=k2), allocatable :: pxpz_dens(:,:) !! momentum density on the pxpz plane
complex(kind=k2), allocatable :: rhor_phi(:,:) !! monoparticular 3D density for each orbitals
!! Density of a single orbital.
 complex(kind=k2), allocatable :: rhor_all(:) !! 3D density for all the electron
!! Radial density for ionization
 real(kind=k1), allocatable :: dens(:) !! Density
 real(kind=k1), allocatable :: dens_hamming(:) !! Density in momentum space using Hamming Window
!! 2R density 
  complex(kind=k2), allocatable:: density2r(:,:) !! Two-Body density in the FE-DVR
  real(kind=k1), allocatable:: density2r_plot(:,:) !! Two-Body density in the FE-DVR
  real(kind=k1) :: double_ion !! yield of the double ionization
  real(kind=k1) :: sum_twobody !! norm of the double ionization

!! Auxiliar variables
Character*24:: rhoall
Character*22:: pxpz_file
Character*35:: rhophi
Character*24:: rhorall
Character*35:: rhorphi
Character*22:: rhok_radial
Character*27:: rho_r1_r2
real(kind=k1) :: weight, node !! auxiliar weights and nodes
real(kind=k1) :: sum_3rout2, sum_rout, sum_rout2, sum !! auxiliar variables for the ionization probability

INTERFACE
  FUNCTION sign_m1(m) 
    IMPLICIT NONE
    Integer, INTENT(IN)  :: m
    REAL*8              :: sign_m1
  END FUNCTION sign_m1
END INTERFACE
 
!! Update the density matrices

  call rho1()
  call rho2()

!
! calculation of the energy 
!

  If (analysis%i_energy==1) then

     call update_ham_onebody()
     !! update one body hamiltonian

     call update_vv_twobody()
     !! update two body hamiltonian
     
     call onebody_expectation(ham_onebody,cden1,system%nptot,onebody_ener)
     !! expectation value of the one body Hamiltonian
     call twobody_expectation(tei_spatial,cden2,system%nptot,twobody_ener)
     !! expectation value of the two body Hamiltonian
     
     total_energy = onebody_ener + twobody_ener
     !! total energy
     
     write(file1,*) time, real(total_energy)
     print*,'time=', time, '    energy (a.u.)=', real(total_energy)
     
  End If

!
! calc. auto function
!
if(analysis%i_auto==0) then


endif 

!! dipole moment in imaginary time propagation

If (.not.laser%tdornot) then
         !!Calculation of the expectation value of the dipole 
      Call dipole_expectation(phi,fedvrx_global,lm_l,lm_m,cden1,xmat_3d,ymat_3d,zmat_3d,fedvr3d%nb_angle,xd,yd,zd) !! calculation of the dipole moment
      
      write(101,*) time, xd, yd, zd !! time, <mu_x>, <mu_y> and <mu_z> in atomic units
      
      !!Calculation of the expectation value of the dipole acceleration.
      
      Call dipole_acceleration_expectation(phi,system%zz,fedvrx_global,lm_l,lm_m,cden1,xmat_3d,ymat_3d,zmat_3d,fedvr3d%nb_angle,xdtt,ydtt,zdtt) !! calculation of the dipole moment
      
      write(108,*) time, xdtt, ydtt, zdtt !! time, <\partial_{tt}mu_x>, <\partial_{tt}mu_y> and <\partial_{tt}mu_z> in atomic units

   Else  

      write(100,*) time, laser%ex_t,laser%ey_t,laser%ez_t
!! write the electric field

      write(123,*) time, laser%ax_t,laser%ay_t,laser%az_t

      !!

   End If

!
! calc. hhg
!
if(analysis%i_hhg==1.and.laser%tdornot.and.icount.le.prop%icount_max_pulse) then

      !!Calculation of the expectation value of the dipole 
      Call dipole_expectation(phi,fedvrx_global,lm_l,lm_m,cden1,xmat_3d,ymat_3d,zmat_3d,fedvr3d%nb_angle,xd,yd,zd) !! calculation of the dipole moment
      
      write(101,*) time, xd, yd, zd !! time, <mu_x>, <mu_y> and <mu_z> in atomic units
      
      !!Calculation of the expectation value of the dipole acceleration.
      
      Call dipole_acceleration_expectation(phi,system%zz,fedvrx_global,lm_l,lm_m,cden1,xmat_3d,ymat_3d,zmat_3d,fedvr3d%nb_angle,xdtt,ydtt,zdtt) !! calculation of the dipole moment
      
      write(108,*) time, xdtt, ydtt, zdtt !! time, <\partial_{tt}mu_x>, <\partial_{tt}mu_y> and <\partial_{tt}mu_z> in atomic units
   
   endif



!
! calc. twoe_density
!
if(analysis%i_twoe_density==1.and.icount.eq.prop%icount_max_pulse.and.laser%tdornot) then

   print*, 
   print*, 'P(r) at the end of the pulse'
   print*, 

   Call plot_phi_fedvr3d(phi,fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_w, cden1)

   print*, 
   print*, 'P(r_1,r_2) at the end of the pulse'
   print*, 'Using the coupled basis'
   print*, 

   call plot_phi2_fedvr3d_coupled(phi,fedvrx_global,which_element,which_basis, fedvr_w,cden2)

End if

!
! calc. i_ionization
!

If (analysis%i_ionization==1.and.laser%tdornot) then

!! prepare the file to store info

    Call radial_density(phi,fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_w, cden1,dens)  !! ionization in the radial coordinate

!! dens is the one body radial density in the coordinate.

    Call phi2_fedvr3d_coupled(phi,fedvrx_global,cden2,density2r) !! calculate the two-body density in FE-DVR basis

    Call double_ionization(fedvr3d%r_out,fedvrx_global,density2r,double_ion) !! yield of the double ionization
    Call double_ionization(fedvr3d%r_end,fedvrx_global,density2r,sum_twobody) !! norm of the double ionization

    sum=zero
!! set the variables to -1.0d0 if they have not been used.
    sum_rout2 = -1.d0 !! r_out/2 to infinity
    sum_rout = -1.d0 !! r_out to infinity
    sum_3rout2 = -1.d0 !! 3*r_out/2 to infinity

    Do idicp=1,size(dens)
!! electrons from 0 to r_out*3/2
       If (fedvrx_global(idicp).ge.fedvr3d%r_out*1.5d0.and.sum_3rout2.lt.zero) sum_3rout2=sum 
!! electrons from 0 to r_out
       If (fedvrx_global(idicp).ge.fedvr3d%r_out.and.sum_rout.lt.zero) sum_rout=sum
!! electrons from 0 to r_out/2
       If (fedvrx_global(idicp).ge.fedvr3d%r_out*.5d0.and.sum_rout2.lt.zero) sum_rout2=sum
       sum=sum+dens(idicp)
    End Do

    write(file3,*) time, sum, dble(system%numelectron)-sum_rout, dble(system%numelectron*(system%numelectron-1))/2.0d0-sum_twobody,double_ion

    deallocate(dens,density2r)

End if

!! Ionization after the pulse
!! Calculate the 3D one body density
print*, icount, prop%icount_max_pulse,prop%icount_max
print*, 'module_analysis.f90'

if(analysis%i_ionization==1.and.laser%tdornot.and.icount.ge.prop%icount_max_pulse) then

   call basis_set_rho3d(which_element,which_basis,fedvr3d%l_max,fedvr3d%m_max,global_to_local_rho) !! prepares the couple representation

   call rho3d_coupled(phi, cden1,rhor_phi, rhor_all)   
   !! calculates the density in 3D

!! Write the results in files

   If (allocated(rhor_phi)) then !! DEBUG 04.03.2015

      Do jdicp=1,size(rhor_phi(1,:))
         write(rhorphi,'(a15,i2.2,a6,i4.4,a4,i4.4)') 'rhor3d_orbital_',jdicp,'_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization

         Open(197+jdicp,file=rhorphi)
         write(197+jdicp,*) '#Time (a. u.) = ', time
         write(197+jdicp,*) 
         write(197+jdicp,*) '#Columns'
         write(197+jdicp,*) '#1: Radial node'
         write(197+jdicp,*) '#2: Weight'
         write(197+jdicp,*) '#3: Total angular momentum quantum number'
         write(197+jdicp,*) '#4: Magnetic quantum number in coupled basis'
         write(197+jdicp,*) '#5: Real part of the coefficient'
         write(197+jdicp,*) '#6: Imaginary part of the coefficient'
         write(197+jdicp,*) 
         write(197+jdicp,*) 

Do idicp=1,size(rhor_phi(:,1))
   weight=fedvr_w(which_element(global_to_local_rho(idicp,1)),which_basis(global_to_local_rho(idicp,1))) !! weight
   node=fedvrx_global(global_to_local_rho(idicp,1)) !! node
   
   write(197+jdicp,*) node, weight,global_to_local_rho(idicp,2),global_to_local_rho(idicp,3),real(rhor_phi(idicp,jdicp)), aimag(rhor_phi(idicp,jdicp))
End Do
close(197+jdicp)
End Do
      
      write(rhorall,'(a12,i4.4,a4,i4.4)') 'rhor3d_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization
      Open(198,file=rhorall)
      write(198,*) '#Time (a. u.) = ', time
write(198,*) 
write(198,*) '#Columns'
write(198,*) '#1: Radial node'
write(198,*) '#2: Weight'
write(198,*) '#3: Total angular momentum quantum number'
write(198,*) '#4: Magnetic quantum number in coupled basis'
write(198,*) '#5: Real part of the coefficient'
write(198,*) '#6: Imaginary part of the coefficient'
write(198,*) 
write(198,*) 
      
      Do idicp=1,size(rhor_phi(:,1))
         weight=fedvr_w(which_element(global_to_local_rho(idicp,1)),which_basis(global_to_local_rho(idicp,1))) !! weight
         node=fedvrx_global(global_to_local_rho(idicp,1)) !! node         
         write(198,*) node, weight,global_to_local_rho(idicp,2),global_to_local_rho(idicp,3),real(rhor_all(idicp)), aimag(rhor_all(idicp))
      End Do
      
      close(198)
      
   End If !! DEBUG 04.03.2015

   deallocate(global_to_local_rho,rhor_phi,rhor_all)

   !! radial distribution

   If (allocated(dens)) deallocate(dens) 

   Call radial_density(phi,fedvrx_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_w, cden1,dens)  !!  ionization in the radial coordinate
   
   !! Open the file for the radial part distribution
   write(rhok_radial,'(a10,i4.4,a4,i4.4)') 'rhor_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization
   Open(309,file=rhok_radial) 

   write(309,*) '# Time (atomic units):', time
   write(309,*) ''
   write(309,*) '# Column #1: Radial position, r, in atomic units'
   write(309,*) '# Column #2: Density'
   write(309,*) ''

   Do idicp=1,size(dens),1
      
      If (idicp.lt.size(dens)) then
         If (which_basis(idicp).lt.which_basis(idicp+1)) then !! basis(k) is not the last
            !! node of the element.
            !! adapt the density to plot   
            
            dens(idicp)=dens(idicp)/fedvr_w(which_element(idicp),which_basis(idicp))       
         else !! otherwise, for a bridge function (which_element(idicp)+1=which_element(idicp+1))
            dens(idicp)=dens(idicp)/(fedvr_w(which_element(idicp),which_basis(idicp))+fedvr_w(which_element(idicp+1),1))
            
         End If
         
         write(309,*) fedvrx_global(idicp), dens(idicp)
         
      End If
   End Do

deallocate(dens)

close(309)

   !! Open file for the two body density

   write(rho_r1_r2,'(a15,i4.4,a4,i4.4)') 'rho_r1_r2_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization
   Open(309,file=rho_r1_r2)

   write(309,*) '# Time (atomic units):', time
   write(309,*) ''
   write(309,*) '# Column #1: Radial position r1 in atomic units'
   write(309,*) '# Column #2: Radial position r2 in atomic units'
   write(309,*) '# Column #3: Two-Body Density'
   write(309,*) ''

   call phi2_fedvr3d_coupled(phi,fedvrx_global,cden2,density2r) !! calculation of the Two-Body density

   call plot_density2(which_element,which_basis, fedvr_w,density2r,density2r_plot) !! calculation of the Two-Body radial density

   !! Plot the Two-Body radial density

   Do idicp=1,size(density2r(:,1))-1
      Do jdicp=1,size(density2r(:,1))-1
         write(309,*) fedvrx_global(idicp),fedvrx_global(jdicp),density2r_plot(idicp,jdicp)
      End Do
      write(309,*)
   End Do
   
close(309)

endif  !! CHECK 03.06.2015

!
!calc. the momentum distribution
!

!! Calculation of the 3D and radial monoparticular density in momentum space.
!!

if(analysis%i_momentum==1.and.icount.ge.prop%icount_max_pulse.and.laser%tdornot) then

print*, 'module_analysis.f90, pre initial_fedvr3d_radial_momentum'

!!$write(*,*) '' 
!!$write(*,*) '#============================================================='
!!$write(*,*) '# CALCULATION OF THE MOMENTUM DISTRIBUTION'
!!$write(*,*) '#============================================================='
!!$write(*,*) '' 
!!$write(*,*) 'Calculation of the wavefunction in the momentum space' 
!!$write(*,*) '' 

!! Calculation of the FE-DVR + SPHERICAL HARMONICS

call initial_fedvr3d_radial_momentum(fedvr3d%k_max,fedvr3d%k_min,which_element,which_basis,fedvr_k_weight,fedvr_k_node, fedvr_k_global) 

print*, 'module_analysis.f90, pre phi_momentum'

!! Orbitals in momentum space

call phi_momentum(phi,fedvr3d%r_out,lm_l,lm_m,fedvrx_global, fedvr_w,which_element,which_basis,fedvr_k_global,fedvr_k_weight,phik) !! calculation of the reduced radial wavefunction in FE-DVR times the spherical harmonics

call phi_momentum_hamming(phi,fedvr3d%r_out,fedvr3d%r_end,20.0d0,lm_l,lm_m,fedvrx_global, fedvr_w,which_element,which_basis,fedvr_k_global,fedvr_k_weight,phik_hamming) !! calculation of the reduced radial wavefunction in FE-DVR times the spherical harmonics

!! Radial density in momentum space (the output is in the FE-DVR). 

print*, 'module_analysis.f90, pre deallocate(dens)'

If (allocated(dens)) deallocate(dens)

!! phik are the orbitals in the momentum space

print*, 'module_analysis.f90, pre radial_density'

Call radial_density(phik,fedvr_k_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_k_weight, cden1,dens)  !! radial density in the momentum space.

Call radial_density(phik_hamming,fedvr_k_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_k_weight, cden1,dens_hamming)  !! radial density in the momentum space using Hamming Window.

!! Open the file for the radial part distribution
write(rhok_radial,'(a10,i4.4,a4,i4.4)') 'rhok_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization
Open(309,file=rhok_radial) 

write(309,*) '# Time (atomic units):', time
write(309,*) ''
write(309,*) '# Column #1: Radial momentum, k_r, in atomic units'
write(309,*) '# Column #2: Density'
write(309,*) '# Column #3: Density using a Hamming Window'
   write(309,*) ''

Do idicp=1,size(dens),1

   If (idicp.lt.size(dens)) then
      If (which_basis(idicp).lt.which_basis(idicp+1)) then !! basis(k) is not the last
         !! node of the element.
         !! adapt the density to plot   
         
         dens(idicp)=dens(idicp)/fedvr_k_weight(which_element(idicp),which_basis(idicp))
         dens_hamming(idicp)=dens_hamming(idicp)/fedvr_k_weight(which_element(idicp),which_basis(idicp))
      else !! otherwise, for a bridge function (which_element(idicp)+1=which_element(idicp+1))
         dens(idicp)=dens(idicp)/(fedvr_k_weight(which_element(idicp),which_basis(idicp))+fedvr_k_weight(which_element(idicp+1),1))
         dens_hamming(idicp)=dens_hamming(idicp)/(fedvr_k_weight(which_element(idicp),which_basis(idicp))+fedvr_k_weight(which_element(idicp+1),1))
      End If
      
      write(309,*) fedvr_k_global(idicp), dens(idicp),dens_hamming(idicp)
   End If
End Do

deallocate(dens,dens_hamming)

close(309)

   !! Open file for the two body density

   write(rho_r1_r2,'(a15,i4.4,a4,i4.4)') 'rho_k1_k2_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization
   Open(309,file=rho_r1_r2)

   write(309,*) '# Time (atomic units):', time
   write(309,*) ''
   write(309,*) '# Column #1: Radial momentum k1 in atomic units'
   write(309,*) '# Column #2: Radial momentum k2 in atomic units'
   write(309,*) '# Column #3: Two-Body Density'
   write(309,*) ''

   if (allocated(density2r)) deallocate(density2r) !! deallocate the variable to store density in FE-DVR.

   call phi2_fedvr3d_coupled(phik,fedvr_k_global,cden2,density2r) !! calculation of the Two-Body density

   if (allocated(density2r_plot)) deallocate(density2r_plot) !! deallocate the variable to store density to plot.

   call plot_density2(which_element,which_basis, fedvr_k_weight,density2r,density2r_plot) !! calculation of the Two-Body radial density

   !! Plot the Two-Body radial density in momentum space

   Do idicp=1,size(density2r(:,1))-1
      Do jdicp=1,size(density2r(:,1))-1
         write(309,*) fedvr_k_global(idicp),fedvr_k_global(jdicp),density2r_plot(idicp,jdicp)
      End Do
      write(309,*)
   End Do
   
close(309)



!! Prepare the two radial density in momentum space

!!$call plot_phik_fedvr3d(phik,fedvr_k_global,fedvr3d%nb_angle,which_element,which_basis,lm_l,lm_m, fedvr_k_weight,cden1) !! This plots the wavefunction and calculates the corresponding

!! Calculation of the spatial basis set for the density in 3D.

call basis_set_rho3d(which_element,which_basis,fedvr3d%l_max,fedvr3d%m_max,global_to_local_rho)

!! Calculation of the 3D monoparticular density.
!! rho_phi is the density of each orbital.
!! rho_all is the total density.

!!$   print*, 
!!$   print*, 'Using the coupled basis'
!!$   print*, 

call rho3d_coupled(phik, cden1,rho_phi, rho_all)   !! without Hamming Window
call rho3d_coupled(phik_hamming, cden1,rho_phi, rho_all_hamming)   !! without Hamming Window

!!$   print*, 
!!$   print*, 'Not using the coupled basis'
!!$   print*, 
!!$
!!$
!!$call rho3d(phik,fedvr_k_global,global_to_local, global_to_local_rho,cden1,rho_phi, rho_all) !! CONSTRUCTION OF THE DENSITY IN 3D TOTAL WAVEFUNCTION AND EACH ORBITAL

!! Write the results in files
If (allocated(rho_phi)) then

Do jdicp=1,size(rho_phi(1,:))

   write(rhophi,'(a15,i2.2,a6,i4.4,a4,i4.4)') 'rhok3d_orbital_',jdicp,'_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization
   Open(197+jdicp,file=rhophi)
write(197+jdicp,*) '#Time (a. u.) = ', time
write(197+jdicp,*) 
write(197+jdicp,*) '#Columns'
write(197+jdicp,*) '#1: Radial node'
write(197+jdicp,*) '#2: Weight'
write(197+jdicp,*) '#3: Total angular momentum quantum number'
write(197+jdicp,*) '#4: Magnetic quantum number in coupled basis'
write(197+jdicp,*) '#5: Real part of the coefficient'
write(197+jdicp,*) '#6: Imaginary part of the coefficient'
write(197+jdicp,*) 
write(197+jdicp,*) 

   Do idicp=1,size(rho_phi(:,1))  !! Run in the orbitals

      weight=fedvr_k_weight(which_element(global_to_local_rho(idicp,1)),which_basis(global_to_local_rho(idicp,1))) !! weight
      node=fedvr_k_global(global_to_local_rho(idicp,1)) !! node

      write(197+jdicp,*) node, weight,global_to_local_rho(idicp,2),global_to_local_rho(idicp,3),real(rho_phi(idicp,jdicp)), aimag(rho_phi(idicp,jdicp))

   End Do
   close(197+jdicp)
End Do

Endif

If (allocated(rho_all)) then

write(rhoall,'(a12,i4.4,a4,i4.4)') 'rhok3d_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization

Open(198,file=rhoall)

write(198,*) '#Time (a. u.) = ', time
write(198,*) 
write(198,*) '#Columns'
write(198,*) '#1: Radial node'
write(198,*) '#2: Weight'
write(198,*) '#3: Total angular momentum quantum number'
write(198,*) '#4: Magnetic quantum number in coupled basis'
write(198,*) '#5: Real part of the coefficient'
write(198,*) '#6: Imaginary part of the coefficient'
write(198,*) '#7: Real part of the coefficient using the Hamming Window'
write(198,*) '#8: Imaginary part of the coefficient usig the Hamming Window'
write(198,*) 
write(198,*) 

   Do idicp=1,size(rho_all(:))
      If (which_element(global_to_local_rho(idicp,1)).eq.size(fedvr_k_weight(:,1))) cycle
      weight=fedvr_k_weight(which_element(global_to_local_rho(idicp,1)),which_basis(global_to_local_rho(idicp,1))) !! weight

      If (which_element(global_to_local_rho(idicp,1)).eq.which_element(global_to_local_rho(idicp,1))+1) then
         weight=fedvr_k_weight(which_element(global_to_local_rho(idicp,1)),which_basis(global_to_local_rho(idicp,1))) !! weight
      else
         weight=fedvr_k_weight(which_element(global_to_local_rho(idicp,1)),which_basis(global_to_local_rho(idicp,1)))+fedvr_k_weight(which_element(global_to_local_rho(idicp,1))+1,1) !! weight of the bridge function
      End If

      node=fedvr_k_global(global_to_local_rho(idicp,1)) !! node
      write(198,'(2E20.10,2I4.2,4E20.10)') node, weight,global_to_local_rho(idicp,2),global_to_local_rho(idicp,3),real(rho_all(idicp)), aimag(rho_all(idicp)),real(rho_all_hamming(idicp)), aimag(rho_all_hamming(idicp))
   End Do

close(198)

!! Now, we obtain the differential probability density for phi=0, and
!! theta= 0, pi/4, pi/2, 3pi/2, pi as a function of the radial momentum.

call rho3d_evaluation(rho_all,rho_all_hamming,global_to_local_rho,fedvr_k_global,diff,diff_hamming)

write(rhoall,'(a12,i4.4,a4,i4.4)') 'diff_k_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization

Open(198,file=rhoall)

write(198,*) '#Time (a. u.) = ', time
write(198,*) 
write(198,*) '#Columns'
write(198,*) '#1: Radial node of momentum'
write(198,*) '#2: Weight'
write(198,*) '#3: Theta=0'
write(198,*) '#4: Theta=pi/4'
write(198,*) '#5: Theta=pi/2'
write(198,*) '#6: Theta=3pi/4'
write(198,*) '#7: Theta=pi'
write(198,*) '#8: Theta=0 using Hamming Window'
write(198,*) '#9: Theta=pi/4 using Hamming Window'
write(198,*) '#10: Theta=pi/2 using Hamming Window'
write(198,*) '#11: Theta=3pi/4 using Hamming Window'
write(198,*) '#12: Theta=pi using Hamming Window'
write(198,*) 
write(198,*) 

   Do idicp=1,size(fedvr_k_global(:))-1
      
      If (which_element(idicp).eq.which_element(idicp+1)) then
         weight=fedvr_k_weight(which_element(idicp),which_basis(idicp)) !! weight
      else
         weight=fedvr_k_weight(which_element(idicp),which_basis(idicp))+fedvr_k_weight(which_element(idicp+1),1) !! weight of the bridge function
      End If

      node=fedvr_k_global(idicp) !! node
      write(198,'(12E20.10)') node, weight,real(diff(global_to_local_rho(idicp,1),1)),real(diff(global_to_local_rho(idicp,1),2)),real(diff(global_to_local_rho(idicp,1),3)),real(diff(global_to_local_rho(idicp,1),4)),real(diff(global_to_local_rho(idicp,1),5)),real(diff_hamming(global_to_local_rho(idicp,1),1)),real(diff_hamming(global_to_local_rho(idicp,1),2)),real(diff_hamming(global_to_local_rho(idicp,1),3)),real(diff_hamming(global_to_local_rho(idicp,1),4)),real(diff_hamming(global_to_local_rho(idicp,1),5))
   End Do
 
   close(198)

!! NOW WE PLOT THE DENSITY IN THE Px AND Pz PLANE ASSUMING THAT Py = 0
!! we set 

call pxpz_density(rho_all_hamming,global_to_local_rho,fedvr_k_global,100,pxpz_dens)

!! storing the results

write(pxpz_file,'(a10,i4.4,a4,i4.4)') 'pxpz_time_',icount*prop%icount_ionization/prop%icount_max,'_of_',prop%icount_ionization

Open(198,file=pxpz_file)

write(198,*) '#Time (a. u.) = ', time
write(198,*) 
write(198,*) '#Columns'
write(198,*) '#1: Radial node of momentum'
write(198,*) '#2: Weight'
write(198,*) '#3: Px'
write(198,*) '#4: Pz'
write(198,*) '#5: Density'
write(198,*) 
write(198,*) 

write(*,*) size(pxpz_dens(1,:)),size(pxpz_dens(:,1))

Do idicp=1,size(fedvr_k_global(:))-1
   
   If (which_element(idicp).eq.which_element(idicp+1)) then
      weight=fedvr_k_weight(which_element(idicp),which_basis(idicp)) !! weight
   else
      weight=fedvr_k_weight(which_element(idicp),which_basis(idicp))+fedvr_k_weight(which_element(idicp+1),1) !! weight of the bridge function
   End If
   
   node=fedvr_k_global(idicp) !! node
   Do iangle = 1, size(pxpz_dens(1,:))
      angle_polar = dble((iangle-1))/dble(size(pxpz_dens(1,:)))*acos(-1.0d0)*2.0d0
      write(198,'(5E20.10)') node, weight, node*sin(angle_polar),node*cos(angle_polar),real(pxpz_dens(idicp,iangle))
   End Do
      write(198,'(5E20.10)') node, weight, 0.0d0,node,real(pxpz_dens(idicp,1)) !! This one out of the loop is used to file momentum pxpz plane.
   write(198,*) ''

End Do
close(198)

End if

!! Deallocate auxiliar variables used

deallocate(fedvr_k_weight,fedvr_k_node,fedvr_k_global,phik,phik_hamming,global_to_local_rho,rho_phi,rho_all,pxpz_dens)


!!$write(*,*) '' 
!!$write(*,*) '#============================================================='
!!$write(*,*) '# DONE CALCULATION OF THE MOMENTUM DISTRIBUTION'
!!$write(*,*) '#============================================================='
!!$write(*,*) '' 

write(*,*) 'Ionization spectrum done!!'
print*, 'Line 1418 in module_analysis.f90'

endif

!
!calc. the atas
!
if(analysis%i_atas==0.and.laser%tdornot) then

endif

!! print the laser envelope

 return
end subroutine drive_analysis_part

!! The function besselj calculates the value of the spherical bessel function of integer index 'nu' at the real value 'x'

Function besselj(nu, x)
implicit none
!!INPUT
integer, intent(in) :: nu !! index of the spherical Bessel function.
real*8, intent(in) :: x !! argument of the spherical Bessel function.
!!OUTPUT
real*8 :: besselj !! spherical Bessel function.
!!Auxiliar variables
integer :: ii, jj, kk 
real*8 :: jn, jn_1,jn_plus_1
real*8 :: pi

pi=acos(-1.0d0)

besselj=0.0d0

jn_1=sqrt(2.0d0/(pi*x))*sin(x) !! J_1/2(x)

If (nu.eq.0) then
   besselj=jn_1*sqrt(pi/(2.0d0*x))
   return
end If

jn=sqrt((2.0d0*x)/pi)*(sin(x)/x**2.0d0-cos(x)/x)!! J_3/2(x)

If (nu.eq.1) then
   besselj=jn*sqrt(pi/(2.0d0*x))
   return
End If

Do ii=3,2*nu-1,2 !! Run to obtain the appropiate index of the Bessel function of the first kind, J_{nu+1/2}(x)
   jn_plus_1=2.0d0/x*(dble(ii)/2.0d0)*jn-jn_1
   jn_1=jn !!update for the next loop
   jn=jn_plus_1 !!update for the next loop
End Do

besselj=jn_plus_1*sqrt(pi/(2.0d0*x))

End Function besselj






