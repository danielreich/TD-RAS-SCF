
!
! module density
! 
  module density
      use global
      use auxiliary
      use wfunction
      use operator_spatial_space
      implicit none
      complex(kind=k2),allocatable,save :: cden1(:,:),inverse_cden1(:,:),cden2(:,:,:,:)
      complex(kind=k2),allocatable,save :: cden3(:,:,:,:),inverse_cden1_fc(:,:)
     
 

  
      contains
!---------------------------------------------------------------------------------------------
!  FI is slater dertermin    
!                                        ____          
!                                        \
! wavefunction is expand : wavefunction = \     C_I * FI = C_1 * F1 + C_2 * F2+ .........
!                                         /
!                                        /___
!                                         I=1
!
!    j                          j
!rho1       =  <wavefunction | E    | wavefunciton> 
!    i                          i 
!            ______   
!            \       
!             \                      j 
!           = /    conjg(C_I)*  <FI|E  |wavefunction>           
!            /                       i 
!            ------- 
!    j
! rho    = cden1(i,j) = 
!    i
!----------------------------------------------------------------------------------------------

  subroutine rho1()
   implicit none
   integer(4) :: idicp,jdicp
   integer(4) :: I_capital,ij,i,j,which,sign_here
   
! initialzaiton of cden1    
  do idicp =1,system%nptot
    do jdicp =1,system%nptot
     cden1 =dcmplx(0.0d0,0.0d0)
    enddo
 enddo 
! accumlate the matrix element  
       do I_capital = 1,n_string_ab
       do ij=1,num_amp1(I_capital)
 
           i = i_amp1(I_capital,ij,1)
           j = i_amp1(I_capital,ij,2)  
           which = i_amp1(I_capital,ij,3)
           sign_here = i_amp1(I_capital,ij,4)

           cden1(i,j) = cden1(i,j) + dcmplx(dble(sign_here),0.0d0)*dconjg(acoff(I_capital))*acoff(which) 

          enddo

       enddo
   return
  end subroutine rho1



  subroutine inverse_rho1
   implicit none
  integer :: idicp,jdicp,kdicp
  integer :: ipiv2(system%nptot)
  integer :: info
  complex(8),dimension(system%nptot,system%nptot) :: den_temp,den_reg
  
  complex(8),dimension(system%nptot,system%nptot) :: den_temp1

!
! regu. the density matrix [1]
!


    do idicp =1, system%nptot
      do jdicp =1, system%nptot
       inverse_cden1(idicp,jdicp) = zzero
       den_temp(idicp,jdicp) = cden1(idicp,jdicp)
   
      enddo
      inverse_cden1(idicp,idicp) = zone
    enddo


    call reg_mat(den_temp,den_reg,system%nptot,eps_q) 

!
! find the inverse of the reg_matrix
!    
  
    call zgesv(system%nptot,system%nptot,den_reg,system%nptot,ipiv2,inverse_cden1,system%nptot,info)
     if(info.ne.0)  then
       write(*,*) '____________________________________________________________'
       write(*,*) 'fail to find the inverse of reg_density_matrix'
       write(*,*) '------------------------------------------------------------'
       stop 
    endif


  return 
  end subroutine inverse_rho1


!
! hf fcore 
!
  subroutine inverse_rho1_fc
   implicit none
  integer :: idicp,jdicp,kdicp
  integer :: ipiv2(system%nptot_fc)
  integer :: info
  complex(8),dimension(system%nptot_fc,system%nptot_fc) :: den_temp,den_reg
  complex(8),dimension(system%nptot_fc,system%nptot_fc) :: den_temp1

!
! regu. the density matrix [1]
!
    do idicp =1, system%nptot_fc
      do jdicp =1, system%nptot_fc
       inverse_cden1(idicp,jdicp) = zzero
       den_temp(idicp,jdicp) = cden1(idicp+system%n_offset,jdicp+system%n_offset)

      enddo
      inverse_cden1_fc(idicp,idicp) = zone
    enddo


    call reg_mat(den_temp,den_reg,system%nptot_fc,eps_q)
!
! find the inverse of the reg_matrix
!    

    call zgesv(system%nptot_fc,system%nptot_fc,den_reg,system%nptot_fc,ipiv2,inverse_cden1_fc,system%nptot_fc,info)
     if(info.ne.0)  then
       write(*,*) '____________________________________________________________'
       write(*,*) 'fail to find the inverse of reg_density_matrix'
       write(*,*) '------------------------------------------------------------'
       stop
    endif

  return
  end subroutine inverse_rho1_fc







!-----------------------------------------------------------------------------------------
! from haru
!
!
!
!
!-----------------------------------------------------------------------------------------
 subroutine rho2()
  implicit none
  integer(4) :: iab,jab,ipqrs,ip,iq,ir,is,iz,iip,iiq,iir,iis
  complex(8) :: cs,cpm
  integer(4) :: mmm
 
 

!----two-electron replacement-----------------------------------------
  do iq=1,system%nptot
  do ip=1,system%nptot
  do ir=1,system%nptot
  do is=1,system%nptot
  cden2(is,ir,ip,iq)=dcmplx(0.D0,0.D0)
  enddo
  enddo
  enddo
  enddo

  do iab=1,n_string_ab
     do ipqrs=1,num_amp2(iab)
        ip   =i_amp2(iab,ipqrs,1)
        iq   =i_amp2(iab,ipqrs,2)
        ir   =i_amp2(iab,ipqrs,3)
        is   =i_amp2(iab,ipqrs,4)
        iz   =i_amp2(iab,ipqrs,5)
        iip  =i_amp2(iab,ipqrs,6)
        iiq  =i_amp2(iab,ipqrs,7)
        iir  =i_amp2(iab,ipqrs,8)
        iis  =i_amp2(iab,ipqrs,9)
        mmm = i_amp2(iab,ipqrs,10)
        jab  =int(abs(iz))
        cpm  =dcmplx(dble(mmm),0.d0)
        cs   =cpm*dconjg(acoff(iab))*acoff(jab)

        if (((iip.eq.0).and.(iir.eq.0).and.(iis.eq.0).and.(iiq.eq.0)).or.&
           &((iip.eq.1).and.(iir.eq.1).and.(iis.eq.1).and.(iiq.eq.1))) then

           cden2(ip,ir,is,iq)=cden2(ip,ir,is,iq)+cs
           cden2(ir,ip,is,iq)=cden2(ir,ip,is,iq)-cs
           cden2(ip,ir,iq,is)=cden2(ip,ir,iq,is)-cs
           cden2(ir,ip,iq,is)=cden2(ir,ip,iq,is)+cs
        else 
           cden2(ip,ir,is,iq)=cden2(ip,ir,is,iq)+cs
           cden2(ir,ip,iq,is)=cden2(ir,ip,iq,is)+cs
        endif

     enddo
  enddo
  
  return
 end subroutine rho2


 subroutine rho3() 
  implicit none
  integer(4) :: ij,il,ik,ii,ip
  complex(kind=k2) :: cs
  do ij=1,system%nptot
  do il=1,system%nptot
  do ik=1,system%nptot
  do ii=1,system%nptot
  cs=zzero
  do ip=1,system%nptot
     cs=cs+inverse_cden1(ii,ip)*cden2(ip,ik,il,ij)
  enddo
  cden3(ii,ik,il,ij)=cs
  enddo
  enddo
  enddo
  enddo
  return
 end subroutine rho3

!
!
!
!
! the same with haru
!
!==================================================================
subroutine zet_mat1
!==================================================================
  
  implicit none
  integer :: iab,jab,iz,ijkl,ij,kl,i,j,k,l
  integer :: NNP1,NNP2,NNA1,NNA2
  complex(kind=k2) :: cs,cpm
  NNP1 = system%np1
  NNP2 = system%np2
  NNA1 = system%np0
  NNA2 = system%np0 + system%np1

  do kl=1,N_ATEN
  do ij=1,N_ATEN
     CZET1(ij,kl)=dcmplx(0.d0,0.d0)
  enddo
  enddo


  do ijkl=1,N_ZMP1
     i  =I_ZMP1(ijkl,1)-NNA1
     j  =I_ZMP1(ijkl,2)-NNA2
     k  =I_ZMP1(ijkl,3)-NNA2
     l  =I_ZMP1(ijkl,4)-NNA1
     iab=I_ZMP1(ijkl,5)
     iz =I_ZMP1(ijkl,6)
     ij=(j-1)*NNP1+i
     kl=(k-1)*NNP1+l
     jab=int(abs(iz))
     cpm=dcmplx(dble(iz/jab),0.d0)
     cs=cpm*dconjg(acoff(iab))*acoff(jab)
     CZET1(ij,kl)=CZET1(ij,kl)+cs
  enddo

  return
end subroutine zet_mat1
!==================================================================
subroutine zet_mat2
!==================================================================
 
  implicit none
  integer :: iab,jab,ij,kl,io1,io2,loop
  integer :: i,j,k,l,m,n,kk,ll,mm,nn,iii,jjj,kkk,lll,mmm,nnn
  integer :: kkkk,llll,mmmm,nnnn,ijkl,ijkl2,iz,ijklmn
  integer :: NNP1,NNP2,NNA1,NNA2
  complex(kind=k2) :: cvv,cs0,cs,cpm
  NNP1= system%np1
  NNP2= system%np2
  NNA1= system%np0
  NNA2= system%np0 + system%np1

  do ij=1,N_ATEN
     CV_ZET2(ij)=cs
  enddo


  do ijklmn=1,N_ZMP2
     iii=I_ZMP2(ijklmn,1)
     jjj=I_ZMP2(ijklmn,2)
     k  =I_ZMP2(ijklmn,3)
     l  =I_ZMP2(ijklmn,4)
     m  =I_ZMP2(ijklmn,5)
     n  =I_ZMP2(ijklmn,6)
     iab=I_ZMP2(ijklmn,7)
     iz =I_ZMP2(ijklmn,8)
     j=jjj-NNA2
     i=iii-NNA1
     ij=(j-1)*NNP1+i

     kk=(k-1)/system%nptot

     kkk=k-kk*system%nptot
     kkkk=0
     if ((NNA1.lt.kkk).and.(kkk.le.NNA2)) kkkk=1
     if (NNA2.lt.kkk) kkkk=2
     mm=(m-1)/system%nptot
     mmm=m-mm*system%nptot
     mmmm=0
     if ((NNA1.lt.mmm).and.(mmm.le.NNA2)) mmmm=1
     if (NNA2.lt.mmm) mmmm=2
     nn=(n-1)/system%nptot
     nnn=n-nn*system%nptot
     nnnn=0
     if ((NNA1.lt.nnn).and.(nnn.le.NNA2)) nnnn=1
     if (NNA2.lt.nnn) nnnn=2
     ll=(l-1)/system%nptot
     lll=l-ll*system%nptot
     llll=0
     if ((NNA1.lt.lll).and.(lll.le.NNA2)) llll=1
     if (NNA2.lt.lll) llll=2
     if ((kk+mm).eq.(ll+nn)) then
     if ((k.lt.m).and.(l.lt.n)) then
     if ((llll+nnnn).le.(kkkk+mmmm)) then
     io1=0
     io2=0
     if ((kk.eq.ll).and.(mm.eq.nn)) then
        io1=1
        if (kk.eq.mm) io2=1
     endif
     jab=int(abs(iz))
     cpm=dcmplx(dble(iz/jab),0.d0)
     cs0=cpm*dconjg(acoff(iab))*acoff(jab)
     cvv=zzero
     if (io1.eq.1) then
        cvv=tei_spatial(kkk,mmm,nnn,lll)
        if (io2.eq.1) cvv=cvv-tei_spatial(kkk,mmm,lll,nnn)
     endif

     CV_ZET2(ij)=CV_ZET2(ij)+cvv*cs0
     endif
     endif
     endif
  enddo

  return
end subroutine zet_mat2






end module density


 subroutine configure_density
  use density
  implicit none
  
  allocate(cden1(system%nptot,system%nptot))
  allocate(cden2(system%nptot,system%nptot,system%nptot,system%nptot))

  if(system%io_fcore/=1) then
    allocate(inverse_cden1(system%nptot,system%nptot))  
    allocate(cden3(system%nptot,system%nptot,system%nptot,system%nptot))
  else
    allocate( inverse_cden1_fc(system%nptot_fc,system%nptot_fc) )
  endif
 



  return
 end subroutine configure_density




!
! auxiary subroutine
!
! 
  subroutine which(iarray,kk)
   use global
   use wfunction
   implicit none
   integer(4),intent(out) :: kk
   integer(4),dimension(2*system%nptot) :: iarray
   integer(4) :: isum,idicp,jdicp
!
! very bad way, should be improved
!
   
   kk =0
   do idicp =1,n_string_ab
     isum=0  
     do jdicp =1,2*system%nptot
	  
	  if(iarray(jdicp) == i_string_ab2(idicp,jdicp)) then
         isum = isum +1
	  endif

	 enddo
     
    
     if(isum==2*system%nptot) then
	   kk = idicp
	   exit
     endif
   enddo

   return
  end subroutine which




