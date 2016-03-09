
!
! the same with haru
!
!==================================================================
subroutine zet_mat1
!==================================================================
  use global
  use auxiliary
  use wfunction
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
  use global
  use auxiliary
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
     cvv=c0
     if (io1.eq.1) then
        cvv=cv(kkk,mmm,nnn,lll)
        if (io2.eq.1) cvv=cvv-tei_spatial(kkk,mmm,lll,nnn)
     endif
 
     CV_ZET2(ij)=CV_ZET2(ij)+cvv*cs0
     endif
     endif
     endif
  enddo
 
  return
end subroutine zet_mat2
