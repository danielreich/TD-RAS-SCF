!==================================================================
subroutine zet_alg2
!==================================================================

  use global
  use auxiliary
  implicit none

  integer :: iab_ij,iab,jab,iz,izz,ijiab_ij,iel,jel,ipm1,ipm2,ijkl
  integer :: k,l,m,n,kk,ll,mm,nn,iii,jjj,kkk,lll,mmm,nnn,ijklmn
  integer :: kel,mel,klmn,ibox(system%numelectron),NNA1,NNA2

  NNA1=system%np0
  NNA2=system%np0+system%np1

  ijkl=0
  do iab_ij=1,N_STRING_AB_Hij
     iii=I_STRING_AB_Hij(iab_ij,system%numelectron+1)
     jjj=I_STRING_AB_Hij(iab_ij,system%numelectron+2)
     iz =I_STRING_AB_Hij(iab_ij,system%numelectron+3)
     iab=int(abs(iz))
     ipm1=iab/iz
     do kel=1,system%numelectron
        k=I_STRING_AB_Hij(iab_ij,kel)
        kk=(k-1)/system%nptot
        kkk=k-kk*system%nptot
     do l=1,system%nptot*2
        ll=(l-1)/system%nptot
        lll=l-ll*system%nptot
     if ((NNA1.lt.lll).and.(lll.le.NNA2)&
   &.and.(NNA2.lt.kkk).and.(kk.eq.ll)) then

        do iel=1,system%numelectron
        jel=I_STRING_AB_Hij(iab_ij,iel)
        if (jel.eq.k) jel=l
        ibox(iel)=jel
        enddo

        call asc_order(ibox,ipm2,jab)

        if (jab.ne.0) then
!=================================================================
           ijkl=ijkl+1        ! counter increment
           if (L_ZMP1.lt.ijkl) stop '** ERROR IN zet_alg2 L_ZMP1 **'
           I_ZMP1(ijkl,1)=iii
           I_ZMP1(ijkl,2)=jjj
           I_ZMP1(ijkl,3)=kkk
           I_ZMP1(ijkl,4)=lll
           I_ZMP1(ijkl,5)=iab
           I_ZMP1(ijkl,6)=ipm1*ipm2*jab
!=================================================================
        endif

     endif
     enddo
     enddo
  enddo
  N_ZMP1=ijkl

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

  ijklmn=0  ! counter reset
  do iab_ij=1,N_STRING_AB_Hij
     iii=I_STRING_AB_Hij(iab_ij,system%numelectron+1)
     jjj=I_STRING_AB_Hij(iab_ij,system%numelectron+2)
     iz =I_STRING_AB_Hij(iab_ij,system%numelectron+3)
     iab=int(abs(iz))
     ipm1=iab/iz
     do mel=2,system%numelectron
     do kel=1,mel-1
        k  =I_STRING_AB_Hij(iab_ij,kel)
        kk =(k-1)/system%nptot
        kkk=k-kk*system%nptot
        m  =I_STRING_AB_Hij(iab_ij,mel)
        mm =(m-1)/system%nptot
        mmm=m-mm*system%nptot
     do n=2,system%nptot*2  ! m and n should be the same spin state
     do l=1,n-1       ! k and l should be the same spin state
        nn =(n-1)/system%nptot
        nnn=n-nn*system%nptot
        ll =(l-1)/system%nptot
        lll=l-ll*system%nptot
     if ((kk+mm).eq.(ll+nn)) then

!---- two-electron replacement (l,n) <-- (k,m) ----
        do iel=1,system%numelectron
        jel=I_STRING_AB_Hij(iab_ij,iel)
        if (jel.eq.k) jel=l
        if (jel.eq.m) jel=n
        ibox(iel)=jel
        enddo

        call asc_order(ibox,ipm2,jab)

        if (jab.ne.0) then
!=================================================================
           ijklmn=ijklmn+1        ! counter increment
           if (L_ZMP2.lt.ijklmn) stop '** ERROR IN zet_alg2 L_ZMP2 **'
           I_ZMP2(ijklmn,1)=iii
           I_ZMP2(ijklmn,2)=jjj
           I_ZMP2(ijklmn,3)=k
           I_ZMP2(ijklmn,4)=l
           I_ZMP2(ijklmn,5)=m
           I_ZMP2(ijklmn,6)=n
           I_ZMP2(ijklmn,7)=iab
           I_ZMP2(ijklmn,8)=ipm1*ipm2*jab
!         write(*,'(50g6.3)') iii,jjj,";",kkk,mmm,nnn,lll,"<---"
!=================================================================
        endif

     endif
     enddo
     enddo
     enddo
     enddo
  enddo
  N_ZMP2=ijklmn
  return
end subroutine zet_alg2
