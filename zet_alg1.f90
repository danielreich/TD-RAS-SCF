!=================================================================
subroutine zet_alg1
!=================================================================

  use global
  use auxiliary
  implicit none
  
  integer :: jab,jab_h,judge,judge1,judge2,jab_ij,ipm,iab_ij,iz,iab
  integer :: iel,jel,i,ii,iii,j,jj,jjj,ie,ibox(system%numelectron),NNA1,NNA2

  NNA1=system%np0
  NNA2=system%np0+system%np1

  jab_h=0
  do jab=1,n_string_ab
     judge=0
     do jel=1,system%numelectron
        j=i_string_ab(jab,jel)
        jj=(j-1)/system%nptot
        jjj=j-jj*system%nptot
        if (jjj.gt.NNA2) judge=judge+1
     enddo

     if (judge.eq.system%io_sdtq) then
        jab_h=jab_h+1
        do ii=1,system%numelectron
           i_string_ab_h(jab_h,ii)=i_string_ab(jab,ii)
        enddo
        i_string_ab_h(jab_h,system%numelectron+1)=jab
        do ii=1,system%nptot*2
           i_string_ab2_H(jab_h,ii)=i_string_ab2(jab,ii)
        enddo
        i_string_ab2_h(jab_h,system%nptot*2+1)=jab
     endif
  enddo
  n_string_ab_h=jab_h



  jab_ij=0
  do jab_h=1,n_string_ab_H
     do ie=1,system%numelectron
     i  =i_string_ab_H(jab_h,ie)
     ii =(i-1)/system%nptot
     iii=i-ii*system%nptot
     if ((ii.eq.0).and.(NNA1.lt.iii).and.(iii.le.NNA2)) then !Only spin UP state
     do jjj=NNA2+1,system%nptot

!--- there should be iii but not jjj ---
     judge1=0
     judge2=0
     do iel=1,system%numelectron/2
        jel=i_string_ab_H(jab_h,iel)
        if (jjj.ne.jel) judge2=judge2+1
     enddo
     if (judge2.eq.system%numelectron/2) then

!---- one-electron replacement (jjj <-- iii) ----
        do iel=1,system%numelectron
        jel=i_string_ab_H(jab_h,iel)
        if (jel.eq.iii) jel=jjj
        ibox(iel)=jel
        enddo

        call asc_order(ibox,ipm,jab)

        if (jab.eq.0) then
!=================================================================
           jab_ij=jab_ij+1        ! counter increment
           do iel=1,system%numelectron
              i_string_ab_Hij(jab_ij,iel)=ibox(iel)
           enddo
           i_string_ab_Hij(jab_ij,system%numelectron+1)=iii
           i_string_ab_Hij(jab_ij,system%numelectron+2)=jjj
           i_string_ab_Hij(jab_ij,system%numelectron+3)=ipm*i_string_ab_H(jab_h,system%numelectron+1)
!=================================================================
        endif

     endif
     enddo
     endif
     enddo
  enddo
  n_string_ab_Hij=jab_ij

!----------------------------------------------------------------------------

  
  return
end subroutine zet_alg1












