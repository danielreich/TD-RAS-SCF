
!
! same with haru
!
!==================================================================
subroutine orb_eqs_p_zet
!==================================================================
  use global
  implicit none
  call orb_eqs_p_zet2(system%NP1,system%NP2,system%NP0,system%NP0+system%NP1)
  return
end subroutine orb_eqs_p_zet
!==================================================================
subroutine orb_eqs_p_zet2(NNP1,NNP2,NNA1,NNA2)
!==================================================================
  use global
  use auxiliary
  use operator_spatial_space

  implicit none
  integer(4) :: NNP1,NNP2,NNP12,NNA1,NNA2
  integer(4) :: iijj,kkll,i,j,ii,jjp,k,l
  integer(4) :: info,ipiv2(NNP1*NNP2)
  complex(kind=k2) :: cs,ce2(NNP1*NNP2),ca2(NNP1*NNP2,NNP1*NNP2)
  complex(kind=k2) :: ca2_reg(NNP1*NNP2,NNP1*NNP2)
  complex(kind=k2) :: here_temp,here_temp1,here_temp2
  NNP12=NNP1*NNP2
!============================================================================
!  calculation of ce1 and ce1
!============================================================================
     iijj=0
     do j=1,NNP2
     do i=1,NNP1
     iijj=iijj+1
     kkll=0
     do k=1,NNP2
     do l=1,NNP1
     kkll=kkll+1
        ca2(iijj,kkll)=czet1(iijj,kkll)
     enddo
     enddo
     enddo
     enddo
!==============================================================================
!==============================================================================
!---- calculataion of ce1 -------------------
     iijj=0
     do j=1,NNP2
     do i=1,NNP1
     iijj=iijj+1
        ce2(iijj)=cv_zet2(iijj)
     enddo
     enddo
     call reg_mat(ca2,ca2_reg,NNP12,eps_Q)
     if (NNP2.gt.0) call zgesv(NNP12,1,ca2_reg,NNP12,ipiv2,ce2,NNP12,info)
!----------------
     iijj=0
     do j=1,NNP2
     jjp=j+NNA2
     do i=1,NNP1
     ii=i+NNA1
        iijj=iijj+1
        cs=-ce2(iijj)
        ch_dummy(jjp,ii)=cs
!        CE(ii,jjp)=dconjg(cs)

        here_temp = ham_onebody(ia( max(ii,jjp)    ) +min(ii,jjp) )

        if(jjp >= ii) then
          here_temp1 = here_temp          
          here_temp2 = dconjg(here_temp)
        else
          here_temp1 = dconjg(here_temp) 
          here_temp2 = here_temp
        endif        
 
!        CE(ii,jjp)=dconjg(cs-CH(jjp,ii))+CH(ii,jjp) ! <-- Be careful.
         ch_dummy(ii,jjp)=dconjg(cs-here_temp1)+ here_temp2 !! ch_dummy is h_i^j-i eta_i^j


      enddo
     enddo
  return
end subroutine orb_eqs_p_zet2
