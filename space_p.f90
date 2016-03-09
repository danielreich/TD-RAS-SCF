
 subroutine orb_eqs_p
 use global
 implicit none
 integer :: io,io_ir_now 

! io =1
! io_ir_now = 1


! if(system%io_fcore==1) io =0
! if((system%io_fcore==2).and.(io_ir_now==2)) io = 0

!if( (system%np0/=0 ).and. (io ==1)) then

!  call orb_eqs_p2(system%np0,system%np1+system%np2,0,system%np0)
! endif


 if((system%io_sdtq==0) .and. (system%np2/=0)) then

  call orb_eqs_p2(system%np1,system%np2,system%np0,system%np0+system%np1)

 endif
end subroutine orb_eqs_p


!==================================================================
  subroutine orb_eqs_p2(NNP1,NNP2,NNA1,NNA2)
! subroutine Pspace_eqs(NNP1,NNP2,NNA1,NNA2)
!==================================================================
   use global
   use density
   use operator_spatial_space

  implicit none

  integer :: NNP1,NNP2,NNP12,NNA1,NNA2
  integer :: iijj,kkll,i,ii,j,jjp,k,kkp,l,ll,m
  integer :: info,ipiv2(NNP1*NNP2)
  complex(kind=k2) :: cs,ce2(NNP1*NNP2)
  complex(kind=k2) :: ca2(NNP1*NNP2,NNP1*NNP2),ca2_reg(NNP1*NNP2,NNP1*NNP2)
  integer :: idicp,jdicp
  complex(kind=k2) :: here_temp,here_temp1,here_temp2
  NNP12=NNP1*NNP2

!============================================================================
!  calculation of ce1 and ce1
!============================================================================  

  do iijj=1,NNP1*NNP2
  j=(iijj-1)/NNP1+1
  i=iijj-(j-1)*NNP1  
  jjp=j+NNA2   ! wenliang ��Ӧ�Ÿ����Ŀռ������   ��ͬ�������NP0--[NP1+NP2]   NP1--NP2 
  ii =i+NNA1   ! wenliang ��Ӧ�ŵ����Ŀռ������ [ NP0 ] --> [NP1 + NP2], or [ NP1 ] --> [ NP2 ] 
 
  do kkll=1,NNP1*NNP2
  k=(kkll-1)/NNP1+1
  l=kkll-(k-1)*NNP1
  kkp=k+NNA2
  ll =l+NNA1
     cs=zzero
        if (jjp.eq.kkp) then !wenliang 
            cs=cs+cden1(ii,ll)
        endif
        if (ii.eq.ll) then
           cs=cs-cden1(kkp,jjp)
        endif
     ca2(iijj,kkll)=cs  !wenliang ,find the matrix, the important one, 

     enddo ! in order to find enta[ij],solve the amplitude equations. 
  enddo
!==============================================================================
!==============================================================================
!---- calculataion of ce1 -------------------
!!$  !$omp do
  do iijj=1,NNP1*NNP2
  j=(iijj-1)/NNP1+1
  i=iijj-(j-1)*NNP1
  jjp=j+NNA2
  ii =i+NNA1
     cs=zzero
        do k=1,system%nptot
        do l=1,system%nptot
        do m=1,system%nptot
           cs=cs+tei_spatial(jjp,m,l,k)*cden2(ii,m,l,k)&
               &-tei_spatial(k,l,m,ii)*cden2(k,l,m,jjp)                    
        enddo
        enddo
        enddo
        ce2(iijj)=-cs
  enddo


!!$  !$omp end do
!!$
!!$  !$omp single
  call reg_mat(ca2,ca2_reg,NNP12,1.0d-10)



  if (NNP2.gt.0) call zgesv(NNP12,1,ca2_reg,NNP12,ipiv2,ce2,NNP12,info)
!!$  !$omp end single 
!!$
!!$!----------------
!!$  !$omp do



  do iijj=1,NNP1*NNP2
  j=(iijj-1)/NNP1+1
  i=iijj-(j-1)*NNP1

  jjp=j+NNA2
  ii =i+NNA1
        cs=ce2(iijj)
        ch_dummy(jjp,ii)=cs
!        ce(ii,jjp)=dconjg(cs) !wenliang ������Haru,��ʽ36
         
        here_temp = ham_onebody(ia( max(ii,jjp)    ) +min(ii,jjp) )

        if(jjp>=ii)  then
          here_temp1 = here_temp
          here_temp2 = dconjg(here_temp)
        else
          here_temp1 = dconjg(here_temp)
          here_temp2 = here_temp
        endif 

     

        ch_dummy(ii,jjp)=dconjg(cs-here_temp1)+here_temp2 ! <-- Be careful.
!       ch_dummy(ii,jjp)=here_temp2 - dconjg(here_temp1)  + dconjg(cs)
!                                   jjp,ii       ii,jjp 

!        write(*,*) here_temp2 ,dconjg(here_temp1)  , dconjg(cs) 


 enddo

 
  return
end subroutine orb_eqs_p2
