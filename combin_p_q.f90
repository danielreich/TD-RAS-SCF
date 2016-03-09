
!------------_____-------------------------------------------------------------
!            \              j 
! |dph_i> =   \   |ph_j> eta     + Q |d_phi>
!             /             i
!            /_____
!------------------------------------------------------------------------------
   subroutine combin_p_q()
    use global
    use operator_spatial_space
    use solver
    implicit none
    integer(4):: ip,iq,ix,idicp
    complex(8):: cs ,here_temp
    complex(8),dimension(system%nptot,system%nptot) :: CE
!-----------------------------------------------------------------------------
! after solve the p space, we get  h - i eta = ce
! we need :   i*(ce - h) = eta
!-----------------------------------------------------------------------------
! from haru
!---- redefine ce ----
     do iq=1,system%nptot
     do ip=1,system%nptot

        here_temp = ham_onebody(ia(max(ip,iq)) + min(ip,iq))
        if(ip<iq) here_temp  = dconjg(here_temp)

        cs=ci*(ch_dummy(ip,iq)-here_temp)
        CE(ip,iq)=cs
      
      enddo
     enddo
    
!---- contribution from the p-space ----
     do iq=1,system%nptot
        do ix=1,fedvr3d%nb_r*fedvr3d%nb_angle
           cs=zzero
           do ip=1,system%nptot
              cs=cs+phi(ix,ip)*CE(ip,iq)
           enddo
           kphi(ix,iq)=kphi(ix,iq)+cs !wenliang, 方程21 贡献来自P和Q空间
           !! kphi(ix,iq) comes from solving space Q
        enddo
     enddo
    return
  end subroutine combin_p_q
