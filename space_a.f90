!----------------------------------------------------------------------------------------
! solve the equation, reference from haru,lars, equation(12)
!
!
!            ______                                         ________
!   .        \         i           i         j              \         i k        j l 
! ci*C_i  =   \      (h  -  ci* eta ) <F_I |E  | wf> + 0.5*  \       v    <F_I |E    | wf>
!             /        j           j         i               /        j l        i k  
!            /__ij__                                        /__ijkl__
!           
!----------------------------------------------------------------------------------------
    subroutine acoff_eqs()
     use global
     use auxiliary
     use density
     use wfunction
     use operator_spatial_space
     use solver
     implicit none

     integer :: iab,jab,ipq,ipqrs,ip,iq,ir,is,iz,iip,iiq,iir,iis,iipr
     complex(kind=k2) :: cs,cpm,cvv,cenergy
     integer :: sign_here
     integer :: idicp,jdicp



     do iab=1,n_string_ab
      cs=zzero
      do ipq=1,num_amp1(iab)    ! wenliang find the matrix element in right hand, the details can be found in the paper
        ip =i_amp1(iab,ipq,1)  
        iq =i_amp1(iab,ipq,2)      ! iq <- ip
        jab =i_amp1(iab,ipq,3)
        sign_here = i_amp1(iab,ipq,4)

        cpm=dcmplx(dble(sign_here),0.d0)
        cs=cs+ch_dummy(ip,iq)*cpm*acoff(jab)
     enddo


      do ipqrs=1,num_amp2(iab)
        ip   =i_amp2(iab,ipqrs,1)
        iq   =i_amp2(iab,ipqrs,2)
        ir   =i_amp2(iab,ipqrs,3)
        is   =i_amp2(iab,ipqrs,4)

        jab  =i_amp2(iab,ipqrs,5)
        iip  =i_amp2(iab,ipqrs,6)
        iiq  =i_amp2(iab,ipqrs,7)
        iir  =i_amp2(iab,ipqrs,8)
        iis  =i_amp2(iab,ipqrs,9)
        sign_here = i_amp2(iab,ipqrs,10)
        
        cpm  =dcmplx(dble(sign_here),0.d0)
        cvv=zzero
        if ((iip.eq.iiq).and.(iir.eq.iis)) then
           cvv=tei_spatial(ip,ir,is,iq)
           if (iip.eq.iis) cvv=cvv-tei_spatial(ip,ir,iq,is)
        endif

        cs=cs+cvv*cpm*acoff(jab)
     enddo

    ka(iab)=-ci*cs
    
    enddo


return
end subroutine acoff_eqs
