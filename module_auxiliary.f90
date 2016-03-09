!===============================================================================
! module auxiliary 
!===============================================================================
 module auxiliary
  use global
  implicit none
!---------------------information for alpha string -----------------------------
   integer,save :: n_string_a_fci         
   integer,allocatable,save :: i_string_a(:,:)
   integer,allocatable,save :: i_string_a2(:,:)
!---------------------information for alpha string + beta string----------------   
   integer,save :: n_string_ab
   integer,allocatable,save :: i_string_ab(:,:)
   integer,allocatable,save :: i_string_ab2(:,:)
!---------------------space and spino ------------------------------------------
   integer,allocatable,save :: space(:)
   integer,allocatable,save :: spino(:)

!  
   integer,allocatable,save :: i_string_ab_h(:,:),i_string_ab2_h(:,:)
   integer :: n_zmp1,n_zmp2
   integer :: n_string_ab_h

!
   integer,allocatable,save :: num_amp1(:),i_amp1(:,:,:),num_amp2(:),i_amp2(:,:,:)
   integer ::n_aten,n_string_ab_hij,l_zmp1,l_zmp2
   integer,allocatable,save :: i_zmp1(:,:),i_zmp2(:,:)
   integer, allocatable,save :: i_string_ab_hij(:,:)
   complex(kind=k2),allocatable,save :: czet1(:,:),cv_zet2(:)

  
!
! set the alphastring and betestring parameter
!
  integer,allocatable,save :: w_vertex(:,:),n_string_ab_inv(:)





  contains
!
! calculate the i_string_a, i_string_a2
!     nume  = numelectron/2
      subroutine config1()
       implicit none
       integer :: ifault,idicp,jdicp,kdicp
       real(kind=k1) :: resultp

       
       call factorialcnm(system%nptot,system%numelectron/2,resultp)
       n_string_a_fci =int(resultp)


!
! allocate the address
!
      allocate(i_string_a(n_string_a_fci,system%numelectron/2))
      allocate(i_string_a2(n_string_a_fci,system%nptot)) 
      i_string_a = 0
      i_string_a2 = 0
!
! get jout            
!     
       call allnr(system%nptot,system%numelectron/2,n_string_a_fci,i_string_a,ifault)
 
      
       if(ifault/=0) then
          write(*,*) 'error happen in subroutine config1'
          write(*,*) '--maybe the input parameter is something wrong--' 
          stop 
       endif

       write(*,*) 'total-spatial orbitals',system%nptot
       write(*,*) 
       write(*,*) '-----------i_string_a----------------' 
       write(*,*)            
       do idicp =1,n_string_a_fci
         write(*,*) (i_string_a(idicp,jdicp),jdicp=1,system%numelectron/2)
       enddo
!
! second quantuzation 
! 4 spatial orbital, 4 electrons | f1:1 f2:2 f3:3
! jout  1,2  jout2  1, 1, 0, 0
!       1,3         1, 0, 1, 0
!       1,4         1, 0, 0, 1
!       2,3         0, 1, 1, 0 
!       2,4         0, 1, 0, 1
!       3,4         0, 0, 1, 1
!
       do idicp = 1,n_string_a_fci
        do jdicp = 1,system%numelectron/2
          kdicp = i_string_a(idicp,jdicp)

          if(kdicp>=1.and.kdicp<=system%nptot) then
             i_string_a2(idicp,kdicp) = 1
           else
              write(*,*) 'error happen in subroutine config1' 
           endif
       enddo
   enddo   

       write(*,*) '--------------i_string_a2---------------'
       write(*,*)
       do idicp=1,n_string_a_fci
         write(*,*) (i_string_a2(idicp,jdicp),jdicp=1,system%nptot)
       enddo
   
   return
  end subroutine config1
!
! definiation the spatial orbital and spin orbital
! 
  subroutine space_plus_spino
   implicit none
   integer :: idicp
!
! allocate the address
!
   allocate(space(system%nptot))
   allocate(spino(2*system%nptot))
   
   space = -1
   spino = 0
!
! initial space
!
  do idicp = 1, system%nptot

    if(idicp<=system%np0) then
      space(idicp) = 0
    elseif(idicp>=(system%np0+1) .and. idicp<=(system%np0+system%np1)) then
      space(idicp) = 1
    elseif(idicp>=(system%np1+1) .and. idicp<=(system%np1+system%np2)) then
      space(idicp) = 2
    else
      write(*,*) 'error happen in subroutine space_plus_spino'
    endif
  enddo
!
! initial spino(:) 
!
  do idicp =1, 2*system%nptot
  
   if(idicp<=system%nptot) then
     spino(idicp) =  0
   else
     spino(idicp) =  1
   endif  
  enddo
!
! ia
!
  do idicp =1,10000
    ia(idicp) = idicp*(idicp-1)/2
  enddo

  return
 end subroutine space_plus_spino

!
! alpha string + beta string
!
! for this method, the N_alpha == N_beta
!
     
  subroutine config2()
  implicit none  
  integer(4) :: idicp,jdicp,kdicp
  integer(4) :: np1p,jab,jab_old,np0p
  character(len=50) :: fname

  allocate(i_string_ab(n_string_a_fci**2,system%numelectron))
  allocate(i_string_ab2(n_string_a_fci**2,system%nptot*2))
!
! combin a + b, the same with Haru
!
       jab = 0  ! sum the combin a+b fullfill the condition !
       jab_old=0
       system%io_orb_eqs_p = 1
       loop1: do idicp =1,n_string_a_fci 
         loop2: do jdicp =1,n_string_a_fci
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
! np0 = 0, it means that, no orbital are frozen in the calcuations
! all the electrons are active
!
           if(system%np0==0) then
               if(system%np2==0) then ! np0= 0 and np2 = 0, it means that, only np1 space contains orbitals == '  MCTDHF '
                    fname='mctdhf'  ! all configurations are included in the calcuations  
                    if(2*system%np1==system%numelectron) fname='tdhf' 
                    jab = jab +1
                    system%io_orb_eqs_p = 0
               else                ! np0= 0, but np2 /=0, it means that, excitation is allowed, s/d/sd are working normally, the sdtq are not debug still.   
                    if(system%io_sdtq==0) then  ! io_sdtq ==0, it means tdrasscf-d method
                        fname ='td-ras-scf-d without core'
                        call howmany(space,i_string_a(idicp,:),i_string_a(jdicp,:),system%numelectron/2,system%nptot,np1p,1)
                   !
                   !if the total numelectron are not allowed occupied in np1, then, there are two situation are allowed, first:the np1 are full with electron,
                   ! second : there are two electron are excited to np2 space
                   !
                        if(2*system%np1<system%numelectron) then
                           if((np1p==2*system%np1).or.(np1p==(2*system%np1-2))) jab = jab +1
                         else
                           if((np1p==system%numelectron).or.(np1p==(system%numelectron-2))) &
                           jab = jab +1
                         endif 
                     else ! io_sdqt/=0, 1,2,3  
                        if(system%io_sdtq==1) fname ='td-ras-scf-s w/o core'
                        if(system%io_sdtq==2) fname ='td-ras-scf-sd w/o core'
                        if(system%io_sdtq==3) fname ='td-ras-scf-sdt w/o core' 
                        call howmany(space,i_string_a(idicp,:),i_string_a(jdicp,:),system%numelectron/2,system%nptot,np1p,1) 

                        if(2*system%np1<system%numelectron) then
                          if(np1p>=2*system%np1-system%io_sdtq) jab = jab +1
                        else
                          if(np1p>=(system%numelectron-system%io_sdtq)) &
                          jab = jab +1
                        endif         
                     endif        
                endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else
             ! we need to figure out how many electron allocate in this np0 space, it need to be alaways full of electrons
             call howmany(space,i_string_a(idicp,:),i_string_a(jdicp,:),system%numelectron/2,system%nptot,np0p,0)
             if (np0p==2*system%np0) then

                if(system%np2==0) then ! td-cas-scf

                   fname='td-cas-scf w/ core'
                   jab = jab+1

                else

                   if(system%io_sdtq==0) then
                       fname = 'td-cas-scf w/ core'
                       call howmany(space,i_string_a(idicp,:),i_string_a(jdicp,:),system%numelectron/2,system%nptot,np1p,1)     
                       if(2*(system%np0+system%np1)<system%numelectron) then
                          if(np1p==2*system%np1 .or. np1p==2*system%np1-2) jab = jab +1
                       else
                          if((2*system%np0+np1p)==system%numelectron .or. ((2*system%np0+np1p) == (system%numelectron-2))) &
                          jab = jab +1    
                       endif
                   else
                                     
                      if(system%io_sdtq==1) fname ='td-ras-scf-s w/ core'
                      if(system%io_sdtq==2) fname ='td-ras-scf-sd w/ core'
                      if(system%io_sdtq==3) fname ='td-ras-scf-sdt w/ core'
                      call howmany(space,i_string_a(idicp,:),i_string_a(jdicp,:),system%numelectron/2,system%nptot,np1p,1) 

                      if(2*(system%np0+system%np1)<system%numelectron) then
                           if(np1p>=(2*system%np1-system%io_sdtq)) jab = jab +1
                      else
                           if((np1p+2*system%np0)>=(system%numelectron-system%io_sdtq)) jab = jab +1
                      endif  
                   endif

                endif
             endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 endif

     if(jab_old < jab) then

          do kdicp=1,system%numelectron

             if(kdicp<=system%numelectron/2) then
             
                i_string_ab(jab,kdicp) = i_string_a(idicp,kdicp)
             else
               i_string_ab(jab,kdicp) = i_string_a(jdicp,kdicp-system%numelectron/2) + system%nptot 
            endif
         enddo

          do kdicp=1,system%nptot*2
            if(kdicp<=system%nptot) then
              i_string_ab2(jab,kdicp) = i_string_a2(idicp,kdicp) 
            else
              i_string_ab2(jab,kdicp) = i_string_a2(jdicp,kdicp-system%nptot) 
            endif
         enddo


      endif

      jab_old = jab


       enddo loop2
   enddo loop1

   n_string_ab = jab
    write(*,*) '--------there are total number of cofficent A---------'
    write(*,*)
    write(*,*) '----------n_string_ab--------- :', n_string_ab
    write(*,*) '----------i_string_ab-----------'


    
    do idicp=1,n_string_ab
         write(*,*) (i_string_ab(idicp,jdicp),jdicp=1,system%numelectron)
    enddo




   

    write(*,*) '----------i_string_ab2-----------'
     

    do idicp=1, n_string_ab
      write(*,*) (i_string_ab2(idicp,jdicp),jdicp=1,2*system%nptot)
      print*, idicp,n_string_ab
    enddo

   return
  end subroutine config2
!===========================================================================
! input sp, is spatial orbital, nume = numelectron/2
! figure out howmany electron occupied in the space
!===========================================================================
  subroutine howmany(space,sp1,sp2,nume,nspatial,npe,which)
  implicit none
  integer nume,nspatial,npe
  integer space(nspatial),sp1(nume),sp2(nume)
  integer idicp,jdicp,kdicp,which
!
! find how many electrons locate in np0,np1,np2 range
!   

  npe = 0
   
  do idicp=1,nume

   jdicp = sp1(idicp)

   if(space(jdicp) == which) then
    npe = npe + 1
   endif
    
    kdicp = sp2(idicp)
   if(space(kdicp) == which) then
    npe = npe + 1
   endif
  enddo

  return
  end subroutine howmany

!-----------------------------------------------------------------------------------------------
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
!
! 
!------------------------------------------------------------------------------------------------
!
!                                            j 
! subroutine amplitude1 calculate the  <FI |E  | wavefunction>  NB: j,i is spatial orbital
!                                            i  
!   
!     
  subroutine amplitude1
    implicit none
    integer,dimension(2*system%nptot) :: iarray,iatemp
    integer,dimension(system%numelectron) :: iarrayp
    integer :: I_capital,m_spin,n_spin,num,method
    integer :: sign_here,iresult,kk
    integer :: idicp,jdicp,kdicp,spin1,ispin1,ispin2
    integer :: ip,jp,ipt,jpt


    allocate(num_amp1(n_string_ab))
    allocate(i_amp1(n_string_ab,100000,10))  !need to revised here


! two method, one is from haru
!------------------------------------------------------------------------------------------------
! another method 
!
!        + 
! < FI| Cm Cn | wavefunction>  NB: m,n is spin orbital,  Cm spin == Cn spin
!                                               +
! for every FI, there are many combination <FI|Cm Cn|wavefunction>
! count howmany is not zero. and remember the sign and the corrsponding the right Slater dertermin
!
!------------------------------------------------------------------------------------------------
   method = 3
   if(method ==0) then

   do I_capital = 1, n_string_ab

    num = 0
!deal with the problem in second representation
           do m_spin = 1, 2*system%nptot   
           do n_spin =1,2*system%nptot

       ! such as |0,1,1,0,1,1,0>
           iarray = i_string_ab2(I_capital,:)
! spin of m is equal to spin of n            
           if(spino(m_spin) == spino(n_spin)) then
              sign_here =1
              if(iarray(m_spin)==1) then
                 sign_here = sign_here*(-1)**sum(iarray(1 : m_spin-1))
                 iarray(m_spin) = 0 
                 
                 if(iarray(n_spin)==0) then
                 sign_here = sign_here*(-1)**sum(iarray(1 : n_spin-1))
                 iarray(n_spin) = 1
                 
                 call which(iarray,kk)
 
                 if(kk/=0) then
                    num = num +1
                    i_amp1(I_capital,num,1) = m_spin - spino(m_spin)*system%nptot
                    i_amp1(I_capital,num,2) = n_spin - spino(n_spin)*system%nptot
                    i_amp1(I_capital,num,3) = kk
                    i_amp1(I_capital,num,4) = sign_here
                  endif         
                 endif
              endif
              endif
            enddo
         enddo 
       num_amp1(I_capital) = num
     enddo
!deal with the problem in normal way
   elseif(method==1) then

   do idicp =1,n_string_ab
        num = 0
        do jdicp =1,system%numelectron
          spin1 = i_string_ab(idicp,jdicp) 
          ispin1 = spino(spin1)
   
         do kdicp = 1,2*system%nptot
          ispin2 = spino(kdicp)   
              
           if(ispin1==ispin2) then

           iarrayp = i_string_ab(idicp,:)
           iarrayp(jdicp) = kdicp

               call  findindex(iarrayp,iresult,kk)

               if(kk/=0) then
                call  getsign(i_string_ab(iresult,:),iarrayp,system%numelectron,sign_here)
                 
                num = num + 1
                 i_amp1(idicp,num,1) = spin1 - spino(spin1)*system%nptot
                 i_amp1(idicp,num,2) = kdicp - spino(kdicp)*system%nptot
                 i_amp1(idicp,num,3) = iresult
                 i_amp1(idicp,num,4) = sign_here
              endif
             endif
           enddo
        enddo
     num_amp1(idicp) = num
 enddo   
   else  

!------------------------------------------------------------------------------
!         jp  
! <F_I | E  | wavefunction>  ! here, jp,ip is spatial orbitals
!         ip
!
!         +          + 
!= <F_I |C   C   +  C   C   | wavefunction>
!         ip+  jp+     ip-  jp-  
!------------------------------------------------------------------------------

   do I_capital = 1, n_string_ab
   
     num = 0

     do ip =1,system%nptot
        do jp =1,system%nptot
!
! spin + + 
!
         iatemp = i_string_ab2(I_capital,:)   
         ipt = ip
         jpt = jp
         sign_here = 1

         if(iatemp(ipt) ==1 ) then
           sign_here = sign_here*(-1)**sum(iatemp(1:ipt-1))
           iatemp(ipt) = 0

           if(iatemp(jpt) ==0) then

           sign_here = sign_here*(-1)**sum(iatemp(1:jpt-1))
              iatemp(jpt) = 1  
           
          call which(iatemp,kk)
 
         if(kk/=0) then
            num = num +1
            i_amp1(I_capital,num,1) = ip
            i_amp1(I_capital,num,2) = jp
            i_amp1(I_capital,num,3) = kk
            i_amp1(I_capital,num,4) = sign_here
         endif         
         endif
       endif
!
! spin - - 
!
       iatemp = i_string_ab2(I_capital,:)   
       ipt = ip + system%nptot
       jpt = jp + system%nptot
        sign_here = 1

      if(iatemp(ipt) ==1 ) then
        sign_here = sign_here*(-1)**sum(iatemp(1:ipt-1))
        iatemp(ipt) = 0

        if(iatemp(jpt) ==0) then

           sign_here = sign_here*(-1)**sum(iatemp(1:jpt-1))
           iatemp(jpt) = 1
  
           call which(iatemp,kk)
 
           if(kk/=0) then
               num = num +1
               i_amp1(I_capital,num,1) = ip 
               i_amp1(I_capital,num,2) = jp 
               i_amp1(I_capital,num,3) = kk
               i_amp1(I_capital,num,4) = sign_here
             endif         
         endif
         endif
   enddo
   enddo
      num_amp1(I_capital) = num
   enddo          
 endif


!--------------------------------------------------------------------
! output the debug information
!--------------------------------------------------------------------
!    write(*,*) 'information about the <F_I | Eji | wavefunction> '
!	write(*,*) 'number of the slater configuration is : ', n_string_ab
!	write(*,*)
!	do idicp =1,n_string_ab 
     
!	   write(*,*) 'there are ', num_amp1(idicp), ' term /= 0, for the--', idicp,'-- slater term'
!       write(*,*)'------------------------------------------------------------------------------'
!       write(*,*)    '        Eji-i', '      Eji-j', '      postion ', '      sign '    
!	   do jdicp =1,num_amp1(idicp)
!         write(*,*)
!		 write(*,*) i_amp1(idicp,jdicp,1),i_amp1(idicp,jdicp,2),i_amp1(idicp,jdicp,3),i_amp1(idicp,jdicp,4)
!		 write(*,*) 

!	   
!	   enddo
!	     write(*,*)
!	     write(*,*)
!   enddo


  
   return
  end subroutine amplitude1

!------------------------------------------------------------------------------------------------------
!find the value :
!
!  
!
! from haru
!------------------------------------------------------------------------------------------------------
 subroutine amplitude2
  implicit none
  integer(4) :: iab,ipqrs,ire,ipe,iir,iip,is,iq,iis,iiq,ii,jj,ir,ip,ipm,kk
  integer(4),dimension(system%numelectron) :: ibox
  integer(4) :: idicp,iresult
! from haru

  allocate(num_amp2(n_string_ab))
  allocate(i_amp2(n_string_ab,100000,10))  !need to revised here
  
  do iab=1,n_string_ab
     ipqrs=0   ! counter reset
     do ire=2,system%numelectron
     do ipe=1,ire-1
        ir=i_string_ab(iab,ire)
        ip=i_string_ab(iab,ipe)  ! ip < ir
        iir=spino(ir) !(ir-1)/ras1%nptot
        iip=spino(ip)  !(ip-1)/ras1%nptot

     do is=2,system%nptot*2
     do iq=1,is-1               ! iq < is
        iis=spino(is)     !(is-1)/ras1%nptot
        iiq=spino(iq)     !(iq-1)/ras1%nptot

!----judge-------------------
     if (iip+iir.eq.iis+iiq) then

        do ii=1,system%numelectron
           jj=i_string_ab(iab,ii)
           if (jj.eq.ip) then
              jj=iq
           elseif (jj.eq.ir) then
              jj=is
           endif
           ibox(ii)=jj
        enddo


         call  findindex(ibox,iresult,kk)

         if(kk/=0) then
          call  getsign(i_string_ab(iresult,:),ibox,system%numelectron,ipm)

!=================================================================
           ipqrs=ipqrs+1        ! counter increment
           i_amp2(iab,ipqrs,1)=ip-system%nptot*iip
           i_amp2(iab,ipqrs,2)=iq-system%nptot*iiq        ! iq <- ip
           i_amp2(iab,ipqrs,3)=ir-system%nptot*iir
           i_amp2(iab,ipqrs,4)=is-system%nptot*iis        ! is <- ir
!           i_amp2(iab,ipqrs,5)=ipm*iresult
           i_amp2(iab,ipqrs,5)=iresult
           i_amp2(iab,ipqrs,6)=iip
           i_amp2(iab,ipqrs,7)=iiq
           i_amp2(iab,ipqrs,8)=iir
           i_amp2(iab,ipqrs,9)=iis 
           i_amp2(iab,ipqrs,10) = ipm 
!=================================================================
        endif

     endif

     enddo
     enddo
     enddo
     enddo
     num_amp2(iab)=ipqrs

  enddo


!
! for debug output
!
!
!   write(*,*) '______________two electron replacement term______________________  '

!  do iab = 1,n_string_ab
!    write(*,*)
!	write(*,*)
!    write(*,*) 'for slater ',iab,' number of term is ', num_amp2(iab)
!	write(*,*)'                          -ip-', '        -iq-', '        -ir-', '       -is-', '   ipm*iresult'
!	do idicp =1, num_amp2(iab)
!     ip = i_amp2(iab,idicp,1)
!	 iq = i_amp2(iab,idicp,2)
!	 ir = i_amp2(iab,idicp,3)
!	 is = i_amp2(iab,idicp,4)
!     ipm = i_amp2(iab,idicp,5)
!    write(*,*) idicp,'    ', ip,iq,ir,is,ipm
!	enddo
!  enddo



  return
 end subroutine amplitude2

!
! find the position of iarray in jout 
!  haru 
!
    subroutine findindex(iarray,iresult,kk)
     implicit none
     integer,intent(out) :: iresult
     integer :: idicp,jdicp,numdiffer
     integer,intent(in),dimension(system%numelectron) :: iarray
     integer,dimension(system%numelectron):: itemp
     integer :: kk,judge0

     kk = 1
     judge0 = 0
     do idicp = 1,system%numelectron-1
         do jdicp =idicp+1,system%numelectron
          if(iarray(idicp) == iarray(jdicp)) then
           judge0 = 1
           kk =0
          endif
         enddo
     enddo

    if(judge0==0) then
    do idicp=1,n_string_ab
         itemp = i_string_ab(idicp,:)
         call get(iarray,itemp,system%numelectron,numdiffer)
         iresult=-1
   
      if(numdiffer==0) then
       iresult=idicp
       exit
      endif
     enddo

    if(iresult==-1) then
      kk = 0
    endif
   endif
  return
 end subroutine findindex


!====================================================================================
! set alphastring and betastring 
! set the weigth of vertex , the details can be found in the paper of Jeppe, Olsen
! j.chem.phys. 89, 1988, or the book of W. Duch, GRMS or Graphical represnetation of 
! model spaces, Lecture notes in Chemistry, Section 1.1, introducing graphical repres-
! entation 
!===================================================================================
 subroutine set_alphastring()
 implicit none
  integer :: ie,io
  integer :: iab,ii_add
  integer :: ibox(system%numelectron)

  allocate(w_vertex(0:system%nptot,0:system%numelectron/2)) 
  allocate(n_string_ab_inv(n_string_a_fci**2))

!
! initialization of the array
! this array store the weight of every vertex, have important relation with the single occupied src
!
  w_vertex = 0 

!
! this array store the index of the path (every path is corrsponding to one configuration)
!
  n_string_ab_inv = 0   

!
! set the value of the first column
!
  do io=0,system%nptot-system%numelectron/2
     w_vertex(io,0)=1
  enddo
!
!set the value of the rightmost position
!
  do ie=0,system%numelectron/2
     w_vertex(ie,ie)=1
  enddo

  w_vertex(1,0)=1
!
! the relation between (e,o) with (e-1,o-1) and (e,o-1)
!
  do ie=1,system%numelectron/2
  do io=1,system%nptot
     w_vertex(io,ie)=w_vertex(io-1,ie)+w_vertex(io-1,ie-1)
  enddo
  enddo

!
!find the index of the path and store 
!
  do iab=1,n_string_ab
   
    ibox = i_string_ab(iab,:)

!
! input : ibox store the spin orbital (configuration)
!
    call address(ii_add,ibox)
    
!
! output: ii_add, is the index of the path(configuration)
!

   n_string_ab_inv(ii_add)=iab
  enddo 

 return
end subroutine set_alphastring

!==================================================================
subroutine amp_alg()
!==================================================================
  use global

  implicit none
  integer :: iab,jab,ipq,ipqrs,ip,iq,ir,is,iip,iiq,iir,iis
  integer :: ipe,ire,ii,jj,kk,ipm,ibox(system%numelectron)
  integer :: l_amp1,l_amp2, ispin_ip, ispin_iq
 
    l_amp1 = system%numelectron/2*(system%nptot-system%numelectron/2+1)*2
     allocate(num_amp1(n_string_ab))
     allocate(i_amp1(n_string_ab,l_amp1,4))  !need to revised here

    l_amp2 = system%numelectron/2*(system%numelectron/2-1)*(system%nptot-system%numelectron/2+2)*&
    (system%nptot-system%numelectron/2+1)/4*2+system%numelectron/2*system%numelectron/2*&
    (system%nptot-system%numelectron/2+1)*(system%nptot-system%numelectron/2+1)

    allocate(num_amp2(n_string_ab)) 
    allocate(i_amp2(n_string_ab,l_amp2,10))

  do iab=1,n_string_ab
     ipq=0  
     do ipe=1,system%numelectron
     ip=i_string_ab(iab,ipe)
  !
  !wenliang, ip is the serial number of spin orbital, ipp is the spin of this spin orbital
  ! for example, if ipp =0, positive, ipp =1, nagative
  !    
     ispin_ip = spino(ip) ! remember the spin of the orbital ip  0, or 1
      
     do iq=1,system%nptot*2
       ispin_iq = spino(iq) 
!        
!----judge-----------------------
     if (ispin_ip.eq.ispin_iq) then !wenliang spin = spin, in the same spin space 
        do ii=1,system%numelectron
        jj=i_string_ab(iab,ii)
        if (jj.eq.ip) jj=iq
        ibox(ii)=jj
        enddo                        !wenliang 2 1 5 6  ! np_tot = 3, nel =4, nel/2 = 2
                                     !   
        call asc_order(ibox,ipm,kk)  ! 
                                     ! 1 2 5 6  ipm is coefficient,  
        if (kk.ne.0) then  !wenliang kk=0, there are two electron in the same spin orbital, not allowed
!=================================================================
           ipq=ipq+1         ! counter increment
           i_amp1(iab,ipq,1)=ip-system%nptot*ispin_ip  !记录哪一个空间轨道
           i_amp1(iab,ipq,2)=iq-system%nptot*ispin_iq        ! iq <- ip
           i_amp1(iab,ipq,3)=kk
           i_amp1(iab,ipq,4)=ipm
!=================================================================
        endif
     endif
     
     enddo
     enddo
     if (l_amp1.lt.ipq) stop  '***** ERROR IN AMP_ALG 1 *****'
     num_amp1(iab)=ipq
  enddo
!!  !$omp end do nowait
!----two-electron replacement------------------------------------------
!!  !$omp do
  do iab=1,n_string_ab
     ipqrs=0   ! counter reset
     do ire=2,system%numelectron
     do ipe=1,ire-1
        ir=i_string_ab(iab,ire)
        ip=i_string_ab(iab,ipe)  ! ip < ir
        iir= spino(ir)
        iip= spino(ip)
     do is=2,system%nptot*2
     do iq=1,is-1               ! iq < is
        iis=spino(is)
        iiq=spino(iq)
!----judge-------------------
     if (iip+iir.eq.iis+iiq) then
        do ii=1,system%numelectron
           jj=i_string_ab(iab,ii)
           if (jj.eq.ip) then
              jj=iq
           elseif (jj.eq.ir) then
              jj=is
           endif
           ibox(ii)=jj
        enddo
        call asc_order(ibox,ipm,kk)
        if (kk.ne.0) then
!=================================================================
           ipqrs=ipqrs+1        ! counter increment
           i_amp2(iab,ipqrs,1) = ip-system%nptot*iip
           i_amp2(iab,ipqrs,2) = iq-system%nptot*iiq        ! iq <- ip
           i_amp2(iab,ipqrs,3) = ir-system%nptot*iir
           i_amp2(iab,ipqrs,4) = is-system%nptot*iis        ! is <- ir
           i_amp2(iab,ipqrs,5) = kk
           i_amp2(iab,ipqrs,6) = iip
           i_amp2(iab,ipqrs,7) = iiq
           i_amp2(iab,ipqrs,8) = iir
           i_amp2(iab,ipqrs,9) = iis
           i_amp2(iab,ipqrs,10) = ipm
!=================================================================
        endif
     endif
     enddo
     enddo
     enddo
     enddo
     if (l_amp2.lt.ipqrs) stop  '***** ERROR IN AMP_ALG 2 *****'
     num_amp2(iab)=ipqrs
  enddo
  
  return
end subroutine amp_alg

!
! construct the zet_alg1 parameter
!
!=================================================================
subroutine zet_alg1
!=================================================================
  implicit none

  integer :: jab,jab_h,judge,judge1,judge2,jab_ij,ipm,iab_ij,iz,iab
  integer :: iel,jel,i,ii,iii,j,jj,jjj,ie,ibox(system%numelectron),NNA1,NNA2

  allocate(i_string_ab_h(n_string_a_fci**2,system%numelectron+1))
  allocate(i_string_ab2_h(n_string_a_fci**2,2*system%nptot+1)) 
 

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



!==================================================================
subroutine zet_alg2
!==================================================================
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

!============================================================================================================
! input: ibox(1:numelectron) is array, store the spin orbital, (A slater Derter.)
! ouput: kk: if kk ==0, this slater are not exist, or it is not a real slater(not satifsy the pauli principle)
!            if kk /=0, kk is the index of the slater, 
!            if kk /=0, ipm store the sign, + / -  
!============================================================================================================
subroutine asc_order(ibox,i_sign,i_index)

!================================================================
  implicit none
  integer,intent(out) :: i_sign,i_index
  integer  :: ibox(system%numelectron)
  integer :: judge0,ie,iep,ii_add,nn
!
! if two electron occupied the same orbital, not real slater, return,kk=0 
!    
  judge0=0
  do ie=2,system%numelectron
    do iep=1,ie-1
      if (ibox(ie)==ibox(iep)) judge0=1
    enddo
  enddo

!----------------------------------------------------------------------------
! if it is a real slater, then we have to try find its position (index)
!
  if (judge0.eq.0) then !-----------------------------------------------------
!----------------------------------------------------------------------------
!
  i_sign=1
!
!exchange the sequence of the orbital, we need ordered the orbitals 
!from small one to large one, the ipm remember the exchange time
! 
  do ie=1,system%numelectron-1
   do iep=ie+1,system%numelectron
     if (ibox(ie).gt.ibox(iep)) then
        nn=ibox(ie)
        ibox(ie)=ibox(iep)
        ibox(iep)=nn
        i_sign=i_sign*(-1)
     endif
   enddo
  enddo
!
!find the index of the slater derter.
!
  call address(ii_add,ibox)
  i_index = n_string_ab_inv(ii_add)

!----------------------------------------------------------------------------
  else !---------------------------------------------------------------------
!----------------------------------------------------------------------------
  i_index = 0
!----------------------------------------------------------------------------
  endif !--------------------------------------------------------------------
!----------------------------------------------------------------------------
  return
end subroutine asc_order

!=======================================================================
! find the index of the path by using alphastring and betastring method
! the details can be found in the paper, j.chem.phys. 89(4),1988
! Jeppe Olsen etc.
!-----------------------------------------------------------------------
! input, ibox(nel), nel is the number of the electron,
!        ibox store the index of the spin orbitals
! output: ii_add is the index of the path
!    0     1    2
! 0  (1)
!     |  \
!     |   \           
! 1  (1)  (1)
!     | \  | \ 
!     |  \ |  \ 
! 2  (1)  (2)  (1)
!       \  | \  |
!        \ |  \ |
! 3       (3)  (3)   
!           \   |
!            \  |
! 4           (6)  
!
! the w_vertex is the weight of the vertex, please read the details in 
! the paper, the above example: ne = 2, nptot = 4
!=======================================================================
subroutine address(ii_add,ibox)
!=======================================================================
  implicit none
  integer :: ibox(system%numelectron),ibox_2nd(system%nptot*2)
  integer :: ii,mm,km,ii_a,ii_b,ii_add,ie
!
! ibox_2nd is the second quantum representation of the slater ibox
! for example, ibox(4) = /1,2,4,5/
! then, ibox_2nd = /1,1,0,1,1/
!
  do mm=1,system%nptot*2
     ii=0
     do ie=1,system%numelectron
        if (mm.eq.ibox(ie)) ii=1
     enddo
     ibox_2nd(mm)=ii
  enddo
!
! find the index of alphastring path
!
  ii_a=0
  do mm=1,system%nptot
     km=0
     do ii=1,mm
        km=km+ibox_2nd(ii)
     enddo
     ii_a=ii_a+ibox_2nd(mm)*W_VERTEX(mm-1,km)
  enddo
!
! find the index of betastring path
!
  ii_b=0
  do mm=1,system%nptot
     km=0
     do ii=1,mm
        km=km+ibox_2nd(ii+system%nptot)
     enddo
     ii_b=ii_b+ibox_2nd(mm+system%nptot)*W_VERTEX(mm-1,km)
  enddo
!
! find the total slater index
!
  ii_add=1+ii_a+ii_b*N_STRING_A_FCI
  return
end subroutine address




 end module auxiliary


 subroutine drive_auxiliary
  use global
  use auxiliary
  implicit none
!
! i_string_a, i_string_a2 
!

  call  config1()
!
! spino() and space()
!

  call space_plus_spino()

!
! n_string_ab, n_string_ab2,n_string_ab
!

  call config2()
!
! set alpha string
!

  call set_alphastring()


!  call amplitude1()
!  call amplitude2()


  call amp_alg()


!
! taken from haru's method
!
  if(system%io_sdtq/=0) then

   write(*,*) '----------modified density matrices------------------- '
   n_aten = system%np1*system%np2

   allocate(i_string_ab_hij(n_string_a_fci**2*system%nptot**2,system%numelectron+3))
   allocate(czet1(n_aten,n_aten))
   allocate(cv_zet2(n_aten)) 

   call zet_alg1()

  l_zmp1 = n_string_ab_hij*(system%io_sdtq+1)*system%np1
  l_zmp2 = n_string_ab_hij*(system%io_sdtq+1)*(system%numelectron-1)*(system%np1+system%np2)

  allocate(i_zmp1(l_zmp1,6))
  allocate(i_zmp2(l_zmp2,8))
 


  call zet_alg2()
 endif

 return
end subroutine drive_auxiliary

