!
! absorb potential , from right sides [  ]
!
    double precision function fun_absorb_rightside(xgrid)
     use global
        implicit none
        real(kind=k1),intent(in):: xgrid
	real(kind=k1) :: capot
      
     if(xgrid<=absorb%r_right)then
          capot=0.d0
        !elseif(absorb1%grid(i)>=absorb1%r_left-absorb1%ramp.and.absorb1%grid(i)<absorb1%r_left)then
        !  capot=absorb1%strength*(1.d0-dcos(abs(absorb1%grid(i)-absorb1%r_left)/absorb1%ramp*pi))/2.d0
        elseif(xgrid>absorb%r_right.and.xgrid<=absorb%r_right+absorb%ramp)then
          capot=absorb%strength*(1.d0-dcos(abs(xgrid-absorb%r_right)/absorb%ramp*pi))/2.d0
        elseif(xgrid>absorb%r_right+absorb%ramp)then
          capot=absorb%strength
        endif
         capot=-capot
       
        fun_absorb_rightside = capot
      return
     end function fun_absorb_rightside




!
! absorb potential , from two sides [  ]
!
    double precision function fun_absorb_twoside(xgrid)
      use global
        implicit none
        real(kind=k1),intent(in):: xgrid
	real(kind=k1) :: capot

         
     if(xgrid>=absorb%r_left.and.xgrid<=absorb%r_right)then
          capot=0.d0
        elseif(xgrid>=absorb%r_left-absorb%ramp.and.xgrid<absorb%r_left)then
          capot=absorb%strength*(1.d0-dcos(abs(xgrid-absorb%r_left)/absorb%ramp*pi))/2.d0
        elseif(xgrid>absorb%r_right.and.xgrid<=absorb%r_right+absorb%ramp)then
          capot=absorb%strength*(1.d0-dcos(abs(xgrid-absorb%r_right)/absorb%ramp*pi))/2.d0
        elseif(xgrid<absorb%r_left-absorb%ramp.or.xgrid>absorb%r_right+absorb%ramp)then
          capot=absorb%strength
        endif
         capot=-capot
       
        fun_absorb_twoside = capot

      return
   end function fun_absorb_twoside
