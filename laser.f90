!
! laser field: electric field
!
 subroutine laserpara(time)
   use global
   implicit none
   real(kind=k1) :: time
   real(kind=k1):: tcycle,t_duration
        
!   q(1)=ex_max    !0    !Ex
!   q(2)=ey_max    !0    !Ey
!   q(3)=ez_max    !0.01d0 !Ez
!   omega= laser1%omega         !5.d0 ! omega
    !/2.418884d-2  ! fs-->as
    tcycle=2.d0*pi/laser%omega
!    t_duration=6.d0*tcycle
    t_duration=laser%ncycle*tcycle  !! changed by Juan

    If (time.ge.t_duration/2) then !! the pulse is finished
       laser%ex_t=zero
       laser%ey_t=zero
       laser%ez_t=zero
       return
    End If

!! This is the pulse   
!!$   laser%ex_t = (dsin(pi*time/t_duration))**2.d0*laser%ex_max*dcos(laser%omega*time)!(dsin(pi*time/t_duration))**2.d0*laser1%ex_max*dcos(omega*time)
!!$   laser%ey_t = (dsin(pi*time/t_duration))**2.d0*laser%ey_max*dcos(laser%omega*time)!(dsin(pi*time/t_duration))**2.d0*laser1%ey_max*dcos(omega*time)
!!$   laser%ez_t =(dsin(pi*time/t_duration))**2.d0*laser%ez_max*dcos(laser%omega*time)


!!$    laser%ey_t = laser%ey_max*(sin(2.0d0*pi*time/t_duration)*sin(laser%omega*time+laser%cep)/(2.0d0*dble(laser%ncycle))-cos(pi*time/t_duration)**2.0d0*cos(laser%omega*time+laser%cep))
!!$    laser%ez_t = laser%ez_max*(sin(2.0d0*pi*time/t_duration)*sin(laser%omega*time+laser%cep)/(2.0d0*dble(laser%ncycle))-cos(pi*time/t_duration)**2.0d0*cos(laser%omega*time+laser%cep))

     laser%ex_t = laser%e_max*laser%ellipticity/sqrt(1.0d0+laser%ellipticity**2.0d0)*(sin(2.0d0*pi*time/t_duration)*cos(laser%omega*time+laser%cep)/(2.0d0*dble(laser%ncycle))+cos(pi*time/t_duration)**2.0d0*sin(laser%omega*time+laser%cep))   
     laser%ez_t = laser%e_max/sqrt(1.0d0+laser%ellipticity**2.0d0)*(sin(2.0d0*pi*time/t_duration)*sin(laser%omega*time+laser%cep)/(2.0d0*dble(laser%ncycle))-cos(pi*time/t_duration)**2.0d0*cos(laser%omega*time+laser%cep))   

   return
  end subroutine laserpara

!
! laser field: vector potential
!
 subroutine laserpara_a(time)
   use global
   implicit none
   real(kind=k1) :: time
   real(kind=k1):: tcycle,t_duration
        
    tcycle=2.d0*pi/laser%omega
    t_duration=laser%ncycle*tcycle

    If (time.ge.t_duration/2) then !! the pulse is finished
       laser%ax_t=zero
       laser%ay_t=zero
       laser%az_t=zero
       return
    End If

     laser%ax_t = (laser%e_max/laser%omega)*laser%ellipticity/sqrt(1.0d0+laser%ellipticity**2.0d0)*cos(pi*time/t_duration)**2.0d0*cos(laser%omega*time+laser%cep) !!A_x(t)
     laser%az_t = (laser%e_max/laser%omega)/sqrt(1.0d0+laser%ellipticity**2.0d0)*cos(pi*time/t_duration)**2.0d0*sin(laser%omega*time+laser%cep) !! A_z(t)

   return
 end subroutine laserpara_a


!
!
! linear laser field
! -->        2                                     ----->
!  E  ==  Sin  ( pi*time / time ) * Cos(w*t + fi)  ex_y_z  
!

subroutine laser_linear_x(time)
 use global
 implicit none
 real(kind=k1) :: cep
 real(kind=k1) :: time,tcycle,t_duration 

 cep = 0.0d0
 tcycle = 2.0d0*pi/laser%omega
 t_duration = laser%ncycle * tcycle

 laser%ex_t = laser%ex_max*dsin(pi*time/t_duration)**2*&
dcos(laser%omega*time + cep)   

 return
end subroutine laser_linear_x

!
! linear laser field along with y direction
!
subroutine laser_linear_y(time)
 use global
 implicit none
 real(kind=k1) :: cep
 real(kind=k1) :: time,tcycle,t_duration
 cep = 0.0d0
 tcycle = 2.0d0*pi/laser%omega
 t_duration = laser%ncycle * tcycle

 laser%ey_t = laser%ey_max*dsin(pi*time/t_duration)**2*&
dcos(laser%omega*time + cep)

 return
end subroutine laser_linear_y

!
! linear laser filed alon with z direction
!
subroutine laser_linear_z(time)
 use global
 implicit none
 real(kind=k1) :: time,tcycle,t_duration
 real(kind=k1) :: cep

 cep = 0.0d0
 tcycle = 2.0d0*pi/laser%omega
 t_duration = laser%ncycle * tcycle

 laser%ez_t = laser%ez_max*dsin(pi*time/t_duration)**2*&
dcos(laser%omega*time + cep)
 return
end subroutine laser_linear_z




!
! laser pulse  ellipitically polarized in the x - y plance
!
! -->      E_max      2        pi * t                          -->                         -->
!  E  ==  -------- Sin  (---------------- ) * { Cos(w*t + fi)*  ex   + e * Sin(w*t + fi) *  ey   }
!          ______          ncycle*T_cycle
!         /     2 
!       \/ 1 + e  
!       
!  
! e: ellipticity
! 
! ncycle: the number of the cycle
!
! T_cycle : the 
! 
! fi: carrier envolpe phase
!
subroutine laser_ellipticity_xy_x(time)
 use global
 implicit none
 real(kind=k1) :: time
 real(kind=k1) :: tcycle,t_duration
 real(kind=k1) :: cep

 cep = 0.0d0
 
 tcycle = 2.0d0*pi/laser%omega
 t_duration = laser%ncycle*tcycle

 laser%ex_t = laser%ex_max * sin(pi*time/t_duration)**2/(dsqrt(1.0 + laser%ellipticity**2 ))*&

 cos(laser%omega*time + cep) 

 return
end subroutine laser_ellipticity_xy_x




subroutine laser_ellipticity_xy_y(time)
 use global
 implicit none
 real(kind=k1) :: time
 real(kind=k1) :: tcycle,t_duration
 real(kind=k1) :: cep
 cep = 0.0d0

 tcycle = 2.0d0*pi/laser%omega
 t_duration = laser%ncycle*tcycle

 laser%ey_t = laser%ey_max * sin(pi*time/t_duration)**2* & 
 laser%ellipticity/(dsqrt(1.0 + laser%ellipticity**2 ))*sin(laser%omega*time + cep)

 return
end subroutine laser_ellipticity_xy_y


!
!
! two color laser pulse, pump + prob
!
!







