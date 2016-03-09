
!
! -----------------------------------------------------------------
! log information file
! -----------------------------------------------------------------
!
  subroutine log_information()
   use global
   implicit none
   character(LEN=50) :: exe_name,in_name,out_name
   logical :: alive
!
! initial the input_channel_file
!

   call getarg(0,exe_name)
   call getarg(1,in_name)

   inquire(file=in_name(1:len_trim(in_name)),exist=alive)

! check if the input file is exist or not   
    if(.not.alive) then
      write(*,*) trim(in_name), ' is absent'
      stop 'check the subroutine openfile'
    endif  
   
   open(unit=in_channel,file=in_name(1:len_trim(in_name)),status='old')
!
! write the log into the output_file
!  
   write(*,*) '#########################################################################'
   write(*,*) "# ------3d-td-ras-scf-calculation for atomic and molecular systems------#"
   write(*,*) '#                                                                       #'
   write(*,*) "#------basis set: spherical coordiante fe-dvr 3d------------------------# " 
   write(*,*) '#------basis set: prolate speridoial coordinate------------------------ #'
   write(*,*) "#-----------------------------------------------------------------------#"
   write(*,*) "#---------------All Rights Reserved-------------------------------------#"
   write(*,*) '#########################################################################'
   write(*,*)                                                        
   return
 end subroutine log_information

