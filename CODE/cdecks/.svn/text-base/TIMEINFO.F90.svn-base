!*****************************

MODULE TIMEINFO

! This module contains the variables describing the start and current time
!   of the model run.

! HISTORY:
!  2010.02.09 -DD- Converted to an f90 module from timeinfo.com

   USE KINDS
   
IMPLICIT NONE
PRIVATE

!*****************************
! Note: Fractional Julian day is specified as the day of year during 
!       iyr (0.0 = 00Z Jan 1)
!   formerly common/day1/

   REAL (KIND=dbl_kind), public ::            &
      rjday0,          &   ! Fractional Julian day at the start of the model simulation
      rjday,           &   ! Fractional Julian day of the current model time step
      utc_time             ! UTC time of the current model time step
      
!*****************************
!   formerly common/day2/

   INTEGER (KIND=int_kind), public ::         &
      iyr,             &   ! Year of current time step
      imonth,          &   ! Month of current time step
      iday                 ! Calendar day of current time step

!*****************************

END MODULE TIMEINFO