module parametersmod

! Simple module defining some types and parameters

use iso_fortran_env, only : int8,int16,int32,real32,real64,output_unit

implicit none

public :: ndaymonth

integer, parameter :: i1 = int8    ! 1 byte integer
integer, parameter :: i2 = int16   ! 2 byte integer
integer, parameter :: i4 = int32   ! 4 byte integer
integer, parameter :: sp = real32  ! 4 byte real
integer, parameter :: dp = real32  ! 8 byte real
! integer, parameter :: dp = real64  ! 8 byte real

integer, parameter :: so = output_unit  ! unit number for standard output

real(sp), parameter :: Tfreeze = 273.15 ! freezing temperature of freshwater (K)

! integer(i2), parameter :: baseyr = 1871

real(sp), parameter :: hsp = huge(sp)    ! largest positive 4-byte real

contains

!-----------------------------------------------------------------------------

integer function ndaymonth(yr,mon)

! Function to find out the number of days in a month, considering leap years and the year given as AD

! Input: Current year and month

integer, intent(in) :: yr
integer, intent(in) :: mon

! Arrays defining standard and leap years

integer, parameter, dimension(12) :: std_year = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
integer, parameter, dimension(12) :: leapyear = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

!---------------------------------

if (mod(yr,400) == 0) then        ! If year can be divided by 400, then it's a leap year

  ndaymonth = leapyear(mon)       ! Choose amount of days from leap year array

else if (mod(yr,100) == 0) then   ! If year can't be divided by 400 but can by 100, then it's a standard year

  ndaymonth = std_year(mon)

else if (mod(yr,4) == 0) then     ! If year can't be divided by 400 or 100, but can be divided by 4, it's a leap year

  ndaymonth = leapyear(mon)

else                              ! Any other case, it's a standard year

  ndaymonth = std_year(mon)

end if

end function ndaymonth

!-----------------------------------------------------------------------------

end module parametersmod
