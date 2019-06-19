module parametersmod

! Simple module defining some types and parameters

use iso_fortran_env, only : int8,int16,int32,real32,real64,output_unit

implicit none

integer, parameter :: dp = real64       ! 8 byte real
integer, parameter :: sp = real32       ! 4 byte real
integer, parameter :: i2 = int16        ! 2 byte integer
integer, parameter :: i4 = int32        ! 4 byte integer

integer, parameter :: so = output_unit  ! unit number for standard output

real(sp), parameter :: Tfreeze = 273.15 ! freezing temperature of freshwater (K)

end module parametersmod
