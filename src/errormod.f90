module errormod

implicit none

public :: netcdf_err

integer :: ncstat

contains

!Internal subroutines for error handline - checks error status after each call,
!prints out text message each time an error code is returned. 

!-------------------------

subroutine netcdf_err(ncstat)

use netcdf, only : nf90_strerror

implicit none

integer, intent(in) :: ncstat

write(0,'(a,i5,a,a)')' NetCDF error ',ncstat,' encountered: ',trim(nf90_strerror(ncstat))
stop

end subroutine netcdf_err

!-------------------------

end module errormod
