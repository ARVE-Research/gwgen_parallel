module getdatamod

implicit none

public :: readdata

contains

!-----------------------------------------------

subroutine readdata(ncid,varname,srt,cnt,values)

use parametersmod, only : sp,i2
use errormod,      only : ncstat,netcdf_err
use netcdf

implicit none

integer,                    intent(in)  :: ncid
character(*),               intent(in)  :: varname
integer, dimension(3),      intent(in)  :: srt
integer, dimension(3),      intent(in)  :: cnt
real(sp), dimension(:,:,:), intent(out) :: values

!---

integer, dimension(3) :: s
integer(i2), allocatable, dimension(:,:,:) :: var_in

real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer :: i
integer :: varid

!---

s = [(size(values,dim=i),i=1,3)]

allocate(var_in(s(1),s(2),s(3)))

ncstat = nf90_inq_varid(ncid,varname,varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,var_in,start=srt,count=cnt)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,"add_offset",add_offset)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) 
  values = real(var_in) * scale_factor + add_offset
elsewhere
  values = -9999.
end where

end subroutine readdata

!-----------------------------------------------

end module getdatamod