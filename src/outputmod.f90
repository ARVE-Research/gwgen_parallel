module outputmod

use parametersmod, only : i2,i4,sp

implicit none

public :: genoutfile
public :: putlonlat

real(sp),    parameter :: missing_sp = -9999.
integer(i2), parameter :: missing_i2 = -32768

contains

!-------------------------------------------------------------------------------------------------

subroutine genoutfile(outfile,id,chunks,ofid)

use parametersmod, only : i4,sp,dp
use netcdf
use errormod,     only : ncstat,netcdf_err
use coordsmod,    only : index

implicit none

character(*),           intent(in)  :: outfile  ! file name
type(index),            intent(in)  :: id
integer, dimension(:),  intent(in)  :: chunks
integer,                intent(out) :: ofid

!local variables

integer(i4) :: dimid
integer(i4) :: varid

real(dp), dimension(2) :: xrange
real(dp), dimension(2) :: yrange

character(8)  :: today
character(10) :: now

integer(i4), allocatable, dimension(:) :: dimids

integer, parameter :: ndims = 3

!---

xrange = [id%minlon,id%maxlon]
yrange = [id%minlat,id%maxlat]

!---

!write(0,'(a,a)')'creating output file: ',trim(outfile)
!write(0,*)'create',id%countx,id%county

ncstat = nf90_create(outfile,nf90_hdf5,ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','weathergen output file')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'created',today//' '//now(1:4))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!write(0,*)'added global atts'

!-----------
!dimensions

allocate(dimids(ndims))

!----
! lon

ncstat = nf90_def_dim(ofid,'lon',id%countx,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(1) = dimid

ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! lat

ncstat = nf90_def_dim(ofid,'lat',id%county,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(2) = dimid

ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! time

ncstat = nf90_def_dim(ofid,'time',nf90_unlimited,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(3) = dimid

ncstat = nf90_def_var(ofid,'time',nf90_float,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','time')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','days since 1950-01-01')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!absolute minimum temperature (lon,lat)

ncstat = nf90_def_var(ofid,'abs_tmin',nf90_short,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','absolute minimum temperature')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degC')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
!absolute maximum temperature (lon,lat)

ncstat = nf90_def_var(ofid,'abs_tmax',nf90_short,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','absolute maximum temperature')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degC')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_enddef(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

end subroutine genoutfile

!-------------------------------------------------------------------------------------------------

subroutine putlonlat(ofid,id,lon,lat)

use parametersmod, only : i4,sp,dp
use netcdf
use coordsmod,     only : index
use errormod,      only : ncstat,netcdf_err

implicit none

integer,                intent(in) :: ofid
type(index),            intent(in) :: id
real(dp), dimension(:), intent(in) :: lon
real(dp), dimension(:), intent(in) :: lat

integer :: varid

!---

ncstat = nf90_inq_varid(ofid,'lon',varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lon)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lat)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

end subroutine putlonlat

!-------------------------------------------------------------------------------------------------

end module outputmod
