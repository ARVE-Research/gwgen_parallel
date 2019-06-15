program gwgen_grid

! use the Makefile to compile this program

! program to run gwgen with gridded input, provide the name of a climate data input file and geographic bounds for the simulation using a xmin/xmax/ymin/ymax string
! JO Kaplan, HKU, 2019

use parametersmod, only : sp,dp,i4,i2,so
use errormod, only : ncstat,netcdf_err
use coordsmod, only : coordstring,bounds,parsecoords,calcpixels
use netcdf

implicit none

! inquire about the dimensions of the input file

character(100) :: infile

integer :: xlen
integer :: ylen
integer :: tlen

integer :: ifid
integer :: dimid
integer :: varid

real(dp), allocatable, dimension(:) :: lon
real(dp), allocatable, dimension(:) :: lat

integer, dimension(2) :: xpos
integer, dimension(2) :: ypos

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty

integer(i2), allocatable, dimension(:,:,:) :: var_in  ! temporary array for input data in i2 format

real(sp), allocatable, dimension(:,:,:) :: tmp  ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: dtr  ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:,:) :: pre  ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:,:) :: wet  ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:,:) :: cld  ! mean monthly cloud cover (percent)
real(sp), allocatable, dimension(:,:,:) :: wnd  ! mean monthly 10m windspeed (m s-1)

real(sp) :: scale_factor
real(sp) :: add_offset
integer(i2) :: missing_value

integer :: i

!----------------------------------------------------

call getarg(1,infile)

ncstat = nf90_open(infile,nf90_nowrite,ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


write(so,*)xlen,ylen

allocate(lon(xlen))
allocate(lat(ylen))

ncstat = nf90_inq_varid(ifid,"lon",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lon)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lat)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

call getarg(2,coordstring)

call parsecoords(coordstring,bounds)

write(so,*)bounds

call calcpixels(lon,lat,bounds,xpos,ypos,srtx,srty,cntx,cnty)

write(so,*)xpos
write(so,*)ypos
write(so,*)cntx,cnty

allocate(var_in(cntx,cnty,tlen))

!---------------------------------------------------------------------
! get the temperature array timeseries

allocate(tmp(cntx,cnty,tlen))

tmp = -9999.

ncstat = nf90_inq_varid(ifid,"tmp",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) tmp = real(var_in) * scale_factor + add_offset

do i = 1,tlen
  write(so,*)i,tmp(1,1,i)
end do

!---------------------------------------------------------------------
! get the diurnal temperature range array timeseries

allocate(dtr(cntx,cnty,tlen))

tmp = -9999.

ncstat = nf90_inq_varid(ifid,"dtr",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) dtr = real(var_in) * scale_factor + add_offset

do i = 1,tlen
  write(so,*)i,tmp(1,1,i)
end do

!---------------------------------------------------------------------
! get the pre array timeseries

! ...








ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end program gwgen_grid