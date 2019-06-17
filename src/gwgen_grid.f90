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

character(100) :: infile									! input file name

integer :: xlen												! length of dimension 'lat'
integer :: ylen												! length of dimension 'long'
integer :: tlen												! length of dimension 'time'

integer :: ifid
integer :: dimid
integer :: varid

real(dp), allocatable, dimension(:) :: lon					! Make allocatable array for longitude
real(dp), allocatable, dimension(:) :: lat					! Make allocatable array for latitude

integer, dimension(2) :: xpos
integer, dimension(2) :: ypos

integer :: srtx												! Start value x for boundary box (upper left??????)
integer :: srty												! Start value y for boundary box
integer :: cntx												! Amount of longitude cells
integer :: cnty												! Amount of latitude cells

integer(i2), allocatable, dimension(:,:,:) :: var_in  		! temporary array for input data in i2 format

real(sp), allocatable, dimension(:,:,:) :: tmp  			! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: dtr  			! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:,:) :: pre  			! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:,:) :: wet  			! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:,:) :: cld  			! mean monthly cloud cover (percent)
real(sp), allocatable, dimension(:,:,:) :: wnd  			! mean monthly 10m windspeed (m s-1)

real(sp) :: scale_factor
real(sp) :: add_offset
integer(i2) :: missing_value

integer :: i,y,m
integer :: nyrs												! Number of years (tlen/12)
integer :: nmos


!----------------------------------------------------
! Read dimension IDs and lengths of dimensions

call getarg(1,infile)
	
ncstat = nf90_open(infile,nf90_nowrite,ifid)				! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)			! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)					! get dimension ID from dimension 'lon' in the input file (file id, dimension name, dimension ID)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)		! get dimension name and length from input file for dimension previously inquired (file id, dimension ID, length will be written in variable xlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)					! Get dimension ID from dimension 'lat' 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)		! Get length of dimension 'lat' and assign it to variable ylen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)					! Get dimension ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)		! Get length of dimension 'time' and assign it to variable tlen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

nyrs = tlen / 12											! Number of years = months / 12

write(so,*)xlen,ylen,nyrs


!----------------------------------------------------
! Read variable IDs and values 

allocate(lon(xlen))											! Allocate length to longitude array
allocate(lat(ylen))											! Allocate length to latitude array

ncstat = nf90_inq_varid(ifid,"lon",varid)					! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
	
ncstat = nf90_get_var(ifid,varid,lon)						! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)					! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)			

ncstat = nf90_get_var(ifid,varid,lat)						! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


!----------------------------------------------------
! In terminal, call the programs coordstring and parsecoords to determine boundaries of area of interest (translates lat/long values into indices of the lat/long arrays)

call getarg(2,coordstring)									

call parsecoords(coordstring,bounds)

write(so,*)bounds											! Print boundaries of area of interest

call calcpixels(lon,lat,bounds,xpos,ypos,srtx,srty,cntx,cnty)

write(so,*)'xpos: ',xpos									! Print start position of longitude
write(so,*)'ypos: ',ypos									! Print start position of latitude
write(so,*)'cntx, cnty: ',cntx,cnty							! Print counts of x and y cells	

allocate(var_in(cntx,cnty,tlen))							! Allocate space in input array 'var_in' (x-range, y-range and temporal range)


!---------------------------------------------------------------------
! get the temperature array timeseries

allocate(tmp(cntx,cnty,tlen))								! Allocate space of area of interest and temporal range to tmp array (mean monthly temperature)

tmp = -9999.												! Set tmp to -9999

ncstat = nf90_inq_varid(ifid,"tmp",varid)					! Get variable ID of variable tmp 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])		! Get values for variable tmp from input file, starting at the starting point and going for cnt x and y cells
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)							! Get attribute 'missing value' in the variable temperature
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)							! Get attribute 'scale factor' 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)								! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) tmp = real(var_in) * scale_factor + add_offset			! Where the temperature attribute is not missing value, calculate the real temperature using the scale factor and the add_offset


!---------------------------------------------------------------------
! get the diurnal temperature range array timeseries

allocate(dtr(cntx,cnty,tlen))								! Allocate space to array dtr 

dtr = -9999.												! Set dtr to -9999								

ncstat = nf90_inq_varid(ifid,"dtr",varid)					! Get variable ID of variable dtr
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])		! Get values for variable dtr from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)							! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)							! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)								! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) dtr = real(var_in) * scale_factor + add_offset			! Where dtr is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the precipitation array timeseries

allocate(pre(cntx,cnty,tlen))															! Allocate space to array pre (precipitation) 

pre = -9999.																			! Set pre to -9999								

ncstat = nf90_inq_varid(ifid,"pre",varid)												! Get variable ID of variable pre
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])		! Get values for variable pre from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)							! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)							! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)								! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) pre = real(var_in) * scale_factor + add_offset			! Where pre is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the wet days array timeseries

allocate(wet(cntx,cnty,tlen))															! Allocate space to array wet (precipitation) 

wet = -9999.																			! Set wet to -9999								

ncstat = nf90_inq_varid(ifid,"wet",varid)												! Get variable ID of variable wet
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])		! Get values for variable wet from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)							! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)							! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)								! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) wet = real(var_in) * scale_factor + add_offset			! Where wet is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the cloud cover array timeseries

allocate(cld(cntx,cnty,tlen))															! Allocate space to array cld (precipitation) 

cld = -9999.																			! Set cld to -9999								

ncstat = nf90_inq_varid(ifid,"cld",varid)												! Get variable ID of variable cld
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])		! Get values for variable cld from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)							! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)							! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)								! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) cld = real(var_in) * scale_factor + add_offset			! Where cld is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the wind speed array timeseries

allocate(wnd(cntx,cnty,tlen))															! Allocate space to array wnd (precipitation) 

wnd = -9999.																			! Set wnd to -9999								

ncstat = nf90_inq_varid(ifid,"wnd",varid)												! Get variable ID of variable wnd
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])		! Get values for variable wnd from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)							! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)							! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)								! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) wnd = real(var_in) * scale_factor + add_offset			! Where wnd is not missing value, calculate real values using add_offset and scale_factor

!---------------------------------------------------------------------
! close the input file

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


!---------------------------------------------------------------------
! data check

write(so,*)'Year, month, tmin, tmean, tmax, precip, wet days, cloud, wind:'

i = 1
do y = 1,nyrs
  do m = 1,12
  
    write(so,'(2i5, 7f7.2)')y,m,&
    tmp(1,1,i) - 0.5 * dtr(1,1,i),tmp(1,1,i),tmp(1,1,i) + 0.5 * dtr(1,1,i),& 
    pre(1,1,i), wet(1,1,i), cld(1,1,i)/100, wnd(1,1,i)
    
    i = i + 1

  end do
end do

!---------------------------------------------------------------------

! '(2i5,3f7.1)'



end program gwgen_grid