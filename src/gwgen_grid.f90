program gwgen_grid

! Use Makefile to compile this program

! Program to run gwgen with gridded input, provide the name of a climate data input file and geographic bounds for the simulation using a xmin/xmax/ymin/ymax string
! JO Kaplan, HKU, 2019

! Terminal command line: ./gwgen_grid ~/path/to/input.file x_coordinate/y_coordinate

use parametersmod, only : sp,dp,i4,i2,so
use errormod, only : ncstat,netcdf_err
use coordsmod, only : coordstring,bounds,parsecoords,calcpixels
use netcdf

implicit none

! Inquire about the dimensions of the input file

character(100) :: infile									! input file name

integer :: xlen												! length of dimension 'lat'
integer :: ylen												! length of dimension 'long'
integer :: tlen												! length of dimension 'time'

! IDs for file, dimensions and variables

integer :: ifid												! Input file ID
integer :: dimid											! Dimension ID
integer :: varid											! Variable ID


! Allocatable arrays for longitude and latitude

real(dp), allocatable, dimension(:) :: lon					
real(dp), allocatable, dimension(:) :: lat					
real(dp), allocatable, dimension(:) :: time

integer, dimension(2) :: xpos
integer, dimension(2) :: ypos


! Start values of x and y (LLC), and counts of cells in both directions: 

integer :: srtx												
integer :: srty												
integer :: cntx												
integer :: cnty												


! Array to store the input attributes

integer(i2), allocatable, dimension(:,:,:) :: var_in  		


! Monthly input attributes

real(sp), allocatable, dimension(:,:,:) :: tmp  			! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: dtr  			! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:,:) :: pre  			! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:,:) :: wet  			! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:,:) :: cld  			! mean monthly cloud cover (percent)
real(sp), allocatable, dimension(:,:,:) :: wnd  			! mean monthly 10m windspeed (m s-1)

! Monthly input attributes calculated here

real(sp), allocatable, dimension(:,:,:) :: mtmin  			! maximum monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: mtmax  			! monthly minimum temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: wetf				! fraction of wet days in a month

real(sp) :: scale_factor									! Value for the calculation of the "real" value of the parameters. Can be found in the netCDF file
real(sp) :: add_offset										! Value for the calculation of the "real" value of the parameters. Can be found in the netCDF file
integer(i2) :: missing_value								! Missing values in the input file


! Elements to calculate current year and amount of days in current month

integer :: i,y,m,j,x,t
integer :: nyrs												! Number of years (tlen/12)
integer, parameter :: nmos = 12								! Number of months

integer, parameter :: startyr = 1871						! Start year set to 1871
integer :: endyr											! End year (2010, not set yet)

integer :: yr												! Variable year 
integer :: mon												! Variable month 

integer, allocatable, dimension(:) :: nd	


! Elements for the smoothing process:

real, dimension(-1:1) :: mtmin3m							! Scalar for monthly minimum temperature
real, dimension(-1:1) :: mtmax3m							! Scalar for monthly maximum temperature
real, dimension(-1:1) :: nd3m								! Scalar for number of days in the month
real, dimension(-1:1) :: cld3m								! Scalar for monthly cloud fractions
real, dimension(-1:1) :: wnd3m								! Scalar for monthly wind speeds	

real(sp), dimension(2) :: bcond_tmin       					! boundary conditions of min temp for smoothing
real(sp), dimension(2) :: bcond_tmax       					! boundary conditions of max temp for smoothing
real(sp), dimension(2) :: bcond_cld 		   				! boundary conditions of cloud for smoothing
real(sp), dimension(2) :: bcond_wnd       					! boundary conditions of wind speed for smoothing				
real(sp), dimension(2) :: bcond_nd 		      				! boundary conditions of number of days for smoothing

real(sp), dimension(31) :: tmin_sm							! smoothed daily values of min temperature
real(sp), dimension(31) :: tmax_sm 							! smoothed daily values of max temperature
real(sp), dimension(31) :: nd_sm 							! smoothed daily values of max temperature
real(sp), dimension(31) :: cld_sm 							! smoothed daily values of cloudiness
real(sp), dimension(31) :: wnd_sm  							! smoothed daily values of wind speed

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

call getarg(1,infile)										! Reads first argument in the command line (path to input file)
	
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

allocate(lon(xlen))														! Allocate length to longitude array
allocate(lat(ylen))														! Allocate length to latitude array
allocate(time(tlen))													! Allocate length to latitude array

ncstat = nf90_inq_varid(ifid,"lon",varid)								! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
	
ncstat = nf90_get_var(ifid,varid,lon)									! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)								! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)			

ncstat = nf90_get_var(ifid,varid,lat)									! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"time",varid)								! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)			

ncstat = nf90_get_var(ifid,varid,time)									! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


!----------------------------------------------------
! In terminal, call the programs coordstring and parsecoords to determine boundaries of area of interest (translates lat/long values into indices of the lat/long arrays)

call getarg(2,coordstring)												! Reads second argument in the command line (coordinates in lat/long format divided by /)

call parsecoords(coordstring,bounds)

write(so,*)bounds														! Print boundaries of area of interest

call calcpixels(lon,lat,bounds,xpos,ypos,srtx,srty,cntx,cnty)

write(so,*)'xpos: ',xpos												! Print start position of longitude
write(so,*)'ypos: ',ypos												! Print start position of latitude
write(so,*)'cntx, cnty: ',cntx,cnty										! Print counts of x and y cells	

allocate(var_in(cntx,cnty,tlen))										! Allocate space in input array 'var_in' (x-range, y-range and temporal range)


!---------------------------------------------------------------------
! get the TEMPERATURE array timeseries

allocate(tmp(cntx,cnty,tlen))											! Allocate space of area of interest and temporal range to tmp array (mean monthly temperature)

tmp = -9999.															! Set tmp to -9999

ncstat = nf90_inq_varid(ifid,"tmp",varid)								! Get variable ID of variable tmp 
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

allocate(dtr(cntx,cnty,tlen))											! Allocate space to array dtr 

dtr = -9999.															! Set dtr to -9999								

ncstat = nf90_inq_varid(ifid,"dtr",varid)								! Get variable ID of variable dtr
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
! get the PRECIPITATION array timeseries

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
! get the NUMBER OF WET DAYS array timeseries

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
! get the CLOUD COVER array timeseries

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
! get the WIND SPEED array timeseries

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
! get the MINIMUM TEMPERATURE array timeseries

allocate(mtmin(cntx,cnty,tlen))															! Allocate space to array tmin (monthly minimum temperature)

mtmin = -9999.

mtmin = tmp - 0.5 * dtr


!---------------------------------------------------------------------
! get the MAXIMUM TEMPERATURE array timeseries

allocate(mtmax(cntx,cnty,tlen))															! Allocate space to array tmax (monthly minimum temperature)

mtmax = -9999.

mtmax = tmp + 0.5 * dtr


!---------------------------------------------------------------------
! get the FRACTION OF WET DAYS array timeseries

allocate(wetf(cntx,cnty,tlen))												

wetf = -9999.																			! Set wetf to -9999

allocate(nd(tlen))																		! Allocate space to array nd (length: time (tlen))

i = 1

endyr = startyr + nyrs - 1																! Calculate endyr (start year + number of years in time series)

! apply function ndaymonth to the time series:
do yr = startyr,endyr																	! Do, for every year from start year to end year
  do mon = 1,12																			! and for every month from 1 to 12
    
    nd(i) = ndaymonth(yr,mon)																

      
     write(*,*)i,yr,mon,nd(i)

    
    i = i + 1 
    
  end do
end do

! use output of step above to calculate fraction of wet days a month
do j = 1,cnty																			! Do, for every year j in amount of years
  do i = 1,cntx																			! and for every month 
    wetf(i,j,:) = wet(i,j,:) / real(nd)													! calculate fraction of wet days and save it in 'wetf"

    write(*,*)i,j,y,wet(i,j,t),nd(t)

  end do
end do


!---------------------------------------------------------------------
! DATA CHECK W/ OUTPUT

write(so,*)'Year, month, tmin (C), tmean (C), tmax (C), precip (mm), wet days, &
cloud cover fraction, wind speed (m/s), fraction wet days:'								! Caption for data check output

i = 1
do y = 1,nyrs																			! Do, for each year
  do m = 1,nmos																			! and within each year, for each month
  	
    write(so,'(2i5, 8f9.2)')y+startyr-1,m,&												! print the attribute values for (lat, long, time step)
    mtmin(1,1,i),tmp(1,1,i),mtmax(1,1,i),& 
    pre(1,1,i), wet(1,1,i), cld(1,1,i)/100, wnd(1,1,i), wetf(1,1,i)
    
    i = i + 1																			! then add one to the time step

  end do
end do


!---------------------------------------------------------------------
! PREPARE data for SMOOTHING

i = 1																		! Two placeholders i,j 
j = 1													


! SMOOTHING DO LOOP: 

do t = 1,tlen																! For very time step in tlen

  if (t == 1) then  														! If t is the first time step (first month), use the current month's values for the preceding month
  
    mtmin3m(-1) = mtmin(i,j,t)												! Preceding month gets current months's values
    mtmin3m(0)  = mtmin(i,j,t)	
    mtmin3m(1)  = mtmin(i,j,t+1)
    mtmax3m(-1) = mtmax(i,j,t)
    mtmax3m(0)  = mtmax(i,j,t)
    mtmax3m(1)  = mtmax(i,j,t+1)
    nd3m(-1) = nd(t)
    nd3m(0)  = nd(t)
    nd3m(1)  = nd(t+1)
    cld3m(-1) = cld(i,j,t)
    cld3m(0)  = cld(i,j,t)
    cld3m(1)  = cld(i,j,t+1)
    wnd3m(-1) = wnd(i,j,t)
    wnd3m(0)  = wnd(i,j,t)
    wnd3m(1)  = wnd(i,j,t+1)
  
  else if (t == tlen) then														! If it's the last time step (last month), use the current month's values for the month after the current one
  
    mtmin3m(-1) = mtmin(i,j,t-1)			
    mtmin3m(0)  = mtmin(i,j,t)
    mtmin3m(1)  = mtmin(i,j,t)												! Month after current one (not in the time frame anymore) gets current month's values
    mtmax3m(-1) = mtmax(i,j,t-1)
    mtmax3m(0)  = mtmax(i,j,t)
    mtmax3m(1)  = mtmax(i,j,t)
    nd3m(-1) = nd(t-1)
    nd3m(0)  = nd(t)
    nd3m(1)  = nd(t)
    cld3m(-1) = cld(i,j,t-1)
    cld3m(0)  = cld(i,j,t)
    cld3m(1)  = cld(i,j,t)
    wnd3m(-1) = wnd(i,j,t-1)
    wnd3m(0)  = wnd(i,j,t)
    wnd3m(1)  = wnd(i,j,t)
    
  else																		! All other cases: time steps include month before current month, current month, 
  																			! and month after the current month
    mtmin3m(-1) = mtmin(i,j,t-1)
    mtmin3m(0)  = mtmin(i,j,t)
    mtmin3m(1)  = mtmin(i,j,t+1)
    mtmax3m(-1) = mtmax(i,j,t-1)
    mtmax3m(0)  = mtmax(i,j,t)
    mtmax3m(1)  = mtmax(i,j,t+1)
	nd3m(-1) = nd(t-1)
    nd3m(0)  = nd(t)
    nd3m(1)  = nd(t+1)
	cld3m(-1) = cld(i,j,t-1)
    cld3m(0)  = cld(i,j,t)
    cld3m(1)  = cld(i,j,t+1)
    wnd3m(-1) = wnd(i,j,t-1)
    wnd3m(0)  = wnd(i,j,t)
    wnd3m(1)  = wnd(i,j,t+1)
    
  end if

  bcond_tmin(1) = mtmin3m(-1)												! Set boundary conditions for variables 
  bcond_tmin(2) = mtmin3m(1)
  bcond_tmax(1) = mtmax3m(-1)
  bcond_tmax(2) = mtmax3m(1)
  bcond_nd(1) = nd3m(-1)
  bcond_nd(2) = nd3m(1)
  bcond_cld(1) = cld3m(-1)
  bcond_cld(2) = cld3m(1)
  bcond_wnd(1) = wnd3m(-1)
  bcond_wnd(2) = wnd3m(1)

  call rmsmooth(mtmin3m,nd3m,bcond_tmin,tmin_sm)							! Smooth minimum variables
  call rmsmooth(mtmax3m,nd3m,bcond_tmax,tmax_sm)
  call rmsmooth(cld3m,nd3m,bcond_cld,cld_sm)
  call rmsmooth(wnd3m,nd3m,bcond_wnd,wnd_sm)

end do																		! End smoothing loop


!---------------------------------------------------------------------
! CLOSE INPUT FILE

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


!---------------------------------------------------------------------
!---------------------------------------------------------------------
! FUNCTION TO DETERMINE AMOUNT OF DAYS IN A MONTH

contains

integer function ndaymonth(yr,mon)													! Function to find out the number of days in a month, considering leap years and the year given as AD

! Input: Current year and month
integer, intent(in) :: yr 															
integer, intent(in) :: mon 															

! Arrays defining standard and leap years
integer, parameter, dimension(12) :: std_year = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
integer, parameter, dimension(12) :: leapyear = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]


! If-statement to find out which kind of year it is:

if (mod(yr,400) == 0) then					! If year can be divided by 400, then it's a leap year

	ndaymonth = leapyear(mon)				! Chose amount of days from leap year array		

  else if (mod(yr,100) == 0) then			! If year can be divided by 100, then it's a standard year 

	ndaymonth = std_year(mon)					

  else if (mod(yr,4) == 0) then				! If year can be divided by 4, it's a leapyear

	ndaymonth = leapyear(mon)

  else										! Else, it's a standard year 

	ndaymonth = std_year(mon)

end if 

end function ndaymonth



!---------------------------------------------------------------------
! END PROGRAM GWGEN_GRID

end program gwgen_grid