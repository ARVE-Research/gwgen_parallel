program gwgen_grid

! Use Makefile to compile this program

! Program to run gwgen with gridded input, provide the name of a climate data input file and geographic bounds for the simulation using a xmin/xmax/ymin/ymax string
! JO Kaplan, HKU, 2019

! Terminal command line: ./gwgen_grid ~/path/to/input.file x_coordinate/y_coordinate

! List of modules that will be used and the variables within these modules that will be used in this script:
use parametersmod, only : sp,dp,i4,i2,so
use errormod, only : ncstat,netcdf_err
use coordsmod, only : coordstring,bounds,parsecoords,calcpixels
use netcdf
use geohashmod, only : geohash
use randomdistmod, only : ran_seed
use weathergenmod, only : metvars_in, metvars_out, weathergen, init_weathergen

implicit none

! Inquire about the dimensions of the input file

character(100) :: infile               ! input file name

integer :: xlen                        ! length of dimension 'lat'
integer :: ylen                        ! length of dimension 'long'
integer :: tlen                        ! length of dimension 'time'

! IDs for file, dimensions and variables

integer :: ifid                        ! Input file ID
integer :: dimid                       ! Dimension ID
integer :: varid                       ! Variable ID

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

real(sp), allocatable, dimension(:,:,:) :: tmp        ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:,:) :: pre        ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:,:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:,:) :: cld        ! mean monthly cloud cover (percent)
real(sp), allocatable, dimension(:,:,:) :: wnd        ! mean monthly 10m windspeed (m s-1)

! Monthly input attributes calculated here

real(sp), allocatable, dimension(:,:,:) :: mtmin        ! maximum monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: mtmax        ! monthly minimum temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: wetf        ! fraction of wet days in a month

real(sp) :: scale_factor                  ! Value for the calculation of the "real" value of the parameters. Can be found in the netCDF file
real(sp) :: add_offset                    ! Value for the calculation of the "real" value of the parameters. Can be found in the netCDF file
integer(i2) :: missing_value                ! Missing values in the input file


! Elements to calculate current year and amount of days in current month

integer :: n,i_count,outd
integer :: i,j,t,d,y,m,s
integer :: nyrs                        ! Number of years (tlen/12)
integer :: d0
integer :: d1
integer :: calyr

integer, parameter :: nmos = 12                ! Number of months

integer, parameter :: startyr = 1871            ! Start year set to 1871
integer :: endyr                      ! End year (2010, not set yet)

integer :: yr                        ! Variable year 
integer :: mon                        ! Variable month 

integer, allocatable, dimension(:) :: nd  


! Elements for the smoothing process:

integer, parameter :: w = 3  !number of months before and after to use in the smoothing

real, dimension(-w:w) :: mtminbuf              ! Scalar for monthly minimum temperature
real, dimension(-w:w) :: mtmaxbuf             ! Scalar for monthly maximum temperature
integer, dimension(-w:w) :: ndbuf              ! Scalar for number of days in the month
real, dimension(-w:w) :: cldbuf                ! Scalar for monthly cloud fractions
real, dimension(-w:w) :: wndbuf                ! Scalar for monthly wind speeds  

real(sp), dimension(2) :: bcond_tmin                 ! boundary conditions of min temp for smoothing
real(sp), dimension(2) :: bcond_tmax                 ! boundary conditions of max temp for smoothing
real(sp), dimension(2) :: bcond_cld                ! boundary conditions of cloud for smoothing
real(sp), dimension(2) :: bcond_wnd                 ! boundary conditions of wind speed for smoothing        
real(sp), dimension(2) :: bcond_nd                   ! boundary conditions of number of days for smoothing

real(sp), dimension(31*(1+2*w)) :: tmin_sm              ! smoothed daily values of min temperature
real(sp), dimension(31*(1+2*w)) :: tmax_sm               ! smoothed daily values of max temperature
real(sp), dimension(31*(1+2*w)) :: cld_sm               ! smoothed daily values of cloudiness
real(sp), dimension(31*(1+2*w)) :: wnd_sm                ! smoothed daily values of wind speed


! For weathergen do loop 

integer  :: mwetd_sim                           ! monthly simulation of wet days  
integer  :: pdaydiff
integer  :: pdaydiff1                       

real(sp) :: mprec_sim                            ! monthly simulation of precipitation
real(sp) :: precdiff
real(sp) :: precdiff1
real(sp) :: prec_t
integer, parameter  :: pday_t = 3                          

type(metvars_in)  :: met_in
type(metvars_out) :: met_out

type(metvars_out), dimension(31) :: month_met



!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

call init_weathergen()

!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

call getarg(1,infile)                    ! Reads first argument in the command line (path to input file)
  
ncstat = nf90_open(infile,nf90_nowrite,ifid)        ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)          ! get dimension ID from dimension 'lon' in the input file (file id, dimension name, dimension ID)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)    ! get dimension name and length from input file for dimension previously inquired (file id, dimension ID, length will be written in variable xlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)          ! Get dimension ID from dimension 'lat' 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)    ! Get length of dimension 'lat' and assign it to variable ylen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)          ! Get dimension ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)    ! Get length of dimension 'time' and assign it to variable tlen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

nyrs = tlen / 12                      ! Number of years = months / 12

! write(so,*)'xlen, ylen, nyears: ',xlen,ylen,nyrs


!----------------------------------------------------
! Read variable IDs and values 

allocate(lon(xlen))                            ! Allocate length to longitude array
allocate(lat(ylen))                            ! Allocate length to latitude array
allocate(time(tlen))                          ! Allocate length to latitude array

ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_get_var(ifid,varid,lon)                  ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      

ncstat = nf90_get_var(ifid,varid,lat)                  ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"time",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      

ncstat = nf90_get_var(ifid,varid,time)                  ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


!----------------------------------------------------
! In terminal, call the programs coordstring and parsecoords to determine boundaries of area of interest (translates lat/long values into indices of the lat/long arrays)

call getarg(2,coordstring)                        ! Reads second argument in the command line (coordinates in lat/long format divided by /)

call parsecoords(coordstring,bounds)

! write(so,*)'Bounds: ',bounds                            ! Print boundaries of area of interest

call calcpixels(lon,lat,bounds,xpos,ypos,srtx,srty,cntx,cnty)

!write(so,*)'xpos: ',xpos                        ! Print start position of longitude
!write(so,*)'ypos: ',ypos                        ! Print start position of latitude
! write(so,*)'cntx, cnty: ',cntx,cnty                    ! Print counts of x and y cells  

allocate(var_in(cntx,cnty,tlen))                    ! Allocate space in input array 'var_in' (x-range, y-range and temporal range)

deallocate(lon)
deallocate(lat)

! reallocate coordinate variables to only hold the selected subset of the grid and get the data

allocate(lon(cntx))
allocate(lat(cnty))

ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_get_var(ifid,varid,lon,start=[srtx],count=[cntx])                  ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      

ncstat = nf90_get_var(ifid,varid,lat,start=[srty],count=[cnty])                  ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


!---------------------------------------------------------------------
! get the TEMPERATURE array timeseries

allocate(tmp(cntx,cnty,tlen))                      ! Allocate space of area of interest and temporal range to tmp array (mean monthly temperature)

tmp = -9999.                              ! Set tmp to -9999

ncstat = nf90_inq_varid(ifid,"tmp",varid)                ! Get variable ID of variable tmp 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    ! Get values for variable tmp from input file, starting at the starting point and going for cnt x and y cells
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)              ! Get attribute 'missing value' in the variable temperature
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)              ! Get attribute 'scale factor' 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)                ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) tmp = real(var_in) * scale_factor + add_offset      ! Where the temperature attribute is not missing value, calculate the real temperature using the scale factor and the add_offset


!---------------------------------------------------------------------
! get the diurnal temperature range array timeseries

allocate(dtr(cntx,cnty,tlen))                      ! Allocate space to array dtr 

dtr = -9999.                              ! Set dtr to -9999                

ncstat = nf90_inq_varid(ifid,"dtr",varid)                ! Get variable ID of variable dtr
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    ! Get values for variable dtr from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)              ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)              ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)                ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) dtr = real(var_in) * scale_factor + add_offset      ! Where dtr is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the PRECIPITATION array timeseries

allocate(pre(cntx,cnty,tlen))                              ! Allocate space to array pre (precipitation) 

pre = -9999.                                      ! Set pre to -9999                

ncstat = nf90_inq_varid(ifid,"pre",varid)                        ! Get variable ID of variable pre
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    ! Get values for variable pre from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)              ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)              ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)                ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) pre = real(var_in) * scale_factor + add_offset      ! Where pre is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the NUMBER OF WET DAYS array timeseries

allocate(wet(cntx,cnty,tlen))                              ! Allocate space to array wet (precipitation) 

wet = -9999.                                      ! Set wet to -9999                

ncstat = nf90_inq_varid(ifid,"wet",varid)                        ! Get variable ID of variable wet
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    ! Get values for variable wet from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)              ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)              ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)                ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) wet = real(var_in) * scale_factor + add_offset      ! Where wet is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the CLOUD COVER array timeseries

allocate(cld(cntx,cnty,tlen))                              ! Allocate space to array cld (precipitation) 

cld = -9999.                                      ! Set cld to -9999                

ncstat = nf90_inq_varid(ifid,"cld",varid)                        ! Get variable ID of variable cld
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    ! Get values for variable cld from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)              ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)              ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)                ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) cld = real(var_in) * scale_factor + add_offset      ! Where cld is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the WIND SPEED array timeseries

allocate(wnd(cntx,cnty,tlen))                              ! Allocate space to array wnd (precipitation) 

wnd = -9999.                                      ! Set wnd to -9999                

ncstat = nf90_inq_varid(ifid,"wnd",varid)                        ! Get variable ID of variable wnd
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    ! Get values for variable wnd from input file, based on area and time scale of interest
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)              ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)              ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)                ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (var_in /= missing_value) wnd = real(var_in) * scale_factor + add_offset      ! Where wnd is not missing value, calculate real values using add_offset and scale_factor


!---------------------------------------------------------------------
! get the MINIMUM TEMPERATURE array timeseries

allocate(mtmin(cntx,cnty,tlen))                              ! Allocate space to array tmin (monthly minimum temperature)

mtmin = -9999.

mtmin = tmp - 0.5 * dtr


!---------------------------------------------------------------------
! get the MAXIMUM TEMPERATURE array timeseries

allocate(mtmax(cntx,cnty,tlen))                              ! Allocate space to array tmax (monthly minimum temperature)

mtmax = -9999.

mtmax = tmp + 0.5 * dtr


!---------------------------------------------------------------------
! get the FRACTION OF WET DAYS array timeseries

allocate(wetf(cntx,cnty,tlen))                        

wetf = -9999.                                      ! Set wetf to -9999

allocate(nd(tlen))                                    ! Allocate space to array nd (length: time (tlen))

i = 1

endyr = startyr + nyrs - 1                                ! Calculate endyr (start year + number of years in time series)

! apply function ndaymonth to the time series:
do yr = startyr,endyr                                  ! Do, for every year from start year to end year
  do mon = 1,12                                      ! and for every month from 1 to 12
    
    nd(i) = ndaymonth(yr,mon)                                

      
!      write(*,*)i,yr,mon,nd(i)

    
    i = i + 1 
    
  end do
end do

t = 1

! use output of step above to calculate fraction of wet days a month
do j = 1,cnty                                      ! Do, for every year j in amount of years
  do i = 1,cntx                                      ! and for every month 
  
    ! enforce reasonable values of prec and wetdays

    where (pre(i,j,:) >= 0.1)
  
      wet(i,j,:) = max(wet(i,j,:),1.)      ! enforce at least one wet day if there is any mesurable precip
      wet(i,j,:) = min(wet(i,j,:),10. * pre(i,j,:))
      wetf(i,j,:) = wet(i,j,:) / real(nd)  ! calculate fraction of wet days and save it in 'wetf"

    elsewhere

      pre(i,j,:) = 0.   ! zero out precip if it is less than 0.1 mm per month
      wet(i,j,:) = 0.
      wetf(i,j,:) = 0.
 
    end where
  end do
end do


!---------------------------------------------------------------------
! DATA CHECK W/ OUTPUT

!write(so,*)'Year, month, tmin (C), tmean (C), tmax (C),precip (mm), wet days, cloud cover fraction, wind speed (m/s),fraction wet days:'                            

! i = 1
! do y = 1,2                                      ! Do, for each year
!   do m = 1,nmos                                      ! and within each year, for each month
!     
! !     write(so,'(2i5, 8f9.2)')y+startyr-1,m,&                        ! print the attribute values for (lat, long, time step)
! !     mtmin(1,1,i),tmp(1,1,i),mtmax(1,1,i),& 
! !     pre(1,1,i), wet(1,1,i), cld(1,1,i)/100, wnd(1,1,i), wetf(1,1,i)
!     
!     i = i + 1                                      ! then add one to the time step
! 
!   end do
! end do


!---------------------------------------------------------------------
! PREPARE data for SMOOTHING

i = 1                                     
j = 1                          

! initialize the random number generator for this gridcell so the stream of random numbers is always the same

call ran_seed(geohash(lon(i),lat(j)),met_in%rndst)



! FINAL OUTPUT WRITE STATEMENT - HEADER:

write(*,*)'Year, Month, Day, ',& 
              'mtmin, mtmax, mtmean, mcloud, mwind, mprecip,',&
              'sm_tmin, sm_tmax, sm_cloud_fr, sm_wind, ',&
              'daily_tmin, daily_tmax, daily_cloud_fr, daily_wind, daily_precip'

! SMOOTHING DO LOOP: 

! write(*,*)'start init loop'

! initialize the smoothing buffer variables with all values from the first month

mtminbuf = mtmin(i,j,1)  !this copies the first value across all of those in the buffer
mtmaxbuf = mtmax(i,j,1)  !this copies the first value across all of those in the buffer
ndbuf    = nd(1)  !this copies the first value across all of those in the buffer
cldbuf   = cld(i,j,1)  !this copies the first value across all of those in the buffer
wndbuf   = wnd(i,j,1)  !this copies the first value across all of those in the buffer

! look ahead the filter width worth of months and preinitialize the buffers

do s = 1,w

  mtminbuf = eoshift(mtminbuf,1,mtmin(i,j,t+s))
  mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(i,j,t+s))
  ndbuf    = eoshift(ndbuf,1,nd(t+s))
  cldbuf   = eoshift(cldbuf,1,cld(i,j,t+s))
  wndbuf   = eoshift(wndbuf,1,wnd(i,j,t+s))

end do

! start time loop

! write(*,*)'start time loop'



met_out%pday(1) = .false.
met_out%pday(2) = .false.
met_out%resid = 0.

do yr = 1,nyrs-1
  do m = 1,12

    t = m + 12 * (yr - 1)


    bcond_tmin = [mtminbuf(-w),mtminbuf(w)]    ! Set boundary conditions for variables 
    bcond_tmax = [mtmaxbuf(-w),mtmaxbuf(w)]    ! Set boundary conditions for variables 
    bcond_nd   = [ndbuf(-w),ndbuf(w)]                                      ! Set boundary conditions for variables 
    bcond_cld  = [cldbuf(-w),cldbuf(w)]                                      ! Set boundary conditions for variables 
    bcond_wnd  = [wndbuf(-w),wndbuf(w)]                                      ! Set boundary conditions for variables 
  
!   write(*,*)'ndmonth buffer',ndbuf
  
  call rmsmooth(mtminbuf,ndbuf,bcond_tmin,tmin_sm(1:sum(ndbuf)))                  ! Smooth minimum variables
  call rmsmooth(mtmaxbuf,ndbuf,bcond_tmax,tmax_sm(1:sum(ndbuf)))
  call rmsmooth(cldbuf,ndbuf,bcond_cld,cld_sm(1:sum(ndbuf)))
  call rmsmooth(wndbuf,ndbuf,bcond_wnd,wnd_sm(1:sum(ndbuf)))


! write(*,*)'Smoothed data: Day,  minimum temp (C), maximum temp (C), cloud fraction, wind speed (m/s)'

d0 = sum(ndbuf(-w:-1))+1
d1 = d0 + ndbuf(0) - 1 

!   s = 1
!   do d = d0,d1
!    
!      write(*,'(2i5,4f9.2)')s,d,tmin_sm(d),tmax_sm(d),cld_sm(d),wnd_sm(d)
!  
!      s = s + 1
!  
!   end do
! 
! read(*,*)

!---------------------------------------------------------------------
! Long loop calling WEATHERGEN

! save the current state of the residuals and pday


prec_t = max(0.5,0.05 * pre(i,j,t))                                   !set quality threshold for preciptation amount

i_count = 1




do    ! quality control loop
! do y = 1,1                                      ! Do, for each year
! do m = 1,3 

  mwetd_sim = 0
  mprec_sim = 0.

outd = 1
do d = d0,d1  ! day loop

!      write(*,'(2i5,4f9.2)')outd,d,tmin_sm(d),tmax_sm(d),cld_sm(d),wnd_sm(d)


! write(*,*)yr,m,d0,d,d1,ndbuf

  met_in%prec = pre(i,j,t)
  met_in%wetd = wet(i,j,t)
  met_in%wetf = wetf(i,j,t)                                     ! OR just use real(wetf(t))
  

  
      !start day loop

          met_in%tmin  = tmin_sm(d)
          met_in%tmax  = tmax_sm(d)
          met_in%cldf  = real(cld_sm(d)) * 0.01
          met_in%wind  = real(wnd_sm(d))
          met_in%pday  = met_out%pday
          met_in%resid = met_out%resid
          
          call weathergen(met_in,met_out)
          
          
!           	    write(so,'(2i5, 8f9.2)')y+startyr-1,m,&                        ! print the attribute values for (lat, long, time step)
! 			    mtmin(1,1,i),tmp(1,1,i),mtmax(1,1,i),& 
! 			    pre(1,1,i), wet(1,1,i), cld(1,1,i)/100, wnd(1,1,i), wetf(1,1,i)
          

          

          met_in%rndst = met_out%rndst
          month_met(outd) = met_out    ! save this day into a month holder

          if (met_out%prec > 0.) then
          
!               write(*,*)met_out%prec,mprec_sim,mwetd_sim
              
              mwetd_sim = mwetd_sim + 1
              mprec_sim = mprec_sim + met_out%prec
              
          end if
          
      !end of month
!       write(*,'(2i5,4f9.2)')outd,d,tmin_sm(d),tmax_sm(d),cld_sm(d),wnd_sm(d)

      outd = outd + 1

  end do  ! day loop            

  ! quality control checks                                               

      if (pre(i,j,t) == 0.) then
          pdaydiff = 0
          precdiff = 0.
          exit

      else if (i_count >= 2) then                                   !enforce at least two times over the month to get initial values ok
          pdaydiff = abs(wet(i,j,t) - mwetd_sim)
          precdiff = abs(pre(i,j,t) - mprec_sim)
          
!           write(*,*)'quality',i_count,pre(i,j,t),mprec_sim,wet(i,j,t),mwetd_sim,wetf(i,j,t)

          ! restrict simulated total monthly precip to +/-5% or 0.5 mm of observed value
          if (pdaydiff <= pday_t .and. precdiff <= prec_t) then
              exit

          else if (pdaydiff < pdaydiff1 .and. precdiff < precdiff1) then
              !save the values you have in a buffer in case you have to leave the loop
              pdaydiff1 = pdaydiff
              precdiff1 = precdiff

          else if (i_count > 1000) then
              write (*,*) "No good solution found after 1000 iterations."
              stop 1

          end if

      end if
      
!        write(0,*)'failed quality check',i_count,' re-running this month'

      i_count = i_count + 1




end do  ! end of quality control loop

! write(*,*)'month completed'



!---------------------------------------------------------------------
! DATA CHECK W/ OUTPUT

! write(*,*)'passed quality check, moving on to next month'

! read(*,*)

    mtminbuf = eoshift(mtminbuf,1,mtmin(i,j,t+w+1))
    mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(i,j,t+w+1))
    ndbuf    = eoshift(ndbuf,1,nd(t+w+1))
    cldbuf   = eoshift(cldbuf,1,cld(i,j,t+w+1))
    wndbuf   = eoshift(wndbuf,1,wnd(i,j,t+w+1))



           calyr = yr+startyr-1
           if (calyr > 1874 .and. calyr <= 1875) then
          
          do outd = 1,ndaymonth(calyr,m)
          
           ! FINAL OUTPUT WRITE STATEMENT
              write(*,'(3i5, 15f9.2)')calyr, m, outd,&
              mtmin(1,1,t), mtmax(1,1,t), tmp(1,1,t), cld(1,1,t)/100, wnd(1,1,t), pre(1,1,t), &
              met_in%tmin, met_in%tmax, met_in%cldf, met_in%wind,&
              month_met(outd)%tmin, month_met(outd)%tmax, month_met(outd)%cldf, month_met(outd)%wind, month_met(outd)%prec

           end do

           end if


end do                                    ! end month loop
end do  ! end year loop

! d-d0+1
!---------------------------------------------------------------------
! CLOSE INPUT FILE

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------


!#########################################
!##                                     ##
!##   F  U  N  C  T  I  O  N  S         ##
!##                                     ##
!#########################################


! FUNCTION TO DETERMINE AMOUNT OF DAYS IN A MONTH

contains

integer function ndaymonth(yr,mon)                          ! Function to find out the number of days in a month, considering leap years and the year given as AD

! Input: Current year and month
integer, intent(in) :: yr                               
integer, intent(in) :: mon                               

! Arrays defining standard and leap years
integer, parameter, dimension(12) :: std_year = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
integer, parameter, dimension(12) :: leapyear = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]


! If-statement to find out which kind of year it is:

if (mod(yr,400) == 0) then          ! If year can be divided by 400, then it's a leap year

  ndaymonth = leapyear(mon)        ! Chose amount of days from leap year array    

  else if (mod(yr,100) == 0) then      ! If year can be divided by 100, then it's a standard year 

  ndaymonth = std_year(mon)          

  else if (mod(yr,4) == 0) then        ! If year can be divided by 4, it's a leapyear

  ndaymonth = leapyear(mon)

  else                    ! Else, it's a standard year 

  ndaymonth = std_year(mon)

end if 

end function ndaymonth



!---------------------------------------------------------------------
!---------------------------------------------------------------------
! FUNCTION TO SMOOTH MONTHLY INPUT DATA

! Iterative, mean preserving method to smoothly interpolate mean data to pseudo-sub-timestep values
! Adapted from Rymes, M.D. and D.R. Myers, 2001. Solar Energy (71) 4, 225-231.

subroutine rmsmooth(mv_sts,dmonth,bcond,r)

use parametersmod, only : sp

implicit none


! Arguments

real(sp), dimension(:), intent(in)  :: mv_sts       ! Vector of mean values at super-time step (e.g., monthly), minimum three values
integer,  dimension(:), intent(in)  :: dmonth       ! Vector of number of intervals for the time step (e.g., days per month)
real(sp), dimension(2), intent(in)  :: bcond        ! Boundary conditions for the result vector (1=left side, 2=right side)
real(sp), dimension(:), intent(out) :: r            ! Result vector of values at chosen time step


! Parameters

real(sp), parameter :: ot = 1. / 3


! Local variables

integer :: n
integer :: ni
integer :: a
integer :: b
integer :: i
integer :: j
integer :: k
integer :: l
integer, dimension(size(r)) :: g

real(sp) :: ck
real(sp), dimension(2) :: bc

n  = size(mv_sts)
ni = size(r)

bc = bcond


! INITIALIZE RESULT VECTOR

i = 1
do a = 1,n

  j = i
  do b = 1,dmonth(a)
    
    r(i) = mv_sts(a)
    g(i) = j
    
    i = i + 1
  
  end do
end do


! Iteratively SMOOTH AND CORRECT THE RESULT to preserve the mean

do i = 1,ni                              ! Iteration loop
  do j = 2,ni-1
    
    r(j) = ot * (r(j-1) + r(j) + r(j+1))                ! Equation 1

  
  end do


  r(1)  = ot * (bc(1)   + r(1)  +  r(2))               ! Equations 2
  r(ni) = ot * (r(ni-1) + r(ni) + bc(2))

  j = 1
  
  
  do k = 1,n                                           ! Calculate one correction factor per super-timestep
  
    a = g(j)                                           ! Index of the first timestep value of the super-timestep
    b = g(j) + dmonth(k) - 1                           ! Index of the last timestep value of the super-timestep

    ck = sum(mv_sts(k) - r(a:b)) / ni                       ! Equation 4

    do l = 1,dmonth(k)                                 ! Apply the correction to all timestep values in the super-timestep
  
      r(j) = r(j) + ck
      j = j + 1
  
    end do

     bc(1) = r(ni)                        ! Correction for circular conditions when using climatology (do not use for transient simulations)
     bc(2) = r(1)

  end do
end do

end subroutine rmsmooth



!---------------------------------------------------------------------
!---------------------------------------------------------------------
! END PROGRAM GWGEN_GRID

end program gwgen_grid