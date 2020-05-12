program gwgen_grid

! Use the Makefile to compile this program

! Program to run gwgen with gridded input, provide the name of a climate data input file and
! geographic bounds for the simulation using a xmin/xmax/ymin/ymax string
! JO Kaplan, HKU, 2019

! Terminal command line: ./gwgen_grid ~/path/to/input.file x_coordinate/y_coordinate

! List of modules that will be used and the variables within these modules that will be used in this script:
use parametersmod, only : sp,dp,i4,i2,so,ndaymonth
use errormod,      only : ncstat,netcdf_err
use coordsmod,     only : coordstring,bounds,parsecoords,calcpixels
use geohashmod,    only : geohash
use randomdistmod, only : ran_seed
use weathergenmod, only : metvars_in, metvars_out, weathergen,rmsmooth
use netcdf

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

real(sp), allocatable, dimension(:,:,:) :: mtmin      ! maximum monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: mtmax      ! monthly minimum temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: wetf      ! fraction of wet days in a month

real(sp) :: scale_factor      ! Value for the calculation of the "real" value of the parameters. Can be found in the netCDF file
real(sp) :: add_offset        ! Value for the calculation of the "real" value of the parameters. Can be found in the netCDF file
integer(i2) :: missing_value  ! Missing values in the input file


! Elements to calculate current year and amount of days in current month

integer :: i_count,outd
integer :: i,j,t,d,m,s
integer :: nyrs                        ! Number of years (tlen/12)
integer :: d0
integer :: d1
integer :: calyr

integer, parameter :: nmos = 12   ! Number of months

integer, parameter :: startyr = 1871  ! Start year set to 1871
integer            :: endyr           ! End year (2010, not set yet)

integer :: yr    ! Variable year 
integer :: mon   ! Variable month 

integer, allocatable, dimension(:) :: nd  


! Variables for the smoothing process

integer, parameter :: w = 3              ! filter half-width for smoothing of monthly mean climate variables to pseudo-daily values (months)
integer, parameter :: wbuf = 31*(1+2*w)  ! length of the buffer in which to hold the smoothed pseudo-daily  meteorological values (days)

integer,  dimension(-w:w) :: ndbuf       ! number of days in the month
real(sp), dimension(-w:w) :: mtminbuf    ! monthly minimum temperature
real(sp), dimension(-w:w) :: mtmaxbuf    ! monthly maximum temperature
real(sp), dimension(-w:w) :: cldbuf      ! monthly cloud fractions
real(sp), dimension(-w:w) :: wndbuf      ! monthly wind speed

real(sp), dimension(2) :: bcond_tmin     ! boundary conditions of min temp for smoothing
real(sp), dimension(2) :: bcond_tmax     ! boundary conditions of max temp for smoothing
real(sp), dimension(2) :: bcond_cld      ! boundary conditions of cloud for smoothing
real(sp), dimension(2) :: bcond_wnd      ! boundary conditions of wind speed for smoothing        
real(sp), dimension(2) :: bcond_nd       ! boundary conditions of number of days for smoothing

real(sp), dimension(wbuf) :: tmin_sm     ! smoothed pseudo-daily values of min temperature
real(sp), dimension(wbuf) :: tmax_sm     ! smoothed pseudo-daily values of max temperature
real(sp), dimension(wbuf) :: cld_sm      ! smoothed pseudo-daily values of cloudiness
real(sp), dimension(wbuf) :: wnd_sm      ! smoothed pseudo-daily values of wind speed

! quality control variables

integer  :: mwetd_sim    ! simulated number of wet days  
real(sp) :: mprec_sim    ! simulated total monthly precipitation (mm)

integer  :: pdaydiff     ! difference between input and simulated wet days
real(sp) :: precdiff     ! difference between input and simulated total monthly precipitation (mm)

real(sp) :: prec_t       ! tolerance for difference between input and simulated total monthly precipitation (mm)
integer, parameter  :: wetd_t = 3  ! tolerance for difference between input and simulated wetdays (days)

integer  :: pdaydiff1 = huge(i4)   ! stored value of the best match difference between input and simulated wet days
real(sp) :: precdiff1 = huge(sp)   ! stored value of the difference between input and simulated total monthly precipitation (mm)

! data structures for meteorology

type(metvars_in)  :: met_in   ! structure containing one day of meteorology input to weathergen
type(metvars_out) :: met_out  ! structure containing one day of meteorology output from weathergen

type(metvars_out), dimension(31) :: month_met  ! buffer containing one month of simulated daily meteorology

!--------------------------------------------------------------------------------------------------------------------------------------------------
! program starts here

!-----------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

call getarg(1,infile)                                  ! Reads first argument in the command line (path to input file)
  
ncstat = nf90_open(infile,nf90_nowrite,ifid)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)              ! get dimension ID from dimension 'lon' in the input file 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! (file id, dimension name, dimension ID)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! (file id, dimension ID, length will be written in variable xlen)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)              ! Get dimension ID from dimension 'lat' 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)   ! Get length of dimension 'lat' and assign it to variable ylen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)             ! Get dimension ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)   ! Get length of dimension 'time' and assign it to variable tlen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

nyrs = tlen / 12                                       ! Calculate the number of years of data in the input file (months / 12)

!----------------------------------------------------
! Read variable IDs and values 

allocate(lon(xlen))       ! Allocate length to longitude array
allocate(lat(ylen))       ! Allocate length to latitude array
allocate(time(tlen))      ! Allocate length to latitude array

ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_get_var(ifid,varid,lon)                    ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      

ncstat = nf90_get_var(ifid,varid,lat)                    ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"time",varid)               ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      

ncstat = nf90_get_var(ifid,varid,time)                   ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!----------------------------------------------------
! Read the coordinates to run from the command line,
! call the programs coordstring and parsecoords 
! to determine boundaries of area of interest 
! (translates lat/long values into indices of the lat/long arrays)

! Read the second argument in the command line (coordinates in lat/long format divided by /)
call getarg(2,coordstring)                        

call parsecoords(coordstring,bounds)

call calcpixels(lon,lat,bounds,xpos,ypos,srtx,srty,cntx,cnty)

! Allocate space in input array 'var_in' (x-range, y-range and temporal range)
allocate(var_in(cntx,cnty,tlen))

! reallocate coordinate variables to only hold the selected subset of the grid and get the data
deallocate(lon)
deallocate(lat)

allocate(lon(cntx))
allocate(lat(cnty))

ncstat = nf90_inq_varid(ifid,"lon",varid)                        ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_get_var(ifid,varid,lon,start=[srtx],count=[cntx])  ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                        ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      

ncstat = nf90_get_var(ifid,varid,lat,start=[srty],count=[cnty])  ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------
! read the timeseries of monthly temperature

! Allocate space of area of interest and temporal range to tmp array (mean monthly temperature)
allocate(tmp(cntx,cnty,tlen)) 

ncstat = nf90_inq_varid(ifid,"tmp",varid)                         ! Get variable ID of variable tmp 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Get values for variable tmp from input file, starting at the starting point and going for cnt x and y cells
ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)   ! Get attribute 'missing value' in the variable temperature
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)     ! Get attribute 'scale factor' 
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)         ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Where the temperature attribute is not missing value, calculate the real temperature using the scale factor and the add_offset
where (var_in /= missing_value) tmp = real(var_in) * scale_factor + add_offset      

!---------------------------------------------------------------------
! read the timeseries of monthly diurnal temperature range

allocate(dtr(cntx,cnty,tlen))                                     ! Allocate space to array dtr

ncstat = nf90_inq_varid(ifid,"dtr",varid)                         ! Get variable ID of variable dtr
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Get values for variable dtr from input file, based on area and time scale of interest
ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)   ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)     ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)         ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Where dtr is not missing value, calculate real values using add_offset and scale_factor
where (var_in /= missing_value) dtr = real(var_in) * scale_factor + add_offset

!---------------------------------------------------------------------
! read the timeseries of monthly total precipitation

allocate(pre(cntx,cnty,tlen))                                     ! Allocate space to array pre (precipitation)

ncstat = nf90_inq_varid(ifid,"pre",varid)                         ! Get variable ID of variable pre
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Get values for variable pre from input file, based on area and time scale of interest
ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)   ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)     ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)         ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Where pre is not missing value, calculate real values using add_offset and scale_factor
where (var_in /= missing_value) pre = real(var_in) * scale_factor + add_offset      

!---------------------------------------------------------------------
! read the timeseries of number of days with precipitation > 0.1 mm (wet days)

allocate(wet(cntx,cnty,tlen))                                     ! Allocate space to array wet (precipitation)

ncstat = nf90_inq_varid(ifid,"wet",varid)                         ! Get variable ID of variable wet
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Get values for variable wet from input file, based on area and time scale of interest
ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)   ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)     ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)         ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Where wet is not missing value, calculate real values using add_offset and scale_factor
where (var_in /= missing_value) wet = real(var_in) * scale_factor + add_offset

!---------------------------------------------------------------------
! read the timeseries of cloud cover percent

allocate(cld(cntx,cnty,tlen))                                     ! Allocate space to array cld (precipitation) 

ncstat = nf90_inq_varid(ifid,"cld",varid)                         ! Get variable ID of variable cld
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Get values for variable cld from input file, based on area and time scale of interest
ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)   ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)     ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)         ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Where cld is not missing value, calculate real values using add_offset and scale_factor
where (var_in /= missing_value) cld = real(var_in) * scale_factor + add_offset      

!---------------------------------------------------------------------
! read the timeseries of windspeed

allocate(wnd(cntx,cnty,tlen))                                     ! Allocate space to array wnd (precipitation) 

ncstat = nf90_inq_varid(ifid,"wnd",varid)                         ! Get variable ID of variable wnd
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Get values for variable wnd from input file, based on area and time scale of interest
ncstat = nf90_get_var(ifid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])    
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)   ! Get attribute 'missing_value'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)     ! Get attribute 'scale_factor'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)         ! Get attribute 'add_offset'
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! Where wnd is not missing value, calculate real values using add_offset and scale_factor
where (var_in /= missing_value) wnd = real(var_in) * scale_factor + add_offset

!---------------------------------------------------------------------
! allocate space and calculate the derivative variables monthly tmin and tmax

allocate(mtmin(cntx,cnty,tlen))
allocate(mtmax(cntx,cnty,tlen))

mtmin = tmp - 0.5 * dtr
mtmax = tmp + 0.5 * dtr

!---------------------------------------------------------------------
! allocate and calculate the derivative variable wetf
! and run some quality control on precip and wet days

allocate(nd(tlen))                                   ! vector containing the number of in every month of the input timeseries

i = 1

endyr = startyr + nyrs - 1                           ! Calculate endyr (start year + number of years in time series)

! apply function ndaymonth to the time series:
do yr = startyr,endyr                                ! Do, for every year from start year to end year
  do mon = 1,12                                      ! and for every month from 1 to 12
    
    nd(i) = ndaymonth(yr,mon)                                

    i = i + 1 
    
  end do
end do

allocate(wetf(cntx,cnty,tlen))                        

t = 1

! use output of step above to calculate fraction of wet days a month
do j = 1,cnty                                        ! Do, for every year j in amount of years
  do i = 1,cntx                                      ! and for every month 
  
    ! enforce reasonable values of prec and wetdays

    where (pre(i,j,:) >= 0.1)
  
      wet(i,j,:)  = max(wet(i,j,:),1.)                ! enforce at least one wet day if there is any mesurable precip
      wet(i,j,:)  = min(wet(i,j,:),10. * pre(i,j,:))  ! do not allow wetdays to imply a mean daily precip smaller than the detection threshold  
      wetf(i,j,:) = wet(i,j,:) / real(nd)             ! calculate fraction of wet days and save it in 'wetf"

    elsewhere
    
      ! monthly total is less than the detection threshold (0.1 mm per month; happens with gridded datasets)
      ! zero out the pre and wet variables

      pre(i,j,:)  = 0.   
      wet(i,j,:)  = 0.
      wetf(i,j,:) = 0.
 
    end where
  end do
end do

!---------------------------------------------------------------------
! grid loop starts here

do j = 1,cnty
  do i = 1,cntx

    ! initialize the random number generator for this gridcell so the stream of random numbers is always the same

    call ran_seed(geohash(lon(i),lat(j)),met_in%rndst)

    ! FINAL OUTPUT WRITE STATEMENT - HEADER:

    ! write(*,'(a)')'X,Y,Year,Month,Day,mtmin,mtmax,mtmean,mcloud,mwind,mprecip,sm_tmin,sm_tmax,sm_cloud_fr,sm_wind,daily_tmin,daily_tmax,daily_cloud_fr,daily_wind,daily_precip'

    !---------------------------------------------------------------------
    ! prepare pseudo-daily smoothed meteorological variables
    ! initialize the smoothing buffer variables with all values from the first month

    ndbuf    = nd(1)         !this copies the first value across all of those in the buffer
    mtminbuf = mtmin(i,j,1)  !this copies the first value across all of those in the buffer
    mtmaxbuf = mtmax(i,j,1)  !this copies the first value across all of those in the buffer
    cldbuf   = cld(i,j,1)    !this copies the first value across all of those in the buffer
    wndbuf   = wnd(i,j,1)    !this copies the first value across all of those in the buffer

    ! look ahead the filter half-width worth of months and preinitialize the buffers

    do s = 1,w

      mtminbuf = eoshift(mtminbuf,1,mtmin(i,j,t+s))
      mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(i,j,t+s))
      ndbuf    = eoshift(ndbuf,1,nd(t+s))
      cldbuf   = eoshift(cldbuf,1,cld(i,j,t+s))
      wndbuf   = eoshift(wndbuf,1,wnd(i,j,t+s))

    end do
 
    met_out%pday(1) = .false.
    met_out%pday(2) = .false.
    met_out%resid = 0.

    ! start time loop

    do yr = 1,nyrs-1
      do m = 1,12

        t = m + 12 * (yr - 1)

        bcond_nd   = [ndbuf(-w),ndbuf(w)]        ! Set boundary conditions for variables
        bcond_tmin = [mtminbuf(-w),mtminbuf(w)]  ! Set boundary conditions for variables
        bcond_tmax = [mtmaxbuf(-w),mtmaxbuf(w)]  ! Set boundary conditions for variables
        bcond_cld  = [cldbuf(-w),cldbuf(w)]      ! Set boundary conditions for variables
        bcond_wnd  = [wndbuf(-w),wndbuf(w)]      ! Set boundary conditions for variables
  
        ! generate pseudo-daily smoothed meteorological variables (using means-preserving algorithm)
  
        call rmsmooth(mtminbuf,ndbuf,bcond_tmin,tmin_sm(1:sum(ndbuf)))  
        call rmsmooth(mtmaxbuf,ndbuf,bcond_tmax,tmax_sm(1:sum(ndbuf)))
        call rmsmooth(cldbuf,ndbuf,bcond_cld,cld_sm(1:sum(ndbuf)))
        call rmsmooth(wndbuf,ndbuf,bcond_wnd,wnd_sm(1:sum(ndbuf)))

        ! calculcate start and end positons of the current month pseudo-daily buffer

        d0 = sum(ndbuf(-w:-1))+1
        d1 = d0 + ndbuf(0) - 1
        
        ! restrict simulated total monthly precip to +/-5% or 0.5 mm of observed value

        prec_t = max(0.5,0.05 * pre(i,j,t))

        i_count = 0

        !---------------------------------------------------------------------------------
        ! quality control loop calling the weathergen -  this loop principally checks that
        ! the number of wet days and total precip stayed close to the input data

        do
          i_count = i_count + 1    ! increment iteration number

          mwetd_sim = 0
          mprec_sim = 0.

          outd = 1

          do d = d0,d1  ! day loop

            !write(*,*)yr,m,d0, d,d1,ndbuf

            met_in%prec  = pre(i,j,t)
            met_in%wetd  = wet(i,j,t)
            met_in%wetf  = wetf(i,j,t)
            met_in%tmin  = tmin_sm(d)
            met_in%tmax  = tmax_sm(d)
            met_in%cldf  = real(cld_sm(d)) * 0.01
            met_in%wind  = real(wnd_sm(d))
            met_in%pday  = met_out%pday
            met_in%resid = met_out%resid

            call weathergen(met_in,met_out)

            met_in%rndst = met_out%rndst
            month_met(outd) = met_out    ! save this day into a month holder

            if (met_out%prec > 0.) then

              mwetd_sim = mwetd_sim + 1
              mprec_sim = mprec_sim + met_out%prec

            end if

            outd = outd + 1

          end do  ! day loop            

          ! quality control checks                                               

          if (pre(i,j,t) == 0.) then ! if there is no precip in this month a single iteration is ok

            pdaydiff = 0
            precdiff = 0.

            exit

          else if (i_count < 2) then
          
            cycle  !enforce at least two times over the month to get initial values ok

          end if

          pdaydiff = abs(wet(i,j,t) - mwetd_sim)
          precdiff = abs(pre(i,j,t) - mprec_sim)

          if (pdaydiff <= wetd_t .and. precdiff <= prec_t) then

            exit

          else if (pdaydiff < pdaydiff1 .and. precdiff < precdiff1) then

            ! save the values you have in a buffer in case you have to leave the loop
            ! should save the entire monthly state so that the "closest" acceptable value
            ! could be used in the event of needing a very large number of iteration cycles

            pdaydiff1 = pdaydiff
            precdiff1 = precdiff

          else if (i_count > 1000) then

            write (*,*) "No good solution found after 1000 iterations."
            stop

          end if
 
        end do  

        ! end of quality control loop
        !---------------------------------------------------------------------------------

        mtminbuf = eoshift(mtminbuf,1,mtmin(i,j,t+w+1))
        mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(i,j,t+w+1))
        ndbuf    = eoshift(ndbuf,1,nd(t+w+1))
        cldbuf   = eoshift(cldbuf,1,cld(i,j,t+w+1))
        wndbuf   = eoshift(wndbuf,1,wnd(i,j,t+w+1))

        !-----------------------------------------------------------
        ! diagnostic output
    
        calyr = yr+startyr-1

        if (calyr > 1874 .and. calyr <= 1875) then
          do outd = 1,ndaymonth(calyr,m)

            ! FINAL OUTPUT WRITE STATEMENT
            write(*,'(5i5, 15f9.2)')i,j,calyr, m, outd,&
            mtmin(1,1,t), mtmax(1,1,t), tmp(1,1,t), cld(1,1,t)/100, wnd(1,1,t), pre(1,1,t), &
            tmin_sm(d0+outd-1), tmax_sm(d0+outd-1), (cld_sm(d0+outd-1)*0.01), wnd_sm(d0+outd-1), &                               ! met_in%cldf, met_in%wind,&
            month_met(outd)%tmin, month_met(outd)%tmax, month_met(outd)%cldf, month_met(outd)%wind, month_met(outd)%prec

          end do
        end if

      end do  ! month loop
    end do    ! year loop
  end do      ! columns
end do        ! rows

!---------------------------------------------------------------------
! CLOSE INPUT FILE

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end program gwgen_grid
