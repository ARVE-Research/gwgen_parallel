program gwgen_grid

! Use the Makefile to compile this program

!! Program to run gwgen with gridded input, provide the name of a climate data input file and
!! geographic bounds for the simulation using a xmin/xmax/ymin/ymax string
! JO Kaplan, HKU, 2019

! Terminal command line: ./gwgen_grid ~/path/to/input.file x_coordinate/y_coordinate

! List of modules that will be used and the variables within these modules that are used in this program:

use parametersmod, only : sp,dp,i4,i2,so,ndaymonth
use errormod,      only : ncstat,netcdf_err
use coordsmod,     only : index,parsecoords,calcpixels
use geohashmod,    only : geohash
use randomdistmod,
use newsplinemod,  only : newspline_all
use weathergenmod, only : metvars_in, metvars_out, weathergen,rmsmooth,roundto
use getdatamod,    only : readdata
use outputmod,     only : genoutfile,putlonlat
use netcdf
use mpi

implicit none

! Inquire about the dimensions of the input file

character(100) :: infile               ! input file name
character(100) :: outfile

integer :: xlen                        ! length of dimension 'lat'
integer :: ylen                        ! length of dimension 'long'
integer :: tlen                        ! length of dimension 'time'

! IDs for file, dimensions and variables

integer :: ifid                        ! Input file ID
integer :: dimid                       ! Dimension ID
integer :: varid                       ! Variable ID
integer :: ofid                        ! output file ID

! Allocatable arrays for longitude and latitude

real(dp), allocatable, dimension(:) :: lon
real(dp), allocatable, dimension(:) :: lat
real(dp), allocatable, dimension(:) :: time

! Start values of x and y (LLC), and counts of cells in both directions:

type(index), target :: id
type(index) :: timevals

integer, pointer :: srtx
integer, pointer :: srty
integer, pointer :: cntx
integer, pointer :: cnty

logical, allocatable, dimension(:,:) :: valid_pixel

! monthly input driver variables

real(sp), allocatable, dimension(:,:,:) :: tmp        ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:,:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:,:) :: pre        ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:,:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:,:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), allocatable, dimension(:,:,:) :: wnd        ! mean monthly 10m windspeed (m s-1)

! monthly derived driver variables

real(sp), allocatable, dimension(:) :: mtmin      ! maximum monthly temperature (degC)
real(sp), allocatable, dimension(:) :: mtmax      ! monthly minimum temperature (degC)
real(sp), allocatable, dimension(:) :: wetf       ! fraction of wet days in a month

! output variable

real(sp), allocatable, dimension(:,:) :: abs_tmin      ! absolute minimum temperature (degC)
real(sp), allocatable, dimension(:,:) :: abs_tmax      ! absolute maximum temperature (degC)

integer(i2), allocatable, dimension(:,:) :: outvar    ! Output variable for ncfile output adjusted by scale factor

real(sp) :: tmin_sim
real(sp) :: tmax_sim

! Elements to calculate current year and amount of days in current month

integer :: i_count,outd
integer :: i,j,t,d,m,s
integer :: nyrs                        ! Number of years (tlen/12)
integer :: d0
integer :: d1
integer :: calyr
integer :: ndm

integer :: yr    ! Variable year
integer :: mon   ! Variable month

! integer :: b

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
integer, parameter  :: wetd_t = 1  ! tolerance for difference between input and simulated wetdays (days)

integer  :: pdaydiff1 = huge(i4)   ! stored value of the best match difference between input and simulated wet days
real(sp) :: precdiff1 = huge(sp)   ! stored value of the difference between input and simulated total monthly precipitation (mm)

! data structures for meteorology

type(metvars_in)  :: met_in   ! structure containing one day of meteorology input to weathergen
type(metvars_out) :: met_out  ! structure containing one day of meteorology output from weathergen

type(metvars_out), dimension(31) :: month_met  ! buffer containing one month of simulated daily meteorology

real(sp) :: mtmin_sim
real(sp) :: mtmax_sim
real(sp) :: mcldf_sim
real(sp) :: mwind_sim

real(sp) :: prec_corr
real(sp) :: tmin_corr
real(sp) :: tmax_corr
real(sp) :: cldf_corr
real(sp) :: wind_corr

character(60) :: basedate

character(45) :: coordstring
character(60) :: timestring

integer :: t0
integer :: t1
integer :: p0
integer :: p1
integer :: cntt

integer :: baseyr
integer :: startyr
integer :: calcyrs

integer :: endyr

integer, dimension(3) :: srt
integer, dimension(3) :: cnt

integer :: nt

! Variables for output write option (Leo)
character(1) :: write
integer :: write_opt
integer :: l

!--- Debugging variables (Leo)
integer :: start_debug, end_debug
integer :: start_debug2, end_debug2
integer :: baddata_check
integer :: k


!-----------------------------------------------------------------------------------------------------------------------------------------
! program starts here

srtx => id%startx
srty => id%starty
cntx => id%countx
cnty => id%county

!-----------------------------------------------------
! Read output writing option (Leo)
! Input 1 as fifth argument on terminal to write out data (Leo)
call getarg(5, write)
write_opt = ichar(write) - 48

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

ncstat = nf90_inq_varid(ifid,"time",varid)               ! Get variable ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,time)                   ! Get variable values for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,'units',basedate)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!----------------------------------------------------
! Read the coordinates to run from the command line,
! call the programs coordstring and parsecoords
! to determine boundaries of area of interest
! (translates lat/long values into indices of the lat/long arrays)

! Read the second argument in the command line (coordinates in lat/long format divided by /)
call getarg(2,coordstring)

call parsecoords(coordstring,id)

call calcpixels(lon,lat,id)

deallocate(lon)
deallocate(lat)

! reallocate coordinate variables to only hold the selected subset of the grid and get the data

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
! to limit memory usage, only get as much data in the time dimension as we need for the run
! because the pseudo-daily smoother algorithm requires a window of 3 months on either side of the current month,
! we need to read in extra data before the first year of run after the last
! to ensure reliable results over the period of interest, always get one additional year of data at either end (if available)
! otherwise, copy the first and last years' data into the padding at either end of the data

baseyr  = 1871

call getarg(3,timestring)

call parsecoords(timestring,timevals)

startyr = timevals%minlon
calcyrs = timevals%minlat

cntt = 12 * (calcyrs + 2)          ! months, includes one-year buffer on either end

allocate(tmp(cntx,cnty,cntt))
allocate(dtr(cntx,cnty,cntt))
allocate(pre(cntx,cnty,cntt))
allocate(wet(cntx,cnty,cntt))
allocate(cld(cntx,cnty,cntt))
allocate(wnd(cntx,cnty,cntt))

! calculate file start and count indices

t0 = 1 + 12 * (startyr - baseyr)     ! index of the first month
t1 = t0 + 12 * calcyrs - 1           ! index of the last month

if (t0 == 1) then         ! there is no additional data to be had at the front end so copy the first year twice
  p0 = 13
else
  p0 = 1                  ! there is additional data before first year so grab it
  t0 = t0 - 12
end if

if (t1 == tlen) then      ! there is no additional data to be had at the back end so copy the last year twice
  p1 = cntt - 12
else
  p1 = cntt
  t1 = t1 + 12
end if

nt = t1 - t0 + 1

srt = [srtx,srty,t0]
cnt = [cntx,cnty,nt]

write(0,*)startyr,calcyrs
write(0,*)cntt,nt
write(0,*)t0,t1
write(0,*)p0,p1

!---------------------------------------------------------------------
! read the timeseries of monthly climate

call readdata(ifid,'tmp',srt,cnt,tmp(:,:,p0:p1))
call readdata(ifid,'dtr',srt,cnt,dtr(:,:,p0:p1))
call readdata(ifid,'pre',srt,cnt,pre(:,:,p0:p1))
call readdata(ifid,'wet',srt,cnt,wet(:,:,p0:p1))
call readdata(ifid,'cld',srt,cnt,cld(:,:,p0:p1))
call readdata(ifid,'wnd',srt,cnt,wnd(:,:,p0:p1))

!---------------------------------------------------------------------
! copy first and last year into buffers at each end

if (p0 == 13) then
  tmp(:,:,1:12) = tmp(:,:,13:24)
  dtr(:,:,1:12) = dtr(:,:,13:24)
  pre(:,:,1:12) = pre(:,:,13:24)
  wet(:,:,1:12) = wet(:,:,13:24)
  cld(:,:,1:12) = cld(:,:,13:24)
  wnd(:,:,1:12) = wnd(:,:,13:24)
end if

if (cntt > p1) then
  tmp(:,:,p1+1:cntt) = tmp(:,:,p1-11:p1)
  dtr(:,:,p1+1:cntt) = dtr(:,:,p1-11:p1)
  pre(:,:,p1+1:cntt) = pre(:,:,p1-11:p1)
  wet(:,:,p1+1:cntt) = wet(:,:,p1-11:p1)
  cld(:,:,p1+1:cntt) = cld(:,:,p1-11:p1)
  wnd(:,:,p1+1:cntt) = wnd(:,:,p1-11:p1)
end if

!---------------------------------------------------------------------
! enforce reasonable values of prec and wetdays

where (pre > 0.1)
  wet = max(wet,1.)
  wet = min(wet,10. * pre)
elsewhere
  pre = 0.
  wet = 0.
end where

!---------------------------------------------------------------------
! find valid pixels

allocate(valid_pixel(cntx,cnty))

where (tmp(:,:,1) /= -9999.)
  valid_pixel = .true.
elsewhere
  valid_pixel = .false.
end where

!---------------------------------------------------------------------
! calculate days per month

endyr = startyr + calcyrs - 1

allocate(nd(cntt))

i = 1
do yr = startyr-1,endyr+1
  do mon = 1,12

    nd(i) = ndaymonth(yr,mon)

    i = i + 1

  end do
end do

! write(0,*)cntx,cnty,cntt ! Debugging dianostics (Leo)

!---------------------------------------------------------------------
! create a netCDF output file with the dimensions of the area of interest

allocate(abs_tmin(cntx,cnty))
allocate(abs_tmax(cntx,cnty))

abs_tmin = 9999.
abs_tmax = -9999.

!---------------------------------------------------------------------
! grid loop starts here

allocate(mtmin(cntt))
allocate(mtmax(cntt))
allocate(wetf(cntt))

do j = 1,cnty

  write(0,*)'working on row',j

  do i = 1,cntx

    ! check that this gridcell has valid meteorology

    if (.not. valid_pixel(i,j)) cycle

    ! initialize the random number generator for this gridcell so the stream of random numbers is always the same

    call ran_seed(geohash(lon(i),lat(j)),met_in%rndst)

    !---------------------------------------------------------------------
    ! calculate derived climate variables

    mtmin = tmp(i,j,:) - 0.5 * dtr(i,j,:)
    mtmax = tmp(i,j,:) + 0.5 * dtr(i,j,:)
    wetf  = wet(i,j,:) / nd

    !--- Checking bad data (Leo)

    do k = 1,cntt

      if(mtmin(k) < -273.15) then

        write(0,*) 'Messed up min temp. at time slice', k, 'with', mtmin(k), 'degC'
        write(0,*) 'tmp: ', tmp(i,j,k) !, 'dtr: ', dtr(i,j,k)
        write(0,*) 'cntx:', i, 'cnty:', j

      end if

    end do

    !--- Cycling bad data (Leo)

    baddata_check = 0

    do k = 1,cntt

      if(mtmin(k) < -273.15) then

        baddata_check = baddata_check + 1

      end if

    end do

    if(baddata_check /= 0) cycle

    !--- Ressign abs_min vector for minval function after monthloop (Leo)
    ! abs_tmin(i,j) = 9999.


    !---------------------------------------------------------------------
    ! prepare pseudo-daily smoothed meteorological variables
    ! initialize the smoothing buffer variables with all values from the first month of the input
    ! since we will always start in year 2, we can do:

    t = 13

    ndbuf    = nd(t-w:t+w)
    mtminbuf = mtmin(t-w:t+w)
    mtmaxbuf = mtmax(t-w:t+w)
    cldbuf   = cld(i,j,t-w:t+w)
    wndbuf   = wnd(i,j,t-w:t+w)

    met_out%pday(1) = .false.
    met_out%pday(2) = .false.
    met_out%resid = 0.

    ! start year loop

    yearloop : do yr = 2,(calcyrs+1) ! Changed to calcyrs+1 (Leo)

      !write(0,*)yr,t0+yr-1

      monthloop : do m = 1,12

        t = m + 12 * (yr - 1)

        bcond_nd   = [ndbuf(-w),ndbuf(w)]        ! Set boundary conditions for variables
        bcond_tmin = [mtminbuf(-w),mtminbuf(w)]  ! Set boundary conditions for variables
        bcond_tmax = [mtmaxbuf(-w),mtmaxbuf(w)]  ! Set boundary conditions for variables
        bcond_cld  = [cldbuf(-w),cldbuf(w)]      ! Set boundary conditions for variables
        bcond_wnd  = [wndbuf(-w),wndbuf(w)]      ! Set boundary conditions for variables

        ! generate pseudo-daily smoothed meteorological variables (using means-preserving algorithm)

        ! call rmsmooth(mtminbuf,ndbuf,bcond_tmin,tmin_sm(1:sum(ndbuf)))
        ! call rmsmooth(mtmaxbuf,ndbuf,bcond_tmax,tmax_sm(1:sum(ndbuf)))
        ! call rmsmooth(cldbuf,ndbuf,bcond_cld,cld_sm(1:sum(ndbuf)))
        ! call rmsmooth(wndbuf,ndbuf,bcond_wnd,wnd_sm(1:sum(ndbuf)))

        call newspline_all(mtminbuf,ndbuf,tmin_sm(1:sum(ndbuf)))
        call newspline_all(mtmaxbuf,ndbuf,tmax_sm(1:sum(ndbuf)))
        call newspline_all(cldbuf,ndbuf,cld_sm(1:sum(ndbuf)))
        call newspline_all(wndbuf,ndbuf,wnd_sm(1:sum(ndbuf)))

        ! calculcate start and end positons of the current month pseudo-daily buffer

        d0 = sum(ndbuf(-w:-1)) + 1
        d1 = d0 + ndbuf(0) - 1

        !print *, d0, d1

        ndm = d1 - d0 + 1

        ! restrict simulated total monthly precip to +/-10% or 1 mm of observed value

        prec_t = max(1.,0.1 * pre(i,j,t))

        i_count = 0

        !---------------------------------------------------------------------------------
        ! quality control loop calling the weathergen - this loop principally checks that
        ! the number of wet days and total precip stayed close to the input data

        qualityloop : do
          i_count = i_count + 1    ! increment iteration number

          mwetd_sim = 0
          mprec_sim = 0.

          outd = 1

          dayloop : do d = d0,d1  ! day loop

            !write(*,*)yr,m,d0, d,d1,ndbuf

            met_in%prec  = pre(i,j,t)
            met_in%wetd  = wet(i,j,t)
            met_in%wetf  = wetf(t)
            met_in%tmin  = tmin_sm(d)
            met_in%tmax  = tmax_sm(d)
            met_in%cldf  = real(cld_sm(d))
            met_in%wind  = real(wnd_sm(d))
            met_in%pday  = met_out%pday
            met_in%resid = met_out%resid

            call weathergen(met_in,met_out)

            ! write(0,*) met_out%cldf ! in fractional form 0 to 1

            met_in%rndst = met_out%rndst
            month_met(outd) = met_out    ! save this day into a month holder

            if (met_out%prec > 0.) then

              mwetd_sim = mwetd_sim + 1
              mprec_sim = mprec_sim + met_out%prec

            end if

            outd = outd + 1

          end do dayloop ! day loop

          ! quality control checks

          if (pre(i,j,t) == 0.) then ! if there is no precip in this month a single iteration is ok

            pdaydiff = 0
            precdiff = 0.

            exit

          else if (i_count < 2) then

            cycle  !enforce at least two times over the month to get initial values for residuals ok

          else if (pre(i,j,t) > 0. .and. mprec_sim == 0.) then

            cycle  ! need to get at least some precip if there is some in the input data

          end if

          pdaydiff = abs(mwetd_sim - wet(i,j,t))

          precdiff = (mprec_sim - pre(i,j,t)) / pre(i,j,t)

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

        end do qualityloop

        ! end of quality control loop
        !---------------------------------------------------------------------------------

        ! adjust meteorological values to match the input means following Richardson & Wright 1984

        mtmin_sim = sum(month_met(1:ndm)%tmin) / ndm
        mtmax_sim = sum(month_met(1:ndm)%tmax) / ndm
        mcldf_sim = sum(month_met(1:ndm)%cldf) / ndm
        mwind_sim = sum(month_met(1:ndm)%wind) / ndm

        if (mprec_sim == 0.) then
          if (pre(i,j,t) > 0.) stop 'simulated monthly prec = 0 but input prec > 0'
          prec_corr = 1.
        else
          prec_corr = pre(i,j,t) / mprec_sim
        end if

        tmin_corr = mtmin(t) - mtmin_sim
        tmax_corr = mtmax(t) - mtmax_sim

        mcldf_sim = 100. * mcldf_sim    ! Convert fraction back to percentage (0 to 100) (Leo)

        if (mcldf_sim == 0.) then
          if (cld(i,j,t) > 0.) stop 'simulated monthly cloud = 0 but input cloud > 0'
          cldf_corr = 1.
        else
          cldf_corr = cld(i,j,t) / mcldf_sim
        end if

        ! if (mwind_sim == 0.) then
        !   if (wnd(i,j,t) > 0.) stop 'simulated monthly wind = 0 but input wind > 0'
        !   wind_corr = 1.
        ! else
        !   wind_corr = wnd(i,j,t) / mwind_sim
        ! end if

        !--- Replaced "stop" command for now since only tmin and tmax is wanted (Leo)
        if (mwind_sim == 0.) then
          if (wnd(i,j,t) > 0.) then
            write(0,*) 'Warning: simulated monthly wind = 0 but input wind > 0'
            wind_corr = 1.
          end if
        else
          wind_corr = wnd(i,j,t) / mwind_sim
        end if

        !---

        month_met(1:ndm)%prec = month_met(1:ndm)%prec * prec_corr
        month_met(1:ndm)%tmin = month_met(1:ndm)%tmin + tmin_corr
        month_met(1:ndm)%tmax = month_met(1:ndm)%tmax + tmax_corr
        month_met(1:ndm)%cldf = month_met(1:ndm)%cldf * cldf_corr
        month_met(1:ndm)%wind = month_met(1:ndm)%wind * wind_corr

        month_met(1:ndm)%cldf = min(max(month_met(1:ndm)%cldf,0.),1.)
        month_met(1:ndm)%wind = max(month_met(1:ndm)%wind,0.)

        month_met(1:ndm)%prec = roundto(month_met(1:ndm)%prec,1)
        month_met(1:ndm)%tmin = roundto(month_met(1:ndm)%tmin,1)
        month_met(1:ndm)%tmax = roundto(month_met(1:ndm)%tmax,1)
        month_met(1:ndm)%cldf = roundto(month_met(1:ndm)%cldf,3)
        month_met(1:ndm)%wind = roundto(month_met(1:ndm)%wind,2)


        !--- Write option #1: write out variables input mean, smoothed daily values and simulated daliy values (Leo)
        if (write_opt == 1) then      ! If argument read from start of program (Leo)
          do k = d0, d1

            l = k - sum(ndbuf(-w:-1))

            write(*,*) mtminbuf(0), tmin_sm(k), month_met(l)%tmin, &
                       mtmaxbuf(0), tmax_sm(k), month_met(l)%tmax, &
                       cldbuf(0), cld_sm(k), 100*month_met(l)%cldf, &
                       wndbuf(0), wnd_sm(k), month_met(l)%wind

          end do
        end if

        ! --- Debuggin tmin series error (Leo)

        ! start_debug = sum(ndbuf(-3:-1)) + 1
        ! end_debug = start_debug + ndbuf(0) - 1
        !
        ! write(0,*) 'yr:',yr,'month:',m,'ndm:',ndm, 'check:',size(tmin_sm(start_debug:end_debug)), size(month_met(1:ndm)%tmin)
        ! write(0,*) 'input mean:',mtminbuf(0),'smoothed:', (sum(tmin_sm(start_debug:end_debug))/ndm), &
        ! 'simulated mean:', (sum(month_met(1:ndm)%tmin)/ndm)
        !
        !
        ! do k = start_debug, end_debug
        !
        !   start_debug2 = k - sum(ndbuf(-3:-1))
        !   ! end_debug2 = start_debug2 + ndbuf(0) - 1
        !
        !   write(*,*) mtminbuf(0), tmin_sm(k), month_met(start_debug2)%tmin
        ! ! write(0,10) tmin_sm(start_debug:end_debug), month_met(1:ndm)%tmin
        ! !
        !   ! 10 format(f6.2, 1x, f6.2, 1x, f6.2)
        ! end do

        ! --- Debuggin tmin series error (Leo)

        ! start_debug = sum(ndbuf(-3:-1)) + 1
        ! end_debug = start_debug + ndbuf(0) - 1
        !
        ! write(0,*) 'yr:',yr,'month:',m,'ndm:',ndm, 'check:',size(tmax_sm(start_debug:end_debug)), size(month_met(1:ndm)%tmax)
        ! write(0,*) 'input mean:',mtmaxbuf(0),'smoothed:', (sum(tmax_sm(start_debug:end_debug))/ndm), &
        ! 'simulated mean:', (sum(month_met(1:ndm)%tmax)/ndm)
        !
        !
        ! do k = start_debug, end_debug
        !
        !   start_debug2 = k - sum(ndbuf(-3:-1))
        !   ! end_debug2 = start_debug2 + ndbuf(0) - 1
        !
        !   write(*,*) mtmaxbuf(0), tmax_sm(k), month_met(start_debug2)%tmax
        ! ! write(0,10) tmin_sm(start_debug:end_debug), month_met(1:ndm)%tmin
        ! !
        !   ! 10 format(f6.2, 1x, f6.2, 1x, f6.2)
        ! end do


        ! --- Debuggin cloud series error (Leo)

        ! start_debug = sum(ndbuf(-3:-1)) + 1
        ! end_debug = start_debug + ndbuf(0) - 1
        !
        ! write(0,*) 'yr:',yr,'month:',m,'ndm:',ndm, 'check:',size(cld_sm(start_debug:end_debug)), size(month_met(1:ndm)%cldf)
        ! write(0,*) 'input mean:',cldbuf(0),'smoothed:', (sum(cld_sm(start_debug:end_debug))/ndm), &
        ! 'simulated mean:', (sum(100*month_met(1:ndm)%cldf)/ndm)


        ! do k = start_debug, end_debug
        !
        !   start_debug2 = k - sum(ndbuf(-3:-1))
        !   ! end_debug2 = start_debug2 + ndbuf(0) - 1
        !
        !   write(*,*) cldbuf(0), cld_sm(k), 100*month_met(start_debug2)%cldf
        ! ! write(0,10) tmin_sm(start_debug:end_debug), month_met(1:ndm)%tmin
        ! !
        !   ! 10 format(f6.2, 1x, f6.2, 1x, f6.2)
        ! end do

        ! --- Debuggin cloud series error (Leo)

        ! start_debug = sum(ndbuf(-3:-1)) + 1
        ! end_debug = start_debug + ndbuf(0) - 1
        !
        ! write(0,*) 'yr:',yr,'month:',m,'ndm:',ndm, 'check:',size(wnd_sm(start_debug:end_debug)), size(month_met(1:ndm)%wind)
        ! write(0,*) 'input mean:',wndbuf(0),'smoothed:', (sum(wnd_sm(start_debug:end_debug))/ndm), &
        ! 'simulated mean:', (sum(month_met(1:ndm)%wind)/ndm)
        !
        !
        ! do k = start_debug, end_debug
        !
        !   start_debug2 = k - sum(ndbuf(-3:-1))
        !   ! end_debug2 = start_debug2 + ndbuf(0) - 1
        !
        !   write(*,*) wndbuf(0), wnd_sm(k), month_met(start_debug2)%wind
        ! ! write(0,10) tmin_sm(start_debug:end_debug), month_met(1:ndm)%tmin
        ! !
        !   ! 10 format(f6.2, 1x, f6.2, 1x, f6.2)
        ! end do


        !---


        !-----------------------------------------------------------
        ! add the current monthly values on to the smoothing buffer
        ! write(0,*)cntt,t,w,t+w+1

        mtminbuf = eoshift(mtminbuf,1,mtmin(t+w+1))
        mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(t+w+1))
        ndbuf    = eoshift(ndbuf,1,nd(t+w+1))
        cldbuf   = eoshift(cldbuf,1,cld(i,j,t+w+1))
        wndbuf   = eoshift(wndbuf,1,wnd(i,j,t+w+1))

        !-----------------------------------------------------------
        ! diagnostic output

        calyr = yr+startyr-1

        !-----------------------------------------------------------
        ! save the min and max temperature of the month and replace original value (Leo)

        tmin_sim = minval(month_met(1:ndm)%tmin)
        tmax_sim = maxval(month_met(1:ndm)%tmax)

        abs_tmin(i,j) = min(abs_tmin(i,j),tmin_sim)
        abs_tmax(i,j) = max(abs_tmax(i,j),tmax_sim)

      end do monthloop ! month loop

      ! tmin_sim = minval(month_met(1:ndm)%tmin) ! This original code will only read Dec... (Leo)
      ! tmax_sim = maxval(month_met(1:ndm)%tmax)

      ! abs_tmin(i,j) = min(abs_tmin(i,j),tmin_sim)
      ! abs_tmax(i,j) = max(abs_tmax(i,j),tmax_sim)


    end do yearloop    ! year loop

    !--- Write option #2: write out all tmin and tmax of individual cells
    if (write_opt == 2) then
      write(*,*) abs_tmin(i,j), abs_tmax(i,j)
    end if

  end do      ! columns
end do        ! rows


!---------------------------------------------------------------------
! write out calculated values

call getarg(4,outfile)

call genoutfile(outfile,id,[cntx,cnty],ofid)

call putlonlat(ofid,id,lon,lat)

ncstat = nf90_inq_varid(ofid,'abs_tmin',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(outvar(cntx,cnty))

where (abs_tmin /= -9999.)
  outvar = nint(abs_tmin / 0.1)
elsewhere
  outvar = -32768
end where

ncstat = nf90_put_var(ofid,varid,outvar)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---

ncstat = nf90_inq_varid(ofid,'abs_tmax',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (abs_tmax /= -9999.)
  outvar = nint(abs_tmax / 0.1)
elsewhere
  outvar = -32768
end where

ncstat = nf90_put_var(ofid,varid,outvar)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------
! close files

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end program gwgen_grid
