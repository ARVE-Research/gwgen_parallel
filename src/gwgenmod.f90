module gwgenmod

! Use the Makefile to compile this program

!! Program to run gwgen with gridded input, provide global lon/lat indexed list output of
!! absolute minimum and maximum temperature
! JO Kaplan, HKU, 2019; Leo O Lai, HKU, 2021

! Terminal command line: mpirun -np 18 ./src/gwgen ~/path/to/input.file startyr/calcyr ~/path/to/outfile

! List of modules that will be used and the variables within these modules that are used in this program:

use parametersmod, only : sp,dp,i4,i2,so,ndaymonth
use errormod,      only : ncstat,netcdf_err
use coordsmod,     only : index,parsecoords,calcpixels
use geohashmod,    only : geohash
use randomdistmod,
use newsplinemod,  only : newspline_all
use weathergenmod, only : metvars_in, metvars_out, weathergen,rmsmooth,roundto
use getdatamod,    only : readdata
use outputmod,     only : putlonlat,infompi
use netcdf
use mpi

implicit none

public :: gwgen

contains

! Inquire about the dimensions of the input file

!-------------------------------------------------------

subroutine gwgen(info,job,rank)

type(infompi), target    , intent(in) :: info
integer(i4), dimension(2), intent(in) :: job
integer                  , intent(in) :: rank

! Pointers for mpi info variables

character(100), pointer :: infile
character(100), pointer :: outfile
character(100), pointer :: timestring
integer(i4)   , pointer :: nproc
integer(i4)   , pointer :: t0
integer(i4)   , pointer :: nt

! IDs for file, dimensions and variables

integer :: ifid                        ! Input file ID
integer :: dimid                       ! Dimension ID
integer :: varid                       ! Variable ID
integer :: ofid                        ! output file ID

integer :: xlen                        ! length of dimension 'lat'
integer :: ylen                        ! length of dimension 'long'
integer :: ilen                        ! length of dimension 'index'
integer :: tlen                        ! length of dimension 'time'

! Allocatable arrays for longitude and latitude

real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat
real(dp),    allocatable, dimension(:)   :: time
integer(i4), allocatable, dimension(:,:) :: indx

integer :: baddata_check

integer, dimension(2) :: ll_loc
integer :: lon_loc
integer :: lat_loc
integer :: ll

! monthly input driver variables

real(sp), allocatable, dimension(:,:) :: tmp        ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:,:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:,:) :: pre        ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:,:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:,:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), allocatable, dimension(:,:) :: wnd        ! mean monthly 10m windspeed (m s-1)

! monthly derived driver variables

real(sp), allocatable, dimension(:) :: mtmin      ! maximum monthly temperature (degC)
real(sp), allocatable, dimension(:) :: mtmax      ! monthly minimum temperature (degC)
real(sp), allocatable, dimension(:) :: wetf       ! fraction of wet days in a month

! output variable

real(sp), allocatable, dimension(:) :: abs_tmin      ! absolute minimum temperature (degC)
real(sp), allocatable, dimension(:) :: abs_tmax      ! absolute maximum temperature (degC)

integer(i2), allocatable, dimension(:) :: outvar    ! Output variable for ncfile output adjusted by scale factor

real(sp) :: tmin_sim
real(sp) :: tmax_sim

! Elements to calculate current year and amount of days in current month

integer :: i_count,outd
integer :: i,k,t,d,m
integer :: d0
integer :: d1
integer :: calyr
integer :: ndm

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

integer :: t1
integer :: p0
integer :: p1
integer :: cntt

integer :: baseyr
integer :: startyr
integer :: calcyrs

integer :: endyr

integer, dimension(2) :: srt
integer, dimension(2) :: cnt

! Variables for output write option (Leo)

! character(1) :: write
! integer :: write_opt
! integer :: l

! Time variables

real :: start_time
real :: end_time


!-----------------------------------------------------------------------------------------------------------------------------------------
! program starts here

call CPU_TIME(start_time)

infile     => info%infile
outfile    => info%outfile
timestring => info%timestring
nproc      => info%nproc
t0         => info%t0
nt         => info%nt

!-----------------------------------------------------
! Read output writing option (Leo)
! Input 1 as fifth argument on terminal to write out data (Leo)
! call getarg(5, write)
! write_opt = ichar(write) - 48

! -----------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

ncstat = nf90_open(infile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)

ncstat = nf90_inq_dimid(ifid,'lon',dimid)              ! get dimension ID from dimension 'lon' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)              ! get dimension ID from dimension 'lon' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'index',dimid)             ! get dimension ID from dimension 'index' in the input file
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)   ! get dimension name and length from input file for dimension previously inquired
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'time',dimid)             ! Get dimension ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)   ! Get length of dimension 'time' and assign it to variable tlen
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read variable IDs and values

allocate(lon(xlen))       ! Allocate length to longitude array
allocate(lat(ylen))       ! Allocate length to latitude array
allocate(time(tlen))      ! Allocate length to time array
allocate(indx(xlen,ylen)) ! Allocate length to index array

ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lon)                    ! Get variable values for longitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lat)                    ! Get variable values for latitude
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"index",varid)               ! Get variable ID for index
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,indx)                    ! Get variable values for index
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"time",varid)               ! Get variable ID for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,time)                   ! Get variable values for time
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ifid,varid,'units',basedate)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)



!----------------------------------------------------
! Get the buffered start and end of the time dimension of each validpixel

calcyrs = nt / 12

cntt = 12 * (calcyrs + 2)          ! months, includes one-year buffer on either end

srt = [job(1), t0]
cnt = [job(2), nt]

allocate(tmp(cnt(1),cntt))
allocate(dtr(cnt(1),cntt))
allocate(pre(cnt(1),cntt))
allocate(wet(cnt(1),cntt))
allocate(cld(cnt(1),cntt))
allocate(wnd(cnt(1),cntt))

! calculate file start and count indices

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

!------

if (rank == 0) then
  write(0,*)startyr,calcyrs
  write(0,*)cntt,nt
  write(0,*)t0,t1
  write(0,*)p0,p1
end if

!---------------------------------------------------------------------
! read the timeseries of monthly climate

call readdata(ifid,'tmp',srt,cnt,tmp(:,p0:p1))
call readdata(ifid,'dtr',srt,cnt,dtr(:,p0:p1))
call readdata(ifid,'pre',srt,cnt,pre(:,p0:p1))
call readdata(ifid,'wet',srt,cnt,wet(:,p0:p1))
call readdata(ifid,'cld',srt,cnt,cld(:,p0:p1))
call readdata(ifid,'wnd',srt,cnt,wnd(:,p0:p1))

!---------------------------------------------------------------------
! copy first and last year into buffers at each end

if (p0 == 13) then
  tmp(:,1:12) = tmp(:,13:24)
  dtr(:,1:12) = dtr(:,13:24)
  pre(:,1:12) = pre(:,13:24)
  wet(:,1:12) = wet(:,13:24)
  cld(:,1:12) = cld(:,13:24)
  wnd(:,1:12) = wnd(:,13:24)
end if

if (cntt > p1) then
  tmp(:,p1+1:cntt) = tmp(:,p1-11:p1)
  dtr(:,p1+1:cntt) = dtr(:,p1-11:p1)
  pre(:,p1+1:cntt) = pre(:,p1-11:p1)
  wet(:,p1+1:cntt) = wet(:,p1-11:p1)
  cld(:,p1+1:cntt) = cld(:,p1-11:p1)
  wnd(:,p1+1:cntt) = wnd(:,p1-11:p1)
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
! calculate days per month
baseyr  = 1871

startyr = baseyr
endyr = startyr + calcyrs - 1

allocate(nd(cntt))

i = 1
do yr = startyr-1,endyr+1
  do mon = 1,12

    nd(i) = ndaymonth(yr,mon)

    i = i + 1

  end do
end do

!---------------------------------------------------------------------
! create a netCDF output file with the dimensions of the area of interest

allocate(abs_tmin(cnt(1)))
allocate(abs_tmax(cnt(1)))

abs_tmin = -9999.
abs_tmax = -9999.

!---------------------------------------------------------------------
! grid loop starts here

allocate(mtmin(cntt))
allocate(mtmax(cntt))
allocate(wetf(cntt))

ll = srt(1)   ! Get the current index value

do i = 1,cnt(1)

  ! Get the value of lon and lat from index dimension

  ll_loc = findloc(indx, ll)

  lon_loc = ll_loc(1)
  lat_loc = ll_loc(2)

  call ran_seed(geohash(lon(lon_loc),lat(lat_loc)),met_in%rndst)

  ll = ll + 1

  !---------------------------------------------------------------------
  ! calculate derived climate variables

  mtmin = tmp(i,:) - 0.5 * dtr(i,:)
  mtmax = tmp(i,:) + 0.5 * dtr(i,:)
  wetf  = wet(i,:) / nd

  !--- Checking bad data (Leo)

  do k = 1,cntt

    if(mtmin(k) < -273.15) then

      write(0,*) 'Messed up min temp. at time slice', k, 'with', mtmin(k), 'degC'
      write(0,*) 'tmp: ', tmp(i,k) !, 'dtr: ', dtr(i,j,k)
      write(0,*) 'cnti:', i, 'cntt:', k

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
  abs_tmin(i) = 9999.


  !---------------------------------------------------------------------
  ! prepare pseudo-daily smoothed meteorological variables
  ! initialize the smoothing buffer variables with all values from the first month of the input
  ! since we will always start in year 2, we can do:

  t = 13

  ndbuf    = nd(t-w:t+w)
  mtminbuf = mtmin(t-w:t+w)
  mtmaxbuf = mtmax(t-w:t+w)
  cldbuf   = cld(i,t-w:t+w)
  wndbuf   = wnd(i,t-w:t+w)

  met_out%pday(1) = .false.
  met_out%pday(2) = .false.
  met_out%resid = 0.

  ! start year loop

  yearloop : do yr = 2,(calcyrs+1) ! Changed to calcyrs+1 (Leo)

    ! start month loop

    monthloop : do m = 1,12

      t = m + 12 * (yr - 1)

      bcond_nd   = [ndbuf(-w),ndbuf(w)]        ! Set boundary conditions for variables
      bcond_tmin = [mtminbuf(-w),mtminbuf(w)]  ! Set boundary conditions for variables
      bcond_tmax = [mtmaxbuf(-w),mtmaxbuf(w)]  ! Set boundary conditions for variables
      bcond_cld  = [cldbuf(-w),cldbuf(w)]      ! Set boundary conditions for variables
      bcond_wnd  = [wndbuf(-w),wndbuf(w)]      ! Set boundary conditions for variables

      ! generate pseudo-daily smoothed meteorological variables (using means-preserving algorithm)

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

      prec_t = max(1.,0.1 * pre(i,t))

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

          met_in%prec  = pre(i,t)
          met_in%wetd  = wet(i,t)
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

        if (pre(i,t) == 0.) then ! if there is no precip in this month a single iteration is ok

          pdaydiff = 0
          precdiff = 0.

          exit

        else if (i_count < 2) then

          cycle  !enforce at least two times over the month to get initial values for residuals ok

        else if (pre(i,t) > 0. .and. mprec_sim == 0.) then

          cycle  ! need to get at least some precip if there is some in the input data

        end if

        pdaydiff = abs(mwetd_sim - wet(i,t))

        precdiff = (mprec_sim - pre(i,t)) / pre(i,t)

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
        if (pre(i,t) > 0.) stop 'simulated monthly prec = 0 but input prec > 0'
        prec_corr = 1.
      else
        prec_corr = pre(i,t) / mprec_sim
      end if

      tmin_corr = mtmin(t) - mtmin_sim
      tmax_corr = mtmax(t) - mtmax_sim

      mcldf_sim = 100. * mcldf_sim    ! Convert fraction back to percentage (0 to 100) (Leo)

      if (mcldf_sim == 0.) then
        if (cld(i,t) > 0.) stop 'simulated monthly cloud = 0 but input cloud > 0'
        cldf_corr = 1.
      else
        cldf_corr = cld(i,t) / mcldf_sim
      end if

      ! if (mwind_sim == 0.) then
      !   if (wnd(i,j,t) > 0.) stop 'simulated monthly wind = 0 but input wind > 0'
      !   wind_corr = 1.
      ! else
      !   wind_corr = wnd(i,j,t) / mwind_sim
      ! end if

      !--- Replaced "stop" command for now since only tmin and tmax is wanted (Leo)
      if (mwind_sim == 0.) then
        if (wnd(i,t) > 0.) then
          ! write(0,*) 'Warning: simulated monthly wind = 0 but input wind > 0'
          wind_corr = 1.
        end if
      else
        wind_corr = wnd(i,t) / mwind_sim
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
      ! if (write_opt == 1) then      ! If argument read from start of program (Leo)
      !   do k = d0, d1
      !
      !     l = k - sum(ndbuf(-w:-1))
      !
      !     write(*,*) mtminbuf(0), tmin_sm(k), month_met(l)%tmin, &
      !                mtmaxbuf(0), tmax_sm(k), month_met(l)%tmax, &
      !                cldbuf(0), cld_sm(k), 100*month_met(l)%cldf, &
      !                wndbuf(0), wnd_sm(k), month_met(l)%wind
      !
      !   end do
      ! end if

      !-----------------------------------------------------------
      ! add the current monthly values on to the smoothing buffer
      ! write(0,*)cntt,t,w,t+w+1

      mtminbuf = eoshift(mtminbuf,1,mtmin(t+w+1))
      mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(t+w+1))
      ndbuf    = eoshift(ndbuf,1,nd(t+w+1))
      cldbuf   = eoshift(cldbuf,1,cld(i,t+w+1))
      wndbuf   = eoshift(wndbuf,1,wnd(i,t+w+1))

      !-----------------------------------------------------------
      ! diagnostic output

      calyr = yr+startyr-1

      !-----------------------------------------------------------
      ! save the min and max temperature of the month and replace original value (Leo)

      tmin_sim = minval(month_met(1:ndm)%tmin)
      tmax_sim = maxval(month_met(1:ndm)%tmax)

      abs_tmin(i) = min(abs_tmin(i),tmin_sim)
      abs_tmax(i) = max(abs_tmax(i),tmax_sim)

    end do monthloop ! month loop

  end do yearloop    ! year loop

  !--- Write option #2: write out all tmin and tmax of individual cells
  ! if (write_opt == 2) then
  !   write(*,*) abs_tmin(i), abs_tmax(i)
  ! end if

  ! write(0,*) abs_tmin(i), abs_tmax(i)

end do      ! columns

!---------------------------------------------------------------------
! write out calculated values

! call putlonlat(ofid,id,lon,lat)

allocate(outvar(cnt(1)))

ncstat = nf90_open(outfile,nf90_write,ofid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'abs_tmin',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (abs_tmin /= -9999.)
  outvar = nint(abs_tmin / 0.1)
elsewhere
  outvar = -32768
end where

ncstat = nf90_put_var(ofid,varid,outvar,start=srt,count=cnt)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---

ncstat = nf90_inq_varid(ofid,'abs_tmax',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

where (abs_tmax /= -9999.)
  outvar = nint(abs_tmax / 0.1)
elsewhere
  outvar = -32768
end where

ncstat = nf90_put_var(ofid,varid,outvar,start=srt,count=cnt)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------
! close files

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!------

call CPU_TIME(end_time)

write(0,*) 'Rank:', rank, 'min/max',minval(abs_tmin), maxval(abs_tmax), 'using seconds', end_time-start_time



end subroutine gwgen


end module gwgenmod
