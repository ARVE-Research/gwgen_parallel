module outputmod

use parametersmod, only : i2,i4,sp,dp
use errormod,      only : ncstat,netcdf_err
use netcdf

implicit none

public :: getoutfile
public :: putlonlat

real(sp),    parameter :: missing_sp = -9999.
integer(i2), parameter :: missing_i2 = -32768

!--------------------

type infompi
  character(100) :: infile
  character(100) :: outfile
  character(100) :: timestring
  integer(i4)    :: nproc
  integer(i4)    :: validcell
  integer(i4)    :: t0
  integer(i4)    :: nt
end type

!--------------------

contains

!-------------------------------------------------------------------------------------------------

subroutine getoutfile(info)

type(infompi), target, intent(in) :: info

character(100), pointer :: infile
character(100), pointer :: outfile
integer(i4) ,   pointer :: validcell

integer :: ifid
integer :: ofid
integer :: ncstat
integer :: varid
integer :: dimid

integer, dimension(3) :: dimids

integer :: xlen
integer :: ylen
integer :: ilen

real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat
integer(i4), allocatable, dimension(:,:) :: indx

integer :: i

character(8) :: today
character(10) :: now

!--------------------

infile    => info%infile
outfile   => info%outfile
validcell => info%validcell

!--------------------

write(0,*)
write(0,*) 'Reading lon, lat and index from infile'

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

ncstat = nf90_inq_dimid(ifid,'index',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!------

allocate(lon(xlen))
allocate(lat(ylen))
allocate(indx(xlen,ylen))

ncstat = nf90_inq_varid(ifid,"lon",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lon)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"lat",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,lat)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ifid,"index",varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ifid,varid,indx)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

write(0,*) 'Got lon, lat and index from infile'
write(0,*)

!--------------------

write(0,*) 'Creating outfile'

ncstat = nf90_create(outfile,nf90_netcdf4,ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_create_par(outfile,ior(nf90_netcdf4,nf90_mpiio),MPI_COMM_WORLD,MPI_INFO_NULL,ofid)
! if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','weathergen parallel output file')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'created',today//' '//now(1:4))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! lon

ncstat = nf90_def_dim(ofid,'lon',xlen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(1) = dimid

ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lon)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! lat

ncstat = nf90_def_dim(ofid,'lat',ylen,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(2) = dimid

ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lat)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----
! index

ncstat = nf90_def_dim(ofid,'index',validcell,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

dimids(3) = dimid

ncstat = nf90_def_var(ofid,'index',nf90_int,dimids(1:2),varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,indx)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','index of lon and lat')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','1 to length')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',-1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)


!----
!absolute minimum temperature

ncstat = nf90_def_var(ofid,'abs_tmin',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','test average temperature')
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
!absolute maximum temperature

ncstat = nf90_def_var(ofid,'abs_tmax',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','test average temperature')
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

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

write(0,*) 'Finished create outfile'
write(0,*)


end subroutine getoutfile

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
