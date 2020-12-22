module coordsmod

use parametersmod, only : dp,i4

implicit none

public :: parsecoords
public :: calcpixels

type index
  real(dp)    :: minlon
  real(dp)    :: maxlon
  real(dp)    :: minlat
  real(dp)    :: maxlat
  integer(i4) :: startx
  integer(i4) :: starty
  integer(i4) :: endx
  integer(i4) :: endy
  integer(i4) :: countx
  integer(i4) :: county
end type index

contains

!-----------------------------------------------------------------------------------------------

subroutine parsecoords(coordstring,id)

!subroutine to parse a coordinate string

use parametersmod, only : dp

implicit none

!arguments

character(*), intent(in)  :: coordstring
type(index),  intent(out) :: id

!local variables

real(dp),      allocatable, dimension(:) :: val
character(50), allocatable, dimension(:) :: cval

integer :: i
integer :: last
integer :: part
integer :: np

!----
! takes the input string minus whitespace, converts it to an array of characters and counts occurrences of "/"
! number of parts = 1 + the number of "/"

np = 1 + count(transfer(coordstring,'a',len_trim(coordstring)) == '/')

allocate(cval(np))
allocate(val(np))

last = 1
part = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(last:i-1)
    last = i + 1
    part = part + 1
  end if
end do

cval(part) = coordstring(last:i-1)

read(cval,*)val

if (np == 2) then
  id%minlon = val(1)
  id%maxlon = val(1)
  id%minlat = val(2)
  id%maxlat = val(2)
else
  id%minlon = val(1)
  id%maxlon = val(2)
  id%minlat = val(3)
  id%maxlat = val(4)
end if

end subroutine parsecoords

!-----------------------------------------------------------------------------------------------

subroutine calcpixels(lon,lat,id)

implicit none

! arguments

real(dp), dimension(:), intent(in) :: lon
real(dp), dimension(:), intent(in) :: lat

type(index),  intent(inout) :: id

integer, dimension(2) :: xpos
integer, dimension(2) :: ypos

! local variables

integer, dimension(1) :: pos

!-----------------------------------------------------

pos = minloc(abs(lon - id%minlon))
xpos(1) = pos(1)

pos = minloc(abs(lon - id%maxlon))
xpos(2) = pos(1)

pos = minloc(abs(lat - id%minlat))
ypos(1) = pos(1)

pos = minloc(abs(lat - id%maxlat))
ypos(2) = pos(1)

id%startx = minval(xpos)
id%starty = minval(ypos)

if (id%minlon == id%maxlon .and. id%minlat == id%maxlat) then  !special case, run just one nearest gridcell even if coords were ambiguous
  
  id%countx = 1
  id%county = 1
  
else

  if (lon(id%startx) < id%minlon) id%startx = id%startx + 1
  id%countx = 1 + abs(maxval(xpos) - id%startx)

  if (lat(id%starty) < id%minlat) id%starty = id%starty + 1
  id%county = 1 + abs(maxval(ypos) - id%starty)

end if

if (xpos(2) - xpos(1) + 1 > id%countx) xpos(2) = xpos(2) - 1
if (ypos(2) - ypos(1) + 1 > id%county) ypos(2) = ypos(2) - 1

end subroutine calcpixels

!-----------------------------------------------------------------------------------------------

end module coordsmod
