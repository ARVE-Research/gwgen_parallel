module coordsmod

use parametersmod, only : dp

implicit none

public :: parsecoords
public :: calcpixels

character(45) :: coordstring
real(dp), dimension(4) :: bounds

contains

!-----------------------------------------------------------------------------------------------

subroutine parsecoords(coordstring,val)

use parametersmod, only : dp

!subroutine to parse a coordinate string

implicit none

character(45),          intent(in)  :: coordstring
real(dp), dimension(4), intent(out) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

end subroutine parsecoords

!-----------------------------------------------------------------------------------------------

subroutine calcpixels(lon,lat,bounds,xpos,ypos,srtx,srty,cntx,cnty)

implicit none

! arguments

real(dp), dimension(:), intent(in) :: lon
real(dp), dimension(:), intent(in) :: lat
real(dp), dimension(4), intent(in) :: bounds

integer, dimension(2), intent(out) :: xpos
integer, dimension(2), intent(out) :: ypos
integer, intent(out) :: srtx
integer, intent(out) :: srty
integer, intent(out) :: cntx
integer, intent(out) :: cnty

! local variables

integer, dimension(1) :: pos

!-----------------------------------------------------

pos = minloc(abs(lon - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lon - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(lat - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(lat - bounds(4)))
ypos(2) = pos(1)

srtx = minval(xpos)
srty = minval(ypos)

if (bounds(1) == bounds(2) .and. bounds(3) == bounds(4)) then  !special case, run just one nearest gridcell even if coords were ambiguous
  
  cntx = 1
  cnty = 1
  
else

  if (lon(srtx) < bounds(1)) srtx = srtx + 1
  cntx = 1 + abs(maxval(xpos) - srtx)

  if (lat(srty) < bounds(3)) srty = srty + 1
  cnty = 1 + abs(maxval(ypos) - srty)

end if

if (xpos(2) - xpos(1) + 1 > cntx) xpos(2) = xpos(2) - 1
if (ypos(2) - ypos(1) + 1 > cnty) ypos(2) = ypos(2) - 1

end subroutine calcpixels

!-----------------------------------------------------------------------------------------------

end module coordsmod
