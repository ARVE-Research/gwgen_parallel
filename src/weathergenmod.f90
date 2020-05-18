module weathergenmod

! This module includes subroutines to calculate daily maximum and minimum temperature and cloud cover fraction
! based on an annual timeseries of monthly values of these variables
! The weather generator is based on the WGEN model (Richardson, 1981) with extension to use monthly summary
! data from Geng et al., 1986, and Geng and Auburn 1986.
! Additional statistical relationships for both temperature, cloudiness and wind speed have been produced
! by P.S. Sommer and J.O. Kaplan using global weather station datasets (GHCN and global synoptic cloud
! reports).
!
! Coded in 2007-2009 by Jed Kaplan and Joe Melton, ARVE Group, EPFL/UVic
! 2011, Shawn Koppenhoefer, EPFL
! 2016, Philipp Sommer, UNIL
! 2020, Jed Kaplan, HKU, jed.kaplan@hku.hk

use parametersmod, only : sp,dp,i4,hsp
use randomdistmod, only : randomstate

implicit none

public  :: metvars_in
public  :: metvars_out
public  :: rmsmooth
public  :: weathergen

private :: daymetvars
private :: meansd

!-------------------------------

type metvars_in
  ! Derived datatype for the monthly weather generator input

  real(sp) :: prec    ! monthly total precipitation amount (mm)
  real(sp) :: wetd    ! number of days in month with precipitation
  real(sp) :: wetf    ! fraction of days in month with precipitation

  real(sp) :: tmin    ! minumum temperture (C)
  real(sp) :: tmax    ! maximum temperture (C)
  real(sp) :: cldf    ! cloud fraction (0=clear sky, 1=overcast) (fraction)
  real(sp) :: wind    ! wind speed (m/s)

  logical, dimension(2)  :: pday    !precipitation status: true if the day was a rain day
  type(randomstate)      :: rndst   !state of the random number generator
  real(sp), dimension(4) :: resid   !previous day's weather residuals

end type metvars_in

type metvars_out
  ! Derived datatype for the daily weather generator output

  real(sp) :: prec    ! 24 hour total precipitation (mm)
  real(sp) :: tmin    ! 24 hour mean minimum temperature (degC)
  real(sp) :: tmax    ! 24 hour mean maximum temperature (degC)
  real(sp) :: cldf    ! 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind    ! wind speed (m s-1)

  logical, dimension(2)  :: pday    !precipitation state
  type(randomstate)      :: rndst   !state of the random number generator, 15 elements
  real(sp), dimension(4) :: resid   !previous day's weather residuals
  real(sp) :: wind_bias
  real(sp) :: wind_intercept_bias
  real(sp) :: tmin_mn
  real(sp) :: tmin_sd
  real(sp) :: wind_mn
  real(sp) :: wind_sd
  real(sp), dimension(4) :: unorm

end type metvars_out

type daymetvars
  ! Derived datatype for monthly climate variables

  real(sp) :: tmax_mn     ! maximum temperature monthly mean (degC)
  real(sp) :: tmin_mn     ! minimum temperature mothly mean (degC)
  real(sp) :: cldf_mn     ! mean cloud fraction (fraction)
  real(sp) :: wind_mn     ! wind speed

  real(sp) :: tmax_sd     ! standard deviation of corresponding variable above
  real(sp) :: tmin_sd     ! "
  real(sp) :: cldf_sd     ! "
  real(sp) :: wind_sd     ! "

end type daymetvars

! -----------------------------------------------------------------------------
! parameters

integer, parameter, dimension(6) :: exponents = [0,1,2,3,4,5]

real(sp), parameter :: p_trans = 5.   ! Threshold for transition from the gamma to gp distribution (mm)

! coefficient to esimate the gamma scale parameter via
! g_scale = g_scale_coeff * mean_monthly_precip / number_of_wet_days
! following Geng et al., 1986
real(sp), parameter :: g_scale_coeff = 1.268022 ! coefficient to esimate the gamma scale parameter

real(sp), parameter :: gp_shape = 1.5 ! shape parameter for the Generalized Pareto distribution

! "A" matrix of present-day cross correlations among met variables (Richardson, 1984; Eqn. 4)
real(sp),  dimension(4,4) :: A = reshape( &
      [0.913437,  0.032532, -0.020658,  0.000573,   &
       0.488761,  0.137304, -0.072510, -0.046058,   &
      -0.001990, -0.045730,  0.591761,  0.026439,   &
       0.010905, -0.044171, -0.018568,  0.666672],  & 
[4,4])

! "B" matrix of lag-1 cross correlations among met variables (Richardson, 1984; Eqn. 4)
real(sp),  dimension(4,4) :: B = reshape( & 
      [0.361854,  0.      ,  0.      ,  0.      ,   &
       0.114410,  0.802987,  0.      ,  0.      ,   &
       0.144862, -0.060622,  0.782791,  0.      ,   &
       0.080593, -0.015829,  0.066186,  0.736713],  &
[4,4])

! transition probability correlations
real(sp), parameter :: p11_1  = 0.254877  ! intercept of p11 best line fit
real(sp), parameter :: p11_2  = 0.745123  ! slope of p11 best line fit
real(sp), parameter :: p101_1 = 0.000000  ! intercept of p101 best line fit
real(sp), parameter :: p101_2 = 0.846326  ! slope of p101 best line fit
real(sp), parameter :: p001_1 = 0.000000  ! intercept of p001 best line fit
real(sp), parameter :: p001_2 = 0.724019  ! slope of p001 best line fit

! temperature and cloud correlation parameters corresponding to wet or dry day

! minimum temperature regression results
real(sp), parameter :: tmin_w1 =  1.164653  ! intercept of best line fit of tmin on wet days
real(sp), parameter :: tmin_w2 =  0.955787  ! slope of best line fit of tmin on wet days
real(sp), parameter :: tmin_d1 = -0.528308  ! intercept of best line fit of tmin on dry days
real(sp), parameter :: tmin_d2 =  1.020964  ! slope of best line fit of tmin on dry days

real(sp), parameter, dimension(4) :: tmin_sd_breaks = [ -40., 0., 25., hsp ]  ! breaks of the tmin sd correlation

! polynomial coefficients for correlating tmin sd on wet days
real(sp), parameter, dimension(4,6) :: tmin_sd_w = reshape(     &
   ! < -40      -40 - 0         0 - 25          > 25
   [9.727157,   3.054988,       3.218742,       0.5570704,      &    
    0.1010504, -0.2115882,     -4.5076340E-02,  2.4431231E-02,  &
    0.0000000,  1.3749480E-02,  2.0944821E-02,  0.0000000,      &
    0.0000000,  1.4053800E-03, -2.6457701E-03,  0.0000000,      &
    0.0000000,  3.6860001E-05,  9.8179997E-05,  0.0000000,      &
    0.0000000,  3.2000000E-07, -1.1300000E-06,  0.0000000],     &
[4,6])

! polynomial coefficients for correlating tmin sd on dry days
real(sp), parameter, dimension(4,6) :: tmin_sd_d = reshape(     &
   ! < -40      -40 - 0         0 - 25          > 25
  [10.89901,    3.567557,       3.794117,      -4.619434,       &
    0.1270989, -0.1154459,      3.2986969E-02,  0.2260560,      &
    0.0000000,  2.8244009E-02, -1.5045540E-02,  0.0000000,      &
    0.0000000,  1.9561199E-03,  1.9034600E-03,  0.0000000,      &
    0.0000000,  4.3140000E-05, -1.1362000E-04,  0.0000000,      &
    0.0000000,  3.2000000E-07,  2.1300000E-06,  0.0000000],     &
[4,6])

! maximum temperature regression results
real(sp), parameter :: tmax_w1 = -0.586296  ! intercept of best line fit of tmax on wet days
real(sp), parameter :: tmax_w2 =  0.948669  ! slope of best line fit of tmax on wet days
real(sp), parameter :: tmax_d1 =  0.386508  ! intercept of best line fit of tmax on dry days
real(sp), parameter :: tmax_d2 =  1.0061    ! slope of best line fit of tmax on dry days

real(sp), parameter, dimension(4) :: tmax_sd_breaks = [ -30., 0., 35., hsp ]  ! breaks of the tmax sd correlation

! polynomial coefficients for correlating tmax sd on wet days
real(sp), parameter, dimension(4,6) :: tmax_sd_w = reshape(      &
! < -30           -30 - 0         0 - 35          > 35
 [6.672004,       3.860109,       3.791932,       5.552928,      &
  3.6439080E-02, -0.2186120,     -3.1260211E-02, -9.7347148E-02, &
  0.0000000,      3.8846501E-03,  1.6114730E-02,  0.0000000,     &
  0.0000000,      1.4617400E-03, -1.2029801E-03,  0.0000000,     &
  0.0000000,      6.0589999E-05,  2.9119999E-05,  0.0000000,     &
  0.0000000,      7.4000002E-07, -2.3999999E-07,  0.0000000],    &
[4,6])

! polynomial coefficients for correlating tmax sd on dry days
real(sp), parameter, dimension(4,6) :: tmax_sd_d = reshape(      &
! < -30           -30 - 0         0 - 35          > 35
 [7.374552,       4.617019,       4.745510,       3.255418,      &
  1.5355260E-02, -0.3387282,     -7.6098159E-02, -2.1786051E-02, &
  0.0000000,     -1.8756600E-02,  1.8930580E-02,  0.0000000,     &
  0.0000000,     -3.1850001E-04, -1.3494299E-03,  0.0000000,     &
  0.0000000,      3.5000000E-06,  3.2090000E-05,  0.0000000,     &
  0.0000000,      1.1000000E-07, -2.5000000E-07,  0.0000000],    &
[4,6])

! cloud regression results
real(sp), parameter :: cldf_w    = -0.738271  ! *a* parameter for cloud fit on wet days
real(sp), parameter :: cldf_d    =  0.420534  ! *a* parameter for cloud fit on dry days
real(sp), parameter :: cld_sd_w1 =  0.981917  ! *a* parameter for std. dev. of cloud fit on wet days
real(sp), parameter :: cld_sd_d1 =  1.041732  ! *a* parameter for std. dev. of cloud fit on dry days

real(sp), parameter :: cldf_w1   =  -cldf_w - 1.
real(sp), parameter :: cldf_w2   =   cldf_w * cldf_w
real(sp), parameter :: cldf_w3   = -(cldf_w * cldf_w) - cldf_w
real(sp), parameter :: cldf_w4   = -1. / cldf_w
real(sp), parameter :: cldf_sd_w = cld_sd_w1**2

real(sp), parameter :: cldf_d1   =  -cldf_d - 1.
real(sp), parameter :: cldf_d2   =   cldf_d * cldf_d
real(sp), parameter :: cldf_d3   = -(cldf_d * cldf_d) - cldf_d
real(sp), parameter :: cldf_d4   = -1./cldf_d
real(sp), parameter :: cldf_sd_d = cld_sd_d1**2

! wind regression results
real(sp), parameter :: wind_w1 = 0.          ! intercept of best line fit of wind on wet days
real(sp), parameter :: wind_w2 = 1.092938     ! slope of best line fit of wind on wet days
real(sp), parameter :: wind_d1 = 0.          ! intercept of best line fit of wind on dry days
real(sp), parameter :: wind_d2 = 0.945229     ! slope of best line fit of wind on wet days

! polygon coefficients for wind standard deviation on wet days
real(sp), parameter, dimension(6) :: wind_sd_w = [ 0., 0.81840997, -0.12633931, 0.00933591, 0., 0. ]

! polygon coefficients for wind standard deviation on dry days
real(sp), parameter, dimension(6) :: wind_sd_d = [ 0., 1.08596114, -0.24073323, 0.02216454, 0., 0. ]

! bias correction for wind speed (section 3.4)
! min. and max range for bias correction (1st and 99th percentile)
real(sp), parameter :: wind_bias_min = -2.3263478740
real(sp), parameter :: wind_bias_max =  2.3263478740

! coefficients for the polynomial slope correction (eqn. 22)
real(sp), parameter, dimension(6) :: wind_bias_coeffs = [ 0.9953539, 0.8507947, 0.027799824, -0.067101441, 0., 0. ]

! coefficients for the exponential intercept correction (eqn. 21)
real(sp), parameter :: wind_intercept_bias_a =  1.158225  ! slope in the exponent
real(sp), parameter :: wind_intercept_bias_b = -1.335892  ! intercept in the exponent

! coefficients for the polynomial intercept correction (shown in figure 14)

real(sp), parameter, dimension(6) :: wind_intercept_bias_coeffs = [0.173,0.437,0.290,0.060,0.,0.]

! -----------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

subroutine weathergen(met_in,met_out)

use parametersmod, only : sp,dp,tfreeze
use randomdistmod, only : ranur,ran_normal,ran_gamma_gp,ran_gamma,gamma_cdf,gamma_pdf

implicit none

!---------------
!arguments

type(metvars_in),  intent(in)  :: met_in
type(metvars_out), intent(out) :: met_out

!---------------
!local variables

integer  :: i

real(sp) :: pre    ! monthly total precipitation amount (mm)
real(sp) :: wetd   ! number of days in month with precipitation (fraction)
real(sp) :: wetf   ! fraction of days in month with precipitation (fraction)
real(sp) :: tmn    ! minumum temperture (C)
real(sp) :: tmx    ! maximum temperture (C)
real(sp) :: cld    ! cloud fraction (0=clear sky, 1=overcast) (fraction)
real(sp) :: wnd    ! wind (m/s)

real(sp), pointer :: tmax_mn
real(sp), pointer :: tmin_mn
real(sp), pointer :: cldf_mn
real(sp), pointer :: wind_mn
real(sp), pointer :: tmax_sd
real(sp), pointer :: tmin_sd
real(sp), pointer :: cldf_sd
real(sp), pointer :: wind_sd

type(randomstate) :: rndst       ! integer state of the random number generator
logical,  dimension(2) :: pday   ! element for yesterday and the day before yesterday
real(sp), dimension(4) :: resid  ! previous day's weather residuals

real(sp) :: prec
real(sp) :: tmin
real(sp) :: tmax
real(sp) :: cldf
real(sp) :: wind

real(sp) :: pbar     ! mean amount of precipitation per wet day (mm)
real(sp) :: pwet     ! probability that today will be wet
real(sp) :: u        ! uniformly distributed random number (0-1)

real(sp) :: g_shape
real(sp) :: g_scale
real(sp) :: gp_scale

! bias correction
real(sp) :: slopecorr      ! slope correction for wind
real(sp) :: intercept_corr ! intercept correction for wind

real(dp) :: cdf_thresh  ! gamma cdf at the threshold
real(dp) :: pdf_thresh  ! gamma pdf at the threshold

type(daymetvars), target :: dmetvars

real(sp), dimension(4) :: unorm  ! vector of uniformly distributed random numbers (0-1)

!---------------------------------------------------------
!input

pre   = met_in%prec
wetd  = met_in%wetd
wetf  = met_in%wetf
tmn   = met_in%tmin
tmx   = met_in%tmax
cld   = met_in%cldf
wnd   = met_in%wind
rndst = met_in%rndst
pday  = met_in%pday
resid = met_in%resid

!shorthand to mean and CV structure

tmin_mn => dmetvars%tmin_mn
tmax_mn => dmetvars%tmax_mn
cldf_mn => dmetvars%cldf_mn
wind_mn => dmetvars%wind_mn
tmin_sd => dmetvars%tmin_sd
tmax_sd => dmetvars%tmax_sd
cldf_sd => dmetvars%cldf_sd
wind_sd => dmetvars%wind_sd

! do i = 1,4
!   write(0,*)a(:,i)
! end do
! write(0,*)
! do i = 1,4
!   write(0,*)b(:,i)
! end do
! read(*,*)

! these values come from one of the old namelists

    a(:, 1) = [0.9161179, 0.48541803, 0.00353546, 0.01237203]
    a(:, 2) = [0.03119116, 0.13461942, -0.04303678, -0.04305549]
    a(:, 3) = [-0.01836789, -0.06865181, 0.59195356, -0.02017241]
    a(:, 4) = [0.00087781, -0.04650604, 0.0234441, 0.67185901]

    b(:, 1) = [0.35827493, 0.11247485, 0.14180909, 0.07727885]
    b(:, 2) = [0.0, 0.80889378, -0.06030167, -0.01563132]
    b(:, 3) = [0.0, 0.0, 0.78463208, 0.06061359]
    b(:, 4) = [0.0, 0.0, 0.0, 0.73288315]



!-----------------------------------------------------------
!1) Precipitation occurrence

!if there is precipitation this month, calculate the precipitation state for today

if (wetf > 0. .and. pre > 0.) then

  !calculate transitional probabilities for dry to wet and wet to wet days
  !Relationships from Geng & Auburn, 1986, Weather simulation models based on summaries of long-term data

  if (pday(1)) then !yesterday was raining, use p11

    pwet = p11_1 + p11_2 * wetf

  else if (pday(2)) then !yesterday was not raining but the day before yesterday was raining, use p101

    pwet = p101_1 + p101_2 * wetf

  else  !both yesterday and the day before were dry, use p001

    pwet = p001_1 + p001_2 * wetf

  end if

  ! -----
  ! determine the precipitation state of the current day using the Markov chain approach
  
  u = ranur(rndst)
  
  if (u <= pwet) then  !today is a rain day

    pday = eoshift(pday,-1,.true.)

  else  !today is dry

    pday = eoshift(pday,-1,.false.)

  end if

  !---------------------------
  !2) precipitation amount

  if (pday(1)) then  !today is a wet day, calculate the rain amount

    !calculate parameters for the distribution function of precipitation amount

    pbar = pre / wetd
    
    g_scale = g_scale_coeff * pbar
    g_shape = pbar / g_scale
        
    call gamma_cdf(p_trans,0.,g_scale,g_shape,cdf_thresh)

    call gamma_pdf(p_trans,0.,g_scale,g_shape,pdf_thresh)

    gp_scale = (1. - cdf_thresh) / pdf_thresh

    i = 1
    do

      !today's precipitation

      prec = ran_gamma_gp(rndst,.true.,g_shape,g_scale,p_trans,gp_shape,gp_scale)

      !simulated precipitation should have no more precision than the input (0.1 mm day-1)

      prec = roundto(prec,1)    

      ! enforce positive precipitation that is not more than 5% greater than the monthly total

      if (prec > 0. .and. prec <= 1.05 * pre) exit

      if (i == 1000) then
        write (0,*)'Could not find good precipitation with ', pre, ' mm and ', wetd, ' wet days'
        stop
      else
        i = i + 1
      end if

    end do

  else

    prec = 0.

  end if

else  ! there was no precipitation in this month, so none on this day either

  pday = .false.
  prec = 0.

end if

!---------------------------
!3) temperature min and max, cloud fraction

! calculate a baseline mean and SD for today's weather dependent on precip status

call meansd(pday(1),tmn,tmx,cld,wnd,dmetvars)

! draw four random numbers from a normal distribution

do i = 1,4
  call ran_normal(rndst,unorm(i))
end do 

! calculate today's residuals for weather variables
! this captures the cross correlation between variables and within and internal temporal autocorrelation

resid = matmul(A,resid) + matmul(B,unorm)  !Richardson 1981, eqn 5; WGEN tech report eqn. 3
! resid = matmul(B,unorm)  !Richardson 1981, eqn 5; WGEN tech report eqn. 3
! resid = matmul(A,unorm)  !Richardson 1981, eqn 5; WGEN tech report eqn. 3

! calculate the weather variables as a function of the mean and standard deviation modified by the residual

!----- 
! minimum temperature

tmin = resid(1) * tmin_sd + tmin_mn

if (abs(tmin - tmin_mn) > 2.* tmin_sd) then

  write(0,'(l3,6f10.4)')pday(1),tmn,tmin_mn,tmin_sd,resid(1),unorm(1),tmin
!   read(*,*)

end if

!----- 
! maximum temperature

tmax = resid(2) * tmax_sd + tmax_mn

!----- 
! cloud fraction

cldf = resid(3) * cldf_sd + cldf_mn

!----- 
! wind speed and wind bias correction

wind = (resid(4) * sqrt(max(0., wind_sd)) + sqrt(max(0., wind_mn)))**2

slopecorr = sum(wind_bias_coeffs * (max(wind_bias_min,min(wind_bias_max,resid(4)))**exponents))

! The intercept can be calculated using either an exponential (eqn 21) or a polynomial (fig. 14) function. Choose one from below

! intercept_corr = exp(wind_intercept_bias_b + wind_intercept_bias_a * max(wind_bias_min,min(wind_bias_max,resid(4))))
intercept_corr = sum(wind_intercept_bias_coeffs * (max(wind_bias_min,min(wind_bias_max,resid(4)))**exponents))

wind = (wind - intercept_corr) / max(slopecorr,9.e-4)

!-----
! check and correct invalid values

if (tmin+Tfreeze < 0.) then
  write(0,*)'Unphysical min. temperature with ',tmin,'K from a monthly mean ',tmin_mn,'degC'
  stop
end if

if (tmax+Tfreeze < 0.) then
  write(0,*)'Unphysical max. temperature with ',tmax,'K from a monthly mean ',tmax_mn,'degC'
  stop
end if

wind = max(0.,wind)

cldf = min(max(cldf,0.),1.)

!----- 
! adjustment to input precision

tmin = roundto(tmin,1)
tmax = roundto(tmax,1)
cldf = roundto(cldf,3)
wind = roundto(wind,2)

!---

met_out%prec  = prec
met_out%tmin  = tmin
met_out%tmax  = tmax
met_out%cldf  = cldf
met_out%wind  = wind
met_out%pday  = pday
met_out%rndst = rndst
met_out%resid = resid
met_out%tmin_mn = tmin_mn
met_out%tmin_sd = tmin_sd
met_out%wind_bias = slopecorr
met_out%wind_intercept_bias = intercept_corr
met_out%wind_mn = wind_mn
met_out%wind_sd = wind_sd
met_out%unorm = unorm

end subroutine weathergen

!------------------------------------------------------------------------------------------------------------

subroutine meansd(pday,tmn,tmx,cld,wind,dm)

! Adjust the monthly means of temperature, cloud and wind corresponding to the wet/dry state
!
! This routine makes the first approximation inside the weather generator to adjust the monthly
! mean according to the wet/dry state using the best fit lines from the parameterization.
!
! Min. and max. temperature, as well as the wind speed, are calculated via
!
! .. math::
!
!     x_{w/d} = x_{w/d1} + x_{w/d2} \cdot \bar{x}
!
! Where :math:`x` stands either for the :math:`T_{min}, T_{max}, T_{min, sd}, T_{max, sd}, wind`
! or :math:`wind_{sd}`. :math:`w/d` stands for the wet dry state deterimined by `pday`.
!
! The cloud fraction is calculated via
!
! .. math::
!
!     c_{w/d} = \frac{-a_{w/d} - 1}{a_{w/d}^2 * \bar{c} - a_{w/d}^2 - a_{w/d}}  - \frac{1}{a_{w/d}}
!
! and it's standard deviation via
!
! .. math::
!
!     c_{sd, w/d} = a_{sd, w/d}^2 \cdot c_{w/d} \cdot (1 - c_{w/d})

implicit none

logical,          intent(in)  :: pday    ! precipitation status (mm/day)
real(sp),         intent(in)  :: tmn     ! smooth interpolation of monthly minimum temperature (degC)
real(sp),         intent(in)  :: tmx     ! smooth interpolation of monthly maximum temperature (degC)
real(sp),         intent(in)  :: cld     ! fraction (0-1)
real(sp),         intent(in)  :: wind    ! wind speed (m/s)
type(daymetvars), intent(out) :: dm      ! the :f:type:`daymetvars` for the first daily approximation

!local variables

!---

if (pday) then  !calculate mean and SD for a wet day

  dm%tmin_mn = tmin_w1 + tmin_w2 * tmn

  dm%tmax_mn = tmax_w1 + tmax_w2 * tmx

  dm%wind_mn = wind_w1 + wind_w2 * wind

  dm%cldf_mn = cldf_w1 / (cldf_w2 * cld + cldf_w3) + cldf_w4

  dm%wind_sd = sum(wind_sd_w * (dm%wind_mn ** [ 0, 1, 2, 3, 4, 5 ]))

  dm%cldf_sd = cldf_sd_w * dm%cldf_mn * (1. - dm%cldf_mn)

else  !dry day

  dm%tmin_mn = tmin_d1 + tmin_d2 * tmn

  dm%tmax_mn = tmax_d1 + tmax_d2 * tmx

  dm%wind_mn = wind_d1 + wind_d2 * wind

  dm%cldf_mn = cldf_d1 / (cldf_d2 * cld + cldf_d3) + cldf_d4

  dm%wind_sd = sum(wind_sd_d * (dm%wind_mn ** [ 0, 1, 2, 3, 4, 5 ]))

  dm%cldf_sd = cldf_sd_d * dm%cldf_mn * (1. - dm%cldf_mn)

end if

call temp_sd(pday, dm)

end subroutine meansd

!----------------------------------------------------------------------------------------------------------------

subroutine rmsmooth(m,dmonth,bcond,r)

! Iterative, mean preserving method to smoothly interpolate mean data to pseudo-sub-timestep values
! From Rymes, M.D. and D.R. Myers, 2001. Solar Energy (71) 4, 225-231

use parametersmod, only : sp

implicit none

!arguments
real(sp), dimension(:), intent(in)  :: m      ! vector of mean values at super-time step (e.g., monthly), minimum three values
integer,  dimension(:), intent(in)  :: dmonth ! vector of number of intervals for the time step (e.g., days per month)
real(sp), dimension(2), intent(in)  :: bcond  ! boundary conditions for the result vector (1=left side, 2=right side)
real(sp), dimension(:), intent(out) :: r      ! result vector of values at chosen time step

!parameters
real(sp), parameter :: ot = 1. / 3

!local variables
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

!----------

n  = size(m)
ni = size(r)

bc = bcond

!initialize the result vector
i = 1
do a = 1,n
  j = i
  do b = 1,dmonth(a)
    r(i) = m(a)
    g(i) = j
    i = i + 1
  end do
end do

!iteratively smooth and correct the result to preserve the mean

!iteration loop
do i = 1,ni

  do j = 2,ni-1
    r(j) = ot * (r(j-1) + r(j) + r(j+1))   !Eqn. 1
  end do

  r(1)  = ot * (bc(1)   + r(1)  +  r(2))   !Eqns. 2
  r(ni) = ot * (r(ni-1) + r(ni) + bc(2))

  j = 1
  do k = 1,n                               !calculate one correction factor per super-timestep

    a = g(j)                               !index of the first timestep value of the super-timestep
    b = g(j) + dmonth(k) - 1               !index of the last timestep value of the super-timestep

    ck = sum(m(k) - r(a:b)) / ni           !Eqn. 4

    do l = 1,dmonth(k)                     !apply the correction to all timestep values in the super-timestep
      r(j) = r(j) + ck
      j = j + 1
    end do

    !correction for circular conditions when using climatology (do not use for transient simulations)
    bc(1) = r(ni)
    bc(2) = r(1)

  end do
end do

end subroutine rmsmooth

!-----------------------------------------------------------------------

real(sp) function roundto(val,precision)

! round a value to the given precision

implicit none

real(sp), intent(in) :: val        ! the input value
integer,  intent(in) :: precision  ! the precision

real(sp) :: scale

!----

scale = 10.**precision

roundto = real(nint(val * scale)) / scale

end function roundto

!---------------------------------------------------------

subroutine temp_sd(pday,dm)

! calculate the standard deviation of tmin and tmax depending on the input value and precipitation state

implicit none

logical,          intent(in)    :: pday    ! precipitation status (mm/day)
type(daymetvars), intent(inout) :: dm

integer :: i

!--------------------

do i=1,4

  if (i == 4 .or. dm%tmin_mn <= tmin_sd_breaks(i)) then

    if (pday) then
    
      dm%tmin_sd = sum(tmin_sd_w(i,:) * (dm%tmin_mn**exponents))  !the vector 'exponents' is a clever way of calculating a polynomial expansion

    else

      dm%tmin_sd = sum(tmin_sd_d(i,:) * (dm%tmin_mn**exponents))

    end if

    exit

  end if
end do

do i=1,4

  if (i == 4 .or. dm%tmax_mn <= tmax_sd_breaks(i)) then

    if (pday) then

      dm%tmax_sd = sum(tmax_sd_w(i, :) * (dm%tmax_mn**exponents))

    else

      dm%tmax_sd = sum(tmax_sd_d(i, :) * (dm%tmax_mn**exponents))

    end if

    exit

  end if
end do

end subroutine temp_sd

!---------------------------------------------------------

end module weathergenmod
