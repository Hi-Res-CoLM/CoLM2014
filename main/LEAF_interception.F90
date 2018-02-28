
#include <define.h>

 subroutine LEAF_interception (deltim,dewmx,chil,sigf,lai,sai,tlsun,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                               ldew,pg_rain,pg_snow,qintr)

!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

  use precision
  use PhysicalConstants, only : tfrz
  implicit none

!-----------------------Arguments---------------------------------------

  real(r8), INTENT(in) :: deltim   ! seconds in a time step [second]
  real(r8), INTENT(in) :: dewmx    ! maximum dew [mm]
  real(r8), INTENT(in) :: chil     ! leaf angle distribution factor
  real(r8), INTENT(in) :: prc_rain ! convective ranfall [mm/s]
  real(r8), INTENT(in) :: prc_snow ! convective snowfall [mm/s]
  real(r8), INTENT(in) :: prl_rain ! large-scale rainfall [mm/s]
  real(r8), INTENT(in) :: prl_snow ! large-scale snowfall [mm/s]
  real(r8), INTENT(in) :: sigf     ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: lai      ! leaf area index [-]
  real(r8), INTENT(in) :: sai      ! stem area index [-]
  real(r8), INTENT(in) :: tlsun    ! sunlit canopy leaf temperature [K]

  real(r8), INTENT(inout) :: ldew  ! depth of water on foliage [mm]
  real(r8), INTENT(out) :: pg_rain ! rainfall onto ground including canopy runoff [kg/(m2 s)]
  real(r8), INTENT(out) :: pg_snow ! snowfall onto ground including canopy runoff [kg/(m2 s)]
  real(r8), INTENT(out) :: qintr   ! interception [kg/(m2 s)]

!-----------------------Local Variables---------------------------------

  real(r8) :: satcap   ! maximum allowed water on canopy [mm]
  real(r8) :: lsai     ! sum of leaf area index and stem area index [-]
  real(r8) :: chiv     ! leaf angle distribution factor
  real(r8) :: ppc      ! convective precipitation in time-step [mm]
  real(r8) :: ppl      ! large-scale precipitation in time-step [mm]
  real(r8) :: p0       ! precipitation in time-step [mm]
  real(r8) :: fpi      ! coefficient of interception
  real(r8) :: pinf     ! interception of precipitation in time step [mm]
  real(r8) :: tti_rain ! direct rain throughfall in time step [mm]
  real(r8) :: tti_snow ! direct snow throughfall in time step [mm]
  real(r8) :: tex_rain ! canopy rain drainage in time step [mm]
  real(r8) :: tex_snow ! canopy snow drainage in time step [mm]
  real(r8) :: vegt     ! sigf*lsai
  real(r8) :: xs       ! proportion of the grid area where the intercepted rainfall 
		       ! plus the preexisting canopy water storage

  real(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
  real(r8) :: thru_rain, thru_snow
  real(r8) :: xsc_rain, xsc_snow
  real pcoefs (2,2)

!-----------------------End Variable List-------------------------------

 if(sigf>=0.001)then

    pcoefs(1,1) = 20.
    pcoefs(1,2) = 0.206e-8
    pcoefs(2,1) = 0.0001 
    pcoefs(2,2) = 0.9999 
    bp = 20. 

    lsai = lai + sai
    vegt = sigf*lsai
    satcap = dewmx*vegt

    p0 = (prc_rain + prc_snow + prl_rain + prl_snow)*deltim
    ppc = (prc_rain+prc_snow)*deltim
    ppl = (prl_rain+prl_snow)*deltim

    w = ldew+p0

    if(tlsun > tfrz)then
       xsc_rain = max(0., ldew-satcap)
       xsc_snow = 0.
    else
       xsc_rain = 0.
       xsc_snow = max(0., ldew-satcap)
    endif

    ldew = ldew - (xsc_rain + xsc_snow)

    ap = pcoefs(2,1)
    cp = pcoefs(2,2)

    if(p0>1.e-8)then
       ap = ppc/p0 * pcoefs(1,1) + ppl/p0 * pcoefs(2,1)
       cp = ppc/p0 * pcoefs(1,2) + ppl/p0 * pcoefs(2,2)

!----------------------------------------------------------------------
!      proportional saturated area (xs) and leaf drainage(tex)
!-----------------------------------------------------------------------

       chiv = chil
       if ( abs(chiv) .le. 0.01 ) chiv = 0.01
       aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
       bb = 0.877 * ( 1. - 2. * aa )
       exrain = aa + bb

       ! coefficient of interception
       ! set fraction of potential interception to max 0.25 (Lawrence et al. 2007)
       alpha = 0.25
       fpi = alpha * ( 1.-exp(-exrain*lsai) ) * sigf
       tti_rain = (prc_rain+prl_rain)*deltim * ( 1.-fpi )
       tti_snow = (prc_snow+prl_snow)*deltim * ( 1.-fpi )

       xs = 1.
       if(p0*fpi>1.e-9)then
          arg = (satcap-ldew)/(p0*fpi*ap) - cp/ap
          if(arg>1.e-9)then
             xs = -1./bp * log( arg )
             xs = min( xs, 1. )
             xs = max( xs, 0. )
          endif
       endif

       ! assume no fall down of the intercepted snowfall in a time step
       tex_rain = (prc_rain+prl_rain)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) &
                - (satcap-ldew) * xs
       tex_rain = max( tex_rain, 0. )
       tex_snow = 0.

#if(defined CLMDEBUG)
       if(tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10) then
         write(6,*) 'tex_ + tti_ > p0 in interception code : '
       endif
#endif

    else
       ! all intercepted by canopy leves for very small precipitation
       tti_rain = 0.
       tti_snow = 0.
       tex_rain = 0.
       tex_snow = 0.
    endif

!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------

    thru_rain = tti_rain + tex_rain
    thru_snow = tti_snow + tex_snow
    pinf = p0 - (thru_rain + thru_snow)
    ldew = ldew + pinf

    pg_rain = (xsc_rain + thru_rain) / deltim
    pg_snow = (xsc_snow + thru_snow) / deltim
    qintr   = pinf / deltim

#if(defined CLMDEBUG)
    w = w - ldew - (pg_rain+pg_snow)*deltim
    if(abs(w)>1.e-6)then
       write(6,*) 'something wrong in interception code : '
       write(6,*) w, ldew, (pg_rain+pg_snow)*deltim, satcap
       call abort
    endif
#endif

 else

    ldew=0.
    pg_rain = prc_rain + prl_rain
    pg_snow = prc_snow + prl_snow
    qintr   = 0.

 endif

 end subroutine LEAF_interception
