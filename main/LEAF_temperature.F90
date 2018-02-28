#include <define.h>

MODULE LEAF_temperature

!-----------------------------------------------------------------------
 use precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: leaftemone
  public :: leaftemtwo


! PRIVATE MEMBER FUNCTIONS:
  private :: dewfraction


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  subroutine  leaftemone (deltim ,csoilc ,dewmx  ,htvp   ,lai    ,&
              sai     ,displa    ,sqrtdi ,z0m    ,effcon ,vmax25 ,&
              slti    ,hlti      ,shti   ,hhti   ,trda   ,trdm   ,&
              trop    ,gradm     ,binter ,extkn  ,extkb  ,extkd  ,&
              hu      ,ht        ,hq     ,us     ,vs     ,thm    ,&
              th      ,thv       ,qm     ,psrf   ,rhoair ,par    ,&
              sabv    ,frl       ,thermk ,rstfac ,po2m   ,pco2m  ,&
              sigf    ,etrc      ,tg     ,qg     ,dqgdT  ,emg    ,&
              tl      ,ldew      ,taux   ,tauy   ,fseng  ,fevpg  ,&
              cgrnd   ,cgrndl    ,cgrnds ,tref   ,qref   ,rst    ,&
              assim   ,respc     ,fsenl  ,fevpl  ,etr    ,dlrad  ,&
              ulrad   ,z0ma      ,zol    ,rib    ,ustar  ,qstar  ,&
              tstar   ,fm        ,fh     ,fq) 
 
!=======================================================================
! Original author : Yongjiu Dai, August 15, 2001
!
! Foliage energy conservation is given by foliage energy budget equation
!                      Rnet - Hf - LEf = 0
! The equation is solved by Newton-Raphson iteration, in which this iteration
! includes the calculation of the photosynthesis and stomatal resistance, and the
! integration of turbulent flux profiles. The sensible and latent heat
! transfer between foliage and atmosphere and ground is linked by the equations:
!                      Ha = Hf + Hg and Ea = Ef + Eg
!
! ________________
! REVISION HISTORY:
! 07/09/2014, Hua Yuan: imbalanced energy due to T/q adjustment is
!                       allocated to sensible heat flux.
!
!=======================================================================

  use precision
  use PhysicalConstants, only : vonkar, grav, hvap, cpair, stefnc
  use FRICTION_VELOCITY
  use ASSIM_STOMATA_conductance
  implicit none
 
!-----------------------Arguments---------------------------------------

  real(r8), INTENT(in) :: &
        deltim,     &! seconds in a time step [second]
        csoilc,     &! drag coefficient for soil under canopy [-]
        dewmx,      &! maximum dew
        htvp         ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  real(r8), INTENT(in) :: &
        lai,        &! adjusted leaf area index for seasonal variation [-]
        sai,        &! stem area index  [-]
        displa,     &! displacement height [m]
        sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]
        z0m,        &! roughness length, momentum [m]

        effcon,     &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,     &! maximum carboxylation rate at 25 C at canopy top
                     ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,       &! slope of high temperature inhibition function     (s1)
        hhti,       &! 1/2 point of high temperature inhibition function (s2)
        slti,       &! slope of low temperature inhibition function      (s3)
        hlti,       &! 1/2 point of low temperature inhibition function  (s4)
        trda,       &! temperature coefficient in gs-a model             (s5)
        trdm,       &! temperature coefficient in gs-a model             (s6)
        trop,       &! temperature coefficient in gs-a model         (273+25)
        gradm,      &! conductance-photosynthesis slope parameter
        binter,     &! conductance-photosynthesis intercept
        extkn        ! coefficient of leaf nitrogen allocation

! input variables
  real(r8), INTENT(in) :: &
        hu,         &! observational height of wind [m]
        ht,         &! observational height of temperature [m]
        hq,         &! observational height of humidity [m]
        us,         &! wind component in eastward direction [m/s]
        vs,         &! wind component in northward direction [m/s]
        thm,        &! intermediate variable (tm+0.0098*ht)
        th,         &! potential temperature (kelvin)
        thv,        &! virtual potential temperature (kelvin)
        qm,         &! specific humidity at reference height [kg/kg]
        psrf,       &! pressure at reference height [pa]
        rhoair,     &! density air [kg/m**3]

        par,        &! par absorbed per unit lai [w/m**2]
        sabv,       &! solar radiation absorbed by vegetation [W/m2]
        frl,        &! atmospheric infrared (longwave) radiation [W/m2]

        extkb,      &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,      &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,     &! canopy gap fraction for tir radiation
        rstfac,     &! factor of soil water stress to plant physiologocal processes

        po2m,       &! atmospheric partial pressure  o2 (pa)
        pco2m,      &! atmospheric partial pressure co2 (pa)

        sigf,       &! fraction of veg cover, excluding snow-covered veg [-]
        etrc,       &! maximum possible transpiration rate (mm/s)
        tg,         &! ground surface temperature [K]
        qg,         &! specific humidity at ground surface [kg/kg]
        dqgdT,      &! temperature derivative of "qg"
        emg          ! vegetation emissivity

  real(r8), INTENT(inout) :: &
        tl,         &! leaf temperature [K]
        ldew,       &! depth of water on foliage [mm]
        taux,       &! wind stress: E-W [kg/m/s**2]
        tauy,       &! wind stress: N-S [kg/m/s**2]
        fseng,      &! sensible heat flux from ground [W/m2]
        fevpg,      &! evaporation heat flux from ground [mm/s]
        cgrnd,      &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,     &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,     &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,       &! 2 m height air temperature (kelvin)
        qref         ! 2 m height air specific humidity

  real(r8), INTENT(out) :: &
        rst,        &! stomatal resistance
        assim,      &! rate of assimilation
        respc,      &! rate of respiration
        fsenl,      &! sensible heat from leaves [W/m2]
        fevpl,      &! evaporation+transpiration from leaves [mm/s]
        etr,        &! transpiration rate [mm/s]
        dlrad,      &! downward longwave radiation blow the canopy [W/m2]
        ulrad,      &! upward longwave radiation above the canopy [W/m2]

        z0ma,       &! effective roughness [m]
        zol,        &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,        &! bulk Richardson number in surface layer
        ustar,      &! friction velocity [m/s]
        tstar,      &! temperature scaling parameter
        qstar,      &! moisture scaling parameter
        fm,         &! integral of profile function for momentum
        fh,         &! integral of profile function for heat
        fq           ! integral of profile function for moisture

!-----------------------Local Variables---------------------------------
! assign iteration parameters
   integer, parameter :: itmax  = 40   ! maximum number of iteration
   integer, parameter :: itmin  = 6    ! minimum number of iteration
   real(r8),parameter :: delmax = 3.0  ! maximum change in leaf temperature [K]
   real(r8),parameter :: dtmin  = 0.01 ! max limit for temperature convergence [K]
   real(r8),parameter :: dlemin = 0.1  ! max limit for energy flux convergence [w/m2]

   real(r8) dtl(0:itmax+1)     ! difference of tl between two iterative step

   real(r8) :: &
        zldis,      &! reference height "minus" zero displacement heght [m]
        zii,        &! convective boundary layer height [m]
        z0mv,       &! roughness length, momentum [m]
        z0hv,       &! roughness length, sensible heat [m]
        z0qv,       &! roughness length, latent heat [m]
        zeta,       &! dimensionless height used in Monin-Obukhov theory
        beta,       &! coefficient of conective velocity [-]
        wc,         &! convective velocity [m/s]
        wc2,        &! wc**2
        dth,        &! diff of virtual temp. between ref. height and surface 
        dthv,       &! diff of vir. poten. temp. between ref. height and surface
        dqh,        &! diff of humidity between ref. height and surface
        obu,        &! monin-obukhov length (m)
        um,         &! wind speed including the stablity effect [m/s]
        ur,         &! wind speed at reference height [m/s]
        uaf,        &! velocity of air within foliage [m/s]
        fh2m,       &! relation for temperature at 2m
        fq2m,       &! relation for specific humidity at 2m
        fm10m,      &! integral of profile function for momentum at 10m
        thvstar,    &! virtual potential temperature scaling parameter
        taf,        &! air temperature within canopy space [K]
        qaf,        &! humidity of canopy air [kg/kg]
        eah,        &! canopy air vapor pressure (pa)
        pco2g,      &! co2 pressure (pa) at ground surface (pa)
        pco2a,      &! canopy air co2 pressure (pa)

        fdry,       &! fraction of foliage that is green and dry [-]
        fwet,       &! fraction of foliage covered by water [-]
        cf,         &! heat transfer coefficient from leaves [-]
        rb,         &! leaf boundary layer resistance [s/m]
        rbone,      &! canopy bulk boundary layer resistance 
        rd,         &! aerodynamical resistance between ground and canopy air
        ram,        &! aerodynamical resistance [s/m]
        rah,        &! thermal resistance [s/m]
        raw,        &! moisture resistance [s/m]
        clai,       &! canopy heat capacity [Jm-2K-1]
        cah,        &! heat conductance for air [m/s]
        cgh,        &! heat conductance for ground [m/s]
        cfh,        &! heat conductance for leaf [m/s]
        caw,        &! latent heat conductance for air [m/s]
        cgw,        &! latent heat conductance for ground [m/s]
        cfw,        &! latent heat conductance for leaf [m/s]
        wtshi,      &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,      &! latent heat resistance for air, grd and leaf [-]
        wta0,       &! normalized heat conductance for air [-]
        wtg0,       &! normalized heat conductance for ground [-]
        wtl0,       &! normalized heat conductance for air and leaf [-]
        wtaq0,      &! normalized latent heat conductance for air [-]
        wtgq0,      &! normalized heat conductance for ground [-]
        wtlq0,      &! normalized latent heat cond. for air and leaf [-]

        ei,         &! vapor pressure on leaf surface [pa]
        deidT,      &! derivative of "ei" on "tl" [pa/K]
        qsatl,      &! leaf specific humidity [kg/kg]
        qsatldT,    &! derivative of "qsatl" on "tlef"

        del,        &! absolute change in leaf temp in current iteration [K]
        del2,       &! change in leaf temperature in previous iteration [K]
        dele,       &! change in heat fluxes from leaf [K]
        dele2,      &! change in heat fluxes from leaf [K]
        det,        &! maximum leaf temp. change in two consecutive iter [K]
 
        obuold,     &! monin-obukhov length from previous iteration
        tlbef,      &! leaf temperature from previous iteration [K]
        ecidif,     &! excess energies [W/m2]
        err,        &! balance error

        rs,         &! leaf stomatal resistance [s/m]
        rsoil,      &! soil respiration
        gah2o,      &! conductance between canopy and atmosphere
        gdh2o,      &! conductance between canopy and ground
        tprcor       ! tf*psur*100./1.013e5

   integer it, nmozsgn 

   real(r8) delta, fac
   real(r8) evplwet, evplwet_dtl, etr_dtl, elwmax, elwdif
   real(r8) irab, dirab_dtl, fsenl_dtl, fevpl_dtl  
   real(r8) w, csoilcn, z0mg, cint(3)
   real(r8) fevpl_bef, fevpl_noadj, dtl_noadj, errt, erre

!-----------------------End Variable List-------------------------------

! initialization of errors and  iteration parameters
       it     = 1    ! counter for leaf temperature iteration
       del    = 0.0  ! change in leaf temperature from previous iteration
       dele   = 0.0  ! latent head flux from leaf for previous iteration

       dtl(0) = 0.
       fevpl_bef = 0.

!-----------------------------------------------------------------------
! scaling-up coefficients from leaf to canopy
!-----------------------------------------------------------------------

! /05/2014/
!      cint(1) = (1.-exp(-extkn*lai))/extkn
       cint(1) = (1.-exp(-0.110*lai))/0.110
       cint(2) = (1.-exp(-extkd*lai))/extkd
       cint(3) = lai

!-----------------------------------------------------------------------
! get fraction of wet and dry canopy surface (fwet & fdry)
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

       !clai = 4.2 * 1000. * 0.2
       clai = 0.0

       call dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)

       call qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

       nmozsgn = 0    ! number of times moz changes sign
       obuold = 0.    ! monin-obukhov length from previous iteration
       zii = 1000.    ! m  (pbl height)
       beta = 1.      ! -  (in computing W_*)

       z0mv = z0m; z0hv = z0m; z0qv = z0m

       taf = 0.5 * (tg + thm)
       qaf = 0.5 * (qm + qg)

       pco2a = pco2m
       tprcor = 44.6*273.16*psrf/1.013e5
       rsoil = 0.  !respiration (mol m-2 s-1)
!      rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13)))
!      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
!      rsoil = 5.22 * 1.e-6
       rsoil = 0.22 * 1.e-6

       ur = max(0.1, sqrt(us*us+vs*vs))    ! limit set to 0.1
       dth = thm - taf
       dqh = qm - qaf
       dthv = dth*(1.+0.61*qm) + 0.61*th*dqh
       zldis = hu - displa

       if(zldis <= 0.0) then
          write(6,*) 'the obs height of u less than the zero displacement heght'
          call abort
       endif

       call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu)

! ======================================================================
!     BEGIN stability iteration 
! ======================================================================

      DO WHILE (it .le. itmax) 

         tlbef = tl

         del2 = del
         dele2 = dele
 
!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration
        call moninobuk(hu,ht,hq,displa,z0mv,z0hv,z0qv,obu,um,&
                       ustar,fh2m,fq2m,fm10m,fm,fh,fq)
 
! Aerodynamic resistance
        ram = 1./(ustar*ustar/um)
        rah = 1./(vonkar/fh*ustar) 
        raw = 1./(vonkar/fq*ustar) 
 
! Bulk boundary layer resistance of leaves
        uaf = ustar
        cf = 0.01*sqrtdi/sqrt(uaf)
        rb = 1/(cf*uaf)

      ! rd = 1./(csoilc*uaf)                 ! BATS legacy
      ! w = exp(-0.5*(lai+sai))              ! Dickinson's modification :
      ! csoilc = ( 1.-w + w*um/uaf)/rah      ! "rah" here is the resistance over
      ! rd = 1./(csoilc*uaf)                 ! bare ground fraction

! modified by Xubin Zeng's suggestion at 08-07-2002
        z0mg = 0.01
        w = exp(-(lai+sai))
        csoilcn = (vonkar/(0.13*(z0mg*uaf/1.5e-5)**0.45))*w + csoilc*(1.-w)
        rd = 1./(csoilcn*uaf)

!-----------------------------------------------------------------------
! stomatal resistances 
!-----------------------------------------------------------------------

        if(lai .gt. 0.001) then
           rbone = rb / lai
           eah = qaf * psrf / ( 0.622 + 0.378 * qaf )    ! pa

           CALL stomata (vmax25 ,effcon ,slti   ,hlti ,&
                shti    ,hhti   ,trda   ,trdm   ,trop ,&
                gradm   ,binter ,thm    ,psrf   ,po2m ,&
                pco2m   ,pco2a  ,eah    ,ei     ,tl   ,&
                par     ,rbone  ,raw    ,rstfac ,cint ,&
                assim   ,respc  ,rs     )
        else
           rs = 2.e4; assim = 0.; respc = 0.
        endif

! above stomatal resistances are for the canopy, the stomatal rsistances 
! and the "rb" in the following calculations are the average for single leaf. thus,
        rs = rs * lai

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

        delta = 0.0
        if(qsatl-qaf .gt. 0.) delta = 1.0
 
        cah = sigf / rah
        cgh = sigf / rd
        cfh = sigf * (lai + sai) / rb

        caw = sigf / raw
        cgw = sigf / rd
        cfw = sigf * ( (1.-delta*(1.-fwet))*(lai+sai)/rb + (1.-fwet)*delta*lai/(rb+rs) )

        wtshi = 1. / ( cah + cgh + cfh )
        wtsqi = 1. / ( caw + cgw + cfw )

        wta0 = cah * wtshi
        wtg0 = cgh * wtshi
        wtl0 = cfh * wtshi

        wtaq0 = caw * wtsqi
        wtgq0 = cgw * wtsqi
        wtlq0 = cfw * wtsqi
 
!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored 
! which cannot be determined analtically
        fac = sigf * (1. - thermk)

! longwave absorption and their derivatives 
        irab = (frl - 2. * stefnc * tl**4 + emg*stefnc*tg**4 ) * fac 
        dirab_dtl = - 8. * stefnc * tl**3                      * fac

! sensible heat fluxes and their derivatives
        fsenl = rhoair * cpair * cfh * ( (wta0 + wtg0)*tl - wta0*thm - wtg0*tg )
        fsenl_dtl = rhoair * cpair * cfh * (wta0 + wtg0)

! latent heat fluxes and their derivatives
        etr = sigf * rhoair * (1.-fwet) * delta * lai / (rb + rs) &
            * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg )
        etr_dtl = sigf * rhoair * (1.-fwet) * delta * lai / (rb + rs) &
            * (wtaq0 + wtgq0)*qsatlDT 
     
        if(etr.ge.sigf*etrc)then
           etr = sigf*etrc
           etr_dtl = 0.
        endif

        evplwet = sigf * rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb &
                * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg )
        evplwet_dtl = sigf * rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb &
                    * (wtaq0 + wtgq0)*qsatlDT 
        if(evplwet.ge.ldew/deltim)then
           evplwet = ldew/deltim
           evplwet_dtl = 0.
        endif
 
        fevpl = etr + evplwet
        fevpl_dtl = etr_dtl + evplwet_dtl

        erre = 0.
        fevpl_noadj = fevpl
        if ( fevpl*fevpl_bef < 0. ) then 
           erre  = -0.9*fevpl
           fevpl =  0.1*fevpl
        endif

!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------

        dtl(it) = (sabv + irab - fsenl - hvap*fevpl) &
            / ((lai+sai)*clai/deltim - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl)
        dtl_noadj = dtl(it)

        ! check magnitude of change in leaf temperature limit to maximum allowed value

        if(it .le. itmax) then

        ! put brakes on large temperature excursions
          if(abs(dtl(it)).gt.delmax)then
              dtl(it) = delmax*dtl(it)/abs(dtl(it))
          endif

          if((it.ge.2) .and. (dtl(it-1)*dtl(it).le.0.))then
              dtl(it) = 0.5*(dtl(it-1) + dtl(it))
          endif

        endif

        tl = tlbef + dtl(it)

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------

        del  = sqrt( dtl(it)*dtl(it) )
        dele = dtl(it) * dtl(it) * ( dirab_dtl**2 + fsenl_dtl**2 + hvap*fevpl_dtl**2 ) 
        dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
        call qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT)

! update vegetation/ground surface temperature, canopy air temperature, 
! canopy air humidity
        taf = wta0*thm + wtg0*tg + wtl0*tl 

        qaf = wtaq0*qm + wtgq0*qg + wtlq0*qsatl

! update co2 partial pressure within canopy air
        gah2o = 1.0/raw * tprcor/thm                     ! mol m-2 s-1
        gdh2o = 1.0/rd  * tprcor/thm                     ! mol m-2 s-1
        pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * (assim - respc - rsoil)

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

        dth = thm - taf       
        dqh = qm - qaf

        tstar = vonkar/fh*dth
        qstar = vonkar/fq*dqh

        thvstar = tstar + 0.61*th*qstar
        zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)
        if(zeta .ge. 0.)then                             !stable
           zeta = min(2.,max(zeta,1.e-6))
        else                                             !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
        endif
        obu = zldis/zeta

        if(zeta .ge. 0.)then
          um = max(ur,.1)
        else
          wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
         wc2 = beta*beta*(wc*wc)
          um = sqrt(ur*ur+wc2)
        endif

        if(obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
        if(nmozsgn .ge. 4) obu = zldis/(-0.01)
        obuold = obu
 
!-----------------------------------------------------------------------
! Test for convergence
!-----------------------------------------------------------------------

      it = it+1

      if(it .gt. itmin) then
         fevpl_bef = fevpl
         det = max(del,del2)
         del = max(dele,dele2)
         if(det .lt. dtmin .AND. del .lt. dlemin) exit 
      endif
 
      END DO 

! ======================================================================
!     END stability iteration 
! ======================================================================

      z0ma = z0mv
      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! canopy fluxes and total assimilation amd respiration

      if(lai .gt. 0.001) then
         rst = rs / lai
      else
        rst = 2.0e4
        assim = 0.
        respc = 0.
      endif
      respc = respc + rsoil

! canopy fluxes and total assimilation amd respiration
      fsenl = fsenl + fsenl_dtl*dtl(it-1) &
              ! add the imbalanced energy below due to T adjustment to sensibel heat
              + (dtl_noadj-dtl(it-1)) * ((lai+sai)*clai/deltim - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl) &
              ! add the imbalanced energy below due to q adjustment to sensibel heat
              + hvap*erre

      etr     = etr     +     etr_dtl*dtl(it-1)
      evplwet = evplwet + evplwet_dtl*dtl(it-1)
      fevpl   = fevpl_noadj
      fevpl   = fevpl   +   fevpl_dtl*dtl(it-1)

      elwmax  = ldew/deltim
      elwdif  = max(0., evplwet-elwmax)
      evplwet = min(evplwet, elwmax)

      fevpl = fevpl - elwdif
      fsenl = fsenl + hvap*elwdif 

      ! wind stresses 
      taux = taux - sigf*rhoair*us/ram
      tauy = tauy - sigf*rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

      fseng = fseng + cpair*rhoair*cgh*(tg-taf)
      fevpg = fevpg +    rhoair*cgw*(qg-qaf)

!-----------------------------------------------------------------------
! downward (upward) longwave radiation below (above) the canopy
!-----------------------------------------------------------------------

      dlrad = sigf * thermk * frl &
              + stefnc * fac * tlbef**3 * (tlbef + 4.*dtl(it-1)) 
      ulrad = stefnc * ( fac * tlbef**3 * (tlbef + 4.*dtl(it-1)) &
              + sigf*thermk*emg*tg**4 )  

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------

      cgrnds = cgrnds + cpair*rhoair*cgh*(1.-wtg0)
      cgrndl = cgrndl + rhoair*cgw*(1.-wtgq0)*dqgdT
      cgrnd  = cgrnds + cgrndl*htvp

!-----------------------------------------------------------------------
! balance check
! (the computational error was created by the assumed 'dtl' in line 406-408) 
!-----------------------------------------------------------------------

      err = sabv + irab + dirab_dtl*dtl(it-1) - fsenl - hvap*fevpl

#if(defined CLMDEBUG)
      if(abs(err) .gt. .2) &
      write(6,*) 'energy imbalance in leaftemone.F90',it-1,err,sabv,irab,fsenl,hvap*fevpl
#endif

!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

      ldew = max(0., ldew-evplwet*deltim)

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------

      tref = tref + sigf*(thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)) 
      qref = qref + sigf*( qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar))


  end subroutine leaftemone
!----------------------------------------------------------------------         



  subroutine leaftemtwo (deltim ,csoilc  ,dewmx  ,htvp   ,lai    ,&
             sai     ,displa    ,sqrtdi  ,z0m    ,effcon ,vmax25 ,&
             slti    ,hlti      ,shti    ,hhti   ,trda   ,trdm   ,&
             trop    ,gradm     ,binter  ,extkn  ,extkb  ,extkd  ,&
             hu      ,ht        ,hq      ,us     ,vs     ,thm    ,&
             th      ,thv       ,qm      ,psrf   ,rhoair ,parsun ,&
             parsha  ,sabvsun   ,sabvsha ,frl    ,fsun   ,thermk ,&
             rstfac  ,po2m      ,pco2m   ,sigf   ,etrc   ,tg     ,&
             qg      ,dqgdT     ,emg     ,tlsun  ,tlsha  ,ldew   ,&
             taux    ,tauy      ,fseng   ,fevpg  ,cgrnd  ,cgrndl ,&
             cgrnds  ,tref      ,qref    ,rst    ,assim  ,respc  ,&
             fsenl   ,fevpl     ,etr     ,dlrad  ,ulrad  ,z0ma   ,&
             zol     ,rib       ,ustar   ,qstar  ,tstar  ,&
             fm      ,fh        ,fq      )

!=======================================================================
! Original author : Yongjiu Dai, August 15, 2001
!
! Foliage energy conservation is given by foliage energy budget equation
!                   sunlit:  [Rnet - Hf - LEf] = 0
!                   shaded:  [Rnet - Hf - LEf] = 0
! The equation is solved by Newton-Raphson iteration, in which this iteration
! includes the calculation of the photosynthesis and stomatal resistance, and the
! integration of turbulent flux profiles. The sensible and latent heat
! transfer between foliage and atmosphere and ground is linked by the equations:
!                   Ha = [Hf]sunlit + [Hf]shaded + Hg 
!                   Ea = [Ef]sunlit + [Ef]shaded + Eg
!
! ________________
! REVISION HISTORY:
! 07/09/2014, Hua Yuan: imbalanced energy due to T/q adjustment is
!                       allocated to sensible heat flux.
!
!=======================================================================

  use precision
  use PhysicalConstants, only : vonkar, grav, hvap, cpair, stefnc
  use FRICTION_VELOCITY
  use ASSIM_STOMATA_conductance
  implicit none
 
!-----------------------Arguments---------------------------------------

  real(r8), INTENT(in) :: &
        deltim,     &! seconds in a time step [second]
        csoilc,     &! drag coefficient for soil under canopy [-]
        dewmx,      &! maximum dew
        htvp         ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  real(r8), INTENT(in) :: &
        lai,        &! adjusted leaf area index for seasonal variation [-]
        sai,        &! stem area index  [-]
        displa,     &! displacement height [m]
        sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]
        z0m,        &! roughness length, momentum [m]

        effcon,     &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,     &! maximum carboxylation rate at 25 C at canopy top
                     ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,       &! slope of high temperature inhibition function     (s1)
        hhti,       &! 1/2 point of high temperature inhibition function (s2)
        slti,       &! slope of low temperature inhibition function      (s3)
        hlti,       &! 1/2 point of low temperature inhibition function  (s4)
        trda,       &! temperature coefficient in gs-a model             (s5)
        trdm,       &! temperature coefficient in gs-a model             (s6)
        trop,       &! temperature coefficient in gs-a model             (273+25)
        gradm,      &! conductance-photosynthesis slope parameter
        binter,     &! conductance-photosynthesis intercept
        extkn        ! coefficient of leaf nitrogen allocation

! input variables
  real(r8), INTENT(in) :: &
        hu,         &! observational height of wind [m]
        ht,         &! observational height of temperature [m]
        hq,         &! observational height of humidity [m]
        us,         &! wind component in eastward direction [m/s]
        vs,         &! wind component in northward direction [m/s]
        thm,        &! intermediate variable (tm+0.0098*ht)
        th,         &! potential temperature (kelvin)
        thv,        &! virtual potential temperature (kelvin)
        qm,         &! specific humidity at reference height [kg/kg]
        psrf,       &! pressure at reference height [pa]
        rhoair,     &! density air [kg/m**3]

        parsun,     &! par absorbed per unit sunlit lai [w/m**2]
        parsha,     &! par absorbed per unit shaded lai [w/m**2]
        sabvsun,    &! solar radiation absorbed by vegetation [W/m2]
        sabvsha,    &! solar radiation absorbed by vegetation [W/m2]
        frl,        &! atmospheric infrared (longwave) radiation [W/m2]

        fsun,       &! sunlit fraction of canopy
        extkb,      &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,      &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,     &! canopy gap fraction for tir radiation
        rstfac,     &! factor of soil water stress to transpiration

        po2m,       &! atmospheric partial pressure  o2 (pa)
        pco2m,      &! atmospheric partial pressure co2 (pa)

        sigf,       &! fraction of veg cover, excluding snow-covered veg [-]
        etrc,       &! maximum possible transpiration rate (mm/s)
        tg,         &! ground surface temperature [K]
        qg,         &! specific humidity at ground surface [kg/kg]
        dqgdT,      &! temperature derivative of "qg"
        emg          ! vegetation emissivity

  real(r8), INTENT(inout) :: &
        tlsun,      &! sunlit leaf temperature [K]
        tlsha,      &! shaded leaf temperature [K]
        ldew,       &! depth of water on foliage [mm]
        taux,       &! wind stress: E-W [kg/m/s**2]
        tauy,       &! wind stress: N-S [kg/m/s**2]
        fseng,      &! sensible heat flux from ground [W/m2]
        fevpg,      &! evaporation heat flux from ground [mm/s]
        cgrnd,      &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,     &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,     &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,       &! 2 m height air temperature (kelvin)
        qref         ! 2 m height air specific humidity

  real(r8), INTENT(out) :: &
        rst,        &! total stomatal resistance (s m-1)
        assim,      &! total assimilation (mol m-2 s-1)
        respc,      &! total respiration (mol m-2 s-1)
        fsenl,      &! sensible heat from leaves [W/m2]
        fevpl,      &! evaporation+transpiration from leaves [mm/s]
        etr,        &! transpiration rate [mm/s]
        dlrad,      &! downward longwave radiation blow the canopy [W/m2]
        ulrad,      &! upward longwave radiation above the canopy [W/m2]

        z0ma,       &! effective roughness [m]
        zol,        &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,        &! bulk Richardson number in surface layer
        ustar,      &! friction velocity [m/s]
        tstar,      &! temperature scaling parameter
        qstar,      &! moisture scaling parameter
        fm,         &! integral of profile function for momentum
        fh,         &! integral of profile function for heat
        fq           ! integral of profile function for moisture

!-----------------------Local Variables---------------------------------

! assign iteration parameters
   integer,  parameter :: itmax  = 40   ! maximum number of iteration
   integer,  parameter :: itmin  = 6    ! minimum number of iteration
   real(r8), parameter :: delmax = 1.0  ! maximum change in leaf temperature [K]
   real(r8), parameter :: dtmin  = 0.01 ! max limit for temperature convergence [K]
   real(r8), parameter :: dlemin = 0.1  ! max limit for energy flux convergence [w/m2]

   real(r8) dtlsun(0:itmax+1)     ! difference of tlsun between two iterative step
   real(r8) dtlsha(0:itmax+1)     ! difference of tlsha between two iterative step

   real(r8) :: &
        zldis,      &! reference height "minus" zero displacement heght [m]
        zii,        &! convective boundary layer height [m]
        z0mv,       &! roughness length, momentum [m]
        z0hv,       &! roughness length, sensible heat [m]
        z0qv,       &! roughness length, latent heat [m]
        zeta,       &! dimensionless height used in Monin-Obukhov theory
        beta,       &! coefficient of conective velocity [-]
        wc,         &! convective velocity [m/s]
        wc2,        &! wc**2
        dth,        &! diff of virtual temp. between ref. height and surface 
        dthv,       &! diff of vir. poten. temp. between ref. height and surface
        dqh,        &! diff of humidity between ref. height and surface
        obu,        &! monin-obukhov length (m)
        um,         &! wind speed including the stablity effect [m/s]
        ur,         &! wind speed at reference height [m/s]
        uaf,        &! velocity of air within foliage [m/s]
        fh2m,       &! relation for temperature at 2m
        fq2m,       &! relation for specific humidity at 2m
        fm10m,      &! integral of profile function for momentum at 10m
        thvstar,    &! virtual potential temperature scaling parameter
        taf,        &! air temperature within canopy space [K]
        qaf,        &! humidity of canopy air [kg/kg]
        eah,        &! canopy air vapor pressure (pa)
        pco2a,      &! canopy air co2 pressure (pa)

        fdry,       &! fraction of foliage that is green and dry [-]
        fwet,       &! fraction of foliage covered by water [-]
        cf,         &! heat transfer coefficient from leaves [-]
        rb,         &! leaf boundary layer resistance [s/m]
        rbsun,      &! bulk boundary layer resistance of sunlit fraction of canopy
        rbsha,      &! bulk boundary layer resistance of shaded fraction of canopy
        rd,         &! aerodynamical resistance between ground and canopy air
        ram,        &! aerodynamical resistance [s/m]
        rah,        &! thermal resistance [s/m]
        raw,        &! moisture resistance [s/m]
        cah,        &! heat conductance for air [m/s]
        cgh,        &! heat conductance for ground [m/s]
        cfsunh,     &! heat conductance for sunlit leaf [m/s]
        cfshah,     &! heat conductance for shaded leaf [m/s]
        caw,        &! latent heat conductance for air [m/s]
        cgw,        &! latent heat conductance for ground [m/s]
        cfsunw,     &! latent heat conductance for sunlit leaf [m/s]
        cfshaw,     &! latent heat conductance for shaded leaf [m/s]
        wtshi,      &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,      &! latent heat resistance for air, grd and leaf [-]
        wta0,       &! normalized heat conductance for air [-]
        wtg0,       &! normalized heat conductance for ground [-]
        wtlsun0,    &! normalized heat conductance for air and sunlit leaf [-]
        wtlsha0,    &! normalized heat conductance for air and shaded leaf [-]
        wtaq0,      &! normalized latent heat conductance for air [-]
        wtgq0,      &! normalized heat conductance for ground [-]
        wtlsunq0,   &! normalized latent heat cond. for air and sunlit leaf [-]
        wtlshaq0,   &! normalized latent heat cond. for air and shaded leaf [-]

        del,        &! absolute change in leaf temp in current iteration [K]
        del2,       &! change in leaf temperature in previous iteration [K]
        dele,       &! change in heat fluxes from leaf [K]
        dele2,      &! change in heat fluxes from leaf [K]
        det,        &! maximum leaf temp. change in two consecutive iter [K]
 
        obuold,     &! monin-obukhov length from previous iteration
        tlsunbef,   &! sunlit leaf temperature from previous iteration [K]
        tlshabef,   &! shaded leaf temperature from previous iteration [K]
        ecidif,     &! excess energies [W/m2]
        err,        &! balance error

        fsha,       &! shaded fraction of canopy
        laisun,     &! sunlit leaf area index, one-sided
        laisha,     &! shaded leaf area index, one-sided
        rssun,      &! sunlit leaf stomatal resistance [s/m]
        rssha,      &! shaded leaf stomatal resistance [s/m]
        assimsun,   &! sunlit leaf assimilation rate [umol co2 /m**2/ s] [+]
        assimsha,   &! shaded leaf assimilation rate [umol co2 /m**2/ s] [+]
        respcsun,   &! sunlit leaf respiration rate [umol co2 /m**2/ s] [+]
        respcsha,   &! shaded leaf respiration rate [umol co2 /m**2/ s] [+]
        rsoil,      &! soil respiration
        gah2o,      &! 
        tprcor       !

   integer  it, nmozsgn
   real(r8) delta1, delta2, fac, fac1, fac2 
   real(r8) cintsun(3), cintsha(3)
   real(r8) etrsun, etrsha, evplwetsun, evplwetsha, evplwet
   real(r8) etrsun_dtlsun, etrsha_dtlsun, evplwetsun_dtlsun, evplwetsha_dtlsun
   real(r8) etrsun_dtlsha, etrsha_dtlsha, evplwetsun_dtlsha, evplwetsha_dtlsha

   real(r8) irabsun, senlsun, evplsun, ftsun
   real(r8) irabsha, senlsha, evplsha, ftsha
   real(r8) dirabsun_dtlsun, senlsun_dtlsun, dirabsun_dtlsha, senlsun_dtlsha  
   real(r8) dirabsha_dtlsun, senlsha_dtlsun, dirabsha_dtlsha, senlsha_dtlsha
   real(r8) evplsun_dtlsun, evplsun_dtlsha, dftsunDTlsun, dftsunDTlsha
   real(r8) evplsha_dtlsun, evplsha_dtlsha, dftshaDTlsun, dftshaDTlsha     

   real(r8) eisun, eisha, deisundT, deishadT
   real(r8) qsatlsun, qsatlsha, qsatlsunDT, qsatlshaDT
   real(r8) z0mg, w, csoilcn
   real(r8) clai, elwmax, elwdif

   real(r8) evplsun_bef, evplsha_bef, dtlsun_noadj, dtlsha_noadj
   real(r8) evplsun_noadj, evplsha_noadj
   real(r8) erresun, erresha
!-----------------------------------------------------------------------
! initialization of errors and  iteration parameters
!-----------------------------------------------------------------------

       it     = 1    ! counter for leaf temperature iteration
       del    = 0.0  ! change in leaf temperature from previous iteration
       dele   = 0.0  ! latent head flux from leaf for previous iteration
 
       dtlsun(0) = 0.
       dtlsha(0) = 0.

       evplsun_bef = 0.
       evplsha_bef = 0.

!-----------------------------------------------------------------------
! leaf area index
!-----------------------------------------------------------------------
! partion visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves
       fsha = 1. - fsun

       laisun = lai*fsun
       laisha = lai*fsha

! scaling-up coefficients from leaf to canopy
! /05/2014/
!      cintsun(1) = (1.-exp(-(extkn+extkb)*lai))/(extkn+extkb)
       cintsun(1) = (1.-exp(-(0.110+extkb)*lai))/(0.110+extkb)
       cintsun(2) = (1.-exp(-(extkb+extkd)*lai))/(extkb+extkd)
       cintsun(3) = (1.-exp(-extkb*lai))/extkb

! /05/2014/
!      cintsha(1) = (1.-exp(-extkn*lai))/extkn - cintsun(1)
       cintsha(1) = (1.-exp(-0.110*lai))/0.110 - cintsun(1)
       cintsha(2) = (1.-exp(-extkd*lai))/extkd - cintsun(2)
       cintsha(3) = lai - cintsun(3)

!-----------------------------------------------------------------------
! get fraction of wet and dry canopy surface (fwet & fdry)
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

     ! clai = 4.2 * 1000. * 0.2
       clai = 0.0
       call dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)

       call qsadv(tlsun,psrf,eisun,deisunDT,qsatlsun,qsatlsunDT)
       call qsadv(tlsha,psrf,eisha,deishaDT,qsatlsha,qsatlshaDT)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------
       nmozsgn = 0    ! number of times moz changes sign
       obuold = 0.    ! monin-obukhov length from previous iteration
       zii = 1000.    ! m  (pbl height)
       beta = 1.      ! -  (in computing W_*)

       z0mv = z0m; z0hv = z0m; z0qv = z0m

       taf = 0.5 * (tg + thm)
       qaf = 0.5 * (qm + qg)

       pco2a = pco2m
       tprcor = 44.6*273.16*psrf/1.013e5
       rsoil = 0. !soil respiration (mol m-2 s-1)
!      rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13))) 
!      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
!      rsoil = 5.22 * 1.e-6
       rsoil = 0.22 * 1.e-6

       ur = max(0.1, sqrt(us*us+vs*vs))    ! limit set to 0.1
       dth = thm - taf
       dqh = qm - qaf
       dthv = dth*(1.+0.61*qm) + 0.61*th*dqh
       zldis = hu - displa

#if(defined CLMDEBUG)
       if(zldis <= 0.0) then
          write(6,*) 'the obs height of u less than the zero displacement heght'
          call abort
       endif
#endif

       call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu)

! ======================================================================
!      BEGIN stability iteration 
! ======================================================================

       DO WHILE (it .le. itmax) 

          tlsunbef = tlsun
          tlshabef = tlsha

          del2 = del
          dele2 = dele
 
!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration
          call moninobuk(hu,ht,hq,displa,z0mv,z0hv,z0qv,obu,um,&
                        ustar,fh2m,fq2m,fm10m,fm,fh,fq)
 
! Aerodynamic resistance
          ram = 1./(ustar*ustar/um) 
          rah = 1./(vonkar/fh*ustar)
          raw = 1./(vonkar/fq*ustar)
 
! Bulk boundary layer resistance of leaves
          uaf = ustar
          cf = 0.01*sqrtdi/sqrt(uaf)
          rb = 1/(cf*uaf) 

!<->      rd = 1./(csoilc*uaf)          ! legacy from BATS
! modified by Xubin Zeng's suggestion at 08-07-2002
          z0mg = 0.01
          w = exp(-(lai+sai))
          csoilcn = (vonkar/(0.13*(z0mg*uaf/1.5e-5)**0.45))*w + csoilc*(1.-w)
          rd = 1./(csoilcn*uaf)
!-----------------------------------------------------------------------
! stomatal resistances for sunlit and shaded fractions of canopy.
! should do each iteration to account for differences in eah, leaf temperatures.
!-----------------------------------------------------------------------
          if(lai .gt. 0.001) then
             rbsun = rb / laisun
             rbsha = rb / laisha

             eah = qaf * psrf / ( 0.622 + 0.378 * qaf )    ! pa

! Sunlit leaves
             CALL stomata  (vmax25   ,effcon ,slti   ,hlti    ,&
                  shti     ,hhti     ,trda   ,trdm   ,trop    ,&
                  gradm    ,binter   ,thm    ,psrf   ,po2m    ,&
                  pco2m    ,pco2a    ,eah    ,eisun  ,tlsun   ,&
                  parsun   ,rbsun    ,raw    ,rstfac ,cintsun ,&
                  assimsun ,respcsun ,rssun  )

! Shaded leaves
             CALL stomata  (vmax25   ,effcon ,slti   ,hlti    ,&
                  shti     ,hhti     ,trda   ,trdm   ,trop    ,&
                  gradm    ,binter   ,thm    ,psrf   ,po2m    ,&
                  pco2m    ,pco2a    ,eah    ,eisha  ,tlsha   ,&
                  parsha   ,rbsha    ,raw    ,rstfac ,cintsha ,&
                  assimsha ,respcsha ,rssha  )

          else
             rssun = 2.0e4 ; rssha = 2.0e4
             assimsun = 0. ; assimsha = 0.
             respcsun = 0. ; respcsha = 0.
          endif

! above stomatal resistances are for the canopy, the stomatal resistance 
! the "rb" in the following calculations are the averages for single leaf. thus,
          rssun = rssun * laisun
          rssha = rssha * laisha

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------
          delta1 = 0.0
          delta2 = 0.0
          if(qsatlsun-qaf .gt. 0.) delta1 = 1.0
          if(qsatlsha-qaf .gt. 0.) delta2 = 1.0
 
          cah = sigf / rah
          cgh = sigf / rd
          cfsunh = sigf * laisun / rb
          cfshah = sigf * (laisha + sai) / rb

          caw = sigf / raw
          cgw = sigf / rd
          cfsunw = sigf * ( (1.-delta1*(1.-fwet)) * laisun / rb &
                                + (1. - fwet) * delta1 * laisun / (rb + rssun) )
          cfshaw = sigf * ( (1.-delta2*(1.-fwet)) * (laisha + sai) / rb &
                                + (1. - fwet) * delta2 * laisha / (rb + rssha) )

          wtshi = 1. / ( cah + cgh + cfsunh + cfshah )
          wtsqi = 1. / ( caw + cgw + cfsunw + cfshaw )

          wta0 = cah * wtshi
          wtg0 = cgh * wtshi
          wtlsun0 = cfsunh * wtshi
          wtlsha0 = cfshah * wtshi

          wtaq0 = caw * wtsqi
          wtgq0 = cgw * wtsqi
          wtlsunq0 = cfsunw * wtsqi
          wtlshaq0 = cfshaw * wtsqi
 
!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored 
! which cannot be determined analtically
          fac = sigf * (1. - thermk)
          fac1 = fac * fsun
          fac2 = fac * fsha

! longwave absorption and their derivatives 
          irabsun = (frl - 2. * stefnc * tlsun**4 + emg*stefnc*tg**4 ) * fac1 
          dirabsun_dtlsun = - 8.* stefnc * tlsun**3                    * fac1
          dirabsun_dtlsha = 0.

          irabsha = (frl - 2. * stefnc * tlsha**4 + emg*stefnc*tg**4 ) * fac2 
          dirabsha_dtlsha = - 8.* stefnc * tlsha**3                    * fac2 
          dirabsha_dtlsun = 0.

! sensible heat fluxes and their derivatives
          senlsun = rhoair * cpair * cfsunh &
                    * ( (wta0 + wtg0 + wtlsha0)*tlsun &
                    - wta0*thm - wtg0*tg - wtlsha0*tlsha )
          senlsun_dtlsun = rhoair * cpair * cfsunh * (wta0 + wtg0 + wtlsha0)
          senlsun_dtlsha = rhoair * cpair * cfsunh * ( - wtlsha0 )

          senlsha = rhoair * cpair * cfshah &
                    * ( (wta0 + wtg0 + wtlsun0)*tlsha &
                    - wta0*thm - wtg0*tg - wtlsun0*tlsun )
          senlsha_dtlsun = rhoair * cpair * cfshah * ( - wtlsun0 ) 
          senlsha_dtlsha = rhoair * cpair * cfshah * ( wta0 + wtg0 + wtlsun0 )

! latent heat fluxes and their derivatives
          etrsun = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun) &
                   * ( (wtaq0 + wtgq0 + wtlshaq0)*qsatlsun &
                   - wtaq0*qm - wtgq0*qg - wtlshaq0*qsatlsha )
          etrsun_dtlsun = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun) &
                                    * (wtaq0 + wtgq0 + wtlshaq0)*qsatlsunDT 
          etrsun_dtlsha = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun) &
                                    * ( - wtlshaq0*qsatlshaDT )

          if(etrsun .ge. sigf*etrc*laisun/lai)then
             etrsun = sigf*etrc*laisun/lai
             etrsun_dtlsun = 0.
             etrsun_dtlsha = 0.
          end if

          evplwetsun = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb &
                            * ( (wtaq0 + wtgq0 + wtlshaq0)*qsatlsun &
                            - wtaq0*qm - wtgq0*qg - wtlshaq0*qsatlsha )
          evplwetsun_dtlsun = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb &
                              * (wtaq0 + wtgq0 + wtlshaq0)*qsatlsunDT 
          evplwetsun_dtlsha = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb &
                              * ( - wtlshaq0*qsatlshaDT )

          if(evplwetsun .ge. ldew/deltim*laisun/(lai+sai))then
             evplwetsun = ldew/deltim*laisun/(lai+sai)
             evplwetsun_dtlsun = 0.
             evplwetsun_dtlsha = 0.
          endif

          evplsun = etrsun + evplwetsun
          evplsun_dtlsun = etrsun_dtlsun + evplwetsun_dtlsun
          evplsun_dtlsha = etrsun_dtlsha + evplwetsun_dtlsha

          ! ---------------

          etrsha = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha) &
                   * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlsha &
                   - wtaq0*qm - wtgq0*qg - wtlsunq0*qsatlsun )
          etrsha_dtlsun = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha) &
                          * ( - wtlsunq0*qsatlsunDT )
          etrsha_dtlsha = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha) &
                          * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlshaDT )

          if(etrsha .ge. sigf*etrc*laisha/lai)then
             etrsha = sigf*etrc*laisha/lai
             etrsha_dtlsun = 0.
             etrsha_dtlsha = 0.
          endif

          evplwetsha = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / (rb) &
                       * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlsha &
                       - wtaq0*qm - wtgq0*qg - wtlsunq0*qsatlsun )
          evplwetsha_dtlsun = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / rb &
                              * ( - wtlsunq0*qsatlsunDT )
          evplwetsha_dtlsha = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / rb &
                              * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlshaDT )

          if(evplwetsha .ge. ldew/deltim*(laisha+sai)/(lai+sai))then
             evplwetsha = ldew/deltim*(laisha+sai)/(lai+sai) 
             evplwetsha_dtlsun = 0.
             evplwetsha_dtlsha = 0.
          endif

          evplsha = etrsha + evplwetsha

          evplsha_dtlsun = etrsha_dtlsun + evplwetsha_dtlsun
          evplsha_dtlsha = etrsha_dtlsha + evplwetsha_dtlsha
          
          erresun = 0.
          erresha = 0.

          evplsun_noadj = evplsun
          evplsha_noadj = evplsha
          
          if ( evplsun*evplsun_bef < 0. ) then 
             erresun = -0.9*evplsun
             evplsun =  0.1*evplsun
          endif
          
          if ( evplsha*evplsha_bef < 0. ) then 
             erresha = -0.9*evplsha
             evplsha =  0.1*evplsha
          endif

! functions and their derivatives with respect to temperatures
          ftsun = sabvsun + irabsun - senlsun - hvap*evplsun 
          ftsha = sabvsha + irabsha - senlsha - hvap*evplsha

          dftsunDTlsun = dirabsun_dtlsun - senlsun_dtlsun - hvap*evplsun_dtlsun
          dftsunDTlsha = dirabsun_dtlsha - senlsun_dtlsha - hvap*evplsun_dtlsha

          dftshaDTlsun = dirabsha_dtlsun - senlsha_dtlsun - hvap*evplsha_dtlsun
          dftshaDTlsha = dirabsha_dtlsha - senlsha_dtlsha - hvap*evplsha_dtlsha

!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------
          dtlsun(it) = - (ftsun * (dftshaDTlsha-(laisha+sai)*clai/deltim) - ftsha * dftsunDTlsha) &
                       / ((dftsunDTlsun-laisun*clai/deltim) * (dftshaDTlsha-(laisha+sai)*clai/deltim) &
                       - dftsunDTlsha * dftshaDTlsun)

          dtlsha(it) = - (ftsun * dftshaDTlsun - ftsha * (dftsunDTlsun-laisun*clai/deltim)) &
                       / (dftsunDTlsha * dftshaDTlsun &
                       - (dftsunDTlsun-laisun*clai/deltim) * (dftshaDTlsha-(laisha+sai)*clai/deltim))

          dtlsun_noadj = dtlsun(it)
          dtlsha_noadj = dtlsha(it)

          if(it .le. itmax)then

           ! put brakes on large temperature excursions
             if(abs(dtlsun(it)).gt.delmax) dtlsun(it) = delmax*dtlsun(it)/abs(dtlsun(it))
             if(abs(dtlsha(it)).gt.delmax) dtlsha(it) = delmax*dtlsha(it)/abs(dtlsha(it))

             if(it.ge.2)then
                if(dtlsun(it-1)*dtlsun(it) .le. 0.) dtlsun(it) = 0.5*(dtlsun(it-1) + dtlsun(it))
                if(dtlsha(it-1)*dtlsha(it) .le. 0.) dtlsha(it) = 0.5*(dtlsha(it-1) + dtlsha(it))
             endif

          endif

          tlsun = tlsunbef + dtlsun(it)
          tlsha = tlshabef + dtlsha(it)

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------
          del  = sqrt( dtlsun(it)*dtlsun(it) + dtlsha(it)*dtlsha(it) )
          dele = (    dirabsun_dtlsun * dtlsun(it)  +     dirabsun_dtlsha * dtlsha(it))**2 &
               + (    dirabsha_dtlsun * dtlsun(it)  +     dirabsha_dtlsha * dtlsha(it))**2 &
               + (     senlsun_dtlsun * dtlsun(it)  +      senlsun_dtlsha * dtlsha(it))**2 &
               + (     senlsha_dtlsun * dtlsun(it)  +      senlsha_dtlsha * dtlsha(it))**2 &
               + ( hvap*evplsun_dtlsun * dtlsun(it) + hvap*evplsun_dtlsha * dtlsha(it))**2 &
               + ( hvap*evplsha_dtlsun * dtlsun(it) + hvap*evplsha_dtlsha * dtlsha(it))**2  
          dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
          call qsadv(tlsun,psrf,eisun,deisunDT,qsatlsun,qsatlsunDT)

          call qsadv(tlsha,psrf,eisha,deishaDT,qsatlsha,qsatlshaDT)
 
! update vegetation/ground surface temperature, canopy air temperature, 
! canopy air humidity
          taf = wta0*thm + wtg0*tg + wtlsun0*tlsun + wtlsha0*tlsha

          qaf = wtaq0*qm + wtgq0*qg + wtlsunq0*qsatlsun + wtlshaq0*qsatlsha

! update co2 partial pressure within canopy air
          gah2o = 1.0/raw * tprcor/thm                     ! mol m-2 s-1
          pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) &
                * (assimsun + assimsha - respcsun - respcsha - rsoil)

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------
          dth = thm - taf       
          dqh = qm - qaf

          tstar = vonkar/fh*dth
          qstar = vonkar/fq*dqh

          thvstar = tstar + 0.61*th*qstar
          zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)

          if(zeta .ge. 0.)then                             !stable
             zeta = min(2.,max(zeta,1.e-6))
          else                                             !unstable
             zeta = max(-100.,min(zeta,-1.e-6))
          endif
          obu = zldis/zeta

          if(zeta .ge. 0.)then
             um = max(ur,0.1)
          else
             wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
             wc2 = beta*beta*wc*wc
             um = sqrt(ur*ur+wc2)
          endif

          if(obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
          if(nmozsgn .ge. 4) obu = zldis/(-0.01)
          obuold = obu
 
!-----------------------------------------------------------------------
! Test for convergence
!-----------------------------------------------------------------------
          it = it+1

          if(it .gt. itmin) then
             evplsun_bef = evplsun
             evplsha_bef = evplsha
             det = max(del,del2)
             del = max(dele,dele2)
             if(det .lt. dtmin .AND. del .lt. dlemin) exit 
          endif
 
       END DO                ! ITERATION

! ======================================================================
!      END stability iteration 
! ======================================================================

       z0ma = z0mv
       zol = zeta
       rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! canopy fluxes and total assimilation amd respiration

       if(lai .gt. 0.001) then
          rst = 1./(laisun/rssun + laisha/rssha)
       else
          rssun = 2.0e4 ; rssha = 2.0e4
          assimsun = 0. ; assimsha = 0.
          respcsun = 0. ; respcsha = 0.
          rst = 2.0e4
       endif
       assim = assimsun + assimsha
       respc = respcsun + respcsha + rsoil

       fsenl = senlsun + senlsha + senlsun_dtlsun*dtlsun(it-1) + senlsun_dtlsha*dtlsha(it-1) &
                                 + senlsha_dtlsun*dtlsun(it-1) + senlsha_dtlsha*dtlsha(it-1) &
                                 ! add the imbalanced energy below due to T adjustment to sensibel heat
                                 + (dtlsun_noadj-dtlsun(it-1)) * (laisun*clai/deltim       -dftsunDTlsun) &
                                 + (dtlsun_noadj-dtlsun(it-1)) * (                         -dftshaDTlsun) &
                                 + (dtlsha_noadj-dtlsha(it-1)) * (                         -dftsunDTlsha) &
                                 + (dtlsha_noadj-dtlsha(it-1)) * ((laisha+sai)*clai/deltim -dftshaDTlsha) &
                                 ! add the imbalanced energy below due to q adjustment to sensibel heat
                                 + hvap*erresun + hvap*erresha
          
       etr = etrsun + etrsha     +  etrsun_dtlsun*dtlsun(it-1) +  etrsun_dtlsha*dtlsha(it-1) &
                                 +  etrsha_dtlsun*dtlsun(it-1) +  etrsha_dtlsha*dtlsha(it-1)

       evplwet = evplwetsun + evplwetsun_dtlsun*dtlsun(it-1) + evplwetsun_dtlsha*dtlsha(it-1) &
               + evplwetsha + evplwetsha_dtlsun*dtlsun(it-1) + evplwetsha_dtlsha*dtlsha(it-1) 

       evplsun = evplsun_noadj
       evplsha = evplsha_noadj
       fevpl   = evplsun + evplsha + evplsun_dtlsun*dtlsun(it-1) + evplsun_dtlsha*dtlsha(it-1) &
                                   + evplsha_dtlsun*dtlsun(it-1) + evplsha_dtlsha*dtlsha(it-1)

       elwmax = ldew/deltim
       elwdif = max(0., evplwet-elwmax)
       evplwet = min (evplwet, elwmax)

       fevpl = fevpl - elwdif
       fsenl = fsenl + hvap*elwdif

       ! wind stresses 
       taux  = taux - sigf*rhoair*us/ram
       tauy  = tauy - sigf*rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

       fseng = fseng + cpair*rhoair*cgh*(tg-taf)
       fevpg = fevpg +    rhoair*cgw*(qg-qaf)

!-----------------------------------------------------------------------
! downward (upward) longwave radiation below (above) the canopy
!-----------------------------------------------------------------------

       dlrad = sigf * thermk * frl  &
             + stefnc * ( fac1*tlsunbef**3*(tlsunbef+4.*dtlsun(it-1)) &
             + fac2*tlshabef**3*(tlshabef+4.*dtlsha(it-1)) )

       ulrad = stefnc * ( fac1*tlsunbef**3*(tlsunbef+4.*dtlsun(it-1)) &
             + fac2*tlshabef**3*(tlshabef+4.*dtlsha(it-1)) &
             + sigf*emg*thermk*tg**4 )  

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------
       cgrnds = cgrnds + cpair*rhoair*cgh*(1.-wtg0)
       cgrndl = cgrndl + rhoair*cgw*(1.-wtgq0)*dqgdT
       cgrnd  = cgrnds + cgrndl*htvp

!-----------------------------------------------------------------------
! balance check
!-----------------------------------------------------------------------
       err = sabvsun+sabvsha &
           + irabsun+dirabsun_dtlsun*dtlsun(it-1)+dirabsun_dtlsha*dtlsha(it-1) &
           + irabsha+dirabsha_dtlsun*dtlsun(it-1)+dirabsha_dtlsha*dtlsha(it-1) &
           - fsenl-hvap*fevpl

#if(defined CLMDEBUG) 
       if(abs(err) .gt. .2) write(6,*) 'leaftemtwo.F90 : energy imbalance',err, it-1
#endif
 
!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

       ldew = max(0.,ldew-evplwet*deltim)

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------
       tref = tref + sigf*(thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar))
       qref = qref + sigf*( qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar))


  end subroutine leaftemtwo



  subroutine dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)
       
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! determine fraction of foliage covered by water and
! fraction of foliage that is dry and transpiring
!
!=======================================================================

 use precision
 implicit none

  real(r8), INTENT(in) :: sigf   ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: lai    ! leaf area index  [-]
  real(r8), INTENT(in) :: sai    ! stem area index  [-]
  real(r8), INTENT(in) :: dewmx  ! maximum allowed dew [0.1 mm]
  real(r8), INTENT(in) :: ldew   ! depth of water on foliage [kg/m2/s]

  real(r8), INTENT(out) :: fwet  ! fraction of foliage covered by water [-]
  real(r8), INTENT(out) :: fdry  ! fraction of foliage that is green and dry [-]

  real(r8) lsai                  ! lai + sai
  real(r8) dewmxi                ! inverse of maximum allowed dew [1/mm]
  real(r8) vegt                  ! sigf*lsai
!
!-----------------------------------------------------------------------
! Fwet is the fraction of all vegetation surfaces which are wet 
! including stem area which contribute to evaporation
      lsai = lai + sai
      dewmxi  = 1.0/dewmx
      vegt   =  sigf*lsai

      fwet = 0
      if(ldew > 0.) then
         fwet = ((dewmxi/vegt)*ldew)**.666666666666

! Check for maximum limit of fwet
         fwet = min(fwet,1.0)

      end if

! fdry is the fraction of lai which is dry because only leaves can 
! transpire. Adjusted for stem area which does not transpire
      fdry = (1.-fwet)*lai/lsai


  end subroutine dewfraction
!----------------------------------------------------------------------         



END MODULE LEAF_temperature
