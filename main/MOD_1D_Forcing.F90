MODULE MOD_1D_Forcing
! -------------------------------
! Meteorogical Forcing
!
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
IMPLICIT NONE
SAVE

! -----------------------------------------------------------------
  real(r8), allocatable :: forc_pco2m (:) ! CO2 concentration in atmos. (pascals)
  real(r8), allocatable :: forc_po2m  (:) ! O2 concentration in atmos. (pascals)
  real(r8), allocatable :: forc_us    (:) ! wind in eastward direction [m/s]
  real(r8), allocatable :: forc_vs    (:) ! wind in northward direction [m/s]
  real(r8), allocatable :: forc_t     (:) ! temperature at reference height [kelvin]
  real(r8), allocatable :: forc_q     (:) ! specific humidity at reference height [kg/kg]
  real(r8), allocatable :: forc_prc   (:) ! convective precipitation [mm/s]
  real(r8), allocatable :: forc_prl   (:) ! large scale precipitation [mm/s]
  real(r8), allocatable :: forc_rain  (:) ! rain [mm/s]
  real(r8), allocatable :: forc_snow  (:) ! snow [mm/s]
  real(r8), allocatable :: forc_psrf  (:) ! atmospheric pressure at the surface [pa]
  real(r8), allocatable :: forc_pbot  (:) ! atm bottom level pressure (or reference height) (pa)
  real(r8), allocatable :: forc_sols  (:) ! atm vis direct beam solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_soll  (:) ! atm nir direct beam solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_solsd (:) ! atm vis diffuse solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_solld (:) ! atm nir diffuse solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_frl   (:) ! atmospheric infrared (longwave) radiation [W/m2]
  real(r8), allocatable :: forc_hgt_u (:) ! observational height of wind [m]
  real(r8), allocatable :: forc_hgt_t (:) ! observational height of temperature [m]
  real(r8), allocatable :: forc_hgt_q (:) ! observational height of humidity [m]
  real(r8), allocatable :: forc_rhoair(:) ! air density [kg/m3]

! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_1D_Forcing
      public :: deallocate_1D_Forcing

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_Forcing (numpatch)
! ------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------
  use precision
  IMPLICIT NONE
  integer, INTENT(in) :: numpatch

      allocate ( forc_pco2m  (numpatch) ) ! CO2 concentration in atmos. (pascals)
      allocate ( forc_po2m   (numpatch) ) ! O2 concentration in atmos. (pascals)
      allocate ( forc_us     (numpatch) ) ! wind in eastward direction [m/s]
      allocate ( forc_vs     (numpatch) ) ! wind in northward direction [m/s]
      allocate ( forc_t      (numpatch) ) ! temperature at reference height [kelvin]
      allocate ( forc_q      (numpatch) ) ! specific humidity at reference height [kg/kg]
      allocate ( forc_prc    (numpatch) ) ! convective precipitation [mm/s]
      allocate ( forc_prl    (numpatch) ) ! large scale precipitation [mm/s]
      allocate ( forc_rain   (numpatch) ) ! rain [mm/s]
      allocate ( forc_snow   (numpatch) ) ! snow [mm/s]
      allocate ( forc_psrf   (numpatch) ) ! atmospheric pressure at the surface [pa]
      allocate ( forc_pbot   (numpatch) ) ! atm bottom level pressure (or reference height) (pa)
      allocate ( forc_sols   (numpatch) ) ! atm vis direct beam solar rad onto srf [W/m2]
      allocate ( forc_soll   (numpatch) ) ! atm nir direct beam solar rad onto srf [W/m2]
      allocate ( forc_solsd  (numpatch) ) ! atm vis diffuse solar rad onto srf [W/m2]
      allocate ( forc_solld  (numpatch) ) ! atm nir diffuse solar rad onto srf [W/m2]
      allocate ( forc_frl    (numpatch) ) ! atmospheric infrared (longwave) radiation [W/m2]
      allocate ( forc_hgt_u  (numpatch) ) ! observational height of wind [m]
      allocate ( forc_hgt_t  (numpatch) ) ! observational height of temperature [m]
      allocate ( forc_hgt_q  (numpatch) ) ! observational height of humidity [m]
      allocate ( forc_rhoair (numpatch) ) ! air density [kg/m3]

  END SUBROUTINE allocate_1D_Forcing


  SUBROUTINE deallocate_1D_Forcing
      deallocate ( forc_pco2m  ) ! CO2 concentration in atmos. (pascals)
      deallocate ( forc_po2m   ) ! O2 concentration in atmos. (pascals)
      deallocate ( forc_us     ) ! wind in eastward direction [m/s]
      deallocate ( forc_vs     ) ! wind in northward direction [m/s]
      deallocate ( forc_t      ) ! temperature at reference height [kelvin]
      deallocate ( forc_q      ) ! specific humidity at reference height [kg/kg]
      deallocate ( forc_prc    ) ! convective precipitation [mm/s]
      deallocate ( forc_prl    ) ! large scale precipitation [mm/s]
      deallocate ( forc_rain   ) ! rain [mm/s]
      deallocate ( forc_snow   ) ! snow [mm/s]
      deallocate ( forc_psrf   ) ! atmospheric pressure at the surface [pa]
      deallocate ( forc_pbot   ) ! atm bottom level pressure (or reference height) (pa)
      deallocate ( forc_sols   ) ! atm vis direct beam solar rad onto srf [W/m2]
      deallocate ( forc_soll   ) ! atm nir direct beam solar rad onto srf [W/m2]
      deallocate ( forc_solsd  ) ! atm vis diffuse solar rad onto srf [W/m2]
      deallocate ( forc_solld  ) ! atm nir diffuse solar rad onto srf [W/m2]
      deallocate ( forc_frl    ) ! atmospheric infrared (longwave) radiation [W/m2]
      deallocate ( forc_hgt_u  ) ! observational height of wind [m]
      deallocate ( forc_hgt_t  ) ! observational height of temperature [m]
      deallocate ( forc_hgt_q  ) ! observational height of humidity [m]
      deallocate ( forc_rhoair ) ! air density [kg/m3]
  END SUBROUTINE deallocate_1D_Forcing

END MODULE MOD_1D_Forcing
! ------ EOP --------
