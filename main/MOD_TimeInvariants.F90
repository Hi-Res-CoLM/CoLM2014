MODULE MOD_TimeInvariants 
! -------------------------------
!
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! surface classification and soil information
      integer,  allocatable :: ixy_patch(:)  ! patch longitude index
      integer,  allocatable :: jxy_patch(:)  ! patch latitude index
      integer,  allocatable :: mxy_patch(:)  ! index of land cover type of the patches at the fraction > 0
      integer,  allocatable :: spatch_xy(:,:)! start patch number of grid
      integer,  allocatable :: epatch_xy(:,:)! end patch number of grid
      real(r8), allocatable :: wtxy_patch(:) ! patch weight
      real(r8), allocatable :: area_patch(:) ! area of grid cell the patch located in

      real(r8), allocatable :: dlat    (:) ! latitude in radians
      real(r8), allocatable :: dlon    (:) ! longitude in radians
      real(r8), allocatable :: lats    (:) ! latitude in degrees
      real(r8), allocatable :: lons    (:) ! longitude in degrees
      integer , allocatable :: itypwat (:) ! land water type

      real(r8), allocatable :: z_soi (:,:) ! node depth [m]
      real(r8), allocatable :: dz_soi(:,:) ! interface depth [m]
      real(r8), allocatable :: lakedepth(:) !
      real(r8), allocatable :: dz_lake(:,:) ! new lake scheme

      real(r8), allocatable :: soil_s_v_alb(:) ! albedo of visible of the saturated soil
      real(r8), allocatable :: soil_d_v_alb(:) ! albedo of visible of the dry soil
      real(r8), allocatable :: soil_s_n_alb(:) ! albedo of near infrared of the saturated soil
      real(r8), allocatable :: soil_d_n_alb(:) ! albedo of near infrared of the dry soil
      real(r8), allocatable :: porsl (:,:) ! fraction of soil that is voids [-]
      real(r8), allocatable :: psi0  (:,:) ! minimum soil suction [mm] (NOTE: "-" valued)
      real(r8), allocatable :: bsw   (:,:) ! clapp and hornbereger "b" parameter [-]
      real(r8), allocatable :: hksati(:,:) ! hydraulic conductivity at saturation [mm h2o/s]
      real(r8), allocatable :: csol  (:,:) ! heat capacity of soil solids [J/(m3 K)]
      real(r8), allocatable :: dksatu(:,:) ! thermal conductivity of saturated soil [W/m-K]
      real(r8), allocatable :: dkdry (:,:) ! thermal conductivity for dry soil  [W/(m-K)]
      real(r8), allocatable :: rootfr(:,:) ! fraction of roots in each soil layer

      real(r8), allocatable :: z0m     (:) ! aerodynamic roughness length [m]
      real(r8), allocatable :: displa  (:) ! displacement height [m]
      real(r8), allocatable :: sqrtdi  (:) ! inverse sqrt of leaf dimension [m**-0.5]
      real(r8), allocatable :: effcon  (:) ! quantum efficiency of RuBP regeneration 
      real(r8), allocatable :: vmax25  (:) ! maximum carboxylation rate at 25 C at canopy top
      real(r8), allocatable :: slti    (:) ! s3: slope of low temperature inhibition function     
      real(r8), allocatable :: hlti    (:) ! s4: 1/2 point of low temperature inhibition function
      real(r8), allocatable :: shti    (:) ! s1: slope of high temperature inhibition function  
      real(r8), allocatable :: hhti    (:) ! s2: 1/2 point of high temperature inhibition function 
      real(r8), allocatable :: trda    (:) ! s5: temperature coefficient in gs-a model            
      real(r8), allocatable :: trdm    (:) ! s6: temperature coefficient in gs-a model           
      real(r8), allocatable :: trop    (:) ! temperature coefficient in gs-a model          
      real(r8), allocatable :: gradm   (:) ! conductance-photosynthesis slope parameter
      real(r8), allocatable :: binter  (:) ! conductance-photosynthesis intercep
      real(r8), allocatable :: extkn   (:) ! coefficient of leaf nitrogen allocation
      real(r8), allocatable :: chil    (:) ! leaf angle distribution factor
      real(r8), allocatable :: ref (:,:,:) ! leaf reflectance (iw=iband, il=life and dead)
      real(r8), allocatable :: tran(:,:,:) ! leaf transmittance (iw=iband, il=life and dead)

      real(r8) zlnd    ! roughness length for soil [m]
      real(r8) zsno    ! roughness length for snow [m]
      real(r8) csoilc  ! drag coefficient for soil under canopy [-]
      real(r8) dewmx   ! maximum dew
      real(r8) wtfact  ! fraction of model area with high water table
      real(r8) capr    ! tuning factor to turn first layer T into surface T
      real(r8) cnfac   ! Crank Nicholson factor between 0 and 1
      real(r8) ssi     ! irreducible water saturation of snow
      real(r8) wimp    ! water impremeable if porosity less than wimp
      real(r8) pondmx  ! ponding depth (mm)
      real(r8) smpmax  ! wilting point potential in mm
      real(r8) smpmin  ! restriction for min of soil poten. (mm)
      real(r8) trsmx0  ! max transpiration for moist soil+100% veg.  [mm/s]
      real(r8) tcrit   ! critical temp. to determine rain or snow

      real(r8), parameter :: spval = -1.e36_r8  ! a special value for missing and filling use

! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_TimeInvariants
      public :: deallocate_TimeInvariants
      public :: WRITE_TimeInvariants
      public :: READ_TimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeInvariants (nl_soil,nl_lake,numpatch, &
                                      lon_points, lat_points)
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------

  use precision
  IMPLICIT NONE
  integer, INTENT(in) :: nl_soil
  integer, INTENT(in) :: nl_lake
  integer, INTENT(in) :: numpatch
  integer, INTENT(in) :: lon_points
  integer, INTENT(in) :: lat_points

      allocate (ixy_patch            (numpatch))
      allocate (jxy_patch            (numpatch))
      allocate (mxy_patch            (numpatch))
      allocate (spatch_xy            (lon_points, lat_points))
      allocate (epatch_xy            (lon_points, lat_points))
      allocate (wtxy_patch           (numpatch))
      allocate (area_patch           (numpatch))

      allocate (itypwat              (numpatch))
      allocate (dlat                 (numpatch))
      allocate (dlon                 (numpatch))
      allocate (lats                 (lat_points))
      allocate (lons                 (lon_points))

      allocate (z_soi        (nl_soil,numpatch))
      allocate (dz_soi       (nl_soil,numpatch))
      allocate (lakedepth            (numpatch))
      allocate (dz_lake      (nl_lake,numpatch))

      allocate (soil_s_v_alb         (numpatch))
      allocate (soil_d_v_alb         (numpatch))
      allocate (soil_s_n_alb         (numpatch))
      allocate (soil_d_n_alb         (numpatch))

      allocate (porsl        (nl_soil,numpatch))
      allocate (psi0         (nl_soil,numpatch))
      allocate (bsw          (nl_soil,numpatch))
      allocate (hksati       (nl_soil,numpatch))
      allocate (csol         (nl_soil,numpatch))
      allocate (dksatu       (nl_soil,numpatch))
      allocate (dkdry        (nl_soil,numpatch))
      allocate (rootfr       (nl_soil,numpatch))

      allocate (z0m                  (numpatch))
      allocate (displa               (numpatch))
      allocate (sqrtdi               (numpatch))
      allocate (effcon               (numpatch))
      allocate (vmax25               (numpatch))
      allocate (slti                 (numpatch))
      allocate (hlti                 (numpatch))
      allocate (shti                 (numpatch))
      allocate (hhti                 (numpatch))
      allocate (trda                 (numpatch))
      allocate (trdm                 (numpatch))
      allocate (trop                 (numpatch))
      allocate (gradm                (numpatch))
      allocate (binter               (numpatch))
      allocate (extkn                (numpatch))
      allocate (chil                 (numpatch))
      allocate (ref              (2,2,numpatch))
      allocate (tran             (2,2,numpatch))

  END SUBROUTINE allocate_TimeInvariants


  SUBROUTINE READ_TimeInvariants(dir_restart_hist,site)
! --------------------------------------------------------------------
! Write out as a restart file [histTimeConst]
! ...............................................
  use precision
  IMPLICIT NONE
      character(LEN=256), INTENT(in) :: site           ! site name
      character(LEN=256), INTENT(in) :: dir_restart_hist

      character(LEN=256) :: fhistTimeConst
      integer :: lhistTimeConst


      lhistTimeConst = 100
      fhistTimeConst = trim(dir_restart_hist)//trim(site)//'-'//'rstTimeConst'
      OPEN(unit=lhistTimeConst,file=trim(fhistTimeConst),status='unknown',&
                               form='unformatted',action='read')
      READ (lhistTimeConst)  & !
            ixy_patch,       & ! longitude index for each patch point
            jxy_patch,       & ! latitude index for each patch point
            mxy_patch,       & ! index of land cover type of the patches at the fraction > 0
            spatch_xy,       & ! start patch number of grid
            epatch_xy,       & ! end patch number of grid
            wtxy_patch,      & ! subgrid weight for each patch point
            area_patch         ! area of grid the patch located in

      READ (lhistTimeConst)  & !
            dlat    ,        & ! latitude in radians
            dlon    ,        & ! longitude in radians
            lats    ,        & ! latitude in degrees
            lons    ,        & ! longitude in degrees
            itypwat ,        & ! land water type

      ! Soil and plant parameters OF CLM
            z_soi   ,        & !
            dz_soi  ,        & !
            lakedepth,       & !
            dz_lake ,        & !

            soil_s_v_alb,    & ! albedo of visible of the saturated soil
            soil_d_v_alb,    & ! albedo of visible of the dry soil
            soil_s_n_alb,    & ! albedo of near infrared of the saturated soil
            soil_d_n_alb,    & ! albedo of near infrared of the dry soil

            porsl   ,        & ! fraction of soil that is voids [-]
            psi0    ,        & ! minimum soil suction [mm] (NOTE: "-" valued)
            bsw     ,        & ! clapp and hornbereger "b" parameter [-]
            hksati  ,        & ! hydraulic conductivity at saturation [mm h2o/s]
            csol    ,        & ! heat capacity of soil solids [J/(m3 K)]
            dksatu  ,        & ! thermal conductivity of saturated soil [W/m-K]
            dkdry   ,        & ! thermal conductivity for dry soil  [W/(m-K)]
            rootfr  ,        & ! fraction of roots in each soil layer

            z0m     ,        & ! aerodynamic roughness length [m]
            displa  ,        & ! displacement height [m]
            sqrtdi  ,        & ! inverse sqrt of leaf dimension [m**-0.5]
            effcon  ,        & ! quantum efficiency of RuBP regeneration
            vmax25  ,        & ! maximum carboxylation rate at 25 C at canopy top
            slti    ,        & ! s3: slope of low temperature inhibition function
            hlti    ,        & ! s4: 1/2 point of low temperature inhibition function
            shti    ,        & ! s1: slope of high temperature inhibition function
            hhti    ,        & ! s2: 1/2 point of high temperature inhibition function
            trda    ,        & ! s5: temperature coefficient in gs-a model
            trdm    ,        & ! s6: temperature coefficient in gs-a model
            trop    ,        & ! temperature coefficient in gs-a model
            gradm   ,        & ! conductance-photosynthesis slope parameter
            binter  ,        & ! conductance-photosynthesis intercep
            extkn   ,        & ! coefficient of leaf nitrogen allocation
            chil    ,        & ! leaf angle distribution factor
            ref     ,        & ! leaf reflectance (iw=iband, il=life and dead)
            tran    ,        & ! leaf transmittance (iw=iband, il=life and dead)

      ! CLM TUNABLE constants
            zlnd    ,        & ! roughness length for soil [m]
            zsno    ,        & ! roughness length for snow [m]
            csoilc  ,        & ! drag coefficient for soil under canopy [-]
            dewmx   ,        & ! maximum dew
            wtfact  ,        & ! fraction of model area with high water table
            capr    ,        & ! tuning factor to turn first layer T into surface T
            cnfac   ,        & ! Crank Nicholson factor between 0 and 1
            ssi     ,        & ! irreducible water saturation of snow
            wimp    ,        & ! water impremeable if porosity less than wimp
            pondmx  ,        & ! ponding depth (mm)
            smpmax  ,        & ! wilting point potential in mm
            smpmin  ,        & ! restriction for min of soil poten. (mm)
            trsmx0  ,        & ! max transpiration for moist soil+100% veg.  [mm/s]
            tcrit              ! critical temp. to determine rain or snow

      CLOSE(lhistTimeConst)

  END SUBROUTINE READ_TimeInvariants


  SUBROUTINE WRITE_TimeInvariants(dir_restart_hist,site)
! --------------------------------------------------------------------
! Write out as a restart file [histTimeConst]
! ...............................................
  use precision
  IMPLICIT NONE
      character(LEN=256), INTENT(in) :: site           ! site name
      character(LEN=256), INTENT(in) :: dir_restart_hist

      character(LEN=256) :: fhistTimeConst
      integer :: lhistTimeConst

      lhistTimeConst = 100
      fhistTimeConst = trim(dir_restart_hist)//trim(site)//'-'//'rstTimeConst'
      OPEN(unit=lhistTimeConst,file=trim(fhistTimeConst),status='unknown',&
                               form='unformatted',action='write')
      WRITE(lhistTimeConst)  & !
            ixy_patch,       & ! longitude index for each patch point
            jxy_patch,       & ! latitude index for each patch point
            mxy_patch,       & ! index of land cover type of the patches at the fraction > 0
            spatch_xy,       & ! start patch number of grid
            epatch_xy,       & ! end patch number of grid
            wtxy_patch,      & ! subgrid weight for each patch point
            area_patch         ! area of grid cell the patch located in

      WRITE(lhistTimeConst)  & !
            dlat    ,        & ! latitude in radians
            dlon    ,        & ! longitude in radians
            lats    ,        & ! latitude in degrees
            lons    ,        & ! longitude in degrees
            itypwat ,        & ! land water type

      ! Soil and plant parameters OF CLM
            z_soi   ,        & !
            dz_soi  ,        & !
            lakedepth,       & !
            dz_lake ,        & !

            soil_s_v_alb,    & ! albedo of visible of the saturated soil
            soil_d_v_alb,    & ! albedo of visible of the dry soil
            soil_s_n_alb,    & ! albedo of near infrared of the saturated soil
            soil_d_n_alb,    & ! albedo of near infrared of the dry soil

            porsl   ,        & ! fraction of soil that is voids [-]
            psi0    ,        & ! minimum soil suction [mm] (NOTE: "-" valued)
            bsw     ,        & ! clapp and hornbereger "b" parameter [-]
            hksati  ,        & ! hydraulic conductivity at saturation [mm h2o/s]
            csol    ,        & ! heat capacity of soil solids [J/(m3 K)]
            dksatu  ,        & ! thermal conductivity of saturated soil [W/m-K]
            dkdry   ,        & ! thermal conductivity for dry soil  [W/(m-K)]
            rootfr  ,        & ! fraction of roots in each soil layer

            z0m     ,        & ! aerodynamic roughness length [m] 
            displa  ,        & ! displacement height [m]
            sqrtdi  ,        & ! inverse sqrt of leaf dimension [m**-0.5]
            effcon  ,        & ! quantum efficiency of RuBP regeneration
            vmax25  ,        & ! maximum carboxylation rate at 25 C at canopy top
            slti    ,        & ! s3: slope of low temperature inhibition function
            hlti    ,        & ! s4: 1/2 point of low temperature inhibition function
            shti    ,        & ! s1: slope of high temperature inhibition function
            hhti    ,        & ! s2: 1/2 point of high temperature inhibition function
            trda    ,        & ! s5: temperature coefficient in gs-a model
            trdm    ,        & ! s6: temperature coefficient in gs-a model
            trop    ,        & ! temperature coefficient in gs-a model
            gradm   ,        & ! conductance-photosynthesis slope parameter
            binter  ,        & ! conductance-photosynthesis intercep
            extkn   ,        & ! coefficient of leaf nitrogen allocation
            chil    ,        & ! leaf angle distribution factor
            ref     ,        & ! leaf reflectance (iw=iband, il=life and dead)
            tran    ,        & ! leaf transmittance (iw=iband, il=life and dead)

      ! CLM TUNABLE constants
            zlnd    ,        & ! roughness length for soil [m]
            zsno    ,        & ! roughness length for snow [m]
            csoilc  ,        & ! drag coefficient for soil under canopy [-]
            dewmx   ,        & ! maximum dew
            wtfact  ,        & ! fraction of model area with high water table
            capr    ,        & ! tuning factor to turn first layer T into surface T
            cnfac   ,        & ! Crank Nicholson factor between 0 and 1
            ssi     ,        & ! irreducible water saturation of snow
            wimp    ,        & ! water impremeable if porosity less than wimp
            pondmx  ,        & ! ponding depth (mm)
            smpmax  ,        & ! wilting point potential in mm
            smpmin  ,        & ! restriction for min of soil poten. (mm)
            trsmx0  ,        & ! max transpiration for moist soil+100% veg.  [mm/s]
            tcrit              ! critical temp. to determine rain or snow

      CLOSE(lhistTimeConst)
  END SUBROUTINE WRITE_TimeInvariants


  SUBROUTINE deallocate_TimeInvariants
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------
      deallocate (ixy_patch )
      deallocate (jxy_patch )
      deallocate (mxy_patch )
      deallocate (spatch_xy )
      deallocate (epatch_xy )
      deallocate (wtxy_patch)
      deallocate (area_patch)

      deallocate (dlat   )
      deallocate (dlon   )
      deallocate (lats   )
      deallocate (lons   )
      deallocate (itypwat)

      deallocate (z_soi  )
      deallocate (dz_soi )
      deallocate (lakedepth)
      deallocate (dz_lake )

      deallocate (soil_s_v_alb)
      deallocate (soil_d_v_alb)
      deallocate (soil_s_n_alb)
      deallocate (soil_d_n_alb)

      deallocate (porsl  )
      deallocate (psi0   )
      deallocate (bsw    )
      deallocate (hksati )
      deallocate (csol   )
      deallocate (dksatu )
      deallocate (dkdry  )
      deallocate (rootfr )

      deallocate (z0m    )
      deallocate (displa )
      deallocate (sqrtdi )
      deallocate (effcon )
      deallocate (vmax25 )
      deallocate (slti   )
      deallocate (hlti   )
      deallocate (shti   )
      deallocate (hhti   )
      deallocate (trda   )
      deallocate (trdm   )
      deallocate (trop   )
      deallocate (gradm  )
      deallocate (binter )
      deallocate (extkn  )
      deallocate (chil   )
      deallocate (ref    )
      deallocate (tran   )

  END SUBROUTINE deallocate_TimeInvariants

END MODULE MOD_TimeInvariants
