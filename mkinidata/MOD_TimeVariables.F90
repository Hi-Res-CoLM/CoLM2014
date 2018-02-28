MODULE MOD_TimeVariables
! -------------------------------
!
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
use timemanager

IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
      real(r8), allocatable :: z_sno (:,:)  ! node depth [m]
      real(r8), allocatable :: dz_sno(:,:)  ! interface depth [m]
      real(r8), allocatable :: t_soisno (:,:) ! soil temperature [K]
      real(r8), allocatable :: wliq_soisno(:,:) ! liquid water in layers [kg/m2]
      real(r8), allocatable :: wice_soisno(:,:) ! ice lens in layers [kg/m2]
      real(r8), allocatable :: h2osoi(:,:)  ! volumetric soil water in layers [m3/m3]
      real(r8), allocatable :: rstfac   (:) ! factor of soil water stress 
      real(r8), allocatable :: t_grnd   (:) ! ground surface temperature [K]
      real(r8), allocatable :: tlsun    (:) ! sunlit leaf temperature [K]
      real(r8), allocatable :: tlsha    (:) ! shaded leaf temperature [K]
      real(r8), allocatable :: ldew     (:) ! depth of water on foliage [mm]
      real(r8), allocatable :: sag      (:) ! non dimensional snow age [-]
      real(r8), allocatable :: scv      (:) ! snow cover, water equivalent [mm]
      real(r8), allocatable :: snowdp   (:) ! snow depth [meter]
      real(r8), allocatable :: fveg     (:) ! fraction of vegetation cover
      real(r8), allocatable :: fsno     (:) ! fraction of snow cover on ground
      real(r8), allocatable :: sigf     (:) ! fraction of veg cover, excluding snow-covered veg [-]
      real(r8), allocatable :: green    (:) ! leaf greenness
      real(r8), allocatable :: lai      (:) ! leaf area index
      real(r8), allocatable :: laisun   (:) ! leaf area index
      real(r8), allocatable :: laisha   (:) ! leaf area index
      real(r8), allocatable :: sai      (:) ! stem area index
      real(r8), allocatable :: coszen   (:) ! cosine of solar zenith angle
      real(r8), allocatable :: albg (:,:,:) ! albedo, ground [-]
      real(r8), allocatable :: albv (:,:,:) ! albedo, vegetation [-]
      real(r8), allocatable :: alb  (:,:,:) ! averaged albedo [-]
      real(r8), allocatable :: ssun (:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: ssha (:,:,:) ! shaded canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: thermk   (:) ! canopy gap fraction for tir radiation
      real(r8), allocatable :: extkb    (:) ! (k, g(mu)/mu) direct solar extinction coefficient
      real(r8), allocatable :: extkd    (:) ! diffuse and scattered diffuse PAR extinction coefficient
      real(r8), allocatable :: zwt      (:) ! the depth to water table [m]
      real(r8), allocatable :: wa       (:) ! water storage in aquifer [mm]
      real(r8), allocatable :: wat      (:) ! total water storage [mm]

      real(r8), allocatable :: t_lake(:,:)      ! lake layer teperature [K]
      real(r8), allocatable :: lake_icefrac(:,:)! lake mass fraction of lake layer that is frozen

      real(r8), allocatable :: trad     (:) ! radiative temperature of surface [K]
      real(r8), allocatable :: tref     (:) ! 2 m height air temperature [kelvin]
      real(r8), allocatable :: qref     (:) ! 2 m height air specific humidity
      real(r8), allocatable :: rst      (:) ! canopy stomatal resistance (s/m)
      real(r8), allocatable :: emis     (:) ! averaged bulk surface emissivity
      real(r8), allocatable :: z0ma     (:) ! effective roughness [m]
      real(r8), allocatable :: zol      (:) ! dimensionless height (z/L) used in Monin-Obukhov theory
      real(r8), allocatable :: rib      (:) ! bulk Richardson number in surface layer
      real(r8), allocatable :: ustar    (:) ! u* in similarity theory [m/s]
      real(r8), allocatable :: qstar    (:) ! q* in similarity theory [kg/kg]
      real(r8), allocatable :: tstar    (:) ! t* in similarity theory [K]
      real(r8), allocatable :: fm       (:) ! integral of profile function for momentum
      real(r8), allocatable :: fh       (:) ! integral of profile function for heat
      real(r8), allocatable :: fq       (:) ! integral of profile function for moisture

! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_TimeVariables
      public :: READ_TimeVariables
      public :: WRITE_TimeVariables
      public :: deallocate_TimeVariables

! PRIVATE MEMBER FUNCTIONS:


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables (nl_soil,nl_lake,maxsnl,numpatch)
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------------

  use precision
  IMPLICIT NONE
      integer, INTENT(in) :: nl_soil
      integer, INTENT(in) :: nl_lake
      integer, INTENT(in) :: maxsnl
      integer, INTENT(in) :: numpatch

      allocate (z_sno  (maxsnl+1:0,numpatch))
      allocate (dz_sno (maxsnl+1:0,numpatch))
      allocate (t_soisno  (maxsnl+1:nl_soil,numpatch))
      allocate (wliq_soisno(maxsnl+1:nl_soil,numpatch))
      allocate (wice_soisno(maxsnl+1:nl_soil,numpatch))
      allocate (h2osoi     (1:nl_soil,numpatch))
      allocate (rstfac               (numpatch))
      allocate (t_grnd               (numpatch))
      allocate (tlsun                (numpatch))
      allocate (tlsha                (numpatch))
      allocate (ldew                 (numpatch))
      allocate (sag                  (numpatch))
      allocate (scv                  (numpatch))
      allocate (snowdp               (numpatch))
      allocate (fveg                 (numpatch))
      allocate (fsno                 (numpatch))
      allocate (sigf                 (numpatch))
      allocate (green                (numpatch))
      allocate (lai                  (numpatch))
      allocate (laisun               (numpatch))
      allocate (laisha               (numpatch))
      allocate (sai                  (numpatch))
      allocate (coszen               (numpatch))
      allocate (albg             (2,2,numpatch))
      allocate (albv             (2,2,numpatch))
      allocate (alb              (2,2,numpatch))
      allocate (ssun             (2,2,numpatch))
      allocate (ssha             (2,2,numpatch))
      allocate (thermk               (numpatch))
      allocate (extkb                (numpatch))
      allocate (extkd                (numpatch))
      allocate (zwt                  (numpatch))
      allocate (wa                   (numpatch))
      allocate (wat                  (numpatch))

      allocate (t_lake       (nl_lake,numpatch))    !new lake scheme
      allocate (lake_icefrac (nl_lake,numpatch))    !new lake scheme

      allocate (trad                 (numpatch))
      allocate (tref                 (numpatch))
      allocate (qref                 (numpatch))
      allocate (rst                  (numpatch))
      allocate (emis                 (numpatch))
      allocate (z0ma                 (numpatch))
      allocate (zol                  (numpatch))
      allocate (rib                  (numpatch))
      allocate (ustar                (numpatch))
      allocate (qstar                (numpatch))
      allocate (tstar                (numpatch))
      allocate (fm                   (numpatch))
      allocate (fh                   (numpatch))
      allocate (fq                   (numpatch))

  END SUBROUTINE allocate_TimeVariables


  SUBROUTINE READ_TimeVariables (idate,dir_restart_hist,site)
! --------------------------------------------------------------------
! Read the model variables for restart run [histTimeVar]
! ...............................................................

  use precision
  IMPLICIT NONE
      character(LEN=255), INTENT(in) :: dir_restart_hist
      character(LEN=256), INTENT(in) :: site
      integer, INTENT(in) :: idate(3)     ! calendar (year, julian day, seconds)

      integer :: lhistTimeVar             ! logical unit number of restart time-varying file
      integer :: id(3)                    ! calendar (year, julian day, seconds)
      character(LEN=255) :: cdate         ! character for date
      character(LEN=256) :: fhistTimeVar  ! file name of time-varying file

    ! the model variables for restart run
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)

      lhistTimeVar = 100
      fhistTimeVar = trim(dir_restart_hist)//trim(site)//'-'//'rstTimeVar'//'-'//trim(cdate)
      print*,trim(fhistTimeVar)
      OPEN(unit=lhistTimeVar,file=trim(fhistTimeVar),status='unknown',&
                             form='unformatted',action='read')
      READ (lhistTimeVar) id, & !
          ! Time-varying state variables which reaquired by restart run
            z_sno   , & !    node depth [m]
            dz_sno  , & !    interface depth [m]
            t_soisno   , & ! soil temperature [K]
            wliq_soisno, & ! liquid water in layers [kg/m2]
            wice_soisno, & ! ice lens in layers [kg/m2]
            t_grnd  , & !    ground surface temperature [K]
            tlsun   , & !    sunlit leaf temperature [K]
            tlsha   , & !    shaded leaf temperature [K]
            ldew    , & !    depth of water on foliage [mm]
            sag     , & !    non dimensional snow age [-]
            scv     , & !    snow cover, water equivalent [mm]
            snowdp  , & !    snow depth [meter]
            fveg    , & !    fraction of vegetation cover
            fsno    , & !    fraction of snow cover on ground
            sigf    , & !    fraction of veg cover, excluding snow-covered veg [-]
            green   , & !    leaf greenness
            lai     , & !    leaf area index
            sai     , & !    stem area index
            coszen  , & !    cosine of solar zenith angle
            albg    , & !    albedo, ground [-]
            albv    , & !    albedo, vegetation [-]
            alb     , & !    averaged albedo [-]
            ssun    , & !    sunlit canopy absorption for solar radiation (0-1)
            ssha    , & !    shaded canopy absorption for solar radiation (0-1)
            thermk  , & !    canopy gap fraction for tir radiation
            extkb   , & !    (k, g(mu)/mu) direct solar extinction coefficient
            extkd   , & !    diffuse and scattered diffuse PAR extinction coefficient
            zwt     , & !    the depth to water table [m]
            wa      , & !    water storage in aquifer [mm]

            t_lake       , & !
            lake_icefrac , & !

          ! Additional variables required by reginal model (such as WRF & RSM)
            trad    , & !    radiative temperature of surface [K]
            tref    , & !    2 m height air temperature [kelvin]
            qref    , & !    2 m height air specific humidity
            rst     , & !    canopy stomatal resistance (s/m)
            emis    , & !    averaged bulk surface emissivity
            z0ma    , & !    effective roughness [m]
            zol     , & !    dimensionless height (z/L) used in Monin-Obukhov theory
            rib     , & !    bulk Richardson number in surface layer
            ustar   , & !    u* in similarity theory [m/s]
            qstar   , & !    q* in similarity theory [kg/kg]
            tstar   , & !    t* in similarity theory [K]
            fm      , & !    integral of profile function for momentum
            fh      , & !    integral of profile function for heat
            fq          !    integral of profile function for moisture

            if(id(1) /= idate(1) .or. id(2) /= idate(2) .or. id(3) /= idate(3))then
               print*, 'id = ', id, 'idate = ', idate
               print*, 'The date of initial data is NOT IDENTICAL TO initial set-up'
               call abort
            endif

      CLOSE(lhistTimeVar)

  END SUBROUTINE READ_TimeVariables


  SUBROUTINE WRITE_TimeVariables (idate,dir_restart_hist,site)
! --------------------------------------------------------------------
! Write out the model variables for restart run [histTimeVar]
! --------------------------------------------------------------------

  use precision
  IMPLICIT NONE
      integer, INTENT(in) :: idate(3)     ! calendar (year, julian day, seconds)
      character(LEN=255), INTENT(in) :: dir_restart_hist
      character(LEN=256), INTENT(in) :: site

      integer :: lhistTimeVar             ! logical unit number of restart time-varying file
      integer :: id(3)                    ! calendar (year, julian day, seconds), temporal
      character(LEN=255) :: cdate         ! character for date
      character(LEN=256) :: fhistTimeVar  ! file name of time-varying file

! ...............................................................

      id(:) = idate(:)
      call adj2begin(id)

    ! the model variables for restart run 
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') id(1), id(2), id(3)

      lhistTimeVar = 100
      fhistTimeVar = trim(dir_restart_hist)//trim(site)//'-'//'rstTimeVar'//'-'//trim(cdate)
      print*,trim(fhistTimeVar)
      OPEN(unit=lhistTimeVar,file=trim(fhistTimeVar),status='unknown',&
                             form='unformatted',action='write')
      WRITE(lhistTimeVar) id, & !
          ! Time-varying state variables which reaquired by restart run
            z_sno   , & !     node depth [m]
            dz_sno  , & !     interface depth [m]
            t_soisno   , & !  soil temperature [K]
            wliq_soisno, & !  liquid water in layers [kg/m2]
            wice_soisno, & !  ice lens in layers [kg/m2]
            t_grnd  , & !     ground surface temperature [K]
            tlsun   , & !     sunlit leaf temperature [K]
            tlsha   , & !     shaded leaf temperature [K]
            ldew    , & !     depth of water on foliage [mm]
            sag     , & !     non dimensional snow age [-]
            scv     , & !     snow cover, water equivalent [mm]
            snowdp  , & !     snow depth [meter]
            fveg    , & !     fraction of vegetation cover
            fsno    , & !     fraction of snow cover on ground
            sigf    , & !     fraction of veg cover, excluding snow-covered veg [-]
            green   , & !     leaf greenness
            lai     , & !     leaf area index
            sai     , & !     stem area index
            coszen  , & !     cosine of solar zenith angle
            albg    , & !     albedo, ground [-]
            albv    , & !     albedo, vegetation [-]
            alb     , & !     averaged albedo [-]
            ssun    , & !     sunlit canopy absorption for solar radiation (0-1)
            ssha    , & !     shaded canopy absorption for solar radiation (0-1)
            thermk  , & !     canopy gap fraction for tir radiation
            extkb   , & !     (k, g(mu)/mu) direct solar extinction coefficient
            extkd   , & !     diffuse and scattered diffuse PAR extinction coefficient
            zwt     , & !     the depth to water table [m]
            wa      , & !     water storage in aquifer [mm]
 
            t_lake       , & !
            lake_icefrac , & !

          ! Additional variables required by reginal model (such as WRF & RSM) 
            trad    , & !     radiative temperature of surface [K]
            tref    , & !     2 m height air temperature [kelvin]
            qref    , & !     2 m height air specific humidity
            rst     , & !     canopy stomatal resistance (s/m)
            emis    , & !     averaged bulk surface emissivity
            z0ma    , & !     effective roughness [m]
            zol     , & !     dimensionless height (z/L) used in Monin-Obukhov theory
            rib     , & !     bulk Richardson number in surface layer
            ustar   , & !     u* in similarity theory [m/s]
            qstar   , & !     q* in similarity theory [kg/kg]
            tstar   , & !     t* in similarity theory [K]
            fm      , & !     integral of profile function for momentum
            fh      , & !     integral of profile function for heat
            fq          !     integral of profile function for moisture

      CLOSE(lhistTimeVar)

  END SUBROUTINE WRITE_TimeVariables

  SUBROUTINE deallocate_TimeVariables
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------
      deallocate (z_sno    )
      deallocate (dz_sno   )
      deallocate (t_soisno    )
      deallocate (wliq_soisno )
      deallocate (wice_soisno )
      deallocate (h2osoi )
      deallocate (rstfac )
      deallocate (t_grnd )
      deallocate (tlsun  )
      deallocate (tlsha  )
      deallocate (ldew   )
      deallocate (sag    )
      deallocate (scv    )
      deallocate (snowdp )
      deallocate (fveg   )
      deallocate (fsno   )
      deallocate (sigf   )
      deallocate (green  )
      deallocate (lai    )
      deallocate (laisun )
      deallocate (laisha )
      deallocate (sai    )
      deallocate (coszen )
      deallocate (albg   )
      deallocate (albv   )
      deallocate (alb    )
      deallocate (ssun   )
      deallocate (ssha   )
      deallocate (thermk )
      deallocate (extkb  )
      deallocate (extkd  )
      deallocate (zwt    )
      deallocate (wa     )
      deallocate (wat    )

      deallocate (t_lake )      ! new lake scheme
      deallocate (lake_icefrac) ! new lake scheme

      deallocate (trad   )
      deallocate (tref   )
      deallocate (qref   )
      deallocate (rst    )
      deallocate (emis   )
      deallocate (z0ma   )
      deallocate (zol    )
      deallocate (rib    )
      deallocate (ustar  )
      deallocate (qstar  )
      deallocate (tstar  )
      deallocate (fm     )
      deallocate (fh     )
      deallocate (fq     )
! -----------------------
  END SUBROUTINE deallocate_TimeVariables

END MODULE MOD_TimeVariables
! ------ EOP --------------
