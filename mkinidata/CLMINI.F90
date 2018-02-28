#include <define.h>

PROGRAM CLMINI
! ======================================================================
! Initialization of Land Characteristic Parameters and Initial State Variables
!
! Reference:
!     [1] Dai et al., 2003: The Common Land Model (CoLM).
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. Journal of Climate
!     [3] Dai et al., 2014: The Terrestrial Modeling System (TMS).
!
!     Created by Yongjiu Dai Februay 2004
!     Revised by Yongjiu Dai Februay 2014
! ======================================================================

      use precision
      implicit none

! ----------------local variables ---------------------------------
      character(LEN=256) :: site ! site name
      character(LEN=256) :: dir_model_landdata
      character(LEN=256) :: dir_restart_hist
      character(LEN=256) :: dir_infolist
      integer :: s_year      ! starting date for run in year
      integer :: s_julian    ! starting date for run in julian day
      integer :: s_seconds   ! starting time of day for run in seconds
      integer :: idate(3)    ! starting date
      logical :: greenwich   ! true: greenwich time, false: local time

      integer :: lon_points  ! number of longitude points on model grid
      integer :: lat_points  ! number of latitude points on model grid
      integer :: numpatch    ! total number of patches of grids

      character(LEN=256) :: finfolist      ! file name of run information

#if(defined UGGS_CLASSIFICATION)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
#endif
#if(defined IGBP_CLASSIFICATION)
      integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
#endif
      integer, parameter :: nl_soil = 10
      integer, parameter :: maxsnl = -5
      integer, parameter :: nl_lake = 10

!  Required by atmospheric models's initialization (such as GRAPES/WRF/RSM/EMSs) 
      real(r8), allocatable :: tg_xy   (:,:)
      real(r8), allocatable :: albvb_xy(:,:)
      real(r8), allocatable :: albvd_xy(:,:)
      real(r8), allocatable :: albnb_xy(:,:)
      real(r8), allocatable :: albnd_xy(:,:)
      real(r8), allocatable :: trad_xy (:,:)
      real(r8), allocatable :: rib_xy  (:,:)
      real(r8), allocatable :: fm_xy   (:,:)
      real(r8), allocatable :: fh_xy   (:,:)
      real(r8), allocatable :: fq_xy   (:,:)

      namelist /clminiexp/ site,dir_model_landdata,&
                           dir_restart_hist,dir_infolist,&
                           lon_points,lat_points,greenwich,&
                           s_year,s_julian,s_seconds
! ----------------------------------------------------------------------
      read (5,clminiexp)

      idate(1) = s_year; idate(2) = s_julian; idate(3) = s_seconds

      allocate ( tg_xy   (lon_points,lat_points) )
      allocate ( albvb_xy(lon_points,lat_points) )
      allocate ( albvd_xy(lon_points,lat_points) )
      allocate ( albnb_xy(lon_points,lat_points) )
      allocate ( albnd_xy(lon_points,lat_points) )
      allocate ( trad_xy (lon_points,lat_points) )
      allocate ( rib_xy  (lon_points,lat_points) )
      allocate ( fm_xy   (lon_points,lat_points) )
      allocate ( fh_xy   (lon_points,lat_points) )
      allocate ( fq_xy   (lon_points,lat_points) )


      CALL initialize (site,dir_model_landdata,dir_restart_hist,&
                       idate,greenwich,lon_points,lat_points,&
                       nl_soil,maxsnl,nl_lake,numpatch,&
                       tg_xy,albvb_xy,albvd_xy,albnb_xy,albnd_xy,&
                       trad_xy,rib_xy,fm_xy,fh_xy,fq_xy)

      finfolist = trim(dir_infolist)//'clmini.infolist'
      OPEN(100,file=trim(finfolist),form='formatted')
      write(100,*) 'numpatch  = ', numpatch   !1
      if ( greenwich ) then
      write(100,*) 'greenwich =       .true.' !2
      else
      write(100,*) 'greenwich =      .false.' !2
      end if
      write(100,*) 's_year    = ', idate(1)   !3
      write(100,*) 's_julian  = ', idate(2)   !4
      write(100,*) 's_seconds = ', idate(3)   !5
      write(100,*) '/'
      CLOSE(100)

      deallocate ( tg_xy    )
      deallocate ( albvb_xy )
      deallocate ( albvd_xy )
      deallocate ( albnb_xy )
      deallocate ( albnd_xy )
      deallocate ( trad_xy  )
      deallocate ( rib_xy   )
      deallocate ( fm_xy    )
      deallocate ( fh_xy    )
      deallocate ( fq_xy    )

      write(6,*) 'CLM Initialization Execution Completed'

END PROGRAM CLMINI
! ----------------------------------------------------------------------
! EOP
