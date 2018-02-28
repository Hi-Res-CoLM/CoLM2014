#include <define.h>

SUBROUTINE initialize (site,dir_model_landdata,dir_restart_hist,&
                       idate,greenwich,lon_points,lat_points,&
                       nl_soil,maxsnl,nl_lake,numpatch,&
                       tg_xy,albvb_xy,albvd_xy,albnb_xy,albnd_xy,&
                       trad_xy,rib_xy,fm_xy,fh_xy,fq_xy)

! ======================================================================
! initialization routine for land surface model.
!
! Created by Yongjiu Dai, 09/15/1999
! Revised by Yongjiu Dai, 08/30/2002
! Revised by Yongjiu Dai, 03/2014
!             
! ======================================================================
   use precision
   use PhysicalConstants
   use MOD_TimeInvariants
   use MOD_TimeVariables
   use timemanager

   IMPLICIT NONE

! ----------------------------------------------------------------------
#if(defined USGS_CLASSIFICATION)
   integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
#endif
#if(defined IGBP_CLASSIFICATION)
   integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
#endif
   character(LEN=256), INTENT(in) :: site           ! site name
   character(LEN=256), INTENT(in) :: dir_model_landdata
   character(LEN=256), INTENT(in) :: dir_restart_hist
   logical, INTENT(in)    :: greenwich  ! true: greenwich time, false: local time
   integer, INTENT(in)    :: lon_points ! number of longitude points on model grid
   integer, INTENT(in)    :: lat_points ! number of latitude points on model grid
   integer, INTENT(in)    :: nl_soil    ! number of soil layers
   integer, INTENT(in)    :: maxsnl     ! max number of snow layers
   integer, INTENT(in)    :: nl_lake    ! number of land water bodies' layers
   integer, INTENT(out)   :: numpatch   ! total number of patches of grids
   integer, INTENT(inout) :: idate(3)   ! year, julian day, seconds of the starting time

! required by atmospheric models initialization (such as GRAPES, RSM, ...)
   real(r8), intent(out) :: tg_xy   (lon_points,lat_points) ! 
   real(r8), intent(out) :: albvb_xy(lon_points,lat_points) ! 
   real(r8), intent(out) :: albvd_xy(lon_points,lat_points) ! 
   real(r8), intent(out) :: albnb_xy(lon_points,lat_points) ! 
   real(r8), intent(out) :: albnd_xy(lon_points,lat_points) ! 
   real(r8), intent(out) :: trad_xy (lon_points,lat_points) ! 
   real(r8), intent(out) :: rib_xy  (lon_points,lat_points) ! 
   real(r8), intent(out) :: fm_xy   (lon_points,lat_points) ! 
   real(r8), intent(out) :: fh_xy   (lon_points,lat_points) ! 
   real(r8), intent(out) :: fq_xy   (lon_points,lat_points) ! 

! ------------------------ local variables -----------------------------
! surface classification and soil information

  real(r8) latixy  (lon_points,lat_points)          ! latitude in radians
  real(r8) longxy  (lon_points,lat_points)          ! longitude in radians
  real(r8) latdeg  (lat_points)                     ! latitude in degree
  real(r8) londeg  (lon_points)                     ! longitude in degree
  real(r8) area_gridcells (lon_points,lat_points)   ! area of gridcells (km^2)

  real(r8), allocatable :: fraction_patches(:,:,:) ! fraction of the patch of landtypes in gridcells

#if(defined SOILINI)
  integer :: lusoil
  integer :: nl_soil_ini
  real(r8), allocatable :: snow_d_grid(:,:)
  real(r8), allocatable :: snow_d(:)

  real(r8), allocatable :: soil_z_grid(:)
  real(r8), allocatable :: soil_t_grid(:,:,:)
  real(r8), allocatable :: soil_w_grid(:,:,:)
  real(r8), allocatable :: soil_z(:)
  real(r8), allocatable :: soil_t(:,:)
  real(r8), allocatable :: soil_w(:,:)
#endif

! CLM soil layer thickiness and depths
  real(r8), allocatable :: zsoi(:)     ! soil layer depth [m]
  real(r8), allocatable :: dzsoi(:)    ! soil node thickness [m]
  real(r8), allocatable :: zsoih(:)    ! interface level below a zsoi level [m]
  real(r8), allocatable :: z_soisno (:,:)
  real(r8), allocatable :: dz_soisno(:,:)

  real(r8) :: calday                   ! Julian cal day (1.xx to 365.xx)
  integer  :: year, jday, msec               ! Julian day and seconds
  real(r8) :: pi                       ! pie
  integer  :: i,j,k,l,m,npatch,np,nsl  ! indices
  integer  :: numpatch_lat(lat_points) ! number of patches of grids at lon. strip

  character(LEN=256) :: c
  character(LEN=255) :: cdate          ! character for date
  character(len=256) :: lndname
  integer iunit
  integer Julian_8day

  real(r8), external :: orb_coszen     ! cosine of the solar zenith angle


! ----------------------------------------------------------------------
! [1] READ IN LAND INFORMATION
! read time-invariant boundary data on [lon_points] x [lat_points] grid.
! ----------------------------------------------------------------------
! Read in the coordinate of the center of the model grids and area of grid cells
      iunit = 100
      lndname = trim(dir_model_landdata)//'model_lonlat_gridcell.bin'
      print*,trim(lndname)
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      READ(iunit) latixy
      READ(iunit) longxy
      READ(iunit) area_gridcells
      close(iunit)

      ! get grid latitudes and longitudes
      latdeg = latixy(1,:)
      londeg = longxy(:,1)
      
      ! convert latitudes and longitudes from degress to radians
      pi = 4.*atan(1.)        
      latixy(:,:) = latixy(:,:)*pi/180. 
      longxy(:,:) = longxy(:,:)*pi/180. 

! Read in the patch fraction of the lantypes of the gridcells
      allocate (fraction_patches(0:N_land_classification,1:lon_points,1:lat_points))
      lndname = trim(dir_model_landdata)//'model_landtypes.bin'
      print*,trim(lndname)
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      READ(iunit,err=100) fraction_patches
      close(iunit)
      print*,'fraction   =', minval(fraction_patches, mask = fraction_patches .gt. -1.0e30), &
                             maxval(fraction_patches ,mask = fraction_patches .gt. -1.0e30)

! ----------------------------------------------------------------------
! [2] MAPPING and ALLOCATE
! Build 1d subgrid patch <-> 2d grid mapping indices and weights
! 
! Build mapping indices and weights: [lon_points]x[lat_points] 2d grid <->
! <-> [numpatch] vector of subgrid patches. 
! The land surface model works by gathering all the land points on a
! [lon_points]x[lat_points] grid into a vector, and then expanded into 
! a vector of [numpatch] subgrid patches, allowing
! for up to [maxpatch=N_land_classification + 1] subgrid patches per land point. 
! [ixy], [jxy], [patch], and [land] are indices for the mapping: 
! [lon_points]x[lat_points] grid <-> [numpatch] vector of subgrid points. 
!
!-----------------------------------------------------------------------
! Find total number of patches [numpatch] allowing for multiple subgrid 
! patches in a grid cell.
! --------------------------------------------------------------------
      npatch = 0
      numpatch_lat(:) = 0

      do j = 1, lat_points
         do i = 1, lon_points
#if(defined LANDONLY)
            do np = 1, N_land_classification
               if(fraction_patches(np,i,j)> 0.)then
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               endif
            end do
#elif(defined LAND_SEA || defined USE_POINT_DATA)
            do np = 0, N_land_classification
               if(fraction_patches(np,i,j)> 0.)then
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               endif
            end do
#endif
         end do
      end do
      numpatch = npatch
      if(numpatch.ne.sum(numpatch_lat))then
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         call abort
      endif
      write(6,*) 'Total land patches = ', numpatch


! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------------------------

      CALL allocate_TimeInvariants (nl_soil,nl_lake,numpatch,lon_points,lat_points)

      CALL allocate_TimeVariables (nl_soil,nl_lake,maxsnl,numpatch)


! set the lat/lon values
      lats = latdeg
      lons = londeg

! --------------------------------------------------------------------
! Build 1d land vector and 1d patch vector mapping components
! --------------------------------------------------------------------

! Determine land vector and patch vector mapping components

      npatch = 0
      wtxy_patch(:)  = 0.
      l = 0; m = 0
      spatch_xy(:,:) = -1
      epatch_xy(:,:) = -1

      do j = 1, lat_points
         do i = 1, lon_points
#if(defined LANDONLY)
            do np = 1, N_land_classification
#elif(defined LAND_SEA || defined USE_POINT_DATA)
            do np = 0, N_land_classification
#endif
               if(fraction_patches(np,i,j)>0.)then                

                  npatch             = npatch+1                      
                  ixy_patch(npatch)  = i            ! patch longitude index
                  jxy_patch(npatch)  = j            ! patch latitude index
                  mxy_patch(npatch)  = np           ! index of land cover type of the patches at the fraction > 0
                  dlat(npatch)       = latixy(i,j)  ! latitude in radians
                  dlon(npatch)       = longxy(i,j)  ! longitude in radians

                  wtxy_patch(npatch) = fraction_patches(np,i,j) ! patch weight
                  area_patch(npatch) = area_gridcells(i,j)      ! grid cell area
                  epatch_xy(i,j)     = npatch

                  if (l.ne.i .OR. m.ne.j) then
                     l = i; m = j; spatch_xy(i,j) = npatch
                  end if

               end if
            end do
         end do
      end do

      if(numpatch.ne.npatch)then
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         call abort
      endif

! ---------------------------------------------------------------
! [3] INITIALIZE TIME INVARIANT VARIABLES
! ---------------------------------------------------------------
! 3.1 Define the soil and lake layers's thickness 
! ...............................................................
      allocate ( zsoi (1:nl_soil) )
      allocate ( dzsoi(1:nl_soil) )
      allocate ( zsoih(0:nl_soil) )
      allocate ( z_soisno (maxsnl+1:nl_soil,numpatch) )
      allocate ( dz_soisno(maxsnl+1:nl_soil,numpatch) )

      do nsl = 1, nl_soil
         zsoi(nsl) = 0.025*(exp(0.5*(nsl-0.5))-1.)  ! node depths
      end do

      dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))              ! =zsoih(1)
      dzsoi(nl_soil) = zsoi(nl_soil)-zsoi(nl_soil-1)
      do nsl = 2, nl_soil-1
         dzsoi(nsl) = 0.5*(zsoi(nsl+1)-zsoi(nsl-1)) ! thickness b/n two interfaces
      end do

      zsoih(0) = 0.
      zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil)
      do nsl = 1, nl_soil-1
         zsoih(nsl) = 0.5*(zsoi(nsl)+zsoi(nsl+1))   ! interface depths
      enddo

      do npatch = 1, numpatch 
         z_soi (1:nl_soil,npatch) = zsoi(1:nl_soil)
         dz_soi(1:nl_soil,npatch) = dzsoi(1:nl_soil)
      enddo

    ! ------------------------------------------
    ! Lake depth and layers' thickness
    ! ------------------------------------------
      CALL lakedepth_readin (lon_points,lat_points,nl_lake,numpatch,dir_model_landdata)

! ...............................................................
! 3.2 Read in the soil parameters of the patches of the gridcells
! ...............................................................

      CALL soil_parameters_readin (lon_points,lat_points,nl_soil,numpatch,dir_model_landdata)

! ...............................................................
! 3.3 Plant time-invariant variables (based on the look-up tables) 
! ...............................................................
      do i = 1, numpatch
      CALL IniTimeConst(nl_soil,zsoih,mxy_patch(i),itypwat(i),&
           z0m(i),displa(i),sqrtdi(i),effcon(i),vmax25(i),slti(i),hlti(i),&
           shti(i),hhti(i),trda(i),trdm(i),trop(i),gradm(i),binter(i),extkn(i),&
           chil(i),ref(1:,1:,i),tran(1:,1:,i),rootfr(1:,i))
      enddo

! ................................
! 3.4 Initialize TUNABLE constants
! ................................
      zlnd   = 0.01    !Roughness length for soil [m]
      zsno   = 0.0024  !Roughness length for snow [m]
      csoilc = 0.004   !Drag coefficient for soil under canopy [-]
      dewmx  = 0.1     !maximum dew
      wtfact = 0.38    !Maximum saturated fraction (global mean; see Niu et al., 2005)
      capr   = 0.34    !Tuning factor to turn first layer T into surface T
      cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
      ssi    = 0.033   !Irreducible water saturation of snow
      wimp   = 0.05    !Water impremeable if porosity less than wimp
      pondmx = 10.0    !Ponding depth (mm)
      smpmax = -1.5e5  !Wilting point potential in mm
      smpmin = -1.e8   !Restriction for min of soil poten. (mm)
      trsmx0 = 2.e-4   !Max transpiration for moist soil+100% veg. [mm/s]
      tcrit  = 2.5     !critical temp. to determine rain or snow

! ...............................................
! 3.5 Write out as a restart file [histTimeConst]
! ...............................................

      CALL WRITE_TimeInvariants(dir_restart_hist,site)

      write (6,*)
      write (6,*) ('Successfully to Initialize the Land Time-Invariants')

! ----------------------------------------------------------------------
! [4] INITIALIZE TIME-VARYING VARIABLES 
! as subgrid vectors of length [numpatch]
! initial run: create the time-varying variables based on :
!              i) observation (NOT CODING CURRENTLY), or
!             ii) some already-known information (NO CODING CURRENTLY), or
!            iii) arbitrarily 
! continuation run: time-varying data read in from restart file 
! ----------------------------------------------------------------------

! 4.1 current time of model run
! ............................

      call initimetype(greenwich)

      if(.not. greenwich)then
         print *, ".........greenwich false", longxy(1,1)
      end if

      year = idate(1)
      jday = idate(2)
      msec = idate(3)

! ................................
! 4.2 cosine of solar zenith angle 
! ................................
      calday = calendarday(idate, lons(1))
      do i = 1, numpatch
         coszen(i) = orb_coszen(calday,dlon(i),dlat(i))
      enddo

#if(defined SOILINI)
! ...........................................
!4.3 READ in or GUSSES land state information
! ...........................................
!!! PLEASE CHANGE
!!! PLEASE CHANGE
!!! PLEASE CHANGE the root of directory when the soil T and W are ready !!!
     lusoil = 100
     fsoildat_name = trim(fsoildat)//'-'//trim(cdate)
     OPEN(lusoil,file=fsoildat_name,status='old',form='unformatted',action='read')

     read(lusoil) nl_soil_ini

     allocate (snow_d_grid(lon_points,lat_points))
     allocate (snow_d(numpatch))

     allocate (soil_z_grid(nl_soil_ini))
     allocate (soil_t_grid(lon_points,lat_points,nl_soil_ini))
     allocate (soil_w_grid(lon_points,lat_points,nl_soil_ini))

     allocate (soil_z(nl_soil_ini))
     allocate (soil_t(nl_soil_ini,numpatch))
     allocate (soil_w(nl_soil_ini,numpatch))

     read(lusoil) soil_z_grid  ! soil layer node depth (m)
     read(lusoil) soil_t_grid  ! soil layer temperature (K)
     read(lusoil) soil_w_grid  ! soil layer wetness (-)
     read(lusoil) snow_d_grid  ! snow depth (m)

     close (lusoil)

     soil_z(:) = soil_z_grid(:) 
     do npatch = 1, numpatch
        i = ixy_patch(npatch)
        j = jxy_patch(npatch)

        snow_d(npatch) = snow_d_grid(i,j)
        do l = 1, nl_soil_ini
           soil_t(l,npatch) = soil_t_grid(i,j,l)
           soil_w(l,npatch) = soil_w_grid(i,j,l)
        enddo
     end do
#endif

! ...................
! 4.4 LEAF area index
! ...................
#if(defined DYN_PHENOLOGY)
    ! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
      lai(:)=0.0; sai(:)=0.0; green(:)=0.0; fveg(:)=0.0
      do i = 1, numpatch
#if(defined SOILINI)
         do l = 1, nl_soil
            t_soisno(l,i) = soil_t(nl_soil_ini,i)
         enddo
#else
         t_soisno(1:,i) = 283.
#endif
       ! Call Ecological Model() 
         if(mxy_patch(i)>0)then
            call lai_empirical(mxy_patch(i),nl_soil,rootfr(1:,i),&
                               t_soisno(1:,i),lai(i),sai(i),fveg(i),green(i))

         endif
      enddo
#else
    ! READ in Leaf area index and stem area index
      Julian_8day = int(calendarday(idate)-1)/8*8 + 1
      CALL LAI_readin (lon_points,lat_points,&
                       Julian_8day,numpatch,dir_model_landdata)
#endif

! ..............................................................................
! 4.5 initialize time-varying variables, as subgrid vectors of length [numpatch]
! ..............................................................................
      do i = 1, numpatch
         z_soisno (1:nl_soil ,i) = z_soi (1:nl_soil ,i)
         dz_soisno(1:nl_soil ,i) = dz_soi(1:nl_soil ,i)
      enddo

      do i = 1, numpatch
      CALL iniTimeVar(nl_soil,maxsnl,itypwat(i)&
          ,porsl(1:,i)&
          ,soil_s_v_alb(i),soil_d_v_alb(i),soil_s_n_alb(i),soil_d_n_alb(i)&
          ,z0m(i),zlnd,chil(i),ref(1:,1:,i),tran(1:,1:,i)&
          ,z_soisno(maxsnl+1:,i),dz_soisno(maxsnl+1:,i)&
          ,t_soisno(maxsnl+1:,i),wliq_soisno(maxsnl+1:,i),wice_soisno(maxsnl+1:,i)&
          ,zwt(i),wa(i)&
          ,t_grnd(i),tlsun(i),tlsha(i),ldew(i),sag(i),scv(i)&
          ,snowdp(i),fveg(i),fsno(i),sigf(i),green(i),lai(i),sai(i),coszen(i)&
          ,albg(1:,1:,i),albv(1:,1:,i),alb(1:,1:,i),ssun(1:,1:,i),ssha(1:,1:,i)&
          ,thermk(i),extkb(i),extkd(i)&
          ,trad(i),tref(i),qref(i),rst(i),emis(i),z0ma(i),zol(i),rib(i)&
          ,ustar(i),qstar(i),tstar(i),fm(i),fh(i),fq(i)&
#if(defined SOILINI)
          ,nl_soil_ini,soil_z,soil_t(1:,i),soil_w(1:,i),snow_d(i))
#else
          )
#endif
      enddo

      do i = 1, numpatch
         z_sno (maxsnl+1:0,i) = z_soisno (maxsnl+1:0,i)
         dz_sno(maxsnl+1:0,i) = dz_soisno(maxsnl+1:0,i)
      end do

    ! ------------------------------------------
    ! PLEASE  
    ! PLEASE UPDATE
    ! PLEASE UPDATE when have the observed lake status
      t_lake      (:,:) = 285.
      lake_icefrac(:,:) = 0.
    ! ------------------------------------------

! ...............................................................
! 4.6 Write out the model variables for restart run [histTimeVar]
! ...............................................................

      CALL WRITE_TimeVariables (idate,dir_restart_hist,site)

      write (6,*)
      write (6,*) ('Successfully to Initialize the Land Time-Vraying Variables')


! ----------------------------------------------------------------------
    ! average subgrid albedos, srf temperature, etc. for atmospheric model
! ----------------------------------------------------------------------
      tg_xy   (:,:) = 0.0
      albvb_xy(:,:) = 0.0
      albvd_xy(:,:) = 0.0
      albnb_xy(:,:) = 0.0
      albnd_xy(:,:) = 0.0
      trad_xy (:,:) = 0.0
      rib_xy  (:,:) = 0.0
      fm_xy   (:,:) = 0.0
      fh_xy   (:,:) = 0.0
      fq_xy   (:,:) = 0.0

      do npatch = 1, numpatch
         i = ixy_patch(npatch)
         j = jxy_patch(npatch)
            tg_xy(i,j) =    tg_xy(i,j) + wtxy_patch(npatch)*t_grnd (npatch)   ! 
         albvb_xy(i,j) = albvb_xy(i,j) + wtxy_patch(npatch)*alb(1,1,npatch)   ! (v=visble, b=direct beam)
         albvd_xy(i,j) = albvd_xy(i,j) + wtxy_patch(npatch)*alb(1,2,npatch)   ! (v=visble, d=diffuse)
         albnb_xy(i,j) = albnb_xy(i,j) + wtxy_patch(npatch)*alb(2,1,npatch)   ! (n=near infrared, b=direct beam)
         albnd_xy(i,j) = albnd_xy(i,j) + wtxy_patch(npatch)*alb(2,2,npatch)   ! (n=near infrared, d=diffuse)
         trad_xy (i,j) =  trad_xy(i,j) + wtxy_patch(npatch)*t_grnd (npatch)   !

         rib_xy(i,j) = -0.1       !
         fm_xy (i,j) = alog(30.)  !
         fh_xy (i,j) = alog(30.)  !
         fq_xy (i,j) = alog(30.)  !
      enddo

    ! --------------------------------------------------
    ! Deallocates memory for CLM 1d [numpatch] variables
    ! --------------------------------------------------

      CALL deallocate_TimeInvariants

      CALL deallocate_TimeVariables

      deallocate (fraction_patches)
      deallocate (zsoi   )
      deallocate (dzsoi  )
      deallocate (zsoih  )
      deallocate (z_soisno  )
      deallocate (dz_soisno )

#if(defined SOILINI)
      deallocate (snow_d_grid)
      deallocate (snow_d)
      deallocate (soil_z_grid)
      deallocate (soil_t_grid)
      deallocate (soil_w_grid)
      deallocate (soil_z)
      deallocate (soil_t)
      deallocate (soil_w)
#endif


      go to 1000
100   print 101,lndname
101   format(' error occured on file: ',a50)
1000  continue


END SUBROUTINE initialize
! --------------------------------------------------
! EOP
