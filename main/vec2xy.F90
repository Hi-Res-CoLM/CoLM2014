
SUBROUTINE vec2xy (lon_points,lat_points,numpatch,nl_soil,maxsnl,&
                   nl_lake,nac,nac_ln,a_rnof)
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------

use precision
use PhysicalConstants, only: vonkar, stefnc, cpair, rgas, grav
use MOD_TimeInvariants
use MOD_TimeVariables
use MOD_1D_Forcing
use MOD_2D_Forcing
use MOD_1D_Fluxes
use MOD_2D_Fluxes
use FRICTION_VELOCITY

IMPLICIT NONE

integer, INTENT(in) :: lon_points
integer, INTENT(in) :: lat_points
integer, INTENT(in) :: numpatch
integer, INTENT(in) :: nl_soil
integer, INTENT(in) :: maxsnl
integer, INTENT(in) :: nl_lake
integer, INTENT(inout) :: nac
integer, INTENT(inout) :: nac_ln(lon_points,lat_points)
      
!---------------------------------------------------------------------
real(r8) a_xy_us  (lon_points,lat_points)  ! wind in eastward direction [m/s]
real(r8) a_xy_vs  (lon_points,lat_points)  ! wind in northward direction [m/s]
real(r8) a_xy_t   (lon_points,lat_points)  ! temperature at reference height [kelvin]
real(r8) a_xy_q   (lon_points,lat_points)  ! specific humidity at reference height [kg/kg]
real(r8) a_xy_prc (lon_points,lat_points)  ! convective precipitation [mm/s]
real(r8) a_xy_prl (lon_points,lat_points)  ! large scale precipitation [mm/s]
real(r8) a_xy_pbot(lon_points,lat_points)  ! atmospheric pressure at the surface [pa]
real(r8) a_xy_frl (lon_points,lat_points)  ! atmospheric infrared (longwave) radiation [W/m2]
real(r8) a_xy_solarin(lon_points,lat_points)  ! downward solar radiation at surface [W/m2]
real(r8) a_xy_rain(lon_points,lat_points)  ! rain [mm/s]
real(r8) a_xy_snow(lon_points,lat_points)  ! snow [mm/s]

!---------------------------------------------------------------------
real(r8) a_taux   (lon_points,lat_points)  ! wind stress: E-W [kg/m/s2]
real(r8) a_tauy   (lon_points,lat_points)  ! wind stress: N-S [kg/m/s2]
real(r8) a_fsena  (lon_points,lat_points)  ! sensible heat from canopy height to atmosphere [W/m2]
real(r8) a_lfevpa (lon_points,lat_points)  ! latent heat flux from canopy height to atmosphere [W/m2]
real(r8) a_fevpa  (lon_points,lat_points)  ! evapotranspiration from canopy to atmosphere [mm/s]
real(r8) a_fsenl  (lon_points,lat_points)  ! sensible heat from leaves [W/m2]
real(r8) a_fevpl  (lon_points,lat_points)  ! evaporation+transpiration from leaves [mm/s]
real(r8) a_etr    (lon_points,lat_points)  ! transpiration rate [mm/s]
real(r8) a_fseng  (lon_points,lat_points)  ! sensible heat flux from ground [W/m2]
real(r8) a_fevpg  (lon_points,lat_points)  ! evaporation heat flux from ground [mm/s]
real(r8) a_fgrnd  (lon_points,lat_points)  ! ground heat flux [W/m2]
real(r8) a_sabvsun(lon_points,lat_points)  ! solar absorbed by sunlit canopy [W/m2]
real(r8) a_sabvsha(lon_points,lat_points)  ! solar absorbed by shaded [W/m2]
real(r8) a_sabg   (lon_points,lat_points)  ! solar absorbed by ground  [W/m2]
real(r8) a_olrg   (lon_points,lat_points)  ! outgoing long-wave radiation from ground+canopy [W/m2]
real(r8) a_rnet   (lon_points,lat_points)  ! net radiation [W/m2]
real(r8) a_xerr   (lon_points,lat_points)  ! the error of water banace [mm/s]
real(r8) a_zerr   (lon_points,lat_points)  ! the error of energy balance [W/m2]
real(r8) a_rsur   (lon_points,lat_points)  ! surface runoff [mm/s]
real(r8) a_rnof   (lon_points,lat_points)  ! total runoff [mm/s]
real(r8) a_qintr  (lon_points,lat_points)  ! interception [mm/s]
real(r8) a_qinfl  (lon_points,lat_points)  ! inflitration [mm/s]
real(r8) a_qdrip  (lon_points,lat_points)  ! throughfall [mm/s]

real(r8) a_assim  (lon_points,lat_points)  ! canopy assimilation rate [mol m-2 s-1]
real(r8) a_respc  (lon_points,lat_points)  ! respiration (plant+soil) [mol m-2 s-1]
real(r8) a_qcharge(lon_points,lat_points)  ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
real(r8) a_t_grnd (lon_points,lat_points)  ! ground surface temperature [K]
real(r8) a_tlsun  (lon_points,lat_points)  ! sunlit leaf temperature [K]
real(r8) a_tlsha  (lon_points,lat_points)  ! shaded leaf temperature [K]
real(r8) a_ldew   (lon_points,lat_points)  ! depth of water on foliage [mm]
real(r8) a_scv    (lon_points,lat_points)  ! snow cover, water equivalent [mm]
real(r8) a_snowdp (lon_points,lat_points)  ! snow depth [meter]
real(r8) a_fsno   (lon_points,lat_points)  ! fraction of snow cover on ground
real(r8) a_sigf   (lon_points,lat_points)  ! fraction of veg cover, excluding snow-covered veg [-]
real(r8) a_green  (lon_points,lat_points)  ! leaf greenness
real(r8) a_lai    (lon_points,lat_points)  ! leaf area index
real(r8) a_laisun (lon_points,lat_points)  ! sunlit leaf area index
real(r8) a_laisha (lon_points,lat_points)  ! shaded leaf area index
real(r8) a_sai    (lon_points,lat_points)  ! stem area index
real(r8) a_alb(2,2,lon_points,lat_points)  ! averaged albedo [visible, direct; direct, diffuse]
real(r8) a_emis   (lon_points,lat_points)  ! averaged bulk surface emissivity
real(r8) a_z0ma   (lon_points,lat_points)  ! effective roughness [m]
real(r8) a_trad   (lon_points,lat_points)  ! radiative temperature of surface [K]
real(r8) a_tref   (lon_points,lat_points)  ! 2 m height air temperature [kelvin]
real(r8) a_qref   (lon_points,lat_points)  ! 2 m height air specific humidity [kg/kg]

!---------------------------------------------------------------------
real(r8) a_t_soisno   (maxsnl+1:nl_soil,lon_points,lat_points)  ! soil temperature [K]
real(r8) a_wliq_soisno(maxsnl+1:nl_soil,lon_points,lat_points)  ! liquid water in soil layers [kg/m2]
real(r8) a_wice_soisno(maxsnl+1:nl_soil,lon_points,lat_points)  ! ice lens in soil layers [kg/m2]
real(r8) a_h2osoi            (1:nl_soil,lon_points,lat_points)  ! volumetric soil water in layers [m3/m3]
real(r8) a_rstfac                      (lon_points,lat_points)  ! factor of soil water stress 
real(r8) a_zwt                         (lon_points,lat_points)  ! the depth to water table [m]
real(r8) a_wa                          (lon_points,lat_points)  ! water storage in aquifer [mm]
real(r8) a_wat                         (lon_points,lat_points)  ! total water storage [mm]

real(r8) a_t_lake      (nl_lake,lon_points,lat_points) ! lake temperature [K]
real(r8) a_lake_icefrac(nl_lake,lon_points,lat_points) ! lake ice fraction cover [0-1]

!---------------------------------------------------------------------
real(r8) a_ustar  (lon_points,lat_points)  ! u* in similarity theory [m/s]
real(r8) a_tstar  (lon_points,lat_points)  ! t* in similarity theory [kg/kg]
real(r8) a_qstar  (lon_points,lat_points)  ! q* in similarity theory [kg/kg]
real(r8) a_zol    (lon_points,lat_points)  ! dimensionless height (z/L) used in Monin-Obukhov theory
real(r8) a_rib    (lon_points,lat_points)  ! bulk Richardson number in surface layer
real(r8) a_fm     (lon_points,lat_points)  ! integral of profile function for momentum
real(r8) a_fh     (lon_points,lat_points)  ! integral of profile function for heat
real(r8) a_fq     (lon_points,lat_points)  ! integral of profile function for moisture

real(r8) a_us10m  (lon_points,lat_points)  ! 10m u-velocity [m/s]
real(r8) a_vs10m  (lon_points,lat_points)  ! 10m v-velocity [m/s]
real(r8) a_fm10m  (lon_points,lat_points)  ! integral of profile function for momentum at 10m [-]

!---------------------------------------------------------------------
real(r8) a_sr     (lon_points,lat_points)  ! total reflected solar radiation (W/m2)
real(r8) a_solvd  (lon_points,lat_points)  ! incident direct beam vis solar radiation (W/m2)
real(r8) a_solvi  (lon_points,lat_points)  ! incident diffuse beam vis solar radiation (W/m2)
real(r8) a_solnd  (lon_points,lat_points)  ! incident direct beam nir solar radiation (W/m2)
real(r8) a_solni  (lon_points,lat_points)  ! incident diffuse beam nir solar radiation (W/m2)
real(r8) a_srvd   (lon_points,lat_points)  ! reflected direct beam vis solar radiation (W/m2)
real(r8) a_srvi   (lon_points,lat_points)  ! reflected diffuse beam vis solar radiation (W/m2)
real(r8) a_srnd   (lon_points,lat_points)  ! reflected direct beam nir solar radiation (W/m2)
real(r8) a_srni   (lon_points,lat_points)  ! reflected diffuse beam nir solar radiation (W/m2)
real(r8) a_solvdln(lon_points,lat_points) ! incident direct beam vis solar radiation at local noon (W/m2)
real(r8) a_solviln(lon_points,lat_points) ! incident diffuse beam vis solar radiation at local noon (W/m2)
real(r8) a_solndln(lon_points,lat_points) ! incident direct beam nir solar radiation at local noon (W/m2)
real(r8) a_solniln(lon_points,lat_points) ! incident diffuse beam nir solar radiation at local noon (W/m2)
real(r8) a_srvdln (lon_points,lat_points) ! reflected direct beam vis solar radiation at local noon (W/m2)
real(r8) a_srviln (lon_points,lat_points) ! reflected diffuse beam vis solar radiation at local noon (W/m2)
real(r8) a_srndln (lon_points,lat_points) ! reflected direct beam nir solar radiation at local noon (W/m2)
real(r8) a_srniln (lon_points,lat_points)  ! reflected diffuse beam nir solar radiation at local noon (W/m2)

!---------------------------------------------------------------------
! local variables

      integer  i,j,np,l
      real(r8) sumwt(lon_points,lat_points)
      real(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
      real(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf
      real(r8) obu,fh2m,fq2m
      real(r8) um,thvstar,beta,zii,wc,wc2

! ---------------------------------------------------
! Meteorological forcing
! ---------------------------------------------------
      a_xy_us     (:,:) = forc_xy_us     (:,:)
      a_xy_vs     (:,:) = forc_xy_vs     (:,:)
      a_xy_t      (:,:) = forc_xy_t      (:,:)
      a_xy_q      (:,:) = forc_xy_q      (:,:)
      a_xy_prc    (:,:) = forc_xy_prc    (:,:)
      a_xy_prl    (:,:) = forc_xy_prl    (:,:)
      a_xy_pbot   (:,:) = forc_xy_pbot   (:,:)
      a_xy_frl    (:,:) = forc_xy_frl    (:,:)

      a_xy_solarin(:,:) = forc_xy_sols (:,:) + forc_xy_soll (:,:) &
                        + forc_xy_solsd(:,:) + forc_xy_solld(:,:)

! ------------------------------------------------------------------------------------------
! Mapping the fluxes and state variables at patch [numpatch] to grid [lon_points,lat_points]
! ------------------------------------------------------------------------------------------
      sumwt    (:,:) = 0.
      a_taux   (:,:) = 0.
      a_tauy   (:,:) = 0.
      a_fsena  (:,:) = 0.
      a_lfevpa (:,:) = 0.
      a_fevpa  (:,:) = 0.
      a_fsenl  (:,:) = 0.
      a_fevpl  (:,:) = 0.
      a_etr    (:,:) = 0.
      a_fseng  (:,:) = 0.
      a_fevpg  (:,:) = 0.
      a_fgrnd  (:,:) = 0.
      a_sabvsun(:,:) = 0.
      a_sabvsha(:,:) = 0.
      a_sabg   (:,:) = 0.
      a_olrg   (:,:) = 0.
      a_rnet   (:,:) = 0.
      a_xerr   (:,:) = 0.
      a_zerr   (:,:) = 0.
      a_rsur   (:,:) = 0.
      a_rnof   (:,:) = 0.
      a_qintr  (:,:) = 0.
      a_qinfl  (:,:) = 0.
      a_qdrip  (:,:) = 0.
      a_wat    (:,:) = 0.
      a_assim  (:,:) = 0.
      a_respc  (:,:) = 0.

      a_qcharge(:,:) = 0.
      a_t_grnd (:,:) = 0.
      a_tlsun  (:,:) = 0.
      a_tlsha  (:,:) = 0.
      a_ldew   (:,:) = 0.
      a_scv    (:,:) = 0.
      a_snowdp (:,:) = 0.
      a_fsno   (:,:) = 0.
      a_sigf   (:,:) = 0.
      a_green  (:,:) = 0.
      a_lai    (:,:) = 0.
      a_laisun (:,:) = 0.
      a_laisha (:,:) = 0.
      a_sai    (:,:) = 0.
      a_alb(:,:,:,:) = 0.
      a_emis   (:,:) = 0.
      a_z0ma   (:,:) = 0.
      a_trad   (:,:) = 0.
      a_tref   (:,:) = 0.
      a_qref   (:,:) = 0.
      a_xy_rain(:,:) = 0.
      a_xy_snow(:,:) = 0.

      a_sr     (:,:) = spval
      a_solvd  (:,:) = spval
      a_solvi  (:,:) = spval
      a_solnd  (:,:) = spval
      a_solni  (:,:) = spval
      a_srvd   (:,:) = spval
      a_srvi   (:,:) = spval
      a_srnd   (:,:) = spval
      a_srni   (:,:) = spval
      a_solvdln(:,:) = spval
      a_solviln(:,:) = spval
      a_solndln(:,:) = spval
      a_solniln(:,:) = spval
      a_srvdln (:,:) = spval
      a_srviln (:,:) = spval
      a_srndln (:,:) = spval
      a_srniln (:,:) = spval

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         do i = 1, lon_points
          
            if (spatch_xy(i,j) .le. 0) cycle
            DO np = spatch_xy(i,j), epatch_xy(i,j)
               sumwt(i,j) = sumwt(i,j) + wtxy_patch(np)
! Fluxes
               a_taux   (i,j) = a_taux   (i,j) + wtxy_patch(np)*taux   (np)
               a_tauy   (i,j) = a_tauy   (i,j) + wtxy_patch(np)*tauy   (np)
               a_fsena  (i,j) = a_fsena  (i,j) + wtxy_patch(np)*fsena  (np)
               a_lfevpa (i,j) = a_lfevpa (i,j) + wtxy_patch(np)*lfevpa (np)
               a_fevpa  (i,j) = a_fevpa  (i,j) + wtxy_patch(np)*fevpa  (np)
               a_fsenl  (i,j) = a_fsenl  (i,j) + wtxy_patch(np)*fsenl  (np)
               a_fevpl  (i,j) = a_fevpl  (i,j) + wtxy_patch(np)*fevpl  (np)
               a_etr    (i,j) = a_etr    (i,j) + wtxy_patch(np)*etr    (np)
               a_fseng  (i,j) = a_fseng  (i,j) + wtxy_patch(np)*fseng  (np)
               a_fevpg  (i,j) = a_fevpg  (i,j) + wtxy_patch(np)*fevpg  (np)
               a_fgrnd  (i,j) = a_fgrnd  (i,j) + wtxy_patch(np)*fgrnd  (np)
               a_sabvsun(i,j) = a_sabvsun(i,j) + wtxy_patch(np)*sabvsun(np)
               a_sabvsha(i,j) = a_sabvsha(i,j) + wtxy_patch(np)*sabvsha(np)
               a_sabg   (i,j) = a_sabg   (i,j) + wtxy_patch(np)*sabg   (np)
               a_olrg   (i,j) = a_olrg   (i,j) + wtxy_patch(np)*olrg   (np)
               a_rnet   (i,j) = a_rnet   (i,j) + wtxy_patch(np)*(sabg(np)+sabvsun(np)+sabvsha(np)-olrg(np))
               a_xerr   (i,j) = a_xerr   (i,j) + wtxy_patch(np)*xerr   (np)
               a_zerr   (i,j) = a_zerr   (i,j) + wtxy_patch(np)*zerr   (np)
               a_rsur   (i,j) = a_rsur   (i,j) + wtxy_patch(np)*rsur   (np)
               a_rnof   (i,j) = a_rnof   (i,j) + wtxy_patch(np)*rnof   (np)
               a_qintr  (i,j) = a_qintr  (i,j) + wtxy_patch(np)*qintr  (np)
               a_qinfl  (i,j) = a_qinfl  (i,j) + wtxy_patch(np)*qinfl  (np)
               a_qdrip  (i,j) = a_qdrip  (i,j) + wtxy_patch(np)*qdrip  (np)
               a_wat    (i,j) = a_wat    (i,j) + wtxy_patch(np)*wat    (np)
               a_assim  (i,j) = a_assim  (i,j) + wtxy_patch(np)*assim  (np)
               a_respc  (i,j) = a_respc  (i,j) + wtxy_patch(np)*respc  (np)

               a_qcharge(i,j) = a_qcharge(i,j) + wtxy_patch(np)*qcharge(np)

! State and other variables
               a_t_grnd (i,j) = a_t_grnd (i,j) + wtxy_patch(np)*t_grnd (np)
               a_tlsun  (i,j) = a_tlsun  (i,j) + wtxy_patch(np)*tlsun  (np)
               a_tlsha  (i,j) = a_tlsha  (i,j) + wtxy_patch(np)*tlsha  (np)
               a_ldew   (i,j) = a_ldew   (i,j) + wtxy_patch(np)*ldew   (np)
               a_scv    (i,j) = a_scv    (i,j) + wtxy_patch(np)*scv    (np)
               a_snowdp (i,j) = a_snowdp (i,j) + wtxy_patch(np)*snowdp (np)
               a_fsno   (i,j) = a_fsno   (i,j) + wtxy_patch(np)*fsno   (np)
               a_sigf   (i,j) = a_sigf   (i,j) + wtxy_patch(np)*sigf   (np)
               a_green  (i,j) = a_green  (i,j) + wtxy_patch(np)*green  (np)
               a_lai    (i,j) = a_lai    (i,j) + wtxy_patch(np)*lai    (np)
               a_laisun (i,j) = a_laisun (i,j) + wtxy_patch(np)*laisun (np)
               a_laisha (i,j) = a_laisha (i,j) + wtxy_patch(np)*laisha (np)
               a_sai    (i,j) = a_sai    (i,j) + wtxy_patch(np)*sai    (np)
               a_alb(1,1,i,j) = a_alb(1,1,i,j) + wtxy_patch(np)*alb(1,1,np)
               a_alb(1,2,i,j) = a_alb(1,2,i,j) + wtxy_patch(np)*alb(1,2,np)
               a_alb(2,1,i,j) = a_alb(2,1,i,j) + wtxy_patch(np)*alb(2,1,np)
               a_alb(2,2,i,j) = a_alb(2,2,i,j) + wtxy_patch(np)*alb(2,2,np)
               a_emis   (i,j) = a_emis   (i,j) + wtxy_patch(np)*emis   (np)
               a_z0ma   (i,j) = a_z0ma   (i,j) + wtxy_patch(np)*z0ma   (np)
               a_tref   (i,j) = a_tref   (i,j) + wtxy_patch(np)*tref   (np)
               a_qref   (i,j) = a_qref   (i,j) + wtxy_patch(np)*qref   (np)
               a_xy_rain(i,j) = a_xy_rain(i,j) + wtxy_patch(np)*forc_rain(np)
               a_xy_snow(i,j) = a_xy_snow(i,j) + wtxy_patch(np)*forc_snow(np)

               ! local noon fluxes 
               call acc(sr     (np), wtxy_patch(np), a_sr     (i,j))
               call acc(solvd  (np), wtxy_patch(np), a_solvd  (i,j))
               call acc(solvi  (np), wtxy_patch(np), a_solvi  (i,j))
               call acc(solnd  (np), wtxy_patch(np), a_solnd  (i,j))
               call acc(solni  (np), wtxy_patch(np), a_solni  (i,j))
               call acc(srvd   (np), wtxy_patch(np), a_srvd   (i,j))
               call acc(srvi   (np), wtxy_patch(np), a_srvi   (i,j))
               call acc(srnd   (np), wtxy_patch(np), a_srnd   (i,j))
               call acc(srni   (np), wtxy_patch(np), a_srni   (i,j))
               call acc(solvdln(np), wtxy_patch(np), a_solvdln(i,j))
               call acc(solviln(np), wtxy_patch(np), a_solviln(i,j))
               call acc(solndln(np), wtxy_patch(np), a_solndln(i,j))
               call acc(solniln(np), wtxy_patch(np), a_solniln(i,j))
               call acc(srvdln (np), wtxy_patch(np), a_srvdln (i,j))
               call acc(srviln (np), wtxy_patch(np), a_srviln (i,j))
               call acc(srndln (np), wtxy_patch(np), a_srndln (i,j))
               call acc(srniln (np), wtxy_patch(np), a_srniln (i,j))

               area(i,j) = area_patch(np)
            ENDDO
         enddo
      ENDDO 
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j)
#endif
      DO j = 1, lat_points
         do i = 1, lon_points
            if(sumwt(i,j).gt.1.00001)then
               print*, 'summation of fraction patches (1) = ', sumwt(i,j), i,j
            !  call abort
            endif
            if(sumwt(i,j).gt.0.00001)then
               a_taux   (i,j) = a_taux   (i,j) / sumwt(i,j)
               a_tauy   (i,j) = a_tauy   (i,j) / sumwt(i,j)
               a_fsena  (i,j) = a_fsena  (i,j) / sumwt(i,j)
               a_lfevpa (i,j) = a_lfevpa (i,j) / sumwt(i,j)
               a_fevpa  (i,j) = a_fevpa  (i,j) / sumwt(i,j)
               a_fsenl  (i,j) = a_fsenl  (i,j) / sumwt(i,j)
               a_fevpl  (i,j) = a_fevpl  (i,j) / sumwt(i,j)
               a_etr    (i,j) = a_etr    (i,j) / sumwt(i,j)
               a_fseng  (i,j) = a_fseng  (i,j) / sumwt(i,j)
               a_fevpg  (i,j) = a_fevpg  (i,j) / sumwt(i,j)
               a_fgrnd  (i,j) = a_fgrnd  (i,j) / sumwt(i,j)
               a_sabvsun(i,j) = a_sabvsun(i,j) / sumwt(i,j)
               a_sabvsha(i,j) = a_sabvsha(i,j) / sumwt(i,j)
               a_sabg   (i,j) = a_sabg   (i,j) / sumwt(i,j)
               a_olrg   (i,j) = a_olrg   (i,j) / sumwt(i,j)
               a_rnet   (i,j) = a_rnet   (i,j) / sumwt(i,j) + a_xy_frl(i,j)
               a_xerr   (i,j) = a_xerr   (i,j) / sumwt(i,j)
               a_zerr   (i,j) = a_zerr   (i,j) / sumwt(i,j)
            
               a_rsur   (i,j) = a_rsur   (i,j) / sumwt(i,j)
               a_rnof   (i,j) = a_rnof   (i,j) / sumwt(i,j)
               a_qintr  (i,j) = a_qintr  (i,j) / sumwt(i,j)
               a_qinfl  (i,j) = a_qinfl  (i,j) / sumwt(i,j)
               a_qdrip  (i,j) = a_qdrip  (i,j) / sumwt(i,j)
               a_wat    (i,j) = a_wat    (i,j) / sumwt(i,j)
               a_assim  (i,j) = a_assim  (i,j) / sumwt(i,j)
               a_respc  (i,j) = a_respc  (i,j) / sumwt(i,j)

               a_qcharge(i,j) = a_qcharge(i,j) / sumwt(i,j)

               a_t_grnd (i,j) = a_t_grnd (i,j) / sumwt(i,j)
               a_tlsun  (i,j) = a_tlsun  (i,j) / sumwt(i,j)
               a_tlsha  (i,j) = a_tlsha  (i,j) / sumwt(i,j)
               a_ldew   (i,j) = a_ldew   (i,j) / sumwt(i,j)
               a_scv    (i,j) = a_scv    (i,j) / sumwt(i,j)
               a_snowdp (i,j) = a_snowdp (i,j) / sumwt(i,j)
               a_fsno   (i,j) = a_fsno   (i,j) / sumwt(i,j)
               a_sigf   (i,j) = a_sigf   (i,j) / sumwt(i,j)
               a_green  (i,j) = a_green  (i,j) / sumwt(i,j)
               a_lai    (i,j) = a_lai    (i,j) / sumwt(i,j)
               a_laisun (i,j) = a_laisun (i,j) / sumwt(i,j)
               a_laisha (i,j) = a_laisha (i,j) / sumwt(i,j)
               a_sai    (i,j) = a_sai    (i,j) / sumwt(i,j)
               a_alb(1,1,i,j) = a_alb(1,1,i,j) / sumwt(i,j)
               a_alb(1,2,i,j) = a_alb(1,2,i,j) / sumwt(i,j)
               a_alb(2,1,i,j) = a_alb(2,1,i,j) / sumwt(i,j)
               a_alb(2,2,i,j) = a_alb(2,2,i,j) / sumwt(i,j)
               a_emis   (i,j) = a_emis   (i,j) / sumwt(i,j)
               a_z0ma   (i,j) = a_z0ma   (i,j) / sumwt(i,j)
               a_trad   (i,j) = (a_olrg(i,j)/stefnc)**0.25
               a_tref   (i,j) = a_tref   (i,j) / sumwt(i,j)
               a_qref   (i,j) = a_qref   (i,j) / sumwt(i,j)
               a_xy_rain(i,j) = a_xy_rain(i,j) / sumwt(i,j)
               a_xy_snow(i,j) = a_xy_snow(i,j) / sumwt(i,j)
               
               if (a_sr     (i,j) /= spval) a_sr     (i,j) = a_sr     (i,j) / sumwt(i,j)
               if (a_solvd  (i,j) /= spval) a_solvd  (i,j) = a_solvd  (i,j) / sumwt(i,j)
               if (a_solvi  (i,j) /= spval) a_solvi  (i,j) = a_solvi  (i,j) / sumwt(i,j)
               if (a_solnd  (i,j) /= spval) a_solnd  (i,j) = a_solnd  (i,j) / sumwt(i,j)
               if (a_solni  (i,j) /= spval) a_solni  (i,j) = a_solni  (i,j) / sumwt(i,j)
               if (a_srvd   (i,j) /= spval) a_srvd   (i,j) = a_srvd   (i,j) / sumwt(i,j)
               if (a_srvi   (i,j) /= spval) a_srvi   (i,j) = a_srvi   (i,j) / sumwt(i,j)
               if (a_srnd   (i,j) /= spval) a_srnd   (i,j) = a_srnd   (i,j) / sumwt(i,j)
               if (a_srni   (i,j) /= spval) a_srni   (i,j) = a_srni   (i,j) / sumwt(i,j)
               if (a_solvdln(i,j) /= spval) a_solvdln(i,j) = a_solvdln(i,j) / sumwt(i,j)
               if (a_solviln(i,j) /= spval) a_solviln(i,j) = a_solviln(i,j) / sumwt(i,j)
               if (a_solndln(i,j) /= spval) a_solndln(i,j) = a_solndln(i,j) / sumwt(i,j)
               if (a_solniln(i,j) /= spval) a_solniln(i,j) = a_solniln(i,j) / sumwt(i,j)
               if (a_srvdln (i,j) /= spval) a_srvdln (i,j) = a_srvdln (i,j) / sumwt(i,j)
               if (a_srviln (i,j) /= spval) a_srviln (i,j) = a_srviln (i,j) / sumwt(i,j)
               if (a_srndln (i,j) /= spval) a_srndln (i,j) = a_srndln (i,j) / sumwt(i,j)
               if (a_srniln (i,j) /= spval) a_srniln (i,j) = a_srniln (i,j) / sumwt(i,j)
               
               mask(i,j) = 1
               frac(i,j) = sumwt(i,j)

            else
               a_taux   (i,j) = spval 
               a_tauy   (i,j) = spval 
               a_fsena  (i,j) = spval 
               a_lfevpa (i,j) = spval 
               a_fevpa  (i,j) = spval 
               a_fsenl  (i,j) = spval 
               a_fevpl  (i,j) = spval 
               a_etr    (i,j) = spval 
               a_fseng  (i,j) = spval 
               a_fevpg  (i,j) = spval 
               a_fgrnd  (i,j) = spval 
               a_sabvsun(i,j) = spval 
               a_sabvsha(i,j) = spval 
               a_sabg   (i,j) = spval 
               a_olrg   (i,j) = spval 
               a_rnet   (i,j) = spval 
               a_xerr   (i,j) = spval 
               a_zerr   (i,j) = spval 

               a_rsur   (i,j) = spval 
               a_rnof   (i,j) = spval 
               a_qintr  (i,j) = spval 
               a_qinfl  (i,j) = spval 
               a_qdrip  (i,j) = spval 
               a_wat    (i,j) = spval 
               a_assim  (i,j) = spval 
               a_respc  (i,j) = spval 

               a_qcharge(i,j) = spval 

               a_t_grnd (i,j) = spval 
               a_tlsun  (i,j) = spval 
               a_tlsha  (i,j) = spval 
               a_ldew   (i,j) = spval 
               a_scv    (i,j) = spval 
               a_snowdp (i,j) = spval 
               a_fsno   (i,j) = spval 
               a_sigf   (i,j) = spval 
               a_green  (i,j) = spval 
               a_lai    (i,j) = spval 
               a_laisun (i,j) = spval 
               a_laisha (i,j) = spval 
               a_sai    (i,j) = spval 
               a_alb(1,1,i,j) = spval 
               a_alb(1,2,i,j) = spval 
               a_alb(2,1,i,j) = spval 
               a_alb(2,2,i,j) = spval 
               a_emis   (i,j) = spval 
               a_z0ma   (i,j) = spval
               a_trad   (i,j) = spval
               a_tref   (i,j) = spval
               a_qref   (i,j) = spval
               a_xy_rain(i,j) = spval
               a_xy_snow(i,j) = spval

               a_sr     (i,j) = spval
               a_solvd  (i,j) = spval
               a_solvi  (i,j) = spval
               a_solnd  (i,j) = spval
               a_solni  (i,j) = spval
               a_srvd   (i,j) = spval
               a_srvi   (i,j) = spval
               a_srnd   (i,j) = spval
               a_srni   (i,j) = spval
               a_solvdln(i,j) = spval
               a_solviln(i,j) = spval
               a_solndln(i,j) = spval
               a_solniln(i,j) = spval
               a_srvdln (i,j) = spval
               a_srviln (i,j) = spval
               a_srndln (i,j) = spval
               a_srniln (i,j) = spval

               mask(i,j) = 0
               frac(i,j) = 0.


            endif

         enddo
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! --------------------------------------------------------------------
! Temperature and water (excluding land water bodies and ocean patches)
! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
! --------------------------------------------------------------------
      sumwt(:,:) = 0.
      a_t_soisno   (:,:,:) = 0.
      a_wliq_soisno(:,:,:) = 0.
      a_wice_soisno(:,:,:) = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         do i = 1, lon_points

            if (spatch_xy(i,j) .le. 0) cycle
            DO np = spatch_xy(i,j), epatch_xy(i,j)
               if(itypwat(np) <= 3)then  ! excluded the land water bodies and ocean patches
                  sumwt(i,j) = sumwt(i,j) + wtxy_patch(np)
                  a_t_soisno   (maxsnl+1:nl_soil,i,j) = a_t_soisno   (maxsnl+1:nl_soil,i,j) &
                     + wtxy_patch(np)*t_soisno   (maxsnl+1:nl_soil,np)
                  a_wliq_soisno(maxsnl+1:nl_soil,i,j) = a_wliq_soisno(maxsnl+1:nl_soil,i,j) &
                     + wtxy_patch(np)*wliq_soisno(maxsnl+1:nl_soil,np)
                  a_wice_soisno(maxsnl+1:nl_soil,i,j) = a_wice_soisno(maxsnl+1:nl_soil,i,j) &
                     + wtxy_patch(np)*wice_soisno(maxsnl+1:nl_soil,np)
               endif
            ENDDO

         enddo
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            if(sumwt(i,j).gt.1.00001)then
               print*, 'summation of fraction patches (2) = ', sumwt(i,j), i,j
            !  call abort
            endif
            if(sumwt(i,j).gt.0.00001)then
               do l = maxsnl+1, nl_soil
                  a_t_soisno   (l,i,j) = a_t_soisno   (l,i,j) / sumwt(i,j)
                  a_wliq_soisno(l,i,j) = a_wliq_soisno(l,i,j) / sumwt(i,j)
                  a_wice_soisno(l,i,j) = a_wice_soisno(l,i,j) / sumwt(i,j)
               enddo
            else
               a_t_soisno   (maxsnl+1:nl_soil,i,j) = spval 
               a_wliq_soisno(maxsnl+1:nl_soil,i,j) = spval 
               a_wice_soisno(maxsnl+1:nl_soil,i,j) = spval 
            endif

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! --------------------------------------------------------------------
! additial diagnostic variables for output (vegetated land only <=2)
! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
! --------------------------------------------------------------------
      sumwt(:,:) = 0.
      a_h2osoi (:,:,:) = 0.
      a_rstfac (:,:)   = 0.
      a_zwt    (:,:)   = 0.
      a_wa     (:,:)   = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         do i = 1, lon_points

            if (spatch_xy(i,j) .le. 0) cycle
            DO np = spatch_xy(i,j), epatch_xy(i,j)
               if(itypwat(np) <= 2)then  ! excluded the land water bodies and ocean patches
                  sumwt(i,j) = sumwt(i,j) + wtxy_patch(np)
                  a_h2osoi(1:nl_soil,i,j) = a_h2osoi(1:nl_soil,i,j) &
                     + wtxy_patch(np)*h2osoi(1:nl_soil,np)
                  a_rstfac(i,j) = a_rstfac(i,j) + wtxy_patch(np)*rstfac(np)
                  a_zwt   (i,j) = a_zwt   (i,j) + wtxy_patch(np)*zwt   (np)
                  a_wa    (i,j) = a_wa    (i,j) + wtxy_patch(np)*wa    (np)
               endif
            ENDDO

         enddo
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            if(sumwt(i,j).gt.1.00001)then
               print*, 'summation of fraction patches (2) = ', sumwt(i,j), i,j
            !  call abort
            endif
            if(sumwt(i,j).gt.0.00001)then
               do l = 1, nl_soil
                  a_h2osoi (l,i,j) = a_h2osoi (l,i,j) / sumwt(i,j)
               enddo
               a_rstfac (i,j) = a_rstfac (i,j) / sumwt(i,j)
               a_zwt    (i,j) = a_zwt    (i,j) / sumwt(i,j)
               a_wa     (i,j) = a_wa     (i,j) / sumwt(i,j)
            else
               a_h2osoi (:,i,j) = spval
               a_rstfac (i,j)   = spval
               a_zwt    (i,j)   = spval
               a_wa     (i,j)   = spval
            endif

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! -----------------------------------------------
! Land water bodies' ice fraction and temperature
! -----------------------------------------------
      sumwt(:,:) = 0.
      a_t_lake(:,:,:) = 0.
      a_lake_icefrac(:,:,:) = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         do i = 1, lon_points

            if (spatch_xy(i,j) .le. 0) cycle
            DO np = spatch_xy(i,j), epatch_xy(i,j)
               if(itypwat(np) == 4)then  ! land water bodies only
                  sumwt(i,j) = sumwt(i,j) + wtxy_patch(np)
                  a_t_lake(1:nl_lake,i,j) = a_t_lake(1:nl_lake,i,j) + wtxy_patch(np)*t_lake(1:nl_lake,np)
                  a_lake_icefrac(1:nl_lake,i,j) = a_lake_icefrac(1:nl_lake,i,j) + wtxy_patch(np)*lake_icefrac(1:nl_lake,np)
               endif
            ENDDO

         enddo
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            if(sumwt(i,j).gt.1.00001)then
               print*, 'summation of fraction patches (3) = ', sumwt(i,j), i,j
            !  call abort
            endif
            if(sumwt(i,j).gt.0.00001)then
               do l = 1, nl_lake
                  a_t_lake(l,i,j) = a_t_lake(l,i,j) / sumwt(i,j)
                  a_lake_icefrac(l,i,j) = a_lake_icefrac(l,i,j) / sumwt(i,j)
               enddo
            else
               a_t_lake(1:nl_lake,i,j) = spval
               a_lake_icefrac(1:nl_lake,i,j) = spval
            endif
         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


! --------------------------------
! Retrieve through averaged fluxes
! --------------------------------
      sumwt(:,:) = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         do i = 1, lon_points

            if (spatch_xy(i,j) .le. 0) cycle
            DO np = spatch_xy(i,j), epatch_xy(i,j)
               sumwt(i,j) = sumwt(i,j) + wtxy_patch(np)
            ENDDO

         enddo
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j) &
!$OMP PRIVATE(rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q) &
!$OMP PRIVATE(z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf) &
!$OMP PRIVATE(obu,fh2m,fq2m) &
!$OMP PRIVATE(um,thvstar,beta,zii,wc,wc2)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            if(sumwt(i,j) > 0.00001)then
               z0m_av = a_z0ma(i,j) 
               z0h_av = a_z0ma(i,j)
               z0q_av = a_z0ma(i,j)
               displa_av = 2./3.*z0m_av/0.07

               hgt_u = max(forc_xy_hgt_u(i,j),5.+displa_av)
               hgt_t = max(forc_xy_hgt_t(i,j),5.+displa_av)
               hgt_q = max(forc_xy_hgt_q(i,j),5.+displa_av)
               zldis = hgt_u-displa_av

               us = forc_xy_us(i,j)
               vs = forc_xy_vs(i,j)
               tm = forc_xy_t(i,j)
               qm = forc_xy_q(i,j)
               psrf = forc_xy_psrf(i,j)
               rhoair = (psrf - 0.378*qm*psrf/(0.622+0.378*qm)) / (rgas*tm)


               a_ustar(i,j) = sqrt(max(1.e-6,sqrt(a_taux(i,j)**2+a_tauy(i,j)**2))/rhoair)

               a_tstar(i,j) = -a_fsena(i,j)/(rhoair*a_ustar(i,j))/cpair
               a_qstar(i,j) = -a_fevpa(i,j)/(rhoair*a_ustar(i,j))
   
               thm = tm + 0.0098*hgt_t
               th = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)

               a_zol(i,j) = zldis*vonkar*grav &
                          * (a_tstar(i,j)+0.61*th*a_qstar(i,j))/(a_ustar(i,j)**2*thv)

               if(a_zol(i,j) >= 0.)then   !stable
                  a_zol(i,j) = min(2.,max(a_zol(i,j),1.e-6))
               else                       !unstable
                  a_zol(i,j) = max(-100.,min(a_zol(i,j),-1.e-6))
               endif

               beta = 1.
               zii = 1000.
               thvstar=a_tstar(i,j)+0.61*th*a_qstar(i,j)
               ur = sqrt(us*us+vs*vs)
               if(a_zol(i,j) >= 0.)then
                  um = max(ur,0.1)
               else
                  wc = (-grav*a_ustar(i,j)*thvstar*zii/thv)**(1./3.)
                 wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               endif

               obu = zldis/a_zol(i,j)
               call moninobuk(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                    obu,um,a_ustar(i,j),fh2m,fq2m,&
                    a_fm10m(i,j),a_fm(i,j),a_fh(i,j),a_fq(i,j))

! bug found by chen qiying 2013/07/01 
               a_rib(i,j) = a_zol(i,j)/vonkar*a_ustar(i,j)**2/(vonkar/a_fh(i,j)*um**2)
               a_rib(i,j) = min(5.,a_rib(i,j))

               a_us10m(i,j) = us/um * a_ustar(i,j)/vonkar * a_fm10m(i,j)
               a_vs10m(i,j) = vs/um * a_ustar(i,j)/vonkar * a_fm10m(i,j)

            else

               a_ustar(i,j) = spval
               a_tstar(i,j) = spval
               a_qstar(i,j) = spval
               a_zol(i,j)   = spval
               a_rib(i,j)   = spval
               a_fm(i,j)    = spval
               a_fh(i,j)    = spval
               a_fq(i,j)    = spval
               a_fm10m(i,j) = spval
               a_us10m(i,j) = spval
               a_vs10m(i,j) = spval

            endif

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


! ---------------------------------------------------
! ACCUMULATION in each time step
! ---------------------------------------------------
      nac = nac + 1
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      do j = 1, lat_points
         do i = 1, lon_points

            call acc(a_xy_us     (i,j), 1., f_xy_us     (i,j))
            call acc(a_xy_vs     (i,j), 1., f_xy_vs     (i,j))
            call acc(a_xy_t      (i,j), 1., f_xy_t      (i,j))
            call acc(a_xy_q      (i,j), 1., f_xy_q      (i,j))
            call acc(a_xy_prc    (i,j), 1., f_xy_prc    (i,j))
            call acc(a_xy_prl    (i,j), 1., f_xy_prl    (i,j))
            call acc(a_xy_pbot   (i,j), 1., f_xy_pbot   (i,j))
            call acc(a_xy_frl    (i,j), 1., f_xy_frl    (i,j))
            call acc(a_xy_solarin(i,j), 1., f_xy_solarin(i,j))

            call acc(a_taux   (i,j), 1., f_taux   (i,j))
            call acc(a_tauy   (i,j), 1., f_tauy   (i,j))
            call acc(a_fsena  (i,j), 1., f_fsena  (i,j))
            call acc(a_lfevpa (i,j), 1., f_lfevpa (i,j))
            call acc(a_fevpa  (i,j), 1., f_fevpa  (i,j))
            call acc(a_fsenl  (i,j), 1., f_fsenl  (i,j))
            call acc(a_fevpl  (i,j), 1., f_fevpl  (i,j))
            call acc(a_etr    (i,j), 1., f_etr    (i,j))
            call acc(a_fseng  (i,j), 1., f_fseng  (i,j))
            call acc(a_fevpg  (i,j), 1., f_fevpg  (i,j))
            call acc(a_fgrnd  (i,j), 1., f_fgrnd  (i,j))
            call acc(a_sabvsun(i,j), 1., f_sabvsun(i,j))
            call acc(a_sabvsha(i,j), 1., f_sabvsha(i,j))
            call acc(a_sabg   (i,j), 1., f_sabg   (i,j))
            call acc(a_olrg   (i,j), 1., f_olrg   (i,j))
            call acc(a_rnet   (i,j), 1., f_rnet   (i,j))
            call acc(a_xerr   (i,j), 1., f_xerr   (i,j))
            call acc(a_zerr   (i,j), 1., f_zerr   (i,j))
            call acc(a_rsur   (i,j), 1., f_rsur   (i,j))
            call acc(a_rnof   (i,j), 1., f_rnof   (i,j))
            call acc(a_qintr  (i,j), 1., f_qintr  (i,j))
            call acc(a_qinfl  (i,j), 1., f_qinfl  (i,j))
            call acc(a_qdrip  (i,j), 1., f_qdrip  (i,j))
            call acc(a_rstfac (i,j), 1., f_rstfac (i,j))
            call acc(a_zwt    (i,j), 1., f_zwt    (i,j))
            call acc(a_wa     (i,j), 1., f_wa     (i,j))
            call acc(a_wat    (i,j), 1., f_wat    (i,j))
            call acc(a_assim  (i,j), 1., f_assim  (i,j))
            call acc(a_respc  (i,j), 1., f_respc  (i,j))

            call acc(a_qcharge(i,j), 1., f_qcharge(i,j))

            call acc(a_t_grnd (i,j), 1., f_t_grnd (i,j))
            call acc(a_tlsun  (i,j), 1., f_tlsun  (i,j))
            call acc(a_tlsha  (i,j), 1., f_tlsha  (i,j))
            call acc(a_ldew   (i,j), 1., f_ldew   (i,j))
            call acc(a_scv    (i,j), 1., f_scv    (i,j))
            call acc(a_snowdp (i,j), 1., f_snowdp (i,j))
            call acc(a_fsno   (i,j), 1., f_fsno   (i,j))
            call acc(a_sigf   (i,j), 1., f_sigf   (i,j))
            call acc(a_green  (i,j), 1., f_green  (i,j))
            call acc(a_lai    (i,j), 1., f_lai    (i,j))
            call acc(a_laisun (i,j), 1., f_laisun (i,j))
            call acc(a_laisha (i,j), 1., f_laisha (i,j))
            call acc(a_sai    (i,j), 1., f_sai    (i,j))
            call acc(a_alb(1,1,i,j), 1., f_alb(1,1,i,j))
            call acc(a_alb(2,1,i,j), 1., f_alb(2,1,i,j))
            call acc(a_alb(1,2,i,j), 1., f_alb(1,2,i,j))
            call acc(a_alb(2,2,i,j), 1., f_alb(2,2,i,j))
            call acc(a_emis   (i,j), 1., f_emis   (i,j))
            call acc(a_z0ma   (i,j), 1., f_z0ma   (i,j))
            call acc(a_trad   (i,j), 1., f_trad   (i,j))
            call acc(a_tref   (i,j), 1., f_tref   (i,j))
            call acc(a_qref   (i,j), 1., f_qref   (i,j))
            call acc(a_xy_rain(i,j), 1., f_xy_rain(i,j))
            call acc(a_xy_snow(i,j), 1., f_xy_snow(i,j))

            do l = maxsnl+1, nl_soil
               call acc(a_t_soisno   (l,i,j), 1., f_t_soisno   (l,i,j))
               call acc(a_wliq_soisno(l,i,j), 1., f_wliq_soisno(l,i,j))
               call acc(a_wice_soisno(l,i,j), 1., f_wice_soisno(l,i,j))
            end do
            
            do l = 1, nl_soil
               call acc(a_h2osoi     (l,i,j), 1., f_h2osoi     (l,i,j))
            end do

            do l = 1, nl_lake
               call acc(a_t_lake(l,i,j), 1., f_t_lake(l,i,j))
               call acc(a_lake_icefrac(l,i,j), 1., f_lake_icefrac(l,i,j))
            end do

            call acc(a_ustar(i,j), 1., f_ustar(i,j))
            call acc(a_tstar(i,j), 1., f_tstar(i,j))
            call acc(a_qstar(i,j), 1., f_qstar(i,j))
            call acc(a_zol  (i,j), 1., f_zol  (i,j))
            call acc(a_rib  (i,j), 1., f_rib  (i,j))
            call acc(a_fm   (i,j), 1., f_fm   (i,j))
            call acc(a_fh   (i,j), 1., f_fh   (i,j))
            call acc(a_fq   (i,j), 1., f_fq   (i,j))

            call acc(a_us10m(i,j), 1., f_us10m(i,j))
            call acc(a_vs10m(i,j), 1., f_vs10m(i,j))
            call acc(a_fm10m(i,j), 1., f_fm10m(i,j))

            call acc(a_sr     (i,j), 1., f_sr     (i,j))
            call acc(a_solvd  (i,j), 1., f_solvd  (i,j))
            call acc(a_solvi  (i,j), 1., f_solvi  (i,j))
            call acc(a_solnd  (i,j), 1., f_solnd  (i,j))
            call acc(a_solni  (i,j), 1., f_solni  (i,j))
            call acc(a_srvd   (i,j), 1., f_srvd   (i,j))
            call acc(a_srvi   (i,j), 1., f_srvi   (i,j))
            call acc(a_srnd   (i,j), 1., f_srnd   (i,j))
            call acc(a_srni   (i,j), 1., f_srni   (i,j))
            call acc(a_solvdln(i,j), 1., f_solvdln(i,j))
            call acc(a_solviln(i,j), 1., f_solviln(i,j))
            call acc(a_solndln(i,j), 1., f_solndln(i,j))
            call acc(a_solniln(i,j), 1., f_solniln(i,j))
            call acc(a_srvdln (i,j), 1., f_srvdln (i,j))
            call acc(a_srviln (i,j), 1., f_srviln (i,j))
            call acc(a_srndln (i,j), 1., f_srndln (i,j))
            call acc(a_srniln (i,j), 1., f_srniln (i,j))

            if (a_solvdln(i,j) /= spval) nac_ln(i,j) = nac_ln(i,j) + 1

         end do
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

END SUBROUTINE vec2xy

SUBROUTINE acc(var, wgt, s)
   
   use precision
   use MOD_TimeInvariants, only: spval

   IMPLICIT NONE

   real(r8), intent(in)  :: var
   real(r8), intent(in)  :: wgt
   real(r8), intent(out) :: s

   if (var /= spval) then
      if (s /= spval) then
         s = s + wgt*var
      else
         s = wgt*var
      end if
   end if

END SUBROUTINE acc

! ----- EOP ---------
