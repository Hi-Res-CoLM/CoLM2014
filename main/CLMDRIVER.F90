#include <define.h>

SUBROUTINE CLMDRIVER (nl_soil,maxsnl,nl_lake,numpatch,idate,deltim,&
                      dolai,doalb,dosst,oro)


!=======================================================================
!
! CLM MODEL DRIVER
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002, 03/2014
!
!=======================================================================

 use precision
 use PhysicalConstants, only : tfrz, rgas, vonkar
 use MOD_TimeInvariants
 use MOD_TimeVariables
 use MOD_1D_Forcing
 use MOD_1D_Fluxes
 use omp_lib

 IMPLICIT NONE

  integer,  INTENT(in) :: nl_soil  ! number of soil layers
  integer,  INTENT(in) :: nl_lake  ! number of lake layers
  integer,  INTENT(in) :: maxsnl   ! max number of snow layers
  integer,  INTENT(in) :: numpatch ! number of clm grid points
  integer,  INTENT(in) :: idate(3) ! model calendar for next time step (year, julian day, seconds)
  real(r8), INTENT(in) :: deltim   ! seconds in a time-step

  logical,  INTENT(in) :: dolai    ! true if time for time-varying vegetation paramter
  logical,  INTENT(in) :: doalb    ! true if time for surface albedo calculation
  logical,  INTENT(in) :: dosst    ! true if time for update sst/ice/snow

  real(r8), INTENT(inout) :: oro(numpatch)  ! ocean(0)/seaice(2)/ flag

! -------------- Local varaibles -------------------
  real(r8) :: z_soisno (maxsnl+1:nl_soil,numpatch)
  real(r8) :: dz_soisno(maxsnl+1:nl_soil,numpatch)

  integer :: i

! ======================================================================

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP)
#endif
      DO i = 1, numpatch

         z_soisno (maxsnl+1:0,i) = z_sno (maxsnl+1:0,i)
         z_soisno (1:nl_soil ,i) = z_soi (1:nl_soil ,i)
         dz_soisno(maxsnl+1:0,i) = dz_sno(maxsnl+1:0,i)
         dz_soisno(1:nl_soil ,i) = dz_soi(1:nl_soil ,i)

         CALL CLMMAIN (deltim,doalb,dolai,dosst,oro(i), &
         nl_soil, nl_lake, maxsnl, i, &

         dlon(i), dlat(i), mxy_patch(i), itypwat(i), &
         lakedepth(i), dz_lake(1:,i), &                     ! new lake scheme

       ! SOIL INFORMATION AND LAKE DEPTH
         soil_s_v_alb(i), soil_d_v_alb(i), soil_s_n_alb(i), soil_d_n_alb(i), &
         porsl(1:,i),     psi0(1:,i),      bsw(1:,i),       hksati(1:,i),    &
         csol(1:,i),      dksatu(1:,i),    dkdry(1:,i),     rootfr(1:,i),    &

       ! VEGETATION INFORMATION
         z0m(i),          displa(i),       sqrtdi(i),                        &
         effcon(i),       vmax25(i),       slti(i),         hlti(i),         &
         shti(i),         hhti(i),         trda(i),         trdm(i),         &
         trop(i),         gradm(i),        binter(i),       extkn(i),        &
         chil(i),         ref(1:,1:,i),    tran(1:,1:,i),                    &

       ! ATMOSPHERIC FORCING
         forc_pco2m(i),   forc_po2m(i),    forc_us(i),      forc_vs(i),      &
         forc_t(i),       forc_q(i),       forc_prc(i),     forc_prl(i),     &
         forc_rain(i),    forc_snow(i),    forc_psrf(i),    forc_pbot(i),    &
         forc_sols(i),    forc_soll(i),    forc_solsd(i),   forc_solld(i),   &
         forc_frl(i),     forc_hgt_u(i),   forc_hgt_t(i),   forc_hgt_q(i),   &
         forc_rhoair(i),                                                     &

       ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
         idate,                                                              &
         z_soisno(maxsnl+1:,i),            dz_soisno(maxsnl+1:,i),           &
         t_soisno(maxsnl+1:,i),            wliq_soisno(maxsnl+1:,i),         &
         wice_soisno(maxsnl+1:,i),                                           &

         t_grnd(i),       tlsun(i),        tlsha(i),        ldew(i),         &
         sag(i),          scv(i),          snowdp(i),       fveg(i),         &
         fsno(i),         sigf(i),         green(i),        lai(i),          &
         sai(i),          coszen(i),       albg(1:,1:,i),   albv(1:,1:,i),   &
         alb(1:,1:,i),    ssun(1:,1:,i),   ssha(1:,1:,i),   thermk(i),       &
         extkb(i),        extkd(i),        &
        
         zwt(i),          wa(i),                                             &
         t_lake(1:,i),    lake_icefrac(1:,i),                                & ! new lake scheme

       ! additional diagnostic variables for output
         laisun(i),       laisha(i),                                         &
         rstfac(i),       h2osoi(1:,i),    wat(i),                           &

       ! FLUXES
         taux(i),         tauy(i),         fsena(i),        fevpa(i),        &
         lfevpa(i),       fsenl(i),        fevpl(i),        etr(i),          &
         fseng(i),        fevpg(i),        olrg(i),         fgrnd(i),        &
         trad(i),         tref(i),         qref(i),         rsur(i),         &
         rnof(i),         qintr(i),        qinfl(i),        qdrip(i),        &
         rst(i),          assim(i),        respc(i),        sabvsun(i),      &
         sabvsha(i),      sabg(i),         sr(i),           solvd(i),        &
         solvi(i),        solnd(i),        solni(i),        srvd(i),         &
         srvi(i),         srnd(i),         srni(i),         solvdln(i),      &
         solviln(i),      solndln(i),      solniln(i),      srvdln(i),       &
         srviln(i),       srndln(i),       srniln(i),       qcharge(i),      &
         xerr(i),         zerr(i),                                           &

       ! TUNABLE modle constants
         zlnd,            zsno,            csoilc,          dewmx,           &
         wtfact,          capr,            cnfac,           ssi,             &
         wimp,            pondmx,          smpmax,          smpmin,          &
         trsmx0,          tcrit,                                             &

       ! additional variables required by coupling with WRF model
         emis(i),         z0ma(i),         zol(i),          rib(i),          &
         ustar(i),        qstar(i),        tstar(i),                         &
         fm(i),           fh(i),           fq(i) )


         z_sno (maxsnl+1:0,i) = z_soisno (maxsnl+1:0,i)
         dz_sno(maxsnl+1:0,i) = dz_soisno(maxsnl+1:0,i)


      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


END SUBROUTINE CLMDRIVER
! --------- EOP --------
