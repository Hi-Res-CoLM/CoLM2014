
#include <define.h>

SUBROUTINE IniTimeVar(nl_soil,maxsnl,itypwat&
	             ,porsl,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb&
                     ,z0m,zlnd,chil,ref,tran,z_soisno,dz_soisno&
                     ,t_soisno,wliq_soisno,wice_soisno,zwt,wa&
                     ,t_grnd,tlsun,tlsha,ldew,sag,scv&   
                     ,snowdp,fveg,fsno,sigf,green,lai,sai,coszen&
                     ,albg,albv,alb,ssun,ssha,thermk,extkb,extkd&
                     ,trad,tref,qref,rst,emis,z0ma,zol,rib&
                     ,ustar,qstar,tstar,fm,fh,fq&
#if(defined SOILINI)
                     ,nl_soil_ini,soil_z,soil_t,soil_w,snow_d)
#else
                     )
#endif

!=======================================================================
! Created by Yongjiu Dai, 09/15/1999
! Revised by Yongjiu Dai, 08/30/2002
!                         03/2014
!=======================================================================

  use precision
  use PhysicalConstants, only : tfrz
  use ALBEDO

  IMPLICIT NONE 

  integer, INTENT(in) ::        &! 
        nl_soil,                &! soil layer number
        maxsnl,                 &! maximum snow layer number
        itypwat                  ! index for land cover type [-]

  real(r8), INTENT(in) ::       &!
        fveg,                   &! fraction of vegetation cover
        green,                  &! leaf greenness
        lai,                    &! leaf area index
        sai,                    &! stem area index
        coszen,                 &! cosine of solar zenith angle
        soil_s_v_alb,           &! albedo of visible of the saturated soil
        soil_d_v_alb,           &! albedo of visible of the dry soil
        soil_s_n_alb,           &! albedo of near infrared of the saturated soil
        soil_d_n_alb,           &! albedo of near infrared of the dry soil
	z0m,                    &! aerodynamic roughness length [m]
        zlnd,                   &! aerodynamic roughness length over soil surface [m] 
        chil,                   &! leaf angle distribution factor
        ref (2,2),              &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2),              &! leaf transmittance (iw=iband, il=life and dead)
        porsl(1:nl_soil)         ! porosity of soil

#if(defined SOILINI)
  integer, INTENT(in) :: nl_soil_ini
  real(r8), INTENT(in) ::       &!
        soil_z(nl_soil_ini),    &! soil layer depth for initial (m)
        soil_t(nl_soil_ini),    &! soil temperature from initial file (K)
        soil_w(nl_soil_ini),    &! soil wetness from initial file (-)
        snow_d                   ! snow depth (m)
#endif

  real(r8), INTENT(inout) ::    &!
        z_soisno (maxsnl+1:nl_soil),   &! node depth [m]
        dz_soisno(maxsnl+1:nl_soil)     ! layer thickness [m]

  real(r8), INTENT(out) ::      &!
        t_soisno (maxsnl+1:nl_soil), &! soil temperature [K]
        wliq_soisno(maxsnl+1:nl_soil), &! liquid water in layers [kg/m2]
        wice_soisno(maxsnl+1:nl_soil), &! ice lens in layers [kg/m2]
        t_grnd,                 &! ground surface temperature [K]
        tlsun,                  &! sunlit leaf temperature [K]
        tlsha,                  &! shaded leaf temperature [K]
        ldew,                   &! depth of water on foliage [mm]
        sag,                    &! non dimensional snow age [-]
        scv,                    &! snow cover, water equivalent [mm]
        snowdp,                 &! snow depth [meter]
        fsno,                   &! fraction of snow cover on ground
        sigf,                   &! fraction of veg cover, excluding snow-covered veg [-]

        albg(2,2),              &! albedo, ground [-]
        albv(2,2),              &! albedo, vegetation [-]
        alb (2,2),              &! averaged albedo [-]
        ssun(2,2),              &! sunlit canopy absorption for solar radiation
        ssha(2,2),              &! shaded canopy absorption for solar radiation
        thermk,                 &! canopy gap fraction for tir radiation
        extkb,                  &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,                  &! diffuse and scattered diffuse PAR extinction coefficient
        wa,                     &! water storage in aquifer [mm]
        zwt,                    &! the depth to water table [m]

                    ! Additional variables required by reginal model (WRF & RSM) 
        trad,                   &! radiative temperature of surface [K]
        tref,                   &! 2 m height air temperature [kelvin]
        qref,                   &! 2 m height air specific humidity
        rst,                    &! canopy stomatal resistance (s/m)
        emis,                   &! averaged bulk surface emissivity
        z0ma,                   &! effective roughness [m]
        zol,                    &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,                    &! bulk Richardson number in surface layer
        ustar,                  &! u* in similarity theory [m/s]
        qstar,                  &! q* in similarity theory [kg/kg]
        tstar,                  &! t* in similarity theory [K]
        fm,                     &! integral of profile function for momentum
        fh,                     &! integral of profile function for heat
        fq                       ! integral of profile function for moisture

        integer j, snl                      
        real(r8) wet(nl_soil), wt, ssw, oro, rhosno_ini, a
!-----------------------------------------------------------------------

  if(itypwat <= 5)then ! land grid
     rhosno_ini = 250.
#if(defined SOILINI)
     do j = 1, nl_soil
        call polint(soil_z,soil_t,nl_soil_ini,z_soisno(j),t_soisno(j))
        call polint(soil_z,soil_w,nl_soil_ini,z_soisno(j),wet(j))
        a = min(soil_t(1),soil_t(2),soil_t(3))-5.
        t_soisno(j) = max(t_soisno(j), a)
        a = max(soil_t(1),soil_t(2),soil_t(3))+5.
        t_soisno(j) = min(t_soisno(j), a)

        a = min(soil_w(1),soil_w(2),soil_w(3))
	wet(j) = max(wet(j), a, 0.1)
        a = max(soil_w(1),soil_w(2),soil_w(3))
	wet(j) = min(wet(j), a, 0.5)

        if(t_soisno(j).ge.tfrz)then
           wliq_soisno(j) = wet(j)*dz_soisno(j)*1000.
!          wliq_soisno(j) = porsl(j)*wet(j)*dz_soisno(j)*1000.
           wice_soisno(j) = 0.
        else
           wliq_soisno(j) = 0.
           wice_soisno(j) = wet(j)*dz_soisno(j)*1000.
!          wliq_soisno(j) = porsl(j)*wet(j)*dz_soisno(j)*1000.
        endif
     enddo

     snowdp = snow_d
     sag    = 0.
     scv    = snowdp*rhosno_ini

     call snowfraction (fveg,z0m,zlnd,scv,snowdp,wt,sigf,fsno)
     call snow_ini (itypwat,maxsnl,snowdp,snl,z_soisno,dz_soisno)

     if(snl.lt.0)then
        do j = snl+1, 0
           t_soisno(j) = min(tfrz-1., t_soisno(1))
           wliq_soisno(j) = 0.
           wice_soisno(j) = dz_soisno(j)*rhosno_ini         ! m * kg m-3 = kg m-2
        enddo
     endif

     if(snl>maxsnl)then
        t_soisno (maxsnl+1:snl) = -999.
        wice_soisno(maxsnl+1:snl) = 0.
        wliq_soisno(maxsnl+1:snl) = 0.
        z_soisno   (maxsnl+1:snl) = 0.
        dz_soisno  (maxsnl+1:snl) = 0.
     endif

     ldew  = 0.
     tlsun = t_soisno(1)
     tlsha = t_soisno(1)
     t_grnd = t_soisno(1)
#else
! soil temperature and water content
     do j = 1, nl_soil
        if(itypwat==3)then ! land ice 
           t_soisno(j) = 253.
           wliq_soisno(j) = 0.
           wice_soisno(j) = dz_soisno(j)*1000.
        else
           t_soisno(j) = 283.
           wliq_soisno(j) = dz_soisno(j)*porsl(j)*1000.
           wice_soisno(j) = 0.
        endif
     enddo

! water table depth (initially at 1.0 m below the model bottom; wa when zwt
!                    is below the model bottom zi(nl_soil)

     wa  = 4800.                             ! assuming aquifer capacity is 5000 mm
     zwt = (25. + z_soisno(nl_soil))+dz_soisno(nl_soil)/2. - wa/1000./0.2 ! to result in zwt = zi(nl_soil) + 1.0 m

! snow temperature and water content
     t_soisno(maxsnl+1:0) = -999.
     wice_soisno(maxsnl+1:0) = 0.
     wliq_soisno(maxsnl+1:0) = 0.
     z_soisno (maxsnl+1:0) = 0.
     dz_soisno(maxsnl+1:0) = 0.

     sigf   = fveg
     fsno   = 0.
     ldew   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.
     tlsun  = t_soisno(1)
     tlsha  = t_soisno(1)
     t_grnd = t_soisno(1)
#endif

! surface albedo
     ssw = min(1.,1.e-3*wliq_soisno(1)/dz_soisno(1))
     call albland (itypwat,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb,&
                   chil,ref,tran,fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,t_grnd,&
                   alb,albg,albv,ssun,ssha,thermk,extkb,extkd)
  else                 ! ocean grid
     t_soisno(:) = 300.
     wice_soisno(:) = 0.
     wliq_soisno(:) = 1000.
     z_soisno (maxsnl+1:0) = 0.
     dz_soisno(maxsnl+1:0) = 0.
     sigf   = 0.
     fsno   = 0.
     ldew   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.
     tlsun  = 300.
     tlsha  = 300.
     t_grnd = 300.

     oro = 0
     call albocean (oro,scv,coszen,alb)
     albg(:,:) = alb(:,:)
     albv(:,:) = 0.0
     ssun(:,:) = 0.0
     ssha(:,:) = 0.0
     thermk = 0.0
     extkb = 0.0
     extkd = 0.0
  endif

! Additional variables required by reginal model (WRF & RSM)
! totally arbitrarily assigned here
  trad  = t_grnd      
  tref  = t_grnd      
  qref  = 0.3     
  rst   = 1.e36   
  emis  = 1.0     
  z0ma  = 0.01    
  zol   = -1.0    
  rib   = -0.1    
  ustar = 0.25    
  qstar = 0.001   
  tstar = -1.5    
  fm    = alog(30.)  
  fh    = alog(30.)  
  fq    = alog(30.)  

END SUBROUTINE IniTimeVar
!-----------------------------------------------------------------------
! EOP


SUBROUTINE snow_ini(itypwat,maxsnl,snowdp,snl,z_soisno,dz_soisno)

! Snow spatial discretization initially

  use precision
  implicit none

  integer, intent(in) :: maxsnl  ! maximum of snow layers
  integer, intent(in) :: itypwat ! index for land cover type [-]
  real(r8), intent(in) :: snowdp ! snow depth [m]
  real(r8), intent(out) :: z_soisno (maxsnl+1:0) ! node depth [m]
  real(r8), intent(out) :: dz_soisno(maxsnl+1:0) ! layer thickness [m]
  integer, intent(out) :: snl ! number of snow layer
  real(r8) zi
  integer i
!-----------------------------------------------------------------------

  dz_soisno(:0) = 0.
  z_soisno(:0) = 0.
  snl = 0
  if(itypwat.le.3)then ! non water bodies

     if(snowdp.lt.0.01)then
        snl = 0
     else
        if(snowdp>=0.01 .and. snowdp<=0.03)then
           snl = -1
           dz_soisno(0)  = snowdp
        else if(snowdp>0.03 .and. snowdp<=0.04)then
           snl = -2
           dz_soisno(-1) = snowdp/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.04 .and. snowdp<=0.07)then
           snl = -2
           dz_soisno(-1) = 0.02
           dz_soisno( 0) = snowdp - dz_soisno(-1)
        else if(snowdp>0.07 .and. snowdp<=0.12)then
           snl = -3
           dz_soisno(-2) = 0.02
           dz_soisno(-1) = (snowdp - 0.02)/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.12 .and. snowdp<=0.18)then
           snl = -3
           dz_soisno(-2) = 0.02
           dz_soisno(-1) = 0.05
           dz_soisno( 0) = snowdp - dz_soisno(-2) - dz_soisno(-1)
        else if(snowdp>0.18 .and. snowdp<=0.29)then
           snl = -4
           dz_soisno(-3) = 0.02
           dz_soisno(-2) = 0.05
           dz_soisno(-1) = (snowdp - dz_soisno(-3) - dz_soisno(-2))/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.29 .and. snowdp<=0.41)then
           snl = -4
           dz_soisno(-3) = 0.02
           dz_soisno(-2) = 0.05
           dz_soisno(-1) = 0.11
           dz_soisno( 0) = snowdp - dz_soisno(-3) - dz_soisno(-2) - dz_soisno(-1)
        else if(snowdp>0.41 .and. snowdp<=0.64)then
           snl = -5
           dz_soisno(-4) = 0.02
           dz_soisno(-3) = 0.05
           dz_soisno(-2) = 0.11
           dz_soisno(-1) = (snowdp - dz_soisno(-4) - dz_soisno(-3) - dz_soisno(-2))/2.
           dz_soisno( 0) = dz_soisno(-1)
        else if(snowdp>0.64)then
           snl = -5
           dz_soisno(-4) = 0.02
           dz_soisno(-3) = 0.05
           dz_soisno(-2) = 0.11
           dz_soisno(-1) = 0.23
           dz_soisno( 0) = snowdp - dz_soisno(-4) - dz_soisno(-3) - dz_soisno(-2) - dz_soisno(-1)
        endif

        zi = 0.
        do i = 0, snl+1, -1
           z_soisno(i) = zi - dz_soisno(i)/2.
           zi = -zi-dz_soisno(i)
        enddo
     endif

  endif

END SUBROUTINE snow_ini
!-----------------------------------------------------------------------
! EOP


SUBROUTINE polint(xa,ya,n,x,y)

! Given arrays xa and ya, each of length n, and gi
! value y, and an error estimate dy. If P (x) is the p
! P (xa(i)) = ya(i), i = 1, . . . , n, then the returned value
! (from: "Numerical Recipes")

  use precision
  implicit none
  integer n,NMAX
  real(r8) dy,x,y,xa(n),ya(n)
  parameter (NMAX=10)      !Largest anticipated val
  integer i,m,ns
  real(r8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))

  do i=1,n       !Here we find the index ns of the closest table entry,
     dift=abs(x-xa(i))
     if(dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)  !and initialize the tableau of c's and d's.
     d(i)=ya(i)
  enddo

  y=ya(ns)       !This is the initial approximation to y.
  ns=ns-1

  do m=1,n-1  !For each column of the tableau,
     do i=1,n-m   !we loop over the current c's and d's and update them.
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) print*, 'failure in polint'  !two input xa's are identical.
        den=w/den
        d(i)=hp*den                                !here the c's and d's are updated.
        c(i)=ho*den
     enddo
     if(2*ns.lt.n-m)then  !After each column in the tableau is completed, we decide
        dy=c(ns+1)        !which correction, c or d, we want to add to our accumulating
     else                 !value of y, i.e., which path to take through
        dy=d(ns)          !the tableau-forking up or down. We do this in such a
        ns=ns-1           !way as to take the most "straight line" route through the
     endif                !tableau to its apex, updating ns accordingly to keep track
     y=y+dy               !of where we are. This route keeps the partial approximations
  enddo                   !centered (insofar as possible) on the target x. T he
                          !last dy added is thus the error indication.
END SUBROUTINE polint
!-----------------------------------------------------------------------
! EOP

