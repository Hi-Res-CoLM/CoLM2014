
 subroutine netsolar (idate,dlon,deltim,&
                      itypwat,sigf,albg,albv,alb,ssun,ssha,&
                      forc_sols,forc_soll,forc_solsd,forc_solld,&
                      parsun,parsha,sabvsun,sabvsha,sabg,sabvg,sr,&
                      solvd,solvi,solnd,solni,srvd,srvi,srnd,srni,&
                      solvdln,solviln,solndln,solniln,srvdln,srviln,srndln,srniln)
!=======================================================================
! Net solar absorbed by surface
! Original author : Yongjiu Dai, 09/15/1999; 09/11/2001
!=======================================================================

  use precision
  use MOD_TimeInvariants, only: spval
  use timemanager, only: isgreenwich
  implicit none

! Dummy argument
  integer,  INTENT(in) :: idate(3) ! model time
  integer,  INTENT(in) :: itypwat  ! land water type (99-sea)
  
  real(r8), INTENT(in) :: dlon     ! logitude in radians
  real(r8), INTENT(in) :: deltim   ! seconds in a time step [second]
 
  real(r8), dimension(1:2,1:2), INTENT(in) :: &
        albg,   &! albedo, ground [-]
        albv,   &! albedo, vegetation [-]
        alb,    &! averaged albedo [-]
        ssun,   &! sunlit canopy absorption for solar radiation
        ssha     ! shaded canopy absorption for solar radiation

  real(r8), INTENT(in) :: &
        sigf,   &! fraction of veg cover, excluding snow-buried veg [-]
        forc_sols,   &! atm vis direct beam solar rad onto srf [W/m2]
        forc_soll,   &! atm nir direct beam solar rad onto srf [W/m2]
        forc_solsd,  &! atm vis diffuse solar rad onto srf [W/m2]
        forc_solld    ! atm nir diffuse solar rad onto srf [W/m2]

  real(r8), INTENT(out) :: &
        parsun, &! PAR absorbed by sunlit vegetation [W/m2]
        parsha, &! PAR absorbed by shaded vegetation [W/m2]
        sabvsun,&! solar absorbed by sunlit vegetation [W/m2]
        sabvsha,&! solar absorbed by shaded vegetation [W/m2]
        sabg,   &! solar absorbed by ground  [W/m2]
        sabvg,  &! solar absorbed by ground + vegetation [W/m2]
        sr,     &! total reflected solar radiation (W/m2)
        solvd,  &! incident direct beam vis solar radiation (W/m2)
        solvi,  &! incident diffuse beam vis solar radiation (W/m2)
        solnd,  &! incident direct beam nir solar radiation (W/m2)
        solni,  &! incident diffuse beam nir solar radiation (W/m2)
        srvd,   &! reflected direct beam vis solar radiation (W/m2)
        srvi,   &! reflected diffuse beam vis solar radiation (W/m2)
        srnd,   &! reflected direct beam nir solar radiation (W/m2)
        srni,   &! reflected diffuse beam nir solar radiation (W/m2)
        solvdln,&! incident direct beam vis solar radiation at local noon(W/m2)
        solviln,&! incident diffuse beam vis solar radiation at local noon(W/m2)
        solndln,&! incident direct beam nir solar radiation at local noon(W/m2)
        solniln,&! incident diffuse beam nir solar radiation at local noon(W/m2)
        srvdln, &! reflected direct beam vis solar radiation at local noon(W/m2)
        srviln, &! reflected diffuse beam vis solar radiation at local noon(W/m2)
        srndln, &! reflected direct beam nir solar radiation at local noon(W/m2)
        srniln   ! reflected diffuse beam nir solar radiation at local noon(W/m2)

! ----------------local variables ---------------------------------
   integer  :: local_secs
   real(r8) :: pi, radpsec

!=======================================================================

        sabvsun = 0.
        sabvsha = 0.
        parsun = 0.
        parsha = 0.

        sabg = 0.
        sabvg = 0.

        if(forc_sols+forc_soll+forc_solsd+forc_solld>0.)then
           if(itypwat<4)then        !non lake and ocean
           ! Radiative fluxes onto surface
              parsun  = ssun(1,1)*forc_sols + ssun(1,2)*forc_solsd
              parsha  = ssha(1,1)*forc_sols + ssha(1,2)*forc_solsd
              sabvsun = forc_sols*ssun(1,1) + forc_soll*ssun(2,1) &
                 + forc_solsd*ssun(1,2) + forc_solld*ssun(2,2)
              sabvsha = forc_sols*ssha(1,1) + forc_soll*ssha(2,1) &
                 + forc_solsd*ssha(1,2) + forc_solld*ssha(2,2)
              sabvsun = sigf*sabvsun
              sabvsha = sigf*sabvsha
              sabvg = forc_sols *(1.-alb(1,1)) + forc_soll *(1.-alb(2,1)) &
                 + forc_solsd*(1.-alb(1,2)) + forc_solld*(1.-alb(2,2))
              sabg = sabvg - sabvsun - sabvsha
           else                     !lake or ocean
              sabvg = forc_sols *(1.-alb(1,1)) + forc_soll *(1.-alb(2,1)) &
                       + forc_solsd*(1.-alb(1,2)) + forc_solld*(1.-alb(2,2))
              sabg = sabvg
           endif
        endif

        solvd = forc_sols
        solvi = forc_solsd
        solnd = forc_soll
        solni = forc_solld
        srvd  = solvd*alb(1,1)
        srvi  = solvi*alb(1,2)
        srnd  = solnd*alb(2,1)
        srni  = solni*alb(2,2)
        sr    = srvd + srvi + srnd + srni

      ! calculate the local secs
        pi = 4.*atan(1.)
        radpsec = pi/12./3600.
        if ( isgreenwich ) then
           local_secs = idate(3) + nint((dlon/radpsec)/deltim)*deltim
           local_secs = mod(local_secs,86400)
        else 
           local_secs = idate(3)
        end if
        
        if (local_secs == 86400/2) then
           solvdln = forc_sols
           solviln = forc_solsd
           solndln = forc_soll
           solniln = forc_solld
           srvdln  = solvdln*alb(1,1)
           srviln  = solviln*alb(1,2)
           srndln  = solndln*alb(2,1)
           srniln  = solniln*alb(2,2)
        else
           solvdln = spval
           solviln = spval
           solndln = spval
           solniln = spval
           srvdln  = spval
           srviln  = spval
           srndln  = spval
           srniln  = spval
        end if


 end subroutine netsolar
