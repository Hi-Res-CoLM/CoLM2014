
 subroutine snowfraction (fveg,z0m,zlnd,scv,snowdp,wt,sigf,fsno)

!=======================================================================
! Provide snow cover fraction
!
! Original author : Yongjiu Dai, /09/1999/, /04/2014/
!=======================================================================

  use precision
  implicit none

! dummy arguments
  real(r8), INTENT(in) :: scv    ! snow water equivalent [mm or kg/m3]
  real(r8), INTENT(in) :: snowdp ! snow depth [m]
  real(r8), INTENT(in) :: z0m    ! aerodynamic roughness length [m]
  real(r8), INTENT(in) :: zlnd   ! aerodynamic roughness length over soil surface [m]
  real(r8), INTENT(in) :: fveg   ! fractional vegetation cover [-]

  real(r8), INTENT(out) :: wt    ! fraction of vegetation covered with snow [-]
  real(r8), INTENT(out) :: sigf  ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(out) :: fsno  ! fraction of soil covered by snow [-]

  real(r8) :: fmelt              ! dimensionless metling factor
  real(r8), parameter :: m = 1.0 ! the value of m used in CLM4.5 is 1.0.
                                 ! while the value of m given by Niu et al (2007) is 1.6 
                                 ! while Niu (2012) suggested 3.0
!-----------------------------------------------------------------------
      if(fveg > 0.001) then

! Fraction of vegetation buried (covered) by snow
         wt = 0.1*snowdp/z0m
         wt = wt/(1.+wt)

! Fraction of vegetation cover free of snow
         sigf = (1.-wt)*fveg
      else
         wt = 0.
         sigf = 0.
      endif

      if(sigf < 0.001) sigf = 0.
      if(sigf > 0.999) sigf = 1.

! Fraction of soil covered by snow
      fsno = 0.0
      if(snowdp > 0.) then
         fmelt = (scv/snowdp/100.) ** m
         fsno  = tanh(snowdp/(2.5 * zlnd * fmelt))
      end if

 end subroutine snowfraction
