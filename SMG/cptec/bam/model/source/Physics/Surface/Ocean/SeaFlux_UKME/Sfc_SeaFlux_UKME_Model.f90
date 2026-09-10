MODULE Sfc_SeaFlux_UKME_Model
    USE PhysicalFunctions , ONLY : fpvs

  IMPLICIT NONE
SAVE

  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  !*L------------------COMDECK C_EPSLON-----------------------------------
  ! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

  REAL(KIND=r8), PARAMETER :: Epsilonk   = 0.62198_r8
  REAL(KIND=r8), PARAMETER :: C_Virtual = 1.0_r8/Epsilonk-1.0_r8
  REAL(kind=r8), PARAMETER :: con_rd      =2.8705e+2_r8 ! gas constant air    (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_rv      =4.6150e+2_r8 ! gas constant H2O    (J/kg/K)
  REAL(kind=r8), PARAMETER:: EPS      =con_rd/con_rv
  REAL(kind=r8), PARAMETER:: EPSM1    =con_rd/con_rv-1.0_r8

  !*----------------------------------------------------------------------
  !*L------------------COMDECK C_O_DG_C-----------------------------------
  ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
  ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
  ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

  REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
  REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
  REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

  !*----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  !  Local parameters and other physical constants
  !-----------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: one_minus_epsilon = 1.0_r8 -epsilonk
  !      One minus the ratio of the molecular weights of water and dry
  !      air
  REAL(KIND=r8)   , PARAMETER   :: hltm   = 2.52e6_r8            !  latent heat of vaporization (J kg^-1)

  REAL(KIND=r8), PARAMETER ::                                                &
       & T_low  = 183.15_r8                                                  &
                                ! Lowest temperature for which look_up table of
                                ! saturation water vapour presssure is valid (K)
       &,T_high = 338.15_r8                                                  &
                                ! Highest temperature for which look_up table of
                                ! saturation water vapour presssure is valid (K)
       &,delta_T = 0.1_r8    ! temperature increment of the look-up table
  ! of saturation vapour pressures.

  INTEGER, PARAMETER ::                                             &
       &  N = ((T_high - T_low + (delta_T*0.5_r8))/delta_T) + 1.0_r8
  !      Gives  N=1551, size of the look-up table of saturation water
  !      vapour pressure.

  INTEGER IES       ! LOOP COUNTER FOR DATA STATEMENT LOOK-UP TABLE

  !  Local dynamic arrays-------------------------------------------------
  REAL(KIND=r8) ES(0:N+1)    ! TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
  !                       - SET BY DATA STATEMENT CALCULATED FROM THE
  !                         GOFF-GRATCH FORMULAE AS TAKEN FROM LANDOLT-
  !                         BORNSTEIN, 1987 NUMERICAL DATA AND FUNCTIONAL
  !                         RELATIONSHIPS IN SCIENCE AND TECHNOLOGY.
  !                         GROUP V/ VOL 4B METEOROLOGY. PHYSICAL AND
  !                         CHEMICAL PROPERTIES OF AIR, P35
  !
  INTEGER,PARAMETER :: NDIV = 6       ! number of particle size divisions
  REAL(KIND=r8), PARAMETER ::  CHARNOCK=.014_r8
  INTEGER , PARAMETER ::ISeaZ0T =1
  !
  !     Options for marine boundary layers
  INTEGER, PARAMETER :: Fixed_Z0T = 0
  !       Stanard flixed value of thermal roughness length over sea
  INTEGER, PARAMETER :: SurfDivZ0T = 1
  !       Thermal roughness length over sea defined from surface
  !       divergence theory

  REAL(KIND=r8), PARAMETER :: SeaSalinityFactor=1.0_r8
  !                                  ! Factor allowing for the
  !                                  ! effect of the salinity of
  !                                  ! sea water on the evaporative
  !                                  ! flux.
  LOGICAL,PARAMETER  :: lq_mix_bl =.FALSE. !TRUE ! IN TRUE if mixing ratios used in boundary layer code
  LOGICAL,PARAMETER  :: L_DUST    =.FALSE.! IN switch for mineral dust
  LOGICAL,PARAMETER  :: SFME      =.TRUE. ! IN STASH flag for wind mixing energy flux
  LOGICAL,PARAMETER  :: QSAT_opt  =.FALSE.

CONTAINS

  SUBROUTINE InitSfc_SeaFlux_UKME_Model()
    IMPLICIT NONE

    CALL qsat_data()

  END SUBROUTINE InitSfc_SeaFlux_UKME_Model




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
!!!   SUBROUTINE SF_EXCH------------------------------------------------
!!!
!!!  Purpose: Calculate coefficients of turbulent exchange between
!!!           the surface and the lowest atmospheric layer, and
!!!           "explicit" fluxes between the surface and this layer.
!!!
!!!  Suitable for Single Column use.
!!!
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version 2, dated 18/1/90.
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Project task:
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!                 See especially sub-section (ix).
!!!
!!!---------------------------------------------------------------------

  ! Arguments :-
  SUBROUTINE SF_EXCH (&
       nCols  , &
       kMAx   , &
       CF     , &
       QCF    , &
       QCL    , &
       Q      , &
       T      , &
       gu     , &
       gv     , &
       prsi   , &
       prsl   , &
       phii   , &
       phil   , &
       sinclt , &
       RADNET , &
       TI     , &
       ZH     ,&
       mskant   , &
       iMask      , &
       TSTAR_LAND , &
       TSTAR_SSI , &
       TSTAR_SEA ,&
       TSTAR_SICE ,&
       rmi       , & 
       rhi       , &
       evap , &
       sens , &
       ustar      , &
       Z0MSEA &
       )
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols ! IN Number of X points?
    INTEGER, INTENT(IN   ) :: kMAx
    !       Factor allowing for the effect of the salinity of
    !       sea water on the evaporative flux.
    REAL(KIND=r8)   , INTENT(IN   ) :: CF    (nCols,kMAx)  ! IN Cloud fraction (decimal).
    REAL(KIND=r8)   , INTENT(IN   ) :: QCF   (nCols,kMAx)    ! IN Cloud ice (kg per kg air)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCL   (nCols,kMAx)    ! IN Cloud liquid water
    REAL(KIND=r8)   , INTENT(IN   ) :: Q     (nCols,kMAx)    ! IN specific humidity
    REAL(KIND=r8)   , INTENT(IN   ) :: T     (nCols,kMAx)    ! IN temperature
    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: prsi     (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8), INTENT(IN   ) :: prsl     (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8), INTENT(IN   ) :: phii     (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8), INTENT(IN   ) :: phil     (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: RADNET(nCols) ! IN Net surface radiation (W/m2) positive downwards
    REAL(KIND=r8)   , INTENT(IN   ) :: TI    (nCols) ! IN Temperature of sea-ice surface layer (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: ZH    (nCols) ! IN Depth of boundary layer (m). boundary layer (metres).

    !       Switch for the treatment of the thermal roughness
    !       length over the sea.
    !       


    INTEGER(KIND=i8)   , INTENT(IN   )  :: mskant  (nCols)  ! IN Fraction of gridbox which is sea-ice.
    INTEGER(KIND=i8)   , INTENT(IN   )  :: iMask     (nCols)  ! IN Land fraction on all tiles.
    !                            !    Note, if the grid used is staggered in
    !                            !    the vertical, Z1_UV and Z1_TQ can be
    !                            !    different.

    REAL(KIND=r8)   , INTENT(IN   )  :: TSTAR_LAND(nCols) ! IN Land mean surface temperature (K).
    REAL(KIND=r8)   , INTENT(IN   )  :: TSTAR_SSI (nCols) ! IN Mean sea surface temperature (K).
    REAL(KIND=r8)   , INTENT(IN   )  :: TSTAR_SEA (nCols) ! IN Open sea surface temperature (K).
    REAL(KIND=r8)   , INTENT(IN   )  :: TSTAR_SICE(nCols) ! IN Sea-ice surface temperature (K).
    REAL(KIND=r8)   , INTENT(OUT  )  :: rmi   (nCols) 
    REAL(KIND=r8)   , INTENT(OUT  )  :: rhi   (nCols) 
    REAL(KIND=r8)   , INTENT(OUT  )  :: evap  (nCols) 
    REAL(KIND=r8)   , INTENT(OUT  )  :: sens  (nCols) 
    REAL(KIND=r8)   , INTENT(OUT  )  :: ustar (nCols) 

    REAL(KIND=r8)   , INTENT(INOUT)  :: Z0MSEA    (nCols) ! INOUT Sea-surface roughness length for momentum (m).  F617.

    REAL(KIND=r8)          :: Z1_TQ      (nCols)  ! IN Height of lowest tq level (m).
    REAL(KIND=r8)     :: Z1_UV     (nCols) ! IN Height of lowest uv level (m).
    REAL(KIND=r8)   :: PSTAR     (nCols) ! IN Surface pressure (Pascals).
    REAL(KIND=r8) :: VSHR_SSI  (nCols) ! IN Mag. of mean sea sfc-to-lowest-level wind shear

    REAL(KIND=r8)    ::  TL_1(nCols) ! IN Liquid/frozen water temperature for lowest atmospheric layer (K).

    REAL(KIND=r8)    :: BQ_1(nCols)
    ! IN A buoyancy parameter for lowest atm
    !                            !    level ("beta-q twiddle").
    REAL(KIND=r8)   :: BT_1(nCols)
    ! IN A buoyancy parameter for lowest atm
    !                            !    level ("beta-T twiddle").
    REAL(KIND=r8)     :: QW_1(nCols) ! IN Total water content of lowest
    !                            !    atmospheric layer (kg per kg air).


    REAL(KIND=r8)     :: RHOSTAR(nCols) ! OUT Surface air density

    REAL(KIND=r8)    :: FLANDG     (nCols)  ! IN Land fraction on all tiles.

    REAL(KIND=r8)     :: ICE_FRACT  (nCols)  ! IN Fraction of gridbox which is sea-ice.




    !   Define local storage.

    !   (a) Workspace.

    REAL(KIND=r8)   :: QS1(nCols)        ! Sat. specific humidity
    !                                  ! qsat(TL_1,PSTAR)
    REAL(KIND=r8)   :: QSTAR_SEA(nCols)
    ! Surface saturated sp humidity
    REAL(KIND=r8)   :: QSTAR_ICE(nCols)! Surface saturated sp humidity
    !      REAL(KIND=r8)   :: PSTAR_LAND(LAND_PTS)! Surface pressure for land points.
    REAL(KIND=r8)   :: Z0H_SEA(nCols)! Roughness length for heat and
    !                                  ! moisture transport

    INTEGER   :: NSICE                         ! Number of sea-ice points.
    INTEGER   :: SICE_INDEX (nCols,2)! Index of sea-ice points
    REAL(KIND=r8)      :: Z1_TQ_SEA(nCols)! Height of lowest model level
    !                            ! relative to sea.
    REAL(KIND=r8)      :: Z0_ICE(nCols) ! Roughness length.
    REAL(KIND=r8)      :: Z0_MIZ(nCols) ! Roughness length.
    REAL(KIND=r8)      :: RIB_SEA(nCols)! Bulk Richardson number
    REAL(KIND=r8)      :: RIB_ICE(nCols)! Bulk Richardson number
    REAL(KIND=r8)      :: DB_SEA(nCols) ! Buoyancy difference for sea points
    REAL(KIND=r8)      :: DB_ICE(nCols) ! Buoyancy difference for sea ice

    REAL(KIND=r8)      :: BT     (nCols,kMAx)
    REAL(KIND=r8)      :: BQ     (nCols,kMAx)
    REAL(KIND=r8)      :: P     (nCols,kMAx)
    REAL(KIND=r8)      :: BT_CLD (nCols,kMAx)
    REAL(KIND=r8)      :: BQ_CLD (nCols,kMAx)
    REAL(KIND=r8)      :: BT_GB  (nCols,kMAx)
    REAL(KIND=r8)      :: BQ_GB  (nCols,kMAx)
    REAL(KIND=r8)      :: A_QS   (nCols,kMAx)
    REAL(KIND=r8)      :: A_DQSDT(nCols,kMAx)
    REAL(KIND=r8)      :: DQSDT  (nCols,kMAx)
    REAL(KIND=r8)      :: QW(nCols, kMAx)! OUT Total water content
    REAL(KIND=r8)      :: TL(nCols, kMAx)! OUT Ice/liquid water temperature
    REAL(KIND=r8)      :: U_0_P(nCols)  ! IN W'ly component of surface current (m/s). P grid
    REAL(KIND=r8)      :: V_0_P(nCols)  ! IN S'ly component of surface current (m/s). P grid
    REAL(KIND=r8)      :: U_1_P(nCols)   ! IN U_1 on P-grid.
    REAL(KIND=r8)      :: V_1_P(nCols)   ! IN V_1 on P-grid.

    REAL(KIND=r8)      :: ALPHA1_SICE(nCols)
    REAL(KIND=r8)      :: ASHTF(nCols)
    REAL(KIND=r8)      :: E_SEA(nCols)
    REAL(KIND=r8)      :: FTL_ICE(nCols)
    REAL(KIND=r8)      :: FTL_1(nCols)
    REAL(KIND=r8)      :: H_SEA(nCols)
    REAL(KIND=r8)      :: RHO_ARESIST(nCols)
    !                            ! OUT RHOSTAR*CD_STD*VSHR  for SCYCLE
    !      REAL(KIND=r8)      :: RHO_ARESIST_LAND(nCols)
    !                            ! Land mean of rho_aresist_tile
    REAL(KIND=r8)      :: ARESIST(nCols)
    !                            ! OUT 1/(CD_STD*VSHR)      for SCYCLE
    REAL(KIND=r8)      :: RESIST_B(nCols)
    !                            ! OUT (1/CH-1/CD_STD)/VSHR for SCYCLE
    !      REAL(KIND=r8)      :: RHOKM_LAND(nCols,nCols)
    !                            ! OUT For land momentum. NB: This is output
    REAL(KIND=r8)      :: RHOKH_1_SICE(nCols)
    !                            ! OUT Surface exchange coefficient for sea
    !                            !     or sea-ice.

    ! OUT "Explicit" surface flux of TL = H/CP.
    !                            !     (sensible heat / CP). grid-box mean
    REAL(KIND=r8)      :: FQW_1(nCols)
    ! OUT "Explicit" surface flux of QW (i.e.
    !                            !     evaporation), on P-grid (kg/m2/s).
    !                            !     for whole grid-box
    REAL(KIND=r8)      :: FQW_ICE(nCols)
    !                            ! OUT GBM FQW_1 for sea-ice.
    REAL(KIND=r8)      :: RHOKM_SSI(nCols)
    !                            ! OUT For mean sea mom. NB: This is output
    !                            !     on UV-grid, but with the first and
    !                            !     last set to "missing data".
    REAL(KIND=r8)      :: RHOKPM_SICE(nCols)
    !                            ! OUT Mixing coefficient for sea-ice.
    REAL(KIND=r8)      :: CD_SSI(nCols)
    !                            ! OUT Bulk transfer coefficient for
    !                            !      momentum over sea mean.
    REAL(KIND=r8)      :: CH_SSI(nCols)
    !                            ! OUT Bulk transfer coefficient for heat
    !                            !    and/or moisture over sea mean.
    REAL(KIND=r8)      :: CD_LAND(nCols)
    ! Bulk transfer coefficient for
    !                                  !      momentum over land.
    REAL(KIND=r8)      :: CD(nCols)
    ! OUT Bulk transfer coefficient for
    !                            !      momentum.
    REAL(KIND=r8)      :: CH(nCols)
    ! OUT Bulk transfer coefficient for heat
    !                            !     and/or moisture.
    REAL(KIND=r8)      :: CD_STD_DUST(nCols)
    !OUT Bulk transfer coef. for
    !                             ! momentum, excluding orographic effects
    !  Workspace for sea and sea-ice leads
    REAL(KIND=r8)      ::    CD_ICE(nCols)
    ! Drag coefficient

    REAL(KIND=r8)      :: CD_SEA(nCols)
    ! Drag coefficient
    REAL(KIND=r8)      :: CH_SEA(nCols)
    ! Transfer coefficient for heat and
    !                                  ! moisture
    REAL(KIND=r8)      :: CD_MIZ(nCols)
    ! Drag coefficient
    REAL(KIND=r8)      :: CH_ICE(nCols)
    ! Transfer coefficient for heat and
    !                                  ! moisture
    REAL(KIND=r8)      :: CH_MIZ(nCols)
    ! Transfer coefficient for heat and
    !                                  ! moisture
    REAL(KIND=r8)      :: RECIP_L_MO_SEA(nCols)
    !                            ! OUT Reciprocal of the Monin-Obukhov
    !                            !     length for sea points (m^-1)
    REAL(KIND=r8)      :: RECIP_L_MO_MIZ(nCols)
    !                                  ! Reciprocal of the Monin-Obukhov
    !                                  ! length for marginal sea ice (m^-1).
    REAL(KIND=r8)      :: RECIP_L_MO_ICE(nCols)
    !                                  ! Reciprocal of the Monin-Obukhov
    !                                  ! length for sea ice (m^-1).

    REAL(KIND=r8)      :: V_S_SEA(nCols)    ! Surface layer scaling velocity
    !                                  ! for sea points (m/s).
    REAL(KIND=r8)      :: V_S_ICE(nCols)
    ! Surface layer scaling velocity
    !                                  ! for sea ice (m/s).
    REAL(KIND=r8)      :: V_S_MIZ(nCols)
    ! Surface layer scaling velocity
    !                                  ! for marginal sea ice (m/s).
    REAL(KIND=r8)      :: Q1_SD(nCols)
    ! OUT Standard deviation of turbulent
    !                            !     fluctuations of surface layer
    !                            !     specific humidity (kg/kg).
    REAL(KIND=r8)      :: T1_SD(nCols)
    ! OUT Standard deviation of turbulent
    !                            !     fluctuations of surface layer
    !                            !     temperature (K).
    REAL(KIND=r8)      :: speedm(nCols)
    REAL(KIND=r8)      :: FME(nCols)
    ! OUT Wind mixing energy flux (Watts/sq m).
    REAL(KIND=r8)      :: VSHR(nCols)
    REAL(KIND=r8)     :: TAU ! Magnitude of surface wind stress over sea.

    REAL(KIND=r8)      :: TOL_USTR_N  ! Tolerance for USTR_N
    REAL(KIND=r8)      :: TOL_USTR_L  ! Tolerance for USTR_L (see below)
    REAL(KIND=r8)      :: USTR_L      ! Low-wind estimate of friction velocity
    REAL(KIND=r8)      :: USTR_N      ! Neutral estimate of friction velocity
    REAL(KIND=r8)      :: R_B_DUST(nCols,NDIV)
    !OUT surf layer res for dust
    REAL(KIND=r8)      :: qss(nCols)

    INTEGER   :: I,K
    INTEGER   :: JITS,IDIV
    !!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
    ! Sea ice parameters
    ! Z0FSEA = roughness length for free convective heat and moisture
    !          transport over the sea (m).
    !          DUMMY VARIABLE - Only used in 7A boundary layer scheme
    ! Z0HSEA = roughness length for heat and moisture transport
    !          over the sea (m).
    ! Z0MIZ  = roughness length for heat, moisture and momentum over
    !          the Marginal Ice Zone (m).
    ! Z0SICE = roughness length for heat, moisture and momentum over
    !          sea-ice (m).
    REAL(KIND=r8), PARAMETER :: Z0HSEA=1.14e-5_r8
    REAL(KIND=r8), PARAMETER :: Z0MIZ=1.5e-5_r8
    REAL(KIND=r8), PARAMETER :: Z0SICE= 0.002e0_r8!!0.0002e0_r8! 

    !!*----------------------------------------------------------------------
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    REAL(KIND=r8), PARAMETER  :: CP     = 1005.0_r8


    REAL(KIND=r8), PARAMETER  :: R      = 287.05_r8

    !!----------------------------------------------------------------------
    ! C_VKMAN start
    REAL(KIND=r8),PARAMETER:: VKMAN=0.4_r8 ! Von Karman's constant
    !-----------------------------------------------------------------------
    !! Calculate total water content, QW and Liquid water temperature, TL
    !-----------------------------------------------------------------------
    ! Derived local parameters.
    !*----------------------------------------------------------------------
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8

    REAL(KIND=r8),PARAMETER:: RHOSEA = 1000.0_r8 ! density of sea water (kg/m3)

    REAL(KIND=r8),PARAMETER :: LCRCP=LC/CP  ! Evaporation-to-dT conversion factor.
    REAL(KIND=r8),PARAMETER :: LS=LF+LC       ! Latent heat of sublimation.
    REAL(KIND=r8),PARAMETER :: LSRCP=LS/CP   ! Sublimation-to-dT conversion factor.

    REAL(KIND=r8) :: USHEAR,VSHEAR,VSHR2
    rmi  =0.0_r8
    rhi  =0.0_r8
    evap =0.0_r8
    sens =0.0_r8
    ustar=0.0_r8
   ! REAL(KIND=r8), INTENT(IN   ) :: prsi     (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
   ! REAL(KIND=r8), INTENT(IN   ) :: prsl     (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa

    DO K=1,kMAx
       DO I=1,nCols
          IF(mskant(i) == 1_i8)THEN
             P(i,k)  = prsl     (i,k) 
             QW(I,K) = Q(I,K) + QCL(I,K)       + QCF(I,K)      ! P243.10
             TL(I,K) = T(I,K) - LCRCP*QCL(I,K) - LSRCP*QCF(I,K)! P243.9
          END IF
       ENDDO
    ENDDO

    ICE_FRACT=0.0_r8
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (TSTAR_SEA(i) < 0.0_r8 .AND. ABS(TSTAR_SEA(i)) < tice+0.01_r8) THEN
             ICE_FRACT(i)=1.0_r8
          ELSE IF (TSTAR_SEA(i) < 0.0_r8 .AND. ABS(TSTAR_SEA(i)) > tice+0.01_r8) THEN
             ICE_FRACT(i)=0.0_r8
          END IF
       END IF
    END DO

    FLANDG=0.0_r8
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN
          IF(IMASK(i) > 1_i8)FLANDG(i)=1.0_r8
          QW_1(i)=QW(I,1)
          TL_1(i)=TL(I,1)
          Z1_TQ(i)=0.5_r8*MAX((phii(i,2) - phii(i,1)),2.0_r8) !rbyg*T(i,1)
          Z1_UV(i)=Z1_TQ(i)
          Pstar(i) =prsi     (i,1)
          !speedm(i)=SQRT(gu(i,1)**2+gv(i,1)**2)*sinclt(i)
          !USHEAR = ((gu(I,1) - 0.01_r8*ABS(gu(I,1)))*sinclt(i))
          !VSHEAR = ((gv(I,1) - 0.01_r8*ABS(gv(I,1)))*sinclt(i))      
        
           U_1_P(I) = ((gu(I,2)/sinclt(i))/MAX(ABS(gu(I,2)/sinclt(i)),1.0e-12_r8)) * ABS(gu(I,1)/sinclt(i))
           V_1_P(I) = ((gv(I,2)/sinclt(i))/MAX(ABS(gv(I,2)/sinclt(i)),1.0e-12_r8)) * ABS(gv(I,1)/sinclt(i))

           U_0_P(I) = ((gu(I,1)/sinclt(i))/MAX(ABS(gu(I,1)/sinclt(i)),1.0e-12_r8)) * ((ABS((gu(I,1)/sinclt(i))))**0.25 )
           V_0_P(I) = ((gv(I,1)/sinclt(i))/MAX(ABS(gv(I,1)/sinclt(i)),1.0e-12_r8)) * ((ABS((gv(I,1)/sinclt(i))))**0.25 )
           USHEAR = U_1_P(I) - U_0_P(I)
           VSHEAR = V_1_P(I) - V_0_P(I)
 
          speedm(i)=SQRT((gu(i,1)/sinclt(i))**2+(gv(i,1)/sinclt(i))**2)
          !USHEAR = (gu(I,1)/sinclt(i)) - 0.01_r8*ABS((gu(I,1)/sinclt(i)))
          !VSHEAR = (gv(I,1)/sinclt(i)) - 0.01_r8*ABS((gv(I,1)/sinclt(i)))

          !USHEAR = ((gu(I,1) - ((gu(I,1)/MAX(ABS(gu(I,1)),0.1_r8))*sqrt(ABS(gu(I,1)))))*sinclt(i))/Z1_UV(i)
          !VSHEAR = ((gv(I,1) - ((gv(I,1)/MAX(ABS(gv(I,1)),0.1_r8))*sqrt(ABS(gv(I,1)))))*sinclt(i))/Z1_UV(i)
 
          !USHEAR = VKMAN*(gu(I,1) - ((gu(I,1)/MAX(ABS(gu(I,1)),0.1_r8))*sqrt(ABS(gu(I,1)))))* &
          !   sinclt(i)/(log(Z1_UV(i)/Z0MSEA(i)))
          !VSHEAR = VKMAN*(gv(I,1) - ((gv(I,1)/MAX(ABS(gv(I,1)),0.1_r8))*sqrt(ABS(gv(I,1)))))* &
          !   sinclt(i)/(log(Z1_UV(i)/Z0MSEA(i)))
          VSHR2 = MAX (1.0E-6_r8 , USHEAR*USHEAR + VSHEAR*VSHEAR)
          VSHR_SSI(I) = SQRT(VSHR2)
       END IF
    END DO

    CALL BOUY_TQ (&
         nCols                         , & !INTEGER, INTENT(IN   ) :: nCols
         kMAx                          , & !INTEGER, INTENT(IN   ) :: kMAx   
         LQ_MIX_BL                     , & !LOGICAL, INTENT(IN   ) :: LQ_MIX_BL   
         mskant   (1:nCols )                    , &
         P        (1:nCols,1:kMAx), & !REAL, INTENT(IN    )   :: P      (nCols,rows,kMAx)
         T        (1:nCols,1:kMAx), & !REAL, INTENT(IN    )   :: T      (nCols,rows,kMAx)
         Q        (1:nCols,1:kMAx), & !REAL, INTENT(IN    )   :: Q      (nCols,rows,kMAx)
         QCF      (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(IN    )   :: QCF    (nCols,rows,kMAx)
         QCL      (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(IN    )   :: QCL    (nCols,rows,kMAx)
         CF       (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(IN    )   :: CF     (nCols,rows,kMAx)
         BT       (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(out   )   :: BT     (nCols,rows,kMAx)
         BQ       (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(out   )   :: BQ     (nCols,rows,kMAx)
         BT_CLD   (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(out   )   :: BT_CLD (nCols,rows,kMAx)
         BQ_CLD   (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(out   )   :: BQ_CLD (nCols,rows,kMAx)
         BT_GB    (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(out   )   :: BT_GB  (nCols,rows,kMAx)
         BQ_GB    (1:nCols,1:kMAx), & !REAL(KIND=r8)(KIND=r8), INTENT(out   )   :: BQ_GB  (nCols,rows,kMAx)
         A_QS     (1:nCols,1:kMAx), & !REAL(KIND=r8), INTENT(out   )   :: A_QS   (nCols,rows,kMAx)
         A_DQSDT  (1:nCols,1:kMAx), & !REAL(KIND=r8), INTENT(out   )   :: A_DQSDT(nCols,rows,kMAx)
         DQSDT    (1:nCols,1:kMAx)  & !REAL(KIND=r8), INTENT(out   )   :: DQSDT  (nCols,rows,kMAx)
         )
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN
          BQ_1(i)=BQ_GB    (i,1)
          BT_1(i)=BT_GB    (I,1)
       END IF
    END DO

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('SFEXCH  ',3)
    !      ENDIF
    !---------     &,VSHR_SSI(nCols,ROWS)                                        &
    !                            ! IN Mag. of mean sea sfc-to-lowest-level
    !                            !    wind shear
    !--------------------------------------------------------------
    !!  1. Index array for sea-ice
    !-----------------------------------------------------------------------

    NSICE = 0
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( ICE_FRACT(I) >  0.0_r8 .AND. FLANDG(I) <  1.0_r8 ) THEN
             NSICE = NSICE + 1
             SICE_INDEX(NSICE,1) = I
          ENDIF
       END IF
    ENDDO
    !-----------------------------------------------------------------------
    !!  1.1 Calculate height of lowest model level relative to sea.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          Z1_TQ_SEA(I)=Z1_TQ(I)
          !IF(L_CTILE.AND.FLANDG(I) >  0.0.AND.FLANDG(I) <  1.0) Z1_TQ_SEA(I)=Z1_TQ(I)+Z_LAND(I)
          !IF(FLANDG(I) >  0.0.AND.FLANDG(I) <  1.0) Z1_TQ_SEA(I)=Z1_TQ(I)+Z_LAND(I)
       END IF
    ENDDO

    !-----------------------------------------------------------------------
    !!  2.  Calculate QSAT values required later.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          RHOSTAR(I) = PSTAR(I) /( R*(FLANDG(I)*TSTAR_LAND(I)+(1.0_r8-FLANDG(I))*TSTAR_SSI(I)) )
          !                        ... surface air density from ideal gas equation
       END IF
    ENDDO
    ! DEPENDS ON: qsat_mix
    IF(QSAT_opt)THEN
       CALL QSAT_mix(  QS1    (1:nCols)   , &!REAL(KIND=r8), intent(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturationeing processed by qSAT scheme.
            TL_1   (1:nCols)   , &!REAL(KIND=r8), intent(in)  :: T(npnts)      !  Temperature (K).
            PSTAR  (1:nCols)   , &!REAL(KIND=r8), intent(in)  :: P(npnts)      !  Pressure (Pa).
            nCols                , &!Integer, intent(in) :: npnts   !Points (=horizontal dimensions) being processed by qSAT scheme.
            lq_mix_bl           ,&!logical, intent(in)  :: lq_mix  .true. return qsat as a mixing ratio
            mskant (1:nCols)     )                       !  .false. return qsat as a specific humidity
    ELSE
       DO i=1,nCols
          IF(mskant(i) == 1_i8)THEN
             qss(i) = fpvs(TL_1(i))
             QS1(I) = EPS * qss(I) / (PSTAR(I) + EPSM1 * qss(I))
          END IF
       END DO
    END IF
    ! DEPENDS ON: qsat_mix
    IF(QSAT_opt)THEN
       CALL QSAT_mix(  QSTAR_SEA (1:nCols), &!REAL(KIND=r8), intent(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturationeing processed by qSAT scheme.
            TSTAR_SEA (1:nCols), &!REAL(KIND=r8), intent(in)  :: T(npnts)!  Temperature (K).
            PSTAR     (1:nCols), &!REAL(KIND=r8), intent(in)  :: P(npnts)!  Pressure (Pa).
            nCols                , &!Integer, intent(in) :: npnts   !Points (=horizontal dimensions) being processed by qSAT scheme.
            lq_mix_bl            ,&!logical, intent(in)  :: lq_mix  .true. return qsat as a mixing ratio
            mskant  (1:nCols)               )       !  .false. return qsat as a specific humidity
    ELSE
       DO i=1,nCols
          IF(mskant(i) == 1_i8)THEN
          qss(i) = fpvs(TSTAR_SEA(i))
          QSTAR_SEA(I) = EPS * qss(I) / (PSTAR(I) + EPSM1 * qss(I))
          END IF
       END DO
    END IF
    ! DEPENDS ON: qsat_mix
    IF(QSAT_opt)THEN
       CALL QSAT_mix  (QSTAR_ICE (1:nCols), &!REAL(KIND=r8), intent(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturationeing processed by qSAT scheme.
            TSTAR_SICE(1:nCols), &!REAL(KIND=r8), intent(in)  :: T(npnts)!  Temperature (K).
            PSTAR     (1:nCols), &!REAL(KIND=r8), intent(in)  :: P(npnts)!  Pressure (Pa).
            nCols                , &!Integer, intent(in) :: npnts   !Points (=horizontal dimensions) being processed by qSAT scheme.
            lq_mix_bl            ,&!logical, intent(in)  :: lq_mix  .true. return qsat as a mixing ratio
            mskant  (1:nCols)              )        !  .false. return qsat as a specific humidity
    ELSE
       DO i=1,nCols
          IF(mskant(i) == 1_i8)THEN
          qss(i) = fpvs(TSTAR_SICE(i))
          QSTAR_ICE(I) = EPS * qss(I) / (PSTAR(I) + EPSM1 * qss(I))
          END IF
       END DO
    END IF
    !-----------------------------------------------------------------------
    !!  3. Calculation of transfer coefficients and surface layer stability
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !!  3.1 Calculate neutral roughness lengths
    !-----------------------------------------------------------------------

    ! Sea, sea-ice leads, sea-ice and marginal ice zone
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          Z0_MIZ(I) = Z0MIZ
          Z0_ICE(I) = Z0SICE
          RIB_SEA(I) = 0.0_r8
          RIB_ICE(I) = 0.0_r8
          DB_SEA(I) = 0.0_r8
          DB_ICE(I) = 0.0_r8
       END IF
    ENDDO

    IF (ISeaZ0T == SurfDivZ0T) THEN
       !
       !       Composite formulation for thermal roughness lengths,
       !       incoporating the smooth aerodynamic limit for low
       !       wind speeds and a value based on surface divergence
       !       theory for higher wind speeds.
       !
       !       The friction velocity is diagnosed in the surface
       !       transfer scheme, using z0m from the previous time-step.
       !       z0[T,q] is also required but depends on u_* in this
       !       scheme. For practical purposes, it is sufficient to
       !       infer it from u_* determined from z0m, but in general
       !       a given value of z0m corresponds to two values of
       !       u_*, so we need to know whether we are on the low or
       !       high wind speed side of the minimum value of z0m.
       !       If we are on the high side, z0[T,q] will be inversely
       !       proportional to z0m, but on the low side it may follow
       !       this relationship, or be aerodynamically smooth. In
       !       the smooth case we iterate u_* from the expression for
       !       the roughness length and then take the maximum of the
       !       smooth and high-wind expressions for z0[T,q]. An
       !       iteration for the low-wind value of u_*, USTR_L is
       !       carried out. This will converge to the correct limit only
       !       on the low-wind side of the minimum, and the standard
       !       criterion that the gradient of a fixed-point iteration
       !       should be less than 1 in modulus gievs a more precise
       !       boundary, TOL_USTR_L. For consistency with earlier versions
       !       of the modset, hard-wired values are retained for the
       !       operational value of Charnock's parameter. An additional
       !       check is required, since z0m can be large at low or at
       !       high wind-speeds. This is less precise and a fixed
       !       value of 0.07 is used to test USTR_N, which was determined
       !       by inspection of a graph of z0m against u_*: it is
       !       unlikely that this will need to be changed unless
       !       Charnock's constant is greatly altered.
       !
       IF (CHARNOCK == 0.018_r8) THEN
          TOL_USTR_L = 0.055_r8
          TOL_USTR_N = 0.07_r8
       ELSE
          TOL_USTR_L = 0.75_r8*(1.54E-6_r8*G/(2.0_r8*CHARNOCK))**0.33333_r8
          TOL_USTR_N = 0.07_r8
       ENDIF
       !
       DO I=1, nCols
          IF(mskant(i) == 1_i8)THEN

             Z0MSEA(I) = MIN(Z0MSEA(I),0.1_r8 * Z1_UV(I))
             Z0H_SEA(I) = Z0MSEA(I)

             ustar(i)=VSHR_SSI(I)

             !           We need to infer u_* from Z0M.
             IF (VSHR_SSI(I)  > 0.0_r8) THEN
                !             Compute u_* using neutral stratification.
                !             stratification.
                USTR_N = VKMAN * VSHR_SSI(I) / &
                     &          LOG(Z1_UV(I) / Z0MSEA(I) )
                ustar(i)=USTR_N
                !             Compute u_* using low wind approximation.
                USTR_L = 1.54E-06_r8 /  Z0MSEA(I) - 1.0E-05_r8
                !             Since Z0M could be large for low and high u_*, we use
                !             USTR_N as an additional check on the regime.
                IF ( (USTR_N < TOL_USTR_N) .AND.                          &
                     &             (USTR_L < TOL_USTR_L) ) THEN
                   !               Iterate u_* for low winds.
                   DO JITS=1, 5
                      USTR_L=1.54E-06_r8/(Z0MSEA(I)-(CHARNOCK/G)*USTR_L**2)  &
                           &              -1.0E-05_r8
                   ENDDO
                   !               Take the maximum of the smooth and high-wind values.
                   !               A lower limit is imposed on the friction velocity to
                   !               allow for the exceptional case of very low winds: the
                   !               value of 10^-5 is the same as the limit for the momentum
                   !               roughness length.
                   ustar(i)=USTR_L
                   Z0H_SEA(I) = MAX( 2.52E-6_r8/(USTR_L+1.0E-05_r8),2.56E-9_r8/Z0MSEA(I) )
                ELSE
                   !               Take the high-wind value, but limit it to the molecular
                   !               mean free path (we should not hit this limit
                   !               in practice).
                   Z0H_SEA(I) = MAX( 2.56E-9_r8/Z0MSEA(I), 7.0E-08_r8 )
                ENDIF
             ELSE
                IF(ustar(i)==0.0_r8 )ustar(i)=0.000001_r8
             ENDIF
             Z0H_SEA(I) = MIN(Z0H_SEA(I),0.1_r8 * Z1_UV(I))
          END IF
       ENDDO
       !
    ELSE IF (ISeaZ0T == Fixed_Z0T) THEN
       !
       !       Use a fixed thermal roughness length.
       Z0MSEA(1:nCols) = MIN(Z0MSEA(1:nCols),0.1_r8 * Z1_UV(1:nCols))
       Z0H_SEA(1:nCols) = Z0MSEA(1:nCols)
       !Z0H_SEA(1:nCols) = Z0HSEA
       !
    ENDIF
    !
    !-----------------------------------------------------------------------
    !!  3.2 Calculate bulk Richardson number for the lowest model level.
    !-----------------------------------------------------------------------

    ! Sea, sea-ice and sea-ice leads
    ! DEPENDS ON: sf_rib_sea
    CALL SF_RIB_SEA (&
         nCols                  , &!INTEGER, INTENT(IN   ) ::  nCols                     ! IN Number of X points?
         mskant     (1:nCols)   , &!INTEGER, INTENT(IN   ) ::  nCols                     ! IN Number of X points?
         FLANDG     (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  FLANDG     (nCols,ROWS)   ! IN Land fraction on all pts.
         NSICE                  , &!INTEGER, INTENT(IN   ) ::  NSICE                     ! IN Number of sea-ice points.
         SICE_INDEX (1:nCols,1:2) , &!INTEGER, INTENT(IN   ) ::  SICE_INDEX (nCols*ROWS,2) ! IN Index of sea-ice points.
         BQ_1       (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  BQ_1       (nCols,ROWS)   ! IN A buoyancy parameter for lowest atm level. ("beta-q twiddle").
         BT_1       (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  BT_1       (nCols,ROWS)   ! IN A buoyancy parameter for lowest atm level. ("beta-T twiddle").
         QSTAR_ICE  (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  QSTAR_ICE  (nCols,ROWS)   ! IN Surface saturated sp humidity over sea-ice.
         QSTAR_SEA  (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  QSTAR_SEA  (nCols,ROWS)   ! IN Surface saturated sp humidity over sea and sea-ice leads.
         QW_1       (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  QW_1       (nCols,ROWS)   ! IN Total water content of lowest atmospheric layer (kg per kg air)
         TL_1       (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  TL_1       (nCols,ROWS)   ! IN Liquid/frozen water temperature for lowest atmospheric layer (K).
         TSTAR_SICE (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  TSTAR_SICE (nCols,ROWS)   ! IN Surface temperature of sea-ice (K).
         TSTAR_SEA  (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  TSTAR_SEA  (nCols,ROWS)   ! IN Surface temperature of sea and sea-ice leads (K).
         VSHR_SSI   (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  VSHR       (nCols,ROWS)   ! IN Magnitude of surface- to-lowest-level wind shear.
         Z0_ICE     (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0H_ICE    (nCols,ROWS)   ! IN Roughness length for heat and moisture transport over sea-ice (m).
         Z0H_SEA    (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0H_SEA    (nCols,ROWS)   ! IN Roughness length for heat and moisture transport over sea or sea-ice 
         Z0_ICE     (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0M_ICE    (nCols,ROWS)   ! IN Roughness length for momentum over sea-ice (m).
         Z0MSEA     (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0M_SEA    (nCols,ROWS)   ! IN Roughness length for momentum over sea or sea-ice leads (m).
         Z1_TQ_SEA  (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z1_TQ      (nCols,ROWS)   ! IN Height of lowest TQ level (m).
         Z1_UV      (1:nCols)   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z1_UV      (nCols,ROWS)   ! IN Height of lowest UV level (m).
         RIB_SEA    (1:nCols)   , &!REAL(KIND=r8)   , INTENT(OUT  ) ::  RIB_SEA    (nCols,ROWS)   ! OUT Bulk Richardson number for lowest layer over sea or sea-ice leads.
         RIB_ICE    (1:nCols)   , &!REAL(KIND=r8)   , INTENT(OUT  ) ::  RIB_ICE    (nCols,ROWS)   ! OUT Bulk Richardson number for lowest layer over sea-ice.
         DB_SEA     (1:nCols)   , &!REAL(KIND=r8)   , INTENT(OUT  ) ::  DB_SEA     (nCols,ROWS)   ! OUT Buoyancy difference between
         DB_ICE     (1:nCols)     &!REAL(KIND=r8)   , INTENT(OUT  ) ::  DB_ICE     (nCols,ROWS)   ! OUT Buoyancy difference between
         )

    !-----------------------------------------------------------------------
    !!  3.4 Calculate CD, CH via routine FCDCH.
    !!      Note that these are returned as the dimensionless surface
    !!      exchange coefficients.
    !-----------------------------------------------------------------------

    ! Sea-ice
    ! DEPENDS ON: fcdch_sea
    CALL FCDCH_SEA( &
         nCols                     , &!INTEGER, INTENT(IN  ) :: nCols                   ! IN Number of X points?
         mskant        (1:nCols), &!INTEGER, INTENT(IN   ) ::  nCols                     ! IN Number of X points?
         FLANDG        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: FLANDG     (nCols,ROWS) ! IN Land fraction
         DB_ICE        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: DB         (nCols,ROWS) !
         VSHR_SSI      (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: VSHR       (nCols,ROWS) 
         Z0_ICE        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z0M        (nCols,ROWS) 
         Z0_ICE        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z0H        (nCols,ROWS) 
         ZH            (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: ZH         (nCols,i)  
         Z1_UV         (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z1_UV      (nCols,ROWS)  
         Z1_TQ_SEA     (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z1_TQ      (nCols,ROWS)  
         CD_ICE        (1:nCols), &!REAL(KIND=r8)   , INTENT(OUT ) :: CDV        (nCols,ROWS) 
         CH_ICE        (1:nCols), &!REAL(KIND=r8)   , INTENT(OUT ) :: CHV        (nCols,ROWS)
         V_S_ICE       (1:nCols), &!REAL(KIND=r8)   , INTENT(OUT ) :: V_S        (nCols,ROWS)
         RECIP_L_MO_ICE(1:nCols) &!REAL(KIND=r8)   , INTENT(OUT ) :: RECIP_L_MO (nCols,ROWS) 
         )


    ! Marginal Ice Zone
    ! DEPENDS ON: fcdch_sea
    CALL FCDCH_SEA(&
         nCols                     , &!INTEGER, INTENT(IN  ) :: nCols      ! IN Number of X points?
         mskant     (1:nCols), &!INTEGER, INTENT(IN   ) ::  nCols                     ! IN Number of X points?
         FLANDG        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: FLANDG     (nCols,ROWS) ! IN Land fraction
         DB_ICE        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: DB         (nCols,ROWS) !
         VSHR_SSI      (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: VSHR       (nCols,ROWS) 
         Z0_MIZ        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z0M        (nCols,ROWS) 
         Z0_MIZ        (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z0H        (nCols,ROWS) 
         ZH            (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: ZH         (nCols,ROWS)  
         Z1_UV         (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z1_UV      (nCols,ROWS)  
         Z1_TQ_SEA     (1:nCols), &!REAL(KIND=r8)   , INTENT(IN  ) :: Z1_TQ      (nCols,ROWS)  
         CD_MIZ        (1:nCols), &!REAL(KIND=r8)   , INTENT(OUT ) :: CDV        (nCols,ROWS) 
         CH_MIZ        (1:nCols), &!REAL(KIND=r8)   , INTENT(OUT ) :: CHV        (nCols,ROWS)	
         V_S_MIZ       (1:nCols), &!REAL(KIND=r8)   , INTENT(OUT ) :: V_S        (nCols,ROWS)
         RECIP_L_MO_MIZ(1:nCols) &!REAL(KIND=r8)   , INTENT(OUT ) :: RECIP_L_MO (nCols,ROWS) 
         )

    ! Sea and sea-ice leads
    ! DEPENDS ON: fcdch_sea
    CALL FCDCH_SEA(&
         nCols                     , & !INTEGER, INTENT(IN  ) :: nCols! IN Number of X points?
         mskant        (1:nCols), &!INTEGER, INTENT(IN   ) ::  nCols                     ! IN Number of X points?
         FLANDG        (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: FLANDG(nCols,ROWS) ! IN Land fraction
         DB_SEA        (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: DB    (nCols,ROWS) !
         VSHR_SSI      (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: VSHR  (nCols,ROWS) 
         Z0MSEA        (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: Z0M   (nCols,ROWS) 
         Z0H_SEA       (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: Z0H   (nCols,ROWS) 
         ZH            (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: ZH   (nCols,ROWS)  
         Z1_UV         (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: Z1_UV(nCols,ROWS)  
         Z1_TQ_SEA     (1:nCols), & !REAL(KIND=r8)   , INTENT(IN  ) :: Z1_TQ(nCols,ROWS)  
         CD_SEA        (1:nCols), & !REAL(KIND=r8)   , INTENT(OUT ) :: CDV   (nCols,ROWS) 
         CH_SEA        (1:nCols), & !REAL(KIND=r8)   , INTENT(OUT ) :: CHV   (nCols,ROWS)     
         V_S_SEA       (1:nCols), & !REAL(KIND=r8)   , INTENT(OUT ) :: V_S(nCols,ROWS)
         RECIP_L_MO_SEA(1:nCols) & !REAL(KIND=r8)   , INTENT(OUT ) :: RECIP_L_MO(nCols,ROWS) 
         ) 
    !-----------------------------------------------------------------------
    ! Calculate gridbox-means of transfer coefficients.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          CD_LAND(I) = 0.0_r8
          CD_SSI(I) = 0.0_r8
          CH_SSI(I) = 0.0_r8
          CD(I) = 0.0_r8
          CH(I) = 0.0_r8
          IF (L_DUST) CD_STD_DUST(I) = 0.0_r8

          ! Sea and sea-ice
          IF (FLANDG(I) <  1.0_r8 ) THEN
             IF ( ICE_FRACT(I)  <   0.7_r8 ) THEN
                CD_SSI(I) = ( ICE_FRACT(I)*CD_MIZ(I) + (0.7_r8-ICE_FRACT(I))*CD_SEA(I) ) / 0.7_r8  ! P2430.5
                CH_SSI(I) = ( ICE_FRACT(I)*CH_MIZ(I) + (0.7_r8-ICE_FRACT(I))*CH_SEA(I) ) / 0.7_r8  ! P2430.4
             ELSE
                CD_SSI(I) = ( (1.0_r8-ICE_FRACT(I))*CD_MIZ(I) + (ICE_FRACT(I)-0.7_r8)*CD_ICE(I) ) / 0.3_r8  ! P2430.7
                CH_SSI(I) = ( (1.0_r8-ICE_FRACT(I))*CH_MIZ(I) + (ICE_FRACT(I)-0.7_r8)*CH_ICE(I) ) / 0.3_r8  ! P2430.7
             ENDIF
             CD(I)=(1.0_r8-FLANDG(I))*CD_SSI(I)
             CH(I)=(1.0_r8-FLANDG(I))*CH_SSI(I)
             IF (L_DUST) CD_STD_DUST(I)=(1.0_r8-FLANDG(I))*CD_SSI(I)
          ENDIF
       END IF
    ENDDO

    ! Land tiles
    !      DO N=1,NTILES
    !CDIR NODEP
    !        DO K=1,TILE_PTS(N)
    !          L = TILE_INDEX(K,N)
    !          J=(LAND_INDEX(L)-1)/nCols + 1
    !          I = LAND_INDEX(L) - (J-1)*nCols
    !          CD_LAND(I,J) = CD_LAND(I,J) + TILE_FRAC(L,N)*CD_TILE(L,N)
    !          CD(I,J) = CD(I,J) + FLANDG(I,J)*TILE_FRAC(L,N)*CD_TILE(L,N)
    !          CH(I,J) = CH(I,J) + FLANDG(I,J)*TILE_FRAC(L,N)*CH_TILE(L,N)
    !          IF (L_DUST) CD_STD_DUST(I,J) = CD_STD_DUST(I,J) +             &
    !     &                  FLANDG(I,J)*TILE_FRAC(L,N)*CD_STD(L,N)!
    !
    !        ENDDO
    !      ENDDO


    !-----------------------------------------------------------------------
    !!  4.3 Calculate the surface exchange coefficients RHOK(*) and
    !       resistances for use in Sulphur Cycle
    !       (Note that CD_STD, CH and VSHR should never = 0)
    !     RHOSTAR * CD * VSHR stored for diagnostic output before
    !     horizontal interpolation.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          RHO_ARESIST(I) = 0.0_r8
          !        RHO_ARESIST_LAND(I)=0.0
          ARESIST(I) = 0.0_r8
          RESIST_B(I) = 0.0_r8
          !        RHOKM_LAND(I) = 0.
          RHOKM_SSI(I) = 0.0_r8

          ! Sea and sea-ice
          IF ( FLANDG(I) <  1.0_r8 ) THEN
             RHOKM_SSI(I)    = RHOSTAR(I)*CD_SSI(I)*VSHR_SSI(I)
             !                                                          ! P243.124
             RHOKH_1_SICE(I) = RHOSTAR(I)*CH_SSI(I)*VSHR_SSI(I)
             !                                                           ! P243.125
             RHO_ARESIST(I) = RHOSTAR(I)*CD_SSI(I)*VSHR_SSI(I)
             ARESIST (I) =  1.0_r8 / (CD_SSI(I) * VSHR_SSI(I))
             RESIST_B(I)= (CD_SSI(I)/CH_SSI(I) - 1.0_r8) * ARESIST(I)
          ENDIF
       END IF
    ENDDO


    !-----------------------------------------------------------------------
    !!  Calculate surface layer resistance for mineral dust
    !-----------------------------------------------------------------------
    IF (L_DUST) THEN

       DO I = 1,nCols
          IF(mskant(i) == 1_i8)THEN

             DO IDIV=1,NDIV
                R_B_DUST(I,IDIV)=0.0_r8
             ENDDO !IDIV
             VSHR(I)=(1.0_r8 - FLANDG(I)) * VSHR_SSI(I) !+ FLANDG(I) * VSHR_LAND(I)
          END IF
       ENDDO !I

       ! DEPENDS ON: dustresb
       CALL DUSTRESB (&
            nCols       , &!INTEGER,INTENT(IN ) :: nCols	   !IN
            mskant, &
            PSTAR       , &!REAL(KIND=r8), INTENT(IN   ) :: PSTAR(nCols,ROWS)  !IN surface pressure
            TSTAR_SEA   , &!REAL(KIND=r8), INTENT(IN   ) :: TSTAR(nCols,ROWS)  !IN surface temperature
            RHOSTAR     , &!REAL(KIND=r8), INTENT(IN   ) :: RHOSTAR(nCols,ROWS)!IN surface air density
            VSHR        , &!REAL(KIND=r8), INTENT(IN   ) :: VSHR(nCols,ROWS)   !IN surface to lowest lev windspeed difference
            CD_STD_DUST , &!REAL(KIND=r8), INTENT(IN   ) :: CD_STD_DUST(nCols,ROWS) !IN surface transfer coeffient for
            R_B_DUST      &! REAL(KIND=r8), INTENT(OUT  ) :: R_B_DUST(nCols,ROWS,NDIV) !OUT surface layer resistance for
            )

    ENDIF !(L_DUST)

    !-----------------------------------------------------------------------
    !!  Calculate local and gridbox-average surface fluxes of heat and
    !!  moisture.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          FTL_1(I) = 0.0_r8
          FQW_1(I) = 0.0_r8
       END IF
    ENDDO

    ! Sea and sea-ice
    ! DEPENDS ON: sf_flux_sea
    CALL SF_FLUX_SEA (  &
         nCols                          , &!INTEGER, INTENT(IN   ) :: nCols                   ! IN Number of X points?
         NSICE                          , &!INTEGER, INTENT(IN   ) :: NSICE                   ! IN Number of sea-ice points.
         mskant                         , &
         SICE_INDEX       (1:nCols,1:2) , &!INTEGER, INTENT(IN   ) :: SICE_INDEX(nCols*ROWS,2)! IN Index of sea-ice points
         FLANDG           (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: FLANDG    (nCols,ROWS)
         ICE_FRACT        (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: ICE_FRACT (nCols,ROWS)  ! IN Fraction of gridbox which is sea-ice.
         QS1              (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: QS1       (nCols,ROWS)  ! IN Sat. specific humidity qsat(TL_1,PSTAR)
         QSTAR_ICE        (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: QSTAR_ICE (nCols,ROWS)  ! IN Surface qsat for sea-ice.
         QSTAR_SEA        (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: QSTAR_SEA (nCols,ROWS)  ! IN Surface qsat for sea or sea-ice leads.
         QW_1             (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: QW_1      (nCols,ROWS)  ! IN Total water content of lowest
         RADNET           (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: RADNET    (nCols,ROWS)  ! IN Net surface radiation (W/m2)positive downwards
         RHOKH_1_SICE     (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: RHOKH_1   (nCols,ROWS)  ! IN Surface exchange coefficient.
         TI               (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: TI        (nCols,ROWS)  ! IN Temperature of sea-ice surfacelayer (K)
         TL_1             (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: TL_1      (nCols,ROWS)  ! IN Liquid/frozen water temperature for lowest atmospheric layer (K)
         TSTAR_SICE       (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: TSTAR_SICE(nCols,ROWS)  ! IN Sea-ice surface temperature (K)
         TSTAR_SEA        (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: TSTAR_SEA (nCols,ROWS)  ! IN Sea surface temperature (K).
         Z0_ICE           (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z0H_ICE   (nCols,ROWS)  ! IN Sea-ice heat and moisture
         Z0_ICE           (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z0M_ICE   (nCols,ROWS)  ! IN Sea-ice momentum roughness
         Z0H_SEA          (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z0H_SEA   (nCols,ROWS)  ! IN Sea and lead heat and moisture
         Z0MSEA           (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z0M_SEA   (nCols,ROWS)  ! IN Sea and lead momentum roughness
         Z1_TQ_SEA        (1:nCols)  , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z1_TQ     (nCols,ROWS)  ! IN Height of lowest atmospheric
         SeaSalinityFactor           , &!REAL(KIND=r8)   , Intent(IN   ) :: SeaSalinityFactor
         ALPHA1_SICE      (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: ALPHA1    (nCols,ROWS)
         ASHTF            (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: ASHTF     (nCols,ROWS)
         E_SEA            (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: E_SEA     (nCols,ROWS)
         FQW_ICE          (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: FQW_ICE   (nCols,ROWS)
         FQW_1            (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: FQW_1     (nCols,ROWS) 
         FTL_ICE          (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: FTL_ICE   (nCols,ROWS)
         FTL_1            (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: FTL_1     (nCols,ROWS)
         H_SEA            (1:nCols)  , &!REAL(KIND=r8)   , INTENT(OUT  ) :: H_SEA     (nCols,ROWS)
         RHOKPM_SICE      (1:nCols)    &!REAL(KIND=r8)   , INTENT(OUT  ) :: RHOKPM    (nCols,ROWS)     ! OUT Modified surface exchange
         )
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN
          rmi   (i)=CD_SSI(I)*VSHR_SSI(I)
          rhi   (i)=CH_SSI(I)*VSHR_SSI(I)
          evap  (i)=FQW_1(i)*hltm
          sens  (i)=FTL_1(i)*CP
       END IF
    ENDDO

    !-----------------------------------------------------------------------
    !!  4.4   Calculate the standard deviations of layer 1 turbulent
    !!        fluctuations of temperature and humidity using approximate
    !!        formulae from first order closure.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          Q1_SD(I) = 0.0_r8
          T1_SD(I) = 0.0_r8
       END IF
    ENDDO

    ! Sea and sea-ice
    ! DEPENDS ON: stdev1_sea
    CALL STDEV1_SEA (&
         nCols                   , & !INTEGER, INTENT(IN	) :: nCols		    ! IN Number of X points?
         mskant , &
         FLANDG      (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: FLANDG   (nCols,ROWS)  ! IN Land fraction.
         BQ_1        (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: BQ_1     (nCols,ROWS)  ! IN Buoyancy parameter.
         BT_1        (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: BT_1     (nCols,ROWS)  ! IN Buoyancy parameter.
         FQW_1       (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: FQW_1    (nCols,ROWS)  ! IN Surface flux of QW.
         FTL_1       (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: FTL_1    (nCols,ROWS)  ! IN Surface flux of TL.
         ICE_FRACT   (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: ICE_FRACT(nCols,ROWS)  ! IN Fraction of gridbox which is sea-ice.
         RHOKM_SSI   (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: RHOKM_1  (nCols,ROWS)  ! IN Surface momentum exchange
         RHOSTAR     (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: RHOSTAR  (nCols,ROWS)  ! IN Surface air density.
         VSHR_SSI    (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: VSHR     (nCols,ROWS)  ! IN Magnitude of surface-
         Z0MSEA      (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: Z0MSEA   (nCols,ROWS)  ! IN Sea roughness length.
         Z0_ICE      (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: Z0_ICE   (nCols,ROWS)  ! IN Sea-ice roughness length.
         Z1_TQ_SEA   (1:nCols), & !REAL(KIND=r8)   , INTENT(IN	) :: Z1_TQ    (nCols,ROWS)  ! IN Height of lowest tq level.
         Q1_SD       (1:nCols), & !REAL(KIND=r8)   , INTENT(out  ) :: Q1_SD    (nCols,ROWS)  ! OUT Standard deviation of
         T1_SD       (1:nCols)  ) !REAL(KIND=r8)   , INTENT(out  ) :: T1_SD    (nCols,ROWS)  ! OUT Standard deviation of


    !-----------------------------------------------------------------------
    !!  4.6 For sea points, calculate the wind mixing energy flux and the
    !!      sea-surface roughness length on the P-grid, using time-level n
    !!      quantities.
    !-----------------------------------------------------------------------

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF (SFME) FME(I) = 0.0_r8
          IF (FLANDG(I) <  1.0_r8) THEN
             TAU = RHOKM_SSI(I) * VSHR_SSI(I)             ! P243.130
             IF (ICE_FRACT(I)  >   0.0_r8)                                  &
                  &      TAU = RHOSTAR(I) * CD_SEA(I) * VSHR_SSI(I) * VSHR_SSI(I)

             IF (SFME)FME(I) = (1.0_r8-ICE_FRACT(I)) * TAU * SQRT(TAU/RHOSEA)
             !                                                            ! P243.96
             ! Limit Z0MSEA to 0.154m for TAU very small
             Z0MSEA(I) = 1.54E-6_r8 / (SQRT(TAU/RHOSTAR(I)) + 1.0E-5_r8)     &
                  &               +  (CHARNOCK/G) * (TAU / RHOSTAR(I))
             Z0MSEA(I) = MAX ( Z0HSEA , Z0MSEA(I) )
             !                                       ... P243.B6 (Charnock formula)
             !                    TAU/RHOSTAR is "mod VS squared", see eqn P243.131
          ENDIF
       END IF
    ENDDO

  END SUBROUTINE SF_EXCH
  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************

  !-----------------------------------------------------------------------
  !
  ! Subroutines SF_FLUX_LAND and SF_FLUX_SEA to calculate explicit surface
  ! fluxes of heat and moisture
  !
  !
  !    Model            Modification history
  !   version  date
  !    5.2   15/11/00   New Deck         M. Best
  !    5.3  25/04/01  Add coastal tiling. Nic Gedney
  !  6.1  01/09/04  Calculate potential evaporation related variables.
  !                                                          Nic Gedney
  !    6.2  07/11/05  Allow for the salinity of sea water in
  !                     the evaporative flux.
  !                     J. M. Edwards
  !
  !    Programming standard:
  !
  !-----------------------------------------------------------------------

  !     SUBROUTINE SF_FLUX_LAND-------------------------------------------
  !
  !     Calculate explicit surface fluxes of heat and moisture over
  !     land tiles
  !
  !     ------------------------------------------------------------------

  !     SUBROUTINE SF_FLUX_SEA--------------------------------------------
  !
  !     Calculate explicit surface fluxes of heat and moisture over sea
  !     and sea-ice
  !
  !     ------------------------------------------------------------------
  SUBROUTINE SF_FLUX_SEA (&
       nCols       , &!INTEGER, INTENT(IN   ) :: nCols ! IN Number of X points?
       NSICE            , &!INTEGER, INTENT(IN   ) :: NSICE      ! IN Number of sea-ice points.
       mskant, &
       SICE_INDEX       , &!INTEGER, INTENT(IN   ) :: SICE_INDEX(nCols*ROWS,2)! IN Index of sea-ice points
       FLANDG           , &!REAL(KIND=r8), INTENT(IN   ) :: FLANDG(nCols,ROWS)
       ICE_FRACT        , &!REAL(KIND=r8), INTENT(IN   ) :: ICE_FRACT(nCols,ROWS) ! IN Fraction of gridbox which is sea-ice.
       QS1              , &!REAL(KIND=r8), INTENT(IN   ) ::QS1(nCols,ROWS)! IN Sat. specific humidity qsat(TL_1,PSTAR)
       QSTAR_ICE        , &!REAL(KIND=r8), INTENT(IN   ) ::QSTAR_ICE(nCols,ROWS)! IN Surface qsat for sea-ice.
       QSTAR_SEA        , &!REAL(KIND=r8), INTENT(IN   ) ::QSTAR_SEA(nCols,ROWS)! IN Surface qsat for sea or sea-ice leads.
       QW_1             , &!REAL(KIND=r8), INTENT(IN   ) ::QW_1(nCols,ROWS)! IN Total water content of lowest
       RADNET           , &!REAL(KIND=r8), INTENT(IN   ) ::RADNET(nCols,ROWS)! IN Net surface radiation (W/m2)positive downwards
       RHOKH_1          , &!REAL(KIND=r8), INTENT(IN   ) ::RHOKH_1(nCols,ROWS)! IN Surface exchange coefficient.
       TI               , &!REAL(KIND=r8), INTENT(IN   ) ::TI(nCols,ROWS)  ! IN Temperature of sea-ice surfacelayer (K)
       TL_1             , &!REAL(KIND=r8), INTENT(IN   ) ::TL_1(nCols,ROWS) ! IN Liquid/frozen water temperature for lowest atmospheric layer (K)
       TSTAR_SICE       , &!REAL(KIND=r8), INTENT(IN   ) ::TSTAR_SICE(nCols,ROWS) ! IN Sea-ice surface temperature (K)
       TSTAR_SEA        , &!REAL(KIND=r8), INTENT(IN   ) ::TSTAR_SEA(nCols,ROWS) ! IN Sea surface temperature (K).
       Z0H_ICE          , &!REAL(KIND=r8), INTENT(IN   ) ::Z0H_ICE(nCols,ROWS)! IN Sea-ice heat and moisture
       Z0M_ICE          , &!REAL(KIND=r8), INTENT(IN   ) ::Z0M_ICE(nCols,ROWS)! IN Sea-ice momentum roughness
       Z0H_SEA          , &!REAL(KIND=r8), INTENT(IN   ) ::Z0H_SEA(nCols,ROWS) ! IN Sea and lead heat and moisture
       Z0M_SEA          , &!REAL(KIND=r8), INTENT(IN   ) ::Z0M_SEA(nCols,ROWS) ! IN Sea and lead momentum roughness
       Z1_TQ            , &!REAL(KIND=r8), INTENT(IN   ) ::Z1_TQ(nCols,ROWS)! IN Height of lowest atmospheric
       SeaSalinityFactor, &!REAL(KIND=r8), Intent(IN) :: SeaSalinityFactor
       ALPHA1           , &!REAL(KIND=r8) , INTENT(OUT  ):: ALPHA1(nCols,ROWS)
       ASHTF            , &!REAL(KIND=r8) , INTENT(OUT  ):: ASHTF(nCols,ROWS)                     
       E_SEA            , &!REAL(KIND=r8) , INTENT(OUT  ):: E_SEA(nCols,ROWS)
       FQW_ICE          , &!REAL(KIND=r8) , INTENT(OUT  ):: FQW_ICE(nCols,ROWS)
       FQW_1            , &!REAL(KIND=r8) , INTENT(OUT  ):: FQW_1(nCols,ROWS) 
       FTL_ICE          , &!REAL(KIND=r8) , INTENT(OUT  ):: FTL_ICE(nCols,ROWS)               
       FTL_1            , &!REAL(KIND=r8) , INTENT(OUT  ):: FTL_1(nCols,ROWS)              
       H_SEA            , &!REAL(KIND=r8) , INTENT(OUT  ):: H_SEA(nCols,ROWS)                
       RHOKPM            &!REAL(KIND=r8) , INTENT(OUT  ):: RHOKPM(nCols,ROWS)     ! OUT Modified surface exchange
       )
    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: nCols ! IN Number of X points?
    INTEGER, INTENT(IN   ) :: NSICE      ! IN Number of sea-ice points.
    INTEGER, INTENT(IN   ) :: SICE_INDEX(nCols,2)! IN Index of sea-ice points
    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: FLANDG(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: ICE_FRACT(nCols) ! IN Fraction of gridbox which is sea-ice.
    REAL(KIND=r8), INTENT(IN   ) ::QS1(nCols)! IN Sat. specific humidity qsat(TL_1,PSTAR)
    REAL(KIND=r8), INTENT(IN   ) ::QSTAR_ICE(nCols)! IN Surface qsat for sea-ice.
    REAL(KIND=r8), INTENT(IN   ) ::QSTAR_SEA(nCols)! IN Surface qsat for sea or sea-ice leads.
    REAL(KIND=r8), INTENT(IN   ) ::QW_1(nCols)! IN Total water content of lowest
    !                                  !    atmospheric layer
    !                                  !    (kg per kg air).
    REAL(KIND=r8), INTENT(IN   ) ::RADNET(nCols)! IN Net surface radiation (W/m2)
    !                                  !    positive downwards
    REAL(KIND=r8), INTENT(IN   ) ::RHOKH_1(nCols)! IN Surface exchange coefficient.
    REAL(KIND=r8), INTENT(IN   ) ::TI(nCols)  ! IN Temperature of sea-ice surface
    !                                  !    layer (K)
    REAL(KIND=r8), INTENT(IN   ) ::TL_1(nCols) ! IN Liquid/frozen water temperature
    !                                  !    for lowest atmospheric layer (K)
    REAL(KIND=r8), INTENT(IN   ) ::TSTAR_SICE(nCols) ! IN Sea-ice surface temperature (K)
    REAL(KIND=r8), INTENT(IN   ) ::TSTAR_SEA(nCols) ! IN Sea surface temperature (K).
    REAL(KIND=r8), INTENT(IN   ) ::Z0H_ICE(nCols)! IN Sea-ice heat and moisture
    !                                  !    roughness length (m).
    REAL(KIND=r8), INTENT(IN   ) ::Z0M_ICE(nCols)! IN Sea-ice momentum roughness
    !                                  !    length (m).
    REAL(KIND=r8), INTENT(IN   ) ::Z0H_SEA(nCols) ! IN Sea and lead heat and moisture
    !                                  !    roughness length (m).
    REAL(KIND=r8), INTENT(IN   ) ::Z0M_SEA(nCols) ! IN Sea and lead momentum roughness
    !                                  !    length.
    REAL(KIND=r8), INTENT(IN   ) ::Z1_TQ(nCols)! IN Height of lowest atmospheric
    !                                  !    level (m).

    !
    !
    REAL(KIND=r8), INTENT(IN) :: SeaSalinityFactor
    !                                  ! Factor allowing for the
    !                                  ! effect of the salinity of
    !                                  ! sea water on the evaporative
    !                                  ! flux.
    REAL(KIND=r8) , INTENT(OUT  ):: ALPHA1(nCols)
    ! OUT Gradient of saturated specific
    !                                  !     humidity with respect to
    !                                  !     temperature between the bottom
    !                                  !     model layer and the surface.
    REAL(KIND=r8) , INTENT(OUT  ):: ASHTF(nCols)
    ! OUT Coefficient to calculate
    !                                  !     surface heat flux into sea-ice
    !                                  !     (W/m2/K).
    REAL(KIND=r8) , INTENT(OUT  ):: E_SEA(nCols) ! OUT Evaporation from sea times
    !                                  !     leads fraction (kg/m2/s).
    REAL(KIND=r8) , INTENT(OUT  ):: FQW_ICE(nCols)
    ! OUT Surface flux of QW for sea-ice.
    REAL(KIND=r8) , INTENT(OUT  ):: FQW_1(nCols)
    ! OUT GBM surface flux of QW
    !                                  !     (kg/m2/s).
    REAL(KIND=r8) , INTENT(OUT  ):: FTL_ICE(nCols)
    ! OUT Surface flux of TL for sea-ice.
    REAL(KIND=r8) , INTENT(OUT  ):: FTL_1(nCols)
    ! OUT GBM surface flux of TL.
    REAL(KIND=r8) , INTENT(OUT  ):: H_SEA(nCols)! OUT Surface sensible heat flux over
    !                                  !     sea times leads fraction (W/m2)
    REAL(KIND=r8) , INTENT(OUT  ):: RHOKPM(nCols)     ! OUT Modified surface exchange
    !                                  !     coefficient.

    !*L------------------COMDECK C_EPSLON-----------------------------------
    ! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

    REAL(KIND=r8), PARAMETER :: Epsilon   = 0.62198_r8
    REAL(KIND=r8), PARAMETER :: C_Virtual = 1.0_r8/Epsilon-1.0_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_O_DG_C-----------------------------------
    ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
    ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
    ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

    REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
    REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
    REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------
    ! C_KAPPAI start

    ! Thermal conductivity of sea-ice (W per m per K).
    REAL(KIND=r8),PARAMETER:: KAPPAI=2.09_r8

    ! Effective thickness of sea-ice surface layer (m).
    REAL(KIND=r8),PARAMETER:: DE = 0.1_r8

    ! C_KAPPAI end
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8

    ! C_LHEAT end
    !*L------------------COMDECK C_R_CP-------------------------------------
    ! History:
    ! Version  Date      Comment.
    !  5.0  07/05/99  Add variable P_zero for consistency with
    !                 conversion to C-P 'C' dynamics grid. R. Rawlins
    !  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

    ! R IS GAS CONSTANT FOR DRY AIR
    ! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
    ! PREF IS REFERENCE SURFACE PRESSURE

    REAL(KIND=r8), PARAMETER  :: R      = 287.05_r8
    REAL(KIND=r8), PARAMETER  :: CP     = 1005.0_r8
    REAL(KIND=r8), PARAMETER  :: Kappa  = R/CP
    REAL(KIND=r8), PARAMETER  :: Pref   = 100000.0_r8

    ! Reference surface pressure = PREF
    REAL(KIND=r8), PARAMETER  :: P_zero = Pref
    REAL(KIND=r8), PARAMETER  :: sclht  = 6.8E+03_r8  ! Scale Height H
    !*----------------------------------------------------------------------
    ! CSIGMA start
    ! Stefan-Boltzmann constant (W/m**2/K**4).
    REAL(KIND=r8), PARAMETER ::  SBCON=5.67E-8_r8
    ! CSIGMA end

    ! Derived local parameters.
    REAL(KIND=r8) GRCP,LS
    PARAMETER (                                                       &
         & GRCP=G/CP                                                        &
         &,LS=LF+LC                                                         &
                                ! Latent heat of sublimation.
         & )

    ! Scalars
    INTEGER                                                           &
         & I                                                              &
                                ! Horizontal field index.
         &,K                   ! Sea-ice field index.
    REAL(KIND=r8)                                                              &
         & DQ1                                                              &
                                ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1
         &,D_T                                                              &
                                ! Temporary in calculation of alpha1.
         &,RAD_REDUC           ! Radiation term required for surface flux
    !                          ! calcs.

    !      EXTERNAL TIMER

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('SF_FLUX ',3)
    !      ENDIF

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          ALPHA1(I) = 0.0_r8
          E_SEA(I) = 0.0_r8
          H_SEA(I) = 0.0_r8
          FQW_ICE(I) = 0.0_r8
          FTL_ICE(I) = 0.0_r8
          RHOKPM(I) = 0.0_r8
       END IF
    ENDDO

    !----------------------------------------------------------------------
    !!  1 Calculate gradient of saturated specific humidity for use in
    !!    calculation of surface fluxes - only required for sea-ice points
    !----------------------------------------------------------------------
    DO K=1,NSICE
       I = SICE_INDEX(K,1)
       IF(mskant(i) == 1_i8)THEN

          D_T = TSTAR_SICE(I) - TL_1(I)
          IF (D_T  >   0.05_r8 .OR. D_T  <   -0.05_r8) THEN
             ALPHA1(I) = (QSTAR_ICE(I) - QS1(I)) / D_T
          ELSEIF (TL_1(I)  >   TM) THEN
             ALPHA1(I) = EPSILON*LC*QS1(I)*( 1.0_r8+C_VIRTUAL*QS1(I) )  &
                  &                                        /(R*TL_1(I)*TL_1(I))
          ELSE
             ALPHA1(I) = EPSILON*LS*QS1(I)*( 1.0_r8+C_VIRTUAL*QS1(I) )  &
                  &                                        /(R*TL_1(I)*TL_1(I))
          ENDIF
          ASHTF(I) = 2 * KAPPAI / DE + 4*SBCON*TI(I)**3
       END IF
    ENDDO

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF (FLANDG(I) <  1.0_r8 ) THEN

             E_SEA(I) = - (1.0_r8 - ICE_FRACT(I))    * RHOKH_1(I)*(QW_1(I) - SeaSalinityFactor*QSTAR_SEA(I))
             H_SEA(I) = - (1.0_r8 - ICE_FRACT(I)) *CP* RHOKH_1(I) * ( TL_1(I) - TSTAR_SEA(I)   &
                  &              + GRCP*(Z1_TQ(I) + Z0M_SEA(I) - Z0H_SEA(I)) )

             IF ( ICE_FRACT(I)  >   0.0_r8 ) THEN
                ! Sea-ice
                RHOKPM(I) = RHOKH_1(I) / ( ASHTF(I) + RHOKH_1(I)*(LS*ALPHA1(I) + CP) )
                RAD_REDUC = RADNET(I) - ICE_FRACT(I) * ASHTF(I) *     &
                     &                  ( TL_1(I) - TI(I) +                         &
                     &                GRCP*(Z1_TQ(I) + Z0M_ICE(I) - Z0H_ICE(I)) )
                DQ1 = QS1(I) - QW_1(I) +                                &
                     &                           GRCP*ALPHA1(I)*                      &
                     &                       (Z1_TQ(I) + Z0M_ICE(I) - Z0H_ICE(I))
                FQW_ICE(I) = RHOKPM(I) * ( ALPHA1(I)*RAD_REDUC +      &
                     &                        (CP*RHOKH_1(I) +                        &
                     &                         ASHTF(I))*DQ1*ICE_FRACT(I) )
                FTL_ICE(I) = RHOKPM(I) * ( RAD_REDUC -                  &
                     &                             ICE_FRACT(I)*LS*RHOKH_1(I)*DQ1 )

             ENDIF

             FTL_1(I) = (1.0_r8-FLANDG(I))*(FTL_ICE(I) + H_SEA(I) / CP)
             FQW_1(I) = (1.0_r8-FLANDG(I))*(FQW_ICE(I) + E_SEA(I))

          ENDIF
       END IF
    ENDDO

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('SF_FLUX ',4)
    !      ENDIF

    RETURN
  END SUBROUTINE SF_FLUX_SEA



  !
  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  ! subroutine DUSTRESB
  !
  !
  ! Purpose:
  !   To calculate the surface layer resistance for mineral dust
  !
  ! Called by sfexch
  !
  ! Current owners of code: S.Woodward
  !
  ! History:
  ! Version     Date     Comment
  ! -------     ----     -------
  !
  !   5.5      12/02/03  Original code   S Woodward
  !   6.2      11/05/06  Fix for bit-non-reproducibility   A. Malcolm
  !   6.4      15/01/07  Malcolm McVean's portability fix   S.Woodward
  !
  ! Code Description:
  !  Language: FORTRAN77 + common extensions
  !  This code is written to UMDP3 v6 programming standards
  !
  ! Documentation: "Modelling the atmospheric lifecycle..."
  !                 Woodward, JGR106, D16, pp18155-18166
  !---------------------------------------------------------------------
  !
  SUBROUTINE DUSTRESB(                                             &
       nCols      , &!INTEGER,INTENT(IN ) :: nCols        !IN
       mskant     , &
       
       PSTAR      , &!REAL(KIND=r8), INTENT(IN   ) :: PSTAR(nCols,ROWS)  !IN surface pressure
       TSTAR      , &!REAL(KIND=r8), INTENT(IN   ) :: TSTAR(nCols,ROWS)  !IN surface temperature
       RHOSTAR    , &!REAL(KIND=r8), INTENT(IN   ) :: RHOSTAR(nCols,ROWS)!IN surface air density
       VSHR       , &!REAL(KIND=r8), INTENT(IN   ) :: VSHR(nCols,ROWS)   !IN surface to lowest lev windspeed difference
       CD_STD_DUST, &!REAL(KIND=r8), INTENT(IN   ) :: CD_STD_DUST(nCols,ROWS) !IN surface transfer coeffient for
       R_B_DUST     &! REAL(KIND=r8), INTENT(OUT  ) :: R_B_DUST(nCols,ROWS,NDIV) !OUT surface layer resistance for
       )
    IMPLICIT NONE

    !C_DUST_NDIV.............................................................
    ! Description: Contains parameters for mineral dust code
    ! Current Code Owner: Stephanie Woodward
    !
    ! History:
    ! Version  Date     Comment
    ! -------  ----     -------
    !  5.5      12/02/03  Original Code.   Stephanie Woodward
    !
    ! Declarations:
    !
    INTEGER NDIV        ! number of particle size divisions
    PARAMETER (NDIV = 6)
    !.....................................................................

    INTEGER,INTENT(IN ) :: nCols        !IN
    INTEGER(KIND=i8),INTENT(IN ) :: mskant(nCols)

    REAL(KIND=r8), INTENT(IN   ) :: PSTAR(nCols)  !IN surface pressure
    REAL(KIND=r8), INTENT(IN   ) :: TSTAR(nCols)  !IN surface temperature
    REAL(KIND=r8), INTENT(IN   ) :: RHOSTAR(nCols)!IN surface air density
    REAL(KIND=r8), INTENT(IN   ) :: VSHR(nCols)   !IN surface to lowest lev windspeed difference
    REAL(KIND=r8), INTENT(IN   ) :: CD_STD_DUST(nCols) !IN surface transfer coeffient for
    !                                                    !   momentum, excluding orographic
    !                                                    !   form drag

    REAL(KIND=r8), INTENT(OUT  ) :: R_B_DUST(nCols,NDIV) !OUT surface layer resistance for
    !                                                      !    mineral dust

    !     local variables

    INTEGER                                                           &
         & IDIV                                                             &
                                !loop counter, dust divisions
         &,I                                                                &
                                !loop counter
         &,LEV1 !number of levels for vstokes calculation

    REAL(KIND=r8)  :: ETAA(nCols)                                            &
                                !dynamic viscosity of air
         &,VSTOKES1(nCols)                                        &
                                !gravitational settling velocity, lev1
         &,NSTOKES(nCols)                                         &
                                !stokes number = VstokesVshrVshr/nu g
         &,NSCHMIDT(nCols)                                        &
                                !schmidt number = nu/diffusivit
         &,CCF(nCols)                                             &
                                !Cunningham correction factor
         &,STOKES_EXP                                                       &
                                !stokes term in R_B_DUST equation
         &,SMALLP                    !small +ve number, negligible compared to 1

    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_PI---------------------------------------
    !LL
    !LL 4.0 19/09/95  New value for PI. Old value incorrect
    !LL               from 12th decimal place. D. Robinson
    !LL 5.1 7/03/00   Fixed/Free format P.Selwood
    !LL

    ! Pi
    REAL(KIND=r8), PARAMETER :: Pi                 = 3.14159265358979323846_r8

    ! Conversion factor degrees to radians
    REAL(KIND=r8), PARAMETER :: Pi_Over_180        = Pi/180.0_r8

    ! Conversion factor radians to degrees
    REAL(KIND=r8), PARAMETER :: Recip_Pi_Over_180  = 180.0_r8/Pi

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_O_DG_C-----------------------------------
    ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
    ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
    ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

    REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
    REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
    REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_R_CP-------------------------------------
    ! History:
    ! Version  Date      Comment.
    !  5.0  07/05/99  Add variable P_zero for consistency with
    !                 conversion to C-P 'C' dynamics grid. R. Rawlins
    !  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

    ! R IS GAS CONSTANT FOR DRY AIR
    ! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
    ! PREF IS REFERENCE SURFACE PRESSURE

    REAL(KIND=r8), PARAMETER  :: R      = 287.05_r8
    REAL(KIND=r8), PARAMETER  :: CP     = 1005.0_r8
    REAL(KIND=r8), PARAMETER  :: Kappa  = R/CP
    REAL(KIND=r8), PARAMETER  :: Pref   = 100000.0_r8

    ! Reference surface pressure = PREF
    REAL(KIND=r8), PARAMETER  :: P_zero = Pref
    REAL(KIND=r8), PARAMETER  :: sclht  = 6.8E+03_r8  ! Scale Height H
    !*----------------------------------------------------------------------
    !C_DUSTGEN..............................................................
    ! Description: Contains parameters for mineral dust generation
    ! Current Code Owner: Stephanie Woodward
    !
    ! History:
    ! Version  Date     Comment
    ! -------  ----     -------
    !  5.5      12/02/03  Original Code.   Stephanie Woodward
    !  6.2      14/02/06  Alternative emissions terms for HadGEM1A.
    !                                                     Stephanie Woodward
    !  6.4      03/01/07  Alternative U*t terms.   Stephanie Woodward
    ! Parameters for mineral dust generation
    !
    REAL(KIND=r8), PARAMETER :: HORIZ_C = 2.61_r8 ! C in horizontal flux calc.     
    REAL(KIND=r8), PARAMETER :: VERT_A = 13.4_r8 ! A in vertical flux calc
    REAL(KIND=r8), PARAMETER :: VERT_B = -6.0_r8 ! B in vertical flux calc

    REAL(KIND=r8),PARAMETER :: RHOP = 2.65E+3_r8  ! density of a dust particle
    REAL(KIND=r8),PARAMETER :: Z0B = 0.0003_r8    !roughness length for bare soil
    REAL(KIND=r8) Z0S(NDIV)  ! smooth roughness len (calc.d from part size)
    REAL(KIND=r8) DREP(NDIV) ! representative particle diameter
    REAL(KIND=r8) DMAX(NDIV) ! max diameter of particles in each div.
    REAL(KIND=r8) DMIN(NDIV) ! min diameter of particles in each div.
    ! note that by using two arrays here we can set
    ! up overlapping divisions, however this means
    ! we have to be careful to make them consistent
    !
    DATA Z0S/ .374894E-08_r8, .118552E-07_r8, .374894E-07_r8, .118552E-06_r8,     &
         &           .374894E-06_r8, .118552E-05_r8/
    DATA DREP/ .112468E-06_r8, .355656E-06_r8, .112468E-05_r8, .355656E-05_r8,    &
         &           .112468E-04_r8, .355656E-04_r8/
    DATA DMAX/2.0E-7_r8,6.32456E-7_r8,2.0E-6_r8,                               &
         &          6.32456E-6_r8,2.0E-5_r8,6.32456E-5_r8/
    DATA DMIN/6.32456E-8_r8,2.0E-7_r8,6.32456E-7_r8,                           &
         &          2.0E-6_r8,6.32456E-6_r8,2.0E-5_r8/
    !.......................................................................
    !C_DUSTGRAV.............................................................
    ! Description:
    ! Contains parameters for mineral dust gravitational settling
    ! Current Code Owner: Stephanie Woodward
    !
    ! History:
    ! Version  Date     Comment
    ! -------  ----     -------
    !  1      12/02/03  Original Code.   Stephanie Woodward
    !
    REAL(KIND=r8),PARAMETER :: ACCF=1.257_r8 ! Cunningham correction factor term A
    REAL(KIND=r8),PARAMETER :: BCCF=0.4_r8 ! Cunningham correction factor term B
    REAL(KIND=r8),PARAMETER :: CCCF=-1.1_r8 ! Cunningham correction factor term C
    !.......................................................................
    !-------------------COMDECK C_SULCHM--------------------------------
    ! Parameters for Sulphur Cycle Chemistry
    REAL(KIND=r8)                                                              &
         &     EVAPTAU,                                                     &
                                ! timescale for dissolved SO4 to evaporate
         &     NUCTAU,                                                      &
                                ! timescale for accumulation mode particles
         !                           to nucleate once they enter a cloud.
         &     DIFFUSE_AIT,                                                 &
                                ! diffusion coefficient of Aitken particles
         &     K_SO2OH_HI,                                                  &
                                ! high pressure reaction rate limit
         &     K_DMS_OH,                                                    &
                                ! reaction rate for DMS+OH  cc/mcl/s
         &      K4_CH3SO2_O3,                                               &
                                ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
         &      K5_CH3SO3_HO2,                                              &
                                ! Rate coeff for CH3SO3+HO2 -> MSA+O2
         &      RMM_O3,                                                     &
                                ! relative molecular mass O3
         &     BRAT_SO2,                                                    &
                                ! branching ratio for SO2 in DMS oxidn
         &     BRAT_MSA,                                                    &
                                ! branching ratio for MSA in DMS oxidn
         &     AVOGADRO,                                                    &
                                ! no. of molecules in 1 mole
         &     RMM_S,                                                       &
                                ! relative molecular mass S kg/mole
         &     RMM_H2O2,                                                    &
                                ! relative molecular mass H2O2 kg/mole
         &     RMM_HO2,                                                     &
                                ! relative molecular mass HO2 kg/mole
         &     RMM_AIR,                                                     &
                                ! relative molecular mass dry air
         &     RMM_W,                                                       &
                                ! relative molecular mass water
         &     RELM_S_H2O2,                                                 &
                                ! rel atomic mass sulphur/RMM_H2O2
         &     RELM_S_2N,                                                   &
                                ! rel atomic mass Sulphur/2*Nitrogen
         &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
         &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
         &     T1,                                                          &
                                !
         &     FC,                                                          &
                                ! parameters for interpolation between
         &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
         &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
         &     T2,T3,T4,                                                    &
                                !
         &     CLOUDTAU,                                                    &
                                ! air parcel lifetime in cloud
         &     CHEMTAU,                                                     &
                                ! chem lifetime in cloud before oxidn
         &     O3_MIN,                                                      &
                                ! min mmr of O3 required for oxidn
         &     THOLD                  ! threshold for cloud liquid water
    !
    !
    PARAMETER (                                                       &
         &           EVAPTAU = 300.0_r8,                                       &
                                ! secs  (=5 mins)
         &             NUCTAU = 30.0_r8,                                       &
                                ! secs
         &       DIFFUSE_AIT = 1.7134E-9_r8,                                   &
                                ! sq m/s
         &        K_SO2OH_HI = 2.0E-12_r8,                                     &
                                ! cc/mcl/s from STOCHEM model
         &           K_DMS_OH = 9.1E-12_r8,                                    &
                                ! cc/mcl/s
         &       K4_CH3SO2_O3 = 1.0E-14_r8,                                    &
                                ! cc/mcl/s
         &      K5_CH3SO3_HO2 = 4.0E-11_r8,                                    &
         &             RMM_O3 = 4.8E-2_r8,                                     &
                                ! kg/mole
         &          BRAT_SO2 = 0.9_r8,                                         &
         &           BRAT_MSA = 1.0_r8-BRAT_SO2,                               &
         &           AVOGADRO = 6.022E23_r8,                                   &
                                ! per mole
         
         &           RMM_S = 3.20E-2_r8,                                       &
                                ! kg/mole
         &           RMM_H2O2 = 3.40E-2_r8,                                    &
                                ! kg/mole
         &           RMM_HO2 = 3.30E-2_r8,                                     &
                                ! kg/mole
         &            RMM_AIR = 2.896E-2_r8,                                   &
                                ! kg/mole
         &              RMM_W = 1.8E-2_r8,                                     &
                                ! kg/mole
         &        RELM_S_H2O2 = 3.206_r8/3.40_r8,                                 &
         &           RELM_S_2N = 3.206_r8/2.80_r8,                                &
         &               PARH = 3.3_r8,                                        &
         &                K1 = 4.0E-31_r8,                                     &
                                ! (cc/mcl)2/s from STOCHEM
         &                 T1 = 300.0_r8,                                      &
                                ! K
         &                FC = 0.45_r8,                                        &
                                ! from STOCHEM model
         &              FAC1 = 1.1904_r8,                                      &
                                ! 0.75-1.27*LOG10(FC) from STOCHEM
         &                 K2 = 2.2E-13_r8,                                    &
                                ! cc/mcl/s
         &                 K3 = 1.9E-33_r8,                                    &
                                ! (cc/mcl)2/s
         &                 K4 = 1.4E-21_r8,                                    &
                                ! cc/mcl
         &                 T2 = 600.0_r8,                                      &
                                ! K
         &                 T3 = 890.0_r8,                                      &
                                ! K
         &                 T4 = 2200.0_r8,                                     &
                                ! K
         &           CLOUDTAU = 1.08E4_r8,                                     &
                                ! secs (=3 hours)
         &            CHEMTAU = 9.0E2_r8,                                      &
                                ! secs (=15 mins)
         &              O3_MIN = 1.6E-8_r8,                                    &
                                !(kg/kg, equiv. 10ppbv)
         &              THOLD = 1.0E-8_r8                                      &
                                ! kg/kg
         &          )
    !
    REAL(KIND=r8) RAD_AIT,                                                     &
                                ! median radius of Aitken mode particles
         &     DIAM_AIT,                                                    &
                                !   "    diameter    "
         &     RAD_ACC,                                                     &
                                ! median radius of acccumulation mode
         &     DIAM_ACC,                                                    &
                                !   "    diameter    "
         &     CHI,                                                         &
                                ! mole fraction of S in particle
         &     RHO_SO4,                                                     &
                                ! density of  SO4 particle
         &     SIGMA,                                                       &
                                ! standard devn of particle size distn
         !                                 for accumulation mode
         &     E_PARM,                                                      &
                                ! param relating size distns of Ait & Acc
         &     NUM_STAR         ! threshold concn of accu mode particles
    !  below which PSI=1
    !
    PARAMETER (                                                       &
         &           RAD_AIT = 6.5E-9_r8,                                      &
                                ! m
         &          DIAM_AIT = 2.0_r8*RAD_AIT,                                 &
         &           RAD_ACC = 95.0E-9_r8,                                     &
                                ! m
         &          DIAM_ACC = 2.0_r8*RAD_ACC,                                 &
         &               CHI = 32.0_r8/132.0_r8,                                  &
         &           RHO_SO4 = 1769.0_r8,                                      &
                                ! kg/m3
         &             SIGMA = 1.4_r8,                                         &
         &            E_PARM = 0.9398_r8,                                      &
         &          NUM_STAR = 1.0E6_r8                                        &
                                ! m-3
         &          )
    !
    REAL(KIND=r8) BOLTZMANN       !Boltzmanns constant.
    REAL(KIND=r8) MFP_REF         !Reference value of mean free path.
    REAL(KIND=r8) TREF_MFP        !Reference temperature for mean free path.
    REAL(KIND=r8) PREF_MFP        !Reference pressure for mean free path.
    REAL(KIND=r8) SIGMA_AIT       !Geometric standard deviation of the Aitken
    !                             mode distribution.
    !
    PARAMETER (BOLTZMANN = 1.3804E-23_r8)  ! J K-1
    PARAMETER (MFP_REF = 6.6E-8_r8                                       &
                                ! m
         &        ,  TREF_MFP = 293.15_r8                                      &
                                ! K
         &        ,  PREF_MFP = 1.01325E5_r8)    ! Pa
    PARAMETER (SIGMA_AIT = 1.30_r8)
    !
    !*---------------------------------------------------------------------
    !
    !       EXTERNAL VGRAV
    !
    !... epsilon() is defined as almost negligible, so eps/100 is negligible
    !
    SMALLP = EPSILON(1.0_r8) / 100.0_r8
    !
    !...calc stokes number, schmidt number and finally resistance
    !
    LEV1=1

    DO IDIV=1,NDIV
       ! DEPENDS ON: vgrav
       CALL VGRAV(                                                     &
            &  nCols,mskant,LEV1,DREP(IDIV),RHOP,PSTAR,TSTAR,               &
            &  VSTOKES1,CCF,ETAA                                               &
            &  )

       !CDIR NOVECTOR
       DO I= 1,nCols
          IF(mskant(i) == 1_i8)THEN

             NSCHMIDT(I)=3.0_r8*PI*ETAA(I)*ETAA(I)*DREP(IDIV)/         &
                  &       (RHOSTAR(I)*BOLTZMANN*TSTAR(I)*CCF(I))
             NSTOKES(I)=VSTOKES1(I)*CD_STD_DUST(I)*RHOSTAR(I)*   &
                  &       VSHR(I)*VSHR(I)/(ETAA(I)*G)
             ! Avoid underflow in Stokes term by setting to zero if 
             ! negligible compared to Schmidt term, i.e., if NSTOKES
             ! is too small.
             IF ( 3.0_r8 / NSTOKES(I) <                                   &
                  - LOG10( SMALLP *NSCHMIDT(I)**(-2.0_r8/3.0_r8) ) ) THEN
                STOKES_EXP = 10.0_r8**(-3.0_r8/NSTOKES(I))
             ELSE
                STOKES_EXP = 0.0_r8
             ENDIF
             R_B_DUST(I,IDIV)=1.0_r8/( SQRT(CD_STD_DUST(I)) *            &
                  &       (NSCHMIDT(I)**(-2.0_r8/3.0_r8)+STOKES_EXP) )
          END IF

       ENDDO !nCols
    ENDDO !NDIV

    RETURN
  END SUBROUTINE DUSTRESB


  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !    Subroutine VGRAV -------------------------------------------------
  !
  ! Purpose: To calculate the gravitational sedimentation velocity of
  !          tracer particles according to Stoke's law, including the
  !          Cunningham correction factor.
  !
  ! Current owners of code:                 S Woodward, M Woodage
  !
  ! History:
  ! Version    Date     Comment
  ! -------    ----     -------
  !   4.4    03/10/97   Original code        S Woodward, M Woodage
  !   5.5    12/02/03   Updated for vn 5.5   S Woodward
  !
  !
  ! Code description:
  !  Language: FORTRAN77 + extensions
  !  Programming standard: UMDP 3 Vn 6
  !
  ! System components covered:
  !
  ! System task:
  !
  !Documentation: Ref. Pruppacher & Klett
  !                    Microphysics of clouds & ppn    1978,1980 edns.
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE VGRAV(                                                 &
       & nCols,mskant,NLEVS,DIAM,RHOP,P,T,                             &
       & VSTOKES,CCF,ETAA)




    !
    IMPLICIT NONE
    !
    !
    INTEGER nCols         !IN row length
    INTEGER NLEVS              !IN number of levels
    !
    INTEGER(KIND=i8) :: mskant(nCols)

    REAL(KIND=r8) DIAM                  !IN particle diameter
    REAL(KIND=r8) RHOP                  !IN particles density
    REAL(KIND=r8) P(nCols,NLEVS)!IN pressure
    REAL(KIND=r8) T(nCols,NLEVS)!IN temperature
    !
    REAL(KIND=r8) VSTOKES(nCols,NLEVS) !OUT sedimentation velocity
    REAL(KIND=r8) ETAA(nCols,NLEVS)!OUT viscosity of air
    REAL(KIND=r8) CCF(nCols,NLEVS) !OUT cunningham correction factor
    !
    !
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_O_DG_C-----------------------------------
    ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
    ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
    ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

    REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
    REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
    REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

    !*----------------------------------------------------------------------
    !-------------------COMDECK C_SULCHM--------------------------------
    ! Parameters for Sulphur Cycle Chemistry
    REAL(KIND=r8)                                                              &
         &     EVAPTAU,                                                     &
                                ! timescale for dissolved SO4 to evaporate
         &     NUCTAU,                                                      &
                                ! timescale for accumulation mode particles
         !                           to nucleate once they enter a cloud.
         &     DIFFUSE_AIT,                                                 &
                                ! diffusion coefficient of Aitken particles
         &     K_SO2OH_HI,                                                  &
                                ! high pressure reaction rate limit
         &     K_DMS_OH,                                                    &
                                ! reaction rate for DMS+OH  cc/mcl/s
         &      K4_CH3SO2_O3,                                               &
                                ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
         &      K5_CH3SO3_HO2,                                              &
                                ! Rate coeff for CH3SO3+HO2 -> MSA+O2
         &      RMM_O3,                                                     &
                                ! relative molecular mass O3
         &     BRAT_SO2,                                                    &
                                ! branching ratio for SO2 in DMS oxidn
         &     BRAT_MSA,                                                    &
                                ! branching ratio for MSA in DMS oxidn
         &     AVOGADRO,                                                    &
                                ! no. of molecules in 1 mole
         &     RMM_S,                                                       &
                                ! relative molecular mass S kg/mole
         &     RMM_H2O2,                                                    &
                                ! relative molecular mass H2O2 kg/mole
         &     RMM_HO2,                                                     &
                                ! relative molecular mass HO2 kg/mole
         &     RMM_AIR,                                                     &
                                ! relative molecular mass dry air
         &     RMM_W,                                                       &
                                ! relative molecular mass water
         &     RELM_S_H2O2,                                                 &
                                ! rel atomic mass sulphur/RMM_H2O2
         &     RELM_S_2N,                                                   &
                                ! rel atomic mass Sulphur/2*Nitrogen
         &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
         &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
         &     T1,                                                          &
                                !
         &     FC,                                                          &
                                ! parameters for interpolation between
         &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
         &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
         &     T2,T3,T4,                                                    &
                                !
         &     CLOUDTAU,                                                    &
                                ! air parcel lifetime in cloud
         &     CHEMTAU,                                                     &
                                ! chem lifetime in cloud before oxidn
         &     O3_MIN,                                                      &
                                ! min mmr of O3 required for oxidn
         &     THOLD                  ! threshold for cloud liquid water
    !
    !
    PARAMETER (                                                       &
         &           EVAPTAU = 300.0_r8,                                       &
                                ! secs  (=5 mins)
         &             NUCTAU = 30.0_r8,                                       &
                                ! secs
         &       DIFFUSE_AIT = 1.7134E-9_r8,                                   &
                                ! sq m/s
         &        K_SO2OH_HI = 2.0E-12_r8,                                     &
                                ! cc/mcl/s from STOCHEM model
         &           K_DMS_OH = 9.1E-12_r8,                                    &
                                ! cc/mcl/s
         &       K4_CH3SO2_O3 = 1.0E-14_r8,                                    &
                                ! cc/mcl/s
         &      K5_CH3SO3_HO2 = 4.0E-11_r8,                                    &
         &             RMM_O3 = 4.8E-2_r8,                                     &
                                ! kg/mole
         &          BRAT_SO2 = 0.9_r8,                                         &
         &           BRAT_MSA = 1.0_r8-BRAT_SO2,                               &
         &           AVOGADRO = 6.022E23_r8,                                   &
                                ! per mole
         
         &           RMM_S = 3.20E-2_r8,                                       &
                                ! kg/mole
         &           RMM_H2O2 = 3.40E-2_r8,                                    &
                                ! kg/mole
         &           RMM_HO2 = 3.30E-2_r8,                                     &
                                ! kg/mole
         &            RMM_AIR = 2.896E-2_r8,                                   &
                                ! kg/mole
         &              RMM_W = 1.8E-2_r8,                                     &
                                ! kg/mole
         &        RELM_S_H2O2 = 3.206_r8/3.40_r8,                                 &
         &           RELM_S_2N = 3.206_r8/2.80_r8,                                &
         &               PARH = 3.3_r8,                                        &
         &                K1 = 4.0E-31_r8,                                     &
                                ! (cc/mcl)2/s from STOCHEM
         &                 T1 = 300.0_r8,                                      &
                                ! K
         &                FC = 0.45_r8,                                        &
                                ! from STOCHEM model
         &              FAC1 = 1.1904_r8,                                      &
                                ! 0.75-1.27*LOG10(FC) from STOCHEM
         &                 K2 = 2.2E-13_r8,                                    &
                                ! cc/mcl/s
         &                 K3 = 1.9E-33_r8,                                    &
                                ! (cc/mcl)2/s
         &                 K4 = 1.4E-21_r8,                                    &
                                ! cc/mcl
         &                 T2 = 600.0_r8,                                      &
                                ! K
         &                 T3 = 890.0_r8,                                      &
                                ! K
         &                 T4 = 2200.0_r8,                                     &
                                ! K
         &           CLOUDTAU = 1.08E4_r8,                                     &
                                ! secs (=3 hours)
         &            CHEMTAU = 9.0E2_r8,                                      &
                                ! secs (=15 mins)
         &              O3_MIN = 1.6E-8_r8,                                    &
                                !(kg/kg, equiv. 10ppbv)
         &              THOLD = 1.0E-8_r8                                      &
                                ! kg/kg
         &          )
    !
    REAL(KIND=r8) RAD_AIT,                                                     &
                                ! median radius of Aitken mode particles
         &     DIAM_AIT,                                                    &
                                !   "    diameter    "
         &     RAD_ACC,                                                     &
                                ! median radius of acccumulation mode
         &     DIAM_ACC,                                                    &
                                !   "    diameter    "
         &     CHI,                                                         &
                                ! mole fraction of S in particle
         &     RHO_SO4,                                                     &
                                ! density of  SO4 particle
         &     SIGMA,                                                       &
                                ! standard devn of particle size distn
         !                                 for accumulation mode
         &     E_PARM,                                                      &
                                ! param relating size distns of Ait & Acc
         &     NUM_STAR         ! threshold concn of accu mode particles
    !  below which PSI=1
    !
    PARAMETER (                                                       &
         &           RAD_AIT = 6.5E-9_r8,                                      &
                                ! m
         &          DIAM_AIT = 2.0_r8*RAD_AIT,                                 &
         &           RAD_ACC = 95.0E-9_r8,                                     &
                                ! m
         &          DIAM_ACC = 2.0_r8*RAD_ACC,                                 &
         &               CHI = 32.0_r8/132.0_r8,                                  &
         &           RHO_SO4 = 1769.0_r8,                                      &
                                ! kg/m3
         &             SIGMA = 1.4_r8,                                         &
         &            E_PARM = 0.9398_r8,                                      &
         &          NUM_STAR = 1.0E6_r8                                        &
                                ! m-3
         &          )
    !
    REAL(KIND=r8) BOLTZMANN       !Boltzmanns constant.
    REAL(KIND=r8) MFP_REF         !Reference value of mean free path.
    REAL(KIND=r8) TREF_MFP        !Reference temperature for mean free path.
    REAL(KIND=r8) PREF_MFP        !Reference pressure for mean free path.
    REAL(KIND=r8) SIGMA_AIT       !Geometric standard deviation of the Aitken
    !                             mode distribution.
    !
    PARAMETER (BOLTZMANN = 1.3804E-23_r8)  ! J K-1
    PARAMETER (MFP_REF = 6.6E-8_r8                                       &
                                ! m
         &        ,  TREF_MFP = 293.15_r8                                      &
                                ! K
         &        ,  PREF_MFP = 1.01325E5_r8)    ! Pa
    PARAMETER (SIGMA_AIT = 1.30_r8)
    !
    !*---------------------------------------------------------------------
    !C_DUSTGRAV.............................................................
    ! Description:
    ! Contains parameters for mineral dust gravitational settling
    ! Current Code Owner: Stephanie Woodward
    !
    ! History:
    ! Version  Date     Comment
    ! -------  ----     -------
    !  1      12/02/03  Original Code.   Stephanie Woodward
    !
    REAL(KIND=r8),PARAMETER :: ACCF=1.257_r8 ! Cunningham correction factor term A
    REAL(KIND=r8),PARAMETER :: BCCF=0.4_r8 ! Cunningham correction factor term B
    REAL(KIND=r8),PARAMETER :: CCCF=-1.1_r8 ! Cunningham correction factor term C
    !.......................................................................
    !
    ! local variables
    !
    INTEGER ILEV               !LOC loop counter for levels
    INTEGER I                  !LOC loop counter
    !
    REAL(KIND=r8) TC(nCols)   !LOC temperature in deg C
    REAL(KIND=r8) LAMDAA(nCols)!LOC mean free path of particle
    REAL(KIND=r8) ALPHACCF(nCols)!LOC
    !
    ! Calculate viscosity of air (Pruppacher & Klett p.323)
    DO ILEV=1,NLEVS
       DO I = 1,nCols
          IF(mskant(i) == 1_i8)THEN

             TC(I)=T(I,ILEV)-ZERODEGC
             IF (TC(I)  >=  0.0_r8) THEN
                ETAA(I,ILEV)=(1.718_r8+0.0049_r8*TC(I))*1.E-5_r8
             ELSE
                ETAA(I,ILEV)=(1.718_r8+0.0049_r8*TC(I)-1.2E-5_r8*TC(I)*TC(I))*1.E-5_r8
             ENDIF
             !
          END IF
       ENDDO !nCols
    ENDDO !NLEVS
    !
    DO ILEV=1,NLEVS
       DO I=1,nCols
          IF(mskant(i) == 1_i8)THEN

             !
             ! Calculate mean free path of particle (Pruppacher & Klett p.323)
             LAMDAA(I)=MFP_REF*PREF_MFP*T(I,ILEV)/                    &
                  &      (P(I,ILEV)*TREF_MFP)
             ! Calculate Cunningham correction factor(Pruppacher & Klett p.361)
             ALPHACCF(I)=ACCF+BCCF*EXP(CCCF*DIAM*.5_r8/LAMDAA(I))
             CCF(I,ILEV)=(1.0_r8+ALPHACCF(I)*LAMDAA(I)/(0.5_r8*DIAM))
             ! Calculate sedimentation velocity (Pruppacher & Klett p.362)
             VSTOKES(I,ILEV)=CCF(I,ILEV)*(DIAM*DIAM*G*RHOP)/          &
                  &             (18.0_r8*ETAA(I,ILEV))
             !
          END IF
       ENDDO !nCols
    ENDDO !NLEV
    !
    RETURN
  END SUBROUTINE VGRAV



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
!!!  SUBROUTINES SFL_INT_SEA AND SFL_INT_LAND--------------------------
!!!
!!!  Purpose: To calculate interpolation coefficients for 10m winds
!!!           and 1.5m temperature/specific humidity diagnostics.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2  15/11/00   New Deck         M. Best
  !    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!  Programming standard:
!!!
!!!  Logical component covered: Part of P243.
!!!
!!!  System Task:
!!!
!!!  External Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
  !*L  Arguments :-
  SUBROUTINE SFL_INT_SEA (                                          &
       & nCols,mskant,FLANDG                                                     &
       &,VSHR,CD,CH,Z0M,Z0H                                               &
       &,RECIP_L_MO,V_S                                                   &
       &,SU10,SV10,ST1P5,SQ1P5,CDR10M,CHR1P5M)
    IMPLICIT NONE

    INTEGER :: nCols      ! IN Number of X points?
    INTEGER(KIND=i8) :: mskant(nCols) 
    REAL(KIND=r8)                                                              &
         & FLANDG(nCols)                                          &
         !                           ! IN Land fraction.
         &,Z0M(nCols)                                             &
                                ! IN Roughness length for momentum (m).
         &,Z0H(nCols)                                             &
                                ! IN Roughness length for heat and
         !                           !    moisture (m).
         &,VSHR(nCols)                                            &
                                ! IN Wind speed difference between the
         !                           !    surface and the lowest wind level in
         !                           !    the atmosphere (m/s).
         &,CD(nCols)                                              &
                                ! IN Surface drag coefficient.
         &,CH(nCols)                                              &
                                ! IN Surface transfer coefficient for heat
         !                           !    and moisture.
         &,RECIP_L_MO(nCols)                                      &
         !                           ! IN Reciprocal of the Monin-Obukhov
         !                           ! length (m)
         &,V_S(nCols) ! IN Surface layer scaling velocity
    !                           !    including orographic form drag (m/s).

    LOGICAL                                                           &
         & SU10                                                             &
                                ! IN 10m U-wind diagnostic flag
         &,SV10                                                             &
                                ! IN 10m V-wind diagnostic flag
         &,ST1P5                                                            &
                                ! IN screen temp diagnostic flag
         &,SQ1P5
    ! IN screen specific humidity
    !                                !    diagnostic flag
    ! Output variables
    !
    REAL(KIND=r8) :: CDR10M(nCols)              &
         !                                ! OUT interpolation coefficicent for
         !                                ! 10m wind
         &,CHR1P5M(nCols)  ! OUT Interpolation coefficient for
    !                                !     1.5m temperature

    !*
    !*L---------------------------------------------------------------------


    !  External routines called :-
    !      EXTERNAL PHI_M_H_SEA
    !      EXTERNAL TIMER

    !*
    !*L---------------------------------------------------------------------
    !    Local and other symbolic constants :-
    ! C_VKMAN start
    REAL(KIND=r8),PARAMETER:: VKMAN=0.4_r8 ! Von Karman's constant
    ! C_VKMAN end
    REAL(KIND=r8) Z_OBS_TQ,Z_OBS_WIND
    PARAMETER (                                                       &
         & Z_OBS_TQ = 1.5_r8                                                   &
                                ! Height of screen observations of temperature
         !                        ! and humidity.
         &,Z_OBS_WIND = 10.0_r8                                                &
                                ! Height of surface wind observations.
         &)
    !
    !  Define local storage.
    !
    !  (a) Local work arrays.
    !
    REAL(KIND=r8)                                                              &
         & Z_WIND(nCols)                                          &
                                ! Height of wind observations.
         &,Z_TEMP(nCols)                                          &
                                ! Height of temperature and humidity
         !                                ! observations.
         &,PHI_M_OBS(nCols)                                       &
                                ! Monin-Obukhov stability function for
         !                                ! momentum integrated to the wind
         !                                ! observation height.
         &,PHI_H_OBS(nCols)! Monin-Obukhov stability function for
    !                                ! scalars integrated to their
    !                                ! observation height.
    !
    !  (b) Scalars.
    !
    INTEGER                                                           &
         & I     ! Loop counter (horizontal field index).
    !*
    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('SFL_INT   ',3)
    !      ENDIF
    !
    !-----------------------------------------------------------------------
    !! 1. If diagnostics required calculate M-O stability functions at
    !!    observation heights.
    !-----------------------------------------------------------------------

    IF (SU10 .OR. SV10 .OR. ST1P5 .OR. SQ1P5) THEN
       DO I=1,nCols
          IF ( FLANDG(I) <  1.0_r8 ) THEN
             Z_WIND(I) = Z_OBS_WIND
             Z_TEMP(I) = Z_OBS_TQ + Z0H(I) - Z0M(I)
          END IF
       ENDDO
       ! DEPENDS ON: phi_m_h_sea 

       CALL PHI_M_H_SEA (nCols,mskant,FLANDG,                       &
            &                    RECIP_L_MO,Z_WIND,Z_TEMP,Z0M,Z0H,             &
            &                    PHI_M_OBS,PHI_H_OBS)
    ENDIF

    !-----------------------------------------------------------------------
    !! 2. If diagnostics required calculate interpolation coefficient
    !!    for 1.5m screen temperature and specific humidity.
    !-----------------------------------------------------------------------
    !
    IF (ST1P5 .OR. SQ1P5) THEN
       DO I=1,nCols
          IF ( FLANDG(I) <  1.0_r8 ) THEN
             CHR1P5M(I) = CH(I) * VSHR(I) *                        &
                  &                          PHI_H_OBS(I)/(VKMAN*V_S(I))
          ENDIF
       ENDDO
    ENDIF
    !
    !-----------------------------------------------------------------------
    !! 3. If diagnostics required calculate interpolation coefficient
    !!    for 10m winds.
    !-----------------------------------------------------------------------
    !
    IF ( SU10 .OR. SV10 ) THEN
       DO I=1,nCols
          IF ( FLANDG(I) <  1.0_r8 ) THEN
             CDR10M(I) = (1.0_r8-FLANDG(I)) * CD(I) * VSHR(I) *     &
                  &                       PHI_M_OBS(I)/(VKMAN*V_S(I))
          ENDIF
       ENDDO
    ENDIF
    !
    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('SFL_INT ',4)
    !      ENDIF
    RETURN
  END SUBROUTINE SFL_INT_SEA

!!!
!!!---------------------------------------------------------------------
  !*L  Arguments :-


  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
!!!  SUBROUTINES STDEV1_SEA and STDEV1_LAND ----------------------------
!!!
!!!  Purpose: Calculate the standard deviations of layer 1 turbulent
!!!           fluctuations of temperature and humidity using approximate
!!!           formulae from first order closure.
!!!
!!!    Model            Modification history
!!!   version  date
!!!    5.2   15/11/00   New Deck         M. Best
  !    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!    Programming standard:
!!!
!!!  -------------------------------------------------------------------
  !

!!!  SUBROUTINE STDEV1_SEA ---------------------------------------------
!!!  Layer 1 standard deviations for sea and sea-ice
!!!  -------------------------------------------------------------------
  SUBROUTINE STDEV1_SEA (&
       nCols   , &!INTEGER, INTENT(IN   ) :: nCols! IN Number of X points?
       mskant       , &
       FLANDG       , &!REAL(KIND=r8)   , INTENT(IN   ) :: FLANDG(nCols,ROWS)! IN Land fraction.
       BQ_1         , &!REAL(KIND=r8)   , INTENT(IN   ) :: BQ_1(nCols,ROWS)  ! IN Buoyancy parameter.
       BT_1         , &!REAL(KIND=r8)   , INTENT(IN   ) :: BT_1(nCols,ROWS)  ! IN Buoyancy parameter.
       FQW_1        , &!REAL(KIND=r8)   , INTENT(IN   ) :: FQW_1(nCols,ROWS) ! IN Surface flux of QW.
       FTL_1        , &!REAL(KIND=r8)   , INTENT(IN   ) :: FTL_1(nCols,ROWS) ! IN Surface flux of TL.
       ICE_FRACT    , &!REAL(KIND=r8)   , INTENT(IN   ) :: ICE_FRACT(nCols,ROWS) ! IN Fraction of gridbox which is sea-ice.
       RHOKM_1      , &!REAL(KIND=r8)   , INTENT(IN   ) :: RHOKM_1(nCols,ROWS) ! IN Surface momentum exchange
       RHOSTAR      , &!REAL(KIND=r8)   , INTENT(IN   ) :: RHOSTAR(nCols,ROWS) ! IN Surface air density.
       VSHR         , &!REAL(KIND=r8)   , INTENT(IN   ) :: VSHR(nCols,ROWS)    ! IN Magnitude of surface-
       Z0MSEA       , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z0MSEA(nCols,ROWS)  ! IN Sea roughness length.
       Z0_ICE       , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z0_ICE(nCols,ROWS)  ! IN Sea-ice roughness length.
       Z1_TQ        , &!REAL(KIND=r8)   , INTENT(IN   ) :: Z1_TQ(nCols,ROWS)     ! IN Height of lowest tq level.
       Q1_SD        , &!REAL(KIND=r8)   , INTENT(out  ) :: Q1_SD(nCols,ROWS)  ! OUT Standard deviation of
       T1_SD         &!REAL(KIND=r8)   , INTENT(out  ) :: T1_SD(nCols,ROWS)     ! OUT Standard deviation of
       )

    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: nCols! IN Number of X points?
    INTEGER(KIND=i8), INTENT(IN   ) ::  mskant(nCols)

    REAL(KIND=r8)   , INTENT(IN   ) :: FLANDG(nCols)! IN Land fraction.
    REAL(KIND=r8)   , INTENT(IN   ) :: BQ_1(nCols)  ! IN Buoyancy parameter.
    REAL(KIND=r8)   , INTENT(IN   ) :: BT_1(nCols)  ! IN Buoyancy parameter.
    REAL(KIND=r8)   , INTENT(IN   ) :: FQW_1(nCols) ! IN Surface flux of QW.
    REAL(KIND=r8)   , INTENT(IN   ) :: FTL_1(nCols) ! IN Surface flux of TL.
    REAL(KIND=r8)   , INTENT(IN   ) :: ICE_FRACT(nCols) ! IN Fraction of gridbox which is sea-ice.
    REAL(KIND=r8)   , INTENT(IN   ) :: RHOKM_1(nCols) ! IN Surface momentum exchange
    !                                 !    coefficient.
    REAL(KIND=r8)   , INTENT(IN   ) :: RHOSTAR(nCols) ! IN Surface air density.
    REAL(KIND=r8)   , INTENT(IN   ) :: VSHR(nCols)    ! IN Magnitude of surface-
    !                                 !    to-lowest-level wind shear.
    REAL(KIND=r8)   , INTENT(IN   ) :: Z0MSEA(nCols)  ! IN Sea roughness length.
    REAL(KIND=r8)   , INTENT(IN   ) :: Z0_ICE(nCols)  ! IN Sea-ice roughness length.
    REAL(KIND=r8)   , INTENT(IN   ) :: Z1_TQ(nCols)     ! IN Height of lowest tq level.

    REAL(KIND=r8)   , INTENT(out  ) :: Q1_SD(nCols)  ! OUT Standard deviation of
    !                                 !     turbulent fluctuations of
    !                                 !     surface layer specific
    !                                 !     humidity (kg/kg).
    REAL(KIND=r8)   , INTENT(out  ) :: T1_SD(nCols)     ! OUT Standard deviation of
    !                                 !     turbulent fluctuations of
    !                                 !     surface layer temperature (K).


    !  External routines called :-
    !      EXTERNAL TIMER


    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------

    !  Workspace --------------------------------------------------------
    INTEGER                                                           &
         & I                   ! Loop counter (horizontal field index).
    REAL(KIND=r8)      :: VS &
                                ! Surface layer friction velocity
         &,VSF1_CUBED                                                       &
                                ! Cube of surface layer free convective
         !                            ! scaling velocity
         &,WS1                                                              &
                                ! Turbulent velocity scale for surface
         !                            ! layer
         &,Z0                    ! Roughness length

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('STDEV1  ',3)
    !      ENDIF

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( FLANDG(I) <  1.0_r8 ) THEN

             Z0 = Z0MSEA(I)
             IF ( ICE_FRACT(I)  >   0.0_r8 ) Z0 = Z0_ICE(I)
             VS = SQRT ( RHOKM_1(I)/RHOSTAR(I) * VSHR(I) )
             VSF1_CUBED = 1.25_r8*G*(Z1_TQ(I) + Z0) *                       &
                  &                ( BT_1(I)*FTL_1(I) +                          &
                  &                   BQ_1(I)*FQW_1(I) )/RHOSTAR(I)
             IF ( VSF1_CUBED  >   0.0_r8 ) THEN
                WS1 = ( VSF1_CUBED + VS*VS*VS ) ** (1.0_r8/3.0_r8)
                T1_SD(I) = MAX ( 0.0_r8 ,                                    &
                     &          (1.0_r8-FLANDG(I))*1.93_r8*FTL_1(I) / (RHOSTAR(I)*WS1) )
                Q1_SD(I) = MAX ( 0.0_r8 ,                                    &
                     &          (1.0_r8-FLANDG(I))*1.93_r8*FQW_1(I) / (RHOSTAR(I)*WS1) )

             ENDIF

          ENDIF
       END IF
    ENDDO

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('STDEV1  ',4)
    !      ENDIF

    RETURN
  END SUBROUTINE STDEV1_SEA

!!!  SUBROUTINE STDEV1_LAND --------------------------------------------
!!!  Layer 1 standard deviations for land tiles
!!!  -------------------------------------------------------------------




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
!!!   SUBROUTINES FCDCH_SEA AND FCDCH_LAND-----------------------------
!!!
!!!  Purpose: Calculate surface transfer coefficients at one or more
!!!           gridpoints.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  5.2   15/11/00   New Deck         M. Best
  !    5.3  25/04/01  Add coastal tiling. Nic Gedney
  !    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
  !  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!

  !     SUBROUTINE FCDCH_SEA---------------------------------------------
  !
  !     Transfer coefficients for sea, sea-ice and leads
  !
  !     -----------------------------------------------------------------
  !!  Arguments:---------------------------------------------------------
  SUBROUTINE FCDCH_SEA(&
       nCols      , & !INTEGER, INTENT(IN  ) :: nCols              ! IN Number of X points?
       mskant     , &
       FLANDG     , & !REAL(KIND=r8)   , INTENT(IN  ) :: FLANDG(nCols,ROWS) ! IN Land fraction
       DB         , & !REAL(KIND=r8)   , INTENT(IN  ) :: DB    (nCols,ROWS) !
       VSHR       , & !REAL(KIND=r8)   , INTENT(IN  ) :: VSHR  (nCols,ROWS) 
       Z0M        , & !REAL(KIND=r8)   , INTENT(IN  ) :: Z0M   (nCols,ROWS) 
       Z0H        , & !REAL(KIND=r8)   , INTENT(IN  ) :: Z0H   (nCols,ROWS) 
       ZH         , & !REAL(KIND=r8)   , INTENT(IN  ) :: ZH   (nCols,ROWS)  
       Z1_UV      , & !REAL(KIND=r8)   , INTENT(IN  ) :: Z1_UV(nCols,ROWS)  
       Z1_TQ      , & !REAL(KIND=r8)   , INTENT(IN  ) :: Z1_TQ(nCols,ROWS)  
       CDV        , & !REAL(KIND=r8)   , INTENT(OUT ) :: CDV   (nCols,ROWS) 
       CHV        , & !REAL(KIND=r8)   , INTENT(OUT ) :: CHV   (nCols,ROWS)     
       V_S        , & !REAL(KIND=r8)   , INTENT(OUT ) :: V_S(nCols,ROWS)       
       RECIP_L_MO  & !REAL(KIND=r8)   , INTENT(OUT ) :: RECIP_L_MO(nCols,ROWS) 
       )
    IMPLICIT NONE

    INTEGER, INTENT(IN  ) :: nCols              ! IN Number of X points?

    INTEGER(KIND=i8), INTENT(IN  ) :: mskant(nCols)
    REAL(KIND=r8)   , INTENT(IN  ) :: FLANDG(nCols) ! IN Land fraction
    REAL(KIND=r8)   , INTENT(IN  ) :: DB    (nCols) ! IN Buoyancy difference between surface
    ! and lowest temperature and humidity
    ! level in the atmosphere (m/s^2).
    REAL(KIND=r8)   , INTENT(IN  ) :: VSHR  (nCols) ! IN Wind speed difference between the
    ! surface and the lowest wind level in
    ! the atmosphere (m/s).
    REAL(KIND=r8)   , INTENT(IN  ) :: Z0M   (nCols) ! IN Roughness length for momentum
    ! transport (m).
    REAL(KIND=r8)   , INTENT(IN  ) :: Z0H   (nCols) ! IN Roughness length for heat and
    ! moisture (m).
    REAL(KIND=r8)   , INTENT(IN  ) :: ZH   (nCols)  ! IN Depth of boundary layer (m).
    REAL(KIND=r8)   , INTENT(IN  ) :: Z1_UV(nCols)  ! IN Height of lowest wind level (m).
    REAL(KIND=r8)   , INTENT(IN  ) :: Z1_TQ(nCols)  ! IN Height of lowest temperature and
    ! humidity level (m).

    REAL(KIND=r8)   , INTENT(OUT ) :: CDV   (nCols)     ! OUT Surface transfer coefficient for
    ! momentum including orographic form drag (m/s).
    REAL(KIND=r8)   , INTENT(OUT ) :: CHV   (nCols)     ! OUT Surface transfer coefficient for
    !     heat, moisture & other scalars (m/s).
    REAL(KIND=r8)   , INTENT(OUT ) :: V_S(nCols)        ! OUT Surface layer scaling velocity
    !                                                          !     including orographic form drag (m/s).
    REAL(KIND=r8)   , INTENT(OUT ) :: RECIP_L_MO(nCols) ! OUT Reciprocal of the Monin-Obukhov length
    !     (m^-1).

    !*L  Workspace usage----------------------------------------------------
    !
    !     Local work arrays.
    !
    REAL(KIND=r8)                                                              &
         & PHI_M(nCols)                                           &
                                ! Monin-Obukhov stability function for
         !                             ! momentum integrated to the model's
         !                             ! lowest wind level.
         &,PHI_H(nCols) ! Monin-Obukhov stability function for
    !                             ! scalars integrated to the model's lowest
    !                             ! temperature and humidity level.
    !
    !*----------------------------------------------------------------------

    !      EXTERNAL PHI_M_H_SEA
    !      EXTERNAL TIMER

    !*----------------------------------------------------------------------
    !  Common and local constants.
    ! C_VKMAN start
    REAL(KIND=r8),PARAMETER:: VKMAN=0.4_r8 ! Von Karman's constant
    ! C_VKMAN end
    REAL(KIND=r8),PARAMETER:: CTRL_TUNABLE=1.0_r8
                                ! Tunable parameter in the surface layer scaling
         !                   ! velocity formula (multiplying the turbulent
         !                   ! convective scaling velocity).

    REAL(KIND=r8) BETA,THIRD
    PARAMETER (                                                       &
         & BETA=0.09_r8,                                                       &
                             ! 0.009 original pkubota |Tunable parameter in the surface layer scaling
         !                   ! velocity formula (multiplying the turbulent
         !                   ! convective scaling velocity).
         & THIRD=1.0_r8/3.0_r8                                                      &
                                ! One third.
         &)
    INTEGER, PARAMETER ::  N_ITS=10  ! 5 Number of iterations for Monin-Obukhov length
    !                   ! and stability functions.
    !
    !  Define local variables
    !
    INTEGER I   ! Loop counters; horizontal field index.
    INTEGER IT    ! Iteration loop counter.

    REAL(KIND=r8)                                                              &
         & B_FLUX                                                           &
                                ! Surface bouyancy flux over air density.
         &,U_S                                                              &
                                ! Surface friction velocity (effective value).
         &,W_S          ! Surface turbulent convective scaling velocity.

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('FCDCH   ',3)
    !      ENDIF

    !
    !-----------------------------------------------------------------------
    !! 1. Set initial values for the iteration.
    !-----------------------------------------------------------------------
    !

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( FLANDG(I) <  1.0_r8 ) THEN
             IF (DB(I)  <   0.0_r8 .AND. VSHR(I)  <   2.0_r8) THEN
                !-----------------------------------------------------------------------
                !         Start the iteration from the convective limit.
                !-----------------------------------------------------------------------
                RECIP_L_MO(I) = -VKMAN/(BETA*BETA*BETA*ZH(I))
             ELSE
                !-----------------------------------------------------------------------
                !         Start the iteration from neutral values.
                !-----------------------------------------------------------------------
                RECIP_L_MO(I) = 0.0_r8
             ENDIF
          ENDIF  ! SEA_MASK
       END IF
    ENDDO

    ! DEPENDS ON: phi_m_h_sea
    CALL PHI_M_H_SEA (nCols,mskant,FLANDG,                         &
         &                  RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,                 &
         &                  PHI_M,PHI_H)
    !
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( FLANDG(I) <  1.0_r8 ) THEN
             IF (DB(I)  <   0.0_r8 .AND. VSHR(I)  <   2.0_r8) THEN
                !-----------------------------------------------------------------------
                !         Start the iteration from the convective limit.
                !-----------------------------------------------------------------------
                V_S(I) = BETA *                                           &
                     &          SQRT( BETA * ( VKMAN / PHI_H(I) ) *                   &
                     &                 ZH(I) * (-DB(I)) )
             ELSE
                !-----------------------------------------------------------------------
                !         Start the iteration from neutral values.
                !-----------------------------------------------------------------------
                V_S(I) = ( VKMAN / PHI_M(I) ) * VSHR(I)
             ENDIF
             CHV(I) = ( VKMAN / PHI_H(I) ) * V_S(I)
             CDV(I) = ( VKMAN / PHI_M(I) ) * V_S(I)
          ENDIF  ! SEA_MASK
       END IF
    ENDDO
    !-----------------------------------------------------------------------
    !! 2. Iterate to obtain sucessively better approximations for CD & CH.
    !-----------------------------------------------------------------------
    DO IT = 1,N_ITS
       !
       DO I=1,nCols
          IF(mskant(i) == 1_i8)THEN

             IF ( FLANDG(I) <  1.0_r8 ) THEN
                B_FLUX = -CHV(I) * DB(I)
                U_S = SQRT( CDV(I) * VSHR(I) )
                IF (DB(I)  <   0.0_r8) THEN
                   W_S = (ZH(I) * B_FLUX)**THIRD
                   V_S(I) = SQRT(U_S*U_S + BETA*BETA*W_S*W_S)
                ELSE
                   V_S(I) = U_S
                ENDIF
                RECIP_L_MO(I) = -VKMAN * B_FLUX /                         &
                     &                       (V_S(I)*V_S(I)*V_S(I))
             ENDIF  ! SEA_MASK
          END IF
       ENDDO
       ! DEPENDS ON: phi_m_h_sea
       CALL PHI_M_H_SEA (nCols,mskant,FLANDG,                       &
            &                    RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,               &
            &                    PHI_M,PHI_H)
       !

       DO I=1,nCols
          IF(mskant(i) == 1_i8)THEN

             IF ( FLANDG(I) <  1.0_r8 ) THEN
                CHV(I) = ( VKMAN / PHI_H(I) ) * V_S(I)
                CDV(I) = ( VKMAN / PHI_M(I) ) * V_S(I)
             ENDIF  ! SEA_MASK
          END IF
       ENDDO
    ENDDO ! Iteration loop

    !-----------------------------------------------------------------------
    !! Set CD's and CH's to be dimensionless paremters
    !-----------------------------------------------------------------------
    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( FLANDG(I) <  1.0_r8 ) THEN
             CDV(I) = (CDV(I) / VSHR(I))*CTRL_TUNABLE
             CHV(I) = (CHV(I) / VSHR(I))*CTRL_TUNABLE
          ENDIF  ! SEA_MASK
       END IF
    ENDDO


    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('FCDCH   ',4)
    !      ENDIF

    RETURN
  END SUBROUTINE FCDCH_SEA


  !!  Arguments:---------------------------------------------------------



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
!!!   SUBROUTINES PHI_M_H_SEA AND PHI_M_H_LAND ------------------------
!!!
!!!  Purpose: Calculate the integrated froms of the Monin-Obukhov
!!!           stability functions for surface exchanges.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2  15/11/00   New Deck         M. Best
  !    5.3  25/04/01  Add coastal tiling. Nic Gedney
  !  6.0  19/08/03  NEC SX-6 optimisation - force vectorisation of loop
  !                 in phi_m_h_land.  R Barnes & J-C Rioual.
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
  !*L  Arguments:---------------------------------------------------------
  SUBROUTINE PHI_M_H_SEA(                                           &
       & nCols,mskant,FLANDG,                                          &
       & RECIP_L_MO,Z_UV,Z_TQ,Z0M,Z0H,PHI_M,PHI_H)
    IMPLICIT NONE

    INTEGER                                                           &
         & nCols                            ! IN Number of X points?

    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(nCols)

    REAL(KIND=r8)                                                              &
         & FLANDG(nCols)                                          &
         !                    ! IN Land fraction
         &,RECIP_L_MO(nCols)                                      &
         !                    ! IN Reciprocal of the Monin-Obukhov length (m^-1).
         &,Z_UV(nCols)                                            &
         !                    ! IN Height of wind level above roughness height (m
         &,Z_TQ(nCols)                                            &
         !                    ! IN Height of temperature, moisture and scalar lev
         !                    !    above the roughness height (m).
         &,Z0M(nCols)                                             &
                                ! IN Roughness length for momentum (m).
         &,Z0H(nCols)
    !                    ! IN Roughness length for heat/moisture/scalars (m)
    !
    REAL(KIND=r8)                                                              &
         & PHI_M(nCols)                                           &
         !                    ! OUT Stability function for momentum.
         &,PHI_H(nCols)
    !                    ! OUT Stability function for heat/moisture/scalars.
    !
    !*L  Workspace usage----------------------------------------------------
    !    No work areas are required.
    !
    !*----------------------------------------------------------------------
    !*L  External subprograms called:

    !      EXTERNAL TIMER

    !*----------------------------------------------------------------------
    !  Common and local physical constants.
    !
    REAL(KIND=r8) A,B,C,D,C_OVER_D
    PARAMETER (                                                       &
         & A=1.0_r8                                                            &
                                !
         &,B=2.0_r8/3.0_r8                                                        &
                                ! Constants used in the Beljaars and
         &,C=5.0_r8                                                            &
                                ! Holtslag stable stability functions
         &,D=0.35_r8                                                           &
                                !
         &,C_OVER_D=C/D                                                     &
                                !
         &)
    !
    !  Define local variables.
    !
    INTEGER I     ! Loop counter; horizontal field index.
    !
    REAL(KIND=r8)                                                              &
         & PHI_MN                                                           &
                                ! Neutral value of stability function for momentum
         &,PHI_HN                                                           &
                                ! Neutral value of stability function for scalars.
         &,ZETA_UV                                                          &
                                ! Temporary in calculation of PHI_M.
         &,ZETA_0M                                                          &
                                ! Temporary in calculation of PHI_M.
         &,ZETA_TQ                                                          &
                                ! Temporary in calculation of PHI_H.
         &,ZETA_0H                                                          &
                                ! Temporary in calculation of PHI_H.
         &,X_UV_SQ                                                          &
                                ! Temporary in calculation of PHI_M.
         &,X_0M_SQ                                                          &
                                ! Temporary in calculation of PHI_M.
         &,X_UV                                                             &
                                ! Temporary in calculation of PHI_M.
         &,X_0M                                                             &
                                ! Temporary in calculation of PHI_M.
         &,Y_TQ                                                             &
                                ! Temporary in calculation of PHI_H.
         &,Y_0H                                                             &
                                ! Temporary in calculation of PHI_H.
         &,PHI_H_FZ1                                                        &
                                ! Temporary in calculation of PHI_H.
         &,PHI_H_FZ0      ! Temporary in calculation of PHI_H.

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('PHI_M_H ',3)
    !      ENDIF

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( FLANDG(I) <  1.0_r8 ) THEN

             !-----------------------------------------------------------------------
             !! 1. Calculate neutral values of PHI_M and PHI_H.
             !-----------------------------------------------------------------------
             !
             PHI_MN = LOG( (Z_UV(I) + Z0M(I)) / Z0M(I) )
             PHI_HN = LOG( (Z_TQ(I) + Z0M(I)) / Z0H(I) )
             !
             !-----------------------------------------------------------------------
             !! 2. Calculate stability parameters.
             !-----------------------------------------------------------------------
             !
             ZETA_UV = (Z_UV(I) + Z0M(I)) * RECIP_L_MO(I)
             ZETA_TQ = (Z_TQ(I) + Z0M(I)) * RECIP_L_MO(I)
             ZETA_0M = Z0M(I) * RECIP_L_MO(I)
             ZETA_0H = Z0H(I) * RECIP_L_MO(I)
             !
             !-----------------------------------------------------------------------
             !! 3. Calculate PHI_M and PHI_H for neutral and stable conditions.
             !!    Formulation of Beljaars and Holtslag (1991).
             !-----------------------------------------------------------------------
             !
             IF (RECIP_L_MO(I)  >=  0.0_r8) THEN
                PHI_M(I) = PHI_MN                                         &
                     &                 + A * (ZETA_UV - ZETA_0M)                        &
                     &                 + B * ( (ZETA_UV - C_OVER_D) * EXP(-D*ZETA_UV)   &
                     &                        -(ZETA_0M - C_OVER_D) * EXP(-D*ZETA_0M) )
                PHI_H_FZ1 = SQRT(1.0_r8 + (2.0_r8/3.0_r8)*A*ZETA_TQ)
                PHI_H_FZ0 = SQRT(1.0_r8 + (2.0_r8/3.0_r8)*A*ZETA_0H)
                PHI_H(I) = PHI_HN +                                       &
                     &                   PHI_H_FZ1*PHI_H_FZ1*PHI_H_FZ1                  &
                     &                 - PHI_H_FZ0*PHI_H_FZ0*PHI_H_FZ0                  &
                     &                 + B * ( (ZETA_TQ - C_OVER_D) * EXP(-D*ZETA_TQ)   &
                     &                        -(ZETA_0H - C_OVER_D) * EXP(-D*ZETA_0H) )
                !
                !-----------------------------------------------------------------------
                !! 4. Calculate PHI_M and PHI_H for unstable conditions.
                !-----------------------------------------------------------------------
                !
             ELSE

                X_UV_SQ = SQRT(1.0_r8 - 16.0_r8*ZETA_UV)
                X_0M_SQ = SQRT(1.0_r8 - 16.0_r8*ZETA_0M)
                X_UV = SQRT(X_UV_SQ)
                X_0M = SQRT(X_0M_SQ)
                PHI_M(I) = PHI_MN - 2.0_r8*LOG( (1.0_r8+X_UV) / (1.0_r8+X_0M) )    &
                     &                      - LOG( (1.0_r8+X_UV_SQ) / (1.0_r8+X_0M_SQ) )      &
                     &                      + 2.0_r8*( ATAN(X_UV) - ATAN(X_0M) )

                Y_TQ = SQRT(1.0_r8 - 16.0_r8*ZETA_TQ)
                Y_0H = SQRT(1.0_r8 - 16.0_r8*ZETA_0H)
                PHI_H(I) = PHI_HN - 2.0_r8*LOG( (1.0_r8+Y_TQ) / (1.0_r8+Y_0H) )

             ENDIF

          ENDIF  ! SEA_MASK
       END IF
    ENDDO

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('PHI_M_H ',4)
    !      ENDIF

    RETURN
  END SUBROUTINE PHI_M_H_SEA

!!!
  !*L  Arguments:---------------------------------------------------------




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
!!!  SUBROUTINES SF_RIB_LAND and SF_RIB_SEA ---------------------------
!!!
!!!  Purpose: Calculate bulk Richardson number for surface layer
!!!
!!!
!!!    Model            Modification history
!!!   version  date
!!!    5.2   15/11/00   New Deck         M. Best
  !    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!    Programming standard:
!!!
!!!  ------------------------------------------------------------------

  !    SUBROUTINE SF_RIB_LAND--------------------------------------------
  !
  !    Calculate RIB for land tiles
  !
  !    ------------------------------------------------------------------

  !    SUBROUTINE SF_RIB_SEA---------------------------------------------
  !
  !    Calculate RIB for sea, sea-ice and sea-ice leads
  !
  !    ------------------------------------------------------------------
  SUBROUTINE SF_RIB_SEA (                                           &
       nCols, &!INTEGER, INTENT(IN   ) ::  nCols ! IN Number of X points?
       mskant,&
       FLANDG    , &!REAL(KIND=r8)   , INTENT(IN   ) ::  FLANDG(nCols,ROWS) ! IN Land fraction on all pts.
       NSICE     , &!INTEGER, INTENT(IN   ) ::  NSICE      ! IN Number of sea-ice points.
       SICE_INDEX, &!INTEGER, INTENT(IN   ) ::  SICE_INDEX(nCols*ROWS,2) ! IN Index of sea-ice points.
       BQ_1      , &!REAL(KIND=r8)   , INTENT(IN   ) ::  BQ_1(nCols,ROWS) ! IN A buoyancy parameter for lowest atm level. ("beta-q twiddle").
       BT_1      , &!REAL(KIND=r8)   , INTENT(IN   ) ::  BT_1(nCols,ROWS)! IN A buoyancy parameter for lowest atm level. ("beta-T twiddle").
       QSTAR_ICE , &!REAL(KIND=r8)   , INTENT(IN   ) ::  QSTAR_ICE(nCols,ROWS)! IN Surface saturated sp humidity over sea-ice.
       QSTAR_SEA , &!REAL(KIND=r8)   , INTENT(IN   ) ::  QSTAR_SEA(nCols,ROWS)! IN Surface saturated sp humidity over sea and sea-ice leads.
       QW_1      , &!REAL(KIND=r8)   , INTENT(IN   ) ::  QW_1(nCols,ROWS)! IN Total water content of lowest atmospheric layer (kg per kg air)
       TL_1      , &!REAL(KIND=r8)   , INTENT(IN   ) ::  TL_1(nCols,ROWS)! IN Liquid/frozen water temperature for lowest atmospheric layer (K).
       TSTAR_SICE, &!REAL(KIND=r8)   , INTENT(IN   ) ::  TSTAR_SICE(nCols,ROWS)! IN Surface temperature of sea-ice (K).
       TSTAR_SEA , &!REAL(KIND=r8)   , INTENT(IN   ) ::  TSTAR_SEA(nCols,ROWS)! IN Surface temperature of sea and sea-ice leads (K).
       VSHR      , &!REAL(KIND=r8)   , INTENT(IN   ) ::  VSHR(nCols,ROWS)! IN Magnitude of surface- to-lowest-level wind shear.
       Z0H_ICE   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0H_ICE(nCols,ROWS)! IN Roughness length for heat and moisture transport over sea-ice (m).
       Z0H_SEA   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0H_SEA(nCols,ROWS)! IN Roughness length for heat and moisture transport over sea or sea-ice leads (m).
       Z0M_ICE   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0M_ICE(nCols,ROWS)! IN Roughness length for momentum over sea-ice (m).
       Z0M_SEA   , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z0M_SEA(nCols,ROWS)! IN Roughness length for momentum over sea or sea-ice leads (m).
       Z1_TQ     , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z1_TQ(nCols,ROWS)! IN Height of lowest TQ level (m).
       Z1_UV     , &!REAL(KIND=r8)   , INTENT(IN   ) ::  Z1_UV(nCols,ROWS)          ! IN Height of lowest UV level (m).
       RIB_SEA   , &!REAL(KIND=r8)   , INTENT(OUT  ) ::  RIB_SEA(nCols,ROWS)! OUT Bulk Richardson number for lowest layer over sea or sea-ice leads.
       RIB_ICE   , &!REAL(KIND=r8)   , INTENT(OUT  ) ::  RIB_ICE(nCols,ROWS)! OUT Bulk Richardson number for lowest layer over sea-ice.
       DB_SEA    , &!REAL(KIND=r8)   , INTENT(OUT  ) ::  DB_SEA(nCols,ROWS)! OUT Buoyancy difference between
       DB_ICE     &!REAL(KIND=r8)   , INTENT(OUT  ) ::  DB_ICE(nCols,ROWS)   ! OUT Buoyancy difference between
       )

    IMPLICIT NONE

    INTEGER, INTENT(IN   ) ::  nCols ! IN Number of X points?
    INTEGER, INTENT(IN   ) ::  NSICE      ! IN Number of sea-ice points.
    INTEGER, INTENT(IN   ) ::  SICE_INDEX(nCols,2) ! IN Index of sea-ice points.
    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(nCols)

    REAL(KIND=r8)   , INTENT(IN   ) ::  FLANDG(nCols) ! IN Land fraction on all pts.
    REAL(KIND=r8)   , INTENT(IN   ) ::  BQ_1(nCols) ! IN A buoyancy parameter for lowest atm level. ("beta-q twiddle").
    REAL(KIND=r8)   , INTENT(IN   ) ::  BT_1(nCols)! IN A buoyancy parameter for lowest atm level. ("beta-T twiddle").
    REAL(KIND=r8)   , INTENT(IN   ) ::  QSTAR_ICE(nCols)! IN Surface saturated sp humidity over sea-ice.
    REAL(KIND=r8)   , INTENT(IN   ) ::  QSTAR_SEA(nCols)! IN Surface saturated sp humidity over sea and sea-ice leads.
    REAL(KIND=r8)   , INTENT(IN   ) ::  QW_1(nCols)! IN Total water content of lowest atmospheric layer (kg per kg air)
    REAL(KIND=r8)   , INTENT(IN   ) ::  TL_1(nCols)! IN Liquid/frozen water temperature for lowest atmospheric layer (K).
    REAL(KIND=r8)   , INTENT(IN   ) ::  TSTAR_SICE(nCols)! IN Surface temperature of sea-ice (K).
    REAL(KIND=r8)   , INTENT(IN   ) ::  TSTAR_SEA(nCols)! IN Surface temperature of sea and sea-ice leads (K).
    REAL(KIND=r8)   , INTENT(IN   ) ::  VSHR(nCols)! IN Magnitude of surface- to-lowest-level wind shear.
    REAL(KIND=r8)   , INTENT(IN   ) ::  Z0H_ICE(nCols)! IN Roughness length for heat and moisture transport over sea-ice (m).
    REAL(KIND=r8)   , INTENT(IN   ) ::  Z0H_SEA(nCols)! IN Roughness length for heat and moisture transport over sea or sea-ice leads (m).
    REAL(KIND=r8)   , INTENT(IN   ) ::  Z0M_ICE(nCols)! IN Roughness length for momentum over sea-ice (m).
    REAL(KIND=r8)   , INTENT(IN   ) ::  Z0M_SEA(nCols)! IN Roughness length for momentum over sea or sea-ice leads (m).
    REAL(KIND=r8)   , INTENT(IN   ) ::  Z1_TQ(nCols)! IN Height of lowest TQ level (m).
    REAL(KIND=r8)   , INTENT(IN   ) ::  Z1_UV(nCols)     ! IN Height of lowest UV level (m).
    REAL(KIND=r8)   , INTENT(OUT  ) ::  RIB_SEA(nCols)! OUT Bulk Richardson number for lowest layer over sea or sea-ice leads.
    REAL(KIND=r8)   , INTENT(OUT  ) ::  RIB_ICE(nCols)! OUT Bulk Richardson number for lowest layer over sea-ice.
    REAL(KIND=r8)   , INTENT(OUT  ) ::  DB_SEA(nCols)! OUT Buoyancy difference between
    !                                !     surface and lowest atmospheric
    !                                !     level over sea or sea-ice leads.
    REAL(KIND=r8)   , INTENT(OUT  ) ::  DB_ICE(nCols)   ! OUT Buoyancy difference between
    !                                !     surface and lowest atmospheric
    !                                !     level over sea-ice.


    !  External routines called :-
    !      EXTERNAL TIMER


    !  Symbolic constants -----------------------------------------------

    !*L------------------COMDECK C_O_DG_C-----------------------------------
    ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
    ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
    ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

    REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
    REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
    REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_R_CP-------------------------------------
    ! History:
    ! Version  Date      Comment.
    !  5.0  07/05/99  Add variable P_zero for consistency with
    !                 conversion to C-P 'C' dynamics grid. R. Rawlins
    !  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

    ! R IS GAS CONSTANT FOR DRY AIR
    ! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
    ! PREF IS REFERENCE SURFACE PRESSURE

    REAL(KIND=r8), PARAMETER  :: R      = 287.05_r8
    REAL(KIND=r8), PARAMETER  :: CP     = 1005.0_r8
    REAL(KIND=r8), PARAMETER  :: Kappa  = R/CP
    REAL(KIND=r8), PARAMETER  :: Pref   = 100000.0_r8

    ! Reference surface pressure = PREF
    REAL(KIND=r8), PARAMETER  :: P_zero = Pref
    REAL(KIND=r8), PARAMETER  :: sclht  = 6.8E+03_r8  ! Scale Height H
    !*----------------------------------------------------------------------

    !  Workspace --------------------------------------------------------
    INTEGER                                                           &
         & I                                                              &
                                ! Horizontal field index.
         &,SI                  !Sea-ice field index.
    REAL(KIND=r8)                                                              &
         & DQ                                                               &
                                ! Sp humidity difference between surface
         !                          ! and lowest atmospheric level (Q1 - Q*).
         &,DTEMP               ! Modified temperature difference between
    !                          ! surface and lowest atmospheric level.

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !        CALL TIMER('SF_RIB  ',3)
    !      ENDIF

    DO I=1,nCols
       IF(mskant(i) == 1_i8)THEN

          IF ( FLANDG(I) <  1.0_r8 ) THEN
             ! Sea and sea-ice leads
             DTEMP = TL_1(I) - TSTAR_SEA(I)                            &
                                ! P243.118
                  &            + (G/CP)*(Z1_TQ(I) + Z0M_SEA(I) - Z0H_SEA(I))
             DQ = QW_1(I) - QSTAR_SEA(I)                  ! P243.119
             DB_SEA(I) = G*( BT_1(I)*DTEMP + BQ_1(I)*DQ )
             RIB_SEA(I) = Z1_UV(I)*DB_SEA(I) /                       &
                  &                     ( VSHR(I)*VSHR(I) )
          ENDIF
       END IF
    ENDDO

    DO SI=1,NSICE
       IF(mskant(i) == 1_i8)THEN

          I = SICE_INDEX(SI,1)
          ! Sea-ice
          DTEMP = TL_1(I) - TSTAR_SICE(I)                             &
               &           + (G/CP)*(Z1_TQ(I) + Z0M_ICE(I) - Z0H_ICE(I))
          DQ = QW_1(I) - QSTAR_ICE(I)
          DB_ICE(I) = G*( BT_1(I)*DTEMP + BQ_1(I)*DQ )
          RIB_ICE(I) = Z1_UV(I)*DB_ICE(I) /                         &
               &                    ( VSHR(I) * VSHR(I) )
       END IF
    ENDDO

    !      IF (LTIMER) THEN
    ! DEPENDS ON: timer
    !         !CALL TIMER('SF_RIB  ',4)
    !      ENDIF

    RETURN
  END SUBROUTINE SF_RIB_SEA

  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  ! SUBROUTINE BOUY_TQ

  ! PURPOSE: To calculate buoyancy parameters on p,T,q-levels
  !
  ! METHOD:
  !
  ! HISTORY:
  ! DATE   VERSION   COMMENT
  ! ----   -------   -------
  ! new deck
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
  !LL   5.2  27/09/00   change from QSAT_2D to QSAT           A.Malcolm
!!!  5.3    Feb. 01  Calculate grid-box mean parameters (APLock)
  !  6.1  17/05/04  Change Q, QCL, QCF dims to enable substepping.
  !                                                       M. Diamantakis
  !
  ! CODE DESCRIPTION:
  !   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
  !   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
  !

  SUBROUTINE BOUY_TQ (&
       nCols, & !INTEGER, INTENT(IN   ) :: nCols
       kMAx , & !INTEGER, INTENT(IN   ) :: kMAx   
       LQ_MIX_BL , & !LOGICAL, INTENT(IN   ) :: LQ_MIX_BL   
       mskant    , &
       P         , & !REAL(KIND=r8), INTENT(IN   ) :: P(nCols,rows,kMAx)
       T         , & !REAL(KIND=r8), INTENT(IN   ) ::T(nCols,rows,kMAx)
       Q         , & !REAL(KIND=r8), INTENT(IN   ) ::Q(nCols,rows,kMAx)
       QCF       , & !REAL(KIND=r8), INTENT(IN   ) ::QCF(nCols,rows,kMAx)
       QCL       , & !REAL(KIND=r8), INTENT(IN   ) ::QCL(nCols,rows,kMAx)
       CF        , & !REAL(KIND=r8), INTENT(IN   ) ::CF(nCols, rows, kMAx)
       BT        , & !REAL(KIND=r8), INTENT(out   ) :: BT(nCols,rows,kMAx)
       BQ        , & !REAL(KIND=r8), INTENT(out   ) :: BQ(nCols,rows,kMAx)
       BT_CLD    , & !REAL(KIND=r8), INTENT(out   ) :: BT_CLD(nCols,rows,kMAx)
       BQ_CLD    , & !REAL(KIND=r8), INTENT(out   ) :: BQ_CLD(nCols,rows,kMAx)
       BT_GB     , & !REAL(KIND=r8), INTENT(out   ) :: BT_GB(nCols,rows,kMAx)
       BQ_GB     , & !REAL(KIND=r8), INTENT(out   ) :: BQ_GB(nCols,rows,kMAx)
       A_QS      , & !REAL(KIND=r8), INTENT(out   ) :: A_QS(nCols,rows,kMAx)
       A_DQSDT   , & !REAL(KIND=r8), INTENT(out   ) :: A_DQSDT(nCols,rows,kMAx)
       DQSDT      & !REAL(KIND=r8), INTENT(out   ) :: DQSDT(nCols,rows,kMAx)
       )

    IMPLICIT NONE

    ! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.



    LOGICAL, INTENT(IN   ) :: LQ_MIX_BL       ! IN switch for using mixing ratios
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMAx              ! IN No. of atmospheric levels for which
    !                                boundary layer fluxes are calculated.
    !                                Assumed  <=30 for dimensioning GAMMA()
    !                                in common deck C_GAMMA
    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(nCols) 

    REAL(KIND=r8), INTENT(IN   ) :: P(nCols,kMAx)
    ! IN Pressure at pressure points.
    REAL(KIND=r8), INTENT(IN   ) ::T(nCols,kMAx)
    ! IN Temperature (K). At P points
    REAL(KIND=r8), INTENT(IN   ) ::Q(nCols,kMAx)
    ! IN Sp humidity (kg water per kg air).
    REAL(KIND=r8), INTENT(IN   ) ::QCL(nCols,kMAx)
    ! IN Cloud liq water (kg per kg air).
    REAL(KIND=r8), INTENT(IN   ) ::QCF(nCols,kMAx)
    ! IN Cloud liq water (kg per kg air).
    REAL(KIND=r8), INTENT(IN   ) ::CF(nCols, kMAx)! IN Cloud fraction (decimal).


    ! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

    REAL(KIND=r8), INTENT(out   ) :: BQ(nCols,kMAx)
    ! OUT A buoyancy parameter for clear air
    REAL(KIND=r8), INTENT(out   ) :: BT(nCols,kMAx)
    ! OUT A buoyancy parameter for clear air
    REAL(KIND=r8), INTENT(out   ) :: BQ_CLD(nCols,kMAx)
    !                             ! OUT A buoyancy parameter for cloudy air
    REAL(KIND=r8), INTENT(out   ) :: BT_CLD(nCols,kMAx)
    !                             ! OUT A buoyancy parameter for cloudy air
    REAL(KIND=r8), INTENT(out   ) :: BQ_GB(nCols,kMAx)
    ! OUT A grid-box mean buoyancy parameter
    REAL(KIND=r8), INTENT(out   ) :: BT_GB(nCols,kMAx)
    ! OUT A grid-box mean buoyancy parameter
    REAL(KIND=r8), INTENT(out   ) :: A_QS(nCols,kMAx)
    !                             ! OUT Saturated lapse rate factor
    REAL(KIND=r8), INTENT(out   ) :: A_DQSDT(nCols,kMAx)
    !                             ! OUT Saturated lapse rate factor
    REAL(KIND=r8), INTENT(out   ) :: DQSDT(nCols,kMAx)
    !                             ! OUT Derivative of q_SAT w.r.t. T

    ! LOCAL VARIABLES.
    REAL(KIND=r8) ::  qss(nCols)
    REAL(KIND=r8)                                                              &
         & QS(nCols)            ! WORK Saturated mixing ratio.

    INTEGER                                                           &
         &  I                                                             &
         &, K

    REAL(KIND=r8)                                                              &
         &  BC

    !     EXTERNAL                                                          &
    !    &  QSAT, TIMER

    !*L------------------COMDECK C_O_DG_C-----------------------------------
    ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
    ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
    ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

    REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
    REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
    REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

    !*----------------------------------------------------------------------
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8

    ! C_LHEAT end
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_R_CP-------------------------------------
    ! History:
    ! Version  Date      Comment.
    !  5.0  07/05/99  Add variable P_zero for consistency with
    !                 conversion to C-P 'C' dynamics grid. R. Rawlins
    !  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

    ! R IS GAS CONSTANT FOR DRY AIR
    ! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
    ! PREF IS REFERENCE SURFACE PRESSURE

    REAL(KIND=r8), PARAMETER  :: R      = 287.05_r8
    REAL(KIND=r8), PARAMETER  :: CP     = 1005.0_r8
    REAL(KIND=r8), PARAMETER  :: Kappa  = R/CP
    REAL(KIND=r8), PARAMETER  :: Pref   = 100000.0_r8

    ! Reference surface pressure = PREF
    REAL(KIND=r8), PARAMETER  :: P_zero = Pref
    REAL(KIND=r8), PARAMETER  :: sclht  = 6.8E+03_r8  ! Scale Height H
    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_EPSLON-----------------------------------
    ! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

    REAL(KIND=r8), PARAMETER :: Epsilon   = 0.62198_r8
    REAL(KIND=r8), PARAMETER :: C_Virtual = 1.0_r8/Epsilon-1.0_r8

    !*----------------------------------------------------------------------
    ! C_VKMAN start
    REAL(KIND=r8),PARAMETER:: VKMAN=0.4_r8 ! Von Karman's constant
    ! C_VKMAN end
    ! C_SOILH start
    ! No. of soil layers (must = NSOIL).
    REAL(KIND=r8),PARAMETER:: PSOIL=4.0_r8

    ! Tunable characteristic freq (rad/s)
    REAL(KIND=r8),PARAMETER:: OMEGA1=3.55088E-4_r8

    ! Density of lying snow (kg per m**3)
    REAL(KIND=r8),PARAMETER:: RHO_SNOW=250.0_r8

    ! Depth of `effective' snow surface layer (m)
    REAL(KIND=r8),PARAMETER:: DEFF_SNOW=0.1_r8

    ! Thermal conductivity of lying snow (Watts per m per K).
    REAL(KIND=r8),PARAMETER:: SNOW_HCON=0.265_r8

    ! Thermal capacity of lying snow (J/K/m3)
    REAL(KIND=r8),PARAMETER:: SNOW_HCAP=0.63E6_r8

    ! C_SOILH end
    ! History:
    ! Version  Date  Comment
    !  3.4   18/5/94 Add PP missing data indicator. J F Thomson
    !  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
    !*L------------------COMDECK C_MDI-------------------------------------
    ! PP missing data indicator (-1.0E+30)
    REAL(KIND=r8), PARAMETER    :: RMDI_PP  = -1.0E+30_r8

    ! Old REAL(KIND=r8) missing data indicator (-32768.0)
    REAL(KIND=r8), PARAMETER    :: RMDI_OLD = -32768.0_r8

    ! New REAL(KIND=r8) missing data indicator (-2**30)
    REAL(KIND=r8), PARAMETER    :: RMDI     = -32768.0_r8*32768.0_r8

    ! Integer missing data indicator
    INTEGER, PARAMETER :: IMDI     = -32768
    !*----------------------------------------------------------------------


    REAL(KIND=r8) ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
    PARAMETER (                                                       &
         & ETAR=1.0_r8/(1.0_r8-EPSILON)                                           &
                                ! Used in buoyancy parameter BETAC.
         &,GRCP=G/CP                                                        &
                                ! Used in DZTL, FTL calculations.
         &,LCRCP=LC/CP                                                      &
                                ! Latent heat of condensation / CP.
         &,LFRCP=LF/CP                                                      &
                                ! Latent heat of fusion / CP.
         &,LS=LC+LF                                                         &
                                ! Latent heat of sublimation.
         &,LSRCP=LS/CP                                                      &
                                ! Latent heat of sublimation / CP.
         &)

    !      IF (LTIMER) THEN
    !! DEPENDS ON: timer
    !        CALL TIMER('BOUY_TQ ',3)
    !      ENDIF
    !-----------------------------------------------------------------------
    !! 1.  Loop round levels.
    !-----------------------------------------------------------------------
    DO K=1,kMAx
       !-----------------------------------------------------------------------
       !! 1.1 Calculate saturated specific humidity at pressure and
       !!     temperature of current level.
       !-----------------------------------------------------------------------
       ! DEPENDS ON: qsat_mix
       IF(QSAT_opt)THEN
          CALL QSAT_mix(QS,T(1,K),P(1,K),nCols,Lq_mix_bl,mskant)
       END IF
       !
       DO I=1,nCols
          IF(mskant(i) == 1_i8)THEN
             IF(.NOT. QSAT_opt)THEN
                qss(i) = fpvs(T(i,k))
                QS (i) = EPS * qss(i) / (P(i,k) + EPSM1 * qss(i))
             END IF

             !ajm        DO I=P1,P1+P_POINTS-1

             !-----------------------------------------------------------------------
             !! 1.2 Calculate buoyancy parameters BT and BQ, required for the
             !!     calculation of stability.
             !-----------------------------------------------------------------------

             BT(I,K) = 1.0_r8/T(I,K)
             BQ(I,K) =                                                   &
                  &      C_VIRTUAL/(1.0_r8+C_VIRTUAL*Q(I,K)-QCL(I,K)-QCF(I,K))
             !
             IF (T(I,K)  >   TM) THEN
                DQSDT(I,K) = (EPSILON * LC * QS(I))                     &
                     &                   / ( R * T(I,K) * T(I,K) )
                !                      ...  (Clausius-Clapeyron) for T above freezing
                !
                A_QS(I,K) = 1.0_r8 / (1.0_r8 + LCRCP*DQSDT(I,K))
                !
                A_DQSDT(I,K) = A_QS(I,K) * DQSDT(I,K)
                !
                BC = LCRCP*BT(I,K) - ETAR*BQ(I,K)
                !
             ELSE
                DQSDT(I,K) = (EPSILON * LS * QS(I))                     &
                     &                   / ( R * T(I,K) * T(I,K) )
                !                      ...  (Clausius-Clapeyron) for T below freezing
                !
                A_QS(I,K) = 1.0_r8 / (1.0_r8 + LSRCP*DQSDT(I,K))
                !
                A_DQSDT(I,K) = A_QS(I,K) * DQSDT(I,K)
                !
                BC = LSRCP*BT(I,K) - ETAR*BQ(I,K)
                !
             ENDIF
             !
             !-----------------------------------------------------------------------
             !! 1.3 Calculate in-cloud buoyancy parameters.
             !-----------------------------------------------------------------------
             !
             BT_CLD(I,K) = BT(I,K) - A_DQSDT(I,K) * BC
             BQ_CLD(I,K) = BQ(I,K) + A_QS(I,K) * BC

             !-----------------------------------------------------------------------
             !! 1.4 Calculate grid-box mean buoyancy parameters.
             !-----------------------------------------------------------------------
             !
             BT_GB(I,K) = BT(I,K) + CF(I,K)*( BT_CLD(I,K) - BT(I,K) )
             BQ_GB(I,K) = BQ(I,K) + CF(I,K)*( BQ_CLD(I,K) - BQ(I,K) )
             !
          END IF

       ENDDO ! p_points,j
    ENDDO ! kMAx

    !      IF (LTIMER) THEN
    !! DEPENDS ON: timer
    !        CALL TIMER('BOUY_TQ ',4)
    !      ENDIF
    RETURN
  END SUBROUTINE BOUY_TQ




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  SUBROUTINE qsat_data()
    !
    !  Local dynamic arrays-------------------------------------------------
    REAL(KIND=r8) ES_Aux(0:N+1)    ! TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
    !                       - SET BY DATA STATEMENT CALCULATED FROM THE
    !                         GOFF-GRATCH FORMULAE AS TAKEN FROM LANDOLT-
    !                         BORNSTEIN, 1987 NUMERICAL DATA AND FUNCTIONAL
    !                         RELATIONSHIPS IN SCIENCE AND TECHNOLOGY.
    !                         GROUP V/ VOL 4B METEOROLOGY. PHYSICAL AND
    !                         CHEMICAL PROPERTIES OF AIR, P35
    !

    ! ----------------------------------------------------------------------
    ! SATURATION WATER VAPOUR PRESSURE
    !
    ! ABOVE 0 DEG C VALUES ARE OVER WATER
    !
    ! BELOW 0 DEC C VALUES ARE OVER ICE
    ! ----------------------------------------------------------------------
    ! Note: 0 element is a repeat of 1st element to cater for special case
    !       of low temperatures ( <= T_LOW) for which the array index is
    !       rounded down due to machine precision.
    DATA (ES_Aux(IES),IES=    0, 95) / 0.966483E-02_r8,                      &
         0.966483E-02_r8,0.984279E-02_r8,0.100240E-01_r8,0.102082E-01_r8,0.103957E-01_r8, &
         0.105865E-01_r8,0.107803E-01_r8,0.109777E-01_r8,0.111784E-01_r8,0.113825E-01_r8, &
         0.115902E-01_r8,0.118016E-01_r8,0.120164E-01_r8,0.122348E-01_r8,0.124572E-01_r8, &
         0.126831E-01_r8,0.129132E-01_r8,0.131470E-01_r8,0.133846E-01_r8,0.136264E-01_r8, &
         0.138724E-01_r8,0.141225E-01_r8,0.143771E-01_r8,0.146356E-01_r8,0.148985E-01_r8, &
         0.151661E-01_r8,0.154379E-01_r8,0.157145E-01_r8,0.159958E-01_r8,0.162817E-01_r8, &
         0.165725E-01_r8,0.168680E-01_r8,0.171684E-01_r8,0.174742E-01_r8,0.177847E-01_r8, &
         0.181008E-01_r8,0.184216E-01_r8,0.187481E-01_r8,0.190801E-01_r8,0.194175E-01_r8, &
         0.197608E-01_r8,0.201094E-01_r8,0.204637E-01_r8,0.208242E-01_r8,0.211906E-01_r8, &
         0.215631E-01_r8,0.219416E-01_r8,0.223263E-01_r8,0.227172E-01_r8,0.231146E-01_r8, &
         0.235188E-01_r8,0.239296E-01_r8,0.243465E-01_r8,0.247708E-01_r8,0.252019E-01_r8, &
         0.256405E-01_r8,0.260857E-01_r8,0.265385E-01_r8,0.269979E-01_r8,0.274656E-01_r8, &
         0.279405E-01_r8,0.284232E-01_r8,0.289142E-01_r8,0.294124E-01_r8,0.299192E-01_r8, &
         0.304341E-01_r8,0.309571E-01_r8,0.314886E-01_r8,0.320285E-01_r8,0.325769E-01_r8, &
         0.331348E-01_r8,0.337014E-01_r8,0.342771E-01_r8,0.348618E-01_r8,0.354557E-01_r8, &
         0.360598E-01_r8,0.366727E-01_r8,0.372958E-01_r8,0.379289E-01_r8,0.385717E-01_r8, &
         0.392248E-01_r8,0.398889E-01_r8,0.405633E-01_r8,0.412474E-01_r8,0.419430E-01_r8, &
         0.426505E-01_r8,0.433678E-01_r8,0.440974E-01_r8,0.448374E-01_r8,0.455896E-01_r8, &
         0.463545E-01_r8,0.471303E-01_r8,0.479191E-01_r8,0.487190E-01_r8,0.495322E-01_r8/
    DATA (ES_Aux(IES),IES= 96,190) /                                      &
         0.503591E-01_r8,0.511977E-01_r8,0.520490E-01_r8,0.529145E-01_r8,0.537931E-01_r8, &
         0.546854E-01_r8,0.555924E-01_r8,0.565119E-01_r8,0.574467E-01_r8,0.583959E-01_r8, &
         0.593592E-01_r8,0.603387E-01_r8,0.613316E-01_r8,0.623409E-01_r8,0.633655E-01_r8, &
         0.644053E-01_r8,0.654624E-01_r8,0.665358E-01_r8,0.676233E-01_r8,0.687302E-01_r8, &
         0.698524E-01_r8,0.709929E-01_r8,0.721490E-01_r8,0.733238E-01_r8,0.745180E-01_r8, &
         0.757281E-01_r8,0.769578E-01_r8,0.782061E-01_r8,0.794728E-01_r8,0.807583E-01_r8, &
         0.820647E-01_r8,0.833905E-01_r8,0.847358E-01_r8,0.861028E-01_r8,0.874882E-01_r8, &
         0.888957E-01_r8,0.903243E-01_r8,0.917736E-01_r8,0.932464E-01_r8,0.947407E-01_r8, &
         0.962571E-01_r8,0.977955E-01_r8,0.993584E-01_r8,0.100942E+00_r8,0.102551E+00_r8, &
         0.104186E+00_r8,0.105842E+00_r8,0.107524E+00_r8,0.109231E+00_r8,0.110963E+00_r8, &
         0.112722E+00_r8,0.114506E+00_r8,0.116317E+00_r8,0.118153E+00_r8,0.120019E+00_r8, &
         0.121911E+00_r8,0.123831E+00_r8,0.125778E+00_r8,0.127755E+00_r8,0.129761E+00_r8, &
         0.131796E+00_r8,0.133863E+00_r8,0.135956E+00_r8,0.138082E+00_r8,0.140241E+00_r8, &
         0.142428E+00_r8,0.144649E+00_r8,0.146902E+00_r8,0.149190E+00_r8,0.151506E+00_r8, &
         0.153859E+00_r8,0.156245E+00_r8,0.158669E+00_r8,0.161126E+00_r8,0.163618E+00_r8, &
         0.166145E+00_r8,0.168711E+00_r8,0.171313E+00_r8,0.173951E+00_r8,0.176626E+00_r8, &
         0.179342E+00_r8,0.182096E+00_r8,0.184893E+00_r8,0.187724E+00_r8,0.190600E+00_r8, &
         0.193518E+00_r8,0.196473E+00_r8,0.199474E+00_r8,0.202516E+00_r8,0.205604E+00_r8, &
         0.208730E+00_r8,0.211905E+00_r8,0.215127E+00_r8,0.218389E+00_r8,0.221701E+00_r8/
    DATA (ES_Aux(IES),IES=191,285) /                                      &
         0.225063E+00_r8,0.228466E+00_r8,0.231920E+00_r8,0.235421E+00_r8,0.238976E+00_r8, &
         0.242580E+00_r8,0.246232E+00_r8,0.249933E+00_r8,0.253691E+00_r8,0.257499E+00_r8, &
         0.261359E+00_r8,0.265278E+00_r8,0.269249E+00_r8,0.273274E+00_r8,0.277358E+00_r8, &
         0.281498E+00_r8,0.285694E+00_r8,0.289952E+00_r8,0.294268E+00_r8,0.298641E+00_r8, &
         0.303078E+00_r8,0.307577E+00_r8,0.312135E+00_r8,0.316753E+00_r8,0.321440E+00_r8, &
         0.326196E+00_r8,0.331009E+00_r8,0.335893E+00_r8,0.340842E+00_r8,0.345863E+00_r8, &
         0.350951E+00_r8,0.356106E+00_r8,0.361337E+00_r8,0.366636E+00_r8,0.372006E+00_r8, &
         0.377447E+00_r8,0.382966E+00_r8,0.388567E+00_r8,0.394233E+00_r8,0.399981E+00_r8, &
         0.405806E+00_r8,0.411714E+00_r8,0.417699E+00_r8,0.423772E+00_r8,0.429914E+00_r8, &
         0.436145E+00_r8,0.442468E+00_r8,0.448862E+00_r8,0.455359E+00_r8,0.461930E+00_r8, &
         0.468596E+00_r8,0.475348E+00_r8,0.482186E+00_r8,0.489124E+00_r8,0.496160E+00_r8, &
         0.503278E+00_r8,0.510497E+00_r8,0.517808E+00_r8,0.525224E+00_r8,0.532737E+00_r8, &
         0.540355E+00_r8,0.548059E+00_r8,0.555886E+00_r8,0.563797E+00_r8,0.571825E+00_r8, &
         0.579952E+00_r8,0.588198E+00_r8,0.596545E+00_r8,0.605000E+00_r8,0.613572E+00_r8, &
         0.622255E+00_r8,0.631059E+00_r8,0.639962E+00_r8,0.649003E+00_r8,0.658144E+00_r8, &
         0.667414E+00_r8,0.676815E+00_r8,0.686317E+00_r8,0.695956E+00_r8,0.705728E+00_r8, &
         0.715622E+00_r8,0.725641E+00_r8,0.735799E+00_r8,0.746082E+00_r8,0.756495E+00_r8, &
         0.767052E+00_r8,0.777741E+00_r8,0.788576E+00_r8,0.799549E+00_r8,0.810656E+00_r8, &
         0.821914E+00_r8,0.833314E+00_r8,0.844854E+00_r8,0.856555E+00_r8,0.868415E+00_r8/
    DATA (ES_Aux(IES),IES=286,380) /                                      &
         0.880404E+00_r8,0.892575E+00_r8,0.904877E+00_r8,0.917350E+00_r8,0.929974E+00_r8, &
         0.942771E+00_r8,0.955724E+00_r8,0.968837E+00_r8,0.982127E+00_r8,0.995600E+00_r8, &
         0.100921E+01_r8,0.102304E+01_r8,0.103700E+01_r8,0.105116E+01_r8,0.106549E+01_r8, &
         0.108002E+01_r8,0.109471E+01_r8,0.110962E+01_r8,0.112469E+01_r8,0.113995E+01_r8, &
         0.115542E+01_r8,0.117107E+01_r8,0.118693E+01_r8,0.120298E+01_r8,0.121923E+01_r8, &
         0.123569E+01_r8,0.125234E+01_r8,0.126923E+01_r8,0.128631E+01_r8,0.130362E+01_r8, &
         0.132114E+01_r8,0.133887E+01_r8,0.135683E+01_r8,0.137500E+01_r8,0.139342E+01_r8, &
         0.141205E+01_r8,0.143091E+01_r8,0.145000E+01_r8,0.146933E+01_r8,0.148892E+01_r8, &
         0.150874E+01_r8,0.152881E+01_r8,0.154912E+01_r8,0.156970E+01_r8,0.159049E+01_r8, &
         0.161159E+01_r8,0.163293E+01_r8,0.165452E+01_r8,0.167640E+01_r8,0.169852E+01_r8, &
         0.172091E+01_r8,0.174359E+01_r8,0.176653E+01_r8,0.178977E+01_r8,0.181332E+01_r8, &
         0.183709E+01_r8,0.186119E+01_r8,0.188559E+01_r8,0.191028E+01_r8,0.193524E+01_r8, &
         0.196054E+01_r8,0.198616E+01_r8,0.201208E+01_r8,0.203829E+01_r8,0.206485E+01_r8, &
         0.209170E+01_r8,0.211885E+01_r8,0.214637E+01_r8,0.217424E+01_r8,0.220242E+01_r8, &
         0.223092E+01_r8,0.225979E+01_r8,0.228899E+01_r8,0.231855E+01_r8,0.234845E+01_r8, &
         0.237874E+01_r8,0.240937E+01_r8,0.244040E+01_r8,0.247176E+01_r8,0.250349E+01_r8, &
         0.253560E+01_r8,0.256814E+01_r8,0.260099E+01_r8,0.263431E+01_r8,0.266800E+01_r8, &
         0.270207E+01_r8,0.273656E+01_r8,0.277145E+01_r8,0.280671E+01_r8,0.284248E+01_r8, &
         0.287859E+01_r8,0.291516E+01_r8,0.295219E+01_r8,0.298962E+01_r8,0.302746E+01_r8/
    DATA (ES_Aux(IES),IES=381,475) /                                      &
         0.306579E+01_r8,0.310454E+01_r8,0.314377E+01_r8,0.318351E+01_r8,0.322360E+01_r8, &
         0.326427E+01_r8,0.330538E+01_r8,0.334694E+01_r8,0.338894E+01_r8,0.343155E+01_r8, &
         0.347456E+01_r8,0.351809E+01_r8,0.356216E+01_r8,0.360673E+01_r8,0.365184E+01_r8, &
         0.369744E+01_r8,0.374352E+01_r8,0.379018E+01_r8,0.383743E+01_r8,0.388518E+01_r8, &
         0.393344E+01_r8,0.398230E+01_r8,0.403177E+01_r8,0.408175E+01_r8,0.413229E+01_r8, &
         0.418343E+01_r8,0.423514E+01_r8,0.428746E+01_r8,0.434034E+01_r8,0.439389E+01_r8, &
         0.444808E+01_r8,0.450276E+01_r8,0.455820E+01_r8,0.461423E+01_r8,0.467084E+01_r8, &
         0.472816E+01_r8,0.478607E+01_r8,0.484468E+01_r8,0.490393E+01_r8,0.496389E+01_r8, &
         0.502446E+01_r8,0.508580E+01_r8,0.514776E+01_r8,0.521047E+01_r8,0.527385E+01_r8, &
         0.533798E+01_r8,0.540279E+01_r8,0.546838E+01_r8,0.553466E+01_r8,0.560173E+01_r8, &
         0.566949E+01_r8,0.573807E+01_r8,0.580750E+01_r8,0.587749E+01_r8,0.594846E+01_r8, &
         0.602017E+01_r8,0.609260E+01_r8,0.616591E+01_r8,0.623995E+01_r8,0.631490E+01_r8, &
         0.639061E+01_r8,0.646723E+01_r8,0.654477E+01_r8,0.662293E+01_r8,0.670220E+01_r8, &
         0.678227E+01_r8,0.686313E+01_r8,0.694495E+01_r8,0.702777E+01_r8,0.711142E+01_r8, &
         0.719592E+01_r8,0.728140E+01_r8,0.736790E+01_r8,0.745527E+01_r8,0.754352E+01_r8, &
         0.763298E+01_r8,0.772316E+01_r8,0.781442E+01_r8,0.790676E+01_r8,0.800001E+01_r8, &
         0.809435E+01_r8,0.818967E+01_r8,0.828606E+01_r8,0.838343E+01_r8,0.848194E+01_r8, &
         0.858144E+01_r8,0.868207E+01_r8,0.878392E+01_r8,0.888673E+01_r8,0.899060E+01_r8, &
         0.909567E+01_r8,0.920172E+01_r8,0.930909E+01_r8,0.941765E+01_r8,0.952730E+01_r8/
    DATA (ES_Aux(IES),IES=476,570) /                                      &
         0.963821E+01_r8,0.975022E+01_r8,0.986352E+01_r8,0.997793E+01_r8,0.100937E+02_r8, &
         0.102105E+02_r8,0.103287E+02_r8,0.104481E+02_r8,0.105688E+02_r8,0.106909E+02_r8, &
         0.108143E+02_r8,0.109387E+02_r8,0.110647E+02_r8,0.111921E+02_r8,0.113207E+02_r8, &
         0.114508E+02_r8,0.115821E+02_r8,0.117149E+02_r8,0.118490E+02_r8,0.119847E+02_r8, &
         0.121216E+02_r8,0.122601E+02_r8,0.124002E+02_r8,0.125416E+02_r8,0.126846E+02_r8, &
         0.128290E+02_r8,0.129747E+02_r8,0.131224E+02_r8,0.132712E+02_r8,0.134220E+02_r8, &
         0.135742E+02_r8,0.137278E+02_r8,0.138831E+02_r8,0.140403E+02_r8,0.141989E+02_r8, &
         0.143589E+02_r8,0.145211E+02_r8,0.146845E+02_r8,0.148501E+02_r8,0.150172E+02_r8, &
         0.151858E+02_r8,0.153564E+02_r8,0.155288E+02_r8,0.157029E+02_r8,0.158786E+02_r8, &
         0.160562E+02_r8,0.162358E+02_r8,0.164174E+02_r8,0.166004E+02_r8,0.167858E+02_r8, &
         0.169728E+02_r8,0.171620E+02_r8,0.173528E+02_r8,0.175455E+02_r8,0.177406E+02_r8, &
         0.179372E+02_r8,0.181363E+02_r8,0.183372E+02_r8,0.185400E+02_r8,0.187453E+02_r8, &
         0.189523E+02_r8,0.191613E+02_r8,0.193728E+02_r8,0.195866E+02_r8,0.198024E+02_r8, &
         0.200200E+02_r8,0.202401E+02_r8,0.204626E+02_r8,0.206871E+02_r8,0.209140E+02_r8, &
         0.211430E+02_r8,0.213744E+02_r8,0.216085E+02_r8,0.218446E+02_r8,0.220828E+02_r8, &
         0.223241E+02_r8,0.225671E+02_r8,0.228132E+02_r8,0.230615E+02_r8,0.233120E+02_r8, &
         0.235651E+02_r8,0.238211E+02_r8,0.240794E+02_r8,0.243404E+02_r8,0.246042E+02_r8, &
         0.248704E+02_r8,0.251390E+02_r8,0.254109E+02_r8,0.256847E+02_r8,0.259620E+02_r8, &
         0.262418E+02_r8,0.265240E+02_r8,0.268092E+02_r8,0.270975E+02_r8,0.273883E+02_r8/
    DATA (ES_Aux(IES),IES=571,665) /                                      &
         0.276822E+02_r8,0.279792E+02_r8,0.282789E+02_r8,0.285812E+02_r8,0.288867E+02_r8, &
         0.291954E+02_r8,0.295075E+02_r8,0.298222E+02_r8,0.301398E+02_r8,0.304606E+02_r8, &
         0.307848E+02_r8,0.311119E+02_r8,0.314424E+02_r8,0.317763E+02_r8,0.321133E+02_r8, &
         0.324536E+02_r8,0.327971E+02_r8,0.331440E+02_r8,0.334940E+02_r8,0.338475E+02_r8, &
         0.342050E+02_r8,0.345654E+02_r8,0.349295E+02_r8,0.352975E+02_r8,0.356687E+02_r8, &
         0.360430E+02_r8,0.364221E+02_r8,0.368042E+02_r8,0.371896E+02_r8,0.375790E+02_r8, &
         0.379725E+02_r8,0.383692E+02_r8,0.387702E+02_r8,0.391744E+02_r8,0.395839E+02_r8, &
         0.399958E+02_r8,0.404118E+02_r8,0.408325E+02_r8,0.412574E+02_r8,0.416858E+02_r8, &
         0.421188E+02_r8,0.425551E+02_r8,0.429962E+02_r8,0.434407E+02_r8,0.438910E+02_r8, &
         0.443439E+02_r8,0.448024E+02_r8,0.452648E+02_r8,0.457308E+02_r8,0.462018E+02_r8, &
         0.466775E+02_r8,0.471582E+02_r8,0.476428E+02_r8,0.481313E+02_r8,0.486249E+02_r8, &
         0.491235E+02_r8,0.496272E+02_r8,0.501349E+02_r8,0.506479E+02_r8,0.511652E+02_r8, &
         0.516876E+02_r8,0.522142E+02_r8,0.527474E+02_r8,0.532836E+02_r8,0.538266E+02_r8, &
         0.543737E+02_r8,0.549254E+02_r8,0.554839E+02_r8,0.560456E+02_r8,0.566142E+02_r8, &
         0.571872E+02_r8,0.577662E+02_r8,0.583498E+02_r8,0.589392E+02_r8,0.595347E+02_r8, &
         0.601346E+02_r8,0.607410E+02_r8,0.613519E+02_r8,0.619689E+02_r8,0.625922E+02_r8, &
         0.632204E+02_r8,0.638550E+02_r8,0.644959E+02_r8,0.651418E+02_r8,0.657942E+02_r8, &
         0.664516E+02_r8,0.671158E+02_r8,0.677864E+02_r8,0.684624E+02_r8,0.691451E+02_r8, &
         0.698345E+02_r8,0.705293E+02_r8,0.712312E+02_r8,0.719398E+02_r8,0.726542E+02_r8/
    DATA (ES_Aux(IES),IES=666,760) /                                      &
         0.733754E+02_r8,0.741022E+02_r8,0.748363E+02_r8,0.755777E+02_r8,0.763247E+02_r8, &
         0.770791E+02_r8,0.778394E+02_r8,0.786088E+02_r8,0.793824E+02_r8,0.801653E+02_r8, &
         0.809542E+02_r8,0.817509E+02_r8,0.825536E+02_r8,0.833643E+02_r8,0.841828E+02_r8, &
         0.850076E+02_r8,0.858405E+02_r8,0.866797E+02_r8,0.875289E+02_r8,0.883827E+02_r8, &
         0.892467E+02_r8,0.901172E+02_r8,0.909962E+02_r8,0.918818E+02_r8,0.927760E+02_r8, &
         0.936790E+02_r8,0.945887E+02_r8,0.955071E+02_r8,0.964346E+02_r8,0.973689E+02_r8, &
         0.983123E+02_r8,0.992648E+02_r8,0.100224E+03_r8,0.101193E+03_r8,0.102169E+03_r8, &
         0.103155E+03_r8,0.104150E+03_r8,0.105152E+03_r8,0.106164E+03_r8,0.107186E+03_r8, &
         0.108217E+03_r8,0.109256E+03_r8,0.110303E+03_r8,0.111362E+03_r8,0.112429E+03_r8, &
         0.113503E+03_r8,0.114588E+03_r8,0.115684E+03_r8,0.116789E+03_r8,0.117903E+03_r8, &
         0.119028E+03_r8,0.120160E+03_r8,0.121306E+03_r8,0.122460E+03_r8,0.123623E+03_r8, &
         0.124796E+03_r8,0.125981E+03_r8,0.127174E+03_r8,0.128381E+03_r8,0.129594E+03_r8, &
         0.130822E+03_r8,0.132058E+03_r8,0.133306E+03_r8,0.134563E+03_r8,0.135828E+03_r8, &
         0.137109E+03_r8,0.138402E+03_r8,0.139700E+03_r8,0.141017E+03_r8,0.142338E+03_r8, &
         0.143676E+03_r8,0.145025E+03_r8,0.146382E+03_r8,0.147753E+03_r8,0.149133E+03_r8, &
         0.150529E+03_r8,0.151935E+03_r8,0.153351E+03_r8,0.154783E+03_r8,0.156222E+03_r8, &
         0.157678E+03_r8,0.159148E+03_r8,0.160624E+03_r8,0.162117E+03_r8,0.163621E+03_r8, &
         0.165142E+03_r8,0.166674E+03_r8,0.168212E+03_r8,0.169772E+03_r8,0.171340E+03_r8, &
         0.172921E+03_r8,0.174522E+03_r8,0.176129E+03_r8,0.177755E+03_r8,0.179388E+03_r8/
    DATA (ES_Aux(IES),IES=761,855) /                                      &
         0.181040E+03_r8,0.182707E+03_r8,0.184382E+03_r8,0.186076E+03_r8,0.187782E+03_r8, &
         0.189503E+03_r8,0.191240E+03_r8,0.192989E+03_r8,0.194758E+03_r8,0.196535E+03_r8, &
         0.198332E+03_r8,0.200141E+03_r8,0.201963E+03_r8,0.203805E+03_r8,0.205656E+03_r8, &
         0.207532E+03_r8,0.209416E+03_r8,0.211317E+03_r8,0.213236E+03_r8,0.215167E+03_r8, &
         0.217121E+03_r8,0.219087E+03_r8,0.221067E+03_r8,0.223064E+03_r8,0.225080E+03_r8, &
         0.227113E+03_r8,0.229160E+03_r8,0.231221E+03_r8,0.233305E+03_r8,0.235403E+03_r8, &
         0.237520E+03_r8,0.239655E+03_r8,0.241805E+03_r8,0.243979E+03_r8,0.246163E+03_r8, &
         0.248365E+03_r8,0.250593E+03_r8,0.252830E+03_r8,0.255093E+03_r8,0.257364E+03_r8, &
         0.259667E+03_r8,0.261979E+03_r8,0.264312E+03_r8,0.266666E+03_r8,0.269034E+03_r8, &
         0.271430E+03_r8,0.273841E+03_r8,0.276268E+03_r8,0.278722E+03_r8,0.281185E+03_r8, &
         0.283677E+03_r8,0.286190E+03_r8,0.288714E+03_r8,0.291266E+03_r8,0.293834E+03_r8, &
         0.296431E+03_r8,0.299045E+03_r8,0.301676E+03_r8,0.304329E+03_r8,0.307006E+03_r8, &
         0.309706E+03_r8,0.312423E+03_r8,0.315165E+03_r8,0.317930E+03_r8,0.320705E+03_r8, &
         0.323519E+03_r8,0.326350E+03_r8,0.329199E+03_r8,0.332073E+03_r8,0.334973E+03_r8, &
         0.337897E+03_r8,0.340839E+03_r8,0.343800E+03_r8,0.346794E+03_r8,0.349806E+03_r8, &
         0.352845E+03_r8,0.355918E+03_r8,0.358994E+03_r8,0.362112E+03_r8,0.365242E+03_r8, &
         0.368407E+03_r8,0.371599E+03_r8,0.374802E+03_r8,0.378042E+03_r8,0.381293E+03_r8, &
         0.384588E+03_r8,0.387904E+03_r8,0.391239E+03_r8,0.394604E+03_r8,0.397988E+03_r8, &
         0.401411E+03_r8,0.404862E+03_r8,0.408326E+03_r8,0.411829E+03_r8,0.415352E+03_r8/
    DATA (ES_Aux(IES),IES=856,950) /                                      &
         0.418906E+03_r8,0.422490E+03_r8,0.426095E+03_r8,0.429740E+03_r8,0.433398E+03_r8, &
         0.437097E+03_r8,0.440827E+03_r8,0.444570E+03_r8,0.448354E+03_r8,0.452160E+03_r8, &
         0.455999E+03_r8,0.459870E+03_r8,0.463765E+03_r8,0.467702E+03_r8,0.471652E+03_r8, &
         0.475646E+03_r8,0.479674E+03_r8,0.483715E+03_r8,0.487811E+03_r8,0.491911E+03_r8, &
         0.496065E+03_r8,0.500244E+03_r8,0.504448E+03_r8,0.508698E+03_r8,0.512961E+03_r8, &
         0.517282E+03_r8,0.521617E+03_r8,0.525989E+03_r8,0.530397E+03_r8,0.534831E+03_r8, &
         0.539313E+03_r8,0.543821E+03_r8,0.548355E+03_r8,0.552938E+03_r8,0.557549E+03_r8, &
         0.562197E+03_r8,0.566884E+03_r8,0.571598E+03_r8,0.576351E+03_r8,0.581131E+03_r8, &
         0.585963E+03_r8,0.590835E+03_r8,0.595722E+03_r8,0.600663E+03_r8,0.605631E+03_r8, &
         0.610641E+03_r8,0.615151E+03_r8,0.619625E+03_r8,0.624140E+03_r8,0.628671E+03_r8, &
         0.633243E+03_r8,0.637845E+03_r8,0.642465E+03_r8,0.647126E+03_r8,0.651806E+03_r8, &
         0.656527E+03_r8,0.661279E+03_r8,0.666049E+03_r8,0.670861E+03_r8,0.675692E+03_r8, &
         0.680566E+03_r8,0.685471E+03_r8,0.690396E+03_r8,0.695363E+03_r8,0.700350E+03_r8, &
         0.705381E+03_r8,0.710444E+03_r8,0.715527E+03_r8,0.720654E+03_r8,0.725801E+03_r8, &
         0.730994E+03_r8,0.736219E+03_r8,0.741465E+03_r8,0.746756E+03_r8,0.752068E+03_r8, &
         0.757426E+03_r8,0.762819E+03_r8,0.768231E+03_r8,0.773692E+03_r8,0.779172E+03_r8, &
         0.784701E+03_r8,0.790265E+03_r8,0.795849E+03_r8,0.801483E+03_r8,0.807137E+03_r8, &
         0.812842E+03_r8,0.818582E+03_r8,0.824343E+03_r8,0.830153E+03_r8,0.835987E+03_r8, &
         0.841871E+03_r8,0.847791E+03_r8,0.853733E+03_r8,0.859727E+03_r8,0.865743E+03_r8/
    DATA (ES_Aux(IES),IES=951,1045) /                                     &
         0.871812E+03_r8,0.877918E+03_r8,0.884046E+03_r8,0.890228E+03_r8,0.896433E+03_r8, &
         0.902690E+03_r8,0.908987E+03_r8,0.915307E+03_r8,0.921681E+03_r8,0.928078E+03_r8, &
         0.934531E+03_r8,0.941023E+03_r8,0.947539E+03_r8,0.954112E+03_r8,0.960708E+03_r8, &
         0.967361E+03_r8,0.974053E+03_r8,0.980771E+03_r8,0.987545E+03_r8,0.994345E+03_r8, &
         0.100120E+04_r8,0.100810E+04_r8,0.101502E+04_r8,0.102201E+04_r8,0.102902E+04_r8, &
         0.103608E+04_r8,0.104320E+04_r8,0.105033E+04_r8,0.105753E+04_r8,0.106475E+04_r8, &
         0.107204E+04_r8,0.107936E+04_r8,0.108672E+04_r8,0.109414E+04_r8,0.110158E+04_r8, &
         0.110908E+04_r8,0.111663E+04_r8,0.112421E+04_r8,0.113185E+04_r8,0.113952E+04_r8, &
         0.114725E+04_r8,0.115503E+04_r8,0.116284E+04_r8,0.117071E+04_r8,0.117861E+04_r8, &
         0.118658E+04_r8,0.119459E+04_r8,0.120264E+04_r8,0.121074E+04_r8,0.121888E+04_r8, &
         0.122709E+04_r8,0.123534E+04_r8,0.124362E+04_r8,0.125198E+04_r8,0.126036E+04_r8, &
         0.126881E+04_r8,0.127731E+04_r8,0.128584E+04_r8,0.129444E+04_r8,0.130307E+04_r8, &
         0.131177E+04_r8,0.132053E+04_r8,0.132931E+04_r8,0.133817E+04_r8,0.134705E+04_r8, &
         0.135602E+04_r8,0.136503E+04_r8,0.137407E+04_r8,0.138319E+04_r8,0.139234E+04_r8, &
         0.140156E+04_r8,0.141084E+04_r8,0.142015E+04_r8,0.142954E+04_r8,0.143896E+04_r8, &
         0.144845E+04_r8,0.145800E+04_r8,0.146759E+04_r8,0.147725E+04_r8,0.148694E+04_r8, &
         0.149672E+04_r8,0.150655E+04_r8,0.151641E+04_r8,0.152635E+04_r8,0.153633E+04_r8, &
         0.154639E+04_r8,0.155650E+04_r8,0.156665E+04_r8,0.157688E+04_r8,0.158715E+04_r8, &
         0.159750E+04_r8,0.160791E+04_r8,0.161836E+04_r8,0.162888E+04_r8,0.163945E+04_r8/
    DATA (ES_Aux(IES),IES=1046,1140) /                                    &
         0.165010E+04_r8,0.166081E+04_r8,0.167155E+04_r8,0.168238E+04_r8,0.169325E+04_r8, &
         0.170420E+04_r8,0.171522E+04_r8,0.172627E+04_r8,0.173741E+04_r8,0.174859E+04_r8, &
         0.175986E+04_r8,0.177119E+04_r8,0.178256E+04_r8,0.179402E+04_r8,0.180552E+04_r8, &
         0.181711E+04_r8,0.182877E+04_r8,0.184046E+04_r8,0.185224E+04_r8,0.186407E+04_r8, &
         0.187599E+04_r8,0.188797E+04_r8,0.190000E+04_r8,0.191212E+04_r8,0.192428E+04_r8, &
         0.193653E+04_r8,0.194886E+04_r8,0.196122E+04_r8,0.197368E+04_r8,0.198618E+04_r8, &
         0.199878E+04_r8,0.201145E+04_r8,0.202416E+04_r8,0.203698E+04_r8,0.204983E+04_r8, &
         0.206278E+04_r8,0.207580E+04_r8,0.208887E+04_r8,0.210204E+04_r8,0.211525E+04_r8, &
         0.212856E+04_r8,0.214195E+04_r8,0.215538E+04_r8,0.216892E+04_r8,0.218249E+04_r8, &
         0.219618E+04_r8,0.220994E+04_r8,0.222375E+04_r8,0.223766E+04_r8,0.225161E+04_r8, &
         0.226567E+04_r8,0.227981E+04_r8,0.229399E+04_r8,0.230829E+04_r8,0.232263E+04_r8, &
         0.233708E+04_r8,0.235161E+04_r8,0.236618E+04_r8,0.238087E+04_r8,0.239560E+04_r8, &
         0.241044E+04_r8,0.242538E+04_r8,0.244035E+04_r8,0.245544E+04_r8,0.247057E+04_r8, &
         0.248583E+04_r8,0.250116E+04_r8,0.251654E+04_r8,0.253204E+04_r8,0.254759E+04_r8, &
         0.256325E+04_r8,0.257901E+04_r8,0.259480E+04_r8,0.261073E+04_r8,0.262670E+04_r8, &
         0.264279E+04_r8,0.265896E+04_r8,0.267519E+04_r8,0.269154E+04_r8,0.270794E+04_r8, &
         0.272447E+04_r8,0.274108E+04_r8,0.275774E+04_r8,0.277453E+04_r8,0.279137E+04_r8, &
         0.280834E+04_r8,0.282540E+04_r8,0.284251E+04_r8,0.285975E+04_r8,0.287704E+04_r8, &
         0.289446E+04_r8,0.291198E+04_r8,0.292954E+04_r8,0.294725E+04_r8,0.296499E+04_r8/
    DATA (ES_Aux(IES),IES=1141,1235) /                                    &
         0.298288E+04_r8,0.300087E+04_r8,0.301890E+04_r8,0.303707E+04_r8,0.305529E+04_r8, &
         0.307365E+04_r8,0.309211E+04_r8,0.311062E+04_r8,0.312927E+04_r8,0.314798E+04_r8, &
         0.316682E+04_r8,0.318577E+04_r8,0.320477E+04_r8,0.322391E+04_r8,0.324310E+04_r8, &
         0.326245E+04_r8,0.328189E+04_r8,0.330138E+04_r8,0.332103E+04_r8,0.334073E+04_r8, &
         0.336058E+04_r8,0.338053E+04_r8,0.340054E+04_r8,0.342069E+04_r8,0.344090E+04_r8, &
         0.346127E+04_r8,0.348174E+04_r8,0.350227E+04_r8,0.352295E+04_r8,0.354369E+04_r8, &
         0.356458E+04_r8,0.358559E+04_r8,0.360664E+04_r8,0.362787E+04_r8,0.364914E+04_r8, &
         0.367058E+04_r8,0.369212E+04_r8,0.371373E+04_r8,0.373548E+04_r8,0.375731E+04_r8, &
         0.377929E+04_r8,0.380139E+04_r8,0.382355E+04_r8,0.384588E+04_r8,0.386826E+04_r8, &
         0.389081E+04_r8,0.391348E+04_r8,0.393620E+04_r8,0.395910E+04_r8,0.398205E+04_r8, &
         0.400518E+04_r8,0.402843E+04_r8,0.405173E+04_r8,0.407520E+04_r8,0.409875E+04_r8, &
         0.412246E+04_r8,0.414630E+04_r8,0.417019E+04_r8,0.419427E+04_r8,0.421840E+04_r8, &
         0.424272E+04_r8,0.426715E+04_r8,0.429165E+04_r8,0.431634E+04_r8,0.434108E+04_r8, &
         0.436602E+04_r8,0.439107E+04_r8,0.441618E+04_r8,0.444149E+04_r8,0.446685E+04_r8, &
         0.449241E+04_r8,0.451810E+04_r8,0.454385E+04_r8,0.456977E+04_r8,0.459578E+04_r8, &
         0.462197E+04_r8,0.464830E+04_r8,0.467468E+04_r8,0.470127E+04_r8,0.472792E+04_r8, &
         0.475477E+04_r8,0.478175E+04_r8,0.480880E+04_r8,0.483605E+04_r8,0.486336E+04_r8, &
         0.489087E+04_r8,0.491853E+04_r8,0.494623E+04_r8,0.497415E+04_r8,0.500215E+04_r8, &
         0.503034E+04_r8,0.505867E+04_r8,0.508707E+04_r8,0.511568E+04_r8,0.514436E+04_r8/
    DATA (ES_Aux(IES),IES=1236,1330) /                                    &
         0.517325E+04_r8,0.520227E+04_r8,0.523137E+04_r8,0.526068E+04_r8,0.529005E+04_r8, &
         0.531965E+04_r8,0.534939E+04_r8,0.537921E+04_r8,0.540923E+04_r8,0.543932E+04_r8, &
         0.546965E+04_r8,0.550011E+04_r8,0.553064E+04_r8,0.556139E+04_r8,0.559223E+04_r8, &
         0.562329E+04_r8,0.565449E+04_r8,0.568577E+04_r8,0.571727E+04_r8,0.574884E+04_r8, &
         0.578064E+04_r8,0.581261E+04_r8,0.584464E+04_r8,0.587692E+04_r8,0.590924E+04_r8, &
         0.594182E+04_r8,0.597455E+04_r8,0.600736E+04_r8,0.604039E+04_r8,0.607350E+04_r8, &
         0.610685E+04_r8,0.614036E+04_r8,0.617394E+04_r8,0.620777E+04_r8,0.624169E+04_r8, &
         0.627584E+04_r8,0.631014E+04_r8,0.634454E+04_r8,0.637918E+04_r8,0.641390E+04_r8, &
         0.644887E+04_r8,0.648400E+04_r8,0.651919E+04_r8,0.655467E+04_r8,0.659021E+04_r8, &
         0.662599E+04_r8,0.666197E+04_r8,0.669800E+04_r8,0.673429E+04_r8,0.677069E+04_r8, &
         0.680735E+04_r8,0.684415E+04_r8,0.688104E+04_r8,0.691819E+04_r8,0.695543E+04_r8, &
         0.699292E+04_r8,0.703061E+04_r8,0.706837E+04_r8,0.710639E+04_r8,0.714451E+04_r8, &
         0.718289E+04_r8,0.722143E+04_r8,0.726009E+04_r8,0.729903E+04_r8,0.733802E+04_r8, &
         0.737729E+04_r8,0.741676E+04_r8,0.745631E+04_r8,0.749612E+04_r8,0.753602E+04_r8, &
         0.757622E+04_r8,0.761659E+04_r8,0.765705E+04_r8,0.769780E+04_r8,0.773863E+04_r8, &
         0.777975E+04_r8,0.782106E+04_r8,0.786246E+04_r8,0.790412E+04_r8,0.794593E+04_r8, &
         0.798802E+04_r8,0.803028E+04_r8,0.807259E+04_r8,0.811525E+04_r8,0.815798E+04_r8, &
         0.820102E+04_r8,0.824427E+04_r8,0.828757E+04_r8,0.833120E+04_r8,0.837493E+04_r8, &
         0.841895E+04_r8,0.846313E+04_r8,0.850744E+04_r8,0.855208E+04_r8,0.859678E+04_r8/
    DATA (ES_Aux(IES),IES=1331,1425) /                                    &
         0.864179E+04_r8,0.868705E+04_r8,0.873237E+04_r8,0.877800E+04_r8,0.882374E+04_r8, &
         0.886979E+04_r8,0.891603E+04_r8,0.896237E+04_r8,0.900904E+04_r8,0.905579E+04_r8, &
         0.910288E+04_r8,0.915018E+04_r8,0.919758E+04_r8,0.924529E+04_r8,0.929310E+04_r8, &
         0.934122E+04_r8,0.938959E+04_r8,0.943804E+04_r8,0.948687E+04_r8,0.953575E+04_r8, &
         0.958494E+04_r8,0.963442E+04_r8,0.968395E+04_r8,0.973384E+04_r8,0.978383E+04_r8, &
         0.983412E+04_r8,0.988468E+04_r8,0.993534E+04_r8,0.998630E+04_r8,0.100374E+05_r8, &
         0.100888E+05_r8,0.101406E+05_r8,0.101923E+05_r8,0.102444E+05_r8,0.102966E+05_r8, &
         0.103492E+05_r8,0.104020E+05_r8,0.104550E+05_r8,0.105082E+05_r8,0.105616E+05_r8, &
         0.106153E+05_r8,0.106693E+05_r8,0.107234E+05_r8,0.107779E+05_r8,0.108325E+05_r8, &
         0.108874E+05_r8,0.109425E+05_r8,0.109978E+05_r8,0.110535E+05_r8,0.111092E+05_r8, &
         0.111653E+05_r8,0.112217E+05_r8,0.112782E+05_r8,0.113350E+05_r8,0.113920E+05_r8, &
         0.114493E+05_r8,0.115070E+05_r8,0.115646E+05_r8,0.116228E+05_r8,0.116809E+05_r8, &
         0.117396E+05_r8,0.117984E+05_r8,0.118574E+05_r8,0.119167E+05_r8,0.119762E+05_r8, &
         0.120360E+05_r8,0.120962E+05_r8,0.121564E+05_r8,0.122170E+05_r8,0.122778E+05_r8, &
         0.123389E+05_r8,0.124004E+05_r8,0.124619E+05_r8,0.125238E+05_r8,0.125859E+05_r8, &
         0.126484E+05_r8,0.127111E+05_r8,0.127739E+05_r8,0.128372E+05_r8,0.129006E+05_r8, &
         0.129644E+05_r8,0.130285E+05_r8,0.130927E+05_r8,0.131573E+05_r8,0.132220E+05_r8, &
         0.132872E+05_r8,0.133526E+05_r8,0.134182E+05_r8,0.134842E+05_r8,0.135503E+05_r8, &
         0.136168E+05_r8,0.136836E+05_r8,0.137505E+05_r8,0.138180E+05_r8,0.138854E+05_r8/
    DATA (ES_Aux(IES),IES=1426,1520) /                                    &
         0.139534E+05_r8,0.140216E+05_r8,0.140900E+05_r8,0.141588E+05_r8,0.142277E+05_r8, &
         0.142971E+05_r8,0.143668E+05_r8,0.144366E+05_r8,0.145069E+05_r8,0.145773E+05_r8, &
         0.146481E+05_r8,0.147192E+05_r8,0.147905E+05_r8,0.148622E+05_r8,0.149341E+05_r8, &
         0.150064E+05_r8,0.150790E+05_r8,0.151517E+05_r8,0.152250E+05_r8,0.152983E+05_r8, &
         0.153721E+05_r8,0.154462E+05_r8,0.155205E+05_r8,0.155952E+05_r8,0.156701E+05_r8, &
         0.157454E+05_r8,0.158211E+05_r8,0.158969E+05_r8,0.159732E+05_r8,0.160496E+05_r8, &
         0.161265E+05_r8,0.162037E+05_r8,0.162811E+05_r8,0.163589E+05_r8,0.164369E+05_r8, &
         0.165154E+05_r8,0.165942E+05_r8,0.166732E+05_r8,0.167526E+05_r8,0.168322E+05_r8, &
         0.169123E+05_r8,0.169927E+05_r8,0.170733E+05_r8,0.171543E+05_r8,0.172356E+05_r8, &
         0.173173E+05_r8,0.173993E+05_r8,0.174815E+05_r8,0.175643E+05_r8,0.176471E+05_r8, &
         0.177305E+05_r8,0.178143E+05_r8,0.178981E+05_r8,0.179826E+05_r8,0.180671E+05_r8, &
         0.181522E+05_r8,0.182377E+05_r8,0.183232E+05_r8,0.184093E+05_r8,0.184955E+05_r8, &
         0.185823E+05_r8,0.186695E+05_r8,0.187568E+05_r8,0.188447E+05_r8,0.189326E+05_r8, &
         0.190212E+05_r8,0.191101E+05_r8,0.191991E+05_r8,0.192887E+05_r8,0.193785E+05_r8, &
         0.194688E+05_r8,0.195595E+05_r8,0.196503E+05_r8,0.197417E+05_r8,0.198332E+05_r8, &
         0.199253E+05_r8,0.200178E+05_r8,0.201105E+05_r8,0.202036E+05_r8,0.202971E+05_r8, &
         0.203910E+05_r8,0.204853E+05_r8,0.205798E+05_r8,0.206749E+05_r8,0.207701E+05_r8, &
         0.208659E+05_r8,0.209621E+05_r8,0.210584E+05_r8,0.211554E+05_r8,0.212524E+05_r8, &
         0.213501E+05_r8,0.214482E+05_r8,0.215465E+05_r8,0.216452E+05_r8,0.217442E+05_r8/
    DATA (ES_Aux(IES),IES=1521,1552) /                                    &
         0.218439E+05_r8,0.219439E+05_r8,0.220440E+05_r8,0.221449E+05_r8,0.222457E+05_r8, &
         0.223473E+05_r8,0.224494E+05_r8,0.225514E+05_r8,0.226542E+05_r8,0.227571E+05_r8, &
         0.228606E+05_r8,0.229646E+05_r8,0.230687E+05_r8,0.231734E+05_r8,0.232783E+05_r8, &
         0.233839E+05_r8,0.234898E+05_r8,0.235960E+05_r8,0.237027E+05_r8,0.238097E+05_r8, &
         0.239173E+05_r8,0.240254E+05_r8,0.241335E+05_r8,0.242424E+05_r8,0.243514E+05_r8, &
         0.244611E+05_r8,0.245712E+05_r8,0.246814E+05_r8,0.247923E+05_r8,0.249034E+05_r8, &
         0.250152E+05_r8,0.250152E+05_r8/
    !

    ES=ES_Aux
  END SUBROUTINE qsat_data
  ! ======================================================================


  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Saturation Specific Humidity/mixing Scheme(Qsat):Vapour to Liquid/Ice
  !
  ! Subroutine Interface:
  SUBROUTINE QSAT_mix (                                             &
       !      Output field
       &  QmixS                                                           &
       !      Input fields
       &, T, P                                                            &
       !      Array dimensions
       &, NPNTS                                                           &
       !      logical control
       &, lq_mix                                                          &
       &, mskant                                                          &
       
       &  )
    IMPLICIT NONE
    !
    ! Purpose:
    !   Returns a saturation specific humidity or mixing ratio given a
    !   temperature and pressure using the saturation vapour pressure
    !   calculated using the Goff-Gratch formulae, adopted by the WMO as
    !   taken from Landolt-Bornstein, 1987 Numerical Data and Functional
    !   Relationships in Science and Technolgy. Group V/vol 4B meteorology.
    !   Phyiscal and Chemical properties or air, P35.
    !
    !   Values in the lookup table are over water above 0 deg C and over ice
    !   below this temperature.
    !
    ! Method:
    !   Uses lookup tables to find eSAT, calculates qSAT directly from that.
    !
    ! Current Owner of Code: OWNER OF LARGE_SCALE CLOUD CODE.
    !
    ! Code Description:
    !   Language: FORTRAN 77  + common extensions also in Fortran90.
    !   This code is written to UMDP3 version 6 programming standards.
    !
    !   Documentation: UMDP No.29
    !
    ! Declarations:
    !
    !  Global Variables:----------------------------------------------------

    !-----------------------------------------------------------------------
    ! Subroutine Arguments
    !-----------------------------------------------------------------------
    !
    ! arguments with intent in. ie: input variables.
    !
    INTEGER, INTENT(in) :: npnts    ! Points (=horizontal dimensions) being processed by qSAT scheme.
    REAL(KIND=r8), INTENT(in)  :: T(npnts)      !  Temperature (K).
    REAL(KIND=r8), INTENT(in)  :: P(npnts)      !  Pressure (Pa).
    INTEGER(KIND=i8), INTENT(in)  :: mskant(npnts)

    LOGICAL, INTENT(in)  :: lq_mix      !  .true. return qsat as a mixing ratio
    !  .false. return qsat as a specific humidity

    !
    ! arguments with intent out
    !
    REAL(KIND=r8), INTENT(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturation
    ! specific humidity at temperature T and pressure
    ! P (kg/kg).

    !-----------------------------------------------------------------------
    !  Local scalars
    !-----------------------------------------------------------------------
    INTEGER :: itable, itable_p1    ! Work variables

    REAL(KIND=r8) :: atable, atable_p1       ! Work variables

    REAL(KIND=r8) :: fsubw, fsubw_p1
    ! FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE IN A PURE
    ! WATER SYSTEM TO SAT VAPOUR PRESSURE IN AIR.

    REAL(KIND=r8) :: TT, TT_P1

    REAL(KIND=r8), PARAMETER :: R_delta_T = 1.0_r8/delta_T

    INTEGER :: I, II
    !
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !
    ! loop over points
    !
    DO i = 1, npnts-1,2
       IF(mskant(i) == 1_i8)THEN

          !      Compute the factor that converts from sat vapour pressure in a
          !      pure water system to sat vapour pressure in air, fsubw.
          !      This formula is taken from equation A4.7 of Adrian Gill's book:
          !      Atmosphere-Ocean Dynamics. Note that his formula works in terms
          !      of pressure in mb and temperature in Celsius, so conversion of
          !      units leads to the slightly different equation used here.

          fsubw    = 1.0_r8 + 1.0E-8_r8 * P(I)   * ( 4.5_r8 +                      &
               6.0E-4_r8*( T(I) - zerodegC )   * ( T(I) - zerodegC ) )
          fsubw_p1 = 1.0_r8 + 1.0E-8_r8 * P(I+1) * ( 4.5_r8 +                      &
               6.0E-4_r8*( T(I+1) - zerodegC ) * ( T(I+1) - zerodegC ) )

          !      Use the lookup table to find saturated vapour pressure, and store
          !      it in qmixs.
          !
          TT = MAX(T_low,T(I))
          TT = MIN(T_high,TT)
          atable = (TT - T_low + delta_T) * R_delta_T
          itable = atable
          atable = atable - itable

          TT_p1 = MAX(T_LOW,T(I+1))
          TT_p1 = MIN(T_HIGH,TT_p1)
          ATABLE_p1 = (TT_p1 - T_LOW + DELTA_T) * R_DELTA_T
          ITABLE_p1 = ATABLE_p1
          ATABLE_p1 = ATABLE_p1 - ITABLE_p1

          QmixS(I)   = (1.0_r8 - atable)   *ES(itable) + atable*ES(itable+1)
          QmixS(I+1) = (1.0_r8 - atable_p1)*ES(itable_p1) +                 &
               atable_p1*ES(itable_p1+1)

          !      Multiply by fsubw to convert to saturated vapour pressure in air
          !      (equation A4.6 of Adrian Gill's book).
          !
          QmixS(I)   = QmixS(I)   * fsubw
          QmixS(I+1) = QmixS(I+1) * fsubw_p1
          !
          !      Now form the accurate expression for qmixs, which is a rearranged
          !      version of equation A4.3 of Gill's book.
          !
          !-----------------------------------------------------------------------
          ! For mixing ratio,  rsat = epsilon *e/(p-e)
          ! e - saturation vapour pressure
          ! Note applying the fix to qsat for specific humidity at low pressures
          ! is not possible, this implies mixing ratio qsat tends to infinity.
          ! If the pressure is very low then the mixing ratio value may become
          ! very large.
          !-----------------------------------------------------------------------

          IF (lq_mix) THEN

             QmixS(I)   = ( epsilonk*QmixS(I) ) /                            &
                  ( MAX(P(I),  1.1_r8*QmixS(I))   - QmixS(I) )
             QmixS(I+1) = ( epsilonk*QmixS(I+1) ) /                          &
                  ( MAX(P(I+1),1.1_r8*QmixS(I+1)) - QmixS(I+1) )

             !-----------------------------------------------------------------------
             ! For specific humidity,   qsat = epsilon*e/(p-(1-epsilon)e)
             !
             ! Note that at very low pressure we apply a fix, to prevent a
             ! singularity (qsat tends to 1. kg/kg).
             !-----------------------------------------------------------------------
          ELSE

             QmixS(I)   = ( epsilonk*QmixS(I) ) /                           &
                  ( MAX(P(I),  QmixS(I))   - one_minus_epsilon*QmixS(I) )
             QmixS(I+1) = ( epsilonk*QmixS(I+1) ) /                         &
                  ( MAX(P(I+1),QmixS(I+1)) - one_minus_epsilon*QmixS(I+1))

          ENDIF  ! test on lq_mix
          !
       END IF
    END DO  ! Npnts_do_1

    II = I
    DO I = II, NPNTS
       IF(mskant(i) == 1_i8)THEN

          fsubw = 1.0_r8 + 1.0E-8_r8*P(I)*( 4.5_r8 +                               &
               &    6.0E-4_r8*( T(I) - zerodegC )*( T(I) - zerodegC ) )
          !
          TT = MAX(T_low,T(I))
          TT = MIN(T_high,TT)
          atable = (TT - T_low + delta_T) * R_delta_T
          itable = atable
          atable = atable - itable

          QmixS(I) = (1.0_r8 - atable)*ES(itable) + atable*ES(itable+1)
          QmixS(I) = QmixS(I) * fsubw

          IF (lq_mix) THEN
             QmixS(I) = ( epsilonk*QmixS(I) ) /                              &
                  &              ( MAX(P(I),1.1_r8*QmixS(i)) - QmixS(I) )
          ELSE
             QmixS(I) = ( epsilonk*QmixS(I) ) /                              &
                  &          ( MAX(P(I),QmixS(I)) - one_minus_epsilon*QmixS(I) )
          ENDIF  ! test on lq_mix
       END IF
    END DO  ! Npnts_do_1

    RETURN
  END SUBROUTINE QSAT_mix
  ! ======================================================================

END MODULE Sfc_SeaFlux_UKME_Model


!PROGRAM Main
!  USE Sfc_SeaFlux_UKME_Model
! IMPLICIT NONE
!
!END PROGRAM Main
