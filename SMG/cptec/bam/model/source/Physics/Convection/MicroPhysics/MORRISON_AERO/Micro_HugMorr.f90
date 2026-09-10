
! This version is tagged as V3.3 of the microphysics code (this is not WRFV3.3!!!!!!)

! NOTE: THIS CODE CONTAINS OPTION FOR PREDICTED DROPLET CONCENTRATION

! THIS MODULE CONTAINS THE TWO-MOMENT MICROPHYSICS CODE DESCRIBED BY
!     MORRISON ET AL. (2009, MWR)

! CHANGES FOR V3.2, RELATIVE TO MOST RECENT (BUG-FIX) CODE FOR V3.1

! 1) ADDED ACCELERATED MELTING OF GRAUPEL/SNOW DUE TO COLLISION WITH RAIN, FOLLOWING LIN ET AL. (1983)
! 2) INCREASED MINIMUM LAMBDA FOR RAIN, AND ADDED RAIN DROP BREAKUP FOLLOWING MODIFIED VERSION
!     OF VERLINDE AND COTTON (1993)
! 3) CHANGE MINIMUM ALLOWED MIXING RATIOS IN DRY CONDITIONS (RH < 90%), THIS IMPROVES RADAR REFLECTIIVITY
!     IN LOW REFLECTIVITY REGIONS
! 4) BUG FIX TO MAXIMUM ALLOWED PARTICLE FALLSPEEDS AS A FUNCTION OF AIR DENSITY
! 5) BUG FIX TO CALCULATION OF LIQUID WATER SATURATION VAPOR PRESSURE (CHANGE IS VERY MINOR)
! 6) INCLUDE WRF CONSTANTS PER SUGGESTION OF JIMY

! changes for consistency with WRFV3.3 microphysics (updated version)
! minor revisions by Andy Ackerman
! 1) replaced kinematic with dynamic viscosity 
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! Additional fixes
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NO SUBTRACTED SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SINCE
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS 
! 7) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 8) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 9) BUG FIX TO IGRAUP SWITCH FOR NO GRAUPEL/HAIL

! NOTE!!! THERE ARE ADDITIONAL CHANGES FOR V3.3 DUE TO COUPLING WITH WRF-CHEM,
! THESE ARE NOT INCLUDED HERE

! hm, fixes 3/4/13 -- for version 3.3

! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted
! 2) removed second initialization of evpms (non-answer-changing)
! 3) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass 
!    collected by rain (very minor change)
! 4) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!    cm-3. This was done to address the problem of excessive and persistent
!    anvil cirrus produced by the scheme, and was found to greatly improve forecasts over
!    at convection-permitting scales over the central U.S. in summertime.
! 5) some changes to comments

! hm, changes 7/25/13 for version 3.4

! 1) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!    Hailong Wang (PNNL)
! 2) very minor fix to immersion freezing rate formulation (negligible impact)
! 3) clarifications to code comments
! 4) minor change to shedding of rain, remove limit so that the number of
!    collected drops can smaller than number of shed drops
! 5) change of specific heat of liquid water from 4218 to 4187 J/kg/K

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.

MODULE Micro_HugMorr
  IMPLICIT NONE
  INTEGER      , PARAMETER :: r8  = SELECTED_REAL_KIND(P=13,R=300)
  REAL(KIND=r8)    , PARAMETER :: r_d          = 287.04_r8
  REAL(KIND=r8)    , PARAMETER :: r_v          = 461.6_r8
  REAL(KIND=r8)    , PARAMETER :: cp           = 1004.6_r8



  REAL(KIND=r8), PARAMETER :: PI = 3.1415926535897932384626434_r8
  REAL(KIND=r8), PARAMETER :: SQRTPI = 0.9189385332046727417803297_r8
  REAL(R8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_R8       ! molecular weight dry air ~ kg/kmole  
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(R8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole  
  REAL(R8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
  REAL(R8),PARAMETER :: SHR_CONST_G      = 9.80616_R8      ! acceleration of gravity ~ m/s^2

  REAL(R8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole

  REAL(R8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg

  REAL(r8),PARAMETER :: rair   = SHR_CONST_RDAIR    ! Gas constant for dry air (J/K/kg)
  REAL(r8),PARAMETER :: gravit = SHR_CONST_G      ! gravitational acceleration
  REAL(r8),PARAMETER :: zvir   = SHR_CONST_RWV/SHR_CONST_RDAIR - 1          ! rh2o/rair - 1
  LOGICAL, PARAMETER :: f_qndrop=.FALSE.

  INTEGER,     PARAMETER  :: PBL=1

  ! amy 
  ! for YSU pbl scheme:
  ! coupling as in WRF2
  !IF (PBL.ne.1) THEN
  ! for MYJ pbl scheme or 3D TKE:
  !    WVAR(I,K)     = (0.667*tke(i,k))**0.5
  !ELSE
  ! for YSU pbl scheme:
  !    WVAR(I,K) = KZH(I,K)/20.
  !END IF 

  !   PUBLIC  ::  MP_MORR_TWO_MOMENT
  PUBLIC  ::  POLYSVP

  PRIVATE :: GAMMA, DERF1
  PRIVATE :: PI, SQRTPI
  PRIVATE :: MORR_TWO_MOMENT_MICRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SWITCHES FOR MICROPHYSICS SCHEME
  ! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
  ! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

  INTEGER, PRIVATE ::  IACT

  ! INUM = 0, PREDICT DROPLET CONCENTRATION
  ! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

  INTEGER, PRIVATE ::  INUM

  ! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
  REAL(KIND=r8), PRIVATE ::      NDCNST

  ! SWITCH FOR LIQUID-ONLY RUN
  ! ILIQ = 0, INCLUDE ICE
  ! ILIQ = 1, LIQUID ONLY, NO ICE

  INTEGER, PRIVATE ::  ILIQ

  ! SWITCH FOR ICE NUCLEATION
  ! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
  !      = 1, USE MPACE OBSERVATIONS

  INTEGER, PRIVATE ::  INUC

  ! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
  !             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
  !             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
  !             NON-EQULIBRIUM SUPERSATURATION, 
  !             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
  ! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
  !             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
  !             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
  !             SUPERSATURATION, BASED ON THE 
  !             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
  !             AT THE GRID POINT

  ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

  INTEGER, PRIVATE ::  IBASE

  ! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
  ! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
  ! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

  INTEGER, PRIVATE ::  ISUB      

  ! SWITCH FOR GRAUPEL/NO GRAUPEL
  ! IGRAUP = 0, INCLUDE GRAUPEL
  ! IGRAUP = 1, NO GRAUPEL

  INTEGER, PRIVATE ::  IGRAUP

  ! HM ADDED NEW OPTION FOR HAIL
  ! SWITCH FOR HAIL/GRAUPEL
  ! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
  ! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

  INTEGER, PRIVATE ::  IHAIL

  ! CLOUD MICROPHYSICS CONSTANTS

  REAL(KIND=r8), PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
  REAL(KIND=r8), PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
  !REAL(KIND=r8), PRIVATE ::      R           ! GAS CONSTANT FOR AIR
  !REAL(KIND=r8), PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
  !REAL(KIND=r8), PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
  REAL(KIND=r8), PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
  REAL(KIND=r8), PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
  REAL(KIND=r8), PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
  REAL(KIND=r8), PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
  REAL(KIND=r8), PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
  REAL(KIND=r8), PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
  REAL(KIND=r8), PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
  REAL(KIND=r8), PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
  REAL(KIND=r8), PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
  REAL(KIND=r8), PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
  REAL(KIND=r8), PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
  REAL(KIND=r8), PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
  REAL(KIND=r8), PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
  REAL(KIND=r8), PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
  REAL(KIND=r8), PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
  !REAL(KIND=r8), PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
  REAL(KIND=r8), PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
  REAL(KIND=r8), PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
  REAL(KIND=r8), PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
  REAL(KIND=r8), PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
  REAL(KIND=r8), PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
  ! hm, add for V3.2
  REAL(KIND=r8), PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

  ! CCN SPECTRA FOR IACT = 1

  REAL(KIND=r8), PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
  REAL(KIND=r8), PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

  ! AEROSOL PARAMETERS FOR IACT = 2

  REAL(KIND=r8), PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
  REAL(KIND=r8), PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
  REAL(KIND=r8), PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
  REAL(KIND=r8), PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
  REAL(KIND=r8), PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
  REAL(KIND=r8), PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
  REAL(KIND=r8), PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
  REAL(KIND=r8), PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
  REAL(KIND=r8), PRIVATE ::      BACT    ! ACTIVATION PARAMETER
  REAL(KIND=r8), PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
  REAL(KIND=r8), PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
  REAL(KIND=r8), PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
  REAL(KIND=r8), PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
  REAL(KIND=r8), PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
  REAL(KIND=r8), PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
  REAL(KIND=r8), PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
  REAL(KIND=r8), PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
  REAL(KIND=r8), PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
  REAL(KIND=r8), PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2     
  REAL(KIND=r8), PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
  REAL(KIND=r8), PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING

  ! CONSTANTS TO IMPROVE EFFICIENCY

  REAL(KIND=r8), PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
  REAL(KIND=r8), PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
  REAL(KIND=r8), PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
  REAL(KIND=r8), PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
  REAL(KIND=r8), PRIVATE :: CONS41


  !+---+-----------------------------------------------------------------+
  !..This set of routines facilitates computing radar reflectivity.
  !.. This module is more library code whereas the individual microphysics
  !.. schemes contains specific details needed for the final computation,
  !.. so refer to location within each schemes calling the routine named
  !.. rayleigh_soak_wetgraupel.
  !.. The bulk of this code originated from Ulrich Blahak (Germany) and
  !.. was adapted to WRF by G. Thompson.  This version of code is only
  !.. intended for use when Rayleigh scattering principles dominate and
  !.. is not intended for wavelengths in which Mie scattering is a
  !.. significant portion.  Therefore, it is well-suited to use with
  !.. 5 or 10 cm wavelength like USA NEXRAD radars.
  !.. This code makes some rather simple assumptions about water
  !.. coating on outside of frozen species (snow/graupel).  Fraction of
  !.. meltwater is simply the ratio of mixing ratio below melting level
  !.. divided by mixing ratio at level just above highest T>0C.  Also,
  !.. immediately 90% of the melted water exists on the ice's surface
  !.. and 10% is embedded within ice.  No water is "shed" at all in these
  !.. assumptions. The code is quite slow because it does the reflectivity
  !.. calculations based on 50 individual size bins of the distributions.
  !+---+-----------------------------------------------------------------+
  PRIVATE :: rayleigh_soak_wetgraupel
  PRIVATE :: radar_init
  PRIVATE :: m_complex_water_ray
  PRIVATE :: m_complex_ice_maetzler
  PRIVATE :: m_complex_maxwellgarnett
  PRIVATE :: get_m_mix_nested
  PRIVATE :: get_m_mix
  PRIVATE :: WGAMMA
  PRIVATE :: GAMMLN


  INTEGER, PARAMETER, PUBLIC:: nrbins = 50
  REAL(KIND=r8), DIMENSION(nrbins+1), PUBLIC:: xxDx
  REAL(KIND=r8), DIMENSION(nrbins), PUBLIC:: xxDs,xdts,xxDg,xdtg
  REAL(KIND=r8), PARAMETER, PUBLIC:: lamda_radar = 0.10_r8           ! in meters
  REAL(KIND=r8), PUBLIC :: K_w, PI5, lamda4
  COMPLEX*16      , PUBLIC :: m_w_0, m_i_0
  REAL(KIND=r8), DIMENSION(nrbins+1), PUBLIC:: simpson
  REAL(KIND=r8), DIMENSION(3), PARAMETER, PUBLIC:: basis =       &
       (/1.0_r8/3.0_r8, 4.0_r8/3.0_r8, 1.0_r8/3.0_r8/)
  REAL(KIND=r8), DIMENSION(4), PUBLIC:: xcre, xcse, xcge, xcrg, xcsg, xcgg
  REAL(KIND=r8), PUBLIC:: xam_r, xbm_r, xmu_r, xobmr
  REAL(KIND=r8), PUBLIC:: xam_s, xbm_s, xmu_s, xoams, xobms, xocms
  REAL(KIND=r8), PUBLIC:: xam_g, xbm_g, xmu_g, xoamg, xobmg, xocmg
  REAL(KIND=r8), PUBLIC:: xorg2, xosg2, xogg2

  INTEGER, PARAMETER, PUBLIC:: slen = 20
  CHARACTER(len=slen), PUBLIC::                                     &
       mixingrulestring_s, matrixstring_s, inclusionstring_s,    &
       hoststring_s, hostmatrixstring_s, hostinclusionstring_s,  &
       mixingrulestring_g, matrixstring_g, inclusionstring_g,    &
       hoststring_g, hostmatrixstring_g, hostinclusionstring_g

  !..Single melting snow/graupel particle 90% meltwater on external sfc
  REAL(KIND=r8), PARAMETER:: melt_outside_s = 0.90_r8
  REAL(KIND=r8), PARAMETER:: melt_outside_g = 0.90_r8

  CHARACTER*256:: radar_debug


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Init_Micro_HugMorr(ibMax,kMax,jbMax,EFFCS,EFFIS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
    ! NEEDED BY THE MICROPHYSICS SCHEME.
    ! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    INTEGER      , INTENT(IN   ) :: ibMax,kMax,jbMax
    REAL(KIND=r8), INTENT(OUT  ) :: EFFCS(ibMax,kMax,jbMax)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFIS(ibMax,kMax,jbMax)
    REAL(KIND=r8), PARAMETER :: reimin  = 10.0_r8 ! Minimum of Ice particle efective radius (microns)  

    !INTEGER n,i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
    ! SET PRIOR TO CODE COMPILATION

    ! INUM = 0, PREDICT DROPLET CONCENTRATION
    ! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

    INUM = 0

    ! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)

    NDCNST = 250.0_r8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
    ! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
    ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

    IACT = 2

    ! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
    !             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
    !             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
    !             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, 
    !             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
    ! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
    !             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
    !             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
    !             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE 
    !             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
    !             AT THE GRID POINT

    ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

    IBASE = 2

    ! INCLUDE SUB-GRID VERTICAL VELOCITY (standard deviation of w) IN DROPLET ACTIVATION
    ! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
    ! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

    ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

    ISUB = 0      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! SWITCH FOR LIQUID-ONLY RUN
    ! ILIQ = 0, INCLUDE ICE
    ! ILIQ = 1, LIQUID ONLY, NO ICE

    ILIQ = 0

    ! SWITCH FOR ICE NUCLEATION
    ! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
    !      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

    INUC = 0

    ! SWITCH FOR GRAUPEL/HAIL NO GRAUPEL/HAIL
    ! IGRAUP = 0, INCLUDE GRAUPEL/HAIL
    ! IGRAUP = 1, NO GRAUPEL/HAIL

    IGRAUP = 0

    ! HM ADDED 11/7/07
    ! SWITCH FOR HAIL/GRAUPEL
    ! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
    ! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL
    ! NOTE ---> RECOMMEND IHAIL = 1 FOR CONTINENTAL DEEP CONVECTION

    IHAIL = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SET PHYSICAL CONSTANTS

    ! FALLSPEED PARAMETERS (V=AD^B)
    AI = 700.0_r8
    AC = 3.E7_r8
    AS = 11.72_r8
    AR = 841.99667_r8
    BI = 1.0_r8
    BC = 2.0_r8
    BS = 0.41_r8
    BR = 0.8_r8
    IF (IHAIL.EQ.0) THEN
       AG = 19.3_r8
       BG = 0.37_r8
    ELSE ! (MATSUN AND HUGGINS 1980)
       AG = 114.5_r8 
       BG = 0.5_r8
    END IF

    ! CONSTANTS AND PARAMETERS
    !R = 287.15_r8
    !RV = 461.5_r8
    !CP = 1005.0_r8
    RHOSU = 85000.0_r8/(287.15_r8*273.15_r8)
    RHOW = 997.0_r8
    RHOI = 500.0_r8
    RHOSN = 100.0_r8
    IF (IHAIL.EQ.0) THEN
       RHOG = 400.0_r8
    ELSE
       RHOG = 900.0_r8
    END IF
    AIMM = 0.66_r8
    BIMM = 100.0_r8
    ECR = 1.0_r8
    DCS = 125.E-6_r8
    MI0 = 4.0_r8/3.0_r8*PI*RHOI*(10.E-6_r8)**3
    MG0 = 1.6E-10_r8
    F1S = 0.86_r8
    F2S = 0.28_r8
    F1R = 0.78_r8
    !         F2R = 0.32_r8
    ! AA revision 4/1/11
    F2R = 0.308_r8
    !G = 9.806_r8
    QSMALL = 1.E-14_r8
    EII = 0.1_r8
    ECI = 0.7_r8
    ! HM, ADD FOR V3.2
    ! hm, 7/23/13
    !         CPW = 4218.0_r8
    CPW = 4187.0_r8


    ! SIZE DISTRIBUTION PARAMETERS

    CI = RHOI*PI/6.0_r8
    DI = 3.0_r8
    CS = RHOSN*PI/6.0_r8
    DS = 3.0_r8
    CG = RHOG*PI/6.0_r8
    DG = 3.0_r8

    ! RADIUS OF CONTACT NUCLEI
    RIN = 0.1E-6_r8

    MMULT = 4.0_r8/3.0_r8*PI*RHOI*(5.E-6_r8)**3

    ! SIZE LIMITS FOR LAMBDA

    LAMMAXI = 1.0_r8/1.E-6_r8
    LAMMINI = 1.0_r8/(2.0_r8*DCS+100.E-6_r8)
    LAMMAXR = 1.0_r8/20.E-6_r8
    !         LAMMINR = 1.0_r8/500.E-6_r8
    LAMMINR = 1.0_r8/2800.E-6_r8

    LAMMAXS = 1.0_r8/10.E-6_r8
    LAMMINS = 1.0_r8/2000.E-6_r8
    LAMMAXG = 1.0_r8/20.E-6_r8
    LAMMING = 1.0_r8/2000.E-6_r8

    ! CCN SPECTRA FOR IACT = 1

    ! MARITIME
    ! MODIFIED FROM RASMUSSEN ET AL. 2002
    ! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

    K1 = 0.4_r8
    C1 = 120.0_r8 

    ! CONTINENTAL

    !              K1 = 0.5_r8
    !              C1 = 1000.0_r8 

    ! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
    ! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

    MW = 0.018_r8
    OSM = 1.0_r8
    VI = 3.0_r8
    EPSM = 0.7_r8
    RHOA = 1777.0_r8
    MAP = 0.132_r8
    MA = 0.0284_r8
    RR = 8.3187_r8
    BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

    ! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE 
    ! (see morrison et al. 2007, JGR)
    ! MODE 1

    RM1 = 0.052E-6_r8
    SIG1 = 2.04_r8
    NANEW1 = 72.2E6_r8
    F11 = 0.5_r8*EXP(2.5_r8*(LOG(SIG1))**2)
    F21 = 1.0_r8+0.25_r8*LOG(SIG1)

    ! MODE 2

    RM2 = 1.3E-6_r8
    SIG2 = 2.5_r8
    NANEW2 = 1.8E6_r8
    F12 = 0.5_r8*EXP(2.5_r8*(LOG(SIG2))**2)
    F22 = 1.0_r8+0.25_r8*LOG(SIG2)

    ! CONSTANTS FOR EFFICIENCY

    CONS1=GAMMA(1.0_r8+DS)*CS
    CONS2=GAMMA(1.0_r8+DG)*CG
    CONS3=GAMMA(4.0_r8+BS)/6.0_r8
    CONS4=GAMMA(4.0_r8+BR)/6.0_r8
    CONS5=GAMMA(1.0_r8+BS)
    CONS6=GAMMA(1.0_r8+BR)
    CONS7=GAMMA(4.0_r8+BG)/6.0_r8
    CONS8=GAMMA(1.0_r8+BG)
    CONS9=GAMMA(5.0_r8/2.0_r8+BR/2.0_r8)
    CONS10=GAMMA(5.0_r8/2.0_r8+BS/2.0_r8)
    CONS11=GAMMA(5.0_r8/2.0_r8+BG/2.0_r8)
    CONS12=GAMMA(1.0_r8+DI)*CI
    CONS13=GAMMA(BS+3.0_r8)*PI/4.0_r8*ECI
    CONS14=GAMMA(BG+3.0_r8)*PI/4.0_r8*ECI
    CONS15=-1108.0_r8*EII*PI**((1.0_r8-BS)/3.0_r8)*RHOSN**((-2.0_r8-BS)/3.0_r8)/(4.0_r8*720.0_r8)
    CONS16=GAMMA(BI+3.0_r8)*PI/4.0_r8*ECI
    CONS17=4.0_r8*2.0_r8*3.0_r8*RHOSU*PI*ECI*ECI*GAMMA(2.0_r8*BS+2.0_r8)/(8.0_r8*(RHOG-RHOSN))
    CONS18=RHOSN*RHOSN
    CONS19=RHOW*RHOW
    CONS20=20.0_r8*PI*PI*RHOW*BIMM
    CONS21=4.0_r8/(DCS*RHOI)
    CONS22=PI*RHOI*DCS**3/6.0_r8
    CONS23=PI/4.0_r8*EII*GAMMA(BS+3.0_r8)
    CONS24=PI/4.0_r8*ECR*GAMMA(BR+3.0_r8)
    CONS25=PI*PI/24.0_r8*RHOW*ECR*GAMMA(BR+6.0_r8)
    CONS26=PI/6.0_r8*RHOW
    CONS27=GAMMA(1.0_r8+BI)
    CONS28=GAMMA(4.0_r8+BI)/6.0_r8
    CONS29=4.0_r8/3.0_r8*PI*RHOW*(25.E-6_r8)**3
    CONS30=4.0_r8/3.0_r8*PI*RHOW
    CONS31=PI*PI*ECR*RHOSN
    CONS32=PI/2.0_r8*ECR
    CONS33=PI*PI*ECR*RHOG
    CONS34=5.0_r8/2.0_r8+BR/2.0_r8
    CONS35=5.0_r8/2.0_r8+BS/2.0_r8
    CONS36=5.0_r8/2.0_r8+BG/2.0_r8
    CONS37=4.0_r8*PI*1.38E-23_r8/(6.0_r8*PI*RIN)
    CONS38=PI*PI/3.0_r8*RHOW
    CONS39=PI*PI/36.0_r8*RHOW*BIMM
    CONS40=PI/6.0_r8*BIMM
    CONS41=PI*PI*ECR*RHOW

    !+---+-----------------------------------------------------------------+
    !..Set these variables needed for computing radar reflectivity.  These
    !.. get used within radar_init to create other variables used in the
    !.. radar module.

    xam_r = PI*RHOW/6.0_r8
    xbm_r = 3.0_r8
    xmu_r = 0.0_r8
    xam_s = CS
    xbm_s = DS
    xmu_s = 0.0_r8
    xam_g = CG
    xbm_g = DG
    xmu_g = 0.0_r8

    CALL radar_init()
    !+---+-----------------------------------------------------------------+
    EFFCS=reimin/2.0_r8
    EFFIS=reimin

  END SUBROUTINE Init_Micro_HugMorr
  


  SUBROUTINE RunMicro_HugMorr( &
       nCols       , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax        , &!INTEGER      , INTENT(IN   ) :: kMax 
       prsi        , &
       prsl        , &
       tc          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       QV          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       QC          , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       QR          , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       QI          , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS          , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG          , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI          , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NS          , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       NR          , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       NG          , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       NC          , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)   
       dTcdt       , &!
       dqvdt       , &!
       dqcdt       , &!
       dqrdt       , &!
       dqidt       , &!
       dqsdt       , &!
       dqgdt       , &!
       dnidt       , &!
       dnsdt       , &!
       dnrdt       , &!
       dNGdt       , &!
       dNCdt       , &!
       TKE         , &!REAL(KIND=r8), INTENT(IN   ) :: TKE (1:nCols, 1:kMax)   
       KZH         , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax)   
       DT_IN       , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       omega       , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       EFFCS       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       EFFIS       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       LSRAIN      , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       LSSNOW        )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)

    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: kMax 
    REAL(KIND=r8), INTENT(IN   ) :: prsi       (1:nCols,1:kMax+1)  
    REAL(KIND=r8), INTENT(IN   ) :: prsl       (1:nCols,1:kMax)    

    ! Temporary changed from INOUT to IN
    REAL(KIND=r8), INTENT(INOUT) :: Tc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qv(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qi(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qs(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qg(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ni(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ns(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: nr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: NG(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: NC(1:nCols, 1:kMax)

    REAL(KIND=r8), INTENT(INOUT) :: dTcdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqvdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqcdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqrdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqidt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqsdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dqgdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dnidt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dnsdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dnrdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dNGdt(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: dNCdt(1:nCols, 1:kMax)

    REAL(KIND=r8), INTENT(IN   ) :: tke(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: kzh(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: dt_in
    REAL(KIND=r8), INTENT(IN   ) :: omega (1:nCols, 1:kMax) ! omega (Pa/s)
    REAL(KIND=r8), INTENT(OUT  ) :: LSRAIN(1:nCols)
    REAL(KIND=r8), INTENT(OUT  ) :: LSSNOW(1:nCols)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
    ! note: effis not currently passed out of microphysics (no coupling of ice eff rad with radiation)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)

    REAL(KIND=r8) :: SR        (1:nCols)
    REAL(KIND=r8) :: refl_10cm (1:nCols, 1:kMax)

    ! add cumulus tendencies

    REAL(KIND=r8) :: qrcuten (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qscuten (1:nCols, 1:kMax)
    REAL(KIND=r8) :: qicuten (1:nCols, 1:kMax)
    REAL(KIND=r8) :: mu      (1:nCols)

    LOGICAL       :: diagflag
    INTEGER       :: do_radar_ref  ! GT added for reflectivity calcs
    LOGICAL       :: F_QNDROP      ! wrf-chem
    REAL(KIND=r8) :: qndrop(1:nCols, 1:kMax) ! hm added, wrf-chem 

    REAL(KIND=r8) :: flip_pint  (nCols,kMax+1)   ! Interface pressures  
    REAL(KIND=r8) :: flip_pmid  (nCols,kMax)! Midpoint pressures 
    REAL(KIND=r8) :: flip_t     (nCols,kMax)! temperature
    REAL(KIND=r8) :: flip_q     (nCols,kMax)! specific humidity
    REAL(KIND=r8) :: flip_pdel  (nCols,kMax)! layer thickness
    REAL(KIND=r8) :: flip_rpdel (nCols,kMax)! inverse of layer thickness
    REAL(KIND=r8) :: flip_lnpmid(nCols,kMax)! Log Midpoint pressures    
    REAL(KIND=r8) :: flip_lnpint(nCols,kMax+1)   ! Log interface pressures
    REAL(KIND=r8) :: flip_zi    (nCols,kMax+1)! Height above surface at interfaces 
    REAL(KIND=r8) :: flip_zm    (nCols,kMax)  ! Geopotential height at mid level

    REAL(KIND=r8) :: zi         (nCols,kMax+1)     ! Height above surface at interfaces
    REAL(KIND=r8) :: zm         (nCols,kMax)        ! Geopotential height at mid level
    REAL(KIND=r8) :: p          (nCols,kMax) ! pressure at all points, on u,v levels (Pa).
    REAL(KIND=r8) :: dz         (nCols,kMax)
    REAL(KIND=r8) :: RAINNC (1:nCols)
    REAL(KIND=r8) :: RAINNCV(1:nCols)
    REAL(KIND=r8) :: SNOW(1:nCols)
    REAL(KIND=r8) :: rho(nCols,kMax)     
    REAL(KIND=r8) :: w  (1:nCols, 1:kMax) !, tke, nctend, nnColsnd,kzh
    REAL(KIND=r8) :: Tc_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: qv_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: qc_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: qr_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: qi_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: qs_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: qg_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: ni_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: ns_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: nr_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: NG_mic(1:nCols, 1:kMax)
    REAL(KIND=r8) :: NC_mic(1:nCols, 1:kMax)

    INTEGER       :: I,K,kflip
    RAINNC =0.0_r8
    refl_10cm=0.0_r8
    RAINNCV=0.0_r8
    SNOW   =0.0_r8
    qrcuten=0.0_r8
    qscuten=0.0_r8
    qicuten=0.0_r8
    mu     =1.0_r8
    F_QNDROP=.FALSE.
    qndrop=0.0_r8
    EFFCS =0.0_r8
    EFFIS =0.0_r8
    diagflag=.TRUE.
    do_radar_ref=1
    DO i=1,nCols
       !flip_pint       (i,kMax+1) = gps(i)*si(1) ! gps --> Pa
       flip_pint       (i,kMax+1) = prsi(i,1)
    END DO
    DO k=kMax,1,-1
       kflip=kMax+2-k
       DO i=1,nCols
         ! flip_pint    (i,k)      = MAX(si(kflip)*gps(i) ,1.0e-12_r8)
          flip_pint    (i,k)      = MAX(prsi(i,kflip)    ,1.0e-12_r8)
       END DO
    END DO
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols


          dTcdt(i,k) =  0.0_r8
          dqvdt(i,k) =  0.0_r8
          dqcdt(i,k) =  0.0_r8
          dqrdt(i,k) =  0.0_r8
          dqidt(i,k) =  0.0_r8
          dqsdt(i,k) =  0.0_r8
          dqgdt(i,k) =  0.0_r8
          dnidt(i,k) =  0.0_r8
          dnsdt(i,k) =  0.0_r8
          dnrdt(i,k) =  0.0_r8
          dNGdt(i,k) =  0.0_r8
          dNCdt(i,k) =  0.0_r8

          Tc_mic(i,k) = Tc(i,k)
          qv_mic(i,k) = qv(i,k)
          qc_mic(i,k) = qc(i,k)
          qr_mic(i,k) = qr(i,k)
          qi_mic(i,k) = qi(i,k)
          qs_mic(i,k) = qs(i,k)
          qg_mic(i,k) = qg(i,k)
          ni_mic(i,k) = ni(i,k)
          ns_mic(i,k) = ns(i,k)
          nr_mic(i,k) = nr(i,k)
          NG_mic(i,k) = NG(i,k)
          NC_mic(i,k) = NC(i,k)

          flip_t   (i,kflip) =  TC (i,k)
          flip_q   (i,kflip) =  qv (i,k)
          !flip_pmid(i,kflip) =  sl(  k)*gps (i)
          flip_pmid(i,kflip) =  prsl(i,k)
       END DO
    END DO
    DO k=1,kMax
       DO i=1,nCols    
          flip_pdel    (i,k) = MAX(flip_pint(i,k+1) - flip_pint(i,k),1.0e-12_r8)
          flip_rpdel   (i,k) = 1.0_r8/MAX((flip_pint(i,k+1) - flip_pint(i,k)),1.0e-12_r8)
          flip_lnpmid  (i,k) = LOG(flip_pmid(i,k))
       END DO
    END DO
    DO k=1,kMax+1
       DO i=1,nCols
          flip_lnpint(i,k) =  LOG(flip_pint  (i,k))
       END DO
    END DO

    !
    !..delsig     k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(2)
    !             k=3/2----si,ric,rf,km,kh,b,l -----------
    !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(1)
    !             k=1/2----si ----------------------------

    ! Derive new temperature and geopotential fields

    CALL geopotential_t(                                 &
         flip_lnpint(1:nCols,1:kMax+1)   , flip_pint  (1:nCols,1:kMax+1)  , &
         flip_pmid  (1:nCols,1:kMax)     , flip_pdel  (1:nCols,1:kMax)   , flip_rpdel(1:nCols,1:kMax)   , &
         flip_t     (1:nCols,1:kMax)     , flip_q     (1:nCols,1:kMax)   , rair   , gravit , zvir   ,&
         flip_zi    (1:nCols,1:kMax+1)   , flip_zm    (1:nCols,1:kMax)   , nCols   ,nCols, kMax)
    DO i=1,nCols
       zi(i,1) = flip_zi    (i,kMax+1)
    END DO
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols
          zi (i,k+1) = flip_zi    (i,kflip)
          zm (i,k  ) = flip_zm    (i,kflip)
          p  (i,k  ) = flip_pmid  (i,kflip)
       END DO
    END DO
    DO k=1,kMax
       DO i=1,nCols
          dz (i,k  ) = MAX(zi(i,k+1)-zi(i,k),1.0e-12)
          !j/kg/kelvin
          !
          ! P = rho * R * T
          !
          !            P
          ! rho  = -------
          !          R * T
          !
          rho   (i,k) =  (p(i,k)/(r_d*tc_mic(i,k)))       ! density
          w     (i,k) = -omega(i,k)/(rho(i,k)*gravit) ! (Pa/s)  - (m/s)
       END DO
    END DO

    CALL MP_MORR_TWO_MOMENT( &
         nCols                        , &!INTEGER      , INTENT(IN   ) :: nCols
         kMax                         , &!INTEGER      , INTENT(IN   ) :: kMax 
         tc_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: Tc_mic (1:nCols, 1:kMax)
         QV_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: qv_mic (1:nCols, 1:kMax)
         QC_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: qc_mic (1:nCols, 1:kMax)
         QR_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: qr_mic (1:nCols, 1:kMax)
         QI_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: qi_mic (1:nCols, 1:kMax)
         QS_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: qs_mic (1:nCols, 1:kMax)
         QG_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: qg_mic (1:nCols, 1:kMax)
         NI_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: ni_mic (1:nCols, 1:kMax)
         NS_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: ns_mic (1:nCols, 1:kMax)
         NR_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: nr_mic (1:nCols, 1:kMax)
         NG_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: NG_mic (1:nCols, 1:kMax)   
         NC_mic      (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: NC_mic (1:nCols, 1:kMax)  
         TKE         (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: TKE (1:nCols, 1:kMax)  
         KZH         (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: KZH (1:nCols, 1:kMax)  
         PBL                          , &!INTEGER(KIND=r8), INTENT(IN   ) :: PBL
         P           (1:nCols, 1:kMax), &!AIR PRESSURE (PA)
         DT_IN                        , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
         DZ          (1:nCols, 1:kMax), &!hm
         W           (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: w  (1:nCols, 1:kMax) !, tke, nctend, nnColsnd,kzh
         RAINNC      (1:nCols)        , &!REAL(KIND=r8), INTENT(INOUT) :: RAINNC (1:nCols)
         RAINNCV     (1:nCols)        , &!REAL(KIND=r8), INTENT(INOUT) :: RAINNCV(1:nCols)
         SNOW        (1:nCols)        , &!REAL(KIND=r8), INTENT(INOUT) :: SNOW(1:nCols)
         SR          (1:nCols)        , &!REAL(KIND=r8), INTENT(INOUT) :: SR     (1:nCols)
         refl_10cm   (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: refl_10cm(1:nCols, 1:kMax)
         qrcuten     (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: qrcuten(1:nCols, 1:kMax)
         qscuten     (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: qscuten(1:nCols, 1:kMax)
         qicuten     (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: qicuten(1:nCols, 1:kMax)
         mu          (1:nCols)        , &!REAL(KIND=r8), INTENT(IN   ) :: mu     (1:nCols)  ! hm added
         diagflag                     , &!LOGICAL      , OPTIONAL, INTENT(IN) :: diagflag
         do_radar_ref                 , &!INTEGER      , OPTIONAL, INTENT(IN) :: do_radar_ref ! GT added for reflectivity calcs
         F_QNDROP                     , &!LOGICAL      , OPTIONAL, INTENT(IN) :: F_QNDROP  ! wrf-chem
         qndrop      (1:nCols, 1:kMax), &!REAL(KIND=r8), OPTIONAL, INTENT(INOUT):: qndrop(1:nCols, 1:kMax) ! hm added, wrf-chem 
         EFFCS       (1:nCols, 1:kMax), &
         EFFIS       (1:nCols, 1:kMax)  )


    DO k=1,kMax
       DO i=1,nCols


          dTcdt(i,k) =  (Tc_mic(i,k) - Tc(i,k))/DT_IN
          dqvdt(i,k) =  (qv_mic(i,k) - qv(i,k))/DT_IN
          dqcdt(i,k) =  (qc_mic(i,k) - qc(i,k))/DT_IN
          dqrdt(i,k) =  (qr_mic(i,k) - qr(i,k))/DT_IN
          dqidt(i,k) =  (qi_mic(i,k) - qi(i,k))/DT_IN
          dqsdt(i,k) =  (qs_mic(i,k) - qs(i,k))/DT_IN
          dqgdt(i,k) =  (qg_mic(i,k) - qg(i,k))/DT_IN
          dnidt(i,k) =  (ni_mic(i,k) - ni(i,k))/DT_IN
          dnsdt(i,k) =  (ns_mic(i,k) - ns(i,k))/DT_IN
          dnrdt(i,k) =  (nr_mic(i,k) - nr(i,k))/DT_IN
          dNGdt(i,k) =  (NG_mic(i,k) - NG(i,k))/DT_IN
          dNCdt(i,k) =  (NC_mic(i,k) - NC(i,k))/DT_IN


          Tc(i,k) = Tc_mic(i,k)
          qv(i,k) = qv_mic(i,k)
          qc(i,k) = qc_mic(i,k)
          qr(i,k) = qr_mic(i,k)
          qi(i,k) = qi_mic(i,k)
          qs(i,k) = qs_mic(i,k)
          qg(i,k) = qg_mic(i,k)
          ni(i,k) = ni_mic(i,k)
          ns(i,k) = ns_mic(i,k)
          nr(i,k) = nr_mic(i,k)
          NG(i,k) = NG_mic(i,k)
          NC(i,k) = NC_mic(i,k)
       END DO
    END DO

    LSRAIN(1:nCols)=0.5_r8*RAINNC(1:nCols)/1000.0_r8  !(mm)->m
    LSSNOW(1:nCols)=0.5_r8*SNOW  (1:nCols)/1000.0_r8  !(mm)->m
  END SUBROUTINE RunMicro_HugMorr



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
  ! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
  ! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE MORR_TWO_MOMENT_MICRO) 
  ! WHICH OPERATES ON 1D VERTICAL COLUMNS.
  ! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
  ! BACK TO DRIVER MODEL USING THIS INTERFACE.
  ! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK TO DRIVER MODEL.

  ! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND SLAVA TATARSKII (GEORGIA TECH).

  ! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MP_MORR_TWO_MOMENT( &
       nCols       , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax        , &!INTEGER      , INTENT(IN   ) :: kMax 
       tc          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       QV          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       QC          , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       QR          , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       QI          , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS          , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG          , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI          , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NS          , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       NR          , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       NG          , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       NC          , &!REAL(KIND=r8), INTENT(INOUT) :: nc (1:nCols, 1:kMax)
       TKE         , &!REAL(KIND=r8), INTENT(INOUT) :: TKE (1:nCols, 1:kMax)  
       KZH         , &!REAL(KIND=r8), INTENT(INOUT) :: KZH (1:nCols, 1:kMax)  
       PBL         , &!INTEGER(KIND=r8), INTENT(IN   ) :: PBL
       P           , &! AIR PRESSURE (PA)
       DT_IN       , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       DZ          , &!* !hm
       W           , &!REAL(KIND=r8), INTENT(IN   ) :: w  (1:nCols, 1:kMax) !, tke, nctend, nnColsnd,kzh
       RAINNC      , &!REAL(KIND=r8), INTENT(INOUT) :: RAINNC    (1:nCols)
       RAINNCV     , &!REAL(KIND=r8), INTENT(INOUT) :: RAINNCV   (1:nCols)
       SNOW        , &!REAL(KIND=r8), INTENT(INOUT) :: SNOW      (1:nCols)
       SR          , &!REAL(KIND=r8), INTENT(INOUT) :: SR        (1:nCols)
       refl_10cm   , &!REAL(KIND=r8), INTENT(INOUT) :: refl_10cm (1:nCols, 1:kMax)
       qrcuten     , &!REAL(KIND=r8), INTENT(IN   ) :: qrcuten   (1:nCols, 1:kMax)
       qscuten     , &!REAL(KIND=r8), INTENT(IN   ) :: qscuten(1:nCols, 1:kMax)
       qicuten     , &!REAL(KIND=r8), INTENT(IN   ) :: qicuten(1:nCols, 1:kMax)
       mu          , &!REAL(KIND=r8), INTENT(IN   ) :: mu     (1:nCols)  ! hm added
       diagflag    , &!LOGICAL      , OPTIONAL, INTENT(IN) :: diagflag
       do_radar_ref, &!INTEGER      , OPTIONAL, INTENT(IN) :: do_radar_ref ! GT added for reflectivity calcs
       F_QNDROP    , &!LOGICAL      , OPTIONAL, INTENT(IN) :: F_QNDROP  ! wrf-chem
       qndrop      , &!REAL(KIND=r8), OPTIONAL, INTENT(INOUT):: qndrop(1:nCols, 1:kMax) ! hm added, wrf-chem 
       EFFCS       , &!
       EFFIS         )!

    ! QV - water vapor mixing ratio (kg/kg)
    ! QC - cloud water mixing ratio (kg/kg)
    ! QR - rain water mixing ratio  (kg/kg)
    ! QI - cloud ice mixing ratio   (kg/kg)
    ! QS - snow mixing ratio        (kg/kg)
    ! QG - graupel mixing ratio     (KG/KG)
    ! NI - cloud ice number concentration (1/kg)
    ! NS - Snow Number concentration (1/kg)
    ! NC - Cloud droplet Number concentration (1/kg)
    ! NR - Rain Number concentration (1/kg)
    ! NG - Graupel number concentration (1/kg)
    ! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
    ! P - AIR PRESSURE (PA)
    ! gps- AIR PRESSURE (PA)
    ! W - VERTICAL AIR VELOCITY (M/S)
    ! Tc -  TEMPERATURE (K)
    ! DZ - difference in height over interface (m)
    ! DT_IN - model time step (sec)
    ! RAINNC - accumulated grid-scale precipitation (mm)
    ! RAINNCV - one time step grid scale precipitation (mm/time step)
    ! SR - one time step mass ratio of snow to total precip
    ! qrcuten, rain tendency from parameterized cumulus convection
    ! qscuten, snow tendency from parameterized cumulus convection
    ! qicuten, cloud ice tendency from parameterized cumulus convection

    ! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
    ! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
    ! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
    ! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
    ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
    ! note: effis not currently passed out of microphysics (no coupling of ice eff rad with radiation)
    ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! reflectivity currently not included!!!!
    ! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
    !................................
    ! GRID_CLOCK, GRID_ALARMS - parameters to limit radar reflectivity calculation only when needed
    ! otherwise radar reflectivity calculation every time step is too slow
    ! only needed for coupling with WRF, see code below for details
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
    ! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
    ! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
    ! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

    ! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

    ! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
    ! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
    ! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
    ! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
    ! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

    ! ADDITIONAL INPUT NEEDED BY MICRO
    ! ********NOTE: WVAR IS SHOULD BE USED IN DROPLET ACTIVATION
    ! FOR CASES WHEN UPDRAFT IS NOT RESOLVED, EITHER BECAUSE OF
    ! LOW MODEL RESOLUTION OR CLOUD TYPE
    ! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

    IMPLICIT NONE

    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: kMax 

    ! Temporary changed from INOUT to IN
    REAL(KIND=r8), INTENT(INOUT) :: Tc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qv(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qi(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qs(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: qg(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ni(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: ns(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: nr(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(INOUT) :: NG(1:nCols, 1:kMax)   
    REAL(KIND=r8), INTENT(INOUT) :: nc(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: TKE(1:nCols, 1:kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: KZH(1:nCols, 1:kMax)  
    INTEGER      , INTENT(IN   ) :: PBL
    REAL(KIND=r8), INTENT(IN   ) :: p          (nCols,kMax) ! pressure at all points, on u,v levels (Pa).

    !    REAL(KIND=r8), INTENT(IN   ) :: rho(1:nCols, 1:kMax)

    REAL(KIND=r8), INTENT(IN   ) :: dt_in
    REAL(KIND=r8), INTENT(IN   ) :: dz         (nCols,kMax)

    !    REAL(KIND=r8), INTENT(IN   ) :: ht( 1:nCols  )
    REAL(KIND=r8), INTENT(IN   ) :: w(1:nCols, 1:kMax) !, tke, nctend, nnColsnd,kzh
    REAL(KIND=r8), INTENT(INOUT) :: RAINNC (1:nCols)
    REAL(KIND=r8), INTENT(INOUT) :: RAINNCV(1:nCols)
    REAL(KIND=r8), INTENT(INOUT) :: SNOW   (1:nCols)
    REAL(KIND=r8), INTENT(INOUT) :: SR     (1:nCols)
    REAL(KIND=r8), INTENT(INOUT) :: refl_10cm(1:nCols, 1:kMax)
    ! add cumulus tendencies

    REAL(KIND=r8), INTENT(IN   ):: qrcuten(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ):: qscuten(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ):: qicuten(1:nCols, 1:kMax)
    REAL(KIND=r8), INTENT(IN   ):: mu     (1:nCols)

    LOGICAL      , OPTIONAL, INTENT(IN) :: diagflag
    INTEGER      , OPTIONAL, INTENT(IN) :: do_radar_ref
    LOGICAL      , OPTIONAL, INTENT(IN) :: F_QNDROP  ! wrf-chem

    !jdf                      qndrop ! hm added, wrf-chem
    REAL(KIND=r8), OPTIONAL,INTENT(INOUT):: qndrop(1:nCols, 1:kMax)

    !jdf  REAL(KIND=r8), DIMENSION(1:nCols, 1:kMax),INTENT(INOUT):: CSED3D, &
    REAL(KIND=r8), OPTIONAL,INTENT(INOUT):: EFFCS (1:nCols, 1:kMax)
    REAL(KIND=r8), OPTIONAL,INTENT(INOUT):: EFFIS (1:nCols, 1:kMax)
    !, effcs, effis








    ! LOCAL VARIABLES

    REAL(KIND=r8) :: effi(1:nCols, 1:kMax)
    REAL(KIND=r8) :: effs(1:nCols, 1:kMax)
    REAL(KIND=r8) :: effr(1:nCols, 1:kMax)
    REAL(KIND=r8) :: EFFG(1:nCols, 1:kMax)

    REAL(KIND=r8) :: WVAR(1:nCols, 1:kMax)
    REAL(KIND=r8) :: EFFC(1:nCols, 1:kMax)

    REAL(KIND=r8) :: QC_TEND1D(1:kMax)
    REAL(KIND=r8) :: QI_TEND1D(1:kMax)
    REAL(KIND=r8) :: QNI_TEND1D(1:kMax)
    REAL(KIND=r8) :: QR_TEND1D(1:kMax)
    REAL(KIND=r8) :: NC_TEND1D(1:kMax)
    REAL(KIND=r8) :: NI_TEND1D(1:kMax)
    REAL(KIND=r8) :: NS_TEND1D(1:kMax)
    REAL(KIND=r8) :: NR_TEND1D(1:kMax)
    REAL(KIND=r8) :: QC1D(1:kMax)
    REAL(KIND=r8) :: QI1D(1:kMax)
    REAL(KIND=r8) :: QR1D(1:kMax)
    REAL(KIND=r8) :: NC1D(1:kMax)
    REAL(KIND=r8) :: NI1D(1:kMax)
    REAL(KIND=r8) :: NS1D(1:kMax)
    REAL(KIND=r8) :: NR1D(1:kMax)
    REAL(KIND=r8) :: QS1D(1:kMax)
    REAL(KIND=r8) :: T_TEND1D(1:kMax)
    REAL(KIND=r8) :: QV_TEND1D(1:kMax)
    REAL(KIND=r8) :: T1D(1:kMax)
    REAL(KIND=r8) :: QV1D(1:kMax)
    REAL(KIND=r8) :: P1D(1:kMax)
    REAL(KIND=r8) :: W1D(1:kMax)
    REAL(KIND=r8) :: WVAR1D(1:kMax)
    REAL(KIND=r8) :: EFFC1D(1:kMax)
    REAL(KIND=r8) :: EFFI1D(1:kMax)
    REAL(KIND=r8) :: EFFS1D(1:kMax)
    REAL(KIND=r8) :: EFFR1D(1:kMax)
    REAL(KIND=r8) :: DZ1D(1:kMax)
    ! HM ADD GRAUPEL
    REAL(KIND=r8) :: QG_TEND1D(1:kMax)
    REAL(KIND=r8) :: NG_TEND1D(1:kMax)
    REAL(KIND=r8) :: QG1D(1:kMax)
    REAL(KIND=r8) :: NG1D(1:kMax)
    REAL(KIND=r8) :: EFFG1D(1:kMax)
    ! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
    REAL(KIND=r8) :: QGSTEN(1:kMax)
    REAL(KIND=r8) :: QRSTEN(1:kMax)
    REAL(KIND=r8) :: QISTEN(1:kMax)
    REAL(KIND=r8) :: QNISTEN(1:kMax)
    REAL(KIND=r8) :: QCSTEN(1:kMax)
    ! ADD CUMULUS TENDENCIES
    REAL(KIND=r8) :: QRCU1D(1:kMax)
    REAL(KIND=r8) :: QSCU1D(1:kMax)
    REAL(KIND=r8) :: QICU1D(1:kMax)

    LOGICAL :: flag_qndrop  ! wrf-chem
    INTEGER :: iinum ! wrf-chem

    ! wrf-chem
    !REAL(KIND=r8) :: nc1d(1:kMax)
    !REAL(KIND=r8) :: nc_tend1d(1:kMax)
    REAL(KIND=r8) :: C2PREC(1:kMax)
    !REAL(KIND=r8) :: CSED(1:kMax)
    REAL(KIND=r8) :: ISED(1:kMax)
    REAL(KIND=r8) :: SSED(1:kMax)
    REAL(KIND=r8) :: GSED(1:kMax)
    REAL(KIND=r8) :: RSED(1:kMax)    
    ! HM add reflectivity      
    REAL(KIND=r8) :: dBZ(1:kMax)

    REAL(KIND=r8) :: PRECPRT1D
    REAL(KIND=r8) :: SNOWRT1D

    INTEGER       :: I,K

    REAL(KIND=r8) :: DT

    ! below for wrf-chem
    iinum=0
    C2PREC = 0.0_r8
    ISED = 0.0_r8;    SSED = 0.0_r8
    GSED = 0.0_r8;    RSED = 0.0_r8    
    dBZ = 0.0_r8;    PRECPRT1D= 0.0_r8    ;SNOWRT1D= 0.0_r8    
    ! LOCAL VARIABLES

    effi = 0.0_r8;    effs = 0.0_r8;    effr = 0.0_r8
    EFFG = 0.0_r8;    WVAR = 0.0_r8;    EFFC = 0.0_r8
    QC_TEND1D = 0.0_r8;    QI_TEND1D = 0.0_r8;    QNI_TEND1D = 0.0_r8;    QR_TEND1D = 0.0_r8
    NI_TEND1D = 0.0_r8;    NS_TEND1D = 0.0_r8;    NR_TEND1D = 0.0_r8;    QC1D = 0.0_r8
    NC_TEND1D = 0.0_r8
    QI1D = 0.0_r8;    QR1D = 0.0_r8; NC1D = 0.0_r8;     NI1D = 0.0_r8;    NS1D = 0.0_r8
    NR1D = 0.0_r8;    QS1D = 0.0_r8;    T_TEND1D = 0.0_r8;    QV_TEND1D = 0.0_r8
    T1D = 0.0_r8;    QV1D = 0.0_r8;    P1D = 0.0_r8;    W1D = 0.0_r8
    WVAR1D = 0.0_r8;    EFFC1D = 0.0_r8;    EFFI1D = 0.0_r8;    EFFS1D = 0.0_r8
    EFFR1D = 0.0_r8;    DZ1D = 0.0_r8
    ! HM ADD GRAUPEL
    QG_TEND1D = 0.0_r8;    NG_TEND1D = 0.0_r8
    QG1D = 0.0_r8;    NG1D = 0.0_r8;    EFFG1D = 0.0_r8;
    ! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
    QGSTEN = 0.0_r8;    QRSTEN = 0.0_r8;    QISTEN = 0.0_r8;    QNISTEN = 0.0_r8;    QCSTEN = 0.0_r8
    ! ADD CUMULUS TENDENCIES
    QRCU1D = 0.0_r8;    QSCU1D = 0.0_r8;    QICU1D = 0.0_r8

    flag_qndrop = .FALSE.
    IF ( PRESENT ( f_qndrop ) ) flag_qndrop = f_qndrop
!!!!!!!!!!!!!!!!!!!!!!

    ! Initialize tendencies (all set to 0) and transfer
    ! array to local variables
    DT = DT_IN   



    DO I=1,nCols
       DO K=1,kMax
          !T(I,K,J)        = temp(i,k,j)!TH(i,k,j)*PII(i,k,j)

          ! wvar is the ST. DEV. OF sub-grid vertical velocity, used for calculating droplet 
          ! activation rates.
          ! WVAR CAN BE DERIVED EITHER FROM PREDICTED TKE (AS IN MYJ PBL SCHEME),
          ! OR FROM EDDY DIFFUSION COEFFICIENT KZH (AS IN YSU PBL SCHEME),
          ! DEPENDING ON THE PARTICULAR pbl SCHEME DRIVER MODEL IS COUPLED WITH
          ! NOTE: IF MODEL HAS HIGH ENOUGH RESOLUTION TO RESOLVE UPDRAFTS, WVAR MAY 
          ! NOT BE NEEDED 


          ! amy WVAR(I,K)     = 0.5

          ! amy 
          ! for YSU pbl scheme:
          ! coupling as in WRF2
          IF (PBL.NE.1) THEN
             ! for MYJ pbl scheme or 3D TKE:
             WVAR(I,K)     = (0.667*tke(i,k))**0.5
          ELSE
             ! for YSU pbl scheme:
             WVAR(I,K) = KZH(I,K)/20.
          END IF

          WVAR(I,K) = MAX(0.1_r8,WVAR(I,K))
          WVAR(I,K) = MIN(4.0_r8,WVAR(I,K))

          ! add tendency from pbl to droplet and cloud ice concentration
          ! NEEDED FOR WRF TEMPORARILY!!!!
          ! OTHER DRIVER MODELS MAY ADD TURBULENT DIFFUSION TENDENCY FOR
          ! SCALARS SOMEWHERE ELSE IN THE MODEL (I.E, NOT IN THE MICROPHYSICS)
          ! IN THIS CASE THESE 2 LINES BELOW MAY BE REMOVED
          !
          !! amy added in physics_addtendc
          !
          !       nc(i,k,j) = nc(i,k,j)+nctend(i,k,j)*dt
          !       ni(i,k,j) = ni(i,k,j)+nitend(i,k,j)*dt

       END DO
    END DO

    DO i=1,nCols      ! i loop (east-west)
       !       DO j=jts,jte      ! j loop (north-south)
       !
       ! Transfer 3D arrays into 1D for microphysical calculations
       !

       ! hm , initialize 1d tendency arrays to zero

       DO k=1,kMax   ! k loop (vertical)

          QC_TEND1D(k)  = 0.0_r8
          QI_TEND1D(k)  = 0.0_r8
          QNI_TEND1D(k) = 0.0_r8
          QR_TEND1D(k)  = 0.0_r8
          NC_TEND1D(k)  = 0.0_r8  ! amy
          NI_TEND1D(k)  = 0.0_r8
          NS_TEND1D(k)  = 0.0_r8
          NR_TEND1D(k)  = 0.0_r8
          T_TEND1D(k)   = 0.0_r8
          QV_TEND1D(k)  = 0.0_r8

          QC1D(k)       = QC(i,k)
          QI1D(k)       = QI(i,k)
          QS1D(k)       = QS(i,k)
          QR1D(k)       = QR(i,k)

          NI1D(k)       = NI(i,k)
          !amy added nc1d
          NC1D(k)       = NC(i,k)
          NS1D(k)       = NS(i,k)
          NR1D(k)       = NR(i,k)
          ! HM ADD GRAUPEL
          QG1D(K)       = QG(I,K)
          NG1D(K)       = NG(I,K)
          QG_TEND1D(K)  = 0.0_r8
          NG_TEND1D(K)  = 0.0_r8

          T1D(k)        = TC(i,k)
          QV1D(k)       = QV(i,k)
          P1D(k)        = P(i,k)
          DZ1D(k)       = DZ(i,k)
          W1D(k)        = W(i,k)
          WVAR1D(k)     = WVAR(i,k)
          ! add cumulus tendencies, decouple from mu
          qrcu1d(k)     = qrcuten(i,k)/mu(i)
          qscu1d(k)     = qscuten(i,k)/mu(i)
          qicu1d(k)     = qicuten(i,k)/mu(i)
       END DO  !jdf added this


  CALL MORR_TWO_MOMENT_MICRO(&
       QC_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: QC3DTEN  (KTS:KTE) ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
       QI_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: QI3DTEN  (KTS:KTE) ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
       QNI_TEND1D , &!REAL(KIND=r8), INTENT(INOUT) :: QNI3DTEN (KTS:KTE) ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
       QR_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: QR3DTEN  (KTS:KTE) ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
       NC_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: NC3DTEN  (KTS:KTE) ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S) amy add
       NI_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: NI3DTEN  (KTS:KTE) ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
       NS_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: NS3DTEN  (KTS:KTE) ! SNOW NUMBER CONCENTRATION (1/KG/S)
       NR_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: NR3DTEN  (KTS:KTE) ! RAIN NUMBER CONCENTRATION (1/KG/S)
       QC1D       , &!REAL(KIND=r8), INTENT(INOUT) :: QC3D(KTS:KTE) ! CLOUD WATER MIXING RATIO (KG/KG)
       QI1D       , &!REAL(KIND=r8), INTENT(INOUT) :: QI3D(KTS:KTE) ! CLOUD ICE MIXING RATIO (KG/KG)
       QS1D       , &!REAL(KIND=r8), INTENT(INOUT) :: QNI3D(KTS:KTE) ! SNOW MIXING RATIO (KG/KG)
       QR1D       , &!REAL(KIND=r8), INTENT(INOUT) :: QR3D(KTS:KTE) ! RAIN MIXING RATIO (KG/KG)
       NI1D       , &!REAL(KIND=r8), INTENT(INOUT) :: NI3D(KTS:KTE) ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
       NS1D       , &!REAL(KIND=r8), INTENT(INOUT) :: NS3D(KTS:KTE) ! SNOW NUMBER CONCENTRATION (1/KG)
       NR1D       , &!REAL(KIND=r8), INTENT(INOUT) :: NR3D(KTS:KTE) ! RAIN NUMBER CONCENTRATION (1/KG)
       NC1D       , &!REAL(KIND=r8), INTENT(INOUT) :: NC3D(KTS:KTE) ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
       T_TEND1D   , &!REAL(KIND=r8), INTENT(INOUT) :: T3DTEN   (KTS:KTE) ! TEMPERATURE TENDENCY (K/S)
       QV_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: QV3DTEN  (KTS:KTE) ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
       T1D        , &!REAL(KIND=r8), INTENT(INOUT) :: T3D      (KTS:KTE) ! TEMPERATURE (K)
       QV1D       , &!REAL(KIND=r8), INTENT(INOUT) :: QV3D     (KTS:KTE) ! WATER VAPOR MIXING RATIO (KG/KG)
       P1D        , &!REAL(KIND=r8), INTENT(IN   ) :: PRES     (KTS:KTE) ! ATMOSPHERIC PRESSURE (PA)
       DZ1D       , &!REAL(KIND=r8), INTENT(IN   ) :: DZQ      (KTS:KTE) ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
       W1D        , &!REAL(KIND=r8), INTENT(IN   ) :: W3D      (KTS:KTE) ! GRID-SCALE VERTICAL VELOCITY (M/S)
       WVAR1D     , &!REAL(KIND=r8), INTENT(IN   ) :: WVAR     (KTS:KTE) ! SUB-GRID VERTICAL VELOCITY (M/S)
       PRECPRT1D  , &!REAL(KIND=r8), INTENT(OUT  ) :: PRECRT ! TOTAL PRECIP PER TIME STEP (mm)
       SNOWRT1D   , &!REAL(KIND=r8), INTENT(OUT  ) :: SNOWRT ! SNOW PER TIME STEP (mm)
       EFFC1D     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFC     (KTS:KTE) ! DROPLET EFFECTIVE RADIUS (MICRON)
       EFFI1D     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFI     (KTS:KTE) ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
       EFFS1D     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFS     (KTS:KTE) ! SNOW EFFECTIVE RADIUS (MICRON)
       EFFR1D     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFR     (KTS:KTE) ! RAIN EFFECTIVE RADIUS (MICRON)
       DT         , &!REAL(KIND=r8), INTENT(IN   ) :: DT     ! MODEL TIME STEP (SEC)
       1          , &!INTEGER	   , INTENT( IN)  :: KTS
       kMax       , &!INTEGER	   , INTENT( IN)  :: KTE
       QG_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: QG3DTEN  (KTS:KTE) ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
       NG_TEND1D  , &!REAL(KIND=r8), INTENT(INOUT) :: NG3DTEN  (KTS:KTE) ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
       QG1D       , &!REAL(KIND=r8), INTENT(INOUT) :: QG3D  (KTS:KTE) ! GRAUPEL MIX RATIO (KG/KG)
       NG1D       , &!REAL(KIND=r8), INTENT(INOUT) :: NG3D  (KTS:KTE) ! GRAUPEL NUMBER CONC (1/KG)
       EFFG1D     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFG  (KTS:KTE) ! GRAUPEL EFFECTIVE RADIUS (MICRON)
       qrcu1d     , &!REAL(KIND=r8), INTENT(IN   ) :: qrcu1d   (KTS:KTE)
       qscu1d     , &!REAL(KIND=r8), INTENT(IN   ) :: qscu1d   (KTS:KTE)
       qicu1d     , &!REAL(KIND=r8), INTENT(IN   ) :: qicu1d   (KTS:KTE)
       QGSTEN     , &!REAL(KIND=r8), INTENT(OUT  ) :: QGSTEN   (KTS:KTE) ! GRAUPEL SED TEND (KG/KG/S)
       QRSTEN     , &!REAL(KIND=r8), INTENT(OUT  ) :: QRSTEN   (KTS:KTE) ! RAIN SED TEND (KG/KG/S)
       QISTEN     , &!REAL(KIND=r8), INTENT(OUT  ) :: QISTEN   (KTS:KTE) ! CLOUD ICE SED TEND (KG/KG/S)
       QNISTEN    , &!REAL(KIND=r8), INTENT(OUT  ) :: QNISTEN  (KTS:KTE) ! SNOW SED TEND (KG/KG/S)
       QCSTEN       )!REAL(KIND=r8), INTENT(OUT  ) :: QCSTEN   (KTS:KTE) ! CLOUD WAT SED TEND (KG/KG/S)   


       !
       ! Transfer 1D arrays back into 3D arrays
       !
       DO k=1,kMax

          ! hm, add tendencies to update global variables 
          ! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
          ! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

          QC(i,k)        = QC1D(k)
          QI(i,k)        = QI1D(k)
          QS(i,k)        = QS1D(k)
          QR(i,k)        = QR1D(k)
          NI(i,k)        = NI1D(k)
          NS(i,k)        = NS1D(k)          
          NR(i,k)        = NR1D(k)
          !amy added nc
          NC(i,k)        = NC1D(k)
          QG(I,K)        = QG1D(K)
          NG(I,K)        = NG1D(K)

          TC(i,k)         = T1D(k)
          !temp(I,K)      = T(i,k)!!T(i,k)/PII(i,k) ! CONVERT TEMP BACK TO POTENTIAL TEMP
          QV(i,k)        = QV1D(k)

          EFFC(i,k)      = EFFC1D(k)
          EFFI(i,k)      = EFFI1D(k)
          EFFS(i,k)      = EFFS1D(k)
          EFFR(i,k)      = EFFR1D(k)
          EFFG(I,K)      = EFFG1D(K)

          ! wrf-chem
          IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
             qndrop(i,k) = nc1d(k)
             !jdf         CSED3D(I,K,J) = CSED(K)
          END IF
          ! EFFECTIVE RADIUS FOR RADIATION CODE (currently not coupled)
          ! amy added back in to cam shortwave
          ! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
          EFFCS(I,K)     = MIN(EFFC (I,K),50.0_r8)
          EFFCS(I,K)     = MAX(EFFCS(I,K),1.0_r8)
          EFFIS(I,K)     = MIN(EFFI (I,K),130.)
          EFFIS(I,K)     = MAX(EFFIS(I,K),13.)

       END DO

       ! hm modified so that m2005 precip variables correctly match wrf precip variables
       RAINNC (i) = RAINNC(I)+ PRECPRT1D
       RAINNCV(i) = PRECPRT1D
       SNOW   (i) = SNOWRT1D
       SR(i) = SNOWRT1D/(PRECPRT1D+1.E-12_r8)
       !+---+-----------------------------------------------------------------+
       IF ( PRESENT (diagflag) ) THEN
          IF (diagflag .AND. do_radar_ref == 1) THEN
             !CALL refl10cm_hm (qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d,   &
             !     t1d, p1d, dBZ,  kMax)
             DO k = 1, kMax
                refl_10cm(i,k) = MAX(-35.0_r8, dBZ(k))
             ENDDO
          ENDIF
       ENDIF
       !+---+-----------------------------------------------------------------+

       !       END DO
    END DO

  END SUBROUTINE MP_MORR_TWO_MOMENT


  !+---+-----------------------------------------------------------------+

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
  ! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
  ! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE MORR_TWO_MOMENT_MICRO) 
  ! WHICH OPERATES ON 1D VERTICAL COLUMNS.
  ! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
  ! BACK TO DRIVER MODEL USING THIS INTERFACE.
  ! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK TO DRIVER MODEL.

  ! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND SLAVA TATARSKII (GEORGIA TECH).

  ! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  SUBROUTINE MORR_TWO_MOMENT_MICRO(&
       QC3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: QC3DTEN  (KTS:KTE) ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
       QI3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: QI3DTEN  (KTS:KTE) ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
       QNI3DTEN , &!REAL(KIND=r8), INTENT(INOUT) :: QNI3DTEN (KTS:KTE) ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
       QR3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: QR3DTEN  (KTS:KTE) ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
       NC3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: NC3DTEN  (KTS:KTE) ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S) amy add
       NI3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: NI3DTEN  (KTS:KTE) ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
       NS3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: NS3DTEN  (KTS:KTE) ! SNOW NUMBER CONCENTRATION (1/KG/S)
       NR3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: NR3DTEN  (KTS:KTE) ! RAIN NUMBER CONCENTRATION (1/KG/S)
       QC3D     , &!REAL(KIND=r8), INTENT(INOUT) :: QC3D     (KTS:KTE) ! CLOUD WATER MIXING RATIO (KG/KG)
       QI3D     , &!REAL(KIND=r8), INTENT(INOUT) :: QI3D     (KTS:KTE) ! CLOUD ICE MIXING RATIO (KG/KG)
       QNI3D    , &!REAL(KIND=r8), INTENT(INOUT) :: QNI3D    (KTS:KTE) ! SNOW MIXING RATIO (KG/KG)
       QR3D     , &!REAL(KIND=r8), INTENT(INOUT) :: QR3D     (KTS:KTE) ! RAIN MIXING RATIO (KG/KG)
       NI3D     , &!REAL(KIND=r8), INTENT(INOUT) :: NI3D     (KTS:KTE) ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
       NS3D     , &!REAL(KIND=r8), INTENT(INOUT) :: NS3D     (KTS:KTE) ! SNOW NUMBER CONCENTRATION (1/KG)
       NR3D     , &!REAL(KIND=r8), INTENT(INOUT) :: NR3D     (KTS:KTE) ! RAIN NUMBER CONCENTRATION (1/KG)
       NC3D     , &!REAL(KIND=r8), INTENT(INOUT) :: NC3D     (KTS:KTE) ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
       T3DTEN   , &!REAL(KIND=r8), INTENT(INOUT) :: T3DTEN   (KTS:KTE) ! TEMPERATURE TENDENCY (K/S)
       QV3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: QV3DTEN  (KTS:KTE) ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
       T3D      , &!REAL(KIND=r8), INTENT(INOUT) :: T3D      (KTS:KTE) ! TEMPERATURE (K)
       QV3D     , &!REAL(KIND=r8), INTENT(INOUT) :: QV3D     (KTS:KTE) ! WATER VAPOR MIXING RATIO (KG/KG)
       PRES     , &!REAL(KIND=r8), INTENT(IN   ) :: PRES     (KTS:KTE) ! ATMOSPHERIC PRESSURE (PA)
       DZQ      , &!REAL(KIND=r8), INTENT(IN   ) :: DZQ      (KTS:KTE) ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
       W3D      , &!REAL(KIND=r8), INTENT(IN   ) :: W3D      (KTS:KTE) ! GRID-SCALE VERTICAL VELOCITY (M/S)
       WVAR     , &!REAL(KIND=r8), INTENT(IN   ) :: WVAR     (KTS:KTE) ! SUB-GRID VERTICAL VELOCITY (M/S)
       PRECRT   , &!REAL(KIND=r8), INTENT(OUT  ) :: PRECRT ! TOTAL PRECIP PER TIME STEP (mm)
       SNOWRT   , &!REAL(KIND=r8), INTENT(OUT  ) :: SNOWRT ! SNOW PER TIME STEP (mm)
       EFFC     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFC     (KTS:KTE) ! DROPLET EFFECTIVE RADIUS (MICRON)
       EFFI     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFI     (KTS:KTE) ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
       EFFS     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFS     (KTS:KTE) ! SNOW EFFECTIVE RADIUS (MICRON)
       EFFR     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFR     (KTS:KTE) ! RAIN EFFECTIVE RADIUS (MICRON)
       DT       , &!REAL(KIND=r8), INTENT(IN   ) :: DT         ! MODEL TIME STEP (SEC)
       KTS      , &!INTEGER      , INTENT( IN)  :: KTS
       KTE      , &!INTEGER      , INTENT( IN)  :: KTE
       QG3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: QG3DTEN  (KTS:KTE) ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
       NG3DTEN  , &!REAL(KIND=r8), INTENT(INOUT) :: NG3DTEN  (KTS:KTE) ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
       QG3D     , &!REAL(KIND=r8), INTENT(INOUT) :: QG3D     (KTS:KTE) ! GRAUPEL MIX RATIO (KG/KG)
       NG3D     , &!REAL(KIND=r8), INTENT(INOUT) :: NG3D     (KTS:KTE) ! GRAUPEL NUMBER CONC (1/KG)
       EFFG     , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFG     (KTS:KTE) ! GRAUPEL EFFECTIVE RADIUS (MICRON)
       qrcu1d   , &!REAL(KIND=r8), INTENT(IN   ) :: qrcu1d   (KTS:KTE)
       qscu1d   , &!REAL(KIND=r8), INTENT(IN   ) :: qscu1d   (KTS:KTE)
       qicu1d   , &!REAL(KIND=r8), INTENT(IN   ) :: qicu1d   (KTS:KTE)
       QGSTEN   , &!REAL(KIND=r8), INTENT(OUT  ) :: QGSTEN   (KTS:KTE) ! GRAUPEL SED TEND (KG/KG/S)
       QRSTEN   , &!REAL(KIND=r8), INTENT(OUT  ) :: QRSTEN   (KTS:KTE) ! RAIN SED TEND (KG/KG/S)
       QISTEN   , &!REAL(KIND=r8), INTENT(OUT  ) :: QISTEN   (KTS:KTE) ! CLOUD ICE SED TEND (KG/KG/S)
       QNISTEN  , &!REAL(KIND=r8), INTENT(OUT  ) :: QNISTEN  (KTS:KTE) ! SNOW SED TEND (KG/KG/S)
       QCSTEN     )!REAL(KIND=r8), INTENT(OUT  ) :: QCSTEN   (KTS:KTE) ! CLOUD WAT SED TEND (KG/KG/S)	

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
    ! MORRISON ET AL. 2005 JAS; MORRISON AND PINTO 2005 JAS.
    ! ADDITIONAL CHANGES ARE DESCRIBED IN DETAIL BY MORRISON, THOMPSON, TATARSKII (MWR, SUBMITTED)

    ! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
    ! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
    ! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.

    ! CODE STRUCTURE: MAIN SUBROUTINE IS 'MORR_TWO_MOMENT'. ALSO INCLUDED IN THIS FILE IS
    ! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
    ! 'FUNCTION GAMMA'.

    ! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! DECLARATIONS

    IMPLICIT NONE

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
    ! DEFINE ARRAY SIZES

    ! INPUT NUMBER OF GRID CELLS

    ! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
    INTEGER, INTENT( IN)  :: KTS,KTE

    REAL(KIND=r8), INTENT(INOUT) :: QC3DTEN  (KTS:KTE) ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: QI3DTEN  (KTS:KTE) ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: QNI3DTEN (KTS:KTE) ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: QR3DTEN  (KTS:KTE) ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: NC3DTEN  (KTS:KTE) ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S) amy add
    REAL(KIND=r8), INTENT(INOUT) :: NI3DTEN  (KTS:KTE) ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: NS3DTEN  (KTS:KTE) ! SNOW NUMBER CONCENTRATION (1/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: NR3DTEN  (KTS:KTE) ! RAIN NUMBER CONCENTRATION (1/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: QC3D     (KTS:KTE) ! CLOUD WATER MIXING RATIO (KG/KG)
    REAL(KIND=r8), INTENT(INOUT) :: QI3D     (KTS:KTE) ! CLOUD ICE MIXING RATIO (KG/KG)
    REAL(KIND=r8), INTENT(INOUT) :: QNI3D    (KTS:KTE) ! SNOW MIXING RATIO (KG/KG)
    REAL(KIND=r8), INTENT(INOUT) :: QR3D     (KTS:KTE) ! RAIN MIXING RATIO (KG/KG)
    REAL(KIND=r8), INTENT(INOUT) :: NC3D     (KTS:KTE) ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
    REAL(KIND=r8), INTENT(INOUT) :: NI3D     (KTS:KTE) ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
    REAL(KIND=r8), INTENT(INOUT) :: NS3D     (KTS:KTE) ! SNOW NUMBER CONCENTRATION (1/KG)
    REAL(KIND=r8), INTENT(INOUT) :: NR3D     (KTS:KTE) ! RAIN NUMBER CONCENTRATION (1/KG)
    REAL(KIND=r8), INTENT(INOUT) :: T3DTEN   (KTS:KTE) ! TEMPERATURE TENDENCY (K/S)
    REAL(KIND=r8), INTENT(INOUT) :: QV3DTEN  (KTS:KTE) ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: T3D      (KTS:KTE) ! TEMPERATURE (K)
    REAL(KIND=r8), INTENT(INOUT) :: QV3D     (KTS:KTE) ! WATER VAPOR MIXING RATIO (KG/KG)
    REAL(KIND=r8), INTENT(IN   ) :: PRES     (KTS:KTE) ! ATMOSPHERIC PRESSURE (PA)
    REAL(KIND=r8), INTENT(IN   ) :: DZQ      (KTS:KTE) ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
    REAL(KIND=r8), INTENT(IN   ) :: W3D      (KTS:KTE) ! GRID-SCALE VERTICAL VELOCITY (M/S)
    REAL(KIND=r8), INTENT(IN   ) :: WVAR     (KTS:KTE) ! SUB-GRID VERTICAL VELOCITY (M/S)

    ! HM ADDED GRAUPEL VARIABLES
    REAL(KIND=r8), INTENT(INOUT) :: QG3DTEN  (KTS:KTE) ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: NG3DTEN  (KTS:KTE) ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
    REAL(KIND=r8), INTENT(INOUT) :: QG3D     (KTS:KTE) ! GRAUPEL MIX RATIO (KG/KG)
    REAL(KIND=r8), INTENT(INOUT) :: NG3D     (KTS:KTE) ! GRAUPEL NUMBER CONC (1/KG)

    ! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

    REAL(KIND=r8), INTENT(OUT  ) :: QGSTEN   (KTS:KTE) ! GRAUPEL SED TEND (KG/KG/S)
    REAL(KIND=r8), INTENT(OUT  ) :: QRSTEN   (KTS:KTE) ! RAIN SED TEND (KG/KG/S)
    REAL(KIND=r8), INTENT(OUT  ) :: QISTEN   (KTS:KTE) ! CLOUD ICE SED TEND (KG/KG/S)
    REAL(KIND=r8), INTENT(OUT  ) :: QNISTEN  (KTS:KTE) ! SNOW SED TEND (KG/KG/S)
    REAL(KIND=r8), INTENT(OUT  ) :: QCSTEN   (KTS:KTE) ! CLOUD WAT SED TEND (KG/KG/S)        

    ! hm add cumulus tendencies for precip
    REAL(KIND=r8), INTENT(IN   ) :: qrcu1d   (KTS:KTE)
    REAL(KIND=r8), INTENT(IN   ) :: qscu1d   (KTS:KTE)
    REAL(KIND=r8), INTENT(IN   ) :: qicu1d   (KTS:KTE)

    ! OUTPUT VARIABLES

    REAL(KIND=r8), INTENT(OUT  ) :: PRECRT ! TOTAL PRECIP PER TIME STEP (mm)
    REAL(KIND=r8), INTENT(OUT  ) :: SNOWRT ! SNOW PER TIME STEP (mm)

    REAL(KIND=r8), INTENT(OUT  ) :: EFFC     (KTS:KTE) ! DROPLET EFFECTIVE RADIUS (MICRON)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFI     (KTS:KTE) ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFS     (KTS:KTE) ! SNOW EFFECTIVE RADIUS (MICRON)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFR     (KTS:KTE) ! RAIN EFFECTIVE RADIUS (MICRON)
    REAL(KIND=r8), INTENT(OUT  ) :: EFFG     (KTS:KTE) ! GRAUPEL EFFECTIVE RADIUS (MICRON)

    ! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

    REAL(KIND=r8), INTENT(IN   ) :: DT         ! MODEL TIME STEP (SEC)

    !.....................................................................................................
    ! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
    ! REST OF THE MODEL.

    ! SIZE PARAMETER VARIABLES

    REAL(KIND=r8) :: LAMC  (KTS:KTE) ! SLOPE PARAMETER FOR DROPLETS (M-1)
    REAL(KIND=r8) :: LAMI  (KTS:KTE) ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
    REAL(KIND=r8) :: LAMS  (KTS:KTE) ! SLOPE PARAMETER FOR SNOW (M-1)
    REAL(KIND=r8) :: LAMR  (KTS:KTE) ! SLOPE PARAMETER FOR RAIN (M-1)
    REAL(KIND=r8) :: LAMG  (KTS:KTE) ! SLOPE PARAMETER FOR GRAUPEL (M-1)
    REAL(KIND=r8) :: CDIST1(KTS:KTE) ! PSD PARAMETER FOR DROPLETS
    REAL(KIND=r8) :: N0I   (KTS:KTE) ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
    REAL(KIND=r8) :: N0S   (KTS:KTE) ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
    REAL(KIND=r8) :: N0RR  (KTS:KTE) ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
    REAL(KIND=r8) :: N0G   (KTS:KTE) ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
    REAL(KIND=r8) :: PGAM  (KTS:KTE) ! SPECTRAL SHAPE PARAMETER FOR DROPLETS

    ! MICROPHYSICAL PROCESSES

    REAL(KIND=r8) :: NSUBC   (KTS:KTE)! LOSS OF NC DURING EVAP
    REAL(KIND=r8) :: NSUBI   (KTS:KTE)! LOSS OF NI DURING SUB.
    REAL(KIND=r8) :: NSUBS   (KTS:KTE)! LOSS OF NS DURING SUB.
    REAL(KIND=r8) :: NSUBR   (KTS:KTE)! LOSS OF NR DURING EVAP
    REAL(KIND=r8) :: PRD     (KTS:KTE)! DEP CLOUD ICE
    REAL(KIND=r8) :: PRE     (KTS:KTE)! EVAP OF RAIN
    REAL(KIND=r8) :: PRDS    (KTS:KTE)! DEP SNOW
    REAL(KIND=r8) :: NNUCCC  (KTS:KTE)! CHANGE N DUE TO CONTACT FREEZ DROPLETS
    REAL(KIND=r8) :: MNUCCC  (KTS:KTE)! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
    REAL(KIND=r8) :: PRA     (KTS:KTE)! ACCRETION DROPLETS BY RAIN
    REAL(KIND=r8) :: PRC     (KTS:KTE)! AUTOCONVERSION DROPLETS
    REAL(KIND=r8) :: PCC     (KTS:KTE)! COND/EVAP DROPLETS
    REAL(KIND=r8) :: NNUCCD  (KTS:KTE)! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
    REAL(KIND=r8) :: MNUCCD  (KTS:KTE)! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
    REAL(KIND=r8) :: MNUCCR  (KTS:KTE)! CHANGE Q DUE TO CONTACT FREEZ RAIN
    REAL(KIND=r8) :: NNUCCR  (KTS:KTE)! CHANGE N DUE TO CONTACT FREEZ RAIN
    REAL(KIND=r8) :: NPRA    (KTS:KTE)! CHANGE IN N DUE TO DROPLET ACC BY RAIN
    REAL(KIND=r8) :: NRAGG   (KTS:KTE)! SELF-COLLECTION/BREAKUP OF RAIN
    REAL(KIND=r8) :: NSAGG   (KTS:KTE)! SELF-COLLECTION OF SNOW
    REAL(KIND=r8) :: NPRC    (KTS:KTE)! CHANGE NC AUTOCONVERSION DROPLETS
    REAL(KIND=r8) :: NPRC1   (KTS:KTE)! CHANGE NR AUTOCONVERSION DROPLETS
    REAL(KIND=r8) :: PRAI    (KTS:KTE)! CHANGE Q ACCRETION CLOUD ICE BY SNOW
    REAL(KIND=r8) :: PRCI    (KTS:KTE)! CHANGE Q AUTOCONVERSION CLOUD ICE
    REAL(KIND=r8) :: PSACWS  (KTS:KTE)! CHANGE Q DROPLET ACCRETION BY SNOW
    REAL(KIND=r8) :: NPSACWS (KTS:KTE)! CHANGE N DROPLET ACCRETION BY SNOW
    REAL(KIND=r8) :: PSACWI  (KTS:KTE)! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
    REAL(KIND=r8) :: NPSACWI (KTS:KTE)! CHANGE N DROPLET ACCRETION BY CLOUD ICE
    REAL(KIND=r8) :: NPRCI   (KTS:KTE)! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
    REAL(KIND=r8) :: NPRAI   (KTS:KTE)! CHANGE N ACCRETION CLOUD ICE
    REAL(KIND=r8) :: NMULTS  (KTS:KTE)! ICE MULT DUE TO RIMING DROPLETS BY SNOW
    REAL(KIND=r8) :: NMULTR  (KTS:KTE)! ICE MULT DUE TO RIMING RAIN BY SNOW
    REAL(KIND=r8) :: QMULTS  (KTS:KTE)! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
    REAL(KIND=r8) :: QMULTR  (KTS:KTE)! CHANGE Q DUE TO ICE RAIN/SNOW
    REAL(KIND=r8) :: PRACS   (KTS:KTE)! CHANGE Q RAIN-SNOW COLLECTION
    REAL(KIND=r8) :: NPRACS  (KTS:KTE)! CHANGE N RAIN-SNOW COLLECTION
    !REAL(KIND=r8) :: PCCN    (KTS:KTE)! CHANGE Q DROPLET ACTIVATION
    REAL(KIND=r8) :: PSMLT   (KTS:KTE)! CHANGE Q MELTING SNOW TO RAIN
    REAL(KIND=r8) :: EVPMS   (KTS:KTE)! CHNAGE Q MELTING SNOW EVAPORATING
    REAL(KIND=r8) :: NSMLTS  (KTS:KTE)! CHANGE N MELTING SNOW
    REAL(KIND=r8) :: NSMLTR  (KTS:KTE)! CHANGE N MELTING SNOW TO RAIN
    ! HM ADDED 12/13/06
    REAL(KIND=r8) :: PIACR   (KTS:KTE)! CHANGE QR, ICE-RAIN COLLECTION
    REAL(KIND=r8) :: NIACR   (KTS:KTE)! CHANGE N, ICE-RAIN COLLECTION
    REAL(KIND=r8) :: PRACI   (KTS:KTE)! CHANGE QI, ICE-RAIN COLLECTION
    REAL(KIND=r8) :: PIACRS  (KTS:KTE)! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
    REAL(KIND=r8) :: NIACRS  (KTS:KTE)! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
    REAL(KIND=r8) :: PRACIS  (KTS:KTE)! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
    REAL(KIND=r8) :: EPRD    (KTS:KTE)! SUBLIMATION CLOUD ICE
    REAL(KIND=r8) :: EPRDS   (KTS:KTE)! SUBLIMATION SNOW
    ! HM ADDED GRAUPEL PROCESSES
    REAL(KIND=r8) :: PRACG   (KTS:KTE)! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
    REAL(KIND=r8) :: PSACWG  (KTS:KTE)! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
    REAL(KIND=r8) :: PGSACW  (KTS:KTE)! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    REAL(KIND=r8) :: PGRACS  (KTS:KTE)! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    REAL(KIND=r8) :: PRDG    (KTS:KTE)! DEP OF GRAUPEL
    REAL(KIND=r8) :: EPRDG   (KTS:KTE)! SUB OF GRAUPEL
    REAL(KIND=r8) :: EVPMG   (KTS:KTE)! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
    REAL(KIND=r8) :: PGMLT   (KTS:KTE)! CHANGE Q MELTING OF GRAUPEL
    REAL(KIND=r8) :: NPRACG  (KTS:KTE)! CHANGE N COLLECTION RAIN BY GRAUPEL
    REAL(KIND=r8) :: NPSACWG (KTS:KTE)! CHANGE N COLLECTION DROPLETS BY GRAUPEL
    REAL(KIND=r8) :: NSCNG   (KTS:KTE)! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    REAL(KIND=r8) :: NGRACS  (KTS:KTE)! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    REAL(KIND=r8) :: NGMLTG  (KTS:KTE)! CHANGE N MELTING GRAUPEL
    REAL(KIND=r8) :: NGMLTR  (KTS:KTE)! CHANGE N MELTING GRAUPEL TO RAIN
    REAL(KIND=r8) :: NSUBG   (KTS:KTE)! CHANGE N SUB/DEP OF GRAUPEL
    REAL(KIND=r8) :: PSACR   (KTS:KTE)! CONVERSION DUE TO COLL OF SNOW BY RAIN
    REAL(KIND=r8) :: NMULTG  (KTS:KTE)! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
    REAL(KIND=r8) :: NMULTRG (KTS:KTE)! ICE MULT DUE TO ACC RAIN BY GRAUPEL
    REAL(KIND=r8) :: QMULTG  (KTS:KTE)! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
    REAL(KIND=r8) :: QMULTRG (KTS:KTE)! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL

    ! TIME-VARYING ATMOSPHERIC PARAMETERS

    REAL(KIND=r8) :: KAP     (KTS:KTE)! THERMAL CONDUCTIVITY OF AIR
    REAL(KIND=r8) :: EVS     (KTS:KTE)! SATURATION VAPOR PRESSURE
    REAL(KIND=r8) :: EIS     (KTS:KTE)! ICE SATURATION VAPOR PRESSURE
    REAL(KIND=r8) :: QVS     (KTS:KTE)! SATURATION MIXING RATIO
    REAL(KIND=r8) :: QVI     (KTS:KTE)! ICE SATURATION MIXING RATIO
    REAL(KIND=r8) :: QVQVS   (KTS:KTE)! SATURATION RATIO
    REAL(KIND=r8) :: QVQVSI  (KTS:KTE)! ICE SATURAION RATIO
    REAL(KIND=r8) :: DV      (KTS:KTE)! DIFFUSIVITY OF WATER VAPOR IN AIR
    REAL(KIND=r8) :: XXLS    (KTS:KTE)! LATENT HEAT OF SUBLIMATION
    REAL(KIND=r8) :: XXLV    (KTS:KTE)! LATENT HEAT OF VAPORIZATION
    REAL(KIND=r8) :: CPM     (KTS:KTE)! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
    REAL(KIND=r8) :: MU      (KTS:KTE)! VISCOCITY OF AIR
    REAL(KIND=r8) :: SC      (KTS:KTE)! SCHMIDT NUMBER
    REAL(KIND=r8) :: XLF     (KTS:KTE)! LATENT HEAT OF FREEZING
    REAL(KIND=r8) :: RHO     (KTS:KTE)! AIR DENSITY
    REAL(KIND=r8) :: AB      (KTS:KTE)! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
    REAL(KIND=r8) :: ABI     (KTS:KTE)! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

    ! TIME-VARYING MICROPHYSICS PARAMETERS

    REAL(KIND=r8) :: DAP   (KTS:KTE)! DIFFUSIVITY OF AEROSOL
    REAL(KIND=r8) :: NACNT          ! NUMBER OF CONTACT IN
    REAL(KIND=r8) :: FMULT          ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
    !REAL(KIND=r8) :: COFFI          ! ICE AUTOCONVERSION PARAMETER

    ! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

    REAL(KIND=r8) :: DUMI    (KTS:KTE)
    REAL(KIND=r8) :: DUMR    (KTS:KTE)
    REAL(KIND=r8) :: DUMFNI  (KTS:KTE)
    REAL(KIND=r8) :: DUMG    (KTS:KTE)
    REAL(KIND=r8) :: DUMFNG  (KTS:KTE)
    REAL(KIND=r8) :: UNI
    REAL(KIND=r8) :: UMI
    REAL(KIND=r8) :: UMR
    REAL(KIND=r8) :: FR      (KTS:KTE)
    REAL(KIND=r8) :: FI      (KTS:KTE)
    REAL(KIND=r8) :: FNI     (KTS:KTE)
    REAL(KIND=r8) :: FG      (KTS:KTE)
    REAL(KIND=r8) :: FNG     (KTS:KTE)
    REAL(KIND=r8) :: RGVM
    REAL(KIND=r8) :: FALOUTR (KTS:KTE)
    REAL(KIND=r8) :: FALOUTI (KTS:KTE)
    REAL(KIND=r8) :: FALOUTNI(KTS:KTE)
    REAL(KIND=r8) :: FALTNDR
    REAL(KIND=r8) :: FALTNDI
    REAL(KIND=r8) :: FALTNDNI
    !REAL(KIND=r8) :: RHO2
    REAL(KIND=r8) :: DUMQS   (KTS:KTE)
    REAL(KIND=r8) :: DUMFNS  (KTS:KTE)
    REAL(KIND=r8) :: UMS
    REAL(KIND=r8) :: UNS
    REAL(KIND=r8) :: FS      (KTS:KTE)
    REAL(KIND=r8) :: FNS     (KTS:KTE)
    REAL(KIND=r8) :: FALOUTS (KTS:KTE)
    REAL(KIND=r8) :: FALOUTNS(KTS:KTE)
    REAL(KIND=r8) :: FALOUTG (KTS:KTE)
    REAL(KIND=r8) :: FALOUTNG(KTS:KTE)
    REAL(KIND=r8) :: FALTNDS
    REAL(KIND=r8) :: FALTNDNS
    REAL(KIND=r8) :: UNR
    REAL(KIND=r8) :: FALTNDG
    REAL(KIND=r8) :: FALTNDNG
    REAL(KIND=r8) :: DUMC    (KTS:KTE)
    REAL(KIND=r8) :: DUMFNC  (KTS:KTE)
    REAL(KIND=r8) :: UNC
    REAL(KIND=r8) :: UMC
    REAL(KIND=r8) :: UNG
    REAL(KIND=r8) :: UMG
    REAL(KIND=r8) :: FC      (KTS:KTE)
    REAL(KIND=r8) :: FALOUTC (KTS:KTE)
    REAL(KIND=r8) :: FALOUTNC(KTS:KTE)
    REAL(KIND=r8) :: FALTNDC
    REAL(KIND=r8) :: FALTNDNC
    REAL(KIND=r8) :: FNC     (KTS:KTE)
    REAL(KIND=r8) :: DUMFNR  (KTS:KTE)
    REAL(KIND=r8) :: FALOUTNR(KTS:KTE)
    REAL(KIND=r8) :: FALTNDNR
    REAL(KIND=r8) :: FNR     (KTS:KTE)

    ! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

    REAL(KIND=r8) :: AIN     (KTS:KTE)
    REAL(KIND=r8) :: ARN     (KTS:KTE)
    REAL(KIND=r8) :: ASN     (KTS:KTE)
    REAL(KIND=r8) :: ACN     (KTS:KTE)
    REAL(KIND=r8) :: AGN     (KTS:KTE)

    ! EXTERNAL FUNCTION CALL RETURN VARIABLES

    !      REAL GAMMA,      ! EULER GAMMA FUNCTION
    !      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
    !      REAL DERF1        ! ERROR FUNCTION

    ! DUMMY VARIABLES

    REAL(KIND=r8) :: DUM
    REAL(KIND=r8) :: DUM1
    REAL(KIND=r8) :: DUM2
    REAL(KIND=r8) :: DUMT
    REAL(KIND=r8) :: DUMQV
    REAL(KIND=r8) :: DUMQSS
    !REAL(KIND=r8) :: DUMQSI
    REAL(KIND=r8) :: DUMS

    ! PROGNOSTIC SUPERSATURATION

    REAL(KIND=r8) :: DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
    REAL(KIND=r8) :: DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
    REAL(KIND=r8) :: EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
    REAL(KIND=r8) :: EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
    REAL(KIND=r8) :: EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
    REAL(KIND=r8) :: EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

    ! NEW DROPLET ACTIVATION VARIABLES
    REAL(KIND=r8) :: TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
    REAL(KIND=r8) :: TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
    REAL(KIND=r8) :: TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
    REAL(KIND=r8) :: TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
    REAL(KIND=r8) :: TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
    REAL(KIND=r8) :: DUMACT
    REAL(KIND=r8) :: DUM3

    ! COUNTING/INDEX VARIABLES

    INTEGER :: K
    INTEGER :: NSTEP
    INTEGER :: N ! ,I

    ! LTRUE IS ONLY USED TO SPEED UP THE CODE !!
    ! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN, 
    !               = 1, HYDROMETEORS IN COLUMN

    INTEGER :: LTRUE

    ! DROPLET ACTIVATION/FREEZING AEROSOL


    !REAL(KIND=r8) :: CT      ! DROPLET ACTIVATION PARAMETER
    !REAL(KIND=r8) :: TEMP1   ! DUMMY TEMPERATURE
    !REAL(KIND=r8) :: SAT1    ! DUMMY SATURATION
    REAL(KIND=r8) :: SIGVL   ! SURFACE TENSION LIQ/VAPOR
    !REAL(KIND=r8) :: KEL     ! KELVIN PARAMETER
    REAL(KIND=r8) :: KC2     ! TOTAL ICE NUCLEATION RATE

    !REAL(KIND=r8) :: CRY   ! AEROSOL ACTIVATION PARAMETERS
    !REAL(KIND=r8) :: KRY   ! AEROSOL ACTIVATION PARAMETERS

    ! MORE WORKING/DUMMY VARIABLES

    !REAL(KIND=r8) :: DUMQI
    !REAL(KIND=r8) :: DUMNI
    !REAL(KIND=r8) :: DC0
    !REAL(KIND=r8) :: DS0
    !REAL(KIND=r8) :: DG0
    REAL(KIND=r8) :: DUMQC
    !REAL(KIND=r8) :: DUMQR
    REAL(KIND=r8) :: RATIO
    REAL(KIND=r8) :: SUM_DEP
    REAL(KIND=r8) :: FUDGEF

    ! EFFECTIVE VERTICAL VELOCITY  (M/S)
    !REAL(KIND=r8) :: WEF

    ! WORKING PARAMETERS FOR ICE NUCLEATION

!    REAL(KIND=r8) :: ANUC
!    REAL(KIND=r8) :: BNUC

    ! WORKING PARAMETERS FOR AEROSOL ACTIVATION

    REAL(KIND=r8) :: AACT
    REAL(KIND=r8) :: GAMM
    REAL(KIND=r8) :: GG
    REAL(KIND=r8) :: PSI
    REAL(KIND=r8) :: ETA1
    REAL(KIND=r8) :: ETA2
    REAL(KIND=r8) :: SM1
    REAL(KIND=r8) :: SM2
    REAL(KIND=r8) :: SMAX
    REAL(KIND=r8) :: UU1
    REAL(KIND=r8) :: UU2
    REAL(KIND=r8) :: ALPHA

    ! DUMMY SIZE DISTRIBUTION PARAMETERS

    REAL(KIND=r8) :: DLAMS
    REAL(KIND=r8) :: DLAMR
    REAL(KIND=r8) :: DLAMI
    REAL(KIND=r8) :: DLAMC
    REAL(KIND=r8) :: DLAMG
    REAL(KIND=r8) :: LAMMAX
    REAL(KIND=r8) :: LAMMIN

    INTEGER :: IDROP

    ! DROPLET CONCENTRATION AND ITS TENDENCY
    ! NOTE: CURRENTLY DROPLET CONCENTRATION IS SPECIFIED !!!!!
    ! TENDENCY OF NC IS CALCULATED BUT IT IS NOT USED !!!
    ! amy       REAL, DIMENSION(KTS:KTE) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! SET LTRUE INITIALLY TO 0

    LTRUE = 0

    ! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
    DO K = KTS,KTE

       ! 
       !.....................................................................
       ! ZERO OUT PROCESS RATES
       ! hm 3/4/13 --  note that this initialization is redundant
       ! this was added by Amy likely for output of process rates, so I'll leave
       ! this in for now. this will result in a small (likely trivial) slowdown of the code.

       PRC(K) = 0.0_r8
       NPRC(K) = 0.0_r8
       NPRC1(K) = 0.0_r8
       PRA(K) = 0.0_r8
       NPRA(K) = 0.0_r8
       NRAGG(K) = 0.0_r8
       NSMLTS(K) = 0.0_r8
       NSMLTR(K) = 0.0_r8
       EVPMS(K) = 0.0_r8
       PCC(K) = 0.0_r8
       PRE(K) = 0.0_r8
       NSUBC(K) = 0.0_r8
       NSUBR(K) = 0.0_r8
       PRACG(K) = 0.0_r8
       NPRACG(K) = 0.0_r8
       PSMLT(K) = 0.0_r8
       PGMLT(K) = 0.0_r8
       EVPMG(K) = 0.0_r8
       PRACS(K) = 0.0_r8
       NPRACS(K) = 0.0_r8
       NGMLTG(K) = 0.0_r8
       NGMLTR(K) = 0.0_r8


       MNUCCC(K) = 0.0_r8
       NNUCCC(K) = 0.0_r8
       PRC(K) = 0.0_r8
       NPRC(K) = 0.0_r8
       NPRC1(K) = 0.0_r8
       NSAGG(K) = 0.0_r8
       PSACWS(K) = 0.0_r8
       NPSACWS(K) = 0.0_r8
       PSACWI(K) = 0.0_r8
       NPSACWI(K) = 0.0_r8
       PRACS(K) = 0.0_r8
       NPRACS(K) = 0.0_r8
       NMULTS(K) = 0.0_r8
       QMULTS(K) = 0.0_r8
       NMULTR(K) = 0.0_r8
       QMULTR(K) = 0.0_r8
       NMULTG(K) = 0.0_r8
       QMULTG(K) = 0.0_r8
       NMULTRG(K) = 0.0_r8
       QMULTRG(K) = 0.0_r8
       MNUCCR(K) = 0.0_r8
       NNUCCR(K) = 0.0_r8
       PRA(K) = 0.0_r8
       NPRA(K) = 0.0_r8
       NRAGG(K) = 0.0_r8
       PRCI(K) = 0.0_r8
       NPRCI(K) = 0.0_r8
       PRAI(K) = 0.0_r8
       NPRAI(K) = 0.0_r8
       NNUCCD(K) = 0.0_r8
       MNUCCD(K) = 0.0_r8
       PCC(K) = 0.0_r8
       PRE(K) = 0.0_r8
       PRD(K) = 0.0_r8
       PRDS(K) = 0.0_r8
       EPRD(K) = 0.0_r8
       EPRDS(K) = 0.0_r8
       NSUBC(K) = 0.0_r8
       NSUBI(K) = 0.0_r8
       NSUBS(K) = 0.0_r8
       NSUBR(K) = 0.0_r8
       PIACR(K) = 0.0_r8
       NIACR(K) = 0.0_r8
       PRACI(K) = 0.0_r8
       PIACRS(K) = 0.0_r8
       NIACRS(K) = 0.0_r8
       PRACIS(K) = 0.0_r8
       ! HM: ADD GRAUPEL PROCESSES
       PRACG(K) = 0.0_r8
       PSACR(K) = 0.0_r8
       PSACWG(K) = 0.0_r8
       PGSACW(K) = 0.0_r8
       PGRACS(K) = 0.0_r8
       PRDG(K) = 0.0_r8
       EPRDG(K) = 0.0_r8
       NPRACG(K) = 0.0_r8
       NPSACWG(K) = 0.0_r8
       NSCNG(K) = 0.0_r8
       NGRACS(K) = 0.0_r8
       NSUBG(K) = 0.0_r8

       ! NC3DTEN LOCAL ARRAY INITIALIZED
       ! amy               NC3DTEN(K) = 0.

       ! LATENT HEAT OF VAPORATION

       XXLV(K) = 3.1484E6_r8 -2370.0_r8*T3D(K)

       ! LATENT HEAT OF SUBLIMATION

       XXLS(K) = 3.15E6_r8-2370.0_r8*T3D(K)+0.3337E6_r8

       CPM(K) = CP*(1.0_r8+0.887_r8*QV3D(K))

       ! SATURATION VAPOR PRESSURE AND MIXING RATIO

       EVS(K) = POLYSVP(T3D(K),0)   ! PA
       EIS(K) = POLYSVP(T3D(K),1)   ! PA

       ! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

       IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

       QVS(K) = 0.622_r8*EVS(K)/(PRES(K)-EVS(K))
       QVI(K) = 0.622_r8*EIS(K)/(PRES(K)-EIS(K))

       QVQVS(K) = QV3D(K)/QVS(K)
       QVQVSI(K) = QV3D(K)/QVI(K)

       ! AIR DENSITY

       RHO(K) = PRES(K)/(r_d*T3D(K))

       ! ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
       ! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
       ! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
       ! FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON

       IF (QRCU1D(K).GE.1.E-10_r8) THEN
          DUM=1.8e5_r8*(QRCU1D(K)*DT/(PI*RHOW*RHO(K)**3))**0.25_r8
          NR3D(K)=NR3D(K)+DUM
       END IF
       IF (QSCU1D(K).GE.1.E-10_r8) THEN
          DUM=3.e5_r8*(QSCU1D(K)*DT/(CONS1*RHO(K)**3))**(1.0_r8/(DS+1.0_r8))
          NS3D(K)=NS3D(K)+DUM
       END IF
       IF (QICU1D(K).GE.1.E-10_r8) THEN
          DUM=QICU1D(K)*DT/(CI*(80.E-6_r8)**DI)
          NI3D(K)=NI3D(K)+DUM
       END IF

       ! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
       ! hm modify 7/0/09 change limit to 1.e-8

       IF (QVQVS(K).LT.0.9_r8) THEN
          IF (QR3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QR3D(K)
             T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
             QR3D(K)=0.0_r8
          END IF
          IF (QC3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QC3D(K)
             T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
             QC3D(K)=0.0_r8
          END IF
       END IF

       IF (QVQVSI(K).LT.0.9_r8) THEN
          IF (QI3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QI3D(K)
             T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
             QI3D(K)=0.0_r8
          END IF
          IF (QNI3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QNI3D(K)
             T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
             QNI3D(K)=0.0_r8
          END IF
          IF (QG3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QG3D(K)
             T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
             QG3D(K)=0.0_r8
          END IF
       END IF

       ! HEAT OF FUSION

       XLF(K) = XXLS(K)-XXLV(K)

       !..................................................................
       ! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
          QC3D(K) = 0.0_r8
          NC3D(K) = 0.0_r8
          EFFC(K) = 0.0_r8
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
          QR3D(K) = 0.0_r8
          NR3D(K) = 0.0_r8
          EFFR(K) = 0.0_r8
       END IF
       IF (QI3D(K).LT.QSMALL) THEN
          QI3D(K) = 0.0_r8
          NI3D(K) = 0.0_r8
          EFFI(K) = 0.0_r8
       END IF
       IF (QNI3D(K).LT.QSMALL) THEN
          QNI3D(K) = 0.0_r8
          NS3D(K) = 0.0_r8
          EFFS(K) = 0.0_r8
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
          QG3D(K) = 0.0_r8
          NG3D(K) = 0.0_r8
          EFFG(K) = 0.0_r8
       END IF

       ! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

       QRSTEN(K) = 0.0_r8
       QISTEN(K) = 0.0_r8
       QNISTEN(K) = 0.0_r8
       QCSTEN(K) = 0.0_r8
       QGSTEN(K) = 0.0_r8

       !..................................................................
       ! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

       ! DYNAMIC VISCOSITY OF AIR

       MU(K) = 1.496E-6_r8*T3D(K)**1.5_r8/(T3D(K)+120.0_r8)

       ! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

       DUM = (RHOSU/RHO(K))**0.54_r8

       !            AIN(K) = DUM*AI
       ! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction 
       ! hm fix 11/18/11
       AIN(K) = (RHOSU/RHO(K))**0.35_r8*AI
       ARN(K) = DUM*AR
       ASN(K) = DUM*AS
       !            ACN(K) = DUM*AC
       ! AA revision 4/1/11: temperature-dependent Stokes fall speed
       ACN(K) = SHR_CONST_G*RHOW/(18.0_r8*MU(K))
       ! HM ADD GRAUPEL 8/28/06
       AGN(K) = DUM*AG

       !hm 4/7/09 bug fix, initialize lami to prevent later division by zero
       LAMI(K)=0.0_r8

       !..................................
       ! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
       ! FOR THIS LEVEL

       IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
            .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) THEN
          IF (T3D(K).LT.273.15_r8.AND.QVQVSI(K).LT.0.999_r8) GOTO 200
          IF (T3D(K).GE.273.15_r8.AND.QVQVS(K).LT.0.999_r8) GOTO 200
       END IF

       ! THERMAL CONDUCTIVITY FOR AIR

       KAP(K) = 1.414E3_r8*MU(K)

       ! DIFFUSIVITY OF WATER VAPOR

       DV(K) = 8.794E-5_r8*T3D(K)**1.81_r8/PRES(K)

       ! SCHMIT NUMBER

       SC(K) = MU(K)/(RHO(K)*DV(K))

       ! PSYCHOMETIC CORRECTIONS

       ! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

       DUM = (r_v*T3D(K)**2)

       DQSDT = XXLV(K)*QVS(K)/DUM
       DQSIDT =  XXLS(K)*QVI(K)/DUM

       ABI(K) = 1.0_r8+DQSIDT*XXLS(K)/CPM(K)
       AB(K) = 1.0_r8+DQSDT*XXLV(K)/CPM(K)

       !.....................................................................
       !.....................................................................
       ! CASE FOR TEMPERATURE ABOVE FREEZING

       IF (T3D(K).GE.273.15_r8) THEN

          !......................................................................
          !HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
          ! INUM = 0, PREDICT DROPLET NUMBER
          ! INUM = 1, SET CONSTANT DROPLET NUMBER

          !!amy 
          IF (INUM.EQ.1) THEN
             ! CONVERT NDCNST FROM CM-3 TO KG-1
             NC3D(K)=NDCNST*1.E6_r8/RHO(K)
          END IF

          ! GET SIZE DISTRIBUTION PARAMETERS

          ! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
          IF (QNI3D(K).LT.1.E-6_r8) THEN
             QR3D(K)=QR3D(K)+QNI3D(K)
             NR3D(K)=NR3D(K)+NS3D(K)
             T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
             QNI3D(K) = 0.0_r8
             NS3D(K) = 0.0_r8
          END IF
          IF (QG3D(K).LT.1.E-6_r8) THEN
             QR3D(K)=QR3D(K)+QG3D(K)
             NR3D(K)=NR3D(K)+NG3D(K)
             T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
             QG3D(K) = 0.0_r8
             NG3D(K) = 0.0_r8
          END IF

          IF (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8_r8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8_r8) GOTO 300

          ! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

          NS3D(K) = MAX(0.0_r8,NS3D(K))
          NC3D(K) = MAX(0.0_r8,NC3D(K))
          NR3D(K) = MAX(0.0_r8,NR3D(K))
          NG3D(K) = MAX(0.0_r8,NG3D(K))

          !......................................................................
          ! RAIN

          IF (QR3D(K).GE.QSMALL) THEN
             LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1.0_r8/3.0_r8)
             N0RR(K) = NR3D(K)*LAMR(K)

             ! CHECK FOR SLOPE

             ! ADJUST VARS

             IF (LAMR(K).LT.LAMMINR) THEN

                LAMR(K) = LAMMINR

                N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                NR3D(K) = N0RR(K)/LAMR(K)
             ELSE IF (LAMR(K).GT.LAMMAXR) THEN
                LAMR(K) = LAMMAXR
                N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                NR3D(K) = N0RR(K)/LAMR(K)
             END IF
          END IF

          !......................................................................
          ! CLOUD DROPLETS

          ! MARTIN ET AL. (1994) FORMULA FOR PGAM

          IF (QC3D(K).GE.QSMALL) THEN

             DUM = PRES(K)/(287.15_r8*T3D(K))
             PGAM(K)=0.0005714_r8*(NC3D(K)/1.E6_r8*DUM)+0.2714_r8
             PGAM(K)=1.0_r8/(PGAM(K)**2)-1.0_r8
             PGAM(K)=MAX(PGAM(K),2.0_r8)
             PGAM(K)=MIN(PGAM(K),10.0_r8)

             ! CALCULATE LAMC

             LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.0_r8)/   &
                  (QC3D(K)*GAMMA(PGAM(K)+1.0_r8)))**(1.0_r8/3.0_r8)

             ! LAMMIN, 60 MICRON DIAMETER
             ! LAMMAX, 1 MICRON

             LAMMIN = (PGAM(K)+1.0_r8)/60.E-6_r8
             LAMMAX = (PGAM(K)+1.0_r8)/1.E-6_r8

             IF (LAMC(K).LT.LAMMIN) THEN
                LAMC(K) = LAMMIN

                NC3D(K) = EXP(3.0_r8*LOG(LAMC(K))+LOG(QC3D(K))+              &
                     LOG(GAMMA(PGAM(K)+1.0_r8))-LOG(GAMMA(PGAM(K)+4.0_r8)))/CONS26
             ELSE IF (LAMC(K).GT.LAMMAX) THEN
                LAMC(K) = LAMMAX

                NC3D(K) = EXP(3.0_r8*LOG(LAMC(K))+LOG(QC3D(K))+              &
                     LOG(GAMMA(PGAM(K)+1.0_r8))-LOG(GAMMA(PGAM(K)+4.0_r8)))/CONS26

             END IF

          END IF

          !......................................................................
          ! SNOW

          IF (QNI3D(K).GE.QSMALL) THEN
             LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1.0_r8/DS)
             N0S(K) = NS3D(K)*LAMS(K)

             ! CHECK FOR SLOPE

             ! ADJUST VARS

             IF (LAMS(K).LT.LAMMINS) THEN
                LAMS(K) = LAMMINS
                N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                NS3D(K) = N0S(K)/LAMS(K)

             ELSE IF (LAMS(K).GT.LAMMAXS) THEN

                LAMS(K) = LAMMAXS
                N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                NS3D(K) = N0S(K)/LAMS(K)
             END IF
          END IF

          !......................................................................
          ! GRAUPEL

          IF (QG3D(K).GE.QSMALL) THEN
             LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1.0_r8/DG)
             N0G(K) = NG3D(K)*LAMG(K)

             ! ADJUST VARS

             IF (LAMG(K).LT.LAMMING) THEN
                LAMG(K) = LAMMING
                N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                NG3D(K) = N0G(K)/LAMG(K)

             ELSE IF (LAMG(K).GT.LAMMAXG) THEN

                LAMG(K) = LAMMAXG
                N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                NG3D(K) = N0G(K)/LAMG(K)
             END IF
          END IF

          !.....................................................................
          ! ZERO OUT PROCESS RATES

          PRC(K) = 0.0_r8
          NPRC(K) = 0.0_r8
          NPRC1(K) = 0.0_r8
          PRA(K) = 0.0_r8
          NPRA(K) = 0.0_r8
          NRAGG(K) = 0.0_r8
          NSMLTS(K) = 0.0_r8
          NSMLTR(K) = 0.0_r8
          EVPMS(K) = 0.0_r8
          PCC(K) = 0.0_r8
          PRE(K) = 0.0_r8
          NSUBC(K) = 0.0_r8
          NSUBR(K) = 0.0_r8
          PRACG(K) = 0.0_r8
          NPRACG(K) = 0.0_r8
          PSMLT(K) = 0.0_r8
          PGMLT(K) = 0.0_r8
          EVPMG(K) = 0.0_r8
          PRACS(K) = 0.0_r8
          NPRACS(K) = 0.0_r8
          NGMLTG(K) = 0.0_r8
          NGMLTR(K) = 0.0_r8

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

          !.................................................................
          !.......................................................................
          ! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
          ! FORMULA FROM BEHENG (1994)
          ! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
          ! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
          ! AS A GAMMA DISTRIBUTION

          ! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

          IF (QC3D(K).GE.1.E-6_r8) THEN

             ! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
             ! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

             PRC(K)=1350.0_r8*QC3D(K)**2.47_r8*  &
                  (NC3D(K)/1.e6_r8*RHO(K))**(-1.79_r8)

             ! note: nprc1 is change in Nr,
             ! nprc is change in Nc

             NPRC1(K) = PRC(K)/CONS29
             NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

             ! hm bug fix 3/4/13
             NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
             NPRC1(K) = MIN(NPRC1(K),NPRC(K))

          END IF

          !.......................................................................
          ! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
          ! FORMULA FROM IKAWA AND SAITO (1991)

          IF (QR3D(K).GE.1.E-8_r8.AND.QNI3D(K).GE.1.E-8_r8) THEN

             UMS = ASN(K)*CONS3/(LAMS(K)**BS)
             UMR = ARN(K)*CONS4/(LAMR(K)**BR)
             UNS = ASN(K)*CONS5/LAMS(K)**BS
             UNR = ARN(K)*CONS6/LAMR(K)**BR

             ! SET REASLISTIC LIMITS ON FALLSPEEDS

             ! bug fix, 10/08/09
             dum=(rhosu/rho(k))**0.54_r8
             UMS=MIN(UMS,1.2_r8*dum)
             UNS=MIN(UNS,1.2_r8*dum)
             UMR=MIN(UMR,9.1_r8*dum)
             UNR=MIN(UNR,9.1_r8*dum)

             ! hm fix, 3/4/13
             ! for above freezing conditions to get accelerated melting of snow,
             ! we need collection of rain by snow (following Lin et al. 1983)
             !            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
             !                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
             !                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
             !                  (5./(LAMS(K)**3*LAMR(K))+                    &
             !                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
             !                  0.5/(LAMS(K)*LAMR(K)**3)))

             PRACS(K) = CONS41*(((1.2_r8*UMR-0.95_r8*UMS)**2+                   &
                  0.08_r8*UMS*UMR)**0.5_r8*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5.0_r8/(LAMR(K)**3*LAMS(K))+                    &
                  2.0_r8/(LAMR(K)**2*LAMS(K)**2)+                  &                                 
                  0.5_r8/(LAMR(k)*LAMS(k)**3)))

             ! v3 5/27/11 npracs no longer used
             !            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
             !                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
             !                (1./(LAMR(K)**3*LAMS(K))+                      &
             !                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
             !                 1./(LAMR(K)*LAMS(K)**3))

          END IF

          ! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
          ! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
          ! ASSUME SHED DROPS ARE 1 MM IN SIZE

          IF (QR3D(K).GE.1.E-8_r8.AND.QG3D(K).GE.1.E-8_r8) THEN

             UMG = AGN(K)*CONS7/(LAMG(K)**BG)
             UMR = ARN(K)*CONS4/(LAMR(K)**BR)
             UNG = AGN(K)*CONS8/LAMG(K)**BG
             UNR = ARN(K)*CONS6/LAMR(K)**BR

             ! SET REASLISTIC LIMITS ON FALLSPEEDS
             ! bug fix, 10/08/09
             dum=(rhosu/rho(k))**0.54_r8
             UMG=MIN(UMG,20.0_r8*dum)
             UNG=MIN(UNG,20.0_r8*dum)
             UMR=MIN(UMR,9.1_r8*dum)
             UNR=MIN(UNR,9.1_r8*dum)

             ! DUM IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
             DUM = CONS41*(((1.2_r8*UMR-0.95_r8*UMG)**2+                   &
                  0.08_r8*UMG*UMR)**0.5_r8*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5.0_r8/(LAMR(K)**3*LAMG(K))+                    &
                  2.0_r8/(LAMR(K)**2*LAMG(K)**2)+                                   &
                  0.5_r8/(LAMR(k)*LAMG(k)**3)))

             ! ASSUME 1 MM DROPS ARE SHED, GET NUMBER SHED PER SEC

             DUM = DUM/5.2E-7_r8

             NPRACG(K) = CONS32*RHO(K)*(1.7_r8*(UNR-UNG)**2+            &
                  0.3_r8*UNR*UNG)**0.5_r8*N0RR(K)*N0G(K)*              &
                  (1.0_r8/(LAMR(K)**3*LAMG(K))+                      &
                  1.0_r8/(LAMR(K)**2*LAMG(K)**2)+                   &
                  1.0_r8/(LAMR(K)*LAMG(K)**3))

             ! hm 7/15/13, remove limit so that the number of collected drops can smaller than
             ! number of shed drops
             !            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
             NPRACG(K)=NPRACG(K)-DUM


          END IF

          !.......................................................................
          ! ACCRETION OF CLOUD LIQUID WATER BY RAIN
          ! CONTINUOUS COLLECTION EQUATION WITH
          ! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

          IF (QR3D(K).GE.1.E-8_r8 .AND. QC3D(K).GE.1.E-8_r8) THEN

             ! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
             ! KHAIROUTDINOV AND KOGAN 2000, MWR

             DUM=(QC3D(K)*QR3D(K))
             PRA(K) = 67.0_r8*(DUM)**1.15_r8
             NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

          END IF
          !.......................................................................
          ! SELF-COLLECTION OF RAIN DROPS
          ! FROM BEHENG(1994)
          ! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
          ! AS DESCRINED ABOVE FOR AUTOCONVERSION

          IF (QR3D(K).GE.1.E-8_r8) THEN
             ! include breakup add 10/09/09
             dum1=300.e-6_r8
             IF (1.0_r8/lamr(k).LT.dum1) THEN
                dum=1.0_r8
             ELSE IF (1.0_r8/lamr(k).GE.dum1) THEN
                dum=2.0_r8-EXP(2300.0_r8*(1.0_r8/lamr(k)-dum1))
             END IF
             !            NRAGG(K) = -8.0_r8*NR3D(K)*QR3D(K)*RHO(K)
             NRAGG(K) = -5.78_r8*dum*NR3D(K)*QR3D(K)*RHO(K)
          END IF

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)

          IF (QR3D(K).GE.QSMALL) THEN
             EPSR = 2.0_r8*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                  (F1R/(LAMR(K)*LAMR(K))+                       &
                  F2R*(ARN(K)*RHO(K)/MU(K))**0.5_r8*                      &
                  SC(K)**(1.0_r8/3.0_r8)*CONS9/                   &
                  (LAMR(K)**CONS34))
          ELSE
             EPSR = 0.0_r8
          END IF

          ! NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED

          IF (QV3D(K).LT.QVS(K)) THEN
             PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
             PRE(K) = MIN(PRE(K),0.0_r8)
          ELSE
             PRE(K) = 0.0_r8
          END IF

          !.......................................................................
          ! MELTING OF SNOW

          ! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
          ! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

          IF (QNI3D(K).GE.1.E-8_r8) THEN

             ! v3 5/27/11 bug fix
             !             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15_r8)*PRACS(K)

             ! hm fix 1/20/15
             !             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
             !                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
             !                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
             !                    SC(K)**(1./3.)*CONS10/                   &
             !                   (LAMS(K)**CONS35))+DUM
             PSMLT(K)=2.0_r8*PI*N0S(K)*KAP(K)*(273.15_r8-T3D(K))/       &
                  XLF(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                  F2S*(ASN(K)*RHO(K)/MU(K))**0.5_r8*                      &
                  SC(K)**(1.0_r8/3.0_r8)*CONS10/                   &
                  (LAMS(K)**CONS35))+DUM

             ! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

             IF (QVQVS(K).LT.1.0_r8) THEN
                EPSS = 2.0_r8*PI*N0S(K)*RHO(K)*DV(K)*                            &
                     (F1S/(LAMS(K)*LAMS(K))+                       &
                     F2S*(ASN(K)*RHO(K)/MU(K))**0.5_r8*                      &
                     SC(K)**(1.0_r8/3.0_r8)*CONS10/                   &
                     (LAMS(K)**CONS35))
                ! bug fix V1.4
                EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)    
                EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
                PSMLT(K) = PSMLT(K)-EVPMS(K)
             END IF
          END IF

          !.......................................................................
          ! MELTING OF GRAUPEL

          ! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
          ! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

          IF (QG3D(K).GE.1.E-8_r8) THEN

             ! v3 5/27/11 bug fix
             !             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15_r8)*PRACG(K)

             ! hm fix 1/20/15
             !             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/                  &
             !                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
             !                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
             !                    SC(K)**(1./3.)*CONS11/                   &
             !                   (LAMG(K)**CONS36))+DUM
             PGMLT(K)=2.0_r8*PI*N0G(K)*KAP(K)*(273.15_r8-T3D(K))/                  &
                  XLF(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                  F2S*(AGN(K)*RHO(K)/MU(K))**0.5_r8*                      &
                  SC(K)**(1.0_r8/3.0_r8)*CONS11/                   &
                  (LAMG(K)**CONS36))+DUM

             ! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

             IF (QVQVS(K).LT.1.0_r8) THEN
                EPSG = 2.0_r8*PI*N0G(K)*RHO(K)*DV(K)*                                &
                     (F1S/(LAMG(K)*LAMG(K))+                               &
                     F2S*(AGN(K)*RHO(K)/MU(K))**0.5_r8*                      &
                     SC(K)**(1.0_r8/3.0_r8)*CONS11/                   &
                     (LAMG(K)**CONS36))
                ! bug fix V1.4
                EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
                EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
                PGMLT(K) = PGMLT(K)-EVPMG(K)
             END IF
          END IF

          ! HM, V2.1
          ! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
          ! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
          ! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

          PRACG(K) = 0.0_r8
          PRACS(K) = 0.0_r8

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
          ! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
          ! CALCULATION

          ! CONSERVATION OF QC

          DUM = (PRC(K)+PRA(K))*DT

          IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN

             RATIO = QC3D(K)/DUM

             PRC(K) = PRC(K)*RATIO
             PRA(K) = PRA(K)*RATIO

          END IF

          ! CONSERVATION OF SNOW

          DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

          IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

             ! NO SOURCE TERMS FOR SNOW AT T > FREEZING
             RATIO = QNI3D(K)/DUM

             PSMLT(K) = PSMLT(K)*RATIO
             EVPMS(K) = EVPMS(K)*RATIO
             PRACS(K) = PRACS(K)*RATIO

          END IF

          ! CONSERVATION OF GRAUPEL

          DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT

          IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

             ! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
             RATIO = QG3D(K)/DUM

             PGMLT(K) = PGMLT(K)*RATIO
             EVPMG(K) = EVPMG(K)*RATIO
             PRACG(K) = PRACG(K)*RATIO

          END IF

          ! CONSERVATION OF QR
          ! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

          DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT

          IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN
             IF(PRE(K) == 0.0_r8)THEN
                RATIO = 0.0_r8
             ELSE
               RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                    (-PRE(K))
             END IF
             PRE(K) = PRE(K)*RATIO

          END IF

          !....................................

          QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

          T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
               (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)

          QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
          QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
          QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
          QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
          ! v3 5/27/11
          !      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
          !      NG3DTEN(K) = NG3DTEN(K)
          NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
          NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))

          IF (PRE(K).LT.0.0_r8) THEN
             DUM = PRE(K)*DT/QR3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSUBR(K) = DUM*NR3D(K)/DT
          END IF

          IF (EVPMS(K)+PSMLT(K).LT.0.0_r8) THEN
             DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSMLTS(K) = DUM*NS3D(K)/DT
          END IF
          IF (PSMLT(K).LT.0.0_r8) THEN
             DUM = PSMLT(K)*DT/QNI3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSMLTR(K) = DUM*NS3D(K)/DT
          END IF
          IF (EVPMG(K)+PGMLT(K).LT.0.0_r8) THEN
             DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NGMLTG(K) = DUM*NG3D(K)/DT
          END IF
          IF (PGMLT(K).LT.0.0_r8) THEN
             DUM = PGMLT(K)*DT/QG3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NGMLTR(K) = DUM*NG3D(K)/DT
          END IF

          NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
          NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
          NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))

300       CONTINUE

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
          ! WATER SATURATION

          DUMT = T3D(K)+DT*T3DTEN(K)
          DUMQV = QV3D(K)+DT*QV3DTEN(K)
          DUMQSS = 0.622_r8*POLYSVP(DUMT,0)/ (PRES(K)-POLYSVP(DUMT,0))
          DUMQC = QC3D(K)+DT*QC3DTEN(K)
          DUMQC = MAX(DUMQC,0.0_r8)

          ! SATURATION ADJUSTMENT FOR LIQUID

          DUMS = DUMQV-DUMQSS
          PCC(K) = DUMS/(1.0_r8+XXLV(K)**2*DUMQSS/(CPM(K)*r_v*DUMT**2))/DT
          IF (PCC(K)*DT+DUMQC.LT.0.0_r8) THEN
             PCC(K) = -DUMQC/DT
          END IF

          QV3DTEN(K) = QV3DTEN(K)-PCC(K)
          T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
          QC3DTEN(K) = QC3DTEN(K)+PCC(K)

          !.......................................................................
          ! ACTIVATION OF CLOUD DROPLETS
          ! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
          ! DROPLET CONCENTRATION IS SPECIFIED !!!!!
          !.......................................................................

          !!amy
          IF (INUM.EQ.0) THEN

             IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) THEN

                ! EFFECTIVE VERTICAL VELOCITY (M/S)

                IF (ISUB.EQ.0) THEN
                   ! ADD SUB-GRID VERTICAL VELOCITY
                   DUM = W3D(K)+WVAR(K)

                   ! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
                   DUM = MAX(DUM,0.10_r8)

                ELSE IF (ISUB.EQ.1) THEN
                   DUM=W3D(K)
                END IF

                ! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
                IF (DUM.GE.0.001_r8) THEN

                   IF (IBASE.EQ.1) THEN

                      ! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
                      ! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

                      IDROP=0

                      IF (QC3D(K)+QC3DTEN(K)*DT.LE.0.05E-3_r8/RHO(K)) THEN
                         IDROP=1
                      END IF
                      IF (K.EQ.1) THEN
                         IDROP=1
                      ELSE IF (K.GE.2) THEN
                         IF (QC3D(K)+QC3DTEN(K)*DT.GT.0.05E-3_r8/RHO(K).AND. &
                              QC3D(K-1)+QC3DTEN(K-1)*DT.LE.0.05E-3_r8/RHO(K-1)) THEN
                            IDROP=1
                         END IF
                      END IF

                      IF (IDROP.EQ.1) THEN
                         ! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

                         IF (IACT.EQ.1) THEN
                            ! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
                            ! BASED ON TWOMEY 1959

                            DUM=DUM*100.0_r8  ! CONVERT FROM M/S TO CM/S
                            DUM2 = 0.88_r8*C1**(2.0_r8/(K1+2.0_r8))*(7.E-2_r8*DUM**1.5_r8)**(K1/(K1+2.0_r8))
                            DUM2=DUM2*1.E6_r8 ! CONVERT FROM CM-3 TO M-3
                            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                            DUM2 = (DUM2-NC3D(K))/DT
                            DUM2 = MAX(0.0_r8,DUM2)
                            NC3DTEN(K) = NC3DTEN(K)+DUM2

                         ELSE IF (IACT.EQ.2) THEN
                            ! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

                            SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                            AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)
                            ALPHA = SHR_CONST_G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-SHR_CONST_G*MA/(RR*T3D(K))
                            GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

                            GG = 1.0_r8/(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                 (T3D(K)*RR)-1.0_r8))

                            PSI = 2.0_r8/3.0_r8*(ALPHA*DUM/GG)**0.5_r8*AACT

                            ETA1 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW1)
                            ETA2 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW2)

                            SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                            SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                            DUM1 = 1.0_r8/SM1**2*(F11*(PSI/ETA1)**1.5_r8+F21*(SM1**2/(ETA1+3.0_r8*PSI))**0.75_r8)
                            DUM2 = 1.0_r8/SM2**2*(F12*(PSI/ETA2)**1.5_r8+F22*(SM2**2/(ETA2+3.0_r8*PSI))**0.75_r8)

                            SMAX = 1.0_r8/(DUM1+DUM2)**0.5_r8

                            UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                            UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                            DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                            DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                            DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                            ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                            DUM2 = (DUM2-NC3D(K))/DT
                            DUM2 = MAX(0.0_r8,DUM2)
                            NC3DTEN(K) = NC3DTEN(K)+DUM2
                         END IF  ! IACT

                         !.............................................................................
                      ELSE IF (IDROP.EQ.0) THEN
                         ! ACTIVATE IN CLOUD INTERIOR
                         ! FIND EQUILIBRIUM SUPERSATURATION

                         TAUC=1.0_r8/(2.0_r8*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.0_r8)/LAMC(K))
                         IF (EPSR.GT.1.E-8_r8) THEN
                            TAUR=1.0_r8/EPSR
                         ELSE
                            TAUR=1.E8_r8
                         END IF

                         ! hm fix 1/20/15
                         !           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
                         DUM3=(-QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*SHR_CONST_G*DUM
                         DUM3=DUM3*TAUC*TAUR/(TAUC+TAUR)

                         IF (DUM3/QVS(K).GE.1.E-6_r8) THEN
                            IF (IACT.EQ.1) THEN

                               ! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQUILIBRIUM SS

                               DUM=DUM*100.0_r8  ! CONVERT FROM M/S TO CM/S
                               DUMACT = 0.88_r8*C1**(2.0_r8/(K1+2.0_r8))*(7.E-2_r8*DUM**1.5_r8)**(K1/(K1+2.0_r8))

                               ! USE POWER LAW CCN SPECTRA

                               ! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
                               DUM3=DUM3/QVS(K)*100.0_r8

                               DUM2=C1*DUM3**K1
                               ! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
                               DUM2=MIN(DUM2,DUMACT)
                               DUM2=DUM2*1.E6_r8 ! CONVERT FROM CM-3 TO M-3
                               DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                               DUM2 = (DUM2-NC3D(K))/DT
                               DUM2 = MAX(0.0_r8,DUM2)
                               NC3DTEN(K) = NC3DTEN(K)+DUM2

                            ELSE IF (IACT.EQ.2) THEN

                               ! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQUILIBRIUM SS

                               SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                               AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)
                               ALPHA = SHR_CONST_G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-SHR_CONST_G*MA/(RR*T3D(K))
                               GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

                               GG = 1.0_r8/(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                    (T3D(K)*RR)-1.0_r8))

                               PSI = 2.0_r8/3.0_r8*(ALPHA*DUM/GG)**0.5_r8*AACT

                               ETA1 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW1)
                               ETA2 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW2)

                               SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                               SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                               DUM1 = 1.0_r8/SM1**2*(F11*(PSI/ETA1)**1.5_r8+F21*(SM1**2/(ETA1+3.0_r8*PSI))**0.75_r8)
                               DUM2 = 1.0_r8/SM2**2*(F12*(PSI/ETA2)**1.5_r8+F22*(SM2**2/(ETA2+3.0_r8*PSI))**0.75_r8)

                               SMAX = 1.0_r8/(DUM1+DUM2)**0.5_r8

                               UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                               UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                               DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                               DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                               DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                               ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                               DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                               ! USE LOGNORMAL AEROSOL
                               SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                               AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)

                               SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                               SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                               ! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
                               SMAX = DUM3/QVS(K)

                               UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                               UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                               DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                               DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                               DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                               ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                               DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                               ! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
                               DUM2=MIN(DUM2,DUMACT)

                               DUM2 = (DUM2-NC3D(K))/DT
                               DUM2 = MAX(0.0_r8,DUM2)
                               NC3DTEN(K) = NC3DTEN(K)+DUM2

                            END IF ! IACT
                         END IF ! DUM3/QVS > 1.E-6
                      END IF  ! IDROP = 1

                      !.......................................................................
                   ELSE IF (IBASE.EQ.2) THEN

                      IF (IACT.EQ.1) THEN
                         ! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
                         ! BASED ON TWOMEY 1959

                         DUM=DUM*100.0_r8  ! CONVERT FROM M/S TO CM/S
                         DUM2 = 0.88_r8*C1**(2.0_r8/(K1+2.0_r8))*(7.E-2_r8*DUM**1.5_r8)**(K1/(K1+2.0_r8))
                         DUM2=DUM2*1.E6_r8 ! CONVERT FROM CM-3 TO M-3
                         DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                         DUM2 = (DUM2-NC3D(K))/DT
                         DUM2 = MAX(0.0_r8,DUM2)
                         NC3DTEN(K) = NC3DTEN(K)+DUM2

                      ELSE IF (IACT.EQ.2) THEN

                         SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                         AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)
                         ALPHA = SHR_CONST_G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-SHR_CONST_G*MA/(RR*T3D(K))
                         GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

                         GG = 1.0_r8/(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                              (T3D(K)*RR)-1.0_r8))

                         PSI = 2.0_r8/3.0_r8*(ALPHA*DUM/GG)**0.5_r8*AACT

                         ETA1 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW1)
                         ETA2 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW2)

                         SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                         SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                         DUM1 = 1.0_r8/SM1**2*(F11*(PSI/ETA1)**1.5_r8+F21*(SM1**2/(ETA1+3.0_r8*PSI))**0.75_r8)
                         DUM2 = 1.0_r8/SM2**2*(F12*(PSI/ETA2)**1.5_r8+F22*(SM2**2/(ETA2+3.0_r8*PSI))**0.75_r8)

                         SMAX = 1.0_r8/(DUM1+DUM2)**0.5_r8

                         UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                         UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                         DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                         DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                         DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                         ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                         DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                         DUM2 = (DUM2-NC3D(K))/DT
                         DUM2 = MAX(0.0_r8,DUM2)
                         NC3DTEN(K) = NC3DTEN(K)+DUM2
                      END IF  ! IACT
                   END IF  ! IBASE
                END IF  ! W > 0.001
             END IF  ! QC3D > QSMALL
          END IF  ! INUM = 0
          !!amy to here

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
          ! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
          ! LOSS OF NUMBER CONCENTRATION

          !     IF (PCC(K).LT.0.) THEN
          !        DUM = PCC(K)*DT/QC3D(K)
          !           DUM = MAX(-1.,DUM)
          !        NSUBC(K) = DUM*NC3D(K)/DT
          !     END IF

          ! UPDATE TENDENCIES

          !        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

          !.....................................................................
          !.....................................................................
       ELSE  ! TEMPERATURE < 273.15

          !......................................................................
          !HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
          ! INUM = 0, PREDICT DROPLET NUMBER
          ! INUM = 1, SET CONSTANT DROPLET NUMBER

          !!amy 
          IF (INUM.EQ.1) THEN
             ! CONVERT NDCNST FROM CM-3 TO KG-1
             NC3D(K)=NDCNST*1.E6_r8/RHO(K)
          END IF

          ! CALCULATE SIZE DISTRIBUTION PARAMETERS
          ! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

          NI3D(K) = MAX(0.0_r8,NI3D(K))
          NS3D(K) = MAX(0.0_r8,NS3D(K))
          NC3D(K) = MAX(0.0_r8,NC3D(K))
          NR3D(K) = MAX(0.0_r8,NR3D(K))
          NG3D(K) = MAX(0.0_r8,NG3D(K))

          !......................................................................
          ! CLOUD ICE

          IF (QI3D(K).GE.QSMALL) THEN
             LAMI(K) = (CONS12*                 &
                  NI3D(K)/QI3D(K))**(1.0_r8/DI)
             N0I(K) = NI3D(K)*LAMI(K)

             ! CHECK FOR SLOPE

             ! ADJUST VARS

             IF (LAMI(K).LT.LAMMINI) THEN

                LAMI(K) = LAMMINI

                N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

                NI3D(K) = N0I(K)/LAMI(K)
             ELSE IF (LAMI(K).GT.LAMMAXI) THEN
                LAMI(K) = LAMMAXI
                N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

                NI3D(K) = N0I(K)/LAMI(K)
             END IF
          END IF

          !......................................................................
          ! RAIN

          IF (QR3D(K).GE.QSMALL) THEN
             LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1.0_r8/3.0_r8)
             N0RR(K) = NR3D(K)*LAMR(K)

             ! CHECK FOR SLOPE

             ! ADJUST VARS

             IF (LAMR(K).LT.LAMMINR) THEN

                LAMR(K) = LAMMINR

                N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                NR3D(K) = N0RR(K)/LAMR(K)
             ELSE IF (LAMR(K).GT.LAMMAXR) THEN
                LAMR(K) = LAMMAXR
                N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                NR3D(K) = N0RR(K)/LAMR(K)
             END IF
          END IF

          !......................................................................
          ! CLOUD DROPLETS

          ! MARTIN ET AL. (1994) FORMULA FOR PGAM

          IF (QC3D(K).GE.QSMALL) THEN

             DUM = PRES(K)/(287.15_r8*T3D(K))
             PGAM(K)=0.0005714_r8*(NC3D(K)/1.E6_r8*DUM)+0.2714
             PGAM(K)=1.0_r8/(PGAM(K)**2)-1.0_r8
             PGAM(K)=MAX(PGAM(K),2.0_r8)
             PGAM(K)=MIN(PGAM(K),10.0_r8)

             ! CALCULATE LAMC

             LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.0_r8)/   &
                  (QC3D(K)*GAMMA(PGAM(K)+1.0_r8)))**(1.0_r8/3.0_r8)

             ! LAMMIN, 60 MICRON DIAMETER
             ! LAMMAX, 1 MICRON

             LAMMIN = (PGAM(K)+1.0_r8)/60.E-6_r8
             LAMMAX = (PGAM(K)+1.0_r8)/1.E-6_r8

             IF (LAMC(K).LT.LAMMIN) THEN
                LAMC(K) = LAMMIN

                NC3D(K) = EXP(3.0_r8*LOG(LAMC(K))+LOG(QC3D(K))+              &
                     LOG(GAMMA(PGAM(K)+1.0_r8))-LOG(GAMMA(PGAM(K)+4.0_r8)))/CONS26
             ELSE IF (LAMC(K).GT.LAMMAX) THEN
                LAMC(K) = LAMMAX
                NC3D(K) = EXP(3.0_r8*LOG(LAMC(K))+LOG(QC3D(K))+              &
                     LOG(GAMMA(PGAM(K)+1.0_r8))-LOG(GAMMA(PGAM(K)+4.0_r8)))/CONS26

             END IF

             ! TO CALCULATE DROPLET FREEZING

             CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.0_r8)

          END IF

          !......................................................................
          ! SNOW

          IF (QNI3D(K).GE.QSMALL) THEN
             LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1.0_r8/DS)
             N0S(K) = NS3D(K)*LAMS(K)

             ! CHECK FOR SLOPE

             ! ADJUST VARS

             IF (LAMS(K).LT.LAMMINS) THEN
                LAMS(K) = LAMMINS
                N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                NS3D(K) = N0S(K)/LAMS(K)

             ELSE IF (LAMS(K).GT.LAMMAXS) THEN

                LAMS(K) = LAMMAXS
                N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                NS3D(K) = N0S(K)/LAMS(K)
             END IF
          END IF

          !......................................................................
          ! GRAUPEL

          IF (QG3D(K).GE.QSMALL) THEN
             LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1.0_r8/DG)
             N0G(K) = NG3D(K)*LAMG(K)

             ! CHECK FOR SLOPE

             ! ADJUST VARS

             IF (LAMG(K).LT.LAMMING) THEN
                LAMG(K) = LAMMING
                N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                NG3D(K) = N0G(K)/LAMG(K)

             ELSE IF (LAMG(K).GT.LAMMAXG) THEN

                LAMG(K) = LAMMAXG
                N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                NG3D(K) = N0G(K)/LAMG(K)
             END IF
          END IF

          !.....................................................................
          ! ZERO OUT PROCESS RATES

          MNUCCC(K) = 0.0_r8
          NNUCCC(K) = 0.0_r8
          PRC(K) = 0.0_r8
          NPRC(K) = 0.0_r8
          NPRC1(K) = 0.0_r8
          NSAGG(K) = 0.0_r8
          PSACWS(K) = 0.0_r8
          NPSACWS(K) = 0.0_r8
          PSACWI(K) = 0.0_r8
          NPSACWI(K) = 0.0_r8
          PRACS(K) = 0.0_r8
          NPRACS(K) = 0.0_r8
          NMULTS(K) = 0.0_r8
          QMULTS(K) = 0.0_r8
          NMULTR(K) = 0.0_r8
          QMULTR(K) = 0.0_r8
          NMULTG(K) = 0.0_r8
          QMULTG(K) = 0.0_r8
          NMULTRG(K) = 0.0_r8
          QMULTRG(K) = 0.0_r8
          MNUCCR(K) = 0.0_r8
          NNUCCR(K) = 0.0_r8
          PRA(K) = 0.0_r8
          NPRA(K) = 0.0_r8
          NRAGG(K) = 0.0_r8
          PRCI(K) = 0.0_r8
          NPRCI(K) = 0.0_r8
          PRAI(K) = 0.0_r8
          NPRAI(K) = 0.0_r8
          NNUCCD(K) = 0.0_r8
          MNUCCD(K) = 0.0_r8
          PCC(K) = 0.0_r8
          PRE(K) = 0.0_r8
          PRD(K) = 0.0_r8
          PRDS(K) = 0.0_r8
          EPRD(K) = 0.0_r8
          EPRDS(K) = 0.0_r8
          NSUBC(K) = 0.0_r8
          NSUBI(K) = 0.0_r8
          NSUBS(K) = 0.0_r8
          NSUBR(K) = 0.0_r8
          PIACR(K) = 0.0_r8
          NIACR(K) = 0.0_r8
          PRACI(K) = 0.0_r8
          PIACRS(K) = 0.0_r8
          NIACRS(K) = 0.0_r8
          PRACIS(K) = 0.0_r8
          ! HM: ADD GRAUPEL PROCESSES
          PRACG(K) = 0.0_r8
          PSACR(K) = 0.0_r8
          PSACWG(K) = 0.0_r8
          PGSACW(K) = 0.0_r8
          PGRACS(K) = 0.0_r8
          PRDG(K) = 0.0_r8
          EPRDG(K) = 0.0_r8
          NPRACG(K) = 0.0_r8
          NPSACWG(K) = 0.0_r8
          NSCNG(K) = 0.0_r8
          NGRACS(K) = 0.0_r8
          NSUBG(K) = 0.0_r8

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! CALCULATION OF MICROPHYSICAL PROCESS RATES
          ! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
          !.......................................................................
          ! FREEZING OF CLOUD DROPLETS
          ! ONLY ALLOWED BELOW -4 C
          IF (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15_r8) THEN

             ! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
             ! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

             ! MEYERS CURVE
             NACNT = EXP(-2.80_r8+0.262_r8*(273.15_r8-T3D(K)))*1000.0_r8

             ! COOPER CURVE
             !        NACNT =  5.*EXP(0.304*(273.15-T3D(K)))

             ! FLECTHER
             !     NACNT = 0.01*EXP(0.6*(273.15-T3D(K)))

             ! CONTACT FREEZING

             ! MEAN FREE PATH

             DUM = 7.37_r8*T3D(K)/(288.0_r8*10.0_r8*PRES(K))/100.0_r8

             ! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
             ! BASED ON BROWNIAN DIFFUSION

             DAP(K) = CONS37*T3D(K)*(1.0_r8+DUM/RIN)/MU(K)

             MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                  LOG(GAMMA(PGAM(K)+5.0_r8))-4.0_r8*LOG(LAMC(K)))
             NNUCCC(K) = 2.0_r8*PI*DAP(K)*NACNT*CDIST1(K)*           &
                  GAMMA(PGAM(K)+2.0_r8)/                         &
                  LAMC(K)

             ! IMMERSION FREEZING (BIGG 1953)

             !           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
             !                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
             !                   EXP(AIMM*(273.15-T3D(K)))

             !           NNUCCC(K) = NNUCCC(K)+                                  &
             !            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
             !                *EXP(AIMM*(273.15-T3D(K)))

             ! hm 7/15/13 fix for consistency w/ original formula
             MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.0_r8+PGAM(K)))-6.0_r8*LOG(LAMC(K)))*             &
                  (EXP(AIMM*(273.15_r8-T3D(K)))-1.0_r8)

             NNUCCC(K) = NNUCCC(K)+                                  &
                  CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.0_r8))-3.0_r8*LOG(LAMC(K)))              &
                  *(EXP(AIMM*(273.15_r8-T3D(K)))-1.0_r8)

             ! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
             ! MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC

             NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

          END IF

          !.................................................................
          !.......................................................................
          ! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
          ! FORMULA FROM BEHENG (1994)
          ! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
          ! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
          ! AS A GAMMA DISTRIBUTION

          ! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

          IF (QC3D(K).GE.1.E-6_r8) THEN

             ! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
             ! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

             PRC(K)=1350.0_r8*QC3D(K)**2.47_r8*  &
                  (NC3D(K)/1.e6_r8*RHO(K))**(-1.79_r8)

             ! note: nprc1 is change in Nr,
             ! nprc is change in Nc

             NPRC1(K) = PRC(K)/CONS29
             NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

             ! hm bug fix 3/4/13
             NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
             NPRC1(K) = MIN(NPRC1(K),NPRC(K))

          END IF

          !.......................................................................
          ! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

          ! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
          ! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

          IF (QNI3D(K).GE.1.E-8_r8) THEN
             NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
                  ((2.0_r8+BS)/3.0_r8)*QNI3D(K)**((2.0_r8+BS)/3.0_r8)*                  &
                  (NS3D(K)*RHO(K))**((4.0_r8-BS)/3.0_r8)/                       &
                  (RHO(K))
          END IF

          !.......................................................................
          ! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
          ! HERE USE CONTINUOUS COLLECTION EQUATION WITH
          ! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

          ! SNOW

          IF (QNI3D(K).GE.1.E-8_r8 .AND. QC3D(K).GE.QSMALL) THEN

             PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.0_r8)
             NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.0_r8)

          END IF

          !............................................................................
          ! COLLECTION OF CLOUD WATER BY GRAUPEL

          IF (QG3D(K).GE.1.E-8_r8 .AND. QC3D(K).GE.QSMALL) THEN

             PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.0_r8)
             NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.0_r8)
          END IF

          !.......................................................................
          ! HM, ADD 12/13/06
          ! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
          ! BEFORE RIMING CAN OCCUR
          ! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
          ! TO HALLET-MOSSOP SPLINTERING

          IF (QI3D(K).GE.1.E-8_r8 .AND. QC3D(K).GE.QSMALL) THEN

             ! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
             ! FROM THOMPSON ET AL. 2004, MWR

             IF (1.0_r8/LAMI(K).GE.100.E-6_r8) THEN

                PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                     N0I(K)/                        &
                     LAMI(K)**(BI+3.0_r8)
                NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                     N0I(K)/                        &
                     LAMI(K)**(BI+3.0_r8)
             END IF
          END IF

          !.......................................................................
          ! ACCRETION OF RAIN WATER BY SNOW
          ! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

          IF (QR3D(K).GE.1.E-8_r8.AND.QNI3D(K).GE.1.E-8_r8) THEN

             UMS = ASN(K)*CONS3/(LAMS(K)**BS)
             UMR = ARN(K)*CONS4/(LAMR(K)**BR)
             UNS = ASN(K)*CONS5/LAMS(K)**BS
             UNR = ARN(K)*CONS6/LAMR(K)**BR

             ! SET REASLISTIC LIMITS ON FALLSPEEDS
             ! bug fix, 10/08/09
             dum=(rhosu/rho(k))**0.54_r8
             UMS=MIN(UMS,1.2_r8*dum)
             UNS=MIN(UNS,1.2_r8*dum)
             UMR=MIN(UMR,9.1_r8*dum)
             UNR=MIN(UNR,9.1_r8*dum)

             PRACS(K) = CONS41*(((1.2_r8*UMR-0.95_r8*UMS)**2+                   &
                  0.08_r8*UMS*UMR)**0.5_r8*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5.0_r8/(LAMR(K)**3*LAMS(K))+                    &
                  2.0_r8/(LAMR(K)**2*LAMS(K)**2)+                  &                                 
                  0.5_r8/(LAMR(k)*LAMS(k)**3)))

             NPRACS(K) = CONS32*RHO(K)*(1.7_r8*(UNR-UNS)**2+            &
                  0.3_r8*UNR*UNS)**0.5_r8*N0RR(K)*N0S(K)*              &
                  (1.0_r8/(LAMR(K)**3*LAMS(K))+                      &
                  1.0_r8/(LAMR(K)**2*LAMS(K)**2)+                   &
                  1.0_r8/(LAMR(K)*LAMS(K)**3))

             ! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
             ! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
             ! RIME-SPLINTERING

             PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

             ! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
             ! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

             ! HM MODIFY FOR WRFV3.1
             !            IF (IHAIL.EQ.0) THEN
             IF (QNI3D(K).GE.0.1E-3_r8.AND.QR3D(K).GE.0.1E-3_r8) THEN
                PSACR(K) = CONS31*(((1.2_r8*UMR-0.95_r8*UMS)**2+              &
                     0.08_r8*UMS*UMR)**0.5_r8*RHO(K)*                     &
                     N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                     (5.0_r8/(LAMS(K)**3*LAMR(K))+                    &
                     2.0_r8/(LAMS(K)**2*LAMR(K)**2)+                  &
                     0.5_r8/(LAMS(K)*LAMR(K)**3)))            
             END IF
             !            END IF

          END IF

          !.......................................................................

          ! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990, 
          ! USED BY REISNER ET AL 1998
          IF (QR3D(K).GE.1.E-8_r8.AND.QG3D(K).GE.1.E-8_r8) THEN

             UMG = AGN(K)*CONS7/(LAMG(K)**BG)
             UMR = ARN(K)*CONS4/(LAMR(K)**BR)
             UNG = AGN(K)*CONS8/LAMG(K)**BG
             UNR = ARN(K)*CONS6/LAMR(K)**BR

             ! SET REASLISTIC LIMITS ON FALLSPEEDS
             ! bug fix, 10/08/09
             dum=(rhosu/rho(k))**0.54_r8
             UMG=MIN(UMG,20.0_r8*dum)
             UNG=MIN(UNG,20.0_r8*dum)
             UMR=MIN(UMR,9.1_r8*dum)
             UNR=MIN(UNR,9.1_r8*dum)

             PRACG(K) = CONS41*(((1.2_r8*UMR-0.95_r8*UMG)**2+                   &
                  0.08_r8*UMG*UMR)**0.5_r8*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5.0_r8/(LAMR(K)**3*LAMG(K))+                    &
                  2.0_r8/(LAMR(K)**2*LAMG(K)**2)+                                   &
                  0.5_r8/(LAMR(k)*LAMG(k)**3)))

             NPRACG(K) = CONS32*RHO(K)*(1.7_r8*(UNR-UNG)**2+            &
                  0.3_r8*UNR*UNG)**0.5_r8*N0RR(K)*N0G(K)*              &
                  (1.0_r8/(LAMR(K)**3*LAMG(K))+                      &
                  1.0_r8/(LAMR(K)**2*LAMG(K)**2)+                   &
                  1.0_r8/(LAMR(K)*LAMG(K)**3))

             ! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
             ! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
             ! RIME-SPLINTERING

             PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

          END IF

          !.......................................................................
          ! RIME-SPLINTERING - SNOW
          ! HALLET-MOSSOP (1974)
          ! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

          ! DUM1 = MASS OF INDIVIDUAL SPLINTERS

          ! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
          ! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
          ! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

          !v1.4
          IF (QNI3D(K).GE.0.1E-3_r8) THEN
             IF (QC3D(K).GE.0.5E-3_r8.OR.QR3D(K).GE.0.1E-3_r8) THEN
                IF (PSACWS(K).GT.0.0_r8.OR.PRACS(K).GT.0.0_r8) THEN
                   IF (T3D(K).LT.270.16_r8 .AND. T3D(K).GT.265.16_r8) THEN

                      IF (T3D(K).GT.270.16_r8) THEN
                         FMULT = 0.0_r8
                      ELSE IF (T3D(K).LE.270.16_r8.AND.T3D(K).GT.268.16_r8)  THEN
                         FMULT = (270.16_r8-T3D(K))/2.0_r8
                      ELSE IF (T3D(K).GE.265.16_r8.AND.T3D(K).LE.268.16_r8)   THEN
                         FMULT = (T3D(K)-265.16_r8)/3.0_r8
                      ELSE IF (T3D(K).LT.265.16_r8) THEN
                         FMULT = 0.0_r8
                      END IF

                      ! 1000 IS TO CONVERT FROM KG TO G

                      ! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

                      IF (PSACWS(K).GT.0.0_r8) THEN
                         NMULTS(K) = 35.E4_r8*PSACWS(K)*FMULT*1000.0_r8
                         QMULTS(K) = NMULTS(K)*MMULT

                         ! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
                         ! THAN WAS RIMED ONTO SNOW

                         QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                         PSACWS(K) = PSACWS(K)-QMULTS(K)

                      END IF

                      ! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

                      IF (PRACS(K).GT.0.0_r8) THEN
                         NMULTR(K) = 35.E4_r8*PRACS(K)*FMULT*1000.0_r8
                         QMULTR(K) = NMULTR(K)*MMULT

                         ! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
                         ! THAN WAS RIMED ONTO SNOW

                         QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                         PRACS(K) = PRACS(K)-QMULTR(K)

                      END IF

                   END IF
                END IF
             END IF
          END IF

          !.......................................................................
          ! RIME-SPLINTERING - GRAUPEL 
          ! HALLET-MOSSOP (1974)
          ! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

          ! DUM1 = MASS OF INDIVIDUAL SPLINTERS

          ! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
          ! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

          !         IF (IHAIL.EQ.0) THEN
          ! v1.4
          IF (QG3D(K).GE.0.1E-3_r8) THEN
             IF (QC3D(K).GE.0.5E-3_r8 .OR.QR3D(K).GE.0.1E-3_r8) THEN
                IF (PSACWG(K).GT.0.0_r8.OR.PRACG(K).GT.0.0_r8) THEN
                   IF (T3D(K).LT.270.16_r8 .AND. T3D(K).GT.265.16_r8) THEN

                      IF (T3D(K).GT.270.16_r8) THEN
                         FMULT = 0.0_r8
                      ELSE IF (T3D(K).LE.270.16_r8.AND.T3D(K).GT.268.16_r8)  THEN
                         FMULT = (270.16_r8-T3D(K))/2.0_r8
                      ELSE IF (T3D(K).GE.265.16_r8.AND.T3D(K).LE.268.16_r8)   THEN
                         FMULT = (T3D(K)-265.16_r8)/3.0_r8
                      ELSE IF (T3D(K).LT.265.16_r8) THEN
                         FMULT = 0.0_r8
                      END IF

                      ! 1000 IS TO CONVERT FROM KG TO G

                      ! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

                      IF (PSACWG(K).GT.0.0_r8) THEN
                         NMULTG(K) = 35.E4_r8*PSACWG(K)*FMULT*1000.0_r8
                         QMULTG(K) = NMULTG(K)*MMULT

                         ! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
                         ! THAN WAS RIMED ONTO GRAUPEL

                         QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                         PSACWG(K) = PSACWG(K)-QMULTG(K)

                      END IF

                      ! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

                      IF (PRACG(K).GT.0.0_r8) THEN
                         NMULTRG(K) = 35.E4_r8*PRACG(K)*FMULT*1000.0_r8
                         QMULTRG(K) = NMULTRG(K)*MMULT

                         ! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
                         ! THAN WAS RIMED ONTO GRAUPEL

                         QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                         PRACG(K) = PRACG(K)-QMULTRG(K)

                      END IF
                   END IF
                END IF
             END IF
          END IF
          !         END IF

          !........................................................................
          ! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL

          !           IF (IHAIL.EQ.0) THEN
          IF (PSACWS(K).GT.0.0_r8) THEN
             ! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
             IF (QNI3D(K).GE.0.1E-3_r8.AND.QC3D(K).GE.0.5E-3_r8) THEN

                ! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
                PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
                     ASN(K)*ASN(K)/ &
                     (RHO(K)*LAMS(K)**(2.0_r8*BS+2.0_r8))) 

                ! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
                DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.0_r8) 

                ! NUMBER CONCENTRATION OF EMBRYO GRAUPEL FROM RIMING OF SNOW
                NSCNG(K) = DUM/MG0*RHO(K)
                ! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
                NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

                ! PORTION OF RIMING LEFT FOR SNOW
                PSACWS(K) = PSACWS(K) - PGSACW(K)
             END IF
          END IF

          ! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

          IF (PRACS(K).GT.0.0_r8) THEN
             ! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
             IF (QNI3D(K).GE.0.1E-3_r8.AND.QR3D(K).GE.0.1E-3_r8) THEN
                ! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
                DUM = CONS18*(4.0_r8/LAMS(K))**3*(4.0_r8/LAMS(K))**3 &    
                     /(CONS18*(4.0_r8/LAMS(K))**3*(4.0_r8/LAMS(K))**3+ &  
                     CONS19*(4.0_r8/LAMR(K))**3*(4.0_r8/LAMR(K))**3)
                DUM=MIN(DUM,1.0_r8)
                DUM=MAX(DUM,0.0_r8)
                PGRACS(K) = (1.0_r8-DUM)*PRACS(K)
                NGRACS(K) = (1.0_r8-DUM)*NPRACS(K)
                ! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
                NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
                NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

                ! AMOUNT LEFT FOR SNOW PRODUCTION
                PRACS(K) = PRACS(K) - PGRACS(K)
                NPRACS(K) = NPRACS(K) - NGRACS(K)
                ! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
                PSACR(K)=PSACR(K)*(1.0_r8-DUM)
             END IF
          END IF
          !           END IF

          !.......................................................................
          ! FREEZING OF RAIN DROPS
          ! FREEZING ALLOWED BELOW -4 C

          IF (T3D(K).LT.269.15_r8.AND.QR3D(K).GE.QSMALL) THEN

             ! IMMERSION FREEZING (BIGG 1953)
             !            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3 &
             !                 /LAMR(K)**3

             !            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3

             ! hm fix 7/15/13 for consistency w/ original formula
             MNUCCR(K) = CONS20*NR3D(K)*(EXP(AIMM*(273.15_r8-T3D(K)))-1.0_r8)/LAMR(K)**3 &
                  /LAMR(K)**3

             NNUCCR(K) = PI*NR3D(K)*BIMM*(EXP(AIMM*(273.15_r8-T3D(K)))-1.0_r8)/LAMR(K)**3

             ! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
             NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

          END IF

          !.......................................................................
          ! ACCRETION OF CLOUD LIQUID WATER BY RAIN
          ! CONTINUOUS COLLECTION EQUATION WITH
          ! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

          IF (QR3D(K).GE.1.E-8_r8 .AND. QC3D(K).GE.1.E-8_r8) THEN

             ! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
             ! KHAIROUTDINOV AND KOGAN 2000, MWR

             DUM=(QC3D(K)*QR3D(K))
             PRA(K) = 67.0_r8*(DUM)**1.15_r8
             NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

          END IF
          !.......................................................................
          ! SELF-COLLECTION OF RAIN DROPS
          ! FROM BEHENG(1994)
          ! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
          ! AS DESCRINED ABOVE FOR AUTOCONVERSION

          IF (QR3D(K).GE.1.E-8_r8) THEN
             ! include breakup add 10/09/09
             dum1=300.e-6_r8
             IF (1.0_r8/lamr(k).LT.dum1) THEN
                dum=1.0_r8
             ELSE IF (1.0_r8/lamr(k).GE.dum1) THEN
                dum=2.0_r8-EXP(2300.0_r8*(1.0_r8/lamr(k)-dum1))
             END IF
             !            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
             NRAGG(K) = -5.78_r8*dum*NR3D(K)*QR3D(K)*RHO(K)
          END IF

          !.......................................................................
          ! AUTOCONVERSION OF CLOUD ICE TO SNOW
          ! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
          ! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
          ! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

          IF (QI3D(K).GE.1.E-8_r8 .AND.QVQVSI(K).GE.1.0_r8) THEN

             !           COFFI = 2.0_r8/LAMI(K)
             !           IF (COFFI.GE.DCS) THEN
             NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                  *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
             PRCI(K) = CONS22*NPRCI(K)
             NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

             !           END IF
          END IF

          !.......................................................................
          ! ACCRETION OF CLOUD ICE BY SNOW
          ! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
          ! AND DS >> DI FOR CONTINUOUS COLLECTION

          IF (QNI3D(K).GE.1.E-8_r8.AND. QI3D(K).GE.QSMALL) THEN
             PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                  LAMS(K)**(BS+3.0_r8)
             NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                  RHO(K)*N0S(K)/                                 &
                  LAMS(K)**(BS+3.0_r8)
             NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
          END IF

          !.......................................................................
          ! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
          ! FOLLOWS REISNER ET AL. 1998
          ! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

          IF (QR3D(K).GE.1.E-8_r8.AND.QI3D(K).GE.1.E-8_r8.AND.T3D(K).LE.273.15_r8) THEN

             ! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
             ! OTHERWISE ADD TO SNOW

             IF (QR3D(K).GE.0.1E-3_r8) THEN
                NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                     /LAMR(K)**(BR+3.0_r8)*RHO(K)
                PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                     /LAMR(K)**(BR+3.0_r8)/LAMR(K)**3*RHO(K)
                PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                     LAMR(K)**(BR+3.0_r8)*RHO(K)
                NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
                NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
             ELSE 
                NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                     /LAMR(K)**(BR+3.0_r8)*RHO(K)
                PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                     /LAMR(K)**(BR+3.0_r8)/LAMR(K)**3*RHO(K)
                PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                     LAMR(K)**(BR+3.0_r8)*RHO(K)
                NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
                NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
             END IF
          END IF

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

          IF (INUC.EQ.0) THEN

             ! FREEZING OF AEROSOL ONLY ALLOWED BELOW -5 C
             ! AND ABOVE DELIQUESCENCE THRESHOLD OF 80%
             ! AND ABOVE ICE SATURATION

             ! add threshold according to Greg Thomspon

             IF ((QVQVS(K).GE.0.999_r8.AND.T3D(K).LE.265.15_r8).OR. &
                  QVQVSI(K).GE.1.08_r8) THEN

                ! hm, modify dec. 5, 2006, replace with cooper curve
                kc2 = 0.005_r8*EXP(0.304_r8*(273.15_r8-T3D(K)))*1000.0_r8 ! convert from L-1 to m-3
                ! limit to 500 L-1
                kc2 = MIN(kc2,500.e3_r8)
                kc2=MAX(kc2/rho(k),0.0_r8)  ! convert to kg-1

                IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
                   NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
                   MNUCCD(K) = NNUCCD(K)*MI0
                END IF

             END IF

          ELSE IF (INUC.EQ.1) THEN
             IF (T3D(K).LT.273.15_r8.AND.QVQVSI(K).GT.1.0_r8) THEN
                KC2 = 0.16_r8*1000.0_r8/RHO(K)  ! CONVERT FROM L-1 TO KG-1
                IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
                   NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
                   MNUCCD(K) = NNUCCD(K)*MI0
                END IF
             END IF

          END IF

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!101       CONTINUE

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

          ! NO VENTILATION FOR CLOUD ICE

          IF (QI3D(K).GE.QSMALL) THEN

             EPSI = 2.0_r8*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

          ELSE
             EPSI = 0.0_r8
          END IF

          IF (QNI3D(K).GE.QSMALL) THEN
             EPSS = 2.0_r8*PI*N0S(K)*RHO(K)*DV(K)*                        &
                  (F1S/(LAMS(K)*LAMS(K))+                       &
                  F2S*(ASN(K)*RHO(K)/MU(K))**0.5*              &
                  SC(K)**(1.0_r8/3.0_r8)*CONS10/                       &
                  (LAMS(K)**CONS35))
          ELSE
             EPSS = 0.0_r8
          END IF

          IF (QG3D(K).GE.QSMALL) THEN
             EPSG = 2.0_r8*PI*N0G(K)*RHO(K)*DV(K)*                        &
                  (F1S/(LAMG(K)*LAMG(K))+                       &
                  F2S*(AGN(K)*RHO(K)/MU(K))**0.5_r8*              &
                  SC(K)**(1.0_r8/3.0_r8)*CONS11/                       &
                  (LAMG(K)**CONS36))

          ELSE
             EPSG = 0.0_r8
          END IF

          IF (QR3D(K).GE.QSMALL) THEN
             EPSR = 2.0_r8*PI*N0RR(K)*RHO(K)*DV(K)*                       &
                  (F1R/(LAMR(K)*LAMR(K))+                       &
                  F2R*(ARN(K)*RHO(K)/MU(K))**0.5_r8*              &
                  SC(K)**(1.0_r8/3.0_r8)*CONS9/                        &
                  (LAMR(K)**CONS34))
          ELSE
             EPSR = 0.0_r8
          END IF

          ! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
          ! DUM IS FRACTION OF D*N(D) < DCS

          ! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
          IF (QI3D(K).GE.QSMALL) THEN              
             DUM=(1.0_r8-EXP(-LAMI(K)*DCS)*(1.0_r8+LAMI(K)*DCS))
             PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
          ELSE
             DUM=0.0_r8
          END IF
          ! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
          IF (QNI3D(K).GE.QSMALL) THEN
             PRDS(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
                  EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.0_r8-DUM)
             ! OTHERWISE ADD TO CLOUD ICE
          ELSE
             PRD(K) = PRD(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.0_r8-DUM)
          END IF
          ! VAPOR DEPOSITION ON GRAUPEL
          PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)

          ! NO CONDENSATION ONTO RAIN, ONLY EVAP

          IF (QV3D(K).LT.QVS(K)) THEN
             PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
             PRE(K) = MIN(PRE(K),0.0_r8)
          ELSE
             PRE(K) = 0.0_r8
          END IF

          ! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
          ! FORMULA FROM REISNER 2 SCHEME

          DUM = (QV3D(K)-QVI(K))/DT

          FUDGEF = 0.9999_r8
          SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)

          IF( (DUM.GT.0.0_r8 .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
               (DUM.LT.0.0_r8 .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
             MNUCCD(K) = FUDGEF*MNUCCD(K)*DUM/SUM_DEP
             PRD(K) = FUDGEF*PRD(K)*DUM/SUM_DEP
             PRDS(K) = FUDGEF*PRDS(K)*DUM/SUM_DEP
             PRDG(K) = FUDGEF*PRDG(K)*DUM/SUM_DEP
          ENDIF

          ! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

          IF (PRD(K).LT.0.0_r8) THEN
             EPRD(K)=PRD(K)
             PRD(K)=0.0_r8
          END IF
          IF (PRDS(K).LT.0.0_r8) THEN
             EPRDS(K)=PRDS(K)
             PRDS(K)=0.0_r8
          END IF
          IF (PRDG(K).LT.0.0_r8) THEN
             EPRDG(K)=PRDG(K)
             PRDG(K)=0.0_r8
          END IF

          !.......................................................................
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ! CONSERVATION OF WATER
          ! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
          ! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.

          ! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
          ! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

          ! NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
          ! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
          ! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
          ! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
          ! STEP

          !****SENSITIVITY - NO ICE

          IF (ILIQ.EQ.1) THEN
             MNUCCC(K)=0.0_r8
             NNUCCC(K)=0.0_r8
             MNUCCR(K)=0.0_r8
             NNUCCR(K)=0.0_r8
             MNUCCD(K)=0.0_r8
             NNUCCD(K)=0.0_r8
          END IF

          ! ****SENSITIVITY - NO GRAUPEL
          IF (IGRAUP.EQ.1) THEN
             PRACG(K) = 0.0_r8
             PSACR(K) = 0.0_r8
             PSACWG(K) = 0.0_r8
             PGSACW(K) = 0.0_r8
             PGRACS(K) = 0.0_r8
             PRDG(K) = 0.0_r8
             EPRDG(K) = 0.0_r8
             EVPMG(K) = 0.0_r8
             PGMLT(K) = 0.0_r8
             NPRACG(K) = 0.0_r8
             NPSACWG(K) = 0.0_r8
             NSCNG(K) = 0.0_r8
             NGRACS(K) = 0.0_r8
             NSUBG(K) = 0.0_r8
             NGMLTG(K) = 0.0_r8
             NGMLTR(K) = 0.0_r8
             ! fix 5/27/11
             PIACRS(K)=PIACRS(K)+PIACR(K)
             PIACR(K) = 0.0_r8
             ! fix 070713
             PRACIS(K)=PRACIS(K)+PRACI(K)
             PRACI(K) = 0.0_r8
             PSACWS(K)=PSACWS(K)+PGSACW(K)
             PGSACW(K) = 0.0_r8
             PRACS(K)=PRACS(K)+PGRACS(K)
             PGRACS(K) = 0.0_r8
          END IF

          ! CONSERVATION OF QC

          DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT

          IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN
             RATIO = (QC3D(K)/DUM)

             PRC(K) = PRC(K)*RATIO
             PRA(K) = PRA(K)*RATIO
             MNUCCC(K) = MNUCCC(K)*RATIO
             PSACWS(K) = PSACWS(K)*RATIO
             PSACWI(K) = PSACWI(K)*RATIO
             QMULTS(K) = QMULTS(K)*RATIO
             QMULTG(K) = QMULTG(K)*RATIO
             PSACWG(K) = PSACWG(K)*RATIO
             PGSACW(K) = PGSACW(K)*RATIO
          END IF

          ! CONSERVATION OF QI

          DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
               -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT

          IF (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) THEN

             RATIO = ((QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
                  MNUCCD(K)+PSACWI(K))/ &
                  (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K)))

             PRCI(K) = PRCI(K)*RATIO
             PRAI(K) = PRAI(K)*RATIO
             PRACI(K) = PRACI(K)*RATIO
             PRACIS(K) = PRACIS(K)*RATIO
             EPRD(K) = EPRD(K)*RATIO

          END IF

          ! CONSERVATION OF QR

          DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
               PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT

          IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

             RATIO = ((QR3D(K)/DT+PRC(K)+PRA(K))/ &
                  (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K)))

             PRE(K) = PRE(K)*RATIO
             PRACS(K) = PRACS(K)*RATIO
             QMULTR(K) = QMULTR(K)*RATIO
             QMULTRG(K) = QMULTRG(K)*RATIO
             MNUCCR(K) = MNUCCR(K)*RATIO
             PIACR(K) = PIACR(K)*RATIO
             PIACRS(K) = PIACRS(K)*RATIO
             PGRACS(K) = PGRACS(K)*RATIO
             PRACG(K) = PRACG(K)*RATIO

          END IF

          ! CONSERVATION OF QNI
          ! CONSERVATION FOR GRAUPEL SCHEME

          IF (IGRAUP.EQ.0) THEN

             DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT

             IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

                RATIO = ((QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K)))

                EPRDS(K) = EPRDS(K)*RATIO
                PSACR(K) = PSACR(K)*RATIO

             END IF

             ! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
          ELSE IF (IGRAUP.EQ.1) THEN

             DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT

             IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

                RATIO = ((QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K) &
                     + PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K)))

                EPRDS(K) = EPRDS(K)*RATIO
                PSACR(K) = PSACR(K)*RATIO

             END IF

          END IF

          ! CONSERVATION OF QG

          DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT

          IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

             RATIO = ((QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
                  PIACR(K)+PRACI(K))/(-EPRDG(K)))

             EPRDG(K) = EPRDG(K)*RATIO

          END IF

          ! TENDENCIES

          QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))

          ! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
          T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
               *XXLV(K)+(PRD(K)+PRDS(K)+                            &
               MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
               (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
               QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
               +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*XLF(K))/CPM(K)

          QC3DTEN(K) = QC3DTEN(K)+                                      &
               (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
               PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
          QI3DTEN(K) = QI3DTEN(K)+                                      &
               (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
               PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
          QR3DTEN(K) = QR3DTEN(K)+                                      &
               (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
               -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))

          IF (IGRAUP.EQ.0) THEN

             QNI3DTEN(K) = QNI3DTEN(K)+                                    &
                  (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
             NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
             QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
                  PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
             NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

             ! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
          ELSE IF (IGRAUP.EQ.1) THEN

             QNI3DTEN(K) = QNI3DTEN(K)+                                    &
                  (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
             NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))

          END IF

          NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
               -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))

          NI3DTEN(K) = NI3DTEN(K)+                                      &
               (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
               NNUCCD(K)-NIACR(K)-NIACRS(K))

          NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
               +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
          ! WATER SATURATION

          DUMT = T3D(K)+DT*T3DTEN(K)
          DUMQV = QV3D(K)+DT*QV3DTEN(K)
          DUMQSS = 0.622_r8*POLYSVP(DUMT,0)/ (PRES(K)-POLYSVP(DUMT,0))
          DUMQC = QC3D(K)+DT*QC3DTEN(K)
          DUMQC = MAX(DUMQC,0.0_r8)

          ! SATURATION ADJUSTMENT FOR LIQUID

          DUMS = DUMQV-DUMQSS
          PCC(K) = DUMS/(1.0_r8+XXLV(K)**2*DUMQSS/(CPM(K)*r_v*DUMT**2))/DT
          IF (PCC(K)*DT+DUMQC.LT.0.0_r8) THEN
             PCC(K) = -DUMQC/DT
          END IF

          QV3DTEN(K) = QV3DTEN(K)-PCC(K)
          T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
          QC3DTEN(K) = QC3DTEN(K)+PCC(K)

          !.......................................................................
          ! ACTIVATION OF CLOUD DROPLETS
          ! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
          ! DROPLET CONCENTRATION IS SPECIFIED !!!!!

          !!amy added code to predict droplet concentration back in 
          IF (INUM.EQ.0) THEN     

             IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) THEN

                ! EFFECTIVE VERTICAL VELOCITY (M/S)

                IF (ISUB.EQ.0) THEN
                   ! ADD SUB-GRID VERTICAL VELOCITY
                   DUM = W3D(K)+WVAR(K)

                   ! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
                   DUM = MAX(DUM,0.10_r8)

                ELSE IF (ISUB.EQ.1) THEN
                   DUM=W3D(K)
                END IF

                ! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
                IF (DUM.GE.0.001_r8) THEN

                   IF (IBASE.EQ.1) THEN

                      ! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
                      ! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

                      IDROP=0

                      IF (QC3D(K)+QC3DTEN(K)*DT.LE.0.05E-3_r8/RHO(K)) THEN
                         IDROP=1
                      END IF
                      IF (K.EQ.1) THEN
                         IDROP=1
                      ELSE IF (K.GE.2) THEN
                         IF (QC3D(K)+QC3DTEN(K)*DT.GT.0.05E-3_r8/RHO(K).AND. &
                              QC3D(K-1)+QC3DTEN(K-1)*DT.LE.0.05E-3_r8/RHO(K-1)) THEN
                            IDROP=1
                         END IF
                      END IF

                      IF (IDROP.EQ.1) THEN
                         ! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

                         IF (IACT.EQ.1) THEN
                            ! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
                            ! BASED ON TWOMEY 1959

                            DUM=DUM*100.0_r8  ! CONVERT FROM M/S TO CM/S
                            DUM2 = 0.88_r8*C1**(2.0_r8/(K1+2.0_r8))*(7.E-2_r8*DUM**1.5_r8)**(K1/(K1+2.0_r8))
                            DUM2=DUM2*1.E6_r8 ! CONVERT FROM CM-3 TO M-3
                            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                            DUM2 = (DUM2-NC3D(K))/DT
                            DUM2 = MAX(0.0_r8,DUM2)
                            NC3DTEN(K) = NC3DTEN(K)+DUM2

                         ELSE IF (IACT.EQ.2) THEN
                            ! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

                            SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                            AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)
                            ALPHA = SHR_CONST_G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-SHR_CONST_G*MA/(RR*T3D(K))
                            GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

                            GG = 1.0_r8/(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                 (T3D(K)*RR)-1.0_r8))

                            PSI = 2.0_r8/3.0_r8*(ALPHA*DUM/GG)**0.5_r8*AACT

                            ETA1 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW1)
                            ETA2 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW2)

                            SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                            SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                            DUM1 = 1.0_r8/SM1**2*(F11*(PSI/ETA1)**1.5_r8+F21*(SM1**2/(ETA1+3.0_r8*PSI))**0.75_r8)
                            DUM2 = 1.0_r8/SM2**2*(F12*(PSI/ETA2)**1.5_r8+F22*(SM2**2/(ETA2+3.0_r8*PSI))**0.75_r8)

                            SMAX = 1.0_r8/(DUM1+DUM2)**0.5_r8

                            UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                            UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                            DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                            DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                            DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                            ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                            DUM2 = (DUM2-NC3D(K))/DT
                            DUM2 = MAX(0.0_r8,DUM2)
                            NC3DTEN(K) = NC3DTEN(K)+DUM2
                         END IF  ! IACT

                         !.............................................................................
                      ELSE IF (IDROP.EQ.0) THEN
                         ! ACTIVATE IN CLOUD INTERIOR
                         ! FIND EQUILIBRIUM SUPERSATURATION

                         TAUC=1.0_r8/(2.0_r8*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.0_r8)/LAMC(K))
                         IF (EPSR.GT.1.E-8_r8) THEN
                            TAUR=1.0_r8/EPSR
                         ELSE
                            TAUR=1.E8_r8
                         END IF
                         !!amy taui,taus,taug lines added in v3
                         IF (EPSI.GT.1.E-8_r8) THEN
                            TAUI=1.0_r8/EPSI
                         ELSE
                            TAUI=1.E8_r8
                         END IF
                         IF (EPSS.GT.1.E-8_r8) THEN
                            TAUS=1.0_r8/EPSS
                         ELSE
                            TAUS=1.E8_r8
                         END IF
                         IF (EPSG.GT.1.E-8_r8) THEN
                            TAUG=1.0_r8/EPSG
                         ELSE
                            TAUG=1.E8_r8
                         END IF

                         ! EQUILIBRIUM SS INCLUDING BERGERON EFFECT
                         !!amy added taui,taus,taug to these lines in v3
                         ! hm fix 1/20/15
                         !           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
                         DUM3=(-QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*SHR_CONST_G*DUM

                         DUM3=(DUM3*TAUC*TAUR*TAUI*TAUS*TAUG- &
                              (QVS(K)-QVI(K))*(TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS))/ &
                              (TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS+ &
                              TAUR*TAUI*TAUS*TAUG+TAUC*TAUI*TAUS*TAUG)

                         IF (DUM3/QVS(K).GE.1.E-6_r8) THEN
                            IF (IACT.EQ.1) THEN

                               ! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQUILIBRIUM SS

                               DUM=DUM*100.0_r8  ! CONVERT FROM M/S TO CM/S
                               DUMACT = 0.88_r8*C1**(2.0_r8/(K1+2.0_r8))*(7.E-2_r8*DUM**1.5_r8)**(K1/(K1+2.0_r8))

                               ! USE POWER LAW CCN SPECTRA

                               ! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
                               DUM3=DUM3/QVS(K)*100.0_r8

                               DUM2=C1*DUM3**K1
                               ! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
                               DUM2=MIN(DUM2,DUMACT)
                               DUM2=DUM2*1.E6_r8 ! CONVERT FROM CM-3 TO M-3
                               DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                               DUM2 = (DUM2-NC3D(K))/DT
                               DUM2 = MAX(0.0_r8,DUM2)
                               NC3DTEN(K) = NC3DTEN(K)+DUM2

                            ELSE IF (IACT.EQ.2) THEN

                               ! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQUILIBRIUM SS

                               SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                               AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)
                               ALPHA = SHR_CONST_G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-SHR_CONST_G*MA/(RR*T3D(K))
                               GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

                               GG = 1.0_r8/(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                                    (T3D(K)*RR)-1.0_r8))

                               PSI = 2.0_r8/3.0_r8*(ALPHA*DUM/GG)**0.5_r8*AACT

                               ETA1 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW1)
                               ETA2 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW2)

                               SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                               SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                               DUM1 = 1.0_r8/SM1**2*(F11*(PSI/ETA1)**1.5_r8+F21*(SM1**2/(ETA1+3.0_r8*PSI))**0.75_r8)
                               DUM2 = 1.0_r8/SM2**2*(F12*(PSI/ETA2)**1.5_r8+F22*(SM2**2/(ETA2+3.0_r8*PSI))**0.75_r8)

                               SMAX = 1.0_r8/(DUM1+DUM2)**0.5_r8

                               UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                               UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                               DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                               DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                               DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                               ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                               DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                               ! USE LOGNORMAL AEROSOL
                               SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                               AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)

                               SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                               SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                               ! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
                               SMAX = DUM3/QVS(K)

                               UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                               UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                               DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                               DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                               DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                               ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                               DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                               ! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
                               DUM2=MIN(DUM2,DUMACT)

                               DUM2 = (DUM2-NC3D(K))/DT
                               DUM2 = MAX(0.0_r8,DUM2)
                               NC3DTEN(K) = NC3DTEN(K)+DUM2

                            END IF ! IACT
                         END IF ! DUM3/QVS > 1.E-6
                      END IF  ! IDROP = 1

                      !.......................................................................
                   ELSE IF (IBASE.EQ.2) THEN

                      IF (IACT.EQ.1) THEN
                         ! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
                         ! BASED ON TWOMEY 1959

                         DUM=DUM*100.0_r8  ! CONVERT FROM M/S TO CM/S
                         DUM2 = 0.88_r8*C1**(2.0_r8/(K1+2.0_r8))*(7.E-2_r8*DUM**1.5_r8)**(K1/(K1+2.0_r8))
                         DUM2=DUM2*1.E6_r8 ! CONVERT FROM CM-3 TO M-3
                         DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
                         DUM2 = (DUM2-NC3D(K))/DT
                         DUM2 = MAX(0.0_r8,DUM2)
                         NC3DTEN(K) = NC3DTEN(K)+DUM2

                      ELSE IF (IACT.EQ.2) THEN

                         SIGVL = 0.0761_r8-1.55E-4_r8*(T3D(K)-273.15_r8)
                         AACT = 2.0_r8*MW/(RHOW*RR)*SIGVL/T3D(K)
                         ALPHA = SHR_CONST_G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-SHR_CONST_G*MA/(RR*T3D(K))
                         GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

                         GG = 1.0_r8/(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
                              (T3D(K)*RR)-1.0_r8))

                         PSI = 2.0_r8/3.0_r8*(ALPHA*DUM/GG)**0.5_r8*AACT

                         ETA1 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW1)
                         ETA2 = (ALPHA*DUM/GG)**1.5_r8/(2.0_r8*PI*RHOW*GAMM*NANEW2)

                         SM1 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM1))**1.5_r8
                         SM2 = 2.0_r8/BACT**0.5_r8*(AACT/(3.0_r8*RM2))**1.5_r8

                         DUM1 = 1.0_r8/SM1**2*(F11*(PSI/ETA1)**1.5_r8+F21*(SM1**2/(ETA1+3.0_r8*PSI))**0.75_r8)
                         DUM2 = 1.0_r8/SM2**2*(F12*(PSI/ETA2)**1.5_r8+F22*(SM2**2/(ETA2+3.0_r8*PSI))**0.75_r8)

                         SMAX = 1.0_r8/(DUM1+DUM2)**0.5_r8

                         UU1 = 2.0_r8*LOG(SM1/SMAX)/(4.242_r8*LOG(SIG1))
                         UU2 = 2.0_r8*LOG(SM2/SMAX)/(4.242_r8*LOG(SIG2))
                         DUM1 = NANEW1/2.0_r8*(1.0_r8-DERF1(UU1))
                         DUM2 = NANEW2/2.0_r8*(1.0_r8-DERF1(UU2))

                         DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

                         ! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

                         DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

                         DUM2 = (DUM2-NC3D(K))/DT
                         DUM2 = MAX(0.0_r8,DUM2)
                         NC3DTEN(K) = NC3DTEN(K)+DUM2
                      END IF  ! IACT
                   END IF  ! IBASE
                END IF  ! W > 0.001
             END IF  ! QC3D > QSMALL
          END IF  ! INUM = 0
          !!amy to here

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
          ! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
          ! LOSS OF NUMBER CONCENTRATION

          !     IF (PCC(K).LT.0.) THEN
          !        DUM = PCC(K)*DT/QC3D(K)
          !           DUM = MAX(-1.,DUM)
          !        NSUBC(K) = DUM*NC3D(K)/DT
          !     END IF

          IF (EPRD(K).LT.0.0_r8) THEN
             DUM = EPRD(K)*DT/QI3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSUBI(K) = DUM*NI3D(K)/DT
          END IF
          IF (EPRDS(K).LT.0.0_r8) THEN
             DUM = EPRDS(K)*DT/QNI3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSUBS(K) = DUM*NS3D(K)/DT
          END IF
          IF (PRE(K).LT.0.0_r8) THEN
             DUM = PRE(K)*DT/QR3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSUBR(K) = DUM*NR3D(K)/DT
          END IF
          IF (EPRDG(K).LT.0.0_r8) THEN
             DUM = EPRDG(K)*DT/QG3D(K)
             DUM = MAX(-1.0_r8,DUM)
             NSUBG(K) = DUM*NG3D(K)/DT
          END IF

          !        nsubr(k)=0.
          !        nsubs(k)=0.
          !        nsubg(k)=0.

          ! UPDATE TENDENCIES

          !        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
          NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
          NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
          NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
          NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)

       END IF !!!!!! TEMPERATURE

       ! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
       LTRUE = 1

200    CONTINUE

    END DO

    ! INITIALIZE PRECIP AND SNOW RATES
    PRECRT = 0.0_r8
    SNOWRT = 0.0_r8

    ! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

    IF (LTRUE.EQ.0) GOTO 400

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !.......................................................................
    ! CALCULATE SEDIMENATION
    ! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
    ! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
    ! STABILITY, I.E. COURANT# < 1

    !.......................................................................

    NSTEP = 1

    DO K = KTE,KTS,-1

       DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
       DUMQS(K) = QNI3D(K)+QNI3DTEN(K)*DT
       DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
       DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
       DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
       DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
       DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
       DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
       DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
       DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT

       !!amy 
       ! SWITCH FOR CONSTANT DROPLET NUMBER
       IF (INUM.EQ.1) THEN
          DUMFNC(K) = NC3D(K)
       END IF

       ! GET DUMMY LAMDA FOR SEDIMENTATION CALCULATIONS

       ! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
       DUMFNI(K) = MAX(0.0_r8,DUMFNI(K))
       DUMFNS(K) = MAX(0.0_r8,DUMFNS(K))
       DUMFNC(K) = MAX(0.0_r8,DUMFNC(K))
       DUMFNR(K) = MAX(0.0_r8,DUMFNR(K))
       DUMFNG(K) = MAX(0.0_r8,DUMFNG(K))

       !......................................................................
       ! CLOUD ICE

       IF (DUMI(K).GE.QSMALL) THEN
          DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1.0_r8/DI)
          DLAMI=MAX(DLAMI,LAMMINI)
          DLAMI=MIN(DLAMI,LAMMAXI)
       END IF
       !......................................................................
       ! RAIN

       IF (DUMR(K).GE.QSMALL) THEN
          DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1.0_r8/3.0_r8)
          DLAMR=MAX(DLAMR,LAMMINR)
          DLAMR=MIN(DLAMR,LAMMAXR)
       END IF
       !......................................................................
       ! CLOUD DROPLETS

       IF (DUMC(K).GE.QSMALL) THEN
          DUM = PRES(K)/(287.15_r8*T3D(K))
          PGAM(K)=0.0005714_r8*(NC3D(K)/1.E6_r8*DUM)+0.2714_r8
          PGAM(K)=1.0_r8/(PGAM(K)**2)-1.0_r8
          PGAM(K)=MAX(PGAM(K),2.0_r8)
          PGAM(K)=MIN(PGAM(K),10.0_r8)

          DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.0_r8)/(DUMC(K)*GAMMA(PGAM(K)+1.0_r8)))**(1.0_r8/3.0_r8)
          LAMMIN = (PGAM(K)+1.0_r8)/60.E-6_r8
          LAMMAX = (PGAM(K)+1.0_r8)/1.E-6_r8
          DLAMC=MAX(DLAMC,LAMMIN)
          DLAMC=MIN(DLAMC,LAMMAX)
       END IF
       !......................................................................
       ! SNOW

       IF (DUMQS(K).GE.QSMALL) THEN
          DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1.0_r8/DS)
          DLAMS=MAX(DLAMS,LAMMINS)
          DLAMS=MIN(DLAMS,LAMMAXS)
       END IF
       !......................................................................
       ! GRAUPEL

       IF (DUMG(K).GE.QSMALL) THEN
          DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1.0_r8/DG)
          DLAMG=MAX(DLAMG,LAMMING)
          DLAMG=MIN(DLAMG,LAMMAXG)
       END IF

       !......................................................................
       ! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

       ! CLOUD WATER

       IF (DUMC(K).GE.QSMALL) THEN
          UNC =  ACN(K)*GAMMA(1.0_r8+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.0_r8))
          UMC = ACN(K)*GAMMA(4.0_r8+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.0_r8))
       ELSE
          UMC = 0.0_r8
          UNC = 0.0_r8
       END IF

       IF (DUMI(K).GE.QSMALL) THEN
          UNI =  AIN(K)*CONS27/DLAMI**BI
          UMI = AIN(K)*CONS28/(DLAMI**BI)
       ELSE
          UMI = 0.0_r8
          UNI = 0.0_r8
       END IF

       IF (DUMR(K).GE.QSMALL) THEN
          UNR = ARN(K)*CONS6/DLAMR**BR
          UMR = ARN(K)*CONS4/(DLAMR**BR)
       ELSE
          UMR = 0.0_r8
          UNR = 0.0_r8
       END IF

       IF (DUMQS(K).GE.QSMALL) THEN
          UMS = ASN(K)*CONS3/(DLAMS**BS)
          UNS = ASN(K)*CONS5/DLAMS**BS
       ELSE
          UMS = 0.0_r8
          UNS = 0.0_r8
       END IF

       IF (DUMG(K).GE.QSMALL) THEN
          UMG = AGN(K)*CONS7/(DLAMG**BG)
          UNG = AGN(K)*CONS8/DLAMG**BG
       ELSE
          UMG = 0.0_r8
          UNG = 0.0_r8
       END IF

       ! SET REALISTIC LIMITS ON FALLSPEED

       ! bug fix, 10/08/09
       dum=(rhosu/rho(k))**0.54_r8
       UMS=MIN(UMS,1.2_r8*dum)
       UNS=MIN(UNS,1.2_r8*dum)
       ! fix for correction by AA 4/6/11
       UMI=MIN(UMI,1.2_r8*(rhosu/rho(k))**0.35_r8)
       UNI=MIN(UNI,1.2_r8*(rhosu/rho(k))**0.35_r8)
       UMR=MIN(UMR,9.1_r8*dum)
       UNR=MIN(UNR,9.1_r8*dum)
       UMG=MIN(UMG,20.0_r8*dum)
       UNG=MIN(UNG,20.0_r8*dum)

       FR(K) = UMR
       FI(K) = UMI
       FNI(K) = UNI
       FS(K) = UMS
       FNS(K) = UNS
       FNR(K) = UNR
       FC(K) = UMC
       FNC(K) = UNC
       FG(K) = UMG
       FNG(K) = UNG

       ! V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

       IF (K.LE.KTE-1) THEN
          IF (FR(K).LT.1.E-10_r8) THEN
             FR(K)=FR(K+1)
          END IF
          IF (FI(K).LT.1.E-10_r8) THEN
             FI(K)=FI(K+1)
          END IF
          IF (FNI(K).LT.1.E-10_r8) THEN
             FNI(K)=FNI(K+1)
          END IF
          IF (FS(K).LT.1.E-10_r8) THEN
             FS(K)=FS(K+1)
          END IF
          IF (FNS(K).LT.1.E-10_r8) THEN
             FNS(K)=FNS(K+1)
          END IF
          IF (FNR(K).LT.1.E-10_r8) THEN
             FNR(K)=FNR(K+1)
          END IF
          IF (FC(K).LT.1.E-10_r8) THEN
             FC(K)=FC(K+1)
          END IF
          IF (FNC(K).LT.1.E-10_r8) THEN
             FNC(K)=FNC(K+1)
          END IF
          IF (FG(K).LT.1.E-10_r8) THEN
             FG(K)=FG(K+1)
          END IF
          IF (FNG(K).LT.1.E-10_r8) THEN
             FNG(K)=FNG(K+1)
          END IF
       END IF ! K LE KTE-1

       ! CALCULATE NUMBER OF SPLIT TIME STEPS

       RGVM = MAX(FR(K),FI(K),FS(K),FC(K),FNI(K),FNR(K),FNS(K),FNC(K),FG(K),FNG(K))
       ! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
       NSTEP = MAX(INT(RGVM*DT/DZQ(K)+1.0_r8),NSTEP)

       ! MULTIPLY VARIABLES BY RHO
       DUMR(k) = DUMR(k)*RHO(K)
       DUMI(k) = DUMI(k)*RHO(K)
       DUMFNI(k) = DUMFNI(K)*RHO(K)
       DUMQS(k) = DUMQS(K)*RHO(K)
       DUMFNS(k) = DUMFNS(K)*RHO(K)
       DUMFNR(k) = DUMFNR(K)*RHO(K)
       DUMC(k) = DUMC(K)*RHO(K)
       DUMFNC(k) = DUMFNC(K)*RHO(K)
       DUMG(k) = DUMG(K)*RHO(K)
       DUMFNG(k) = DUMFNG(K)*RHO(K)

    END DO

    DO N = 1,NSTEP

       DO K = KTS,KTE
          FALOUTR(K) = FR(K)*DUMR(K)
          FALOUTI(K) = FI(K)*DUMI(K)
          FALOUTNI(K) = FNI(K)*DUMFNI(K)
          FALOUTS(K) = FS(K)*DUMQS(K)
          FALOUTNS(K) = FNS(K)*DUMFNS(K)
          FALOUTNR(K) = FNR(K)*DUMFNR(K)
          FALOUTC(K) = FC(K)*DUMC(K)
          FALOUTNC(K) = FNC(K)*DUMFNC(K)
          FALOUTG(K) = FG(K)*DUMG(K)
          FALOUTNG(K) = FNG(K)*DUMFNG(K)
       END DO

       ! TOP OF MODEL

       K = KTE
       FALTNDR = FALOUTR(K)/DZQ(k)
       FALTNDI = FALOUTI(K)/DZQ(k)
       FALTNDNI = FALOUTNI(K)/DZQ(k)
       FALTNDS = FALOUTS(K)/DZQ(k)
       FALTNDNS = FALOUTNS(K)/DZQ(k)
       FALTNDNR = FALOUTNR(K)/DZQ(k)
       FALTNDC = FALOUTC(K)/DZQ(k)
       FALTNDNC = FALOUTNC(K)/DZQ(k)
       FALTNDG = FALOUTG(K)/DZQ(k)
       FALTNDNG = FALOUTNG(K)/DZQ(k)
       ! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

       QRSTEN(K) = QRSTEN(K)-FALTNDR/NSTEP/RHO(k)
       QISTEN(K) = QISTEN(K)-FALTNDI/NSTEP/RHO(k)
       NI3DTEN(K) = NI3DTEN(K)-FALTNDNI/NSTEP/RHO(k)
       QNISTEN(K) = QNISTEN(K)-FALTNDS/NSTEP/RHO(k)
       NS3DTEN(K) = NS3DTEN(K)-FALTNDNS/NSTEP/RHO(k)
       NR3DTEN(K) = NR3DTEN(K)-FALTNDNR/NSTEP/RHO(k)
       QCSTEN(K) = QCSTEN(K)-FALTNDC/NSTEP/RHO(k)
       NC3DTEN(K) = NC3DTEN(K)-FALTNDNC/NSTEP/RHO(k)
       QGSTEN(K) = QGSTEN(K)-FALTNDG/NSTEP/RHO(k)
       NG3DTEN(K) = NG3DTEN(K)-FALTNDNG/NSTEP/RHO(k)

       DUMR(K) = DUMR(K)-FALTNDR*DT/NSTEP
       DUMI(K) = DUMI(K)-FALTNDI*DT/NSTEP
       DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT/NSTEP
       DUMQS(K) = DUMQS(K)-FALTNDS*DT/NSTEP
       DUMFNS(K) = DUMFNS(K)-FALTNDNS*DT/NSTEP
       DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT/NSTEP
       DUMC(K) = DUMC(K)-FALTNDC*DT/NSTEP
       DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT/NSTEP
       DUMG(K) = DUMG(K)-FALTNDG*DT/NSTEP
       DUMFNG(K) = DUMFNG(K)-FALTNDNG*DT/NSTEP

       DO K = KTE-1,KTS,-1
          FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)
          FALTNDI = (FALOUTI(K+1)-FALOUTI(K))/DZQ(K)
          FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))/DZQ(K)
          FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
          FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
          FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)
          FALTNDC = (FALOUTC(K+1)-FALOUTC(K))/DZQ(K)
          FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))/DZQ(K)
          FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
          FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)

          ! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

          QRSTEN(K) = QRSTEN(K)+FALTNDR/NSTEP/RHO(k)
          QISTEN(K) = QISTEN(K)+FALTNDI/NSTEP/RHO(k)
          NI3DTEN(K) = NI3DTEN(K)+FALTNDNI/NSTEP/RHO(k)
          QNISTEN(K) = QNISTEN(K)+FALTNDS/NSTEP/RHO(k)
          NS3DTEN(K) = NS3DTEN(K)+FALTNDNS/NSTEP/RHO(k)
          NR3DTEN(K) = NR3DTEN(K)+FALTNDNR/NSTEP/RHO(k)
          QCSTEN(K) = QCSTEN(K)+FALTNDC/NSTEP/RHO(k)
          NC3DTEN(K) = NC3DTEN(K)+FALTNDNC/NSTEP/RHO(k)
          QGSTEN(K) = QGSTEN(K)+FALTNDG/NSTEP/RHO(k)
          NG3DTEN(K) = NG3DTEN(K)+FALTNDNG/NSTEP/RHO(k)

          DUMR(K) = DUMR(K)+FALTNDR*DT/NSTEP
          DUMI(K) = DUMI(K)+FALTNDI*DT/NSTEP
          DUMFNI(K) = DUMFNI(K)+FALTNDNI*DT/NSTEP
          DUMQS(K) = DUMQS(K)+FALTNDS*DT/NSTEP
          DUMFNS(K) = DUMFNS(K)+FALTNDNS*DT/NSTEP
          DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT/NSTEP
          DUMC(K) = DUMC(K)+FALTNDC*DT/NSTEP
          DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT/NSTEP
          DUMG(K) = DUMG(K)+FALTNDG*DT/NSTEP
          DUMFNG(K) = DUMFNG(K)+FALTNDNG*DT/NSTEP

       END DO

       ! GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
       ! FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
       ! OF LIQUID WATER CANCELS THIS FACTOR OF 1000

       PRECRT = PRECRT+(FALOUTR(KTS)+FALOUTC(KTS)+FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))  &
            *DT/NSTEP
       SNOWRT = SNOWRT+(FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))*DT/NSTEP

    END DO

    DO K=KTS,KTE

       ! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

       QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
       QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
       QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
       QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
       QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)

       ! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

       !hm 4/7/09 bug fix
       !        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
       IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15_r8.AND.LAMI(K).GE.1.E-10_r8) THEN
          IF (1.0_r8/LAMI(K).GE.2.0_r8*DCS) THEN
             QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
             NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
             QI3DTEN(K) = -QI3D(K)/DT
             NI3DTEN(K) = -NI3D(K)/DT
          END IF
       END IF

       ! hm add tendencies here, then call sizeparameter
       ! to ensure consisitency between mixing ratio and number concentration

       QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
       QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
       QNI3D(k)       = QNI3D(k)+QNI3DTEN(k)*DT
       QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
       NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
       NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
       NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
       NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT

       IF (IGRAUP.EQ.0) THEN
          QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
          NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
       END IF

       ! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
       T3D(K)         = T3D(K)+T3DTEN(k)*DT
       QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT

       ! SATURATION VAPOR PRESSURE AND MIXING RATIO

       EVS(K) = POLYSVP(T3D(K),0)   ! PA
       EIS(K) = POLYSVP(T3D(K),1)   ! PA

       ! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

       IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

       QVS(K) = 0.622_r8*EVS(K)/(PRES(K)-EVS(K))
       QVI(K) = 0.622_r8*EIS(K)/(PRES(K)-EIS(K))

       QVQVS(K) = QV3D(K)/QVS(K)
       QVQVSI(K) = QV3D(K)/QVI(K)

       ! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
       ! hm 7/9/09 change limit to 1.e-8

       IF (QVQVS(K).LT.0.9_r8) THEN
          IF (QR3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QR3D(K)
             T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
             QR3D(K)=0.0_r8
          END IF
          IF (QC3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QC3D(K)
             T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
             QC3D(K)=0.0_r8
          END IF
       END IF

       IF (QVQVSI(K).LT.0.9_r8) THEN
          IF (QI3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QI3D(K)
             T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
             QI3D(K)=0.0_r8
          END IF
          IF (QNI3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QNI3D(K)
             T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
             QNI3D(K)=0.0_r8
          END IF
          IF (QG3D(K).LT.1.E-8_r8) THEN
             QV3D(K)=QV3D(K)+QG3D(K)
             T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
             QG3D(K)=0.0_r8
          END IF
       END IF

       !..................................................................
       ! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
          QC3D(K) = 0.0_r8
          NC3D(K) = 0.0_r8
          EFFC(K) = 0.0_r8
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
          QR3D(K) = 0.0_r8
          NR3D(K) = 0.0_r8
          EFFR(K) = 0.0_r8
       END IF
       IF (QI3D(K).LT.QSMALL) THEN
          QI3D(K) = 0.0_r8
          NI3D(K) = 0.0_r8
          EFFI(K) = 0.0_r8
       END IF
       IF (QNI3D(K).LT.QSMALL) THEN
          QNI3D(K) = 0.0_r8
          NS3D(K) = 0.0_r8
          EFFS(K) = 0.0_r8
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
          QG3D(K) = 0.0_r8
          NG3D(K) = 0.0_r8
          EFFG(K) = 0.0_r8
       END IF

       !..................................
       ! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

       IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
            .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) GOTO 500

       !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ! CALCULATE INSTANTANEOUS PROCESSES

       ! ADD MELTING OF CLOUD ICE TO FORM RAIN

       IF (QI3D(K).GE.QSMALL.AND.T3D(K).GE.273.15_r8) THEN
          QR3D(K) = QR3D(K)+QI3D(K)
          T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
          QI3D(K) = 0.0_r8
          NR3D(K) = NR3D(K)+NI3D(K)
          NI3D(K) = 0.0_r8
       END IF

       ! ****SENSITIVITY - NO ICE
       IF (ILIQ.EQ.1) GOTO 778

       ! HOMOGENEOUS FREEZING OF CLOUD WATER

       IF (T3D(K).LE.233.15_r8.AND.QC3D(K).GE.QSMALL) THEN
          QI3D(K)=QI3D(K)+QC3D(K)
          T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
          QC3D(K)=0.0_r8
          NI3D(K)=NI3D(K)+NC3D(K)
          NC3D(K)=0.0_r8
       END IF

       ! HOMOGENEOUS FREEZING OF RAIN

       IF (IGRAUP.EQ.0) THEN

          IF (T3D(K).LE.233.15_r8.AND.QR3D(K).GE.QSMALL) THEN
             QG3D(K) = QG3D(K)+QR3D(K)
             T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
             QR3D(K) = 0.0_r8
             NG3D(K) = NG3D(K)+ NR3D(K)
             NR3D(K) = 0.0_r8
          END IF

       ELSE IF (IGRAUP.EQ.1) THEN

          IF (T3D(K).LE.233.15_r8.AND.QR3D(K).GE.QSMALL) THEN
             QNI3D(K) = QNI3D(K)+QR3D(K)
             T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
             QR3D(K) = 0.0_r8
             NS3D(K) = NS3D(K)+NR3D(K)
             NR3D(K) = 0.0_r8
          END IF

       END IF

778    CONTINUE

       ! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

       NI3D(K) = MAX(0.0_r8,NI3D(K))
       NS3D(K) = MAX(0.0_r8,NS3D(K))
       NC3D(K) = MAX(0.0_r8,NC3D(K))
       NR3D(K) = MAX(0.0_r8,NR3D(K))
       NG3D(K) = MAX(0.0_r8,NG3D(K))

       !......................................................................
       ! CLOUD ICE

       IF (QI3D(K).GE.QSMALL) THEN
          LAMI(K) = (CONS12*                 &
               NI3D(K)/QI3D(K))**(1.0_r8/DI)

          ! CHECK FOR SLOPE

          ! ADJUST VARS

          IF (LAMI(K).LT.LAMMINI) THEN

             LAMI(K) = LAMMINI

             N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

             NI3D(K) = N0I(K)/LAMI(K)
          ELSE IF (LAMI(K).GT.LAMMAXI) THEN
             LAMI(K) = LAMMAXI
             N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

             NI3D(K) = N0I(K)/LAMI(K)
          END IF
       END IF

       !......................................................................
       ! RAIN

       IF (QR3D(K).GE.QSMALL) THEN
          LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1.0_r8/3.0_r8)

          ! CHECK FOR SLOPE

          ! ADJUST VARS

          IF (LAMR(K).LT.LAMMINR) THEN

             LAMR(K) = LAMMINR

             N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

             NR3D(K) = N0RR(K)/LAMR(K)
          ELSE IF (LAMR(K).GT.LAMMAXR) THEN
             LAMR(K) = LAMMAXR
             N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

             NR3D(K) = N0RR(K)/LAMR(K)
          END IF

       END IF

       !......................................................................
       ! CLOUD DROPLETS

       ! MARTIN ET AL. (1994) FORMULA FOR PGAM

       IF (QC3D(K).GE.QSMALL) THEN

          DUM = PRES(K)/(287.15_r8*T3D(K))
          PGAM(K)=0.0005714_r8*(NC3D(K)/1.E6_r8*DUM)+0.2714_r8
          PGAM(K)=1.0_r8/(PGAM(K)**2)-1.0_r8
          PGAM(K)=MAX(PGAM(K),2.0_r8)
          PGAM(K)=MIN(PGAM(K),10.0_r8)

          ! CALCULATE LAMC

          LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.0_r8)/   &
               (QC3D(K)*GAMMA(PGAM(K)+1.0_r8)))**(1.0_r8/3.0_r8)

          ! LAMMIN, 60 MICRON DIAMETER
          ! LAMMAX, 1 MICRON

          LAMMIN = (PGAM(K)+1.0_r8)/60.E-6_r8
          LAMMAX = (PGAM(K)+1.0_r8)/1.E-6_r8

          IF (LAMC(K).LT.LAMMIN) THEN
             LAMC(K) = LAMMIN
             NC3D(K) = EXP(3.0_r8*LOG(LAMC(K))+LOG(QC3D(K))+              &
                  LOG(GAMMA(PGAM(K)+1.0_r8))-LOG(GAMMA(PGAM(K)+4.0_r8)))/CONS26

          ELSE IF (LAMC(K).GT.LAMMAX) THEN
             LAMC(K) = LAMMAX
             NC3D(K) = EXP(3.0_r8*LOG(LAMC(K))+LOG(QC3D(K))+              &
                  LOG(GAMMA(PGAM(K)+1.0_r8))-LOG(GAMMA(PGAM(K)+4.0_r8)))/CONS26

          END IF

       END IF

       !......................................................................
       ! SNOW

       IF (QNI3D(K).GE.QSMALL) THEN
          LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1.0_r8/DS)

          ! CHECK FOR SLOPE

          ! ADJUST VARS

          IF (LAMS(K).LT.LAMMINS) THEN
             LAMS(K) = LAMMINS
             N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

             NS3D(K) = N0S(K)/LAMS(K)

          ELSE IF (LAMS(K).GT.LAMMAXS) THEN

             LAMS(K) = LAMMAXS
             N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1
             NS3D(K) = N0S(K)/LAMS(K)
          END IF

       END IF

       !......................................................................
       ! GRAUPEL

       IF (QG3D(K).GE.QSMALL) THEN
          LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1.0_r8/DG)

          ! CHECK FOR SLOPE

          ! ADJUST VARS

          IF (LAMG(K).LT.LAMMING) THEN
             LAMG(K) = LAMMING
             N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

             NG3D(K) = N0G(K)/LAMG(K)

          ELSE IF (LAMG(K).GT.LAMMAXG) THEN

             LAMG(K) = LAMMAXG
             N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

             NG3D(K) = N0G(K)/LAMG(K)
          END IF

       END IF

500    CONTINUE

       ! CALCULATE EFFECTIVE RADIUS

       IF (QI3D(K).GE.QSMALL) THEN
          EFFI(K) = 3.0_r8/LAMI(K)/2.0_r8*1.E6_r8
       ELSE
          EFFI(K) = 25.0_r8
       END IF

       IF (QNI3D(K).GE.QSMALL) THEN
          EFFS(K) = 3.0_r8/LAMS(K)/2.0_r8*1.E6_r8
       ELSE
          EFFS(K) = 25.0_r8
       END IF

       IF (QR3D(K).GE.QSMALL) THEN
          EFFR(K) = 3.0_r8/LAMR(K)/2.0_r8*1.E6_r8
       ELSE
          EFFR(K) = 25.0_r8
       END IF

       IF (QC3D(K).GE.QSMALL) THEN
          EFFC(K) = GAMMA(PGAM(K)+4.0_r8)/                        &
               GAMMA(PGAM(K)+3.0_r8)/LAMC(K)/2.0_r8*1.E6_r8
       ELSE
          EFFC(K) = 25.0_r8
       END IF

       IF (QG3D(K).GE.QSMALL) THEN
          EFFG(K) = 3.0_r8/LAMG(K)/2.0_r8*1.E6_r8
       ELSE
          EFFG(K) = 25.0_r8
       END IF

       ! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
       ! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
       ! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
       !          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
       ! HM, 3/4/13, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
       ! OF EXCESSIVE AND PERSISTENT ANVIL
       ! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
       NI3D(K) = MIN(NI3D(K),0.3E6_r8/RHO(K))

       ! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
       IF (INUM.EQ.0.AND.IACT.EQ.2) THEN
          NC3D(K) = MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K))
       END IF
       !!amy 
       ! SWITCH FOR CONSTANT DROPLET NUMBER
       IF (INUM.EQ.1) THEN
          ! CHANGE NDCNST FROM CM-3 TO KG-1
          NC3D(K) = NDCNST*1.E6_r8/RHO(K)
       END IF

    END DO !!! K LOOP

400 CONTINUE

    ! ALL DONE !!!!!!!!!!!
    RETURN
  END SUBROUTINE MORR_TWO_MOMENT_MICRO
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  REAL(KIND=r8) FUNCTION POLYSVP (T,TYPE)

    !-------------------------------------------

    !  COMPUTE SATURATION VAPOR PRESSURE

    !  POLYSVP RETURNED IN UNITS OF PA.
    !  T IS INPUT IN UNITS OF K.
    !  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

    ! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM MARAT KHROUTDINOV

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN   ) :: T
    INTEGER      , INTENT(IN   ) :: TYPE

    ! ice
    REAL(KIND=r8) :: a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
    DATA a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
         6.11147274_r8    , 0.503160820_r8    , 0.188439774e-1_r8, &
         0.420895665e-3_r8, 0.615021634e-5_r8 , 0.602588177e-7_r8, &
         0.385852041e-9_r8, 0.146898966e-11_r8, 0.252751365e-14_r8/

    ! liquid
    REAL(KIND=r8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8 

    ! V1.7
    DATA a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
         6.11239921_r8     , 0.443987641_r8    , 0.142986287e-1_r8, &
         0.264847430e-3_r8 , 0.302950461e-5_r8 , 0.206739458e-7_r8, &
         0.640689451e-10_r8,-0.952447341e-13_r8,-0.976195544e-15_r8/
    REAL(KIND=r8) :: dt

    ! ICE

    IF (TYPE.EQ.1) THEN

       !         POLYSVP = 10.0_r8**(-9.09718_r8*(273.16_r8/T-1.0_r8)-3.56654_r8*&
       !          LOG10(273.16_r8/T)+0.876793_r8*(1.0_r8-T/273.16_r8)+&
       !          LOG10(6.1071_r8))*100.0_r8


       dt = MAX(-80.0_r8,t-273.16_r8)
       polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
       polysvp = polysvp*100.0_r8

    END IF

    ! LIQUID

    IF (TYPE.EQ.0) THEN

       dt = MAX(-80.0_r8,t-273.16_r8)
       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp = polysvp*100.0_r8

       !         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
       !             5.02808*LOG10(373.16/T)-&
       !             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+&
       !             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+&
       !             LOG10(1013.246))*100.

    END IF


  END FUNCTION POLYSVP


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  REAL(KIND=r8) FUNCTION GAMMA(X)
    !----------------------------------------------------------------------
    !
    ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
    !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
    !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
    !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
    !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
    !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
    !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
    !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
    !   MACHINE-DEPENDENT CONSTANTS.
    !
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
    !
    ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
    ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
    ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
    !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
    !                  GAMMA(XBIG) = BETA**MAXEXP
    ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
    !          APPROXIMATELY BETA**MAXEXP
    ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    !          1.0+EPS .GT. 1.0
    ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    !          1/XMININ IS MACHINE REPRESENTABLE
    !
    !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
    !
    !                            BETA       MAXEXP        XBIG
    !
    ! CRAY-1         (S.P.)        2         8191        966.961
    ! CYBER 180/855
    !   UNDER NOS    (S.P.)        2         1070        177.803
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (S.P.)        2          128        35.040
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (D.P.)        2         1024        171.624
    ! IBM 3033       (D.P.)       16           63        57.574
    ! VAX D-FORMAT   (D.P.)        2          127        34.844
    ! VAX G-FORMAT   (D.P.)        2         1023        171.489
    !
    !                            XINF         EPS        XMININ
    !
    ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
    ! CYBER 180/855
    !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
    ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
    ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
    ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! ERROR RETURNS
    !
    !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
    !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
    !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
    !
    !
    !  INTRINSIC FUNCTIONS REQUIRED ARE:
    !
    !     INT, DBLE, EXP, LOG, REAL, SIN
    !
    !
    ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
    !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
    !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
    !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
    !
    !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
    !              SONS, NEW YORK, 1968.
    !
    !  LATEST MODIFICATION: OCTOBER 12, 1989
    !
    !  AUTHORS: W. J. CODY AND L. STOLTZ
    !           APPLIED MATHEMATICS DIVISION
    !           ARGONNE NATIONAL LABORATORY
    !           ARGONNE, IL 60439
    !
    !----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: X  
    INTEGER I,N
    LOGICAL PARITY
    REAL(KIND=r8)                                                          &
         CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
         TWO,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
    REAL(KIND=r8), DIMENSION(7) :: C
    REAL(KIND=r8), DIMENSION(8) :: P
    REAL(KIND=r8), DIMENSION(8) :: Q
    !----------------------------------------------------------------------
    !  MATHEMATICAL CONSTANTS
    !----------------------------------------------------------------------
    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/


    !----------------------------------------------------------------------
    !  MACHINE DEPENDENT PARAMETERS
    !----------------------------------------------------------------------
    DATA XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/,XINF/3.4E38_r8/
    !----------------------------------------------------------------------
    !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
    !     APPROXIMATION OVER (1,2).
    !----------------------------------------------------------------------
    DATA P/-1.71618513886549492533811E+0_r8, 2.47656508055759199108314E+1_r8,  &
         -3.79804256470945635097577E+2_r8, 6.29331155312818442661052E+2_r8,  &
         8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,  &
         -3.61444134186911729807069E+4_r8, 6.64561438202405440627855E+4_r8/
    DATA Q/-3.08402300119738975254353E+1_r8, 3.15350626979604161529144E+2_r8,  &
         -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,  &
         2.25381184209801510330112E+4_r8, 4.75584627752788110767815E+3_r8,  &
         -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
    !----------------------------------------------------------------------
    !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
    !----------------------------------------------------------------------
    DATA C/-1.910444077728E-03_r8,         8.4171387781295E-04_r8,       &
         -5.952379913043012E-04_r8,      7.93650793500350248E-04_r8,   &
         -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8, &
         5.7083835261E-03_r8/
    !----------------------------------------------------------------------
    !  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
    !----------------------------------------------------------------------
    CONV(I) = REAL(I,KIND=r8)
    PARITY=.FALSE.
    FACT=ONE
    N=0
    Y=X
    IF(Y.LE.ZERO)THEN
       !----------------------------------------------------------------------
       !  ARGUMENT IS NEGATIVE
       !----------------------------------------------------------------------
       Y=-X
       Y1=AINT(Y)
       RES=Y-Y1
       IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
       ELSE
          RES=XINF
          GOTO 900
       ENDIF
    ENDIF
    !----------------------------------------------------------------------
    !  ARGUMENT IS POSITIVE
    !----------------------------------------------------------------------
    IF(Y.LT.EPS)THEN
       !----------------------------------------------------------------------
       !  ARGUMENT .LT. EPS
       !----------------------------------------------------------------------
       IF(Y.GE.XMININ)THEN
          RES=ONE/Y
       ELSE
          RES=XINF
          GOTO 900
       ENDIF
    ELSEIF(Y.LT.TWELVE)THEN
       Y1=Y
       IF(Y.LT.ONE)THEN
          !----------------------------------------------------------------------
          !  0.0 .LT. ARGUMENT .LT. 1.0
          !----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
       ELSE
          !----------------------------------------------------------------------
          !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
          !----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
       ENDIF
       !----------------------------------------------------------------------
       !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
       !----------------------------------------------------------------------
       XNUM=ZERO
       XDEN=ONE
       DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
       END DO
       RES=XNUM/XDEN+ONE
       IF(Y1.LT.Y)THEN
          !----------------------------------------------------------------------
          !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
          !----------------------------------------------------------------------
          RES=RES/Y1
       ELSEIF(Y1.GT.Y)THEN
          !----------------------------------------------------------------------
          !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
          !----------------------------------------------------------------------
          DO I=1,N
             RES=RES*Y
             Y=Y+ONE
          END DO
       ENDIF
    ELSE
       !----------------------------------------------------------------------
       !  EVALUATE FOR ARGUMENT .GE. 12.0,
       !----------------------------------------------------------------------
       IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
             SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
       ELSE
          RES=XINF
          GOTO 900
       ENDIF
    ENDIF
    !----------------------------------------------------------------------
    !  FINAL ADJUSTMENTS AND RETURN
    !----------------------------------------------------------------------
    IF(PARITY)RES=-RES
    IF(FACT.NE.ONE)RES=FACT/RES
900 GAMMA=RES
    RETURN
    ! ---------- LAST LINE OF GAMMA ----------
  END FUNCTION GAMMA

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------


  REAL(KIND=r8) FUNCTION DERF1(X)
    IMPLICIT NONE
    REAL(KIND=r8) ::  X
    REAL(KIND=r8) :: A(0 : 64), B(0 : 64)
    REAL(KIND=r8) :: W,T,Y
    INTEGER K,I
    DATA A/                                                 &
         0.00000000005958930743E0_r8, -0.00000000113739022964E0_r8, &
         0.00000001466005199839E0_r8, -0.00000016350354461960E0_r8, &
         0.00000164610044809620E0_r8, -0.00001492559551950604E0_r8, &
         0.00012055331122299265E0_r8, -0.00085483269811296660E0_r8, &
         0.00522397762482322257E0_r8, -0.02686617064507733420E0_r8, &
         0.11283791670954881569E0_r8, -0.37612638903183748117E0_r8, &
         1.12837916709551257377E0_r8,                               &
         0.00000000002372510631E0_r8, -0.00000000045493253732E0_r8, &
         0.00000000590362766598E0_r8, -0.00000006642090827576E0_r8, &
         0.00000067595634268133E0_r8, -0.00000621188515924000E0_r8, &
         0.00005103883009709690E0_r8, -0.00037015410692956173E0_r8, &
         0.00233307631218880978E0_r8, -0.01254988477182192210E0_r8, &
         0.05657061146827041994E0_r8, -0.21379664776456006580E0_r8, &
         0.84270079294971486929E0_r8,                               &
         0.00000000000949905026E0_r8, -0.00000000018310229805E0_r8, &
         0.00000000239463074000E0_r8, -0.00000002721444369609E0_r8, &
         0.00000028045522331686E0_r8, -0.00000261830022482897E0_r8, &
         0.00002195455056768781E0_r8, -0.00016358986921372656E0_r8, &
         0.00107052153564110318E0_r8, -0.00608284718113590151E0_r8, &
         0.02986978465246258244E0_r8, -0.13055593046562267625E0_r8, &
         0.67493323603965504676E0_r8,                               &
         0.00000000000382722073E0_r8, -0.00000000007421598602E0_r8, &
         0.00000000097930574080E0_r8, -0.00000001126008898854E0_r8, &
         0.00000011775134830784E0_r8, -0.00000111992758382650E0_r8, &
         0.00000962023443095201E0_r8, -0.00007404402135070773E0_r8, &
         0.00050689993654144881E0_r8, -0.00307553051439272889E0_r8, &
         0.01668977892553165586E0_r8, -0.08548534594781312114E0_r8, &
         0.56909076642393639985E0_r8,                               &
         0.00000000000155296588E0_r8, -0.00000000003032205868E0_r8, &
         0.00000000040424830707E0_r8, -0.00000000471135111493E0_r8, &
         0.00000005011915876293E0_r8, -0.00000048722516178974E0_r8, &
         0.00000430683284629395E0_r8, -0.00003445026145385764E0_r8, &
         0.00024879276133931664E0_r8, -0.00162940941748079288E0_r8, &
         0.00988786373932350462E0_r8, -0.05962426839442303805E0_r8, &
         0.49766113250947636708E0_r8 /
    DATA (B(I), I = 0, 12) /                                         &
         -0.00000000029734388465E0_r8,  0.00000000269776334046E0_r8, &
         -0.00000000640788827665E0_r8, -0.00000001667820132100E0_r8, &
         -0.00000021854388148686E0_r8,  0.00000266246030457984E0_r8, &
         0.00001612722157047886E0_r8, -0.00025616361025506629E0_r8, &
         0.00015380842432375365E0_r8,  0.00815533022524927908E0_r8, &
         -0.01402283663896319337E0_r8, -0.19746892495383021487E0_r8, &
         0.71511720328842845913E0_r8 /
    DATA (B(I), I = 13, 25) /                                        &
         -0.00000000001951073787E0_r8, -0.00000000032302692214E0_r8, &
         0.00000000522461866919E0_r8,  0.00000000342940918551E0_r8, &
         -0.00000035772874310272E0_r8,  0.00000019999935792654E0_r8, &
         0.00002687044575042908E0_r8, -0.00011843240273775776E0_r8, &
         -0.00080991728956032271E0_r8,  0.00661062970502241174E0_r8, &
         0.00909530922354827295E0_r8, -0.20160072778491013140E0_r8, &
         0.51169696718727644908E0_r8 /
    DATA (B(I), I = 26, 38) /                                        &
         0.00000000003147682272E0_r8, -0.00000000048465972408E0_r8, &
         0.00000000063675740242E0_r8,  0.00000003377623323271E0_r8, &
         -0.00000015451139637086E0_r8, -0.00000203340624738438E0_r8, &
         0.00001947204525295057E0_r8,  0.00002854147231653228E0_r8, &
         -0.00101565063152200272E0_r8,  0.00271187003520095655E0_r8, &
         0.02328095035422810727E0_r8, -0.16725021123116877197E0_r8, &
         0.32490054966649436974E0_r8 /
    DATA (B(I), I = 39, 51) /                                        &
         0.00000000002319363370E0_r8, -0.00000000006303206648E0_r8, &
         -0.00000000264888267434E0_r8,  0.00000002050708040581E0_r8, &
         0.00000011371857327578E0_r8, -0.00000211211337219663E0_r8, &
         0.00000368797328322935E0_r8,  0.00009823686253424796E0_r8, &
         -0.00065860243990455368E0_r8, -0.00075285814895230877E0_r8, &
         0.02585434424202960464E0_r8, -0.11637092784486193258E0_r8, &
         0.18267336775296612024E0_r8 /
    DATA (B(I), I = 52, 64) /                                        &
         -0.00000000000367789363E0_r8,  0.00000000020876046746E0_r8, &
         -0.00000000193319027226E0_r8, -0.00000000435953392472E0_r8, &
         0.00000018006992266137E0_r8, -0.00000078441223763969E0_r8, &
         -0.00000675407647949153E0_r8,  0.00008428418334440096E0_r8, &
         -0.00017604388937031815E0_r8, -0.00239729611435071610E0_r8, &
         0.02064129023876022970E0_r8, -0.06905562880005864105E0_r8, &
         0.09084526782065478489E0_r8 /
    W = ABS(X)
    IF (W .LT. 2.20_r8) THEN
       T = W * W
       K = INT(T)
       T = T - K
       K = K * 13
       Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
            A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
            A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
            A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T +   &
            A(K + 11)) * T + A(K + 12)) * W
    ELSE IF (W .LT. 6.90_r8) THEN
       K = INT(W)
       T = W - K
       K = 13 * (K - 2)
       Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
            B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T +   &
            B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T +   &
            B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T +   &
            B(K + 11)) * T + B(K + 12)
       Y = Y * Y
       Y = Y * Y
       Y = Y * Y
       Y = 1 - Y * Y
    ELSE
       Y = 1
    END IF
    IF (X .LT. 00_r8) Y = -Y
    DERF1 = Y
  END FUNCTION DERF1

  !+---+-----------------------------------------------------------------+

  !+---+-----------------------------------------------------------------+

  SUBROUTINE refl10cm_hm (qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d, &
       t1d, p1d, dBZ,  kMax)

    IMPLICIT NONE

    !..Sub arguments
    INTEGER, INTENT(IN):: kMax
    REAL(KIND=r8), INTENT(IN):: qv1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: qr1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: nr1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: qs1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: ns1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: qg1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: ng1d(1:kMax)
    REAL(KIND=r8), INTENT(IN):: t1d (1:kMax)
    REAL(KIND=r8), INTENT(IN):: p1d (1:kMax)

    REAL(KIND=r8), INTENT(INOUT):: dBZ(1:kMax)

    !..Local variables
    REAL(KIND=r8) :: temp (1:kMax)
    REAL(KIND=r8) :: pres (1:kMax)
    REAL(KIND=r8) :: qv   (1:kMax)
    REAL(KIND=r8) :: rho  (1:kMax)
    REAL(KIND=r8) :: rr   (1:kMax)
    REAL(KIND=r8) :: nr   (1:kMax)
    REAL(KIND=r8) :: rs   (1:kMax)
    REAL(KIND=r8) :: ns   (1:kMax)
    REAL(KIND=r8) :: rg   (1:kMax)
    REAL(KIND=r8) :: ng   (1:kMax)

    REAL(KIND=r8) :: ilamr(1:kMax)
    REAL(KIND=r8) :: ilamg(1:kMax)
    REAL(KIND=r8) :: ilams(1:kMax)
    REAL(KIND=r8) :: N0_r (1:kMax)
    REAL(KIND=r8) :: N0_g (1:kMax)
    REAL(KIND=r8) :: N0_s (1:kMax)
    REAL(KIND=r8) :: lamr, lamg, lams
    LOGICAL       :: L_qr(1:kMax)
    LOGICAL       :: L_qs(1:kMax)
    LOGICAL       :: L_qg(1:kMax)

    REAL(KIND=r8) :: ze_rain   (1:kMax)
    REAL(KIND=r8) :: ze_snow   (1:kMax)
    REAL(KIND=r8) :: ze_graupel(1:kMax)
    REAL(KIND=r8) :: fmelt_s, fmelt_g
    REAL(KIND=r8) :: cback, x, eta, f_d

    INTEGER:: k, k_0,  n
    LOGICAL:: melti

    !+---+
    temp = 0.0_r8;   pres = 0.0_r8;   qv = 0.0_r8;   rho  = 0.0_r8;
    rr = 0.0_r8;   nr = 0.0_r8;   rs = 0.0_r8;   ns = 0.0_r8;
    rg = 0.0_r8;   ng = 0.0_r8;   ilamr = 0.0_r8;   ilamg = 0.0_r8;
    ilams = 0.0_r8;   N0_r = 0.0_r8;   N0_g = 0.0_r8;   N0_s = 0.0_r8;
    lamr = 0.0_r8; lamg = 0.0_r8; lams = 0.0_r8;
    ze_rain   = 0.0_r8;   ze_snow   = 0.0_r8;
    ze_graupel = 0.0_r8;   fmelt_s= 0.0_r8; fmelt_g= 0.0_r8;
    cback= 0.0_r8; x= 0.0_r8; eta= 0.0_r8; f_d= 0.0_r8;

    DO k = 1, kMax
       dBZ(k) = -35.0_r8
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Put column of data into local arrays.
    !+---+-----------------------------------------------------------------+
    DO k = 1, kMax
       temp(k) = t1d(k)
       qv(k) = MAX(1.E-10_r8, qv1d(k))
       pres(k) = p1d(k)
       rho(k) = 0.622_r8*pres(k)/(R_d*temp(k)*(qv(k)+0.622_r8))

       IF (qr1d(k) .GT. 1.E-9_r8) THEN
          rr(k) = qr1d(k)*rho(k)
          nr(k) = nr1d(k)*rho(k)
          lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
          ilamr(k) = 1.0_r8/lamr
          N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
          L_qr(k) = .TRUE.
       ELSE
          rr(k) = 1.E-12_r8
          nr(k) = 1.E-12_r8
          L_qr(k) = .FALSE.
       ENDIF

       IF (qs1d(k) .GT. 1.E-9_r8) THEN
          rs(k) = qs1d(k)*rho(k)
          ns(k) = ns1d(k)*rho(k)
          lams = (xam_s*xcsg(3)*xosg2*ns(k)/rs(k))**xobms
          ilams(k) = 1.0_r8/lams
          N0_s(k) = ns(k)*xosg2*lams**xcse(2)
          L_qs(k) = .TRUE.
       ELSE
          rs(k) = 1.E-12_r8
          ns(k) = 1.E-12_r8
          L_qs(k) = .FALSE.
       ENDIF

       IF (qg1d(k) .GT. 1.E-9_r8) THEN
          rg(k) = qg1d(k)*rho(k)
          ng(k) = ng1d(k)*rho(k)
          lamg = (xam_g*xcgg(3)*xogg2*ng(k)/rg(k))**xobmg
          ilamg(k) = 1.0_r8/lamg
          N0_g(k) = ng(k)*xogg2*lamg**xcge(2)
          L_qg(k) = .TRUE.
       ELSE
          rg(k) = 1.E-12_r8
          ng(k) = 1.E-12_r8
          L_qg(k) = .FALSE.
       ENDIF
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Locate K-level of start of melting (k_0 is level above).
    !+---+-----------------------------------------------------------------+
    melti = .FALSE.
    k_0 = 1
    DO k = kMax-1, 1, -1
       IF ( (temp(k).GT.273.15_r8) .AND. L_qr(k)                         &
            .AND. (L_qs(k+1).OR.L_qg(k+1)) ) THEN
          k_0 = MAX(k+1, k_0)
          melti=.TRUE.
          GOTO 195
       ENDIF
    ENDDO
195 CONTINUE

    !+---+-----------------------------------------------------------------+
    !..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
    !.. and non-water-coated snow and graupel when below freezing are
    !.. simple. Integrations of m(D)*m(D)*N(D)*dD.
    !+---+-----------------------------------------------------------------+

    DO k = 1, kMax
       ze_rain(k) = 1.e-22_r8
       ze_snow(k) = 1.e-22_r8
       ze_graupel(k) = 1.e-22_r8
       IF (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
       IF (L_qs(k)) ze_snow(k) = (0.176_r8/0.93_r8) * (6.0_r8/PI)*(6.0_r8/PI)     &
            * (xam_s/900.0_r8)*(xam_s/900.0_r8)          &
            * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
       IF (L_qg(k)) ze_graupel(k) = (0.176_r8/0.93_r8) * (6.0_r8/PI)*(6.0_r8/PI)  &
            * (xam_g/900.0_r8)*(xam_g/900.0_r8)       &
            * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
    ENDDO

    !+---+-----------------------------------------------------------------+
    !..Special case of melting ice (snow/graupel) particles.  Assume the
    !.. ice is surrounded by the liquid water.  Fraction of meltwater is
    !.. extremely simple based on amount found above the melting level.
    !.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
    !.. routines).
    !+---+-----------------------------------------------------------------+

    IF (melti .AND. k_0.GE.1+1) THEN
       DO k = k_0-1, 1, -1

          !..Reflectivity contributed by melting snow
          IF (L_qs(k) .AND. L_qs(k_0) ) THEN
             fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
             eta = 0.d0
             lams = 1./ilams(k)
             DO n = 1, nrbins
                x = xam_s * xxDs(n)**xbm_s
                CALL rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                     fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                     CBACK, mixingrulestring_s, matrixstring_s,          &
                     inclusionstring_s, hoststring_s,                    &
                     hostmatrixstring_s, hostinclusionstring_s)
                f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
                eta = eta + f_d * CBACK * simpson(n) * xdts(n)
             ENDDO
             ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          ENDIF


          !..Reflectivity contributed by melting graupel

          IF (L_qg(k) .AND. L_qg(k_0) ) THEN
             fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
             eta = 0.d0
             lamg = 1./ilamg(k)
             DO n = 1, nrbins
                x = xam_g * xxDg(n)**xbm_g
                CALL rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
                     fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                     CBACK, mixingrulestring_g, matrixstring_g,          &
                     inclusionstring_g, hoststring_g,                    &
                     hostmatrixstring_g, hostinclusionstring_g)
                f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
                eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
             ENDDO
             ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          ENDIF

       ENDDO
    ENDIF

    DO k = kMax, 1, -1
       dBZ(k) = 10.0_r8*LOG10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
    ENDDO


  END SUBROUTINE refl10cm_hm

  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+
  !+---+-----------------------------------------------------------------+

  SUBROUTINE radar_init()

    IMPLICIT NONE
    INTEGER:: n
    PI5    = 3.14159_r8*3.14159_r8*3.14159_r8*3.14159_r8*3.14159_r8
    lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
    m_w_0  = m_complex_water_ray (lamda_radar, 0.0d0)
    m_i_0  = m_complex_ice_maetzler (lamda_radar, 0.0d0)
    K_w    = (ABS( (m_w_0*m_w_0 - 1.0_r8) /(m_w_0*m_w_0 + 2.0_r8) ))**2

    DO n = 1, nrbins+1
       simpson(n) = 0.0d0
    ENDDO
    DO n = 1, nrbins-1, 2
       simpson(n) = simpson(n) + basis(1)
       simpson(n+1) = simpson(n+1) + basis(2)
       simpson(n+2) = simpson(n+2) + basis(3)
    ENDDO

    DO n = 1, slen
       mixingrulestring_s(n:n) = CHAR(0)
       matrixstring_s(n:n) = CHAR(0)
       inclusionstring_s(n:n) = CHAR(0)
       hoststring_s(n:n) = CHAR(0)
       hostmatrixstring_s(n:n) = CHAR(0)
       hostinclusionstring_s(n:n) = CHAR(0)
       mixingrulestring_g(n:n) = CHAR(0)
       matrixstring_g(n:n) = CHAR(0)
       inclusionstring_g(n:n) = CHAR(0)
       hoststring_g(n:n) = CHAR(0)
       hostmatrixstring_g(n:n) = CHAR(0)
       hostinclusionstring_g(n:n) = CHAR(0)
    ENDDO

    mixingrulestring_s = 'maxwellgarnett'
    hoststring_s = 'air'
    matrixstring_s = 'water'
    inclusionstring_s = 'spheroidal'
    hostmatrixstring_s = 'icewater'
    hostinclusionstring_s = 'spheroidal'

    mixingrulestring_g = 'maxwellgarnett'
    hoststring_g = 'air'
    matrixstring_g = 'water'
    inclusionstring_g = 'spheroidal'
    hostmatrixstring_g = 'icewater'
    hostinclusionstring_g = 'spheroidal'

    !..Create bins of snow (from 100 microns up to 2 cm).
    xxDx(1) = 100.D-6
    xxDx(nrbins+1) = 0.02d0
    DO n = 2, nrbins
       xxDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nrbins) &
            *DLOG(xxDx(nrbins+1)/xxDx(1)) +DLOG(xxDx(1)))
    ENDDO
    DO n = 1, nrbins
       xxDs(n) = DSQRT(xxDx(n)*xxDx(n+1))
       xdts(n) = xxDx(n+1) - xxDx(n)
    ENDDO

    !..Create bins of graupel (from 100 microns up to 5 cm).
    xxDx(1) = 100.D-6
    xxDx(nrbins+1) = 0.05d0
    DO n = 2, nrbins
       xxDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nrbins) &
            *DLOG(xxDx(nrbins+1)/xxDx(1)) +DLOG(xxDx(1)))
    ENDDO
    DO n = 1, nrbins
       xxDg(n) = DSQRT(xxDx(n)*xxDx(n+1))
       xdtg(n) = xxDx(n+1) - xxDx(n)
    ENDDO


    !..The calling program must set the m(D) relations and gamma shape
    !.. parameter mu for rain, snow, and graupel.  Easily add other types
    !.. based on the template here.  For majority of schemes with simpler
    !.. exponential number distribution, mu=0.

    xcre(1) = 1.0_r8 + xbm_r
    xcre(2) = 1.0_r8 + xmu_r
    xcre(3) = 4.0_r8 + xmu_r
    xcre(4) = 7.0_r8 + xmu_r
    DO n = 1, 4
       xcrg(n) = WGAMMA(xcre(n))
    ENDDO
    xorg2 = 1.0_r8/xcrg(2)

    xcse(1) = 1.0_r8 + xbm_s
    xcse(2) = 1.0_r8 + xmu_s
    xcse(3) = 4.0_r8 + xmu_s
    xcse(4) = 7.0_r8 + xmu_s
    DO n = 1, 4
       xcsg(n) = WGAMMA(xcse(n))
    ENDDO

    xosg2 = 1.0_r8/xcsg(2)

    xcge(1) = 1.0_r8 + xbm_g
    xcge(2) = 1.0_r8 + xmu_g
    xcge(3) = 4.0_r8 + xmu_g
    xcge(4) = 7.0_r8 + xmu_g
    DO n = 1, 4
       xcgg(n) = WGAMMA(xcge(n))
    ENDDO
    xogg2 = 1.0_r8/xcgg(2)

    xobmr = 1.0_r8/xbm_r
    xoams = 1.0_r8/xam_s
    xobms = 1.0_r8/xbm_s
    xocms = xoams**xobms
    xoamg = 1.0_r8/xam_g
    xobmg = 1.0_r8/xbm_g
    xocmg = xoamg**xobmg


  END SUBROUTINE radar_init



  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8) FUNCTION GAMMLN(XX)
    !     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: XX
    REAL(KIND=r8), PARAMETER:: STP = 2.5066282746310005D0
    REAL(KIND=r8), DIMENSION(6), PARAMETER:: &
         COF = (/76.18009172947146D0, -86.50532032941677D0, &
         24.01409824083091D0, -1.231739572450155D0, &
         .1208650973866179D-2, -.5395239384953D-5/)
    REAL(KIND=r8):: SER,TMP,X,Y
    INTEGER:: J

    X=XX
    Y=X
    TMP=X+5.5D0
    TMP=(X+0.5D0)*LOG(TMP)-TMP
    SER=1.000000000190015D0
    DO J=1,6
       Y=Y+1.D0
       SER=SER+COF(J)/Y
    END DO
    GAMMLN=TMP+LOG(STP*SER/X)
  END FUNCTION GAMMLN


  !+---+-----------------------------------------------------------------+

  SUBROUTINE rayleigh_soak_wetgraupel (x_g, a_geo, b_geo, fmelt,    &
       meltratio_outside, m_w, m_i, lambda, C_back,       &
       mixingrule,matrix,inclusion,                       &
       host,hostmatrix,hostinclusion)

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in):: x_g, a_geo, b_geo, fmelt, lambda,  &
         meltratio_outside
    REAL(KIND=r8), INTENT(out):: C_back
    COMPLEX*16, INTENT(in):: m_w, m_i
    CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion,     &
         host, hostmatrix, hostinclusion

    COMPLEX*16:: m_core, m_air
    REAL(KIND=r8):: D_large, D_g, rhog, x_w, fm, fmgrenz,    &
         volg, vg, volair, volice, volwater,            &
         meltratio_outside_grenz, mra,aa
    INTEGER:: error
    REAL(KIND=r8), PARAMETER:: PIx=3.1415926535897932384626434d0
    aa=lambda
    !     refractive index of air:
    m_air = (1.0d0,0.0d0)

    !     Limiting the degree of melting --- for safety: 
    fm = DMAX1(DMIN1(fmelt, 1.0d0), 0.0d0)
    !     Limiting the ratio of (melting on outside)/(melting on inside):
    mra = DMAX1(DMIN1(meltratio_outside, 1.0d0), 0.0d0)

    !    ! The relative portion of meltwater melting at outside should increase
    !    ! from the given input value (between 0 and 1)
    !    ! to 1 as the degree of melting approaches 1,
    !    ! so that the melting particle "converges" to a water drop.
    !    ! Simplest assumption is linear:
    mra = mra + (1.0d0-mra)*fm

    x_w = x_g * fm

    D_g = a_geo * x_g**b_geo

    IF (D_g .GE. 1d-12) THEN

       vg = PIx/6.0_r8 * D_g**3
       rhog = DMAX1(DMIN1(x_g / vg, 900.0d0), 10.0d0)
       vg = x_g / rhog

       meltratio_outside_grenz = 1.0d0 - rhog / 1000.0_r8

       IF (mra .LE. meltratio_outside_grenz) THEN
          !..In this case, it cannot happen that, during melting, all the
          !.. air inclusions within the ice particle get filled with
          !.. meltwater. This only happens at the end of all melting.
          volg = vg * (1.0d0 - mra * fm)

       ELSE
          !..In this case, at some melting degree fm, all the air
          !.. inclusions get filled with meltwater.
          fmgrenz=(900.0_r8-rhog)/(mra*900.0_r8-rhog+900.0_r8*rhog/1000.0_r8)

          IF (fm .LE. fmgrenz) THEN
             !.. not all air pockets are filled:
             volg = (1.0_r8 - mra * fm) * vg
          ELSE
             !..all air pockets are filled with meltwater, now the
             !.. entire ice sceleton melts homogeneously:
             volg = (x_g - x_w) / 900.0_r8 + x_w / 1000.0_r8
          ENDIF

       ENDIF

       D_large  = (6.0_r8 / PIx * volg) ** (1.0_r8/3.0_r8)
       volice = (x_g - x_w) / (volg * 900.0_r8)
       volwater = x_w / (1000.0_r8 * volg)
       volair = 1.0_r8 - volice - volwater

       !..complex index of refraction for the ice-air-water mixture
       !.. of the particle:
       m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,      &
            volwater, mixingrule, host, matrix, inclusion, &
            hostmatrix, hostinclusion, error)
       IF (error .NE. 0) THEN
          C_back = 0.0d0
          RETURN
       ENDIF

       !..Rayleigh-backscattering coefficient of melting particle: 
       C_back = (ABS((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2           &
            * PI5 * D_large**6 / lamda4

    ELSE
       C_back = 0.0d0
    ENDIF

  END SUBROUTINE rayleigh_soak_wetgraupel

  !+---+-----------------------------------------------------------------+

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION m_complex_water_ray(lambda,T)

    !      Complex refractive Index of Water as function of Temperature T
    !      [deg C] and radar wavelength lambda [m]; valid for
    !      lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
    !      after Ray (1972)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: T,lambda
    REAL(KIND=r8):: epsinf,epss,epsr,epsi
    REAL(KIND=r8):: alpha,lambdas,nenner
    !COMPLEX*16, PARAMETER:: i = (0e0_r8,1e0_r8)
    REAL(KIND=r8), PARAMETER:: PIx=3.1415926535897932384626434d0

    epsinf  = 5.27137d0 + 0.02164740d0 * T - 0.00131198d0 * T*T
    epss    = 78.54d+0 * (1.0_r8 - 4.579d-3 * (T - 25.0_r8)                 &
         + 1.190d-5 * (T - 25.0_r8)*(T - 25.0_r8)                        &
         - 2.800d-8 * (T - 25.0_r8)*(T - 25.0_r8)*(T - 25.0_r8))
    alpha   = -16.8129d0/(T+273.16_r8) + 0.0609265d0
    lambdas = 0.00033836d0 * EXP(2513.98d0/(T+273.16_r8)) * 1e-2_r8

    nenner = 1.d0+2.d0*(lambdas/lambda)**(1d0-alpha)*SIN(alpha*PIx*0.5_r8) &
         + (lambdas/lambda)**(2d0-2d0*alpha)
    epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)   &
         * SIN(alpha*PIx*0.5_r8)+1d0)) / nenner
    epsi = ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)            &
         * COS(alpha*PIx*0.5_r8)+0d0)) / nenner                           &
         + lambda*1.25664_r8/1.88496_r8

    m_complex_water_ray = SQRT(CMPLX(epsr,-epsi))

  END FUNCTION m_complex_water_ray

  !+---+-----------------------------------------------------------------+


  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION m_complex_ice_maetzler(lambda,T)

    !      complex refractive index of ice as function of Temperature T
    !      [deg C] and radar wavelength lambda [m]; valid for
    !      lambda in [0.0001,30] m; T in [-250.0,0.0] C
    !      Original comment from the Matlab-routine of Prof. Maetzler:
    !      Function for calculating the relative permittivity of pure ice in
    !      the microwave region, according to C. Maetzler, "Microwave
    !      properties of ice and snow", in B. Schmitt et al. (eds.) Solar
    !      System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
    !      Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
    !      TK = temperature (K), range 20 to 273.15
    !      f = frequency in GHz, range 0.01 to 3000

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: T,lambda
    REAL(KIND=r8):: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa

    c = 2.99d8
    TK = T + 273.16_r8
    f = c / lambda * 1d-9

    B1 = 0.0207_r8
    B2 = 1.16d-11
    b = 335.0d0
    deltabeta = EXP(-10.02_r8 + 0.0364_r8*(TK-273.16_r8))
    betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1)**2) ) + B2*f*f
    beta = betam + deltabeta
    theta = 300.0_r8 / TK - 1.0_r8
    alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
    m_complex_ice_maetzler = 3.1884_r8 + 9.1e-4_r8*(TK-273.16_r8)
    m_complex_ice_maetzler = m_complex_ice_maetzler                   &
         + CMPLX(0.0d0, (alfa/f + beta*f)) 
    m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))

  END FUNCTION m_complex_ice_maetzler

  !+---+-----------------------------------------------------------------+

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION get_m_mix_nested (m_a, m_i, m_w, volair,      &
       volice, volwater, mixingrule, host, matrix,        &
       inclusion, hostmatrix, hostinclusion, cumulerror)

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in):: volice, volair, volwater
    COMPLEX*16, INTENT(in):: m_a, m_i, m_w
    CHARACTER(len=*), INTENT(in):: mixingrule, host, matrix,          &
         inclusion, hostmatrix, hostinclusion
    INTEGER, INTENT(out):: cumulerror

    REAL(KIND=r8):: vol1, vol2
    COMPLEX*16:: mtmp
    INTEGER:: error

    !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
    !.. air, ice, or water

    cumulerror = 0
    get_m_mix_nested = CMPLX(1.0d0,0.0d0)

    IF (host .EQ. 'air') THEN

       IF (matrix .EQ. 'air') THEN
          WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
          CALL wrf_debug(150, radar_debug)
          cumulerror = cumulerror + 1
       ELSE
          vol1 = volice / MAX(volice+volwater,1d-10)
          vol2 = 1.0d0 - vol1
          mtmp = get_m_mix (m_a, m_i, m_w, 0.0d0, vol1, vol2,             &
               mixingrule, matrix, inclusion, error)
          cumulerror = cumulerror + error

          IF (hostmatrix .EQ. 'air') THEN
             get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                  volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                  hostmatrix, hostinclusion, error)
             cumulerror = cumulerror + error
          ELSEIF (hostmatrix .EQ. 'icewater') THEN
             get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                  volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                  'ice', hostinclusion, error)
             cumulerror = cumulerror + error
          ELSE
             WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                  hostmatrix
             CALL wrf_debug(150, radar_debug)
             cumulerror = cumulerror + 1
          ENDIF
       ENDIF

    ELSEIF (host .EQ. 'ice') THEN

       IF (matrix .EQ. 'ice') THEN
          WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
          CALL wrf_debug(150, radar_debug)
          cumulerror = cumulerror + 1
       ELSE
          vol1 = volair / MAX(volair+volwater,1d-10)
          vol2 = 1.0d0 - vol1
          mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0d0, vol2,             &
               mixingrule, matrix, inclusion, error)
          cumulerror = cumulerror + error

          IF (hostmatrix .EQ. 'ice') THEN
             get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                  (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                  hostmatrix, hostinclusion, error)
             cumulerror = cumulerror + error
          ELSEIF (hostmatrix .EQ. 'airwater') THEN
             get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                  (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                  'air', hostinclusion, error)
             cumulerror = cumulerror + error          
          ELSE
             WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                  hostmatrix
             CALL wrf_debug(150, radar_debug)
             cumulerror = cumulerror + 1
          ENDIF
       ENDIF

    ELSEIF (host .EQ. 'water') THEN

       IF (matrix .EQ. 'water') THEN
          WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
          CALL wrf_debug(150, radar_debug)
          cumulerror = cumulerror + 1
       ELSE
          vol1 = volair / MAX(volice+volair,1d-10)
          vol2 = 1.0d0 - vol1
          mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0d0,             &
               mixingrule, matrix, inclusion, error)
          cumulerror = cumulerror + error

          IF (hostmatrix .EQ. 'water') THEN
             get_m_mix_nested = get_m_mix (2*m_a, mtmp, m_w,                &
                  0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                  hostmatrix, hostinclusion, error)
             cumulerror = cumulerror + error
          ELSEIF (hostmatrix .EQ. 'airice') THEN
             get_m_mix_nested = get_m_mix (2*m_a, mtmp, m_w,                &
                  0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                  'ice', hostinclusion, error)
             cumulerror = cumulerror + error          
          ELSE
             WRITE(radar_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
                  hostmatrix
             CALL wrf_debug(150, radar_debug)
             cumulerror = cumulerror + 1
          ENDIF
       ENDIF

    ELSEIF (host .EQ. 'none') THEN

       get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                     &
            volair, volice, volwater, mixingrule,            &
            matrix, inclusion, error)
       cumulerror = cumulerror + error

    ELSE
       WRITE(radar_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
       CALL wrf_debug(150, radar_debug)
       cumulerror = cumulerror + 1
    ENDIF

    IF (cumulerror .NE. 0) THEN
       WRITE(radar_debug,*) 'GET_M_MIX_NESTED: error encountered'
       CALL wrf_debug(150, radar_debug)
       get_m_mix_nested = CMPLX(1.0d0,0.0d0)    
    ENDIF

  END FUNCTION get_m_mix_nested

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION get_m_mix (m_a, m_i, m_w, volair, volice,     &
       volwater, mixingrule, matrix, inclusion, error)

    IMPLICIT NONE

    REAL(KIND=r8), INTENT(in):: volice, volair, volwater
    COMPLEX*16, INTENT(in):: m_a, m_i, m_w
    CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion
    INTEGER, INTENT(out):: error

    error = 0
    get_m_mix = CMPLX(1.0d0,0.0d0)

    IF (mixingrule .EQ. 'maxwellgarnett') THEN
       IF (matrix .EQ. 'ice') THEN
          get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,  &
               m_i, m_a, m_w, inclusion, error)
       ELSEIF (matrix .EQ. 'water') THEN
          get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,  &
               m_w, m_a, m_i, inclusion, error)
       ELSEIF (matrix .EQ. 'air') THEN
          get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,  &
               m_a, m_w, m_i, inclusion, error)
       ELSE
          WRITE(radar_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
          CALL wrf_debug(150, radar_debug)
          error = 1
       ENDIF

    ELSE
       WRITE(radar_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
       CALL wrf_debug(150, radar_debug)
       error = 2
    ENDIF

    IF (error .NE. 0) THEN
       WRITE(radar_debug,*) 'GET_M_MIX: error encountered'
       CALL wrf_debug(150, radar_debug)
    ENDIF

  END FUNCTION get_m_mix

  !+---+-----------------------------------------------------------------+

  COMPLEX*16 FUNCTION m_complex_maxwellgarnett(vol1, vol2, vol3,    &
       m1, m2, m3, inclusion, error)

    IMPLICIT NONE

    COMPLEX*16 :: m1, m2, m3
    REAL(KIND=r8) :: vol1, vol2, vol3
    CHARACTER(len=*) :: inclusion

    COMPLEX*16 :: beta2, beta3, m1t, m2t, m3t
    INTEGER, INTENT(out) :: error

    error = 0

    IF (DABS(vol1+vol2+vol3-1.0d0) .GT. 1d-6) THEN
       WRITE(radar_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',       &
            'partial volume fractions is not 1...ERROR'
       CALL wrf_debug(150, radar_debug)
       m_complex_maxwellgarnett=CMPLX(-999.99d0,-999.99d0)
       error = 1
       RETURN
    ENDIF

    m1t = m1**2
    m2t = m2**2
    m3t = m3**2

    IF (inclusion .EQ. 'spherical') THEN
       beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
       beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
    ELSEIF (inclusion .EQ. 'spheroidal') THEN
       beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
       beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
    ELSE
       WRITE(radar_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ',                  &
            'unknown inclusion: ', inclusion
       CALL wrf_debug(150, radar_debug)
       m_complex_maxwellgarnett=DCMPLX(-999.99d0,-999.99d0)
       error = 1
       RETURN
    ENDIF

    m_complex_maxwellgarnett = &
         SQRT(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
         (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

  END FUNCTION m_complex_maxwellgarnett
  !  (C) Copr. 1986-92 Numerical Recipes Software 2.02
  !+---+-----------------------------------------------------------------+
  REAL(KIND=r8) FUNCTION WGAMMA(y)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN):: y

    WGAMMA = EXP(GAMMLN(y))

  END FUNCTION WGAMMA

  !+---+-----------------------------------------------------------------+


  !+---+-----------------------------------------------------------------+

  SUBROUTINE wrf_debug( level , str ) 
    IMPLICIT NONE 
    CHARACTER(Len=*), INTENT (IN) :: str 
    INTEGER         , INTENT (IN) :: level
    WRITE(0,'(A)')str ,level
    RETURN 
  END SUBROUTINE wrf_debug

  !+---+-----------------------------------------------------------------+


  !===============================================================================
  SUBROUTINE geopotential_t(                                 &
       piln   ,  pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   ,nCols, kMax)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute the geopotential height (above the surface) at the midpoints and 
    ! interfaces using the input temperatures and pressures.
    !
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! Number of longitudes
    INTEGER, INTENT(in) :: nCols
    INTEGER, INTENT(in) :: kMax
    REAL(r8), INTENT(in) :: piln (nCols,kMax+1)   ! Log interface pressures
    REAL(r8), INTENT(in) :: pint (nCols,kMax+1)   ! Interface pressures
    REAL(r8), INTENT(in) :: pmid (nCols,kMax)    ! Midpoint pressures
    REAL(r8), INTENT(in) :: pdel (nCols,kMax)    ! layer thickness
    REAL(r8), INTENT(in) :: rpdel(nCols,kMax)    ! inverse of layer thickness
    REAL(r8), INTENT(in) :: t    (nCols,kMax)    ! temperature
    REAL(r8), INTENT(in) :: q    (nCols,kMax)    ! specific humidity
    REAL(r8), INTENT(in) :: rair                 ! Gas constant for dry air
    REAL(r8), INTENT(in) :: gravit               ! Acceleration of gravity
    REAL(r8), INTENT(in) :: zvir                 ! rh2o/rair - 1

    ! Output arguments

    REAL(r8), INTENT(out) :: zi(nCols,kMax+1)     ! Height above surface at interfaces
    REAL(r8), INTENT(out) :: zm(nCols,kMax)      ! Geopotential height at mid level
    !
    !---------------------------Local variables-----------------------------
    !
    LOGICAL  :: fvdyn              ! finite volume dynamics
    INTEGER  :: i,k                ! Lon, level indices
    REAL(r8) :: hkk(nCols)         ! diagonal element of hydrostatic matrix
    REAL(r8) :: hkl(nCols)         ! off-diagonal element
    REAL(r8) :: rog                ! Rair / gravit
    REAL(r8) :: tv                 ! virtual temperature
    REAL(r8) :: tvfac              ! Tv/T
    zi= 0.0_r8;    zm= 0.0_r8;    hkk= 0.0_r8;hkl= 0.0_r8
    rog= 0.0_r8;tv = 0.0_r8;tvfac= 0.0_r8
    !
    !-----------------------------------------------------------------------
    !
    rog = rair/gravit

    ! Set dynamics flag

    fvdyn = .FALSE.!dycore_is ('LR')

    ! The surface height is zero by definition.

    DO i = 1,ncol
       zi(i,kMax+1) = 0.0_r8
    END DO

    ! Compute zi, zm from bottom up. 
    ! Note, zi(i,k) is the interface above zm(i,k)

    DO k = kMax, 1, -1

       ! First set hydrostatic elements consistent with dynamics

       IF (fvdyn) THEN
          DO i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1.0_r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          END DO
       ELSE
          DO i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF

       ! Now compute tv, zm, zi

       DO i = 1,ncol
          tvfac   = 1.0_r8 + zvir * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv * hkl(i)
       END DO
    END DO

    RETURN
  END SUBROUTINE geopotential_t


END MODULE Micro_HugMorr
