Module CalCloudCover
IMPLICIT NONE
SAVE
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  !INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(6) ! Kind for 64-bits Real Numbers

CONTAINS








! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Ozone set in lower
!                                               levels.
!                                               (J. M. Edwards)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!       4.5             18-05-98                Provision for treating
!                                               extra (H)(C)FCs
!                                               included.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Move HCFCs to a more
!                                               natural place in the
!                                               code.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Remove the explicit
!                                               limit on the
!                                               concentration of
!                                               water vapour.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             29-05-02                Add diagnostic call
!                                               for column-integrated
!                                               cloud droplet number.
!                                               (A. Jones)
!       5.4             25-04-02                Replace land/sea mask
!                                               with land fraction in
!                                               call to NUMBER_DROPLET
!                                               (A. Jones)
!       5.5             17-02-03                Change I_CLIM_POINTER
!                                               for hp compilation
!                                               (M.Hughes)
!  6.1   20/08/03  Code for STOCHEM feedback.  C. Johnson
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!  6.2   03-11-05   Enable HadGEM1 climatological aerosols. C. F. Durman
!                   Reworked to use switch instead of #defined. R Barnes
!  6.2   25/05/05  Convert compilation into a more universally portable
!                  form. Tom Edwards
!  6.2   15/12/05  Set negative specific humidities to zero.
!                                               (J. Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Old formulation over
!                                               sea-ice removed.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Ground temperature
!                                               set equal to that
!                                               in the middle of the
!                                               bottom layer.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.3             25-04-01                Alter the specification
!                                               of temperature on rho
!                                               levels (layer
!                                               boundaries).  S.Cusack
!       5.3             25-04-01   Gather land, sea and
!                                  sea-ice temperatures and
!                                  land fraction. Replace TFS
!                                  with general sea temperature.
!                                       (N. Gedney)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New flag L_AEROSOL_CCN
!                                               introduced to allow
!                                               inclusion of indirect
!                                               aerosol forcing alone.
!                                               Correction of comments
!                                               for LCCWC1 and LCCWC2.
!                                               Correction of level at
!                                               which temperature for
!                                               partitioning
!                                               convective homogeneously
!                                               mixed cloud is taken.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             15-09-97                A parametrization of
!                                               ice crystals with a
!                                               temperature dependedence
!                                               of the size has been
!                                               added.
!                                               Explicit checking of
!                                               the sizes of particles
!                                               for the domain of
!                                               validity of the para-
!                                               metrization has been
!                                               added.
!                                               (J. M. Edwards)
!       5.0             15-04-98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  A.C.Bushell
!       4.5             18-05-98                New option for
!                                               partitioning between
!                                               ice and water in
!                                               convective cloud
!                                               included.
!                                               (J. M. Edwards)
!       4.5             13/05/98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  S. Cusack
!       5.1             04-04-00                Remove obsolete tests
!                                               for convective cloud
!                                               and removal of very
!                                               thin cloud (no longer
!                                               required with current
!                                               solvers, but affects
!                                               bit-comparison).
!                                               (J. M. Edwards)
!       5.1             06-04-00                Correct some comments
!                                               and error messages.
!                                               (J. M. Edwards)
!       5.2             10-11-00                With local partitioning
!                                               of convective cloud
!                                               between water and ice
!                                               force homogeneous
!                                               nucleation at
!                                               -40 Celsius.
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Pass sea-salt variables
!                                               down to R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       5.3             11-10-01                Convert returned
!                                               diagnostics to 2-D
!                                               arrays.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check on small cloud
!                                               fractions and liquid for
!                                               contents for PC2 scheme.
!                                               (D. Wilson)
!       5.5             24-02-03                Addition of new ice
!                                               aggregate
!                                               parametrization.
!                                               (J. M. Edwards)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       6.1             07-04-04                Add variables for
!                                               column-droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Pass through PC2 logical
!                                               (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_CLOUD_FIELD( &
       N_PROFILE               , &
       NLEVS                   , &
       N_LAYER                 , &
       NCLDS                   , &
       I_GATHER                , &
       P                       , &
       T                       , &
       D_MASS                  , &
       CCB                     , &
       CCT                     , &
       CCA                     , &
       CCCWP                   , &
       LCCWC1                  , &
       LCCWC2                  , &
       LCA_AREA                , &
       LCA_BULK                , &
       L_PC2                   , &
       L_MICROPHYSICS          , &
       L_AEROSOL_CCN           , &
       SEA_SALT_FILM           , &
       SEA_SALT_JET            , &
       L_SEASALT_CCN           , &
       SALT_DIM_A              , &
       SALT_DIM_B              , &
       L_USE_BIOGENIC          , &
       BIOGENIC                , &
       BIOGENIC_DIM1           , &
       BIOGENIC_DIM2           , &
       SULP_DIM1               , &
       SULP_DIM2               , &
       ACCUM_SULPHATE          , &
       DISS_SULPHATE           , &
       AITKEN_SULPHATE         , &
       L_BIOMASS_CCN           , &
       BMASS_DIM1              , &
       BMASS_DIM2              , &
       AGED_BMASS              , &
       CLOUD_BMASS             , &
       L_OCFF_CCN              , &
       OCFF_DIM1               , &
       OCFF_DIM2               , &
       AGED_OCFF               , &
       CLOUD_OCFF              , &
       LYING_SNOW              , &
       L_CLOUD_WATER_PARTITION , &
       LAND_G                  , &
       FLANDG_G                , &
       I_CLOUD_REPRESENTATION  , &
       I_CONDENSED_PARAM       , &
       CONDENSED_MIN_DIM       , &
       CONDENSED_MAX_DIM       , &
       N_CONDENSED             , &
       TYPE_CONDENSED          , &
       W_CLOUD                 , &
       FRAC_CLOUD              , &
       L_LOCAL_CNV_PARTITION   , &
       CONDENSED_MIX_RAT_AREA  , &
       CONDENSED_DIM_CHAR      , &
       RE_CONV                 , &
       RE_CONV_FLAG            , &
       RE_STRAT                , &
       RE_STRAT_FLAG           , &
       WGT_CONV                , &
       WGT_CONV_FLAG           , &
       WGT_STRAT               , &
       WGT_STRAT_FLAG          , &
       LWP_STRAT               , &
       LWP_STRAT_FLAG          , &
       NTOT_DIAG               , &
       NTOT_DIAG_FLAG          , &
       STRAT_LWC_DIAG          , &
       STRAT_LWC_DIAG_FLAG     , &
       SO4_CCN_DIAG            , &
       SO4_CCN_DIAG_FLAG       , &
       COND_SAMP_WGT           , &
       COND_SAMP_WGT_FLAG      , &
       NC_DIAG                 , &
       NC_DIAG_FLAG            , &
       NC_WEIGHT               , &
       NC_WEIGHT_FLAG          , &
       col_list                , &
       row_list                , &
       row_length              , &
       rows                    , &
       NPD_FIELD               , &
       NPD_PROFILE             , &
       NPD_LAYER               , &
       NPD_AEROSOL_SPECIES     , &
       N_CCA_LEV               , &
       L_3D_CCA                , &
       Ntot_land               , &
       Ntot_sea                  )
!
!
!
      IMPLICIT NONE

!
!     COMDECKS INCLUDED.
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
! CLDCMP3A sets components of clouds for two-stream radiation code.

      ! stratiform water droplets
      INTEGER,PARAMETER:: IP_CLCMP_ST_WATER=1

      ! stratiform ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_ST_ICE=2

      ! convective water droplets
      INTEGER,PARAMETER:: IP_CLCMP_CNV_WATER=3

      ! convective ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_CNV_ICE=4

! CLDCMP3A end
! CLDTYP3A defines cloud types for TWO-STREAM RADIATION CODE.

      INTEGER,PARAMETER:: IP_CLOUD_TYPE_HOMOGEN=1 ! water and ice
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_WATER=1   ! Water only
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_ICE=2     ! Ice only

      ! mixed-phase stratiform cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_STRAT=1

      ! mixed-phase convective cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CONV=2

      INTEGER,PARAMETER:: IP_CLOUD_TYPE_SW=1 ! stratiform water cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_SI=2 ! stratiform ice cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CW=3 ! convective water cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CI=4 ! convective ice cloud

! CLDTYP3A end
! CLREPP3A defines representations of clouds in two-stream radiation
! code.

      ! all components are mixed homogeneously
      INTEGER,PARAMETER:: IP_CLOUD_HOMOGEN     = 1

      ! ice and water clouds are treated separately
      INTEGER,PARAMETER:: IP_CLOUD_ICE_WATER   = 2

      ! clouds are divided into homogeneously mixed stratiform and
      ! convective parts
      INTEGER,PARAMETER:: IP_CLOUD_CONV_STRAT  = 3

      ! clouds divided into ice and water phases and into stratiform and
      ! convective components.
      INTEGER,PARAMETER:: IP_CLOUD_CSIW        = 4

! Types of clouds (values in CLREPD3A)

      ! number of type of clouds in representation
      INTEGER :: NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)

      ! map of components contributing to types of clouds
      INTEGER :: IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT,NPD_CLOUD_REPRESENTATION)

! CLREPP3A end
! ICLPRM3A defines numbers for ice cloud schemes in two-stream radiation
! code.
!
!   5.5   Feb 2003     Addition of the aggregate parametrization.
!                                                 John Edwards
!   6.2   Jan 2006     Various options for radiation code
!                      3Z added.   (J.-C. Thelen)
!

      ! number of cloud fitting schemes
      INTEGER,PARAMETER:: NPD_ICE_CLOUD_FIT=10

      ! parametrization of slingo and schrecker.
      INTEGER,PARAMETER:: IP_SLINGO_SCHRECKER_ICE=1

      ! unparametrized ice crystal data
       INTEGER,PARAMETER:: IP_ICE_UNPARAMETRIZED=3

      ! sun and shine's parametrization in the visible (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_VIS=4

      ! sun and shine's parametrization in the ir (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_IR=5

      ! scheme based on anomalous diffraction theory for ice crystals
      INTEGER,PARAMETER:: IP_ICE_ADT=6

!           ADT-based scheme for ice crystals using 10th order
!           polynomials
      INTEGER,PARAMETER:: IP_ICE_ADT_10=7
      ! Provisional agregate parametrization.
      INTEGER,PARAMETER:: IP_ICE_AGG_DE=12
!           Fu's parametrization in the solar region of the spectrum
      INTEGER,PARAMETER:: IP_ICE_FU_SOLAR=9
!           Fu's parametrization in the infra-red region of the spectrum
      INTEGER,PARAMETER:: IP_ICE_FU_IR=10
!           Parametrization of Slingo and Schrecker
!           (Moments of phase function).
      INTEGER,PARAMETER:: IP_SLINGO_SCHR_ICE_PHF=11
!           Parametrization like Fu
!           (Moments of phase function).
      INTEGER,PARAMETER:: IP_ICE_FU_PHF=13

! ICLPRM3A end
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

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

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
      Real, Parameter  :: sclht  = 6.8E+03  ! Scale Height H

!
!*----------------------------------------------------------------------
!
!
!     DIMENSIONS OF ARRAYS:
      Integer, Intent(IN) :: row_length !Number of grid-points in EW-direction
!                              in the local domain
      Integer, Intent(IN) :: rows!Number of grid-points in NS-direction
!                              in the local domain
      INTEGER, INTENT(IN) :: NPD_FIELD
!             FIELD SIZE IN CALLING PROGRAM
      INTEGER, INTENT(IN) :: NPD_PROFILE
!             SIZE OF ARRAY OF PROFILES
      INTEGER, INTENT(IN) :: NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
      INTEGER, INTENT(IN) :: NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOL_SPECIES
      INTEGER, INTENT(IN) :: SULP_DIM1
!             1ST DIMENSION OF ARRAYS OF SULPHATE
      INTEGER, INTENT(IN) :: SULP_DIM2
!             2ND DIMENSION OF ARRAYS OF SULPHATE
      INTEGER, INTENT(IN) :: BMASS_DIM1
!             1ST DIMENSION OF ARRAYS OF BIOMASS SMOKE
      INTEGER, INTENT(IN) :: BMASS_DIM2
!             2ND DIMENSION OF ARRAYS OF BIOMASS SMOKE
      INTEGER, INTENT(IN) :: OCFF_DIM1
!             1ST DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
      INTEGER, INTENT(IN) :: OCFF_DIM2
!             2ND DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
      INTEGER, INTENT(IN) :: SALT_DIM_A
!             1ST DIMENSION OF ARRAYS OF SEA-SALT
      INTEGER, INTENT(IN) :: SALT_DIM_B
!             2ND DIMENSION OF ARRAYS OF SEA-SALT
      INTEGER, INTENT(IN) :: BIOGENIC_DIM1
!             1ST DIMENSION OF BIOGENIC AEROSOL ARRAY
      INTEGER, INTENT(IN) :: BIOGENIC_DIM2
!             2ND DIMENSION OF BIOGENIC AEROSOL ARRAY
      INTEGER, INTENT(IN) :: N_CCA_LEV
!             NUMBER OF LEVELS FOR CONVECTIVE CLOUD AMOUNT
!
!     ACTUAL SIZES USED:
      INTEGER, INTENT(IN) :: N_PROFILE
!             NUMBER OF PROFILES
      INTEGER, INTENT(IN) :: NLEVS
!             Number of layers used outside the radiation scheme
      INTEGER, INTENT(IN) :: N_LAYER
!             Number of layers seen by the radiation code
      INTEGER, INTENT(IN) :: NCLDS
!             NUMBER OF CLOUDY LEVELS
!
!     GATHERING ARRAY:
      INTEGER, INTENT(IN) :: I_GATHER(NPD_FIELD) !LIST OF POINTS TO BE GATHERED
      Integer, Intent(IN) :: col_list(npd_field) !EW indices of gathered points in the 2-D domain
      Integer, Intent(IN) :: row_list(npd_field) !NS indices of gathered points in the 2-D domain
!
!     THERMODYNAMIC FIELDS:
      REAL   , INTENT(IN) :: P     (NPD_PROFILE, 0: NPD_LAYER)!PRESSURES
      REAL   , INTENT(IN) :: T     (NPD_PROFILE, 0: NPD_LAYER)!TEMPERATURES
      REAL   , INTENT(IN) :: D_MASS(NPD_PROFILE,    NPD_LAYER)!MASS THICKNESSES OF LAYERS
!
!     CONVECTIVE CLOUDS:
      INTEGER, INTENT(IN) :: CCB(NPD_FIELD)          !BASE OF CONVECTIVE CLOUD
      INTEGER, INTENT(IN) :: CCT(NPD_FIELD)          !TOP OF CONVECTIVE CLOUD
      REAL   , INTENT(IN) :: CCA(NPD_FIELD,N_CCA_LEV)!FRACTION OF CONVECTIVE CLOUD
      REAL   , INTENT(IN) :: CCCWP(NPD_FIELD)        !WATER PATH OF CONVECTIVE CLOUD
      LOGICAL, INTENT(IN) :: L_3D_CCA 
      LOGICAL, INTENT(IN) :: L_LOCAL_CNV_PARTITION

!             FLAG TO CARRY OUT THE PARTITIONING BETWEEN ICE
!             AND WATER IN CONVECTIVE CLOUDS AS A FUNCTION OF
!             THE LOCAL TEMPERATURE

      LOGICAL, INTENT(IN) :: L_SEASALT_CCN
!              FLAG FOR SEA-SALT PARAMETRIZATION FOR CCN
      LOGICAL, INTENT(IN) :: L_USE_BIOGENIC
!              FLAG TO USE BIOGENIC AEROSOLS AS CCN
      LOGICAL, INTENT(IN) :: L_BIOMASS_CCN
!              FLAG FOR SEA-SALT PARAMETRIZATION FOR CCN
      LOGICAL, INTENT(IN) :: L_OCFF_CCN
!              FLAG FOR FOSSIL-FUEL ORG CARB PARAMETRIZATION FOR CCN
!
!     LAYER CLOUDS:
      REAL, INTENT(IN) ::  LCCWC1(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             LIQUID WATER CONTENTS
      REAL, INTENT(IN) ::  LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             ICE WATER CONTENTS
      REAL, INTENT(IN) ::  LCA_AREA(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             AREA COVERAGE FRACTIONS OF LAYER CLOUDS
      REAL, INTENT(IN) ::  LCA_BULK(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             BULK COVERAGE FRACTIONS OF LAYER CLOUDS
!
!     ARRAYS FOR MICROPHYSICS:
      LOGICAL , INTENT(IN) :: L_MICROPHYSICS 
!             MICROPHYSICAL FLAG
      LOGICAL , INTENT(IN) :: L_PC2
!             PC2 cloud scheme is in use
      LOGICAL , INTENT(IN) :: L_AEROSOL_CCN
!             FLAG TO USE AEROSOLS TO FIND CCN
      LOGICAL , INTENT(IN) :: L_CLOUD_WATER_PARTITION
!             FLAG TO USE PROGNOSTIC CLOUD ICE CONTENTS
      LOGICAL , INTENT(IN) :: LAND_G(NPD_PROFILE)
!             FLAG FOR LAND POINTS
      REAL   , INTENT(IN)  ::  ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF ACCUMULATION-MODE SULPHATE
      REAL   , INTENT(IN)  ::  AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)
!             Mixing ratios of Aitken-mode sulphate
      REAL   , INTENT(IN)  ::  DISS_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MIXING RATIOS OF DISSOLVED SULPHATE
      REAL   , INTENT(IN)  ::  AGED_BMASS(BMASS_DIM1, BMASS_DIM2)
!             MIXING RATIOS OF AGED BIOMASS SMOKE
      REAL   , INTENT(IN)  ::  CLOUD_BMASS(BMASS_DIM1, BMASS_DIM2)
!             MIXING RATIOS OF IN-CLOUD BIOMASS SMOKE
      REAL   , INTENT(IN)  ::  AGED_OCFF(OCFF_DIM1, OCFF_DIM2)
!             MIXING RATIOS OF AGED FOSSIL-FUEL ORGANIC CARBON
      REAL   , INTENT(IN)  ::  CLOUD_OCFF(OCFF_DIM1, OCFF_DIM2)
!             MIXING RATIOS OF IN-CLOUD FOSSIL-FUEL ORGANIC CARBON
      REAL   , INTENT(IN)  ::  SEA_SALT_FILM(SALT_DIM_A, SALT_DIM_B)
!             NUMBER CONCENTRATION OF FILM-MODE SEA-SALT AEROSOL
      REAL   , INTENT(IN)  ::  SEA_SALT_JET(SALT_DIM_A, SALT_DIM_B)
!             NUMBER CONCENTRATION OF JET-MODE SEA-SALT AEROSOL
      REAL   , INTENT(IN)  ::  BIOGENIC(BIOGENIC_DIM1, BIOGENIC_DIM2)
!             M.M.R. OF BIOGENIC AEROSOL
!
!     REPRESENTATION OF CLOUDS
      INTEGER, INTENT(IN) ::  I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS
!
!     PARAMETRIZATIONS FOR CLOUDS:
      INTEGER, INTENT(IN) ::  I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
!     LIMITS ON SIZES OF PARTICLES
      REAL   , INTENT(IN) :: CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
      REAL   , INTENT(IN) :: CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT
!
      Real   , Intent(IN) :: Ntot_land ! Number of droplets over land / m-3
      Real   , Intent(IN) :: Ntot_sea  ! Number of droplets over sea / m-3
!
!     ASSIGNED CLOUD FIELDS:
      INTEGER , INTENT(OUT) ::  N_CONDENSED!NUMBER OF CONDENSED COMPONENTS
      INTEGER , INTENT(OUT) ::  TYPE_CONDENSED(NPD_CLOUD_COMPONENT)
!             TYPES OF CONDENSED COMPONENTS
      REAL    , INTENT(OUT) :: W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             TOTAL AMOUNTS OF CLOUD
      REAL    , INTENT(OUT) :: FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTION OF EACH TYPE OF CLOUD
      REAL    , INTENT(OUT) :: CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT)
!             CHARACTERISTIC DIMENSIONS OF CLOUDY COMPONENTS
      REAL    , INTENT(OUT) :: CONDENSED_MIX_RAT_AREA(NPD_PROFILE, 0: NPD_LAYER , NPD_CLOUD_COMPONENT)
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING AREA CLD
      REAL     :: NTOT_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
      REAL     :: STRAT_LWC_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
      REAL     :: SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
!
!
      REAL, INTENT(IN   ) ::  LYING_SNOW(NPD_FIELD)     !  SNOW DEPTH (>5000m = LAND ICE SHEET)
      REAL ::  LYING_SNOW_G(NPD_PROFILE) !GATHERED VERSION OF THE ABOVE
      REAL, INTENT(IN   ) ::  FLANDG_G(NPD_PROFILE)     !GATHERED global LAND FRACTION FIELD
!
!     MICROPHYSICAL DIAGNOSTICS:
      LOGICAL, INTENT(IN   ) ::    RE_CONV_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
      LOGICAL, INTENT(IN   ) ::  RE_STRAT_FLAG
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
      LOGICAL, INTENT(IN   ) ::  WGT_CONV_FLAG
!             DIAGNOSE WEIGHT FOR CONVECTIVE CLOUD
      LOGICAL, INTENT(IN   ) ::  WGT_STRAT_FLAG
!             DIAGNOSE WEIGHT FOR STRATIFORM CLOUD
      LOGICAL, INTENT(IN   ) ::  LWP_STRAT_FLAG
!             DIAGNOSE LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
      LOGICAL, INTENT(IN   ) ::  NTOT_DIAG_FLAG
!             DIAGNOSE DROPLET CONCENTRATION*WEIGHT
      LOGICAL, INTENT(IN   ) ::  STRAT_LWC_DIAG_FLAG
!             DIAGNOSE STRATIFORM LWC*WEIGHT
      LOGICAL, INTENT(IN   ) ::  SO4_CCN_DIAG_FLAG
!             DIAGNOSE SO4 CCN MASS CONC*COND. SAMP. WEIGHT
      LOGICAL, INTENT(IN   ) ::  COND_SAMP_WGT_FLAG
!             DIAGNOSE CONDITIONAL SAMPLING WEIGHT
      LOGICAL, INTENT(IN   ) ::  NC_DIAG_FLAG
!             DIAGNOSE COLUMN DROPLET CONCENTRATION * SAMP. WEIGHT
      LOGICAL, INTENT(IN   ) ::  NC_WEIGHT_FLAG
!             DIAGNOSE COLUMN DROPLET SAMPLING WEIGHT
!
      REAL , INTENT(OUT  ) :: RE_CONV(row_length, rows, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
      REAL , INTENT(OUT   ) :: RE_STRAT(row_length, rows, NCLDS)
!             EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
      REAL , INTENT(OUT  ) :: WGT_CONV(row_length, rows, NCLDS)
!             WEIGHT FOR CONVECTIVE CLOUD
      REAL , INTENT(OUT   ) :: WGT_STRAT(row_length, rows, NCLDS)
!             WEIGHT FOR STRATIFORM CLOUD
      REAL , INTENT(OUT   ) :: LWP_STRAT(row_length, rows, NCLDS)
!             LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
      REAL , INTENT(OUT   ) :: NTOT_DIAG(row_length, rows, NCLDS)
!             DROPLET CONCENTRATION*WEIGHT
      REAL , INTENT(OUT   ) :: STRAT_LWC_DIAG(row_length, rows, NCLDS)
!             STRATIFORM LWC*WEIGHT
      REAL , INTENT(OUT   ) :: SO4_CCN_DIAG(row_length, rows, NCLDS)
!             SO4 CCN MASS CONC*COND. SAMP. WEIGHT
      REAL , INTENT(OUT   ) :: COND_SAMP_WGT(row_length, rows, NCLDS)
!             CONDITIONAL SAMPLING WEIGHT
      REAL , INTENT(OUT   ) :: NC_DIAG(row_length, rows)
!             COLUMN DROPLET CONCENTRATION * SAMPLING WEIGHT
      REAL , INTENT(OUT   ) :: NC_WEIGHT(row_length, rows)
!             COLUMN DROPLET CONCENTRATION SAMPLING WEIGHT
!
!
!
!     LOCAL VARIABLES:
      INTEGER :: I
!             LOOP VARIABLE
      INTEGER :: J
!             LOOP VARIABLE
      INTEGER :: L
!             LOOP VARIABLE
      INTEGER :: LG
!             INDEX TO GATHER
!
      LOGICAL :: L_GLACIATED_TOP(NPD_PROFILE)
!             LOGICAL FOR GLACIATED TOPS IN CONVECTIVE CLOUD.
!
      REAL :: LIQ_FRAC(NPD_PROFILE)
!             FRACTION OF LIQUID CLOUD WATER
      REAL :: LIQ_FRAC_CONV(NPD_PROFILE)
!             FRACTION OF LIQUID WATER IN CONVECTIVE CLOUD
      REAL :: T_GATHER(NPD_PROFILE)
!             GATHERED TEMPERATURE FOR LSP_FOCWWIL
      REAL :: T_LIST(NPD_PROFILE)
!             LIST OF TEMPERATURES
      REAL :: TOTAL_MASS(NPD_PROFILE)
!             TOTAL MASS IN CONVECTIVE CLOUD
      REAL :: CC_DEPTH(NPD_PROFILE)
!             DEPTH OF CONVECTIVE CLOUD
      REAL :: CONDENSED_MIX_RAT_BULK(NPD_PROFILE, 0: NPD_LAYER,NPD_CLOUD_COMPONENT)
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING BULK CLD
      REAL :: DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             DENSITY OF AIR
      REAL :: CONVECTIVE_CLOUD_LAYER(NPD_PROFILE)
!             AMOUNT OF CONVECTIVE CLOUD IN TH CURRENT LAYER
      REAL :: STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             STRATIFORM LIQUID CLOUD FRACTION (T>273K)
      REAL :: CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             CONVECTIVE LIQUID CLOUD FRACTION (T>273K)
      REAL :: NC_DIAG_G(NPD_PROFILE)
!             DIAGNOSTIC ARRAY FOR COLUMN DROPLET NUMBER (GATHERED)
      REAL :: NC_WEIGHT_G(NPD_PROFILE)
!             DIAGNOSTIC ARRAY FOR COL. DROP NO. SAMPLING WGT (GATHERED)
!
!
!     Parameters for the aggregate parametrization.
      REAL, Parameter :: a0_agg_cold = 7.5094588E-04
      REAL, Parameter :: b0_agg_cold = 5.0830326E-07
      REAL, Parameter :: a0_agg_warm = 1.3505403E-04
      REAL, Parameter :: b0_agg_warm = 2.6517429E-05
      REAL, Parameter :: t_switch    = 216.208
      REAL, Parameter :: t0_agg      = 279.5
      REAL, Parameter :: s0_agg      = 0.05
!
!
!     SET THE COMPONENTS WITHIN THE CLOUDS. IN THE UNIFIED MODEL WE
!     HAVE FOUR COMPONENTS: STRATIFORM ICE AND WATER AND CONVECTIVE
!     ICE AND WATER.
      N_CONDENSED=4
      TYPE_CONDENSED(1)=IP_CLCMP_ST_WATER
      TYPE_CONDENSED(2)=IP_CLCMP_ST_ICE
      TYPE_CONDENSED(3)=IP_CLCMP_CNV_WATER
      TYPE_CONDENSED(4)=IP_CLCMP_CNV_ICE
!
!
!
!     SET THE TOTAL AMOUNTS OF CLOUD AND THE FRACTIONS COMPRISED BY
!     CONVECTIVE AND STRATIFORM COMPONENTS.
!
!     ZERO THE AMOUNTS OF CLOUD IN THE UPPER LAYERS.
      DO I=1, N_LAYER-NCLDS
         DO L=1, N_PROFILE
            W_CLOUD(L, I)=0.0E+00
         ENDDO
      ENDDO
!
      IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT .AND.           &
          .NOT. L_CLOUD_WATER_PARTITION) THEN
!  This cloud representation not available with new cloud microphysics
!
!        THE CLOUDS ARE DIVIDED INTO MIXED-PHASE STRATIFORM AND
!        CONVECTIVE CLOUDS: LSP_FOCWWIL GIVES THE PARTITIONING BETWEEN
!        ICE AND WATER IN STRATIFORM CLOUDS AND IN CONVECTIVE CLOUD,
!        UNLESS THE OPTION TO PARTITION AS A FUNCTION OF THE LOCAL
!        TEMPERATURE IS SELECTED. WITHIN CONVECTIVE CLOUD THE LIQUID
!        WATER CONTENT IS DISTRIBUTED UNIFORMLY THROUGHOUT THE CLOUD.
!
!        CONVECTIVE CLOUD:
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
            ENDDO
         ENDDO
!
!
         IF (L_LOCAL_CNV_PARTITION) THEN
!
!           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
!           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
!           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
!           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.
!
            DO L=1, N_PROFILE
!              MIN is required since CCT may be 0 if there is no
!              convective cloud.
               L_GLACIATED_TOP(L)                                       &
                  =(T(L, MIN(N_LAYER+2-CCT(I_GATHER(L))                 &
                  , N_LAYER-NCLDS+1)) <  TM)
            ENDDO

         ELSE
!
!           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
!           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
!           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.
!
            DO L=1, N_PROFILE
               T_GATHER(L)=T(L, MIN(N_LAYER+2-CCT(I_GATHER(L))          &
                  , N_LAYER-NCLDS+1))
            ENDDO
! DEPENDS ON: lsp_focwwil
            CALL LSP_FOCWWIL(   &
                 T_GATHER     , & !INTEGER, INTENT(IN        ) :: POINTS             ! IN Number of points to be processed.
                 N_PROFILE    , & !REAL   , INTENT(IN        ) :: T(POINTS)       ! IN Temperature at this level (K).
                 LIQ_FRAC_CONV  ) !REAL   , INTENT(OUT  ) :: ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.
!
         ENDIF
!
!
         DO L=1, N_PROFILE
            TOTAL_MASS(L)=0.0E+00
         ENDDO
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
                    (CCB(LG) <= N_LAYER+1-I) ) THEN
                  TOTAL_MASS(L)=TOTAL_MASS(L)+D_MASS(L, I)
               ENDIF
            ENDDO
         ENDDO
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
                    (CCB(LG) <= N_LAYER+1-I) ) THEN
                  IF (L_LOCAL_CNV_PARTITION) THEN
!                    THE PARTITIONING IS RECALCULATED FOR EACH LAYER
!                    OTHERWISE A GENERIC VALUE IS USED.
                     LIQ_FRAC_CONV(L)=MAX(0.0E+00, MIN(1.0E+00          &
                        , 1.61E-02*(T(L, I)-TM)+8.9E-01))
!                    Do not allow ice above 0 Celsius unless the top
!                    of the cloud is glaciated and force homogeneous
!                    nucleation at -40 Celsius.
                     IF ( (T(L, I) >  TM).AND.                          &
                          (.NOT.L_GLACIATED_TOP(L)) ) THEN
                       LIQ_FRAC_CONV(L)=1.0E+00
                     ELSE IF (T(L, I) <  TM-4.0E+01) THEN
                       LIQ_FRAC_CONV(L)=0.0E+00
                     ENDIF
                  ENDIF
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)      &
                     =CCCWP(LG)*LIQ_FRAC_CONV(L)                        &
                     /(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)        &
                     =CCCWP(LG)*(1.0E+00-LIQ_FRAC_CONV(L))              &
                     /(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)      &
                     =CCCWP(LG)*LIQ_FRAC_CONV(L)                        &
                     /(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)        &
                     =CCCWP(LG)*(1.0E+00-LIQ_FRAC_CONV(L))              &
                     /(TOTAL_MASS(L)+TINY(CCCWP))
               ENDIF
            ENDDO
         ENDDO
!
!
!        STRATIFORM CLOUDS:
!
!        PARTITION BETWEEN ICE AND WATER DEPENDING ON THE
!        LOCAL TEMPERATURE.
!
         DO I=1, NCLDS
! DEPENDS ON: lsp_focwwil
            CALL LSP_FOCWWIL( &
                 T(L, N_LAYER+1-I)    , &!INTEGER, INTENT(IN   ) :: POINTS            ! IN Number of points to be processed.
                 N_PROFILE            , &!REAL   , INTENT(IN   ) :: T(POINTS)            ! IN Temperature at this level (K).
                 LIQ_FRAC               )!REAL   , INTENT(OUT  ) :: ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF (LCA_AREA(LG, I) >  EPSILON(LCA_AREA)) THEN
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)                                &
                     =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
                     *LIQ_FRAC(L)/LCA_AREA(LG, I)
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)                                  &
                     =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
                     *(1.0E+00-LIQ_FRAC(L))/LCA_AREA(LG, I)
               ELSE
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
!
               IF (LCA_BULK(LG, I) >  EPSILON(LCA_BULK)) THEN
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)                                &
                     =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
                     *LIQ_FRAC(L)/LCA_BULK(LG, I)
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)                                  &
                     =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
                     *(1.0E+00-LIQ_FRAC(L))/LCA_BULK(LG, I)
               ELSE
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
!        CLOUD FRACTIONS:
!
       IF (L_3D_CCA) THEN
         DO I=1, NCLDS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               W_CLOUD(L, N_LAYER+1-I)                                  &
                  =CCA(LG,I)+(1.0E+00-CCA(LG,I))*LCA_AREA(LG, I)
               FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)           &
                  =CCA(LG,I)/(W_CLOUD(L, N_LAYER+1-I)+TINY(CCA))
               FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_STRAT)          &
                  =1.0E+00-FRAC_CLOUD(L, N_LAYER+1-I                    &
                  , IP_CLOUD_TYPE_CONV)
            ENDDO
         ENDDO
       ELSE
         DO I=1, NCLDS
            DO L=1, N_PROFILE
              LG=I_GATHER(L)
               IF ( (I <= CCT(LG)-1).AND.(I >= CCB(LG)) ) THEN
                  W_CLOUD(L, N_LAYER+1-I)                               &
                     =CCA(LG,1)+(1.0E+00-CCA(LG,1))*LCA_AREA(LG, I)
                  FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)        &
                     =CCA(LG,1)/(W_CLOUD(L, N_LAYER+1-I)+TINY(CCA))
               ELSE
                  W_CLOUD(L, N_LAYER+1-I)=LCA_AREA(LG, I)
                  FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)        &
                     =0.0E+00
               ENDIF
               FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_STRAT)          &
                  =1.0E+00-FRAC_CLOUD(L, N_LAYER+1-I                    &
                  , IP_CLOUD_TYPE_CONV)
            ENDDO
         ENDDO
       ENDIF
!
!
!
!
      ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
!
!        HERE THE CLOUDS ARE SPLIT INTO FOUR SEPARATE TYPES.
!        THE PARTITIONING BETWEEN ICE AND WATER IS REGARDED AS
!        DETERMINING THE AREAS WITHIN THE GRID_BOX COVERED BY
!        ICE OR WATER CLOUD, RATHER THAN AS DETERMINING THE IN-CLOUD
!        MIXING RATIOS. THE GRID-BOX MEAN ICE WATER CONTENTS IN
!        STRATIFORM CLOUDS MAY BE PREDICTED BY THE ICE MICROPHYSICS
!        SCHEME OR MAY BE DETERMINED AS A FUNCTION OF THE TEMPERATURE
!        (LSP_FOCWWIL). IN CONVECTIVE CLOUDS THE PARTITIONING MAY BE
!        DONE USING THE SAME FUNCTION, LSP_FOCWWIL, BASED ON A SINGLE
!        TEMPERATURE, OR USING A PARTITION BASED ON THE LOCAL
!        TEMPERATURE.
!
!        CONVECTIVE CLOUD:
!
          DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
            TOTAL_MASS(L)=0.0E+00
         ENDDO
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
                    (CCB(LG) <= N_LAYER+1-I) ) THEN
                  TOTAL_MASS(L)=TOTAL_MASS(L)+D_MASS(L, I)
               ENDIF
            ENDDO
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
                    (CCB(LG) <= N_LAYER+1-I) ) THEN
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)      &
                     =CCCWP(LG)/(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)        &
                     =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)      &
                     =CCCWP(LG)/(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)        &
                     =CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)
               ENDIF
            ENDDO
         ENDDO
!
!        STRATIFORM CLOUDS:
!
         DO I=1, NCLDS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF (LCA_AREA(LG, I) >  EPSILON(LCA_AREA)) THEN
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)                                &
                    =(LCCWC1(LG, I)+LCCWC2(LG, I))/LCA_AREA(LG, I)
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)                                  &
                    =CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I              &
                    , IP_CLCMP_ST_WATER)
               ELSE
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
!
               IF (LCA_BULK(LG, I) >  EPSILON(LCA_BULK)) THEN
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)                                &
                    =(LCCWC1(LG, I)+LCCWC2(LG, I))/LCA_BULK(LG, I)
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)                                  &
                    =CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I              &
                    , IP_CLCMP_ST_WATER)
               ELSE
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
                    , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
!        CLOUD FRACTIONS:
!
         IF (L_LOCAL_CNV_PARTITION) THEN
!
!           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
!           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
!           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
!           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.
!
            DO L=1, N_PROFILE
!              MIN is required since CCT may be 0 if there is no
!              convective cloud.
               L_GLACIATED_TOP(L)                                       &
                  =(T(L, MIN(N_LAYER+2-CCT(I_GATHER(L))                 &
                  , N_LAYER-NCLDS+1)) <  TM)
            ENDDO

         ELSE
!
!           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
!           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
!           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.
!
            DO L=1, N_PROFILE
               T_GATHER(L)=T(L, MIN(N_LAYER+2-CCT(I_GATHER(L))          &
                  , N_LAYER-NCLDS+1))
            ENDDO
! DEPENDS ON: lsp_focwwil
            CALL LSP_FOCWWIL(&
               T_GATHER       , & !INTEGER, INTENT(IN        ) :: POINTS             ! IN Number of points to be processed.
               N_PROFILE      , & !REAL   , INTENT(IN        ) :: T(POINTS)       ! IN Temperature at this level (K).
               LIQ_FRAC_CONV    ) !REAL   , INTENT(OUT  ) :: ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.
!
         ENDIF
!
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
!
            IF (.NOT. L_CLOUD_WATER_PARTITION)                          &
!           PARTITION STRATIFORM CLOUDS USING THE LOCAL TEMPERATURE.
! DEPENDS ON: lsp_focwwil
              CALL LSP_FOCWWIL( &
               T(1, I)        , & !INTEGER, INTENT(IN        ) :: POINTS             ! IN Number of points to be processed.
               N_PROFILE      , & !REAL   , INTENT(IN        ) :: T(POINTS)       ! IN Temperature at this level (K).
               LIQ_FRAC         ) !REAL   , INTENT(OUT  ) :: ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.
!
          IF (L_3D_CCA) THEN
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
              CONVECTIVE_CLOUD_LAYER(L)=CCA(LG,N_LAYER+1-I)
            ENDDO
          ELSE
            DO L=1, N_PROFILE
            LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
                    (CCB(LG) <= N_LAYER+1-I) ) THEN
                CONVECTIVE_CLOUD_LAYER(L)=CCA(LG,1)
               ELSE
                  CONVECTIVE_CLOUD_LAYER(L)=0.0E+00
               ENDIF
            ENDDO
          ENDIF
!
            DO L=1, N_PROFILE
            LG=I_GATHER(L)
               W_CLOUD(L, I)                                            &
                  =CONVECTIVE_CLOUD_LAYER(L)                            &
                  +(1.0E+00-CONVECTIVE_CLOUD_LAYER(L))                  &
                  *LCA_AREA(LG, N_LAYER+1-I)
!
               IF (L_CLOUD_WATER_PARTITION) THEN
!  PARTITION STRATIFORM CLOUDS USING THE RATIO OF CLOUD WATER CONTENTS.
                 IF (LCA_AREA(LG, N_LAYER+1-I) >   EPSILON(LCA_AREA)    &
                    .AND.                                               &
                    (LCCWC1(LG, N_LAYER+1-I) + LCCWC2(LG, N_LAYER+1-I)) &
                     >   0.0) THEN
                   LIQ_FRAC(L) = LCCWC1(LG, N_LAYER+1-I) /              &
                    (LCCWC1(LG, N_LAYER+1-I) + LCCWC2(LG, N_LAYER+1-I))
                 ELSE
                   LIQ_FRAC(L) = 0.0E+00
                 ENDIF
               ENDIF
!
               IF (L_LOCAL_CNV_PARTITION) THEN
!
!                THE PARTITIONING BETWEEN ICE AND WATER MUST BE
!                RECALCULATED FOR THIS LAYER AS A FUNCTION OF THE
!                LOCAL TEMPERATURE, BUT ICE IS ALLOWED ABOVE THE
!                FREEZING POINT ONLY IF THE TOP OF THE CLOUD IS
!                GLACIATED.
                 LIQ_FRAC_CONV(L)=MAX(0.0E+00, MIN(1.0E+00              &
                    , 1.61E-02*(T(L, I)-TM)+8.9E-01))
!                Do not allow ice above 0 Celsius unless the top
!                of the cloud is glaciated and force homogeneous
!                nucleation at -40 Celsius.
                 IF ( (T(L, I) >  TM).AND.                              &
                      (.NOT.L_GLACIATED_TOP(L)) ) THEN
                    LIQ_FRAC_CONV(L)=1.0E+00
                 ELSE IF (T(L, I) <  TM-4.0E+01) THEN
                    LIQ_FRAC_CONV=0.0E+00
                 ENDIF

               ENDIF
!
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)                       &
                  =LIQ_FRAC(L)*(1.0E+00-CONVECTIVE_CLOUD_LAYER(L))      &
                  *LCA_AREA(LG, N_LAYER+1-I)                            &
                  /(W_CLOUD(L, I)+TINY(W_CLOUD))
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)                       &
                  =(1.0E+00-LIQ_FRAC(L))                                &
                  *(1.0E+00-CONVECTIVE_CLOUD_LAYER(L))                  &
                  *LCA_AREA(LG, N_LAYER+1-I)                            &
                  /(W_CLOUD(L, I)+TINY(W_CLOUD))
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                       &
                  =LIQ_FRAC_CONV(L)*CONVECTIVE_CLOUD_LAYER(L)           &
                  /(W_CLOUD(L, I)+TINY(W_CLOUD))
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)                       &
                  =(1.0E+00-LIQ_FRAC_CONV(L))*CONVECTIVE_CLOUD_LAYER(L) &
                  /(W_CLOUD(L, I)+TINY(W_CLOUD))
!
            ENDDO
         ENDDO
!
!
      ENDIF
!
!
!
!     EFFECTIVE RADII OF WATER CLOUDS: A MICROPHYSICAL PARAMETRIZATION
!     IS AVAILABLE; OTHERWISE STANDARD VALUES ARE USED.
!
      IF (L_MICROPHYSICS) THEN
!
!        STANDARD VALUES ARE USED FOR ICE CRYSTALS, BUT
!        A PARAMETRIZATION PROVIDED BY UMIST AND MRF
!        IS USED FOR DROPLETS.
!
!        CALCULATE THE DENSITY OF AIR.
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               DENSITY_AIR(L, I)=P(L, I)/(R*T(L, I))
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
            CC_DEPTH(L)=0.0E+00
         ENDDO
!
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
!           This loop should be safe even when convective
!           cloud is not present, since CCB should not exceed CCT.
            DO I=N_LAYER+2-CCT(LG), N_LAYER+1-CCB(LG)
               CC_DEPTH(L)=CC_DEPTH(L)+D_MASS(L, I)/DENSITY_AIR(L, I)
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
            LYING_SNOW_G(L)=LYING_SNOW(I_GATHER(L))
         ENDDO
!
         IF (NC_DIAG_FLAG) THEN
            DO I=N_LAYER+1-NCLDS, N_LAYER
               DO L=1, N_PROFILE
                  IF (T(L, I)  >=  TM) THEN
                     STRAT_LIQ_CLOUD_FRACTION(L, I)=W_CLOUD(L, I)       &
                                 *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)
                     CONV_LIQ_CLOUD_FRACTION(L, I)=W_CLOUD(L, I)        &
                                 *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)
                  ELSE
                     STRAT_LIQ_CLOUD_FRACTION(L, I)=0.0
                     CONV_LIQ_CLOUD_FRACTION(L, I)=0.0
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
! DEPENDS ON: r2_re_mrf_umist

         CALL R2_RE_MRF_UMIST(&
             N_PROFILE                , &
             N_LAYER                  , &
             NCLDS                    , &
             I_GATHER                 , &
             L_PC2                    , &
             L_AEROSOL_CCN            , &
             L_BIOMASS_CCN            , &
             L_OCFF_CCN               , &
             SEA_SALT_FILM            , &
             SEA_SALT_JET             , &
             L_SEASALT_CCN            , &
             SALT_DIM_A               , &
             SALT_DIM_B               , &
             L_USE_BIOGENIC           , &
             BIOGENIC                 , &
             BIOGENIC_DIM1            , &
             BIOGENIC_DIM2            , &
             ACCUM_SULPHATE           , &
             DISS_SULPHATE            , &
             AITKEN_SULPHATE          , &
             AGED_BMASS               , &
             CLOUD_BMASS              , &
             AGED_OCFF                , &
             CLOUD_OCFF               , &
             LYING_SNOW_G             , &
             I_CLOUD_REPRESENTATION   , &
             LAND_G                   , &
             FLANDG_G                 , &
             DENSITY_AIR              , &
             CONDENSED_MIX_RAT_BULK   , &
             CC_DEPTH                 , &
             CONDENSED_DIM_CHAR       , &
             D_MASS                   , &
             STRAT_LIQ_CLOUD_FRACTION , &
             CONV_LIQ_CLOUD_FRACTION  , &
             NC_DIAG_FLAG             , &
             NC_DIAG_G                , &
             NC_WEIGHT_G              , &
             NTOT_DIAG_G              , &
             Ntot_land                , &
             Ntot_sea                 , &
             STRAT_LWC_DIAG_G         , &
             SO4_CCN_DIAG_G           , &
             SULP_DIM1                , &
             SULP_DIM2                , &
             BMASS_DIM1               , &
             BMASS_DIM2               , &
             OCFF_DIM1                , &
             OCFF_DIM2                , &
             NPD_FIELD                , &
             NPD_PROFILE              , &
             NPD_LAYER                , &
             NPD_AEROSOL_SPECIES        )
!
!        CONSTRAIN THE SIZES OF DROPLETS TO LIE WITHIN THE RANGE OF
!        VALIDITY OF THE PARAMETRIZATION SCHEME.
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)              &
                  =MAX(CONDENSED_MIN_DIM(IP_CLCMP_ST_WATER)             &
                  , MIN(CONDENSED_MAX_DIM(IP_CLCMP_ST_WATER)            &
                  , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)))
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)             &
                  =MAX(CONDENSED_MIN_DIM(IP_CLCMP_CNV_WATER)            &
                  , MIN(CONDENSED_MAX_DIM(IP_CLCMP_CNV_WATER)           &
                  , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)))
            ENDDO
         ENDDO
!
!
!        SET MICROPHYSICAL DIAGNOSTICS. WEIGHTS FOR CLOUD CALCULATED
!        HERE ARE USED SOLELY FOR THE MICROPHYSICS AND DO NOT HAVE
!        AN INDEPENDENT MEANING.
!
         IF (WGT_CONV_FLAG) THEN
            IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_CONV(col_list(l), row_list(l), I)              &
                        =W_CLOUD(L, N_LAYER+1-I)                        &
                        *FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)
                  ENDDO
               ENDDO
            ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_CONV(col_list(l), row_list(l), I)              &
                        =W_CLOUD(L, N_LAYER+1-I)                        &
                        *FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CW)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         IF (RE_CONV_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
!                 EFFECTIVE RADII ARE GIVEN IN MICRONS.
                  RE_CONV(col_list(l), row_list(l), I)                  &
                     =CONDENSED_DIM_CHAR(L, N_LAYER+1-I                 &
                     , IP_CLCMP_CNV_WATER)                              &
                     *WGT_CONV(col_list(l), row_list(l), I)*1.0E+06
               ENDDO
            ENDDO
         ENDIF
!
         IF (WGT_STRAT_FLAG) THEN
            IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_STRAT(col_list(l), row_list(l), I)             &
                        =W_CLOUD(L, N_LAYER+1-I)                        &
                        *FRAC_CLOUD(L, N_LAYER+1-I                      &
                        , IP_CLOUD_TYPE_STRAT)
                  ENDDO
               ENDDO
            ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_STRAT(col_list(l), row_list(l), I)             &
                        =W_CLOUD(L, N_LAYER+1-I)                        &
                        *FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_SW)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         IF (RE_STRAT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
!                 EFFECTIVE RADII ARE GIVEN IN MICRONS.
                  RE_STRAT(col_list(l), row_list(l), I)                 &
                     =CONDENSED_DIM_CHAR(L, N_LAYER+1-I                 &
                     , IP_CLCMP_ST_WATER)                               &
                     *WGT_STRAT(col_list(l), row_list(l), I)*1.0E+06
               ENDDO
            ENDDO
         ENDIF

         IF (LWP_STRAT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LWP_STRAT(col_list(l), row_list(l), I)                &
                     =CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I             &
                     , IP_CLCMP_ST_WATER)*D_MASS(L, N_LAYER+1-I)        &
                     *WGT_STRAT(col_list(l), row_list(l), I)
               ENDDO
            ENDDO
         ENDIF

         IF (NC_DIAG_FLAG .AND. NC_WEIGHT_FLAG) THEN
            DO L=1, N_PROFILE
               NC_DIAG(col_list(L), row_list(L))                        &
                  =NC_DIAG_G(L)*NC_WEIGHT_G(L)
            ENDDO
         ENDIF

         IF (NC_WEIGHT_FLAG) THEN
            DO L=1, N_PROFILE
               NC_WEIGHT(col_list(L), row_list(L))=NC_WEIGHT_G(L)
            ENDDO
         ENDIF

         IF (NTOT_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  NTOT_DIAG(col_list(l), row_list(l), I)                &
                     =NTOT_DIAG_G(L, N_LAYER+1-I)                       &
                     *WGT_STRAT(col_list(l), row_list(l), I)
               ENDDO
            ENDDO
         ENDIF

         IF (STRAT_LWC_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  STRAT_LWC_DIAG(col_list(l), row_list(l), I)           &
                     =STRAT_LWC_DIAG_G(L, N_LAYER+1-I)                  &
                     *WGT_STRAT(col_list(l), row_list(l), I)
               ENDDO
            ENDDO
         ENDIF

! Non-cloud diagnostics are "weighted" by the conditional sampling
! weight COND_SAMP_WGT, but as this is 1.0 if the SW radiation is
! active, and 0.0 if it is not, there is no need to actually
! multiply by it.

         IF (COND_SAMP_WGT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  COND_SAMP_WGT(col_list(l), row_list(l), I)=1.0
               ENDDO
            ENDDO
         ENDIF

         IF (SO4_CCN_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  SO4_CCN_DIAG(col_list(l), row_list(l), I)             &
                          =SO4_CCN_DIAG_G(L, N_LAYER+1-I)
               ENDDO
            ENDDO
         ENDIF
!
!
      ELSE
!
!        ALL EFFECTIVE RADII ARE SET TO STANDARD VALUES.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)=7.E-6
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)=7.E-6
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     SET THE CHARACTERISTIC DIMENSIONS OF ICE CRYSTALS:
!
!     ICE CRYSTALS IN STRATIFORM CLOUDS:
!
      IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE) ==                         &
         IP_SLINGO_SCHRECKER_ICE) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
!        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)=30.E-6
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE) ==                    &
         IP_ICE_ADT) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
!        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
!        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
!        AT THE FREEZING LEVEL.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                &
                  =MIN(7.198755E-04                                     &
                  , EXP(5.522E-02*(T(L, I)-2.7965E+02))/9.702E+02)
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE) ==                    &
         IP_ICE_AGG_DE) THEN
!
!      Aggregate parametrization based on effective dimension.
!      In the initial form, the same approach is used for stratiform
!      and convective cloud.
!
!      The fit provided here is based on Stephan Havemann's fit of
!      Dge with temperature, consistent with David Mitchell's treatment
!      of the variation of the size distribution with temperature. The
!      parametrization of the optical properties is based on De
!      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
!      (=(2*SQRT(3)/3)*volume/projected area), which explains the
!      conversion factor. The fit to Dge is in two sections, because
!      Mitchell's relationship predicts a cusp at 216.208 K. Limits
!      of 8 and 124 microns are imposed on Dge: these are based on this
!      relationship and should be reviewed if it is changed. Note also
!      that the relationship given here is for polycrystals only.
       DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
!          Preliminary calculation of Dge.
           IF (T(L, I) < t_switch) THEN
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                  &
               = a0_agg_cold*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_cold
           ELSE
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                  &
               = a0_agg_warm*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_warm
           ENDIF
!          Limit and convert to De.
           CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                    &
             = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
               MIN(1.24E-04, MAX(8.0E-06,                               &
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)))
         ENDDO
       ENDDO
!
      ENDIF
!
!
!     ICE CRYSTALS IN CONVECTIVE CLOUDS:
!
      IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE) ==                        &
         IP_SLINGO_SCHRECKER_ICE) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
!        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)=30.E-6
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE) ==                   &
         IP_ICE_ADT) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
!        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
!        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
!        AT THE FREEZING LEVEL.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)               &
                  =MIN(7.198755E-04                                     &
                  , EXP(5.522E-02*(T(L, I)-2.7965E+02))/9.702E+02)
            ENDDO
         ENDDO
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE) ==                   &
         IP_ICE_AGG_DE) THEN
!
!      Aggregate parametrization based on effective dimension.
!      In the initial form, the same approach is used for stratiform
!      and convective cloud.
!
!      The fit provided here is based on Stephan Havemann's fit of
!      Dge with temperature, consistent with David Mitchell's treatment
!      of the variation of the size distribution with temperature. The
!      parametrization of the optical properties is based on De
!      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
!      (=(2*SQRT(3)/3)*volume/projected area), which explains the
!      conversion factor. The fit to Dge is in two sections, because
!      Mitchell's relationship predicts a cusp at 216.208 K. Limits
!      of 8 and 124 microns are imposed on Dge: these are based on this
!      relationship and should be reviewed if it is changed. Note also
!      that the relationship given here is for polycrystals only.
       DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
!          Preliminary calculation of Dge.
           IF (T(L, I) < t_switch) THEN
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                 &
               = a0_agg_cold*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_cold
           ELSE
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                 &
               = a0_agg_warm*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_warm
           ENDIF
!          Limit and convert to De.
           CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                   &
             = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
               MIN(1.24E-04, MAX(8.0E-06,                               &
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)))
         ENDDO
       ENDDO
!
!
      ENDIF
!
!
!
!     CONSTRAIN THE SIZES OF ICE CRYSTALS TO LIE WITHIN THE RANGE
!     OF VALIDITY OF THE PARAMETRIZATION SCHEME.
      DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
            CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                   &
               =MAX(CONDENSED_MIN_DIM(IP_CLCMP_ST_ICE)                  &
               , MIN(CONDENSED_MAX_DIM(IP_CLCMP_ST_ICE)                 &
               , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)))
            CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                  &
               =MAX(CONDENSED_MIN_DIM(IP_CLCMP_CNV_ICE)                 &
               , MIN(CONDENSED_MAX_DIM(IP_CLCMP_CNV_ICE)                &
               , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)))
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_SET_CLOUD_FIELD
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code to check the
!                                               range of validity of
!                                               parametrizations
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Error message for
!                                               ice corrected.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             12-06-96                Code rewritten to
!                                               include two types
!                                               of sulphate provided
!                                               by the sulphur cycle.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Climatological aerosol
!                                               model added.
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code for aerosols
!                                               generalized to allow
!                                               arbitrary combinations.
!                                               (J. M. Edwards)
!       4.5   April 1998   Option to use interactive soot in place
!                          of climatological soot.     Luke Robinson.
!                          (Repositioned more logically at 5.1)
!                                                      J. M. Edwards
!       5.1             11-04-00                The boundary layer
!                                               depth is passed to
!                                               the routine setting
!                                               the climatological
!                                               aerosol.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Set up sea-salt aerosol
!                                               (and deactivate climat-
!                                               ological sea-salt) if
!                                               required.
!                                               (A. Jones)
!       5.3             16-10-01                Switch off the
!                                               climatological
!                                               water soluble aerosol
!                                               when the sulphur
!                                               cycle is on.
!                                               (J. M. Edwards)
!       5.3             04-04-01                Include mesoscale
!                                               aerosols if required.
!                                                            S. Cusack
!       5.4             09-05-02                Use L_USE_SOOT_DIRECT
!                                               to govern the extra-
!                                               polation of aerosol
!                                               soot to the extra top
!                                               layer if required.
!                                               (A. Jones)
!       5.5             05-02-03                Include biomass aerosol
!                                               if required.  P Davison
!       5.5             10-01-03                Revision to sea-salt
!                                               density parameter.
!                                               (A. Jones)
!       5.5             21-02-03                Include mineral dust
!                                               aerosol if required.
!                                               S Woodward
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             29-09-97                Original Code
!                                               very closely based on
!                                               previous versions of
!                                               this scheme.
!                                               (J. M. Edwards)
!  4.5  12/05/98  Swap loop order in final nest of loops to
!                 improve vectorization.  RBarnes@ecmwf.int
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Impose a minimum
!                                               thickness on the
!                                               layer filled by the
!                                               boundary layer
!                                               aerosol and correct
!                                               the dimensioning of T.
!                                               (J. M. Edwards)
!       5.3             17-10-01                Restrict the height
!                                               of the BL to ensure
!                                               that at least one layer
!                                               lies in the free
!                                               troposphere in all
!                                               cases. This is needed
!                                               to allow for changes
!                                               in the range of the
!                                               tropopause.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code added for coherent
!                                               convective cloud.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check that cloud
!                                               fraction is between 0
!                                               and 1 for PC2 scheme.
!                                               (D. Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete bounds on
!                                               effective radius
!                                               removed.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Subroutine for droplet
!                                               number concentration
!                                               replaced by new function
!                                               with option to use sea-
!                                               salt to supplement
!                                               sulphate aerosol.
!                                               Treatment of convective
!                                               cloud effective radii
!                                               updated.
!                                               (A. Jones)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               NUMBER_DROPLET.
!                                               (A. Jones)
!       6.1             07-04-04                Add new variables for
!                                               column cloud droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Protect calculations
!                                               from failure in PC2
!                                               (Damian Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate column-integrated cloud droplet number.
!
! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.
!
! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       5.4             29-05-02                Original Code
!                                               (A. Jones)
!       6.1             07-04-04                Modified in accordance
!                                               with AVHRR retrievals:
!                                               only clouds >273K used,
!                                               convective clouds also
!                                               included.
!                                               (A. Jones)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             04-03-96                Original Code
!                                               (J. M. Edwards)
!                                               Parts of this code are
!                                               rather redundant. The
!                                               form of writing is for
!                                               near consistency with
!                                               HADAM3.
!
!       4.5   April 1998   Check for inconsistencies between soot
!                          spectral file and options used. L Robinson.
!       5.3     04/04/01   Include mesoscale aerosol switch when
!                          checking if aerosols are required.  S. Cusack
!       5.4     09/05/02   Include logical flag for sea-salt aerosol.
!                                                              A. Jones
!       5.5     05/02/03   Include logical for biomass smoke
!                          aerosol.               P Davison
! Description of Code:
!   5.5    21/02/03 Add logical for d mineral dust
!                                                S Woodward
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------




! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE LSP_FOCWWIL--------------------------------------------
!LL
!LL  Purpose: Calculate from temperature the Fraction Of Cloud Water
!LL           Which Is Liquid.
!LL     NOTE: Operates within range 0 to -9 deg.C based upon MRF
!LL           observational analysis. Not robust to changes in TM or T0C
!LL
!LL A.Bushell   <- programmer of some or all of previous code or changes
!LL
!LL  Model
!LL version  Date     Modification history from model version 4.0:
!LL
!LL   4.0    27/09/95 Subroutine created from in-line COMDECK.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!LL
!LL
!LL  Programming standard: Unified Model Documentation Paper No 4,
!LL                        Version 1, dated 12/9/89.
!LL
!LL  Logical component covered: Part of P26.
!LL
!LL  System task:
!LL
!LL  Documentation: Unified Model Documentation Paper No 26: Eq 26.50.
!LL
!LL  Called by components P26, P23.
!*
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_FOCWWIL( &
      T         , &!INTEGER, INTENT(IN   ) :: POINTS          ! IN Number of points to be processed.
      POINTS    , &!REAL   , INTENT(IN   ) :: T(POINTS)       ! IN Temperature at this level (K).
      ROCWWIL     )!REAL   , INTENT(OUT  ) :: ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.

      IMPLICIT NONE
          ! Input integer scalar :-
      INTEGER, INTENT(IN   ) :: POINTS          ! IN Number of points to be processed.
          ! Input real arrays :-
      REAL   , INTENT(IN   ) :: T(POINTS)       ! IN Temperature at this level (K).
          ! Updated real arrays :-
      REAL   , INTENT(OUT  ) :: ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.
!*L   External subprogram called :-
!     EXTERNAL None.
!*
!-----------------------------------------------------------------------
!  Common, then local, physical constants.
!-----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
      REAL ,PARAMETER     :: TSTART=TM ! Temperature at which ROCWWIL reaches 1.
      REAL ,PARAMETER     :: TRANGE=9.0! Temperature range over which 0 < ROCWWIL < 1.
!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
                        ! Real workspace. At end of DO loop, contains :-
      REAL :: TFOC      ! T(I) within DO loop. Allows routines to call
!                        LSP_FOCWWIL(WORK1, POINTS, WORK1) to save space
!  (b) Others.
      INTEGER :: I       ! Loop counter (horizontal field index).
!
      DO  I = 1, POINTS
!
        TFOC = T(I)
!-----------------------------------------------------------------------
!L 0. Calculate fraction of cloud water which is liquid (FL),
!L    according to equation P26.50.
!-----------------------------------------------------------------------
        IF (TFOC  <=  (TSTART - TRANGE)) THEN
!       Low temperatures, cloud water all frozen------------------------
          ROCWWIL(I) = 0.0
!
        ELSE IF (TFOC  <   TSTART) THEN
!       Intermediate temperatures---------------------------------------
          ROCWWIL(I) = (TFOC - TSTART + TRANGE) / TRANGE
!
        ELSE
!       High temperatures, cloud water all liquid-----------------------
          ROCWWIL(I) = 1.0
!
        END IF
!
      END DO ! Loop over points
!
      RETURN
      END SUBROUTINE LSP_FOCWWIL




! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Ozone set in lower
!                                               levels.
!                                               (J. M. Edwards)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!       4.5             18-05-98                Provision for treating
!                                               extra (H)(C)FCs
!                                               included.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Move HCFCs to a more
!                                               natural place in the
!                                               code.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Remove the explicit
!                                               limit on the
!                                               concentration of
!                                               water vapour.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             29-05-02                Add diagnostic call
!                                               for column-integrated
!                                               cloud droplet number.
!                                               (A. Jones)
!       5.4             25-04-02                Replace land/sea mask
!                                               with land fraction in
!                                               call to NUMBER_DROPLET
!                                               (A. Jones)
!       5.5             17-02-03                Change I_CLIM_POINTER
!                                               for hp compilation
!                                               (M.Hughes)
!  6.1   20/08/03  Code for STOCHEM feedback.  C. Johnson
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!  6.2   03-11-05   Enable HadGEM1 climatological aerosols. C. F. Durman
!                   Reworked to use switch instead of #defined. R Barnes
!  6.2   25/05/05  Convert compilation into a more universally portable
!                  form. Tom Edwards
!  6.2   15/12/05  Set negative specific humidities to zero.
!                                               (J. Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Old formulation over
!                                               sea-ice removed.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Ground temperature
!                                               set equal to that
!                                               in the middle of the
!                                               bottom layer.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.3             25-04-01                Alter the specification
!                                               of temperature on rho
!                                               levels (layer
!                                               boundaries).  S.Cusack
!       5.3             25-04-01   Gather land, sea and
!                                  sea-ice temperatures and
!                                  land fraction. Replace TFS
!                                  with general sea temperature.
!                                       (N. Gedney)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New flag L_AEROSOL_CCN
!                                               introduced to allow
!                                               inclusion of indirect
!                                               aerosol forcing alone.
!                                               Correction of comments
!                                               for LCCWC1 and LCCWC2.
!                                               Correction of level at
!                                               which temperature for
!                                               partitioning
!                                               convective homogeneously
!                                               mixed cloud is taken.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             15-09-97                A parametrization of
!                                               ice crystals with a
!                                               temperature dependedence
!                                               of the size has been
!                                               added.
!                                               Explicit checking of
!                                               the sizes of particles
!                                               for the domain of
!                                               validity of the para-
!                                               metrization has been
!                                               added.
!                                               (J. M. Edwards)
!       5.0             15-04-98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  A.C.Bushell
!       4.5             18-05-98                New option for
!                                               partitioning between
!                                               ice and water in
!                                               convective cloud
!                                               included.
!                                               (J. M. Edwards)
!       4.5             13/05/98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  S. Cusack
!       5.1             04-04-00                Remove obsolete tests
!                                               for convective cloud
!                                               and removal of very
!                                               thin cloud (no longer
!                                               required with current
!                                               solvers, but affects
!                                               bit-comparison).
!                                               (J. M. Edwards)
!       5.1             06-04-00                Correct some comments
!                                               and error messages.
!                                               (J. M. Edwards)
!       5.2             10-11-00                With local partitioning
!                                               of convective cloud
!                                               between water and ice
!                                               force homogeneous
!                                               nucleation at
!                                               -40 Celsius.
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Pass sea-salt variables
!                                               down to R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       5.3             11-10-01                Convert returned
!                                               diagnostics to 2-D
!                                               arrays.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check on small cloud
!                                               fractions and liquid for
!                                               contents for PC2 scheme.
!                                               (D. Wilson)
!       5.5             24-02-03                Addition of new ice
!                                               aggregate
!                                               parametrization.
!                                               (J. M. Edwards)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       6.1             07-04-04                Add variables for
!                                               column-droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Pass through PC2 logical
!                                               (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code to check the
!                                               range of validity of
!                                               parametrizations
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Error message for
!                                               ice corrected.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             12-06-96                Code rewritten to
!                                               include two types
!                                               of sulphate provided
!                                               by the sulphur cycle.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Climatological aerosol
!                                               model added.
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code for aerosols
!                                               generalized to allow
!                                               arbitrary combinations.
!                                               (J. M. Edwards)
!       4.5   April 1998   Option to use interactive soot in place
!                          of climatological soot.     Luke Robinson.
!                          (Repositioned more logically at 5.1)
!                                                      J. M. Edwards
!       5.1             11-04-00                The boundary layer
!                                               depth is passed to
!                                               the routine setting
!                                               the climatological
!                                               aerosol.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Set up sea-salt aerosol
!                                               (and deactivate climat-
!                                               ological sea-salt) if
!                                               required.
!                                               (A. Jones)
!       5.3             16-10-01                Switch off the
!                                               climatological
!                                               water soluble aerosol
!                                               when the sulphur
!                                               cycle is on.
!                                               (J. M. Edwards)
!       5.3             04-04-01                Include mesoscale
!                                               aerosols if required.
!                                                            S. Cusack
!       5.4             09-05-02                Use L_USE_SOOT_DIRECT
!                                               to govern the extra-
!                                               polation of aerosol
!                                               soot to the extra top
!                                               layer if required.
!                                               (A. Jones)
!       5.5             05-02-03                Include biomass aerosol
!                                               if required.  P Davison
!       5.5             10-01-03                Revision to sea-salt
!                                               density parameter.
!                                               (A. Jones)
!       5.5             21-02-03                Include mineral dust
!                                               aerosol if required.
!                                               S Woodward
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             29-09-97                Original Code
!                                               very closely based on
!                                               previous versions of
!                                               this scheme.
!                                               (J. M. Edwards)
!  4.5  12/05/98  Swap loop order in final nest of loops to
!                 improve vectorization.  RBarnes@ecmwf.int
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Impose a minimum
!                                               thickness on the
!                                               layer filled by the
!                                               boundary layer
!                                               aerosol and correct
!                                               the dimensioning of T.
!                                               (J. M. Edwards)
!       5.3             17-10-01                Restrict the height
!                                               of the BL to ensure
!                                               that at least one layer
!                                               lies in the free
!                                               troposphere in all
!                                               cases. This is needed
!                                               to allow for changes
!                                               in the range of the
!                                               tropopause.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code added for coherent
!                                               convective cloud.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check that cloud
!                                               fraction is between 0
!                                               and 1 for PC2 scheme.
!                                               (D. Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete bounds on
!                                               effective radius
!                                               removed.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Subroutine for droplet
!                                               number concentration
!                                               replaced by new function
!                                               with option to use sea-
!                                               salt to supplement
!                                               sulphate aerosol.
!                                               Treatment of convective
!                                               cloud effective radii
!                                               updated.
!                                               (A. Jones)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               NUMBER_DROPLET.
!                                               (A. Jones)
!       6.1             07-04-04                Add new variables for
!                                               column cloud droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Protect calculations
!                                               from failure in PC2
!                                               (Damian Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_RE_MRF_UMIST( &
          N_PROFILE                , &
          N_LAYER                  , &
          NCLDS                    , &
          I_GATHER                 , &
          L_PC2                    , &
          L_AEROSOL_CCN            , &
          L_BIOMASS_CCN            , &
          L_OCFF_CCN               , &
          SEA_SALT_FILM            , &
          SEA_SALT_JET             , &
          L_SEASALT_CCN            , &
          SALT_DIM_A               , &
          SALT_DIM_B               , &
          L_USE_BIOGENIC           , &
          BIOGENIC                 , &
          BIOGENIC_DIM1            , &
          BIOGENIC_DIM2            , &
          ACCUM_SULPHATE           , &
          DISS_SULPHATE            , &
          AITKEN_SULPHATE          , &
          AGED_BMASS               , &
          CLOUD_BMASS              , &
          AGED_OCFF                , &
          CLOUD_OCFF               , &
          LYING_SNOW_G             , &
          I_CLOUD_REPRESENTATION   , &
          LAND_G                   , &
          FLANDG_G                 , &
          DENSITY_AIR              , &
          CONDENSED_MIX_RATIO      , &
          CC_DEPTH                 , &
          CONDENSED_RE             , &
          D_MASS                   , &
          STRAT_LIQ_CLOUD_FRACTION , &
          CONV_LIQ_CLOUD_FRACTION  , &
          NC_DIAG_FLAG             , &
          NC_DIAG_G                , &
          NC_WEIGHT_G              , &
          NTOT_DIAG_G              , &
          Ntot_land                , &
          Ntot_sea                 , &
          STRAT_LWC_DIAG_G         , &
          SO4_CCN_DIAG_G           , &
          SULP_DIM1                , &
          SULP_DIM2                , &
          BMASS_DIM1               , &
          BMASS_DIM2               , &
          OCFF_DIM1                , &
          OCFF_DIM2                , &
          NPD_FIELD                , &
          NPD_PROFILE              , &
          NPD_LAYER                , &
          NPD_AEROSOL_SPECIES        )
!
!
!
      IMPLICIT NONE
!
!     COMDECKS INCLUDED:
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_DENSTY for subroutine SF_EXCH
! Gawd knows why this routine expects the sea to be fresh water. JRB
      REAL,PARAMETER:: RHOSEA = 1000.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
!
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  NTOT_* is the total number concentration (m-3) of cloud droplets;
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!  5.2     111000   Updated in line with Bower et al. 1994 (J. Atmos.
!                   Sci., 51, 2722-2732) and subsequent pers. comms.
!                   Droplet concentrations now as used in HadAM4.
!                                     Andy Jones
!  5.4     02/09/02 Moved THOMO here from C_LSPMIC.      Damian Wilson
!  6.2     17/11/05 Remove variables that are now in UMUI. D. Wilson
!
!     REAL,PARAMETER:: NTOT_LAND is set in UMUI
!     REAL,PARAMETER:: NTOT_SEA is set in UMUI
      REAL,PARAMETER:: KPARAM_LAND = 0.67
      REAL,PARAMETER:: KPARAM_SEA = 0.80
      REAL,PARAMETER:: DCONRE_LAND = 9.5E-06
      REAL,PARAMETER:: DCONRE_SEA = 16.0E-06
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
      REAL,PARAMETER:: THOMO = -40.0
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
! CLDCMP3A sets components of clouds for two-stream radiation code.

      ! stratiform water droplets
      INTEGER,PARAMETER:: IP_CLCMP_ST_WATER=1

      ! stratiform ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_ST_ICE=2

      ! convective water droplets
      INTEGER,PARAMETER:: IP_CLCMP_CNV_WATER=3

      ! convective ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_CNV_ICE=4

! CLDCMP3A end
! CLREPP3A defines representations of clouds in two-stream radiation
! code.

      ! all components are mixed homogeneously
      INTEGER,PARAMETER:: IP_CLOUD_HOMOGEN     = 1

      ! ice and water clouds are treated separately
      INTEGER,PARAMETER:: IP_CLOUD_ICE_WATER   = 2

      ! clouds are divided into homogeneously mixed stratiform and
      ! convective parts
      INTEGER,PARAMETER:: IP_CLOUD_CONV_STRAT  = 3

      ! clouds divided into ice and water phases and into stratiform and
      ! convective components.
      INTEGER,PARAMETER:: IP_CLOUD_CSIW        = 4

! Types of clouds (values in CLREPD3A)

      ! number of type of clouds in representation
      INTEGER :: NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)

      ! map of components contributing to types of clouds
      INTEGER :: IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT,NPD_CLOUD_REPRESENTATION)

!
!
!     DUMMY ARGUMENTS:
!
!     SIZES OF ARRAYS:
      INTEGER, INTENT(IN) :: NPD_FIELD!SIZE OF INPUT FIELDS TO THE RADIATION
      INTEGER, INTENT(IN) :: NPD_PROFILE!MAXIMUM NUMBER OF PROFILES
      INTEGER, INTENT(IN) :: NPD_LAYER!MAXIMUM NUMBER OF LAYERS
      INTEGER, INTENT(IN) :: NPD_AEROSOL_SPECIES!MAXIMUM NUMBER OF AEROSOL SPECIES
      INTEGER, INTENT(IN) :: SULP_DIM1!1ST DIMENSION OF ARRAYS OF SULPHATE
      INTEGER, INTENT(IN) :: SULP_DIM2!2ND DIMENSION OF ARRAYS OF SULPHATE
      INTEGER, INTENT(IN) :: BMASS_DIM1!1ST DIMENSION OF ARRAYS OF BIOMASS SMOKE
      INTEGER, INTENT(IN) :: BMASS_DIM2!2ND DIMENSION OF ARRAYS OF BIOMASS SMOKE
      INTEGER, INTENT(IN) :: OCFF_DIM1!1ST DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
      INTEGER, INTENT(IN) :: OCFF_DIM2!2ND DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
      INTEGER, INTENT(IN) :: SALT_DIM_A!1ST DIMENSION OF ARRAYS OF SEA-SALT
      INTEGER, INTENT(IN) :: SALT_DIM_B!2ND DIMENSION OF ARRAYS OF SEA-SALT
      INTEGER, INTENT(IN) :: BIOGENIC_DIM1!1ST DIMENSION OF BIOGENIC AEROSOL ARRAY
      INTEGER, INTENT(IN) :: BIOGENIC_DIM2!2ND DIMENSION OF BIOGENIC AEROSOL ARRAY
!
      INTEGER, INTENT(IN) :: N_PROFILE!NUMBER OF ATMOSPHERIC PROFILES
      INTEGER, INTENT(IN) :: N_LAYER!Number of layers seen in radiation
      INTEGER, INTENT(IN) :: NCLDS!NUMBER OF CLOUDY LEVELS
!
      INTEGER, INTENT(IN) :: I_GATHER(NPD_FIELD)!LIST OF POINTS TO BE GATHERED
      LOGICAL, INTENT(IN) :: LAND_G(NPD_PROFILE)!GATHERED MASK FOR LAND POINTS
      INTEGER, INTENT(IN) :: I_CLOUD_REPRESENTATION!REPRESENTATION OF CLOUDS
!
!     VARIABLES FOR PC2
      LOGICAL, INTENT(IN) :: L_PC2 !PC2 cloud scheme is in use
!
!     VARIABLES FOR AEROSOLS
      LOGICAL, INTENT(IN) :: L_AEROSOL_CCN!FLAG TO USE AEROSOLS TO FIND CCN.
      LOGICAL, INTENT(IN) :: L_SEASALT_CCN!FLAG TO USE SEA-SALT AEROSOL FOR CCN
      LOGICAL, INTENT(IN) :: L_USE_BIOGENIC!FLAG TO USE BIOGENIC AEROSOL FOR CCN
      LOGICAL, INTENT(IN) :: L_BIOMASS_CCN!FLAG TO USE BIOMASS SMOKE AEROSOL FOR CCN
      LOGICAL, INTENT(IN) :: L_OCFF_CCN!FLAG TO USE FOSSIL-FUEL ORGANIC CARBON AEROSOL FOR CCN
      LOGICAL, INTENT(IN) :: NC_DIAG_FLAG!FLAG TO DIAGNOSE COLUMN-INTEGRATED DROPLET NUMBER
!
      REAL   , INTENT(IN) :: ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)!MIXING RATIOS OF ACCUMULATION MODE SULPHATE
      REAL   , INTENT(IN) :: AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)!Mixing ratios of Aitken-mode sulphate
      REAL   , INTENT(IN) :: DISS_SULPHATE(SULP_DIM1, SULP_DIM2)!MIXING RATIOS OF DISSOLVED SULPHATE
      REAL   , INTENT(IN) :: AGED_BMASS(BMASS_DIM1, BMASS_DIM2)!MIXING RATIOS OF AGED BIOMASS SMOKE
      REAL   , INTENT(IN) :: CLOUD_BMASS(BMASS_DIM1, BMASS_DIM2)!MIXING RATIOS OF IN-CLOUD BIOMASS SMOKE
      REAL   , INTENT(IN) :: AGED_OCFF(OCFF_DIM1, OCFF_DIM2)!MIXING RATIOS OF AGED FOSSIL-FUEL ORGANIC CARBON
      REAL   , INTENT(IN) :: CLOUD_OCFF(OCFF_DIM1, OCFF_DIM2)!MIXING RATIOS OF IN-CLOUD FOSSIL-FUEL ORGANIC CARBON
      REAL   , INTENT(IN) :: SEA_SALT_FILM(SALT_DIM_A, SALT_DIM_B)!NUMBER CONCENTRATION OF FILM-MODE SEA-SALT AEROSOL
      REAL   , INTENT(IN) :: SEA_SALT_JET(SALT_DIM_A, SALT_DIM_B)!NUMBER CONCENTRATION OF JET-MODE SEA-SALT AEROSOL
      REAL   , INTENT(IN) :: BIOGENIC(BIOGENIC_DIM1, BIOGENIC_DIM2)!M.M.R. OF BIOGENIC AEROSOL
!
      REAL   , INTENT(IN) :: DENSITY_AIR(NPD_PROFILE, NPD_LAYER)!DENSITY OF AIR
!
      REAL   , INTENT(IN) :: CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT)
!             MIXING RATIOS OF CONDENSED SPECIES
      REAL   , INTENT(IN) :: CC_DEPTH(NPD_PROFILE)
!             DEPTH OF CONVECTIVE CLOUD
      REAL   , INTENT(IN) :: D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF LAYER
      REAL   , INTENT(IN) :: STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             STRATIFORM LIQUID CLOUD COVER IN LAYERS (T>273K)
      REAL   , INTENT(IN) :: CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             CONVECTIVE LIQUID CLOUD COVER IN LAYERS (T>273K)
      Real   , Intent(IN) :: Ntot_land! Number of droplets over land / m-3
      Real   , Intent(IN) :: Ntot_sea ! Number of droplets over sea / m-3

!
      REAL   , INTENT(OUT) :: CONDENSED_RE(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CONDENSED COMPONENTS OF CLOUDS
!
      REAL   , INTENT(OUT) :: NTOT_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
      REAL   , INTENT(OUT) :: STRAT_LWC_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
      REAL   , INTENT(OUT) :: SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
      REAL   , INTENT(OUT) :: NC_DIAG_G(NPD_PROFILE)
!             DIAGNOSTIC ARRAY FOR INTEGRATED DROPLET NUMBER (GATHERED)
      REAL   , INTENT(OUT) :: NC_WEIGHT_G(NPD_PROFILE)
!             DIAGNOSTIC ARRAY FOR INT DROP NO. SAMPLING WGT (GATHERED)
!
!
      REAL   , Intent(IN) :: LYING_SNOW_G(NPD_PROFILE)! GATHERED SNOW DEPTH (>5000m = LAND ICE SHEET)
      REAL   , Intent(IN) :: FLANDG_G(NPD_PROFILE)!             GATHERED global LAND FRACTION
!

!

! CLREPP3A end
!     LOCAL VARIABLES:
      INTEGER  :: I!             LOOP VARIABLE
      INTEGER  ::  L
!             LOOP VARIABLE
      INTEGER  ::  SULPHATE_PTR_A
      INTEGER  ::  SULPHATE_PTR_B
!             POINTERS FOR SULPHATE ARRAYS
      INTEGER  ::  SEASALT_PTR_A
      INTEGER  ::  SEASALT_PTR_B
!             POINTERS FOR SEA-SALT ARRAYS
      INTEGER  ::  BIOMASS_PTR_A
      INTEGER  ::  BIOMASS_PTR_B
!             POINTERS FOR BIOMASS SMOKE ARRAYS
      INTEGER  ::  OCFF_PTR_A
      INTEGER  ::  OCFF_PTR_B
!             POINTERS FOR FOSSIL-FUEL ORGANIC CARBON ARRAYS
      INTEGER  ::  BIOGENIC_PTR_A
      INTEGER  ::  BIOGENIC_PTR_B
!             POINTERS FOR BIOGENIC ARRAY
!
      REAL     ::  TOTAL_MIX_RATIO_ST(NPD_PROFILE)
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
      REAL     ::  TOTAL_MIX_RATIO_CNV(NPD_PROFILE)
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
      REAL     ::   TOTAL_STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE)
!             TOTAL STRATIFORM LIQUID CLOUD COVER (T>273K)
      REAL     ::  TOTAL_CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE)
!             TOTAL CONVECTIVE LIQUID CLOUD COVER (T>273K)
!
      REAL     ::  N_DROP(NPD_PROFILE, NPD_LAYER)
!             NUMBER DENSITY OF DROPLETS
      REAL     ::  KPARAM
!             RATIO OF CUBES OF VOLUME RADIUS TO EFFECTIVE RADIUS
      REAL     ::  TEMP
!             Temporary in calculation of effective radius
!
!     FIXED CONSTANTS OF THE PARAMETRIZATION:
      REAL  , PARAMETER   ::  DEEP_CONVECTIVE_CLOUD=5.0E+02
!             THRESHOLD VALUE FOR DEEP CONVECTIVE CLOUD
!
!
!     Subroutines called:
!      EXTERNAL                                                          &
!     &     R2_CALC_TOTAL_CLOUD_COVER                                    &
!     &   , R2_COLUMN_DROPLET_CONC
!
!     Functions called:
!      REAL                                                              &
!     &     NUMBER_DROPLET
!             Function to calculate the number of clouds droplets
!      EXTERNAL                                                          &
!     &     NUMBER_DROPLET
!
!
!
!
!     CALCULATE THE NUMBER DENSITY OF DROPLETS
!
      DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
            IF (L_AEROSOL_CCN) THEN
               SULPHATE_PTR_A=I_GATHER(L)
               SULPHATE_PTR_B=N_LAYER+1-I
            ELSE
               SULPHATE_PTR_A=1
               SULPHATE_PTR_B=1
            ENDIF
            IF (L_SEASALT_CCN) THEN
               SEASALT_PTR_A=I_GATHER(L)
               SEASALT_PTR_B=N_LAYER+1-I
            ELSE
               SEASALT_PTR_A=1
               SEASALT_PTR_B=1
            ENDIF
            IF (L_BIOMASS_CCN) THEN
               BIOMASS_PTR_A=I_GATHER(L)
               BIOMASS_PTR_B=N_LAYER+1-I
            ELSE
               BIOMASS_PTR_A=1
               BIOMASS_PTR_B=1
            ENDIF
            IF (L_OCFF_CCN) THEN
               OCFF_PTR_A=I_GATHER(L)
               OCFF_PTR_B=N_LAYER+1-I
            ELSE
               OCFF_PTR_A=1
               OCFF_PTR_B=1
            ENDIF
            IF (L_USE_BIOGENIC) THEN
               BIOGENIC_PTR_A=I_GATHER(L)
               BIOGENIC_PTR_B=N_LAYER+1-I
            ELSE
               BIOGENIC_PTR_A=1
               BIOGENIC_PTR_B=1
            ENDIF
! DEPENDS ON: number_droplet
            N_DROP(L, I)=NUMBER_DROPLET(&                         !      REAL(KIND=r8) FUNCTION NUMBER_DROPLET( &
               L_AEROSOL_CCN                                   , &!           L_AEROSOL_DROPLET, &
               .TRUE.                                          , &!           L_NH42SO4        , &
!               AITKEN_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B) , &!          !AITKEN_SULPHATE  , &
               ACCUM_SULPHATE (SULPHATE_PTR_A, SULPHATE_PTR_B) , &!           ACCUM_SULPHATE   , &
               DISS_SULPHATE  (SULPHATE_PTR_A, SULPHATE_PTR_B) , &!           DISS_SULPHATE    , &
               L_SEASALT_CCN                                   , &!           L_SEASALT_CCN    , &
               SEA_SALT_FILM  (SEASALT_PTR_A, SEASALT_PTR_B)   , &!           SEA_SALT_FILM    , &
               SEA_SALT_JET   (SEASALT_PTR_A, SEASALT_PTR_B)   , &!           SEA_SALT_JET     , &
               L_USE_BIOGENIC                                  , &!           L_BIOGENIC_CCN   , &
               BIOGENIC       (BIOGENIC_PTR_A, BIOGENIC_PTR_B) , &!           BIOGENIC         , &
               L_BIOMASS_CCN                                   , &!           L_BIOMASS_CCN    , &
               AGED_BMASS     (BIOMASS_PTR_A, BIOMASS_PTR_B)   , &!           BIOMASS_AGED     , &
               CLOUD_BMASS    (BIOMASS_PTR_A, BIOMASS_PTR_B)   , &!           BIOMASS_CLOUD    , &
               L_OCFF_CCN                                      , &!           L_OCFF_CCN       , &
               AGED_OCFF      (OCFF_PTR_A, OCFF_PTR_B)         , &!           OCFF_AGED        , &
               CLOUD_OCFF     (OCFF_PTR_A, OCFF_PTR_B)         , &!           OCFF_CLOUD       , &
               DENSITY_AIR    (L, I)                           , &!           DENSITY_AIR      , &
               LYING_SNOW_G   (L)                              , &!           SNOW_DEPTH       , &
               FLANDG_G       (L)                              , &!           LAND_FRACT       , &
               Ntot_land                                       , &!           NTOT_LAND        , &
               Ntot_sea                                          )!           NTOT_SEA           )
         ENDDO
      ENDDO

!  Diagnose column-integrated cloud droplet number if required.

      IF (NC_DIAG_FLAG) THEN

! DEPENDS ON: r2_calc_total_cloud_cover
         CALL R2_CALC_TOTAL_CLOUD_COVER(   &
           N_PROFILE                     , &!INTEGER , INTENT(IN ):: N_PROFILE !NUMBER OF PROFILES
           N_LAYER                       , &!INTEGER , INTENT(IN ):: N_LAYER   !Number of layers seen in radiation
           NCLDS                         , &!INTEGER , INTENT(IN ):: NCLDS     !NUMBER OF CLOUDY LAYERS
           I_CLOUD_REPRESENTATION        , &!INTEGER , INTENT(IN ):: I_CLOUD   !CLOUD SCHEME EMPLOYED
           STRAT_LIQ_CLOUD_FRACTION      , &!REAL    , INTENT(IN ):: W_CLOUD_IN(NPD_PROFILE, NPD_LAYER)!CLOUD AMOUNTS
           TOTAL_STRAT_LIQ_CLOUD_FRACTION, &!REAL    , INTENT(OUT):: TOTAL_CLOUD_COVER(NPD_PROFILE) !TOTAL CLOUD COVER
           NPD_PROFILE                   , &!INTEGER , INTENT(IN )::  NPD_PROFILE!MAXIMUM NUMBER OF PROFILES
           NPD_LAYER                       )!INTEGER , INTENT(IN )::  NPD_LAYER  !MAXIMUM NUMBER OF LAYERS

! DEPENDS ON: r2_calc_total_cloud_cover
         CALL R2_CALC_TOTAL_CLOUD_COVER(   &
           N_PROFILE                     , &!INTEGER , INTENT(IN ):: N_PROFILE !NUMBER OF PROFILES
           N_LAYER                       , &!INTEGER , INTENT(IN ):: N_LAYER   !Number of layers seen in radiation
           NCLDS                         , &!INTEGER , INTENT(IN ):: NCLDS     !NUMBER OF CLOUDY LAYERS
           I_CLOUD_REPRESENTATION        , &!INTEGER , INTENT(IN ):: I_CLOUD   !CLOUD SCHEME EMPLOYED
           CONV_LIQ_CLOUD_FRACTION       , &!REAL    , INTENT(IN ):: W_CLOUD_IN(NPD_PROFILE, NPD_LAYER)!CLOUD AMOUNTS
           TOTAL_CONV_LIQ_CLOUD_FRACTION , &!REAL    , INTENT(OUT):: TOTAL_CLOUD_COVER(NPD_PROFILE) !TOTAL CLOUD COVER
           NPD_PROFILE                   , &!INTEGER , INTENT(IN )::  NPD_PROFILE!MAXIMUM NUMBER OF PROFILES
           NPD_LAYER                       )!INTEGER , INTENT(IN )::  NPD_LAYER  !MAXIMUM NUMBER OF LAYERS

! DEPENDS ON: r2_column_droplet_conc

         CALL R2_COLUMN_DROPLET_CONC(      & !
           NPD_PROFILE                   , & ! INTEGER, INTENT(IN   ) :: NPD_PROFILE!Maximum number of profiles
           NPD_LAYER                     , & ! INTEGER, INTENT(IN   ) :: NPD_LAYER  !Maximum number of layers
           N_PROFILE                     , & ! INTEGER, INTENT(IN   ) :: N_PROFILE  !Number of atmospheric profiles
           N_LAYER                       , & ! INTEGER, INTENT(IN   ) :: N_LAYER    !Number of layers seen in radiation
           NCLDS                         , & ! INTEGER, INTENT(IN   ) :: NCLDS      !Number of cloudy layers
           STRAT_LIQ_CLOUD_FRACTION      , & ! REAL   , INTENT(IN   ) :: STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
           TOTAL_STRAT_LIQ_CLOUD_FRACTION, & ! REAL   , INTENT(IN   ) :: TOTAL_STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE)
           CONV_LIQ_CLOUD_FRACTION       , & ! REAL   , INTENT(IN   ) :: CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
           TOTAL_CONV_LIQ_CLOUD_FRACTION , & ! REAL   , INTENT(IN   ) :: TOTAL_CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE)
           N_DROP                        , & ! REAL   , INTENT(IN   ) :: N_DROP(NPD_PROFILE, NPD_LAYER)
           D_MASS                        , & ! REAL   , INTENT(IN   ) :: D_MASS(NPD_PROFILE, NPD_LAYER)
           DENSITY_AIR                   , & ! REAL   , INTENT(IN   ) :: DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
           NC_DIAG_G                     , & ! REAL   , INTENT(OUT  ) :: NC_DIAG(NPD_PROFILE)
           NC_WEIGHT_G                      )! REAL   , INTENT(OUT  ) :: NC_WEIGHT(NPD_PROFILE)

      ENDIF

!  Diagnose SO4 aerosol concentrations. Mass mixing ratio of ammonium
!  sulphate is converted to microgrammes of the sulphate ion per m3
!  for diagnostic purposes.

      IF (L_AEROSOL_CCN) THEN
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               SO4_CCN_DIAG_G(L, I)=                                    &
                       (AITKEN_SULPHATE(I_GATHER(L), N_LAYER+1-I)      &
                        +ACCUM_SULPHATE(I_GATHER(L), N_LAYER+1-I)      &
                         +DISS_SULPHATE(I_GATHER(L), N_LAYER+1-I))     &
                      * DENSITY_AIR(L, I) * (96./132.) * 1.0E+09
            ENDDO
         ENDDO
      ENDIF
!
      DO I=N_LAYER+1-NCLDS, N_LAYER
!
!        FIND THE TOTAL MIXING RATIO OF WATER SUBSTANCE IN THE CLOUD
!        AS IMPLIED BY THE REPRESENTATION.
         IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
            DO L=1, N_PROFILE
               TOTAL_MIX_RATIO_ST(L)                                    &
                 =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_WATER)         &
                 +CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_ICE)
               TOTAL_MIX_RATIO_CNV(L)                                   &
                 =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_WATER)        &
                 +CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_ICE)
            ENDDO
         ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
            DO L=1, N_PROFILE
               TOTAL_MIX_RATIO_ST(L)                                    &
                 =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_WATER)
               TOTAL_MIX_RATIO_CNV(L)                                   &
                 =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_WATER)
            ENDDO
         ENDIF
         DO L=1, N_PROFILE
            IF (LAND_G(L)) THEN
               KPARAM=KPARAM_LAND
            ELSE
               KPARAM=KPARAM_SEA
            ENDIF
            IF (.NOT. L_PC2) THEN
              CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)                    &
              =(3.0E+00*TOTAL_MIX_RATIO_CNV(L)*DENSITY_AIR(L, I)       &
              /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))             &
              **(1.0E+00/3.0E+00)
              CONDENSED_RE(L, I, IP_CLCMP_ST_WATER)                     &
              =(3.0E+00*TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)        &
              /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))             &
              **(1.0E+00/3.0E+00)
            ELSE
              TEMP                                                      &
              =(3.0E+00*TOTAL_MIX_RATIO_CNV(L)*DENSITY_AIR(L, I)       &
              /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))
              IF (TEMP  >=  0.0) THEN
                CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)                  &
               = TEMP**(1.0E+00/3.0E+00)
              ELSE
                CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER) = 0.0
              END IF
              TEMP                                                      &
              =(3.0E+00*TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)        &
              /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))
              IF (TEMP  >=  0.0) THEN
                CONDENSED_RE(L, I, IP_CLCMP_ST_WATER)                   &
               = TEMP**(1.0E+00/3.0E+00)
              ELSE
                CONDENSED_RE(L, I, IP_CLCMP_ST_WATER) = 0.0
              END IF
            END IF
         ENDDO
         DO L=1, N_PROFILE
            NTOT_DIAG_G(L, I)=N_DROP(L, I)*1.0E-06
            STRAT_LWC_DIAG_G(L, I)                                      &
              =TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)*1.0E03
         ENDDO
      ENDDO
!
!     RESET THE EFFECTIVE RADII FOR DEEP CONVECTIVE CLOUDS.
      DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
            IF (LAND_G(L)) THEN
               IF (CC_DEPTH(L)  >   DEEP_CONVECTION_LIMIT_LAND) THEN
                  CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)=DCONRE_LAND
               ELSE
                  IF (CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)            &
                                                   >   DCONRE_LAND)    &
                 CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)=DCONRE_LAND
               ENDIF
            ELSE
               IF (CC_DEPTH(L)  >   DEEP_CONVECTION_LIMIT_SEA) THEN
                  CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)=DCONRE_SEA
               ELSE
                  IF (CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)            &
                                                   >   DCONRE_SEA)     &
                 CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)=DCONRE_SEA
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_RE_MRF_UMIST
!+ Subroutine to calculate column-integrated cloud droplet number.
!
! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.
!
! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       5.4             29-05-02                Original Code
!                                               (A. Jones)
!       6.1             07-04-04                Modified in accordance
!                                               with AVHRR retrievals:
!                                               only clouds >273K used,
!                                               convective clouds also
!                                               included.
!                                               (A. Jones)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             04-03-96                Original Code
!                                               (J. M. Edwards)
!                                               Parts of this code are
!                                               rather redundant. The
!                                               form of writing is for
!                                               near consistency with
!                                               HADAM3.
!
!       4.5   April 1998   Check for inconsistencies between soot
!                          spectral file and options used. L Robinson.
!       5.3     04/04/01   Include mesoscale aerosol switch when
!                          checking if aerosols are required.  S. Cusack
!       5.4     09/05/02   Include logical flag for sea-salt aerosol.
!                                                              A. Jones
!       5.5     05/02/03   Include logical for biomass smoke
!                          aerosol.               P Davison
! Description of Code:
!   5.5    21/02/03 Add logical for d mineral dust
!                                                S Woodward
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------





! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Ozone set in lower
!                                               levels.
!                                               (J. M. Edwards)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!       4.5             18-05-98                Provision for treating
!                                               extra (H)(C)FCs
!                                               included.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Move HCFCs to a more
!                                               natural place in the
!                                               code.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Remove the explicit
!                                               limit on the
!                                               concentration of
!                                               water vapour.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             29-05-02                Add diagnostic call
!                                               for column-integrated
!                                               cloud droplet number.
!                                               (A. Jones)
!       5.4             25-04-02                Replace land/sea mask
!                                               with land fraction in
!                                               call to NUMBER_DROPLET
!                                               (A. Jones)
!       5.5             17-02-03                Change I_CLIM_POINTER
!                                               for hp compilation
!                                               (M.Hughes)
!  6.1   20/08/03  Code for STOCHEM feedback.  C. Johnson
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!  6.2   03-11-05   Enable HadGEM1 climatological aerosols. C. F. Durman
!                   Reworked to use switch instead of #defined. R Barnes
!  6.2   25/05/05  Convert compilation into a more universally portable
!                  form. Tom Edwards
!  6.2   15/12/05  Set negative specific humidities to zero.
!                                               (J. Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Old formulation over
!                                               sea-ice removed.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Ground temperature
!                                               set equal to that
!                                               in the middle of the
!                                               bottom layer.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.3             25-04-01                Alter the specification
!                                               of temperature on rho
!                                               levels (layer
!                                               boundaries).  S.Cusack
!       5.3             25-04-01   Gather land, sea and
!                                  sea-ice temperatures and
!                                  land fraction. Replace TFS
!                                  with general sea temperature.
!                                       (N. Gedney)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New flag L_AEROSOL_CCN
!                                               introduced to allow
!                                               inclusion of indirect
!                                               aerosol forcing alone.
!                                               Correction of comments
!                                               for LCCWC1 and LCCWC2.
!                                               Correction of level at
!                                               which temperature for
!                                               partitioning
!                                               convective homogeneously
!                                               mixed cloud is taken.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             15-09-97                A parametrization of
!                                               ice crystals with a
!                                               temperature dependedence
!                                               of the size has been
!                                               added.
!                                               Explicit checking of
!                                               the sizes of particles
!                                               for the domain of
!                                               validity of the para-
!                                               metrization has been
!                                               added.
!                                               (J. M. Edwards)
!       5.0             15-04-98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  A.C.Bushell
!       4.5             18-05-98                New option for
!                                               partitioning between
!                                               ice and water in
!                                               convective cloud
!                                               included.
!                                               (J. M. Edwards)
!       4.5             13/05/98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  S. Cusack
!       5.1             04-04-00                Remove obsolete tests
!                                               for convective cloud
!                                               and removal of very
!                                               thin cloud (no longer
!                                               required with current
!                                               solvers, but affects
!                                               bit-comparison).
!                                               (J. M. Edwards)
!       5.1             06-04-00                Correct some comments
!                                               and error messages.
!                                               (J. M. Edwards)
!       5.2             10-11-00                With local partitioning
!                                               of convective cloud
!                                               between water and ice
!                                               force homogeneous
!                                               nucleation at
!                                               -40 Celsius.
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Pass sea-salt variables
!                                               down to R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       5.3             11-10-01                Convert returned
!                                               diagnostics to 2-D
!                                               arrays.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check on small cloud
!                                               fractions and liquid for
!                                               contents for PC2 scheme.
!                                               (D. Wilson)
!       5.5             24-02-03                Addition of new ice
!                                               aggregate
!                                               parametrization.
!                                               (J. M. Edwards)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       6.1             07-04-04                Add variables for
!                                               column-droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Pass through PC2 logical
!                                               (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code to check the
!                                               range of validity of
!                                               parametrizations
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Error message for
!                                               ice corrected.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             12-06-96                Code rewritten to
!                                               include two types
!                                               of sulphate provided
!                                               by the sulphur cycle.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Climatological aerosol
!                                               model added.
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code for aerosols
!                                               generalized to allow
!                                               arbitrary combinations.
!                                               (J. M. Edwards)
!       4.5   April 1998   Option to use interactive soot in place
!                          of climatological soot.     Luke Robinson.
!                          (Repositioned more logically at 5.1)
!                                                      J. M. Edwards
!       5.1             11-04-00                The boundary layer
!                                               depth is passed to
!                                               the routine setting
!                                               the climatological
!                                               aerosol.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Set up sea-salt aerosol
!                                               (and deactivate climat-
!                                               ological sea-salt) if
!                                               required.
!                                               (A. Jones)
!       5.3             16-10-01                Switch off the
!                                               climatological
!                                               water soluble aerosol
!                                               when the sulphur
!                                               cycle is on.
!                                               (J. M. Edwards)
!       5.3             04-04-01                Include mesoscale
!                                               aerosols if required.
!                                                            S. Cusack
!       5.4             09-05-02                Use L_USE_SOOT_DIRECT
!                                               to govern the extra-
!                                               polation of aerosol
!                                               soot to the extra top
!                                               layer if required.
!                                               (A. Jones)
!       5.5             05-02-03                Include biomass aerosol
!                                               if required.  P Davison
!       5.5             10-01-03                Revision to sea-salt
!                                               density parameter.
!                                               (A. Jones)
!       5.5             21-02-03                Include mineral dust
!                                               aerosol if required.
!                                               S Woodward
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             29-09-97                Original Code
!                                               very closely based on
!                                               previous versions of
!                                               this scheme.
!                                               (J. M. Edwards)
!  4.5  12/05/98  Swap loop order in final nest of loops to
!                 improve vectorization.  RBarnes@ecmwf.int
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Impose a minimum
!                                               thickness on the
!                                               layer filled by the
!                                               boundary layer
!                                               aerosol and correct
!                                               the dimensioning of T.
!                                               (J. M. Edwards)
!       5.3             17-10-01                Restrict the height
!                                               of the BL to ensure
!                                               that at least one layer
!                                               lies in the free
!                                               troposphere in all
!                                               cases. This is needed
!                                               to allow for changes
!                                               in the range of the
!                                               tropopause.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code added for coherent
!                                               convective cloud.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check that cloud
!                                               fraction is between 0
!                                               and 1 for PC2 scheme.
!                                               (D. Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete bounds on
!                                               effective radius
!                                               removed.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Subroutine for droplet
!                                               number concentration
!                                               replaced by new function
!                                               with option to use sea-
!                                               salt to supplement
!                                               sulphate aerosol.
!                                               Treatment of convective
!                                               cloud effective radii
!                                               updated.
!                                               (A. Jones)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               NUMBER_DROPLET.
!                                               (A. Jones)
!       6.1             07-04-04                Add new variables for
!                                               column cloud droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Protect calculations
!                                               from failure in PC2
!                                               (Damian Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate column-integrated cloud droplet number.
!
! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.
!
! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       5.4             29-05-02                Original Code
!                                               (A. Jones)
!       6.1             07-04-04                Modified in accordance
!                                               with AVHRR retrievals:
!                                               only clouds >273K used,
!                                               convective clouds also
!                                               included.
!                                               (A. Jones)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_COLUMN_DROPLET_CONC(&!
       NPD_PROFILE                    , &! INTEGER, INTENT(IN        ) :: NPD_PROFILE!Maximum number of profiles
       NPD_LAYER                      , &! INTEGER, INTENT(IN        ) :: NPD_LAYER  !Maximum number of layers
       N_PROFILE                      , &! INTEGER, INTENT(IN        ) :: N_PROFILE  !Number of atmospheric profiles
       N_LAYER                        , &! INTEGER, INTENT(IN        ) :: N_LAYER        !Number of layers seen in radiation
       NCLDS                          , &! INTEGER, INTENT(IN        ) :: NCLDS        !Number of cloudy layers
       STRAT_LIQ_CLOUD_FRACTION       , &! REAL   , INTENT(IN        ) :: STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
       TOTAL_STRAT_LIQ_CLOUD_FRACTION , &! REAL   , INTENT(IN        ) :: TOTAL_STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE)
       CONV_LIQ_CLOUD_FRACTION        , &! REAL   , INTENT(IN        ) :: CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
       TOTAL_CONV_LIQ_CLOUD_FRACTION  , &! REAL   , INTENT(IN        ) :: TOTAL_CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE)
       N_DROP                         , &! REAL   , INTENT(IN        ) :: N_DROP(NPD_PROFILE, NPD_LAYER)
       D_MASS                         , &! REAL   , INTENT(IN        ) :: D_MASS(NPD_PROFILE, NPD_LAYER)
       DENSITY_AIR                    , &! REAL   , INTENT(IN        ) :: DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
       NC_DIAG                        , &! REAL   , INTENT(OUT  ) :: NC_DIAG(NPD_PROFILE)
       NC_WEIGHT                        )! REAL   , INTENT(OUT  ) :: NC_WEIGHT(NPD_PROFILE)
!
!
!
      IMPLICIT NONE
!
!
!
!  Input variables:
!
      INTEGER, INTENT(IN   ) :: NPD_PROFILE!Maximum number of profiles
      INTEGER, INTENT(IN   ) :: NPD_LAYER  !Maximum number of layers
      INTEGER, INTENT(IN   ) :: N_PROFILE  !Number of atmospheric profiles
      INTEGER, INTENT(IN   ) :: N_LAYER    !Number of layers seen in radiation
      INTEGER, INTENT(IN   ) :: NCLDS      !Number of cloudy layers
!
      REAL   , INTENT(IN   ) :: STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             Stratiform liquid (T>273K) cloud cover in layers
      REAL   , INTENT(IN   ) :: TOTAL_STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE)
!             Total liquid (T>273K) stratiform cloud cover
      REAL   , INTENT(IN   ) :: CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             Convective liquid (T>273K) cloud cover in layers
      REAL   , INTENT(IN   ) :: TOTAL_CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE)
!             Total liquid (T>273K) convective cloud cover
      REAL   , INTENT(IN   ) :: N_DROP(NPD_PROFILE, NPD_LAYER)
!             Number concentration of cloud droplets (m-3)
      REAL   , INTENT(IN   ) :: D_MASS(NPD_PROFILE, NPD_LAYER)
!             Mass thickness of layer (kg m-2)
      REAL   , INTENT(IN   ) :: DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             Air density (kg m-3)
!
!
!  Output variables:
!
      REAL   , INTENT(OUT  ) :: NC_DIAG(NPD_PROFILE)
!             Column-integrated droplet number diagnostic (m-2)
      REAL   , INTENT(OUT  ) :: NC_WEIGHT(NPD_PROFILE)
!             Weighting factor for column droplet number
!
!
!  Local variables:
!
      INTEGER :: I, L !Loop counters
!
      REAL    :: N_SUM_S
      REAL    :: N_SUM_C
!             Temporary sums
      REAL    :: NCOL_S
      REAL    :: WGT_S
      REAL    :: NCOL_C
      REAL    :: WGT_C
!             N-column values and weights for stratiform
!             and convective clouds separately.
!
!
!
      DO L=1, N_PROFILE

         N_SUM_S=0.0
         DO I=N_LAYER+1-NCLDS, N_LAYER
            N_SUM_S=N_SUM_S+(STRAT_LIQ_CLOUD_FRACTION(L, I)             &
                     *N_DROP(L, I)*D_MASS(L, I)/DENSITY_AIR(L, I))
         ENDDO

         N_SUM_C=0.0
         DO I=N_LAYER+1-NCLDS, N_LAYER
            N_SUM_C=N_SUM_C+(CONV_LIQ_CLOUD_FRACTION(L, I)              &
                     *N_DROP(L, I)*D_MASS(L, I)/DENSITY_AIR(L, I))
         ENDDO

         IF (TOTAL_STRAT_LIQ_CLOUD_FRACTION(L)  >   0.0) THEN
            NCOL_S=N_SUM_S/TOTAL_STRAT_LIQ_CLOUD_FRACTION(L)
            WGT_S=TOTAL_STRAT_LIQ_CLOUD_FRACTION(L)
         ELSE
            NCOL_S=0.0
            WGT_S=0.0
         ENDIF

         IF (TOTAL_CONV_LIQ_CLOUD_FRACTION(L)  >   0.0) THEN
            NCOL_C=N_SUM_C/TOTAL_CONV_LIQ_CLOUD_FRACTION(L)
            WGT_C=TOTAL_CONV_LIQ_CLOUD_FRACTION(L)
         ELSE
            NCOL_C=0.0
            WGT_C=0.0
         ENDIF

         IF ((WGT_S+WGT_C)  >   0.0) THEN
            NC_DIAG(L)=((NCOL_S*WGT_S)+(NCOL_C*WGT_C))/(WGT_S+WGT_C)
            NC_WEIGHT(L)=WGT_S+WGT_C
         ELSE
            NC_DIAG(L)=0.0
            NC_WEIGHT(L)=0.0
         ENDIF

      ENDDO
!
!     Note: weighting is done later in R2_SET_CLOUD_FIELD
!
!
!
      RETURN
      END SUBROUTINE R2_COLUMN_DROPLET_CONC
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             04-03-96                Original Code
!                                               (J. M. Edwards)
!                                               Parts of this code are
!                                               rather redundant. The
!                                               form of writing is for
!                                               near consistency with
!                                               HADAM3.
!
!       4.5   April 1998   Check for inconsistencies between soot
!                          spectral file and options used. L Robinson.
!       5.3     04/04/01   Include mesoscale aerosol switch when
!                          checking if aerosols are required.  S. Cusack
!       5.4     09/05/02   Include logical flag for sea-salt aerosol.
!                                                              A. Jones
!       5.5     05/02/03   Include logical for biomass smoke
!                          aerosol.               P Davison
! Description of Code:
!   5.5    21/02/03 Add logical for d mineral dust
!                                                S Woodward
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------





! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Ozone set in lower
!                                               levels.
!                                               (J. M. Edwards)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!       4.5             18-05-98                Provision for treating
!                                               extra (H)(C)FCs
!                                               included.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Move HCFCs to a more
!                                               natural place in the
!                                               code.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Remove the explicit
!                                               limit on the
!                                               concentration of
!                                               water vapour.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             29-05-02                Add diagnostic call
!                                               for column-integrated
!                                               cloud droplet number.
!                                               (A. Jones)
!       5.4             25-04-02                Replace land/sea mask
!                                               with land fraction in
!                                               call to NUMBER_DROPLET
!                                               (A. Jones)
!       5.5             17-02-03                Change I_CLIM_POINTER
!                                               for hp compilation
!                                               (M.Hughes)
!  6.1   20/08/03  Code for STOCHEM feedback.  C. Johnson
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!  6.2   03-11-05   Enable HadGEM1 climatological aerosols. C. F. Durman
!                   Reworked to use switch instead of #defined. R Barnes
!  6.2   25/05/05  Convert compilation into a more universally portable
!                  form. Tom Edwards
!  6.2   15/12/05  Set negative specific humidities to zero.
!                                               (J. Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Old formulation over
!                                               sea-ice removed.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Ground temperature
!                                               set equal to that
!                                               in the middle of the
!                                               bottom layer.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.3             25-04-01                Alter the specification
!                                               of temperature on rho
!                                               levels (layer
!                                               boundaries).  S.Cusack
!       5.3             25-04-01   Gather land, sea and
!                                  sea-ice temperatures and
!                                  land fraction. Replace TFS
!                                  with general sea temperature.
!                                       (N. Gedney)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New flag L_AEROSOL_CCN
!                                               introduced to allow
!                                               inclusion of indirect
!                                               aerosol forcing alone.
!                                               Correction of comments
!                                               for LCCWC1 and LCCWC2.
!                                               Correction of level at
!                                               which temperature for
!                                               partitioning
!                                               convective homogeneously
!                                               mixed cloud is taken.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             15-09-97                A parametrization of
!                                               ice crystals with a
!                                               temperature dependedence
!                                               of the size has been
!                                               added.
!                                               Explicit checking of
!                                               the sizes of particles
!                                               for the domain of
!                                               validity of the para-
!                                               metrization has been
!                                               added.
!                                               (J. M. Edwards)
!       5.0             15-04-98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  A.C.Bushell
!       4.5             18-05-98                New option for
!                                               partitioning between
!                                               ice and water in
!                                               convective cloud
!                                               included.
!                                               (J. M. Edwards)
!       4.5             13/05/98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  S. Cusack
!       5.1             04-04-00                Remove obsolete tests
!                                               for convective cloud
!                                               and removal of very
!                                               thin cloud (no longer
!                                               required with current
!                                               solvers, but affects
!                                               bit-comparison).
!                                               (J. M. Edwards)
!       5.1             06-04-00                Correct some comments
!                                               and error messages.
!                                               (J. M. Edwards)
!       5.2             10-11-00                With local partitioning
!                                               of convective cloud
!                                               between water and ice
!                                               force homogeneous
!                                               nucleation at
!                                               -40 Celsius.
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Pass sea-salt variables
!                                               down to R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       5.3             11-10-01                Convert returned
!                                               diagnostics to 2-D
!                                               arrays.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check on small cloud
!                                               fractions and liquid for
!                                               contents for PC2 scheme.
!                                               (D. Wilson)
!       5.5             24-02-03                Addition of new ice
!                                               aggregate
!                                               parametrization.
!                                               (J. M. Edwards)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       6.1             07-04-04                Add variables for
!                                               column-droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Pass through PC2 logical
!                                               (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code to check the
!                                               range of validity of
!                                               parametrizations
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Error message for
!                                               ice corrected.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             12-06-96                Code rewritten to
!                                               include two types
!                                               of sulphate provided
!                                               by the sulphur cycle.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Climatological aerosol
!                                               model added.
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code for aerosols
!                                               generalized to allow
!                                               arbitrary combinations.
!                                               (J. M. Edwards)
!       4.5   April 1998   Option to use interactive soot in place
!                          of climatological soot.     Luke Robinson.
!                          (Repositioned more logically at 5.1)
!                                                      J. M. Edwards
!       5.1             11-04-00                The boundary layer
!                                               depth is passed to
!                                               the routine setting
!                                               the climatological
!                                               aerosol.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Set up sea-salt aerosol
!                                               (and deactivate climat-
!                                               ological sea-salt) if
!                                               required.
!                                               (A. Jones)
!       5.3             16-10-01                Switch off the
!                                               climatological
!                                               water soluble aerosol
!                                               when the sulphur
!                                               cycle is on.
!                                               (J. M. Edwards)
!       5.3             04-04-01                Include mesoscale
!                                               aerosols if required.
!                                                            S. Cusack
!       5.4             09-05-02                Use L_USE_SOOT_DIRECT
!                                               to govern the extra-
!                                               polation of aerosol
!                                               soot to the extra top
!                                               layer if required.
!                                               (A. Jones)
!       5.5             05-02-03                Include biomass aerosol
!                                               if required.  P Davison
!       5.5             10-01-03                Revision to sea-salt
!                                               density parameter.
!                                               (A. Jones)
!       5.5             21-02-03                Include mineral dust
!                                               aerosol if required.
!                                               S Woodward
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             29-09-97                Original Code
!                                               very closely based on
!                                               previous versions of
!                                               this scheme.
!                                               (J. M. Edwards)
!  4.5  12/05/98  Swap loop order in final nest of loops to
!                 improve vectorization.  RBarnes@ecmwf.int
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Impose a minimum
!                                               thickness on the
!                                               layer filled by the
!                                               boundary layer
!                                               aerosol and correct
!                                               the dimensioning of T.
!                                               (J. M. Edwards)
!       5.3             17-10-01                Restrict the height
!                                               of the BL to ensure
!                                               that at least one layer
!                                               lies in the free
!                                               troposphere in all
!                                               cases. This is needed
!                                               to allow for changes
!                                               in the range of the
!                                               tropopause.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code added for coherent
!                                               convective cloud.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check that cloud
!                                               fraction is between 0
!                                               and 1 for PC2 scheme.
!                                               (D. Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_CALC_TOTAL_CLOUD_COVER( &
          N_PROFILE         , &!INTEGER , INTENT(IN ):: N_PROFILE !NUMBER OF PROFILES
          N_LAYER           , &!INTEGER , INTENT(IN ):: N_LAYER   !Number of layers seen in radiation
          NCLDS             , &!INTEGER , INTENT(IN ):: NCLDS          !NUMBER OF CLOUDY LAYERS
          I_CLOUD           , &!INTEGER , INTENT(IN ):: I_CLOUD   !CLOUD SCHEME EMPLOYED
          W_CLOUD_IN        , &!REAL        , INTENT(IN ):: W_CLOUD_IN(NPD_PROFILE, NPD_LAYER)!CLOUD AMOUNTS
          TOTAL_CLOUD_COVER , &!REAL        , INTENT(OUT):: TOTAL_CLOUD_COVER(NPD_PROFILE) !TOTAL CLOUD COVER
          NPD_PROFILE       , &!INTEGER , INTENT(IN )::  NPD_PROFILE!MAXIMUM NUMBER OF PROFILES
          NPD_LAYER           )!INTEGER , INTENT(IN )::  NPD_LAYER  !MAXIMUM NUMBER OF LAYERS
!
!
!
      IMPLICIT NONE
!
!
!     DECLARATION OF ARRAY SIZES.
      INTEGER , INTENT(IN )::  NPD_PROFILE!MAXIMUM NUMBER OF PROFILES
      INTEGER , INTENT(IN )::  NPD_LAYER  !MAXIMUM NUMBER OF LAYERS
! CLSCHM3A end
!
!
!     DUMMY ARGUMENTS.
      INTEGER , INTENT(IN ):: N_PROFILE !NUMBER OF PROFILES
      INTEGER , INTENT(IN ):: N_LAYER   !Number of layers seen in radiation
      INTEGER , INTENT(IN ):: NCLDS     !NUMBER OF CLOUDY LAYERS
      INTEGER , INTENT(IN ):: I_CLOUD   !CLOUD SCHEME EMPLOYED
      REAL    , INTENT(IN ):: W_CLOUD_IN(NPD_PROFILE, NPD_LAYER)!CLOUD AMOUNTS
      REAL    , INTENT(OUT):: TOTAL_CLOUD_COVER(NPD_PROFILE) !TOTAL CLOUD COVER
!
!

!
!     COMDECKS INCLUDED
! CLSCHM3A defines reference numbers for cloud schemes in two-stream
! radiation code.

      ! maximum/random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_MAX=2

      ! random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_RANDOM=4

      ! maximum overlap in a column model
      INTEGER,PARAMETER:: IP_CLOUD_COLUMN_MAX=3

      ! clear column
      INTEGER,PARAMETER:: IP_CLOUD_CLEAR=5

      ! mixed column with split between  convective and layer cloud.
      INTEGER,PARAMETER:: IP_CLOUD_TRIPLE=6

      ! Coupled overlap with partial correlation of cloud
      INTEGER,Parameter:: IP_cloud_part_corr=7

      ! Coupled overlap with partial correlation of cloud
      ! with a separate treatment of convective cloud
      INTEGER,Parameter:: IP_cloud_part_corr_cnv=8

!
!
!     LOCAL VARIABLES.
      INTEGER :: L!             LOOP VARIABLE
      INTEGER :: I!             LOOP VARIABLE
      REAL    :: W_CLOUD(NPD_PROFILE, NPD_LAYER)
!
!
!
!     COPY W_CLOUD_IN TO W_CLOUD AND THEN CHECK THAT VALUES ARE
!     BETWEEN 0 AND 1.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
              W_CLOUD(L,I) = MIN( MAX(W_CLOUD_IN(L,I),0.0) , 1.0 )
            ENDDO
         ENDDO
!
!     DIFFERENT OVERLAP ASSUMPTIONS ARE CODED INTO EACH SOLVER.
!
      IF (I_CLOUD == IP_CLOUD_MIX_MAX) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-W_CLOUD(L, N_LAYER+1-NCLDS)
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER-1
            DO L=1, N_PROFILE
               IF (W_CLOUD(L, I+1) >  W_CLOUD(L, I)) THEN
                  TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)             &
                     *(1.0E+00-W_CLOUD(L, I+1))/(1.0E+00-W_CLOUD(L, I))
               ENDIF
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_MIX_RANDOM) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)                &
                  *(1.0E+00-W_CLOUD(L, I))
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=0.0E+00
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               TOTAL_CLOUD_COVER(L)=MAX(TOTAL_CLOUD_COVER(L)            &
                 , W_CLOUD(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_TRIPLE) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-W_CLOUD(L, N_LAYER+1-NCLDS)
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER-1
            DO L=1, N_PROFILE
               IF (W_CLOUD(L, I+1) >  W_CLOUD(L, I)) THEN
                  TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)             &
                     *(1.0E+00-W_CLOUD(L, I+1))/(1.0E+00-W_CLOUD(L, I))
               ENDIF
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_CLEAR) THEN
!
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=0.0E+00
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE R2_CALC_TOTAL_CLOUD_COVER
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete bounds on
!                                               effective radius
!                                               removed.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Subroutine for droplet
!                                               number concentration
!                                               replaced by new function
!                                               with option to use sea-
!                                               salt to supplement
!                                               sulphate aerosol.
!                                               Treatment of convective
!                                               cloud effective radii
!                                               updated.
!                                               (A. Jones)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               NUMBER_DROPLET.
!                                               (A. Jones)
!       6.1             07-04-04                Add new variables for
!                                               column cloud droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Protect calculations
!                                               from failure in PC2
!                                               (Damian Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate column-integrated cloud droplet number.
!
! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.
!
! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       5.4             29-05-02                Original Code
!                                               (A. Jones)
!       6.1             07-04-04                Modified in accordance
!                                               with AVHRR retrievals:
!                                               only clouds >273K used,
!                                               convective clouds also
!                                               included.
!                                               (A. Jones)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             04-03-96                Original Code
!                                               (J. M. Edwards)
!                                               Parts of this code are
!                                               rather redundant. The
!                                               form of writing is for
!                                               near consistency with
!                                               HADAM3.
!
!       4.5   April 1998   Check for inconsistencies between soot
!                          spectral file and options used. L Robinson.
!       5.3     04/04/01   Include mesoscale aerosol switch when
!                          checking if aerosols are required.  S. Cusack
!       5.4     09/05/02   Include logical flag for sea-salt aerosol.
!                                                              A. Jones
!       5.5     05/02/03   Include logical for biomass smoke
!                          aerosol.               P Davison
! Description of Code:
!   5.5    21/02/03 Add logical for d mineral dust
!                                                S Woodward
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !+ Function to calculate cloud droplet number concentration.
  !
  ! Purpose:
  !   Cloud droplet number concentration is calculated from aerosol
  !   concentration or else fixed values are assigned.
  !
  ! Method:
  !   Sulphate aerosol mass concentration is converted to number
  !   concentration by assuming a log-normal size distribution.
  !   Sea-salt and/or biomass-burning aerosols may then
  !   be added if required. The total is then converted to cloud
  !   droplet concentration following the parametrization of
  !   Jones et al. (1994) and lower limits are imposed.
  !   Alternatively, fixed droplet values are assigned if the
  !   parametrization is not required.
  !
  ! Current Owner of Code: A. Jones
  !
  ! History:
  !       Version         Date                    Comment
  !       6.2             20-10-05                Based on NDROP1.
  !                                               Aitken-mode SO4 no
  !                                               longer used & aerosol
  !                                               parameters revised.
  !                                               (A. Jones)
  !
  !
  ! Description of Code:
  !   FORTRAN 77  with extensions listed in documentation.
  !
  !- ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION NUMBER_DROPLET( &
       L_AEROSOL_DROPLET, &
       L_NH42SO4        , &
      !AITKEN_SULPHATE  , &
       ACCUM_SULPHATE   , &
       DISS_SULPHATE    , &
       L_SEASALT_CCN    , &
       SEA_SALT_FILM    , &
       SEA_SALT_JET     , &
       L_BIOGENIC_CCN   , &
       BIOGENIC         , &
       L_BIOMASS_CCN    , &
       BIOMASS_AGED     , &
       BIOMASS_CLOUD    , &
       L_OCFF_CCN       , &
       OCFF_AGED        , &
       OCFF_CLOUD       , &
       DENSITY_AIR      , &
       SNOW_DEPTH       , &
       LAND_FRACT       , &
       NTOT_LAND        , &
       NTOT_SEA           )



    IMPLICIT NONE


    !     Comdecks included:
    !*L------------------COMDECK C_PI---------------------------------------
    !LL
    !LL 4.0 19/09/95  New value for PI. Old value incorrect
    !LL               from 12th decimal place. D. Robinson
    !LL 5.1 7/03/00   Fixed/Free format P.Selwood
    !LL

    ! Pi
    REAL(KIND=r8), PARAMETER :: Pi = 3.14159265358979323846_r8

    ! Conversion factor degrees to radians
    REAL(KIND=r8), PARAMETER :: Pi_Over_180        = Pi/180.0_r8

    ! Conversion factor radians to degrees
    REAL(KIND=r8), PARAMETER :: Recip_Pi_Over_180  = 180.0_r8/Pi

    !*----------------------------------------------------------------------


    LOGICAL, INTENT(IN   ) ::  L_AEROSOL_DROPLET !             Flag to use aerosols to find droplet number
    LOGICAL, INTENT(IN   ) ::  L_NH42SO4         !             Is the input "sulphate" aerosol in the form of
    !             ammonium sulphate (T) or just sulphur (F)?
    LOGICAL, INTENT(IN   ) ::  L_SEASALT_CCN     !             Is sea-salt aerosol to be used?
    LOGICAL, INTENT(IN   ) ::  L_BIOMASS_CCN     !             Is biomass smoke aerosol to be used?
    LOGICAL, INTENT(IN   ) ::  L_BIOGENIC_CCN    !             Is biogenic aerosol to be used?
    LOGICAL, INTENT(IN   ) ::  L_OCFF_CCN        !             Is fossil-fuel organic carbon aerosol to be used?

    !      REAL    ::  AITKEN_SULPHATE                  !not use      Dummy in version 2 of the routine
    REAL(KIND=r8)   , INTENT(IN   ) ::  ACCUM_SULPHATE    !             Mixing ratio of accumulation-mode sulphate aerosol
    REAL(KIND=r8)   , INTENT(IN   ) ::  DISS_SULPHATE     !             Mixing ratio of dissolved sulphate aerosol
    REAL(KIND=r8)   , INTENT(IN   ) ::  BIOMASS_AGED      !             Mixing ratio of aged biomass smoke
    REAL(KIND=r8)   , INTENT(IN   ) ::  BIOMASS_CLOUD     !             Mixing ratio of in-cloud biomass smoke
    REAL(KIND=r8)   , INTENT(IN   ) ::  SEA_SALT_FILM     !             Number concentration of film-mode sea salt aerosol (m-3)
    REAL(KIND=r8)   , INTENT(IN   ) ::  SEA_SALT_JET      !             Number concentration of jet-mode sea salt aerosol (m-3)
    REAL(KIND=r8)   , INTENT(IN   ) ::  BIOGENIC          !             Mixing ratio of biogenic aerosol
    REAL(KIND=r8)   , INTENT(IN   ) ::  OCFF_AGED         !             Mixing ratio of aged fossil-fuel organic carbon
    REAL(KIND=r8)   , INTENT(IN   ) ::  OCFF_CLOUD        !             Mixing ratio of in-cloud fossil-fuel organic carbon
    REAL(KIND=r8)   , INTENT(IN   ) ::  DENSITY_AIR       !             Density of air (kg m-3)
    REAL(KIND=r8)   , INTENT(IN   ) ::  SNOW_DEPTH        !             Snow depth (m; >5000 is flag for ice-sheets)
    REAL(KIND=r8)   , INTENT(IN   ) ::  LAND_FRACT        !             Land fraction
    REAL(KIND=r8)   , INTENT(IN   ) ::  NTOT_LAND         !             Droplet number over land if parameterization is off (m-3)
    REAL(KIND=r8)   , INTENT(IN   ) ::  NTOT_SEA          !             Droplet number over sea if parameterization is off (m-3)

    !REAL(KIND=r8)    :: NUMBER_DROPLET     !             Returned number concentration of cloud droplets (m-3)



    !     Local variables:

    REAL(KIND=r8)    ::  PARTICLE_VOLUME   !             Mean volume of aerosol particle
    REAL(KIND=r8)    ::  N_CCN             !             Number density of CCN

    REAL(KIND=r8) ,PARAMETER :: RADIUS_0_NH42SO4=9.5E-8_r8  ! Median radius of log-normal distribution for (NH4)2SO4
    REAL(KIND=r8) ,PARAMETER :: SIGMA_0_NH42SO4=1.4_r8      ! Geometric standard deviation of same
    REAL(KIND=r8) ,PARAMETER :: DENSITY_NH42SO4=1.769E+03_r8! Density of ammonium sulphate aerosol
    REAL(KIND=r8) ,PARAMETER :: RADIUS_0_BIOMASS=1.2E-07_r8 ! Median radius of log-normal distribution for biomass smoke
    REAL(KIND=r8) ,PARAMETER :: SIGMA_0_BIOMASS=1.30_r8     ! Geometric standard deviation of same
    REAL(KIND=r8) ,PARAMETER :: DENSITY_BIOMASS=1.35E+03_r8 ! Density of biomass smoke aerosol
    REAL(KIND=r8) ,PARAMETER :: RADIUS_0_BIOGENIC=9.5E-08_r8! Median radius of log-normal dist. for biogenic aerosol
    REAL(KIND=r8) ,PARAMETER :: SIGMA_0_BIOGENIC=1.50_r8    ! Geometric standard deviation of same
    REAL(KIND=r8) ,PARAMETER :: DENSITY_BIOGENIC=1.3E+03_r8 ! Density of biogenic aerosol
    REAL(KIND=r8) ,PARAMETER :: RADIUS_0_OCFF=0.12E-06_r8   ! Median radius of log-normal dist. for OCFF aerosol
    REAL(KIND=r8) ,PARAMETER :: SIGMA_0_OCFF=1.30_r8        ! Geometric standard deviation of same
    REAL(KIND=r8) ,PARAMETER :: DENSITY_OCFF=1350.0_r8      ! Density of OCFF aerosol




    IF (L_AEROSOL_DROPLET) THEN

       !        If active, aerosol concentrations are used to calculate the
       !        number of CCN, which is then used to determine the number
       !        concentration of cloud droplets (m-3).

       PARTICLE_VOLUME=(4.0E+00_r8*PI/3.0E+00_r8)*RADIUS_0_NH42SO4**3       &
            &                         *EXP(4.5E+00_r8*(LOG(SIGMA_0_NH42SO4))**2)

       IF (L_NH42SO4) THEN
          !           Input data have already been converted to ammonium sulphate.
          N_CCN=(ACCUM_SULPHATE+DISS_SULPHATE)                        &
               &                  *DENSITY_AIR/(DENSITY_NH42SO4*PARTICLE_VOLUME)
       ELSE
          !           Convert m.m.r. of sulphur to ammonium sulphate by
          !           multiplying by ratio of molecular weights:
          N_CCN=(ACCUM_SULPHATE+DISS_SULPHATE)*4.125_r8                  &
               &                  *DENSITY_AIR/(DENSITY_NH42SO4*PARTICLE_VOLUME)
       ENDIF

       IF (L_SEASALT_CCN) THEN
          N_CCN=N_CCN+SEA_SALT_FILM+SEA_SALT_JET
       ENDIF

       IF (L_BIOMASS_CCN) THEN
          PARTICLE_VOLUME=(4.0E+00_r8*PI/3.0E+00_r8)*RADIUS_0_BIOMASS**3    &
               &                         *EXP(4.5E+00_r8*(LOG(SIGMA_0_BIOMASS))**2)
          N_CCN=N_CCN+((BIOMASS_AGED+BIOMASS_CLOUD)*DENSITY_AIR       &
               &                             /(DENSITY_BIOMASS*PARTICLE_VOLUME))
       ENDIF

       IF (L_BIOGENIC_CCN) THEN
          PARTICLE_VOLUME=(4.0E+00_r8*PI/3.0E+00_r8)*RADIUS_0_BIOGENIC**3   &
               &                         *EXP(4.5E+00_r8*(LOG(SIGMA_0_BIOGENIC))**2)
          N_CCN=N_CCN+(BIOGENIC*DENSITY_AIR                           &
               &                             /(DENSITY_BIOGENIC*PARTICLE_VOLUME))
       ENDIF

       IF (L_OCFF_CCN) THEN
          PARTICLE_VOLUME=(4.0E+00_r8*PI/3.0E+00_r8)*RADIUS_0_OCFF**3       &
               &                         *EXP(4.5E+00_r8*(LOG(SIGMA_0_OCFF))**2)
          N_CCN=N_CCN+((OCFF_AGED+OCFF_CLOUD)*DENSITY_AIR             &
               &                             /(DENSITY_OCFF*PARTICLE_VOLUME))
       ENDIF

       !        Apply relation of Jones et al. (1994) to get droplet number
       !        and apply minimum value (equivalent to 5 cm-3):

       NUMBER_DROPLET=3.75E+08_r8*(1.0E+00_r8-EXP(-2.5E-9_r8*N_CCN))

       IF (NUMBER_DROPLET  <   5.0E+06_r8) THEN
          NUMBER_DROPLET=5.0E+06_r8
       ENDIF

       !        If gridbox is more than 20% land AND this land is not covered
       !        by an ice-sheet, use larger minimum droplet number (=35 cm-3):

       IF (LAND_FRACT  >   0.2_r8 .AND. SNOW_DEPTH  <   5000.0_r8           &
            &                        .AND. NUMBER_DROPLET  <   35.0E+06_r8) THEN
          NUMBER_DROPLET=35.0E+06_r8
       ENDIF

    ELSE

       !        Without aerosols, the number of droplets is fixed; a simple
       !        50% criterion is used for land or sea in this case.

       IF (LAND_FRACT  >=  0.5_r8) THEN
          NUMBER_DROPLET=NTOT_LAND
       ELSE
          NUMBER_DROPLET=NTOT_SEA
       ENDIF
       !
    ENDIF
    !
    !
    !
    RETURN
  END FUNCTION NUMBER_DROPLET



END Module CalCloudCover
PROGRAM Main
  USE CalCloudCover
END PROGRAM Main
