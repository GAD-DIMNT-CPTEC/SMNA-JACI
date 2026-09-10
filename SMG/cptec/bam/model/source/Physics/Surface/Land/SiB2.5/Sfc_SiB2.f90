
! Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
!  Fine aerosol mode is included
!  Modifications (5) are marked by 
!  !tar begin   and  !tar end

MODULE SFC_SiB2
  ! InitSSiB2
  !3947-4548
  !sib --------|-- rnload 
  !            | 
  !            |-- balan 
  !            | 
  !            |-- begtem 
  !            | 
  !            |-- netrad 
  !            | 
  !            |-- vnqsat 
  !            | 
  !            |-- vntlat --- vmfcalz 
  !            |           |
  !            |           |- vmfcal
  !            |           | 
  !            |           |- rbrd
  !            |           | 
  !            |           |- respsib
  !            |           |
  !            |           |- phosib --- sortin
  !            |                      |
  !            |                      |--cycalc
  !            |
  !            | 
  !            |-- vmfcalzo
  !            | 
  !            |-- vmfcalo
  !            | 
  !            |-- dellwf
  !            |
  !            |-- delhf
  !            |
  !            |-- delef
  !            |
  !            |-- soilprop
  !            |
  !            |-- sibslv --- gauss
  !            |
  !            |-- updat2
  !            !
  !            |-- soiltherm
  !            |
  !            |-- addinc 
  !            |
  !            |-- inter2 -- adjust
  USE Constants, ONLY :     &
       r8,r4,i4,i8,pie,gasr,ityp, imon, icg, iwv, idp, ibd,oceald,icealv,icealn

  USE Options, ONLY: &
       reducedGrid,nfsibt,nfsoiltb,nfmorftb,nfbioctb,nfbioctb,nfsib2mk,nfsand,&
       nfclay,nftext,nfsibd,record_type,isimp,nftgz0,nfzol,initlz,ifalb,ifsst,ifco2flx,ifndvi,ifslm,&
       ifslmsib2,ifsnw,ifozone,sstlag,intsst,intndvi,intsoilm,fint,iglsm_w,mxiter,nfsibi,ifozone,iftracer,&
!tar begin
!climate aerosol selection parameter
        ifaeros, &
!tar end
       nfsoiltp,nfvegtp,nfslmtp,nfsoiltp,OCFLUX,omlmodel,oml_hml0,&
       nfprt, nfctrl,  nfalb,filta,epsflt,istrt,yrl   ,monl,sfcpbl,atmpbl,&
       fNameSoilType,fNameVegType,fNameSoilMoist,fNameRouLen, &
       fNameSibmsk,fNameTg3zrl  ,fNameSoilTab,fNameMorfTab,fNameBioCTab,&
       fNameAeroTab,fNameSiB2Mask,fNameSandMask,fNameClayMask,fNameTextMask,SLABOCEAN,ICEMODEL
  USE FieldsPhysics, ONLY: &
      npatches     , &
      npatches_actual, &
      nzg           , &
      sheleg       , &
      imask           , &
      SoilMask        , &
      AlbVisDiff   , &
      gtsea           , &
      gco2flx,       &
      gndvi           , &
      soilm           , &
      o3mix           , &
!tar begin
!climate aerosol optical parameters of coarse mode
      aod             , &
      asy             , &
      ssa             , &
      z_aer           , &  
!tar end
!
!tar begin
!climate aerosol optical parameters of fine mode
      aodF             , &
      asyF             , &
      ssaF             , &
      z_aerF           , &  
!tar end
      tg1           , &
      tg2           , &
      tg3           , &
      rVisDiff     , &
      ssib           , &
      wsib3d       , &
      ppli           , &
      ppci           , &
      capac0       , &
      gl0           , &
      Mmlen           , &
      tseam           , &
      ndvim           , &
      w0           , &
      tg0           , &
      tc0           , &
      tm0          , &
      qm0          , &
      tmm          , &
      qmm          , &
      tcm          , &  
      tgm          , &  
      wm           , &   
      capacm       , &
      capacc       , &
      qsfc0           , &
      tsfc0           , &
      qsfcm           , &
      tsfcm           , &
      tkemyj       , &
      z0           , &
      zorl           , &
      MskAnt        , &
      tracermix     , &
      HML            , &
      HUML           , &
      HVML           , &
      TSK            , &
      z0sea          , &
      mlsi           , &  ! add solange 13-11-2012
      sm0            , &  ! add solange 13-11-2012
      laymld,       hbath,     tdeep,sdeep, &
      sfc,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 

  USE Utils, ONLY: &
       IJtoIBJB, &
       LinearIJtoIBJB, &
       NearestIJtoIBJB, &
       !SeaMaskIJtoIBJB, &
       !SplineIJtoIBJB, &
       !AveBoxIJtoIBJB, &
       FreqBoxIJtoIBJB, &
       vfirec

  USE InputOutput, ONLY: &
       getsbc

  USE IOLowLevel, ONLY: &
       ReadGetNFTGZ
 
  USE Parallelism, ONLY: &
       MsgOne, FatalError

  USE Sfc_SeaFlux_Interface   , Only :  seasfc

  USE SlabOceanModel  , Only : GetOceanAlb

  USE Sfc_SeaIceFlux_WRF_Model  , Only : GetIceOceanAlb,&
      TC_SeaIce     ,&
      TGS_SeaIce ,&
      TD_SeaIce  ,&
      TA_SeaIce  ,&
      SNOA_SeaIce,&
      SNOB_SeaIce

  USE PhysicalFunctions, ONLY : fpvs,fpvs2es5

  IMPLICIT NONE
SAVE

  PRIVATE
  INTEGER      ,PUBLIC  :: nsoil         !  number of soil levels
  REAL(KIND=r8)           :: ztemp             !  height of temperature measurement (m)
  REAL(KIND=r8)           :: zwind             !  height of wind measurement (m)
  INTEGER :: idatec(4)
  REAL (KIND=r8)  , PARAMETER   :: z0ice  = 0.001e0_r8! 
  REAL(KIND=r8)   , PARAMETER   :: rgas   = 287.0_r8!  dry air gas constant (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: kapa   = 0.2861328125_r8!  rgas/cp (unitless)
  REAL(KIND=r8)   , PARAMETER   :: pi     = 3.1415926_r8 !  pi 
  REAL(KIND=r8)   , PARAMETER   :: grav   = 9.81_r8                !  gravitational constant (m sec^-2)
  REAL(KIND=r8)   , PARAMETER   :: cp     = rgas/kapa!  specific of heat of dry air at constant pressure (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: cv     = 1952.0_r8!  specific heat of water vapor at constant pressure (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: hltm   = 2.52e6_r8            !  latent heat of vaporization (J kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: delta  = 0.608_r8               !  UNKNOWN  
  REAL(KIND=r8)   , PARAMETER   :: asnow  = 16.7_r8                !  UNKNOWN
  REAL(KIND=r8)   , PARAMETER   :: snomel  = 333624.2e0_r8 * 1000.0e0       !  latent heat of fusion of ice (J m^-3) 
  REAL(KIND=r8)   , PARAMETER   :: clai   = 4.2e0_r8*1000.0_r8*0.2_r8    !  leaf heat capacity  (J m^-3 deg^-1)
  REAL(KIND=r8)   , PARAMETER   :: cww    = 4.2e0_r8*1000.0_r8*1000.0_r8 !  water heat capacity (J m^-3 deg^-1)
  REAL(KIND=r8)   , PARAMETER   :: pr0    = 0.74_r8                !  turb Prandtl Number at neutral stblty
  REAL(KIND=r8)   , PARAMETER   :: ribc   = 3.05_r8                !  critical Richardson Number (unitless)
  REAL(KIND=r8)   , PARAMETER   :: vkrmn  = 0.35_r8                !  Von Karmann's constant (unitless)
  REAL(KIND=r8)   , PARAMETER   :: po2m   = 20900.0_r8              !  mixed layer O2 concentration
  REAL(KIND=r8)   , PARAMETER   :: stefan = 5.67e-8_r8           !  Stefan-Boltzmann constant
  REAL(KIND=r8)   , PARAMETER   :: grav2  = grav *0.01_r8          !  grav/100 (Pa/hPa)
  REAL(KIND=r8)   , PARAMETER   :: tice   = 273.1_r8               !  freezing temperature (K)
 
  LOGICAL, PARAMETER   :: forcerestore  = .FALSE.    !
  INTEGER, PARAMETER   :: SchemeDifus   = 1 !=1 Original Sib2
                                            !=2 Original SSib
  
  LOGICAL, PARAMETER   :: sibdrv   = .TRUE.      !  kept because it is in the SiB code...
  LOGICAL, PARAMETER   :: dosibco2 = .TRUE.     !  flag-calculate CO2 flux?
  LOGICAL              :: fixday   = .FALSE.    !  perpetual day of year
  LOGICAL              :: dotkef   = .FALSE.     !  use changan instead of deardorff flux
  INTEGER, PARAMETER   :: louis    =   3          !  use Louis (1850) vmf calculations
  !INTEGER, PARAMETER   :: louis =  2          !  use  VENTILATION MASS FLUX, BASED ON DEARDORFF, MWR, 1972
  !INTEGER, PARAMETER   :: louis =  3          !  use Deardorff, mwr, 1972? BASED ON SSiB.
   REAL(KIND=r8), PARAMETER :: asat0 =  6.1078000_r8
  REAL(KIND=r8), PARAMETER :: asat1 =  4.4365185e-1_r8
  REAL(KIND=r8), PARAMETER :: asat2 =  1.4289458e-2_r8
  REAL(KIND=r8), PARAMETER :: asat3 =  2.6506485e-4_r8
  REAL(KIND=r8), PARAMETER :: asat4 =  3.0312404e-6_r8
  REAL(KIND=r8), PARAMETER :: asat5 =  2.0340809e-8_r8
  REAL(KIND=r8), PARAMETER :: asat6 =  6.1368209e-11_r8
  !
  !
  REAL(KIND=r8), PARAMETER :: bsat0 =  6.1091780_r8
  REAL(KIND=r8), PARAMETER :: bsat1 =  5.0346990e-1_r8
  REAL(KIND=r8), PARAMETER :: bsat2 =  1.8860134e-2_r8
  REAL(KIND=r8), PARAMETER :: bsat3 =  4.1762237e-4_r8
  REAL(KIND=r8), PARAMETER :: bsat4 =  5.8247203e-6_r8
  REAL(KIND=r8), PARAMETER :: bsat5 =  4.8388032e-8_r8
  REAL(KIND=r8), PARAMETER :: bsat6 =  1.8388269e-10_r8
  !
  REAL(KIND=r8), PARAMETER :: csat0 =  4.4381000e-1_r8
  REAL(KIND=r8), PARAMETER :: csat1 =  2.8570026e-2_r8
  REAL(KIND=r8), PARAMETER :: csat2 =  7.9380540e-4_r8
  REAL(KIND=r8), PARAMETER :: csat3 =  1.2152151e-5_r8
  REAL(KIND=r8), PARAMETER :: csat4 =  1.0365614e-7_r8
  REAL(KIND=r8), PARAMETER :: csat5 =  3.5324218e-10_r8
  REAL(KIND=r8), PARAMETER :: csat6 = -7.0902448e-13_r8
  !
  !
  REAL(KIND=r8), PARAMETER :: dsat0 =  5.0303052e-1_r8
  REAL(KIND=r8), PARAMETER :: dsat1 =  3.7732550e-2_r8
  REAL(KIND=r8), PARAMETER :: dsat2 =  1.2679954e-3_r8
  REAL(KIND=r8), PARAMETER :: dsat3 =  2.4775631e-5_r8
  REAL(KIND=r8), PARAMETER :: dsat4 =  3.0056931e-7_r8
  REAL(KIND=r8), PARAMETER :: dsat5 =  2.1585425e-9_r8
  REAL(KIND=r8), PARAMETER :: dsat6 =  7.1310977e-12_r8
 
  LOGICAL              :: isotope = .FALSE.      ! invoke isotope calculations?
  REAL(KIND=r8), PARAMETER :: Epsilonk   = 0.62198_r8
  REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
  REAL(KIND=r8), PARAMETER :: one_minus_epsilon = 1.0_r8 -epsilonk

  !      One minus the ratio of the molecular weights of water and dry

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
       &  N = INT(((T_high - T_low + (delta_T*0.5_r8))/delta_T) + 1.0_r8,kind=i4)
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

  !  begin input data variables
  CHARACTER(LEN=70) :: CenturyLai      ! century lai input file
  REAL(KIND=r8)     :: NDVIoffset
  REAL(KIND=r8)     :: NDVIscale
  REAL(KIND=r8)     :: LatMin
  REAL(KIND=r8)     :: LonMin
  INTEGER           :: use100     ! switch for using Century LAI (1=yes)
  REAL(KIND=r8)     :: Dlat
  REAL(KIND=r8)     :: Dlon
  INTEGER           :: imm  
  INTEGER           :: jmm  
  INTEGER SoilNum   ! soil type number
  !
  ! begin soil dependant variables
  !
  !
  ! begin biome dependant, physical morphology variables
  REAL(KIND=r8), PARAMETER :: fPARmax=0.95_r8  ! Max possible FPAR for 98th percentile
  REAL(KIND=r8), PARAMETER :: fPARmin=0.01_r8  ! Min possible FPAR for 2nd percentile

  !
  ! begin aerodynamic interpolation tables
  !
  REAL (KIND=r8) :: LAIgrid    (50)   ! grid of LAI values for lookup table
  REAL (KIND=r8) :: fVCovergrid(50)! grid of fVCover values for interpolation table
  !

  !
  ! begin Input program control variables
  INTEGER :: nm         ! total number of months in a year
  INTEGER :: ndays      ! total number of days in a year
  INTEGER :: nv         ! total number of possible vegetation types
  INTEGER :: ns         ! total number of possible soil types
  INTEGER :: minBiome   ! minimum biome number for biome subset option
  INTEGER :: maxBiome   ! maximum biome number for biome subset option


  !
  ! Begin Input program option flags
  TYPE Flags
     INTEGER Map      ! execute mapper subroutine
     INTEGER Mode     ! 'grid' or 'single' point mode
     INTEGER Sib      ! generate SiB BC input file
     INTEGER single   ! generate SiB BC input file ascii format single pt
     INTEGER EzPlot   ! generate generic ezplot output
     INTEGER Stats    ! generate output for statistics
     INTEGER monthly  ! generate monthly output files
     INTEGER GridKey  ! generate ascci file of lat/lon for each SiB point
     INTEGER PRINT    ! print statements to screen
     INTEGER SoRefTab ! use soil relectance look up table
     INTEGER SoilMap  ! use %clay/sand or soil type input maps
     INTEGER SoilProp ! soil properties from look up table or % sand/clay
     INTEGER fVCov    ! calculate fVCover or use input map
  END TYPE Flags
  TYPE(Flags) Flagsib

  !
  !
  ! Begin filenames for input lookup tables and maps
  TYPE FileNames
     CHARACTER(LEN=255) :: BioMap   ! vegetation type map
     CHARACTER(LEN=255) :: BioTab   ! veg. type char. lookup table
     CHARACTER(LEN=255) :: MorphTab ! veg. morphilogical lookup table
     CHARACTER(LEN=255) :: SoilTab  ! soil type lookup table
     CHARACTER(LEN=255) :: SoRefVis ! soil ref. map visible
     CHARACTER(LEN=255) :: SoRefNIR ! soil ref map near IR
     CHARACTER(LEN=255) :: SoilMap  ! soil characteristic map
     CHARACTER(LEN=255) :: AeroVar  ! aerodynamic interpolation tables
     CHARACTER(LEN=255) :: PercClay ! percentage soil that is clay map
     CHARACTER(LEN=255) :: PercSand ! percentage soil that is sand map
     CHARACTER(LEN=255) :: fVCovMap ! fraction veg. cover map
     CHARACTER(LEN=255) :: sib_grid ! SiB gridmap file
     CHARACTER(LEN=255) :: sib_bc   ! SiB boundary condition file
     CHARACTER(LEN=255) :: sib_biom ! GCM biome file
     CHARACTER(LEN=255) :: stats    ! generic statistics file
     CHARACTER(LEN=255) :: Monthly  ! monthly files
     CHARACTER(LEN=255) :: GridKey  ! lat/lon grid key for SiB Points
     CHARACTER(LEN=255) :: asciiNDVI! ascii ndvi file for single point option
  END TYPE FileNames
  TYPE(FileNames) FileName    ! file of filenames for input tables and maps

  TYPE biome_morph_var
     REAL(KIND=r8) :: zc         ! Canopy inflection height (m)
     REAL(KIND=r8) :: LWidth         ! Leaf width
     REAL(KIND=r8) :: LLength         ! Leaf length 
     REAL(KIND=r8) :: LAImax         ! Maximum LAI
     REAL(KIND=r8) :: stems         ! Stem area index
     REAL(KIND=r8) :: NDVImax         ! Maximum NDVI
     REAL(KIND=r8) :: NDVImin         ! Minimum NDVI
     REAL(KIND=r8) :: SRmax         ! Maximum simple ratio
     REAL(KIND=r8) :: SRmin         ! Minimum simple ratio
  END TYPE biome_morph_var
  TYPE(biome_morph_var), ALLOCATABLE :: MorphTab(:)! Lookup table of biome dependant morphology

  TYPE soil_Physical
     INTEGER(KIND=i8):: SoilNum ! soil type number
     REAL(KIND=r8)   :: BEE     ! Soil wetness exponent
     REAL(KIND=r8)   :: PhiSat  ! Soil tension at saturation
     REAL(KIND=r8)   :: SatCo   ! Hydraulic conductivity at saturation
     REAL(KIND=r8)   :: poros   ! Soil porosity
     REAL(KIND=r8)   :: Slope   ! Cosine of mean slope
     REAL(KIND=r8)   :: Wopt    ! optimal soil wetness for soil respiration
     REAL(KIND=r8)   :: Skew    ! skewness exponent of soil respiration vs. wetness curve
     REAL(KIND=r8)   :: RespSat ! assures soil respiration is 60-80% of max @ saturation
  END TYPE soil_Physical
  TYPE(soil_Physical) SoilVar  ! time ind., soil dependant variables
  TYPE(soil_Physical), ALLOCATABLE :: SoilTab(:)
  TYPE(soil_Physical), ALLOCATABLE :: bSoilGrd(:,:)


  TYPE aero_var
     REAL(KIND=r8) :: zo                ! Canopy roughness coeff 
     REAL(KIND=r8) :: zp_disp    ! Zero plane displacement
     REAL(KIND=r8) :: RbC        ! RB Coefficient
     REAL(KIND=r8) :: RdC        ! RC Coefficient
  END TYPE aero_var
  TYPE(aero_var), ALLOCATABLE :: AeroVar(:,:,:) ! interpolation tables for aero variables

  !
  !
  ! begin Biome-dependent variables
  ! same structure applies to input lookup tables and output files
  !
  TYPE Biome_dep_var
     INTEGER bioNum  ! biome or vegetation cover type
     REAL(KIND=r8) :: z2         ! Canopy top height (m)
     REAL(KIND=r8) :: z1         ! Canopy base height (m)
     REAL(KIND=r8) :: fVCover    ! Canopy cover fraction
     REAL(KIND=r8) :: ChiL       ! Leaf angle distribution factor
     REAL(KIND=r8) :: SoDep      ! Total depth of 3 soil layers (m)
     REAL(KIND=r8) :: RootD      ! Rooting depth (m)
     REAL(KIND=r8) :: Phi_half   ! 1/2 Critical leaf water potential limit (m)
     REAL(KIND=r8) :: LTran(2,2) ! Leaf transmittance for green/brown plants
     REAL(KIND=r8) :: LRef(2,2)  ! Leaf reflectance for green/brown plants
     ! For LTran and LRef:
     ! (1,1)=shortwave, green plants
     ! (2,1)=longwave, green plants
     ! (1,2)=shortwave, brown plants
     ! (2,2)=longwave, brown plants
     REAL(KIND=r8) :: vmax0      ! Rubisco velocity of sun-leaf (Mol m^-2 s^-1)
     REAL(KIND=r8) :: EffCon     ! Quantum efficiency (Mol Mol^-1)
     REAL(KIND=r8) :: gsSlope    ! Conductance-Photosynthesis Slope Parameter
     REAL(KIND=r8) :: gsMin      ! Conductance-Photosynthesis Intercept
     REAL(KIND=r8) :: Atheta     ! WC WE Coupling Parameter
     REAL(KIND=r8) :: Btheta     ! WC & WE, WS Coupling Parameter
     REAL(KIND=r8) :: TRDA       ! Temperature Coefficient in GS-A Model (K^-1)
     REAL(KIND=r8) :: TRDM       ! "" (K)
     REAL(KIND=r8) :: TROP       ! "" (K)
     REAL(KIND=r8) :: respcp     ! Respiration Fraction of Vmax
     REAL(KIND=r8) :: SLTI       ! Slope of low-temp inhibition function (K^-1)
     REAL(KIND=r8) :: HLTI       ! Slope of high-temp inhibition function (K^-1)
     REAL(KIND=r8) :: SHTI       ! 1/2 Point of low-temp inhibition function (K)
     REAL(KIND=r8) :: HHTI       ! 1/2 Point of high-temp inhibition function (K)
     REAL(KIND=r8) :: SoRef(2)   ! 2-stream soil and litter reflectivity
     ! Soref(1)=visible soil and litter reflectivity
     ! Soref(2)=Near IR soil and litter reflectivity
  END TYPE Biome_dep_var
  !
  TYPE(Biome_dep_var) BioVar ! time ind., biome dependant variables
  TYPE(Biome_dep_var), ALLOCATABLE :: BioTab(:)! Lookup table of biome dependant variables

  !
  ! begin Soil texture variables (Only clay and sand need to be specified)
  !
  TYPE soil_texture
     REAL(KIND=r8)    :: clay       ! Percent clay content
     REAL(KIND=r8)    :: silt       ! Percent silt content
     REAL(KIND=r8)    :: sand       ! Percent sand content
     INTEGER :: class      ! Soil texture class
  END TYPE soil_texture
  TYPE(soil_texture) text ! soil texture: percent sand, silt, and clay

  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE :: iMaskSiB2 (:,:)   !  vegetation mask

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: NDVI(:,:,:)! array ofFASIR NDVI values for single pt
  REAL(KIND=r8), TARGET   , ALLOCATABLE :: vcoverg (:,:)           !  vegetation cover fraction

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_zc_sib2     (:,:)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_lwidth_sib2 (:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_llength_sib2(:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_laimax_sib2 (:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_stems_sib2  (:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_ndvimax_sib2(:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_ndvimin_sib2(:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_srmax_sib2  (:,:) 
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: morphtab_srmin_sib2  (:,:) 

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: aerovar_zo_sib2      (:,:,:,:)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: aerovar_zp_disp_sib2 (:,:,:,:)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: aerovar_rbc_sib2     (:,:,:,:)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: aerovar_rdc_sib2     (:,:,:,:)

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: chil_sib2 (:,:)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: tran_sib2 (:,:,:,:)
  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: ref_sib2  (:,:,:,:)

  REAL(KIND=r8), PUBLIC   , ALLOCATABLE :: vcover2g (:,:)           !  vegetation cover fraction

  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_fpar      (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_lai       (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_green     (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_zo             (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_zp_disp   (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_rbc      (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_rdc      (:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: timevar_gmudmu    (:,:)

  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: sandfrac(:,:)           !  fraction of sand in the soil
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: clayfrac(:,:)           !  fraction of clay in the soil

  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: respfactor(:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: pco2m2(:,:)         !  mixed layer CO2 partial pressure (Pa)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: pco2ap(:,:)        !  canopy air space CO2 partial pressure (Pa)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: d13cca(:,:)           !  del 13C of canopy air space
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tmgc       (:,:)     !  ground temperature (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: qmgc       (:,:)     !  ground temperature (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tcasc     (:,:)     !  canopy air space (CAS) temp (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: qcasc     (:,:)     !  (CAS) moisture mixing ratio (kg/kg)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tcalc     (:,:)     !  canopy leaves temp (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tgrdc     (:,:)     !  ground temperature (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tdgc    (:,:,:)   !  soil temp (nsib,nsoil) (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: stores  (:,:)     !  stomatal resistance (1/gs) (m/sec)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: snowg   (:,:,:)    ! snow depth  (m)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: wwwgc     (:,:,:)    ! soil moisture (% of saturation)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: ventmf2   (:,:)      !  output-ventilation mass flux
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: thvgm2   (:,:)            !  output-sfc air moisture deficit (theta)  

  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: c4fractg  (:,:)  ! fraction of C4 vegetation  
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: d13crespg (:,:)  ! del 13C of respiration (per mil vs PDB)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: d13ccag   (:,:)  ! del 13C of canopy CO2 (per mil vs PDB)

  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: sensf     (:,:)    ! soil moisture (% of saturation)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: latef     (:,:)    ! soil moisture (% of saturation)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: co2fx     (:,:)   
  INTEGER(KIND=r8), PUBLIC, ALLOCATABLE :: MskAntSib2(:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: radfac_gbl(:,:,:,:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: AlbGblSiB2(:,:,:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: thermk_gbl(:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tgeff4_sib_gbl(:,:)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: glsm_w      (:,:,:     ) 
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: veg_type (:,:,:)! SIB veg type
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: zlwup_SiB2(:,:) 
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tcasm     (:,:)     !  canopy air space (CAS) temp (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tcas0     (:,:)     !  canopy air space (CAS) temp (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: qcasm     (:,:)     !  (CAS) moisture mixing ratio (kg/kg)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: qcas0     (:,:)     !  (CAS) moisture mixing ratio (kg/kg)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: frac_occ(:,:,:) ! fractional area
  ! coverage
  REAL(KIND=r8), PUBLIC, ALLOCATABLE :: soil_type(:,:  )! FAO/USDA soil texture
  REAL(KIND=r8), PUBLIC, ALLOCATABLE :: tdgm      (:,:,:)   !  soil temp (nsib,nsoil) (K)
  REAL(KIND=r8),PUBLIC   , ALLOCATABLE :: tdg0      (:,:,:)   !  soil temp (nsib,nsoil) (K)

  REAL (KIND=r4), ALLOCATABLE :: rstpar_r4(:,:,:)
  REAL (KIND=r4), ALLOCATABLE :: chil_r4  (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: topt_r4  (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: tll_r4   (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: tu_r4    (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: defac_r4 (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: ph1_r4   (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: ph2_r4   (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: rootd_r4 (:,:) 
  REAL (KIND=r4), ALLOCATABLE :: bee_r4   (:)         
  REAL (KIND=r4), ALLOCATABLE :: phsat_r4 (:)         
  REAL (KIND=r4), ALLOCATABLE :: satco_r4 (:)         
  REAL (KIND=r4), ALLOCATABLE :: poros_r4 (:)         
  REAL (KIND=r4), ALLOCATABLE :: zdepth_r4(:,:)        
  REAL (KIND=r4), ALLOCATABLE :: green_r4 (:,:,:)
  REAL (KIND=r4), ALLOCATABLE :: xcover_r4(:,:,:)
  REAL (KIND=r4), ALLOCATABLE :: zlt_r4   (:,:,:)
  REAL (KIND=r4), ALLOCATABLE :: x0x_r4   (:,:)  
  REAL (KIND=r4), ALLOCATABLE :: xd_r4    (:,:)  
  REAL (KIND=r4), ALLOCATABLE :: z2_r4    (:,:)  
  REAL (KIND=r4), ALLOCATABLE :: z1_r4    (:,:)  
  REAL (KIND=r4), ALLOCATABLE :: xdc_r4   (:,:)  
  REAL (KIND=r4), ALLOCATABLE :: xbc_r4   (:,:)  

  CHARACTER(LEN=255) :: path
  CHARACTER(LEN=255) :: fNameSibVeg
  CHARACTER(LEN=255) :: fNameSibAlb
  real, parameter, public :: MAPL_AIRMW  = 28.97                  ! kg/Kmole
  real, parameter, public :: MAPL_H2OMW  = 18.01                  ! kg/Kmole

  real, parameter, public :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0   ! --

  PUBLIC :: InitSfcSib2
  PUBLIC :: SiB2_Driver
  PUBLIC :: Albedo_sib2
  PUBLIC :: Phenology_sib2
  PUBLIC :: InitCheckSiB2File
  PUBLIC :: InitSurfTempSiB2
  PUBLIC :: ReStartSiB2
  PUBLIC :: Finalize_SiB2
CONTAINS
  SUBROUTINE InitSfcSib2(ibMax         ,jbMax         ,iMax      ,jMax   , &
                         kMax          ,path_in       ,fNameSibVeg_in   , &
                         fNameSibAlb_in,ifdy              ,ids         ,                &
                         idc           ,ifday         ,tod       ,todsib        , &
                         idate         ,idatec        , &
                         ibMaxPerJB  )
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: iMax
    INTEGER, INTENT(IN   ) :: jMax
    INTEGER, INTENT(IN   ) :: kMax
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibVeg_in
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibAlb_in
        
    INTEGER         , INTENT(OUT  ) :: ifdy
    INTEGER         , INTENT(OUT  ) :: ids(:)
    INTEGER         , INTENT(OUT  ) :: idc(:)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    REAL(KIND=r8)   , INTENT(OUT  ) :: todsib
    INTEGER         , INTENT(IN   ) :: idate(:)
    INTEGER         , INTENT(IN   ) :: idatec(:)
    !REAL(KIND=r8)   , INTENT(IN   ) :: si(:)
    !REAL(KIND=r8)   , INTENT(IN   ) :: sl(:)
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(:)
   

    path=path_in
    fNameSibVeg = fNameSibVeg_in
    fNameSibAlb = fNameSibAlb_in

    IF(forcerestore)THEN
       nsoil          = 2
    ELSE
       nsoil          = 6
    END IF

    ALLOCATE(iMaskSiB2  (ibMax,jbMax))
    iMaskSiB2=0
    ALLOCATE(bSoilGrd(ibMax,jbMax)) 
    !bSoilGrd=0.0_r8
    ALLOCATE(NDVI    (ibMax,jbMax,12))
    NDVI=0.0_r8
    ALLOCATE(vcoverg(ibMax,jbMax)) 
    vcoverg=0.0_r8
    sfc%vcover=>vcoverg
    ALLOCATE(morphtab_zc_sib2   (ibMax,jbMax))               
    morphtab_zc_sib2=0.0_r8
    ALLOCATE(morphtab_lwidth_sib2 (ibMax,jbMax))           
    morphtab_lwidth_sib2 =0.0_r8
    ALLOCATE(morphtab_llength_sib2(ibMax,jbMax))            
    morphtab_llength_sib2=0.0_r8
    ALLOCATE(morphtab_laimax_sib2 (ibMax,jbMax))            
    morphtab_laimax_sib2 =0.0_r8
    ALLOCATE(morphtab_stems_sib2  (ibMax,jbMax))            
    morphtab_stems_sib2  =0.0_r8
    ALLOCATE(morphtab_ndvimax_sib2(ibMax,jbMax))            
    morphtab_ndvimax_sib2=0.0_r8
    ALLOCATE(morphtab_ndvimin_sib2(ibMax,jbMax))            
    morphtab_ndvimin_sib2=0.0_r8
    ALLOCATE(morphtab_srmax_sib2  (ibMax,jbMax))            
    morphtab_srmax_sib2  =0.0_r8
    ALLOCATE(morphtab_srmin_sib2  (ibMax,jbMax))             
    morphtab_srmin_sib2=0.0_r8  
    ALLOCATE(aerovar_zo_sib2      (ibMax,50,50,jbMax))
    aerovar_zo_sib2         =0.0_r8    
    ALLOCATE(aerovar_zp_disp_sib2 (ibMax,50,50,jbMax))
    aerovar_zp_disp_sib2 =0.0_r8  
    ALLOCATE(aerovar_rbc_sib2     (ibMax,50,50,jbMax))
    aerovar_rbc_sib2     =0.0_r8  
    ALLOCATE(aerovar_rdc_sib2     (ibMax,50,50,jbMax))
    aerovar_rdc_sib2  =0.0_r8     
    ALLOCATE(chil_sib2             (ibMax,jbMax))         
    chil_sib2         =0.0_r8       
    ALLOCATE(tran_sib2             (ibMax,2,2,jbMax)) 
    tran_sib2         =0.0_r8       
    ALLOCATE(ref_sib2             (ibMax,2,2,jbMax)) 
    ref_sib2         =0.0_r8       
    ALLOCATE(vcover2g(ibMax,jbMax)) 
    vcover2g=0.0_r8
    ALLOCATE(timevar_fpar         (ibMax,jbMax)) 
    timevar_fpar=0.0_r8                 
    ALLOCATE(timevar_lai         (ibMax,jbMax)) 
    timevar_lai         =0.0_r8        
    ALLOCATE(timevar_green         (ibMax,jbMax)) 
    timevar_green=0.0_r8                 
    ALLOCATE(timevar_zo                 (ibMax,jbMax)) 
    timevar_zo          =0.0_r8       
    ALLOCATE(timevar_zp_disp     (ibMax,jbMax)) 
    timevar_zp_disp    =0.0_r8      
    ALLOCATE(timevar_rbc         (ibMax,jbMax)) 
    timevar_rbc           =0.0_r8      
    ALLOCATE(timevar_rdc         (ibMax,jbMax)) 
    timevar_rdc        =0.0_r8         
    ALLOCATE(timevar_gmudmu      (ibMax,jbMax)) 
    timevar_gmudmu   =0.0_r8        
    ALLOCATE(sandfrac  (ibMax,jbMax))
    sandfrac=0.0_r8   
    ALLOCATE(clayfrac  (ibMax,jbMax))
    clayfrac=0.0_r8       
    ALLOCATE(respfactor(ibMax,nsoil+1))
    respfactor=0.0_r8        
    ALLOCATE(pco2m2     (ibMax,jbMax))
    pco2m2=0.0_r8       
    ALLOCATE(pco2ap    (ibMax,jbMax))
    pco2ap=0.0_r8       
    ALLOCATE(d13cca    (ibMax,jbMax))
    d13cca=0.0_r8  
    ALLOCATE(tmgc      (ibMax,jbMax))
    tmgc=0.0_r8  
    ALLOCATE(qmgc      (ibMax,jbMax))
    qmgc=0.0_r8  
    ALLOCATE(tcasc       (ibMax,jbMax))
    tcasc=0.0_r8  
    ALLOCATE(qcasc       (ibMax,jbMax))
    qcasc=0.0_r8  
    ALLOCATE(tcalc       (ibMax,jbMax))
    tcalc=0.0_r8  
    ALLOCATE(tgrdc       (ibMax,jbMax))
    tgrdc=0.0_r8  
    ALLOCATE(tdgc        (ibMax,nsoil,jbMax))
    tdgc =0.0_r8  
    ALLOCATE(stores     (ibMax,jbMax))
    stores=0.0_r8  
    ALLOCATE(snowg       (ibMax,jbMax,2))
    snowg=0.0_r8  
    ALLOCATE(wwwgc       (ibMax,3,jbMax))
    wwwgc=0.0_r8  
    ALLOCATE(ventmf2 (ibMax,jbMax)) 
    ventmf2=0.0_r8  
    ALLOCATE(thvgm2  (ibMax,jbMax)) 
    thvgm2 =0.0_r8  
    ALLOCATE(c4fractg (ibMax,jbMax)) 
    c4fractg =0.0_r8  
    ALLOCATE(d13crespg(ibMax,jbMax)) 
    d13crespg=0.0_r8  
    ALLOCATE(d13ccag  (ibMax,jbMax)) 
    d13ccag  =0.0_r8  
    ALLOCATE(sensf (ibMax,jbMax)) 
    sensf=0.0_r8  
    ALLOCATE(latef (ibMax,jbMax)) 
    latef=0.0_r8  
    ALLOCATE(co2fx(ibMax,jbMax))
    co2fx=0.0_r8 
    ALLOCATE(MskAntSib2(ibMax,jbMax))
    MskAntSib2 =0
    ALLOCATE(radfac_gbl(ibMax,2,2,2,jbMax))
    radfac_gbl=0.0_r8 
    ALLOCATE(AlbGblSiB2(ibMax,2,2,jbMax))
    AlbGblSiB2=0.0_r8 
    ALLOCATE(thermk_gbl(ibMax,jbMax))
    thermk_gbl=0.0_r8 
    ALLOCATE(tgeff4_sib_gbl(ibMax,jbMax))
    tgeff4_sib_gbl=0.0_r8 
    ALLOCATE(glsm_w   (ibMax,jbMax,nzg     ))
    glsm_w=0.0_r8
    ALLOCATE(veg_type (ibMax,jbMax,npatches))
    veg_type=0.0_r8
    ALLOCATE(zlwup_SiB2(ibMax,jbMax))
    zlwup_SiB2=0.0_r8    
    ALLOCATE(tcasm(ibMax,jbMax))
    tcasm=0.0_r8    
    ALLOCATE(tcas0(ibMax,jbMax))
    tcas0=0.0_r8    
    ALLOCATE(qcasm(ibMax,jbMax))
    qcasm=0.0_r8    
    ALLOCATE(qcas0(ibMax,jbMax))
    qcas0=0.0_r8 
    ALLOCATE(frac_occ (ibMax,jbMax,npatches))
    frac_occ=0.0_r8
    ALLOCATE(soil_type(ibMax,jbMax         ))
    soil_type=0.0_r8
    ALLOCATE(tdgm(ibMax,6,jbMax))
    tdgm=0.0_r8 
    ALLOCATE(tdg0(ibMax,6,jbMax))
    tdg0=0.0_r8 
    CALL qsat_data()

    !----------------------------------------------------------------------
    !   !itb...set sandfrac to 0.0 right now, not being used...
    !
    sandfrac = 0.0_r8

    respfactor = 3.0E-6_r8
    !Bio...set values for pco2m and pco2ap here
    pco2m2  = 35.0_r8
    pco2ap = 35.0_r8
    d13cca = -7.8_r8
    !Bio...set initial CAS temp value to canopy veg temp
    !    tcas(:) = tcal(:)

    IF(TRIM(isimp).NE.'YES') THEN
    

    CALL InitBoundCond(iMax,jMax,ibMax,jbMax,kMax,ifdy,&
         ids,idc,ifday,&
       tod,todsib,idate,idatec,record_type,&
       fNameSoilType,fNameVegType,fNameSoilMoist, &
       fNameTg3zrl  ,fNameSoilTab,fNameMorfTab,fNameBioCTab,&
       fNameAeroTab,fNameSiB2Mask,fNameSandMask,fNameClayMask,fNameTextMask,&
       ibMaxPerJB)
    END IF

      
  END SUBROUTINE InitSfcSib2

  SUBROUTINE InitBoundCond(iMax,jMax,&
       ibMax,jbMax,kMax,ifdy,ids,idc,ifday, &
       tod,todsib,idate,idatec,record_type,&
       fNameSoilType,fNameVegType,fNameSoilMoist, &
       fNameTg3zrl  ,fNameSoilTab,fNameMorfTab,fNameBioCTab,&
       fNameAeroTab,fNameSiB2Mask,fNameSandMask,fNameClayMask,fNameTextMask,&
       ibMaxPerJB)
    INTEGER         , INTENT(IN   ) :: iMax
    INTEGER         , INTENT(IN   ) :: jMax
    INTEGER         , INTENT(IN   ) :: ibMax
    INTEGER         , INTENT(IN   ) :: jbMax
    INTEGER         , INTENT(IN   ) :: kMax
    INTEGER         , INTENT(OUT  ) :: ifdy
    INTEGER         , INTENT(OUT  ) :: ids(:)
    INTEGER         , INTENT(OUT  ) :: idc(:)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod

    REAL(KIND=r8)   , INTENT(OUT  ) :: todsib
    INTEGER         , INTENT(IN   ) :: idate(:)
    INTEGER         , INTENT(IN   ) :: idatec(:)

    !REAL(KIND=r8)   , INTENT(IN   ) :: si(:)
    !REAL(KIND=r8)   , INTENT(IN   ) :: sl(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: record_type
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSoilType
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameVegType
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSoilMoist
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSoilTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMorfTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameBioCTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameAeroTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSiB2Mask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSandMask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameClayMask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameTextMask
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(:)
    !CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibmsk
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameTg3zrl

    REAL(KIND=r8)                            :: tice
    REAL(KIND=r8)                            :: sinmax
    INTEGER                         :: j
    INTEGER                         :: i,k,ierr
    INTEGER                         :: ncount,LRecIN,irec
    REAL(KIND=r8)                            :: wsib  (ibMax,jbMax)
    REAL(KIND=r8)                            :: zero
    REAL(KIND=r8)                            :: thousd
    REAL(KIND=r8)            , PARAMETER     :: xl0   =10.0_r8
    REAL(KIND=r8)            , PARAMETER     :: t0 =271.17_r8
    REAL(KIND=r8) :: VegType(iMax,jMax,npatches) ! SIB veg type
    REAL(KIND=r8) ::   buf (iMax,jMax,4)
    REAL(KIND=r4) ::   brf (iMax,jMax)
    INTEGER :: ier(iMax,jMax)
    CHARACTER(LEN=*), PARAMETER :: h='**(InitBoundCondSib2)**'
    tice  =271.16e0_r8
    zero  =0.0e3_r8
    thousd=1.0e3_r8
    brf=0.0_r4
    buf=0.0_r8
    VegType=0.0_r8
    ier=0
    IF(TRIM(isimp).NE.'YES') THEN

       sheleg=0.0_r8
       CALL vegin_sib2(fNameSoilTab,fNameMorfTab,fNameBioCTab,fNameAeroTab)
       CALL ReadSurfaceMaskSib2(ibMax,jbMax,iMax,jMax,fNameSiB2Mask,fNameSandMask,&
                              fNameClayMask,fNameTextMask)

       INQUIRE (IOLENGTH=LRecIN) brf
       OPEN (UNIT=nftgz0,FILE=TRIM(fNameTg3zrl), FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN, &
         ACTION='read', STATUS='OLD', IOSTAT=ierr) 
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameTg3zrl), ierr
          STOP "**(ERROR)**"
       END IF
 
       INQUIRE (IOLENGTH=LRecIN) brf
       OPEN (UNIT=nfzol,FILE=TRIM(fNameRouLen),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN, &
            ACTION='READ', STATUS='OLD', IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameRouLen), ierr
          STOP "**(ERROR)**"
       END IF

       !
       !
       !     initialize sib variables
       !
       IF(ifday == 0 .AND. tod == zero .AND. initlz >= 0 ) THEN

          call MsgOne(h,'Cold start SSib variables')

          CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin  
!climate aerosol parameters
               aod,asy,ssa,z_aer,ifaeros,&
!tar end
!
!tar begin  
!climate aerosol parameters
               aodF,asyF,ssaF,z_aerF, &
!tar end  
               ifday , tod  ,idate ,idatec, &
               ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
               sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
               yrl  ,monl,ibMax,jbMax,ibMaxPerJB)
               
          irec=1
          CALL ReadGetNFTGZ(nftgz0,irec,buf(:,:,1),buf(:,:,2),buf(:,:,3))
          READ (UNIT=nfzol, REC=1) brf
          
          buf(1:iMax,1:jMax,4)=brf(1:iMax,1:jMax)
          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(buf(:,:,1),tg1)
            !CALL AveBoxIJtoIBJB(buf(:,:,1),tg1)
          ELSE
             CALL IJtoIBJB(buf(:,:,1) ,tg1 )
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(buf(:,:,2),tg2)
            !CALL AveBoxIJtoIBJB(buf(:,:,2),tg2)
          ELSE
             CALL IJtoIBJB(buf(:,:,2) ,tg2 )
          END IF

          IF (reducedGrid) THEN
           CALL LinearIJtoIBJB(buf(:,:,3),tg3)
            !CALL AveBoxIJtoIBJB(buf(:,:,3),tg3)
           ELSE
             CALL IJtoIBJB(buf(:,:,3) ,tg3 )
          END IF
 
          IF (reducedGrid) THEN
            CALL LinearIJtoIBJB(buf(:,:,4),zorl)
            !CALL AveBoxIJtoIBJB(buf(:,:,4),zorl)
          ELSE
             CALL IJtoIBJB(buf(:,:,4),zorl )
          END IF
          z0=zorl
          sinmax=150.0_r8
          !
          !     use rvisd as temporary for abs(soilm)
          !
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                rVisDiff(i,j)=ABS(soilm(i,j))
             END DO
          END DO
          !-srf--------------------------------
          IF(iglsm_w == 0) THEN
             CALL sibwet_sib2(ibmax,jbmax,rVisDiff,sinmax,wsib,ssib,mxiter,ibMaxPerJB)
          ELSE
             !
             !- rotina para chamar leitura da umidade do solo
             !
             CALL read_gl_sm_bc(imax           , & !   IN
                  jmax           , & !   IN
                  jbMax          , & !   IN
                  ibMaxPerJB     , & !   IN
                  record_type    , & !   IN
                  fNameSoilType  , & !   IN
                  fNameVegType   , & !   IN
                  fNameSoilMoist , & !   IN
                  VegType          ) !   INOUT

             CALL re_assign_sib_soil_prop(imax            , & ! IN
                  jmax            , & ! IN
                  npatches        , & ! IN
                  VegType         ) ! IN

             !IF (reducedGrid) THEN
             !   CALL FreqBoxIJtoIBJB(imask_in,imask)
             !ELSE
             !   CALL IJtoIBJB( imask_in,imask)
             !END IF

             !
             !- for output isurf, use rlsm array
             !
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   !  rlsm(i,j)=REAL(imask(i,j),r8)
                END DO
             END DO

             CALL sibwet_GLSM (ibMax          , & ! IN
                  jbMax          , & ! IN
                  !imask          , & ! IN
                  wsib           , & ! OUT
                  ssib           , & ! IN
                  mxiter         , & ! IN
                  ibMaxPerJB     , & ! IN
                  soilm          , & ! OUT
                  nzg            , & ! in
                  wsib3d         , & ! OUT
                  glsm_w)            ! IN

          END IF
          !-srf--------------------------------
          ppli=0.0_r8
          ppci=0.0_r8
          capac0=0.0_r8
          capacm=0.0_r8
          !capacc=0.0_r8
          !
          !     td0 (deep soil temp) is temporarily defined as tg3
          !
          !$OMP PARALLEL DO PRIVATE(ncount,i)
          DO j=1,jbMax
             ncount=0
             DO i=1,ibMaxPerJB(j)
                gl0(i,j)=xl0
                Mmlen(i,j)=xl0
                IF(iMaskSiB2(i,j) >= 1_i8)gtsea(i,j)=290.0_r8
                IF(iMaskSiB2(i,j) >= 1_i8)gco2flx(i,j)=0.0_r8
                IF(iMaskSiB2(i,j) == 0_i8)gndvi(i,j)=0.0_r8
                
                tseam(i,j)=gtsea(i,j)
                ndvim(i,j)=gndvi(i,j)
                TSK (I,J)=ABS(gtsea(i,j))
                IF (omlmodel) THEN
                   HML  (i,J)= oml_hml0 - 13.5_r8*log(MAX(ABS(TSK(i,j))-tice+0.01_r8,1.0_r8))
                   HUML (I,J)=0.0_r8
                   HVML (I,J)=0.0_r8
                END IF

                IF(iMaskSiB2(i,j) == 0_i8) THEN
                   IF(-gtsea(i,j).LT.t0) THEN
                      iMaskSiB2(i,j)=-1_i8 
                      iMask    (i,j)=-1_i8
                   END IF
                ELSE
                   ncount=ncount+1
                   IF(iglsm_w == 0) THEN
                      w0 (ncount,1,j)=wsib(i,j)
                      w0 (ncount,2,j)=wsib(i,j)
                      w0 (ncount,3,j)=wsib(i,j)
                      w0    (ncount,1,j)=wsib(i,j)
                      w0    (ncount,2,j)=wsib(i,j)
                      w0    (ncount,3,j)=wsib(i,j)
                   ELSE
                      !-srf--------------------------------
                      w0 (ncount,1,j)=wsib3d(i,j,1)
                      w0 (ncount,2,j)=wsib3d(i,j,2)
                      w0 (ncount,3,j)=wsib3d(i,j,3)
                      w0    (ncount,1,j)=wsib3d(i,j,1)
                      w0    (ncount,2,j)=wsib3d(i,j,2)
                      w0    (ncount,3,j)=wsib3d(i,j,3)
                      !-srf--------------------------------
                   END IF
                   DO k=1,nsoil
                      tdg0   (ncount,k,j)=tg3 (i,j)
                      tdgc   (ncount,k,j)=tg3 (i,j)
                   END DO   
                   sm0(ncount,1,j)=w0(ncount,1,j)*bSoilGrd(i,j)%poros
                   sm0(ncount,2,j)=w0(ncount,2,j)*bSoilGrd(i,j)%poros
                   sm0(ncount,3,j)=w0(ncount,3,j)*bSoilGrd(i,j)%poros

                   IF(iglsm_w == 0) THEN
                      wwwgc(ncount,1,j)=wsib(i,j)
                      wwwgc(ncount,2,j)=wsib(i,j)
                      wwwgc(ncount,3,j)=wsib(i,j)
                      !wm(ncount,1,j)=wsib(i,j)
                      !wm(ncount,2,j)=wsib(i,j)
                      !wm(ncount,3,j)=wsib(i,j)
                      wm   (ncount,1,j)=wsib(i,j)
                      wm   (ncount,2,j)=wsib(i,j)
                      wm   (ncount,3,j)=wsib(i,j)
                   ELSE
                      !-srf--------------------------------
                      wwwgc(ncount,1,j)=wsib3d(i,j,1)
                      wwwgc(ncount,2,j)=wsib3d(i,j,2)
                      wwwgc(ncount,3,j)=wsib3d(i,j,3)                         
                            wm(ncount,1,j)=wsib3d(i,j,1)
                      wm(ncount,2,j)=wsib3d(i,j,2)
                      wm(ncount,3,j)=wsib3d(i,j,3)
                            wm   (ncount,1,j)=wsib3d(i,j,1)
                      wm   (ncount,2,j)=wsib3d(i,j,2)
                      wm   (ncount,3,j)=wsib3d(i,j,3)

                      !-srf--------------------------------
                   END IF
                   DO k=1,nsoil
                      tdgm   (ncount, k, j)=tg3 (i,j)
                   END DO   
                   tgm   (ncount,  j)=tg3 (i,j)

                   tgrdc   (ncount,  j)=tg3 (i,j)
                   tg0   (ncount,  j)=tg3 (i,j)
                   tcm   (ncount,  j)=tg3 (i,j)
                   tc0   (ncount,  j)=tg3 (i,j)
                   
                   tcalc   (ncount,  j)=tg3 (i,j)
                   ssib  (ncount,j  )=0.0_r8
                   IF(soilm(i,j).LT.0.0_r8)ssib(ncount,j)=wsib(i,j)
                  !IF(sheleg(i,j).GT.zero) THEN
                  !    capac0(ncount,2,j)=sheleg(i,j)/thousd
                  !    capacm(ncount,2,j)=sheleg(i,j)/thousd
                  !    capacc(ncount,2,j)=sheleg(i,j)/thousd
                  ! END IF
                END IF
             END DO
          END DO
          !$OMP END PARALLEL DO

       ELSE

          call MsgOne(h,'Warm start SSib variables')

          READ(UNIT=nfsibi)ifdy,todsib,ids,idc
          READ(UNIT=nfsibi) tm0,tcasm,tcasc,tcas0,tmm
          READ(UNIT=nfsibi) qm0,qcas0,qcasc,qcasm,qmm
          READ(UNIT=nfsibi) tdg0 ,tdgc ,tdgm
          READ(UNIT=nfsibi) tg0,tgrdc,tgm
          READ(UNIT=nfsibi) tc0,tcalc,tcm
          READ(UNIT=nfsibi) w0 ,wwwgc,wm   ,wm
          READ(UNIT=nfsibi) capac0,capacc,capacm
          READ(UNIT=nfsibi) ppci  ,ppli,tkemyj
          READ(UNIT=nfsibi) gl0   ,zorl ,gtsea,gco2flx,gndvi,tseam,ndvim,qsfc0,&
          TSfc0,QSfcm,TSfcm,co2fx,pco2ap,c4fractg,& 
          d13crespg,d13ccag,stores,HML,HUML,HVML,TSK,z0sea,&
          TC_SeaIce,TGS_SeaIce,TD_SeaIce,TA_SeaIce,SNOA_SeaIce,SNOB_SeaIce,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 
          READ(UNIT=nfsibi) laymld,       hbath,     tdeep,sdeep
          Mmlen=gl0
          REWIND nfsibi



          IF(initlz < 0.AND. initlz >= -3 )THEN

             IF(initlz == -2 .or. initlz == -3 )ifsst=-1
             IF(ifco2flx == -2 .or. ifco2flx == -3 )ifco2flx=-1

             CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin  
!climate aerosol parameters
                  aod,asy,ssa,z_aer,ifaeros,&
!tar end 
!
!tar begin  
!climate aerosol parameters
                  aodF,asyF,ssaF,z_aerF, &
!tar end   
                  ifday , tod  ,idate ,idatec, &
                  ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
                  sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
                  yrl  ,monl,ibMax,jbMax,ibMaxPerJB)

             IF( initlz == -2  .or. initlz == -3 ) THEN
                !$OMP PARALLEL DO PRIVATE(i)
                DO j=1,jbMax
                   DO i=1,ibMaxPerJB(j)
                      IF(iMaskSiB2(i,j) >= 1_i8) gtsea(i,j)=290.0_r8
                      IF(iMaskSiB2(i,j) >= 1_i8) gco2flx(i,j)=0.0_r8
                      !tseam(i,j) = gtsea(i,j)
                      TSK  (I,J) = ABS(gtsea(i,j))
                      IF (omlmodel) THEN
                         HML  (i,j) = oml_hml0 - 13.5_r8*log(MAX(ABS(TSK(i,j))-tice+0.01_r8,1.0_r8))
                         HUML (I,J)=0.0_r8
                         HVML (I,J)=0.0_r8
                      END IF
                      IF(iMaskSiB2(i,j) == 0_i8) THEN
                         IF(-gtsea(i,j) < t0) THEN
                            iMaskSiB2(i,j)=-1_i8 
                            iMask    (i,j)=-1_i8
                         END IF
                      END IF
                   END DO
                END DO
                !$OMP END PARALLEL DO
             END IF

          END IF

          !$OMP PARALLEL DO PRIVATE(ncount,i)
          DO j=1,jbMax
             ncount=0
             DO i=1,ibMaxPerJB(j)
                IF(iMaskSiB2(i,j) >=1_i8)THEN
                   ncount=ncount+1
                   ssib(ncount,j)=0.0_r8
                   IF(w0(ncount,1,j).LT.0.0_r8)THEN
                      ssib(ncount,j)=ABS(w0(ncount,1,j))
                      w0(ncount,1,j)=ABS(w0(ncount,1,j))
                      w0(ncount,2,j)=ABS(w0(ncount,2,j))
                      w0(ncount,3,j)=ABS(w0(ncount,3,j))
                      wm(ncount,1,j)=ABS(wm(ncount,1,j))
                      wm(ncount,2,j)=ABS(wm(ncount,2,j))
                      wm(ncount,3,j)=ABS(wm(ncount,3,j))
                   END IF
                END IF
 
                IF(iMaskSiB2(i,j) >=1_i8)THEN
                   ncount=ncount+1
                   ssib(ncount,j)=0.0_r8
                   IF(w0(ncount,1,j).LT.0.0_r8)THEN
                      ssib(ncount,j)=ABS(w0(ncount,1,j))
                      w0(ncount,1,j)=ABS(w0(ncount,1,j))
                      w0(ncount,2,j)=ABS(w0(ncount,2,j))
                      w0(ncount,3,j)=ABS(w0(ncount,3,j))
                      wwwgc(ncount,1,j)=ABS(wwwgc(ncount,1,j))
                      wwwgc(ncount,2,j)=ABS(wwwgc(ncount,2,j))
                      wwwgc(ncount,3,j)=ABS(wwwgc(ncount,3,j))
                      
                      wm(ncount,1,j)=ABS(wm(ncount,1,j))
                      wm(ncount,2,j)=ABS(wm(ncount,2,j))
                      wm(ncount,3,j)=ABS(wm(ncount,3,j))
                   END IF
                   sm0(ncount,1,j)=w0(ncount,1,j)*bSoilGrd(i,j)%poros
                   sm0(ncount,2,j)=w0(ncount,2,j)*bSoilGrd(i,j)%poros
                   sm0(ncount,3,j)=w0(ncount,3,j)*bSoilGrd(i,j)%poros
                END IF
 
             END DO
          END DO
          !$OMP END PARALLEL DO

          IF(nfctrl(5).GE.1)WRITE(UNIT=nfprt,FMT=444) ifdy,todsib,ids,idc

       END IF
    END IF

444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
!555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
!         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
  END SUBROUTINE InitBoundCond




  SUBROUTINE InitCheckSiB2File(iMax,jMax,ibMax,&
       jbMax  ,kMax, ifdy ,ids   ,idc   ,ifday , &
       tod   ,idate ,idatec   ,todsib,ibMaxPerJB )   
    INTEGER, INTENT(IN   ) :: iMax
    INTEGER, INTENT(IN   ) :: jMax
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(OUT  ) :: ifdy
    INTEGER, INTENT(OUT  ) :: ids   (:)
    INTEGER, INTENT(OUT  ) :: idc   (:)
    INTEGER, INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    INTEGER, INTENT(IN   ) :: idate (:)
    INTEGER, INTENT(IN   ) :: idatec(:)
    REAL(KIND=r8)   , INTENT(OUT  ) :: todsib
    INTEGER, INTENT(IN   ) :: ibMaxPerJB(:)

    INTEGER                :: j
    INTEGER                :: ncount
    INTEGER                :: i
    REAL(KIND=r8)                   :: tice  =271.16e0_r8
    CHARACTER(LEN=*), PARAMETER :: h="**(InitCheckSiB2File)**"

    !
    !     read cloud dataset for cold start
    !
    IF(initlz < 0)THEN

       CALL MsgOne(h,'Read SSib variables from warm-start file')

       READ(UNIT=nfsibi) ifdy,todsib,ids,idc
       READ(UNIT=nfsibi) tm0,tcasm,tcasc,tcas0,tmm
       READ(UNIT=nfsibi) qm0,qcas0,qcasc,qcasm,qmm
       READ(UNIT=nfsibi) tdg0,tdgc   ,tdgm
       READ(UNIT=nfsibi) tg0,tgrdc, tgm
       READ(UNIT=nfsibi) tc0,tcalc  ,tcm
       READ(UNIT=nfsibi) w0 ,wwwgc   ,wm,wm
       READ(UNIT=nfsibi) capac0,capacc,capacm
       READ(UNIT=nfsibi) ppci  ,ppli,tkemyj
       READ(UNIT=nfsibi) gl0   ,zorl  ,gtsea,gco2flx,gndvi,tseam,ndvim,QSfc0,&
       TSfc0,QSfcm,TSfcm,co2fx,pco2ap,c4fractg,& 
       d13crespg,d13ccag,stores,HML,HUML,HVML,TSK,z0sea,&
       TC_SeaIce,TGS_SeaIce,TD_SeaIce,TA_SeaIce,SNOA_SeaIce,SNOB_SeaIce,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 
       READ(UNIT=nfsibi)  laymld,       hbath,     tdeep,sdeep
       Mmlen=gl0
       REWIND nfsibi

       CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin  
! climateaerosol parameters
            aod,asy,ssa,z_aer,ifaeros,&
!tar end 
!
!tar begin  
! climateaerosol parameters
            aodF,asyF,ssaF,z_aerF,&
!tar end       
            ifday , tod  ,idate ,idatec, &
            ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
            sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
            yrl  ,monl,ibMax,jbMax,ibMaxPerJB)

       DO j=1,jbMax
          ncount=0
          DO i=1,ibMaxPerJB(j)
          
             IF(iMaskSiB2(i,j) >=1_i8)THEN
                ncount=ncount+1
                ssib(ncount,j)=0.0_r8
                IF(w0(ncount,1,j).LT.0.0_r8)THEN
                   ssib(ncount,  j)=ABS(w0(ncount,1,j))
                   w0  (ncount,1,j)=ABS(w0(ncount,1,j))
                   w0  (ncount,2,j)=ABS(w0(ncount,2,j))
                   w0  (ncount,3,j)=ABS(w0(ncount,3,j))                   
                   wm  (ncount,1,j)=ABS(wm(ncount,1,j))
                   wm  (ncount,2,j)=ABS(wm(ncount,2,j))
                   wm  (ncount,3,j)=ABS(wm(ncount,3,j))
                END IF
             END IF

             IF(iMaskSiB2(i,j) >= 1_i8)THEN
                ncount=ncount+1
                ssib(ncount,j)=0.0_r8
                IF(w0(ncount,1,j).LT.0.0_r8)THEN
                   ssib(ncount,  j)=ABS(w0(ncount,1,j))
                   w0  (ncount,1,j)=ABS(w0(ncount,1,j))
                   w0  (ncount,2,j)=ABS(w0(ncount,2,j))
                   w0  (ncount,3,j)=ABS(w0(ncount,3,j))
                   wwwgc  (ncount,1,j)=ABS(wwwgc(ncount,1,j))
                   wwwgc  (ncount,2,j)=ABS(wwwgc(ncount,2,j))
                   wwwgc  (ncount,3,j)=ABS(wwwgc(ncount,3,j))                   
                   wm  (ncount,1,j)=ABS(wm(ncount,1,j))
                   wm  (ncount,2,j)=ABS(wm(ncount,2,j))
                   wm  (ncount,3,j)=ABS(wm(ncount,3,j))
                END IF
             END IF
             
          END DO
       END DO
       IF(nfctrl(5).GE.1)WRITE(UNIT=nfprt,FMT=444) ifdy,todsib,ids,idc
    END IF

444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
!555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
!         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)

  END SUBROUTINE InitCheckSiB2File

  SUBROUTINE Finalize_SiB2()
    DEALLOCATE(AlbGblSiB2)
    DEALLOCATE(MskAntSiB2)   
  END SUBROUTINE Finalize_SiB2


  SUBROUTINE ReStartSiB2 (jbMax,ifday,tod,idate ,idatec, &
       nfsibo,ibMaxPerJB)

    INTEGER           ,INTENT(IN   ) :: jbMax
    INTEGER           ,INTENT(IN   ) :: ifday
    REAL(KIND=r8)              ,INTENT(IN   ) :: tod
    INTEGER           ,INTENT(IN   ) :: idate(:)
    INTEGER           ,INTENT(IN   ) :: idatec(:)
    INTEGER           ,INTENT(IN   ) :: nfsibo
    INTEGER           ,INTENT(IN   ) :: ibMaxPerJB(:)
    INTEGER                         :: i
    INTEGER                         :: j
    INTEGER                         :: ncount

    IF(TRIM(isimp).NE.'YES') THEN
  
       CALL MsgOne('**(restartphyscs)**','Saving physics state for restart')

       !$OMP DO PRIVATE(ncount, i)
       DO j=1,jbMax
          ncount=0
          DO i=1,ibMaxPerJB(j)
             IF(iMaskSiB2(i,j) >= 1_i8)THEN
                ncount=ncount+1
                IF(ssib(ncount,j).GT.0.0_r8)THEN
                   w0  (ncount,1,j)=-ssib(ncount,j)
                   w0  (ncount,2,j)=-ssib(ncount,j)
                   w0  (ncount,3,j)=-ssib(ncount,j)
                   wwwgc(ncount,1,j)=-ssib(ncount,j)
                   wwwgc(ncount,2,j)=-ssib(ncount,j)
                   wwwgc(ncount,3,j)=-ssib(ncount,j)
                   wm(ncount,1,j)=-ssib(ncount,j)
                   wm(ncount,2,j)=-ssib(ncount,j)
                   wm(ncount,3,j)=-ssib(ncount,j)


                   w0  (ncount,1,j)=-ssib(ncount,j)
                   w0  (ncount,2,j)=-ssib(ncount,j)
                   w0  (ncount,3,j)=-ssib(ncount,j)
                   wm  (ncount,1,j)=-ssib(ncount,j)
                   wm  (ncount,2,j)=-ssib(ncount,j)
                   wm  (ncount,3,j)=-ssib(ncount,j)
                END IF
             END IF
          END DO
       END DO


       !$OMP SINGLE
       WRITE(UNIT=nfsibo) ifday,tod,idate,idatec
       WRITE(UNIT=nfsibo) tm0,tcasm,tcasc,tcas0,tmm
       WRITE(UNIT=nfsibo) qm0,qcas0,qcasc,qcasm,qmm
       WRITE(UNIT=nfsibo) tdg0,tdgc,tdgm
       WRITE(UNIT=nfsibo) tg0,tgrdc,tgm
       WRITE(UNIT=nfsibo) tc0,tcalc,tcm
       WRITE(UNIT=nfsibo) w0 ,wwwgc,wm,wm
       WRITE(UNIT=nfsibo) capac0,capacc,capacm
       WRITE(UNIT=nfsibo) ppci,ppli,tkemyj
       WRITE(UNIT=nfsibo) gl0 ,zorl,gtsea,gco2flx,gndvi,tseam,ndvim,QSfc0,&
       TSfc0,QSfcm,TSfcm,co2fx,pco2ap,c4fractg,& 
       d13crespg,d13ccag,stores,HML,HUML,HVML,TSK,z0sea,&
       TC_SeaIce,TGS_SeaIce,TD_SeaIce,TA_SeaIce,SNOA_SeaIce,SNOB_SeaIce,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 
       WRITE(UNIT=nfsibo)  laymld,       hbath,     tdeep,sdeep
       !$OMP END SINGLE
    END IF

  END SUBROUTINE ReStartSiB2



  SUBROUTINE InitSurfTempSiB2 (jbMax ,ibMaxPerJB)

    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: ibMaxPerJB(:)
    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: ncount
    REAL(KIND=r8)                   :: zero  =0.0e3_r8
    REAL(KIND=r8)                   :: thousd=1.0e3_r8
    REAL(KIND=r8)                   :: tf    =273.16e0_r8
    capacm=0.0_r8
    capac0=0.0_r8
                !  capacc=0.0_r8
    DO j=1,jbMax
       ncount=0
       DO i=1,ibMaxPerJB(j)
          IF(iMaskSiB2(i,j) >= 1_i8) THEN
             ncount=ncount+1
             IF(sheleg(i,j).GT.zero) THEN
                capac0(ncount,2,j) = sheleg(i,j)/thousd
                !capacc(ncount,2,j) = sheleg(i,j)/thousd
                capacm  (ncount,2,j) = sheleg(i,j)/thousd
                tgrdc   (ncount,  j) = MIN(tgrdc(ncount,j),tf-0.01e0_r8)
                                
                capac0(ncount,2,j) = sheleg(i,j)/thousd
                tg0   (ncount,  j) = MIN(tg0(ncount,j),tf-0.01e0_r8)
                tgm   (ncount,  j) = MIN(tgm(ncount,j),tf-0.01e0_r8)
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE InitSurfTempSiB2






  ! sibwet :transform mintz-serafini and national meteoroLOGICAL center fields
  !         of soil moisture into sib compatible fields of soil moisture.
  SUBROUTINE sibwet_GLSM (ibMax          , & ! IN
       jbMax          , & ! IN
       !imask          , & ! IN
       wsib           , & ! IN
       ssib           , & ! IN
       mxiter         , & ! OUT
       ibMaxPerJB     , & ! OUT
       soilm          , & ! in
       nzg         , & ! in
       wsib3d         , & ! OUT
       glsm_w)            ! IN

    !
    ! $Author: pkubota $
    ! $Date: 2008/04/09 12:42:57 $
    ! $Revision: 1.12 $
    !
    ! sibwet :transform mintz-serafini and national meteoroLOGICAL center fields
    !         of soil moisture into sib compatible fields of soil moisture.
    !
    !     piers sellers : 29 april 1987
    !
    INTEGER, INTENT(IN   )            :: ibMax
    INTEGER, INTENT(IN   )            :: jbMax
    INTEGER, INTENT(IN   )            :: mxiter
    REAL(KIND=r8)   , INTENT(OUT  )            :: soilm          (ibMax,jbMax)
    !INTEGER(KIND=i8), INTENT(IN   )            :: imask          (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  )            :: wsib           (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  )            :: ssib           (ibMax,jbMax)
    INTEGER, INTENT(in   )            :: ibMaxPerJB     (:)
    INTEGER, INTENT(in   )            :: nzg
    REAL(KIND=r8)   , INTENT(OUT  )            :: wsib3d    (ibMax,jbMax,8       )
    REAL(KIND=r8)   , INTENT(IN   )            :: glsm_w    (ibMax,jbMax,nzg     )

    REAL(KIND=r8)               :: sm  (ityp,mxiter)
    REAL(KIND=r8)               :: time(ityp,mxiter)
    REAL(KIND=r8)               :: fact(ityp,mxiter)
    !
    !-srf
    !
    INTEGER, PARAMETER :: nzgmax=20
    REAL(KIND=r8)               :: glsm_w1d  (0:nzgmax)     ! dummy 1d initial soil  wetness
    REAL(KIND=r8)               :: glsm_tzdep(0:3)          ! sib soil levels
    REAL(KIND=r8)               :: glsm_w_sib(0:3)          ! SIB dummy 1d initial and interpolated soil  wetness
    !
    !-srf
    !
    REAL(KIND=r8)               :: tzdep (3)
    REAL(KIND=r8)               :: tzltm (2)
    REAL(KIND=r8)               :: sibmax(ityp)
    INTEGER            :: k
    REAL(KIND=r8)               :: fx
    INTEGER            :: lonmax
    INTEGER            :: latmax
    INTEGER            :: is
    REAL(KIND=r8)               :: tphsat
    REAL(KIND=r8)               :: tbee
    REAL(KIND=r8)               :: tporos
    INTEGER            :: imm1
    INTEGER            :: imm2
    INTEGER            :: im
    INTEGER            :: imm
    INTEGER            :: ivegm
    REAL(KIND=r8)               :: cover
    REAL(KIND=r8)               :: tph1
    REAL(KIND=r8)               :: tph2
    REAL(KIND=r8)               :: sref
    REAL(KIND=r8)               :: smin
    REAL(KIND=r8)               :: dssib
    REAL(KIND=r8)               :: dw
    REAL(KIND=r8)               :: times
    REAL(KIND=r8)               :: soilmo
    REAL(KIND=r8)               :: w
    REAL(KIND=r8)               :: rsoilm
    INTEGER            :: iter
    REAL(KIND=r8)               :: psit
    REAL(KIND=r8)               :: factor
    REAL(KIND=r8)               :: dt
    INTEGER            :: lat
    INTEGER            :: lon

    !
    !              wsinp    = m-s or nmc fractional wetness
    !              ms       = 1, mintz-serafini
    !              nmc      = 1, national meteoroLOGICAL center
    !              bee      = sib : soil moisture potential factor
    !              phsat     = sib : soil potential at saturation (m)
    !              zdepth(3)= sib : depth of 3 soil layers (m)
    !              poros    = sib : soil porosity
    !              ph1      = sib : leaf potential, stress onset (m)
    !              ph2      = sib : leaf potential, no e-t (m)
    !
    !   output :   wsibt    = sib : fractional wetness
    !              ssibt    = sib : soil moisture content (m)
    !              psit     = sib : soil moisture potential (m)
    !              factor   = sib : extraction factor
    !
    REAL(KIND=r8), PARAMETER :: xph1(13,2) = RESHAPE( &
         (/-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8,-120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8, &
         -200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8,-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8, &
         -120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8,-200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8/), &
         (/13,2/))
    REAL(KIND=r8), PARAMETER :: xph2(13,2) = RESHAPE( &
         (/-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8, &
         -400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8,-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8, &
         -230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8,-400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8/) , &
         (/13,2/))

    !-srf
    !hmjb
    !    REAL, PARAMETER :: glsm_slz(0:nzgmax) = (/  0., 0.1, 0.25, 0.5, 1., 2., 3.,& !7  values
    !         0., 0.,  0.,   0.,  0., 0., 0., 0., 0., 0.,& !10 values
    !         0., 0.,  0.,   0.                         /) !4  values
    !versao para NZG=8 => 9 niveis no MCGA
    REAL(KIND=r8), PARAMETER :: glsm_slz(0:nzgmax) = (/  0.0_r8, 0.05_r8, 0.13_r8, 0.25_r8, 0.5_r8, 1.0_r8, 1.75_r8,& !9  values
         2.5_r8, 4.5_r8,  0.0_r8,   0.0_r8,  0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,& !10 values
         0.0_r8, 0.0_r8,  0.0_r8,   0.0_r8                         /) !4  values
    !-srf
    sm =0.0_r8
    time  =0.0_r8
    fact  =0.0_r8
    ssib  =0.0_r8
    wsib  =0.0_r8

    lonmax=ibMax
    latmax=jbMax

    DO is = 1,ityp
       tzdep (1)= zdepth_r4(is,1)
       tzdep (2)= zdepth_r4(is,2)
       tzdep (3)= zdepth_r4(is,3)
       tphsat   = phsat_r4 (is)
       tbee     = bee_r4   (is)
       tporos   = poros_r4 (is)
       imm1=1
       imm2=1
       tzltm(1)=zlt_r4(is,1,1)
       tzltm(2)=zlt_r4(is,1,2)
       DO im=2,12
          IF (tzltm(1).le.zlt_r4(is,im,1) ) THEN
             imm1=im
             tzltm(1)=zlt_r4(is,im,1)
          END IF

          IF (tzltm(2).le.zlt_r4(is,im,2) ) THEN
             imm2=im
             tzltm(2)=zlt_r4(is,im,2)
          END IF
       END DO

       imm=imm1
       ivegm=1

       IF (tzltm(1).le.tzltm(2)) THEN
          imm=imm2
          ivegm=2
       END IF
       cover=xcover_r4(is,imm,ivegm)
       tph1=xph1(is,ivegm)
       tph2=xph2(is,ivegm)
       !
       !srf- max water content
       !
       sibmax(is) = ( tzdep(1) + tzdep(2) + tzdep(3) ) * tporos
       IF (nfctrl(83).ge.1) WRITE(UNIT=nfprt,FMT=999)is,sibmax(is),tzdep(1), &
            tzdep(2),tzdep(3),tporos
       sref = sibmax(is) * exp( log(tphsat /(-1.0e0_r8)) /tbee )
       smin = sibmax(is) * exp( log(tphsat /(-1.0e10_r8)) /tbee )
       dssib= (sref - smin) / REAL(mxiter,r8)
       dw   = dssib / sibmax(is)
       times  = 0.0e0_r8
       soilmo = sref
       w      = soilmo / sibmax(is)
       rsoilm = 101840.0_r8 * (1.0_r8 - w**0.0027_r8)

       DO iter = 1, mxiter
          CALL extrak      ( w     , &  ! IN
               dw    , &  ! IN
               tbee  , &  ! IN
               tphsat, &  ! IN
               rsoilm, &  ! IN
               cover , &  ! IN
               tph1  , &  ! IN
               tph2  , &  ! IN
               psit  , &  ! OUT
               factor  )  ! OUT
          dt            = dssib / factor
          soilmo        = soilmo - dssib
          w             = soilmo / sibmax(is)
          times         = times + dt
          sm  (is,iter) = soilmo
          time(is,iter) = times
          fact(is,iter) = factor
       END DO
    END DO
    !
    !     input soil moisture map is now transformed to sib fields.
    !
    DO lat = 1, latmax
       DO lon = 1, ibMaxPerJB(lat)

          wsib3d(lon,lat,:) = 0.e0_r8

          is=INT(iMaskSiB2(lon,lat),kind=i4)
          IF (is.ne.0) THEN

             tzdep (1)= zdepth_r4(is,1)
             tzdep (2)= zdepth_r4(is,2)
             tzdep (3)= zdepth_r4(is,3)
             tphsat   = phsat_r4 (is)
             tbee     = bee_r4   (is)
             tporos   = poros_r4 (is)
             !
             !-sib soil levels
             !
             glsm_tzdep(0) = 0.e0_r8
             glsm_w_sib(0) = 0.e0_r8

             DO k=1,3
                glsm_tzdep (k) = zdepth_r4(is,k) + glsm_tzdep (k-1)
                glsm_w_sib (k) = 0.e0_r8
             END DO
             !
             !- copy 3d soil moisture array to 1d column array
             !
             DO k=1,nzg
                glsm_w1d(k)=glsm_w(lon,lat,k)
             END DO
             !
             !- performs vertical interpolation from soil moisture
             !  levels to sib levels
             !
             CALL vert_interp(4               , &  ! IN
                  nzg+1           , &  ! IN
                  glsm_tzdep(0:3) , &  ! IN
                  glsm_slz(0:nzg) , &  ! IN
                  glsm_w1d(0:nzg) , &  ! IN
                  glsm_w_sib(0:3)  )   ! OUT


             !endif
             !
             !- stores 1d sib soil moisture at 3d array
             !
             DO k=1,3
                wsib3d(lon,lat,k) = glsm_w_sib(k)
             END DO
             !
             !------------------------- remove this later--------------------------------X
             !- for now fix zero soil moisture inside the land
             !- latter fix this at soil moisture original data
             !
             !IF (imask(lon,lat) > 0 ) THEN
             !   ssm=0.
             !   DO k=1,3
             !      ssm=ssm+wsib3d(lon,lat,k)
             !   END DO
             !
             !   IF (ssm < 0.15) THEN
             !      !
             !      !print*,'SM null inside land portion', imask(lon,lat)
             !      !print*,'1',lon,lat,wsib3d(lon,lat,:)
             !      !
             !      ssm1d(:) = 0.
             !      ncount = 0
             !      DO i=max(1,lon-4),min(lonmax,lon+4)
             !         DO j=max(1,lat-4),min(latmax,lat+4)
             !    IF (imask(i,j) > 0) THEN !only points inside the land
             !       ssm=0.
             !       DO k=1,3
             !  ssm=ssm+wsib3d(i,j,k)
             !       END DO
             !
             !       IF (ssm > 0.15) THEN
             !  ncount=ncount  + 1
             !  ssm1d(:) = ssm1d(:) + wsib3d(i,j,:)
             !       END IF
             !    END IF
             !         END DO
             !      END DO
             !
             !     IF (ncount > 1) THEN
             !         wsib3d(lon,lat,:)=ssm1d(:)/float(ncount)
             !      ELSE
             !         wsib3d(lon,lat,:)=0.5
             !      END IF
             !      !
             !      !print*,'2',lon,lat,wsib3d(lon,lat,:)
             !      !
             !   END IF
             !END IF
             !
             !-----------------------------------------------------------------------------X
             !
             ssib(lon,lat) = 0.0_r8
             wsib(lon,lat) = 0.0_r8

             DO k=1,3

                fx            = ( glsm_tzdep(k)-glsm_tzdep(k-1) ) / glsm_tzdep(3)
                wsib(lon,lat) = wsib(lon,lat) + glsm_w_sib(k) * fx
                ssib(lon,lat) = ssib(lon,lat) + glsm_w_sib(k) * fx * tporos

             END DO
             !
             ! total water in mm
             !
             soilm(lon,lat) = ( tzdep(1)*wsib3d(lon,lat,1) + &
                  tzdep(2)*wsib3d(lon,lat,2) + &
                  tzdep(3)*wsib3d(lon,lat,3) ) * tporos
             !
          END IF
       END DO
    END DO
999 FORMAT(' IS,MAX,D1,D2,D3,POROS=',I2,1X,5E12.5)
  END SUBROUTINE sibwet_GLSM


 !
  !------------------------------------------------------------
  !
  SUBROUTINE re_assign_sib_soil_prop(iMax            , & ! IN
       jMax            , & ! IN
       npatches        , & ! IN
       veg_type          ) ! IN
    INTEGER, INTENT(IN   )  :: iMax
    INTEGER, INTENT(IN   )  :: jMax
    INTEGER, INTENT(IN   )  :: npatches

    REAL(KIND=r8)   , INTENT(IN   )  :: veg_type (imax,jmax,npatches)
    REAL(KIND=r8)    :: GSWP_soil_input_data(10,12  )
    INTEGER(KIND=i8)  :: imask_in    (imax,jmax         )
    INTEGER :: nnn
    INTEGER :: i
    INTEGER :: j
    !
    !-------------------------------Soil data from GSWP-2 -------------------------------------
    !
    DATA GSWP_soil_input_data/  &
                                !1     2    3     4        5       6      7      8      9     10
                                !SAND(%) SILT CLAY QUARTZ  Wfc    Wwilt  Wsat    b    PHIsat  Ksat
         92.0_r8, 5.0_r8, 3.0_r8,0.92_r8,0.132_r8,0.033_r8,0.373_r8, 3.30_r8,-0.05_r8,2.45E-05_r8,&!1  Sand
         82.0_r8,12.0_r8, 6.0_r8,0.82_r8,0.156_r8,0.051_r8,0.386_r8, 3.80_r8,-0.07_r8,1.75E-05_r8,&!2  Loamy Sand
         58.0_r8,32.0_r8,10.0_r8,0.60_r8,0.196_r8,0.086_r8,0.419_r8, 4.34_r8,-0.16_r8,8.35E-06_r8,&!3  Sandy Loam
         10.0_r8,85.0_r8, 5.0_r8,0.25_r8,0.361_r8,0.045_r8,0.471_r8, 3.63_r8,-0.84_r8,1.10E-06_r8,&!4  Silt Loam
         17.0_r8,70.0_r8,13.0_r8,0.40_r8,0.270_r8,0.169_r8,0.476_r8, 5.25_r8,-0.65_r8,2.36E-06_r8,&!5  Loam
         58.0_r8,15.0_r8,27.0_r8,0.60_r8,0.253_r8,0.156_r8,0.412_r8, 7.32_r8,-0.12_r8,6.31E-06_r8,&!6  Sandy Clay Loam
         32.0_r8,34.0_r8,34.0_r8,0.10_r8,0.301_r8,0.211_r8,0.447_r8, 8.34_r8,-0.28_r8,2.72E-06_r8,&!7  Silty Clay Loam
         10.0_r8,56.0_r8,34.0_r8,0.35_r8,0.334_r8,0.249_r8,0.478_r8, 8.41_r8,-0.63_r8,1.44E-06_r8,&!8  Clay Loam
         52.0_r8, 6.0_r8,42.0_r8,0.52_r8,0.288_r8,0.199_r8,0.415_r8, 9.70_r8,-0.12_r8,4.25E-06_r8,&!9  Sandy Clay
         6.0_r8,47.0_r8,47.0_r8,0.10_r8,0.363_r8,0.286_r8,0.478_r8,10.78_r8,-0.58_r8,1.02E-06_r8,&!10 Silty Clay
         22.0_r8,20.0_r8,58.0_r8,0.25_r8,0.353_r8,0.276_r8,0.450_r8,12.93_r8,-0.27_r8,1.33E-06_r8,&!11 Clay
         43.0_r8,39.0_r8,18.0_r8,0.10_r8,0.250_r8,0.148_r8,0.437_r8, 5.96_r8,-0.24_r8,4.66E-06_r8 /!12 Silt
    !
    !-srf: avoid this for now, only use it when all arrays above are used like:
    ! bee(int(soil_type(lon,lat))) and not the usual way: bee(isurf(lon,lat))),
    ! where isurf is the vegetation index
    !
    !GO TO 332
    DO nnn = 1,12
       !
       !   sslfc(nnn)  = GSWP_soil_input_data(5,nnn)        !not in use
       !   sswlts(nnn) = GSWP_soil_input_data(6,nnn)        !not in use
       !   sswlts(nnn) = max(0.06_r8,GSWP_soil_input_data(6,n) !not in use nn)
       !
       ! print*,nnn,'poros bee phsat satco'
       ! print*,poros(nnn) , GSWP_soil_input_data(7,nnn)
       ! print*,bee(nnn)   ,GSWP_soil_input_data(8,nnn)
       ! print*,phsat(nnn) ,GSWP_soil_input_data(9,nnn)
       ! print*,satco(nnn) ,GSWP_soil_input_data(10,nnn)

       ! poros(nnn) = GSWP_soil_input_data(7,nnn)
       ! bee  (nnn) = GSWP_soil_input_data(8,nnn)
       ! phsat(nnn) = GSWP_soil_input_data(9,nnn)
       ! satco(nnn) = GSWP_soil_input_data(10,nnn)
    END DO

    !332 continue
    !
    !- for now, set isurf(:,:) as the veg data of the predominant biome:
    !
    DO j=1,jMax
       DO i= 1,iMax
          !imask_in(i,j) = 0 => ocean  / imask_in(i,j) = 13 => ice
          !IF (imask_in(i,j) > 0_i8 .and. imask_in(i,j) < 13_i8) THEN
             !print*,'1',i,j,int(veg_type(i,j,2)),imask_in(i,j)
             imask_in(i,j) = int(veg_type(i,j,2))
          !END IF
       END DO
    END DO
    
    IF (reducedGrid) THEN
       CALL FreqBoxIJtoIBJB(imask_in,iMaskSiB2)
    ELSE
       CALL IJtoIBJB( imask_in,iMaskSiB2)
    END IF
    iMask = iMaskSiB2
    !
    !stop 44433
    !srf- original SSIB from MCGA requires 13 soil classes, while USDA/GSWP2 has only 12
    !srf- the soil class 13 is not changed here (see vegin.f90)
    !  bee(13) = 4.8_r8
    !  phsat(13) = -0.167_r8
    !  satco(13) = 0.762e-4_r8
    !  poros(13) = 0.4352_r8
    !  zdepth(13,1) = 1.0_r8
    !  zdepth(13,2) = 1.0_r8
    !  zdepth(13,3) = 1.0_r8
    !
    RETURN
  END SUBROUTINE re_assign_sib_soil_prop



  !------------------------------------------------------------
  SUBROUTINE read_gl_sm_bc(iMax           , &!   IN
       jMax           , &!   IN
       jbMax          , &!   IN
       ibMaxPerJB     , &!   IN
       record_type   , &! IN
       fNameSoilType  , &! IN
       fNameVegType   , &! IN
       fNameSoilMoist , &
       VegType  )! IN

    INTEGER, INTENT(IN   )            :: iMax
    INTEGER, INTENT(IN   )            :: jMax
    INTEGER, INTENT(IN   )            :: jbMax
    INTEGER, INTENT(IN   )            :: ibMaxPerJB(:)
    CHARACTER(LEN=*), INTENT(IN   )   :: record_type
    CHARACTER(LEN=*), INTENT(IN   )   :: fNameSoilType
    CHARACTER(LEN=*), INTENT(IN   )   :: fNameVegType
    CHARACTER(LEN=*), INTENT(IN   )   :: fNameSoilMoist
    REAL(KIND=r8) , INTENT(INOUT   ):: VegType(iMax,jMax,npatches) ! SIB veg type
    REAL(KIND=r8) :: FracOcc(iMax,jMax,npatches) ! fractional area
    REAL(KIND=r8) :: glsm(iMax,jMax,nzg     )   ! initial soil wetness data
    REAL(KIND=r8) :: SoilType(iMax,jMax         )! FAO/USDA soil texture
    !
    ! Local
    !
    INTEGER            :: i
    INTEGER            :: j
    INTEGER            :: k
    INTEGER            :: ipatch
    INTEGER            :: ierr
    REAL(KIND=r8)               :: fractx
    !
    !------------------------- soil type initialization ------------
    !
    call MsgOne('**(read_gl_sm_bc)**','Opening GL soil file='//TRIM(fNameSoilType))
    FracOcc=0.0_r8
    glsm=0.0_r8
    SoilType=0.0_r8
    IF (record_type == 'seq') THEN      !sequential mode

       OPEN(UNIT=nfsoiltp,FILE=TRIM(fNameSoilType),FORM='unformatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilType), ierr
          STOP "**(ERROR)**"
       END IF
       READ(UNIT=nfsoiltp) ((SoilType(i,j),i=1,iMax),j=1,jMax)

       IF (reducedGrid) THEN
          CALL NearestIJtoIBJB(SoilType,soil_type)
       ELSE
          CALL IJtoIBJB( SoilType,soil_type)
       END IF

       CLOSE(UNIT=nfsoiltp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfsoiltp,FILE=TRIM(fNameSoilType),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilType), ierr
          STOP "**(ERROR)**"
       END IF

       CALL  vfirec(nfsoiltp,SoilType,imax*jmax,'LIN')

       IF (reducedGrid) THEN
          CALL NearestIJtoIBJB(SoilType,soil_type)
       ELSE
          CALL IJtoIBJB( SoilType,soil_type)
       END IF

       CLOSE(UNIT=nfsoiltp)

    END IF
    DO i=1,iMax
       DO j=1,jMax
          SoilType(i,j)=REAL(INT(SoilType(i,j)+0.1_r8),r8)
       END DO
    END DO
    DO j=1,jbMax
       DO i=1,ibMaxPerJB(j)
          soil_type(i,j)=REAL(INT(soil_type(i,j)+0.1_r8),r8)
       END DO
    END DO
    !
    !-------------------veg type and fractional area initialization ------------
    !
    call MsgOne('**(read_gl_sm_bc)**','Opening GL veg file='//TRIM(fNameVegType))

    IF (record_type == 'seq') THEN !sequential mode

       OPEN(UNIT=nfvegtp,FILE=TRIM(fNameVegType),FORM='unformatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameVegType), ierr
          STOP "**(ERROR)**"
       END IF

       DO ipatch=1,npatches_actual

          READ(UNIT=nfvegtp) ((VegType(i,j,ipatch),i=1,iMax),j=1,jMax) !veg dominante no patch
          READ(UNIT=nfvegtp) ((FracOcc(i,j,ipatch),i=1,iMax),j=1,jMax) !fracao ocupada pelo patch

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(VegType(:,:,ipatch) ,veg_type(:,:,ipatch) )
          ELSE
             CALL IJtoIBJB( VegType(:,:,ipatch) ,veg_type(:,:,ipatch) )
          END IF

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch) )
          ELSE
             CALL IJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch)  )
          END IF

       END DO

       CLOSE(UNIT=nfvegtp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfvegtp,FILE=TRIM(fNameVegType),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameVegType), ierr
          STOP "**(ERROR)**"
       END IF

       DO ipatch=1,npatches_actual
          !
          !print*,'=======================VEGET =======================',ipatch
          !
          CALL vfirec(nfvegtp,VegType(1,1,ipatch),iMax*jMax,'LIN') !veg dominante no patch
          !
          !print*,'=======================FRACA =======================',ipatch
          !
          CALL vfirec(nfvegtp,FracOcc(1,1,ipatch),iMax*jMax,'LIN') !fracao ocupada pelo patch

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(VegType(:,:,ipatch) ,veg_type(:,:,ipatch))
          ELSE
             CALL IJtoIBJB(VegType(:,:,ipatch) ,veg_type(:,:,ipatch) )
          END IF

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch) )
          ELSE
             CALL IJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch)  )
          END IF

       END DO

       CLOSE(UNIT=nfvegtp)

    END IF
    !
    ! fractional area normalization
    !
    DO j=1,jbMax
       DO i=1,ibMaxPerJB(j)
          IF(frac_occ(i,j,1) < 0.99999_r8) THEN
             fractx=0.0_r8

             DO ipatch=1,npatches_actual-1
                fractx=fractx+frac_occ(i,j,ipatch)
             END DO

             frac_occ(i,j,npatches_actual)= 1.0_r8 - fractx
          END IF

       END DO
    END DO

    !IF (reducedGrid) THEN
    !   CALL NearestIBJBtoIJ(frac_occ(:,:,npatches_actual),FracOcc(:,:,npatches_actual))
    !ELSE
    !   CALL IBJBtoIJ(frac_occ(:,:,npatches_actual),FracOcc(:,:,npatches_actual))
    !END IF
    !!
    !!-
    !!
    !DO ipatch=1,npatches_actual
    !   DO j=1,jbMax
    !      DO i=1,ibMaxPerJB(j)
    !         veg_type(i,j,ipatch)=REAL(INT(veg_type(i,j,ipatch)+0.1_r8),r8)
    !      END DO
    !   END DO
    ! 
    !   IF (reducedGrid) THEN
    !      CALL NearestIBJBtoIJ(veg_type(:,:,npatches_actual),VegType(:,:,npatches_actual))
    !   ELSE
    !      CALL IBJBtoIJ(veg_type(:,:,npatches_actual),VegType(:,:,npatches_actual))
    !   END IF
    !
    !END DO
    !
    !------------------------- soil moisture initialization ------------
    !
    call MsgOne('**(read_gl_sm_bc)**','Opening GL_SM file='//TRIM(fNameSoilMoist))

    IF (record_type == 'seq') THEN !sequential mode

       OPEN(UNIT=nfslmtp,FILE=TRIM(fNameSoilMoist),FORM='unformatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilMoist), ierr
          STOP "**(ERROR)**"
       END IF

       ! do k=1,nzg   ! direct order

       DO k=nzg,1,-1 ! revert reading order
          READ(UNIT=nfslmtp) ((glsm(i,j,k),i=1,iMax),j=1,jMax) ! wetness

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          ELSE
             CALL IJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          END IF

       END DO

       CLOSE(UNIT=nfslmtp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfslmtp,FILE=TRIM(fNameSoilMoist),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilMoist), ierr
          STOP "**(ERROR)**"
       END IF


       ! do k=1,nzg   ! direct order

       DO k=nzg,1,-1 ! revert reading order
          !
          !print*,'================== GLSM for k====================',k
          !
          CALL vfirec(nfslmtp,glsm(1,1,k),iMax*jMax,'LIN')

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          ELSE
             CALL IJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          END IF

       END DO

       CLOSE(UNIT=nfslmtp)

    ELSE

       call FatalError('**(read_gl_sm_bc)** unknown record type')

    END IF

    call MsgOne('**(read_gl_sm_bc)**','DONE')

    RETURN
  END SUBROUTINE read_gl_sm_bc


  SUBROUTINE sibwet_sib2 &
       (ibmax,jbmax,sinp,sinmax,wsib,ssib,mxiter,ibMaxPerJB)
    !
    !
    !     piers sellers : 29 april 1987
    !
    !
    !   input  :   sinp     = mintz-serafini or national meteoroLOGICAL
    !                         center soil moisture (mm)
    !              sinmax   = maximum value of sinp (mm)
    !              wsinp    = m-s or nmc fractional wetness
    !              ms       = 1, mintz-serafini
    !              nmc      = 1, national meteoroLOGICAL center
    !              bee      = sib : soil moisture potential factor
    !              phsat    = sib : soil potential at saturation (m)
    !              zdepth(3)= sib : depth of 3 soil layers (m)
    !              poros    = Porosidade do solo (m"3/m"3)
    !
    !   output :   wsibt    = sib : fractional wetness
    !              ssibt    = sib : soil moisture content (m)
    !              psit     = sib : soil moisture potential (m)
    !              factor   = sib : extraction factor
    !
    INTEGER, INTENT(in   ) :: ibmax
    INTEGER, INTENT(in   ) :: jbmax
    INTEGER, INTENT(in   ) :: mxiter
    REAL(KIND=r8)   , INTENT(in   ) :: sinp(ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(in   ) :: sinmax
    !

    REAL(KIND=r8)   , INTENT(inout  ) :: wsib  (ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(inout  ) :: ssib  (ibmax,jbmax)
    INTEGER, INTENT(in   ) :: ibMaxPerJB(:)

    REAL(KIND=r8) :: sm(ityp,mxiter)
    REAL(KIND=r8) :: time(ityp,mxiter)
    REAL(KIND=r8) :: fact(ityp,mxiter)

    REAL(KIND=r8), PARAMETER :: xph1(13,2) = RESHAPE( &
         (/-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8,-120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8, &
         -200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8,-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8, &
         -120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8,-200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8/), &
         (/13,2/))
    REAL(KIND=r8), PARAMETER :: xph2(13,2) = RESHAPE( &
         (/-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8, &
         -400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8,-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8, &
         -230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8,-400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8/) , &
         (/13,2/))

    REAL(KIND=r8)    :: tzdep(3)
    REAL(KIND=r8)    :: tzltm(2)
    REAL(KIND=r8)    :: sibmax(ityp)
    REAL(KIND=r8)    :: tphsat
    REAL(KIND=r8)    :: tbee
    REAL(KIND=r8)    :: tporos
    INTEGER :: imm1
    INTEGER :: imm2
    INTEGER :: is
    INTEGER :: im
    INTEGER :: imm
    INTEGER :: ivegm
    REAL(KIND=r8)    :: cover
    REAL(KIND=r8)    :: tph1
    REAL(KIND=r8)    :: tph2
    REAL(KIND=r8)    :: sref
    REAL(KIND=r8)    :: smin
    REAL(KIND=r8)    :: dssib
    REAL(KIND=r8)    :: dw
    REAL(KIND=r8)    :: times
    REAL(KIND=r8)    :: soilmo
    REAL(KIND=r8)    :: w
    REAL(KIND=r8)    :: rsoilm
    INTEGER :: iter
    INTEGER :: latmax
    INTEGER :: lonmax
    INTEGER :: lat
    INTEGER :: lon
    REAL(KIND=r8)    :: tsinp
    REAL(KIND=r8)    :: etp
    REAL(KIND=r8)    :: facmod
    REAL(KIND=r8)    :: ssibt
    REAL(KIND=r8)    :: psit
    REAL(KIND=r8)    :: factor
    REAL(KIND=r8)    :: dt
    INTEGER :: itsoil
    INTEGER :: itfac

    sm  =0.0_r8
    time=0.0_r8
    fact=0.0_r8
    ssib=0.0_r8
    wsib=0.0_r8

    lonmax=ibmax
    latmax=jbmax

    DO is = 1,ityp
       !zdepth(3)= sib : depth of 3 soil layers (m)
       tzdep (1)= zdepth_r4(is,1)
       tzdep (2)= zdepth_r4(is,2)
       tzdep (3)= zdepth_r4(is,3)
       tphsat   = phsat_r4 (is)
       tbee     = bee_r4   (is)
       tporos   = poros_r4 (is)
       imm1=1
       imm2=1
       tzltm(1)=zlt_r4(is,1,1)
       tzltm(2)=zlt_r4(is,1,2)
       DO im=2,12
          IF(tzltm(1).LE.zlt_r4(is,im,1) ) THEN
             imm1=im
             tzltm(1)=zlt_r4(is,im,1)
          END IF
          IF(tzltm(2).LE.zlt_r4(is,im,2) )THEN
             imm2=im
             tzltm(2)=zlt_r4(is,im,2)
          END IF
       END DO
       imm=imm1
       ivegm=1
       IF(tzltm(1).LE.tzltm(2)) THEN
          imm=imm2
          ivegm=2
       END IF
       !
       !     xcover......Fracao de cobertura vegetal icg=1 topo
       !     xcover......Fracao de cobertura vegetal icg=2 base
       !
       cover=xcover_r4   (is,imm,ivegm)
       tph1=xph1         (is,ivegm)
       tph2=xph2         (is,ivegm)
       !
       !                                                     m^3
       ! sibmax(is) =(Z1 + Z2 + Z3) * poros = [m + m + m] * ----- = m = Os
       !                                                     m^3
       !
       sibmax(is) = ( tzdep(1) + tzdep(2) + tzdep(3) ) * tporos
       !
       IF(nfctrl(83).GE.1)WRITE(UNIT=nfprt,FMT=999)is,sibmax(is),tzdep(1), &
            tzdep(2),tzdep(3),tporos
       !
       !            bee      = soil moisture potential factor
       !            phsat    = soil potential at saturation   (m)
       !
       !                   --              --
       !                  | log ( - tphsat/1)|
       !  O  = Os * EXP * | -----------------|
       !                  |        b         |
       !                   --              --
       !
       sref = sibmax(is) * EXP( LOG(tphsat /(-1.0e0_r8)) /tbee)
       !                   --                          --
       !                  | log ( - tphsat/(-1.0e10) )   |
       !Omin = Os * EXP * | -----------------------------|
       !                  |              b               |
       !                   --                          --
       !
       smin    = sibmax(is) * EXP( LOG(tphsat /(-1.0e10_r8)) / tbee)
       !
       !             O - Omin
       !dssib  = ------------------
       !              mxiter
       !
       dssib   = (sref - smin) / REAL(mxiter,r8)
       !
       !              O - Omin
       ! dw    =  ------------------
       !             mxiter*Os
       !
       dw      = dssib / sibmax(is)
       !
       times   = 0.0e0_r8
       soilmo  = sref
       !
       !       O
       ! w = -----
       !       Os
       !
       w = soilmo / sibmax(is)
       !
       !                      --             --
       !                     |       0.0027    |
       !rsoilm  = 101840.0 * |1.0 - w          |
       !                     |                 |
       !                      --             --
       !
       rsoilm  = 101840.0_r8 * (1.0_r8 - w**0.0027_r8)
       DO iter = 1, mxiter
          CALL extrak_sib2( w   ,dw  ,tbee,tphsat, rsoilm, cover, &
               tph1,tph2,psit,factor )
          !
          !       dssib
          !dt = ----------
          !       factor
          !
          dt            = dssib  / factor
          !
          soilmo        = soilmo - dssib
          !
          !       O
          ! w = -----
          !       Os
          !
          w             = soilmo / sibmax(is)
          times         = times  + dt
          sm  (is,iter) = soilmo
          time(is,iter) = times
          fact(is,iter) = factor
       END DO

    END DO
    !
    !     input soil moisture map is now transformed to sib fields.
    !
    DO lat = 1, latmax
       DO lon = 1, ibMaxPerJB(lat)
          is=INT(iMaskSiB2(lon,lat),kind=i4)
          IF(is.NE.0)THEN
             tsinp = sinp(lon,lat)
             tsinp = MAX (sinmax/100.0e3_r8 , tsinp )
             tsinp = MIN (sinmax,tsinp)
             IF (tsinp .GT. 0.75e0_r8*sinmax ) etp = sinmax - tsinp
             facmod=MIN(1.0e0_r8,tsinp/(0.75e0_r8*sinmax) )
             IF (tsinp .LE. 0.75e0_r8*sinmax ) THEN
                etp = 0.75e0_r8*sinmax*LOG(0.75e0_r8*sinmax/tsinp ) + 0.25e0_r8*sinmax
             END IF
             etp = etp / 1000.0e0_r8
             DO iter = 1, mxiter
                itsoil=iter
                IF ( time(is,iter) - etp .GT. 0.0e0_r8  ) EXIT
             END DO
             DO iter=1,mxiter
                itfac=iter
                IF( fact(is,iter)-facmod-0.01e0_r8.LT.0.0e0_r8)EXIT
             END DO
             ssibt=MIN(sm(is,itsoil),sm(is,itfac))
             DO iter=1,mxiter
                IF(ssibt.GT.sm(is,iter))EXIT
             END DO
             ssib(lon,lat) = sm(is,iter)
             !
             !          O
             ! wsib = -----
             !         Os
             !
             wsib(lon,lat) = (sm(is,iter) / sibmax(is))
          END IF
       END DO
    END DO
999 FORMAT(' IS,MAX,D1,D2,D3,POROS=',I2,1X,5E12.5)
  END SUBROUTINE sibwet_sib2
  
  
  SUBROUTINE extrak_sib2( w, dw, tbee, tphsat, rsoilm, cover, tph1, tph2, &
       psit, factor )
    REAL(KIND=r8), INTENT(in   ) :: w
    REAL(KIND=r8), INTENT(in   ) :: dw
    REAL(KIND=r8), INTENT(in   ) :: tbee
    REAL(KIND=r8), INTENT(in   ) :: tphsat
    REAL(KIND=r8), INTENT(in   ) :: rsoilm
    REAL(KIND=r8), INTENT(in   ) :: cover
    REAL(KIND=r8), INTENT(in   ) :: tph1
    REAL(KIND=r8), INTENT(in   ) :: tph2
    REAL(KIND=r8), INTENT(inout  ) :: psit
    REAL(KIND=r8), INTENT(inout  ) :: factor
    REAL(KIND=r8) :: rsoil
    REAL(KIND=r8) :: argg
    REAL(KIND=r8) :: hr
    REAL(KIND=r8) :: rplant
    !                --     -- (-b)
    !               |      dw |                  0
    ! psit = PHYs * | w - --- |      where w = -----
    !               |      2  |                  0s
    !                --     --
    psit   = tphsat * ( w-dw/2.0e0_r8 ) ** (-tbee)
    !
    !                      --                        --
    !                     |       --     -- (0.0027)   |
    !                     |      |      dw |           |
    !rsoil   = 101840.0 * |1.0 - | w - --- |           |
    !                     |      |      2  |           |
    !                     |       --     --            |
    !                      --                        --
    !
    rsoil  = 101840.0_r8 * (1.0_r8-( w-dw/2.0_r8) ** 0.0027_r8)
    !
    !                9.81       1
    !argg = psit * -------- * -------
    !               461.50     310.0
    !
    argg   = MAX ( -10.0e0_r8 , ((psit * 9.81e0_r8 / 461.5e0_r8) / 310.e0_r8))
    !
    !            --                       --
    !           |         9.81       1      |
    !hr   = EXP |psit * -------- * -------  |
    !           |        461.50     310.0   |
    !            --                       --
    !
    hr     = EXP ( argg )
    !
    !         rsoilm
    ! rsoil =--------- * hr
    !         rsoil
    !
    rsoil  = rsoilm /rsoil * hr
    !
    !          ( psit - tph2 - 50.0)
    !rplant = -------------------------
    !             ( tph1 - tph2 )
    !
    rplant = ( psit - tph2 -50.0_r8) / ( tph1 - tph2 )
    rplant = MAX ( 0.0e0_r8, MIN ( 1.0e0_r8, rplant ) )
    !                                                                     --                   --
    !                  --                 --                             |     --     -- (0.0027)|
    !                 |( psit - tph2 - 50)  |                            |    |      dw |        |
    !factor = cover * |---------------------| + (1 - cover) * 101840.0 * |1 - | w - --- |        |
    !                 |   ( tph1 - tph2 )   |                            |    |      2  |        |
    !                  --                 --                             |     --     --         |
    !                                                                     --                   --
    factor = cover * rplant + ( 1.0e0_r8 - cover ) * rsoil
    factor = MAX ( 1.e-6_r8, factor )
  END SUBROUTINE extrak_sib2

  SUBROUTINE Phenology_sib2(latco,nCols,idatec,ndvi,ndvim,colrad,itype)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: latco
    INTEGER, INTENT(IN   ) :: idatec(4)
    REAL(KIND=r8),    INTENT(INOUT) :: ndvi  (ncols) 
    REAL(KIND=r8),    INTENT(INOUT) :: ndvim  (ncols) 
    REAL(KIND=r8),    INTENT(IN   ) :: colrad(ncols) 
    !INTEGER, INTENT(IN   ) :: mon  (ncols) 
    INTEGER, INTENT(IN   ) :: itype(ncols) 
    INTEGER       :: ncount
    !INTEGER       :: imes
    !INTEGER       :: irec
    !INTEGER       :: j
    INTEGER       :: i
    INTEGER       :: k
    INTEGER       :: kk    
    INTEGER       :: doy,nsib
    REAL(KIND=r8) ::   pndvi_sib2 (ncols)
    REAL(KIND=r8) ::   cndvi_sib2 (ncols)
    REAL(KIND=r8) ::   latitude    (ncols)
    ncount=0
    DO i=1,nCols       
       IF(iMaskSiB2(i,latco)>=1_i8)THEN
          ncount=ncount+1
          IF(ndvi (i) > 1.0_r8)WRITE(*,*)ndvi (i)
          IF(ndvi (i) < 0.0_r8)WRITE(*,*)ndvi (i)
              pndvi_sib2           (ncount)            =  MAX(MIN(ndvim(i),0.99_r8),0.001_r8)
              cndvi_sib2           (ncount)            =  MAX(MIN(ndvi (i),0.99_r8),0.001_r8)
              IF(iMaskSiB2(i,latco)>=13_i8)pndvi_sib2(ncount)=  0.000_r8
              IF(iMaskSiB2(i,latco)>=11_i8)pndvi_sib2(ncount)=  0.000_r8
              IF(iMaskSiB2(i,latco)>=13_i8)cndvi_sib2(ncount)=  0.000_r8
              IF(iMaskSiB2(i,latco)>=11_i8)cndvi_sib2(ncount)=  0.000_r8
          !
          ! old ndvi conversion: NDVI(k)=(real(value)-511.)/512.
          !
          !  NDVI(k)=(value+NDVIoffset)/NDVIscale

              latitude           (ncount)            =  90.0_r8-180.0e0_r8/pie * colrad(ncount)!lati(jbMax+1-j)
       END IF
    END DO
    nsib=ncount
    DO i=1,nsib
       !            old ndvi conversion: NDVI(k)=(real(value)-511.)/512.
       IF(Flagsib%fVcov.EQ.1) THEN
             CALL FractionVegCover( cndvi_sib2(i),fPARmax,fPARmin,&
                  MorphTab(itype(i) )%SRmax,&
                  MorphTab(itype(i) )%SRmin,BioVar%fVCover)
             vcoverg(i,latco)=   BioVar%fVCover
       ELSEIF(Flagsib%fVcov.EQ.2) THEN
          WRITE(*,*)'read(16,rec=RecNum) BioVar%fVCover' 
       ELSEIF(Flagsib%fVcov.EQ.3) THEN
             vcoverg(i,latco)=   BioTab(itype(i))%fVCover!vcoverg(i,j)
       ENDIF
    END DO

    doy=julday(idatec(3),idatec(2),idatec(4))

    DO i=1,nsib    
              morphtab_zc_sib2     (i,latco) = MorphTab(itype(i))%zc   
              morphtab_lwidth_sib2 (i,latco) = MorphTab(itype(i))%LWidth 
              morphtab_llength_sib2(i,latco) = MorphTab(itype(i))%LLength
              morphtab_laimax_sib2 (i,latco) = MorphTab(itype(i))%LAImax 
              morphtab_stems_sib2  (i,latco) = MorphTab(itype(i))%stems  
              morphtab_ndvimax_sib2(i,latco) = MorphTab(itype(i))%NDVImax
              morphtab_ndvimin_sib2(i,latco) = MorphTab(itype(i))%NDVImin
              morphtab_srmax_sib2  (i,latco) = MorphTab(itype(i))%SRmax  
              morphtab_srmin_sib2  (i,latco) = MorphTab(itype(i))%SRmin
              DO k=1,50
                 DO kk=1,50
                    aerovar_zo_sib2      (i,k,kk,latco) = AeroVar(itype(i),k,kk)%zo
                    aerovar_zp_disp_sib2 (i,k,kk,latco) = AeroVar(itype(i),k,kk)%zp_disp
                    aerovar_rbc_sib2     (i,k,kk,latco) = AeroVar(itype(i),k,kk)%RbC  
                    aerovar_rdc_sib2     (i,k,kk,latco) = AeroVar(itype(i),k,kk)%RdC
                 END DO
              END DO
              chil_sib2 (i,latco) = BioTab(itype(i))%ChiL 
              DO k=1,2
                 DO kk=1,2
                    tran_sib2(i,k,kk,latco) = BioTab(itype(i))%LTran(k,kk)
                    ref_sib2 (i,k,kk,latco) = BioTab(itype(i))%LRef (k,kk) 
                 END DO
              END DO
              vcover2g           (i,latco)            = vcoverg(i,latco)! BioTab(itype(i))%fVCover!vcoverg(i,j)
    END DO
          !          Calculate time dependant variables

          CALL mapper(                                       &
               latitude             (1:nsib) , &  !IN
               doy                                  , &  !IN
               pndvi_sib2           (1:nsib) , &  !IN
               cndvi_sib2           (1:nsib) , &  !IN
               vcover2g             (1:nsib,latco) , &  !IN
               chil_sib2            (1:nsib,latco) , &  !IN
               tran_sib2            (1:nsib,1:2,1:2,latco) , &  !IN
               ref_sib2             (1:nsib,1:2,1:2,latco) , &  !IN
               !morphtab_zc_sib2     (1:nsib,latco) , &  !IN 
               !morphtab_lwidth_sib2 (1:nsib,latco) , &  !IN 
               !morphtab_llength_sib2(1:nsib,latco) , &  !IN 
               morphtab_laimax_sib2 (1:nsib,latco) , &  !IN 
               morphtab_stems_sib2  (1:nsib,latco) , &  !IN 
               morphtab_ndvimax_sib2(1:nsib,latco) , &  !IN 
               morphtab_ndvimin_sib2(1:nsib,latco) , &  !IN 
               morphtab_srmax_sib2  (1:nsib,latco) , &  !IN 
               morphtab_srmin_sib2  (1:nsib,latco) , &  !IN 
               aerovar_zo_sib2            (1:nsib,1:50,1:50,latco) , &  !IN 
               aerovar_zp_disp_sib2 (1:nsib,1:50,1:50,latco) , &  !IN 
               aerovar_rbc_sib2            (1:nsib,1:50,1:50,latco) , &  !IN 
               aerovar_rdc_sib2            (1:nsib,1:50,1:50,latco) , &  !IN 
               LAIgrid              (              1:50) , &  !IN
               fVCovergrid          (              1:50) , &  !IN
               timevar_fpar            (1:nsib,latco) , &  !INOUT
               timevar_lai             (1:nsib,latco) , &  !INOUT zlt_r4  (ityp,imon,icg)
               timevar_green            (1:nsib,latco) , &  !INOUT green_r4(ityp,imon,icg)
               timevar_zo              (1:nsib,latco) , &  !OUT   x0x_r4  (ityp,imon)
               timevar_zp_disp            (1:nsib,latco) , &  !OUT   xd_r4   (ityp,imon)
               timevar_rbc             (1:nsib,latco) , &  !OUT   xbc_r4  (ityp,imon)
               timevar_rdc             (1:nsib,latco) , &  !OUT   xdc_r4  (ityp,imon)
               timevar_gmudmu            (1:nsib,latco) , &  !INOUT
               nsib                                )  
  END SUBROUTINE Phenology_sib2




 SUBROUTINE ReadSurfaceMaskSib2(ibMax,jbMax,iMax,jMax,fNameSiB2Mask,fNameSandMask,&
            fNameClayMask,fNameTextMask)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: iMax
    INTEGER, INTENT(IN   ) :: jMax  
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMax  
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSiB2Mask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSandMask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameClayMask
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameTextMask

    REAL(KIND=r4) :: SoilNum
    INTEGER (KIND=i8)      :: imask_in (iMax,jMax)
    INTEGER(KIND=i8)       :: mskant_in(iMax,jMax)
    INTEGER       :: ier(iMax,jMax)
    !REAL(KIND=r4) ::   brf (iMax,jMax)
   ! INTEGER       ::   bif (iMax,jMax)
    INTEGER(KIND=i4) :: iBioNum
    REAL(KIND=r8) :: TextClay  (iMax,jMax)        
    REAL(KIND=r8) :: TextSand  (iMax,jMax)                
    TYPE(soil_Physical) :: SoilGrd(iMax,jMax)

    INTEGER       :: irec,ierr,LRecIN
    INTEGER       :: j
    INTEGER       :: i
    imask_in =0_i8     
    !
    ! open soil reflectivity maps not used
    !IF(Flagsib%mode.EQ.1) THEN
    !   IF(Flagsib%SoRefTab.EQ.2) THEN
    !      OPEN(12, file=TRIM(FileName%SoRefVis), form='unformatted', &
    !            access='direct',recl=1, status='old', action='read')
    !      OPEN(13, file=TRIM(FileName%SoRefNIR), form='unformatted', &
    !            access='direct',recl=1, status='old', action='read')
    !   ENDIF
    !ENDIF
    !
    ! open fractional vegetation cover map not used
    !IF(Flagsib%fVCov.EQ.2) THEN
    !   OPEN(16, file=TRIM(FileName%fVCovMap), form='unformatted', &
    !            access='direct',recl=1, status='old', action='read')
    !ENDIF
    !
    !--------------------------------------------------------------------------
    ! Execute Mapper on grid specified in map.in
    !--------------------------------------------------------------------------
    !IF (Flagsib%Map.EQ.1.AND.Flagsib%mode.EQ.1) THEN
       !
       ! open vegetation/biome type map
       iBioNum=0
       INQUIRE (IOLENGTH=LRecIN) iBioNum
       OPEN(nfsib2mk, file=TRIM(fNameSiB2Mask), form='unformatted', &
            access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSiB2Mask), ierr
          STOP "**(ERROR)**"
       END IF
 !      READ (UNIT=nfsib2mk, REC=1) brf

       irec=0
       DO j=1,jMax
          DO i=1,iMax
             irec=irec+1
             READ(nfsib2mk,rec=irec) iBioNum
             !IF(BioNum2==13.0 .or.BioNum2==11.0 )BioNum2=6.0
             !IF(BioNum2.GE.REAL(minBiome,kind=4).AND.BioNum2.LE.REAL(maxBiome,KIND=4)) THEN
             !                 Map lookup table values onto grid cell values
                     imask_in(i,j)=INT(iBioNum,kind=i8)
             !END IF
          ENDDO
       ENDDO
       CLOSE(nfsib2mk)

       !DO j=1,jMax
       !   DO i=1,iMax
!                imask_in  (i,j) = INT(brf (i,j),kind=i8)
                !iMaskSiB2 (i,j) = imask_in(i,j)
               ! WRITE(*,*)'imask_in(i,j)',i,j,imask_in(i,j),iMax,jMAx,ibMax,jbMax
!          ENDDO
!       ENDDO

       IF (reducedGrid) THEN
          CALL FreqBoxIJtoIBJB(imask_in,iMaskSiB2)
       ELSE
          CALL IJtoIBJB( imask_in,iMaskSiB2)
       END IF
       imask=iMaskSiB2
       DO j=1,jMax
          DO i=1,iMax            
             IF (imask_in(i,j) >= 1_i8) THEN
                 ier(i,j) = 0
             ELSE
                  ier(i,j) = 1
             END IF 
          END DO
          IF (ANY( ier(1:iMax,j) /= 0)) THEN
             DO i=1,iMax          
                mskant_in(i,j) = 1_i8
             END DO        
          ELSE
             DO i=1,iMax          
                mskant_in(i,j) = 0_i8
             END DO                     
          END IF          
       END DO
       IF (reducedGrid) THEN
          CALL FreqBoxIJtoIBJB(mskant_in,MskAntSib2)
       ELSE
          CALL IJtoIBJB( mskant_in,MskAntSib2)
       END IF
       MskAnt=MskAntSib2
       ! set soil characteristics

       IF(Flagsib%SoilMap.EQ.1) THEN ! read %sand/clay from map
          INQUIRE (IOLENGTH=LRecIN) iBioNum
          OPEN(nfclay, file=TRIM(fNameClayMask), form='unformatted', &
               access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameClayMask), ierr
             STOP "**(ERROR)**"
          END IF
          
          INQUIRE (IOLENGTH=LRecIN) iBioNum
          OPEN(nfsand, file=TRIM(fNameSandMask), form='unformatted', &
               access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSandMask), ierr
             STOP "**(ERROR)**"
          END IF

          irec=0
          DO j=1,jMax
             DO i=1,iMax
                irec=irec+1

                READ(nfclay,rec=irec) iBioNum
                TextClay(i,j) = iBioNum
                READ(nfsand,rec=irec) iBioNum
                TextSand(i,j) = iBioNum
                IF     (Flagsib%SoilProp.EQ.1) THEN ! calc soil prop. from %sand/clay
                   CALL SoilProperties(TextClay(i,j),TextSand(i,j))                   
                   SoilGrd(i,j)%SoilNum =SoilVar%SoilNum
                   SoilGrd(i,j)%BEE     =SoilVar%BEE        
                   SoilGrd(i,j)%PhiSat  =SoilVar%PhiSat 
                   SoilGrd(i,j)%SatCo   =SoilVar%SatCo  
                   SoilGrd(i,j)%poros   =SoilVar%poros  
                   SoilGrd(i,j)%Slope   =SoilVar%Slope  
                   SoilGrd(i,j)%Wopt    =SoilVar%Wopt        
                   SoilGrd(i,j)%Skew    =SoilVar%Skew        
                   SoilGrd(i,j)%RespSat =SoilVar%RespSat
                   CALL textclass(TextClay(i,j),TextSand(i,j))
                   SoilNum              =Text%Class
                   SoilGrd(i,j)%SoilNum=SoilNum
                ELSE IF(Flagsib%SoilProp.EQ.2) THEN ! calc class & use soil table
                   CALL textclass(TextClay(i,j),TextSand(i,j))
                   IF(imask_in(i,j) >= 1_i8)THEN
                     ! IF(imask_in(i,j)==13 .AND.Text%Class ==  0) Text%Class=6
                      SoilNum              =Text%Class
                      SoilVar              =SoilTab(INT(SoilNum))
                      SoilGrd(i,j)%SoilNum =SoilVar%SoilNum
                      SoilGrd(i,j)%BEE     =SoilVar%BEE        
                      SoilGrd(i,j)%PhiSat  =SoilVar%PhiSat 
                      SoilGrd(i,j)%SatCo   =SoilVar%SatCo   
                      SoilGrd(i,j)%poros   =SoilVar%poros  
                      SoilGrd(i,j)%Slope   =SoilVar%Slope  
                      SoilGrd(i,j)%Wopt           =SoilVar%Wopt        
                      SoilGrd(i,j)%Skew           =SoilVar%Skew        
                      SoilGrd(i,j)%RespSat =SoilVar%RespSat
                   END IF
                ENDIF
             END DO
          END DO
         ! DO j=1,jbMax
         !    DO i=1,ibMax
         !       bSoilGrd(i,j)%SoilNum=SoilGrd(i,j)%SoilNum
        !        bSoilGrd(i,j)%BEE    =SoilGrd(i,j)%BEE         
        !        bSoilGrd(i,j)%PhiSat =SoilGrd(i,j)%PhiSat 
        !        bSoilGrd(i,j)%SatCo  =SoilGrd(i,j)%SatCo  
        !        bSoilGrd(i,j)%poros  =SoilGrd(i,j)%poros  
        !        bSoilGrd(i,j)%Slope  =SoilGrd(i,j)%Slope  
        !        bSoilGrd(i,j)%Wopt   =SoilGrd(i,j)%Wopt   
        !        bSoilGrd(i,j)%Skew   =SoilGrd(i,j)%Skew   
        !        bSoilGrd(i,j)%RespSat=SoilGrd(i,j)%RespSat
        !     END DO
        !  END DO
          IF (reducedGrid) THEN
             CALL FreqBoxIJtoIBJB(TextSand,sandfrac)
             sandfrac=sandfrac/100.0_r8   
          ELSE
             CALL IJtoIBJB( TextSand/100.0_r8,sandfrac)
          END IF
          IF (reducedGrid) THEN
             CALL FreqBoxIJtoIBJB(TextClay,clayfrac)
             clayfrac =clayfrac/100.0_r8    
          ELSE
             CALL IJtoIBJB( TextClay/100.0_r8,clayfrac)
          END IF

          IF (reducedGrid) THEN
             CALL FreqBoxIJtoIBJB(SoilGrd%SoilNum,bSoilGrd%SoilNum)
          ELSE
             CALL IJtoIBJB( SoilGrd%SoilNum,bSoilGrd%SoilNum)
          END IF
          
          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%BEE,bSoilGrd%BEE)
          ELSE
             CALL IJtoIBJB( SoilGrd%BEE,bSoilGrd%BEE)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%PhiSat,bSoilGrd%PhiSat)
          ELSE
             CALL IJtoIBJB( SoilGrd%PhiSat,bSoilGrd%PhiSat)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%SatCo,bSoilGrd%SatCo)
          ELSE
             CALL IJtoIBJB( SoilGrd%SatCo,bSoilGrd%SatCo)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%poros,bSoilGrd%poros)
          ELSE
             CALL IJtoIBJB( SoilGrd%poros,bSoilGrd%poros)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%Slope,bSoilGrd%Slope)
          ELSE
             CALL IJtoIBJB( SoilGrd%Slope,bSoilGrd%Slope)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%Wopt,bSoilGrd%Wopt)
          ELSE
             CALL IJtoIBJB( SoilGrd%Wopt,bSoilGrd%Wopt)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%Skew,bSoilGrd%Skew)
          ELSE
             CALL IJtoIBJB( SoilGrd%Skew,bSoilGrd%Skew)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%RespSat,bSoilGrd%RespSat)
          ELSE
             CALL IJtoIBJB( SoilGrd%RespSat,bSoilGrd%RespSat)
          END IF

         CLOSE(nfclay)
         CLOSE(nfsand)
       ELSE IF(Flagsib%SoilMap.EQ.2) THEN ! read class from map, use soil table
          INQUIRE (IOLENGTH=LRecIN) iBioNum
          OPEN(nftext, file=TRIM(fNameTextMask), form='unformatted', &
               access='direct',recl=LRecIN, status='old', action='read', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameTextMask), ierr
             STOP "**(ERROR)**"
          END IF

          irec=0
          DO j=1,jMax
             DO i=1,iMax
                     irec=irec+1
                READ(nftext,rec=irec) iBioNum
                IF(iBioNum==6.0)SoilNum=13.0
                SoilVar=SoilTab(INT(iBioNum))
                SoilGrd(i,j)%SoilNum =SoilVar%SoilNum
                SoilGrd(i,j)%BEE        =SoilVar%BEE        
                SoilGrd(i,j)%PhiSat  =SoilVar%PhiSat 
                SoilGrd(i,j)%SatCo        =SoilVar%SatCo  
                SoilGrd(i,j)%poros        =SoilVar%poros  
                SoilGrd(i,j)%Slope        =SoilVar%Slope  
                SoilGrd(i,j)%Wopt        =SoilVar%Wopt        
                SoilGrd(i,j)%Skew        =SoilVar%Skew        
                SoilGrd(i,j)%RespSat =SoilVar%RespSat
             END DO
          END DO
        !  DO j=1,jbMax
        !     DO i=1,ibMax
        !        bSoilGrd(i,j)%SoilNum=SoilGrd(i,j)%SoilNum
        !        bSoilGrd(i,j)%BEE    =SoilGrd(i,j)%BEE         
        !        bSoilGrd(i,j)%PhiSat =SoilGrd(i,j)%PhiSat 
        !        bSoilGrd(i,j)%SatCo  =SoilGrd(i,j)%SatCo  
        !        bSoilGrd(i,j)%poros  =SoilGrd(i,j)%poros  
        !        bSoilGrd(i,j)%Slope  =SoilGrd(i,j)%Slope  
        !        bSoilGrd(i,j)%Wopt   =SoilGrd(i,j)%Wopt   
        !        bSoilGrd(i,j)%Skew   =SoilGrd(i,j)%Skew   
        !        bSoilGrd(i,j)%RespSat=SoilGrd(i,j)%RespSat
          !   END DO
          !END DO

          IF (reducedGrid) THEN
             CALL FreqBoxIJtoIBJB(SoilGrd%SoilNum,bSoilGrd%SoilNum)
          ELSE
             CALL IJtoIBJB( SoilGrd%SoilNum,bSoilGrd%SoilNum)
          END IF
          
          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%BEE,bSoilGrd%BEE)
          ELSE
             CALL IJtoIBJB( SoilGrd%BEE,bSoilGrd%BEE)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%PhiSat,bSoilGrd%PhiSat)
          ELSE
             CALL IJtoIBJB( SoilGrd%PhiSat,bSoilGrd%PhiSat)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%SatCo,bSoilGrd%SatCo)
          ELSE
             CALL IJtoIBJB( SoilGrd%SatCo,bSoilGrd%SatCo)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%poros,bSoilGrd%poros)
          ELSE
             CALL IJtoIBJB( SoilGrd%poros,bSoilGrd%poros)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%Slope,bSoilGrd%Slope)
          ELSE
             CALL IJtoIBJB( SoilGrd%Slope,bSoilGrd%Slope)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%Wopt,bSoilGrd%Wopt)
          ELSE
             CALL IJtoIBJB( SoilGrd%Wopt,bSoilGrd%Wopt)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%Skew,bSoilGrd%Skew)
          ELSE
             CALL IJtoIBJB( SoilGrd%Skew,bSoilGrd%Skew)
          END IF

          IF (reducedGrid) THEN
             CALL LinearIJtoIBJB(SoilGrd%RespSat,bSoilGrd%RespSat)
          ELSE
             CALL IJtoIBJB( SoilGrd%RespSat,bSoilGrd%RespSat)
          END IF          
          CLOSE(nftext)
       ENDIF
       DO j=1,jbMax
          DO i=1,ibMax
             SoilMask(i,j)=bSoilGrd(i,j)%SoilNum
          END DO
       END DO

  END SUBROUTINE ReadSurfaceMaskSib2

  SUBROUTINE vegin_sib2(fNameSoilTab,fNameMorfTab,fNameBioCTab,fNameAeroTab)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSoilTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMorfTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameBioCTab
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameAeroTab
    INTEGER :: ierr
    ALLOCATE(rstpar_r4(ityp,icg,iwv))
    ALLOCATE(chil_r4  (ityp,icg))
    ALLOCATE(topt_r4  (ityp,icg))
    ALLOCATE(tll_r4   (ityp,icg))
    ALLOCATE(tu_r4    (ityp,icg))
    ALLOCATE(defac_r4 (ityp,icg))
    ALLOCATE(ph1_r4   (ityp,icg))
    ALLOCATE(ph2_r4   (ityp,icg))
    ALLOCATE(rootd_r4 (ityp,icg))
    ALLOCATE(bee_r4   (ityp))
    ALLOCATE(phsat_r4 (ityp))
    ALLOCATE(satco_r4 (ityp))
    ALLOCATE(poros_r4 (ityp))
    ALLOCATE(zdepth_r4(ityp,idp))
    ALLOCATE(green_r4 (ityp,imon,icg))
    ALLOCATE(xcover_r4(ityp,imon,icg))
    ALLOCATE(zlt_r4   (ityp,imon,icg))
    ALLOCATE(x0x_r4   (ityp,imon))
    ALLOCATE(xd_r4    (ityp,imon))
    ALLOCATE(z2_r4    (ityp,imon))
    ALLOCATE(z1_r4    (ityp,imon))
    ALLOCATE(xdc_r4   (ityp,imon))
    ALLOCATE(xbc_r4   (ityp,imon))
    OPEN(UNIT=nfsibd, FILE=TRIM(fNameSibVeg),FORM='UNFORMATTED', ACCESS='SEQUENTIAL',&
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSibVeg), ierr
       STOP "**(ERROR)**"
    END IF
    READ (UNIT=nfsibd) rstpar_r4, chil_r4, topt_r4, tll_r4, tu_r4, defac_r4, ph1_r4, ph2_r4, &
                       rootd_r4, bee_r4, phsat_r4, satco_r4, poros_r4, zdepth_r4
    READ (UNIT=nfsibd) green_r4, xcover_r4, zlt_r4, x0x_r4, xd_r4, z2_r4, z1_r4, xdc_r4, xbc_r4

    !
    !--------------------------------------------------------------------------
    ! input control variables for various wrapper options
    !--------------------------------------------------------------------------
    !
    CALL ReadControl
    !
    ! allocate variables based on user input from file map.in
    !ALLOCATE(TimeVar ))
    ALLOCATE(SoilTab (ns))
    ALLOCATE(MorphTab(nv))
    ALLOCATE(BioTab  (nv))
    ALLOCATE(AeroVar (nv,50,50))
    !
    !--------------------------------------------------------------------------
    ! Read all lookup tables, open input and output files
    !--------------------------------------------------------------------------
    IF(Flagsib%Print.EQ.1) PRINT *, 'Reading Lookup tables'
    use100=0
    !PRINT *,'wrapper use100=',use100
    !
    ! read Century generated LAI's
    !      call ReadCenturyTab    
    !
    ! read table containing soil types 
    !
    CALL ReadSoilTable(fNameSoilTab)
    !PRINT *,'1'
    !
    !
    ! read table containing morphology characteristics for each biome type
    CALL ReadMorphTable(fNameMorfTab)
    !PRINT *,'2'
    !
    ! read table of biome dependant variables
    CALL ReadBioTable(fNameBioCTab)
    !PRINT *,'3'

    !
    ! read in interpolation tables of aerodynamic variables
    CALL ReadAeroTables(fNameAeroTab)
    !PRINT *,'4'
    !
    !

  END SUBROUTINE vegin_sib2

   !
  !=======================================================================
  SUBROUTINE ReadControl
    !=======================================================================
    ! This subroutine reads in all the control variables, flags, 
    ! and input files for various wrapper options.
    !
    !
    IMPLICIT NONE
    !
    !
    ! Read input flags                       start input in column 46--------------------|
    Flagsib%Map     =1! Run mapper?                     (1=yes 2=no)=1
    Flagsib%Mode    =1! Mapper grid or single point(1=grid 2=single)=2
    Flagsib%EzPlot  =2! Make Ezplot file?             (1=yes 2=no)=2
    Flagsib%Sib     =2! Make Sib_bc file?             (1=yes 2=no)=2
    Flagsib%single  =1! Make single pt Sib_bc file?    (1=yes 2=no)=1
    Flagsib%Stats   =2! Make stats file?             (1=yes 2=no)=2
    Flagsib%Monthly =2! Monthly output files?             (1=yes 2=no)=2
    Flagsib%GridKey =2! lat/lon key for sib points?    (1=yes 2=no)=2
    Flagsib%Print   =1! Print messages?              (1=yes 2=no)=1
    Flagsib%SoRefTab=1! Soil Reflectances?          (1=Table 2=Map)=1
    Flagsib%SoilMap =1! Soil input/Map?     (1=%clay/sand 2=class)=1
    Flagsib%SoilProp=1! Soil properties?    (1=%clay/sand 2=table)=1
    Flagsib%fVCov   =3! Fraction Veg. Cover?  (1=calc 2=Map 3=table)=1
    !                                     --------------------------------------------|
    ! Read header line
    !
    ! Read input variables                    --------------------------------------------|
    !      year        =98   ! Calculation Year (last 2 digits)           =98
    nm        =12   ! Number of months in a year                   =12
    ndays        =1  ! Number of days in a year                   =365
    !      DOYstart  =15.2 ! Day of Year (DOY) first NDVI input map      =15.2
    !      DDoy        =30.4 ! Time between NDVI input maps (day)           =30.4
    NDVIoffset=0.0_r8   ! ndvi offset (ndvi=(value+offset)/scale)     =0.
    NDVIscale =1.0_r8! ndvi scale factor(ndvi=(value+offset)/scale)=1000.
    ns        =13   ! Total number of soil types                   =12
    nv        =13   ! Total number of biome types                   =13
    minBiome  =1    ! min biome type number (1-13)                   =1
    maxBiome  =13   ! max biome type number (1-13)                   =12
    LatMin        =-90.0_r8 ! Min Latitude (deg) lower left Domain corner =-90.
    LonMin        =-180.0_r8! Min Longitude (deg) lower left Domain corner=-180.
    Dlat        =1.   ! Domain grid spacing: Delta Latitude (deg)   =1.
    Dlon        =1.   ! Domain grid spacing: Delta Longitude (deg)  =1.
    imm        =360  ! Total Number Longitude pts for Domain            =360
    jmm        =180  ! Total Number Latitude pts for Domain           =180
    !      SubLatMin =-90. ! Min Latitude (deg) for sub-grid option      =-90.
    !      SubLatMax =90.  ! Max Latitude (deg) for sub-grid option      =90.
    !      SubLonMin =-180.! Min Longitude (deg) for sub-grid option     =-180.
    !      SubLonMax =180. ! Max Longitude (deg) for sub-grid option     =180.
    !                                    --------------------------------------------|
    ! Read input file names
    !
    FileName%BioMap        ='veg192x96.T062.bin' !Lookup map of Vegetation types    
    FileName%BioTab        ='BioChar.Tab'        !Lookup table of Veg. type char.   
    FileName%MorphTab ='Morph.Tab.orig' !Lookup table veg. morphilogical char.
    FileName%SoilTab  ='SoilChar.Tab'        !Lookup table of soil physical char.
    FileName%SoRefVis ='rhoVis.Map'        !Lookup map soil visible reflectances
    FileName%SoRefNIR ='rhoNIR.Map'        !Lookup map soil near IR reflectances
    FileName%SoilMap  ='SoilType.Map'        !Lookup map of soil types           
    FileName%PercSand ='SoilSand192x96.T062.bin'        !Lookup map of soil percent sand   
    FileName%PercClay ='SoilClay192x96.T062.bin'        !Lookup map of soil percent clay   
    FileName%fVCovMap ='fVCover.Map'        !Lookup map of fractional vegetation cover
    FileName%AeroVar  ='AeroVar.Tab'        !Interpolation tables for aero variables
    FileName%asciiNDVI ='ndvi_ascii'                     !ascii NDVI file for single point option  
    ! this file contains all the century generated lai's    
    CenturyLAI='CenturyLAI.Tab.1998.10c'!Century generated LAI's         
    !
    !                                       --------------------------------------------
    ! Read output file names                ! file names for output files (70 char max)---
    FileName%sib_bc         ='sib_bc_'     !Sib BC filename (yy attached)          
    FileName%sib_grid  ='sib_gridmap' !Sib gridmap filename                    
    FileName%sib_biom  ='sib_biome'   !GCM biome type filename            
    FileName%stats         ='StatsDat_'   !Stats data filename (yy attached)  
    FileName%Monthly   ='sib_bc__'    !monthly data filename (yy attached)
    FileName%GridKey   ='Grid_key'    !Sib point lat/lon grid key file    
    !                                                                 
    ! Read NDVI file names                        -------------------------------------------
    !      do k = 1, nm
    !           read (9,17) NDVIfiles(k)
    !      enddo
    RETURN
  END SUBROUTINE ReadControl

  !
  !======================================================================
  SUBROUTINE ReadSoilTable(fNameSoilTab)
    !======================================================================
    !     read in soil properties for each texture class
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSoilTab
    CHARACTER(LEN=100) :: junk  ! junk variable to read character text from table
    INTEGER            :: s    ! soil type index
    OPEN (nfsoiltb,file=TRIM(fNameSoilTab),form='formatted')
    READ(nfsoiltb,'(a100)') junk
    DO s=1,ns
       READ(nfsoiltb,'(i2,1x,2(f7.3,1x),e9.3, 5(f7.3,1x))') SoilTab(s)%SoilNum,&
            SoilTab(s)%BEE,    &
            SoilTab(s)%PhiSat, &
            SoilTab(s)%SatCo,  &
            SoilTab(s)%poros,  &
            SoilTab(s)%Slope,  &
            SoilTab(s)%Wopt,   &
            SoilTab(s)%Skew,   &
            SoilTab(s)%RespSat 
    ENDDO
    CLOSE(nfsoiltb)
    RETURN
  END SUBROUTINE ReadSoilTable
  !
  !
  !======================================================================
  SUBROUTINE ReadMorphTable(fNameMorfTab)
    !======================================================================
    !     read in morphological and bounding ndvi values for each vegtype
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    ! 
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMorfTab
    INTEGER :: iv       ! biome type do-loop counter
    INTEGER :: ivv      ! biome type number
    !
    OPEN(nfmorftb,file=TRIM(fNameMorfTab),form='formatted')
    !
    DO iv=1, nv
       !
       READ(nfmorftb,*) ivv, MorphTab(iv)%zc, MorphTab(iv)%LWidth, &
            MorphTab(iv)%LLength, &
            MorphTab(iv)%LAImax, &
            MorphTab(iv)%stems, &
            MorphTab(iv)%NDVImax, & 
            MorphTab(iv)%NDVImin
       !
       ! Convert maximum/minimum NDVI values to simple ratios 
       MorphTab(iv)%SRmax=(1.0_r8+MorphTab(iv)%NDVImax)/(1.0_r8-MorphTab(iv)%NDVImax)
       MorphTab(iv)%SRmin=(1.0_r8+MorphTab(iv)%NDVImin)/(1.0_r8-MorphTab(iv)%NDVImin)
       !
    ENDDO
    CLOSE(nfmorftb)
    RETURN
  END SUBROUTINE ReadMorphTable
  !=======================================================================
  SUBROUTINE ReadBioTable(fNameBioCTab)
    !=======================================================================       
    !     read in vegetation parameter lookup table                                
    !-----------------------------------------------------------------------
    !      
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameBioCTab
    CHARACTER(LEN=100) ::  junk  ! junk variable to read character text from table
    INTEGER :: iv        ! counter for biome type do-loop
    INTEGER :: ivv       ! temporary variable for reading vegetation type

    !                                                                               
    OPEN(nfbioctb,file=TRIM(fNameBioCTab),form='formatted')
    !
    READ(nfbioctb,'(a100)') junk
    DO iv=1,nv
       READ(nfbioctb,13) ivv,                    &
            BioTab(iv)%z2,            &
            BioTab(iv)%z1,            &
            BioTab(iv)%fVcover,     &
            BioTab(iv)%Chil,            &
            BioTab(iv)%SoDep,            &
            BioTab(iv)%RootD,            &
            BioTab(iv)%Phi_half
       !
       !       assign biome type to BioTab
       BioTab(iv)%BioNum=ivv
    ENDDO

    READ(nfbioctb,'(a100)') junk
    READ(nfbioctb,'(a100)') junk
    DO iv=1,nv
       READ(nfbioctb,14) ivv,                    &
            BioTab(iv)%LTran(1,1),  &
            BioTab(iv)%LTran(2,1),  &
            BioTab(iv)%LTran(1,2),  &
            BioTab(iv)%LTran(2,2),  &
            BioTab(iv)%LRef(1,1),   &
            BioTab(iv)%LRef(2,1),   &
            BioTab(iv)%LRef(1,2),   &
            BioTab(iv)%LRef(2,2)
    ENDDO
    READ(nfbioctb,'(a100)') junk
    DO iv=1,nv
       READ(nfbioctb,12) ivv,                &
            BioTab(iv)%vmax0,         &
            BioTab(iv)%EffCon,  &
            BioTab(iv)%gsSlope, &
            BioTab(iv)%gsMin,        &
            BioTab(iv)%Atheta,  &
            BioTab(iv)%Btheta
    ENDDO
    READ(nfbioctb,'(a100)') junk
    DO iv=1,nv
       READ(nfbioctb,14) ivv,               &
            BioTab(iv)%TRDA,   &
            BioTab(iv)%TRDM,   &
            BioTab(iv)%TROP,   &
            BioTab(iv)%respcp, &               
            BioTab(iv)%SLTI,   &
            BioTab(iv)%HLTI,   &
            BioTab(iv)%SHTI,   &
            BioTab(iv)%HHTI
    ENDDO
    READ(nfbioctb,'(a100)') junk
    DO iv=1,nv
       READ(nfbioctb,14) ivv,                 &
            BioTab(iv)%SoRef(1), &
            BioTab(iv)%SoRef(2)
    ENDDO
    CLOSE(nfbioctb)
    !

12  FORMAT (i2,4x,e9.3,5(f8.3,1x))
13  FORMAT (i2,2x,6(f8.3,1x),f10.3)
14  FORMAT (i2,2x,8(f8.3,1x))
    !
    RETURN                                                                    
  END SUBROUTINE ReadBioTable



  !
  !=======================================================================        
  SUBROUTINE ReadAeroTables(fNameAeroTab)
    !=======================================================================
    ! This subroutine reads in interpolation tables of previously 
    ! calculated aerodynamic variables. 
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameAeroTab
    REAL(KIND=r4)   :: aero_zo_sib2         (nv,50,50)
    REAL(KIND=r4)   :: aero_zp_disp_sib2 (nv,50,50)
    REAL(KIND=r4)   :: aero_rbc_sib2         (nv,50,50)
    REAL(KIND=r4)   :: aero_rdc_sib2         (nv,50,50)
    REAL(KIND=r4)   :: aero_LAIgrid_sib2     (50)
    REAL(KIND=r4)   :: aero_fVCovergrid_sib2 (50)
    !
    ! open interpolation tables
    OPEN(nfbioctb,file=TRIM(fNameAeroTab), form='unformatted')
    !
    !   Read in interpolation tables
    READ (nfbioctb) aero_LAIgrid_sib2
    READ (nfbioctb) aero_fVCovergrid_sib2
    READ (nfbioctb) aero_zo_sib2      
    READ (nfbioctb) aero_zp_disp_sib2 
    READ (nfbioctb) aero_rbc_sib2     
    READ (nfbioctb) aero_rdc_sib2     
    LAIgrid            = REAL(aero_LAIgrid_sib2,KIND=r8)
    fVCovergrid     = REAL(aero_fVCovergrid_sib2,KIND=r8)
    AeroVar%zo            = REAL(aero_zo_sib2     ,KIND=r8)
    AeroVar%zp_disp = REAL(aero_zp_disp_sib2,KIND=r8)
    AeroVar%RbC     = REAL(aero_rbc_sib2    ,KIND=r8)
    AeroVar%RdC     = REAL(aero_rdc_sib2    ,KIND=r8)
    !,KIND=r8)
    ! close interpolation tables
    CLOSE (nfbioctb)
    !
    RETURN
  END SUBROUTINE ReadAeroTables

  !=======================================================================
  SUBROUTINE SoilProperties(TextClay,TextSand)
    !=======================================================================
    ! calculates soil physical properties given sand and clay content
    !
    ! Modifications
    !  Kevin Schaefer created subroutine for soil hydraulic properties (4/22/00)
    !  Kevin Schaefer resp. variable curve fits from Raich et al., 1991 (6/19/00)
    !  Kevin Schaefer combine code for hydraulic & respiration variables (3/30/01)
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: TextClay
    REAL(KIND=r8), INTENT(IN   ) :: TextSand
    !
    ! begin local variables
    REAL(KIND=r8) :: fclay   ! fraction of clay in soil
    REAL(KIND=r8) :: fsand   ! fraction of sand in soil

    Text%Clay=TextClay
    Text%Sand=TextSand
    !
    ! calculate Soil hydraulic and thermal variables based on Klapp and Hornberger
    SoilVar%PhiSat=-10.0_r8*10**(1.88_r8-0.0131_r8*Text%Sand)/1000.0_r8
    SoilVar%poros=0.489_r8-0.00126_r8*Text%Sand
    SoilVar%SatCo=0.0070556_r8*10**(-0.884_r8+0.0153_r8*Text%sand)/1000.0_r8
    SoilVar%bee=2.91_r8+0.159_r8*Text%Clay
    !
    ! Calculate clay and sand fractions from percentages
    fclay=Text%clay/100.0_r8
    fsand=Text%sand/100.0_r8
    !
    ! Calculate soil respiration variables based on curve fits to 
    ! data shown in Raich et al. (1991)
    SoilVar%Wopt=(-0.08_r8*fclay**2+0.22_r8*fclay+0.59_r8)*100.0_r8
    SoilVar%Skew=-2*fclay**3-0.4491_r8*fclay**2+0.2101_r8*fclay+0.3478_r8
    SoilVar%RespSat=0.25_r8*fclay+0.5_r8
    !
    ! assign value for mean slope of terrain
    SoilVar%Slope=0.176_r8
    !
    SoilVar%SoilNum=1
    RETURN
  END  SUBROUTINE SoilProperties
  !******************************************************************
  SUBROUTINE textclass(TextClay,TextSand)
    !******************************************************************
    ! Assigns soil texture classes based on the USDA texture triangle
    ! using subroutines developed by aris gerakis
    !******************************************************************
    !* +-----------------------------------------------------------------------
    !* |                         T R I A N G L E
    !* | Main program that calls WHAT_TEXTURE, a function that classifies soil
    !* | in the USDA textural triangle using sand and clay %
    !* +-----------------------------------------------------------------------
    !* | Created by: aris gerakis, apr. 98 with help from brian baer
    !* | Modified by: aris gerakis, july 99: now all borderline cases are valid
    !* | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
    !* |              main program
    !* +-----------------------------------------------------------------------
    !* | COMMENTS
    !* | Supply a data file with two columns, in free format: 1st column sand,
    !* |   2nd column clay %, no header.  The output is a file with the classes.
    !* +-----------------------------------------------------------------------
    !* | You may use, distribute and modify this code provided you maintain
    !* ! this header and give appropriate credit.
    !* +-----------------------------------------------------------------------
    !
    ! Modifications:
    !   Lara Prihodko customized triangle program for mapper (1/31/01)
    !
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: TextClay
    REAL(KIND=r8), INTENT(IN   ) :: TextSand
    !INTEGER    :: texture
    !integer    :: what_texture
    REAL(KIND=r8)       :: sand, clay
    REAL(KIND=r8)       :: silty_loam(1:7,1:2), sandy(1:7,1:2),&
         silty_clay_loam(1:7,1:2),&
         loam(1:7,1:2), clay_loam(1:7,1:2), sandy_loam(1:7,1:2),&
         silty_clay(1:7,1:2), sandy_clay_loam(1:7,1:2),&
         loamy_sand(1:7,1:2), clayey(1:7,1:2), silt(1:7,1:2),&
         sandy_clay(1:7,1:2)
    !LOGICAL :: inpoly

    !Initalize polygon coordinates:

    DATA sandy           /85.0_r8, 90.0_r8, 100.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,&
                          0.0_r8, 10.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA loamy_sand         /70.0_r8, 85.0_r8,  90.0_r8,85.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,  & 
                              15.0_r8, 10.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA sandy_loam         /50.0_r8, 43.0_r8,  52.0_r8,52.0_r8,80.0_r8,85.0_r8, 70.0_r8,  &
                               0.0_r8,  7.0_r8,   7.0_r8,20.0_r8,20.0_r8,15.0_r8,  0.0_r8/
    DATA loam                 /43.0_r8, 23.0_r8,  45.0_r8,52.0_r8,52.0_r8, 0.0_r8,  0.0_r8,    &
                               7.0_r8, 27.0_r8,  27.0_r8,20.0_r8, 7.0_r8, 0.0_r8,  0.0_r8/
    DATA silty_loam         / 0.0_r8,  0.0_r8,  23.0_r8,50.0_r8, 0.0_r8, 0.0_r8,  0.0_r8, 0.0_r8,    &
                              27.0_r8, 27.0_r8,   0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/ 
!   DATA silt                 /0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
    DATA sandy_clay_loam /52.0_r8, 45.0_r8, 45.0_r8, 65.0_r8, 80.0_r8, 0.0_r8, 0.0_r8,    & 
                              20.0_r8, 27.0_r8, 35.0_r8, 35.0_r8, 20.0_r8, 0.0_r8, 0.0_r8/
    DATA clay_loam         /20.0_r8, 20.0_r8, 45.0_r8, 45.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,          &
                              27.0_r8, 40.0_r8, 40.0_r8, 27.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA silty_clay_loam /0.0_r8, 0.0_r8, 20.0_r8, 20.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 27.0_r8,   &
                             40.0_r8, 40.0_r8, 27.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA sandy_clay         /45.0_r8, 45.0_r8, 65.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,          &
                              35.0_r8, 55.0_r8, 35.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA silty_clay         /0.0_r8, 0.0_r8, 20.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 40.0_r8,    &
                             60.0_r8, 40.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8/
    DATA clayey          /20.0_r8, 0.0_r8, 0.0_r8, 45.0_r8, 45.0_r8, 0.0_r8, 0.0_r8,          &
                              40.0_r8, 60.0_r8, 100.0_r8, 55.0_r8, 40.0_r8, 0.0_r8, 0.0_r8/

    !DATA silty_loam     /0 ,  0,  23, 50, 20, 8, 0, 12, 27, 27,  0, 0, 12, 0/
    !DATA silty_clay_loam/0 ,  0,  20, 20,  0, 0, 0, 27, 40, 40, 27, 0,  0, 0/
    !DATA loam           /43, 23,  45, 52, 52, 0, 0,  7, 27, 27, 20,  7, 0, 0/
    !DATA clay_loam      /20, 20,  45, 45,  0, 0, 0, 27, 40, 40, 27,  0, 0, 0/
    !DATA sandy_loam     /50, 43,  52, 52, 80,85,70,  0,  7,  7, 20, 20,15, 0/
    !DATA silty_clay     / 0,  0,  20,  0,  0, 0, 0, 40, 60, 40,  0,  0, 0, 0/
    !DATA sandy_clay_loam/52, 45,  45, 65, 80, 0, 0, 20, 27, 35, 35, 20, 0, 0/
    !DATA loamy_sand     /70, 85,  90, 85,  0, 0, 0,  0, 15, 10,  0,  0, 0, 0/
    !DATA clayey         /20,  0,   0, 45, 45, 0, 0, 40, 60,100, 55, 40, 0, 0/
    !DATA silt           / 0,  0,   8, 20,  0, 0, 0,  0, 12, 12,  0,  0, 0, 0/
    !DATA sandy_clay     /45, 45,  65,  0,  0, 0, 0, 35, 55, 35,  0,  0, 0, 0/

    text%sand=TextSand
    text%clay=TextClay
    sand = 0.0_r8
    clay = 0.0_r8

    !Read input:

    sand = text%sand
    clay = text%clay

    !Call function that estimates texture and put into structure:
    text%class = what_texture (sand, clay, silty_loam, sandy, &
         silty_clay_loam,                                                &
         loam, clay_loam, sandy_loam, silty_clay,                        &
         sandy_clay_loam, loamy_sand, clayey, silt,                &
         sandy_clay)

    RETURN
  END SUBROUTINE textclass
  !
  !******************************************************************
  !* +-----------------------------------------------------------------------
  !* | WHAT TEXTURE?
  !* | Function to classify a soil in the triangle based on sand and clay %
  !* +-----------------------------------------------------------------------
  !* | Created by: aris gerakis, apr. 98
  !* | Modified by: aris gerakis, june 99.  Now check all polygons instead of
  !* | stopping when a right solution is found.  This to cover all borderline
  !* | cases.
  !* +-----------------------------------------------------------------------

  FUNCTION what_texture (sand, clay, silty_loam, sandy,  &
       silty_clay_loam, loam, clay_loam,                      &
       sandy_loam, silty_clay, sandy_clay_loam,                     &
       loamy_sand, clayey, silt, sandy_clay)

    IMPLICIT NONE

    !Declare arguments:

    REAL(KIND=r8), INTENT(in) :: clay, sand, silty_loam(1:7,1:2),   &
         sandy(1:7,1:2),                                             &
         silty_clay_loam(1:7,1:2), loam(1:7,1:2),                     &
         clay_loam(1:7,1:2), sandy_loam(1:7,1:2),                     &
         silty_clay(1:7,1:2), sandy_clay_loam(1:7,1:2),             &
         loamy_sand(1:7,1:2), clayey(1:7,1:2), silt(1:7,1:2),   &
         sandy_clay(1:7,1:2)

    !Declare local variables:

    !logical :: inpoly
    INTEGER :: texture, what_texture

    !Find polygon(s) where the point is.

    texture = 0

    IF (sand .GT. 0.0 .AND. clay .GT. 0.0) THEN

       IF (inpoly(sandy, 3, sand, clay)) THEN
          texture = 1   ! sand
       ENDIF
       IF (inpoly(loamy_sand, 4, sand, clay)) THEN
          texture = 2   ! loamy sand
       ENDIF
       IF (inpoly(sandy_loam, 7, sand, clay)) THEN
          texture = 3 ! sandy loam
       ENDIF
       IF (inpoly(loam, 5, sand, clay)) THEN
          texture = 4   ! loam
       ENDIF
       IF (inpoly(silty_loam, 4, sand, clay)) THEN
          texture = 5   ! silt loam
       ENDIF
       IF (inpoly(sandy_clay_loam, 5, sand, clay)) THEN
          texture = 6   ! sandy clay loam
       ENDIF
       IF (inpoly(clay_loam, 4, sand, clay)) THEN
          texture = 7    ! clay loam
       ENDIF      
       IF (inpoly(silty_clay_loam, 4, sand, clay)) THEN
          texture = 8   ! silty clay loam
       ENDIF
       IF (inpoly(sandy_clay, 3, sand, clay)) THEN
          texture = 9  ! sandy clay
       ENDIF       
       IF (inpoly(silty_clay, 3, sand, clay)) THEN
          texture = 10 ! silty clay
       ENDIF
       IF (inpoly(clayey, 5, sand, clay)) THEN
          texture = 11 ! clay
       ENDIF
       IF (inpoly(silt, 4, sand, clay)) THEN
          texture = 5 ! silt loam
       ENDIF

    ENDIF
    IF (texture == 0) THEN
         texture = 5         ! silt loam
!
!        write (*, 1000) msand, mclay
! 1000   format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
    END IF

    !IF (sand == 100) THEN
    !   texture = 1
    !ENDIF

    !IF (clay == 100) THEN
    !   texture = 12
    !ENDIF

    !IF (texture == 6 ) THEN
    !   texture = 6
    !ENDIF

    what_texture = texture


  END FUNCTION what_texture
  !******************************************************************
  !--------------------------------------------------------------------------
  !                            INPOLY
  !   Function to tell if a point is inside a polygon or not.
  !--------------------------------------------------------------------------
  !   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
  !
  !   Please feel free to use this source code for any purpose, commercial
  !   or otherwise, as long as you don't restrict anyone else's use of
  !   this source code.  Please give credit where credit is due.
  !
  !   Point-in-polygon algorithm, created especially for World-Wide Web
  !   servers to process image maps with mouse-clickable regions.
  !
  !   Home for this file:  http://www.gcomm.com/develop/inpoly.c
  !
  !                                       6/19/95 - Bob Stein & Craig Yap
  !                                       stein@gcomm.com
  !                                       craig@cse.fau.edu
  !--------------------------------------------------------------------------
  !   Modified by:
  !   Aris Gerakis, apr. 1998: 1.  translated to Fortran
  !                            2.  made it work with real coordinates
  !                            3.  now resolves the case where point falls
  !                                on polygon border.
  !   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
  !   Aris Gerakis, july 1999: Now all borderline cases are valid
  !--------------------------------------------------------------------------
  !   Glossary:
  !   function inpoly: true=inside, false=outside (is target point inside
  !                    a 2D polygon?)
  !   poly(*,2):  polygon points, [0]=x, [1]=y
  !   npoints: number of points in polygon
  !   xt: x (horizontal) of target point
  !   yt: y (vertical) of target point
  !--------------------------------------------------------------------------

    LOGICAL FUNCTION inpoly (poly, npoints, xt, yt)

    IMPLICIT NONE

    !Declare arguments:

    INTEGER :: npoints
    REAL(KIND=r8), INTENT(in)    :: poly(7, 2), xt, yt

    !Declare local variables:

    REAL(KIND=r8)    :: xnew, ynew, xold, yold, x1, y1, x2, y2
    INTEGER :: i
    LOGICAL :: inside, on_border

    inside    = .FALSE.
    on_border = .FALSE.

    IF (npoints < 3)  THEN
       inpoly = .FALSE.
       RETURN
    END IF

    xold = poly(npoints,1)
    yold = poly(npoints,2)

    DO i = 1 , npoints
       xnew = poly(i,1)
       ynew = poly(i,2)

       IF (xnew > xold)  THEN
          x1 = xold
          x2 = xnew
          y1 = yold
          y2 = ynew
       ELSE
          x1 = xnew
          x2 = xold
          y1 = ynew
          y2 = yold
       END IF

       !The outer IF is the 'straddle' test and the 'vertical border' test.
       !The inner IF is the 'non-vertical border' test and the 'north' test.

       !The first statement checks whether a north pointing vector crosses
       !(stradles) the straight segment.  There are two possibilities, depe-
       !nding on whether xnew < xold or xnew > xold.  The '<' is because edge
       !must be "open" at left, which is necessary to keep correct count when
       !vector 'licks' a vertix of a polygon.

       IF ((xnew < xt .AND. xt <= xold) .OR. (.NOT. xnew < xt .AND. &
            .NOT. xt <= xold)) THEN
          !The test point lies on a non-vertical border:
          IF ((yt-y1)*(x2-x1) == (y2-y1)*(xt-x1)) THEN
             on_border = .TRUE.
             !Check if segment is north of test point.  If yes, reverse the
             !value of INSIDE.  The +0.001 was necessary to avoid errors due
             !arithmetic (e.g., when clay = 98.87 and sand = 1.13):
          ELSEIF ((yt-y1)*(x2-x1) < (y2-y1)*(xt-x1) + 0.001_r8) THEN
             inside = .NOT.inside ! cross a segment
          ENDIF
          !This is the rare case when test point falls on vertical border or
          !left edge of non-vertical border. The left x-coordinate must be
          !common.  The slope requirement must be met, but also point must be
          !between the lower and upper y-coordinate of border segment.  There
          !are two possibilities,  depending on whether ynew < yold or ynew >
          !yold:
       ELSEIF ((xnew == xt .OR. xold == xt) .AND. (yt-y1)*(x2-x1) == &
            (y2-y1)*(xt-x1) .AND. ((ynew <= yt .AND. yt <= yold) .OR.   &
            (.NOT. ynew < yt .AND. .NOT. yt < yold))) THEN
          on_border = .TRUE.
       ENDIF

       xold = xnew
       yold = ynew

    ENDDO

    !If test point is not on a border, the function result is the last state
    !of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
    !inside the polygon if it falls on any of its borders:

    IF (.NOT. on_border) THEN
       inpoly = inside
    ELSE
       inpoly = .TRUE.
    ENDIF
    !
  END FUNCTION inpoly
  !

  !***************************************************************************
  !                      (imonth,iday,iyear)
  INTEGER FUNCTION julday (imonth,iday,iyear)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: imonth
    INTEGER, INTENT(IN   ) :: iday
    INTEGER, INTENT(IN   ) :: iyear
    !
    ! compute the julian day from a normal date
    !
    julday= iday  &
         + MIN(1,MAX(0,imonth-1))*31  &
         + MIN(1,MAX(0,imonth-2))*(28+(1-MIN(1,MOD(iyear,4))))  &
         + MIN(1,MAX(0,imonth-3))*31  &
         + MIN(1,MAX(0,imonth-4))*30  &
         + MIN(1,MAX(0,imonth-5))*31  &
         + MIN(1,MAX(0,imonth-6))*30  &
         + MIN(1,MAX(0,imonth-7))*31  &
         + MIN(1,MAX(0,imonth-8))*31  &
         + MIN(1,MAX(0,imonth-9))*30  &
         + MIN(1,MAX(0,imonth-10))*31  &
         + MIN(1,MAX(0,imonth-11))*30  &
         + MIN(1,MAX(0,imonth-12))*31

  END FUNCTION julday

  !
  !=======================================================================
  SUBROUTINE FractionVegCover(NDVI,fPARmax,fPARmin, &
       SRmax,SRmin,fVCover)
    !=======================================================================
    ! calculates the vegetation cover fraction for a single pixel.
    ! The maximum fPAR for pixel during entire year determines vegetation
    ! cover fraction.  Maximum yearly NDVI corresponds to maximum fPAR.
    ! Calculate fPAR from Simple Ratio using an empirical linear relationship.
    !
    IMPLICIT NONE
    !
    ! begin input variables
    REAL(KIND=r8), INTENT(IN   ) :: NDVI      ! maximum FASIR NDVI value for a grid cell
    REAL(KIND=r8), INTENT(IN   ) :: fPARmax   ! Maximum possible FPAR corresponding to 98th percentile
    REAL(KIND=r8), INTENT(IN   ) :: fPARmin   ! Minimum possible FPAR corresponding to 2nd percentile
    REAL(KIND=r8), INTENT(IN   ) :: SRmax     ! Maximum simple ratio for biome type
    REAL(KIND=r8), INTENT(IN   ) :: SRmin     ! Minimum simple ratio for biome type
    !
    ! begin output variables
    REAL(KIND=r8), INTENT(OUT  ) ::  fVCover   ! fractional vegetation cover
    !
    ! begin internal variables
    REAL(KIND=r8) :: fPAR      ! maximum fPAR associated with maximum ndvi
    !
    ! The maximum fPAR for pixel during entire year determines vegetation
    ! cover fraction.  Maximum yearly NDVI corresponds to maximum fPAR.
    ! Calculate fPAR from Simple Ratio using an empirical linear relationship.
    !
    CALL srapar (NDVI,    &
         SRmin,    &
         SRmax,    &
         fPAR,     &
         fPARmax,  &
         fPARmin)
    !
    ! calculate fractional vegetation cover
    fVCover=fPAR/fPARmax
    !
    RETURN
  END SUBROUTINE FractionVegCover
  !
  !=======================================================================
  SUBROUTINE srapar (ndvi, SRmin, SRmax, fPAR, fPARmax, fParmin)
    !=======================================================================      
    ! calculates Canopy absorbed fraction of Photosynthetically 
    ! Active Radiation (fPAR) using the Simple Ratio (sr) method 
    ! (Los et al. (1998), eqn 6). This empirical method assumes a linear
    ! relationship between fPAR and sr.
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! begin input variables
    !
    REAL(KIND=r8), INTENT(IN   ) :: ndvi    ! normalized difference vegetation index
    REAL(KIND=r8), INTENT(IN   ) :: SRmin   ! minimum simple ratio for vegetation type
    REAL(KIND=r8), INTENT(IN   ) :: SRmax   ! maximum simple ratio for vegetation type
    REAL(KIND=r8), INTENT(IN   ) :: fPARmax ! Maximum possible FPAR corresponding to 98th percentile
    REAL(KIND=r8), INTENT(IN   ) :: fPARmin ! Minimum possible FPAR corresponding to 2nd percentile
    !
    ! begin output variables
    !
    REAL(KIND=r8), INTENT(OUT  ) :: fPAR    ! Canopy absorbed fraction of PAR
    !
    ! begin internal variables
    !
    REAL(KIND=r8)                :: sr      ! simple ratio of near IR and visible radiances
    !
    ! Calculate simple ratio (SR)
    !
    sr=(1.0_r8+ndvi)/(1.0_r8-ndvi)
    !
    ! Insure calculated SR value falls within physical limits for veg. type
    !
    sr=MAX(sr,SRmin)
    sr=MIN(sr,SRmax)
    !
    ! Calculate fPAR using SR method (Los et al. (1998), eqn 6)
    !
    fPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
    !
    RETURN                                                                    
  END SUBROUTINE srapar

  SUBROUTINE mapper(                    &   
       mapper_lat               , &   
       DOY               , &   
       mapper_prevNDVI          , &   
       mapper_curNDVI           , &   
       mapper_fVCover           , &   
       mapper_ChiL              , &   
       mapper_LTran             , &   
       mapper_LRef              , & 
       !mphtab_zc            , & 
       !mphtab_lwidth   , & 
       !mphtab_llength  , & 
       mphtab_laimax   , & 
       mphtab_stems    , & 
       mphtab_ndvimax  , & 
       mphtab_ndvimin  , & 
       mphtab_srmax    , & 
       mphtab_srmin    , & 
       aerotab_zo            , &
       aerotab_zp_disp   , &
       aerotab_rbc            , &
       aerotab_rdc            , &
       LAIgrid           , & 
       fVCovergrid       , & 
       timetab_fpar      , &    
       timetab_lai            , &   
       timetab_green     , &  
       timetab_zo            , & 
       timetab_zp_disp   , & 
       timetab_rbc            , & 
       timetab_rdc            , & 
       timetab_gmudmu    , & 
       ijmax               ) 
    !=======================================================================
    ! calculates time dependant boundary condition variables for SiB.
    !
    IMPLICIT NONE
    !
    ! begin input variables
    !
    INTEGER, INTENT(IN   ) :: ijmax
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_lat     (ijmax)     ! center latitude of grid cell
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_curNDVI (ijmax)     ! FASIR NDVI values for a grid cell
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_prevNDVI(ijmax)     ! previous month's NDVI value
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_fVCover (ijmax)     !
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_ChiL    (ijmax)     !
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_LTran   (ijmax,2,2) !
    REAL(KIND=r8)   , INTENT(IN   ) :: mapper_LRef    (ijmax,2,2) !
    INTEGER, INTENT(IN   ) :: DOY               ! Day of Year (DOY) of ndvi input map
    !
    ! begin input biome dependant, physical morphology variables
    !
    !REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_zc      (ijmax)
    !REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_LWidth  (ijmax)
    !REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_llength (ijmax)
    REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_LAImax  (ijmax)
    REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_stems   (ijmax)
    REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_NDVImax (ijmax)
    REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_NDVImin (ijmax)
    REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_SRmax   (ijmax)
    REAL(KIND=r8)   , INTENT(IN   ) :: mphtab_SRmin   (ijmax)
    !
    ! begin input aerodynamic parameters
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: aerotab_zo       (ijmax,50,50) 
    REAL(KIND=r8)   , INTENT(IN   ) :: aerotab_zp_disp  (ijmax,50,50) 
    REAL(KIND=r8)   , INTENT(IN   ) :: aerotab_rbc      (ijmax,50,50) 
    REAL(KIND=r8)   , INTENT(IN   ) :: aerotab_rdc      (ijmax,50,50) 
    !
    !  interpolation tables
    !
    REAL(KIND=r8)    , INTENT(IN  ) :: LAIgrid(50)    ! grid of LAI values for lookup table
    REAL(KIND=r8)    , INTENT(IN  ) :: fVCovergrid(50)! grid of fVCover values for interpolation table
    !
    ! begin time dependant, output variables
    !
    REAL(KIND=r8)    , INTENT(INOUT) :: timetab_fpar    (ijmax)  ! Canopy absorbed fraction of PAR
    REAL(KIND=r8)    , INTENT(INOUT) :: timetab_lai     (ijmax)  ! Leaf-area index
    REAL(KIND=r8)    , INTENT(INOUT) :: timetab_green   (ijmax)  ! Canopy greeness fraction of LAI
    REAL(KIND=r8)    , INTENT(OUT  ) :: timetab_zo      (ijmax)  ! Canopy roughness coeff
    REAL(KIND=r8)    , INTENT(OUT  ) :: timetab_zp_disp (ijmax)  ! Zero plane displacement
    REAL(KIND=r8)    , INTENT(OUT  ) :: timetab_rbc     (ijmax)  ! RB Coefficient (c1)
    REAL(KIND=r8)    , INTENT(OUT  ) :: timetab_rdc     (ijmax)  ! RC Coefficient (c2)
    REAL(KIND=r8)    , INTENT(INOUT) :: timetab_gmudmu  (ijmax)  ! Time-mean leaf projection    
    !
    ! begin internal variables
    !
    REAL(KIND=r8)            :: prevfPAR        (ijmax) ! previous month's fPAR value
    REAL(KIND=r8), PARAMETER :: fPARmax=0.95_r8
    !                   ! Maximum possible FPAR corresponding to 98th percentile
    REAL(KIND=r8), PARAMETER :: fPARmin=0.001_r8
    !                   ! Minimum possible FPAR corresponding to 2nd percentile
    !     For more information on fPARmin and fPARmax, see
    !     Sellers et al. (1994a, pg. 3532); Los (1998, pg. 29, 37-39)
    !

    !-----------------------------------------------------------------------
    ! Calculate time dependent variables
    !-----------------------------------------------------------------------

    !
    ! Calculate first guess fPAR 
    ! use average of Simple Ratio (SR) and NDVI methods.
    !
    !  print*,'call avgapar:',prevndvi,mphtab%ndvimin,mphtab%ndvimax
    !
    CALL AverageAPAR (mapper_prevNDVI         (1:ijmax) , & !IN
         mphtab_NDVImin (1:ijmax) , & !IN
         mphtab_NDVImax (1:ijmax) , & !IN
         mphtab_SRmin   (1:ijmax) , & !IN
         mphtab_SRmax   (1:ijmax) , & !IN
         fPARmax                    , & !IN
         fParmin                       , & !IN
         prevfPAR         (1:ijmax) , & !OUT
         ijmax                               ) !in

    CALL AverageAPAR (mapper_curNDVI                  (1:ijmax)   , & !IN
         mphtab_NDVImin         (1:ijmax)   , & !IN
         mphtab_NDVImax         (1:ijmax)   , & !IN
         mphtab_SRmin           (1:ijmax)   , & !IN
         mphtab_SRmax           (1:ijmax)   , & !IN
         fPARmax                              , & !IN
         fParmin                              , & !IN
         timetab_fPAR             (1:ijmax)   , & !OUT
         ijmax                                    ) !in
    !
    !
    ! Calculate leaf area index (LAI) and greeness fraction (Green)
    !   See S. Los et al 1998 section 4.2.
    !
    !   Select previous month
    !
    !
    CALL laigrn (timetab_fPAR              , &!IN
         prevfPAR                  , &!IN
         fPARmax                   , &!IN
         mapper_fVCover                   , &!IN
         mphtab_stems            , &!IN
         mphtab_LAImax           , &!IN
         timetab_Green             , &!OUT
         timetab_LAI               , &!OUT
         ijmax                       )!IN
    !
    ! Interpolate to calculate aerodynamic, time varying variables
    !
    !  PRINT*,'call aeroint:',timetab%lai,fvcover
    !
    CALL AeroInterpolate (                   &
         timetab_LAI                , &!IN 
         mapper_fVCover                    , &!IN
         LAIgrid                    , &!IN
         fVCovergrid                , &!IN
         aerotab_zo                 , &!IN
         aerotab_zp_disp            , &!IN
         aerotab_rbc                , &!IN
         aerotab_rdc                , &!IN
         timetab_zo                 , &!OUT
         timetab_zp_disp            , &!OUT
         timetab_RbC                , &!OUT
         timetab_RdC                , &!OUT
         ijmax                        )!IN
    !
    ! Calculate mean leaf orientation to par flux (gmudmu)
    !
    CALL gmuder (mapper_lat                       , &!IN
         DOY                       , &!IN
         mapper_ChiL                      , &!IN
         timetab_gmudmu            , &!OUT
         ijmax                         )!IN
    !
    ! recalculate fPAR adjusting for Sun angle, vegetation cover fraction,
    ! and greeness fraction, and LAI
    !
    CALL aparnew (timetab_LAI               , &!IN
         timetab_Green             , &!IN
         mapper_LTran                     , &!IN
         mapper_LRef                      , &!IN
         timetab_gmudmu            , &!IN
         mapper_fVCover                   , &!IN
         timetab_fPAR              , &!OUT
         fPARmax                   , &!IN
         fPARmin                   , &!IN
         ijmax                         )!IN

    !
    RETURN
  END SUBROUTINE mapper
  SUBROUTINE AverageAPAR (ndvi                          , &
       NDVImin                       , &
       NDVImax                       , &
       SRmin                         , &
       SRmax                         , &
       fPARmax                       , &
       fParmin                       , &
       fPAR                          , &
       ijmax                           )
    !=======================================================================      
    ! calculates Canopy absorbed fraction of Photosynthetically 
    ! Active Radiation (fPAR) using an average of the Simple Ratio (sr) 
    ! and NDVI methods (Los et al. (1999), eqn 5-6).  The empirical
    ! SR method assumes a linear relationship between fPAR and SR.
    ! The NDVI method assumes a linear relationship between fPAR and NDVI.
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! begin input variables
    !
    INTEGER         , INTENT(IN   ) :: ijmax
    REAL(KIND=r8)   , INTENT(IN   ) :: ndvi   (ijmax)  ! normalized difference vegetation index
    REAL(KIND=r8)   , INTENT(IN   ) :: NDVImin(ijmax)  ! minimum NDVI for vegetation type
    REAL(KIND=r8)   , INTENT(IN   ) :: NDVImax(ijmax)  ! maximum NDVI for vegetation type
    REAL(KIND=r8)   , INTENT(IN   ) :: SRmin  (ijmax)  ! minimum NDVI for vegetation type
    REAL(KIND=r8)   , INTENT(IN   ) :: SRmax  (ijmax)  ! maximum NDVI for vegetation type
    REAL(KIND=r8)   , INTENT(IN   ) :: fPARmax       ! Maximum possible FPAR corresponding to 98th percentile
    REAL(KIND=r8)   , INTENT(IN   ) :: fPARmin       ! Minimum possible FPAR corresponding to 2nd percentile
    !
    ! begin output variables
    !
    REAL(KIND=r8)   , INTENT(OUT  ) :: fPAR    (ijmax) ! Canopy absorbed fraction of PAR
    !
    ! begin internal variables
    !
    REAL(KIND=r8)                :: LocNDVI  ! local value of NDVI to prevent changes in input value
    REAL(KIND=r8)                :: sr          !  simple ratio of near IR and visible radiances
    REAL(KIND=r8)                :: NDVIfPAR ! fPAR from NDVI method
    REAL(KIND=r8)                :: SRfPAR   ! fPAR from SR method
    INTEGER             :: i
    DO i=1,ijmax
       !
       ! switch to local value of ndvi to prevent any changes going back to main
       !
       LocNDVI=NDVI(i)
       !
       ! Insure calculated NDVI value falls within physical limits for veg. type
       !
       !WRITE(*,*)NDVI(i),fPARmax,fPARmin,SRmax(i),SRmin(i)
       LocNDVI=MAX(LocNDVI,NDVImin(i))
       LocNDVI=MIN(LocNDVI,NDVImax(i))
       !  print*,'fpar1:',locndvi,ndvimin,ndvimax
       !
       ! Calculate simple ratio (SR)
       !
       sr=(1.0_r8+LocNDVI)/(1.0_r8-LocNDVI)
       !
       ! Calculate fPAR using SR method (Los et al. (1999), eqn 5)
       !
       SRfPAR=(sr-SRmin(i))*(fPARmax-fPARmin)/(SRmax(i)-SRmin(i))+fPARmin
       !  print*,'fpar2:',sr,srmin,fparmax,fparmin,srmin,srmax,fparmin
       !
       ! Calculate fPAR using NDVI method (Los et al. (1999), eqn 6)
       !
       NDVIfPAR=(LocNDVI-NDVImin(i))*(fPARmax-fPARmin)/(NDVImax(i)-NDVImin(i))+fPARmin
       !
       ! take average of two methods
       !
       fPAR(i)=0.50_r8*(SRfPAR+NDVIfPAR)
       !  print*,'fpar3:',fpar,srfpar,ndvifpar
       !
    END DO
    RETURN  
  END SUBROUTINE AverageAPAR


  !
  !======================================================================= 
  SUBROUTINE aparnew (LAI                              , &            
       Green                            , &            
       LTran                            , &            
       LRef                             , &            
       gmudmu                           , &            
       fVCover                          , &            
       fPAR                             , &            
       fPARmax                          , &            
       fPARmin                          , & 
       ijmax                              )                   
    !=======================================================================
    ! recomputes the Canopy absorbed fraction of Photosynthetically
    ! Active Radiation (fPAR), adjusting for solar zenith angle and the 
    ! vegetation cover fraction (fVCover) using a modified form of Beer's law..
    ! See Sellers et al. Part II (1996), eqns. 9-13.
    !
    IMPLICIT NONE
    !
    ! begin input variables
    !
    INTEGER, INTENT(IN   ) :: ijmax
    REAL(KIND=r8)   , INTENT(IN   ) :: LAI  (ijmax)    ! Leaf Area Index
    REAL(KIND=r8)   , INTENT(IN   ) :: Green(ijmax)    ! Greeness fraction of Leaf Area Index
    REAL(KIND=r8)   , INTENT(IN   ) :: LTran(ijmax,2,2)! Leaf transmittance for green/brown plants
    REAL(KIND=r8)   , INTENT(IN   ) :: LRef (ijmax,2,2)! Leaf reflectance for green/brown plants
    ! For LTran and LRef:
    ! (1,1)=shortwave, green plants
    ! (2,1)=longwave, green plants
    ! (1,2)=shortwave, brown plants
    ! (2,2)=longwave, brown plants
    REAL(KIND=r8)   , INTENT(IN   ) :: gmudmu (ijmax)  ! daily Time-mean canopy optical depth
    REAL(KIND=r8)   , INTENT(IN   ) :: fVCover(ijmax)  ! Canopy cover fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: fPARmax       ! Maximum possible FPAR corresponding to 98th percentile
    REAL(KIND=r8)   , INTENT(IN   ) :: fPARmin       ! Minimum possible FPAR corresponding to 2nd percentile
    !
    ! begin output variables
    !
    REAL(KIND=r8)   , INTENT(OUT  ) :: fPAR   (ijmax)  ! area average Canopy absorbed fraction of PAR
    !
    ! begin internal variables
    !
    REAL(KIND=r8)                   :: scatp         ! Canopy transmittance + reflectance coefficient wrt PAR
    REAL(KIND=r8)                   :: PARk          ! mean canopy absorption optical depth wrt PAR
    INTEGER                :: i
    DO i=1,ijmax
       !
       ! Calculate canopy transmittance + reflectance coefficient wrt PAR
       ! transmittance + reflectance coefficient=green plants + brown plants
       !
       scatp=Green(i)*(LTran(i,1,1)+LRef(i,1,1))+        &
            (1.0_r8-Green(i))*(LTran(i,1,2)+LRef(i,1,2))
       !
       ! Calculate PAR absorption optical depth in canopy adjusting for 
       ! variance in projected leaf area wrt solar zenith angle
       ! (Sellers et al. Part II (1996), eqn. 13b)
       ! PAR absorption coefficient=(1-scatp)
       !
       PARk=SQRT(1.0_r8-scatp)*gmudmu(i)
       !
       ! Calculate the new fPAR (Sellers et al. Part II (1996), eqn. 9)
       !
       fPAR(i)=fVCover(i)*(1.0_r8-EXP(-PARk*LAI(i)/fVCover(i)))
       !
       ! Ensure calculated fPAR falls within physical limits
       !
       fPAR(i)=MAX(fPARmin,fPAR(i))
       fPAR(i)=MIN(fPARmax,fPAR(i))
       !
    END DO
    RETURN
  END SUBROUTINE aparnew
  !
  !
  !=======================================================================      


  SUBROUTINE gmuder (Lat               , &
       DOY               , &
       ChiL              , &
       gmudmu            , &
       ijmax                 )
    !=======================================================================      
    ! calculates daily time mean optical depth of canopy relative to the Sun.
    !
    IMPLICIT NONE
    !
    ! begin input variables
    !
    INTEGER , INTENT(IN   ) :: ijmax
    REAL(KIND=r8)    , INTENT(IN   ) :: Lat   (ijmax)   ! latitude in degrees
    INTEGER , INTENT(IN   ) :: DOY             ! day-of-year (typically middle day of the month)
    REAL(KIND=r8)    , INTENT(IN   ) :: ChiL  (ijmax)   ! leaf angle distribution factor
    !
    ! begin output variables
    !
    REAL(KIND=r8)    , INTENT(OUT  ) :: gmudmu(ijmax)   ! daily time mean canopy optical depth relative to Sun
    !
    ! begin internal variables
    !
    REAL(KIND=r8)                    :: mumax    ! max cosine of the Solar zenith angle (noon)
    REAL(KIND=r8)                    :: mumin    ! min cosine of the Solar zenith angle (rise/set)
    REAL(KIND=r8)                    :: dec      ! declination of the Sun (Solar Declination)
    REAL(KIND=r8)                    :: pi180    ! conversion factor from degrees to radians
    REAL(KIND=r8)                    :: aa       ! minimum possible LAI projection vs. cosine Sun angle
    REAL(KIND=r8)                    :: bb       ! slope leaf area projection vs. cosine Sun angle
    INTEGER                 :: i
    !
    ! Calculate conversion factor from degrees to radians
    !
    !
    pi180=3.141590_r8/180.0_r8 
    !
    ! Calculate solar declination in degrees
    !
    !dec=23.50_r8*SIN(1.72e-20_r8*(DOY-80))
    dec=23.50_r8*SIN((360.0_r8/365.0_r8)*(284.0_r8+ DOY))
    !
    ! Calculate maximum cosine of zenith angle corresponding to noon
    !
    DO i=1,ijmax  
       mumax=COS((dec-lat(i))*pi180)
       mumax=MAX(0.020_r8, mumax)
       !
       ! Assign min cosine zenith angle corresponding to start disc set (cos(89.4))
       !
       mumin=0.010_r8
       !
       ! The projected leaf area relative to the Sun is G(mu)=aa+bb*mu
       ! Calculate minimum projected leaf area
       !
       aa=0.50_r8-0.6330_r8*ChiL(i)-0.330_r8*ChiL(i)*ChiL(i)
       !
       ! Calculate slope of projected leaf area wrt cosine sun angle
       !
       bb=0.8770_r8*(1.0_r8-2.0_r8*aa) 
       !
       ! Calculate mean optical depth of canopy by integrating G(mu)/mu over
       ! all values of mu.  Since G(mu) has an analytical form, this comes to
       !
       gmudmu(i)=aa*log(mumax/mumin)/(mumax-mumin)+bb
       !
    END DO
    RETURN                                                                    
  END SUBROUTINE gmuder

  SUBROUTINE laigrn (fPAR                    , & 
       fPARm                   , & 
       fPARmax                 , & 
       fVCover                 , & 
       stems                   , & 
       LAImax                  , & 
       Green                   , & 
       LAI                     , & 
       ijmax                       )
    !=======================================================================
    ! calculate leaf area index (LAI) and greenness fraction (Green) from fPAR. 
    ! LAI is linear with vegetation fraction and exponential with fPAR.
    ! See Sellers et al (1994), Equations 7 through 13.
    !                                                                               
    IMPLICIT NONE
    !
    ! begin input variables
    INTEGER , INTENT(IN   ) :: ijmax
    REAL(KIND=r8)    , INTENT(IN   ) :: fPAR    (ijmax)  ! fraction of PAR absorbed by plants at current time
    REAL(KIND=r8)    , INTENT(IN   ) :: fPARm   (ijmax)  ! fraction of PAR absorbed by plants at previous time
    REAL(KIND=r8)    , INTENT(IN   ) :: fPARmax          ! maximum possible FPAR corresponding to 98th percentile
    REAL(KIND=r8)    , INTENT(IN   ) :: fVCover (ijmax)  ! vegetation cover fraction
    REAL(KIND=r8)    , INTENT(IN   ) :: stems   (ijmax)  ! stem area index for the specific biome type
    REAL(KIND=r8)    , INTENT(IN   ) :: LAImax  (ijmax)  ! maximum total leaf area index for specific biome type
    !
    ! begin output variables
    !
    REAL(KIND=r8)    , INTENT(OUT  ) :: Green   (ijmax)  ! greeness fraction of the total leaf area index
    REAL(KIND=r8)    , INTENT(OUT  ) :: LAI     (ijmax)  ! area average total leaf area index
    !
    ! begin internal variables
    !
    REAL(KIND=r8)                    :: LAIg     ! green leaf area index at current time
    REAL(KIND=r8)                    :: LAIgm    ! green leaf area index at previous time
    REAL(KIND=r8)                    :: LAId     ! dead leaf area index at current time
    INTEGER                 :: i
    !
    ! Calculate current and previous green leaf area index (LAIg and LAIgm):
    ! LAIg is log-linear with fPAR.  Since measured fPAR is an area average, 
    ! divide by fVCover to get portion due to vegetation.  Since fVCover can
    ! be specified, check to assure that calculated fPAR does not exceed fPARMax.
    !
    !  print*,'lai_1',fpar,fvcover,fparmax
    DO i=1,ijmax
       IF(fPAR(i)/fVCover(i).GE.fPARmax) THEN
          LAIg=LAImax(i)
       ELSE
          LAIg=log(1.0_r8-fPAR(i)/fVCover(i))*LAImax(i)/log(1.0_r8-fPARmax)
       ENDIF
       !
       IF(fPARm(i)/fVCover(i).GE.fPARmax) THEN
          LAIgm=LAImax(i)
       ELSE
          LAIgm=log(1.0_r8-fPARm(i)/fVCover(i))*LAImax(i)/log(1.0_r8-fPARmax)
       ENDIF
       !
       ! Calculate dead leaf area index (LAId):
       ! If LAIg is increasing or unchanged, the vegetation is in growth mode.
       ! LAId is then very small (very little dead matter).
       ! If LAIg is decreasing, the peak in vegetation growth has passed and
       ! leaves have begun to die off.  LAId is then half the change in LAIg,
       ! assuming half the dead leaves fall off.
       !
       !     Growth mode dead leaf area index:
       IF (LAIg.GE.LAIgm) LAId=0.00010_r8
       !
       !     die-off (post peak growth) dead leaf area index:
       IF (LAIg.LT.LAIgm) LAId=0.50_r8*(LAIgm-LAIg)
       !
       ! Calculate area average, total leaf area index (LAI):
       LAI(i)=(LAIg+LAId+stems(i))*fVCover(i)
       !  print*,'laigrn1',laig,laid,stems,fvcover
       !
       ! Calculate greeness fraction (Green):
       ! Greeness fraction=(green leaf area index)/(total leaf area index)
       Green(i)=LAIg/(LAIg+LAId+stems(i))
       !  PRINT*,'end laigrn',LAI,Green,laimax
    END DO
    RETURN                                                                    
  END SUBROUTINE laigrn
  !
  !=======================================================================
  SUBROUTINE AeroInterpolate (LAI                     , & 
       fVCover                 , & 
       LAIgrid                 , & 
       fVCovergrid             , & 
       AeroVar_zo              , & 
       AeroVar_zp_disp         , & 
       AeroVar_rbc             , & 
       AeroVar_rdc             , & 
       zo                      , & 
       zp_disp                 , & 
       RbC                     , & 
       RdC                     , & 
       ijmax                     ) 
    !=======================================================================
    ! This subroutine calculates the aerodynamic parameters by bi-linear 
    ! interpolation from a lookup table of previously calculated values.  
    ! The interpolation table is a numpts x numpts LAI/fVCover grid with
    ! LAI ranging from 0.02 to 10 and fVCover ranging from 0.01 to 1.
    !
    IMPLICIT NONE
    !
    ! begin input variables
    INTEGER , INTENT(IN   ) :: ijmax
    REAL(KIND=r8)    , INTENT(IN   ) :: LAI        (ijmax)    ! actual area averaged LAI for interpolation
    REAL(KIND=r8)    , INTENT(IN   ) :: fVCover    (ijmax)    ! vegetation cover fraction for interpolation
    REAL(KIND=r8)    , INTENT(IN   ) :: LAIgrid    (50)       ! grid of LAI values for lookup table 
    REAL(KIND=r8)    , INTENT(IN   ) :: fVCovergrid(50)       ! grid of fVCover values for interpolation table

    REAL(KIND=r8)   , INTENT(IN   ) :: AeroVar_zo           (ijmax,50,50)
    REAL(KIND=r8)   , INTENT(IN   ) :: AeroVar_zp_disp     (ijmax,50,50)
    REAL(KIND=r8)   , INTENT(IN   ) :: AeroVar_rbc           (ijmax,50,50)
    REAL(KIND=r8)   , INTENT(IN   ) :: AeroVar_rdc           (ijmax,50,50)
    !
    ! begin output variables
    !
    REAL(KIND=r8)    , INTENT(OUT  ) ::  RbC     (ijmax)      ! interpolated Rb coefficient
    REAL(KIND=r8)    , INTENT(OUT  ) ::  RdC     (ijmax)      ! interpolated Rd coefficient
    REAL(KIND=r8)    , INTENT(OUT  ) ::  zo      (ijmax)      ! interpolated roughness length
    REAL(KIND=r8)    , INTENT(OUT  ) ::  zp_disp (ijmax)      ! interpolated zero plane displacement
    !
    ! begin internal variables
    !
    INTEGER                          :: i              ! index for LAI grid location
    INTEGER                          :: j              ! index for fVCover grid location
    REAL(KIND=r8)                    :: LocLAI         ! local LAI var. to prevent changing main LAI value
    REAL(KIND=r8)                    :: LocfVCover     ! local fVCover var. to prevent changing fVCover value
    REAL(KIND=r8)                    :: DLAI           ! grid spacing between LAI values in tables
    REAL(KIND=r8)                    :: DfVCover       ! grid spacing between fVCover values in tables
    INTEGER                 :: ii
    !
    !  !print*,'aerointerp:lai,fvc=',lai,fvcover
    !
    ! calculate grid spacing (assumed fixed)
    !
    DLAI=LAIgrid(2)-LAIgrid(1)
    DfVCover=fVCovergrid(2)-fVCovergrid(1)
    !
    ! Assign input LAI and fVCover to local variables and make sure
    ! they lie within the limits of the interpolation tables, assuring 
    ! the LAI and fVCover values returned from the subroutine are not modified.
    !
    DO ii=1,ijmax
       LocLAI    =MAX(LAI(ii),0.020_r8)
       LocfVCover=MAX(fVCover(ii),0.010_r8)
       !
       ! determine the nearest array location for the desired LAI and fVCover
       !
       i=INT(LocLAI/DLAI+1)
       j=INT(LocfVCover/DfVCover+1)
       j=MIN(j,49)
       !
       ! interpolate RbC variable
       !
       CALL interpolate(                                      &
            LAIgrid(i)                                      , & !IN
            LocLAI                                          , & !IN
            DLAI                                            , & !IN
            fVCovergrid(j)                                  , & !IN
            LocfVCover                                      , & !IN
            DfVCover                                        , & !IN
            AeroVar_RbC (ii,i,j    )                               , & !IN
            AeroVar_RbC (ii,i+1,j  )                               , & !IN
            AeroVar_RbC (ii,i,j+1  )                               , & !IN
            AeroVar_RbC (ii,i+1,j+1)                               , & !IN
            RbC         (ii)                                  ) !OUT
       !
       ! interpolate RdC variable
       !
       CALL interpolate(                                      &
            LAIgrid(i)                                      , & !IN
            LocLAI                                          , & !IN
            DLAI                                            , & !IN
            fVCovergrid(j)                                  , & !IN
            LocfVCover                                      , & !IN
            DfVCover                                        , & !IN
            AeroVar_RdC (ii,i,j    )                               , & !IN
            AeroVar_RdC (ii,i+1,j  )                               , & !IN
            AeroVar_RdC (ii,i,j+1  )                               , & !IN
            AeroVar_RdC (ii,i+1,j+1)                               , & !IN
            RdC         (ii)                                  ) !OUT
       !
       ! interpolate roughness length'
       !
       CALL interpolate(                                      &
            LAIgrid(i)                                      , & !IN
            LocLAI                                          , & !IN
            DLAI                                            , & !IN
            fVCovergrid(j)                                  , & !IN
            LocfVCover                                      , & !IN
            DfVCover                                        , & !IN
            AeroVar_zo(ii,i,j    )                            , & !IN
            AeroVar_zo(ii,i+1,j  )                            , & !IN
            AeroVar_zo(ii,i,j+1  )                            , & !IN
            AeroVar_zo(ii,i+1,j+1)                            , & !IN
            zo        (ii)                                    ) !OUT
       !
       ! interpolate zero plane displacement
       !
       CALL interpolate(                                      &
            LAIgrid(i)                                      , &!IN
            LocLAI                                          , &!IN
            DLAI                                            , &!IN
            fVCovergrid(j)                                  , &!IN
            LocfVCover                                      , &!IN
            DfVCover                                        , &!IN
            AeroVar_zp_disp (ii,i,j    )                       , &!IN
            AeroVar_zp_disp (ii,i+1,j  )                       , &!IN
            AeroVar_zp_disp (ii,i,j+1  )                       , &!IN
            AeroVar_zp_disp (ii,i+1,j+1)                       , &!IN
            zp_disp         (ii)                              )!OUT
       !
    END DO
    RETURN
  END SUBROUTINE AeroInterpolate

  !
  !=======================================================================
  !=======================================================================
  SUBROUTINE interpolate(x1                        , & 
       x                         , & 
       Dx                        , & 
       y1                        , & 
       y                         , & 
       Dy                        , & 
       z11                       , & 
       z21                       , & 
       z12                       , & 
       z22                       , & 
       z                           ) 
    !=======================================================================

    IMPLICIT NONE

    ! calculates the value of z=f(x,y) by linearly interpolating
    ! between the 4 closest data points on a uniform grid.  The subroutine
    ! requires a grid point (x1, y1), the grid spacing (Dx and Dy), and the 
    ! 4 closest data points (z11, z21, z12, and z22).
    !
    ! begin input variables
    !
    REAL(KIND=r8)    , INTENT(IN   ) :: x1  ! the x grid location of z11
    REAL(KIND=r8)    , INTENT(IN   ) :: x   ! x-value at which you will interpolate z=f(x,y)
    REAL(KIND=r8)    , INTENT(IN   ) :: Dx  ! grid spacing in the x direction
    REAL(KIND=r8)    , INTENT(IN   ) :: y1  ! the y grid location of z11
    REAL(KIND=r8)    , INTENT(IN   ) :: y   ! y-value at which you will interpolate z=f(x,y)
    REAL(KIND=r8)    , INTENT(IN   ) :: Dy  ! grid spacing in the y direction
    REAL(KIND=r8)    , INTENT(IN   ) :: z11 ! f(x1, y1)
    REAL(KIND=r8)    , INTENT(IN   ) :: z21 ! f(x1+Dx, y1)
    REAL(KIND=r8)    , INTENT(IN   ) :: z12 ! f(x1, y1+Dy)
    REAL(KIND=r8)    , INTENT(IN   ) :: z22 ! f(x1+Dx, y1+Dy)
    !
    ! begin output variables
    !
    REAL(KIND=r8)    , INTENT(OUT  ) :: z   ! f(x,y), the desired interpolated value
    !
    ! begin internal variables
    !
    REAL(KIND=r8)                    :: zp  ! z'=first interpolated value at (x, y1)
    REAL(KIND=r8)                    :: zpp ! z''=second interpolated value at (x, Y1+Dy)
    !
    ! interpolate between z11 and z21 to calculate z' (zp) at (x, y1)
    !
    zp=z11+(x-x1)*(z21-z11)/Dx
    !
    ! interpolate between z12 and z22 to calculate z'' (zpp) at (x, Y1+Dy)
    !
    zpp=z12+(x-x1)*(z22-z12)/Dx
    !
    ! interpolate between zp and zpp to calculate z at (x,y)
    !
    z=zp+(y-y1)*(zpp-zp)/Dy
    !
    RETURN
  END SUBROUTINE interpolate

  SUBROUTINE SiB2_Driver(&
         nCols              ,nmax                 ,kMax                ,&
         latco              ,ktm                  ,initlz              ,&
         kt                 ,iswrad               ,ilwrad               ,dtc3x               ,&
         intg               ,tkes_sib             ,t_sib                ,sh_sib              ,&
         gu                 ,gv                   ,pl2g_sib             ,imask               ,&
         prsi               ,prsl                 ,phii                 ,phil,&
         zenith             ,beam_visb            ,beam_visd            ,beam_nirb           ,&
         beam_nird          ,cos2                 ,dlwbot_sib           ,radvbc_sib          ,&
         radvdc_sib         ,radnbc_sib           ,radndc_sib           ,dlspr_sib           ,&
         dcupr_sib          ,itype                ,slrad               ,&
         qsurf              ,colrad                           ,&
         MskAnt             ,tsea                 ,tseam               ,&
         tsurf              ,tmtx                 ,qmtx                 ,umtx                ,&
         cu                 ,ustar                ,cosz_sib             ,hr                  ,&
         ect                ,eci                  ,egt                  ,egi                 ,&
         egs                ,ec                   ,eg                   ,hc                  ,&
         hg                 ,chf                  ,shf                  ,roff                ,&
         drag               ,ra                   ,rb                   ,rd                  ,&
         rc                 ,rg                   ,ta                   ,ea                  ,&
         etc                ,etg                  ,&
         rsoil              ,tg                   ,ndvi                 ,ndvim               ,&
         sens               ,evap                 ,&
         umom               ,vmom                 ,zorl                 ,rmi                 ,&
         rhi                ,cond                 ,stor                 ,z0d                 ,&
         spdm_sib           ,Ustarm               ,z0sea                ,rho                 ,&
         dd                 ,qsfc                 ,tsfc                 ,bstar               ,&
         HML                ,HUML                 ,HVML                 ,TSK                 ,&
         cldtot             ,ySwSfcNet            ,LwSfcNet             ,pblh                ,&
         QCF                ,QCL                  ,sm0                  ,mlsi                ,&
         LwSfcDown          ,month                ,Mmlen                ,idatec,dump)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: nmax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: latco
    INTEGER, INTENT(IN   ) :: ktm
    INTEGER, INTENT(IN   ) :: initlz
    INTEGER, INTENT(IN   ) :: kt
    INTEGER, INTENT(IN   ) :: intg
    INTEGER      , INTENT(IN   ) :: idatec(1:4) 
    CHARACTER(len=*), INTENT(IN   ) :: iswrad
    CHARACTER(len=*), INTENT(IN   ) :: ilwrad
    REAL(KIND=r8)   , INTENT(IN   ) :: dtc3x
    REAL(KIND=r8)   , INTENT(IN   ) :: tkes_sib   (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: t_sib      (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: sh_sib     (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gu         (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gv         (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: pl2g_sib   (1:nCols)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask      (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: prsi       (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8)   , INTENT(IN   ) :: prsl       (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8)   , INTENT(IN   ) :: phii       (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: phil       (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: zenith     (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: beam_visb(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: beam_visd(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: beam_nirb(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: beam_nird(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cos2      (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: dlwbot_sib (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: radvbc_sib (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: radvdc_sib (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: radnbc_sib (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: radndc_sib (1:nCols)
    INTEGER         , INTENT(IN   ) :: itype      (1:nCols)
    !REAL(KIND=r8)   , INTENT(IN   ) ::  ssib       (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: dcupr_sib  (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: dlspr_sib  (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad     (1:nCols)
    INTEGER(KIND=i8), INTENT(IN   ) :: MskAnt     (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: spdm_sib   (1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cosz_sib   (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tseam      (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: zorl(1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: qsfc     (ncols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsfc     (ncols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsea     (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: slrad    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsurf    (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: qsurf    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: umom     (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: vmom     (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rmi      (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rhi      (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: cond     (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: stor     (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: Ustarm   (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: z0sea    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rho      (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tmtx           (1:nCols,1:kMax,1:3)
    REAL(KIND=r8)   , INTENT(INOUT) :: qmtx           (1:nCols,1:kMax,1:3)
    REAL(KIND=r8)   , INTENT(INOUT) :: umtx           (1:nCols,1:kMax,1:4)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: cu        (1:nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ustar     (1:nCols)   
    REAL(KIND=r8)    ,INTENT(IN OUT) ::dump(1:nCols,1:kMax )
    REAL(KIND=r8)    ,INTENT(IN OUT) :: hr    (nCols)!,INTENT(INOUT        )
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ect   (nCols)     ! transpiration (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: eci   (nCols)     ! canopy interception evaporation (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: egt   (nCols)   
    REAL(KIND=r8)    ,INTENT(IN OUT) :: egi   (nCols)     ! ground interception evaporation (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: egs   (nCols)     ! soil evaporation (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: eg    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ec    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: hc    (nCols)     
    REAL(KIND=r8)    ,INTENT(IN OUT) :: hg    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: chf   (nCols)     ! canopy heat flux (W/m^2)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: shf   (nCols)     ! soil heat flux (W/m^2)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: roff  (nCols)   ! total runoff (surface and subsurface)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ra    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rb    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rd    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rc    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rg    (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ea    (nCols)      ! canopy airspace water vapor pressure (hPa)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: etc   (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: etg   (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rsoil (nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: drag  (1:nCols,2)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ta    (1:nCols)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: tg    (1:nCols)
    REAL(KIND=r8)  :: tc    (1:nCols)  
    REAL(KIND=r8)  :: td1   (1:nCols)  
    REAL(KIND=r8)  :: capac (1:nCols,2)
    REAL(KIND=r8)  :: w     (1:nCols,3)
    REAL(KIND=r8)    ,INTENT(INOUT) :: ndvi        (1:nCols)             
    REAL(KIND=r8)    ,INTENT(INOUT) :: ndvim    (1:nCols)
    REAL(KIND=r8)    ,INTENT(INOUT) :: sens     (nCols) !sensible heat flux
    REAL(KIND=r8)    ,INTENT(INOUT) :: evap     (nCols) !latent heat flux
    REAL(KIND=r8)    ,INTENT(INOUT) :: z0d        (1:nCols)
    REAL(KIND=r8)    ,INTENT(INOUT) :: dd        (1:nCols)                
    REAL(KIND=r8)    ,INTENT(OUT) :: bstar      (1:nCols)                
    REAL(KIND=r8),    INTENT(INOUT) :: HML  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HUML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HVML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: TSK  (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cldtot (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ySwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: LwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCF(nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCL(nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: sm0      (ncols,3)
    INTEGER(KIND=i8), INTENT(INOUT) :: mlsi  (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: LwSfcDown(1:nCols )
    INTEGER         , INTENT(IN   ) :: month(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: Mmlen(1:nCols )


    REAL(KIND=r8)                   :: xsea       (1:nCols) 
    REAL(KIND=r8)                   :: xndvi      (1:nCols)
    REAL(KIND=r8)                   :: ps         (1:nCols)
    REAL(KIND=r8)                   :: psb        (1:nCols)
    REAL(KIND=r8)                   :: tkes       (1:nCols)
    REAL(KIND=r8)                   :: bpeg       (1:nCols,2)
    REAL(KIND=r8)                   :: bps        (1:nCols)
    REAL(KIND=r8)                   :: ts         (1:nCols)
    REAL(KIND=r8)                   :: rhoair     (1:nCols)
    REAL(KIND=r8)                   :: cupr         (1:nCols)
    REAL(KIND=r8)                   :: lspr         (1:nCols)
    REAL(KIND=r8)                   :: radvbc     (1:nCols)
    REAL(KIND=r8)                   :: radnbc     (1:nCols)
    REAL(KIND=r8)                   :: radvdc     (1:nCols)
    REAL(KIND=r8)                   :: radndc     (1:nCols)
    REAL(KIND=r8)                   :: zb         (1:nCols)
    REAL(KIND=r8)                   :: spdm       (1:nCols)
    REAL(KIND=r8)                   :: dlwbot     (1:nCols)
    REAL(KIND=r8)                   :: tm         (1:nCols)
    REAL(KIND=r8)                   :: qm         (1:nCols)
    REAL(KIND=r8)                   :: sh         (1:nCols)
    REAL(KIND=r8)                   :: radt  (nCols,3)
    REAL(KIND=r8)                   :: qa         (1:nCols)
    REAL(KIND=r8)                   :: rst         (1:nCols)
    REAL(KIND=r8)                   :: snow         (1:nCols,2)
    REAL(KIND=r8)                   :: zlt        (1:nCols)
    REAL(KIND=r8)                   :: z1        (1:nCols)
    REAL(KIND=r8)                   :: z2        (1:nCols)
    REAL(KIND=r8)                   :: cc1        (1:nCols)
    REAL(KIND=r8)                   :: cc2        (1:nCols)
    REAL(KIND=r8)                   :: poros        (1:nCols)         
    REAL(KIND=r8)                   :: rootd        (1:nCols)
    REAL(KIND=r8)                   :: SoDep        (1:nCols)    
    REAL(KIND=r8)                   :: zdepth    (1:nCols,3)
    REAL(KIND=r8)                   :: phsat        (1:nCols)
    REAL(KIND=r8)                   :: bee        (1:nCols)
    REAL(KIND=r8)                   :: respcp    (1:nCols)
    REAL(KIND=r8)                   :: vmax0        (1:nCols)
    REAL(KIND=r8)                   :: green        (1:nCols)
    REAL(KIND=r8)                   :: tran      (1:nCols,2,2)
    REAL(KIND=r8)                   :: ref       (1:nCols,2,2)
    REAL(KIND=r8)                   :: gmudmu    (1:nCols)
    REAL(KIND=r8)                   :: trop        (1:nCols)
    REAL(KIND=r8)                   :: phc        (1:nCols)
    REAL(KIND=r8)                   :: trda        (1:nCols)
    REAL(KIND=r8)                   :: trdm        (1:nCols)
    REAL(KIND=r8)                   :: slti        (1:nCols)
    REAL(KIND=r8)                   :: shti        (1:nCols)
    REAL(KIND=r8)                   :: hltii        (1:nCols)
    REAL(KIND=r8)                   :: hhti        (1:nCols)
    REAL(KIND=r8)                   :: effcon    (1:nCols)
    REAL(KIND=r8)                   :: binter    (1:nCols)
    REAL(KIND=r8)                   :: gradm        (1:nCols)
    REAL(KIND=r8)                   :: atheta    (1:nCols)
    REAL(KIND=r8)                   :: btheta    (1:nCols)
    REAL(KIND=r8)                   :: aparc        (1:nCols)
    REAL(KIND=r8)                   :: wopt        (1:nCols)
    REAL(KIND=r8)                   :: zm        (1:nCols)
    REAL(KIND=r8)                   :: wsat        (1:nCols)
    REAL(KIND=r8)                   :: vcover    (1:nCols)
    REAL(KIND=r8)                   :: radfac    (1:nCols,2,2,2)
    REAL(KIND=r8)                   :: thermk    (1:nCols)
    REAL(KIND=r8)                   :: satco     (1:nCols)
    REAL(KIND=r8)                   :: slope     (1:nCols)
    REAL(KIND=r8)                   :: chil      (1:nCols)
    REAL(KIND=r8)                   :: sandfrac2 (1:nCols)
    REAL(KIND=r8)                   :: clayfrac2 (1:nCols)
    REAL(KIND=r8)                   :: ztdep     (1:nCols,nsoil)
    REAL(KIND=r8)                   :: ventmf    (1:nCols)
    REAL(KIND=r8)                   :: thvgm     (1:nCols)
    REAL(KIND=r8)                   :: xgpp      (1:nCols)
    REAL(KIND=r8)                   :: pco2ap2   (1:nCols)
    REAL(KIND=r8)                   :: dttsib
    REAL(KIND=r8)                   :: c4fract   (1:nCols)
    REAL(KIND=r8)                   :: d13cresp  (1:nCols)
    REAL(KIND=r8)                   :: d13cca    (1:nCols)
    REAL(KIND=r8)                   :: td0       (1:nCols,nsoil)
    REAL(KIND=r8)                   :: td        (1:nCols,nsoil)
    REAL(KIND=r8)                   :: tdm       (1:nCols,nsoil)
    REAL(KIND=r8)                   :: pco2m     (1:nCols)
    INTEGER                         :: ncount,i,k
    REAL(KIND=r8)                   :: fss       (1:nCols)
    REAL(KIND=r8)                   :: fws       (1:nCols)
    REAL(KIND=r8)                   :: es       (1:nCols)
    REAL(KIND=r8)                   :: hs       (1:nCols)
    REAL(KIND=r8)                   :: co2flx    (1:nCols)
    REAL(KIND=r8)                   :: ct        (1:nCols)
    REAL(KIND=r8)                   :: cosz      (1:nCols)
    REAL(KIND=r8)                   :: SoRef     (1:nCols,2)
    REAL(KIND=r8)                   :: salb      (1:nCols,2,2)
    REAL(KIND=r8)                   :: tgeff4_sib(1:nCols)
    REAL(KIND=r8)                   :: r100,cpdgrv,tau,snofac
    REAL(KIND=r8)                   :: smal2 
    REAL(KIND=r8)                   :: tmin    (1:ncols)
    REAL(KIND=r8)                   :: tmax    (1:ncols)
    REAL(KIND=r8)                   :: thm     (1:ncols)
    REAL(KIND=r8)                   :: areas   (1:ncols)     
    REAL(KIND=r8)                   :: bstar1   (1:ncols)    
    REAL(KIND=r8)                   :: zlwup  (1:ncols)    
    REAL(KIND=r8) :: um    (nCols)
    REAL(KIND=r8) :: vm    (nCols)
    REAL(KIND=r8) :: gmt   (ncols,3)
    REAL(KIND=r8) :: gmq   (ncols,3)
    REAL(KIND=r8) :: gmu   (ncols,4)
    REAL(KIND=r8) :: sinclt (nCols)
    REAL(KIND=r8)   :: GSW (nCols)
    REAL(KIND=r8)   :: GLW (nCols)
    REAL(KIND=r8)    :: sigki      (nCols)           
    REAL(KIND=r8)    :: delsig     (nCols)
    REAL(KIND=r8)    :: rbyg
    LOGICAL                         :: InitMod
    INTEGER                         :: ioffset   ! subdomain offset
    INTEGER                         :: nsib
    INTEGER                         :: nint
    INTEGER                         :: IntSib
    INTEGER                         :: itr,ind
    smal2 = 1.0e-3_r8
    ioffset = 1
    r100=100.0e0_r8 /rgas
    dttsib  = dtc3x
    tau      = dttsib  
    cpdgrv = cp/grav  
          !rbyg=gasr/grav*delsig(1)*0.5e0_r8
          !ztn   (ncount)=MAX((rbyg * tvland(ncount) ),0.5_r8)
          !ztn   (ncount)=MAX((phii(i,2) - phii(i,1)),0.5_r8)

   ! rbyg  =rgas/grav*delsig(1)*0.5_r8
    radt=0.0_r8
    ncount=0
    InitMod = (initlz >= 0 .AND. ktm == -1 .AND. kt == 0 .AND. nmax >= 1)
    DO i=1,nCols
       GSW(I) = radvbc_sib  (i)+radvdc_sib  (i)+radnbc_sib  (i)+radndc_sib  (i)
       GLW(I) =  dlwbot_sib(i) 
    END DO
    DO i=1,nCols
       spdm_sib(i) = SQRT((gu(i,1) /SIN(colrad(i)))**2  + (gv(i,1) /SIN(colrad(i)))**2)
       spdm_sib(i) = MAX(.1_r8 ,spdm_sib(i))
       IF(iMask(i) >= 1_i8)THEN 
          ncount=ncount+1
          um      (ncount)  = gu(i,1) /SIN(colrad(i))
          vm      (ncount)  = gv(i,1) /SIN(colrad(i))
          sinclt  (ncount)  = SIN(colrad(i))
          !bps     (ncount)  = sigki(1)
          delsig(ncount)    =((prsi(i,1)/prsi(i,1)) - (prsi(i,2)/prsi(i,1)) )
          rbyg  =rgas/grav*delsig(ncount)*0.5_r8

          sigki   (ncount)= (prsi(i,1)/(prsi(i,2)))**(gasr/cp)
          bps     (ncount)= (prsi(i,1)/(prsi(i,2)))**(gasr/cp)
          thm     (ncount)  = t_sib    (i,1) * bps     (ncount) 
          tm      (ncount)  = t_sib    (i,1)    
          ps      (ncount)  = pl2g_sib (i) ! driver pressure (mb)
          qm      (ncount)  = sh_sib   (i,1)
          sh      (ncount)  = sh_sib   (i,1)
          !psb     (ncount)  = pl2g_sib(i)*delsig
          psb     (ncount)  = pl2g_sib(i)*(grav/(rgas*t_sib(i,1)))*rbyg*(tm(ncount)*(1.0_r8+0.608_r8*qm(ncount)))
          !psb     (ncount)  = psb_g    (ncount)           ! pbl depth (mb)
          tkes    (ncount)  = 0.00000001*( 1 + 0.0_r8*tkes_sib (i))
          bpeg    (ncount,1)= bps     (ncount) ! (0.001_r8* ps(ncount))**kapa! (ps/1000)**kapa--> bps(ncount)=sigki(1)
          bpeg    (ncount,2)= ((ps(ncount)-psb(ncount))/ps(ncount))**(-kapa)
          ts      (ncount)  = t_sib    (i,1) * bps     (ncount) !tm(ncount) / bpeg(ncount,1)
          rhoair  (ncount)  = r100*ps(ncount)/t_sib    (i,1)
          cupr    (ncount)  = dcupr_sib (i)
          lspr    (ncount)  = dlspr_sib (i)
          radvbc  (ncount)  = radvbc_sib(i)
          radnbc  (ncount)  = radnbc_sib(i)
          radvdc  (ncount)  = radvdc_sib(i)
          radndc  (ncount)  = radndc_sib(i)
          zb      (ncount)  = rbyg*(tm(ncount)*(1.0_r8+0.608_r8*qm(ncount)))
          !zb      (ncount)  = cpdgrv * ts(ncount) * (1.0_r8 - bpeg(ncount,2) / bpeg(ncount,1) )!! boundary layer thickness (m)
          spdm    (ncount)  = spdm_sib(i) 
          dlwbot  (ncount)  = dlwbot_sib(i)
          cosz    (ncount)  = cosz_sib  (ncount) 
          gmt        (ncount,1)=tmtx(i,1,1)
          gmt        (ncount,2)=tmtx(i,1,2)
          gmt        (ncount,3)=tmtx(i,1,3)
          gmq        (ncount,1)=qmtx(i,1,1)
          gmq        (ncount,2)=qmtx(i,1,2)
          gmq        (ncount,3)=qmtx(i,1,3)
          gmu        (ncount,1)=umtx(i,1,1)
          gmu        (ncount,2)=umtx(i,1,2)
          gmu        (ncount,3)=umtx(i,1,3)
          gmu        (ncount,4)=umtx(i,1,4)
       END IF
    END DO
    IF(InitMod)THEN
       ncount=0
       DO i=1,nCols
       IF(iMask(i) >= 1_i8)THEN 
          ncount=ncount+1
          capacc(ncount,1,latco) = capacm(ncount,1,latco)
          capacc(ncount,2,latco) = capacm(ncount,2,latco)

          tcalc (ncount,latco) = tcm (ncount,latco)

          tgrdc (ncount,latco) =  tgm (ncount,latco)    

          tmm  (ncount,latco) = t_sib    (i,1)
          tmgc (ncount,latco) = t_sib    (i,1) 
          tm0  (ncount,latco) = t_sib    (i,1) 
          
          tcasm (ncount,latco) = t_sib         (i,1)
          tcasc (ncount,latco) = t_sib         (i,1)
          tcas0 (ncount,latco) = t_sib         (i,1)
          
          qcasm (ncount,latco) = sh_sib   (i,1)
          qcasc (ncount,latco) = sh_sib   (i,1)
          qcas0 (ncount,latco) = sh_sib   (i,1)
          
          
          qmm  (ncount,latco) = sh_sib   (i,1)
          qmgc (ncount,latco) = sh_sib   (i,1)
          qm0  (ncount,latco) = sh_sib   (i,1)
       END IF
       END DO
       DO k=1,nsoil
          ncount=0
          DO i=1,nCols
             IF(iMask(i) >=1_i8)THEN 
                ncount=ncount+1
                    td0 (ncount,k) = tdgm(ncount,k,latco)
                td  (ncount,k) = tdgm(ncount,k,latco)
                tdm (ncount,k) = tdgm(ncount,k,latco)
             END IF
          END DO
       END DO
    END IF
    
    ncount=0
    DO i=1,nCols
       IF(iMask(i) >=1_i8)THEN 
          ncount=ncount+1
          mlsi    (i) = 1_i8   !add solange 13-11-2012
          ta      (ncount)  = tcasc (ncount,latco)          
          qa      (ncount)  = qcasc (ncount,latco)
          tg      (ncount)  = tgrdc (ncount,latco) 
          tc      (ncount)  = tcalc (ncount,latco) 
          w       (ncount,1)= wwwgc (ncount,1,latco)
          w       (ncount,2)= wwwgc (ncount,2,latco)
          w       (ncount,3)= wwwgc (ncount,3,latco)
          capac   (ncount,1)= capacc(ncount,1,latco)
          capac   (ncount,2)= capacc(ncount,2,latco)

          rst     (ncount)  = stores(i,latco)
          snow    (ncount,1)= snowg (i,latco,1)
          snow    (ncount,2)= snowg (i,latco,2)       
          z0d     (ncount)  = 0.64_r8*timevar_zo (nCount,latco)
          zlt     (ncount)  = timevar_lai(nCount,latco)
          z1      (ncount)  = BioTab(iMask(i))%z1! surface parameters
          z2      (ncount)  = BioTab(iMask(i))%z2! surface parameters
          cc1     (ncount)  = timevar_rbc(nCount,latco)! surface parameters
          cc2     (ncount)  = timevar_rdc(nCount,latco)! surface parameters
          dd      (ncount)  = timevar_zp_disp(nCount,latco)! surface parameters
          !                                           soil moisture layers
          poros   (ncount)  =   bSoilGrd(i,latco)%poros!, &! surface parameters

          rootd   (ncount) = MIN(BioTab(iMask(i))%RootD,(BioTab(iMask(i))%SoDep*0.75_r8))
          SoDep   (ncount) = BioTab(iMask(i))%SoDep
          zdepth  (ncount,1) = 0.02_r8   
          zdepth  (ncount,2) = rootd(ncount) + 0.02_r8*SoDep (ncount)                      
          zdepth  (ncount,3) = SoDep(ncount) - zdepth(ncount,1) - zdepth(ncount,2)   
          

          !zdepth  (ncount,1) = 0.02   
          !zdepth  (ncount,2) = rootd - 0.02                    
          !zdepth  (ncount,3) = sodep - zdepth(ncount,1) - zdepth(ncount,2)   

          !zdepth  (ncount,1) = 0.02_r8 * poros(ncount)
          !zdepth  (ncount,2) = (rootd(ncount)  - 0.02_r8)*poros(ncount)
          !zdepth  (ncount,3) = poros(ncount)*BioTab(iMask(i))%SoDep - &
          !     ( zdepth(ncount,1)+zdepth(ncount,2) )

          phsat    (ncount)  = bSoilGrd(i,latco)%PhiSat! surface parameters
          bee      (ncount)  = bSoilGrd(i,latco)%BEE ! surface parameters
          respcp   (ncount)  = bSoilGrd(i,latco)%RespSat! surface parameters
          vmax0    (ncount)  = BioTab(iMask(i))%vmax0
          green    (ncount)  = timevar_green(nCount,latco)
          tran     (ncount,1:2,1:2) = BioTab(iMask(i))%LTran(1:2,1:2)
          ref      (ncount,1:2,1:2) = BioTab(iMask(i))%LRef (1:2,1:2)
          gmudmu   (ncount)  = timevar_gmudmu            (nCount,latco) 
          trop     (ncount)  = BioTab(iMask(i))%TROP
          phc      (ncount)  = BioTab(iMask(i))%Phi_half
          trda     (ncount)  = BioTab(iMask(i))%TRDA
          trdm     (ncount)  = BioTab(iMask(i))%TRDM
          slti     (ncount)  = BioTab(iMask(i))%SLTI
          shti     (ncount)  = BioTab(iMask(i))%SHTI 
          hltii    (ncount)  = BioTab(iMask(i))%HLTI
          hhti     (ncount)  = BioTab(iMask(i))%HHTI  
          effcon   (ncount)  = BioTab(iMask(i))%EffCon
          binter   (ncount)  = BioTab(iMask(i))%gsMin  
          gradm    (ncount)  = BioTab(iMask(i))%gsSlope
          atheta   (ncount)  = BioTab(iMask(i))%Atheta
          btheta   (ncount)  = BioTab(iMask(i))%Btheta
          aparc    (ncount)  = timevar_fpar(nCount,latco) 
          wopt     (ncount)  = bSoilGrd(i,latco)%Wopt    
          zm       (ncount)  = bSoilGrd(i,latco)%Skew    
          wsat     (ncount)  = bSoilGrd(i,latco)%RespSat  
          vcover   (ncount)  = BioTab(iMask(i))%fVCover
          SoRef    (ncount,1) = BioTab(iMask(i))%SoRef(1)
          SoRef    (ncount,2) = BioTab(iMask(i))%SoRef(2)
          radfac    (ncount,1:2,1:2,1:2)  = radfac_gbl     (ncount,1:2,1:2,1:2,latco)
          salb      (ncount,1:2,1:2)      = AlbGblSiB2     (ncount,1:2,1:2,latco)
          thermk    (ncount)              = thermk_gbl     (ncount,latco) 
          tgeff4_sib(ncount)              = tgeff4_sib_gbl (ncount,latco) 

          !IF(vcover(ncount) < zlt(ncount)/10.0_r8)THEN
          !   vcover(ncount) = vcover(ncount) * 10.0_r8
          !END IF
          pco2m  (ncount)  =  pco2m2  (i,latco) 
          satco  (ncount)  =  bSoilGrd(i,latco)%SatCo  
          slope  (ncount)  =  bSoilGrd(i,latco)%Slope 
          chil   (ncount)  =  BioTab(iMask(i))%ChiL
          sandfrac2(ncount)  =  sandfrac(i,latco)            
          clayfrac2(ncount)  =  clayfrac(i,latco)     
          IF(.NOT.forcerestore) THEN
             ztdep(ncount,1) = 6.0_r8 - sodep(ncount)                            !   0.10_r8
             ztdep(ncount,2) = sodep(ncount) - rootd(ncount)                    !   0.15_r8
             ztdep(ncount,3) = 8.0_r8 * (rootd(ncount)-0.02_r8) / 15.0_r8   !   0.25_r8
             ztdep(ncount,4) = 4. * (rootd(ncount)-0.02) / 15.0_r8            !   0.50_r8
             ztdep(ncount,5) = 2. * (rootd(ncount)-0.02) / 15.0_r8            !   1.00_r8
             ztdep(ncount,6) = 1. * (rootd(ncount)-0.02) / 15.0_r8            !   2.00_r8
          ELSE
             ztdep(ncount,1) = zdepth  (ncount,3)
             ztdep(ncount,2) = zdepth  (ncount,1)
          END IF
          IF(itype(ncount) == 13)THEN
              bee  (ncount)  =  2.9100_r8
              phsat(ncount)  = -0.758_r8
              satco(ncount)  =  9.216E-07_r8
              poros(ncount)  =  0.489_r8
              green(ncount)  =  0.1_r8              
              !zdepth  (ncount,1) = 0.02_r8 * poros(ncount)
              !zdepth  (ncount,2) = (rootd(ncount)  - 0.02_r8)*poros(ncount)
              !zdepth  (ncount,3) = poros(ncount)*BioTab(iMask(i))%SoDep - &
              ! ( zdepth(ncount,1)+zdepth(ncount,2) )
              zlt     (ncount)   = 0.0001_r8
              vcover  (ncount)   = 0.0001_r8
              chil    (ncount)   = 0.01_r8
              z1      (ncount)   = 0.0001_r8
              z2      (ncount)   = 0.1_r8
              z0d     (ncount)   = 0.01_r8
              dd      (ncount)   = 0.0004_r8
              cc1     (ncount)   = 35461.0_r8!timevar_rbc    rbc(ncount)=  
              cc2     (ncount)   = 28.5_r8   !timevar_rdc    rdc(ncount)=  
              !timetab_fpar    (ijmax)  ! Canopy absorbed fraction of PAR
              !timetab_lai     (ijmax)  ! Leaf-area index
              !timetab_green   (ijmax)  ! Canopy greeness fraction of LAI
              !timetab_zo      (ijmax)  ! Canopy roughness coeff
              !timetab_zp_disp (ijmax)  ! Zero plane displacement
              !timetab_rbc     (ijmax)  ! RB Coefficient (c1)
              !timetab_rdc     (ijmax)  ! RC Coefficient (c2)
              !timetab_gmudmu  (ijmax)  ! Time-mean leaf projection    

          END IF
!          ventmf(ncount)= ventmf2(i,latco) 
!          thvgm (ncount)= thvgm2 (i,latco) 
          pco2ap2 (ncount)= pco2ap   (i,latco)
          c4fract (ncount)= c4fractg (i,latco)
          d13cresp(ncount)= d13crespg(i,latco)
          d13cca  (ncount)= d13ccag  (i,latco)
       END IF
    END DO
    nsib=ncount

    DO k=1,nsoil
       ncount=0
       DO i=1,nCols
          IF(iMask(i) >=1_i8)THEN 
             ncount=ncount+1
             td0 (ncount,k) = tdg0(ncount,k,latco)
             td  (ncount,k) = tdgc(ncount,k,latco)
             tdm (ncount,k) = tdgm(ncount,k,latco)
          END IF
       END DO
    END DO

    !IF(InitMod)CALL vnqsat(1, tg(1:nsib), ps(1:nsib), qa(1:nsib), nsib ) 
    IF(InitMod)THEN
       nint=2
       IntSib=5
    ELSE
       nint=1
       IntSib=1
    END IF
    IF(TRIM(iswrad).NE.'NON'.AND.TRIM(ilwrad).NE.'NON') THEN
       IF(InitMod .and. nsib >= 1)THEN
       
          DO ind=1,nint
            ncount=0
             DO i=1,nCols
                IF(imask(i) >= 1_i8) THEN
                   ncount=ncount+1
                   !IF(ind.EQ.1) THEN
                   !   !
                   !   !     night
                   !   !
                   !   radvbc(ncount)    =0.0e0_r8
                   !   radnbc(ncount)    =0.0e0_r8
                   !   radvdc(ncount)    =0.0e0_r8
                   !   radndc(ncount)    =0.0e0_r8
                   !   cosz(ncount)      =0.0e0_r8
                   !ELSE
                      !
                      !     noon
                      !
                      radvbc(ncount) =beam_visb (i)
                      radnbc(ncount) =beam_nirb (i)
                      radvdc(ncount) =beam_visd (i)!
                      radndc(ncount) =beam_nird (i)
                      cosz(ncount)   =cos2(i)
                   !END IF
                   !
                   !     precipitation
                   !
                   lspr (ncount)=0.0e0_r8
                   cupr (ncount)=0.0e0_r8

                END IF
             END DO
             DO itr=1,IntSib
             
                CALL rada2(snow(1:nsib,1:2),zlt(1:nsib),z1(1:nsib),z2(1:nsib), &
                       asnow,tg(1:nsib),cosz(1:nsib),tice,ref(1:nsib,1:2,1:2),&
                       tran(1:nsib,1:2,1:2),chil(1:nsib) ,&
                       green(1:nsib),vcover(1:nsib),soref(1:nsib,1:2),radfac(1:nsib,1:2,1:2,1:2),&
                       salb(1:nsib,1:2,1:2),thermk(1:nsib),tgeff4_sib(1:nsib),&
                       tc(1:nsib),nsib)

                 CALL SIB(nsib                     , &! array and loop bounds  , INTENT(IN   ) :: nsib
                & nsib                            , &! array and loop bounds  , INTENT(IN   ) :: len
                & nsoil                           , &! array and loop bounds  , INTENT(IN   ) :: nsoil! array and loop bounds
                & tau                             , &! time independent variab, INTENT(IN   ) :: tau! time independent variable (hr)
                & ta       (1:nsib)           , &! CAS temperature (K) CAS, INTENT(INOUT) :: ta   (len)! CAS temperature (K)
                & qa       (1:nsib)           , &! CAS prognostic variable, INTENT(INOUT) :: sha  (len)! CAS water vapor mixing ratio (kg/kg)
                & tc       (1:nsib)           , &! prognostic variables   , INTENT(INOUT) :: tc   (len)! canopy (vegetation) temperature (K)
                & tg       (1:nsib)           , &! prognostic variables   , INTENT(INOUT) :: tg   (len)! surface boundary temperature (K)
                & td       (1:nsib,1:nsoil)   , &! prognostic variables   , INTENT(INOUT) :: td   (nsib,nsoil)! deep soil temperature (K)
                & w        (1:nsib,1:3)       , &! prognostic variables   , INTENT(INOUT) :: www  (nsib,3)! soil wetness 
                & snow     (1:nsib,1:2)       , &! prognostic variables   , INTENT(INOUT) :: snow (nsib,2)! snow cover ( kg /m^2)
                & capac    (1:nsib,1:2)       , &! prognostic variables   , INTENT(INOUT) :: capac(nsib,2)! liquid interception store (kg/m^2)
                & rst      (1:nsib)           , &! prognostic variables   , INTENT(INOUT) :: rst  (len)! stomatal resistance
                & tkes     (1:nsib)           , &! prognostic variables,INTENT(IN) :: tke(len)! turbulent kinetic energy
                & thm      (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: thm(len)! mixed layer potential temperature (K)
                & tm       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: thm(len)! mixed layer potential temperature (K)
                & qm       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: sh(len)! mixed layer water vapor mixing ratio (kg/kg)
                & um       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: thm(len)! mixed layer potential temperature (K)
                & vm       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: sh(len)! mixed layer water vapor mixing ratio (kg/kg)
                & sinclt   (1:nsib)           , &
                & ps       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: ps(len)! surface pressure (hPa)
                & bps      (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: bps(len)! (ps/1000)**kapa
                & rhoair   (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: ros (len)! surface air density (kg/m^3)
                & ts       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: ts  (len)! surface mixed layer air temperature
                & psb      (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: psb (len)! boundary layer mass depth (hPa)
                & cupr     (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: cupr(len)! convective precipitation rate (mm/s)
                & lspr     (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: lspr(len)! stratiform precipitation rate (mm/s)
                & radvbc   (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: radvbc(len)! surface incident visible direct beam (W/m^2)
                & radnbc   (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: radnbc(len)! surface incident near IR direct beam (W/m^2)
                & radvdc   (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: radvdc(len)! surface incident visible diffuse beam (W/m^2)
                & radndc   (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: radndc(len)! surface incident near IR diffuse beam (W/m^2)
                & zb       (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: zb(len)! boundary layer thickness (m)
                & spdm     (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: spdm(len)! boundary layer wind speed (m/s) 
                & dlwbot   (1:nsib)           , &! atmospheric forcing, INTENT(IN) :: dlwbot(len)! surface incident longwave radiation (W/m^2)
                & z0d      (1:nsib)           , &! surface parameters , INTENT(IN) :: z0d (len)! surface roughness length (m)
                & zlt      (1:nsib)           , &! surface parameters , INTENT(IN) :: zlt (len)! leaf area index
                & z1       (1:nsib)           , &! surface parameters , INTENT(IN) :: z1  (len)
                & z2       (1:nsib)           , &! surface parameters , INTENT(IN) :: z2  (len)
                & cc1      (1:nsib)           , &! surface parameters , INTENT(IN) :: cc1 (len)
                & cc2      (1:nsib)           , &! surface parameters , INTENT(IN) :: cc2 (len)
                & dd       (1:nsib)           , &! surface parameters , INTENT(IN) :: dd  (len)
                & zdepth   (1:nsib,1:3)       , &! surface parameters , INTENT(IN) :: zdepth(nsib,3) ! porosity * soil hydrology model layer depths (
                & poros    (1:nsib)           , &! surface parameters , INTENT(IN) :: poros(len)! soil porosity
                & phsat    (1:nsib)           , &! surface parameters , INTENT(IN) :: phsat(len)
                & bee      (1:nsib)           , &! surface parameters , INTENT(IN) :: bee(len)
                & respcp   (1:nsib)           , &! surface parameters , INTENT(IN) :: respcp(len)
                & vmax0    (1:nsib)           , &! surface parameters , INTENT(IN) :: vmax0(len)
                & green    (1:nsib)           , &! surface parameters , INTENT(IN) :: green(len)
                & tran     (1:nsib,1:2,1:2)   , &! surface parameters , INTENT(IN) :: tran(nsib,2,2)
                & ref      (1:nsib,1:2,1:2)   , &! surface parameters , INTENT(IN) :: ref(nsib,2,2)
                & gmudmu   (1:nsib)           , &! surface parameters , INTENT(IN) :: gmudmu(len)
                & trop     (1:nsib)           , &! surface parameters , INTENT(IN) :: trop(len)
                & phc      (1:nsib)           , &! surface parameters , INTENT(IN) :: phc(len)
                & trda     (1:nsib)           , &! surface parameters , INTENT(IN) :: trda(len)
                & trdm     (1:nsib)           , &! surface parameters , INTENT(IN) :: trdm(len)
                & slti     (1:nsib)           , &! surface parameters , INTENT(IN) :: slti(len)
                & shti     (1:nsib)           , &! surface parameters , INTENT(IN) :: shti(len)
                & hltii    (1:nsib)           , &! surface parameters , INTENT(IN) :: hltii(len)
                & hhti     (1:nsib)           , &! surface parameters , INTENT(IN) :: hhti(len)
                & effcon   (1:nsib)           , &! surface parameters , INTENT(IN) :: effcon(len)
                & binter   (1:nsib)           , &! surface parameters , INTENT(IN) :: binter(len)
                & gradm    (1:nsib)           , &! surface parameters , INTENT(IN) :: gradm(len)
                & atheta   (1:nsib)           , &! surface parameters , INTENT(IN) :: atheta(len)
                & btheta   (1:nsib)           , &! surface parameters , INTENT(IN) :: btheta(len)
                & aparc    (1:nsib)           , &! surface parameters , INTENT(IN) :: aparc(len)
                & wopt     (1:nsib)           , &! surface parameters , INTENT(IN) :: wopt(len)
                & zm       (1:nsib)           , &! surface parameters , INTENT(IN) :: zm(len)
                & wsat     (1:nsib)           , &! surface parameters , INTENT(IN) :: wsat(len)
                & vcover   (1:nsib)           , &! surface parameters , INTENT(IN) :: vcover(len)   ! vegetation cover fraction
                & radfac   (1:nsib,1:2,1:2,1:2), &! surface paramet, INTENT(IN):: radfac(nsib,2,2,2)
                & thermk   (1:nsib)           , &! surface parameters , INTENT(IN) :: thermk(len)
                & satco    (1:nsib)           , &! surface parameters , INTENT(IN) :: satco(len)
                & slope    (1:nsib)           , &! surface parameters , INTENT(IN) :: slope(len)
                & chil     (1:nsib)           , &! surface parameters , INTENT(IN) :: chil(len)
                & ztdep    (1:nsib,1:nSoil)   , &! surface parameters , INTENT(IN) :: ztdep(nsib,nsoil) ! soil thermal model layer depths (m)
                & sandfrac2(1:nsib)           , &! surface parameters , INTENT(IN) :: sandfrac(len)
                & clayfrac2(1:nsib)           , &! surface parameters , INTENT(IN) :: soil texture clay fraction
                & pi                          , &! constants , INTENT(IN  ) :: pi ! 3.1415926....
                & grav                        , &! constants , INTENT(IN  ) :: grav! gravitational acceleration (m/s^2)
                & dtc3x                       , &! constants , INTENT(IN  ) :: dt! time step (s)
                & cp                          , &! constants , INTENT(IN  ) :: cp! heat capacity of dry air (J/(kg K) )
                & cv                          , &! constants , INTENT(IN  ) :: cv ! heat capacity of water vapor (J/(kg K))
                & hltm                        , &! constants , INTENT(IN  ) :: hltm! latent heat of vaporization (J/kg)
                & delta                       , &! constants , INTENT(IN  ) :: delta
                & asnow                       , &! constants , INTENT(IN  ) :: asnow
                & kapa                        , &! constants , INTENT(IN  ) :: kapa
                & snomel                      , &! constants , INTENT(IN  ) :: snomel
                & clai                        , &! constants , INTENT(IN  ) :: clai   ! leaf heat capacity
                & cww                         , &! constants , INTENT(IN  ) :: cww! water heat capacity
                & pr0                         , &! constants , INTENT(IN  ) :: pr0
                & ribc                        , &! constants , INTENT(IN  ) :: ribc
                & vkrmn                       , &! constants , INTENT(IN  ) :: vkrmn! von karmann's constant
                & pco2m    (1:nsib)           , &! constants , INTENT(IN  ) :: pco2m(len)! CO2 partial pressure (Pa)
                & po2m                        , &! constants , INTENT(IN  ) :: po2m! O2 partial pressure (Pa)
                & stefan                      , &! constants , INTENT(IN  ) :: stefan! stefan-boltzman constant
                & grav2                       , &! constants , INTENT(IN  ) :: grav2! grav / 100 Pa/hPa
                & tice                        , &! constants , INTENT(IN  ) :: tice! freezing temperature (KP
                & gmt      (1:nsib,:)         , &
                & gmq      (1:nsib,:)         , &
                & gmu      (1:nsib,:)         , &
                & snofac                      , &! constants , INTENT(IN) :: snofac
                & fss      (1:nsib)           , &! INTENT(OUT)::fss(len)! surface sensible heat flux (W/m^2)
                & fws      (1:nsib)           , &! INTENT(OUT)::fws(len)! surface evaporation (kg/m^2/s)
                & drag     (1:nsib,1:2)       , &! INTENT(OUT)::drag(nsib,2)! surface drag coefficient
                & co2flx   (1:nsib)           , &! INTENT(OUT)::co2flx(len)! surface CO2 flux
                & cu       (1:nsib)           , &! INTENT(OUT)::cu(len)
                & ct       (1:nsib)           , &! INTENT(OUT)::ct(len)
                & ustar    (1:nsib)           , &! INTENT(OUT)::ustar (len)! friction velocity (m/s)
                & hr       (1:nsib)           , &! INTENT(OUT)
                & ect      (1:nsib)           , &! INTENT(OUT)
                & eci      (1:nsib)           , &! INTENT(OUT)
                & egt      (1:nsib)           , &! INTENT(OUT)
                & egi      (1:nsib)           , &! INTENT(OUT)
                & egs      (1:nsib)           , &! INTENT(OUT)
                & eg       (1:nsib)           , &! INTENT(OUT)
                & ec       (1:nsib)           , &! INTENT(OUT)
                & es       (1:nsib)           , &! INTENT(OUT)
                & hc       (1:nsib)           , &! INTENT(OUT)
                & hg       (1:nsib)           , &! INTENT(OUT)
                & hs       (1:nsib)           , &! INTENT(OUT)
                & chf      (1:nsib)           , &! INTENT(OUT)
                & shf      (1:nsib)           , &! INTENT(OUT)
                & roff     (1:nsib)           , &! INTENT(OUT)
                & ra       (1:nsib)           , &! INTENT(OUT)
                & rb       (1:nsib)           , &! INTENT(OUT)
                & rd       (1:nsib)           , &! INTENT(OUT)
                & rc       (1:nsib)           , &! INTENT(OUT)
                & rg       (1:nsib)           , &! INTENT(OUT)
                & ea       (1:nsib)           , &! INTENT(OUT)
                & etc      (1:nsib)           , &! INTENT(OUT)
                & etg      (1:nsib)           , &! INTENT(OUT)
                & radt     (1:nsib,1:3)       , &! INTENT(OUT)
                & rsoil    (1:nsib)           , &! INTENT(OUT)
                & ioffset                     , &! logical switches
                & sibdrv                      , &! logical switches
                & forcerestore                , &! logical switches
                & doSiBco2                    , &! logical switches
                & fixday                      , &! logical switches
                & dotkef                      , &! logical switches
                & ventmf   (1:nsib)           , &! logical switches INTENT(OUT) :: ventmf(len) ! ventilation mass flux
                & thvgm    (1:nsib)           , &! logical switches INTENT(OUT) :: thvgm(len)
                & louis                       , &! logical switches
                & respfactor(1:nsib,1:nsoil+1), &! logical switches
                & xgpp      (1:nsib)          , &! logical switches
                & pco2ap2   (1:nsib)          , &! logical switches  INTENT(IN OUT) :: pco2ap(len)! canopy air space pCO2 (Pa)
                & c4fract   (1:nsib)          , &! logical switches  INTENT(IN OUT) :: c4fract(len)! fraction of C4 vegetation  
                & d13cresp  (1:nsib)          , &! logical switches  INTENT(IN OUT) :: d13cresp(len)! del 13C of respiration (per mil vs PDB)
                & d13cca    (1:nsib)          , &! logical switches  INTENT(IN OUT) :: d13cca(len)! del 13C of canopy CO2 (per mil vs PDB)
                & isotope                     , &! logical switches
                & itype     (1:nsib)          , &
                & delsig    (1:nsib)             , &
                & sigki     (1:nsib)               , &
                & bstar1    (1:nsib)          , &
                & zlwup     (1:nsib)          , &
                & areas                         )
                ncount=0
                DO i=1,nCols
                   IF(iMask(i) >= 1_i8) THEN
                      ncount=ncount+1
                      tm (ncount  )= t_sib  (i,1)
                      qm (ncount  )= sh_sib (i,1)
                      gmt(ncount,1)=tmtx(i,1,1)
                      gmt(ncount,2)=tmtx(i,1,2)
                      gmt(ncount,3)=tmtx(i,1,3)
                      gmq(ncount,1)=qmtx(i,1,1)
                      gmq(ncount,2)=qmtx(i,1,2)
                      gmq(ncount,3)=qmtx(i,1,3)
                      gmu(ncount,1)=umtx(i,1,1)
                      gmu(ncount,2)=umtx(i,1,2)
                      gmu(ncount,3)=umtx(i,1,3)
                      gmu(ncount,4)=umtx(i,1,4)
                   END IF
                END DO
             END DO
             DO k=1,nsoil
                DO i=1,nsib
                     td   (i,k)  =tdm   (i,k)
                END DO
             END DO
             DO i=1,nsib
                capac(i,1)=capacm(i,1,latco)
                capac(i,2)=capacm(i,2,latco)
                w    (i,1)=wm    (i,1,latco)
                w    (i,2)=wm    (i,2,latco)
                w    (i,3)=wm    (i,3,latco)
                tc   (i)  =tcm   (i,latco)
                ta   (i)  =tcasm (i,latco)
                qa   (i)  =qcasm(i,latco)
                IF(ind.EQ.1) THEN
                   tmin (i) =tg (i)
                ELSE
                   tmax (i) =tg (i)
                END IF
                tg   (i) =tgm(i,latco)
             END DO             
          END DO
          DO k=1,nsoil
             DO i=1,nsib
                  td   (i,k) =tdm   (i,k)
                  td   (i,k) =0.9_r8*0.5_r8*(tmax(i)+tmin(i))+0.1_r8*tdm(i,k)
                  tdm  (i,k) =td(i,k)
                  td0  (i,k) =td(i,k)
             END DO
          END DO
          !
          !     this is a start of equilibrium tg,tc comp.
          !
          ncount=0
          DO i=1,nCols
             IF(imask(i) >= 1_i8) THEN
                ncount=ncount+1
                cosz(ncount)    =zenith(i)
             END IF
          END DO

          DO i=1,nsib
             IF(cosz(i).LT.0.0e0_r8) THEN
                tgm  (i,latco)  =tmin(i)
                tg0  (i,latco)  =tmin(i)
             END IF
          END DO
          CALL rada2(snow(1:nsib,1:2),zlt(1:nsib),z1(1:nsib),z2(1:nsib), &
                   asnow,tg(1:nsib),cosz(1:nsib),tice,ref(1:nsib,1:2,1:2),&
                   tran(1:nsib,1:2,1:2),chil(1:nsib) ,&
                   green(1:nsib),vcover(1:nsib),soref(1:nsib,1:2),radfac(1:nsib,1:2,1:2,1:2),&
                   salb(1:nsib,1:2,1:2),thermk(1:nsib),tgeff4_sib(1:nsib),&
                   tc(1:nsib),nsib)
       END IF
    END IF
    IF(nsib.GE.1) THEN
       ncount=0
       DO i=1,nCols
          IF(imask(i) >= 1_i8) THEN
             ncount=ncount+1
             !
             !     this is for radiation interpolation
             !
             !IF(cosz(ncount).GE.0.01746e0_r8 ) THEN
                radvbc(ncount)  = radvbc_sib(i)
                radnbc(ncount)  = radnbc_sib(i)
                radvdc(ncount)  = radvdc_sib(i)
                radndc(ncount)  = radndc_sib(i)
             !ELSE
             !   radvbc(ncount)=0.0e0_r8
             !   radnbc(ncount)=0.0e0_r8
             !   radvdc(ncount)=0.0e0_r8
             !   radndc(ncount)=0.0e0_r8
             !END IF
             !
             !     precipitation
             !
             lspr (ncount)=dlspr_sib (i)/(tau) !convert mm --> mm/s ! tau 2 DT
             cupr (ncount)=dcupr_sib (i)/(tau) !convert mm --> mm/s 
          END IF
       END DO

       CALL SIB(nsib               , &! array and loop bounds  , INTENT(IN   ) :: nsib
      & nsib                      , &! array and loop bounds  , INTENT(IN   ) :: len
        nsoil                     , &! array and loop bounds  , INTENT(IN   ) :: nsoil         ! array and loop bounds
        tau                       , &! time independent variab, INTENT(IN   ) :: tau         ! time independent variable (hr)
      & ta       (1:nsib)         , &! CAS temperature (K) CAS, INTENT(INOUT) :: ta   (len)! CAS temperature (K)
        qa       (1:nsib)         , &! CAS prognostic variable, INTENT(INOUT) :: sha  (len)! CAS water vapor mixing ratio (kg/kg)
        tc       (1:nsib)         , &! prognostic variables   , INTENT(INOUT) :: tc   (len)    ! canopy (vegetation) temperature (K)
        tg       (1:nsib)         , &! prognostic variables   , INTENT(INOUT) :: tg   (len)           ! surface boundary temperature (K)
        td       (1:nsib,1:nsoil) , &! prognostic variables   , INTENT(INOUT) :: td   (nsib,nsoil)  ! deep soil temperature (K)
        w        (1:nsib,1:3)     , &! prognostic variables   , INTENT(INOUT) :: www  (nsib,3)           ! soil wetness 
        snow     (1:nsib,1:2)     , &! prognostic variables   , INTENT(INOUT) :: snow (nsib,2)           ! snow cover ( kg /m^2)
        capac    (1:nsib,1:2)     , &! prognostic variables   , INTENT(INOUT) :: capac(nsib,2)  ! liquid interception store (kg/m^2)
        rst      (1:nsib)         , &! prognostic variables   , INTENT(INOUT) :: rst  (len)     ! stomatal resistance
        tkes     (1:nsib)         , &! prognostic variables,INTENT(IN   ) :: tke(len)  ! turbulent kinetic energy
        thm      (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: thm(len)           ! mixed layer potential temperature (K)
        tm       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: thm(len)           ! mixed layer potential temperature (K)
        qm       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: sh(len)           ! mixed layer water vapor mixing ratio (kg/kg)
        um       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: thm(len)           ! mixed layer potential temperature (K)
        vm       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: sh(len)           ! mixed layer water vapor mixing ratio (kg/kg)
        sinclt   (1:nsib)         , &
        ps       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: ps(len)           ! surface pressure (hPa)
        bps      (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: bps(len)           ! (ps/1000)**kapa
        rhoair   (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: ros (len)           ! surface air density (kg/m^3)
        ts       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: ts  (len)           ! surface mixed layer air temperature
        psb      (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: psb (len)           ! boundary layer mass depth (hPa)
        cupr     (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: cupr(len)           ! convective precipitation rate (mm/s)
        lspr     (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: lspr(len)           ! stratiform precipitation rate (mm/s)
        radvbc   (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: radvbc(len)    ! surface incident visible direct beam (W/m^2)
        radnbc   (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: radnbc(len)    ! surface incident near IR direct beam (W/m^2)
        radvdc   (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: radvdc(len)    ! surface incident visible diffuse beam (W/m^2)
        radndc   (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: radndc(len)    ! surface incident near IR diffuse beam (W/m^2)
        zb       (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: zb(len)           ! boundary layer thickness (m)
        spdm     (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: spdm(len)           ! boundary layer wind speed (m/s) 
        dlwbot   (1:nsib)         , &! atmospheric forcing, INTENT(IN   ) :: dlwbot(len)    ! surface incident longwave radiation (W/m^2)
        z0d      (1:nsib)         , &! surface parameters , INTENT(IN   ) :: z0d (len)           ! surface roughness length (m)
        zlt      (1:nsib)         , &! surface parameters , INTENT(IN   ) :: zlt (len)           ! leaf area index
        z1       (1:nsib)         , &! surface parameters , INTENT(IN   ) :: z1  (len)
        z2       (1:nsib)         , &! surface parameters , INTENT(IN   ) :: z2  (len)
        cc1      (1:nsib)         , &! surface parameters , INTENT(IN   ) :: cc1 (len)
        cc2      (1:nsib)         , &! surface parameters , INTENT(IN   ) :: cc2 (len)
        dd       (1:nsib)         , &! surface parameters , INTENT(IN   ) :: dd  (len)
        zdepth   (1:nsib,1:3)     , &! surface parameters , INTENT(IN   ) :: zdepth(nsib,3) ! porosity * soil hydrology model layer depths (
        poros    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: poros(len)     ! soil porosity
        phsat    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: phsat(len)
        bee      (1:nsib)         , &! surface parameters , INTENT(IN   ) :: bee(len)
        respcp   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: respcp(len)
        vmax0    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: vmax0(len)
        green    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: green(len)
        tran     (1:nsib,1:2,1:2) , &! surface parameters , INTENT(IN   ) :: tran(nsib,2,2)
        ref      (1:nsib,1:2,1:2) , &! surface parameters , INTENT(IN   ) :: ref(nsib,2,2)
        gmudmu   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: gmudmu(len)
        trop     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: trop(len)
        phc      (1:nsib)         , &! surface parameters , INTENT(IN   ) :: phc(len)
        trda     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: trda(len)
        trdm     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: trdm(len)
        slti     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: slti(len)
        shti     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: shti(len)
        hltii    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: hltii(len)
        hhti     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: hhti(len)
        effcon   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: effcon(len)
        binter   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: binter(len)
        gradm    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: gradm(len)
        atheta   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: atheta(len)
        btheta   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: btheta(len)
        aparc    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: aparc(len)
        wopt     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: wopt(len)
        zm       (1:nsib)         , &! surface parameters , INTENT(IN   ) :: zm(len)
        wsat     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: wsat(len)
        vcover   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: vcover(len)   ! vegetation cover fraction
        radfac   (1:nsib,1:2,1:2,1:2), &! surface paramet, INTENT(IN   ) :: radfac(nsib,2,2,2)
        thermk   (1:nsib)         , &! surface parameters , INTENT(IN   ) :: thermk(len)
        satco    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: satco(len)
        slope    (1:nsib)         , &! surface parameters , INTENT(IN   ) :: slope(len)
        chil     (1:nsib)         , &! surface parameters , INTENT(IN   ) :: chil(len)
        ztdep    (1:nsib,1:nSoil) , &! surface parameters , INTENT(IN   ) :: ztdep(nsib,nsoil) ! soil thermal model layer depths (m)
        sandfrac2(1:nsib)         , &! surface parameters , INTENT(IN   ) :: sandfrac(len)     ! soil texture sand fraction
        clayfrac2(1:nsib)         , &! surface parameters , INTENT(IN   ) :: clayfrac(len)     ! soil texture sand fraction
        pi                        , &! constants , INTENT(IN   ) :: pi         ! 3.1415926....
        grav                      , &! constants , INTENT(IN   ) :: grav          ! gravitational acceleration (m/s^2)
        dtc3x                     , &! constants , INTENT(IN   ) :: dt         ! time step (s)
        cp                        , &! constants , INTENT(IN   ) :: cp         ! heat capacity of dry air (J/(kg K) )
        cv                        , &! constants , INTENT(IN   ) :: cv         ! heat capacity of water vapor (J/(kg K))
        hltm                      , &! constants , INTENT(IN   ) :: hltm          ! latent heat of vaporization (J/kg)
        delta                     , &! constants , INTENT(IN   ) :: delta
        asnow                     , &! constants , INTENT(IN   ) :: asnow
        kapa                      , &! constants , INTENT(IN   ) :: kapa
        snomel                    , &! constants , INTENT(IN   ) :: snomel
        clai                      , &! constants , INTENT(IN   ) :: clai          ! leaf heat capacity
        cww                       , &! constants , INTENT(IN   ) :: cww   ! water heat capacity
        pr0                       , &! constants , INTENT(IN   ) :: pr0
        ribc                      , &! constants , INTENT(IN   ) :: ribc
        vkrmn                     , &! constants , INTENT(IN   ) :: vkrmn         ! von karmann's constant
        pco2m    (1:nsib)         , &! constants , INTENT(IN   ) :: pco2m(len)! CO2 partial pressure (Pa)
        po2m                      , &! constants , INTENT(IN   ) :: po2m          ! O2 partial pressure (Pa)
        stefan                    , &! constants , INTENT(IN   ) :: stefan    ! stefan-boltzman constant
        grav2                     , &! constants , INTENT(IN   ) :: grav2         ! grav / 100 Pa/hPa
        tice                      , &! constants , INTENT(IN   ) :: tice          ! freezing temperature (KP
        gmt      (1:nsib,:)       , &
        gmq      (1:nsib,:)       , &
        gmu      (1:nsib,:)       , &
        snofac                    , &! constants , INTENT(IN   ) :: snofac
        fss      (1:nsib)         , &! output ,INTENT(OUT   ) :: fss   (len)   ! surface sensible heat flux (W/m^2)
        fws      (1:nsib)         , &! output ,INTENT(OUT   ) :: fws   (len)   ! surface evaporation (kg/m^2/s)
        drag     (1:nsib,1:2)     , &! output ,INTENT(OUT   ) :: drag  (nsib,2)! surface drag coefficient
        co2flx   (1:nsib)         , &! output ,INTENT(OUT   ) :: co2flx(len)   ! surface CO2 flux
        cu       (1:nsib)         , &! output ,INTENT(OUT   ) :: cu    (len)
        ct       (1:nsib)         , &! output ,INTENT(OUT   ) :: ct    (len)
        ustar    (1:nsib)         , &! output ,INTENT(OUT   ) :: ustar (len)   ! friction velocity (m/s)
        hr       (1:nsib)         , &! output ,INTENT(OUT   ) 
        ect      (1:nsib)         , &! output ,INTENT(OUT   ) 
        eci      (1:nsib)         , &! output ,INTENT(OUT   ) 
        egt      (1:nsib)         , &! output ,INTENT(OUT   ) 
        egi      (1:nsib)         , &! output ,INTENT(OUT   ) 
        egs      (1:nsib)         , &! output ,INTENT(OUT   ) 
        eg       (1:nsib)         , &! output ,INTENT(OUT   ) 
        ec       (1:nsib)         , &! output ,INTENT(OUT   ) 
        es       (1:nsib)         , &! output ,INTENT(OUT   )  
        hc       (1:nsib)         , &! output ,INTENT(OUT   ) 
        hg       (1:nsib)         , &! output ,INTENT(OUT   ) 
        hs       (1:nsib)         , &! output ,INTENT(OUT   ) 
        chf      (1:nsib)         , &! output ,INTENT(OUT   ) 
        shf      (1:nsib)         , &! output ,INTENT(OUT   ) 
        roff     (1:nsib)         , &! output ,INTENT(OUT   ) 
        ra       (1:nsib)         , &! output ,INTENT(OUT   ) 
        rb       (1:nsib)         , &! output ,INTENT(OUT   ) 
        rd       (1:nsib)         , &! output ,INTENT(OUT   ) 
        rc       (1:nsib)         , &! output ,INTENT(OUT   ) 
        rg       (1:nsib)         , &! output ,INTENT(OUT   ) 
        ea       (1:nsib)         , &! output ,INTENT(OUT   ) 
        etc      (1:nsib)         , &! output ,INTENT(OUT   ) 
        etg      (1:nsib)         , &! output ,INTENT(OUT   ) 
        radt     (1:nsib,1:3)     , &! output ,INTENT(OUT   ) 
        rsoil    (1:nsib)         , &! output ,INTENT(OUT   ) 
        ioffset                   , &! logical switches
        sibdrv                    , &! logical switches
        forcerestore              , &! logical switches
        doSiBco2                  , &! logical switches
        fixday                    , &! logical switches
        dotkef                    , &! logical switches
        ventmf   (1:nsib)         , &! logical switches INTENT(OUT) :: ventmf(len) ! ventilation mass flux
        thvgm    (1:nsib)         , &! logical switches INTENT(OUT) :: thvgm(len)
        louis                     , &! logical switches
        respfactor(1:nsib,1:nsoil+1), &! logical switches    
        xgpp      (1:nsib)        , &! logical switches        
        pco2ap2   (1:nsib)        , &! logical switches  INTENT(IN OUT) :: pco2ap(len)         ! canopy air space pCO2 (Pa)
        c4fract   (1:nsib)        , &! logical switches  INTENT(IN OUT) :: c4fract(len)   ! fraction of C4 vegetation  
        d13cresp  (1:nsib)        , &! logical switches  INTENT(IN OUT) :: d13cresp(len)  ! del 13C of respiration (per mil vs PDB)
        d13cca    (1:nsib)        , &! logical switches  INTENT(IN OUT) :: d13cca(len)         ! del 13C of canopy CO2 (per mil vs PDB)
        isotope                   , &! logical switches
        itype     (1:nsib)        , &
        delsig    (1:nsib)           , &
        sigki     (1:nsib)             , &
        bstar1    (1:nsib)        , &
        zlwup     (1:nsib)        , &
        areas                       )


   !       CALL rada2(snow(1:nsib,1:2),zlt(1:nsib),z1(1:nsib),z2(1:nsib), &
   !                asnow,tg(1:nsib),cosz(1:nsib),tice,ref(1:nsib,1:2,1:2),&
   !                tran(1:nsib,1:2,1:2),chil(1:nsib) ,&
   !                green(1:nsib),vcover(1:nsib),soref(1:nsib,1:2),radfac(1:nsib,1:2,1:2,1:2),&
   !                salb(1:nsib,1:2,1:2),thermk(1:nsib),tgeff4_sib(1:nsib),&
   !                tc(1:nsib),nsib)

    END IF
    !     
    !     temperature and snow depths in Antarctica and Groenland
    !     
    ncount=0
    DO i=1,nCols
       IF(imask(i).GE.1_i8 ) THEN
          ncount=ncount+1
          IF ( imask(i).EQ.13_i8 ) THEN
            w (ncount,1)  = 1.0_r8
            w (ncount,2)  = 1.0_r8
            w (ncount,3)  = 1.0_r8
            DO k=1,nsoil
               TD(ncount ,k )  = MAX(MIN(TD(ncount,k) ,273.15_r8),218.15_r8)
            END DO 
            TC(ncount  )  = MAX(MIN(TC(ncount) ,273.15_r8),218.15_r8)
            tg(ncount  )  = MAX(MIN(tg(ncount) ,273.15_r8),218.15_r8)
          END IF
       END IF
    END DO

    !
    !     sib time integaration and time filter
    !
    DO i=1,nsib
       qm(i)=MAX(1.0e-12_r8,qm(i)) 
       !tm(i)=tm(i) / bps     (i) 
    END DO

    CALL sextrp &
         (td     ,qa     ,ta     ,tg     , &
         tc     ,w      ,capac  ,td0    , &
         qcas0(:,latco)    ,tcas0 (:,latco)   ,tg0(:,latco)    ,tc0(:,latco)    , &
         w0(:,:,latco)     ,capac0(:,:,latco) ,tdm    ,qcasm (:,latco)   , &
         tcasm(:,latco)    ,tgm(:,latco)    ,tcm(:,latco)   ,wm(:,:,latco)     , &
         capacm(:,:,latco) ,istrt  ,ncols  ,nsib   , &
         epsflt ,intg      ,tm0(:,latco)    , &
         qm0(:,latco)    ,tm     ,qm     ,tmm(:,latco)    , &
         qmm(:,latco)    ,nsoil    )


    !
    !     fix soil moisture at selected locations
    !
!    DO i=1,nsib
!       IF(ssib(i).GT.0.0_r8)THEN
!          qm(i)=MAX(1.0e-12_r8,qm(i))
!          w0(i,1)=ssib(i)
!          w0(i,2)=ssib(i)
!          w0(i,3)=ssib(i)
!          wm(i,1)=ssib(i)
!          wm(i,2)=ssib(i)
!          wm(i,3)=ssib(i)
!       END IF
!    END DO
    ncount=0
    DO i=1,nCols
       IF(imask(i) >=1_i8) THEN
          ncount=ncount+1
          tmtx(i,1,3)=gmt(ncount,3)
          qmtx(i,1,3)=gmq(ncount,3)
          umtx(i,1,3)=gmu(ncount,3)
          umtx(i,1,4)=gmu(ncount,4)
          tsea(i)    =tgeff4_sib(ncount)!/ bps     (ncount) 
          IF(omlmodel)TSK (i)=tgeff4_sib(ncount)
       END IF
    END DO

    ncount=0
    DO i=1,nCols
       IF(iMask(i) >=1_i8)THEN 
          ncount=ncount+1
          co2fx    (i,latco) = co2flx  (ncount)
          pco2ap   (i,latco) = pco2ap2 (ncount)  
          c4fractg (i,latco) = c4fract (ncount)  
          d13crespg(i,latco) = d13cresp(ncount)  
          d13ccag  (i,latco) = d13cca  (ncount)  
          tmgc     (ncount,latco) = tm      (ncount)
          qmgc     (ncount,latco)  = qm      (ncount)
          tcasc    (ncount,latco) = ta      (ncount  )  
          qcasc    (ncount,latco) = qa      (ncount  ) 
          tgrdc  (ncount,latco  ) = tg      (ncount  )
          tcalc  (ncount,latco  ) = tc      (ncount  ) 
          stores (i,latco  ) = rst     (ncount  )  
          snowg  (i,latco,1) = snow    (ncount,1)
          snowg  (i,latco,2) = snow    (ncount,2)         
          wwwgc  (ncount,1,latco) = w       (ncount,1) 
          wwwgc  (ncount,2,latco) = w       (ncount,2) 
          wwwgc  (ncount,3,latco) = w       (ncount,3) 
          w0     (ncount,1,latco) = w       (ncount,1) 
          w0     (ncount,2,latco) = w       (ncount,2) 
          w0     (ncount,3,latco) = w       (ncount,3) 
          wm     (ncount,1,latco) = w       (ncount,1) 
          wm     (ncount,2,latco) = w       (ncount,2) 
          wm     (ncount,3,latco) = w       (ncount,3) 
          
          tm0     (ncount,latco) = tm            (ncount) 
          tmm     (ncount,latco) = tm            (ncount) 

          qm0     (ncount,latco) = qm            (ncount) 
          qmm     (ncount,latco) = qm            (ncount) 
          tg0     (ncount,latco) = tg            (ncount) 

          tgm     (ncount,latco) = tg       (ncount) 
          tc0     (ncount,latco) = tc       (ncount) 
          tcm     (ncount,latco) = tc       (ncount) 

          capacc (ncount,1,latco) = capac  (ncount,1) 
          capacc (ncount,2,latco) = capac  (ncount,2) 
          capac0 (ncount,1,latco) = capac  (ncount,1) 
          capac0 (ncount,2,latco) = capac  (ncount,2) 
          capacm (ncount,1,latco) = capac  (ncount,1) 
          capacm (ncount,2,latco) = capac  (ncount,2) 
       END IF
    END DO
    ncount=0
    DO i=1,nCols
       IF(imask(i) >= 1_i8 ) THEN
          ncount=ncount+1
          IF ( imask(i) == 13_i8 ) THEN
             sm0 (ncount,1)    = poros (ncount)
          ELSE
             sm0 (ncount,1)    = w0(ncount,1,latco)* poros (ncount)
             sm0 (ncount,2)    = w0(ncount,2,latco)* poros (ncount)
             sm0 (ncount,3)    = w0(ncount,3,latco)* poros (ncount)
          END IF
       END IF
    END DO
    DO k=1,nsoil
       ncount=0
       DO i=1,nCols
          IF(iMask(i) >=1_i8)THEN 
             ncount=ncount+1
             tdgm(ncount,k,latco) = td  (ncount,k)
             tdm(ncount,k) = td  (ncount,1)
             td0(ncount,k) = td  (ncount,1)
             tdgc(ncount,k,latco) = td  (ncount,k)
             td1 (ncount)         = td  (ncount,1)
             tdg0(ncount,k,latco) = td  (ncount,k)
          END IF
       END DO
    END DO
    !
    !     sea or sea ice
    ! gu gv gps colrad sigki delsig sens evap umom vmom rmi rhi cond stor zorl rnet ztn2 THETA_2M VELC_2m MIXQ_2M
    ! THETA_10M VELC_10M MIXQ_10M
    ! mmax=ncols-nmax+1
    ! including case 1D physics

!    IF(InitMod)THEN
!      tsfc0(1:nCols,latco)=t_sib  (1:nCols)
!      QSfc0SiB2(1:nCols,latco)=sh_sib (1:nCols)
!      TSfcmSiB2(1:nCols,latco)=t_sib  (1:nCols)
!      QSfcmSib2(1:nCols,latco)=sh_sib (1:nCols)
!      tsfc (1:nCols)=t_sib  (1:nCols)
!      qsfc (1:nCols)=sh_sib (1:nCols)
!    END IF
    DO i=1,nCols
       xndvi(i) = ndvim(i)
       IF(mskant(i) == 1_i8)THEN
          xsea (i) = tseam(i)
          tsfc (i) = TSfcm(i,latco)
          qsfc (i) = QSfcm(i,latco)
       END IF   
    END DO
   CALL seasfc( &
           tmtx  (1:nCols,1:kMax,1:3)  ,umtx  (1:nCols,1:kMax,1:4),qmtx  (1:nCols,1:kMax,1:3)  ,&
           kmax                        ,kmax                      ,slrad (1:nCols)             ,&
           tsurf (1:nCols)             ,qsurf (1:nCols)           ,gu    (1:nCols,1:kMax)      ,&
           gv    (1:nCols,1:kMax)      ,t_sib (1:nCols,1:kMax)    ,sh_sib(1:nCols,1:kMax)      ,&
           prsi(1:ncols,1:kMax+1),prsl(1:ncols,1:kMax),phii(1:nCols,1:kMax+1),phil(1:nCols,1:kMax),&
           xsea  (1:nCols)           ,dtc3x                       ,SIN(colrad(1:nCols))        ,&
           sensf (1:nCols,latco)       ,latef (1:nCols,latco)     ,umom  (1:nCols)             ,&
           vmom  (1:nCols)             ,rmi   (1:nCols)           ,rhi   (1:nCols)             ,&
           cond  (1:nCols)             ,stor  (1:nCols)           ,zorl  (1:nCols)             ,&
           nCols                       ,spdm_sib(1:nCols)         ,bstar (1:nCols)             ,&
           Ustarm(1:nCols)             ,z0sea (1:nCols)           ,rho   (1:nCols)             ,&
           qsfc  (1:nCols)             ,tsfc  (1:nCols)           ,MskAnt(1:nCols)             ,&
           iMask (1:nCols)             ,zenith (1:nCols)          ,dlspr_sib (1:nCols)         ,&
           dcupr_sib(1:nCols)          ,LwSfcDown(1:nCols)        ,radvbc_sib(1:nCols)         ,&
           radvdc_sib(1:nCols)         ,radnbc_sib(1:nCols)       ,radndc_sib (1:nCols)        ,&
           HML   (1:nCols)             ,HUML (1:nCols)            ,HVML (1:nCols)              ,&
           TSK   (1:nCols)             ,GSW(1:nCols)              ,GLW(1:nCols)                ,&
           cldtot(1:nCols,1:kMax)      ,ySwSfcNet(1:nCols)        ,month(1:nCols)             ,&
           LwSfcNet(1:nCols)           ,pblh  (1:nCols)           ,QCF (1:nCols,1:kMax)        ,&
           QCL  (1:nCols,1:kMax)       ,mlsi  (1:nCols)           ,latco                       ,&
           Mmlen(1:nCols)              ,colrad(1:nCols)           ,idatec,dump(1:nCols,1:kMax ))

    DO i=1,nCols
       IF(mskant(i) == 1_i8 .and. tsea(i) <= 0.0e0_r8 .AND. tsurf(i) < tice+0.01e0_r8 ) THEN
              IF(intg.EQ.2) THEN
                 IF(istrt.EQ.0) THEN
                    tseam(i)=filta*tsea (i) + epsflt*(tseam(i)+xsea(i))
                    qsfc (i)=MAX(1.0e-12_r8,qsfc(i))
                    TSfcm(i,latco)=filta*TSfc0 (i,latco) + epsflt*(TSfcm(i,latco)+tsfc(i))
                    QSfcm(i,latco)=filta*QSfc0 (i,latco) + epsflt*(QSfcm(i,latco)+qsfc(i))
                 END IF
                 tsea (i) = xsea(i)
                 qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
                 TSfc0(i,latco) = tsfc  (i)
                 QSfc0(i,latco) = qsfc  (i)
              ELSE
                 tsea (i) = xsea(i)
                 tseam(i) = xsea(i)
                 qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
                 TSfc0(i,latco) = tsfc(i)
                 QSfc0(i,latco) = qsfc(i)
                 TSfcm(i,latco) = tsfc(i)
                 QSfcm(i,latco) = qsfc(i)
              END IF
       END IF
       IF(mskant(i) == 1_i8 .and. tsea(i).LT.0.0e0_r8.AND.tsurf(i).GE.tice+0.01e0_r8) THEN
              tseam(i)       = tsea (i)
              TSfcm(i,latco) = TSfc0(i,latco)
              QSfcm(i,latco) = QSfc0(i,latco)
       END IF
    END DO

    DO i=1,nCols
       IF(intg.EQ.2) THEN
          IF(istrt.EQ.0) THEN
                 ndvim (i      )=filta*ndvi (i) + epsflt*(ndvim(i)+xndvi(i))
          END IF
          ndvi (i)=ndvim (i      )
       ELSE
          ndvi (i)=ndvi (i      )
          ndvim(i)=ndvi (i      )
       END IF
    END DO

    ncount=0
    DO i=1,nCols
       sens        (i)          =  sensf    (i,latco)
       evap        (i)          =  latef    (i,latco)
       IF(iMask(i) >=1_i8)THEN 
          ncount=ncount+1
          bstar(i)        = bstar1(ncount)
          zlwup_SiB2(ncount,latco)        =zlwup(ncount)
          !TGS(i) = TSNOW*AREAS(i) + TG(i)*(1.0_r8-AREAS(i))
          sensf    (i,latco) = (hc   (ncount) + hg(ncount)*(1.0_r8-AREAS(ncount)) + hs(ncount)*AREAS(ncount))*(1.0_r8/dtc3x)!fss           (ncount)
          latef    (i,latco) = (ec   (ncount) + eg(ncount)*(1.0_r8-AREAS(ncount)) + es(ncount)*AREAS(ncount))*(1.0_r8/dtc3x)!fws           (ncount)*hltm
          sens     (i)       = (hc   (ncount) + hg(ncount)*(1.0_r8-AREAS(ncount)) + hs(ncount)*AREAS(ncount))*(1.0_r8/dtc3x) !fss     (ncount)
          evap     (i)       = (ec   (ncount) + eg(ncount)*(1.0_r8-AREAS(ncount)) + es(ncount)*AREAS(ncount))*(1.0_r8/dtc3x) !fws     (ncount)*hltm
          TSfc0    (i,latco) = tm0     (ncount,latco)
          QSfc0    (i,latco) = MAX(1.0e-12_r8,qm0 (ncount,latco))
          TSfcm    (i,latco) = tmm     (ncount,latco)
          QSfcm    (i,latco) = MAX(1.0e-12_r8,qmm (ncount,latco))
       END IF
    END DO    

  END SUBROUTINE SiB2_Driver


  !
  !------------------------------------------------------------
  !
  SUBROUTINE vert_interp(nsib      , & ! IN
       nzg       , & ! IN
       tzdep     , & ! IN
       glsm_slz  , & ! IN
       gl_sm     , & ! IN
       glsm_w_sib  ) ! OUT

    INTEGER, INTENT(IN   ) :: nsib
    INTEGER, INTENT(IN   ) :: nzg
    REAL(KIND=r8)   , INTENT(IN   ) :: tzdep     (nsib)
    REAL(KIND=r8)   , INTENT(IN   ) :: glsm_slz  (: )
    REAL(KIND=r8)   , INTENT(IN   ) :: gl_sm     (: )
    REAL(KIND=r8)   , INTENT(OUT  ) :: glsm_w_sib(nsib)

    REAL(KIND=r8)    :: zm        (nsib)
    REAL(KIND=r8)    :: wf        (nsib)
    REAL(KIND=r8)    :: zc        (nzg )
    REAL(KIND=r8)    :: wi        (nzg )
    REAL(KIND=r8)    :: dzlft
    INTEGER :: ZDM
    INTEGER :: k
    INTEGER :: kstart
    INTEGER :: L

    DO k=1,nzg
       zc(k)=glsm_slz(k)
    END DO

    DO k=1,nsib
       zm(k)=tzdep(k)
    END DO

    zdm=nsib
    KSTART=3
    !
    ! Transfere valores da grade de MAIOR resolucao (WI)
    !                     para a grade de MENOR resolucao (WF)
    !
    ! OS valores de WI devem estar definidos nos pontos de grade ZCS=zc/2
    ! OS valores de WF saem nos niveis ZMS = ZM/2
    !
    !
    !
    !    Dados da grade de maior resolucao
    !
    DO K=1,NZG
       WI(K) = gl_sm(k)
       !print*,'wi=',k,wi(k)
    END DO
    !
    !     Dado interpolado
    !
    wf(:)=0.0_r8
    !
    !     Valor de superficie:
    !
    WF(1)=WI(2)
    WF(2)=WI(2)
    !
    !
    DZLFT=0.0_r8
    L=2
    DO K=KSTART,ZDM
       !
       !    if(k==4) print*,'0',l,WF(K),WI(L),DZLFT
       !
       IF(DZLFT.NE.0.0_r8) THEN

          WF(K)=WF(K)+WI(L)*DZLFT
          !    if(k==4) print*,'1',l,WF(K),WI(L),DZLFT
          L=L+1

       END IF

70     CONTINUE

       IF(ZC(L).LE.ZM(K)) THEN

          WF(K)=WF(K)+WI(L)*(ZC(L)-ZC(L-1))

          !   if(k==4) print*,'2',l,WF(K),WI(L),ZC(L),zm(k)

          L=L+1
          DZLFT=0.0_r8
          IF (L>nzg) GO TO 1000
          GO TO 70
       ELSE

          WF(K)=WF(K)+WI(L)*(ZM(K)-ZC(L-1))
          DZLFT=ZC(L)-ZM(K)
       ENDIF
    ENDDO

1000 CONTINUE

    DO K=KSTART,ZDM
       !
       !   WF(K) =WF(K)/(ZM(K)-ZM(K-1))
       !         if(k==4)print*,zm(k),zc(nzg),ZM(K-1),WF(K)
       !
       IF (ZM(K) > ZC(nzg)) THEN
          WF(K) = WF(K)/(ZC(NZG)-ZM(K-1))
       ELSE
          WF(K)  = WF(K)/(ZM(K)-ZM(K-1))
       END IF
    END DO
    !
    !valores na grade do SIB
    !
    DO k=1,nsib
       glsm_w_sib(k)=WF(k)
       !print*,'SIB',k,glsm_w_sib(k)
    END DO
    !
    !check conservacao
    !srf - verifique se a integral de ambos calculos percorrem
    !srf - o mesmo intervalo
    !      print*,'        '
    !      sumf=0.0_r8
    !      DO K=2,ZDM
    !       sumf=sumf+wf(k)*(ZM(K)-ZM(K-1))
    !       print*,sumf,wf(k),zm(k),ZM(K)-ZM(K-1)
    !      ENDDO
    !      print*,'--------sumf-----',sumf
    !      sumi=0.0_r8
    !      DO K=2,nzg
    !       sumi=sumi+wi(k)*(glsm_slz(K)-glsm_slz(K-1))
    !       print*,k,sumi,wi(k),glsm_slz(K),(glsm_slz(K)-glsm_slz(K-1))
    !      ENDDO
    !      print*,'--------sumi-----',sumi, 100*(sumf-sumi)/sumi
    !
    RETURN
  END SUBROUTINE vert_interp
 
   SUBROUTINE extrak( w, dw, tbee, tphsat, rsoilm, cover, tph1, tph2, &
       psit, factor )
    REAL(KIND=r8), INTENT(in   ) :: w
    REAL(KIND=r8), INTENT(in   ) :: dw
    REAL(KIND=r8), INTENT(in   ) :: tbee
    REAL(KIND=r8), INTENT(in   ) :: tphsat
    REAL(KIND=r8), INTENT(in   ) :: rsoilm
    REAL(KIND=r8), INTENT(in   ) :: cover
    REAL(KIND=r8), INTENT(in   ) :: tph1
    REAL(KIND=r8), INTENT(in   ) :: tph2
    REAL(KIND=r8), INTENT(inout  ) :: psit
    REAL(KIND=r8), INTENT(inout  ) :: factor
    REAL(KIND=r8) :: rsoil
    REAL(KIND=r8) :: argg
    REAL(KIND=r8) :: hr
    REAL(KIND=r8) :: rplant
    !                --     -- (-b)
    !               |      dw |                  0
    ! psit = PHYs * | w - --- |      where w = -----
    !               |      2  |                  0s
    !                --     --
    psit   = tphsat * ( w-dw/2.0e0_r8 ) ** (-tbee)
    !
    !                      --                        --
    !                     |       --     -- (0.0027)   |
    !                     |      |      dw |           |
    !rsoil   = 101840.0 * |1.0 - | w - --- |           |
    !                     |      |      2  |           |
    !                     |       --     --            |
    !                      --                        --
    !
    rsoil  = 101840.0_r8 * (1.0_r8-( w-dw/2.0_r8) ** 0.0027_r8)
    !
    !                9.81       1
    !argg = psit * -------- * -------
    !               461.50     310.0
    !
    argg   = MAX ( -10.0e0_r8 , ((psit * 9.81e0_r8 / 461.5e0_r8) / 310.e0_r8))
    !
    !            --                       --
    !           |         9.81       1      |
    !hr   = EXP |psit * -------- * -------  |
    !           |        461.50     310.0   |
    !            --                       --
    !
    hr     = EXP ( argg )
    !
    !         rsoilm
    ! rsoil =--------- * hr
    !         rsoil
    !
    rsoil  = rsoilm /rsoil * hr
    !
    !          ( psit - tph2 - 50.0)
    !rplant = -------------------------
    !             ( tph1 - tph2 )
    !
    rplant = ( psit - tph2 -50.0_r8) / ( tph1 - tph2 )
    rplant = MAX ( 0.0e0_r8, MIN ( 1.0e0_r8, rplant ) )
    !                                                                     --                   --
    !                  --                 --                             |     --     -- (0.0027)|
    !                 |( psit - tph2 - 50)  |                            |    |      dw |        |
    !factor = cover * |---------------------| + (1 - cover) * 101840.0 * |1 - | w - --- |        |
    !                 |   ( tph1 - tph2 )   |                            |    |      2  |        |
    !                  --                 --                             |     --     --         |
    !                                                                     --                   --
    factor = cover * rplant + ( 1.0e0_r8 - cover ) * rsoil
    factor = MAX ( 1.e-6_r8, factor )
  END SUBROUTINE extrak 
  SUBROUTINE sextrp &
       (td    ,qa,ta,tg    ,tc    ,w     ,capac ,td0   ,qa0,ta0,tg0   ,tc0   ,w0    , &
       capac0,tdm   ,qam,tam,tgm   ,tcm   ,wm    ,capacm,istrt ,ncols ,nmax   , &
       epsflt,intg  ,tm0   ,qm0   ,tm    ,qm    ,tmm    ,qmm   , &
       nsoil  )
    INTEGER, INTENT(in   ) :: istrt
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: nmax
    INTEGER, INTENT(in   ) :: nsoil
    REAL(KIND=r8)   , INTENT(in   ) :: epsflt
    INTEGER, INTENT(in   ) :: intg
!    INTEGER, INTENT(in   ) :: latitu
    REAL(KIND=r8),    INTENT(in   ) :: tm    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: qm    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: td    (ncols,nsoil)
    REAL(KIND=r8),    INTENT(in   ) :: qa    (ncols)    
    REAL(KIND=r8),    INTENT(in   ) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: w     (ncols,3)
    REAL(KIND=r8),    INTENT(in   ) :: capac (ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: td0   (ncols,nsoil)
    REAL(KIND=r8),    INTENT(inout) :: qa0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ta0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg0   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tc0   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: w0    (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: capac0(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: tdm   (ncols,nsoil)
    REAL(KIND=r8),    INTENT(inout) :: qam    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tam    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tgm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tcm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wm    (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: capacm(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: tm0 (ncols)
    REAL(KIND=r8),    INTENT(inout) :: qm0 (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tmm (ncols)
    REAL(KIND=r8),    INTENT(inout) :: qmm (ncols)
    INTEGER :: i, nc, ii,k

    IF (intg == 2) THEN
       IF (istrt >= 1) THEN
          DO k=1, nsoil 
             DO i = 1, nmax
                td0   (i,k)  =td   (i,k)
             END DO
          END DO
          DO i = 1, nmax
             tm0   (i)  =tm   (i)
             qm0   (i)  =qm   (i)
             qa0   (i)  =qa   (i)
             ta0   (i)  =ta   (i)
             tg0   (i)  =tg   (i)
             tc0   (i)  =tc   (i)
             w0    (i,1)=w    (i,1)
             w0    (i,2)=w    (i,2)
             w0    (i,3)=w    (i,3)
             capac0(i,1)=capac(i,1)
             capac0(i,2)=capac(i,2)
             IF (capac0(i,2) > 0.0_r8 .AND. tg0(i) > 273.16_r8) THEN
                nc=0
                ii=0
                !DO ii = 1, ncols
                !   IF (imask(ii) >= 1) nc=nc+1
                !   IF (nc == i) EXIT
                !END DO
                !WRITE(UNIT=*,FMT=200)ii,latitu,i,capac0(i,2),tg0(i)
             END IF
          END DO
       ELSE
          DO k=1, nsoil 
             DO i = 1, nmax
                td0(i,k)=td0(i,k)+epsflt*(td(i,k)+tdm(i,k)-2.0_r8  *td0(i,k))
             END DO
          END DO

          DO i = 1, nmax
             qa0(i)=qa0(i)+epsflt*(qa(i)+qam(i)-2.0_r8  *qa0(i))
             ta0(i)=ta0(i)+epsflt*(ta(i)+tam(i)-2.0_r8  *ta0(i))
             tg0(i)=tg0(i)+epsflt*(tg(i)+tgm(i)-2.0_r8  *tg0(i))
             tc0(i)=tc0(i)+epsflt*(tc(i)+tcm(i)-2.0_r8  *tc0(i))

             tm0(i)=tm0(i)+epsflt*(tm(i)+tmm(i)-2.0_r8  *tm0(i))
             qm0(i)=qm0(i)+epsflt*(qm(i)+qmm(i)-2.0_r8  *qm0(i))

             IF(w0    (i,1) > 0.0_r8 ) THEN
                w0(i,1)=w0(i,1)+epsflt*(w(i,1)+wm(i,1)-2.0_r8  *w0(i,1))
             END IF
             IF(w0    (i,2) > 0.0_r8 ) THEN
                w0(i,2)=w0(i,2)+epsflt*(w(i,2)+wm(i,2)-2.0_r8  *w0(i,2))
             END IF
             IF(w0    (i,3) > 0.0_r8 ) THEN
                w0(i,3)=w0(i,3)+epsflt*(w(i,3)+wm(i,3)-2.0_r8  *w0(i,3))
             END IF
             IF(capac0(i,1) > 0.0_r8 ) THEN
                capac0(i,1)=capac0(i,1) &
                     +epsflt*(capac(i,1)+capacm(i,1)-2.0_r8*capac0(i,1))
             END IF
             IF(capac0(i,2) > 0.0_r8 ) THEN
                capac0(i,2)=capac0(i,2) &
                     +epsflt*(capac(i,2)+capacm(i,2)-2.0_r8*capac0(i,2))
             END IF
          END DO
          DO k=1, nsoil 
             DO i = 1, nmax
                tdm   (i,k)  =td0   (i,k)
             END DO
          END DO

          DO i = 1, nmax
             qam   (i)  =qa0   (i)
             tam   (i)  =ta0   (i)
             tgm   (i)  =tg0   (i)
             tcm   (i)  =tc0   (i)
             tmm   (i)  =tm0   (i)
             qmm   (i)  =qm0   (i)
             wm    (i,1)=w0    (i,1)
             wm    (i,2)=w0    (i,2)
             wm    (i,3)=w0    (i,3)
             capacm(i,1)=capac0(i,1)
             capacm(i,2)=capac0(i,2)
             IF (capacm(i,2) > 0.0_r8) tgm(i)=MIN(tgm(i),273.06_r8)
          END DO
          DO k=1, nsoil 
             DO i = 1, nmax
                td0   (i,k)  =td    (i,k)
             END DO
          END DO

          DO i = 1, nmax
             qa0   (i)  =qa    (i)
             ta0   (i)  =ta    (i)
             tg0   (i)  =tg    (i)
             tc0   (i)  =tc    (i)
             tm0   (i)  =tm    (i)
             qm0   (i)  =qm    (i)
             w0    (i,1)=w     (i,1)
             w0    (i,2)=w     (i,2)
             w0    (i,3)=w     (i,3)
             capac0(i,1)=capac (i,1)
             capac0(i,2)=capac (i,2)
             IF (capac0(i,2) > 0.0_r8 .AND. tg0(i) > 273.16_r8) THEN
                nc=0
                !DO ii = 1, ncols
                !   IF (imask(ii) >= 1) nc=nc+1
                !   IF (nc == i) EXIT
                !END DO
                !WRITE(UNIT=*,FMT=200)ii,latitu,i,capac0(i,2),tg0(i)
             END IF
          END DO
       END IF
    ELSE
       DO k=1, nsoil 
          DO i = 1, nmax
             tdm   (i,k)  =td   (i,k)
          END DO
       END DO

       DO i = 1, nmax
          qam   (i)  =qa   (i)
          tam   (i)  =ta   (i)
          tgm   (i)  =tg   (i)
          tcm   (i)  =tc   (i)
          tmm   (i)  =tm   (i)
          qmm   (i)  =qm   (i)
          wm    (i,1)=w    (i,1)
          wm    (i,2)=w    (i,2)
          wm    (i,3)=w    (i,3)
          capacm(i,1)=capac(i,1)
          capacm(i,2)=capac(i,2)
          IF (capacm(i,2) > 0.0_r8 .AND. tgm(i) > 273.16_r8) THEN
             nc=0
             !DO ii = 1, ncols
             !   IF (imask(ii) >= 1) nc=nc+1
             !   IF (nc == i) EXIT
             !END DO
             !WRITE(UNIT=*,FMT=650)ii,latitu,i,capacm(i,2),tgm(i)
          END IF
       END DO
       DO k=1, nsoil 
          DO i = 1, nmax
             td0   (i,k)  =td   (i,k)
          END DO
       END DO
       DO i = 1, nmax
          qa0   (i)  =qa   (i)
          ta0   (i)  =ta   (i)
          tg0   (i)  =tg   (i)
          tc0   (i)  =tc   (i)
          tm0   (i)  =tm   (i)
          qm0   (i)  =qm   (i)
          w0    (i,1)=w    (i,1)
          w0    (i,2)=w    (i,2)
          w0    (i,3)=w    (i,3)
          capac0(i,1)=capac(i,1)
          capac0(i,2)=capac(i,2)
       END DO
    END IF
!200 FORMAT(' CAPAC0 AND TG0 NOT CONSISTENT AT I,J,IS=',3I4, &
!         ' CAPAC=',G16.8,' TG=',G16.8)
!650 FORMAT(' CAPACM AND TGM NOT CONSISTENT AT I,J,IS=',3I4, &
!         ' CAPAC=',G16.8,' TG=',G16.8)
  END SUBROUTINE sextrp



  !-------------------------------------------------------------------

  SUBROUTINE SIB &
       ( nsib        ,len  , nsoil  , &                          ! array and loop bounds
       tau , &                                                  ! time independent variable
       ta, sha , &                                          ! CAS prognostic variables
       tc, tg, td, www, snow, capac, rst, tke  , &          ! prognostic variables
       thm,tm, qm, um,vm,sinclt,ps, bps, ros, ts, psb  , &                  ! atmospheric forcing
       cupr, lspr, radvbc, radnbc, radvdc, radndc , &          ! atmospheric forcing
       zb, spdm, dlwbot         , &                                  ! atmospheric forcing
       z0d, zlt, z1, z2, cc1, cc2, dd, zdepth, poros , &  ! surface parameters
       phsat, bee, respcp, vmax0, green, tran, ref , &          ! surface parameters
       gmudmu, trop, phc, trda, trdm, slti, shti  , &          ! surface parameters
       hltii, hhti, effcon, binter, gradm, atheta  , &          ! surface parameters
       btheta, aparc, wopt, zm, wsat, vcover, radfac , &  ! surface parameters
       thermk, satco, slope, chil, ztdep, sandfrac ,clayfrac, &          ! surface parameters
       pi, grav, dt, cp,cv,  hltm , &                          ! constants
       delta, asnow, kapa, snomel, clai, cww, pr0  , &          ! constants
       ribc, vkrmn, pco2m, po2m, stefan, grav2 , &          ! constants
       tice,gmt,gmq,gmu, snofac, &                                          ! constants
       fss, fws, drag, co2flx, cu, ct,ustar , &                  ! output
       hr , ect,eci ,egt,egi,egs,eg ,ec,es,hc,hg ,hs,chf,shf,roff,&  ! output
       ra ,rb ,rd ,rc ,rg ,ea ,etc,etg,radt,rsoil,& ! output
       ioffset , &        
       sibdrv, forcerestore, doSiBco2, fixday, dotkef  , &         ! logical switches
       ventmf, thvgm, louis,  respfactor , &        
       xgpp ,pco2ap,c4fract,d13cresp,d13cca,isotope,itype,delsig,sigki,bstar,zlwup,areas)

    IMPLICIT NONE

    !-------------------------------------------------------------------
    !     Subroutine  bldif conducts implicit time differencing for 
    !          turbulence kinetic energy,ventilation mass flux calculations,
    !          drag, and boundary layer processes, and calls subroutines 
    !          calculating PBL and surface fluxes. 

    !     REFERENCES: Sato, N., P. J. Sellers, D. A. Randall, E. K. Schneider, 
    !          J. Shukla, J. L Kinter III, Y-T, Hou, and Albertazzi (1989) 
    !          "Effects of implementing the simple biosphere model in a general
    !          circulation model. J. Atmos. Sci., 46, 2767-2782.

    !                 Sellers, P. J., D. A. Randall, C. J. Collatz, J. A. Berry,
    !          C. B. Field, D. A. Dazlich, C. Zhang, G. Collelo (1996) A revise 
    !          land-surface parameterization (SiB2) for atmospheric GCMs. Part 1:
    !          Model formulation. (accepted by JCL)


    !     MODIFICATIONS:
    !           - changed VQSAT call to VNQSAT.  kwitt 10/23
    !          - added in the prognostic stomatal conductance in addinc. changan 
    !          - moved sib diagnostics accumulation from dcontrol to bldif
    !                dd 950202

    !     SUBROUTINES CALLED:  VNQSAT, SNOW1, balan, VNTLAT
    !          DELHF, DELEF, NETRAD, SIBSLV, endtem, updat2, addinc
    !          inter2, balan, soilprop, soiltherm, begtem, rnload

    !     FUNCTIONS CALLED:
    !          none

    !     Argument list variables
    !       Intent: in
    INTEGER, INTENT(IN   ) :: nsib
    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: nsoil    ! array and loop bounds
    REAL(KIND=r8)   , INTENT(IN   ) :: tau      ! time independent variable (hr)
    REAL(KIND=r8)   , INTENT(INOUT) :: ta   (len)! CAS temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: sha  (len)! CAS water vapor mixing ratio (kg/kg)
    REAL(KIND=r8)   , INTENT(INOUT) :: tc   (len)    ! canopy (vegetation) temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: tg   (len)         ! surface boundary temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: td   (nsib,nsoil)  ! deep soil temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: www  (nsib,3)      ! soil wetness 
    REAL(KIND=r8)   , INTENT(INOUT) :: snow (nsib,2)      ! snow cover ( kg /m^2)
    REAL(KIND=r8)   , INTENT(INOUT) :: capac(nsib,2)  ! liquid interception store (kg/m^2)
    REAL(KIND=r8)   , INTENT(INOUT) :: rst  (len)     ! stomatal resistance
    !
    !  atmospheric forcing
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: tke(len)       ! turbulent kinetic energy
    REAL(KIND=r8)   , INTENT(IN   ) :: thm(len)       ! mixed layer potential temperature (K)
    REAL(KIND=r8)   , INTENT(IN OUT) :: tm(len)       ! mixed layer potential temperature (K)
    REAL(KIND=r8)   , INTENT(IN OUT) :: qm(len)

    REAL(KIND=r8)   , INTENT(IN   ) :: um (len) ! atmospheric forcing, INTENT(IN   ) :: thm(len)        ! mixed layer potential temperature (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: vm (len) ! atmospheric forcing, INTENT(IN   ) :: sh(len) ! mixed layer water vapor mixing ratio (kg/kg)
    REAL(KIND=r8)   , INTENT(IN   ) ::sinclt (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ps(len)        ! surface pressure (hPa)
    REAL(KIND=r8)   , INTENT(IN   ) :: bps(len)       ! (ps/1000)**kapa
    REAL(KIND=r8)   , INTENT(IN   ) :: ros (len)      ! surface air density (kg/m^3)
    REAL(KIND=r8)   , INTENT(IN   ) :: ts  (len)      ! surface mixed layer air temperature
    REAL(KIND=r8)   , INTENT(IN   ) :: psb (len)      ! boundary layer mass depth (hPa)
    REAL(KIND=r8)   , INTENT(IN   ) :: cupr(len)      ! convective precipitation rate (mm/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: lspr(len)      ! stratiform precipitation rate (mm/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: radvbc(len)    ! surface incident visible direct beam (W/m^2)
    REAL(KIND=r8)   , INTENT(IN   ) :: radnbc(len)    ! surface incident near IR direct beam (W/m^2)
    REAL(KIND=r8)   , INTENT(IN   ) :: radvdc(len)    ! surface incident visible diffuse beam (W/m^2)
    REAL(KIND=r8)   , INTENT(IN   ) :: radndc(len)    ! surface incident near IR diffuse beam (W/m^2)
    REAL(KIND=r8)   , INTENT(IN   ) :: zb(len)        ! boundary layer thickness (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: spdm(len)      ! boundary layer wind speed (m/s) 
    REAL(KIND=r8)   , INTENT(IN   ) :: dlwbot(len)    ! surface incident longwave radiation (W/m^2)                

    !  surface parameters

    REAL(KIND=r8)   , INTENT(IN   ) :: z0d (len)      ! surface roughness length (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: zlt (len)      ! leaf area index
    REAL(KIND=r8)   , INTENT(IN   ) :: z1  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: z2  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: cc1 (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: cc2 (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: dd  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zdepth(nsib,3) ! porosity * soil hydrology model layer depths (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros(len)     ! soil porosity
    REAL(KIND=r8)   , INTENT(IN   ) :: phsat(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: bee(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: respcp(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: vmax0(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: green(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tran(nsib,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: ref(nsib,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: gmudmu(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: trop(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: phc(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: trda(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: trdm(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: slti(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: shti(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hltii(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hhti(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: effcon(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: binter(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: gradm(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: atheta(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: btheta(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: aparc(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: wopt(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zm(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsat(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: vcover(len)   ! vegetation cover fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: radfac(nsib,2,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: thermk(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: satco(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: slope(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: chil(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ztdep(nsib,nsoil) ! soil thermal model layer depths (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: sandfrac(len)     ! soil texture sand fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: clayfrac (len) ! soil texture sand fraction

    REAL(KIND=r8)   , INTENT(IN   ) :: pi          ! 3.1415926....
    REAL(KIND=r8)   , INTENT(IN   ) :: grav          ! gravitational acceleration (m/s^2)
    REAL(KIND=r8)   , INTENT(IN   ) :: dt          ! time step (s)
    REAL(KIND=r8)   , INTENT(IN   ) :: cp          ! heat capacity of dry air (J/(kg K) )
    REAL(KIND=r8)   , INTENT(IN   ) :: cv          ! heat capacity of water vapor (J/(kg K))
    REAL(KIND=r8)   , INTENT(IN   ) :: hltm          ! latent heat of vaporization (J/kg)
    REAL(KIND=r8)   , INTENT(IN   ) :: delta
    REAL(KIND=r8)   , INTENT(IN   ) :: asnow
    REAL(KIND=r8)   , INTENT(IN   ) :: kapa
    REAL(KIND=r8)   , INTENT(IN   ) :: snomel
    REAL(KIND=r8)   , INTENT(IN   ) :: clai          ! leaf heat capacity
    REAL(KIND=r8)   , INTENT(IN   ) :: cww          ! water heat capacity
    REAL(KIND=r8)   , INTENT(IN   ) :: pr0
    REAL(KIND=r8)   , INTENT(IN   ) :: ribc
    REAL(KIND=r8)   , INTENT(IN   ) :: vkrmn          ! von karmann's constant
    REAL(KIND=r8)   , INTENT(IN   ) :: pco2m(len)! CO2 partial pressure (Pa)
    REAL(KIND=r8)   , INTENT(IN   ) :: po2m          ! O2 partial pressure (Pa)
    REAL(KIND=r8)   , INTENT(IN   ) :: stefan    ! stefan-boltzman constant
    REAL(KIND=r8)   , INTENT(IN   ) :: grav2          ! grav / 100 Pa/hPa
    REAL(KIND=r8)   , INTENT(IN   ) :: tice          ! freezing temperature (KP
    REAL(KIND=r8)   , INTENT(INOUT) :: gmt(len,3)
    REAL(KIND=r8)   , INTENT(INOUT) :: gmq(len,3)
    REAL(KIND=r8)   , INTENT(INOUT) :: gmu(len,4)
    REAL(KIND=r8)   , INTENT(INOUT) :: snofac
    REAL(KIND=r8)   , INTENT(IN   ) :: delsig(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: sigki(len)
    !     Intent:out
    !fss, fws, drag, co2flx, cu, ct,ustar , &        
    REAL(KIND=r8)    ,INTENT(IN OUT) :: fss(len)    ! surface sensible heat flux (W/m^2)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: fws(len)    ! surface evaporation (kg/m^2/s)
    REAL(KIND=r8)    ,INTENT(OUT   ) :: co2flx(len) ! surface CO2 flux
    REAL(KIND=r8)    ,INTENT(OUT   ) :: ct(len)
    REAL(KIND=r8)    ,INTENT(OUT   ) :: drag(nsib,2)   ! surface drag coefficient

    REAL(KIND=r8)    ,INTENT(OUT   ) :: cu    (len)
    REAL(KIND=r8)    ,INTENT(OUT   ) :: ustar (len)  ! friction velocity (m/s)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: hr    (len)!,INTENT(INOUT        )
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ect   (len)     ! transpiration (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: eci   (len)     ! canopy interception evaporation (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: egt   (len)   
    REAL(KIND=r8)    ,INTENT(IN OUT) :: egi   (len)     ! ground interception evaporation (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: egs   (len)     ! soil evaporation (J)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: eg    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ec    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: es(len)

    REAL(KIND=r8)    ,INTENT(IN OUT) :: hc    (len)     
    REAL(KIND=r8)    ,INTENT(IN OUT) :: hg    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: hs    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: chf   (len)     ! canopy heat flux (W/m^2)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: shf   (len)     ! soil heat flux (W/m^2)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: roff  (len)   ! total runoff (surface and subsurface)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ra    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rb    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rd    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rc    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rg    (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: ea    (len)      ! canopy airspace water vapor pressure (hPa)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: etc   (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: etg   (len)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: radt  (len,3)
    REAL(KIND=r8)    ,INTENT(IN OUT) :: rsoil (len)
    REAL(KIND=r8)    ,INTENT(OUT   ) ::  bstar(len)
    REAL(KIND=r8)    ,INTENT(OUT   ) ::  zlwup(len)
    INTEGER, INTENT(IN    ) :: ioffset   ! subdomain offset
    LOGICAL, INTENT(IN    ) :: sibdrv           ! run SiB with prescribe meteorology
    LOGICAL, INTENT(IN    ) :: forcerestore          ! do force-restore soil thermodynamics
    LOGICAL, INTENT(IN    ) :: doSiBco2          ! prognostic co2
    LOGICAL, INTENT(IN    ) :: fixday           ! perpetual day of year conditions
    LOGICAL, INTENT(IN    ) :: dotkef           ! use changan instead of deardorff flux 
    REAL(KIND=r8)   , INTENT(IN OUT) :: ventmf(len) ! ventilation mass flux
    REAL(KIND=r8)   , INTENT(IN OUT) :: thvgm(len)
    INTEGER, INTENT(IN    ) :: louis            
    REAL(KIND=r8)   , INTENT(IN    ) :: respfactor(nsib,nsoil+1)
    REAL(KIND=r8)   , INTENT(OUT   ) :: xgpp(len)   ! gross primary productivity (micromoles/m**2/s)
    REAL(KIND=r8)   , INTENT(IN OUT) :: pco2ap(len)    ! canopy air space pCO2 (Pa)
    REAL(KIND=r8)   , INTENT(IN OUT) :: c4fract(len)   ! fraction of C4 vegetation  
    REAL(KIND=r8)   , INTENT(IN OUT) :: d13cresp(len)  ! del 13C of respiration (per mil vs PDB)
    REAL(KIND=r8)   , INTENT(IN OUT) :: d13cca(len)    ! del 13C of canopy CO2 (per mil vs PDB)
    LOGICAL         , INTENT(IN    ) :: isotope   
    INTEGER         , INTENT(IN   ) :: itype(len)
    REAL(KIND=r8)   , INTENT(IN OUT)  :: areas(len)

    ! will be moles air / m^2 in phosib
    !     !REAL(KIND=r8)     :: qp2(nsib,nnqp2) ! time average diagnostic fields (all points)
    !     REAL(KIND=r8)     :: pbp1(npbp1+1,ijtlen) ! time series diagnostic fields (select points)
    !     REAL(KIND=r8)     :: qp3(nsib, nsoil, nnqp3)
    !     REAL(KIND=r8)     :: pbp2(nsoil, npbp2, ijtlen)

    !      integer :: nnqp2mx  ! diagnostics array bounds
    !integer :: imultpbp(ijtlen)    ! grid points to save pbp at
    !integer :: indxqp2(nnqp2mx)! index of saved diagnostics
    !integer :: indxpbp1(npbp1mx)  ! index of saved diagnostics
    !integer :: nnqp3mx
    !      integer :: nnqp3
    !integer :: indxqp3(nnqp3mx)         
    !      integer :: npbp2mx
    !      integer :: npbp2
    !      integer :: indxpbp2(npbp2mx)
    !logical :: qpintp
    !logical :: histpp
    !logical :: doqp2(nnqp2mx)
    !logical :: doqp3(nnqp3mx)
    !  parameterization


    !     Local variables: 
    REAL(KIND=r8)     :: drag2(nsib)   ! surface drag coefficient
    REAL(KIND=r8)     :: cflux(len)  ! new formulation of CO2 flux (phosib)
    REAL(KIND=r8)     :: cas_cap_heat(len)  ! CAS heat capacity (J/K m^2)
    REAL(KIND=r8)     :: cas_cap_vap(len)   ! CAS 
    REAL(KIND=r8)     :: cas_cap_co2(len)   ! CAS CO2 capacity (m / m^2)
    REAL(KIND=r8)     :: dtinv
    !     REAL(KIND=r8)     :: AUXeadem
    REAL(KIND=r8)     :: psy(len)     ! psychrometric 'constant'
    REAL(KIND=r8)     :: tha(len)     ! canopy airspace potential temperature (K)
    REAL(KIND=r8)     :: em(len)      ! mixed layer water vapor pressure (hPa)
    REAL(KIND=r8)     :: d(len)       ! dd corrected for snow covered canopy
    REAL(KIND=r8)     :: rbc(len)     ! cc1 corrected for snow covered canopy
    REAL(KIND=r8)     :: rdc(len)     ! cc2 corrected for snow covered canopy
    REAL(KIND=r8)     :: etmass(len)  ! evapotranspiration
    REAL(KIND=r8)     :: totwb(len)   ! total surface and soil water at beginning of timestep
!    REAL(KIND=r8)     :: chf(len)     ! canopy heat flux (W/m^2)
    !     REAL(KIND=r8)     :: ahf(len)     ! CAS heat flux (W/m^2)
!    REAL(KIND=r8)     :: shf(len)     ! soil heat flux (W/m^2)
!    REAL(KIND=r8)     :: egs(len)     ! soil evaporation (J)
!    REAL(KIND=r8)     :: egi(len)     ! ground interception evaporation (J)
!    REAL(KIND=r8)     :: hc(len)     
!    REAL(KIND=r8)     :: hg(len)
    !     REAL(KIND=r8)     :: esdif(len)
    REAL(KIND=r8)     :: heaten(len)
    REAL(KIND=r8)     :: hflux(len)
    REAL(KIND=r8)     :: tgs(len)     ! bare ground and snow surface mean temperature
    REAL(KIND=r8)     :: tsnow(len)
    REAL(KIND=r8)     :: czc(len)     ! canopy heat capacity
    REAL(KIND=r8)     :: btc(len)
    REAL(KIND=r8)     :: btg(len)
    !     REAL(KIND=r8)     :: btci(len)
    !     REAL(KIND=r8)     :: btct(len)
    !     REAL(KIND=r8)     :: btgs(len)
    REAL(KIND=r8)     :: bts(len)
    REAL(KIND=r8)     :: rstfac(len,4)
    REAL(KIND=r8)     :: wc(len)
    REAL(KIND=r8)     :: wg(len)
    REAL(KIND=r8)     :: satcap(len,2)
    REAL(KIND=r8)     :: csoil(len)
    REAL(KIND=r8)     :: gect(len)
    REAL(KIND=r8)     :: geci(len)
    REAL(KIND=r8)     :: gegs(len)
    REAL(KIND=r8)     :: gegi(len)
    REAL(KIND=r8)     :: zmelt(len)
    REAL(KIND=r8)     :: rds(len)
    REAL(KIND=r8)     :: rib(len)     ! soil resistance
    REAL(KIND=r8)     :: ecmass(len)
    REAL(KIND=r8)     :: hrr(len)
    REAL(KIND=r8)     :: bintc(len) 
    REAL(KIND=r8)     :: aparkk(len)
    REAL(KIND=r8)     :: wsfws(len)
    REAL(KIND=r8)     :: wsfht(len)
    REAL(KIND=r8)     :: wsflt(len)
    REAL(KIND=r8)     :: wci(len)
    REAL(KIND=r8)     :: whs(len)
    REAL(KIND=r8)     :: omepot(len)
    REAL(KIND=r8)     :: assimpot(len)
    REAL(KIND=r8)     :: assimci(len)
    REAL(KIND=r8)     :: antemp(len)
    REAL(KIND=r8)     :: assimnp(len)
    REAL(KIND=r8)     :: wags(len)
    REAL(KIND=r8)     :: wegs(len)
    REAL(KIND=r8)     :: pfd(len)
    REAL(KIND=r8)     :: assim(len)
    REAL(KIND=r8)     :: zmstscale(len,2)
    REAL(KIND=r8)     :: zltrscale(len)
    REAL(KIND=r8)     :: zmlscale(len)
    REAL(KIND=r8)     :: drst(len)    ! stomatal resistance increment 
    REAL(KIND=r8)     :: soilq10(len,nsoil+1)
    REAL(KIND=r8)     :: ansqr(len)
    REAL(KIND=r8)     :: soilscaleold(len)
    REAL(KIND=r8)     :: anwtc(len)
    REAL(KIND=r8)     :: hgdtg(len)
    !     REAL(KIND=r8)     :: hgdth(len)
    REAL(KIND=r8)     :: hgdta(len)
    REAL(KIND=r8)     :: hsdts(len)
    !     REAL(KIND=r8)     :: hsdtc(len)
    !     REAL(KIND=r8)     :: hsdth(len)
    REAL(KIND=r8)     :: hsdta(len)
    !     REAL(KIND=r8)     :: hcdtg(len)
    REAL(KIND=r8)     :: hcdtc(len)
    !     REAL(KIND=r8)     :: hcdth(len)
    REAL(KIND=r8)     :: hcdta(len)
    !     REAL(KIND=r8)     :: hadtg(len)
    !     REAL(KIND=r8)     :: hadtc(len)
    REAL(KIND=r8)     :: hadth(len)
    REAL(KIND=r8)     :: hadta(len)
    !     REAL(KIND=r8)     :: aag(len)
    !     REAL(KIND=r8)     :: aac(len)
    !     REAL(KIND=r8)     :: aam(len)
    REAL(KIND=r8)     :: fc(len)
    REAL(KIND=r8)     :: fg(len)
    !Bio these are the derivatives of the LW fluxes needed for 
    !Bio the prognostic canopy
    REAL(KIND=r8)     :: lcdtc(len)
    REAL(KIND=r8)     :: lcdtg(len)
    REAL(KIND=r8)     :: lcdts(len)
    REAL(KIND=r8)     :: lgdtg(len)
    REAL(KIND=r8)     :: lgdtc(len)
    REAL(KIND=r8)     :: lsdts(len)
    REAL(KIND=r8)     :: lsdtc(len)
    REAL(KIND=r8)     :: egdtg(len)
    !     REAL(KIND=r8)     :: egdqm(len)
    !     REAL(KIND=r8)     :: ecdtg(len)
    REAL(KIND=r8)     :: ecdtc(len)
    !     REAL(KIND=r8)     :: ecdqm(len)
    !     REAL(KIND=r8)     :: egdtc(len)
    !     REAL(KIND=r8)     :: deadtg(len)
    !     REAL(KIND=r8)     :: deadtc(len)
    !     REAL(KIND=r8)     :: deadqm(len)
    REAL(KIND=r8)     :: ecdea(len)
    REAL(KIND=r8)     :: egdea(len)
    REAL(KIND=r8)     :: esdts(len)
    REAL(KIND=r8)     :: esdea(len)
    REAL(KIND=r8)     :: eadea(len)
    REAL(KIND=r8)     :: eadem(len)
    !     REAL(KIND=r8)     :: bbg(len)
    !     REAL(KIND=r8)     :: bbc(len)
    !     REAL(KIND=r8)     :: bbm(len)
    REAL(KIND=r8)     :: dtg(len,2)   ! surface ground and snow temperature increments (K)
    REAL(KIND=r8)     :: dtc(len)     ! canopy temperature increment (K)
    REAL(KIND=r8)     :: dta(len)     ! CAS temperature increment (K)
    REAL(KIND=r8)     :: dea(len)     ! CAS moisture increment (Pa)
    REAL(KIND=r8)     :: dtd(len,nsoil)   ! deep soil temperature increments (K)
    REAL(KIND=r8)     :: q3l(len)
    REAL(KIND=r8)     :: q3o(len)
    REAL(KIND=r8)     :: exo(len)
    REAL(KIND=r8)     :: qqq(len,3)
    REAL(KIND=r8)     :: zmelt1(len)
    REAL(KIND=r8)     :: evt(len)
    REAL(KIND=r8)     :: eastar(len)
    REAL(KIND=r8)     :: roffo(len)
    REAL(KIND=r8)     :: zmelt2(len)
    REAL(KIND=r8)     :: rha(len)     ! canopy airspace relative humidity
    REAL(KIND=r8)     :: ggl(len)
    REAL(KIND=r8)     :: egmass(len)
    !     REAL(KIND=r8)     :: etci(len)
    !     REAL(KIND=r8)     :: etgi(len)
    !     REAL(KIND=r8)     :: etct(len)
    !     REAL(KIND=r8)     :: etgs(len)
    REAL(KIND=r8)     :: ets(len)
    REAL(KIND=r8)     :: slamda(len,nsoil)  ! soil thermal conductivities
    REAL(KIND=r8)     :: shcap(len,nsoil)   ! soil heat capacities
    REAL(KIND=r8)     :: fac1(len)
    REAL(KIND=r8)     :: dth(nsib)   ! mixed layer potential temperature increment (K)
    REAL(KIND=r8)     :: dqm(nsib)   ! mixed layer water vapor mixing ratio increment (kg/kg)
    REAL(KIND=r8)     :: assimn(len)
    REAL(KIND=r8)     :: soilscale(len,nsoil+1)
    REAL(KIND=r8)     :: czh(len)    ! surface layer heat capacity
    REAL(KIND=r8)     :: radc3(len,2)
    REAL(KIND=r8)     :: radn(len,2,2)
    REAL(KIND=r8)     :: ustaro(len)  ! friction velocity (m/s) ( for oceanic z0)
    REAL(KIND=r8)     :: cuo(len)   ! ( for oceanic z0)
    REAL(KIND=r8)     :: z0(len)    ! surface roughness length corrected for canopy snow (m)
    REAL(KIND=r8)     :: wwwtem(len,3)  ! soil wetness copy
    REAL(KIND=r8)     :: thgeff(len), tgeff(len)
    REAL(KIND=r8)     :: shgeff(len), canex(len)
    REAL(KIND=r8)     :: cuprt(len), lsprt(len) ! copies of cupr and lspr
    REAL(KIND=r8)     :: thmtem(len), shtem(len)  ! copies of thm and sh
    REAL(KIND=r8)     :: zzwind(len), zztemp(len)
    REAL(KIND=r8)     :: respg(len)
    REAL(KIND=r8)     :: discrim(len)
    REAL(KIND=r8)     :: discrim2(len)
    !     REAL(KIND=r8)     :: discrim3(len)
    REAL(KIND=r8)     :: closs(len)     ! vegetation IR loss
    REAL(KIND=r8)     :: gloss(len)     ! ground IR loss
    REAL(KIND=r8)     :: sloss(len)     ! snow IR loss
    REAL(KIND=r8)     :: dtc4(len)      ! 1st derivative of vegetation T^4
    REAL(KIND=r8)     :: dtg4(len)      ! 1st derivative of ground T^4
    REAL(KIND=r8)     :: dts4(len)      ! 1st derivative of snow T^4
    !itb    ! discrim factors-neil suits' programs
    REAL(KIND=r8)     :: pco2i(len)     ! leaf internal pCO2 (Pa)
    REAL(KIND=r8)     :: pco2ap_old(len)  !previous time step pCO2 for cfrax eqns
    REAL(KIND=r8)     :: pco2c(len)     ! chloroplast pCO2 (Pa)
    REAL(KIND=r8)     :: pco2s(len)     ! leaf surface pCO2 (Pa)
    REAL(KIND=r8)     :: d13cm(len)     ! del 13C of mixed layer (constant, per mil vs PDB)
    REAL(KIND=r8)     :: co2cap(len)    ! moles of air in the canopy (moles/canopy air space)
    REAL(KIND=r8)     :: kiepsc3(len)   ! kinetic isotope effect during C3 photosynthesis
    REAL(KIND=r8)     :: kiepsc4(len)   ! kinetic isotope effect during C4 photosynthesis
    REAL(KIND=r8)     :: d13cassc3(len)  ! del 13C of CO2 assimilated by C3 plants
    REAL(KIND=r8)     :: d13cassc4(len)  ! del 13C of CO2 assimilated by C4 plants
    REAL(KIND=r8)     :: d13casstot(len) ! del 13C of CO2 assimilated by all plants
    REAL(KIND=r8)     :: flux13c(len)   ! turbulent flux of 13CO2 out of the canopy
    REAL(KIND=r8)     :: flux12c(len)   ! turbulent flux of 12CO2 out of the canopy
    REAL(KIND=r8)     :: d13cflux(len)  ! del 13C of the turb flux of CO2 out of the canopy 
    REAL(KIND=r8)     :: conductra(len) ! 1/ra, for neil's program
    REAL(KIND=r8)     :: flux_turb(len) ! nsuits -- turb flux of CO2 out of canopy
    REAL(KIND=r8)     :: umom(len)
    REAL(KIND=r8)     :: vmom(len)
    REAL(KIND=r8)     :: gb100  
    REAL(KIND=r8)     :: fac
    REAL(KIND=r8)     :: gby100
    !REAL(KIND=r8)     :: dtmdt
    !REAL(KIND=r8)     :: dqmdt
    !REAL(KIND=r8)     :: dtm 
    !REAL(KIND=r8)     :: dqms
    REAL(KIND=r8)     :: aaa(len)
    !REAL(KIND=r8)     :: sensf(len)
    !REAL(KIND=r8)     :: latef(len)
    !itb   added in for neil's programs
    REAL(KIND=r8)     ::  gbycp
    REAL(KIND=r8)     ::  gbyhl,esat12
    REAL(KIND=r8)     :: rmi   (len)
    REAL(KIND=r8)     :: rhi   (len)
    REAL(KIND=r8)     :: am   (len)  
    REAL(KIND=r8)     :: cti (len)
    REAL(KIND=r8)     :: cui (len)
    INTEGER :: i
    REAL(KIND=r8)     :: ah  (len)  
    REAL(KIND=r8)     :: al    (len) 
    REAL(KIND=r8)    :: cuni(len)
    REAL(KIND=r8)    :: ctni(len)
    REAL(KIND=r8)    :: u2(len)
    REAL(KIND=r8)    :: qh(len)
    LOGICAL          :: jstneu

    !      integer :: i, j, k, n, l

    !-------------------------------------------------------------------


    dtinv = 1.0_r8 / dt   ! inverse time step

    !      print*,'SiB.F: tau=',tau

    !     some initialization, copy of soil wetness
    DO i = 1,len
       tsnow(i)      = MIN(tg(i),tice)
       cuprt(i)      = cupr(i) * 0.001_r8 !convert mm/s ---> m/s
       lsprt(i)      = lspr(i) * 0.001_r8 !convert mm/s ---> m/s
       zmelt(i)      = 0.0_r8
       roff (i)      = 0.0_r8
       wwwtem(i,1)   = www(i,1)
       wwwtem(i,2)   = www(i,2)
       wwwtem(i,3)   = www(i,3)
       pco2ap_old(i) = pco2ap(i)
       heaten(i)     = 0.0_r8
       rg(i)         = 0.0_r8
       egt           = 0.0_r8
    ENDDO

    !     first guesses for ta and ea (see temrec 120)

    DO I=1,len
       THA(I) = TA(I) 
    ENDDO

    DO I=1,len
       !EA(I) = SHA(I) * PS(I) / (0.622_r8 + SHA(I))
       EA(I)=fpvs2es5(TA(I) )/100.0_r8
       !EM(I) = qm (I) * PS(I) / (0.622_r8 + qm (I))
       EM(I) =fpvs2es5(TM(I) )/100.0_r8
    ENDDO

    !    distribute incident radiation between canopy and surface
    CALL rnload( &
         len    , &! INTENT(IN   ) :: nlen
         radvbc , &! INTENT(IN   ) :: radvbc(nlen)
         radvdc , &! INTENT(IN   ) :: radvdc(nlen)
         radnbc , &! INTENT(IN   ) :: radnbc(nlen)
         radndc , &! INTENT(IN   ) :: radndc(nlen)
         dlwbot , &! INTENT(IN   ) :: dlwbot(nlen)
         VCOVER , &! INTENT(IN   ) :: VCOVER(nlen)
         thermk , &! INTENT(IN   ) :: thermk(nlen)
         radfac , &! INTENT(IN   ) :: radfac(nlen,2,2,2)
         radn   , &! INTENT(OUT  ) :: radn  (nlen,2,2)
         radc3    )! INTENT(OUT  ) :: radc3 (nlen,2)

    DO i = 1,len
       CANEX (i)  = 1.0_r8-( SNOW(i,1)*5.0_r8-Z1(i))/(Z2(i)-Z1(i))
       CANEX (i)  = MAX( 0.1E0_r8, CANEX(i) )
       CANEX (i)  = MIN( 1.0E0_r8, CANEX(i) )
       D     (i)  = Z2(i) - ( Z2(i)-DD(i) ) * CANEX(i)
       Z0    (i)  = (Z0D(i)/( Z2(i)-DD(i) )) * ( Z2(i)-D(i) )
       RBC   (i)  = CC1(i)/CANEX(i)
       RDC   (i)  = CC2(i)*CANEX(i)
       AREAS (i)  = MIN(1.0E0_r8, ASNOW*SNOW(i,2))
       SATCAP(i,1)= ZLT(i)*0.0001_r8 * CANEX(i)
       !
       ! Collatz-Bounoua change satcap(2) to 0.0002
       !
       !bl      SATCAP(i,2) = 0.002_r8
       !SATCAP(i,2) = 0.0002_r8             ! lahouari
       SATCAP(i,2) =  ZLT(i)*0.0001_r8 * CANEX(i)            ! lahouari
    END DO

    !    initialize energy and water budgets
    CALL balan(  1     , & ! INTENT(IN   ) :: iplace
         tau   , & ! INTENT(IN   ) :: tau
         zdepth, & ! INTENT(IN   ) :: zdepth(nlen,3)
         wwwtem, & ! INTENT(IN   ) :: www   (nlen,3)
         capac , & ! INTENT(IN   ) :: capac (nlen,2)
         cupr  , & ! INTENT(IN   ) :: ppc   (nlen)
         lspr  , & ! INTENT(IN   ) :: ppl   (nlen)
         roff  , & ! INTENT(INOUT) :: roff  (nlen)
         etmass, & ! INTENT(INOUT) :: etmass(nlen)
         totwb , & ! INTENT(INOUT) :: totwb (nlen)
         radt  , & ! INTENT(IN   ) :: radt  (nlen,2)
         chf   , & ! INTENT(IN   ) :: chf   (nlen)
         shf   , & ! INTENT(IN   ) :: shf   (nlen)
         dt    , & ! INTENT(IN   ) :: dtt
         ect   , & ! INTENT(IN   ) :: ect   (nlen)
         eci   , & ! INTENT(IN   ) :: eci   (nlen)
         egs   , & ! INTENT(IN   ) :: egs   (nlen)
         egi   , & ! INTENT(IN   ) :: egi   (nlen)
         hc    , & ! INTENT(IN   ) :: hc    (nlen)
         hg    , & ! INTENT(IN   ) :: hg    (nlen)
         heaten, & ! INTENT(IN   ) :: heaten(nlen)
         hflux , & ! INTENT(IN   ) :: hflux (nlen)
         snow  , & ! INTENT(IN   ) :: snoww (nlen,2)    
         thm   , & ! INTENT(IN   ) :: thm   (nlen)
         tc    , & ! INTENT(IN   ) :: tc    (nlen)
         tg    , & ! INTENT(IN   ) :: tg    (nlen)
         tgs   , & ! INTENT(IN   ) :: tgs   (nlen)
         td    , & ! INTENT(IN   ) :: td    (nlen, nsoil)
         ps    , & ! INTENT(IN   ) :: ps    (nlen)  
         kapa  , & ! INTENT(IN   ) :: kapa
         !nsib , & ! INTENT(IN   ) :: nsib
         len   , & ! INTENT(IN   ) :: nlen
         ioffset,& ! INTENT(IN   ) :: ioffset
         nsoil )   ! INTENT(IN   ) :: nsoil

    CALL begtem(tc      , &         ! INTENT(IN   ) :: tc    (nlen)
         tg      , &         ! INTENT(IN   ) :: tg    (nlen)
         cp      , &         ! INTENT(IN   ) :: cpair    
         hltm    , &         ! INTENT(IN   ) :: hlat
         ps      , &         ! INTENT(IN   ) :: psur  (nlen)
         snomel  , &         ! INTENT(IN   ) :: snomel
         zlt     , &         ! INTENT(IN   ) :: zlt   (nlen)
         clai    , &         ! INTENT(IN   ) :: clai
         cww     , &         ! INTENT(IN   ) :: cw
         wwwtem  , &         ! INTENT(IN   ) :: www   (len,3)
         poros   , &         ! INTENT(IN   ) :: poros (nlen)
         pi      , &         ! INTENT(IN   ) :: pie
         psy     , &         ! INTENT(OUT  ) :: psy   (len) 
         phsat   , &         ! INTENT(IN   ) :: phsat (nlen)
         bee     , &         ! INTENT(IN   ) :: bee   (nlen)
         czc     , &         ! INTENT(OUT  ) :: ccx   (len)
         czh     , &         ! INTENT(OUT  ) :: cg    (nlen)
         phc     , &         ! INTENT(IN   ) :: phc   (nlen)
         tgs     , &         ! INTENT(OUT  ) :: tgs   (len)
         etc     , &         ! INTENT(OUT  ) :: etc   (len)
         etg     , &         ! INTENT(OUT  ) :: etgs  (len)
         btc     , &         ! INTENT(OUT  ) :: getc  (len)
         btg     , &         ! INTENT(OUT  ) :: getgs (len)
         rstfac  , &         ! INTENT(OUT  ) :: rstfac(len,4)
         rsoil   , &         ! INTENT(OUT  ) :: rsoil (len)
         hr      , &         ! INTENT(OUT  ) :: hr    (len)
         wc      , &         ! INTENT(OUT  ) :: wc    (len)
         wg      , &         ! INTENT(OUT  ) :: wg    (len)
         snow    , &         ! INTENT(IN   ) :: snoww (nlen,2)
         capac   , &         ! INTENT(IN   ) :: capac (nlen,2)
         areas   , &         ! INTENT(IN   ) :: areas (len)
         satcap  , &         ! INTENT(IN   ) :: satcap(len,2)   
         csoil   , &         ! INTENT(OUT  ) :: csoil (len)
         tice    , &         ! INTENT(IN   ) :: tf
         grav    , &         ! INTENT(IN   ) :: g
         snofac  , &         ! INTENT(OUT  ) :: snofac
         len     , &         ! INTENT(IN   ) :: len
         sandfrac, &         ! INTENT(IN   ) :: 
         clayfrac, &         ! INTENT(IN   ) :: 
         nsib    , &         ! INTENT(IN   ) :: nlen
         forcerestore )      ! INTENT(IN   ) :: forcerestore


    !pl now that we have the new psy, calculate the new CAS capacities
    !itb...PL made this max(4.,z2(:)), but we might have to boost the
    !itb...min value upwards...
    cas_cap_heat(:) = ros(:) * cp * MAX(4.0_r8,z2(:))
    !itb...I think cas_cap_vap should use cv instead of cp...
    !           cas_cap_vap(:)  = ros(:) * cp * max(4.,z2(:)) / psy(:)
    cas_cap_vap(:)  = ros(:) * cv * MAX(4.0_r8,z2(:)) / psy(:)
    cas_cap_co2(:)  =               MAX(4.0_r8,z2(:))  ! this goes 
    ! out to phosib

    !Bio approximate snow sfc vapor pressure and d(esnow)/dt with 
    !Bio ground surface values

    ets(:) = etg(:)
    bts(:) = btg(:)


    !     CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT
    !     LOSSES FROM CANOPY AND GROUND

    CALL NETRAD(radc3, &! INTENT(IN   ) :: radc3 (len,2)
         radt      , &! INTENT(OUT  ) :: radt  (len,3)
         stefan    , &! INTENT(IN   ) :: stefan
         fac1      , &! INTENT(OUT  ) :: fac1  (len)
         vcover    , &! INTENT(IN   ) :: vcover(len)
         thermk    , &! INTENT(IN   ) :: thermk(len)
         tc        , &! INTENT(IN   ) :: tc        (len)
         tg        , &! INTENT(IN   ) :: tg        (len)
         tice      , &! INTENT(IN   ) :: tf
         dtc4      , &! INTENT(OUT  ) :: dtc4  (len)
         dtg4      , &! INTENT(OUT  ) :: dtg4  (len)
         dts4      , &! INTENT(OUT  ) :: dts4  (len)
         closs     , &! INTENT(OUT  ) :: closs (len)
         gloss     , &! INTENT(OUT  ) :: gloss (len)
         sloss     , &! INTENT(OUT  ) :: sloss (len)
         tgeff     , &! INTENT(OUT  ) :: tgeff (len)
         areas     , &! INTENT(IN   ) :: areas (len)
         zlwup     , &! INTENT(OUT  ) :: tgeff (len)
         len         )! INTENT(IN   ) :: len

    DO I=1,len
       THgeff(I) = Tgeff(I) !/ BPS(I)
       esat12    = esat (tgeff(i))
       SHgeff(i) = qsat (esat12, PS(i)*100.0_r8)
    ENDDO
!    CALL VNQSAT(1       , &! INTENT(IN   ) :: iflag
!         tgeff   , &! INTENT(IN   ) :: TQS(IM)
!         PS      , &! INTENT(IN   ) :: PQS(IM)
!         SHgeff  , &! INTENT(OUT  ) :: QSS(IM)
!         len       )! INTENT(IN   ) :: IM
!    CALL QSAT_mix(SHgeff   , &!REAL(KIND=r8), intent(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturationeing processed by qSAT scheme.
!            tgeff      , &!REAL(KIND=r8), intent(in)  :: T(npnts)      !  Temperature (K).
!            PS*100.0_r8, &!REAL(KIND=r8), intent(in)  :: P(npnts)      !  Pressure (Pa).
!            len     , &!Integer, intent(in) :: npnts   !Points (=horizontal dimensions) being processed by qSAT scheme.
!            .FALSE.  )!logical, intent(in)  :: lq_mix  .true. return qsat as a mixing ratio
!       GET RESISTANCES FOR SIB
 
    CALL VNTLAT(grav , &  ! INTENT(IN        ) :: grav
         tice       , &  ! INTENT(IN        ) :: tice
         pr0        , &  ! INTENT(IN        ) :: pr0
         ribc       , &  ! INTENT(IN        ) :: ribc
         vkrmn      , &  ! INTENT(IN        ) :: vkrmn
         delta      , &  ! INTENT(IN        ) :: delta
         dt         , &  ! INTENT(IN        ) :: dtt
         tc         , &  ! INTENT(IN        ) :: tc   (len)
         tg         , &  ! INTENT(IN        ) :: tg   (len)
         ts         , &  ! INTENT(IN        ) :: ts   (len)
         ps         , &  ! INTENT(IN        ) :: ps   (len)
         zlt        , &  ! INTENT(IN        ) :: zlt  (len)
         wwwtem     , &  ! INTENT(IN        ) :: www  (len,3)
         tgs        , &  ! INTENT(IN        ) :: tgs  (len)
         etc        , &  ! INTENT(IN        ) :: etc  (len)
         etg        , &  ! INTENT(IN        ) :: etg  (len)
         snow       , &  ! INTENT(IN        ) :: snoww(len,2) ! snow cover (veg and ground) (m)
         rstfac     , &  ! INTENT(INOUT) :: rstfac(len,4)
         rsoil      , &  ! INTENT(IN        ) :: rsoil(len)
         hr         , &  ! INTENT(IN        ) :: hr   (len)
         wc         , &  ! INTENT(IN        ) :: wc   (len)
         wg         , &  ! INTENT(IN        ) :: wg   (len)
         !areas     , &  ! INTENT(IN        ) :: areas(len)
         snofac     , &  ! INTENT(IN        ) :: snofac
         qm         , &  ! INTENT(IN        ) :: sh   (len)
         z0         , &  ! INTENT(IN        ) :: z0   (len)
         spdm       , &  ! INTENT(IN        ) :: spdm (len)
         sha        , &  ! INTENT(IN        ) :: sha  (len)
         zb         , &  ! INTENT(IN        ) :: zb   (len)
         ros        , &  ! INTENT(IN        ) :: ros  (len) 
         cas_cap_co2, &  ! INTENT(IN        ) :: cas_cap_co2(len) 
         cu         , &  ! INTENT(OUT  ) :: cu   (len)
         ra         , &  ! INTENT(INOUT) :: ra   (len)
         thvgm      , &  ! INTENT(OUT  ) :: thvgm(len)
         rib        , &  ! INTENT(OUT  ) :: rib  (len)
         ustar      , &  ! INTENT(OUT  ) :: ustar(len) 
         ventmf     , &  ! INTENT(OUT  ) :: ventmf(len)
         thm        , &  ! INTENT(IN        ) :: thm  (len)
         tha        , &  ! INTENT(IN        ) :: tha  (len)
         z2         , &  ! INTENT(IN        ) :: z2   (len)
         d          , &  ! INTENT(IN        ) :: d    (len) 
         fc         , &  ! INTENT(OUT  ) :: fc   (len)
         fg         , &  ! INTENT(INOUT) :: fg   (len)
         rbc        , &  ! INTENT(IN        ) :: rbc  (len)
         rdc        , &  ! INTENT(IN        ) :: rdc  (len)
         gect       , &  ! INTENT(OUT  ) :: gect (len) 
         geci       , &  ! INTENT(OUT  ) :: geci (len) 
         gegs       , &  ! INTENT(OUT  ) :: gegs (len) 
         gegi       , &  ! INTENT(OUT  ) :: gegi (len)
         respcp     , &  ! INTENT(IN        ) :: respcp(len)
         rb         , &  ! INTENT(INOUT) :: rb   (len)
         rd         , &  ! INTENT(OUT  ) :: rd   (len) 
         rds        , &  ! INTENT(OUT  ) :: rds  (len) 
         !bps        , &  ! INTENT(IN        ) :: bps  (len)
         rst        , &  ! INTENT(IN        ) :: rst  (len)
         rc         , &  ! INTENT(OUT  ) :: rc   (len)
         ecmass     , &  ! INTENT(INOUT) :: ecmass(len)
         ea         , &  ! INTENT(IN        ) :: ea   (len) 
         !em         , &  ! INTENT(IN        ) :: em   (len)
         hrr        , &  ! INTENT(OUT  ) :: hrr  (len)
         assimn     , &  ! INTENT(OUT  ) :: assimn(len)
         bintc      , &  ! INTENT(OUT  ) :: bintc(len)
         ta         , &  ! INTENT(IN        ) :: ta   (len) 
         pco2m      , &  ! INTENT(IN        ) :: pco2m(len)
         po2m       , &  ! INTENT(IN        ) :: po2m
         vmax0      , &  ! INTENT(IN        ) :: vmax0(len)
         green      , &  ! INTENT(IN        ) :: green(len)
         tran       , &  ! INTENT(IN        ) :: tran (len,2,2) 
         ref        , &  ! INTENT(IN        ) :: ref  (len,2,2)
         gmudmu     , &  ! INTENT(IN        ) :: gmudmu(len)
         trop       , &  ! INTENT(IN        ) :: trop (len)
         trda       , &  ! INTENT(IN        ) :: trda (len) 
         trdm       , &  ! INTENT(IN        ) :: trdm (len)
         slti       , &  ! INTENT(IN        ) :: slti (len)
         shti       , &  ! INTENT(IN        ) :: shti (len)
         hltii      , &  ! INTENT(IN        ) :: hltii(len) 
         hhti       , &  ! INTENT(IN        ) :: hhti (len)
         radn       , &  ! INTENT(IN        ) :: radn (len,2,2)
         effcon     , &  ! INTENT(IN        ) :: effcon(len)
         binter     , &  ! INTENT(IN        ) :: binter(len) 
         gradm      , &  ! INTENT(IN        ) :: gradm (len)
         atheta     , &  ! INTENT(IN        ) :: atheta(len)
         btheta     , &  ! INTENT(IN        ) :: btheta(len)
         aparkk     , &  ! INTENT(OUT  ) :: aparkk(len)
         wsfws      , &  ! INTENT(OUT  ) :: wsfws (len)
         wsfht      , &  ! INTENT(OUT  ) :: wsfht (len)
         wsflt      , &  ! INTENT(OUT  ) :: wsflt (len) 
         wci        , &  ! INTENT(OUT  ) :: wci   (len)
         whs        , &  ! INTENT(OUT  ) :: whs   (len)
         omepot     , &  ! INTENT(OUT  ) :: omepot(len)
         assimpot   , &  ! INTENT(OUT  ) :: assimpot(len) 
         assimci    , &  ! INTENT(OUT  ) :: assimci(len)
         antemp     , &  ! INTENT(OUT  ) :: antemp (len)
         assimnp    , &  ! INTENT(OUT  ) :: assimnp(len)
         wags       , &  ! INTENT(OUT  ) :: wags(len) 
         wegs       , &  ! INTENT(OUT  ) :: wegs(len)
         aparc      , &  ! INTENT(IN        ) :: aparc(len)
         pfd        , &  ! INTENT(OUT  ) :: pfd(len)
         assim      , &  ! INTENT(OUT  ) :: assim(len)
         td         , &  ! INTENT(IN        ) :: td(len,nsoil) 
         wopt       , &  ! INTENT(IN        ) :: wopt(len)
         zm         , &  ! INTENT(IN        ) :: zm  (len)
         wsat       , &  ! INTENT(IN        ) :: wsat(len)
         soilscale  , &  ! INTENT(INOUT) :: soilscale(len,nsoil+1) 
         zmstscale  , &  ! INTENT(INOUT) :: zmstscale(len,2)
         zltrscale  , &  ! INTENT(OUT  ) :: zltrscale(len)
         zmlscale   , &  ! INTENT(OUT  ) :: zmlscale(len) 
         drst       , &  ! INTENT(OUT  ) :: drst(len)
         soilq10    , &  ! INTENT(INOUT) :: soilq10(len,nsoil+1)
         ansqr      , &  ! INTENT(OUT  ) :: ansqr(len) 
         soilscaleold,&  ! INTENT(OUT  ) :: soilscaleold(len)
         !nsib       , &  ! INTENT(IN        ) :: nsib
         len        , &  ! INTENT(IN        ) :: len
         nsoil      , &  ! INTENT(IN        ) :: nsoil
         forcerestore,&  ! INTENT(IN        ) :: forcerestore
         dotkef     , &  ! INTENT(IN        ) :: dotkef
         thgeff     , &  ! INTENT(IN        ) :: thgeff(len)
         shgeff     , &  ! INTENT(IN        ) :: shgeff(len)
         tke        , &  ! INTENT(IN        ) :: tke(len)
         ct         , &  ! INTENT(OUT  ) :: ct(len) 
         louis      , &  ! INTENT(IN        ) :: louis
         !zwind      , &  ! INTENT(IN        ) :: zwind
         !ztemp      , &  ! INTENT(IN        ) :: ztemp
         respg      , &  ! INTENT(INOUT) :: respg(len) 
         respfactor , &  ! INTENT(IN        ) :: respfactor(len, nsoil+1)
         pco2ap     , &  ! INTENT(INOUT) :: pco2ap(len)
         pco2i      , &  ! INTENT(OUT  ) :: pco2i(len)   ! added for neil suits' programs
         pco2c      , &  ! INTENT(OUT  ) :: pco2c(len)  ! more added nsuits vars
         pco2s      , &  ! INTENT(OUT  ) :: pco2s(len)
         co2cap     , &  ! INTENT(OUT  ) :: co2cap(len)  ! moles of air in canopy
         cflux        )  ! INTENT(OUT  ) :: cflux(len)


    !   this call for ustar, cu for oceanic value of z0 
    IF(louis == 1) THEN
       DO i = 1,len
!          zzwind(i) = z2(i) - d(i) + zb (i)
!          zztemp(i) = z2(i) - d(i) + zb (i)
          zzwind(i) = d(i) + zb (i)
          zztemp(i) = d(i) + zb (i)
       ENDDO
       CALL VMFCALZO(&
            !z0       , & ! INTENT(IN   ) :: pr0 ! is not used
            !PR0      , & ! INTENT(IN   ) :: pr0 ! is not used
            !RIBC     , & ! INTENT(IN   ) :: ribc! is not used
            VKRMN    , & ! INTENT(IN   ) :: vkrmn
            !DELTA    , & ! INTENT(IN   ) :: delta! is not used
            GRAV     , & ! INTENT(IN   ) :: grav
            !PS       , & ! INTENT(IN   ) :: PS        (len)! is not used   
            tha      , & ! INTENT(IN   ) :: THa        (len)
            SPDM     , & ! INTENT(IN   ) :: SPDM        (len)
            !ROS      , & ! INTENT(IN   ) :: ROS        (len)! is not used   
            CUo      , & ! INTENT(OUT  ) :: CU        (len)
            THVGM    , & ! INTENT(IN   ) :: THVGM        (len)
            RIB      , & ! INTENT(OUT  ) :: RIB        (len) 
            USTARo   , & ! INTENT(OUT  ) :: USTAR        (len)
            zzwind   , & ! INTENT(IN   ) :: zzwind        (len)
            zztemp   , & ! INTENT(IN   ) :: zztemp        (len)
            len        ) ! INTENT(IN   ) :: len

      DO I=1,len
         DRAG (I,1) = ROS(I) * CU(I) * USTAR(I)
         DRAG (I,2) = ROS(I) * CUo(I) * USTARo(I)
         ANWTC(i)   = ANTEMP(i)* TC(i)
      END DO


    ELSE  IF(louis == 2) THEN                                 
       CALL VMFCALo(&
            !z0       , & ! INTENT(IN   ) :: pr0 ! is not used
            !PR0      , &! INTENT(IN   ) :: pr0
            RIBC     , &! INTENT(IN   ) :: ribc
            VKRMN    , &! INTENT(IN   ) :: vkrmn
            GRAV     , &! INTENT(IN   ) :: grav               
            SPDM     , &! INTENT(IN   ) :: SPDM   (len)
            ZB       , &! INTENT(IN   ) :: ZB     (len)
            !ROS      , &! INTENT(IN   ) :: ROS    (len)     
            CUo      , &! INTENT(OUT  ) :: CU     (len)
            THVGM    , &! INTENT(IN   ) :: THVGM  (len)
            USTARo   , &! INTENT(OUT  ) :: USTAR  (len)
            tha      , &! INTENT(IN   ) :: tha    (len)
            len      , &! INTENT(IN   ) :: len
            thgeff   , &! INTENT(IN   ) :: thgeff (len)
            tke      , &! INTENT(IN   ) :: tke    (len)
            dotkef )    ! INTENT(IN   ) :: dotkef   
    DO I=1,len
       DRAG(I,1) = ROS(I) * CU(I) * USTAR(I)
       DRAG(I,2) = ROS(I) * CUo(I) * USTARo(I)
       ANWTC (i) = ANTEMP(i)* TC(i)
    END DO


    ELSE  IF(louis == 3) THEN    
       !
       !     the first call to vntlat just gets the neutral values of ustar
       !     and ventmf.
       !
       jstneu=.TRUE.
       CALL vntlax(        &
                   USTARo, &!INTENT(inout) :: ustarn(ncols)
                   !bps   , &!INTENT(in   ) :: bps   (ncols)
                   zb    , &!INTENT(in   ) :: dzm   (ncols)
                   ROS   , &!INTENT(IN   ) :: ROS   (len)
                   THVGM , &!INTENT(in   ) :: dzm   (ncols)
                   VENTMF, &!INTENT(OUT  ) :: VENTMF(len)
                   RIB   , &!INTENT(OUT  ) :: RIB   (len)
                   cu    , &!INTENT(inout) :: cu    (ncols)
                   ct    , &!INTENT(inout) :: ct    (ncols)
                   cti   , &!INTENT(inout) :: cti    (ncols)
                   cui   , &!INTENT(OUT  ) :: CUI   (len)
                   cuni  , &!INTENT(inout) :: cuni  (ncols)
                   ctni  , &!INTENT(inout) :: ctni  (ncols)
                   ustar , &!INTENT(inout) :: ustar (ncols)
                   ra    , &!INTENT(inout) :: ra    (ncols)
                   ta    , &!INTENT(in   ) :: ta    (ncols)
                   u2    , &!INTENT(inout) :: u2    (ncols)
                   THM   , &!INTENT(in   ) :: tm  (ncols)
                   SPDM  , &!INTENT(in   ) :: um  (ncols)
                   d     , &!INTENT(in   ) :: d     (ncols)
                   z0    , &!INTENT(in   ) :: z0    (ncols)
                   z2    , &!INTENT(in   ) :: z2    (ncols)
                   len   , &!INTENT(in   ) :: nmax
                   jstneu, &!INTENT(in   ) :: jstneu
                   len     )!INTENT(in   )
        jstneu=.FALSE.
       CALL vntlax(        &
                   USTARo, &!INTENT(inout) :: ustarn(ncols)
                   !bps   , &!INTENT(in   ) :: bps   (ncols)
                   zb    , &!INTENT(in   ) :: dzm   (ncols)
                   ROS   , &!INTENT(IN   ) :: ROS   (len)
                   THVGM , &!INTENT(in   ) :: dzm   (ncols)
                   VENTMF, &!INTENT(OUT  ) :: VENTMF(len)
                   RIB   , &!INTENT(OUT  ) :: RIB   (len)
                   cu    , &!INTENT(inout) :: cu    (ncols) 
                   ct    , &!INTENT(inout) :: cu    (ncols)
                   cti   , &!INTENT(inout) :: cti    (ncols)
                   cui   , &!INTENT(OUT  ) :: CUI   (len)
                   cuni  , &!INTENT(inout) :: cuni  (ncols)
                   ctni  , &!INTENT(inout) :: ctni  (ncols)
                   ustar , &!INTENT(inout) :: ustar (ncols)
                   ra    , &!INTENT(inout) :: ra    (ncols)
                   ta    , &!INTENT(in   ) :: ta    (ncols)
                   u2    , &!INTENT(inout) :: u2    (ncols)
                   THM   , &!INTENT(in   ) :: tm  (ncols)
                   SPDM  , &!INTENT(in   ) :: um  (ncols)
                   d     , &!INTENT(in   ) :: d     (ncols)
                   z0    , &!INTENT(in   ) :: z0    (ncols)
                   z2    , &!INTENT(in   ) :: z2    (ncols)
                   len   , &!INTENT(in   ) :: nmax
                   jstneu, &!INTENT(in   ) :: jstneu
                   len     )!INTENT(in   )
    DO I=1,len
       DRAG(I,1) = ROS(I) * CU(I) * USTAR(I)
       DRAG(I,2) = ROS(I) * CU(I) * USTARo(I)
       ANWTC (i) = ANTEMP(i)* TC(i)
    END DO


    END IF

    !itb...calculate partial derivatives of the various heat fluxes
    !itb...with respect to ground/canopy/snow temp, as well as
    !itb...some other derivatives.


    CALL DELLWF(&
         !DT, & ! INTENT(IN   ) :: dta ! is not used
         dtc4    , & ! INTENT(IN   ) :: dtc4 (len)
         dtg4    , & ! INTENT(IN   ) :: dtg4 (len)
         dts4    , & ! INTENT(IN   ) :: dts4 (len)
         fac1    , & ! INTENT(IN   ) :: fac1 (len)
         areas   , & ! INTENT(IN   ) :: areas(len)
         lcdtc   , & ! INTENT(OUT  ) :: lcdtc(len)
         lcdtg   , & ! INTENT(OUT  ) :: lcdtg(len)
         lcdts   , & ! INTENT(OUT  ) :: lcdts(len)
         lgdtg   , & ! INTENT(OUT  ) :: lgdtg(len)
         lgdtc   , & ! INTENT(OUT  ) :: lgdtc(len)
         lsdts   , & ! INTENT(OUT  ) :: lsdts(len)
         lsdtc   , & ! INTENT(OUT  ) :: lsdtc(len)
         len       ) ! INTENT(IN   ) :: len


    CALL DELHF( DT      , & ! INTENT(IN   ) :: dta
         CP      , & ! INTENT(IN   ) :: cp
         !bps     , & ! INTENT(IN   ) :: bps   (len)
         ts      , & ! INTENT(IN   ) :: tm    (len)
         tgs     , & ! INTENT(IN   ) :: tg    (len)
         !tsnow   , & ! INTENT(IN   ) :: ts    (len) !is not used
         tc      , & ! INTENT(IN   ) :: tc    (len)
         ta      , & ! INTENT(IN   ) :: ta    (len)
         ros     , & ! INTENT(IN   ) :: ros   (len)
         ra      , & ! INTENT(IN   ) :: ra    (len)
         rb      , & ! INTENT(IN   ) :: rb    (len)
         rd      , & ! INTENT(IN   ) :: rd    (len)
         HCDTC   , & ! INTENT(OUT  ) :: HCDTC (len)
         HCDTA   , & ! INTENT(OUT  ) :: HCDTA (len)
         HGDTG   , & ! INTENT(OUT  ) :: HGDTG (len)
         HGDTA   , & ! INTENT(OUT  ) :: HGDTA (len)
         HSDTS   , & ! INTENT(OUT  ) :: HSDTS (len)
         HSDTA   , & ! INTENT(OUT  ) :: HSDTA (len)
         HADTA   , & ! INTENT(OUT  ) :: HADTA (len)
         HADTH   , & ! INTENT(OUT  ) :: HADTH (len)
         hc      , & ! INTENT(OUT  ) :: hc    (len)
         hg      , & ! INTENT(OUT  ) :: hg    (len) 
         hs      , & ! INTENT(OUT  ) :: fss   (len)
         fss     , & ! INTENT(OUT  ) :: hs    (len)
         len       ) ! INTENT(IN   ) :: len


    CALL DELEF( DT         , & ! INTENT(IN   ) :: dta
         CP         , & ! INTENT(IN   ) :: cp
         !ps         , & ! INTENT(IN   ) :: ps    (len)! is not used
         em         , & ! INTENT(IN   ) :: em    (len)
         ea         , & ! INTENT(IN   ) :: ea    (len)
         ros        , & ! INTENT(IN   ) :: ros   (len)
         HRr        , & ! INTENT(OUT  ) :: hrr   (len)
         !fc         , & ! INTENT(IN   ) :: fc    (len)! is not used
         fg         , & ! INTENT(IN   ) :: fg    (len)
         ra         , & ! INTENT(IN   ) :: ra    (len)
         !rb         , & ! INTENT(IN   ) :: rb    (len)
         rd         , & ! INTENT(IN   ) :: rd    (len)
         !rc         , & ! INTENT(IN   ) :: rc    (len)! is not used
         rsoil      , & ! INTENT(IN   ) :: rsoil (len)
         !snow       , & ! INTENT(IN   ) :: snow  (len,2)! is not used
         !capac      , & ! INTENT(IN   ) :: capac (len,2)! is not used
         !www        , & ! INTENT(IN   ) :: www   (len,3) ! is not used
         ECDTC      , & ! INTENT(OUT  ) :: ECDTC (len)
         ECDEA      , & ! INTENT(OUT  ) :: ECDEA (len)
         EGDTG      , & ! INTENT(OUT  ) :: EGDTG (len)
         EGDEA      , & ! INTENT(OUT  ) :: EGDEA (len)
         ESDTS      , & ! INTENT(OUT  ) :: ESDTS (len)
         ESDEA      , & ! INTENT(OUT  ) :: ESDEA (len)
         EADEA      , & ! INTENT(OUT  ) :: EADEA (len)
         EADEM      , & ! INTENT(OUT  ) :: EADEM (len)
         ec         , & ! INTENT(OUT  ) :: ec    (len)
         eg         , & ! INTENT(OUT  ) :: eg    (len)
         es         , & ! INTENT(OUT  ) :: es    (len)
         fws        , & ! INTENT(OUT  ) :: fws   (len)
       !  hltm       , & ! INTENT(IN   ) :: hltm            ! is not used
         etc        , & ! INTENT(IN   ) :: etc   (len)
         etg        , & ! INTENT(IN   ) :: etg   (len)
         btc        , & ! INTENT(IN   ) :: btc   (len)
         btg        , & ! INTENT(IN   ) :: btg   (len)
        ! bts        , & ! INTENT(IN   ) :: bts   (len) ! is not used   
        ! areas      , & ! INTENT(IN   ) :: areas (len)! is not used   
         gect       , & ! INTENT(IN   ) :: gect  (len)
         geci       , & ! INTENT(IN   ) :: geci  (len)
         gegs       , & ! INTENT(IN   ) :: gegs  (len)
         gegi       , & ! INTENT(IN   ) :: gegi  (len)
         psy        , & ! INTENT(IN   ) :: psy   (len)
         snofac     , & ! INTENT(IN   ) :: snofac
         hr         , & ! INTENT(IN   ) :: hr    (len)
         len          ) ! INTENT(IN          ) :: len

    !     get soil thermal properties
    IF(.NOT.forcerestore)THEN
         CALL soilprop( td       , &! INTENT(IN   ) :: td        (len,nsoil)
         tgs      , &! INTENT(IN   ) :: tg        (len)
         slamda   , &! INTENT(OUT  ) :: slamda        (len,nsoil)
         shcap    , &! INTENT(OUT  ) :: shcap        (len,nsoil)
         wwwtem   , &! INTENT(IN   ) :: www        (len,3)
         poros    , &! INTENT(IN   ) :: poros        (len)
         ztdep    , &! INTENT(IN   ) :: ztdep        (len,nsoil)
         asnow    , &! INTENT(IN   ) :: asnow
         snow(1:len,2), &! INTENT(IN   ) :: snoww        (len)
         areas    , &! INTENT(IN   ) :: areas        (len)
         tice     , &! INTENT(IN   ) :: tf
         snomel   , &! INTENT(IN   ) :: snomel
         sandfrac , &! INTENT(IN   ) :: sandfrac (len)!nor used
         clayfrac , &! INTENT(IN   ) :: sandfrac (len)!nor used
         !nsib     , &! INTENT(IN   ) :: nsib!nor used
         len      , &! INTENT(IN   ) :: len
         nsoil      )! INTENT(IN   ) :: nsoil
    ELSE
       CALL soilprop( td       , &! INTENT(IN   ) :: td        (len,nsoil)
         tgs      , &! INTENT(IN   ) :: tg        (len)
         slamda   , &! INTENT(OUT  ) :: slamda        (len,nsoil)
         shcap    , &! INTENT(OUT  ) :: shcap        (len,nsoil)
         wwwtem   , &! INTENT(IN   ) :: www        (len,3)
         poros    , &! INTENT(IN   ) :: poros        (len)
         ztdep    , &! INTENT(IN   ) :: ztdep        (len,nsoil)
         asnow    , &! INTENT(IN   ) :: asnow
         snow(1:len,2), &! INTENT(IN   ) :: snoww        (len)
         areas    , &! INTENT(IN   ) :: areas        (len)
         tice     , &! INTENT(IN   ) :: tf
         snomel   , &! INTENT(IN   ) :: snomel
         sandfrac , &! INTENT(IN   ) :: sandfrac (len)!nor used
         clayfrac , &! INTENT(IN   ) :: sandfrac (len)!nor used
         !nsib     , &! INTENT(IN   ) :: nsib
         len      , &! INTENT(IN   ) :: len
         nsoil      )! INTENT(IN   ) :: nsoil
    END IF

    !     check against new code, rn derivatives may be zero
    CALL SIBSLV( DT             , & ! INTENT(IN   ) :: dt
         GRAV2          , & ! INTENT(IN   ) :: grav2
         CP             , & ! INTENT(IN   ) :: cp
!         HLTM           , & ! INTENT(IN   ) :: hltm       !is not used
         tgs            , & ! INTENT(IN   ) :: tg    (len)
!         tsnow          , & ! INTENT(IN   ) :: tsnow (len)!is not used
         td(1:len,nsoil)    , & ! INTENT(IN   ) :: td    (len)
         slamda(1:len,nsoil), & ! INTENT(IN   ) :: slamda(len)
         pi             , & ! INTENT(IN   ) :: pi
         areas          , & ! INTENT(IN   ) :: areas (len)
!         fac1           , & ! INTENT(IN   ) :: fac1  (len)!is not used
         VENTMF         , & ! INTENT(IN   ) :: VENTMF(len)
         PSB            , & ! INTENT(IN   ) :: PSB   (len)
!         BPS            , & ! INTENT(IN   ) :: BPS   (len)!is not used
!         ros            , & ! INTENT(IN   ) :: ros   (len)!is not used
         psy            , & ! INTENT(IN   ) :: psy   (len)    
         czh            , & ! INTENT(IN   ) :: cg    (len)
         czc            , & ! INTENT(IN   ) :: ccx   (len)
         cas_cap_heat   , & ! INTENT(IN   ) :: cas_cap_heat(len)
         cas_cap_vap    , & ! INTENT(IN   ) :: cas_cap_vap (len)
         lcdtc          , & ! INTENT(IN   ) :: lcdtc  (len)
         lcdtg          , & ! INTENT(IN   ) :: lcdtg  (len)
         lcdts          , & ! INTENT(IN   ) :: lcdts  (len)
         lgdtg          , & ! INTENT(IN   ) :: lgdtg  (len)
         lgdtc          , & ! INTENT(IN   ) :: lgdtc  (len)
         lsdts          , & ! INTENT(IN   ) :: lsdts  (len)
         lsdtc          , & ! INTENT(IN   ) :: lsdtc  (len)
         HCDTC          , & ! INTENT(IN   ) :: HCDTC  (len)
         HCDTA          , & ! INTENT(IN   ) :: HCDTA  (len)
         HGDTG          , & ! INTENT(IN   ) :: HGDTG  (len)
         HGDTA          , & ! INTENT(IN   ) :: HGDTA  (len)
         HSDTS          , & ! INTENT(IN   ) :: HSDTS  (len)
         HSDTA          , & ! INTENT(IN   ) :: HSDTA  (len)
         HADTA          , & ! INTENT(IN   ) :: HADTA  (len)
         HADTH          , & ! INTENT(IN   ) :: HADTH  (len)
         hc             , & ! INTENT(IN   ) :: hc     (len)
         hg             , & ! INTENT(IN   ) :: hg     (len)
         hs             , & ! INTENT(IN   ) :: hs     (len)
         fss            , & ! INTENT(IN   ) :: FSS    (len)
         ECDTC          , & ! INTENT(IN   ) :: ECDTC  (len)
         ECDEA          , & ! INTENT(IN   ) :: ECDEA  (len)
         EGDTG          , & ! INTENT(IN   ) :: EGDTG  (len)
         EGDEA          , & ! INTENT(IN   ) :: EGDEA  (len)
         ESDTS          , & ! INTENT(IN   ) :: ESDTS  (len)
         ESDEA          , & ! INTENT(IN   ) :: ESDEA  (len)
         EADEA          , & ! INTENT(IN   ) :: EADEA  (len)
         EADEM          , & ! INTENT(IN   ) :: EADEM  (len)
         ec             , & ! INTENT(IN   ) :: ec     (len)
         eg             , & ! INTENT(IN   ) :: eg     (len)
         es             , & ! INTENT(IN   ) :: es     (len)
         fws            , & ! INTENT(IN   ) :: FWS    (len)
!         etc            , & ! INTENT(IN   ) :: etc    (len) ! is not used
!         etg            , & ! INTENT(IN   ) :: etg    (len) ! is not used
!         ets            , & ! INTENT(IN   ) :: ets    (len)! is not used
!         btc            , & ! INTENT(IN   ) :: btc    (len)! is not used
!         btg            , & ! INTENT(IN   ) :: btg    (len)! is not used
!         bts            , & ! INTENT(IN   ) :: bts    (len)! is not used
         RADT           , & ! INTENT(IN   ) :: RADT   (len,3)
         dtc            , & ! INTENT(OUT  ) :: dtc    (len)
         dtg            , & ! INTENT(OUT  ) :: dtg    (len,2)
         dth            , & ! INTENT(OUT  ) :: dth    (len)
         dqm            , & ! INTENT(OUT  ) :: dqm    (len)
         dta            , & ! INTENT(OUT  ) :: dta    (len)
         dea            , & ! INTENT(OUT  ) :: dea    (len)
         len            , & ! INTENT(IN   ) :: len
         sibdrv           ) ! INTENT(IN   ) :: sibdrv
!         forcerestore     ) ! INTENT(IN   ) :: forcerestore! is not used


    DO i = 1,len
       radt(i,2) = (1.0_r8-areas(i))*radt(i,2)+areas(i)*radt(i,3)
    ENDDO

    !         call endtem25(etc,etg, ts, ta, ea, em, 
    !     *      btc, btg, bts, 
    !     *      deadtg, deadtc, dtc, dtg, dta, dea, 
    !     *      deadqm, dqm, wc, wg, ra, rb, rd, dt, ros, psy, snow, capac,
    !     *      hltm, areas, asnow, snofac, csoil, hr, rst, fc, fg,
    !     *      rsoil, zdepth, wwwtem, czc, czh, cas_cap_heat,cas_cap_vap,
    !     *      hc, hg, fss, fws, ec, eg, eci, ect,
    !     *      egi, egs, chf, shf, heaten, tgs, td(1,nsoil), 
    !     *      hcdtc, hcdtg, hcdta,
    !     *      hcdth, hgdtc, hgdtg, hgdta, hgdth, dth, tg, tice, cp, 
    !     *      slamda(1,nsoil), len, nsib, geci, gect, gegi, gegs, 
    !     *      forcerestore )


    !        print*,'call updat2'
    CALL updat2(&
         snow        , &!INTENT(inout) :: snoww(nsib,2) ! snow-interception (1-veg, 2-ground) (meters)
         capac       , &!INTENT(inout) :: capac(nsib,2) ! liquid interception
         snofac      , &!INTENT(in   ) :: snofac        !  ___(lat ht of vap)___     (unitless)
         ect         , &!INTENT(out  ) :: ect  (len)    ! transpiration flux (J m^-2 for the timestep)
         eci         , &!INTENT(out  ) :: eci  (len)    ! interception flux (veg - CAS) (J m^-2)
         egi         , &!INTENT(out  ) :: egi  (len)    ! ground interception flux (J m^-2)
         egs         , &!INTENT(out  ) :: egs  (len)    ! ground evaporative flux (J m^-2)
         hltm        , &!INTENT(IN   ) :: hltm       !is not used
         wwwtem      , &!INTENT(inout) :: www  (len,3)  ! soil moisture (% of saturation)
         pi          , &!INTENT(in   ) :: pi              ! 3.1415...
         czh         , &!INTENT(in   ) :: cg   (len) ! surface layer heat capacity (J m^-2 deg^-1)
         dtd         , &!INTENT(out  ) :: dtd(len,nsoil)  ! deep soil temperature increment (K)
         dtg         , &!INTENT(inout) :: dtg  (len,2)  ! delta ground surface temperature (K)
         dtc         , &!INTENT(inout) :: dtc  (len)    ! delta vegetation temperature (K)
         ta          , &!INTENT(in   ) :: ta   (len)     ! CAS temperature (K)
         dta         , &!INTENT(in   ) :: dta  (len) ! delta canopy air space (CAS) temperature (K)
         dea         , &!INTENT(in   ) :: dea  (len) ! delta CAS vapor pressure (Pa)
         dt          , &!INTENT(in) :: dtt              ! timestep (seconds)
         roff        , &!INTENT(inout) :: roff (len)    ! runoff (mm)
         tc          , &!INTENT(in) :: tc        (len) ! canopy temperature (K)
         td          , &!INTENT(in) :: td   (nsib,nsoil) ! deep soil temperature (K))
         tg          , &!INTENT(in) :: tg        (len) ! ground surface temperature (K)
         bee         , &!INTENT(in) :: bee     (len) ! Clapp&Hornberger 'b' exponent
         poros       , &!INTENT(in) :: poros (len) ! soil porosity (fraction)
         satco       , &!INTENT(in) :: satco (len) ! hydraulic conductivity at saturation (UNITS?)
         slope       , &!INTENT(in) :: slope (len) ! cosine of mean slope
         phsat       , &!INTENT(in) :: phsat (len) ! soil tension at saturation  (UNITS?)
         zdepth      , &!INTENT(in) :: zdepth(nsib,3) ! soil layer depth * porosity (meters)
         ecmass      , &!INTENT(out) :: ecmass (len)    ! canopy evaporation (mm)
         egmass      , &!INTENT(out) :: egmass (len)    ! ground evaporation (mm)
         shf         , &!INTENT(out) :: shf         (len)    ! soil heat flux (W m^-2)
         tice        , &!INTENT(in) :: tf              ! freezing temperature (273.16 K)
         snomel      , &!INTENT(in) :: snomel      ! latent heat of fusion for ice (J m^-3)
         asnow       , &!INTENT(in) :: asnow       ! conversion factor for kg water to depth of snow (16.7)
         czc         , &!INTENT(in) :: ccx        (len) ! canopy heat capacity (J m^-2 deg^-1)
         csoil       , &!INTENT(in) :: csoil (len) ! soil heat capacity (J m^-2 deg^-1)
         chf         , &!INTENT(out) :: chf         (len)    ! canopy heat flux (W m^-2)
         hc          , &!INTENT(inout) :: hc   (len)    ! canopy sensible heat flux (J m^-2)
         hg          , &!INTENT(inout) :: hg   (len)    ! ground sensible heat flux (J m^-2)
         areas       , &!INTENT(inout) :: areas(len)    ! fractional snow coverage (0 to 1)
         q3l         , &!INTENT(out) :: q3l         (len)    ! 'Liston' drainage from the bottom of soil layer 3 (mm)
         q3o         , &!INTENT(out) :: q3o         (len)    ! gravitational drainage out of soil layer 3 (mm)
         qqq         , &!INTENT(out) :: QQQ(len,3)      ! soil layer drainage (mm m^-2 timestep)
         zmelt1      , &!INTENT(out) :: zmelt  (len)    ! depth of melted water (m)
!         cww         , &!INTENT(in) :: cw              ! water heat capacity (J m^-3 deg^-1)
         len         , &!INTENT(IN   ) :: len
         nsib        , &!INTENT(IN   ) :: nsib
         nsoil       , &!INTENT(IN   ) :: nsoil
         forcerestore, &!INTENT(IN   ) :: forcerestore
         etc         , &!INTENT(in) :: etc        (len) ! 'E-star' of the canopy - vapor pressure (Pa)
         ea          , &!INTENT(in) :: ea        (len) ! CAS vapor pressure (Pa)
         btc         , &!INTENT(in) :: btc        (len) ! d(E(Tc))/d(Tc) - Clausius-Clapyron
         geci        , &!INTENT(in) :: gegi  (len) ! wet fraction of ground/Rd
         ros         , &!INTENT(in) :: ros        (len) ! air density (kg m^-3)
         cp          , &!INTENT(in) :: cp              ! specific heat of air at const pres (J kg-1 deg-1)
         psy         , &!INTENT(in) :: psy        (len) ! 
         gect        , &!INTENT(in) :: gect  (len) ! dry fraction of canopy/(Rst + 2Rb)
         etg         , &!INTENT(in) :: etg        (len) ! 'E-star' of the ground surface  (Pa)
         btg         , &!INTENT(in) :: btg        (len) ! d(E(tg))/d(Tg) - Clausius-Clapyron
         gegs        , &!INTENT(in) :: gegs  (len) ! dry fraction of ground/(fg*rsoil + Rd)
         hr          , &!INTENT(in) :: hr        (len) ! soil surface relative humidity
         fg          , &!INTENT(in) :: fg        (len) ! flag indicating direction of vapor pressure
         !                           deficit between CAS and ground: 0 => ea>e(Tg)
         !                                                           1 => ea<e(Tg)
         gegi        , &!INTENT(in) :: gegi  (len) ! wet fraction of ground/Rd
         rd          , &!INTENT(in) :: rd      (len) ! ground-CAS resistance
         rb          , &!INTENT(in) :: rb      (len) ! leaf-CAS resistance
!         hcdtc       , &!INTENT(in) dHc/dTc
!         hcdta       , &!INTENT(in) dHc/dTa
!         hgdta       , &!INTENT(in) dHg/dTa
         slamda(1:len,nsoil))!INTENT(in) :: slamda(len) !  

    DO i = 1, len
       EVT(i) = (55.56_r8 *dtinv) * (ECMASS(i) + EGMASS(i))
       wegs(i) = wegs(i) * evt(i)
    ENDDO

    !     get soil temperature increments
    IF(.NOT.forcerestore) &
         CALL soiltherm(&
         td        ,& ! INTENT(IN  ) :: td    (nsib,nsoil)
         dtd       ,& ! INTENT(OUT ) :: dtd   (len,nsoil)
         tgs       ,& ! INTENT(IN  ) :: tgs   (len)
         dtg       ,& ! INTENT(IN  ) :: dtg   (len)
         slamda    ,& ! INTENT(IN  ) :: slamda(len,nsoil)
         shcap     ,& ! INTENT(IN  ) :: shcap (len,nsoil)
         !ztdep     ,& ! INTENT(IN  ) :: ztdep (len,nsoil) !is not used
         dt        ,& ! INTENT(IN  ) :: dt
         nsib      ,& ! INTENT(IN  ) :: nsib
         len       ,& ! INTENT(IN  ) :: len
         nsoil      ) ! INTENT(IN  ) :: nsoil


    !    update prognostic variables, get total latent and sensible fluxes
    DO i = 1,len
       thmtem(i) = thm(i)
       shtem(i) = qm(i)
    ENDDO

    CALL addinc(&
!        grav2       , & ! INTENT(IN   ) :: grav2!is not used
         cp          , & ! INTENT(IN   ) :: cp
!         dt          , & ! INTENT(IN   ) :: dtt  !is not used
!         hc          , & ! INTENT(IN   ) :: hc(len)!is not used
!         hg          , & ! INTENT(IN   ) :: hg(len)!is not used
         ps          , & ! INTENT(IN   ) :: ps(len)
!         bps         , & ! INTENT(IN   ) :: bps(len)!is not used
!         ecmass      , & ! INTENT(IN   ) :: ecmass(len)!is not used
         psy         , & ! INTENT(IN   ) :: psy(len)
         ros         , & ! INTENT(IN   ) :: rho(len)
         hltm        , & ! INTENT(IN   ) :: hltm
!         egmass      , & ! INTENT(IN   ) :: egmass(len)!is not used
         fss         , & ! INTENT(OUT  ) :: fss(len)
         fws         , & ! INTENT(OUT  ) :: fws(len)
         hflux       , & ! INTENT(OUT  ) :: hflux(len)
         etmass      , & ! INTENT(OUT  ) :: etmass(len)
!         psb         , & ! INTENT(IN   ) :: psb(len)!is not used
         td          , & ! INTENT(INOUT) :: td(nsib,nsoil)
!         thmtem      , & ! INTENT(IN   ) :: thm(len)!is not used
         ts          , & ! INTENT(IN   ) :: ts(len)
         shtem       , & ! INTENT(OUT  ) :: sh(len)
         tc          , & ! INTENT(INOUT) :: tc(len)
         tg          , & ! INTENT(INOUT) :: tg(len)
         ta          , & ! INTENT(INOUT) :: ta(len)
         ea          , & ! INTENT(INOUT) :: ea(len)
         ra          , & ! INTENT(IN   ) :: ra(len)
         em          , & ! INTENT(IN   ) :: em(len)
         sha         , & ! INTENT(INOUT) :: sha(len)
         dtd         , & ! INTENT(IN   ) :: dtd(len,nsoil)
         dtc         , & ! INTENT(IN   ) :: dtc(len)
         dtg         , & ! INTENT(IN   ) :: dtg(len)
         dta         , & ! INTENT(IN   ) :: dta(len)
         dea         , & ! INTENT(IN   ) :: dea(len)
         drst        , & ! INTENT(IN   ) :: DRST(len)
         rst         , & ! INTENT(INOUT) :: rst(len)
         bintc       , & ! INTENT(IN   ) :: bintc(len)
         itype       , & ! INTENT(IN   ) :: bintc(len)
         len         , & ! INTENT(IN   ) :: len
         nsib        , & ! INTENT(IN   ) :: nsib
         nsoil         ) ! INTENT(IN   ) :: nsoil


    CALL vnqsat(2         , &! INTENT(IN        ) :: iflag
         ta         , &! INTENT(IN        ) :: TQS(IM)
         ps         , &! INTENT(IN        ) :: PQS(IM)
         eastar     , &! INTENT(OUT  ) :: QSS(IM)
         len          )! INTENT(IN        ) :: IM
   !     inter2 replaces interc

    CALL inter2(cuprt     , &! INTENT(INOUT) :: ppc   (len)
         lsprt     , &! INTENT(INOUT) :: ppl   (len)
         snow      , &! INTENT(INOUT) :: snoww (nsib,2)  
         capac     , &! INTENT(INOUT) :: capac (nsib,2)
         wwwtem    , &! INTENT(INOUT) :: WWW   (len,3)
         satcap    , &! INTENT(IN        ) :: satcap(len,2)
         cww       , &! INTENT(IN        ) :: cw
         tc        , &! INTENT(INOUT) :: tc(len)
         tg        , &! INTENT(INOUT) :: tg(len)
         clai      , &! INTENT(IN        ) :: clai
         zlt       , &! INTENT(IN        ) :: zlt(len)
         chil      , &! INTENT(IN        ) :: chil(len)
         roff      , &! INTENT(INOUT) :: roff(len)
         snomel    , &! INTENT(IN        ) :: snomel
         zdepth    , &! INTENT(IN        ) :: ZDEPTH(nsib,3)
         ts        , &! INTENT(IN        ) :: tm(len)
         tice      , &! INTENT(IN        ) :: tf
         asnow     , &! INTENT(IN        ) :: asnow
         csoil     , &! INTENT(IN        ) :: csoil(len)
         satco     , &! INTENT(IN        ) :: satco(len)
         dt        , &! INTENT(IN        ) :: dtt
         vcover    , &! INTENT(IN        ) :: vcover(len)
         roffo     , &! INTENT(OUT  ) :: roffo (len)
         zmelt2    , &! INTENT(OUT  ) :: zmelt(len)
         len       , &! INTENT(IN   ) :: len
         nsib      , &! INTENT(IN   ) :: nsib
         exo         )! INTENT(OUT  ) :: exo(len)
         
    gby100 = 0.01_r8  * grav
    !gbyhl =grav/(hltm*delsig*100.0_r8 )

    DO i = 1, len
       gbyhl =grav/hltm*delsig(i)*100.0_r8 
       gbycp =grav/(cp*delsig(i)*100.0_r8 *sigki(i))
       rmi   (i)=cu(i)*ustar(i)
       rhi   (i)=ct(i)*ustar(i)
       ah    (i)=gbycp/ps(i)
       al    (i)=gbyhl/ps(i)
       fac     =grav/(100.0_r8 *psb(i)*dt)
       !
       !     total mass of water and total sensible heat lost from the veggies.
       !
       IF(sfcpbl == 2)THEN       
          !sensf    (i) =(hc   (i) + hg(i)*(1.0_r8-AREAS(i)) + hs(i)*AREAS(i))*(1.0_r8/dt)!fss            (ncount)
          !latef    (i) =(ec   (i) + eg(i)*(1.0_r8-AREAS(i)) + es(i)*AREAS(i))*(1.0_r8/dt)!fws            (ncount)*hltm
          !!dtmdt=((ah(i)*sensf(i))/(dt*ah(i)*ros(i)*cp  *rhi(i)))/bps(i)
          !dqmdt=((al(i)*latef(i))/(dt*al(i)*ros(i)*hltm*rhi(i)))
          !dtm =dtmdt*dt
          !dqms=dqmdt*dt
          gmt(i,3)=(dth(i)/ BPS(I))/dt
          gmq(i,3)=dqm(i)/ dt
          tm   (i)=tg(i)/BPS(i)
          qm   (i)=qm(i)
       ELSE
          IF(atmpbl == 1)THEN
             !dtmdt   =(gmt(i,3) + ((hc(i)+hg(i))*(1.0_r8-AREAS(i)) + hs(i)*AREAS(i)) * fac   /(cp*bps(i)))/gmt(i,2)
             !dqmdt   =(gmq(i,3) + ((ec(i)+eg(i))*(1.0_r8-AREAS(i)) + es(i)*AREAS(i)) * hltm * fac)/ gmq(i,2)
             !dtm     =dtmdt * dt
             !dqms    =dqmdt * dt
             gmt(i,3)= gmt(i,3) + (dth(i) / BPS(I))/dt
             gmq(i,3)= gmq(i,3) + dqm(i)/ dt
             tm   (i)=tm(i)+(dth(i)/BPS(i))
             qm   (i)=qm(i)+ dqm(i)
          ELSE
             !sensf    (i) =(hc   (i) + hg(i)*(1.0_r8-AREAS(i)) + hs(i)*AREAS(i))*(1.0_r8/dt)!fss            (ncount)
             !latef    (i) =(ec   (i) + eg(i)*(1.0_r8-AREAS(i)) + es(i)*AREAS(i))*(1.0_r8/dt)!fws            (ncount)*hltm
             !dtmdt=((ah(i)*sensf(i))/(dt*ah(i)*ros(i)*cp  *rhi(i)))/bps(i)
             !dqmdt=((al(i)*latef(i))/(dt*al(i)*ros(i)*hltm*rhi(i)))
             !dtm     =dtmdt * dt
             !dqms    =dqmdt * dt
             !gmt(i,3)=dtmdt
             !gmq(i,3)=dqmdt
             gmt(i,3)=(dth(i) / BPS(I))/dt
             gmq(i,3)=dqm(i)/ dt
             !tm   (i)=tgeff(i)+dtm
             !qm   (i)=sha  (i)+dqms
             tm   (i)=tg(i)/BPS(i)
             qm   (i)=qm(i)
          END IF
       END IF
       !
       !     solve implicit system for winds
       !
       ! psb(i) = psur(i) * ( si(k) - si(k+1) )
       !
       !DRAG(I,1) = ROS(I) * CU(I) * USTAR(I)
       !DRAG(I,2) = ROS(I) * CUo(I) * USTARo(I)
       drag2(i)  =ROS(i)*cu(i)*ustar(i)
       !
       ! P=rho*G*Z ===> DP=rho*G*DZ
       !
       ! D                D
       !---- = rho * g * ----
       ! DZ               DP
       !                                D
       ! aaa = cu * ustar *  rho * g * ----
       !                                DP
       !
       !                                                   g
       ! aaa (i)  = rhoair(i)*cu(i)*ustar(i) * -------------------------------
       !                                       100*psur(i) * ( si(k) - si(k+1) )
       !
       aaa (i)  =drag2  (i)*gby100/psb(i)
       IF(sfcpbl == 2)THEN
          umom(i) =ROS(i)*(um(i)*sinclt(i))*rmi(i)
          vmom(i) =ROS(i)*(vm(i)*sinclt(i))*rmi(i)
          gb100   =grav/(delsig(i)*100.0_r8 )
          am   (i)=gb100/ps(i)
          gmu(i,3)=(am(i)*umom(i))/(dt*am(i)*ROS(i)*rmi(i))
          gmu(i,4)=(am(i)*vmom(i))/(dt*am(i)*ROS(i)*rmi(i))
       ELSE
          IF(atmpbl == 1)THEN
             gmu (i,2) =  gmu(i,2) + dt*aaa(i)
             gmu (i,3) = (gmu(i,3) - aaa(i) * um(i)*sinclt(i) ) / gmu(i,2)
             gmu (i,4) = (gmu(i,4) - aaa(i) * vm(i)*sinclt(i) ) / gmu(i,2)
          ELSE
            umom (i  ) = ROS(i)*(um(i)*sinclt(i))*rmi(i)
            vmom (i  ) = ROS(i)*(vm(i)*sinclt(i))*rmi(i)
            gb100      = grav/(delsig(i)*100.0_r8 )
            am   (i  ) = gb100/ps(i)
            gmu  (i,3) =(am(i)*umom(i))/(dt*am(i)*ROS(i)*rmi(i))
            gmu  (i,4) =(am(i)*vmom(i))/(dt*am(i)*ROS(i)*rmi(i))
          END IF
       END IF
!ssib2        umom(i)=ROS(i)*(um(i)*sinclt(i))*rmi(i)
!ssib2        vmom(i)=ROS(i)*(vm(i)*sinclt(i))*rmi(i)
!ssib2        gb100  =grav/(delsig*100.0_r8 )
!ssib2        am   (i)=gb100/ps(i)
!ssib2        gmu(i,3)=(gmu(i,3)-am(i)*umom(i))/(gmu(i,2)+dt*am(i)*ROS(i)*rmi(i))
!ssib2        gmu(i,4)=(gmu(i,4)-am(i)*vmom(i))/(gmu(i,2)+dt*am(i)*ROS(i)*rmi(i))
        qh(i)=0.622e0_r8*EXP(21.65605e0_r8 -5418.0e0_r8 /ta(i))/psb(i)
        bstar(i)=cu(i)*grav*(ct(i)*(ta(i)/BPS(I)-tm(i))/tm(i)*BPS(I))!+mapl_vireps*ct(i)*(qh(i)-qm(i)))

    END DO

    !     for perpetual conditions, do not update soil moisture
    IF(.NOT.fixday) THEN
       DO i = 1,len
          www(i,1) = wwwtem(i,1)
          www(i,2) = wwwtem(i,2)
          www(i,3) = wwwtem(i,3)
       ENDDO
    ENDIF

    DO i = 1,len
       xgpp(i) = assim(i) * 1.e6_r8
    ENDDO
    !     Calculate the surface flux of CO2 due to SiB (tracer T18)

    !     The net flux to the atmosphere is given by the release of CO2
    !     by soil respiration minus the uptake of CO2 by photosynthesis

    !     ASSIMN is the net assimilation of CO2 by the plants (from SiB2)
    !     respFactor*soilScale is the rate of release of CO2 by the soil

    !     soilScale is a diagnostic of the instantaneous rate of 
    !        soil respiration (derived by Jim Collatz, similar to TEM)
    !     respFactor is the annual total accumulation of carbon in the
    !        previous year at each grid cell (annual total ASSIMN) 
    !        divided by the annual total of soilScale at the same grid pt.

    !     Surface flux of CO2 used to be merely Assimn-Respg. With the
    !     prognostic CAS, the calculation becomes
    !
    !     co2flux =  (CO2A - CO2M)/ra
    !
    !     with a temperature correction thrown in. This calculation is 
    !     performed in phosib.


    IF(doSiBco2) THEN
       DO i = 1,len

          co2flx(i) = cflux(i)

       ENDDO
    ELSE
       DO i= 1,len
          co2flx(i) = 0.0_r8
       ENDDO
    ENDIF

    !    some quantities for diagnostic output
    DO i = 1,len
       rha(i) = ea(i)/eastar(i)
       zmelt(i) = zmelt1(i)+zmelt2(i)

       !           Calculate an overall leaf conductance, which is QP2(162) 
       !           and is used as a weighting function in QP2(162 and 163)
       ggl(i) = 1.0_r8 / (rst(i) * (rb(i) + rc(i)))
    ENDDO


    !itb...call neil suits' program cfrax2, which calculates 13C and 12C 
    !itb...fluxes and concentrations.

    !itb...will probably need to calculate some pre-subroutine variables...
    conductra = 1.0_r8/ra !units become meter/seconds

    IF(isotope)THEN
       CALL cfrax(dt          , &!INTENT(IN   ) :: dtt !the time  step in seconds
            d13cca      , &!INTENT(INOUT) :: d13Cca   (len) !del13C of canopy CO2 (per mil vs PDB)
            d13cm       , &!INTENT(IN   ) :: d13Cm    (len) !del13C of mixed layer CO2 (per mil vs PDB)
            d13cresp    , &!INTENT(IN   ) :: d13Cresp (len) !del13C of respiration CO2 (per mil vs PDB)   
            pco2ap_old  , &!INTENT(IN   ) :: pco2a    (len) !CO2 pressure in the canopy (pascals)      
            pco2m       , &!INTENT(IN   ) :: pco2m    (len) !CO2 pressure in the mixed layer (pascals)
            pco2s       , &!INTENT(IN   ) :: pco2s    (len) !CO2 pressure at the leaf surface (pascals)
            pco2c       , &!INTENT(IN   ) :: pco2c    (len) !CO2 pressure in the chloroplast (pascals)
            pco2i       , &!INTENT(IN   ) :: pco2i    (len) !CO2 pressure in the stoma (pascals)
            respg       , &!INTENT(IN   ) :: respg    (len) !rate of ground respiration (moles/m2/sec)         
            assimn      , &!INTENT(IN   ) :: assimn   (len) !rate of assimilation by plants (moles/m2/sec)
            conductra   , &!INTENT(INOUT) :: ga       (len) !1/resistance factor (ra) to mixing between the canopy 
            !and the overlying mixed layer. ga units:(m/sec)
       kiepsc3     , &!INTENT(OUT        ) :: KIECpsC3(len) !Kinetic Isotope Effect during C3 photosynthesis
            kiepsc4     , &!INTENT(OUT        ) :: KIECpsC4(len) !Kinetic Isotope Effect during C4 photosynthesis
            d13cassc3   , &!INTENT(OUT        ) :: d13CassimnC3  (len) !del13C of CO2 assimilated by C3 plants      
            d13cassc4   , &!INTENT(OUT        ) :: d13CassimnC4  (len) !del13C of CO2 assimilated by C3 plants           
            d13casstot  , &!INTENT(OUT        ) :: d13Cassimntot (len) !del13C of CO2 assimilated by all plants
            flux13c     , &!INTENT(OUT        ) :: Flux13C           (len) !turbulent flux of 13CO2 out of the canopy
            flux12c     , &!INTENT(OUT        ) :: Flux12C           (len) !turbulent flux of 12CO2 out of the canopy
            len         , &!INTENT(IN    ) :: len
            c4fract     , &!INTENT(IN        ) :: c4fract (len)         !fraction of vegetation that is C4: C. Still, pers.comm.
            z2          , &!INTENT(IN        ) :: ca_depth (len)         !depth of canopy
            ta          , &!INTENT(IN        ) :: t_air(len)          ! canopy temperature (K)
            flux_turb     )!INTENT(OUT        ) :: d13cflx_turb(len)   ! del 13C of turbulent flux of
       !                                        CO2 out of the canopy
       DO i=1,len
          discrim(i)=kiepsc3(i)*(1.0_r8-c4fract(i))+kiepsc4(i)*c4fract(i)
          discrim2(i)= discrim(i) * antemp(i)
       ENDDO
    ELSE
       d13cca(:)    = 0.0_r8
       d13cm(:)     = 0.0_r8
       d13cresp(:)  = 0.0_r8
       kiepsc3(:)   = 0.0_r8
       kiepsc4(:)   = 0.0_r8
       d13cassc3(:) = 0.0_r8
       d13cassc4(:) = 0.0_r8
       d13casstot(:)= 0.0_r8
       flux13c(:)   = 0.0_r8
       flux12c(:)   = 0.0_r8
       d13cflux(:)  = 0.0_r8
       c4fract(:)   = 0.0_r8
       discrim(:)   = 0.0_r8
       discrim2(:)  = 0.0_r8
       flux_turb(:) = 0.0_r8
    ENDIF
    RETURN
  END SUBROUTINE SIB

  SUBROUTINE Albedo_sib2( &
            ! Model information
            ncols     ,kmax      ,latco     , &
            ! Model Geometry
            cosz      ,zenith    , &
            ! Time info
            month2    ,month     , &
            ! Atmospheric fields
            wind      ,tsea      , &
            ! Microphysics
            taud      , &
            ! LW Radiation fields at last integer hour
            LwSfcDown , &
            ! Radiation field (Interpolated) at time = tod
            xVisBeam  , xVisDiff ,xNirBeam   , &
            xNirDiff  , &
            ! Surface Albedo
            avisb     ,avisd     ,anirb      , &
            anird             )
    IMPLICIT NONE 
   ! Model information
   INTEGER         , INTENT(IN   ) :: ncols
   INTEGER         , INTENT(IN   ) :: kmax
   INTEGER         , INTENT(IN   ) :: latco
   ! Model Geometry
   REAL   (KIND=r8), INTENT(IN   ) :: cosz  (ncols)
   REAL   (KIND=r8), INTENT(IN   ) :: zenith(ncols)
   ! Time info
   INTEGER         , INTENT(INOUT) :: month2(ncols)
   INTEGER         , INTENT(IN   ) :: month (ncols)
   ! Atmospheric fields
   REAL(KIND=r8) ,   INTENT(IN   ) :: wind  (ncols)!wind speed in m/s
   REAL(KIND=r8) ,   INTENT(IN   ) :: tsea  (ncols)
   ! Microphysics
   REAL(KIND=r8) ,   INTENT(IN   ) :: taud  (ncols,kMax)
   ! LW Radiation fields at last integer hour
   REAL(KIND=r8) ,   INTENT(IN   ) :: LwSfcDown(1:nCols)
   ! Radiation field (Interpolated) at time = tod
   REAL(KIND=r8) ,   INTENT(IN   ) :: xVisBeam (1:nCols)
   REAL(KIND=r8) ,   INTENT(IN   ) :: xVisDiff (1:nCols)
   REAL(KIND=r8) ,   INTENT(IN   ) :: xNirBeam (1:nCols)
   REAL(KIND=r8) ,   INTENT(IN   ) :: xNirDiff (1:nCols)
   ! Surface Albedo
   REAL(KIND=r8) ,   INTENT(OUT  ) :: avisb (ncols)
   REAL(KIND=r8) ,   INTENT(OUT  ) :: avisd (ncols)
   REAL(KIND=r8) ,   INTENT(OUT  ) :: anirb (ncols)
   REAL(KIND=r8) ,   INTENT(OUT  ) :: anird (ncols)


    REAL(KIND=r8) :: tc      (ncols) 
    REAL(KIND=r8) :: tg      (ncols)
    REAL(KIND=r8) :: snow    (ncols,2)
    REAL(KIND=r8) :: zlt     (ncols)  
    REAL(KIND=r8) :: z1      (ncols)
    REAL(KIND=r8) :: z2      (ncols)
    REAL(KIND=r8) :: tran    (ncols,2,2)
    REAL(KIND=r8) :: ref     (ncols,2,2)
    REAL(KIND=r8) :: chil    (ncols)
    REAL(KIND=r8) :: green   (ncols)
    REAL(KIND=r8) :: vcover  (ncols)
    REAL(KIND=r8) :: SoRef      (ncols,2)
    REAL(KIND=r8) :: capac     (ncols,2)
    REAL(KIND=r8) :: radfac    (ncols,2,2,2)
    REAL(KIND=r8) :: salb      (ncols,2,2)    
    REAL(KIND=r8) :: thermk    (ncols)            
    REAL(KIND=r8) :: tgeff4_sib(ncols)            
    REAL(KIND=r8) :: f
    REAL(KIND=r8) :: ocealb
    REAL(KIND=r8) :: IceOceanAlb (nCols,2,2)  

!     --------------------------- INPUT ---------------------------------------
!   
!    specify the parameters for albedo here:

    REAL(KIND=r8) ::         tau  (ncols)             !aerosol/cloud optical depth
    REAL(KIND=r8) ::         chl   (ncols)            !chlorophyll concentration in mg/m3
    
    INTEGER :: ncount
    INTEGER :: i,k
    INTEGER :: nsib

    ncount=0
    DO i=1,nCols
       IF(iMaskSiB2(i,latco) >=1_i8)THEN 
          ncount=ncount+1
          tc      (ncount)         = tcalc (ncount,latco) 
          tg      (ncount)         = tgrdc (ncount,latco)
          snow    (ncount,1)       = snowg (i,latco,1)
          snow    (ncount,2)       = snowg (i,latco,2)
          zlt     (ncount)         = timevar_lai(nCount,latco)
          z1      (ncount)         = BioTab(iMaskSiB2(i,latco))%z1! surface parameters
          z2      (ncount)         = BioTab(iMaskSiB2(i,latco))%z2! surface parameters
          tran    (ncount,1:2,1:2) = BioTab(iMaskSiB2(i,latco))%LTran(1:2,1:2)
          ref     (ncount,1:2,1:2) = BioTab(iMaskSiB2(i,latco))%LRef (1:2,1:2)
          chil    (ncount)         = BioTab(iMaskSiB2(i,latco))%ChiL
          green   (ncount)         = timevar_green(nCount,latco)
          vcover  (ncount)         = BioTab(iMaskSiB2(i,latco))%fVCover
          !IF(vcover(ncount) < zlt(ncount)/10.0_r8)THEN
          !   vcover(ncount) = vcover(ncount) * 10.0_r8
          !END IF
          SoRef     (ncount,1)  = BioTab(iMaskSiB2(i,latco))%SoRef(1)
          SoRef     (ncount,2)  = BioTab(iMaskSiB2(i,latco))%SoRef(2)
          capac            (ncount,1:2) = capacc(ncount,1:2,latco) 

          IF(iMaskSiB2(i,latco) == 13_i8)THEN
              zlt     (ncount)   = 0.0001_r8
              z1      (ncount)   = 0.0001_r8
              z2      (ncount)   = 0.1_r8
              vcover  (ncount)   = 0.0001_r8
              chil    (ncount)   = 0.01_r8
              green   (ncount)   = 0.1_r8         
          END IF
          radfac    (ncount,1:2,1:2,1:2)   = radfac_gbl(ncount,1:2,1:2,1:2,latco)
          salb      (ncount,1:2,1:2)       = AlbGblSiB2  (ncount,1:2,1:2,latco)
          thermk    (ncount)               = thermk_gbl (ncount,latco) 
          tgeff4_sib(ncount)               = tgeff4_sib_gbl (ncount,latco) 
       END IF
    END DO
    
    nsib=ncount
    
    CALL rada2(snow(1:nsib,1:2),zlt(1:nsib),z1(1:nsib),z2(1:nsib), &
            asnow,tg(1:nsib),cosz(1:nsib),tice,ref(1:nsib,1:2,1:2),&
            tran(1:nsib,1:2,1:2),chil(1:nsib) ,&
            green(1:nsib),vcover(1:nsib),soref(1:nsib,1:2),radfac(1:nsib,1:2,1:2,1:2),&
            salb(1:nsib,1:2,1:2),thermk(1:nsib),tgeff4_sib(1:nsib),&
            tc(1:nsib),nsib)
        
    ncount=0
    DO i=1,nCols
       IF(iMaskSiB2(i,latco) >=1_i8)THEN 
          ncount=ncount+1
          radfac_gbl     (ncount,1:2,1:2,1:2,latco)   =  radfac    (ncount,1:2,1:2,1:2)
          AlbGblSiB2     (ncount,1:2,1:2,latco)       =  salb      (ncount,1:2,1:2)
          thermk_gbl     (ncount,latco)               =  thermk    (ncount)            
          tgeff4_sib_gbl (ncount,latco)               =  tgeff4_sib(ncount)      
       END IF
    END DO
!     --------------------------- INPUT ---------------------------------------
!   
!    specify the parameters for albedo here:
    IceOceanAlb=0.0_r8
    chl = 0.10_r8              !chlorophyll concentration in mg/m3
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !
    tau=0.0_r8
    DO k=1,kMax
       DO i=1,ncols
          tau(i)=tau(i)+ taud(i,k) ! tau = SUM(taud(i,1:kMax))              !aerosol/cloud optical depth
       END DO
    END DO   
    DO i=1,ncols
        tau(i)=MIN(MAX(tau(i),0.0_r8),25.0_r8)      !aerosol/cloud optical depth
    END DO
  
    
    
    IF(TRIM(SLABOCEAN) == 'SLAB')THEN
       ncount=0
       DO i=1,ncols
          IF(iMaskSiB2(i,latco) >= 1_i8) THEN
             ncount=ncount+1
             avisb(i)=salb(ncount,1,1)
             avisd(i)=salb(ncount,1,2)
             anirb(i)=salb(ncount,2,1)
             anird(i)=salb(ncount,2,2)
          ELSE IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
             f=MAX(zenith(i),0.0e0_r8 )
             avisb(i)=GetOceanAlb(i,tau(i),f,wind(i),chl(i),0.2_r8,0.7_r8) !   s (shortwave)  = 0.2-0.7 micro-meters
             avisd(i)=oceald
             anirb(i)=GetOceanAlb(i,tau(i),f,wind(i),chl(i),0.7_r8,5.0_r8) !   l (longwave)   = 0.7-5.0 micro-meters
             anird(i)=oceald
          ELSE
              IF (TRIM(ICEMODEL)=='SSIB')THEN
                f=MAX(zenith(i),0.0e0_r8 )
                IceOceanAlb(i,:,:)=GetIceOceanAlb(i,latco,month(i),xVisBeam(i),xVisDiff(i),&
                                       xNirBeam(i),xNirDiff(i),f,LwSfcDown(i))
                avisb(i)= IceOceanAlb(i,1,1)
                avisd(i)= IceOceanAlb(i,1,2)
                anirb(i)= IceOceanAlb(i,2,1)
                anird(i)= IceOceanAlb(i,2,2)
             ELSE IF (TRIM(ICEMODEL)=='COLA')THEN
                avisb(i)=icealv
                avisd(i)=icealv
                anirb(i)=icealn
                anird(i)=icealn
             ELSE
                STOP "ICEMODEL ->OPTIONS"
             END IF
          END IF
       END DO

    ELSE IF(TRIM(SLABOCEAN) == 'COLA')THEN
    
       ncount=0
       DO i=1,ncols
          IF(iMaskSiB2(i,latco) >= 1_i8) THEN
             ncount=ncount+1
             avisb(i)=salb(ncount,1,1)
             avisd(i)=salb(ncount,1,2)
             anirb(i)=salb(ncount,2,1)
             anird(i)=salb(ncount,2,2)
          ELSE IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
             f=MAX(zenith(i),0.0e0_r8 )
             ocealb=0.12347e0_r8 +f*(0.34667e0_r8+f*(-1.7485e0_r8 + &
                  f*(2.04630e0_r8 -0.74839e0_r8 *f)))
             avisb(i)=ocealb
             avisd(i)=oceald
             anirb(i)=ocealb
             anird(i)=oceald
          ELSE
              IF (TRIM(ICEMODEL)=='SSIB')THEN
                f=MAX(zenith(i),0.0e0_r8 )
                IceOceanAlb(i,:,:)=GetIceOceanAlb(i,latco,month(i),xVisBeam(i),xVisDiff(i),&
                                       xNirBeam(i),xNirDiff(i),f,LwSfcDown(i))
                avisb(i)= IceOceanAlb(i,1,1)
                avisd(i)= IceOceanAlb(i,1,2)
                anirb(i)= IceOceanAlb(i,2,1)
                anird(i)= IceOceanAlb(i,2,2)
             ELSE IF (TRIM(ICEMODEL)=='COLA')THEN
                avisb(i)=icealv
                avisd(i)=icealv
                anirb(i)=icealn
                anird(i)=icealn
             ELSE
                STOP "ICEMODEL ->OPTIONS"
             END IF
          END IF
       END DO
    ELSE
       WRITE(0,*)"ERRO SLABOCEAN",TRIM(SLABOCEAN)
       STOP 
    END IF
  END SUBROUTINE Albedo_sib2

  !
  !=================SUBROUTINE RADA2====================================== 
  !
  SUBROUTINE RADA2(snoww,zlt,z1,z2                 &
       ,asnow,tg,sunang,tf,ref,tran,chil      &
       ,green,vcover,soref,radfac,salb,thermk &
       ,tgeff4,tc,nlen)

    IMPLICIT NONE                                            

    !=======================================================================
    !
    !     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT
    !     AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY
    !
    !-----------------------------------------------------------------------


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       SALB(2,2)      SURFACE ALBEDOS 
    !       TGEFF4         EFFECTIVE SURFACE RADIATIVE TEMPERATURE (K) 
    !       RADFAC(2,2,2)  RADIATION ABSORPTION FACTORS 
    !       THERMK         CANOPY GAP FRACTION FOR TIR RADIATION 
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       ALBEDO(2,2,2)  COMPONENT REFLECTANCES 
    !
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !     arg list declarations
    !
    INTEGER      , INTENT(IN   ) :: nlen
    REAL(KIND=r8)         , INTENT(IN   ) :: snoww (nlen,2)
    REAL(KIND=r8), INTENT(IN   ) :: ref   (nlen,2,2) !  leaf reflectance -- 4 components
    !    (:,1,1)  shortwave, green plants
    !    (:,1,2)  longwave, green plants
    !    (:,2,1)  shortwave, brown plants
    !    (:,2,2)  longwave, brown plants
    REAL(KIND=r8), INTENT(IN   ) :: tran  (nlen,2,2) !  leaf transmittance -- 4 components
    !    (:,1,1)  shortwave, green plants
    !    (:,1,2)  longwave, green plants
    !    (:,2,1)  shortwave, brown plants
    !    (:,2,2)  longwave, brown plants
    REAL(KIND=r8), INTENT(IN   ) :: soref (nlen,2)   ! soil reflectance (shortwave and longwave)
    REAL(KIND=r8), INTENT(OUT  ) :: salb  (nlen,2,2)
    REAL(KIND=r8), INTENT(OUT  ) :: radfac(nlen,2,2,2)
    REAL(KIND=r8), INTENT(IN   ) :: asnow
    REAL(KIND=r8), INTENT(IN   ) :: zlt   (nlen)     ! leaf area index
    REAL(KIND=r8), INTENT(IN   ) :: z1    (nlen)     ! height of bottom of canopy (m)
    REAL(KIND=r8), INTENT(IN   ) :: z2    (nlen)     ! height of top of canopy (m)
    REAL(KIND=r8), INTENT(IN   ) :: tg    (nlen)     ! ground temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: sunang(nlen)     !  
    REAL(KIND=r8), INTENT(IN   ) :: chil  (nlen)     ! leaf angle distribution factor
    REAL(KIND=r8), INTENT(IN   ) :: green (nlen)     ! green leaf fraction
    REAL(KIND=r8), INTENT(IN   ) :: vcover(nlen)     ! vegetation cover fraction
    REAL(KIND=r8), INTENT(OUT  ) :: tgeff4(nlen)
    REAL(KIND=r8), INTENT(IN   ) :: tc    (nlen)     ! canopy temperature (K)
    REAL(KIND=r8), INTENT(OUT  ) :: thermk(nlen)
    REAL(KIND=r8), INTENT(IN   ) :: tf

    !     local variables

    REAL(KIND=r8) :: TRANC1(2)
    REAL(KIND=r8) :: TRANC2(2)
    REAL(KIND=r8) :: TRANC3(2)
    REAL(KIND=r8) :: satcap(nlen,2)
    REAL(KIND=r8) :: f     (nlen)
    REAL(KIND=r8) :: fmelt (nlen)
    REAL(KIND=r8) :: zmew  (nlen)
    REAL(KIND=r8) :: albedo(nlen,2,2,2)
    REAL(KIND=r8) :: canex (nlen)
    REAL(KIND=r8) :: areas (nlen)
    !REAL(KIND=r8) :: deltg(nlen)
    REAL(KIND=r8) :: facs
    REAL(KIND=r8) :: scov
    REAL(KIND=r8) :: scat
    REAL(KIND=r8) :: chiv
    REAL(KIND=r8) :: aa
    REAL(KIND=r8) :: bb
    REAL(KIND=r8) :: fac2
    REAL(KIND=r8) :: fac1
    REAL(KIND=r8) :: zkat
    REAL(KIND=r8) :: tg4
    REAL(KIND=r8) :: tc4
    REAL(KIND=r8) :: tgs
    REAL(KIND=r8) :: hh6
    REAL(KIND=r8) :: hh5
    REAL(KIND=r8) :: hh3
    REAL(KIND=r8) :: hh2
    REAL(KIND=r8) :: zmk
    REAL(KIND=r8) :: hh10
    REAL(KIND=r8) :: hh9
    REAL(KIND=r8) :: hh8
    REAL(KIND=r8) :: hh7
    REAL(KIND=r8) :: den
    REAL(KIND=r8) :: zp
    REAL(KIND=r8) :: f1
    REAL(KIND=r8) :: ge
    REAL(KIND=r8) :: ek
    REAL(KIND=r8) :: epsi
    REAL(KIND=r8) :: power2
    REAL(KIND=r8) :: power1
    REAL(KIND=r8) :: zat
    REAL(KIND=r8) :: psi
    REAL(KIND=r8) :: hh4
    REAL(KIND=r8) :: hh1
    REAL(KIND=r8) :: fe
    REAL(KIND=r8) :: de
    REAL(KIND=r8) :: bot
    REAL(KIND=r8) :: ce
    REAL(KIND=r8) :: be
    REAL(KIND=r8) :: betao
    REAL(KIND=r8) :: upscat
    REAL(KIND=r8) :: acss
    REAL(KIND=r8) :: extkb
    REAL(KIND=r8) :: proj
    REAL(KIND=r8) :: tran2
    REAL(KIND=r8) :: tran1
    REAL(KIND=r8) :: reff1
    REAL(KIND=r8) :: reff2
    INTEGER       :: iwave
    INTEGER       :: i
    INTEGER       :: irad
    !
    !----------------------------------------------------------------------
    !
    !
    !     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO
    !         SNOW REFLECTANCE   = 0.80_r8, 0.40_r8_r8 . MULTIPLY BY 0.6 IF MELTING
    !         SNOW TRANSMITTANCE = 0.20_r8, 0.54_r8
    !
    !
    !-----------------------------------------------------------------------
    !
    !
    DO i = 1,nlen                      ! loop over gridpoint
       !        this portion is snow1 inlined
       CANEX (i)   = 1.0_r8 - ( SNOWw(i,1)*5.0_r8-Z1(i))/(Z2(i)-Z1(i))
       CANEX (i)   = MAX( 0.1E0_r8, CANEX(i) )
       CANEX (i)   = MIN( 1.0E0_r8, CANEX(i) )
       AREAS (i)   = MIN( 1.0E0_r8, ASNOW*SNOWw(i,2))
       SATCAP(i,1) = ZLT(i)*0.0001_r8 * CANEX(i)
       !
       ! Collatz-Bounoua change satcap(2) to 0.0002
       !
       !bl      SATCAP(i,2) = 0.002_r8
       !SATCAP(i,2) = 0.0002_r8             ! lahouari
       SATCAP(i,2) =  ZLT(i)*0.0001_r8 * CANEX(i)            ! lahouari

       !    end old snow1
       F(i) = MAX(0.01746E0_r8,SUNANG(i))
       FACS  = ( TG(i)-TF ) * 0.04_r8
       FACS  = MAX( 0.0E0_r8 , FACS)
       FACS  = MIN( 0.4E0_r8, FACS)
       FMELT(i) = 1.0_r8 - FACS
      !  deltg(i)=tf-tg(i)
      ! fmelt(i)=1.0_r8
      ! IF (ABS(deltg(i)) < 0.5_r8 .AND. deltg(i) > 0.0_r8) THEN
      !    fmelt(i)=0.6_r8
      ! END IF

    END DO
    !
    !-----------------------------------------------------------------------
    !
    DO IWAVE = 1, 2
       !
       !----------------------------------------------------------------------
       DO i = 1,nlen                      ! loop over gridpoint
          SCOV =  MIN( 0.5E0_r8, SNOWw(i,1)/SATCAP(i,1) )
          !IF (tc(i) <= tf) THEN
          !   SCOV= MIN( 0.5_r8  , capac(i,1)/satcap(i,1))
          !END IF
          REFF1 = ( 1.0_r8 - SCOV ) * REF(i,IWAVE,1) + SCOV * ( 1.2_r8 - &
               IWAVE * 0.4_r8 ) * FMELT(i)
          REFF2 = ( 1.0_r8 - SCOV ) * REF(i,IWAVE,2) + SCOV * ( 1.2_r8 - &
               IWAVE * 0.4_r8 ) * FMELT(i)
          TRAN1 = TRAN(i,IWAVE,1) * ( 1.0_r8 - SCOV) + SCOV * ( 1.0_r8 - &
               ( 1.2_r8 - IWAVE * 0.4_r8 ) * FMELT(i) ) * TRAN(i,IWAVE,1)
          TRAN2 = TRAN(i,IWAVE,2) * ( 1.0_r8 - SCOV ) &
               + SCOV * ( 1.0_r8 - ( 1.2_r8 - IWAVE * 0.4_r8 ) * FMELT(i) )&
               * 0.9_r8 * TRAN(i,IWAVE,2)
          !
          !-----------------------------------------------------------------------
          !
          !     CALCULATE AVERAGE SCATTERING COEFFICIENT, LEAF PROJECTION AND
          !     OTHER COEFFICIENTS FOR TWO-STREAM MODEL.
          !
          !      SCAT  (OMEGA)        : EQUATION (1,2) , SE-85
          !      PROJ  (G(MU))        : EQUATION (13)  , SE-85
          !      EXTKB (K, G(MU)/MU)  : EQUATION (1,2) , SE-85
          !      ZMEW  (INT(MU/G(MU)) : EQUATION (1,2) , SE-85
          !      ACSS  (A-S(MU))      : EQUATION (5)   , SE-85
          !      EXTK  (K, VARIOUS)   : EQUATION (13)  , SE-85
          !      UPSCAT(OMEGA*BETA)   : EQUATION (3)   , SE-85
          !      BETAO (BETA-0)       : EQUATION (4)   , SE-85 
          !
          !-----------------------------------------------------------------------
          !
          SCAT = GREEN(i)*( TRAN1 + REFF1 ) +( 1.0_r8 - GREEN(i) ) * &
               ( TRAN2 + REFF2)
          CHIV = CHIL(i)
          !
          IF ( ABS(CHIV) .LE. 0.01_r8 ) CHIV = 0.01_r8
          AA = 0.5_r8 - 0.633_r8 * CHIV - 0.33_r8 * CHIV * CHIV
          BB = 0.877_r8 * ( 1.0_r8 - 2.0_r8 * AA )
          !
          PROJ = AA + BB * F(i)
          EXTKB = ( AA + BB * F(i) ) / F(i)
          ZMEW(i) = 1.0_r8 / BB * ( 1.0_r8 - AA / BB  &
               * LOG ( ( AA + BB ) / AA ) )
          ACSS = SCAT / 2.0_r8 * PROJ / ( PROJ + F(i) * BB )

          ACSS = ACSS * ( 1.0_r8 - F(i) * AA           &
               / ( PROJ + F(i) * BB ) * LOG ( ( PROJ     &
               +   F(i) * BB + F(i) * AA ) / ( F(i) * AA ) ) )
          !
          UPSCAT = GREEN(i) * TRAN1 + ( 1.0_r8-GREEN(i) ) * TRAN2
          UPSCAT = 0.5_r8 * ( SCAT + ( SCAT - 2.0_r8 * UPSCAT ) * &
               (( 1.0_r8 - CHIV ) / 2.0_r8 ) ** 2 )
          BETAO = ( 1.0_r8 + ZMEW(i) * EXTKB )   &
               / ( SCAT * ZMEW(i) * EXTKB ) * ACSS
          !
          !-----------------------------------------------------------------------
          !
          !     Intermediate variables identified in appendix of SE-85.
          !
          !      BE          (B)     : APPENDIX      , SE-85
          !      CE          (C)     : APPENDIX      , SE-85
          !      BOT         (SIGMA) : APPENDIX      , SE-85
          !      HH1         (H1)    : APPENDIX      , SE-85
          !      HH2         (H2)    : APPENDIX      , SE-85
          !      HH3         (H3)    : APPENDIX      , SE-85
          !      HH4         (H4)    : APPENDIX      , SE-85
          !      HH5         (H5)    : APPENDIX      , SE-85
          !      HH6         (H6)    : APPENDIX      , SE-85
          !      HH7         (H7)    : APPENDIX      , SE-85
          !      HH8         (H8)    : APPENDIX      , SE-85
          !      HH9         (H9)    : APPENDIX      , SE-85
          !      HH10        (H10)   : APPENDIX      , SE-85
          !      PSI         (H)     : APPENDIX      , SE-85
          !      ZAT         (L-T)   : APPENDIX      , SE-85
          !      EPSI        (S1)    : APPENDIX      , SE-85
          !      EK          (S2)    : APPENDIX      , SE-85
          !-----------------------------------------------------------------------
          !
          BE = 1.0_r8 - SCAT + UPSCAT
          CE = UPSCAT
          BOT = ( ZMEW(i) * EXTKB ) ** 2 + ( CE**2 - BE**2 )
          IF ( ABS(BOT) .LE. 1.0E-10_r8 ) THEN
             SCAT = SCAT* 0.98_r8
             BE   = 1.0_r8 - SCAT + UPSCAT
             BOT  = ( ZMEW(i) * EXTKB ) ** 2 + ( CE**2 - BE**2 )
          END IF
          DE = SCAT * ZMEW(i) * EXTKB * BETAO
          FE = SCAT * ZMEW(i) * EXTKB * ( 1.0_r8 - BETAO )
          HH1 = -DE * BE + ZMEW(i) * DE * EXTKB - CE * FE
          HH4 = -BE * FE - ZMEW(i) * FE * EXTKB - CE * DE
          !
          PSI = SQRT(BE**2 - CE**2)/ZMEW(i)
          !
          ZAT = ZLT(i)/VCOVER(i)*CANEX(i)
          !
          POWER1 = MIN( PSI*ZAT, 50.0E0_r8 )
          POWER2 = MIN( EXTKB*ZAT, 50.E0_r8 )
          EPSI = EXP( - POWER1 )
          EK   = EXP ( - POWER2 )
          !
          ALBEDO(i,2,IWAVE,1) = SOREF(i,IWAVE)*(1.0_r8-AREAS(i)) &
               + ( 1.2_r8 - IWAVE*0.4_r8 )*FMELT(i) * AREAS(i)
          ALBEDO(i,2,IWAVE,2) = SOREF(i,IWAVE)*(1.0_r8-AREAS(i))   &
               + ( 1.2_r8-IWAVE*0.4_r8 )*FMELT(i) * AREAS(i)
          GE = ALBEDO(i,2,IWAVE,1)/ALBEDO(i,2,IWAVE,2)
          !
          !-----------------------------------------------------------------------
          !     CALCULATION OF DIFFUSE ALBEDOS
          !
          !      ALBEDO(1,IWAVE,2) ( I-UP ) : APPENDIX , SE-85
          !-----------------------------------------------------------------------
          !
          F1 = BE - CE / ALBEDO(i,2,IWAVE,2)
          ZP = ZMEW(i) * PSI
          !
          DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI - &
               ( BE - ZP ) * ( F1 + ZP ) * EPSI
          HH7 = CE * ( F1 - ZP ) / EPSI / DEN
          HH8 = -CE * ( F1 + ZP ) * EPSI / DEN
          F1 = BE - CE * ALBEDO(i,2,IWAVE,2)
          DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI
          !
          HH9 = ( F1 + ZP ) / EPSI / DEN
          HH10 = - ( F1 - ZP ) * EPSI / DEN
          TRANC2(IWAVE) = HH9 * EPSI + HH10 / EPSI
          !
          ALBEDO(i,1,IWAVE,2) =  HH7 + HH8
          !
          !-----------------------------------------------------------------------
          !     CALCULATION OF DIRECT ALBEDOS AND CANOPY TRANSMITTANCES.
          !
          !      ALBEDO(1,IWAVE,1) ( I-UP )   : EQUATION(11)   , SE-85
          !      TRANC(IWAVE)      ( I-DOWN ) : EQUATION(10)   , SE-85
          !
          !-----------------------------------------------------------------------
          !
          F1 = BE - CE / ALBEDO(i,2,IWAVE,2)
          ZMK = ZMEW(i) * EXTKB
          !
          DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI - &
               ( BE - ZP ) * ( F1 + ZP ) * EPSI
          HH2 = ( DE - HH1 / BOT * ( BE + ZMK ) ) &
               * ( F1 - ZP ) / EPSI - &
               ( BE - ZP ) * ( DE - CE*GE - HH1 / BOT &
               * ( F1 + ZMK ) ) * EK
          HH2 = HH2 / DEN
          HH3 = ( BE + ZP ) * (DE - CE*GE -  &
               HH1 / BOT * ( F1 + ZMK ))* EK - &
               ( DE - HH1 / BOT * ( BE + ZMK ) ) *  &
               ( F1 + ZP ) * EPSI
          HH3 = HH3 / DEN
          F1 = BE - CE * ALBEDO(i,2,IWAVE,2)
          DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI
          HH5 = - HH4 / BOT * ( F1 + ZP ) / EPSI - &
               ( FE + CE*GE*ALBEDO(i,2,IWAVE,2) + &
               HH4 / BOT*( ZMK-F1 ) ) * EK
          HH5 = HH5 / DEN
          HH6 =   HH4 / BOT * ( F1 - ZP ) * EPSI + &
               ( FE + CE*GE*ALBEDO(i,2,IWAVE,2) + &
               HH4 / BOT*( ZMK-F1 ) ) * EK
          HH6 = HH6 / DEN
          TRANC1(IWAVE) = EK
          TRANC3(IWAVE) = HH4 / BOT * EK + HH5 * EPSI + HH6 / EPSI
          !
          ALBEDO(i,1,IWAVE,1) = HH1 / BOT + HH2 + HH3
          !
          !----------------------------------------------------------------------
          !
          !
          !----------------------------------------------------------------------
          !     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES
          !     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND
          !
          !      RADFAC      (F(IL,IMU,IV)) : EQUATION (19,20) , SE-86
          !----------------------------------------------------------------------
          !
          RADFAC(i,2,IWAVE,1) = ( 1.0_r8-VCOVER(i) ) &
               * ( 1.0_r8-ALBEDO(i,2,IWAVE,1) ) + VCOVER(i) &
               * ( TRANC1(IWAVE) * ( 1.0_r8-ALBEDO(i,2,IWAVE,1) ) &
               + TRANC3(IWAVE) * ( 1.0_r8-ALBEDO(i,2,IWAVE,2) ) )
          !
          RADFAC(i,2,IWAVE,2) = ( 1.0_r8-VCOVER(i) )      &
               * ( 1.0_r8-ALBEDO(i,2,IWAVE,2) ) + VCOVER(i) &
               *  TRANC2(IWAVE) * ( 1.0_r8-ALBEDO(i,2,IWAVE,2) ) 
          !
          RADFAC(i,1,IWAVE,1) = VCOVER(i) &
               * ( ( 1.0_r8-ALBEDO(i,1,IWAVE,1) )&
               - TRANC1(IWAVE) * ( 1.0_r8-ALBEDO(i,2,IWAVE,1) )&
               - TRANC3(IWAVE) * ( 1.0_r8-ALBEDO(i,2,IWAVE,2) ) )
          !
          RADFAC(i,1,IWAVE,2) = VCOVER(i) &
               * ( ( 1.0_r8-ALBEDO(i,1,IWAVE,2) )&
               - TRANC2(IWAVE) * ( 1.0_r8-ALBEDO(i,2,IWAVE,2) ) )
       ENDDO
       !
       !----------------------------------------------------------------------
       !     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB ) WITH WEIGHTING
       !     FOR COVER FRACTIONS.
       !----------------------------------------------------------------------
       !
       DO  IRAD = 1, 2
          DO i = 1,nlen        !  loop over gridpoint
             SALB(i,IWAVE,IRAD) = ( 1.0_r8-VCOVER(i) )   &
                  * ALBEDO(i,2,IWAVE,IRAD) +     &
                  VCOVER(i) * ALBEDO(i,1,IWAVE,IRAD)
          END DO
       END DO
       !
       !----------------------------------------------------------------------
       !
    END DO
    !
    !----------------------------------------------------------------------
    !
    !     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND
    !
    !----------------------------------------------------------------------
    !
    DO i = 1,nlen                  !  loop over gridpoint
       TGS = MIN(TF,TG(i))*AREAS(i)  &
            + TG(i)*(1.0_r8-AREAS(i))
       TC4 = TC(i)**4
       TG4 = TGS**4
       !
       ZKAT = 1.0_r8/ZMEW(i) * ZLT(i) / VCOVER(i)
       ZKAT = MIN( 50.0E0_r8 , ZKAT )
       ZKAT = MAX( 1.0E-5_r8, ZKAT )
       THERMK(i) = EXP(-ZKAT)
       !
       FAC1 =  VCOVER(i) * ( 1.0_r8-THERMK(i) )
       FAC2 =  1.0_r8
       TGEFF4(i) =  FAC1 * TC4 &
            + (1.0_r8 - FAC1 ) * FAC2 * TG4
       TGEFF4(i) = SQRT ( SQRT ( TGEFF4(i)   ))
    END DO              
    RETURN
  END  SUBROUTINE RADA2
  !
  !==================SUBROUTINE RNLOAD====================================
  !
  SUBROUTINE rnload(nlen     , &! INTENT(IN   ) :: nlen
                 radvbc    , &! INTENT(IN   ) :: radvbc(nlen)
                 radvdc    , &! INTENT(IN   ) :: radndc(nlen)
                 radnbc    , &! INTENT(IN   ) :: radnbc(nlen)
                 radndc    , &! INTENT(IN   ) :: radvdc(nlen)
                 dlwbot    , &! INTENT(IN   ) :: dlwbot(nlen)
                 VCOVER    , &! INTENT(IN   ) :: VCOVER(nlen)
                 thermk    , &! INTENT(IN   ) :: thermk(nlen)
                 radfac    , &! INTENT(IN   ) :: radfac(nlen,2,2,2)
                 radn      , &! INTENT(OUT  ) :: radn  (nlen,2,2)
                 radc3       )! INTENT(OUT  ) :: radc3 (nlen,2)

    IMPLICIT NONE
    !
    !=======================================================================
    !
    !    calculation of absorption of radiation by surface.  Note that
    !       output from this calculation (radc3) only accounts for the 
    !       absorption of incident longwave and shortwave fluxes.  The
    !       total net radiation calculation is performed in subroutine
    !       netrad.
    !
    !=======================================================================
    !

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       RADN(2,3)      INCIDENT RADIATION FLUXES (W M-2)
    !       RADC3(2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2) 
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER, INTENT(IN   ) :: nlen
    REAL(KIND=r8)   , INTENT(OUT  ) :: radc3 (nlen,2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: radn  (nlen,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: radfac(nlen,2,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: radvbc(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: radvdc(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: radnbc(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: radndc(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: dlwbot(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: VCOVER(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: thermk(nlen)

    INTEGER i, iveg, iwave, irad

    !-----------------------------------------------------------------------
    !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
    !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
    !     SOURCE FOR TRANSPIRATION.
    !
    !      RADN        (F(IW,IMU,O)) : EQUATION (19-22) , SE-86
    !      RADC3       (FC,FGS)      : EQUATION (21,22) , SE-86
    !-----------------------------------------------------------------------

    DO i = 1,nlen
       RADc3(i,1) = 0.0_r8
       RADc3(i,2) = 0.0_r8
       radn(i,1,1) = radvbc(i)
       radn(i,1,2) = radvdc(i)
       radn(i,2,1) = radnbc(i)
       radn(i,2,2) = radndc(i)
    ENDDO

    DO iveg=1,2
       DO iwave=1,2
          DO irad=1,2
             DO i = 1,nlen
                radc3(i,iveg) = radc3(i,iveg) + &
                     radfac(i,iveg,iwave,irad) * &
                     radn(i,iwave,irad)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !     absorbed downwelling radiation only

    DO i = 1,nlen
       RADc3(i,1) = RADc3(i,1) + dlwbot(i) * &
            VCOVER(i) * (1.0_r8- THERMK(i))
       RADc3(i,2) = RADc3(i,2) + dlwbot(i) * &
            (1.0_r8-VCOVER(i) * (1.0_r8-THERMK(i)))
    ENDDO

    RETURN
  END SUBROUTINE rnload
  !
  !===================SUBROUTINE BALAN====================================
  !
  SUBROUTINE BALAN (IPLACE     , &! INTENT(IN   ) :: iplace
                  tau        , &! INTENT(IN   ) :: tau
                  zdepth     , &! INTENT(IN   ) :: zdepth(nlen,3)
                  www        , &! INTENT(IN   ) :: www   (nlen,3)
                  capac      , &! INTENT(IN   ) :: capac (nlen,2)
                  ppc        , &! INTENT(IN   ) :: ppc   (nlen)
                  ppl        , &! INTENT(IN   ) :: ppl   (nlen)
                  roff       , &! INTENT(INOUT) :: roff  (nlen)
                  etmass     , &! INTENT(INOUT) :: etmass(nlen)
                  totwb      , &! INTENT(INOUT) :: totwb (nlen)
                  radt       , &! INTENT(IN   ) :: radt  (nlen,2)
                  chf        , &! INTENT(IN   ) :: chf   (nlen)
                  shf        , &! INTENT(IN   ) :: shf   (nlen)
                  dtt        , &! INTENT(IN   ) :: dtt
                  ect        , &! INTENT(IN   ) :: ect   (nlen)
                  eci        , &! INTENT(IN   ) :: eci   (nlen)
                  egs        , &! INTENT(IN   ) :: egs   (nlen)
                  egi        , &! INTENT(IN   ) :: egi   (nlen)
                  hc         , &! INTENT(IN   ) :: hc    (nlen)
                  hg         , &! INTENT(IN   ) :: hg    (nlen)
                  heaten     , &! INTENT(IN   ) :: heaten(nlen)
                  hflux      , &! INTENT(IN   ) :: hflux (nlen)
                  snoww      , &! INTENT(IN   ) :: snoww (nlen,2)
                  thm        , &! INTENT(IN   ) :: thm   (nlen)
                  tc         , &! INTENT(IN   ) :: tc    (nlen)
                  tg         , &! INTENT(IN   ) :: tg    (nlen)
                  tgs        , &! INTENT(IN   ) :: tgs   (nlen)
                  td         , &! INTENT(IN   ) :: td    (nlen,nsoil)
                  ps         , &! INTENT(IN   ) :: ps    (nlen)
                  kapa       , &! INTENT(IN   ) :: kapa
                  !nsib       , &! INTENT(IN   ) :: nsib
                  nlen       , &! INTENT(IN   ) :: nlen
                  ioffset    , &! INTENT(IN   ) :: ioffset
                  nsoil        )! INTENT(IN   ) :: nsoil

    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: nlen
    !INTEGER, INTENT(IN   ) :: nsib

    INTEGER, INTENT(IN   ) :: iplace
    REAL(KIND=r8)   , INTENT(IN   ) :: tau
    REAL(KIND=r8)   , INTENT(IN   ) :: zdepth(nlen,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: www   (nlen,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: capac (nlen,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: ppc   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: ppl   (nlen)
    REAL(KIND=r8)   , INTENT(INOUT) :: roff  (nlen)
    REAL(KIND=r8)   , INTENT(INOUT) :: etmass(nlen)
    REAL(KIND=r8)   , INTENT(INOUT) :: totwb (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: radt  (nlen,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: chf   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: shf   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: dtt
    REAL(KIND=r8)   , INTENT(IN   ) :: ect   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: eci   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: egs   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: egi   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: hc    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: hg    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: heaten(nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: hflux (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: snoww (nlen,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: thm   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: tc    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: tgs   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: td    (nlen, nsoil)
    REAL(KIND=r8)   , INTENT(IN   ) :: ps    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: kapa
    !   INTEGER, INTENT(IN   ) :: nsib
    !    INTEGER, INTENT(IN   ) :: nlen
    INTEGER, INTENT(IN   ) :: ioffset
    INTEGER, INTENT(IN   ) :: nsoil

    !     local variables

    INTEGER  :: i
    INTEGER  :: j
    INTEGER  :: igp
    INTEGER  :: nerror
    INTEGER  :: indxerr(nlen)
    INTEGER  :: n
    REAL(KIND=r8)     :: endwb  (nlen)
    REAL(KIND=r8)     :: errorw (nlen)
    REAL(KIND=r8)     :: pcmeter(nlen)
    REAL(KIND=r8)     :: plmeter(nlen)
    REAL(KIND=r8)     :: emeter (nlen)
    REAL(KIND=r8)     :: cbal   (nlen)
    REAL(KIND=r8)     :: gbal   (nlen)
    REAL(KIND=r8)     :: errore (nlen)
    REAL(KIND=r8)     :: zlhs   (nlen)
    REAL(KIND=r8)     :: zrhs   (nlen)
    REAL(KIND=r8)     :: tm
    !
    !=======================================================================
    !
    !     ENERGY AND WATER BALANCE CHECK.
    !
    !-----------------------------------------------------------------------
    !
    IF( IPLACE .EQ. 1 ) THEN
       ! 
       DO i = 1,nlen
          ETMASS(i) = 0.0_r8
          ROFF(i)   = 0.0_r8
          !
          TOTWB(i) = WWW(i,1) * ZDEPTH(i,1)     &
               + WWW(i,2) * ZDEPTH(i,2)  &
               + WWW(i,3) * ZDEPTH(i,3)   &
               + CAPAC(i,1) + CAPAC(i,2) + snoww(i,1) + snoww(i,2) 
       ENDDO
       !
    ELSE
       !   
       nerror = 0     
       DO i = 1,nlen
          ENDWB(i) = WWW(i,1) * ZDEPTH(i,1)     &
               + WWW(i,2) * ZDEPTH(i,2)         &
               + WWW(i,3) * ZDEPTH(i,3)           &
               + CAPAC(i,1) + CAPAC(i,2) + snoww(i,1) + snoww(i,2)    &
               - (PPL(i)+PPC(i))*0.001_r8*dtt &
               + ETMASS(i)*0.001_r8 + ROFF(i)
          ERRORW(i)= TOTWB(i) - ENDWB(i)
          pcmeter(i) = ppc(i) * 0.001_r8*dtt
          plmeter(i) = ppl(i) * 0.001_r8*dtt
          EMETER(i)= ETMASS(i) * 0.001_r8
          !itb...trying a different error check, 1.e-6 is in the
          !itb...noise in 32-bit arithmetic, works fine in 64-bit
          !         if(abs(errorw(i)).gt.1.e-6) then
          IF(ABS(errorw(i)).GT.1.e-5_r8*totwb(i))THEN
             nerror = nerror + 1
             indxerr(nerror) = i
          ENDIF
       ENDDO
       !
       DO j = 1,nerror
          i = indxerr(j)
          igp = i+ioffset
          WRITE(6,900) tau,IGP, TOTWB(i), ENDWB(i), ERRORW(i),  &
               WWW(i,1), WWW(i,2), WWW(i,3),               &
               CAPAC(i,1), CAPAC(i,2),snoww(i,1),snoww(i,2), &
               pcmeter(i),plmeter(i), EMETER(i), ROFF(i) 
       ENDDO
       !
       !-----------------------------------------------------------------------
       !     
       nerror = 0
       DO i = 1,nlen
          CBAL(i) = RADT(i,1) - CHF(i) -  &
               (ECT(i)+HC(i)+ECI(i) )/DTT
          GBAL(i) = RADT(i,2) - SHF(i) - (EGS(i)+HG(i)+EGI(i) &
               - HEATEN(i) )/DTT             
          ZLHS(i) = RADT(i,1)+RADT(i,2) - CHF(i) - SHF(i)
          ZRHS(i) = HFLUX(i) + (ECT(i) + ECI(i) + EGI(i) + EGS(i) &
               + HEATEN(i) ) /DTT
          !
          ERRORE(i)= ZLHS(i) - ZRHS(i)
          IF(ABS(errore(i)).GT.1.0_r8) THEN
             nerror = nerror + 1
             indxerr(nerror) = i
          ENDIF
       ENDDO
       !
       DO j = 1,nerror
          i = indxerr(j)
          tm = thm(i) * (0.001_r8*ps(i))**kapa
          igp = i+ioffset
          WRITE(6,910) tau,IGP, ZLHS(i), ZRHS(i), &
               RADT(i,1), RADT(i,2),CHF(i), SHF(i), &   
               HFLUX(i), ECT(i), ECI(i), EGI(i), EGS(i), &
               tm,tc(i),tg(i),tgs(i)
          WRITE(6,911)(td(i,n),n=nsoil,1,-1)
          WRITE(6,912) HC(i), HG(i), HEATEN(i), CBAL(i), GBAL(i) 
          WRITE(6,901)  TOTWB(i), ENDWB(i), ERRORW(i),     &
               WWW(i,1), WWW(i,2), WWW(i,3),        &
               CAPAC(i,1), CAPAC(i,2),snoww(i,1),snoww(i,2), &
               pcmeter(i),plmeter(i), EMETER(i), ROFF(i)
       ENDDO
    ENDIF

900 FORMAT(//,10X,'** WARNING: WATER BALANCE VIOLATION **  ',//,   &  
         & /,1X,'TAU ', F10.2,' AT SIB POINT (I) = ',I5,                  &
         &/,1X,'BEGIN, END, DIFF  ', 2(F10.7,1X),g13.5,                   &
         &/,1X,'WWW,1-3           ', 3(F10.8,1X),                         &
         &/,1X,'CAPAC1-2,snow 1-2 ', 4(g13.5,1X),                         &
         &/,1X,'PPc,PPl, ET, ROFF ', 4(g13.5,1X) )                          
910 FORMAT(//,10X,'** WARNING: ENERGY BALANCE VIOLATION **',//,& 
         & /,1X,'TAU ', F10.2,' AT SIB POINT (I) = ',I5,              &
         & /,1X,'RHS, LHS              ', 2G13.5,                     &
         & /,1X,'RN1, RN2, CHF, SHF, H ', 5G13.5,                     &
         & /,1X,'ECT, ECI, EGI, EGS    ', 4G13.5,                     &
         & /,1X,'TM,  TC,  TG, TGS     ', 4G13.5,                     &
         & /,1X,'HC        HG          ',  G12.5, 12X, G12.5,         &
         & /,1X,'HEATEN, C-BAL, G-BAL  ', 3G13.5 )
911 FORMAT(1x,'TD ', 5G13.5)                          
912 FORMAT(1X,'HC        HG          ',  G12.5, 12X, G12.5,  &
         & /,1X,'HEATEN, C-BAL, G-BAL  ', 3G13.5 )
901 FORMAT(10X,'WATER BALANCE'                      &
         &/,1X,'BEGIN, END, DIFF  ', 2(F10.7,1X),g13.5, &
         &/,1X,'WWW,1-3           ', 3(F10.8,1X),       &
         &/,1X,'CAPAC1-2,snow 1-2 ', 4(g13.5,1X),       &
         &/,1X,'PPc,PPl, ET, ROFF ', 4(g13.5,1X) )
    !
    RETURN
  END SUBROUTINE BALAN
  !
  !================SUBROUTINE BEGTEM=======================================
  !
  SUBROUTINE begtem(tc      , &
                  tg      , &
                  cpair   , &
                  hlat    , &
                  psur    , &
                  snomel  , &
                  zlt     , &
                  clai    , &
                  cw      , &
                  www     , &
                  poros   , &
                  pie     , &
                  psy     , &
                  phsat   , &
                  bee     , &
                  ccx     , &
                  cg      , &
                  phc     , &
                  tgs     , &
                  etc     , &
                  etgs    , &
                  getc    , &
                  getgs   , &
                  rstfac  , &
                  rsoil   , &
                  hr      , &
                  wc      , &
                  wg      , &
                  snoww   , &
                  capac   , &
                  areas   , &
                  satcap  , &
                  csoil   , &
                  tf      , &
                  g       , &
                  snofac  , &
                  len     , &
                  sandfrac, &
                  clayfrac, &
                  nlen    , &
                  forcerestore )

    IMPLICIT NONE

    !========================================================================
    !
    !     Calculation of flux potentials and constants prior to heat 
    !         flux calculations.  Corresponds to first half of TEMREC
    !         in 1D model.
    !
    !========================================================================


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       RSTFAC(2)      SOIL MOISTURE STRESS FACTOR 
    !       RSOIL          SOIL SURFACE RESISTANCE (S M-1)
    !       HR             SOIL SURFACE RELATIVE HUMIDITY
    !       WC             CANOPY WETNESS FRACTION
    !       WG             GROUND WETNESS FRACTION
    !       CCX            CANOPY HEAT CAPACITY (J M-2 K-1)
    !       CG             GROUND HEAT CAPACITY (J M-2 K-1)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: nlen
    REAL(KIND=r8)   , INTENT(IN   ) :: tc    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: cpair    
    REAL(KIND=r8)   , INTENT(IN   ) :: hlat
    REAL(KIND=r8)   , INTENT(IN   ) :: psur  (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: snomel
    REAL(KIND=r8)   , INTENT(IN   ) :: zlt   (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: clai
    REAL(KIND=r8)   , INTENT(IN   ) :: cw
    REAL(KIND=r8)   , INTENT(IN   ) :: www   (len,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: pie
    REAL(KIND=r8)   , INTENT(OUT  ) :: psy   (len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: phsat (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: bee   (nlen)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ccx   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cg    (nlen)
    REAL(KIND=r8)   , INTENT(IN   ) :: phc   (nlen)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tgs   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: etc   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: etgs  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: getc  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: getgs (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rstfac(len,4)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rsoil (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hr    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wc    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wg    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: snoww (nlen,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: capac (nlen,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: areas (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: satcap(len,2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: csoil (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tf
    REAL(KIND=r8)   , INTENT(IN   ) :: g
    REAL(KIND=r8)   , INTENT(IN   ) :: sandfrac(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: clayfrac(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: snofac
    !   INTEGER, INTENT(IN   ) :: len
    !   INTEGER, INTENT(IN   ) :: nlen
    LOGICAL, INTENT(IN   ) :: forcerestore

    !     local variables
    INTEGER :: i  
    REAL(KIND=r8)    :: e
    REAL(KIND=r8)    :: csol
    REAL(KIND=r8)    :: ge
    REAL(KIND=r8)    :: slamda
    REAL(KIND=r8)    :: shcap
    REAL(KIND=r8)    :: x
    REAL(KIND=r8)    :: tsnow
    REAL(KIND=r8)    :: rsnow
    REAL(KIND=r8)    :: fac
    REAL(KIND=r8)    :: psit
    REAL(KIND=r8)    :: argg 
    REAL(KIND=r8)    :: phroot(len)
    REAL(KIND=r8)    :: chisl(len)
    REAL(KIND=r8)    :: siltfrac(len)
    REAL(KIND=r8)    :: shcapf
    REAL(KIND=r8)    :: shcapu
    REAL(KIND=r8)    :: one
    REAL(KIND=r8)    :: theta
    REAL(KIND=r8)    :: powliq 
    REAL(KIND=r8)    :: zcondry
    REAL(KIND=r8)    :: d1x
    REAL(KIND=r8)    :: spwet(len)
    REAL(KIND=r8)    :: ksoil(len)

    !----------------------------------------------------------------------
    !     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNCTION OF TEMPERATURE
    !     GE(X) IS D E(X) / D ( TEMP )
    !----------------------------------------------------------------------

    E (X) = EXP( 21.18123_r8 - 5418.0_r8 / X ) / .622_r8
    GE(X) = EXP( 21.18123_r8 - 5418.0_r8 / X ) * 5418.0_r8/ (X*X) / .622_r8

    one = 1.0_r8

    !----------------------------------------------------------------------
!hltm   = 2.52e6_r8
    SNOFAC   = HLAT/ ( HLAT + SNOMEL * 1.E-3_r8 )
    !----------------------------------------------------------------------
    !     CALCULATION OF CANOPY AND GROUND HEAT CAPACITIES.
    !     N.B. THIS SPECIFICATION DOES NOT NECESSARILY CONSERVE ENERGY WHEN
    !     DEALING WITH VERY LARGE SNOWPACKS.
    !     HEAT CAPACITY OF THE SOIL, AS USED IN FORCE-RESTORE HEAT FLUX
    !     DESCRIPTION. DEPENDENCE OF CSOIL ON POROSITY AND WETNESS.
    !
    !      CG          (cg)    : EQUATION (?) , CS-81
    !----------------------------------------------------------------------
    IF(forcerestore) THEN
       IF(SchemeDifus == 1 ) THEN
          DO i = 1,len
             SLAMDA = ( 1.5_r8*(1.0_r8-POROS(i)) + 1.3_r8*WWW(i,1)*POROS(i) ) / &
                      ( 0.75_r8 + 0.65_r8*POROS(i) -  0.4_r8*WWW(i,1)*POROS(i) ) &
                  * 0.4186_r8
             SHCAP  = ( 0.5_r8*(1.0_r8-POROS(i)) + WWW(i,1) * POROS(i)) * 4.186E6_r8
             CSOIL(i)  = SQRT( SLAMDA * SHCAP * 86400.0_r8/PIE ) * 0.5_r8

             !950511 adjust for different heat capacity of snow
             !       CCX(i) = ZLT(i) * CLAI +
             !    &             (0.5_r8*SNOWw(i,1)+CAPAC(i,1))*CW
             CCX(i) = ZLT(i)*CLAI+ (SNOWw(i,1)+CAPAC(i,1))*CW
             !950511 adjust for different heat capacity of snow
             !       CG(i)  = CSOIL(i) + 
             !    *       MIN ( 0.025_r8*one, (0.5_r8 *SNOWw(i,2)+CAPAC(i,2))) * CW
             CG(i)  = CSOIL(i) +  &
                    MIN ( 0.05_r8*one, (SNOWw(i,2)+CAPAC(i,2))) * CW
          ENDDO
       ELSE IF(SchemeDifus == 2 ) THEN
          DO i = 1,len
             siltfrac(i) = 0.01_r8 * max((100.0_r8 - 100.0_r8*sandfrac(i) - 100.0_r8*clayfrac(i)),0.0_r8)
             powliq = poros(i) * (www(i,1)+www(i,2)+www(i,3))/3.0_r8
             zcondry = sandfrac(i) * 0.300_r8 +  &
                       siltfrac(i) * 0.265_r8 +  &
                       clayfrac(i) * 0.250_r8 ! +
             ksoil(i) = zcondry * ((0.56_r8*100.0_r8)**powliq)

             !
             !     diffusivity of the soil
             !            --          --
             !           |    86400.0   |
             !d1x   =SQRT|--------------|*0.5
             !           |  (pie*difsl  |
             !            --          --
             d1x   =SQRT(86400.0_r8 /(3.14159_r8*5.0e-7_r8))*0.5_r8
             !
       
             theta   = ((www(i,1)+www(i,2)+www(i,3))/3.0_r8)*POROS(i)
             chisl(i) = ( 9.8e-4_r8 + 1.2e-3_r8 *theta ) / ( 1.1_r8 - 0.4_r8 *theta )
             chisl(i) = chisl(i)*4.186e2_r8
             !csoil(i)=chisl(i)*d1x       
             csoil(i)=ksoil(i)
             ccx(i)=zlt(i)*clai+capac(i,1)*cww
             spwet(i)=MIN( 0.05_r8 , capac(i,2))*cww
             cg(i)=csoil(i)+spwet(i)
          END DO
       ELSE
          WRITE(nfprt,'(A)')'ERRO DIFFUSIVITY OF THE SOIL SchemeDifus'
       END IF
    ELSE
       IF(SchemeDifus == 1 ) THEN
          DO i = 1,len
             !
             !--------------------
             !   new calculation for ground heat capacity cg - no longer force-restore
             !      now for the top 1cm of snow and soil, with phase change in the soil
             !      incorporated into the heat capacity from +0.5C to -0.5C
             csol    = (2.128_r8*sandfrac(i)+2.385_r8*clayfrac(i)) / (sandfrac(i)+clayfrac(i))*1.e6_r8  ! J/(m3 K)
             SHCAPu  = ( 0.5_r8 * (1.0_r8 - POROS(i)) + WWW(i,1) * POROS(i)) * 4.186E6_r8
             shcapf  =   0.5_r8 * (1.0_r8 + poros(i) * (www(i,1) - 1.0_r8) ) * 4.186E6_r8
             SHCAPu  = 0.2_r8*SHCAPu + 0.8_r8*csol
             shcapf  = 0.2_r8*shcapf + 0.8_r8*csol

             !SHCAPu  = ( 0.5_r8*(1.0_r8-POROS(i)) + WWW(i,1)* POROS(i)) * 4.186E6_r8
             !shcapf =  0.5_r8 * (1.0_r8 + poros(i) * (www(i,1)-1.0_r8)) * 4.186E6_r8
             IF(tg(i).GE.tf+0.5_r8) THEN
                csoil(i) = shcapu * 0.02_r8
             ELSE IF(tg(i).LE.tf-0.5_r8) THEN
                csoil(i) = shcapf * 0.02_r8
             ELSE
                csoil(i) = (0.5_r8*(shcapu+shcapf) +  snomel*poros(i)*www(i,1) ) * 0.02_r8
             ENDIF

             CCX(i) = ZLT(i)*CLAI + (0.5_r8*SNOWw(i,1)+CAPAC(i,1))*CW
             CG(i)  = (1.0_r8-areas(i))*CSOIL(i) + cw * (capac(i,2) + 0.01_r8 * areas(i))
             !
          ENDDO
       ELSE IF(SchemeDifus == 2 ) THEN
          DO i = 1,len
             siltfrac(i) = 0.01_r8 * max((100.0_r8 - 100.0_r8*sandfrac(i) - 100.0_r8*clayfrac(i)),0.0_r8)
             powliq = poros(i) * ((www(i,1)+www(i,2)+www(i,3))/3.0_r8)
             zcondry = sandfrac(i) * 0.300_r8 +  &
                      siltfrac(i) * 0.265_r8 +  &
                       clayfrac(i) * 0.250_r8 ! +
             ksoil(i) = zcondry * ((0.56_r8*100.0_r8)**powliq)
             
            ! ksoil(i)= (2.128*sandfrac+2.385*clayfrac) / (sandfrac+clayfrac)*1.e6  ! J/(m3 K)

             SHCAPu  = ( 0.5_r8*(1.0_r8-POROS(i)) + ((www(i,1)+www(i,2)+www(i,3))/3.0_r8)* POROS(i)) * 4.186E6_r8
             shcapf =  0.5_r8 * (1.0_r8 + poros(i) * (((www(i,1)+www(i,2)+www(i,3))/3.0_r8)-1.0_r8)) * 4.186E6_r8

             !
             !     diffusivity of the soil
             !            --          --
             !           |    86400.0   |
             !d1x   =SQRT|--------------|*0.5
             !           |  (pie*difsl  |
             !            --          --
              d1x   =SQRT(86400.0_r8 /(3.14159_r8*5.0e-7_r8))*0.5_r8
             !
       
             theta    = ((www(i,1)+www(i,2)+www(i,3))/3.0_r8)*POROS(i)
             chisl(i) = ( 9.8e-4_r8 + 1.2e-3_r8 *theta ) / ( 1.1_r8 - 0.4_r8 *theta )
             chisl(i) = chisl(i)*4.186e2_r8
             IF(tg(i).GE.tf+0.5_r8) THEN
                !csoil(i) = chisl(i)*d1x   
                csoil(i) = SHCAPu * ksoil(i) * 0.002_r8   
             ELSE IF(tg(i).LE.tf-0.5_r8) THEN
                csoil(i) = shcapf * 0.02_r8
             ELSE
                csoil(i) = (0.5_r8*(shcapu+shcapf) +  snomel*poros(i)*((www(i,1)+www(i,2)+www(i,3))/3.0_r8) ) * 0.02_r8
             ENDIF


             CCX(i) = ZLT(i)*CLAI + (0.5_r8*SNOWw(i,1)+CAPAC(i,1))*CW
             CG(i)  = (1.0_r8-areas(i))*CSOIL(i) + cw * (capac(i,2) + 0.01_r8 * areas(i))
          END DO
       ELSE
          WRITE(nfprt,'(A)')'ERRO DIFFUSIVITY OF THE SOIL SchemeDifus'
       END IF
    ENDIF
    DO i = 1,len
       !       HLAT(i)  = ( 3150.19_r8 - 2.378_r8 * TM(i) ) * 1000._r8 !use constant passed in
       PSY(i)      = CPAIR / HLAT * PSUR(i) / .622_r8
       !
       !----------------------------------------------------------------------
       !      Calculation of ground surface temperature and wetness fractions
       !        
       !----------------------------------------------------------------------
       !
       TSNOW = MIN ( TF-0.01_r8, TG(i) )
       RSNOW = SNOWw(i,2) / (SNOWw(i,2)+CAPAC(i,2)+1.E-10_r8)
       !
       TGS(i) = TSNOW*AREAS(i) + TG(i)*(1.0_r8-AREAS(i))
       IF(tgs(i).LT.0.0_r8)tgs(i) = SQRT(tgs(i))
       !
       etc(i)  = esat(tc(i))/100.0_r8
       GETC(i) = desat(tc(i))/100.0_r8

       ETGS(i)  = esat (TGS(i))/100.0_r8
       GETGS(i) = desat(TGS(i))/100.0_r8


       !ETC(i)   = E(TC(i))
       !ETGS(i)  = E(TGS(i))
       !GETC(i)  = GE(TC(i))
       !GETGS(i) = GE(TGS(i))
       !
       WC(i) = MIN( one,( CAPAC(i,1) + SNOWw(i,1))/SATCAP(i,1) )
       WG(i) = MAX( 0.0_r8*one,  CAPAC(i,2)/SATCAP(i,2) )*0.25_r8

       !-----------------------------------------------------------------------
       !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
       !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
       !     SOURCE FOR TRANSPIRATION.
       !
       !      PHROOT      (PSI-R) : EQUATION (48) , SE-86
       !      RSTFAC(2)  F(PSI-L) : MODIFICATION OF EQUATION (12), SE-89
       !-----------------------------------------------------------------------
       !
       PHROOT(i) = PHSAT(i) * MAX( 0.02_r8*one, WWW(i,2) ) ** ( - BEE(i) )
       PHROOT(i) = MAX ( PHROOT(i), -2.E3_r8*one )
       RSTFAC(i,2) = 1.0_r8/( 1.0_r8 + EXP( 0.02_r8*( PHC(i)-PHROOT(i)) ))
       RSTFAC(i,2) = MAX( 0.0001_r8*one, RSTFAC(i,2) )
       RSTFAC(i,2) = MIN( one,     RSTFAC(i,2) )
       !
       !----------------------------------------------------------------------
       !
       !      RSOIL FUNCTION FROM FIT TO FIFE-87 DATA.  Soil surface layer
       !         relative humidity.
       !
       !      RSOIL      (RSOIL) : HEISER 1992 (PERSONAL COMMUNICATION)
       !      HR         (Fh)    : EQUATION (66) , SE-86
       !----------------------------------------------------------------------
       !
       FAC = MIN( WWW(i,1), one )
       FAC = MAX( FAC, 0.02_r8*one  )
       !
       ! Collatz-Bounoua change rsoil to  FIFE rsoil formulation from eq(19) SE-92
       !
       !cbl     RSOIL(i) =  MAX (0.1*one, 694. - FAC*1500.) + 23.6 
       rsoil(i) =  EXP(8.206_r8 - 4.255_r8 * fac)      ! lahouari 

       !
       PSIT = PHSAT(i) * FAC ** (- BEE(i) )
       ARGG = MAX(-10.0_r8*one,(PSIT*G/ (461.5_r8*TGS(i)) ))
       HR(i) = EXP(ARGG)
    ENDDO
    RETURN
  END SUBROUTINE begtem
  !
  !==================SUBROUTINE NETRAD====================================
  !
  SUBROUTINE NETRAD(radc3 , &! INTENT(IN   ) :: radc3 (len,2)
                  radt  , &! INTENT(OUT  ) :: radt  (len,3)
                  stefan, &! INTENT(IN   ) :: stefan
                  fac1  , &! INTENT(OUT  ) :: fac1  (len)
                  vcover, &! INTENT(IN   ) :: vcover(len)
                  thermk, &! INTENT(IN   ) :: thermk(len)
                  tc    , &! INTENT(IN   ) :: tc    (len)
                  tg    , &! INTENT(IN   ) :: tg    (len)
                  tf    , &! INTENT(IN   ) :: tf
                  dtc4  , &! INTENT(OUT  ) :: dtc4  (len)
                  dtg4  , &! INTENT(OUT  ) :: dtg4  (len)
                  dts4  , &! INTENT(OUT  ) :: dts4  (len)
                  closs , &! INTENT(OUT  ) :: closs (len)
                  gloss , &! INTENT(OUT  ) :: gloss (len)
                  sloss , &! INTENT(OUT  ) :: sloss (len)
                  tgeff , &! INTENT(OUT  ) :: tgeff (len)
                  areas , &! INTENT(IN   ) :: areas (len)
                  zlwup , &! INTENT(OUT  ) :: zlwup(len)
                  len     )! INTENT(IN   ) :: len
    IMPLICIT NONE
    !
    !=======================================================================
    !
    !
    !        CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT
    !        LOSSES FROM CANOPY AND GROUND
    !
    !
    !=======================================================================
    !
    !pl bands in sib: 0.2 to 0.7 micromets are VIS, then 
    !pl               0.7 to 4.0 is NIR, above 4.0 it is thermal
    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       RADt (2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2) 
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER, INTENT(IN   ) :: len
    REAL(KIND=r8)   , INTENT(IN   ) :: radc3 (len,2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: radt  (len,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: stefan
    REAL(KIND=r8)   , INTENT(OUT  ) :: fac1  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: vcover(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: thermk(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tc    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tf
    REAL(KIND=r8)   , INTENT(OUT  ) :: dtc4  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dtg4  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dts4  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: closs (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gloss (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sloss (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tgeff (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: areas (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: zlwup(len)
    !   INTEGER, INTENT(IN   ) :: len

    !     local variables

    INTEGER :: i
    REAL(KIND=r8)    :: tc4
    REAL(KIND=r8)    :: tg4
    REAL(KIND=r8)    :: ts4
    REAL(KIND=r8)    :: tsnow
    REAL(KIND=r8)    :: radtbar

    REAL(KIND=r8)  ,PARAMETER :: feedfac=1.0_r8!  feedback factor 

    DO I=1,len
       tsnow = MIN(tf, tg(i))
       TC4 = TC(i)**4
       TG4 = TG(i)**4
       ts4 = tsnow**4

       !pl effective ground cover for thermal radiation 
       FAC1(i) =  VCOVER(i) * ( 1.0_r8-THERMK(i) )

       !pl derivatives
       DTC4(i) = 4*STEFAN * TC(i)**3
       DTG4(i) = 4*STEFAN * TG(i)**3
       Dts4(i) = 4*STEFAN * tsnow**3

       CLOSS(i) =  2.0_r8 * FAC1(i) * STEFAN * TC4
       !pl canopy leaves thermal radiation loss
       CLOSS(i) =  CLOSS(i) - FAC1(i) * STEFAN * &
            ( (1.0_r8-areas(i))*TG4+areas(i)*ts4)
       !pl ground loss
       GLOSS(i) =  STEFAN * TG4 - FAC1(i) * STEFAN * TC4
       !pl snow loss
       SLOSS(i) =  STEFAN * Ts4 - FAC1(i) * STEFAN * TC4
       !pl canopy leaves net radiation
       RADT(i,1) = RADC3(i,1) - closs(i)
       !pl ground net radiation
       RADT(i,2) = RADC3(i,2) - gloss(i)
       !pl snow net radiation 
       RADT(i,3) = RADC3(i,2) - sloss(i)
       !pl bulk, weighted net radiation from combined ground and snow
       radtbar = areas(i)*radt(i,3) + (1.0_r8-areas(i))*radt(i,2)
       !pl this is the exchange meant to help out exchanges between
       !pl ground and snow

       radt(i,2) = radtbar + (1.0_r8+feedfac)*(radt(i,2)-radtbar)
       radt(i,3) = radtbar + (1.0_r8+feedfac)*(radt(i,3)-radtbar)
       !pl total thermal radiation up from surface
       zlwup(i) = fac1(i) * tc4 + &
            (1.0_r8-fac1(i)) * (areas(i)*ts4+(1.0_r8-areas(i))*tg4)
       !pl effective (combined) skin temperature from surface thermal radiation
       TGEFF(i) =  ZLWUP(i) ** 0.25_r8
    ENDDO

    RETURN
  END SUBROUTINE NETRAD
  !--------------------

  SUBROUTINE VNQSAT (iflag   , &! INTENT(IN   ) :: iflag
                   TQS     , &! INTENT(IN   ) :: TQS(IM)
                   PQS     , &! INTENT(IN   ) :: PQS(IM)
                   QSS     , &! INTENT(OUT  ) :: QSS(IM)
                   IM        )! INTENT(IN   ) :: IM

    !     Computes saturation mixing ratio or saturation vapour pressure
    !     as a function of temperature and pressure for water vapour.

    !     INPUT VARIABLES     TQS,PQS,iflag,im
    !     OUTPUT VARIABLES    QSS
    !     SUBROUTINES CALLED  (AMAX1,AMIN1)
    !     iflag = 1 for saturation mixing ratio, otherwise for saturation 
    !             vapor pressure

    !        Modifications:
    !                - removed routine VHQSAT and made the call to it in c3vint.F
    !                  a call to VNQSAT adding the iflag=1 argument.  changed
    !                  subroutine name from VQSAT to VNQSAT and all calls to
    !                  VQSAT (c1subs.F, c3subs.F, comp3.F, hstatc.F, ocean.F) to
    !                  VNQSAT with the iflag=1 argument.  kwitt  10/23/91

    !     converted to fortran 90 syntax - dd 6/17/93

    !--------------------

    !     argument declarations
    INTEGER, INTENT(IN   ) :: IM
    INTEGER, INTENT(IN   ) :: iflag
    REAL(KIND=r8)   , INTENT(IN   ) :: TQS(IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: PQS(IM)
    REAL(KIND=r8)   , INTENT(OUT  ) :: QSS(IM)
    !   INTEGER, INTENT(IN   ) :: IM

    !     local declarations
    INTEGER :: ic(IM)
    REAL(KIND=r8)    :: tq1(IM)
    REAL(KIND=r8)    :: es1(IM)
    REAL(KIND=r8)    :: es2(IM)
    REAL(KIND=r8)    :: epsinv

    REAL(KIND=r8), PARAMETER :: EST(1:139)=(/ &
         0.0031195_r8,   0.0036135_r8,   0.0041800_r8,   0.0048227_r8,   0.0055571_r8, &
         0.0063934_r8,   0.0073433_r8,   0.0084286_r8,   0.0096407_r8,   0.0110140_r8, &
         0.0125820_r8,   0.0143530_r8,   0.0163410_r8,   0.0185740_r8,   0.0210950_r8, &
         0.0239260_r8,   0.0270960_r8,   0.0306520_r8,   0.0346290_r8,   0.0390730_r8, &
         0.0440280_r8,   0.0495460_r8,   0.0556910_r8,   0.0625080_r8,   0.0700770_r8, &
         0.0787000_r8,   0.0881280_r8,   0.0984770_r8,   0.1098300_r8,   0.1223300_r8, &
         0.1360800_r8,   0.1512100_r8,   0.1678400_r8,   0.1861500_r8,   0.2062700_r8, &
         0.2283700_r8,   0.2526300_r8,   0.2792300_r8,   0.3083800_r8,   0.3403000_r8, &
         0.3752000_r8,   0.4133400_r8,   0.4549700_r8,   0.5003700_r8,   0.5498400_r8, &
         0.6036900_r8,   0.6622500_r8,   0.7258900_r8,   0.7949700_r8,   0.8699100_r8, &
         0.9511300_r8,   1.0391000_r8,   1.1343000_r8,   1.2372000_r8,   1.3484000_r8, &
         1.4684000_r8,   1.5979000_r8,   1.7375000_r8,   1.8879000_r8,   2.0499000_r8, &
         2.2241000_r8,   2.4113000_r8,   2.6126000_r8,   2.8286000_r8,   3.0604000_r8, &
         3.3091000_r8,   3.5755000_r8,   3.8608000_r8,   4.1663000_r8,   4.4930000_r8, &
         4.8423000_r8,   5.2155000_r8,   5.6140000_r8,   6.0394000_r8,   6.4930000_r8, &
         6.9767000_r8,   7.4919000_r8,   8.0406000_r8,   8.6246000_r8,   9.2457000_r8, &
         9.9061000_r8,  10.6080000_r8,  11.3530000_r8,  12.1440000_r8,  12.9830000_r8, &
         13.8730000_r8,  14.8160000_r8,  15.8150000_r8,  16.8720000_r8,  17.9920000_r8, &
         19.1760000_r8,  20.4280000_r8,  21.7500000_r8,  23.1480000_r8,  24.6230000_r8, &
         26.1800000_r8,  27.8220000_r8,  29.5530000_r8,  31.3780000_r8,  33.3000000_r8, &
         35.3240000_r8,  37.4540000_r8,  39.6960000_r8,  42.0530000_r8,  44.5310000_r8, &
         47.1340000_r8,  49.8690000_r8,  52.7410000_r8,  55.7540000_r8,  58.9160000_r8, &
         62.2320000_r8,  65.7080000_r8,  69.3510000_r8,  73.1680000_r8,  77.1640000_r8, &
         81.3480000_r8,  85.7250000_r8,  90.3050000_r8,  95.0940000_r8, 100.1000000_r8, &
         105.3300000_r8, 110.8000000_r8, 116.5000000_r8, 122.4600000_r8, 128.6800000_r8, &
         135.1700000_r8, 141.9300000_r8, 148.9900000_r8, 156.3400000_r8, 164.0000000_r8, &
         171.9900000_r8, 180.3000000_r8, 188.9500000_r8, 197.9600000_r8, 207.3300000_r8, &
         217.0800000_r8, 227.2200000_r8, 237.7600000_r8, 248.71_r8/)

    EPSINV = 1.0_r8/1.622_r8

    tq1(:) = MAX(1.00001E0_r8,(tqs(:)-198.99999_r8))
    tq1(:) = MIN(138.900001E0_r8,tq1(:))
    ic (:) = INT(tq1(:),kind=i4)

    es1(:) = est(ic(:))
    es2(:) = est(ic(:)+1)

    qss(:) = ic(:)
    qss(:) = es1(:) + (es2(:)-es1(:)) * (tq1(:)-qss(:))
    tq1(:) = pqs(:) * epsinv
    qss(:) = MIN(tq1(:),qss(:))

    IF (IFLAG .EQ. 1) qss(:) = 0.622_r8 * qss(:) / (pqs(:)-qss(:))

    RETURN
  END SUBROUTINE VNQSAT


  SUBROUTINE VMFCALZ(&
                   !PR0    , &! INTENT(IN   ) :: pr0 !not used
                   !RIBC   , &! INTENT(IN   ) :: ribc!not used
                   VKRMN  , &! INTENT(IN   ) :: vkrmn
                   DELTA  , &! INTENT(IN   ) :: delta
                   GRAV   , &! INTENT(IN   ) :: grav
                   SH     , &! INTENT(IN   ) :: SH    (len)
                   !PS     , &! INTENT(IN   ) :: PS    (len)            !not used
                   z0     , &! INTENT(IN   ) :: z0    (len)
                   SPDM   , &! INTENT(IN   ) :: SPDM  (len)
                   SHA    , &! INTENT(IN   ) :: SHA   (len)
                   ROS    , &! INTENT(IN   ) :: ROS   (len)
                   CU     , &! INTENT(OUT  ) :: CU    (len)      
                   ct     , &! INTENT(OUT  ) :: CT    (len)
                   THVGM  , &! INTENT(OUT  ) :: THVGM (len)
                   RIB    , &! INTENT(OUT  ) :: RIB   (len)
                   USTAR  , &! INTENT(OUT  ) :: USTAR (len)
                   VENTMF , &! INTENT(OUT  ) :: VENTMF(len)
                   THM    , &! INTENT(IN   ) :: THM   (len)
                   tha    , &! INTENT(IN   ) :: tha   (len) 
                   zzwind , &! INTENT(IN   ) :: zzwind (len)
                   zztemp , &! INTENT(IN   ) :: zztemp (len)
                   cuni   , &! INTENT(OUT  ) :: cuni  (len)
                   cun    , &! INTENT(OUT  ) :: cun    (len)
                   ctn    , &! INTENT(OUT  ) :: ctn    (len)
                   z1z0Urt, &! INTENT(OUT  ) :: z1z0Urt(len)
                   z1z0Trt, &! INTENT(OUT  ) :: z1z0Trt(len)
                   len      )! INTENT(IN   ) :: len

    IMPLICIT NONE

    !*****************************************************************************
    !     VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum 
    !     and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
    !     and Boville(1993), and by Beljaars and Holtslag (1991).              
    !  
    !     Rerences:
    !       Beljars and Holtslag (1991): Flux parameterization over land surfaces
    !              for atmospheric models. J. Appl. Meteo., 30, 327-341.
    !       Holtslag and Boville (1993): Local versus nonlocal boundary-layer 
    !              diffusion in a global climate model. J. of Climate, 6, 1825-
    !              1842.
    !       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
    !              atmosphere. Boundary-Layer Meteo., 17, 187-202.
    !       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
    !              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
    !              Layer parameterization, ECMWF, 59-79.
    !
    !     General formulation:
    !        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.) 
    !     Transfer coefficients for mommentum and heat fluxes are:
    !        CU = CUN*Fm, and
    !        CT = CTN*Fh
    !        where  CUN and CTN are nutral values of momentum and heat transfers,
    !           and Fm and Fh are stability functions derived from surface
    !           similarity relationships.     
    !*****************************************************************************

    INTEGER, INTENT(IN   ) :: len
    !REAL(KIND=r8)   , INTENT(IN   ) :: pr0 !not used
    !REAL(KIND=r8)   , INTENT(IN   ) :: ribc!not used
    REAL(KIND=r8)   , INTENT(IN   ) :: vkrmn
    REAL(KIND=r8)   , INTENT(IN   ) :: delta
    REAL(KIND=r8)   , INTENT(IN   ) :: grav
    REAL(KIND=r8)   , INTENT(IN   ) :: z0    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CU    (len)      
    REAL(KIND=r8)   , INTENT(IN   ) :: ROS   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: SH    (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: PS    (len)        !not used
    REAL(KIND=r8)   , INTENT(IN   ) :: THM   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: THVGM (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: RIB   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: SHA   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: SPDM  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: USTAR (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CT    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tha   (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: VENTMF(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cuni  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zzwind (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zztemp (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: z1z0Urt(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: z1z0Trt(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cun    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ctn    (len)

    !     local variables
    REAL(KIND=r8)    :: CTI   (len)
    REAL(KIND=r8)    :: CUI   (len)
    REAL(KIND=r8)    :: THGM  (len)  
    REAL(KIND=r8)    :: TEMV     (len)
    REAL(KIND=r8)    :: wgm      (len)
    REAL(KIND=r8)    :: bunstablM
    REAL(KIND=r8)    :: bunstablT
    REAL(KIND=r8)    :: cunstablM
    REAL(KIND=r8)    :: cunstablT
    REAL(KIND=r8)    :: bstabl
    REAL(KIND=r8)    :: cstabl
    REAL(KIND=r8)    :: zrib     (len)
    REAL(KIND=r8)    :: z1z0U    (len)
    REAL(KIND=r8)    :: z1z0T    (len)
    REAL(KIND=r8)    :: fmomn    (len)
    REAL(KIND=r8)    :: fheat    (len)
    REAL(KIND=r8)    :: ribtemp
    REAL(KIND=r8)    :: dm
    REAL(KIND=r8)    :: dh

    INTEGER :: i    
    !      
    !  condtants for surface flux functions, according to Holtslag and
    !      Boville (1993, J. Climate)

    bunstablM = 10.0_r8       ! constants for unstable function
    bunstablT = 15.0_r8
    cunstablM = 75.0_r8
    cunstablT = 75.0_r8
    bstabl =  8.0_r8           ! constants for stable function
    cstabl = 10.0_r8
    !
    DO I=1,len
       zrib(i) = zzwind(i) ** 2 / zztemp(i)
       WGM(i)  = SHA(i) - SH(i)
       !
       !        SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE
       !        WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.
       !
       THGM(i)  = THA(i)  - THM(i)
       THVGM(i) = THGM(i) + THA(i) * DELTA * WGM(i)
    END DO

    !   Ratio of reference height (zwind/ztemp) and roughness length:
    DO i = 1, len                !for all grid points
       z1z0U(i) = zzwind(i)/ z0(i)
       z1z0Urt(i) = SQRT( z1z0U(i) )
       z1z0U(i) = LOG( z1z0U(i) )
       z1z0T(i) = zzwind(i)/ z0(i)
       z1z0Trt(i) = SQRT( z1z0T(i) )
       z1z0T(i) = LOG( z1z0T(i) )
    ENDDO

    !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

    DO i = 1, len       
       cun(i) = VKRMN*VKRMN / (z1z0U(i)*z1z0U(i) )   !neutral Cm & Ct
       ctn(i) = VKRMN*VKRMN / (z1z0T(i)*z1z0T(i) )
       cuni(i) = z1z0u(i) / vkrmn
    ENDDO
    !
    !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.            
    !   RIB IS THE BULK RICHARDSON NUMBER, between reference height and surface.
    !
    DO I=1,len
       TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
       temv(i) = MAX(0.000001E0_r8,temv(i))
       RIB(I) = -THVGM(I) * GRAV * zrib(i) / TEMV(i) 
    END DO

    !   The stability functions for momentum and heat/moisture fluxes as
    !   derived from the surface-similarity theory by Luis (1079, 1982), and
    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag 
    !   (1991).

    DO I=1,len 
       IF(rib(i).GE.0.0_r8) THEN
          !
          !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
          !
          rib(i) = MIN( rib(i), 0.5E0_r8)
          fmomn(i) = (1.0_r8 + cstabl * rib(i) * (1.0_r8+ bstabl * rib(i)))
          fmomn(i) = 1.0_r8 / fmomn(i)
          fmomn(i) = MAX(0.0001E0_r8,fmomn(i))
          fheat(i) = fmomn(i)

       ELSE
          !
          !        THE UNSTABLE CASE.
          !            
          ribtemp = ABS(rib(i))
          ribtemp = SQRT( ribtemp )
          dm = 1.0_r8 + cunstablM * cun(i) * z1z0Urt(i) * ribtemp
          dh = 1.0_r8 + cunstablT * ctn(i) * z1z0Trt(i) * ribtemp
          fmomn(i) = 1.0_r8 - (bunstablM * rib(i) ) / dm
          fheat(i) = 1.0_r8 - (bunstablT * rib(i) ) / dh

       END IF
    END DO

    !   surface-air transfer coefficients for momentum CU, for heat and 
    !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

    DO i = 1, len
       CU(i) = CUN(i) * fmomn(i) 
       CT(i) = CTN(i) * fheat(i)
       CUI(i) = 1.0_r8 / CU(i)
       CTI(i) = 1.0_r8 / CT(i)
    END DO

    !   Ustar and ventlation mass flux: note that the ustar and ventlation 
    !   are calculated differently from the Deardoff's methods due to their
    !   differences in define the CU and CT.

    DO I=1,len 
       USTAR(i) = SPDM(i)*SPDM(i)*CU(i) 
       USTAR(i) = SQRT( USTAR(i) ) 
       VENTMF(i)= ROS(i)*CT(i)* SPDM(i)
    END DO
    !
    !   Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION  
    !   VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
    !   free convection conditions.
    !  

    RETURN
  END SUBROUTINE VMFCALZ


  SUBROUTINE VMFCAL(PR0     , &! INTENT(IN   ) :: pr0
                  RIBC    , &! INTENT(IN   ) :: ribc
                  VKRMN   , &! INTENT(IN   ) :: vkrmn
                  DELTA   , &! INTENT(IN   ) :: delta
                  GRAV    , &! INTENT(IN   ) :: grav
                  SH      , &! INTENT(IN   ) :: SH    (len)
                  !PS      , &! INTENT(IN   ) :: PS(len) ! not used
                  z0      , &! INTENT(IN   ) :: z0    (len)
                  SPDM    , &! INTENT(IN   ) :: SPDM  (len)
                  SHA     , &! INTENT(IN   ) :: SHA   (len)
                  ZB      , &! INTENT(IN   ) :: ZB    (len)
                  ROS     , &! INTENT(IN   ) :: ROS   (len)
                  CU      , &! INTENT(OUT  ) :: CU    (len)
                  THVGM   , &! INTENT(OUT  ) :: THVGM (len)
                  RIB     , &! INTENT(OUT  ) :: RIB   (len)
                  USTAR   , &! INTENT(OUT  ) :: USTAR (len)
                  VENTMF  , &! INTENT(OUT  ) :: VENTMF(len)
                  THM     , &! INTENT(IN   ) :: THM   (len)
                  tha     , &! INTENT(IN   ) :: tha   (len) 
                  cni     , &! INTENT(OUT  ) :: CNI   (len)  
                  cuni    , &! INTENT(OUT  ) :: CUNI  (len)
                  ctni    , &! INTENT(OUT  ) :: CTNI  (len)
                  ctni3   , &! INTENT(OUT  ) :: CTNI3 (len) 
                  cti     , &! INTENT(OUT  ) :: CTI   (len)
                  cui     , &! INTENT(OUT  ) :: CUI   (len)
                  len     , &! INTENT(IN   ) :: len
                  thgeff  , &! INTENT(IN   ) :: thgeff(len)
                  shgeff  , &! INTENT(IN   ) :: shgeff(len)
                  tke     , &! INTENT(IN   ) :: tke   (len)
                  ct      , &! INTENT(OUT  ) :: CT    (len)
                  dotkef    )! INTENT(IN   ) :: dotkef
    IMPLICIT NONE
    !
    !***  VENTILATION MASS FLUX, BASED ON DEARDORFF, MWR, 1972
    !
    INTEGER, INTENT(IN   ) :: len
    LOGICAL, INTENT(IN   ) :: dotkef
    REAL(KIND=r8)   , INTENT(IN   ) :: pr0
    REAL(KIND=r8)   , INTENT(IN   ) :: ribc
    REAL(KIND=r8)   , INTENT(IN   ) :: vkrmn
    REAL(KIND=r8)   , INTENT(IN   ) :: delta
    REAL(KIND=r8)   , INTENT(IN   ) :: grav
    REAL(KIND=r8)   , INTENT(IN   ) :: z0    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CU    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ZB    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ROS   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: SH    (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: PS    (len) ! not used
    REAL(KIND=r8)   , INTENT(IN   ) :: THM   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: THVGM (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: RIB   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CNI   (len)  
    REAL(KIND=r8)   , INTENT(OUT  ) :: CUNI  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: SHA   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: SPDM  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: USTAR (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CTNI  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CUI   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CTI   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CT    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tha   (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: VENTMF(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: CTNI3 (len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: thgeff(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: shgeff(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tke   (len)

    !     local variables

    REAL(KIND=r8)    :: ribmax
    REAL(KIND=r8)    :: vkrinv
    REAL(KIND=r8)    :: athird
    REAL(KIND=r8)    :: TEMV   (len)
    REAL(KIND=r8)    :: ZDRDRF (len)
    REAL(KIND=r8)    :: CHOKE  (len)
    REAL(KIND=r8)    :: wgm    (len)
    REAL(KIND=r8)    :: THGM   (len)
    REAL(KIND=r8)    :: sqrtke
    INTEGER i          
    !
    RIBMAX = 0.9_r8*RIBC
    VKRINV = 1.0_r8/VKRMN
    ATHIRD = 1.0_r8/3.0_r8
    !
    DO I=1,len
       WGM(i)  = SHA(i) - SH(i)
       !
       !        SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE
       !        WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.
       !
       THGM(i)  = THA(i)  - THM(i)
    END DO
    IF (dotkef) THEN
       DO i = 1,len
          thvgm(i) = (thgeff(i)-thm(i)) + thgeff(i) * delta *((shgeff(i)-sh(i)))
       ENDDO
    ELSE
       DO i = 1,len
          THVGM(i) = THGM(i) + THA(i) * DELTA * WGM(i)
       ENDDO
    ENDIF
    !
    !        CUNI AND CTN1 ARE INVERSES OF THE NEUTRAL TRANSFER COEFFICIENTS
    !        DEARDORFF EQS(33) AND (34).
    !        PR0 IS THE TURBULENT PRANDTL NUMBER AT NEUTRAL STABILITY.
    !
    !
    DO I=1,len
       CNI(i) = 0.025_r8*ZB(i)/z0(i)
       cni(i) = LOG(cni(i))
       CNI(i)  = CNI(i)  * VKRINV
       CUNI(i) = CNI(i)  + 8.4_r8
       CTNI(i) = CNI(i)  * PR0 + 7.3_r8
       CTNI3(i)= 0.3_r8     * CTNI(i)
       CNI(i)  = GRAV    * ZB(i)
       !
       !        SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
       !        RIB IS THE BULK RICHARDSON NUMBER, DEARDORFF EQ (25)
       !
       TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
       temv(i) = MAX(0.000001E0_r8,temv(i))
    END DO
    IF (dotkef) THEN
       DO I=1,len
          rib(i) = -thvgm(i)*grav*zb(i)/(thgeff(i)*tke(i)) 
       ENDDO
    ELSE
       DO I=1,len
          RIB(I) = -THVGM(I) * CNI(I) / TEMV(i) 
       ENDDO
    ENDIF
    !   
    !
    IF (dotkef) THEN
       DO I=1,len 
          IF(rib(i).GE.0.0_r8) THEN
             cu(i) = ((-0.0307_r8*rib(i)**2)/(61.5_r8+rib(i)**2))+0.044_r8
             ct(i) = ((-0.0152_r8*rib(i)**2)/(190.4_r8+rib(i)**2))+0.025_r8
          ELSE
             cu(i) = ((-0.016_r8*rib(i)**2)/(4.2e4_r8+rib(i)**2))+0.044_r8
             ct(i) = ((-0.0195_r8*rib(i)**2)/(2.1e4_r8+rib(i)**2))+0.025_r8
          ENDIF
       ENDDO
       DO I=1,len 
          sqrtke = SQRT(tke(i))
          USTAR(I) = SQRT(sqrtke*spdm(I)*CU(I))
          VENTMF(I)=ROS(I)*CT(I)*sqrtke
       ENDDO
    ELSE
       DO I=1,len 
          IF(rib(i).GE.0.0_r8) THEN
             !
             !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
             !
             CHOKE(i)=1.0_r8-MIN(RIB(I),RIBMAX)/RIBC
             CU(I)=CHOKE(i)/CUNI(I)
             CT(I)=CHOKE(i)/CTNI(I)
          ELSE
             !
             !        FIRST, THE UNSTABLE CASE. DEARDORFF EQS(28), (29), AND (30)
             !
             ZDRDRF(i) = LOG10(-RIB(I)) - 3.5_r8
             CUI(I)    = CUNI(I) - 25.0_r8 * EXP(0.26_r8 * ZDRDRF(i) &
                  - 0.03_r8 * ZDRDRF(i) * ZDRDRF(i))
             CTI(I)    = CUI(I) + CTNI(I) - CUNI(I)
             CU(I)     = 1.0_r8 / MAX(CUI(I),0.5_r8 * CUNI(I))
             CT(I)     = 1.0_r8 / MAX(CTI(I),CTNI3(I))
          END IF
       END DO
       DO I=1,len
          USTAR(i) =SPDM(i)*CU(i)
          VENTMF(i)=ROS(i)*CT(i)*USTAR(i)
       END DO
    ENDIF

    !
    !     CHECK THAT VENTFC EXCEEDS TOWNSENDS(1962) FREE CONVECTION VALUE, 
    !     DEARDORFF EQ(40B)
    !  
    IF (.NOT.dotkef) THEN
       DO I=1,len 
          IF( rib(i).LT.0.0_r8) THEN
             IF( cti(i).LT.ctni3(i) ) &
                  VENTMF(I)=MAX(VENTMF(I),ROS(I)*0.00186_r8*(THVGM(I))**ATHIRD) 
          ENDIF
       END DO
    ENDIF
    RETURN
    !
  END SUBROUTINE VMFCAL


  !
  !==================SUBROUTINE RBRD=======================================
  !
  SUBROUTINE RBRD(&
       tc,&
       !tm,&
       rbc,&
       zlt,&
       !tg,&
       z2,&
       u2,   &
       rd,&
       rb,&
       ta,&
       g,&
       rdc,&
       !ra,&
       tgs,&
       !wg,    &
       !rsoil,&
       !fg,&
       !hr,&
       !areas,&
       !respcp,&
       len )

    IMPLICIT NONE

    !========================================================================
    !
    !      CALCULATION OF RB AND RD AS FUNCTIONS OF U2 AND TEMPERATURES
    !
    !======================================================================== 


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       RB (GRB)       CANOPY TO CAS AERODYNAMIC RESISTANCE (S M-1)
    !       RD (GRD)       GROUND TO CAS AERODYNAMIC RESISTANCE (S M-1)
    !       TA (GTA)       CAS TEMPERATURE (K)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER, INTENT(IN   ) :: len
    REAL(KIND=r8)   , INTENT(IN   ) :: tc    (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: tm    (len)!isnot used
    REAL(KIND=r8)   , INTENT(IN   ) :: rbc   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zlt   (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: tg    (len)!isnot used
    REAL(KIND=r8)   , INTENT(IN   ) :: z2    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: u2    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rd    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rb    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ta    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rdc   (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: ra    (len)!isnot used
    REAL(KIND=r8)   , INTENT(IN   ) :: tgs   (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: wg    (len)!isnot used
    !REAL(KIND=r8)   , INTENT(IN   ) :: rsoil (len)!isnot used
    !REAL(KIND=r8)   , INTENT(IN   ) :: fg    (len)!isnot used
    !REAL(KIND=r8)   , INTENT(IN   ) :: hr    (len)!isnot used
    !REAL(KIND=r8)   , INTENT(IN   ) :: areas (len)!isnot used
    !REAL(KIND=r8)   , INTENT(IN   ) :: respcp(len)!isnot used
    REAL(KIND=r8)   , INTENT(IN   ) :: g

    REAL(KIND=r8)    :: temdif(len)
    INTEGER :: i
    REAL(KIND=r8)    :: fac
    REAL(KIND=r8)    :: fih
    !
    !-----------------------------------------------------------------------
    !      RB       (RB)       : EQUATION (A9), SE-86
    !-----------------------------------------------------------------------
    !
    DO i = 1,len                      !  loop over gridpoint
       TEMDIF(i) = MAX( 0.1E0_r8,  TC(i)-TA(i) )
       FAC = ZLT(i)/890.0_r8* (TEMDIF(i)*20.0_r8)**0.25_r8
       !RB(i)  = 1.0_r8/(SQRT(U2(i))/RBC(i)+FAC)
       rb (i)=1.0_r8  /(SQRT(u2(i))/rbc(i)+zlt(i)*0.004_r8 )

       !
       !-----------------------------------------------------------------------
       !      RD       (RD)       : EQUATION (A15), SE-86
       !-----------------------------------------------------------------------
       !
       TEMDIF(i) = MAX( 0.1E0_r8, TGs(i)-TA(i) )
       FIH = SQRT( 1.0_r8+9.0_r8*G*TEMDIF(i)*Z2(i)/(TGS(i)*U2(i)*U2(i)) )
       RD(i)  = RDC(i) / (U2(i) * FIH) 

    ENDDO
    !
    RETURN
  END SUBROUTINE RBRD



  SUBROUTINE respsib(len,  nsoil, wopt, zm, www, wsat, &
       tg, td, forcerestore, respfactor, respg, soilscale, &
       zmstscale, soilq10)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: len
    !INTEGER, INTENT(IN   ) :: nsib!is not used
    INTEGER, INTENT(IN   ) :: nsoil
    REAL(KIND=r8)   , INTENT(IN   ) :: wopt       (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zm         (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: www        (len,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsat       (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg         (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: td         (len,nsoil)
    REAL(KIND=r8)   , INTENT(IN   ) :: respfactor (len,nsoil+1)
    ! output      
    REAL(KIND=r8)   , INTENT(OUT  ) :: respg      (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: soilscale  (len,nsoil+1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: zmstscale  (len,2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: soilq10    (len,nsoil+1)
    LOGICAL, INTENT(IN   ) :: forcerestore

    ! local arrays

    REAL(KIND=r8)    :: b(len,2)
    REAL(KIND=r8)    :: woptzm
    INTEGER :: i,l

    !    Calculates the rate of CO2 efflux from soils, according to the "R-star"
    !    approach of Denning et al (1996), adapted for use with the Bonan
    !    6-layer soil thermodynamics module

    !    Changed soil Q10 value for soil respiration from 2.0 to 2.4
    !    following Raich and Schelsinger (1992, Tellus 44B, 81-89),
    !    Scott Denning, 9/14/95

    IF(.NOT.forcerestore) THEN

       DO i = 1,len

          ! Moisture effect on soil respiration, shallow and root zone

          woptzm = wopt(i)**zm(i)
          b(i,1) = (((100.0_r8*www(i,1))**zm(i)-woptzm)/ &
                   (woptzm - 100.0_r8**zm(i)))**2
          b(i,2) = (((100.0_r8*www(i,2))**zm(i)-woptzm)/ &
               (woptzm - 100.0_r8**zm(i)))**2
          b(i,1) = MIN(b(i,1),10.0E0_r8)
          b(i,2) = MIN(b(i,2),10.0E0_r8)
          zmstscale(i,1) = 0.8_r8*wsat(i)**b(i,1) + 0.2_r8
          zmstscale(i,2) = 0.8_r8*wsat(i)**b(i,2) + 0.2_r8

          ! Temperature effect is a simple Q10, with a reference T of 25 C

          ! Deepest soil layers do not respire (no carbon below root zone)
          DO L = 1, 2
             soilscale(i,L) = 0.0_r8
             soilq10(i,L) = 0.0_r8
          ENDDO

          ! Layers 3 through nsoil (root zone) use WWW(2) and TD(3:nSoil)
          DO l = 3, nsoil
             soilQ10(i,L) = EXP(0.087547_r8 * (td(i,L) - 298.15_r8))
             soilscale(i,L) = soilQ10(i,L) * zmstscale(i,2)
          ENDDO

          ! Surface soil uses TG and water layer 1
          soilQ10(i,nsoil+1) = EXP(0.087547_r8 * (tg(i) - 298.15_r8))
          soilscale(i,nsoil+1) = soilQ10(i,nsoil+1) * zmstscale(i,1)

          ! Dimensionalize soil resp flux to balance annual budget
          respg(i) = respfactor(i,1) * soilscale(i,1)

          DO L = 2, nSoil+1
             respg(i) = respg(i) + respfactor(i,L) * soilscale(i,L)
          ENDDO

       ENDDO

    ELSE  ! (FORCERESTORE CASE ... only two soil T levels available)

       DO i = 1,len

          ! Moisture effect from TEM
          woptzm = wopt(i)**zm(i)
          b(i,2) = (((100.0_r8*www(i,2))**zm(i)-woptzm)/ &
               (woptzm - 100.0_r8**zm(i)))**2
          b(i,2) = MIN(b(i,2),10.0E0_r8)
          zmstscale(i,1) = 0.8_r8*wsat(i)**b(i,2) +0.2_r8

          ! Temperature effect is Q10 =2.4 from ref T of 25 C
          soilQ10(i,1) = EXP(0.087547_r8 * (td(i,nsoil) - 298.15_r8))
          soilscale(i,1) = soilQ10(i,1) * zmstscale(i,1)

          ! Dimensionalize soil resp flux to balance annual budget
          respg(i) = respfactor(i,1) * soilscale(i,1)

       ENDDO


    ENDIF
    RETURN
  END SUBROUTINE respsib

  !
  !===================SUBROUTINE SORTIN===================================
  !
  SUBROUTINE SORTIN( EYY, PCO2Y, RANGE, GAMMAS, IC,len )
    !
    !=======================================================================
    !
    !     ARRANGES SUCCESSIVE PCO2/ERROR PAIRS IN ORDER OF INCREASING PCO2.
    !       ESTIMATES NEXT GUESS FOR PCO2 USING COMBINATION OF LINEAR AND
    !       QUADRATIC FITS.
    !
    !=======================================================================
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: ic
    REAL(KIND=r8)   , INTENT(INOUT) :: EYY     (len,6)
    REAL(KIND=r8)   , INTENT(INOUT) :: PCO2Y   (len,6)
    REAL(KIND=r8)   , INTENT(IN   ) :: RANGE   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: gammas  (len)

    !     work arrays

    REAL(KIND=r8)    :: eyyi1   (len)
    REAL(KIND=r8)    :: eyyi2   (len)
    REAL(KIND=r8)    :: eyyi3   (len)
    REAL(KIND=r8)    :: eyyis   (len)
    REAL(KIND=r8)    :: eyyisp  (len)
    REAL(KIND=r8)    :: pco2yis (len)
    REAL(KIND=r8)    :: pco2yisp(len)
    REAL(KIND=r8)    :: pco2b   (len)
    REAL(KIND=r8)    :: pco2yi1 (len)
    REAL(KIND=r8)    :: pco2yi2 (len)
    REAL(KIND=r8)    :: pco2yi3 (len)
    REAL(KIND=r8)    :: aterm
    REAL(KIND=r8)    :: bterm
    REAL(KIND=r8)    :: cterm
    REAL(KIND=r8)    :: pco2yq
    REAL(KIND=r8)    :: cc1
    REAL(KIND=r8)    :: cc2
    REAL(KIND=r8)    :: bc1
    REAL(KIND=r8)    :: bc2
    REAL(KIND=r8)    :: ac1
    REAL(KIND=r8)    :: ac2
    REAL(KIND=r8)    :: pco2yl
    REAL(KIND=r8)    :: a
    REAL(KIND=r8)    :: b
    REAL(KIND=r8)    :: pmin
    REAL(KIND=r8)    :: emin
    REAL(KIND=r8)    :: one
    INTEGER :: is(len)
    LOGICAL :: bitx(len)
    INTEGER :: i
    INTEGER :: ix
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: isp
    INTEGER :: n
    INTEGER :: l
    INTEGER :: j

    !
    one = 1.0_r8
    IF( IC .LT. 4 ) THEN
       DO i = 1,len
          PCO2Y(i,1) = GAMMAS(i) + 0.5_r8*RANGE(i)
          PCO2Y(i,2) = GAMMAS(i)  &
               + RANGE(i)*( 0.5_r8 - 0.3_r8*SIGN(one,EYY(i,1)) )
          PCO2Y(i,3) = PCO2Y(i,1)- (PCO2Y(i,1)-PCO2Y(i,2)) &
               /(EYY(i,1)-EYY(i,2)+1.E-10_r8)*EYY(i,1)
          !
          PMIN = MIN( PCO2Y(i,1), PCO2Y(i,2) )
          EMIN = MIN(   EYY(i,1),   EYY(i,2) )
          IF ( EMIN .GT. 0.0_r8 .AND. PCO2Y(i,3) .GT. PMIN ) &
               PCO2Y(i,3) = GAMMAS(i)
       ENDDO
    ELSE
       !
       N = IC - 1
       DO l = 1,len
          bitx(l) = ABS(eyy(l,n)).GT.0.1_r8
          IF(.NOT.bitx(l)) pco2y(l,ic) = pco2y(l,n)
       ENDDO
       DO l = 1,len
          IF(bitx(l)) THEN
             DO J = 2, N
                A = EYY(l,J)
                B = PCO2Y(l,J)
                DO I = J-1,1,-1
                   IF(EYY(l,I) .LE. A ) go to 100
                   EYY(l,I+1) = EYY(l,I)
                   PCO2Y(l,I+1) = PCO2Y(l,I)
                END DO
                i = 0
100             CONTINUE
                EYY(l,I+1) = A
                PCO2Y(l,I+1) = B
             END DO
          ENDIF
       ENDDO
       !
       !-----------------------------------------------------------------------
       !
       DO l = 1,len
          IF(bitx(l)) THEN
             PCO2B(l) = 0.0_r8
             IS(l)    = 1
          ENDIF
       ENDDO

       DO IX = 1, N
          DO l = 1,len
             IF(bitx(l)) THEN
                IF( EYY(l,IX) .LT. 0.0_r8 )  THEN
                   PCO2B(l) = PCO2Y(l,IX)
                   IS(l) = IX
                ENDIF
             ENDIF
          ENDDO
       END DO

       DO l = 1,len
          IF(bitx(l)) THEN
             I1 = IS(l)-1
             I1 = MAX(1, I1)
             I1 = MIN(N-2, I1)
             I2 = I1 + 1
             I3 = I1 + 2
             ISP   = IS(l) + 1
             ISP = MIN0( ISP, N )
             IS(l) = ISP - 1
             eyyisp(l) = eyy(l,isp)
             eyyis(l) = eyy(l,is(l))
             eyyi1(l) = eyy(l,i1)
             eyyi2(l) = eyy(l,i2)
             eyyi3(l) = eyy(l,i3)
             pco2yisp(l) = pco2y(l,isp)
             pco2yis(l) = pco2y(l,is(l))
             pco2yi1(l) = pco2y(l,i1)
             pco2yi2(l) = pco2y(l,i2)
             pco2yi3(l) = pco2y(l,i3)
          ENDIF
       ENDDO
       !
       DO l = 1,len
          IF(bitx(l)) THEN

             !itb...Neil Suits' patch to check for zero in the denominator...
             IF(EYYis(l) /= EYYisp(l))THEN
                PCO2YL=PCO2Yis(l)              &
                     - (PCO2Yis(l)-PCO2Yisp(l))  &
                     /(EYYis(l)-EYYisp(l))*EYYis(l)
             ELSE
                PCO2YL = PCO2Yis(l) * 1.01_r8
             ENDIF
             !
             !   METHOD USING A QUADRATIC FIT
             !
             AC1 = EYYi1(l)*EYYi1(l) - EYYi2(l)*EYYi2(l)
             AC2 = EYYi2(l)*EYYi2(l) - EYYi3(l)*EYYi3(l)
             BC1 = EYYi1(l) - EYYi2(l)
             BC2 = EYYi2(l) - EYYi3(l)
             CC1 = PCO2Yi1(l) - PCO2Yi2(l)
             CC2 = PCO2Yi2(l) - PCO2Yi3(l)

             !itb...Neil Suits' patch to prevent zero in denominator...
             IF(BC1*AC2-AC1*BC2 /= 0.0_r8 .AND. AC1 /= 0.0_r8)THEN
                BTERM = (CC1*AC2-CC2*AC1)/(BC1*AC2-AC1*BC2)
                ATERM = (CC1-BC1*BTERM)/AC1
                CTERM = PCO2Yi2(l)   &
                     -ATERM*EYYi2(l)*EYYi2(l)-BTERM*EYYi2(l)
                PCO2YQ= CTERM
                PCO2YQ= MAX( PCO2YQ, PCO2B(l) )
                PCO2Y(l,IC) = ( PCO2YL+PCO2YQ)/2.0_r8
             ELSE
                PCO2Y(l,IC) = PCO2Y(l,IC) * 1.01_r8
             ENDIF

          ENDIF
       ENDDO
       !
    ENDIF

    DO i = 1,len
       pco2y(i,ic) = MAX(pco2y(i,ic),0.01E0_r8)
    ENDDO
    !
    RETURN
  END SUBROUTINE SORTIN
  !
  !=====================SUBROUTINE CYCALC=================================
  !

  SUBROUTINE CYCALC( APARKK, VM, ATHETA, BTHETA, par,   &
       GAMMAS, RESPC, RRKK, OMSS, C3, C4, &
       PCO2I, ASSIMN, assim, len )
    IMPLICIT NONE
    !
    !=======================================================================
    !
    !     CALCULATION EQUIVALENT TO STEPS IN FIGURE 4 OF SE-92A
    !     C4 CALCULATION BASED ON CO-92.
    ! 
    !=======================================================================
    ! 

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
    !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
    !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       OMC            RUBISCO LIMITED ASSIMILATION (MOL M-2 S-1)
    !       OME            LIGHT LIMITED ASSIMILATION (MOL M-2 S-1)
    !       OMS            SINK LIMITED ASSIMILATION (MOL M-2 S-1)
    !       CO2S           CANOPY SURFACE CO2 CONCENTRATION (MOL MOL-1)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER, INTENT(IN   ) :: len
    REAL(KIND=r8)   , INTENT(IN   ) :: aparkk(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: vm    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: atheta(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: btheta(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: gammas(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: par   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: respc (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rrkk  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: omss  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: c3    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: c4    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: pco2i (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimn(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assim (len)

    !    local variables
    REAL(KIND=r8)      :: ome
    REAL(KIND=r8)      :: omc
    REAL(KIND=r8)      :: omp
    REAL(KIND=r8)      :: oms
    REAL(KIND=r8)      :: sqrtin
    INTEGER :: i

    !-----------------------------------------------------------------------
    !     CALCULATE ASSIMILATION RATE
    !
    !      OMC         (OMEGA-C): EQUATION (11) , SE-92A
    !      OME         (OMEGA-E): EQUATION (12) , SE-92A
    !      OMS         (OMEGA-S): EQUATION (13) , SE-92A
    !      ASSIMN      (A-N)    : EQUATION (14,15), SE-92A
    !-----------------------------------------------------------------------

    DO i = 1,len
       OMC = VM(i) *(PCO2I(i)-GAMMAS(i))/(PCO2I(i) + RRKK(i))*C3(i) &
            + VM(i) * C4(i)
       OME = PAR(i)*(PCO2I(i)-GAMMAS(i))/(PCO2I(i)+2.0_r8*GAMMAS(i))*C3(i)  &
            + PAR(i) * C4(i)
       SQRTIN= MAX( 0.0E0_r8, ( (OME+OMC)**2 - 4.0_r8*ATHETA(i)*OME*OMC ) )
       OMP  = ( ( OME+OMC ) - SQRT( SQRTIN ) ) / ( 2.0_r8*ATHETA(i) )
       OMS  = OMSS(i) * C3(i) + OMSS(i)*PCO2I(i) * C4(i)
       SQRTIN= MAX( 0.E0_r8, ( (OMP+OMS)**2 - 4.0_r8*BTHETA(i)*OMP*OMS ) )
       ASSIM(i) = ( ( OMS+OMP ) - SQRT( SQRTIN ) ) / &
            ( 2.0_r8*BTHETA(i) )
       ASSIMN(i)= ( ASSIM(i) - RESPC(i)) * APARKK(i)

    ENDDO
    !
    RETURN
  END SUBROUTINE CYCALC

  !
  !==================SUBROUTINE PHOSIB===================================
  !
  SUBROUTINE PHOSIB(pco2m,pco2ap,po2m,vmax0,tf,psur,green  ,    &
       tran,&
       ref,&
       gmudmu,&
       zlt,&
       cas_cap_co2,&
       tc,&
       ta,&
       trop,&
       trda ,  &
       trdm,&
       slti,&
       shti,&
       hltii,&
       hhti,&
       radn,&
       etc,&
       !etgs,&
       !wc      ,  &
       ea,&
       !em,&
       rb,&
       ra,&
       tm                                       ,  &
       effcon,rstfac,binter,gradm,assimn               ,  &
       rst,atheta,btheta,tgs,respcp                       ,  &
       aparkk,len,                                   &
       omepot,assimpot,assimci,antemp,assimnp,                    &
       wsfws,wsfht,wsflt,wci,whs,                                    &
       wags,&
       wegs,&
       aparc,&
       pfd,&
       assim,&
       td,&
       !www,                            &
       !wopt,&
       !zm,&
       !wsat,&
       !tg,&
       !soilscale,&
       !zmstscale,&
       zltrscale,&
       zmlscale,   &
       drst,&
       pdamp,&
       qdamp,&
       ecmass,&
       dtt,&
       bintc,&
       tprcor,&
       !soilq10,&
       ansqr,   &
       soilscaleold,&
       nsoil, &
       !forcerestore, &
       respg, &
       pco2c, &
       pco2i,   &
       pco2s,&
       co2cap,&
       cflux)

    IMPLICIT NONE

    !
    !
    !=======================================================================
    !
    !     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED   
    !     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
    !     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
    !     BASE REFERENCE IS SE-92A
    !
    !                          UNITS
    !                         -------
    !
    !      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
    !      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
    !      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
    !      EFFCON                                   : MOL CO2 MOL QUANTA-1
    !      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
    !      EVAPKG                                   : KG M-2 S-1
    !      Q                                        : KG KG-1
    !
    !                       CONVERSIONS
    !                      -------------
    !
    !      1 MOL H2O           = 0.018 KG
    !      1 MOL CO2           = 0.044 KG
    !      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
    !      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
    !pl the next line applies to the Ci to Cs pathway
    !      GS  (CO2)           = GS (H2O) * 1./1.6
    !pl 44.6 is the number of moles of air per cubic meter
    !      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
    !      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
    !      MM  (MOLAIR/MOLH2O) = 1.611
    !
    !
    !                         OUTPUT
    !                      -------------
    !
    !      ASSIMN              = CANOPY NET ASSIMILATION RATE
    !      EA                  = CANOPY AIR SPACE VAPOR PRESSURE
    !      1/RST               = CANOPY CONDUCTANCE
    !      PCO2I               = INTERNAL CO2 CONCENTRATION
    !      RESPC               = CANOPY RESPIRATION
    !      RESPG               = GROUND RESPIRATION
    !
    !----------------------------------------------------------------------
    !
    !         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
    !         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
    !         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
    !         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
    !
    !-----------------------------------------------------------------------
    !

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1) 
    !       RST            CANOPY RESISTANCE (S M-1)
    !       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS 
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
    !       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
    !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
    !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
    !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !     Modifications:
    !       - gs (stomatal conductance reduced for freezing soils per Jim Collatz
    !         dd 950221      
    !
    !      Modified for multitasking - introduced gather/scatter indices
    !          - DD 951206
    !
    !itb   Added in pco2c (chloroplast partial co2) for neil's fractionation
    !itb   calculations
    !itb       - IB Sep99

    !
    !     input arrays:
    INTEGER , INTENT(IN   ) :: len
    !INTEGER , INTENT(IN   ):: nsib
    INTEGER , INTENT(IN   ) :: nsoil
    !LOGICAL , INTENT(IN   ) :: forcerestore!is not used

    REAL(KIND=r8)    , INTENT(IN) :: vmax0(len)
    REAL(KIND=r8)    , INTENT(IN) :: psur       (len)
    REAL(KIND=r8)    , INTENT(IN) :: green      (len)
    REAL(KIND=r8)    , INTENT(IN) :: gmudmu     (len)
    REAL(KIND=r8)    , INTENT(IN) :: zlt        (len)
    REAL(KIND=r8)    , INTENT(IN) :: cas_cap_co2(len)
    REAL(KIND=r8)    , INTENT(IN) :: tc         (len)
    REAL(KIND=r8)    , INTENT(IN) :: ta         (len)
    REAL(KIND=r8)    , INTENT(IN) :: trop       (len)
    REAL(KIND=r8)    , INTENT(IN) :: trda       (len)
    REAL(KIND=r8)    , INTENT(IN) :: slti       (len)
    REAL(KIND=r8)    , INTENT(IN) :: shti       (len)
    REAL(KIND=r8)    , INTENT(IN) :: hltii      (len)
    REAL(KIND=r8)    , INTENT(IN) :: hhti       (len)
    !REAL(KIND=r8)    , INTENT(IN) :: etgs       (len)!is not used
    !REAL(KIND=r8)    , INTENT(IN) :: wc         (len)!is not used
    !REAL(KIND=r8)    , INTENT(IN) :: em         (len)!is not used
    REAL(KIND=r8)    , INTENT(IN) :: ra         (len)
    REAL(KIND=r8)    , INTENT(IN) :: rb         (len)
    REAL(KIND=r8)    , INTENT(IN) :: tm         (len)
    REAL(KIND=r8)    , INTENT(IN) :: effcon     (len)
    REAL(KIND=r8)    , INTENT(IN) :: binter     (len)
    REAL(KIND=r8)    , INTENT(IN) :: gradm      (len)
    REAL(KIND=r8)    , INTENT(IN) :: atheta     (len)
    REAL(KIND=r8)    , INTENT(IN) :: btheta     (len)
    REAL(KIND=r8)    , INTENT(IN) :: tgs        (len)
    REAL(KIND=r8)    , INTENT(IN) :: respcp     (len)
    REAL(KIND=r8)    , INTENT(IN) :: radn       (len,2,2)
    REAL(KIND=r8)    , INTENT(IN) :: ecmass     (len)
    REAL(KIND=r8)    , INTENT(IN) :: trdm       (len)
    REAL(KIND=r8)    , INTENT(IN) :: etc        (len)
    REAL(KIND=r8)    , INTENT(IN) :: aparc      (len)
    !REAL(KIND=r8)    , INTENT(IN) :: www        (len,2)!is not used
    !REAL(KIND=r8)    , INTENT(IN) :: tg         (len)!is not used
    REAL(KIND=r8)    , INTENT(INOUT) :: rst        (len)
    REAL(KIND=r8)    , INTENT(OUT  ) :: cflux      (len)
    REAL(KIND=r8)    , INTENT(IN) :: pdamp
    REAL(KIND=r8)    , INTENT(IN) :: qdamp
    REAL(KIND=r8)    , INTENT(IN) :: dtt
    REAL(KIND=r8)    , INTENT(IN   ) :: pco2m      (len)
    REAL(KIND=r8)    , INTENT(INOUT) :: pco2ap     (len)
    REAL(KIND=r8)    , INTENT(IN   ) :: tf
    REAL(KIND=r8)    , INTENT(IN   ) :: po2m

    !     output arrays:


    REAL(KIND=r8)    , INTENT(OUT) :: assimn     (len)
    REAL(KIND=r8)    , INTENT(IN) :: ea         (len)
    REAL(KIND=r8)    , INTENT(INOUT):: rstfac     (len,4)
    REAL(KIND=r8)    , INTENT(OUT) :: pco2i      (len)
    REAL(KIND=r8)    , INTENT(IN) :: respg      (len)
    REAL(KIND=r8)    , INTENT(OUT) :: drst       (len)

    !zz new diagostics 10/14/92 
    !
    ! output arrays

    REAL(KIND=r8)   , INTENT(OUT  ) :: omepot   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimpot (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimci  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimnp  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: whs      (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: antemp   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsfws    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsflt    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wci      (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wags     (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wegs     (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pfd      (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: td       (len,nsoil)
    REAL(KIND=r8)   , INTENT(OUT  ) :: zmlscale (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assim    (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: wopt     (len)!is not used
    !REAL(KIND=r8)   , INTENT(IN   ) :: zm       (len)!is not used
    !REAL(KIND=r8)   , INTENT(IN   ) :: wsat     (len)!is not used
    !REAL(KIND=r8)   , INTENT(IN   ) :: soilscale(len,nsoil+1)!is not used
    !REAL(KIND=r8)   , INTENT(IN   ) :: zmstscale(len,2)!is not used
    REAL(KIND=r8)   , INTENT(OUT  ) :: zltrscale(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tprcor   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: bintc    (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: soilq10  (len,nsoil+1)!is not used
    REAL(KIND=r8)   , INTENT(OUT  ) :: ansqr    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: soilscaleold(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pco2c    (len) !chloroplast pco2
    REAL(KIND=r8)   , INTENT(INOUT) :: aparkk   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pco2s    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: co2cap   (len)

    REAL(KIND=r8)    :: xgah2o   (len) 
    REAL(KIND=r8)    :: xgco2m   (len)
    REAL(KIND=r8)    :: blinder2   (len)

    !     work arrays:

    REAL(KIND=r8)    :: tran     (len,2,2)
    REAL(KIND=r8)    :: ref      (len,2,2)
    REAL(KIND=r8)    :: respc    (len)
    REAL(KIND=r8)    :: wsfht    (len)
    REAL(KIND=r8)    :: PCO2Y    (len,6)
    REAL(KIND=r8)    :: EYY      (len,6)
    REAL(KIND=r8)    :: assimny  (len,6)
    REAL(KIND=r8)    :: assimy   (len,6)
    REAL(KIND=r8)    :: c3       (len)
    REAL(KIND=r8)    :: c4       (len)
    REAL(KIND=r8)    :: RANGE    (len)
    REAL(KIND=r8)    :: gammas   (len)
    REAL(KIND=r8)    :: gah2o    (len)
    REAL(KIND=r8)    :: gbh2o    (len)
    REAL(KIND=r8)    :: par      (len)
    REAL(KIND=r8)    :: rrkk     (len)
    REAL(KIND=r8)    :: omss     (len)
    REAL(KIND=r8)    :: vm       (len)
    REAL(KIND=r8)    :: gsh2o    (len)
    REAL(KIND=r8)    :: templ    (len)
    REAL(KIND=r8)    :: temph    (len)
    REAL(KIND=r8)    :: qt       (len)
    REAL(KIND=r8)    :: co2s     (len)
    REAL(KIND=r8)    :: scatp    (len)
    REAL(KIND=r8)    :: scatg    (len)
    REAL(KIND=r8)    :: park     (len)
    REAL(KIND=r8)    :: respn    (len)
    REAL(KIND=r8)    :: zkc      (len)
    REAL(KIND=r8)    :: zko      (len)
    REAL(KIND=r8)    :: spfy     (len)
    REAL(KIND=r8)    :: co2a     (len)
    REAL(KIND=r8)    :: co2m     (len)


    INTEGER :: icconv(len)
    INTEGER :: igath(len)

    REAL(KIND=r8)    :: soilfrz(len)
    REAL(KIND=r8)    :: cwsflt
    REAL(KIND=r8)    :: cwsfht
    REAL(KIND=r8)    :: cwsfws
    REAL(KIND=r8)    :: ccoms
    REAL(KIND=r8)    :: ccomc
    REAL(KIND=r8)    :: ascitemp
    REAL(KIND=r8)    :: dompdomc
    REAL(KIND=r8)    :: omsci
    REAL(KIND=r8)    :: ompci
    REAL(KIND=r8)    :: omcci
    REAL(KIND=r8)    :: omcpot
    REAL(KIND=r8)    :: omppot
    REAL(KIND=r8)    :: sqrtin
    REAL(KIND=r8)    :: omspot
    REAL(KIND=r8)    :: pco2ipot
    REAL(KIND=r8)    :: ohtp2
    REAL(KIND=r8)    :: sttp2
    REAL(KIND=r8)    :: gsh2oinf(len)
    REAL(KIND=r8)    :: h2osrh
    REAL(KIND=r8)    :: h2os
    REAL(KIND=r8)    :: ecmole
    REAL(KIND=r8)    :: h2oa
    REAL(KIND=r8)    :: h2oi
    REAL(KIND=r8)    :: dtti
    REAL(KIND=r8)    :: pco2in
    REAL(KIND=r8)    :: pco2a
    REAL(KIND=r8)    :: soilfrztd
    REAL(KIND=r8)    :: soilfrztg
    REAL(KIND=r8)    :: Min_assimn
    REAL(KIND=r8)    :: Min_gsmin

    INTEGER i, ic1, ic



    !pl introduce a co2 capacity 
    !pl this will basically be the mass of air under the top of the canopy (in
    !pl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.

    !pl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
    !pl carbon are in m2 s / mol air

    !pl one mole of gas occupies 22.4 cubic dm
    !pl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
    !pl the units for the carbon capacity are mol air /m2. 
    !pl (e.g. here 893 moles if thickness of the layer is 20m)
    !pl this means that the units for pc02ap should be mol co2 / mol air, but
    !pl it is also possible to keep just co2 pressure and convert

   icconv=0
    DO i = 1,len                !  LOOP OVER GRIDPOINT


       TPRCOR(i) = TF*PSUR(i)*100.0_r8/1.013E5_r8
       co2cap(i) = cas_cap_co2(i) * 44.6_r8 * tprcor(i)/ta(i)     ! moles air / m2

       !pl this needs to be modified as in sibslv3 to automatically use the 
       !pl thickness of the canopy air space. 

       !
       !----------------------------------------------------------------------
       !
       !pl        RESPG(i) = 0. E -6 ! fixed respiration at 5 micromoles
       !pl   no longe rused since we now have respsib
       !
       !----------------------------------------------------------------------
       !
       IF( EFFCON(i) .GT. 0.07_r8 ) THEN
          C3(i) = 1.0_r8
       ELSE
          C3(i) = 0.0_r8
       ENDIF
       C4(i)     = 1.0_r8 - C3(i)

       !
       !-----------------------------------------------------------------------
       !
       !     CALCULATION OF CANOPY PAR USE PARAMETER.
       !
       !      APARKK      (PI)     : EQUATION (31) , SE-92A
       !-----------------------------------------------------------------------
       !  
       SCATP(I) =     GREEN(i)   *            &
            ( TRAN(i,1,1) + REF(i,1,1) )   &
            +( 1.0_r8-GREEN(i) ) *             &
            ( TRAN(i,1,2) + REF(i,1,2) )
       SCATG(i) = TRAN(i,1,1) + REF(i,1,1)
       PARK(i) = SQRT(1.0_r8-SCATP(i)) * GMUDMU(i)
       !
       ! Collatz-Bounoua commented the calculation of  aparc
       ! replaced it with theone calculated in new_mapper.
       !
       !b        APARC(i) = 1. - EXP ( -PARK(i)*ZLT(i) )   ! lahouari
       !
       APARKK(i)   = APARC(i) / PARK(i) * GREEN(i)
       !-----------------------------------------------------------------------
       !
       !     Q-10 AND STRESS TEMPERATURE EFFECTS
       !
       !      QT          (QT)    : TABLE (2)     , SE-92A
       !-----------------------------------------------------------------------
       !
       qt(i) = 0.1_r8*( TC(i) - TROP(i) )
       RESPN(i) = RESPCP(i) * VMAX0(i) * RSTFAC(i,2)

       !itb...patch to prevent underflow if temp is too cool...
       IF(TC(i) >= TRDM(i))THEN
          RESPC(i) = RESPN(i) * 2.0_r8**qt(i) &
               /( 1.0_r8 + EXP( TRDA(i)*(TC(i)-TRDM(i))))
       ELSE
          RESPC(i) = RESPN(i) * 2.0_r8**qt(i)
       ENDIF

       VM(i) = VMAX0(i) * 2.1_r8**qt(i)
       TEMPL(i) = 1.0_r8 + EXP(SLTI(i)*(HLTIi(i)-TC(i)))
       TEMPH(i) = 1.0_r8 + EXP(SHTI(i)*(TC(i)-HHTI(i)))
       RSTFAC(i,3) = 1.0_r8/( TEMPL(i)*TEMPH(i))
       VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i) &
            + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)
       !
       !-----------------------------------------------------------------------
       !
       !     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
       !     COMPENSATION POINT       
       !
       !      ZKC          (KC)     : TABLE (2)     , SE-92A
       !      ZKO          (KO)     : TABLE (2)     , SE-92A
       !      SPFY         (S)      : TABLE (2)     , SE-92A
       !      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
       !      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
       !      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
       !-----------------------------------------------------------------------
       !
       ZKC(i) = 30.0_r8 * 2.1_r8**qt(i)
       ZKO(i) = 30000.0_r8 * 1.2_r8**qt(i)
       SPFY(i) = 2600.0_r8 * 0.57_r8**qt(i)
       GAMMAS(i) = 0.5_r8 * PO2M/SPFY(i) * C3(i)
       PFD(i)    = 4.6E-6_r8 * GMUDMU(i)* &
            ( RADN(i,1,1)+RADN(i,1,2) )

       !
       !pl these here all go from being m/s to being mol/ (m2 sec)
       !PK       GSH2O(i)  = 1.0_r8/RST(i) * 44.6_r8*TPRCOR(i)/TC(i)
       GBH2O(i)  = 0.5_r8/RB(i) * 44.6_r8*TPRCOR(i)/TC(i)
       GAH2O(i)  = 1.0_r8/RA(i) * 44.6_r8*TPRCOR(i)/TM(i)

       xgah2o(i) = MAX(0.466E0_r8, gah2o(i) )
       xgco2m(i) = 4000.0_r8 * vmax0(i)
       !
       !itb...this is changed slightly from older version of code...
       RRKK(i)   = ZKC(i)*( 1.0_r8 + PO2M/ZKO(i) ) * C3(i) &
            + VMAX0(i)/5._r8* ( 1.8_r8**qt(i)) * C4(i)
       !        RRKK(i)   = ZKC(i)*( 1.0_r8 + PO2M/ZKO(i) ) * C3(i)
       !     &               + VMAX0(i)*200.0_r8/psur(i)* ( 1.8**qt(i)) * C4(i)
       PAR(i)    = pfd(i)*EFFCON(i)*( 1.0_r8-SCATG(i) )
       soilfrztg = 1.0_r8+EXP(-1.5_r8 * (MAX(270.0E0_r8,tgs(i))-273.16_r8) )
       soilfrztd = 1.0_r8+EXP(-1.5_r8 * (MAX(270.0E0_r8,td (i,nsoil))-273.16_r8) )
       soilfrz(i) = MAX(1.0_r8/soilfrztg, 1.0_r8/soilfrztd)
       soilfrz(i) = MAX( soilfrz(i), 0.05E0_r8)
       
       ! bintc(i)- smallest canopy stomatal conductance needs to be passed in here.

       BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)*RSTFAC(i,2)*soilfrz(i)

       blinder2(i)= bintc(i) * tc(i) / ( 44.6_r8 * tprcor(i)) 
       IF(rst(i)<=0.0_r8)THEN
          ! rst(i)=MIN(MAX( ten, 1.0_r8/blinder2(i) ),1.0e5_r8)
          rst(i)=MIN(MAX( 10.0_r8, 1.0_r8/blinder2(i) ),1.0e5_r8)
       END IF
       GSH2O(i)  = 1.0_r8/RST(i) * 44.6_r8*TPRCOR(i)/TC(i)


       !itb...this is changed slightly from older version of code...
       OMSS(i)   = ( VMAX0(i)/2.0_r8 ) * ( 1.8_r8**qt(i) )    & 
            /TEMPL(i) * RSTFAC(i,2) * C3(i)  &
            + RRKK(i) *RSTFAC(i,2) * C4(i)
       !        OMSS(i)   = ( VMAX0(i)/2.0_r8 ) * ( 1.8_r8**qt(i) )
       !     &                  /TEMPH(i) * RSTFAC(i,2) * C3(i)
       !     &                  + RRKK(i) * RSTFAC(i,3)*RSTFAC(i,2) * C4(i)
       !
       !-----------------------------------------------------------------------
       !
       !     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
       !     ASSIMILATION RATE.
       !
       !-----------------------------------------------------------------------


       RANGE(i)    = PCO2M(i) * ( 1.0_r8 - 1.6_r8/GRADM(i) ) - GAMMAS(i)
       icconv(i) = 1

    ENDDO

    !
    DO IC = 1, 6
       DO i = 1,len        ! LOOP OVER GRIDPOINT
          PCO2Y(i,IC) = 0.0_r8
          EYY(i,IC) = 0.0_r8
       ENDDO
    END DO
    !

    !pl beginning of PL's setup

    !      do i=1,len
    !       gah2o(i) =  1.0_r8 / MAX(0.446_r8 E 0,GAH2O(i))
    !      enddo

    !Bio...HERE IS THE ITERATION LOOP.
    !Bio...
    !Bio...We iterate on PCO2C-sortin makes a 'first guess' at
    !Bio...then orders PCO2C/Error pairs on increasing error size,
    !Bio...then uses a combination of linear/quadratic fit to obtain 
    !Bio...the 'next best guess' as iteration count increases.
    !Bio...CYCALC uses that value of PCO2C to get the next value 
    !Bio...of ASSIMN. CO2A and CO2S follow.

    DO IC = 1, 6
       !
       CALL       SORTIN( EYY, PCO2Y, RANGE, GAMMAS, ic,len )

       CALL       CYCALC( APARKK, VM, ATHETA, BTHETA,par,     &
            GAMMAS, RESPC, RRKK, OMSS, C3, C4,  &
            PCO2Y(1:len,ic), assimny(1:len,ic), assimy(1:len,ic),&
            len  )
       !

       DO i = 1,len

          !pl first diagnose the current CO2 flux in mol / (m2 s)
          !pl this is a modified ra that will get us the right units
          !pl in the conservation equation. its units are m2 s / (mol_air)

          !        resa(i) =   1. / MAX(0.446 E 0,GAH2O(i))

          !        resb(i) =   1.4/GBH2O(i)

          !        resc(i) =   1.6/GSH2O(i) 

          !pl now prognose the new CAS CO2 according to flux divergence
          !pl we are going to do this in mol C / mol air (same as PaC/PaAir)

          CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.0_r8)
          co2m(i)    = pco2m(i)  /   (PSUR(i)*100.0_r8)   

          CO2A(i)   = (  CO2A(i) + (dtt/co2cap(i)) *  &
               (respg(i) - assimny(i,ic)       &
               +co2m(i)*gah2o(i)        ) )    &
               / (1+dtt*gah2o(i)/ co2cap(i) ) 


          !        PCO2A   =  PCO2Ap(i) + (dtt/co2cap) 
          !     &            * (  (PSUR(i) *100.0 * respg(i))
          !     &            +    (PCO2Y(i,IC)/(resb+resc))
          !     &            +    (pco2m(i)   /       resa)  )
          !     &            / (1.0_r8+ (dtt/co2cap(i)) * (1.0_r8/resa + 1.0_r8/(resb+resc)) )


          !pl        PCO2A = PCO2M(i) - (1.4_r8/MAX(0.446_r8 E 0,GAH2O(i)) * 
          !pl     >             (ASSIMNy(i,ic) - RESPG(i))* PSUR(i)*100.0_r8)


          pco2a = co2a(i) * psur(i) * 100.0_r8

          PCO2S(i) = PCO2A - (1.4_r8/GBH2O(i) * ASSIMNy(i,ic) &
               * PSUR(i)*100.0_r8)
          !        PCO2IN   = PCO2S(i) - (1.6_r8/GSH2O(i) * ASSIMNy(i,ic) &
          !                  * PSUR(i)*100.0_r8)
          !hanging to iterate on pco2c instead of pco2i
          PCO2IN = PCO2S(i) - ASSIMNy(i,ic)*psur(i)*100.0_r8* &
               (1.6_r8/gsh2o(i) +1.0_r8/xgco2m(i))
          EYY(i,IC) = PCO2Y(i,IC) - PCO2IN
       ENDDO
       !
       
       icconv = 1
       IF(ic.GE.2) THEN
          ic1 = ic-1
          DO i = 1,len        ! LOOP OVER GRIDPOINT
             IF(ABS(eyy(i,ic1)).GE.0.1_r8)THEN
                icconv(i) = ic
             ELSE
                icconv(i) = ic1
                eyy  (i,ic) = eyy(i,ic1)
                pco2y(i,ic) = pco2y(i,ic1)
             ENDIF
          ENDDO
       ENDIF
       !
    END DO
    !
    !

    !pl end of PL's setup

    DO i = 1,len        ! LOOP OVER GRIDPOINT
       icconv(i) = MIN(icconv(i),6)
       igath(i) = i+(icconv(i)-1)*len
    ENDDO


    DO i = 1,len         ! LOOP OVER GRIDPOINT

       pco2c(i)  = pco2y  (i,icconv(i))
       assimn(i) = assimny(i,icconv(i))
       assim(i)  = assimy (i,icconv(i))



       !        pco2i(i)  = pco2y  (i,icconv(i))
       !        assimn(i) = assimny(i,icconv(i))
       !        assim(i)  = assimy (i,icconv(i))

       pco2i(i) = pco2c(i) + assimn(i)/xgco2m(i)*psur(i)*100.0_r8
       pco2s(i) = pco2i(i) + assimn(i)/gsh2o(i)*psur(i)*100.0_r8


       !        pco2c(i) = pco2i(i) - assimn(i)/xgco2m(i)*psur(i)*100.0_r8

       !pl now do the real C_A forecast with the iterated fluxes_r8

       CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.0_r8)
       co2m(i)    = pco2m(i)  /   (PSUR(i)*100.0_r8)   

       CO2A(i) = (CO2A(i) + (dtt/co2cap(i)) *  &
            (respg(i) - assimn(i)       &
            +co2m(i)*gah2o(i) ) )       &
            / (1+dtt*gah2o(i)/co2cap(i))
       !pl go back from molC / mol air to Pascals 

       pco2ap(i) = co2a(i) * psur(i) * 100.0_r8

       !itb...carbon flux between CAS and reference level
       cflux(i) = gah2o(i)*(co2a(i)-co2m(i))

       !        PCO2Ap(i) =  PCO2Ap(i) + (dtt/co2cap) 
       !     &            * (  (PSUR(i) *100.0_r8 * respg(i))
       !     &            +    (PCO2i(i)   /(resb+resc))
       !     &            +    (pco2m(i)   /       resa)  )
       !     &            / (1.0_r8+ (dtt/co2cap) * (1.0_r8/resa + 1.0_r8/(resb+resc)) )
    ENDDO

    !
    dtti = 1.0_r8/dtt
    GSH2OINF=-10.0_r8
    Min_assimn=-1.0_r8
    Min_gsmin=-0.007_r8
    DO 
       DO i = 1,len        ! LOOP OVER GRIDPOINT
          IF(GSH2OINF(i) < Min_gsmin  ) THEN
             !czzggrst5 - new code
             H2OI   = ETC(i)/PSUR(i)
             H2OA   =  EA(i)/PSUR(i)
             ECMOLE = 55.56_r8 * ECMASS(i) * dtti  ! ecmass must be computed and passed in
             H2OS = H2OA + ECMOLE / GBH2O(i)
             H2OS  = MIN( H2OS, H2OI )
             H2OS  = MAX( H2OS, 1.0E-7_r8)
             H2OSRH = H2OS/H2OI
             !  need qdamp and pdamp calculated and passed to here!
             !pl        CO2S(i) = MAX(PCO2S(I),PCO2M*0.5_r8) / (PSUR(i)*100.0_r8)

             !pl I have relaxed this condition to 1/10 of previous. The old way made
             !pl the CO2 on top of the leaves always at least 1/2 of the value at the
             !pl reference level.

             CO2S(i) = MAX(PCO2S(I),PCO2M(i)*0.05_r8) / (PSUR(i)*100.0_r8)

             !pl Ball-Berry equation right here !

             GSH2OINF(i) = (GRADM(i) * MAX(ASSIMN(i),Min_assimn) * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)

             !pl this is the change in stomatal resistance

             DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF(i))/ &
                (PDAMP*GSH2O(i)+QDAMP*GSH2OINF(i)))


             !pl this is the 'would be change' if we did not use the damping factor.

             !        rstnew = (1./gsh2oinf(i)) * 44.6 * tprcor(i)/tc(i)
             !        DRST(i) = rstnew - RST(i)

             !
             RSTFAC(i,1) = H2OS/H2OI
             RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)
          END IF
       END DO
       IF (ANY(GSH2OINF < Min_gsmin)) THEN
          Min_assimn = Min_assimn/10.0_r8
       ELSE
         EXIT
       END IF
    END DO
    !
    !Z CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
    !-----------------------------------------------------------------------
    !  INPUTS: PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
    !    C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
    !    TEMPL,RSTFAC(i,1),VM(i),ASSIM,GSH20(i),EFFCON(i),QT,GAMMAS(i),
    !    PFD(i)
    !

    sttp2 = 73.0_r8**0.2_r8
    ohtp2 = 100.0_r8**0.2_r8
    DO i = 1,len
       !-----------------------------------------------------------------------
       ! CALCULATION OF POTENTIAL ASSIMILATION
       !-----------------------------------------------------------------------
       !
       ! Make assimn a top leaf, not the canopy.
       ASSIMNp(i) = ASSIMN(i) / APARKK(i)
       !
       ! Bottom stopped assim.
       ANTEMP(i) = MAX(0.E0_r8,ASSIMNp(i))
       !
       ! Potential intercellular co2.
       PCO2IPOT = PSUR(i)*100.0_r8*(co2s(i)-(1.6_r8*ASSIMNp(i)/ &
            ((GRADM(i)*ANTEMP(i)/co2s(i))+BINTC(i))))
       !
       ! Potential rubisco limitation.
       OMCPOT = VMAX0(i)*2.1_r8**qt(i)*((PCO2IPOT-GAMMAS(i))/ &
            (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
       !
       ! Potential light limitation.
       OMEPOT(i) = PAR(i)*((PCO2IPOT-GAMMAS(i))/ &
            (PCO2IPOT+2.0_r8*GAMMAS(i))*C3(i) + C4(i))
       !
       ! Quad 1.
       SQRTIN = MAX(0.E0_r8,((OMEPOT(i)+OMCPOT)**2- &
            4.0_r8*ATHETA(i)*OMEPOT(i)*OMCPOT))
       !
       ! Quad 1. Intermediate  top leaf photosynthesis.
       OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.0_r8*ATHETA(i))
       !
       ! Potential sink or pep limitation.
       OMSPOT = (VMAX0(i)/2.0_r8)*(1.8_r8**qt(i))*C3(i) &
            + RRKK(i)*PCO2IPOT*C4(i)
       !
       ! Quad 2.
       SQRTIN=MAX(0.0E0_r8,((OMPPOT+OMSPOT)**2-4.0_r8*BTHETA(i)* &
            OMPPOT*OMSPOT))
       !
       ! Quad 2. Final Potential top leaf photosynthesis.
       ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.0_r8*BTHETA(i))
       !
       !-----------------------------------------------------------------------
       ! CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
       !-----------------------------------------------------------------------
       !
       ! Stressed rubisco limitation.
       OMCCI = VM(i)*((PCO2IPOT-GAMMAS(i))/(PCO2IPOT+RRKK(i))*C3(i) &
            + C4(i))
       !
       ! Quad 1.
       SQRTIN = MAX(0.0E0_r8,(OMEPOT(i)+OMCCI)**2- &
            4.0_r8*ATHETA(i)*OMEPOT(i)*OMCCI)
       !
       ! Quad 1. Intermediate stress limited top leaf photosynthesis.
       OMPCI = ((OMEPOT(i)+OMCCI)-SQRT(SQRTIN))/(2.0_r8*ATHETA(i))
       !
       ! Stressed sink or pep limitation.
       OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
       !
       ! Quad 2.
       SQRTIN = MAX(0.0E0_r8,(OMPCI+OMSCI)**2-4.0_r8*BTHETA(i)*OMPCI*OMSCI)
       ! 
       ! Quad 2. Final stress limited top leaf photosynthesis.
       ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.0_r8*BTHETA(i))
       !
       !-----------------------------------------------------------------------
       ! CALCULATION OF CONTROL COEFFICIENTS
       !-----------------------------------------------------------------------
       !
       ! Intermediate.
       DOMPDOMC = (OMPCI-OMEPOT(i))/ &
            (2.0_r8*ATHETA(i)*OMPCI-OMCCI-OMEPOT(i))
       !
       ! Bottom stopped final stress limited top leaf photosynthesis.
       ASCITEMP = MAX(ASSIMCI(i),1.0E-12_r8)
       !
       ! Rubisco control coefficient.
       CCOMC = (DOMPDOMC*(ASSIMCI(i)-OMSCI)/ &
            (2.0_r8*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMCCI/ASCITEMP
       !
       ! Sink or pep control coefficient.
       CCOMS = ((ASSIMCI(i)-OMPCI)/ &
            (2.0_r8*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMSCI/ASCITEMP
       !
       !-----------------------------------------------------------------------
       !  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
       !-----------------------------------------------------------------------
       ! Canopy values (overwrites top leaf).
       ! 
       OMEPOT(i) = OMEPOT(i)*APARKK(i)
       ASSIMPOT(i) = ASSIMPOT(i)*APARKK(i)
       ASSIMCI(i) = ASSIMCI(i)*APARKK(i)
       ASSIM(i) = ASSIM(i)*APARKK(i)
       ANTEMP(i) = ANTEMP(i)*APARKK(i)
       ANSQR(i) = ANTEMP(i)*ANTEMP(i)
       ASSIMNp(i) = ASSIMNp(i)*APARKK(i)
       !
       !-----------------------------------------------------------------------
       ! OUTPUT:  WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
       !-----------------------------------------------------------------------
       !
       ! Water stress.
       WSFWS(i) = ASSIMPOT(i)*(1.0_r8-RSTFAC(i,2))*(CCOMC+CCOMS)
       !
       ! High temperature stress.
       WSFHT(i) = ASSIMPOT(i)*(1.0_r8-1.0_r8/TEMPH(i))*CCOMC
       !
       ! Low temperature stress.
       WSFLT(i) = ASSIMPOT(i)*(1.0_r8-1.0_r8/TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
       !
       !  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
       cwsfws = (1.0_r8-RSTFAC(i,2))*(CCOMC+CCOMS)
       IF(cwsfws.GT.1.0_r8 .OR. cwsfws.LT.0.) wsfws(i)=0.0_r8
       cwsfht = (1.0_r8-1.0_r8/TEMPH(i))*CCOMC
       IF(cwsfht.GT.1.0_r8 .OR. cwsfht.LT.0.0_r8) wsfht(i)=0.0_r8
       cwsflt = (1.0_r8-1.0_r8/TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
       IF(cwsflt.GT.1.0_r8 .OR. cwsflt.LT.0.0_r8) wsflt(i)=0.0_r8

       !
       ! Intermediate assimilation weighted Ci.
       WCI(i) = ANTEMP(i)*PCO2I(i)
       !
       ! Intermediate assimilation weighted relative humidty stress factor.
       WHS(i) = ANTEMP(i)*RSTFAC(i,1)
       !
       ! Intermediate assimilation weighted stomatal conductance.
       WAGS(i) = GSH2O(i)*ANTEMP(i)
       !
       ! Intermediate evaporation weighted stomatal conductance.(Step 1.
       !   Step 2 after subroutine updat2)
       WEGS(i) = GSH2O(i)
       !
       !      bl = (((100.* www(i,1))**0.2 - sttp2)/(sttp2 - ohtp2))**2
       !      bl = min(bl,10. E 0)
       !      zmlscale(i) = 0.8*0.75**bl + 0.2
       !
       !      zltrscale(i) = (exp(0.0693*(tgs(i)-298.15)))
       !      if (zltrscale(i) .lt. 0.17 ) zltrscale(i)= 0.
       !      zltrscale(i) = zltrscale(i) * zmlscale(i)

       !      make these zero until canned from the qp diagnostic, then scrap
       soilscaleold(i) = 0.0_r8
       zmlscale(i) = 0.0_r8
       zltrscale(i) = 0.0_r8
    ENDDO
    !

    RETURN
  END SUBROUTINE PHOSIB

  SUBROUTINE VNTLAT                                               &
       (grav, tice,                                              &
       pr0, ribc, vkrmn, delta, dtt, tc, tg, ts, ps, zlt,   &
       www, tgs, etc, etg,snoww,                              &
       rstfac, rsoil, hr, wc, wg, snofac,              &
       sh, z0, spdm, sha, zb, ros,cas_cap_co2,              &
       cu, ra, thvgm, rib, ustar, ventmf, thm, tha, z2, d,  &
       fc, fg, rbc, rdc,gect,geci,gegs,gegi,                      &
       respcp, rb, rd, rds,  rst, rc, ecmass,              &
       ea,  hrr, assimn, bintc, ta, pco2m, po2m, vmax0,  &
       green, tran, ref, gmudmu, trop, trda, trdm, slti,    &
       shti, hltii, hhti, radn, effcon, binter, gradm,      &
       atheta, btheta, aparkk, wsfws, wsfht, wsflt, wci,    &
       whs, omepot, assimpot, assimci, antemp, assimnp,     &
       wags, wegs, aparc, pfd, assim, td, wopt, zm, wsat,   &
       soilscale, zmstscale, zltrscale, zmlscale, drst,     &
       soilq10, ansqr, soilscaleold,                        &
       len, nsoil, forcerestore, dotkef,              &
       thgeff, shgeff, tke, ct, louis,          &
       respg, respfactor, pco2ap, pco2i, pco2c, pco2s,      &
       co2cap,cflux)

    !
    !
    !        - optimized subroutine phosib
    !          dd 92.06.10
    !
    !
    IMPLICIT NONE

    !    argument list declarations  


    !INTEGER, INTENT(IN   ) :: nsib
    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: nsoil

    REAL(KIND=r8)   , INTENT(IN   ) :: grav
    REAL(KIND=r8)   , INTENT(IN   ) :: tice
    REAL(KIND=r8)   , INTENT(IN   ) :: pr0
    REAL(KIND=r8)   , INTENT(IN   ) :: ribc
    REAL(KIND=r8)   , INTENT(IN   ) :: vkrmn
    REAL(KIND=r8)   , INTENT(IN   ) :: delta
    REAL(KIND=r8)   , INTENT(IN   ) :: dtt
    REAL(KIND=r8)   , INTENT(IN   ) :: tc   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ts   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ps   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zlt  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: www  (len,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: tgs  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: etc  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: etg  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: snoww(len,2) ! snow cover (veg and ground) (m)
    REAL(KIND=r8)   , INTENT(INOUT) :: rstfac(len,4)
    REAL(KIND=r8)   , INTENT(IN   ) :: rsoil(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hr   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: wc   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: wg   (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: areas(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: snofac
    REAL(KIND=r8)   , INTENT(IN   ) :: sh   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: z0   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: spdm (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: sha  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zb   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ros  (len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: cas_cap_co2(len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: cu   (len)
    REAL(KIND=r8)   , INTENT(INOUT) :: ra   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: thvgm(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rib  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ustar(len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: ventmf(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: thm  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tha  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: z2   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: d    (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: fc   (len)
    REAL(KIND=r8)   , INTENT(INOUT) :: fg   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rbc  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rdc  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: gect (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: geci (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: gegs (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: gegi (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: respcp(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: rb   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rd   (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: rds  (len) 
    !REAL(KIND=r8)   , INTENT(IN   ) :: bps  (len)
    REAL(KIND=r8)   , INTENT(INOUT) :: rst  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: rc   (len)
    REAL(KIND=r8)   , INTENT(INOUT) :: ecmass(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ea   (len) 
    !REAL(KIND=r8)   , INTENT(IN   ) :: em   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hrr  (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimn(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: bintc(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ta   (len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: pco2m(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: po2m
    REAL(KIND=r8)   , INTENT(IN   ) :: vmax0(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: green(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tran (len,2,2) 
    REAL(KIND=r8)   , INTENT(IN   ) :: ref  (len,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: gmudmu(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: trop (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: trda (len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: trdm (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: slti (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: shti (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hltii(len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: hhti (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: radn (len,2,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: effcon(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: binter(len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: gradm (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: atheta(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: btheta(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: aparkk(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsfws (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsfht (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wsflt (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: wci   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: whs   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: omepot(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimpot(len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimci(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: antemp (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assimnp(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: wags(len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: wegs(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: aparc(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pfd(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: assim(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: td(len,nsoil) 
    REAL(KIND=r8)   , INTENT(IN   ) :: wopt(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: zm  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: wsat(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: soilscale(len,nsoil+1) 
    REAL(KIND=r8)   , INTENT(INOUT) :: zmstscale(len,2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: zltrscale(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: zmlscale(len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: drst(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: soilq10(len,nsoil+1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ansqr(len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: soilscaleold(len)
    !     integer, INTENT(IN   ) :: nsib
    !     integer, INTENT(IN   ) :: len
    !     integer, INTENT(IN   ) :: nsoil
    LOGICAL         , INTENT(IN   ) :: forcerestore
    LOGICAL         , INTENT(IN   ) :: dotkef
    REAL(KIND=r8)   , INTENT(IN   ) :: thgeff(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: shgeff(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tke(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ct(len) 
    INTEGER, INTENT(IN   ) :: louis
    !REAL(KIND=r8)   , INTENT(IN   ) :: zwind
    !REAL(KIND=r8)   , INTENT(IN   ) :: ztemp
    REAL(KIND=r8)   , INTENT(INOUT) :: respg(len) 
    REAL(KIND=r8)   , INTENT(IN   ) :: respfactor(len, nsoil+1)
    REAL(KIND=r8)   , INTENT(INOUT) :: pco2ap(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pco2i(len)   ! added for neil suits' programs
    REAL(KIND=r8)   , INTENT(OUT  ) :: pco2c(len)  ! more added nsuits vars
    REAL(KIND=r8)   , INTENT(OUT  ) :: pco2s(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: co2cap(len)  ! moles of air in canopy
    REAL(KIND=r8)   , INTENT(OUT  ) :: cflux(len)





    !    local variables

    INTEGER :: i 
    !      LOGICAL  FRSTVM 
    REAL(KIND=r8)    :: gcon
    REAL(KIND=r8)    :: zln2
    REAL(KIND=r8)    :: ghalf
    REAL(KIND=r8)    :: dmin
    REAL(KIND=r8)    :: pdamp
    REAL(KIND=r8)    :: qdamp
    REAL(KIND=r8)    :: dttin
    REAL(KIND=r8)    :: eps

    REAL(KIND=r8)    :: u2(len)
    REAL(KIND=r8)    :: vsib(len)
    REAL(KIND=r8)    :: thsib(len)
    REAL(KIND=r8)    :: coc(len)
    REAL(KIND=r8)    :: tprcor(len)
    REAL(KIND=r8)    :: cni(len)
    REAL(KIND=r8)    :: cuni(len)
    REAL(KIND=r8)    :: ctni(len)
    REAL(KIND=r8)    :: ctni3(len)
    REAL(KIND=r8)    :: cti(len)
    REAL(KIND=r8)    :: cui(len)
    REAL(KIND=r8)    :: temv(len)
    REAL(KIND=r8)    :: cun(len)
    REAL(KIND=r8)    :: ctn(len) 
    REAL(KIND=r8)    :: z1z0Urt(len)
    REAL(KIND=r8)    :: z1z0Trt(len)
    REAL(KIND=r8)    :: zzwind(len)
    REAL(KIND=r8)    :: zztemp(len) 
    REAL(KIND=r8)    :: epsc(len)
    REAL(KIND=r8)    :: epsg(len) 
    REAL(KIND=r8)    :: USTARo(len) 
    LOGICAL :: jstneu

    !
    INTEGER, PARAMETER :: ITMAX=4
    !  
    !      tbeg = second()
    !      sbeg = timef()

    GCON   = GRAV / 461.5_r8
    EPS    = 1.0_r8 / SNOFAC
    !
    !czzggrst   calculate damping factors
    zln2 = 6.9314718e-1_r8
    ghalf = 1.0257068e1_r8
    dttin = 3.6e3_r8 
    dmin = 6.0e1_r8
    pdamp = EXP (-1.0_r8 * zln2*(dtt*dmin)/(dttin*ghalf))
    qdamp = 1.0_r8 - pdamp


    !czzggrst
    !     FOR OCEAN POINTS, THIS IS THE ONLY CALL TO VMFCAL.  FOR VEGETATED 
    !     POINTS, ITS JUST THE FIRST CALL BEFORE AN ITERATION.
    !
    !      FRSTVM = .TRUE.
    !xx
    IF(louis == 1) THEN
       DO i = 1,len
!          zzwind(i) = z2(i) - d(i) + ZB(i)!zwind
!          zztemp(i) = z2(i) - d(i) + ZB(i)!ztemp
          zzwind(i) =d(i) + ZB(i)!zwind
          zztemp(i) =d(i) + ZB(i)!ztemp

       ENDDO
       CALL VMFCALZ(&
            !PR0      ,&! INTENT(IN) :: pr0 !not used
            !RIBC     ,&! INTENT(IN) :: ribc!not used
            VKRMN    ,&! INTENT(IN) :: vkrmn
            DELTA    ,&! INTENT(IN) :: delta
            GRAV     ,&! INTENT(IN) :: grav
            SH       ,&! INTENT(IN) :: SH    (len)
            !PS       ,&! INTENT(IN) :: PS    (len)    !not used
            z0       ,&! INTENT(IN) :: z0    (len)
            SPDM     ,&! INTENT(IN) :: SPDM  (len)
            SHA      ,&! INTENT(IN) :: SHA   (len)
            ROS      ,&! INTENT(IN) :: ROS   (len)
            CU       ,&! INTENT(OUT) :: CU    (len)
            ct       ,&! INTENT(OUT) :: CT    (len)
            THVGM    ,&! INTENT(OUT) :: THVGM (len)
            RIB      ,&! INTENT(OUT) :: RIB   (len)
            USTAR    ,&! INTENT(OUT) :: USTAR (len)
            VENTMF   ,&! INTENT(OUT) :: VENTMF(len)
            THM      ,&! INTENT(IN) :: THM   (len)
            tha      ,&! INTENT(IN) :: tha   (len) 
            zzwind   ,&! INTENT(IN) :: zzwind (len)
            zztemp   ,&! INTENT(IN) :: zztemp (len)
            cuni     ,&! INTENT(OUT):: cuni  (len)
            cun      ,&! INTENT(OUT):: cun    (len)
            ctn      ,&! INTENT(OUT):: ctn    (len)
            z1z0Urt  ,&! INTENT(OUT):: z1z0Urt(len)
            z1z0Trt  ,&! INTENT(OUT):: z1z0Trt(len)
            len       )! INTENT(IN ) :: len
            !
            !     AERODYNAMIC RESISTANCE
            !
            DO I=1,len
               RA(I)    = ROS(i) / VENTMF(i)
               TEMV(I)  = (Z2(i) - D(i)) / z0(i)
               U2(I)    = SPDM(i) / (CUNI(i) * VKRMN)
               temv(i)  = LOG(temv(i))
               U2(I)    = U2(I) * TEMV(I)
            END DO
    ELSE IF(louis == 2) THEN
       CALL VMFCAL(&
            PR0,&! INTENT(IN) :: pr0
            RIBC     ,&! INTENT(IN) :: ribc
            VKRMN    ,&! INTENT(IN) :: vkrmn
            DELTA    ,&! INTENT(IN) :: delta
            GRAV     ,&! INTENT(IN) :: grav
            SH       ,&! INTENT(IN) :: SH    (len)
            !PS       ,&! INTENT(IN) :: PS    (len) ! not used
            Z0       ,&! INTENT(IN) :: z0    (len)
            SPDM     ,&! INTENT(IN) :: SPDM  (len)
            SHA      ,&! INTENT(IN) :: SHA   (len)
            ZB       ,&! INTENT(IN) :: ZB    (len)
            ROS      ,&! INTENT(IN) :: ROS   (len)
            CU       ,&! INTENT(OUT) :: CU    (len)
            THVGM    ,&! INTENT(OUT) :: THVGM (len)
            RIB      ,&! INTENT(OUT) :: RIB   (len)
            USTAR    ,&! INTENT(OUT) :: USTAR (len)
            VENTMF   ,&! INTENT(OUT) :: VENTMF(len)
            THM      ,&! INTENT(IN):: THM   (len)
            tha      ,&! INTENT(IN):: tha   (len)
            cni      ,&! INTENT(OUT) :: CNI   (len)
            cuni     ,&! INTENT(OUT) :: CUNI  (len)
            ctni     ,&! INTENT(OUT) :: CTNI  (len)
            ctni3    ,&! INTENT(OUT) :: CTNI3 (len)
            cti      ,&! INTENT(OUT) :: CTI   (len)
            cui      ,&! INTENT(OUT) :: CUI   (len)
            len      ,&! INTENT(IN)::len
            thgeff   ,&! INTENT(IN):: thgeff(len)
            shgeff   ,&! INTENT(IN):: shgeff(len)
            tke      ,&! INTENT(IN):: tke   (len)
            ct       ,&! INTENT(OUT)::CT    (len)
            dotkef    )! INTENT(IN):: dotkef 
            !
            !     AERODYNAMIC RESISTANCE
            !
            DO I=1,len
               RA(I)    = ROS(i) / VENTMF(i)
               TEMV(I)  = (Z2(i) - D(i)) / z0(i)
               U2(I)    = SPDM(i) / (CUNI(i) * VKRMN)
               temv(i)  = LOG(temv(i))
               U2(I)    = U2(I) * TEMV(I)
            END DO
    ELSE IF(louis == 3) THEN
       !
       !     the first call to vntlat just gets the neutral values of ustar
       !     and ventmf.
       !
       jstneu=.TRUE.
       CALL vntlax(        &
                   USTARo, &!INTENT(inout) :: ustarn(ncols)
                   !bps   , &!INTENT(in   ) :: bps   (ncols)
                   zb    , &!INTENT(in   ) :: dzm   (ncols)
                   ROS   , &! INTENT(IN        ) :: ROS   (len)
                   THVGM , &!INTENT(in   ) :: dzm   (ncols)
                   VENTMF, &! INTENT(OUT  ) :: VENTMF(len)
                   RIB   , &! INTENT(OUT  ) :: RIB   (len)
                   cu    , &!INTENT(inout) :: cu    (ncols)
                   ct    , &!INTENT(inout) :: ct    (ncols)
                   cti   , &!INTENT(inout) :: cti    (ncols)
                   cui   , &! INTENT(OUT  ) :: CUI   (len)
                   cuni  , &!INTENT(inout) :: cuni  (ncols)
                   ctni  , &!INTENT(inout) :: ctni  (ncols)
                   ustar , &!INTENT(inout) :: ustar (ncols)
                   ra    , &!INTENT(inout) :: ra    (ncols)
                   ta    , &!INTENT(in   ) :: ta    (ncols)
                   u2    , &!INTENT(inout) :: u2    (ncols)
                   THM   , &!INTENT(in   ) :: tm  (ncols)
                   SPDM  , &!INTENT(in   ) :: um  (ncols)
                   d     , &!INTENT(in   ) :: d     (ncols)
                   z0    , &!INTENT(in   ) :: z0    (ncols)
                   z2    , &!INTENT(in   ) :: z2    (ncols)
                   len   , &!INTENT(in   ) :: nmax
                   jstneu, &!INTENT(in   ) :: jstneu
                   len     )!INTENT(in   )
        jstneu=.FALSE.
       CALL vntlax(        &
                   USTARo  , &!INTENT(inout) :: ustarn(ncols)
                   !bps    , &!INTENT(in   ) :: bps   (ncols)
                   zb      , &!INTENT(in   ) :: dzm   (ncols)
                   ROS     , &! INTENT(IN        ) :: ROS   (len)
                   THVGM   , &!INTENT(in   ) :: dzm   (ncols)
                   VENTMF  , &! INTENT(OUT  ) :: VENTMF(len)
                   RIB     , &! INTENT(OUT  ) :: RIB   (len)
                   cu      , &!INTENT(inout) :: cu    (ncols) 
                   ct      , &!INTENT(inout) :: cu    (ncols)
                   cti     , &!INTENT(inout) :: cti    (ncols)
                   cui     , &! INTENT(OUT  ) :: CUI   (len)
                   cuni    , &!INTENT(inout) :: cuni  (ncols)
                   ctni  , &!INTENT(inout) :: ctni  (ncols)
                   ustar , &!INTENT(inout) :: ustar (ncols)
                   ra    , &!INTENT(inout) :: ra    (ncols)
                   ta    , &!INTENT(in   ) :: ta    (ncols)
                   u2    , &!INTENT(inout) :: u2    (ncols)
                   THM    , &!INTENT(in   ) :: tm  (ncols)
                   SPDM    , &!INTENT(in   ) :: um  (ncols)
                   d     , &!INTENT(in   ) :: d     (ncols)
                   z0    , &!INTENT(in   ) :: z0    (ncols)
                   z2    , &!INTENT(in   ) :: z2    (ncols)
                   len   , &!INTENT(in   ) :: nmax
                   jstneu, &!INTENT(in   ) :: jstneu
                   len     )!INTENT(in   )
    ENDIF
    !      FRSTVM = .FALSE.
    !
    !      DO 100 I=1,len
    !         FC(I)=1.0
    !         FG(I)=1.0
    !  100 CONTINUE

    FC(:) = 1.0_r8
    FG(:) = 1.0_r8
    ! 

    CALL RBRD(&
         tc     , &! INTENT(IN        ) :: tc    (len)
         !ts     , &! INTENT(IN        ) :: tm    (len)!isnot used
         rbc    , &! INTENT(IN        ) :: rbc   (len)
         zlt    , &! INTENT(IN        ) :: zlt   (len)
         !tg     , &! INTENT(IN        ) :: tg    (len)!isnot used
         z2     , &! INTENT(IN        ) :: z2    (len)
         u2     , &! INTENT(IN        ) :: u2    (len)
         rd     , &! INTENT(OUT  ) :: rd    (len)
         rb     , &! INTENT(OUT  ) :: rb    (len)
         ta     , &! INTENT(IN        ) :: ta    (len)
         grav   , &! INTENT(IN        ) :: g
         rdc    , &! INTENT(IN        ) :: rdc   (len)
         !ra     , &! INTENT(IN        ) :: ra    (len)!isnot used
         tgs    , &! INTENT(IN        ) :: tgs   (len)
         !wg     , &! INTENT(IN        ) :: wg    (len)!isnot used
         !rsoil  , &! INTENT(IN        ) :: rsoil (len)!isnot used
         !fg     , &! INTENT(IN        ) :: fg    (len)!isnot used
         !hr     , &! INTENT(IN        ) :: hr    (len)!isnot used
         !areas  , &! INTENT(IN        ) :: areas (len)!isnot used
         !respcp , &! INTENT(IN        ) :: respcp(len)!isnot used
         len      )! INTENT(IN   ) :: len
    !


    DO I=1,len 

       !itb...here is inserted some PL prog CAS stuff...
       epsc(i) = 1.0_r8
       epsg(i) = 1.0_r8 
       !itb...pl says " this only makes sense for canopy leaves, since
       !itb...there can only be water OR snow, not both. switching epsc
       !itb...epsc to eps makes the hltm adapt to freezing/fusion.

       IF(snoww(i,1) .GT. 0.0_r8) epsc(i) = eps
       IF(snoww(i,2) .GT. 0.0_r8) epsg(i) = eps

       RC(i) = RST(i) + RB(i) + RB(i)

       RDS(i) = RSOIL(i) * FG(i) + RD(i)

       GECT(i) =  (1.0_r8 - WC(i)) / RC(i)
       GECI(i) = epsc(i) * WC(i) / (RB(i) + RB(i))

       GEGS(i) =  (1.0_r8 - WG(i)) / RDS(i)
       GEGI(i) = epsg(i) * WG(i) / RD(i)

       COC(i) = GECT(i) + GECI(i)


       VSIB(I)  = 1.0_r8/RB(I) + 1.0_r8/RD(I)
       THSIB(I) = (tg(i)/RD(I) + TC(i)/RB(I)) / (VSIB(I)) 
    END DO
    !

    !czzggrst calculate ecmass -- canopy evapotranspiration
    !
    DO i=1,len
       ecmass(i) = (etc(i) - ea(i)) * coc(i) * ros (i) * 0.622e0_r8 /ps(i) * dtt
    ENDDO


    !pl include here a call to respsib
    !pl pass it 

    CALL respsib(len    ,&! INTENT(IN   ) :: len
         !nsib          ,&! INTENT(IN   ) :: nsib!is not used
         nsoil         ,&! INTENT(IN   ) :: nsoil
         wopt          ,&! INTENT(IN   ) :: wopt         (len)
         zm            ,&! INTENT(IN   ) :: zm         (len)
         www           ,&! INTENT(IN   ) :: www         (len,3)
         wsat          ,&! INTENT(IN   ) :: wsat         (len)
         tgs           ,&! INTENT(IN   ) :: tg         (len)
         td            ,&! INTENT(IN   ) :: td         (len,nsoil)
         forcerestore  ,&! INTENT(IN   ) :: forcerestore
         respfactor    ,&! INTENT(IN   ) :: respfactor (len,nsoil+1)
         respg         ,&! INTENT(OUT  ) :: respg         (len)
         soilscale     ,&! INTENT(OUT  ) :: soilscale  (len,nsoil+1)
         zmstscale     ,&! INTENT(OUT  ) :: zmstscale  (len,2)
         soilq10        )! INTENT(OUT  ) :: soilq10         (len,nsoil+1)


    !czzggrst
    !
    !     calculation of canopy conductance and photosynthesis
    !
    CALL phosib(pco2m         ,& ! INTENT(IN   ) :: pco2m         (len)
         pco2ap        ,& ! INTENT(INOUT) :: pco2ap         (len)
         po2m          ,& ! INTENT(IN   ) :: po2m
         vmax0         ,& ! INTENT(IN   ) :: vmax0         (len)
         tice          ,& ! INTENT(IN   ) :: tf
         ps            ,& ! INTENT(IN   ) :: psur         (len)
         green         ,& ! INTENT(IN   ) :: green         (len)
         tran          ,& ! INTENT(IN   ) :: tran         (len,2,2)
         ref           ,& ! INTENT(IN   ) :: ref         (len,2,2)
         gmudmu        ,& ! INTENT(IN   ) :: gmudmu         (len)
         zlt           ,& ! INTENT(IN   ) :: zlt         (len)
         cas_cap_co2   ,& ! INTENT(IN   ) :: cas_cap_co2(len)
         tc            ,& ! INTENT(IN   ) :: tc         (len)
         ta            ,& ! INTENT(IN   ) :: ta         (len) 
         trop          ,& ! INTENT(IN   ) :: trop         (len)
         trda          ,& ! INTENT(IN   ) :: trda         (len)
         trdm          ,& ! INTENT(IN   ) :: trdm         (len)
         slti          ,& ! INTENT(IN   ) :: slti         (len)
         shti          ,& ! INTENT(IN   ) :: shti         (len)
         hltii         ,& ! INTENT(IN   ) :: hltii         (len)
         hhti          ,& ! INTENT(IN   ) :: hhti         (len)
         radn          ,& ! INTENT(IN   ) :: radn         (len,2,2)
         etc           ,& ! INTENT(IN   ) :: etc         (len)
         !etg           ,& ! INTENT(IN   ) :: etgs         (len)!is not used
         !wc            ,& ! INTENT(IN   ) :: wc         (len)!is not used
         ea            ,& ! INTENT(IN   ) :: ea         (len)
         !em            ,& ! INTENT(IN   ) :: em         (len)!is not used
         rb            ,& ! INTENT(IN   ) :: rb         (len)
         ra            ,& ! INTENT(IN   ) :: ra         (len)
         ts            ,& ! INTENT(IN   ) :: tm         (len)
         effcon        ,& ! INTENT(IN   ) :: effcon         (len)
         rstfac        ,& ! INTENT(INOUT) :: rstfac         (len,4)
         binter        ,& ! INTENT(IN   ) :: binter         (len)
         gradm         ,& ! INTENT(IN   ) :: gradm         (len)
         assimn        ,& ! INTENT(OUT  ) :: assimn         (len)
         rst           ,& ! INTENT(IN   ) :: rst         (len)
         atheta        ,& ! INTENT(IN   ) :: atheta         (len)
         btheta        ,& ! INTENT(IN   ) :: btheta         (len)
         tgs           ,& ! INTENT(IN   ) :: tgs         (len)
         respcp        ,& ! INTENT(IN   ) :: respcp         (len)
         aparkk        ,& ! INTENT(INOUT) :: aparkk   (len)
         len           ,& ! INTENT(IN   ) :: len
         !nsib          ,& ! INTENT(IN   ) :: nsib
         omepot        ,& ! INTENT(OUT  ) :: omepot   (len)
         assimpot      ,& ! INTENT(OUT  ) :: assimpot (len)
         assimci       ,& ! INTENT(OUT  ) :: assimci  (len)
         antemp        ,& ! INTENT(OUT  ) :: antemp   (len)
         assimnp       ,& ! INTENT(OUT  ) :: assimnp  (len)
         wsfws         ,& ! INTENT(OUT  ) :: wsfws    (len)
         wsfht         ,& ! INTENT(OUT  ) :: wsfht    (len)
         wsflt         ,& ! INTENT(OUT  ) :: wsflt    (len)
         wci           ,& ! INTENT(OUT  ) :: wci      (len)
         whs           ,& ! INTENT(OUT  ) :: whs      (len)
         wags          ,& ! INTENT(OUT  ) :: wags     (len)
         wegs          ,& ! INTENT(OUT  ) :: wegs     (len)
         aparc         ,& ! INTENT(IN   ) :: aparc    (len)
         pfd           ,& ! INTENT(OUT  ) :: pfd      (len)
         assim         ,& ! INTENT(OUT  ) :: assim    (len)
         td            ,& ! INTENT(IN   ) :: td       (len,nsoil)
         !www           ,& ! INTENT(IN   ) :: www         (len,2)!is not used
         !wopt          ,& ! INTENT(IN   ) :: wopt     (len)!is not used
         !zm            ,& ! INTENT(IN   ) :: zm       (len)!is not used
         !wsat          ,& ! INTENT(IN   ) :: wsat     (len)!is not used
         !tg            ,& ! INTENT(IN   ) :: tg         (len)!is not used
         !soilscale     ,& ! INTENT(IN   ) :: soilscale(len,nsoil+1)!is not used
         !zmstscale     ,& ! INTENT(IN   ) :: zmstscale(len,2)!is not used
         zltrscale     ,& ! INTENT(OUT  ) :: zltrscale(len)
         zmlscale      ,& ! INTENT(OUT  ) :: zmlscale (len)
         drst          ,& ! INTENT(OUT  ) :: drst         (len)
         pdamp         ,& ! INTENT(IN   ) :: pdamp
         qdamp         ,& ! INTENT(IN   ) :: qdamp
         ecmass        ,& ! INTENT(IN   ) :: ecmass         (len)
         dtt           ,& ! INTENT(IN   ) :: dtt
         bintc         ,& ! INTENT(OUT  ) :: bintc    (len)
         tprcor        ,& ! INTENT(OUT  ) :: tprcor   (len)
         !soilq10       ,& ! INTENT(IN   ) :: soilq10  (len,nsoil+1)!is not used
         ansqr         ,& ! INTENT(OUT  ) :: ansqr    (len)
         soilscaleold  ,& ! INTENT(OUT  ) :: soilscaleold(len)
         nsoil         ,& ! INTENT(IN   ) :: nsoil
         !forcerestore  ,& ! INTENT(IN   ) :: forcerestore!is not used
         respg         ,& ! INTENT(IN   ) :: respg         (len)
         pco2c         ,& ! INTENT(OUT  ) :: pco2c    (len) !chloroplast pco2
         pco2i         ,& ! INTENT(OUT  ) :: pco2i         (len)
         pco2s         ,& ! INTENT(OUT  ) :: pco2s    (len)
         co2cap        ,& ! INTENT(OUT  ) :: co2cap   (len)
         cflux          ) ! INTENT(OUT  ) :: cflux         (len)


    !czzggrst block moved up (cxx-140-2000 loop)
    !

    DO i = 1,len
       bintc(i) = bintc(i) * tc(i) / ( 44.6_r8 * tprcor(i)) 
       IF(ea(i).GT.etc(i)) fc(i) = 0.0_r8
       IF(ea(i).GT.etg(i)) fg(i) = 0.0_r8
       HRR(I) = HR(I)
       IF (FG(I) .LT. 0.5_r8) HRR(I) = 1.0_r8
    ENDDO
    !


    RETURN
  END SUBROUTINE VNTLAT



  ! vntlax :performs ventilation mass flux, based on deardorff, mwr, 1972?.


  SUBROUTINE vntlax(&
                  ustarn, &! INTENT(inout) :: ustarn(ncols)
                  !bps   , &! INTENT(in   ) :: bps   (ncols)
                  dzm   , &! INTENT(in   ) :: dzm   (ncols)
                  ROS   , &! INTENT(in   ) :: ros(ncols)!**(JP)** unused
                  THVGM , &! INTENT(out  ) :: thvgm(ncols)!**(JP)** scalar
                  VENTMF, &! INTENT(out  ) ::VENTMF(ncols)!**(JP)** scalar
                  RIB   , &! INTENT(out  ) :: rib(ncols)!**(JP)** scalar
                  cu    , &! INTENT(inout) :: cu    (ncols)
                  ct    , &! INTENT(inout) :: ct    (ncols)
                  cti   , &! INTENT(out  ) :: cti(ncols)!**(JP)** scalar
                  cui   , &! INTENT(out  ) :: cui(ncols)!**(JP)** scalar
                  cuni  , &! INTENT(inout) :: cuni  (ncols)
                  ctni  , &! INTENT(inout) :: ctni  (ncols)
                  ustar , &! INTENT(inout) :: ustar (ncols)
                  ra    , &! INTENT(inout) :: ra    (ncols)
                  ta    , &! INTENT(in   ) :: ta    (ncols)
                  u2    , &! INTENT(inout) :: u2    (ncols)
                  thm   , &! INTENT(in   ) :: thm  (ncols)
                  speedm, &! INTENT(in   ) :: speedm  (ncols)
                  d     , &! INTENT(in   ) :: d     (ncols)
                  z0    , &! INTENT(in   ) :: z0    (ncols)
                  z2    , &! INTENT(in   ) :: z2    (ncols)
                  nmax  , &! INTENT(in   ) :: nmax
                  jstneu, &! INTENT(in   ) :: jstneu
                  ncols   )! INTENT(in   ) :: ncols
    !
    !
    !-----------------------------------------------------------------------
    !         input parameters
    !-----------------------------------------------------------------------
    !
    !   ea..........Pressao de vapor
    !   ta..........Temperatura no nivel de fonte de calor do dossel (K)
    !   um..........Razao entre zonal pseudo-wind (fourier) e seno da
    !               colatitude
    !   vm..........Razao entre meridional pseudo-wind (fourier) e seno da
    !               colatitude
    !   qm..........specific humidity of reference (fourier)
    !   tm..........Temperature of reference (fourier)
    !   dzm  .......Altura media de referencia  para o vento para o calculo
    !               da estabilidade do escoamento
    !   grav........gravity constant      (m/s**2)
    !   cpair.......specific heat of air (j/kg/k)
    !   gasr........gas constant of dry air      (j/kg/k)
    !   bps ........
    !   z2..........height of canopy top
    !   d...........displacement height                        (m)
    !   epsfac......parametro para o gas 0.622
    !
    !
    !
    !-----------------------------------------------------------------------
    !        output parameters
    !-----------------------------------------------------------------------
    !
    !   ustar.........surface friction velocity  (m/s)
    !   ra............Resistencia Aerodinamica (s/m)
    !   u2............wind speed at top of canopy                (m/s)
    !   ventmf........ventilation mass flux
    !-----------------------------------------------------------------------
    !=======================================================================
    !   ncols........Numero de ponto por faixa de latitude
    !   ityp.........Numero do tipo de solo
    !   jstneu.......The first call to vntlat just gets the neutral values
    !                of ustar and ventmf para jstneu=.TRUE..
    !   nmax.........
    !   z0...........roughness length
    !   bps..........bps   (i)=sigki(1)=1.0e0/EXP(akappa*LOG(sig(k)))
    !   cu...........friction  transfer coefficients.
    !   ct...........heat transfer coefficients.
    !   cuni.........neutral friction transfer  coefficients.
    !   ctni.........neutral heat transfer coefficients.
    !=======================================================================
    INTEGER, INTENT(in   ) :: ncols
    LOGICAL, INTENT(in   ) :: jstneu
    INTEGER, INTENT(in   ) :: nmax
    !
    !     vegetation and soil parameters
    !
    REAL(KIND=r8),    INTENT(in   ) :: z2    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: d     (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: z0    (ncols)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(in   ) :: thm   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: speedm(ncols)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout) :: ra    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: u2    (ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    !REAL(KIND=r8),    INTENT(in   ) :: bps   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dzm   (ncols)    
    REAL(KIND=r8),    INTENT(in   ) :: ros(ncols)      !**(JP)** unused

    REAL(KIND=r8),    INTENT(out  ) :: thvgm(ncols)    !**(JP)** scalar
    REAL(KIND=r8),    INTENT(out  ) :: VENTMF(ncols)    !**(JP)** scalar
    REAL(KIND=r8),    INTENT(out  ) :: rib(ncols)      !**(JP)** scalar
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)

    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(out  ) :: cti(ncols)      !**(JP)** scalar
    REAL(KIND=r8),    INTENT(out  ) :: cui(ncols)      !**(JP)** scalar

    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustarn(ncols)


    !REAL(KIND=r8) :: thm(ncols)      !**(JP)** scalar
    !REAL(KIND=r8) :: thvgm(ncols)    !**(JP)** scalar
    !REAL(KIND=r8) :: thvgm
    !REAL(KIND=r8) :: rib
    !REAL(KIND=r8) :: cui
    !REAL(KIND=r8) :: ran(ncols)      !**(JP)** unused
    !REAL(KIND=r8), PARAMETER ::  fsc=66.85_r8
    !REAL(KIND=r8), PARAMETER ::  ftc=0.904_r8
    !REAL(KIND=r8), PARAMETER ::  fvc=0.315_r8

    REAL(KIND=r8), PARAMETER ::  vkrmn=0.40_r8
    REAL(KIND=r8), PARAMETER ::  fsc  =66.85_r8
    REAL(KIND=r8), PARAMETER ::  ftc  =0.904_r8
    REAL(KIND=r8), PARAMETER ::  fvc  =0.915_r8
    REAL(KIND=r8) :: rfac
    REAL(KIND=r8) :: vkrmni
    REAL(KIND=r8) :: g2
    REAL(KIND=r8) :: g3
    REAL(KIND=r8) :: zl
    REAL(KIND=r8) :: xct1
    REAL(KIND=r8) :: xct2
    REAL(KIND=r8) :: xctu1
    REAL(KIND=r8) :: xctu2
    REAL(KIND=r8) :: grib
    REAL(KIND=r8) :: grzl
    REAL(KIND=r8) :: grz2
    REAL(KIND=r8) :: fvv
    REAL(KIND=r8) :: ftt
    REAL(KIND=r8) :: rzl
    REAL(KIND=r8) :: rz2
    INTEGER :: i
    !INTEGER :: ntyp


    rfac  =1.0e2_r8 /gasr

    vkrmni=1.0_r8  /vkrmn
    g2 = 0.75_r8
    g3 = 1.0_r8
    !
    !     cu and ct are the friction and heat transfer coefficients.
    !     cun and ctn are the neutral friction and heat transfer
    !     coefficients.
    !
    IF (jstneu) THEN
       DO i = 1, nmax
          zl = z2(i) + 11.785_r8  * z0(i)
          cuni(i)=LOG((dzm(i)-d(i))/z0(i))*vkrmni
          ustarn(i)=speedm(i)/cuni(i)
          IF (zl < dzm(i)) THEN
             xct1 = LOG((dzm(i)-d(i))/(zl-d(i)))
             xct2 = LOG((zl-d(i))/z0(i))
             xctu1 = xct1
             xctu2 = LOG((zl-d(i))/(z2(i)-d(i)))
             ctni(i) = (xct1 + g2 * xct2) *vkrmni
          ELSE
             xct2 =  LOG((dzm(i)-d(i))/z0(i))
             xctu1 =  0.0_r8
             xctu2 =  LOG((dzm(i)-d(i))/(z2(i)-d(i)))
             ctni(i) = g2 * xct2 *vkrmni
          END IF
          !
          !     neutral values of ustar and ventmf
          !
          u2(i) = speedm(i) - ustarn(i)*vkrmni*(xctu1 + g2*xctu2)
       END DO
       RETURN
    END IF
    !
    !     stability branch based on bulk richardson number.
    !
    DO i = 1, nmax
          !
          !     freelm(i)=.false.
          !
         
          zl = z2(i) + 11.785_r8  * z0(i)
          thvgm(i)   = ta(i)- thm(i)
          rib(i)     =-thvgm(i)   *grav*(dzm(i)-d(i)) &
               /( thm(i)*(speedm(i)-u2(i))**2)
          ! Manzi Suggestion:
          ! rib   (i)=max(-10.0_r8  ,rib(i))
          rib(i)      =MAX(-1.5_r8    ,rib (i)  )
          rib(i)      =MIN( 0.165_r8  ,rib (i)  )
          IF (rib(i)    < 0.0_r8) THEN
             grib = -rib(i)
             grzl = -rib(i)   * (zl-d(i))/(dzm(i)-d(i))
             grz2 = -rib(i)   * z0(i)/(dzm(i)-d(i))
             fvv = fvc*grib
             IF (zl < dzm(i)) THEN
                ftt = (ftc*grib) + (g2-1.0_r8) * (ftc*grzl) - g2 * (ftc*grz2)
             ELSE
                ftt = g2*((ftc*grib) - (ftc*grz2))
             END IF
             cui(i)    = g3*(cuni(i) - fvv)
             cti(i)    = g3*(ctni(i) - ftt)
          ELSE
             rzl = rib(i)   /(dzm(i)-d(i))*(zl-d(i))
             rz2 = rib(i)   /(dzm(i)-d(i))*z0(i)
             fvv = fsc*rib(i)
             IF (zl < dzm(i)) THEN
                ftt = (fsc*rib(i)) + (g2-1) * (fsc*rzl) - g2 * (fsc*rz2)
             ELSE
                ftt = g2 * ((fsc*rib(i)) - (fsc*rz2))
             END IF
             cui(i)    = g3*( cuni(i) + fvv)
             cti(i)    = g3*( ctni(i) + ftt)
          ENDIF
          cu    (i)=1.0_r8/cui(i)
          !**(JP)** ct is not used anywhere else
          ct    (i)=1.0_r8/cti(i)
          !
          !
          !     surface friction velocity and ventilation mass flux
          !
          ustar (i)=MAX(speedm(i)*cu(i),0.0000000001_r8)
          ra(i) = cti(i)    / ustar(i)
          !**(JP)** ran is not used anywhere else
          !ran(i) = ctni(i) / ustarn(i)
          !ran(i) = MAX(ran(i), 0.8_r8 )
          ra(i) = MAX(ra(i), 0.8_r8 )
          VENTMF(i)= ROS(i)*CT(i)* speedm(i)

    END DO
  END SUBROUTINE vntlax


  SUBROUTINE VMFCALZO(&
                    !z0       , &! INTENT(IN   ) :: z0    (len)! is not used   
                    !PR0      , &! INTENT(IN   ) :: pr0 ! is not used
                    !RIBC     , &! INTENT(IN   ) :: ribc! is not used
                    VKRMN    , &! INTENT(IN   ) :: vkrmn
                    !DELTA    , &! INTENT(IN   ) :: delta! is not used
                    GRAV     , &! INTENT(IN   ) :: grav
                    !PS       , &! INTENT(IN   ) :: PS       (len)! is not used   
                    tha      , &! INTENT(IN   ) :: THa      (len)
                    SPDM     , &! INTENT(IN   ) :: SPDM     (len)
                    !ROS      , &! INTENT(IN   ) :: ROS      (len)! is not used   
                    CU       , &! INTENT(OUT  ) :: CU       (len)      
                    THVGM    , &! INTENT(IN   ) :: THVGM    (len)
                    RIB      , &! INTENT(OUT  ) :: RIB      (len) 
                    USTAR    , &! INTENT(OUT  ) :: USTAR    (len)
                    zzwind   , &! INTENT(IN   ) :: zzwind   (len)
                    zztemp   , &! INTENT(IN   ) :: zztemp   (len)
                    len        )! INTENT(IN   ) :: len

    IMPLICIT NONE

    !*****************************************************************************
    !     VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum 
    !     and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
    !     and Boville(1993), and by Beljaars and Holtslag (1991).
    !  
    !     Rerences:
    !       Beljars and Holtslag (1991): Flux parameterization over land surfaces
    !              for atmospheric models. J. Appl. Meteo., 30, 327-341.
    !       Holtslag and Boville (1993): Local versus nonlocal boundary-layer 
    !              diffusion in a global climate model. J. of Climate, 6, 1825-
    !              1842.
    !       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
    !              atmosphere. Boundary-Layer Meteo., 17, 187-202.
    !       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
    !              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
    !              Layer parameterization, ECMWF, 59-79.
    !
    !     General formulation:
    !        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.) 
    !     Transfer coefficients for mommentum and heat fluxes are:
    !        CU = CUN*Fm, and
    !        CT = CTN*Fh
    !        where  CUN and CTN are nutral values of momentum and heat transfers,
    !           and Fm and Fh are stability functions derived from surface
    !           similarity relationships.     
    !*****************************************************************************

    INTEGER, INTENT(IN   ) :: len
    !REAL(KIND=r8), INTENT(IN   ) :: z0    (len)! is not used   
    !REAL(KIND=r8), INTENT(IN   ) :: pr0 ! is not used
    !REAL(KIND=r8), INTENT(IN   ) :: ribc! is not used
    REAL(KIND=r8), INTENT(IN   ) :: vkrmn
    !REAL(KIND=r8), INTENT(IN   ) :: delta! is not used
    REAL(KIND=r8), INTENT(IN   ) :: grav
    !REAL(KIND=r8), INTENT(IN   ) :: PS       (len)! is not used   
    REAL(KIND=r8), INTENT(IN   ) :: THa      (len)
    REAL(KIND=r8), INTENT(IN   ) :: SPDM     (len)
    !REAL(KIND=r8), INTENT(IN   ) :: ROS      (len)! is not used   
    REAL(KIND=r8), INTENT(OUT  ) :: CU       (len)      
    REAL(KIND=r8), INTENT(IN   ) :: THVGM    (len)
    REAL(KIND=r8), INTENT(OUT  ) :: RIB      (len) 
    REAL(KIND=r8), INTENT(OUT  ) :: USTAR    (len)
    REAL(KIND=r8), INTENT(IN   ) :: zzwind   (len)
    REAL(KIND=r8), INTENT(IN   ) :: zztemp   (len)
    !     integer, INTENT(IN   ) :: len

    !     local variables
    !REAL(KIND=r8):: eve      (len)
    !REAL(KIND=r8):: SH       (len)
    REAL(KIND=r8):: CUI      (len)
    REAL(KIND=r8):: TEMV     (len)
    !REAL(KIND=r8):: wgm      (len)
    REAL(KIND=r8):: bunstablM
    REAL(KIND=r8):: cunstablM
    REAL(KIND=r8):: bstabl
    REAL(KIND=r8):: cstabl
    REAL(KIND=r8):: zrib     (len)
    REAL(KIND=r8):: cun      (len)
    REAL(KIND=r8):: z1z0U    (len)
    REAL(KIND=r8):: z1z0Urt  (len)
    REAL(KIND=r8):: fmomn    (len)
    REAL(KIND=r8):: ribtemp
    REAL(KIND=r8):: dm

    INTEGER :: i    
    !      
    !  condtants for surface flux functions, according to Holtslag and
    !      Boville (1993, J. Climate)

    bunstablM = 10.0_r8! constants for unstable function
    cunstablM = 75.0_r8
    bstabl    =  8.0_r8! constants for stable function
    cstabl    = 10.0_r8
    !
    DO I=1,len

       zrib(i) = zzwind(i) **2 / zztemp(i)

       !   Ratio of reference height (zwind/ztemp) and roughness length:
       
       !z1z0U  (i) = zzwind(i)/ z0(i)!0.0002_r8   ! oceanic roughness length
       z1z0U  (i) = zzwind(i)/ 0.0002_r8   ! oceanic roughness length
       
       z1z0Urt(i) = SQRT( z1z0U(i) )
       z1z0U  (i) = LOG( z1z0U(i) )

       !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

       cun(i) = VKRMN*VKRMN / (z1z0U(i)*z1z0U(i) )   !neutral Cm & Ct
       !
       !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
       !   RIB IS THE BULK RICHARDSON NUMBER, between reference height and surface.
       !
       TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
       temv(i) = MAX(0.000001E0_r8,temv(i))
       RIB(I) = -THVGM(I) * GRAV * zrib(i) / TEMV(i) 
    ENDDO

    !   The stability functions for momentum and heat/moisture fluxes as
    !   derived from the surface-similarity theory by Luis (1079, 1982), and
    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag 
    !   (1991).

    DO I=1,len 
       IF(rib(i).GE.0.0_r8) THEN
          !
          !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
          !
          rib(i) = MIN( rib(i), 0.5E0_r8)
          fmomn(i) = (1.0_r8 + cstabl * rib(i) * (1.0_r8+ bstabl * rib(i)))
          fmomn(i) = 1.0_r8 / fmomn(i)
          fmomn(i) = MAX(0.0001E0_r8,fmomn(i))

       ELSE
          !
          !        THE UNSTABLE CASE.
          !            
          ribtemp = ABS(rib(i))
          ribtemp = SQRT( ribtemp )
          dm = 1.0_r8 + cunstablM * cun(i) * z1z0Urt(i) * ribtemp
          fmomn(i) = 1.0_r8 - (bunstablM * rib(i) ) / dm

       END IF
    END DO

    !   surface-air transfer coefficients for momentum CU, for heat and 
    !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

    DO i = 1, len
       CU(i) = CUN(i) * fmomn(i) 
       CUI(i) = 1.0_r8 / CU(i)

       !   Ustar and ventlation mass flux: note that the ustar and ventlation 
       !   are calculated differently from the Deardoff's methods due to their
       !   differences in define the CU and CT.

       USTAR(i) = SPDM(i)*SPDM(i)*CU(i) 
       USTAR(i) = SQRT( USTAR(i) ) 
    ENDDO
    !
    !   Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION  
    !   VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
    !   free convection conditions.
    !  

    RETURN
  END SUBROUTINE VMFCALZO



  SUBROUTINE VMFCALo(&
                  !z0     , &! INTENT(IN   ) :: z0   (len)
                  !PR0     , &! INTENT(IN   ) :: pr0
                  RIBC    , &! INTENT(IN   ) :: ribc
                  VKRMN   , &! INTENT(IN   ) :: vkrmn
                  GRAV    , &! INTENT(IN   ) :: grav
                  SPDM    , &! INTENT(IN   ) :: SPDM   (len)
                  ZB      , &! INTENT(IN   ) :: ZB     (len)
                  !ROS     , &! INTENT(IN   ) :: ROS    (len)
                  CU      , &! INTENT(OUT  ) :: CU     (len)
                  THVGM   , &! INTENT(IN   ) :: THVGM  (len)
                  USTAR   , &! INTENT(OUT  ) :: USTAR  (len)
                  tha     , &! INTENT(IN   ) :: tha    (len)
                  len     , &! INTENT(IN   ) :: len
                  thgeff  , &! INTENT(IN   ) :: thgeff (len)
                  tke     , &! INTENT(IN   ) :: tke    (len)
                  dotkef    )! INTENT(IN   ) :: dotkef
    IMPLICIT NONE
    ! subroutine is stripped from vmfcal.F and is used only to calculate
    !   ustart and cu for oceanic values of the surface roughness length
    !
    !***  VENTILATION MASS FLUX, BASED ON DEARDORFF, MWR, 1972
    !
    INTEGER      , INTENT(IN   ) :: len
    !REAL(KIND=r8), INTENT(IN   ) :: pr0
    REAL(KIND=r8), INTENT(IN   ) :: ribc
    REAL(KIND=r8), INTENT(IN   ) :: vkrmn
    REAL(KIND=r8), INTENT(IN   ) :: grav
    !REAL(KIND=r8), INTENT(IN   ) :: z0   (len)
    REAL(KIND=r8), INTENT(IN   ) :: SPDM   (len)
    REAL(KIND=r8), INTENT(IN   ) :: ZB     (len)
    !REAL(KIND=r8), INTENT(IN   ) :: ROS    (len)
    REAL(KIND=r8), INTENT(OUT  ) :: CU     (len)
    REAL(KIND=r8), INTENT(IN   ) :: THVGM  (len)
    REAL(KIND=r8), INTENT(OUT  ) :: USTAR  (len)
    REAL(KIND=r8), INTENT(IN   ) :: tha    (len)
    !     integer,INTENT(IN   ) :: len
    REAL(KIND=r8), INTENT(IN   ) :: thgeff (len)
    REAL(KIND=r8), INTENT(IN   ) :: tke    (len)
    LOGICAL       , INTENT(IN   ) :: dotkef

    !     local variables
    REAL(KIND=r8):: RIB    (len) 
    REAL(KIND=r8):: ribmax
    REAL(KIND=r8):: vkrinv
    REAL(KIND=r8):: TEMV   (len)
    REAL(KIND=r8):: ZDRDRF (len)
    REAL(KIND=r8):: CHOKE  (len)
    REAL(KIND=r8):: sqrtke
    REAL(KIND=r8):: cui    (len)
    REAL(KIND=r8):: cuni   (len)
    REAL(KIND=r8):: cni    (len)
    INTEGER :: i
    !
    RIBMAX = 0.9_r8*RIBC
    VKRINV = 1.0_r8/VKRMN
    !
    !        CUNI AND CTN1 ARE INVERSES OF THE NEUTRAL TRANSFER COEFFICIENTS
    !        DEARDORFF EQS(33) AND (34).
    !        PR0 IS THE TURBULENT PRANDTL NUMBER AT NEUTRAL STABILITY.
    !
    DO I=1,len
       !CNI(i) = 0.025_r8*ZB(i)/z0(i)!0.0002_r8    ! oceanic surface roughness length
       CNI(i) = 0.025_r8*ZB(i)/0.0002_r8    ! oceanic surface roughness length

       cni(i) = LOG(cni(i))
       CNI(i)  = CNI(i)  * VKRINV
       CUNI(i) = CNI(i)  + 8.4_r8
       CNI(i)  = GRAV    * ZB(i)
       !
       !        SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
       !        RIB IS THE BULK RICHARDSON NUMBER, DEARDORFF EQ (25)
       !
    END DO

    IF (dotkef) THEN
       DO I=1,len
          rib(i) = -thvgm(i)*grav*zb(i)/(thgeff(i)*tke(i)) 
       ENDDO
       DO I=1,len 
          IF(rib(i).GE.0.0_r8) THEN
             cu(i) = ((-0.0307_r8*rib(i)**2)/(61.5_r8+rib(i)**2))+0.044_r8
          ELSE
             cu(i) = ((-0.016_r8*rib(i)**2)/(4.2e4_r8+rib(i)**2))+0.044_r8
          ENDIF
       ENDDO
       DO I=1,len 
          sqrtke = SQRT(tke(i))
          USTAR(I) = SQRT(sqrtke*spdm(I)*CU(I))
       ENDDO
    ELSE
       DO I=1,len
          TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
          temv(i) = MAX(0.000001E0_r8,temv(i))
          RIB(I) = -THVGM(I) * CNI(I) / TEMV(i) 
       END DO

       DO I=1,len 
          IF(rib(i).GE.0.0_r8) THEN
             !
             !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
             !
             CHOKE(i)=1.0_r8-MIN(RIB(I),RIBMAX)/RIBC
             CU(I)=CHOKE(i)/CUNI(I)
          ELSE
             !
             !        FIRST, THE UNSTABLE CASE. DEARDORFF EQS(28), (29), AND (30)
             !
             ZDRDRF(i) = LOG10(-RIB(I)) - 3.5_r8
             CUI(I)    = CUNI(I) - 25.0_r8 * EXP(0.26_r8 * ZDRDRF(i)  &
                  - 0.03_r8 * ZDRDRF(i) * ZDRDRF(i))
             CU(I)     = 1.0_r8 / MAX(CUI(I),0.5_r8 * CUNI(I))
          END IF
       END DO

       DO I=1,len
          USTAR(i) =SPDM(i)*CU(i)
       END DO
    ENDIF
    !   
    RETURN
  END SUBROUTINE VMFCALo


  !
  !===================SUBROUTINE DELRF=====================================
  !

  SUBROUTINE DELLWF(&
                  !DTA    , &! INTENT(IN   ) :: dta ! is not used
                  dtc4   , &! INTENT(IN   ) :: dtc4 (len)
                  dtg4   , &! INTENT(IN   ) :: dtg4 (len)
                  dts4   , &! INTENT(IN   ) :: dts4 (len)
                  fac1   , &! INTENT(IN   ) :: fac1 (len)
                  areas  , &! INTENT(IN   ) :: areas(len)
                  lcdtc  , &! INTENT(OUT  ) :: lcdtc(len)
                  lcdtg  , &! INTENT(OUT  ) :: lcdtg(len)
                  lcdts  , &! INTENT(OUT  ) :: lcdts(len)
                  lgdtg  , &! INTENT(OUT  ) :: lgdtg(len)
                  lgdtc  , &! INTENT(OUT  ) :: lgdtc(len)
                  lsdts  , &! INTENT(OUT  ) :: lsdts(len)
                  lsdtc  , &! INTENT(OUT  ) :: lsdtc(len)
                  len      )! INTENT(IN   ) :: len

    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground radiative
    !        heat fluxes with respect to Tc, Tgs
    !pl   Here we are doing only the long wave radiative loss, which is the
    !pl   only radiative quantity we are trying to bring to the next time step.
    !
    !======================================================================== 

    !------------------------------INPUT is coming from Netrad-------------
    !
    !       dtc4, dtg4, dts4, which are the derivatives of the LW loss
    !
    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       LCDTC          dLC/dTC 
    !       LCDTG          dLC/dTG
    !       LCDTS          dLC/dTS
    !       LGDTG          dLG/dTG
    !       LGDTC          dLG/dTC
    !       LSDTS          dLS/dTS
    !       LSDTC          dLS/dTC
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE
    !

    INTEGER, INTENT(IN   ) :: len

    !REAL(KIND=r8)   , INTENT(IN   ) :: dta ! is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: dtc4 (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: dtg4 (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: dts4 (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: fac1 (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: areas(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lcdtc(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lcdtg(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lcdts(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lgdtg(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lgdtc(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lsdts(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: lsdtc(len)
    !   INTEGER, INTENT(IN   ) :: len



    !     local variables

    INTEGER :: i

    DO I=1,len  

       !pl canopy leaves:
       LCDTC(I) =   2 * DTC4(i) * fac1(i)
       LCDTG(I) =     - DTG4(i) * fac1(i) * (1.0_r8-areas(i))
       LCDTS(I) =     - DTS4(i) * fac1(i) * (   areas(i))

       !pl ground:
       LGDTG(I) =   DTG4(i)
       LGDTC(I) = - DTC4(i) * fac1(i)

       !pl snow:
       LSDTS(I) =   DTS4(i)
       LSDTC(I) = - DTC4(i) * fac1(i)
       !
    ENDDO

    RETURN
  END SUBROUTINE DELLWF


  !
  !===================SUBROUTINE DELHF25=====================================
  !
  !pl ATTENTION receiving tgs instead of tg !!!
  SUBROUTINE DELHF(&
       DTA,&
       CP,&
      ! bps,&
       tm,&
       tg,&
       !ts,&
       tc,&
       ta,&
       ros,&
       ra,&
       rb,&
       rd  ,& 
       HCDTC,&
       HCDTA,&
       HGDTG,&
       HGDTA,&
       HSDTS,&
       HSDTA,&
       HADTA,&
       HADTH,&
       hc, &
       hg, &
       hs, &
       fss, &
       len ) 

    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground sensible
    !        heat fluxes with respect to Tc, Tgs, and Theta-m.
    !     Calculation of initial sensible heat fluxes.
    !
    !========================================================================


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       HC             CANOPY SENSIBLE HEAT FLUX (J M-2)
    !       HG             GROUND SENSIBLE HEAT FLUX (J M-2)
    !       HS             SNOW   SENSIBLE HEAT FLUX (J M-2)
    !       HA             CAS    SENSIBLE HEAT FLUX (J M-2)
    !       HCDTC          dHC/dTC 
    !       HCDTA          dHC/dTA
    !       HGDTG          dHG/dTG
    !       HGDTA          dHG/dTA
    !       HSDTS          dHS/dTS
    !       HSDTA          dHS/dTA
    !       HADTA          dHA/dTA
    !       HADTH          dHA/dTH
    !       AAC            dH/dTC
    !       AAG            dH/dTGS
    !       AAM            dH/dTH
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: len
    REAL(KIND=r8)   , INTENT(IN   ) :: dta
    REAL(KIND=r8)   , INTENT(IN   ) :: cp
    !REAL(KIND=r8)   , INTENT(IN   ) :: bps   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tm    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg    (len)
    !REAL(KIND=r8)   , INTENT(IN   ) :: ts    (len) !is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: tc    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ta    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ros   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ra    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rb    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rd    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HCDTC (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HCDTA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HGDTG (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HGDTA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HSDTS (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HSDTA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HADTA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: HADTH (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hc    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hg    (len) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: fss   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hs    (len)
    !     integer, INTENT(IN   ) :: len

    !     local variables

    REAL(KIND=r8)    :: rai
    REAL(KIND=r8)    :: rbi
    REAL(KIND=r8)    :: rdi
    INTEGER :: i

    !-----------------------------------------------------------------------
    !
    !     FLUXES EXPRESSED IN JOULES M-2, although in SIBSLV WE THEN WANT W/m2
    !pl WHY ????
    !
    !pl if we were to keep things simple, there is no need to separate
    !pl HG and HS, but it helps the derivatives keep clean.
    !
    !      HC          (HC)    : EQUATION (63) , SE-86
    !      HG          (HG)    : EQUATION (65) , SE-86
    !      HS          (HS)    : EQUATION (65) , SE-86
    !      HA          (HA)    : EQUATION ???
    !-----------------------------------------------------------------------

    DO I=1,len  
       rai = 1.0_r8 / ra(I)
       rbi = 1.0_r8 / rb(I)
       rdi = 1.0_r8 / rd(I)
       !        D1 = rai + rbi + rdi
       !        d1i = 1.0_r8 / d1

       !pl these are the current time step fluxes in J/m2
       !pl can we change this to W/m2 ???

       HC(I)   = CP * ros(i) * (tc(i) - ta(i)) * rbi * DTA
       HG(I)   = CP * ros(i) * (tg(I) - ta(i)) * rdi * DTA
       HS(I)   = CP * ros(i) * (tg(I) - ta(i)) * rdi * DTA
       fss(I)  = CP * ros(i) * (ta(I) - tm(i)) * rai * DTA

       !pl now we do the partial derivatives

       !pl these are done assuming the fluxes in W/m2

       !pl for canopy leaves sensible heat flux: W/(m2 * K) 
       HCDTC(I) =   CP * ros(i) * rbi
       HCDTA(I) = - HCDTC(I)
       !pl for ground and snow sensible heat fluxes: W/(m2 * K)
       HGDTG(I) =   CP * ros(i) * rdi  
       HSDTS(I) =   HGDTG(I)  
       HGDTA(I) = - HGDTG(I)  
       HSDTA(I) = - HGDTG(I)  
       !pl for the canopy air space (CAS) sensible heat flux: W/(m2 * K)
       HADTA(I) =   CP * ros(i) * rai
       HADTH(I) = - HADTA(I)!/bps(i)

       !pl ATTENTION !!!! DANGER !!!!! THIS WILL NOT WORK WITHOUT sibdrv = true
       !pl for mixed layer (ref temp if not sibdrv): YET TO BE DONE
       !        AAG(I) = rdi * d1i
       !        AAC(I) = rbi * d1i
       !        AAM(I) = rai * d1i * bps(i)
       !
    ENDDO

    !itb...
    !      print*,'delhf: rbi,rdi,rai,dta=',rbi,rdi,rai,dta
    !      print*,'delhf: cp, ros=',cp,ros
    !      print*,'delhf: tm,tc,tg,ta=',tm,tc,tg,ta
    !      print*,'delhf: hc,hg,hs,fss=',hc,hg,hs,fss

    RETURN
  END SUBROUTINE DELHF


  !
  !==================SUBROUTINE DELEF25======================================
  !
  SUBROUTINE DELEF(&
       DTA,&
       CP,&
       !ps,&
       em,&
       ea,&
       ros,&
       HRr,&
       !fc,&
       fg,&
       ra,&
!       rb,&
       rd,&
       !rc,&
       rsoil,&
!       snow,&
!       capac,&
!       www,&
       ECDTC,&
       ECDEA,&
       EGDTG,&
       EGDEA,&
       ESDTS,&
       ESDEA,&
       EADEA,&
       EADEM,&
       ec,&
       eg,&
       es,&
       fws,&
       !hltm,&
       etc,&
       etg,&
       btc,&
       btg,&
!       bts, &
!       areas, &
       gect,&
       geci,&
       gegs,&
       gegi, &
       psy, &
       snofac, &
       hr,    &
       len)    
    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground latent
    !        heat fluxes with respect to Tc, Tgs, Theta-m, and Qm.
    !     Calculation of initial latent heat fluxes.

    !pl the ETC, ETG and so on are the vapor pressure at temps TC, TG and so on
    !pl the BTC, BTG are the derivatives of ETC, ETG with relation to TC, TG etc.
    !
    !======================================================================== 

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       EC             ECT + ECI
    !       EG             EGS + EGI
    !       ECDTC          dEC/dTC
    !       ECDTG          dEC/dTGS
    !       ECDQM          dEC/dQM
    !       EGDTC          dEG/dTC
    !       EGDTG          dEG/dTGS
    !       EGDQM          dEG/dQM
    !       BBC            dE/dTC
    !       BBG            dE/dTGS
    !       BBM            dE/dQM
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE
    !         
    INTEGER, INTENT(IN   ) :: len

    REAL(KIND=r8)   , INTENT(IN   ) :: dta
    REAL(KIND=r8)   , INTENT(IN   ) :: cp
!    REAL(KIND=r8)   , INTENT(IN   ) :: ps    (len)! is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: em    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ea    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ros   (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hrr   (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: fc    (len)! is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: fg    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ra    (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: rb    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rd    (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: rc    (len)! is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: rsoil (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: snow  (len,2)! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: capac (len,2)! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: www   (len,3) ! is not used
    REAL(KIND=r8)   , INTENT(OUT  ) :: ECDTC (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ECDEA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: EGDTG (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: EGDEA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ESDTS (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ESDEA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: EADEA (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: EADEM (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ec    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: eg    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: es    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fws   (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: hltm         ! is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: etc   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: etg   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: btc   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: btg   (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: bts   (len) ! is not used   
!    REAL(KIND=r8)   , INTENT(IN   ) :: areas (len)! is not used   
    REAL(KIND=r8)   , INTENT(IN   ) :: gect  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: geci  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: gegs  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: gegi  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: psy   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: snofac
    REAL(KIND=r8)   , INTENT(IN   ) :: hr    (len)
    !     integer, INTENT(IN   ) :: len


!    REAL(KIND=r8)    :: rst   (len)
    REAL(KIND=r8)    :: rds   (len)
!    REAL(KIND=r8)    :: sh    (len)
    REAL(KIND=r8)    :: COG1  (len)
    REAL(KIND=r8)    :: COC   (len)
 !   REAL(KIND=r8)    :: D2    (len) 
!    REAL(KIND=r8)    :: tc_0  (len)
!    REAL(KIND=r8)    :: tgs_0 (len)
!    REAL(KIND=r8)    :: intg  (len)
!    REAL(KIND=r8)    :: intc  (len)
    REAL(KIND=r8)    :: COG2  (len)
    REAL(KIND=r8)    :: rcp   (len)
    REAL(KIND=r8)    :: rcpg  (len)    
!    REAL(KIND=r8)    :: limci (len)
!    REAL(KIND=r8)    :: limgi (len)
!    REAL(KIND=r8)    :: limct (len)
!    REAL(KIND=r8)    :: limgs (len)


!    REAL(KIND=r8)    :: cogr   (len)
!    REAL(KIND=r8)    :: cogs   (len)
    REAL(KIND=r8)    :: cpdpsy (len)
    REAL(KIND=r8)    :: resrat

    INTEGER :: i
    !                                                                       
    !     MODIFICATION FOR SOIL DRYNESS : HR=REL. HUMIDITY IN TOP LAYER     
    !                                                                       
    resrat = 0.5_r8

    DO I=1,len                                                   
       hrr(I) = HR(I)               
       IF(fg(i).LT.0.5_r8) hrr(i) = 1.0_r8
    ENDDO

    DO I=1,len
       rcp(i)  = ros(i) * cp   
       rcpg(i) = rcp(i)/psy(i)              ! this is rho * cp / gamma
       cpdpsy(i) = cp / psy(i)
       rds(i) = rsoil(i) + rd(i)

       !-----------------------------------------------------------------------
       !                                                                       
       !     CALCULATION OF SURFACE RESISTANCE COMPONENTS, SEE EQUATIONS (64,66)
       !       OF SE-86                                   
       !pl the ge?? coefficients come all the way from VNTLAT and are common to
       !pl all subroutines:
       !pl     gect(i)  =      (1.0_r8 -wc(i)) /  rc(i)
       !pl        geci(i)  = epsc(i) * wc(i)  / (RB(I) + RB(I))
       !pl        gegs(i)  =        (1-wg(i)) / (rds(i))
       !pl        gegi(i)  = epsg(i) * wg(i)  /  rd(i)
       !                                                                       
       !-----------------------------------------------------------------------

       COC (I) =  gect(i) + geci(i)
       COG1(i) = (gegi(i) + gegs(i)*HRR(i))
       COG2(i) = (gegi(i) + gegs(i)       )

       !            D2(I)   = 1.0_r8 / RA(I) + COC(I) + COG2(I)
       !-----------------------------------------------------------------------
       !                                                                       
       !     FLUXES EXPRESSED IN JOULES M-2   CPL WHY ?????
       !                                                                       
       !      ec         (EC)    : EQUATION (64) , SE-86
       !      eg         (EG)    : EQUATION (66) , SE-86
       !      es         (ES)    : EQUATION (66) , SE-86
       !      ea         (EA)    : EQUATION ????
       !-----------------------------------------------------------------------

       !pl these are the current time step fluxes in J/m2  WHY ?????

       !pl notice that the fluxes are already limited by the altered e*(T) values

       ec(I)  = ( etc(I) - ea(i)) * COC(I) * ros(i)  * dta * cpdpsy(i)

       eg(I)  = ( etg(I) * COG1(I) - ea(i) * COG2(I)) * ros(i) * dta * cpdpsy(i)

       es(I)  = ((etg(I) - ea(i))/rd(i) )* ros(i) * dta * cpdpsy(i)/snofac
       fws(I) = ((ea(I)  - em(i) ) / ra(i)) * ros(i) * dta * cpdpsy(i)

       !pl now we do the partial derivatives  these assume W/m2

       !pl for the canopy leaves vapor pressure: W/ (m2* K)
       ECDTC(I) =    btc(I) * COC(I) * ros(i) * CPDPSY(i)   
       ECDEA(I) = - COC(I) * ros(i) * CPDPSY(i)

       !pl for ground latent heat fluxes: W/ (m2* K)
       EGDTG(I) =   btg(I) * COG1(I) * ros(i) * CPDPSY(i)
       EGDEA(I) = - (COG2(I)) * ros(i) * CPDPSY(i)               

       !pl for snow latent heat fluxes: W/ (m2* K)
       ESDTS(I) =   btg(I) * ros(i) * CPDPSY(i)/RD(i)

       !pl for snow latent heat fluxes: W/ (m2 * Pa)
       ESDEA(I) = - ros(i) * CPDPSY(i)/RD(i)              

       !pl for CAS latent heat fluxes: W/ (m2* Pa)
       EADEA(I) = ros(i) * CPDPSY(i) / ra(i)              
       EADEM(I) = - EADEA(I)            

       !PL ATTENTION !!!! DANGER !!! do not use without sibdrv = true
       !pl these all need to be re-done for the GCM (no sibdrv)
       !-----------------------------------------------------------------------
       !      BBC       (dE/dTC)  : EQUATION (13) , SA-89B
       !      BBG       (dE/dTGS) : EQUATION (13) , SA-89B
       !      BBM       (dE/dQM)  : EQUATION (13) , SA-89B
       !-----------------------------------------------------------------------
       !        BBG(I) = (COG1(I) / D2(i))
       !     *          * btg(I) * 0.622_r8 * ps(i)           
       !     *       / ((ps(i) - etg(I)) * (ps(i) - etg(I)))                
       !        BBC(I) = (COC(I)  / D2(i))
       !     *            * btc(I) * 0.622_r8 * ps(i)           
       !     *       / ((ps(i) - etc(I)) * (ps(i) - etc(I)))                
       !        BBM(I) = 1.0_r8   / (ra(I)  * D2(i))                             
       !                                                                       
    ENDDO
    !           

    !      print*,'delef: em,ea,etg,etc=',em,ea,etg,etc
    !      print*,'delef: coc,cog1,cog2,dta=',coc,cog1,cog2,dta
    !      print*,'delef: ros, cpdpsy=',ros,cpdpsy
    !      print*,'delef: ec,eg,es,fws=',ec,eg,es,fws
    RETURN
  END SUBROUTINE DELEF


  SUBROUTINE soilprop( td, tg, slamda, shcap, www, poros, ztdep, &
       asnow, snoww, areas, tf, snomel, sandfrac,&
       clayfrac, &
       !nsib, &
       len, &
       nsoil )

    IMPLICIT NONE

    !     this subroutine calculates the soil thermal properties heat capacity
    !         (shcap) and conductivity (slamda). Phase changes are incorporated
    !         into the heat capacity over a range of -0.5C to 0.5C.
    !     slamda(n) is the conductivity between layers n and n+1 / delta z
    !     layer 1 is the deepest layer, layer nsoil is immediately below the 
    !         surface
    !     treatment based on Bonan

    !     argument list variables
    INTEGER, INTENT(IN   ) :: len
    !INTEGER, INTENT(IN   ) :: nsib
    INTEGER, INTENT(IN   ) :: nsoil
    REAL(KIND=r8)   , INTENT(IN   ) :: td       (len,nsoil)
    REAL(KIND=r8)   , INTENT(IN   ) :: tg       (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: slamda   (len,nsoil)
    REAL(KIND=r8)   , INTENT(OUT  ) :: shcap    (len,nsoil)
    REAL(KIND=r8)   , INTENT(IN   ) :: www      (len,3)
    REAL(KIND=r8)   , INTENT(IN   ) :: poros    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ztdep    (len,nsoil)
    REAL(KIND=r8)   , INTENT(IN   ) :: asnow
    REAL(KIND=r8)   , INTENT(IN   ) :: snoww    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: areas    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: tf
    REAL(KIND=r8)   , INTENT(IN   ) :: snomel
    REAL(KIND=r8)   , INTENT(IN   ) :: sandfrac (len)!nor used
    REAL(KIND=r8)   , INTENT(IN   ) :: clayfrac(len)!nor used
    !     integer, INTENT(IN   ) :: nsib
    !     integer, INTENT(IN   ) :: len
    !     integer, INTENT(IN   ) :: nsoil





    !     local variables

    REAL(KIND=r8)    :: shcapu
    REAL(KIND=r8)    :: shcapf
    REAL(KIND=r8)    :: shsnow
    REAL(KIND=r8)    :: klamda(len,nsoil)
    REAL(KIND=r8)    :: zsnow(len)
    REAL(KIND=r8)    :: kf
    REAL(KIND=r8)    :: ku
    REAL(KIND=r8)    :: ksnow
    REAL(KIND=r8)    :: ksoil(len)
    REAL(KIND=r8)    :: siltfrac (len)
!    REAL(KIND=r8)    :: consoi (len)

    REAL(KIND=r8)    :: kiw(len,3)
    REAL(KIND=r8)    :: kii(len,3)
    REAL(KIND=r8)    :: shcap1
    REAL(KIND=r8)    :: klamda1
    REAL(KIND=r8)    :: bd
    REAL(KIND=r8)    :: csol
    REAL(KIND=r8)    :: klamda2
    INTEGER :: i
    INTEGER :: n

    !     snow density assumed 0.25 that of water to convert from Bonan's constants
    !        for ksnow and max effective snow depth
    ksnow = 0.085_r8   
    DO i = 1,len
       siltfrac(i) = 0.01_r8 * max((100.0_r8 - 100.0_r8*sandfrac(i) - 100.0_r8*clayfrac(i)),0.0_r8)
       !powliq = poros(i) * (www(i,1)+www(i,2)+www(i,3))/3.0_r8
       !zcondry = sandfrac(i) * 0.300_r8 +  &
       !          siltfrac(i) * 0.265_r8 +  &
       !          clayfrac(i) * 0.250_r8 ! +
       !ksoil(i) = zcondry !* ((0.56_r8*100.0_r8)**powliq)

!
       !         ksoil(i) = (8.8_r8 * sandfrac(i) + 
       !     *                    2.92_r8 * (1.0_r8-sandfrac(i)) )
       !     *                           **(1.0_r8-poros(i))
       bd    = (1.0_r8-poros(i))*2.7e3_r8
       ksoil(i) = (0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd)! (W/m/Kelvin) (nlevsoi)

       !ksoil(i) = 6.0_r8**(1.0_r8-poros(i))
    ENDDO
    DO n = 1,3
       DO i = 1,len
          kiw(i,n) = 0.6_r8**(www(i,n)*poros(i))
          kii(i,n) = 2.2_r8**(www(i,n)*poros(i))
       ENDDO
    ENDDO
    !     heat capacity calculation and conductivity calculation

    !     soil layers 1 and 2 ( deepest layers, use www(3) )
    DO n = 1,2
       DO i = 1,len
          csol    = (2.128_r8*sandfrac(i)+2.385_r8*clayfrac(i)) / (sandfrac(i)+clayfrac(i))*1.e6_r8  ! J/(m3 K)
          SHCAPu  = ( 0.5_r8 * (1.0_r8 - POROS(i)) + WWW(i,3) * POROS(i)) * 4.186E6_r8
          shcapf  =    0.5_r8 * (1.0_r8 + poros(i) * (www(i,3) - 1.0_r8)) * 4.186E6_r8
          SHCAPu  = 0.2_r8*SHCAPu + 0.8_r8*csol
          shcapf  = 0.2_r8*shcapf + 0.8_r8*csol
          
          ku = (ksoil(i)*kiw(i,3)-0.15_r8)*www(i,3) + 0.15_r8 
          kf = (ksoil(i)*kii(i,3)-0.15_r8)*www(i,3) + 0.15_r8
          IF(td(i,n).GE.tf+0.5_r8) THEN
             shcap(i,n) = shcapu * ztdep(i,n)
             klamda(i,n) = ku
          ELSE IF(td(i,n).LE.tf-0.5_r8) THEN
             shcap(i,n) = shcapf * ztdep(i,n)
             klamda(i,n) = kf
          ELSE
             shcap(i,n) = (0.5_r8*(shcapu+shcapf) + &
                  poros(i)*www(i,3)*snomel) * ztdep(i,n)
             klamda(i,n) = kf + (ku-kf)*(td(i,n)+0.5_r8-tf)
          ENDIF
       ENDDO
    ENDDO

    !     soil layers 3,4,5 and nsoil ( intermediate layers, use www(2) )
    DO n = 3,nsoil
       DO i = 1,len
          csol    = (2.128_r8*sandfrac(i)+2.385_r8*clayfrac(i)) / (sandfrac(i)+clayfrac(i))*1.e6_r8  ! J/(m3 K)
          SHCAPu  = ( 0.5_r8*(1.0_r8-POROS(i)) + WWW(i,2)*POROS(i)) * 4.186E6_r8
          shcapf =  0.5_r8 * (1.0_r8 + poros(i) * (www(i,2)-1.0_r8)) * 4.186E6_r8
          SHCAPu  = 0.2_r8*SHCAPu + 0.8_r8*csol
          shcapf  = 0.2_r8*shcapf + 0.8_r8*csol
          ku = (ksoil(i)*kiw(i,2)-0.15_r8)*www(i,2) + 0.15_r8 
          kf = (ksoil(i)*kii(i,2)-0.15_r8)*www(i,2) + 0.15_r8
          IF(td(i,n).GE.tf+0.5_r8) THEN
             shcap(i,n) = shcapu * ztdep(i,n)
             klamda(i,n) = ku
          ELSE IF(td(i,n).LE.tf-0.5_r8) THEN
             shcap(i,n) = shcapf * ztdep(i,n)
             klamda(i,n) = kf
          ELSE
             shcap(i,n) = (0.5_r8*(shcapu+shcapf) + &
                  poros(i)*www(i,2)*snomel) * ztdep(i,n)
             klamda(i,n) = kf + (ku-kf)*(td(i,n)+0.5_r8-tf)
          ENDIF
       ENDDO
    ENDDO

    !     soil layer nsoil ( top layer, use www(1) )
    !     if snow covered add additional heat capacity due to snow
    DO i = 1,len
       csol    = (2.128_r8*sandfrac(i)+2.385_r8*clayfrac(i)) / (sandfrac(i)+clayfrac(i))*1.e6_r8  ! J/(m3 K)
       SHCAPu  = ( 0.5_r8 * (1.0_r8 - POROS(i)) + WWW(i,1) * POROS(i)) * 4.186E6_r8
       shcapf  =   0.5_r8 * (1.0_r8 + poros(i) * (www(i,1) - 1.0_r8) ) * 4.186E6_r8
       SHCAPu  = 0.2_r8*SHCAPu + 0.8_r8*csol
       shcapf  = 0.2_r8*shcapf + 0.8_r8*csol
       ku = (ksoil(i)*kiw(i,1)-0.15_r8)*www(i,1) + 0.15_r8 
       kf = (ksoil(i)*kii(i,1)-0.15_r8)*www(i,1) + 0.15_r8
       IF(td(i,nsoil).GE.tf+0.5_r8) THEN
          shcap1 = shcapu 
          klamda1 = ku
       ELSE IF(td(i,nsoil).LE.tf-0.5_r8) THEN
          shcap1 = shcapf
          klamda1 = kf
       ELSE
          shcap1 = 0.5_r8*(shcapu+shcapf) + poros(i)*www(i,1)*snomel
          klamda1 = kf + (ku-kf)*(td(i,nsoil)+0.5_r8-tf)
       ENDIF
       IF(tg(i).GE.tf+0.5_r8) THEN
          klamda2 = ku
       ELSE IF(tg(i).LE.tf-0.5_r8) THEN
          klamda2 = kf
       ELSE
          klamda2 = kf + (ku-kf)*(tg(i)+0.5_r8-tf)
       ENDIF
       zsnow(i) = MIN( MAX(1.0_r8/asnow,snoww(i)) - 0.02_r8, 0.25E0_r8 )
       shsnow = (0.5_r8 * 4.186e6_r8 * zsnow(i) + shcap1 * 0.02_r8) / (zsnow(i)+0.02_r8)
       shcap(i,nsoil) =  shcap(i,nsoil) * (1.0_r8-areas(i)) + areas(i) *   &
            shcap(i,nsoil)*shsnow*(ztdep(i,nsoil)+zsnow(i)+0.02_r8) / &
            (shsnow*ztdep(i,nsoil)+shcap(i,nsoil) * &
            (zsnow(i)+0.02_r8)) * (ztdep(i,nsoil)+zsnow(i)+0.02_r8)
       klamda1 = ksnow*klamda1 * (0.02_r8+zsnow(i)) / &
            (ksnow*0.02_r8 + klamda1*zsnow(i))

       !    soil conductivities / delta z

       slamda(i,nsoil) = (1.0_r8-areas(i)) * 2.0_r8*klamda(i,nsoil)*klamda2 / &
            (klamda2*ztdep(i,nsoil)+0.02_r8*klamda(i,nsoil)) + &
            areas(i) * klamda(i,nsoil)*klamda1 / &
            (klamda1*ztdep(i,nsoil)*0.05_r8 + &
            klamda(i,nsoil)*(zsnow(i)+0.02_r8))

    ENDDO

    DO n = 1,nsoil-1
       DO i = 1,len
          slamda(i,n) = 2.0_r8 * klamda(i,n)*klamda(i,n+1) /   &
               (klamda(i,n)*ztdep(i,n+1)+           &
               klamda(i,n+1)*ztdep(i,n))
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE soilprop
  !
  !===================SUBROUTINE GAUSS=====================================
  !
  SUBROUTINE GAUSS(len,WORK,N,NP1,X)                         
    IMPLICIT NONE
    !========================================================================
    !
    !     SOLVE A LINEAR SYSTEM BY GAUSSIAN ELIMINATION.  DEVELOPED BY      
    !     DR. CHIN-HOH MOENG.  A IS THE MATRIX OF COEFFICIENTS, WITH THE    
    !     VECTOR OF CONSTANTS APPENDED AS AN EXTRA COLUMN.  X IS THE VECTOR 
    !     CONTAINING THE RESULTS.  THE INPUT MATRIX IS NOT DESTROYED.       
    !
    !======================================================================== 
    !          
    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: n
    INTEGER, INTENT(IN   ) :: np1                                                             
    REAL(KIND=r8)   , INTENT(INOUT) :: WORK(len,N,NP1)
    REAL(KIND=r8)   , INTENT(OUT  ) :: X   (len,N)
    !
    !     local variables
    !
    REAL(KIND=r8)    :: TEMV(len,2)
    INTEGER :: ii
    INTEGER :: j
    INTEGER :: i
    INTEGER :: k
    INTEGER :: kk
    INTEGER :: l                
    !                                                                       
    DO II=2,N                                                      
       DO J=II,N                                                      
          DO I=1,len                                                   
             TEMV(I,1) = WORK(I,J,II-1) / WORK(I,II-1,II-1)                    
          END DO
          DO K=1,NP1                                                     
             DO I=1,len                                                   
                WORK(I,J,K) = WORK(I,J,K) - TEMV(I,1) * WORK(I,II-1,K)            
             END DO
          END DO
       END DO
    END DO

    DO K=N,2,-1                                                    
       DO I=1,len                                                   
          TEMV(I,1) = WORK(I,K,NP1) / WORK(I,K,K)                           
       END DO

       KK = K-1                                                          
       DO L=KK,1,-1                                                   
          DO I=1,len                                                   
             TEMV(I,2) = TEMV(I,1) * WORK(I,L,K)                               
          END DO

          DO  I=1,len
             WORK(I,L,NP1) = WORK(I,L,NP1) - TEMV(I,2)
          END DO
       END DO
    END DO
    !
    DO II=1,N                                                      
       DO I=1,len                                                   
          X(I,II) = WORK(I,II,NP1) / WORK(I,II,II)     
       END DO
    END DO
    !                                                                       
    RETURN                                                            
  END SUBROUTINE GAUSS


  !
  !======================SUBROUTINE SIBSLV25=================================
  !
  !pl ATTENTION !!! I am calling this wiht the TGS temp instead of tg
  SUBROUTINE SIBSLV(&
       DT,&
       GRAV2,&
       CP,&
!       HLTM,&
       tg,&
!       tsnow,&
       td,&
       slamda,&
       pi     ,&
       areas  ,&
!       fac1   ,&
       VENTMF ,&
       PSB,&
!       BPS,&
!       ros,&
       psy                        ,&
       cg,ccx,cas_cap_heat,cas_cap_vap                ,&
       lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc ,&
       HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA        ,&
       HADTA,HADTH                                ,&
       hc, hg, hs, fss                                ,&
       ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                ,&
       ESDEA,EADEA,EADEM                         ,&
       ec,eg,es,fws                                ,&
!       etc,&
!       etg,&
!       ets                                ,&
!       btc,&
!       btg,&
!       bts,&
       RADT                                        ,&
       dtc, dtg, dth, dqm, dta, dea                 ,&
       len, sibdrv  )   
    !========================================================================
    !
    !     Calculation of time increments in Tc, Tgs, Theta-m and Qm using an
    !        implicit backwards method with explicit coefficients.  
    !pl   Similar to equations (10-15), SA-92B. 
    !
    !     Longwave feedbacks are now really included
    !
    !pl ATTENTION !!! this is hardwired to work only with Bonan soil.
    !pl if force-restore is wanted, we need to include approporiate if loops like
    !pl in sibslv.F
    !======================================================================== 

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       DTC            CANOPY TEMPERATURE INCREMENT (K)
    !       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
    !       DTH            MIXED LAYER POTENTIAL TEMPERATURE INCREMENT (K)
    !       DQM            MIXED LAYER MIXING RATIO INCREMENT (KG KG-1)
    !       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
    !       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE
    !           
    ! Declare parameters
    !
    INTEGER, PARAMETER :: sgl = SELECTED_REAL_KIND(p=6)   ! Single
    INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=13)  ! Double


    INTEGER, INTENT(IN   ) :: len


    REAL(KIND=r8)   , INTENT(IN   ) :: dt
    REAL(KIND=r8)   , INTENT(IN   ) :: grav2
    REAL(KIND=r8)   , INTENT(IN   ) :: cp
!    REAL(KIND=r8)   , INTENT(IN   ) :: hltm       !is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: tg    (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: tsnow (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: td    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: slamda(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: pi
    REAL(KIND=r8)   , INTENT(IN   ) :: areas (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: fac1  (len)!is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: VENTMF(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: PSB   (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: BPS   (len)!is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: ros   (len)!is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: psy   (len)    
    REAL(KIND=r8)   , INTENT(IN   ) :: cg    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ccx   (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: cas_cap_heat(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: cas_cap_vap (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lcdtc  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lcdtg  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lcdts  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lgdtg  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lgdtc  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lsdts  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: lsdtc  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HCDTC  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HCDTA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HGDTG  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HGDTA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HSDTS  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HSDTA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HADTA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: HADTH  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hc     (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hg     (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hs     (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: FSS    (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ECDTC  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ECDEA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: EGDTG  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: EGDEA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ESDTS  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ESDEA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: EADEA  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: EADEM  (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ec     (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: eg     (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: es     (len)
    REAL(KIND=r8)   , INTENT(IN   ) :: FWS    (len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: etc    (len) ! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: etg    (len) ! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: ets    (len)! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: btc    (len)! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: btg    (len)! is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: bts    (len)! is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: RADT   (len,3)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dtc    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dtg    (len,2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dth    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqm    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dta    (len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dea    (len)
    !     INTEGER, INTENT(IN   ) :: len
    LOGICAL, INTENT(IN   ) :: sibdrv
!    LOGICAL, INTENT(IN   ) :: forcerestore! is not used


    !logical :: ipvt
    !REAL(KIND=r8)    :: dtc4 (len)
    !REAL(KIND=r8)    :: dtg4 (len)
    !REAL(KIND=r8)    :: dts4(len)
    !REAL(KIND=r8)    :: z2(len)
    !REAL(KIND=r8)    :: z1(len)
    REAL(KIND=r8)    :: AAG(len)
    REAL(KIND=r8)    :: AAC(len)          
    REAL(KIND=r8)    :: AAM(len)
    REAL(KIND=r8)    :: BBG(len)
    REAL(KIND=r8)    :: BBC(len)
    REAL(KIND=r8)    :: BBM(len)
    !REAL(KIND=r8)    :: tc(len)
    !REAL(KIND=r8)    :: ta(len)
    !REAL(KIND=r8)    :: ea(len)

    !     local variables

    REAL(KIND=r8)    :: timcon
    !REAL(KIND=r8)    :: tmcn2
    REAL(KIND=r8)    :: dti
    !REAL(KIND=r8)    :: ddthl
    !REAL(KIND=r8)    :: cpdt
    REAL(KIND=r8)    :: TEMV(len)
    REAL(KIND=r8)    :: cpdpsy(len)
    REAL(KIND=r8)    :: AVEC(len,7,8)
    REAL(KIND=r8)    :: BVEC(len,7)
    !REAL(KIND=r8)*8  :: dAVEC(len,7,8)
    !REAL(KIND=r8)*8  :: dBVEC(len,7)


    INTEGER :: i
    !integer :: j
    !integer :: k
    !integer :: error_flag     

    !pl  this routine sets up the coupled system of partial 
    !pl  differential equations described in Sato et al., 
    !pl  with the exception that now Ta and ea are prognostic
    !pl  variables, and so we have two more equations, reflected
    !pl  in two more lines and two more columns as compared 
    !pl  to the old sibslv.F (used for no prognistic CAS calculations)


    !pl          J: /variables
    !pl J: equation/  1     2     3     4      5     6     7     8
    !pl              TC,   TG,   TS, TREF,  EREF,   TA,   EA,  FORCING past t.l.
    !pl 1: TC       
    !pl 2: TG
    !pl 3: TS
    !pl 4: TREF
    !pl 5: EREF
    !pl 6: TA
    !pl 7: EA


    TIMCON = PI  / 86400.0_r8     
    DTI  = 1.0_r8 / DT  
    !
    DO I=1,len 
       AAC(I) = 0.0_r8
       BBC(I) = 0.0_r8
       AAG(I) = 0.0_r8
       BBG(I) = 0.0_r8
       AAM(I) = 0.0_r8
       BBM(I) = 0.0_r8
       TEMV(I) = GRAV2 * VENTMF(i)
       cpdpsy(i) = cp / psy(i)
       !        
       !     DTC EQUATION        
       !       
       AVEC(I,1,1) = ccx(I) * DTI + HCDTC(I) + ECDTC(I) + lcdtc(i)
       AVEC(I,1,2) = LCDTG(i)
       AVEC(I,1,3) = LCDTS(i)
       AVEC(I,1,4) = 0.0_r8                                            
       AVEC(I,1,5) = 0.0_r8 
       AVEC(I,1,6) = HCDTA(I)       
       AVEC(I,1,7) = ECDEA(I)
       AVEC(I,1,8) = RADT(i,1) - HC(I) * DTI - ec(I) * DTI 
       !     
       !     DTG EQUATION      
       !      
       AVEC(i,2,1) = lgdtc(i) 
       AVEC(i,2,2) = cg(i)   * DTI + HGDTG(I) + EGDTG(I) &
            + lgdtg(i) + slamda(i)
       AVEC(i,2,3) = 0.0_r8                                    
       AVEC(I,2,4) = 0.0_r8
       AVEC(I,2,5) = 0.0_r8                                        
       AVEC(I,2,6) = hgdta(i)      
       AVEC(i,2,7) = egdea(i)                                    
       AVEC(I,2,8) = RADT(i,2) - HG(I) * DTI - eg(I) * DTI &
            - slamda(i) * (tg(I) - td(i))
       !         
       !     DTS EQUATION        
       !        
       AVEC(i,3,1) = lsdtc(i) 
       AVEC(i,3,2) = 0.0_r8
       AVEC(i,3,3) = cg(i)   * DTI + HSDTS(I) + ESDTS(I) &
            + lsdts(i) + slamda(i)                    
       AVEC(I,3,4) = 0.0_r8
       AVEC(I,3,5) = 0.0_r8                                        
       AVEC(I,3,6) = hsdta(i)      
       AVEC(i,3,7) = esdea(i)                                    
       AVEC(I,3,8) = RADT(i,3) - HS(I) * DTI - es(I) * DTI &
            - slamda(i) * (tg(I) - td(i))   
       !
       !     DTA EQUATION        
       !  
       AVEC(i,6,1) = - HCDTC(i)
       AVEC(i,6,2) = - HGDTG(i) * (1.0_r8-areas(i))
       AVEC(i,6,3) = - HSDTS(i) * (   areas(i))
       AVEC(I,6,4) =   HADTH(i)
       AVEC(I,6,5) =   0.0_r8                                   
       AVEC(I,6,6) = cas_cap_heat(i)   * DTI  &
            + HADTA(I)  - HCDTA(i)  &
            - (1.0_r8-areas(i))*HGDTA(I) - areas(i)*HSDTA(I)
       AVEC(i,6,7) = 0.0_r8          
       AVEC(I,6,8) = HC(I) * DTI - FSS(I) * DTI &
            + (1.0_r8-areas(i))*HG(I) * DTI &
            +     areas(i) *HS(I) * DTI
       !
       !     DEA EQUATION        
       !        
       AVEC(i,7,1) = - ECDTC(i)
       AVEC(i,7,2) = - EGDTG(i) * (1.0_r8-areas(i))
       AVEC(i,7,3) = - ESDTS(i) * (   areas(i))
       AVEC(I,7,4) = 0.0_r8  
       AVEC(I,7,5) = EADEM(i)
       AVEC(I,7,6) = 0.0_r8
       AVEC(i,7,7) = cas_cap_vap(i)   * DTI         &
            + (EADEA(I)  - ECDEA(I)       & 
            -  (1.0_r8-areas(i))*EGDEA(I)     &
            -      areas(i) *ESDEA(I))
       AVEC(I,7,8) = (EC(I) * DTI  -  FWS(I) * DTI  &
            +  (1.0_r8-areas(i))*EG(I) * DTI  &
            +      areas(i) *ES(I) * DTI)
    ENDDO
    !

    IF(sibdrv) THEN             
       DO i=1,len
          !       
          !     DTHETA EQUATION        
          !       
          AVEC(I,4,1) = 0.0_r8                                 
          AVEC(I,4,2) = 0.0_r8                                 
          AVEC(I,4,3) = 0.0_r8                                 
          AVEC(I,4,4) = 1.0_r8         
          AVEC(I,4,5) =  0.0_r8             
          AVEC(I,4,6) =  0.0_r8             
          AVEC(I,4,7) =  0.0_r8             
          AVEC(I,4,8) =  0.0_r8  
          !       
          !     DSH EQUATION        
          !        
          AVEC(I,5,1) = 0.0_r8
          AVEC(I,5,2) = 0.0_r8
          AVEC(I,5,3) =  0.0_r8
          AVEC(I,5,4) =  0.0_r8
          AVEC(I,5,5) = 1.0_r8
          AVEC(I,5,6) = 0.0_r8
          AVEC(I,5,7) = 0.0_r8
          AVEC(I,5,8) = 0.0_r8
       ENDDO
    ELSE
       DO i=1,len
          !        
          !     DTHETA EQUATION       
          !        
          AVEC(I,4,1) = -temv(I) * AAG(I) * (1.0_r8-areas(i))
          AVEC(I,4,2) = -temv(I) * AAG(I) * areas(i)
          AVEC(I,4,3) = -TEMV(I) * AAC(I)          
          AVEC(I,4,4) = -TEMV(I) * (AAM(I) - 1.0_r8) + &
               PSB(i) * DTI
          AVEC(I,4,5) =  0.0_r8             
          !
          !     DSH EQUATION        
          !         
          AVEC(I,5,1) = -TEMV(I)*BBG(I) * (1.0_r8-areas(i))
          AVEC(I,5,2) = -TEMV(I)*BBG(I) * areas(i)
          AVEC(I,5,3) = -TEMV(I)*BBC(I)           
          AVEC(I,5,4) =  0.0_r8  
          AVEC(I,5,5) = -TEMV(I) * (BBM(I) - 1.0_r8) + &
               PSB(i) * DTI 
          AVEC(I,4,6) = GRAV2 * FSS(i)  
          AVEC(I,5,6) = GRAV2 * FWS(i) 
       ENDDO
    ENDIF

    !      
    !     SOLVE 7 X 8 MATRIX EQUATION       
    !

    CALL GAUSS(len,AVEC,7,8,BVEC)      
    !
    DO I=1,len  
       !
       DTC(I)   = BVEC(I,1)           
       DTG(I,1) = BVEC(I,2)   ! this is DTG
       DTG(I,2) = BVEC(I,3)   ! this is DTS 
       DTH(i)   = BVEC(I,4)           
       DQM(i)   = BVEC(I,5)
       DTA(i)   = BVEC(I,6)
       DEA(i)   = BVEC(I,7)
    ENDDO


    RETURN       
  END SUBROUTINE SIBSLV


  !
  !=====================SUBROUTINE UPDAT25=================================
  !
  SUBROUTINE UPDAT2(snoww,capac,snofac,ect,eci,egi            ,&
       egs,hlat,www,pi,cg,dtd,dtg,dtc,ta,dta,dea ,&
       dtt                                            ,&
       roff,tc,td,tg,bee, poros                    ,&
       satco,slope,phsat,zdepth,ecmass            ,&
       egmass,shf,tf,snomel,asnow                    ,&
       ccx,csoil,chf,hc,hg,areas                    ,&
       q3l  ,&
       q3o  ,&
       qqq  ,&
       zmelt,&
!       cw   , &
       len  , &
       nsib ,&
       nsoil, forcerestore,etc,ea,btc             ,&
       geci,ros,cp,psy,gect,etg,btg                    ,&
       gegs,&
       hr,&
       fg,&
       gegi,&
       rd,&
       rb,&
!       hcdtc ,&
!       hcdta ,&
!       hgdta ,&
       slamda )

    IMPLICIT NONE



    !itb====================================================================
    !itb
    !itb  MOVING STORAGE TERM UPDATES (SNOW, CAPAC) HERE FROM ENDTEM, WHICH
    !itb     NO LONGER EXISTS. FLUXES PREVIOUSLY CALCULATED IN ENDTEM ARE
    !itb     TAKEN CARE OF IN THE PROGNOSTIC C.A.S. CALCULATIONS, SO WE
    !itb     MERELY NEED TO TAKE CARE OF STORAGE TERMS NOW.
    !itb
    !itb      November 2000
    !
    !=======================================================================
    !
    !     UPDATING OF ALL HYDROLOGICAL PROGNOSTIC VARIABLES.  SNOW AND
    !        RUNOFF CALCULATIONS (SEE ALSO INTER2).  SUBROUTINES SNOW2 AND
    !        RUN2 OF 1D MODEL ARE INLINED IN THIS CODE.
    !
    !=======================================================================
    !

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       DTC            CANOPY TEMPERATURE INCREMENT (K)
    !       DTD            DEEP SOIL TEMPERATURE INCREMENT (K)
    !       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
    !       WWW(3)         GROUND WETNESS 
    !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
    !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
    !       ROFF           RUNOFF (MM)
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       ECMASS         CANOPY EVAPOTRANSPIRATION (MM)
    !       EGMASS         GROUND EVAPOTRANSPIRATION (MM)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER , INTENT(IN   ) :: len
    INTEGER , INTENT(IN   ) :: nsib
    INTEGER , INTENT(IN   ) :: nsoil
    INTEGER  :: i
    INTEGER  :: l
    INTEGER  :: iveg
    INTEGER  :: ksoil
!    INTEGER  :: j

    LOGICAL , INTENT(IN   ) :: forcerestore

    !
    !-----------------------------------------------------------
    !
    !              INTENT = IN VARIABLES
    !
    !----------------------------------------------------------

    REAL(KIND=r8), INTENT(in) :: cg        (len) ! surface layer heat capacity (J m^-2 deg^-1)
    REAL(KIND=r8), INTENT(in) :: dta        (len) ! delta canopy air space (CAS) temperature (K)
    REAL(KIND=r8), INTENT(in) :: dea        (len) ! delta CAS vapor pressure (Pa)
    REAL(KIND=r8), INTENT(in) :: tc        (len) ! canopy temperature (K)
    REAL(KIND=r8), INTENT(in) :: ta        (len) ! CAS temperature (K)
    REAL(KIND=r8), INTENT(in) :: tg        (len) ! ground surface temperature (K)
    REAL(KIND=r8), INTENT(in) :: bee        (len) ! Clapp&Hornberger 'b' exponent
    REAL(KIND=r8), INTENT(in) :: poros (len) ! soil porosity (fraction)
    REAL(KIND=r8), INTENT(in) :: satco (len) ! hydraulic conductivity at saturation (UNITS?)
    REAL(KIND=r8), INTENT(in) :: slope (len) ! cosine of mean slope
    REAL(KIND=r8), INTENT(in) :: phsat (len) ! soil tension at saturation  (UNITS?)
    REAL(KIND=r8), INTENT(in) :: ccx        (len) ! canopy heat capacity (J m^-2 deg^-1)
    REAL(KIND=r8), INTENT(in) :: csoil (len) ! soil heat capacity (J m^-2 deg^-1)
    REAL(KIND=r8), INTENT(in) :: etc        (len) ! 'E-star' of the canopy - vapor pressure (Pa)
    REAL(KIND=r8), INTENT(in) :: ea        (len) ! CAS vapor pressure (Pa)
    REAL(KIND=r8), INTENT(in) :: btc        (len) ! d(E(Tc))/d(Tc) - Clausius-Clapyron
    REAL(KIND=r8), INTENT(in) :: gegs  (len) ! dry fraction of ground/(fg*rsoil + Rd)
    REAL(KIND=r8), INTENT(in) :: hr        (len) ! soil surface relative humidity
    REAL(KIND=r8), INTENT(in) :: fg        (len) ! flag indicating direction of vapor pressure
    ! deficit between CAS and ground: 0 => ea>e(Tg)
    !                                 1 => ea<e(Tg)
    REAL(KIND=r8), INTENT(in) :: gegi  (len) ! wet fraction of ground/Rd
    REAL(KIND=r8), INTENT(in) :: rd        (len) ! ground-CAS resistance
    REAL(KIND=r8), INTENT(in) :: rb        (len) ! leaf-CAS resistance
!    REAL(KIND=r8), INTENT(in) :: hcdtc (len) ! dHc/dTc
!    REAL(KIND=r8), INTENT(in) :: hcdta (len) ! dHc/dTa
!    REAL(KIND=r8), INTENT(in) :: hgdta (len) ! dHg/dTa
    REAL(KIND=r8), INTENT(in) :: slamda(len) ! 
    REAL(KIND=r8), INTENT(in) :: etg        (len) ! 'E-star' of the ground surface  (Pa)
    REAL(KIND=r8), INTENT(in) :: btg        (len) ! d(E(tg))/d(Tg) - Clausius-Clapyron
    REAL(KIND=r8), INTENT(in) :: geci  (len) ! wetted fraction of canopy/2Rb
    REAL(KIND=r8), INTENT(in) :: gect  (len) ! dry fraction of canopy/(Rst + 2Rb)
    REAL(KIND=r8), INTENT(in) :: ros        (len) ! air density (kg m^-3)
    REAL(KIND=r8), INTENT(in) :: psy        (len) ! 


    REAL(KIND=r8), INTENT(in) :: snofac      !  ___(lat ht of vap)___     (unitless)
    ! (lat ht vap + lat ht ice)
    REAL(KIND=r8), INTENT(in) :: hlat        ! latent heat of vaporization of water (J kg^-1)
    REAL(KIND=r8), INTENT(in) :: pi              ! 3.1415...
    REAL(KIND=r8), INTENT(in) :: dtt              ! timestep (seconds)
    REAL(KIND=r8), INTENT(in) :: tf              ! freezing temperature (273.16 K)
    REAL(KIND=r8), INTENT(in) :: snomel      ! latent heat of fusion for ice (J m^-3)
    REAL(KIND=r8), INTENT(in) :: asnow       ! conversion factor for kg water to depth
    !   of snow (16.7)
!    REAL(KIND=r8), INTENT(in) :: cw              ! water heat capacity (J m^-3 deg^-1)
    REAL(KIND=r8), INTENT(in) :: cp              ! specific heat of air at const pres (J kg-1 deg-1)


    REAL(KIND=r8), INTENT(in) :: zdepth(nsib,3) ! soil layer depth * porosity (meters)

    REAL(KIND=r8), INTENT(in) :: td(nsib,nsoil) ! deep soil temperature (K)


    !
    !----------------------------------------------------------------------
    !
    !              INTENT = OUT VARIABLES
    !
    !----------------------------------------------------------------------

    REAL(KIND=r8),INTENT(out) :: ect         (len)    ! transpiration flux (J m^-2 for the timestep)
    REAL(KIND=r8),INTENT(out) :: eci         (len)    ! interception flux (veg - CAS) (J m^-2)
    REAL(KIND=r8),INTENT(out) :: egi         (len)    ! ground interception flux (J m^-2)
    REAL(KIND=r8),INTENT(out) :: egs         (len)    ! ground evaporative flux (J m^-2)
    REAL(KIND=r8),INTENT(out) :: ecmass (len)    ! canopy evaporation (mm)
    REAL(KIND=r8),INTENT(out) :: egmass (len)    ! ground evaporation (mm)
    REAL(KIND=r8),INTENT(out) :: shf         (len)    ! soil heat flux (W m^-2)
    REAL(KIND=r8),INTENT(out) :: chf         (len)    ! canopy heat flux (W m^-2)
    REAL(KIND=r8),INTENT(out) :: q3l         (len)    ! 'Liston' drainage from the bottom of
    ! soil layer 3 (mm)
    REAL(KIND=r8),INTENT(out) :: q3o         (len)    ! gravitational drainage out of soil
    ! layer 3 (mm)
    REAL(KIND=r8),INTENT(out) :: zmelt  (len)    ! depth of melted water (m) 

    REAL(KIND=r8),INTENT(out) :: dtd(len,nsoil)  ! deep soil temperature increment (K)
    REAL(KIND=r8),INTENT(out) :: QQQ(len,3)      ! soil layer drainage (mm m^-2 timestep)


    !
    !----------------------------------------------------------------------
    !
    !              INTENT = IN/OUT VARIABLES
    !
    !----------------------------------------------------------------------

    REAL(KIND=r8),INTENT(inout) :: www  (len,3)  ! soil moisture (% of saturation)
    REAL(KIND=r8),INTENT(inout) :: snoww(nsib,2) ! snow-interception (1-veg, 2-ground) (meters)
    REAL(KIND=r8),INTENT(inout) :: capac(nsib,2) ! liquid interception
    REAL(KIND=r8),INTENT(inout) :: dtg  (len,2)  ! delta ground surface temperature (K)
    REAL(KIND=r8),INTENT(inout) :: roff (len)    ! runoff (mm)
    REAL(KIND=r8),INTENT(inout) :: hc   (len)    ! canopy sensible heat flux (J m^-2)
    REAL(KIND=r8),INTENT(inout) :: hg   (len)    ! ground sensible heat flux (J m^-2)
    REAL(KIND=r8),INTENT(inout) :: areas(len)    ! fractional snow coverage (0 to 1)
    REAL(KIND=r8),INTENT(inout) :: dtc  (len)    ! delta vegetation temperature (K)


    !
    !----------------------------------------------------------------------
    !
    !              LOCAL VARIABLES
    !
    !----------------------------------------------------------------------

    REAL(KIND=r8) :: dtgs          (len)  ! snow/dry soil averaged temp increment (K)
    REAL(KIND=r8) :: zmelt_total (len)  ! total melting
    REAL(KIND=r8) :: cctt          (len)
    REAL(KIND=r8) :: cct          (len)
    REAL(KIND=r8) :: ts          (len)
    REAL(KIND=r8) :: dts          (len)
    REAL(KIND=r8) :: flux          (len)
    REAL(KIND=r8) :: fluxef          (len)
    REAL(KIND=r8) :: tsnow          (len)
    REAL(KIND=r8) :: dpdw          (len)
    REAL(KIND=r8) :: tgs          (len)
    REAL(KIND=r8) :: dts2          (len)
    REAL(KIND=r8) :: cpdpsy          (len)
    REAL(KIND=r8) :: heaten          (len)   ! energy to heat snow to ground temp (J m^-2)


    REAL(KIND=r8) :: TEMW   (len,3)
    REAL(KIND=r8) :: TEMWP  (len,3)
    REAL(KIND=r8) :: TEMWPP (len,3)

    REAL(KIND=r8) :: AAA    (len,2)
    REAL(KIND=r8) :: BBB    (len,2)
    REAL(KIND=r8) :: CCC    (len,2)

    REAL(KIND=r8) :: hlati
    REAL(KIND=r8) :: rsnow
    REAL(KIND=r8) :: facks
    REAL(KIND=r8) :: realc
    REAL(KIND=r8) :: realg
    REAL(KIND=r8) :: dpdwdz
    REAL(KIND=r8) :: qmax
    REAL(KIND=r8) :: qmin 
    REAL(KIND=r8) :: rdenom
    REAL(KIND=r8) :: denom
    REAL(KIND=r8) :: props
    REAL(KIND=r8) :: avkmax
    REAL(KIND=r8) :: avkmin
    REAL(KIND=r8) :: div
    REAL(KIND=r8) :: rsame
    REAL(KIND=r8) :: pmin
    REAL(KIND=r8) :: wmin
    REAL(KIND=r8) :: pmax
    REAL(KIND=r8) :: wmax
    REAL(KIND=r8) :: egsdif
    REAL(KIND=r8) :: ectdif
    REAL(KIND=r8) :: extrak
    REAL(KIND=r8) :: facl
    REAL(KIND=r8) :: dtsg2
    REAL(KIND=r8) :: dtsg3
    REAL(KIND=r8) :: cool
    REAL(KIND=r8) :: zmelt2
    REAL(KIND=r8) :: dtsg
    REAL(KIND=r8) :: heat
    REAL(KIND=r8) :: exmelt
    REAL(KIND=r8) :: exheat
    REAL(KIND=r8) :: safe
    REAL(KIND=r8) :: avheat
    REAL(KIND=r8) :: avmelt
    REAL(KIND=r8) :: avk
    REAL(KIND=r8) :: pows
    REAL(KIND=r8) :: avex
    REAL(KIND=r8) :: freeze
    REAL(KIND=r8) :: ecpot
    REAL(KIND=r8) :: egpot
    REAL(KIND=r8) :: hrr
    REAL(KIND=r8) :: ecidif         ! actual amount of canopy moisture put into CAS (J m^-2)
    REAL(KIND=r8) :: egidif         ! actual amount of ground interception moisture
    ! put into CAS (J m^-2))
    REAL(KIND=r8)  :: egit          ! temporary ground heat flux holder (J m^-2)
    REAL(KIND=r8)  :: t1,t2          ! snow depth measures (intermediate)
    REAL(KIND=r8)  :: aven          ! energy difference between actual snow and areas=1
    REAL(KIND=r8)  :: darea          ! adjustment to areas
    REAL(KIND=r8)  :: arean          ! adjustment to areas
    REAL(KIND=r8)  :: ectmax          ! upper bound for transpiratoin (J m^-2)
    REAL(KIND=r8)  :: egsmax          ! upper bound for soil evaporation (J m^-2)
    REAL(KIND=r8)  :: dti          ! 1/timestep
    REAL(KIND=r8)  :: timcon          ! pi/seconds per day
    REAL(KIND=r8)  :: cogs1          ! non snowcovered fraction * soil humidity
    REAL(KIND=r8)  :: cogs2          ! non snowcovered fraction
    INTEGER :: k

    !
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    !

    dti = 1.0_r8/DTT
    timcon = 3.1415926_r8 / 86400.0_r8

    DO i=1,len
       tsnow(i) = MIN(tf - 0.01_r8, tg(i))
       cpdpsy(i) = cp/psy(i)
       tgs(i) = (1.0_r8 - areas(i))*tg(i) + areas(i)*tsnow(i)
       dtgs(i) = (1.0_r8-areas(i))*dtg(i,1) + areas(i)*dtg(i,2)
       rsnow = snoww(i,2) / (snoww(i,2) + capac(i,2) + 1.0E-10_r8)
       !pl this is the potential gradient in Pa
       !pl this WAS realized in sibslv

       ECPOT =     (etc(i) + btc(i)*DTC(i)) - (ea(I) + DEA(i))

       !pl and this is the  INTERCEPTION flux in J/m2
       ECI(i) = ECPOT * geci(i) * ros(i) * CPDPSY(i) * DTT

       !pl and this is the TRANSPIRATION flux in J/m2
       ECT(i) = ECPOT * gect(i) * ros(i) * CPDPSY(i) * DTT


       !pl this is the potential gradient in Pa
       !pl this WAS realized in sibslv25

       EGPOT =   (etg(i) + btg(i)*DTGS(i)) - ( ea(i) + DEA(i))

       !pl and this is the  INTERCEPTION flux in J/m2
       EGI(i) = EGPOT * (gegi(i) * (1.0_r8-areas(i)) + areas(i)/rd(i)) &  
            * ros(i) * cpdpsy(i) * DTT

       HRR = HR(i)
       IF ( FG(i) .LT. 0.5_r8 ) HRR = 1.0_r8
       COGS1    =  gegs(i) * HRR * (1.0_r8-AREAS(i))                
       COGS2    =  gegs(i)       * (1.0_r8-AREAS(i))

       !pl and this is the EVAPORATIVE flux in J/m2
       EGS(i) =  (etg(i) + btg(i) * dtgs(i)) * COGS1  &
            -(ea(i) +           dea(i)) * COGS2
       EGS(i) = EGS(i) * ros(i) * CPDPSY(i) * DTT

       !itb...make sure you don't evap more than you have...
       EGSMAX = WWW(i,1) * ZDEPTH(i,1) * hlat * 1.e3_r8 * 0.5_r8
       EGS(i) = MIN ( EGS(i), EGSMAX )
       !itb...make sure you don't transpire more water than is in the soil
       ECTMAX = WWW(i,2) * ZDEPTH(i,2) * hlat * 1.e3_r8 * 0.5_r8
       ECT(i) = MIN ( ECT(i), ECTMAX )


       !itb...these fluxes were all realized in sibslv. If positive, they
       !itb...imply transfer of water vapor INTO the CAS. If negative,
       !itb...they imply transfer OUT OF the CAS. We need to adjust
       !itb...the various reserviors as well as the CAS vapor capacity, 
       !itb...making sure that none go to negative values.

       !itb...the actual movement of vapor is taken care of in the
       !itb...equations in sibslv. All we do now is adjust the surface
       !itb...and vegetation storage reservoirs to reflect what we've
       !itb...already added or taken out.

       !pl this is the limitation to the ECI flux in J/m2

       ECIdif=MAX(0.0E0_r8,(ECI(i)-(SNOWw(i,1)+CAPAC(i,1)) &
            *1.E3_r8*hlat))

       ECI(i)   =MIN(ECI(i), &
            ( (SNOWw(i,1)+CAPAC(i,1))*1.E3_r8*hlat))


       !pl this is the EGI flux in J/m2

       EGIdif= &
            MAX(0.0E0_r8,EGI(i)-(SNOWw(i,2)+CAPAC(i,2))*1.E3_r8*hlat) &
            *(1.0_r8-RSNOW)


       EGIT  = &
            MIN(EGI(i), (SNOWw(i,2)+CAPAC(i,2))*1.E3_r8*hlat ) &
            *(1.0_r8-RSNOW)


       !itb...print this stuff out, for grins
       !        print*,'updat2: eci,ect,egi,egs,ecidif,egidif'
       !        print*,eci(i),ect(i),egi(i),egs(i),ecidif,egidif


       !
       !----------------------------------------------------------------------
       !     CALCULATION OF INTERCEPTION LOSS FROM GROUND-SNOW. IF SNOW PATCH
       !     SHRINKS, ENERGY IS TAKEN FROM EGI TO WARM EXPOSED SOIL TO TGS.
       !----------------------------------------------------------------------
       !
       T1 = SNOWw(i,2) - 1.0_r8/ASNOW
       T2 = MAX( 0.E0_r8, T1 )
       AVEN = EGI(i) - T2*hlat*1.E3_r8/SNOFAC
       IF ( (T1-T2)*EGI(i) .GT. 0.0_r8 ) AVEN = EGI(i)
       DAREA = AVEN/( (TSNOW(i)-TG(i))*CSOIL(i) &
            - 1.0_r8/ASNOW*hlat*1.E3_r8/SNOFAC)
       AREAN = AREAS(i) + DAREA
       EGIdif = EGIdif - MIN( 0.E0_r8, AREAN ) &
            *hlat*1.E3_r8/(asnow*SNOFAC)*RSNOW
       DAREA = MAX( DAREA, -AREAS(i) )
       DAREA = MIN( 1.0_r8-AREAS(i), DAREA )
       HEATEN(i) = (TSNOW(i)-TG(i))*CSOIL(i)*DAREA*RSNOW
       EGIT = EGIT + ( EGI(i) - HEATEN(i) - EGIdif )*RSNOW
       EGI(i) = EGIT


       !---------------------------------------------------------------------
       !     CALCULATION OF SENSIBLE HEAT FLUXES FOR THE END OF THE TIMESTEP.
       !        SEE FIGURE (2) OF SE-86.  NOTE THAT INTERCEPTION LOSS EXCESS
       !        ENERGIES (ECIDIF, EGIDIF) ARE ADDED.
       !
       !      HC          (HC)    : EQUATION (63) , SE-86
       !      HG          (HGS)   : EQUATION (65) , SE-86
       !----------------------------------------------------------------------
       !
       !        HC(i) = HC(i) + (HCDTC(i)*DTC(i) 
       !     &                +  HCDTA(i)*dta(i))*DTT + ECIdif
       !        HG(i) = HG(i) + (HGDTC(i)*DTC(i) 
       !     &                +  HGDTA(i)*dta(i))*DTT + EGIdif

       !itb...i've left the leaf one-sided, for now...
       HC(i) = ( (tc(i)+dtc(i)) - (ta(i)+dta(i)) ) /rb(i) &
            * ros(i) * cp * DTT + ECIdif

       !itb...ground sensible heat flux includes soil and snow by using
       !itb...dtgs
       HG(i) = ( (tg(i)+dtgs(i)) - (ta(i)+dta(i)) ) /rd(i) &
            * ros(i) * cp * DTT + EGIdif


       CHF(i) = CCX(i) * dti * DTC(i)

    ENDDO
    !----------------------------------------------------------------------
    !     CALCULATION OF STORAGE HEAT FLUXES
    !
    !---------------------------------------------------------------------- 
    !

    IF(forcerestore) THEN
       DO i = 1,len
          SHF(i) = dti * ( (1.0_r8-areas(i))*dtg(i,1) + &
               areas(i)*dtg(i,2) ) * cg(i)&
               + TIMCON*csoil(i)*2.0_r8 *( TGS(i)+dtgs(i) - TD(i,nsoil) )
       ENDDO
    ELSE
       DO i = 1,len
          !  new soil thermodynamic model         
          SHF(i) = dti * ( (1.0_r8-areas(i))*dtg(i,1) + &
               areas(i)*dtg(i,2) ) * cg(i) &
               + slamda(i) *( TGS(i)+dtgs(i) - TD(i,nsoil) )
          !           print*,'soil heat flux'
          !           print*,'dti,areas,dtg1,2,slamda,tgs,dtgs,td1'
          !           print*,dti,areas(i),dtg(i,1),dtg(i,2),slamda(i),tgs(i)
          !     &       ,dtgs(i),td(i,nsoil)
          !           print*,'SOIL HEAT FLUX = ',shf(i)
          !           print*,'soil temps'
          !           print*,'*********'
       ENDDO
    ENDIF

    IF(forcerestore) THEN
       ksoil = nsoil
    ELSE
       ksoil = 3
    ENDIF
    !
    !----------------------------------------------------------------------
    !    INTERCEPTION LOSSES APPLIED TO SURFACE WATER STORES.                      
    !    EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY
    !    ( HLAT*1000.) LOSS IS IN M M-2. MASS TERMS ARE IN KG M-2 DT-1
    !    INTERCEPTION AND DRAINAGE TREATED IN INTER2.
    !
    !      CAPAC/SNOWW(1) (M-C)   : EQUATION (3)  , SE-86
    !      CAPAC/SNOWW(2) (M-G)   : EQUATION (4)  , SE-86
    !----------------------------------------------------------------------
    ! 
    hlati = 1.0_r8 / hlat
    !PL HERE WE DO A CHECK FOR CONDENSATION AND MAKE SURE THAT IT ONLY
    !PL HAPPENS TRHOUGH ECI AND EGI

    DO i = 1,len      
       RSNOW = SNOWW(i,1)/(SNOWW(i,1)+CAPAC(i,1)+1.E-10_r8)
       FACKS = 1.0_r8 + RSNOW * ( SNOFAC-1.0_r8 )
       IF ( (ECT(i)+ECI(i)) .LE. 0.0_r8) THEN
          ECI(i) = ECT(i)+ECI(i)
          ECT(i) = 0.0_r8
          FACKS = 1.0_r8 / FACKS
       ENDIF
       CAPAC(i,1) = CAPAC(i,1)-( 1.0_r8-RSNOW )*ECI(i)*FACKS*hlati*0.001_r8
       SNOWW(i,1) = SNOWW(i,1)-     RSNOW  *ECI(i)*FACKS*hlati*0.001_r8
       snoww(i,1) = MAX(snoww(i,1),0.0E0_r8)
       capac(i,1) = MAX(capac(i,1),0.0E0_r8)
       ECMASS(i) = ECI(i)*FACKS *hlati
       zmelt_total(i) = 0.0_r8
    ENDDO
    !      do i = 1,len
    !         if(snoww(i,1).lt.0.0)
    !     *      print *,'snoww1 after updat2 100 ',i,snoww(i,1)
    !         if(capac(i,1).lt.0.0)
    !     *      print *,'capac1 after updat2 100 ',i,capac(i,1)
    !      enddo
    !
    DO i = 1,len
       RSNOW = SNOWW(i,2)/(SNOWW(i,2)+CAPAC(i,2)+1.e-10_r8)
       FACKS = 1.0_r8 + RSNOW * ( SNOFAC-1.0_r8 )
       IF ( (EGS(i)+EGI(i)) .LE. 0.0_r8 ) THEN
          EGI(i) = EGS(i)+EGI(i)
          EGS(i)= 0.0_r8
          FACKS = 1.0_r8 / FACKS
       ENDIF
       CAPAC(i,2) = CAPAC(i,2)-( 1.0_r8-RSNOW )*EGI(i)*FACKS*hlati*0.001_r8
       SNOWW(i,2) = SNOWW(i,2)-     RSNOW  *EGI(i)*FACKS*hlati*0.001_r8
       snoww(i,2) = MAX(snoww(i,2),0.0E0_r8)
       capac(i,2) = MAX(capac(i,2),0.0E0_r8)
       EGMASS(i) = EGI(i)*FACKS *hlati
    ENDDO
    !      do i = 1,len
    !         if(snoww(i,2).lt.0.0_r8)
    !     *      print *,'snoww2 after updat2 200 ',i,snoww(i,2)
    !         if(capac(i,2).lt.0.0)
    !     *      print *,'capac2 after updat2 200 ',i,capac(i,2)
    !      enddo
    !
    !----------------------------------------------------------------------
    !    DUMPING OF SMALL CAPAC VALUES ONTO SOIL SURFACE STORE
    !----------------------------------------------------------------------
    !
    DO IVEG = 1, 2
       DO i = 1,len
          IF ( (SNOWW(i,iveg)+CAPAC(i,IVEG)) .LE. 0.00001_r8 ) THEN
             WWW(i,1) = WWW(i,1) + (SNOWW(i,IVEG)+CAPAC(i,IVEG)) / &
                  ZDEPTH(i,1) 
             CAPAC(i,IVEG) = 0.0_r8
             SNOWW(i,IVEG) = 0.0_r8
          ENDIF
       ENDDO
    END DO
    DO i = 1,len
       IF(www(i,1).LT.0.0_r8)PRINT *,'www after updat2 1000 ',i,www(i,1)
    ENDDO
    !
    !
    !=======================================================================
    !------------------SNOW2 INLINED-------------------------------------
    !----------------------------------------------------------------------
    !    SNOWMELT / REFREEZE CALCULATION                                  
    !----------------------------------------------------------------------
    !                                                                     
    !     CALCULATION OF SNOWMELT AND MODIFICATION OF TEMPERATURES       
    !                                                                   
    !     MODIFICATION DEALS WITH SNOW PATCHES:                        
    !          TS < TF, TSNOW = TS                                    
    !          TS > TF, TSNOW = TF                                 
    !                                                             
    !=======================================================================
    !                                                                      
    DO IVEG = 1, 2
       !
       REALC = (2 - IVEG)*1.0_r8
       REALG = (IVEG - 1)*1.0_r8
       !
       DO i = 1,len
          CCTT(i) = REALC*CCX (i) +  REALG*CG(i)
          CCT(i)  = REALC*CCX(i)  +  REALG*CSOIL(i)
          TS(i)   = REALC*TC(i)   +  REALG*TG(i)
          DTS(i)  = REALC*DTC(i)  +  REALG*DTG(i,1)
          DTS2(i)  = REALC*DTC(i)  +  REALG*DTG(i,2)
          FLUX(i) = REALC*CHF(i)  +  REALG* &
               ( (1.0_r8-areas(i))*DTG(i,1)+ &
               areas(i)*dtg(i,2)  )*cg(i) /DTT
          !  fluxef moved up here to conserve energy
          fluxef(i) = ( shf(i) - flux(i)) * realg
          TSNOW(i) = MIN ( TF-0.01_r8, TS(i) )
          ZMELT(i) = 0.0_r8
       ENDDO
       !
       DO i = 1,len  ! this scalar loop needs vector optimization
          !itb       print*,'updat2:ts,dts,ts+dts,tf',ts(i),dts(i),ts(i)+dts(i),tf
          IF ( SNOWW(i,IVEG) .GT. 0.0_r8 ) GO TO 102
          IF ( ( TS(i)+DTS(i)) .GT. TF ) GO TO 502
          !-----------------------------------------------------------------------
          !
          !     NO SNOW  PRESENT, SIMPLE THERMAL BALANCE WITH POSSIBLE FREEZING.
          !
          !-----------------------------------------------------------------------
          FREEZE = MIN ( 0.0E0_r8, (FLUX(i)*DTT-( TF-0.01_r8 - TS(i)) &
               *CCTT(i)))
          SNOWW(i,IVEG) = MIN( CAPAC(i,IVEG), - FREEZE/SNOMEL )
          ZMELT(i) = CAPAC(i,IVEG) - SNOWW(i,IVEG)
          CAPAC(i,IVEG) = 0.0_r8
          DTS(i) = DTS(i) + SNOWW(i,IVEG)*SNOMEL/CCTT(i)
          GO TO 502
          !
          !-----------------------------------------------------------------------
          !
          !     SNOW PRESENT                                                      
          !                                                                      
          !---------------------------------------------------------------------
102       CONTINUE                                                       
          !                                                                   
          !itb      IF ( TS(i) .LT. TF .AND. (TS(i)+DTS(i)) .LT. TF ) GO TO 502
          IF ( TS(i) .LT. TF .AND. (TS(i)+DTS2(i)) .LT. TF ) GO TO 502           
          IF ( ts(i) .GT. TF ) GO TO 202                                 
          !----------------------------------------------------------------
          !                                                               
          !     SNOW PRESENT : TS < TF,  TS+DTS > TF                     
          !                                                             
          !------------------------------------------------------------
          AVEX = FLUX(i) - ( TF-0.01_r8 - TS(i) ) * CCTT(i)/DTT             
          AVMELT = ( AVEX/SNOMEL * (AREAS(i)*REALG + REALC ) )*DTT
          ZMELT(i) = MIN( AVMELT, SNOWW(i,IVEG) )   
          SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT(i)                   
          AVHEAT = AVEX*( 1.0_r8-AREAS(i) )*REALG + &
               ( AVMELT-ZMELT(i) )*SNOMEL/DTT
          AREAS(i) = MIN( 0.999E0_r8, ASNOW*SNOWW(i,2) )                        
          SAFE = MAX( ( 1.0_r8-AREAS(i)*REALG ), 1.0E-8_r8 )                      
          DTS(i) = TF-0.01_r8 - TS(i) + AVHEAT / ( CCTT(i)*SAFE )*DTT             
          GO TO 502                                                         
          !----------------------------------------------------------------------
          !                                                                     
          !     SNOW PRESENT AND TS > TF : GROUND ONLY.                        
          !                                                                   
          !------------------------------------------------------------------
202       CONTINUE                                                    
          !                                                                
          EXHEAT = CCT(i)*( 1.001_r8-MAX(0.1E0_r8,AREAS(i))) * DTS(i)              
          EXMELT = FLUX(i)*DTT - EXHEAT                             
          HEAT = EXHEAT                                         
          DTSG = EXHEAT / ( CCT(i)*(1.001_r8-AREAS(i) ))                
          IF ( (TS(i)+DTSG) .GT. TF ) GO TO 302                  
          HEAT = ( TF-0.01 - TS(i) ) * ( CCT(i)*(1.0_r8-AREAS(i)) )       
          DTSG = TF-0.01_r8 - TS(i)                               
          !                                                      
302       EXMELT = EXMELT + EXHEAT - HEAT                 
          !                                                    
          IF( EXMELT .LT. 0.0_r8 ) GO TO 402                
          ZMELT(i) = EXMELT/SNOMEL  
          IF( ASNOW*(SNOWW(i,IVEG)-ZMELT(i)) .LT. 1.0_r8 )  &                       
               ZMELT(i) = MAX( 0.0E0_r8, SNOWW(i,IVEG) - 1.0_r8/ASNOW )           
          SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT(i)                             
          EXMELT = EXMELT - ZMELT(i)*SNOMEL                               
          ZMELT2 = EXMELT/ ( CCT(i)*( TS(i)-TF )*ASNOW + SNOMEL )           
          ZMELT2 = MIN( ZMELT2, SNOWW(i,IVEG) )                        
          ZMELT(i) = ZMELT(i) + ZMELT2                                    
          SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT2                       
          EXMELT = EXMELT - ZMELT2*( CCT(i)*( TS(i)-TF )*ASNOW + SNOMEL )
          DTS(i)  = DTSG + EXMELT/CCT(i)                                
          GO TO 502                                              
          !                                                           
402       COOL = MIN( 0.0E0_r8, TF-0.01_r8 -(TS(i)+DTSG)) * CCT(i) &
               *(1.0_r8-AREAS(i))
          DTSG2 = MAX ( COOL, EXMELT ) / ( CCT(i)*( 1.001_r8-AREAS(i) ) )       
          EXMELT = EXMELT - DTSG2*CCT(i)*(1.0_r8-AREAS(i))                      
          DTSG3 =EXMELT/CCTT(i)                                         
          DTS(i) = DTSG + DTSG2 + DTSG3                                       
          !                                                                     
502       CONTINUE  
       ENDDO
       !
       DO i = 1,len
          !itb...patch
          IF(ZMELT(i) < 0.0_r8 ) ZMELT(i) = 0.0_r8
          !itb...patch
          WWW(i,1) = WWW(i,1) + ZMELT(i) / &
               ZDEPTH(i,1)  

          IF(www(i,1).LT.0.0_r8)zmelt_total(i) = SQRT(www(i,1))
          !
          DTC(i) = DTC(i)*REALG + DTS(i)*REALC
          DTG(i,1) = DTG(i,1)*REALC + DTS(i)*REALG
          zmelt_total(i) = zmelt_total(i) + zmelt(i)
       ENDDO
       !
    END DO
    !

    !itb...put zmelt_total into zmelt
    zmelt(:) = zmelt_total(:)

    IF(forcerestore) THEN
       DO k=1,nsoil
          DO i = 1,len
             !960320 fluxef calculation moved up to conserve energy
             !950511 changed cg to csoil          
             DTD(i,k) = FLUXEF(i) / &
               ( csoil(i) * 2.0_r8*SQRT( PI*365.0_r8 ) ) * DTT
          END DO
       END DO
    END IF
    DO i = 1,len
       !
       !------------------END SNOW2  -------------------------------------
       !----------------------------------------------------------------------
       !    EVAPOTRANSPIRATION LOSSES APPLIED TO SOIL MOISTURE STORE.
       !    EXTRACTION OF TRANSPIRATION LOSS FROM ROOT ZONE, SOIL EVAPORATION.
       !
       !      ECT         (E-DC)  : EQUATION (5,6), SE-86
       !      EGS         (E-S)   : EQUATION (5)  , SE-86
       !----------------------------------------------------------------------
       !
       !PL STEP THREE part II
       !pl we have done the potential ECT and EGS inside ENDTEM25
       !pl now we limit these fluxes according to 1/2 of what is 
       !pl in the soil reservoirs, WWW(i,1) for EGS and WWW(i,2) for ECT
       !pl we 'donate' the excess to HC and HG, if any.

       FACL   = hlati*0.001_r8/ZDEPTH(i,2)
       EXTRAK = ECT(i)*FACL
       EXTRAK = MIN( EXTRAK, WWW(i,2)*0.5_r8 )
       ECTDIF = ECT(i) - EXTRAK/FACL
       ECT(i)    = EXTRAK/FACL
       HC(i)     = HC(i) + ECTDIF
       ECMASS(i) = ECMASS(i) + ECT(i)*hlati
       WWW(i,2) = WWW(i,2) - ECT(i)*FACL
       !
       FACL   = 0.001_r8*hlati/ZDEPTH(i,1)
       EXTRAK = EGS(i)*FACL
       EXTRAK = MIN( EXTRAK, WWW(i,1) *0.5_r8 )
       EGSDIF = EGS(i) - EXTRAK/FACL
       EGS(i)    = EXTRAK/FACL
       HG(i)     = HG(i) + EGSDIF
       EGMASS(i) = EGMASS(i) + EGS(i)*hlati
       WWW(i,1) = WWW(i,1) - EGS(i)*FACL


    ENDDO
    !
    !========================================================================
    !------------------RUN2 INLINED-------------------------------------
    !========================================================================
    !----------------------------------------------------------------------
    !    CALCULATION OF INTERFLOW, INFILTRATION EXCESS AND LOSS TO
    !    GROUNDWATER .  ALL LOSSES ARE ASSIGNED TO VARIABLE 'ROFF' .
    !----------------------------------------------------------------------
    !
    !
    DO I = 1, 3
       !
       DO l = 1,len
          TEMW(l,I)   = MAX( 0.03E0_r8, WWW(l,I) )
          TEMWP(l,I)  = TEMW(l,I) ** ( -BEE(l) )
          TEMWPP(l,I) = MIN( 1.0E0_r8, TEMW(l,I))**( 2.0_r8*BEE(l)+ 3.0_r8 )
       ENDDO
    END DO
    !
    !-----------------------------------------------------------------------
    !
    !    CALCULATION OF GRAVITATIONALLY DRIVEN DRAINAGE FROM W(3) : TAKEN
    !    AS AN INTEGRAL OF TIME VARYING CONDUCTIVITY.
    !
    !     qqq(3) (Q3) : EQUATION (62) , SE-86
    !
    !    QQQ(3) IS AUGMENTED BY A LINEAR LOSS TERM RECOMMENDED BY LISTON (1992)
    !-----------------------------------------------------------------------
    !
    DO i = 1,len
       POWS = 2.0_r8*BEE(i)+2.0_r8
       qqq(i,3) = TEMW(i,3)**(-POWS) + SATCO(i)/ &
            (ZDEPTH(i,3) )* &
            SLOPE(i)*POWS*DTT
       qqq(i,3) = qqq(i,3) ** ( 1.0_r8 / POWS )
       qqq(i,3) = - ( 1.0_r8 / qqq(i,3) - WWW(i,3) ) * &
            ZDEPTH(i,3) / DTT
       qqq(i,3) = MAX( 0.0E0_r8, qqq(i,3) )
       q3o(i) = qqq(i,3) * dtt
       qqq(i,3) = MIN( qqq(i,3), WWW(i,3)* &
            ZDEPTH(i,3)/DTT )
       !
       Q3l(i) = 0.002_r8*ZDEPTH(i,3)*0.5_r8 / 86400.0_r8*  &
            MAX(0.0E0_r8,(www(i,3)-0.01_r8)/0.99_r8 )
       qqq(i,3) = qqq(i,3) + q3l(i)
       q3l(i) = q3l(i) * dtt
       !
       !----------------------------------------------------------------------
       !
       !    CALCULATION OF INTER-LAYER EXCHANGES OF WATER DUE TO GRAVITATION
       !    AND HYDRAULIC GRADIENT. THE VALUES OF W(X) + DW(X) ARE USED TO
       !    CALCULATE THE POTENTIAL GRADIENTS BETWEEN LAYERS.
       !    MODIFIED CALCULATION OF MEAN CONDUCTIVITIES FOLLOWS MILLY AND
       !    EAGLESON (1982 ), REDUCES RECHARGE FLUX TO TOP LAYER.
       !
       !      DPDW           : ESTIMATED DERIVATIVE OF SOIL MOISTURE POTENTIAL
       !                       WITH RESPECT TO SOIL WETNESS. ASSUMPTION OF
       !                       GRAVITATIONAL DRAINAGE USED TO ESTIMATE LIKELY
       !                       MINIMUM WETNESS OVER THE TIME STEP.
       !
       !      QQQ  (Q     )  : EQUATION (61) , SE-86
       !             I,I+1
       !            -
       !      AVK  (K     )  : EQUATION (4.14) , ME-82
       !             I,I+1
       !
       !----------------------------------------------------------------------
       !
       WMAX = MAX( WWW(i,1), WWW(i,2), WWW(i,3), 0.05E0_r8 )
       WMAX = MIN( WMAX, 1.0E0_r8 )
       PMAX = WMAX**(-BEE(i))
       WMIN = (PMAX-2.0_r8*poros(i)/( PHSAT(i)* &
            (ZDEPTH(i,1)+2.0_r8*ZDEPTH(i,2)+ZDEPTH(i,3)))) &
            **(-1.0_r8/BEE(i))
       WMIN = MIN( WWW(i,1), WWW(i,2), WWW(i,3), WMIN )
       WMIN = MAX( WMIN, 0.02E0_r8 )
       PMIN = WMIN**(-BEE(i))
       DPDW(i) = PHSAT(i)*( PMAX-PMIN )/( WMAX-WMIN )
    ENDDO
    !
    DO I = 1, 2
       !
       DO l = 1,len
          RSAME = 0.0_r8
          AVK  = TEMWP(l,I)*TEMWPP(l,I) - TEMWP(l,I+1)*TEMWPP(l,I+1)
          DIV  = TEMWP(l,I+1) - TEMWP(l,I)
          IF ( ABS(DIV) .LT. 1.E-6_r8 ) RSAME = 1.0_r8
          AVK = SATCO(l)*AVK / &
               ( ( 1.0_r8 + 3.0_r8/BEE(l) ) * DIV + RSAME )
          AVKMIN = SATCO(l) * MIN( TEMWPP(l,I), TEMWPP(l,I+1) )
          AVKMAX = SATCO(l) * MAX( TEMWPP(l,I), TEMWPP(l,I+1) )*1.01_r8
          AVK = MAX( AVK, AVKMIN )
          AVK = MIN( AVK, AVKMAX )
          !
          !c  Collatz-Bounoua change to effective hydraulic conductivity making 
          !   it 10x harder for water to move up than down if the upper soil layer
          !   is wetter than lower soil layer.
          !----------------------------------------------------------------------
          IF (www(l,i) .LT. www(l,i+1)) avk = 0.1_r8 * avk     ! lahouari
          !
          !-----------------------------------------------------------------------
          !     CONDUCTIVITIES AND BASE FLOW REDUCED WHEN TEMPERATURE DROPS BELOW
          !        FREEZING.
          !-----------------------------------------------------------------------

          TSNOW(l) = MIN ( TF-0.01_r8, TG(l) )
          TGS  (l) = TSNOW(l)*AREAS(l) + TG(l)*(1.0_r8-AREAS(l))
          TS   (l)    = TGS(l)*(2-I) + TD(l,ksoil)*(I-1)
          PROPS = ( TS(l)-(TF-10.0_r8) ) / 10.0_r8
          PROPS = MAX( 0.05E0_r8, MIN( 1.0E0_r8, PROPS ) )
          AVK   = AVK * PROPS
          qqq  (l,3)  = qqq(l,3) * PROPS

          !-----------------------------------------------------------------------
          !     BACKWARD IMPLICIT CALCULATION OF FLOWS BETWEEN SOIL LAYERS.
          !-----------------------------------------------------------------------

          DPDWDZ = DPDW(l)*2.0_r8*poros(l)/    &
               ( ZDEPTH(l,I) + ZDEPTH(l,I+1) )
          AAA(l,I) = 1.0_r8 + AVK*DPDWDZ*                &
               ( 1.0_r8/ZDEPTH(l,I)+1.0_r8/ZDEPTH(l,I+1) ) &
               *DTT
          BBB(l,I) =-AVK *   DPDWDZ * 1.0_r8/ZDEPTH(l,2)*DTT
          CCC(l,I) = AVK * ( DPDWDZ * ( WWW(l,I)-WWW(l,I+1) ) + 1.0_r8 + &
               (I-1)*DPDWDZ*qqq(l,3)*1.0_r8/ZDEPTH(l,3)* &
               DTT )
       ENDDO
    END DO
    !
    DO i = 1,len
       DENOM  = ( AAA(i,1)*AAA(i,2) - BBB(i,1)*BBB(i,2) )
       RDENOM = 0.0_r8
       IF ( ABS(DENOM) .LT. 1.E-6_r8 ) RDENOM = 1.0_r8
       RDENOM = ( 1.0_r8-RDENOM)/( DENOM + RDENOM )
       QQQ(i,1)  = ( AAA(i,2)*CCC(i,1) - BBB(i,1)*CCC(i,2) ) * RDENOM
       QQQ(i,2)  = ( AAA(i,1)*CCC(i,2) - BBB(i,2)*CCC(i,1) ) * RDENOM

       !-----------------------------------------------------------------------
       !     UPDATE WETNESS OF EACH SOIL MOISTURE LAYER DUE TO LAYER INTERFLOW
       !        AND BASE FLOW.
       !-----------------------------------------------------------------------

       WWW(i,3) = WWW(i,3) - qqq(i,3)*DTT/ZDEPTH(i,3)
       ROFF(i) = ROFF(i) + qqq(i,3) * DTT
    ENDDO

    DO I = 1, 2
       !
       DO l = 1,len
          QMAX   =  WWW(l,I)   *    &
               (ZDEPTH(l,I)  /DTT) * 0.5_r8
          QMIN   = -WWW(l,I+1) *    &
               (ZDEPTH(l,I+1)/DTT) * 0.5_r8
          QQQ(l,I) = MIN( QQQ(l,I),QMAX)
          QQQ(l,I) = MAX( QQQ(l,I),QMIN)
          WWW(l,I)   =   WWW(l,I)   - &
               QQQ(l,I)/ZDEPTH(l,I) *DTT
          WWW(l,I+1) =   WWW(l,I+1) + &
               QQQ(l,I)/ZDEPTH(l,I+1)*DTT
       ENDDO
    END DO


    DO i = 1,len
       IF(www(i,1).LT.0.0_r8)PRINT *,'www after updat2 ',i,www(i,1)
    ENDDO
    RETURN
  END SUBROUTINE UPDAT2


  SUBROUTINE SOILTHERM(&
       td, &       !td        ,& ! INTENT(IN  ) :: td(nsib,nsoil)
       dtd, &      !dtd       ,& ! INTENT(OUT ) :: dtd(len,nsoil)
       tgs, &      !tgs       ,& ! INTENT(IN  ) :: tgs(len)
       dtg, &      !dtg       ,& ! INTENT(IN  ) :: dtg(len)
       slamda, &   !slamda    ,& ! INTENT(IN  ) :: slamda(len,nsoil)
       shcap, &    !shcap     ,& ! INTENT(IN  ) :: shcap (len,nsoil)
!       ztdep, &    !ztdep     ,& ! INTENT(IN  ) :: ztdep (len,nsoil) !is not used
       dt, &       !dt        ,& ! INTENT(IN  ) :: dt
       nsib, &     !nsib      ,& ! INTENT(IN  ) :: nsib
       len, &      !len       ,& ! INTENT(IN  ) :: len
       nsoil )     !nsoil      ) ! INTENT(IN  ) :: nsoil

    IMPLICIT NONE

    !     this subroutine calculates the temperature increments dtd for
    !         the soil, based on the soil thermodynamic model of Bonan
    !     layer 1 is the deepest layer, layer nsoil is immediately below the 
    !         surface
    !     the time step is crank-nicholson. a tridiagonal matrix system is solved

    !     argument list variables
    INTEGER, INTENT(IN  ) :: len
    INTEGER, INTENT(IN  ) :: nsib
    INTEGER, INTENT(IN  ) :: nsoil

    REAL(KIND=r8)   , INTENT(IN  ) :: td    (nsib,nsoil)
    REAL(KIND=r8)   , INTENT(OUT ) :: dtd   (len,nsoil)
    REAL(KIND=r8)   , INTENT(IN  ) :: tgs   (len)
    REAL(KIND=r8)   , INTENT(IN  ) :: dtg   (len,2)
    REAL(KIND=r8)   , INTENT(IN  ) :: slamda(len,nsoil)
    REAL(KIND=r8)   , INTENT(IN  ) :: shcap (len,nsoil)
!    REAL(KIND=r8)   , INTENT(IN  ) :: ztdep (len,nsoil) !is not used
    REAL(KIND=r8)   , INTENT(IN  ) :: dt
    !     integer, INTENT(IN  ) :: nsib
    !     integer, INTENT(IN  ) :: len
    !     integer, INTENT(IN  ) :: nsoil

    !     local variables    a(n)*dtd(n-1)+b(n)*dtd(n)+c(n)*dtd(n+1) = r(n)

    REAL(KIND=r8)    :: a       (len,nsoil) ! lower sub-diagonal
    REAL(KIND=r8)    :: b       (len,nsoil) ! diagonal
    REAL(KIND=r8)    :: c       (len,nsoil) ! upper sib-diagonal
    REAL(KIND=r8)    :: r       (len,nsoil) ! right hand side
    REAL(KIND=r8)    :: lamtem
    REAL(KIND=r8)    :: rtem 
    REAL(KIND=r8)    :: dti
    REAL(KIND=r8)    :: fac
    INTEGER i, n

    dti = 1.0_r8 / dt   ! inverse time step

    !     construct matrix
    DO n = 1,nsoil
       DO i = 1,len
          b(i,n) = shcap(i,n) * dti
          r(i,n) = 0.0_r8
       ENDDO
    ENDDO

    DO n = 1,nsoil-1
       DO i = 1,len
          lamtem = -0.5_r8 * slamda(i,n)
          rtem = slamda(i,n) * (td(i,n+1) - td(i,n))
          a(i,n+1) = lamtem
          c(i,n) = lamtem
          b(i,n) = b(i,n) - c(i,n)
          b(i,n+1) = b(i,n+1) - a(i,n+1)
          r(i,n) = r(i,n) + rtem
          r(i,n+1) = r(i,n+1) - rtem
       ENDDO
    ENDDO

    DO i = 1,len
       r(i,nsoil) = r(i,nsoil) + slamda(i,nsoil) * (tgs(i)+dtg(i,1) &
            - td(i,nsoil))
    ENDDO

    !     eliminate lower diagonal
    DO n = 1,nsoil - 1
       DO i = 1,len
          fac = a(i,n+1) / b(i,n)
          b(i,n+1) = b(i,n+1) - c(i,n) * fac
          r(i,n+1) = r(i,n+1) - r(i,n) * fac
       ENDDO
    ENDDO
    !     back-substitution
    DO i = 1,len
       dtd(i,nsoil) = r(i,nsoil) / b(i,nsoil)
    ENDDO

    DO n = nsoil-1,1,-1
       DO i = 1,len
          dtd(i,n) = (r(i,n) - c(i,n) * dtd(i,n+1)) / b(i,n)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE SOILTHERM

  !
  !==================SUBROUTINE ADDINC25====================================
  !
  SUBROUTINE ADDINC(&           
!       grav2,&             !( grav2     , & ! INTENT(IN   ) :: grav2!is not used
       cp,&                !         cp             , & ! INTENT(IN   ) :: cp
!       dtt,&               !          dt              , & ! INTENT(IN        ) :: dtt  !is not used
!       hc,&                !          hc              , & ! INTENT(IN        ) :: hc(len)!is not used
!       hg, &               !          hg              , & ! INTENT(IN        ) :: hg(len)!is not used
       ps, &               !         ps             , & ! INTENT(IN   ) :: ps(len)
!       bps,&              !          bps              , & ! INTENT(IN        ) :: bps(len)!is not used
!       ecmass,&            !          ecmass      , & ! INTENT(IN        ) :: ecmass(len)!is not used
       psy,&               !         psy             , & ! INTENT(IN   ) :: psy(len)
       rho, &              !         ros             , & ! INTENT(IN   ) :: rho(len)
       hltm, &             !         hltm             , & ! INTENT(IN   ) :: hltm
!       egmass,&            !          egmass      , & ! INTENT(IN        ) :: egmass(len)!is not used
       fss,&               !         fss             , & ! INTENT(OUT  ) :: fss(len)
       fws,&               !         fws             , & ! INTENT(OUT  ) :: fws(len)
       hflux,&             !         hflux       , & ! INTENT(OUT  ) :: hflux(len)
       etmass,&            !         etmass      , & ! INTENT(OUT  ) :: etmass(len)
!       psb,&               !          psb              , & ! INTENT(IN        ) :: psb(len)!is not used
       td,&                !         td             , & ! INTENT(INOUT) :: td(nsib,nsoil)
!       thm,&               !          thmtem      , & ! INTENT(IN        ) :: thm(len)!is not used
       ts,&                !         ts             , & ! INTENT(IN   ) :: ts(len)
       sh,&                !         shtem       , & ! INTENT(OUT  ) :: sh(len)
       tc,&                !         tc             , & ! INTENT(INOUT) :: tc(len)
       tg,&                !         tg             , & ! INTENT(INOUT) :: tg(len)
       ta,&                !         ta             , & ! INTENT(INOUT) :: ta(len)
       ea,&                !         ea             , & ! INTENT(INOUT) :: ea(len)
       ra,&                !         ra             , & ! INTENT(IN   ) :: ra(len)
       em,&                !         em             , & ! INTENT(IN   ) :: em(len)
       sha,&               !         sha             , & ! INTENT(INOUT) :: sha(len)
       dtd,&               !         dtd             , & ! INTENT(IN   ) :: dtd(len,nsoil)
       dtc,&               !         dtc             , & ! INTENT(IN   ) :: dtc(len)
       dtg,&               !         dtg             , & ! INTENT(IN   ) :: dtg(len)
       dta,&               !         dta             , & ! INTENT(IN   ) :: dta(len)
       dea,&               !         dea             , & ! INTENT(IN   ) :: dea(len)
       drst,&              !         drst             , & ! INTENT(IN   ) :: DRST(len)
       rst,&               !         rst             , & ! INTENT(INOUT) :: rst(len)
       bintc,&             !         bintc       , & ! INTENT(IN   ) :: bintc(len)
       itype,&             !         itype       , & ! INTENT(IN   ) :: bintc(len)
       len,&               !         len             , & ! INTENT(IN   ) :: len
       nsib,&              !         nsib             , & ! INTENT(IN   ) :: nsib
       nsoil )             !         nsoil         ) ! INTENT(IN   ) :: nsoil

    IMPLICIT NONE
    !
    !=======================================================================
    !
    !        Add prognostic variable increments to prognostic variables
    !           and diagnose evapotranspiration and sensible heat flux.
    !
    !        Modified for multitasking - introduced gather/scatter indices
    !           - DD 951206                                  
    !
    !=======================================================================
    !

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       TC             CANOPY TEMPERATURE (K)
    !       TG             GROUND SURFACE TEMPERATURE (K)
    !       TD             DEEP SOIL TEMPERATURE (K)
    !       THM            MIXED LAYER POTENTIAL TEMPERATURE (K)
    !       QM (gsh)       MIXED LAYER MIXING RATIO (KG KG-1)
    !       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
    !pl now FWS mm/s, not mm !
    !       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
    !       rst (FSS)      STOMATAL RESISTANCE (S M-1)
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       DTH            MIXED LAYER POTENTIAL TEMPERATURE INCREMENT (K)
    !       DQM            MIXED LAYER MIXING RATIO INCREMENT (KG KG-1)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: nsib
    INTEGER, INTENT(IN   ) :: nsoil
!    REAL(KIND=r8)   , INTENT(IN   ) :: grav2!is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: cp
!    REAL(KIND=r8)   , INTENT(IN   ) :: dtt  !is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: hc(len)!is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: hg(len)!is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: ps(len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: bps(len)!is not used
!    REAL(KIND=r8)   , INTENT(IN   ) :: ecmass(len)!is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: psy(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: rho(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: hltm
!    REAL(KIND=r8)   , INTENT(IN   ) :: egmass(len)!is not used
    REAL(KIND=r8)   , INTENT(OUT  ) :: fss(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: fws(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hflux(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: etmass(len)
!    REAL(KIND=r8)   , INTENT(IN   ) :: psb(len)!is not used
    REAL(KIND=r8)   , INTENT(INOUT) :: td(nsib,nsoil)
!    REAL(KIND=r8)   , INTENT(IN   ) :: thm(len)!is not used
    REAL(KIND=r8)   , INTENT(IN   ) :: ts(len)
    REAL(KIND=r8)   , INTENT(OUT  ) :: sh(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: tc(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: tg(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: ta(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: ea(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: ra(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: em(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: sha(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: dtd(len,nsoil)
    REAL(KIND=r8)   , INTENT(IN   ) :: dtc(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: dtg(len,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: dta(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: dea(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: DRST(len)
    REAL(KIND=r8)   , INTENT(INOUT) :: rst(len)
    REAL(KIND=r8)   , INTENT(IN   ) :: bintc(len)
    INTEGER         , INTENT(IN   ) :: itype(len)
    !     integer, INTENT(IN   ) :: len
    !     integer, INTENT(IN   ) :: nsib
    !     integer, INTENT(IN   ) :: nsoil
    REAL(KIND=r8)  :: qa(len)


    !     local variables

    INTEGER :: i
    INTEGER :: n
    INTEGER :: ntyp
    REAL(KIND=r8)    :: ten,esat12

    ten = 10.0_r8
    DO I=1,len

       tc(i)   = tc(i) + dtc(i)
       tg(i)   = tg(i) + dtg(i,1)
       ta(i)   = ta(i) + dta(i) 
       ea(i)   = ea(i) + dea(i) 

       esat12 = esat (ta(i))
       qa(i)  = qsat (esat12, PS(i)*100.0_r8)

    END DO
!    CALL vnqsat(1        , & ! INTENT(IN   ) :: iflag
!               ta(1:len), & ! INTENT(IN   ) :: TQS(IM)
!               ps(1:len), & ! INTENT(IN   ) :: PQS(IM)
!               qa(1:len), & ! INTENT(OUT  ) :: QSS(IM)
!               len        ) ! INTENT(IN   ) :: IM
!    CALL QSAT_mix(qa    , &!REAL(KIND=r8), intent(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturationeing processed by qSAT scheme.
!            ta          , &!REAL(KIND=r8), intent(in)  :: T(npnts)      !  Temperature (K).
!            ps*100.0_r8 , &!REAL(KIND=r8), intent(in)  :: P(npnts)      !  Pressure (Pa).
!            len         , &!Integer, intent(in) :: npnts   !Points (=horizontal dimensions) being processed by qSAT scheme.
!            .false.     )!logical, intent(in)  :: lq_mix  .true. return qsat as a mixing ratio
    DO I=1,len

       !tc(i)   = tc(i) + dtc(i)
       !tg(i)   = tg(i) + dtg(i,1)
       !ta(i)   = ta(i) + dta(i) 
       !ea(i)   = ea(i) + dea(i) 

      ! IF(tg(i).LT.0.0_r8 .OR. ta(i).LT.0.0_r8 .OR. &
      !      ea(i).LT.0.0_r8 .OR. tc(i).LT.0.0_r8)THEN
      !    PRINT*,'BAD Ta OR ea VALUE:'
      !    PRINT*,'SiB point:',i,'itype(i)=',itype(i)
      !    PRINT*,'ea:',ea(i)-dea(i),dea(i),ea(i)
      !    PRINT*,'Ta:',ta(i)-dta(i),dta(i),ta(i)
      !    PRINT*,'Tg:',tg(i)-dtg(i,1),dtg(i,1),tg(i)
      !    PRINT*,'Tc:',tc(i)-dtc(i),dtc(i),tc(i)
      !    PRINT*,' '
      ! ENDIF


       !           print*,'addinc: new tc,tg,ta,ea=',tc,tg,ta,ea

       !pl now do the flux from the CAS to the ref level
       !pl vidale et al. (1999) equation ??
       !pl here we are using W/m2

       FSS(i) = rho(i)*cp * (ta(i)-ts(i)) / ra(i)
       hflux(i) = fss(i)

       SH(i)  = MAX(1.0e-12_r8,0.622_r8 / ( ps(i)/em(i) -1.0_r8))
       IF(ea(i) < 0.0_r8) THEN
          SHa(i) = qa(i)
       ELSE
          SHa(i) = MAX(1.0e-12_r8,0.622_r8 / ( ps(i)/ea(i) -1.0_r8))
       END IF

       !pl this is the latent heat flux from the CAS
       !pl instead of using W/m2 we stick to kg/m^2/s
       !pl the conversion is then done at the output, for
       !pl instance in RAMS module scontrol
       !pl vidale et al. (1999) equation ??

       !pl so, here we have want W/m2
       !pl in the next equation we need to multiply by cpdpsy

       fws(i) = (ea(i) - em(i)) / ra(i) * cp * rho(i) / psy(i) 

       !pl but now let us go back to mm/s (or kg/(m2 s)) for the water flux,
       !pl in order to keep to the (confusing) system we had before

       fws(i) = fws(i) / hltm
       etmass(i) = fws(i)
       ntyp=itype(i)
       IF ((ntyp == 11) .OR. (ntyp == 13)) THEN
          rst(i) = 1.0e5_r8
          CYCLE
       ELSE
          rst(i) = rst(i) + drst(i)       
       END IF
       ! bintc(i)- smallest canopy stomatal conductance needs to be passed in here.
       ! ---- c.zhang, 2/3/93
       rst(i)=MIN( 1.0_r8/bintc(i), rst(i) )
       rst(i)=MAX( ten, rst(i) )                                  
       !          rst(i)=MAX( 1.0_r8, rst(i) )                                  
    ENDDO
  
!    CALL vnqsat(1, ta(1:len), ps(1:len), qa(1:len), len ) 

    DO n = 1,nsoil
       DO I=1,len
          TD(i,n)    = TD(i,n) + dtd(i,n)      
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE ADDINC


  !
  !======================SUBROUTINE ADJUST================================
  !
  SUBROUTINE ADJUST ( TS, SPECHC, CAPACP, SNOWWP, IVEG ,&
       capac,snoww,tm,tf,snomel,www,zdepth,&
       satcap,cw,nsib, len)

    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: len
    INTEGER, INTENT(IN   ) :: nsib
    REAL(KIND=r8)   , INTENT(INOUT) :: ts
    REAL(KIND=r8)   , INTENT(IN   ) :: spechc
    REAL(KIND=r8)   , INTENT(IN   ) :: capacp
    REAL(KIND=r8)   , INTENT(IN   ) :: snowwp
    INTEGER, INTENT(IN   ) :: iveg
    REAL(KIND=r8)   , INTENT(INOUT) :: capac(nsib,2)
    REAL(KIND=r8)   , INTENT(INOUT) :: snoww(nsib,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: tm
    REAL(KIND=r8)   , INTENT(IN   ) :: tf
    REAL(KIND=r8)   , INTENT(IN   ) :: snomel
    REAL(KIND=r8)   , INTENT(INOUT) :: www
    REAL(KIND=r8)   , INTENT(IN   ) :: zdepth
    REAL(KIND=r8)   , INTENT(IN   ) :: satcap(len,2)
    REAL(KIND=r8)   , INTENT(IN   ) :: cw
    !     integer, INTENT(IN   ) :: nsib
    !     integer, INTENT(IN   ) :: len

    !
    !=======================================================================
    !
    !     TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION
    !
    !=======================================================================
    !
    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       TC             CANOPY TEMPERATURE (K)
    !       WWW(1)         GROUND WETNESS OF SURFACE LAYER 
    !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
    !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !     local variables

    REAL(KIND=r8) :: freeze
    REAL(KIND=r8) :: diff
    REAL(KIND=r8) :: ccp
    REAL(KIND=r8) :: cct
    REAL(KIND=r8) :: tsd
    REAL(KIND=r8) :: tta
    REAL(KIND=r8) :: ttb
    REAL(KIND=r8) :: cca
    REAL(KIND=r8) :: ccb
    REAL(KIND=r8) :: ccc
    REAL(KIND=r8) :: xs

    FREEZE = 0.0_r8
    DIFF = ( CAPAC(1,IVEG)+SNOWW(1,IVEG) - CAPACP-SNOWWP )*CW
    CCP = SPECHC
    CCT = SPECHC + DIFF
    !
    TSD = ( TS * CCP + TM * DIFF ) / CCT
    !
    IF ( ( TS .GT. TF .AND. TM .LE. TF ) .OR. &
         ( TS .LE. TF .AND. TM .GT. TF ) )THEN
       !
       TTA = TS
       TTB = TM
       CCA = CCP
       CCB = DIFF
       IF ( TSD .LE. TF ) THEN
          !
          !----------------------------------------------------------------------
          !    FREEZING OF WATER ON CANOPY OR GROUND
          !----------------------------------------------------------------------
          !
          CCC = CAPACP * SNOMEL
          IF ( TS .LT. TM ) CCC = DIFF * SNOMEL / CW
          TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT
          !
          FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )
          FREEZE = (MIN ( CCC, FREEZE )) / SNOMEL
          IF(TSD .GT. TF)TSD = TF - 0.01_r8
          !
       ELSE
          !
          !----------------------------------------------------------------------
          !    MELTING OF SNOW ON CANOPY OR GROUND, WATER INFILTRATES.
          !----------------------------------------------------------------------
          !
          CCC = - SNOWW(1,IVEG) * SNOMEL
          IF ( TS .GT. TM ) CCC = - DIFF * SNOMEL / CW
          !
          TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT
          !
          FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )
          FREEZE = (MAX( CCC, FREEZE )) / SNOMEL
          IF(TSD .LE. TF)TSD = TF - 0.01_r8
          !
       ENDIF
    ENDIF
    SNOWW(1,IVEG) = SNOWW(1,IVEG) + FREEZE
    CAPAC(1,IVEG) = CAPAC(1,IVEG) - FREEZE
    snoww(1,IVEG) = MAX(snoww(1,IVEG),0.0E0_r8)
    capac(1,IVEG) = MAX(capac(1,IVEG),0.0E0_r8)

    !
    XS = MAX( 0.0E0_r8, ( CAPAC(1,IVEG) - SATCAP(1,IVEG) ) )
    IF( SNOWW(1,IVEG) .GE. 0.0000001_r8 ) XS = CAPAC(1,IVEG)
    WWW = WWW + XS / ZDEPTH 
    CAPAC(1,IVEG) = CAPAC(1,IVEG) - XS
    TS = TSD

    !
    RETURN
  END SUBROUTINE ADJUST

  !
  !=====================SUBROUTINE INTER2=================================
  !
  SUBROUTINE INTER2(ppc,ppl,snoww,capac,www ,&
       satcap,cw,tc,tg,clai,zlt,chil,roff,&
       snomel,zdepTH,tm,tf,asnow,csoil,satco,dtt,vcover,roffo,&
       zmelt, len, nsib, exo)

    IMPLICIT NONE

    INTEGER , INTENT(IN   ) :: len
    INTEGER , INTENT(IN   ) :: nsib
    REAL(KIND=r8)    , INTENT(INOUT) :: ppc   (len)
    REAL(KIND=r8)    , INTENT(INOUT) :: ppl   (len)
    REAL(KIND=r8)    , INTENT(INOUT) :: snoww (nsib,2)  
    REAL(KIND=r8)    , INTENT(INOUT) :: capac (nsib,2)
    REAL(KIND=r8)    , INTENT(INOUT) :: WWW   (len,3)
    REAL(KIND=r8)    , INTENT(IN   ) :: satcap(len,2)
    REAL(KIND=r8)    , INTENT(IN   ) :: cw
    REAL(KIND=r8)    , INTENT(INOUT) :: tc(len)
    REAL(KIND=r8)    , INTENT(INOUT) :: tg(len)
    REAL(KIND=r8)    , INTENT(IN   ) :: clai
    REAL(KIND=r8)    , INTENT(IN   ) :: zlt(len)
    REAL(KIND=r8)    , INTENT(IN   ) :: chil(len)
    REAL(KIND=r8)    , INTENT(INOUT) :: roff(len)
    REAL(KIND=r8)    , INTENT(IN   ) :: snomel
    REAL(KIND=r8)    , INTENT(IN   ) :: ZDEPTH(nsib,3)
    REAL(KIND=r8)    , INTENT(IN   ) :: tm(len)
    REAL(KIND=r8)    , INTENT(IN   ) :: tf
    REAL(KIND=r8)    , INTENT(IN   ) :: asnow
    REAL(KIND=r8)    , INTENT(IN   ) :: csoil(len)
    REAL(KIND=r8)    , INTENT(IN   ) :: satco(len)
    REAL(KIND=r8)    , INTENT(IN   ) :: dtt
    REAL(KIND=r8)    , INTENT(IN   ) :: vcover(len)
    REAL(KIND=r8)    , INTENT(OUT  ) :: roffo (len)
    REAL(KIND=r8)    , INTENT(OUT  ) :: zmelt(len)
    ! INTEGER , INTENT(IN   ) :: len
    ! INTEGER , INTENT(IN   ) :: nsib
    REAL(KIND=r8)    , INTENT(OUT  ) :: exo(len)

    !=======================================================================
    !
    !     CALCULATION OF  INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW
    !     INCORPORATING EFFECTS OF PATCHY SNOW COVER AND TEMPERATURE
    !     ADJUSTMENTS.
    !
    !----------------------------------------------------------------------
    !
    !     (1) NON-UNIFORM PRECIPITATION
    !         CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY
    !         RELATIONSHIP :-
    !
    !                   F(X) = A*EXP(-B*X)+C
    !
    !         THROUGHFALL, INTERCEPTION AND INFILTRATION
    !         EXCESS ARE FUNCTIONAL ON THIS RELATIONSHIP
    !         AND PROPORTION OF LARGE-SCALE PPN.
    !         REFERENCE: SA-89B, APPENDIX.
    !
    !     (2) REORGANISATION OF SNOWMELT AND RAIN-FREEZE PROCEDURES.
    !               SUBROUTINE ADJUST
    !
    !     (3) ADDITIONAL CALCULATION FOR PARTIAL SNOW-COVER CASE.
    !               SUBROUTINE PATCHS
    !
    !     (4) REORGANISATION OF OVERLAND FLOW.
    !         REFERENCE: SA-89B, APPENDIX.
    !
    !     (5) MODIFIED CALCULATION OF SOIL HEAT CAPACITY AND
    !         CONDUCTIVITY.
    !
    !=======================================================================
    !      1D MODEL SUBROUTINE PATCHS INLINED.
    !----------------------------------------------------------------------
    !

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       ROFF           RUNOFF (MM)
    !       TC             CANOPY TEMPERATURE (K)
    !       TG             GROUND SURFACE TEMPERATURE (K)
    !       WWW(1)         GROUND WETNESS OF SURFACE LAYER 
    !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
    !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !     local variables

    INTEGER :: l

    REAL(KIND=r8)   :: excess
    REAL(KIND=r8)   :: PCOEFS(2,2)
    REAL(KIND=r8)   :: bp
    REAL(KIND=r8)   :: totalp
    REAL(KIND=r8)   :: ap(len)
    REAL(KIND=r8)   :: cp(len)
    REAL(KIND=r8)   :: thru(len)
    REAL(KIND=r8)   :: fpi(len)
    REAL(KIND=r8)   :: realc
    REAL(KIND=r8)   :: realg
    REAL(KIND=r8)   :: xss
    REAL(KIND=r8)   :: xsc
    REAL(KIND=r8)   :: chiv
    REAL(KIND=r8)   :: aa
    REAL(KIND=r8)   :: bb
    REAL(KIND=r8)   :: exrain
    REAL(KIND=r8)   :: zload
    REAL(KIND=r8)   :: xs
    REAL(KIND=r8)   :: arg
    REAL(KIND=r8)   :: tex(len)
    REAL(KIND=r8)   :: tti(len)
    REAL(KIND=r8)   :: snowwp(len)
    REAL(KIND=r8)   :: capacp(len)
    REAL(KIND=r8)   :: spechc(len)
    REAL(KIND=r8)   :: pinf(len)
    REAL(KIND=r8)   :: ts(len)
    REAL(KIND=r8)   :: equdep
    REAL(KIND=r8)   :: dcap
    REAL(KIND=r8)   :: tsd
    REAL(KIND=r8)   :: ex
    REAL(KIND=r8)   :: dareas
    REAL(KIND=r8)   :: rhs
    REAL(KIND=r8)   :: areas
    REAL(KIND=r8)   :: snowhc
    REAL(KIND=r8)   :: p0(len)
    INTEGER :: i
    INTEGER :: iveg

    DATA PCOEFS(1,1)/ 20.0_r8 /, PCOEFS(1,2)/ 0.206E-8_r8 /, &
         PCOEFS(2,1)/ 0.0001_r8 /, PCOEFS(2,2)/ 0.9999_r8 /, BP /20.0_r8 /
    !
    !-----------------------------------------------------------------------
    !
    !     PREC ( PI-X )   : EQUATION (C.3), SA-89B
    !
    !-----------------------------------------------------------------------
    !
    DO i = 1,len
       roffo(i) = 0.0_r8
       zmelt(i) = 0.0_r8
       AP(i) = PCOEFS(2,1)
       CP(i) = PCOEFS(2,2)
       TOTALP = (PPC(i) + PPL(i)) * dtt 
       IF( SNOWW(i,1) .GT. 0.0_r8 .OR. SNOWW(i,2) .GT. 0.0_r8 &
            .OR. TM(i) .LT. TF ) PPC(i) = 0.0_r8
       PPL(i) = TOTALP/dtt - PPC(i)
       IF(TOTALP.GE.1.E-8_r8) THEN
          AP(i) = PPC(i)*dtt/TOTALP * PCOEFS(1,1) + &
               PPL(i)*dtt/TOTALP * PCOEFS(2,1)
          CP(i) = PPC(i)*dtt/TOTALP * PCOEFS(1,2) + &
               PPL(i)*dtt/TOTALP * PCOEFS(2,2)
       ENDIF
       !
       THRU(i) = 0.0_r8
       FPI(i)  = 0.0_r8
       !
       !----------------------------------------------------------------------
       !     PRECIP INPUT INTO INTER2 IN M/SEC; TOTALP IS IN METERS
       !----------------------------------------------------------------------
       !
       P0(i) = TOTALP
    ENDDO

    DO IVEG = 1, 2
       REALC = 2.0_r8 - IVEG                                                
       REALG = IVEG - 1.0_r8                                                
       ! 
       DO i = 1,len                                                                     
          !
          XSC = MAX(0.0E0_r8, CAPAC(i,IVEG) - SATCAP(i,IVEG) )
          CAPAC(i,IVEG) = CAPAC(i,IVEG) - XSC
          XSS = MAX(0.0E0_r8, SNOWW(i,IVEG) - SATCAP(i,IVEG) ) * REALC
          SNOWW(i,IVEG) = SNOWW(i,IVEG) - XSS
          ROFF(i) = ROFF(i) + XSC + XSS

          CAPACP(i) = CAPAC(i,IVEG)                                             
          SNOWWP(i) = SNOWW(i,IVEG)
          !
          SPECHC(i) = &
               MIN( 0.05E0_r8, ( CAPAC(i,IVEG) + SNOWW(i,IVEG) ) ) * CW &
               + REALC * ZLT(i) * CLAI + REALG * CSOIL(i)
          TS(i) = TC(i) * REALC + TG(i) * REALG 
          !
          !----------------------------------------------------------------------
          !    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)
          !
          !     TTI ( D-D )     : EQUATION (C.4), SA-89B
          !     XS  ( X-S )     : EQUATION (C.7), SA-89B
          !     TEX ( D-C )     : EQUATION (C.8), SA-89B
          !
          !-----------------------------------------------------------------------
          !
          CHIV = CHIL(i)
          IF ( ABS(CHIV) .LE. 0.01_r8 ) CHIV = 0.01_r8
          AA = 0.5_r8 - 0.633_r8 * CHIV - 0.33_r8 * CHIV * CHIV
          BB = 0.877_r8 * ( 1.0_r8 - 2.0_r8 * AA )
          EXRAIN = AA + BB
          !
          ZLOAD = CAPAC(i,IVEG) + SNOWW(i,IVEG)
          FPI(i)=( 1.0_r8-EXP( - EXRAIN*ZLT(i)/VCOVER(i) ) )* &
               VCOVER(i)*REALC + REALG
          TTI(i) = P0(i) * ( 1.0_r8-FPI(i) )
          XS = 1.0_r8
          IF ( P0(i) .GE. 1.E-9_r8 ) THEN                                    
             ARG =  ( SATCAP(i,IVEG)-ZLOAD )/ &
                  ( P0(i)*FPI(i)*AP(i) ) -CP(i)/AP(i)               
             IF ( ARG .GE. 1.E-9_r8 ) THEN                                 
                XS = -1.0_r8/BP * LOG( ARG )                                      
                XS = MIN( XS, 1.0E0_r8 )                                            
                XS = MAX( XS, 0.0E0_r8 )                                           
             ENDIF
          ENDIF
          TEX(i) = P0(i)*FPI(i) * &
               ( AP(i)/BP*( 1.0_r8- EXP( -BP*XS )) + CP(i)*XS ) -     &
               ( SATCAP(i,IVEG) - ZLOAD ) * XS                        
          TEX(i) = MAX( TEX(i), 0.0E0_r8 ) 
       ENDDO
       !                                                              
       !-------------------------------------------------------------
       !    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION         
       !-----------------------------------------------------------
       !                                                          
       IF ( IVEG .EQ. 1 ) THEN                        
          DO i = 1,len                                           
             thru(i) = TTI(i) + TEX(i)                                  
             PINF(i) = P0(i) - THRU(i)                                 
             IF( TM(i).GT.TF ) THEN
                CAPAC(i,IVEG) = CAPAC(i,IVEG) + PINF(i) 
             ELSE
                SNOWW(i,IVEG) = SNOWW(i,IVEG) + PINF(i)               
             ENDIF
             !                                                                    
             CALL ADJUST(Tc(i)      , &! INTENT(INOUT) :: ts
                  SPECHC(i)  , &! INTENT(IN   ) :: spechc
                  CAPACP(i)  , &! INTENT(IN   ) :: capacp
                  SNOWWP(i)  , &! INTENT(IN   ) :: snowwp
                  IVEG       , &! INTENT(IN   ) :: iveg
                  capac(i,1) , &! INTENT(INOUT) :: capac(nsib,2)
                  snoww(i,1) , &! INTENT(INOUT) :: snoww(nsib,2)
                  tm(i)      , &! INTENT(IN   ) :: tm
                  tf         , &! INTENT(IN   ) :: tf
                  snomel     , &! INTENT(IN   ) :: snomel
                  www(i,1)   , &! INTENT(INOUT) :: www
                  zdepth(i,1), &! INTENT(IN   ) :: zdepth
                  satcap(i,1), &! INTENT(IN   ) :: satcap(len,2)
                  cw         , &! INTENT(IN   ) :: cw
                  nsib       , &! integer, INTENT(IN   ) :: nsib
                  len )              ! integer, INTENT(IN   ) :: len

             !                                                                   
             P0(i) = THRU(i)                                                    
          ENDDO

       ELSE
          ! 
          !DIR$ INLINE
          DO i = 1,len                                                             
             IF ( TG(i) .GT. TF .AND. SNOWW(i,2) .GT. 0.0_r8 ) THEN           
                !                                                            
                !========================================================================
                !------------------PATCHS INLINED-------------------------------------
                !========================================================================
                !
                !     CALCULATION OF EFFECT OF INTERCEPTED SNOW AND RAINFALL ON GROUND.
                !     PATCHY SNOWCOVER SITUATION INVOLVES COMPLEX TREATMENT TO KEEP
                !     ENERGY CONSERVED.
                !
                !=======================================================================

                !                                                                      
                !     MARGINAL SITUATION: SNOW EXISTS IN PATCHES AT TEMPERATURE TF    
                !     WITH REMAINING AREA AT TEMPERATURE TG > TF.                    
                !                                                                   
                !------------------------------------------------------------------
                !
                PINF(i) = P0(i)                                                   
                THRU(i) = 0.0_r8                                                  
                SNOWHC = MIN( 0.05E0_r8, SNOWW(i,2) ) * CW                       
                areas = MIN( 1.0E0_r8,(ASNOW*SNOWW(i,2)) )                       
                IF( TM(i) .LE. TF ) THEN                              
                   !                                                            
                   !-----------------------------------------------------------
                   !     SNOW FALLING ONTO AREA                               
                   !---------------------------------------------------------
                   !                                                        
                   RHS = TM(i)*PINF(i)*CW + TF*(SNOWHC + &
                        CSOIL(i)*areas)      &
                        + TG(i)*CSOIL(i)*(1.0_r8-areas)                        
                   DAREAS = MIN( ASNOW*PINF(i), ( 1.0_r8-areas ) )                       
                   EX = RHS - TF*PINF(i)*CW -  &
                        TF*(SNOWHC + CSOIL(i)*(areas + DAREAS))   &
                        - TG(i)*CSOIL(i)*(1.0_r8-areas-DAREAS)                              
                   IF( (areas+DAREAS) .GE. 0.999_r8 ) &
                        TG(i) = TF - 0.01_r8              
                   IF( EX .GE. 0.0_r8 ) THEN                                 
                      !                                                               
                      !----------------------------------------------------------------------
                      !     EXCESS ENERGY IS POSITIVE, SOME SNOW MELTS AND INFILTRATES. 
                      !---------------------------------------------------------------------
                      !                                                                    
                      ZMELT(i) = EX/SNOMEL                                             
                      IF( ASNOW*(SNOWW(i,2) + PINF(i) - ZMELT(i)) &
                           .LE. 1.0_r8 ) THEN      
                         ZMELT(i) = 0.0_r8                                                  
                         IF( ASNOW*(SNOWW(i,2) + PINF(i)) .GE. 1.0_r8 ) &
                              ZMELT(i) = ( ASNOW*(SNOWW(i,2) + &
                              PINF(i)) - 1.0_r8 ) / ASNOW      
                         ZMELT(i) = ( EX - ZMELT(i)*SNOMEL )/ &
                              ( SNOMEL + ASNOW*CSOIL(i)*   &
                              (TG(i)-TF) ) + ZMELT(i)                                                   
                      ENDIF
                      SNOWW(i,2) =  SNOWW(i,2) + PINF(i) - ZMELT(i)                            
                      WWW(i,1) = WWW(i,1) + ZMELT(i)/ZDEPTH(i,1)
                   ELSE                                                         
                      !
                      !----------------------------------------------------------------------
                      !     EXCESS ENERGY IS NEGATIVE, BARE GROUND COOLS TO TF, THEN WHOLE
                      !     AREA COOLS TOGETHER TO LOWER TEMPERATURE.
                      !----------------------------------------------------------------------
                      !
                      TSD = 0.0_r8
                      IF( (areas+DAREAS) .LE. 0.999_r8 ) &
                           TSD = EX/(CSOIL(i)*( 1.0_r8-areas-DAREAS)) &
                           + TG(i)
                      IF( TSD .LE. TF )  &
                           TSD = TF + ( EX - (TF-TG(i))* &
                           CSOIL(i)*(1.0_r8-areas-DAREAS) ) &
                           /(SNOWHC+PINF(i)*CW+CSOIL(i))
                      TG(i) = TSD
                      SNOWW(i,2) = SNOWW(i,2) + PINF(i)
                   ENDIF
                ELSE
                   !
                   !----------------------------------------------------------------------
                   !     RAIN FALLING ONTO AREA
                   !----------------------------------------------------------------------
                   !
                   !----------------------------------------------------------------------
                   !     RAIN FALLS ONTO SNOW-FREE SECTOR FIRST.
                   !----------------------------------------------------------------------
                   TSD = TF - 0.01_r8
                   IF ( areas .LT. 0.999_r8 ) TSD =  &
                        ( TM(i)*PINF(i)*CW + &
                        TG(i)*CSOIL(i) )   &
                        /  ( PINF(i)*CW + CSOIL(i) )
                   TG(i) = TSD
                   WWW(i,1)= WWW(i,1)+PINF(i)*(1.0_r8-areas)/ &
                        ZDEPTH(i,1)
                   !----------------------------------------------------------------------
                   !     RAIN FALLS ONTO SNOW-COVERED SECTOR NEXT.
                   !----------------------------------------------------------------------
                   EX = ( TM(i) - TF )*PINF(i)*CW*areas
                   DCAP = -EX / ( SNOMEL + ( TG(i)-TF )* &
                        CSOIL(i)*ASNOW )
                   IF( (SNOWW(i,2) + DCAP) .GE. 0.0_r8 ) THEN
                      WWW(i,1) = WWW(i,1)+(PINF(i)*areas-DCAP)/ &
                           ZDEPTH(i,1)
                      SNOWW(i,2) = SNOWW(i,2) + DCAP
                   ELSE
                      TG(i) = ( EX - SNOMEL*SNOWW(i,2) - &
                           ( TG(i)-TF )*CSOIL(i)*areas ) / &
                           CSOIL(i) + TG(i)
                      WWW(i,1)=WWW(i,1)+(SNOWW(i,2)+PINF(i)* &
                           areas)/zdepth(i,1)
                      CAPAC(i,2) = 0.0_r8
                      SNOWW(i,2) = 0.0_r8
                   ENDIF
                   !
                ENDIF
                !
                !-----------------------------------------------------------------------
                !---------------------END OF PATCHS ------------------------------------
                !-----------------------------------------------------------------------

             ELSE
                !                                                                   
                THRU(i) = TTI(i) + TEX(i)                                             
                IF ( TG(i) .LE. TF .OR. TM(i) .LE. TF ) &
                     THRU(i) = 0.0_r8                 
                PINF(i) = P0(i) - THRU(i)                                           
                IF( TM(i) .GT. TF )THEN                                        
                   CAPAC(i,IVEG) = CAPAC(i,IVEG) + PINF(i)         
                   !                                                                      
                   !---------------------------------------------------------------------
                   !                                                                    
                   !    INSTANTANEOUS OVERLAND FLOW CONTRIBUTION ( ROFF )              
                   !                                                                  
                   !     ROFF( R-I )     : EQUATION (C.13), SA-89B                   
                   !                                                                
                   !---------------------------------------------------------------
                   !                                                              
                   EQUDEP = SATCO(i) * DTT                                    
                   !                                                            
                   XS = 1.0_r8                                               
                   IF ( THRU(i) .GE. 1.E-9_r8 ) THEN                     
                      ARG = EQUDEP / ( THRU(i) * AP(i) ) -CP(i)/AP(i)                 
                      IF ( ARG .GE. 1.E-9_r8 ) THEN                    
                         XS = -1.0_r8/BP * LOG( ARG )                         
                         XS = MIN( XS, 1.0E0_r8 )                               
                         XS = MAX( XS, 0.0E0_r8 )                              
                      ENDIF
                   ENDIF
                   ROFFO(i) = THRU(i) * ( AP(i)/BP * ( 1.0_r8-EXP( -BP*XS )) + CP(i)*XS ) -EQUDEP*XS
                   ROFFO(i) = MAX ( ROFFO(i), 0.0E0_r8 )
                   ROFF(i) = ROFF(i) + ROFFO(i)
                   WWW(i,1) = WWW(i,1) + (THRU(i) - ROFFO(i)) / ZDEPTH(i,1) 
                ELSE
                   SNOWW(i,IVEG) = SNOWW(i,IVEG) + PINF(i)        
                ENDIF
                !
                !----------------------------------------------------------------------
                !
                CALL ADJUST ( Tg(i)      , & ! INTENT(INOUT) :: ts
                     SPECHC(i)  , & ! INTENT(IN   ) :: spechc
                     CAPACP(i)  , & ! INTENT(IN   ) :: capacp
                     SNOWWP(i)  , & ! INTENT(IN   ) :: snowwp
                     IVEG       , & ! INTENT(IN   ) :: iveg
                     capac(i,1) , & ! INTENT(INOUT) :: capac(nsib,2)
                     snoww(i,1) , & ! INTENT(INOUT) :: snoww(nsib,2)
                     tm(i)      , & ! INTENT(IN   ) :: tm
                     tf         , & ! INTENT(IN   ) :: tf
                     snomel     , & ! INTENT(IN   ) :: snomel
                     www(i,1)   , & ! INTENT(INOUT) :: www
                     zdepth(i,1), & ! INTENT(IN   ) :: zdepth
                     satcap(i,1), & ! INTENT(IN   ) :: satcap(len,2)
                     cw         , & ! INTENT(IN   ) :: cw
                     nsib       , & ! INTENT(IN   ) :: nsib
                     len          ) ! INTENT(IN   ) :: len

                !
                !----------------------------------------------------------------------
                !
             ENDIF
          ENDDO
       ENDIF   ! if(iveg.eq.1)

       !     make either all capac or all snow

       DO i = 1,len
          IF(capac(i,iveg).GT.snoww(i,iveg)) THEN
             capac(i,iveg) = capac(i,iveg) + snoww(i,iveg)
             snoww(i,iveg) = 0.0_r8
          ELSE
             snoww(i,iveg) = snoww(i,iveg) + capac(i,iveg)
             capac(i,iveg) = 0.0_r8
          ENDIF
       ENDDO
       !
    END DO
    !
    DO i = 1,len
       exo(i) = 0.0_r8
    ENDDO

    DO I = 1, 3
       DO l = 1,len
          EXCESS = MAX(0.0E0_r8,(WWW(l,I) - 1.0_r8))
          WWW(l,I) = WWW(l,I) - EXCESS
          exo(l) = exo(l) + EXCESS * ZDEPTH(l,I)
          !
          ! Collatz-Bounoua put excess water into runoff according to
          ! original sib2 offline code .
          !
          roff(l) = roff(l) + EXCESS * ZDEPTH(l,I)    ! lahouari
       ENDDO

    END DO

    RETURN
  END SUBROUTINE INTER2




  SUBROUTINE CFRAX(dtt,d13Cca,d13Cm,d13Cresp,pco2a,pco2m,pco2s, &
       pco2c,pco2i,respg,assimn,ga,KIECpsC3,KIECpsC4, &
       d13CassimnC3,d13CassimnC4,d13Cassimntot,Flux13C,Flux12C, &
       len,c4fract,ca_depth,t_air,flux_turb)

    !  CFRAX calculates 13C and 12C fluxes and concentrations in the canopy,
    !  mixed layer, carbon assimilation (photosynthate), respired soil carbon,
    !  assuming that discrimination against 13C during photosynthesis is a 
    !  function of discrimination during diffusion and the pCO2i/pCO2c ratio.
    !  C4 discrimination against 13C only results from diffusion. 


    !itb...modified 01 Oct 99 to try to settle some mass-balance problems
    !itb...we've been having - new implicit eqn's for canopy C12 and C13.

    IMPLICIT NONE

    ! input variables

    INTEGER :: i
    INTEGER , INTENT(IN   ) :: len
    REAL(KIND=r8)    , INTENT(IN   ) :: dtt      !the time  step in seconds
    REAL(KIND=r8)    , INTENT(INOUT) :: d13Cca   (len)       !del13C of canopy CO2 (per mil vs PDB)
    REAL(KIND=r8)    , INTENT(IN   ) :: d13Cm    (len)       !del13C of mixed layer CO2 (per mil vs PDB)
    REAL(KIND=r8)    , INTENT(IN   ) :: d13Cresp (len)       !del13C of respiration CO2 (per mil vs PDB)   
    REAL(KIND=r8)    , INTENT(IN   ) :: pco2a    (len)       !CO2 pressure in the canopy (pascals)      
    REAL(KIND=r8)    , INTENT(IN   ) :: pco2m    (len)       !CO2 pressure in the mixed layer (pascals)
    REAL(KIND=r8)    , INTENT(IN   ) :: pco2s    (len)              !CO2 pressure at the leaf surface (pascals)
    REAL(KIND=r8)    , INTENT(IN   ) :: pco2c    (len)       !CO2 pressure in the chloroplast (pascals)
    REAL(KIND=r8)    , INTENT(IN   ) :: pco2i          (len)       !CO2 pressure in the stoma (pascals)
    REAL(KIND=r8)    , INTENT(IN   ) :: respg    (len)       !rate of ground respiration (moles/m2/sec)      
    REAL(KIND=r8)    , INTENT(IN   ) :: assimn   (len)       !rate of assimilation by plants (moles/m2/sec)
    REAL(KIND=r8)    , INTENT(INOUT) :: ga       (len)       !1/resistance factor (ra) to mixing between the canopy 
    !and the overlying mixed layer. ga units:(m/sec)
    REAL(KIND=r8)    , INTENT(OUT   ) :: KIECpsC3 (len)      !Kinetic Isotope Effect during C3 photosynthesis
    REAL(KIND=r8)    , INTENT(OUT   ) :: KIECpsC4 (len)      !Kinetic Isotope Effect during C4 photosynthesis
    REAL(KIND=r8)    , INTENT(OUT   ) :: d13CassimnC3  (len) !del13C of CO2 assimilated by C3 plants      
    REAL(KIND=r8)    , INTENT(OUT   ) :: d13CassimnC4  (len) !del13C of CO2 assimilated by C3 plants
    REAL(KIND=r8)    , INTENT(OUT   ) :: d13Cassimntot (len) !del13C of CO2 assimilated by all plants
    REAL(KIND=r8)    , INTENT(OUT   ) :: Flux13C       (len) !turbulent flux of 13CO2 out of the canopy
    REAL(KIND=r8)    , INTENT(OUT   ) :: Flux12C       (len) !turbulent flux of 12CO2 out of the canopy
    !the canopy 
    !      INTEGER , INTENT(IN   ) :: len
    REAL(KIND=r8)    , INTENT(IN    ) ::c4fract (len)        !fraction of vegetation that is C4: C. Still, pers.comm.

    REAL(KIND=r8)    , INTENT(IN    ) :: ca_depth (len)      !depth of canopy
    REAL(KIND=r8)    , INTENT(IN    ) :: t_air(len)          ! canopy temperature (K)
    REAL(KIND=r8)    , INTENT(OUT   ) :: flux_turb(len)      ! turbulent flux of CO2 out of the
    ! CO2 out of the canopy    



    ! Temporarily the division of C3 and C4 assimilation rates will be calculated in CFRAX.
    ! Therefore the following two variables are local.

    REAL(KIND=r8)    :: assimnC3(len)!rate of assimilation by C3 plants (moles/m2/sec)    
    REAL(KIND=r8)    :: assimnC4(len)!rate of assimilation by C4 plants(moles/m2/sec)


    ! local variables  (concentration units: moles/m3; flux units:moles/m2/sec)

    REAL(KIND=r8)    :: Rca             (len)!C13/C12 ratio of canopy CO2 
    REAL(KIND=r8)    :: c13ca       (len)!concentration of C13 in canopy CO2 
    REAL(KIND=r8)    :: c12ca       (len)!concentration of C12 in canopy CO2
    REAL(KIND=r8)    :: Rcm             (len)!C13/C12 ratio of mixed layer CO2
    REAL(KIND=r8)    :: c13cm       (len)!concentration of C13 in mixed layer CO2
    REAL(KIND=r8)    :: c12cm       (len)!concentration of C12 in mixed layer CO2
    REAL(KIND=r8)    :: Rcresp      (len)!C13/C12 ratio of respiratory CO2
    REAL(KIND=r8)    :: c13resp     (len)!flux of C13 in respiratory CO2
    REAL(KIND=r8)    :: c12resp     (len)!flux of C12 in respiratory CO2
    REAL(KIND=r8)    :: RcassimnC3  (len)!C13/C12 ratio of CO2 assimilated by C3 plants
    REAL(KIND=r8)    :: c13assimnC3 (len)!flux of C13 in CO2 assimilated by C3 plants
    REAL(KIND=r8)    :: c12assimnC3 (len)!flux of C13 in CO2 assimilated by C3 plants
    REAL(KIND=r8)    :: RcassimnC4  (len)!C13/C12 ratio of CO2 assimilated by C4 plants
    REAL(KIND=r8)    :: c13assimnC4 (len)!flux of C13 in CO2 assimilated by C3 plants
    REAL(KIND=r8)    :: c12assimnC4 (len)!flux of C12 in CO2 assimilated by C3 plants
    REAL(KIND=r8)    :: c13assimntot(len)!total flux of C13 in CO2 assimilated by plants     
    REAL(KIND=r8)    :: c12assimntot(len)!total flux of C12 in CO2 assimilated by plants      
    REAL(KIND=r8)    :: Rcassimntot (len)!C13/C12 ratio of CO2 assimilated by all plants

    REAL(KIND=r8)    :: co2a_conc (len)
    REAL(KIND=r8)    :: co2m_conc (len) ! canopy and mixed layer CO2 concentrations 
    ! in moles per meter cubed rather
    ! than partial pressure (pascals)

    REAL(KIND=r8)    :: ga_temp (len)              ! local copy of aero conductance
    ! store real value here while we mess
    ! with the ga to prevent excessive stability



    !NOTE: IT MAY BE NECESSARY TO ALSO DIMENSION BY MONTH

    ! Carbon isotopic fractionation constants (units = per mil)         

    ! KIEC refers to Kinetic Isotope Effect (KIE) for Carbon (C), and can be converted   
    ! to alpha notation by alpha = (1 - KIEC/1000).  For a chemical reaction,
    ! alpha = Rreactant/Rproduct.  KIEs are sometimes referred to as epsilon factors.  

    REAL(KIND=r8), PARAMETER :: PDB = 0.0112372_r8   ! 13C/12C ratio of Pee Dee Belemnite (no units)

    REAL(KIND=r8), PARAMETER :: KIEClfbl  = - 2.9_r8    ! canopy air space to leaf boundary layer
    REAL(KIND=r8), PARAMETER :: KIECstom  = - 4.4_r8    ! leaf boundary layer to stomatal cavity
    REAL(KIND=r8), PARAMETER :: KIEClphas =  -0.7_r8    ! liquid phase fractionation 
    REAL(KIND=r8), PARAMETER :: KIECdis   =  -1.1_r8    ! dissolution

    REAL(KIND=r8), PARAMETER :: KIECrbsco = -28.2_r8    ! C3 C-fixation enzyme rubisco

    REAL(KIND=r8), PARAMETER :: TREF = 298.16_r8       ! standard temperature (K)
    REAL(KIND=r8), PARAMETER :: PREF = 101325.0_r8     ! standard pressure (Pa)

    !  output variables 




    !  ***********************************************************************************************

    ! I am temporarily hardwiring the carbon isotopic ratios of the mixed layer
    ! and respireation into SiBDRIVE

    !        d13Cm = -7.8_r8
    !        d13Cresp = -28.0_r8

    DO i = 1,len 



       !itb...converts canopy and mixed layer CO2 pressures from Pa to moles/m^3 
       ! uses PV = nRT; at STP and 44.6 moles of gas per m3.  


       co2a_conc(i) = pco2a(i) * (TREF/t_air(i))*(44.6_r8/PREF)

       co2m_conc(i) = pco2m(i) * (TREF/t_air(i))*(44.6_r8/PREF)


       !  the conductance is given a lower limit of 0.01 m/sec.  This is the same
       !  as the limit provided in phosib.F, where gah2o = max(0.446, gah2o),
       !  because the units are different, i.e., moles/m2/sec vs. m/s.             


       ga_temp(i) = ga(i)            ! remember the real ga
       !             ga(i) = max(ga(i), 1./1000.)   ! bottom-stop ga at 1/100


    END DO

    !  d13Cca and d13Cm are converted to concentrations (moles/m3) of 13C 
    !  and 12C by first calculating isotope ratios (13C/12C) of the canopy         
    ! (Ca) and mixed layer (m). 

    DO i = 1,len            
       Rca(i)  = ((d13Cca(i) * PDB) / 1000.0_r8) + PDB

       c13ca(i)= (Rca(i) * co2a_conc(i)) / (1.0_r8 + Rca(i))
       c12ca(i) = co2a_conc(i) / (1.0_r8 + Rca(i))

       Rcm(i)  = ((d13Cm(i) * PDB) / 1000.0_r8) + PDB


       c13cm(i)= (Rcm(i) * co2m_conc(i)) / (1.0_r8 + Rcm(i))
       c12cm(i) = co2m_conc(i) / (1.0_r8 + Rcm(i))

       !            print*,'co2a_local,c13ca,c12ca',
       !     +            co2a_local(i),c13ca(i),c12ca(i)
       !            print*,'co2m_local,c13cm,c12cm',
       !     +            co2m_local(i),c13cm(i),c12cm(i) 
       !            print*,' '

    END DO

    !  13c and 12c fluxes (moles/m2/sec) arising from respiration are        
    !  calculated using conversions between delta notation and epsilon 
    !  notation and 13C/12C ratios.

    DO i = 1,len
       Rcresp(i)  = ((d13Cresp(i) * PDB) / 1000.0_r8) + PDB
       c13resp(i) = (Rcresp(i) * respg(i)) / (1.0_r8 + Rcresp(i))
       c12resp(i) = respg(i) / (1.0_r8 + Rcresp(i))
    END DO

    ! C13/C12 discrimination for C3 plants.  The isotope effect during C3
    ! photosynthesis is a function of a combination of the isotope effects
    ! associated with molecular transport of CO2 across the leaf boundary
    ! layer (lfbl), into the stoma (stom), dissolution to in mesophyll H2O
    ! (dis), and transport in the liquid phase (lphas).  The isotope effect 
    ! during C4 photosynthesis is only a function (for now) of transport into
    ! the stoma.    note: IECpsC3 is the isotope effect for carbon isotopic  
    ! discrimination during photosynthesis of C3 plants.  Similarly for IECpsC4,
    ! but for C4 plants. 


    DO i = 1,len

       IF(assimn(i).GT. 0.0_r8) THEN

          KIECpsC3(i)=(KIEClfbl*pco2a(i)+(KIECstom-KIEClfbl)*pco2s(i)+ &
               (KIECdis+KIEClphas-KIECstom)*pco2i(i)+ &
               (KIECrbsco-KIECdis-KIEClphas)*pco2c(i))/pco2a(i)

          KIECpsC4(i)= KIECstom 

       ELSE 

          ! Since we are temporarily using net assimilation (Assimn), we set
          ! del13C value of leaf respiration to the del13C of the soild organic
          ! matter.  The atmospheric effect on carbon isotope ratios of plants
          ! is deducted here, and added back as the carbon isotope ratios of the
          ! canopy later on.

          KIECpsC3(i)= d13Cresp(i) - d13cm(i)
          KIECpsC4(i)= d13Cresp(i) - d13cm(i)

       END IF


    END DO

    ! Total net assimilation is divided between C3 and C4 plants using monthly
    ! maps developed by Chris Still.  C4fract ranges from 0.00 to 1.00

    DO i = 1,len 

       assimnC4(i) = Assimn(i) * c4fract(i)
       assimnC3(i) = Assimn(i) * (1.0_r8 - c4fract(i))

    END DO


    ! calculates del13C value of carbon assimilated by C3 plants as well as
    ! the fluxes of 13C and 12C assimilated by C3 plants (moles/m2/sec)


    DO i = 1,len 

       RcassimnC3(i)  =   Rca(i) / ((-KIECpsC3(i)/1000.0_r8) + 1.0_r8)       

       d13CassimnC3(i) =  ((RcassimnC3(i) - PDB)/ PDB) * 1000.0_r8

       !        d13CassimnC3(i) = d13Cca(i) + KIECpsC3(i) 

       !        RcassimnC3(i)  = (d13CassimnC3(i) * PDB)/1000.0_r8  + PDB 

       c13assimnC3(i) = (RcassimnC3(i) * assimnC3(i))/ &
            (1.0_r8 + RcassimnC3(i))

       c12assimnC3(i) =  assimnC3(i) / (1.0_r8+ RcassimnC3(i))


       ! calculates del13C value of carbon assimilated by C4 plants as well as
       ! the fluxes of 13C and 12C assimilated by C4 plants (moles/m2/sec)

       RcassimnC4(i)  =   Rca(i) /  ((-KIECpsC4(i)/1000.0_r8) + 1.0_r8)      

       d13CassimnC4(i) =  ((RcassimnC4(i) - PDB)/ PDB) * 1000.0_r8

       !        d13CassimnC4(i) = d13Cca(i) + KIECpsC4(i) 

       !        RcassimnC4(i)  = (d13CassimnC4(i) * PDB)/1000.  + PDB

       c13assimnC4(i) = (RcassimnC4(i) * assimnC4(i))/ &
            (1.0_r8 + RcassimnC4(i))

       c12assimnC4(i) =  assimnC4(i) / (1.0_r8+ RcassimnC4(i))


       ! Total assimilated fluxes of 13C and 12C, as well as d13C of the
       ! assimilated flux are calculated (moles/m2/sec and 
       ! (per mil moles/m2/sec)) 

       c13assimntot(i) = c13assimnC3(i) + c13assimnC4(i)
       c12assimntot(i) = c12assimnC3(i) + c12assimnC4(i)

       d13Cassimntot(i) = (((c13assimntot(i) /  &
            c12assimntot(i)) - PDB) / PDB) * 1000.0_r8

       Rcassimntot(i)  = (d13Cassimntot(i) * PDB)/1000.0_r8  + PDB


       !calculates retrodiffused flux of co2 using the flux-gradient 
       !proportionality found in keeling (1996)
       !  I believe that this is only used for oxygen isotopes, so I have 
       !  commented it out for the time being.  nss

       !      retroflux(i) = ASSIM(i)*(pco2c(i)/(pco2a(i)-
       !     & pco2c(i)))       

    END DO


    !  Canopy concentrations at time n+1 is calculated using an implicit 
    !  scheme.

    DO i = 1,len 


       c13ca(i) = (c13ca(i) + (dtt / ca_depth(i)) * (c13resp(i) - &
            c13assimntot(i) + (c13cm(i) * ga(i))))  &
            / (1.0_r8 + (dtt*ga(i)) / ca_depth(i))

       c12ca(i) = (c12ca(i) + (dtt / ca_depth(i)) * (c12resp(i) - &
            c12assimntot(i) + (c12cm(i) * ga(i)))) &
            / (1.0_r8 + (dtt*ga(i)) / ca_depth(i))


    END DO


    !  del13C of canopy is recalculated using concentrations of 13Cca and 
    !  12Cca.  The fluxes (moles/m2/sec) of 13C and 12C out of the canopy 
    !  (the turbulent flux), and the del13C value (per mil vs PDB) of this 
    !  flux are calculated. 

    DO i = 1,len

       d13Cca(i) = ((c13ca(i)/c12ca(i) - PDB) / PDB) *1000.0_r8


       !  Use the following if you want the net flux from the canopy to 
       !  be based on differences in 12C and 13C net fluxes from respiration
       !  and photosynthesis

       !           Flux13C(i) = respg(i) * Rcresp(i) / (1. + Rcresp(i)) -
       !     &  (assimn(i) * Rcassimntot(i) / (1. + Rcassimntot(i)))

       !           Flux12C(i) = respg(i) / (1. + Rcresp(i)) -
       !     &  (assimn(i) / (1. + Rcassimntot(i)))

       !  Use the following if you want the net flux from the canopy to 
       !  be based on differences  in concentration gradients between the 
       !  canopy and overlying atmosphere.

       Flux13C(i) = ( c13ca(i) - c13cm(i)) * ga(i)
       Flux12C(i) = ( c12ca(i) - c12cm(i)) * ga(i)

       ! To prevent generation of NaNs by dividing by zero: del13C of the 
       ! turbulent flux is set the del13c of respiration.   nss

       !            IF   (Flux12C(i) .ne. 0.0) THEN
       !               d13CFlux(i) = (((Flux13C(i) / Flux12C(i))  
       !     &                                   - PDB)  / PDB) * 1000.          
       !            ELSE  
       !               d13CFlux(i) =  d13Cresp(i)                 
       !            END IF     

       flux_turb(i) = Flux13C(i) + Flux12C(i)

       !        d13cflx_turb(i) = d13CFlux(i) * flux_turb(i)


    END DO

    DO i = 1, len
       ga(i) = ga_temp(i)  ! restore ga to its original value
    END DO

    RETURN
  END SUBROUTINE cfrax
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
          itable = INT(atable,kind=i4)
          atable = atable - itable

          TT_p1 = MAX(T_LOW,T(I+1))
          TT_p1 = MIN(T_HIGH,TT_p1)
          ATABLE_p1 = (TT_p1 - T_LOW + DELTA_T) * R_DELTA_T
          ITABLE_p1 = INT(ATABLE_p1,kind=i4)
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
    END DO  ! Npnts_do_1


    II = I
    DO I = II, NPNTS

          fsubw = 1.0_r8 + 1.0E-8_r8*P(I)*( 4.5_r8 +                               &
               &    6.0E-4_r8*( T(I) - zerodegC )*( T(I) - zerodegC ) )
          !
          TT = MAX(T_low,T(I))
          TT = MIN(T_high,TT)
          atable = (TT - T_low + delta_T) * R_delta_T
          itable = INT(atable,kind=i4)
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
    END DO  ! Npnts_do_1

    RETURN
  END SUBROUTINE QSAT_mix
  ! ======================================================================

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION dqsat (t, q1)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE      
    REAL(KIND=r8), INTENT(IN   ) :: t
    REAL(KIND=r8), INTENT(IN   ) :: q1

    !
    !
    ! statement function dqsat is d(qsat)/dt, with t in deg k and q1
    ! in kg/kg (q1 is *saturation* specific humidity)
    !
    dqsat = desat(t) * q1 * (1.0_r8 + q1*(1.0_r8/0.622_r8 - 1.0_r8)) / &
         esat(t)
    !
    RETURN
  END  FUNCTION dqsat

  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION qsat(e1, p1)
    ! ---------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: e1 
    REAL(KIND=r8), INTENT(IN   ) :: p1 

    ! statement function qsat is saturation specific humidity,
    ! with svp e1 and ambient pressure p in n/m**2. impose an upper
    ! limit of 1 to avoid spurious values for very high svp
    ! and/or small p1
    !
    qsat = 0.622_r8 * e1 /  &
         MAX ( p1 - (1.0_r8 - 0.622_r8) * e1, 0.622_r8 * e1 )

  END  FUNCTION qsat
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION desat(t)
    ! ---------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: t 

    desat  = 100.0_r8*( cvmgt (csat0, dsat0, t.GE.273.16_r8)             &
         + tsatl(t)*(csat1 + tsatl(t)*(csat2 + tsatl(t)*(csat3             &
         + tsatl(t)*(csat4 + tsatl(t)*(csat5 + tsatl(t)* csat6)))))  &
         + tsati(t)*(dsat1 + tsati(t)*(dsat2 + tsati(t)*(dsat3             &
         + tsati(t)*(dsat4 + tsati(t)*(dsat5 + tsati(t)* dsat6)))))  &
         )
  END  FUNCTION desat



  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION esat(t)
    ! ---------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN   ) :: t
    
    esat=fpvs2es5(t)

 !   esat = 100.0_r8*(cvmgt (asat0, bsat0, t.GE.273.16_r8)             &
 !        + tsatl(t)*(asat1 + tsatl(t)*(asat2 + tsatl(t)*(asat3             &
 !        + tsatl(t)*(asat4 + tsatl(t)*(asat5 + tsatl(t)* asat6)))))  &
 !        + tsati(t)*(bsat1 + tsati(t)*(bsat2 + tsati(t)*(bsat3             &
 !        + tsati(t)*(bsat4 + tsati(t)*(bsat5 + tsati(t)* bsat6)))))  &
 !        )
  END  FUNCTION esat
  !
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION tsatl (t)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: t
    !
    tsatl = MIN (100.0_r8, MAX (t-273.16_r8, 0.0_r8))
    !
    RETURN
  END  FUNCTION tsatl

  REAL(KIND=r8) FUNCTION tsati (t)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    REAL(KIND=r8), INTENT(IN   ) :: t
    !
    tsati = MAX (-60.0_r8, MIN (t-273.16_r8, 0.0_r8))
    !
    RETURN
  END  FUNCTION tsati
  ! ---------------------------------------------------------------------
  REAL(KIND=r8) FUNCTION cvmgt (x,y,l)
    ! ---------------------------------------------------------------------
    !
    ! chooses between two things.  Used in canopy.f
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN   ) :: l
    REAL(KIND=r8), INTENT(IN   ) :: x
    REAL(KIND=r8), INTENT(IN   ) :: y
    !
    IF (l) THEN
       cvmgt = x
    ELSE
       cvmgt = y
    END IF
    !
    RETURN
  END  FUNCTION cvmgt

END MODULE SFC_SiB2
