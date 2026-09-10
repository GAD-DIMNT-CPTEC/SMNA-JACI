!
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.17 $
!
MODULE GridHistory

  USE Parallelism, ONLY: &
       myid,             &
       maxNodes

  USE Constants, ONLY: &
       ngfld,&
       grav,i8,r8

  USE InputOutput, ONLY: &
       cnvray

  USE Communications, ONLY: &
       Collect_Grid_His

  USE Options, ONLY: &
       yrl       , &
       monl      , &
       nfprt     , &
       nferr     , &
       nfctrl    , &
       fNameGHLoc, &
       fNameGHTable, &
       start       , &
       TruncLev

  USE Utils, ONLY: &
       tmstmp2   , &
       lati

  USE Sizes, ONLY: &
       myfirstlat, &
       mylastlat,  &
       myfirstlon, &
       mylastlon,  &
       gridmap,    &
       nodehasj


  USE IOLowLevel, ONLY: &
       WriteGrdHist,&
       WrTopoGrdHist


   IMPLICIT NONE
  SAVE       


  INCLUDE 'mpif.h'

  PRIVATE

  PUBLIC :: InitGridHistory
  PUBLIC :: StoreGridHistory
  PUBLIC :: StoreMaskedGridHistory
  PUBLIC :: WriteGridHistory
  PUBLIC :: TurnOnGridHistory
  PUBLIC :: TurnOffGridHistory
  PUBLIC :: IsGridHistoryOn
  PUBLIC :: WriteGridHistoryTopo

  INTERFACE StoreGridHistory
     MODULE PROCEDURE Store2D, Store1D, Store2DV, Store1DV
  END INTERFACE

  !------------------------------
  ! Module GridHistory
  !------------------------------
  !
  ! The GridHistory module  has the function of to extract
  ! atmospheric fields during the simulation with the GCM.
  !
  ! The available atmospheric fields to be extracted during
  ! the simulation had been placed on Type "AvailFields",
  ! that it contains variables that describe the characteristics
  ! of each field, described bellow:
  !
  ! FldsName => Name of the meteorological available fields
  ! SurfFlds => It defines if the variable is of surface
  !             or has other layers
  ! UnGraFld => Code of field units used at conversion
  ! tfrmgf   => Time frame grads field
  ! GradsFld => Name of the used meteorological fields in
  !             the visualization, with software grads
  ! LyrsFlds => It stores the number of levels for each
  !             meteorological field
  ! LocGrFld => It locates in the second dimension of the
  !             diagnostic variable "dignos", the initial
  !             position of each meteorological field
  ! IsReqFld => It indicates which meteorological fields is
  !             required to compose grid history point
  !
  TYPE AvailFields
     CHARACTER(LEN=40 ) :: FldsName ! field name                   FldsName
     LOGICAL            :: SurfFlds ! if surface field or not      LyrsFlds
     INTEGER            :: UnGraFld ! code of field units          UnGraFld
     CHARACTER(LEN=1)   :: tfrmgf   ! time frame grads field       tfrmgf
     CHARACTER(LEN=4)   :: GradsFld ! field name for grads output  GradsFld
     INTEGER            :: LyrsFlds ! number of vertical levels
     INTEGER            :: LocGrFld ! total number of verticals in
     LOGICAL            :: IsReqFld ! if an available field is required
  END TYPE AvailFields
  TYPE(AvailFields) , ALLOCATABLE :: GHAF(:)
  !------------------------------------
  ! Grid points to collect and Store Grid History
  !------------------------------------
  !
  ! Required fields are collected on selected grid points.
  ! Each selected grid point has a name, longitude and
  ! latitude.
  !
  ! Required fields at selected grid points are copied
  ! into a data structure for a single time step. The
  ! data structure is a rank two array; first dimension
  ! is addressed by the selected grid point number;
  ! second dimension is addressed by a combination of
  ! vertical level and required field number.
  !
  ! PtCty   - grid point name        (contain the city name)
  ! ptLon   - grid point longitude   (  0:360 deg.)    (private)
  ! ptLat   - grid point latitude    (-90:90  deg.)    (private)
  ! ptCoor  - grid point coordinates (character)       (private)
  ! dignos  - It stores for each grid history point the
  !           required meteorological fields
  ! InReFd  - it stores the index for the required meteorological fields
  !
  !
  TYPE GridHistPoint
     CHARACTER(LEN=40)              :: PtCty
     INTEGER                        :: ibLoc
     INTEGER                        :: jbLoc
     REAL(KIND=r8)                           :: ptLon
     REAL(KIND=r8)                           :: ptLat
     CHARACTER(LEN=11)              :: ptCoor
     REAL(KIND=r8)   , POINTER, DIMENSION(:) :: dignos
     INTEGER, POINTER, DIMENSION(:) :: InReFd
  END TYPE GridHistPoint
  TYPE(GridHistPoint) , ALLOCATABLE :: GPt(:)

  !
  ! Mapping grid points into selected grid history points is
  ! done by rank two array *MPt*, indexed by grid point
  ! and block number. Array values are zero if the grid point
  ! is not selected and the number of grid point history
  ! (first dimension index of dignos) otherwise.
  !
  TYPE MapGrHist
     INTEGER                                    :: SumPt
     INTEGER,        POINTER, DIMENSION(:)      :: AlcPt
  END TYPE MapGrHist
  TYPE(MapGrHist), ALLOCATABLE:: MPt(:,:)
  !
  ! Available Grid History fields are computed for a set of
  ! grid points concurrently. The set of computed grid points
  ! may have some or none selected grid points for Grid History.
  !
  ! Rank two array *DoGrH*, indexed by available field and
  ! block of grid point number, indicates if the available field
  ! is required and if the block of computed grid points has
  ! selected grid history points
  !
  LOGICAL, PUBLIC, ALLOCATABLE :: DoGrH(:,:)!(ngfld, jbMax)
  !
  ! Available Grid History Indexes
  !
  INTEGER, PUBLIC, PARAMETER :: nGHis_presfc =  1 ! surface pressure
  INTEGER, PUBLIC, PARAMETER :: nGHis_tcanop =  2 ! canopy temperature
  INTEGER, PUBLIC, PARAMETER :: nGHis_tgfccv =  3 ! ground/surface cover temperature
  INTEGER, PUBLIC, PARAMETER :: nGHis_tgdeep =  4 ! deep soil temperature
  INTEGER, PUBLIC, PARAMETER :: nGHis_swtsfz =  5 ! soil wetness of surface zone
  INTEGER, PUBLIC, PARAMETER :: nGHis_swtrtz =  6 ! soil wetness of root zone
  INTEGER, PUBLIC, PARAMETER :: nGHis_swtrcz =  7 ! soil wetness of recharge zone
  INTEGER, PUBLIC, PARAMETER :: nGHis_mostca =  8 ! moisture store on canopy
  INTEGER, PUBLIC, PARAMETER :: nGHis_mostgc =  9 ! moisture store on ground cover
  INTEGER, PUBLIC, PARAMETER :: nGHis_snowdp = 10 ! snow depth
  INTEGER, PUBLIC, PARAMETER :: nGHis_snowfl = 11 ! snowfall
  INTEGER, PUBLIC, PARAMETER :: nGHis_rouglg = 12 ! roughness length
  INTEGER, PUBLIC, PARAMETER :: nGHis_ustres = 13 ! surface zonal wind stress
  INTEGER, PUBLIC, PARAMETER :: nGHis_vstres = 14 ! surface meridional wind stress
  INTEGER, PUBLIC, PARAMETER :: nGHis_sheatf = 15 ! sensible heat flux from surface
  INTEGER, PUBLIC, PARAMETER :: nGHis_lheatf = 16 ! latent heat flux from surface
  INTEGER, PUBLIC, PARAMETER :: nGHis_toprec = 17 ! total precipitation
  INTEGER, PUBLIC, PARAMETER :: nGHis_cvprec = 18 ! convective precipitation
  INTEGER, PUBLIC, PARAMETER :: nGHis_swdtop = 19 ! incident short wave flux
  INTEGER, PUBLIC, PARAMETER :: nGHis_lwutop = 20 ! outgoing long wave at top
  INTEGER, PUBLIC, PARAMETER :: nGHis_lwdbot = 21 ! downward long wave at ground
  INTEGER, PUBLIC, PARAMETER :: nGHis_lwubot = 22 ! upward long wave flux at ground
  INTEGER, PUBLIC, PARAMETER :: nGHis_swutop = 23 ! upward short wave at top
  INTEGER, PUBLIC, PARAMETER :: nGHis_swdbvb = 24 ! downward short wave flux at ground (vb)
  INTEGER, PUBLIC, PARAMETER :: nGHis_swdbvd = 25 ! downward short wave flux at ground (vd)
  INTEGER, PUBLIC, PARAMETER :: nGHis_swdbnb = 26 ! downward short wave flux at ground (nb)
  INTEGER, PUBLIC, PARAMETER :: nGHis_swdbnd = 27 ! downward short wave flux at ground (nd)
  INTEGER, PUBLIC, PARAMETER :: nGHis_vibalb = 28 ! visible beam albedo
  INTEGER, PUBLIC, PARAMETER :: nGHis_vidalb = 29 ! visible diffuse albedo
  INTEGER, PUBLIC, PARAMETER :: nGHis_nibalb = 30 ! near infrared beam albedo
  INTEGER, PUBLIC, PARAMETER :: nGHis_nidalb = 31 ! near infrared diffuse albedo
  INTEGER, PUBLIC, PARAMETER :: nGHis_vegtyp = 32 ! vegetation type
  INTEGER, PUBLIC, PARAMETER :: nGHis_nrdcan = 33 ! net radiation of canopy
  INTEGER, PUBLIC, PARAMETER :: nGHis_nrdgsc = 34 ! net radiation of ground surface/cover
  INTEGER, PUBLIC, PARAMETER :: nGHis_coszen = 35 ! cosine of zenith angle
  INTEGER, PUBLIC, PARAMETER :: nGHis_dragcf = 36 ! drag
  INTEGER, PUBLIC, PARAMETER :: nGHis_mofres = 37 ! momentum flux resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_casrrs = 38 ! canopy air spc to ref. lvl resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_cascrs = 39 ! canopy air spc to canopy resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_casgrs = 40 ! canopy air spc to ground resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_canres = 41 ! canopy resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_gcovrs = 42 ! ground cover resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_bssfrs = 43 ! bare soil surface resistance
  INTEGER, PUBLIC, PARAMETER :: nGHis_ecairs = 44 ! vapor pressure of canopy air space
  INTEGER, PUBLIC, PARAMETER :: nGHis_tcairs = 45 ! temperature of canopy air space
  INTEGER, PUBLIC, PARAMETER :: nGHis_tracan = 46 ! transpiration from canopy
  INTEGER, PUBLIC, PARAMETER :: nGHis_inlocp = 47 ! interception loss from canopy
  INTEGER, PUBLIC, PARAMETER :: nGHis_tragcv = 48 ! transpiration from ground cover
  INTEGER, PUBLIC, PARAMETER :: nGHis_inlogc = 49 ! interception loss from ground cover
  INTEGER, PUBLIC, PARAMETER :: nGHis_bsevap = 50 ! bare soil evaporation
  INTEGER, PUBLIC, PARAMETER :: nGHis_shfcan = 51 ! sensible heat flux from canopy
  INTEGER, PUBLIC, PARAMETER :: nGHis_shfgnd = 52 ! sensible heat flux from ground
  INTEGER, PUBLIC, PARAMETER :: nGHis_canhea = 53 ! canopy heating rate
  INTEGER, PUBLIC, PARAMETER :: nGHis_gcheat = 54 ! ground/surface cover heating rate
  INTEGER, PUBLIC, PARAMETER :: nGHis_runoff = 55 ! runoff
  INTEGER, PUBLIC, PARAMETER :: nGHis_hcseai = 56 ! heat conduction through sea ice
  INTEGER, PUBLIC, PARAMETER :: nGHis_hsseai = 57 ! heat storage tendency over sea ice
  INTEGER, PUBLIC, PARAMETER :: nGHis_uzonal = 58 ! zonal wind (u)
  INTEGER, PUBLIC, PARAMETER :: nGHis_vmerid = 59 ! meridional wind (v)
  INTEGER, PUBLIC, PARAMETER :: nGHis_temper = 60 ! virtual temperature
  INTEGER, PUBLIC, PARAMETER :: nGHis_spchum = 61 ! specific humidity
  INTEGER, PUBLIC, PARAMETER :: nGHis_swheat = 62 ! short wave radiative heating
  INTEGER, PUBLIC, PARAMETER :: nGHis_lwheat = 63 ! long wave radiative heating
  INTEGER, PUBLIC, PARAMETER :: nGHis_sslaht = 64 ! supersaturation latent heating
  INTEGER, PUBLIC, PARAMETER :: nGHis_clheat = 65 ! convective latent heating
  INTEGER, PUBLIC, PARAMETER :: nGHis_sclhea = 66 ! shallow convective heating
  INTEGER, PUBLIC, PARAMETER :: nGHis_vdheat = 67 ! vertical diffusion heating
  INTEGER, PUBLIC, PARAMETER :: nGHis_spstms = 68 ! supersaturation moisture source
  INTEGER, PUBLIC, PARAMETER :: nGHis_cvmosr = 69 ! convective moisture source
  INTEGER, PUBLIC, PARAMETER :: nGHis_shcvmo = 70 ! shallow convective moistening
  INTEGER, PUBLIC, PARAMETER :: nGHis_vdmois = 71 ! vertical diffusion moistening
  INTEGER, PUBLIC, PARAMETER :: nGHis_vduzon = 72 ! vertical diffusion du/dt
  INTEGER, PUBLIC, PARAMETER :: nGHis_vdvmer = 73 ! vertical diffusion dv/dt
  INTEGER, PUBLIC, PARAMETER :: nGHis_cloudc = 74 ! cloud cover
  INTEGER, PUBLIC, PARAMETER :: nGHis_vdtclc = 75 ! vertical dist total cloud cover
  INTEGER, PUBLIC, PARAMETER :: nGHis_uzonsf = 76 ! surface zonal wind (u)
  INTEGER, PUBLIC, PARAMETER :: nGHis_vmersf = 77 ! surface meridional wind (v)
  INTEGER, PUBLIC, PARAMETER :: nGHis_tvirsf = 78 ! surface virtual temperature
  INTEGER, PUBLIC, PARAMETER :: nGHis_sphusf = 79 ! surface specific humidity
  INTEGER, PUBLIC, PARAMETER :: nGHis_tep02m = 80 ! temp at 2-m from sfc layer
  INTEGER, PUBLIC, PARAMETER :: nGHis_mxr02m = 81 ! especific humid at 2-m from sfc layer
  INTEGER, PUBLIC, PARAMETER :: nGHis_zwn10m = 82 ! Zonal Wind at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nGHis_mwn10m = 83 ! Meridional wind at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nGHis_swdgrd = 84 ! shortwave downward at ground
  INTEGER, PUBLIC, PARAMETER :: nGHis_swugrd = 85 ! shortwave upward at bottom
  INTEGER, PUBLIC, PARAMETER :: nGHis_tsgrnd = 86 ! Temperatura da superficie do solo  (K)
  INTEGER, PUBLIC, PARAMETER :: nGHis_tspres = 87 ! Tendency Surface Pressure (Pa/s)
  INTEGER, PUBLIC, PARAMETER :: nGHis_cldlow = 88 ! fraction of clouds for Low (%)
  INTEGER, PUBLIC, PARAMETER :: nGHis_cldmed = 89 ! fraction of clouds for medium (%)
  INTEGER, PUBLIC, PARAMETER :: nGHis_cldHig = 90 ! fraction of clouds for High (%)
  !
  !
  ! it is possible to select all available grid history
  ! fields:
  !
  ! allghf - if all available fields are required    (private)
  !
  LOGICAL :: allghf
  !
  !----------------------------
  ! Turning Grid History On/Off
  !----------------------------
  !
  ! The input flag *grhflg* defines if this run will have
  ! grid history turned on or off. (private)
  !
  LOGICAL :: grhflg
  !
  ! Whenever *grhflg* is turned on, there are timesteps
  ! where Grid History will be collected, and there are
  ! timesteps where grid history will not be collected.
  ! These are controlled by the (private) variable
  ! *grhOn*, that is turned on/off by TurnOnGridHistory
  ! and TurnOffGridHistory. It is inquired by IsGridHistoryOn.
  ! Default is GridHistory turned off.
  !
  LOGICAL              :: grhOn

  INTEGER              :: ngrfld    ! Number required fields                       (private)
  INTEGER              :: nghsl     ! total number of verticals in *dignos*        (private)
  INTEGER              :: ngpts     ! There are *ngpts* selected points            (private)
  INTEGER              :: ngptslocal
  INTEGER, ALLOCATABLE :: gptslocal(:)
  INTEGER, ALLOCATABLE :: ngptsperjb(:)
  INTEGER, ALLOCATABLE :: map(:)
  INTEGER, ALLOCATABLE :: mapGlobal(:)
  INTEGER, ALLOCATABLE :: procmap(:)
  INTEGER, ALLOCATABLE :: iniperjb(:)
  INTEGER, ALLOCATABLE :: ngptsperproc(:)
  INTEGER              :: ngpts_new

  INTEGER              :: iMax
  INTEGER              :: jMax
  INTEGER              :: ibMax
  INTEGER              :: jbMax
  INTEGER,ALLOCATABLE  :: ibMaxPerJB(:)
  INTEGER,ALLOCATABLE  :: iMaxPerJ  (:)
  INTEGER, ALLOCATABLE :: ibPerIJ (:,:)
  INTEGER, ALLOCATABLE :: jbPerIJ (:,:)
  INTEGER              :: kMax

  CHARACTER(LEN=40), ALLOCATABLE :: PtCtyGlobal(:)
  CHARACTER(LEN=18), ALLOCATABLE :: PtCoorGlobal(:)
  REAL(KIND=r8)    , ALLOCATABLE :: PtLonGlobal(:)
  REAL(KIND=r8)    , ALLOCATABLE :: PtLatGlobal(:)
  INTEGER           :: jbLocMin(1)
  INTEGER           :: IbLocMin(1)

CONTAINS
  SUBROUTINE InitGridHistFields
    ALLOCATE(GHAF(ngfld))
    ! field name

    GHAF(1:ngfld)%FldsName= (/ "SURFACE PRESSURE                        ", &
         "CANOPY TEMPERATURE                      ", "GROUND/SURFACE COVER TEMPERATURE        ", &
         "DEEP SOIL TEMPERATURE                   ", "SOIL WETNESS OF SURFACE ZONE            ", &
         "SOIL WETNESS OF ROOT ZONE               ", "SOIL WETNESS OF RECHARGE ZONE           ", &
         "MOISTURE STORE ON CANOPY                ", "MOISTURE STORE ON GROUND COVER          ", &
         "SNOW DEPTH                              ", "SNOWFALL                                ", &
         "ROUGHNESS LENGTH                        ", "SURFACE ZONAL WIND STRESS               ", &
         "SURFACE MERIDIONAL WIND STRESS          ", "SENSIBLE HEAT FLUX FROM SURFACE         ", &
         "LATENT HEAT FLUX FROM SURFACE           ", "TOTAL PRECIPITATION                     ", &
         "CONVECTIVE PRECIPITATION                ", "INCIDENT SHORT WAVE FLUX                ", &
         "OUTGOING LONG WAVE AT TOP               ", "DOWNWARD LONG WAVE AT GROUND            ", &
         "UPWARD LONG WAVE FLUX AT GROUND         ", "UPWARD SHORT WAVE AT TOP                ", &
         "DOWNWARD SHORT WAVE FLUX AT GROUND (VB) ", "DOWNWARD SHORT WAVE FLUX AT GROUND (VD) ", &
         "DOWNWARD SHORT WAVE FLUX AT GROUND (NB) ", "DOWNWARD SHORT WAVE FLUX AT GROUND (ND) ", &
         "VISIBLE BEAM ALBEDO                     ", "VISIBLE DIFFUSE ALBEDO                  ", &
         "NEAR INFRARED BEAM ALBEDO               ", "NEAR INFRARED DIFFUSE ALBEDO            ", &
         "VEGETATION TYPE                         ", "NET RADIATION OF CANOPY                 ", &
         "NET RADIATION OF GROUND SURFACE/COVER   ", "COSINE OF ZENITH ANGLE                  ", &
         "DRAG                                    ", "MOMENTUM FLUX RESISTANCE                ", &
         "CANOPY AIR SPC TO REF. LVL RESISTANCE   ", "CANOPY AIR SPC TO CANOPY RESISTANCE     ", &
         "CANOPY AIR SPC TO GROUND RESISTANCE     ", "CANOPY RESISTANCE                       ", &
         "GROUND COVER RESISTANCE                 ", "BARE SOIL SURFACE RESISTANCE            ", &
         "VAPOR PRESSURE OF CANOPY AIR SPACE      ", "TEMPERATURE OF CANOPY AIR SPACE         ", &
         "TRANSPIRATION FROM CANOPY               ", "INTERCEPTION LOSS FROM CANOPY           ", &
         "TRANSPIRATION FROM GROUND COVER         ", "INTERCEPTION LOSS FROM GROUND COVER     ", &
         "BARE SOIL EVAPORATION                   ", "SENSIBLE HEAT FLUX FROM CANOPY          ", &
         "SENSIBLE HEAT FLUX FROM GROUND          ", "CANOPY HEATING RATE                     ", &
         "GROUND/SURFACE COVER HEATING RATE       ", "RUNOFF                                  ", &
         "HEAT CONDUCTION THROUGH SEA ICE         ", "HEAT STORAGE TENDENCY OVER SEA ICE      ", &
         "ZONAL WIND (U)                          ", "MERIDIONAL WIND (V)                     ", &
         "VIRTUAL TEMPERATURE                     ", "SPECIFIC HUMIDITY                       ", &
         "SHORT WAVE RADIATIVE HEATING            ", "LONG WAVE RADIATIVE HEATING             ", &
         "SUPERSATURATION LATENT HEATING          ", "CONVECTIVE LATENT HEATING               ", &
         "SHALLOW CONVECTIVE HEATING              ", "VERTICAL DIFFUSION HEATING              ", &
         "SUPERSATURATION MOISTURE SOURCE         ", "CONVECTIVE MOISTURE SOURCE              ", &
         "SHALLOW CONVECTIVE MOISTENING           ", "VERTICAL DIFFUSION MOISTENING           ", &
         "VERTICAL DIFFUSION DU/DT                ", "VERTICAL DIFFUSION DV/DT                ", &
         "CLOUD COVER                             ", "VERTICAL DIST TOTAL CLOUD COVER         ", &
         "SURFACE ZONAL WIND (U)                  ", "SURFACE MERIDIONAL WIND (V)             ", &
         "SURFACE VIRTUAL TEMPERATURE             ", "SURFACE SPECIFIC HUMIDITY               ", &
         "TEMP AT 2-M FROM SFC                    ", "SPECIFIC HUMID AT 2-M FROM SFC          ", &
         "ZONAL WIND AT 10-M FROM SFC             ", "MERIDIONAL WIND AT 10-M FROM SFC        ", &
         "SHORTWAVE DOWNWARD AT GROUND            ", "SHORTWAVE UPWARD AT BOTTOM              ", &
         "TEMPERATURA DA SUPERFICIE DO SOLO       ", "TENDENCY SURFACE PRESSURE               ", &
         "FRACTION OF CLOUDS FOR LOW              ", "FRACTION OF CLOUDS FOR MEDIUM           ", &
         "FRACTION OF CLOUDS FOR HIGH             "/)


    ! surface field or not

    GHAF(1:ngfld)%SurfFlds=    (/ &
         .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
         .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
         .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
         .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
         .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
         .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
         .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .TRUE. , .FALSE., .TRUE. , .TRUE. , &
         .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE., &
         .TRUE. , .TRUE./)

    ! units

    GHAF(1:ngfld)%UnGraFld=  (/ &
         131,  40,  40,  40,   0,   0,   0, 110, 110, 110, &
         120,  10, 130, 130, 170, 170, 120, 120, 170, 170, &
         170, 170, 170, 170, 170, 170, 170,   0,   0,   0, &
           0,   0, 170, 170,   0, 200, 190, 190, 190, 190, &
         190, 190, 190, 131,  40, 170, 170, 170, 170, 170, &
         170, 170, 170, 170, 120, 170, 170,  60,  60,  40, &
           0,  70,  70,  70,  70,  70, 70 ,  50,  50,  50, &
          50, 100, 100,   0,   0,  60,  60,  40,   0,  40, &
           0,  60,  60, 170, 170,  40, 150,   0,   0,   0/)

    ! time frame

    GHAF(1:ngfld)%tfrmgf =     (/ &
         "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "C", &
         "I", "C", "C", "P", "P", "C", "C", "C", "C", "C", "C", &
         "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", &
         "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "P", &
         "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", &
         "C", "C", "I", "I", "I", "I", "C", "C", "C", "C", "C", &
         "C", "C", "C", "C", "C", "C", "C", "C", "C", "I", "I", &
         "I", "I", "I", "I", "I", "I", "C", "C", "I", "I", "I", &
         "I", "I"/)

    ! grads names

    GHAF(1:ngfld)%GradsFld =        (/&
         "PSLC", "TDSL", "TGSC", "TGDP", "USSL", "UZRS", "UZDS", "AUDL", "AUGC", "PNEV", "NEVE", &
         "ZORL", "USST", "VSST", "CSSF", "CLSF", "PREC", "PRCV", "FOCI", "ROLE", "OLIS", "OLES", &
         "ROCE", "CIVB", "CIVD", "CINB", "CIND", "ALVB", "ALVD", "ALIB", "ALID", "TIVG", "SRDL", &
         "SRGC", "CAZL", "DRAG", "RFLM", "RDNR", "RDAD", "RDAG", "RDSL", "RSGC", "RSBS", "PVDL", &
         "TADL", "TRDL", "PIDL", "TRGC", "PIGC", "EVBS", "CSDL", "CSGR", "TAQD", "TAQG", "RNOF", &
         "CAGM", "ACGM", "UVEL", "VVEL", "TEMV", "UMES", "AROC", "AROL", "CLSS", "CLCV", "ACVR", &
         "DVAQ", "FUSS", "FUCV", "UCVR", "DVUM", "DVTU", "DVTV", "CBNV", "VDCC", "UVES", "VVES", &
         "TEVS", "UESS", "TP2M", "QQ2M", "US2M", "VS2M", "OCIS", "OCES", "TGRD", "TNPS", "CLLW", &
         "CLMD", "CLHI"/)



  END SUBROUTINE InitGridHistFields



  SUBROUTINE InitGridHistory (idate, idatec,iov, &
       allghf_in, grhflg_in, nfghds, nfghloc, nfghdr,iMax_in,jMax_in,&
       ibMax_in,jbMax_in,ibPerIJ_in,jbPerIJ_in,kMax_in,ibMaxPerJB_in,iMaxPerJ_in,start)
    INTEGER,           INTENT(in) :: iMax_in
    INTEGER,           INTENT(in) :: jMax_in
    INTEGER,           INTENT(in) :: ibMax_in
    INTEGER,           INTENT(in) :: jbMax_in
    INTEGER,           INTENT(in) :: iov
    INTEGER,           INTENT(in) :: ibPerIJ_in(1-iov:iMax_in+iov,-1:jMax_in+2)
    INTEGER,           INTENT(in) :: jbPerIJ_in(1-iov:iMax_in+iov,-1:jMax_in+2)
    INTEGER,           INTENT(in) :: kMax_in
    INTEGER,           INTENT(in) :: idate(4)
    INTEGER,           INTENT(in) :: idatec(4)
    LOGICAL,           INTENT(in) :: allghf_in
    LOGICAL,           INTENT(in) :: grhflg_in
    INTEGER,           INTENT(in) :: nfghds
    INTEGER,           INTENT(in) :: nfghloc
    INTEGER,           INTENT(in) :: nfghdr
    INTEGER,           INTENT(in) :: ibMaxPerJB_in(:)
    INTEGER,           INTENT(in) :: iMaxPerJ_in(:)
    CHARACTER(len=*) , INTENT(in) ::  start

    REAL(KIND=r8)                          :: Lat
    INTEGER, ALLOCATABLE          :: jslocal(:)
    INTEGER, ALLOCATABLE          :: islocal(:)

    CHARACTER(len= *), PARAMETER :: h="**(InitGridHistory)**"
    CHARACTER(len=20), PARAMETER :: typgh='GRID POINT HISTORY  '
    CHARACTER(len= 4), PARAMETER :: iacc='SEQU'
    CHARACTER(len= 4), PARAMETER :: idev='TAPE'
    CHARACTER(LEN= 4), PARAMETER :: nexp='0003'
    CHARACTER(LEN=40)             :: rdesc
    INTEGER                       :: UniGraFld

    REAL(KIND=r8)    :: pi
    INTEGER :: n, nloc
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: ib
    INTEGER :: jb
    INTEGER :: nn
    INTEGER :: ierr
    LOGICAL :: notfound
    REAL(KIND=r8)              :: dlon   (iMax_in)
    REAL(KIND=r8)              :: dlat   (jMax_in)

    CALL InitGridHistFields

    pi = 4.0_r8 * ATAN(1.0_r8)
    ALLOCATE (ibPerIJ(iMax_in ,jMax_in ))
    ibPerIJ=-1
    ALLOCATE (jbPerIJ(iMax_in ,jMax_in ))
    jbPerIJ=-1
    ALLOCATE (ibMaxPerJB(jbMax_in))
    ibMaxPerJB=-1
    ALLOCATE (iMaxPerJ(jMax_in))
    iMaxPerJ = -1

    ! store input data

    allghf = allghf_in
    grhflg = grhflg_in
    grhOn  = .FALSE.

    iMax    = iMax_in
    jMax    = jMax_in
    ibMax   = ibMax_in
    jbMax   = jbMax_in
    ibPerIJ = ibPerIJ_in(1:iMax,1:jMax)
    jbPerIJ = jbPerIJ_in(1:iMax,1:jMax)
    iMaxPerJ= iMaxPerJ_in
    kMax    = kMax_in
    ibMaxPerJB=ibMaxPerJB_in

    ! if grid history not required, fix data structure and return

    IF (.NOT. grhflg) THEN
       ALLOCATE (DoGrH(ngfld, jbMax))
       DoGrH = .FALSE.
       RETURN
    END IF

    ! # vertical levels of available grid history fields

    DO n = 1, ngfld
       IF (GHAF(n)%SurfFlds) THEN
          GHAF(n)%LyrsFlds = 1
       ELSE
          GHAF(n)%LyrsFlds = kMax
       END IF
    END DO

    ! dump available fields

    IF (nfctrl(55) >= 1) THEN
       WRITE(UNIT=nfprt,FMT=110)
       DO n = 1, ngfld
          WRITE(UNIT=nfprt,FMT=140) n,  GHAF(n)%FldsName, &
               GHAF(n)%LyrsFlds, &
               GHAF(n)%UnGraFld, &
               GHAF(n)%tfrmgf  , &
               GHAF(n)%GradsFld
       END DO
    END IF

    ! # required fields

    OPEN(UNIT=nfghds,FILE=TRIM(fNameGHTable),FORM='formatted',ACCESS='sequential',&
         ACTION='read',STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameGHTable), ierr
       STOP "**(ERROR)**"
    END IF

    IF (allghf) THEN
       ngrfld = ngfld
    ELSE IF (nfghds  > 0) THEN
       ngrfld = 0
       DO
          READ (UNIT=nfghds, FMT=225, IOSTAT=ierr)
          IF (ierr == 0) THEN
             ngrfld = ngrfld + 1
          ELSE IF (ierr > 0) THEN
             WRITE(UNIT=nfprt, FMT="(a,' error reading unit ',i4)") h, nfghds
             WRITE(UNIT=nferr, FMT="(a,' error reading unit ',i4)") h, nfghds
             STOP h
          ELSE
             EXIT
          END IF
       END DO
    END IF

    IF (ngrfld == 0 .AND. .NOT. allghf) THEN
       WRITE(UNIT=nfprt,FMT=2135)
       WRITE(UNIT=nferr,FMT=2135)
       STOP h
    END IF

    ! initialize required field data structure


    IF (allghf) THEN
       DO n=1,ngfld
          GHAF(n)%IsReqFld = .TRUE.
       END DO
    ELSE
       GHAF(1:ngfld)%IsReqFld = .FALSE.
       REWIND nfghds
       IF (nfctrl(55) > 0) WRITE(UNIT=nfprt,FMT=210)
       DO n = 1, ngrfld
          READ (UNIT=nfghds, FMT=225) rdesc, UniGraFld
          IF (nfctrl(55) > 0) WRITE(UNIT=nfprt,FMT=240) n, rdesc, UniGraFld
          notfound = .TRUE.
          DO nn = 1, ngfld
             IF ( rdesc == GHAF(nn)%FldsName) THEN
                GHAF(nn)%IsReqFld   = .TRUE.
                Notfound  = .FALSE.
                EXIT
             END IF
          END DO
          IF (notfound) THEN
             WRITE(UNIT=nfprt, FMT="(a,' required field ',a,' not available')") h, TRIM(rdesc)
             WRITE(UNIT=nferr, FMT="(a,' required field ',a,' not available')") h, TRIM(rdesc)
             STOP h
          END IF
       END DO
    END IF

    ! where to store each required field and total # of verticals

    nghsl = 0
    DO n = 1, ngfld
       IF (GHAF(n)%IsReqFld) THEN
          GHAF(n)%LocGrFld = nghsl + 1
          nghsl    = nghsl + GHAF(n)%LyrsFlds
       ELSE
          GHAF(n)%LocGrFld = 0
       END IF
    END DO

    ! # grid points to collect grid history

    OPEN(UNIT=nfghloc,FILE=TRIM(fNameGHLoc),FORM='formatted',ACCESS='sequential',&
         ACTION='read',STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameGHLoc), ierr
       STOP "**(ERROR)**"
    END IF
    ngpts = 0
    DO
       READ (UNIT=nfghloc, FMT='(A40,2F11.5,1X,A11)', IOSTAT=ierr)
       IF (ierr == 0) THEN
          ngpts = ngpts + 1
       ELSE IF (ierr > 0) THEN
          WRITE(UNIT=nfprt, FMT="(a,' error reading unit ',i4)") h, nfghloc
          WRITE(UNIT=nferr, FMT="(a,' error reading unit ',i4)") h, nfghloc
          STOP h
       ELSE
          EXIT
       END IF
    END DO
    !
    !  detect which are the points to handle in this process
    !
    ALLOCATE(gptslocal(ngpts))
    ALLOCATE(ngptsperproc(0:maxnodes-1))
    ALLOCATE(procmap(ngpts))
    ALLOCATE(jslocal(ngpts))
    ALLOCATE(islocal(ngpts))
    ALLOCATE(PtCtyGlobal(ngpts))
    ALLOCATE(PtCoorGlobal(ngpts))
    ALLOCATE(PtLonGlobal(ngpts))
    ALLOCATE(PtLatGlobal(ngpts))

    ! fill grid point data structure

    REWIND nfghloc
    DO n = 1, ngpts
       READ (nfghloc,*) PtCtyGlobal(n),PtLonGlobal(n),PtLatGlobal(n)
    END DO

    CALL GridHistLabel (ngpts, iMax,jMax,iMaxPerJ,PtLonGlobal, PtLatGlobal, PtCoorGlobal)
    dlon=0.0_r8
    dlat=0.0_r8
    ngptslocal = 0
    ngptsperproc = 0
    DO n = 1, ngpts
       Lat=-(PtLatGlobal (n) - 90.0_r8)*(pi/180.0_r8)
       DO j=1,jMax
          dlat(j) = ABS(lati(j) - Lat)
       END DO
       jblocmin    = MINLOC (dlat)
       j           = jblocmin(1)
       DO i=1,iMaxPerJ(j)
          dlon(i) = ABS(((pi*((i-1)*360.0_r8/REAL(iMaxPerJ(j),r8)))/180.0_r8) &
               &         -((pi*PtLonGlobal(n))/180.0_r8))
       END DO
       IbLocMin=MINLOC (dlon(1:iMaxPerJ(j)))
       i = IbLocMin(1)
       procmap(n) = gridmap(i,j)
       ngptsperproc(gridmap(i,j)) = ngptsperproc(gridmap(i,j)) + 1

       IF (gridmap(i,j).eq.myid) THEN
          ngptslocal            = ngptslocal + 1
          gptslocal(ngptslocal) = n
          jslocal(ngptslocal)   = j
          islocal(ngptslocal)   = i
       ENDIF
    END DO
    ALLOCATE(mapglobal(ngpts))
    i = 1
    DO k=0,maxnodes-1
       DO n=1,ngpts
          IF(procmap(n).eq.k) THEN
             mapglobal(i) = n
             i = i + 1
          ENDIF
       ENDDO
    ENDDO

    ! sort grid points per latitude

    ALLOCATE(map(ngptslocal))
    nloc = 0
    DO j=myfirstlat,mylastlat
       DO i=myfirstlon(j),mylastlon(j)
          DO n=1,ngptslocal
             IF (jslocal(n).eq.j.and.islocal(n).eq.i) THEN
                nloc = nloc + 1
                !compacta os pontos no array local
                map(n)              = nloc
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! allocate grid point data structure

    ALLOCATE(GPt(ngptslocal))
    ALLOCATE(ngptsperjb(jbmax))
    ALLOCATE(iniperjb(jbmax))

    ! fill grid point data structure

    DO n=1,ngptslocal
       nloc = map(n)
       GPt(nloc)%PtCty  = PtCtyGlobal (gptslocal(n))
       GPt(nloc)%ptLon  = PtLonGlobal (gptslocal(n))
       GPt(nloc)%ptLat  =-(PtLatGlobal(gptslocal(n)) - 90.0_r8)*(pi/180.0_r8)! (0-pi)--(N-S)
       GPt(nloc)%ptCoor = PtCoorGlobal(gptslocal(n))(1:11)
    ENDDO


    IF (nfctrl(55) > 0) THEN
       WRITE(UNIT=nfprt,FMT=520)
       DO nloc=1,ngptslocal
          n = map(nloc)
          WRITE (UNIT=nfprt, FMT='(1X,2I4,1X,A40,2F11.5,1X,A11)')&
               &        myid,n, GPt(n)%PtCty, GPt(n)%ptLon, &
               &                90.0_r8- (180.0_r8/pi)*GPt(n)%ptLat, GPt(n)%ptCoor
       ENDDO
    END IF

    ! allocate and initialize grid history buffer

    DO n = 1, ngptslocal
       CALL NULLIFY_dignos(GPt( n ),nghsl  )
       CALL NULLIFY_InReFd(GPt( n ),ngrfld )
    END DO

    DO nloc = 1, ngptslocal
       n = map(nloc)
       ib=0
       DO nn = 1, ngfld
          IF (  GHAF(nn)%IsReqFld ) THEN
             ib=ib+1
             GPt( n )%InReFd(ib) = nn
          END IF
       END DO
    END DO
    !
    ! computes the indices for grid history points
    !
    ngptsperjb = 0
    DO nloc = 1, ngptslocal
       n = map(nloc)
       GPt(n)%ibLoc =  ibPerIJ(  islocal(nloc) , jslocal(nloc)  )
       GPt(n)%jbLoc =  jbPerIJ(  islocal(nloc) , jslocal(nloc)  )
       ngptsperjb( GPt(n)%jbLoc ) = ngptsperjb( GPt(n)%jbLoc ) + 1
    END DO

    iniperjb(1)=0
    DO jb = 2, jbMax
       iniperjb(jb) = iniperjb(jb-1)+ngptsperjb(jb-1)
    END DO

    ! allocate and initialize where to do grid history

    ALLOCATE (DoGrH(ngfld, jbMax))
    DO jb = 1, jbMax
       IF (ngptsperjb(jb).gt.0) THEN
          DoGrH(:, jb) = GHAF(:)%IsReqFld
       ELSE
          DoGrH(:, jb) = .FALSE.
       END IF
    END DO
    IF(myid ==0 ) THEN
       ! dumping for the file the information of the  grid history points
       WRITE(UNIT=nfghdr,FMT=700)typgh
       IF( TRIM(start) /= "warm" )THEN
          WRITE(UNIT=nfghdr,FMT=720) nexp, iacc, ibMax, jbMax, kMax, kMax, ngpts, nghsl, nghsl, &
               idate, idev
       ELSE
           WRITE(UNIT=nfghdr,FMT=720) nexp, iacc, ibMax, jbMax, kMax, kMax, ngpts, nghsl, nghsl, &
               idatec, idev      
       END IF
       WRITE(UNIT=nfghdr,FMT=710) 'CPTEC AGCM R1.2 2001  '//TRIM(TruncLev)//'  '//TRIM(start)
!      WRITE(UNIT=nfghdr,FMT=740)(del(k),k=1,kMax)
       DO n=1,ngfld
          IF(GHAF(n)%IsReqFld )THEN
             WRITE(UNIT=nfghdr,FMT=126) GHAF(n)%FldsName, &
                  &                GHAF(n)%LyrsFlds, &
                  &                GHAF(n)%UnGraFld, &
                  &                GHAF(n)%GradsFld
          END IF
       END DO
       DO n=1,ngpts
          WRITE(UNIT=nfghdr,FMT=550) PtCtyGlobal(n),PtCoorGlobal(n)
       ENDDO
       ENDFILE nfghdr
    END IF
110 FORMAT(' NO.',T10,'AVAILABLE GRID HISTORY FIELD DESCRIPTION', &
         T47,'NUM. OF LAYERS',T66,'UNITS',T75,'TIME FRAME')
126 FORMAT(A40,I5,2X,I5,1X,A4)
140 FORMAT(' ',I4,' ',A40,I5,2X,I5,2X,A1,1X,A4)
210 FORMAT(' NO.',T10,'REQUESTED GRID HISTORY FIELD DESCRIPTION', &
         T55,'REQUESTED UNITS')
225 FORMAT(A40,I5)
240 FORMAT(' ',I4,' ',A40,I5)
520 FORMAT(' NO.',T17,'POINT DESCRIPTION',T51,'I PT.',T59,'J PT.')
550 FORMAT(A40,1X,A11)
!570 FORMAT(' ',I4,' ',A40,I5,I5,1X,A11)
700 FORMAT(A20)
710 FORMAT(A60)
720 FORMAT(A4,1X,A4,11I5,1X,A4)
740 FORMAT(5E16.8)
2135 FORMAT(' REQUESTED GRID POINT HISTORY FIELD TABLE EMPTY OR', &
         ' NOT FOUND'/' WITH ALLGHF=F')
  END SUBROUTINE InitGridHistory

  SUBROUTINE GridHistLabel (Nmax,iMax,jMax,iMaxPerJ,GHLon, GHLat, Label)

    IMPLICIT NONE

    ! This Procedure Generates a Label for the Grid History Points
    ! with the Following Layout: DloMloXDlatMlaY_N, Where:

    ! Dlo: Int  Part of Model Longitude in Degree (3 Char)
    ! Mlo: Frac Part of Model Longitude in Minutes (2 Char)
    !   X: Label for Model Longitude Hemisphere
    !      W for Weast and E for East (2 Char)
    ! Dla: Int  Part of Model Latitude in Degree (2 Char)
    ! Mla: Frac Part of Model Latitude in Minutes (2 Char)
    !   Y: Label for Model Latitude Hemisphere
    !      S for South and N for North (1 Char)
    !   N: Record Number at the Input List (6 Char)

    INTEGER, INTENT (IN) :: Nmax ! Number of Grid History Points
    INTEGER, INTENT (IN) :: iMax
    INTEGER, INTENT (IN) :: jMax
    INTEGER, INTENT (IN) :: iMaxPerJ(:)
    REAL (KIND=r8), INTENT (INOUT) :: GHLon(Nmax) ! Model Longitudes of Grid History Points

    REAL (KIND=r8), INTENT (INOUT) :: GHLat(Nmax) ! Model Latitudes  of Grid History Points

    CHARACTER (LEN=*), INTENT (OUT) :: Label(Nmax) ! Output Label Described Above

    INTEGER :: n, &       ! Index for Grid History Points
         LonDeg, &  ! Degree  Part of Model Longitude
         LonSig, &  ! Signal  Part of Model Longitude
         LonMin, &  ! Minutes Part of Model Longitude
         LatDeg, &  ! Degree  Part of Model Latitude
         LatSig, &  ! Signal  Part of Model Latitude
         LatMin     ! Minutes Part of Model Latitude

    CHARACTER (LEN=1) :: LonHem, & ! Hemisphere of Longitude: W or E
         LatHem    ! Hemisphere of Latitude:  S or N

    INTEGER :: i, &     ! Model Zonal      Index
         j        ! Model Meridional Index

    REAL (KIND=r8) :: DegConv, & ! Convertion Factor from Radian to Degree
         RadConv, & ! Convertion Factor from Degree to Radian
         DLon(iMax), & ! Delta Lon to Obtain Nearest Model Longitude
         DLat(jMax)    ! Delta Lat to Obtain Nearest Model Latitude


    INTEGER           :: ILocMin(1)
    INTEGER           :: JLocMin(1)
    REAL(KIND=r8)     :: pi

    pi = 4.0_r8*ATAN(1.0_r8)
    DegConv = 45.0_r8/ATAN(1.0_r8)
    RadConv = ATAN(1.0_r8)/45.0_r8
    GHLon = GHLon  * RadConv
    GHLat=(GHLat+90.0_r8)*RadConv
    DO n=1,Nmax
       DLat=1000.0_r8
       DO j=1,jMax
          dlat(j) = ABS(lati(j) - GHLat(n))
       END DO
       jlocmin=MINLOC (dlat)

       dlon=1000.0_r8
       DO i=1,iMaxPerJ(jlocmin(1))
          dlon(i) = ABS( ((pi*((i-1)*360.0_r8/REAL(iMaxPerJ(jlocmin(1)),r8)))/180.0_r8) &
               &         -  GHLon(n) )
       END DO
       ILocMin=MINLOC (dlon(1:iMaxPerJ(jlocmin(1))))
       GHLat(n)=lati(jlocmin(1))*DegConv
       GHLon(n)= ((pi*((ILocMin(1)-1)*360.0_r8/REAL(iMaxPerJ(jlocmin(1)),r8)))/180.0_r8)*DegConv
    END DO
    WHERE (GHLon > 180.0_r8)
       GHLon=GHLon-360.0_r8
    END WHERE
    GHLat=GHLat-90_r8

    DO n=1,Nmax
       ! Longitudes
       LonDeg=INT(ABS(GHLon(n)))
       LonSig=INT(SIGN(1.0_r8,GHLon(n)))
       LonMin=NINT(60.0_r8*(ABS(GHLon(n))-REAL(LonDeg,r8)))
       IF (LonMin > 99) LonMin=99
       IF (LonSig < 0) THEN
          LonHem='W'
       ELSE
          LonHem='E'
       ENDIF
       ! Latitudes
       LatDeg=INT(ABS(GHLat(n)))
       LatSig=INT(SIGN(1.1_r8,GHLat(n)))
       LatMin=NINT(60.0_r8*(ABS(GHLat(n))-REAL(LatDeg,r8)))
       IF (LatMin > 99) LatMin=99
       IF (LatSig < 0) THEN
          LatHem='S'
       ELSE
          LatHem='N'
       ENDIF
       ! Label
       WRITE (Label(n), FMT='(I3.3,I2.2,A1,2I2.2,2A1,I6.6)') &
            LonDeg, LonMin, LonHem, LatDeg, LatMin, LatHem, '_', n
    END DO
    WHERE (GHLon < 0.0_r8)
       GHLon=GHLon+360.0_r8
    END WHERE

  END SUBROUTINE GridHistLabel

  !----------------------------------------------------------------------------
  SUBROUTINE NULLIFY_grid(a,nlyr )
    TYPE(MapGrHist) :: a
    INTEGER, INTENT(in):: nlyr
    INTEGER :: i
    IF ( ASSOCIATED( a%AlcPt) )  NULLIFY ( a%AlcPt )
    ALLOCATE( a%AlcPt(nlyr) )
    DO i=1, nlyr
       a%AlcPt(i)=0
    END DO
  END SUBROUTINE NULLIFY_grid
  !----------------------------------------------------------------------------
  SUBROUTINE NULLIFY_dignos(a,nlyr )
    TYPE(GridHistPoint) :: a
    INTEGER, INTENT(in):: nlyr
    IF ( ASSOCIATED( a%dignos) )  NULLIFY ( a%dignos )
    ALLOCATE( a%dignos(nlyr) )
  END SUBROUTINE NULLIFY_dignos
  !----------------------------------------------------------------------------
  SUBROUTINE NULLIFY_InReFd(a,nlyr )
    TYPE(GridHistPoint) :: a
    INTEGER, INTENT(in):: nlyr
    IF ( ASSOCIATED( a%InReFd) )  NULLIFY ( a%InReFd )
    ALLOCATE( a%InReFd(nlyr) )
  END SUBROUTINE NULLIFY_InReFd
  !----------------------------------------------------------------------------
  SUBROUTINE WriteGridHistoryTopo (fgzs,TopoGridH,nfghtop)
    REAL(KIND=r8)            , INTENT(IN   ) :: fgzs(:,:)
    CHARACTER(LEN=*), INTENT(IN   ) :: TopoGridH
    INTEGER         , INTENT(IN   ) :: nfghtop
    REAL(KIND=r8)    :: grhtop(ngpts)
    REAL(KIND=r8)    :: grhloc(ngptslocal)
    INTEGER :: n,nloc
    INTEGER :: ierr
    IF (.NOT. grhflg) THEN
       RETURN
    END IF
    IF (myid.eq.0)THEN
       OPEN(UNIT=nfghtop,FILE=TRIM(TopoGridH),FORM='unformatted',ACCESS='sequential',&
            ACTION='write', STATUS='replace',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(TopoGridH), ierr
          STOP "**(ERROR)**"
       END IF
    END IF


    DO nloc=1,ngptslocal
       n = map(nloc)
       grhloc(nloc)=fgzs(GPt(n)%ibLoc,GPt(n)%jbLoc)/grav
    END DO

    IF (maxnodes.gt.1) THEN
       CALL Collect_Grid_His(grhloc, grhtop, ngpts, ngptslocal, 0, 1, &
                             ngptsperproc, mapglobal)
       IF (myid.eq.0) CALL WrTopoGrdHist  (nfghtop,grhtop)
     ELSE
       grhtop=grhloc
       CALL WrTopoGrdHist  (nfghtop,grhtop)
    ENDIF

  END SUBROUTINE WriteGridHistoryTopo



  SUBROUTINE Store2D (field, fId, jb, cf)
    REAL(KIND=r8),    INTENT(IN) :: field(:,:)
    INTEGER, INTENT(IN) :: fId
    INTEGER, INTENT(IN) :: jb
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: cf

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: kfirst
    INTEGER :: k
    INTEGER :: i
    INTEGER :: n
    CHARACTER(LEN=*), PARAMETER :: h = "**(StoreGridHistory)**"

    IF (.NOT. IsGridHistoryOn()) THEN
       RETURN
    END IF
    dim1 = SIZE(field,1)
    dim2 = SIZE(field,2)
    IF (fId < 1 .OR. fId > ngfld) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' fId out of range =', i10)") h, fId
       WRITE(UNIT=nferr, FMT="(a, ' fId out of range =', i10)") h, fId
       STOP h
    ELSE IF (jb < 1 .OR. jb > jbMax) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' jb out of range =', i10)") h, jb
       WRITE(UNIT=nferr, FMT="(a, ' jb out of range =', i10)") h, jb
       STOP h
    ELSE IF (dim1 /= ibMaxPerJB(jb)) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       WRITE(UNIT=nferr, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       STOP h
    ELSE IF (dim2 /= kMax) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field second dimension out of range =', i10)") h, dim2
       WRITE(UNIT=nferr, FMT="(a, ' field second dimension out of range =', i10)") h, dim2
       STOP h
    END IF

    kfirst = GHAF(fId)%LocGrFld
    IF (PRESENT(cf)) THEN
       DO k = 1, kMax
          DO n = iniperjb(jb)+1, iniperjb(jb)+ngptsperjb(jb)
             i = GPt( n )%ibloc
             GPt( n )%dignos(k+kfirst-1) = field(i, k) * cf
          END DO
       END DO
    ELSE
       DO k = 1, kMax
          DO n = iniperjb(jb)+1, iniperjb(jb)+ngptsperjb(jb)
             i = GPt( n )%ibloc
             GPt( n )%dignos(k+kfirst-1) = field(i, k)
          END DO
       END DO
    END IF
  END SUBROUTINE Store2D






  SUBROUTINE Store1D (field, fId, jb, cf)
    REAL(KIND=r8),    INTENT(IN) :: field(:)
    INTEGER, INTENT(IN) :: fId
    INTEGER, INTENT(IN) :: jb
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: cf

    INTEGER :: dim1
    INTEGER :: kfirst
    INTEGER :: i
    INTEGER :: n
    CHARACTER(LEN=*), PARAMETER :: h = "**(StoreGridHistory)**"

    IF (.NOT. IsGridHistoryOn()) THEN
       RETURN
    END IF
    dim1 = SIZE(field,1)
    IF (fId < 1 .OR. fId > ngfld) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' fId out of range =', i10)") h, fId
       WRITE(UNIT=nferr, FMT="(a, ' fId out of range =', i10)") h, fId
       STOP h
    ELSE IF (jb < 1 .OR. jb > jbMax) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' jb out of range =', i10)") h, jb
       WRITE(UNIT=nferr, FMT="(a, ' jb out of range =', i10)") h, jb
       STOP h
    ELSE IF (dim1 /= ibMaxPerJB(jb)) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       WRITE(UNIT=nferr, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       STOP h
    END IF

    kfirst = GHAF(fId)%LocGrFld
    IF (PRESENT(cf)) THEN
       DO n = iniperjb(jb)+1, iniperjb(jb)+ngptsperjb(jb)
          i = GPt( n )%ibloc
          GPt( n )%dignos(kfirst) = field(i) * cf
       END DO
    ELSE
       DO n = iniperjb(jb)+1, iniperjb(jb)+ngptsperjb(jb)
          i = GPt( n )%ibloc
          GPt( n )%dignos(kfirst) = field(i)
       END DO
    END IF
  END SUBROUTINE Store1D






  SUBROUTINE WriteGridHistory (nfghou, ifday, tod, idate)
    INTEGER, INTENT(IN) :: nfghou
    INTEGER, INTENT(IN) :: ifday
    REAL(KIND=r8),    INTENT(IN) :: tod
    INTEGER, INTENT(IN) :: idate(4)

    INTEGER :: j
    INTEGER :: m
    INTEGER :: n,nloc
    INTEGER :: iqstmp(6)
    INTEGER :: isg
    REAL(KIND=r8)    :: sg
    REAL(KIND=r8)    :: stmp(6)
    REAL(KIND=r8)    :: qwork(ngptslocal,nghsl)
    REAL(KIND=r8)    :: work(ngpts,nghsl)

    IF (.NOT. IsGridHistoryOn()) THEN
       RETURN
    END IF
    CALL tmstmp2 (idate, ifday, tod, iqstmp(3), iqstmp(4), iqstmp(5), iqstmp(6))
    sg = MOD(tod+0.03125_r8,3600.0_r8)-0.03125_r8
    isg=INT(sg)
    iqstmp(2)=isg/60
    iqstmp(1)=MOD(isg,60)
    DO j=1,6
       stmp(j)=iqstmp(j)
    END DO
    DO j = 1, nghsl
       DO nloc=1,ngptslocal!ngpts
          n = map(nloc)
          qwork(nloc,j)=GPt( n )%dignos(j)
       END DO
    END DO
    IF(.NOT.allghf.and.ngptslocal.gt.0)THEN
       m = 1
       DO n = 1, ngfld
          IF (GHAF(n)%IsReqFld) THEN
             CALL cnvray(qwork(1,m) , GHAF(n)%LyrsFlds*ngpts,GHAF(n)%UnGraFld,&
                  &           GHAF(n)%UnGraFld  )
             m=m+GHAF(n)%LyrsFlds
          END IF
       END DO
    ENDIF

    IF (maxnodes.gt.1) THEN
       CALL Collect_Grid_His(qwork, work, ngpts, ngptslocal, 0, nghsl, &
                             ngptsperproc, mapglobal)
       IF (myid.eq.0) CALL WriteGrdHist(nfghou,stmp,work)
    ELSE
       CALL WriteGrdHist(nfghou,stmp,qwork)
    ENDIF
    DO n=1,ngptslocal
       GPt(n)%dignos(:)=0.0_r8
    END DO
  END SUBROUTINE WriteGridHistory






  SUBROUTINE StoreMaskedGridHistory (field, imask, fId, jb, cf)
    REAL(KIND=r8),    INTENT(IN) :: field(:)
    INTEGER(KIND=i8), INTENT(IN) :: imask(:)
    INTEGER, INTENT(IN) :: fId
    INTEGER, INTENT(IN) :: jb
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: cf

    INTEGER :: i
    INTEGER :: ncount
    REAL(KIND=r8) :: bfr(ibMax)
    INTEGER :: dim1
    dim1 = SIZE(imask,1)
    IF (.NOT. IsGridHistoryOn()) THEN
       RETURN
    END IF
    ncount=0
    DO i = 1, dim1
       IF (imask(i) >= 1_i8) THEN
          ncount=ncount+1
          bfr(i) = field(ncount)
       ELSE
          bfr(i) = 0.0_r8
       END IF
    END DO
    IF (PRESENT(cf)) THEN
       CALL StoreGridHistory (bfr(1:dim1), fId, jb, cf)
    ELSE
       CALL StoreGridHistory (bfr(1:dim1), fId, jb)
    END IF
  END SUBROUTINE StoreMaskedGridHistory




  SUBROUTINE TurnOnGridHistory()
    grhOn = grhflg
  END SUBROUTINE TurnOnGridHistory




  SUBROUTINE TurnOffGridHistory()
    grhOn = .FALSE.
  END SUBROUTINE TurnOffGridHistory



  FUNCTION IsGridHistoryOn()
    LOGICAL :: IsGridHistoryOn
    IsGridHistoryOn = grhflg .AND. grhOn
  END FUNCTION IsGridHistoryOn







  SUBROUTINE Store2DV (field, fId, jb, cf)
    REAL(KIND=r8),    INTENT(IN) :: field(:,:)
    INTEGER, INTENT(IN) :: fId
    INTEGER, INTENT(IN) :: jb
    REAL(KIND=r8),    INTENT(IN) :: cf(:)

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dimcf
    INTEGER :: kfirst
    INTEGER :: k
    INTEGER :: i
    INTEGER :: n
    CHARACTER(LEN=*), PARAMETER :: h = "**(StoreGridHistory)**"

    IF (.NOT. IsGridHistoryOn()) THEN
       RETURN
    END IF
    dim1 = SIZE(field,1)
    dim2 = SIZE(field,2)
    dimcf = SIZE(cf,1)
    IF (fId < 1 .OR. fId > ngfld) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' fId out of range =', i10)") h, fId
       WRITE(UNIT=nferr, FMT="(a, ' fId out of range =', i10)") h, fId
       STOP h
    ELSE IF (jb < 1 .OR. jb > jbMax) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' jb out of range =', i10)") h, jb
       WRITE(UNIT=nferr, FMT="(a, ' jb out of range =', i10)") h, jb
       STOP h
    ELSE IF (dim1 /= ibMaxPerJB(jb)) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       WRITE(UNIT=nferr, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       STOP h
    ELSE IF (dim2 /= kMax) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field second dimension out of range =', i10)") h, dim2
       WRITE(UNIT=nferr, FMT="(a, ' field second dimension out of range =', i10)") h, dim2
       STOP h
    ELSE IF (dim1 /= dimcf) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field first dimension and cf dimension do not match',2i10)") h, dim1, dimcf
       WRITE(UNIT=nferr, FMT="(a, ' field first dimension and cf dimension do not match',2i10)") h, dim1, dimcf
       STOP h
    END IF

    kfirst = GHAF(fId)%LocGrFld
    DO k = 1, kMax
       DO n = iniperjb(jb)+1, iniperjb(jb)+ngptsperjb(jb)
          i = GPt( n )%ibloc
          GPt( n )%dignos(k+kfirst-1) = field(i, k) * cf(i)
       END DO
    END DO
  END SUBROUTINE Store2DV






  SUBROUTINE Store1DV (field, fId, jb, cf)
    REAL(KIND=r8),    INTENT(IN) :: field(:)
    INTEGER, INTENT(IN) :: fId
    INTEGER, INTENT(IN) :: jb
    REAL(KIND=r8),    INTENT(IN) :: cf(:)

    INTEGER :: dim1
    INTEGER :: dimcf
    INTEGER :: kfirst
    INTEGER :: i
    INTEGER :: n
    CHARACTER(LEN=*), PARAMETER :: h = "**(StoreGridHistory)**"
    IF (.NOT. IsGridHistoryOn()) THEN
       RETURN
    END IF
    dim1 = SIZE(field,1)
    dimcf = SIZE(cf,1)
    IF (fId < 1 .OR. fId > ngfld) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' fId out of range =', i10)") h, fId
       WRITE(UNIT=nferr, FMT="(a, ' fId out of range =', i10)") h, fId
       STOP h
    ELSE IF (jb < 1 .OR. jb > jbMax) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' jb out of range =', i10)") h, jb
       WRITE(UNIT=nferr, FMT="(a, ' jb out of range =', i10)") h, jb
       STOP h
    ELSE IF (dim1 /= ibMaxPerJB(jb)) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       WRITE(UNIT=nferr, FMT="(a, ' field first dimension out of range =', 3i10)") h, dim1,ibMaxPerJB(jb),ibMax
       STOP h
    ELSE IF (dim1 /= dimcf) THEN
       WRITE(UNIT=nfprt, FMT="(a, ' field first dimension and cf dimension do not match',2i10)") h, dim1, dimcf
       WRITE(UNIT=nferr, FMT="(a, ' field first dimension and cf dimension do not match',2i10)") h, dim1, dimcf
       STOP h
    END IF

    kfirst = GHAF(fId)%LocGrFld
    DO n = iniperjb(jb)+1, iniperjb(jb)+ngptsperjb(jb)
       i = GPt( n )%ibloc
       GPt( n )%dignos(kfirst) = field(i) * cf(i)
    END DO
  END SUBROUTINE Store1DV

END MODULE GridHistory
