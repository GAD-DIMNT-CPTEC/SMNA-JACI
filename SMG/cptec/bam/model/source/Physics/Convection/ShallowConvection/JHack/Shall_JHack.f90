MODULE Shall_JHack
  USE uwshcu, Only : init_uwshcu,fqsatd,compute_uwshcu_inv
  IMPLICIT NONE
SAVE

  SAVE
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers
  !----------------------------------------------------------------------------
  ! physical constants (all data public)
  !----------------------------------------------------------------------------
  REAL(r8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number ~ molecules/kmole
  REAL(r8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant ~ J/K/molecule
  REAL(r8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(r8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air ~ kg/kmole
  REAL(r8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(r8),PARAMETER :: SHR_CONST_G      = 9.80616_r8      ! acceleration of gravity ~ m/s^2
  REAL(r8),PARAMETER :: SHR_CONST_TKFRZ  = 273.16_r8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(r8),PARAMETER :: SHR_CONST_LATVAP = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_LATICE = 3.337e5_r8      ! latent heat of fusion ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_RHOFW  = 1.000e3_r8      ! density of fresh water ~ kg/m^3

  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg

  INTEGER ,PARAMETER :: ixcldliq=2
  INTEGER ,PARAMETER :: ixcldice=3

  ! Constants for air

  REAL(r8), PUBLIC, PARAMETER :: rair = shr_const_rdair    ! Gas constant for dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: zvir = SHR_CONST_RWV/rair - 1          ! rh2o/rair - 1

  ! Constants for Earth
  REAL(r8), PUBLIC, PARAMETER :: gravit = shr_const_g      ! gravitational acceleration
  ! Constants for water

  REAL(r8), PUBLIC, PARAMETER :: tmelt = shr_const_tkfrz   ! Freezing point of water
  REAL(r8), PUBLIC, PARAMETER :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights 
  REAL(r8), PUBLIC, PARAMETER :: latvap = shr_const_latvap ! Latent heat of vaporization
  REAL(r8), PUBLIC, PARAMETER :: latice = shr_const_latice ! Latent heat of fusion
  REAL(r8), PRIVATE :: rhoh2o = shr_const_rhofw  ! Density of liquid water (STP)
  REAL(r8), PUBLIC, PARAMETER :: rh2o    =SHR_CONST_RWV   !! Gas constant for water vapor
  LOGICAL :: PERGRO=.FALSE.
  LOGICAL :: cnst_need_pdeldry=.FALSE.

  REAL(r8) :: ke                     ! Tunable evaporation efficiency

  !------------wv_saturation-------------
  !
  ! Data
  !
  INTEGER, PARAMETER :: plenest = 250! length of saturation vapor pressure table

  REAL(r8)         :: estbl(plenest)      ! table values of saturation vapor pressure
  REAL(r8),PARAMETER :: tmn  = 173.16_r8          ! Minimum temperature entry in table
  REAL(r8),PARAMETER :: tmx  = 375.16_r8          ! Maximum temperature entry in table
  REAL(r8),PARAMETER :: trice  =  20.00_r8         ! Trans range from es over h2o to es over ice
  REAL(r8),PARAMETER :: tmin=tmn       ! min temperature (K) for table
  REAL(r8),PARAMETER :: tmax= tmx      ! max temperature (K) for table
  LOGICAL  ,PARAMETER :: icephs=.TRUE.  ! false => saturation vapor press over water only
  INTEGER,PARAMETER   ::  iterp =2             ! #iterations for precipitation calculation
  REAL(r8) :: pcf(6)     ! polynomial coeffs -> es transition water to ice
  REAL(r8) ,PARAMETER  :: ttrice=trice
  REAL(r8), PUBLIC, PARAMETER :: tmax_fice = tmelt - 10.0_r8       ! max temperature for cloud ice formation
  REAL(r8), PUBLIC, PARAMETER :: tmin_fice = tmax_fice - 30.0_r8   ! min temperature for cloud ice formation
  REAL(r8), PUBLIC, PARAMETER :: tmax_fsnow = tmelt            ! max temperature for transition to convective snow
  REAL(r8), PUBLIC, PARAMETER :: tmin_fsnow = tmelt-5.0_r8         ! min temperature for transition to convective snow

  REAL(r8) :: epsqs=epsilo

  !------------wv_saturation end -------------
  REAL(r8), PRIVATE :: hlatv  = latvap
  REAL(r8), PRIVATE :: hlatf  = latice
  REAL(r8), PRIVATE :: rgasv  = SHR_CONST_RWV    ! Gas constant for water vapor
  !------------  module physics_types-------------
  !-------------------------------------------------------------------------------
  TYPE physics_state
     INTEGER , POINTER        :: ncol(:) ! number of active columns
     REAL(r8), POINTER :: lat           (:,:) !(pcols)     ! latitude (radians)
     REAL(r8), POINTER :: lon           (:,:) !(pcols)     ! longitude (radians)
     REAL(r8), POINTER :: ps           (:,:) !(pcols)     ! surface pressure
     !     REAL(r8), POINTER :: psdry           (:,:) !(pcols)     ! dry surface pressure
     REAL(r8), POINTER :: phis           (:,:) !(pcols)     ! surface geopotential
     REAL(r8), POINTER :: t           (:,:,:) !(pcols,pver)! temperature (K)
     REAL(r8), POINTER :: u           (:,:,:) !(pcols,pver)! zonal wind (m/s)
     REAL(r8), POINTER :: v           (:,:,:) !(pcols,pver)! meridional wind (m/s)
     REAL(r8), POINTER :: s           (:,:,:) !(pcols,pver)! dry static energy
     REAL(r8), POINTER :: omega           (:,:,:) !(pcols,pver)! vertical pressure velocity (Pa/s) 
     REAL(r8), POINTER :: pmid           (:,:,:) !(pcols,pver)! midpoint pressure (Pa) 
     !     REAL(r8), POINTER :: pmiddry  (:,:,:) !(pcols,pver)! midpoint pressure dry (Pa) 
     REAL(r8), POINTER :: pdel           (:,:,:) !(pcols,pver)! layer thickness (Pa)
     !     REAL(r8), POINTER :: pdeldry  (:,:,:) !(pcols,pver)! layer thickness dry (Pa)
     REAL(r8), POINTER :: rpdel           (:,:,:) !(pcols,pver)! reciprocal of layer thickness (Pa)
     !     REAL(r8), POINTER :: rpdeldry (:,:,:) !(pcols,pver)! recipricol layer thickness dry (Pa)
     REAL(r8), POINTER :: lnpmid   (:,:,:) !(pcols,pver)! ln(pmid)
     !     REAL(r8), POINTER :: lnpmiddry(:,:,:) !(pcols,pver)! log midpoint pressure dry (Pa) 
     !     REAL(r8), POINTER :: exner           (:,:,:) !(pcols,pver)! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
     REAL(r8), POINTER :: zm           (:,:,:) !(pcols,pver)! geopotential height above surface at midpoints (m)
     REAL(r8), POINTER :: q        (:,:,:,:) !(pcols,pver,ppcnst)! constituent mixing ratio (kg/kg moist or dry air depending on type)
     REAL(r8), POINTER :: pint           (:,:,:) !(pcols,pver+1)! interface pressure (Pa)
     !     REAL(r8), POINTER :: pintdry  (:,:,:) !(pcols,pver+1)! interface pressure dry (Pa) 
     REAL(r8), POINTER :: lnpint   (:,:,:) !(pcols,pver+1)! ln(pint)
     !     REAL(r8), POINTER :: lnpintdry(:,:,:) !(pcols,pver+1)! log interface pressure dry (Pa) 
     REAL(r8), POINTER :: zi           (:,:,:) !(pcols,pver+1)! geopotential height above surface at interfaces (m)
     REAL(r8), POINTER :: te_ini   (:,:) !(pcols)  ! vertically integrated total (kinetic + static) energy of initial state
     REAL(r8), POINTER :: te_cur   (:,:) !(pcols)  ! vertically integrated total (kinetic + static) energy of current state
     REAL(r8), POINTER :: tw_ini   (:,:) !(pcols)  ! vertically integrated total water of initial state
     REAL(r8), POINTER :: tw_cur   (:,:) !(pcols)  ! vertically integrated total water of new state
     INTEGER , POINTER  :: COUNT (:)              ! count of values with significant energy or water imbalances     
  END TYPE physics_state
  TYPE(physics_state):: state
  TYPE(physics_state):: state1
  !-------------------------------------------------------------------------------
  TYPE physics_tend
     REAL(r8), POINTER :: dtdt   (:,:,:) !(pcols,pver) 
     REAL(r8), POINTER :: dudt   (:,:,:) !(pcols,pver) 
     REAL(r8), POINTER :: dvdt   (:,:,:) !(pcols,pver) 
     REAL(r8), POINTER :: flx_net(:,:) !(pcols          ) 
     REAL(r8), POINTER :: te_tnd (:,:) !(pcols)      ! cumulative boundary flux of total energy
     REAL(r8), POINTER :: tw_tnd (:,:) !(pcols)      ! cumulative boundary flux of total water
  END TYPE physics_tend
  TYPE(physics_tend )  :: tend        ! Physics tendencies (empty, needed for physics_update call)

  !-------------------------------------------------------------------------------
  ! This is for tendencies returned from individual parameterizations
  TYPE physics_ptend
     CHARACTER(LEN=24), POINTER :: name(:)   ! name of parameterization which produced tendencies.
     LOGICAL   , POINTER:: ls(:)                ! true if dsdt is returned
     LOGICAL   , POINTER:: lu(:)                ! true if dudt is returned
     LOGICAL   , POINTER:: lv(:)                ! true if dvdt is returned
     LOGICAL  , POINTER :: lq (:,:)  !(ppcnst)       ! true if dqdt() is returned
     INTEGER  , POINTER :: top_level(:)        ! top level index for which nonzero tendencies have been set
     INTEGER  , POINTER :: bot_level(:)        ! bottom level index for which nonzero tendencies have been set
     REAL(r8), POINTER :: s(:,:,:) !(pcols,pver)! heating rate (J/kg/s)
     REAL(r8), POINTER :: u(:,:,:) !(pcols,pver)! u momentum tendency (m/s/s)
     REAL(r8), POINTER :: v(:,:,:) !(pcols,pver)! v momentum tendency (m/s/s)
     REAL(r8), POINTER :: q(:,:,:,:)!(pcols,pver,ppcnst)                 ! consituent tendencies (kg/kg/s)

     ! boundary fluxes

     REAL(r8), POINTER :: hflux_srf(:,:) !(pcols)! net heat flux at surface (W/m2)
     REAL(r8), POINTER :: hflux_top(:,:) !(pcols)! net heat flux at top of model (W/m2)
     REAL(r8), POINTER :: taux_srf (:,:) !(pcols)! net zonal stress at surface (Pa)
     REAL(r8), POINTER :: taux_top (:,:) !(pcols)! net zonal stress at top of model (Pa)
     REAL(r8), POINTER :: tauy_srf (:,:) !(pcols)! net meridional stress at surface (Pa)
     REAL(r8), POINTER :: tauy_top (:,:) !(pcols)! net meridional stress at top of model (Pa)
     REAL(r8), POINTER :: cflx_srf (:,:,:) !(pcols,ppcnst)! constituent flux at surface (kg/m2/s)
     REAL(r8), POINTER :: cflx_top (:,:,:) !(pcols,ppcnst)! constituent flux top of model (kg/m2/s)
  END TYPE physics_ptend
  TYPE(physics_ptend)  :: ptend_loc   ! package tendencies
  TYPE(physics_ptend)        :: ptend_all   ! package tendencies
  !REAL(r8),POINTER ::  state_concld (:,:,:) !(1:pcols,1:pver,latco)
  !REAL(r8),POINTER ::  state_cld    (:,:,:) !(1:pcols,1:pver,latco)
 ! REAL(r8),POINTER ::  state_icwmr  (:,:,:) !(1:pcols,1:pver,latco)
 ! REAL(r8),POINTER ::  state_rprddp (:,:,:) !(1:pcols,1:pver,latco)
 ! REAL(r8),POINTER ::  state_rprdsh (:,:,:) !(1:pcols,1:pver,latco)
  !REAL(r8),POINTER ::  state_RPRDTOT(:,:,:) !(1:pcols,1:pver,latco)
!  REAL(r8),POINTER ::  state_cnt    (:,:) !(1:pcols,latco)
!  REAL(r8),POINTER ::  state_cnb    (:,:) !(1:pcols,latco)
 ! REAL(r8),POINTER ::  state_shfrc  (:,:,:) 
  !REAL(r8),POINTER ::  state_evapcsh(:,:,:) 
!  REAL(r8),POINTER ::  state_cush     (:,:) !(1:pcols,latco)=> pbuf(pbuf_get_fld_idx('cush'))%fld_ptr(1,1:pcols,1,lchnk,itim)

  !---------------------------------------------------------------------------------
  ! Purpose:
  !
  ! CAM interface to the Hack shallow convection scheme
  !
  ! Author: D.B. Coleman
  !
  !---------------------------------------------------------------------------------
  PRIVATE


  !
  ! Private data used for Hack shallow convection
  !
  REAL(r8) :: hlat        ! latent heat of vaporization
  REAL(r8) :: c0          ! rain water autoconversion coefficient
  REAL(r8) :: betamn      ! minimum overshoot parameter
  REAL(r8) :: rhlat       ! reciprocal of hlat
  REAL(r8) :: rcp         ! reciprocal of cp
  REAL(r8) :: cmftau      ! characteristic adjustment time scale
  REAL(r8) :: dzmin       ! minimum convective depth for precipitation
  REAL(r8) :: tiny        ! arbitrary small num used in transport estimates
  REAL(r8) :: eps         ! convergence criteria (machine dependent)
  REAL(r8) :: tpmax       ! maximum acceptable t perturbation (degrees C)
  REAL(r8) :: shpmax      ! maximum acceptable q perturbation (g/g)           

  INTEGER :: iloc         ! longitude location for diagnostics
  INTEGER :: jloc         ! latitude  location for diagnostics
  INTEGER :: nsloc        ! nstep for which to produce diagnostics
  !
  LOGICAL :: rlxclm       ! logical to relax column versus cloud triplet

  REAL(r8) :: cp          ! specific heat of dry air
  REAL(r8) :: grav        ! gravitational constant       
  REAL(r8) :: rgrav       ! reciprocal of grav
  REAL(r8) :: rgas        ! gas constant for dry air
  INTEGER   :: limcnv          ! top interface level limit for convection

  INTEGER, PARAMETER  :: ppcnst=3
  REAL(r8) :: qmin(3)
  CHARACTER*3, PUBLIC :: cnst_type(1:ppcnst)          ! wet or dry mixing ratio
  CHARACTER*3, PARAMETER :: mixtype(1:ppcnst)=(/'wet','wet', 'wet'/) ! mixing ratio type (dry, wet)
  ! The following namelist variable controls which shallow convection package is used.
  !           'Hack'   = Hack shallow convection (default)
  !           'UW'     = UW shallow convection by Sungsu Park and Christopher S. Bretherton
  !           'off'    = No shallow convection

  character(len=16) :: shallow_scheme      ! Default set in phys_control.F90, use namelist to change

  ! Public methods

  PUBLIC ::&
       convect_shallow_init,               &! initialize donner_shallow module
       convect_shallow_tend                 ! return tendencies

  !=========================================================================================
CONTAINS
  !=========================================================================================

  SUBROUTINE convect_shallow_init(ISCON,kMax,jMax,ibMax,jbMax,ppcnst,&
                                  a_hybr,b_hybr)

    !----------------------------------------
    ! Purpose:  declare output fields, initialize variables needed by convection
    !----------------------------------------

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN   ) :: ISCON
    INTEGER, INTENT(in   ) :: kMax                  ! number of vertical levels
    INTEGER, INTENT(IN   ) :: jMax
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: ppcnst
    REAL(KIND=r8), INTENT(IN   ) :: a_hybr   (kMax+1)
    REAL(KIND=r8), INTENT(IN   ) :: b_hybr   (kMax+1)
    LOGICAL ip           ! Ice phase (true or false)
    REAL(KIND=r8) :: hypi (kMax+1)! reference pressures at interfaces



    INTEGER limcnv          ! top interface level limit for convection
    INTEGER k
    INTEGER ind


    PRINT*,'convect_shallow_init',jbMax,kMax,ibMax
    ALLOCATE(state%ncol     (1:jbMax))                  ;state%ncol     =0     ! number of active columns
    ALLOCATE(state%lat      (1:ibMax,1:jbMax))            ;state%lat      =0.0_r8!(pcols)     ! latitude (radians)
    ALLOCATE(state%lon      (1:ibMax,1:jbMax))            ;state%lon      =0.0_r8!(pcols)     ! longitude (radians)
    ALLOCATE(state%ps             (1:ibMax,1:jbMax))            ;state%ps       =0.0_r8!(pcols)     ! surface pressure
    !    ALLOCATE(state%psdry    (1:ibMax,1:jbMax))            ;state%psdry    =0.0_r8!(pcols)     ! dry surface pressure
    ALLOCATE(state%phis     (1:ibMax,1:jbMax))            ;state%phis     =0.0_r8!(pcols)     ! surface geopotential
    ALLOCATE(state%t             (1:ibMax,1:kMax,1:jbMax))       ;state%t        =0.0_r8!(pcols,pver)! temperature (K)
    ALLOCATE(state%u             (1:ibMax,1:kMax,1:jbMax))       ;state%u        =0.0_r8!(pcols,pver)! zonal wind (m/s)
    ALLOCATE(state%v             (1:ibMax,1:kMax,1:jbMax))       ;state%v        =0.0_r8!(pcols,pver)! meridional wind (m/s)
    ALLOCATE(state%s             (1:ibMax,1:kMax,1:jbMax))       ;state%s        =0.0_r8!(pcols,pver)! dry static energy
    ALLOCATE(state%omega    (1:ibMax,1:kMax,1:jbMax))       ;state%omega    =0.0_r8!(pcols,pver)! vertical pressure velocity (Pa/s) 
    ALLOCATE(state%pmid     (1:ibMax,1:kMax,1:jbMax))       ;state%pmid     =0.0_r8!(pcols,pver)! midpoint pressure (Pa) 
    !    ALLOCATE(state%pmiddry  (1:ibMax,1:kMax,1:jbMax))       ;state%pmiddry  =0.0_r8!(pcols,pver)! midpoint pressure dry (Pa) 
    ALLOCATE(state%pdel     (1:ibMax,1:kMax,1:jbMax))       ;state%pdel     =0.0_r8!(pcols,pver)! layer thickness (Pa)
    !    ALLOCATE(state%pdeldry  (1:ibMax,1:kMax,1:jbMax))       ;state%pdeldry  =0.0_r8!(pcols,pver)! layer thickness dry (Pa)
    ALLOCATE(state%rpdel    (1:ibMax,1:kMax,1:jbMax))       ;state%rpdel    =0.0_r8!(pcols,pver)! reciprocal of layer thickness (Pa)
    !    ALLOCATE(state%rpdeldry (1:ibMax,1:kMax,1:jbMax))       ;state%rpdeldry =0.0_r8!(pcols,pver)! recipricol layer thickness dry (Pa)
    ALLOCATE(state%lnpmid   (1:ibMax,1:kMax,1:jbMax))       ;state%lnpmid   =0.0_r8! (pcols,pver)! ln(pmid)
    !    ALLOCATE(state%lnpmiddry(1:ibMax,1:kMax,1:jbMax))       ;state%lnpmiddry=0.0_r8!(pcols,pver)! log midpoint pressure dry (Pa) 
    !    ALLOCATE(state%exner    (1:ibMax,1:kMax,1:jbMax))       ;state%exner    =0.0_r8!(pcols,pver)! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
    ALLOCATE(state%zm             (1:ibMax,1:kMax,1:jbMax))       ;state%zm       =0.0_r8!(pcols,pver)! geopotential height above surface at midpoints (m)
    ALLOCATE(state%q         (1:ibMax,1:kMax,1:jbMax,1:ppcnst));state%q        =0.0_r8!(pcols,pver,ppcnst)! constituent mixing ratio 
    !(kg/kg moist or dry air depending on type)
    ALLOCATE(state%pint     (1:ibMax,1:kMax+1,1:jbMax))     ;state%pint     =0.0_r8!(pcols,pver+1)! interface pressure (Pa)
    !    ALLOCATE(state%pintdry  (1:ibMax,1:kMax+1,1:jbMax))     ;state%pintdry  =0.0_r8!(pcols,pver+1)! interface pressure dry (Pa) 
    ALLOCATE(state%lnpint   (1:ibMax,1:kMax+1,1:jbMax))     ;state%lnpint   =0.0_r8!(pcols,pver+1)! ln(pint)
    !    ALLOCATE(state%lnpintdry(1:ibMax,1:kMax+1,1:jbMax))     ;state%lnpintdry=0.0_r8!(pcols,pver+1)! log interface pressure dry (Pa) 
    ALLOCATE(state%zi             (1:ibMax,1:kMax+1,1:jbMax))     ;state%zi       =0.0_r8!(pcols,pver+1)! geopotential height above surface at interfaces (m)
    ALLOCATE(state%te_ini   (1:ibMax,1:jbMax))            ;state%te_ini   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of initial state
    ALLOCATE(state%te_cur   (1:ibMax,1:jbMax))            ;state%te_cur   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of current state
    ALLOCATE(state%tw_ini   (1:ibMax,1:jbMax))            ;state%tw_ini   =0.0_r8!(pcols)  ! vertically integrated total water of initial state
    ALLOCATE(state%tw_cur   (1:ibMax,1:jbMax))            ;state%tw_cur   =0.0_r8!(pcols)  ! vertically integrated total water of new state
    ALLOCATE(state%count    (1:jbMax))                  ;state%count    =0! count of values with significant energy or water imbalances          

    ALLOCATE(state1%ncol     (1:jbMax))                  ;state1%ncol     =0     ! number of active columns
    ALLOCATE(state1%lat      (1:ibMax,1:jbMax))            ;state1%lat      =0.0_r8!(pcols)     ! latitude (radians)
    ALLOCATE(state1%lon      (1:ibMax,1:jbMax))            ;state1%lon      =0.0_r8!(pcols)     ! longitude (radians)
    ALLOCATE(state1%ps             (1:ibMax,1:jbMax))            ;state1%ps       =0.0_r8!(pcols)     ! surface pressure
    !    ALLOCATE(state1%psdry    (1:ibMax,1:jbMax))            ;state1%psdry    =0.0_r8!(pcols)     ! dry surface pressure
    ALLOCATE(state1%phis     (1:ibMax,1:jbMax))            ;state1%phis     =0.0_r8!(pcols)     ! surface geopotential
    ALLOCATE(state1%t             (1:ibMax,1:kMax,1:jbMax))       ;state1%t        =0.0_r8!(pcols,pver)! temperature (K)
    ALLOCATE(state1%u             (1:ibMax,1:kMax,1:jbMax))       ;state1%u        =0.0_r8!(pcols,pver)! zonal wind (m/s)
    ALLOCATE(state1%v             (1:ibMax,1:kMax,1:jbMax))       ;state1%v        =0.0_r8!(pcols,pver)! meridional wind (m/s)
    ALLOCATE(state1%s             (1:ibMax,1:kMax,1:jbMax))       ;state1%s        =0.0_r8!(pcols,pver)! dry static energy
    ALLOCATE(state1%omega    (1:ibMax,1:kMax,1:jbMax))       ;state1%omega    =0.0_r8!(pcols,pver)! vertical pressure velocity (Pa/s) 
    ALLOCATE(state1%pmid     (1:ibMax,1:kMax,1:jbMax))       ;state1%pmid     =0.0_r8!(pcols,pver)! midpoint pressure (Pa) 
    !    ALLOCATE(state1%pmiddry  (1:ibMax,1:kMax,1:jbMax))       ;state1%pmiddry  =0.0_r8!(pcols,pver)! midpoint pressure dry (Pa) 
    ALLOCATE(state1%pdel     (1:ibMax,1:kMax,1:jbMax))       ;state1%pdel     =0.0_r8!(pcols,pver)! layer thickness (Pa)
    !    ALLOCATE(state1%pdeldry  (1:ibMax,1:kMax,1:jbMax))       ;state1%pdeldry  =0.0_r8!(pcols,pver)! layer thickness dry (Pa)
    ALLOCATE(state1%rpdel    (1:ibMax,1:kMax,1:jbMax))       ;state1%rpdel    =0.0_r8!(pcols,pver)! reciprocal of layer thickness (Pa)
    !    ALLOCATE(state1%rpdeldry (1:ibMax,1:kMax,1:jbMax))       ;state1%rpdeldry =0.0_r8!(pcols,pver)! recipricol layer thickness dry (Pa)
    ALLOCATE(state1%lnpmid   (1:ibMax,1:kMax,1:jbMax))       ;state1%lnpmid   =0.0_r8! (pcols,pver)! ln(pmid)
    !    ALLOCATE(state1%lnpmiddry(1:ibMax,1:kMax,1:jbMax))       ;state1%lnpmiddry=0.0_r8!(pcols,pver)! log midpoint pressure dry (Pa) 
    !    ALLOCATE(state1%exner    (1:ibMax,1:kMax,1:jbMax))       ;state1%exner    =0.0_r8!(pcols,pver)! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
    ALLOCATE(state1%zm             (1:ibMax,1:kMax,1:jbMax))       ;state1%zm       =0.0_r8!(pcols,pver)! geopotential height above surface at midpoints (m)
    ALLOCATE(state1%q        (1:ibMax,1:kMax,1:jbMax,1:ppcnst));state1%q        =0.0_r8!(pcols,pver,ppcnst)! constituent mixing ratio 
    !(kg/kg moist or dry air depending on type)
    ALLOCATE(state1%pint     (1:ibMax,1:kMax+1,1:jbMax))     ;state1%pint     =0.0_r8!(pcols,pver+1)! interface pressure (Pa)
    !    ALLOCATE(state1%pintdry  (1:ibMax,1:kMax+1,1:jbMax))     ;state1%pintdry  =0.0_r8!(pcols,pver+1)! interface pressure dry (Pa) 
    ALLOCATE(state1%lnpint   (1:ibMax,1:kMax+1,1:jbMax))     ;state1%lnpint   =0.0_r8!(pcols,pver+1)! ln(pint)
    !    ALLOCATE(state1%lnpintdry(1:ibMax,1:kMax+1,1:jbMax))     ;state1%lnpintdry=0.0_r8!(pcols,pver+1)! log interface pressure dry (Pa) 
    ALLOCATE(state1%zi             (1:ibMax,1:kMax+1,1:jbMax))     ;state1%zi       =0.0_r8!(pcols,pver+1)! geopotential height above surface at interfaces (m)
    ALLOCATE(state1%te_ini   (1:ibMax,1:jbMax))            ;state1%te_ini   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of initial state
    ALLOCATE(state1%te_cur   (1:ibMax,1:jbMax))            ;state1%te_cur   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of current state
    ALLOCATE(state1%tw_ini   (1:ibMax,1:jbMax))            ;state1%tw_ini   =0.0_r8!(pcols)  ! vertically integrated total water of initial state
    ALLOCATE(state1%tw_cur   (1:ibMax,1:jbMax))            ;state1%tw_cur   =0.0_r8!(pcols)  ! vertically integrated total water of new state
    ALLOCATE(state1%count    (1:jbMax))                  ;state1%count    =0! count of values with significant energy or water imbalances          
    ALLOCATE(tend%dtdt   (1:ibMax,1:kMax,1:jbMax));tend%dtdt        =0.0_r8!(pcols,pver) 
    ALLOCATE(tend%dudt   (1:ibMax,1:kMax,1:jbMax));tend%dudt        =0.0_r8!(pcols,pver) 
    ALLOCATE(tend%dvdt   (1:ibMax,1:kMax,1:jbMax));tend%dvdt        =0.0_r8!(pcols,pver) 
    ALLOCATE(tend%flx_net(1:ibMax,1:jbMax))     ;tend%flx_net=0.0_r8!(pcols      ) 
    ALLOCATE(tend%te_tnd (1:ibMax,1:jbMax))     ;tend%te_tnd =0.0_r8!(pcols)      ! cumulative boundary flux of total energy
    ALLOCATE(tend%tw_tnd (1:ibMax,1:jbMax))     ;tend%tw_tnd =0.0_r8!(pcols)      ! cumulative boundary flux of total water


    !-------------------------------------------------------------------------------
    ! This is for tendencies returned from individual parameterizations
    ALLOCATE(ptend_loc%name(1:jbMax))                  ;ptend_loc%name=''!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_loc%ls  (1:jbMax))                  ;ptend_loc%ls=.TRUE.!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_loc%lu  (1:jbMax))                  ;ptend_loc%lu=.TRUE.!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_loc%lv  (1:jbMax))                  ;ptend_loc%lv=.TRUE.!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_loc%top_level (1:jbMax))            ;ptend_loc%top_level=-1
    ALLOCATE(ptend_loc%bot_level (1:jbMax))            ;ptend_loc%bot_level=-1
    ALLOCATE(ptend_loc%lq  (1:jbMax,1:ppcnst))           ;ptend_loc%lq=.TRUE.!(ppcnst)         ! true if dqdt() is returned
    ALLOCATE(ptend_loc%s   (1:ibMax,1:kMax,1:jbMax))       ;ptend_loc%s =0.0_r8!(pcols,pver)! heating rate (J/kg/s)
    ALLOCATE(ptend_loc%u   (1:ibMax,1:kMax,1:jbMax))       ;ptend_loc%u =0.0_r8!(pcols,pver)! u momentum tendency (m/s/s)
    ALLOCATE(ptend_loc%v   (1:ibMax,1:kMax,1:jbMax))       ;ptend_loc%v =0.0_r8!(pcols,pver)! v momentum tendency (m/s/s)
    ALLOCATE(ptend_loc%q   (1:ibMax,1:kMax,1:jbMax,1:ppcnst));ptend_loc%q =0.0_r8!(pcols,pver,ppcnst)! consituent tendencies (kg/kg/s)

    ! boundary fluxes

    ALLOCATE(ptend_loc%hflux_srf  (1:ibMax,1:jbMax))       ;ptend_loc%hflux_srf =0.0_r8!(pcols)! net heat flux at surface (W/m2)
    ALLOCATE(ptend_loc%hflux_top  (1:ibMax,1:jbMax))       ;ptend_loc%hflux_top =0.0_r8!(pcols)! net heat flux at top of model (W/m2)
    ALLOCATE(ptend_loc%taux_srf   (1:ibMax,1:jbMax))       ;ptend_loc%taux_srf  =0.0_r8!(pcols)! net zonal stress at surface (Pa)
    ALLOCATE(ptend_loc%taux_top   (1:ibMax,1:jbMax))       ;ptend_loc%taux_top  =0.0_r8!(pcols)! net zonal stress at top of model (Pa)
    ALLOCATE(ptend_loc%tauy_srf   (1:ibMax,1:jbMax))       ;ptend_loc%tauy_srf  =0.0_r8!(pcols)! net meridional stress at surface (Pa)
    ALLOCATE(ptend_loc%tauy_top   (1:ibMax,1:jbMax))       ;ptend_loc%tauy_top  =0.0_r8!(pcols)! net meridional stress at top of model (Pa)
    ALLOCATE(ptend_loc%cflx_srf   (1:ibMax,1:jbMax,1:ppcnst));ptend_loc%cflx_srf  =0.0_r8!(pcols,ppcnst)! constituent flux at surface (kg/m2/s)
    ALLOCATE(ptend_loc%cflx_top   (1:ibMax,1:jbMax,1:ppcnst));ptend_loc%cflx_top  =0.0_r8!(pcols,ppcnst)! constituent flux top of model (kg/m2/s)
    !-------------------------------------------------------------------------------
    ! This is for tendencies returned from individual parameterizations
    ALLOCATE(ptend_all%name(1:jbMax))                  ;ptend_all%name=''!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_all%ls(1:jbMax))                    ;ptend_all%ls=.TRUE.!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_all%lu(1:jbMax))                    ;ptend_all%lu=.TRUE.!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_all%lv(1:jbMax))                    ;ptend_all%lv=.TRUE.!(ppcnst)         ! true if dqdt() is returned    
    ALLOCATE(ptend_all%top_level (1:jbMax))            ;ptend_all%top_level=-1
    ALLOCATE(ptend_all%bot_level (1:jbMax))            ;ptend_all%bot_level=-1
    ALLOCATE(ptend_all%lq  (1:jbMax,1:ppcnst))           ;ptend_all%lq =.TRUE.! (ppcnst)          ! true if dqdt() is returned
    ALLOCATE(ptend_all%s   (1:ibMax,1:kMax,1:jbMax))       ;ptend_all%s  =0.0_r8!(pcols,pver)! heating rate (J/kg/s)
    ALLOCATE(ptend_all%u   (1:ibMax,1:kMax,1:jbMax))       ;ptend_all%u  =0.0_r8!(pcols,pver)! u momentum tendency (m/s/s)
    ALLOCATE(ptend_all%v   (1:ibMax,1:kMax,1:jbMax))       ;ptend_all%v  =0.0_r8!(pcols,pver)! v momentum tendency (m/s/s)
    ALLOCATE(ptend_all%q   (1:ibMax,1:kMax,1:jbMax,1:ppcnst));ptend_all%q  =0.0_r8!(pcols,pver,ppcnst)! consituent tendencies (kg/kg/s)

    ! boundary fluxes

    ALLOCATE(ptend_all%hflux_srf  (1:ibMax,1:jbMax))       ;ptend_all%hflux_srf=0.0_r8 !(pcols)! net heat flux at surface (W/m2)
    ALLOCATE(ptend_all%hflux_top  (1:ibMax,1:jbMax))       ;ptend_all%hflux_top =0.0_r8!(pcols)! net heat flux at top of model (W/m2)
    ALLOCATE(ptend_all%taux_srf   (1:ibMax,1:jbMax))       ;ptend_all%taux_srf =0.0_r8 !(pcols)! net zonal stress at surface (Pa)
    ALLOCATE(ptend_all%taux_top   (1:ibMax,1:jbMax))       ;ptend_all%taux_top =0.0_r8 !(pcols)! net zonal stress at top of model (Pa)
    ALLOCATE(ptend_all%tauy_srf   (1:ibMax,1:jbMax))       ;ptend_all%tauy_srf =0.0_r8 !(pcols)! net meridional stress at surface (Pa)
    ALLOCATE(ptend_all%tauy_top   (1:ibMax,1:jbMax))       ;ptend_all%tauy_top  =0.0_r8!(pcols)! net meridional stress at top of model (Pa)
    ALLOCATE(ptend_all%cflx_srf   (1:ibMax,1:jbMax,1:ppcnst));ptend_all%cflx_srf  =0.0_r8!(pcols,ppcnst)! constituent flux at surface (kg/m2/s)
    ALLOCATE(ptend_all%cflx_top   (1:ibMax,1:jbMax,1:ppcnst));ptend_all%cflx_top  =0.0_r8!(pcols,ppcnst)! constituent flux top of model (kg/m2/s)

!    ALLOCATE(state_concld (1:ibMax,1:kMax,1:jbMax))  ;state_concld =0.0_r8
!    ALLOCATE(state_cld    (1:ibMax,1:kMax,1:jbMax))  ;state_cld    =0.0_r8
   ! ALLOCATE(state_icwmr  (1:ibMax,1:kMax,1:jbMax))  ;state_icwmr  =0.0_r8
   ! ALLOCATE(state_rprddp (1:ibMax,1:kMax,1:jbMax))  ;state_rprddp =0.0_r8
   ! ALLOCATE(state_rprdsh (1:ibMax,1:kMax,1:jbMax))  ;state_rprdsh =0.0_r8
   ! ALLOCATE(state_RPRDTOT(1:ibMax,1:kMax,1:jbMax))  ;state_RPRDTOT =0.0_r8
   ! ALLOCATE(state_cnt    (1:ibMax,1:jbMax))       ;state_cnt    =0.0_r8
   ! ALLOCATE(state_cnb    (1:ibMax,1:jbMax))       ;state_cnb    =0.0_r8
   ! ALLOCATE(state_shfrc(1:ibMax,1:kMax,1:jbMax))    ;state_shfrc =0.0_r8
   ! ALLOCATE(state_evapcsh(1:ibMax,1:kMax,1:jbMax))    ;state_evapcsh =0.0_r8
   ! ALLOCATE(state_cush   (1:ibMax,1:jbMax))       ;state_cush   =0.0_r8
    IF(TRIM(ISCON).EQ.'JHK')THEN
       shallow_scheme='Hack'
    ELSE IF (TRIM(ISCON).EQ.'UW')THEN
       shallow_scheme='UW'
    ELSE
       shallow_scheme='off'
    END IF
    
    DO k=1,kMax+1
!      hypi(k) =  si_in         (kMax+2-k)
!      SB  changed to hybrid (already top to bottom)
       hypi(k) =  a_hybr(k) / 1.e5_r8 + b_hybr(k) 
    END DO
    qmin=1.0e-12_r8
    DO ind=1,ppcnst
       ! set constituent mixing ratio type
       !if ( present(mixtype) )then
       cnst_type(ind) = mixtype(ind) 
       !else
       !   cnst_type(ind) = 'wet'
       !end if
    END DO

    !
    ! Limit shallow convection to regions below 40 mb or 0.040 sigma   40/1000
    ! Note this calculation is repeated in the deep convection interface
    !4.e3 4000
    IF (hypi(1) >= 0.040_r8) THEN
       limcnv = 1
    ELSE
       DO k=1,kMax
          IF (hypi(k) < 0.040 .AND. hypi(k+1) >= 0.040) THEN
             limcnv = k
             GOTO 10
          END IF
       END DO
       limcnv = kMax+1
    END IF

10  CONTINUE

    !if (masterproc) then
    !   write(6,*)'MFINTI: Convection will be capped at intfc ',limcnv, &
    !             ' which is ',hypi(limcnv),' pascals'
    !end if



    CALL mfinti(rair    ,cpair   ,gravit  ,latvap  ,rhoh2o, limcnv,jMax) ! get args from inti.F90
    !
    !-----------------------------------------------------------------------
    !
    ! Specify control parameters first
    !
    ip    = .TRUE.

    CALL gestbl(tmn     ,tmx     ,trice   ,ip      ,epsilo  , &
         latvap  ,latice  ,rh2o    ,cpair   ,tmelt )

    CALL init_uwshcu( r8, latvap, cpair, latice, zvir, rair, gravit, epsilo)

  END SUBROUTINE convect_shallow_init

  !=========================================================================================

  SUBROUTINE convect_shallow_tend( &
       latco   , &!INTEGER, INTENT(IN   )  :: latco
       pcols   , &!INTEGER, INTENT(IN   )  :: pcols
       pver    , &!INTEGER, INTENT(IN   )  :: pver
       pverp   , &!INTEGER, INTENT(IN   )  :: pverp
       nstep   , &!INTEGER, INTENT(IN   )  :: nstep
       pcnst   , &!INTEGER, INTENT(IN   )  :: pcnst
       pnats   , &!INTEGER, INTENT(IN   )  :: pnats
       ztodt   , &!real(r8), intent(in) :: ztodt               ! 2 delta-t (seconds)
       prsi        , &
       prsl        , &
       state_phis  , &! REAL(r8), INTENT(in)  :: state_phis   (pcols)    !(pcols)          ! surface geopotential
       state_t     , &! REAL(r8), INTENT(in)  :: state_t          (pcols,pver)!(pcols,pver)! temperature (K)
       state_qv    , &! REAL(r8), INTENT(in)  :: state_qv          (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       state_ql    , &! REAL(r8), INTENT(in)  :: state_ql          (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       state_qi    , &! REAL(r8), INTENT(in)  :: state_qi          (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       state_u    , &! REAL(r8), INTENT(in)  :: state_ql          (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       state_v    , &! REAL(r8), INTENT(in)  :: state_qi          (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)       
       state_omega , &! REAL(r8), INTENT(in)  :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
       state_concld, &! 
       state_cld   , &! 
       qpert_in   , &!real(r8), intent(inout) :: qpert(pcols)  ! PBL perturbation specific humidity
       pblht   , &!real(r8), intent(in) :: pblht(pcols)        ! PBL height (provided by PBL routine)
       state_cmfmc   , &!real(r8), intent(inout) :: cmfmc(pcols,pverp)  ! moist deep + shallow convection cloud mass flux
       state_cmfmc2  , &!real(r8), intent(out) :: cmfmc2(pcols,pverp)  ! moist shallow convection cloud mass flux
       precc   , &!real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
       state_qc      , &!real(r8), intent(inout) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
       rliq    , &!real(r8), intent(inout) :: rliq(pcols) 
       rliq2   , &!real(r8), intent(out)  :: rliq2(pcols) 
       snow    , &!real(r8), intent(out) :: snow(pcols)  ! snow from this convection
       tke_in  , &  
       ktop    , &!(1:iMax)             , &
       KUO     , &!(1:iMax)             , &
       dtdt   , &
       dqdt   , &
       dqldt  , &
       dqidt  , &
       kctop1  , &
       kcbot1  , &
       noshal1   )
    ! Arguments
    INTEGER, INTENT(IN   )  :: latco
    INTEGER, INTENT(IN   )  :: pcols
    INTEGER, INTENT(IN   )  :: pver
    INTEGER, INTENT(IN   )  :: pverp
    INTEGER, INTENT(IN   )  :: nstep
    INTEGER, INTENT(IN   )  :: pcnst
    INTEGER, INTENT(IN   )  :: pnats
    REAL(r8), INTENT(in) :: ztodt               ! 2 delta-t (seconds)
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in   ) :: prsi    (1:pcols,1:pver+1)  
    REAL(r8), INTENT(in   ) :: prsl    (1:pcols,1:pver)    
    REAL(r8), INTENT(in   ) :: state_phis(pcols)    ! REAL(r8), INTENT(in)  :: state_ps          (pcols)    !(pcols)          ! surface pressure(Pa)
    REAL(r8), INTENT(inout) :: state_t          (pcols,pver)
    REAL(r8), INTENT(inout) :: state_u          (pcols,pver)
    REAL(r8), INTENT(inout) :: state_v          (pcols,pver)
    REAL(r8), INTENT(inout) :: state_qv     (pcols,pver)
    REAL(r8), INTENT(inout) :: state_ql     (pcols,pver)
    REAL(r8), INTENT(inout) :: state_qi     (pcols,pver)
    REAL(r8), INTENT(in)    :: state_omega  (pcols,pver)
    REAL(r8), INTENT(inout) :: state_concld(pcols,pver)
    REAL(r8), INTENT(inout) :: state_cld   (pcols,pver)
    REAL(r8), INTENT(in) :: pblht(pcols)        ! PBL height (provided by PBL routine)
    REAL(r8), INTENT(in) :: qpert_in(pcols)  ! PBL perturbation specific humidity
    REAL(r8) :: qpert(pcols,pcnst+pnats)  ! PBL perturbation specific humidity

    !
    ! fields which combine deep convection with shallow
    !
    REAL(r8), INTENT(inout) :: state_cmfmc(pcols,pverp)  ! moist deep + shallow convection cloud mass flux
    REAL(r8), INTENT(inout) :: state_qc(pcols,pver)      ! dq/dt due to export of cloud water
    REAL(r8), INTENT(inout) :: rliq(pcols) 

    !
    ! Output arguments
    !
    REAL(r8), INTENT(in) :: tke_in(pcols,pver)
    REAL(r8), INTENT(out) :: state_cmfmc2(pcols,pverp)  ! moist shallow convection cloud mass flux
    REAL(r8), INTENT(out) :: precc(pcols)        ! convective precipitation rate
    REAL(r8), INTENT(out)  :: rliq2(pcols) 
    REAL(r8), INTENT(out) :: snow(pcols)  ! snow from this convection   
    INTEGER, INTENT(in   ) :: ktop   ( pcols )
    INTEGER, INTENT(in   ) :: kuo    ( pcols )
    INTEGER, INTENT(OUT   ) :: kctop1(pcols)
    INTEGER, INTENT(OUT   ) :: kcbot1(pcols)
    INTEGER, INTENT(OUT   ) :: noshal1(pcols)
    REAL(KINd=r8)   , INTENT(OUT  ) :: dtdt (pcols,pver)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqdt (pcols,pver)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqldt(pcols,pver)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqidt(pcols,pver)

    ! Local variables

    INTEGER :: i,k,m
    INTEGER :: n,x
    INTEGER :: ilon                      ! global longitude index of a column
    INTEGER :: ilat                      ! global latitude index of a column
    INTEGER :: ncol                 ! number of atmospheric columns
    REAL(r8) :: cmfsl(pcols,pver+1 )  ! convective lw static energy flux
    REAL(r8) :: cmflq(pcols,pver+1 )  ! convective total water flux
    REAL(r8), TARGET :: hkk(pcols)
    REAL(r8), TARGET :: hkl(pcols)
    REAL(r8), TARGET :: tvfac
    REAL(r8) :: cmfmc(pcols,pverp)  ! moist deep + shallow convection cloud mass flux
    REAL(r8) :: cmfmc2(pcols,pverp)  ! moist shallow convection cloud mass flux
    REAL(r8) :: qc(pcols,pver)      ! dq/dt due to export of cloud water

    LOGICAL :: fvdyn
    
    REAL(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
    REAL(r8) :: qc2(pcols,pver)      ! dq/dt due to export of cloud water
    REAL(r8) :: cnt2(pcols)          ! top level of convective activity
    REAL(r8) :: cnb2(pcols)          ! bottom level of convective activity
    REAL(r8) :: tpert(pcols)        ! PBL perturbation theta
    REAL(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
    REAL(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
    REAL(r8) :: flxprec(pcols,pverp)   ! evap outfld: Convective-scale flux of precip at interfaces (kg/m2/s)
    REAL(r8) :: flxsnow(pcols,pverp)   ! evap outfld: Convective-scale flux of snow   at interfaces (kg/m2/s)
    REAL(r8) :: zero(pcols)            ! Array of zeros
    REAL(r8) :: cmfdqs(pcols, pver)    ! Shallow convective snow production
    REAL(r8) :: slflx(pcols,pverp)     ! Shallow convective liquid water static energy flux
    REAL(r8) :: qtflx(pcols,pverp)     ! Shallow convective total water flux
    real(r8) :: iccmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC+IWC [ kg/m2 ]
    real(r8) :: icwmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC     [ kg/m2 ]
    real(r8) :: icimr_UW(pcols,pver)                                      ! In-cloud Cumulus IWC     [ kg/m2 ]
    real(r8) :: ptend_tracer(pcols,pver,pcnst+pnats)                            ! Tendencies of tracers
    REAL(r8) :: cush(pcols)    
    REAL(r8) :: state_pmiddry (pcols,pver)
    REAL(r8) :: state_pdeldry (pcols,pver)
    REAL(r8) :: state_rpdeldry(pcols,pver)
    REAL(r8) :: state_buf_t    (pcols,pver)
    REAL(r8) :: state_buf_u    (pcols,pver)
    REAL(r8) :: state_buf_v    (pcols,pver)
    REAL(r8) :: state_buf_omega(pcols,pver)
    REAL(r8) :: state_buf_pmid  (pcols,pver)  
    REAL(r8) :: state_buf_pint(pcols,pverp)  
    REAL(r8) :: state_buf_lnpint(pcols,pverp)  
    REAL(r8) :: state_buf_pdel(pcols,pver)
    REAL(r8) :: state_buf_rpdel (pcols,pver)
    REAL(r8) :: state_buf_lnpmid  (pcols,pver)
    ! physics types
    !   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
    !   type(physics_ptend) :: ptend_loc     ! local tend from processes, added up to return as ptend_all
    !   type(physics_tend ) :: tend          ! Physics tendencies (empty, needed for physics_update call)

    ! physics buffer fields 
    INTEGER itim, ifld
    REAL(r8) :: concld(pcols,pver)
    REAL(r8) :: cld(pcols,pver)
    REAL(r8) :: icwmr(pcols,pver)    ! in cloud water mixing ratio
    REAL(r8) :: rprddp (pcols,pver)  ! dq/dt due to deep convective rainout
    REAL(r8) :: rprdsh (pcols,pver)  ! dq/dt due to deep and shallow convective rainout
    !REAL(r8) :: RPRDTOT(pcols,pver)
    REAL(r8) :: shfrc(pcols,pver)
    REAL(r8) :: evapcsh(pcols,pver)
    REAL(r8) :: tke(pcols,pverp)  
    REAL(r8) :: cbmf(pcols)                                               ! Shallow cloud base mass flux [ kg/s/m2 ]
    INTEGER :: icheck(pcols)
    INTEGER :: noshal(pcols)

    REAL(r8) :: cnt(pcols) 
    REAL(r8) :: cnb(pcols) 

    IF(pcols < 1) RETURN

   ! ----------------------- !
   ! Main Computation Begins ! 
   ! ----------------------- !
   zero  = 0._r8
   cnt2= 0._r8
   cnb2= 0._r8
   cnt= 0._r8
   cnb= 0._r8

    !---------------------------------------------------------------------- 
    state%ncol(latco)         = pcols
    !state%ps  (1:pcols,latco) = state_ps  (1:pcols) 
    state%ps  (1:pcols,latco) = prsi(1:pcols,1) 
    state%phis(1:pcols,latco) = state_phis(1:pcols) 

    qpert(1:pcols,1)  = qpert_in(1:pcols)  
    !
    !     zero  arrays
    !
    DO i=1,pcols
       noshal(i)=0
       icheck(i)=0
    END DO

    !
    ! initialize 
    !
    DO i=1,pcols
       IF(kuo(i) .EQ. 1) THEN
          noshal(i)=1
       END IF
    END DO

    ncol = state%ncol(latco)

    DO i=1,pcols
       !state_buf_pint       (i,pver+1) = state_ps(i)*si(1)
       state_buf_pint       (i,pver+1)  = prsi(i,1)
    END DO
    DO k=pver,1,-1
       DO i=1,pcols
          !state_buf_pint    (i,k) = MAX(si(pver+2-k)*state_ps(i) ,0.0001_r8)
          state_buf_pint     (i,k) = MAX(prsi(i,pver+2-k) ,0.0001_r8)
       END DO
    END DO

    DO k=1,pver+1
       DO i=1,pcols
          state_buf_lnpint(i,k) =  LOG(state_buf_pint  (i,k))
          cmfmc (i,pver+2-k)     =  state_cmfmc (i,k)
          !cmfmc2(i,pver+2-k)     =  state_cmfmc2(i,k)
       END DO
    END DO

    state%pint     (1:pcols,1:pver+1,latco) = state_buf_pint(1:pcols,1:pver+1)
    state%lnpint   (1:pcols,1:pver+1,latco) = state_buf_lnpint(1:pcols,1:pver+1)

    DO k=1,pver
       DO i=1,pcols

          dtdt  (i,pver+1-k) = 0.0_r8
          dqdt  (i,pver+1-k) = 0.0_r8
          dqldt (i,pver+1-k) = 0.0_r8
          dqidt (i,pver+1-k) = 0.0_r8

          state_buf_t       (i,pver+1-k) = state_t          (i,k)
          state_buf_u       (i,pver+1-k) = state_u          (i,k)
          state_buf_v       (i,pver+1-k) = state_v          (i,k)
          state_buf_omega   (i,pver+1-k) = state_omega      (i,k)
          !state_buf_pmid    (i,pver+1-k) = sl(k)*state_ps (i)
          state_buf_pmid    (i,pver+1-k) = prsl(i,k)
          qc                (i,pver+1-k) = state_qc (i,k) 
       END DO
    END DO

    DO k=1,pver
       DO i=1,pcols          
          state_buf_pdel    (i,k) = MAX(state%pint(i,k+1,latco) - state%pint(i,k,latco),0.5_r8)
          state_buf_rpdel   (i,k) = 1.0_r8/MAX((state%pint(i,k+1,latco) - state%pint(i,k,latco)),0.5_r8)
          state_buf_lnpmid  (i,k) = LOG(state_buf_pmid(i,k))        
       END DO
    END DO
    state%t           (1:pcols,1:pver,latco)= state_buf_t     (1:pcols,1:pver) 
    state%u           (1:pcols,1:pver,latco)= state_buf_u     (1:pcols,1:pver) 
    state%v           (1:pcols,1:pver,latco)= state_buf_v     (1:pcols,1:pver) 
    state%omega    (1:pcols,1:pver,latco)= state_buf_omega (1:pcols,1:pver) 
    state%pmid     (1:pcols,1:pver,latco)= state_buf_pmid  (1:pcols,1:pver) 
    state%pdel     (1:pcols,1:pver,latco)= state_buf_pdel  (1:pcols,1:pver) 
    state%rpdel    (1:pcols,1:pver,latco)= state_buf_rpdel (1:pcols,1:pver) 
    state%lnpmid   (1:pcols,1:pver,latco)= state_buf_lnpmid(1:pcols,1:pver) 

    state_pmiddry (1:pcols,1:pver) = state_buf_pmid  (1:pcols,1:pver) 
    state_pdeldry (1:pcols,1:pver) = state_buf_pdel  (1:pcols,1:pver) 
    state_rpdeldry(1:pcols,1:pver) = state_buf_rpdel (1:pcols,1:pver)

    DO k=1,pver
       DO i=1,pcols
          state%q    (i,pver+1-k,latco,1)         =  state_qv(i,k)
          state%q    (i,pver+1-k,latco,ixcldliq)  =  state_ql(i,k)
          state%q    (i,pver+1-k,latco,ixcldice)  =  state_qi(i,k)
       END DO
    END DO


    ! Derive new temperature and geopotential fields

    CALL geopotential_t(                                 &
         state%lnpint(1:pcols,1:pver+1,latco)   , state%lnpmid(1:pcols,1:pver,latco)   , state%pint (1:pcols,1:pver+1,latco)   , &
         state%pmid  (1:pcols,1:pver,latco)     , state%pdel  (1:pcols,1:pver,latco)   , state%rpdel(1:pcols,1:pver,latco)   , &
         state%t     (1:pcols,1:pver,latco)     , state%q     (1:pcols,1:pver,latco,1) , rair   , gravit , zvir   ,          &
         state%zi    (1:pcols,1:pver+1,latco)   , state%zm    (1:pcols,1:pver,latco)   , ncol   ,pcols, pver, pverp)

    fvdyn = dycore_is ('LR')
    DO k = pver, 1, -1
       ! First set hydrostatic elements consistent with dynamics
       IF (fvdyn) THEN
          DO i = 1,ncol
             hkl(i) = state%lnpmid(i,k+1,latco) - state%lnpmid(i,k,latco)
             hkk(i) = 1.0_r8 - state%pint (i,k,latco)* hkl(i)* state%rpdel(i,k,latco)
          END DO
       ELSE
          DO i = 1,ncol
             hkl(i) = state%pdel(i,k,latco)/state%pmid  (i,k,latco)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF
       ! Now compute s
       DO i = 1,ncol
          tvfac   = 1.0_r8 + zvir * state%q(i,k,latco,1) 
          state%s(i,k,latco) =  (state%t(i,k,latco)* cpair) + (state%t(i,k,latco) * tvfac * rair*hkk(i))  +  &
               ( state%phis(i,latco) + gravit*state%zi(i,k+1,latco))           
       END DO
    END DO

    CALL physics_state_copy(pcols,latco,pver, pverp,ppcnst, cnst_need_pdeldry)   ! copy state to local state1.
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)  ! initialize local ptend type
    CALL physics_ptend_init(ptend_all,latco,ppcnst,pver)  ! initialize output ptend type
    CALL physics_tend_init(tend,latco)                    ! tend here is just a null place holder

    !
    ! Associate pointers with physics buffer fields
    !

    !itim = pbuf_old_tim_idx()
    !ifld = pbuf_get_fld_idx('CLD')
    !cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    !cld (1:pcols,1:pver)=state_cld (1:pcols,1:pver,latco)

    !ifld = pbuf_get_fld_idx('ICWMRSH')
    !icwmr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    !icwmr(1:pcols,1:pver)=state_icwmr (1:pcols,1:pver,latco)

    !ifld = pbuf_get_fld_idx('RPRDDP')
    !rprddp => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    !rprddp(1:pcols,1:pver)=state_rprddp (1:pcols,1:pver,latco)

    !ifld = pbuf_get_fld_idx('RPRDSH')
    !rprdsh => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    !rprdsh(1:pcols,1:pver)=state_rprdsh (1:pcols,1:pver,latco)

    !ifld = pbuf_get_fld_idx('CLDTOP')
    !cnt => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,1)
    !cnt(1:pcols)=state_cnt (1:pcols,latco)

    !ifld = pbuf_get_fld_idx('CLDBOT')
    !cnb => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,1)
    !cnt(1:pcols)=state_cnb (1:pcols,latco)
    !itim   =  pbuf_old_tim_idx()
    !ifld   =  pbuf_get_fld_idx('CONCLD')
    !concld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

    DO k=1,pver
       DO i=1,pcols
          concld(i,pver+1-k)=  state_concld (i,k)  
          cld   (i,pver+1-k)=  state_cld    (i,k)  
          !icwmr (i,pver+1-k)=  state_icwmr  (i,k,latco)  
          !rprddp(i,pver+1-k)=  state_rprddp (i,k,latco)  
          !rprdsh(i,pver+1-k)=  state_rprdsh (i,k,latco)  
          !RPRDTOT(i,pver+1-k)= state_RPRDTOT(i,k,latco)  
          !shfrc  (i,pver+1-k)= state_shfrc  (i,k,latco)  
          !evapcsh(i,pver+1-k)= state_evapcsh(i,k,latco)  
       END DO
    END DO
    !DO i=1,pcols
    !   cnt(i) = state_cnt (i,latco) 
    !   cnb(i) = state_cnb (i,latco) 
    !END DO

    ! Initialize

    tpert(1:ncol  ) =0.0_r8
    qpert(1:ncol,2:pcnst+pnats) = 0.0_r8

   select case (shallow_scheme)
   case('off') ! None
      cmfmc2      = 0._r8
      ptend_loc%q = 0._r8
      ptend_loc%s = 0._r8
      rprdsh      = 0._r8
      rprddp      = 0._r8
      cmfdqs      = 0._r8
      precc       = 0._r8
      slflx       = 0._r8
      qtflx       = 0._r8
      icwmr       = 0._r8
      rliq2       = 0._r8
      qc2         = 0._r8
      cmfsl       = 0._r8
      cmflq       = 0._r8
      cnt2        = 0._r8
      cnb2        = 0._r8
      evapcsh     = 0._r8   
   case('Hack') ! Hack scheme
    !
    !  Call Hack convection 
    !
    CALL cmfmca (&
         ncol                                            , &!integer, intent(in) :: ncol          ! number of atmospheric columns
         pcols                                                , &!integer, intent(in) :: pcols
         pver                                                , &!integer, intent(in) :: pver
         pverp                                                , &!integer, intent(in) :: pverp
         pcnst                                                , &!integer, intent(in) :: pcnst
         pnats                                                , &!integer, intent(in) :: pnats
         nstep                                           , &!integer, intent(in) :: nstep         ! current time step index
         ztodt                                           , &!real(r8), intent(in) :: ztodt               ! 2 delta-t (seconds)
         state%pmid (1:pcols,1:pver,latco)               , &!real(r8), intent(in) :: pmid(pcols,pver)    ! pressure
         state%pdel (1:pcols,1:pver,latco)               , &!real(r8), intent(in) :: pdel(pcols,pver)    ! delta-p
         state%rpdel(1:pcols,1:pver,latco)               , &!real(r8), intent(in) :: rpdel(pcols,pver)   ! 1./pdel
         state%zm   (1:pcols,1:pver,latco)               , &!real(r8), intent(in) :: zm(pcols,pver)      ! height abv sfc at midpoints
         tpert      (1:ncol  )                           , &!real(r8), intent(in) :: tpert(pcols)      ! PBL perturbation theta
         qpert      (1:ncol,1:pcnst+pnats)               , &!real(r8), intent(in) :: qpert(pcols,pcnst+pnats)  ! PBL perturbation specific humidity
         state%phis (1:pcols,latco)                      , &!real(r8), intent(in) :: phis(pcols)            ! surface geopotential
         pblht      (1:ncol  )                           , &!real(r8), intent(in) :: pblht(pcols)      ! PBL height (provided by PBL routine)
         state%t    (1:pcols,1:pver,latco)               , &!real(r8), intent(in) :: t(pcols,pver)     ! temperature (t bar)
         state%q    (1:pcols,1:pver,latco,1:pcnst+pnats) , &!real(r8), intent(in) :: q(pcols,pver,pcnst+pnats) ! specific humidity (sh bar)
         ptend_loc%s(1:pcols,1:pver,latco)               , &!real(r8), intent(out) :: cmfdt(pcols,pver)   ! dt/dt due to moist convection
         ptend_loc%q(1:pcols,1:pver,latco,1:pcnst+pnats) , &!real(r8), intent(out) :: dq(pcols,pver,pcnst+pnats) ! constituent tendencies
         cmfmc2     (1:pcols,1:pverp)                    , &!real(r8), intent(out) :: cmfmc(pcols,pverp)  ! moist convection cloud mass flux
         rprdsh     (1:pcols,1:pver)                     , &!real(r8), intent(out) :: cmfdqr(pcols,pver)  ! dq/dt due to convective rainout
         cmfsl      (1:pcols,1:pver+1)                     , &!real(r8), intent(out) :: cmfsl(pcols,pver )  ! convective lw static energy flux
         cmflq      (1:pcols,1:pver+1)                     , &!real(r8), intent(out) :: cmflq(pcols,pver )  ! convective total water flux
         precc      (1:pcols)                            , &!real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
         qc2        (1:pcols,1:pver)                     , &!real(r8), intent(out) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
         cnt2       (1:pcols)                            , &!real(r8), intent(out) :: cnt(pcols)               ! top level of convective activity
         cnb2       (1:pcols)                            , &!real(r8), intent(out) :: cnb(pcols)               ! bottom level of convective activity
         icwmr      (1:pcols,1:pver)                     , &!real(r8), intent(out) :: icwmr(pcols,pver)
         rliq2      (1:pcols)                            , &!real(r8), intent(out) :: rliq(pcols) 
         state_pmiddry(1:pcols,1:pver)                   , &!real(r8), intent(in) :: pmiddry(pcols,pver)         ! pressure
         state_pdeldry(1:pcols,1:pver)                   , &!real(r8), intent(in) :: pdeldry(pcols,pver)         ! delta-p
         state_rpdeldry(1:pcols,1:pver)                  , &!real(r8), intent(in) :: rpdeldry(pcols,pver)   ! 1./pdel 
         kctop1     (1:pcols)                            , &
         kcbot1     (1:pcols)                            , &
         noshal     (1:pcols)) 
   case('UW')   ! UW shallow convection scheme

      !cush(1:pcols) =state_cush(1:pcols,latco)! => pbuf(pbuf_get_fld_idx('cush'))%fld_ptr(1,1:pcols,1,lchnk,itim)
      tke(1:pcols,pver+1) = 0.01_r8*0.8_r8   +  0.2_r8*tke_in(1:pcols,1)

      DO k=1,pver
         DO i=1,pcols
            tke(i,pver+1-k)=tke_in(i,k) ! => pbuf(pbuf_get_fld_idx('tke'))%fld_ptr(1,1:pcols,1:pverp,lchnk,itim)
         END DO
      END DO            
      if( nstep .le. 1 ) then
          cush(:)  = 1.e3_r8
          tke(:,:) = 0.01_r8
      end if

      call compute_uwshcu_inv( pcols                                         , &!  mix                , &
                               pver                                          , &!  mkx                , &
                               ncol                                          , &!  iend          , &
                               pcnst+pnats                                   , &!  ncnst         , &
                               ztodt                                         , &!  dt                , & 
                               state%pint   (1:pcols,1:pver+1,latco)         , &!  ps0_inv        , &
                               state%zi     (1:pcols,1:pver+1,latco)         , &!  zs0_inv        , &
                               state%pmid   (1:pcols,1:pver  ,latco)         , &!  p0_inv        , &
                               state%zm     (1:pcols,1:pver  ,latco)         , &!  z0_inv        , &
                               state%pdel   (1:pcols,1:pver  ,latco)         , &!  dp0_inv        , &
                               state%u      (1:pcols,1:pver  ,latco)         , &!  u0_inv        , &
                               state%v      (1:pcols,1:pver  ,latco)         , &!  v0_inv        , &
                               state%q      (1:pcols,1:pver  ,latco,1       ), &!  qv0_inv        , &
                               state%q      (1:pcols,1:pver  ,latco,ixcldliq), &!  ql0_inv        , &
                               state%q      (1:pcols,1:pver  ,latco,ixcldice), &!  qi0_inv        , &
                               state%t      (1:pcols,1:pver  ,latco)         , &!  t0_inv        , &
                               state%s      (1:pcols,1:pver  ,latco)         , &!  s0_inv        , &
                               state%q      (1:pcols,1:pver  ,latco,1:pcnst+pnats), &!  tr0_inv        , &
                               tke          (1:pcols,1:pver+1               ), &!  tke_inv        , &
                               cld          (1:pcols,1:pver                 ), &!  cldfrct_inv        , &
                               concld       (1:pcols,1:pver                 ), &!  concldfrct_inv, &
                               pblht        (1:ncol  )                       , &!  pblh          , &
                               cush         (1:ncol  )                       , &!  cush          , & 
                               cmfmc2       (1:pcols,1:pverp)                , &!  umf_inv        , &
                               slflx        (1:pcols,1:pverp)                , &!  slflx_inv        , &
                               qtflx        (1:pcols,1:pverp)                , &!  qtflx_inv        , & 
                               ptend_loc%q  (1:pcols,1:pver  ,latco,1       ), &!  qvten_inv        , &
                               ptend_loc%q  (1:pcols,1:pver  ,latco,ixcldliq), &!  qlten_inv        , &
                               ptend_loc%q  (1:pcols,1:pver  ,latco,ixcldice), &!  qiten_inv        , &
                               ptend_loc%s  (1:ncol ,1:pver  ,latco         ), &!  sten_inv        , &
                               ptend_loc%u  (1:pcols,1:pver  ,latco         ), &!  uten_inv        , &
                               ptend_loc%v  (1:pcols,1:pver  ,latco         ), &!  vten_inv        , &
                               ptend_tracer (1:pcols,1:pver  ,1:pcnst+pnats ), &!  trten_inv        , &  
                               rprdsh       (1:pcols,1:pver                 ), &!  qrten_inv        , &
                               cmfdqs       (1:pcols,1:pver)                 , &!  qsten_inv        , &
                               precc        (1:pcols)                        , &!  precip        , &
                               snow         (1:pcols)                        , &!  snow          , &
                               evapcsh      (1:pcols,1:pver)                 , &!  evapc_inv        , &
                               shfrc        (1:pcols,1:pver)                 , &!  cufrc_inv        , &
                               iccmr_UW     (1:pcols,1:pver)                 , &!  qcu_inv        , &
                               icwmr_UW     (1:pcols,1:pver)                 , &!  qlu_inv        , &
                               icimr_UW     (1:pcols,1:pver)                 , &!  qiu_inv        , &   
                               cbmf         (1:pcols)                        , &!  cbmf          , &
                               qc2          (1:pcols,1:pver)                 , &!  qc_inv        , &
                               rliq2        (1:pcols)                        , &!  rliq          , &
                               cnt2         (1:pcols)                             , &!  cnt_inv        , &
                               cnb2         (1:pcols)                             , &!  cnb_inv        , &
                               fqsatd                                        , &!  qsat          , &
                               state_pdeldry(1:pcols,1:pver)                   )!  dpdry0_inv           ) 
      DO i=1,pcols
         precc(i)=precc(i)-snow(i)
         kctop1(i) = pver-cnt2(i)
         kcbot1(i) = pver-cnb2(i)
         IF(precc(i)<=0.0_r8)THEN
            noshal(i) = 1
         END IF
      END DO


      ! state_cush(1:pcols,latco)=cush(1:pcols) ! => pbuf(pbuf_get_fld_idx('cush'))%fld_ptr(1,1:pcols,1,lchnk,itim)

      ! --------------------------------------------------------------------- !
      ! Here, 'rprdsh = qrten', 'cmfdqs = qsten' both in unit of [ kg/kg/s ]  !
      ! In addition, define 'icwmr' which includes both liquid and ice.       !
      ! --------------------------------------------------------------------- !

      icwmr(1:pcols,1:pver)  = iccmr_UW(1:pcols,1:pver) 
      rprdsh(1:pcols,1:pver) = rprdsh(1:pcols,1:pver) + cmfdqs(1:pcols,1:pver)
      do m = 4, pcnst
         IF(m > 3)THEN
            ptend_loc%q(1:pcols,1:pver,latco,m) = ptend_tracer(1:ncol,1:pver,m)
         END IF
      enddo


      ! Conservation check
      
      !  do i = 1, ncol
      !  do m = 1, pcnst
      !     sum1 = 0._r8
      !     sum2 = 0._r8
      !     sum3 = 0._r8
      !  do k = 1, pver
      !       if(cnst_get_type_byind(m).eq.'wet') then
      !          pdelx = state%pdel(i,k)
      !       else
      !          pdelx = state%pdeldry(i,k)
      !       endif
      !       sum1 = sum1 + state%q(i,k,m)*pdelx
      !       sum2 = sum2 +(state%q(i,k,m)+ptend_loc%q(i,k,m)*ztodt)*pdelx  
      !       sum3 = sum3 + ptend_loc%q(i,k,m)*pdelx 
      !  enddo
      !  if( m .gt. 3 .and. abs(sum1) .gt. 1.e-13_r8 .and. abs(sum2-sum1)/sum1 .gt. 1.e-12_r8 ) then
      !! if( m .gt. 3 .and. abs(sum3) .gt. 1.e-13_r8 ) then
      !      write(iulog,*) 'Sungsu : convect_shallow.F90 does not conserve tracers : ', m, sum1, sum2, abs(sum2-sum1)/sum1
      !!     write(iulog,*) 'Sungsu : convect_shallow.F90 does not conserve tracers : ', m, sum3
      !  endif
      !  enddo
      !  enddo

      ! ------------------------------------------------- !
      ! Convective fluxes of 'sl' and 'qt' in energy unit !
      ! ------------------------------------------------- !

      cmfsl(1:ncol,1:pverp) = slflx(1:ncol,1:pverp)
      cmflq(1:ncol,1:pverp) = qtflx(1:ncol,1:pverp) * latvap



      ! -------------------------------------- !
      ! uwshcu does momentum transport as well !
      ! -------------------------------------- !

      ptend_loc%lu = .TRUE.
      ptend_loc%lv = .TRUE.

      !call outfld( 'PRECSH' , precc  , pcols, lchnk )

   end select

    ptend_loc%name(latco)  = 'cmfmca'
    ptend_loc%ls(latco)    = .TRUE.
    ptend_loc%lq(latco,:)  = .TRUE.

    !
    ! Merge shallow/mid-level output with prior results from deep
    !

    !combine cmfmc2 (from shallow)  with cmfmc (from deep)
    DO k=1,pver+1
       DO i=1,pcols
         IF(noshal(i) == 0)THEN
            cmfmc (i,k) = cmfmc (i,k) + cmfmc2 (i,k)
         ELSE
            cmfmc2 (i,k)=0.0_r8
         END IF
      END DO
   END DO         
  ! -------------------------------------------------------------- !
   ! 'cnt2' & 'cnb2' are from shallow, 'cnt' & 'cnb' are from deep  !
   ! 'cnt2' & 'cnb2' are the interface indices of cloud top & base: ! 
   !        cnt2 = float(kpen)                                      !
   !        cnb2 = float(krel - 1)                                  !
   ! Note that indices decreases with height.                       !
   ! -------------------------------------------------------------- !
    DO i=1,ncol
       IF (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
       IF (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
    END DO

    ! This quantity was previously known as CMFDQR.  Now CMFDQR is the shallow rain production only.
    !ifld = pbuf_get_fld_idx( 'RPRDTOT' )
    !pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = rprdsh(:ncol,:pver) + rprddp(:ncol,:pver)
    !RPRDTOT(1:ncol,1:pver)= rprdsh(1:ncol,1:pver) + rprddp(1:ncol,1:pver)
    ! Add shallow cloud water detrainment to cloud water detrained from deep convection
    DO k=1,pver
       DO i=1,pcols
          IF(noshal(i) == 0)THEN
              qc(i,k) = qc(i,k) + qc2(i,k)
          ELSE
              qc2(i,k)=0.0_r8
         END IF
      END DO
   END DO         
    DO i=1,pcols
       IF(noshal(i) == 0)THEN
          rliq(i) = rliq(i) + rliq2(i)    
       ELSE
          rliq2(i)=0.0_r8
       END IF
   END DO

    ftem(1:ncol,1:pver) = ptend_loc%s(1:ncol,:pver,latco)/cpair

    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, state,ppcnst,latco)

    ! update physics state type state1 with ptend_loc 
    CALL physics_update(state1, tend, ptend_loc,ztodt,ppcnst,pcols,pver,latco)

    ! initialize ptend for next process
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)
 
   ! ------------------------------------------------------------------------ !
   ! UW-Shallow Cumulus scheme includes                                       !
   ! evaporation physics inside in it. So when 'shallow_scheme = UW', we must !
   ! NOT perform below 'zm_conv_evap'.                                        !
   ! ------------------------------------------------------------------------ !

    if( shallow_scheme .eq. 'Hack' ) then

    !
    ! Determine the phase of the precipitation produced and add latent heat of fusion
    ! Evaporate some of the precip directly into the environment (Sundqvist)
    ! Allow this to use the updated state1 and a fresh ptend_loc type
    ! Heating and specific humidity tendencies produced
    !
    ptend_loc%name(latco)  = 'zm_conv_evap'
    ptend_loc%ls(latco)    = .TRUE.
    ptend_loc%lq(latco,1)  = .TRUE.

    CALL zm_conv_evap( &
         ncol                               , &!INTENT(in) :: ncol               ! number of columns and chunk index
         pcols                              , &!INTENT(in) :: pcols                    ! number of columns (max)
         pver                               , &!INTENT(in) :: pver                    ! number of vertical levels
         pverp                              , &!INTENT(in) :: pverp                    ! number of vertical levels + 1        
         state1%t   (1:pcols,1:pver,latco)  , &!INTENT(in), DIMENSION(pcols,pver) :: t    ! temperature (K)
         state1%pmid(1:pcols,1:pver,latco)  , &!INTENT(in), DIMENSION(pcols,pver) :: pmid         ! midpoint pressure (Pa) 
         state1%pdel(1:pcols,1:pver,latco)  , &!INTENT(in), DIMENSION(pcols,pver) :: pdel         ! layer thickness (Pa)
         state1%q   (1:pcols,1:pver,latco,1), &!INTENT(in), DIMENSION(pcols,pver) :: q          ! water vapor (kg/kg)
         ptend_loc%s(1:pcols,1:pver,latco)  , &!INTENT(inout), DIMENSION(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
         ptend_loc%q(1:pcols,1:pver,latco,1), &!INTENT(inout), DIMENSION(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
         rprdsh     (1:pcols,1:pver)        , &!INTENT(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
         cld        (1:pcols,1:pver)        , &!INTENT(in   ) :: cldfrc(pcols,pver) ! cloud fraction
         ztodt                              , &!INTENT(in   ) :: deltat             ! time step
         precc      (1:pcols)               , &!INTENT(inout) :: prec(pcols)             ! Convective-scale preciptn rate
         snow       (1:pcols)               , &!INTENT(out)   :: snow(pcols)             ! Convective-scale snowfall rate
         ntprprd    (1:pcols,1:pver)        , &!INTENT(out) :: ntprprd(pcols,pver)    ! net precip production in layer
         ntsnprd    (1:pcols,1:pver)        , &!INTENT(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
         flxprec    (1:pcols,1:pver+1)      , &!INTENT(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
         flxsnow    (1:pcols,1:pver+1  )      )!INTENT(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)

    ! record history variables from zm_conv_evap

    ftem(1:ncol,1:pver) = ptend_loc%s(1:ncol,1:pver,latco)/cpair

    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, state,ppcnst,latco)

    ! -------------------------------------------- !
    ! Do not perform evaporation process for UW-Cu !
    ! -------------------------------------------- !

 
    end if

    ! update name of parameterization tendencies to send to tphysbc
    ptend_all%name(latco) = 'convect_shallow'

    CALL physics_update   (state, tend, ptend_all, ztodt,ppcnst,pcols,pver,latco)
    noshal1=noshal
    DO k=1,pver
       DO i=1,pcols
         IF(noshal(i) == 0)THEN

            dtdt  (i,pver+1-k) =  (state%t(i,k,latco         ) - state_t   (i,pver+1-k))/ztodt
            dqdt  (i,pver+1-k) =  (state%q(i,k,latco,1       ) - state_qv  (i,pver+1-k))/ztodt
            dqldt (i,pver+1-k) =  (state%q(i,k,latco,ixcldliq) - state_ql  (i,pver+1-k))/ztodt
            dqidt (i,pver+1-k) =  (state%q(i,k,latco,ixcldice) - state_qi  (i,pver+1-k))/ztodt


            state_t   (i,pver+1-k)= state%t(i,k,latco)
            state_qv  (i,pver+1-k)= state%q(i,k,latco,1)
            state_ql  (i,pver+1-k)= state%q(i,k,latco,ixcldliq)
            state_qi  (i,pver+1-k)= state%q(i,k,latco,ixcldice)
            state_qc  (i,k)       = qc(i,pver+1-k)  
         END IF
            state_concld (i,k) =  concld (i,pver+1-k)
            state_cld    (i,k) =  cld   (i,pver+1-k)
            !state_icwmr  (i,k,latco) =  icwmr (i,pver+1-k)
            !state_rprddp (i,k,latco) =  rprddp(i,pver+1-k)
            !state_rprdsh (i,k,latco) =  rprdsh(i,pver+1-k)
            !state_RPRDTOT(i,k,latco) =  RPRDTOT(i,pver+1-k)
            !state_shfrc  (i,k,latco) =  shfrc(i,pver+1-k)
            !state_evapcsh(i,k,latco) =  evapcsh(i,pver+1-k)
       END DO
    END DO
    DO k=1,pver+1
       DO i=1,pcols
         IF(noshal(i) == 0)THEN
            state_cmfmc (i,k) = cmfmc (i,pver+2-k)
            state_cmfmc2(i,k) = cmfmc2(i,pver+2-k)
         END IF
       END DO
    END DO
    DO i=1,pcols
       IF(noshal(i) == 1)THEN
         precc(pcols)=0.0_r8        ! convective precipitation rate
         rliq2(pcols) =0.0_r8
         snow(pcols) =0.0_r8 ! snow from this convection   
        END IF
    END DO

    !DO i=1,pcols
    !     state_cnt (i,latco) = cnt(i)
    !     state_cnb (i,latco) = cnb(i)
    !END DO
  END SUBROUTINE convect_shallow_tend
  !=========================================================================================

  SUBROUTINE cmfmca( &
       ncol    , &!integer, intent(in) :: ncol                 ! number of atmospheric columns
       pcols   , &!integer, intent(in) :: pcols
       pver    , &!integer, intent(in) :: pver
       pverp   , &!integer, intent(in) :: pverp
       pcnst   , &!integer, intent(in) :: pcnst
       pnats   , &!integer, intent(in) :: pnats
       nstep   , &!integer, intent(in) :: nstep                ! current time step index
       ztodt   , &!real(r8), intent(in) :: ztodt               ! 2 delta-t (seconds)
       pmid    , &!real(r8), intent(in) :: pmid(pcols,pver)    ! pressure
       pdel    , &!real(r8), intent(in) :: pdel(pcols,pver)    ! delta-p
       rpdel   , &!real(r8), intent(in) :: rpdel(pcols,pver)   ! 1./pdel
       zm      , &!real(r8), intent(in) :: zm(pcols,pver)      ! height abv sfc at midpoints
       tpert   , &!real(r8), intent(in) :: tpert(pcols)            ! PBL perturbation theta
       qpert   , &!real(r8), intent(in) :: qpert(pcols,pcnst+pnats)  ! PBL perturbation specific humidity
       phis    , &!real(r8), intent(in) :: phis(pcols)            ! surface geopotential
       pblht   , &!real(r8), intent(in) :: pblht(pcols)            ! PBL height (provided by PBL routine)
       t       , &!real(r8), intent(in) :: t(pcols,pver)            ! temperature (t bar)
       q       , &!real(r8), intent(in) :: q(pcols,pver,pcnst+pnats) ! specific humidity (sh bar)
       cmfdt   , &!real(r8), intent(out) :: cmfdt(pcols,pver)   ! dt/dt due to moist convection
       dq      , &!real(r8), intent(out) :: dq(pcols,pver,pcnst+pnats) ! constituent tendencies
       cmfmc   , &!real(r8), intent(out) :: cmfmc(pcols,pverp)  ! moist convection cloud mass flux
       cmfdqr  , &!real(r8), intent(out) :: cmfdqr(pcols,pver)  ! dq/dt due to convective rainout
       cmfsl   , &!real(r8), intent(out) :: cmfsl(pcols,pver )  ! convective lw static energy flux
       cmflq   , &!real(r8), intent(out) :: cmflq(pcols,pver )  ! convective total water flux
       precc   , &!real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
       qc      , &!real(r8), intent(out) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
       cnt     , &!real(r8), intent(out) :: cnt(pcols)          ! top level of convective activity
       cnb     , &!real(r8), intent(out) :: cnb(pcols)          ! bottom level of convective activity
       icwmr   , &!real(r8), intent(out) :: icwmr(pcols,pver)
       rliq    , &! real(r8), intent(out) :: rliq(pcols) 
       pmiddry , &!real(r8), intent(in) :: pmiddry(pcols,pver)    ! pressure
       pdeldry , &!real(r8), intent(in) :: pdeldry(pcols,pver)    ! delta-p
       rpdeldry, &!real(r8), intent(in) :: rpdeldry(pcols,pver)   ! 1./pdel
       kctop1  , &
       kcbot1  , &
       noshal   )

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Moist convective mass flux procedure:
    ! 
    ! Method: 
    ! If stratification is unstable to nonentraining parcel ascent,
    ! complete an adjustment making successive use of a simple cloud model
    ! consisting of three layers (sometimes referred to as a triplet)
    !
    ! Code generalized to allow specification of parcel ("updraft")
    ! properties, as well as convective transport of an arbitrary
    ! number of passive constituents (see q array).  The code
    ! is written so the water vapor field is passed independently
    ! in the calling list from the block of other transported
    ! constituents, even though as currently designed, it is the
    ! first component in the constituents field.
    ! 
    ! Author: J. Hack
    !
    ! BAB: changed code to report tendencies in cmfdt and dq, instead of
    ! updating profiles. Cmfdq contains water only, made it a local variable
    ! made dq (all constituents) the argument.
    ! 
    !-----------------------------------------------------------------------

    !#######################################################################
    !#                                                                     #
    !# Debugging blocks are marked this way for easy identification        #
    !#                                                                     #
    !#######################################################################

    REAL(r8) ,PARAMETER :: ssfac = 1.001              ! supersaturation bound (detrained air)

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                 ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    INTEGER, INTENT(in) :: pverp
    INTEGER, INTENT(in) :: pcnst
    INTEGER, INTENT(in) :: pnats
    INTEGER, INTENT(in) :: nstep                ! current time step index

    REAL(r8), INTENT(in) :: ztodt               ! 2 delta-t (seconds)
    REAL(r8), INTENT(in) :: pmid(pcols,pver)    ! pressure
    REAL(r8), INTENT(in) :: pdel(pcols,pver)    ! delta-p
    REAL(r8), INTENT(in) :: pmiddry(pcols,pver)    ! pressure
    REAL(r8), INTENT(in) :: pdeldry(pcols,pver)    ! delta-p
    REAL(r8), INTENT(in) :: rpdel(pcols,pver)   ! 1./pdel
    REAL(r8), INTENT(in) :: rpdeldry(pcols,pver)   ! 1./pdel
    REAL(r8), INTENT(in) :: zm(pcols,pver)      ! height abv sfc at midpoints
    REAL(r8), INTENT(in) :: tpert(pcols)        ! PBL perturbation theta
    REAL(r8), INTENT(in) :: qpert(pcols,pcnst+pnats)  ! PBL perturbation specific humidity
    REAL(r8), INTENT(in) :: phis(pcols)         ! surface geopotential
    REAL(r8), INTENT(in) :: pblht(pcols)        ! PBL height (provided by PBL routine)
    REAL(r8), INTENT(in) :: t(pcols,pver)       ! temperature (t bar)
    REAL(r8), INTENT(in) :: q(pcols,pver,pcnst+pnats) ! specific humidity (sh bar)
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: cmfdt(pcols,pver)   ! dt/dt due to moist convection
    REAL(r8), INTENT(out) :: cmfmc(pcols,pverp)  ! moist convection cloud mass flux
    REAL(r8), INTENT(out) :: cmfdqr(pcols,pver)  ! dq/dt due to convective rainout
    REAL(r8), INTENT(out) :: cmfsl(pcols,pver+1 )  ! convective lw static energy flux
    REAL(r8), INTENT(out) :: cmflq(pcols,pver+1 )  ! convective total water flux
    REAL(r8), INTENT(out) :: precc(pcols)        ! convective precipitation rate
    ! JJH mod to explicitly export cloud water
    REAL(r8), INTENT(out) :: qc(pcols,pver)      ! dq/dt due to export of cloud water
    REAL(r8), INTENT(out) :: cnt(pcols)          ! top level of convective activity
    REAL(r8), INTENT(out) :: cnb(pcols)          ! bottom level of convective activity
    REAL(r8), INTENT(out) :: dq(pcols,pver,pcnst+pnats) ! constituent tendencies
    REAL(r8), INTENT(out) :: icwmr(pcols,pver)
    REAL(r8), INTENT(out) :: rliq(pcols) 
    INTEGER, INTENT(OUT  ) :: kctop1(pcols) 
    INTEGER, INTENT(OUT  ) :: kcbot1(pcols) 
    INTEGER, INTENT(INOUT) :: noshal(pcols) 
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(r8) pm(pcols,pver)    ! pressure
    REAL(r8) pd(pcols,pver)    ! delta-p
    REAL(r8) rpd(pcols,pver)   ! 1./pdel

    REAL(r8) cmfdq(pcols,pver)   ! dq/dt due to moist convection
    REAL(r8) gam(pcols,pver)     ! 1/cp (d(qsat)/dT)
    REAL(r8) sb(pcols,pver)      ! dry static energy (s bar)
    REAL(r8) hb(pcols,pver)      ! moist static energy (h bar)
    REAL(r8) shbs(pcols,pver)    ! sat. specific humidity (sh bar star)
    REAL(r8) hbs(pcols,pver)     ! sat. moist static energy (h bar star)
    REAL(r8) shbh(pcols,pverp)   ! specific humidity on interfaces
    REAL(r8) sbh(pcols,pverp)    ! s bar on interfaces
    REAL(r8) hbh(pcols,pverp)    ! h bar on interfaces
    REAL(r8) cmrh(pcols,pverp)   ! interface constituent mixing ratio
    REAL(r8) prec(pcols)         ! instantaneous total precipitation
    REAL(r8) dzcld(pcols)        ! depth of convective layer (m)
    REAL(r8) beta(pcols)         ! overshoot parameter (fraction)
    REAL(r8) betamx(pcols)       ! local maximum on overshoot
    REAL(r8) eta(pcols)          ! convective mass flux (kg/m^2 s)
    REAL(r8) etagdt(pcols)       ! eta*grav*dt
    REAL(r8) cldwtr(pcols)       ! cloud water (mass)
    REAL(r8) rnwtr(pcols)        ! rain water  (mass)
    !  JJH extension to facilitate export of cloud liquid water
    REAL(r8) totcond(pcols)        ! total condensate; mix of precip and cloud water (mass)
    REAL(r8) sc  (pcols)         ! dry static energy   ("in-cloud")
    REAL(r8) shc (pcols)         ! specific humidity   ("in-cloud")
    REAL(r8) hc  (pcols)         ! moist static energy ("in-cloud")
    REAL(r8) cmrc(pcols)         ! constituent mix rat ("in-cloud")
    REAL(r8) dq1(pcols)          ! shb  convective change (lower lvl)
    REAL(r8) dq2(pcols)          ! shb  convective change (mid level)
    REAL(r8) dq3(pcols)          ! shb  convective change (upper lvl)
    REAL(r8) ds1(pcols)          ! sb   convective change (lower lvl)
    REAL(r8) ds2(pcols)          ! sb   convective change (mid level)
    REAL(r8) ds3(pcols)          ! sb   convective change (upper lvl)
    REAL(r8) dcmr1(pcols)        ! q convective change (lower lvl)
    REAL(r8) dcmr2(pcols)        ! q convective change (mid level)
    REAL(r8) dcmr3(pcols)        ! q convective change (upper lvl)
    REAL(r8) estemp(pcols,pver)  ! saturation vapor pressure (scratch)
    REAL(r8) vtemp1(2*pcols)     ! intermediate scratch vector
    REAL(r8) vtemp2(2*pcols)     ! intermediate scratch vector
    REAL(r8) vtemp3(2*pcols)     ! intermediate scratch vector
    REAL(r8) vtemp4(2*pcols)     ! intermediate scratch vector
    INTEGER indx1(pcols)     ! longitude indices for condition true
    LOGICAL etagt0           ! true if eta > 0.0
    REAL(r8) sh1                 ! dummy arg in qhalf statement func.
    REAL(r8) sh2                 ! dummy arg in qhalf statement func.
    REAL(r8) shbs1               ! dummy arg in qhalf statement func.
    REAL(r8) shbs2               ! dummy arg in qhalf statement func.
    REAL(r8) cats                ! modified characteristic adj. time
    REAL(r8) rtdt                ! 1./ztodt
    REAL(r8) qprime              ! modified specific humidity pert.
    REAL(r8) tprime              ! modified thermal perturbation
    REAL(r8) pblhgt              ! bounded pbl height (max[pblh,1m])
    REAL(r8) fac1                ! intermediate scratch variable
    REAL(r8) shprme              ! intermediate specific humidity pert.
    REAL(r8) qsattp              ! sat mix rat for thermally pert PBL parcels
    REAL(r8) dz                  ! local layer depth
    REAL(r8) temp1               ! intermediate scratch variable
    REAL(r8) b1                  ! bouyancy measure in detrainment lvl
    REAL(r8) b2                  ! bouyancy measure in condensation lvl
    REAL(r8) temp2               ! intermediate scratch variable
    REAL(r8) temp3               ! intermediate scratch variable
    REAL(r8) g                   ! bounded vertical gradient of hb
    REAL(r8) tmass               ! total mass available for convective exch
    REAL(r8) denom               ! intermediate scratch variable
    REAL(r8) qtest1              ! used in negative q test (middle lvl)
    REAL(r8) qtest2              ! used in negative q test (lower lvl)
    REAL(r8) fslkp               ! flux lw static energy (bot interface)
    REAL(r8) fslkm               ! flux lw static energy (top interface)
    REAL(r8) fqlkp               ! flux total water (bottom interface)
    REAL(r8) fqlkm               ! flux total water (top interface)
    REAL(r8) botflx              ! bottom constituent mixing ratio flux
    REAL(r8) topflx              ! top constituent mixing ratio flux
    REAL(r8) efac1               ! ratio q to convectively induced chg (btm lvl)
    REAL(r8) efac2               ! ratio q to convectively induced chg (mid lvl)
    REAL(r8) efac3               ! ratio q to convectively induced chg (top lvl)
    REAL(r8) tb(pcols,pver)      ! working storage for temp (t bar)
    REAL(r8) shb(pcols,pver)     ! working storage for spec hum (sh bar)
    REAL(r8) adjfac              ! adjustment factor (relaxation related)
    REAL(r8) rktp
    REAL(r8) rk
    INTEGER i,k              ! longitude, level indices
    INTEGER ii               ! index on "gathered" vectors
    INTEGER len1             ! vector length of "gathered" vectors
    INTEGER m                ! constituent index
    INTEGER ktp              ! tmp indx used to track top of convective layer
    !
    !---------------------------Statement functions-------------------------
    !
    !   real(r8) qhalf
    !   qhalf(sh1,sh2,shbs1,shbs2) = min(max(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
    !
    !-----------------------------------------------------------------------

    !** BAB initialize output tendencies here
    !       copy q to dq; use dq below for passive tracer transport
    cmfdt(1:ncol,:)  = 0.0_r8
    cmfdq(1:ncol,:)  = 0.0_r8
    dq(1:ncol,:,2:)  = q(1:ncol,:,2:)
    cmfmc(1:ncol,:)  = 0.0_r8
    cmfdqr(1:ncol,:) = 0.0_r8
    cmfsl(1:ncol,:)  = 0.0_r8
    cmflq(1:ncol,:)  = 0.0_r8
    qc(1:ncol,:)     = 0.0_r8
    rliq(1:ncol)     = 0.0_r8
    !
    ! Ensure that characteristic adjustment time scale (cmftau) assumed
    ! in estimate of eta isn't smaller than model time scale (ztodt)
    ! The time over which the convection is assumed to act (the adjustment
    ! time scale) can be applied with each application of the three-level
    ! cloud model, or applied to the column tendencies after a "hard"
    ! adjustment (i.e., on a 2-delta t time scale) is evaluated
    !
    IF (rlxclm) THEN
       cats   = ztodt             ! relaxation applied to column
       adjfac = ztodt/(MAX(ztodt,cmftau))
    ELSE
       cats   = MAX(ztodt,cmftau) ! relaxation applied to triplet
       adjfac = 1.0_r8
    ENDIF
    rtdt = 1.0_r8/ztodt
    !
    ! Move temperature and moisture into working storage
    !
    DO k=limcnv,pver
       DO i=1,ncol
          tb (i,k) = t(i,k)
          shb(i,k) = q(i,k,1)
       END DO
    END DO
    DO k=1,pver
       DO i=1,ncol
          icwmr(i,k) = 0.0_r8
       END DO
    END DO
    !
    ! Compute sb,hb,shbs,hbs
    !
    CALL aqsatd(tb      ,pmid    ,estemp ,shbs    ,gam     , &
         pcols   ,ncol    ,pver   ,limcnv  ,pver    )
    !
    DO k=limcnv,pver
       DO i=1,ncol
          sb (i,k) = cp*tb(i,k) + zm(i,k)*grav + phis(i)
          hb (i,k) = sb(i,k) + hlat*shb(i,k)
          hbs(i,k) = sb(i,k) + hlat*shbs(i,k)
       END DO
    END DO
    !
    ! Compute sbh, shbh
    !
    DO k=limcnv+1,pver
       DO i=1,ncol
          sbh (i,k) = 0.5_r8*(sb(i,k-1) + sb(i,k))
          shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
          hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
       END DO
    END DO
    !
    ! Specify properties at top of model (not used, but filling anyway)
    !
    DO i=1,ncol
       sbh (i,limcnv) = sb(i,limcnv)
       shbh(i,limcnv) = shb(i,limcnv)
       hbh (i,limcnv) = hb(i,limcnv)
    END DO
    !
    ! Zero vertically independent control, tendency & diagnostic arrays
    !
    DO i=1,ncol
       prec(i)  = 0.0_r8
       dzcld(i) = 0.0_r8
       cnb(i)   = 0.0_r8
       cnt(i)   = real(pver+1,r8)
       kctop1(i)   = 0
       kcbot1(i)   = 0
    END DO
    !
    ! Begin moist convective mass flux adjustment procedure.
    ! Formalism ensures that negative cloud liquid water can never occur
    !
    DO  k=pver-1,limcnv+1,-1
       DO  i=1,ncol
          etagdt(i) = 0.0_r8
          eta   (i) = 0.0_r8
          beta  (i) = 0.0_r8
          ds1   (i) = 0.0_r8
          ds2   (i) = 0.0_r8
          ds3   (i) = 0.0_r8
          dq1   (i) = 0.0_r8
          dq2   (i) = 0.0_r8
          dq3   (i) = 0.0_r8
          !
          ! Specification of "cloud base" conditions
          !
          qprime    = 0.0_r8
          tprime    = 0.0_r8
          !
          ! Assign tprime within the PBL to be proportional to the quantity
          ! tpert (which will be bounded by tpmax), passed to this routine by
          ! the PBL routine.  Don't allow perturbation to produce a dry
          ! adiabatically unstable parcel.  Assign qprime within the PBL to be
          ! an appropriately modified value of the quantity qpert (which will be
          ! bounded by shpmax) passed to this routine by the PBL routine.  The
          ! quantity qprime should be less than the local saturation value
          ! (qsattp=qsat[t+tprime,p]).  In both cases, tpert and qpert are
          ! linearly reduced toward zero as the PBL top is approached.
          !
          pblhgt = MAX(pblht(i),1.0_r8)
          IF ( (zm(i,k+1) <= pblhgt) .AND. dzcld(i) == 0.0_r8 ) THEN
             fac1   = MAX(0.0_r8,1.0_r8-zm(i,k+1)/pblhgt)
             tprime = MIN(tpert(i),tpmax)*fac1
             qsattp = shbs(i,k+1) + cp*rhlat*gam(i,k+1)*tprime
             shprme = MIN(MIN(qpert(i,1),shpmax)*fac1,MAX(qsattp-shb(i,k+1),0.0_r8))
             qprime = MAX(qprime,shprme)
          ELSE
             tprime = 0.0_r8
             qprime = 0.0_r8
          END IF
          !
          ! Specify "updraft" (in-cloud) thermodynamic properties
          !
          sc (i)    = sb (i,k+1) + cp*tprime
          shc(i)    = shb(i,k+1) + qprime
          hc (i)    = sc (i    ) + hlat*shc(i)
          vtemp4(i) = hc(i) - hbs(i,k)
          dz        = pdel(i,k)*rgas*tb(i,k)*rgrav/pmid(i,k)
          IF (vtemp4(i) > 0.0_r8) THEN
             dzcld(i) = dzcld(i) + dz
          ELSE
             dzcld(i) = 0.0_r8
          END IF
10        CONTINUE
       END DO
       !
       ! Check on moist convective instability
       ! Build index vector of points where instability exists
       !
       len1 = 0
       DO i=1,ncol
          IF (vtemp4(i) > 0.0_r8) THEN
             len1 = len1 + 1
             indx1(len1) = i
          END IF
       END DO
       IF (len1 <= 0) go to 70
       !
       ! Current level just below top level => no overshoot
       !
       IF (k <= limcnv+1) THEN
          DO ii=1,len1
             i = indx1(ii)
             temp1     = vtemp4(i)/(1.0_r8 + gam(i,k))
             cldwtr(i) = MAX(0.0_r8,(sb(i,k) - sc(i) + temp1))
             beta(i)   = 0.0_r8
             vtemp3(i) = (1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k))
          END DO
       ELSE
          !
          ! First guess at overshoot parameter using crude buoyancy closure
          ! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
          ! If pre-existing supersaturation in detrainment layer, beta=0
          ! cldwtr is temporarily equal to hlat*l (l=> liquid water)
          !
          !cdir nodep
          !DIR$ CONCURRENT
          DO ii=1,len1
             i = indx1(ii)
             temp1     = vtemp4(i)/(1.0_r8 + gam(i,k))
             cldwtr(i) = MAX(0.0_r8,(sb(i,k)-sc(i)+temp1))
             betamx(i) = 1.0_r8 - c0*MAX(0.0_r8,(dzcld(i)-dzmin))
             b1        = (hc(i) - hbs(i,k-1))*pdel(i,k-1)
             b2        = (hc(i) - hbs(i,k  ))*pdel(i,k  )
             beta(i)   = MAX(betamn,MIN(betamx(i), 1.0_r8 + b1/b2))
             IF (hbs(i,k-1) <= hb(i,k-1)) beta(i) = 0.0_r8
             !
             ! Bound maximum beta to ensure physically realistic solutions
             !
             ! First check constrains beta so that eta remains positive
             ! (assuming that eta is already positive for beta equal zero)
             !
             vtemp1(i) = -(hbh(i,k+1) - hc(i))*pdel(i,k)*rpdel(i,k+1)+ &
                  (1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i))
             vtemp2(i) = (1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k))
             vtemp3(i) = vtemp2(i)
             IF ((beta(i)*vtemp2(i) - vtemp1(i)) > 0.0_r8) THEN
                betamx(i) = 0.99_r8*(vtemp1(i)/vtemp2(i))
                beta(i)   = MAX(0.0_r8,MIN(betamx(i),beta(i)))
             END IF
          END DO
          !
          ! Second check involves supersaturation of "detrainment layer"
          ! small amount of supersaturation acceptable (by ssfac factor)
          !
          !cdir nodep
          !DIR$ CONCURRENT
          DO ii=1,len1
             i = indx1(ii)
             IF (hb(i,k-1) < hbs(i,k-1)) THEN
                vtemp1(i) = vtemp1(i)*rpdel(i,k)
                temp2 = gam(i,k-1)*(sbh(i,k) - sc(i) + cldwtr(i)) -  &
                     hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
                temp3 = vtemp3(i)*rpdel(i,k)
                vtemp2(i) = (ztodt/cats)*(hc(i) - hbs(i,k))*temp2/ &
                     (pdel(i,k-1)*(hbs(i,k-1) - hb(i,k-1))) + temp3
                IF ((beta(i)*vtemp2(i) - vtemp1(i)) > 0.0_r8) THEN
                   betamx(i) = ssfac*(vtemp1(i)/vtemp2(i))
                   beta(i)   = MAX(0.0_r8,MIN(betamx(i),beta(i)))
                END IF
             ELSE
                beta(i) = 0.0_r8
             END IF
          END DO
          !
          ! Third check to avoid introducing 2 delta x thermodynamic
          ! noise in the vertical ... constrain adjusted h (or theta e)
          ! so that the adjustment doesn't contribute to "kinks" in h
          !
          !cdir nodep
          !DIR$ CONCURRENT
          DO ii=1,len1
             i = indx1(ii)
             g = MIN(0.0_r8,hb(i,k) - hb(i,k-1))
             temp1 = (hb(i,k) - hb(i,k-1) - g)*(cats/ztodt)/(hc(i) - hbs(i,k))
             vtemp1(i) = temp1*vtemp1(i) + (hc(i) - hbh(i,k+1))*rpdel(i,k)
             vtemp2(i) = temp1*vtemp3(i)*rpdel(i,k) + (hc(i) - hbh(i,k) - cldwtr(i))* &
                  (rpdel(i,k) + rpdel(i,k+1))
             IF ((beta(i)*vtemp2(i) - vtemp1(i)) > 0.0_r8) THEN
                IF (vtemp2(i) /= 0.0_r8) THEN
                   betamx(i) = vtemp1(i)/vtemp2(i)
                ELSE
                   betamx(i) = 0.0_r8
                END IF
                beta(i) = MAX(0.0_r8,MIN(betamx(i),beta(i)))
             END IF
          END DO
       END IF
       !
       ! Calculate mass flux required for stabilization.
       !
       ! Ensure that the convective mass flux, eta, is positive by
       ! setting negative values of eta to zero..
       ! Ensure that estimated mass flux cannot move more than the
       ! minimum of total mass contained in either layer k or layer k+1.
       ! Also test for other pathological cases that result in non-
       ! physical states and adjust eta accordingly.
       !
       !cdir nodep
       !DIR$ CONCURRENT
       DO ii=1,len1
          i = indx1(ii)
          beta(i) = MAX(0.0_r8,beta(i))
          temp1 = hc(i) - hbs(i,k)
          temp2 = ((1.0_r8 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i)) - &
               beta(i)*vtemp3(i))*rpdel(i,k) - (hbh(i,k+1) - hc(i))*rpdel(i,k+1)
          eta(i) = temp1/(temp2*grav*cats)
          tmass = MIN(pdel(i,k),pdel(i,k+1))*rgrav
          IF (eta(i) > tmass*rtdt .OR. eta(i) <= 0.0_r8) eta(i) = 0.0_r8
          !
          ! Check on negative q in top layer (bound beta)
          !
          IF (shc(i)-shbh(i,k) < 0.0_r8 .AND. beta(i)*eta(i) /= 0.0_r8) THEN
             denom = eta(i)*grav*ztodt*(shc(i) - shbh(i,k))*rpdel(i,k-1)
             beta(i) = MAX(0.0_r8,MIN(-0.999_r8*shb(i,k-1)/denom,beta(i)))
          END IF
          !
          ! Check on negative q in middle layer (zero eta)
          !
          qtest1 = shb(i,k) + eta(i)*grav*ztodt*((shc(i) - shbh(i,k+1)) - &
               (1.0_r8 - beta(i))*cldwtr(i)*rhlat - beta(i)*(shc(i) - shbh(i,k)))* &
               rpdel(i,k)
          IF (qtest1 <= 0.0_r8) eta(i) = 0.0_r8
          !
          ! Check on negative q in lower layer (bound eta)
          !
          fac1 = -(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
          qtest2 = shb(i,k+1) - eta(i)*grav*ztodt*fac1
          IF (qtest2 < 0.0_r8) THEN
             eta(i) = 0.99_r8*shb(i,k+1)/(grav*ztodt*fac1)
          END IF
          etagdt(i) = eta(i)*grav*ztodt
       END DO
       !
       ! Calculate cloud water, rain water, and thermodynamic changes
       !
       !cdir nodep
       !DIR$ CONCURRENT
       DO  ii=1,len1
          i = indx1(ii)
          icwmr(i,k) = cldwtr(i)*rhlat
          cldwtr(i) = etagdt(i)*cldwtr(i)*rhlat*rgrav
          ! JJH changes to facilitate export of cloud liquid water --------------------------------
          totcond(i) = (1.0_r8 - beta(i))*cldwtr(i)
          rnwtr(i) = MIN(totcond(i),c0*(dzcld(i)-dzmin)*cldwtr(i))
          ds1(i) = etagdt(i)*(sbh(i,k+1) - sc(i))*rpdel(i,k+1)
          dq1(i) = etagdt(i)*(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
          ds2(i) = (etagdt(i)*(sc(i) - sbh(i,k+1)) +  &
               hlat*grav*cldwtr(i) - beta(i)*etagdt(i)*(sc(i) - sbh(i,k)))*rpdel(i,k)
          ! JJH change for export of cloud liquid water; must use total condensate 
          ! since rainwater no longer represents total condensate
          dq2(i) = (etagdt(i)*(shc(i) - shbh(i,k+1)) - grav*totcond(i) - beta(i)* &
               etagdt(i)*(shc(i) - shbh(i,k)))*rpdel(i,k)
          ds3(i) = beta(i)*(etagdt(i)*(sc(i) - sbh(i,k)) - hlat*grav*cldwtr(i))* &
               rpdel(i,k-1)
          dq3(i) = beta(i)*etagdt(i)*(shc(i) - shbh(i,k))*rpdel(i,k-1)
          !
          ! Isolate convective fluxes for later diagnostics
          !
          fslkp = eta(i)*(sc(i) - sbh(i,k+1))
          fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) - hlat*cldwtr(i)*rtdt)
          fqlkp = eta(i)*(shc(i) - shbh(i,k+1))
          fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))
          !
          ! Update thermodynamic profile (update sb, hb, & hbs later)
          !
          tb (i,k+1) = tb(i,k+1)  + ds1(i)*rcp
          tb (i,k  ) = tb(i,k  )  + ds2(i)*rcp
          tb (i,k-1) = tb(i,k-1)  + ds3(i)*rcp
          shb(i,k+1) = shb(i,k+1) + dq1(i)
          shb(i,k  ) = shb(i,k  ) + dq2(i)
          shb(i,k-1) = shb(i,k-1) + dq3(i)
          !
          ! ** Update diagnostic information for final budget **
          ! Tracking precipitation, temperature & specific humidity tendencies,
          ! rainout term, convective mass flux, convective liquid
          ! water static energy flux, and convective total water flux
          ! The variable afac makes the necessary adjustment to the
          ! diagnostic fluxes to account for adjustment time scale based on
          ! how relaxation time scale is to be applied (column vs. triplet)
          !
          prec(i)    = prec(i) + (rnwtr(i)/rhoh2o)*adjfac
          !
          ! The following variables have units of "units"/second
          !
          cmfdt (i,k+1) = cmfdt (i,k+1) + ds1(i)*rtdt*adjfac
          cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)*rtdt*adjfac
          cmfdt (i,k-1) = cmfdt (i,k-1) + ds3(i)*rtdt*adjfac
          cmfdq (i,k+1) = cmfdq (i,k+1) + dq1(i)*rtdt*adjfac
          cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rtdt*adjfac
          cmfdq (i,k-1) = cmfdq (i,k-1) + dq3(i)*rtdt*adjfac
          ! JJH changes to export cloud liquid water --------------------------------
          qc    (i,k  ) = (grav*(totcond(i)-rnwtr(i))*rpdel(i,k))*rtdt*adjfac
          cmfdqr(i,k  ) = cmfdqr(i,k  ) + (grav*rnwtr(i)*rpdel(i,k))*rtdt*adjfac
          cmfmc (i,k+1) = cmfmc (i,k+1) + eta(i)*adjfac
          cmfmc (i,k  ) = cmfmc (i,k  ) + beta(i)*eta(i)*adjfac
          !
          ! The following variables have units of w/m**2
          !
          cmfsl (i,k+1) = cmfsl (i,k+1) + fslkp*adjfac
          cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm*adjfac
          cmflq (i,k+1) = cmflq (i,k+1) + hlat*fqlkp*adjfac
          cmflq (i,k  ) = cmflq (i,k  ) + hlat*fqlkm*adjfac
30        CONTINUE
       END DO! end of ii=1,len1
       !
       ! Next, convectively modify passive constituents
       ! For now, when applying relaxation time scale to thermal fields after
       ! entire column has undergone convective overturning, constituents will
       ! be mixed using a "relaxed" value of the mass flux determined above
       ! Although this will be inconsistant with the treatment of the thermal
       ! fields, it's computationally much cheaper, no more-or-less justifiable,
       ! and consistent with how the history tape mass fluxes would be used in
       ! an off-line mode (i.e., using an off-line transport model)
       !
       DO  m=2,pcnst+pnats    ! note: indexing assumes water is first field
          IF (cnst_get_type_byind(m,pcnst+pnats).EQ.'dry') THEN
             pd(:ncol,:) = pdeldry(:ncol,:)
             rpd(:ncol,:) = rpdeldry(:ncol,:)
             pm(:ncol,:) = pmiddry(:ncol,:)
          ELSE
             pd(:ncol,:) = pdel(:ncol,:)
             rpd(:ncol,:) = rpdel(:ncol,:)
             pm(:ncol,:) = pmid(:ncol,:)
          ENDIF
          !cdir nodep
          !DIR$ CONCURRENT
          DO  ii=1,len1
             i = indx1(ii)
             !
             ! If any of the reported values of the constituent is negative in
             ! the three adjacent levels, nothing will be done to the profile
             !
             IF ((dq(i,k+1,m) < 0.0_r8) .OR. (dq(i,k,m) < 0.0_r8) .OR. (dq(i,k-1,m) < 0.0_r8)) go to 40
             !
             ! Specify constituent interface values (linear interpolation)
             !
             cmrh(i,k  ) = 0.5_r8*(dq(i,k-1,m) + dq(i,k  ,m))
             cmrh(i,k+1) = 0.5_r8*(dq(i,k  ,m) + dq(i,k+1,m))
             !
             ! Specify perturbation properties of constituents in PBL
             !
             pblhgt = MAX(pblht(i),1.0_r8)
             IF ( (zm(i,k+1) <= pblhgt) .AND. dzcld(i) == 0.0_r8 ) THEN
                fac1 = MAX(0.0_r8,1.0_r8-zm(i,k+1)/pblhgt)
                cmrc(i) = dq(i,k+1,m) + qpert(i,m)*fac1
             ELSE
                cmrc(i) = dq(i,k+1,m)
             END IF
             !
             ! Determine fluxes, flux divergence => changes due to convection
             ! Logic must be included to avoid producing negative values. A bit
             ! messy since there are no a priori assumptions about profiles.
             ! Tendency is modified (reduced) when pending disaster detected.
             !
             botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
             topflx   = beta(i)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
             dcmr1(i) = -botflx*rpd(i,k+1)
             efac1    = 1.0_r8
             efac2    = 1.0_r8
             efac3    = 1.0_r8
             !
             if (dq(i,k+1,m)+dcmr1(i) < 0.0_r8) then
                     if ( abs(dcmr1(i)) > 1.e-300_r8 ) then
                        efac1 = max(tiny,abs(dq(i,k+1,m)/dcmr1(i)) - eps)
                     else
                        efac1 = tiny
                     endif
             end if
             !
             IF (efac1 == tiny .OR. efac1 > 1.0_r8) efac1 = 0.0_r8
             dcmr1(i) = -efac1*botflx*rpd(i,k+1)
             dcmr2(i) = (efac1*botflx - topflx)*rpd(i,k)
             !
             if (dq(i,k,m)+dcmr2(i) < 0.0_r8) then
                     if ( abs(dcmr2(i)) > 1.e-300_r8 ) then
                        efac2 = max(tiny,abs(dq(i,k  ,m)/dcmr2(i)) - eps)
                     else
                        efac2 = tiny
                     endif
             end if
             !
             IF (efac2 == tiny .OR. efac2 > 1.0_r8) efac2 = 0.0_r8
             dcmr2(i) = (efac1*botflx - efac2*topflx)*rpd(i,k)
             dcmr3(i) = efac2*topflx*rpdel(i,k-1)
             !
             if ( (dq(i,k-1,m)+dcmr3(i) < 0.0_r8 ) ) then
                     if  ( abs(dcmr3(i)) > 1.e-300_r8 ) then
                        efac3 = max(tiny,abs(dq(i,k-1,m)/dcmr3(i)) - eps)
                     else
                        efac3 = tiny
                     endif
             end if
             !
             IF (efac3 == tiny .OR. efac3 > 1.0_r8) efac3 = 0.0_r8
             efac3    = MIN(efac2,efac3)
             dcmr2(i) = (efac1*botflx - efac3*topflx)*rpd(i,k)
             dcmr3(i) = efac3*topflx*rpd(i,k-1)
             !
             dq(i,k+1,m) = dq(i,k+1,m) + dcmr1(i)
             dq(i,k  ,m) = dq(i,k  ,m) + dcmr2(i)
             dq(i,k-1,m) = dq(i,k-1,m) + dcmr3(i)
40           CONTINUE
          END DO !  end of ii=1,len1
50        CONTINUE                ! end of m=2,pcnst+pnats loop
       END DO !end of m=2,pcnst+pnats loop
       !
       ! Constituent modifications complete
       !
       IF (k == limcnv+1) go to 60
       !
       ! Complete update of thermodynamic structure at integer levels
       ! gather/scatter points that need new values of shbs and gamma
       !
       DO ii=1,len1
          i = indx1(ii)
          vtemp1(ii     ) = tb(i,k)
          vtemp1(ii+len1) = tb(i,k-1)
          vtemp2(ii     ) = pmid(i,k)
          vtemp2(ii+len1) = pmid(i,k-1)
       END DO
       CALL vqsatd (vtemp1  ,vtemp2  ,estemp  ,vtemp3  , vtemp4  , &
            2*len1   )    ! using estemp as extra long vector
       !cdir nodep
       !DIR$ CONCURRENT
       DO ii=1,len1
          i = indx1(ii)
          shbs(i,k  ) = vtemp3(ii     )
          shbs(i,k-1) = vtemp3(ii+len1)
          gam(i,k  ) = vtemp4(ii     )
          gam(i,k-1) = vtemp4(ii+len1)
          sb (i,k  ) = sb(i,k  ) + ds2(i)
          sb (i,k-1) = sb(i,k-1) + ds3(i)
          hb (i,k  ) = sb(i,k  ) + hlat*shb(i,k  )
          hb (i,k-1) = sb(i,k-1) + hlat*shb(i,k-1)
          hbs(i,k  ) = sb(i,k  ) + hlat*shbs(i,k  )
          hbs(i,k-1) = sb(i,k-1) + hlat*shbs(i,k-1)
       END DO
       !
       ! Update thermodynamic information at half (i.e., interface) levels
       !
       !DIR$ CONCURRENT
       DO ii=1,len1
          i = indx1(ii)
          sbh (i,k) = 0.5_r8*(sb(i,k) + sb(i,k-1))
          shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
          hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
          sbh (i,k-1) = 0.5_r8*(sb(i,k-1) + sb(i,k-2))
          shbh(i,k-1) = qhalf(shb(i,k-2),shb(i,k-1),shbs(i,k-2),shbs(i,k-1))
          hbh (i,k-1) = sbh(i,k-1) + hlat*shbh(i,k-1)
       END DO
       !
       ! Ensure that dzcld is reset if convective mass flux zero
       ! specify the current vertical extent of the convective activity
       ! top of convective layer determined by size of overshoot param.
       !
60     DO i=1,ncol
          etagt0 = eta(i).GT.0.0_r8
          IF ( .NOT. etagt0) dzcld(i) = 0.0_r8
          IF (etagt0 .AND. beta(i) > betamn) THEN
             ktp = k-1
          ELSE
             ktp = k
          END IF
          IF (etagt0) THEN
             rk=k
             rktp=ktp
             cnt(i) = MIN(cnt(i),rktp)
             cnb(i) = MAX(cnb(i),rk)
          END IF
       END DO
70     CONTINUE                  ! end of k loop
    END DO ! end of k loop
    !
    ! ** apply final thermodynamic tendencies **
    !
    !**BAB don't update input profiles
!!$                  do k=limcnv,pver
!!$                     do i=1,ncol
!!$                        t (i,k) = t (i,k) + cmfdt(i,k)*ztodt
!!$                        q(i,k,1) = q(i,k,1) + cmfdq(i,k)*ztodt
!!$                     end do
!!$                  end do
    ! Set output q tendencies 
    dq(1:ncol,:,1 ) = cmfdq(1:ncol,:)
    dq(1:ncol,:,2:) = (dq(1:ncol,:,2:) - q(1:ncol,:,2:))/ztodt
    !
    ! Kludge to prevent cnb-cnt from being zero (in the event
    ! someone decides that they want to divide by this quantity)
    !
    DO i=1,ncol
       IF (cnb(i) /= 0.0_r8 .AND. cnb(i) == cnt(i)) THEN
          cnt(i) = cnt(i) - 1.0_r8
       END IF
    END DO
    !
    DO i=1,ncol
       precc(i) = prec(i)*rtdt
       kctop1(i) = pver-cnt(i)
       kcbot1(i) = pver-cnb(i)
       IF(precc(i)<=0.0_r8)THEN
          noshal(i) = 1
       END IF
    END DO
    !
    ! Compute reserved liquid (not yet in cldliq) for energy integrals.
    ! Treat rliq as flux out bottom, to be added back later.
    DO k = 1, pver
       DO i = 1, ncol
          IF( noshal(i) == 0)THEN
             rliq(i) = rliq(i) + qc(i,k)*pdel(i,k)/grav
          END IF
       END DO
    END DO
    rliq(1:ncol) = rliq(1:ncol) /1000.0_r8

    RETURN                 ! we're all done ... return to calling procedure
  END SUBROUTINE cmfmca
  !=========================================================================================
  !===============================================================================
  !===============================================================================
  SUBROUTINE physics_update(state, tend, ptend, dt,ppcnst,pcols,pver,latco)
    !-----------------------------------------------------------------------
    ! Update the state and or tendency structure with the parameterization tendencies
    !-----------------------------------------------------------------------
    !    use geopotential, only: geopotential_dse
    !    use physconst,    only: cpair, gravit, rair, zvir
    !    use constituents, only: cnst_get_ind
    !#if ( defined SCAM )
    !#include <max.h>
    !    use scamMod, only: switch
    !#endif
    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN) :: ppcnst,latco
    INTEGER, INTENT(IN) :: pcols
    INTEGER, INTENT(IN) :: pver
    TYPE(physics_ptend), INTENT(inout)  :: ptend   ! Parameterization tendencies

    TYPE(physics_state), INTENT(inout)  :: state   ! Physics state variables
    TYPE(physics_tend ), INTENT(inout)  :: tend    ! Physics tendencies

    REAL(r8), INTENT(in) :: dt                     ! time step
    REAL(r8) :: qq  (pcols,pver,ppcnst)                     ! time step
    REAL(r8) :: qq2 (pcols,pver)                     ! time step

    !
    !---------------------------Local storage-------------------------------
    INTEGER :: i,k,m                               ! column,level,constituent indices
    !INTEGER :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    INTEGER :: ncol                                ! number of columns
    CHARACTER*40 :: name    ! param and tracer name for qneg3
    !-----------------------------------------------------------------------
    !#if ( defined SCAM )
    !    ! The column radiation model does not update the state
    !    if(switch(CRM_SW+1)) return
    !#endif
    ncol = state%ncol(latco)

    ! Update u,v fields
    IF(ptend%lu(latco)) THEN
       DO k = ptend%top_level(latco), ptend%bot_level(latco)
          DO i = 1, ncol
             state%u  (i,k,latco) = state%u  (i,k,latco) + ptend%u(i,k,latco) * dt
             tend%dudt(i,k,latco) = tend%dudt(i,k,latco) + ptend%u(i,k,latco)
          END DO
       END DO
    END IF

    IF(ptend%lv(latco)) THEN
       DO k = ptend%top_level(latco), ptend%bot_level(latco)
          DO i = 1, ncol
             state%v  (i,k,latco) = state%v  (i,k,latco) + ptend%v(i,k,latco) * dt
             tend%dvdt(i,k,latco) = tend%dvdt(i,k,latco) + ptend%v(i,k,latco)
          END DO
       END DO
    END IF

    ! Update dry static energy
    IF(ptend%ls(latco)) THEN
       DO k = ptend%top_level(latco), ptend%bot_level(latco)
          DO i = 1, ncol
             state%s(i,k,latco)   = state%s(i,k,latco)   + ptend%s(i,k,latco) * dt
             tend%dtdt(i,k,latco) = tend%dtdt(i,k,latco) + ptend%s(i,k,latco)/cpair
          END DO
       END DO
    END IF

    ! Update constituents, all schemes use time split q: no tendency kept
    !call cnst_get_ind('CLDICE', ixcldice)
    !call cnst_get_ind('CLDLIQ', ixcldliq)
    ! ixcldliq=2
    !ixcldice=3
    DO m = 1, ppcnst
       IF(ptend%lq(latco,m)) THEN
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             DO i = 1,ncol
                !PRINT*,state%q(i,k,latco,m) , ptend%q(i,k,latco,m) , dt
                state%q(i,k,latco,m) = state%q(i,k,latco,m) + ptend%q(i,k,latco,m) * dt
             END DO
          END DO
          ! now test for mixing ratios which are too small
          name = TRIM(ptend%name(latco)) !// '/' // trim(cnst_name(m))
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             DO i = 1,ncol
                qq(i,k,m)=state%q(i,k,latco,m)
             END DO
          END DO
          CALL qneg3(TRIM(name),  ncol, pcols, pver, m, m, qmin(m), qq(1,1,m))
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             DO i = 1,ncol
                state%q(i,k,latco,m)=qq(i,k,m)
             END DO
          END DO
       END IF
    END DO

    ! special test for cloud water
    IF(ptend%lq(latco,ixcldliq)) THEN
       IF (ptend%name(latco) == 'stratiform' .OR. ptend%name(latco) == 'cldwat'  ) THEN
          IF(PERGRO)THEN
             WHERE (state%q(1:ncol,1:pver,latco,ixcldliq) < 1.0e-12_r8)
                state%q(1:ncol,1:pver,latco,ixcldliq) = 0.0_r8
             END WHERE
          ENDIF
       ELSE IF (ptend%name(latco) == 'convect_deep' .OR. ptend%name(latco) == 'convect_shallow' &
       .OR. ptend%name(latco) == 'cmfmca') THEN
          WHERE (state%q(1:ncol,1:pver,latco,ixcldliq) < 1.0e-36_r8)
             state%q(1:ncol,1:pver,latco,ixcldliq) = 0.0_r8
          END WHERE
       END IF
    END IF
    IF(ptend%lq(latco,ixcldice)) THEN
       IF (ptend%name(latco) == 'stratiform' .OR. ptend%name(latco) == 'cldwat'  ) THEN
          IF( PERGRO)THEN
             WHERE (state%q(1:ncol,1:pver,latco,ixcldice) < 1.0e-12_r8)
                state%q(1:ncol,1:pver,latco,ixcldice) = 0.0_r8
             END WHERE
          ENDIF
       ELSE IF (ptend%name(latco) == 'convect_deep' .OR. ptend%name(latco) == 'convect_shallow' &
       .OR. ptend%name(latco) == 'cmfmca') THEN
          WHERE (state%q(1:ncol,1:pver,latco,ixcldice) < 1.0e-36_r8)
             state%q(1:ncol,1:pver,latco,ixcldice) = 0.0_r8
          END WHERE
       END IF
    END IF

    ! Derive new temperature and geopotential fields if heating or water tendency not 0.
    IF (ptend%ls(latco) .OR. ptend%lq(latco,1)) THEN
       CALL geopotential_dse(                                                                    &
            state%lnpint(1:pcols,1:pver+1,latco), state%lnpmid(1:pcols,1:pver,latco) , state%pint (1:pcols,1:pver+1,latco)  , &
            state%pmid  (1:pcols,1:pver,latco), state%pdel  (1:pcols,1:pver,latco) , state%rpdel(1:pcols,1:pver,latco)  , &
            state%s     (1:pcols,1:pver,latco), state%q     (1:ncol,1:pver,latco,1), state%phis (1:pcols,latco) , rair  , &
            gravit      , cpair        ,zvir        , &
            state%t(1:pcols,1:pver,latco)     , state%zi(1:pcols,1:pver+1,latco)    , state%zm(1:pcols,1:pver,latco)       ,&
            ncol      ,pcols, pver, pver+1          )

    END IF

    ! Reset all parameterization tendency flags to false
    CALL physics_ptend_reset(ptend,latco,ppcnst,pver)

  END SUBROUTINE physics_update




  !===============================================================================
  SUBROUTINE physics_ptend_sum(ptend, ptend_sum, state,ppcnst,latco)
    !-----------------------------------------------------------------------
    ! Add ptend fields to ptend_sum for ptend logical flags = .true.
    ! Where ptend logical flags = .false, don't change ptend_sum
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN   ) ::ppcnst,latco
    TYPE(physics_ptend), INTENT(in)     :: ptend   ! New parameterization tendencies
    TYPE(physics_ptend), INTENT(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    TYPE(physics_state), INTENT(in)     :: state   ! New parameterization tendencies

    !---------------------------Local storage-------------------------------
    INTEGER :: i,k,m                               ! column,level,constituent indices
    INTEGER :: ncol                                ! number of columns

    !-----------------------------------------------------------------------
    ncol = state%ncol(latco)


    ! Update u,v fields
    IF(ptend%lu(latco)) THEN
       ptend_sum%lu(latco) = .TRUE.
       DO i = 1, ncol
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             ptend_sum%u(i,k,latco) = ptend_sum%u(i,k,latco) + ptend%u(i,k,latco)
          END DO
          ptend_sum%taux_srf(i,latco) = ptend_sum%taux_srf(i,latco) + ptend%taux_srf(i,latco)
          ptend_sum%taux_top(i,latco) = ptend_sum%taux_top(i,latco) + ptend%taux_top(i,latco)
       END DO
    END IF

    IF(ptend%lv(latco)) THEN
       ptend_sum%lv(latco) = .TRUE.
       DO i = 1, ncol
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             ptend_sum%v(i,k,latco) = ptend_sum%v(i,k,latco) + ptend%v(i,k,latco)
          END DO
          ptend_sum%tauy_srf(i,latco) = ptend_sum%tauy_srf(i,latco) + ptend%tauy_srf(i,latco)
          ptend_sum%tauy_top(i,latco) = ptend_sum%tauy_top(i,latco) + ptend%tauy_top(i,latco)
       END DO
    END IF


    IF(ptend%ls(latco)) THEN
       ptend_sum%ls(latco) = .TRUE.
       DO i = 1, ncol
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             ptend_sum%s(i,k,latco) = ptend_sum%s(i,k,latco) + ptend%s(i,k,latco)
          END DO
          ptend_sum%hflux_srf(i,latco) = ptend_sum%hflux_srf(i,latco) + ptend%hflux_srf(i,latco)
          ptend_sum%hflux_top(i,latco) = ptend_sum%hflux_top(i,latco) + ptend%hflux_top(i,latco)
       END DO
    END IF

    ! Update constituents
    DO m = 1, ppcnst
       IF(ptend%lq(latco,m)) THEN
          ptend_sum%lq(latco,m) = .TRUE.
          DO i = 1,ncol
             DO k = ptend%top_level(latco), ptend%bot_level(latco)
                ptend_sum%q(i,k,latco,m) = ptend_sum%q(i,k,latco,m) + ptend%q(i,k,latco,m)
             END DO
             ptend_sum%cflx_srf(i,latco,m) = ptend_sum%cflx_srf(i,latco,m) + ptend%cflx_srf(i,latco,m)
             ptend_sum%cflx_top(i,latco,m) = ptend_sum%cflx_top(i,latco,m) + ptend%cflx_top(i,latco,m)
          END DO
       END IF
    END DO


  END SUBROUTINE physics_ptend_sum

  !===============================================================================
  SUBROUTINE physics_ptend_init(ptend,latco,ppcnst,pver)
    !-----------------------------------------------------------------------
    ! Initialize the parameterization tendency structure to "empty"
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN   ) :: ppcnst,latco
    INTEGER, INTENT(IN   ) :: pver
    TYPE(physics_ptend), INTENT(inout)  :: ptend   ! Parameterization tendencies
    !-----------------------------------------------------------------------
    ptend%name(latco)  = "none"
    ptend%lq(latco,:) = .TRUE.
    ptend%ls(latco)    = .TRUE.
    ptend%lu(latco)    = .TRUE.
    ptend%lv(latco)    = .TRUE.

    CALL physics_ptend_reset(ptend,latco,ppcnst,pver)

    RETURN
  END SUBROUTINE physics_ptend_init
  !===============================================================================
  SUBROUTINE physics_ptend_reset(ptend,latco,ppcnst,pver)
    !-----------------------------------------------------------------------
    ! Reset the parameterization tendency structure to "empty"
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN   ) :: ppcnst,latco
    INTEGER, INTENT(IN   ) :: pver
    TYPE(physics_ptend), INTENT(inout)  :: ptend   ! Parameterization tendencies
    !-----------------------------------------------------------------------
    INTEGER :: m             ! Index for constiuent
    !-----------------------------------------------------------------------

    IF(ptend%ls(latco)) THEN
       ptend%s(:,:,latco) = 0.0_r8
       ptend%hflux_srf(:,latco) = 0.0_r8
       ptend%hflux_top(:,latco) = 0.0_r8
    ENDIF
    IF(ptend%lu(latco)) THEN
       ptend%u(:,:,latco) = 0.0_r8
       ptend%taux_srf(:,latco) = 0.0_r8
       ptend%taux_top(:,latco) = 0.0_r8
    ENDIF
    IF(ptend%lv(latco)) THEN
       ptend%v(:,:,latco) = 0.0_r8
       ptend%tauy_srf(:,latco) = 0.0_r8
       ptend%tauy_top(:,latco) = 0.0_r8
    ENDIF
    DO m = 1, ppcnst
       IF(ptend%lq(latco,m)) THEN
          ptend%q(:,:,latco,m) = 0.0_r8
          ptend%cflx_srf(:,latco,m) = 0.0_r8
          ptend%cflx_top(:,latco,m) = 0.0_r8
       ENDIF
    END DO

    ptend%name(latco)  = "none"
    ptend%lq(latco,:) = .FALSE.
    ptend%ls(latco)    = .FALSE.
    ptend%lu(latco)    = .FALSE.
    ptend%lv(latco)    = .FALSE.

    ptend%top_level(latco) = 1
    ptend%bot_level(latco) = pver

    RETURN
  END SUBROUTINE physics_ptend_reset


  !===============================================================================
  SUBROUTINE physics_state_copy(pcols,latco,&
       pver, pverp,ppcnst, cnst_need_pdeldry)

    !use ppgrid,       only: pver, pverp
    !use constituents, only: ppcnst, cnst_need_pdeldry

    IMPLICIT NONE

    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: pver, pverp,ppcnst,latco,pcols
    LOGICAL, INTENT(IN   ) ::  cnst_need_pdeldry
    !TYPE(physics_state), INTENT(in) :: state
    !TYPE(physics_state), INTENT(out) :: state1

    !
    ! Local variables
    !
    INTEGER i, k, m


    state1%ncol(latco)  = pcols
    state1%count(latco) = state%count (latco)

    DO i = 1, pcols
       state1%lat   (i,latco)    = state%lat(i,latco)
       state1%lon   (i,latco)    = state%lon(i,latco)
       state1%ps    (i,latco)     = state%ps(i,latco)
       state1%phis  (i,latco)   = state%phis(i,latco)
       state1%te_ini(i,latco) = state%te_ini(i,latco) 
       state1%te_cur(i,latco) = state%te_cur(i,latco) 
       state1%tw_ini(i,latco) = state%tw_ini(i,latco) 
       state1%tw_cur(i,latco) = state%tw_cur(i,latco) 
    END DO

    DO k = 1, pver
       DO i = 1, pcols
          state1%t(i,k,latco)         = state%t(i,k,latco) 
          state1%u(i,k,latco)         = state%u(i,k,latco) 
          state1%v(i,k,latco)         = state%v(i,k,latco) 
          state1%s(i,k,latco)         = state%s(i,k,latco) 
          state1%omega(i,k,latco)     = state%omega(i,k,latco) 
          state1%pmid(i,k,latco)      = state%pmid(i,k,latco) 
          state1%pdel(i,k,latco)      = state%pdel(i,k,latco) 
          state1%rpdel(i,k,latco)     = state%rpdel(i,k,latco) 
          state1%lnpmid(i,k,latco)    = state%lnpmid(i,k,latco) 
          !          state1%exner(i,k,latco)     = state%exner(i,k,latco) 
          state1%zm(i,k,latco)        = state%zm(i,k,latco)
       END DO
    END DO

    DO k = 1, pverp
       DO i = 1, pcols
          state1%pint(i,k,latco)      = state%pint(i,k,latco) 
          state1%lnpint(i,k,latco)    = state%lnpint(i,k,latco) 
          state1%zi(i,k,latco)        = state% zi(i,k,latco) 
       END DO
    END DO


    !    IF ( cnst_need_pdeldry ) THEN
    !       DO i = 1, pcols
    !          state1%psdry(i,latco)  = state%psdry(i,latco) 
    !       END DO
    !       DO k = 1, pver
    !          DO i = 1, pcols
    !             state1%lnpmiddry(i,k,latco)  = state%lnpmiddry(i,k,latco) 
    !             state1%pmiddry  (i,k,latco)  = state%pmiddry  (i,k,latco) 
    !             state1%pdeldry  (i,k,latco)  = state%pdeldry  (i,k,latco) 
    !             state1%rpdeldry (i,k,latco)  = state%rpdeldry (i,k,latco)
    !          END DO
    !       END DO
    !       DO k = 1, pverp
    !          DO i = 1, pcols
    !             state1%pintdry  (i,k,latco)   = state%pintdry(i,k,latco)
    !             state1%lnpintdry(i,k,latco) = state%lnpintdry(i,k,latco) 
    !          END DO
    !       END DO
    !    ENDIF !cnst_need_pdeldry

    DO m = 1, ppcnst
       DO k = 1, pver
          DO i = 1, pcols
             state1%q(i,k,latco,m) = state%q(i,k,latco,m) 
          END DO
       END DO
    END DO

  END  SUBROUTINE physics_state_copy
  !-------------module physics_types
  !===============================================================================

  SUBROUTINE physics_tend_init(tend,latco)


    IMPLICIT NONE

    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: latco
    TYPE(physics_tend), INTENT(inout) :: tend

    !
    ! Local variables
    !

    tend%dtdt   (:,:,latco)     = 0.0_r8
    tend%dudt   (:,:,latco)    = 0.0_r8
    tend%dvdt   (:,:,latco)    = 0.0_r8
    tend%flx_net(:,latco)   = 0.0_r8
    tend%te_tnd (:,latco)   = 0.0_r8
    tend%tw_tnd (:,latco)   = 0.0_r8

  END SUBROUTINE physics_tend_init


  !===============================================================================
  SUBROUTINE zm_conv_evap( &
       ncol        , &!INTENT(in) :: ncol               ! number of columns and chunk index
       pcols       , &!INTENT(in) :: pcols                    ! number of columns (max)
       pver        , &!INTENT(in) :: pver                    ! number of vertical levels
       pverp       , &!INTENT(in) :: pverp                    ! number of vertical levels + 1
       t           , &!INTENT(in), DIMENSION(pcols,pver) :: t          ! temperature (K)
       pmid        , &!INTENT(in), DIMENSION(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
       pdel        , &!INTENT(in), DIMENSION(pcols,pver) :: pdel       ! layer thickness (Pa)
       q           , &!INTENT(in), DIMENSION(pcols,pver) :: q          ! water vapor (kg/kg)
       tend_s      , &!INTENT(inout), DIMENSION(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
       tend_q      , &!INTENT(inout), DIMENSION(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
       prdprec     , &!INTENT(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
       cldfrc      , &!INTENT(in   ) :: cldfrc(pcols,pver) ! cloud fraction
       deltat      , &!INTENT(in   ) :: deltat             ! time step
       prec        , &!INTENT(inout) :: prec(pcols)        ! Convective-scale preciptn rate
       snow        , &!INTENT(out)   :: snow(pcols)        ! Convective-scale snowfall rate
       ntprprd     , &!INTENT(out) :: ntprprd(pcols,pver)    ! net precip production in layer
       ntsnprd     , &!INTENT(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
       flxprec     , &!INTENT(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
       flxsnow       )!INTENT(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)

    !-----------------------------------------------------------------------
    ! Compute tendencies due to evaporation of rain from ZM scheme
    !--
    ! Compute the total precipitation and snow fluxes at the surface.
    ! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
    ! in the Zhang-MacFarlane parameterization.
    ! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
    !-----------------------------------------------------------------------


    !------------------------------Arguments--------------------------------
    INTEGER,INTENT(in) :: ncol             ! number of columns and chunk index
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: pverp                 ! number of vertical levels + 1

    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: t          ! temperature (K)
    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: pdel       ! layer thickness (Pa)
    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: q          ! water vapor (kg/kg)
    REAL(r8),INTENT(inout), DIMENSION(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
    REAL(r8),INTENT(inout), DIMENSION(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)



    REAL(r8), INTENT(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
    REAL(r8), INTENT(in   ) :: cldfrc(pcols,pver) ! cloud fraction
    REAL(r8), INTENT(in   ) :: deltat             ! time step

    REAL(r8), INTENT(inout) :: prec(pcols)        ! Convective-scale preciptn rate
    REAL(r8), INTENT(out)   :: snow(pcols)        ! Convective-scale snowfall rate
    !
    !---------------------------Local storage-------------------------------

    REAL(r8) :: est    (pcols,pver)    ! Saturation vapor pressure
    REAL(r8) :: fice   (pcols,pver)    ! ice fraction in precip production
    REAL(r8) :: fsnow_conv(pcols,pver) ! snow fraction in precip production
    REAL(r8) :: qsat   (pcols,pver)    ! saturation specific humidity
    REAL(r8),INTENT(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
    REAL(r8),INTENT(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    REAL(r8),INTENT(out) :: ntprprd(pcols,pver)    ! net precip production in layer
    REAL(r8),INTENT(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
    REAL(r8) :: work1                  ! temp variable (pjr)
    REAL(r8) :: work2                  ! temp variable (pjr)

    REAL(r8) :: evpvint(pcols)         ! vertical integral of evaporation
    REAL(r8) :: evpprec(pcols)         ! evaporation of precipitation (kg/kg/s)
    REAL(r8) :: evpsnow(pcols)         ! evaporation of snowfall (kg/kg/s)
    REAL(r8) :: snowmlt(pcols)         ! snow melt tendency in layer
    REAL(r8) :: flxsntm(pcols)         ! flux of snow into layer, after melting

    REAL(r8) :: evplimit               ! temp variable for evaporation limits
    REAL(r8) :: rlat(pcols)

    INTEGER :: i,k                     ! longitude,level indices


    !-----------------------------------------------------------------------
!!$    ke = conke
!!$    ke = 2.0E-6


    ! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000.0_r8

    ! determine saturation vapor pressure
    CALL aqsat (t    ,pmid  ,est    ,qsat    ,pcols   , &
         ncol ,pver  ,1       ,pver    )

    ! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    CALL cldwat_fice(ncol,pcols,pver ,t, fice, fsnow_conv)

    ! zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0.0_r8
    flxsnow(:ncol,1) = 0.0_r8
    evpvint(:ncol)   = 0.0_r8

    DO k = 1, pver
       DO i = 1, ncol

          ! Melt snow falling into layer, if necessary. 
          IF (t(i,k) > tmelt) THEN
             flxsntm(i) = 0.0_r8
             snowmlt(i) = flxsnow(i,k) * gravit/ pdel(i,k)
          ELSE
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0.0_r8
          END IF

          ! relative humidity depression must be > 0 for evaporation
          evplimit = MAX(1.0_r8 - q(i,k)/qsat(i,k), 0.0_r8)

          ! total evaporation depends on flux in the top of the layer
          ! flux prec is the net production above layer minus evaporation into environmet
          evpprec(i) = ke * (1.0_r8 - cldfrc(i,k)) * evplimit * SQRT(flxprec(i,k))
          !**********************************************************
!!$           evpprec(i) = 0.    ! turn off evaporation for now
          !**********************************************************

          ! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
          ! Currently does not include heating/cooling change to qsat
          evplimit   = MAX(0.0_r8, (qsat(i,k)-q(i,k)) / deltat)

          ! Don't evaporate more than is falling into the layer - do not evaporate rain formed
          ! in this layer but if precip production is negative, remove from the available precip
          ! Negative precip production occurs because of evaporation in downdrafts.
!!$          evplimit   = flxprec(i,k) * gravit / pdel(i,k) + min(prdprec(i,k), 0.)
          evplimit   = MIN(evplimit, flxprec(i,k) * gravit / pdel(i,k))

          ! Total evaporation cannot exceed input precipitation
          evplimit   = MIN(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

          evpprec(i) = MIN(evplimit, evpprec(i))

          ! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          IF (flxprec(i,k) > 0.0_r8) THEN
             !            evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(i,k)
             !            prevent roundoff problems
             work1 = MIN(MAX(0.0_r8,flxsntm(i)/flxprec(i,k)),1.0_r8)
             evpsnow(i) = evpprec(i) * work1
          ELSE
             evpsnow(i) = 0.0_r8
          END IF

          ! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

          ! net precip production is production - evaporation
          ntprprd(i,k) = prdprec(i,k) - evpprec(i)
          ! net snow production is precip production * ice fraction - evaporation - melting
          !pjrworks ntsnprd(i,k) = prdprec(i,k)*fice(i,k) - evpsnow(i) - snowmlt(i)
          !pjrwrks2 ntsnprd(i,k) = prdprec(i,k)*fsnow_conv(i,k) - evpsnow(i) - snowmlt(i)
          ! the small amount added to flxprec in the work1 expression has been increased from 
          ! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
          ! scheme to be used for small flxprec amounts.  This is to address error growth problems.
          !#ifdef PERGRO
          IF(PERGRO)THEN
             work1 = MIN(MAX(0.0_r8,flxsnow(i,k)/(flxprec(i,k)+8.64e-11_r8)),1.0_r8)
             !#else
          ELSE
             IF (flxprec(i,k).GT.0.0_r8) THEN
                work1 = MIN(MAX(0.0_r8,flxsnow(i,k)/flxprec(i,k)),1.0_r8)
             ELSE
                work1 = 0.0_r8
             ENDIF
             !#endif
          ENDIF
          work2 = MAX(fsnow_conv(i,k), work1)
          IF (snowmlt(i).GT.0.0_r8) work2 = 0.0_r8
          !         work2 = fsnow_conv(i,k)
          ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)

          ! precipitation fluxes
          flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/gravit
          flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/gravit

          ! protect against rounding error
          flxprec(i,k+1) = MAX(flxprec(i,k+1), 0.0_r8)
          flxsnow(i,k+1) = MAX(flxsnow(i,k+1), 0.0_r8)
          ! more protection (pjr)
          !         flxsnow(i,k+1) = min(flxsnow(i,k+1), flxprec(i,k+1))

          ! heating (cooling) and moistening due to evaporation 
          ! - latent heat of vaporization for precip production has already been accounted for
          ! - snow is contained in prec
          tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice
          tend_q(i,k) = evpprec(i)
       END DO
    END DO

    ! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000.0_r8
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000.0_r8

    !**********************************************************
!!$    tend_s(:ncol,:)   = 0.      ! turn heating off
    !**********************************************************

  END SUBROUTINE zm_conv_evap





  SUBROUTINE mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhowtr,limcnv_in ,jlat)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Initialize moist convective mass flux procedure common block, cmfmca
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: rair              ! gas constant for dry air
    REAL(r8), INTENT(in) :: cpair             ! specific heat of dry air
    REAL(r8), INTENT(in) :: gravit            ! acceleration due to gravity
    REAL(r8), INTENT(in) :: latvap            ! latent heat of vaporization
    REAL(r8), INTENT(in) :: rhowtr            ! density of liquid water (STP)
    INTEGER,INTENT(in)   :: limcnv_in       ! top interface level limit for convection
    INTEGER,INTENT(in)   :: jlat
    !
    !-----------------------------------------------------------------------
    !
    ! Initialize physical constants for moist convective mass flux procedure
    !
    cp     = cpair         ! specific heat of dry air
    hlat   = latvap        ! latent heat of vaporization
    grav   = gravit        ! gravitational constant
    rgas   = rair          ! gas constant for dry air
    rhoh2o = rhowtr        ! density of liquid water (STP)

    limcnv = limcnv_in

    !
    ! Initialize free parameters for moist convective mass flux procedure
    !
    IF (dycore_is('LR'))THEN
       IF ( get_resolution(jlat) == '1x1.25' ) THEN
          cmftau = 3600.0_r8
          c0 = 3.5E-3_r8
          ke = 1.0E-6_r8
       ELSE IF ( get_resolution(jlat) == '4x5' ) THEN
          cmftau = 1800.0_r8         ! characteristic adjustment time scale
          c0     = 2.0e-4_r8        ! rain water autoconversion coeff (1/m)
          ke = 1.0E-6_r8
       ELSE
          cmftau = 1800.0_r8         ! characteristic adjustment time scale
          c0     = 1.0e-4_r8        ! rain water autoconversion coeff (1/m)
          ke = 1.0E-6_r8
       ENDIF
    ELSE
       IF(get_resolution(jlat) == 'T85')THEN
          cmftau = 1800.0_r8         ! characteristic adjustment time scale
          c0     = 1.0e-4_r8        ! rain water autoconversion coeff (1/m)
          ke = 1.0E-6_r8
       ELSEIF(get_resolution(jlat) == 'T31')THEN
          cmftau = 1800.0_r8         ! characteristic adjustment time scale
          c0     = 5.0e-4_r8        ! rain water autoconversion coeff (1/m)
          ke = 1.0E-6_r8
       ELSE
          cmftau = 1800.0_r8         ! characteristic adjustment time scale
          c0     = 2.0e-4_r8        ! rain water autoconversion coeff (1/m)
          ke = 3.0E-6_r8
       ENDIF
    ENDIF
    dzmin  = 0.0_r8           ! minimum cloud depth to precipitate (m)
    betamn = 0.10_r8          ! minimum overshoot parameter


    tpmax  = 1.50_r8          ! maximum acceptable t perturbation (deg C)
    shpmax = 1.50e-3_r8       ! maximum acceptable q perturbation (g/g)
    rlxclm = .TRUE.        ! logical variable to specify that relaxation
    !                                time scale should applied to column as
    !                                opposed to triplets individually
    !
    ! Initialize miscellaneous (frequently used) constants
    !
    rhlat  = 1.0_r8/hlat      ! reciprocal latent heat of vaporization
    rcp    = 1.0_r8/cp        ! reciprocal specific heat of dry air
    rgrav  = 1.0_r8/grav      ! reciprocal gravitational constant
    !
    ! Initialize diagnostic location information for moist convection scheme
    !
    iloc   = 1             ! longitude point for diagnostic info
    jloc   = 1             ! latitude  point for diagnostic info
    nsloc  = 1             ! nstep value at which to begin diagnostics
    !
    ! Initialize other miscellaneous parameters
    !
    tiny   = 1.0e-36_r8       ! arbitrary small number (scalar transport)
    eps    = 1.0e-13_r8       ! convergence criteria (machine dependent)
    !
    RETURN
  END SUBROUTINE mfinti

  SUBROUTINE gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   , &
       latvap  ,latice  ,rh2o    ,cpair   ,tmeltx   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Builds saturation vapor pressure table for later lookup procedure.
    ! 
    ! Method: 
    ! Uses Goff & Gratch (1946) relationships to generate the table
    ! according to a set of free parameters defined below.  Auxiliary
    ! routines are also included for making rapid estimates (well with 1%)
    ! of both es and d(es)/dt for the particular table configuration.
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use pmgrid, only: masterproc
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: tmn           ! Minimum temperature entry in es lookup table
    REAL(r8), INTENT(in) :: tmx           ! Maximum temperature entry in es lookup table
    REAL(r8), INTENT(in) :: epsil         ! Ratio of h2o to dry air molecular weights
    REAL(r8), INTENT(in) :: trice         ! Transition range from es over range to es over ice
    REAL(r8), INTENT(in) :: latvap        ! Latent heat of vaporization
    REAL(r8), INTENT(in) :: latice        ! Latent heat of fusion
    REAL(r8), INTENT(in) :: rh2o          ! Gas constant for water vapor
    REAL(r8), INTENT(in) :: cpair         ! Specific heat of dry air
    REAL(r8), INTENT(in) :: tmeltx        ! Melting point of water (K)
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(r8) t             ! Temperature
    INTEGER n          ! Increment counter
    INTEGER lentbl     ! Calculated length of lookup table
    INTEGER itype      ! Ice phase: 0 -> no ice phase
    !            1 -> ice phase, no transitiong
    !           -x -> ice phase, x degree transition
    LOGICAL ip         ! Ice phase logical flag
    !
    !-----------------------------------------------------------------------
    !
    ! Set es table parameters
    !
    !   tmin   = tmn       ! Minimum temperature entry in table
    !   tmax   = tmx       ! Maximum temperature entry in table
    !   ttrice = trice     ! Trans. range from es over h2o to es over ice
    !   icephs = ip        ! Ice phase (true or false)
    !
    ! Set physical constants required for es calculation
    !
    epsqs  = epsil
    hlatv  = latvap
    hlatf  = latice
    rgasv  = rh2o
    cp     = cpair
    !   tmelt  = tmeltx
    !
    lentbl = INT(tmax-tmin+2.000001_r8)
    IF (lentbl .GT. plenest) THEN
       WRITE(6,9000) tmax, tmin, plenest
       CALL endrun ('GESTBL')    ! Abnormal termination
    END IF
    !
    ! Begin building es table.
    ! Check whether ice phase requested.
    ! If so, set appropriate transition range for temperature
    !
    IF (icephs) THEN
       IF (ttrice /= 0.0_r8) THEN
          itype = -ttrice
       ELSE
          itype = 1
       END IF
    ELSE
       itype = 0
    END IF
    !
    t = tmin - 1.0_r8
    DO n=1,lentbl
       t = t + 1.0_r8
       CALL gffgch(t,estbl(n),itype)
    END DO
    !
    DO n=lentbl+1,plenest
       estbl(n) = -99999.0_r8
    END DO
    !
    ! Table complete -- Set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
    ! is valid in the range -40 < t < 0 (degrees C).
    !
    !                  --- Degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_r8
    pcf(2) = -5.47288442819e+00_r8
    pcf(3) = -3.67471858735e-01_r8
    pcf(4) = -8.95963532403e-03_r8
    pcf(5) = -7.78053686625e-05_r8
    !
    !                  --- Degree 6 approximation ---
    !
    !-----pcf(1) =  7.63285250063e-02
    !-----pcf(2) = -5.86048427932e+00
    !-----pcf(3) = -4.38660831780e-01
    !-----pcf(4) = -1.37898276415e-02
    !-----pcf(5) = -2.14444472424e-04
    !-----pcf(6) = -1.36639103771e-06
    !
    !   if (masterproc) then
    !      write(6,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
    !   end if

    RETURN
    !
9000 FORMAT('GESTBL: FATAL ERROR *********************************',/, &
         ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
         ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
         ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
    !
  END SUBROUTINE gestbl



  REAL(r8) FUNCTION estblf( td )
    !
    ! Saturation vapor pressure table lookup
    !
    REAL(r8), INTENT(in) :: td         ! Temperature for saturation lookup
    !
    REAL(r8) :: e       ! intermediate variable for es look-up
    REAL(r8) :: ai
    INTEGER  :: i
    !
    e = MAX(MIN(td,tmax),tmin)   ! partial pressure
    i = INT(e-tmin)+1
    ai = AINT(e-tmin)
    estblf = (tmin+ai-e+1.0_r8)* &
         estbl(i)-(tmin+ai-e)* &
         estbl(i+1)


  END FUNCTION estblf
  !---------------------------Statement functions-------------------------
  !   
  REAL(r8) FUNCTION qhalf(sh1,sh2,shbs1,shbs2)

    REAL(r8), INTENT(IN) :: sh1,sh2,shbs1,shbs2
    !qhalf(sh1,sh2,shbs1,shbs2) = min(max(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
    qhalf = MIN(MAX(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))

  END FUNCTION qhalf
  !
  !-----------------------------------------------------------------------

  SUBROUTINE gffgch(t       ,es      ,itype   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Computes saturation vapor pressure over water and/or over ice using
    ! Goff & Gratch (1946) relationships. 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! T (temperature), and itype are input parameters, while es (saturation
    ! vapor pressure) is an output parameter.  The input parameter itype
    ! serves two purposes: a value of zero indicates that saturation vapor
    ! pressures over water are to be returned (regardless of temperature),
    ! while a value of one indicates that saturation vapor pressures over
    ! ice should be returned when t is less than freezing degrees.  If itype
    ! is negative, its absolute value is interpreted to define a temperature
    ! transition region below freezing in which the returned
    ! saturation vapor pressure is a weighted average of the respective ice
    ! and water value.  That is, in the temperature range 0 => -itype
    ! degrees c, the saturation vapor pressures are assumed to be a weighted
    ! average of the vapor pressure over supercooled water and ice (all
    ! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use shr_kind_mod, only: r8 => shr_kind_r8
    !   use physconst, only: tmelt
    !   use abortutils, only: endrun

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: t          ! Temperature
    !
    ! Output arguments
    !
    INTEGER, INTENT(inout) :: itype   ! Flag for ice phase and associated transition

    REAL(r8), INTENT(out) :: es         ! Saturation vapor pressure
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(r8) e1         ! Intermediate scratch variable for es over water
    REAL(r8) e2         ! Intermediate scratch variable for es over water
    REAL(r8) eswtr      ! Saturation vapor pressure over water
    REAL(r8) f          ! Intermediate scratch variable for es over water
    REAL(r8) f1         ! Intermediate scratch variable for es over water
    REAL(r8) f2         ! Intermediate scratch variable for es over water
    REAL(r8) f3         ! Intermediate scratch variable for es over water
    REAL(r8) f4         ! Intermediate scratch variable for es over water
    REAL(r8) f5         ! Intermediate scratch variable for es over water
    REAL(r8) ps         ! Reference pressure (mb)
    REAL(r8) t0         ! Reference temperature (freezing point of water)
    REAL(r8) term1      ! Intermediate scratch variable for es over ice
    REAL(r8) term2      ! Intermediate scratch variable for es over ice
    REAL(r8) term3      ! Intermediate scratch variable for es over ice
    REAL(r8) tr         ! Transition range for es over water to es over ice
    REAL(r8) ts         ! Reference temperature (boiling point of water)
    REAL(r8) weight     ! Intermediate scratch variable for es transition
    INTEGER itypo   ! Intermediate scratch variable for holding itype
    !
    !-----------------------------------------------------------------------
    !
    ! Check on whether there is to be a transition region for es
    !
    IF (itype < 0) THEN
       tr    = ABS(float(itype))
       itypo = itype
       itype = 1
    ELSE
       tr    = 0.0_r8
       itypo = itype
    END IF
    IF (tr > 40.0_r8) THEN
       WRITE(6,900) tr
       CALL endrun ('GFFGCH')                ! Abnormal termination
    END IF
    !
    IF(t < (tmelt - tr) .AND. itype == 1) go to 10
    !
    ! Water
    !
    ps = 1013.246_r8
    ts = 373.16_r8
    e1 = 11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*LOG10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = LOG10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8
    eswtr = es
    !
    IF(t >= tmelt .OR. itype == 0) go to 20
    !
    ! Ice
    !
10  CONTINUE
    t0    = tmelt
    term1 = 2.01889049_r8/(t0/t)
    term2 = 3.56654_r8*LOG(t0/t)
    term3 = 20.947031_r8*(t0/t)
    es    = 575.185606e10_r8*EXP(-(term1 + term2 + term3))
    !
    IF (t < (tmelt - tr)) go to 20
    !
    ! Weighted transition between water and ice
    !
    weight = MIN((tmelt - t)/tr,1.0_r8)
    es = weight*es + (1.0_r8 - weight)*eswtr
    !
20  CONTINUE
    itype = itypo
    RETURN
    !
900 FORMAT('GFFGCH: FATAL ERROR ******************************',/, &
         'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
         ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
         ' 40.0 DEGREES C',/, ' TR = ',f7.2)
    !
  END SUBROUTINE gffgch


  SUBROUTINE vqsatd(t       ,p       ,es      ,qs      ,gam      , &
       len     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg        ! true if in temperature transition region
    !
    INTEGER i           ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(r8) tc        ! temperature (in degrees C)
    REAL(r8) weight    ! weight for es transition from water to ice
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) tterm     ! account for d(es)/dT in transition region
    REAL(r8) desdt     ! d(es)/dT
    REAL(r8) epsqs
    !
    !-----------------------------------------------------------------------
    !
    epsqs = epsilo

    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = estblf(t(i))
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    DO i=1,len
       !
       ! Weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition
       ! range from ice to water also accounting for change of hlatv with t
       ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
       weight = MIN(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       IF (lflg) THEN
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       ELSE
          tterm = 0.0_r8
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    RETURN
    !
    ! No icephs or water to ice transition
    !
10  DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       IF (icephs) THEN
          hlatsb = hlatv + hlatf
       ELSE
          hlatsb = hlatv
       END IF
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd

  SUBROUTINE cldwat_fice(ncol,pcols,pver ,t, fice, fsnow)
    !
    ! Compute the fraction of the total cloud water which is in ice phase.
    ! The fraction depends on temperature only. 
    ! This is the form that was used for radiation, the code came from cldefr originally
    ! 
    ! Author: B. A. Boville Sept 10, 2002
    !  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    ! Arguments
    INTEGER,  INTENT(in)  :: ncol                 ! number of active columns
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in)  :: t(pcols,pver)        ! temperature

    REAL(r8), INTENT(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    REAL(r8), INTENT(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

    ! Local variables
    INTEGER :: i,k                                   ! loop indexes

    !-----------------------------------------------------------------------

    ! Define fractional amount of cloud that is ice
    DO k=1,pver
       DO i=1,ncol

          ! If warmer than tmax then water phase
          IF (t(i,k) > tmax_fice) THEN
             fice(i,k) = 0.0_r8

             ! If colder than tmin then ice phase
          ELSE IF (t(i,k) < tmin_fice) THEN
             fice(i,k) = 1.0_r8

             ! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          ELSE 
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          END IF

          ! snow fraction partitioning

          ! If warmer than tmax then water phase
          IF (t(i,k) > tmax_fsnow) THEN
             fsnow(i,k) = 0.0_r8

             ! If colder than tmin then ice phase
          ELSE IF (t(i,k) < tmin_fsnow) THEN
             fsnow(i,k) = 1.0_r8

             ! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          ELSE 
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          END IF

       END DO
    END DO

    RETURN
  END SUBROUTINE cldwat_fice

  SUBROUTINE aqsat(t       ,p       ,es      ,qs        ,ii      , &
       ILEN    ,kk      ,kstart  ,kend      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
    ! This routine is useful for evaluating only a selected region in the
    ! vertical.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER  , INTENT(in) :: ii             ! I dimension of arrays t, p, es, qs
    INTEGER  , INTENT(in) :: kk             ! K dimension of arrays t, p, es, qs
    INTEGER  , INTENT(in) :: ILEN           ! Length of vectors in I direction which
    INTEGER  , INTENT(in) :: kstart         ! Starting location in K direction
    INTEGER  , INTENT(in) :: kend           ! Ending location in K direction
    REAL(r8), INTENT(in) :: t(ii,kk)          ! Temperature
    REAL(r8), INTENT(in) :: p(ii,kk)          ! Pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(ii,kk)         ! Saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(ii,kk)         ! Saturation specific humidity
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(r8) epsqs      ! Ratio of h2o to dry air molecular weights 
    REAL(r8) omeps             ! 1 - 0.622
    INTEGER i, k           ! Indices
    !
    !-----------------------------------------------------------------------
    !
    epsqs = epsilo
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative values
          ! that can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE aqsat

  SUBROUTINE aqsatd(t       ,p       ,es      ,qs      ,gam     , &
       ii      ,ILEN    ,kk      ,kstart  ,kend    )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g).   
    ! 
    ! Method: 
    ! Differs from aqsat by also calculating and returning
    ! gamma (l/cp)*(d(qsat)/dT)
    ! Input arrays temperature and pressure (dimensioned ii,kk).
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ii            ! I dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: ILEN          ! Vector length in I direction
    INTEGER, INTENT(in) :: kk            ! K dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: kstart        ! Starting location in K direction
    INTEGER, INTENT(in) :: kend          ! Ending location in K direction

    REAL(r8), INTENT(in) :: t(ii,kk)         ! Temperature
    REAL(r8), INTENT(in) :: p(ii,kk)         ! Pressure

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(ii,kk)        ! Saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(ii,kk)        ! Saturation specific humidity
    REAL(r8), INTENT(out) :: gam(ii,kk)       ! (l/cp)*(d(qs)/dt)
    !
    !---------------------------Local workspace-----------------------------
    !
    LOGICAL lflg          ! True if in temperature transition region
    INTEGER i             ! i index for vector calculations
    INTEGER k             ! k index
    REAL(r8) omeps            ! 1. - 0.622
    REAL(r8) trinv            ! Reciprocal of ttrice (transition range)
    REAL(r8) tc               ! Temperature (in degrees C)
    REAL(r8) weight           ! Weight for es transition from water to ice
    REAL(r8) hltalt           ! Appropriately modified hlat for T derivatives
    REAL(r8) hlatsb           ! hlat weighted in transition region
    REAL(r8) hlatvp           ! hlat modified for t changes above freezing
    REAL(r8) tterm            ! Account for d(es)/dT in transition region
    REAL(r8) desdt            ! d(es)/dT
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative qs
          ! values which can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    !
    DO k=kstart,kend
       DO i=1,ILEN
          !
          ! Weighting of hlat accounts for transition from water to ice
          ! polynomial expression approximates difference between es over
          ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
          ! -40): required for accurate estimate of es derivative in transition
          ! range from ice to water also accounting for change of hlatv with t
          ! above freezing where constant slope is given by -2369 j/(kg c) =cpv - cw
          !
          tc     = t(i,k) - tmelt
          lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
          weight = MIN(-tc*trinv,1.0_r8)
          hlatsb = hlatv + weight*hlatf
          hlatvp = hlatv - 2369.0_r8*tc
          IF (t(i,k) < tmelt) THEN
             hltalt = hlatsb
          ELSE
             hltalt = hlatvp
          END IF
          IF (lflg) THEN
             tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
          ELSE
             tterm = 0.0_r8
          END IF
          desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
          gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
          IF (qs(i,k) == 1.0_r8) gam(i,k) = 0.0_r8
       END DO
    END DO
    !
    go to 20
    !
    ! No icephs or water to ice transition
    !
10  DO k=kstart,kend
       DO i=1,ILEN
          !
          ! Account for change of hlatv with t above freezing where
          ! constant slope is given by -2369 j/(kg c) = cpv - cw
          !
          hlatvp = hlatv - 2369.0_r8*(t(i,k)-tmelt)
          IF (icephs) THEN
             hlatsb = hlatv + hlatf
          ELSE
             hlatsb = hlatv
          END IF
          IF (t(i,k) < tmelt) THEN
             hltalt = hlatsb
          ELSE
             hltalt = hlatvp
          END IF
          desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
          gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
          IF (qs(i,k) == 1.0_r8) gam(i,k) = 0.0_r8
       END DO
    END DO
    !
20  RETURN
  END SUBROUTINE aqsatd
  SUBROUTINE geopotential_dse(                                &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , ncol   ,pcols, pver, pverp          )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute the temperature  and geopotential height (above the surface) at the
    ! midpoints and interfaces from the input dry static energy and pressures.
    !
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    INTEGER, INTENT(in) :: ncol                  ! Number of longitudes
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    INTEGER, INTENT(in) :: pverp
    REAL(r8), INTENT(in) :: piln (pcols,pverp)   ! Log interface pressures
    REAL(r8), INTENT(in) :: pmln (pcols,pver)    ! Log midpoint pressures
    REAL(r8), INTENT(in) :: pint (pcols,pverp)   ! Interface pressures
    REAL(r8), INTENT(in) :: pmid (pcols,pver)    ! Midpoint pressures
    REAL(r8), INTENT(in) :: pdel (pcols,pver)    ! layer thickness
    REAL(r8), INTENT(in) :: rpdel(pcols,pver)    ! inverse of layer thickness
    REAL(r8), INTENT(in) :: dse  (pcols,pver)    ! dry static energy
    REAL(r8), INTENT(in) :: q    (pcols,pver)    ! specific humidity
    REAL(r8), INTENT(in) :: phis (pcols)         ! surface geopotential
    REAL(r8), INTENT(in) :: rair                 ! Gas constant for dry air
    REAL(r8), INTENT(in) :: gravit               ! Acceleration of gravity
    REAL(r8), INTENT(in) :: cpair                ! specific heat at constant p for dry air
    REAL(r8), INTENT(in) :: zvir                 ! rh2o/rair - 1

    ! Output arguments

    REAL(r8), INTENT(out) :: t(pcols,pver)       ! temperature
    REAL(r8), INTENT(out) :: zi(pcols,pverp)     ! Height above surface at interfaces
    REAL(r8), INTENT(out) :: zm(pcols,pver)      ! Geopotential height at mid level
    !
    !---------------------------Local variables-----------------------------
    !
    LOGICAL  :: fvdyn              ! finite volume dynamics
    INTEGER  :: i,k                ! Lon, level, level indices
    REAL(r8) :: hkk(pcols)         ! diagonal element of hydrostatic matrix
    REAL(r8) :: hkl(pcols)         ! off-diagonal element
    REAL(r8) :: rog                ! Rair / gravit
    REAL(r8) :: tv                 ! virtual temperature
    REAL(r8) :: tvfac              ! Tv/T
    !
    !-----------------------------------------------------------------------
    rog = rair/gravit

    ! Set dynamics flag
    fvdyn = dycore_is ('LR')

    ! The surface height is zero by definition.
    DO i = 1,ncol
       zi(i,pverp) = 0.0_r8
    END DO

    ! Compute the virtual temperature, zi, zm from bottom up
    ! Note, zi(i,k) is the interface above zm(i,k)
    DO k = pver, 1, -1

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

       ! Now compute tv, t, zm, zi
       DO i = 1,ncol
          tvfac   = 1.0_r8 + zvir * q(i,k)
          tv      = (dse(i,k) - phis(i) - gravit*zi(i,k+1)) / ((cpair / tvfac) + rair*hkk(i))

          t (i,k) = tv / tvfac

          zm(i,k) = zi(i,k+1) + rog * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv * hkl(i)
       END DO
    END DO

    RETURN
  END SUBROUTINE geopotential_dse


  !===============================================================================
  SUBROUTINE geopotential_t(                                 &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   ,pcols, pver, pverp)

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
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    INTEGER, INTENT(in) :: pverp
    REAL(r8), INTENT(in) :: piln (pcols,pverp)   ! Log interface pressures
    REAL(r8), INTENT(in) :: pmln (pcols,pver)    ! Log midpoint pressures
    REAL(r8), INTENT(in) :: pint (pcols,pverp)   ! Interface pressures
    REAL(r8), INTENT(in) :: pmid (pcols,pver)    ! Midpoint pressures
    REAL(r8), INTENT(in) :: pdel (pcols,pver)    ! layer thickness
    REAL(r8), INTENT(in) :: rpdel(pcols,pver)    ! inverse of layer thickness
    REAL(r8), INTENT(in) :: t    (pcols,pver)    ! temperature
    REAL(r8), INTENT(in) :: q    (pcols,pver)    ! specific humidity
    REAL(r8), INTENT(in) :: rair                 ! Gas constant for dry air
    REAL(r8), INTENT(in) :: gravit               ! Acceleration of gravity
    REAL(r8), INTENT(in) :: zvir                 ! rh2o/rair - 1

    ! Output arguments

    REAL(r8), INTENT(out) :: zi(pcols,pverp)     ! Height above surface at interfaces
    REAL(r8), INTENT(out) :: zm(pcols,pver)      ! Geopotential height at mid level
    !
    !---------------------------Local variables-----------------------------
    !
    LOGICAL  :: fvdyn              ! finite volume dynamics
    INTEGER  :: i,k                ! Lon, level indices
    REAL(r8) :: hkk(pcols)         ! diagonal element of hydrostatic matrix
    REAL(r8) :: hkl(pcols)         ! off-diagonal element
    REAL(r8) :: rog                ! Rair / gravit
    REAL(r8) :: tv                 ! virtual temperature
    REAL(r8) :: tvfac              ! Tv/T
    !
    !-----------------------------------------------------------------------
    !
    rog = rair/gravit

    ! Set dynamics flag

    fvdyn = dycore_is ('LR')

    ! The surface height is zero by definition.

    DO i = 1,ncol
       zi(i,pverp) = 0.0_r8
    END DO

    ! Compute zi, zm from bottom up. 
    ! Note, zi(i,k) is the interface above zm(i,k)

    DO k = pver, 1, -1

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
  !end module geopotential
  SUBROUTINE qneg3 (subnam  ,ncol    ,ncold   ,lver    ,lconst_beg  , &
       lconst_end       ,qmin    ,q       )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Check moisture and tracers for minimum value, reset any below
    ! minimum value to minimum value and return information to allow
    ! warning message to be printed. The global average is NOT preserved.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Rosinski
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    CHARACTER*(*), INTENT(in) :: subnam ! name of calling routine

    INTEGER, INTENT(in) :: ncol         ! number of atmospheric columns
    INTEGER, INTENT(in) :: ncold        ! declared number of atmospheric columns
    INTEGER, INTENT(in) :: lver         ! number of vertical levels in column
    INTEGER, INTENT(in) :: lconst_beg   ! beginning constituent
    INTEGER, INTENT(in) :: lconst_end   ! ending    constituent

    REAL(r8), INTENT(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

    !
    ! Input/Output arguments
    !
    REAL(r8), INTENT(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER indx(ncol,lver)  ! array of indices of points < qmin
    INTEGER nval(lver)       ! number of points < qmin for 1 level
    INTEGER nvals            ! number of values found < qmin
    INTEGER i,ii,k           ! longitude, level indices
    INTEGER m                ! constituent index
    INTEGER iw,kw            ! i,k indices of worst violator

    LOGICAL found            ! true => at least 1 minimum violator found

    REAL(r8) worst           ! biggest violator
    !
    !-----------------------------------------------------------------------
    !
    DO m=lconst_beg,lconst_end
       nvals = 0
       found = .FALSE.
       worst = 1.0e35_r8
       !
       ! Test all field values for being less than minimum value. Set q = qmin
       ! for all such points. Trace offenders and identify worst one.
       !
       !CDIR$ preferstream
       DO k=1,lver
          nval(k) = 0
          !CDIR$ prefervector
          DO i=1,ncol
             IF (q(i,k,m) < qmin(m)) THEN
                nval(k) = nval(k) + 1
                indx(nval(k),k) = i
             END IF
          END DO
       END DO

       DO k=1,lver
          IF (nval(k) > 0) THEN
             found = .TRUE.
             nvals = nvals + nval(k)
             DO ii=1,nval(k)
                i = indx(ii,k)
                IF (q(i,k,m) < worst) THEN
                   worst = q(i,k,m)
                   kw = k
                   iw = i
                END IF
                q(i,k,m) = qmin(m)
             END DO
          END IF
       END DO
       !#if ( defined WACCM_GHG || defined WACCM_MOZART )
       !      if (found .and. abs(worst)>1.e-12) then
       !#else               
       !      if (found .and. abs(worst)>1.e-16) then
       !#endif
       !         write(6,9000)subnam,m,idx,nvals,qmin(m),worst,iw,kw
       !      end if
    END DO
    !
    RETURN
9000 FORMAT(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i3, &
         ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
         1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
  END SUBROUTINE qneg3
  !==============================================================================================

  CHARACTER*3 FUNCTION cnst_get_type_byind (ind,ppcnst)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Get the type of a constituent 
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author:  P. J. Rasch
    ! 
    !-----------------------------Arguments---------------------------------
    !
    INTEGER, INTENT(in)   :: ind    ! global constituent index (in q array)
    INTEGER, INTENT(in)   :: ppcnst
    !---------------------------Local workspace-----------------------------
    INTEGER :: m                                   ! tracer index

    !-----------------------------------------------------------------------

    IF (ind.LE.ppcnst) THEN
       cnst_get_type_byind = cnst_type(ind)
    ELSE
       ! Unrecognized name
       WRITE(6,*) 'CNST_GET_TYPE_BYIND, ind:', ind
       CALL endrun
    ENDIF


  END FUNCTION cnst_get_type_byind

  LOGICAL FUNCTION dycore_is (name)
    !
    ! Input arguments
    !
    CHARACTER(len=*), INTENT(in) :: name      
    IF (name == 'eul' .OR. name == 'EUL') THEN
       dycore_is = .TRUE.
    ELSE
       dycore_is = .FALSE.
    END IF
    RETURN
  END FUNCTION dycore_is

  CHARACTER(len=7) FUNCTION get_resolution(plat)

    INTEGER, INTENT(IN   ) ::  plat! number of latitudes

    SELECT CASE ( plat )
    CASE ( 8 )
       get_resolution = 'T5'
    CASE ( 32 )
       get_resolution = 'T21'
    CASE ( 48 )
       get_resolution = 'T31'
    CASE ( 64 )
       get_resolution = 'T42'
    CASE ( 96 )
       get_resolution = 'T62'
    CASE ( 128 )
       get_resolution = 'T85'
    CASE ( 256 )
       get_resolution = 'T170'
    CASE DEFAULT
       get_resolution = 'UNKNOWN'
    END SELECT

    RETURN
  END FUNCTION get_resolution
  SUBROUTINE endrun (msg)
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Abort the model for abnormal termination
    !
    ! Author: CCM Core group
    !
    !-----------------------------------------------------------------------
    ! $Id: abortutils.F90,v 1.1.2.4 2004/09/17 16:59:22 mvr Exp $
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(in), OPTIONAL :: msg    ! string to be printed

    IF (PRESENT (msg)) THEN
       WRITE(6,*)'ENDRUN:', msg
    ELSE
       WRITE(6,*)'ENDRUN: called without a message string'
    END IF
    STOP
  END SUBROUTINE endrun

END MODULE Shall_JHack
!PROGRAM MAIN
!  USE Shall_JHack
!END PROGRAM Main
