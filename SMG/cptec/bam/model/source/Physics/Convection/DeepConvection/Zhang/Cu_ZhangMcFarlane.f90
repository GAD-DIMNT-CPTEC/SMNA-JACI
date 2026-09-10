MODULE Cu_ZhangMcFarlane
  IMPLICIT NONE
SAVE

  PRIVATE
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
  REAL(r8),PARAMETER :: SHR_CONST_TKFRZ  = 273.16_r8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(r8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air ~ kg/kmole
  REAL(r8),PARAMETER :: SHR_CONST_LATVAP = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_LATICE = 3.337e5_r8      ! latent heat of fusion ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(r8),PARAMETER :: SHR_CONST_G      = 9.80616_r8      ! acceleration of gravity ~ m/s^2
  REAL(r8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(r8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
  real(R8),parameter :: SHR_CONST_CPFW    = 4.188e3_R8      
  real(R8),parameter :: SHR_CONST_CPWV    = 1.810e3_R8      

  ! Constants for Earth

  REAL(r8), PUBLIC, PARAMETER :: gravit = shr_const_g      ! gravitational acceleration

  ! Constants for air
  real(r8), public, parameter :: cpliq       = shr_const_cpfw       
  real(r8), public, parameter :: cpwv        = shr_const_cpwv  

  REAL(r8), PUBLIC, PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: rair = shr_const_rdair    ! Gas constant for dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: zvir = SHR_CONST_RWV/rair - 1          ! rh2o/rair - 1

  ! Constants for water

  REAL(r8), PUBLIC, PARAMETER :: tmelt = shr_const_tkfrz   ! Freezing point of water
  REAL(r8), PUBLIC, PARAMETER :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights 
  REAL(r8), PUBLIC, PARAMETER :: latvap = shr_const_latvap ! Latent heat of vaporization
  REAL(r8), PUBLIC, PARAMETER :: latice = shr_const_latice ! Latent heat of fusion
  REAL(r8), PUBLIC, PARAMETER :: rh2o  =SHR_CONST_RWV   !! Gas constant for water vapor
  REAL(r8), PUBLIC, PARAMETER :: trice  =  20.00_r8         ! Trans range from es over h2o to es over ice
  REAL(r8), PUBLIC, PARAMETER :: ttrice=trice
  REAL(r8), PUBLIC, PARAMETER :: tmax_fice = tmelt - 10.0_r8       ! max temperature for cloud ice formation
  REAL(r8), PUBLIC, PARAMETER :: tmin_fice = tmax_fice - 30.0_r8   ! min temperature for cloud ice formation
  REAL(r8), PUBLIC, PARAMETER :: tmax_fsnow = tmelt            ! max temperature for transition to convective snow
  REAL(r8), PUBLIC, PARAMETER :: tmin_fsnow = tmelt-5.0_r8         ! min temperature for transition to convective snow
  ! if (cnst_get_type_byind(m).eq.'dry') then
  !character(LEN=3) function cnst_get_type_byind (ind)
  LOGICAL, PARAMETER :: cnst_need_pdeldry=.FALSE.

  LOGICAL, PARAMETER :: masterproc=.FALSE.
  LOGICAL, PARAMETER :: PERGRO=.FALSE.
  !------------wv_saturation-------------
  !
  ! Data
  !
  INTEGER, PARAMETER :: plenest = 250! length of saturation vapor pressure table

  REAL(r8)           :: estbl(plenest)      ! table values of saturation vapor pressure
  REAL(r8),PARAMETER :: tmn  = 173.16_r8          ! Minimum temperature entry in table
  REAL(r8),PARAMETER :: tmx  = 375.16_r8          ! Maximum temperature entry in table
  REAL(r8),PARAMETER :: tmin=tmn       ! min temperature (K) for table
  REAL(r8),PARAMETER :: tmax= tmx      ! max temperature (K) for table
  LOGICAL  ,PARAMETER :: icephs=.TRUE.  ! false => saturation vapor press over water only
  REAL(r8)            :: pcf(6)     ! polynomial coeffs -> es transition water to ice
  INTEGER ,PARAMETER  :: ixcldliq=2
  INTEGER ,PARAMETER  :: ixcldice=3
  
  logical :: no_deep_pbl=.false. ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL 

  real(r8),parameter ::  tiedke_add = 0.5_r8   

  PUBLIC :: Init_Cu_ZhangMcFarlane
  PUBLIC :: RunCu_ZhangMcFarlane
  !
  ! Private data
  !
  !bundy no one uses these   public rl, cpres, capelmt
  REAL(r8) rl         ! wg latent heat of vaporization.
  REAL(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
  REAL(r8), PARAMETER :: capelmt = 70.0_r8  ! threshold value for cape for deep convection.
  REAL(r8) :: ke                        ! Tunable evaporation efficiency
  REAL(r8) c0
  REAL(r8) :: tau   ! convective time scale
  REAL(r8),PARAMETER :: a = 21.656_r8 
  REAL(r8),PARAMETER :: b = 5418.0_r8 
  REAL(r8),PARAMETER :: c1 = 6.112_r8 
  REAL(r8),PARAMETER :: c2 = 17.67_r8 
  REAL(r8),PARAMETER :: c3 = 243.5_r8 

  REAL(r8) :: tfreez
  REAL(r8) :: eps1
  !bundy the following used to live in moistconvection.F90
  REAL(r8) :: rgrav       ! reciprocal of grav
  REAL(r8) :: rgas        ! gas constant for dry air
  REAL(r8) :: grav        ! = gravit
  REAL(r8) :: cp          ! = cpres = cpair

  INTEGER  :: limcnv       ! top interface level limit for convection
  INTEGER, PARAMETER  :: ppcnst=3
  REAL(r8) :: qmin(3)
  CHARACTER*3, PUBLIC :: cnst_type(ppcnst)          ! wet or dry mixing ratio
  CHARACTER*3, PARAMETER :: mixtype(ppcnst)=(/'wet','wet', 'wet'/) ! mixing ratio type (dry, wet)
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

!  REAL(r8), POINTER, DIMENSION(:,:,:) :: state_cld
!  REAL(r8), POINTER, DIMENSION(:,:,:) :: state_qliq           ! wg grid slice of cloud liquid water.
!  REAL(r8), POINTER, DIMENSION(:,:,:) :: state_rprd         ! rain production rate
  REAL(r8), POINTER, DIMENSION(:,:,:,:) :: state_fracis  ! fraction of transported species that are insoluble

 ! REAL(r8), POINTER, DIMENSION(:,:)  :: state_jctop
 ! REAL(r8), POINTER, DIMENSION(:,:)  :: state_jcbot

  !
  ! Private module data
  !

  !REAL(r8), POINTER, DIMENSION(:,:,:) :: mu  !(pcols,pver,begchunk:endchunk)
  !REAL(r8), POINTER, DIMENSION(:,:,:) :: eu  !(pcols,pver,begchunk:endchunk)
  !REAL(r8), POINTER, DIMENSION(:,:,:) :: du  !(pcols,pver,begchunk:endchunk)
  !REAL(r8), POINTER, DIMENSION(:,:,:) :: md  !(pcols,pver,begchunk:endchunk)
  !REAL(r8), POINTER, DIMENSION(:,:,:) :: ed  !(pcols,pver,begchunk:endchunk)
  !REAL(r8), POINTER, DIMENSION(:,:,:) :: dp  !(pcols,pver,begchunk:endchunk) 
  ! wg layer thickness in mbs (between upper/lower interface).
  !REAL(r8), POINTER, DIMENSION(:,:)   :: dsubcld  !(pcols,begchunk:endchunk)
  ! wg layer thickness in mbs between lcl and maxi.

  !INTEGER, POINTER, DIMENSION(:,:) :: jt   !(pcols,begchunk:endchunk)
  ! wg top  level index of deep cumulus convection.
  !INTEGER, POINTER, DIMENSION(:,:) :: maxg !(pcols,begchunk:endchunk)
  ! wg gathered values of maxi.
  !INTEGER, POINTER, DIMENSION(:,:) :: ideep !(pcols,begchunk:endchunk)               
  ! w holds position of gathered points vs longitude index

  !INTEGER, POINTER, DIMENSION(:) :: lengath !(begchunk:endchunk)

  REAL(KIND=r8)            :: ae  (2)
  REAL(KIND=r8)            :: be  (2)
  REAL(KIND=r8)            :: ht  (2)


CONTAINS
  SUBROUTINE Init_Cu_ZhangMcFarlane(dt,a_hybr,b_hybr,ibMax,kMax,jbMax,jMax)

    INTEGER,INTENT(IN   ) :: jMax           ! number of latitudes
    INTEGER,INTENT(IN   ) :: ibMax
    INTEGER,INTENT(IN   ) :: kMax
    INTEGER,INTENT(IN   ) :: jbMax
    INTEGER              :: k 
    REAL(KIND=r8), INTENT(IN   ) :: dt
    REAL(KIND=r8), INTENT(IN   ) :: a_hybr(kMax+1)
    REAL(KIND=r8), INTENT(IN   ) :: b_hybr(kMax+1)

    REAL(KIND=r8) :: hypi (kMax+1)
    INTEGER :: limcnv_in       ! top interface level limit for convection

    !---------------------------Local workspace-----------------------------

    LOGICAL ip           ! Ice phase (true or false)
    INTEGER :: ind    
    qmin=1.0e-12_r8

    DO k=1,kMax+1
!      hypi(k) =  si_in    (kMax+2-k)
!      SB  change to hybrid (already from top to bottom)
       hypi(k) =  a_hybr(k) / 1.e5_r8 + b_hybr(k)
    END DO
    limcnv_in=1
    CALL convect_deep_init(dt,hypi,limcnv_in,ibMax,kMax,jbMax,jMax)

    !-----------------------------------------------------------------------
    !
    ! Specify control parameters first
    !
    ip    = .TRUE.

    CALL gestbl(tmn     ,tmx     ,trice   ,ip      ,epsilo  , &
         latvap  ,latice  ,rh2o    ,cpair   ,tmelt )

    DO ind=1,ppcnst
       ! set constituent mixing ratio type
       !if ( present(mixtype) )then
       cnst_type(ind) = mixtype(ind) 
       !else
       !   cnst_type(ind) = 'wet'
       !end if
    END DO
    !
    ht(1)=latvap/cpair
    
    ht(2)=2.834e6_r8/cpair
    
    be(1)=0.622_r8*ht(1)/0.286_r8
    
    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)
    
    be(2)=0.622_r8*ht(2)/0.286_r8
    
    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)

  END SUBROUTINE Init_Cu_ZhangMcFarlane


  !=========================================================================================

  SUBROUTINE convect_deep_init(dt,hypi,limcnv_in,ibMax,kMax,jbMax,jMax)

    !----------------------------------------
    ! Purpose:  declare output fields, initialize variables needed by convection
    !----------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: jMax
    INTEGER, INTENT(INOUT) :: limcnv_in! top interface level limit for convection
    REAL(r8),INTENT(in) :: hypi(kMax+1)        ! reference pressures at interfaces
    REAL(r8),INTENT(in) :: dt
    INTEGER k, istat
    PRINT*,'convect_deep_init',jbMax,kMax,ibMax

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



!    ALLOCATE(state_cld    (1:ibMax,1:kMax,1:jbMax));state_cld=0.0_r8
    !ALLOCATE(state_qliq   (1:ibMax,1:kMax,1:jbMax));state_qliq  =0.0_r8    ! wg grid slice of cloud liquid water.
    !ALLOCATE(state_rprd   (1:ibMax,1:kMax,1:jbMax));state_rprd   =0.0_r8    ! rain production rate
    ALLOCATE(state_fracis (1:ibMax,1:kMax,1:jbMax,1:ppcnst)); state_fracis=0.0_r8 ! fraction of transported species that are insoluble

!    ALLOCATE(state_jctop  (1:ibMax,1:jbMax));state_jctop=0.0_r8
!    ALLOCATE(state_jcbot  (1:ibMax,1:jbMax));state_jcbot=0.0_r8

    !
    ! Allocate space for arrays private to this module
    !
!    ALLOCATE( mu     (1:ibMax,1:kMax,1:jbMax), stat=istat )  ;mu   =0.0_r8    
!    ALLOCATE( eu     (1:ibMax,1:kMax,1:jbMax), stat=istat )  ;eu  =0.0_r8     
!    ALLOCATE( du     (1:ibMax,1:kMax,1:jbMax), stat=istat )  ;du   =0.0_r8    
!    ALLOCATE( md     (1:ibMax,1:kMax,1:jbMax), stat=istat )  ;md   =0.0_r8    
!    ALLOCATE( ed     (1:ibMax,1:kMax,1:jbMax), stat=istat )  ;ed     =0.0_r8  
!    ALLOCATE( dp     (1:ibMax,1:kMax,1:jbMax), stat=istat )  ;dp     =0.0_r8  
!    ALLOCATE( dsubcld(1:ibMax,1:jbMax), stat=istat )       ;dsubcld =0.0_r8 
!    ALLOCATE( jt     (1:ibMax,1:jbMax), stat=istat )       ;jt      =0
!    ALLOCATE( maxg   (1:ibMax,1:jbMax), stat=istat )       ;maxg    =0
!    ALLOCATE( ideep  (1:ibMax,1:jbMax), stat=istat )       ;ideep   =0
!    ALLOCATE( lengath(1:jbMax), stat=istat )               ;lengath =0

    !
    ! Limit deep convection to regions below 40 mb or 40/1000 =  0.040 sigma
    ! Note this calculation is repeated in the shallow convection interface
    !
    limcnv_in = 0   ! null value to check against below
    IF (hypi(1) >= 0.04_r8) THEN
       limcnv_in = 1
    ELSE
       DO k=1,kMax
          IF (hypi(k) < 0.04_r8 .AND. hypi(k+1) >= 0.04_r8) THEN
             limcnv_in = k
             EXIT
          END IF
       END DO
       IF ( limcnv_in == 0 ) limcnv_in = kMax+1
    END IF

    IF (masterproc) THEN
       WRITE(6,*)'CONVECT_DEEP_INIT: Deep convection will be capped at intfc ',limcnv_in, &
            ' which is ',hypi(limcnv_in),' pascals'
    END IF

    CALL zm_convi(limcnv_in,jMax,dt)


  END SUBROUTINE convect_deep_init


  !=========================================================================================
  !subroutine RunCu_ZhangMcFarlane(state, ptend, tdt, pbuf)

  SUBROUTINE RunCu_ZhangMcFarlane( &
       pcols       , &! INTEGER , INTENT(in ) :: pcols                 ! number of columns (max)
       pverp       , &! INTEGER , INTENT(in ) :: pverp                 ! number of vertical levels + 1
       pver        , &! INTEGER , INTENT(in ) :: pver                  ! number of vertical levels
       latco       , &! INTEGER , INTENT(in ) :: latco                  ! number of latitudes
       pcnst       , &! INTEGER , INTENT(in ) :: pcnst        ! number of advected constituents (including water vapor)
       pnats       , &! INTEGER , INTENT(in ) :: pnats        ! number of non-advected constituents
       prsi        , &
       prsl        , &
       state_phis  , &! REAL(r8), INTENT(in ) :: state_phis   (pcols)    !(pcols)          ! surface geopotential
       state_t     , &! REAL(r8), INTENT(in ) :: state_t          (pcols,pver)!(pcols,pver)! temperature (K)
       state_qv    , &! REAL(r8), INTENT(in ) :: state_qv          (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       state_ql    , &! REAL(r8), INTENT(in ) :: state_ql          (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       state_qi    , &! REAL(r8), INTENT(in ) :: state_qi          (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       state_omega , &! REAL(r8), INTENT(in ) :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
       state_cld   , &! REAL(r8), INTENT(out) :: state_cld(ibMax,kMax)          !  cloud fraction
       prec        , &! real(r8), intent(out) :: prec(pcols)   ! total precipitation
       pblh        , &! real(r8), intent(in ) :: pblh(pcols)
       state_mcon  , &! real(r8), intent(out) :: mcon(pcols,pverp)
       tpert       , &! real(r8), intent(in ) :: tpert(pcols)
       state_dlf   , &! real(r8), intent(out) :: dlf(pcols,pver)! scattrd version of the detraining cld h2o tend
       state_zdu   , &! real(r8), intent(out) :: zdu(pcols,pver)
       rliq        , &! real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
       ztodt       , &! real(r8), intent(in ) :: ztodt                          ! 2 delta t (model time increment)
       snow        , &! real(r8), intent(out) :: snow(pcols)   ! snow from ZM convection 
       kctop1      , &
       kcbot1      , &
       kuo         , &
       dtdt        , &
       dqdt        , &
       dqldt       , &
       dqidt         )
         !,&
         !state   ,ptend_all   ,pbuf  )


    !use history,       only: outfld
    !use physics_types, only: physics_state, physics_ptend, physics_tend
    !use physics_types, only: physics_ptend_init,  physics_tend_init,physics_update
    !use physics_types, only: physics_state_copy
    !use physics_types, only: physics_ptend_sum

    !use phys_grid,     only: get_lat_p, get_lon_p
    !use time_manager,  only: get_nstep, is_first_step
    !use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
    !use constituents, only: pcnst, pnats,ppcnst, cnst_get_ind
    !use check_energy,    only: check_energy_chng
    !use physconst,       only: gravit

    ! Arguments
    !   type(physics_state), intent(in ) :: state          ! Physics state variables
    !   type(physics_ptend), intent(out) :: ptend_all          ! indivdual parameterization tendencies
    !   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer

    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                 ! number of vertical levels + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: latco
    INTEGER, INTENT(in) :: pcnst
    INTEGER, INTENT(in) :: pnats
    REAL(r8), INTENT(in)  ::    prsi    (1:pcols,1:pver+1)   ! 
    REAL(r8), INTENT(in)  ::    prsl    (1:pcols,1:pver)     ! 

    REAL(r8), INTENT(in)  :: state_phis   (pcols) 
    REAL(r8), INTENT(inout)  :: state_t          (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_qv     (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_ql     (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_qi     (pcols,pver)
    REAL(r8), INTENT(in)  :: state_omega  (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_cld    (pcols,pver)
    REAL(r8), INTENT(in) :: ztodt                          ! 2 delta t (model time increment)
    REAL(r8), INTENT(in) :: pblh(pcols)
    REAL(r8), INTENT(in) :: tpert(pcols)

    REAL(r8), INTENT(out) :: state_mcon(pcols,pverp)
    REAL(r8), INTENT(out) :: state_dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
    REAL(r8) :: pflx(pcols,pverp)  ! scattered precip flux at each level
    REAL(r8) :: cme(pcols,pver)
    REAL(r8), INTENT(out) :: state_zdu(pcols,pver)

    REAL(r8), INTENT(out) :: prec(pcols)   ! total precipitation
    REAL(r8), INTENT(out) :: snow(pcols)   ! snow from ZM convection 
    REAL(r8), INTENT(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
    INTEGER  , INTENT(OUT) ::  kctop1 (pcols)
    INTEGER  , INTENT(OUT) ::  kcbot1 (pcols)
    INTEGER  , INTENT(OUT) ::  kuo (pcols)

    REAL(KINd=r8), INTENT(OUT)  :: dtdt (pcols,pver)
    REAL(KIND=r8), INTENT(OUT)  :: dqdt (pcols,pver)
    REAL(KIND=r8), INTENT(OUT)  :: dqldt(pcols,pver)
    REAL(KIND=r8), INTENT(OUT)  :: dqidt(pcols,pver)

    ! Local variables

    INTEGER :: i,k,m,ll
    INTEGER :: ilon                      ! global longitude index of a column
    INTEGER :: ilat                      ! global latitude index of a column
    !integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
    INTEGER  :: itim, ifld  ! for physics buffer fields 
    REAL(r8) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
    REAL(r8) :: zdu(pcols,pver)

    REAL(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
    REAL(r8) ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
    REAL(r8) ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
    REAL(r8) flxprec(pcols,pverp)   ! evap outfld: Convective-scale flux of precip at interfaces (kg/m2/s)
    REAL(r8) flxsnow(pcols,pverp)   ! evap outfld: Convective-scale flux of snow   at interfaces (kg/m2/s)
    REAL(r8) :: mcon(pcols,pverp)

    ! physics types
    !    type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
    !    type(physics_tend ) :: tend          ! Physics tendencies (empty, needed for physics_update call)
    !    type(physics_ptend)  :: ptend_loc   ! package tendencies

    ! physics buffer fields 

    REAL(r8) :: cld    (pcols,pver)
    REAL(r8) :: ql     (pcols,pver)          ! wg grid slice of cloud liquid water.
    REAL(r8) :: rprd   (pcols,pver)          ! rain production rate
    REAL(r8) :: fracis(pcols,pver,ppcnst)  ! fraction of transported species that are insoluble

    REAL(r8) :: jctop(pcols)
    REAL(r8) :: jcbot(pcols)
    REAL(r8):: state_buf_t         (pcols,pver)  
    REAL(r8):: state_buf_omega   (pcols,pver) 
    REAL(r8):: state_buf_pmid         (pcols,pver) 
    REAL(r8):: state_buf_pdel         (pcols,pver) 
    REAL(r8):: state_buf_rpdel   (pcols,pver)
    REAL(r8):: state_buf_lnpmid  (pcols,pver)
    REAL(r8):: state_buf_pint         (pcols,pver+1)
    REAL(r8):: state_buf_lnpint  (pcols,pver+1) 
    REAL(r8):: state_buf_q         (pcols,pver,ppcnst)   
    REAL(r8):: ptend_loc_buf_q  (pcols,pver,ppcnst)
    REAL(r8):: state_buf_zm(pcols,pver)
    REAL(r8):: state_buf_zi(pcols,pver+1)
    REAL(r8):: state_buf_s(pcols,pver)
    REAL(r8):: ptend_loc_buf_cflx_srf(pcols,ppcnst)    
    REAL(r8):: ptend_loc_buf_cflx_top(pcols,ppcnst)    
    REAL(r8) :: tv(pcols,pver)
    REAL(r8) :: press(pcols,pver)
    REAL(r8) :: delz(pcols,pver)
    REAL(r8) :: diffpl1(pcols,pver)
    REAL(r8), TARGET :: rbyg
    LOGICAL , TARGET  :: aux2lq(1:ppcnst)
    REAL(r8), TARGET :: qnew(pcols,pver)
    REAL(r8), TARGET :: hkk(pcols)
    REAL(r8), TARGET :: hkl(pcols)
    REAL(r8), TARGET :: tvfac
    REAL(r8) :: mu  (pcols,pver)
    REAL(r8) :: eu  (pcols,pver)
    REAL(r8) :: du  (pcols,pver)
    REAL(r8) :: md  (pcols,pver)
    REAL(r8) :: ed  (pcols,pver)
    REAL(r8) :: dp  (pcols,pver) 
    REAL(r8) :: dsubcld  (pcols)
    INTEGER  :: jt   (pcols)
    INTEGER  :: maxg (pcols)
    INTEGER  :: ideep (pcols)               
    INTEGER  :: lengath 

    LOGICAL :: fvdyn

    IF(pcols < 1) RETURN
    state%ncol(latco)         = pcols
    
    state%ps  (1:pcols,latco) = prsi(1:pcols,1) !state_ps(1:pcols) 
    state%phis(1:pcols,latco) = state_phis(1:pcols) 
    state_mcon                = 0.0_r8
    state_dlf                 = 0.0_r8
    state_zdu                 = 0.0_r8

    !
    ! initialize 
    !
    ftem = 0.0_r8   
    DO i=1,pcols
       !state_buf_pint       (i,pver+1) = state_ps(i)*si(1)
       state_buf_pint       (i,pver+1) = prsi(i,1)
    END DO
    DO k=pver,1,-1
       DO i=1,pcols
          state_buf_pint    (i,k) = MAX(prsi    (i,pver+2-k)  ,0.0001_r8)
          !state_buf_pint    (i,k) = MAX(si(pver+2-k)*state_ps(i) ,0.0001_r8)
       END DO
    END DO

    DO k=1,pver+1
       DO i=1,pcols
          state_buf_lnpint(i,k) =  LOG(state_buf_pint  (i,k))
          mcon (i,pver+2-k)     =  state_mcon (i,k)
       END DO
    END DO

    state%pint     (1:pcols,1:pver+1,latco) = state_buf_pint(1:pcols,1:pver+1)
    state%lnpint   (1:pcols,1:pver+1,latco) = state_buf_lnpint(1:pcols,1:pver+1)

    DO k=1,pver
       DO i=1,pcols

          dtdt (i,pver+1-k)=0.0_r8
          dqdt (i,pver+1-k)=0.0_r8
          dqldt(i,pver+1-k)=0.0_r8
          dqidt(i,pver+1-k)=0.0_r8

          state_buf_t       (i,pver+1-k) = state_t          (i,k)
          state_buf_omega   (i,pver+1-k) = state_omega    (i,k)
          !state_buf_pmid    (i,pver+1-k) = sl(k)*state_ps (i)
          state_buf_pmid    (i,pver+1-k) = prsl(i,k)
          dlf               (i,pver+1-k) = state_dlf (i,k) 
          zdu               (i,pver+1-k) = state_zdu (i,k) 
       END DO
    END DO

    DO k=1,pver
       DO i=1,pcols          
          state_buf_pdel    (i,k) = MAX(state%pint(i,k+1,latco) - state%pint(i,k,latco),0.5_r8)
          state_buf_rpdel   (i,k) = 1.0_r8/MAX((state%pint(i,k+1,latco) - state%pint(i,k,latco)),0.5_r8)
          state_buf_lnpmid  (i,k) = LOG(state_buf_pmid(i,k))        
       END DO
    END DO

    state%t        (1:pcols,1:pver,latco)= state_buf_t     (1:pcols,1:pver) 
    state%omega    (1:pcols,1:pver,latco)= state_buf_omega (1:pcols,1:pver) 
    state%pmid     (1:pcols,1:pver,latco)= state_buf_pmid  (1:pcols,1:pver) 
    state%pdel     (1:pcols,1:pver,latco)= state_buf_pdel  (1:pcols,1:pver) 
    state%rpdel    (1:pcols,1:pver,latco)= state_buf_rpdel (1:pcols,1:pver) 
    state%lnpmid   (1:pcols,1:pver,latco)= state_buf_lnpmid(1:pcols,1:pver) 

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
       state%pmid  (1:pcols,1:pver,latco)     , state%pdel  (1:pcols,1:pver,latco)   , state%rpdel(1:pcols,1:pver,latco)     , &
       state%t     (1:pcols,1:pver,latco)     , state%q     (1:pcols,1:pver,latco,1) , rair   , gravit , zvir   ,          &
       state%zi    (1:pcols,1:pver+1,latco)   , state%zm    (1:pcols,1:pver,latco)   , pcols   , pver, pverp)

    fvdyn = dycore_is ('LR')
    DO k = pver, 1, -1
       ! First set hydrostatic elements consistent with dynamics
       IF (fvdyn) THEN
          DO i = 1,pcols
             hkl(i) = state%lnpmid(i,k+1,latco) - state%lnpmid(i,k,latco)
             hkk(i) = 1.0_r8  - state%pint (i,k,latco)* hkl(i)* state%rpdel(i,k,latco)
          END DO
       ELSE
          DO i = 1,pcols
             hkl(i) = state%pdel(i,k,latco)/state%pmid  (i,k,latco)
             hkk(i) = 0.5_r8  * hkl(i)
          END DO
       END IF
       ! Now compute s
       DO i = 1,pcols
          tvfac   = 1.0_r8 + zvir * state%q(i,k,latco,1) 
          state%s(i,k,latco) =  (state%t(i,k,latco)* cpair) + (state%t(i,k,latco) * tvfac * rair*hkk(i))  +  &
                                ( state%phis(i,latco) + gravit*state%zi(i,k+1,latco))           
       END DO
    END DO

    !-------------------------------------------------------------------------------
    ! This is for tendencies returned from individual parameterizations

    CALL physics_state_copy(pcols,latco,pver, pverp,ppcnst, cnst_need_pdeldry)   ! copy state to local state1.  
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)  ! initialize local ptend type
    CALL physics_ptend_init(ptend_all,latco,ppcnst,pver)  ! initialize output ptend type
    CALL physics_tend_init (tend,latco)                    ! tend here is just a null place holder
    !
    ! Associate pointers with physics buffer fields
    !
    DO k=1,pver
       DO i=1,pcols
          cld    (i,pver+1-k) =   state_cld    (i,k)
          ql     (i,pver+1-k) =   0.0_r8
          !rprd   (i,pver+1-k) =   state_rprd   (i,k,latco)
       END DO
    END DO
!    DO i=1,pcols
!  !     jctop  (i)=state_jctop  (i,latco)
!       jcbot  (i)=state_jcbot  (i,latco)
!    END DO
    DO ll=1,ppcnst
       DO k=1,pver
          DO i=1,pcols
             fracis (i,pver+1-k,ll) = state_fracis (i,k,latco,ll)
          END DO
       END DO
    END DO


    !
    ! Begin with Zhang-McFarlane (1996) convection parameterization
    !
    !  call t_startf ('zm_convr')
    CALL zm_convr(  pcols,pver, pcnst,pnats,pverp,&
         state%t      (1:pcols,1:pver,latco)   , &
         state%q      (1:pcols,1:pver,latco,1:pcnst+pnats)  , &
         prec         (1:pcols)  , &
         jctop        (1:pcols) , & ! REAL(r8), INTENT(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
         jcbot        (1:pcols) , & ! REAL(r8), INTENT(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
         pblh         (1:pcols) , &
         state%zm     (1:pcols,1:pver,latco)  , &
         state%phis   (1:pcols,latco)  , &
         state%zi     (1:pcols,1:pver+1,latco)  , &
         ptend_loc%q  (1:pcols,1:pver,latco,1)    , &
         ptend_loc%s  (1:pcols,1:pver,latco)           , &
         state%pmid   (1:pcols,1:pver,latco)            , &
         state%pint   (1:pcols,1:pver+1,latco)            , &
         state%pdel   (1:pcols,1:pver,latco)            , &
         !PK 0.5_r8*ztodt             , &
         ztodt                        ,&
         mcon         (1:pcols,1:pver+1)          , &
         cme          (1:pcols,1:pver)          , &
         tpert        (1:pcols)          , &
         dlf          (1:pcols,1:pver)          , &
         pflx         (1:pcols,1:pver+1)          , &
         zdu          (1:pcols,1:pver)          , &
         rprd         (1:pcols,1:pver)          , &!out
         mu           (1:pcols,1:pver)         , &!   REAL(r8), INTENT(out) :: mu(pcols,pver)
         md           (1:pcols,1:pver)         , &!   REAL(r8), INTENT(out) :: eu(pcols,pver)
         du           (1:pcols,1:pver)         , &!   REAL(r8), INTENT(out) :: du(pcols,pver)
         eu           (1:pcols,1:pver)         , &!   REAL(r8), INTENT(out) :: md(pcols,pver)
         ed           (1:pcols,1:pver)         , &!   REAL(r8), INTENT(out) :: ed(pcols,pver)
         dp           (1:pcols,1:pver)         , &!   REAL(r8), INTENT(out) :: dp(pcols,pver)           ! wg layer thickness in mbs (between upper/lower interface).
         dsubcld      (1:pcols)                , &! REAL(r8), INTENT(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
         jt           (1:pcols)           , &
         maxg         (1:pcols)         , &
         ideep        (1:pcols)        , &
         lengath                       , &
         ql           (1:pcols,1:pver)               , &!out
         rliq         (1:pcols) ,&
         kctop1 (1:pcols),&
         kcbot1 (1:pcols),&
         kuo (1:pcols))

    !
    ! Convert mass flux from reported mb/s to kg/m^2/s
    !
    mcon(1:pcols,1:pver)   = mcon(1:pcols,1:pver) * 100.0_r8/gravit

    ptend_loc%name(latco)  = 'zm_convr'
    ptend_loc%ls(latco)    = .TRUE.
    ptend_loc%lq(latco,1)  = .TRUE.

    ftem(1:pcols,1:pver) = ptend_loc%s(1:pcols,1:pver,latco)/cpair

    ! add tendency from this process to tendencies from other processes
    CALL physics_ptend_sum(ptend_loc,ptend_all, state,ppcnst,latco,pcols)

    ! update physics state type state1 with ptend_loc 
    CALL physics_update(state1, tend, ptend_loc, ztodt,ppcnst,pcols,pver,latco)
    ! initialize ptend for next process
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)
    !
    ! Determine the phase of the precipitation produced and add latent heat of fusion
    ! Evaporate some of the precip directly into the environment (Sundqvist)
    ! Allow this to use the updated state1 and the fresh ptend_loc type
    ! heating and specific humidity tendencies produced
    !
    ptend_loc%name(latco)  = 'zm_conv_evap'
    ptend_loc%ls(latco)    = .TRUE.
    ptend_loc%lq(latco,1)  = .TRUE.


    CALL zm_conv_evap(pcols, pver,pverp      , &
         state1%t    (1:pcols,1:pver,latco)      , &
         state1%pmid (1:pcols,1:pver,latco)      , &
         state1%pdel (1:pcols,1:pver,latco)      , &
         state1%q    (1:pcols,1:pver,latco,1)    , &
         ptend_loc%s (1:pcols,1:pver,latco)      , &
         ptend_loc%q (1:pcols,1:pver,latco,1)    , &
         rprd        (1:pcols,1:pver)            , &
         cld         (1:pcols,1:pver)            , &
         0.5_r8*ztodt                                   , &
         prec        (1:pcols )                  , &
         snow        (1:pcols)                   , &
         ntprprd     (1:pcols,1:pver  )          , &
         ntsnprd     (1:pcols,1:pver  )          , &
         flxprec     (1:pcols,1:pver+1  )        , &
         flxsnow     (1:pcols,1:pver+1  ))

    ftem(1:pcols,1:pver) = ptend_loc%s(1:pcols,1:pver,latco)/cpair

    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, state,ppcnst,latco,pcols)

    ! update physics state type state1 with ptend_loc 
    CALL physics_update(state1, tend, ptend_loc, ztodt,ppcnst,pcols,pver,latco)

    ! initialize ptend for next process
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)

    !
    ! Transport cloud water and ice only
    !
    ptend_loc%name(latco) = 'convtran1'
    !PK so para escalaes  ptend_loc%lq(latco,ixcldice) = .TRUE.
    !PK so para escalaes  ptend_loc%lq(latco,ixcldliq) = .TRUE.
    ptend_loc%lq(latco,ixcldice) = .FALSE.
    ptend_loc%lq(latco,ixcldliq) = .FALSE.

    CALL convtran (pcols                         , &
         pver                                   , &
         ptend_loc%lq(latco,:)                  , &
         state1%q(1:pcols,1:pver,latco,1:ppcnst), &
         ppcnst                                 , &
         mu          (1:pcols,1:pver)          , &
         md          (1:pcols,1:pver)          , &
         du          (1:pcols,1:pver)          , &
         eu          (1:pcols,1:pver)          , &
         ed          (1:pcols,1:pver)          , &
         dp          (1:pcols,1:pver)          , &
         dsubcld     (1:pcols)    , &
         jt          (1:pcols)    , &
         maxg        (1:pcols)    , &
         ideep       (1:pcols)    , &
         1                              , &
         lengath                        , &
         fracis      (1:pcols,1:pver,1:ppcnst)            , &
         ptend_loc%q (1:pcols,1:pver,latco,1:ppcnst) )  


    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, state,ppcnst,latco,pcols)

    ! ptend_all will be applied to original state on return to tphysbc
    ! This name triggers a special case in physics_types.F90:physics_update()
    ptend_all%name(latco) = 'convect_deep'


    call physics_update (state, tend, ptend_all, ztodt,ppcnst,pcols,pver,latco)

    !
    ! Associate pointers with physics buffer fields
    !

    DO k=1,pver
       DO i=1,pcols          
          dtdt (i,pver+1-k)=(state%t(i,k,latco         )-state_t   (i,pver+1-k))/ztodt
          dqdt (i,pver+1-k)=(state%q(i,k,latco,1       )-state_qv  (i,pver+1-k))/ztodt
          dqldt(i,pver+1-k)=(state%q(i,k,latco,ixcldliq)-state_ql  (i,pver+1-k))/ztodt
          dqidt(i,pver+1-k)=(state%q(i,k,latco,ixcldice)-state_qi  (i,pver+1-k))/ztodt
          
          state_t   (i,pver+1-k)= state%t(i,k,latco)
          state_qv  (i,pver+1-k)= state%q(i,k,latco,1)
          state_ql  (i,pver+1-k)= state%q(i,k,latco,ixcldliq)
          state_qi  (i,pver+1-k)= state%q(i,k,latco,ixcldice)          
          state_dlf (i,k)       = dlf    (i,pver+1-k)  
          state_zdu (i,k)       = zdu    (i,pver+1-k)

          state_cld (i,k)       = cld    (i,pver+1-k) 
          !state_qliq   (i,k,latco)= ql     (i,pver+1-k)
          !state_rprd   (i,k,latco)= rprd   (i,pver+1-k)
       END DO
    END DO

    DO k=1,pver+1
       DO i=1,pcols
           state_mcon (i,k) =mcon (i,pver+2-k)
       END DO
    END DO

  !  DO i=1,pcols
  !     !state_jctop  (i,latco)=jctop  (i)
  !     state_jcbot  (i,latco)=jcbot  (i)
  !  END DO

    DO ll=1,ppcnst
       DO k=1,pver
          DO i=1,pcols
             state_fracis (i,k,latco,ll)=fracis (i,pver+1-k,ll) 
          END DO
       END DO

    END DO
  END SUBROUTINE RunCu_ZhangMcFarlane

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
    REAL(r8) :: qq(pcols,pver,ppcnst)                     ! time step
    !
    !---------------------------Local storage-------------------------------
    INTEGER :: i,k,m                               ! column,level,constituent indices
    !INTEGER :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    CHARACTER*40 :: name    ! param and tracer name for qneg3
    !-----------------------------------------------------------------------
    !#if ( defined SCAM )
    !    ! The column radiation model does not update the state
    !    if(switch(CRM_SW+1)) return
    !#endif

    ! Update u,v fields
    IF(ptend%lu(latco)) THEN
       DO k = ptend%top_level(latco), ptend%bot_level(latco)
          DO i = 1, pcols
             state%u  (i,k,latco) = state%u  (i,k,latco) + ptend%u(i,k,latco) * dt
             tend%dudt(i,k,latco) = tend%dudt(i,k,latco) + ptend%u(i,k,latco)
          END DO
       END DO
    END IF

    IF(ptend%lv(latco)) THEN
       DO k = ptend%top_level(latco), ptend%bot_level(latco)
          DO i = 1, pcols
             state%v  (i,k,latco) = state%v  (i,k,latco) + ptend%v(i,k,latco) * dt
             tend%dvdt(i,k,latco) = tend%dvdt(i,k,latco) + ptend%v(i,k,latco)
          END DO
       END DO
    END IF

    ! Update dry static energy
    IF(ptend%ls(latco)) THEN
       DO k = ptend%top_level(latco), ptend%bot_level(latco)
          DO i = 1, pcols
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
             DO i = 1,pcols
                !PRINT*,state%q(i,k,latco,m) , ptend%q(i,k,latco,m) , dt
                state%q(i,k,latco,m) = state%q(i,k,latco,m) + ptend%q(i,k,latco,m) * dt
             END DO
          END DO
          ! now test for mixing ratios which are too small
          name = TRIM(ptend%name(latco)) !// '/' // trim(cnst_name(m))
          qq(1:pcols,1:pver,m)=state%q(1:pcols,1:pver,latco,m)         
          CALL qneg3(TRIM(name),  pcols, pcols, pver, m, m, qmin(m), qq(1,1,m))
          state%q(1:pcols,1:pver,latco,m)=qq(1:pcols,1:pver,m)
       END IF
    END DO

    ! special test for cloud water
    IF(ptend%lq(latco,ixcldliq)) THEN
       IF (ptend%name(latco) == 'stratiform' .OR. ptend%name(latco) == 'cldwat'  ) THEN
          IF(PERGRO)THEN
             WHERE (state%q(1:pcols,1:pver,latco,ixcldliq) < 1.0e-12_r8)
                state%q(1:pcols,1:pver,latco,ixcldliq) = 0.0_r8
             END WHERE
          ENDIF
       ELSE IF (ptend%name(latco) == 'convect_deep') THEN
          WHERE (state%q(1:pcols,1:pver,latco,ixcldliq) < 1.0e-36_r8)
             state%q(1:pcols,1:pver,latco,ixcldliq) = 0.0_r8
          END WHERE
       END IF
    END IF
    IF(ptend%lq(latco,ixcldice)) THEN
       IF (ptend%name(latco) == 'stratiform' .OR. ptend%name(latco) == 'cldwat'  ) THEN
          IF( PERGRO)THEN
             WHERE (state%q(1:pcols,1:pver,latco,ixcldice) < 1.0e-12_r8)
                state%q(1:pcols,1:pver,latco,ixcldice) = 0.0_r8
             END WHERE
          ENDIF
       ELSE IF (ptend%name(latco) == 'convect_deep') THEN
          WHERE (state%q(1:pcols,1:pver,latco,ixcldice) < 1.0e-36_r8)
             state%q(1:pcols,1:pver,latco,ixcldice) = 0.0_r8
          END WHERE
       END IF
    END IF

    ! Derive new temperature and geopotential fields if heating or water tendency not 0.
    IF (ptend%ls(latco) .OR. ptend%lq(latco,1)) THEN
       CALL geopotential_dse(                                                                    &
            state%lnpint(1:pcols,1:pver+1,latco), state%lnpmid(1:pcols,1:pver,latco) , state%pint (1:pcols,1:pver+1,latco)  , &
            state%pmid  (1:pcols,1:pver,latco), state%pdel  (1:pcols,1:pver,latco) , state%rpdel(1:pcols,1:pver,latco)  , &
            state%s     (1:pcols,1:pver,latco), state%q     (1:pcols,1:pver,latco,1), state%phis (1:pcols,latco) , rair  , &
            gravit      , cpair        ,zvir        , &
            state%t(1:pcols,1:pver,latco)     , state%zi(1:pcols,1:pver+1,latco)    , state%zm(1:pcols,1:pver,latco)       ,&
            pcols, pver, pver+1          )

    END IF

    ! Reset all parameterization tendency flags to false
    CALL physics_ptend_reset(ptend,latco,ppcnst,pver)

  END SUBROUTINE physics_update




  !===============================================================================
  SUBROUTINE physics_ptend_sum(ptend, ptend_sum, state,ppcnst,latco,pcols)
    !-----------------------------------------------------------------------
    ! Add ptend fields to ptend_sum for ptend logical flags = .true.
    ! Where ptend logical flags = .false, don't change ptend_sum
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN   ) ::ppcnst,latco,pcols
    TYPE(physics_ptend), INTENT(in)     :: ptend   ! New parameterization tendencies
    TYPE(physics_ptend), INTENT(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    TYPE(physics_state), INTENT(in)     :: state   ! New parameterization tendencies

    !---------------------------Local storage-------------------------------
    INTEGER :: i,k,m                               ! column,level,constituent indices

    !-----------------------------------------------------------------------


    ! Update u,v fields
    IF(ptend%lu(latco)) THEN
       ptend_sum%lu(latco) = .TRUE.
       DO i = 1, pcols
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             ptend_sum%u(i,k,latco) = ptend_sum%u(i,k,latco) + ptend%u(i,k,latco)
          END DO
          ptend_sum%taux_srf(i,latco) = ptend_sum%taux_srf(i,latco) + ptend%taux_srf(i,latco)
          ptend_sum%taux_top(i,latco) = ptend_sum%taux_top(i,latco) + ptend%taux_top(i,latco)
       END DO
    END IF

    IF(ptend%lv(latco)) THEN
       ptend_sum%lv(latco) = .TRUE.
       DO i = 1, pcols
          DO k = ptend%top_level(latco), ptend%bot_level(latco)
             ptend_sum%v(i,k,latco) = ptend_sum%v(i,k,latco) + ptend%v(i,k,latco)
          END DO
          ptend_sum%tauy_srf(i,latco) = ptend_sum%tauy_srf(i,latco) + ptend%tauy_srf(i,latco)
          ptend_sum%tauy_top(i,latco) = ptend_sum%tauy_top(i,latco) + ptend%tauy_top(i,latco)
       END DO
    END IF


    IF(ptend%ls(latco)) THEN
       ptend_sum%ls(latco) = .TRUE.
       DO i = 1, pcols
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
          DO i = 1,pcols
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
!    TYPE(physics_state), INTENT(in) :: state
!    TYPE(physics_state), INTENT(out) :: state1

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

  SUBROUTINE zm_convr(pcols,pver, pcnst,pnats,pverp,&
       t       ,qh      ,prec    ,jctop   ,jcbot   , &
       pblh    ,zm      ,geos    ,zi      ,qtnd    , &
       heat    ,pap     ,paph    ,dpp     , &
       delt    ,mcon    ,cme     ,          &
       tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
       mu      ,md      ,du      ,eu      ,ed      , &
       dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
       lengath ,ql      ,rliq    ,kctop1   ,kcbot1,kuo)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Main driver for zhang-mcfarlane convection scheme 
    ! 
    ! Method: 
    ! performs deep convective adjustment based on mass-flux closure
    ! algorithm.
    ! 
    ! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
    !
    ! This is contributed code not fully standardized by the CAM core group.
    ! All variables have been typed, where most are identified in comments
    ! The current procedure will be reimplemented in a subsequent version
    ! of the CAM where it will include a more straightforward formulation
    ! and will make use of the standard CAM nomenclature
    ! 
    !-----------------------------------------------------------------------
    !   use constituents, only: pcnst, pnats

    !
    ! ************************ index of variables **********************
    !
    !  wg * alpha    array of vertical differencing used (=1. for upstream).
    !  w  * cape     convective available potential energy.
    !  wg * capeg    gathered convective available potential energy.
    !  c  * capelmt  threshold value for cape for deep convection.
    !  ic  * cpres    specific heat at constant pressure in j/kg-degk.
    !  i  * dpp      
    !  ic  * delt     length of model time-step in seconds.
    !  wg * dp       layer thickness in mbs (between upper/lower interface).
    !  wg * dqdt     mixing ratio tendency at gathered points.
    !  wg * dsdt     dry static energy ("temp") tendency at gathered points.
    !  wg * dudt     u-wind tendency at gathered points.
    !  wg * dvdt     v-wind tendency at gathered points.
    !  wg * dsubcld  layer thickness in mbs between lcl and maxi.
    !  ic  * grav     acceleration due to gravity in m/sec2.
    !  wg * du       detrainment in updraft. specified in mid-layer
    !  wg * ed       entrainment in downdraft.
    !  wg * eu       entrainment in updraft.
    !  wg * hmn      moist static energy.
    !  wg * hsat     saturated moist static energy.
    !  w  * ideep    holds position of gathered points vs longitude index.
    !  ic  * pver     number of model levels.
    !  wg * j0       detrainment initiation level index.
    !  wg * jd       downdraft   initiation level index.
    !  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
    !  wg * jt       top  level index of deep cumulus convection.
    !  w  * lcl      base level index of deep cumulus convection.
    !  wg * lclg     gathered values of lcl.
    !  w  * lel      index of highest theoretical convective plume.
    !  wg * lelg     gathered values of lel.
    !  w  * lon      index of onset level for deep convection.
    !  w  * maxi     index of level with largest moist static energy.
    !  wg * maxg     gathered values of maxi.
    !  wg * mb       cloud base mass flux.
    !  wg * mc       net upward (scaled by mb) cloud mass flux.
    !  wg * md       downward cloud mass flux (positive up).
    !  wg * mu       upward   cloud mass flux (positive up). specified
    !                at interface
    !  ic  * msg      number of missing moisture levels at the top of model.
    !  w  * p        grid slice of ambient mid-layer pressure in mbs.
    !  i  * pblt     row of pbl top indices.
    !  w  * pcpdh    scaled surface pressure.
    !  w  * pf       grid slice of ambient interface pressure in mbs.
    !  wg * pg       grid slice of gathered values of p.
    !  w  * q        grid slice of mixing ratio.
    !  wg * qd       grid slice of mixing ratio in downdraft.
    !  wg * qg       grid slice of gathered values of q.
    !  i/o * qh       grid slice of specific humidity.
    !  w  * qh0      grid slice of initial specific humidity.
    !  wg * qhat     grid slice of upper interface mixing ratio.
    !  wg * ql       grid slice of cloud liquid water.
    !  wg * qs       grid slice of saturation mixing ratio.
    !  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
    !  wg * qstpg    grid slice of gathered values of qstp.
    !  wg * qu       grid slice of mixing ratio in updraft.
    !  ic  * rgas     dry air gas constant.
    !  wg * rl       latent heat of vaporization.
    !  w  * s        grid slice of scaled dry static energy (t+gz/cp).
    !  wg * sd       grid slice of dry static energy in downdraft.
    !  wg * sg       grid slice of gathered values of s.
    !  wg * shat     grid slice of upper interface dry static energy.
    !  wg * su       grid slice of dry static energy in updraft.
    !  i/o * t       
    !  o  * jctop    row of top-of-deep-convection indices passed out.
    !  O  * jcbot    row of base of cloud indices passed out.
    !  wg * tg       grid slice of gathered values of t.
    !  w  * tl       row of parcel temperature at lcl.
    !  wg * tlg      grid slice of gathered values of tl.
    !  w  * tp       grid slice of parcel temperatures.
    !  wg * tpg      grid slice of gathered values of tp.
    !  i/o * u        grid slice of u-wind (real).
    !  wg * ug       grid slice of gathered values of u.
    !  i/o * utg      grid slice of u-wind tendency (real).
    !  i/o * v        grid slice of v-wind (real).
    !  w  * va       work array re-used by called subroutines.
    !  wg * vg       grid slice of gathered values of v.
    !  i/o * vtg      grid slice of v-wind tendency (real).
    !  i  * w        grid slice of diagnosed large-scale vertical velocity.
    !  w  * z        grid slice of ambient mid-layer height in metres.
    !  w  * zf       grid slice of ambient interface height in metres.
    !  wg * zfg      grid slice of gathered values of zf.
    !  wg * zg       grid slice of gathered values of z.
    !
    !-----------------------------------------------------------------------
    !
    ! multi-level i/o fields:
    !  i      => input arrays.
    !  i/o    => input/output arrays.
    !  w      => work arrays.
    !  wg     => work arrays operating only on gathered points.
    !  ic     => input data constants.
    !  c      => data constants pertaining to subroutine itself.
    !
    ! input arguments
    !
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: pverp
    INTEGER, INTENT(in) :: pcnst     ! number of advected constituents (including water vapor)
    INTEGER, INTENT(in) :: pnats     ! number of non-advected constituents
    !integer, parameter, public :: ppcnst = pcnst+pnats! total number of constituents

    REAL(r8), INTENT(in) :: t(pcols,pver)              ! grid slice of temperature at mid-layer.
    REAL(r8), INTENT(in) :: qh(pcols,pver,pcnst+pnats) ! grid slice of specific humidity.
    REAL(r8), INTENT(in) :: pap(pcols,pver)     
    REAL(r8), INTENT(in) :: paph(pcols,pver+1)
    REAL(r8), INTENT(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
    REAL(r8), INTENT(in) :: zm(pcols,pver)
    REAL(r8), INTENT(in) :: geos(pcols)
    REAL(r8), INTENT(in) :: zi(pcols,pver+1)
    REAL(r8), INTENT(in) :: pblh(pcols)
    REAL(r8), INTENT(in) :: tpert(pcols)
    !
    ! output arguments
    !
    REAL(r8), INTENT(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
    REAL(r8), INTENT(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
    REAL(r8), INTENT(out) :: mcon(pcols,pverp)
    REAL(r8), INTENT(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
    REAL(r8), INTENT(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
    REAL(r8), INTENT(out) :: cme(pcols,pver)
    REAL(r8), INTENT(out) :: zdu(pcols,pver)
    REAL(r8), INTENT(out) :: rprd(pcols,pver)     ! rain production rate
    ! move these vars from local storage to output so that convective
    ! transports can be done in outside of conv_cam.
    REAL(r8), INTENT(out) :: mu(pcols,pver)
    REAL(r8), INTENT(out) :: eu(pcols,pver)
    REAL(r8), INTENT(out) :: du(pcols,pver)
    REAL(r8), INTENT(out) :: md(pcols,pver)
    REAL(r8), INTENT(out) :: ed(pcols,pver)
    REAL(r8), INTENT(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
    REAL(r8), INTENT(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
    REAL(r8), INTENT(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
    REAL(r8), INTENT(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
    REAL(r8), INTENT(out) :: prec(pcols)
    REAL(r8), INTENT(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
    INTEGER, INTENT(out) :: kctop1(pcols)
    INTEGER, INTENT(out) :: kcbot1(pcols)
    INTEGER, INTENT(out) :: kuo   (pcols)
    
    REAL(r8) zs(pcols)
    REAL(r8) dlg(pcols,pver)    ! gathrd version of the detraining cld h2o tend
    REAL(r8) pflxg(pcols,pverp) ! gather precip flux at each level
    REAL(r8) cug(pcols,pver)    ! gathered condensation rate
    REAL(r8) evpg(pcols,pver)   ! gathered evap rate of rain in downdraft
    REAL(r8) mumax(pcols)
    INTEGER jt(pcols)                          ! wg top  level index of deep cumulus convection.
    INTEGER maxg(pcols)                        ! wg gathered values of maxi.
    INTEGER ideep(pcols)                       ! w holds position of gathered points vs longitude index.
    INTEGER lengath
    !     diagnostic field used by chem/wetdep codes
    REAL(r8) ql(pcols,pver)                    ! wg grid slice of cloud liquid water.
    !
    REAL(r8) pblt(pcols)           ! i row of pbl top indices.
    !
    !-----------------------------------------------------------------------
    !
    ! general work fields (local variables):
    !
    REAL(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
    REAL(r8) p(pcols,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
    REAL(r8) z(pcols,pver)              ! w  grid slice of ambient mid-layer height in metres.
    REAL(r8) s(pcols,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
    REAL(r8) tp(pcols,pver)             ! w  grid slice of parcel temperatures.
    REAL(r8) zf(pcols,pver+1)           ! w  grid slice of ambient interface height in metres.
    REAL(r8) pf(pcols,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
    REAL(r8) qstp(pcols,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

    REAL(r8) cape(pcols)                ! w  convective available potential energy.
    REAL(r8) tl(pcols)                  ! w  row of parcel temperature at lcl.

    INTEGER lcl(pcols)                  ! w  base level index of deep cumulus convection.
    INTEGER lel(pcols)                  ! w  index of highest theoretical convective plume.
    INTEGER lon(pcols)                  ! w  index of onset level for deep convection.
    INTEGER maxi(pcols)                 ! w  index of level with largest moist static energy.
    INTEGER INDEX(pcols)
    REAL(r8) precip
    !
    ! gathered work fields:
    !
    REAL(r8) qg(pcols,pver)             ! wg grid slice of gathered values of q.
    REAL(r8) tg(pcols,pver)             ! w  grid slice of temperature at interface.
    REAL(r8) pg(pcols,pver)             ! wg grid slice of gathered values of p.
    REAL(r8) zg(pcols,pver)             ! wg grid slice of gathered values of z.
    REAL(r8) sg(pcols,pver)             ! wg grid slice of gathered values of s.
    REAL(r8) tpg(pcols,pver)            ! wg grid slice of gathered values of tp.
    REAL(r8) zfg(pcols,pver+1)          ! wg grid slice of gathered values of zf.
    REAL(r8) qstpg(pcols,pver)          ! wg grid slice of gathered values of qstp.
    REAL(r8) ug(pcols,pver)             ! wg grid slice of gathered values of u.
    REAL(r8) vg(pcols,pver)             ! wg grid slice of gathered values of v.
    REAL(r8) cmeg(pcols,pver)

    REAL(r8) rprdg(pcols,pver)           ! wg gathered rain production rate
    REAL(r8) capeg(pcols)               ! wg gathered convective available potential energy.
    REAL(r8) tlg(pcols)                 ! wg grid slice of gathered values of tl.
    INTEGER lclg(pcols)                 ! wg gathered values of lcl.
    INTEGER lelg(pcols)
    !
    ! work fields arising from gathered calculations.
    !
    REAL(r8) dqdt(pcols,pver)           ! wg mixing ratio tendency at gathered points.
    REAL(r8) dsdt(pcols,pver)           ! wg dry static energy ("temp") tendency at gathered points.
    !      real(r8) alpha(pcols,pver)      ! array of vertical differencing used (=1. for upstream).
    REAL(r8) sd(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
    REAL(r8) qd(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
    REAL(r8) mc(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
    REAL(r8) qhat(pcols,pver)           ! wg grid slice of upper interface mixing ratio.
    REAL(r8) qu(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
    REAL(r8) su(pcols,pver)             ! wg grid slice of dry static energy in updraft.
    REAL(r8) qs(pcols,pver)             ! wg grid slice of saturation mixing ratio.
    REAL(r8) shat(pcols,pver)           ! wg grid slice of upper interface dry static energy.
    REAL(r8) hmn(pcols,pver)            ! wg moist static energy.
    REAL(r8) hsat(pcols,pver)           ! wg saturated moist static energy.
    REAL(r8) qlg(pcols,pver)
    REAL(r8) dudt(pcols,pver)           ! wg u-wind tendency at gathered points.
    REAL(r8) dvdt(pcols,pver)           ! wg v-wind tendency at gathered points.
    !      real(r8) ud(pcols,pver)
    !      real(r8) vd(pcols,pver)

    REAL(r8) mb(pcols)                  ! wg cloud base mass flux.

    INTEGER jlcl(pcols)
    INTEGER j0(pcols)                 ! wg detrainment initiation level index.
    INTEGER jd(pcols)                 ! wg downdraft initiation level index.

    REAL(r8) delt                     ! length of model time-step in seconds.

    INTEGER i
    INTEGER ii
    INTEGER k
    INTEGER msg                      !  ic number of missing moisture levels at the top of model.
    REAL(r8) qdifr
    REAL(r8) sdifr
    !
    !--------------------------Data statements------------------------------
    !
    ! Set internal variable "msg" (convection limit) to "limcnv-1"
    !
    msg = limcnv - 1
    !
    ! initialize necessary arrays.
    ! zero out variables not used in cam
    !
    qtnd(:,:) = 0.0_r8
    heat(:,:) = 0.0_r8
    mcon(:,:) = 0.0_r8
    rliq(:pcols)   = 0.0_r8
    !
    ! initialize convective tendencies
    !
    prec(:pcols) = 0.0_r8
    DO k = 1,pver
       DO i = 1,pcols
          dqdt(i,k)  = 0.0_r8
          dsdt(i,k)  = 0.0_r8
          dudt(i,k)  = 0.0_r8
          dvdt(i,k)  = 0.0_r8
          pflx(i,k)  = 0.0_r8
          pflxg(i,k) = 0.0_r8
          cme(i,k)   = 0.0_r8
          rprd(i,k)  = 0.0_r8
          zdu(i,k)   = 0.0_r8
          ql(i,k)    = 0.0_r8
          qlg(i,k)   = 0.0_r8
          dlf(i,k)   = 0.0_r8
          dlg(i,k)   = 0.0_r8
       END DO
    END DO
    DO i = 1,pcols
       pflx(i,pverp) = 0.0_r8
       pflxg(i,pverp) = 0.0_r8
    END DO
    !
    DO i = 1,pcols
       pblt(i) = pver
       dsubcld(i) = 0.0_r8
       kctop1(i)=1
       kcbot1(i)=1
       kuo  (i)=0
       jctop(i) = pver
       jcbot(i) = 1
    END DO
    !
    ! calculate local pressure (mbs) and height (m) for both interface
    ! and mid-layer locations.
    !
    DO i = 1,pcols
       zs(i) = geos(i)*rgrav
       pf(i,pver+1) = paph(i,pver+1)*0.01_r8
       zf(i,pver+1) = zi(i,pver+1) + zs(i)
    END DO
    DO k = 1,pver
       DO i = 1,pcols
          p(i,k) = pap(i,k)*0.01_r8
          pf(i,k) = paph(i,k)*0.01_r8
          z(i,k) = zm(i,k) + zs(i)
          zf(i,k) = zi(i,k) + zs(i)
       END DO
    END DO
    !
    DO k = pver - 1,msg + 1,-1
       DO i = 1,pcols
          IF (ABS(z(i,k)-zs(i)-pblh(i)) < (zf(i,k)-zf(i,k+1))*0.5_r8) pblt(i) = k
       END DO
    END DO
    !
    ! store incoming specific humidity field for subsequent calculation
    ! of precipitation (through change in storage).
    ! define dry static energy (normalized by cp).
    !
    DO k = 1,pver
       DO i = 1,pcols
          q(i,k) = qh(i,k,1)
          s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
          tp(i,k)=0.0_r8
          shat(i,k) = s(i,k)
          qhat(i,k) = q(i,k)
       END DO
    END DO
    DO i = 1,pcols
       capeg(i) = 0.0_r8
       lclg(i) = 1
       lelg(i) = pver
       maxg(i) = 1
       tlg(i) = 400.0_r8
       dsubcld(i) = 0.0_r8
    END DO
    !
    ! evaluate covective available potential energy (cape).
    !
    !CALL buoyan( pcols,pver,&
    !     q       ,t       ,p       ,z       ,pf       , &
    !     tp      ,qstp    ,tl      ,rl      ,cape     , &
    !     pblt    ,lcl     ,lel     ,lon     ,maxi     , &
    !     rgas    ,grav    ,cpres   ,msg     , &
    !     tpert   )
    call buoyan_dilute(1        ,pcols,pcols,pver    , &
               q       ,t       ,p       ,z       ,pf       , &
               tp      ,qstp    ,tl      ,rl      ,cape     , &
               pblt    ,lcl     ,lel     ,lon     ,maxi     , &
               rgas    ,grav    ,cpres   ,msg     , &
               tpert   )

    !
    ! determine whether grid points will undergo some deep convection
    ! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
    ! (require cape.gt. 0 and lel<lcl as minimum conditions).
    !
    lengath = 0
    DO i=1,pcols
       IF (cape(i) > capelmt) THEN
          lengath = lengath + 1
          INDEX(lengath) = i
       END IF
    END DO

    IF (lengath.EQ.0) RETURN
    DO ii=1,lengath
       i=INDEX(ii)
       ideep(ii)=i
    END DO
    !
    ! obtain gathered arrays necessary for ensuing calculations.
    !
    DO k = 1,pver
       DO i = 1,lengath
          dp(i,k) = 0.01_r8*dpp(ideep(i),k)
          qg(i,k) = q(ideep(i),k)
          tg(i,k) = t(ideep(i),k)
          pg(i,k) = p(ideep(i),k)
          zg(i,k) = z(ideep(i),k)
          sg(i,k) = s(ideep(i),k)
          tpg(i,k) = tp(ideep(i),k)
          zfg(i,k) = zf(ideep(i),k)
          qstpg(i,k) = qstp(ideep(i),k)
          ug(i,k) = 0.0_r8
          vg(i,k) = 0.0_r8
       END DO
    END DO
    !
    DO i = 1,lengath
       zfg(i,pver+1) = zf(ideep(i),pver+1)
    END DO
    DO i = 1,lengath
       capeg(i) = cape(ideep(i))
       lclg(i) = lcl(ideep(i))
       lelg(i) = lel(ideep(i))
       maxg(i) = maxi(ideep(i))
       tlg(i) = tl(ideep(i))
    END DO
    !
    ! calculate sub-cloud layer pressure "thickness" for use in
    ! closure and tendency routines.
    !
    DO k = msg + 1,pver
       DO i = 1,lengath
          IF (k >= maxg(i)) THEN
             dsubcld(i) = dsubcld(i) + dp(i,k)
          END IF
       END DO
    END DO
    !
    ! define array of factors (alpha) which defines interfacial
    ! values, as well as interfacial values for (q,s) used in
    ! subsequent routines.
    !
    DO k = msg + 2,pver
       DO i = 1,lengath
          !            alpha(i,k) = 0.5
          sdifr = 0.0_r8
          qdifr = 0.0_r8
          IF (sg(i,k) > 0.0_r8 .OR. sg(i,k-1) > 0.0_r8) &
               sdifr = ABS((sg(i,k)-sg(i,k-1))/MAX(sg(i,k-1),sg(i,k)))
          IF (qg(i,k) > 0.0_r8 .OR. qg(i,k-1) > 0.0_r8) &
               qdifr = ABS((qg(i,k)-qg(i,k-1))/MAX(qg(i,k-1),qg(i,k)))
          IF (sdifr > 1.0E-6_r8) THEN
             shat(i,k) = LOG(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
          ELSE
             shat(i,k) = 0.5_r8* (sg(i,k)+sg(i,k-1))
          END IF
          IF (qdifr > 1.0E-6_r8) THEN
             qhat(i,k) = LOG(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
          ELSE
             qhat(i,k) = 0.5_r8* (qg(i,k)+qg(i,k-1))
          END IF
       END DO
    END DO
    !
    ! obtain cloud properties.
    !
    CALL cldprp(pcols,pver,pverp ,&
         qg      ,tg      ,ug      ,vg      ,pg      , &
         zg      ,sg      ,mu      ,eu      ,du      , &
         md      ,ed      ,sd      ,qd      ,mc      , &
         qu      ,su      ,zfg     ,qs      ,hmn     , &
         hsat    ,shat    ,qlg     , &
         cmeg    ,maxg    ,lelg    ,jt      ,jlcl    , &
         maxg    ,j0      ,jd      ,rl      ,lengath , &
         rgas    ,grav    ,cpres   ,msg     , &
         pflxg   ,evpg    ,cug     ,rprdg   ,limcnv  )
    !
    ! convert detrainment from units of "1/m" to "1/mb".
    !
    DO k = msg + 1,pver
       DO i = 1,lengath
          du   (i,k) = du   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
          eu   (i,k) = eu   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
          ed   (i,k) = ed   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
          cug  (i,k) = cug  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
          cmeg (i,k) = cmeg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
          rprdg(i,k) = rprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
          evpg (i,k) = evpg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
       END DO
    END DO

    CALL closure(pcols,pver, &
         qg      ,tg      ,pg      ,zg      ,sg      , &
         tpg     ,qs      ,qu      ,su      ,mc      , &
         du      ,mu      ,md      ,qd      ,sd      , &
         qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
         qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
         lclg    ,lelg    ,jt      ,maxg    ,1       , &
         lengath ,rgas    ,grav    ,cpres   ,rl      , &
         msg     ,capelmt    )
    !
    ! limit cloud base mass flux to theoretical upper bound.
    !
    DO i=1,lengath
       mumax(i) = 0
    END DO
    DO k=msg + 2,pver
       DO i=1,lengath
          mumax(i) = MAX(mumax(i), mu(i,k)/dp(i,k))
       END DO
    END DO

    DO i=1,lengath
       IF (mumax(i) > 0.0_r8) THEN
          mb(i) = MIN(mb(i),0.5_r8/(delt*mumax(i)))
       ELSE
          mb(i) = 0.0_r8
       ENDIF
    END DO
    ! If no_deep_pbl = .true., don't allow convection entirely 
    ! within PBL (suggestion of Bjorn Stevens, 8-2000)

   if (no_deep_pbl) then
      do i=1,lengath
         if (zm(ideep(i),jt(i)) < pblh(ideep(i))) mb(i) = 0
      end do
   end if


    DO k=msg+1,pver
       DO i=1,lengath
          mu   (i,k)  = mu   (i,k)*mb(i)
          md   (i,k)  = md   (i,k)*mb(i)
          mc   (i,k)  = mc   (i,k)*mb(i)
          du   (i,k)  = du   (i,k)*mb(i)
          eu   (i,k)  = eu   (i,k)*mb(i)
          ed   (i,k)  = ed   (i,k)*mb(i)
          cmeg (i,k)  = cmeg (i,k)*mb(i)
          rprdg(i,k)  = rprdg(i,k)*mb(i)
          cug  (i,k)  = cug  (i,k)*mb(i)
          evpg (i,k)  = evpg (i,k)*mb(i)
          pflxg(i,k+1)= pflxg(i,k+1)*mb(i)*100.0_r8/grav
       END DO
    END DO
    !
    ! compute temperature and moisture changes due to convection.
    !
    CALL q1q2_pjr(pcols,pver, &
         dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
         su      ,du      ,qhat    ,shat    ,dp      , &
         mu      ,md      ,sd      ,qd      ,qlg     , &
         dsubcld ,jt      ,maxg    ,1       ,lengath , &
         cpres   ,rl      ,msg     ,          &
         dlg     ,evpg    ,cug     )
    !
    ! gather back temperature and mixing ratio.
    !
    DO k = msg + 1,pver
       !DIR$ CONCURRENT
       DO i = 1,lengath
          !
          ! q is updated to compute net precip.
          !
          q(ideep(i),k) = qh(ideep(i),k,1) + 2.0_r8*delt*dqdt(i,k)
          qtnd(ideep(i),k) = dqdt (i,k)
          cme (ideep(i),k) = cmeg (i,k)
          rprd(ideep(i),k) = rprdg(i,k)
          zdu (ideep(i),k) = du   (i,k)
          mcon(ideep(i),k) = mc   (i,k)
          heat(ideep(i),k) = dsdt (i,k)*cpres
          dlf (ideep(i),k) = dlg  (i,k)
          pflx(ideep(i),k) = pflxg(i,k)
          ql  (ideep(i),k) = qlg  (i,k)
       END DO
    END DO
    !
    !DIR$ CONCURRENT
    DO i = 1,lengath
       jctop(ideep(i)) = jt(i)
       !++bee
       jcbot(ideep(i)) = maxg(i)
       !--bee
       kctop1(ideep(i))=pver-jctop(ideep(i))
       kcbot1(ideep(i))=pver-jcbot(ideep(i))
       kuo  (ideep(i))=1

       !--bee
       pflx(ideep(i),pverp) = pflxg(i,pverp)
    END DO

    ! Compute precip by integrating change in water vapor minus detrained cloud water
    DO k = pver,msg + 1,-1
       DO i = 1,pcols
          prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k,1)) - dpp(i,k)*dlf(i,k)*2*delt
       END DO
    END DO

    ! obtain final precipitation rate in m/s.
    DO i = 1,pcols
       prec(i) = rgrav*MAX(prec(i),0.0_r8)/ (2.0_r8*delt)/1000.0_r8
    END DO

    ! Compute reserved liquid (not yet in cldliq) for energy integrals.
    ! Treat rliq as flux out bottom, to be added back later.
    DO k = 1, pver
       DO i = 1, pcols
          rliq(i) = rliq(i) + dlf(i,k)*dpp(i,k)/gravit
       END DO
    END DO
    rliq(:pcols) = rliq(:pcols) /1000.0_r8

    RETURN
  END SUBROUTINE zm_convr

  SUBROUTINE q1q2_pjr(pcols,pver, &
       dqdt    ,dsdt    ,q       ,qs      ,qu      , &
       su      ,du      ,qhat    ,shat    ,dp      , &
       mu      ,md      ,sd      ,qd      ,ql      , &
       dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
       cp      ,rl      ,msg     ,          &
       dl      ,evp     ,cu      )


    IMPLICIT NONE

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: phil rasch dec 19 1995
    ! 
    !-----------------------------------------------------------------------


    REAL(r8), INTENT(in) :: cp
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    INTEGER, INTENT(in) :: il1g
    INTEGER, INTENT(in) :: il2g
    INTEGER, INTENT(in) :: msg

    REAL(r8), INTENT(in) :: q(pcols,pver)
    REAL(r8), INTENT(in) :: qs(pcols,pver)
    REAL(r8), INTENT(in) :: qu(pcols,pver)
    REAL(r8), INTENT(in) :: su(pcols,pver)
    REAL(r8), INTENT(in) :: du(pcols,pver)
    REAL(r8), INTENT(in) :: qhat(pcols,pver)
    REAL(r8), INTENT(in) :: shat(pcols,pver)
    REAL(r8), INTENT(in) :: dp(pcols,pver)
    REAL(r8), INTENT(in) :: mu(pcols,pver)
    REAL(r8), INTENT(in) :: md(pcols,pver)
    REAL(r8), INTENT(in) :: sd(pcols,pver)
    REAL(r8), INTENT(in) :: qd(pcols,pver)
    REAL(r8), INTENT(in) :: ql(pcols,pver)
    REAL(r8), INTENT(in) :: evp(pcols,pver)
    REAL(r8), INTENT(in) :: cu(pcols,pver)
    REAL(r8), INTENT(in) :: dsubcld(pcols)

    REAL(r8),INTENT(out) :: dqdt(pcols,pver),dsdt(pcols,pver)
    REAL(r8),INTENT(out) :: dl(pcols,pver)
    INTEGER kbm
    INTEGER ktm
    INTEGER jt(pcols)
    INTEGER mx(pcols)
    !
    ! work fields:
    !
    INTEGER i
    INTEGER k

    REAL(r8) emc
    REAL(r8) rl
    !-------------------------------------------------------------------
    DO k = msg + 1,pver
       DO i = il1g,il2g
          dsdt(i,k) = 0.0_r8
          dqdt(i,k) = 0.0_r8
          dl(i,k) = 0.0_r8
       END DO
    END DO
    !
    ! find the highest level top and bottom levels of convection
    !
    ktm = pver
    kbm = pver
    DO i = il1g, il2g
       ktm = MIN(ktm,jt(i))
       kbm = MIN(kbm,mx(i))
    END DO

    DO k = ktm,pver-1
       DO i = il1g,il2g
          emc = -cu (i,k)               &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft

          dsdt(i,k) = -rl/cp*emc &
               + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
               -mu(i,k)*   (su(i,k)-shat(i,k)) &
               +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
               -md(i,k)*   (sd(i,k)-shat(i,k)) &
               )/dp(i,k)

          dqdt(i,k) = emc + &
               (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
               -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
               +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
               -md(i,k)*   (qd(i,k)-qhat(i,k)) &
               )/dp(i,k)

          dl(i,k) = du(i,k)*ql(i,k+1)

       END DO
    END DO

    !
    !DIR$ NOINTERCHANGE!
    DO k = kbm,pver
       DO i = il1g,il2g
          IF (k == mx(i)) THEN
             dsdt(i,k) = (1.0_r8/dsubcld(i))* &
                  (-mu(i,k)* (su(i,k)-shat(i,k)) &
                  -md(i,k)* (sd(i,k)-shat(i,k)) &
                  )
             dqdt(i,k) = (1.0_r8/dsubcld(i))* &
                  (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                  -md(i,k)*(qd(i,k)-qhat(i,k)) &
                  )
          ELSE IF (k > mx(i)) THEN
             dsdt(i,k) = dsdt(i,k-1)
             dqdt(i,k) = dqdt(i,k-1)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE q1q2_pjr



  SUBROUTINE closure(pcols,pver, &
       q       ,t       ,p       ,z       ,s       , &
       tp      ,qs      ,qu      ,su      ,mc      , &
       du      ,mu      ,md      ,qd      ,sd      , &
       qhat    ,shat    ,dp      ,qstp    ,zf      , &
       ql      ,dsubcld ,mb      ,cape    ,tl      , &
       lcl     ,lel     ,jt      ,mx      ,il1g    , &
       il2g    ,rd      ,grav    ,cp      ,rl      , &
       msg     ,capelmt )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: G. Zhang and collaborators. CCM contact:P. Rasch
    ! This is contributed code not fully standardized by the CCM core group.
    !
    ! this code is very much rougher than virtually anything else in the CCM
    ! We expect to release cleaner code in a future release
    !
    ! the documentation has been enhanced to the degree that we are able
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !
    !-----------------------------Arguments---------------------------------
    !
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(inout) :: q(pcols,pver)        ! spec humidity
    REAL(r8), INTENT(inout) :: t(pcols,pver)        ! temperature
    REAL(r8), INTENT(inout) :: p(pcols,pver)        ! pressure (mb)
    REAL(r8), INTENT(inout) :: mb(pcols)            ! cloud base mass flux
    REAL(r8), INTENT(in) :: z(pcols,pver)        ! height (m)
    REAL(r8), INTENT(in) :: s(pcols,pver)        ! normalized dry static energy
    REAL(r8), INTENT(in) :: tp(pcols,pver)       ! parcel temp
    REAL(r8), INTENT(in) :: qs(pcols,pver)       ! sat spec humidity
    REAL(r8), INTENT(in) :: qu(pcols,pver)       ! updraft spec. humidity
    REAL(r8), INTENT(in) :: su(pcols,pver)       ! normalized dry stat energy of updraft
    REAL(r8), INTENT(in) :: mc(pcols,pver)       ! net convective mass flux
    REAL(r8), INTENT(in) :: du(pcols,pver)       ! detrainment from updraft
    REAL(r8), INTENT(in) :: mu(pcols,pver)       ! mass flux of updraft
    REAL(r8), INTENT(in) :: md(pcols,pver)       ! mass flux of downdraft
    REAL(r8), INTENT(in) :: qd(pcols,pver)       ! spec. humidity of downdraft
    REAL(r8), INTENT(in) :: sd(pcols,pver)       ! dry static energy of downdraft
    REAL(r8), INTENT(in) :: qhat(pcols,pver)     ! environment spec humidity at interfaces
    REAL(r8), INTENT(in) :: shat(pcols,pver)     ! env. normalized dry static energy at intrfcs
    REAL(r8), INTENT(in) :: dp(pcols,pver)       ! pressure thickness of layers
    REAL(r8), INTENT(in) :: qstp(pcols,pver)     ! spec humidity of parcel
    REAL(r8), INTENT(in) :: zf(pcols,pver+1)     ! height of interface levels
    REAL(r8), INTENT(in) :: ql(pcols,pver)       ! liquid water mixing ratio

    REAL(r8), INTENT(in) :: cape(pcols)          ! available pot. energy of column
    REAL(r8), INTENT(in) :: tl(pcols)
    REAL(r8), INTENT(in) :: dsubcld(pcols)       ! thickness of subcloud layer

    INTEGER, INTENT(in) :: lcl(pcols)        ! index of lcl
    INTEGER, INTENT(in) :: lel(pcols)        ! index of launch leve
    INTEGER, INTENT(in) :: jt(pcols)         ! top of updraft
    INTEGER, INTENT(in) :: mx(pcols)         ! base of updraft
    !
    !--------------------------Local variables------------------------------
    !
    REAL(r8) dtpdt(pcols,pver)
    REAL(r8) dqsdtp(pcols,pver)
    REAL(r8) dtmdt(pcols,pver)
    REAL(r8) dqmdt(pcols,pver)
    REAL(r8) dboydt(pcols,pver)
    REAL(r8) thetavp(pcols,pver)
    REAL(r8) thetavm(pcols,pver)

    REAL(r8) dtbdt(pcols),dqbdt(pcols),dtldt(pcols)
    REAL(r8) beta
    REAL(r8) capelmt
    REAL(r8) cp
    REAL(r8) dadt(pcols)
    REAL(r8) debdt
    REAL(r8) dltaa
    REAL(r8) eb
    REAL(r8) grav

    INTEGER i
    INTEGER il1g
    INTEGER il2g
    INTEGER k, kmin, kmax
    INTEGER msg

    REAL(r8) rd
    REAL(r8) rl
    ! change of subcloud layer properties due to convection is
    ! related to cumulus updrafts and downdrafts.
    ! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
    ! to define betau, betad and f(z).
    ! note that this implies all time derivatives are in effect
    ! time derivatives per unit cloud-base mass flux, i.e. they
    ! have units of 1/mb instead of 1/sec.
    !
    DO i = il1g,il2g
       mb(i) = 0.0_r8
       eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
       dtbdt(i) = (1.0_r8/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
            md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
       dqbdt(i) = (1.0_r8/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
            md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
       debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
       dtldt(i) = -2840.0_r8* (3.5_r8/t(i,mx(i))*dtbdt(i)-debdt/eb)/ &
            (3.5_r8*LOG(t(i,mx(i)))-LOG(eb)-4.805_r8)**2
    END DO
    !
    !   dtmdt and dqmdt are cumulus heating and drying.
    !
    DO k = msg + 1,pver
       DO i = il1g,il2g
          dtmdt(i,k) = 0.0_r8
          dqmdt(i,k) = 0.0_r8
       END DO
    END DO
    !
    DO k = msg + 1,pver - 1
       DO i = il1g,il2g
          IF (k == jt(i)) THEN
             dtmdt(i,k) = (1.0_r8/dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- &
                  rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
             dqmdt(i,k) = (1.0_r8/dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                  qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
          END IF
       END DO
    END DO
    !
    beta = 0.0_r8
    DO k = msg + 1,pver - 1
       DO i = il1g,il2g
          IF (k > jt(i) .AND. k < mx(i)) THEN
             dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                  dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
             !          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
             !     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
             !     2                +du(i,k)*(qs(i,k)-q(i,k))
             !     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

             dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                  mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                  (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                  (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                  du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
          END IF
       END DO
    END DO
    !
    DO k = msg + 1,pver
       DO i = il1g,il2g
          IF (k >= lel(i) .AND. k <= lcl(i)) THEN
             thetavp(i,k) = tp(i,k)* (1000.0_r8/p(i,k))** (rd/cp)*(1.0_r8+1.608_r8*qstp(i,k)-q(i,mx(i)))
             thetavm(i,k) = t(i,k)* (1000.0_r8/p(i,k))** (rd/cp)*(1.0_r8+0.608_r8*q(i,k))
             dqsdtp(i,k) = qstp(i,k)* (1.0_r8+qstp(i,k)/eps1)*eps1*rl/(rd*tp(i,k)**2)
             !
             ! dtpdt is the parcel temperature change due to change of
             ! subcloud layer properties during convection.
             !
             dtpdt(i,k) = tp(i,k)/ (1.0_r8+rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                  (dtbdt(i)/t(i,mx(i))+rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/ &
                  tl(i)**2*dtldt(i)))
             !
             ! dboydt is the integrand of cape change.
             !
             dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1.0_r8/(1.0_r8+1.608_r8*qstp(i,k)-q(i,mx(i)))* &
                  (1.608_r8 * dqsdtp(i,k) * dtpdt(i,k) -dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608_r8/ &
                  (1.0_r8+0.608_r8*q(i,k))*dqmdt(i,k)))*grav*thetavp(i,k)/thetavm(i,k)
          END IF
       END DO
    END DO
    !
    DO k = msg + 1,pver
       DO i = il1g,il2g
          IF (k > lcl(i) .AND. k < mx(i)) THEN
             thetavp(i,k) = tp(i,k)* (1000.0_r8/p(i,k))** (rd/cp)*(1.0_r8+0.608_r8*q(i,mx(i)))
             thetavm(i,k) = t(i,k)* (1000.0_r8/p(i,k))** (rd/cp)*(1.0_r8+0.608_r8*q(i,k))
             !
             ! dboydt is the integrand of cape change.
             !
             dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+0.608_r8/ (1.0_r8+0.608_r8*q(i,mx(i)))*dqbdt(i)- &
                  dtmdt(i,k)/t(i,k)-0.608_r8/ (1.0_r8+0.608_r8*q(i,k))*dqmdt(i,k))* &
                  grav*thetavp(i,k)/thetavm(i,k)
          END IF
       END DO
    END DO

    !
    ! buoyant energy change is set to 2/3*excess cape per 3 hours
    !
    dadt(il1g:il2g)  = 0.0_r8
    kmin = MINVAL(lel(il1g:il2g))
    kmax = MAXVAL(mx(il1g:il2g)) - 1
    DO k = kmin, kmax
       DO i = il1g,il2g
          IF ( k >= lel(i) .AND. k <= mx(i) - 1) THEN
             dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
          ENDIF
       END DO
    END DO
    DO i = il1g,il2g
       dltaa = -1.0_r8* (cape(i)-capelmt)
       IF (dadt(i) /= 0.0_r8) mb(i) = MAX(dltaa/tau/dadt(i),0.0_r8)
    END DO
    !
    RETURN
  END SUBROUTINE closure

  SUBROUTINE cldprp(pcols,pver,pverp ,&
       q       ,t       ,u       ,v       ,p       , &
       z       ,s       ,mu      ,eu      ,du      , &
       md      ,ed      ,sd      ,qd      ,mc      , &
       qu      ,su      ,zf      ,qst     ,hmn     , &
       hsat    ,shat    ,ql      , &
       cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
       mx      ,j0      ,jd      ,rl      ,il2g    , &
       rd      ,grav    ,cp      ,msg     , &
       pflx    ,evp     ,cu      ,rprd    ,limcnv  )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
    !             original version cldprop.
    ! 
    ! Author: See above, modified by P. Rasch
    ! This is contributed code not fully standardized by the CCM core group.
    !
    ! this code is very much rougher than virtually anything else in the CCM
    ! there are debug statements left strewn about and code segments disabled
    ! these are to facilitate future development. We expect to release a
    ! cleaner code in a future release
    !
    ! the documentation has been enhanced to the degree that we are able
    !
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                 ! number of vertical levels + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: q(pcols,pver)         ! spec. humidity of env
    REAL(r8), INTENT(in) :: t(pcols,pver)         ! temp of env
    REAL(r8), INTENT(in) :: p(pcols,pver)         ! pressure of env
    REAL(r8), INTENT(in) :: z(pcols,pver)         ! height of env
    REAL(r8), INTENT(in) :: s(pcols,pver)         ! normalized dry static energy of env
    REAL(r8), INTENT(in) :: zf(pcols,pverp)       ! height of interfaces
    REAL(r8), INTENT(in) :: u(pcols,pver)         ! zonal velocity of env
    REAL(r8), INTENT(in) :: v(pcols,pver)         ! merid. velocity of env

    INTEGER, INTENT(in) :: jb(pcols)              ! updraft base level
    INTEGER, INTENT(in) :: lel(pcols)             ! updraft launch level
    INTEGER, INTENT(out) :: jt(pcols)              ! updraft plume top
    INTEGER, INTENT(out) :: jlcl(pcols)            ! updraft lifting cond level
    INTEGER, INTENT(in) :: mx(pcols)              ! updraft base level (same is jb)
    INTEGER, INTENT(out) :: j0(pcols)              ! level where updraft begins detraining
    INTEGER, INTENT(out) :: jd(pcols)              ! level of downdraft
    INTEGER, INTENT(in) :: limcnv                 ! convection limiting level
    INTEGER, INTENT(in) :: il2g                   !CORE GROUP REMOVE
    INTEGER, INTENT(in) :: msg                    ! missing moisture vals (always 0)
    REAL(r8), INTENT(in) :: rl                    ! latent heat of vap
    REAL(r8), INTENT(in) :: shat(pcols,pver)      ! interface values of dry stat energy
    !
    ! output
    !
    REAL(r8), INTENT(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
    REAL(r8), INTENT(out) :: du(pcols,pver)       ! detrainement rate of updraft
    REAL(r8), INTENT(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
    REAL(r8), INTENT(out) :: eu(pcols,pver)       ! entrainment rate of updraft
    REAL(r8), INTENT(out) :: hmn(pcols,pver)      ! moist stat energy of env
    REAL(r8), INTENT(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
    REAL(r8), INTENT(out) :: mc(pcols,pver)       ! net mass flux
    REAL(r8), INTENT(out) :: md(pcols,pver)       ! downdraft mass flux
    REAL(r8), INTENT(out) :: mu(pcols,pver)       ! updraft mass flux
    REAL(r8), INTENT(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
    REAL(r8), INTENT(out) :: qd(pcols,pver)       ! spec humidity of downdraft
    REAL(r8), INTENT(out) :: ql(pcols,pver)       ! liq water of updraft
    REAL(r8), INTENT(out) :: qst(pcols,pver)      ! saturation spec humidity of env.
    REAL(r8), INTENT(out) :: qu(pcols,pver)       ! spec hum of updraft
    REAL(r8), INTENT(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
    REAL(r8), INTENT(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft


    REAL(r8) rd                   ! gas constant for dry air
    REAL(r8) grav                 ! gravity
    REAL(r8) cp                   ! heat capacity of dry air

    !
    ! Local workspace
    !
    REAL(r8) gamma(pcols,pver)
    REAL(r8) dz(pcols,pver)
    REAL(r8) iprm(pcols,pver)
    REAL(r8) hu(pcols,pver)
    REAL(r8) hd(pcols,pver)
    REAL(r8) eps(pcols,pver)
    REAL(r8) f(pcols,pver)
    REAL(r8) k1(pcols,pver)
    REAL(r8) i2(pcols,pver)
    REAL(r8) ihat(pcols,pver)
    REAL(r8) i3(pcols,pver)
    REAL(r8) idag(pcols,pver)
    REAL(r8) i4(pcols,pver)
    REAL(r8) qsthat(pcols,pver)
    REAL(r8) hsthat(pcols,pver)
    REAL(r8) gamhat(pcols,pver)
    REAL(r8) cu(pcols,pver)
    REAL(r8) evp(pcols,pver)
    REAL(r8) cmeg(pcols,pver)
    REAL(r8) qds(pcols,pver)
    REAL(r8) hmin(pcols)
    REAL(r8) expdif(pcols)
    REAL(r8) expnum(pcols)
    REAL(r8) ftemp(pcols)
    REAL(r8) eps0(pcols)
    REAL(r8) rmue(pcols)
    REAL(r8) zuef(pcols)
    REAL(r8) zdef(pcols)
    REAL(r8) epsm(pcols)
    REAL(r8) ratmjb(pcols)
    REAL(r8) est(pcols)
    REAL(r8) totpcp(pcols)
    REAL(r8) totevp(pcols)
    REAL(r8) alfa(pcols)
    REAL(r8) ql1
    REAL(r8) tu
    REAL(r8) estu
    REAL(r8) qstu

    REAL(r8) small
    REAL(r8) mdt

    INTEGER khighest
    INTEGER klowest
    INTEGER kount
    INTEGER i,k

    LOGICAL doit(pcols)
    LOGICAL done(pcols)
    !
    !------------------------------------------------------------------------------
    !
    DO i = 1,il2g
       ftemp(i) = 0.0_r8
       expnum(i) = 0.0_r8
       expdif(i) = 0.0_r8
    END DO
    !
    !jr Change from msg+1 to 1 to prevent blowup
    !
    DO k = 1,pver
       DO i = 1,il2g
          dz(i,k) = zf(i,k) - zf(i,k+1)
       END DO
    END DO

    !
    ! initialize many output and work variables to zero
    !
    pflx(:il2g,1) = 0

    DO k = 1,pver
       DO i = 1,il2g
          k1(i,k) = 0.0_r8
          i2(i,k) = 0.0_r8
          i3(i,k) = 0.0_r8
          i4(i,k) = 0.0_r8
          mu(i,k) = 0.0_r8
          f(i,k) = 0.0_r8
          eps(i,k) = 0.0_r8
          eu(i,k) = 0.0_r8
          du(i,k) = 0.0_r8
          ql(i,k) = 0.0_r8
          cu(i,k) = 0.0_r8
          evp(i,k) = 0.0_r8
          cmeg(i,k) = 0.0_r8
          qds(i,k) = q(i,k)
          md(i,k) = 0.0_r8
          ed(i,k) = 0.0_r8
          sd(i,k) = s(i,k)
          qd(i,k) = q(i,k)
          mc(i,k) = 0.0_r8
          qu(i,k) = q(i,k)
          su(i,k) = s(i,k)
          !        est(i)=exp(a-b/t(i,k))
          est(i) = c1*EXP((c2* (t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
          !++bee
          IF ( p(i,k)-est(i) > 0.0_r8 ) THEN
             qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
          ELSE
             qst(i,k) = 1.0_r8
          END IF
          !--bee
          gamma(i,k) = qst(i,k)*(1.0_r8 + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
          hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
          hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
          hu(i,k) = hmn(i,k)
          hd(i,k) = hmn(i,k)
          rprd(i,k) = 0.0_r8
       END DO
    END DO
    !
    !jr Set to zero things which make this routine blow up
    !
    DO k=1,msg
       DO i=1,il2g
          rprd(i,k) = 0.0_r8
       END DO
    END DO
    !
    ! interpolate the layer values of qst, hsat and gamma to
    ! layer interfaces
    !
    DO i = 1,il2g
       hsthat(i,msg+1) = hsat(i,msg+1)
       qsthat(i,msg+1) = qst(i,msg+1)
       gamhat(i,msg+1) = gamma(i,msg+1)
       totpcp(i) = 0.0_r8
       totevp(i) = 0.0_r8
    END DO
    DO k = msg + 2,pver
       DO i = 1,il2g
          IF (ABS(qst(i,k-1)-qst(i,k)) > 1.0E-6_r8) THEN
             qsthat(i,k) = LOG(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
          ELSE
             qsthat(i,k) = qst(i,k)
          END IF
          hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
          IF (ABS(gamma(i,k-1)-gamma(i,k)) > 1.0E-6_r8) THEN
             gamhat(i,k) = LOG(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                  (gamma(i,k-1)-gamma(i,k))
          ELSE
             gamhat(i,k) = gamma(i,k)
          END IF
       END DO
    END DO
    !
    ! initialize cloud top to highest plume top.
    !jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
    !
   jt(:) = pver
    DO i = 1,il2g
       jt(i) = MAX(lel(i),limcnv+1)
       jt(i) = MIN(jt(i),pver)
       jd(i) = pver
       jlcl(i) = lel(i)
       hmin(i) = 1.0E6_r8
    END DO
    !
    ! find the level of minimum hsat, where detrainment starts
    !
    DO k = msg + 1,pver
       DO i = 1,il2g
          IF (hsat(i,k) <= hmin(i) .AND. k >= jt(i) .AND. k <= jb(i)) THEN
             hmin(i) = hsat(i,k)
             j0(i) = k
          END IF
       END DO
    END DO
    DO i = 1,il2g
       j0(i) = MIN(j0(i),jb(i)-2)
       j0(i) = MAX(j0(i),jt(i)+2)
       !
       ! Fix from Guang Zhang to address out of bounds array reference
       !
       j0(i) = MIN(j0(i),pver)
    END DO
    !
    ! Initialize certain arrays inside cloud
    !
    DO k = msg + 1,pver
       DO i = 1,il2g
          IF (k >= jt(i) .AND. k <= jb(i)) THEN
             hu(i,k) = hmn(i,mx(i)) + cp*0.5_r8
             su(i,k) = s(i,mx(i)) + 0.5_r8
          END IF
       END DO
    END DO
    !
    ! *********************************************************
    ! compute taylor series for approximate eps(z) below
    ! *********************************************************
    !
    DO k = pver - 1,msg + 1,-1
       DO i = 1,il2g
          IF (k < jb(i) .AND. k >= jt(i)) THEN
             k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
             ihat(i,k) = 0.5_r8* (k1(i,k+1)+k1(i,k))
             i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
             idag(i,k) = 0.5_r8* (i2(i,k+1)+i2(i,k))
             i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
             iprm(i,k) = 0.5_r8* (i3(i,k+1)+i3(i,k))
             i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
          END IF
       END DO
    END DO
    !
    ! re-initialize hmin array for ensuing calculation.
    !
    DO i = 1,il2g
       hmin(i) = 1.0E6_r8
    END DO
    DO k = msg + 1,pver
       DO i = 1,il2g
          IF (k >= j0(i) .AND. k <= jb(i) .AND. hmn(i,k) <= hmin(i)) THEN
             hmin(i) = hmn(i,k)
             expdif(i) = hmn(i,mx(i)) - hmin(i)
          END IF
       END DO
    END DO
    !
    ! *********************************************************
    ! compute approximate eps(z) using above taylor series
    ! *********************************************************
    !
    DO k = msg + 2,pver
       DO i = 1,il2g
          expnum(i) = 0.0_r8
          ftemp(i) = 0.0_r8
          IF (k < jt(i) .OR. k >= jb(i)) THEN
             k1(i,k) = 0.0_r8
             expnum(i) = 0.0_r8
          ELSE
             expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) + &
                  hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
          END IF
          IF ((expdif(i) > 100.0_r8 .AND. expnum(i) > 0.0_r8) .AND. &
               k1(i,k) > expnum(i)*dz(i,k)) THEN
             ftemp(i) = expnum(i)/k1(i,k)
             f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
                  (2.0_r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
                  ftemp(i)**3 + (-5.0_r8*k1(i,k)*i2(i,k)*i3(i,k)+ &
                  5.0_r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
                  k1(i,k)**3*ftemp(i)**4
             f(i,k) = MAX(f(i,k),0.0_r8)
             f(i,k) = MIN(f(i,k),0.0002_r8)
          END IF
       END DO
    END DO
    DO i = 1,il2g
       IF (j0(i) < jb(i)) THEN
          IF (f(i,j0(i)) < 1.0E-6_r8 .AND. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
       END IF
    END DO
    DO k = msg + 2,pver
       DO i = 1,il2g
          IF (k >= jt(i) .AND. k <= j0(i)) THEN
             f(i,k) = MAX(f(i,k),f(i,k-1))
          END IF
       END DO
    END DO
    DO i = 1,il2g
       eps0(i) = f(i,j0(i))
       eps(i,jb(i)) = eps0(i)
    END DO
    !
    ! This is set to match the Rasch and Kristjansson paper
    !
    DO k = pver,msg + 1,-1
       DO i = 1,il2g
          IF (k >= j0(i) .AND. k <= jb(i)) THEN
             eps(i,k) = f(i,j0(i))
          END IF
       END DO
    END DO
    DO k = pver,msg + 1,-1
       DO i = 1,il2g
          IF (k < j0(i) .AND. k >= jt(i)) eps(i,k) = f(i,k)
       END DO
    END DO
    !
    ! specify the updraft mass flux mu, entrainment eu, detrainment du
    ! and moist static energy hu.
    ! here and below mu, eu,du, md and ed are all normalized by mb
    !
    DO i = 1,il2g
       IF (eps0(i) > 0.0_r8) THEN
          mu(i,jb(i)) = 1.0_r8
          eu(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
       END IF
    END DO
    DO k = pver,msg + 1,-1
       DO i = 1,il2g
          IF (eps0(i) > 0.0_r8 .AND. (k >= jt(i) .AND. k < jb(i))) THEN
             zuef(i) = zf(i,k) - zf(i,jb(i))
             rmue(i) = (1.0_r8/eps0(i))* (EXP(eps(i,k+1)*zuef(i))-1.0_r8)/zuef(i)
             mu(i,k) = (1.0_r8/eps0(i))* (EXP(eps(i,k  )*zuef(i))-1.0_r8)/zuef(i)
             eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
             du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
          END IF
       END DO
    END DO
    !
    khighest = pverp
    klowest = 1
    DO i=1,il2g
       khighest = MIN(khighest,lel(i))
       klowest = MAX(klowest,jb(i))
    END DO
    DO k = klowest-1,khighest,-1
       !cdir$ ivdep
       DO i = 1,il2g
          IF (k <= jb(i)-1 .AND. k >= lel(i) .AND. eps0(i) > 0.0_r8) THEN
             IF (mu(i,k) < 0.01_r8) THEN
                hu(i,k) = hu(i,jb(i))
                mu(i,k) = 0.0_r8
                eu(i,k) = 0.0_r8
                du(i,k) = mu(i,k+1)/dz(i,k)
             ELSE
                hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                     dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)- du(i,k)*hsat(i,k))
             END IF
          END IF
       END DO
    END DO
    !
    ! reset cloud top index beginning from two layers above the
    ! cloud base (i.e. if cloud is only one layer thick, top is not reset
    !
    DO i=1,il2g
       doit(i) = .TRUE.
    END DO
    DO k=klowest-2,khighest-1,-1
       DO i=1,il2g
          IF (doit(i) .AND. k <= jb(i)-2 .AND. k >= lel(i)-1) THEN
             IF (hu(i,k) <= hsthat(i,k) .AND. hu(i,k+1) > hsthat(i,k+1) &
                  .AND. mu(i,k) >= 0.02_r8) THEN
                IF (hu(i,k)-hsthat(i,k) < -2000.0_r8) THEN
                   jt(i) = k + 1
                   doit(i) = .FALSE.
                ELSE
                   jt(i) = k
                   IF (eps0(i) <= 0.0_r8) doit(i) = .FALSE.
                END IF
             ELSE IF (hu(i,k) > hu(i,jb(i)) .OR. mu(i,k) < 0.01_r8) THEN
                jt(i) = k + 1
                doit(i) = .FALSE.
             END IF
          END IF
       END DO
    END DO
    DO k = pver,msg + 1,-1
       !cdir$ ivdep
       DO i = 1,il2g
          IF (k >= lel(i) .AND. k <= jt(i) .AND. eps0(i) > 0.0_r8) THEN
             mu(i,k) = 0.0_r8
             eu(i,k) = 0.0_r8
             du(i,k) = 0.0_r8
             hu(i,k) = hu(i,jb(i))
          END IF
          IF (k == jt(i) .AND. eps0(i) > 0.0_r8) THEN
             du(i,k) = mu(i,k+1)/dz(i,k)
             eu(i,k) = 0.0_r8
             mu(i,k) = 0.0_r8
          END IF
       END DO
    END DO
    !
    ! specify downdraft properties (no downdrafts if jd.ge.jb).
    ! scale down downward mass flux profile so that net flux
    ! (up-down) at cloud base in not negative.
    !
    DO i = 1,il2g
       !
       ! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
       !
       alfa(i) = 0.1_r8
       jt(i) = MIN(jt(i),jb(i)-1)
       jd(i) = MAX(j0(i),jt(i)+1)
       jd(i) = MIN(jd(i),jb(i))
       hd(i,jd(i)) = hmn(i,jd(i)-1)
       IF (jd(i) < jb(i) .AND. eps0(i) > 0.0_r8) THEN
          epsm(i) = eps0(i)
          md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
       END IF
    END DO
    DO k = msg + 1,pver
       DO i = 1,il2g
          IF ((k > jd(i) .AND. k <= jb(i)) .AND. eps0(i) > 0.0_r8) THEN
             zdef(i) = zf(i,jd(i)) - zf(i,k)
             md(i,k) = -alfa(i)/ (2.0_r8*eps0(i))*(EXP(2.0_r8*epsm(i)*zdef(i))-1.0_r8)/zdef(i)
          END IF
       END DO
    END DO
    DO k = msg + 1,pver
       DO i = 1,il2g
          IF ((k >= jt(i) .AND. k <= jb(i)) .AND. eps0(i) > 0.0_r8 .AND. jd(i) < jb(i)) THEN
             ratmjb(i) = MIN(ABS(mu(i,jb(i))/md(i,jb(i))),1.0_r8)
             md(i,k) = md(i,k)*ratmjb(i)
          END IF
       END DO
    END DO

    small = 1.0e-20_r8
    DO k = msg + 1,pver
       DO i = 1,il2g
          IF ((k >= jt(i) .AND. k <= pver) .AND. eps0(i) > 0.0_r8) THEN
             ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
             mdt = MIN(md(i,k),-small)
             hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
          END IF
       END DO
    END DO
    !
    ! calculate updraft and downdraft properties.
    !
    DO k = msg + 2,pver
       DO i = 1,il2g
          IF ((k >= jd(i) .AND. k <= jb(i)) .AND. eps0(i) > 0. .AND. jd(i) < jb(i)) THEN
             qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/ &
                  (rl*(1.0_r8 + gamhat(i,k)))
          END IF
       END DO
    END DO
    !
    DO i = 1,il2g
       done(i) = .FALSE.
    END DO
    kount = 0
    DO k = pver,msg + 2,-1
       DO i = 1,il2g
          IF (( .NOT. done(i) .AND. k > jt(i) .AND. k < jb(i)) .AND. eps0(i) > 0.0_r8) THEN
             su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                  dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
             qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)- &
                  du(i,k)*qst(i,k))
             tu = su(i,k) - grav/cp*zf(i,k)
             estu = c1*EXP((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
             qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2.0_r8-estu)
             IF (qu(i,k) >= qstu) THEN
                jlcl(i) = k
                kount = kount + 1
                done(i) = .TRUE.
             END IF
          END IF
       END DO
       IF (kount >= il2g) GOTO 690
    END DO
690 CONTINUE
    DO k = msg + 2,pver
       DO i = 1,il2g
          IF (k == jb(i) .AND. eps0(i) > 0.0_r8) THEN
             qu(i,k) = q(i,mx(i))
             su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
          END IF
          IF ((k > jt(i) .AND. k <= jlcl(i)) .AND. eps0(i) > 0.0_r8) THEN
             su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/(cp* (1.0_r8+gamhat(i,k)))
             qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k))/ &
                  (rl* (1.0_r8+gamhat(i,k)))
          END IF
       END DO
    END DO

    ! compute condensation in updraft
    DO k = pver,msg + 2,-1
       DO i = 1,il2g
          IF (k >= jt(i) .AND. k < jb(i) .AND. eps0(i) > 0.0_r8) THEN
             cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                  dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/(rl/cp)
             IF (k == jt(i)) cu(i,k) = 0.0_r8
             cu(i,k) = MAX(0.0_r8,cu(i,k))
          END IF
       END DO
    END DO

    ! compute condensed liquid, rain production rate
    ! accumulate total precipitation (condensation - detrainment of liquid)
    ! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
    ! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
    ! consistently applied.
    !    mu, ql are interface quantities
    !    cu, du, eu, rprd are midpoint quantites
    DO k = pver,msg + 2,-1
       DO i = 1,il2g
          rprd(i,k) = 0.0_r8
          IF (k >= jt(i) .AND. k < jb(i) .AND. eps0(i) > 0.0_r8 .AND. mu(i,k) >= 0.0_r8) THEN
             IF (mu(i,k) > 0.0_r8) THEN
                ql1 = 1.0_r8/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
                     dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
                ql(i,k) = ql1/ (1.0_r8+dz(i,k)*c0)
             ELSE
                ql(i,k) = 0.0_r8
             END IF
             totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))
             rprd(i,k) = c0*mu(i,k)*ql(i,k)
          END IF
       END DO
    END DO
    !
    DO i = 1,il2g
       qd(i,jd(i)) = qds(i,jd(i))
       sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
    END DO
    !
    DO k = msg + 2,pver
       DO i = 1,il2g
          IF (k >= jd(i) .AND. k < jb(i) .AND. eps0(i) > 0.0_r8) THEN
             qd(i,k+1) = qds(i,k+1)
             evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
             evp(i,k) = MAX(evp(i,k),0.0_r8)
             mdt = MIN(md(i,k+1),-small)
             sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
             totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          END IF
       END DO
    END DO
    DO i = 1,il2g
       !*guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
       totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
    END DO
!!$   if (.true.) then
    IF (.FALSE.) THEN
       DO i = 1,il2g
          k = jb(i)
          IF (eps0(i) > 0.0_r8) THEN
             evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
             evp(i,k) = MAX(evp(i,k),0.0_r8)
             totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          END IF
       END DO
    ENDIF

    DO i = 1,il2g
       totpcp(i) = MAX(totpcp(i),0.0_r8)
       totevp(i) = MAX(totevp(i),0.0_r8)
    END DO
    !
    DO k = msg + 2,pver
       DO i = 1,il2g
          IF (totevp(i) > 0.0_r8 .AND. totpcp(i) > 0.0_r8) THEN
             md(i,k)  = md (i,k)*MIN(1.0_r8, totpcp(i)/(totevp(i)+totpcp(i)))
             ed(i,k)  = ed (i,k)*MIN(1.0_r8, totpcp(i)/(totevp(i)+totpcp(i)))
             evp(i,k) = evp(i,k)*MIN(1.0_r8, totpcp(i)/(totevp(i)+totpcp(i)))
          ELSE
             md(i,k) = 0.0_r8
             ed(i,k) = 0.0_r8
             evp(i,k) = 0.0_r8
          END IF
          ! cmeg is the cloud water condensed - rain water evaporated
          ! rprd is the cloud water converted to rain - (rain evaporated)
          cmeg(i,k) = cu(i,k) - evp(i,k)
          rprd(i,k) = rprd(i,k)-evp(i,k)
       END DO
    END DO

    ! compute the net precipitation flux across interfaces
    pflx(:il2g,1) = 0.0_r8
    DO k = 2,pverp
       DO i = 1,il2g
          pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
       END DO
    END DO
    !
    DO k = msg + 1,pver
       DO i = 1,il2g
          mc(i,k) = mu(i,k) + md(i,k)
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE cldprp

  !=========================================================================================
subroutine buoyan_dilute(lchnk   ,ncol    ,pcols,pver, &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  rd      ,grav    ,cp      ,msg     , &
                  tpert   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
! 
   implicit none




   integer, intent(in) :: lchnk                 
   integer, intent(in) :: ncol   ! number of atmospheric columns               
   integer, intent(in) :: pcols
   integer, intent(in) :: pver
   real(r8), intent(in) :: q(pcols,pver)        ! spec. humidity
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: p(pcols,pver)        ! pressure
   real(r8), intent(in) :: z(pcols,pver)        ! height
   real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(pcols)          ! index of pbl depth
   real(r8), intent(in) :: tpert(pcols)   ! perturbation temperature by pbl processes      
!
! output arguments
!
   real(r8), intent(out) :: tp(pcols,pver)    ! parcel temperature   
   real(r8), intent(out) :: qstp(pcols,pver)   ! saturation mixing ratio of parcel (only above lcl, just q below).  
   real(r8), intent(out) :: tl(pcols)          ! parcel temperature at lcl  
   real(r8), intent(out) :: cape(pcols)         ! convective aval. pot. energy. 
   integer lcl(pcols)        
   integer lel(pcols)        
   integer lon(pcols)! level of onset of deep convection        
   integer mx(pcols) ! level of max moist static energy        

!
!--------------------------Local Variables------------------------------
!


   real(r8) capeten(pcols,5)    ! provisional value of cape  
   real(r8) tv(pcols,pver)       
   real(r8) tpv(pcols,pver)      
   real(r8) buoy(pcols,pver)

   real(r8) a1(pcols)
   real(r8) a2(pcols)
   real(r8) estp(pcols)
   real(r8) pl(pcols)
   real(r8) plexp(pcols)
   real(r8) hmax(pcols)
   real(r8) hmn(pcols)
   real(r8) y(pcols)

   logical plge600(pcols)
   integer knt(pcols)
   integer lelten(pcols,5)

   real(r8) cp
   real(r8) e
   real(r8) grav

   integer i
   integer k
   integer msg
   integer n

   real(r8) rd
   real(r8) rl

   real(r8) rhd



   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = pver
         capeten(i,n) = 0._r8
      end do
   end do

   do i = 1,ncol
      lon(i) = pver
      knt(i) = 0
      lel(i) = pver
      mx(i) = lon(i)
      cape(i) = 0._r8
      hmax(i) = 0._r8
   end do

   tp(:ncol,:) = t(:ncol,:)
   qstp(:ncol,:) = q(:ncol,:)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.

   tv(:ncol,:) = t(:ncol,:) *(1._r8+1.608_r8*q(:ncol,:))/ (1._r8+q(:ncol,:))
   tpv(:ncol,:) = tv(:ncol,:)
   buoy(:ncol,:) = 0._r8


!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!

IF(PERGRO)THEN
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
!
! Reset max moist static energy level when relative difference exceeds 1.e-4
!
         rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.e-4_r8) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
ELSE
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
END IF

! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.
   do i = 1,ncol ! Initialise LCL variables.
      lcl(i) = mx(i)
      tl(i) = t(i,mx(i))
      pl(i) = p(i,mx(i))
   end do
!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!

   call parcel_dilute(lchnk, ncol,pcols,pver, msg, mx, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)

! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._r8 ! Just change to always allow buoy calculation.
   end do
!
! Main buoyancy calculation.
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
            tv(i,k) = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add  ! +0.5K or not?
         else
            qstp(i,k) = q(i,k)
            tp(i,k)   = t(i,k)            
            tpv(i,k)  = tv(i,k)
         endif
      end do
   end do



!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
   do k = msg + 2,pver
      do i = 1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1) > 0. .and. buoy(i,k) <= 0._r8) then
               knt(i) = min(5,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do k = msg + 1,pver
         do i = 1,ncol
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do

   return
end subroutine buoyan_dilute

  !=========================================================================================

  SUBROUTINE buoyan( pcols,pver,&
       q       ,t       ,p       ,z       ,pf      , &
       tp      ,qstp    ,tl      ,rl      ,cape    , &
       pblt    ,lcl     ,lel     ,lon     ,mx      , &
       rd      ,grav    ,cp      ,msg     , &
       tpert   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author:
    ! This is contributed code not fully standardized by the CCM core group.
    ! The documentation has been enhanced to the degree that we are able.
    ! Reviewed:          P. Rasch, April 1996
    ! 
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !
    ! input arguments
    !
    INTEGER, INTENT(in) :: pcols                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: q(pcols,pver)        ! spec. humidity
    REAL(r8), INTENT(in) :: t(pcols,pver)        ! temperature
    REAL(r8), INTENT(in) :: p(pcols,pver)        ! pressure
    REAL(r8), INTENT(in) :: z(pcols,pver)        ! height
    REAL(r8), INTENT(in) :: pf(pcols,pver+1)     ! pressure at interfaces
    REAL(r8), INTENT(in) :: pblt(pcols)          ! index of pbl depth
    REAL(r8), INTENT(in) :: tpert(pcols)         ! perturbation temperature by pbl processes

    !
    ! output arguments
    !
    REAL(r8), INTENT(out) :: tp(pcols,pver)       ! parcel temperature
    REAL(r8), INTENT(out) :: qstp(pcols,pver)     ! saturation mixing ratio of parcel
    REAL(r8), INTENT(out) :: tl(pcols)            ! parcel temperature at lcl
    REAL(r8), INTENT(out) :: cape(pcols)          ! convective aval. pot. energy.
    INTEGER lcl(pcols)        !
    INTEGER lel(pcols)        !
    INTEGER lon(pcols)        ! level of onset of deep convection
    INTEGER mx(pcols)         ! level of max moist static energy
    !
    !--------------------------Local Variables------------------------------
    !
    REAL(r8) capeten(pcols,5)     ! provisional value of cape
    REAL(r8) tv(pcols,pver)       !
    REAL(r8) tpv(pcols,pver)      !
    REAL(r8) buoy(pcols,pver)

    REAL(r8) a1(pcols)
    REAL(r8) a2(pcols)
    REAL(r8) estp(pcols)
    REAL(r8) pl(pcols)
    REAL(r8) plexp(pcols)
    REAL(r8) hmax(pcols)
    REAL(r8) hmn(pcols)
    REAL(r8) y(pcols)

    LOGICAL plge600(pcols)
    INTEGER knt(pcols)
    INTEGER lelten(pcols,5)

    REAL(r8) cp
    REAL(r8) e
    REAL(r8) grav

    INTEGER i
    INTEGER k
    INTEGER msg
    INTEGER n

    REAL(r8) rd
    REAL(r8) rl
    REAL(r8) rhd
    !
    !-----------------------------------------------------------------------
    !
    DO n = 1,5
       DO i = 1,pcols
          lelten(i,n) = pver
          capeten(i,n) = 0.0_r8
       END DO
    END DO
    !
    DO i = 1,pcols
       lon(i) = pver
       knt(i) = 0
       lel(i) = pver
       mx(i) = lon(i)
       cape(i) = 0.0_r8
       hmax(i) = 0.0_r8
    END DO

    tp(:pcols,:) = t(:pcols,:)
    qstp(:pcols,:) = q(:pcols,:)
    !
    ! set "launching" level(mx) to be at maximum moist static energy.
    ! search for this level stops at planetary boundary layer top.
    !
    !#ifdef PERGRO
    IF(PERGRO)THEN
       DO k = pver,msg + 1,-1
          DO i = 1,pcols
             hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
             !
             ! Reset max moist static energy level when relative difference exceeds 1.e-4
             !
             rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
             IF (k >= NINT(pblt(i)) .AND. k <= lon(i) .AND. rhd > -1.0e-4_r8) THEN
                hmax(i) = hmn(i)
                mx(i) = k
             END IF
          END DO
       END DO
    ELSE
       !#else
       DO k = pver,msg + 1,-1
          DO i = 1,pcols
             hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
             IF (k >= NINT(pblt(i)) .AND. k <= lon(i) .AND. hmn(i) > hmax(i)) THEN
                hmax(i) = hmn(i)
                mx(i) = k
             END IF
          END DO
       END DO
       !#endif
    ENDIF
    !
    DO i = 1,pcols
       lcl(i) = mx(i)
       e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
       tl(i) = 2840.0_r8/ (3.5_r8*LOG(t(i,mx(i)))-LOG(e)-4.805_r8) + 55.0_r8
       IF (tl(i) < t(i,mx(i))) THEN
          plexp(i) = (1.0_r8/ (0.2854_r8 * (1.0_r8-0.28_r8*q(i,mx(i)))))
          pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)
       ELSE
          tl(i) = t(i,mx(i))
          pl(i) = p(i,mx(i))
       END IF
    END DO
    !
    ! calculate lifting condensation level (lcl).
    !
    DO k = pver,msg + 2,-1
       DO i = 1,pcols
          IF (k <= mx(i) .AND. (p(i,k) > pl(i) .AND. p(i,k-1) <= pl(i))) THEN
             lcl(i) = k - 1
          END IF
       END DO
    END DO
    !
    ! if lcl is above the nominal level of non-divergence (600 mbs),
    ! no deep convection is permitted (ensuing calculations
    ! skipped and cape retains initialized value of zero).
    !
    DO i = 1,pcols
       plge600(i) = pl(i).GE.600.0_r8
    END DO
    !
    ! initialize parcel properties in sub-cloud layer below lcl.
    !
    DO k = pver,msg + 1,-1
       DO i=1,pcols
          IF (k > lcl(i) .AND. k <= mx(i) .AND. plge600(i)) THEN
             tv(i,k) = t(i,k)* (1.0_r8+1.608_r8*q(i,k))/ (1.0_r8+q(i,k))
             qstp(i,k) = q(i,mx(i))
             tp(i,k) = t(i,mx(i))* (p(i,k)/p(i,mx(i)))**(0.2854_r8* (1.0_r8-0.28_r8*q(i,mx(i))))
             !
             ! buoyancy is increased by 0.5_r8 k as in tiedtke
             !
             !-jjh          tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))    / (1.+q(i,mx(i)))
             !-kubota tpv(i,k) = (tp(i,k)+tpert(i))*(1.0_r8+1.608_r8*q(i,mx(i)))/ (1.0_r8+q(i,mx(i)))
             tpv(i,k) = (tp(i,k))*(1.0_r8+1.608_r8*q(i,mx(i)))/ (1.0_r8+q(i,mx(i)))

             buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5_r8
          END IF
       END DO
    END DO
    !
    ! define parcel properties at lcl (i.e. level immediately above pl).
    !
    DO k = pver,msg + 1,-1
       DO i=1,pcols
          IF (k == lcl(i) .AND. plge600(i)) THEN
             tv(i,k) = t(i,k)* (1.0_r8+1.608_r8*q(i,k))/ (1.0_r8+q(i,k))
             qstp(i,k) = q(i,mx(i))
             tp(i,k) = tl(i)* (p(i,k)/pl(i))**(0.2854_r8* (1.0_r8-0.28_r8*qstp(i,k)))
             !              estp(i)  =exp(a-b/tp(i,k))
             ! use of different formulas for est has about 1 g/kg difference
             ! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
             ! above giving larger qs.
             !
             estp(i) = c1*EXP((c2* (tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))

             qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
             a1(i) = cp / rl + qstp(i,k) * (1.0_r8+ qstp(i,k) / eps1) * rl * eps1 / &
                  (rd * tp(i,k) ** 2)
             a2(i) = 0.5_r8* (qstp(i,k)* (1.0_r8+2.0_r8/eps1*qstp(i,k))* &
                  (1.0_r8+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                  (rd**2*tp(i,k)**4)-qstp(i,k)* &
                  (1.0_r8+qstp(i,k)/eps1)*2.0_r8*eps1*rl/ &
                  (rd*tp(i,k)**3))
             a1(i) = 1.0_r8/a1(i)
             a2(i) = -a2(i)*a1(i)**3
             y(i) = q(i,mx(i)) - qstp(i,k)
             tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
             !          estp(i)  =exp(a-b/tp(i,k))
             estp(i) = c1*EXP((c2* (tp(i,k)-tfreez))/ ((tp(i,k)-tfreez)+c3))

             qstp(i,k) = eps1*estp(i) / (p(i,k)-estp(i))
             !
             ! buoyancy is increased by 0.5_r8 k in cape calculation.
             ! dec. 9, 1994
             !-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))
             !
             !kubota tpv(i,k) = (tp(i,k)+tpert(i))* (1.0_r8+1.608_r8*qstp(i,k)) / (1.0_r8+q(i,mx(i)))
             tpv(i,k) = (tp(i,k))* (1.0_r8+1.608_r8*qstp(i,k)) / (1.0_r8+q(i,mx(i)))
             buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5_r8
          END IF
       END DO
    END DO
    !
    ! main buoyancy calculation.
    !
    DO k = pver - 1,msg + 1,-1
       DO i=1,pcols
          IF (k < lcl(i) .AND. plge600(i)) THEN
             tv(i,k) = t(i,k)* (1.0_r8+1.608_r8*q(i,k))/ (1.0_r8+q(i,k))
             qstp(i,k) = qstp(i,k+1)
             tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**(0.2854_r8* (1.0_r8-0.28_r8*qstp(i,k)))
             !          estp(i) = exp(a-b/tp(i,k))
             estp(i) = c1*EXP((c2* (tp(i,k)-tfreez))/((tp(i,k)-tfreez)+c3))

             qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
             a1(i) = cp/rl + qstp(i,k)* (1.0_r8+qstp(i,k)/eps1)*rl*eps1/ (rd*tp(i,k)**2)
             a2(i) = 0.5_r8* (qstp(i,k)* (1.0_r8+2.0_r8/eps1*qstp(i,k))* &
                  (1.0_r8+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                  (rd**2*tp(i,k)**4)-qstp(i,k)* &
                  (1.0_r8+qstp(i,k)/eps1)*2.0_r8*eps1*rl/ &
                  (rd*tp(i,k)**3))
             a1(i) = 1.0_r8/a1(i)
             a2(i) = -a2(i)*a1(i)**3
             y(i) = qstp(i,k+1) - qstp(i,k)
             tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
             !          estp(i)  =exp(a-b/tp(i,k))
             estp(i) = c1*EXP((c2* (tp(i,k)-tfreez))/ ((tp(i,k)-tfreez)+c3))

             qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
             !-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
             !jt            (1.+q(i,mx(i)))
             !kubota tpv(i,k) = (tp(i,k)+tpert(i))* (1.0_r8+1.608_r8*qstp(i,k))/(1.0_r8+q(i,mx(i)))
             tpv(i,k) = (tp(i,k))* (1.0_r8+1.608_r8*qstp(i,k))/(1.0_r8+q(i,mx(i)))
             buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5_r8
          END IF
       END DO
    END DO
    !
    DO k = msg + 2,pver
       DO i = 1,pcols
          IF (k < lcl(i) .AND. plge600(i)) THEN
             IF (buoy(i,k+1) > 0.0_r8 .AND. buoy(i,k) <= 0.0_r8) THEN
                knt(i) = MIN(5,knt(i) + 1)
                lelten(i,knt(i)) = k
             END IF
          END IF
       END DO
    END DO
    !
    ! calculate convective available potential energy (cape).
    !
    DO n = 1,5
       DO k = msg + 1,pver
          DO i = 1,pcols
             IF (plge600(i) .AND. k <= mx(i) .AND. k > lelten(i,n)) THEN
                capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*LOG(pf(i,k+1)/pf(i,k))
             END IF
          END DO
       END DO
    END DO
    !
    ! find maximum cape from all possible tentative capes from
    ! one sounding,
    ! and use it as the final cape, april 26, 1995
    !
    DO n = 1,5
       DO i = 1,pcols
          IF (capeten(i,n) > cape(i)) THEN
             cape(i) = capeten(i,n)
             lel(i) = lelten(i,n)
          END IF
       END DO
    END DO
    !
    ! put lower bound on cape for diagnostic purposes.
    !
    DO i = 1,pcols
       cape(i) = MAX(cape(i), 0.0_r8)
    END DO
    !
    RETURN
  END SUBROUTINE buoyan

  SUBROUTINE zm_convi(limcnv_in,plat,dt)

    INTEGER      , INTENT(in   ) :: limcnv_in     ! top interface level limit for convection
    INTEGER      , INTENT(IN   ) :: plat          ! number of latitudes
    REAL(KIND=r8), INTENT(IN   ) :: dt
    !
    ! Initialization of ZM constants
    !
    limcnv = limcnv_in
    tfreez = tmelt
    eps1   = epsilo
    rl     = latvap
    cpres  = cpair
    rgrav = 1.0_r8/gravit
    rgas = rair
    grav = gravit
    cp = cpres
    !
    ! tau=4800. were used in canadian climate center. however, in echam3 t42, 
    ! convection is too weak, thus adjusted to 2400.
    !
    IF ( dycore_is ('LR') ) THEN
       IF ( get_resolution(plat) == '1x1.25' ) THEN
          IF(masterproc) WRITE(6,*) 'ZM: Found LR dycore at 1x1.25'
          tau = 3600.0_r8
          c0 = 3.5E-3_r8
          ke = 1.0E-6_r8
       ELSEIF ( get_resolution(plat) == '4x5' ) THEN
          IF(masterproc) WRITE(6,*) 'ZM: Found LR dycore at 4x5'
          tau = 3600.0_r8
          c0 = 3.5E-3_r8
          ke = 1.0E-6_r8
       ELSE
          IF(masterproc) WRITE(6,*) 'ZM: Found LR dycore at default resolution'
          tau = 3600.0_r8
          c0 = 3.5E-3_r8
          ke = 1.0E-6_r8
       ENDIF
    ELSE
       IF(get_resolution(plat) == 'T85')THEN
          tau = 3600.0_r8
          c0 = 4.E-3_r8
          ke = 1.0E-6_r8
          IF(masterproc) WRITE(6,*) 'ZM: Found spectral dycore at T85 resolution'
       ELSEIF(get_resolution(plat) == 'T31')THEN
          tau  = 3600.0_r8
          c0 = 2.E-3_r8
          ke = 3.0E-6_r8
          IF(masterproc) WRITE(6,*) 'ZM: Found spectral dycore at non-T85 resolution'
       ELSE
          tau  = 3600.0_r8
          c0 = 3.E-3_r8
          ke = 3.0E-6_r8
          IF(masterproc) WRITE(6,*) 'ZM: Found spectral dycore at default resolution'
       ENDIF
    ENDIF

    tau = 4800.0_r8
    c0  =  0.5E-3_r8
    ke  = 1.0E-6_r8

  END SUBROUTINE zm_convi


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
    REAL(r8)  epsqs
    REAL(r8)  hlatv
    REAL(r8)  hlatf
    REAL(r8)  rgasv
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
       tr    = ABS(real(itype,kind=r8))
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




  !===============================================================================
  SUBROUTINE cldwat_fice(pcols,pver ,t, fice, fsnow)
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
    INTEGER,  INTENT(in)  :: pcols                 ! number of active columns
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in)  :: t(pcols,pver)        ! temperature

    REAL(r8), INTENT(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    REAL(r8), INTENT(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

    ! Local variables
    INTEGER :: i,k                                   ! loop indexes

    !-----------------------------------------------------------------------

    ! Define fractional amount of cloud that is ice
    DO k=1,pver
       DO i=1,pcols

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



  SUBROUTINE aqsat(t       ,q    ,p       ,es      ,qs        ,ii      , &
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
    REAL(r8), INTENT(inout) :: q(ii,kk)          ! Temperature
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

          !IF(qs(i,k) <= 1.0e-08_r8  )      qs(i,k)=1.0e-08_r8
          !IF(q(i,k)   >  qs(i,k))          q(i,k)=qs(i,k)

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

  REAL(r8) FUNCTION estblf( td )
    !
    ! Saturation vapor pressure table lookup
    !
    REAL(r8), INTENT(in) :: td         ! Temperature for saturation lookup
    !
    REAL(r8) :: e       ! intermediate variable for es look-up
    REAL(r8) :: tmin       ! min temperature (K) for table
    REAL(r8) :: tmax       ! max temperature (K) for table

    REAL(r8) :: ai
    INTEGER  :: i

    tmin=tmn  
    tmax=tmx  

    !
    e = MAX(MIN(td,tmax),tmin)   ! partial pressure
    i = INT(e-tmin)+1
    ai = AINT(e-tmin)
    estblf = (tmin+ai-e+1.0_r8)* &
         estbl(i)-(tmin+ai-e)* &
         estbl(i+1)
  END FUNCTION estblf


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
  !end module geopotential

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



  SUBROUTINE convtran(pcols,pver, &
       doconvtran,q       ,ncnst   ,mu      ,md      , &
       du      ,eu      ,ed      ,dp      ,dsubcld , &
       jt      ,mx      ,ideep   ,il1g    ,il2g    , &
       fracis  ,dqdt      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Convective transport of trace species
    !
    ! Mixing ratios may be with respect to either dry or moist air
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------
    ! use shr_kind_mod, only: r8 => shr_kind_r8
    ! use constituents,    only: cnst_get_type_byind
    ! use ppgrid
    ! use abortutils, only: endrun

    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    INTEGER, INTENT(in) :: ncnst                 ! number of tracers to transport
    LOGICAL  , INTENT(in) :: doconvtran(ncnst)     ! flag for doing convective transport
    REAL(r8), INTENT(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
    REAL(r8), INTENT(in) :: mu(pcols,pver)       ! Mass flux up
    REAL(r8), INTENT(in) :: md(pcols,pver)       ! Mass flux down
    REAL(r8), INTENT(in) :: du(pcols,pver)       ! Mass detraining from updraft
    REAL(r8), INTENT(in) :: eu(pcols,pver)       ! Mass entraining from updraft
    REAL(r8), INTENT(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
    REAL(r8), INTENT(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
    REAL(r8), INTENT(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
    REAL(r8), INTENT(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

    INTEGER, INTENT(in) :: jt(pcols)         ! Index of cloud top for each column
    INTEGER, INTENT(in) :: mx(pcols)         ! Index of cloud top for each column
    INTEGER, INTENT(in) :: ideep(pcols)      ! Gathering array
    INTEGER, INTENT(in) :: il1g              ! Gathered min lon indices over which to operate
    INTEGER, INTENT(in) :: il2g              ! Gathered max lon indices over which to operate

    REAL(r8) :: dpdry(pcols,pver)       ! Delta pressure between interfaces


    ! input/output

    REAL(r8), INTENT(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

    !--------------------------Local Variables------------------------------

    INTEGER i                 ! Work index
    INTEGER k                 ! Work index
    INTEGER kbm               ! Highest altitude index of cloud base
    INTEGER kk                ! Work index
    INTEGER kkp1              ! Work index
    INTEGER km1               ! Work index
    INTEGER kp1               ! Work index
    INTEGER ktm               ! Highest altitude index of cloud top
    INTEGER m                 ! Work index

    REAL(r8) cabv                 ! Mix ratio of constituent above
    REAL(r8) cbel                 ! Mix ratio of constituent below
    REAL(r8) cdifr                ! Normalized diff between cabv and cbel
    REAL(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
    REAL(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
    REAL(r8) const(pcols,pver)    ! Gathered tracer array
    REAL(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
    REAL(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
    REAL(r8) dcondt(pcols,pver)   ! Gathered tend array
    REAL(r8) small                ! A small number
    REAL(r8) mbsth                ! Threshold for mass fluxes
    REAL(r8) mupdudp              ! A work variable
    REAL(r8) minc                 ! A work variable
    REAL(r8) maxc                 ! A work variable
    REAL(r8) fluxin               ! A work variable
    REAL(r8) fluxout              ! A work variable
    REAL(r8) netflux              ! A work variable

    REAL(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
    REAL(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
    REAL(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
    REAL(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces
    !-----------------------------------------------------------------------
    !
    !call t_startf ('convtran')




    small = 1.0e-36_r8
    ! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
    mbsth = 1.0e-15_r8
    dpdry=0.0_r8

    ! Find the highest level top and bottom levels of convection
    ktm = pver
    kbm = pver
    DO i = il1g, il2g
       ktm = MIN(ktm,jt(i))
       kbm = MIN(kbm,mx(i))
    END DO

    ! Loop ever each constituent
    DO m = 2, ncnst
       IF (doconvtran(m)) THEN

          IF (cnst_get_type_byind(m,ncnst).EQ.'dry') THEN
             !if (  .not. present(dpdry) ) then
             !   write(6,*)'convtran was asked to do dry tracers but was called without dpdry argument'
             !   call endrun('convtran')
             !endif
             DO k = 1,pver
                DO i =il1g,il2g
                   dptmp(i,k) = dpdry(i,k)
                   dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                   eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                   edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
                END DO
             END DO
          ELSE
             DO k = 1,pver
                DO i =il1g,il2g
                   dptmp(i,k) = dp(i,k)
                   dutmp(i,k) = du(i,k)
                   eutmp(i,k) = eu(i,k)
                   edtmp(i,k) = ed(i,k)
                END DO
             END DO
          ENDIF
          !        dptmp = dp

          ! Gather up the constituent and set tend to zero
          DO k = 1,pver
             DO i =il1g,il2g
                const(i,k) = q(ideep(i),k,m)
                fisg(i,k) = fracis(ideep(i),k,m)
             END DO
          END DO

          ! From now on work only with gathered data

          ! Interpolate environment tracer values to interfaces
          DO k = 1,pver
             km1 = MAX(1,k-1)
             DO i = il1g, il2g
                minc = MIN(const(i,km1),const(i,k))
                maxc = MAX(const(i,km1),const(i,k))
                IF (minc < 0) THEN
                   cdifr = 0.0_r8
                ELSE
                   cdifr = ABS(const(i,k)-const(i,km1))/MAX(maxc,small)
                ENDIF

                ! If the two layers differ significantly use a geometric averaging
                ! procedure
                IF (cdifr > 1.0E-6_r8) THEN
                   cabv = MAX(const(i,km1),maxc*1.0e-12_r8)
                   cbel = MAX(const(i,k),maxc*1.0e-12_r8)
                   chat(i,k) = LOG(cabv/cbel)/(cabv-cbel)*cabv*cbel

                ELSE             ! Small diff, so just arithmetic mean
                   chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
                END IF

                ! Provisional up and down draft values
                conu(i,k) = chat(i,k)
                cond(i,k) = chat(i,k)

                !              provisional tends
                dcondt(i,k) = 0.0_r8

             END DO
          END DO

          ! Do levels adjacent to top and bottom
          k = 2
          km1 = 1
          kk = pver
          DO i = il1g,il2g
             mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
             IF (mupdudp > mbsth) THEN
                conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
             ENDIF
             IF (md(i,k) < -mbsth) THEN
                cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
             ENDIF
          END DO

          ! Updraft from bottom to top
          DO kk = pver-1,1,-1
             kkp1 = MIN(pver,kk+1)
             DO i = il1g,il2g
                mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
                IF (mupdudp > mbsth) THEN
                   conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                        const(i,kk)*dptmp(i,kk) )/mupdudp
                ENDIF
             END DO
          END DO

          ! Downdraft from top to bottom
          DO k = 3,pver
             km1 = MAX(1,k-1)
             DO i = il1g,il2g
                IF (md(i,k) < -mbsth) THEN
                   cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                        *dptmp(i,km1) )/md(i,k)
                ENDIF
             END DO
          END DO


          DO k = ktm,pver
             km1 = MAX(1,k-1)
             kp1 = MIN(pver,k+1)
             DO i = il1g,il2g

                ! version 1 hard to check for roundoff errors
                !               dcondt(i,k) =
                !     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
                !     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
                !     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
                !     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
                !     $                   )/dp(i,k)

                ! version 2 hard to limit fluxes
                !               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
                !     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
                !               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
                !     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

                ! version 3 limit fluxes outside convection to mass in appropriate layer
                ! these limiters are probably only safe for positive definite quantitities
                ! it assumes that mu and md already satify a courant number limit of 1
                fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*MIN(chat(i,k),const(i,km1)) &
                     -(md(i,k)  *cond(i,k) + md(i,kp1)*MIN(chat(i,kp1),const(i,kp1)))
                fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*MIN(chat(i,kp1),const(i,k)) &
                     -(md(i,kp1)*cond(i,kp1) + md(i,k)*MIN(chat(i,k),const(i,k)))

                netflux = fluxin - fluxout
                IF (ABS(netflux) < MAX(fluxin,fluxout)*1.0e-12_r8) THEN
                   netflux = 0.0_r8
                ENDIF
                dcondt(i,k) = netflux/dptmp(i,k)
             END DO
          END DO
          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !
          !DIR$ NOINTERCHANGE
          DO k = kbm,pver
             km1 = MAX(1,k-1)
             DO i = il1g,il2g
                IF (k == mx(i)) THEN

                   ! version 1
                   !                  dcondt(i,k) = (1./dsubcld(i))*
                   !     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
                   !     $               -md(i,k)*(cond(i,k)-chat(i,k))
                   !     $              )

                   ! version 2
                   !                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
                   !                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
                   ! version 3
                   fluxin =  mu(i,k)*MIN(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                   fluxout = mu(i,k)*conu(i,k) - md(i,k)*MIN(chat(i,k),const(i,k))

                   netflux = fluxin - fluxout
                   IF (ABS(netflux) < MAX(fluxin,fluxout)*1.0e-12_r8) THEN
                      netflux = 0.0_r8
                   ENDIF
                   !                  dcondt(i,k) = netflux/dsubcld(i)
                   dcondt(i,k) = netflux/dptmp(i,k)
                ELSE IF (k > mx(i)) THEN
                   !                  dcondt(i,k) = dcondt(i,k-1)
                   dcondt(i,k) = 0.0_r8
                END IF
             END DO
          END DO

          ! Initialize to zero everywhere, then scatter tendency back to full array
          dqdt(:,:,m) = 0.0_r8
          DO k = 1,pver
             kp1 = MIN(pver,k+1)
             !DIR$ CONCURRENT
             DO i = il1g,il2g
                dqdt(ideep(i),k,m) = dcondt(i,k)
             END DO
          END DO

       END IF      ! for doconvtran

    END DO

    !call t_stopf ('convtran')
    RETURN
  END SUBROUTINE convtran

  !===============================================================================
  SUBROUTINE zm_conv_evap(pcols, pver,pverp,&
       t,pmid,pdel,q, &
       tend_s, tend_q, &
       prdprec, cldfrc, deltat,  &
       prec, snow, ntprprd, ntsnprd, flxprec, flxsnow )

    !-----------------------------------------------------------------------
    ! Compute tendencies due to evaporation of rain from ZM scheme
    !--
    ! Compute the total precipitation and snow fluxes at the surface.
    ! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
    ! in the Zhang-MacFarlane parameterization.
    ! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
    !-----------------------------------------------------------------------


    !------------------------------Arguments--------------------------------
    INTEGER,INTENT(in) :: pcols             ! number of columns and chunk index
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: pverp                 ! number of vertical levels + 1

    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: t          ! temperature (K)
    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    REAL(r8),INTENT(in), DIMENSION(pcols,pver) :: pdel       ! layer thickness (Pa)
    REAL(r8),INTENT(inout), DIMENSION(pcols,pver) :: q          ! water vapor (kg/kg)
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
    prec(:pcols) = prec(:pcols)*1000.0_r8

    ! determine saturation vapor pressure
    CALL aqsat (t    ,q      ,pmid  ,est    ,qsat    ,pcols   , &
         pcols ,pver  ,1       ,pver    )

    ! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    CALL cldwat_fice(pcols,pver ,t, fice, fsnow_conv)

    ! zero the flux integrals on the top boundary
    flxprec(:pcols,1) = 0.0_r8
    flxsnow(:pcols,1) = 0.0_r8
    evpvint(:pcols)   = 0.0_r8

    DO k = 1, pver
       DO i = 1, pcols

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
    prec(:pcols) = flxprec(:pcols,pver+1) / 1000.0_r8
    snow(:pcols) = flxsnow(:pcols,pver+1) / 1000.0_r8

    !**********************************************************
!!$    tend_s(:pcols,:)   = 0.      ! turn heating off
    !**********************************************************

  END SUBROUTINE zm_conv_evap

  !#include <misc.h>
  !#include <params.h>

  !module geopotential

  !---------------------------------------------------------------------------------
  ! Compute geopotential from temperature or
  ! compute geopotential and temperature from dry static energy.
  !
  ! The hydrostatic matrix elements must be consistent with the dynamics algorithm.
  ! The diagonal element is the itegration weight from interface k+1 to midpoint k.
  ! The offdiagonal element is the weight between interfaces.
  ! 
  ! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
  !---------------------------------------------------------------------------------

  !  use shr_kind_mod, only: r8 => shr_kind_r8
  !use ppgrid, only: pcols, pver, pverp
  !  use dycore, only: dycore_is

  !contains
  !===============================================================================
  SUBROUTINE geopotential_dse(                                &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , pcols   ,pver, pverp          )
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
    INTEGER, INTENT(in) :: pcols                  ! Number of longitudes
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
    DO i = 1,pcols
       zi(i,pverp) = 0.0_r8
    END DO

    ! Compute the virtual temperature, zi, zm from bottom up
    ! Note, zi(i,k) is the interface above zm(i,k)
    DO k = pver, 1, -1

       ! First set hydrostatic elements consistent with dynamics
       IF (fvdyn) THEN
          DO i = 1,pcols
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1.0_r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          END DO
       ELSE
          DO i = 1,pcols
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF

       ! Now compute tv, t, zm, zi
       DO i = 1,pcols
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
       zi     , zm     , pcols   , pver, pverp)

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
    INTEGER, INTENT(in) :: pcols                  ! Number of longitudes
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

    DO i = 1,pcols
       zi(i,pverp) = 0.0_r8
    END DO

    ! Compute zi, zm from bottom up. 
    ! Note, zi(i,k) is the interface above zm(i,k)

    DO k = pver, 1, -1

       ! First set hydrostatic elements consistent with dynamics

       IF (fvdyn) THEN
          DO i = 1,pcols
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1.0_r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          END DO
       ELSE
          DO i = 1,pcols
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF

       ! Now compute tv, zm, zi

       DO i = 1,pcols
          tvfac   = 1.0_r8 + zvir * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv * hkl(i)
       END DO
    END DO

    RETURN
  END SUBROUTINE geopotential_t
  !end module geopotential
  SUBROUTINE qneg3 (subnam  ,pcols    ,ncold   ,lver    ,lconst_beg  , &
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

    INTEGER, INTENT(in) :: pcols         ! number of atmospheric columns
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
    INTEGER indx(pcols,lver)  ! array of indices of points < qmin
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
!!       !DIR$ preferstream
       DO k=1,lver
          nval(k) = 0
!!          !DIR$ prefervector
          DO i=1,pcols
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


real(r8) function entropy(TK,p,qtot)





     real(r8), intent(in) :: p,qtot,TK
     real(r8) :: qv,qsat,e,esat,L,pref

pref = 1000.0_r8           

L = rl - (cpliq - cpwv)*(TK-tfreez)         



esat = c1*exp(c2*(TK-tfreez)/(c3+TK-tfreez))       
qsat=eps1*esat/(p-esat)                      

qv = min(qtot,qsat)                         
e = qv*p / (eps1 +qv)

entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
        L*qv/TK - qv*rh2o*log(qv/qsat)

return
end FUNCTION entropy


subroutine parcel_dilute (lchnk, ncol, pcols,pver,msg, klaunch, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)


! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.


implicit none


integer, intent(in) :: lchnk
integer, intent(in) :: ncol
integer, intent(in) :: pcols
integer, intent(in) :: pver
integer, intent(in) :: msg

integer, intent(in), dimension(pcols) :: klaunch(pcols)

real(r8), intent(in), dimension(pcols,pver) :: p
real(r8), intent(in), dimension(pcols,pver) :: t
real(r8), intent(in), dimension(pcols,pver) :: q
real(r8), intent(in), dimension(pcols) :: tpert  ! PBL temperature perturbation.

real(r8), intent(inout), dimension(pcols,pver) :: tp     ! Parcel temp.
real(r8), intent(inout), dimension(pcols,pver) :: qstp   ! Parcel water vapour (sat value above lcl).
real(r8), intent(inout), dimension(pcols) :: tl    ! Actual temp of LCL.     
real(r8), intent(inout), dimension(pcols) :: pl   ! Actual pressure of LCL.         

integer, intent(inout), dimension(pcols) :: lcl ! Lifting condesation level (first model level with saturation).

real(r8), intent(out), dimension(pcols,pver) :: tpv   ! Define tpv within this routine.

!--------------------

! Have to be careful as s is also dry static energy.


! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(r8) tmix(pcols,pver)        ! Tempertaure of the entraining parcel.
real(r8) qtmix(pcols,pver)       ! Total water of the entraining parcel.
real(r8) qsmix(pcols,pver)       ! Saturated mixing ratio at the tmix.
real(r8) smix(pcols,pver)        ! Entropy of the entraining parcel.
real(r8) xsh2o(pcols,pver)       ! Precipitate lost from parcel.
real(r8) ds_xsh2o(pcols,pver)    ! Entropy change due to loss of condensate.
real(r8) ds_freeze(pcols,pver)   ! Entropy change sue to freezing of precip.

real(r8) mp(pcols)    ! Parcel mass flux.
real(r8) qtp(pcols)   ! Parcel total water.
real(r8) sp(pcols)    ! Parcel entropy.

real(r8) sp0(pcols)    ! Parcel launch entropy.
real(r8) qtp0(pcols)   ! Parcel launch total water.
real(r8) mp0(pcols)    ! Parcel launch relative mass flux.

real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
real(r8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(r8) senv       ! Environmental entropy at each grid point.
real(r8) qtenv      ! Environmental total water "   "   ".
real(r8) penv       ! Environmental total pressure "   "   ".
real(r8) tenv       ! Environmental total temperature "   "   ".
real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(r8) dp         ! Layer thickness (center to center)
real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(r8) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(r8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(r8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat     ! Number of iterations for condensation/freezing loop.
integer i,k,ii   ! Loop counters.


!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!


nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.
dmpdz=-1.e-3_r8        ! Entrainment rate. (-ve for /m)
!dmpdpc = 3.e-2_r8   ! In cloud entrainment rate (/mb).
lwmax = 1.e-3_r8    ! Need to put formula in for this.
tscool = 0.0_r8   ! Temp at which water loading freezes in the cloud.

qtmix=0._r8
smix=0._r8

qtenv = 0._r8
senv = 0._r8
tenv = 0._r8
penv = 0._r8

qtp0 = 0._r8
sp0  = 0._r8
mp0 = 0._r8

qtp = 0._r8
sp = 0._r8
mp = 0._r8

new_q = 0._r8
new_s = 0._r8

! **** Begin loops ****


do k = pver, msg+1, -1
   do i=1,ncol 

! Initialize parcel values at launch level.

      if (k == klaunch(i)) then 
         qtp0(i) = q(i,k)   ! Parcel launch total water (assuming subsaturated) - OK????.
         sp0(i)  = entropy(t(i,k),p(i,k),qtp0(i))  ! Parcel launch entropy.
         mp0(i)  = 1._r8       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 
         smix(i,k)  = sp0(i)
         qtmix(i,k) = qtp0(i)
         tfguess = t(i,k)
         rcall = 1
         call ientropy ( &
                        rcall     , & !integer , intent(in)  :: rcall
                        i         , & !integer , intent(in)  :: icol
                        lchnk     , & !integer , intent(in)  :: lchnk
                        smix(i,k) , & !real(r8), intent(in)  :: s 
                        p(i,k)    , & !real(r8), intent(in)  :: p 
                        qtmix(i,k), & !real(r8), intent(in)  :: qt
                        tmix(i,k) , & !real(r8), intent(out) :: T
                        qsmix(i,k), & !real(r8), intent(out) :: qsat
                        tfguess     ) !real(r8), intent(in)  :: Tfg 
      end if

! Entraining levels
      
      if (k < klaunch(i)) then 

! Set environmental values for this level.                 
         
         dp = (p(i,k)-p(i,k+1)) ! In -ve mb as p decreasing with height - difference between center of layers.
         qtenv = 0.5_r8*(q(i,k)+q(i,k+1))   ! Total water of environment.       
         tenv  = 0.5_r8*(t(i,k)+t(i,k+1)) 
         penv  = 0.5_r8*(p(i,k)+p(i,k+1))

         senv  = entropy(tenv,penv,qtenv)   ! Entropy of environment.   

! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
         dzdp = 1._r8/dpdz                  ! in m/mb
         dmpdp = dmpdz*dzdp               ! /mb Fractional entrainment

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp(i)  = sp(i)  - dmpdp*dp*senv 
         qtp(i) = qtp(i) - dmpdp*dp*qtenv 
         mp(i)  = mp(i)  - dmpdp*dp
            
! Entrain s and qt to next level.

         smix(i,k)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
         qtmix(i,k) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(i,k) used as a first guess so that it converges faster.

         tfguess = tmix(i,k+1)
         rcall = 2
         call ientropy( &
                        rcall      , &!integer , intent(in)  :: rcall
                        i          , &!integer , intent(in)  :: icol
                        lchnk      , &!integer , intent(in)  :: lchnk
                        smix(i,k)  , &!real(r8), intent(in)  :: s 
                        p(i,k)     , &!real(r8), intent(in)  :: p 
                        qtmix(i,k) , &!real(r8), intent(in)  :: qt
                        tmix(i,k)  , &!real(r8), intent(out) :: T
                        qsmix(i,k) , &!real(r8), intent(out) :: qsat
                        tfguess      )!real(r8), intent(in)  :: Tfg 


!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.



         if (qsmix(i,k) <= qtmix(i,k) .and. qsmix(i,k+1) > qtmix(i,k+1)) then
            lcl(i) = k
            qxsk   = qtmix(i,k) - qsmix(i,k)
            qxskp1 = qtmix(i,k+1) - qsmix(i,k+1)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(i,k+1) - qxskp1/dqxsdp     ! pressure level of actual lcl.
            dsdp   = (smix(i,k)  - smix(i,k+1))/dp
            dqtdp  = (qtmix(i,k) - qtmix(i,k+1))/dp
            slcl   = smix(i,k+1)  +  dsdp* (pl(i)-p(i,k+1))  
            qtlcl  = qtmix(i,k+1) +  dqtdp*(pl(i)-p(i,k+1))

            tfguess = tmix(i,k)
            rcall = 3
            call ientropy (   &
                            rcall        , & !integer , intent(in)  :: rcall
                            i            , & !integer , intent(in)  :: icol
                            lchnk        , & !integer , intent(in)  :: lchnk
                            slcl         , & !real(r8), intent(in)  :: s 
                            pl(i)        , & !real(r8), intent(in)  :: p 
                            qtlcl        , & !real(r8), intent(in)  :: qt
                            tl(i)        , & !real(r8), intent(out) :: T
                            qslcl        , & !real(r8), intent(out) :: qsat
                            tfguess         )!real(r8), intent(in)  :: Tfg 








         endif

      end if 

 
   end do 
end do 


!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 
 

xsh2o = 0._r8
ds_xsh2o = 0._r8
ds_freeze = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



do k = pver, msg+1, -1
   do i=1,ncol    
      
! Initialize variables at k=klaunch

      
      if (k == klaunch(i)) then

! Set parcel values at launch level assume no liquid water.            


         tp(i,k)    = tmix(i,k)
         qstp(i,k)  = q(i,k) 
         tpv(i,k)   =  (tp(i,k) + tpert(i)) * (1._r8+1.608_r8*qstp(i,k)) / (1._r8+qstp(i,k))
         
      end if

      if (k < klaunch(i)) then
            
! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

! Iterate nit_lheat times for s,qt changes.

         do ii=0,nit_lheat-1            

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).

            xsh2o(i,k) = max (0._r8, qtmix(i,k) - qsmix(i,k) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)                     

                     
            ds_xsh2o(i,k) = ds_xsh2o(i,k+1) - cpliq * log (tmix(i,k)/tfreez) * max(0._r8,(xsh2o(i,k)-xsh2o(i,k+1)))

!
! Entropy of freezing: latice times amount of water involved divided by T.
!
 
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) == 0._r8) then 
               ds_freeze(i,k) = (latice/tmix(i,k)) * max(0._r8,qtmix(i,k)-qsmix(i,k)-xsh2o(i,k)) 
            end if
            
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) /= 0._r8) then 
               ds_freeze(i,k) = ds_freeze(i,k+1)+(latice/tmix(i,k)) * max(0._r8,(qsmix(i,k+1)-qsmix(i,k)))
            end if
            
! Adjust entropy and accordingly to sum of ds (be careful of signs).

            new_s = smix(i,k) + ds_xsh2o(i,k) + ds_freeze(i,k) 

! Adjust liquid water and accordingly to xsh2o.

            new_q = qtmix(i,k) - xsh2o(i,k)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(i,k)
            rcall =4
            call ientropy ( &
                           rcall     , &
                           i         , &
                           lchnk     , &
                           new_s     , &
                           p(i,k)    , &
                           new_q     , &
                           tmix(i,k) , &
                           qsmix(i,k),&
                           tfguess    )
            
         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water. 

         tp(i,k)    = tmix(i,k)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         if (new_q > qsmix(i,k)) then  ! Super-saturated so condensate present - reduces buoyancy. 
            qstp(i,k) = qsmix(i,k)
         else                      ! Just saturated/sub-saturated - no condensate virtual effects.    
            qstp(i,k) = new_q
         end if

         tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+1.608_r8*qstp(i,k)) / (1._r8+ new_q) 

      end if 
      
   end do 
   
end do  


return
end subroutine parcel_dilute

   SUBROUTINE ientropy ( &
                         rcall        , &!integer , intent(in)  :: rcall
                         icol         , &!integer , intent(in)  :: icol
                         lchnk        , &!integer , intent(in)  :: lchnk
                         s            , &!real(r8), intent(in)  :: s        (J/kg)
                         p            , &!real(r8), intent(in)  :: p       p(mb)
                         qt           , &!real(r8), intent(in)  :: qt     (kg/kg)
                         T            , &!real(r8), intent(out) :: T      T(K)
                         qsat         , &!real(r8), intent(out) :: qsat
                         Tfg            )!real(r8), intent(in)  :: Tfg Tfg(K)




!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 



     integer , intent(in)  :: icol
     integer , intent(in)  :: lchnk
     integer , intent(in)  :: rcall
     real(r8), intent(in)  :: s 
     real(r8), intent(in)  :: p 
     real(r8), intent(in)  :: Tfg 
     real(r8), intent(in)  :: qt
     real(r8), intent(out) :: qsat
     real(r8), intent(out) :: T
     real(r8) :: qv,Ts,dTs,fs1,fs2,esat     
     real(r8) :: pref,L,e
     real(r8) :: this_lat,this_lon,esft
     integer :: LOOPMAX,i

LOOPMAX = 500      !* max number of iteration loops                


pref = 1000.0_r8           
!pref = ps           



Ts = Tfg                  

converge: do i=0, LOOPMAX

   L = rl - (cpliq - cpwv)*(Ts-tfreez) 

   esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez)) ! Bolton (eq. 10)
   qsat = eps1*esat/(p-esat)     
   qv = min(qt,qsat) 
   e = qv*p / (eps1 +qv)  ! Bolton (eq. 16)
   fs1 = (cpres + qt*cpliq)*log( Ts/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/Ts - qv*rh2o*log(qv/qsat) - s
   
   L = rl - (cpliq - cpwv)*(Ts-1._r8-tfreez)         

   esat = c1*exp(c2*(Ts-1._r8-tfreez)/(c3+Ts-1._r8-tfreez))
   qsat = eps1*esat/(p-esat)  
   qv = min(qt,qsat) 
   e = qv*p / (eps1 +qv)
   fs2 = (cpres + qt*cpliq)*log( (Ts-1._r8)/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/(Ts-1._r8) - qv*rh2o*log(qv/qsat) - s 
   
   dTs = fs1/(fs2 - fs1)
   Ts  = Ts+dTs
   if (abs(dTs).lt.0.001_r8) exit converge
   if (i .eq. LOOPMAX - 1) then

      this_lat = 0.
      this_lon = 0.
      write(0,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      !call wrf_message(0)
      write(0,100) 'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
       ' lat: ',this_lat,' lon: ',this_lon, &
       ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
       ' qsat(g/kg) = ', 1000._r8*qsat,', s(J/kg) = ',s
      !call wrf_message(0)
      call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
   end if
enddo converge

! Replace call to satmixutils.


esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez))
qsat=eps1*esat/(p-esat)

qv = min(qt,qsat)                             
T = Ts 

 100    format (A,I1,I4,I4,7(A,e12.5))

return
end SUBROUTINE ientropy
  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    REAL(KIND=r8), INTENT(IN) :: t

    IF (t <= tmelt) THEN
       es5 = EXP(ae(2)-be(2)/t)
    ELSE
       es5 = EXP(ae(1)-be(1)/t)
    END IF
  END FUNCTION es5

END MODULE Cu_ZhangMcFarlane

!PROGRAM MAIN
!  USE Cu_ZhangMcFarlane
!
!END PROGRAM MAIN
