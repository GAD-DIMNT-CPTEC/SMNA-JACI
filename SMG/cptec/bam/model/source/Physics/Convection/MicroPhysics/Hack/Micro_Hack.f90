MODULE Micro_Hack
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
  REAL(r8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air ~ kg/kmole
  REAL(r8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(r8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_G      = 9.80616_r8      ! acceleration of gravity ~ m/s^2
  REAL(r8),PARAMETER :: SHR_CONST_TKFRZ  = 273.16_r8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(r8),PARAMETER :: SHR_CONST_LATVAP = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_LATICE = 3.337e5_r8      ! latent heat of fusion ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_RHOFW  = 1.000e3_r8      ! density of fresh water ~ kg/m^3
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_PSTD   = 101325.0_r8     ! standard pressure ~ pascals

  REAL(r8),PARAMETER :: SHR_CONST_RHODAIR=SHR_CONST_PSTD/ &
       (SHR_CONST_RDAIR*SHR_CONST_TKFRZ)         ! density of dry air at STP   ~ kg/m^3

  ! Constants for air

  REAL(r8), PUBLIC, PARAMETER :: rair = shr_const_rdair    ! Gas constant for dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: cappa = rair/cpair        ! R/Cp
  REAL(r8), PUBLIC, PARAMETER :: rhodair = shr_const_rhodair ! density of dry air at STP (kg/m3)

  ! Constants for Earth
  REAL(r8), PUBLIC, PARAMETER :: gravit = shr_const_g      ! gravitational acceleration
  ! Constants for water

  REAL(r8), PUBLIC, PARAMETER :: tmelt = shr_const_tkfrz   ! Freezing point of water
  REAL(r8), PUBLIC, PARAMETER :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights 
  REAL(r8), PUBLIC, PARAMETER :: latvap = shr_const_latvap ! Latent heat of vaporization
  REAL(r8), PUBLIC, PARAMETER :: latice = shr_const_latice ! Latent heat of fusion
  REAL(r8), PUBLIC, PARAMETER :: rhoh2o = shr_const_rhofw  ! Density of liquid water (STP)

  !------------cloud_fraction ----------------
  !
  ! Private data
  !
  REAL(KIND=r8) :: rhminl                ! minimum rh for low stable clouds
  REAL(KIND=r8) :: rhminh                ! minimum rh for high stable clouds
  REAL(KIND=r8) :: sh1,sh2               ! parameters for shallow convection cloud fraction
  REAL(KIND=r8) :: dp1,dp2               ! parameters for deep convection cloud fraction
  REAL(KIND=r8) :: premit                ! top pressure bound for mid level cloud

  INTEGER        :: k700                  ! cldconst -- model level nearest 700 mb


  PUBLIC :: cldfrc_init ! Inititialization of cloud_fraction run-time parameters
  PUBLIC :: cldfrc      !  Computation of cloud fraction

  !------------end cloud_fraction-------------
  ! 
  LOGICAL :: PERGRO=.FALSE.
  LOGICAL :: OLDLIQSED=.FALSE.
  LOGICAL :: UNICOSMP=.FALSE.
  LOGICAL :: HEAVYNEW=.TRUE.
  !LOGICAL :: cnst_need_pdeldry=.FALSE.
  !------------wv_saturation-------------
  !
  ! Data
  !
  INTEGER, PARAMETER :: plenest = 250! length of saturation vapor pressure table

  REAL(r8)           :: estbl(plenest)      ! table values of saturation vapor pressure
  REAL(r8),PARAMETER :: tmn  = 173.16_r8          ! Minimum temperature entry in table
  REAL(r8),PARAMETER :: tmx  = 375.16_r8          ! Maximum temperature entry in table
  REAL(r8),PARAMETER :: trice  =  20.00_r8         ! Trans range from es over h2o to es over ice
  REAL(r8),PARAMETER :: tmin=tmn       ! min temperature (K) for table
  REAL(r8),PARAMETER :: tmax= tmx      ! max temperature (K) for table
  LOGICAL ,PARAMETER :: icephs=.TRUE.  ! false => saturation vapor press over water only
  INTEGER,PARAMETER  :: iterp =50             ! #iterations for precipitation calculation

  !------------wv_saturation end -------------

  !------------module pkg_cld_sediment-------------
  REAL (r8), PARAMETER :: mxsedfac   = 0.99_r8       ! maximum sedimentation flux factor
  REAL (r8), PARAMETER :: vland  = 1.5_r8            ! liquid fall velocity over land  (cm/s)
  REAL (r8), PARAMETER :: vocean = 2.8_r8            ! liquid fall velocity over ocean (cm/s)
  LOGICAL,    PARAMETER :: stokes = .TRUE.            ! use Stokes velocity instead of McFarquhar and Heymsfield
  ! parameter for modified McFarquhar and Heymsfield
  REAL (r8), PARAMETER :: vice_small = 1.0_r8         ! ice fall velocity for small concentration (cm/s)

  ! parameters for Stokes velocity
  REAL (r8), PARAMETER :: eta =  1.7e-5_r8           ! viscosity of air (kg m / s)
  REAL (r8), PARAMETER :: r40 =  40.0_r8             !  40 micron radius
  REAL (r8), PARAMETER :: r400= 400.0_r8             ! 400 micron radius
  REAL (r8), PARAMETER :: v40 = (2.0_r8/9.0_r8) * rhoh2o * gravit/eta * r40**2 * 1.0e-12_r8  
  ! Stokes fall velocity of 40 micron sphere (m/s)
  REAL (r8), PARAMETER :: v400= 1.00_r8              ! fall velocity of 400 micron sphere (m/s)
  REAL (r8), PARAMETER :: vslope = (v400 - v40)/(r400 -r40) ! linear slope for large particles m/s/micron

  !------------end module pkg_cld_sediment-------------

  !------------   module cldwat-------------

  REAL(r8), PRIVATE, PARAMETER :: tmax_fice = tmelt - 10.0_r8       ! max temperature for cloud ice formation
!!$   real(r8), private, parameter :: tmax_fice = tmelt          ! max temperature for cloud ice formation
!!$   real(r8), private, parameter :: tmin_fice = tmax_fice - 20.! min temperature for cloud ice formation
  REAL(r8), PRIVATE, PARAMETER :: tmin_fice = tmax_fice - 30.0_r8   ! min temperature for cloud ice formation
  !  pjr
  REAL(r8), PRIVATE, PARAMETER :: tmax_fsnow = tmelt            ! max temperature for transition to convective snow
  REAL(r8), PRIVATE, PARAMETER :: tmin_fsnow = tmelt-5.0_r8         ! min temperature for transition to convective snow

  REAL(r8),PUBLIC ::  icritc   =   5.e-6_r8           ! threshold for autoconversion of cold ice
  REAL(r8),PUBLIC ::  icritw               ! threshold for autoconversion of warm ice
!!$   real(r8),public,parameter::  conke  = 1.e-6    ! tunable constant for evaporation of precip
!!$   real(r8),public,parameter::  conke  =  2.e-6    ! tunable constant for evaporation of precip
  REAL(r8),PUBLIC ::  conke                          ! tunable constant for evaporation of precip
  REAL(r8),PUBLIC ::  conke_land                     ! tunable constant for evaporation of precip

  REAL(r8), PRIVATE :: rhonot   ! air density at surface

  REAL(r8), PRIVATE :: t0       ! Freezing temperature
  REAL(r8), PRIVATE :: cldmin   ! assumed minimum cloud amount
  REAL(r8), PRIVATE :: small    ! small number compared to unity
  REAL(r8), PRIVATE :: c        ! constant for graupel like snow cm**(1-d)/s
  REAL(r8), PRIVATE :: d        ! constant for graupel like snow
  REAL(r8), PRIVATE :: esi      ! collection efficient for ice by snow
  REAL(r8), PRIVATE :: esw      ! collection efficient for water by snow
  REAL(r8), PRIVATE :: nos      ! particles snow / cm**4
  REAL(r8), PRIVATE :: pi       ! Mathematical constant
  !real(r8), private :: gravit   ! Gravitational acceleration at surface  
  REAL(r8), PRIVATE :: rh2o  =SHR_CONST_RWV   !! Gas constant for water vapor
  REAL(r8), PRIVATE :: prhonos
  REAL(r8), PRIVATE :: thrpd    ! numerical three added to d  
  REAL(r8), PRIVATE :: gam3pd   ! gamma function on (3+d)
  REAL(r8), PRIVATE :: gam4pd   ! gamma function on (4+d)
  REAL(r8), PRIVATE :: rhoi     ! ice density
  REAL(r8), PRIVATE :: rhos     ! snow density
  REAL(r8), PRIVATE :: rhow     ! water density
  REAL(r8), PRIVATE :: mcon01   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon02   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon03   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon04   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon05   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon06   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon07   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: mcon08   ! constants used in cloud microphysics
  REAL(r8), PRIVATE :: ttrice=trice
  REAL(r8), PRIVATE :: cp    =cpair
  REAL(r8), PRIVATE :: hlatv  = latvap
  REAL(r8), PRIVATE :: hlatf  = latice
  REAL(r8), PRIVATE :: rgasv  = SHR_CONST_RWV    ! Gas constant for water vapor
  REAL(r8), PRIVATE :: zvir = SHR_CONST_RWV/rair - 1          ! rh2o/rair - 1

  INTEGER  , PRIVATE ::  k1mb  =1  ! index of the eta level near 1 mb
  REAL(r8) :: pcf(6)     ! polynomial coeffs -> es transition water to ice

  !------------  END module cldwat-------------

  !------------  module physics_types-------------
  !-------------------------------------------------------------------------------
  TYPE physics_state
     INTEGER , POINTER :: ncol    (:) ! number of active columns
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
     INTEGER , POINTER  :: COUNT   (:)              ! count of values with significant energy or water imbalances     
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

     LOGICAL , POINTER  :: ls (:)                 ! true if dsdt is returned
     LOGICAL , POINTER  :: lu (:)                 ! true if dudt is returned
     LOGICAL , POINTER  :: lv (:)                 ! true if dvdt is returned
     LOGICAL  , POINTER :: lq  (:,:)  !(ppcnst)       ! true if dqdt() is returned
     INTEGER , POINTER  :: top_level  (:)        ! top level index for which nonzero tendencies have been set
     INTEGER , POINTER  :: bot_level  (:)       ! bottom level index for which nonzero tendencies have been set
     REAL(r8), POINTER :: s(:,:,:) !(pcols,pver)! heating rate (J/kg/s)
     REAL(r8), POINTER :: u(:,:,:) !(pcols,pver)! u momentum tendency (m/s/s)
     REAL(r8), POINTER :: v(:,:,:) !(pcols,pver)! v momentum tendency (m/s/s)
     REAL(r8), POINTER :: q(:,:,:,:)!(pcols,pver,ppcnst)                 ! consituent tendencies (kg/kg/s)
     REAL(r8), POINTER :: qt(:,:,:) !(pcols,pver)! heating rate (J/kg/s)

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



  !REAL(r8), POINTER :: state_qcwat (:,:,:)!(1:pcols,1:pver)  ! cloud water old q
  !REAL(r8), POINTER :: state_tcwat (:,:,:)!(1:pcols,1:pver)  ! cloud water old temperature
  !REAL(r8), POINTER :: state_lcwat (:,:,:)!(1:pcols,1:pver)  ! cloud liquid water old q
  !REAL(r8), POINTER :: state_cld   (:,:,:)!(1:pcols,1:pver)  ! cloud fraction
  !REAL(r8), POINTER :: state_qme   (:,:,:)!(1:pcols,1:pver) 
  !REAL(r8), POINTER :: state_prain (:,:,:)!(1:pcols,1:pver) 
  !REAL(r8), POINTER :: state_nevapr(:,:,:)!(1:pcols,1:pver) 
  !REAL(r8), POINTER :: state_rel   (:,:,:)!(1:pcols,1:pver)   ! liquid effective drop radius (microns)
  !REAL(r8), POINTER :: state_rei   (:,:,:)!(1:pcols,1:pver)   ! ice effective drop size (microns)

  !------------  module physics_types-------------

  REAL(r8) :: qmin(3)
  INTEGER :: ixcldliq=2
  INTEGER :: ixcldice=3


  PUBLIC :: Init_Micro_Hack
  PUBLIC :: RunMicro_Hack
  PUBLIC :: DestroyMicro_Hack

CONTAINS
  SUBROUTINE Init_Micro_Hack(kMax,jMax,ibMax,jbMax,ppcnst,a_hybr,b_hybr)
    IMPLICIT NONE
    INTEGER, INTENT(in   ) :: kMax                  ! number of vertical levels
    INTEGER, INTENT(IN   ) :: jMax
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: ppcnst
    REAL(KIND=r8), INTENT(IN   ) :: a_hybr   (kMax+1)
    REAL(KIND=r8), INTENT(IN   ) :: b_hybr   (kMax+1)
    REAL(KIND=r8) :: hypi (kMax+1)
    INTEGER :: k
    DO k=1,kMax+1
!      hypi(k) =  si_in    (kMax+2-k)
!      SB   changed to hybrid (already top to bottom)
       hypi(k) =  a_hybr(k) / 1.0e5_r8 + b_hybr(k)
    END DO

    qmin=1e-12

    ALLOCATE(state%ncol     (jbMax))                  ;state%ncol     =0     ! number of active columns
    ALLOCATE(state%lat      (ibMax,jbMax))            ;state%lat      =0.0_r8!(pcols)     ! latitude (radians)
    ALLOCATE(state%lon      (ibMax,jbMax))            ;state%lon      =0.0_r8!(pcols)     ! longitude (radians)
    ALLOCATE(state%ps            (ibMax,jbMax))            ;state%ps       =0.0_r8!(pcols)     ! surface pressure
!    ALLOCATE(state%psdry    (ibMax,jbMax))            ;state%psdry    =0.0_r8!(pcols)     ! dry surface pressure
    ALLOCATE(state%phis     (ibMax,jbMax))            ;state%phis     =0.0_r8!(pcols)     ! surface geopotential
    ALLOCATE(state%t            (ibMax,kMax,jbMax))       ;state%t        =0.0_r8!(pcols,pver)! temperature (K)
    ALLOCATE(state%u            (ibMax,kMax,jbMax))       ;state%u        =0.0_r8!(pcols,pver)! zonal wind (m/s)
    ALLOCATE(state%v            (ibMax,kMax,jbMax))       ;state%v        =0.0_r8!(pcols,pver)! meridional wind (m/s)
    ALLOCATE(state%s            (ibMax,kMax,jbMax))       ;state%s        =0.0_r8!(pcols,pver)! dry static energy
    ALLOCATE(state%omega    (ibMax,kMax,jbMax))       ;state%omega    =0.0_r8!(pcols,pver)! vertical pressure velocity (Pa/s) 
    ALLOCATE(state%pmid     (ibMax,kMax,jbMax))       ;state%pmid     =0.0_r8!(pcols,pver)! midpoint pressure (Pa) 
!    ALLOCATE(state%pmiddry  (ibMax,kMax,jbMax))       ;state%pmiddry  =0.0_r8!(pcols,pver)! midpoint pressure dry (Pa) 
    ALLOCATE(state%pdel     (ibMax,kMax,jbMax))       ;state%pdel     =0.0_r8!(pcols,pver)! layer thickness (Pa)
!    ALLOCATE(state%pdeldry  (ibMax,kMax,jbMax))       ;state%pdeldry  =0.0_r8!(pcols,pver)! layer thickness dry (Pa)
    ALLOCATE(state%rpdel    (ibMax,kMax,jbMax))       ;state%rpdel    =0.0_r8!(pcols,pver)! reciprocal of layer thickness (Pa)
!    ALLOCATE(state%rpdeldry (ibMax,kMax,jbMax))       ;state%rpdeldry =0.0_r8!(pcols,pver)! recipricol layer thickness dry (Pa)
    ALLOCATE(state%lnpmid   (ibMax,kMax,jbMax))       ;state%lnpmid   =0.0_r8! (pcols,pver)! ln(pmid)
!    ALLOCATE(state%lnpmiddry(ibMax,kMax,jbMax))       ;state%lnpmiddry=0.0_r8!(pcols,pver)! log midpoint pressure dry (Pa) 
!    ALLOCATE(state%exner    (ibMax,kMax,jbMax))       ;state%exner    =0.0_r8!(pcols,pver)! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
    ALLOCATE(state%zm             (ibMax,kMax,jbMax))       ;state%zm       =0.0_r8!(pcols,pver)! geopotential height above surface at midpoints (m)
    ALLOCATE(state%q         (ibMax,kMax,jbMax,ppcnst));state%q        =0.0_r8!(pcols,pver,ppcnst)! constituent mixing ratio 
    !(kg/kg moist or dry air depending on type)
    ALLOCATE(state%pint     (ibMax,kMax+1,jbMax))     ;state%pint     =0.0_r8!(pcols,pver+1)! interface pressure (Pa)
!    ALLOCATE(state%pintdry  (ibMax,kMax+1,jbMax))     ;state%pintdry  =0.0_r8!(pcols,pver+1)! interface pressure dry (Pa) 
    ALLOCATE(state%lnpint   (ibMax,kMax+1,jbMax))     ;state%lnpint   =0.0_r8!(pcols,pver+1)! ln(pint)
!    ALLOCATE(state%lnpintdry(ibMax,kMax+1,jbMax))     ;state%lnpintdry=0.0_r8!(pcols,pver+1)! log interface pressure dry (Pa) 
    ALLOCATE(state%zi             (ibMax,kMax+1,jbMax))     ;state%zi       =0.0_r8!(pcols,pver+1)! geopotential height above surface at interfaces (m)
    ALLOCATE(state%te_ini   (ibMax,jbMax))            ;state%te_ini   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of initial state
    ALLOCATE(state%te_cur   (ibMax,jbMax))            ;state%te_cur   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of current state
    ALLOCATE(state%tw_ini   (ibMax,jbMax))            ;state%tw_ini   =0.0_r8!(pcols)  ! vertically integrated total water of initial state
    ALLOCATE(state%tw_cur   (ibMax,jbMax))            ;state%tw_cur   =0.0_r8!(pcols)  ! vertically integrated total water of new state
    ALLOCATE(state%count    (jbMax))                  ;state%count    =0! count of values with significant energy or water imbalances          

    ALLOCATE(state1%ncol     (jbMax))                  ;state1%ncol     =0     ! number of active columns
    ALLOCATE(state1%lat      (ibMax,jbMax))            ;state1%lat      =0.0_r8!(pcols)     ! latitude (radians)
    ALLOCATE(state1%lon      (ibMax,jbMax))            ;state1%lon      =0.0_r8!(pcols)     ! longitude (radians)
    ALLOCATE(state1%ps             (ibMax,jbMax))            ;state1%ps       =0.0_r8!(pcols)     ! surface pressure
!    ALLOCATE(state1%psdry    (ibMax,jbMax))            ;state1%psdry    =0.0_r8!(pcols)     ! dry surface pressure
    ALLOCATE(state1%phis     (ibMax,jbMax))            ;state1%phis     =0.0_r8!(pcols)     ! surface geopotential
    ALLOCATE(state1%t             (ibMax,kMax,jbMax))       ;state1%t        =0.0_r8!(pcols,pver)! temperature (K)
    ALLOCATE(state1%u             (ibMax,kMax,jbMax))       ;state1%u        =0.0_r8!(pcols,pver)! zonal wind (m/s)
    ALLOCATE(state1%v             (ibMax,kMax,jbMax))       ;state1%v        =0.0_r8!(pcols,pver)! meridional wind (m/s)
    ALLOCATE(state1%s             (ibMax,kMax,jbMax))       ;state1%s        =0.0_r8!(pcols,pver)! dry static energy
    ALLOCATE(state1%omega    (ibMax,kMax,jbMax))       ;state1%omega    =0.0_r8!(pcols,pver)! vertical pressure velocity (Pa/s) 
    ALLOCATE(state1%pmid     (ibMax,kMax,jbMax))       ;state1%pmid     =0.0_r8!(pcols,pver)! midpoint pressure (Pa) 
!    ALLOCATE(state1%pmiddry  (ibMax,kMax,jbMax))       ;state1%pmiddry  =0.0_r8!(pcols,pver)! midpoint pressure dry (Pa) 
    ALLOCATE(state1%pdel     (ibMax,kMax,jbMax))       ;state1%pdel     =0.0_r8!(pcols,pver)! layer thickness (Pa)
!    ALLOCATE(state1%pdeldry  (ibMax,kMax,jbMax))       ;state1%pdeldry  =0.0_r8!(pcols,pver)! layer thickness dry (Pa)
    ALLOCATE(state1%rpdel    (ibMax,kMax,jbMax))       ;state1%rpdel    =0.0_r8!(pcols,pver)! reciprocal of layer thickness (Pa)
!    ALLOCATE(state1%rpdeldry (ibMax,kMax,jbMax))       ;state1%rpdeldry =0.0_r8!(pcols,pver)! recipricol layer thickness dry (Pa)
    ALLOCATE(state1%lnpmid   (ibMax,kMax,jbMax))       ;state1%lnpmid   =0.0_r8! (pcols,pver)! ln(pmid)
!    ALLOCATE(state1%lnpmiddry(ibMax,kMax,jbMax))       ;state1%lnpmiddry=0.0_r8!(pcols,pver)! log midpoint pressure dry (Pa) 
!    ALLOCATE(state1%exner    (ibMax,kMax,jbMax))       ;state1%exner    =0.0_r8!(pcols,pver)! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
    ALLOCATE(state1%zm             (ibMax,kMax,jbMax))       ;state1%zm       =0.0_r8!(pcols,pver)! geopotential height above surface at midpoints (m)
    ALLOCATE(state1%q        (ibMax,kMax,jbMax,ppcnst));state1%q        =0.0_r8!(pcols,pver,ppcnst)! constituent mixing ratio 
    !(kg/kg moist or dry air depending on type)
    ALLOCATE(state1%pint     (ibMax,kMax+1,jbMax))     ;state1%pint     =0.0_r8!(pcols,pver+1)! interface pressure (Pa)
!    ALLOCATE(state1%pintdry  (ibMax,kMax+1,jbMax))     ;state1%pintdry  =0.0_r8!(pcols,pver+1)! interface pressure dry (Pa) 
    ALLOCATE(state1%lnpint   (ibMax,kMax+1,jbMax))     ;state1%lnpint   =0.0_r8!(pcols,pver+1)! ln(pint)
!    ALLOCATE(state1%lnpintdry(ibMax,kMax+1,jbMax))     ;state1%lnpintdry=0.0_r8!(pcols,pver+1)! log interface pressure dry (Pa) 
    ALLOCATE(state1%zi             (ibMax,kMax+1,jbMax))     ;state1%zi       =0.0_r8!(pcols,pver+1)! geopotential height above surface at interfaces (m)
    ALLOCATE(state1%te_ini   (ibMax,jbMax))            ;state1%te_ini   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of initial state
    ALLOCATE(state1%te_cur   (ibMax,jbMax))            ;state1%te_cur   =0.0_r8!(pcols)  ! vertically integrated total (kinetic + static) energy of current state
    ALLOCATE(state1%tw_ini   (ibMax,jbMax))            ;state1%tw_ini   =0.0_r8!(pcols)  ! vertically integrated total water of initial state
    ALLOCATE(state1%tw_cur   (ibMax,jbMax))            ;state1%tw_cur   =0.0_r8!(pcols)  ! vertically integrated total water of new state
    ALLOCATE(state1%count    (jbMax))                  ;state1%count    =0! count of values with significant energy or water imbalances          


    ALLOCATE(tend%dtdt   (ibMax,kMax,jbMax));tend%dtdt        =0.0_r8!(pcols,pver) 
    ALLOCATE(tend%dudt   (ibMax,kMax,jbMax));tend%dudt        =0.0_r8!(pcols,pver) 
    ALLOCATE(tend%dvdt   (ibMax,kMax,jbMax));tend%dvdt        =0.0_r8!(pcols,pver) 
    ALLOCATE(tend%flx_net(ibMax,jbMax))     ;tend%flx_net=0.0_r8!(pcols      ) 
    ALLOCATE(tend%te_tnd (ibMax,jbMax))     ;tend%te_tnd =0.0_r8!(pcols)      ! cumulative boundary flux of total energy
    ALLOCATE(tend%tw_tnd (ibMax,jbMax))     ;tend%tw_tnd =0.0_r8!(pcols)      ! cumulative boundary flux of total water


    !-------------------------------------------------------------------------------
    ! This is for tendencies returned from individual parameterizations
    ALLOCATE(ptend_loc%name       (jbMax))            ! true if dsdt is returned
    ALLOCATE(ptend_loc%ls         (jbMax))            ! true if dsdt is returned
    ALLOCATE(ptend_loc%lu         (jbMax))            ! true if dudt is returned
    ALLOCATE(ptend_loc%lv         (jbMax))            ! true if dvdt is returned
    ALLOCATE(ptend_loc%top_level  (jbMax))       ! top level index for which nonzero tendencies have been set
    ALLOCATE(ptend_loc%bot_level  (jbMax))       ! bottom level index for which nonzero tendencies have been set

    ALLOCATE(ptend_loc%lq  (jbMax,ppcnst))                 ;ptend_loc%lq=.TRUE.!(ppcnst)         ! true if dqdt() is returned
    ALLOCATE(ptend_loc%s   (ibMax,kMax,jbMax))       ;ptend_loc%s =0.0_r8!(pcols,pver)! heating rate (J/kg/s)
    ALLOCATE(ptend_loc%u   (ibMax,kMax,jbMax))       ;ptend_loc%u =0.0_r8!(pcols,pver)! u momentum tendency (m/s/s)
    ALLOCATE(ptend_loc%v   (ibMax,kMax,jbMax))       ;ptend_loc%v =0.0_r8!(pcols,pver)! v momentum tendency (m/s/s)
    ALLOCATE(ptend_loc%q   (ibMax,kMax,jbMax,ppcnst));ptend_loc%q =0.0_r8!(pcols,pver,ppcnst)! consituent tendencies (kg/kg/s)
    ALLOCATE(ptend_loc%qt  (ibMax,kMax,jbMax))       ;ptend_loc%qt =0.0_r8!(pcols,pver)! heating rate (J/kg/s)

    ! boundary fluxes

    ALLOCATE(ptend_loc%hflux_srf  (ibMax,jbMax))       ;ptend_loc%hflux_srf =0.0_r8!(pcols)! net heat flux at surface (W/m2)
    ALLOCATE(ptend_loc%hflux_top  (ibMax,jbMax))       ;ptend_loc%hflux_top =0.0_r8!(pcols)! net heat flux at top of model (W/m2)
    ALLOCATE(ptend_loc%taux_srf   (ibMax,jbMax))       ;ptend_loc%taux_srf  =0.0_r8!(pcols)! net zonal stress at surface (Pa)
    ALLOCATE(ptend_loc%taux_top   (ibMax,jbMax))       ;ptend_loc%taux_top  =0.0_r8!(pcols)! net zonal stress at top of model (Pa)
    ALLOCATE(ptend_loc%tauy_srf   (ibMax,jbMax))       ;ptend_loc%tauy_srf  =0.0_r8!(pcols)! net meridional stress at surface (Pa)
    ALLOCATE(ptend_loc%tauy_top   (ibMax,jbMax))       ;ptend_loc%tauy_top  =0.0_r8!(pcols)! net meridional stress at top of model (Pa)
    ALLOCATE(ptend_loc%cflx_srf   (ibMax,jbMax,ppcnst));ptend_loc%cflx_srf  =0.0_r8!(pcols,ppcnst)! constituent flux at surface (kg/m2/s)
    ALLOCATE(ptend_loc%cflx_top   (ibMax,jbMax,ppcnst));ptend_loc%cflx_top  =0.0_r8!(pcols,ppcnst)! constituent flux top of model (kg/m2/s)
    !-------------------------------------------------------------------------------
    ! This is for tendencies returned from individual parameterizations
    ALLOCATE(ptend_all%name      (jbMax))            ! true if dsdt is returned
    ALLOCATE(ptend_all%ls         (jbMax))            ! true if dsdt is returned
    ALLOCATE(ptend_all%lu         (jbMax))            ! true if dudt is returned
    ALLOCATE(ptend_all%lv         (jbMax))            ! true if dvdt is returned
    ALLOCATE(ptend_all%top_level  (jbMax))       ! top level index for which nonzero tendencies have been set
    ALLOCATE(ptend_all%bot_level  (jbMax))       ! bottom level index for which nonzero tendencies have been set
    ALLOCATE(ptend_all%lq  (jbMax,ppcnst))                 ;ptend_all%lq =.TRUE.! (ppcnst)          ! true if dqdt() is returned
    ALLOCATE(ptend_all%s   (ibMax,kMax,jbMax))       ;ptend_all%s  =0.0_r8!(pcols,pver)! heating rate (J/kg/s)
    ALLOCATE(ptend_all%u   (ibMax,kMax,jbMax))       ;ptend_all%u  =0.0_r8!(pcols,pver)! u momentum tendency (m/s/s)
    ALLOCATE(ptend_all%v   (ibMax,kMax,jbMax))       ;ptend_all%v  =0.0_r8!(pcols,pver)! v momentum tendency (m/s/s)
    ALLOCATE(ptend_all%q   (ibMax,kMax,jbMax,ppcnst));ptend_all%q  =0.0_r8!(pcols,pver,ppcnst)! consituent tendencies (kg/kg/s)
    ALLOCATE(ptend_all%qt   (ibMax,kMax,jbMax));ptend_all%qt  =0.0_r8!(pcols,pver,ppcnst)! consituent tendencies (kg/kg/s)

    ! boundary fluxes

    ALLOCATE(ptend_all%hflux_srf  (ibMax,jbMax))       ;ptend_all%hflux_srf=0.0_r8 !(pcols)! net heat flux at surface (W/m2)
    ALLOCATE(ptend_all%hflux_top  (ibMax,jbMax))       ;ptend_all%hflux_top =0.0_r8!(pcols)! net heat flux at top of model (W/m2)
    ALLOCATE(ptend_all%taux_srf   (ibMax,jbMax))       ;ptend_all%taux_srf =0.0_r8 !(pcols)! net zonal stress at surface (Pa)
    ALLOCATE(ptend_all%taux_top   (ibMax,jbMax))       ;ptend_all%taux_top =0.0_r8 !(pcols)! net zonal stress at top of model (Pa)
    ALLOCATE(ptend_all%tauy_srf   (ibMax,jbMax))       ;ptend_all%tauy_srf =0.0_r8 !(pcols)! net meridional stress at surface (Pa)
    ALLOCATE(ptend_all%tauy_top   (ibMax,jbMax))       ;ptend_all%tauy_top  =0.0_r8!(pcols)! net meridional stress at top of model (Pa)
    ALLOCATE(ptend_all%cflx_srf   (ibMax,jbMax,ppcnst));ptend_all%cflx_srf  =0.0_r8!(pcols,ppcnst)! constituent flux at surface (kg/m2/s)
    ALLOCATE(ptend_all%cflx_top   (ibMax,jbMax,ppcnst));ptend_all%cflx_top  =0.0_r8!(pcols,ppcnst)! constituent flux top of model (kg/m2/s)

    !ALLOCATE(state_qcwat (ibMax,kMax,jbMax));state_qcwat =0.0_r8!(1:pcols,1:pver)  ! cloud water old q
    !ALLOCATE(state_tcwat (ibMax,kMax,jbMax));state_tcwat =0.0_r8!(1:pcols,1:pver)  ! cloud water old temperature
    !ALLOCATE(state_lcwat (ibMax,kMax,jbMax));state_lcwat =0.0_r8!(1:pcols,1:pver)  ! cloud liquid water old q
    !ALLOCATE(state_cld   (ibMax,kMax,jbMax));state_cld   =0.0_r8!(1:pcols,1:pver)  ! cloud fraction
    !ALLOCATE(state_qme   (ibMax,kMax,jbMax));state_qme   =0.0_r8!(1:pcols,1:pver) 
    !ALLOCATE(state_prain (ibMax,kMax,jbMax));state_prain =0.0_r8!(1:pcols,1:pver) 
    !ALLOCATE(state_nevapr(ibMax,kMax,jbMax));state_nevapr=0.0_r8!(1:pcols,1:pver) 
    !ALLOCATE(state_rel   (ibMax,kMax,jbMax));state_rel   =0.0_r8!(1:pcols,1:pver) ! liquid effective drop radius (microns)
    !ALLOCATE(state_rei   (ibMax,kMax,jbMax));state_rei   =0.0_r8!(1:pcols,1:pver) ! ice effective drop size (microns)

    CALL stratiform_init( kMax,jMax,hypi)
    CALL cldinti (kMax,hypi)
  END SUBROUTINE Init_Micro_Hack





  !===============================================================================
  SUBROUTINE stratiform_init( pver,plat,hypi)
    !----------------------------------------------------------------------
    !
    ! Initialize the cloud water parameterization
    ! 
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: pver                  ! number of vertical levels
    INTEGER, INTENT(IN   ) :: plat   
    REAL(KIND=r8), INTENT(IN) :: hypi(pver+1)

    !---------------------------Local workspace-----------------------------

    !LOGICAL ip           ! Ice phase (true or false)

    ! initialization routine for prognostic cloud water
    CALL inimc (pver,plat ,tmelt, rhodair/1000.0_r8, rh2o, hypi)
    !
    !-----------------------------------------------------------------------
    !
    ! Specify control parameters first
    !
    !ip    = .TRUE.
    CALL cldfrc_init(pver,plat)
    CALL gestbl(epsilo  , &
         latvap  ,latice  ,rh2o    ,cpair    )

    RETURN
  END SUBROUTINE stratiform_init


subroutine cldinti (pver,hypi)
   INTEGER, INTENT(IN   ) :: pver
   REAL(KIND=r8), INTENT(IN) :: hypi(pver+1)
   INTEGER ::k
!
! Find vertical level nearest 700 mb 0r 0.7 sigma
!
   k700 = 1
   do k=1,pver-1
      if (hypi(k) < 0.7_r8 .and. hypi(k+1) >= 0.7_r8) then
         if (0.7_r8-hypi(k) < hypi(k+1)-0.7_r8) then
            k700 = k
         else
            k700 = k + 1
         end if
         goto 20
      end if
   end do

!   call endrun ('CLDINTI: model levels bracketing 700 mb not found')

20 continue
    PRINT*,' k700 =' ,k700
!   if (masterproc) then
!      write(6,*)'CLDINTI: model level nearest 700 mb is',k700,'which is',hypm(k700),'pascals'
!   end if

   return
end subroutine cldinti


  SUBROUTINE RunMicro_Hack(&
       nstep       , &
       pcols       , &! INTEGER , INTENT(in) :: pcols                 ! number of columns (max)
       pver        , &! INTEGER , INTENT(in) :: pver                  ! number of vertical levels
       pverp       , &! INTEGER , INTENT(in) :: pverp                 ! number of vertical levels + 1
       ppcnst      , &! INTEGER , INTENT(in) :: ppcnst          ! number of constituent
       latco       , &! INTEGER , INTENT(in) :: ppcnst          ! number of constituent
       dtime       , &! REAL(r8), INTENT(in)  :: dtime                   ! timestep
       prsi        , &
       prsl        , &
       state_t2    , &! REAL(r8), INTENT(in)  :: state_t      (pcols,pver)!(pcols,pver)! temperature (K)
       state_t     , &! REAL(r8), INTENT(in)  :: state_t      (pcols,pver)!(pcols,pver)! temperature (K)
       state_qv2   , &! REAL(r8), INTENT(in)  :: state_qv     (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       state_qv    , &! REAL(r8), INTENT(in)  :: state_qv     (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       state_ql2   , &! REAL(r8), INTENT(in)  :: state_ql     (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       state_ql    , &! REAL(r8), INTENT(in)  :: state_ql     (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       state_qi2    , &! REAL(r8), INTENT(in)  :: state_qi     (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       state_qi    , &! REAL(r8), INTENT(in)  :: state_qi     (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       state_omega , &! REAL(r8), INTENT(in)  :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
       icefrac     , &! REAL(r8), INTENT(in)  :: icefrac      (pcols)             ! sea ice fraction (fraction)
       landfrac    , &! REAL(r8), INTENT(in)  :: landfrac     (pcols)             ! land fraction (fraction)
       ocnfrac     , &! REAL(r8), INTENT(in)  :: ocnfrac      (pcols)             ! ocean fraction (fraction)
       landm       , &! REAL(r8), INTENT(in)  :: landm              (pcols)             ! land fraction ramped over water
       snowh       , &! REAL(r8), INTENT(in)  :: snowh              (pcols)             ! Snow depth over land, water equivalent (m)
       state_dlf   , &! REAL(r8), INTENT(in)  :: state_dlf    (pcols,pver)    ! detrained water from ZM
       rliq        , &! REAL(r8), INTENT(in)  :: rliq         (pcols)         ! vertical integral of liquid not yet in q(ixcldliq)
       state_cmfmc , &! REAL(r8), INTENT(in)  :: state_cmfmc  (pcols,pverp)   ! convective mass flux--m sub c
       state_cmfmc2, &! REAL(r8), INTENT(in)  :: state_cmfmc2 (pcols,pverp)   ! shallow convective mass flux--m sub c
       state_concld, &! REAL(r8), INTENT(out) :: state_concld (pcols,pver)    ! convective cloud cover
       state_cld   , &! REAL(r8), INTENT(out) :: state_concld (pcols,pver)    !  cloud fraction
       sst         , &! REAL(r8), INTENT(in)  :: sst              (pcols)              ! sea surface temperature
       !state_zdu   , &! REAL(r8), INTENT(in)  :: state_zdu    (pcols,pver)           ! detrainment rate from deep convection
       prec_str    , &! REAL(r8), INTENT(out)  :: prec_str    (pcols)         ! [Total] sfc flux of precip from stratiform (m/s) 
       snow_str    , &! REAL(r8), INTENT(out)  :: snow_str    (pcols)  ! [Total] sfc flux of snow from stratiform   (m/s)
       prec_sed    , &! REAL(r8), INTENT(out)  :: prec_sed    (pcols)  ! surface flux of total cloud water from sedimentation
       snow_sed    , &! REAL(r8), INTENT(out)  :: snow_sed    (pcols)  ! surface flux of cloud ice from sedimentation
       prec_pcw    , &! REAL(r8), INTENT(out)  :: prec_pcw    (pcols)  ! sfc flux of precip from microphysics(m/s)
       snow_pcw    , &! REAL(r8), INTENT(out)  :: snow_pcw    (pcols)  ! sfc flux of snow from microphysics (m/s)
       dtdt     , &
       dqdt     , &
       dqldt    , &
       dqidt    )
    !       pbuf )
    !
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: nstep
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                 ! number of vertical levels + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: ppcnst
    INTEGER, INTENT(in) :: latco

    REAL(r8), INTENT(in)  :: dtime                ! timestep
    REAL(KIND=r8)   , INTENT(in   ) :: prsi    (pcols,pver+1)  
    REAL(KIND=r8)   , INTENT(in   ) :: prsl    (pcols,pver  ) 
    REAL(r8), INTENT(in   )  :: state_t2  (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_t   (pcols,pver)
    REAL(r8), INTENT(in   )  :: state_qv2 (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_qv  (pcols,pver)
    REAL(r8), INTENT(in   )  :: state_ql2 (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_ql  (pcols,pver)
    REAL(r8), INTENT(in   )  :: state_qi2 (pcols,pver)
    REAL(r8), INTENT(inout)  :: state_qi  (pcols,pver)

    REAL(r8), INTENT(in)  :: state_omega  (pcols,pver)

    REAL(r8), INTENT(in)  :: icefrac      (pcols)         ! sea ice fraction (fraction)
    REAL(r8), INTENT(in)  :: landfrac     (pcols)         ! land fraction (fraction)
    REAL(r8), INTENT(in)  :: ocnfrac      (pcols)         ! ocean fraction (fraction)
    REAL(r8), INTENT(in)  :: landm        (pcols)         ! land fraction ramped over water
    REAL(r8), INTENT(in)  :: snowh        (pcols)         ! Snow depth over land, water equivalent (m)

    REAL(r8), INTENT(in) :: rliq          (pcols)         ! vertical integral of liquid not yet in q(ixcldliq)
    REAL(r8), INTENT(in) :: state_dlf     (pcols,pver)    ! detrained water from ZM
    REAL(r8), INTENT(in) :: state_cmfmc   (pcols,pverp)   ! convective mass flux--m sub c
    REAL(r8), INTENT(in) :: state_cmfmc2  (pcols,pverp)   ! shallow convective mass flux--m sub c

    REAL(r8), INTENT(in) :: sst           (pcols)               ! sea surface temperature
    !REAL(r8), INTENT(in) :: state_zdu     (pcols,pver)          ! detrainment rate from deep convection

    REAL(r8), INTENT(out)  :: prec_str(pcols)  ! [Total] sfc flux of precip from stratiform (m/s) 
    REAL(r8), INTENT(out)  :: snow_str(pcols)  ! [Total] sfc flux of snow from stratiform   (m/s)
    REAL(r8), INTENT(out)  :: prec_sed(pcols)  ! surface flux of total cloud water from sedimentation
    REAL(r8), INTENT(out)  :: snow_sed(pcols)  ! surface flux of cloud ice from sedimentation
    REAL(r8), INTENT(out)  :: prec_pcw(pcols)  ! sfc flux of precip from microphysics(m/s)
    REAL(r8), INTENT(out)  :: snow_pcw(pcols)  ! sfc flux of snow from microphysics (m/s)
    REAL(r8), INTENT(out)  :: state_concld(pcols,pver)      ! convective cloud cover
    REAL(r8), INTENT(out)  :: state_cld  (pcols,pver)
    REAL(KINd=r8), INTENT(OUT  ) :: dtdt (pcols,pver)
    REAL(KIND=r8), INTENT(OUT  ) :: dqdt (pcols,pver)
    REAL(KIND=r8), INTENT(OUT  ) :: dqldt(pcols,pver)
    REAL(KIND=r8), INTENT(OUT  ) :: dqidt(pcols,pver)

    ! physics buffer fields
    REAL(r8) :: concld (pcols,pver)     ! convective cloud cover
    REAL(r8) :: cld   (1:pcols,1:pver)  ! cloud fraction
    REAL(r8) :: qcwat (1:pcols,1:pver)  ! cloud water old q
    REAL(r8) :: tcwat (1:pcols,1:pver)  ! cloud water old temperature
    REAL(r8) :: lcwat (1:pcols,1:pver)  ! cloud liquid water old q
    REAL(r8) :: qme   (1:pcols,1:pver) 
    REAL(r8) :: prain (1:pcols,1:pver) 
    REAL(r8) :: nevapr(1:pcols,1:pver) 
    REAL(r8) :: rel   (1:pcols,1:pver)   ! liquid effective drop radius (microns)
    REAL(r8) :: rei   (1:pcols,1:pver)   ! ice effective drop size (microns)

    REAL(r8) :: dlf   (pcols,pver)        ! detrained water from ZM
    REAL(r8) :: cmfmc (pcols,pverp)        ! convective mass flux--m sub c
    REAL(r8) :: cmfmc2(pcols,pverp)        ! shallow convective mass flux--m sub c
    !REAL(r8) :: zdu   (pcols,pver)        ! detrainment rate from deep convection

    ! Parameters
    !
    REAL(r8) :: pnot  (pcols)                ! reference pressure
    !PARAMETER (pnot = 1.e5_r8)

    !
    ! Local variables
    !



    INTEGER i,k


    ! local variables for stratiform_sediment
    REAL(r8) :: rain(pcols)                       ! surface flux of cloud liquid
    REAL(r8) :: pvliq(pcols,pver+1)               ! vertical velocity of cloud liquid drops (Pa/s)
    REAL(r8) :: pvice(pcols,pver+1)               ! vertical velocity of cloud ice particles (Pa/s)

    ! local variables for cldfrc
    REAL(r8) :: cldst(pcols,pver)     ! cloud fraction
    REAL(r8) :: clc(pcols)            ! column convective cloud amount
    REAL(r8) :: rhdfda(pcols,pver)    ! d_RH/d_cloud_fraction    ====wlin
    REAL(r8) :: rhu00(pcols,pver)     ! RH limit, U00             ====wlin
    REAL(r8) :: relhum(pcols,pver)         ! RH, output to determine drh/da
    REAL(r8) :: rhu002(pcols,pver)         ! same as rhu00 but for perturbed rh 
    REAL(r8) :: cld2(pcols,pver)          ! same as cld but for perturbed rh  
    REAL(r8) :: concld2(pcols,pver)        ! same as concld but for perturbed rh 
    REAL(r8) :: cldst2(pcols,pver)         ! same as cldst but for perturbed rh 
    REAL(r8) :: relhum2(pcols,pver)        ! RH after  perturbation            
    REAL(r8) :: pmid(pcols,pver)      ! midpoint pressures
    REAL(r8) :: t(pcols,pver)         ! temperature
    REAL(r8) :: q(pcols,pver)         ! specific humidity
    REAL(r8) :: omga(pcols,pver)      ! vertical pressure velocity
    REAL(r8) :: phis(pcols)           ! surface geopotential
    REAL(r8) :: pdel(pcols,pver)      ! pressure depth of layer
    REAL(r8) :: ps(pcols)             ! surface pressure

    ! local variables for microphysics
    REAL(r8) :: rdtime                          ! 1./dtime
    REAL(r8) :: qtend (pcols,pver)            ! moisture tendencies
    REAL(r8) :: ttend (pcols,pver)            ! temperature tendencies
    REAL(r8) :: ltend (pcols,pver)            ! cloud liquid water tendencies
    REAL(r8) :: evapheat(pcols,pver)          ! heating rate due to evaporation of precip
    REAL(r8) :: evapsnow(pcols,pver)          ! local evaporation of snow
    REAL(r8) :: prfzheat(pcols,pver)          ! heating rate due to freezing of precip (W/kg)
    REAL(r8) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip
    REAL(r8) :: cmeheat (pcols,pver)          ! heating rate due to phase change of precip
    REAL(r8) :: prodsnow(pcols,pver)          ! local production of snow
    REAL(r8) :: totcw   (pcols,pver)          ! total cloud water mixing ratio
    REAL(r8) :: fice    (pcols,pver)          ! Fractional ice content within cloud
    REAL(r8) :: fsnow   (pcols,pver)          ! Fractional snow production
    REAL(r8) :: repartht(pcols,pver)          ! heating rate due to phase repartition of input precip
    REAL(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
    REAL(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
    REAL(r8) :: fwaut(pcols,pver)              
    REAL(r8) :: fsaut(pcols,pver)              
    REAL(r8) :: fracw(pcols,pver)              
    REAL(r8) :: fsacw(pcols,pver)              
    REAL(r8) :: fsaci(pcols,pver)              
    REAL(r8) :: cmeice(pcols,pver)   ! Rate of cond-evap of ice within the cloud
    REAL(r8) :: cmeliq(pcols,pver)   ! Rate of cond-evap of liq within the cloud
    REAL(r8) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
    REAL(r8) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
    REAL(r8) :: liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
    REAL(r8), TARGET :: state_buf_t       (pcols,pver)  
    REAL(r8), TARGET :: state_buf_omega   (pcols,pver) 
    REAL(r8), TARGET :: state_buf_pmid    (pcols,pver) 
    REAL(r8), TARGET :: state_buf_pdel    (pcols,pver) 
    REAL(r8), TARGET :: state_buf_rpdel   (pcols,pver)
    REAL(r8), TARGET :: state_buf_lnpmid  (pcols,pver)
    REAL(r8), TARGET :: state_buf_pint    (pcols,pver+1)
    REAL(r8), TARGET :: state_buf_lnpint  (pcols,pver+1) 
    REAL(r8), TARGET :: hkk(pcols)
    REAL(r8), TARGET :: hkl(pcols)
    REAL(r8), TARGET :: tvfac
    LOGICAL :: fvdyn

    DO i=1,pcols
       prec_str(i)=0.0_r8  ! [Total] sfc flux of precip from stratiform (m/s)
       snow_str(i)=0.0_r8  ! [Total] sfc flux of snow from stratiform   (m/s)
       prec_sed(i)=0.0_r8  ! surface flux of total cloud water from sedimentation
       snow_sed(i)=0.0_r8  ! surface flux of cloud ice from sedimentation
       prec_pcw(i)=0.0_r8  ! sfc flux of precip from microphysics(m/s)
       snow_pcw(i)=0.0_r8  ! sfc flux of snow from microphysics (m/s)
       pnot(i)=0.0_r8
       rain(i)=0.0_r8
       clc(i)=0.0_r8            ! column convective cloud amount
       phis(i)=0.0_r8           ! surface geopotential
       ps(i)=0.0_r8             ! surface pressure
    END DO
    DO k=pver+1,1,-1
       DO i=1,pcols
          pvliq(i,k)=0.0_r8               ! vertical velocity of cloud liquid drops (Pa/s)
          pvice(i,k)=0.0_r8               ! vertical velocity of cloud ice particles (Pa/s)
       END DO
    END DO
    DO k=pver,1,-1
       DO i=1,pcols
          dtdt (i,k)=0.0_r8 
          dqdt (i,k)=0.0_r8 
          dqldt(i,k)=0.0_r8 
          dqidt(i,k)=0.0_r8 
          
          state_concld(i,k)=0.0_r8      ! convective cloud cover
          state_cld  (i,k)=0.0_r8
          concld (i,k)=0.0_r8     ! convective cloud cover
          cld   (i,k)=0.0_r8  ! cloud fraction
          qcwat (i,k)=0.0_r8  ! cloud water old q
          tcwat (i,k)=0.0_r8  ! cloud water old temperature
          lcwat (i,k)=0.0_r8  ! cloud liquid water old q
          qme   (i,k)=0.0_r8
          prain (i,k)=0.0_r8
          nevapr(i,k)=0.0_r8
          rel   (i,k)=0.0_r8   ! liquid effective drop radius (microns)
          rei   (i,k)=0.0_r8   ! ice effective drop size (microns)
          dlf   (i,k)=0.0_r8        ! detrained water from ZM
          cmfmc (i,k)=0.0_r8        ! convective mass flux--m sub c
          cmfmc2(i,k)=0.0_r8        ! shallow convective mass flux--m sub c
    ! local variables for cldfrc
          cldst(i,k)=0.0_r8     ! cloud fraction
          rhdfda(i,k)=0.0_r8    ! d_RH/d_cloud_fraction    ====wlin
          rhu00(i,k) =0.0_r8    ! RH limit, U00             ====wlin
          relhum(i,k)=0.0_r8         ! RH, output to determine drh/da
          rhu002(i,k)=0.0_r8         ! same as rhu00 but for perturbed rh
          cld2(i,k)  =0.0_r8        ! same as cld but for perturbed rh
          concld2(i,k)=0.0_r8        ! same as concld but for perturbed rh
          cldst2(i,k) =0.0_r8        ! same as cldst but for perturbed rh
          relhum2(i,k)=0.0_r8        ! RH after  perturbation
          pmid(i,k)   =0.0_r8   ! midpoint pressures
          t(i,k)      =0.0_r8   ! temperature
          q(i,k)      =0.0_r8   ! specific humidity
          omga(i,k)   =0.0_r8   ! vertical pressure velocity
          pdel(i,k)   =0.0_r8   ! pressure depth of layer
          qtend (i,k) =0.0_r8            ! moisture tendencies
          ttend (i,k) =0.0_r8            ! temperature tendencies
          ltend (i,k) =0.0_r8            ! cloud liquid water tendencies
          evapheat(i,k)=0.0_r8           ! heating rate due to evaporation of precip
          evapsnow(i,k)=0.0_r8           ! local evaporation of snow
          prfzheat(i,k)=0.0_r8           ! heating rate due to freezing of precip (W/kg)
          meltheat(i,k)=0.0_r8           ! heating rate due to phase change of precip
          cmeheat (i,k)=0.0_r8           ! heating rate due to phase change of precip
          prodsnow(i,k)=0.0_r8           ! local production of snow
          totcw   (i,k)=0.0_r8           ! total cloud water mixing ratio
          fice    (i,k)=0.0_r8           ! Fractional ice content within cloud
          fsnow   (i,k)=0.0_r8           ! Fractional snow production
          repartht(i,k)=0.0_r8           ! heating rate due to phase repartition of input precip
          icimr(i,k)=0.0_r8              ! in cloud ice mixing ratio
          icwmr(i,k)=0.0_r8              ! in cloud water mixing ratio
          fwaut(i,k)=0.0_r8 
          fsaut(i,k)=0.0_r8 
          fracw(i,k)=0.0_r8 
          fsacw(i,k)=0.0_r8 
          fsaci(i,k)=0.0_r8 
          cmeice(i,k)=0.0_r8    ! Rate of cond-evap of ice within the cloud
          cmeliq(i,k)=0.0_r8    ! Rate of cond-evap of liq within the cloud
          ice2pr(i,k)=0.0_r8    ! rate of conversion of ice to precip
          liq2pr(i,k)=0.0_r8    ! rate of conversion of liquid to precip
          liq2snow(i,k)=0.0_r8    ! rate of conversion of liquid to snow
       END DO
    END DO
    !======================================================================
    !======================================================================
    !======================================================================
    DO i=1,pcols
       !state_buf_pint       (i,pver+1) = state_ps(i)*si(1)
       state_buf_pint       (i,pver+1) = prsi(i,1)
    END DO
    DO k=pver,1,-1
       DO i=1,pcols
          !state_buf_pint    (i,k)      = MAX(si(pver+2-k)*state_ps(i) ,0.0001_r8)
          state_buf_pint    (i,k)      = MAX( prsi(i,pver+2-k) ,0.0001_r8)
       END DO
    END DO

    DO k=1,pver+1
       DO i=1,pcols
          state_buf_lnpint(i,k) =  LOG(state_buf_pint  (i,k))
          cmfmc (i,pver+2-k)          =  state_cmfmc (i,k)
          cmfmc2(i,pver+2-k)          =  state_cmfmc2(i,k)
       END DO
    END DO
    DO k=1,pver+1
       DO i=1,pcols
          state%pint     (i,k,latco) = state_buf_pint  (i,k)
          state%lnpint   (i,k,latco) = state_buf_lnpint(i,k)
       END DO
    END DO
    DO k=1,pver
       DO i=1,pcols
          state_buf_t       (i,pver+1-k) = state_t          (i,k)
          state_buf_omega   (i,pver+1-k) = state_omega    (i,k)
          !state_buf_pmid    (i,pver+1-k) = sl(k)*state_ps (i)
          state_buf_pmid    (i,pver+1-k) = prsl(i,k) !sl(k)*state_ps (i)
          dlf               (i,pver+1-k) = state_dlf (i,k) 
          !zdu               (i,pver+1-k) = state_zdu (i,k) 
       END DO
    END DO
    DO k=1,pver
       DO i=1,pcols          
          state_buf_pdel    (i,k) = MAX(state%pint(i,k+1,latco) - state%pint(i,k,latco),0.000000005_r8)
          state_buf_rpdel   (i,k) = 1.0_r8/MAX((state%pint(i,k+1,latco) - state%pint(i,k,latco)),0.00000000005_r8)
          state_buf_lnpmid  (i,k) = LOG(state_buf_pmid(i,k))        
       END DO
    END DO
    DO k=1,pver
       DO i=1,pcols          
          state%t        (i,k,latco)= state_buf_t     (i,k) 
          state%omega    (i,k,latco)= state_buf_omega (i,k) 
          state%pmid     (i,k,latco)= state_buf_pmid  (i,k) 
          state%pdel     (i,k,latco)= state_buf_pdel  (i,k) 
          state%rpdel    (i,k,latco)= state_buf_rpdel (i,k) 
          state%lnpmid   (i,k,latco)= state_buf_lnpmid(i,k) 
       END DO
    END DO

    DO k=1,pver
       DO i=1,pcols
          state%q    (i,pver+1-k,latco,1)         =  state_qv(i,k)
          state%q    (i,pver+1-k,latco,ixcldliq)  =  state_ql(i,k)
          state%q    (i,pver+1-k,latco,ixcldice)  =  state_qi(i,k)
       END DO
    END DO
    ! Derive new temperature and geopotential fields

    CALL geopotential_t(                                 &
       state%lnpint(1:pcols,1:pver+1,latco)   , state%pint (1:pcols,1:pver+1,latco)   , &
       state%pmid  (1:pcols,1:pver,latco)     , state%pdel  (1:pcols,1:pver,latco)   , state%rpdel(1:pcols,1:pver,latco)   , &
       state%t     (1:pcols,1:pver,latco)     , state%q     (1:pcols,1:pver,latco,1) , rair   , gravit , zvir   ,          &
       state%zi    (1:pcols,1:pver+1,latco)   , state%zm    (1:pcols,1:pver,latco)   , pcols, pver, pverp)


    fvdyn = dycore_is ('LR')
    DO k = pver, 1, -1
       ! First set hydrostatic elements consistent with dynamics
       IF (fvdyn) THEN
          DO i = 1,pcols
             hkl(i) = state%lnpmid(i,k+1,latco) - state%lnpmid(i,k,latco)
             hkk(i) = 1.0_r8 - state%pint (i,k,latco)* hkl(i)* state%rpdel(i,k,latco)
          END DO
       ELSE
          DO i = 1,pcols
             hkl(i) = state%pdel(i,k,latco)/state%pmid  (i,k,latco)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF
       ! Now compute s
       DO i = 1,pcols
          tvfac   = 1.0_r8 + zvir * state%q(i,k,latco,1) 
          state%s(i,k,latco) =  (state%t(i,k,latco)* cpair) + (state%t(i,k,latco) * tvfac * rair*hkk(i))  +  &
                                ( state%phis(i,latco) + gravit*state%zi(i,k+1,latco))
       END DO
    END DO

    IF(nstep == 1)THEN
    DO k=1,pver
       DO i=1,pcols          
          qcwat(i,pver+1-k) =  state_qv2 (i,k)! cloud water old q
          tcwat(i,pver+1-k) =  state_t2  (i,k)  ! cloud water old temperature
          lcwat(i,pver+1-k) =  state_qi2 (i,k) + state_ql2(i,k)! cloud liquid water old q
       END DO
    END DO    
    ELSE 

    DO k=1,pver
       DO i=1,pcols
          !concld(i,pver+1-k)= state_concld(i,k)          ! cloud water old q
          qcwat(i,k) = state_qv2(i,pver+1-k)              ! cloud water old q
          !PK qcwat (i,pver+1-k)= state_qcwat (i,k,latco) ! cloud water old q
          tcwat(i,k) = state_t2 (i,pver+1-k)              ! cloud water old temperature
          !PK tcwat (i,pver+1-k)= state_tcwat (i,k,latco)! cloud water old temperature
          lcwat(i,k) = state_qi2 (i,pver+1-k) + state_ql2(i,pver+1-k)! cloud liquid water old q
          !PK lcwat (i,pver+1-k)= state_lcwat (i,k,latco)! cloud liquid water old q
          !cld   (i,pver+1-k)= state_cld   (i,k,latco)! cloud fraction
          !qme   (i,pver+1-k)= state_qme   (i,k,latco)
          !prain (i,pver+1-k)= state_prain (i,k,latco)
          !nevapr(i,pver+1-k)= state_nevapr(i,k,latco)
          !rel   (i,pver+1-k)= state_rel   (i,k,latco) ! liquid effective drop radius (microns)
          !rei   (i,pver+1-k)= state_rei   (i,k,latco) ! ice effective drop size (microns)
       END DO
    END DO
    END IF


    !-------------------------------------------------------------------------------
    ! This is for tendencies returned from individual parameterizations

    CALL physics_state_copy(pcols,latco,pver, pverp,ppcnst)   ! copy state to local state1.
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)  ! initialize local ptend type
    CALL physics_ptend_init(ptend_all,latco,ppcnst,pver)  ! initialize output ptend type
    CALL physics_tend_init (tend,latco)                    ! tend here is just a null place holder
    !+++sediment ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Allow the cloud liquid drops and ice particles to sediment
    ! Occurs before adding convectively detrained cloud water, because the phase of the
    ! of the detrained water is unknown.
    !call t_startf('stratiform_sediment')


    ptend_loc%name(latco)        = 'pcwsediment'
    ptend_loc%ls(latco)          = .TRUE.
    ptend_loc%lq(latco,1)        = .TRUE.
    ptend_loc%lq(latco,ixcldice) = .TRUE.
    ptend_loc%lq(latco,ixcldliq) = .TRUE.
    !
    ! cloud fraction after transport and convection,
    ! derive the relationship between rh and cld from 
    ! the employed c7loud scheme
    !
    DO k=1,pver
       DO i=1,pcols
          pmid(i,k) = state1%pmid  (i,k,latco)
          t   (i,k) = state1%t     (i,k,latco)
          q   (i,k) = state1%q     (i,k,latco,1)
          omga(i,k) = state1%omega (i,k,latco)
          pdel(i,k) = state1%pdel  (i,k,latco)
       END DO
    END DO
    DO i=1,pcols
       ps  (i)        = state1%pint  (i,pverp,latco)
       phis(i)        = state1%phis  (i,latco)
    END DO

    CALL cldfrc(&
         pcols                    , &! !INTENT(in ) :: pcols           ! number of columns  (max)
         pver                     , &! !INTENT(in ) :: pver        ! number of vertical levels
         pverp                    , &! !INTENT(in ) :: pverp           ! pver + 1
         pcols                    , &! !INTENT(in ) :: pcols                    ! number of atmospheric columns
         pnot    (1:pcols)        , &! !INTENT(in ) :: pnot    (pcols)                ! reference pressure
         pmid    (1:pcols,1:pver) , &! !INTENT(in ) :: pmid    (pcols,pver)      ! midpoint pressures
         t       (1:pcols,1:pver) , &! !INTENT(in ) :: temp    (pcols,pver)      ! temperature
         q       (1:pcols,1:pver) , &! !INTENT(in ) :: q       (pcols,pver)      ! specific humidity
         phis    (1:pcols)        , &! !INTENT(in ) :: phis    (pcols)               ! surface geopotential 
         cld     (1:pcols,1:pver) , &! !INTENT(out) :: cloud   (pcols,pver)      ! cloud fraction
         clc     (1:pcols)        , &! !INTENT(out) :: clc     (pcols)               ! column convective cloud amount 
         cmfmc   (1:pcols,1:pverp), &! !INTENT(in ) :: cmfmc   (pcols,pverp)     ! convective mass flux--m sub c
         cmfmc2  (1:pcols,1:pverp), &! !INTENT(in ) :: cmfmc2  (pcols,pverp)     ! shallow convective mass flux--m sub c
         landfrac(1:pcols)        , &! !INTENT(in ) :: landfrac(pcols)               ! Land fraction 
         snowh   (1:pcols)        , &! !INTENT(in ) :: snowh   (pcols)               ! snow depth (liquid water equivalent) ! sea surface temperature
         concld  (1:pcols,1:pver) , &! !INTENT(out) :: concld  (pcols,pver)   ! convective cloud cover
         cldst   (1:pcols,1:pver) , &! !INTENT(in ) :: cldst   (pcols,pver)      ! not used ! detrainment rate from deep convection
         sst     (1:pcols)        , &! !INTENT(in ) :: sst     (pcols)             ! sea surface temperature
         ps      (1:pcols)        , &! !INTENT(in ) :: ps      (pcols)             ! surface pressure
         !zdu     (1:pcols,1:pver), &! !INTENT(in ) :: zdu     (pcols,pver)      ! not used ! detrainment rate from deep convection
         ocnfrac (1:pcols)        , &! !INTENT(in ) :: ocnfrac (pcols)               ! Ocean fraction
         rhu00   (1:pcols,1:pver) , &! !INTENT(out) :: rhu00   (pcols,pver)      ! RH threshold for cloud
         relhum  (1:pcols,1:pver) , &! !INTENT(out) :: relhum  (pcols,pver)      ! RH 
         0                          ) ! !INTENT(in ) :: dindex                     ! 0 or 1 to perturb rh

    CALL cld_sediment_vel ( &
         pcols                                  , &!INTENT(in) :: pcols                   ! number of colums to process
         pcols                                  , &!INTENT(in) :: pcols                 ! number of columns (max)
         pver                                   , &!INTENT(in) :: pver                  ! number of vertical levels
         pverp                                  , &!INTENT(in) :: pverp                  ! number of vertical levels + 1
         icefrac    (1:pcols)                    , &!INTENT(in) :: icefrac (pcols)          ! sea ice fraction (fraction)
         landfrac   (1:pcols)                    , &!INTENT(in) :: landfrac(pcols)          ! land fraction (fraction)
         state1%pmid(1:pcols,1:pver,latco)      , &!INTENT(in) :: pmid  (pcols,pver)          ! pressure of midpoint levels (Pa)
         state1%t   (1:pcols,1:pver,latco)      , &!INTENT(in) :: t         (pcols,pver)          ! temperature (K)
         cld        (1:pcols,1:pver)            , &!INTENT(in) :: cloud (pcols,pver)          ! cloud fraction (fraction)
         state1%q(1:pcols,1:pver,latco,ixcldliq), &!INTENT(in) :: cldliq(pcols,pver)          ! cloud water, liquid (kg/kg)
         state1%q(1:pcols,1:pver,latco,ixcldice), &!INTENT(in) :: cldice(pcols,pver)          ! cloud water, ice        (kg/kg)
         pvliq   (1:pcols,1:pver+1)             , &!INTENT(out):: pvliq (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)
         pvice   (1:pcols,1:pver+1)             , &!INTENT(out):: pvice (pcols,pverp)    ! vertical velocity of cloud ice particles (Pa/s)
         landm   (1:pcols)                       , &!INTENT(in) :: landm(pcols)            ! land fraction ramped over water
         snowh   (1:pcols)                         )!INTENT(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)


  
    CALL cld_sediment_tend (&
         pcols                                      , &!INTENT(in)  :: ncol                         ! number of colums to process
         pcols                                      , &!INTENT(in)  :: pcols                    ! number of columns (max)
         pver                                       , &!INTENT(in)  :: pver                    ! number of vertical levels
         pverp                                      , &!INTENT(in)  :: pverp                     ! number of vertical levels + 1
         dtime                                      , &!INTENT(in)  :: dtime                         ! time step
         state1%pint(1:pcols,1:pver+1,latco)        , &!INTENT(in)  :: pint  (pcols,pverp)         ! interfaces pressure (Pa)
         state1%pdel(1:pcols,1:pver,latco)          , &!INTENT(in)  :: pdel  (pcols,pver)         ! pressure diff across layer (Pa)
         cld        (1:pcols,1:pver)                , &!INTENT(in)  :: cloud (pcols,pver)         ! cloud fraction (fraction)
         state1%q   (1:pcols,1:pver,latco,ixcldliq) , &!INTENT(in)  :: cldliq(pcols,pver)         ! cloud liquid water (kg/kg)
         state1%q   (1:pcols,1:pver,latco,ixcldice) , &!INTENT(in)  :: cldice(pcols,pver)         ! cloud ice water    (kg/kg)
         pvliq      (1:pcols,1:pver+1)              , &!INTENT(in ) :: pvliq (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)
         pvice      (1:pcols,1:pver+1)              , &!INTENT(in ) :: pvice (pcols,pverp)    ! vertical velocity of cloud ice particles (Pa/s)
         ptend_loc%q(1:pcols,1:pver,latco,ixcldliq) , &!INTENT(out) :: liqtend(pcols,pver)         ! liquid condensate tend
         ptend_loc%q(1:pcols,1:pver,latco,ixcldice) , &!INTENT(out) :: icetend(pcols,pver)         ! ice condensate tend
         ptend_loc%q(1:pcols,1:pver,latco,1)        , &!INTENT(out) :: wvtend (pcols,pver)         ! water vapor tend
         ptend_loc%s(1:pcols,1:pver,latco)          , &!INTENT(out) :: htend  (pcols,pver)         ! heating rate
         rain       (1:pcols)                        , &!INTENT(out) :: sfliq  (pcols)                 ! surface flux of liquid (rain, kg/m/s)
         snow_sed   (1:pcols)                          )!INTENT(out) :: sfice  (pcols)                 ! surface flux of ice    (snow, kg/m/s)

    DO i=1,pcols

       ! convert rain and snow from kg/m2 to m/s
       snow_sed(i) = snow_sed(i)/1000.0_r8
       rain    (i) = rain(i)/1000.0_r8
       ! compute total precip (m/s)
       prec_sed(i) = rain(i) + snow_sed(i)
    END DO


    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, ppcnst,latco,pcols)
    ! update physics state type state1 with ptend_loc 
    CALL physics_update(state1, tend, ptend_loc, dtime,ppcnst,pcols,pver,latco)
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)



    ! accumulate prec and snow
    DO i=1,pcols
       prec_str(i) = prec_sed(i)
       snow_str(i) = snow_sed(i)
    END DO
    !++detrain ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Put the detraining cloud water from convection into the cloud and environment. 

    ptend_loc%name(latco)         = 'pcwdetrain'
!!$    ptend_loc%ls           = .TRUE.
!!$    ptend_loc%lq(1)        = .TRUE.
!!$    ptend_loc%lq(ixcldice) = .TRUE.
    ptend_loc%lq(latco,ixcldliq) = .TRUE.
    !
    ! Put all of the detraining cloud water from convection into the large scale cloud.
    ! It all goes in liquid for the moment.
    DO k = 1,pver
       DO i = 1,pcols
!!$          ptend_loc%q(i,k,1)        = dlf(i,k) * (1.-cld(i,k))
!!$          ptend_loc%s(i,k)          =-dlf(i,k) * (1.-cld(i,k))*latvap
!!$          ptend_loc%q(i,k,ixcldice) = 0.
!!$          ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * cld(i,k)
          ptend_loc%q(i,k,latco,ixcldliq) = dlf(i,k)
       END DO
    END DO


    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, ppcnst,latco,pcols)
    CALL physics_update(state1, tend, ptend_loc, dtime,ppcnst,pcols,pver,latco)
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)


    ! accumulate prec and snow, reserved liquid has now been used
    DO i=1,pcols
       prec_str(i) = prec_str(i) - rliq(i)  ! ( snow contribution is zero )
    ENDDO

    !++++ cldfrc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !
    ! cloud fraction after transport and convection,
    ! derive the relationship between rh and cld from 
    ! the employed c7loud scheme
    !
    DO k=1,pver
       DO i=1,pcols
          pmid(i,k) = state1%pmid  (i,k,latco)
          t   (i,k) = state1%t     (i,k,latco)
          q   (i,k) = state1%q     (i,k,latco,1)
          omga(i,k) = state1%omega (i,k,latco)
          pdel(i,k) = state1%pdel  (i,k,latco)
       ENDDO
    ENDDO
    DO i=1,pcols
       ps  (i)        = state1%pint  (i,pverp,latco)
       phis(i)        = state1%phis  (i,latco)
    ENDDO

    !

    CALL cldfrc(&
         pcols                    , &! !INTENT(in) ::  pcols           ! number of columns  (max)
         pver                     , &! !INTENT(in) ::  pver ! number of vertical levels
         pverp                    , &! !INTENT(in) ::  pverp           ! pver + 1
         pcols                     , &! !INTENT(in) :: pcols                    ! number of atmospheric columns
         pnot    (1:pcols)         , &! !INTENT(in) :: pnot    (pcols)                ! reference pressure
         pmid    (1:pcols,1:pver) , &! !INTENT(in) :: pmid    (pcols,pver)      ! midpoint pressures
         t       (1:pcols,1:pver) , &! !INTENT(in) :: temp    (pcols,pver)      ! temperature
         q       (1:pcols,1:pver) , &! !INTENT(in) :: q       (pcols,pver)      ! specific humidity
         phis    (1:pcols)         , &! !INTENT(in) :: phis    (pcols)               ! surface geopotential 
         cld     (1:pcols,1:pver) , &! !INTENT(out) :: cloud  (pcols,pver)      ! cloud fraction
         clc     (1:pcols)         , &! !INTENT(out) :: clc    (pcols)               ! column convective cloud amount 
         cmfmc   (1:pcols,1:pverp), &! !INTENT(in) :: cmfmc   (pcols,pverp)     ! convective mass flux--m sub c
         cmfmc2  (1:pcols,1:pverp), &! !INTENT(in) :: cmfmc2  (pcols,pverp)     ! shallow convective mass flux--m sub c
         landfrac(1:pcols)         , &! !INTENT(in) :: landfrac(pcols)               ! Land fraction 
         snowh   (1:pcols)         , &! !INTENT(in) :: snowh   (pcols)               ! snow depth (liquid water equivalent) ! sea surface temperature
         concld  (1:pcols,1:pver) , &!  INTENT(out) concld (pcols,pver)   ! convective cloud cover
         cldst   (1:pcols,1:pver) , &! !INTENT(in) :: cldst     (pcols,pver)      ! not used ! detrainment rate from deep convection
         sst     (1:pcols)         , &! !INTENT(in) :: sst   (pcols)             ! sea surface temperature
         ps      (1:pcols)         , &! !INTENT(in) :: ps    (pcols)             ! surface pressure
         !zdu     (1:pcols,1:pver) , &! !INTENT(in) :: zdu   (pcols,pver)      ! not used ! detrainment rate from deep convection
         ocnfrac (1:pcols)         , &! !INTENT(in) :: ocnfrac (pcols)               ! Ocean fraction
         rhu00   (1:pcols,1:pver) , &! !INTENT(out) :: rhu00  (pcols,pver)      ! RH threshold for cloud
         relhum  (1:pcols,1:pver) , &! !INTENT(out) :: relhum (pcols,pver)      ! RH 
         0                          )! !INTENT(in) :: dindex                     ! 0 or 1 to perturb rh

    ! re-calculate cloud with perturbed rh             add call cldfrc  
    
    CALL cldfrc(                    &
         pcols                   , & !INTENT(in) ::  pcols          ! number of columns  (max)
         pver                    , & !INTENT(in) ::  pver ! number of vertical levels
         pverp                   , & !INTENT(in) ::  pverp          ! pver + 1
         pcols                   , & !INTENT(in) :: pcols                   ! number of atmospheric columns
         pnot    (1:pcols)        , & !INTENT(in) :: pnot    (pcols)               ! reference pressure
         pmid    (1:pcols,1:pver) , & !INTENT(in) :: pmid    (pcols,pver)      ! midpoint pressures
         t       (1:pcols,1:pver) , & !INTENT(in) :: temp    (pcols,pver)      ! temperature
         q       (1:pcols,1:pver) , & !INTENT(in) :: q            (pcols,pver)      ! specific humidity
         phis    (1:pcols)        , & !INTENT(in) :: phis    (pcols)              ! surface geopotential 
         cld2    (1:pcols,1:pver) , & !INTENT(out) :: cloud  (pcols,pver)      ! cloud fraction
         clc     (1:pcols)        , & !INTENT(out) :: clc    (pcols)              ! column convective cloud amount 
         cmfmc   (1:pcols,1:pverp), & !INTENT(in) :: cmfmc   (pcols,pverp)     ! convective mass flux--m sub c
         cmfmc2  (1:pcols,1:pverp), & !INTENT(in) :: cmfmc2  (pcols,pverp)     ! shallow convective mass flux--m sub c
         landfrac(1:pcols)        , & !INTENT(in) :: landfrac(pcols)              ! Land fraction 
         snowh   (1:pcols)        , & !INTENT(in) :: snowh   (pcols)              ! snow depth (liquid water equivalent) ! sea surface temperature
         concld2 (1:pcols,1:pver) , & !INTENT(out) concld (pcols,pver)   ! convective cloud cover
         cldst2  (1:pcols,1:pver) , & !INTENT(in) :: cldst2  (pcols,pver)      ! not used ! detrainment rate from deep convection
         sst     (1:pcols)        , & !INTENT(in) :: sst   (pcols)            ! sea surface temperature
         ps      (1:pcols)        , & !INTENT(in) :: ps    (pcols)            ! surface pressure
         !zdu     (1:pcols,1:pver) , & !INTENT(in) :: zdu   (pcols,pver)      ! not used ! detrainment rate from deep convection
         ocnfrac (1:pcols)        , & !INTENT(in) :: ocnfrac (pcols)              ! Ocean fraction
         rhu002  (1:pcols,1:pver) , & !INTENT(out) :: rhu00  (pcols,pver)      ! RH threshold for cloud
         relhum2 (1:pcols,1:pver) , & !INTENT(out) :: relhum (pcols,pver)      ! RH 
         1                         ) !INTENT(in) :: dindex                   ! 0 or 1 to perturb rh

    ! cldfrc does not define layer cloud for model layer at k=1
    ! so set rhu00(k=1)=2.0 to not calculate cme for this layer
    DO i=1,pcols
       rhu00(i,1) = 2.0_r8 
    END DO
    ! Add following to estimate rhdfda                       
 
    DO k=1,pver
       DO i=1,pcols
          IF(relhum(i,k) < rhu00(i,k) ) THEN
             rhdfda(i,k)=0.0_r8
          ELSE IF (relhum(i,k) >= 1.0_r8 ) THEN
             rhdfda(i,k)=0.0_r8
          ELSE
             !under certain circumstances, rh+ cause cld not to changed
             !when at an upper limit, or w/ strong subsidence
             !need to further check whether this if-block is necessary

             IF((cld2(i,k) - cld(i,k) ) < 1.0e-4_r8 ) THEN
                rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8   !instead of 0.0
             ELSE
                rhdfda(i,k)=0.01_r8*relhum(i,k)/(cld2(i,k)-cld(i,k))
             ENDIF
          ENDIF
       ENDDO
    ENDDO


    !+ mp +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! cloud water and ice parameterizations
    !drb   prec_pcw = 0.
    !drb   snow_pcw = 0.


    ! Initialize chunk id and size

    ! associate local pointers with fields in the physics buffer
!!$    buffld1 => pbuf(ixbuffld1)%fld_ptr(1,1:pcols,1,     lchnk,1)
!!$    buffld2 => pbuf(ixbuffld2)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

    ! Define fractional amount of cloud condensate in ice phase
    
    CALL cldwat_fice(pcols                          , &!INTENT(in)  :: pcols                  ! number of active columns
                    pcols                         , &!INTENT(in) :: pcols              ! number of columns (max)
                    pver                          , &!INTENT(in) :: pver              ! number of vertical levels
                    state1%t(1:pcols,1:pver,latco) , &!INTENT(in)  :: t(pcols,pver)          ! temperature
                    fice    (1:pcols,1:pver)       , &!INTENT(out) :: fice(pcols,pver)          ! Fractional ice content within cloud
                    fsnow   (1:pcols,1:pver)         )!INTENT(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection
    DO i = 1,pcols
       DO k = 1,pver

          ! compute total cloud water
          totcw(i,k) = state1%q(i,k,latco,ixcldice) + state1%q(i,k,latco,ixcldliq)

          ! save input cloud ice
          repartht(i,k) = state1%q(i,k,latco,ixcldice)

          ! Repartition ice and liquid
          !  state1%q(i,k,ixcldice) = totcw(i,k) * fice(i,k)
          !  state1%q(i,k,ixcldliq) = totcw(i,k) * (1.0_r8 - fice(i,k))
          rdtime = 1.0_r8/dtime
          ptend_loc%q(i,k,latco,ixcldice) = &
               ( totcw(i,k) * fice(i,k)- state1%q(i,k,latco,ixcldice) ) * rdtime
          ptend_loc%q(i,k,latco,ixcldliq) = &
               ( totcw(i,k) * (1.0_r8 - fice(i,k)) - state1%q(i,k,latco,ixcldliq) ) * rdtime
       END DO
    END DO


    ! Set output flags
    ptend_loc%name(latco)         = 'cldwat-repartition'
    ptend_loc%lq(latco,ixcldice) = .TRUE.
    ptend_loc%lq(latco,ixcldliq) = .TRUE.

    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, ppcnst,latco,pcols)

    ! update state for use below
    CALL physics_update (state1, tend, ptend_loc, dtime,ppcnst,pcols,pver,latco)
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)


    ! Set output flags for final cldwat update 
    ptend_loc%name(latco)         = 'cldwat'
    ptend_loc%ls(latco)           = .TRUE.
    ptend_loc%lq(latco,1)        = .TRUE.
    ptend_loc%lq(latco,ixcldice) = .TRUE.
    ptend_loc%lq(latco,ixcldliq) = .TRUE.

    DO i = 1,pcols
       DO k = 1,pver

          ! Determine heating from change in cloud ice
          repartht(i,k) = latice/dtime * (state1%q(i,k,latco,ixcldice)&
                                    -repartht(i,k))

          ! calculate the tendencies for moisture, temperature and cloud fraction
          qtend(i,k) = (state1%q(i,k,latco,1) - qcwat(i,k))*rdtime
          ttend(i,k) = (state1%t(i,k,latco)   - tcwat(i,k))*rdtime
          ltend(i,k) = (totcw   (i,k)         - lcwat(i,k))*rdtime
       END DO
    END DO


    ! call microphysics package to calculate tendencies
    CALL pcond (&
         pcols    ,&                                  !INTEGER, INTENT(in) :: pcols                  ! number of atmospheric columns
         pcols  ,&                                  !INTEGER, INTENT(in) :: pcols                ! number of columns (max)
         pverp  ,&                                  !INTEGER, INTENT(in) :: pverp                ! number of vertical levels + 1
         pver    , &                                  !INTEGER, INTENT(in) :: pver                 ! number of vertical levels
         state1%t    (1:pcols,1:pver,latco)  , &   !REAL(r8), INTENT(in ) :: tn(pcols,pver)         ! new temperature    (K)
         ttend       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: ttend(pcols,pver)         ! temp tendencies    (K/s)
         state1%q    (1:pcols,1:pver,latco,1), &   !REAL(r8), INTENT(in ) :: qn (pcols,pver)        ! new water vapor    (kg/kg)
         qtend       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: qtend(pcols,pver)        ! mixing ratio tend  (kg/kg/s)
         totcw       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: cwat(pcols,pver)        ! cloud water (kg/kg)
         state1%pmid (1:pcols,1:pver,latco)  , &   !REAL(r8), INTENT(in ) :: p(pcols,pver)        ! pressure           (K)
         state1%pdel (1:pcols,1:pver,latco)  , &   !REAL(r8), INTENT(in ) :: pdel(pcols,pver)        ! pressure thickness (Pa)
         cld         (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: cldn(pcols,pver)        ! new value of cloud fraction         (fraction)
         fice        (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: fice(pcols,pver)        ! fraction of cwat that is ice
         fsnow       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: fsnow(pcols,pver)        ! fraction of rain that freezes to snow
         qme         (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: cme     (pcols,pver) ! rate of cond-evap of condensate (1/s)
         prain       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: prodprec(pcols,pver) ! rate of conversion of condensate to precip (1/s)
         prodsnow    (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: prodsnow(pcols,pver) ! rate of production of snow
         nevapr      (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: evapprec(pcols,pver) ! rate of evaporation of falling precip (1/s)
         evapsnow    (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: evapsnow(pcols,pver) ! rate of evaporation of falling snow (1/s)
         evapheat    (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: evapheat(pcols,pver) ! heating rate due to evaporation of precip (W/kg)
         prfzheat    (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: prfzheat(pcols,pver) ! heating rate due to freezing of precip (W/kg)
         meltheat    (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: meltheat(pcols,pver) ! heating rate due to snow melt (W/kg)
         prec_pcw    (1:pcols)               , &   !REAL(r8), INTENT(out) :: precip(pcols)         ! rate of precipitation (kg / (m**2 * s))
         snow_pcw    (1:pcols)               , &   !REAL(r8), INTENT(out) :: snowab(pcols)         ! rate of snow (kg / (m**2 * s))
         dtime                               , &   !REAL(r8), INTENT(in ) :: deltat                ! time step to advance solution over
         fwaut       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: fwaut(pcols,pver)         ! relative importance of warm cloud autoconversion    fsaci(1:pcols,:) = 0.0_r8
         fsaut       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: fsaut(pcols,pver)         ! relative importance of ice auto conversion            fsacw(1:pcols,:) = 0.0_r8
         fracw       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: fracw(pcols,pver)         ! relative importance of collection of liquid by rain fwaut(1:pcols,:) = 0.0_r8
         fsacw       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: fsacw(pcols,pver)         ! relative importance of collection of liquid by snow fracw(1:pcols,:) = 0.0_r8
         fsaci       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: fsaci(pcols,pver)         ! relative importance of collection of ice by snow    fsaut(1:pcols,:) = 0.0_r8
         ltend       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: lctend(pcols,pver)        ! cloud liquid water tendencies   ====wlin
         rhdfda      (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: rhdfda(pcols,pver)        ! dG(a)/da, rh=G(a), when rh>u00  ====wlin
         rhu00       (1:pcols,1:pver)        , &   !REAL(r8), INTENT(in ) :: rhu00 (pcols,pver)        ! Rhlim for cloud                  ====wlin
         icefrac     (1:pcols)               , &   !REAL(r8), INTENT(in ) :: seaicef(pcols)        ! sea ice fraction  (fraction)
         state1%zi   (1:pcols,1:pver+1,latco), &   !REAL(r8), INTENT(in ) :: zi(pcols,pverp)        ! layer interfaces (m)
         ice2pr      (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
         liq2pr      (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
         liq2snow    (1:pcols,1:pver)        , &   !REAL(r8), INTENT(out) :: liq2snow(pcols,pver) ! rate of conversion of liquid to snow
         snowh       (1:pcols)               , &   !REAL(r8), INTENT(in ) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
         landm       (1:pcols)               , &   !REAL(r8), INTENT(in ) :: landm(pcols)         ! land fraction ramped over water
         landfrac    (1:pcols)              )    !INTENT(in) :: landfrac(pcols)               ! Land fraction 
!RETURN
    ! make it interactive
    DO i = 1,pcols
       DO k = 1,pver
          ptend_loc%s(i,k,latco) = qme(i,k) * (latvap + latice*fice(i,k)) &
               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)

          ptend_loc%q(i,k,latco,1) =-qme(i,k) + nevapr(i,k)

          ptend_loc%q(i,k,latco,ixcldice) =qme(i,k)*fice(i,k) - ice2pr(i,k)
          ptend_loc%q(i,k,latco,ixcldliq) =qme(i,k)*(1.0_r8-fice(i,k)) - liq2pr(i,k)
       END DO
    END DO

    !#ifdef DEBUG
    !  if (lchnk.eq.248) then
    !     i = 12
    !1     do k = 1,pver
    !        call debug_microphys_1(state1,ptend_loc,i,k, &
    !                dtime,qme,fice,snow_pcw,prec_pcw, &
    !                prain,nevapr,prodsnow, evapsnow, &
    !                ice2pr,liq2pr,liq2snow)
    !     end do
    !  endif
    !  call debug_microphys_2(state1,&
    !       snow_pcw,fsaut,fsacw ,fsaci, meltheat)
    !#endif

    ! Compute in cloud ice and liquid mixing ratios

    DO k=1,pver
       DO i = 1,pcols
          icimr(i,k) = (state1%q(i,k,latco,ixcldice) + dtime*ptend_loc%q(i,k,latco,ixcldice)) / MAX(0.01_r8,cld(i,k))
          icwmr(i,k) = (state1%q(i,k,latco,ixcldliq) + dtime*ptend_loc%q(i,k,latco,ixcldliq)) / MAX(0.01_r8,cld(i,k))
       END DO
    END DO


    ! convert precipitation from kg/m2 to m/s
    DO i = 1,pcols
       snow_pcw  (i) = snow_pcw  (i)/1000.0_r8
       prec_pcw  (i) = prec_pcw  (i)/1000.0_r8
    END DO
    !DO i = 1,pcols
    !   PRINT*, snow_pcw  (i) ,prec_pcw  (i)
    !END DO   
    DO i = 1,pcols
       DO k = 1,pver
          cmeheat(i,k) = qme(i,k) * (latvap + latice*fice(i,k))
          cmeice (i,k) = qme(i,k) * fice(i,k)
          cmeliq (i,k) = qme(i,k) * (1.0_r8 - fice(i,k))
       END DO
    END DO

    ! update boundary quantities
!!$    ptend_loc%hflx_srf = 0.
!!$    ptend_loc%hflx_top = 0.
!!$    ptend_loc%cflx_srf = 0.
!!$    ptend_loc%cflx_top = 0.


    ! add tendency from this process to tend from other processes here
    CALL physics_ptend_sum(ptend_loc,ptend_all, ppcnst,latco,pcols)

    ! Set the name of the final package tendencies. Note that this
    ! is a special case in physics_update, so a change here must be 
    ! matched there.
    ptend_all%name(latco) = 'stratiform'


    ! used below
    CALL physics_update (state1, tend, ptend_loc, dtime,ppcnst,pcols,pver,latco)
    CALL physics_ptend_init(ptend_loc,latco,ppcnst,pver)

    !call t_stopf('stratiform_microphys')

    ! accumulate prec and snow
    DO i = 1,pcols
       prec_str(i) = prec_str(i) + prec_pcw(i)
       snow_str(i) = snow_str(i) + snow_pcw(i)
    END DO
    ! Save off q and t after cloud water
    !call cnst_get_ind('CLDLIQ', ixcldliq)
    !call cnst_get_ind('CLDICE', ixcldice)

    DO k=1,pver
       DO i=1,pcols
          qcwat(i,k) = state1%q(i,k,latco,1)
          tcwat(i,k) = state1%t(i,k,latco)
          lcwat(i,k) = state1%q(i,k,latco,ixcldice) + state1%q(i,k,latco,ixcldliq)
       END DO 
    END DO

    !
    ! Cloud water and ice particle sizes, saved in physics buffer for radiation
    !
    CALL cldefr(&
         pcols                            , & !INTEGER, INTENT(in) :: pcols                    ! number of atmospheric columns
         pcols                            , & !INTEGER, INTENT(in) :: pcols                    ! number of columns (max)
         pver                             , & !INTEGER, INTENT(in) :: pver                    ! number of vertical levels
         state1%t   (1:pcols,1:pver,latco), & !REAL(r8), INTENT(in) :: t       (pcols,pver) ! Temperature
         rel        (1:pcols,1:pver)      , & !REAL(r8), INTENT(out) :: rel(pcols,pver)     ! Liquid effective drop size (microns)
         rei        (1:pcols,1:pver)      , & !REAL(r8), INTENT(out) :: rei(pcols,pver)     ! Ice effective drop size (microns)
         landm      (1:pcols)             , & !REAL(r8), INTENT(in) :: landm   (pcols)
         icefrac    (1:pcols)             , & !REAL(r8), INTENT(in) :: icefrac (pcols)      ! Ice fraction
         snowh      (1:pcols)               ) !REAL(r8), INTENT(in) :: snowh   (pcols)      ! Snow depth over land, water equivalent (m)

   call physics_update (state, tend, ptend_all, dtime,ppcnst,pcols,pver,latco)

    DO k=1,pver
       DO i=1,pcols          

          dtdt      (i,pver+1-k) = ( state1%t(i,k,latco)         - state_t   (i,pver+1-k))/dtime
          dqdt      (i,pver+1-k) = ( state1%q(i,k,latco,1)       - state_qv  (i,pver+1-k))/dtime
          dqldt     (i,pver+1-k) = ( state1%q(i,k,latco,ixcldliq)- state_ql  (i,pver+1-k))/dtime
          dqidt     (i,pver+1-k) = ( state1%q(i,k,latco,ixcldice)- state_qi  (i,pver+1-k))/dtime

          state_t   (i,pver+1-k)= state1%t(i,k,latco)
          state_qv  (i,pver+1-k)= state1%q(i,k,latco,1)
          state_ql  (i,pver+1-k)= state1%q(i,k,latco,ixcldliq)
          state_qi  (i,pver+1-k)= state1%q(i,k,latco,ixcldice)
          
          state_concld(i,k      )=concld(i,pver+1-k) 
          !state_qcwat (i,k,latco)=qcwat (i,pver+1-k) ! cloud water old q
          !state_tcwat (i,k,latco)=tcwat (i,pver+1-k) ! cloud water old temperature
          !state_lcwat (i,k,latco)=lcwat (i,pver+1-k) ! cloud liquid water old q
          state_cld   (i,k      )=cld   (i,pver+1-k) ! cloud fraction
          !state_qme   (i,k,latco)=qme   (i,pver+1-k) 
          !state_prain (i,k,latco)=prain (i,pver+1-k) 
          !state_nevapr(i,k,latco)=nevapr(i,pver+1-k) 
          !state_rel   (i,k,latco)=rel   (i,pver+1-k) ! liquid effective drop radius (microns)
          !state_rei   (i,k,latco)=rei   (i,pver+1-k) ! ice effective drop size (microns)
       END DO
    END DO

  END SUBROUTINE RunMicro_Hack
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

    tend%dtdt   (:,:,latco)   = 0.0_r8
    tend%dudt   (:,:,latco)   = 0.0_r8
    tend%dvdt   (:,:,latco)   = 0.0_r8
    tend%flx_net(:  ,latco)   = 0.0_r8
    tend%te_tnd (:  ,latco)   = 0.0_r8
    tend%tw_tnd (:  ,latco)   = 0.0_r8

  END SUBROUTINE physics_tend_init

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
          CALL qneg3(  pcols, pcols, pver, m, m, qmin(m), qq(1:pcols,1:pver,m))
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
       CALL geopotential_dse(&
            state%lnpint(1:pcols,1:pver+1,latco),  state%pint (1:pcols,1:pver+1,latco)  , &
            state%pmid  (1:pcols,1:pver,latco)  ,  state%pdel  (1:pcols,1:pver,latco)   , &
            state%rpdel(1:pcols,1:pver,latco)   ,  state%s     (1:pcols,1:pver,latco)   , &
            state%q     (1:pcols,1:pver,latco,1), state%phis (1:pcols,latco) , rair  , &
            gravit      , cpair        ,zvir        , &
            state%t(1:pcols,1:pver,latco)     , state%zi(1:pcols,1:pver+1,latco)    ,&
            state%zm(1:pcols,1:pver,latco),&
            pcols, pver, pver+1          )
    END IF

    ! Reset all parameterization tendency flags to false
    CALL physics_ptend_reset(ptend,latco,ppcnst,pver)

  END SUBROUTINE physics_update


  SUBROUTINE qneg3 (pcols    ,ncold   ,lver    ,lconst_beg  , &
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
    !CHARACTER(LEN=*), INTENT(in) :: subnam ! name of calling routine

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
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!DIR$ preferstream
       DO k=1,lver
          nval(k) = 0
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!DIR$ prefervector
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
!9000 FORMAT(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i3, &
!         ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
!         1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
  END SUBROUTINE qneg3


  !===============================================================================
  SUBROUTINE physics_ptend_sum(ptend, ptend_sum, ppcnst,latco,pcols)
    !-----------------------------------------------------------------------
    ! Add ptend fields to ptend_sum for ptend logical flags = .true.
    ! Where ptend logical flags = .false, don't change ptend_sum
    !-----------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN   ) ::ppcnst,latco,pcols
    TYPE(physics_ptend), INTENT(in)     :: ptend   ! New parameterization tendencies
    TYPE(physics_ptend), INTENT(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    !TYPE(physics_state), INTENT(in)     :: state   ! New parameterization tendencies

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
                ptend_sum%qt(i,k,latco)  = ptend_sum%qt(i,k,latco)  + ptend%q(i,k,latco,m)
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
       ptend%s        (:,:,latco) = 0.0_r8
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
          ptend%qt        (:,:,latco) = 0.0_r8
          ptend%q      (:,:,latco,m) = 0.0_r8
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

  SUBROUTINE physics_state_copy(pcols,latco, pver, pverp,ppcnst)

    !use ppgrid,       only: pver, pverp
    !use constituents, only: ppcnst, cnst_need_pdeldry

    IMPLICIT NONE

    !
    ! Arguments
    !
    INTEGER, INTENT(IN   ) :: pver, pverp,ppcnst,latco,pcols
!    LOGICAL, INTENT(IN   ) ::  cnst_need_pdeldry
!    TYPE(physics_state), INTENT(in) :: state
!    TYPE(physics_state), INTENT(out) :: state1

    !
    ! Local variables
    !
    INTEGER i, k, m 

    state1%ncol(latco)  = pcols
    state1%count(latco) = state%count (latco)

    DO i = 1, pcols
       state1%lat(i,latco)    = state%lat(i,latco)
       state1%lon(i,latco)    = state%lon(i,latco)
       state1%ps(i,latco)     = state%ps(i,latco)
       state1%phis(i,latco)   = state%phis(i,latco)
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

  !MODULE CLDWAT

  SUBROUTINE pcond (&
       ncol    , &!INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
       pcols   , &!INTEGER, INTENT(in) :: pcols                ! number of columns (max)
       pverp   , &!INTEGER, INTENT(in) :: pverp                ! number of vertical levels + 1
       pver    , &!INTEGER, INTENT(in) :: pver                 ! number of vertical levels
       tn      , &!REAL(r8), INTENT(in) :: tn(pcols,pver)         ! new temperature    (K)
       ttend   , &!REAL(r8), INTENT(in) :: ttend(pcols,pver)         ! temp tendencies    (K/s)
       qn      , &!REAL(r8), INTENT(in) :: qn (pcols,pver)        ! new water vapor    (kg/kg)
       qtend   , &!REAL(r8), INTENT(in) :: qtend(pcols,pver)    ! mixing ratio tend  (kg/kg/s)
       cwat    , &!REAL(r8), INTENT(in) :: cwat(pcols,pver)     ! cloud water (kg/kg)
       p       , &!REAL(r8), INTENT(in) :: p(pcols,pver)        ! pressure           (K)
       pdel    , &!REAL(r8), INTENT(in) :: pdel(pcols,pver)     ! pressure thickness (Pa)
       cldn    , &!REAL(r8), INTENT(in) :: cldn(pcols,pver)     ! new value of cloud fraction    (fraction)
       fice    , &!REAL(r8), INTENT(in) :: fice(pcols,pver)        ! fraction of cwat that is ice
       fsnow   , &!REAL(r8), INTENT(in) :: fsnow(pcols,pver)        ! fraction of rain that freezes to snow
       cme     , &!REAL(r8), INTENT(out) :: cme     (pcols,pver) ! rate of cond-evap of condensate (1/s)
       prodprec, &!REAL(r8), INTENT(out) :: prodprec(pcols,pver) ! rate of conversion of condensate to precip (1/s)
       prodsnow, &!REAL(r8), INTENT(out) :: prodsnow(pcols,pver) ! rate of production of snow
       evapprec, &!REAL(r8), INTENT(out) :: evapprec(pcols,pver) ! rate of evaporation of falling precip (1/s)
       evapsnow, &!REAL(r8), INTENT(out) :: evapsnow(pcols,pver) ! rate of evaporation of falling snow (1/s)
       evapheat, &!REAL(r8), INTENT(out) :: evapheat(pcols,pver) ! heating rate due to evaporation of precip (W/kg)
       prfzheat, &!REAL(r8), INTENT(out) :: prfzheat(pcols,pver) ! heating rate due to freezing of precip (W/kg)
       meltheat, &!REAL(r8), INTENT(out) :: meltheat(pcols,pver) ! heating rate due to snow melt (W/kg)
       precip  , &!REAL(r8), INTENT(out) :: precip(pcols)        ! rate of precipitation (kg / (m**2 * s))
       snowab  , &!REAL(r8), INTENT(out) :: snowab(pcols)        ! rate of snow (kg / (m**2 * s))
       deltat  , &!REAL(r8), INTENT(in) :: deltat               ! time step to advance solution over
       fwaut   , &!REAL(r8), INTENT(out) :: fwaut(pcols,pver)    ! relative importance of warm cloud autoconversion    fsaci(1:ncol,:) = 0.0_r8
       fsaut   , &!REAL(r8), INTENT(out) :: fsaut(pcols,pver)    ! relative importance of ice auto conversion            fsacw(1:ncol,:) = 0.0_r8
       fracw   , &!REAL(r8), INTENT(out) :: fracw(pcols,pver)    ! relative importance of collection of liquid by rain fwaut(1:ncol,:) = 0.0_r8
       fsacw   , &!REAL(r8), INTENT(out) :: fsacw(pcols,pver)    ! relative importance of collection of liquid by snow fracw(1:ncol,:) = 0.0_r8
       fsaci   , &!REAL(r8), INTENT(out) :: fsaci(pcols,pver)    ! relative importance of collection of ice by snow    fsaut(1:ncol,:) = 0.0_r8
       lctend  , &!REAL(r8), INTENT(in) :: lctend(pcols,pver)   ! cloud liquid water tendencies   ====wlin
       rhdfda  , &!REAL(r8), INTENT(in) :: rhdfda(pcols,pver)   ! dG(a)/da, rh=G(a), when rh>u00  ====wlin
       rhu00   , &!REAL(r8), INTENT(in) :: rhu00 (pcols,pver)   ! Rhlim for cloud                 ====wlin
       seaicef , &!REAL(r8), INTENT(in) :: seaicef(pcols)       ! sea ice fraction  (fraction)
       zi      , &!REAL(r8), INTENT(in) :: zi(pcols,pverp)      ! layer interfaces (m)
       ice2pr  , &!REAL(r8), INTENT(out) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
       liq2pr  , &!REAL(r8), INTENT(out) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
       liq2snow, &!REAL(r8), INTENT(out) :: liq2snow(pcols,pver) ! rate of conversion of liquid to snow
       snowh   , &!REAL(r8), INTENT(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
       landm   , &!REAL(r8), INTENT(in) :: landm(pcols)         ! land fraction ramped over water
       landfrac)! !INTENT(in) :: landfrac(pcols)               ! Land fraction 
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! The public interface to the cloud water parameterization
    ! returns tendencies to water vapor, temperature and cloud water variables
    ! 
    ! For basic method 
    !  See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
    !  model climate using diagnosed and 
    !  predicted condensate parameterizations, 1998, J. Clim., 11,
    !  pp1587---1614.
    ! 
    ! For important modifications to improve the method of determining
    ! condensation/evaporation see Zhang et al (2001, in preparation)
    !
    ! Authors: M. Zhang, W. Lin, P. Rasch and J.E. Kristjansson
    !          B. A. Boville (latent heat of fusion)
    !-----------------------------------------------------------------------
    !   use wv_saturation, only: vqsatd
    !
    !---------------------------------------------------------------------
    !
    ! Input Arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                 ! number of vertical levels + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels


    REAL(r8), INTENT(in) :: fice(pcols,pver)     ! fraction of cwat that is ice
    REAL(r8), INTENT(in) :: fsnow(pcols,pver)    ! fraction of rain that freezes to snow
    REAL(r8), INTENT(in) :: cldn(pcols,pver)     ! new value of cloud fraction    (fraction)
    REAL(r8), INTENT(in) :: cwat(pcols,pver)     ! cloud water (kg/kg)
    REAL(r8), INTENT(in) :: p(pcols,pver)        ! pressure          (K)
    REAL(r8), INTENT(in) :: pdel(pcols,pver)     ! pressure thickness (Pa)
    REAL(r8), INTENT(in) :: qn (pcols,pver)       ! new water vapor    (kg/kg)
    REAL(r8), INTENT(in) :: qtend(pcols,pver)    ! mixing ratio tend  (kg/kg/s)
    REAL(r8), INTENT(in) :: tn(pcols,pver)       ! new temperature    (K)
    REAL(r8), INTENT(in) :: ttend(pcols,pver)    ! temp tendencies    (K/s)
    REAL(r8), INTENT(in) :: deltat               ! time step to advance solution over
    REAL(r8), INTENT(in) :: lctend(pcols,pver)   ! cloud liquid water tendencies   ====wlin
    REAL(r8), INTENT(in) :: rhdfda(pcols,pver)   ! dG(a)/da, rh=G(a), when rh>u00  ====wlin
    REAL(r8), INTENT(in) :: rhu00 (pcols,pver)   ! Rhlim for cloud                 ====wlin
    REAL(r8), INTENT(in) :: seaicef(pcols)       ! sea ice fraction  (fraction)
    REAL(r8), INTENT(in) :: zi(pcols,pverp)      ! layer interfaces (m)
    REAL(r8), INTENT(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    REAL(r8), INTENT(in) :: landm(pcols)         ! land fraction ramped over water
    REAL(r8), INTENT(in) :: landfrac(pcols)         ! land fraction ramped over water

    ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    !
    ! Output Arguments
    !
    REAL(r8), INTENT(out) :: cme     (pcols,pver) ! rate of cond-evap of condensate (1/s)
    REAL(r8), INTENT(out) :: prodprec(pcols,pver) ! rate of conversion of condensate to precip (1/s)
    REAL(r8), INTENT(out) :: evapprec(pcols,pver) ! rate of evaporation of falling precip (1/s)
    REAL(r8), INTENT(out) :: evapsnow(pcols,pver) ! rate of evaporation of falling snow (1/s)
    REAL(r8), INTENT(out) :: evapheat(pcols,pver) ! heating rate due to evaporation of precip (W/kg)
    REAL(r8), INTENT(out) :: prfzheat(pcols,pver) ! heating rate due to freezing of precip (W/kg)
    REAL(r8), INTENT(out) :: meltheat(pcols,pver) ! heating rate due to snow melt (W/kg)
    REAL(r8), INTENT(out) :: precip  (pcols)        ! rate of precipitation (kg / (m**2 * s))
    REAL(r8), INTENT(out) :: snowab  (pcols)        ! rate of snow (kg / (m**2 * s))
    !++pjr
    REAL(r8), INTENT(out) :: ice2pr  (pcols,pver)   ! rate of conversion of ice to precip
    REAL(r8), INTENT(out) :: liq2pr  (pcols,pver)   ! rate of conversion of liquid to precip
    REAL(r8), INTENT(out) :: liq2snow(pcols,pver) ! rate of conversion of liquid to snow

    REAL(r8) nice2pr     ! rate of conversion of ice to snow
    REAL(r8) nliq2pr     ! rate of conversion of liquid to precip
    REAL(r8) nliq2snow   ! rate of conversion of liquid to snow
    REAL(r8), INTENT(out) :: prodsnow(pcols,pver) ! rate of production of snow
    !--pjr

    !
    ! Local workspace
    !
    REAL(r8) :: precab(pcols)        ! rate of precipitation (kg / (m**2 * s))
    INTEGER i                 ! work variable
    INTEGER iter              ! #iterations for precipitation calculation
    INTEGER k                 ! work variable
    INTEGER l                 ! work variable

    REAL(r8) cldm  (pcols)          ! mean cloud fraction over the time step
    REAL(r8) cldmax(pcols)        ! max cloud fraction above
    REAL(r8) coef  (pcols)          ! conversion time scale for condensate to rain
    REAL(r8) cwm   (pcols)           ! cwat mixing ratio at midpoint of time step
    REAL(r8) cwn   (pcols)           ! cwat mixing ratio at end
    REAL(r8) denom                ! work variable
    REAL(r8) dqsdt                ! change in sat spec. hum. wrt temperature
    REAL(r8) es    (pcols)            ! sat. vapor pressure
    REAL(r8), INTENT(out) :: fracw(pcols,pver)    ! relative importance of collection of liquid by rain fwaut(1:ncol,:) = 0.0_r8
    REAL(r8), INTENT(out) :: fsaci(pcols,pver)    ! relative importance of collection of ice by snow        fsaut(1:ncol,:) = 0.0_r8
    REAL(r8), INTENT(out) :: fsacw(pcols,pver)    ! relative importance of collection of liquid by snow fracw(1:ncol,:) = 0.0_r8
    REAL(r8), INTENT(out) :: fsaut(pcols,pver)    ! relative importance of ice auto conversion                fsacw(1:ncol,:) = 0.0_r8
    REAL(r8), INTENT(out) :: fwaut(pcols,pver)    ! relative importance of warm cloud autoconversion        fsaci(1:ncol,:) = 0.0_r8
    REAL(r8) gamma (pcols)         ! d qs / dT
    REAL(r8) icwc  (pcols)          ! in-cloud water content (kg/kg)
    REAL(r8) mincld               ! a small cloud fraction to avoid / zero
    REAL(r8) omeps                ! 1 minus epsilon
    REAL(r8),PARAMETER ::omsm=0.99999_r8                 ! a number just less than unity (for rounding)
    REAL(r8) prprov(pcols)        ! provisional value of precip at btm of layer
    REAL(r8) prtmp                ! work variable
    REAL(r8) q     (pcols,pver)        ! mixing ratio before time step ignoring condensate
    REAL(r8) qs    (pcols)            ! spec. hum. of water vapor
    REAL(r8) qsn, esn             ! work variable
    REAL(r8) qsp   (pcols,pver)      ! sat pt mixing ratio
    REAL(r8) qtl   (pcols)           ! tendency which would saturate the grid box in deltat
    REAL(r8) qtmp, ttmp           ! work variable
    REAL(r8) relhum1(pcols)        ! relative humidity
    REAL(r8) relhum(pcols)        ! relative humidity
!!$   real(r8) tc                   ! crit temp of transition to ice
    REAL(r8) t(pcols,pver)        ! temp before time step ignoring condensate
    REAL(r8) tsp(pcols,pver)      ! sat pt temperature
    REAL(r8) pol                  ! work variable
    REAL(r8) cdt                  ! work variable
    REAL(r8) wtthick              ! work variable

    ! Extra local work space for cloud scheme modification       

    REAL(r8) cpohl                !Cp/Hlatv
    REAL(r8) hlocp                !Hlatv/Cp
    REAL(r8) dto2                 !0.5*deltat (delta=2.0*dt)
    REAL(r8) calpha(pcols)        !alpha of new C - E scheme formulation
    REAL(r8) cbeta (pcols)        !beta  of new C - E scheme formulation
    REAL(r8) cbetah(pcols)        !beta_hat at saturation portion 
    REAL(r8) cgamma(pcols)        !gamma of new C - E scheme formulation
    REAL(r8) cgamah(pcols)        !gamma_hat at saturation portion
    REAL(r8) rcgama(pcols)        !gamma/gamma_hat
    REAL(r8) csigma(pcols)        !sigma of new C - E scheme formulation
    REAL(r8) cmec1 (pcols)        !c1    of new C - E scheme formulation
    REAL(r8) cmec2 (pcols)        !c2    of new C - E scheme formulation
    REAL(r8) cmec3 (pcols)        !c3    of new C - E scheme formulation
    REAL(r8) cmec4 (pcols)        !c4    of new C - E scheme formulation
    REAL(r8) cmeres(pcols)        !residual cond of over-sat after cme and evapprec
    REAL(r8) ctmp                 !a scalar representation of cmeres
    REAL(r8) clrh2o               ! Ratio of latvap to water vapor gas const
    !++pjr
    !REAL(r8) ice(pcols,pver)    ! ice mixing ratio
    !REAL(r8) liq(pcols,pver)    ! liquid mixing ratio
    REAL(r8) rcwn(pcols,iterp,pver), rliq(pcols,iterp,pver), rice(pcols,iterp,pver)
    REAL(r8) cwnsave(pcols,iterp,pver), cmesave(pcols,iterp,pver)
    REAL(r8) prodprecsave(pcols,iterp,pver)
    !--pjr
    LOGICAL error_found
    REAL(r8) epsqs

    cme=0.0_r8;prodprec=0.0_r8;evapprec=0.0_r8;evapsnow=0.0_r8;evapheat=0.0_r8
    prfzheat=0.0_r8;meltheat=0.0_r8;precip  =0.0_r8;snowab  =0.0_r8;prodsnow=0.0_r8
    ice2pr  =0.0_r8;liq2pr  =0.0_r8;liq2snow=0.0_r8;precab=0.0_r8;cldm  =0.0_r8
    cldmax=0.0_r8;coef  =0.0_r8;cwm   =0.0_r8;cwn   =0.0_r8;denom =0.0_r8;dqsdt =0.0_r8
    es  =0.0_r8;fracw=0.0_r8;fsaci=0.0_r8;fsacw=0.0_r8;fsaut=0.0_r8;fwaut =0.0_r8 
    gamma =0.0_r8;icwc =0.0_r8;prprov =0.0_r8;prtmp =0.0_r8;q     =0.0_r8;qs   =0.0_r8   
    qsn=0.0_r8;esn =0.0_r8;qsp =0.0_r8;qtl =0.0_r8;qtmp=0.0_r8;ttmp =0.0_r8   
    relhum1 =0.0_r8;relhum =0.0_r8 ;t    =0.0_r8   ;tsp  =0.0_r8   ;pol=0.0_r8;cdt=0.0_r8
    wtthick =0.0_r8;cpohl=0.0_r8;hlocp=0.0_r8;dto2=0.0_r8;calpha =0.0_r8;cbeta =0.0_r8
    cbetah=0.0_r8;cgamma=0.0_r8;cgamah=0.0_r8;rcgama=0.0_r8;csigma=0.0_r8;cmec1 =0.0_r8;
    cmec2 =0.0_r8;cmec3 =0.0_r8;cmec4 =0.0_r8;cmeres=0.0_r8;ctmp=0.0_r8;clrh2o =0.0_r8;
    rcwn =0.0_r8;rliq =0.0_r8;rice =0.0_r8;cwnsave=0.0_r8;cmesave=0.0_r8;prodprecsave=0.0_r8   
    epsqs=0.0_r8

    !
    !------------------------------------------------------------
!!!!!!!#include <comadj.h>              
    !------------------------------------------------------------
    !
    epsqs = epsilo
    clrh2o = hlatv/rh2o   ! Ratio of latvap to water vapor gas const
    omeps = 1.0_r8 - epsqs
    !!#ifdef PERGRO
    IF(PERGRO)THEN
       mincld = 1.0e-4_r8
       iter = iterp   ! number of times to iterate the precipitation calculation
    ELSE
       mincld = 1.0e-4_r8
       iter = iterp
    END IF
    !   omsm = 0.99999_r8
    cpohl = cp/hlatv
    hlocp = hlatv/cp
    dto2=0.5_r8*deltat
    !
    ! Constant for computing rate of evaporation of precipitation:
    !
!!$   conke = 1.e-5
!!$   conke = 1.e-6
    !
    ! initialize a few single level fields
    !
    DO i = 1,ncol
       precip(i) = 0.0_r8
       precab(i) = 0.0_r8
       snowab(i) = 0.0_r8
       cldmax(i) = 0.0_r8
    END DO
    !
    ! initialize multi-level fields 
    !
    DO k = 1,pver
       DO i = 1,ncol
          q(i,k) = qn(i,k) 
          t(i,k) = tn(i,k)
       END DO
    END DO
    cme     (1:ncol,:) = 0.0_r8
    evapprec(1:ncol,:) = 0.0_r8
    prodprec(1:ncol,:) = 0.0_r8
    evapsnow(1:ncol,:) = 0.0_r8
    prodsnow(1:ncol,:) = 0.0_r8
    evapheat(1:ncol,:) = 0.0_r8
    meltheat(1:ncol,:) = 0.0_r8
    prfzheat(1:ncol,:) = 0.0_r8
    ice2pr(1:ncol,:)   = 0.0_r8
    liq2pr(1:ncol,:)   = 0.0_r8
    liq2snow(1:ncol,:) = 0.0_r8
    fwaut(1:ncol,:) = 0.0_r8
    fsaut(1:ncol,:) = 0.0_r8
    fracw(1:ncol,:) = 0.0_r8
    fsacw(1:ncol,:) = 0.0_r8
    fsaci(1:ncol,:) = 0.0_r8
    !
    ! find the wet bulb temp and saturation value
    ! for the provisional t and q without condensation
    !
    CALL findsp (ncol,pcols,pver, qn, tn, p, tsp, qsp)
    DO k = k1mb,pver
       CALL vqsatd (t(1:ncol,k), p(1:ncol,k), es, qs, gamma, ncol)
       DO i = 1,ncol
          relhum(i) = q(i,k)/MAX(qs(i),0.0e-12)
          !
          cldm(i) = MAX(cldn(i,k),mincld)
          !
          ! the max cloud fraction above this level
          !
          cldmax(i) = MAX(cldmax(i), cldm(i))

          ! define the coefficients for C - E calculation

          calpha(i) = 1.0_r8/qs(i)
          cbeta (i) = q(i,k)/qs(i)**2*gamma(i)*cpohl
          cbetah(i) = 1.0_r8/qs(i)*gamma(i)*cpohl
          cgamma(i) = calpha(i)+hlatv*cbeta(i)/cp
          cgamah(i) = calpha(i)+hlatv*cbetah(i)/cp
          rcgama(i) = cgamma(i)/cgamah(i)

          IF(cldm(i) > mincld) THEN
             icwc(i) = MAX(0.0_r8,cwat(i,k)/cldm(i))
          ELSE
             icwc(i) = 0.0_r8
          ENDIF

          !
          ! initial guess of evaporation, will be updated within iteration
          !
          IF(landfrac(i) > 0.0_r8 )THEN 
              !LAND 
              evapprec(i,k) = conke_land*(1.0_r8 - cldm(i))*SQRT(precab(i)) &
               *(1.0_r8 - MIN(relhum(i),1.0_r8))
          ELSE
               !OCEAN
               evapprec(i,k) = conke*(1.0_r8 - cldm(i))*SQRT(precab(i)) &
               *(1.0_r8 - MIN(relhum(i),1.0_r8))
          END IF
          !
          ! zero cmeres before iteration for each level
          !
          cmeres(i)=0.0_r8

       END DO
       DO i = 1,ncol
          !
          ! fractions of ice at this level
          !
!!$         tc = t(i,k) - t0
!!$         fice(i,k) = max(0._r8,min(-tc*0.05,1.0_r8))
          !
          ! calculate the cooling due to a phase change of the rainwater
          ! from above
          !
          IF (t(i,k) >= t0) THEN
             meltheat(i,k) =  -hlatf * snowab(i) * gravit/pdel(i,k)
             snowab(i) = 0.0_r8
          ELSE
             meltheat(i,k) = 0.0_r8
          ENDIF
       END DO

       !
       ! calculate cme and formation of precip. 
       !
       ! The cloud microphysics is highly nonlinear and coupled with cme
       ! Both rain processes and cme are calculated iteratively.
       ! 
       DO  l = 1,iter

          DO i = 1,ncol

             !
             ! calculation of cme has 4 scenarios
             ! ==================================
             !
             IF(relhum(i) > rhu00(i,k)) THEN

                ! 1. whole grid saturation
                ! ========================
                IF(relhum(i) >= 0.999_r8 .OR. cldm(i) >= 0.999_r8 ) THEN
                   cme(i,k)=(calpha(i)*qtend(i,k)-cbetah(i)*ttend(i,k))/cgamah(i)

                   ! 2. fractional saturation
                   ! ========================
                ELSE
                   csigma(i) = 1.0_r8/(rhdfda(i,k)+cgamma(i)*icwc(i))
                   cmec1(i) = (1.0_r8-cldm(i))*csigma(i)*rhdfda(i,k)
                   cmec2(i) = cldm(i)*calpha(i)/cgamah(i)+(1.0_r8-rcgama(i)*cldm(i))*   &
                        csigma(i)*calpha(i)*icwc(i)
                   cmec3(i) = cldm(i)*cbetah(i)/cgamah(i) +  &
                        (cbeta(i)-rcgama(i)*cldm(i)*cbetah(i))*csigma(i)*icwc(i)
                   cmec4(i) = csigma(i)*cgamma(i)*icwc(i)

                   ! Q=C-E=-C1*Al + C2*Aq - C3* At + C4*Er

                   cme(i,k) = -cmec1(i)*lctend(i,k) + cmec2(i)*qtend(i,k)  &
                        -cmec3(i)*ttend(i,k) + cmec4(i)*evapprec(i,k)
                ENDIF

                ! 3. when rh < rhu00, evaporate existing cloud water
                ! ================================================== 
             ELSE IF(cwat(i,k) > 0.0)THEN
                ! liquid water should be evaporated but not to exceed 
                ! saturation point. if qn > qsp, not to evaporate cwat
                cme(i,k)=-MIN(MAX(0.0_r8,qsp(i,k)-qn(i,k)),cwat(i,k))/deltat 

                ! 4. no condensation nor evaporation
                ! ==================================
             ELSE
                cme(i,k)=0.0_r8
             ENDIF


          END DO    !end loop for cme update

          ! Because of the finite time step, 
          ! place a bound here not to exceed wet bulb point
          ! and not to evaporate more than available water
          !
          DO i = 1, ncol
             qtmp = qn(i,k) - cme(i,k)*deltat

             ! possibilities to have qtmp > qsp
             !
             !   1. if qn > qs(tn), it condenses; 
             !      if after applying cme,  qtmp > qsp,  more condensation is applied. 
             !      
             !   2. if qn < qs, evaporation should not exceed qsp,

             IF(qtmp > qsp(i,k)) THEN
                cme(i,k) = cme(i,k) + (qtmp-qsp(i,k))/deltat
             ENDIF

             !
             ! if net evaporation, it should not exceed available cwat
             !
             IF(cme(i,k) < -cwat(i,k)/deltat)  &
                  cme(i,k) = -cwat(i,k)/deltat
             !
             ! addition of residual condensation from previous step of iteration
             !
             cme(i,k) = cme(i,k) + cmeres(i)

          END DO

          !      limit cme for roundoff errors
          DO i = 1, ncol
             cme(i,k) = cme(i,k)*omsm
          END DO

          DO i = 1,ncol
             !
             ! as a safe limit, condensation should not reduce grid mean rh below rhu00
             ! 
             IF(cme(i,k) > 0.0_r8 .AND. relhum(i) > rhu00(i,k) )  &
                  cme(i,k) = MIN(cme(i,k), (qn(i,k)-qs(i)*rhu00(i,k))/deltat)
             !
             ! initial guess for cwm (mean cloud water over time step) if 1st iteration
             !
             IF(l < 2) THEN
                cwm(i) = MAX(cwat(i,k)+cme(i,k)*dto2,  0.0_r8)
             ENDIF

          ENDDO

          ! provisional precipitation falling through model layer
          DO i = 1,ncol
!!$            prprov(i) =  precab(i) + prodprec(i,k)*pdel(i,k)/gravit
             ! rain produced in this layer not too effective in collection process
             wtthick = MAX(0.0_r8,MIN(0.5_r8,((zi(i,k)-zi(i,k+1))/1000.0_r8)**2))
             prprov(i) =  precab(i) + wtthick*prodprec(i,k)*pdel(i,k)/gravit
          END DO

          ! calculate conversion of condensate to precipitation by cloud microphysics 
          CALL findmcnew (ncol    , pcols,pver, &
               k       ,prprov  ,  t       ,p        , &
               cwm     ,cldm    ,cldmax  ,fice(1,k),coef    , &
               fwaut(1,k),fsaut(1,k),fracw(1,k),fsacw(1,k),fsaci(1,k), &
               seaicef, snowh,landm)

          !
          ! calculate the precip rate
          ! 
          error_found = .FALSE.
          DO i = 1,ncol
             IF (cldm(i) > 0.0_r8) THEN  
                !
                ! first predict the cloud water
                !
                cdt = coef(i)*deltat
                IF (cdt > 0.01_r8) THEN
                   pol = cme(i,k)/coef(i) ! production over loss
                   cwn(i) = MAX(0.0_r8,(cwat(i,k)-pol)*EXP(-cdt)+ pol)
                ELSE
                   cwn(i) = MAX(0.0_r8,(cwat(i,k) + cme(i,k)*deltat)/(1+cdt))
                ENDIF
                !
                ! now back out the tendency of net rain production
                !
                prodprec(i,k) =  MAX(0.0_r8,cme(i,k)-(cwn(i)-cwat(i,k))/deltat)
             ELSE
                prodprec(i,k) = 0.0_r8
                cwn(i) = 0.0_r8
             ENDIF

             ! provisional calculation of conversion terms
             ice2pr(i,k) = prodprec(i,k)*(fsaut(i,k)+fsaci(i,k))
             liq2pr(i,k) = prodprec(i,k)*(fwaut(i,k)+fsacw(i,k)+fracw(i,k))
             !old        liq2snow(i,k) = prodprec(i,k)*fsacw(i,k)

             !           revision suggested by Jim McCaa
             !           it controls the amount of snow hitting the sfc 
             !           by forcing a lot of conversion of cloud liquid to snow phase
             !           it might be better done later by an explicit representation of 
             !           rain accreting ice (and freezing), or by an explicit freezing of raindrops
             liq2snow(i,k) = MAX(prodprec(i,k)*fsacw(i,k), fsnow(i,k)*liq2pr(i,k))

             ! bounds
             nice2pr = MIN(ice2pr(i,k),(cwat(i,k)+cme(i,k)*deltat)*fice(i,k)/deltat)
             nliq2pr = MIN(liq2pr(i,k),(cwat(i,k)+cme(i,k)*deltat)*(1.0_r8-fice(i,k))/deltat)
             !            write (6,*) ' prodprec ', i, k, prodprec(i,k)
             !            write (6,*) ' nliq2pr, nice2pr ', nliq2pr, nice2pr
             IF (liq2pr(i,k).NE.0.0_r8) THEN
                nliq2snow = liq2snow(i,k)*nliq2pr/liq2pr(i,k)   ! correction
             ELSE
                nliq2snow = liq2snow(i,k)
             ENDIF

             !           avoid roundoff problems generating negatives
             nliq2snow = nliq2snow*omsm
             nliq2pr = nliq2pr*omsm
             nice2pr = nice2pr*omsm

             !           final estimates of conversion to precip and snow
             prodprec(i,k) = (nliq2pr + nice2pr)
             prodsnow(i,k) = (nice2pr + nliq2snow)

             rcwn(i,l,k) =  cwat(i,k) + (cme(i,k)-   prodprec(i,k))*deltat
             rliq(i,l,k) = (cwat(i,k) + cme(i,k)*deltat)*(1.0_r8-fice(i,k)) - nliq2pr * deltat
             rice(i,l,k) = (cwat(i,k) + cme(i,k)*deltat)* fice(i,k)- nice2pr*deltat

             !           Save for sanity check later...  
             !           Putting sanity checks inside loops 100 and 800 screws up the 
             !           IBM compiler for reasons as yet unknown.  TBH
             cwnsave(i,l,k)      = cwn(i)
             cmesave(i,l,k)      = cme(i,k)
             prodprecsave(i,l,k) = prodprec(i,k)
             !           End of save for sanity check later...  

             !           final version of condensate to precip terms
             liq2pr(i,k) = nliq2pr
             liq2snow(i,k) = nliq2snow
             ice2pr(i,k) = nice2pr

             cwn(i) = rcwn(i,l,k)
             !
             ! update any remaining  provisional values
             !
             cwm(i) = (cwn(i) + cwat(i,k))*0.5_r8
             !
             ! update in cloud water
             !
             IF(cldm(i) > mincld) THEN
                icwc(i) = cwm(i)/cldm(i)
             ELSE
                icwc(i) = 0.0_r8
             ENDIF

          END DO              ! end of do i = 1,ncol
!RETURN
          !
          ! calculate provisional value of cloud water for
          ! evaporation of precipitate (evapprec) calculation
          !
          DO i = 1,ncol
             qtmp = qn(i,k) - cme(i,k)*deltat
             ttmp = tn(i,k) + deltat/cp * ( meltheat(i,k)       &
                  + (hlatv + hlatf*fice(i,k)) * cme(i,k) )
             esn = estblf(ttmp)
             qsn = MIN(epsqs*esn/(p(i,k) - omeps*esn),1.0_r8)
             qtl(i) = MAX((qsn - qtmp)/deltat,0.0_r8)
             relhum1(i) = qtmp/qsn
          END DO
          !
          DO i = 1,ncol

             !#ifdef PERGRO
             IF( landfrac(i) > 0.0_r8)THEN
                 !LAND       
                IF(PERGRO)THEN
                   evapprec(i,k) = conke_land*(1.0_r8 - MAX(cldm(i),mincld))* &
                        SQRT(precab(i))*(1.0_r8 - MIN(relhum1(i),1.0_r8))
                 ELSE
                   evapprec(i,k) = conke_land*(1.0_r8 - cldm(i))*SQRT(precab(i)) &
                        *(1.0_r8 - MIN(relhum1(i),1.0_r8))
                 ENDIF
             ELSE
                IF(PERGRO)THEN
                   evapprec(i,k) = conke*(1.0_r8 - MAX(cldm(i),mincld))* &
                        SQRT(precab(i))*(1.0_r8 - MIN(relhum1(i),1.0_r8))
                 ELSE
                   evapprec(i,k) = conke*(1.0_r8 - cldm(i))*SQRT(precab(i)) &
                        *(1.0_r8 - MIN(relhum1(i),1.0_r8))
                 ENDIF
             END IF
             !
             ! limit the evaporation to the amount which is entering the box
             ! or saturates the box
             !
             prtmp = precab(i)*gravit/pdel(i,k)
             evapprec(i,k) = MIN(evapprec(i,k), prtmp, qtl(i))*omsm
             !#ifdef PERGRO
             IF(PERGRO)THEN

                !           zeroing needed for pert growth
                evapprec(i,k) = 0.0_r8
             ENDIF

             !#endif
             !
             ! Partition evaporation of precipitate between rain and snow using
             ! the fraction of snow falling into the box. Determine the heating
             ! due to evaporation. Note that evaporation is positive (loss of precip,
             ! gain of vapor) and that heating is negative.
             IF (evapprec(i,k) > 0.0_r8) THEN
                evapsnow(i,k) = evapprec(i,k) * snowab(i) / precab(i)
                evapheat(i,k) = -hlatv * evapprec(i,k) - hlatf * evapsnow(i,k)
             ELSE 
                evapsnow(i,k) = 0.0_r8
                evapheat(i,k) = 0.0_r8
             END IF
             ! Account for the latent heat of fusion for liquid drops collected by falling snow
             prfzheat(i,k) = hlatf * liq2snow(i,k)
          END DO
!RETURN
          ! now remove the residual of any over-saturation. Normally,
          ! the oversaturated water vapor should have been removed by 
          ! cme formulation plus constraints by wet bulb tsp/qsp
          ! as computed above. However, because of non-linearity,
          ! addition of (cme-evapprec) to update t and q may still cause
          ! a very small amount of over saturation. It is called a
          ! residual of over-saturation because theoretically, cme
          ! should have taken care of all of large scale condensation.
          ! 

          DO i = 1,ncol
             qtmp = qn(i,k)-(cme(i,k)-evapprec(i,k))*deltat
             ttmp = tn(i,k) + deltat/cp * ( meltheat(i,k) + evapheat(i,k) + prfzheat(i,k)&
                  + (hlatv + hlatf*fice(i,k)) * cme(i,k) )
             esn = estblf(ttmp)
             qsn = MIN(epsqs*esn/(p(i,k) - omeps*esn),1.0_r8)
             !
             !Upper stratosphere and mesosphere, qsn calculated
             !above may be negative. Here just to skip it instead
             !of resetting it to 1 as in aqsat
             !
             IF(qtmp > qsn .AND. qsn > 0) THEN
                !calculate dqsdt, a more precise calculation
                !which taking into account different range of T 
                !can be found in aqsatd.F. Here follows
                !cond.F to calculate it.
                !
                denom = (p(i,k)-omeps*esn)*ttmp*ttmp
                dqsdt = clrh2o*qsn*p(i,k)/denom
                !
                !now extra condensation to bring air to just saturation
                !
                ctmp = (qtmp-qsn)/(1.0_r8+hlocp*dqsdt)/deltat
                cme(i,k) = cme(i,k)+ctmp
                !
                ! save residual on cmeres to addtion to cme on entering next iteration
                ! cme exit here contain the residual but overrided if back to iteration
                !
                cmeres(i) = ctmp
             ELSE
                cmeres(i) = 0.0_r8
             ENDIF
          END DO

       END DO ! end of do l = 1,iter
!RETURN
       !
       ! precipitation
       !
       DO i = 1,ncol
          precip(i) = precip(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
          precab(i) = precab(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
          IF(precab(i).LT.0.0_r8) precab(i)=0.0_r8
          !         snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodprec(i,k)*fice(i,k) - evapsnow(i,k))
          snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodsnow(i,k) - evapsnow(i,k))


!!$         if ((precab(i)) < 1.e-10_r8) then      
!!$            precab(i) = 0.
!!$            snowab(i) = 0.
!!$         endif
       END DO
    END DO ! level loop (k=1,pver)
    ! begin sanity checks
    error_found = .FALSE.
    DO k = k1mb,pver
       DO l = 1,iter
          DO i = 1,ncol
             IF (rcwn(i,l,k).LT.0.0_r8) error_found = .TRUE.
             IF (rliq(i,l,k).LT.0.0_r8) error_found = .TRUE.
             IF (rice(i,l,k).LT.0.0_r8) error_found = .TRUE.
          ENDDO
       ENDDO
    ENDDO
    IF (error_found) THEN
       DO k = k1mb,pver
          DO l = 1,iter
             DO i = 1,ncol
                IF (rcwn(i,l,k).LT.0.0_r8) THEN
                   WRITE (6,*) ' prob with neg rcwn1 ', rcwn(i,l,k),  &
                        cwnsave(i,l,k)
                   WRITE (6,*) ' cwat, cme*deltat, prodprec*deltat ', &
                        cwat(i,k), cmesave(i,l,k)*deltat,               &
                        prodprecsave(i,l,k)*deltat,                     &
                        (cmesave(i,l,k)-prodprecsave(i,l,k))*deltat
                   CALL endrun('PCOND')
                ENDIF
                IF (rliq(i,l,k).LT.0.0_r8) THEN
                   WRITE (6,*) ' prob with neg rliq1 ', rliq(i,l,k)
                   CALL endrun('PCOND')
                ENDIF
                IF (rice(i,l,k).LT.0.0_r8) THEN
                   WRITE (6,*) ' prob with neg rice ', rice(i,l,k)
                   CALL endrun('PCOND')
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    END IF
    ! end sanity checks

    RETURN
  END SUBROUTINE pcond





  !===============================================================================
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


  !##############################################################################

  SUBROUTINE findmcnew (ncol    ,pcols,pver, &
       k       ,precab  ,  t       ,p       , &
       cwm     ,cldm    ,cldmax  ,fice    ,coef    , &
       fwaut   ,fsaut   ,fracw   ,fsacw   ,fsaci   , &
       seaicef ,snowh   ,landm)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! calculate the conversion of condensate to precipitate
    ! 
    ! Method: 
    ! See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
    !  model climate using diagnosed and 
    !  predicted condensate parameterizations, 1998, J. Clim., 11,
    !  pp1587---1614.
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------
    !   use phys_grid, only: get_rlat_all_p
    !
    ! input args
    !
    INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: k                     ! level index

    REAL(r8), INTENT(in) :: precab(pcols)        ! rate of precipitation from above (kg / (m**2 * s))
    REAL(r8), INTENT(in) :: t(pcols,pver)        ! temperature       (K)
    REAL(r8), INTENT(in) :: p(pcols,pver)        ! pressure          (Pa)
    REAL(r8), INTENT(in) :: cldm(pcols)          ! cloud fraction
    REAL(r8), INTENT(in) :: cldmax(pcols)        ! max cloud fraction above this level
    REAL(r8), INTENT(in) :: cwm(pcols)           ! condensate mixing ratio (kg/kg)
    REAL(r8), INTENT(in) :: fice(pcols)          ! fraction of cwat that is ice
    REAL(r8), INTENT(in) :: seaicef(pcols)       ! sea ice fraction 
    !REAL(r8), INTENT(in) :: snowab(pcols)        ! rate of snow from above (kg / (m**2 * s))
    REAL(r8), INTENT(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    REAL(r8), INTENT(in) :: landm(pcols)    
    ! output arguments
    REAL(r8), INTENT(out) :: coef(pcols)          ! conversion rate (1/s)
    REAL(r8), INTENT(out) :: fwaut(pcols)         ! relative importance of liquid autoconversion (a diagnostic)
    REAL(r8), INTENT(out) :: fsaut(pcols)         ! relative importance of ice autoconversion    (a diagnostic)
    REAL(r8), INTENT(out) :: fracw(pcols)         ! relative importance of rain accreting liquid (a diagnostic)
    REAL(r8), INTENT(out) :: fsacw(pcols)         ! relative importance of snow accreting liquid (a diagnostic)
    REAL(r8), INTENT(out) :: fsaci(pcols)         ! relative importance of snow accreting ice    (a diagnostic)

    ! work variables

    INTEGER i
    INTEGER ii
    INTEGER ind(pcols)
    INTEGER ncols

    REAL(r8), PARAMETER :: degrad = 57.296_r8 ! divide by this to convert degrees to radians
    REAL(r8) alpha                ! ratio of 3rd moment radius to 2nd
    REAL(r8) capc                 ! constant for autoconversion
    REAL(r8) capn                 ! local cloud particles / cm3
    !REAL(r8) capnoice             ! local cloud particles when not over sea ice / cm3
    REAL(r8) capnsi               ! sea ice cloud particles / cm3
    REAL(r8) capnc                ! cold and oceanic cloud particles / cm3
    REAL(r8) capnw                ! warm continental cloud particles / cm3
    REAL(r8) ciaut                ! coefficient of autoconversion of ice (1/s)
    REAL(r8) ciautb               ! coefficient of autoconversion of ice (1/s)
    REAL(r8) cldloc(pcols)        ! non-zero amount of cloud
    REAL(r8) cldpr(pcols)         ! assumed cloudy volume occupied by rain and cloud
    REAL(r8) con1                 ! work constant
    REAL(r8) con2                 ! work constant
    REAL(r8) convfw               ! constant used for fall velocity calculation
    REAL(r8) cracw                ! constant used for rain accreting water
    REAL(r8) critpr               ! critical precip rate collection efficiency changes
    REAL(r8) csacx                ! constant used for snow accreting liquid or ice
!!$   real(r8) dtice                ! interval for transition from liquid to ice
    REAL(r8) effc                 ! collection efficiency
    REAL(r8) icemr(pcols)         ! in-cloud ice mixing ratio
    REAL(r8) icrit                ! threshold for autoconversion of ice
    REAL(r8) kconst               ! const for terminal velocity (stokes regime)
    REAL(r8) liqmr(pcols)         ! in-cloud liquid water mixing ratio
    REAL(r8) pracw                ! rate of rain accreting water
    REAL(r8) prlloc(pcols)        ! local rain flux in mm/day
    REAL(r8) prscgs(pcols)        ! local snow amount in cgs units
    REAL(r8) psaci                ! rate of collection of ice by snow (lin et al 1983)
    REAL(r8) psacw                ! rate of collection of liquid by snow (lin et al 1983)
    REAL(r8) psaut                ! rate of autoconversion of ice condensate
    REAL(r8) ptot                 ! total rate of conversion
    REAL(r8) pwaut                ! rate of autoconversion of liquid condensate
    REAL(r8) r3l                  ! volume radius
    REAL(r8) r3lcrit              ! critical radius at which autoconversion become efficient
    REAL(r8) rainmr(pcols)        ! in-cloud rain mixing ratio
    REAL(r8) rat1                 ! work constant
    REAL(r8) rat2                 ! work constant
!!$   real(r8) rdtice               ! recipricol of dtice
    REAL(r8) rho(pcols)           ! density (mks units)
    REAL(r8) rhocgs               ! density (cgs units)
    !REAL(r8) rlat(pcols)          ! latitude (radians)
    REAL(r8) snowfr               ! fraction of precipate existing as snow
    REAL(r8) totmr(pcols)         ! in-cloud total condensate mixing ratio
    REAL(r8) vfallw               ! fall speed of precipitate as liquid
    !REAL(r8) wp                   ! weight factor used in calculating pressure dep of autoconversion
    !REAL(r8) wsi                  ! weight factor for sea ice
    REAL(r8) wt                   ! fraction of ice
    !REAL(r8) wland                ! fraction of land

    !      real(r8) csaci
    !      real(r8) csacw
    !      real(r8) cwaut
    !      real(r8) efact
    !      real(r8) lamdas
    !      real(r8) lcrit
    !      real(r8) rcwm
    !      real(r8) r3lc2
    !      real(r8) snowmr(pcols)
    !      real(r8) vfalls

    REAL(8) ftot

    !     inline statement functions
    REAL(r8) heavy, heavym, a1, a2, heavyp, heavymp
    heavy(a1,a2) = MAX(0.0_r8,SIGN(1.0_r8,a1-a2))  ! heavyside function
    heavym(a1,a2) = MAX(0.01_r8,SIGN(1.0_r8,a1-a2))  ! modified heavyside function
    !
    ! New heavyside functions to perhaps address error growth problems
    !
    heavyp(a1,a2) = a1/(a2+a1+1.0e-36_r8)
    heavymp(a1,a2) = (a1+0.01_r8*a2)/(a2+a1+1.0e-36_r8)

    ! critical precip rate at which we assume the collector drops can change the
    ! drop size enough to enhance the auto-conversion process (mm/day)
    critpr = 0.5_r8

    convfw = 1.94_r8*2.13_r8*SQRT(rhow*1000.0_r8*9.81_r8*2.7e-4_r8)

    ! liquid microphysics
    !      cracw = 6                 ! beheng
    cracw = 0.884_r8*SQRT(9.81_r8/(rhow*1000.0_r8*2.7e-4_r8)) ! tripoli and cotton

    ! ice microphysics
    !      ciautb = 6.e-4_r8
    !      ciautb = 1.e-3_r8
    ciautb = 5.e-4_r8
    !      icritw = 1.e-5_r8
    !      icritw = 5.e-5_r8
!!$   icritw = 4.e-4_r8
    !      icritc = 4.e-6_r8
    !      icritc = 6.e-6_r8
!!$   icritc = 5.e-6_r8

!!$   dtice = 20.0_r8
!!$   rdtice = 1.0_r8/dtice

    capnw = 400.0_r8              ! warm continental cloud particles / cm3
    capnc = 150.0_r8              ! cold and oceanic cloud particles / cm3
    !  capnsi = 40.0_r8              ! sea ice cloud particles density  / cm3
    capnsi = 75.0_r8              ! sea ice cloud particles density  / cm3

    kconst = 1.18e6_r8           ! const for terminal velocity

    !      effc = 1.0_r8                 ! autoconv collection efficiency following boucher 96
    !      effc = 0.55_r8*0.05_r8           ! autoconv collection efficiency following baker 93
    effc = 0.55_r8                ! autoconv collection efficiency following tripoli and cotton
    !   effc = 0.0_r8    ! turn off warm-cloud autoconv
    alpha = 1.1_r8**4
    capc = pi**(-0.333_r8)*kconst*effc *(0.75_r8)**(1.333_r8)*alpha  ! constant for autoconversion

    r3lcrit = 10.0e-6_r8         ! 10.0u  crit radius where liq conversion begins
    !
    ! find all the points where we need to do the microphysics
    ! and set the output variables to zero
    !
    ncols = 0
    DO i = 1,ncol
       coef(i) = 0.0_r8
       fwaut(i) = 0.0_r8
       fsaut(i) = 0.0_r8
       fracw(i) = 0.0_r8
       fsacw(i) = 0.0_r8
       fsaci(i) = 0.0_r8
       liqmr(i) = 0.0_r8
       rainmr(i) = 0.0_r8
       IF (cwm(i) > 1.e-20_r8) THEN
          ncols = ncols + 1
          ind(ncols) = i
       ENDIF
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cdir nodep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DIR$ CONCURRENT
    DO ii = 1,ncols
       i = ind(ii)
       !
       ! the local cloudiness at this level
       !
       cldloc(i) = MAX(cldmin,cldm(i))
       !
       ! a weighted mean between max cloudiness above, and this layer
       !
       cldpr(i) = MAX(cldmin,(cldmax(i)+cldm(i))*0.5_r8)
       !
       ! decompose the suspended condensate into
       ! an incloud liquid and ice phase component
       !
       totmr(i) = cwm(i)/cldloc(i)
       icemr(i) = totmr(i)*fice(i)
       liqmr(i) = totmr(i)*(1.0_r8-fice(i))
       !
       ! density
       !
       rho(i) = p(i,k)/(287.0_r8*t(i,k))
       rhocgs = rho(i)*1.e-3_r8     ! density in cgs units
       !
       ! decompose the precipitate into a liquid and ice phase
       !
       IF (t(i,k) > t0) THEN
          vfallw = convfw/SQRT(rho(i))
          rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
          snowfr = 0
          !        snowmr(i)
       ELSE
          snowfr = 1
          rainmr(i) = 0.0_r8
       ENDIF
       !     rainmr(i) = (precab(i)-snowab(i))/(rho(i)*vfallw*cldpr(i))
       !
       ! local snow amount in cgs units
       !
       prscgs(i) = precab(i)/cldpr(i)*0.1_r8*snowfr
       !     prscgs(i) = snowab(i)/cldpr(i)*0.1_r8
       !
       ! local rain amount in mm/day
       !
       prlloc(i) = precab(i)*86400.0_r8/cldpr(i)
    END DO

    con1 = 1.0_r8/(1.333_r8*pi)**0.333_r8 * 0.01_r8 ! meters
    !
    ! calculate the conversion terms
    !
    !   call get_rlat_all_p(lchnk, ncol, rlat)

    !!!!!!!!!!!!!!!!!!!cdir nodep
    !!!!!!!!!!!!!!!!!!!!DIR$ CONCURRENT
    DO ii = 1,ncols
       i = ind(ii)
       rhocgs = rho(i)*1.e-3_r8     ! density in cgs units
       !
       ! exponential temperature factor
       !
       !        efact = exp(0.025_r8*(t(i,k)-t0))
       !
       ! some temperature dependent constants
       !
!!$      wt = min(1.0_r8,max(0.0_r8,(t0-t(i,k))*rdtice))
       wt = fice(i)
       icrit = icritc*wt + icritw*(1-wt)
       !
       ! jrm Reworked droplet number concentration algorithm
       ! Start with pressure-dependent value appropriate for continental air
       ! Note: reltab has a temperature dependence here
       capn = capnw + (capnc-capnw) * MIN(1.0_r8,MAX(0.0_r8,1.0_r8-(p(i,k)-0.8_r8*p(i,pver))/&
             (0.2_r8*p(i,pver))))
       ! Modify for snow depth over land
       capn = capn + (capnc-capn) * MIN(1.0_r8,MAX(0.0_r8,snowh(i)*10.0_r8))
       ! Ramp between polluted value over land to clean value over ocean.
       capn = capn + (capnc-capn) * MIN(1.0_r8,MAX(0.0_r8,1.0_r8-landm(i)))
       ! Ramp between the resultant value and a sea ice value in the presence of ice.
       capn = capn + (capnsi-capn) * MIN(1.0_r8,MAX(0.0_r8,seaicef(i)))
       ! end jrm
       !      
       !#ifdef DEBUG2
       !      if ( (lat(i) == latlook(1)) .or. (lat(i) == latlook(2)) ) then
       !         if (i == ilook(1)) then
       !            write (6,*) ' findmcnew: lat, k, seaicef, landm, wp, capnoice, capn ', &
       !                 lat(i), k, seaicef(i), landm(i,lat(i)), wp, capnoice, capn
       !         endif
       !      endif
       !#endif

       !
       ! useful terms in following calculations
       !
       rat1 = rhocgs/rhow
       rat2 = liqmr(i)/capn
       con2 = (rat1*rat2)**0.333_r8
       !
       ! volume radius
       !
       !        r3l = (rhocgs*liqmr(i)/(1.333*pi*capn*rhow))**0.333 * 0.01 ! meters
       r3l = con1*con2
       !
       ! critical threshold for autoconversion if modified for mixed phase
       ! clouds to mimic a bergeron findeisen process
       ! r3lc2 = r3lcrit*(1.-0.5*fice(i)*(1-fice(i)))
       !
       ! autoconversion of liquid
       !
       !        cwaut = 2.e-4
       !        cwaut = 1.e-3
       !        lcrit = 2.e-4
       !        lcrit = 5.e-4
       !        pwaut = max(0._r8,liqmr(i)-lcrit)*cwaut
       !
       ! pwaut is following tripoli and cotton (and many others)
       ! we reduce the autoconversion below critpr, because these are regions where
       ! the drop size distribution is likely to imply much smaller collector drops than
       ! those relevant for a cloud distribution corresponding to the value of effc = 0.55
       ! suggested by cotton (see austin 1995 JAS, baker 1993)

       ! easy to follow form
       !        pwaut = capc*liqmr(i)**2*rhocgs/rhow
       !    $           *(liqmr(i)*rhocgs/(rhow*capn))**(.333)
       !    $           *heavy(r3l,r3lcrit)
       !    $           *max(0.10_r8,min(1._r8,prlloc(i)/critpr))
       ! somewhat faster form
       !#define HEAVYNEW
       IF(HEAVYNEW)THEN
          pwaut = capc*liqmr(i)**2*rat1*con2*heavymp(r3l,r3lcrit) * &
               MAX(0.10_r8,MIN(1.0_r8,prlloc(i)/critpr))
       ELSE
          pwaut = capc*liqmr(i)**2*rat1*con2*heavym(r3l,r3lcrit)* &
               MAX(0.10_r8,MIN(1.0_r8,prlloc(i)/critpr))
       ENDIF
       !#endif
       !
       ! autoconversion of ice
       !
       !        ciaut = ciautb*efact
       ciaut = ciautb
       !        psaut = capc*totmr(i)**2*rhocgs/rhoi
       !     $           *(totmr(i)*rhocgs/(rhoi*capn))**(.333)
       !
       ! autoconversion of ice condensate
       !
       !#ifdef PERGRO
       IF(PERGRO)THEN
          psaut = heavyp(icemr(i),icrit)*icemr(i)*ciaut
       ELSE 
          !#else
          psaut = MAX(0.0_r8,icemr(i)-icrit)*ciaut
       END IF
       !#endif
       !
       ! collection of liquid by rain
       !
       !        pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
       pracw = cracw*rho(i)*SQRT(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)
       !!      pracw = 0.
       !
       ! the following lines calculate the slope parameter and snow mixing ratio
       ! from the precip rate using the equations found in lin et al 83
       ! in the most natural form, but it is expensive, so after some tedious
       ! algebraic manipulation you can use the cheaper form found below
       !            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
       !     $               *0.01   ! convert from cm/s to m/s
       !            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
       !            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
       !            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25
       !            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd)
       !
       ! coefficient for collection by snow independent of phase
       !
       csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05

       !
       ! collection of liquid by snow (lin et al 1983)
       !
       psacw = csacx*liqmr(i)*esw
       !#ifdef PERGRO
       IF(PERGRO)THEN
          ! this is necessary for pergro
          psacw = 0.0_r8
       END IF
       ! #endif
       !
       ! collection of ice by snow (lin et al 1983)
       !
       psaci = csacx*icemr(i)*esi
       !
       ! total conversion of condensate to precipitate
       !
       ptot = pwaut + psaut + pracw + psacw + psaci
       !
       ! the recipricol of cloud water amnt (or zero if no cloud water)
       !
       !         rcwm =  totmr(i)/(max(totmr(i),small)**2)
       !
       ! turn the tendency back into a loss rate (1/seconds)
       !
       IF (totmr(i) > 0.0_r8) THEN
          coef(i) = ptot/totmr(i)
       ELSE
          coef(i) = 0.0_r8
       ENDIF

       IF (ptot.GT.0.0_r8) THEN
          fwaut(i) = pwaut/ptot
          fsaut(i) = psaut/ptot
          fracw(i) = pracw/ptot
          fsacw(i) = psacw/ptot
          fsaci(i) = psaci/ptot
       ELSE
          fwaut(i) = 0.0_r8
          fsaut(i) = 0.0_r8
          fracw(i) = 0.0_r8
          fsacw(i) = 0.0_r8
          fsaci(i) = 0.0_r8
       ENDIF

       ftot = fwaut(i)+fsaut(i)+fracw(i)+fsacw(i)+fsaci(i)
       !      if (abs(ftot-1._r8).gt.1.e-14_r8.and.ftot.ne.0._r8) then
       !         write (6,*) ' something is wrong in findmcnew ', ftot, &
       !              fwaut(i),fsaut(i),fracw(i),fsacw(i),fsaci(i)
       !         write (6,*) ' unscaled ', ptot, &
       !              pwaut,psaut,pracw,psacw,psaci
       !         write (6,*) ' totmr, liqmr, icemr ', totmr(i), liqmr(i), icemr(i)
       !         call endrun()
       !      endif
    END DO
    !#ifdef DEBUG
    !   i = icollook(nlook)
    !   if (lchnk == lchnklook(nlook) ) then
    !      write (6,*)
    !      write (6,*) '------', k, i, lchnk
    !      write (6,*) ' liqmr, rainmr,precab ', liqmr(i), rainmr(i), precab(i)*8.64e4
    !      write (6,*) ' frac: waut,saut,racw,sacw,saci ', &
    !           fwaut(i), fsaut(i), fracw(i), fsacw(i), fsaci(i)
    !   endif
    !#endif

    RETURN
  END SUBROUTINE findmcnew

  !##############################################################################

  SUBROUTINE findsp (ncol,pcols,pver, q, t, p, tsp, qsp)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !     find the wet bulb temperature for a given t and q
    !     in a longitude height section
    !     wet bulb temp is the temperature and spec humidity that is 
    !     just saturated and has the same enthalpy
    !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
    !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
    !
    ! Method: 
    ! a Newton method is used
    ! first guess uses an algorithm provided by John Petch from the UKMO
    ! we exclude points where the physical situation is unrealistic
    ! e.g. where the temperature is outside the range of validity for the
    !      saturation vapor pressure, or where the water vapor pressure
    !      exceeds the ambient pressure, or the saturation specific humidity is 
    !      unrealistic
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------
    !
    !     input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: q(pcols,pver)        ! water vapor (kg/kg)
    REAL(r8), INTENT(in) :: t(pcols,pver)        ! temperature (K)
    REAL(r8), INTENT(in) :: p(pcols,pver)        ! pressure    (Pa)
    !
    ! output arguments
    !
    REAL(r8), INTENT(out) :: tsp(pcols,pver)      ! saturation temp (K)
    REAL(r8), INTENT(out) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
    !
    ! local variables
    !
    INTEGER i                 ! work variable
    INTEGER k                 ! work variable
    LOGICAL lflg              ! work variable
    INTEGER iter              ! work variable
    INTEGER l                 ! work variable
    LOGICAL :: error_found

    REAL(r8) omeps                ! 1 minus epsilon
    REAL(r8) trinv                ! work variable
    REAL(r8) es                   ! sat. vapor pressure
    REAL(r8) desdt                ! change in sat vap pressure wrt temperature
    !     real(r8) desdp                ! change in sat vap pressure wrt pressure
    REAL(r8) dqsdt                ! change in sat spec. hum. wrt temperature
    REAL(r8) dgdt                 ! work variable
    REAL(r8) g                    ! work variable
    REAL(r8) weight(pcols)        ! work variable
    REAL(r8) hlatsb               ! (sublimation)
    REAL(r8) hlatvp               ! (vaporization)
    REAL(r8) hltalt(pcols,pver)   ! lat. heat. of vap.
    REAL(r8) tterm                ! work var.
    REAL(r8) qs                   ! spec. hum. of water vapor
    REAL(r8) tc                   ! crit temp of transition to ice

    ! work variables
    REAL(r8) t1, q1, dt, dq
    REAL(r8) dtm, dqm
    REAL(r8) qvd, a1, tmp
    REAL(r8) rair
    REAL(r8) r1b, c1, c2, c3
    REAL(r8) denom
    REAL(r8) dttol
    REAL(r8) dqtol
    INTEGER  doit(pcols) 
    REAL(r8) enin(pcols)
    REAL(r8) enout(pcols)
    REAL(r8) tlim(pcols)
    REAL(r8) epsqs

    tsp= 0.0_r8
    qsp= 0.0_r8
    tterm= 0.0_r8
    t1= 0.0_r8; q1= 0.0_r8; dt= 0.0_r8; dq= 0.0_r8
    qvd= 0.0_r8; a1= 0.0_r8; tmp= 0.0_r8;
    qs   = 0.0_r8
    tc   = 0.0_r8
    weight= 0.0_r8
    hltalt= 0.0_r8
    doit=0
    hlatsb= 0.0_r8
    hlatvp= 0.0_r8
    enin= 0.0_r8
    enout= 0.0_r8
    tlim= 0.0_r8
    es   = 0.0_r8
    g    = 0.0_r8
    desdt= 0.0_r8
    dqsdt= 0.0_r8
    dgdt = 0.0_r8
    epsqs =  epsilo

    omeps = 1.0_r8 - epsqs
    trinv = 1.0_r8/ttrice
    a1 = 7.5_r8*LOG(10.0_r8)
    rair =  287.04_r8
    c3 = rair*a1/cp
    dtm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dqm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
    dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
    !  tmin = 173.16_r8 ! the coldest temperature we can deal with
    !
    ! max number of times to iterate the calculation
    iter = iterp
    !
    DO k = k1mb,pver

       !
       ! first guess on the wet bulb temperature
       !
       DO i = 1,ncol

          !#ifdef DEBUG
          !         if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
          !            write (6,*) ' '
          !            write (6,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
          !         endif
          !#endif
          ! limit the temperature range to that relevant to the sat vap pres tables
          tlim(i) = MIN(MAX(t(i,k),173.0_r8),373.0_r8)
          es      = estblf(tlim(i))
          denom   = p(i,k) - omeps*es
          qs      = epsqs*es/denom
          doit(i) = 0
          enout(i)= 1.0_r8
          ! make sure a meaningful calculation is possible
          IF (p(i,k) > 5.0_r8*es .AND. qs > 0.0_r8 .AND. qs < 0.5_r8) THEN
             !
             ! Saturation specific humidity
             !
             qs = MIN(epsqs*es/denom,1.0_r8)
             !
             ! "generalized" analytic expression for t derivative of es
             !  accurate to within 1 percent for 173.16 < t < 373.16
             !
             ! Weighting of hlat accounts for transition from water to ice
             ! polynomial expression approximates difference between es over
             ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
             ! -40): required for accurate estimate of es derivative in transition
             ! range from ice to water also accounting for change of hlatv with t
             ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
             !
             tc        = tlim(i) - t0
             lflg      = (tc >= -ttrice .AND. tc < 0.0_r8)
             weight(i) = MIN(-tc*trinv,1.0_r8)
             hlatsb    = hlatv + weight(i)*hlatf
             hlatvp    = hlatv - 2369.0_r8*tc
             IF (tlim(i) < t0) THEN
                hltalt(i,k) = hlatsb
             ELSE
                hltalt(i,k) = hlatvp
             END IF
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

             ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp   = q(i,k) - qs
             c1    = hltalt(i,k)*c3
             c2    = (tlim(i) + 36.0_r8)**2
             r1b   = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
             !#ifdef DEBUG
             !             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
             !                write (6,*) ' relative humidity ', q(i,k)/qs
             !                write (6,*) ' first guess ', tsp(i,k)
             !             endif
             !#endif
             es        = estblf(tsp(i,k))
             qsp(i,k)  = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
          ELSE
             doit(i)   = 1
             tsp (i,k) = tlim(i)
             qsp (i,k) = q   (i,k)
             enin(i)   = 1.0_r8
          ENDIF
       END DO   ! end do i
       !
       ! now iterate on first guess
       !
       DO l = 1, iter
          dtm = 0.0_r8
          dqm = 0.0_r8
          DO i = 1,ncol
             IF (doit(i) == 0) THEN
                es = estblf(tsp(i,k))
                !
                ! Saturation specific humidity
                !
                qs = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                !
                ! "generalized" analytic expression for t derivative of es
                ! accurate to within 1 percent for 173.16 < t < 373.16
                !
                ! Weighting of hlat accounts for transition from water to ice
                ! polynomial expression approximates difference between es over
                ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                ! -40): required for accurate estimate of es derivative in transition
                ! range from ice to water also accounting for change of hlatv with t
                ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                !
                tc        = tsp(i,k) - t0
                lflg      = (tc >= -ttrice .AND. tc < 0.0_r8)
                weight(i) = MIN(-tc*trinv,1.0_r8)
                hlatsb    = hlatv + weight(i)*hlatf
                hlatvp    = hlatv - 2369.0_r8*tc
                IF (tsp(i,k) < t0) THEN
                   hltalt(i,k) = hlatsb
                ELSE
                   hltalt(i,k) = hlatvp
                END IF
                IF (lflg) THEN
                   tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
                ELSE
                   tterm = 0.0_r8
                END IF
                desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
                !g    = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
                g     = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
                dgdt  = -(cp + hltalt(i,k)*dqsdt)
                t1    = tsp(i,k) - g/dgdt
                dt    = ABS(t1 - tsp(i,k))/t1
                tsp(i,k) = MAX(t1,tmin)
                es = estblf(tsp(i,k))
                q1 = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                q1 = MAX(q1,1.e-12_r8)
                dq = ABS(q1 - qsp(i,k))/MAX(q1,1.e-12_r8)
                qsp(i,k) = q1
                !#ifdef DEBUG
                !               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                !                  write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
                !               endif
                !#endif
                dtm = MAX(dtm,dt)
                dqm = MAX(dqm,dq)
                ! if converged at this point, exclude it from more iterations
                IF (dt < dttol .AND. dq < dqtol) THEN
                   doit(i) = 2
                ENDIF
                enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
                ! bail out if we are too near the end of temp range
                IF (tsp(i,k) < 174.16_r8) THEN
                   doit(i) = 4
                ENDIF
             ELSE
             ENDIF
          END DO              ! do i = 1,ncol

          IF (dtm < dttol .AND. dqm < dqtol) THEN
             go to 10
          ENDIF

       END DO                 ! do l = 1,iter
10     CONTINUE

       error_found = .FALSE.
       IF (dtm > dttol .OR. dqm > dqtol) THEN
          DO i = 1,ncol
             IF (doit(i) == 0) error_found = .TRUE.
          END DO
          IF (error_found) THEN
             DO i = 1,ncol
                IF (doit(i) == 0) THEN
                   WRITE (6,*) ' findsp not converging at point i, k ', i, k
                   WRITE (6,*) ' t, q, p, enin   ', t(i,k)  , q(i,k)  , p(i,k), enin(i)
                   WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                   CALL endrun ('FINDSP')
                ENDIF
             END DO
          ENDIF
       ENDIF
       DO i = 1,ncol
          IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
             error_found = .TRUE.
          ENDIF
       END DO
       IF (error_found) THEN
          DO i = 1,ncol
             IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
                WRITE (6,*) ' the enthalpy is not conserved for point ', &
                     i, k, enin(i), enout(i)
                WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                CALL endrun ('FINDSP')
             ENDIF
          END DO
       ENDIF

    END DO                    ! level loop (k=1,pver)

    RETURN
  END SUBROUTINE findsp




  SUBROUTINE inimc( pver, plat,tmeltx, rhonotx, rh2ox ,hypi)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! initialize constants for the prognostic condensate
    ! 
    ! Author: P. Rasch, April 1997
    ! 
    !-----------------------------------------------------------------------
    !use phys_grid, only: get_chunk_coord_p
    !use pmgrid, only: plev, plevp
    !use dycore, only: dycore_is, get_resolution
    INTEGER, INTENT(IN   ) :: pver
    INTEGER, INTENT(IN   ) :: plat
    REAL(r8), INTENT(in) :: tmeltx
    REAL(r8), INTENT(in) :: rhonotx
    !REAL(r8), INTENT(in) :: gravitx
    REAL(r8), INTENT(in) :: rh2ox
    REAL(KIND=r8), INTENT(IN) :: hypi(pver+1)

    !CHARACTER(len=7) :: get_resolution
    !external get_resolution
    integer k

!!!!!!!#include <comhyb.h>

    !REAL(r8) signgam              ! variable required by cray gamma function
    !external gamma
    rhonot = rhonotx          ! air density at surface (gm/cm3)
    !gravit = gravitx
    rh2o   = rh2ox!! Gas constant for water vapor
    rhos = 0.1_r8                 ! assumed snow density (gm/cm3)
    rhow = 1.0_r8                 ! water density
    rhoi = 1.0_r8                 ! ice density
    esi = 1.0_r8                 ! collection efficient for ice by snow
    esw = 0.1_r8                 ! collection efficient for water by snow
    t0 = tmeltx               ! approximate freezing temp
    cldmin = 0.02_r8             ! assumed minimum cloud amount
    small = 1.e-22_r8            ! a small number compared to unity
    c = 152.93_r8                ! constant for graupel like snow cm**(1-d)/s
    d = 0.25_r8                  ! constant for graupel like snow
    nos = 3.0e-2_r8               ! particles snow / cm**4
    pi = 4.0_r8*ATAN(1.0_r8)
    prhonos = pi*rhos*nos
    thrpd = 3.0_r8 + d
    IF (d == 0.25_r8) THEN
       gam3pd = 2.549256966718531_r8 ! only right for d = 0.25
       gam4pd = 8.285085141835282_r8
    ELSE
       IF( UNICOSMP)THEN
          !    CALL gamma(3.0_r8+d, signgam, gam3pd)
          !    gam3pd = SIGN(EXP(gam3pd),signgam)
          !    CALL gamma(4.0_r8+d, signgam, gam4pd)
          !    gam4pd = SIGN(EXP(gam4pd),signgam)
          WRITE (6,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
       ELSE
          WRITE (6,*) ' can only use d ne 0.25 on a cray '
          STOP
       ENDIF
    ENDIF

    mcon01 = pi*nos*c*gam3pd/4.0_r8
    mcon02 = 1.0_r8/(c*gam4pd*SQRT(rhonot)/(6*prhonos**(d/4.0_r8)))
    mcon03 = -(0.5_r8+d/4.0_r8)
    mcon04 = 4.0_r8/(4.0_r8+d)
    mcon05 = (3+d)/(4+d)
    mcon06 = (3+d)/4.0_r8
    mcon07 = mcon01*SQRT(rhonot)*mcon02**mcon05/prhonos**mcon06
    mcon08 = -0.5_r8/(4.0_r8+d)

    !  find the level about 1mb, we wont do the microphysics above this level 0.001
 
    k1mb = 1
    do k=1,pver-1
      if (hypi(k) < 0.001_r8 .and. hypi(k+1) >= 0.001_r8) then
         if (0.001_r8 - hypi(k) < hypi(k+1) - 0.001_r8) then
            k1mb = k
         else
            k1mb = k + 1
         end if
    !     goto 20
      end if
    end do
    !if (masterproc) then
    !   write(6,*)'inimc: model levels bracketing 1 mb not found'
    !end if
    !  call endrun
    !k1mb = 1
    !20 if( masterproc ) write(6,*)'inimc: model level nearest 1 mb is',k1mb,'which is',hypm(k1mb),'pascals'

    !#ifdef DEBUG
    !!
    !! Set indicies of the point to examine for debugging
    !!
    !   latlook(:) = (/64, 32/)   ! Latitude indices to examine
    !   ilook(:)   = (/1,   1/)   ! Longitude indicex to examine
    !   call get_chunk_coord_p( nlook, ilook, latlook, icollook, lchnklook )
    !#endif
    !   if( masterproc ) write (6,*) 'cloud water initialization by inimc complete '
!!$   real(r8),public,parameter::  conke  = 1.e-6    ! tunable constant for evaporation of precip
!!$   real(r8),public,parameter::  conke  =  2.e-6    ! tunable constant for evaporation of precip
!!$   conke = 1.e-5
!!$   conke = 1.e-6
    !      icritw = 1.e-5_r8
    !      icritw = 5.e-5_r8
!!$   icritw = 4.e-4_r8
    !      icritc = 4.e-6_r8
    !      icritc = 6.e-6_r8
!!$   icritc = 5.e-6_r8

    IF ( dycore_is ('LR') ) THEN
       icritw = 2.e-4_r8 !! threshold for autoconversion of warm ice
       icritc = 9.5e-6_r8    ! threshold for autoconversion of cold ice
       conke  =  5.e-6_r8    ! tunable constant for evaporation of precip
    ELSE
       IF( get_resolution(plat) == 'T85' ) THEN
          icritw = 4.e-4_r8! threshold for autoconversion of warm ice
          icritc = 16.0e-6_r8    ! threshold for autoconversion of cold ice
          conke  = 5.e-6_r8    ! tunable constant for evaporation of precip
       ELSEIF( get_resolution(plat) == 'T62' ) THEN
          icritw = 4.e-4_r8! threshold for autoconversion of warm ice
          icritc =  5.e-6_r8    ! threshold for autoconversion of cold ice
          conke  = 10.e-6_r8    ! tunable constant for evaporation of precip
       ELSEIF( get_resolution(plat) == 'T31' ) THEN
          icritw = 4.e-4_r8! threshold for autoconversion of warm ice
          icritc =  3.e-6_r8    ! threshold for autoconversion of cold ice
          conke  = 10.e-6_r8    ! tunable constant for evaporation of precip
       ELSE
          icritw = 4.e-4_r8! threshold for autoconversion of warm ice
          icritc =  5.e-6_r8    ! threshold for autoconversion of cold ice
          conke  = 10.e-6_r8    ! tunable constant for evaporation of precip
       ENDIF
    ENDIF
    icritw =   8.e-4_r8    ! threshold for autoconversion of warm ice
    icritc =  10.e-6_r8    ! threshold for autoconversion of cold ice
    conke  = 180.e-5_r8    ! tunable constant for evaporation of precip
    conke_land  = 180.e-5_r8!50.e-6_r8    ! tunable constant for evaporation of precip

!    icritw = 4.e-4_r8! threshold for autoconversion of warm ice
!    icritc = 16.0e-6_r8    ! threshold for autoconversion of cold ice
!    conke  = 50.e-6_r8    ! tunable constant for evaporation of precip

    RETURN
  END SUBROUTINE inimc
  !END MODULE CLDWAT

  !MODULE PKG_CLD_SEDIMENT

  !===============================================================================
  SUBROUTINE cld_sediment_vel (ncol, pcols,  pver, pverp,      &
       icefrac , landfrac , pmid    , t       , &
       cloud   , cldliq  , cldice  , pvliq   , pvice   , landm, snowh)

    !----------------------------------------------------------------------

    ! Compute gravitational sedimentation velocities for cloud liquid water
    ! and ice, based on Lawrence and Crutzen (1998).

    ! LIQUID

    ! The fall velocities assume that droplets have a gamma distribution
    ! with effective radii for land and ocean as assessed by Han et al.;
    ! see Lawrence and Crutzen (1998) for a derivation.

    ! ICE

    ! The fall velocities are based on data from McFarquhar and Heymsfield
    ! or on Stokes terminal velocity for spheres and the effective radius.

    ! NEED TO BE CAREFUL - VELOCITIES SHOULD BE AT THE *LOWER* INTERFACE
    ! (THAT IS, FOR K+1), FLUXES ARE ALSO AT THE LOWER INTERFACE (K+1), 
    ! BUT MIXING RATIOS ARE AT THE MIDPOINTS (K)...

    ! NOTE THAT PVEL IS ON PVERP (INTERFACES), WHEREAS VFALL IS FOR THE CELL
    ! AVERAGES (I.E., MIDPOINTS); ASSUME THE FALL VELOCITY APPLICABLE TO THE 
    ! LOWER INTERFACE (K+1) IS THE SAME AS THAT APPLICABLE FOR THE CELL (V(K))

    !-----------------------------------------------------------------------
    !     MATCH-MPIC version 2.0, Author: mgl, March 1998
    ! adapted by P. J. Rasch
    !            B. A. Boville, September 19, 2002
    !            P. J. Rasch    May 22, 2003 (added stokes flow calc for liquid
    !                                         drops based on effect radii)
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(in) :: ncol                     ! number of colums to process
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                  ! number of vertical levels + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in)  :: icefrac (pcols)        ! sea ice fraction (fraction)
    REAL(r8), INTENT(in)  :: landfrac(pcols)        ! land fraction (fraction)
    !REAL(r8), INTENT(in)  :: ocnfrac (pcols)        ! ocean fraction (fraction)
    REAL(r8), INTENT(in)  :: pmid  (pcols,pver)     ! pressure of midpoint levels (Pa)
    !REAL(r8), INTENT(in)  :: pdel  (pcols,pver)     ! pressure diff across layer (Pa)
    REAL(r8), INTENT(in)  :: cloud (pcols,pver)     ! cloud fraction (fraction)
    REAL(r8), INTENT(in)  :: t     (pcols,pver)     ! temperature (K)
    REAL(r8), INTENT(in)  :: cldliq(pcols,pver)     ! cloud water, liquid (kg/kg)
    REAL(r8), INTENT(in)  :: cldice(pcols,pver)     ! cloud water, ice    (kg/kg)
    REAL(r8), INTENT(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)

    REAL(r8), INTENT(out) :: pvliq (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)
    REAL(r8), INTENT(out) :: pvice (pcols,pverp)    ! vertical velocity of cloud ice particles (Pa/s)
    REAL(r8), INTENT(in) :: landm(pcols)            ! land fraction ramped over water
    ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    ! Local variables
    REAL (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    REAL (r8) :: vfall                              ! settling velocity of cloud particles (m/s)
    REAL (r8) :: icice                              ! in cloud ice water content (kg/kg)
    REAL (r8) :: iciwc                              ! in cloud ice water content in g/m3
    REAL (r8) :: icefac
    REAL (r8) :: logiwc

    REAL (r8) :: rei(pcols,pver)                    ! effective radius of ice particles (microns)
    REAL (r8) :: rel(pcols,pver)                    ! effective radius of liq particles (microns)
    !REAL(r8)  pvliq2 (pcols,pverp)    ! vertical velocity of cloud liquid drops (Pa/s)

    INTEGER i,k

    REAL (r8) :: lbound, ac, bc, cc

    !-----------------------------------------------------------------------
    !--------------------- liquid fall velocity ----------------------------
    !-----------------------------------------------------------------------
    vfall = 0.0_r8 
    icice = 0.0_r8 
    iciwc = 0.0_r8 
    icefac= 0.0_r8 
    logiwc= 0.0_r8 
    ! compute air density
    DO k = 1,pver
       DO i = 1,ncol
          rei(i,k) = 0.0_r8      ! effective radius of ice particles (microns)
          rel(i,k) = 0.0_r8      ! effective radius of liq particles (microns)

          rho(i,k) = pmid(i,k) / (rair * t(i,k))
          pvliq(i,k) = 0.0_r8
       END DO
    END DO

    ! get effective radius of liquid drop
    CALL reltab(ncol,pcols,pver,t, landm, icefrac, rel, snowh)

    DO k = 1,pver
       DO i = 1,ncol
          IF (cloud(i,k) > 0.0_r8 .AND. cldliq(i,k) > 0.0_r8) THEN

             IF(OLDLIQSED)THEN
                ! oldway
                ! merge the liquid fall velocities for land and ocean (cm/s)
                ! SHOULD ALSO ACCOUNT FOR ICEFRAC
                vfall = vland*landfrac(i) + vocean*(1.0_r8-landfrac(i))
!!$          vfall = vland*landfrac(i) + vocean*ocnfrac(i) + vseaice*icefrac(i)

                ! convert the fall speed to pressure units, but do not apply the traditional
                ! negative convention for pvel.
                pvliq(i,k+1) = vfall     &
                     * 0.01_r8                 & ! cm to meters
                     * rho(i,k)*gravit        ! meters/sec to pascals/sec
             ELSE

                ! newway
                IF (rel(i,k) < 40.0_r8 ) THEN
                   vfall = 2.0_r8/9.0_r8 * rhoh2o * gravit * rel(i,k)**2 / eta  * 1.0e-12_r8  ! micons^2 -> m^2
                ELSE
                   vfall = v40 + vslope * (rel(i,k)-r40)      ! linear above 40 microns
                END IF
                ! convert the fall speed to pressure units
                ! but do not apply the traditional
                ! negative convention for pvel.
                !             pvliq2(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
                pvliq(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
             ENDIF
          END IF
       END DO
    END DO

    !-----------------------------------------------------------------------
    !--------------------- ice fall velocity -------------------------------
    !--------------------- stokes terminal velocity < 40 microns -----------
    !-----------------------------------------------------------------------
    DO k = 1,pver
       DO i = 1,ncol
         pvice(i,k) = 0.0_r8
       END DO
    END DO

    ! stokes terminal velocity
    IF (stokes) THEN

       ! get effective radius
       CALL reitab(ncol,pcols,pver, t, rei)

       DO k = 1,pver
          DO i = 1,ncol
             IF (cloud(i,k) > 0.0_r8 .AND. cldice(i,k) > 0.0_r8) THEN
                IF (rei(i,k) < 40.0_r8 ) THEN
                   vfall = 2.0_r8/9.0_r8 * rhoh2o * gravit * rei(i,k)**2 / eta  * 1.0e-12_r8  ! micons^2 -> m^2
                ELSE
                   vfall = v40 + vslope * (rei(i,k)-r40)      ! linear above 40 microns
                END IF
!!$                print *, t(i,k), rei(i,k), vfall*100.
                ! convert the fall speed to pressure units, but do not apply the traditional
                ! negative convention for pvel.
                pvice(i,k+1) = vfall * rho(i,k)*gravit        ! meters/sec to pascals/sec
             END IF
          END DO
       END DO

       RETURN
    END IF

    !-----------------------------------------------------------------------
    !--------------------- ice fall velocity -------------------------------
    !--------------------- McFarquhar and Heymsfield > icritc --------------
    !-----------------------------------------------------------------------

    ! lower bound for iciwc

    cc = 128.64_r8 
    bc = 53.242_r8
    ac = 5.4795_r8
    lbound = REAL((-bc + SQRT(bc*bc-4.0_r8*ac*cc))/(2.0_r8*ac),KIND=r8)
    lbound = 10.0_r8**lbound

    DO k = 1,pver
       DO i = 1,ncol
          IF (cloud(i,k) > 0.0_r8 .AND. cldice(i,k) > 0.0_r8) THEN

             ! compute the in-cloud ice concentration (kg/kg)
             icice = cldice(i,k) / cloud(i,k)

             ! compute the ice water content in g/m3
             iciwc = icice * rho(i,k) * 1.0e3_r8

             ! set the fall velocity (cm/s) to depend on the ice water content in g/m3,
             IF (iciwc > lbound) THEN ! need this because of log10
                logiwc = LOG10(iciwc)
                !          Median - 
                vfall = 128.64_r8 + 53.242_r8*logiwc + 5.4795_r8*logiwc**2
                !          Average - 
!!$             vfall = 122.63 + 44.111*logiwc + 4.2144*logiwc**2
             ELSE
                vfall = 0.0_r8
             END IF

             ! set ice velocity to 1 cm/s if ice mixing ratio < icritc, ramp to value
             ! calculated above at 2*icritc
             IF (icice <= icritc) THEN
                vfall = vice_small
             ELSE IF(icice < 2*icritc) THEN
                icefac = (icice-icritc)/icritc
                vfall = vice_small * (1.0_r8-icefac) + vfall * icefac
             END IF

             ! bound the terminal velocity of ice particles at high concentration
             vfall = MIN(100.0_r8, vfall)

             ! convert the fall speed to pressure units, but do not apply the traditional
             ! negative convention for pvel.
             pvice(i,k+1) = vfall     &
                  * 0.01_r8                 & ! cm to meters
                  * rho(i,k)*gravit        ! meters/sec to pascals/sec
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cld_sediment_vel




  !===============================================================================
  SUBROUTINE cld_sediment_tend (ncol, pcols,pver,pverp, dtime  ,               &
       pint   , pdel   ,                    &
       cloud  , cldliq , cldice , pvliq  , pvice  ,          &
       liqtend, icetend, wvtend , htend  , sfliq  , sfice   )

    !----------------------------------------------------------------------
    !     Apply Cloud Particle Gravitational Sedimentation to Condensate
    !----------------------------------------------------------------------

    IMPLICIT NONE

    ! Arguments
    INTEGER,  INTENT(in)  :: ncol                      ! number of colums to process
    INTEGER,  INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER,  INTENT(in) :: pverp                  ! number of vertical levels + 1
    INTEGER,  INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in)  :: dtime                     ! time step
    REAL(r8), INTENT(in)  :: pint  (pcols,pverp)       ! interfaces pressure (Pa)
    !REAL(r8), INTENT(in)  :: pmid  (pcols,pver)        ! not used !midpoint pressures (Pa)
    REAL(r8), INTENT(in)  :: pdel  (pcols,pver)        ! pressure diff across layer (Pa)
    REAL(r8), INTENT(in)  :: cloud (pcols,pver)        ! cloud fraction (fraction)
    !REAL(r8), INTENT(in)  :: t     (pcols,pver)        ! not used ! temperature (K)
    REAL(r8), INTENT(in)  :: cldliq(pcols,pver)        ! cloud liquid water (kg/kg)
    REAL(r8), INTENT(in)  :: cldice(pcols,pver)        ! cloud ice water    (kg/kg)
    REAL(r8), INTENT(in)  :: pvliq (pcols,pverp)       ! vertical velocity of liquid drops  (Pa/s)
    REAL(r8), INTENT(in)  :: pvice (pcols,pverp)       ! vertical velocity of ice particles (Pa/s)

    ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    REAL(r8), INTENT(out) :: liqtend(pcols,pver)       ! liquid condensate tend
    REAL(r8), INTENT(out) :: icetend(pcols,pver)       ! ice condensate tend
    REAL(r8), INTENT(out) :: wvtend (pcols,pver)       ! water vapor tend
    REAL(r8), INTENT(out) :: htend  (pcols,pver)       ! heating rate
    REAL(r8), INTENT(out) :: sfliq  (pcols)            ! surface flux of liquid (rain, kg/m/s)
    REAL(r8), INTENT(out) :: sfice  (pcols)            ! surface flux of ice    (snow, kg/m/s)

    ! Local variables
    REAL(r8) :: fxliq(pcols,pverp)                     ! fluxes at the interfaces, liquid (positive = down)
    REAL(r8) :: fxice(pcols,pverp)                     ! fluxes at the interfaces, ice    (positive = down)
    REAL(r8) :: cldab(pcols)                           ! cloud in layer above
    REAL(r8) :: evapliq                                ! evaporation of cloud liquid into environment
    REAL(r8) :: evapice                                ! evaporation of cloud ice into environment
    REAL(r8) :: cldovrl                                ! cloud overlap factor

    INTEGER :: i,k
    !----------------------------------------------------------------------

    ! initialize variables
    DO k = 1,pver
       DO i = 1,ncol
           fxliq  (i,k) = 0.0_r8 ! flux at interfaces (liquid)
           fxice  (i,k) = 0.0_r8 ! flux at interfaces (ice)
           liqtend(i,k) = 0.0_r8 ! condensate tend (liquid)
           icetend(i,k) = 0.0_r8 ! condensate tend (ice)
           wvtend (i,k) = 0.0_r8 ! environmental moistening
           htend  (i,k) = 0.0_r8 ! evaporative cooling
           sfliq  (i)   = 0.0_r8 ! condensate sedimentation flux out bot of column (liquid)
           sfice  (i)   = 0.0_r8 ! condensate sedimentation flux out bot of column (ice)
       END DO
    END DO

    ! fluxes at interior points
    !PRINT*,pint
    CALL getflx(ncol,pcols ,pverp,pver, pint, cldliq, pvliq, dtime, fxliq)
    CALL getflx(ncol,pcols ,pverp,pver, pint, cldice, pvice, dtime, fxice)

    ! calculate fluxes at boundaries
    DO i = 1,ncol
       fxliq(i,1) = 0.0_r8
       fxice(i,1) = 0.0_r8
       ! surface flux by upstream scheme
       fxliq(i,pverp) = cldliq(i,pver) * pvliq(i,pverp) * dtime
       fxice(i,pverp) = cldice(i,pver) * pvice(i,pverp) * dtime
    END DO

    ! filter out any negative fluxes from the getflx routine
    ! (typical fluxes are of order > 1e-3 when clouds are present)
    DO k = 2,pver
       DO i = 1,ncol
          fxliq(i,k) = MAX(0.0_r8, fxliq(i,k))
          fxice(i,k) = MAX(0.0_r8, fxice(i,k))
       END DO
    END DO

    ! Limit the flux out of the bottom of each cell to the water content in each phase.
    ! Apply mxsedfac to prevent generating very small negative cloud water/ice
    ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
    ! ***Should we include the flux in the top, to allow for thin surface layers?
    ! ***Requires simple treatment of cloud overlap, already included below.
    DO k = 1,pver
       DO i = 1,ncol
          fxliq(i,k+1) = MIN( fxliq(i,k+1), mxsedfac * cldliq(i,k) * pdel(i,k) )
          fxice(i,k+1) = MIN( fxice(i,k+1), mxsedfac * cldice(i,k) * pdel(i,k) )
!!$        fxliq(i,k+1) = min( fxliq(i,k+1), cldliq(i,k) * pdel(i,k) + fxliq(i,k))
!!$        fxice(i,k+1) = min( fxice(i,k+1), cldice(i,k) * pdel(i,k) + fxice(i,k))
!!$        fxliq(i,k+1) = min( fxliq(i,k+1), cloud(i,k) * cldliq(i,k) * pdel(i,k) )
!!$        fxice(i,k+1) = min( fxice(i,k+1), cloud(i,k) * cldice(i,k) * pdel(i,k) )
       END DO
    END DO

    ! Now calculate the tendencies assuming that condensate evaporates when
    ! it falls into environment, and does not when it falls into cloud.
    ! All flux out of the layer comes from the cloudy part.
    ! Assume maximum overlap for stratiform clouds
    !  if cloud above < cloud,  all water falls into cloud below
    !  if cloud above > cloud,  water split between cloud  and environment
    DO k = 1,pver
       cldab(1:ncol) = 0.0_r8
       DO i = 1,ncol
          ! cloud overlap cloud factor
          IF(cloud(i,k) > 0.0_r8)THEN
             cldovrl  = MIN( cloud(i,k) / (cldab(i)+0.0001_r8), 1.0_r8 )
          ELSE
             cldovrl  = 0.0_r8
          END IF   
          cldab(i) = cloud(i,k)
          ! evaporation into environment cause moistening and cooling
          evapliq = fxliq(i,k) * (1.0_r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s)
          evapice = fxice(i,k) * (1.0_r8-cldovrl) / (dtime * pdel(i,k))  ! into env (kg/kg/s)
          wvtend(i,k) = evapliq + evapice                          ! evaporation into environment (kg/kg/s)
          htend (i,k) = -latvap*evapliq -(latvap+latice)*evapice   ! evaporation (W/kg)
          ! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)
          liqtend(i,k)  = (fxliq(i,k)*cldovrl - fxliq(i,k+1)) / (dtime * pdel(i,k))
          icetend(i,k)  = (fxice(i,k)*cldovrl - fxice(i,k+1)) / (dtime * pdel(i,k))
       END DO
    END DO

    ! convert flux out the bottom to mass units Pa -> kg/m2/s
    DO i = 1,ncol
       sfliq(i) = fxliq(i,pverp) / (dtime*gravit)
       sfice(i) = fxice(i,pverp) / (dtime*gravit)
    END DO
    RETURN
  END SUBROUTINE cld_sediment_tend


  !===============================================================================
  SUBROUTINE getflx(ncol,pcols ,pverp,pver,xw, phi, vel, deltat, flux)

    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......


    IMPLICIT NONE

    INTEGER, INTENT(in) :: ncol                      ! number of colums to process
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                  ! number of vertical levels + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    INTEGER i
    INTEGER k

    REAL (r8), INTENT(in) :: vel(pcols,pverp)
    REAL (r8) flux(pcols,pverp)
    REAL (r8) xw(pcols,pverp)
    REAL (r8) psi(pcols,pverp)
    REAL (r8), INTENT(in) :: phi(pcols,pverp-1)
    REAL (r8) fdot(pcols,pverp)
    !REAL (r8) xx(pcols)
    REAL (r8) fxdot(pcols)
    REAL (r8) fxdd(pcols)

    REAL (r8) psistar(pcols)
    REAL (r8) deltat

    REAL (r8) xxk(pcols,pver)

    DO i = 1,ncol
       !        integral of phi
       psi(i,1) = 0.0_r8
       !        fluxes at boundaries
       flux(i,1) = 0.0_r8
       flux(i,pverp) = 0.0_r8
    END DO

    !     integral function
    DO k = 2,pverp
       DO i = 1,ncol
          !PRINT*,xw(i,k),xw(i,k-1),phi(i,k-1)
          psi(i,k) = phi(i,k-1)*(xw(i,k)-xw(i,k-1)) + psi(i,k-1)
       END DO
    END DO


    !     calculate the derivatives for the interpolating polynomial
    CALL cfdotmc_pro (ncol,pcols,pverp,pver, xw, psi, fdot)

    !  NEW WAY
    !     calculate fluxes at interior pts
    DO k = 2,pver
       DO i = 1,ncol
          xxk(i,k) = xw(i,k)-vel(i,k)*deltat
       END DO
    END DO
    DO k = 2,pver
       CALL cfint2(ncol, pcols,pverp, xw, psi, fdot, xxk(1,k), fxdot, fxdd, psistar)
       DO i = 1,ncol
          flux(i,k) = (psi(i,k)-psistar(i))
       END DO
    END DO


    RETURN
  END SUBROUTINE getflx


  REAL(r8) FUNCTION minmod( a,b )
    REAL(KIND=r8), INTENT(IN   ) ::  a,b 
    minmod = 0.5_r8*(SIGN(1.0_r8,a) + SIGN(1.0_r8,b))*MIN(ABS(a),ABS(b))
  END FUNCTION minmod

  REAL(r8) FUNCTION medan( a,b,c )
    REAL(KIND=r8), INTENT(IN   ) ::  a,b,c
    !REAL(KIND=r8) :: minmod
    medan = a + minmod(b-a,c-a)
  END FUNCTION medan

  !###############################################################################

  SUBROUTINE cfint2 (ncol,pcols,pverp, x, f, fdot, xin, fxdot, fxdd, psistar)


    IMPLICIT NONE

    ! input
    INTEGER ncol                      ! number of colums to process
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                  ! number of vertical levels + 1

    REAL (r8), INTENT(IN   ) ::  x(pcols, pverp)
    REAL (r8), INTENT(IN   ) ::  f(pcols, pverp)
    REAL (r8), INTENT(IN   ) ::  fdot(pcols, pverp)
    REAL (r8), INTENT(IN   ) ::  xin(pcols)

    ! output
    REAL (r8), INTENT(OUT  ) ::  fxdot(pcols)
    REAL (r8), INTENT(OUT  ) ::  fxdd(pcols)
    REAL (r8), INTENT(OUT  ) ::  psistar(pcols)

    INTEGER :: i
    INTEGER :: k
    INTEGER :: intz(pcols)
    REAL (r8) :: dx
    REAL (r8) :: s
    REAL (r8) :: c2
    REAL (r8) :: c3
    REAL (r8) :: xx
    !real (r8) xinf
    REAL (r8) :: psi1, psi2, psi3, psim
    REAL (r8) :: cfint
    REAL (r8) :: cfnew
    REAL (r8) :: xins(pcols)

    !     the minmod function 
    !REAL (r8) :: a, b, c
    !REAL (r8) :: minmod
    !REAL (r8) :: medan
    LOGICAL :: found_error

    !minmod(a,b) = 0.5_r8*(SIGN(1.0_r8,a) + SIGN(1.0_r8,b))*MIN(ABS(a),ABS(b))
    !medan(a,b,c) = a + minmod(b-a,c-a)

    DO i = 1,ncol
       xins(i) = medan(x(i,1), xin(i), x(i,pverp))
       intz(i) = 0
    END DO

    ! first find the interval 
    DO k =  1,pverp-1
       DO i = 1,ncol
          IF ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).GE.0.0_r8) THEN
             intz(i) = k
          ENDIF
       END DO
    END DO

    found_error=.FALSE.
    DO i = 1,ncol
       IF (intz(i).EQ.0) found_error=.TRUE.
    END DO
    IF(found_error) THEN
       DO i = 1,ncol
          IF (intz(i).EQ.0) THEN
             WRITE (6,*) ' interval was not found for col i ', i
             CALL endrun('CFINT2')
          ENDIF
       END DO
    ENDIF

    ! now interpolate
    DO i = 1,ncol
       k = intz(i)
       dx = (x(i,k+1)-x(i,k))
       s = (f(i,k+1)-f(i,k))/dx
       c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
       c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
       xx = (xins(i)-x(i,k))
       fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
       fxdd(i) = 6*c3*xx + 2*c2
       cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)

       !        limit the interpolant
       psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
       IF (k.EQ.1) THEN
          psi2 = f(i,1)
       ELSE
          psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
       ENDIF
       IF (k+1.EQ.pverp) THEN
          psi3 = f(i,pverp)
       ELSE
          psi3 = f(i,k+1) - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
       ENDIF
       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)
       IF (ABS(cfnew-cfint)/(ABS(cfnew)+ABS(cfint)+1.e-36_r8)  .GT.0.03_r8) THEN
          !     CHANGE THIS BACK LATER!!!
          !     $        .gt..1) then


          !     UNCOMMENT THIS LATER!!!
          !            write (6,*) ' cfint2 limiting important ', cfint, cfnew


       ENDIF
       psistar(i) = cfnew
    END DO

    RETURN
  END SUBROUTINE cfint2


  !##############################################################################

  SUBROUTINE cfdotmc_pro (ncol,pcols,pverp,pver, x, f, fdot)

    !     prototype version; eventually replace with final SPITFIRE scheme

    !     calculate the derivative for the interpolating polynomial
    !     multi column version


    IMPLICIT NONE

    ! input
    INTEGER, INTENT(in) :: ncol                 ! number of colums to process
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: pverp                  ! number of vertical levels + 1

    REAL (r8) x(pcols, pverp)
    REAL (r8) f(pcols, pverp)
    ! output
    REAL (r8) fdot(pcols, pverp)          ! derivative at nodes

    ! assumed variable distribution
    !     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
    !     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
    !     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
    !     .........d2.......d3.......d4.......d5.........     2,pver points
    !     .........s2.......s3.......s4.......s5.........     2,pver points
    !     .............dh2......dh3......dh4.............     2,pver-1 points
    !     .............eh2......eh3......eh4.............     2,pver-1 points
    !     ..................e3.......e4..................     3,pver-1 points
    !     .................ppl3......ppl4................     3,pver-1 points
    !     .................ppr3......ppr4................     3,pver-1 points
    !     .................t3........t4..................     3,pver-1 points
    !     ................fdot3.....fdot4................     3,pver-1 points


    ! work variables


    INTEGER i
    INTEGER k

    !REAL (r8) a                    ! work var
    !REAL (r8) b                    ! work var
    !REAL (r8) c                    ! work var
    REAL (r8) s(pcols,pverp)             ! first divided differences at nodes
    REAL (r8) sh(pcols,pverp)            ! first divided differences between nodes
    REAL (r8) d(pcols,pverp)             ! second divided differences at nodes
    REAL (r8) dh(pcols,pverp)            ! second divided differences between nodes
    REAL (r8) e(pcols,pverp)             ! third divided differences at nodes
    REAL (r8) eh(pcols,pverp)            ! third divided differences between nodes
    REAL (r8) pp                   ! p prime
    REAL (r8) ppl(pcols,pverp)           ! p prime on left
    REAL (r8) ppr(pcols,pverp)           ! p prime on right
    REAL (r8) qpl
    REAL (r8) qpr
    REAL (r8) ttt
    REAL (r8) t
    REAL (r8) tmin
    REAL (r8) tmax
    REAL (r8) delxh(pcols,pverp)


    !     the minmod function 
    !REAL (r8) minmod
    !REAL (r8) medan
    !minmod(a,b) = 0.5_r8*(SIGN(1.0_r8,a) + SIGN(1.0_r8,b))*MIN(ABS(a),ABS(b))
    !medan(a,b,c) = a + minmod(b-a,c-a)

    DO k = 1,pver


       !        first divided differences between nodes
       DO i = 1, ncol
          delxh(i,k) = (x(i,k+1)-x(i,k))
          sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
       END DO

       !        first and second divided differences at nodes
       IF (k.GE.2) THEN
          DO i = 1,ncol
             d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
             s(i,k) = minmod(sh(i,k),sh(i,k-1))
          END DO
       ENDIF
    END DO

    !     second and third divided diffs between nodes
    DO k = 2,pver-1
       DO i = 1, ncol
          eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
          dh(i,k) = minmod(d(i,k),d(i,k+1))
       END DO
    END DO

    !     treat the boundaries
    DO i = 1,ncol
       e(i,2) = eh(i,2)
       e(i,pver) = eh(i,pver-1)
       !        outside level
       fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)  &
            - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
       fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
       fdot(i,pverp) = sh(i,pver) + d(i,pver)*delxh(i,pver)  &
            + eh(i,pver-1)*delxh(i,pver)*(x(i,pverp)-x(i,pver-1))
       fdot(i,pverp) = minmod(fdot(i,pverp),3*sh(i,pver))
       !        one in from boundary
       fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
       fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
       fdot(i,pver) = sh(i,pver) - d(i,pver)*delxh(i,pver)   &
            - eh(i,pver-1)*delxh(i,pver)*delxh(i,pver-1)
       fdot(i,pver) = minmod(fdot(i,pver),3*s(i,pver))
    END DO


    DO k = 3,pver-1
       DO i = 1,ncol
          e(i,k) = minmod(eh(i,k),eh(i,k-1))
       END DO
    END DO



    DO k = 3,pver-1

       DO i = 1,ncol

          !           p prime at k-0.5_r8
          ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  
          !           p prime at k+0.5_r8
          ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)

          t = minmod(ppl(i,k),ppr(i,k))

          !           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
          pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 

          !           quartic estimate of fdot
          fdot(i,k) = pp                            &
               - delxh(i,k-1)*delxh(i,k)            &
               *(  eh(i,k-1)*(x(i,k+2)-x(i,k  ))    &
               + eh(i,k  )*(x(i,k  )-x(i,k-2))      &
               )/(x(i,k+2)-x(i,k-2))

          !           now limit it
          qpl = sh(i,k-1)       &
               + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)*(x(i,k)-x(i,k-2)), &
               d(i,k)  -e(i,k)*delxh(i,k))
          qpr = sh(i,k)         &
               + delxh(i,k  )*minmod(d(i,k)  +e(i,k)*delxh(i,k-1),        &
               d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

          fdot(i,k) = medan(fdot(i,k), qpl, qpr)

          ttt = minmod(qpl, qpr)
          tmin = MIN(0.0_r8,3*s(i,k),1.5_r8*t,ttt)
          tmax = MAX(0.0_r8,3*s(i,k),1.5_r8*t,ttt)

          fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))

       END DO

    END DO

    RETURN
  END SUBROUTINE cfdotmc_pro
  !END MODULE PKG_CLD_SEDIMENT




  ! MODULE PKG_CLDOPTICS
  !===============================================================================
  SUBROUTINE cldefr( &
       ncol    , &!INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
       pcols   , &!INTEGER, INTENT(in) :: pcols                  ! number of columns (max)
       pver    , &!INTEGER, INTENT(in) :: pver                  ! number of vertical levels
       t       , &!REAL(r8), INTENT(in) :: t       (pcols,pver)        ! Temperature
       rel     , &!REAL(r8), INTENT(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
       rei     , &!REAL(r8), INTENT(out) :: rei(pcols,pver)      ! Ice effective drop size (microns)
       landm   , &!REAL(r8), INTENT(in) :: landm   (pcols)
       icefrac , &!REAL(r8), INTENT(in) :: icefrac (pcols)       ! Ice fraction
       snowh     )!REAL(r8), INTENT(in) :: snowh   (pcols)         ! Snow depth over land, water equivalent (m)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud water and ice particle size 
    ! 
    ! Method: 
    ! use empirical formulas to construct effective radii
    ! 
    ! Author: J.T. Kiehl, B. A. Boville, P. Rasch
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                  ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: icefrac (pcols)       ! Ice fraction
    REAL(r8), INTENT(in) :: t       (pcols,pver)        ! Temperature
    REAL(r8), INTENT(in) :: landm   (pcols)
    REAL(r8), INTENT(in) :: snowh   (pcols)         ! Snow depth over land, water equivalent (m)
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
    REAL(r8), INTENT(out) :: rei(pcols,pver)      ! Ice effective drop size (microns)
    !
    rel=0.0_r8
    rei=0.0_r8
    !++pjr
    ! following Kiehl
    CALL reltab(ncol,pcols,pver,  t,  landm, icefrac, rel, snowh)

    ! following Kristjansson and Mitchell
    CALL reitab(ncol,pcols,pver, t, rei)
    !--pjr
    !
    !
    RETURN
  END SUBROUTINE cldefr


  !===============================================================================
  SUBROUTINE cldems(ncol ,pcols ,pver  ,clwp    ,fice    ,rei     ,emis    )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud emissivity using cloud liquid water path (g/m**2)
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J.T. Kiehl
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE
    !------------------------------Parameters-------------------------------
    !
    REAL(r8),PARAMETER :: kabsl  = 0.090361_r8                ! longwave liquid absorption coeff (m**2/g)
    !
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                    ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                  ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: clwp(pcols,pver)       ! cloud liquid water path (g/m**2)
    REAL(r8), INTENT(in) :: rei(pcols,pver)        ! ice effective drop size (microns)
    REAL(r8), INTENT(in) :: fice(pcols,pver)       ! fractional ice content within cloud
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: emis(pcols,pver)       ! cloud emissivity (fraction)
    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER i,k                 ! longitude, level indices
    REAL(r8) kabs                   ! longwave absorption coeff (m**2/g)
    REAL(r8) kabsi                  ! ice absorption coefficient
    !
    !-----------------------------------------------------------------------
    !
    DO k=1,pver
       DO i=1,ncol
          kabsi = 0.005_r8 + 1.0_r8/rei(i,k)
          kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k)
          emis(i,k) = 1.0_r8 - EXP(-1.66*kabs*clwp(i,k))
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE cldems



  !===============================================================================
  SUBROUTINE cldovrlap(ncol    ,pcols,pverp,pver,pint    ,cld     ,nmxrgn  ,pmxrgn  )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Partitions each column into regions with clouds in neighboring layers.
    ! This information is used to implement maximum overlap in these regions
    ! with random overlap between them.
    ! On output,
    !    nmxrgn contains the number of regions in each column
    !    pmxrgn contains the interface pressures for the lower boundaries of
    !           each region! 
    ! Method: 

    ! 
    ! Author: W. Collins
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                 ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(in) :: pverp                 ! pver + 1

    REAL(r8), INTENT(in) :: pint(pcols,pverp)   ! Interface pressure
    REAL(r8), INTENT(in) :: cld(pcols,pver)     ! Fractional cloud cover
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: pmxrgn(pcols,pverp)! Maximum values of pressure for each
    !    maximally overlapped region.
    !    0->pmxrgn(i,1) is range of pressure for
    !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
    !    2nd region, etc
    INTEGER nmxrgn(pcols)                    ! Number of maximally overlapped regions
    !
    !---------------------------Local variables-----------------------------
    !
    INTEGER i                    ! Longitude index
    INTEGER k                    ! Level index
    INTEGER n                    ! Max-overlap region counter

    REAL(r8) pnm(pcols,pverp)    ! Interface pressure

    LOGICAL cld_found            ! Flag for detection of cloud
    LOGICAL cld_layer(pver)      ! Flag for cloud in layer
    !
    !------------------------------------------------------------------------
    !

    DO i = 1, ncol
       cld_found = .FALSE.
       cld_layer(:) = cld(i,:) > 0.0_r8
       pmxrgn(i,:) = 0.0_r8
       pnm(i,:)=pint(i,:)*10._r8
       n = 1
       DO k = 1, pver
          IF (cld_layer(k) .AND.  .NOT. cld_found) THEN
             cld_found = .TRUE.
          ELSE IF ( .NOT. cld_layer(k) .AND. cld_found) THEN
             cld_found = .FALSE.
             IF (COUNT(cld_layer(k:pver)) == 0) THEN
                EXIT
             ENDIF
             pmxrgn(i,n) = pnm(i,k)
             n = n + 1
          ENDIF
       END DO
       pmxrgn(i,n) = pnm(i,pverp)
       nmxrgn(i) = n
    END DO

    RETURN
  END SUBROUTINE cldovrlap

  !===============================================================================
  SUBROUTINE cldclw(ncol    ,pcols,pverp,pver, zi      ,clwp    ,tpw     ,hl      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Evaluate cloud liquid water path clwp (g/m**2)
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J.T. Kiehl
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pverp                 ! pver + 1
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: zi(pcols,pverp)      ! height at layer interfaces(m)
    REAL(r8), INTENT(in) :: tpw(pcols)           ! total precipitable water (mm)
    !
    ! Output arguments
    !
    REAL(r8) clwp(pcols,pver)     ! cloud liquid water path (g/m**2)
    REAL(r8) hl(pcols)            ! liquid water scale height
    REAL(r8) rhl(pcols)           ! 1/hl

    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER i,k               ! longitude, level indices
    REAL(r8) clwc0                ! reference liquid water concentration (g/m**3)
    REAL(r8) emziohl(pcols,pverp) ! exp(-zi/hl)
    !
    !-----------------------------------------------------------------------
    !
    ! Set reference liquid water concentration
    !
    clwc0 = 0.21_r8
    !
    ! Diagnose liquid water scale height from precipitable water
    !
    DO i=1,ncol
       hl(i)  = 700.0_r8*LOG(MAX(tpw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !
    ! Evaluate cloud liquid water path (vertical integral of exponential fn)
    !
    DO k=1,pverp
       DO i=1,ncol
          emziohl(i,k) = EXP(-zi(i,k)*rhl(i))
       END DO
    END DO
    DO k=1,pver
       DO i=1,ncol
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k+1) - emziohl(i,k))
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE cldclw


  !===============================================================================
  SUBROUTINE reltab(ncol,pcols,pver, t, landm, icefrac, rel, snowh)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud water size
    ! 
    ! Method: 
    ! analytic formula following the formulation originally developed by J. T. Kiehl
    ! 
    ! Author: Phil Rasch
    ! 
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER  , INTENT(in) :: ncol
    INTEGER, INTENT(in) :: pcols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(r8), INTENT(in) :: icefrac(pcols)       ! Ice fraction
    REAL(r8), INTENT(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    REAL(r8), INTENT(in) :: landm(pcols)         ! Land fraction ramping to zero over ocean
    REAL(r8), INTENT(in) :: t(pcols,pver)        ! Temperature

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER i,k               ! Lon, lev indices
    REAL(r8) rliqland         ! liquid drop size if over land
    REAL(r8) rliqocean        ! liquid drop size if over ocean
    REAL(r8) rliqice          ! liquid drop size if over sea ice
    !
    !-----------------------------------------------------------------------
    !
    rliqocean = 14.0_r8
    rliqice   = 14.0_r8
    rliqland  = 8.0_r8
    DO k=1,pver
       DO i=1,ncol
          ! jrm Reworked effective radius algorithm
          ! Start with temperature-dependent value appropriate for continental air
          ! Note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * MIN(1.0_r8,MAX(0.0_r8,(tmelt-t(i,k))*0.05))
          ! Modify for snow depth over land
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * MIN(1.0_r8,MAX(0.0_r8,snowh(i)*10.0_r8))
          ! Ramp between polluted value over land to clean value over ocean.
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * MIN(1.0_r8,MAX(0.0_r8,1.0_r8-landm(i)))
          ! Ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k) = rel(i,k) + (rliqice-rel(i,k)) * MIN(1.0_r8,MAX(0.0_r8,icefrac(i)))
          ! end jrm
       END DO
    END DO
  END SUBROUTINE reltab



  !===============================================================================
  SUBROUTINE reitab(ncol,pcols,pver, t, re)
    !

    INTEGER  , INTENT(in) :: ncol
    INTEGER  , INTENT(in) :: pcols
    INTEGER  , INTENT(in) :: pver
    REAL(r8), INTENT(in ) :: t(pcols,pver)
    REAL(r8), INTENT(out) :: re(pcols,pver)
    REAL(r8) :: corr
    INTEGER :: i
    INTEGER :: k
    INTEGER :: index
    !
    !       Tabulated values of re(T) in the temperature interval
    !       180 K -- 274 K; hexagonal columns assumed:
    !
    REAL(KIND=r8), PARAMETER :: retab(95)=(/                                                 &
         5.92779_r8, 6.26422_r8, 6.61973_r8, 6.99539_r8, 7.39234_r8,        &
         7.81177_r8, 8.25496_r8, 8.72323_r8, 9.21800_r8, 9.74075_r8, 10.2930_r8,        &
         10.8765_r8, 11.4929_r8, 12.1440_r8, 12.8317_r8, 13.5581_r8, 14.2319_r8,         &
         15.0351_r8, 15.8799_r8, 16.7674_r8, 17.6986_r8, 18.6744_r8, 19.6955_r8,        &
         20.7623_r8, 21.8757_r8, 23.0364_r8, 24.2452_r8, 25.5034_r8, 26.8125_r8,        &
         27.7895_r8, 28.6450_r8, 29.4167_r8, 30.1088_r8, 30.7306_r8, 31.2943_r8,         &
         31.8151_r8, 32.3077_r8, 32.7870_r8, 33.2657_r8, 33.7540_r8, 34.2601_r8,         &
         34.7892_r8, 35.3442_r8, 35.9255_r8, 36.5316_r8, 37.1602_r8, 37.8078_r8,        &
         38.4720_r8, 39.1508_r8, 39.8442_r8, 40.5552_r8, 41.2912_r8, 42.0635_r8,        &
         42.8876_r8, 43.7863_r8, 44.7853_r8, 45.9170_r8, 47.2165_r8, 48.7221_r8,        &
         50.4710_r8, 52.4980_r8, 54.8315_r8, 57.4898_r8, 60.4785_r8, 63.7898_r8,        &
         65.5604_r8, 71.2885_r8, 75.4113_r8, 79.7368_r8, 84.2351_r8, 88.8833_r8,        &
         93.6658_r8, 98.5739_r8, 103.603_r8, 108.752_r8, 114.025_r8, 119.424_r8,         &
         124.954_r8, 130.630_r8, 136.457_r8, 142.446_r8, 148.608_r8, 154.956_r8,        &
         161.503_r8, 168.262_r8, 175.248_r8, 182.473_r8, 189.952_r8, 197.699_r8,        &
         205.728_r8, 214.055_r8, 222.694_r8, 231.661_r8, 240.971_r8, 250.639_r8/)        
    !
    !
    DO k=1,pver
       DO i=1,ncol
          index = INT(t(i,k)-179.0_r8)
          index = MIN(MAX(index,1),94)
          corr = t(i,k) - INT(t(i,k))
          re(i,k) = retab(index)*(1.0_r8-corr)                &
               +retab(index+1)*corr
          !           re(i,k) = amax1(amin1(re(i,k),30.0_r8),10.0_r8)
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE reitab

  ! END MODULE PKG_CLDOPTICS



  !MODULE cloud_fraction

  SUBROUTINE cldfrc_init(pver,plat)
    !
    ! Purpose:
    ! Initialize cloud fraction run-time parameters
    !
    ! Author: J. McCaa
    !    
    INTEGER, INTENT(IN   ) :: pver  ! number of vertical levels        
    INTEGER, INTENT(IN   ) ::  plat! number of latitudes

    IF ( dycore_is ('LR') ) THEN
       IF ( get_resolution(plat) == '1x1.25' ) THEN
          rhminl = .88_r8     ! minimum rh for low stable clouds
          rhminh = .77_r8     ! minimum rh for high stable clouds
          sh1 = 0.04_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.10_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud
       ELSEIF ( get_resolution(plat) == '4x5' .AND. pver == 66 ) THEN
          rhminl = .90_r8     ! minimum rh for low stable clouds
          rhminh = .90_r8     ! minimum rh for high stable clouds
          sh1 = 0.04_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.10_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud
       ELSE
          rhminl = .90_r8     ! minimum rh for low stable clouds
          rhminh = .80_r8     ! minimum rh for high stable clouds
          sh1 = 0.04_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.10_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud
       ENDIF
    ELSE
       IF ( get_resolution(plat) == 'T85' ) THEN
          rhminl = .91_r8     ! minimum rh for low stable clouds
          rhminh = .70_r8     ! minimum rh for high stable clouds
          sh1 = 0.07_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.14_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 250.e2_r8  ! top of area defined to be mid-level cloud
       ELSEIF ( get_resolution(plat) == 'T62' ) THEN
          rhminl = .91_r8     ! minimum rh for low stable clouds
          rhminh = .70_r8     ! minimum rh for high stable clouds
          sh1 = 0.07_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.14_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 250.e2_r8  ! top of area defined to be mid-level cloud
       ELSEIF ( get_resolution(plat) == 'T31' ) THEN
          rhminl = .88_r8     ! minimum rh for low stable clouds
          rhminh = .80_r8     ! minimum rh for high stable clouds
          sh1 = 0.07_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.14_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud
       ELSE
          rhminl = .90_r8     ! minimum rh for low stable clouds
          rhminh = .80_r8     ! minimum rh for high stable clouds
          sh1 = 0.07_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.14_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud
       ENDIF
    ENDIF
          rhminl = .97_r8     ! minimum rh for low stable clouds
          rhminh = .93_r8     ! minimum rh for high stable clouds
          sh1 = 0.07_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.10_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud

          rhminl = 0.95_r8     ! minimum rh for low stable clouds
          rhminh = 0.70_r8     ! minimum rh for high stable clouds
          sh1 = 0.07_r8       ! parameters for shallow convection cloud fraction
          sh2 = 500.0_r8      ! parameters for shallow convection cloud fraction
          dp1 = 0.14_r8       ! parameters for deep convection cloud fraction
          dp2 = 500.0_r8      ! parameters for deep convection cloud fraction
          premit = 250.e2_r8  ! top of area defined to be mid-level cloud

  END SUBROUTINE cldfrc_init

  SUBROUTINE cldfrc( &
       pcols   , &!INTEGER, INTENT(in) ::  pcols        ! number of columns  (max)
       pver    , &!INTEGER, INTENT(in) ::  pver        ! number of vertical levels
       pverp   , &!INTEGER, INTENT(in) ::  pverp        ! pver + 1
       ncol    , &!INTEGER, INTENT(in) :: ncol                   ! number of atmospheric columns
       pnot    , &!REAL(r8), INTENT(in) :: pnot    (pcols)              ! reference pressure
       pmid    , &!REAL(r8), INTENT(in) :: pmid    (pcols,pver)      ! midpoint pressures
       temp    , &!REAL(r8), INTENT(in) :: temp    (pcols,pver)      ! temperature
       q       , &!REAL(r8), INTENT(in) :: q           (pcols,pver)      ! specific humidity
       phis    , &!REAL(r8), INTENT(in) :: phis    (pcols)             ! surface geopotential 
       cloud   , &!REAL(r8), INTENT(out) :: cloud  (pcols,pver)      ! cloud fraction
       clc     , &!REAL(r8), INTENT(out) :: clc    (pcols)           ! column convective cloud amount 
       cmfmc   , &!REAL(r8), INTENT(in) :: cmfmc   (pcols,pverp)     ! convective mass flux--m sub c
       cmfmc2  , &!REAL(r8), INTENT(in) :: cmfmc2  (pcols,pverp)     ! shallow convective mass flux--m sub c
       landfrac, &!REAL(r8), INTENT(in) :: landfrac(pcols)             ! Land fraction 
       snowh   , &!REAL(r8), INTENT(in) :: snowh   (pcols)             ! snow depth (liquid water equivalent) ! sea surface temperature
       concld  , &!REAL(r8), INTENT(out) concld (pcols,pver)   ! convective cloud cover
       cldst   , &!REAL(r8), INTENT(in) :: cldst     (pcols,pver)      ! not used ! detrainment rate from deep convection
       sst     , &!REAL(r8), INTENT(in) :: sst         (pcols)           ! sea surface temperature
       ps      , &!REAL(r8), INTENT(in) :: ps         (pcols)           ! surface pressure
       !zdu     , &!REAL(r8), INTENT(in) :: zdu         (pcols,pver)           ! not used ! detrainment rate from deep convection
       ocnfrac , &!REAL(r8), INTENT(in) :: ocnfrac (pcols)             ! Ocean fraction
       rhu00   , &!REAL(r8), INTENT(out) :: rhu00  (pcols,pver)      ! RH threshold for cloud
       relhum  , &!REAL(r8), INTENT(out) :: relhum (pcols,pver)      ! RH 
       dindex    )!INTEGER, INTENT(in) :: dindex                 ! 0 or 1 to perturb rh
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud fraction 
    ! 
    ! 
    ! Method: 
    ! This calculate cloud fraction using a relative humidity threshold
    ! The threshold depends upon pressure, and upon the presence or absence 
    ! of convection as defined by a reasonably large vertical mass flux 
    ! entering that layer from below.
    ! 
    ! Author: Many. Last modified by Jim McCaa
    ! 
    !-----------------------------------------------------------------------
    !use ppgrid
    !use wv_saturation, only: aqsat
    !use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
    !use dycore,        only: dycore_is, get_resolution

    IMPLICIT NONE


    !REAL(r8), PARAMETER :: pnot = 1.e5_r8       ! reference pressure
    REAL(r8), PARAMETER :: lapse = 6.5e-3_r8    ! U.S. Standard Atmsophere lapse rate
    REAL(r8), PARAMETER :: premib = 750.e2_r8   ! bottom pressure bound of middle cloud
    REAL(r8), PARAMETER :: pretop = 1.0e2_r8    ! pressure bounding high cloud
    !
    ! Arguments
    !
    INTEGER, INTENT(in) ::  pcols      ! number of columns  (max)
    INTEGER, INTENT(in) ::  pver       ! number of vertical levels
    INTEGER, INTENT(in) ::  pverp      ! pver + 1

    INTEGER, INTENT(in) :: ncol                   ! number of atmospheric columns
    INTEGER, INTENT(in) :: dindex                 ! 0 or 1 to perturb rh
    REAL(r8), INTENT(in) :: pnot    (pcols)           ! reference pressure
    REAL(r8), INTENT(in) :: pmid    (pcols,pver)      ! midpoint pressures
    REAL(r8), INTENT(in) :: temp    (pcols,pver)      ! temperature
    REAL(r8), INTENT(in) :: q       (pcols,pver)      ! specific humidity
    REAL(r8), INTENT(in) :: cmfmc   (pcols,pverp)     ! convective mass flux--m sub c
    REAL(r8), INTENT(in) :: cmfmc2  (pcols,pverp)     ! shallow convective mass flux--m sub c
    REAL(r8), INTENT(in) :: snowh   (pcols)           ! snow depth (liquid water equivalent)
    REAL(r8), INTENT(in) :: landfrac(pcols)           ! Land fraction
    REAL(r8), INTENT(in) :: ocnfrac (pcols)           ! Ocean fraction
    REAL(r8), INTENT(in) :: sst     (pcols)           ! sea surface temperature
    REAL(r8), INTENT(in) :: ps      (pcols)           ! surface pressure
    !REAL(r8), INTENT(in) :: zdu     (pcols,pver)      ! not used ! detrainment rate from deep convection
    REAL(r8), INTENT(in) :: phis    (pcols)           ! surface geopotential
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: cloud  (pcols,pver)      ! cloud fraction
    REAL(r8), INTENT(out) :: clc    (pcols)           ! column convective cloud amount
    REAL(r8), INTENT(out) :: cldst  (pcols,pver)      ! cloud fraction
    REAL(r8), INTENT(out) :: rhu00  (pcols,pver)      ! RH threshold for cloud
    REAL(r8), INTENT(out) :: relhum (pcols,pver)      ! RH 
    !      real(r8) dmudp                 ! measure of mass detraining in a layer
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(r8), INTENT(out) :: concld (pcols,pver)   ! convective cloud cover
    !real(r8) cld                   ! not used ! intermediate scratch variable (low cld)
    REAL(r8) dthdpmn(pcols)        ! most stable lapse rate below 750 mb
    REAL(r8) dthdp                 ! lapse rate (intermediate variable)
    REAL(r8) es     (pcols,pver)   ! saturation vapor pressure
    REAL(r8) qs     (pcols,pver)   ! saturation specific humidity
    REAL(r8) rhwght                ! weighting function for rhlim transition
    REAL(r8) rh     (pcols,pver)   ! relative humidity
    REAL(r8) rhdif                 ! intermediate scratch variable
    REAL(r8) strat                 ! intermediate scratch variable
    REAL(r8) theta(pcols,pver)     ! potential temperature
    REAL(r8) rhlim                 ! local rel. humidity threshold estimate
    REAL(r8) coef1                 ! coefficient to convert mass flux to mb/d
    !real(r8) clrsky(pcols)         ! not used ! temporary used in random overlap calc
    REAL(r8) rpdeli(pcols,pver-1)  ! 1./(pmid(k+1)-pmid(k))
    REAL(r8) rhpert                ! the specified perturbation to rh
    REAL(r8) deepcu                ! deep convection cloud fraction
    REAL(r8) shallowcu             ! shallow convection cloud fraction

    LOGICAL cldbnd(pcols)          ! region below high cloud boundary

    INTEGER i, ierror, k           ! column, level indices
    INTEGER kp1
    INTEGER kdthdp(pcols)
    INTEGER numkcld                ! number of levels in which to allow clouds

    REAL(r8) thetas(pcols)                    ! ocean surface potential temperature
    !real(r8) :: clat(pcols)                   ! not used ! current latitudes(radians)
    !real(r8) :: clon(pcols)                   ! not used ! current longitudes(radians)
    !
    ! Statement functions
    !
    LOGICAL land(pcols)


    !call get_rlat_all_p(lchnk, ncol, clat)
    !call get_rlon_all_p(lchnk, ncol, clon)
    DO k=1,pver
       DO i=1,ncol
          cloud  (i,k)=0.0_r8      ! cloud fraction
          cldst  (i,k)=0.0_r8      ! cloud fraction
          rhu00  (i,k)=0.0_r8      ! RH threshold for cloud
          relhum (i,k)=0.0_r8      ! RH 
          concld (i,k)=0.0_r8   ! convective cloud cover
          es     (i,k)=0.0_r8   ! saturation vapor pressure
          qs     (i,k)=0.0_r8   ! saturation specific humidity
          rh     (i,k)=0.0_r8   ! relative humidity
          theta  (i,k)=0.0_r8     ! potential temperature
       END DO
    END DO

    DO k=1,pver-1
       DO i=1,ncol
          rpdeli(i,k)=0.0_r8   ! 1./(pmid(k+1)-pmid(k))
       END DO
    END DO

    !      real(r8) dmudp                 ! measure of mass detraining in a layer
    !
    !---------------------------Local workspace-----------------------------
    !
    !real(r8) cld                   ! not used ! intermediate scratch variable (low cld)
    DO i=1,ncol
       land(i) = NINT(landfrac(i)) == 1
       dthdpmn(i) =0.0_r8        ! most stable lapse rate below 750 mb
       clc   (i) =0.0_r8 ! column convective cloud amount
       kdthdp(i) =0
       thetas(i) =0.0_r8                    ! ocean surface potential temperature
    END DO 
    dthdp=0.0_r8! lapse rate (intermediate variable)
    rhwght=0.0_r8! weighting function for rhlim transition
    rhdif=0.0_r8! intermediate scratch variable
    strat=0.0_r8! intermediate scratch variable
    rhlim=0.0_r8! local rel. humidity threshold estimate
    coef1=0.0_r8! coefficient to convert mass flux to mb/d
    !real(r8) clrsky(pcols)         ! not used ! temporary used in random overlap calc
    rhpert=0.0_r8! the specified perturbation to rh
    deepcu=0.0_r8! deep convection cloud fraction
    shallowcu=0.0_r8! shallow convection cloud fraction
    kp1=0
    numkcld=0! number of levels in which to allow clouds


    !==================================================================================
    ! PHILOSOPHY OF PRESENT IMPLEMENTATION
    !
    ! There are three co-existing cloud types: convective, inversion related low-level
    ! stratocumulus, and layered cloud (based on relative humidity).  Layered and 
    ! stratocumulus clouds do not compete with convective cloud for which one creates 
    ! the most cloud.  They contribute collectively to the total grid-box average cloud 
    ! amount.  This is reflected in the way in which the total cloud amount is evaluated 
    ! (a sum as opposed to a logical "or" operation)
    !
    !==================================================================================
    ! set defaults for rhu00
    rhu00(:,:) = 2.0_r8
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01_r8 

    !
    ! Evaluate potential temperature and relative humidity
    !
    CALL aqsat(temp    ,pmid    ,es      ,qs      ,pcols   , &
         ncol    ,pver    ,1       ,pver    )
    DO k=1,pver
       DO i=1,ncol
          theta(i,k)  = temp(i,k)*(pnot(i)/pmid(i,k))**cappa
          rh(i,k)     = q(i,k)/qs(i,k)*(1.0_r8+real(dindex,kind=r8)*rhpert)
          !
          !  record relhum, rh itself will later be modified related with concld
          !
          relhum(i,k) = rh(i,k)
          cloud(i,k)  = 0.0_r8
          cldst(i,k)  = 0.0_r8
          concld(i,k) = 0.0_r8
       END DO
    END DO
    !
    ! Initialize other temporary variables
    !
    ierror = 0
    DO i=1,ncol
       ! Adjust thetas(i) in the presence of non-zero ocean heights.
       ! This reduces the temperature for positive heights according to a standard lapse rate.
       IF(ocnfrac(i).GT.0.01_r8) thetas(i)  = &
            ( sst(i) - lapse * phis(i) / gravit) * (pnot(i)/ps(i))**cappa
       IF(ocnfrac(i).GT.0.01_r8.AND.sst(i).LT.260.0_r8) ierror = i
       clc(i) = 0.0_r8
    END DO
    coef1 = gravit*864.0_r8    ! conversion to millibars/day

    IF (ierror > 0) THEN
       WRITE(6,*) 'COLDSST: encountered in cldfrc:', ierror,ocnfrac(ierror),sst(ierror)
    ENDIF

    DO k=1,pver-1
       DO i=1,ncol
          rpdeli(i,k) = 1.0_r8/(pmid(i,k+1) - pmid(i,k))
       END DO
    END DO

    !
    ! Estimate of local convective cloud cover based on convective mass flux
    ! Modify local large-scale relative humidity to account for presence of 
    ! convective cloud when evaluating relative humidity based layered cloud amount
    !
    DO k=1,pver
       DO i=1,ncol
          concld(i,k) = 0.0_r8
       END DO
    END DO
    !
    ! cloud mass flux in SI units of kg/m2/s; should produce typical numbers of 20%
    ! shallow and deep convective cloudiness are evaluated separately (since processes
    ! are evaluated separately) and summed
    !   
    !#ifndef PERGRO
    IF(PERGRO)  THEN
       DO k=1,pver-1
          DO i=1,ncol
             shallowcu = MAX(0.0_r8,MIN(sh1*LOG(1.0_r8+sh2* cmfmc2(i,k+1)),0.30_r8))
             deepcu    = MAX(0.0_r8,MIN(dp1*LOG(1.0_r8+dp2*(cmfmc (i,k+1)-cmfmc2(i,k+1))),0.60_r8))
             concld(i,k) = MIN(shallowcu + deepcu,0.80_r8)
             rh(i,k) = (rh(i,k) - concld(i,k))/(1.0_r8 - concld(i,k))
          END DO
       END DO
    END IF
!!!!!#endif
    !==================================================================================
    !
    !          ****** Compute layer cloudiness ******
    !
    !====================================================================
    ! Begin the evaluation of layered cloud amount based on (modified) RH 
    !====================================================================
    !
    numkcld = pver
    DO k=2,numkcld
       kp1 = MIN(k + 1,pver)
       DO i=1,ncol
          !
          cldbnd(i) = pmid(i,k).GE.pretop
          !WRITE(0,*)i,k, pmid(i,k),premit 
          !premib = 750.e2_r8   ! bottom pressure bound of middle cloud
          IF ( pmid(i,k).GE.premib ) THEN
             !==============================================================
             ! This is the low cloud (below premib) block
             !==============================================================
             ! enhance low cloud activation over land with no snow cover
             IF (land(i) .AND. (snowh(i) <= 0.000001_r8)) THEN
                rhlim = rhminl - 0.10_r8
             ELSE
                rhlim = rhminl
             ENDIF
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = MIN(0.999_r8,(MAX(rhdif,0.0_r8))**2)
             ! premit = 250.e2_r8  ! top of area defined to be mid-level cloud
          ELSE IF ( pmid(i,k).LT.premit ) THEN
             !==============================================================
             ! This is the high cloud (above premit) block
             !==============================================================
             !
             rhlim = rhminh
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = MIN(0.999_r8,(MAX(rhdif,0.0_r8))**2)
          ELSE
             !==============================================================
             ! This is the middle cloud block
             !==============================================================
             !
             !       linear rh threshold transition between thresholds for low & high cloud
             !
             rhwght = (premib-(MAX(pmid(i,k),premit)))/(premib-premit)

             IF (land(i) .AND. (snowh(i) <= 0.000001_r8)) THEN
                rhlim = rhminh*rhwght + (rhminl - 0.10_r8)*(1.0_r8-rhwght)
             ELSE
                rhlim = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
             ENDIF
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = MIN(0.999_r8,(MAX(rhdif,0.0_r8))**2)
          END IF
          !==================================================================================
          ! WE NEED TO DOCUMENT THE PURPOSE OF THIS TYPE OF CODE (ASSOCIATED WITH 2ND CALL)
          !==================================================================================
          !      !
          !      ! save rhlim to rhu00, it handles well by itself for low/high cloud
          !      !
          rhu00(i,k)=rhlim
          !==================================================================================
       END DO
       !
       ! Final evaluation of layered cloud fraction
       !
    END DO
    !
    ! Add in the marine strat
    ! MARINE STRATUS SHOULD BE A SPECIAL CASE OF LAYERED CLOUD
    ! CLOUD CURRENTLY CONTAINS LAYERED CLOUD DETERMINED BY RH CRITERIA
    ! TAKE THE MAXIMUM OF THE DIAGNOSED LAYERED CLOUD OR STRATOCUMULUS
    !
    !===================================================================================
    !
    !  SOME OBSERVATIONS ABOUT THE FOLLOWING SECTION OF CODE (missed in earlier look)
    !  K700 IS SET AS A CONSTANT BASED ON HYBRID COORDINATE: IT DOES NOT DEPEND ON 
    !  LOCAL PRESSURE; THERE IS NO PRESSURE RAMP => LOOKS LEVEL DEPENDENT AND 
    !  DISCONTINUOUS IN SPACE (I.E., STRATUS WILL END SUDDENLY WITH NO TRANSITION)
    !
    !  IT APPEARS THAT STRAT IS EVALUATED ACCORDING TO KLEIN AND HARTMANN; HOWEVER,
    !  THE ACTUAL STRATUS AMOUNT (CLDST) APPEARS TO DEPEND DIRECTLY ON THE RH BELOW
    !  THE STRONGEST PART OF THE LOW LEVEL INVERSION.  
    !
    !==================================================================================
    !
    ! Find most stable level below 750 mb for evaluating stratus regimes
    !
    DO i=1,ncol
       ! Nothing triggers unless a stability greater than this minimum threshold is found
       dthdpmn(i) = -0.125_r8
       kdthdp(i) = 0
    END DO
    !
    DO k=2,pver
       DO i=1,ncol
          IF (pmid(i,k) >= premib .AND. ocnfrac(i).GT. 0.01_r8) THEN
             ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
             dthdp = 100.0_r8*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
             IF (dthdp < dthdpmn(i)) THEN
                dthdpmn(i) = dthdp
                kdthdp(i) = k     ! index of interface of max inversion
             END IF
          END IF
       END DO
    END DO

    ! Also check between the bottom layer and the surface
    ! Only perform this check if the criteria were not met above

    DO i = 1,ncol
       IF ( kdthdp(i) .EQ. 0 .AND. ocnfrac(i).GT.0.01_r8) THEN
          dthdp = 100.0_r8 * (thetas(i) - theta(i,pver)) / (ps(i)-pmid(i,pver))
          IF (dthdp < dthdpmn(i)) THEN
             dthdpmn(i) = dthdp
             kdthdp(i) = pver     ! index of interface of max inversion
          ENDIF
       ENDIF
    ENDDO

    DO i=1,ncol
       IF (kdthdp(i) /= 0) THEN
          k = kdthdp(i)
          kp1 = MIN(k+1,pver)
          ! Note: strat will be zero unless ocnfrac > 0.01_r8
          strat = MIN(1.0_r8,MAX(0.0_r8, ocnfrac(i) * ((theta(i,k700)-thetas(i))*.057_r8-.5573_r8) ) )
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
          cldst(i,k) = MIN(strat,MAX(rh(i,k),rh(i,kp1)))
       END IF
    END DO
    !
    ! AGGREGATE CLOUD CONTRIBUTIONS (cldst should be zero everywhere except at level kdthdp(i))
    !
    DO k=1,pver
       DO i=1,ncol
          !
          !       which is greater; standard layered cloud amount or stratocumulus diagnosis
          !
          cloud(i,k) = MAX(cloud(i,k),cldst(i,k))
          !
          !       add in the contributions of convective cloud (determined separately and accounted
          !       for by modifications to the large-scale relative humidity.
          !
          cloud(i,k) = MIN(cloud(i,k)+concld(i,k), 1.0_r8)
       END DO
    END DO

    !
    RETURN
  END SUBROUTINE cldfrc
  !end module cloud_fraction


  SUBROUTINE gestbl(epsil   , &
       latvap  ,latice  ,rh2o    ,cpair      )
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
    !REAL(r8), INTENT(in) :: tmn           ! Minimum temperature entry in es lookup table
    !REAL(r8), INTENT(in) :: tmx           ! Maximum temperature entry in es lookup table
    REAL(r8), INTENT(in) :: epsil         ! Ratio of h2o to dry air molecular weights
    REAL(r8), INTENT(in) :: latvap        ! Latent heat of vaporization
    REAL(r8), INTENT(in) :: latice        ! Latent heat of fusion
    REAL(r8), INTENT(in) :: rh2o          ! Gas constant for water vapor
    REAL(r8), INTENT(in) :: cpair         ! Specific heat of dry air
    !REAL(r8), INTENT(in) :: tmeltx        ! Melting point of water (K)
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(r8)  epsqs
    REAL(r8) t             ! Temperature
    INTEGER n          ! Increment counter
    INTEGER lentbl     ! Calculated length of lookup table
    INTEGER itype      ! Ice phase: 0 -> no ice phase
    !            1 -> ice phase, no transitiong
    !           -x -> ice phase, x degree transition
    !LOGICAL ip         ! Ice phase logical flag
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
          itype = -INT(ttrice,kind=i4)
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
    INTEGER , INTENT(in) :: len       ! vector length
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
       piln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , pcols, pver, pverp          )
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
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    INTEGER, INTENT(in) :: pverp
    REAL(r8), INTENT(in) :: piln (pcols,pverp)   ! Log interface pressures
    !REAL(r8), INTENT(in) :: pmln (pcols,pver)    ! Log midpoint pressures
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
          !
          !  S = CpT + HGT
          !  T = (S - HGT)/cp
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
       piln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     ,pcols, pver, pverp)

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
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    INTEGER, INTENT(in) :: pverp
    REAL(r8), INTENT(in) :: piln (pcols,pverp)   ! Log interface pressures
    !REAL(r8), INTENT(in) :: pmln (pcols,pver)    ! Log midpoint pressures
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
    CHARACTER(len=*), INTENT(in) :: msg    ! string to be printed

    !IF (PRESENT (msg)) THEN
       WRITE(6,*)'ENDRUN:', msg
    !ELSE
    !   WRITE(6,*)'ENDRUN: called without a message string'
    !END IF
    STOP
  END SUBROUTINE endrun


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


  SUBROUTINE DestroyMicro_Hack()
    IMPLICIT NONE

  END SUBROUTINE DestroyMicro_Hack

END MODULE Micro_Hack



!PROGRAM MAIN
!
! USE Micro_Hack
! INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
! INTEGER , PARAMETER:: kMax=28
! INTEGER, PARAMETER :: iMax=1  
! INTEGER, PARAMETER :: jMax=1
! INTEGER, PARAMETER :: ibMax=1
! INTEGER, PARAMETER :: jbMax=1
! INTEGER, PARAMETER :: ppcnst=3
! INTEGER :: si    (kmax+1)
! INTEGER :: sl    (kmax)
! REAL(r8) :: dtime                      ! timestep
! REAL(r8) :: state_ps     (ibMax)    !(ibMax)     ! surface pressure(Pa)
! REAL(r8) :: state_phis   (ibMax)    !(ibMax)     ! surface geopotential
! REAL(r8) :: state_t           (ibMax,kMax)!(ibMax,kMax)! temperature (K)
! REAL(r8) :: state_qv     (ibMax,kMax)!(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
! REAL(r8) :: state_ql     (ibMax,kMax)!(ibMax,kMax,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
! REAL(r8) :: state_qi     (ibMax,kMax)!(ibMax,kMax,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
! REAL(r8) :: state_omega  (ibMax,kMax)!(ibMax,kMax)! vertical pressure velocity (Pa/s) 
! REAL(r8) :: icefrac           (ibMax)        ! sea ice fraction (fraction)
! REAL(r8) :: landfrac     (ibMax)        ! land fraction (fraction)
! REAL(r8) :: ocnfrac           (ibMax)        ! ocean fraction (fraction)
! REAL(r8) :: landm           (ibMax)        ! land fraction ramped over water
! REAL(r8) :: snowh           (ibMax)        ! Snow depth over land, water equivalent (m)
! REAL(r8) :: state_dlf    (ibMax,kMax)    ! detrained water from ZM
! REAL(r8) :: rliq           (ibMax)           ! vertical integral of liquid not yet in q(ixcldliq)
! REAL(r8) :: state_cmfmc  (ibMax,kMax+1)   ! convective mass flux--m sub c
! REAL(r8) :: state_cmfmc2 (ibMax,kMax+1)   ! shallow convective mass flux--m sub c
! REAL(r8) :: state_concld (ibMax,kMax)    ! convective cloud cover
! REAL(r8) :: ts           (ibMax)         ! surface temperature
! REAL(r8) :: sst           (ibMax)         ! sea surface temperature
! REAL(r8) :: state_zdu    (ibMax,kMax)       ! detrainment rate from deep convection
! REAL(r8) :: prec_str     (ibMax)           ! [Total] sfc flux of precip from stratiform (m/s) 
! REAL(r8) :: snow_str     (ibMax)  ! [Total] sfc flux of snow from stratiform   (m/s)
! REAL(r8) :: prec_sed     (ibMax)  ! surface flux of total cloud water from sedimentation
! REAL(r8) :: snow_sed     (ibMax)  ! surface flux of cloud ice from sedimentation
! REAL(r8) :: prec_pcw     (ibMax)  ! sfc flux of precip from Micro_Hack(m/s)
! REAL(r8) :: snow_pcw     (ibMax)  ! sfc flux of snow from Micro_Hack (m/s)
! REAL(r8) :: state_qcwat  (1:ibMax,1:kMax)  ! cloud water old q
! REAL(r8) :: state_tcwat  (1:ibMax,1:kMax)  ! cloud water old temperature
! REAL(r8) :: state_lcwat  (1:ibMax,1:kMax)  ! cloud liquid water old q
! REAL(r8) :: state_cld    (1:ibMax,1:kMax)  ! cloud fraction
! REAL(r8) :: state_qme    (1:ibMax,1:kMax) 
! REAL(r8) :: state_prain  (1:ibMax,1:kMax) 
! REAL(r8) :: state_nevapr (1:ibMax,1:kMax) 
! REAL(r8) :: state_rel    (1:ibMax,1:kMax)   ! liquid effective drop radius (microns)
! REAL(r8) :: state_rei    (1:ibMax,1:kMax)   ! ice effective drop size (microns)
! !qmin=1e-12
! !k700 = 3

! CALL  Init_Micro_Hack(kMax,jMax,ibMax,jbMax,ppcnst,si,sl)


! CALL RunMicro_Hack(&
!      ibMax       , &! INTEGER , INTENT(in) :: ibMax                 ! number of columns (max)
!      kMax        , &! INTEGER , INTENT(in) :: kMax                  ! number of vertical levels
!      kMax+1      , &! INTEGER , INTENT(in) :: kMax+1                 ! number of vertical levels + 1
!      ppcnst      , &! INTEGER , INTENT(in) :: ppcnst          ! number of constituent
!      dtime       , &! REAL(r8), INTENT(in)  :: dtime                   ! timestep
!      state_ps    , &! REAL(r8), INTENT(in)  :: state_ps     (ibMax)    !(ibMax)     ! surface pressure(Pa)
!      state_phis  , &! REAL(r8), INTENT(in)  :: state_phis   (ibMax)    !(ibMax)     ! surface geopotential
!      state_t     , &! REAL(r8), INTENT(in)  :: state_t      (ibMax,kMax)!(ibMax,kMax)! temperature (K)
!      state_qv    , &! REAL(r8), INTENT(in)  :: state_qv     (ibMax,kMax)!(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
!      state_ql    , &! REAL(r8), INTENT(in)  :: state_ql     (ibMax,kMax)!(ibMax,kMax,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
!      state_qi    , &! REAL(r8), INTENT(in)  :: state_qi     (ibMax,kMax)!(ibMax,kMax,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
!      state_omega , &! REAL(r8), INTENT(in)  :: state_omega  (ibMax,kMax)!(ibMax,kMax)! vertical pressure velocity (Pa/s) 
!      icefrac     , &! REAL(r8), INTENT(in)  :: icefrac      (ibMax)             ! sea ice fraction (fraction)
!      landfrac    , &! REAL(r8), INTENT(in)  :: landfrac     (ibMax)             ! land fraction (fraction)
!      ocnfrac     , &! REAL(r8), INTENT(in)  :: ocnfrac      (ibMax)             ! ocean fraction (fraction)
!      landm       , &! REAL(r8), INTENT(in)  :: landm              (ibMax)             ! land fraction ramped over water
!      snowh       , &! REAL(r8), INTENT(in)  :: snowh              (ibMax)             ! Snow depth over land, water equivalent (m)
!      state_dlf   , &! REAL(r8), INTENT(in)  :: state_dlf    (ibMax,kMax)    ! detrained water from ZM
!      rliq        , &! REAL(r8), INTENT(in)  :: rliq         (ibMax)         ! vertical integral of liquid not yet in q(ixcldliq)
!      state_cmfmc , &! REAL(r8), INTENT(in)  :: state_cmfmc  (ibMax,kMax+1)   ! convective mass flux--m sub c
!      state_cmfmc2, &! REAL(r8), INTENT(in)  :: state_cmfmc2 (ibMax,kMax+1)   ! shallow convective mass flux--m sub c
!      state_concld, &! REAL(r8), INTENT(out) :: state_concld (ibMax,kMax)    ! convective cloud cover
!      ts          , &! REAL(r8), INTENT(in)  :: ts              (ibMax)              ! surface temperature
!      sst         , &! REAL(r8), INTENT(in)  :: sst              (ibMax)              ! sea surface temperature
!      state_zdu   , &! REAL(r8), INTENT(in)  :: state_zdu    (ibMax,kMax)           ! detrainment rate from deep convection
!      prec_str    , &! REAL(r8), INTENT(out)  :: prec_str    (ibMax)         ! [Total] sfc flux of precip from stratiform (m/s) 
!      snow_str    , &! REAL(r8), INTENT(out)  :: snow_str    (ibMax)  ! [Total] sfc flux of snow from stratiform   (m/s)
!      prec_sed    , &! REAL(r8), INTENT(out)  :: prec_sed    (ibMax)  ! surface flux of total cloud water from sedimentation
!      snow_sed    , &! REAL(r8), INTENT(out)  :: snow_sed    (ibMax)  ! surface flux of cloud ice from sedimentation
!      prec_pcw    , &! REAL(r8), INTENT(out)  :: prec_pcw    (ibMax)  ! sfc flux of precip from Micro_Hack(m/s)
!      snow_pcw    , &! REAL(r8), INTENT(out)  :: snow_pcw    (ibMax)  ! sfc flux of snow from Micro_Hack (m/s)
!      state_qcwat , &! REAL(r8), INTENT(inout) :: state_qcwat(1:ibMax,1:kMax)  ! cloud water old q
!      state_tcwat , &! REAL(r8), INTENT(inout) :: state_tcwat(1:ibMax,1:kMax)  ! cloud water old temperature
!      state_lcwat , &! REAL(r8), INTENT(inout) :: state_lcwat(1:ibMax,1:kMax)  ! cloud liquid water old q
!      state_cld   , &! REAL(r8), INTENT(inout) :: state_cld  (1:ibMax,1:kMax)  ! cloud fraction
!      state_qme   , &! REAL(r8), INTENT(inout) :: state_qme  (1:ibMax,1:kMax) 
!      state_prain , &! REAL(r8), INTENT(inout) :: state_prain(1:ibMax,1:kMax) 
!      state_nevapr, &! REAL(r8), INTENT(inout) :: state_nevapr(1:ibMax,1:kMax) 
!      state_rel   , &! REAL(r8), INTENT(inout) :: state_rel  (1:ibMax,1:kMax)   ! liquid effective drop radius (microns)
!      state_rei     )! REAL(r8), INTENT(inout) :: state_rei  (1:ibMax,1:kMax)   ! ice effective drop size (microns)
!
!END PROGRAM MAIN
