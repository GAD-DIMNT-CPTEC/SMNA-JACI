MODULE GwddSchemeCAM
    IMPLICIT NONE
  SAVE

  PRIVATE                         ! Make default type private to the module

  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers
  !----------------------------------------------------------------------------
  ! physical constants (all data public)
  !----------------------------------------------------------------------------
  REAL(R8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
  REAL(R8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
  REAL(R8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_R8       ! molecular weight dry air ~ kg/kmole
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(R8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(R8),PARAMETER :: SHR_CONST_G      = 9.80616_R8      ! acceleration of gravity ~ m/s^2
  REAL(R8),PARAMETER :: SHR_CONST_CPWV   = 1.810e3_R8      ! specific heat of water vap ~ J/kg/K
  REAL(R8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg

  REAL(r8), PUBLIC, PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: rair = shr_const_rdair    ! Gas constant for dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: cpwv  = shr_const_cpwv
  REAL(r8), PUBLIC, PARAMETER :: zvir = SHR_CONST_RWV/rair - 1          ! rh2o/rair - 1

  ! Constants for Earth
  REAL(r8), PUBLIC, PARAMETER :: gravit = shr_const_g      ! gravitational acceleration
  REAL(r8), PUBLIC, PARAMETER :: rga = 1.0_r8/gravit           ! reciprocal of gravit



  !
  ! PUBLIC: interfaces
  !
    public gw_inti                  ! Initialization
    public gw_intr                  ! interface to actual parameterization
  !
  ! PRIVATE: Rest of the data and interfaces are private to this module
  !
  INTEGER, PARAMETER :: pgwv = 0  ! number of waves allowed

  INTEGER ::  kbotoro      ! interface of gwd source
  INTEGER :: ktopbg, ktoporo      ! top interface of gwd region

  REAL(r8), ALLOCATABLE ::alpha(:)    ! alpha(0:pver)       ! newtonian cooling coefficients
  REAL(r8) :: c(-pgwv:pgwv)       ! list of wave phase speeds
  REAL(r8) :: dback               ! background diffusivity
  REAL(r8) :: effkwv              ! effective wavenumber (fcrit2*kwv)
  REAL(r8) :: effgw_oro           ! tendency efficiency for orographic gw
  REAL(r8) :: effgw_spec=.125_r8  ! tendency efficiency for internal gw
  REAL(r8) :: fracldv             ! fraction of stress deposited in low level region
  REAL(r8) :: g                   ! acceleration of gravity
  REAL(r8) :: kwv                 ! effective horizontal wave number
  REAL(r8) :: lzmax               ! maximum vertical wavelength at source

  REAL(r8) :: mxasym              ! max asymmetry between tau(c) and tau(-c)
  REAL(r8) :: mxrange             ! max range of tau for all c
  REAL(r8) :: n2min               ! min value of bouyancy frequency
  REAL(r8) :: fcrit2              ! critical froude number
  REAL(r8) :: oroko2              ! 1/2 * horizontal wavenumber
  REAL(r8) :: orohmin             ! min surface displacment height for orographic waves
  REAL(r8) :: orovmin             ! min wind speed for orographic waves
  REAL(r8) :: r                   ! gas constant for dry air
  REAL(r8) :: rog                 ! r / g
  REAL(r8) :: taubgnd             ! background source strength (/tauscal)
  REAL(r8) :: taumin              ! minimum (nonzero) stress
  REAL(r8) :: tauscal             ! scale factor for background stress source
  REAL(r8) :: tndmax              ! maximum wind tendency
  REAL(r8) :: umcfac              ! factor to limit tendency to prevent reversing u-c
  REAL(r8) :: ubmc2mn             ! min (u-c)**2
  REAL(r8) :: zldvcon             ! constant for determining zldv from tau0
  !REAL(KIND=r8),ALLOCATABLE :: si (:)
  !REAL(KIND=r8),ALLOCATABLE :: sl (:)
  !REAL(KIND=r8),ALLOCATABLE :: delsig(:)

CONTAINS
  !===============================================================================
  ! call gw_inti (   ,hypi    )

  SUBROUTINE gw_inti (kMax)
    !-----------------------------------------------------------------------
    ! Time independent initialization for multiple gravity wave parameterization.
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(IN  ) :: kMax
    !REAL(KIND=r8), INTENT(IN   ) :: si_in    (kMax+1)
    !REAL(KIND=r8), INTENT(IN   ) :: sl_in    (kMax)
    !REAL(KIND=r8), INTENT(IN   ) :: delsig_in(kMax)

    !REAL(r8) :: hypi(kMax+1)          ! reference interface pressures

    !---------------------------Local storage-------------------------------
    INTEGER :: k

    !ALLOCATE(si    (kmax+1))
    !ALLOCATE(sl    (kmax  ))
    !ALLOCATE(delsig(kmax  ))
    ALLOCATE(alpha(0:kMax)) 

    !DO k=1,kMax+1
    !   hypi(k) =  si_in    (kMax+2-k)
    !END DO
    !si    =si_in
    !sl    =sl_in
    !delsig=delsig_in
    !-----------------------------------------------------------------------
    ! Copy model constants

    g      = gravit
    r      = rair

    ! Set MGWD constants
    kwv    = 6.28e-5_r8          ! 100 km wave length
    dback  = 0.05_r8             ! background diffusivity
    fcrit2 = 0.5_r8              ! critical froude number squared
    tauscal= 0.001_r8            ! scale factor for background stress
    taubgnd= 6.4_r8              ! background stress amplitude

    zldvcon= 10.0_r8              ! constant for determining zldv
    lzmax  = 7.0E3_r8             ! maximum vertical wavelength at source (m)
    IF ( dycore_is ('LR') ) THEN
       effgw_oro = 0.125_r8
       fracldv= 0.0_r8
       ! Restore these to work with turulent mountain stress
              effgw_oro = 1.0_r8
              fracldv= 0.7_r8           ! fraction of tau0 diverged in low level region
    ELSE
       !effgw_oro = 0.125_r8
       effgw_oro = 1.0_r8
       fracldv= 0.0_r8
       ! Restore these to work with turulent mountain stress
       effgw_oro = 0.125_r8
       !       fracldv= 0.7_r8           ! fraction of tau0 diverged in low level region
    ENDIF
    ! Set phase speeds 
    DO k = -pgwv, pgwv
       c(k)   = 10.0_r8 * k       ! 0, +/- 10, +/- 20, ... m/s
    END DO

    !    if (masterproc) then
    !       write(6,*) ' '
    !       write(6,*) 'GW_INTI: pgwv = ', pgwv
    !       write(6,*) 'GW_INTI: c(l) = ', c
    !       write(6,*) ' '
    !    end if

    ! Set radiative damping times
    DO k = 0, kMax
       alpha(k) = 1.0e-6_r8       ! about 10 days.
    END DO

    ! Min and max values to keep things reasonable
    mxasym =  0.1_r8              ! max factor of 10 from |tau(c)| to |tau(-c)|
    mxrange=  0.001_r8            ! factor of 100 from max to min |tau(c)|
    n2min  =  1.e-8_r8            ! min value of Brunt-Vaisalla freq squared
    orohmin= 10.0_r8              ! min surface displacement for orographic wave drag
    orovmin=  2.0_r8              ! min wind speed for orographic wave drag
    taumin =  1.e-10_r8           ! min stress considered > 0
    tndmax = 500.0_r8 / 86400.0_r8    ! max permitted tendency (500 m/s/day)
    umcfac = 0.5_r8              ! max permitted reduction in u-c
    ubmc2mn= 0.01_r8             ! min value of (u-c)^2

    ! Determine other derived constants
    oroko2 = 0.5_r8 * kwv
    effkwv = fcrit2 * kwv
    rog    = r/g

    ! Determine the bounds of the background and orographic stress regions
    ktopbg  = 0
    kbotoro = kMax
    !DO k = 0, kMax
    !   IF (hypi(k+1) < 0.1_r8) kbotbg  = k    ! spectrum source at 100 mb 100/1000 = 0.1
    !END DO
    ktoporo = 0
    RETURN
  END  SUBROUTINE gw_inti


  !===============================================================================

  SUBROUTINE gw_intr ( &
      pcols    , &!INTEGER, INTENT(IN   ) :: pcols
      pver     , &!INTEGER, INTENT(IN   ) :: pver
      pverp    , &!INTEGER, INTENT(IN   ) :: pverp
      prsi ,prsl  ,phii ,phil    ,&
      gt       , &!REAL(r8), INTENT(in) :: gt (pcols,pver)  
      gq       , &!REAL(r8), INTENT(in) :: gq (pcols,pver)  
      gu       , &!REAL(r8), INTENT(in) :: gu (pcols,pver)  
      gv       , &!REAL(r8), INTENT(in) :: gv (pcols,pver)  
      sgh      , &!REAL(r8), INTENT(in) :: sgh(pcols)                ! standard deviation of orography
      dt       , &!REAL(r8), INTENT(in) :: dt                        ! time step
      landfrac , &!REAL(r8), INTENT(in) :: landfrac(pcols)        ! Land fraction
      rlat     , &!REAL(r8), INTENT(in) :: rlat(pcols)             ! latitude in radians for columns
      psi      , &
      chug     , &
      chvg     , &
      chtg       )
    !-----------------------------------------------------------------------
    ! Interface for multiple gravity wave drag parameterization.
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    INTEGER      , INTENT(IN   ) :: pcols
    INTEGER      , INTENT(IN   ) :: pver
    INTEGER      , INTENT(IN   ) :: pverp
    REAL(KIND=r8), INTENT(in   ) :: prsi   (pcols,pver+1)  !     prsi     - real, pressure at layer interfaces [Pa]
    REAL(KIND=r8), INTENT(in   ) :: prsl   (pcols,pver)    !     prsl     - real, mean layer presure [Pa]
    REAL(KIND=r8), INTENT(in   ) :: phii   (pcols,pver+1) !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8), INTENT(in   ) :: phil   (pcols,pver)   !===>  PHIL(K) INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(r8)     , INTENT(in   ) :: gt (pcols,pver)  
    REAL(r8)     , INTENT(in   ) :: gq (pcols,pver)  
    REAL(r8), INTENT(in) :: gu (pcols,pver)  
    REAL(r8), INTENT(in) :: gv (pcols,pver)  
    REAL(r8), INTENT(in) :: sgh(pcols)            ! standard deviation of orography
    REAL(r8), INTENT(in) :: dt                    ! time step
    REAL(r8), INTENT(in) :: landfrac(pcols)        ! Land fraction
    REAL(r8), INTENT(in) :: rlat(pcols)             ! latitude in radians for columns
    REAL(r8), INTENT(in) :: psi(pcols)        
    REAL(r8), INTENT(inOUT) :: chug  (pcols,pver)  
    REAL(r8), INTENT(inOUT) :: chvg  (pcols,pver)  
    REAL(r8), INTENT(inOUT) :: chtg  (pcols,pver)  
    
    REAL(r8) :: ttgw  (pcols,pver)                  ! temperature tendency

    !    type(physics_state), intent(in) :: state      ! physics state structure
    !    type(physics_ptend), intent(inout):: ptend    ! parameterization tendency structure
    REAL(r8) :: state_t (pcols,pver)  
    REAL(r8) :: state_t2 (pcols,pver)  
    REAL(r8) :: state_q (pcols,pver)  
    REAL(r8) :: state_u (pcols,pver)  
    REAL(r8) :: state_v (pcols,pver)  
    REAL(r8) :: state_pmid (pcols,pver)  
    REAL(r8) :: state_pint(pcols,pver+1)  
    REAL(r8) :: state_lnpint(pcols,pver+1)  
    REAL(r8) :: state_pdel(pcols,pver)  
    REAL(r8) :: state_rpdel(pcols,pver)
    REAL(r8) :: state_lnpmid(pcols,pver)
    REAL(r8) :: ptend_u(pcols,pver)
    REAL(r8) :: ptend_v(pcols,pver)
    REAL(r8) :: ptend_s(pcols,pver)

    REAL(r8) :: state_zi    (1:pcols,1:pver+1)   
    REAL(r8) :: state_zm    (1:pcols,1:pver)
    REAL(r8) :: state_s     (1:pcols,1:pver)
    REAL(r8) :: state_phis  (1:pcols)
    !---------------------------Local storage-------------------------------
    INTEGER :: ncol                               ! number of atmospheric columns

    INTEGER :: i,k                                ! loop indexes
    INTEGER :: kldv(pcols)                        ! top interface of low level stress divergence region
    INTEGER :: kldvmn                             ! min value of kldv
    INTEGER :: ksrc(pcols)                        ! index of top interface of source region
    INTEGER :: ksrcmn                             ! min value of ksrc

    REAL(r8) :: utgw(pcols,pver)                  ! zonal wind tendency
    REAL(r8) :: vtgw(pcols,pver)                  ! meridional wind tendency

    REAL(r8) :: ni(pcols,0:pver)                  ! interface Brunt-Vaisalla frequency
    REAL(r8) :: nm(pcols,pver)                    ! midpoint Brunt-Vaisalla frequency
    REAL(r8) :: rdpldv(pcols)                     ! 1/dp across low level divergence region
    REAL(r8) :: rhoi(pcols,0:pver)                ! interface density
    REAL(r8) :: tau(pcols,-pgwv:pgwv,0:pver)      ! wave Reynolds stress
    REAL(r8) :: tau0x(pcols)                      ! c=0 sfc. stress (zonal)
    REAL(r8) :: tau0y(pcols)                      ! c=0 sfc. stress (meridional)
    REAL(r8) :: ti(pcols,0:pver)                  ! interface temperature
    REAL(r8) :: ubi(pcols,0:pver)                 ! projection of wind at interfaces
    REAL(r8) :: ubm(pcols,pver)                   ! projection of wind at midpoints
    REAL(r8) :: xv(pcols)                         ! unit vectors of source wind (x)
    REAL(r8) :: yv(pcols)                         ! unit vectors of source wind (y)
    REAL(r8) :: hkl(pcols)
    REAL(r8) :: hypi(pcols,pver+1)          ! reference interface pressures
    INTEGER  :: kbot(pcols) 
    INTEGER  :: kbotorov(pcols)
    REAL(r8) :: hkk (pcols) 
    REAL(r8) :: tvfac  
    LOGICAL  :: fvdyn
    !-----------------------------------------------------------------------------
    ncol  = pcols
    DO i=1,pcols
       state_phis(i) = psi(i)
       !state_pint       (i,pver+1) = gps(i)*si(1)
    END DO
    DO k=pver+1,1,-1
       DO i=1,pcols
         ! state_pint    (i,k) = MAX(si(pver+2-k)*gps(i) ,0.0001_r8)
          state_pint    (i,k) = prsi(i,pver+2-k)
          state_zi      (i,k) = phii(i,pver+2-k)
       END DO
    END DO

       DO i=1,pcols
          state_pint    (i,1) = prsl(i,pver)
       END DO

    DO k=1,pver
       DO i=1,pcols
          state_t (i,pver+1-k) =  gt (i,k)
          state_q (i,pver+1-k) =  gq (i,k)
          state_u (i,pver+1-k) =  gu (i,k)
          state_v (i,pver+1-k) =  gv (i,k)
          !state_pmid(i,pver+1-k) = sl(k)*gps (i)
          state_pmid(i,pver+1-k) = prsl(i,k)
          state_zm  (i,pver+1-k) = phil(i,k)  

       END DO
    END DO
    DO k=1,pver+1
       DO i=1,pcols
          !hypi(i,k) =  si_in    (kMax+2-k)
          hypi(i,k) =  prsi(i,pver+2-k)/prsi(i,1)
       END DO
    END DO
    DO k = 0, pver
       DO i=1,pcols
          !IF (hypi(k+1)      < 0.1_r8) kbotbg  = k    ! spectrum source at 100 mb 100/1000 = 0.1
          IF (hypi    (i,k+1) < 0.1_r8) THEN
             ! ! spectrum source at 100 mb 100/1000 = 0.1
             kbot(i)= k
          END IF
          kbotorov(i)=pver
       END DO
    END DO


    DO k=1,pver
       DO i=1,pcols    
          state_pdel    (i,k) =        MAX (state_pint(i,k+1) - state_pint(i,k) ,0.00000000005_r8)
          state_rpdel   (i,k) = 1.0_r8/MAX((state_pint(i,k+1) - state_pint(i,k)),0.00000000005_r8)
          state_lnpmid  (i,k) = LOG(state_pmid(i,k))        
       END DO
    END DO
    DO k=1,pver+1
       DO i=1,pcols
          state_lnpint(i,k) =  LOG(state_pint  (i,k))
       END DO
    END DO

    ! Derive new temperature and geopotential fields

    CALL geopotential_t(                                 &
         state_lnpint(1:pcols,1:pver+1)   , state_pint (1:pcols,1:pver+1)   , &
         state_pmid  (1:pcols,1:pver)     , state_pdel  (1:pcols,1:pver)   , state_rpdel(1:pcols,1:pver)   , &
         state_t     (1:pcols,1:pver)     , state_q     (1:pcols,1:pver)   , rair   , gravit , zvir   ,          &
         state_zi    (1:pcols,1:pver+1)   , state_zm    (1:pcols,1:pver)   , ncol   ,pcols, pver, pverp)

    fvdyn = dycore_is ('LR')
    DO k = pver, 1, -1
       ! First set hydrostatic elements consistent with dynamics
       IF (fvdyn) THEN
          DO i = 1,ncol
             hkl(i) = state_lnpmid(i,k+1) - state_lnpmid(i,k)
             hkk(i) = 1.0_r8 - state_pint (i,k)* hkl(i)* state_rpdel(i,k)
          END DO
       ELSE
          DO i = 1,ncol
             hkl(i) = state_pdel(i,k)/state_pmid  (i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF
       ! Now compute s
       DO i = 1,ncol
          tvfac   = 1.0_r8 + zvir * state_q(i,k) 
          state_s(i,k) =  (state_t(i,k)* cpair) + (state_t(i,k) * tvfac * rair*hkk(i))  +  &
                                ( state_phis(i) + gravit*state_zi(i,k+1))           
       END DO
    END DO

    ! Profiles of background state variables
    CALL gw_prof( &
         ncol                         , & !integer , intent(in ) :: ncol                      ! number of atmospheric columns
         pcols                        , & !integer , intent(in ) :: pcols
         pver                         , & !integer , intent(in ) :: pver
         state_t    (1:pcols,1:pver)  , & !real(r8), intent(in ) :: t(pcols,pver)             ! midpoint temperatures
         state_pmid (1:pcols,1:pver)  , & !real(r8), intent(in ) :: pm(pcols,pver)             ! midpoint pressures
         state_pint (1:pcols,1:pver+1), & !real(r8), intent(in ) :: pi(pcols,0:pver)      ! interface pressures
         rhoi       (1:pcols,0:pver)  , & !real(r8), intent(out) :: rhoi(pcols,0:pver)   ! interface density
         ni         (1:pcols,0:pver)  , & !real(r8), intent(out) :: ni(pcols,0:pver)     ! interface Brunt-Vaisalla frequency
         ti         (1:pcols,0:pver)  , & !real(r8), intent(out) :: ti(pcols,0:pver)     ! interface temperature
         nm         (1:pcols,1:pver)    ) !real(r8), intent(out) :: nm(pcols,pver)       ! midpoint Brunt-Vaisalla frequency

    !-----------------------------------------------------------------------------
    ! Non-orographic backgound gravity wave spectrum
    !-----------------------------------------------------------------------------
    IF (pgwv >0) THEN

       ! Determine the wave source for a background spectrum at ~100 mb

       CALL gw_bgnd ( &
            ncol                                    , &!integer, intent(in) :: ncol! number of atmospheric columns
            pcols                                   , &!integer, intent(in) :: pcols
            pver                                    , &!integer, intent(in) :: pver
            state_u     (1:pcols,1:pver)            , &!real(r8), intent(in) :: u(pcols,pver)! midpoint zonal wind
            state_v     (1:pcols,1:pver)            , &!real(r8), intent(in) :: v(pcols,pver)! midpoint meridional wind
            rlat        (1:pcols)                   , &!real(r8), intent(in) :: rlat(pcols)! latitude in radians for columns
            kldv        (1:pcols)                   , &!integer, intent(out) :: kldv(pcols)! top interface of low level stress 
            kldvmn                                  , &!integer, intent(out) :: kldvmn! min value of kldv
            ksrc        (1:pcols)                   , &!integer, intent(out) :: ksrc(pcols)! index of top interface of source region
            ksrcmn                                  , &!integer, intent(out) :: ksrcmn! min value of ksrc
            tau         (1:pcols,-pgwv:pgwv,0:pver) , &!real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
            ubi         (1:pcols,0:pver)            , &!real(r8), intent(out) :: ubi(pcols,0:pver)! projection of wind at interfaces
            ubm         (1:pcols,1:pver)            , &!real(r8), intent(out) :: ubm(pcols,pver)! projection of wind at midpoints
            xv          (1:pcols)                   , &!real(r8), intent(out) :: xv(pcols)! unit vectors of source wind (x)
            yv          (1:pcols)                   , &!real(r8), intent(out) :: yv(pcols)! unit vectors of source wind (y)
            PGWV                                    , &!integer, intent(in) :: ngwv! number of gravity waves to use
            kbot       (1:pcols)                       )!integer, intent(in) :: kbot! index of bottom (source) interface

       ! Solve for the drag profile

       CALL gw_drag_prof (                            &
            ncol                                   , &!integer, intent(in) :: ncol! number of atmospheric columns
            pcols                                  , &!integer, intent(in) :: pcols 
            pver                                   , &!integer, intent(in) :: pver        
            PGWV                                   , &!integer, intent(in) :: ngwv! number of gravity waves to use
            kbot        (1:pcols)                  , &!integer, intent(in) :: kbot! index of bottom (source) interface
            ktopbg                                 , &!integer, intent(in) :: ktop! index of top interface of gwd region
            state_t     (1:pcols,1:pver)           , &!real(r8), intent(in) :: t(pcols,pver)! midpoint temperatures
            state_pint  (1:pcols,1:pver+1)         , &!real(r8), intent(in) :: pi(pcols,0:pver)! interface pressures
            state_pdel  (1:pcols,1:pver)           , &!real(r8), intent(in) :: dpm(pcols,pver)! midpoint delta p (pi(k)-pi(k-1))
            state_rpdel (1:pcols,1:pver)           , &!real(r8), intent(in) :: rdpm(pcols,pver)! 1. / (pi(k)-pi(k-1))
            state_lnpint(1:pcols,1:pver+1)         , &!real(r8), intent(in) :: piln(pcols,0:pver)! ln(interface pressures)
            rhoi        (1:pcols,0:pver)           , &!real(r8), intent(in) :: rhoi(pcols,0:pver)! interface density
            ni          (1:pcols,0:pver)           , &!real(r8), intent(in) :: ni(pcols,0:pver)! interface Brunt-Vaisalla frequency
            ti          (1:pcols,0:pver)           , &!real(r8), intent(in) :: ti(pcols,0:pver)! interface temperature
            nm          (1:pcols,1:pver)           , &!real(r8), intent(in) :: nm(pcols,pver)! midpoint Brunt-Vaisalla frequency
            dt                                     , &!real(r8), intent(in) :: dt! time step
            kldv        (1:pcols)                  , &!integer, intent(in) :: kldv(pcols)! top interface of low level stress  
            kldvmn                                 , &!integer, intent(in) :: kldvmn! min value of kldv
            ksrc        (1:pcols)                  , &!integer, intent(in) :: ksrc(pcols)! index of top interface of source region
            ksrcmn                                 , &!integer, intent(in) :: ksrcmn! min value of ksrc
            rdpldv      (1:pcols)                  , &!real(r8), intent(in) :: rdpldv(pcols)! 1/dp across low level divergence region
            tau         (1:pcols,-pgwv:pgwv,0:pver), &!real(r8), intent(inout) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
            ubi         (1:pcols,0:pver)           , &!real(r8), intent(in) :: ubi(pcols,0:pver)! projection of wind at interfaces
            ubm         (1:pcols,1:pver)           , &!real(r8), intent(in) :: ubm(pcols,pver)! projection of wind at midpoints
            xv          (1:pcols)                  , &!real(r8), intent(in) :: xv(pcols)! unit vectors of source wind (x)
            yv          (1:pcols)                  , &!real(r8), intent(in) :: yv(pcols)! unit vectors of source wind (y)
            effgw_spec                             , &!real(r8), intent(in) :: effgw! tendency efficiency
            utgw        (1:pcols,1:pver)           , &!real(r8), intent(out) :: ut(pcols,pver)! zonal wind tendency
            vtgw        (1:pcols,1:pver)           , &!real(r8), intent(out) :: vt(pcols,pver)! meridional wind tendency
            tau0x       (1:pcols)                  , &!real(r8), intent(out) :: tau0x(pcols)! c=0 sfc. stress (zonal)
            tau0y       (1:pcols)                    )!real(r8), intent(out) :: tau0y(pcols)! c=0 sfc. stress (meridional)

       ! Add the momentum tendencies to the output tendency arrays
       DO k = 1, pver
          DO i = 1, ncol
             ptend_u(i,k) = utgw(i,k)
             ptend_v(i,k) = vtgw(i,k)
          END DO
       END DO

       ! zero net tendencies if no spectrum computed

    ELSE
       ptend_u = 0.0_r8
       ptend_v = 0.0_r8
    END IF
    !-----------------------------------------------------------------------------
    ! Orographic stationary gravity wave
    !-----------------------------------------------------------------------------

    ! Determine the orographic wave source

    CALL gw_oro ( &
         ncol                                  , &!integer, intent(in) :: ncol! number of atmospheric columns
         pcols                                 , &!integer, intent(in) :: pcols
         pver                                  , &!integer, intent(in) :: pver
         state_u    (1:pcols,1:pver)           , &!real(r8), intent(in) :: u(pcols,pver)! midpoint zonal wind
         state_v    (1:pcols,1:pver)           , &!real(r8), intent(in) :: v(pcols,pver)! midpoint meridional wind
         state_t    (1:pcols,1:pver)           , &!real(r8), intent(in) :: t(pcols,pver)! midpoint temperatures
         sgh        (1:pcols)                  , &!real(r8), intent(in) :: sgh(pcols)! standard deviation of orography
         state_pmid (1:pcols,1:pver)           , &!real(r8), intent(in) :: pm(pcols,pver)! midpoint pressures
         state_pint (1:pcols,1:pver+1)         , &!real(r8), intent(in) :: pi(pcols,0:pver)! interface pressures
         state_pdel (1:pcols,1:pver)           , &!real(r8), intent(in) :: dpm(pcols,pver)! midpoint delta p (pi(k)-pi(k-1))
         state_zm   (1:pcols,1:pver)           , &!real(r8), intent(in) :: zm(pcols,pver)! midpoint heights
         nm         (1:pcols,1:pver)           , &!real(r8), intent(in) :: nm(pcols,pver)! midpoint Brunt-Vaisalla frequency
         kldv       (1:pcols)                  , &!integer, intent(out) :: kldv(pcols)! top interface of low level stress div region
         kldvmn                                , &!integer, intent(out) :: kldvmn! min value of kldv
         ksrc       (1:pcols)                  , &!integer, intent(out) :: ksrc(pcols)! index of top interface of source region
         ksrcmn                                , &!integer, intent(out) :: ksrcmn! min value of ksrc
         rdpldv     (1:pcols)                  , &!real(r8), intent(out) :: rdpldv(pcols)! 1/dp across low level divergence region
         tau        (1:pcols,-pgwv:pgwv,0:pver), &!real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
         ubi        (1:pcols,0:pver)           , &!real(r8), intent(out) :: ubi(pcols,0:pver)! projection of wind at interfaces
         ubm        (1:pcols,1:pver)           , &!real(r8), intent(out) :: ubm(pcols,pver)! projection of wind at midpoints
         xv         (1:pcols)                  , &!real(r8), intent(out) :: xv(pcols)! unit vectors of source wind (x)
         yv         (1:pcols)                    )!real(r8), intent(out) :: yv(pcols)! unit vectors of source wind (y)

    ! Solve for the drag profile

    CALL gw_drag_prof ( &
         ncol                                 , &!integer, intent(in) :: ncol                   ! number of atmospheric columns
         pcols                                , &!integer, intent(in) :: pcols 
         pver                                 , &!integer, intent(in) :: pver      
         0                                    , &!integer, intent(in) :: ngwv                        ! number of gravity waves to use
         kbotorov   (1:pcols)                            , &!integer, intent(in) :: kbot                     ! index of bottom (source) interface
         ktoporo                              , &!integer, intent(in) :: ktop                     ! index of top interface of gwd region
         state_t     (1:pcols,1:pver)              , &!real(r8), intent(in) :: t(pcols,pver)        ! midpoint temperatures
         state_pint  (1:pcols,1:pver+1)              , &!real(r8), intent(in) :: pi(pcols,0:pver)     ! interface pressures
         state_pdel  (1:pcols,1:pver)              , &!real(r8), intent(in) :: dpm(pcols,pver)      ! midpoint delta p (pi(k)-pi(k-1))
         state_rpdel (1:pcols,1:pver)              , &!real(r8), intent(in) :: rdpm(pcols,pver)     ! 1. / (pi(k)-pi(k-1))
         state_lnpint(1:pcols,1:pver+1)              , &!real(r8), intent(in) :: piln(pcols,0:pver)   ! ln(interface pressures)
         rhoi             (1:pcols,0:pver)              , &!real(r8), intent(in) :: rhoi(pcols,0:pver)   ! interface density
         ni             (1:pcols,0:pver)              , &!real(r8), intent(in) :: ni(pcols,0:pver)     ! interface Brunt-Vaisalla frequency
         ti             (1:pcols,0:pver)              , &!real(r8), intent(in) :: ti(pcols,0:pver)     ! interface temperature
         nm             (1:pcols,1:pver)              , &!real(r8), intent(in) :: nm(pcols,pver)       ! midpoint Brunt-Vaisalla frequency
         dt                                      , &!real(r8), intent(in) :: dt                        ! time step
         kldv           (1:pcols)                      , &!integer, intent(in) :: kldv(pcols)               ! top interface of low level stress  
         kldvmn                                      , &!integer, intent(in) :: kldvmn                    ! min value of kldv
         ksrc           (1:pcols)                      , &!integer, intent(in) :: ksrc(pcols)           ! index of top interface of source region
         ksrcmn                                      , &!integer, intent(in) :: ksrcmn                    ! min value of ksrc
         rdpldv           (1:pcols)                      , &!real(r8), intent(in) :: rdpldv(pcols)         ! 1/dp across low level divergence region
         tau           (1:pcols,-pgwv:pgwv,0:pver), &!real(r8), intent(inout) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
         ubi           (1:pcols,0:pver)              , &!real(r8), intent(in) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
         ubm           (1:pcols,1:pver)              , &!real(r8), intent(in) :: ubm(pcols,pver)      ! projection of wind at midpoints
         xv           (1:pcols)                      , &!real(r8), intent(in) :: xv(pcols)               ! unit vectors of source wind (x)
         yv           (1:pcols)                      , &!real(r8), intent(in) :: yv(pcols)               ! unit vectors of source wind (y)
         effgw_oro                                , &!real(r8), intent(in) :: effgw                ! tendency efficiency
         utgw           (1:pcols,1:pver)              , &!real(r8), intent(out) :: ut(pcols,pver)      ! zonal wind tendency
         vtgw           (1:pcols,1:pver)              , &!real(r8), intent(out) :: vt(pcols,pver)      ! meridional wind tendency
         tau0x           (1:pcols)                      , &!real(r8), intent(out) :: tau0x(pcols)        ! c=0 sfc. stress (zonal)
         tau0y           (1:pcols)                        )!real(r8), intent(out) :: tau0y(pcols)        ! c=0 sfc. stress (meridional)

    ! Add the orographic tendencies to the spectrum tendencies
    ! Compute the temperature tendency from energy conservation (includes spectrum).
    DO k = 1, pver
       DO i = 1, ncol
          ptend_u(i,k) = ptend_u(i,k) + utgw(i,k) * landfrac(i)
          ptend_v(i,k) = ptend_v(i,k) + vtgw(i,k) * landfrac(i)
          ptend_s(i,k) = -(ptend_u(i,k) * (state_u(i,k) + ptend_u(i,k)*0.5_r8*dt) &
                         + ptend_v(i,k) * (state_v(i,k) + ptend_v(i,k)*0.5_r8*dt))
          ttgw(i,k) = ptend_s(i,k) / cpair
          state_s(i,k)=state_s(i,k)+ptend_s(i,k)*dt
       END DO
    END DO
    CALL geopotential_dse(                                                                    &
            state_lnpint(1:pcols,1:pver+1), state_pint (1:pcols,1:pver+1)  , &
            state_pmid  (1:pcols,1:pver), state_pdel  (1:pcols,1:pver) , state_rpdel(1:pcols,1:pver)  , &
            state_s     (1:pcols,1:pver), state_q     (1:ncol,1:pver), state_phis (1:pcols) , rair  , &
            gravit      , cpair        ,zvir        , &
            state_t2(1:pcols,1:pver)     , state_zi(1:pcols,1:pver+1)    , state_zm(1:pcols,1:pver)       ,&
            ncol      ,pcols, pver, pver+1          )

    ! Set flags for nonzero tendencies, q not yet affected by gwd
    !    ptend%name  = "vertical diffusion"
    !    ptend%lq(:) = .FALSE.
    !    ptend%ls    = .TRUE.
    !    ptend%lu    = .TRUE.
    !    ptend%lv    = .TRUE.
    !
    !   call physics_update (state, tend, ptend, dt)
    DO k = 1, pver
       DO i = 1, ncol
         chug(i,pver+1-k)= ptend_u(i,k)
         chvg(i,pver+1-k)= ptend_v(i,k)
         chtg(i,pver+1-k)= (state_t2(i,k)-state_t(i,k))/dt
       END DO
    END DO         
    RETURN
  END  SUBROUTINE gw_intr

  !===============================================================================
  SUBROUTINE gw_prof ( &
       ncol       , &!integer , intent(in) :: ncol                   ! number of atmospheric columns
       pcols      , &!integer , intent(in) :: pcols
       pver       , &!integer , intent(in) :: pver
       t          , &!real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
       pm         , &!real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
       pi         , &!real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
       rhoi       , &!real(r8), intent(out) :: rhoi(pcols,0:pver)   ! interface density
       ni         , &!real(r8), intent(out) :: ni(pcols,0:pver)     ! interface Brunt-Vaisalla frequency
       ti         , &!real(r8), intent(out) :: ti(pcols,0:pver)     ! interface temperature
       nm           )!real(r8), intent(out) :: nm(pcols,pver)       ! midpoint Brunt-Vaisalla frequency
    !-----------------------------------------------------------------------
    ! Compute profiles of background state quantities for the multiple
    ! gravity wave drag parameterization.
    ! 
    ! The parameterization is assumed to operate only where water vapor 
    ! concentrations are negligible in determining the density.
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(in) :: ncol                   ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    REAL(r8), INTENT(in) :: t(pcols,pver)         ! midpoint temperatures
    REAL(r8), INTENT(in) :: pm(pcols,pver)        ! midpoint pressures
    REAL(r8), INTENT(in) :: pi(pcols,0:pver)      ! interface pressures

    REAL(r8), INTENT(out) :: rhoi(pcols,0:pver)   ! interface density
    REAL(r8), INTENT(out) :: ni(pcols,0:pver)     ! interface Brunt-Vaisalla frequency
    REAL(r8), INTENT(out) :: ti(pcols,0:pver)     ! interface temperature
    REAL(r8), INTENT(out) :: nm(pcols,pver)       ! midpoint Brunt-Vaisalla frequency

    !---------------------------Local storage-------------------------------
    INTEGER :: i,k                                ! loop indexes

    REAL(r8) :: dtdp
    REAL(r8) :: n2                                ! Brunt-Vaisalla frequency squared

    !-----------------------------------------------------------------------------
    ! Determine the interface densities and Brunt-Vaisala frequencies.
    !-----------------------------------------------------------------------------

    ! The top interface values are calculated assuming an isothermal atmosphere 
    ! above the top level.
    k = 0
    DO i = 1, ncol
       ti(i,k) = t(i,k+1)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k) = SQRT (g*g / (cpair*ti(i,k)))
    END DO

    ! Interior points use centered differences
    DO k = 1, pver-1
       DO i = 1, ncol
          ti(i,k) = 0.5_r8 * (t(i,k) + t(i,k+1))
          rhoi(i,k) = pi(i,k) / (r*ti(i,k))
          dtdp = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
          n2 = g*g/ti(i,k) * (1.0_r8/cpair - rhoi(i,k)*dtdp)
          ni(i,k) = SQRT (MAX (n2min, n2))
       END DO
    END DO

    ! Bottom interface uses bottom level temperature, density; next interface
    ! B-V frequency.
    k = pver
    DO i = 1, ncol
       ti(i,k) = t(i,k)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k) = ni(i,k-1)
    END DO

    !-----------------------------------------------------------------------------
    ! Determine the midpoint Brunt-Vaisala frequencies.
    !-----------------------------------------------------------------------------
    DO k=1,pver
       DO i=1,ncol
          nm(i,k) = 0.5_r8 * (ni(i,k-1) + ni(i,k))
       END DO
    END DO

    RETURN
  END SUBROUTINE gw_prof

  !===============================================================================

  SUBROUTINE gw_oro ( &
       ncol     , &!integer, intent(in) :: ncol                   ! number of atmospheric columns
       pcols    , &!integer, intent(in) :: pcols
       pver     , &!integer, intent(in) :: pver
       u        , &!real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
       v        , &!real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
       t        , &!real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
       sgh      , &!real(r8), intent(in) :: sgh(pcols)                 ! standard deviation of orography
       pm       , &!real(r8), intent(in) :: pm(pcols,pver)         ! midpoint pressures
       pi       , &!real(r8), intent(in) :: pi(pcols,0:pver)         ! interface pressures
       dpm      , &!real(r8), intent(in) :: dpm(pcols,pver)         ! midpoint delta p (pi(k)-pi(k-1))
       zm       , &!real(r8), intent(in) :: zm(pcols,pver)         ! midpoint heights
       nm       , &!real(r8), intent(in) :: nm(pcols,pver)         ! midpoint Brunt-Vaisalla frequency
       kldv     , &!integer, intent(out) :: kldv(pcols)           ! top interface of low level stress div region
       kldvmn   , &!integer, intent(out) :: kldvmn                 ! min value of kldv
       ksrc     , &!integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
       ksrcmn   , &!integer, intent(out) :: ksrcmn                 ! min value of ksrc
       rdpldv   , &!real(r8), intent(out) :: rdpldv(pcols)         ! 1/dp across low level divergence region
       tau      , &!real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
       ubi      , &!real(r8), intent(out) :: ubi(pcols,0:pver)         ! projection of wind at interfaces
       ubm      , &!real(r8), intent(out) :: ubm(pcols,pver)         ! projection of wind at midpoints
       xv       , &!real(r8), intent(out) :: xv(pcols)                 ! unit vectors of source wind (x)
       yv         )!real(r8), intent(out) :: yv(pcols)                 ! unit vectors of source wind (y)
    !-----------------------------------------------------------------------
    ! Orographic source for multiple gravity wave drag parameterization.
    ! 
    ! The stress is returned for a single wave with c=0, over orography.
    ! For points where the orographic variance is small (including ocean),
    ! the returned stress is zero. 
    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(in) :: ncol                   ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    REAL(r8), INTENT(in) :: u(pcols,pver)         ! midpoint zonal wind
    REAL(r8), INTENT(in) :: v(pcols,pver)         ! midpoint meridional wind
    REAL(r8), INTENT(in) :: t(pcols,pver)         ! midpoint temperatures
    REAL(r8), INTENT(in) :: sgh(pcols)            ! standard deviation of orography
    REAL(r8), INTENT(in) :: pm(pcols,pver)        ! midpoint pressures
    REAL(r8), INTENT(in) :: pi(pcols,0:pver)      ! interface pressures
    REAL(r8), INTENT(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    REAL(r8), INTENT(in) :: zm(pcols,pver)        ! midpoint heights
    REAL(r8), INTENT(in) :: nm(pcols,pver)        ! midpoint Brunt-Vaisalla frequency

    INTEGER, INTENT(out) :: kldv(pcols)           ! top interface of low level stress div region
    INTEGER, INTENT(out) :: kldvmn                ! min value of kldv
    INTEGER, INTENT(out) :: ksrc(pcols)           ! index of top interface of source region
    INTEGER, INTENT(out) :: ksrcmn                ! min value of ksrc

    REAL(r8), INTENT(out) :: rdpldv(pcols)        ! 1/dp across low level divergence region
    REAL(r8), INTENT(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
    REAL(r8), INTENT(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    REAL(r8), INTENT(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    REAL(r8), INTENT(out) :: xv(pcols)            ! unit vectors of source wind (x)
    REAL(r8), INTENT(out) :: yv(pcols)            ! unit vectors of source wind (y)

    !---------------------------Local storage-------------------------------
    INTEGER :: i,k                                ! loop indexes

    REAL(r8) :: pil                               ! don't we have pi somewhere?
    REAL(r8) :: lzsrc                             ! vertical wavelength at source
    REAL(r8) :: hdsp(pcols)                       ! surface streamline displacment height (2*sgh)
    REAL(r8) :: sghmax                            ! max orographic sdv to use
    REAL(r8) :: tauoro(pcols)                     ! c=0 stress from orography
    REAL(r8) :: zldv(pcols)                       ! top of the low level stress divergence region
    REAL(r8) :: nsrc(pcols)                       ! b-f frequency averaged over source region
    REAL(r8) :: psrc(pcols)                       ! interface pressure at top of source region
    REAL(r8) :: rsrc(pcols)                       ! density averaged over source region
    REAL(r8) :: usrc(pcols)                       ! u wind averaged over source region
    REAL(r8) :: vsrc(pcols)                       ! v wind averaged over source region

    !---------------------------------------------------------------------------
    ! Average the basic state variables for the wave source over the depth of
    ! the orographic standard deviation. Here we assume that the apropiate
    ! values of wind, stability, etc. for determining the wave source are 
    ! averages over the depth of the atmosphere pentrated by the typical mountain.
    ! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
    !---------------------------------------------------------------------------
    k = pver
    DO i = 1, ncol
       ksrc(i) = k-1
       psrc(i) = pi(i,k-1)
       rsrc(i) = pm(i,k)/(r*t(i,k)) * dpm(i,k)
       usrc(i) = u(i,k) * dpm(i,k)
       vsrc(i) = v(i,k) * dpm(i,k)
       nsrc(i) = nm(i,k)* dpm(i,k)
       hdsp(i) = 2.0_r8 * sgh(i)
    END DO
    DO k = pver-1, pver/2, -1
       DO i = 1, ncol
          IF (hdsp(i) > SQRT(zm(i,k)*zm(i,k+1))) THEN
             ksrc(i) = k-1
             psrc(i) = pi(i,k-1)
             rsrc(i) = rsrc(i) + pm(i,k) / (r*t(i,k))* dpm(i,k)
             usrc(i) = usrc(i) + u(i,k) * dpm(i,k)
             vsrc(i) = vsrc(i) + v(i,k) * dpm(i,k)
             nsrc(i) = nsrc(i) + nm(i,k)* dpm(i,k)
          END IF
       END DO
    END DO
    DO i = 1, ncol
       rsrc(i) = rsrc(i) / (pi(i,pver) - psrc(i))
       usrc(i) = usrc(i) / (pi(i,pver) - psrc(i))
       vsrc(i) = vsrc(i) / (pi(i,pver) - psrc(i))
       nsrc(i) = nsrc(i) / (pi(i,pver) - psrc(i))

       !#if ( defined SCAM )
       !! needed the following fix when winds are identically 0
       !! orig ->  ubi(i,pver) = sqrt (usrc(i)**2 + vsrc(i)**2)
       !
       !       ubi(i,pver) = max(sqrt (usrc(i)**2 + vsrc(i)**2),orovmin)
       !
       !#else

       ubi(i,pver) = SQRT (usrc(i)**2 + vsrc(i)**2)

       !#endif
       xv(i) = usrc(i) / ubi(i,pver)
       yv(i) = vsrc(i) / ubi(i,pver)
    END DO

    ! Project the local wind at midpoints onto the source wind.
    DO k = 1, pver
       DO i = 1, ncol
          ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
       END DO
    END DO

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    DO i = 1, ncol
       ubi(i,0) = ubm(i,1)
    END DO
    DO k = 1, pver-1
       DO i = 1, ncol
          ubi(i,k) = 0.5_r8 * (ubm(i,k) + ubm(i,k+1))
       END DO
    END DO

    !---------------------------------------------------------------------------
    ! Determine the depth of the low level stress divergence region, as
    ! the max of the source region depth and 1/2 the vertical wavelength at the
    ! source. 
    !---------------------------------------------------------------------------
    pil = ACOS(-1.0_r8)
    DO i = 1, ncol
       lzsrc   = MIN(2.0_r8 * pil * usrc(i) / nsrc(i), lzmax)
       zldv(i) = MAX(hdsp(i), 0.5_r8 * lzsrc)
    END DO

    ! find the index of the top of the low level divergence region

    kldv(:) = pver-1
    DO k = pver-1, pver/2, -1
       DO i = 1, ncol
          IF (zldv(i) .GT. SQRT(zm(i,k)*zm(i,k+1))) THEN
             kldv(i)  = k-1
          END IF
       END DO
    END DO

    ! Determine the orographic c=0 source term following McFarlane (1987).
    ! Set the source top interface index to pver, if the orographic term is zero.
    DO i = 1, ncol
       IF ((ubi(i,pver) .GT. orovmin) .AND. (hdsp(i) .GT. orohmin)) THEN
          sghmax = fcrit2 * (ubi(i,pver) / nsrc(i))**2
          tauoro(i) = oroko2 * MIN(hdsp(i)**2, sghmax) * rsrc(i) * nsrc(i) * ubi(i,pver)
       ELSE
          tauoro(i) = 0.0_r8
          ksrc(i) = pver
          kldv(i) = pver
       END IF
    END DO

    ! Set the phase speeds and wave numbers in the direction of the source wind.
    ! Set the source stress magnitude (positive only, note that the sign of the 
    ! stress is the same as (c-u).
    DO i = 1, ncol
       tau(i,0,pver) = tauoro(i)
    END DO

    ! Determine the min value of kldv and ksrc for limiting later loops
    ! and the pressure at the top interface of the low level stress divergence
    ! region.

    ksrcmn = pver
    kldvmn = pver
    DO i = 1, ncol
       ksrcmn = MIN(ksrcmn, ksrc(i))
       kldvmn = MIN(kldvmn, kldv(i))
       IF (kldv(i) .NE. pver) THEN
          rdpldv(i) = 1.0_r8 / (pi(i,kldv(i)) - pi(i,pver))
       END IF
    END DO
    IF (fracldv .LE. 0.0_r8) kldvmn = pver

    RETURN
  END SUBROUTINE gw_oro


  !===============================================================================
  SUBROUTINE gw_bgnd ( &
       ncol     , &!integer, intent(in) :: ncol                   ! number of atmospheric columns
       pcols    , &!integer, intent(in) :: pcols
       pver     , &!integer, intent(in) :: pver
       u        , &!real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
       v        , &!real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
       rlat     , &!real(r8), intent(in) :: rlat(pcols)           ! latitude in radians for columns
       kldv     , &!integer, intent(out) :: kldv(pcols)              ! top interface of low level stress divergence region
       kldvmn   , &!integer, intent(out) :: kldvmn                ! min value of kldv
       ksrc     , &!integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
       ksrcmn   , &!integer, intent(out) :: ksrcmn                ! min value of ksrc
       tau      , &!real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
       ubi      , &!real(r8), intent(out) :: ubi(pcols,0:pver)         ! projection of wind at interfaces
       ubm      , &!real(r8), intent(out) :: ubm(pcols,pver)         ! projection of wind at midpoints
       xv       , &!real(r8), intent(out) :: xv(pcols)                 ! unit vectors of source wind (x)
       yv       , &!real(r8), intent(out) :: yv(pcols)                 ! unit vectors of source wind (y)
       ngwv     , &!integer, intent(in) :: ngwv                   ! number of gravity waves to use
       kbot      ) !integer, intent(in) :: kbot                   ! index of bottom (source) interface
    !-----------------------------------------------------------------------
    ! Driver for multiple gravity wave drag parameterization.
    ! 
    ! The parameterization is assumed to operate only where water vapor 
    ! concentrations are negligible in determining the density.
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(in) :: ncol                   ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols
    INTEGER, INTENT(in) :: pver
    INTEGER, INTENT(in) :: kbot (pcols)           ! index of bottom (source) interface
    INTEGER, INTENT(in) :: ngwv                   ! number of gravity waves to use

    REAL(r8), INTENT(in) :: u(pcols,pver)         ! midpoint zonal wind
    REAL(r8), INTENT(in) :: v(pcols,pver)         ! midpoint meridional wind
    REAL(r8), INTENT(in) :: rlat(pcols)           ! latitude in radians for columns

    INTEGER, INTENT(out) :: kldv(pcols)           ! top interface of low level stress divergence region
    INTEGER, INTENT(out) :: kldvmn                ! min value of kldv
    INTEGER, INTENT(out) :: ksrc(pcols)           ! index of top interface of source region
    INTEGER, INTENT(out) :: ksrcmn                ! min value of ksrc

    REAL(r8), INTENT(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
    REAL(r8), INTENT(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    REAL(r8), INTENT(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    REAL(r8), INTENT(out) :: xv(pcols)            ! unit vectors of source wind (x)
    REAL(r8), INTENT(out) :: yv(pcols)            ! unit vectors of source wind (y)

    !---------------------------Local storage-------------------------------
    INTEGER :: i,k,l                              ! loop indexes

    REAL(r8) :: tauback(pcols)                    ! background stress at c=0
    REAL(r8) :: usrc(pcols)                       ! u wind averaged over source region
    REAL(r8) :: vsrc(pcols)                       ! v wind averaged over source region
    REAL(r8) :: al0                               ! Used in lat dependence of GW spec. 
    REAL(r8) :: dlat0                             ! Used in lat dependence of GW spec.
    REAL(r8) :: flat_gw                           ! The actual lat dependence of GW spec.
    REAL(r8) :: pi_g                              ! 3.14........

    !---------------------------------------------------------------------------
    ! Determine the source layer wind and unit vectors, then project winds.
    !---------------------------------------------------------------------------

    ! Just use the source level interface values for the source
    ! wind speed and direction (unit vector).
    DO i = 1, ncol
       k = kbot(i)
       ksrc(i) = k
       kldv(i) = k
       usrc(i) = 0.5_r8*(u(i,k+1)+u(i,k))
       vsrc(i) = 0.5_r8*(v(i,k+1)+v(i,k))
       ubi(i,kbot(i)) = SQRT (usrc(i)**2 + vsrc(i)**2)
       xv(i) = usrc(i) / ubi(i,k)
       yv(i) = vsrc(i) / ubi(i,k)
    END DO

    ! Project the local wind at midpoints onto the source wind.
    DO k = 1, pver
       DO i = 1, ncol
          IF(k>=1  .and. k<=kbot(i))THEN
             ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
          END IF
       END DO
    END DO

    ! Compute the interface wind projection by averaging the midpoint winds.
    ! Use the top level wind at the top interface.
    DO i = 1, ncol
       ubi(i,0) = ubm(i,1)
    END DO

    DO k = 1, pver
       DO i = 1, ncol
          IF(k>=1  .and. k<=kbot(i)-1)THEN
             ubi(i,k) = 0.5_r8 * (ubm(i,k) + ubm(i,k+1))
          END IF
       END DO
    END DO

    !-----------------------------------------------------------------------
    ! Gravity wave sources
    !-----------------------------------------------------------------------

    ! Determine the background stress at c=0
    DO i=1,ncol
       tauback(i) = taubgnd * tauscal
    ENDDO

    !        Include dependence on latitude:
    !         The lat function was obtained by RR Garcia as 
    !        currently used in his 2D model
    !       [Added by F. Sassi on May 30, 1997]

    pi_g=4.0_r8*ATAN(1.0_r8)  ! 3.14
    al0=40.0_r8*pi_g/180.0_r8  ! 0.69
    dlat0=10.0_r8*pi_g/180.0_r8! 0.175
    !
    DO i=1,ncol
       flat_gw= 0.5_r8*(1.0_r8+TANH((rlat(i)-al0)/dlat0)) + 0.5_r8*(1.0_r8+TANH(-(rlat(i)+al0)/dlat0)) 
       flat_gw=MAX(flat_gw,0.2_r8)
       tauback(i)=tauback(i)*flat_gw
    ENDDO

    ! Set the phase speeds and wave numbers in the direction of the source wind.
    ! Set the source stress magnitude (positive only, note that the sign of the 
    ! stress is the same as (c-u).

    DO l = 1, ngwv
       DO i = 1, ncol
          tau(i, l,kbot(i)) = tauback(i) * EXP(-(c(l)/30.0_r8)**2)
          tau(i,-l,kbot(i)) = tau(i, l,kbot(i))
       END DO
    END DO
    DO i = 1, ncol
       tau(i,0,kbot(i)) = tauback(i)
    END DO

    ! Determine the min value of kldv and ksrc for limiting later loops
    ! and the pressure at the top interface of the low level stress divergence
    ! region.

    ksrcmn = pver
    kldvmn = pver

    RETURN
  END  SUBROUTINE gw_bgnd

  !===============================================================================
  SUBROUTINE gw_drag_prof ( &
       ncol     , &!integer, intent(in) :: ncol                   ! number of atmospheric columns
       pcols    , &!integer, intent(in) :: pcols 
       pver     , &!integer, intent(in) :: pver      
       ngwv     , &!integer, intent(in) :: ngwv                   ! number of gravity waves to use
       kbot     , &!integer, intent(in) :: kbot                   ! index of bottom (source) interface
       ktop     , &!integer, intent(in) :: ktop                   ! index of top interface of gwd region
       t        , &!real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
       pi       , &!real(r8), intent(in) :: pi(pcols,0:pver)         ! interface pressures
       dpm      , &!real(r8), intent(in) :: dpm(pcols,pver)         ! midpoint delta p (pi(k)-pi(k-1))
       rdpm     , &!real(r8), intent(in) :: rdpm(pcols,pver)         ! 1. / (pi(k)-pi(k-1))
       piln     , &!real(r8), intent(in) :: piln(pcols,0:pver)         ! ln(interface pressures)
       rhoi     , &!real(r8), intent(in) :: rhoi(pcols,0:pver)         ! interface density
       ni       , &!real(r8), intent(in) :: ni(pcols,0:pver)         ! interface Brunt-Vaisalla frequency
       ti       , &!real(r8), intent(in) :: ti(pcols,0:pver)         ! interface temperature
       nm       , &!real(r8), intent(in) :: nm(pcols,pver)         ! midpoint Brunt-Vaisalla frequency
       dt       , &!real(r8), intent(in) :: dt                    ! time step
       kldv     , &!integer, intent(in) :: kldv(pcols)                 ! top interface of low level stress  divergence region
       kldvmn   , &!integer, intent(in) :: kldvmn                     ! min value of kldv
       ksrc     , &!integer, intent(in) :: ksrc(pcols)             ! index of top interface of source region
       ksrcmn   , &!integer, intent(in) :: ksrcmn                     ! min value of ksrc
       rdpldv   , &!real(r8), intent(in) :: rdpldv(pcols)         ! 1/dp across low level divergence region
       tau      , &!real(r8), intent(inout) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
       ubi      , &!real(r8), intent(in) :: ubi(pcols,0:pver)         ! projection of wind at interfaces
       ubm      , &!real(r8), intent(in) :: ubm(pcols,pver)         ! projection of wind at midpoints
       xv       , &!real(r8), intent(in) :: xv(pcols)                 ! unit vectors of source wind (x)
       yv       , &!real(r8), intent(in) :: yv(pcols)                 ! unit vectors of source wind (y)
       effgw    , &!real(r8), intent(in) :: effgw                 ! tendency efficiency
       ut       , &!real(r8), intent(out) :: ut(pcols,pver)         ! zonal wind tendency
       vt       , &!real(r8), intent(out) :: vt(pcols,pver)         ! meridional wind tendency
       tau0x    , &!real(r8), intent(out) :: tau0x(pcols)         ! c=0 sfc. stress (zonal)
       tau0y      )!real(r8), intent(out) :: tau0y(pcols)         ! c=0 sfc. stress (meridional)
    !-----------------------------------------------------------------------
    ! Solve for the drag profile from the multiple gravity wave drag
    ! parameterization.
    ! 1. scan up from the wave source to determine the stress profile
    ! 2. scan down the stress profile to determine the tendencies
    !     => apply bounds to the tendency
    !          a. from wkb solution
    !          b. from computational stability constraints
    !     => adjust stress on interface below to reflect actual bounded tendency
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    INTEGER, INTENT(in) :: ncol                   ! number of atmospheric columns
    INTEGER, INTENT(in) :: pcols 
    INTEGER, INTENT(in) :: pver      
    INTEGER, INTENT(in) :: kbot (pcols)                   ! index of bottom (source) interface
    INTEGER, INTENT(in) :: ktop                   ! index of top interface of gwd region
    INTEGER, INTENT(in) :: ngwv                   ! number of gravity waves to use
    INTEGER, INTENT(in) :: kldv(pcols)            ! top interface of low level stress  divergence region
    INTEGER, INTENT(in) :: kldvmn                 ! min value of kldv
    INTEGER, INTENT(in) :: ksrc(pcols)            ! index of top interface of source region
    INTEGER, INTENT(in) :: ksrcmn                 ! min value of ksrc

    REAL(r8), INTENT(in) :: t(pcols,pver)         ! midpoint temperatures
    REAL(r8), INTENT(in) :: pi(pcols,0:pver)      ! interface pressures
    REAL(r8), INTENT(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    REAL(r8), INTENT(in) :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    REAL(r8), INTENT(in) :: piln(pcols,0:pver)    ! ln(interface pressures)
    REAL(r8), INTENT(in) :: rhoi(pcols,0:pver)    ! interface density
    REAL(r8), INTENT(in) :: ni(pcols,0:pver)      ! interface Brunt-Vaisalla frequency
    REAL(r8), INTENT(in) :: ti(pcols,0:pver)      ! interface temperature
    REAL(r8), INTENT(in) :: nm(pcols,pver)        ! midpoint Brunt-Vaisalla frequency
    REAL(r8), INTENT(in) :: dt                    ! time step
    REAL(r8), INTENT(in) :: rdpldv(pcols)         ! 1/dp across low level divergence region
    REAL(r8), INTENT(in) :: ubi(pcols,0:pver)     ! projection of wind at interfaces
    REAL(r8), INTENT(in) :: ubm(pcols,pver)       ! projection of wind at midpoints
    REAL(r8), INTENT(in) :: xv(pcols)             ! unit vectors of source wind (x)
    REAL(r8), INTENT(in) :: yv(pcols)             ! unit vectors of source wind (y)
    REAL(r8), INTENT(in) :: effgw                 ! tendency efficiency

    REAL(r8), INTENT(inout) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress

    REAL(r8), INTENT(out) :: ut(pcols,pver)       ! zonal wind tendency
    REAL(r8), INTENT(out) :: vt(pcols,pver)       ! meridional wind tendency
    REAL(r8), INTENT(out) :: tau0x(pcols)         ! c=0 sfc. stress (zonal)
    REAL(r8), INTENT(out) :: tau0y(pcols)         ! c=0 sfc. stress (meridional)

    !---------------------------Local storage-------------------------------
    INTEGER :: i,k,l                              ! loop indexes

    REAL(r8) :: d(pcols)                          ! "total" diffusivity 
    REAL(r8) :: dsat(pcols,-pgwv:pgwv)            ! saturation diffusivity
    REAL(r8) :: dscal                             ! fraction of dsat to use
    REAL(r8) :: mi                                ! imaginary part of vertical wavenumber
    REAL(r8) :: taudmp                            ! stress after damping
    REAL(r8) :: taumax(pcols)                     ! max(tau) for any l
    REAL(r8) :: tausat(pcols,-pgwv:pgwv)          ! saturation stress
    REAL(r8) :: ubmc(pcols,-pgwv:pgwv)            ! (ub-c)
    REAL(r8) :: ubmc2                             ! (ub-c)**2
    REAL(r8) :: ubt(pcols,pver)                   ! ubar tendency
    REAL(r8) :: ubtl                              ! ubar tendency from wave l
    REAL(r8) :: ubtlsat                           ! saturation tendency

    ! Initialize gravity wave drag tendencies to zero

    DO k=1,pver
       DO i=1,pcols
          ut(i,k) = 0.0_r8
          vt(i,k) = 0.0_r8
       END DO
    END DO

    !---------------------------------------------------------------------------
    ! Compute the stress profiles and diffusivities
    !---------------------------------------------------------------------------

    ! Loop from bottom to top to get stress profiles      

    !DO k = kbot-1, ktop, -1
    DO k = pver  , ktop, -1

       ! Determine the absolute value of the saturation stress and the diffusivity
       ! for each wave.
       ! Define critical levels where the sign of (u-c) changes between interfaces.

       d(:ncol) = dback
       DO l = -ngwv, ngwv
          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                ubmc(i,l) = ubi(i,k) - c(l)
                tausat(i,l) = ABS (effkwv * rhoi(i,k) * ubmc(i,l)**3 / (2.0_r8*ni(i,k)) )
                IF (tausat(i,l) .LE. taumin) tausat(i,l) = 0.0_r8
                IF (ubmc(i,l) / (ubi(i,k+1) - c(l)) .LE. 0.0_r8) tausat(i,l) = 0.0_r8
                dsat(i,l) = (ubmc(i,l) / ni(i,k))**2 * &
                     (effkwv * ubmc(i,l)**2 / (rog * ti(i,k) * ni(i,k)) - alpha(k))
                dscal = MIN (1.0_r8, tau(i,l,k+1) / (tausat(i,l)+taumin))
                d(i) = MAX( d(i), dscal * dsat(i,l))
             END IF
          END DO
       END DO

       ! Compute stress for each wave. The stress at this level is the min of 
       ! the saturation stress and the stress at the level below reduced by damping.
       ! The sign of the stress must be the same as at the level below.

       DO l = -ngwv, ngwv
          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                ubmc2 = MAX(ubmc(i,l)**2, ubmc2mn)
                mi = ni(i,k) / (2.0_r8 * kwv * ubmc2) * (alpha(k) + ni(i,k)**2/ubmc2 * d(i))
                taudmp = tau(i,l,k+1) * EXP(-2.0_r8*mi*rog*t(i,k+1)*(piln(i,k+1)-piln(i,k)))
                IF (taudmp .LE. taumin) taudmp = 0.0_r8
                tau(i,l,k) = MIN (taudmp, tausat(i,l))
             END IF
          END DO
       END DO

       ! The orographic stress term does not change across the source region
       ! Note that k ge ksrcmn cannot occur without an orographic source term

       IF (k .GE. ksrcmn) THEN
          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                IF (k .GE. ksrc(i)) THEN
                   tau(i,0,k) = tau(i,0,pver) 
                END IF
             END IF
          END DO
       END IF

       ! Require that the orographic term decrease linearly (with pressure) 
       ! within the low level stress divergence region. This supersedes the above
       ! requirment of constant stress within the source region.
       ! Note that k ge kldvmn cannot occur without an orographic source term, since
       ! kldvmn=pver then and k<=pver-1

       IF (k .GE. kldvmn) THEN
          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                IF (k .GE. kldv(i)) THEN
                   tau(i,0,k) = MIN (tau(i,0,k), tau(i,0,pver)  * &
                        (1.0_r8 - fracldv * (pi(i,k)-pi(i,pver)) * rdpldv(i)))
                END IF
             END IF
          END DO
       END IF

       ! Apply lower bounds to the stress if ngwv > 0.

       IF (ngwv .GE. 1) THEN

          ! Determine the max value of tau for any l

          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                 taumax(i) = tau(i,-ngwv,k)
             END IF
          END DO
          DO l = -ngwv+1, ngwv
             DO i = 1, ncol
                IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                   taumax(i) = MAX(taumax(i), tau(i,l,k))
                END IF
             END DO
          END DO
          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                taumax(i) = mxrange * taumax(i)
             END IF
          END DO

          ! Set the min value of tau for each wave to the max of mxrange*taumax or
          ! mxasym*tau(-c)

          DO l = 1, ngwv
             DO i = 1, ncol
                IF(k <= kbot(i)-1 .and. k>= ktop)THEN

                   tau(i, l,k) = MAX(tau(i, l,k), taumax(i))
                   tau(i, l,k) = MAX(tau(i, l,k), mxasym*tau(i,-l,k))
                   tau(i,-l,k) = MAX(tau(i,-l,k), taumax(i))
                   tau(i,-l,k) = MAX(tau(i,-l,k), mxasym*tau(i, l,k))
                END IF
             END DO
          END DO
          l = 0
          DO i = 1, ncol
             IF(k <= kbot(i)-1 .and. k>= ktop)THEN
                tau(i,l,k) = MAX(tau(i,l,k), mxasym * 0.5_r8 * (tau(i,l-1,k) + tau(i,l+1,k)))
             END IF 
          END DO
       END IF

    END DO

    ! Put an upper bound on the stress at the top interface to force some stress
    ! divergence in the top layer. This prevents the easterlies from running
    ! away in the summer mesosphere, since most of the gravity wave activity
    ! will pass through a top interface at 75--80 km under strong easterly
    ! conditions. 
    !++BAB fix to match ccm3.10
!!$    do l = -ngwv, ngwv
!!$       do i = 1, ncol
!!$          tau(i,l,0) = min(tau(i,l,0), 0.5*tau(i,l,1))
!!$       end do
!!$    end do
    !--BAB fix to match ccm3.10

    !---------------------------------------------------------------------------
    ! Compute the tendencies from the stress divergence.
    !---------------------------------------------------------------------------

    ! Loop over levels from top to bottom
    !DO k = ktop+1, kbot
    DO k = ktop+1, pver

       ! Accumulate the mean wind tendency over wavenumber.
       DO i = 1, ncol
          IF(k >= ktop+1 .AND. k<= kbot(i) )THEN
             ubt (i,k) = 0.0_r8
          END IF
       END DO
       DO l = -ngwv, ngwv
          DO i = 1, ncol
             IF(k >= ktop+1 .AND. k<= kbot(i) )THEN

                ! Determine the wind tendency including excess stress carried down from above.
                ubtl = g * (tau(i,l,k)-tau(i,l,k-1)) * rdpm(i,k)

                ! Require that the tendency be no larger than the analytic solution for
                ! a saturated region [proportional to (u-c)^3].
                ubtlsat = effkwv * ABS((c(l)-ubm(i,k))**3) /(2.0_r8*rog*t(i,k)*nm(i,k))
                ubtl = MIN(ubtl, ubtlsat)

                ! Apply tendency limits to maintain numerical stability.
                ! 1. du/dt < |c-u|/dt  so u-c cannot change sign (u^n+1 = u^n + du/dt * dt)
                ! 2. du/dt < tndmax    so that ridicuously large tendency are not permitted
                ubtl = MIN(ubtl, umcfac * ABS(c(l)-ubm(i,k)) / dt)
                ubtl = MIN(ubtl, tndmax)

                ! Accumulate the mean wind tendency over wavenumber.
                ubt (i,k) = ubt (i,k) + SIGN(ubtl, c(l)-ubm(i,k))

                ! Redetermine the effective stress on the interface below from the wind 
                ! tendency. If the wind tendency was limited above, then the new stress
                ! will be small than the old stress and will cause stress divergence in
                ! the next layer down. This has the effect of smoothing large stress 
                ! divergences downward while conserving total stress.
                tau(i,l,k) = tau(i,l,k-1) + ubtl * dpm(i,k) / g
             END IF
          END DO
       END DO

       ! Project the mean wind tendency onto the components and scale by "efficiency".
       DO i = 1, ncol
          IF(k >= ktop+1 .AND. k<= kbot(i) )THEN
             ut(i,k) = ubt(i,k) * xv(i) * effgw
             vt(i,k) = ubt(i,k) * yv(i) * effgw
          END IF
       END DO

       ! End of level loop
    END DO

    !-----------------------------------------------------------------------
    ! Project the c=0 stress (scaled) in the direction of the source wind for recording
    ! on the output file.
    !-----------------------------------------------------------------------
    DO i = 1, ncol
       tau0x(i) = tau(i,0,kbot(i)) * xv(i) * effgw
       tau0y(i) = tau(i,0,kbot(i)) * yv(i) * effgw
    END DO

    RETURN
  END SUBROUTINE gw_drag_prof

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
       piln   ,  pint   , pmid   , pdel   , rpdel  , &
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


END MODULE GwddSchemeCAM

!PROGRAM main
! USE GwddSchemeCAM
! IMPLICIT NONE
!
!END PROGRAM main
