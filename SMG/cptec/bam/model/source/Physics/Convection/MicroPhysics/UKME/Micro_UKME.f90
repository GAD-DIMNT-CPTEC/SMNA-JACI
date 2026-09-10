MODULE Micro_UKME
  USE StratCloudFraction, ONLY:       &
      Init_StratCloudFraction,Run_StratCloudFraction

  IMPLICIT NONE
SAVE

  PRIVATE
  !LS_PPN____LSPCON_GAMMAF
  !      |
  !      |___LS_PPNC_______LSP_ICE_______QSAT_mix 
  !                |              |         
  !                |__LSP_SCAV    |______qsat_wat_mix
  !                               |
  !                               |______QSAT_WAT
  
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers
  !#########################################################################################################
  REAL(KIND=r8), PARAMETER :: Epsilonk  = 0.62198_r8
  REAL(KIND=r8), PARAMETER :: Epsilon   = 0.62198_r8
  REAL(KIND=r8), PARAMETER :: C_Virtual = 1.0_r8/Epsilon-1.0_r8
  
  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(R8),PARAMETER :: SHR_CONST_CPWV   = 1.810e3_R8      ! specific heat of water vap ~ J/kg/K
  REAL(R8),PARAMETER :: SHR_CONST_G      = 9.80616_R8      ! acceleration of gravity ~ m/s^2
  REAL(R8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_R8       ! molecular weight dry air ~ kg/kmole  
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(R8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
  REAL(R8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole  
  REAL(R8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(R8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg

  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg

  REAL(r8),PARAMETER :: SHR_CONST_Pi                 = 3.14159265358979323846_r8
  REAL(r8),PARAMETER :: rair   = SHR_CONST_RDAIR    ! Gas constant for dry air (J/K/kg)
  REAL(r8),PARAMETER :: gravit = SHR_CONST_G      ! gravitational acceleration
  REAL(r8),PARAMETER :: zvir   = SHR_CONST_RWV/SHR_CONST_RDAIR - 1          ! rh2o/rair - 1
  REAL(KIND=r8), PARAMETER  :: Earth_Radius = 6371229.0_r8
  REAL(KIND=r8), PARAMETER :: recip_a2  =1.0_r8/(earth_radius*earth_radius)     ! 1/(radius of earth)^2
  REAL (KIND=r8), PARAMETER   :: gasr  =                  287.05_r8! gas constant of dry air        (j/kg/k)

  !*----------------------------------------------------------------------
  !*L------------------COMDECK C_O_DG_C-----------------------------------
  ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
  ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
  ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

  REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
  REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
  REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

  !ftn -g  -fbounds-check -Waliasing -Wall -fbacktrace -ffpe-trap=invalid,overflow,zero  -finit-real=nan -finit-integer=nan -fconvert=big-endian  -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2 
  REAL(KIND=r8), PARAMETER :: one_minus_epsilon = 1.0_r8 -epsilonk
  !      One minus the ratio of the molecular weights of water and dry
  !      air

  REAL(KIND=r8), PARAMETER :: T_low  = 183.15_r8 ! Lowest temperature for which look_up table of
  ! saturation water vapour presssure is valid (K)
  REAL(KIND=r8), PARAMETER :: T_high = 338.15_r8 ! Highest temperature for which look_up table of
  ! saturation water vapour presssure is valid (K)
  REAL(KIND=r8), PARAMETER :: delta_T = 0.1_r8    ! temperature increment of the look-up table
  ! of saturation vapour pressures.

  INTEGER, PARAMETER :: N = INT(((T_high - T_low + (delta_T*0.5_r8))/delta_T) + 1.0_r8)
  !      Gives  N=1551, size of the look-up table of saturation water
  !      vapour pressure.

!  REAL(KIND=r8), PARAMETER :: RHCRIT0(1:38) = (/&
!             0.950_r8,0.925_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.850_r8,0.850_r8,0.850_r8,0.800_r8, &
!             0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8, &
!             0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8, &
!             0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8/)
  REAL(KIND=r8), PARAMETER :: RHCRIT0(1:38) = (/&
         0.990_r8,0.955_r8,0.940_r8,0.940_r8,0.940_r8,0.940_r8,0.920_r8,0.920_r8,0.920_r8,0.900_r8, &
         0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8, &
         0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8, &
         0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8/)


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
  !#########################################################################################################
  INTEGER :: IESW     ! loop counter for data statement look-up table

  REAL(KIND=r8) :: ESW2(0:N+1)
  !                         TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
  !                       - SET BY DATA STATEMENT CALCULATED FROM THE
  !                         GOFF-GRATCH FORMULAE AS TAKEN FROM LANDOLT-
  !                         BORNSTEIN, 1987 NUMERICAL DATA AND FUNCTIONAL
  !                         RELATIONSHIPS IN SCIENCE AND TECHNOLOGY.
  !                         GROUP V/ VOL 4B METEOROLOGY. PHYSICAL AND
  !                         CHEMICAL PROPERTIES OF AIR, P35
  PUBLIC :: Init_Micro_UKME
  PUBLIC :: RunMicro_UKME
CONTAINS
   SUBROUTINE Init_Micro_UKME(iMax)
    IMPLICIT NONE
    INTEGER , INTENT(IN   ) :: iMax

      CALL qsat_wat_data()
      CALL qsat_data()
      CALL Init_StratCloudFraction(iMax)
   END SUBROUTINE Init_Micro_UKME
   
   INTEGER FUNCTION idx(xMax,yMax,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: xMax
    INTEGER, INTENT(IN   ) :: yMax
    INTEGER, INTENT(IN   ) :: x
    REAL(KIND=r8) :: tag_alfa,a

    tag_alfa=REAL(yMax - 1 )/ REAL(xMax - 1)
    ! x=1 ; y=1
    a=1.0_r8 - tag_alfa
    idx = INT(a + tag_alfa*x)
   END FUNCTION idx

   SUBROUTINE RunMicro_UKME(&
                        kMax           , &
                        nCols          , &
                        prsi           , &
                        prsl           , &
                        TIMESTEP       , &
                        pblh           , &
                        colrad         , &
                        kuo            , &
                        q              , &
                        QCF            , &
                        QCL            , &
                        qcf2           , &
                        T              , &
                        CF             , &
                        CFL            , &
                        CFF            , &
                        dtdt           , &
                        dqdt           , &
                        dqldt          , &
                        dqidt          , &
                        dqrdt          , &
                        SNOW_DEPTH     , &
                        LAND_FRACT     , &
                        topo           , &
                        LSRAIN         , &
                        LSSNOW           )
    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: kMax                                 ! Number of model levels
    INTEGER, INTENT(IN   ) :: nCols           
    REAL(KIND=r8)   , INTENT(in   ) :: prsi      (nCols,kMax+1)  
    REAL(KIND=r8)   , INTENT(in   ) :: prsl      (nCols,kMax  ) 
    REAL(KIND=r8)   , INTENT(IN   ) :: TIMESTEP                 ! IN Timestep (sec).
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh      (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad    (nCols)
    INTEGER         , INTENT(IN   ) :: kuo       (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: Q         (nCols,kMax)   ! Specific humidity  (kg water
    REAL(KIND=r8)   , INTENT(INOUT) :: QCF       (nCols,kMax)   ! Cloud ice          (kg per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCL       (nCols,kMax)   ! Cloud liquid water (kg per
    REAL(KIND=r8)   , INTENT(INOUT) :: T         (nCols,kMax)   ! Temperature        (K).
    REAL(KIND=r8)   , INTENT(INOUT) :: qcf2      (nCols,kMax)   ! Ice                (kg per kg air)
    REAL(KIND=r8)   , INTENT(INOUT) :: CF        (nCols,kMax)   ! IN Cloud fraction.
    REAL(KIND=r8)   , INTENT(INOUT) :: CFL       (nCols,kMax)   ! IN Cloud liquid fraction.
    REAL(KIND=r8)   , INTENT(INOUT) :: CFF       (nCols,kMax)   ! IN Cloud ice fraction.
    REAL(KIND=r8)   , INTENT(IN   ) :: SNOW_DEPTH(nCols)                 !, INTENT(IN   )IN Snow depth (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: LAND_FRACT(nCols)                 !, INTENT(IN   )IN Land fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: topo      (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: LSRAIN    (nCols)        ! OUT Surface rainfall rate (kg / sq m /
    REAL(KIND=r8)   , INTENT(OUT  ) :: LSSNOW    (nCols)        ! OUT Surface snowfall rate (kg / sq m /
    REAL(KINd=r8)   , INTENT(OUT  ) :: dtdt (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqdt (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqldt(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqidt(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqrdt(nCols,kMax)

    
    REAL(KIND=r8)    :: LSRAIN3D  (nCols,kMax)   ! OUT Rain rate out of each model layer
    REAL(KIND=r8)    :: LSSNOW3D  (nCols,kMax)        ! OUT Snow rate out of each model layer
    REAL(KIND=r8)    :: RAINFRAC3D(nCols,kMax)        ! OUT rain fraction out of each model layer




    LOGICAL        :: L_MURK                               !, INTENT(IN   )IN Aerosol needs scavenging.
    REAL(KIND=r8)  :: AEROSOL   (nCols,kMax)   !, INTENT(INOUT)  Aerosol            (K).

    LOGICAL :: L_USE_SULPHATE_AUTOCONV              !, INTENT(IN   )IN True if sulphate aerosol used in autoconversion (i.e. second indirect effect)
    LOGICAL :: L_SEASALT_CCN                        !, INTENT(IN   )IN True if sea-salt required for second indirect effect.
    LOGICAL :: L_BIOMASS_CCN                        !, INTENT(IN   )IN True if biomass smoke used for 2nd indirect effect.
    LOGICAL :: L_OCFF_CCN                           !, INTENT(IN   )IN True if fossil-fuel organic carbon aerosol is used for 2nd indirect effect
    LOGICAL :: L_use_biogenic                       !, INTENT(IN   )IN True if biogenic aerosol used for 2nd indirect effect.
    LOGICAL :: L_AUTO_DEBIAS                        !, INTENT(IN   )IN True if autoconversion de-biasing scheme selected.
    LOGICAL :: L_mcr_qcf2                           !, INTENT(IN   )IN Use second prognostic ice if T
    LOGICAL :: L_it_melting                         !, INTENT(IN   )IN Use iterative melting
    LOGICAL :: L_pc2                                !, INTENT(IN   )IN Use the PC2 cloud and condenstn scheme
    LOGICAL :: L_mixing_ratio                       !, INTENT(IN   )   Use mixing ratios removed by deposition depending on
                                                    !                   amount of ice in each category
    LOGICAL :: L_autoconv_murk                      !, INTENT(IN   ) Use murk aerosol to calc. drop number


    REAL(KIND=r8)    :: SO4_AIT      (nCols,kMax)            !, INTENT(IN   )Sulphur Cycle tracers (mmr kg/kg)
    REAL(KIND=r8)    :: SO4_ACC      (nCols,kMax)            !, INTENT(IN   )Sulphur Cycle tracers (mmr kg/kg)
    REAL(KIND=r8)    :: SO4_DIS      (nCols,kMax)            !, INTENT(IN   )Sulphur Cycle tracers (mmr kg/kg)
    REAL(KIND=r8)    :: BMASS_AGD    (nCols,kMax)            !, INTENT(IN   )Biomass smoke tracers
    REAL(KIND=r8)    :: BMASS_CLD    (nCols,kMax)            !, INTENT(IN   )Biomass smoke tracers
    REAL(KIND=r8)    :: OCFF_AGD     (nCols,kMax)            !, INTENT(IN   )Fossil-fuel organic carbon tracers
    REAL(KIND=r8)    :: OCFF_CLD     (nCols,kMax)            !, INTENT(IN   )Fossil-fuel organic carbon tracers
    REAL(KIND=r8)    :: SEA_SALT_FILM(nCols,kMax)            !, INTENT(IN   )(m-3)
    REAL(KIND=r8)    :: SEA_SALT_JET (nCols,kMax)            !, INTENT(IN   )(m-3)
    REAL(KIND=r8)    :: biogenic     (nCols,kMax)            !, INTENT(IN   )(m.m.r.)
    REAL(KIND=r8)    :: RHCRIT       (nCols,kMax)            !, INTENT(IN   )IN Critical humidity for cloud formation.
    REAL(KIND=r8)    :: ec_auto                              !, INTENT(IN   )Collision coalescence efficiency
    REAL(KIND=r8)    :: Ntot_land                            !, INTENT(IN   )Number of droplets over land / m-3
    REAL(KIND=r8)    :: Ntot_sea                             !, INTENT(IN   )Number of droplets over sea / m-3
    REAL(KIND=r8)    :: x1i                                  !, INTENT(IN   )Intercept of aggregate size distribution
    REAL(KIND=r8)    :: x1ic                                 !, INTENT(IN   )Intercept of crystal size distribution
    REAL(KIND=r8)    :: x1r                                  !, INTENT(IN   )Intercept of raindrop size distribution
    REAL(KIND=r8)    :: x2r                                  !, INTENT(IN   )Scaling parameter of raindrop size distribn
    REAL(KIND=r8)    :: x4r                                  !, INTENT(IN   )Shape parameter of raindrop size distribution
    INTEGER          :: ERROR                                ! OUT Return code - 0 if OK, 1 if bad arguments.
    REAL(KIND=r8)    :: rho(nCols,kMax)     

    REAL(KIND=r8)    :: flip_pint  (nCols,kMax+1)   ! Interface pressures  
    REAL(KIND=r8)    :: flip_pmid  (nCols,kMax)     ! Midpoint pressures 
    REAL(KIND=r8)    :: flip_t     (nCols,kMax)     ! temperature
    REAL(KIND=r8)    :: flip_q     (nCols,kMax)     ! specific humidity
    REAL(KIND=r8)    :: flip_pdel  (nCols,kMax)     ! layer thickness
    REAL(KIND=r8)    :: flip_rpdel (nCols,kMax)     ! inverse of layer thickness
    REAL(KIND=r8)    :: flip_lnpmid(nCols,kMax)     ! Log Midpoint pressures    
    REAL(KIND=r8)    :: flip_lnpint(nCols,kMax+1)   ! Log interface pressures
    REAL(KIND=r8)    :: flip_zi   (1:nCols,1:kMax+1)! Height above surface at interfaces 
    REAL(KIND=r8)    :: flip_zm   (1:nCols,1:kMax)  ! Geopotential height at mid level

    REAL(KIND=r8)    :: zi   (nCols,kMax+1)     ! Height above surface at interfaces
    REAL(KIND=r8)    :: zm   (nCols,kMax)        ! Geopotential height at mid level
    REAL(KIND=r8)    :: r_rho_levels       (1:nCols,1:kMax) ! height of rho level above earth's centre (m)
    REAL(KIND=r8)    :: r_theta_levels     (1:nCols,0:kMax) ! height of theta level above earth's centre (m)
    REAL(KIND=r8)    :: p_layer_boundaries (1:nCols,0:kMax)    !, INTENT(IN   )
    REAL(KIND=r8)    :: rho_r2             (1:nCols,1:kMax)    !, INTENT(IN   )IN Air density * earth radius**2
    REAL(KIND=r8)    :: p                  (1:nCols,1:kMax+1)  ! pressure at all points, on u,v levels (Pa).
    REAL(KIND=r8)    :: p_theta_levels     (1:nCols,1:kMax)    ! pressure at all points (Pa).

    REAL(KIND=r8)    :: T_mic         (nCols,kMax)   ! Temperature        (K).
    REAL(KIND=r8)    :: Q_mic         (nCols,kMax)   ! Specific humidity  (kg water
    REAL(KIND=r8)    :: QCF_mic       (nCols,kMax)   ! Cloud ice          (kg per kg air).
    REAL(KIND=r8)    :: QCL_mic       (nCols,kMax)   ! Cloud liquid water (kg per
    REAL(KIND=r8)    :: qcf2_mic      (nCols,kMax)   ! Ice                (kg per kg air)

    INTEGER :: i,k
    INTEGER :: kflip


    DO i=1,nCols
       !flip_pint       (i,kMax+1) = gps(i)*si(1) ! gps --> Pa
       flip_pint       (i,kMax+1) = prsi(i,1) 
    END DO

    DO k=kMax,1,-1
       kflip=kMax+2-k
       DO i=1,nCols
          !flip_pint    (i,k)      = MAX(si(kflip)*gps(i) ,1.0e-12_r8)
          flip_pint    (i,k)      = MAX(prsi(i,kflip),1.0e-12_r8)
       END DO
    END DO
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols
          flip_t   (i,kflip) =  T (i,k)
          flip_q   (i,kflip) =  Q (i,k)
          !flip_pmid(i,kflip) =  sl(  k)*gps (i)
          flip_pmid(i,kflip) =  prsl(i,k)
       END DO
    END DO

    DO k=1,kMax
       DO i=1,nCols    
          dtdt (i,k)=0.0_r8
          dqdt (i,k)=0.0_r8
          dqldt(i,k)=0.0_r8
          dqidt(i,k)=0.0_r8
          dqrdt(i,k)=0.0_r8

          T_mic         (i,k) = T         (i,k)
          Q_mic         (i,k) = Q         (i,k)
          QCF_mic       (i,k) = QCF       (i,k)
          QCL_mic       (i,k) = QCL       (i,k)
          qcf2_mic      (i,k) = qcf2      (i,k)

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
         flip_lnpint(1:nCols,1:kMax+1)   , flip_pint (1:nCols,1:kMax+1)  , &
         flip_pmid  (1:nCols,1:kMax)     , flip_pdel  (1:nCols,1:kMax)   , flip_rpdel(1:nCols,1:kMax)   , &
         flip_t     (1:nCols,1:kMax)     , flip_q     (1:nCols,1:kMax)   , rair   , gravit , zvir   ,&
         flip_zi    (1:nCols,1:kMax+1)   , flip_zm    (1:nCols,1:kMax)   , nCols   ,nCols, kMax)

    DO i=1,nCols
          zi(i,1) = flip_zi    (i,kMax+1)
          p (i,1) = flip_pint    (i,kMax+1)
    END DO
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols
          zi (i,k+1) = flip_zi    (i,kflip)
          zm (i,k  ) = flip_zm    (i,kflip)
          p  (i,k+1) = flip_pint  (i,kflip)
          p_theta_levels(i,k) = flip_pmid(i,kflip)
       END DO
    END DO
    !
    DO i=1,nCols
       R_THETA_LEVELS(i,0)      = EARTH_RADIUS + MAX(topo(i),0.0_r8) + zi(i,1)
       p_layer_boundaries (i,0) = flip_pint(i,kMax+1)
    END DO
    
    DO k=1,kMax
       kflip=kMax+1-k
       DO i=1,nCols
          R_THETA_LEVELS     (i,k) = (EARTH_RADIUS + MAX(topo(i),0.0_r8)) + zi    (i,k+1)
          r_rho_levels       (i,k) = (EARTH_RADIUS + MAX(topo(i),0.0_r8)) + zm    (i,k)
          p_theta_levels     (i,k) = flip_pmid(i,kflip)    !, INTENT(IN)
          p_layer_boundaries (i,k) = flip_pint(i,kflip)    !, INTENT(IN)
          !PRINT*,R_THETA_LEVELS       (i,k),state_zi    (i,kMax+1-k)

          !j/kg/kelvin
          !
          ! P = rho * R * T
          !
          !            P
          ! rho  = -------
          !          R * T
          !
          !           1             R * T
          !  rrho = ----- =        -------
          !          rho              P
          !
          !
          rho   (i,k) = (flip_pmid(i,kflip)/(gasr*flip_t   (i,kflip)))         ! density

          rho_r2(i,k) = rho(i,k) *( r_rho_levels(i,k) * r_rho_levels(i,k) )    !, INTENT(IN)IN Air density * earth radius**2 !         Reference density of air (kg m-3)*m**2

       END DO
    END DO

    L_MURK                 =.FALSE.    !, INTENT(IN   )IN Aerosol needs scavenging.
    L_autoconv_murk        =.FALSE.     !, INTENT(IN   )IN Use murk aerosol to calc. drop number
    AEROSOL=0.0_r8

    L_USE_SULPHATE_AUTOCONV=.FALSE.    !, INTENT(IN   )IN True if sulphate aerosol used in autoconversion (i.e. second indirect effect)
    L_SEASALT_CCN          =.FALSE.    !, INTENT(IN   )IN True if sea-salt required for second indirect effect.
    L_BIOMASS_CCN          =.FALSE.    !, INTENT(IN   )IN True if biomass smoke used for 2nd indirect effect.
    L_OCFF_CCN             =.FALSE.    !, INTENT(IN   )IN True if fossil-fuel organic carbon aerosol is used for 2nd indirect effect
    L_use_biogenic         =.FALSE.    !, INTENT(IN   )IN True if biogenic aerosol used for 2nd indirect effect.
    L_AUTO_DEBIAS          =.TRUE.     !, INTENT(IN   )IN True if autoconversion de-biasing scheme selected.
    L_mcr_qcf2             =.TRUE.     !, INTENT(IN   )IN Use second prognostic ice if T
    L_it_melting           =.TRUE.     !, INTENT(IN   )IN Use iterative melting
    L_pc2                  =.FALSE.    !, INTENT(IN   )IN Use the PC2 cloud and condenstn scheme
    L_mixing_ratio         =.FALSE.    !, INTENT(IN   )IN Use mixing ratios removed by deposition depending on

    EC_AUTO=0.425_r8             !Collision coalescence efficiency
    NTOT_LAND=3.0E8_r8           ! Number of droplets over land / m-3
    NTOT_SEA=1.0E8_r8            ! Number of droplets over sea / m-3
                                       !                  amount of ice in each category
    X1I=4.0E6_r8    !, INTENT(IN)Intercept of aggregate size distribution
    X1IC=80.0E6_r8  !, INTENT(IN)Intercept of crystal size distribution
    X1R  =26.2_r8   !, INTENT(IN)Intercept of raindrop size distribution
    x2r  =1.57_r8   !, INTENT(IN)Scaling parameter of raindrop size distribn
    x4r  = 1.0_r8   !, INTENT(IN)Shape parameter of raindrop size distribution

    DO k=1,kMax
       DO i=1,nCols
          SO4_AIT      (i,k)=0.0_r8
          SO4_ACC      (i,k)=0.0_r8
          SO4_DIS      (i,k)=0.0_r8
          BMASS_AGD    (i,k)=0.0_r8
          BMASS_CLD    (i,k)=0.0_r8
          OCFF_AGD     (i,k)=0.0_r8
          OCFF_CLD     (i,k)=0.0_r8
          SEA_SALT_FILM(i,k)=0.0_r8
          SEA_SALT_JET (i,k)=0.0_r8
          biogenic     (i,k)=0.0_r8
          RHCRIT       (i,k)=RHCRIT0(idx(kMax,38,k))
       END DO
    END DO    

    CALL Run_StratCloudFraction(&
                        nCols                             , &!INTEGER   , INTENT(IN   ) :: nCols
                        kMax                              , &!INTEGER   , INTENT(IN   ) :: LEVELS
                        L_mcr_qcf2                        , &!LOGICAL   , INTENT(IN   ) :: L_mcr_qcf2              ! IN Use second prognostic ice if T
                        l_mixing_ratio                    , &!LOGICAL   ,INTENT(IN   ) :: l_mixing_ratio! true if using mixing ratio formulation
                        RHCRIT         (1:nCols,1:kMax  ) , &!REAL(KIND=r8), INTENT(IN   ) :: rhcrit  (nCols,LEVELS) 
                        zi             (1:nCols,1:kMax+1) , &!REAL(KIND=r8), INTENT(IN   ) :: zi      (1:nCols,1:LEVELS+1)
                        zm             (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(IN   ) :: zm      (1:nCols,1:LEVELS) 
                        pblh           (1:nCols)          , &!REAL(KIND=r8), INTENT(IN   ) :: pblh    (1:nCols)   !
                        colrad         (1:nCols)          , &!REAL(KIND=r8), INTENT(IN   ) :: colrad  (1:nCols)   !
                        p              (1:nCols,1:kMax+1) , &!REAL(KIND=r8), INTENT(IN   ) :: p    (1:nCols,1:LEVELS+1)   ! pressure at all points, on theta levels (Pa).
                        p_theta_levels (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(IN   ) :: p_theta_levels(nCols,levels)    ! pressure at all points (Pa).
                        R_THETA_LEVELS (1:nCols,0:kMax)   , &!REAL(KIND=r8), INTENT(IN   ) :: r_theta_levels(1:nCols,0:LEVELS) ! height of theta level above earth's centre (m)
                        kuo            (1:nCols)          , &!INTEGER      , INTENT(IN   ) :: kuo     (1:nCols) 
                        T_mic          (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(INOUT) :: T3      (1:nCols,LEVELS)
                        Q_mic          (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(INOUT) :: q       (1:nCols,LEVELS)
                        QCL_mic        (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(INOUT) :: qcl     (1:nCols,LEVELS)
                        QCF_mic        (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(IN   ) :: qcf     (1:nCols,LEVELS)
                        qcf2_mic       (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(IN   ) :: qcf2    (1:nCols,LEVELS)
                        CF             (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(INOUT) :: CF      (nCols,LEVELS) !Cloud fraction at processed levels (decimal fraction).
                        CFL            (1:nCols,1:kMax)   , &!REAL(KIND=r8), INTENT(INOUT) :: CFL     (nCols,LEVELS) !Liquid cloud fraction at processed levels (decimal fraction).
                        CFF            (1:nCols,1:kMax)     )!REAL(KIND=r8), INTENT(INOUT) :: CFF     (nCols,LEVELS) !Frozen cloud fraction at processed levels (decimal fraction).

      CALL LS_PPN(&
       p_layer_boundaries     , &!REAL   , INTENT(IN   ) :: p_layer_boundaries (nCols,0:kMax)    !
       p_theta_levels         , &!REAL   , INTENT(IN   ) :: p_theta_levels     (nCols,kMax)      !
       TIMESTEP               , &!REAL   , INTENT(IN   ) :: TIMESTEP                ! IN Timestep (sec).
       L_autoconv_murk        , &!LOGICAL, INTENT(IN   ) :: L_autoconv_murk! Use murk aerosol to calc. drop number
       ec_auto                , &!REAL   , INTENT(IN   ) :: ec_auto        ! Collision coalescence efficiency
       Ntot_land              , &!REAL   , INTENT(IN   ) :: Ntot_land      ! Number of droplets over land / m-3
       Ntot_sea               , &!REAL   , INTENT(IN   ) :: Ntot_sea       ! Number of droplets over sea / m-3
       x1i                    , &!REAL   , INTENT(IN   ) :: x1i            ! Intercept of aggregate size distribution
       x1ic                   , &!REAL   , INTENT(IN   ) :: x1ic           ! Intercept of crystal size distribution
       x1r                    , &!REAL   , INTENT(IN   ) :: x1r            ! Intercept of raindrop size distribution
       x2r                    , &!REAL   , INTENT(IN   ) :: x2r            ! Scaling parameter of raindrop size distribn
       x4r                    , &!REAL   , INTENT(IN   ) :: x4r            ! Shape parameter of raindrop size distribution
       CF                     , &!REAL   , INTENT(IN   ) :: CF  (nCols,kMax)      ! IN Cloud fraction.
       CFL                    , &!REAL   , INTENT(IN   ) :: CFL                (nCols,kMax)      ! IN Cloud liquid fraction.
       CFF                    , &!REAL   , INTENT(IN   ) :: CFF                (nCols,kMax)      ! IN Cloud ice fraction.
       RHCRIT                 , &!REAL   , INTENT(IN   ) :: RHCRIT             (nCols,kMax)      ! IN Critical humidity for cloud formation.
       kMax                   , &!INTEGER, INTENT(IN   ) :: kMax                                 ! Number of model levels
       rho_r2                 , &!REAL   , INTENT(IN   ) :: rho_r2             (1:nCols,kMax)    ! IN Air density * earth radius**2
       r_rho_levels           , &!REAL   , INTENT(IN   ) :: r_rho_levels       (1:nCols,kMax)    ! IN Earths radius at each rho level
       r_theta_levels         , &!REAL   , INTENT(IN   ) :: r_theta_levels     (1:nCols,0:kMax)  ! IN Earths radius at each theta level
       Q_mic                  , &!REAL   , INTENT(INOUT) :: Q      (nCols,kMax)   ! Specific humidity  (kg water
       QCF_mic                , &!REAL   , INTENT(INOUT) :: QCF    (nCols,kMax)   ! Cloud ice          (kg per kg air).
       QCL_mic                , &!REAL   , INTENT(INOUT) :: QCL    (nCols,kMax)   ! Cloud liquid water (kg per
       T_mic                  , &!REAL   , INTENT(INOUT) :: T      (nCols,kMax)   ! Temperature        (K).
       qcf2_mic               , &!REAL   , INTENT(INOUT) :: qcf2   (nCols,kMax)   ! Ice                (kg per kg air)
       L_mcr_qcf2             , &!LOGICAL, INTENT(IN   ) :: L_mcr_qcf2              ! IN Use second prognostic ice if T
       L_it_melting           , &!LOGICAL, INTENT(IN   ) :: L_it_melting            ! IN Use iterative melting
       L_USE_SULPHATE_AUTOCONV, &!LOGICAL, INTENT(IN   ) :: L_USE_SULPHATE_AUTOCONV ! IN True if sulphate aerosol used in autoconversion (i.e. second indirect effect)
       L_AUTO_DEBIAS          , &!LOGICAL, INTENT(IN   ) :: L_AUTO_DEBIAS           ! IN True if autoconversion de-biasing scheme selected.
       l_mixing_ratio         , &!LOGICAL, INTENT(IN   ) :: L_mixing_ratio ! Use mixing ratios removed by deposition depending on
       SEA_SALT_FILM          , &!REAL   , INTENT(IN   ) :: SEA_SALT_FILM(nCols,kMax)! (m-3)
       SEA_SALT_JET           , &!REAL   , INTENT(IN   ) :: SEA_SALT_JET (nCols,kMax)! (m-3)
       L_SEASALT_CCN          , &!LOGICAL, INTENT(IN   ) :: L_SEASALT_CCN           ! IN True if sea-salt required for second indirect effect.
       L_use_biogenic         , &!LOGICAL, INTENT(IN   ) :: L_use_biogenic          ! IN True if biogenic aerosol used for 2nd indirect effect.
       biogenic               , &!REAL   , INTENT(IN   ) :: biogenic     (nCols,kMax)           ! (m.m.r.)
       SNOW_DEPTH             , &!REAL   , INTENT(IN   ) :: SNOW_DEPTH(nCols)        ! IN Snow depth (m)
       LAND_FRACT             , &!REAL   , INTENT(IN   ) :: LAND_FRACT(nCols)        ! IN Land fraction
       SO4_AIT                , &!REAL   , INTENT(IN   ) :: SO4_AIT   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
       SO4_ACC                , &!REAL   , INTENT(IN   ) :: SO4_ACC   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
       SO4_DIS                , &!REAL   , INTENT(IN   ) :: SO4_DIS   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
       BMASS_AGD              , &!REAL   , INTENT(IN   ) :: BMASS_AGD (nCols,kMax)   !Biomass smoke tracers
       BMASS_CLD              , &!REAL   , INTENT(IN   ) :: BMASS_CLD (nCols,kMax)   !Biomass smoke tracers
       L_BIOMASS_CCN          , &!LOGICAL, INTENT(IN   ) :: L_BIOMASS_CCN           ! IN True if biomass smoke used for 2nd indirect effect.
       OCFF_AGD               , &!REAL   , INTENT(IN   ) :: OCFF_AGD  (nCols,kMax)   !Fossil-fuel organic carbon tracers
       OCFF_CLD               , &!REAL   , INTENT(IN   ) :: OCFF_CLD  (nCols,kMax)   !Fossil-fuel organic carbon tracers
       L_OCFF_CCN             , &!LOGICAL, INTENT(IN   ) :: L_OCFF_CCN              ! IN True if fossil-fuel organic carbon aerosol is used for 2nd indirect effect
       AEROSOL                , &!REAL   , INTENT(INOUT) :: AEROSOL(nCols,kMax)   ! Aerosol            (K).
       L_MURK                 , &!LOGICAL, INTENT(IN   ) :: L_MURK                  ! IN Aerosol needs scavenging.
       L_pc2                  , &!LOGICAL, INTENT(IN   ) :: L_pc2                   ! IN Use the PC2 cloud and condenstn scheme
       LSRAIN                 , &!REAL   , INTENT(OUT  ) :: LSRAIN(nCols)       ! OUT Surface rainfall rate (kg / sq m /
       LSSNOW                 , &!REAL   , INTENT(OUT  ) :: LSSNOW(nCols)       ! OUT Surface snowfall rate (kg / sq m /
       LSRAIN3D               , &!REAL   , INTENT(OUT  ) :: LSRAIN3D(nCols,kMax)! OUT Rain rate out of each model layer
       LSSNOW3D               , &!REAL   , INTENT(OUT  ) :: LSSNOW3D(nCols,kMax)! OUT Snow rate out of each model layer
       RAINFRAC3D             , &!REAL   , INTENT(OUT  ) :: RAINFRAC3D(nCols,kMax)! OUT rain fraction out of each model layer
       nCols                  , &!INTEGER, INTENT(IN   ) :: nCols           
       ERROR                    &!INTEGER, INTENT(OUT  ) ::  ERROR          ! OUT Return code - 0 if OK, 1 if bad arguments.
       )

    DO k=1,kMax
       DO i=1,nCols    
          dtdt (i,k)=( T_mic   (i,k) - T    (i,k) )/TIMESTEP
          dqdt (i,k)=( Q_mic   (i,k) - Q    (i,k) )/TIMESTEP
          dqldt(i,k)=( QCF_mic (i,k) - QCF  (i,k) )/TIMESTEP
          dqidt(i,k)=( QCL_mic (i,k) - QCL  (i,k) )/TIMESTEP
          dqrdt(i,k)=( qcf2_mic(i,k) - qcf2 (i,k) )/TIMESTEP



          T    (i,k)  = T_mic   (i,k)
          Q    (i,k)  = Q_mic   (i,k)
          QCF  (i,k)  = QCF_mic (i,k)
          QCL  (i,k)  = QCL_mic (i,k)
          qcf2 (i,k)  = qcf2_mic(i,k)
       END DO
    END DO

    ! (kg m^-2 s^-1).  0.001m/s
    !
    ! after cumulus parameterization
    ! out  tn1, qn1, prec, kuo,ktop, kbot
    !
    DO i = 1,nCols
      LSRAIN(i) = 0.5_r8*LSRAIN(i)*TIMESTEP/1000.0_r8! convert mm/s to meters
      LSSNOW(i) = 0.5_r8*LSSNOW(i)*TIMESTEP/1000.0_r8! convert mm/s to meters
    ! RAINCV(i)=dtime*pre1(i)             !in mm/sec(ditme),by 0.5_r8(if leap-frog or 2dt)
    ! RAINCV(i)=RAINCV(i)/1000.0_r8       !in m for gcm
    END DO


   END SUBROUTINE RunMicro_UKME


  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !*LL  SUBROUTINES LS_PPN and LS_PPNC------------------------------------
  !LL  Purpose:
  !LL          LS_PPN and LS_PPNC:
  !LL           Calculate large-scale (dynamical) precipitation.
  !LL           LS_PPNC is the gather/scatter routine which then
  !LL           calls LSP_ICE.
  !LL  Note: in all cases, level counters (incl subscripts) run from 1
  !LL        (lowest model layer) to kMax (topmost "wet" model
  !LL        layer) - it is assumed that the bottom kMax layers are
  !LL        the "wet" layers.
  !LL
  !LL  Put through fpp on Cray.  Activate *IF definition CRAY if running
  !LL  on the Cray.
  !LL
  !LL  Modification History from Version 4.4
  !LL     Version    Date
  !LL       4.5      March 98          New Deck        Damian Wilson
  !
  !         5.1      16-12-99   Change to allow 3D RHcrit diagnostic
  !                             in place of 1D parameter. AC Bushell
  !         5.2      25/10/00   Reintroduction of tracers. P.Selwood
  !         5.2      21-11-00   Allow interfacting with 3C scheme.
  !                             Return diagnostics. Damian Wilson
  !LL
  !         5.2      28-11-00   Pass down sea-salt arrays.
  !                                                       A. Jones
  !         5.3      29-08-01   Introduce logical to control effect of
  !                             sulphate aerosol on autoconversion.
  !                                                       A. Jones
  !         5.3      24-09-01   Tidy up code (mostly redundant) relating
  !                             to the sulphur cycle.     A. Jones
  !         5.4      30-08-02   Remove lscav_agedsoot, which is now
  !                             calculated in microphys_ctl.  P. Davison
  !         5.4      27-07-02   Change cloud fractions to prognostics
  !                             for PC2                   D. Wilson
  !
  !         5.4      25-04-02   Pass land fraction down to LSP_ICE.
  !                                                       A. Jones
  !         5.5      03-02-03   Pass extra microphysics variables
  !                             down to LSP_ICE.           R.M.Forbes
  !         5.5      03-02-03   Include extra microphysics variables
  !                             (qcf2,qrain,qgraup)          R.M.Forbes
  !         6.1      01-08-04   Include variables for prognostic rain
  !                             and graupel.             R.M.Forbes
  !         6.1      07-04-04   Add biomass smoke to call to LSP_ICE.
  !                                                          A. Jones
  !         6.1      07-04-04   Pass switch for autoconversion de-biasing
  !                             down to LSP_ICE.              A. Jones
  !         6.2      22-08-05   Remove commented out code. P.Selwood.
  !         6.2      23-11-05   Pass through precip variables from
  !                             UMUI. Damian Wilson
  !         6.2      11-01-06   Include non-hydrostatic calculation
  !                             of air density and model thickness.
  !                                                       D. Wilson
  !         6.2      31-01-06   Pass process rate diags through. R.Forbes
  !         6.4      10-01-07   Provide mixing ratio logical control
  !                                                       D. Wilson
  !LL  Programming standard: Unified Model Documentation Paper No 3
  !LL
  !LL  Documentation: UM Documentation Paper 26.
  !LL
  !*L  Arguments:---------------------------------------------------------
  SUBROUTINE LS_PPN(&
       p_layer_boundaries     , &!REAL   , INTENT(IN   ) :: p_layer_boundaries (nCols,0:kMax)    !
       p_theta_levels         , &!REAL   , INTENT(IN   ) :: p_theta_levels     (nCols,kMax)      !
       TIMESTEP               , &!REAL   , INTENT(IN   ) :: TIMESTEP                ! IN Timestep (sec).
       L_autoconv_murk        , &!LOGICAL, INTENT(IN   ) :: L_autoconv_murk! Use murk aerosol to calc. drop number
       ec_auto                , &!REAL   , INTENT(IN   ) :: ec_auto        ! Collision coalescence efficiency
       Ntot_land              , &!REAL   , INTENT(IN   ) :: Ntot_land      ! Number of droplets over land / m-3
       Ntot_sea               , &!REAL   , INTENT(IN   ) :: Ntot_sea       ! Number of droplets over sea / m-3
       x1i                    , &!REAL   , INTENT(IN   ) :: x1i            ! Intercept of aggregate size distribution
       x1ic                   , &!REAL   , INTENT(IN   ) :: x1ic           ! Intercept of crystal size distribution
       x1r                    , &!REAL   , INTENT(IN   ) :: x1r            ! Intercept of raindrop size distribution
       x2r                    , &!REAL   , INTENT(IN   ) :: x2r            ! Scaling parameter of raindrop size distribn
       x4r                    , &!REAL   , INTENT(IN   ) :: x4r            ! Shape parameter of raindrop size distribution
       CF                     , &!REAL   , INTENT(IN   ) :: CF  (nCols,kMax)      ! IN Cloud fraction.
       CFL                    , &!REAL   , INTENT(IN   ) :: CFL                (nCols,kMax)      ! IN Cloud liquid fraction.
       CFF                    , &!REAL   , INTENT(IN   ) :: CFF                (nCols,kMax)      ! IN Cloud ice fraction.
       RHCRIT                 , &!REAL   , INTENT(IN   ) :: RHCRIT             (nCols,kMax)      ! IN Critical humidity for cloud formation.
       kMax                   , &!INTEGER, INTENT(IN   ) :: kMax                                 ! Number of model levels
       rho_r2                 , &!REAL   , INTENT(IN   ) :: rho_r2             (1:nCols,kMax)    ! IN Air density * earth radius**2
       r_rho_levels           , &!REAL   , INTENT(IN   ) :: r_rho_levels       (1:nCols,kMax)    ! IN Earths radius at each rho level
       r_theta_levels         , &!REAL   , INTENT(IN   ) :: r_theta_levels     (1:nCols,0:kMax)  ! IN Earths radius at each theta level
       q                      , &!REAL   , INTENT(INOUT) :: Q      (nCols,kMax)   ! Specific humidity  (kg water
       qcf                    , &!REAL   , INTENT(INOUT) :: QCF    (nCols,kMax)   ! Cloud ice          (kg per kg air).
       qcl                    , &!REAL   , INTENT(INOUT) :: QCL    (nCols,kMax)   ! Cloud liquid water (kg per
       T                      , &!REAL   , INTENT(INOUT) :: T      (nCols,kMax)   ! Temperature        (K).
       qcf2                   , &!REAL   , INTENT(INOUT) :: qcf2   (nCols,kMax)   ! Ice                (kg per kg air)
       L_mcr_qcf2             , &!LOGICAL, INTENT(IN   ) :: L_mcr_qcf2              ! IN Use second prognostic ice if T
       L_it_melting           , &!LOGICAL, INTENT(IN   ) :: L_it_melting            ! IN Use iterative melting
       L_USE_SULPHATE_AUTOCONV, &!LOGICAL, INTENT(IN   ) :: L_USE_SULPHATE_AUTOCONV ! IN True if sulphate aerosol used in autoconversion (i.e. second indirect effect)
       L_AUTO_DEBIAS          , &!LOGICAL, INTENT(IN   ) :: L_AUTO_DEBIAS           ! IN True if autoconversion de-biasing scheme selected.
       l_mixing_ratio         , &!LOGICAL, INTENT(IN   ) :: L_mixing_ratio ! Use mixing ratios removed by deposition depending on
       SEA_SALT_FILM          , &!REAL   , INTENT(IN   ) :: SEA_SALT_FILM(nCols,kMax)! (m-3)
       SEA_SALT_JET           , &!REAL   , INTENT(IN   ) :: SEA_SALT_JET (nCols,kMax)! (m-3)
       L_SEASALT_CCN          , &!LOGICAL, INTENT(IN   ) :: L_SEASALT_CCN           ! IN True if sea-salt required for second indirect effect.
       L_use_biogenic         , &!LOGICAL, INTENT(IN   ) :: L_use_biogenic          ! IN True if biogenic aerosol used for 2nd indirect effect.
       biogenic               , &!REAL   , INTENT(IN   ) :: biogenic     (nCols,kMax)           ! (m.m.r.)
       SNOW_DEPTH             , &!REAL   , INTENT(IN   ) :: SNOW_DEPTH(nCols)        ! IN Snow depth (m)
       LAND_FRACT             , &!REAL   , INTENT(IN   ) :: LAND_FRACT(nCols)        ! IN Land fraction
       SO4_AIT                , &!REAL   , INTENT(IN   ) :: SO4_AIT   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
       SO4_ACC                , &!REAL   , INTENT(IN   ) :: SO4_ACC   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
       SO4_DIS                , &!REAL   , INTENT(IN   ) :: SO4_DIS   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
       BMASS_AGD              , &!REAL   , INTENT(IN   ) :: BMASS_AGD (nCols,kMax)   !Biomass smoke tracers
       BMASS_CLD              , &!REAL   , INTENT(IN   ) :: BMASS_CLD (nCols,kMax)   !Biomass smoke tracers
       L_BIOMASS_CCN          , &!LOGICAL, INTENT(IN   ) :: L_BIOMASS_CCN           ! IN True if biomass smoke used for 2nd indirect effect.
       OCFF_AGD               , &!REAL   , INTENT(IN   ) :: OCFF_AGD  (nCols,kMax)   !Fossil-fuel organic carbon tracers
       OCFF_CLD               , &!REAL   , INTENT(IN   ) :: OCFF_CLD  (nCols,kMax)   !Fossil-fuel organic carbon tracers
       L_OCFF_CCN             , &!LOGICAL, INTENT(IN   ) :: L_OCFF_CCN              ! IN True if fossil-fuel organic carbon aerosol is used for 2nd indirect effect
       AEROSOL                , &!REAL   , INTENT(INOUT) :: AEROSOL(nCols,kMax)   ! Aerosol            (K).
       L_MURK                 , &!LOGICAL, INTENT(IN   ) :: L_MURK                  ! IN Aerosol needs scavenging.
       L_pc2                  , &!LOGICAL, INTENT(IN   ) :: L_pc2                   ! IN Use the PC2 cloud and condenstn scheme
       LSRAIN                 , &!REAL   , INTENT(OUT  ) :: LSRAIN(nCols)       ! OUT Surface rainfall rate (kg / sq m /
       LSSNOW                 , &!REAL   , INTENT(OUT  ) :: LSSNOW(nCols)       ! OUT Surface snowfall rate (kg / sq m /
       LSRAIN3D               , &!REAL   , INTENT(OUT  ) :: LSRAIN3D(nCols,kMax)! OUT Rain rate out of each model layer
       LSSNOW3D               , &!REAL   , INTENT(OUT  ) :: LSSNOW3D(nCols,kMax)! OUT Snow rate out of each model layer
       RAINFRAC3D             , &!REAL   , INTENT(OUT  ) :: RAINFRAC3D(nCols,kMax)! OUT rain fraction out of each model layer
       nCols                  , &!INTEGER, INTENT(IN   ) :: nCols           
       ERROR                    &!INTEGER, INTENT(OUT  ) ::  ERROR          ! OUT Return code - 0 if OK, 1 if bad arguments.
       )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kMax                                 ! Number of model levels
    INTEGER, INTENT(IN) :: nCols           
    REAL(KIND=r8)   , INTENT(INOUT) :: CF                 (nCols,kMax)      ! IN Cloud fraction.
    REAL(KIND=r8)   , INTENT(IN) :: p_theta_levels     (nCols,kMax)      !
    REAL(KIND=r8)   , INTENT(IN) :: p_layer_boundaries (nCols,0:kMax)    !
    REAL(KIND=r8)   , INTENT(IN) :: RHCRIT             (nCols,kMax)      ! IN Critical humidity for cloud formation.
    REAL(KIND=r8)   , INTENT(INOUT) :: CFL                (nCols,kMax)      ! IN Cloud liquid fraction.
    REAL(KIND=r8)   , INTENT(INOUT) :: CFF                (nCols,kMax)      ! IN Cloud ice fraction.
    REAL(KIND=r8)   , INTENT(IN) :: rho_r2             (1:nCols,kMax)    ! IN Air density * earth radius**2
    REAL(KIND=r8)   , INTENT(IN) :: r_rho_levels       (1:nCols,kMax)    ! IN Earths radius at each rho level
    REAL(KIND=r8)   , INTENT(IN) :: r_theta_levels     (1:nCols,0:kMax)  ! IN Earths radius at each theta level
    REAL(KIND=r8)   , INTENT(IN) :: TIMESTEP                ! IN Timestep (sec).
    LOGICAL, INTENT(IN) :: L_MURK                  ! IN Aerosol needs scavenging.
    LOGICAL, INTENT(IN) :: L_USE_SULPHATE_AUTOCONV ! IN True if sulphate aerosol used in autoconversion (i.e. second indirect effect)
    LOGICAL, INTENT(IN) :: L_SEASALT_CCN           ! IN True if sea-salt required for second indirect effect.
    LOGICAL, INTENT(IN) :: L_BIOMASS_CCN           ! IN True if biomass smoke used for 2nd indirect effect.
    LOGICAL, INTENT(IN) :: L_OCFF_CCN              ! IN True if fossil-fuel organic carbon aerosol is used for 2nd indirect effect
    LOGICAL, INTENT(IN) :: L_use_biogenic          ! IN True if biogenic aerosol used for 2nd indirect effect.
    LOGICAL, INTENT(IN) :: L_AUTO_DEBIAS           ! IN True if autoconversion de-biasing scheme selected.
    LOGICAL, INTENT(IN) :: L_mcr_qcf2              ! IN Use second prognostic ice if T
    LOGICAL, INTENT(IN) :: L_it_melting            ! IN Use iterative melting
    LOGICAL, INTENT(IN) :: L_pc2                   ! IN Use the PC2 cloud and condenstn scheme

    REAL(KIND=r8), INTENT(INOUT) :: Q      (nCols,kMax)   ! Specific humidity  (kg water
    REAL(KIND=r8), INTENT(INOUT) :: QCF    (nCols,kMax)   ! Cloud ice          (kg per kg air).
    REAL(KIND=r8), INTENT(INOUT) :: QCL    (nCols,kMax)   ! Cloud liquid water (kg per
    REAL(KIND=r8), INTENT(INOUT) :: T      (nCols,kMax)   ! Temperature        (K).
    REAL(KIND=r8), INTENT(INOUT) :: qcf2   (nCols,kMax)   ! Ice                (kg per kg air)
    REAL(KIND=r8), INTENT(INOUT) :: AEROSOL(nCols,kMax)   ! Aerosol            (K).

    REAL(KIND=r8),INTENT(INOUT) :: SO4_AIT   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
    REAL(KIND=r8),INTENT(INOUT) :: SO4_ACC   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
    REAL(KIND=r8),INTENT(INOUT) :: SO4_DIS   (nCols,kMax)   !Sulphur Cycle tracers (mmr kg/kg)
    REAL(KIND=r8),INTENT(INOUT) :: BMASS_AGD (nCols,kMax)   !Biomass smoke tracers
    REAL(KIND=r8),INTENT(INOUT) :: BMASS_CLD (nCols,kMax)   !Biomass smoke tracers
    REAL(KIND=r8),INTENT(INOUT) :: OCFF_AGD  (nCols,kMax)   !Fossil-fuel organic carbon tracers
    REAL(KIND=r8),INTENT(INOUT) :: OCFF_CLD  (nCols,kMax)   !Fossil-fuel organic carbon tracers
    REAL(KIND=r8),INTENT(IN   ) :: SNOW_DEPTH(nCols)        ! IN Snow depth (m)
    REAL(KIND=r8),INTENT(IN   ) :: LAND_FRACT(nCols)        ! IN Land fraction
    REAL(KIND=r8),INTENT(IN   ) :: SEA_SALT_FILM(nCols,kMax)! (m-3)
    REAL(KIND=r8),INTENT(IN   ) :: SEA_SALT_JET (nCols,kMax)! (m-3)
    REAL(KIND=r8),INTENT(IN   ) :: biogenic     (nCols,kMax)           ! (m.m.r.)
    LOGICAL, INTENT(IN) :: L_mixing_ratio ! Use mixing ratios removed by deposition depending on
                                          ! amount of ice in each category
    LOGICAL, INTENT(IN) :: L_autoconv_murk! Use murk aerosol to calc. drop number

    REAL(KIND=r8)   , INTENT(IN) :: ec_auto        ! Collision coalescence efficiency
    REAL(KIND=r8)   , INTENT(IN) :: Ntot_land      ! Number of droplets over land / m-3
    REAL(KIND=r8)   , INTENT(IN) :: Ntot_sea       ! Number of droplets over sea / m-3
    REAL(KIND=r8)   , INTENT(IN) :: x1i            ! Intercept of aggregate size distribution
    REAL(KIND=r8)   , INTENT(IN) :: x1ic           ! Intercept of crystal size distribution
    REAL(KIND=r8)   , INTENT(IN) :: x1r            ! Intercept of raindrop size distribution
    REAL(KIND=r8)   , INTENT(IN) :: x2r            ! Scaling parameter of raindrop size distribn
    REAL(KIND=r8)   , INTENT(IN) :: x4r            ! Shape parameter of raindrop size distribution

    REAL(KIND=r8)   , INTENT(OUT  ) :: LSRAIN(nCols)       ! OUT Surface rainfall rate (kg / sq m /
    REAL(KIND=r8)   , INTENT(OUT  ) :: LSSNOW(nCols)       ! OUT Surface snowfall rate (kg / sq m /
    REAL(KIND=r8)   , INTENT(OUT  ) :: LSRAIN3D(nCols,kMax)! OUT Rain rate out of each model layer
    REAL(KIND=r8)   , INTENT(OUT  ) :: LSSNOW3D(nCols,kMax)! OUT Snow rate out of each model layer
    REAL(KIND=r8)   , INTENT(OUT  ) :: RAINFRAC3D(nCols,kMax)! OUT rain fraction out of each model layer
    INTEGER, INTENT(OUT  ) ::  ERROR          ! OUT Return code - 0 if OK, 1 if bad arguments.


    !*L  Workspace usage ---------------------------------------------------
    !
    LOGICAL :: L_non_hydrostatic ! Use non-hydrostatic formulation of layer thicknesses

    !ajm      LOGICAL
    !ajm     & H(PFIELD)      ! Used as "logical" in compression.
    INTEGER :: IX(nCols)    ! Index for compress/expand.
    REAL(KIND=r8)    :: VFALL(nCols)        ! snow fall velocity (m per s).
    REAL(KIND=r8)    :: VFALL2(nCols)       ! fall velocity for qcf2 (m/s)
    REAL(KIND=r8)    :: LSSNOW2(nCols)      ! snowfall rate for qcf2
    REAL(KIND=r8)    :: VFALL_RAIN(nCols)   ! fall velocity for rain (m/s)
    REAL(KIND=r8)    :: VFALL_GRAUP(nCols)  ! fall vel. for graupel (m/s)
    REAL(KIND=r8)    :: CTTEMP(nCols)
    REAL(KIND=r8)    :: RAINFRAC(nCols)
    REAL(KIND=r8)    :: FRAC_ICE_ABOVE(nCols) ! Cloud ice fraction
    !                                            in layer above
    REAL(KIND=r8)    :: layer_thickness(nCols)
    REAL(KIND=r8)    :: rho1(nCols)
    REAL(KIND=r8)    :: rho2(nCols)
    REAL(KIND=r8)    :: deltaz(nCols,kMax)
    REAL(KIND=r8)    :: rhodz_dry(nCols,kMax)
    REAL(KIND=r8)    :: rhodz_moist(nCols,kMax)
    REAL(KIND=r8)    :: q_total(nCols)
    REAL(KIND=r8)    :: LSGRAUP(nCols)      ! Graupel fall rate (kg/m2/s)

    INTEGER :: biog_dim_ice ! local ! Array dimension for passing to LSP_ICE.
    INTEGER :: biogenic_ptr ! local ! Pointer for biogenic aerosol array.
    INTEGER :: salt_dim_ice    ! Array dimension for passing to LSP_ICE.
    INTEGER :: sea_salt_ptr    ! Pointer for sea-salt arrays.

    ! Allocate CX and CONSTP arrays
    ! Start C_LSPSIZ
    ! Description: Include file containing idealised forcing options
    ! Author:      R. Forbes
    !
    ! History:
    ! Version  Date      Comment
    ! -------  ----      -------
    !   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
    !   6.2    22/08/05  Include the step size between ice categories.
    !                                                   Damian Wilson

    ! Sets up the size of arrays for CX and CONSTP
    REAL(KIND=r8)    :: CX(80)
    REAL(KIND=r8)    :: CONSTP(80)
    INTEGER,PARAMETER:: ice_type_offset=20

    ! End C_LSPSIZ
    ! --------------------------COMDECK C_LSPMIC----------------------------
    ! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
    ! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
    ! THE MICROPHYSICS AND ICE CLOUD FRACTION METHOD
    ! ----------------------------------------------------------------------
    !
    ! History:
    !
    ! Version    Date     Comment
    ! -------    ----     -------
    !   5.4    16/08/02   Correct comment line, add PC2 parameters and
    !                     move THOMO to c_micro     Damian Wilson
    !   6.0    11/08/03   Correct value of wind_shear_factor for PC2
    !                                                          Damian Wilson
    !   6.2    17/11/05   Remove variables that are now in UMUI. D. Wilson
    !   6.2    03/02/06   Include droplet settling logical. Damian Wilson
    !
    ! ----------------------------------------------------------------------
    !      AUTOCONVERSION TERMS
    ! ----------------------------------------------------------------------
    !
    !     LOGICAL, PARAMETER :: L_AUTOCONV_MURK is set in UMUI
    ! Set to .TRUE. to calculate droplet concentration from MURK aerosol,
    ! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
    ! of sulphate aerosol). If both are .FALSE., droplet concentrations
    ! from comdeck C_MICRO are used to be consistent with the values
    ! used in the radiation scheme.

    ! This next set of parameters is to allow the 3B scheme to
    ! be replicated at 3C/3D
    ! Inhomogeneity factor for autoconversion rate
    REAL(KIND=r8),PARAMETER:: INHOMOG_RATE=1.0_r8

    ! Inhomogeneity factor for autoconversion limit
    REAL(KIND=r8),PARAMETER:: INHOMOG_LIM=1.0_r8

    ! Threshold droplet radius for autoconversion
    REAL(KIND=r8),PARAMETER:: R_THRESH=7.0E-6_r8
    ! End of 3B repeated code

    !Do not alter R_AUTO and N_AUTO since these values are effectively
    ! hard wired into a numerical approximation in the autoconversion
    ! code. EC_AUTO will be multiplied by CONSTS_AUTO

    ! Threshold radius for autoconversion
    REAL(KIND=r8), PARAMETER :: R_AUTO=20.0E-6_r8

    ! Critical droplet number for autoconversion
    REAL(KIND=r8), PARAMETER :: N_AUTO=1000.0_r8

    ! Collision coalesence efficiency for autoconversion
    !      REAL(KIND=r8), PARAMETER :: EC_AUTO is set in UMUI

    ! The autoconversion powers define the variation of the rate with
    ! liquid water content and droplet concentration. The following are
    ! from Tripoli and Cotton

    !  Dependency of autoconversion rate on droplet concentration
    REAL(KIND=r8), PARAMETER :: POWER_DROPLET_AUTO=-0.33333_r8

    ! Dependency of autoconversion rate on water content
    REAL(KIND=r8), PARAMETER :: POWER_QCL_AUTO=2.33333_r8

    ! Dependency of autoconversion rate on air density
    REAL(KIND=r8), PARAMETER :: power_rho_auto=1.33333_r8

    ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
    ! See UM documentation paper 26, equation P26.132_r8

    ! Combination of physical constants
    REAL(KIND=r8), PARAMETER :: CONSTS_AUTO=5907.24_r8

    ! Quantites for calculation of drop number by aerosols.
    ! Need only set if L_AUTOCONV_MURK=.TRUE.  See file C_VISBTY

    ! Scaling concentration (m-3) in droplet number concentration
    REAL(KIND=r8), PARAMETER :: N0_MURK=500.0E6_r8

    ! Scaling mass (kg/kg) in droplet number calculation from aerosols
    REAL(KIND=r8), PARAMETER :: M0_MURK=1.458E-8_r8

    ! Power in droplet number calculation from aerosols
    REAL(KIND=r8), PARAMETER :: POWER_MURK=0.5_r8

    ! Ice water content threshold for graupel autoconversion (kg/m^3)
    REAL(KIND=r8), PARAMETER :: AUTO_GRAUP_QCF_THRESH = 3.E-4_r8

    ! Temperature threshold for graupel autoconversion (degC)
    REAL(KIND=r8), PARAMETER :: AUTO_GRAUP_T_THRESH = -4.0_r8

    ! Temperature threshold for graupel autoconversion
    REAL(KIND=r8), PARAMETER :: AUTO_GRAUP_COEFF = 0.5_r8

    !-----------------------------------------------------------------
    ! Iterations of microphysics
    !-----------------------------------------------------------------

    ! Number of iterations in microphysics.
    INTEGER,PARAMETER :: LSITER=8
    ! Advise 1 iteration for every 10 minutes or less of timestep.

    !-----------------------------------------------------------------
    ! Nucleation of ice
    !-----------------------------------------------------------------

    ! Note that the assimilation scheme uses temperature thresholds
    ! in its calculation of qsat.

    ! Nucleation mass
    REAL(KIND=r8), PARAMETER :: M0=1.0E-12_r8

    ! Maximum Temp for ice nuclei nucleation (deg C)
    REAL(KIND=r8), PARAMETER :: TNUC=-10.0_r8

    ! Maximum temperature for homogenous nucleation is now in c_micro
    ! so that it is available to code outside of section A04.

    !  1.0/Scaling quantity for ice in crystals
    REAL(KIND=r8), PARAMETER :: QCF0=1.0E4_r8       ! This is an inverse quantity

    ! Minimum allowed QCF after microphysics
    REAL(KIND=r8),PARAMETER:: QCFMIN=1.0E-8_r8

    ! 1/scaling temperature in aggregate fraction calculation
    REAL(KIND=r8), PARAMETER :: T_SCALING=0.0384_r8

    !  Minimum temperature limit in calculation  of N0 for ice (deg C)
    REAL(KIND=r8), PARAMETER :: T_AGG_MIN=-45.0_r8

    !-----------------------------------------------------------------
    ! Hallett Mossop process
    !-----------------------------------------------------------------

    ! Switch off Hallett Mossop in this version but allow
    ! functionality

    ! Min temp for production of Hallett Mossop splinters (deg C)
    REAL(KIND=r8), PARAMETER :: HM_T_MIN=-8.0_r8

    ! Max temp for production of Hallett Mossop splinters (deg C)
    REAL(KIND=r8), PARAMETER :: HM_T_MAX=-273.0_r8
    ! REAL(KIND=r8), PARAMETER :: HM_T_MAX=-3.0

    !  Residence distance for Hallett Mossop splinters (1/deg C)
    REAL(KIND=r8), PARAMETER :: HM_DECAY=1.0_r8/7.0_r8

    ! Reciprocal of scaling liquid water content for HM process
    REAL(KIND=r8), PARAMETER :: HM_RQCL=1.0_r8/0.1E-3_r8

    !-----------------------------------------------------------------
    ! PC2 Cloud Scheme Terms
    !-----------------------------------------------------------------

    ! Specifies the ice content (in terms of a fraction of qsat_liq)
    ! that corresponds to a factor of two reduction in the width of
    ! the vapour distribution in the liquid-free part of the gridbox.
    REAL(KIND=r8), PARAMETER :: ICE_WIDTH=0.04_r8

    ! Parameter that governs the rate of spread of ice cloud fraction
    ! due to windshear
    REAL(KIND=r8), PARAMETER :: WIND_SHEAR_FACTOR = 1.5E-4_r8

    !-----------------------------------------------------------------
    ! Droplet settling
    !-----------------------------------------------------------------
    ! Use the droplet settling code (not available for 3C).
    LOGICAL, PARAMETER :: l_droplet_settle = .FALSE.

    !  External subroutines called -----------------------------------------
    !EXTERNAL LS_PPNC,LSPCON
    !*----------------------------------------------------------------------
    !  Physical constants -------------------------------------------------
    REAL(KIND=r8)    ,PARAMETER :: CFMIN=1.0E-3_r8    ! Used for LS_PPNC  compress.
    !  Define local variables ----------------------------------------------
    INTEGER :: I,K ! Loop counters: I - horizontal field index;
    !                                      K - vertical level index.
    INTEGER :: N              ! "nval" for WHEN routine.
    !
    REAL(KIND=r8)    :: work  ! work variable

    !-----------------------------------------------------------------------
    ERROR=0

    ! Define l_non_hydrostatic and l_mixing_ratio. These would be best
    ! passed from the UMUI, but we will keep them here for the moment.
    !IF (l_mixing_ratio .OR. (h_sect(4)  /=  '03B' .AND. h_sect(4)  /=  '03C') ) THEN
    IF (l_mixing_ratio) THEN
       l_non_hydrostatic = .TRUE.   ! This uses most physically accurate code.
    ELSE
       l_non_hydrostatic = .FALSE.  ! This ensures bit reproducibility
    END IF
    CX       =0.0_r8
    CONSTP   =0.0_r8
    ERROR =0

    ! Define CX and CONSTP values
    ! DEPENDS ON: lspcon
    CALL LSPCON( &
         CX        , & !REAL(KIND=r8), INTENT(OUT  ) :: CX(80)
         CONSTP    , & !REAL, INTENT(OUT  ) :: CONSTP(80)
         x1i       , & !REAL, INTENT(IN   ) ::  x1i  ! Intercept of aggregate size distribution
         x1ic      , & !REAL, INTENT(IN   ) ::  x1ic ! Intercept of crystal size distribution
         x1r       , & !REAL, INTENT(IN   ) ::  x1r  ! Intercept of raindrop size distribution
         x2r       , & !REAL, INTENT(IN   ) ::  x2r  ! Scaling parameter of raindrop size distribution
         x4r         ) !REAL, INTENT(IN   ) ::  x4r  ! Shape parameter of raindrop size distribution
    !-----------------------------------------------------------------------
    !L Internal structure.
    !L 1. Initialise rain and snow to zero.
    !   Initialise scavenged amounts of S Cycle tracers to 0 for full field
    !-----------------------------------------------------------------------
    DO I=1,nCols
       IX      (I)=0
       LSRAIN  (I)=0.0_r8
       LSSNOW  (I)=0.0_r8
       LSSNOW2 (i)=0.0_r8
       LSGRAUP (i)=0.0_r8
       CTTEMP      (I)=0.0_r8
       RAINFRAC    (I)=0.0_r8
       FRAC_ICE_ABOVE(I)=0.0_r8
       VFALL       (I)=0.0_r8
       VFALL2      (i)=0.0_r8
       VFALL_RAIN  (i)=0.0_r8
       VFALL_GRAUP (i)=0.0_r8
       q_total     (i)=0.0_r8
       rho1(i)=0.0_r8
       rho2(i)=0.0_r8
       layer_thickness(i)=0.0_r8
    END DO ! Loop over points,i
    DO k = 1, kMax
       DO i = 1, nCols
          LSRAIN3D   (i,k)=0.0_r8! OUT Rain rate out of each model layer
          LSSNOW3D   (i,k)=0.0_r8! OUT Snow rate out of each model layer
          RAINFRAC3D (i,k)=0.0_r8! OUT rain fraction out of each model layer
          deltaz     (i,k)=0.0_r8
          rhodz_dry  (i,k)=0.0_r8
          rhodz_moist(i,k)=0.0_r8
       END DO ! Loop over points,i
    END DO ! Loop over points,i

    IF (l_non_hydrostatic) THEN
       ! ----------------------------------------------------------------------
       ! Calculate the (non-hydrostatic) layer thicknesses (deltaz) and air
       ! densities multiplied by deltaz (rhodz_moist and rhodz_dry).
       ! ----------------------------------------------------------------------
       ! We should note that this formulation, although better than the
       ! hydrostatic formulation, is still not entirely conservative. To ensure
       ! conservation we would need to rewrite the large-scale precipitation
       ! scheme to consider masses in terms of rho<q>, and
       ! not the current <rho>q formulation.

       ! We only need to calculate averages for the moist levels
       DO k = 1, kMax
          DO i = 1, nCols

             ! Calculate densities at the boundaries of the layer
             ! by removing the r**2 term from rho_r2.
             ! Rho1 is the density at the lower boundary.
             rho1(i)= rho_r2(i,k)/( r_rho_levels(i,k) * r_rho_levels(i,k) )

             ! Check whether there is a rho level above the current
             ! moist level.
             IF (k  <   kMax) THEN
                ! Rho2 is the density at the upper boundary.
                rho2(i)= rho_r2(i,k+1)/( r_rho_levels(i,k+1) * r_rho_levels(i,k+1) )

                ! Calculate the average value of rho across the layer
                ! multiplied by the layer thickness and the layer
                ! thickness.
                rhodz_moist(i,k) = rho2(i) * ( r_theta_levels(i,k) - r_rho_levels(i,k) ) &
                     &                          +  rho1(i) * ( r_rho_levels(i,k+1) - r_theta_levels(i,k) )
                deltaz(i,k) = r_rho_levels(i,k+1) - r_rho_levels(i,k)

                IF (k  ==  1) THEN
                   ! For the lowest layer we need to extend the lower
                   ! boundary from the first rho level to the surface.
                   ! The surface is the 0'th theta level.
                   deltaz     (i,1) = r_rho_levels(i,2) - r_theta_levels(i,0)
                   rhodz_moist(i,1) = rhodz_moist(i,1)*deltaz(i,1) &
                        &                   / (r_rho_levels(i,2)-r_rho_levels(i,1))
                END IF  ! k  ==  1

             ELSE
                ! For a top layer higher than the highest rho level
                ! we can calculate a pseudo rho level. We will assume
                ! it has a similar density to the rho level below
                ! and that the intervening theta level is in the centre
                ! of the layer.
                deltaz     (i,k) = 2.0_r8*(r_theta_levels(i,k) - r_rho_levels(i,k))
                rhodz_moist(i,k) = rho1(i) * deltaz(i,k)
             END IF  ! k  <   kMax

             ! Calculate total moisture
             q_total(i) = q(i,k) + qcl(i,k) + qcf(i,k)

             IF (l_mcr_qcf2) THEN
                q_total(i) = q_total(i) + qcf2(i,k)
             END IF  ! l_mcr_qcf2




             ! Rho_r2 uses the moist density of air. If the mixing
             ! ratio framework is in place then we need to also know
             ! the dry density of air.
             IF (l_mixing_ratio) THEN
                rhodz_dry(i,k) = rhodz_moist(i,k)/ (1.0_r8 + q_total(i))
             ELSE
                rhodz_dry(i,k) = rhodz_moist(i,k)* (1.0_r8 - q_total(i))
             END IF  ! l_mixing_ratio

          END DO  ! i
       END DO  ! k

    END IF  ! L_non_hydrostatic
    !
    !-----------------------------------------------------------------------
    !L 2. Loop round levels from top down (counting bottom level as level 1,
    !L    as is standard in the Unified model).
    !-----------------------------------------------------------------------
    !
    DO K=kMax,1,-1


       !-----------------------------------------------------------------------
       !L 2.5 Form INDEX IX to gather/scatter variables in LS_PPNC
       !-----------------------------------------------------------------------
       !
       !  Set index where cloud fraction > CFMIN or where non-zero pptn
       !  Note: whenimd is functionally equivalent to WHENILE (but autotasks).
       !
       !
       N=0
       DO i = 1,nCols
          layer_thickness(i) = p_layer_boundaries(i,k) - p_layer_boundaries(i,k-1)

          ! Set up IF statement to determine whether to call the
          ! microphysics code for this grid box (i.e. if there is
          ! already condensate in the grid box or there is
          ! precipitation about to fall into the grid box)
          work = QCF(i,k)
          ! Include extra microphysics variables if in use
          IF (L_mcr_qcf2 )  work = work + QCF2(i,k) + LSSNOW2(i)
          IF (CFL(i,k) > CFMIN  .OR. (LSRAIN(i)+LSSNOW(i)) > 0.0_r8 .OR. work > 0.0_r8) THEN
             ! include this grid box.
             ! Strictly speaking the CFL > CFMIN clause is too
             ! restrictive since ice nucleation does not require
             ! liquid water, but the code would be very messy.
             N = N + 1
             IX(N) = i
          END IF
       END DO ! Loop over points,i
       !
       IF(N >  0)THEN

          IF (L_SEASALT_CCN) THEN
             sea_salt_ptr=K
             salt_dim_ice=N
          ELSE
             sea_salt_ptr=1
             salt_dim_ice=1
          ENDIF

          IF (L_use_biogenic) THEN
             biogenic_ptr=K
             biog_dim_ice=N
          ELSE
             biogenic_ptr=1
             biog_dim_ice=1
          ENDIF

          ! DEPENDS ON: ls_ppnc
          CALL LS_PPNC( &
               IX                             , & !INTEGER ,INTENT(IN   ) :: IX(nCols*rows,2)! IN gather/scatter index
               N                              , & !INTEGER ,INTENT(IN   ) :: N               ! IN Number of points where pptn non-zero from above or where CF>CFMIN
               TIMESTEP                       , & !REAL    ,INTENT(IN   ) :: TIMESTEP        ! IN Timestep (sec).
               LSRAIN                         , & !REAL    ,INTENT(INOUT) :: LSRAIN(nCols,rows)             ! INOUT Surface rainfall rate (kg m^-2 s^-1).
               LSSNOW                         , & !REAL    ,INTENT(INOUT) :: LSSNOW(nCols,rows)             ! INOUT Surface snowfall rate (kg m^-2 s^-1).
               LSSNOW2                        , & !REAL    ,INTENT(INOUT) :: LSSNOW2(nCols,rows)       ! INOUT layer snowfall rate (kg m^-2 s^-1).
               LSGRAUP                        , & !REAL    ,INTENT(INOUT) :: LSGRAUP(nCols,rows)       ! INOUT layer graupelfall rate (kg m^-2 s^-1)
               CF    (1,K)                    , & !REAL    ,INTENT(INOUT) :: CF (nCols,rows) ! IN Cloud fraction.
               CFL   (1,K)                    , & !REAL    ,INTENT(INOUT) :: CFL(nCols,rows) ! IN Cloud liquid fraction.
               CFF   (1,K)                    , & !REAL    ,INTENT(INOUT) :: CFF(nCols,rows) ! IN Cloud ice fraction.
               QCF   (1,K)                    , & !REAL    ,INTENT(INOUT) :: QCF(nCols,rows) ! INOUT Cloud ice (kg per kg air).
               QCL   (1,K)                    , & !REAL    ,INTENT(INOUT) :: QCL(nCols,rows)            ! INOUT Cloud liquid water (kg per kg air).
               T     (1,K)                    , & !REAL    ,INTENT(INOUT) :: T(nCols,rows)             ! INOUT Temperature (K).
               QCF2  (1,K)                    , & !REAL    ,INTENT(INOUT) :: QCF2(nCols,rows)            ! INOUT Cloud ice2 (kg per kg air).
               L_mcr_qcf2                     , & !LOGICAL ,INTENT(IN   ) :: L_mcr_qcf2                         ! IN Use second prognostic ice if T
               L_it_melting                   , & !LOGICAL ,INTENT(IN   ) :: L_it_melting                   ! IN Use iterative melting
               l_non_hydrostatic              , & !LOGICAL ,INTENT(IN   ) :: L_non_hydrostatic                 ! IN Use non-hydrostatic calculation of model layer depths
               l_mixing_ratio                 , & !LOGICAL ,INTENT(IN   ) :: L_mixing_ratio                 ! IN q is a mixing ratio
               SO4_AIT(1,K)                   , & !REAL    ,INTENT(INOUT) :: SO4_AIT     (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               SO4_ACC(1,K)                   , & !REAL    ,INTENT(INOUT) :: SO4_ACC     (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               SO4_DIS(1,K)                   , & !REAL    ,INTENT(INOUT) :: SO4_DIS     (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               BMASS_AGD(1,K)                 , & !REAL    ,INTENT(INOUT) :: BMASS_AGD   (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               BMASS_CLD(1,K)                 , & !REAL    ,INTENT(INOUT) :: BMASS_CLD   (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               L_BIOMASS_CCN                  , & !LOGICAL ,INTENT(IN   ) :: L_BIOMASS_CCN                  ! IN Biomass smoke aerosols used for
               OCFF_AGD(1,K)                  , & !REAL    ,INTENT(INOUT) :: OCFF_AGD    (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               OCFF_CLD(1,K)                  , & !REAL    ,INTENT(INOUT) :: OCFF_CLD    (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
               L_OCFF_CCN                     , & !LOGICAL ,INTENT(IN   ) :: L_OCFF_CCN                         ! IN Fossil-fuel organic carbon
               AEROSOL(1,K)                   , & !REAL    ,INTENT(INOUT) :: AEROSOL(nCols,rows)            ! INOUT Aerosol (K).
               L_MURK                         , & !LOGICAL ,INTENT(IN   ) :: L_MURK                         ! IN Aerosol needs scavenging.
               l_pc2                          , & !LOGICAL ,INTENT(IN   ) :: L_pc2                          ! IN Use PC2 cloud and condensation
               SNOW_DEPTH                     , & !REAL    ,INTENT(IN   ) :: SNOW_DEPTH(nCols, rows)
               LAND_FRACT                     , & !REAL    ,INTENT(IN   ) :: LAND_FRACT(nCols, rows)
               L_USE_SULPHATE_AUTOCONV        , & !LOGICAL ,INTENT(IN   ) :: L_USE_SULPHATE_AUTOCONV         ! IN Sulphate aerosols used in
               L_AUTO_DEBIAS                  , & !LOGICAL ,INTENT(IN   ) :: L_AUTO_DEBIAS                  ! IN Use autoconversion de-biasing
               SEA_SALT_FILM(1,sea_salt_ptr), & !REAL    ,INTENT(IN   ) :: SEA_SALT_FILM(nCols, salt_dim2)
               SEA_SALT_JET (1,sea_salt_ptr), & ! REAL   ,INTENT(IN   ) :: SEA_SALT_JET(nCols, salt_dim2)
               L_SEASALT_CCN                  , & !LOGICAL ,INTENT(IN   ) :: L_SEASALT_CCN                  ! IN Sea-salt aerosols used for second indirect effect if T
               salt_dim_ice                   , & !INTEGER ,INTENT(IN   ) :: salt_dim_ice  ! Dimension to use for call to LSP_ICE.
               L_use_biogenic                 , & !LOGICAL ,INTENT(IN   ) :: L_biogenic_CCN                 ! IN Biogenic aerosol used for second indirect effect if T
               biogenic(1, biogenic_ptr)      , & !REAL    ,INTENT(IN   ) :: biogenic(nCols, rows)
               biog_dim_ice                   , & !INTEGER ,INTENT(IN   ) :: biog_dim_ice  !     "     "        "   "        "   "          "   .
               Q(1,K)                         , & !REAL    ,INTENT(INOUT) :: Q(nCols,rows)             ! INOUT Specific humidity (kg water/kg air).
               p_theta_levels(1,K)            , & !REAL    ,INTENT(IN   ) :: p_theta_levels(nCols,rows)
               layer_thickness                , & !REAL    ,INTENT(IN   ) :: layer_thickness(nCols,rows)! IN thickness of layer (Pa)
               deltaz(1,k)                    , & !REAL    ,INTENT(IN   ) :: deltaz(nCols,rows)            ! IN thickness of layer (m)
               rhodz_dry(1,k)                 , & !REAL    ,INTENT(IN   ) :: rhodz_dry(nCols,rows)     ! Dry air density layer thickness (kg m-2)
               rhodz_moist(1,k)               , & !REAL    ,INTENT(IN   ) :: rhodz_moist(nCols,rows)   ! Moist air density layer thickness (kg m-2)
               nCols                          , & !INTEGER ,INTENT(IN   ) :: nCols    ! IN gather/scatter index
               RHCRIT(1,K)                    , & !REAL    ,INTENT(IN   ) :: RHCRIT(nCols,rhc_rows)! IN Critical humidity for cloud formation.
               VFALL                          , & !REAL    ,INTENT(INOUT) :: VFALL(nCols,rows)            ! INOUT fall velocity of ice (m per
               VFALL2                         , & !REAL    ,INTENT(INOUT) :: VFALL2(nCols,rows)            ! INOUT fall vel. of rain (m/s)
               VFALL_RAIN                     , & !REAL    ,INTENT(INOUT) :: VFALL_RAIN(nCols,rows)    ! INOUT fall vel. of rain (m/s)
               VFALL_GRAUP                    , & !REAL    ,INTENT(INOUT) :: VFALL_GRAUP(nCols,rows)   ! INOUT fall vel. of graupel (m/s)
               FRAC_ICE_ABOVE                 , & !REAL    ,INTENT(INOUT) :: FRAC_ICE_ABOVE(nCols,rows)! INOUT Ice fraction from layer above water.
               CTTEMP                         , & !REAL    ,INTENT(INOUT) :: CTTEMP(nCols,rows)            ! INOUT Ice cloud top temperature (K)
               RAINFRAC                       , & !REAL    ,INTENT(INOUT) :: RAINFRAC(nCols,rows)      ! INOUT Rain fraction.
               CX                             , & !REAL    ,INTENT(IN   ) ::  CX(80)
               CONSTP                         , & !REAL    ,INTENT(IN   ) ::  CONSTP(80)
               L_autoconv_murk                , & !LOGICAL ,INTENT(IN   ) :: L_autoconv_murk                 ! Use murk aerosol to calc. drop number
               ec_auto                        , & !REAL    ,INTENT(IN   ) :: ec_auto                         ! Collision coalescence efficiency
               Ntot_land                      , & !REAL    ,INTENT(IN   ) :: Ntot_land                         ! Number of droplets over land / m-3
               Ntot_sea                         & !REAL    ,INTENT(IN   ) :: Ntot_sea                         ! Number of droplets over sea / m-3
               )
       ENDIF
       !
       ! Copy rainfall and snowfall rates to 3D fields for diagnostic output
       !
       IF (nCols  ==  nCols .AND. kMax  ==  kMax) THEN
          ! Only copy rain and snow to 3D fields if arrays are dimensionalized.
          DO i=1,nCols
             LSRAIN3D  (I,K) = LSRAIN(I) 
             LSSNOW3D  (I,K) = LSSNOW(I) + LSSNOW2(I) + LSGRAUP(I)
             RAINFRAC3D(I,K) = RAINFRAC(I)
          END DO
       ENDIF
       !
    END DO ! Loop over K

    ! Add together ice crystals, snow aggregates and graupel
    ! for surface snow rate (kg/m2/s)

    IF (L_mcr_qcf2) THEN
       DO i = 1,nCols
          LSSNOW(i) = LSSNOW(i) + LSSNOW2(i) + LSGRAUP(i)
       END DO
    END IF

    IF (L_droplet_settle) THEN
       ! Add droplet settling to the large-scale precip rate in order
       ! to be able to have a full water budget with the precipitation
       ! diagnostics.
       DO i = 1,nCols
          LSRAIN(i) = LSRAIN(i) 
       END DO
    END IF

!20  CONTINUE                  ! Branch for error exit
    RETURN
  END SUBROUTINE LS_PPN
  !*LL  SUBROUTINE LS_PPNC------------------------------------------------
  !*L  Arguments:---------------------------------------------------------



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !*LL  SUBROUTINES LS_PPN and LS_PPNC------------------------------------
  !LL  Purpose:
  !LL          LS_PPN and LS_PPNC:
  !LL           Calculate large-scale (dynamical) precipitation.
  !LL           LS_PPNC is the gather/scatter routine which then
  !LL           calls LSP_ICE.
  !LL  Note: in all cases, level counters (incl subscripts) run from 1
  !LL        (lowest model layer) to kMax (topmost "wet" model
  !LL        layer) - it is assumed that the bottom kMax layers are
  !LL        the "wet" layers.
  !LL
  !LL  Put through fpp on Cray.  Activate *IF definition CRAY if running
  !LL  on the Cray.
  !LL
  !LL  Modification History from Version 4.4
  !LL     Version    Date
  !LL       4.5      March 98          New Deck        Damian Wilson
  !
  !         5.1      16-12-99   Change to allow 3D RHcrit diagnostic
  !                             in place of 1D parameter. AC Bushell
  !         5.2      25/10/00   Reintroduction of tracers. P.Selwood
  !         5.2      21-11-00   Allow interfacting with 3C scheme.
  !                             Return diagnostics. Damian Wilson
  !LL
  !         5.2      28-11-00   Pass down sea-salt arrays.
  !                                                       A. Jones
  !         5.3      29-08-01   Introduce logical to control effect of
  !                             sulphate aerosol on autoconversion.
  !                                                       A. Jones
  !         5.3      24-09-01   Tidy up code (mostly redundant) relating
  !                             to the sulphur cycle.     A. Jones
  !         5.4      30-08-02   Remove lscav_agedsoot, which is now
  !                             calculated in microphys_ctl.  P. Davison
  !         5.4      27-07-02   Change cloud fractions to prognostics
  !                             for PC2                   D. Wilson
  !
  !         5.4      25-04-02   Pass land fraction down to LSP_ICE.
  !                                                       A. Jones
  !         5.5      03-02-03   Pass extra microphysics variables
  !                             down to LSP_ICE.           R.M.Forbes
  !         5.5      03-02-03   Include extra microphysics variables
  !                             (qcf2,qrain,qgraup)          R.M.Forbes
  !         6.1      01-08-04   Include variables for prognostic rain
  !                             and graupel.             R.M.Forbes
  !         6.1      07-04-04   Add biomass smoke to call to LSP_ICE.
  !                                                          A. Jones
  !         6.1      07-04-04   Pass switch for autoconversion de-biasing
  !                             down to LSP_ICE.              A. Jones
  !         6.2      22-08-05   Remove commented out code. P.Selwood.
  !         6.2      23-11-05   Pass through precip variables from
  !                             UMUI. Damian Wilson
  !         6.2      11-01-06   Include non-hydrostatic calculation
  !                             of air density and model thickness.
  !                                                       D. Wilson
  !         6.2      31-01-06   Pass process rate diags through. R.Forbes
  !LL  Programming standard: Unified Model Documentation Paper No 3
  !LL
  !LL  Documentation: UM Documentation Paper 26.
  !LL
  !*L  Arguments:---------------------------------------------------------
  !*LL  SUBROUTINE LS_PPNC------------------------------------------------
  !*L  Arguments:---------------------------------------------------------
  SUBROUTINE LS_PPNC(                                               &
       IX                      , &!INTEGER ,INTENT(IN   ) :: IX(nCols*rows,2)! IN gather/scatter index
       N                       , &!INTEGER ,INTENT(IN   ) :: N               ! IN Number of points where pptn non-zero from above or where CF>CFMIN
       TIMESTEP                , &!REAL    ,INTENT(IN   ) :: TIMESTEP        ! IN Timestep (sec).
       LSRAIN                  , &!REAL    , INTENT(INOUT) :: LSRAIN(nCols,rows)        ! INOUT Surface rainfall rate (kg m^-2 s^-1).
       LSSNOW                  , &!REAL    , INTENT(INOUT) :: LSSNOW(nCols,rows)        ! INOUT Surface snowfall rate (kg m^-2 s^-1).
       LSSNOW2                 , &!REAL    , INTENT(INOUT) :: LSSNOW2(nCols,rows)       ! INOUT layer snowfall rate (kg m^-2 s^-1).
       LSGRAUP                 , &!REAL    , INTENT(INOUT) :: LSGRAUP(nCols,rows)       ! INOUT layer graupelfall rate (kg m^-2 s^-1)
       CF                      , &!REAL    ,INTENT(INOUT) :: CF (nCols,rows) ! IN Cloud fraction.
       CFL                     , &!REAL    ,INTENT(INOUT) :: CFL(nCols,rows) ! IN Cloud liquid fraction.
       CFF                     , &!REAL    ,INTENT(INOUT) :: CFF(nCols,rows) ! IN Cloud ice fraction.
       QCF                     , &!REAL   , INTENT(INOUT) :: QCF(nCols,rows) ! INOUT Cloud ice (kg per kg air).
       QCL                     , &!REAL   , INTENT(INOUT) :: QCL(nCols,rows)           ! INOUT Cloud liquid water (kg per kg air).
       T                       , &!REAL   , INTENT(INOUT) :: T(nCols,rows)             ! INOUT Temperature (K).
       QCF2                    , &!REAL   , INTENT(INOUT) :: QCF2(nCols,rows)          ! INOUT Cloud ice2 (kg per kg air).
       L_mcr_qcf2              , &!LOGICAL , INTENT(IN  ) :: L_mcr_qcf2                     ! IN Use second prognostic ice if T
       L_it_melting            , &!LOGICAL , INTENT(IN  ) :: L_it_melting                   ! IN Use iterative melting
       l_non_hydrostatic       , &!LOGICAL , INTENT(IN  ) :: L_non_hydrostatic              ! IN Use non-hydrostatic calculation of model layer depths
       l_mixing_ratio          , &!LOGICAL , INTENT(IN  ) :: L_mixing_ratio                 ! IN q is a mixing ratio
       SO4_AIT                 , &!REAL   , INTENT(INOUT) :: SO4_AIT     (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       SO4_ACC                 , &!REAL   , INTENT(INOUT) :: SO4_ACC     (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       SO4_DIS                 , &!REAL   , INTENT(INOUT) :: SO4_DIS     (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       BMASS_AGD               , &!REAL   , INTENT(INOUT) :: BMASS_AGD   (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       BMASS_CLD               , &!REAL   , INTENT(INOUT) :: BMASS_CLD   (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       L_BIOMASS_CCN           , &!LOGICAL , INTENT(IN  ) :: L_BIOMASS_CCN                  ! IN Biomass smoke aerosols used for
       OCFF_AGD                , &!REAL   , INTENT(INOUT) :: OCFF_AGD    (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       OCFF_CLD                , &!REAL   , INTENT(INOUT) :: OCFF_CLD    (nCols, rows) !INOUT S Cycle tracers & scavngd amoun
       L_OCFF_CCN              , &!LOGICAL , INTENT(IN  ) :: L_OCFF_CCN                     ! IN Fossil-fuel organic carbon
       AEROSOL                 , &!REAL   , INTENT(INOUT) :: AEROSOL(nCols,rows)       ! INOUT Aerosol (K).
       L_MURK                  , &!LOGICAL , INTENT(IN  ) :: L_MURK                         ! IN Aerosol needs scavenging.
       l_pc2                   , &!LOGICAL , INTENT(IN  ) :: L_pc2                          ! IN Use PC2 cloud and condensation
       SNOW_DEPTH              , &!REAL    ,INTENT(IN   ) :: SNOW_DEPTH(nCols, rows)
       LAND_FRACT              , &!REAL    ,INTENT(IN   ) :: LAND_FRACT(nCols, rows)
       L_USE_SULPHATE_AUTOCONV , &!LOGICAL , INTENT(IN  ) :: L_USE_SULPHATE_AUTOCONV        ! IN Sulphate aerosols used in
       L_AUTO_DEBIAS           , &!LOGICAL , INTENT(IN  ) :: L_AUTO_DEBIAS                  ! IN Use autoconversion de-biasing
       SEA_SALT_FILM           , &!REAL    ,INTENT(IN   ) :: SEA_SALT_FILM(nCols, salt_dim2)
       SEA_SALT_JET            , &!REAL    ,INTENT(IN   ) :: SEA_SALT_JET(nCols, salt_dim2)
       L_SEASALT_CCN           , &!LOGICAL , INTENT(IN  ) :: L_SEASALT_CCN                  ! IN Sea-salt aerosols used for second indirect effect if T
       salt_dim_ice            , &!INTEGER ,INTENT(IN   ) :: salt_dim_ice  ! Dimension to use for call to LSP_ICE.
       L_biogenic_CCN          , &!LOGICAL , INTENT(IN  ) :: L_biogenic_CCN                 ! IN Biogenic aerosol used for second indirect effect if T
       biogenic                , &!REAL    ,INTENT(IN   ) :: biogenic(nCols, rows)
       biog_dim_ice            , &!INTEGER ,INTENT(IN   ) :: biog_dim_ice  !     "     "   "   "   "   "     "   .
       Q                       , &!REAL   , INTENT(INOUT) :: Q(nCols,rows)             ! INOUT Specific humidity (kg water/kg air).
       p_theta_levels          , &!REAL    ,INTENT(IN   ) :: p_theta_levels(nCols,rows)
       layer_thickness         , &!REAL    ,INTENT(IN   ) :: layer_thickness(nCols,rows)! IN thickness of layer (Pa)
       deltaz                  , &!REAL    ,INTENT(IN   ) :: deltaz(nCols,rows)        ! IN thickness of layer (m)
       rhodz_dry               , &!REAL    ,INTENT(IN   ) :: rhodz_dry(nCols,rows)     ! Dry air density layer thickness (kg m-2)
       rhodz_moist             , &!REAL    ,INTENT(IN   ) :: rhodz_moist(nCols,rows)   ! Moist air density layer thickness (kg m-2)
       nCols                   , &!INTEGER ,INTENT(IN   ) :: nCols    ! IN gather/scatter index
       RHCRIT                  , &!REAL    ,INTENT(IN   ) :: RHCRIT(nCols,rhc_rows)! IN Critical humidity for cloud formation.
       VFALL                   , &!REAL   , INTENT(INOUT) :: VFALL(nCols,rows)         ! INOUT fall velocity of ice (m per
       VFALL2                  , &!REAL   , INTENT(INOUT) :: VFALL2(nCols,rows)        ! INOUT fall vel. of rain (m/s)
       VFALL_RAIN              , &!REAL   , INTENT(INOUT) :: VFALL_RAIN(nCols,rows)    ! INOUT fall vel. of rain (m/s)
       VFALL_GRAUP             , &!REAL   , INTENT(INOUT) :: VFALL_GRAUP(nCols,rows)   ! INOUT fall vel. of graupel (m/s)
       FRAC_ICE_ABOVE          , &!REAL   , INTENT(INOUT) :: FRAC_ICE_ABOVE(nCols,rows)! INOUT Ice fraction from layer above water.
       CTTEMP                  , &!REAL   , INTENT(INOUT) :: CTTEMP(nCols,rows)        ! INOUT Ice cloud top temperature (K)
       RAINFRAC                , &!REAL   , INTENT(INOUT) :: RAINFRAC(nCols,rows)      ! INOUT Rain fraction.
       CX                      , &!REAL    ,INTENT(IN   ) ::  CX(80)
       CONSTP                  , &!REAL    ,INTENT(IN   ) ::  CONSTP(80)
       L_autoconv_murk         , &!LOGICAL ,INTENT(IN   ) :: L_autoconv_murk                ! Use murk aerosol to calc. drop number
       ec_auto                 , &!REAL    ,INTENT(IN   ) :: ec_auto                        ! Collision coalescence efficiency
       Ntot_land               , &!REAL    ,INTENT(IN   ) :: Ntot_land                      ! Number of droplets over land / m-3
       Ntot_sea                  &!REAL    ,INTENT(IN   ) :: Ntot_sea                       ! Number of droplets over sea / m-3
       )

    IMPLICIT NONE

    INTEGER ,INTENT(IN   ) :: N             ! IN Number of points where pptn non-zero from above or where CF>CFMIN
    INTEGER ,INTENT(IN   ) :: nCols    ! IN gather/scatter index
    INTEGER ,INTENT(IN   ) :: IX(nCols)! IN gather/scatter index
    INTEGER ,INTENT(IN   ) :: salt_dim_ice  ! Dimension to use for call to LSP_ICE.
    INTEGER ,INTENT(IN   ) :: biog_dim_ice  !     "     "   "   "   "   "     "   .
    ! 
    ! CF, CFL and CFF are IN/OUT if the PC2 cloud scheme is in use.
    !
    REAL(KIND=r8)    ,INTENT(INOUT) :: CF (nCols) ! IN Cloud fraction.
    REAL(KIND=r8)    ,INTENT(INOUT) :: CFL(nCols) ! IN Cloud liquid fraction.
    REAL(KIND=r8)    ,INTENT(INOUT) :: CFF(nCols) ! IN Cloud ice fraction.
    REAL(KIND=r8)    ,INTENT(IN   ) :: p_theta_levels(nCols)
    REAL(KIND=r8)    ,INTENT(IN   ) :: layer_thickness(nCols)! IN thickness of layer (Pa)
    REAL(KIND=r8)    ,INTENT(IN   ) :: deltaz(nCols)        ! IN thickness of layer (m)
    REAL(KIND=r8)    ,INTENT(IN   ) :: rhodz_dry(nCols)     ! Dry air density layer thickness (kg m-2)
    REAL(KIND=r8)    ,INTENT(IN   ) :: rhodz_moist(nCols)   ! Moist air density layer thickness (kg m-2)

    REAL(KIND=r8)    ,INTENT(IN   ) :: RHCRIT(nCols)! IN Critical humidity for cloud formation.
    REAL(KIND=r8)    ,INTENT(IN   ) :: TIMESTEP                       ! IN Timestep (sec).
    ! for conversion to ppn (kg water per m**3).

    LOGICAL , INTENT(IN  ) :: L_MURK                         ! IN Aerosol needs scavenging.
    LOGICAL , INTENT(IN  ) :: L_USE_SULPHATE_AUTOCONV        ! IN Sulphate aerosols used in
    ! autoconversion if T (i.e.
    ! second indirect effect)
    LOGICAL , INTENT(IN  ) :: L_SEASALT_CCN                  ! IN Sea-salt aerosols used for second indirect effect if T
    LOGICAL , INTENT(IN  ) :: L_BIOMASS_CCN                  ! IN Biomass smoke aerosols used for
    !                                                              ! second indirect effect if T
    LOGICAL , INTENT(IN  ) :: L_OCFF_CCN                     ! IN Fossil-fuel organic carbon
    ! aerosol used for second indirect effect if T
    LOGICAL , INTENT(IN  ) :: L_biogenic_CCN                 ! IN Biogenic aerosol used for second indirect effect if T
    LOGICAL , INTENT(IN  ) :: L_AUTO_DEBIAS                  ! IN Use autoconversion de-biasing
    LOGICAL , INTENT(IN  ) :: L_pc2                          ! IN Use PC2 cloud and condensation
    LOGICAL , INTENT(IN  ) :: L_mcr_qcf2                     ! IN Use second prognostic ice if T
    LOGICAL , INTENT(IN  ) :: L_it_melting                   ! IN Use iterative melting
    LOGICAL , INTENT(IN  ) :: L_non_hydrostatic              ! IN Use non-hydrostatic calculation of model layer depths
    LOGICAL , INTENT(IN  ) :: L_mixing_ratio                 ! IN q is a mixing ratio
    ! removed by deposition depending on amount of ice in each category
    !

    REAL(KIND=r8)   , INTENT(INOUT) :: Q(nCols)             ! INOUT Specific humidity (kg water/kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCF(nCols)           ! INOUT Cloud ice (kg per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCL(nCols)           ! INOUT Cloud liquid water (kg per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCF2(nCols)          ! INOUT Cloud ice2 (kg per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: T(nCols)             ! INOUT Temperature (K).
    REAL(KIND=r8)   , INTENT(INOUT) :: AEROSOL(nCols)       ! INOUT Aerosol (K).
    REAL(KIND=r8)   , INTENT(INOUT) :: LSRAIN(nCols)        ! INOUT Surface rainfall rate (kg m^-2 s^-1).
    REAL(KIND=r8)   , INTENT(INOUT) :: LSSNOW(nCols)        ! INOUT Surface snowfall rate (kg m^-2 s^-1).
    REAL(KIND=r8)   , INTENT(INOUT) :: LSSNOW2(nCols)       ! INOUT layer snowfall rate (kg m^-2 s^-1).
    REAL(KIND=r8)   , INTENT(INOUT) :: LSGRAUP(nCols)       ! INOUT layer graupelfall rate (kg m^-2 s^-1)
    REAL(KIND=r8)   , INTENT(INOUT) :: CTTEMP(nCols)        ! INOUT Ice cloud top temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: RAINFRAC(nCols)      ! INOUT Rain fraction.
    REAL(KIND=r8)   , INTENT(INOUT) :: FRAC_ICE_ABOVE(nCols)! INOUT Ice fraction from layer above water.
    REAL(KIND=r8)   , INTENT(INOUT) :: VFALL(nCols)         ! INOUT fall velocity of ice (m per
    REAL(KIND=r8)   , INTENT(INOUT) :: VFALL2(nCols)        ! INOUT fall vel. of rain (m/s)
    REAL(KIND=r8)   , INTENT(INOUT) :: VFALL_RAIN(nCols)    ! INOUT fall vel. of rain (m/s)
    REAL(KIND=r8)   , INTENT(INOUT) :: VFALL_GRAUP(nCols)   ! INOUT fall vel. of graupel (m/s)

    REAL(KIND=r8)   , INTENT(INOUT) :: SO4_AIT     (nCols) !INOUT S Cycle tracers & scavngd amoun
    REAL(KIND=r8)   , INTENT(INOUT) :: SO4_ACC     (nCols) !INOUT S Cycle tracers & scavngd amoun
    REAL(KIND=r8)   , INTENT(INOUT) :: SO4_DIS     (nCols) !INOUT S Cycle tracers & scavngd amoun
    REAL(KIND=r8)   , INTENT(INOUT) :: BMASS_AGD   (nCols) !INOUT S Cycle tracers & scavngd amoun
    REAL(KIND=r8)   , INTENT(INOUT) :: BMASS_CLD   (nCols) !INOUT S Cycle tracers & scavngd amoun
    REAL(KIND=r8)   , INTENT(INOUT) :: OCFF_AGD    (nCols) !INOUT S Cycle tracers & scavngd amoun
    REAL(KIND=r8)   , INTENT(INOUT) :: OCFF_CLD    (nCols) !INOUT S Cycle tracers & scavngd amoun


    !
    REAL(KIND=r8)    ,INTENT(IN   ) :: SNOW_DEPTH(nCols)
    REAL(KIND=r8)    ,INTENT(IN   ) :: LAND_FRACT(nCols)
    ! Sets up the size of arrays for CX and CONSTP
    REAL(KIND=r8)    ,INTENT(IN   ) ::  CX(80)
    REAL(KIND=r8)    ,INTENT(IN   ) ::  CONSTP(80)

    !
    REAL(KIND=r8)    ,INTENT(IN   ) :: SEA_SALT_FILM(nCols)
    REAL(KIND=r8)    ,INTENT(IN   ) :: SEA_SALT_JET (nCols)
    !
    REAL(KIND=r8)    ,INTENT(IN   ) :: biogenic(nCols)
    !
    LOGICAL ,INTENT(IN   ) :: L_autoconv_murk                ! Use murk aerosol to calc. drop number

    REAL(KIND=r8)    ,INTENT(IN   ) :: ec_auto                        ! Collision coalescence efficiency
    REAL(KIND=r8)    ,INTENT(IN   ) :: Ntot_land                      ! Number of droplets over land / m-3
    REAL(KIND=r8)    ,INTENT(IN   ) :: Ntot_sea                       ! Number of droplets over sea / m-3

    ! Microphysical process rate diagnostics (2D arrays on one level)
    ! Note: These arrays will only increase memory usage and are
    !       only referenced if the particular diagnostic is active



    !*L  Workspace usage ---------------------------------------------------
    !
    !ajm      REAL(KIND=r8)    ::      & PSTAR_C(N)         ! gathered Surface pressure (Pa).
    REAL(KIND=r8)    :: CF_C(N)                       ! gathered Cloud fraction.
    REAL(KIND=r8)    :: Q_C(N)                        ! gathered Specific humidity (kg water/kg air).
    REAL(KIND=r8)    :: QCF_C(N)                      ! gathered Cloud ice (kg per kg air).
    REAL(KIND=r8)    :: QCL_C(N)                      ! gathered Cloud liquid water (kg per kg air).
    REAL(KIND=r8)    :: QCF2_C(N)                     ! gathered cloud ice2 (kg per kg air).
    REAL(KIND=r8)    :: T_C(N)                        ! gathered Temperature (K).
    REAL(KIND=r8)    :: AERO_C(N)                     ! gathered Aerosol.
    REAL(KIND=r8)    :: LSRAIN_C(N)                   ! gathered Surface rainfall rate (kg per sq m per s).
    REAL(KIND=r8)    :: LSSNOW_C(N)                   ! gathered Surface snowfall rate (kg per sq m per s).
    REAL(KIND=r8)    :: LSSNOW2_C(N)                  ! gathered layer snowfall rate (kg per sq m per s).
    REAL(KIND=r8)    :: LSGRAUP_C(N)                  ! gathered layer graupel fall rate (kg/sq m/s)
    REAL(KIND=r8)    :: CTTEMP_C(N)                   ! gathered ice cloud top temperature.
    REAL(KIND=r8)    :: RAINFRAC_C(N)                 ! gathered rain fraction.
    REAL(KIND=r8)    :: FRAC_ICE_ABOVE_C(N)           ! gathered fraction of ice in layer above
    REAL(KIND=r8)    :: CFL_C(N)                      ! gathered Cloud liquid fraction.
    REAL(KIND=r8)    :: CFF_C(N)                      ! gathered Cloud ice fraction.
    REAL(KIND=r8)    :: VFALL_C(N)                    ! gathered fall velocity (m per s).
    REAL(KIND=r8)    :: VFALL2_C(N)                   ! gathered fall velocity for qcf2 (m per s).
    REAL(KIND=r8)    :: VFALL_RAIN_C(N)               ! gathered fall velocity for qcf2 (m per s).
    REAL(KIND=r8)    :: VFALL_GRAUP_C(N)              ! gathered fall velocity for qcf2 (m per s).
    REAL(KIND=r8)    :: SEA_SALT_FILM_C(salt_dim_ice) ! gathered film-mode sea-salt (m-3)
    REAL(KIND=r8)    :: SEA_SALT_JET_C (salt_dim_ice) ! gathered jet-mode sea-salt (m-3)
    REAL(KIND=r8)    :: biogenic_C(biog_dim_ice)      ! gathered biogenic aerosol (m.m.r.)
    REAL(KIND=r8)    :: RHC_C(N)                      ! gathered RH_crit value at points.

    REAL(KIND=r8)    :: SO4_AIT_C(N)                  ! gathered sulphate aerosol arrays
    REAL(KIND=r8)    :: SO4_ACC_C(N)                  ! gathered sulphate aerosol arrays
    REAL(KIND=r8)    :: SO4_DIS_C(N)                  ! gathered sulphate aerosol arrays
    !
    REAL(KIND=r8)    :: BMASS_AGD_C(N)                ! gathered biomass aerosol arrays
    REAL(KIND=r8)    :: BMASS_CLD_C(N)                ! gathered biomass aerosol arrays
    !
    REAL(KIND=r8)    :: OCFF_AGD_C(N)                 ! gathered fossil-fuel org carb aerosol arrays
    REAL(KIND=r8)    :: OCFF_CLD_C(N)                 ! gathered fossil-fuel org carb aerosol arrays

    !
    REAL(KIND=r8)    :: SNOW_DEPTH_C(N)
    REAL(KIND=r8)    :: LAND_FRACT_C(N)

    !
    REAL(KIND=r8)    :: RHODZ(N)                      ! WORK Used for air mass p.u.a. in successive layers.
    REAL(KIND=r8)    :: deltaz_c(n)                   ! Thickness of layer (m)
    REAL(KIND=r8)    :: rhodz_dry_c(n)                ! Dry air density * layer thickness (kg m-2)
    REAL(KIND=r8)    :: rhodz_moist_c(n)              ! Moist air density * layer thickness (kg m-2)

    REAL(KIND=r8)    :: P(N)                          ! WORK Used for pressure at successive levels.
    !
    ! Microphysical process rate diagnostics (compressed arrays)

    ! Call size of CX and CONSTP
    ! Start C_LSPSIZ
    ! Description: Include file containing idealised forcing options
    ! Author:      R. Forbes
    !
    ! History:
    ! Version  Date      Comment
    ! -------  ----      -------
    !   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
    !   6.2    22/08/05  Include the step size between ice categories.
    !                                                   Damian Wilson

    INTEGER,PARAMETER:: ice_type_offset=20

    ! End C_LSPSIZ
    ! Call comdeck containing ls ppn scavenging coeffs for Sulphur Cycle
    !
    !  External subroutines called -----------------------------------------
    !EXTERNAL LSP_ICE,LSP_SCAV
    !     &        ,SLSPSCV
    !*----------------------------------------------------------------------
    !  Physical constants -------------------------------------------------
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------
    REAL(KIND=r8) , PARAMETER :: P1UPONG=1.0_r8/G           ! One upon g (sq seconds per m).
    !  Define local variables ----------------------------------------------
    INTEGER :: MULTRHC         ! Zero if (nCols*rhc_rows) le 1, else 1
    !
    INTEGER :: I,II
    ! Loop counters: I - horizontal field index;
    INTEGER :: IRHI   !           : IRHI,IRHJ-indices for RHcrit.
    CF_C=0.0_r8;      Q_C=0.0_r8;     QCF_C=0.0_r8;    QCL_C=0.0_r8;    QCF2_C=0.0_r8;
    T_C=0.0_r8 ;      AERO_C=0.0_r8;  LSRAIN_C=0.0_r8; LSSNOW_C=0.0_r8; LSSNOW2_C=0.0_r8;
    LSGRAUP_C=0.0_r8; CTTEMP_C=0.0_r8;RAINFRAC_C=0.0_r8;FRAC_ICE_ABOVE_C=0.0_r8;    &
    CFL_C=0.0_r8;    CFF_C=0.0_r8;    VFALL_C=0.0_r8;
    VFALL2_C=0.0_r8;    VFALL_RAIN_C=0.0_r8;    VFALL_GRAUP_C=0.0_r8
    SEA_SALT_FILM_C=0.0_r8;    SEA_SALT_JET_C =0.0_r8;    biogenic_C=0.0_r8
    RHC_C=0.0_r8;    SO4_AIT_C=0.0_r8;    SO4_ACC_C=0.0_r8;
    SO4_DIS_C=0.0_r8;        BMASS_AGD_C=0.0_r8;    BMASS_CLD_C=0.0_r8;
    OCFF_AGD_C=0.0_r8;    OCFF_CLD_C=0.0_r8;
    SNOW_DEPTH_C=0.0_r8;    LAND_FRACT_C=0.0_r8;
    RHODZ=0.0_r8;    deltaz_c=0.0_r8;    rhodz_dry_c=0.0_r8;    rhodz_moist_c=0.0_r8;
    P=0.0_r8;                          
    !
    !-----------------------------------------------------------------------
    !L Internal structure.
    !L 1. gather variables using index
    !-----------------------------------------------------------------------
    IF ( (nCols )  >   1) THEN
       MULTRHC = 1
    ELSE
       MULTRHC = 0
    END IF
    !
    DO I=1,N
       II=IX(I)
       ! IJ=IX(I,2)
       IRHI = (MULTRHC * (II - 1)) + 1
       !IRHJ = (MULTRHC * (IJ - 1)) + 1
       LSRAIN_C(I)        = LSRAIN(II)
       LSSNOW_C(I)        = LSSNOW(II)
       LSSNOW2_C(I)       = LSSNOW2(II)
       LSGRAUP_C(I)       = LSGRAUP(II)
       CTTEMP_C(I)        = CTTEMP(II)
       RAINFRAC_C(I)      = RAINFRAC(II)
       FRAC_ICE_ABOVE_C(I)= FRAC_ICE_ABOVE(II)
       !ajm        PSTAR_C(I) =PSTAR(II)
       P(I)               = p_theta_levels(II)
       RHODZ(I)           = -P1UPONG*layer_thickness(II)
       deltaz_c(i)        = deltaz(ii)
       rhodz_dry_c(i)     = rhodz_dry(ii)
       rhodz_moist_c(i)   = rhodz_moist(ii)
       CF_C (I)           = CF(II)
       CFL_C(I)           = CFL(II)
       CFF_C(I)           = CFF(II)
       QCF_C(I)           = QCF(II)
       QCL_C(I)           = QCL(II)
       Q_C(I)=Q(II)
       T_C(I)=T(II)
       IF (L_mcr_qcf2) QCF2_C(I) = QCF2(II)
       IF (L_MURK) AERO_C(I)=AEROSOL(II)
       VFALL_C(I)=VFALL(II)
       VFALL2_C(I) = VFALL2(II)
       VFALL_RAIN_C(I) = VFALL_RAIN(II)
       VFALL_GRAUP_C(I) = VFALL_GRAUP(II)
       RHC_C(I)=RHCRIT(IRHI)
       IF (L_USE_SULPHATE_AUTOCONV) THEN
          SO4_AIT_C(I)=SO4_AIT(II)
          SO4_ACC_C(I)=SO4_ACC(II)
          SO4_DIS_C(I)=SO4_DIS(II)
       END IF
       IF (L_SEASALT_CCN) THEN
          SEA_SALT_FILM_C(I)=SEA_SALT_FILM(II)
          SEA_SALT_JET_C(I)=SEA_SALT_JET(II)
       ELSE
          SEA_SALT_FILM_C(1)=SEA_SALT_FILM(1)
          SEA_SALT_JET_C(1)=SEA_SALT_JET(1)
       ENDIF
       IF (L_BIOMASS_CCN) THEN
          BMASS_AGD_C(I)=BMASS_AGD(II)
          BMASS_CLD_C(I)=BMASS_CLD(II)
       ELSE
          BMASS_AGD_C(I)=BMASS_AGD(1)
          BMASS_CLD_C(I)=BMASS_CLD(1)
       ENDIF
       IF (L_biogenic_CCN) THEN
          biogenic_C(I)=biogenic(II)
       ELSE
          biogenic_C(1)=biogenic(1)
       ENDIF
       IF (L_OCFF_CCN) THEN
          OCFF_AGD_C(I)=OCFF_AGD(II)
          OCFF_CLD_C(I)=OCFF_CLD(II)
       ELSE
          OCFF_AGD_C(I)=OCFF_AGD(1)
          OCFF_CLD_C(I)=OCFF_CLD(1)
       ENDIF
       SNOW_DEPTH_C(I)=SNOW_DEPTH(II)
       LAND_FRACT_C(I)=LAND_FRACT(II)

       ! Process diagnostic arrays are initialised to zero in LSP_ICE
    END DO ! Loop over points
    !
    !-----------------------------------------------------------------------
    ! ICE FORMATION/EVAPORATION/MELTING
    ! WATER CLOUD AND RAIN FORMATION/EVAPORATION
    !-----------------------------------------------------------------------
    ! The call to LSP_ICE replaces the calls to LSP_EVAP, LSPFRMT
    ! and LSP_FORM.
    ! CFL_C contains cloud fraction for ice
    ! CFF_C contains cloud fraction for water
    ! DEPENDS ON: lsp_ice
    CALL LSP_ICE( &
         P                       , &!REAL(KIND=r8)   , INTENT(IN   ) :: P(POINTS)          ! Air pressure at this level (Pa).
         RHODZ                   , &!REAL   , INTENT(IN   ) :: RHODZ(POINTS)      ! Air mass p.u.a. in this layer (kg per sq m).
         deltaz_c                , &!REAL   , INTENT(IN   ) :: deltaz(points)     ! Depth of layer / m
         rhodz_dry_c             , &!REAL   , INTENT(IN   ) :: rhodz_dry(points)  ! Density of dry air / kg m-3
         rhodz_moist_c           , &!REAL   , INTENT(IN   ) :: rhodz_moist(points)! Density of moist air / kg m-3
         TIMESTEP                , &!REAL   , INTENT(IN   ) :: TIMESTEPFIXED      ! Timestep of physics in model (s).
         N                       , &!INTEGER, INTENT(IN   ) :: POINTS       ! Number of points to be processed.
         L_MURK                  , &!LOGICAL, INTENT(IN   ) :: L_MURK       ! Murk aerosol is a valid quantity
         RHC_C                   , &!REAL   , INTENT(IN   ) :: RHCPT(POINTS)      ! Critical relative humidity of all points for cloud formation.
         L_USE_SULPHATE_AUTOCONV , &!LOGICAL, INTENT(IN   ) :: L_USE_SULPHATE_AUTOCONV!Switch to use sulphate aerosol in the calculation of cloud
         L_AUTO_DEBIAS           , &!LOGICAL, INTENT(IN   ) :: L_AUTO_DEBIAS  ! Switch to apply de-biasing correction to autoconversion rate.
         L_mcr_qcf2              , &!LOGICAL, INTENT(IN   ) :: L_mcr_qcf2   ! true if using 2nd cloud ice prognostic
         L_it_melting            , &!LOGICAL, INTENT(IN   ) :: L_it_melting           ! true if using iterative melting
         l_non_hydrostatic       , &!LOGICAL, INTENT(IN   ) :: L_non_hydrostatic ! Use non hydrostatic layer masses
         l_mixing_ratio          , &!LOGICAL, INTENT(IN   ) :: L_mixing_ratio    ! Use mixing ratio formulation
         SO4_ACC_C               , &!REAL   , INTENT(IN   ) :: SO4_ACC(POINTS)    ! Sulphur cycle variable
         SO4_DIS_C               , &!REAL   , INTENT(IN   ) :: SO4_DIS(POINTS)    ! Sulphur cycle variable
         !SO4_AIT_C               , &!REAL   , INTENT(IN   ) :: SO4_AIT(POINTS)    ! Sulphur cycle variable
         L_BIOMASS_CCN           , &!LOGICAL, INTENT(IN   ) :: L_BIOMASS_CCN! Switch to supplement sulphate aerosol with biomass smoke
         BMASS_AGD_C             , &!REAL   , INTENT(IN   ) :: BMASS_AGD(POINTS)  ! Aged biomass smoke mass mixing ratio
         BMASS_CLD_C             , &!REAL   , INTENT(IN   ) :: BMASS_CLD(POINTS)  ! In-cloud biomass smoke mass mixing ratio
         L_OCFF_CCN              , &!LOGICAL, INTENT(IN   ) :: L_OCFF_CCN   ! Switch to supplement sulphate aerosol with fossil-fuel organic
         OCFF_AGD_C              , &!REAL   , INTENT(IN   ) :: OCFF_AGD(POINTS)   ! Aged fossil-fuel organic carbon mass mixing ratio
         OCFF_CLD_C              , &!REAL   , INTENT(IN   ) :: OCFF_CLD(POINTS)   ! In-cloud fossil-fuel organic carbon mass mixing ratio
         L_SEASALT_CCN           , &!LOGICAL, INTENT(IN   ) :: L_SEASALT_CCN! Switch to supplement sulphate aerosol with sea-salt if the
         SEA_SALT_FILM_C         , &!REAL   , INTENT(IN   ) :: SEA_SALT_FILM(SALT_DIM_ICE)! Film-mode sea-salt aerosol number concentration (m-3)
         SEA_SALT_JET_C          , &!REAL   , INTENT(IN   ) :: SEA_SALT_JET(SALT_DIM_ICE) ! Jet-mode sea-salt aerosol number concentration (m-3)
         salt_dim_ice            , &!INTEGER, INTENT(IN   ) :: SALT_DIM_ICE ! Number of points for sea-salt arrays (either POINTS or 1)
         L_biogenic_CCN          , &!LOGICAL, INTENT(IN   ) :: L_biogenic_CCN !Switch to supplement sulphate aerosol with biogenic aerosol
         biogenic_C              , &!REAL   , INTENT(IN   ) :: biogenic(biog_dim_ice)     ! Biogenic aerosol m.m.r.
         biog_dim_ice            , &!INTEGER, INTENT(IN   ) :: biog_dim_ice ! Number of points for biogenic array  (either POINTS or 1)
         SNOW_DEPTH_C            , &!REAL   , INTENT(IN   ) :: SNOW_DEPTH(POINTS) ! Snow depth for aerosol amount (m)
         LAND_FRACT_C            , &!REAL   , INTENT(IN   ) :: LAND_FRACT(POINTS) ! Land fraction
         AERO_C                  , &!REAL   , INTENT(IN   ) :: AEROSOL(POINTS)    ! Aerosol mass (ug/kg)
         L_autoconv_murk         , &!LOGICAL, INTENT(IN   ) :: L_autoconv_murk! Use murk aerosol to calc. drop number
         ec_auto                 , &!Real   , INTENT(IN   ) :: ec_auto       ! Collision coalescence efficiency
         Ntot_land               , &!Real   , INTENT(IN   ) :: Ntot_land     ! Number of droplets over land / m-3
         Ntot_sea                , &!Real   , INTENT(IN   ) :: Ntot_sea      ! Number of droplets over sea  / m-3
         QCF_C                   , &!REAL   , INTENT(INOUT) :: QCF(POINTS)   ! Cloud ice (kg water per kg air).
         QCL_C                   , &!REAL   , INTENT(INOUT) :: QCL(POINTS)   ! Cloud liquid water (kg water per kg air).
         Q_C                     , &!REAL   , INTENT(INOUT) :: Q(POINTS)     ! Specific humidity at this level (kg water per kg air).
         QCF2_C                  , &!REAL   , INTENT(INOUT) :: QCF2(POINTS)  ! Second cloud ice (kg water per kg air).
         LSRAIN_C                , &!REAL   , INTENT(INOUT) :: RAIN(POINTS)  ! On input: Rate of rainfall entering this layer from above.
         LSSNOW_C                , &!REAL   , INTENT(INOUT) :: SNOW_AGG(POINTS) ! Dummy in 3C version of the code
         VFALL_C                 , &!REAL   , INTENT(INOUT) :: VF(POINTS)      ! On input: Fall velocity of snow aggregates entering layer.
         LSSNOW2_C               , &!REAL   , INTENT(INOUT) :: SNOW_CRY(POINTS)! On input: Rate of snow crystals entering layer from above.
         VFALL2_C                , &!REAL   , INTENT(INOUT) :: VF_CRY(POINTS)  ! On input: Fall velocity of snow crystals entering this layer.
         FRAC_ICE_ABOVE_C        , &!REAL   , INTENT(INOUT) :: FRAC_ICE_ABOVE(POINTS)! Fraction of ice in layer above (no units)
         CTTEMP_C                , &!REAL   , INTENT(INOUT) :: CTTEMP(POINTS)        ! Ice cloud top temperature (K)
         RAINFRAC_C              , &!REAL   , INTENT(INOUT) :: RAINFRAC(POINTS)      ! Rain fraction (no units)
         T_C                     , &!REAL   , INTENT(INOUT) :: T(POINTS)     ! Temperature at this level (K).
         CF_C                    , &!REAL   , INTENT(IN   ) :: CFKEEP(POINTS)     ! Total cloud fraction in this layer (no units).
         CFL_C                   , &!REAL   , INTENT(IN   ) :: CFLIQ(POINTS)      ! Liquid cloud fraction in this layer (no units).
         CFF_C                   , &!REAL   , INTENT(IN   ) :: CFICEKEEP(POINTS)  ! Frozen cloud fraction in this layer (no units).
         CX                      , &!REAL   , INTENT(IN   ) :: CX(80)
         CONSTP                    &!REAL   , INTENT(IN   ) :: CONSTP(80)
         )
    !-----------------------------------------------------------------------
    !L 3.4 Lose aerosol by scavenging: call LSP_SCAV
    !-----------------------------------------------------------------------
    !
    IF (L_MURK)  THEN
       ! DEPENDS ON: lsp_scav
       CALL LSP_SCAV( &
            TIMESTEP, & !REAL   , INTENT(IN   ) :: TIMESTEP! Input real scalar :- IN Timestep (s).
            N       , & !INTEGER, INTENT(IN   ) :: POINTS  ! Input integer scalar :- IN Number of points to be processed.
            LSRAIN_C, & !REAL   , INTENT(IN   ) :: RAIN(POINTS)  ! Input real arrays :- IN Rate of rainfall in this layer from above! (kg per sq m per s).
            LSSNOW_C, & !REAL   , INTENT(IN   ) :: SNOW(POINTS)  ! IN Rate of snowfall in this layer from above (kg per sq m per s).
            AERO_C    ) !REAL   , INTENT(INOUT) :: AEROSOL(POINTS) ! Updated real arrays :-INOUT Aerosol mixing ratio
    ENDIF
    !
    !-----------------------------------------------------------------------
    !L 4  Scatter back arrays which will have been changed.
    !L
    !-----------------------------------------------------------------------
    !

    DO I=1,N
       II=IX(I)
       !IJ=IX(I,2)
       T(II)=T_C(I)
       Q(II)=Q_C(I)
       QCF(II)=QCF_C(I)
       QCL(II)=QCL_C(I)
       IF (L_mcr_qcf2) QCF2(II) = QCF2_C(I)
       IF (L_MURK) AEROSOL(II)=AERO_C(I)
       LSRAIN(II)=LSRAIN_C(I)
       LSSNOW(II)=LSSNOW_C(I)
       LSSNOW2(II)=LSSNOW2_C(I)
       LSGRAUP(II)=LSGRAUP_C(I)
       CTTEMP(II)=CTTEMP_C(I)
       RAINFRAC(II)=RAINFRAC_C(I)
       FRAC_ICE_ABOVE(II)=FRAC_ICE_ABOVE_C(I)
       IF (L_pc2) THEN
          CFF(II)=CFF_C(I)
          CFL(II)=CFL_C(I)
          CF (II)=CF_C(I)
       END IF  ! L_pc2

       VFALL(II)=VFALL_C(I)
       VFALL2(II)=VFALL2_C(I)
       VFALL_RAIN(II)=VFALL_RAIN_C(I)
       VFALL_GRAUP(II)=VFALL_GRAUP_C(I)
       ! Only store process rates in array for diagnostic
       ! if a particular diagnostic is requested,
       ! otherwise overwriting will occur
       ! (space for the 3D array in MCR_CTL is only allocated
       ! if the diagnostic is active, to save memory)
    END DO ! Loop over points
    !
    !
    RETURN
  END SUBROUTINE LS_PPNC



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !+ Precipitation microphysics calculations.
  ! Subroutine Interface:
  SUBROUTINE LSP_ICE(&
       P                      , &!REAL   , INTENT(IN   ) :: P(POINTS)          ! Air pressure at this level (Pa).
       RHODZ                  , &!REAL   , INTENT(IN   ) :: RHODZ(POINTS)      ! Air mass p.u.a. in this layer (kg per sq m).
       deltaz                 , &!REAL   , INTENT(IN   ) :: deltaz(points)     ! Depth of layer / m
       rhodz_dry              , &!REAL   , INTENT(IN   ) :: rhodz_dry(points)  ! Density of dry air / kg m-3
       rhodz_moist            , &!REAL   , INTENT(IN   ) :: rhodz_moist(points)! Density of moist air / kg m-3
       TIMESTEPFIXED          , &!REAL   , INTENT(IN   ) :: TIMESTEPFIXED      ! Timestep of physics in model (s).
       POINTS                 , &!INTEGER, INTENT(IN   ) :: POINTS       ! Number of points to be processed.
       L_MURK                 , &!LOGICAL, INTENT(IN   ) :: L_MURK       ! Murk aerosol is a valid quantity
       RHCPT                  , &!REAL   , INTENT(IN   ) :: RHCPT(POINTS)      ! Critical relative humidity of all points for cloud formation.
       L_USE_SULPHATE_AUTOCONV, &!LOGICAL, INTENT(IN   ) :: L_USE_SULPHATE_AUTOCONV!Switch to use sulphate aerosol in the calculation of cloud
       L_AUTO_DEBIAS          , &!LOGICAL, INTENT(IN   ) :: L_AUTO_DEBIAS  ! Switch to apply de-biasing correction to autoconversion rate.
       L_mcr_qcf2             , &!LOGICAL, INTENT(IN   ) :: L_mcr_qcf2   ! true if using 2nd cloud ice prognostic
       L_it_melting           , &!LOGICAL, INTENT(IN   ) :: L_it_melting      ! true if using iterative melting
       l_non_hydrostatic      , &!LOGICAL, INTENT(IN   ) :: L_non_hydrostatic ! Use non hydrostatic layer masses
       l_mixing_ratio         , &!LOGICAL, INTENT(IN   ) :: L_mixing_ratio    ! Use mixing ratio formulation
       SO4_ACC                , &!REAL   , INTENT(IN   ) :: SO4_ACC(POINTS)    ! Sulphur cycle variable
       SO4_DIS                , &!REAL   , INTENT(IN   ) :: SO4_DIS(POINTS)    ! Sulphur cycle variable
!       SO4_AIT                , &!REAL   , INTENT(IN   ) :: SO4_AIT(POINTS)    ! Sulphur cycle variable
       L_BIOMASS_CCN          , &!LOGICAL, INTENT(IN   ) :: L_BIOMASS_CCN! Switch to supplement sulphate aerosol with biomass smoke
       BMASS_AGD              , &!REAL   , INTENT(IN   ) :: BMASS_AGD(POINTS)  ! Aged biomass smoke mass mixing ratio
       BMASS_CLD              , &!REAL   , INTENT(IN   ) :: BMASS_CLD(POINTS)  ! In-cloud biomass smoke mass mixing ratio
       L_OCFF_CCN             , &!LOGICAL, INTENT(IN   ) :: L_OCFF_CCN   ! Switch to supplement sulphate aerosol with fossil-fuel organic
       OCFF_AGD               , &!REAL   , INTENT(IN   ) :: OCFF_AGD(POINTS)   ! Aged fossil-fuel organic carbon mass mixing ratio
       OCFF_CLD               , &!REAL   , INTENT(IN   ) :: OCFF_CLD(POINTS)   ! In-cloud fossil-fuel organic carbon mass mixing ratio
       L_SEASALT_CCN          , &!LOGICAL, INTENT(IN   ) :: L_SEASALT_CCN! Switch to supplement sulphate aerosol with sea-salt if the
       SEA_SALT_FILM          , &!REAL   , INTENT(IN   ) :: SEA_SALT_FILM(SALT_DIM_ICE)! Film-mode sea-salt aerosol number concentration (m-3)
       SEA_SALT_JET           , &!REAL   , INTENT(IN   ) :: SEA_SALT_JET(SALT_DIM_ICE) ! Jet-mode sea-salt aerosol number concentration (m-3)
       SALT_DIM_ICE           , &!INTEGER, INTENT(IN   ) :: SALT_DIM_ICE ! Number of points for sea-salt arrays (either POINTS or 1)
       L_biogenic_CCN         , &!LOGICAL, INTENT(IN   ) :: L_biogenic_CCN !Switch to supplement sulphate aerosol with biogenic aerosol
       biogenic               , &!REAL   , INTENT(IN   ) :: biogenic(biog_dim_ice)     ! Biogenic aerosol m.m.r.
       biog_dim_ice           , &!INTEGER, INTENT(IN   ) :: biog_dim_ice ! Number of points for biogenic array  (either POINTS or 1)
       SNOW_DEPTH             , &!REAL   , INTENT(IN   ) :: SNOW_DEPTH(POINTS) ! Snow depth for aerosol amount (m)
       LAND_FRACT             , &!REAL   , INTENT(IN   ) :: LAND_FRACT(POINTS) ! Land fraction
       AEROSOL                , &!REAL   , INTENT(IN   ) :: AEROSOL(POINTS)    ! Aerosol mass (ug/kg)
       L_autoconv_murk        , &!LOGICAL, INTENT(IN   ) :: L_autoconv_murk! Use murk aerosol to calc. drop number
       ec_auto                , &!Real   , INTENT(IN   ) :: ec_auto       ! Collision coalescence efficiency
       Ntot_land              , &!Real   , INTENT(IN   ) :: Ntot_land     ! Number of droplets over land / m-3
       Ntot_sea               , &!Real   , INTENT(IN   ) :: Ntot_sea      ! Number of droplets over sea  / m-3
       QCF                    , &!REAL   , INTENT(INOUT) :: QCF(POINTS)   !         Cloud ice (kg water per kg air).
       QCL                    , &!REAL   , INTENT(INOUT) :: QCL(POINTS)   !         Cloud liquid water (kg water per kg air).
       Q                      , &!REAL   , INTENT(INOUT) :: Q(POINTS)     !         Specific humidity at this level (kg water per kg air).
       QCF2                   , &!REAL   , INTENT(INOUT) :: QCF2(POINTS)  !         Second cloud ice (kg water per kg air).
       RAIN                   , &!REAL   , INTENT(INOUT) :: RAIN(POINTS)  !         On input: Rate of rainfall entering this layer from above.
       SNOW_AGG               , &!REAL   , INTENT(INOUT) :: SNOW_AGG (POINTS) ! Dummy in 3C version of the code
       VF                     , &!REAL   , INTENT(INOUT) :: VF(POINTS)      ! On input: Fall velocity of snow aggregates entering layer.
       SNOW_CRY               , &!REAL   , INTENT(INOUT) :: SNOW_CRY(POINTS)! On input: Rate of snow crystals entering layer from above.
       VF_CRY                 , &!REAL   , INTENT(INOUT) :: VF_CRY(POINTS)  ! On input: Fall velocity of snow crystals entering this layer.
       FRAC_ICE_ABOVE         , &!REAL   , INTENT(INOUT) :: FRAC_ICE_ABOVE(POINTS)! Fraction of ice in layer above (no units)
       CTTEMP                 , &!REAL   , INTENT(INOUT) :: CTTEMP(POINTS)        ! Ice cloud top temperature (K)
       RAINFRAC               , &!REAL   , INTENT(INOUT) :: RAINFRAC(POINTS)      ! Rain fraction (no units)
       T                      , &!REAL   , INTENT(INOUT) :: T(POINTS)     !         Temperature at this level (K).
       CFKEEP                 , &!REAL   , INTENT(IN   ) :: CFKEEP(POINTS)     ! Total cloud fraction in this layer (no units).
       CFLIQ                  , &!REAL   , INTENT(IN   ) :: CFLIQ(POINTS)      ! Liquid cloud fraction in this layer (no units).
       CFICEKEEP              , &!REAL   , INTENT(IN   ) :: CFICEKEEP(POINTS)  ! Frozen cloud fraction in this layer (no units).
       CX                     , &!REAL   , INTENT(IN   ) :: CX(80)
       CONSTP                   &!REAL   , INTENT(IN   ) :: CONSTP(80)
       )


    IMPLICIT NONE
    !
    ! Description:
    !   Updates ice, liquid and vapour contents and temperature as a
    !   result of microphysical processes.
    !
    ! Method:
    !   Calculates transfers of water between vapour, ice, cloud liquid
    !   and rain. Advects ice downwards. Processes included are:
    !   Fall of ice into and out of the layer;
    !   Homogenous and heterogenous nucleation of ice;
    !   Deposition and sublimation of ice;
    !   Riming; riming by supercooled raindrops;
    !   Melting of ice; Evaporation of rain; accretion;
    !   Autoconversion of liquid water to rain.
    !   This is described in Unified Model Documentation Paper 26.
    !
    !   Microphysics options:
    !   - Second prognostic cloud ice variables
    !      Active if L_mcr_qcf2=.True.
    !      The code supports the use of a second cloud ice prognostic
    !      variable so that both cloud ice aggregates (QCF and QCF_AGG)
    !      and cloud ice pristine crystals (QCF2 and QCF_CRY) can be
    !      represented and advected.
    !      At UM5.5 this code is still experimental.
    !   - Prognostic rain (controlled by L_mcr_qrain)
    !     At UM5.5 no microphysical process/sedimentation code is present.
    !   - Prognostic graupel (controlled by L_mcr_qrain)
    !     At UM5.5 no microphysical process/sedimentation code is present.
    !
    ! Current Code Owner: Damian Wilson
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    ! 5.2       21/11/00 Original code. Damian Wilson
    !LL  5.3  24/09/01  Portability changes.    Z. Gardner
    ! 5.3       06/08/01 Modify fall out algorithm to be more numerically
    !                    stable. Adjust the layer density calculation to
    !                    be more numerically robust.  Damian Wilson
    !           03/12/01 Correct wet bulb temperature calculation.
    !                                                 Damian Wilson
    !
    !
    ! 5.3       29/08/01 Introduce logical to control effect of sulphate
    !                    aerosol on autoconversion.  Andy Jones
    ! 5.4       02/09/02 Include c_micro for homogenous freezing
    !                    temperature                 Damian Wilson
    !
    ! 5.4       25/04/02 Replace land/sea mask with land fraction field
    !                    in call to NUMBER_DROPLET.
    !                                             Andy Jones
    ! 5.5       03/02/03 Include second ice, rain and graupel prognostic
    !                    variables.                            R.M.Forbes
    ! 5.5       24/02/03 Correct sign error in evaporation of melting snow
    !                    term.                    Damian Wilson
    ! 6.1       01/08/04 Include dummy variables in argument list. R.Forbes
    ! 6.1       30/07/04 Alter the subgrid-scale treatment in the
    !                    deposition/sublimation term.  Damian Wilson
    ! 6.1       05/05/04 Correct evaporation of rain error in the
    !                    non VECTLIB branch of the code. Damian Wilson
    ! 6.1       07/04/04 Add biomass smoke to call to NUMBER_DROPLET.
    !                                             Andy Jones
    ! 6.1       07/04/04 Add option to apply bias-removal scheme to the
    !                    autoconversion.          Andy Jones
    ! 6.2       02/11/05 Provide a minimum limit to width of PDF
    !                    in deposition term.  Damian Wilson
    ! 6.2       02/11/05 Add iterative melting    Damian Wilson
    ! 6.2       23/11/05 Pass through precip variables from UMUI. D. Wilson
    ! 6.2       12/01/06 Pass in dummy variables for non-hydrostatic code
    !                                                      Damian Wilson
    ! 6.2       20/10/05 Update comments for bias-removal scheme.
    !                                             Andy Jones
    ! 6.2       31/01/06 Pass in dummy rate diagnostics. Richard Forbes
    ! 6.4       10/01/07 Include mixing ratio control logical. D. Wilson
    !
    ! Code Description:
    !   Language: FORTRAN 77 + common extensions.
    !   This code is written to UMDP3 v6 programming standards.
    !
    INTEGER, INTENT(IN   ) :: POINTS       ! Number of points to be processed.
    INTEGER, INTENT(IN   ) :: SALT_DIM_ICE ! Number of points for sea-salt arrays (either POINTS or 1)
    INTEGER, INTENT(IN   ) :: biog_dim_ice ! Number of points for biogenic array  (either POINTS or 1)
    !
    LOGICAL, INTENT(IN   ) :: L_MURK       ! Murk aerosol is a valid quantity
    LOGICAL, INTENT(IN   ) :: L_USE_SULPHATE_AUTOCONV!Switch to use sulphate aerosol in the calculation of cloud
    !droplet number concentration in the autoconversion section,
    !i.e. activate the second indirect ("lifetime") effect.
    LOGICAL, INTENT(IN   ) :: L_SEASALT_CCN! Switch to supplement sulphate aerosol with sea-salt if the
    ! second indirect effect is active.
    LOGICAL, INTENT(IN   ) :: L_BIOMASS_CCN! Switch to supplement sulphate aerosol with biomass smoke
    ! aerosol if the second indirect effect is active.
    LOGICAL, INTENT(IN   ) :: L_OCFF_CCN   ! Switch to supplement sulphate aerosol with fossil-fuel organic
    ! carbon aerosol if the second indirect effect is active.
    LOGICAL, INTENT(IN   ) :: L_biogenic_CCN !Switch to supplement sulphate aerosol with biogenic aerosol
    ! if the second indirect effect is active.
    LOGICAL, INTENT(IN   ) :: L_AUTO_DEBIAS  ! Switch to apply de-biasing correction to autoconversion rate.
    LOGICAL, INTENT(IN   ) :: L_mcr_qcf2        ! true if using 2nd cloud ice prognostic
    LOGICAL, INTENT(IN   ) :: L_it_melting      ! true if using iterative melting
    LOGICAL, INTENT(IN   ) :: L_non_hydrostatic ! Use non hydrostatic layer masses
    LOGICAL, INTENT(IN   ) :: L_mixing_ratio    ! Use mixing ratio formulation
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: TIMESTEPFIXED      ! Timestep of physics in model (s).
    REAL(KIND=r8)   , INTENT(IN   ) :: CFLIQ(POINTS)      ! Liquid cloud fraction in this layer (no units).
    REAL(KIND=r8)   , INTENT(IN   ) :: CFICEKEEP(POINTS)  ! Frozen cloud fraction in this layer (no units).
    REAL(KIND=r8)   , INTENT(IN   ) :: CFKEEP(POINTS)     ! Total cloud fraction in this layer (no units).
    REAL(KIND=r8)   , INTENT(IN   ) :: P(POINTS)          ! Air pressure at this level (Pa).
    REAL(KIND=r8)   , INTENT(IN   ) :: RHODZ(POINTS)      ! Air mass p.u.a. in this layer (kg per sq m).
    REAL(KIND=r8)   , INTENT(IN   ) :: deltaz(points)     ! Depth of layer / m
    REAL(KIND=r8)   , INTENT(IN   ) :: rhodz_dry(points)  ! Density of dry air / kg m-3
    REAL(KIND=r8)   , INTENT(IN   ) :: rhodz_moist(points)! Density of moist air / kg m-3
    REAL(KIND=r8)   , INTENT(IN   ) :: SO4_ACC(POINTS)    ! Sulphur cycle variable
    REAL(KIND=r8)   , INTENT(IN   ) :: SO4_DIS(POINTS)    ! Sulphur cycle variable
    !REAL(KIND=r8)   , INTENT(IN   ) :: SO4_AIT(POINTS)    ! Sulphur cycle variable
    REAL(KIND=r8)   , INTENT(IN   ) :: BMASS_AGD(POINTS)  ! Aged biomass smoke mass mixing ratio
    REAL(KIND=r8)   , INTENT(IN   ) :: BMASS_CLD(POINTS)  ! In-cloud biomass smoke mass mixing ratio
    REAL(KIND=r8)   , INTENT(IN   ) :: OCFF_AGD(POINTS)   ! Aged fossil-fuel organic carbon mass mixing ratio
    REAL(KIND=r8)   , INTENT(IN   ) :: OCFF_CLD(POINTS)   ! In-cloud fossil-fuel organic carbon mass mixing ratio
    REAL(KIND=r8)   , INTENT(IN   ) :: SNOW_DEPTH(POINTS) ! Snow depth for aerosol amount (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: LAND_FRACT(POINTS) ! Land fraction
    REAL(KIND=r8)   , INTENT(IN   ) :: AEROSOL(POINTS)    ! Aerosol mass (ug/kg)
    REAL(KIND=r8)   , INTENT(IN   ) :: RHCPT(POINTS)      ! Critical relative humidity of all points for cloud formation.
    REAL(KIND=r8)   , INTENT(IN   ) :: SEA_SALT_FILM(SALT_DIM_ICE)! Film-mode sea-salt aerosol number concentration (m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: SEA_SALT_JET(SALT_DIM_ICE) ! Jet-mode sea-salt aerosol number concentration (m-3)
    REAL(KIND=r8)   , INTENT(IN   ) :: biogenic(biog_dim_ice)     ! Biogenic aerosol m.m.r.
    !
    ! End C_LSPSIZ
    !
    REAL(KIND=r8)   , INTENT(INOUT) :: Q(POINTS)     !         Specific humidity at this level (kg water per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCF(POINTS)   !         Cloud ice (kg water per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCL(POINTS)   !         Cloud liquid water (kg water per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: QCF2(POINTS)  !         Second cloud ice (kg water per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) :: T(POINTS)     !         Temperature at this level (K).
    REAL(KIND=r8)   , INTENT(INOUT) :: RAIN(POINTS)  !         On input: Rate of rainfall entering this layer from above.
    !         On output: Rate of rainfall leaving this layer.
    !                   (kg m-2 s-1).
    REAL(KIND=r8)   , INTENT(INOUT) :: SNOW_AGG(POINTS)! On input: Rate of snow aggregates entering layer from above.
    ! On output: Rate of snow aggregates leaving this layer.
    ! (kg m-2 s-1). If only one ice prognostic is active
    ! then this variable contains all the snow.
    REAL(KIND=r8)   , INTENT(INOUT) :: SNOW_CRY(POINTS)! On input: Rate of snow crystals entering layer from above.
    ! On Output: Rate of snow crystals leaving this layer.
    ! (kg m-2 s-1). Only non-zero if two ice
    ! prognostics in use.
    REAL(KIND=r8)   , INTENT(INOUT) :: VF(POINTS)      ! On input: Fall velocity of snow aggregates entering layer.
    ! On Output: Fall velocity of snow aggregates leaving layer.
    ! (m s-1). If only one ice prognostic is active
    ! then this is the velocity of all falling snow.
    REAL(KIND=r8)   , INTENT(INOUT) :: VF_CRY(POINTS)  ! On input: Fall velocity of snow crystals entering this layer.
    ! On Output: Fall velocity of snow crystals leaving this layer.
    ! (m s-1). Only used if two ice prognostics in use.
    REAL(KIND=r8)   , INTENT(INOUT) :: CTTEMP(POINTS)        ! Ice cloud top temperature (K)
    REAL(KIND=r8)   , INTENT(INOUT) :: RAINFRAC(POINTS)      ! Rain fraction (no units)
    REAL(KIND=r8)   , INTENT(INOUT) :: FRAC_ICE_ABOVE(POINTS)! Fraction of ice in layer above (no units)
    !
    ! Process rate diagnostics (Dummy arrays in 3C version of code)
    !
    ! Sets up the size of arrays for CX and CONSTP
    REAL(KIND=r8)   , INTENT(IN   ) ::  CX(80)
    REAL(KIND=r8)   , INTENT(IN   ) ::  CONSTP(80)

    !
    LOGICAL, INTENT(IN   ) :: L_autoconv_murk! Use murk aerosol to calc. drop number

    REAL(KIND=r8)   , INTENT(IN   ) :: ec_auto       ! Collision coalescence efficiency
    REAL(KIND=r8)   , INTENT(IN   ) :: Ntot_land     ! Number of droplets over land / m-3
    REAL(KIND=r8)   , INTENT(IN   ) :: Ntot_sea      ! Number of droplets over sea  / m-3

    ! Declarations:
    !
    ! Global variables:
    !*L------------------COMDECK C_O_DG_C-----------------------------------
    ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
    ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
    ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

    REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
    REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
    REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8

    !*----------------------------------------------------------------------
    ! --------------------------COMDECK C_LSPMIC----------------------------
    ! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
    ! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
    ! THE MICROPHYSICS AND ICE CLOUD FRACTION METHOD
    ! ----------------------------------------------------------------------
    !
    ! History:
    !
    ! Version    Date     Comment
    ! -------    ----     -------
    !   5.4    16/08/02   Correct comment line, add PC2 parameters and
    !                     move THOMO to c_micro     Damian Wilson
    !   6.0    11/08/03   Correct value of wind_shear_factor for PC2
    !                                                          Damian Wilson
    !   6.2    17/11/05   Remove variables that are now in UMUI. D. Wilson
    !   6.2    03/02/06   Include droplet settling logical. Damian Wilson
    !
    ! ----------------------------------------------------------------------
    !      AUTOCONVERSION TERMS
    ! ----------------------------------------------------------------------
    !
    !     LOGICAL, PARAMETER :: L_AUTOCONV_MURK is set in UMUI
    ! Set to .TRUE. to calculate droplet concentration from MURK aerosol,
    ! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
    ! of sulphate aerosol). If both are .FALSE., droplet concentrations
    ! from comdeck C_MICRO are used to be consistent with the values
    ! used in the radiation scheme.

    ! This next set of parameters is to allow the 3B scheme to
    ! be replicated at 3C/3D
    ! Inhomogeneity factor for autoconversion rate
    REAL(KIND=r8),PARAMETER:: INHOMOG_RATE=1.0_r8

    ! Inhomogeneity factor for autoconversion limit
    REAL(KIND=r8),PARAMETER:: INHOMOG_LIM=1.0_r8

    ! Threshold droplet radius for autoconversion
    REAL(KIND=r8),PARAMETER:: R_THRESH=7.0E-6_r8
    ! End of 3B repeated code

    !Do not alter R_AUTO and N_AUTO since these values are effectively
    ! hard wired into a numerical approximation in the autoconversion
    ! code. EC_AUTO will be multiplied by CONSTS_AUTO

    ! Threshold radius for autoconversion
    REAL(KIND=r8), PARAMETER :: R_AUTO=20.0E-6_r8

    ! Critical droplet number for autoconversion
    REAL(KIND=r8), PARAMETER :: N_AUTO=1000.0_r8

    ! Collision coalesence efficiency for autoconversion
    !      REAL(KIND=r8), PARAMETER :: EC_AUTO is set in UMUI

    ! The autoconversion powers define the variation of the rate with
    ! liquid water content and droplet concentration. The following are
    ! from Tripoli and Cotton

    !  Dependency of autoconversion rate on droplet concentration
    REAL(KIND=r8), PARAMETER :: POWER_DROPLET_AUTO=-0.33333_r8

    ! Dependency of autoconversion rate on water content
    REAL(KIND=r8), PARAMETER :: POWER_QCL_AUTO=2.33333_r8

    ! Dependency of autoconversion rate on air density
    REAL(KIND=r8), PARAMETER :: power_rho_auto=1.33333_r8

    ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
    ! See UM documentation paper 26, equation P26.132

    ! Combination of physical constants
    REAL(KIND=r8), PARAMETER :: CONSTS_AUTO=5907.24_r8

    ! Quantites for calculation of drop number by aerosols.
    ! Need only set if L_AUTOCONV_MURK=.TRUE.  See file C_VISBTY

    ! Scaling concentration (m-3) in droplet number concentration
    REAL(KIND=r8), PARAMETER :: N0_MURK=500.0E6_r8

    ! Scaling mass (kg/kg) in droplet number calculation from aerosols
    REAL(KIND=r8), PARAMETER :: M0_MURK=1.458E-8_r8

    ! Power in droplet number calculation from aerosols
    REAL(KIND=r8), PARAMETER :: POWER_MURK=0.5_r8

    ! Ice water content threshold for graupel autoconversion (kg/m^3)
    REAL(KIND=r8), PARAMETER :: AUTO_GRAUP_QCF_THRESH = 3.E-4_r8

    ! Temperature threshold for graupel autoconversion (degC)
    REAL(KIND=r8), PARAMETER :: AUTO_GRAUP_T_THRESH = -4.0_r8

    ! Temperature threshold for graupel autoconversion
    REAL(KIND=r8), PARAMETER :: AUTO_GRAUP_COEFF = 0.5_r8

    !-----------------------------------------------------------------
    ! Iterations of microphysics
    !-----------------------------------------------------------------

    ! Number of iterations in microphysics.
    INTEGER,PARAMETER :: LSITER=8
    ! Advise 1 iteration for every 10 minutes or less of timestep.

    !-----------------------------------------------------------------
    ! Nucleation of ice
    !-----------------------------------------------------------------

    ! Note that the assimilation scheme uses temperature thresholds
    ! in its calculation of qsat.

    ! Nucleation mass
    REAL(KIND=r8), PARAMETER :: M0=1.0E-12_r8

    ! Maximum Temp for ice nuclei nucleation (deg C)
    REAL(KIND=r8), PARAMETER :: TNUC=-10.0_r8

    ! Maximum temperature for homogenous nucleation is now in c_micro
    ! so that it is available to code outside of section A04.

    !  1.0/Scaling quantity for ice in crystals
    REAL(KIND=r8), PARAMETER :: QCF0=1.0E4_r8       ! This is an inverse quantity

    ! Minimum allowed QCF after microphysics
    REAL(KIND=r8),PARAMETER:: QCFMIN=1.0E-8_r8

    ! 1/scaling temperature in aggregate fraction calculation
    REAL(KIND=r8), PARAMETER :: T_SCALING=0.0384_r8

    !  Minimum temperature limit in calculation  of N0 for ice (deg C)
    REAL(KIND=r8), PARAMETER :: T_AGG_MIN=-45.0_r8

    !-----------------------------------------------------------------
    ! Hallett Mossop process
    !-----------------------------------------------------------------

    ! Switch off Hallett Mossop in this version but allow
    ! functionality

    ! Min temp for production of Hallett Mossop splinters (deg C)
    REAL(KIND=r8), PARAMETER :: HM_T_MIN=-8.0_r8

    ! Max temp for production of Hallett Mossop splinters (deg C)
    REAL(KIND=r8), PARAMETER :: HM_T_MAX=-273.0_r8
    ! REAL(KIND=r8), PARAMETER :: HM_T_MAX=-3.0_r8

    !  Residence distance for Hallett Mossop splinters (1/deg C)
    REAL(KIND=r8), PARAMETER :: HM_DECAY=1.0_r8/7.0_r8

    ! Reciprocal of scaling liquid water content for HM process
    REAL(KIND=r8), PARAMETER :: HM_RQCL=1.0_r8/0.1E-3_r8

    !-----------------------------------------------------------------
    ! PC2 Cloud Scheme Terms
    !-----------------------------------------------------------------

    ! Specifies the ice content (in terms of a fraction of qsat_liq)
    ! that corresponds to a factor of two reduction in the width of
    ! the vapour distribution in the liquid-free part of the gridbox.
    REAL(KIND=r8), PARAMETER :: ICE_WIDTH=0.04_r8

    ! Parameter that governs the rate of spread of ice cloud fraction
    ! due to windshear
    REAL(KIND=r8), PARAMETER :: WIND_SHEAR_FACTOR = 1.5E-4_r8

    !-----------------------------------------------------------------
    ! Droplet settling
    !-----------------------------------------------------------------
    ! Use the droplet settling code (not available for 3C).
    LOGICAL, PARAMETER :: l_droplet_settle = .FALSE.

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
    !     REAL(KIND=r8),PARAMETER:: NTOT_LAND is set in UMUI
    !     REAL(KIND=r8),PARAMETER:: NTOT_SEA is set in UMUI
    REAL(KIND=r8),PARAMETER:: KPARAM_LAND = 0.67_r8
    REAL(KIND=r8),PARAMETER:: KPARAM_SEA = 0.80_r8
    REAL(KIND=r8),PARAMETER:: DCONRE_LAND = 9.5E-06_r8
    REAL(KIND=r8),PARAMETER:: DCONRE_SEA = 16.0E-06_r8
    REAL(KIND=r8),PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0_r8
    REAL(KIND=r8),PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0_r8
    !
    ! Maximum Temp for homogenous nucleation (deg C)
    REAL(KIND=r8),PARAMETER:: THOMO = -40.0_r8
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
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8

    ! C_LHEAT end
    ! --------------------------COMDECK C_LSPDIF---------------------------
    ! input variables
    !
    ! Values of reference variables
    REAL(KIND=r8),PARAMETER:: AIR_DENSITY0=1.0_r8             ! kg m-3
    REAL(KIND=r8),PARAMETER:: AIR_VISCOSITY0=1.717E-5_r8      ! kg m-1 s-1
    REAL(KIND=r8),PARAMETER:: AIR_CONDUCTIVITY0=2.40E-2_r8    ! J m-1 s-1 K-1
    REAL(KIND=r8),PARAMETER:: AIR_DIFFUSIVITY0=2.21E-5_r8     ! m2 s-1
    REAL(KIND=r8),PARAMETER:: AIR_PRESSURE0=1.0E5_r8          ! Pa

    ! Values of diffusional growth parameters
    ! Terms in deposition and sublimation
    REAL(KIND=r8),PARAMETER:: APB1=(LC+LF)**2 * EPSILON /(R*AIR_CONDUCTIVITY0)
    REAL(KIND=r8),PARAMETER:: APB2=(LC+LF) / AIR_CONDUCTIVITY0
    REAL(KIND=r8),PARAMETER:: APB3=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)
    ! Terms in evap of melting snow and rain
    REAL(KIND=r8),PARAMETER:: APB4=LC**2*EPSILON/(R*AIR_CONDUCTIVITY0)
    REAL(KIND=r8),PARAMETER:: APB5=LC /AIR_CONDUCTIVITY0
    REAL(KIND=r8),PARAMETER:: APB6=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)

    ! Values of numerical approximation to wet bulb temperature
    ! Numerical fit to wet bulb temperature
    REAL(KIND=r8),PARAMETER:: TW1=1329.31_r8
    REAL(KIND=r8),PARAMETER:: TW2=0.0074615_r8
    REAL(KIND=r8),PARAMETER:: TW3=0.85E5_r8
    ! Numerical fit to wet bulb temperature
    REAL(KIND=r8),PARAMETER:: TW4=40.637_r8
    REAL(KIND=r8),PARAMETER:: TW5=275.0_r8

    ! Ventilation parameters
    REAL(KIND=r8),PARAMETER:: SC=0.6_r8
    ! f(v)  =  VENT_ICE1 + VENT_ICE2  Sc**(1/3) * Re**(1/2)
    REAL(KIND=r8),PARAMETER:: VENT_ICE1=0.65_r8
    REAL(KIND=r8),PARAMETER:: VENT_ICE2=0.44_r8
    ! f(v)  =  VENT_RAIN1 + VENT_RAIN2  Sc**(1/3) * Re**(1/2)
    REAL(KIND=r8),PARAMETER:: VENT_RAIN1=0.78_r8
    REAL(KIND=r8),PARAMETER:: VENT_RAIN2=0.31_r8
    ! c_lspdif will call c_r_cp, c_epslon and c_lheat
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
    !
    ! Subroutine arguments
    !
    ! Obtain the size for CX and CONSTP
    ! Start C_LSPSIZ
    ! Description: Include file containing idealised forcing options
    ! Author:      R. Forbes
    !
    ! History:
    ! Version  Date      Comment
    ! -------  ----      -------
    !   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
    !   6.2    22/08/05  Include the step size between ice categories.
    !                                                   Damian Wilson

    ! Sets up the size of arrays for CX and CONSTP
    !      REAL(KIND=r8) CX(80),CONSTP(80)
    INTEGER,PARAMETER:: ice_type_offset=20


    !  Local parameters
    REAL(KIND=r8),PARAMETER:: LCRCP=LC/CP!         Latent heat of condensation / cp (K).
    REAL(KIND=r8),PARAMETER:: LFRCP=LF/CP!         Latent heat of fusion / cp (K).
    REAL(KIND=r8),PARAMETER:: LSRCP=LCRCP+LFRCP!         Sum of the above (S for Sublimation).
    REAL(KIND=r8),PARAMETER:: RHO1=1.0_r8!         Reference density of air (kg m-3)
    REAL(KIND=r8),PARAMETER:: ONE_OVER_EPSILON=1.0_r8/EPSILON!        Inverse of epsilon to speed up calculations
    REAL(KIND=r8),PARAMETER:: ONE_OVER_ZERODEGC=1.0_r8/ZERODEGC!        Inverse of zero degrees Celsius to speed up code (K-1)
    !
    !  Local scalars and dynamic arrays
    !
    INTEGER :: I            !         Loop counter (horizontal field index).
    INTEGER :: J            !         Counter for the iterations
    INTEGER :: K            !         Counter for the LSITER2 loop
    INTEGER :: LSITER2      !         Number of times advection and melting sections are iterated
    INTEGER :: KK           !         Variable for condensed points compression
    INTEGER :: SEA_SALT_PTR !         Pointer for sea-salt arrays
    INTEGER :: BIOMASS_PTR  !         Pointer for biomass smoke arrays
    INTEGER :: OCFF_PTR     !         Pointer for fossil-fuel organic carbon arrays
    INTEGER :: biogenic_ptr !         Pointer for biogenic aerosol arrays
    !
    REAL(KIND=r8) :: QS(POINTS)         !         Saturated sp humidity for (T,p) in layer (kg kg-1)
    REAL(KIND=r8) :: QSL(POINTS)        !         Saturated sp humidity for (T,p) in layer
                                        !         wrt water at all temps (kg kg-1)
    REAL(KIND=r8) :: SNOWT_AGG(POINTS)  !         Cumulative fall out of snow aggregates within iterations.
                                        !         (kg m-2 s-1). If only one ice prognostic is active then
                                        !         this variable contains all the snow.
    REAL(KIND=r8) :: SNOWT_CRY(POINTS)  !         Cumulative fall out of snow crystals within iterations.
                                        !         (kg m-2 s-1). Only non-zero if two ice prognostics in use.
                                        !REAL(KIND=r8) :: NUMBER_DROPLET     !         Droplet concentration calculated using a function (m-3)
    REAL(KIND=r8) :: RHO(POINTS)        !         Density of air in the layer (kg m-3).
    REAL(KIND=r8) :: RHOR(POINTS)       !         1.0/RHO to speed up calculations (kg-1 m3).
    REAL(KIND=r8) :: ESI(POINTS)        !         saturation vapour pressure (wrt ice below zero Celsius)(Pa)
    REAL(KIND=r8) :: ESW(POINTS)        !         saturation vapour pressure (wrt water at all temperatures)(Pa)
    REAL(KIND=r8) :: DQI                !         increment to/from ice/snow (kg kg-1)
    REAL(KIND=r8) :: DQIL               !         increment to/from cloud water (kg kg-1)
    REAL(KIND=r8) :: DPR                !         increment to/from rain (kg m-2 s-1)
    REAL(KIND=r8) :: CFICE(POINTS)      !         fraction of ice inferred for the microphysics (no units).
    REAL(KIND=r8) :: CFICEI(POINTS)     !         inverse of CFICE (no units)
    REAL(KIND=r8) :: CF(POINTS)         !         total cloud fraction for the microphysics (no units)
    REAL(KIND=r8) :: FQI_AGG(POINTS)    !         fallspeed for aggregates (m s-1)
    REAL(KIND=r8) :: FQI_CRY(POINTS)    !         fallspeed for aggregates (m s-1)
    REAL(KIND=r8) :: DHI(POINTS)        !         CFL limit (s m-1)
    REAL(KIND=r8) :: DHIR(POINTS)       !         1.0/DHI (m s-1)
    REAL(KIND=r8) :: DHILSITERR(POINTS) !         1.0/(DHI*LSITER) (m s-1)
    REAL(KIND=r8) :: FQIRQI_AGG         !         saved flux of ice out of layer (kg m-2 s-1)
    REAL(KIND=r8) :: FQIRQI2_AGG(POINTS)!         saved flux of ice out of layer from layer above (kg m-2 s-1)
    REAL(KIND=r8) :: FQIRQI_CRY         !         saved flux of ice out of layer (kg m-2 s-1)
    REAL(KIND=r8) :: FQIRQI2_CRY(POINTS)!         saved flux of ice out of layer from layer above (kg m-2 s-1)
    REAL(KIND=r8) :: QCLNEW             !         updated liquid cloud in implicit calculations (kg kg-1)
    REAL(KIND=r8) :: TEMP7              !         temporary variable
    REAL(KIND=r8) :: TEMPW              !         temporary for vapour calculations
    REAL(KIND=r8) :: t_rapid_melt       !         Rapid melting temperature
    REAL(KIND=r8) :: PR02               !         term in evaporation of rain
    REAL(KIND=r8) :: PR04               !         square of pr02
    REAL(KIND=r8) :: QC                 !         term in autoconversion of cloud to rain (kg kg-1)
    REAL(KIND=r8) :: APLUSB             !         denominator in deposition or evaporation of ice
    REAL(KIND=r8) :: CORR(POINTS)       !         density correction for fall speed (no units)
    REAL(KIND=r8) :: ROCOR(POINTS)      !         density correction for fall speed (no units)
    REAL(KIND=r8) :: VR1                !         Mean fall speed of rain (m s-1)
    REAL(KIND=r8) :: VS1                !         Mean fall speed of snow (m s-1)
    REAL(KIND=r8) :: LAMR1              !         Inverse lambda in rain exponential distribution (m)
    REAL(KIND=r8) :: LAMR2              !         Inverse lambda in rain exponential distribution (m)
    REAL(KIND=r8) :: LAMFAC1            !         Expression containing calculations with lambda
    REAL(KIND=r8) :: LAMS1              !         Inverse lambda in snow exponential distribution (m)
    REAL(KIND=r8) :: FV1                !         Mean velocity difference between rain and snow (m s-1)
    REAL(KIND=r8) :: TIMESTEP           !         Timestep of each iteration (s)
    REAL(KIND=r8) :: CORR2(POINTS)      !         Temperature correction of viscosity etc. (no units)
    REAL(KIND=r8) :: RHNUC              !         Relative humidity required for nucleation (no units)
    REAL(KIND=r8) :: TCG(POINTS)        !         Temperature Factor for aggregate size distribution (no units)
    REAL(KIND=r8) :: TCGI(POINTS)       !         Inverse of TCG (no units)
    REAL(KIND=r8) :: TCGC(POINTS)       !         Temperature Factor for crystal size distribution (no units)
    REAL(KIND=r8) :: TCGCI(POINTS)      !         Inverse of TCGC (no units)
    REAL(KIND=r8) :: RATEQS(POINTS)     !         Sub grid model variable (no units)
    REAL(KIND=r8) :: HM_NORMALIZE       !         Normalization for Hallett Mossop process (no units)
    REAL(KIND=r8) :: HM_RATE            !         Increase in deposition due to Hallett Mossop process(no units)
    REAL(KIND=r8) :: AREA_LIQ(POINTS)   !         Liquid only area of gridbox (no units)
    REAL(KIND=r8) :: AREA_MIX(POINTS)   !         Mixed phase area of gridbox (no units)
    REAL(KIND=r8) :: AREA_ICE(POINTS)   !         Ice only area of gridbox (no units)
    REAL(KIND=r8) :: AREA_CLEAR(POINTS) !         Cloud free area of gridbox (no units)
    REAL(KIND=r8) :: RAIN_LIQ(POINTS)   !         Overlap fraction of gridbox between rain and liquid cloud
    REAL(KIND=r8) :: RAIN_MIX(POINTS)   !         Overlap fraction of gridbox between rain and mixed phase cloud
    REAL(KIND=r8) :: RAIN_ICE(POINTS)   !         Overlap fraction of gridbox between rain and ice cloud
    REAL(KIND=r8) :: RAIN_CLEAR(POINTS) !         Overlap fraction of gridbox between rain and no cloud
    REAL(KIND=r8) :: Q_ICE(POINTS)      !         Vapour content in the ice only part of the grid box (kg kg-1)
    REAL(KIND=r8) :: Q_CLEAR(POINTS)    !         Vapour content in the cloud free part of the grid box(kg kg-1)
    REAL(KIND=r8) :: QCF_AGG(POINTS)    !         QCF in the form of aggregates (kg kg-1)
    REAL(KIND=r8) :: QCF_CRY(POINTS)    !         QCF in the form of crystals (kg kg-1)
    REAL(KIND=r8) :: FRAC_AGG(POINTS)   !         Fraction of aggregates (no units)
    REAL(KIND=r8) :: TEMP1(POINTS)      !         Temporary arrays for T3E vector functions
    REAL(KIND=r8) :: TEMP3(POINTS)      !         Temporary arrays for T3E vector functions
    REAL(KIND=r8) :: N_DROP(POINTS)     !         Droplet concentration (m-3)
    REAL(KIND=r8) :: A_FACTOR(POINTS)   !         Numerical factors in autoconversion calculation
    REAL(KIND=r8) :: B_FACTOR(POINTS)   !         Numerical factors in autoconversion calculation
    REAL(KIND=r8) :: R_MEAN(POINTS)     !         Mean droplet radius (m)
    REAL(KIND=r8) :: R_MEAN0            !         Constant in the calculation of R_MEAN
    REAL(KIND=r8) :: N_GT_20            !         Concentration of particles greater than a threshold radius in size (m-3)
    REAL(KIND=r8) :: AUTOLIM            !         Minimum water content for autoconversion (kg m-3)
    REAL(KIND=r8) :: AUTORATE           !         Rate constant for autoconversion
    REAL(KIND=r8) :: AC_FACTOR          !         Autoconversion de-biasing factor
    REAL(KIND=r8) :: T_L(POINTS)        !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: QSL_TL(POINTS)     !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: ALPHA_L            !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: A_L                !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: SIGMA_S            !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: G_L                !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: GACB               !         Subsidiary terms used in the autoconversion de-biasing
    REAL(KIND=r8) :: QCFAUTOLIM         !         Minimum ice content for autoconversion (kg/kg)
    REAL(KIND=r8) :: QCFAUTORATE        !         Rate constant for ice autoconversion
    REAL(KIND=r8) :: QCF_TOT(POINTS)    !         Total amount of ice (crystals+aggregates)
    REAL(KIND=r8) :: LHEAT_CORREC_LIQ(POINTS)!         Reduction factor in evaporation limits because of latent heat
    REAL(KIND=r8) :: LHEAT_CORREC_ICE(POINTS)!         Reduction factor in evaporation limits because of latent heat
    REAL(KIND=r8) :: TEMPR              !         Temporary for optimization
    REAL(KIND=r8) :: DQI_DEP
    REAL(KIND=r8) :: DQI_SUB
    REAL(KIND=r8) :: TEMPW_DEP
    REAL(KIND=r8) :: TEMPW_SUB
    REAL(KIND=r8) :: WIDTH(POINTS)
    REAL(KIND=r8) :: Q_ICE_1(POINTS)
    REAL(KIND=r8) :: Q_ICE_2(POINTS)
    REAL(KIND=r8) :: AREA_ICE_1(POINTS)
    REAL(KIND=r8) :: AREA_ICE_2(POINTS)!         Subgrid splitting for deposition term
    LSITER2     =0;    KK =0;    SEA_SALT_PTR=0;    BIOMASS_PTR =0
    OCFF_PTR    =0;    biogenic_ptr=0
    QS=0.0_r8;QSL=0.0_r8;SNOWT_AGG=0.0_r8;SNOWT_CRY=0.0_r8;RHO=0.0_r8;RHOR=0.0_r8;
    ESI=0.0_r8;ESW=0.0_r8;DQI=0.0_r8;DQIL=0.0_r8;DPR=0.0_r8;CFICE=0.0_r8
    CFICEI=0.0_r8;CF=0.0_r8;FQI_AGG=0.0_r8;FQI_CRY=0.0_r8;DHI=0.0_r8;DHIR=0.0_r8;
    DHILSITERR=0.0_r8; FQIRQI_AGG=0.0_r8;FQIRQI2_AGG=0.0_r8;FQIRQI_CRY=0.0_r8;
    FQIRQI2_CRY=0.0_r8;QCLNEW=0.0_r8;TEMP7=0.0_r8;TEMPW=0.0_r8;t_rapid_melt=0.0_r8;
    PR02=0.0_r8;PR04=0.0_r8;QC=0.0_r8;APLUSB=0.0_r8;CORR=0.0_r8;ROCOR=0.0_r8;
    VR1=0.0_r8;VS1=0.0_r8;LAMR1=0.0_r8;LAMR2=0.0_r8;LAMFAC1=0.0_r8;LAMS1=0.0_r8;
    FV1=0.0_r8;TIMESTEP=0.0_r8;CORR2=0.0_r8;RHNUC=0.0_r8;TCG=0.0_r8;
    TCGI=0.0_r8;TCGC=0.0_r8;TCGCI=0.0_r8;RATEQS=0.0_r8;HM_NORMALIZE=0.0_r8;
    HM_RATE=0.0_r8;AREA_LIQ=0.0_r8;AREA_MIX=0.0_r8;AREA_ICE=0.0_r8;AREA_CLEAR=0.0_r8; 
    RAIN_LIQ=0.0_r8;RAIN_MIX=0.0_r8;RAIN_ICE=0.0_r8;RAIN_CLEAR=0.0_r8;Q_ICE=0.0_r8;
    Q_CLEAR=0.0_r8;QCF_AGG=0.0_r8;QCF_CRY=0.0_r8;FRAC_AGG=0.0_r8;TEMP1=0.0_r8;TEMP3=0.0_r8;
    N_DROP=0.0_r8;A_FACTOR=0.0_r8;B_FACTOR=0.0_r8;R_MEAN=0.0_r8;R_MEAN0=0.0_r8;
    N_GT_20=0.0_r8;AUTOLIM=0.0_r8;AUTORATE=0.0_r8;AC_FACTOR=0.0_r8;T_L=0.0_r8;
    QSL_TL=0.0_r8;ALPHA_L=0.0_r8;A_L=0.0_r8;SIGMA_S=0.0_r8;G_L=0.0_r8;
    GACB=0.0_r8;QCFAUTOLIM=0.0_r8;QCFAUTORATE=0.0_r8;QCF_TOT=0.0_r8;
    LHEAT_CORREC_LIQ=0.0_r8;LHEAT_CORREC_ICE=0.0_r8;TEMPR=0.0_r8;DQI_DEP=0.0_r8;
    DQI_SUB=0.0_r8;TEMPW_DEP=0.0_r8;TEMPW_SUB=0.0_r8;WIDTH=0.0_r8;Q_ICE_1=0.0_r8;
    Q_ICE_2=0.0_r8;AREA_ICE_1=0.0_r8;AREA_ICE_2=0.0_r8;
    !
    !  Function and Subroutine calls
    !EXTERNAL NUMBER_DROPLET
    !
    !- End of header
    !
    ! ----------------------------------------------------------------------
    !  2.1 Set up some variables
    ! ----------------------------------------------------------------------
    ! Set up the iterations
    TIMESTEP=TIMESTEPFIXED/LSITER

    DO I=1,POINTS
       ! Set up SNOWT to be zero for all I
       SNOWT_CRY(I) = 0.0_r8
       SNOWT_AGG(I) = 0.0_r8
    END DO

    ! Set qcf for crystals and aggregates depending on
    ! whether there are one or two ice prognostics

    IF (L_mcr_qcf2) THEN ! two ice prognostics

       DO I=1,POINTS
          QCF_CRY(I) = QCF2(I)
          QCF_AGG(I) = QCF(I)
       END DO

    ELSE ! only one ice prognostic, split diagnostically

       DO I=1,POINTS ! Points 0
          ! Work out fraction of ice in aggregates
          FRAC_AGG(I)=MAX(1.0_r8-EXP(-T_SCALING*MAX((T(I)-CTTEMP(I)),0.0_r8)   &
               &               *MAX(QCF(I),0.0_r8)*QCF0) , 0.0_r8)
       END DO ! Points 0
       ! Points_do1:
       DO I=1,POINTS
          ! Allocate ice content to crystals and aggregates
          QCF_CRY(I)=QCF(I)*(1.0_r8-FRAC_AGG(I))
          QCF_AGG(I)=QCF(I)*FRAC_AGG(I)
          ! Assume falling snow is partitioned into crystals and
          ! aggregates, SNOW_AGG contains total snow on input
          SNOW_CRY(I)=SNOW_AGG(I)*(1.0_r8-FRAC_AGG(I))
          SNOW_AGG(I)=SNOW_AGG(I)*FRAC_AGG(I)

          ! The compiler should unroll the above loop so break it here
       END DO
    ENDIF ! on L_mcr_qcf2 (number of ice prognostics)

    DO I=1,POINTS
       QCF_TOT(I) = QCF_CRY(I) + QCF_AGG(I)
    END DO

    ! Set up Hallett Mossop calculation
    HM_NORMALIZE=1.0_r8/(1.0_r8-EXP((HM_T_MIN-HM_T_MAX)*HM_DECAY))
    ! Set up autoconversion factor
    R_MEAN0=(27.0_r8/(80.0_r8*PI*1000.0_r8))**(-1.0_r8/3.0_r8)
    ! ----------------------------------------------------------------------
    !  2.2 Start iterating.
    ! ----------------------------------------------------------------------
    ! Iters_do1:
    DO J=1,LSITER
       ! ----------------------------------------------------------------------
       !  2.3  Calculate saturation specific humidities
       ! ----------------------------------------------------------------------
       ! Qsat with respect to ice
       CALL QSAT_mix ( &
                                !      Output field
            &  qs &             ! REAL(KIND=r8), INTENT(out) :: QmixS(npnts)  ! Output Saturation mixing ratio or saturation
                                ! specific humidity at temperature T and pressure
                                ! P (kg/kg).
                                !      Input fields
            &, T &              ! REAL(KIND=r8), INTENT(in)  :: T(npnts)      !  Temperature (K).
            &, p         &      ! REAL(KIND=r8), INTENT(in)  :: P(npnts)      !  Pressure (Pa).
                                !      Array dimensions
            &, points        &  !INTEGER, INTENT(in) :: npnts    ! Points (=horizontal dimensions) being processed by qSAT scheme
                                !      logical control
            &, l_mixing_ratio & !LOGICAL, INTENT(in)  :: lq_mix      !  .true. return qsat as a mixing ratio
                                !  .false. return qsat as a specific humidity
            &  )

       !
       ! Qsat with respect to liquid water
       ! DEPENDS ON: qsat_wat_mix
       CALL qsat_wat_mix ( &
            !      Output field
            &  qsl              &! Real, intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
            !      Input fields
            &, T                &! Real, intent(in)  :: T(npnts) !  Temperature (K).
            &, P                &! Real, intent(in)  :: P(npnts) !  Pressure (Pa).  
            !      Array dimensions
            &, points           &! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
            !      logical control
            &, l_mixing_ratio   &! logical, intent(in)  :: lq_mix
                                ! .true. return qsat as a mixing ratio
                                ! .false. return qsat as a specific humidity
            &  )

       ! ----------------------------------------------------------------------
       !  2.4 Start loop over points.
       ! ----------------------------------------------------------------------
       ! Points_do2:
       DO I=1,POINTS
          !
          ! ----------------------------------------------------------------------
          !  3.1 Calculate density of air and density dependent correction factors
          ! ----------------------------------------------------------------------
          ESI(I)=QS(I)*P(I)*ONE_OVER_EPSILON
          ESW(I)=QSL(I)*P(I)*ONE_OVER_EPSILON

          !-----------------------------------------------
          ! Calculate density of air
          !-----------------------------------------------
          IF (l_non_hydrostatic) THEN
             IF (l_mixing_ratio) THEN
                ! rho is the dry density
                rho(i) = rhodz_dry(i) / deltaz(i)
             ELSE
                ! rho is the moist density
                rho(i) = rhodz_moist(i) / deltaz(i)
             END IF  ! l_mixing_ratio
          ELSE
             ! Use the pressure to retrieve the moist density.
             ! An exact expression for the air density is
             ! rho = p / (1+c_virtual q - qcl - qcf) / RT but
             ! the approximation below is more numerically
             ! stable.
             rho(i) = p(i)*(1.0_r8-c_virtual*q(i)+qcl(i)+qcf_tot(i))           &
                  &           /(R*T(i))
          END IF  ! l_non_hydrostatic

          RHOR(I)=1.0_r8/MAX(RHO(I),0.00001_r8)
          ! Estimate latent heat correction to rate of evaporation etc.
          LHEAT_CORREC_LIQ(I)=1.0_r8/(1.0_r8+EPSILON*LC**2*QSL(I)              &
               &                            /(CP*R*T(I)**2))
          LHEAT_CORREC_ICE(I)=1.0_r8/(1.0_r8+EPSILON*(LC+LF)**2*QS(I)          &
               &                            /(CP*R*T(I)**2))
       END DO
       DO I=1,POINTS
          ! Correction factor of fall speeds etc. due to density.
          IF (l_mixing_ratio .AND. l_non_hydrostatic) THEN
             corr(i) = (rho1*deltaz(i) / rhodz_moist(i))**0.4_r8
          ELSE
             corr(i) = (rho1*rhor(i))**0.4_r8
          END IF
          ! Correction factor in viscosity etc. due to temperature.
          CORR2(I)=(T(I)*ONE_OVER_ZERODEGC)**1.5_r8 * (393.0_r8/(T(I)+120.0_r8))
       ENDDO
       DO I=1,POINTS
          ! ----------------------------------------------------------------------
          !  3.2 Calculate particle size distributions. Use vector functions
          !       if available because they are much faster
          ! ----------------------------------------------------------------------
          ! Combined correction factor
          IF (l_mixing_ratio .AND. l_non_hydrostatic) THEN
             rocor(i) = SQRT(rhodz_moist(i)/deltaz(i)*corr(i)*corr2(i))
          ELSE
             rocor(i) = SQRT(rho(i)*corr(i)*corr2(i))
          END IF
          ! Calculate a temperature factor for N0aggregates. CX(32)>0.0
          ! if there is a temperature dependence, and 0.0 if there is not.
          TCG(I)=EXP(-CX(32)*MAX(T(I)-ZERODEGC,T_AGG_MIN))
          ! Temperature dependence for N0crystals
          TCGC(I)=EXP(-CX(12)*MAX(T(I)-ZERODEGC,T_AGG_MIN))
          ! Define inverse of TCG values to speed up calculations
          TCGI(I)=1.0_r8/TCG(I)
          TCGCI(I)=1.0_r8/TCGC(I)
       ENDDO
       DO I=1,POINTS
          !
          ! ----------------------------------------------------------------------
          !  4.1 Check that ice cloud fraction is sensible.
          ! ----------------------------------------------------------------------
          CFICE(I)=MAX(MAX(CFICEKEEP(I),0.01_r8),FRAC_ICE_ABOVE(I))
          CFICEI(I)=1.0_r8/CFICE(I)
          CF(I)=CFKEEP(I)
          CF(I)=MAX(CF(I),CFICE(I))
          ! Break loop to aid efficient pipelining of calculations
       END DO
       DO I=1,POINTS
          RATEQS(I)=RHCPT(I)   ! Sub grid parameter
          FRAC_ICE_ABOVE(I)=CFICEKEEP(I)
          !
          ! ----------------------------------------------------------------------
          !  4.2 Calculate overlaps of liquid, ice and rain fractions
          ! ----------------------------------------------------------------------
          AREA_LIQ(I)=CF(I)-CFICE(I)
          AREA_MIX(I)=CFICE(I)+CFLIQ(I)-CF(I)
          AREA_ICE(I)=CF(I)-CFLIQ(I)
          AREA_CLEAR(I)=1.0_r8-CF(I)
          ! Break loop again since the above can be pipelined and the next loop
          ! the compiler should unroll
       END DO
       DO I=1,POINTS
          RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
          RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
          RAIN_ICE(I)=                                                   &
               &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
          RAIN_CLEAR(I)=RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
       END DO
       DO I=1,POINTS
          ! ----------------------------------------------------------------------
          !  4.3 Calculate vapour contents in ice only and clear regions
          ! ----------------------------------------------------------------------
          IF (CFLIQ(I)  <   1.0_r8) THEN
             ! TEMPW is the mean vapour content in the ice only and clear regions
             TEMPW=(Q(I) - CFLIQ(I)*QSL(I)) / (1.0_r8 - CFLIQ(I))
             TEMP7 = ICE_WIDTH * QSL(I)
             ! 0.001 is to avoid divide by zero problems
             WIDTH(I) = 2.0_r8 *(1.0_r8-RATEQS(I))*QSL(I)                       &
                  &                    *MAX(  (1.0_r8-0.5_r8*QCF(I)/TEMP7),0.001_r8  )
             ! The full width cannot be greater than 2q because otherwise part of
             ! gridbox would have negative q.
             ! However, we must also ensure that a q of zero cannot produce an error.
             WIDTH(I) = MIN(WIDTH(I) ,    MAX(2.0_r8*Q(I),0.001_r8*QS(I) )    )
             IF (AREA_ICE(I)  >   0.0_r8) THEN
                Q_CLEAR(I)= TEMPW - 0.5_r8*WIDTH(I)*AREA_ICE(I)
                Q_ICE(I)=(Q(I)-CFLIQ(I)*QSL(I)-AREA_CLEAR(I)*Q_CLEAR(I))  &
                     &                / AREA_ICE(I)
             ELSE
                Q_CLEAR(I) = TEMPW
                Q_ICE(I)=0.0_r8
             END IF
          ELSE
             ! Q_CLEAR and Q_ICE are undefined. Set them to zero
             WIDTH(I)=1.0_r8
             Q_CLEAR(I)=0.0_r8
             Q_ICE(I)=0.0_r8
          ENDIF
          ! ----------------------------------------------------------------------
          !  4.4 Update ice cloud top temperature if no ice falling in
          ! ----------------------------------------------------------------------
          IF (SNOW_CRY(I)+SNOW_AGG(I)  <=  0.0_r8) THEN
             CTTEMP(I)=T(I)
          ENDIF
          ! ----------------------------------------------------------------------
          !  4.5  Remove any small amount of ice to be tidy.
          !      If QCF is less than QCFMIN and isn't growing
          !      by deposition (assumed to be given by RHCPT) then remove it.
          ! ----------------------------------------------------------------------
          IF (QCF_TOT(I)  <   QCFMIN) THEN
             IF (T(I) >  ZERODEGC .OR.                                    &
                  &       (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0_r8)           &
                  &       .OR. QCF_TOT(I)  <   0.0_r8)  THEN
                Q(I)=Q(I)+QCF_TOT(I)
                T(I)=T(I)-LSRCP*QCF_TOT(I)
                QCF(I)     = 0.0_r8
                QCF_TOT(I) = 0.0_r8
                QCF_AGG(I) = 0.0_r8
                QCF_CRY(I) = 0.0_r8
             END IF
          END IF
       END DO
       ! ----------------------------------------------------------------------
       !  5.1   Falling ice is advected downwards
       ! ----------------------------------------------------------------------
       ! Estimate fall speed out of this layer. We want to avoid advecting
       ! very small amounts of snow between layers, as this can cause numerical
       ! problems in other routines, so if QCF is smaller than a single
       ! nucleation mass per metre cubed don't advect it.
       DO I=1,POINTS
          IF (QCF_AGG(I) >  M0) THEN
             FQI_AGG(I)=CONSTP(24)*CORR(I)*                               &
                  &     (RHO(I)*QCF_AGG(I)*CONSTP(25)*TCGI(I)*CFICEI(I))**CX(23)
          ELSE
             ! QCF is smaller than zero so set fall speed to zero
             FQI_AGG(I)=0.0_r8
             ! Endif for calculation of fall speed
          END IF
          IF (QCF_CRY(I) >  M0) THEN
             FQI_CRY(I)=CONSTP(4)*CORR(I)*                                &
                  &       (RHO(I)*QCF_CRY(I)*CONSTP(5)*TCGCI(I)*CFICEI(I))**CX(3)
          ELSE
             ! QCF is smaller than zero so set fall speed to zero
             FQI_CRY(I)=0.0_r8
             ! Endif for calculation of fall speed
          END IF
          !
          ! ----------------------------------------------------------------------
          !  5.2 Calculate CFL quantity of timestep over level separation.
          ! ----------------------------------------------------------------------
          IF (l_non_hydrostatic) THEN
             ! Use the formulation based on the heights of the levels
             dhi(i) = timestep/deltaz(i)
          ELSE
             ! Use the formulation based on the pressure difference
             ! across a layer.
             dhi(i)        = timestep * rho(i)/rhodz(i)
          END IF

          ! Define DHIR and DHILSITERR(I) to speed up calculations.
          DHIR(I)=1.0_r8/DHI(I)
          DHILSITERR(I)=1.0_r8/(DHI(I)*LSITER)
          !
          ! ----------------------------------------------------------------------
          !  5.2a Adjust fall speeds to make a linear combination of the fall
          !       speed from the layer above. This ensures that the fall speed of
          !       ice in this layer will not be calculated as zero even though
          !       there is ice falling into it from above.
          ! ----------------------------------------------------------------------
          ! If using two cloud ice prognostics, calculate fallspeeds
          !  separately for crystals and aggregates, in which case
          !  VF represents aggregates, VF_CRY represents crystals
          ! Otherwise calculate an average fallspeed represented by VF

          IF (L_mcr_qcf2) THEN

             ! Crystals
             TEMP7 = SNOW_CRY(I)*DHI(I)*RHOR(I)
             IF (QCF_CRY(I)+TEMP7  >   M0) THEN
                TEMPR = 1.0_r8/(QCF_CRY(I)+TEMP7)
                FQI_CRY(I)=(FQI_CRY(I)*QCF_CRY(I)+VF_CRY(I)*TEMP7)*TEMPR
                VF_CRY(I) = FQI_CRY(I)
             ELSE
                VF_CRY(I)=0.0_r8
             ENDIF

             ! Aggregates
             TEMP7 = SNOW_AGG(I)*DHI(I)*RHOR(I)
             IF (QCF_AGG(I)+TEMP7  >   M0) THEN
                TEMPR = 1.0_r8/(QCF_AGG(I)+TEMP7)
                FQI_AGG(I)=(FQI_AGG(I)*QCF_AGG(I)+VF(I)*TEMP7)*TEMPR
                VF(I) = FQI_AGG(I)
             ELSE
                VF(I)=0.0_r8
             ENDIF

          ELSE ! L_mcr_qcf2=.F.-> use totals to estimate fallspeed

             TEMP7 = (SNOW_AGG(I) + SNOW_CRY(I)) *DHI(I)*RHOR(I)
             IF((QCF_CRY(I)+QCF_AGG(I)+TEMP7)  >   M0) THEN
                TEMPR=1.0_r8 / (QCF_CRY(I)+QCF_AGG(I)+TEMP7)
                FQI_CRY(I)=(FQI_CRY(I)*(QCF_CRY(I)+QCF_AGG(I))+VF(I)*TEMP7)  &
                     &                * TEMPR
                FQI_AGG(I)=(FQI_AGG(I)*(QCF_CRY(I)+QCF_AGG(I))+VF(I)*TEMP7)  &
                     &                * TEMPR
                VF(I)=FQI_CRY(I)*(1.0_r8-FRAC_AGG(I))+FQI_AGG(I)*FRAC_AGG(I)
             ELSE
                VF(I)=0.0_r8
             ENDIF

          ENDIF ! on L_mcr_qcf2
          !
          ! ----------------------------------------------------------------------
          !  5.2b Calculate the additional iterations required by the advection
          ! and melting schemes.
          ! ----------------------------------------------------------------------
          !
          IF ((QCF_CRY(I)+QCF_AGG(I)) >  M0.AND.T(I) >  ZERODEGC.AND.    &
               &       T(I) <  (ZERODEGC+2.0_r8) .AND. l_it_melting) THEN
             LSITER2 = INT(MAX(FQI_CRY(I),FQI_AGG(I))*DHI(I) + 1)
             IF (LSITER2  >   5) THEN
                LSITER2 = 5
             ENDIF
          ELSE
             LSITER2 = 1
          ENDIF
          !
          ! Now start the additional iterations loop
          DO K=1,LSITER2
             !
             ! ----------------------------------------------------------------------
             !  5.3 Analytical Solution
             ! ----------------------------------------------------------------------
             ! FQIRQI2 is used as a temporary for the exponential function
             ! Calculate values point by point
             FQIRQI2_AGG(I)=EXP(-FQI_AGG(I)*DHI(I)/LSITER2)
             FQIRQI2_CRY(I)=EXP(-FQI_CRY(I)*DHI(I)/LSITER2)
             ! Advect aggregates
             IF (FQI_AGG(I)  >   0.0_r8) THEN
                ! Assume falling snow is partitioned into crystals and
                ! aggregates
                FQIRQI_AGG = SNOW_AGG(I)+DHIR(I)*LSITER2*                    &
                     &                  (RHO(I)*QCF_AGG(I)-SNOW_AGG(I)/FQI_AGG(I))      &
                     &                  * (1.0_r8-FQIRQI2_AGG(I))
                QCF_AGG(I) = SNOW_AGG(I)*RHOR(I)/FQI_AGG(I)                  &
                     &                  * (1.0_r8-FQIRQI2_AGG(I))                          &
                     &                  + QCF_AGG(I)*FQIRQI2_AGG(I)
             ELSE  ! No fall of QCF out of the layer
                FQIRQI_AGG = 0.0_r8
                QCF_AGG(I) = SNOW_AGG(I)*RHOR(I)*DHI(I)/LSITER2
             END IF

             ! Advect crystals
             IF (FQI_CRY(I)  >   0.0_r8) THEN
                FQIRQI_CRY = SNOW_CRY(I)+DHIR(I)*LSITER2*                    &
                     &                  (RHO(I)*QCF_CRY(I)-SNOW_CRY(I)/FQI_CRY(I))      &
                     &                  * (1.0_r8-FQIRQI2_CRY(I))
                QCF_CRY(I) = SNOW_CRY(I)*RHOR(I)/FQI_CRY(I)                  &
                     &                  * (1.0_r8-FQIRQI2_CRY(I))                          &
                     &                  + QCF_CRY(I)*FQIRQI2_CRY(I)
             ELSE  ! No fall of QCF out of the layer
                FQIRQI_CRY=0.0_r8
                QCF_CRY(I)=SNOW_CRY(I)*RHOR(I)*DHI(I)/LSITER2
             END IF
             ! No need to compute fall speed out of the layer in this method.
             ! ----------------------------------------------------------------------
             !  5.4 Snow is used to save fall out of layer
             !      for calculation of fall into next layer
             ! ----------------------------------------------------------------------
             SNOWT_CRY(I) = SNOWT_CRY(I) + (FQIRQI_CRY)/(LSITER*LSITER2)
             SNOWT_AGG(I) = SNOWT_AGG(I) + (FQIRQI_AGG)/(LSITER*LSITER2)
             !
             ! ----------------------------------------------------------------------
             !  11  Melting of snow. Uses a numerical approximation to the wet bulb
             !      temperature.
             ! ----------------------------------------------------------------------
             IF (L_it_melting) THEN
                ! Do iterative melting

                IF(QCF_CRY(I) >  M0.AND.T(I) >  ZERODEGC)THEN
                   ! ----------------------------------------------------------------------
                   !  11.1 Crystals first
                   ! ----------------------------------------------------------------------
                   ! An approximate calculation of wet bulb temperature
                   ! TEMPW represents AVERAGE supersaturation in the ice
                   ! Strictly speaking, we need to do two melting calculations,
                   ! for the two different wet bulb temperatures.
                   TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0_r8)*CFICEI(I)
                   TEMP7=T(I)-ZERODEGC-TEMPW                                    &
                        &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
                   TEMP7=MAX(TEMP7,0.0_r8)
                   ! End of wet bulb temp formulations.
                   PR02=RHO(I)*MAX(QCF_CRY(I),0.0_r8)*CFICEI(I)*CONSTP(5)*TCGCI(I)
                   DPR=TCGC(I)*CONSTP(14)*TIMESTEP/LSITER2*                     &
                        &            (CONSTP(7)*CORR2(I)*PR02**CX(4)                       &
                        &         + CONSTP(8)*ROCOR(I)*PR02**CX(5))*RHOR(I)
                   ! Solve implicitly in terms of temperature
                   DPR=TEMP7*(1.0_r8-1.0_r8/(1.0_r8+DPR*LFRCP))/LFRCP
                   DPR=MIN(DPR,QCF_CRY(I))
                   ! Update values of ice and Rain
                   QCF_CRY(I)=QCF_CRY(I)-DPR
                   RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
                   T(I)=T(I)-LFRCP*DPR
                   IF (DPR >  0.0_r8) THEN
                      RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
                      ! Update rain fractions
                      RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                      RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                      RAIN_ICE(I)=                                               &
                           &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                      RAIN_CLEAR(I)=                                             &
                           &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
                   ENDIF
                   ! ENDIF for melting snow
                END IF
                !
                ! ----------------------------------------------------------------------
                !  11.1 Aggregates next
                ! ----------------------------------------------------------------------
                IF(QCF_AGG(I) >  M0.AND.T(I) >  ZERODEGC)THEN
                   ! An approximate calculation of wet bulb temperature
                   ! TEMPW represents AVERAGE supersaturation in the ice
                   ! Strictly speaking, we need to do two melting calculations,
                   ! for the two different wet bulb temperatures.
                   TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0_r8)*CFICEI(I)
                   TEMP7=T(I)-ZERODEGC-TEMPW                                    &
                        &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
                   TEMP7=MAX(TEMP7,0.0_r8)
                   ! End of wet bulb temp formulations.
                   PR02=RHO(I)*MAX(QCF_AGG(I),0.0_r8)*CFICEI(I)*CONSTP(25)*TCGI(I)
                   DPR=TCG(I)*CONSTP(34)*TIMESTEP/LSITER2*                      &
                        &            (CONSTP(27)*CORR2(I)*PR02**CX(24)                     &
                        &         + CONSTP(28)*ROCOR(I)*PR02**CX(25))*RHOR(I)
                   ! Solve implicitly in terms of temperature
                   DPR=TEMP7*(1.0_r8-1.0_r8/(1.0_r8+DPR*LFRCP))/LFRCP
                   DPR=MIN(DPR,QCF_AGG(I))
                   ! Update values of ice and Rain
                   QCF_AGG(I)=QCF_AGG(I)-DPR
                   RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
                   T(I)=T(I)-LFRCP*DPR
                   IF (DPR >  0.0_r8) THEN
                      RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
                      ! Update rain fractions
                      RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                      RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                      RAIN_ICE(I)=                                               &
                           &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                      RAIN_CLEAR(I)=                                             &
                           &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
                   ENDIF
                   ! ENDIF for melting snow
                END IF

             END IF  ! l_it_melting
             !
             ! End of loop over additional iterations
          ENDDO
          !
          ! ----------------------------------------------------------------------
          !  6.1 Homogenous nucleation takes place at temperatures less than THOMO
          ! ----------------------------------------------------------------------
          IF (T(I) <  (ZERODEGC+THOMO)) THEN
             ! Turn all liquid to ice
             QCF_CRY(I)=QCF_CRY(I)+QCL(I)
             T(I)=T(I)+LFRCP*QCL(I)
             QCL(I)=0.0_r8
          END IF
          ! ----------------------------------------------------------------------
          !  6.2 Heteorgenous nucleation occurs for temps less than TNUC deg C
          ! ----------------------------------------------------------------------
          IF (T(I) <  (ZERODEGC+TNUC).AND.                             &
               &       QCF_TOT(I) <  (1.0E5_r8*M0*RHOR(I))) THEN
             ! Calculate number of active ice nucleii.
             DQI=MIN(0.01_r8*EXP(-0.6_r8*(T(I)-ZERODEGC)),1.0E5_r8)
             ! Each nucleus can grow to arbitary mass of M0 kg
             DQI=M0*DQI*RHOR(I)
             ! DQI is amount of nucleation
             DQI=MAX(DQI-QCF_TOT(I),0.0_r8)
             IF (DQI >  0.0_r8) THEN
                ! RHNUC represents how much moisture is available for ice formation.
                RHNUC=(188.92_r8+2.81_r8*(T(I)-ZERODEGC)                         &
                     &       +0.013336_r8*(T(I)-ZERODEGC)**2)*0.01_r8
                RHNUC=MIN(RHNUC,1.0_r8)-0.1_r8
                ! Predict transfer of mass to ice.
                RHNUC=MAX(QSL(I)*RHNUC,QS(I))
                ! DQIL is amount of moisture available
                DQIL=(QCL(I)+QSL(I)-RHNUC)*CFLIQ(I)                        &
                     &            +MAX(Q_ICE(I)-RHNUC,0.0_r8)*AREA_ICE(I)                  &
                     &            +MAX(Q_CLEAR(I)-RHNUC,0.0_r8)*AREA_CLEAR(I)
                DQI=MAX(MIN(DQI,DQIL*LHEAT_CORREC_ICE(I)),0.0_r8)
                QCF_CRY(I)=QCF_CRY(I)+DQI
                ! This comes initially from liquid water
                DQIL=MIN(DQI,QCL(I))
                QCL(I)=QCL(I)-DQIL
                T(I)=T(I)+LFRCP*DQIL
                ! If more moisture is required then nucleation removes from vapour.
                DQI=DQI-DQIL
                T(I)=T(I)+LSRCP*DQI
                Q(I)=Q(I)-DQI
                ! END IFs for nucleation
             END IF
          END IF
          !
          ! ----------------------------------------------------------------------
          !  7   Deposition/Sublimation of snow.
          !      Hallett Mossop process enhances growth.
          ! ----------------------------------------------------------------------
          ! Assume we can't calculate a meaningful width
          IF (Q_ICE(I)  >   QS(I)) THEN
             ! Deposition
             AREA_ICE_1(I) = AREA_ICE(I)
             AREA_ICE_2(I) = 0.0_r8
             Q_ICE_1(I)    = Q_ICE(I)
             Q_ICE_2(I)    = QS(I)       ! Dummy value
          ELSE
             ! Sublimation
             AREA_ICE_1(I) = 0.0_r8
             AREA_ICE_2(I) = AREA_ICE(I)
             Q_ICE_1(I)    = QS(I)       ! Dummy value
             Q_ICE_2(I)    = Q_ICE(I)
          END IF
          IF (AREA_ICE(I)  >   0.0_r8) THEN
             TEMP7 = 0.5_r8*AREA_ICE(I) + (Q_ICE(I)-QS(I)) / WIDTH(I)
             ! Temp7 is now the estimate of the proportion of the gridbox which
             ! contains ice and has local q greater than saturation wrt ice
             IF (TEMP7  >   0.0_r8 .AND. TEMP7  <   AREA_ICE(I)) THEN
                ! Calculate values of q in each region.
                ! These overwrite previous estimates.
                AREA_ICE_1(I) = TEMP7
                AREA_ICE_2(I) = AREA_ICE(I) - AREA_ICE_1(I)
                Q_ICE_1(I)=QS(I) + 0.5_r8 * AREA_ICE_1(I) * WIDTH(I)
                Q_ICE_2(I)=QS(I) - 0.5_r8 * AREA_ICE_2(I) * WIDTH(I)
             END IF
          END IF
       ENDDO
       !
       DO I=1,POINTS
          IF (QCF_CRY(I) >  M0.AND.T(I) <  ZERODEGC) THEN
             TEMP3(I)=QCF_CRY(I)+QCF_AGG(I)
             ! Diffusional parameters
             APLUSB=(APB1-APB2*T(I))*ESI(I)
             APLUSB=APLUSB+(T(I)**3)*P(I)*APB3
             ! Moisture available from subgrid scale calculation
             TEMPW_DEP = QSL(I)*AREA_MIX(I)                             &
                  &                 + MIN(Q_ICE_1(I),QSL(I)) * AREA_ICE_1(I)         &
                  &                 - QS(I) * (AREA_MIX(I) + AREA_ICE_1(I))
             TEMPW_SUB = (Q_ICE_2(I) - QS(I)) * AREA_ICE_2(I)
             PR02=RHO(I)*MAX(QCF_CRY(I),0.0_r8)                            &
                  &                  *CFICEI(I)*CONSTP(5)*TCGCI(I)
             LAMR1=PR02**CX(4)
             LAMR2=PR02**CX(5)
             ! Transfer rate
             DQI=TCGC(I)*CONSTP(6)*T(I)**2*ESI(I)*                      &
                  &       (CONSTP(7)*CORR2(I)*LAMR1+CONSTP(8)*ROCOR(I)*              &
                  &       LAMR2)/(QS(I)*APLUSB*RHO(I))
             DQI_DEP = DQI * TEMPW_DEP
             DQI_SUB = DQI * TEMPW_SUB
             ! Limits depend on whether deposition or sublimation occurs
             IF (DQI_DEP >  0.0_r8) THEN
                ! Deposition is occuring.
                ! Hallett Mossop Enhancement
                IF ( (T(I)-ZERODEGC)  >=  HM_T_MAX) THEN
                   ! Temperature is greater than maximum threshold for HM.
                   HM_RATE=0.0_r8
                ELSEIF ((T(I)-ZERODEGC)  <   HM_T_MAX                    &
                     ! Temperature is between HM thresholds
                     &           .AND. (T(I)-ZERODEGC)  >   HM_T_MIN) THEN
                   HM_RATE=(1.0_r8-EXP( (T(I)-ZERODEGC-HM_T_MAX)*HM_DECAY) ) &
                        &             *HM_NORMALIZE
                ELSE
                   ! Temperature is less than minimum threshold for HM.
                   HM_RATE=EXP( (T(I)-ZERODEGC-HM_T_MIN)*HM_DECAY)
                ENDIF
                ! Calculate enhancement factor for HM process.
                HM_RATE=1.0_r8+HM_RATE*QCL(I)*HM_RQCL
                ! The molecular diffusion to or from the surface is more efficient
                ! when a particle is at a molecular step. This is more likely when
                ! a particle is subliming. For growth, reduce the rate by 10 percent.
                HM_RATE=0.9_r8*HM_RATE
                TEMPW_DEP=TEMPW_DEP*LHEAT_CORREC_ICE(I)
                ! Calculate Transfer. Limit is available moisture.
                IF (CFLIQ(I) >  0.0_r8) THEN
                   ! Add on liquid water contribution to available moisture
                   ! The latent heat correction at this stage becomes very tedious to
                   ! calculate. Freezing the liquid should raise qsat and hence allow
                   ! for some evaporation of liquid itself. It is best to assume that
                   ! no latent heat effect from the freezing should be employed.
                   TEMPW_DEP=TEMPW_DEP+QCL(I)*AREA_MIX(I)/CFLIQ(I)        &
                        &                    +MAX((Q_ICE_1(I)-QSL(I)),0.0_r8)*AREA_ICE_1(I)
                ENDIF
                DQI_DEP=MIN(DQI_DEP*TIMESTEP*HM_RATE,TEMPW_DEP)
             END IF
             IF (DQI_SUB <  0.0_r8) THEN
                ! Sublimation is occuring. Limits are spare moisture capacity and QCF
                ! outside the liquid cloud
                DQI_SUB=MAX(MAX                                          &
                     &                 (DQI_SUB*TIMESTEP,TEMPW_SUB*LHEAT_CORREC_ICE(I)) &
                     &                 ,-(QCF_CRY(I) * AREA_ICE_2(I) * CFICEI(I) ))
             END IF
             ! Adjust ice content
             QCF_CRY(I)=QCF_CRY(I)+DQI_DEP+DQI_SUB
             !
             IF (CFLIQ(I) >  0.0_r8 .AND. AREA_MIX(I) >  0.0_r8               &
                  &            .AND. QCL(I) >  0.0_r8) THEN
                ! Deposition removes some liquid water content.
                !
                ! First estimate of the liquid water removed is explicit
                DQIL=MAX(MIN( DQI_DEP*AREA_MIX(I)                        &
                     &                /(AREA_MIX(I)+AREA_ICE_1(I)),                     &
                     &                QCL(I)*AREA_MIX(I)/CFLIQ(I)),0.0_r8)
             ELSE
                ! Deposition does not remove any liquid water content
                DQIL=0.0_r8
             ENDIF
             ! Adjust liquid content (deposits before vapour by Bergeron Findeison
             !  process).
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
             DQI=DQI_DEP+DQI_SUB-DQIL
             ! Adjust vapour content
             Q(I)=Q(I)-DQI
             T(I)=T(I)+LSRCP*DQI
             ! END IF for QCF >  M0.
          END IF
          !
       ENDDO
       ! ----------------------------------------------------------------------
       !  7.2 Now aggregates
       ! ----------------------------------------------------------------------
       DO I=1,POINTS
          IF (QCF_AGG(I) >  M0.AND.T(I) <  ZERODEGC) THEN
             TEMP3(I)=QCF_CRY(I)+QCF_AGG(I)
             ! Diffusional parameters
             APLUSB=(APB1-APB2*T(I))*ESI(I)
             APLUSB=APLUSB+(T(I)**3)*P(I)*APB3
             ! Moisture available from subgrid scale calculation
             TEMPW_DEP = QSL(I)*AREA_MIX(I)                             &
                  &                 + MIN(Q_ICE_1(I),QSL(I)) * AREA_ICE_1(I)         &
                  &                 - QS(I) * (AREA_MIX(I) + AREA_ICE_1(I))
             TEMPW_SUB = (Q_ICE_2(I) - QS(I)) * AREA_ICE_2(I)
             PR02=RHO(I)*MAX(QCF_AGG(I),0.0_r8)                            &
                  &                  *CFICEI(I)*CONSTP(25)*TCGI(I)
             LAMR1=PR02**CX(24)
             LAMR2=PR02**CX(25)
             ! Transfer rate
             DQI=TCGC(I)*CONSTP(26)*T(I)**2*ESI(I)*                     &
                  &       (CONSTP(27)*CORR2(I)*LAMR1+CONSTP(28)*ROCOR(I)*            &
                  &       LAMR2)/(QS(I)*APLUSB*RHO(I))
             DQI_DEP = DQI * TEMPW_DEP
             DQI_SUB = DQI * TEMPW_SUB
             ! Limits depend on whether deposition or sublimation occurs
             IF (DQI_DEP >  0.0_r8) THEN
                ! Deposition is occuring.
                !
                ! The molecular diffusion to or from the surface is more efficient
                ! when a particle is at a molecular step. This is more likely when
                ! a particle is subliming. For growth, reduce the rate by 10 percent.
                DQI_DEP=0.9_r8*DQI_DEP
                TEMPW_DEP=TEMPW_DEP*LHEAT_CORREC_ICE(I)
                ! Calculate Transfer. Limit is available moisture.
                IF (CFLIQ(I) >  0.0_r8) THEN
                   ! Add on liquid water contribution to available moisture
                   TEMPW_DEP=TEMPW_DEP+QCL(I)*AREA_MIX(I)/CFLIQ(I)        &
                        &                    +MAX((Q_ICE_1(I)-QSL(I)),0.0_r8)*AREA_ICE_1(I)
                ENDIF
                DQI_DEP=MIN(DQI_DEP*TIMESTEP,TEMPW_DEP)
             END IF
             IF (DQI_SUB  <   0.0_r8) THEN
                ! Sublimation is occuring. Limits are spare moisture capacity and QCF
                ! outside liquid cloud
                DQI_SUB=MAX(MAX(                                         &
                     &                 DQI_SUB*TIMESTEP,TEMPW_SUB*LHEAT_CORREC_ICE(I))  &
                     &                    ,-(QCF_AGG(I) * AREA_ICE_2(I) * CFICEI(I) ))
             END IF
             ! Adjust ice content
             QCF_AGG(I)=QCF_AGG(I)+DQI_SUB+DQI_DEP
             !
             IF (CFLIQ(I) >  0.0_r8 .AND. AREA_MIX(I) >  0.0_r8               &
                  &            .AND. QCL(I) >  0.0_r8) THEN
                ! Deposition removes some liquid water content.
                !
                ! First estimate of the liquid water removed is explicit
                DQIL=MAX(MIN( DQI_DEP*AREA_MIX(I)                        &
                     &                /(AREA_MIX(I)+AREA_ICE_1(I)),                     &
                     &                QCL(I)*AREA_MIX(I)/CFLIQ(I)),0.0_r8)
             ELSE
                ! Deposition does not remove any liquid water content
                DQIL=0.0_r8
             ENDIF
             ! Adjust liquid content (deposits before vapour by Bergeron Findeison
             !  process).
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
             DQI=DQI_DEP+DQI_SUB-DQIL
             ! Adjust vapour content
             Q(I)=Q(I)-DQI
             T(I)=T(I)+LSRCP*DQI
             ! END IF for QCF >  M0.
          END IF

          ! ----------------------------------------------------------------------
          !  7.5  Autoconversion of ice crystals to aggregates
          ! ----------------------------------------------------------------------
          IF (L_mcr_qcf2 .AND. QCF_CRY(I) >  M0) THEN

             ! Simple explicit Kessler type param. of autoconversion
             ! Autoconversion rate from Lin et al. (1983)
             ! QCFAUTORATE = 0.005*EXP(0.025*(T(I)-273.16))

             ! Set autoconversion limit to emulate split-ice scheme
             QCFAUTOLIM = (QCF_AGG(I)+QCF_CRY(I))                       &
                  &               *MAX(EXP(-T_SCALING*MAX((T(I)-CTTEMP(I)),0.0_r8)      &
                  &               *MAX(QCF_AGG(I)+QCF_CRY(I),0.0_r8)*QCF0) , 0.0_r8)

             ! Set rate to emulate spilt-ice scheme, i.e. infinite
             QCFAUTORATE = 1.0_r8/TIMESTEP

             QC  = MIN(QCFAUTOLIM,QCF_CRY(I))
             DPR = MIN(QCFAUTORATE*TIMESTEP*(QCF_CRY(I)-QC)             &
                  &                 ,QCF_CRY(I)-QC)

             ! End of calculation of autoconversion amount DPR
             QCF_CRY(I) = QCF_CRY(I)-DPR
             QCF_AGG(I) = QCF_AGG(I)+DPR

          END IF  ! on autoconversion of crystals to agg


          ! ----------------------------------------------------------------------
          !      Transfer processes only active at T less than 0 deg C
          ! ----------------------------------------------------------------------
          IF(T(I) <  ZERODEGC) THEN
             !
             ! ----------------------------------------------------------------------
             !  8   Riming of snow by cloud water -implicit in QCL
             ! ----------------------------------------------------------------------
             IF (QCF_CRY(I) >  M0.AND.QCL(I) >  0.0_r8                       &
                  &         .AND.AREA_MIX(I) >  0.0_r8.AND.CFLIQ(I) >  0.0_r8) THEN
                ! ----------------------------------------------------------------------
                !  8.1 Crystals first
                ! ----------------------------------------------------------------------
                ! Calculate water content of mixed phase region
                QCLNEW=QCL(I)/(CFLIQ(I)+CFLIQ(I)*CONSTP(9)*TCGC(I)*CORR(I) &
                     &                *TIMESTEP*(RHO(I)*QCF_CRY(I)*CFICEI(I)            &
                     &                *CONSTP(5)*TCGCI(I))**CX(6))
                ! Convert to new grid box total water content
                QCLNEW=QCL(I)*AREA_LIQ(I)/CFLIQ(I)+QCLNEW*AREA_MIX(I)
                ! Recalculate water contents
                QCF_CRY(I)=QCF_CRY(I)+(QCL(I)-QCLNEW)
                T(I)=T(I)+LFRCP*(QCL(I)-QCLNEW)
                QCL(I)=QCLNEW
                ! END IF for QCF >  M0.AND.QCL(I) >  0.0
             END IF
             !
             IF (QCF_AGG(I) >  M0.AND.QCL(I) >  0.0_r8                       &
                  &         .AND.AREA_MIX(I) >  0.0_r8.AND.CFLIQ(I) >  0.0_r8) THEN
                ! ----------------------------------------------------------------------
                !  8.2 Aggregates next
                ! ----------------------------------------------------------------------
                ! Calculate water content of mixed phase region
                QCLNEW=QCL(I)/(CFLIQ(I)+CFLIQ(I)*CONSTP(29)*TCG(I)*CORR(I)&
                     &                *TIMESTEP*(RHO(I)*QCF_AGG(I)*CFICEI(I)            &
                     &                *CONSTP(25)*TCGI(I))**CX(26))
                ! Convert to new grid box total water content
                QCLNEW=QCL(I)*AREA_LIQ(I)/CFLIQ(I)+QCLNEW*AREA_MIX(I)
                ! Recalculate water contents
                QCF_AGG(I)=QCF_AGG(I)+(QCL(I)-QCLNEW)
                T(I)=T(I)+LFRCP*(QCL(I)-QCLNEW)
                QCL(I)=QCLNEW
                ! END IF for QCF >  M0.AND.QCL(I) >  0.0
             END IF
             !
             ! ----------------------------------------------------------------------
             !  9   Capture of rain by snow
             ! ----------------------------------------------------------------------
             IF (QCF_CRY(I) >  M0 .AND. RAIN(I)  >   0.0_r8                  &
                  &          .AND. (RAIN_MIX(I)+RAIN_ICE(I)) >  0.0_r8) THEN
                ! ----------------------------------------------------------------------
                !  9.1 Crystals first
                ! ----------------------------------------------------------------------
                ! Calculate velocities
                VR1=CORR(I)*CONSTP(41)/6.0_r8*                                &
                     &              (RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**CX(41)
                VS1=CONSTP(4)*CORR(I)*(RHO(I)*QCF_CRY(I)*CFICEI(I)         &
                     &           *CONSTP(5)*TCGCI(I))**CX(3)
                ! Estimate the mean absolute differences in velocities.
                FV1=MAX(ABS(VR1-VS1),(VR1+VS1)/8.0_r8)
                ! Calculate functions of slope parameter lambda
                LAMR1=(RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**(CX(42))
                LAMS1=(RHO(I)*QCF_CRY(I)*CFICEI(I)                         &
                     &               *CONSTP(5)*TCGCI(I))**(-CX(7))
                LAMFAC1=CONSTP(10)*CONSTP(43)*                             &
                     &                 (LAMR1**CX(43)*LAMS1**CX(8)) +                   &
                     &               CONSTP(11)*CONSTP(44)*                             &
                     &                 (LAMR1**CX(44)*LAMS1**CX(9)) +                   &
                     &               CONSTP(12)*CONSTP(45)*                             &
                     &                 (LAMR1**CX(45)*LAMS1**CX(10))
                ! Calculate transfer
                DPR=TCGC(I)*CONSTP(13)*LAMS1**(-CX(11))*LAMR1**(-CX(46))*FV1* &
                     &       LAMFAC1*TIMESTEP*RHOR(I)*(RAIN_MIX(I)+RAIN_ICE(I))
                DPR=MIN(DPR,RAIN(I)*(DHI(I)*LSITER)*RHOR(I)                &
                     &                   *(RAIN_MIX(I)+RAIN_ICE(I))/RAINFRAC(I))
                ! Adjust ice and rain contents
                QCF_CRY(I)=QCF_CRY(I)+DPR
                RAIN(I)=RAIN(I)-DPR*RHO(I)*DHILSITERR(I)
                T(I)=T(I)+LFRCP*DPR
                RAINFRAC(I)=RAINFRAC(I)*RAIN(I)/                           &
                     &                   (RAIN(I)+DPR*RHO(I)*DHILSITERR(I))
                ! Update rain fractions
                RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                RAIN_ICE(I)=                                               &
                     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                RAIN_CLEAR(I)=                                             &
                     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
                ! Endif for RAIN >  0.0 in capture term
             END IF
             !
             ! ----------------------------------------------------------------------
             !  9.2  Aggregates next
             ! ----------------------------------------------------------------------
             IF (QCF_AGG(I) >  M0 .AND. RAIN(I)  >   0.0_r8                  &
                  &          .AND. (RAIN_MIX(I)+RAIN_ICE(I)) >  0.0_r8) THEN
                ! Calculate velocities
                VR1=CORR(I)*CONSTP(41)/6.0_r8*                                &
                     &              (RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**CX(41)
                VS1=CONSTP(24)*CORR(I)*(RHO(I)*QCF_AGG(I)*CFICEI(I)        &
                     &           *CONSTP(25)*TCGI(I))**CX(23)
                ! Estimate the mean absolute differences in velocities.
                FV1=MAX(ABS(VR1-VS1),(VR1+VS1)/8.0_r8)
                ! Calculate functions of slope parameter lambda
                LAMR1=(RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**(CX(42))
                LAMS1=(RHO(I)*QCF_AGG(I)*CFICEI(I)                         &
                     &             *CONSTP(25)*TCGI(I))**(-CX(27))
                LAMFAC1=CONSTP(30)*CONSTP(43)*                             &
                     &                 (LAMR1**CX(43)*LAMS1**CX(28)) +                  &
                     &               CONSTP(31)*CONSTP(44)*                             &
                     &                 (LAMR1**CX(44)*LAMS1**CX(29)) +                  &
                     &               CONSTP(32)*CONSTP(45)*                             &
                     &                 (LAMR1**CX(45)*LAMS1**CX(30))
                ! Calculate transfer
                DPR=TCG(I)*CONSTP(33)*LAMS1**(-CX(31))*LAMR1**(-CX(46))*FV1*  &
                     &       LAMFAC1*TIMESTEP*RHOR(I)*(RAIN_MIX(I)+RAIN_ICE(I))
                DPR=MIN(DPR,RAIN(I)*(DHI(I)*LSITER)*RHOR(I)                &
                     &                   *(RAIN_MIX(I)+RAIN_ICE(I))/RAINFRAC(I))
                ! Adjust ice and rain contents
                QCF_AGG(I)=QCF_AGG(I)+DPR
                RAIN(I)=RAIN(I)-DPR*RHO(I)*DHILSITERR(I)
                T(I)=T(I)+LFRCP*DPR
                RAINFRAC(I)=RAINFRAC(I)*RAIN(I)/                           &
                     &                   (RAIN(I)+DPR*RHO(I)*DHILSITERR(I))
                ! Update rain fractions
                RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                RAIN_ICE(I)=                                               &
                     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                RAIN_CLEAR(I)=                                             &
                     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
                !      Endif for RAIN >  0.0_r8 in capture term
             END IF
             ! ----------------------------------------------------------------------
             !      End of transfer processes only active at T less than 0 deg C
             ! ----------------------------------------------------------------------
          END IF
          ! ----------------------------------------------------------------------
          !  10  Evaporate melting snow
          ! ----------------------------------------------------------------------
          IF(QCF_CRY(I) >  M0.AND.T(I) >  ZERODEGC)THEN
             ! ----------------------------------------------------------------------
             !  10.1 Crystals first
             ! ----------------------------------------------------------------------
             ! Calculate transfer as a function of QCF, T and specific humidity
             PR02=RHO(I)*MAX(QCF_CRY(I),0.0_r8)*CFICEI(I)*CONSTP(5)*TCGCI(I)
             PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
             DPR=TCGC(I)*CONSTP(6)*T(I)**2*ESW(I)*TIMESTEP*               &
                  &     (CONSTP(7)*CORR2(I)*PR02**CX(4)                              &
                  &      +CONSTP(8)*ROCOR(I)*PR02**CX(5))/(QSL(I)*RHO(I)*PR04)
             ! TEMPW is the subsaturation in the ice region
             TEMPW=AREA_ICE(I)*(QSL(I)-Q_ICE(I))
             DPR=DPR*TEMPW
             DPR=MAX(MIN(DPR,TEMPW*LHEAT_CORREC_LIQ(I)),0.0_r8)
             ! Extra check to see we don't get a negative QCF
             DPR=MIN(DPR,QCF_CRY(I))
             ! Update values of ice and vapour
             QCF_CRY(I)=QCF_CRY(I)-DPR
             Q(I)=Q(I)+DPR
             T(I)=T(I)-DPR*LSRCP
          END IF
          !
          IF(QCF_AGG(I) >  M0.AND.T(I) >  ZERODEGC)THEN
             ! ----------------------------------------------------------------------
             !  10.2 Aggregates next
             ! ----------------------------------------------------------------------
             ! Calculate transfer as a function of QCF, T and specific humidity
             PR02=RHO(I)*MAX(QCF_AGG(I),0.0_r8)*CFICEI(I)*CONSTP(25)*TCGI(I)
             PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
             DPR=TCG(I)*CONSTP(26)*T(I)**2*ESW(I)*TIMESTEP*               &
                  &     (CONSTP(27)*CORR2(I)*PR02**CX(24)                            &
                  &      +CONSTP(28)*ROCOR(I)*PR02**CX(25))/(QSL(I)*RHO(I)*PR04)
             ! TEMPW is the subsaturation in the ice region
             TEMPW=AREA_ICE(I)*(QSL(I)-Q_ICE(I))
             DPR=DPR*TEMPW
             DPR=MAX(MIN(DPR,TEMPW*LHEAT_CORREC_LIQ(I)),0.0_r8)
             ! Extra check to see we don't get a negative QCF
             DPR=MIN(DPR,QCF_AGG(I))
             ! Update values of ice and vapour
             QCF_AGG(I)=QCF_AGG(I)-DPR
             Q(I)=Q(I)+DPR
             T(I)=T(I)-DPR*LSRCP
          END IF
          !
          !
          ! If iterative melting is active then we have already done this term.
          !
          IF (.NOT. l_it_melting) THEN
             ! ----------------------------------------------------------------------
             !  11  Melting of snow. Uses a numerical approximation to the wet bulb
             !      temperature.
             ! ----------------------------------------------------------------------
             IF(QCF_CRY(I) >  M0.AND.T(I) >  ZERODEGC)THEN
                ! ----------------------------------------------------------------------
                !  11.1 Crystals first
                ! ----------------------------------------------------------------------
                ! An approximate calculation of wet bulb temperature
                ! TEMPW represents AVERAGE supersaturation in the ice
                ! Strictly speaking, we need to do two melting calculations,
                ! for the two different wet bulb temperatures.
                TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0_r8)*CFICEI(I)
                TEMP7=T(I)-ZERODEGC-TEMPW                                    &
                     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
                TEMP7=MAX(TEMP7,0.0_r8)
                ! End of wet bulb temp formulations.
                PR02=RHO(I)*QCF_CRY(I)*CFICEI(I)*CONSTP(5)*TCGCI(I)
                DPR=TCGC(I)*CONSTP(14)*TIMESTEP*                             &
                     &            (CONSTP(7)*CORR2(I)*PR02**CX(4)                       &
                     &         + CONSTP(8)*ROCOR(I)*PR02**CX(5))*RHOR(I)
                ! Solve implicitly in terms of temperature
                DPR=TEMP7*(1.0_r8-1.0_r8/(1.0_r8+DPR*LFRCP))/LFRCP
                DPR=MIN(DPR,QCF_CRY(I))
                ! Update values of ice and Rain
                QCF_CRY(I)=QCF_CRY(I)-DPR
                RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
                T(I)=T(I)-LFRCP*DPR
                IF (DPR >  0.0_r8) THEN
                   RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
                   ! Update rain fractions
                   RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                   RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                   RAIN_ICE(I)=                                               &
                        &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                   RAIN_CLEAR(I)=                                             &
                        &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
                ENDIF
                ! ENDIF for melting snow
             END IF
             !
             ! ----------------------------------------------------------------------
             !  11.1 Aggregates next
             ! ----------------------------------------------------------------------
             IF(QCF_AGG(I) >  M0.AND.T(I) >  ZERODEGC)THEN
                ! An approximate calculation of wet bulb temperature
                ! TEMPW represents AVERAGE supersaturation in the ice
                ! Strictly speaking, we need to do two melting calculations,
                ! for the two different wet bulb temperatures.
                TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0_r8)*CFICEI(I)
                TEMP7=T(I)-ZERODEGC-TEMPW                                    &
                     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
                TEMP7=MAX(TEMP7,0.0_r8)
                ! End of wet bulb temp formulations.
                PR02=RHO(I)*QCF_AGG(I)*CFICEI(I)*CONSTP(25)*TCGI(I)
                DPR=TCG(I)*CONSTP(34)*TIMESTEP*                              &
                     &            (CONSTP(27)*CORR2(I)*PR02**CX(24)                     &
                     &         + CONSTP(28)*ROCOR(I)*PR02**CX(25))*RHOR(I)
                ! Solve implicitly in terms of temperature
                DPR=TEMP7*(1.0_r8-1.0_r8/(1.0_r8+DPR*LFRCP))/LFRCP
                DPR=MIN(DPR,QCF_AGG(I))
                ! Update values of ice and Rain
                QCF_AGG(I)=QCF_AGG(I)-DPR
                RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
                T(I)=T(I)-LFRCP*DPR
                IF (DPR >  0.0_r8) THEN
                   RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
                   ! Update rain fractions
                   RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                   RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                   RAIN_ICE(I)=                                               &
                        &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                   RAIN_CLEAR(I)=                                             &
                        &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
                ENDIF
                ! ENDIF for melting snow
             END IF
             !
          END IF  ! l_it_melting
          !
          ! ----------------------------------------------------------------------
          !  12.-1 Evaporate all small rain contents to avoid numerical problems
          ! ----------------------------------------------------------------------
          IF (RAIN(I)  <   ( QCFMIN*DHILSITERR(I)*RHO(I) ) ) THEN
             ! Evaporate all this rain
             DPR=RAIN(I)*RHOR(I)*DHI(I)*LSITER
             T(I)=T(I)-LCRCP*DPR
             Q(I)=Q(I)+DPR
             RAIN(I)=0.0_r8
             ! Update rain fractions
             RAINFRAC(I)=0.0_r8
             RAIN_LIQ(I)=0.0_r8
             RAIN_MIX(I)=0.0_r8
             RAIN_ICE(I)=0.0_r8
             RAIN_CLEAR(I)=0.0_r8
          ENDIF
          ! ----------------------------------------------------------------------
          !      Break loop at this point to use vector multiplication
          ! ----------------------------------------------------------------------
       ENDDO
       !
       ! ----------------------------------------------------------------------
       !  12  Evaporation of rain - implicit in subsaturation
       ! ----------------------------------------------------------------------
       DO I=1,POINTS
          IF(RAIN(I) >  0.0_r8.AND.RAINFRAC(I) >  0.0_r8)THEN
             PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
             ! Define LAMR1 and LAMR2
             LAMR1=RAIN(I)/(CONSTP(42)*CORR(I)*RAINFRAC(I))
             LAMR2=LAMR1**(CX(47)*CX(48))
             LAMR1=LAMR1**(CX(49)*CX(48))
             ! New, consistent evaporation method, with rain fall speed relationship.
             DPR=CONSTP(46)*T(I)**2*ESW(I)*TIMESTEP
             DPR=DPR*( (CONSTP(47)*CORR2(I)*LAMR2)                        &
                  &               + (CONSTP(48)*ROCOR(I)*LAMR1) )
             ! Calculate transfers.
             ! TEMP7 is grid box mean supersaturation. This provides a limit.
             TEMP7=(Q_ICE(I)-QSL(I))*RAIN_ICE(I)                          &
                  &           +(Q_CLEAR(I)-QSL(I))*RAIN_CLEAR(I)
             DPR=DPR*MAX(-TEMP7*LHEAT_CORREC_LIQ(I),0.0_r8)                  &
                  &           /(QSL(I)*RHO(I)*PR04+DPR)
             DPR=DPR*RHO(I)*DHILSITERR(I)
             ! Another limit is on the amount of rain available
             DPR=MIN(DPR,RAIN(I)*                                         &
                  &            (RAIN_ICE(I)+RAIN_CLEAR(I))/RAINFRAC(I))
             ! Update values of rain
             RAIN(I)=RAIN(I)-DPR
             ! Update vapour and temperature
             Q(I)=Q(I)+DPR*DHI(I)*LSITER*RHOR(I)
             T(I)=T(I)-DPR*LCRCP*DHI(I)*LSITER*RHOR(I)
             RAINFRAC(I)=RAINFRAC(I)*RAIN(I)/(RAIN(I)+DPR)
             ! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                                 &
                  &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
             ! END IF for evaporation of rain.
          END IF
       END DO
       !
       ! ----------------------------------------------------------------------
       !  13  Accretion of cloud on rain - implicit in liquid water content
       ! ----------------------------------------------------------------------
       KK=0
       DO I=1,POINTS
          IF(RAINFRAC(I) >  0.0_r8.AND.RAIN(I) >  0.0_r8                       &
               &      .AND.CFLIQ(I) >  0.0_r8)THEN
             KK=KK+1
             ! New accretion formulation.
             TEMP1(KK)=RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I))
             ! QCLNEW is the in cloud water content left
          END IF
       END DO
       KK=0
       DO I=1,POINTS
          IF(RAINFRAC(I) >  0.0_r8.AND.RAIN(I) >  0.0_r8                       &
               &      .AND.CFLIQ(I) >  0.0_r8)THEN
             KK=KK+1
             QCLNEW=QCL(I)/((CFLIQ(I)+CFLIQ(I)*CONSTP(49)*CORR(I)         &
                  &                             *TIMESTEP*TEMP1(KK)**CX(50)))
             ! Convert QCLNEW to a gridbox mean
             ! TEMP7 is overlap of liquid cloud but no rain
             TEMP7=MAX(CFLIQ(I)-RAIN_LIQ(I)-RAIN_MIX(I),0.0_r8)
             QCLNEW=QCLNEW*(RAIN_LIQ(I)+RAIN_MIX(I))+QCL(I)/CFLIQ(I)*TEMP7
             ! Now calculate increments to rain.
             RAIN(I)=RAIN(I)+(QCL(I)-QCLNEW)*RHO(I)*DHILSITERR(I)
             QCL(I)=QCL(I)-(QCL(I)-QCLNEW)
             ! END IF for accretion of cloud on rain.
          END IF
       END DO
       !
       ! ----------------------------------------------------------------------
       !  14  Autoconversion of cloud to rain
       ! ----------------------------------------------------------------------
       ! This looks complicated because there are several routes through
       ! depending on whether we are using a second indirect effect,
       ! have a land or sea point or are using the T3E optimizations. Look
       ! at UMDP 26 to see what is going on.
       !
       KK=0 ! Index counter
       ! ----------------------------------------------------------------------
       !  14.1 Calculate the droplet concentration.
       ! ----------------------------------------------------------------------
       ! There are several possible methods.
       !
       IF (L_AUTOCONV_MURK .AND. L_MURK) THEN ! Use murk aerosol
          DO I=1,POINTS ! Loop b2
             IF (QCL(I) >  0.0_r8.AND.CFLIQ(I) >  0.0_r8) THEN ! Proceed
                KK=KK+1
                ! Convert aerosol mass to droplet number. See subroutine VISBTY
                N_DROP(KK) = MAX(AEROSOL(I)/M0_MURK*1.0E-9_r8, 0.0001_r8)
                !              1.0E-9 converts from ug/kg to kg/kg
                N_DROP(KK)=N0_MURK*N_DROP(KK)**POWER_MURK
             END IF ! Proceed
          END DO ! Loop b2
          !
       ELSE ! Calculate number of droplets either from the sulphate
          !             (and optionally other) aerosol(s) or from a fixed value
          !
          KK=0
          DO I=1, POINTS ! Loop x
             IF (QCL(I)  >   0.0_r8 .AND. CFLIQ(I)  >   0.0_r8) THEN ! Proceed
                IF (L_SEASALT_CCN) THEN
                   SEA_SALT_PTR=I
                ELSE
                   SEA_SALT_PTR=1
                ENDIF
                IF (L_BIOMASS_CCN) THEN
                   BIOMASS_PTR=I
                ELSE
                   BIOMASS_PTR=1
                ENDIF
                IF (L_OCFF_CCN) THEN
                   OCFF_PTR=I
                ELSE
                   OCFF_PTR=1
                ENDIF
                IF (L_biogenic_CCN) THEN
                   biogenic_ptr=I
                ELSE
                   biogenic_ptr=1
                ENDIF
                KK=KK+1
                ! DEPENDS ON: number_droplet
                N_DROP(KK)=NUMBER_DROPLET(&
                     L_USE_SULPHATE_AUTOCONV      , & !  L_AEROSOL_DROPLET, &
                     .FALSE.                      , & !  L_NH42SO4, &
                     !SO4_AIT(I)                  , & ! !AITKEN_SULPHATE  , &
                     SO4_ACC(I)                   , & !  ACCUM_SULPHATE, &
                     SO4_DIS(I)                   , & !  DISS_SULPHATE, &
                     L_SEASALT_CCN                , & !  L_SEASALT_CCN, &
                     SEA_SALT_FILM(SEA_SALT_PTR)  , & !  SEA_SALT_FILM, &
                     SEA_SALT_JET(SEA_SALT_PTR)   , & !  SEA_SALT_JET, &
                     L_biogenic_CCN               , & !  L_BIOGENIC_CCN, &
                     biogenic(biogenic_ptr)       , & !  BIOGENIC , &
                     L_BIOMASS_CCN                , & !  L_BIOMASS_CCN, &
                     BMASS_AGD(BIOMASS_PTR)       , & !  BIOMASS_AGED, &
                     BMASS_CLD(BIOMASS_PTR)       , & !  BIOMASS_CLOUD, &
                     L_OCFF_CCN                   , & !  L_OCFF_CCN, &
                     OCFF_AGD(OCFF_PTR)           , & !  OCFF_AGED, &
                     OCFF_CLD(OCFF_PTR)           , & !  OCFF_CLOUD, &
                     RHO(I)                       , & !  DENSITY_AIR, &
                     SNOW_DEPTH(I)                , & !  SNOW_DEPTH, &
                     LAND_FRACT(I)                , & !  LAND_FRACT, &
                     Ntot_land                    , & !  NTOT_LAND, &
                     Ntot_sea                       ) !  NTOT_SEA   )
             ENDIF ! Proceed
          END DO ! Loop x
          !
       ENDIF
       ! ----------------------------------------------------------------------
       !  14.2 Calculate whether a sufficient concentration of large droplets
       !      is present to allow autoconversion. Calculate its rate.
       ! ----------------------------------------------------------------------
       IF (L_AUTO_DEBIAS) THEN
          DO I=1,POINTS
             T_L(I)=T(I)-(LCRCP*QCL(I))
          END DO
          ! DEPENDS ON: qsat_wat
          !          CALL QSAT_WAT(QSL_TL,T_L,P,POINTS)
          CALL QSAT_WAT (&
               !      Output field
               &  QSL_TL      &!REAL   , INTENT(OUT)  ::  QS(NPNTS) ! SATURATION MIXING RATIO AT TEMPERATURE T AND PRESSURE P (KG/KG)
               !      Input fields
               &, T_L         &!REAL   , INTENT(IN   ) :: T(NPNTS)  !Temperature (K).
               &, P           &!REAL   , INTENT(IN   ) :: P(NPNTS)  !Pressure (Pa).
               !      Array dimensions
               &, POINTS      &!INTEGER, INTENT(IN   ) :: NPNTS  !, INTENT(IN)
               &  )


       ENDIF
       KK=0
       DO I=1,POINTS ! Loop e
          IF (QCL(I) >  0.0_r8.AND.CFLIQ(I) >  0.0_r8) THEN ! Proceed
             KK=KK+1
             ! Non T3E code. This is rather easier to follow.
             !
             ! Calculate inverse of mean droplet size
             R_MEAN(KK)=R_MEAN0*(RHO(I)*QCL(I)                            &
                  &                /(CFLIQ(I)*N_DROP(KK)))**(-1.0_r8/3.0_r8)
             ! Calculate numerical factors
             B_FACTOR(KK)=3.0_r8*R_MEAN(KK)
             A_FACTOR(KK)=N_DROP(KK)*0.5_r8
             ! Calculate droplet number concentration greater than threshold radius
             N_GT_20=(B_FACTOR(KK)**2*A_FACTOR(KK))*R_AUTO**2             &
                  &                                   *EXP(-B_FACTOR(KK)*R_AUTO)     &
                  &            +(2.0_r8*B_FACTOR(KK)*A_FACTOR(KK))*R_AUTO                &
                  &                                   *EXP(-B_FACTOR(KK)*R_AUTO)     &
                  &            +(2.0_r8*A_FACTOR(KK))                                    &
                  &                                   *EXP(-B_FACTOR(KK)*R_AUTO)
             ! Test to see if there is a sufficient concentration of
             ! droplets >20um for autoconversion to proceed (i.e. > N_AUTO):
             IF (N_GT_20  >=  N_AUTO) THEN ! Number
                AUTORATE=CONSTS_AUTO*EC_AUTO                              &
                     &                 *N_DROP(KK)**POWER_DROPLET_AUTO
             ELSE ! Number
                AUTORATE=0.0_r8
             END IF ! Number
             !
             ! ----------------------------------------------------------------------
             !  14.2.1 Optionally de-bias the autoconversion rate.
             ! ----------------------------------------------------------------------
             ! If this option is active, the autoconversion rate is corrected ("de-
             ! biased") based on Wood et al. (Atmos. Res., 65, 109-128, 2002).
             ! This correction should only be used if the second indirect (or "life-
             ! time") effect has been selected (i.e. cloud droplet number is calcul-
             ! ated interactively) and if the diagnostic RHcrit scheme is on (i.e. an
             ! interactive measure of cloud inhomogeneity is available).
             !
             IF (L_AUTO_DEBIAS) THEN

                IF (QCL(I)  <   1.0E-15_r8) THEN  ! Don't bother for very
                   ! small LWCs - avoids
                   AC_FACTOR=1.0_r8                ! AC_FACTOR --> infinity.

                ELSE

                   ALPHA_L=EPSILON*LC*QSL_TL(I)/(R*T_L(I)**2)
                   A_L=1.0_r8/(1.0_r8+(LCRCP*ALPHA_L))
                   SIGMA_S=(1.0_r8-RHCPT(I))*A_L*QSL_TL(I)/SQRT(6.0_r8)
                   G_L=1.15_r8*(POWER_QCL_AUTO-1.0_r8)*SIGMA_S
                   GACB=EXP(-1.0_r8*QCL(I)/G_L)

                   AC_FACTOR=MAX(                                           &
                        &                  (CFLIQ(I)**(POWER_QCL_AUTO-1.0_r8))*1.0_r8/(1.0_r8-GACB) &
                        &                  ,1.0_r8)

                ENDIF

                AUTORATE=AC_FACTOR*AUTORATE

             ENDIF
             !
             ! ----------------------------------------------------------------------
             !  14.3 How much water content can autoconversion remove?
             ! ----------------------------------------------------------------------
             ! Calculate value of local liquid water content at which the droplet
             ! concentration with radii greater than 20um will
             ! fall below a threshold (1000 m-3), and so determine the minimum
             ! liquid water content AUTOLIM: This is a numerical approximation
             ! which ideally could be expressed in terms of N_AUTO and R_AUTO
             ! but isn't.
             AUTOLIM=(6.20E-31_r8*N_DROP(KK)**3)-(5.53E-22_r8*N_DROP(KK)**2)    &
                  &               +(4.54E-13_r8*N_DROP(KK))+(3.71E-6_r8)-(7.59_r8/N_DROP(KK))
             ! Calculate maximum amount of liquid that can be removed from the
             ! grid box
             QC=MIN(AUTOLIM*CFLIQ(I)*RHOR(I),QCL(I))
             ! Calculate autoconversion amount (finally!)
             DPR=MIN(AUTORATE                                             &
                  &              *(RHO(I)*QCL(I)/CFLIQ(I))**(POWER_QCL_AUTO-1.0_r8)     &
                  &              *TIMESTEP*QCL(I)/CORR2(I),QCL(I)-QC)
             ! Update liquid water content and rain
             QCL(I)=QCL(I)-DPR
             RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
             ! Update rain fractions
             IF (DPR >  0.0_r8) THEN ! Rain fraction
                RAINFRAC(I)=MAX(RAINFRAC(I),CFLIQ(I))
                RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
                RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
                RAIN_ICE(I)=                                               &
                     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
                RAIN_CLEAR(I)=                                             &
                     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
             ENDIF ! Rain fraction
          END IF ! Proceed
          ! ----------------------------------------------------------------------
          !  14.4 Remove small rain amounts.
          ! ----------------------------------------------------------------------
          IF (RAIN(I) <= (QCFMIN*DHILSITERR(I)*RHO(I))                   &
                                ! Tidy up rain
               &       .AND.QCL(I) <= 0.0_r8) THEN
             RAINFRAC(I)=0.0_r8
             Q(I)=Q(I)+RAIN(I)*DHI(I)*RHOR(I)*LSITER
             T(I)=T(I)-RAIN(I)*LCRCP*DHI(I)*LSITER*RHOR(I)
             RAIN(I)=0.0_r8
             ! Update rain fractions
             RAINFRAC(I)=0.0_r8
             RAIN_LIQ(I)=0.0_r8
             RAIN_MIX(I)=0.0_r8
             RAIN_ICE(I)=0.0_r8
             RAIN_CLEAR(I)=0.0_r8
          ENDIF ! Tidy up rain
          ! ----------------------------------------------------------------------
          !  15  Now continue the loops over points and iterations.
          ! ----------------------------------------------------------------------
          ! Continue DO loop over points
       END DO ! Loop e or e2
       ! Continue DO loop over iterations
    END DO ! Iters_do1
    !
    ! Copy contents of SNOWT to SNOW, to fall into next layer down
    ! Points_do3
    DO I=1,POINTS
       IF (L_mcr_qcf2) THEN ! two ice prognostics
          QCF(I)  = QCF_AGG(I)
          QCF2(I) = QCF_CRY(I)
          SNOW_CRY(I) = SNOWT_CRY(I)
          SNOW_AGG(I) = SNOWT_AGG(I)
       ELSE ! only one ice prognostic, put all snow in to snow_cry
          QCF(I) = QCF_CRY(I) + QCF_AGG(I)
          SNOW_CRY(I) = 0.0_r8  ! Redundant variable
          SNOW_AGG(I) = SNOWT_CRY(I) + SNOWT_AGG(I)
       ENDIF ! on L_mcr_qcf2
       ! ----------------------------------------------------------------------
       !  16.1 Remove any small amount of QCF which is left over to be tidy.
       !       If QCF is less than QCFMIN and isn't growing
       !       by deposition (assumed to be given by RHCPT) then remove it.
       ! ----------------------------------------------------------------------
       IF (L_mcr_qcf2) THEN

          ! Aggregates
          IF (QCF(I) <  QCFMIN) THEN
             IF (T(I) >  ZERODEGC .OR.                                   &
                  &        (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0_r8)          &
                  &        .OR. QCF(I) <  0.0_r8) THEN
                Q(I)=Q(I)+QCF(I)
                T(I)=T(I)-LSRCP*QCF(I)
                QCF(I)=0.0_r8
             END IF
          END IF
          ! Crystals
          IF (QCF2(I) <  QCFMIN) THEN
             IF (T(I) >  ZERODEGC .OR.                                   &
                  &        (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0_r8)          &
                  &        .OR. QCF2(I) <  0.0_r8) THEN
                Q(I)=Q(I)+QCF2(I)
                T(I)=T(I)-LSRCP*QCF2(I)
                QCF2(I)=0.0_r8
                ! Update ice cloud top temperature
                IF (QCF(I) == 0.0_r8) CTTEMP(I)=T(I)
             END IF
          END IF

       ELSE  ! only one prognostic is active

          IF (QCF(I) <  QCFMIN) THEN
             IF (T(I) >  ZERODEGC .OR.                                   &
                  &        (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0_r8)          &
                  &        .OR. QCF(I) <  0.0_r8) THEN
                Q(I)=Q(I)+QCF(I)
                T(I)=T(I)-LSRCP*QCF(I)
                QCF(I)=0.0_r8
                ! Update ice cloud top temperature
                CTTEMP(I)=T(I)
             END IF
          END IF

       END IF  ! on L_mcr_qcf2

       ! ----------------------------------------------------------------------
       !  16.2 Emergency melting of any excess snow to avoid surface snowfall
       !       at high temperatures.
       ! ----------------------------------------------------------------------
       ! Define rapid melting temperature
       IF (l_it_melting) THEN
          t_rapid_melt = zerodegc + 2.0_r8
       ELSE
          t_rapid_melt = zerodegc
       END IF
       !
       ! Melt SNOW_AGG first
       IF (SNOW_AGG(I)  >   0.0_r8 .AND. T(I)  >   (t_rapid_melt)) THEN
          ! Numerical approximation of wet bulb temperature.
          ! Similar to the melting calculation in section 11.
          TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0_r8)*CFICEI(I)
          TEMP7=T(I)-ZERODEGC-TEMPW*(TW1+TW2*(P(I)-TW3)-TW4*(T(I)-TW5))
          TEMP7=MAX(TEMP7,0.0_r8)
          ! End of wet bulb calculation
          ! Remember that DHI uses the shortened timestep which is why LSITER
          ! appears in the following statements.
          DPR=TEMP7/(LFRCP*LSITER)
          DPR = MIN(DPR,SNOW_AGG(I)*DHI(I)*RHOR(I))
          ! Update values of snow and rain
          SNOW_AGG(I) = SNOW_AGG(I) - DPR*RHO(I)*DHIR(I)
          RAIN(I)=RAIN(I)+DPR*RHO(I)*DHIR(I)
          T(I)=T(I)-LFRCP*DPR*LSITER
          ! END IF for emergency melting of snow
       END IF

       ! If ice crystals prognostic is active, also melt snow_cry
       IF (L_mcr_qcf2 .AND.                                            &
            &      SNOW_CRY(I)  >   0.0_r8 .AND. T(I)  >   (t_rapid_melt)) THEN
          ! Numerical approximation of wet bulb temperature.
          ! Similar to the melting calculation in section 11.
          TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0_r8)*CFICEI(I)
          TEMP7=T(I)-ZERODEGC-TEMPW*(TW1+TW2*(P(I)-TW3)-TW4*(T(I)-TW5))
          TEMP7=MAX(TEMP7,0.0_r8)
          ! End of wet bulb calculation
          ! Remember that DHI uses the shortened timestep which is
          ! why LSITER appears in the following statements.
          DPR = TEMP7/(LFRCP*LSITER)
          DPR = MIN(DPR, SNOW_CRY(I)*DHI(I)*RHOR(I))
          ! Update values of snow and rain
          SNOW_CRY(I) = SNOW_CRY(I) - DPR*RHO(I)*DHIR(I)
          RAIN(I)     = RAIN(I)     + DPR*RHO(I)*DHIR(I)
          T(I)        = T(I)        - LFRCP*DPR*LSITER
       END IF  ! on emergency melting of SNOW_CRY
       !
       ! END DO for tidying up small amounts of ice and snow.
    END DO ! Points_do3
    ! ----------------------------------------------------------------------
    !  17  End of the LSP_ICE subroutine
    ! ----------------------------------------------------------------------
    RETURN
  END SUBROUTINE LSP_ICE



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !*LL  SUBROUTINE LSP_SCAV-----------------------------------------------
  !LL
  !LL  Purpose: Scavenge aerosol by large scale precipitation.
  !LL
  !LL  Model            Modification history from model version 3.0:
  !LL version  Date
  !LL  3.4  15/08/94  New routine. Pete Clark.
  !LL
  !LL  Programming standard: Unified Model Documentation Paper No 3,
  !LL                        Version 7, dated 11/3/93.
  !LL
  !LL  Logical component covered: Part of P26.
  !LL
  !LL  System task:
  !LL
  !LL  Documentation: Unified Model Documentation Paper No 26.
  !*
  !*L  Arguments:---------------------------------------------------------
  SUBROUTINE LSP_SCAV(&
       TIMESTEP   , &!REAL   , INTENT(IN   ) :: TIMESTEP! Input real scalar :- IN Timestep (s).
       POINTS     , &!INTEGER, INTENT(IN   ) :: POINTS  ! Input integer scalar :- IN Number of points to be processed.
       RAIN       , &!REAL   , INTENT(IN   ) :: RAIN(POINTS)  ! Input real arrays :- IN Rate of rainfall in this layer from above! (kg per sq m per s).
       SNOW       , &!REAL   , INTENT(IN   ) :: SNOW(POINTS)  ! IN Rate of snowfall in this layer from above (kg per sq m per s).
       AEROSOL      &!REAL   , INTENT(INOUT) :: AEROSOL(POINTS) ! Updated real arrays :-INOUT Aerosol mixing ratio
       )
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: POINTS  ! Input integer scalar :-
    ! IN Number of points to be processed.
    REAL(KIND=r8)   , INTENT(IN   ) :: TIMESTEP! Input real scalar :- IN Timestep (s).
    REAL(KIND=r8)   , INTENT(IN   ) :: RAIN(POINTS)  ! Input real arrays :-
    ! IN Rate of rainfall in this layer from
    ! above
    ! (kg per sq m per s).
    REAL(KIND=r8)   , INTENT(IN   ) :: SNOW(POINTS)  ! IN Rate of snowfall in this layer from
    !       above
    !       (kg per sq m per s).
    REAL(KIND=r8)    , INTENT(INOUT) :: AEROSOL(POINTS) ! Updated real arrays :-
    ! INOUT Aerosol mixing ratio
    !*L   External subprogram called :-
    !     EXTERNAL None
    !-----------------------------------------------------------------------
    !  Define local scalars.
    !-----------------------------------------------------------------------
    !  (a) Reals effectively expanded to workspace by the Cray (using
    !      vector registers).
    REAL(KIND=r8) :: RRAIN
    REAL(KIND=r8) :: RSNOW
    REAL(KIND=r8),PARAMETER  :: KRAIN  =1.0E-4_r8     ! REAL(KIND=r8) workspace.
    REAL(KIND=r8),PARAMETER  :: KSNOW  =1.0E-4_r8     ! REAL(KIND=r8) workspace.
    !  (b) Others.
    INTEGER :: I       ! Loop counter (horizontal field index).
    !
    ! Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
    RRAIN=KRAIN*TIMESTEP*3600.0_r8
    RSNOW=KSNOW*TIMESTEP*3600.0_r8
    DO I=1,POINTS
       AEROSOL(I)=AEROSOL(I)/(1.0_r8+RRAIN*RAIN(I)+RSNOW*SNOW(I))
    END DO
    RETURN
  END SUBROUTINE LSP_SCAV





  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !+ Calculates constants used in large-scale precipitation scheme.
  SUBROUTINE LSPCON( &
       CX        , &!REAL, INTENT(OUT   ) :: CX(80)
       CONSTP    , &!REAL, INTENT(OUT   ) :: CONSTP(80)
       x1i       , &!REAL, INTENT(IN   ) ::  x1i  ! Intercept of aggregate size distribution
       x1ic      , &!REAL, INTENT(IN   ) ::  x1ic ! Intercept of crystal size distribution
       x1r       , &!REAL, INTENT(IN   ) ::  x1r  ! Intercept of raindrop size distribution
       x2r       , &!REAL, INTENT(IN   ) ::  x2r  ! Scaling parameter of raindrop size distribution
       x4r         )!REAL, INTENT(IN   ) ::  x4r  ! Shape parameter of raindrop size distribution

    IMPLICIT NONE
    !
    !  Description:
    !     Calculates constants used within the LSP_ICE routine.
    !
    !  Method:
    !     Calculate powers, gamma functions and constants from combinations
    !     of physical parameters selected by the comdecks.
    !
    !  Current Code Owner:  Damian Wilson
    !
    !  History:
    !  Version   Date     Comment
    !   5.2      21/11/00 Original code.  Damian Wilson
    !   5.3      21/11/00 Removes duplicate declarations.  A van der Wal
    !   6.1      01/08/04 Include additional constants for ice crystal,
    !                     aggregate, graupel and rain prognostics R.Forbes
    !   6.2      15/08/05 Free format fixes. P.Selwood.
    !   6.2      23/11/05 Pass through precip variables from UMUI
    !                                               Damian Wilson
    !
    !  Code description:
    !   Language: FORTRAN 77 + common extensions.
    !   This code is written to UMDP3 v6 programming standards.
    !
    ! Declarations:
    !
    ! Global variables:
    ! End C_LSPSIZ
    !
    !, Intent(IN)
    REAL(KIND=r8), INTENT(IN   ) ::  x1i  ! Intercept of aggregate size distribution
    REAL(KIND=r8), INTENT(IN   ) ::  x1ic ! Intercept of crystal size distribution
    REAL(KIND=r8), INTENT(IN   ) ::  x1r  ! Intercept of raindrop size distribution
    REAL(KIND=r8), INTENT(IN   ) ::  x2r  ! Scaling parameter of raindrop size distribution
    REAL(KIND=r8), INTENT(IN   ) ::  x4r  ! Shape parameter of raindrop size distribution
    ! Sets up the size of arrays for CX and CONSTP
    REAL(KIND=r8), INTENT(OUT   ) :: CX(80)
    REAL(KIND=r8), INTENT(OUT   ) :: CONSTP(80)


    ! C_LSPDRP start

    ! Microphysics parameters

    ! Drop size distribution for rain: N(D) =  N0 D^m exp(-lambda D)
    ! where N0 = X1R lambda^X2R  and m=X4R
    !     REAL, PARAMETER :: X1R is set in the UMUI
    !     REAL, PARAMETER :: X2R is set in the UMUI
    !     REAL, PARAMETER :: X4R is set in the UMUI

    ! Drop size distribution for graupel: N(D) =  N0 D^m exp(-lambda D)
    ! where N0 = X1G lambda^X2G  and m=X4G
    REAL(KIND=r8), PARAMETER :: X1G=5.E25_r8
    REAL(KIND=r8), PARAMETER :: X2G=-4.0_r8
    REAL(KIND=r8), PARAMETER :: X4G=2.5_r8

    ! Particle size distribution for ice: N(D) = N0 D^m exp(-lambda D)
    ! where N0 = X1I TCG lambda^X2I, m=X4I and TCG=exp(- X3I T[deg C])
    !     REAL(KIND=r8), PARAMETER :: X1I is set in the UMUI
    REAL(KIND=r8), PARAMETER :: X2I=0.0_r8
    REAL(KIND=r8), PARAMETER :: X3I=0.1222_r8
    REAL(KIND=r8), PARAMETER :: X4I=0.0_r8
    !     REAL(KIND=r8), PARAMETER :: X1IC is set in the UMUI
    REAL(KIND=r8), PARAMETER :: X2IC=0.0_r8
    REAL(KIND=r8), PARAMETER :: X3IC=0.1222_r8
    REAL(KIND=r8), PARAMETER :: X4IC=0.0_r8

    ! Mass diameter relationship for graupel:  m(D) = AG D^BG
    REAL(KIND=r8), PARAMETER :: AG=261.8_r8
    REAL(KIND=r8), PARAMETER :: BG=3.0_r8

    ! Mass diameter relationship for ice:  m(D) = AI D^BI
    REAL(KIND=r8), PARAMETER :: AI=0.0444_r8
    REAL(KIND=r8), PARAMETER :: BI=2.1_r8
    REAL(KIND=r8), PARAMETER :: AIC=0.587_r8
    REAL(KIND=r8), PARAMETER :: BIC=2.45_r8

    ! The area diameter relationships are only used if
    ! L_CALCFALL=.TRUE.
    ! Area diameter relationship for ice:  Area(D) = RI D^SI
    REAL(KIND=r8), PARAMETER :: RI=0.131_r8
    REAL(KIND=r8), PARAMETER :: SI=1.88_r8
    REAL(KIND=r8), PARAMETER :: RIC=0.131_r8
    REAL(KIND=r8), PARAMETER :: SIC=1.88_r8

    ! The Best/Reynolds relationships are only used if
    ! L_CALCFALL=.TRUE.
    ! Relationship between Best number and Reynolds number:
    ! Re(D) =EI Be^FI
    REAL(KIND=r8), PARAMETER :: EI=0.2072_r8
    REAL(KIND=r8), PARAMETER :: FI=0.638_r8
    REAL(KIND=r8), PARAMETER :: EIC=0.2072_r8
    REAL(KIND=r8), PARAMETER :: FIC=0.638_r8

    ! The fall speeds of ice particles are only used if
    ! L_CALCFALL=.FALSE.
    ! Fall speed diameter relationships for ice:
    ! vt(D) = CI D^DI
    REAL(KIND=r8), PARAMETER :: CI0=14.3_r8
    REAL(KIND=r8), PARAMETER :: DI0=0.416_r8
    REAL(KIND=r8), PARAMETER :: CIC0=74.5_r8
    REAL(KIND=r8), PARAMETER :: DIC0=0.640_r8

    ! Axial ratio (c-axis divided by a-axis) ESTIMATES. These are not
    ! consistent with those from the area diameter relationships.
    REAL(KIND=r8), PARAMETER :: AR=1.0_r8
    REAL(KIND=r8), PARAMETER :: ARC=1.0_r8

    ! Fall speed diameter relationship for rain: vt(D) = CR D^DR
    REAL(KIND=r8), PARAMETER :: CR=386.8_r8
    REAL(KIND=r8), PARAMETER :: DR=0.67_r8

    ! Fall speed diameter relationship for graupel: vt(D) = CG D^DG
    REAL(KIND=r8), PARAMETER :: CG=253.0_r8
    REAL(KIND=r8), PARAMETER :: DG=0.734_r8

    ! Do we wish to calculate the ice fall velocities?
    ! TRUE if calculate speeds, FALSE if specify speeds
    LOGICAL, PARAMETER :: L_CALCFALL=.TRUE.
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
    REAL(KIND=r8),PARAMETER :: RHO_WATER = 1000.0_r8! DENSITY OF WATER (KG/M3)
    !
    !
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
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8

    ! C_LHEAT end
    ! --------------------------COMDECK C_LSPDIF---------------------------
    ! input variables
    !
    ! Values of reference variables
    REAL(KIND=r8),PARAMETER:: AIR_DENSITY0=1.0_r8             ! kg m-3
    REAL(KIND=r8),PARAMETER:: AIR_VISCOSITY0=1.717E-5_r8      ! kg m-1 s-1
    REAL(KIND=r8),PARAMETER:: AIR_CONDUCTIVITY0=2.40E-2_r8    ! J m-1 s-1 K-1
    REAL(KIND=r8),PARAMETER:: AIR_DIFFUSIVITY0=2.21E-5_r8     ! m2 s-1
    REAL(KIND=r8),PARAMETER:: AIR_PRESSURE0=1.0E5_r8          ! Pa

    ! Values of diffusional growth parameters
    ! Terms in deposition and sublimation
    REAL(KIND=r8),PARAMETER:: APB1=(LC+LF)**2 * EPSILON /(R*AIR_CONDUCTIVITY0)
    REAL(KIND=r8),PARAMETER:: APB2=(LC+LF) / AIR_CONDUCTIVITY0
    REAL(KIND=r8),PARAMETER:: APB3=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)
    ! Terms in evap of melting snow and rain
    REAL(KIND=r8),PARAMETER:: APB4=LC**2*EPSILON/(R*AIR_CONDUCTIVITY0)
    REAL(KIND=r8),PARAMETER:: APB5=LC /AIR_CONDUCTIVITY0
    REAL(KIND=r8),PARAMETER:: APB6=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)

    ! Values of numerical approximation to wet bulb temperature
    ! Numerical fit to wet bulb temperature
    REAL(KIND=r8),PARAMETER:: TW1=1329.31_r8
    REAL(KIND=r8),PARAMETER:: TW2=0.0074615_r8
    REAL(KIND=r8),PARAMETER:: TW3=0.85E5_r8
    ! Numerical fit to wet bulb temperature
    REAL(KIND=r8),PARAMETER:: TW4=40.637_r8
    REAL(KIND=r8),PARAMETER:: TW5=275.0_r8

    ! Ventilation parameters
    REAL(KIND=r8),PARAMETER:: SC=0.6_r8
    ! f(v)  =  VENT_ICE1 + VENT_ICE2  Sc**(1/3) * Re**(1/2)
    REAL(KIND=r8),PARAMETER:: VENT_ICE1=0.65_r8
    REAL(KIND=r8),PARAMETER:: VENT_ICE2=0.44_r8
    ! f(v)  =  VENT_RAIN1 + VENT_RAIN2  Sc**(1/3) * Re**(1/2)
    REAL(KIND=r8),PARAMETER:: VENT_RAIN1=0.78_r8
    REAL(KIND=r8),PARAMETER:: VENT_RAIN2=0.31_r8
    !
    ! Subroutine arguments
    !  Obtain the size of CONSTP and CX
    ! Start C_LSPSIZ
    ! Description: Include file containing idealised forcing options
    ! Author:      R. Forbes
    !
    ! History:
    ! Version  Date      Comment
    ! -------  ----      -------
    !   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
    !   6.2    22/08/05  Include the step size between ice categories.
    !                                                   Damian Wilson

    INTEGER,PARAMETER:: ice_type_offset=20


    ! Counter to print out the contents of CX and CONSTP
    REAL(KIND=r8) :: TEMP 
    ! Forms input to the GAMMAF routine which calculates gamma functions.
    REAL(KIND=r8) :: G1, G2, G3
    REAL(KIND=r8) :: GB1, GB2, GB3
    REAL(KIND=r8) :: GBC1, GBD1, GBDC1
    REAL(KIND=r8) :: GC1, GC2, GC3
    REAL(KIND=r8) :: GD3, GDC3, GD52, GDC52
    REAL(KIND=r8) :: GDR3, GDR4, GDR52
    REAL(KIND=r8) :: GR2, GR4, GR5, GR6
    REAL(KIND=r8) :: GG1, GGB1, GGBD1
    REAL(KIND=r8) :: GG2, GDG52, GDG3, GG3
    ! Represents the gamma function of BI+DI+1 etc.
    ! Fall speed of ice particles parameters
    REAL(KIND=r8) :: CI,DI,CIC,DIC
    !
    ! Function and subroutine calls:
    !EXTERNAL GAMMAF
    !
    !- End of header
    !
    ! Do we need to calculate fall speeds?
    IF (L_CALCFALL) THEN
       ! Define fall speeds
       CI=EI*AIR_VISCOSITY0**(1.0_r8-2.0_r8*FI)                       &
            &            *AIR_DENSITY0**(FI-1.0_r8)                               &
            &            *(2.0_r8*G)**FI*(AI/RI)**FI
       DI=FI*(BI+2.0_r8-SI)-1.0_r8
       CIC=EIC*AIR_VISCOSITY0**(1.0_r8-2.0_r8*FIC)                    &
            &            *AIR_DENSITY0**(FIC-1.0_r8)                              &
            &            *(2.0_r8*G)**FIC*(AIC/RIC)**FIC
       DIC=FIC*(BIC+2.0_r8-SIC)-1.0_r8
    ELSE
       ! Use preset parameters
       CI=CI0
       DI=DI0
       CIC=CIC0
       DIC=DIC0
    ENDIF  ! Calculation of fall speeds
    !
    ! CX values. 1-20 are for the crystal population. 21-40 are for the
    ! aggregate population. 41-60 are for rain.
    ! Crystals
    CX(1)=(BIC+1.0_r8+X4IC-X2IC)/BIC
    CX(2)=-(X4IC+1.0_r8-X2IC)/BIC
    CX(3)=DIC/(BIC+1.0_r8+X4IC-X2IC)
    CX(4)=(2.0_r8+X4IC-X2IC)/(BIC+1.0_r8+X4IC-X2IC)
    CX(5)=(5.0_r8+DIC+2.0_r8*X4IC-2.0_r8*X2IC)*0.5_r8/(BIC+1.0_r8+X4IC-X2IC)
    CX(6)=(3.0_r8+DIC+X4IC-X2IC)/(BIC+1.0_r8+X4IC-X2IC)
    CX(7)=1.0_r8/(X2IC-X4IC-1.0_r8-BIC)
    CX(8)=1.0_r8+X4IC
    CX(9)=2.0_r8+X4IC
    CX(10)=3.0_r8+X4IC
    CX(11)=X2IC
    CX(12)=X3IC
    CX(13)=1.0_r8+X4IC+BIC
    CX(14)=BIC

    ! Aggregates
    CX(23)=DI/(BI+1.0_r8+X4I-X2I)
    CX(24)=(2.0_r8+X4I-X2I)/(BI+1.0_r8+X4I-X2I)
    CX(25)=(5.0_r8+DI+2.0_r8*X4I-2.0_r8*X2I)*0.5_r8/(BI+1.0_r8+X4I-X2I)
    CX(26)=(3.0_r8+DI+X4I-X2I)/(BI+1.0_r8+X4I-X2I)
    CX(27)=1.0_r8/(X2I-X4I-1.0_r8-BI)
    CX(28)=1.0_r8+X4I
    CX(29)=2.0_r8+X4I
    CX(30)=3.0_r8+X4I
    CX(31)=X2I
    CX(32)=X3I
    CX(33)=3.0_r8+X4I+BI
    CX(34)=2.0_r8+X4I+BI
    CX(35)=1.0_r8+X4I+BI
    ! Rain
    CX(41)=DR/(4.0_r8+DR-X2R+X4R)
    CX(42)=1.0_r8/(4.0_r8+DR-X2R+X4R)
    CX(43)=6.0_r8+X4R
    CX(44)=5.0_r8+X4R
    CX(45)=4.0_r8+X4R
    CX(46)=X2R
    CX(47)=2.0_r8+X4R-X2R
    CX(48)=1.0_r8/(4.0_r8+X4R+DR-X2R)
    CX(49)=(DR+5.0_r8)*0.5_r8-X2R+X4R
    CX(50)=(3.0_r8+DR-X2R+X4R)/(4.0_r8+DR-X2R+X4R)
    ! Rain mixing ratio
    CX(51)=DR/(4.0_r8-X2R+X4R)
    CX(52)=1.0_r8/(4.0_r8-X2R+X4R)
    CX(53)=(3.0_r8+DR-X2R+X4R)/(4.0_r8-X2R+X4R)
    !
    ! Graupel

    CX(63)=DG/(BG+1.0_r8+X4G-X2G)
    CX(64)=(2.0_r8+X4G-X2G)/(BG+1.0_r8+X4G-X2G)
    CX(65)=(5.0_r8+DG+2.0_r8*X4G-2.0_r8*X2G)*0.5_r8/(BG+1.0_r8+X4G-X2G)
    CX(66)=(3.0_r8+DG+X4G-X2G)/(BG+1.0_r8+X4G-X2G)
    CX(67)=1.0_r8/(X2G-X4G-1.0_r8-BG)
    CX(68)=1.0_r8+X4G
    CX(69)=2.0_r8+X4G
    CX(70)=3.0_r8+X4G
    CX(71)=X2G

    CX(73)=3.0_r8+X4G+BG
    CX(74)=2.0_r8+X4G+BG
    CX(75)=1.0_r8+X4G+BG

    ! Values 69 to 80 are unused
    ! Define gamma values
    ! Crystals
    TEMP=1.0_r8+X4IC
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GC1)
    TEMP=BIC+1.0_r8+X4IC
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GBC1)
    TEMP=BIC+DIC+1.0_r8+X4IC
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GBDC1)
    TEMP=2.0_r8+X4IC
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GC2)
    TEMP=(DIC+5.0_r8+2.0_r8*X4IC)*0.5_r8
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDC52)
    TEMP=DIC+3.0_r8+X4IC
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDC3)
    TEMP=3.0_r8+X4IC
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GC3)
    ! Aggregates
    TEMP=1.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,G1)
    TEMP=BI+1.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GB1)
    TEMP=BI+2.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GB2)
    TEMP=BI+3.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GB3)
    TEMP=BI+DI+1.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GBD1)
    TEMP=2.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,G2)
    TEMP=(DI+5.0_r8+2.0_r8*X4I)*0.5_r8
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GD52)
    TEMP=DI+3.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GD3)
    TEMP=3.0_r8+X4I
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,G3)
    ! Rain
    TEMP=DR+4.0_r8+X4R
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDR4)
    TEMP=4.0_r8+X4R
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GR4)
    TEMP=6.0_r8+X4R
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GR6)
    TEMP=5.0_r8+X4R
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GR5)
    TEMP=2.0_r8+X4R
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GR2)
    TEMP=(DR+5.0_r8+2.0_r8*X4R)*0.5_r8
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDR52)
    TEMP=DR+3.0_r8+X4R
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDR3)
    ! Graupel
    TEMP=1.0_r8+X4G
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GG1)
    TEMP=BG+1.0_r8+X4G
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GGB1)
    TEMP=BG+DG+1.0_r8+X4G
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GGBD1)
    TEMP=2.0_r8+X4G
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GG2)
    TEMP=(DG+5.0_r8+2.0_r8*X4G)*0.5_r8
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDG52)
    TEMP=DG+3.0_r8+X4G
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GDG3)
    TEMP=3.0_r8+X4G
    ! DEPENDS ON: gammaf
    CALL GAMMAF(TEMP,GG3)

    !
    ! CONSTP values. 1-20 are for the crystal population. 21-40 are for the
    ! aggregate population. 41-60 are for rain.
    ! Crystals
    CONSTP(1)=X1IC
    CONSTP(2)=1.0_r8/GC1
    CONSTP(3)=1.0_r8/(AIC*GBC1)
    CONSTP(4)=CIC*GBDC1/GBC1
    CONSTP(5)=1.0_r8/(AIC*X1IC*GBC1)
    CONSTP(6)=2.0_r8*PI*X1IC
    CONSTP(7)=VENT_ICE1*GC2
    CONSTP(8)=VENT_ICE2*SC**(1.0_r8/3.0_r8)/AIR_VISCOSITY0**0.5_r8      &
         &                 *SQRT(CIC)*GDC52
    CONSTP(9)=PI*0.25_r8*X1IC*CIC*GDC3
    CONSTP(10)=5.0_r8*GC1
    CONSTP(11)=2.0_r8*GC2
    CONSTP(12)=0.25_r8*GC3
    CONSTP(13)=PI**2*RHO_WATER*X1IC*X1R
    CONSTP(14)=2.0_r8*PI*AIR_CONDUCTIVITY0/LF*X1IC
    ! Capacitance relative to spheres of same maximum dimension
    ! Formula depends on value of axial ratio
    CONSTP(15)=ARC
    IF (ARC  >   1.0_r8) THEN
       ! Prolate
       CONSTP(15)=(1.0_r8-(1.0_r8/CONSTP(15))**2)**0.5_r8 /&
            &     LOG( CONSTP(15) + (CONSTP(15)**2-1.0_r8)**0.5_r8 )
    ELSEIF (ARC  ==  1.0_r8) THEN
       ! Spherical
       CONSTP(15)=1.0_r8
    ELSE
       ! Oblate
       CONSTP(15)=(1.0_r8-CONSTP(15)**2)**0.5_r8                      &
            &                    /ASIN((1.0_r8-CONSTP(15)**2)**0.5_r8)
    ENDIF
    ! Now adjust diffusional growth constants for capacitance
    CONSTP(6)=CONSTP(6)*CONSTP(15)
    CONSTP(14)=CONSTP(14)*CONSTP(15)
    !
    ! Values 16 to 23 are unused
    ! Aggregates
    CONSTP(24)=CI*GBD1/GB1
    CONSTP(25)=1.0_r8/(AI*X1I*GB1)
    CONSTP(26)=2.0_r8*PI*X1I
    CONSTP(27)=VENT_ICE1*G2
    CONSTP(28)=VENT_ICE2*SC**(1.0_r8/3.0_r8)/AIR_VISCOSITY0**0.5_r8     &
         &                  *SQRT(CI)*GD52
    CONSTP(29)=PI*0.25_r8*X1I*CI*GD3
    CONSTP(30)=5.0_r8*G1
    CONSTP(31)=2.0_r8*G2
    CONSTP(32)=0.25_r8*G3
    CONSTP(33)=PI**2*RHO_WATER*X1I*X1R
    CONSTP(34)=2.0_r8*PI*AIR_CONDUCTIVITY0/LF*X1I
    ! Capacitance relative to spheres of same maximum dimension
    ! Formula depends on value of axial ratio
    CONSTP(35)=AR
    IF (AR  >   1.0_r8) THEN
       ! Prolate
       CONSTP(35)=(1.0_r8-(1.0_r8/CONSTP(35))**2)**0.5_r8 /              &
            &                    LOG( CONSTP(35) + (CONSTP(35)**2-1.0_r8)**0.5_r8 )
    ELSEIF (AR  ==  1.0_r8) THEN
       ! Spherical
       CONSTP(35)=1.0_r8
    ELSE
       ! Oblate
       CONSTP(35)=(1.0_r8-CONSTP(35)**2)**0.5_r8                      &
            &                    / ASIN((1.0_r8-CONSTP(35)**2)**0.5_r8)
    ENDIF
    ! Now adjust diffusional growth constants for capacitance
    CONSTP(26)=CONSTP(26)*CONSTP(35)
    CONSTP(34)=CONSTP(34)*CONSTP(35)
    !
    CONSTP(36) = GC1*GB3
    CONSTP(37) = 2.0_r8*GC2*GB2
    CONSTP(38) = GC3*GB1
    CONSTP(39) = AI*X1IC*X1I*PI/4.0_r8

    ! Rain
    CONSTP(41)=6.0_r8*CR*GDR4/GR4
    CONSTP(42)=PI*RHO_WATER/6.0_r8*X1R*GDR4*CR
    CONSTP(43)=1.0_r8/120.0_r8*GR6
    CONSTP(44)=1.0_r8/24.0_r8*GR5
    CONSTP(45)=1.0_r8/6.0_r8*GR4
    CONSTP(46)=2.0_r8*PI*X1R
    CONSTP(47)=VENT_RAIN1*GR2
    CONSTP(48)=VENT_RAIN2*SC**(1.0_r8/3.0_r8)/AIR_VISCOSITY0**0.5_r8    &
         &                  *GDR52*SQRT(CR)
    CONSTP(49)=PI*0.25_r8*X1R*CR*GDR3
    CONSTP(50)=1.0_r8/(PI*RHO_WATER*X1R*GR4/6.0_r8)
    ! Values 51 to 60 are unused
    !
    ! Graupel
    CONSTP(64) = CG*GGBD1/GGB1
    CONSTP(65) = 1.0_r8/(AG*X1G*GGB1)
    CONSTP(67) = VENT_RAIN1*GG2
    CONSTP(68) = VENT_RAIN2*SC**(1.0_r8/3.0_r8)/AIR_VISCOSITY0**0.5_r8  &
         &                    *SQRT(CG)*GDG52
    CONSTP(69) = PI*0.25_r8*X1G*CG*GDG3
    CONSTP(74) = 2.0_r8*PI*AIR_CONDUCTIVITY0*X1G/LF
    CONSTP(76) = GG1*GB3
    CONSTP(77) = 2.0_r8*GG2*GB2
    CONSTP(78) = GG3*GB1
    CONSTP(79) = AI*X1G*X1I*PI/4.0_r8

    ! End the subroutine
    RETURN
  END SUBROUTINE LSPCON
  !
  !  SUBROUTINE GAMMAF--------------------------------------------------
  !   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
  !   A POLYNOMIAL APPROXIMATION
  ! --------------------------------------------------------------------
  !





  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !+ Calculates constants used in large-scale precipitation scheme.
  !
  !  SUBROUTINE GAMMAF--------------------------------------------------
  !   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
  !   A POLYNOMIAL APPROXIMATION
  ! --------------------------------------------------------------------
  SUBROUTINE GAMMAF(Y,GAM)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: Y  !, INTENT(IN)
    REAL(KIND=r8), INTENT(OUT  ) :: GAM
    ! Gamma function of Y
    !
    ! LOCAL VARIABLE
    INTEGER :: I,M
    REAL(KIND=r8)    :: GG,G,PARE,X
    ! --------------------------------------------------------------------
    GG=1.0_r8
    M=INT(Y)
    X=Y-M
    IF (M /= 1) THEN
       DO I=1,M-1
          G=Y-I
          GG=GG*G
       END DO
    END IF
    PARE=-0.5748646_r8*X+0.9512363_r8*X*X-0.6998588_r8*X*X*X              &
         +0.4245549_r8*X*X*X*X-0.1010678_r8*X*X*X*X*X+1.0_r8
    GAM=PARE*GG
    RETURN
  END SUBROUTINE GAMMAF
  !



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




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  SUBROUTINE qsat_data()
    !
    !  Local dynamic arrays-------------------------------------------------
    REAL(KIND=r8) :: ES_Aux(0:N+1)    ! TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
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
  SUBROUTINE qsat_wat_data()
    !
    !  Local dynamic arrays-------------------------------------------------
    REAL(KIND=r8) :: ESW_Aux(0:N+1)    ! TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
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
    ! ----------------------------------------------------------------------
    ! SATURATION WATER VAPOUR PRESSURE
    !
    ! VALUES ABOVE AND BELOW 0 DEG C ARE OVER WATER
    !
    ! VALUES BELOW -50 DEC C ARE OUTSIDE FORMAL RANGE OF WMO ADOPTED FORMULA
    ! (SEE ABOVE).  HOWEVER STANDARD PRACTICE TO USE THESE VALUES
    ! ----------------------------------------------------------------------
    ! Note: 0 element is a repeat of 1st element to cater for special case
    !       of low temperatures ( <= T_LOW) for which the array index is
    !       rounded down due to machine precision.
    DATA (ESW_Aux(IESW),IESW=    0, 95) / 0.186905E-01_r8,                      &
         0.186905E-01_r8,0.190449E-01_r8,0.194059E-01_r8,0.197727E-01_r8,0.201462E-01_r8, &
         0.205261E-01_r8,0.209122E-01_r8,0.213052E-01_r8,0.217050E-01_r8,0.221116E-01_r8, &
         0.225252E-01_r8,0.229463E-01_r8,0.233740E-01_r8,0.238090E-01_r8,0.242518E-01_r8, &
         0.247017E-01_r8,0.251595E-01_r8,0.256252E-01_r8,0.260981E-01_r8,0.265795E-01_r8, &
         0.270691E-01_r8,0.275667E-01_r8,0.280733E-01_r8,0.285876E-01_r8,0.291105E-01_r8, &
         0.296429E-01_r8,0.301835E-01_r8,0.307336E-01_r8,0.312927E-01_r8,0.318611E-01_r8, &
         0.324390E-01_r8,0.330262E-01_r8,0.336232E-01_r8,0.342306E-01_r8,0.348472E-01_r8, &
         0.354748E-01_r8,0.361117E-01_r8,0.367599E-01_r8,0.374185E-01_r8,0.380879E-01_r8, &
         0.387689E-01_r8,0.394602E-01_r8,0.401626E-01_r8,0.408771E-01_r8,0.416033E-01_r8, &
         0.423411E-01_r8,0.430908E-01_r8,0.438524E-01_r8,0.446263E-01_r8,0.454124E-01_r8, &
         0.462122E-01_r8,0.470247E-01_r8,0.478491E-01_r8,0.486874E-01_r8,0.495393E-01_r8, &
         0.504057E-01_r8,0.512847E-01_r8,0.521784E-01_r8,0.530853E-01_r8,0.540076E-01_r8, &
         0.549444E-01_r8,0.558959E-01_r8,0.568633E-01_r8,0.578448E-01_r8,0.588428E-01_r8, &
         0.598566E-01_r8,0.608858E-01_r8,0.619313E-01_r8,0.629926E-01_r8,0.640706E-01_r8, &
         0.651665E-01_r8,0.662795E-01_r8,0.674095E-01_r8,0.685570E-01_r8,0.697219E-01_r8, &
         0.709063E-01_r8,0.721076E-01_r8,0.733284E-01_r8,0.745679E-01_r8,0.758265E-01_r8, &
         0.771039E-01_r8,0.784026E-01_r8,0.797212E-01_r8,0.810577E-01_r8,0.824164E-01_r8, &
         0.837971E-01_r8,0.851970E-01_r8,0.866198E-01_r8,0.880620E-01_r8,0.895281E-01_r8, &
         0.910178E-01_r8,0.925278E-01_r8,0.940622E-01_r8,0.956177E-01_r8,0.971984E-01_r8/
    DATA (ESW_Aux(IESW),IESW= 96,190) /                                      &
         0.988051E-01_r8,0.100433E+00_r8,0.102085E+00_r8,0.103764E+00_r8,0.105467E+00_r8, &
         0.107196E+00_r8,0.108953E+00_r8,0.110732E+00_r8,0.112541E+00_r8,0.114376E+00_r8, &
         0.116238E+00_r8,0.118130E+00_r8,0.120046E+00_r8,0.121993E+00_r8,0.123969E+00_r8, &
         0.125973E+00_r8,0.128009E+00_r8,0.130075E+00_r8,0.132167E+00_r8,0.134296E+00_r8, &
         0.136452E+00_r8,0.138642E+00_r8,0.140861E+00_r8,0.143115E+00_r8,0.145404E+00_r8, &
         0.147723E+00_r8,0.150078E+00_r8,0.152466E+00_r8,0.154889E+00_r8,0.157346E+00_r8, &
         0.159841E+00_r8,0.162372E+00_r8,0.164939E+00_r8,0.167545E+00_r8,0.170185E+00_r8, &
         0.172866E+00_r8,0.175584E+00_r8,0.178340E+00_r8,0.181139E+00_r8,0.183977E+00_r8, &
         0.186855E+00_r8,0.189773E+00_r8,0.192737E+00_r8,0.195736E+00_r8,0.198783E+00_r8, &
         0.201875E+00_r8,0.205007E+00_r8,0.208186E+00_r8,0.211409E+00_r8,0.214676E+00_r8, &
         0.217993E+00_r8,0.221355E+00_r8,0.224764E+00_r8,0.228220E+00_r8,0.231728E+00_r8, &
         0.235284E+00_r8,0.238888E+00_r8,0.242542E+00_r8,0.246251E+00_r8,0.250010E+00_r8, &
         0.253821E+00_r8,0.257688E+00_r8,0.261602E+00_r8,0.265575E+00_r8,0.269607E+00_r8, &
         0.273689E+00_r8,0.277830E+00_r8,0.282027E+00_r8,0.286287E+00_r8,0.290598E+00_r8, &
         0.294972E+00_r8,0.299405E+00_r8,0.303904E+00_r8,0.308462E+00_r8,0.313082E+00_r8, &
         0.317763E+00_r8,0.322512E+00_r8,0.327324E+00_r8,0.332201E+00_r8,0.337141E+00_r8, &
         0.342154E+00_r8,0.347234E+00_r8,0.352387E+00_r8,0.357601E+00_r8,0.362889E+00_r8, &
         0.368257E+00_r8,0.373685E+00_r8,0.379194E+00_r8,0.384773E+00_r8,0.390433E+00_r8, &
         0.396159E+00_r8,0.401968E+00_r8,0.407861E+00_r8,0.413820E+00_r8,0.419866E+00_r8/
    DATA (ESW_Aux(IESW),IESW=191,285) /                                      &
         0.425999E+00_r8,0.432203E+00_r8,0.438494E+00_r8,0.444867E+00_r8,0.451332E+00_r8, &
         0.457879E+00_r8,0.464510E+00_r8,0.471226E+00_r8,0.478037E+00_r8,0.484935E+00_r8, &
         0.491920E+00_r8,0.499005E+00_r8,0.506181E+00_r8,0.513447E+00_r8,0.520816E+00_r8, &
         0.528279E+00_r8,0.535835E+00_r8,0.543497E+00_r8,0.551256E+00_r8,0.559113E+00_r8, &
         0.567081E+00_r8,0.575147E+00_r8,0.583315E+00_r8,0.591585E+00_r8,0.599970E+00_r8, &
         0.608472E+00_r8,0.617069E+00_r8,0.625785E+00_r8,0.634609E+00_r8,0.643556E+00_r8, &
         0.652611E+00_r8,0.661782E+00_r8,0.671077E+00_r8,0.680487E+00_r8,0.690015E+00_r8, &
         0.699656E+00_r8,0.709433E+00_r8,0.719344E+00_r8,0.729363E+00_r8,0.739518E+00_r8, &
         0.749795E+00_r8,0.760217E+00_r8,0.770763E+00_r8,0.781454E+00_r8,0.792258E+00_r8, &
         0.803208E+00_r8,0.814309E+00_r8,0.825528E+00_r8,0.836914E+00_r8,0.848422E+00_r8, &
         0.860086E+00_r8,0.871891E+00_r8,0.883837E+00_r8,0.895944E+00_r8,0.908214E+00_r8, &
         0.920611E+00_r8,0.933175E+00_r8,0.945890E+00_r8,0.958776E+00_r8,0.971812E+00_r8, &
         0.985027E+00_r8,0.998379E+00_r8,0.101193E+01_r8,0.102561E+01_r8,0.103949E+01_r8, &
         0.105352E+01_r8,0.106774E+01_r8,0.108213E+01_r8,0.109669E+01_r8,0.111144E+01_r8, &
         0.112636E+01_r8,0.114148E+01_r8,0.115676E+01_r8,0.117226E+01_r8,0.118791E+01_r8, &
         0.120377E+01_r8,0.121984E+01_r8,0.123608E+01_r8,0.125252E+01_r8,0.126919E+01_r8, &
         0.128604E+01_r8,0.130309E+01_r8,0.132036E+01_r8,0.133782E+01_r8,0.135549E+01_r8, &
         0.137339E+01_r8,0.139150E+01_r8,0.140984E+01_r8,0.142839E+01_r8,0.144715E+01_r8, &
         0.146616E+01_r8,0.148538E+01_r8,0.150482E+01_r8,0.152450E+01_r8,0.154445E+01_r8/
    DATA (ESW_Aux(IESW),IESW=286,380) /                                      &
         0.156459E+01_r8,0.158502E+01_r8,0.160564E+01_r8,0.162654E+01_r8,0.164766E+01_r8, &
         0.166906E+01_r8,0.169070E+01_r8,0.171257E+01_r8,0.173473E+01_r8,0.175718E+01_r8, &
         0.177984E+01_r8,0.180282E+01_r8,0.182602E+01_r8,0.184951E+01_r8,0.187327E+01_r8, &
         0.189733E+01_r8,0.192165E+01_r8,0.194629E+01_r8,0.197118E+01_r8,0.199636E+01_r8, &
         0.202185E+01_r8,0.204762E+01_r8,0.207372E+01_r8,0.210010E+01_r8,0.212678E+01_r8, &
         0.215379E+01_r8,0.218109E+01_r8,0.220873E+01_r8,0.223668E+01_r8,0.226497E+01_r8, &
         0.229357E+01_r8,0.232249E+01_r8,0.235176E+01_r8,0.238134E+01_r8,0.241129E+01_r8, &
         0.244157E+01_r8,0.247217E+01_r8,0.250316E+01_r8,0.253447E+01_r8,0.256617E+01_r8, &
         0.259821E+01_r8,0.263064E+01_r8,0.266341E+01_r8,0.269661E+01_r8,0.273009E+01_r8, &
         0.276403E+01_r8,0.279834E+01_r8,0.283302E+01_r8,0.286811E+01_r8,0.290358E+01_r8, &
         0.293943E+01_r8,0.297571E+01_r8,0.301236E+01_r8,0.304946E+01_r8,0.308702E+01_r8, &
         0.312491E+01_r8,0.316326E+01_r8,0.320208E+01_r8,0.324130E+01_r8,0.328092E+01_r8, &
         0.332102E+01_r8,0.336162E+01_r8,0.340264E+01_r8,0.344407E+01_r8,0.348601E+01_r8, &
         0.352838E+01_r8,0.357118E+01_r8,0.361449E+01_r8,0.365834E+01_r8,0.370264E+01_r8, &
         0.374737E+01_r8,0.379265E+01_r8,0.383839E+01_r8,0.388469E+01_r8,0.393144E+01_r8, &
         0.397876E+01_r8,0.402656E+01_r8,0.407492E+01_r8,0.412378E+01_r8,0.417313E+01_r8, &
         0.422306E+01_r8,0.427359E+01_r8,0.432454E+01_r8,0.437617E+01_r8,0.442834E+01_r8, &
         0.448102E+01_r8,0.453433E+01_r8,0.458816E+01_r8,0.464253E+01_r8,0.469764E+01_r8, &
         0.475321E+01_r8,0.480942E+01_r8,0.486629E+01_r8,0.492372E+01_r8,0.498173E+01_r8/
    DATA (ESW_Aux(IESW),IESW=381,475) /                                      &
         0.504041E+01_r8,0.509967E+01_r8,0.515962E+01_r8,0.522029E+01_r8,0.528142E+01_r8, &
         0.534337E+01_r8,0.540595E+01_r8,0.546912E+01_r8,0.553292E+01_r8,0.559757E+01_r8, &
         0.566273E+01_r8,0.572864E+01_r8,0.579532E+01_r8,0.586266E+01_r8,0.593075E+01_r8, &
         0.599952E+01_r8,0.606895E+01_r8,0.613918E+01_r8,0.621021E+01_r8,0.628191E+01_r8, &
         0.635433E+01_r8,0.642755E+01_r8,0.650162E+01_r8,0.657639E+01_r8,0.665188E+01_r8, &
         0.672823E+01_r8,0.680532E+01_r8,0.688329E+01_r8,0.696198E+01_r8,0.704157E+01_r8, &
         0.712206E+01_r8,0.720319E+01_r8,0.728534E+01_r8,0.736829E+01_r8,0.745204E+01_r8, &
         0.753671E+01_r8,0.762218E+01_r8,0.770860E+01_r8,0.779588E+01_r8,0.788408E+01_r8, &
         0.797314E+01_r8,0.806318E+01_r8,0.815408E+01_r8,0.824599E+01_r8,0.833874E+01_r8, &
         0.843254E+01_r8,0.852721E+01_r8,0.862293E+01_r8,0.871954E+01_r8,0.881724E+01_r8, &
         0.891579E+01_r8,0.901547E+01_r8,0.911624E+01_r8,0.921778E+01_r8,0.932061E+01_r8, &
         0.942438E+01_r8,0.952910E+01_r8,0.963497E+01_r8,0.974181E+01_r8,0.984982E+01_r8, &
         0.995887E+01_r8,0.100690E+02_r8,0.101804E+02_r8,0.102926E+02_r8,0.104063E+02_r8, &
         0.105210E+02_r8,0.106367E+02_r8,0.107536E+02_r8,0.108719E+02_r8,0.109912E+02_r8, &
         0.111116E+02_r8,0.112333E+02_r8,0.113563E+02_r8,0.114804E+02_r8,0.116056E+02_r8, &
         0.117325E+02_r8,0.118602E+02_r8,0.119892E+02_r8,0.121197E+02_r8,0.122513E+02_r8, &
         0.123844E+02_r8,0.125186E+02_r8,0.126543E+02_r8,0.127912E+02_r8,0.129295E+02_r8, &
         0.130691E+02_r8,0.132101E+02_r8,0.133527E+02_r8,0.134965E+02_r8,0.136415E+02_r8, &
         0.137882E+02_r8,0.139361E+02_r8,0.140855E+02_r8,0.142366E+02_r8,0.143889E+02_r8/
    DATA (ESW_Aux(IESW),IESW=476,570) /                                      &
         0.145429E+02_r8,0.146982E+02_r8,0.148552E+02_r8,0.150135E+02_r8,0.151735E+02_r8, &
         0.153349E+02_r8,0.154979E+02_r8,0.156624E+02_r8,0.158286E+02_r8,0.159965E+02_r8, &
         0.161659E+02_r8,0.163367E+02_r8,0.165094E+02_r8,0.166838E+02_r8,0.168597E+02_r8, &
         0.170375E+02_r8,0.172168E+02_r8,0.173979E+02_r8,0.175806E+02_r8,0.177651E+02_r8, &
         0.179513E+02_r8,0.181394E+02_r8,0.183293E+02_r8,0.185210E+02_r8,0.187146E+02_r8, &
         0.189098E+02_r8,0.191066E+02_r8,0.193059E+02_r8,0.195065E+02_r8,0.197095E+02_r8, &
         0.199142E+02_r8,0.201206E+02_r8,0.203291E+02_r8,0.205397E+02_r8,0.207522E+02_r8, &
         0.209664E+02_r8,0.211831E+02_r8,0.214013E+02_r8,0.216221E+02_r8,0.218448E+02_r8, &
         0.220692E+02_r8,0.222959E+02_r8,0.225250E+02_r8,0.227559E+02_r8,0.229887E+02_r8, &
         0.232239E+02_r8,0.234614E+02_r8,0.237014E+02_r8,0.239428E+02_r8,0.241872E+02_r8, &
         0.244335E+02_r8,0.246824E+02_r8,0.249332E+02_r8,0.251860E+02_r8,0.254419E+02_r8, &
         0.256993E+02_r8,0.259600E+02_r8,0.262225E+02_r8,0.264873E+02_r8,0.267552E+02_r8, &
         0.270248E+02_r8,0.272970E+02_r8,0.275719E+02_r8,0.278497E+02_r8,0.281295E+02_r8, &
         0.284117E+02_r8,0.286965E+02_r8,0.289843E+02_r8,0.292743E+02_r8,0.295671E+02_r8, &
         0.298624E+02_r8,0.301605E+02_r8,0.304616E+02_r8,0.307650E+02_r8,0.310708E+02_r8, &
         0.313803E+02_r8,0.316915E+02_r8,0.320064E+02_r8,0.323238E+02_r8,0.326437E+02_r8, &
         0.329666E+02_r8,0.332928E+02_r8,0.336215E+02_r8,0.339534E+02_r8,0.342885E+02_r8, &
         0.346263E+02_r8,0.349666E+02_r8,0.353109E+02_r8,0.356572E+02_r8,0.360076E+02_r8, &
         0.363606E+02_r8,0.367164E+02_r8,0.370757E+02_r8,0.374383E+02_r8,0.378038E+02_r8/
    DATA (ESW_Aux(IESW),IESW=571,665) /                                      &
         0.381727E+02_r8,0.385453E+02_r8,0.389206E+02_r8,0.392989E+02_r8,0.396807E+02_r8, &
         0.400663E+02_r8,0.404555E+02_r8,0.408478E+02_r8,0.412428E+02_r8,0.416417E+02_r8, &
         0.420445E+02_r8,0.424502E+02_r8,0.428600E+02_r8,0.432733E+02_r8,0.436900E+02_r8, &
         0.441106E+02_r8,0.445343E+02_r8,0.449620E+02_r8,0.453930E+02_r8,0.458280E+02_r8, &
         0.462672E+02_r8,0.467096E+02_r8,0.471561E+02_r8,0.476070E+02_r8,0.480610E+02_r8, &
         0.485186E+02_r8,0.489813E+02_r8,0.494474E+02_r8,0.499170E+02_r8,0.503909E+02_r8, &
         0.508693E+02_r8,0.513511E+02_r8,0.518376E+02_r8,0.523277E+02_r8,0.528232E+02_r8, &
         0.533213E+02_r8,0.538240E+02_r8,0.543315E+02_r8,0.548437E+02_r8,0.553596E+02_r8, &
         0.558802E+02_r8,0.564046E+02_r8,0.569340E+02_r8,0.574672E+02_r8,0.580061E+02_r8, &
         0.585481E+02_r8,0.590963E+02_r8,0.596482E+02_r8,0.602041E+02_r8,0.607649E+02_r8, &
         0.613311E+02_r8,0.619025E+02_r8,0.624779E+02_r8,0.630574E+02_r8,0.636422E+02_r8, &
         0.642324E+02_r8,0.648280E+02_r8,0.654278E+02_r8,0.660332E+02_r8,0.666426E+02_r8, &
         0.672577E+02_r8,0.678771E+02_r8,0.685034E+02_r8,0.691328E+02_r8,0.697694E+02_r8, &
         0.704103E+02_r8,0.710556E+02_r8,0.717081E+02_r8,0.723639E+02_r8,0.730269E+02_r8, &
         0.736945E+02_r8,0.743681E+02_r8,0.750463E+02_r8,0.757309E+02_r8,0.764214E+02_r8, &
         0.771167E+02_r8,0.778182E+02_r8,0.785246E+02_r8,0.792373E+02_r8,0.799564E+02_r8, &
         0.806804E+02_r8,0.814109E+02_r8,0.821479E+02_r8,0.828898E+02_r8,0.836384E+02_r8, &
         0.843922E+02_r8,0.851525E+02_r8,0.859198E+02_r8,0.866920E+02_r8,0.874712E+02_r8, &
         0.882574E+02_r8,0.890486E+02_r8,0.898470E+02_r8,0.906525E+02_r8,0.914634E+02_r8/
    DATA (ESW_Aux(IESW),IESW=666,760) /                                      &
         0.922814E+02_r8,0.931048E+02_r8,0.939356E+02_r8,0.947736E+02_r8,0.956171E+02_r8, &
         0.964681E+02_r8,0.973246E+02_r8,0.981907E+02_r8,0.990605E+02_r8,0.999399E+02_r8, &
         0.100825E+03_r8,0.101718E+03_r8,0.102617E+03_r8,0.103523E+03_r8,0.104438E+03_r8, &
         0.105358E+03_r8,0.106287E+03_r8,0.107221E+03_r8,0.108166E+03_r8,0.109115E+03_r8, &
         0.110074E+03_r8,0.111039E+03_r8,0.112012E+03_r8,0.112992E+03_r8,0.113981E+03_r8, &
         0.114978E+03_r8,0.115981E+03_r8,0.116993E+03_r8,0.118013E+03_r8,0.119041E+03_r8, &
         0.120077E+03_r8,0.121122E+03_r8,0.122173E+03_r8,0.123234E+03_r8,0.124301E+03_r8, &
         0.125377E+03_r8,0.126463E+03_r8,0.127556E+03_r8,0.128657E+03_r8,0.129769E+03_r8, &
         0.130889E+03_r8,0.132017E+03_r8,0.133152E+03_r8,0.134299E+03_r8,0.135453E+03_r8, &
         0.136614E+03_r8,0.137786E+03_r8,0.138967E+03_r8,0.140158E+03_r8,0.141356E+03_r8, &
         0.142565E+03_r8,0.143781E+03_r8,0.145010E+03_r8,0.146247E+03_r8,0.147491E+03_r8, &
         0.148746E+03_r8,0.150011E+03_r8,0.151284E+03_r8,0.152571E+03_r8,0.153862E+03_r8, &
         0.155168E+03_r8,0.156481E+03_r8,0.157805E+03_r8,0.159137E+03_r8,0.160478E+03_r8, &
         0.161832E+03_r8,0.163198E+03_r8,0.164569E+03_r8,0.165958E+03_r8,0.167348E+03_r8, &
         0.168757E+03_r8,0.170174E+03_r8,0.171599E+03_r8,0.173037E+03_r8,0.174483E+03_r8, &
         0.175944E+03_r8,0.177414E+03_r8,0.178892E+03_r8,0.180387E+03_r8,0.181886E+03_r8, &
         0.183402E+03_r8,0.184930E+03_r8,0.186463E+03_r8,0.188012E+03_r8,0.189571E+03_r8, &
         0.191146E+03_r8,0.192730E+03_r8,0.194320E+03_r8,0.195930E+03_r8,0.197546E+03_r8, &
         0.199175E+03_r8,0.200821E+03_r8,0.202473E+03_r8,0.204142E+03_r8,0.205817E+03_r8/
    DATA (ESW_Aux(IESW),IESW=761,855) /                                      &
         0.207510E+03_r8,0.209216E+03_r8,0.210928E+03_r8,0.212658E+03_r8,0.214398E+03_r8, &
         0.216152E+03_r8,0.217920E+03_r8,0.219698E+03_r8,0.221495E+03_r8,0.223297E+03_r8, &
         0.225119E+03_r8,0.226951E+03_r8,0.228793E+03_r8,0.230654E+03_r8,0.232522E+03_r8, &
         0.234413E+03_r8,0.236311E+03_r8,0.238223E+03_r8,0.240151E+03_r8,0.242090E+03_r8, &
         0.244049E+03_r8,0.246019E+03_r8,0.248000E+03_r8,0.249996E+03_r8,0.252009E+03_r8, &
         0.254037E+03_r8,0.256077E+03_r8,0.258128E+03_r8,0.260200E+03_r8,0.262284E+03_r8, &
         0.264384E+03_r8,0.266500E+03_r8,0.268629E+03_r8,0.270779E+03_r8,0.272936E+03_r8, &
         0.275110E+03_r8,0.277306E+03_r8,0.279509E+03_r8,0.281734E+03_r8,0.283966E+03_r8, &
         0.286227E+03_r8,0.288494E+03_r8,0.290780E+03_r8,0.293083E+03_r8,0.295398E+03_r8, &
         0.297737E+03_r8,0.300089E+03_r8,0.302453E+03_r8,0.304841E+03_r8,0.307237E+03_r8, &
         0.309656E+03_r8,0.312095E+03_r8,0.314541E+03_r8,0.317012E+03_r8,0.319496E+03_r8, &
         0.322005E+03_r8,0.324527E+03_r8,0.327063E+03_r8,0.329618E+03_r8,0.332193E+03_r8, &
         0.334788E+03_r8,0.337396E+03_r8,0.340025E+03_r8,0.342673E+03_r8,0.345329E+03_r8, &
         0.348019E+03_r8,0.350722E+03_r8,0.353440E+03_r8,0.356178E+03_r8,0.358938E+03_r8, &
         0.361718E+03_r8,0.364513E+03_r8,0.367322E+03_r8,0.370160E+03_r8,0.373012E+03_r8, &
         0.375885E+03_r8,0.378788E+03_r8,0.381691E+03_r8,0.384631E+03_r8,0.387579E+03_r8, &
         0.390556E+03_r8,0.393556E+03_r8,0.396563E+03_r8,0.399601E+03_r8,0.402646E+03_r8, &
         0.405730E+03_r8,0.408829E+03_r8,0.411944E+03_r8,0.415083E+03_r8,0.418236E+03_r8, &
         0.421422E+03_r8,0.424632E+03_r8,0.427849E+03_r8,0.431099E+03_r8,0.434365E+03_r8/
    DATA (ESW_Aux(IESW),IESW=856,950) /                                      &
         0.437655E+03_r8,0.440970E+03_r8,0.444301E+03_r8,0.447666E+03_r8,0.451038E+03_r8, &
         0.454445E+03_r8,0.457876E+03_r8,0.461316E+03_r8,0.464790E+03_r8,0.468281E+03_r8, &
         0.471798E+03_r8,0.475342E+03_r8,0.478902E+03_r8,0.482497E+03_r8,0.486101E+03_r8, &
         0.489741E+03_r8,0.493408E+03_r8,0.497083E+03_r8,0.500804E+03_r8,0.504524E+03_r8, &
         0.508290E+03_r8,0.512074E+03_r8,0.515877E+03_r8,0.519717E+03_r8,0.523566E+03_r8, &
         0.527462E+03_r8,0.531367E+03_r8,0.535301E+03_r8,0.539264E+03_r8,0.543245E+03_r8, &
         0.547265E+03_r8,0.551305E+03_r8,0.555363E+03_r8,0.559462E+03_r8,0.563579E+03_r8, &
         0.567727E+03_r8,0.571905E+03_r8,0.576102E+03_r8,0.580329E+03_r8,0.584576E+03_r8, &
         0.588865E+03_r8,0.593185E+03_r8,0.597514E+03_r8,0.601885E+03_r8,0.606276E+03_r8, &
         0.610699E+03_r8,0.615151E+03_r8,0.619625E+03_r8,0.624140E+03_r8,0.628671E+03_r8, &
         0.633243E+03_r8,0.637845E+03_r8,0.642465E+03_r8,0.647126E+03_r8,0.651806E+03_r8, &
         0.656527E+03_r8,0.661279E+03_r8,0.666049E+03_r8,0.670861E+03_r8,0.675692E+03_r8, &
         0.680566E+03_r8,0.685471E+03_r8,0.690396E+03_r8,0.695363E+03_r8,0.700350E+03_r8, &
         0.705381E+03_r8,0.710444E+03_r8,0.715527E+03_r8,0.720654E+03_r8,0.725801E+03_r8, &
         0.730994E+03_r8,0.736219E+03_r8,0.741465E+03_r8,0.746756E+03_r8,0.752068E+03_r8, &
         0.757426E+03_r8,0.762819E+03_r8,0.768231E+03_r8,0.773692E+03_r8,0.779172E+03_r8, &
         0.784701E+03_r8,0.790265E+03_r8,0.795849E+03_r8,0.801483E+03_r8,0.807137E+03_r8, &
         0.812842E+03_r8,0.818582E+03_r8,0.824343E+03_r8,0.830153E+03_r8,0.835987E+03_r8, &
         0.841871E+03_r8,0.847791E+03_r8,0.853733E+03_r8,0.859727E+03_r8,0.865743E+03_r8/
    DATA (ESW_Aux(IESW),IESW=951,1045) /                                     &
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
    DATA (ESW_Aux(IESW),IESW=1046,1140) /                                    &
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
    DATA (ESW_Aux(IESW),IESW=1141,1235) /                                    &
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
    DATA (ESW_Aux(IESW),IESW=1236,1330) /                                    &
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
    DATA (ESW_Aux(IESW),IESW=1331,1425) /                                    &
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
    DATA (ESW_Aux(IESW),IESW=1426,1520) /                                    &
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
    DATA (ESW_Aux(IESW),IESW=1521,1552) /                                    &
         0.218439E+05_r8,0.219439E+05_r8,0.220440E+05_r8,0.221449E+05_r8,0.222457E+05_r8, &
         0.223473E+05_r8,0.224494E+05_r8,0.225514E+05_r8,0.226542E+05_r8,0.227571E+05_r8, &
         0.228606E+05_r8,0.229646E+05_r8,0.230687E+05_r8,0.231734E+05_r8,0.232783E+05_r8, &
         0.233839E+05_r8,0.234898E+05_r8,0.235960E+05_r8,0.237027E+05_r8,0.238097E+05_r8, &
         0.239173E+05_r8,0.240254E+05_r8,0.241335E+05_r8,0.242424E+05_r8,0.243514E+05_r8, &
         0.244611E+05_r8,0.245712E+05_r8,0.246814E+05_r8,0.247923E+05_r8,0.249034E+05_r8, &
         0.250152E+05_r8,0.250152E+05_r8/
    ESW2=ESW_Aux
  END SUBROUTINE qsat_wat_data

  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Saturation Specific Humidity/mixing Scheme(Qsat):Vapour to Liquid/Ice
  !
  ! Subroutine Interface:
  SUBROUTINE QSAT_mix ( &
                                !      Output field
       &  QmixS &  ! REAL(KIND=r8), INTENT(out)   :: QmixS(npnts)  ! Output Saturation mixing ratio or saturation
                                ! specific humidity at temperature T and pressure
                                ! P (kg/kg).
                                !      Input fields
       &, T &      ! REAL(KIND=r8), INTENT(in)  :: T(npnts)      !  Temperature (K).
       &, P &      ! REAL(KIND=r8), INTENT(in)  :: P(npnts)      !  Pressure (Pa).
                                !      Array dimensions
       &, NPNTS &  !INTEGER, INTENT(in) :: npnts    ! Points (=horizontal dimensions) being processed by qSAT scheme
                                !      logical control
       &, lq_mix & !LOGICAL, INTENT(in)  :: lq_mix      !  .true. return qsat as a mixing ratio
                                !  .false. return qsat as a specific humidity
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
    INTEGER      , INTENT(in)  :: npnts    ! Points (=horizontal dimensions) being processed by qSAT scheme.
    REAL(KIND=r8), INTENT(in)  :: T(npnts)      !  Temperature (K).
    REAL(KIND=r8), INTENT(in)  :: P(npnts)      !  Pressure (Pa).
    !    INTEGER(KIND=i8), INTENT(in)  :: mskant(npnts)

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
    QmixS=0.0_r8
    !-----------------------------------------------------------------------
    !
    ! loop over points
    !
    DO i = 1, npnts-1,2
       !       IF(mskant(i) == 1_i8)THEN

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
       itable = INT(atable)
       atable = atable - itable

       TT_p1 = MAX(T_LOW,T(I+1))
       TT_p1 = MIN(T_HIGH,TT_p1)
       ATABLE_p1 = (TT_p1 - T_LOW + DELTA_T) * R_DELTA_T
       ITABLE_p1 = INT(ATABLE_p1)
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
       !       END IF
    END DO  ! Npnts_do_1

    II = I
    DO I = II, NPNTS
       !       IF(mskant(i) == 1_i8)THEN

       fsubw = 1.0_r8 + 1.0E-8_r8*P(I)*( 4.5_r8 +                               &
            &    6.0E-4_r8*( T(I) - zerodegC )*( T(I) - zerodegC ) )
       !
       TT = MAX(T_low,T(I))
       TT = MIN(T_high,TT)
       atable = (TT - T_low + delta_T) * R_delta_T
       itable = INT(atable)
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
       !       END IF
    END DO  ! Npnts_do_1

    RETURN
  END SUBROUTINE QSAT_mix
  ! ======================================================================



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Saturation Specific Humidity Scheme (Qsat_Wat): Vapour to Liquid.
  ! Subroutine Interface:
  SUBROUTINE qsat_wat_mix ( &
       !      Output field
       &  QmixS            &! Real, intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
       !      Input fields
       &, T                &! Real, intent(in)  :: T(npnts) !  Temperature (K).
       &, P                &! Real, intent(in)  :: P(npnts) !  Pressure (Pa).  
       !      Array dimensions
       &, NPNTS            &! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
       !      logical control
       &, lq_mix           &! logical, intent(in)  :: lq_mix
                                ! .true. return qsat as a mixing ratio
                                ! .false. return qsat as a specific humidity
       &  )

    IMPLICIT NONE

    ! Purpose:
    !   Returns a saturation specific humidity or mixing ratio given a
    !   temperature and pressure using the saturation vapour pressure
    !   calculated using the Goff-Gratch formulae, adopted by the WMO as
    !   taken from Landolt-Bornstein, 1987 Numerical Data and Functional
    !   Relationships in Science and Technolgy. Group V/vol 4B meteorology.
    !   Phyiscal and Chemical properties or air, P35.
    !
    !   Values in the lookup table are over water above and below 0 deg C.
    !
    !   Note : For vapour pressure over water this formula is valid for
    !   temperatures between 373K and 223K. The values for saturated vapour
    !   over water in the lookup table below are out of the lower end of
    !   this range. However it is standard WMO practice to use the formula
    !   below its accepted range for use with the calculation of dew points
    !   in the upper atmosphere
    !
    ! Method:
    !   Uses lookup tables to find eSAT, calculates qSAT directly from that.
    !
    ! Current Owner of Code: Owner of large-scale cloud code.
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
    !
    !  Subroutine Arguments:------------------------------------------------
    !
    ! arguments with intent in. ie: input variables.

    INTEGER, INTENT(in) :: npnts
    ! Points (=horizontal dimensions) being processed by qSAT scheme.

    REAL(KIND=r8), INTENT(in)  :: T(npnts) !  Temperature (K).
    REAL(KIND=r8), INTENT(in)  :: P(npnts) !  Pressure (Pa).  

    LOGICAL, INTENT(in)  :: lq_mix
    !  .true. return qsat as a mixing ratio
    !  .false. return qsat as a specific humidity

    ! arguments with intent out

    REAL(KIND=r8), INTENT(out)   ::  QmixS(npnts)
    ! Output Saturation mixing ratio or saturation specific
    ! humidity at temperature T and pressure P (kg/kg).

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

    !-----------------------------------------------------------------------
    ! loop over points
    QmixS=0.0_r8; fsubw=0.0_r8; fsubw_p1 =0.0_r8 ; atable=0.0_r8; atable_p1 =0.0_r8
    DO i = 1, npnts-1, 2

       !      Compute the factor that converts from sat vapour pressure in a
       !      pure water system to sat vapour pressure in air, fsubw.
       !      This formula is taken from equation A4.7 of Adrian Gill's book:
       !      Atmosphere-Ocean Dynamics. Note that his formula works in terms
       !      of pressure in mb and temperature in Celsius, so conversion of
       !      units leads to the slightly different equation used here.

       fsubw    = 1.0_r8 + 1.0E-8_r8*P(I)   * ( 4.5_r8 +                        &
            6.0E-4_r8*( T(I)   - zerodegC ) * ( T(I) - zerodegC ) )
       fsubw_p1 = 1.0_r8 + 1.0E-8_r8*P(I+1) * ( 4.5_r8 +                        &
            6.0E-4_r8*( T(I+1) - zerodegC ) * ( T(I+1) - zerodegC ) )

       !      Use the lookup table to find saturated vapour pressure, and store
       !      it in qmixs.

       TT = MAX(T_low,T(I))
       TT = MIN(T_high,TT)
       atable = (TT - T_low + delta_T) * R_delta_T
       itable = INT(atable)
       atable = atable - itable

       TT_p1 = MAX(T_low,T(I+1))
       TT_p1 = MIN(T_high,TT_p1)
       atable_p1 = (TT_p1 - T_low + delta_T) * R_delta_T
       itable_p1 = INT(atable_p1)
       atable_p1 = atable_p1 - itable_p1

       QmixS(I)   = (1.0_r8 - atable)    * ESW2(itable)    +                  &
            atable*ESW2(itable+1)
       QmixS(I+1) = (1.0_r8 - atable_p1) * ESW2(itable_p1) +                  &
            atable_p1*ESW2(itable_p1+1)

       !      Multiply by fsubw to convert to saturated vapour pressure in air
       !      (equation A4.6 of Adrian Gill's book).

       QmixS(I)   = QmixS(I)   * fsubw
       QmixS(I+1) = QmixS(I+1) * fsubw_p1

       !      Now form the accurate expression for qmixs, which is a rearranged
       !      version of equation A4.3 of Gill's book.

       !-----------------------------------------------------------------------
       ! For mixing ratio,  rsat = epsilon *e/(p-e)
       ! e - saturation vapour pressure
       ! Note applying the fix to qsat for specific humidity at low pressures
       ! is not possible, this implies mixing ratio qsat tends to infinity.
       ! If the pressure is very low then the mixing ratio value may become
       ! very large.
       !-----------------------------------------------------------------------
       IF (lq_mix) THEN
          QmixS(I)   = ( epsilon*QmixS(I) ) /                            &
               ( MAX(P(I),  1.1_r8*QmixS(I))   - QmixS(I) )
          QmixS(I+1) = ( epsilon*QmixS(I+1) ) /                          &
               ( MAX(P(I+1),1.1_r8*QmixS(I+1)) - QmixS(I+1) )
          !-----------------------------------------------------------------------
          ! For specific humidity,   qsat = epsilon*e/(p-(1-epsilon)e)
          !
          ! Note that at very low pressure we apply a fix, to prevent a
          ! singularity (qsat tends to 1. kg/kg).
          !-----------------------------------------------------------------------
       ELSE
          QmixS(I)   = ( epsilon*QmixS(I) ) /                            &
               ( MAX(P(I),  QmixS(I))   - one_minus_epsilon*QmixS(I) )
          QmixS(I+1) = ( epsilon*QmixS(I+1) ) /                          &
               ( MAX(P(I+1),QmixS(I+1)) - one_minus_epsilon*QmixS(I+1))
       ENDIF  ! test on lq_mix
    END DO  ! Npnts_do_1

    ii = i
    DO i = ii, npnts
       fsubw = 1.0_r8 + 1.0E-8_r8*P(I)*( 4.5_r8 +                               &
            &    6.0E-4_r8*( T(I) - zerodegC )*( T(I) - zerodegC ) )

       TT = MAX(T_low,T(I))
       TT = MIN(T_high,TT)
       atable = (TT - T_low + delta_T) * R_delta_T
       itable = INT(atable)
       atable = atable - itable

       QmixS(I) = (1.0_r8 - atable)*ESW2(itable) + atable*ESW2(itable+1)
       QmixS(I) = QmixS(I) * fsubw

       IF (lq_mix) THEN
          QmixS(I) = ( epsilon*QmixS(I) ) /                              &
               &              ( MAX(P(I),1.1_r8*QmixS(i)) - QmixS(I) )
       ELSE
          QmixS(I) = ( epsilon*QmixS(I) ) /                              &
               &               ( MAX(P(I),QmixS(I)) - one_minus_epsilon*QmixS(I) )
       ENDIF  ! test on lq_mix
    END DO  ! Npnts_do_1

    RETURN
  END SUBROUTINE qsat_wat_mix
  ! ======================================================================




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Saturation Specific Humidity Scheme (Qsat_Wat): Vapour to Liquid.
  ! Subroutine Interface:
  SUBROUTINE QSAT_WAT (&
       !      Output field
       &  QS      &!REAL   , INTENT(OUT)  ::  QS(NPNTS) ! SATURATION MIXING RATIO AT TEMPERATURE T AND PRESSURE P (KG/KG)
       !      Input fields
       &, T       &!REAL   , INTENT(IN   ) :: T(NPNTS)  !Temperature (K).
       &, P       &!REAL   , INTENT(IN   ) :: P(NPNTS)  !Pressure (Pa).
       !      Array dimensions
       &, NPNTS   &!INTEGER, INTENT(IN   ) :: NPNTS  !, INTENT(IN)
       &  )
    !
    IMPLICIT NONE
    !
    ! Purpose:
    !   RETURNS A SATURATION MIXING RATIO GIVEN A TEMPERATURE AND PRESSURE
    !   USING SATURATION VAPOUR PRESSURES CALCULATED USING THE GOFF-GRATCH
    !   FORMULAE, ADOPTED BY THE WMO AS TAKEN FROM LANDOLT-BORNSTEIN, 1987
    !   NUMERICAL DATA AND FUNCTIONAL RELATIONSHIPS IN SCIENCE AND
    !   TECHNOLOGY. GROUP V/VOL 4B METEOROLOGY. PHYSICAL AND CHEMICAL
    !   PROPERTIES OF AIR, P35
    !
    !   VALUES IN THE LOOKUP TABLE ARE OVER WATER ABOVE AND BELOW 0 DEG C.
    !
    !   NOTE : FOR VAPOUR PRESSURE OVER WATER THIS FORMULA IS VALID FOR
    !   TEMPERATURES BETWEEN 373K AND 223K.  THE VALUES FOR SATURATED VAPOUR
    !   OVER WATER IN THE LOOKUP TABLE BELOW ARE OUT OF THE LOWER END OF
    !   THIS RANGE.  HOWEVER IT IS STANDARD WMO PRACTICE TO USE THE FORMULA
    !   BELOW ITS ACCEPTED RANGE FOR USE WITH THE CALCULATION OF DEW POINTS
    !   IN THE UPPER ATMOSPHERE
    !
    ! Method:
    !   Uses lookup tables to find eSAT, calculates qSAT directly from that.
    !
    ! Current Owner of Code: OWNER OF LARGE_SCALE CLOUD CODE.
    !
    ! History:
    !   Version  Date    Comment
    !    5.1   14-03-00  Correct 2B (optimized) version for use in UM5.x .
    !                                                        A.C.Bushell
    !    6.2   03-02-06  Moved to a71_1a.  P.Selwood
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
    !
    !  Subroutine Arguments:------------------------------------------------
    !
    ! arguments with intent in. ie: input variables.
    !
    INTEGER, INTENT(IN   ) :: NPNTS                   !, INTENT(IN)
    !       Points (=horizontal dimensions) being processed by qSAT scheme.
    !
    REAL(KIND=r8)   , INTENT(IN   ) :: T(NPNTS)  !       Temperature (K).
    REAL(KIND=r8)   , INTENT(IN   ) :: P(NPNTS)  !       Pressure (Pa).
    !
    ! arguments with intent out. ie: output variables.
    !
    REAL(KIND=r8)   , INTENT(OUT)  ::  QS(NPNTS) ! SATURATION MIXING RATIO AT TEMPERATURE T AND PRESSURE P (KG/KG)
    !
    !  Local parameters and other physical constants------------------------
    !REAL(KIND=r8) ONE_MINUS_EPSILON  ! ONE MINUS THE RATIO OF THE MOLECULAR
    !                               WEIGHTS OF WATER AND DRY AIR
    !
    !REAL(KIND=r8) T_LOW        ! LOWEST TEMPERATURE FOR WHICH LOOK-UP TABLE OF
    !                         SATURATION WATER VAPOUR PRESSURE IS VALID (K)
    !
    !REAL(KIND=r8) T_HIGH       ! HIGHEST TEMPERATURE FOR WHICH LOOK-UP TABLE OF
    !                         SATURATION WATER VAPOUR PRESSURE IS VALID (K)
    !
    !REAL(KIND=r8) DELTA_T      ! TEMPERATURE INCREMENT OF THE LOOK-UP
    ! TABLE OF SATURATION VAPOUR PRESSURES
    !
    !INTEGER N         ! SIZE OF LOOK-UP TABLE OF SATURATION
    ! WATER VAPOUR PRESSURES
    !
    REAL(KIND=r8), PARAMETER :: ONE_MINUS_EPSILON = 1.0_r8 - EPSILON
    REAL(KIND=r8), PARAMETER :: T_LOW = 183.15_r8
    REAL(KIND=r8), PARAMETER :: T_HIGH = 338.15_r8
    REAL(KIND=r8), PARAMETER :: DELTA_T = 0.1_r8
    INTEGER, PARAMETER :: N = INT(((T_HIGH - T_LOW + (DELTA_T*0.5_r8))/DELTA_T) + 1.0_r8)
    ! gives N=1551
    !
    !  Local scalars--------------------------------------------------------
    !
    !  (a) Scalars effectively expanded to workspace by the Cray (using
    !      vector registers).
    !
    INTEGER ITABLE    ! WORK VARIABLES
    !
    REAL(KIND=r8) ATABLE       ! WORK VARIABLES
    !
    !     VARIABLES INTRODUCED BY DLR.
    !
    REAL(KIND=r8) FSUBW        ! FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE
    !                         IN A PURE WATER SYSTEM TO SAT VAPOUR PRESSURE
    !                         IN AIR.
    REAL(KIND=r8) TT
    !
    !  (b) Others.
    INTEGER I         ! LOOP COUNTERS
    !
    INTEGER IES       ! LOOP COUNTER FOR DATA STATEMENT LOOK-UP TABLE
    !
    !  Local dynamic arrays-------------------------------------------------
    REAL(KIND=r8) ES(0:N+1)    ! TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
    !                       - SET BY DATA STATEMENT CALCULATED FROM THE
    !                         GOFF-GRATCH FORMULAE AS TAKEN FROM LANDOLT-
    !                         BORNSTEIN, 1987 NUMERICAL DATA AND FUNCTIONAL
    !                         RELATIONSHIPS IN SCIENCE AND TECHNOLOGY.
    !                         GROUP V/ VOL 4B METEOROLOGY. PHYSICAL AND
    !                         CHEMICAL PROPERTIES OF AIR, P35
    !
    !  External subroutine calls: ------------------------------------------
    !     EXTERNAL None
    !- End of Header
    !
    ! ==Main Block==--------------------------------------------------------
    ! Subroutine structure :
    ! NO SIGNIFICANT STRUCTURE
    ! ----------------------------------------------------------------------
    ! SATURATION WATER VAPOUR PRESSURE
    !
    ! VALUES ABOVE AND BELOW 0 DEG C ARE OVER WATER
    !
    ! VALUES BELOW -50 DEC C ARE OUTSIDE FORMAL RANGE OF WMO ADOPTED FORMULA
    ! (SEE ABOVE).  HOWEVER STANDARD PRACTICE TO USE THESE VALUES
    ! ----------------------------------------------------------------------
    ! Note: 0 element is a repeat of 1st element to cater for special case
    !       of low temperatures ( <= T_LOW) for which the array index is
    !       rounded down due to machine precision.
    DATA (ES(IES),IES=    0, 95) / 0.186905E-01_r8,                      &
         &0.186905E-01_r8,0.190449E-01_r8,0.194059E-01_r8,0.197727E-01_r8,0.201462E-01_r8, &
         &0.205261E-01_r8,0.209122E-01_r8,0.213052E-01_r8,0.217050E-01_r8,0.221116E-01_r8, &
         &0.225252E-01_r8,0.229463E-01_r8,0.233740E-01_r8,0.238090E-01_r8,0.242518E-01_r8, &
         &0.247017E-01_r8,0.251595E-01_r8,0.256252E-01_r8,0.260981E-01_r8,0.265795E-01_r8, &
         &0.270691E-01_r8,0.275667E-01_r8,0.280733E-01_r8,0.285876E-01_r8,0.291105E-01_r8, &
         &0.296429E-01_r8,0.301835E-01_r8,0.307336E-01_r8,0.312927E-01_r8,0.318611E-01_r8, &
         &0.324390E-01_r8,0.330262E-01_r8,0.336232E-01_r8,0.342306E-01_r8,0.348472E-01_r8, &
         &0.354748E-01_r8,0.361117E-01_r8,0.367599E-01_r8,0.374185E-01_r8,0.380879E-01_r8, &
         &0.387689E-01_r8,0.394602E-01_r8,0.401626E-01_r8,0.408771E-01_r8,0.416033E-01_r8, &
         &0.423411E-01_r8,0.430908E-01_r8,0.438524E-01_r8,0.446263E-01_r8,0.454124E-01_r8, &
         &0.462122E-01_r8,0.470247E-01_r8,0.478491E-01_r8,0.486874E-01_r8,0.495393E-01_r8, &
         &0.504057E-01_r8,0.512847E-01_r8,0.521784E-01_r8,0.530853E-01_r8,0.540076E-01_r8, &
         &0.549444E-01_r8,0.558959E-01_r8,0.568633E-01_r8,0.578448E-01_r8,0.588428E-01_r8, &
         &0.598566E-01_r8,0.608858E-01_r8,0.619313E-01_r8,0.629926E-01_r8,0.640706E-01_r8, &
         &0.651665E-01_r8,0.662795E-01_r8,0.674095E-01_r8,0.685570E-01_r8,0.697219E-01_r8, &
         &0.709063E-01_r8,0.721076E-01_r8,0.733284E-01_r8,0.745679E-01_r8,0.758265E-01_r8, &
         &0.771039E-01_r8,0.784026E-01_r8,0.797212E-01_r8,0.810577E-01_r8,0.824164E-01_r8, &
         &0.837971E-01_r8,0.851970E-01_r8,0.866198E-01_r8,0.880620E-01_r8,0.895281E-01_r8, &
         &0.910178E-01_r8,0.925278E-01_r8,0.940622E-01_r8,0.956177E-01_r8,0.971984E-01_r8/
    DATA (ES(IES),IES= 96,190) /                                      &
         &0.988051E-01_r8,0.100433E+00_r8,0.102085E+00_r8,0.103764E+00_r8,0.105467E+00_r8, &
         &0.107196E+00_r8,0.108953E+00_r8,0.110732E+00_r8,0.112541E+00_r8,0.114376E+00_r8, &
         &0.116238E+00_r8,0.118130E+00_r8,0.120046E+00_r8,0.121993E+00_r8,0.123969E+00_r8, &
         &0.125973E+00_r8,0.128009E+00_r8,0.130075E+00_r8,0.132167E+00_r8,0.134296E+00_r8, &
         &0.136452E+00_r8,0.138642E+00_r8,0.140861E+00_r8,0.143115E+00_r8,0.145404E+00_r8, &
         &0.147723E+00_r8,0.150078E+00_r8,0.152466E+00_r8,0.154889E+00_r8,0.157346E+00_r8, &
         &0.159841E+00_r8,0.162372E+00_r8,0.164939E+00_r8,0.167545E+00_r8,0.170185E+00_r8, &
         &0.172866E+00_r8,0.175584E+00_r8,0.178340E+00_r8,0.181139E+00_r8,0.183977E+00_r8, &
         &0.186855E+00_r8,0.189773E+00_r8,0.192737E+00_r8,0.195736E+00_r8,0.198783E+00_r8, &
         &0.201875E+00_r8,0.205007E+00_r8,0.208186E+00_r8,0.211409E+00_r8,0.214676E+00_r8, &
         &0.217993E+00_r8,0.221355E+00_r8,0.224764E+00_r8,0.228220E+00_r8,0.231728E+00_r8, &
         &0.235284E+00_r8,0.238888E+00_r8,0.242542E+00_r8,0.246251E+00_r8,0.250010E+00_r8, &
         &0.253821E+00_r8,0.257688E+00_r8,0.261602E+00_r8,0.265575E+00_r8,0.269607E+00_r8, &
         &0.273689E+00_r8,0.277830E+00_r8,0.282027E+00_r8,0.286287E+00_r8,0.290598E+00_r8, &
         &0.294972E+00_r8,0.299405E+00_r8,0.303904E+00_r8,0.308462E+00_r8,0.313082E+00_r8, &
         &0.317763E+00_r8,0.322512E+00_r8,0.327324E+00_r8,0.332201E+00_r8,0.337141E+00_r8, &
         &0.342154E+00_r8,0.347234E+00_r8,0.352387E+00_r8,0.357601E+00_r8,0.362889E+00_r8, &
         &0.368257E+00_r8,0.373685E+00_r8,0.379194E+00_r8,0.384773E+00_r8,0.390433E+00_r8, &
         &0.396159E+00_r8,0.401968E+00_r8,0.407861E+00_r8,0.413820E+00_r8,0.419866E+00_r8/
    DATA (ES(IES),IES=191,285) /                                      &
         &0.425999E+00_r8,0.432203E+00_r8,0.438494E+00_r8,0.444867E+00_r8,0.451332E+00_r8, &
         &0.457879E+00_r8,0.464510E+00_r8,0.471226E+00_r8,0.478037E+00_r8,0.484935E+00_r8, &
         &0.491920E+00_r8,0.499005E+00_r8,0.506181E+00_r8,0.513447E+00_r8,0.520816E+00_r8, &
         &0.528279E+00_r8,0.535835E+00_r8,0.543497E+00_r8,0.551256E+00_r8,0.559113E+00_r8, &
         &0.567081E+00_r8,0.575147E+00_r8,0.583315E+00_r8,0.591585E+00_r8,0.599970E+00_r8, &
         &0.608472E+00_r8,0.617069E+00_r8,0.625785E+00_r8,0.634609E+00_r8,0.643556E+00_r8, &
         &0.652611E+00_r8,0.661782E+00_r8,0.671077E+00_r8,0.680487E+00_r8,0.690015E+00_r8, &
         &0.699656E+00_r8,0.709433E+00_r8,0.719344E+00_r8,0.729363E+00_r8,0.739518E+00_r8, &
         &0.749795E+00_r8,0.760217E+00_r8,0.770763E+00_r8,0.781454E+00_r8,0.792258E+00_r8, &
         &0.803208E+00_r8,0.814309E+00_r8,0.825528E+00_r8,0.836914E+00_r8,0.848422E+00_r8, &
         &0.860086E+00_r8,0.871891E+00_r8,0.883837E+00_r8,0.895944E+00_r8,0.908214E+00_r8, &
         &0.920611E+00_r8,0.933175E+00_r8,0.945890E+00_r8,0.958776E+00_r8,0.971812E+00_r8, &
         &0.985027E+00_r8,0.998379E+00_r8,0.101193E+01_r8,0.102561E+01_r8,0.103949E+01_r8, &
         &0.105352E+01_r8,0.106774E+01_r8,0.108213E+01_r8,0.109669E+01_r8,0.111144E+01_r8, &
         &0.112636E+01_r8,0.114148E+01_r8,0.115676E+01_r8,0.117226E+01_r8,0.118791E+01_r8, &
         &0.120377E+01_r8,0.121984E+01_r8,0.123608E+01_r8,0.125252E+01_r8,0.126919E+01_r8, &
         &0.128604E+01_r8,0.130309E+01_r8,0.132036E+01_r8,0.133782E+01_r8,0.135549E+01_r8, &
         &0.137339E+01_r8,0.139150E+01_r8,0.140984E+01_r8,0.142839E+01_r8,0.144715E+01_r8, &
         &0.146616E+01_r8,0.148538E+01_r8,0.150482E+01_r8,0.152450E+01_r8,0.154445E+01_r8/
    DATA (ES(IES),IES=286,380) /                                      &
         &0.156459E+01_r8,0.158502E+01_r8,0.160564E+01_r8,0.162654E+01_r8,0.164766E+01_r8, &
         &0.166906E+01_r8,0.169070E+01_r8,0.171257E+01_r8,0.173473E+01_r8,0.175718E+01_r8, &
         &0.177984E+01_r8,0.180282E+01_r8,0.182602E+01_r8,0.184951E+01_r8,0.187327E+01_r8, &
         &0.189733E+01_r8,0.192165E+01_r8,0.194629E+01_r8,0.197118E+01_r8,0.199636E+01_r8, &
         &0.202185E+01_r8,0.204762E+01_r8,0.207372E+01_r8,0.210010E+01_r8,0.212678E+01_r8, &
         &0.215379E+01_r8,0.218109E+01_r8,0.220873E+01_r8,0.223668E+01_r8,0.226497E+01_r8, &
         &0.229357E+01_r8,0.232249E+01_r8,0.235176E+01_r8,0.238134E+01_r8,0.241129E+01_r8, &
         &0.244157E+01_r8,0.247217E+01_r8,0.250316E+01_r8,0.253447E+01_r8,0.256617E+01_r8, &
         &0.259821E+01_r8,0.263064E+01_r8,0.266341E+01_r8,0.269661E+01_r8,0.273009E+01_r8, &
         &0.276403E+01_r8,0.279834E+01_r8,0.283302E+01_r8,0.286811E+01_r8,0.290358E+01_r8, &
         &0.293943E+01_r8,0.297571E+01_r8,0.301236E+01_r8,0.304946E+01_r8,0.308702E+01_r8, &
         &0.312491E+01_r8,0.316326E+01_r8,0.320208E+01_r8,0.324130E+01_r8,0.328092E+01_r8, &
         &0.332102E+01_r8,0.336162E+01_r8,0.340264E+01_r8,0.344407E+01_r8,0.348601E+01_r8, &
         &0.352838E+01_r8,0.357118E+01_r8,0.361449E+01_r8,0.365834E+01_r8,0.370264E+01_r8, &
         &0.374737E+01_r8,0.379265E+01_r8,0.383839E+01_r8,0.388469E+01_r8,0.393144E+01_r8, &
         &0.397876E+01_r8,0.402656E+01_r8,0.407492E+01_r8,0.412378E+01_r8,0.417313E+01_r8, &
         &0.422306E+01_r8,0.427359E+01_r8,0.432454E+01_r8,0.437617E+01_r8,0.442834E+01_r8, &
         &0.448102E+01_r8,0.453433E+01_r8,0.458816E+01_r8,0.464253E+01_r8,0.469764E+01_r8, &
         &0.475321E+01_r8,0.480942E+01_r8,0.486629E+01_r8,0.492372E+01_r8,0.498173E+01_r8/
    DATA (ES(IES),IES=381,475) /                                      &
         &0.504041E+01_r8,0.509967E+01_r8,0.515962E+01_r8,0.522029E+01_r8,0.528142E+01_r8, &
         &0.534337E+01_r8,0.540595E+01_r8,0.546912E+01_r8,0.553292E+01_r8,0.559757E+01_r8, &
         &0.566273E+01_r8,0.572864E+01_r8,0.579532E+01_r8,0.586266E+01_r8,0.593075E+01_r8, &
         &0.599952E+01_r8,0.606895E+01_r8,0.613918E+01_r8,0.621021E+01_r8,0.628191E+01_r8, &
         &0.635433E+01_r8,0.642755E+01_r8,0.650162E+01_r8,0.657639E+01_r8,0.665188E+01_r8, &
         &0.672823E+01_r8,0.680532E+01_r8,0.688329E+01_r8,0.696198E+01_r8,0.704157E+01_r8, &
         &0.712206E+01_r8,0.720319E+01_r8,0.728534E+01_r8,0.736829E+01_r8,0.745204E+01_r8, &
         &0.753671E+01_r8,0.762218E+01_r8,0.770860E+01_r8,0.779588E+01_r8,0.788408E+01_r8, &
         &0.797314E+01_r8,0.806318E+01_r8,0.815408E+01_r8,0.824599E+01_r8,0.833874E+01_r8, &
         &0.843254E+01_r8,0.852721E+01_r8,0.862293E+01_r8,0.871954E+01_r8,0.881724E+01_r8, &
         &0.891579E+01_r8,0.901547E+01_r8,0.911624E+01_r8,0.921778E+01_r8,0.932061E+01_r8, &
         &0.942438E+01_r8,0.952910E+01_r8,0.963497E+01_r8,0.974181E+01_r8,0.984982E+01_r8, &
         &0.995887E+01_r8,0.100690E+02_r8,0.101804E+02_r8,0.102926E+02_r8,0.104063E+02_r8, &
         &0.105210E+02_r8,0.106367E+02_r8,0.107536E+02_r8,0.108719E+02_r8,0.109912E+02_r8, &
         &0.111116E+02_r8,0.112333E+02_r8,0.113563E+02_r8,0.114804E+02_r8,0.116056E+02_r8, &
         &0.117325E+02_r8,0.118602E+02_r8,0.119892E+02_r8,0.121197E+02_r8,0.122513E+02_r8, &
         &0.123844E+02_r8,0.125186E+02_r8,0.126543E+02_r8,0.127912E+02_r8,0.129295E+02_r8, &
         &0.130691E+02_r8,0.132101E+02_r8,0.133527E+02_r8,0.134965E+02_r8,0.136415E+02_r8, &
         &0.137882E+02_r8,0.139361E+02_r8,0.140855E+02_r8,0.142366E+02_r8,0.143889E+02_r8/
    DATA (ES(IES),IES=476,570) /                                      &
         &0.145429E+02_r8,0.146982E+02_r8,0.148552E+02_r8,0.150135E+02_r8,0.151735E+02_r8, &
         &0.153349E+02_r8,0.154979E+02_r8,0.156624E+02_r8,0.158286E+02_r8,0.159965E+02_r8, &
         &0.161659E+02_r8,0.163367E+02_r8,0.165094E+02_r8,0.166838E+02_r8,0.168597E+02_r8, &
         &0.170375E+02_r8,0.172168E+02_r8,0.173979E+02_r8,0.175806E+02_r8,0.177651E+02_r8, &
         &0.179513E+02_r8,0.181394E+02_r8,0.183293E+02_r8,0.185210E+02_r8,0.187146E+02_r8, &
         &0.189098E+02_r8,0.191066E+02_r8,0.193059E+02_r8,0.195065E+02_r8,0.197095E+02_r8, &
         &0.199142E+02_r8,0.201206E+02_r8,0.203291E+02_r8,0.205397E+02_r8,0.207522E+02_r8, &
         &0.209664E+02_r8,0.211831E+02_r8,0.214013E+02_r8,0.216221E+02_r8,0.218448E+02_r8, &
         &0.220692E+02_r8,0.222959E+02_r8,0.225250E+02_r8,0.227559E+02_r8,0.229887E+02_r8, &
         &0.232239E+02_r8,0.234614E+02_r8,0.237014E+02_r8,0.239428E+02_r8,0.241872E+02_r8, &
         &0.244335E+02_r8,0.246824E+02_r8,0.249332E+02_r8,0.251860E+02_r8,0.254419E+02_r8, &
         &0.256993E+02_r8,0.259600E+02_r8,0.262225E+02_r8,0.264873E+02_r8,0.267552E+02_r8, &
         &0.270248E+02_r8,0.272970E+02_r8,0.275719E+02_r8,0.278497E+02_r8,0.281295E+02_r8, &
         &0.284117E+02_r8,0.286965E+02_r8,0.289843E+02_r8,0.292743E+02_r8,0.295671E+02_r8, &
         &0.298624E+02_r8,0.301605E+02_r8,0.304616E+02_r8,0.307650E+02_r8,0.310708E+02_r8, &
         &0.313803E+02_r8,0.316915E+02_r8,0.320064E+02_r8,0.323238E+02_r8,0.326437E+02_r8, &
         &0.329666E+02_r8,0.332928E+02_r8,0.336215E+02_r8,0.339534E+02_r8,0.342885E+02_r8, &
         &0.346263E+02_r8,0.349666E+02_r8,0.353109E+02_r8,0.356572E+02_r8,0.360076E+02_r8, &
         &0.363606E+02_r8,0.367164E+02_r8,0.370757E+02_r8,0.374383E+02_r8,0.378038E+02_r8/
    DATA (ES(IES),IES=571,665) /                                      &
         &0.381727E+02_r8,0.385453E+02_r8,0.389206E+02_r8,0.392989E+02_r8,0.396807E+02_r8, &
         &0.400663E+02_r8,0.404555E+02_r8,0.408478E+02_r8,0.412428E+02_r8,0.416417E+02_r8, &
         &0.420445E+02_r8,0.424502E+02_r8,0.428600E+02_r8,0.432733E+02_r8,0.436900E+02_r8, &
         &0.441106E+02_r8,0.445343E+02_r8,0.449620E+02_r8,0.453930E+02_r8,0.458280E+02_r8, &
         &0.462672E+02_r8,0.467096E+02_r8,0.471561E+02_r8,0.476070E+02_r8,0.480610E+02_r8, &
         &0.485186E+02_r8,0.489813E+02_r8,0.494474E+02_r8,0.499170E+02_r8,0.503909E+02_r8, &
         &0.508693E+02_r8,0.513511E+02_r8,0.518376E+02_r8,0.523277E+02_r8,0.528232E+02_r8, &
         &0.533213E+02_r8,0.538240E+02_r8,0.543315E+02_r8,0.548437E+02_r8,0.553596E+02_r8, &
         &0.558802E+02_r8,0.564046E+02_r8,0.569340E+02_r8,0.574672E+02_r8,0.580061E+02_r8, &
         &0.585481E+02_r8,0.590963E+02_r8,0.596482E+02_r8,0.602041E+02_r8,0.607649E+02_r8, &
         &0.613311E+02_r8,0.619025E+02_r8,0.624779E+02_r8,0.630574E+02_r8,0.636422E+02_r8, &
         &0.642324E+02_r8,0.648280E+02_r8,0.654278E+02_r8,0.660332E+02_r8,0.666426E+02_r8, &
         &0.672577E+02_r8,0.678771E+02_r8,0.685034E+02_r8,0.691328E+02_r8,0.697694E+02_r8, &
         &0.704103E+02_r8,0.710556E+02_r8,0.717081E+02_r8,0.723639E+02_r8,0.730269E+02_r8, &
         &0.736945E+02_r8,0.743681E+02_r8,0.750463E+02_r8,0.757309E+02_r8,0.764214E+02_r8, &
         &0.771167E+02_r8,0.778182E+02_r8,0.785246E+02_r8,0.792373E+02_r8,0.799564E+02_r8, &
         &0.806804E+02_r8,0.814109E+02_r8,0.821479E+02_r8,0.828898E+02_r8,0.836384E+02_r8, &
         &0.843922E+02_r8,0.851525E+02_r8,0.859198E+02_r8,0.866920E+02_r8,0.874712E+02_r8, &
         &0.882574E+02_r8,0.890486E+02_r8,0.898470E+02_r8,0.906525E+02_r8,0.914634E+02_r8/
    DATA (ES(IES),IES=666,760) /                                      &
         &0.922814E+02_r8,0.931048E+02_r8,0.939356E+02_r8,0.947736E+02_r8,0.956171E+02_r8, &
         &0.964681E+02_r8,0.973246E+02_r8,0.981907E+02_r8,0.990605E+02_r8,0.999399E+02_r8, &
         &0.100825E+03_r8,0.101718E+03_r8,0.102617E+03_r8,0.103523E+03_r8,0.104438E+03_r8, &
         &0.105358E+03_r8,0.106287E+03_r8,0.107221E+03_r8,0.108166E+03_r8,0.109115E+03_r8, &
         &0.110074E+03_r8,0.111039E+03_r8,0.112012E+03_r8,0.112992E+03_r8,0.113981E+03_r8, &
         &0.114978E+03_r8,0.115981E+03_r8,0.116993E+03_r8,0.118013E+03_r8,0.119041E+03_r8, &
         &0.120077E+03_r8,0.121122E+03_r8,0.122173E+03_r8,0.123234E+03_r8,0.124301E+03_r8, &
         &0.125377E+03_r8,0.126463E+03_r8,0.127556E+03_r8,0.128657E+03_r8,0.129769E+03_r8, &
         &0.130889E+03_r8,0.132017E+03_r8,0.133152E+03_r8,0.134299E+03_r8,0.135453E+03_r8, &
         &0.136614E+03_r8,0.137786E+03_r8,0.138967E+03_r8,0.140158E+03_r8,0.141356E+03_r8, &
         &0.142565E+03_r8,0.143781E+03_r8,0.145010E+03_r8,0.146247E+03_r8,0.147491E+03_r8, &
         &0.148746E+03_r8,0.150011E+03_r8,0.151284E+03_r8,0.152571E+03_r8,0.153862E+03_r8, &
         &0.155168E+03_r8,0.156481E+03_r8,0.157805E+03_r8,0.159137E+03_r8,0.160478E+03_r8, &
         &0.161832E+03_r8,0.163198E+03_r8,0.164569E+03_r8,0.165958E+03_r8,0.167348E+03_r8, &
         &0.168757E+03_r8,0.170174E+03_r8,0.171599E+03_r8,0.173037E+03_r8,0.174483E+03_r8, &
         &0.175944E+03_r8,0.177414E+03_r8,0.178892E+03_r8,0.180387E+03_r8,0.181886E+03_r8, &
         &0.183402E+03_r8,0.184930E+03_r8,0.186463E+03_r8,0.188012E+03_r8,0.189571E+03_r8, &
         &0.191146E+03_r8,0.192730E+03_r8,0.194320E+03_r8,0.195930E+03_r8,0.197546E+03_r8, &
         &0.199175E+03_r8,0.200821E+03_r8,0.202473E+03_r8,0.204142E+03_r8,0.205817E+03_r8/
    DATA (ES(IES),IES=761,855) /                                      &
         &0.207510E+03_r8,0.209216E+03_r8,0.210928E+03_r8,0.212658E+03_r8,0.214398E+03_r8, &
         &0.216152E+03_r8,0.217920E+03_r8,0.219698E+03_r8,0.221495E+03_r8,0.223297E+03_r8, &
         &0.225119E+03_r8,0.226951E+03_r8,0.228793E+03_r8,0.230654E+03_r8,0.232522E+03_r8, &
         &0.234413E+03_r8,0.236311E+03_r8,0.238223E+03_r8,0.240151E+03_r8,0.242090E+03_r8, &
         &0.244049E+03_r8,0.246019E+03_r8,0.248000E+03_r8,0.249996E+03_r8,0.252009E+03_r8, &
         &0.254037E+03_r8,0.256077E+03_r8,0.258128E+03_r8,0.260200E+03_r8,0.262284E+03_r8, &
         &0.264384E+03_r8,0.266500E+03_r8,0.268629E+03_r8,0.270779E+03_r8,0.272936E+03_r8, &
         &0.275110E+03_r8,0.277306E+03_r8,0.279509E+03_r8,0.281734E+03_r8,0.283966E+03_r8, &
         &0.286227E+03_r8,0.288494E+03_r8,0.290780E+03_r8,0.293083E+03_r8,0.295398E+03_r8, &
         &0.297737E+03_r8,0.300089E+03_r8,0.302453E+03_r8,0.304841E+03_r8,0.307237E+03_r8, &
         &0.309656E+03_r8,0.312095E+03_r8,0.314541E+03_r8,0.317012E+03_r8,0.319496E+03_r8, &
         &0.322005E+03_r8,0.324527E+03_r8,0.327063E+03_r8,0.329618E+03_r8,0.332193E+03_r8, &
         &0.334788E+03_r8,0.337396E+03_r8,0.340025E+03_r8,0.342673E+03_r8,0.345329E+03_r8, &
         &0.348019E+03_r8,0.350722E+03_r8,0.353440E+03_r8,0.356178E+03_r8,0.358938E+03_r8, &
         &0.361718E+03_r8,0.364513E+03_r8,0.367322E+03_r8,0.370160E+03_r8,0.373012E+03_r8, &
         &0.375885E+03_r8,0.378788E+03_r8,0.381691E+03_r8,0.384631E+03_r8,0.387579E+03_r8, &
         &0.390556E+03_r8,0.393556E+03_r8,0.396563E+03_r8,0.399601E+03_r8,0.402646E+03_r8, &
         &0.405730E+03_r8,0.408829E+03_r8,0.411944E+03_r8,0.415083E+03_r8,0.418236E+03_r8, &
         &0.421422E+03_r8,0.424632E+03_r8,0.427849E+03_r8,0.431099E+03_r8,0.434365E+03_r8/
    DATA (ES(IES),IES=856,950) /                                      &
         &0.437655E+03_r8,0.440970E+03_r8,0.444301E+03_r8,0.447666E+03_r8,0.451038E+03_r8, &
         &0.454445E+03_r8,0.457876E+03_r8,0.461316E+03_r8,0.464790E+03_r8,0.468281E+03_r8, &
         &0.471798E+03_r8,0.475342E+03_r8,0.478902E+03_r8,0.482497E+03_r8,0.486101E+03_r8, &
         &0.489741E+03_r8,0.493408E+03_r8,0.497083E+03_r8,0.500804E+03_r8,0.504524E+03_r8, &
         &0.508290E+03_r8,0.512074E+03_r8,0.515877E+03_r8,0.519717E+03_r8,0.523566E+03_r8, &
         &0.527462E+03_r8,0.531367E+03_r8,0.535301E+03_r8,0.539264E+03_r8,0.543245E+03_r8, &
         &0.547265E+03_r8,0.551305E+03_r8,0.555363E+03_r8,0.559462E+03_r8,0.563579E+03_r8, &
         &0.567727E+03_r8,0.571905E+03_r8,0.576102E+03_r8,0.580329E+03_r8,0.584576E+03_r8, &
         &0.588865E+03_r8,0.593185E+03_r8,0.597514E+03_r8,0.601885E+03_r8,0.606276E+03_r8, &
         &0.610699E+03_r8,0.615151E+03_r8,0.619625E+03_r8,0.624140E+03_r8,0.628671E+03_r8, &
         &0.633243E+03_r8,0.637845E+03_r8,0.642465E+03_r8,0.647126E+03_r8,0.651806E+03_r8, &
         &0.656527E+03_r8,0.661279E+03_r8,0.666049E+03_r8,0.670861E+03_r8,0.675692E+03_r8, &
         &0.680566E+03_r8,0.685471E+03_r8,0.690396E+03_r8,0.695363E+03_r8,0.700350E+03_r8, &
         &0.705381E+03_r8,0.710444E+03_r8,0.715527E+03_r8,0.720654E+03_r8,0.725801E+03_r8, &
         &0.730994E+03_r8,0.736219E+03_r8,0.741465E+03_r8,0.746756E+03_r8,0.752068E+03_r8, &
         &0.757426E+03_r8,0.762819E+03_r8,0.768231E+03_r8,0.773692E+03_r8,0.779172E+03_r8, &
         &0.784701E+03_r8,0.790265E+03_r8,0.795849E+03_r8,0.801483E+03_r8,0.807137E+03_r8, &
         &0.812842E+03_r8,0.818582E+03_r8,0.824343E+03_r8,0.830153E+03_r8,0.835987E+03_r8, &
         &0.841871E+03_r8,0.847791E+03_r8,0.853733E+03_r8,0.859727E+03_r8,0.865743E+03_r8/
    DATA (ES(IES),IES=951,1045) /                                     &
         &0.871812E+03_r8,0.877918E+03_r8,0.884046E+03_r8,0.890228E+03_r8,0.896433E+03_r8, &
         &0.902690E+03_r8,0.908987E+03_r8,0.915307E+03_r8,0.921681E+03_r8,0.928078E+03_r8, &
         &0.934531E+03_r8,0.941023E+03_r8,0.947539E+03_r8,0.954112E+03_r8,0.960708E+03_r8, &
         &0.967361E+03_r8,0.974053E+03_r8,0.980771E+03_r8,0.987545E+03_r8,0.994345E+03_r8, &
         &0.100120E+04_r8,0.100810E+04_r8,0.101502E+04_r8,0.102201E+04_r8,0.102902E+04_r8, &
         &0.103608E+04_r8,0.104320E+04_r8,0.105033E+04_r8,0.105753E+04_r8,0.106475E+04_r8, &
         &0.107204E+04_r8,0.107936E+04_r8,0.108672E+04_r8,0.109414E+04_r8,0.110158E+04_r8, &
         &0.110908E+04_r8,0.111663E+04_r8,0.112421E+04_r8,0.113185E+04_r8,0.113952E+04_r8, &
         &0.114725E+04_r8,0.115503E+04_r8,0.116284E+04_r8,0.117071E+04_r8,0.117861E+04_r8, &
         &0.118658E+04_r8,0.119459E+04_r8,0.120264E+04_r8,0.121074E+04_r8,0.121888E+04_r8, &
         &0.122709E+04_r8,0.123534E+04_r8,0.124362E+04_r8,0.125198E+04_r8,0.126036E+04_r8, &
         &0.126881E+04_r8,0.127731E+04_r8,0.128584E+04_r8,0.129444E+04_r8,0.130307E+04_r8, &
         &0.131177E+04_r8,0.132053E+04_r8,0.132931E+04_r8,0.133817E+04_r8,0.134705E+04_r8, &
         &0.135602E+04_r8,0.136503E+04_r8,0.137407E+04_r8,0.138319E+04_r8,0.139234E+04_r8, &
         &0.140156E+04_r8,0.141084E+04_r8,0.142015E+04_r8,0.142954E+04_r8,0.143896E+04_r8, &
         &0.144845E+04_r8,0.145800E+04_r8,0.146759E+04_r8,0.147725E+04_r8,0.148694E+04_r8, &
         &0.149672E+04_r8,0.150655E+04_r8,0.151641E+04_r8,0.152635E+04_r8,0.153633E+04_r8, &
         &0.154639E+04_r8,0.155650E+04_r8,0.156665E+04_r8,0.157688E+04_r8,0.158715E+04_r8, &
         &0.159750E+04_r8,0.160791E+04_r8,0.161836E+04_r8,0.162888E+04_r8,0.163945E+04_r8/
    DATA (ES(IES),IES=1046,1140) /                                    &
         &0.165010E+04_r8,0.166081E+04_r8,0.167155E+04_r8,0.168238E+04_r8,0.169325E+04_r8, &
         &0.170420E+04_r8,0.171522E+04_r8,0.172627E+04_r8,0.173741E+04_r8,0.174859E+04_r8, &
         &0.175986E+04_r8,0.177119E+04_r8,0.178256E+04_r8,0.179402E+04_r8,0.180552E+04_r8, &
         &0.181711E+04_r8,0.182877E+04_r8,0.184046E+04_r8,0.185224E+04_r8,0.186407E+04_r8, &
         &0.187599E+04_r8,0.188797E+04_r8,0.190000E+04_r8,0.191212E+04_r8,0.192428E+04_r8, &
         &0.193653E+04_r8,0.194886E+04_r8,0.196122E+04_r8,0.197368E+04_r8,0.198618E+04_r8, &
         &0.199878E+04_r8,0.201145E+04_r8,0.202416E+04_r8,0.203698E+04_r8,0.204983E+04_r8, &
         &0.206278E+04_r8,0.207580E+04_r8,0.208887E+04_r8,0.210204E+04_r8,0.211525E+04_r8, &
         &0.212856E+04_r8,0.214195E+04_r8,0.215538E+04_r8,0.216892E+04_r8,0.218249E+04_r8, &
         &0.219618E+04_r8,0.220994E+04_r8,0.222375E+04_r8,0.223766E+04_r8,0.225161E+04_r8, &
         &0.226567E+04_r8,0.227981E+04_r8,0.229399E+04_r8,0.230829E+04_r8,0.232263E+04_r8, &
         &0.233708E+04_r8,0.235161E+04_r8,0.236618E+04_r8,0.238087E+04_r8,0.239560E+04_r8, &
         &0.241044E+04_r8,0.242538E+04_r8,0.244035E+04_r8,0.245544E+04_r8,0.247057E+04_r8, &
         &0.248583E+04_r8,0.250116E+04_r8,0.251654E+04_r8,0.253204E+04_r8,0.254759E+04_r8, &
         &0.256325E+04_r8,0.257901E+04_r8,0.259480E+04_r8,0.261073E+04_r8,0.262670E+04_r8, &
         &0.264279E+04_r8,0.265896E+04_r8,0.267519E+04_r8,0.269154E+04_r8,0.270794E+04_r8, &
         &0.272447E+04_r8,0.274108E+04_r8,0.275774E+04_r8,0.277453E+04_r8,0.279137E+04_r8, &
         &0.280834E+04_r8,0.282540E+04_r8,0.284251E+04_r8,0.285975E+04_r8,0.287704E+04_r8, &
         &0.289446E+04_r8,0.291198E+04_r8,0.292954E+04_r8,0.294725E+04_r8,0.296499E+04_r8/
    DATA (ES(IES),IES=1141,1235) /                                    &
         &0.298288E+04_r8,0.300087E+04_r8,0.301890E+04_r8,0.303707E+04_r8,0.305529E+04_r8, &
         &0.307365E+04_r8,0.309211E+04_r8,0.311062E+04_r8,0.312927E+04_r8,0.314798E+04_r8, &
         &0.316682E+04_r8,0.318577E+04_r8,0.320477E+04_r8,0.322391E+04_r8,0.324310E+04_r8, &
         &0.326245E+04_r8,0.328189E+04_r8,0.330138E+04_r8,0.332103E+04_r8,0.334073E+04_r8, &
         &0.336058E+04_r8,0.338053E+04_r8,0.340054E+04_r8,0.342069E+04_r8,0.344090E+04_r8, &
         &0.346127E+04_r8,0.348174E+04_r8,0.350227E+04_r8,0.352295E+04_r8,0.354369E+04_r8, &
         &0.356458E+04_r8,0.358559E+04_r8,0.360664E+04_r8,0.362787E+04_r8,0.364914E+04_r8, &
         &0.367058E+04_r8,0.369212E+04_r8,0.371373E+04_r8,0.373548E+04_r8,0.375731E+04_r8, &
         &0.377929E+04_r8,0.380139E+04_r8,0.382355E+04_r8,0.384588E+04_r8,0.386826E+04_r8, &
         &0.389081E+04_r8,0.391348E+04_r8,0.393620E+04_r8,0.395910E+04_r8,0.398205E+04_r8, &
         &0.400518E+04_r8,0.402843E+04_r8,0.405173E+04_r8,0.407520E+04_r8,0.409875E+04_r8, &
         &0.412246E+04_r8,0.414630E+04_r8,0.417019E+04_r8,0.419427E+04_r8,0.421840E+04_r8, &
         &0.424272E+04_r8,0.426715E+04_r8,0.429165E+04_r8,0.431634E+04_r8,0.434108E+04_r8, &
         &0.436602E+04_r8,0.439107E+04_r8,0.441618E+04_r8,0.444149E+04_r8,0.446685E+04_r8, &
         &0.449241E+04_r8,0.451810E+04_r8,0.454385E+04_r8,0.456977E+04_r8,0.459578E+04_r8, &
         &0.462197E+04_r8,0.464830E+04_r8,0.467468E+04_r8,0.470127E+04_r8,0.472792E+04_r8, &
         &0.475477E+04_r8,0.478175E+04_r8,0.480880E+04_r8,0.483605E+04_r8,0.486336E+04_r8, &
         &0.489087E+04_r8,0.491853E+04_r8,0.494623E+04_r8,0.497415E+04_r8,0.500215E+04_r8, &
         &0.503034E+04_r8,0.505867E+04_r8,0.508707E+04_r8,0.511568E+04_r8,0.514436E+04_r8/
    DATA (ES(IES),IES=1236,1330) /                                    &
         &0.517325E+04_r8,0.520227E+04_r8,0.523137E+04_r8,0.526068E+04_r8,0.529005E+04_r8, &
         &0.531965E+04_r8,0.534939E+04_r8,0.537921E+04_r8,0.540923E+04_r8,0.543932E+04_r8, &
         &0.546965E+04_r8,0.550011E+04_r8,0.553064E+04_r8,0.556139E+04_r8,0.559223E+04_r8, &
         &0.562329E+04_r8,0.565449E+04_r8,0.568577E+04_r8,0.571727E+04_r8,0.574884E+04_r8, &
         &0.578064E+04_r8,0.581261E+04_r8,0.584464E+04_r8,0.587692E+04_r8,0.590924E+04_r8, &
         &0.594182E+04_r8,0.597455E+04_r8,0.600736E+04_r8,0.604039E+04_r8,0.607350E+04_r8, &
         &0.610685E+04_r8,0.614036E+04_r8,0.617394E+04_r8,0.620777E+04_r8,0.624169E+04_r8, &
         &0.627584E+04_r8,0.631014E+04_r8,0.634454E+04_r8,0.637918E+04_r8,0.641390E+04_r8, &
         &0.644887E+04_r8,0.648400E+04_r8,0.651919E+04_r8,0.655467E+04_r8,0.659021E+04_r8, &
         &0.662599E+04_r8,0.666197E+04_r8,0.669800E+04_r8,0.673429E+04_r8,0.677069E+04_r8, &
         &0.680735E+04_r8,0.684415E+04_r8,0.688104E+04_r8,0.691819E+04_r8,0.695543E+04_r8, &
         &0.699292E+04_r8,0.703061E+04_r8,0.706837E+04_r8,0.710639E+04_r8,0.714451E+04_r8, &
         &0.718289E+04_r8,0.722143E+04_r8,0.726009E+04_r8,0.729903E+04_r8,0.733802E+04_r8, &
         &0.737729E+04_r8,0.741676E+04_r8,0.745631E+04_r8,0.749612E+04_r8,0.753602E+04_r8, &
         &0.757622E+04_r8,0.761659E+04_r8,0.765705E+04_r8,0.769780E+04_r8,0.773863E+04_r8, &
         &0.777975E+04_r8,0.782106E+04_r8,0.786246E+04_r8,0.790412E+04_r8,0.794593E+04_r8, &
         &0.798802E+04_r8,0.803028E+04_r8,0.807259E+04_r8,0.811525E+04_r8,0.815798E+04_r8, &
         &0.820102E+04_r8,0.824427E+04_r8,0.828757E+04_r8,0.833120E+04_r8,0.837493E+04_r8, &
         &0.841895E+04_r8,0.846313E+04_r8,0.850744E+04_r8,0.855208E+04_r8,0.859678E+04_r8/
    DATA (ES(IES),IES=1331,1425) /                                    &
         &0.864179E+04_r8,0.868705E+04_r8,0.873237E+04_r8,0.877800E+04_r8,0.882374E+04_r8, &
         &0.886979E+04_r8,0.891603E+04_r8,0.896237E+04_r8,0.900904E+04_r8,0.905579E+04_r8, &
         &0.910288E+04_r8,0.915018E+04_r8,0.919758E+04_r8,0.924529E+04_r8,0.929310E+04_r8, &
         &0.934122E+04_r8,0.938959E+04_r8,0.943804E+04_r8,0.948687E+04_r8,0.953575E+04_r8, &
         &0.958494E+04_r8,0.963442E+04_r8,0.968395E+04_r8,0.973384E+04_r8,0.978383E+04_r8, &
         &0.983412E+04_r8,0.988468E+04_r8,0.993534E+04_r8,0.998630E+04_r8,0.100374E+05_r8, &
         &0.100888E+05_r8,0.101406E+05_r8,0.101923E+05_r8,0.102444E+05_r8,0.102966E+05_r8, &
         &0.103492E+05_r8,0.104020E+05_r8,0.104550E+05_r8,0.105082E+05_r8,0.105616E+05_r8, &
         &0.106153E+05_r8,0.106693E+05_r8,0.107234E+05_r8,0.107779E+05_r8,0.108325E+05_r8, &
         &0.108874E+05_r8,0.109425E+05_r8,0.109978E+05_r8,0.110535E+05_r8,0.111092E+05_r8, &
         &0.111653E+05_r8,0.112217E+05_r8,0.112782E+05_r8,0.113350E+05_r8,0.113920E+05_r8, &
         &0.114493E+05_r8,0.115070E+05_r8,0.115646E+05_r8,0.116228E+05_r8,0.116809E+05_r8, &
         &0.117396E+05_r8,0.117984E+05_r8,0.118574E+05_r8,0.119167E+05_r8,0.119762E+05_r8, &
         &0.120360E+05_r8,0.120962E+05_r8,0.121564E+05_r8,0.122170E+05_r8,0.122778E+05_r8, &
         &0.123389E+05_r8,0.124004E+05_r8,0.124619E+05_r8,0.125238E+05_r8,0.125859E+05_r8, &
         &0.126484E+05_r8,0.127111E+05_r8,0.127739E+05_r8,0.128372E+05_r8,0.129006E+05_r8, &
         &0.129644E+05_r8,0.130285E+05_r8,0.130927E+05_r8,0.131573E+05_r8,0.132220E+05_r8, &
         &0.132872E+05_r8,0.133526E+05_r8,0.134182E+05_r8,0.134842E+05_r8,0.135503E+05_r8, &
         &0.136168E+05_r8,0.136836E+05_r8,0.137505E+05_r8,0.138180E+05_r8,0.138854E+05_r8/
    DATA (ES(IES),IES=1426,1520) /                                    &
         &0.139534E+05_r8,0.140216E+05_r8,0.140900E+05_r8,0.141588E+05_r8,0.142277E+05_r8, &
         &0.142971E+05_r8,0.143668E+05_r8,0.144366E+05_r8,0.145069E+05_r8,0.145773E+05_r8, &
         &0.146481E+05_r8,0.147192E+05_r8,0.147905E+05_r8,0.148622E+05_r8,0.149341E+05_r8, &
         &0.150064E+05_r8,0.150790E+05_r8,0.151517E+05_r8,0.152250E+05_r8,0.152983E+05_r8, &
         &0.153721E+05_r8,0.154462E+05_r8,0.155205E+05_r8,0.155952E+05_r8,0.156701E+05_r8, &
         &0.157454E+05_r8,0.158211E+05_r8,0.158969E+05_r8,0.159732E+05_r8,0.160496E+05_r8, &
         &0.161265E+05_r8,0.162037E+05_r8,0.162811E+05_r8,0.163589E+05_r8,0.164369E+05_r8, &
         &0.165154E+05_r8,0.165942E+05_r8,0.166732E+05_r8,0.167526E+05_r8,0.168322E+05_r8, &
         &0.169123E+05_r8,0.169927E+05_r8,0.170733E+05_r8,0.171543E+05_r8,0.172356E+05_r8, &
         &0.173173E+05_r8,0.173993E+05_r8,0.174815E+05_r8,0.175643E+05_r8,0.176471E+05_r8, &
         &0.177305E+05_r8,0.178143E+05_r8,0.178981E+05_r8,0.179826E+05_r8,0.180671E+05_r8, &
         &0.181522E+05_r8,0.182377E+05_r8,0.183232E+05_r8,0.184093E+05_r8,0.184955E+05_r8, &
         &0.185823E+05_r8,0.186695E+05_r8,0.187568E+05_r8,0.188447E+05_r8,0.189326E+05_r8, &
         &0.190212E+05_r8,0.191101E+05_r8,0.191991E+05_r8,0.192887E+05_r8,0.193785E+05_r8, &
         &0.194688E+05_r8,0.195595E+05_r8,0.196503E+05_r8,0.197417E+05_r8,0.198332E+05_r8, &
         &0.199253E+05_r8,0.200178E+05_r8,0.201105E+05_r8,0.202036E+05_r8,0.202971E+05_r8, &
         &0.203910E+05_r8,0.204853E+05_r8,0.205798E+05_r8,0.206749E+05_r8,0.207701E+05_r8, &
         &0.208659E+05_r8,0.209621E+05_r8,0.210584E+05_r8,0.211554E+05_r8,0.212524E+05_r8, &
         &0.213501E+05_r8,0.214482E+05_r8,0.215465E+05_r8,0.216452E+05_r8,0.217442E+05_r8/
    DATA (ES(IES),IES=1521,1552) /                                    &
         &0.218439E+05_r8,0.219439E+05_r8,0.220440E+05_r8,0.221449E+05_r8,0.222457E+05_r8, &
         &0.223473E+05_r8,0.224494E+05_r8,0.225514E+05_r8,0.226542E+05_r8,0.227571E+05_r8, &
         &0.228606E+05_r8,0.229646E+05_r8,0.230687E+05_r8,0.231734E+05_r8,0.232783E+05_r8, &
         &0.233839E+05_r8,0.234898E+05_r8,0.235960E+05_r8,0.237027E+05_r8,0.238097E+05_r8, &
         &0.239173E+05_r8,0.240254E+05_r8,0.241335E+05_r8,0.242424E+05_r8,0.243514E+05_r8, &
         &0.244611E+05_r8,0.245712E+05_r8,0.246814E+05_r8,0.247923E+05_r8,0.249034E+05_r8, &
         &0.250152E+05_r8,0.250152E+05_r8/
    !
    QS=0.0_r8
    ! Npnts_do_1:
    DO I=1,NPNTS
       !      COMPUTE THE FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE IN A
       !      PURE WATER SYSTEM TO SAT VAPOUR PRESSURE IN AIR, FSUBW.
       !      THIS FORMULA IS TAKEN FROM EQUATION A4.7 OF ADRIAN GILL'S BOOK:
       !      ATMOSPHERE-ocean DYNAMICS. NOTE THAT HIS FORMULA WORKS IN TERMS
       !      OF PRESSURE IN MB AND TEMPERATURE IN CELSIUS, SO CONVERSION OF
       !      UNITS LEADS TO THE SLIGHTLY DIFFERENT EQUATION USED HERE.
       !
       FSUBW = 1.0_r8 + 1.0E-8_r8*P(I)*( 4.5_r8 +                               &
            &    6.0E-4_r8*( T(I) - ZERODEGC )*( T(I) - ZERODEGC ) )
       !
       !      USE THE LOOKUP TABLE TO FIND SATURATED VAPOUR PRESSURE, AND STORE
       !      IT IN QS.
       !
       TT = MAX(T_LOW,T(I))
       TT = MIN(T_HIGH,TT)
       !
       ATABLE = (TT - T_LOW + DELTA_T) / DELTA_T
       ITABLE = INT(ATABLE)
       ATABLE = ATABLE - ITABLE
       !
       QS(I) = (1.0_r8 - ATABLE)*ES(ITABLE) + ATABLE*ES(ITABLE+1)
       !
       !      MULTIPLY BY FSUBW TO CONVERT TO SATURATED VAPOUR PRESSURE IN AIR
       !      (EQUATION A4.6 OF ADRIAN GILL'S BOOK).
       !
       QS(I) = QS(I) * FSUBW
       !
       !      NOW FORM THE ACCURATE EXPRESSION FOR QS, WHICH IS A REARRANGED
       !      VERSION OF EQUATION A4.3 OF GILL'S BOOK.
       !
       !      NOTE THAT AT VERY LOW PRESSURES WE APPLY A FIX, TO PREVENT A
       !      SINGULARITY (Qsat tends to 1. kg/kg).
       !
       QS(I) = ( EPSILON*QS(I) ) /                                     &
            &          ( MAX(P(I),QS(I)) - ONE_MINUS_EPSILON*QS(I) )
       !
    END DO ! Npnts_do_1
    !
    RETURN
  END SUBROUTINE QSAT_WAT
  ! ======================================================================
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


END MODULE Micro_UKME

!PROGRAM Main
!  Use Micro_UKME
!END PROGRAM Main
