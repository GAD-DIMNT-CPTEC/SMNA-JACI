MODULE CloudFraction_UKME
  IMPLICIT NONE
SAVE

  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers

  !ftn -g  -fbounds-check -Waliasing -Wall -fbacktrace -ffpe-trap=invalid,overflow,zero  -finit-real=nan -finit-integer=nan -fconvert=big-endian  -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2 

  !
  !
  !  qt_bal_cld_____ls_arcld______ls_cld
  !                            |
  !                            |__qsat_wat_mix
  !                            |
  !                            |__qsat_wat_mix
  !                            |
  !                            |__qsat_wat_mix
  !                            |
  !                            |__ls_cld
  !                            |
  !                            |__ls_cld
  !                            |
  !                            |__LS_ACF_Brooks
  !
  !
  !*----------------------------------------------------------------------
  !*L------------------COMDECK C_O_DG_C-----------------------------------
  ! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
  ! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
  ! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

  REAL(KIND=r8), PARAMETER :: ZeroDegC = 273.15_r8
  REAL(KIND=r8), PARAMETER :: TFS      = 271.35_r8
  REAL(KIND=r8), PARAMETER :: TM       = 273.15_r8
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

  REAL(KIND=r8), PARAMETER :: Epsilonk  = 0.62198_r8

  REAL(KIND=r8), PARAMETER :: Epsilon   = 0.62198_r8
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

  INTEGER :: IESW     ! loop counter for data statement look-up table

  REAL(KIND=r8) :: ESW2(0:N+1)
  !                         TABLE OF SATURATION WATER VAPOUR PRESSURE (PA)
  !                       - SET BY DATA STATEMENT CALCULATED FROM THE
  !                         GOFF-GRATCH FORMULAE AS TAKEN FROM LANDOLT-
  !                         BORNSTEIN, 1987 NUMERICAL DATA AND FUNCTIONAL
  !                         RELATIONSHIPS IN SCIENCE AND TECHNOLOGY.
  !                         GROUP V/ VOL 4B METEOROLOGY. PHYSICAL AND
  !                         CHEMICAL PROPERTIES OF AIR, P35
  ! Large scale cloud:

  LOGICAL :: L_eacf                ! Use empirically adjusted
  ! cloud fraction

  INTEGER :: cloud_fraction_method ! Selects total cloud fraction
  ! calculation method
  INTEGER :: ice_fraction_method   ! Selects ice cloud fraction
  ! calculation method

  REAL(KIND=r8)    :: overlap_ice_liquid    ! Generic overlap parameter
  ! between ice and liquid phases
  REAL(KIND=r8)    :: ctt_weight            ! Cloud top temperature weight
  REAL(KIND=r8)    :: t_weight              ! Local temperature weight
  REAL(KIND=r8)    :: qsat_fixed            ! Prescribed qsat value
  REAL(KIND=r8)    :: sub_cld               ! Scaling factor
  !  REAL(KIND=r8)    :: dbsdtbs_turb_0        ! PC2 erosion rate / s-1

  !      NAMELIST/RUN_Cloud/L_eacf,cloud_fraction_method,                  &
  !     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
  !     &  qsat_fixed,sub_cld,dbsdtbs_turb_0
  INTEGER :: iMax

 ! REAL(KIND=r8), PARAMETER :: RHCRIT0(1:38) = (/&
 !      0.950_r8,0.925_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.850_r8,0.850_r8,0.850_r8,0.800_r8, &
 !      0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8, &
 !      0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8, &
 !      0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8,0.800_r8/)
!
 REAL(KIND=r8), PARAMETER :: RHCRIT0(1:38) = (/&
       0.990_r8,0.955_r8,0.940_r8,0.940_r8,0.940_r8,0.940_r8,0.920_r8,0.920_r8,0.920_r8,0.900_r8, &
       0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8, &
       0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8, &
       0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8,0.900_r8/)


  PUBLIC :: Init_CloudFraction_UKME
  PUBLIC :: Run_CloudFraction_UKME

CONTAINS

  SUBROUTINE Init_CloudFraction_UKME(iMax_in)
    IMPLICIT NONE
    INTEGER , INTENT(IN   )  ::  iMax_in
    iMax=iMax_in
    CALL qsat_wat_data()

  END SUBROUTINE Init_CloudFraction_UKME

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

  SUBROUTINE Run_CloudFraction_UKME(nCols,LEVELS,si,sl,pblh,colrad,topo,ps,kuo,T3,q ,qcl ,qcf ,qcf2,CF ,CFL,CFF)
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: LEVELS
    REAL(KIND=r8), INTENT(IN   ) :: si      (LEVELS+1)   !
    REAL(KIND=r8), INTENT(IN   ) :: sl      (LEVELS)   !
    REAL(KIND=r8), INTENT(IN   ) :: pblh    (1:nCols)   !
    REAL(KIND=r8), INTENT(IN   ) :: colrad  (1:nCols)   !
    REAL(KIND=r8), INTENT(IN   ) :: topo    (1:nCols)   !
    REAL(KIND=r8), INTENT(IN   ) :: ps      (1:nCols)   ! pressure at all points, on theta levels (Pa).
    INTEGER      , INTENT(IN   ) :: kuo     (1:nCols) 
    REAL(KIND=r8), INTENT(INOUT) :: T3      (1:nCols,LEVELS)
    REAL(KIND=r8), INTENT(INOUT) :: q       (1:nCols,LEVELS)
    REAL(KIND=r8), INTENT(INOUT) :: qcl     (1:nCols,LEVELS)
    REAL(KIND=r8), INTENT(IN   ) :: qcf     (1:nCols,LEVELS)
    REAL(KIND=r8), INTENT(IN   ) :: qcf2    (1:nCols,LEVELS)
    REAL(KIND=r8), INTENT(OUT  ) :: CF      (nCols,LEVELS) !Cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8), INTENT(OUT  ) :: CFL     (nCols,LEVELS) !Liquid cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8), INTENT(OUT  ) :: CFF     (nCols,LEVELS) !Frozen cloud fraction at processed levels (decimal fraction).

    INTEGER       :: ntml          (nCols)       
    LOGICAL       :: cumulus       (nCols)           ! bl convection flag
    REAL(KIND=r8) :: p_theta_levels(nCols,levels)    ! pressure at all points (Pa).
    REAL(KIND=r8) :: p             (nCols,LEVELS+1)  ! pressure at all points, on u,v levels (Pa).
    REAL(KIND=r8) :: rhcrit        (nCols,LEVELS)    ! Critical relative humidity.  See the the paragraph incorporating
                      ! eqs P292.11 to P292.14; the values need to be tuned for a given
    ! set of levels.
    INTEGER       :: bl_levels                       ! No. of boundary layer levels
    REAL(KIND=r8) :: delta_lambda
    REAL(KIND=r8) :: delta_phi
    REAL(KIND=r8) :: FV_cos_theta_latitude(nCols)     ! Finite volume cos(lat)
    LOGICAL       :: L_CLD_AREA    ! true if using area cloud fraction (ACF)
    LOGICAL       :: L_ACF_Cusack  ! ... and selected Cusack and PC2 off
    LOGICAL       :: L_ACF_Brooks  ! ... and selected Brooks
    LOGICAL       :: L_eacf        ! true if using empirically adjusted cloud fraction
    LOGICAL       :: L_conv4a      ! true if using 4A convection scheme
    LOGICAL       :: L_mcr_qcf2    ! true if second cloud ice variable in use
    LOGICAL       :: l_mixing_ratio! true if using mixing ratio formulation

    REAL(KIND=r8) :: area_cloud_fraction  (1:nCols, LEVELS)
    REAL(KIND=r8) :: bulk_cloud_fraction  (1:nCols,LEVELS)
    REAL(KIND=r8) :: cloud_fraction_liquid(1:nCols, LEVELS)
    REAL(KIND=r8) :: cloud_fraction_frozen(1:nCols, LEVELS)
    REAL(KIND=r8) :: r_theta_levels     (1:nCols,0:LEVELS) ! height of theta level above earth's centre (m)
    REAL(KIND=r8) :: zi   (1:nCols,1:LEVELS+1)! Height above surface at interfaces 
    REAL(KIND=r8) :: zm   (1:nCols,1:LEVELS)  ! Geopotential height at mid level

    REAL(KIND=r8) :: flip_pint  (nCols,LEVELS+1)   ! Interface pressures  
    REAL(KIND=r8) :: flip_pmid  (nCols,LEVELS)     ! Midpoint pressures 
    REAL(KIND=r8) :: flip_t     (nCols,LEVELS)     ! temperature
    REAL(KIND=r8) :: flip_q     (nCols,LEVELS)     ! specific humidity
    REAL(KIND=r8) :: flip_pdel  (nCols,LEVELS)     ! layer thickness
    REAL(KIND=r8) :: flip_rpdel (nCols,LEVELS)     ! inverse of layer thickness
    REAL(KIND=r8) :: flip_lnpmid(nCols,LEVELS)     ! Log Midpoint pressures    
    REAL(KIND=r8) :: flip_lnpint(nCols,LEVELS+1)   ! Log interface pressures
    REAL(KIND=r8) :: flip_zi   (1:nCols,1:LEVELS+1)! Height above surface at interfaces 
    REAL(KIND=r8) :: flip_zm   (1:nCols,1:LEVELS)  ! Geopotential height at mid level
    INTEGER :: kflip
    INTEGER :: i
    INTEGER :: k

    !*----------------------------------------------------------------------
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8
    ! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
    REAL(KIND=r8), PARAMETER  :: CP     = 1005.0_r8



    DO i=1,nCols
       flip_pint       (i,LEVELS+1) = ps(i)*si(1) ! gps --> Pa
    END DO
    DO k=LEVELS,1,-1
       kflip=LEVELS+2-k
       DO i=1,nCols
          flip_pint    (i,k)      = MAX(si(kflip)*ps(i) ,1.0e-12_r8)
       END DO
    END DO
    DO k=1,LEVELS
       kflip=LEVELS+1-k
       DO i=1,nCols
          flip_t   (i,kflip) =  T3 (i,k)
          flip_q   (i,kflip) =  Q (i,k)
          flip_pmid(i,kflip) =  sl(  k)*ps (i)
       END DO
    END DO
    DO k=1,LEVELS
       DO i=1,nCols    
          flip_pdel    (i,k) = MAX(flip_pint(i,k+1) - flip_pint(i,k),1.0e-12_r8)
          flip_rpdel   (i,k) = 1.0_r8/MAX((flip_pint(i,k+1) - flip_pint(i,k)),1.0e-12_r8)
          flip_lnpmid  (i,k) = LOG(flip_pmid(i,k))
       END DO
    END DO
    DO k=1,LEVELS+1
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
         flip_lnpint(1:nCols,1:LEVELS+1)   , flip_pint  (1:nCols,1:LEVELS+1) , &
         flip_pmid  (1:nCols,1:LEVELS)     , flip_pdel  (1:nCols,1:LEVELS)   , flip_rpdel(1:nCols,1:LEVELS)   , &
         flip_t     (1:nCols,1:LEVELS)     , flip_q     (1:nCols,1:LEVELS)   , rair   , gravit , zvir   ,&
         flip_zi    (1:nCols,1:LEVELS+1)   , flip_zm    (1:nCols,1:LEVELS)   , nCols   ,nCols, LEVELS)

    DO i=1,nCols
       zi(i,1) = flip_zi      (i,LEVELS+1)
       p (i,1) = flip_pint    (i,LEVELS+1)

    END DO
    DO k=1,LEVELS
       kflip=LEVELS+1-k
       DO i=1,nCols
          zi (i,k+1) = flip_zi    (i,kflip)
          zm (i,k  ) = flip_zm    (i,kflip)
          p  (i,k+1) = flip_pint  (i,kflip)
          p_theta_levels(i,k) = flip_pmid(i,kflip)
       END DO
    END DO
    DO i=1,nCols
       R_THETA_LEVELS(i,0)      = Earth_Radius + MAX(topo(i),0.0_r8) + zi(i,1)
    END DO
    DO k=1,LEVELS
       DO i=1,nCols
          R_THETA_LEVELS(i,k)   = (Earth_Radius + MAX(topo(i),0.0_r8)) + zi(i,k+1)
       END DO
    END DO

    DO k=1,LEVELS
       DO i=1,nCols
          RHCRIT       (i,k)=RHCRIT0(idx(LEVELS,38,k))
       END DO
    END DO
    ntml=1
    DO k=1,LEVELS
       DO i=1,nCols
          IF(zm(i,k) <= pblh(i))THEN
             ntml(i)=k
          END IF
       END DO
    END DO

    bl_levels=LEVELS
    !cloud_fraction_method  =  1!  Use minimum overlap condition
    cloud_fraction_method  =  2! Calculate possible overlaps between ice and liquid in THIS layer
    ! =\ 1 and 2 No total cloud fraction method defined
    !ice_fraction_method    =  1
    ice_fraction_method    =  2! Use cloud top temperature and a fixed qsat to give QCFRBS
    ! No ice cloud fraction method defined
    overlap_ice_liquid     =  0.5! Generic overlap parameter between ice and liquid phases
    !IF (ICE_FRACTION_METHOD  ==  2) THEN
    CTT_WEIGHT  =0.20_r8 ! Cloud top temperature weight
    T_WEIGHT    =0.60_r8 ! Local temperature weight
    qsat_fixed  =0.062198_r8     ! Prescribed qsat value Saturation mixing ratio or saturation specific humidity at temperature(kg/kg).
    sub_cld     =0.06_r8     ! Scaling factor
    !SUBGRID = SUB_CLD ** (1.0_r8-T_WEIGHT)/ QSAT_FIXED ** (1.0_r8-T_WEIGHT-CTT_WEIGHT) ice_fraction_method eq 2
    !          CALL qsat_wat_mix(QSL_CTT,CTT,p_theta_levels(1,1,k),          &
    !                        nCols*rows,l_mixing_ratio)
    DO i = 1, nCols
       ! colrad.....colatitude  colrad=0 - 3.14 (0-180)from np to sp in radians
       !IF((((colrad(i)*180.0_r8)/3.1415926e0_r8)-90.0_r8)  > 0.0_r8 ) THEN
       FV_cos_theta_latitude(i)   = COS(((colrad(i)))-(3.1415926e0_r8/2.0_r8))
       IF(kuo(i) > 0) cumulus(i)=.TRUE.
    ENDDO
    !180-----pi
    !  x     y
    !  y = x*pi/180
    delta_lambda=(360.0_r8/iMax)*(3.1415926e0_r8/180.0_r8)
    delta_phi   =(360.0_r8/iMax)*(3.1415926e0_r8/180.0_r8)


    L_CLD_AREA    = .TRUE.! true if using area cloud fraction (ACF)
    L_ACF_Cusack  = .TRUE.! ... and selected Cusack and PC2 off
    L_ACF_Brooks  = .TRUE.! ... and selected Brooks
    L_eacf        = .TRUE.! true if using empirically adjusted cloud fraction
    L_conv4a      = .FALSE.! true if using 4A convection scheme
    L_mcr_qcf2    = .FALSE.     !, INTENT(IN   )IN Use second prognostic ice if T
    L_mixing_ratio= .TRUE.     !, INTENT(IN   )IN Use mixing ratios removed by deposition depending on

    CALL qt_bal_cld ( &
         ps                    , &!Real   ,INTENT(IN   ) :: p_star        (1:nCols)
         p_theta_levels        , &!Real   ,INTENT(IN   ) :: p_theta_levels(1:nCols,LEVELS)
         p                     , &!Real   ,INTENT(IN   ) :: p             (1:nCols,LEVELS+1)
         T3                    , &!Real   ,INTENT(INOUT) :: T3            (1:nCols,LEVELS)
         q                     , &!Real   ,INTENT(INOUT) :: q             (1:nCols,LEVELS)
         qcl                   , &!Real   ,INTENT(INOUT) :: qcl           (1:nCols,LEVELS)
         qcf                   , &!Real   ,INTENT(IN   ) :: qcf           (1:nCols,LEVELS)
         qcf2                  , &!Real   ,INTENT(IN   ) :: qcf2          (1:nCols,LEVELS)
         rhcrit                , &!Real   ,INTENT(IN   ) :: rhcpt         (1:nCols,LEVELS)
         nCols                 , &!Integer,INTENT(IN   ) :: nCols
         bl_levels             , &!Integer,INTENT(IN   ) :: bl_levels
         cloud_fraction_method , &!Integer,Intent(IN   ) :: cloud_fraction_method  ! Method for calculating total cloud fract
         overlap_ice_liquid    , &!Real   ,Intent(IN   ) :: overlap_ice_liquid   ! Overlap between ice and liquid phases
         ice_fraction_method   , &!Integer,Intent(IN   ) :: ice_fraction_method        ! Method for calculating ice cloud frac.
         ctt_weight            , &!Real   ,Intent(IN   ) :: ctt_weight              ! Weighting of cloud top temperature
         t_weight              , &!Real   ,Intent(IN   ) :: t_weight              ! Weighting of local temperature
         qsat_fixed            , &!Real   ,Intent(IN   ) :: qsat_fixed              ! Fixed value of saturation humidity
         sub_cld               , &!Real   ,Intent(IN   ) :: sub_cld              ! Scaling parameter
         LEVELS                , &!Integer,INTENT(IN   ) :: LEVELS                    ! number of model levels where moisture
         delta_lambda          , &!Real   ,INTENT(IN   ) :: delta_lambda
         delta_phi             , &!Real   ,INTENT(IN   ) :: delta_phi
         r_theta_levels        , &!Real   ,INTENT(IN   ) :: r_theta_levels(1:nCols,0:LEVELS)
         FV_cos_theta_latitude , &!Real   ,INTENT(IN   ) :: FV_cos_theta_latitude (1:nCols)
         lc                    , &!Real   ,INTENT(IN   ) :: lc
         cp                    , &!Real   ,INTENT(IN   ) :: cp
         L_cld_area            , &!LOGICAL,INTENT(IN   ) :: L_CLD_AREA        ! true if using area cloud fraction (ACF)
         L_ACF_Cusack          , &!LOGICAL,INTENT(IN   ) :: L_ACF_Cusack  ! ... and selected Cusack and PC2 off
         L_ACF_Brooks          , &!LOGICAL,INTENT(IN   ) :: L_ACF_Brooks  ! ... and selected Brooks
         L_eacf                , &!LOGICAL,INTENT(IN   ) :: L_eacf        ! true if using empirically adjusted cloud fraction
         L_mcr_qcf2            , &!LOGICAL,INTENT(IN   ) :: L_mcr_qcf2        ! true if second cloud ice variable in use
         l_mixing_ratio        , &!LOGICAL,INTENT(IN   ) :: l_mixing_ratio! true if using mixing ratio formulation
         ntml                  , &!Integer,Intent(IN   ) :: ntml (nCols)
         cumulus               , &!Logical,Intent(IN   ) :: cumulus (nCols) ! bl convection flag
         L_conv4a              , &!LOGICAL,INTENT(IN   ) :: L_conv4a        ! true if using 4A convection scheme
         area_cloud_fraction   , &!Real   ,INTENT(OUT  ) :: area_cloud_fraction  (1:nCols, LEVELS)
         bulk_cloud_fraction   , &!Real   ,INTENT(OUT  ) :: bulk_cloud_fraction  (1:nCols, LEVELS)
         cloud_fraction_liquid , &!Real   ,INTENT(OUT  ) :: cloud_fraction_liquid(1:nCols, LEVELS)
         cloud_fraction_frozen   )!Real   ,INTENT(OUT  ) :: cloud_fraction_frozen(1:nCols, LEVELS)
    DO k=1,LEVELS
       DO i=1,nCols
          CF (i,k)= bulk_cloud_fraction  (i, k)! Cloud fraction at processed levels (decimal fraction).
          CFL(i,k)= cloud_fraction_liquid(i, k)! Liquid cloud fraction at processed levels (decimal fraction).
          CFF(i,k)= cloud_fraction_frozen(i, k)! Frozen cloud fraction at processed levels (decimal fraction).
       END DO
    END DO

  END SUBROUTINE Run_CloudFraction_UKME

  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  ! Subroutine qt_bal_cld
  SUBROUTINE qt_bal_cld (&!
       p_star                , &!Real  , INTENT(IN   ) ::  p_star        (1:nCols)
       p_theta_levels        , &!Real  , INTENT(IN   ) ::  p_theta_levels(1:nCols,LEVELS)
       p                     , &!Real  , INTENT(IN   ) ::  p             (1:nCols,LEVELS+1)
       T3                    , &!Real  , INTENT(INOUT) ::  T3            (1:nCols,LEVELS)
       q                     , &!Real  , INTENT(INOUT) ::  q             (1:nCols,LEVELS)
       qcl                   , &!Real  , INTENT(INOUT) ::  qcl           (1:nCols,LEVELS)
       qcf                   , &!Real  , INTENT(IN   ) ::  qcf           (1:nCols,LEVELS)
       qcf2                  , &!Real  , INTENT(IN   ) ::  qcf2          (1:nCols,LEVELS)
       rhcpt                 , &!Real  , INTENT(IN   ) ::  rhcpt         (1:nCols,LEVELS)
       nCols                 , &!Integer, INTENT(IN   ) ::  nCols
       bl_levels             , &!Integer, INTENT(IN   ) ::  bl_levels
       cloud_fraction_method , &!Integer, Intent(IN  ) ::   cloud_fraction_method  ! Method for calculating total cloud fraction
       overlap_ice_liquid    , &!Real    , Intent(IN ) :: overlap_ice_liquid   ! Overlap between ice and liquid phases
       ice_fraction_method   , &!Integer, Intent(IN  ) :: ice_fraction_method    ! Method for calculating ice cloud frac.
       ctt_weight            , &!Real    , Intent(IN ) ::  ctt_weight          ! Weighting of cloud top temperature
       t_weight              , &!Real    , Intent(IN ) ::  t_weight            ! Weighting of local temperature
       qsat_fixed            , &!Real    , Intent(IN ) ::  qsat_fixed          ! Fixed value of saturation humidity
       sub_cld               , &!Real    , Intent(IN ) ::  sub_cld             ! Scaling parameter
       LEVELS                , &!Integer, INTENT(IN   ) ::  LEVELS                   ! number of model levels where moisture
       delta_lambda          , &!Real  , INTENT(IN   ) ::  delta_lambda
       delta_phi             , &!Real  , INTENT(IN   ) ::  delta_phi
       r_theta_levels        , &!Real  , INTENT(IN   ) ::  r_theta_levels(1:nCols,0:LEVELS)
       FV_cos_theta_latitude , &!Real  , INTENT(IN   ) ::  FV_cos_theta_latitude (1:nCols)
       lc                    , &!Real , INTENT(IN   ) ::  lc
       cp                    , &!Real , INTENT(IN   ) ::  cp
       L_cld_area            , &!LOGICAL, INTENT(IN   ) :: L_CLD_AREA    ! true if using area cloud fraction (ACF)
       L_ACF_Cusack          , &!LOGICAL, INTENT(IN   ) :: L_ACF_Cusack  ! ... and selected Cusack and PC2 off
       L_ACF_Brooks          , &!LOGICAL, INTENT(IN   ) :: L_ACF_Brooks  ! ... and selected Brooks
       L_eacf                , &!LOGICAL, INTENT(IN   ) :: L_eacf        ! true if using empirically adjusted cloud fraction
       L_mcr_qcf2            , &!LOGICAL, INTENT(IN   ) :: L_mcr_qcf2    ! true if second cloud ice variable in use
       l_mixing_ratio        , &!LOGICAL, INTENT(IN   ) :: l_mixing_ratio! true if using mixing ratio formulation
       ntml                  , &!Integer     , Intent(IN) ::  ntml (nCols)
       cumulus               , &!Logical      , Intent(IN) :: cumulus (nCols) ! bl convection flag
       L_conv4a              , &!LOGICAL, INTENT(IN   ) :: L_conv4a      ! true if using 4A convection scheme
       area_cloud_fraction   , &!Real   , INTENT(OUT) ::   area_cloud_fraction  (1:nCols, LEVELS)
       bulk_cloud_fraction   , &!Real   , INTENT(OUT) ::   bulk_cloud_fraction  (1:nCols,LEVELS)
       cloud_fraction_liquid , &!Real   , INTENT(OUT) ::   cloud_fraction_liquid(1:nCols, LEVELS)
       cloud_fraction_frozen   )!Real   , INTENT(OUT) ::   cloud_fraction_frozen(1:nCols, LEVELS)

    ! Purpose:
    !        reset q, t and the cloud fields to be consistent at the
    !        end of the timestep
    !
    ! Method:
    !          Is described in ;
    !
    ! Original Progammer: Andy Malcolm
    ! Current code owner: Andy Malcolm
    !
    ! History:
    ! Version   Date       Comment
    ! -------  -------     -------
    !  5.4     28/08/02    New Deck                   Andy Malcolm
    !  5.5     03/02/03    Add qcf2 to qcf if active.      Richard Forbes
    !  5.5     20/02/03    Inserted missing #endif    P.Dando
    !  6.2     07/11/05   Pass cloud variables from UMUI. D. Wilson
    !  6.4     08/11/06    Include Brooks Area Cloud Fraction D. Wilson
    !  6.4     10/01/07    Add logical control for mixing ratios. D. Wilson
    !
    ! Code Description:
    !   Language: FORTRAN 77 + CRAY extensions
    !   This code is written to UMDP3 programming standards.
    !

    IMPLICIT NONE

    INTEGER, INTENT(IN   ) ::  LEVELS                   ! number of model levels where moisture
    ! variables are held
    INTEGER, INTENT(IN   ) ::  nCols
    INTEGER, INTENT(IN   ) ::  bl_levels
    ! physical constants
    REAL(KIND=r8) , INTENT(IN   ) ::  lc
    REAL(KIND=r8) , INTENT(IN   ) ::  cp

    LOGICAL, INTENT(IN   ) :: L_CLD_AREA    ! true if using area cloud fraction (ACF)
    LOGICAL, INTENT(IN   ) :: L_ACF_Cusack  ! ... and selected Cusack and PC2 off
    LOGICAL, INTENT(IN   ) :: L_ACF_Brooks  ! ... and selected Brooks
    LOGICAL, INTENT(IN   ) :: L_eacf        ! true if using empirically adjusted cloud fraction
    LOGICAL, INTENT(IN   ) :: L_conv4a      ! true if using 4A convection scheme
    LOGICAL, INTENT(IN   ) :: L_mcr_qcf2    ! true if second cloud ice variable in use
    LOGICAL, INTENT(IN   ) :: l_mixing_ratio! true if using mixing ratio formulation

    REAL(KIND=r8)  , INTENT(IN   ) ::  p             (1:nCols,LEVELS+1)
    REAL(KIND=r8)  , INTENT(IN   ) ::  p_theta_levels(1:nCols,LEVELS)
    REAL(KIND=r8)  , INTENT(IN   ) ::  p_star        (1:nCols)
    REAL(KIND=r8)  , INTENT(INOUT) ::  T3            (1:nCols,LEVELS)
    REAL(KIND=r8)  , INTENT(INOUT) ::  q             (1:nCols,LEVELS)
    REAL(KIND=r8)  , INTENT(INOUT) ::  qcl           (1:nCols,LEVELS)
    REAL(KIND=r8)  , INTENT(IN   ) ::  qcf           (1:nCols,LEVELS)
    REAL(KIND=r8)  , INTENT(IN   ) ::  qcf2          (1:nCols,LEVELS)
    REAL(KIND=r8)  , INTENT(IN   ) ::  rhcpt         (1:nCols,LEVELS)
    ! coordinate arrays
    REAL(KIND=r8)  , INTENT(IN   ) ::  r_theta_levels(1:nCols,0:LEVELS)
    REAL(KIND=r8)  , INTENT(IN   ) ::  delta_lambda
    REAL(KIND=r8)  , INTENT(IN   ) ::  delta_phi
    ! trig arrays
    REAL(KIND=r8)  , INTENT(IN   ) ::  FV_cos_theta_latitude (1:nCols)
    !
    INTEGER, INTENT(IN  ) :: cloud_fraction_method  ! Method for calculating total cloud fraction
    INTEGER, INTENT(IN  ) :: ice_fraction_method    ! Method for calculating ice cloud frac.

    REAL(KIND=r8)    , INTENT(IN ) :: overlap_ice_liquid   ! Overlap between ice and liquid phases
    REAL(KIND=r8)    , INTENT(IN ) :: ctt_weight          ! Weighting of cloud top temperature
    REAL(KIND=r8)    , INTENT(IN ) :: t_weight            ! Weighting of local temperature
    REAL(KIND=r8)    , INTENT(IN ) :: qsat_fixed          ! Fixed value of saturation humidity
    REAL(KIND=r8)    , INTENT(IN ) :: sub_cld             ! Scaling parameter

    ! Diagnostic variables

    REAL(KIND=r8)   , INTENT(OUT) :: area_cloud_fraction  (1:nCols, LEVELS)
    REAL(KIND=r8)   , INTENT(OUT) :: bulk_cloud_fraction  (1:nCols, LEVELS)!CF (nCols,LEVELS) Cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)   , INTENT(OUT) :: cloud_fraction_liquid(1:nCols, LEVELS)!CFL(nCols,LEVELS) Liquid cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)   , INTENT(OUT) :: cloud_fraction_frozen(1:nCols, LEVELS)!CFF(nCols,LEVELS) Frozen cloud fraction at processed levels (decimal fraction).

    INTEGER     , INTENT(IN) ::  ntml (nCols)

    LOGICAL      , INTENT(IN) :: cumulus (nCols) ! bl convection flag

    ! Local variables
    REAL(KIND=r8) ::  p_layer_centres(nCols, 0:LEVELS) 
    ! pressure at layer centres. Same as p_theta_levels
    !except bottom level = p_star, and at top = 0.
    REAL(KIND=r8) ::  p_layer_boundaries(nCols, 0:LEVELS) 
    REAL(KIND=r8) ::  tl       (nCols, LEVELS) 
    REAL(KIND=r8) ::  qt       (nCols, LEVELS) 
    REAL(KIND=r8) ::  qcf_in   (nCols, LEVELS) 
    REAL(KIND=r8) ::  qcl_out  (nCols, LEVELS)
    REAL(KIND=r8) ::  cf_inout (nCols, LEVELS)
    REAL(KIND=r8) ::  cfl_inout(nCols, LEVELS)
    REAL(KIND=r8) ::  cff_inout(nCols, LEVELS)


    INTEGER ::  large_levels
    INTEGER ::  levels_per_level

    INTEGER :: i,k,errorstatus     ! loop variables

    DO i = 1, nCols
       p_layer_centres(i,0)    = p_star(i)
       p_layer_boundaries(i,0) = p_star(i)
    END DO

    DO k = 1, LEVELS-1
       DO i = 1, nCols
          p_layer_centres   (i,k) = p_theta_levels(i,k)
          p_layer_boundaries(i,k) = p(i,k+1)
       END DO
    END DO
    k=LEVELS
    DO i = 1, nCols
       p_layer_centres(i,k) = p_theta_levels(i,k)
       p_layer_boundaries(i,k) = 0.0_r8
    END DO

    ! ----------------------------------------------------------------------
    ! Section  Convert qT and Tl for input to cloud scheme.
    ! ----------------------------------------------------------------------

    ! Create Tl and qT
    DO k = 1, LEVELS
       DO i = 1, nCols
          ! Tl(i,k)       = (theta(i,k) *exner_theta_levels(i,k)) - (lc * qcl(i,k)) / cp
          Tl(i,k)       = T3(i,k) - (lc * qcl(i,k)) / cp

          qt(i,k)       = q(i,k) + qcl(i,k)
          qcf_in(i,k)   = qcf(i,k)
          cf_inout(i,k) = bulk_cloud_fraction(i,k)
          cfl_inout(i,k)= cloud_fraction_liquid(i,k)
          cff_inout(i,k)= cloud_fraction_frozen(i,k)
       END DO
    END DO

    ! If second cloud ice variable in use then add to qcf
    ! for the cloud scheme call
    IF (L_mcr_qcf2) qcf_in(:,:) = qcf_in(:,:) + qcf2(1:nCols,:)

    ! ----------------------------------------------------------------------
    ! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
    !              calculate bulk_cloud fields from qT and qcf
    !               and calculate area_cloud fields.
    ! ----------------------------------------------------------------------
    !

    ! Determine number of sublevels for vertical gradient area cloud
    ! Want an odd number of sublevels per level: 3 is hardwired in do loops
    levels_per_level = 3
    large_levels = ((LEVELS - 2)*levels_per_level) + 2
    !
    ! DEPENDS ON: ls_arcld
    CALL ls_arcld( &
                                !      Pressure related fields
         p_layer_centres       , &!REAL    , INTENT(IN) :: p_layer_centres(nCols,0:LEVELS)
         RHCPT                 , &!REAL    , INTENT(IN) :: rhcrit(nCols,LEVELS)
         p_layer_boundaries    , &!REAL    , INTENT(IN) :: p_layer_boundaries(nCols,0:LEVELS)
                                !      Array dimensions
         LEVELS                , &!INTEGER , INTENT(IN)  ::  LEVELS
         nCols                 , &!INTEGER , INTENT(IN)  ::  nCols
         bl_levels             , &!INTEGER , INTENT(IN)  ::  bl_levels
         cloud_fraction_method , &!Integer , Intent(IN) :: cloud_fraction_method ! Method for calculating total cloud fraction
         overlap_ice_liquid    , &!Real    , Intent(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
         ice_fraction_method   , &!Integer , Intent(IN) :: ice_fraction_method      ! Method for calculating ice cloud frac.
         ctt_weight            , &!Real    , Intent(IN) :: ctt_weight          ! Weighting of cloud top temperature
         t_weight              , &!Real    , Intent(IN) :: t_weight          ! Weighting of local temperature
         qsat_fixed            , &!Real    , Intent(IN) :: qsat_fixed          ! Fixed value of saturation humidity
         sub_cld               , &!Real    , Intent(IN) :: sub_cld           ! Scaling parameter
         levels_per_level      , &!INTEGER , INTENT(IN)  ::  levels_per_level
         large_levels          , &!INTEGER , INTENT(IN)  ::  large_levels
                                !      Switch on area cloud calcuculation and select which to use
         L_cld_area            , &!LOGICAL , INTENT(IN) :: L_AREA_CLOUD
         L_ACF_Cusack          , &!LOGICAL , INTENT(IN) :: L_ACF_Cusack
         L_ACF_Brooks          , &!LOGICAL , INTENT(IN) :: L_ACF_Brooks
         L_eacf                , &!LOGICAL , INTENT(IN) :: L_eacf
                                !      Needed for LS_ACF_Brooks
         delta_lambda          , &!REAL    , INTENT(IN) ::  delta_lambda     ! EW (x) grid spacing in radians
         delta_phi             , &!REAL    , INTENT(IN) ::  delta_phi         ! NS (y) grid spacing in radians
         r_theta_levels        , &!REAL    , INTENT(IN) ::  r_theta_levels (1:nCols,0:LEVELS) ! height of theta levels (from centre of earth)
         FV_cos_theta_latitude , &!REAL    , INTENT(IN) ::  FV_cos_theta_latitude (1:nCols)
                                !      Convection diagnosis inforormation (only used for a05_4a)
         ntml                  , &!INTEGER , INTENT(IN)  ::  ntml(nCols)
         cumulus               , &!LOGICAL , INTENT(IN) :: cumulus(nCols)
         L_conv4a              , &!LOGICAL , INTENT(IN) :: L_conv4a 
         l_mixing_ratio        , &!LOGICAL , INTENT(IN) :: L_mixing_ratio
                                !      Prognostic Fields
         qcf_in                , &!REAL    , INTENT(IN) :: qcf_latest(nCols,LEVELS)
         Tl                    , &!REAL    , INTENT(INOUT) :: T_latest(nCols,LEVELS)
         qt                    , &!REAL    , INTENT(INOUT) :: q_latest(nCols,LEVELS)
         qcl_out               , &!REAL  , INTENT(OUT)  ::  qcl_latest(nCols,LEVELS)
                                !      Various cloud fractions
         area_cloud_fraction   , &!REAL  , INTENT(OUT)  ::  area_cloud_fraction(nCols,LEVELS)
         cf_inout              , &!REAL  , INTENT(OUT)  ::  bulk_cloud_fraction(nCols,LEVELS)
         cfl_inout             , &!REAL  , INTENT(OUT)  ::  cloud_fraction_liquid(nCols,LEVELS)
         cff_inout             , &!REAL  , INTENT(OUT)  ::  cloud_fraction_frozen(nCols,LEVELS)
         errorstatus             )!INTEGER,INTENT(OUT)  ::  error_code

    ! qt holds q (no halos), tl holds t(no halos),
    ! qcl_out holds qcl(no halos)
    DO k = 1, LEVELS
       DO i = 1, nCols
          !theta(i,k)                 = tl(i,k)/exner_theta_levels(i,k)
          T3(i,k)                    = tl(i,k)
          q(i,k)                     = qt(i,k)
          qcl(i,k)                   = qcl_out(i,k)
          bulk_cloud_fraction(i,k)   = cf_inout(i,k)
          cloud_fraction_liquid(i,k) = cfl_inout(i,k)
          cloud_fraction_frozen(i,k) = cff_inout(i,k)
       END DO
    END DO

    RETURN
  END SUBROUTINE qt_bal_cld




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Large-scale Area (Vertical Gradient) Cloud Scheme.
  ! Subroutine Interface:
  SUBROUTINE LS_ARCLD(    &
                                !      Pressure related fields
       p_layer_centres       ,&!REAL    , INTENT(IN) :: p_layer_centres(nCols,0:LEVELS)
       rhcrit                ,&!REAL    , INTENT(IN) :: rhcrit(nCols,LEVELS)
       p_layer_boundaries    ,&!REAL    , INTENT(IN) :: p_layer_boundaries(nCols,0:LEVELS)
                                !      Array dimensions
       LEVELS                ,&!INTEGER , INTENT(IN)  ::  LEVELS
       nCols                 ,&!INTEGER , INTENT(IN)  ::  nCols
       bl_levels             ,&!INTEGER , INTENT(IN)  ::  bl_levels
       cloud_fraction_method ,&!Integer , Intent(IN) :: cloud_fraction_method ! Method for calculating total cloud fraction
       overlap_ice_liquid    ,&!Real    , Intent(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
       ice_fraction_method   ,&!Integer , Intent(IN) :: ice_fraction_method    ! Method for calculating ice cloud frac.
       ctt_weight            ,&!Real    , Intent(IN) :: ctt_weight         ! Weighting of cloud top temperature
       t_weight              ,&!Real    , Intent(IN) :: t_weight           ! Weighting of local temperature
       qsat_fixed            ,&!Real    , Intent(IN) :: qsat_fixed         ! Fixed value of saturation humidity
       sub_cld               ,&!Real    , Intent(IN) :: sub_cld                ! Scaling parameter
       levels_per_level      ,&!INTEGER , INTENT(IN)  ::  levels_per_level
       large_levels          ,&!INTEGER , INTENT(IN)  ::  large_levels
                                !      Switch on area cloud calculation and select which to use
       L_AREA_CLOUD          ,&!LOGICAL , INTENT(IN) :: L_AREA_CLOUD
       L_ACF_Cusack          ,&!LOGICAL , INTENT(IN) :: L_ACF_Cusack
       L_ACF_Brooks          ,&!LOGICAL , INTENT(IN) :: L_ACF_Brooks
       L_eacf                ,&!LOGICAL , INTENT(IN) :: L_eacf
                                !      Needed for LS_ACF_Brooks
       delta_lambda          ,&!REAL    , INTENT(IN) ::  delta_lambda     ! EW (x) grid spacing in radians
       delta_phi             ,&!REAL    , INTENT(IN) ::  delta_phi        ! NS (y) grid spacing in radians
       r_theta_levels        ,&!REAL    , INTENT(IN) ::  r_theta_levels (1:nCols,0:LEVELS) ! height of theta levels (from centre of earth)
       FV_cos_theta_latitude ,&!REAL    , INTENT(IN) ::  FV_cos_theta_latitude (1:nCols)
                                !      Convection diagnosis information (only used for a05_4a)
       ntml                  ,&!INTEGER , INTENT(IN)  ::  ntml(nCols)
       cumulus               ,&!LOGICAL , INTENT(IN) :: cumulus(nCols)
       L_conv4a              ,&!LOGICAL , INTENT(IN) :: L_conv4a 
       L_mixing_ratio        ,&!LOGICAL , INTENT(IN) :: L_mixing_ratio
                                !      Prognostic Fields
       qcf_latest            ,&!REAL    , INTENT(IN) :: qcf_latest(nCols,LEVELS)
       T_latest              ,&!REAL    , INTENT(INOUT) :: T_latest(nCols,LEVELS)
       q_latest              ,&!REAL    , INTENT(INOUT) :: q_latest(nCols,LEVELS)
       qcl_latest            ,&!REAL  , INTENT(OUT)  ::  qcl_latest(nCols,LEVELS)
                                !      Various cloud fractions
       area_cloud_fraction   ,&!REAL  , INTENT(OUT)  ::  area_cloud_fraction(nCols,LEVELS)
       bulk_cloud_fraction   ,&!REAL  , INTENT(OUT)  ::  bulk_cloud_fraction(nCols,LEVELS)
       cloud_fraction_liquid ,&!REAL  , INTENT(OUT)  ::  cloud_fraction_liquid(nCols,LEVELS)
       cloud_fraction_frozen ,&!REAL  , INTENT(OUT)  ::  cloud_fraction_frozen(nCols,LEVELS)
       error_code             )!INTEGER,INTENT(OUT)  ::  error_code
    !
    IMPLICIT NONE
    !
    ! Purpose:
    !   This subroutine calculates liquid and ice cloud fractional cover
    !   for use with the enhanced precipitation microphysics scheme. It
    !   also returns area and bulk cloud fractions for use in radiation.
    !
    ! Method:
    !   Statistical cloud scheme separates input moisture into specific
    !   humidity and cloud liquid water. Temperature calculated from liquid
    !   water temperature. Cloud fractions calculated from statistical
    !   relation between cloud fraction and cloud liquid/ice water content.
    !   Area cloud fraction calculated by subdividing layers, calculating
    !   cloud on each sublayer and taking maximum (use mean for bulk cloud).
    !
    ! Current Owner of Code: A. C. Bushell
    !
    ! History:
    ! Version   Date     Comment
    !  4.5    14-05-98   Original Code (Stephen Cusack for HadAM4)
    !
    !  5.0    14-05-98   Rewritten for New Dynamics with simplified
    !                    interpolation.  AC Bushell
    !
    !  5.1    15-12-99   Change to allow 3D RHcrit specification. AC Bushell
    !
    !  5.2    01-11-00   Remove temporary test diagnostics.       AC Bushell
    !
    !
    !  5.4   10-09-02    If cumulus convection is diagnosed the large
    !                    scale cloud is not allowed in the levels around
    !                    cumulus cloud-base.           Gill Martin
    !  6.2   07-11-05    Pass cloud variables from UMUI.  D. Wilson
    !  6.4   08/11/06    Include Brooks Area Cloud Fraction D. Wilson
    !  6.4   07-08-06    Enable mixing ratios to be used.  D. Wilson
    ! Description of Code:
    !   FORTRAN 77  + common extensions also in Fortran90.
    !   This code is written to UMDP3 version 6 programming standards.
    !
    !   Documentation: UMDP No.29
    !
    !  Global Variables:----------------------------------------------------
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
    ! C_LHEAT start

    ! latent heat of condensation of water at 0degc
    REAL(KIND=r8),PARAMETER:: LC=2.501E6_r8

    ! latent heat of fusion at 0degc
    REAL(KIND=r8),PARAMETER:: LF=0.334E6_r8

    ! C_LHEAT end
    !
    !  Subroutine Arguments:------------------------------------------------
    !
    ! arguments with intent in. ie: input variables.
    !
    INTEGER , INTENT(IN)  ::  LEVELS
    !       No. of levels being processed by cloud scheme.
    INTEGER , INTENT(IN)  ::  bl_levels
    !       No. of boundary layer levels
    INTEGER , INTENT(IN)  ::  nCols
    !       Horizontal dimensions of rh_crit diagnostic.
    INTEGER , INTENT(IN)  ::  levels_per_level
    !       No. of sub-levels being processed by area cloud scheme.
    !       Want an odd number of sublevels per level.
    !       NB: levels_per_level = 3 is currently hardwired in the do loops
    INTEGER , INTENT(IN)  ::  large_levels! Total no. of sub-levels being processed by cloud scheme.
    !       Currently ((LEVELS - 2)*levels_per_level) + 2
    !
    INTEGER , INTENT(IN)  ::  ntml(nCols)!  Height of diagnosed BL top
    LOGICAL , INTENT(IN) :: L_AREA_CLOUD ! true if using area cloud fraction (ACF)
    LOGICAL , INTENT(IN) :: L_ACF_Cusack ! ... and selected Cusack and PC2 off
    LOGICAL , INTENT(IN) :: L_ACF_Brooks ! ... and selected Brooks
    LOGICAL , INTENT(IN) :: L_eacf       ! true if using empirically adjusted cloud fraction
    LOGICAL , INTENT(IN) :: L_mixing_ratio ! true if using mixing rations
    !
    LOGICAL , INTENT(IN) :: L_conv4a        ! true if using 4A convection scheme
    LOGICAL , INTENT(IN) :: cumulus(nCols)   ! Logical indicator of convection

    REAL(KIND=r8)    , INTENT(IN) :: qcf_latest(nCols,LEVELS)
    !       Cloud ice content at processed levels (kg water per kg air).
    REAL(KIND=r8)    , INTENT(IN) :: p_layer_centres(nCols,0:LEVELS)
    !       pressure at all points, on theta levels (Pa).
    REAL(KIND=r8)    , INTENT(IN) :: rhcrit(nCols,LEVELS)
    !       Critical relative humidity.  See the the paragraph incorporating
    !       eqs P292.11 to P292.14; the values need to be tuned for a given
    !       set of levels.
    REAL(KIND=r8)    , INTENT(IN) :: p_layer_boundaries(nCols,0:LEVELS)
    !       pressure at all points, on u,v levels (Pa).
    REAL(KIND=r8)    , INTENT(IN) ::  r_theta_levels (1:nCols,0:LEVELS) ! height of theta levels (from centre of earth)
    REAL(KIND=r8)    , INTENT(IN) ::  FV_cos_theta_latitude (1:nCols)
    ! Finite volume cos(lat)
    REAL(KIND=r8)    , INTENT(IN) ::  delta_lambda     ! EW (x) grid spacing in radians
    REAL(KIND=r8)    , INTENT(IN) ::  delta_phi        ! NS (y) grid spacing in radians

    !
    INTEGER , INTENT(IN) :: cloud_fraction_method ! Method for calculating total cloud fraction
    INTEGER , INTENT(IN) :: ice_fraction_method    ! Method for calculating ice cloud frac.

    REAL(KIND=r8)    , INTENT(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
    REAL(KIND=r8)    , INTENT(IN) :: ctt_weight         ! Weighting of cloud top temperature
    REAL(KIND=r8)    , INTENT(IN) :: t_weight           ! Weighting of local temperature
    REAL(KIND=r8)    , INTENT(IN) :: qsat_fixed         ! Fixed value of saturation humidity
    REAL(KIND=r8)    , INTENT(IN) :: sub_cld                ! Scaling parameter
    !
    ! arguments with intent in/out. ie: input variables changed on output.
    !
    REAL(KIND=r8)    , INTENT(INOUT) :: q_latest(nCols,LEVELS)
    !       On input : Vapour + liquid water content (QW) (kg per kg air).
    !       On output: Specific humidity at processed levels
    !                   (kg water per kg air).
    REAL(KIND=r8)    , INTENT(INOUT) :: T_latest(nCols,LEVELS)
    !       On input : Liquid water temperature (TL) (K).
    !       On output: Temperature at processed levels (K).
    !
    ! arguments with intent out. ie: output variables.
    !
    !     Error Status:
    INTEGER,INTENT(OUT)  ::  error_code  !, 0 if OK; 1 if bad arguments.
    !
    REAL(KIND=r8)  , INTENT(OUT)  ::  qcl_latest(nCols,LEVELS)
    !       Cloud liquid content at processed levels (kg water per kg air).
    REAL(KIND=r8)  , INTENT(OUT)  ::  area_cloud_fraction(nCols,LEVELS)
    !       Area cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)  , INTENT(OUT)  ::  bulk_cloud_fraction(nCols,LEVELS)
    !       Cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)  , INTENT(OUT)  ::  cloud_fraction_liquid(nCols,LEVELS)
    !       Liquid cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)  , INTENT(OUT)  ::  cloud_fraction_frozen(nCols,LEVELS)
    !       Frozen cloud fraction at processed levels (decimal fraction).
    !
    !  Local parameters and other physical constants------------------------
    REAL(KIND=r8),PARAMETER  :: LCRCP =LC/CP! Derived parameters. ! Lat ht of condensation/Cp.
    REAL(KIND=r8),PARAMETER  :: drat_thresh =3.0E-1_r8           ! Test for continuity of sub-levels
    REAL(KIND=r8),PARAMETER  :: tol_test    =1.0E-11_r8           ! Tolerance for non-zero humidities
    !
    !  Local scalars--------------------------------------------------------
    !
    !  (a) Scalars effectively expanded to workspace by the Cray (using
    !      vector registers).
    REAL(KIND=r8) ::  qt_norm_next
    ! Temporary space for qT_norm
    REAL(KIND=r8) ::   stretcher
    !
    REAL(KIND=r8) ::   inverse_level
    ! Set to (1._r8 / levels_per_level)
    REAL(KIND=r8) ::   delta_p        ! Layer pressure thickness * inverse_level
    !
    !  (b) Others.
    INTEGER i,k      ! Loop counters: k - vertical level index.
    !                                       i,j - horizontal field index.
    INTEGER k_index    ! Extra loop counter for large arrays.
    !
    !  Local dynamic arrays-------------------------------------------------
    !    11 blocks of real workspace are required.
    REAL(KIND=r8)  :: qsl(nCols)
    !        Saturated specific humidity for temp TL or T.
    REAL(KIND=r8) ::   qt_norm(nCols)
    !        Total water content normalized to qSAT_WAT.
    REAL(KIND=r8) ::   rhcrit_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   p_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   T_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   q_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   qcf_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   cloud_fraction_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   qcl_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   cloud_fraction_liquid_large(nCols,large_levels)
    !
    REAL(KIND=r8) ::   cloud_fraction_frozen_large(nCols,large_levels)
    !
    INTEGER :: ntml_large(nCols)

    LOGICAL :: linked(nCols,LEVELS)
    !       True for sub-layers that have similar supersaturation properties
    !
    !  External subroutine calls: ------------------------------------------
    !      EXTERNAL QSAT_WAT,LS_CLD,LS_ACF_BROOKS
    !- End of Header
    ! ----------------------------------------------------------------------
    !  Check input arguments for potential over-writing problems.
    ! ----------------------------------------------------------------------
    error_code=0
    !
    ! ==Main Block==--------------------------------------------------------
    ! Subroutine structure :
    ! Loop round levels to be processed.
    ! ----------------------------------------------------------------------
    ! Lsarc_if1:
    IF (.NOT. L_AREA_CLOUD) THEN
       !     As before
       !
       ! DEPENDS ON: ls_cld
       CALL ls_cld(&
                                !Pressure related fields
            p_layer_centres(1,1)  , &  !REAL    , INTENT(IN) :: p_theta_levels(nCols,levels)! pressure at all points (Pa).
            rhcrit                , &  !REAL    , INTENT(IN) :: RHCRIT(nCols,LEVELS)
                                !      Array dimensions
            LEVELS                , &  !INTEGER , INTENT(IN) :: LEVELS!          No. of levels being processed.
            bl_levels             , &  !INTEGER , INTENT(IN) :: BL_LEVELS!       No. of boundary layer levels
            nCols                 , &  !INTEGER , INTENT(IN) :: nCols
            cloud_fraction_method , &  !Integer , Intent(IN) :: cloud_fraction_method ! Method for calculating
            overlap_ice_liquid    , &  !Real    , Intent(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
            ice_fraction_method   , &  !Integer , Intent(IN) :: ice_fraction_method   ! Method for calculating ice 
            ctt_weight            , &  !Real    , Intent(IN) :: ctt_weight! Weighting of cloud top temperature
            t_weight              , &  !Real    , Intent(IN) :: t_weight  ! Weighting of local temperature
            qsat_fixed            , &  !Real    , Intent(IN) :: qsat_fixed ! Fixed value of saturation humidity
            sub_cld               , &  !Real    , Intent(IN) :: sub_cld                   ! Scaling parameter
                                !      From convection diagnosis (s (only used if a05_4a)
            ntml                  , &  !INTEGER, INTENT(IN   ) :: NTML(nCols)     ! IN Height of diagnosed BL top
            cumulus               , &  !LOGICAL, INTENT(IN   ) ::  CUMULUS(nCols)  ! IN Logical indicator of convection
            L_conv4a              , &  !LOGICAL, INTENT(IN   ) ::  L_conv4a! IN true if using 4A convection scheme
            L_eacf                , &  !LOGICAL, INTENT(IN   ) ::  L_eacf! IN true if using empirically adjusted cloud
            L_mixing_ratio        , &  !LOGICAL, INTENT(IN   ) ::  L_mixing_ratio  ! IN true if using mixing ratios
                                !      Prognostic Fields
            T_latest              , &  !REAL   , INTENT(INOUT) ::        T(nCols,LEVELS)
            bulk_cloud_fraction   , &  !REAL, INTENT(OUT)        ::  CF
            q_latest              , &  !REAL   , INTENT(INOUT) ::  Q(nCols,LEVELS)
            qcf_latest            , &  !REAL    , INTENT(IN) :: QCF(nCols,LEVELS)!Cloud ice content at processed levels 
            qcl_latest            , &  !REAL, INTENT(OUT)        ::  QCL(nCols,LEVELS)
                                !      Liquid and frozen ice cloudoud fractions
            cloud_fraction_liquid , &  !REAL, INTENT(OUT)        ::  CFL(nCols,LEVELS)
            cloud_fraction_frozen , &  !REAL, INTENT(OUT)        ::  CFF(nCols,LEVELS)
            error_code               ) !INTEGER, INTENT(OUT) :: ERROR     !  0 if OK; 1 if bad arguments.
       !
       ! Lsarc_do1:
       DO k = 1, LEVELS
          DO i = 1, nCols
             area_cloud_fraction(i,k)=bulk_cloud_fraction(i,k)
          END DO
       END DO ! Lsarc_do1
       !
    ELSE ! L_area_cloud
       IF (L_ACF_CUSACK) THEN
          !       Vertical gradient area cloud option
          !
          inverse_level = 1.0_r8 / levels_per_level
          !
          ! Test for continuity between adjacent layers based on supersaturation
          ! (qt - qsl) / qsl : as we take differences the - qsl term drops out.
          ! DEPENDS ON: qsat_wat_mix
          CALL  qsat_wat_mix ( &
                                !      Output field
               qsl                        ,&! Real, intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
               T_latest(1,1)            ,&! Real, intent(in)  :: T(npnts) !  Temperature (K).
               p_layer_centres(1,1)     ,&! Real, intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
               nCols           ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
               l_mixing_ratio              )! logical, intent(in)  :: lq_mix
          ! .true. return qsat as a mixing ratio
          ! .false. return qsat as a specific humidity

          !        CALL qsat_wat_mix(  &
          !             qsl                    ,&
          !             T_latest(1,1,1)        ,&
          !             p_layer_centres(1,1,1) ,&
          !             nCols        ,&
          !             l_mixing_ratio          )

          DO i = 1, nCols
             qt_norm(i) =(q_latest(i,1)+qcf_latest(i,1))/qsl(i)
          END DO
          ! DEPENDS ON: qsat_wat_mix
          CALL  qsat_wat_mix ( &
                                !      Output field
               qsl                        ,&! Real, intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
               T_latest(1,2)            ,&! Real, intent(in)  :: T(npnts) !  Temperature (K).
               p_layer_centres(1,2)     ,&! Real, intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
               nCols           ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
               l_mixing_ratio              )! logical, intent(in)  :: lq_mix
          ! .true. return qsat as a mixing ratio
          ! .false. return qsat as a specific humidity

          !        CALL qsat_wat_mix(           &
          !             qsl                    ,&
          !             T_latest(1,1,2)        ,&
          !             p_layer_centres(1,1,2) ,&
          !             nCols        ,&
          !             l_mixing_ratio          )
          !
          ! Do nothing to top and bottom layers
          DO i = 1, nCols
             rhcrit_large(i,1) = rhcrit(i,1)
             rhcrit_large(i,large_levels)= rhcrit(i,LEVELS)
          END DO
          !
          ! Lsarc_do2:
          DO i = 1, nCols
             p_large(i,1) = p_layer_centres(i,1)
             T_large(i,1) = T_latest(i,1)
             q_large(i,1) = q_latest(i,1)
             qcf_large(i,1) = qcf_latest(i,1)
             !
             p_large(i,large_levels) = p_layer_centres(i,LEVELS)
             T_large(i,large_levels) = T_latest(i,LEVELS)
             q_large(i,large_levels) = q_latest(i,LEVELS)
             qcf_large(i,large_levels) = qcf_latest(i,LEVELS)
             ! Test for continuity (assumed if linked is .true.)
             qt_norm_next=(q_latest(i,2)+qcf_latest(i,2))/qsl(i)
             linked(i,1) =                                             &
                  (drat_thresh  >=  ABS(qt_norm(i) - qt_norm_next))
             qt_norm(i) = qt_norm_next
          END DO
          !
          ! Lsarc_do3:
          DO k = 2, (LEVELS - 1)
             k_index = 3 + (levels_per_level * (k-2))
             !
             ! DEPENDS ON: qsat_wat_mix
             CALL  qsat_wat_mix ( &
                                !      Output field
                  QSL                        ,&! Real, intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
                  T_latest(1,(k+1))          ,&! Real, intent(in)  :: T(npnts) !  Temperature (K).
                  p_layer_centres(1,(k+1))   ,&! Real, intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
                  nCols           ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
                  l_mixing_ratio              )! logical, intent(in)  :: lq_mix
             ! .true. return qsat as a mixing ratio
             ! .false. return qsat as a specific humidity

             !          CALL qsat_wat_mix(          &
             !           qsl                       ,&
             !           T_latest(1,1,(k+1))       ,&
             !           p_layer_centres(1,1,(k+1)),&
             !           nCols           ,&
             !           l_mixing_ratio             )
             !
             ! Select associated rhcrit values
             DO i = 1, nCols
                rhcrit_large(i,k_index-1) = rhcrit(i,k)
                rhcrit_large(i,k_index)   = rhcrit(i,k)
                rhcrit_large(i,k_index+1) = rhcrit(i,k)
             END DO
             !
             ! Lsarc_do3ij:
             DO i = 1, nCols
                ! Test for continuity (assumed if linked = .true.)
                qt_norm_next=(q_latest(i,(k+1))+qcf_latest(i,(k+1)))  &
                     / qsl(i)
                linked(i,k) =                                           &
                     (drat_thresh  >=  ABS(qt_norm(i) - qt_norm_next))
                qt_norm(i) = qt_norm_next
                !
                ! Select interpolated pressure levels
                delta_p = (p_layer_boundaries(i,(k-1)) -                &
                     p_layer_boundaries(i,k))    * inverse_level
                IF (p_layer_centres(i,k)  >=                            &
                     (p_layer_boundaries(i,k) + delta_p)) THEN
                   p_large(i,k_index) = p_layer_centres(i,k)
                ELSE
                   p_large(i,k_index) = 0.5_r8*(p_layer_boundaries(i,k) + &
                        p_layer_boundaries(i,(k-1)))
                ENDIF
                p_large(i,(k_index-1)) = p_large(i,k_index)+delta_p
                p_large(i,(k_index+1)) = p_large(i,k_index)-delta_p
                !
                ! Select variable values at layer centres
                T_large(i,k_index) = T_latest(i,k)
                q_large(i,k_index) = q_latest(i,k)
                qcf_large(i,k_index) = qcf_latest(i,k)
                !
                ! Calculate increment in variable values, pressure interpolation
                ! NB: Using X_large(i,(k_index+1)) as store for X increments
                ! Lsarc_if2:
                IF ( linked(i,(k-1)) ) THEN
                   IF ( linked(i,k) ) THEN
                      !               Interpolate from level k-1 to k+1
                      stretcher = delta_p /                                 &
                           (p_layer_centres(i,k-1)-p_layer_centres(i,k+1))
                      T_large(i,(k_index+1)) = stretcher *                &
                           (T_latest(i,(k+1)) - T_latest(i,(k-1)))
                      q_large(i,(k_index+1)) = stretcher *                &
                           (q_latest(i,(k+1)) - q_latest(i,(k-1)))
                      qcf_large(i,(k_index+1)) = stretcher *              &
                           (qcf_latest(i,(k+1)) - qcf_latest(i,(k-1)))
                   ELSE
                      !               Interpolate from level k-1 to k
                      stretcher = delta_p /                                 &
                           (p_layer_centres(i,k-1) - p_large(i,k_index))
                      T_large(i,(k_index+1)) = stretcher *                &
                           (T_large(i,k_index) - T_latest(i,(k-1)))
                      q_large(i,(k_index+1)) = stretcher *                &
                           (q_large(i,k_index) - q_latest(i,(k-1)))
                      qcf_large(i,(k_index+1)) = stretcher *              &
                           (qcf_large(i,k_index) - qcf_latest(i,(k-1)))
                   ENDIF
                   !
                ELSE
                   IF ( linked(i,k) ) THEN
                      !               Interpolate from level k to k+1
                      stretcher = delta_p /                                 &
                           (p_large(i,k_index) - p_layer_centres(i,k+1))
                      T_large(i,(k_index+1)) = stretcher *               &
                           (T_latest(i,(k+1)) - T_large(i,k_index))
                      q_large(i,(k_index+1)) = stretcher *               &
                           (q_latest(i,(k+1)) - q_large(i,k_index))
                      qcf_large(i,(k_index+1)) = stretcher *             &
                           (qcf_latest(i,(k+1)) - qcf_large(i,k_index))
                   ELSE
                      !               No interpolation, freeze at level k
                      T_large(i,(k_index+1)) = 0.0_r8
                      q_large(i,(k_index+1)) = 0.0_r8
                      qcf_large(i,(k_index+1)) = 0.0_r8
                   ENDIF
                   !
                END IF ! Lsarc_if2
                !
                ! Protect against q or qcf going negative (T would imply blow-up anyway)
                IF (q_large(i,k_index)  <                               &
                     (ABS(q_large(i,(k_index+1)))+tol_test))      &
                     q_large(i,(k_index+1)) = 0.0_r8
                IF (qcf_large(i,k_index)  <                             &
                     (ABS(qcf_large(i,(k_index+1)))+tol_test))  &
                     qcf_large(i,(k_index+1)) = 0.0_r8
                !
                ! Select variable values at level below layer centre
                T_large(i,(k_index-1)) = T_large(i,k_index) -         &
                     T_large(i,(k_index+1))
                q_large(i,(k_index-1)) = q_large(i,k_index) -         &
                     q_large(i,(k_index+1))
                qcf_large(i,(k_index-1))=qcf_large(i,k_index)-        &
                     qcf_large(i,(k_index+1))
                !
                ! Select variable values at level above layer centre
                ! NB: CEASE using X_large(i,(k_index+1)) as store for X increments.
                T_large(i,(k_index+1)) = T_large(i,(k_index+1)) +     &
                     T_large(i,k_index)
                q_large(i,(k_index+1)) = q_large(i,(k_index+1)) +     &
                     q_large(i,k_index)
                qcf_large(i,(k_index+1))=qcf_large(i,(k_index+1)) +   &
                     qcf_large(i,k_index)
             END DO
          END DO ! Lsarc_do3ij
          !
          !
          !
          ! Create an array of NTML values adjusted to the large levels

          DO i = 1, nCols
             ntml_large(i) = 3+(levels_per_level*(ntml(i)-1))
          END DO

          ! DEPENDS ON: ls_cld
          CALL ls_cld( &
                                !      Pressure related fields
               p_large                     , & !REAL    , INTENT(IN) :: p_theta_levels(nCols,levels)! pressure at all points (Pa).
               rhcrit_large                , & !REAL    , INTENT(IN) :: RHCRIT(nCols,LEVELS)
                                !      Array dimensions
               large_levels                , & !INTEGER , INTENT(IN) :: LEVELS!         No. of levels being processed.
               bl_levels                   , & !INTEGER , INTENT(IN) :: BL_LEVELS!            No. of boundary layer levels
               nCols                  , & !INTEGER , INTENT(IN) :: nCols
               cloud_fraction_method       , & !Integer , Intent(IN) :: cloud_fraction_method ! Method for calculating
               overlap_ice_liquid          , & !Real    , Intent(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
               ice_fraction_method         , & !Integer , Intent(IN) :: ice_fraction_method   ! Method for calculating ice 
               ctt_weight                  , & !Real    , Intent(IN) :: ctt_weight! Weighting of cloud top temperature
               t_weight                    , & !Real    , Intent(IN) :: t_weight  ! Weighting of local temperature
               qsat_fixed                  , & !Real    , Intent(IN) :: qsat_fixed ! Fixed value of saturation humidity
               sub_cld                     , & !Real    , Intent(IN) :: sub_cld                  ! Scaling parameter
                                !      From convection diagnosis (only uss (only used if a05_4a)
               ntml_large                  , & !INTEGER, INTENT(IN        ) :: NTML(nCols)     ! IN Height of diagnosed BL top
               cumulus                     , & !LOGICAL, INTENT(IN        ) ::  CUMULUS(nCols)  ! IN Logical indicator of convection
               L_conv4a                    , & !LOGICAL, INTENT(IN        ) ::  L_conv4a! IN true if using 4A convection scheme
               L_eacf                      , & !LOGICAL, INTENT(IN        ) ::  L_eacf! IN true if using empirically adjusted cloud
               l_mixing_ratio              , & !LOGICAL, INTENT(IN        ) ::  L_mixing_ratio  ! IN true if using mixing ratios
                                !      Prognostic Fields
               T_large                     , & !REAL   , INTENT(INOUT) ::   T(nCols,LEVELS)
               cloud_fraction_large        , & !REAL, INTENT(OUT)   ::  CF
               q_large                     , & !REAL   , INTENT(INOUT) ::  Q(nCols,LEVELS)
               qcf_large                   , & !REAL    , INTENT(IN) :: QCF(nCols,LEVELS)!Cloud ice content at processed levels 
               qcl_large                   , & !REAL, INTENT(OUT)   ::  QCL(nCols,LEVELS)
                                !      Liquid and frozen ice cloud fractioud fractions
               cloud_fraction_liquid_large , & !REAL, INTENT(OUT)   ::  CFL(nCols,LEVELS)
               cloud_fraction_frozen_large , & !REAL, INTENT(OUT)   ::  CFF(nCols,LEVELS)
               error_code                    ) !INTEGER, INTENT(OUT) :: ERROR     !  0 if OK; 1 if bad arguments.

          ! Lsarc_do4:
          DO i = 1, nCols
             T_latest(i,1) = T_large(i,1)
             q_latest(i,1) = q_large(i,1)

             area_cloud_fraction(i,1) = cloud_fraction_large(i,1)
             bulk_cloud_fraction(i,1) = cloud_fraction_large(i,1)

             qcl_latest(i,1) = qcl_large(i,1)
             cloud_fraction_liquid(i,1) =                              &
                  cloud_fraction_liquid_large(i,1)
             cloud_fraction_frozen(i,1) =                              &
                  cloud_fraction_frozen_large(i,1)

             T_latest(i,LEVELS) = T_large(i,large_levels)
             q_latest(i,LEVELS) = q_large(i,large_levels)

             area_cloud_fraction(i,LEVELS) =                 &
                  cloud_fraction_large(i,large_levels)
             bulk_cloud_fraction(i,LEVELS) =                 &
                  cloud_fraction_large(i,large_levels)

             qcl_latest(i,LEVELS)=qcl_large(i,large_levels)
             cloud_fraction_liquid(i,LEVELS) =               &
                  cloud_fraction_liquid_large(i,large_levels)
             cloud_fraction_frozen(i,LEVELS) =               &
                  cloud_fraction_frozen_large(i,large_levels)
          END DO
          !
          ! Output variables for remaining layers
          ! Lsarc_do5:
          DO k = 2, (LEVELS - 1)
             k_index = 3 + (levels_per_level * (k-2))
             ! Lsarc_do5ij:
             DO i = 1, nCols
                ! Area cloud fraction is maximum of sub-layer cloud fractions
                area_cloud_fraction(i,k) =                              &
                     MAX( cloud_fraction_large(i,k_index),                   &
                     (MAX(cloud_fraction_large(i,(k_index+1)),          &
                     cloud_fraction_large(i,(k_index-1)))) )
                !
                ! Bulk cloud fraction is mean of sub-layer cloud fractions : strictly
                ! this is a pressure weighted mean being used to approximate a volume
                ! mean. Over the depth of a layer the difference should not be large.
                bulk_cloud_fraction(i,k) = inverse_level *              &
                     ( cloud_fraction_large(i,(k_index-1)) +              &
                     cloud_fraction_large(i, k_index)    +              &
                     cloud_fraction_large(i,(k_index+1)) )
                !
                ! The pressure weighted mean of qcf is the input qcf: do not update.
                !
                ! Qcl is the pressure weighted mean of qcl from each sub-layer.
                qcl_latest(i,k) = inverse_level *                       &
                     ( qcl_large(i,(k_index-1)) +                              &
                     qcl_large(i,k_index) + qcl_large(i,(k_index+1)) )
                !
                ! Liq. cloud fraction is mean of sub-layer cloud fractions : strictly
                ! this is a pressure weighted mean being used to approximate a volume
                ! mean. Over the depth of a layer the difference should not be large.
                cloud_fraction_liquid(i,k) = inverse_level *            &
                     ( cloud_fraction_liquid_large(i,(k_index-1)) +            &
                     cloud_fraction_liquid_large(i,k_index)     +            &
                     cloud_fraction_liquid_large(i,(k_index+1)) )
                !
                ! Froz cloud fraction is mean of sub-layer cloud fractions : strictly
                ! this is a pressure weighted mean being used to approximate a volume
                ! mean. Over the depth of a layer the difference should not be large.
                cloud_fraction_frozen(i,k) = inverse_level *            &
                     ( cloud_fraction_frozen_large(i,(k_index-1)) +            &
                     cloud_fraction_frozen_large(i,k_index)     +            &
                     cloud_fraction_frozen_large(i,(k_index+1)) )
                !
                ! Transform q_latest from qT(vapour + liquid) to specific humidity.
                ! Transform T_latest from TL(vapour + liquid) to temperature.
                q_latest(i,k) = q_latest(i,k) - qcl_latest(i,k)
                T_latest(i,k) = T_latest(i,k) +                       &
                     (qcl_latest(i,k) * LCRCP)
             END DO
             !
          END DO ! Lsarc_do5
          !
       ELSE IF (L_ACF_Brooks) THEN  ! L_ACF_Cusack or Brooks
          !
          !       As before, update variables that would have been updated
          !       without area cloud fraction on
          !
          ! DEPENDS ON: ls_cld
          CALL ls_cld(&
                                !      Pressure related fields
               p_layer_centres(1,1)    ,&!REAL    , INTENT(IN) :: p_theta_levels(nCols,levels)! pressure at all points (Pa).
               rhcrit                    ,&!REAL    , INTENT(IN) :: RHCRIT(nCols,LEVELS)
                                !      Array dimensions
               LEVELS          ,&!INTEGER , INTENT(IN) :: LEVELS!            No. of levels being processed.
               bl_levels                 ,&!INTEGER , INTENT(IN) :: BL_LEVELS!       No. of boundary layer levels
               nCols                ,&!INTEGER , INTENT(IN) :: nCols
               cloud_fraction_method     ,&!Integer , Intent(IN) :: cloud_fraction_method ! Method for calculating
               overlap_ice_liquid        ,&!Real    , Intent(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
               ice_fraction_method       ,&!Integer , Intent(IN) :: ice_fraction_method   ! Method for calculating ice 
               ctt_weight                ,&!Real    , Intent(IN) :: ctt_weight! Weighting of cloud top temperature
               t_weight                  ,&!Real    , Intent(IN) :: t_weight  ! Weighting of local temperature
               qsat_fixed                ,&!Real    , Intent(IN) :: qsat_fixed ! Fixed value of saturation humidity
               sub_cld                   ,&!Real    , Intent(IN) :: sub_cld                     ! Scaling parameter
                                !      From convection diagnosis (ons (only used if a05_4a)
               ntml                      ,&!INTEGER, INTENT(IN   ) :: NTML(nCols)        ! IN Height of diagnosed BL top
               cumulus                   ,&!LOGICAL, INTENT(IN   ) ::  CUMULUS(nCols)  ! IN Logical indicator of convection
               L_conv4a                  ,&!LOGICAL, INTENT(IN   ) ::  L_conv4a! IN true if using 4A convection scheme
               L_eacf                    ,&!LOGICAL, INTENT(IN   ) ::  L_eacf! IN true if using empirically adusted cloud
               L_mixing_ratio            ,&!LOGICAL, INTENT(IN   ) ::  L_mixing_ratio  ! IN true if using mixing ratios
                                !      Prognostic Fields
               T_latest                  ,&!REAL   , INTENT(INOUT) ::   T(nCols,LEVELS)
               bulk_cloud_fraction       ,&!REAL, INTENT(OUT)   ::  CF
               q_latest                  ,&!REAL   , INTENT(INOUT) ::  Q(nCols,LEVELS)
               qcf_latest                ,&!REAL    , INTENT(IN) :: QCF(nCols,LEVELS)!Cloud ice content at processed levels 
               qcl_latest                ,&!REAL, INTENT(OUT)   ::  QCL(nCols,LEVELS)
                                !      Liquid and frozen ice cloud oud fractions
               cloud_fraction_liquid     ,&!REAL, INTENT(OUT)   ::  CFL(nCols,LEVELS)
               cloud_fraction_frozen     ,&!REAL, INTENT(OUT)   ::  CFF(nCols,LEVELS)
               error_code                 )!INTEGER, INTENT(OUT) :: ERROR        !  0 if OK; 1 if bad arguments.
          !
          !      Calculate the area cloud fraction

          ! DEPENDS ON: ls_acf_brooks
          CALL LS_ACF_Brooks ( &
                                ! model dimensions
               nCols            , &!INTEGER, INTENT(IN        ) ::  nCols
               LEVELS      , &!INTEGER, INTENT(IN        ) :: LEVELS
                                ! in coordinate information
               r_theta_levels        , &!Real ,INTENT(IN)::  r_theta_levels (nCols,0:LEVELS)
               delta_lambda          , &!Real ,INTENT(IN):: delta_lambda ! EW (x) grid spacing in radians
               delta_phi             , &!Real ,INTENT(IN):: delta_phi    ! NS (y) grid spacing in radians
                                ! trig arrays
               FV_cos_theta_latitude , &!Real ,INTENT(IN):: FV_cos_theta_latitude (1:nCols)! Finite volume cos(lat)
                                ! in data fields
               bulk_cloud_fraction   , &!Real ,INTENT(IN)::  bulk_cloud_fraction (nCols, LEVELS)
               cloud_fraction_liquid , &!Real  ,INTENT(IN):: cloud_fraction_liquid (nCols, LEVELS)
               cloud_fraction_frozen , &!Real  ,INTENT(IN):: cloud_fraction_frozen (nCols, LEVELS)
                                ! in logical control
               cumulus               , &!Logical, intent(in):: cumulus(nCols)
                                ! out data fields
               area_cloud_fraction     )!Real::   area_cloud_fraction (nCols, LEVELS)

          !
       END IF ! L_ACF_Brooks
    END IF ! Lsarc_if1  L_area_cloud
    !
    RETURN
  END SUBROUTINE LS_ARCLD
  ! ======================================================================




  !+    Brooks area cloud fraction parametrisation
  !
  !     *********************** COPYRIGHT *************************
  !     Crown Copyright <year>, The Met. Office. All rights reserved.
  !     *********************** COPYRIGHT *************************
  !
  !     Subroutine Interface:

  SUBROUTINE LS_ACF_Brooks ( &
                                ! model dimensions
       nCols                  ,&!INTEGER, INTENT(IN   ) ::  nCols
       LEVELS                 ,&!INTEGER, INTENT(IN   ) :: LEVELS
                                ! in coordinate information
       r_theta_levels         ,&!Real ,INTENT(IN)::  r_theta_levels (nCols,0:LEVELS)
       delta_lambda           ,&!Real ,INTENT(IN):: delta_lambda ! EW (x) grid spacing in radians
       delta_phi              ,&!Real ,INTENT(IN):: delta_phi    ! NS (y) grid spacing in radians
                                ! trig arrays
       FV_cos_theta_latitude  ,&!Real ,INTENT(IN):: FV_cos_theta_latitude (1:nCols)! Finite volume cos(lat)
                                ! in data fields
       bulk_cloud_fraction    ,&!Real ,INTENT(IN)::  bulk_cloud_fraction (nCols, LEVELS)
       cloud_fraction_liquid  ,&!Real  ,INTENT(IN):: cloud_fraction_liquid (nCols, LEVELS)
       cloud_fraction_frozen  ,&!Real  ,INTENT(IN):: cloud_fraction_frozen (nCols, LEVELS)
                                ! in logical control
       cumulus                ,&!Logical, intent(in):: cumulus(nCols)
                                ! out data fields
       area_cloud_fraction     )!Real::   area_cloud_fraction (nCols, LEVELS)

    !     Description:
    !       Calculates area_cloud_fraction from bulk_cloud_fraction
    !
    !     Method:
    !       The calculation is  based on the parametrisation in
    !       Brooks 2005 equations 2-3.
    !       (Brooks et al, July 2005, JAS vol 62 pp 2248-2260)
    !       The initial parametrisation uses the values in equations
    !       4,5,7 and 8 of the paper for ice and liquid cloud without
    !       wind shear.  For mixed phase clouds the maximum of the two
    !       area_cloud_fractions resulting will be used.
    !       Only area_cloud_fraction will be updated.
    !       Grid box size is needed to be known.
    !
    !     Current Code Owner: LS cloud scheme owner
    !
    !     History:
    !     Version  Date      Comment
    !
    !       6.4    19/09/05  Original code. (Amanda Kerr-Munslow, placed
    !                                        into UM by D Wilson)
    !
    !     Code Description:
    !       FORTRAN 77 with extensions recommended in the Met. Office
    !       F77 Standard.
    !
    IMPLICIT NONE

    ! Scalar arguments with INTENT(IN):

    ! Model dimensions
    INTEGER, INTENT(IN   ) ::  nCols
    ! number of points on a row
    INTEGER, INTENT(IN   ) :: LEVELS
    ! number of model levels where moisture
    ! variables are held

    ! Array Arguments with INTENT(IN)
    ! Co-ordinate arrays:
    REAL(KIND=r8) ,INTENT(IN)::  r_theta_levels (nCols,0:LEVELS)
    ! height of theta levels (from centre of earth)
    REAL(KIND=r8) ,INTENT(IN):: FV_cos_theta_latitude (1:nCols)! Finite volume cos(lat)
    REAL(KIND=r8) ,INTENT(IN):: delta_lambda ! EW (x) grid spacing in radians
    REAL(KIND=r8) ,INTENT(IN):: delta_phi    ! NS (y) grid spacing in radians

    ! Data arrays:
    REAL(KIND=r8) ,INTENT(IN)::  bulk_cloud_fraction (nCols, LEVELS)
    !       Cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)  ,INTENT(IN):: cloud_fraction_liquid (nCols, LEVELS)
    !       Liquid cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8)  ,INTENT(IN):: cloud_fraction_frozen (nCols, LEVELS)
    !       Frozen cloud fraction at processed levels (decimal fraction).

    LOGICAL, INTENT(in):: cumulus(nCols)

    ! Arguments with INTENT(OUT):
    ! Data arrays:
    REAL(KIND=r8)::   area_cloud_fraction (nCols, LEVELS)
    !       Cloud fraction at processed levels (decimal fraction).

    ! Local Parameters:

    ! Parameters for liquid clouds, from Brooks 2005 equations 7 and 8
    REAL(KIND=r8) ,PARAMETER :: power_law_gradient_liquid = 0.1635_r8 ! A
    REAL(KIND=r8) ,PARAMETER :: vert_fit_liquid = 0.6694_r8   ! alpha
    REAL(KIND=r8) ,PARAMETER :: horiz_fit_liquid = -0.1882_r8 ! beta

    ! Parameters for frozen clouds, from Brooks 2005 equations 4 and 5
    REAL(KIND=r8),PARAMETER ::  power_law_gradient_frozen = 0.0880_r8! A
    REAL(KIND=r8),PARAMETER ::  vert_fit_frozen = 0.7679_r8  ! alpha
    REAL(KIND=r8),PARAMETER :: horiz_fit_frozen= -0.2254  ! beta

    ! Local Scalars:
    ! Loop counters
    INTEGER ::   i, k

    REAL(KIND=r8) ::  symmetric_adjustment_liquid
    !    function f in eqn 7 in Brooks 2005
    REAL(KIND=r8) ::   symmetric_adjustment_frozen
    !    function f in eqn 4 in Brooks 2005
    REAL(KIND=r8) ::  horiz_scale
    !    horizontal scale size of the grid box (m)
    REAL(KIND=r8) ::   vert_scale
    !    vertical scale size of the grid box (m)


    !  Local Arrays:
    REAL(KIND=r8) ::  acf_liquid (nCols, LEVELS)
    !    area cloud fraction based on liquid parameters
    REAL(KIND=r8) ::  acf_frozen (nCols, LEVELS)
    !    area cloud fraction based on frozen parameters

    !-    End of header
    ! ----------------------------------------------------------------------

    ! ==Main Block==--------------------------------------------------------

    ! Initialise arrays and local variables to zero
    DO k = 1, LEVELS
       DO i = 1, nCols
          area_cloud_fraction(i,k) = 0.0_r8
          acf_liquid(i,k) = 0.0_r8
          acf_frozen(i,k) = 0.0_r8
       END DO
    END DO
    horiz_scale = 0.0_r8
    vert_scale = 0.0_r8
    symmetric_adjustment_liquid = 0.0_r8
    symmetric_adjustment_frozen = 0.0_r8

    DO k = 1, LEVELS
       DO i = 1, nCols
          !
          ! Test if bulk_cloud_fraction is within bounds of possibility
          !
          IF ( bulk_cloud_fraction(i,k) <= 0.0_r8 ) THEN

             area_cloud_fraction(i,k) = 0.0_r8

          ELSE IF ( bulk_cloud_fraction(i,k) >= 1.0_r8 ) THEN

             area_cloud_fraction(i,k) = 1.0_r8

          ELSE IF (cumulus(i)) THEN
             ! This is a convective point so do not apply the
             ! area cloud representation
             area_cloud_fraction(i,k) = bulk_cloud_fraction(i,k)

          ELSE

             ! Only calculate area_cloud_fraction if the bulk_cloud_fraction
             ! is between (not equal to) 0.0_r8 and 1.0_r8
             !
             ! Calculate horizontal and vertical scales.
             ! The horizontal scale is taken as the square root of the
             ! area of the grid box.
             ! The vertical scale is taken as the difference in radius
             ! from the centre of the Earth between the upper and lower
             ! boundaries of the grid box.

             horiz_scale = SQRT (                                      &
                  r_theta_levels(i,k)                   &
                  * r_theta_levels(i,k)                 &
                  * delta_lambda * delta_phi              &
                  * FV_cos_theta_latitude(i) )
             IF (k .EQ. LEVELS) THEN
                ! Assume top layer thickness is the same as the
                ! thickness of the layer below
                vert_scale =  r_theta_levels(i,k)                     &
                     - r_theta_levels(i,k-1)
             ELSE 
                vert_scale =  r_theta_levels(i,k+1)                   &
                     - r_theta_levels(i,k)
             END IF  ! k eq LEVELS

             ! Calculate the symmetric_adjustment (f).
             ! This parameter controls the extent to which the area cloud fraction
             ! is greater than the bulk cloud fraction.  If f = 0, they are equal.

             symmetric_adjustment_liquid =                             &
                  power_law_gradient_liquid                      &
                  * ( vert_scale ** vert_fit_liquid )            &
                  * ( horiz_scale ** horiz_fit_liquid )
             symmetric_adjustment_frozen =                             &
                  power_law_gradient_frozen                      &
                  * ( vert_scale ** vert_fit_frozen )            &
                  * ( horiz_scale ** horiz_fit_frozen )

             ! Calculate the area cloud fractions for liquid and frozen cloud
             ! Calculate the liquid and frozen fractions separately to
             ! allow for greatest flexibility in future choice of decisions
             ! regarding mixed phase cloud.

             acf_liquid(i,k) = 1.0_r8/                                   &
                  ( 1.0_r8 + ( EXP(-1.0_r8*symmetric_adjustment_liquid)         &
                  * ( 1.0_r8/bulk_cloud_fraction(i,k) - 1.0_r8) ) )
             acf_frozen(i,k) = 1.0_r8/                                   &
                  ( 1.0_r8 + ( EXP(-1.0_r8*symmetric_adjustment_frozen)         &
                  * ( 1.0_r8/bulk_cloud_fraction(i,k) - 1.0_r8) ) )

             ! Calculate the final area cloud fraction for each grid box
             ! Currently this is based on which there is more of, ice or liquid.

             IF ( cloud_fraction_frozen(i,k) == 0.0_r8 ) THEN
                IF ( cloud_fraction_liquid(i,k) == 0.0_r8 ) THEN

                   ! If there is no liquid or frozen cloud, there should be no area cloud
                   area_cloud_fraction(i,k) = 0.0_r8

                ELSE

                   ! If there is no frozen cloud but there is liquid cloud,
                   ! then the area cloud fraction is given by the liquid parametrisation
                   ! 0 no cloud, 1 either, 2 liq, 3 ice'
                   area_cloud_fraction(i,k) = acf_liquid(i,k)
                END IF

             ELSE ! cloud_fraction_frozen

                IF ( cloud_fraction_liquid(i,k) == 0.0_r8 ) THEN

                   ! If there is frozen cloud but there is no liquid cloud,
                   ! then the area cloud fraction is given by the frozen parametrisation
                   area_cloud_fraction(i,k) = acf_frozen(i,k)

                ELSE

                   ! If there is frozen cloud and there is liquid cloud,
                   ! then the area cloud fraction is given by the maximum of the two
                   ! parametrisations
                   area_cloud_fraction(i,k) =                          &
                        MAX( acf_liquid(i,k),acf_frozen(i,k) )

                END IF

             END IF ! cloud_fraction_frozen

          END IF ! bulk_cloud_fraction between 0.0_r8 and 1.0_r8

       END DO
    END DO

    RETURN
  END SUBROUTINE LS_ACF_Brooks
  ! ======================================================================


  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Large-scale Cloud Scheme.
  ! Subroutine Interface:
  SUBROUTINE LS_CLD( &
                                !      Pressure related fields
       p_theta_levels        ,& !REAL    , INTENT(IN) :: p_theta_levels(nCols,levels)! pressure at all points (Pa).
       RHCRIT                ,& !REAL    , INTENT(IN) :: RHCRIT(nCols,LEVELS)
                                !      Array dimensions
       LEVELS                ,& !INTEGER , INTENT(IN) :: LEVELS!       No. of levels being processed.
       BL_LEVELS             ,& !INTEGER , INTENT(IN) :: BL_LEVELS!       No. of boundary layer levels
       nCols                 ,& !INTEGER , INTENT(IN) :: nCols
       cloud_fraction_method ,& !Integer , Intent(IN) :: cloud_fraction_method ! Method for calculating
       overlap_ice_liquid    ,& !Real    , Intent(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
       ice_fraction_method   ,& !Integer , Intent(IN) :: ice_fraction_method   ! Method for calculating ice 
       ctt_weight            ,& !Real    , Intent(IN) :: ctt_weight! Weighting of cloud top temperature
       t_weight              ,& !Real    , Intent(IN) :: t_weight  ! Weighting of local temperature
       qsat_fixed            ,& !Real    , Intent(IN) :: qsat_fixed ! Fixed value of saturation humidity
       sub_cld               ,& !Real    , Intent(IN) :: sub_cld                ! Scaling parameter
                                !      From convection diagnosis (only used if a05_4a)
       ntml                  ,& !INTEGER, INTENT(IN   ) :: NTML(nCols)     ! IN Height of diagnosed BL top
       cumulus               ,& !LOGICAL, INTENT(IN   ) ::  CUMULUS(nCols)  ! IN Logical indicator of convection
       L_conv4a              ,& !LOGICAL, INTENT(IN   ) ::  L_conv4a! IN true if using 4A convection scheme
       L_eacf                ,& !LOGICAL, INTENT(IN   ) ::  L_eacf! IN true if using empirically adjusted cloud
       L_mixing_ratio        ,& !LOGICAL, INTENT(IN   ) ::  L_mixing_ratio  ! IN true if using mixing ratios
                                !      Prognostic Fields
       T                     ,& !REAL   , INTENT(INOUT) ::   T(nCols,LEVELS)
       CF                    ,& !REAL, INTENT(OUT)   ::  CF
       Q                     ,& !REAL   , INTENT(INOUT) ::  Q(nCols,LEVELS)
       QCF                   ,& !REAL    , INTENT(IN) :: QCF(nCols,LEVELS)!Cloud ice content at processed levels 
       QCL                   ,& !REAL, INTENT(OUT)   ::  QCL(nCols,LEVELS)
                                !      Liquid and frozen ice cloud fractions
       CFL                   ,& !REAL, INTENT(OUT)   ::  CFL(nCols,LEVELS)
       CFF                   ,& !REAL, INTENT(OUT)   ::  CFF(nCols,LEVELS)
       ERROR                  ) !INTEGER, INTENT(OUT) :: ERROR     !  0 if OK; 1 if bad arguments.
    !

    IMPLICIT NONE
    !
    ! Purpose:
    !   This subroutine calculates liquid and ice cloud fractional cover
    !   for use with the enhanced precipitation microphysics scheme.
    !
    ! Method:
    !   Statistical cloud scheme separates input moisture into specific
    !   humidity and cloud liquid water. Temperature calculated from liquid
    !   water temperature. Cloud fractions calculated from statistical
    !   relation between cloud fraction and cloud liquid/ice water content.
    !
    ! Current Owner of Code: A. C. Bushell
    !
    ! History:
    ! Version   Date     Comment
    !  4.4    14-11-96   Original Code (A. C. Bushell from Wilson/Ballard)
    !
    !  4.5    16-01-98   Change to estimated bs in ice cloud fraction
    !                    calculation (use QSat_Wat).  AC Bushell
    !
    !  5.1    09-12-99   Change to allow 3D RHcrit specification. AC Bushell
    !
    !  5.2    21-11-00   Changes to ice cloud fraction. Damian Wilson
    !
    !LL  5.3  24/09/01  Portability changes.    Z. Gardner
    !
    !  5.4   10-09-02    If cumulus convection is diagnosed the large
    !                    scale cloud is not allowed in the levels around
    !                    cumulus cloud-base.           Gill Martin
    !
    !  5.5   11-02-03    Fix to allow RHcrit = 1.   R.M. Forbes
    !  6.0  19/08/03  NEC SX-6 optimisation - Avoid bank conflicts by
    !                 replacing repeated access to RHCRIT(1,1,K) by
    !                 temporary scalar.  R Barnes & J-C Rioual.
    !  6.2  20/03/06  Revised SX-6 optimisation by J-C Rioual.
    !                 Lodged by M Saunby.
    !  6.2  07-11-05  Include the EACF parametrization for cloud
    !                    fraction.                             D. Wilson
    !  6.2  22-11-05  Include cloud parametrization quantities. D. Wilson
    !  6.4  19-12-06  NEC vectorisation optimisation and tidy-up comments.
    !                 P.Selwood.
    !  6.4  07-08-06  Include mixing ratio formulation. D. Wilson
    !
    ! Description of Code:
    !   FORTRAN 77  + common extensions also in Fortran90.
    !   This code is written to UMDP3 version 6 programming standards.
    !
    !   Documentation: UMDP No. 29
    !

    !  Global Variables:----------------------------------------------------
    ! History:
    ! Version  Date  Comment
    !  3.4   18/5/94 Add PP missing data indicator. J F Thomson
    !  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
    !*L------------------COMDECK C_MDI-------------------------------------
    ! PP missing data indicator (-1.0E+30_r8)
    REAL(KIND=r8), PARAMETER    :: RMDI_PP  = -1.0E+30_r8

    ! Old REAL(KIND=r8) missing data indicator (-32768.0_r8)
    REAL(KIND=r8), PARAMETER    :: RMDI_OLD = -32768.0_r8

    ! New REAL(KIND=r8) missing data indicator (-2**30)
    REAL(KIND=r8), PARAMETER    :: RMDI     = -32768.0_r8*32768.0_r8

    ! Integer missing data indicator
    INTEGER, PARAMETER :: IMDI     = -32768
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
    !
    !  History
    !
    !  6.2   17/11/05  Remove variables that are now set in UMUI. D. Wilson
    ! start C_CLDSGS


    ! Minimum ice content to perform calculations
    REAL(KIND=r8),PARAMETER :: QCFMIN=1.0E-8_r8

    ! end C_CLDSGS
    !
    !  Subroutine Arguments:------------------------------------------------
    INTEGER , INTENT(IN) :: LEVELS!       No. of levels being processed.
    INTEGER , INTENT(IN) :: BL_LEVELS!       No. of boundary layer levels
    INTEGER , INTENT(IN) :: nCols
    !
    REAL(KIND=r8)    , INTENT(IN) :: QCF(nCols,LEVELS)!Cloud ice content at processed levels (kg water per kg air).
    REAL(KIND=r8)    , INTENT(IN) :: p_theta_levels(nCols,levels)! pressure at all points (Pa).
    REAL(KIND=r8)    , INTENT(IN) :: RHCRIT(nCols,LEVELS)
    !       Critical relative humidity.  See the the paragraph incorporating
    !       eqs P292.11 to P292.14; the values need to be tuned for the give
    !       set of levels.
    !
    INTEGER , INTENT(IN) :: cloud_fraction_method ! Method for calculating
    ! total cloud fraction
    INTEGER , INTENT(IN) :: ice_fraction_method   ! Method for calculating ice cloud frac.

    REAL(KIND=r8)    , INTENT(IN) :: overlap_ice_liquid ! Overlap between ice and liquid phases
    REAL(KIND=r8)    , INTENT(IN) :: ctt_weight! Weighting of cloud top temperature
    REAL(KIND=r8)    , INTENT(IN) :: t_weight  ! Weighting of local temperature
    REAL(KIND=r8)    , INTENT(IN) :: qsat_fixed ! Fixed value of saturation humidity
    REAL(KIND=r8)    , INTENT(IN) :: sub_cld                ! Scaling parameter
    !
    INTEGER, INTENT(IN   ) :: NTML(nCols)     ! IN Height of diagnosed BL top

    LOGICAL, INTENT(IN   ) ::  L_conv4a! IN true if using 4A convection scheme
    LOGICAL, INTENT(IN   ) ::  L_eacf! IN true if using empirically adjusted cloud
    !    fraction parametrization
    LOGICAL, INTENT(IN   ) ::  L_mixing_ratio  ! IN true if using mixing ratios

    LOGICAL, INTENT(IN   ) ::  CUMULUS(nCols)  ! IN Logical indicator of convection

    REAL(KIND=r8)   , INTENT(INOUT) ::  Q(nCols,LEVELS)
    !       On input : Total water content (QW) (kg per kg air).
    !       On output: Specific humidity at processed levels
    !                   (kg water per kg air).
    REAL(KIND=r8)   , INTENT(INOUT) ::   T(nCols,LEVELS)
    !       On input : Liquid/frozen water temperature (TL) (K).
    !       On output: Temperature at processed levels (K).
    !
    REAL(KIND=r8), INTENT(OUT)   ::  CF(nCols,LEVELS) 
    !       Cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8), INTENT(OUT)   ::  QCL(nCols,LEVELS)
    !       Cloud liquid water content at processed levels (kg per kg air).
    REAL(KIND=r8), INTENT(OUT)   ::  CFL(nCols,LEVELS)
    !       Liquid cloud fraction at processed levels (decimal fraction).
    REAL(KIND=r8), INTENT(OUT)   ::  CFF(nCols,LEVELS)
    !       Frozen cloud fraction at processed levels (decimal fraction).
    !
    !     Error Status:
    INTEGER, INTENT(OUT) :: ERROR     !  0 if OK; 1 if bad arguments.
    !
    !  Local parameters and other physical constants------------------------
    REAL(KIND=r8) :: ROOTWO       ! Sqrt(2.)
    REAL(KIND=r8) :: SUBGRID      ! Subgrid parameter in ice cloud calculation
    !
    !  Local scalars--------------------------------------------------------
    !
    !  (a) Scalars effectively expanded to workspace by the Cray (using
    !      vector registers).
    REAL(KIND=r8) :: PHIQCF    ! Arc-cosine term in Cloud ice fraction calc.
    REAL(KIND=r8) :: COSQCF    ! Cosine term in Cloud ice fraction calc.
    REAL(KIND=r8) :: OVERLAP_MAX ! Maximum possible overlap
    REAL(KIND=r8) :: OVERLAP_MIN ! Minimum possible overlap
    REAL(KIND=r8) :: OVERLAP_RANDOM! Random overlap
    REAL(KIND=r8) :: TEMP0 
    REAL(KIND=r8) :: TEMP1
    REAL(KIND=r8) :: TEMP2! Temporaries for combinations of the
    REAL(KIND=r8) :: QN_IMP!! overlap parameters
    REAL(KIND=r8) :: QN_ADJ!! overlap parameters
    !
    !  (b) Others.
    INTEGER ::  K,I       ! Loop counters: K - vertical level index.
    !                                        I,J - horizontal field indices.
    !
    INTEGER :: QC_POINTS  ! No. points with non-zero cloud
    INTEGER :: MULTRHC   ! Zero if (nCols*nCols) le 1, else 1
    !
    !  Local dynamic arrays-------------------------------------------------
    !    6 blocks of REAL(KIND=r8) workspace are required.
    REAL(KIND=r8) :: QCFRBS(nCols)! qCF / bs
    REAL(KIND=r8) :: QSL(nCols)!Saturated specific humidity for temp TL or T.
    REAL(KIND=r8) :: QSL_CTT(nCols)!Saturated specific humidity wrt liquid at cloud top temperature
    REAL(KIND=r8) :: QN(nCols)!Cloud water normalised with BS.
    REAL(KIND=r8) :: GRID_QC(nCols,LEVELS)!Gridbox mean saturation excess at processed levels
    !        (kg per kg air). Set to RMDI when cloud is absent.
    REAL(KIND=r8) :: BS(nCols,LEVELS)!Maximum moisture fluctuation /6*sigma at processed levels
    !        (kg per kg air). Set to RMDI when cloud is absent.
    REAL(KIND=r8) :: CTT(nCols)!Ice cloud top temperature (K) - as coded it is REAL(KIND=r8)ly TL
    LOGICAL :: LQC(nCols)      ! True for points with non-zero cloud
    INTEGER :: INDEX(nCols,1)  ! Index for points with non-zero cloud
    INTEGER :: LLWIC(nCols)!       Last Level With Ice Cloud
    REAL(KIND=r8)    :: RHCRITX              ! scalar copy of RHCRIT(I,J,K)
    !
    !  External subroutine calls: ------------------------------------------
    !      EXTERNAL QSAT,QSAT_WAT,LS_CLD_C
    !- End of Header
    ! ----------------------------------------------------------------------
    !  Check input arguments for potential over-writing problems.
    ! ----------------------------------------------------------------------
    ERROR=0
    !
    IF ( (nCols * nCols)  >   1) THEN
       MULTRHC = 1
    ELSE
       MULTRHC = 0
    END IF
    !
    ! ==Main Block==--------------------------------------------------------
    ! Subroutine structure :
    ! Loop round levels to be processed.
    ! ----------------------------------------------------------------------
    ! Initialize cloud-top-temperature and last-level-with-ice-cloud arrays
    !CDIR COLLAPSE
    DO i=1,nCols
       CTT(I)=0.0_r8
       LLWIC(I)=0
    END DO
    ! Levels_do1:
    DO K=LEVELS,1,-1
       !
       ! ----------------------------------------------------------------------
       ! 1. Calculate QSAT at liquid/ice water temperature, TL, and initialize
       !    cloud water, sub-grid distribution and fraction arrays.
       !    This requires a preliminary calculation of the pressure.
       !    NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
       ! ----------------------------------------------------------------------
       ! Rows_do1:
       !CDIR COLLAPSE
       ! nCols_do1:
       DO I=1,nCols
          QCL(I,K) = 0.0_r8
          CFL(I,K) = 0.0_r8
          GRID_QC(I,K) = RMDI
          BS(I,K) = RMDI
       END DO ! nCols_do1
       !
       ! DEPENDS ON: qsat_wat_mix
       CALL  qsat_wat_mix ( &
                                !      Output field
            QSL                   ,&! REAL(KIND=r8), intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
            T             (1,K) ,&! REAL(KIND=r8), intent(in)  :: T(npnts) !  Temperature (K).
            P_theta_levels(1,k) ,&! REAL(KIND=r8), intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
            nCols      ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
            l_mixing_ratio         )! logical, intent(in)  :: lq_mix
       ! .true. return qsat as a mixing ratio
       ! .false. return qsat as a specific humidity

       !
       ! nCols_do2:
       DO I=1,nCols
          IF (MULTRHC==1) THEN
             RHCRITX = RHCRIT(I,K)
          ELSE
             RHCRITX = RHCRIT(1,K)
          END IF

          ! Omit CUMULUS points below (and including) NTML+1

          IF ( .NOT. L_conv4a .OR. (L_conv4a .AND.                     &
               (.NOT. CUMULUS(I) .OR. ( CUMULUS(I)                   &
               .AND. (K  >   NTML(I)+1) ))) ) THEN

             ! Rhcrit_if:
             IF (RHCRITX  <   1.0_r8) THEN
                ! ----------------------------------------------------------------------
                ! 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
                !    if RHcrit is less than 1
                ! ----------------------------------------------------------------------
                !
                QN(I) = (Q(I,K) / QSL(I) - 1.0_r8) /                    &
                     (1.0_r8 - RHCRITX)
                !
                ! ----------------------------------------------------------------------
                ! 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
                !    where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
                ! ----------------------------------------------------------------------
                !
                LQC(I) = (QN(I)  >   -1.)
             ELSE
                ! ----------------------------------------------------------------------
                ! 2.a Calculate QN = QW - QSL if RHcrit equals 1
                ! ----------------------------------------------------------------------
                !
                QN(I) = Q(I,K) - QSL(I)
                !
                ! ----------------------------------------------------------------------
                ! 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
                !     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
                ! ----------------------------------------------------------------------
                !
                LQC(I) = (QN(I)  >   0.0_r8)
             END IF ! Rhcrit_if
          ELSEIF (L_conv4a) THEN
             LQC(I) = .FALSE.
          END IF  ! Test on CUMULUS and NTML for a05_4a only
       END DO ! nCols_do2
       !
       ! ----------------------------------------------------------------------
       ! 4. Form index of points where non-zero liquid cloud fraction
       ! ----------------------------------------------------------------------
       !
       QC_POINTS=0
       ! nCols_do3:
       DO I=1,nCols
          IF (LQC(I)) THEN
             QC_POINTS = QC_POINTS + 1
             INDEX(QC_POINTS,1) = I
          END IF
       END DO ! nCols_do3
       !
       ! ----------------------------------------------------------------------
       ! 5. Call LS_CLD_C to calculate cloud water content, specific humidity,
       !                  water cloud fraction and determine temperature.
       ! ----------------------------------------------------------------------
       ! Qc_points_if:
       IF (QC_POINTS  >   0) THEN
          ! DEPENDS ON: ls_cld_c
          CALL LS_CLD_C( &
               P_theta_levels  (1,K) , & !REAL(KIND=r8), INTENT(IN) :: P_F(nCols)!  pressure (Pa).
               RHCRIT          (1,K) , & !REAL(KIND=r8), INTENT(IN) :: RHCRIT(nCols)! Critical relative humidity.
               QSL                   , & !REAL(KIND=r8), INTENT(IN) :: QSL_F(nCols)!    saturated humidity at 
               QN                    , & !REAL(KIND=r8), INTENT(IN) :: QN_F(nCols)!   Normalised super/subsaturation
               Q               (1,K) , & !REAL(KIND=r8), INTENT(INOUT) :: Q_F(nCols)!       On input : Vapour + liquid water content (QW)
               T               (1,K) , & !REAL(KIND=r8), INTENT(INOUT) :: T_F(nCols)!       On input : Liquid water temperature (TL)
               QCL             (1,K) , & !REAL(KIND=r8), INTENT(OUT  ) :: QCL_F(nCols)!       Cloud liquid water content 
               CFL             (1,K) , & !REAL(KIND=r8), INTENT(OUT  ) :: CF_F(nCols)!       Liquid cloud fraction at processed levels.
               GRID_QC         (1,K) , & !REAL(KIND=r8), INTENT(OUT  ) :: GRID_QC_F(nCols)!   Super/subsaturation on processed levels
               BS              (1,K) , & !REAL(KIND=r8), INTENT(OUT  ) :: BS_F(nCols)!       Value of bs at processed levels
               INDEX                 , & !INTEGER , INTENT(IN) :: INDEX(nCols,1)! INDEX for  points with non-zero 
               QC_POINTS             , & !INTEGER , INTENT(IN) :: POINTS !       No. of gridpoints with non-zero cloud
               nCols                 , & !INTEGER , INTENT(IN) :: nCols !       No. of gridpoints being processed.
               BL_LEVELS             , & !INTEGER , INTENT(IN) :: BL_LEVELS! No. of boundary layer levels
               K                     , & !INTEGER , INTENT(IN) :: K !  Level no.
               L_eacf                , & !LOGICAL , INTENT(IN) :: L_eacf  !  Use empirically adjusted cloud fraction
               l_mixing_ratio           )!LOGICAL , INTENT(IN) :: L_mixing_ratio   !  Use mixing ratio formulation
       END IF ! Qc_points_if
       !
       ! ----------------------------------------------------------------------
       ! 6. Calculate cloud fractions for ice clouds.
       !    THIS IS STILL HIGHLY EXPERIMENTAL.
       !    Begin by calculating Qsat_wat(T,P), at Temp. T, for estimate of bs.
       ! ----------------------------------------------------------------------
       DO i=1,nCols
          ! Check for last level with cloud and update cloud top temperature
          IF (LLWIC(I)  /=  K+1) THEN
             CTT(I)=T(I,K)
          END IF
          ! Check for significant ice content and update last level with ice cloud
          IF (QCF(I,K)  >   QCFMIN) THEN
             LLWIC(I)=K
          ENDIF
       END DO
       ! DEPENDS ON: qsat_wat_mix
       CALL  qsat_wat_mix ( &
                                !      Output field
            QSL                   ,&! REAL(KIND=r8), intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
            T             (1,K) ,&! REAL(KIND=r8), intent(in)  :: T(npnts) !  Temperature (K).
            P_theta_levels(1,k) ,&! REAL(KIND=r8), intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
            nCols              ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
            l_mixing_ratio         )! logical, intent(in)  :: lq_mix
       ! .true. return qsat as a mixing ratio
       ! .false. return qsat as a specific humidity

       !        CALL qsat_wat_mix(QSL,T(1,1,K),P_theta_levels(1,1,k),           &
       !                      nCols*rows,l_mixing_ratio)
       ROOTWO = SQRT(2.0_r8)
       IF (ICE_FRACTION_METHOD  ==  2) THEN
          ! Use cloud top temperature and a fixed qsat to give QCFRBS
          ! DEPENDS ON: qsat_wat_mix
          CALL  qsat_wat_mix ( &
                                !      Output field
               QSL_CTT               ,&! REAL(KIND=r8), intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
               CTT                   ,&! REAL(KIND=r8), intent(in)  :: T(npnts) !  Temperature (K).
               p_theta_levels(1,k) ,&! REAL(KIND=r8), intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
               nCols      ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
               l_mixing_ratio         )! logical, intent(in)  :: lq_mix
          ! .true. return qsat as a mixing ratio
          ! .false. return qsat as a specific humidity

          !          CALL qsat_wat_mix(QSL_CTT,CTT,p_theta_levels(1,1,k),          &
          !                        nCols*rows,l_mixing_ratio)
          SUBGRID = SUB_CLD ** (1.0_r8-T_WEIGHT)                           &
               / QSAT_FIXED ** (1.0_r8-T_WEIGHT-CTT_WEIGHT)
       END IF ! ice_fraction_method eq 2
       !
       ! nCols_do4:
       DO I=1,nCols
          IF ( MULTRHC== 1) THEN
             RHCRITX = RHCRIT(I,K)
          ELSE
             RHCRITX = RHCRIT(1,K)
          ENDIF
          ! ----------------------------------------------------------------------
          ! 6a Calculate qCF/bs.
          ! ----------------------------------------------------------------------
          ! Rhcrit_if2:
          IF (RHCRITX  <   1.0_r8) THEN
             !
             IF (ICE_FRACTION_METHOD  ==  1) THEN
                QCFRBS(i)=  QCF(I,K) / ((1.0_r8-RHCRITX) * QSL(I))
             ELSEIF (ICE_FRACTION_METHOD  ==  2) THEN
                QCFRBS(i) = SUBGRID * QCF(I,K) / ((1.0_r8-RHCRITX)          &
                     * QSL_CTT(I)**CTT_WEIGHT*QSL(I)**T_WEIGHT)
             ELSE
                ! No ice cloud fraction method defined
             END IF ! ice_fraction_method
             !
             ! ----------------------------------------------------------------------
             ! 6b Calculate frozen cloud fraction from frozen cloud water content.
             ! ----------------------------------------------------------------------
             IF (QCFRBS(i)  <=  0.0_r8) THEN
                CFF(I,K) = 0.0_r8
             ELSEIF (0.0_r8  <   QCFRBS(i)  .AND. (6.0_r8 * QCFRBS(i))  <=  1.0_r8) THEN
                CFF(I,K) = 0.5_r8 * ((6.0_r8 * QCFRBS(i))**(2.0_r8/3.0_r8))
             ELSEIF (1.0_r8  <   (6.0_r8*QCFRBS(i)) .AND. QCFRBS(i)  <   1.0_r8) THEN
                PHIQCF = ACOS(ROOTWO * 0.75_r8 * (1.0_r8 - QCFRBS(i)))
                COSQCF = COS((PHIQCF + (4.0_r8 * PI)) / 3.0_r8)
                CFF(I,K) = 1.0_r8 - (4.0_r8 * COSQCF * COSQCF)
             ELSEIF (QCFRBS(i)  >=  1.0_r8) THEN
                CFF(I,K) = 1.0_r8
             END IF
             IF (L_eacf) THEN  ! Empirically adjusted cloud fraction
                ! Back out QN
                IF (0.0_r8 <  QCFRBS(i)  .AND. (6.0_r8 * QCFRBS(i))  <=  1.0_r8) THEN
                   QN_IMP=SQRT(2.0_r8*CFF(I,K))-1.0_r8
                ELSEIF (1.0_r8  <   (6.0_r8*QCFRBS(i)) .AND. QCFRBS(i) <  1.0_r8) THEN
                   QN_IMP=1.0_r8-SQRT((1.0_r8-CFF(I,K))*2.0_r8)
                ELSE
                   QN_IMP = 1.0_r8
                ENDIF

                ! Modify QN with EACF relationship
                IF (K >  BL_LEVELS) THEN
                   QN_ADJ=(QN_IMP+0.0955_r8)/(1.0_r8-0.0955_r8)
                ELSE
                   QN_ADJ=(QN_IMP+0.184_r8)/(1.0_r8-0.184_r8)
                ENDIF

                ! Recalculate ice cloud fraction with modified QN
                IF (QCFRBS(i)  <=  0.0_r8) THEN
                   CFF(I,K) = 0.0_r8
                ELSEIF (QN_ADJ  <=  0.0_r8) THEN
                   CFF(I,K) = 0.5_r8 * (1.0_r8 + QN_ADJ) * (1.0_r8 + QN_ADJ)
                ELSEIF (QN_ADJ  <   1.0_r8) THEN
                   CFF(I,K) = 1.0_r8 - 0.5_r8 * (1.0_r8-QN_ADJ) * (1.0_r8-QN_ADJ)
                ELSE
                   CFF(I,K) = 1.0_r8
                ENDIF

             END IF  ! L_eacf
             !

          ELSE ! RHcrit = 1, set cloud fraction to 1 or 0
             !
             IF (QCF(I,K)  >   0.0_r8) THEN
                CFF(I,K) = 1.0_r8
             ELSE
                CFF(I,K) = 0.0_r8
             ENDIF
             !
          ENDIF
       END DO ! nCols_do4

       ! ----------------------------------------------------------------------
       ! 6c Calculate combined cloud fraction.
       ! ----------------------------------------------------------------------
       ! nCols_do5:
       DO I=1,nCols
          IF (CLOUD_FRACTION_METHOD  ==  1) THEN
             !           Use minimum overlap condition
             CF(I,K) = MIN(CFL(I,K)+CFF(I,K), 1.0_r8)
          ELSEIF (CLOUD_FRACTION_METHOD  ==  2) THEN
             ! Calculate possible overlaps between ice and liquid in THIS layer
             OVERLAP_MAX=MIN(CFL(I,K),CFF(I,K))
             OVERLAP_MIN=MAX(CFL(I,K)+CFF(I,K)-1.0_r8,0.0_r8)
             OVERLAP_RANDOM=CFL(I,K)*CFF(I,K)
             ! Now use the specified degree of overlap to calculate the total
             ! cloud fraction (= cfice + cfliq - overlap)
             TEMP0=OVERLAP_RANDOM
             TEMP1=0.5_r8*(OVERLAP_MAX-OVERLAP_MIN)
             TEMP2=0.5_r8*(OVERLAP_MAX+OVERLAP_MIN)-OVERLAP_RANDOM
             CF(I,K)=CFL(I,K)+CFF(I,K)                                &
                  -(TEMP0+TEMP1*OVERLAP_ICE_LIQUID                    &
                  +TEMP2*OVERLAP_ICE_LIQUID*OVERLAP_ICE_LIQUID)
             ! Check that the overlap wasnt negative
             CF(I,K)=MIN(CF(I,K),CFL(I,K)+CFF(I,K))
          ELSE
             ! No total cloud fraction method defined
          END IF ! cloud_fraction_method
          !
       END DO ! nCols_do5
       !
    END DO ! Levels_do
    !
    RETURN
  END SUBROUTINE LS_CLD


  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Large-scale Cloud Scheme.
  ! Subroutine Interface:
  ! ======================================================================
  !
  !+ Large-scale Cloud Scheme Compression routine (Cloud points only).
  ! Subroutine Interface:
  SUBROUTINE LS_CLD_C( &
       P_F           , &!REAL(KIND=r8)    , INTENT(IN) :: P_F(nCols)!       pressure (Pa).
       RHCRIT        , &!REAL(KIND=r8)    , INTENT(IN) :: RHCRIT(nCols)! Critical relative humidity.
       QSL_F         , &!REAL(KIND=r8)    , INTENT(IN) :: QSL_F(nCols)!       saturated humidity at 
       QN_F          , &!REAL(KIND=r8)    , INTENT(IN) :: QN_F(nCols)!       Normalised super/subsaturation
       Q_F           , &!REAL(KIND=r8)    , INTENT(INOUT) :: Q_F(nCols)!       On input : Vapour + liquid water content (QW)
       T_F           , &!REAL(KIND=r8)    , INTENT(INOUT) :: T_F(nCols)!       On input : Liquid water temperature (TL)
       QCL_F         , &!REAL(KIND=r8)    , INTENT(OUT  ) :: QCL_F(nCols)!       Cloud liquid water content 
       CF_F          , &!REAL(KIND=r8)    , INTENT(OUT  ) :: CF_F(nCols)!       Liquid cloud fraction at processed levels.
       GRID_QC_F     , &!REAL(KIND=r8)    , INTENT(OUT  ) :: GRID_QC_F(nCols)!       Super/subsaturation on processed levels
       BS_F          , &!REAL(KIND=r8)    , INTENT(OUT  ) :: BS_F(nCols)!       Value of bs at processed levels
       INDEX         , &!INTEGER , INTENT(IN) :: INDEX   (nCols,1)! INDEX for  points with non-zero 
       POINTS        , &!INTEGER , INTENT(IN) :: POINTS !       No. of gridpoints with non-zero cloud
       nCols         , &!INTEGER , INTENT(IN) :: nCols !       No. of gridpoints being processed.
       BL_LEVELS     , &!INTEGER , INTENT(IN) :: BL_LEVELS!       No. of boundary layer levels
       K             , &!INTEGER , INTENT(IN) :: K        !        Level no.
       L_eacf        , &!LOGICAL , INTENT(IN) :: L_eacf  !  Use empirically adjusted cloud fraction
       l_mixing_ratio  )!LOGICAL , INTENT(IN) :: L_mixing_ratio   !  Use mixing ratio formulation
    !
    IMPLICIT NONE
    !  Subroutine Arguments:------------------------------------------------
    INTEGER , INTENT(IN) :: nCols !       No. of gridpoints being processed.
    INTEGER , INTENT(IN) :: BL_LEVELS!       No. of boundary layer levels
    INTEGER , INTENT(IN) :: K        !        Level no.
    INTEGER , INTENT(IN) :: POINTS !       No. of gridpoints with non-zero cloud
    INTEGER , INTENT(IN) :: INDEX   (nCols,1)! INDEX for  points with non-zero cloud from lowest model level.
    !
    REAL(KIND=r8)    , INTENT(IN) :: RHCRIT(nCols)! Critical relative humidity.  See the paragraph incorporating
    !       eqs P292.11 to P292.14.
    REAL(KIND=r8)    , INTENT(IN) :: P_F(nCols)!       pressure (Pa).
    REAL(KIND=r8)    , INTENT(IN) :: QSL_F(nCols)!       saturated humidity at temperature TL, and pressure P_F
    REAL(KIND=r8)    , INTENT(IN) :: QN_F(nCols)!       Normalised super/subsaturation ( = QC/BS).
    !
    LOGICAL , INTENT(IN) :: L_eacf  !  Use empirically adjusted cloud fraction
    LOGICAL , INTENT(IN) :: L_mixing_ratio   !  Use mixing ratio formulation
    !
    REAL(KIND=r8)    , INTENT(INOUT) :: Q_F(nCols)!       On input : Vapour + liquid water content (QW) (kg per kg air).
    !       On output: Specific humidity at processed levels
    !                   (kg water per kg air).
    REAL(KIND=r8)    , INTENT(INOUT) :: T_F(nCols)!       On input : Liquid water temperature (TL) (K).
    !       On output: Temperature at processed levels (K).
    !
    REAL(KIND=r8)    , INTENT(OUT  ) :: QCL_F(nCols)!       Cloud liquid water content at processed levels (kg per kg air).
    REAL(KIND=r8)    , INTENT(OUT  ) :: CF_F(nCols)!       Liquid cloud fraction at processed levels.
    REAL(KIND=r8)    , INTENT(OUT  ) :: GRID_QC_F(nCols)!       Super/subsaturation on processed levels. Input initially RMDI.
    REAL(KIND=r8)    , INTENT(OUT  ) :: BS_F(nCols)!       Value of bs at processed levels. Input initialized to RMDI.

    !
    ! Purpose: Calculates liquid cloud water amounts and cloud amounts,
    !          temperature and specific humidity from cloud-conserved and
    !          other model variables. This is done for one model level.
    !
    ! Current Owner of Code: A. C. Bushell
    !
    ! History:
    ! Version   Date     Comment
    !  4.4    14-11-96   Original Code (A. C. Bushell from Wilson/Ballard)
    !  6.2    07-11-05   Include the EACF parametrization for cloud
    !                    fraction.                            D. Wilson
    !  6.4    08-08-06   Include mixing ratio formulation. D. Wilson
    !
    ! Description of Code:
    !   FORTRAN 77  + common extensions also in Fortran90.
    !   This code is written to UMDP3 version 6 programming standards.
    !
    !   Documentation: UMDP No.29
    !
    !  Global Variables:----------------------------------------------------
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
    !
    !
    !  Local parameters and other physical constants------------------------
    REAL(KIND=r8), PARAMETER :: ALPHL=EPSILON*LC/R ! For liquid AlphaL calculation.
    REAL(KIND=r8), PARAMETER :: LCRCP=LC/CP    ! Lat ht of condensation/Cp.
    INTEGER , PARAMETER :: ITS=5    ! Total number of iterations
    INTEGER , PARAMETER :: WTN=0.75_r8  ! Weighting for ALPHAL iteration
    !
    !  Local scalars--------------------------------------------------------
    !
    !  (a) Scalars effectively expanded to workspace by the Cray (using
    !      vector registers).
    REAL(KIND=r8)  :: AL                         ! LOCAL AL (see equation P292.6).
    REAL(KIND=r8)  :: ALPHAL
    ! LOCAL ALPHAL (see equation P292.5).
    REAL(KIND=r8)  :: QN_ADJ
    REAL(KIND=r8)  :: RHCRITX          ! scalar copy of RHCRIT(I,J)
    INTEGER   :: MULTRHC          ! Zero if (nCols*nCols) le 1, else 1
    !
    !  (b) Others.
    INTEGER   I,II,N   ! Loop counters:I,II-horizontal field index.
    !                                       : N - iteration counter.
    !
    !  Local dynamic arrays-------------------------------------------------
    !    7 blocks of REAL(KIND=r8) workspace are required.
    REAL(KIND=r8)  :: P(POINTS)
    !       Pressure  (Pa).
    REAL(KIND=r8)  :: QS(POINTS)
    !       Saturated spec humidity for temp T.
    REAL(KIND=r8)  :: QCN(POINTS)
    !       Cloud water normalised with BS.
    REAL(KIND=r8)  :: T(POINTS)
    !       temperature.
    REAL(KIND=r8)  :: Q(POINTS)
    !       specific humidity.
    REAL(KIND=r8)  :: BS(POINTS)
    !       Sigmas*sqrt(6): sigmas the parametric standard deviation of
    !       local cloud water content fluctuations.
    REAL(KIND=r8)  :: ALPHAL_NM1(POINTS)
    !       ALPHAL at previous iteration.
    !
    !  External subroutine calls: ------------------------------------------
    !      EXTERNAL QSAT_WAT
    !
    !- End of Header
    !
    ! ==Main Block==--------------------------------------------------------
    ! Operate on INDEXed points with non-zero cloud fraction.
    ! ----------------------------------------------------------------------
    IF ( (nCols * nCols)  >   1) THEN
       MULTRHC = 1
    ELSE
       MULTRHC = 0
    END IF
    !
    !        RHCRITX = RHCRIT(1,1)
    ! Points_do1:
    !CDIR NODEP
    DO I=1, POINTS
       II = INDEX(I,1)
       IF ( MULTRHC== 1) THEN
          RHCRITX = RHCRIT(II)
       ELSE
          RHCRITX = RHCRIT(1)
       ENDIF
       P(I)  = P_F(II)
       QCN(I)= QN_F(II)
       ! ----------------------------------------------------------------------
       ! 1. Calculate ALPHAL (eq P292.5) and AL (P292.6).
       !    CAUTION: T_F acts as TL (input value) until update in final section
       !    CAUTION: Q_F acts as QW (input value) until update in final section
       ! ----------------------------------------------------------------------
       !
       ALPHAL = ALPHL * QSL_F(II) / (T_F(II) * T_F(II)) !P292.5
       AL = 1.0_r8 / (1.0_r8 + (LCRCP * ALPHAL))                    ! P292.6
       ALPHAL_NM1(I) = ALPHAL
       !
       ! Rhcrit_if1:
       IF (RHCRITX  <   1.0_r8) THEN
          ! ----------------------------------------------------------------------
          ! 2. Calculate BS (ie. sigma*sqrt(6), where sigma is
          !    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs
          !    P292.15 & 16 if RHcrit < 1.
          ! N.B. QN (input) is initially in QCN
          ! N.B. QN does not depend on AL and so CF and QCN can be calculated
          !      outside the iteration (which is performed in LS_CLD_C).
          !      QN is > -1 for all points processed so CF > 0.
          ! ----------------------------------------------------------------------
          !
          BS(I) = (1.0_r8-RHCRITX) * AL * QSL_F(II)  ! P292.14
          IF (QCN(I)  <=  0.0_r8) THEN
             CF_F(II) = 0.5_r8 * (1.0_r8 + QCN(I)) * (1.0_r8 + QCN(I))
             QCN(I)= (1.0_r8 + QCN(I)) * (1.0_r8 + QCN(I)) * (1.0_r8 + QCN(I)) / 6.0_r8
          ELSEIF (QCN(I)  <   1.0_r8) THEN
             CF_F(II) = 1.0_r8 - 0.5_r8 * (1.0_r8 - QCN(I)) * (1.0_r8 - QCN(I))
             QCN(I)=QCN(I) + (1.0_r8-QCN(I)) * (1.0_r8-QCN(I)) * (1.0_r8-QCN(I))/6.0_r8
          ELSE ! QN  >=  1
             CF_F(II) = 1.0_r8
          END IF ! QCN_if
          !
          ! ----------------------------------------------------------------------
          ! 3.b If necessary, modify cloud fraction using empirically adjusted
          !     cloud fraction parametrization, but keep liquid content the same.
          ! ----------------------------------------------------------------------
          IF (L_eacf) THEN
             ! Adjust QN according to EACF parametrization

             IF (K <= BL_LEVELS) THEN
                QN_ADJ=(QN_F(II)+0.184_r8)/(1.0_r8-0.184_r8)
             ELSE
                QN_ADJ=(QN_F(II)+0.0955_r8)/(1.0_r8-0.0955_r8)
             ENDIF

             !         Calculate cloud fraction using adjusted QN
             IF (QN_ADJ  <=  0.0_r8) THEN
                CF_F(II) = 0.5_r8 * (1.0_r8 + QN_ADJ) * (1.0_r8 + QN_ADJ)
             ELSEIF (QN_ADJ  <   1.0_r8) THEN
                CF_F(II) = 1.0_r8 - 0.5_r8 * (1.0_r8 - QN_ADJ) * (1.0_r8 - QN_ADJ)
             ELSE ! QN_ADJ  >=  1
                CF_F(II) = 1.0_r8
             ENDIF ! QN_ADJ_if

          ENDIF  ! l_eacf
          !
       ELSE ! i.e. if RHcrit = 1
          ! ----------------------------------------------------------------------
          ! 3.a If RHcrit = 1., all points processed have QN > 0 and CF = 1.
          ! ----------------------------------------------------------------------
          BS(I) = AL
          CF_F(II) = 1.0_r8
       END IF ! Rhcrit_if1
       !
       ! ----------------------------------------------------------------------
       ! 3.1 Calculate 1st approx. to qc (store in QCL)
       ! ----------------------------------------------------------------------
       !
       QCL_F(II) = QCN(I) * BS(I)
       !
       ! ----------------------------------------------------------------------
       ! 3.2 Calculate 1st approx. specific humidity (total minus cloud water)
       ! ----------------------------------------------------------------------
       !
       Q(I) = Q_F(II) - QCL_F(II)
       !
       ! ----------------------------------------------------------------------
       ! 3.3 Calculate 1st approx. to temperature, adjusting for latent heating
       ! ----------------------------------------------------------------------
       !
       T(I) = T_F(II) + LCRCP*QCL_F(II)
    END DO ! Points_do1
    !
    ! ----------------------------------------------------------------------
    ! 4. Iteration to find better cloud water values.
    ! ----------------------------------------------------------------------
    ! Its_if:
    IF (ITS  >=  2) THEN
       ! Its_do:
       DO N=2, ITS
          !
          ! DEPENDS ON: qsat_wat_mix
          CALL  qsat_wat_mix ( &
                                !      Output field
               QS                   ,&! REAL(KIND=r8), intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
               T                    ,&! REAL(KIND=r8), intent(in)  :: T(npnts) !  Temperature (K).
               P                    ,&! REAL(KIND=r8), intent(in)  :: P(npnts) !  Pressure (Pa).  
                                !      Array dimensions
               POINTS               ,&! Integer, intent(in) :: npnts  ! Points (=horizontal dimensions) being processed by qSAT scheme.
                                !      logical control
               l_mixing_ratio         )! logical, intent(in)  :: lq_mix
          ! .true. return qsat as a mixing ratio
          ! .false. return qsat as a specific humidity

          !          CALL qsat_wat_mix(QS,T,P,POINTS,l_mixing_ratio)
          ! Points_do2:
          RHCRITX = RHCRIT(1)
          DO I=1, POINTS
             II = INDEX(I,1)
             IF ( MULTRHC== 1) THEN
                RHCRITX = RHCRIT(II)
             ELSE
                RHCRITX = RHCRIT(1)
             ENDIF
             ! T_if:
             IF (T(I)  >   T_F(II)) THEN
                !           NB. T > TL implies cloud fraction > 0.
                ALPHAL = (QS(I) - QSL_F(II)) / (T(I) - T_F(II))
                ALPHAL = WTN * ALPHAL + (1.0_r8 - WTN) * ALPHAL_NM1(I)
                ALPHAL_NM1(I) = ALPHAL
                AL = 1.0_r8 / (1.0_r8 + (LCRCP * ALPHAL))
                ! Rhcrit_if2:
                IF (RHCRITX  <   1.0_r8) THEN
                   BS(I) = (1.0_r8-RHCRITX) * AL * QSL_F(II)
                   !                                                             P292.14
                ELSE
                   BS(I) = AL
                END IF  ! Rhcrit_if2
                !
                ! ----------------------------------------------------------------------
                ! 4.1 Calculate Nth approx. to qc (store in QCL).
                ! ----------------------------------------------------------------------
                !
                QCL_F(II) = QCN(I) * BS(I)
                !
                ! ----------------------------------------------------------------------
                ! 4.2 Calculate Nth approx. spec. humidity (total minus cloud water).
                ! ----------------------------------------------------------------------
                !
                Q(I) = Q_F(II) - QCL_F(II)
                !
                ! ----------------------------------------------------------------------
                ! 4.3 Calculate Nth approx. to temperature, adjusting for latent heating
                ! ----------------------------------------------------------------------
                !
                T(I) = T_F(II) + LCRCP * QCL_F(II)
                !
             END IF ! T_if
          END DO ! Points_do2
       END DO ! Its_do
    END IF ! Its_if
    !
    ! ----------------------------------------------------------------------
    ! 5. Finally scatter back cloud point results to full field arrays.
    !    CAUTION: T_F updated from TL (input) to T (output)
    !    CAUTION: Q_F updated from QW (input) to Q (output)
    ! ----------------------------------------------------------------------
    !
    !DIR$ IVDEP
    ! Points_do3:
    DO I=1,POINTS
       II = INDEX(I,1)
       Q_F(II) = Q(I)
       T_F(II) = T(I)
       GRID_QC_F(II) = BS(I) * QN_F(II)
       BS_F(II) = BS(I)
    END DO ! Points_do3
    !
    RETURN
  END SUBROUTINE LS_CLD_C
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
       &  QmixS            &! REAL(KIND=r8), intent(out)   ::  QmixS(npnts)
                                ! Output Saturation mixing ratio or saturation specific
                                ! humidity at temperature T and pressure P (kg/kg).
                                !      Input fields
       &, T                &! REAL(KIND=r8), intent(in)  :: T(npnts) !  Temperature (K).
       &, P                &! REAL(KIND=r8), intent(in)  :: P(npnts) !  Pressure (Pa).  
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



  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !+ Saturation Specific Humidity Scheme (Qsat): Vapour to Liquid/Ice.
  ! Subroutine Interface:
  SUBROUTINE QSAT (                                                 &
                                !      Output field
       &  QS                                                              &
                                !      Input fields
       &, T, P                                                            &
                                !      Array dimensions
       &, NPNTS                                                           &
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
    !   VALUES IN THE LOOKUP TABLE ARE OVER WATER ABOVE 0 DEG C AND OVER ICE
    !   BELOW THIS TEMPERATURE
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
    !    6.2   03-02-06  Moved to a71_1a. P.Selwood
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
    INTEGER                                                           &
                                !, INTENT(IN)
         &  NPNTS
    !       Points (=horizontal dimensions) being processed by qSAT scheme.
    !
    REAL(KIND=r8)                                                              &
                                !, INTENT(IN)
         &  T(NPNTS)                                                        &
                                !       Temperature (K).
         &, P(NPNTS)
    !       Pressure (Pa).
    !
    ! arguments with intent out. ie: output variables.
    !
    REAL(KIND=r8)                                                              &
                                !, INTENT(OUT)
         &  QS(NPNTS)
    !       SATURATION MIXING RATIO AT TEMPERATURE T AND PRESSURE P (KG/KG)
    !
    !  Local parameters and other physical constants------------------------
    REAL(KIND=r8) ONE_MINUS_EPSILON  ! ONE MINUS THE RATIO OF THE MOLECULAR
    !                               WEIGHTS OF WATER AND DRY AIR
    !
    REAL(KIND=r8) T_LOW        ! LOWEST TEMPERATURE FOR WHICH LOOK-UP TABLE OF
    !                         SATURATION WATER VAPOUR PRESSURE IS VALID (K)
    !
    REAL(KIND=r8) T_HIGH       ! HIGHEST TEMPERATURE FOR WHICH LOOK-UP TABLE OF
    !                         SATURATION WATER VAPOUR PRESSURE IS VALID (K)
    !
    REAL(KIND=r8) DELTA_T      ! TEMPERATURE INCREMENT OF THE LOOK-UP
    ! TABLE OF SATURATION VAPOUR PRESSURES
    !
    INTEGER N         ! SIZE OF LOOK-UP TABLE OF SATURATION
    ! WATER VAPOUR PRESSURES
    !
    PARAMETER ( ONE_MINUS_EPSILON = 1.0_r8 - EPSILON,                    &
         &            T_LOW = 183.15_r8,                                       &
         &            T_HIGH = 338.15_r8,                                      &
         &            DELTA_T = 0.1_r8,                                        &
         &            N = ((T_HIGH - T_LOW + (DELTA_T*0.5_r8))/DELTA_T) + 1.0_r8  &
         &          )    ! gives N=1551
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
    ! ABOVE 0 DEG C VALUES ARE OVER WATER
    !
    ! BELOW 0 DEC C VALUES ARE OVER ICE
    ! ----------------------------------------------------------------------
    ! Note: 0 element is a repeat of 1st element to cater for special case
    !       of low temperatures ( <= T_LOW) for which the array index is
    !       rounded down due to machine precision.
    DATA (ES(IES),IES=    0, 95) / 0.966483E-02_r8,                      &
         &0.966483E-02_r8,0.984279E-02_r8,0.100240E-01_r8,0.102082E-01_r8,0.103957E-01_r8, &
         &0.105865E-01_r8,0.107803E-01_r8,0.109777E-01_r8,0.111784E-01_r8,0.113825E-01_r8, &
         &0.115902E-01_r8,0.118016E-01_r8,0.120164E-01_r8,0.122348E-01_r8,0.124572E-01_r8, &
         &0.126831E-01_r8,0.129132E-01_r8,0.131470E-01_r8,0.133846E-01_r8,0.136264E-01_r8, &
         &0.138724E-01_r8,0.141225E-01_r8,0.143771E-01_r8,0.146356E-01_r8,0.148985E-01_r8, &
         &0.151661E-01_r8,0.154379E-01_r8,0.157145E-01_r8,0.159958E-01_r8,0.162817E-01_r8, &
         &0.165725E-01_r8,0.168680E-01_r8,0.171684E-01_r8,0.174742E-01_r8,0.177847E-01_r8, &
         &0.181008E-01_r8,0.184216E-01_r8,0.187481E-01_r8,0.190801E-01_r8,0.194175E-01_r8, &
         &0.197608E-01_r8,0.201094E-01_r8,0.204637E-01_r8,0.208242E-01_r8,0.211906E-01_r8, &
         &0.215631E-01_r8,0.219416E-01_r8,0.223263E-01_r8,0.227172E-01_r8,0.231146E-01_r8, &
         &0.235188E-01_r8,0.239296E-01_r8,0.243465E-01_r8,0.247708E-01_r8,0.252019E-01_r8, &
         &0.256405E-01_r8,0.260857E-01_r8,0.265385E-01_r8,0.269979E-01_r8,0.274656E-01_r8, &
         &0.279405E-01_r8,0.284232E-01_r8,0.289142E-01_r8,0.294124E-01_r8,0.299192E-01_r8, &
         &0.304341E-01_r8,0.309571E-01_r8,0.314886E-01_r8,0.320285E-01_r8,0.325769E-01_r8, &
         &0.331348E-01_r8,0.337014E-01_r8,0.342771E-01_r8,0.348618E-01_r8,0.354557E-01_r8, &
         &0.360598E-01_r8,0.366727E-01_r8,0.372958E-01_r8,0.379289E-01_r8,0.385717E-01_r8, &
         &0.392248E-01_r8,0.398889E-01_r8,0.405633E-01_r8,0.412474E-01_r8,0.419430E-01_r8, &
         &0.426505E-01_r8,0.433678E-01_r8,0.440974E-01_r8,0.448374E-01_r8,0.455896E-01_r8, &
         &0.463545E-01_r8,0.471303E-01_r8,0.479191E-01_r8,0.487190E-01_r8,0.495322E-01_r8/
    DATA (ES(IES),IES= 96,190) /                                      &
         &0.503591E-01_r8,0.511977E-01_r8,0.520490E-01_r8,0.529145E-01_r8,0.537931E-01_r8, &
         &0.546854E-01_r8,0.555924E-01_r8,0.565119E-01_r8,0.574467E-01_r8,0.583959E-01_r8, &
         &0.593592E-01_r8,0.603387E-01_r8,0.613316E-01_r8,0.623409E-01_r8,0.633655E-01_r8, &
         &0.644053E-01_r8,0.654624E-01_r8,0.665358E-01_r8,0.676233E-01_r8,0.687302E-01_r8, &
         &0.698524E-01_r8,0.709929E-01_r8,0.721490E-01_r8,0.733238E-01_r8,0.745180E-01_r8, &
         &0.757281E-01_r8,0.769578E-01_r8,0.782061E-01_r8,0.794728E-01_r8,0.807583E-01_r8, &
         &0.820647E-01_r8,0.833905E-01_r8,0.847358E-01_r8,0.861028E-01_r8,0.874882E-01_r8, &
         &0.888957E-01_r8,0.903243E-01_r8,0.917736E-01_r8,0.932464E-01_r8,0.947407E-01_r8, &
         &0.962571E-01_r8,0.977955E-01_r8,0.993584E-01_r8,0.100942E+00_r8,0.102551E+00_r8, &
         &0.104186E+00_r8,0.105842E+00_r8,0.107524E+00_r8,0.109231E+00_r8,0.110963E+00_r8, &
         &0.112722E+00_r8,0.114506E+00_r8,0.116317E+00_r8,0.118153E+00_r8,0.120019E+00_r8, &
         &0.121911E+00_r8,0.123831E+00_r8,0.125778E+00_r8,0.127755E+00_r8,0.129761E+00_r8, &
         &0.131796E+00_r8,0.133863E+00_r8,0.135956E+00_r8,0.138082E+00_r8,0.140241E+00_r8, &
         &0.142428E+00_r8,0.144649E+00_r8,0.146902E+00_r8,0.149190E+00_r8,0.151506E+00_r8, &
         &0.153859E+00_r8,0.156245E+00_r8,0.158669E+00_r8,0.161126E+00_r8,0.163618E+00_r8, &
         &0.166145E+00_r8,0.168711E+00_r8,0.171313E+00_r8,0.173951E+00_r8,0.176626E+00_r8, &
         &0.179342E+00_r8,0.182096E+00_r8,0.184893E+00_r8,0.187724E+00_r8,0.190600E+00_r8, &
         &0.193518E+00_r8,0.196473E+00_r8,0.199474E+00_r8,0.202516E+00_r8,0.205604E+00_r8, &
         &0.208730E+00_r8,0.211905E+00_r8,0.215127E+00_r8,0.218389E+00_r8,0.221701E+00_r8/
    DATA (ES(IES),IES=191,285) /                                      &
         &0.225063E+00_r8,0.228466E+00_r8,0.231920E+00_r8,0.235421E+00_r8,0.238976E+00_r8, &
         &0.242580E+00_r8,0.246232E+00_r8,0.249933E+00_r8,0.253691E+00_r8,0.257499E+00_r8, &
         &0.261359E+00_r8,0.265278E+00_r8,0.269249E+00_r8,0.273274E+00_r8,0.277358E+00_r8, &
         &0.281498E+00_r8,0.285694E+00_r8,0.289952E+00_r8,0.294268E+00_r8,0.298641E+00_r8, &
         &0.303078E+00_r8,0.307577E+00_r8,0.312135E+00_r8,0.316753E+00_r8,0.321440E+00_r8, &
         &0.326196E+00_r8,0.331009E+00_r8,0.335893E+00_r8,0.340842E+00_r8,0.345863E+00_r8, &
         &0.350951E+00_r8,0.356106E+00_r8,0.361337E+00_r8,0.366636E+00_r8,0.372006E+00_r8, &
         &0.377447E+00_r8,0.382966E+00_r8,0.388567E+00_r8,0.394233E+00_r8,0.399981E+00_r8, &
         &0.405806E+00_r8,0.411714E+00_r8,0.417699E+00_r8,0.423772E+00_r8,0.429914E+00_r8, &
         &0.436145E+00_r8,0.442468E+00_r8,0.448862E+00_r8,0.455359E+00_r8,0.461930E+00_r8, &
         &0.468596E+00_r8,0.475348E+00_r8,0.482186E+00_r8,0.489124E+00_r8,0.496160E+00_r8, &
         &0.503278E+00_r8,0.510497E+00_r8,0.517808E+00_r8,0.525224E+00_r8,0.532737E+00_r8, &
         &0.540355E+00_r8,0.548059E+00_r8,0.555886E+00_r8,0.563797E+00_r8,0.571825E+00_r8, &
         &0.579952E+00_r8,0.588198E+00_r8,0.596545E+00_r8,0.605000E+00_r8,0.613572E+00_r8, &
         &0.622255E+00_r8,0.631059E+00_r8,0.639962E+00_r8,0.649003E+00_r8,0.658144E+00_r8, &
         &0.667414E+00_r8,0.676815E+00_r8,0.686317E+00_r8,0.695956E+00_r8,0.705728E+00_r8, &
         &0.715622E+00_r8,0.725641E+00_r8,0.735799E+00_r8,0.746082E+00_r8,0.756495E+00_r8, &
         &0.767052E+00_r8,0.777741E+00_r8,0.788576E+00_r8,0.799549E+00_r8,0.810656E+00_r8, &
         &0.821914E+00_r8,0.833314E+00_r8,0.844854E+00_r8,0.856555E+00_r8,0.868415E+00_r8/
    DATA (ES(IES),IES=286,380) /                                      &
         &0.880404E+00_r8,0.892575E+00_r8,0.904877E+00_r8,0.917350E+00_r8,0.929974E+00_r8, &
         &0.942771E+00_r8,0.955724E+00_r8,0.968837E+00_r8,0.982127E+00_r8,0.995600E+00_r8, &
         &0.100921E+01_r8,0.102304E+01_r8,0.103700E+01_r8,0.105116E+01_r8,0.106549E+01_r8, &
         &0.108002E+01_r8,0.109471E+01_r8,0.110962E+01_r8,0.112469E+01_r8,0.113995E+01_r8, &
         &0.115542E+01_r8,0.117107E+01_r8,0.118693E+01_r8,0.120298E+01_r8,0.121923E+01_r8, &
         &0.123569E+01_r8,0.125234E+01_r8,0.126923E+01_r8,0.128631E+01_r8,0.130362E+01_r8, &
         &0.132114E+01_r8,0.133887E+01_r8,0.135683E+01_r8,0.137500E+01_r8,0.139342E+01_r8, &
         &0.141205E+01_r8,0.143091E+01_r8,0.145000E+01_r8,0.146933E+01_r8,0.148892E+01_r8, &
         &0.150874E+01_r8,0.152881E+01_r8,0.154912E+01_r8,0.156970E+01_r8,0.159049E+01_r8, &
         &0.161159E+01_r8,0.163293E+01_r8,0.165452E+01_r8,0.167640E+01_r8,0.169852E+01_r8, &
         &0.172091E+01_r8,0.174359E+01_r8,0.176653E+01_r8,0.178977E+01_r8,0.181332E+01_r8, &
         &0.183709E+01_r8,0.186119E+01_r8,0.188559E+01_r8,0.191028E+01_r8,0.193524E+01_r8, &
         &0.196054E+01_r8,0.198616E+01_r8,0.201208E+01_r8,0.203829E+01_r8,0.206485E+01_r8, &
         &0.209170E+01_r8,0.211885E+01_r8,0.214637E+01_r8,0.217424E+01_r8,0.220242E+01_r8, &
         &0.223092E+01_r8,0.225979E+01_r8,0.228899E+01_r8,0.231855E+01_r8,0.234845E+01_r8, &
         &0.237874E+01_r8,0.240937E+01_r8,0.244040E+01_r8,0.247176E+01_r8,0.250349E+01_r8, &
         &0.253560E+01_r8,0.256814E+01_r8,0.260099E+01_r8,0.263431E+01_r8,0.266800E+01_r8, &
         &0.270207E+01_r8,0.273656E+01_r8,0.277145E+01_r8,0.280671E+01_r8,0.284248E+01_r8, &
         &0.287859E+01_r8,0.291516E+01_r8,0.295219E+01_r8,0.298962E+01_r8,0.302746E+01_r8/
    DATA (ES(IES),IES=381,475) /                                      &
         &0.306579E+01_r8,0.310454E+01_r8,0.314377E+01_r8,0.318351E+01_r8,0.322360E+01_r8, &
         &0.326427E+01_r8,0.330538E+01_r8,0.334694E+01_r8,0.338894E+01_r8,0.343155E+01_r8, &
         &0.347456E+01_r8,0.351809E+01_r8,0.356216E+01_r8,0.360673E+01_r8,0.365184E+01_r8, &
         &0.369744E+01_r8,0.374352E+01_r8,0.379018E+01_r8,0.383743E+01_r8,0.388518E+01_r8, &
         &0.393344E+01_r8,0.398230E+01_r8,0.403177E+01_r8,0.408175E+01_r8,0.413229E+01_r8, &
         &0.418343E+01_r8,0.423514E+01_r8,0.428746E+01_r8,0.434034E+01_r8,0.439389E+01_r8, &
         &0.444808E+01_r8,0.450276E+01_r8,0.455820E+01_r8,0.461423E+01_r8,0.467084E+01_r8, &
         &0.472816E+01_r8,0.478607E+01_r8,0.484468E+01_r8,0.490393E+01_r8,0.496389E+01_r8, &
         &0.502446E+01_r8,0.508580E+01_r8,0.514776E+01_r8,0.521047E+01_r8,0.527385E+01_r8, &
         &0.533798E+01_r8,0.540279E+01_r8,0.546838E+01_r8,0.553466E+01_r8,0.560173E+01_r8, &
         &0.566949E+01_r8,0.573807E+01_r8,0.580750E+01_r8,0.587749E+01_r8,0.594846E+01_r8, &
         &0.602017E+01_r8,0.609260E+01_r8,0.616591E+01_r8,0.623995E+01_r8,0.631490E+01_r8, &
         &0.639061E+01_r8,0.646723E+01_r8,0.654477E+01_r8,0.662293E+01_r8,0.670220E+01_r8, &
         &0.678227E+01_r8,0.686313E+01_r8,0.694495E+01_r8,0.702777E+01_r8,0.711142E+01_r8, &
         &0.719592E+01_r8,0.728140E+01_r8,0.736790E+01_r8,0.745527E+01_r8,0.754352E+01_r8, &
         &0.763298E+01_r8,0.772316E+01_r8,0.781442E+01_r8,0.790676E+01_r8,0.800001E+01_r8, &
         &0.809435E+01_r8,0.818967E+01_r8,0.828606E+01_r8,0.838343E+01_r8,0.848194E+01_r8, &
         &0.858144E+01_r8,0.868207E+01_r8,0.878392E+01_r8,0.888673E+01_r8,0.899060E+01_r8, &
         &0.909567E+01_r8,0.920172E+01_r8,0.930909E+01_r8,0.941765E+01_r8,0.952730E+01_r8/
    DATA (ES(IES),IES=476,570) /                                      &
         &0.963821E+01_r8,0.975022E+01_r8,0.986352E+01_r8,0.997793E+01_r8,0.100937E+02_r8, &
         &0.102105E+02_r8,0.103287E+02_r8,0.104481E+02_r8,0.105688E+02_r8,0.106909E+02_r8, &
         &0.108143E+02_r8,0.109387E+02_r8,0.110647E+02_r8,0.111921E+02_r8,0.113207E+02_r8, &
         &0.114508E+02_r8,0.115821E+02_r8,0.117149E+02_r8,0.118490E+02_r8,0.119847E+02_r8, &
         &0.121216E+02_r8,0.122601E+02_r8,0.124002E+02_r8,0.125416E+02_r8,0.126846E+02_r8, &
         &0.128290E+02_r8,0.129747E+02_r8,0.131224E+02_r8,0.132712E+02_r8,0.134220E+02_r8, &
         &0.135742E+02_r8,0.137278E+02_r8,0.138831E+02_r8,0.140403E+02_r8,0.141989E+02_r8, &
         &0.143589E+02_r8,0.145211E+02_r8,0.146845E+02_r8,0.148501E+02_r8,0.150172E+02_r8, &
         &0.151858E+02_r8,0.153564E+02_r8,0.155288E+02_r8,0.157029E+02_r8,0.158786E+02_r8, &
         &0.160562E+02_r8,0.162358E+02_r8,0.164174E+02_r8,0.166004E+02_r8,0.167858E+02_r8, &
         &0.169728E+02_r8,0.171620E+02_r8,0.173528E+02_r8,0.175455E+02_r8,0.177406E+02_r8, &
         &0.179372E+02_r8,0.181363E+02_r8,0.183372E+02_r8,0.185400E+02_r8,0.187453E+02_r8, &
         &0.189523E+02_r8,0.191613E+02_r8,0.193728E+02_r8,0.195866E+02_r8,0.198024E+02_r8, &
         &0.200200E+02_r8,0.202401E+02_r8,0.204626E+02_r8,0.206871E+02_r8,0.209140E+02_r8, &
         &0.211430E+02_r8,0.213744E+02_r8,0.216085E+02_r8,0.218446E+02_r8,0.220828E+02_r8, &
         &0.223241E+02_r8,0.225671E+02_r8,0.228132E+02_r8,0.230615E+02_r8,0.233120E+02_r8, &
         &0.235651E+02_r8,0.238211E+02_r8,0.240794E+02_r8,0.243404E+02_r8,0.246042E+02_r8, &
         &0.248704E+02_r8,0.251390E+02_r8,0.254109E+02_r8,0.256847E+02_r8,0.259620E+02_r8, &
         &0.262418E+02_r8,0.265240E+02_r8,0.268092E+02_r8,0.270975E+02_r8,0.273883E+02_r8/
    DATA (ES(IES),IES=571,665) /                                      &
         &0.276822E+02_r8,0.279792E+02_r8,0.282789E+02_r8,0.285812E+02_r8,0.288867E+02_r8, &
         &0.291954E+02_r8,0.295075E+02_r8,0.298222E+02_r8,0.301398E+02_r8,0.304606E+02_r8, &
         &0.307848E+02_r8,0.311119E+02_r8,0.314424E+02_r8,0.317763E+02_r8,0.321133E+02_r8, &
         &0.324536E+02_r8,0.327971E+02_r8,0.331440E+02_r8,0.334940E+02_r8,0.338475E+02_r8, &
         &0.342050E+02_r8,0.345654E+02_r8,0.349295E+02_r8,0.352975E+02_r8,0.356687E+02_r8, &
         &0.360430E+02_r8,0.364221E+02_r8,0.368042E+02_r8,0.371896E+02_r8,0.375790E+02_r8, &
         &0.379725E+02_r8,0.383692E+02_r8,0.387702E+02_r8,0.391744E+02_r8,0.395839E+02_r8, &
         &0.399958E+02_r8,0.404118E+02_r8,0.408325E+02_r8,0.412574E+02_r8,0.416858E+02_r8, &
         &0.421188E+02_r8,0.425551E+02_r8,0.429962E+02_r8,0.434407E+02_r8,0.438910E+02_r8, &
         &0.443439E+02_r8,0.448024E+02_r8,0.452648E+02_r8,0.457308E+02_r8,0.462018E+02_r8, &
         &0.466775E+02_r8,0.471582E+02_r8,0.476428E+02_r8,0.481313E+02_r8,0.486249E+02_r8, &
         &0.491235E+02_r8,0.496272E+02_r8,0.501349E+02_r8,0.506479E+02_r8,0.511652E+02_r8, &
         &0.516876E+02_r8,0.522142E+02_r8,0.527474E+02_r8,0.532836E+02_r8,0.538266E+02_r8, &
         &0.543737E+02_r8,0.549254E+02_r8,0.554839E+02_r8,0.560456E+02_r8,0.566142E+02_r8, &
         &0.571872E+02_r8,0.577662E+02_r8,0.583498E+02_r8,0.589392E+02_r8,0.595347E+02_r8, &
         &0.601346E+02_r8,0.607410E+02_r8,0.613519E+02_r8,0.619689E+02_r8,0.625922E+02_r8, &
         &0.632204E+02_r8,0.638550E+02_r8,0.644959E+02_r8,0.651418E+02_r8,0.657942E+02_r8, &
         &0.664516E+02_r8,0.671158E+02_r8,0.677864E+02_r8,0.684624E+02_r8,0.691451E+02_r8, &
         &0.698345E+02_r8,0.705293E+02_r8,0.712312E+02_r8,0.719398E+02_r8,0.726542E+02_r8/
    DATA (ES(IES),IES=666,760) /                                      &
         &0.733754E+02_r8,0.741022E+02_r8,0.748363E+02_r8,0.755777E+02_r8,0.763247E+02_r8, &
         &0.770791E+02_r8,0.778394E+02_r8,0.786088E+02_r8,0.793824E+02_r8,0.801653E+02_r8, &
         &0.809542E+02_r8,0.817509E+02_r8,0.825536E+02_r8,0.833643E+02_r8,0.841828E+02_r8, &
         &0.850076E+02_r8,0.858405E+02_r8,0.866797E+02_r8,0.875289E+02_r8,0.883827E+02_r8, &
         &0.892467E+02_r8,0.901172E+02_r8,0.909962E+02_r8,0.918818E+02_r8,0.927760E+02_r8, &
         &0.936790E+02_r8,0.945887E+02_r8,0.955071E+02_r8,0.964346E+02_r8,0.973689E+02_r8, &
         &0.983123E+02_r8,0.992648E+02_r8,0.100224E+03_r8,0.101193E+03_r8,0.102169E+03_r8, &
         &0.103155E+03_r8,0.104150E+03_r8,0.105152E+03_r8,0.106164E+03_r8,0.107186E+03_r8, &
         &0.108217E+03_r8,0.109256E+03_r8,0.110303E+03_r8,0.111362E+03_r8,0.112429E+03_r8, &
         &0.113503E+03_r8,0.114588E+03_r8,0.115684E+03_r8,0.116789E+03_r8,0.117903E+03_r8, &
         &0.119028E+03_r8,0.120160E+03_r8,0.121306E+03_r8,0.122460E+03_r8,0.123623E+03_r8, &
         &0.124796E+03_r8,0.125981E+03_r8,0.127174E+03_r8,0.128381E+03_r8,0.129594E+03_r8, &
         &0.130822E+03_r8,0.132058E+03_r8,0.133306E+03_r8,0.134563E+03_r8,0.135828E+03_r8, &
         &0.137109E+03_r8,0.138402E+03_r8,0.139700E+03_r8,0.141017E+03_r8,0.142338E+03_r8, &
         &0.143676E+03_r8,0.145025E+03_r8,0.146382E+03_r8,0.147753E+03_r8,0.149133E+03_r8, &
         &0.150529E+03_r8,0.151935E+03_r8,0.153351E+03_r8,0.154783E+03_r8,0.156222E+03_r8, &
         &0.157678E+03_r8,0.159148E+03_r8,0.160624E+03_r8,0.162117E+03_r8,0.163621E+03_r8, &
         &0.165142E+03_r8,0.166674E+03_r8,0.168212E+03_r8,0.169772E+03_r8,0.171340E+03_r8, &
         &0.172921E+03_r8,0.174522E+03_r8,0.176129E+03_r8,0.177755E+03_r8,0.179388E+03_r8/
    DATA (ES(IES),IES=761,855) /                                      &
         &0.181040E+03_r8,0.182707E+03_r8,0.184382E+03_r8,0.186076E+03_r8,0.187782E+03_r8, &
         &0.189503E+03_r8,0.191240E+03_r8,0.192989E+03_r8,0.194758E+03_r8,0.196535E+03_r8, &
         &0.198332E+03_r8,0.200141E+03_r8,0.201963E+03_r8,0.203805E+03_r8,0.205656E+03_r8, &
         &0.207532E+03_r8,0.209416E+03_r8,0.211317E+03_r8,0.213236E+03_r8,0.215167E+03_r8, &
         &0.217121E+03_r8,0.219087E+03_r8,0.221067E+03_r8,0.223064E+03_r8,0.225080E+03_r8, &
         &0.227113E+03_r8,0.229160E+03_r8,0.231221E+03_r8,0.233305E+03_r8,0.235403E+03_r8, &
         &0.237520E+03_r8,0.239655E+03_r8,0.241805E+03_r8,0.243979E+03_r8,0.246163E+03_r8, &
         &0.248365E+03_r8,0.250593E+03_r8,0.252830E+03_r8,0.255093E+03_r8,0.257364E+03_r8, &
         &0.259667E+03_r8,0.261979E+03_r8,0.264312E+03_r8,0.266666E+03_r8,0.269034E+03_r8, &
         &0.271430E+03_r8,0.273841E+03_r8,0.276268E+03_r8,0.278722E+03_r8,0.281185E+03_r8, &
         &0.283677E+03_r8,0.286190E+03_r8,0.288714E+03_r8,0.291266E+03_r8,0.293834E+03_r8, &
         &0.296431E+03_r8,0.299045E+03_r8,0.301676E+03_r8,0.304329E+03_r8,0.307006E+03_r8, &
         &0.309706E+03_r8,0.312423E+03_r8,0.315165E+03_r8,0.317930E+03_r8,0.320705E+03_r8, &
         &0.323519E+03_r8,0.326350E+03_r8,0.329199E+03_r8,0.332073E+03_r8,0.334973E+03_r8, &
         &0.337897E+03_r8,0.340839E+03_r8,0.343800E+03_r8,0.346794E+03_r8,0.349806E+03_r8, &
         &0.352845E+03_r8,0.355918E+03_r8,0.358994E+03_r8,0.362112E+03_r8,0.365242E+03_r8, &
         &0.368407E+03_r8,0.371599E+03_r8,0.374802E+03_r8,0.378042E+03_r8,0.381293E+03_r8, &
         &0.384588E+03_r8,0.387904E+03_r8,0.391239E+03_r8,0.394604E+03_r8,0.397988E+03_r8, &
         &0.401411E+03_r8,0.404862E+03_r8,0.408326E+03_r8,0.411829E+03_r8,0.415352E+03_r8/
    DATA (ES(IES),IES=856,950) /                                      &
         &0.418906E+03_r8,0.422490E+03_r8,0.426095E+03_r8,0.429740E+03_r8,0.433398E+03_r8, &
         &0.437097E+03_r8,0.440827E+03_r8,0.444570E+03_r8,0.448354E+03_r8,0.452160E+03_r8, &
         &0.455999E+03_r8,0.459870E+03_r8,0.463765E+03_r8,0.467702E+03_r8,0.471652E+03_r8, &
         &0.475646E+03_r8,0.479674E+03_r8,0.483715E+03_r8,0.487811E+03_r8,0.491911E+03_r8, &
         &0.496065E+03_r8,0.500244E+03_r8,0.504448E+03_r8,0.508698E+03_r8,0.512961E+03_r8, &
         &0.517282E+03_r8,0.521617E+03_r8,0.525989E+03_r8,0.530397E+03_r8,0.534831E+03_r8, &
         &0.539313E+03_r8,0.543821E+03_r8,0.548355E+03_r8,0.552938E+03_r8,0.557549E+03_r8, &
         &0.562197E+03_r8,0.566884E+03_r8,0.571598E+03_r8,0.576351E+03_r8,0.581131E+03_r8, &
         &0.585963E+03_r8,0.590835E+03_r8,0.595722E+03_r8,0.600663E+03_r8,0.605631E+03_r8, &
         &0.610641E+03_r8,0.615151E+03_r8,0.619625E+03_r8,0.624140E+03_r8,0.628671E+03_r8, &
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
       ITABLE = ATABLE
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
            ( MAX(P(I),QS(I)) - ONE_MINUS_EPSILON*QS(I) )
       !
    END DO ! Npnts_do_1
    !
    RETURN
  END SUBROUTINE QSAT
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
    INTEGER , INTENT(in) :: ncol                 ! Number of longitudes
    INTEGER , INTENT(in) :: nCols
    INTEGER , INTENT(in) :: kMax
    REAL(r8), INTENT(in) :: piln (nCols,kMax+1)  ! Log interface pressures
    REAL(r8), INTENT(in) :: pint (nCols,kMax+1)  ! Interface pressures
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

END MODULE CloudFraction_UKME

