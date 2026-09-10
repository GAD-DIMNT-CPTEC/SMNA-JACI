!
!-----------------------------------------------------------------------
!
MODULE Pbl_MellorYamada1

  USE Constants, ONLY :     &
       r8

    IMPLICIT NONE
  SAVE

  PRIVATE
  !
  !-----------------------------------------------------------------------
  !

  !  2. Following are constants for use in defining REAL(KIND=r8) number bounds.

  !  A really small number.

  REAL(KIND=r8)    , PARAMETER :: epsilon  = 1.E-15_r8
  !
  !  4. Following is information related to the physical constants.
  !
  !  These are the physical constants used within the model.
  !
  ! JM NOTE -- can we name this grav instead?
  !
  REAL(KIND=r8)    , PARAMETER :: g = 9.81_r8  ! acceleration due to gravity (m {s}^-2)
  !
  REAL(KIND=r8)    , PARAMETER :: r_d          = 287.0_r8
  REAL(KIND=r8)    , PARAMETER :: cp           = 7.0_r8*r_d/2.0_r8

  REAL(KIND=r8)    , PARAMETER :: r_v          = 461.6_r8
  REAL(KIND=r8)    , PARAMETER :: cv           = cp-r_d
  REAL(KIND=r8)    , PARAMETER :: cpv          = 4.0_r8*r_v
  REAL(KIND=r8)    , PARAMETER :: cvv          = cpv-r_v
  REAL(KIND=r8)    , PARAMETER :: cvpm         = -cv/cp
  REAL(KIND=r8)    , PARAMETER :: cliq         = 4190.0_r8
  REAL(KIND=r8)    , PARAMETER :: cice         = 2106.0_r8
  REAL(KIND=r8)    , PARAMETER :: psat         = 610.78_r8
  REAL(KIND=r8)    , PARAMETER :: rcv          = r_d/cv
  REAL(KIND=r8)    , PARAMETER :: rcp          = r_d/cp
  REAL(KIND=r8)    , PARAMETER :: rovg         = r_d/g
  REAL(KIND=r8)    , PARAMETER :: c2           = cp * rcv

  REAL(KIND=r8)    , PARAMETER :: p1000mb      = 100000.0_r8
  REAL(KIND=r8)    , PARAMETER :: t0           = 300.0_r8
  REAL(KIND=r8)    , PARAMETER :: p0           = p1000mb
  REAL(KIND=r8)    , PARAMETER :: cpovcv       = cp/(cp-r_d)
  REAL(KIND=r8)    , PARAMETER :: cvovcp       = 1.0_r8/cpovcv
  REAL(KIND=r8)    , PARAMETER :: rvovrd       = r_v/r_d

  REAL(KIND=r8)    , PARAMETER :: reradius     = 1.0_r8/6370.e03_r8 

  REAL(KIND=r8)    , PARAMETER :: asselin      = 0.025_r8
  !   REAL(KIND=r8)    , PARAMETER :: asselin      = 0.0_r8
  REAL(KIND=r8)    , PARAMETER :: cb           = 25.0_r8

  REAL(KIND=r8)    , PARAMETER :: XLV0         = 3.15E6_r8
  REAL(KIND=r8)    , PARAMETER :: XLV1         = 2370.0_r8
  REAL(KIND=r8)    , PARAMETER :: XLS0         = 2.905E6_r8
  REAL(KIND=r8)    , PARAMETER :: XLS1         = 259.532_r8

  REAL(KIND=r8)    , PARAMETER :: XLS          = 2.85E6_r8
  REAL(KIND=r8)    , PARAMETER :: XLV          = 2.5E6_r8
  REAL(KIND=r8)    , PARAMETER :: XLF          = 3.50E5_r8

  REAL(KIND=r8)    , PARAMETER :: rhowater     = 1000.0_r8
  REAL(KIND=r8)    , PARAMETER :: rhosnow      = 100.0_r8
  REAL(KIND=r8)    , PARAMETER :: rhoair0      = 1.28_r8

  REAL(KIND=r8)    , PARAMETER :: DEGRAD       = 3.1415926_r8/180.0_r8
  REAL(KIND=r8)    , PARAMETER :: DPD          = 360.0_r8/365.0_r8

  REAL(KIND=r8)    , PARAMETER ::  SVP1=0.6112_r8
  REAL(KIND=r8)    , PARAMETER ::  SVP2=17.67_r8
  REAL(KIND=r8)    , PARAMETER ::  SVP3=29.65_r8
  REAL(KIND=r8)    , PARAMETER ::  SVPT0=273.15_r8
  REAL(KIND=r8)    , PARAMETER ::  EP_1=R_v/R_d-1.0_r8
  REAL(KIND=r8)    , PARAMETER ::  EP_2=R_d/R_v
  REAL(KIND=r8)    , PARAMETER ::  KARMAN=0.4_r8
  REAL(KIND=r8)    , PARAMETER ::  EOMEG=7.2921E-5_r8
  REAL(KIND=r8)    , PARAMETER ::  STBOLT=5.67051E-8_r8

  ! proportionality constants for eddy viscosity coefficient calc
  REAL(KIND=r8)    , PARAMETER ::  c_s = 0.25_r8  ! turbulence parameterization constant, for smagorinsky
  REAL(KIND=r8)    , PARAMETER ::  c_k = 0.15_r8  ! turbulence parameterization constant, for TKE
  REAL(KIND=r8)    , PARAMETER ::  prandtl = 1.0_r8/3.0_r8
  ! constants for w-damping option
  REAL(KIND=r8)    , PARAMETER ::  w_alpha = 0.3_r8 ! strength m/s/s
  REAL(KIND=r8)    , PARAMETER ::  w_beta  = 1.0_r8 ! activation cfl number

  REAL(KIND=r8) , PARAMETER ::  pq0=379.90516_r8
  !      REAL(KIND=r8) , PARAMETER ::  epsq2=0.2_r8
  REAL(KIND=r8) , PARAMETER ::  epsq2=0.02_r8
  REAL(KIND=r8) , PARAMETER ::  a2=17.2693882_r8
  REAL(KIND=r8) , PARAMETER ::  a3=273.16_r8
  REAL(KIND=r8) , PARAMETER ::  a4=35.86_r8
  REAL(KIND=r8) , PARAMETER ::  epsq=1.e-12_r8
  REAL(KIND=r8) , PARAMETER ::  p608=rvovrd-1.0_r8
  REAL(KIND=r8) , PARAMETER ::  climit=1.e-20_r8
  REAL(KIND=r8) , PARAMETER ::  cm1=2937.4_r8
  REAL(KIND=r8) , PARAMETER ::  cm2=4.9283_r8
  REAL(KIND=r8) , PARAMETER ::  cm3=23.5518_r8
  !       REAL(KIND=r8) , PARAMETER ::  defc=8.0_r8
  !       REAL(KIND=r8) , PARAMETER ::  defm=32.0_r8
  REAL(KIND=r8) , PARAMETER ::  defc=0.0_r8
  REAL(KIND=r8) , PARAMETER ::  defm=99999.0_r8
  REAL(KIND=r8) , PARAMETER ::  epsfc=1.0_r8/1.05_r8
  REAL(KIND=r8) , PARAMETER ::  epswet=0.0_r8
  REAL(KIND=r8) , PARAMETER ::  fcdif=1.0_r8/3.0_r8
  !       REAL(KIND=r8) , PARAMETER ::  fcm=0.003_r8
  REAL(KIND=r8) , PARAMETER ::  fcm=0.0_r8
  REAL(KIND=r8) , PARAMETER ::  gma=-r_d*(1.0_r8-rcp)*0.5_r8
  REAL(KIND=r8) , PARAMETER ::  p400=40000.0_r8
  REAL(KIND=r8) , PARAMETER ::  phitp=15000.0_r8
  REAL(KIND=r8) , PARAMETER ::  pi2=2.0_r8*3.1415926_r8
  REAL(KIND=r8) , PARAMETER ::  plbtm=105000.0_r8
  REAL(KIND=r8) , PARAMETER ::  plomd=64200.0_r8
  REAL(KIND=r8) , PARAMETER ::  pmdhi=35000.0_r8
  REAL(KIND=r8) , PARAMETER ::  q2ini=0.50_r8
  REAL(KIND=r8) , PARAMETER ::  rfcp=0.25_r8/cp
  REAL(KIND=r8) , PARAMETER ::  rhcrit_land=0.75_r8
  REAL(KIND=r8) , PARAMETER ::  rhcrit_sea=0.80_r8
  REAL(KIND=r8) , PARAMETER ::  rlag=14.8125_r8
  REAL(KIND=r8) , PARAMETER ::  rlx=0.90_r8
  REAL(KIND=r8) , PARAMETER ::  scq2=50.0_r8
  REAL(KIND=r8) , PARAMETER ::  slopht=0.001_r8
  REAL(KIND=r8) , PARAMETER ::  tlc=2.0_r8*0.703972477_r8
  REAL(KIND=r8) , PARAMETER ::  wa=0.15_r8
  REAL(KIND=r8) , PARAMETER ::  wght=0.35_r8
  REAL(KIND=r8) , PARAMETER ::  wpc=0.075_r8
  REAL(KIND=r8) , PARAMETER ::  z0land=0.10_r8
  REAL(KIND=r8) , PARAMETER ::  z0max=0.01_r8
  REAL(KIND=r8) , PARAMETER ::  z0sea=0.001_r8
  !#endif
  !
  !-----------------------------------------------------------------------
  !
  ! REFERENCES:  Janjic (2001), NCEP Office Note 437
  !              Mellor and Yamada (1982), Rev. Geophys. Space Phys.
  !
  ! ABSTRACT:
  !     MYJ UPDATES THE TURBULENT KINETIC ENERGY WITH THE PRODUCTION/
  !     DISSIPATION TERM AND THE VERTICAL DIFFUSION TERM
  !     (USING AN IMPLICIT FORMULATION) FROM MELLOR-YAMADA
  !     LEVEL 2.5 AS EXTENDED BY JANJIC.  EXCHANGE COEFFICIENTS FOR
  !     THE SURFACE AND FOR ALL LAYER INTERFACES ARE COMPUTED FROM
  !     MONIN-OBUKHOV THEORY.
  !     THE TURBULENT VERTICAL EXCHANGE IS THEN EXECUTED.
  !
  !-----------------------------------------------------------------------
  !
  INTEGER :: ITRMX=10 ! Iteration count for mixing length computation
  !
  !     REAL(KIND=r8),PARAMETER :: G=9.81_r8,PI=3.1415926_r8,R_D=287.04_r8,R_V=461.6_r8        &
  !    &                 ,VKARMAN=0.4_r8
  REAL(KIND=r8),PARAMETER :: PI=3.1415926_r8,VKARMAN=0.4_r8
  !     REAL(KIND=r8),PARAMETER :: CP=7.0_r8*R_D/2.0_r8
  REAL(KIND=r8),PARAMETER :: CAPA=R_D/CP
  REAL(KIND=r8),PARAMETER :: RLIVWV=XLS/XLV,ELOCP=2.72E6_r8/CP
  REAL(KIND=r8),PARAMETER :: EPS1=1.E-12_r8,EPS2=0.0_r8
  REAL(KIND=r8),PARAMETER :: EPSL=0.10_r8,EPSRU=1.E-7_r8,EPSRS=1.E-7_r8               &
       &                 ,EPSTRB=1.E-24_r8
  REAL(KIND=r8),PARAMETER :: EPSA=1.E-8_r8,EPSIT=1.E-4_r8,EPSU2=1.E-4_r8,EPSUST=0.07_r8  &
       &                 ,FH=1.01_r8
  REAL(KIND=r8),PARAMETER :: ALPH=0.30_r8,BETA=1.0_r8/273.0_r8,EL0MAX=1000.0_r8,EL0MIN=1.0_r8   &
       &                 ,ELFC=0.23_r8*0.5_r8,GAM1=0.2222222222222222222_r8        &
       &                 ,PRT=1.0_r8
  REAL(KIND=r8),PARAMETER :: A1=0.659888514560862645_r8                         &
       &                 ,A2x=0.6574209922667784586_r8                       &
       &                 ,B1=11.87799326209552761_r8                         &
       &                 ,B2=7.226971804046074028_r8                         &
       &                 ,C1=0.000830955950095854396_r8
  REAL(KIND=r8),PARAMETER :: A2S=17.2693882_r8,A3S=273.16_r8,A4S=35.86_r8
  REAL(KIND=r8),PARAMETER :: ELZ0=0.0_r8,ESQ=5.0_r8,EXCM=0.001_r8                      &
       &                 ,FHNEU=0.8_r8,GLKBR=10.0_r8,GLKBS=30.0_r8                   &
       &                 ,QVISC=2.1E-5_r8,RFC=0.191_r8,RIC=0.505_r8,SMALL=0.35_r8     &
       &                 ,SQPR=0.84_r8,SQSC=0.84_r8,SQVISC=258.2_r8,TVISC=2.1E-5_r8   &
       &                 ,USTC=0.7_r8,USTR=0.225_r8,VISC=1.5E-5_r8                 &
       &                 ,WOLD=0.15_r8,WWST=1.2_r8,ZTMAX=1.0_r8,ZTFC=1.0_r8,ZTMIN=-5.0_r8
  !
  REAL(KIND=r8),PARAMETER :: SEAFC=0.98_r8,PQ0SEA=PQ0*SEAFC
  !
  REAL(KIND=r8),PARAMETER :: BTG=BETA*G,CZIV=SMALL*GLKBS                     &
       !    &                 ,EP_1=R_V/R_D-1.0_r8,ESQHF=0.5_r8*5.0_r8,GRRS=GLKBR/GLKBS  &
  &                 ,ESQHF=0.5_r8*5.0_r8,GRRS=GLKBR/GLKBS                  &
       &                 ,RB1=1.0_r8/B1,RTVISC=1.0_r8/TVISC,RVISC=1.0_r8/VISC         &
       &                 ,ZQRZT=SQSC/SQPR
  !
  REAL(KIND=r8),PARAMETER :: ADNH= 9.0_r8*A1*A2x*A2x*(12.0_r8*A1+3.0_r8*B2)*BTG*BTG      &                  
       &                 ,ADNM=18.0_r8*A1*A1*A2x*(B2-3.0_r8*A2x)*BTG              & 
       &                 ,ANMH=-9.0_r8*A1*A2x*A2x*BTG*BTG                     &
       &                 ,ANMM=-3.0_r8*A1*A2x*(3.0_r8*A2x+3.0_r8*B2*C1+18.0_r8*A1*C1-B2)  &
       &                                *BTG                              &   
       &                 ,BDNH= 3.0_r8*A2x*(7.0_r8*A1+B2)*BTG                     &
       &                 ,BDNM= 6.0_r8*A1*A1                                  &
       &                 ,BEQH= A2x*B1*BTG+3.0_r8*A2x*(7.0_r8*A1+B2)*BTG          &
       &                 ,BEQM=-A1*B1*(1.0_r8-3.0_r8*C1)+6.0_r8*A1*A1                 &
       &                 ,BNMH=-A2x*BTG                                   &     
       &                 ,BNMM=A1*(1.0_r8-3.0_r8*C1)                              &
       &                 ,BSHH=9.0_r8*A1*A2x*A2x*BTG                          &
       &                 ,BSHM=18.0_r8*A1*A1*A2x*C1                           &
       &                 ,BSMH=-3.0_r8*A1*A2x*(3.0_r8*A2x+3.0_r8*B2*C1+12.0_r8*A1*C1-B2)  &
       &                                *BTG                              &
       &                 ,CESH=A2x                                        &
       &                 ,CESM=A1*(1.0_r8-3.0_r8*C1)                              &
       &                 ,CNV=EP_1*G/BTG                                  &
       &                 ,ELFCS=VKARMAN*BTG                               &
       &                 ,FZQ1=RTVISC*QVISC*ZQRZT                         &
       &                 ,FZQ2=RTVISC*QVISC*ZQRZT                         &
       &                 ,FZT1=RVISC *TVISC*SQPR                          &
       &                 ,FZT2=CZIV*GRRS*TVISC*SQPR                       &
       &                 ,FZU1=CZIV*VISC                                  &
       &                 ,PIHF=0.5_r8*PI                                     &
       &                 ,RFAC=RIC/(FHNEU*RFC*RFC)                        &
       &                 ,RQVISC=1.0_r8/QVISC                                 &
       &                 ,RRIC=1.0_r8/RIC                                     &
       &                 ,USTFC=0.018_r8/G                                   &
       &                 ,WNEW=1.0_r8-WOLD                                    &
       &                 ,WWST2=WWST*WWST
  !
  !-----------------------------------------------------------------------
  !***  FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2
  !-----------------------------------------------------------------------
  !
  REAL(KIND=r8),PARAMETER :: AEQH=9.0_r8*A1*A2x*A2x*B1*BTG*BTG                   &
       &                      +9.0_r8*A1*A2x*A2x*(12.0_r8*A1+3.0_r8*B2)*BTG*BTG       &
       &                 ,AEQM=3.0_r8*A1*A2x*B1*(3.0_r8*A2x+3.0_r8*B2*C1+18.0_r8*A1*C1-B2)&
       &                      *BTG+18.0_r8*A1*A1*A2x*(B2-3.0_r8*A2x)*BTG
  !
  !-----------------------------------------------------------------------
  !***  FORBIDDEN TURBULENCE AREA
  !-----------------------------------------------------------------------
  !
  REAL(KIND=r8),PARAMETER :: REQU=-AEQH/AEQM                                 &
       &                 ,EPSGH=1.E-9_r8,EPSGM=REQU*EPSGH
  !
  !-----------------------------------------------------------------------
  !***  NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT
  !-----------------------------------------------------------------------
  ! 
  REAL(KIND=r8),PARAMETER :: UBRYL=(18.0_r8*REQU*A1*A1*A2x*B2*C1*BTG             &
       &                         +9.0_r8*A1*A2x*A2x*B2*BTG*BTG)               &
       &                        /(REQU*ADNM+ADNH)                         &
       &                 ,UBRY=(1.0_r8+EPSRS)*UBRYL,UBRY3=3.0_r8*UBRY
  !
  REAL(KIND=r8),PARAMETER :: AUBH=27.0_r8*A1*A2x*A2x*B2*BTG*BTG-ADNH*UBRY3       &
       &                 ,AUBM=54.0_r8*A1*A1*A2x*B2*C1*BTG -ADNM*UBRY3        &
       &                 ,BUBH=(9.0_r8*A1*A2x+3.0_r8*A2x*B2)*BTG-BDNH*UBRY3       &
       &                 ,BUBM=18.0_r8*A1*A1*C1           -BDNM*UBRY3         &
       &                 ,CUBR=1.0_r8                     -     UBRY3         &
       &                 ,RCUBR=1.0_r8/CUBR
  !
  !-----------------------------------------------------------------------
  !



  PUBLIC :: InitPbl_MellorYamada1
  PUBLIC :: MellorYamada1
CONTAINS
  SUBROUTINE MellorYamada1(iMax    , & !INTENT(IN   ) 
    kMax    , & !INTENT(IN   ) 
                    DT      , & !INTENT(IN   ) :: DT! time step (second)
                    STEPBL  , & !INTENT(IN   ) :: STEPBL! bldt (max_dom)= 0,; minutes between boundary-layer physics calls
    HT      , & !INTENT(IN   ) :: HT   (1:iMax)!"HGT" "Terrain Height"   "m"
    DZ      , & !INTENT(IN   ) :: DZ (1:iMax,1:kMax)!-- dz8w dz between full levels (m)
    PMID    , & !INTENT(IN   ) :: PMID(1:iMax,1:kMax)!-- p_phy pressure (Pa)
    PINT    , & !INTENT(IN   ) :: PINT  (1:iMax,1:kMax)!-- p8w pressure at full levels (Pa)
    TH      , & !INTENT(IN   ) :: TH(1:iMax,1:kMax)! potential temperature (K)
    T       , & !INTENT(IN   ) :: T(1:iMax,1:kMax) !t_phy         temperature (K)
    QV      , & !INTENT(IN   ) :: QV(1:iMax,1:kMax)! Qv         water vapor mixing ratio (kg/kg)
    CWM     , & !INTENT(IN   ) :: CWM(1:iMax,1:kMax)! cloud water mixing ratio (kg/kg)
    U       , & !INTENT(IN   ) :: U(1:iMax,1:kMax)! u-velocity interpolated to theta points (m/s)
    V       , & !INTENT(IN   ) :: V(1:iMax,1:kMax)! v-velocity interpolated to theta points (m/s)
    TSK     , & !INTENT(IN   ) :: TSK    (1:iMax)!-- TSK surface temperature (K)
    CHKLOWQ , & !INTENT(IN   ) :: CHKLOWQ(1:iMax)!-- CHKLOWQ - is either 0 or 1 (so far set equal to 1).
                    LOWLYR  , & !INTENT(IN   ) :: LOWLYR (1:iMax)!-- lowlyr index of lowest model layer above ground
    XLAND   , & !INTENT(IN   ) :: XLAND  (1:iMax)!-- XLAND         land mask (1 for land, 2 for water)
    SICE    , & !INTENT(IN   ) :: SICE (1:iMax) !-- SICE liquid water-equivalent ice  depth (m)
    SNOW    , & !INTENT(IN   ) :: SNOW (1:iMax) !-- SNOW liquid water-equivalent snow depth (m)
    ELFLX   , & !INTENT(IN   ) :: ELFLX(1:iMax)!-- ELFLX--LH net upward latent heat flux at surface (W/m^2)
    bps     , & 
    colrad  , &     
    QSFC    , & !INTENT(INOUT) :: QSFC   (1:iMax)!-- qsfc specific humidity at lower boundary (kg/kg)
    THZ0    , & !INTENT(INOUT) :: THZ0   (1:iMax)!-- thz0 potential temperature at roughness length (K)
    QZ0     , & !INTENT(INOUT) :: QZ0    (1:iMax)!-- QZ0 specific humidity at roughness length (kg/kg)
    UZ0     , & !INTENT(INOUT) :: UZ0    (1:iMax)!-- uz0 u wind component at roughness length (m/s)
    VZ0     , & !INTENT(INOUT) :: VZ0    (1:iMax)!-- vz0 v wind component at roughness length (m/s)
    TKE_MYJ , & !INTENT(INOUT) :: TKE_MYJ(1:iMax,1:kMax)!-- tke_myjturbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
    EXCH_H  , & !INTENT(INOUT) :: EXCH_H (1:iMax,1:kMax)
    USTAR   , & !INTENT(INOUT) :: USTAR  (1:iMax)!-- UST           u* in similarity theory (m/s)
    CT      , & !INTENT(INOUT) :: CT     (1:iMax)
                    AKHS    , & !INTENT(INOUT) :: AKHS(1:iMax)!-- akhssfc exchange coefficient of heat/moisture from MYJ
    AKMS    , & !INTENT(INOUT) :: AKMS(1:iMax)!-- akmssfc exchange coefficient of momentum from MYJ
    EL_MYJ  , & !INTENT(OUT  ) :: EL_MYJ   (1:iMax,1:kMax)! mixing length from Mellor-Yamada-Janjic (MYJ) (m)
    PBLH    , & !INTENT(OUT  ) :: PBLH(1:iMax)! PBL height (m)
    KPBL    , & !INTENT(OUT  ) :: KPBL   (1:iMax)!-- KPBL layer index of the PBL
                    PBL_CoefKm,&
            PBL_CoefKh,&
                    RUBLTEN , & !INTENT(OUT  ) :: RUBLTEN  (1:iMax,1:kMax)! U tendency due to PBL parameterization (m/s^2)
    RVBLTEN , & !INTENT(OUT  ) :: RVBLTEN  (1:iMax,1:kMax)! V tendency due to PBL parameterization (m/s^2)
    RTHBLTEN, & !INTENT(OUT  ) :: RTHBLTEN (1:iMax,1:kMax)! Theta tendency due to PBL parameterization (K/s)
    RQVBLTEN, & !INTENT(OUT  ) :: RQVBLTEN (1:iMax,1:kMax)! Qv tendency due to PBL parameterization (kg/kg/s)
    RQCBLTEN  ) !INTENT(OUT  ) :: RQCBLTEN (1:iMax,1:kMax)! Qc tendency due to PBL parameterization (kg/kg/s)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    INTEGER      ,INTENT(IN) :: iMax
    INTEGER      ,INTENT(IN) :: kMax
    REAL(KIND=r8),INTENT(IN) :: DT    ! time step (second)
    INTEGER      ,INTENT(IN) :: STEPBL! bldt (max_dom)= 0,; minutes between boundary-layer physics calls
                                      !-- calculate pbl time step
                                      !   STEPBL = nint(BLDT*60./DT)
                                      !   STEPBL = max(STEPBL,1)
    REAL(KIND=r8),INTENT(IN) :: HT  (1:iMax)!"HGT" "Terrain Height"   "m"
    REAL(KIND=r8),INTENT(IN) :: DZ  (1:iMax,1:kMax)!-- dz8w dz between full levels (m)
    REAL(KIND=r8),INTENT(IN) :: PMID(1:iMax,1:kMax)!-- p_phy pressure (Pa)
    REAL(KIND=r8),INTENT(IN) :: PINT(1:iMax,1:kMax)!-- p8w pressure at full levels (Pa)
    REAL(KIND=r8),INTENT(IN) :: TH  (1:iMax,1:kMax)! potential temperature (K)
    REAL(KIND=r8),INTENT(IN) :: T   (1:iMax,1:kMax)!t_phy         temperature (K)
    REAL(KIND=r8),INTENT(IN) :: QV   (1:iMax,1:kMax)! Qv         water vapor mixing ratio (kg/kg)
    REAL(KIND=r8),INTENT(IN) :: CWM  (1:iMax,1:kMax)! cloud water mixing ratio (kg/kg)
    REAL(KIND=r8),INTENT(IN) :: U    (1:iMax,1:kMax)! u-velocity interpolated to theta points (m/s)
    REAL(KIND=r8),INTENT(IN) :: V    (1:iMax,1:kMax)! v-velocity interpolated to theta points (m/s)
    REAL(KIND=r8),INTENT(IN) :: TSK  (1:iMax) !-- TSK           surface temperature (K)
    REAL(KIND=r8),INTENT(IN) :: CHKLOWQ(1:iMax)!-- CHKLOWQ - is either 0 or 1 (so far set equal to 1).
                                               !-- used only in MYJPBL. 
                                               ! For RUC LSM CHKLOWQ needed for MYJPBL should 
                                               ! 1 because is actual specific humidity at the surface, and
                                               ! not the saturation value
    INTEGER      ,INTENT(IN) :: LOWLYR  (1:iMax)!-- lowlyr index of lowest model layer above ground
    REAL(KIND=r8),INTENT(IN) :: XLAND(1:iMax) !-- XLAND         land mask (1 for land, 2 for water)
    REAL(KIND=r8),INTENT(IN) :: SICE (1:iMax)
    REAL(KIND=r8),INTENT(IN) :: SNOW (1:iMax) !-- SNOW        liquid water-equivalent snow depth (m)
    REAL(KIND=r8),INTENT(IN) :: ELFLX  (1:iMax)!-- ELFLX--LH net upward latent heat flux at surface (W/m^2)
    REAL(KIND=r8),INTENT(IN) :: bps    (1:iMax,1:kMax)  
    REAL(KIND=r8),INTENT(IN) :: colrad (1:iMax) 


    REAL(KIND=r8),INTENT(INOUT) :: QSFC   (1:iMax)!-- qsfc          specific humidity at lower boundary (kg/kg)
    REAL(KIND=r8),INTENT(INOUT) :: THZ0   (1:iMax)!-- thz0          potential temperature at roughness length (K)
    REAL(KIND=r8),INTENT(INOUT) :: QZ0    (1:iMax)!-- QZ0           specific humidity at roughness length (kg/kg)
    REAL(KIND=r8),INTENT(INOUT) :: UZ0    (1:iMax)!-- uz0       u wind component at roughness length (m/s)
    REAL(KIND=r8),INTENT(INOUT) :: VZ0    (1:iMax)!-- vz0       v wind component at roughness length (m/s)
    REAL(KIND=r8),INTENT(INOUT) :: TKE_MYJ(1:iMax,1:kMax)!-- tke_myj       turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
    REAL(KIND=r8),INTENT(INOUT) :: EXCH_H (1:iMax,1:kMax)
    REAL(KIND=r8),INTENT(INOUT) :: USTAR  (1:iMax)!-- UST           u* in similarity theory (m/s)
    REAL(KIND=r8),INTENT(INOUT) :: CT     (1:iMax)
    REAL(KIND=r8),INTENT(INOUT) :: AKHS   (1:iMax)!-- akhs   sfc exchange coefficient of heat/moisture from MYJ
    REAL(KIND=r8),INTENT(INOUT) :: AKMS   (1:iMax)!-- akms   sfc exchange coefficient of momentum from MYJ
    REAL(KIND=r8),    INTENT(INOUT) :: PBL_CoefKm(iMax, kmax)
    REAL(KIND=r8),    INTENT(INOUT) :: PBL_CoefKh(iMax, kmax)

    !
    REAL(KIND=r8),INTENT(OUT  ) :: PBLH(1:iMax)! PBL height (m)
    INTEGER      ,INTENT(OUT  ) :: KPBL   (1:iMax)!-- KPBL layer index of the PBL
    !
    REAL(KIND=r8),INTENT(OUT  ) :: EL_MYJ   (1:iMax,1:kMax)! mixing length from Mellor-Yamada-Janjic (MYJ) (m)
    REAL(KIND=r8),INTENT(OUT  ) :: RQCBLTEN (1:iMax,1:kMax)! Qc tendency due to PBL parameterization (kg/kg/s)
    REAL(KIND=r8),INTENT(OUT  ) :: RQVBLTEN (1:iMax,1:kMax)! Qv tendency due to PBL parameterization (kg/kg/s)
    REAL(KIND=r8),INTENT(OUT  ) :: RTHBLTEN (1:iMax,1:kMax)! Theta tendency due to PBL parameterization (K/s)
    REAL(KIND=r8),INTENT(OUT  ) :: RUBLTEN  (1:iMax,1:kMax)! U tendency due to PBL parameterization (m/s^2)
    REAL(KIND=r8),INTENT(OUT  ) :: RVBLTEN  (1:iMax,1:kMax)! V tendency due to PBL parameterization (m/s^2)
    !
    !
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: I,K,KFLIP,LLOW,LMH,LMXL
    !
    INTEGER,DIMENSION(1:iMax) :: LPBL
    !
    REAL(KIND=r8) :: APEX,DCDT,DELTAZ,DQDT,DTDIF,DTDT,DTTURBL,DUDT,DVDT       &
         &       ,FFSK,PLOW,PSFC,PTOP,QFC1,QLOW,QOLD                &
         &       ,RATIOMX,RDTTURBL,RG,SEAMASK,THNEW,THOLD,TX  
    !
    REAL(KIND=r8),DIMENSION(1:kMax) :: CWMK,PK,Q2K,QK,THEK,TK,UK,VK
    !
    REAL(KIND=r8),DIMENSION(1:kMax-1) :: AKHK,AKMK,EL,GH,GM
    !
    REAL(KIND=r8),DIMENSION(1:kMax+1) :: ZHK
    !
    REAL(KIND=r8),DIMENSION(1:iMax) :: THSK
    !
    REAL(KIND=r8),DIMENSION(1:iMax,1:kMax) :: RHOD
    !
    REAL(KIND=r8),DIMENSION(1:iMax,1:kMax) :: APE,THE
    !
    REAL(KIND=r8),DIMENSION(1:iMax,1:kMax-1) :: AKH,AKM
    !
    REAL(KIND=r8),DIMENSION(1:iMax,1:kMax+1) :: ZINT
    !
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    !
    !***  MAKE PREPARATIONS
    !
    !----------------------------------------------------------------------
    DTTURBL  = DT*STEPBL
    RDTTURBL = 1.0_r8/DTTURBL
    DTDIF    = DTTURBL
    RG       = 1.0_r8/G
    !
    !----------------------------------------------------------------------
    EL_MYJ = 0.0_r8
    !setup_integration:  DO J=jMax0,jMax
    !----------------------------------------------------------------------
    !

    DO K=1,kMax-1
       DO I=1,iMax
          AKM(I,K)=0.0_r8
       ENDDO
    ENDDO
    !
    DO K=1,kMax+1
       DO I=1,iMax
          ZINT(I,K)=0.0_r8
       ENDDO
    ENDDO
    !
    DO I=1,iMax
       ZINT(I,kMax+1)=HT(I)     ! Z at bottom of lowest sigma layer
       !
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!!!!
!!!!!!  ZINT(I,kMax+1)=1.E-4_r8       ! Z of bottom of lowest eta layer
!!!!!!  ZHK(kMax+1)=1.E-4_r8            ! Z of bottom of lowest eta layer
       !
    ENDDO
    !
    DO K=kMax,1,-1
       KFLIP=kMax+1-K
       DO I=1,iMax
          ZINT(I,K)=ZINT(I,K+1)+DZ(I,KFLIP)
          APEX=bps(i,k)             ! 1.0_r8/EXNER(I,K)
          APE(I,K)=APEX
          TX=T(I,K)
          THE(I,K)=(CWM(I,K)*(-ELOCP/TX) + 1.0_r8)*TH(I,K)
       ENDDO
    ENDDO
    !
    EL_MYJ(1:iMax,1:kMax) =0

    DO I=1,iMax
       !
       !***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
       !
       LMH=kMax-LOWLYR(I)+1
       !
       PTOP=PINT(I,kMax)      ! kMax+1=kMax
       PSFC=PINT(I,LOWLYR(I))
       !
       !***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
       !
       SEAMASK=XLAND(I)-1.0_r8
       !
       !***  FILL 1-D VERTICAL ARRAYS
       !***  AND FLIP DIRECTION SINCE MYJ SCHEME
       !***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
       !
       DO K=kMax,1,-1
          KFLIP=kMax+1-K
          TK(K)=T(I,KFLIP)
          THEK(K)=THE(I,KFLIP)
          RATIOMX=QV(I,KFLIP)
          QK(K)=RATIOMX/(1.0_r8+RATIOMX)
          CWMK(K)=CWM(I,KFLIP)
          PK(K)=PMID(I,KFLIP)
          UK(K)=U(I,KFLIP)
          VK(K)=V(I,KFLIP)
          !
          !***  TKE=0.5_r8*(q**2) ==> q**2=2.0_r8*TKE
          !
          Q2K(K)=2.0_r8*TKE_MYJ(I,KFLIP)
          !
          !***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
          !
          ZHK(K)=ZINT(I,K)
          !
       ENDDO
       ZHK(kMax+1)=HT(I)          ! Z at bottom of lowest sigma layer
       !----------------------------------------------------------------------
       !***
       !***  FIND THE MIXING LENGTH
       !***
       CALL MIXLEN(LMH,UK,VK,TK,THEK,QK,CWMK                        &
            &               ,Q2K,ZHK,GM,GH,EL                                 &
            &               ,PBLH(I),LPBL(I),LMXL,CT(I)                 &
            &               ,kMax)

       !
       !----------------------------------------------------------------------
       !***
       !***  SOLVE FOR THE PRODUCTION/DISSIPATION OF
       !***  THE TURBULENT KINETIC ENERGY
       !***
       !
       CALL PRODQ2(LMH,DTTURBL,USTAR(I),GM,GH,EL,Q2K              &
            &                ,kMax)

       !
       !----------------------------------------------------------------------
       !*** THE MODEL LAYER (COUNTING UPWARD) CONTAINING THE TOP OF THE PBL
       !----------------------------------------------------------------------
       !
       KPBL(I)=kMax-LPBL(I)+1
       !
       !----------------------------------------------------------------------
       !***
       !***  FIND THE EXCHANGE COEFFICIENTS IN THE FREE ATMOSPHERE
       !***
       CALL DIFCOF(LMH,GM,GH,EL,Q2K,ZHK,AKMK,AKHK           &
            &,kMax)
       !
       !***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH 
       !***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS 1 TO kMax-1.  COUNTING 
       !***  COUNTING UPWARD FROM THE BOTTOM, THOSE SAME COEFFICIENTS EXCH_H
       !***  ARE DEFINED ON THE TOPS OF THE LAYERS 1 TO kMax-1.
       !
       DO K=1,kMax-1
          KFLIP=kMax-K
          AKH(I,K)=AKHK(K)
          AKM(I,K)=AKMK(K)
          PBL_CoefKm(I, k)  =AKHK(KFLIP)
          PBL_CoefKh(I, k)  =AKMK(KFLIP)
          DELTAZ=0.5_r8*(ZHK(KFLIP)-ZHK(KFLIP+2))
          EXCH_H(I,K)=AKHK(KFLIP)*DELTAZ
       ENDDO
       PBL_CoefKm(I, kMax)=PBL_CoefKm(I, kMax-1) 
       PBL_CoefKh(I, kMax)=PBL_CoefKh(I, kMax-1) 


       !
       !----------------------------------------------------------------------
       !***
       !***  CARRY OUT THE VERTICAL DIFFUSION OF
       !***  TURBULENT KINETIC ENERGY
       !***
       !
       CALL VDIFQ(LMH,DTDIF,Q2K,EL,ZHK,kMax)
       !
       !***  SAVE THE NEW TKE AND MIXING LENGTH.
       !
       DO K=1,kMax
          KFLIP=kMax+1-K
          Q2K(KFLIP)=MAX(Q2K(KFLIP),EPSQ2)
          TKE_MYJ(I,K)=0.5_r8*Q2K(KFLIP)
          IF ( K .LT. kMax ) EL_MYJ(I,K)=EL(K)   ! EL IS NOT DEFINED AT kMax
       ENDDO
       !
       !
    ENDDO
    !
    !----------------------------------------------------------------------
    !ENDDO setup_integration
    !----------------------------------------------------------------------
    !
    !***  CONVERT SURFACE SENSIBLE TEMPERATURE TO POTENTIAL TEMPERATURE.
    !
    DO I=1,iMax
       PSFC=PINT(I,LOWLYR(I))
       THSK(I)=TSK(I)*bps(i,1)!(1.E5_r8/PSFC)**CAPA
    ENDDO
    !ENDDO
    !
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !
    DO I=1,iMax
       !
       !***  FILL 1-D VERTICAL ARRAYS
       !***  AND FLIP DIRECTION SINCE MYJ SCHEME
       !***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
       !
       DO K=kMax-1,1,-1
          KFLIP=kMax-K
          AKHK(K)=AKH(I,K)
       ENDDO
       !
       DO K=kMax,1,-1
          KFLIP=kMax+1-K
          THEK(K)=THE(I,KFLIP)
          RATIOMX=QV(I,KFLIP)
          QK(K)=RATIOMX/(1.0_r8+RATIOMX)
          CWMK(K)=CWM(I,KFLIP)
          ZHK(K)=ZINT(I,K)
       ENDDO
       !
       ZHK(kMax+1)=ZINT(I,kMax+1)
       !
       SEAMASK=XLAND(I)-1.0_r8
       !THZ0(I)=(1.0_r8-SEAMASK)*THSK(I)+SEAMASK*THZ0(I)
       !!!!!!!
       THZ0(I)=THSK(I)
       LLOW=LOWLYR(I)
       PLOW=PMID(I,LLOW)
       QLOW=QK(kMax+1-LLOW)
       FFSK=AKHS(I)*PLOW/((QLOW*P608-CWM(I,LLOW)+1.0_r8)*T(I,LLOW)*R_D)
       QFC1=FFSK*XLV
       QFC1=QFC1*CHKLOWQ(I)
       !
       IF(SNOW(I)>0.0_r8.OR.SICE(I)>0.5_r8)THEN
          QFC1=QFC1*RLIVWV
       ENDIF
       !
       IF(QFC1>0.0_r8)THEN
          QSFC(I)=QLOW+ELFLX(I)/QFC1
       ENDIF
       
       QZ0 (I)=QSFC(I)

       !!*** SEAMASK (1 FOR SEA; 0 FOR LAND)

       !IF(SEAMASK<0.5_r8)THEN
       !   LLOW=LOWLYR(I)
       !   PLOW=PMID(I,LLOW)
       !   QLOW=QK(kMax+1-LLOW)
       !   FFSK=AKHS(I)*PLOW/((QLOW*P608-CWM(I,LLOW)+1.0_r8)*T(I,LLOW)*R_D)
       !   QFC1=FFSK*XLV
       !   QFC1=QFC1*CHKLOWQ(I)
       !   !
       !   IF(SNOW(I)>0.0_r8.OR.SICE(I)>0.5_r8)THEN
       !      QFC1=QFC1*RLIVWV
       !   ENDIF
       !   !
       !   IF(QFC1>0.0_r8)THEN
       !      QSFC(I)=QLOW+ELFLX(I)/QFC1
       !   ENDIF
       !   !
       !ELSE
       !   PSFC=PINT(I,LOWLYR(I))
       !   EXNSFC=bps(i,1)!(1.E5_r8/PSFC)**CAPA
       !  !
       !  !
       !   QSFC (I) = PQ0SEA / PSFC * EXP(A2*(THSK(I) - A3*EXNSFC)   / (THSK(I)-A4*EXNSFC))
       !ENDIF
       !
       !QZ0 (I)=(1.0_r8-SEAMASK)*QSFC(I)+SEAMASK*QZ0 (I)
       !
       !***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
       !
       LMH=kMax-LOWLYR(I)+1
       !
       !----------------------------------------------------------------------
       !***  CARRY OUT THE VERTICAL DIFFUSION OF
       !***  TEMPERATURE AND WATER VAPOR
       !----------------------------------------------------------------------
       !
       CALL VDIFH(DTDIF,LMH,THZ0(I),QZ0(I)                      &
            &              ,AKHS(I),CHKLOWQ(I),CT(I)                    &
            &              ,THEK,QK,CWMK,AKHK,ZHK                             &
            &              ,kMax)

       !----------------------------------------------------------------------
       !***
       !***  COMPUTE PRIMARY VARIABLE TENDENCIES
       !***
       DO K=1,kMax
          KFLIP=kMax+1-K
          THOLD=TH(I,K)
          THNEW=THEK(KFLIP)+CWMK(KFLIP)*ELOCP*APE(I,K)
          DTDT=((THNEW)-(THOLD))*RDTTURBL
          QOLD=QV(I,K)/(1.0_r8+QV(I,K))
          DQDT=(QK(KFLIP) - QOLD)*RDTTURBL
          DCDT=(CWMK(KFLIP)-CWM(I,K))*RDTTURBL
          !
          RHOD(I,K)=PMID(I,K)/(R_D*T(I,K))
          RTHBLTEN(I,K)=DTDT
          RQVBLTEN(I,K)=DQDT/(1.0_r8-QK(KFLIP))**2
          RQCBLTEN(I,K)=DCDT
       ENDDO
       !----------------------------------------------------------------------
    ENDDO
    !----------------------------------------------------------------------
    DO I=1,iMax
       !
       !***  FILL 1-D VERTICAL ARRAYS
       !***  AND FLIP DIRECTION SINCE MYJ SCHEME
       !***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
       !
       DO K=kMax-1,1,-1
          AKMK(K)=AKM(I,K)
       ENDDO
       !
       DO K=kMax,1,-1
          KFLIP=kMax+1-K
          UK(K)=U(I,KFLIP)
          VK(K)=V(I,KFLIP)
          ZHK(K)=ZINT(I,K)
       ENDDO
       ZHK(kMax+1)=ZINT(I,kMax+1)
       !
       !----------------------------------------------------------------------
       !***  CARRY OUT THE VERTICAL DIFFUSION OF
       !***  VELOCITY COMPONENTS
       !----------------------------------------------------------------------
       !
       CALL VDIFV(LMH,DTDIF,UZ0(I),VZ0(I)                       &
            &              ,AKMS(I),UK,VK,AKMK,ZHK                          &
            &              ,kMax)
       !
       !----------------------------------------------------------------------
       !***
       !***  COMPUTE PRIMARY VARIABLE TENDENCIES
       !***
       DO K=1,kMax
          KFLIP=kMax+1-K
          DUDT=(UK(KFLIP)*SIN( colrad(i)) - U(I,K)*SIN( colrad(i)))*RDTTURBL
          DVDT=(VK(KFLIP)*SIN( colrad(i)) - V(I,K)*SIN( colrad(i)))*RDTTURBL
          RUBLTEN(I,K)=DUDT
          RVBLTEN(I,K)=DVDT
       ENDDO
       !
    ENDDO
    !----------------------------------------------------------------------
    !
    !ENDDO setup_integration
    !
    !----------------------------------------------------------------------
    !
  END SUBROUTINE MellorYamada1
  !

  !
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE MIXLEN                            &
       !----------------------------------------------------------------------
       !   ******************************************************************
       !   *                                                                *
       !   *                   LEVEL 2.5 MIXING LENGTH                      *
       !   *                                                                *
       !   ******************************************************************
       !
    &(LMH,U,V,T,THE,Q,CWM,Q2,Z,GM,GH,EL,PBLH,LPBL,LMXL,CT             &
         &,kMax)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    INTEGER,INTENT(IN) :: kMax
    !
    INTEGER,INTENT(IN) :: LMH
    !
    INTEGER,INTENT(OUT) :: LPBL
    !
    REAL(KIND=r8),INTENT(IN)  :: CWM  (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: Q    (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: Q2   (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: T    (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: THE  (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: U    (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: V    (1:kMax)
    REAL(KIND=r8),INTENT(IN)  :: Z    (1:kMax+1)
    REAL(KIND=r8),INTENT(OUT) :: PBLH 
    REAL(KIND=r8),INTENT(OUT) :: EL   (1:kMax-1)
    REAL(KIND=r8),INTENT(OUT) :: GH   (1:kMax-1)
    REAL(KIND=r8),INTENT(OUT) :: GM   (1:kMax-1)
    !
    REAL(KIND=r8),INTENT(INOUT) :: CT
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: K,LMXL,LPBLM
    !
    REAL(KIND=r8) :: A,ADEN,B,BDEN,AUBR,BLMX,BUBR,EL0,ELOQ2X,GHL,GML          &
         &       ,QOL2ST,QOL2UN,QDZL,RDZ,SQ,SREL,SZQ,TEM,THM,VKRMZ
    !
    REAL(KIND=r8),DIMENSION(1:kMax) :: Q1
    !
    REAL(KIND=r8),DIMENSION(1:kMax-1) :: DTH,ELM,REL
    !
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    !--------------FIND THE HEIGHT OF THE PBL-------------------------------
    !----------------------------------------------------------------------
    LPBL=LMH
    DTH=0.0_r8
    !
    DO K=LMH-1,1,-1
       IF(Q2(K)<=EPSQ2*FH)THEN
          LPBL=K
          GO TO 110
       ENDIF
    ENDDO
    !
    LPBL=1
    !----------------------------------------------------------------------
    !--------------THE HEIGHT OF THE PBL------------------------------------
    !----------------------------------------------------------------------
110 PBLH=Z(LPBL)-Z(LMH+1)
    !-----------------------------------------------------------------------
    !
    DO K=1,LMH
       Q1(K)=0.0_r8
    ENDDO
    !
    DO K=1,LMH-1
       DTH(K)=THE(K)-THE(K+1)
    ENDDO
    !
    DO K=LMH-2,1,-1
       IF(DTH(K)>0.0_r8.AND.DTH(K+1)<=0.0_r8)THEN
          DTH(K)=DTH(K)+CT
          EXIT
       ENDIF
    ENDDO
    !
    CT=0.0_r8
    !----------------------------------------------------------------------
    DO K=1,LMH-1
       RDZ=2.0_r8/(Z(K)-Z(K+2))
       GML=((U(K)-U(K+1))**2+(V(K)-V(K+1))**2)*RDZ*RDZ
       GM(K)=MAX(GML,EPSGM)
       !
       TEM=(T(K)+T(K+1))*0.5_r8
       THM=(THE(K)+THE(K+1))*0.5_r8
       !
       A=THM*P608
       B=(ELOCP/TEM-1.0_r8-P608)*THM
       !
       GHL=(DTH(K)*((Q(K)+Q(K+1)+CWM(K)+CWM(K+1))*(0.5_r8*P608)+1.0_r8)      &
            &     +(Q(K)-Q(K+1)+CWM(K)-CWM(K+1))*A                            &
            &     +(CWM(K)-CWM(K+1))*B)*RDZ
       !
       IF(ABS(GHL)<=EPSGH)GHL=EPSGH
       GH(K)=GHL
    ENDDO
    !
    !----------------------------------------------------------------------
    !***  FIND MAXIMUM MIXING LENGTHS AND THE LEVEL OF THE PBL TOP
    !----------------------------------------------------------------------
    !
    LMXL=LMH
    !
    DO K=1,LMH-1
       GML=GM(K)
       GHL=GH(K)
       !
       IF(GHL>=EPSGH)THEN
          IF(GML/GHL<=REQU)THEN
             ELM(K)=EPSL
             LMXL=K
          ELSE
             AUBR=(AUBM*GML+AUBH*GHL)*GHL
             BUBR= BUBM*GML+BUBH*GHL
             QOL2ST=(-0.5_r8*BUBR+SQRT(BUBR*BUBR*0.25_r8-AUBR*CUBR))*RCUBR
             ELOQ2X=1.0_r8/QOL2ST
             ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
          ENDIF
       ELSE
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          QOL2UN=-0.5_r8*BDEN+SQRT(BDEN*BDEN*0.25_r8-ADEN)
          ELOQ2X=1.0_r8/(QOL2UN+EPSRU)       ! repsr1/qol2un
          ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
       ENDIF
    ENDDO
    !
    IF(ELM(LMH-1)==EPSL)LMXL=LMH
    !
    !----------------------------------------------------------------------
    !***  THE HEIGHT OF THE MIXED LAYER
    !----------------------------------------------------------------------
    !
    BLMX=Z(LMXL)-Z(LMH+1)
    !
    !----------------------------------------------------------------------
    DO K=LPBL,LMH
       Q1(K)=SQRT(Q2(K))
    ENDDO
    !----------------------------------------------------------------------
    SZQ=0.0_r8
    SQ =0.0_r8
    !
    DO K=1,LMH-1
       QDZL=(Q1(K)+Q1(K+1))*(Z(K+1)-Z(K+2))
       SZQ=(Z(K+1)+Z(K+2)-Z(LMH+1)-Z(LMH+1))*QDZL+SZQ
       SQ=QDZL+SQ
    ENDDO
    !
    !----------------------------------------------------------------------
    !***  COMPUTATION OF ASYMPTOTIC L IN BLACKADAR FORMULA
    !----------------------------------------------------------------------
    !
    EL0=MIN(ALPH*SZQ*0.5_r8/SQ,EL0MAX)
    EL0=MAX(EL0            ,EL0MIN)
    !
    !----------------------------------------------------------------------
    !***  ABOVE THE PBL TOP
    !----------------------------------------------------------------------
    !
    LPBLM=MAX(LPBL-1,1)
    DO K=1,LPBLM
       EL(K)=MIN((Z(K)-Z(K+2))*ELFC,ELM(K))
       REL(K)=EL(K)/ELM(K)
    ENDDO
    !
    !----------------------------------------------------------------------
    !***  INSIDE THE PBL
    !----------------------------------------------------------------------
    !
    IF(LPBL<LMH)THEN
       DO K=LPBL,LMH-1
          VKRMZ=(Z(K+1)-Z(LMH+1))*VKARMAN
          EL(K)=MIN(VKRMZ/(VKRMZ/EL0+1.0_r8),ELM(K))
          REL(K)=EL(K)/ELM(K)
       ENDDO
    ENDIF
    !
    DO K=LPBL+1,LMH-2
       SREL=MIN(((REL(K-1)+REL(K+1))*0.5_r8+REL(K))*0.5_r8,REL(K))
       EL(K)=MAX(SREL*ELM(K),EPSL)
    ENDDO
    !
    !----------------------------------------------------------------------
  END SUBROUTINE MIXLEN
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE PRODQ2                            &
       !----------------------------------------------------------------------
       !   ******************************************************************
       !   *                                                                *
       !   *            LEVEL 2.5 Q2 PRODUCTION/DISSIPATION                 *
       !   *                                                                *
       !   ******************************************************************
       !
    &(LMH,DTTURBL,USTAR,GM,GH,EL,Q2 &
         &,kMax)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    INTEGER,INTENT(IN) :: kMax
    !
    INTEGER,INTENT(IN) :: LMH
    !
    REAL(KIND=r8)   ,INTENT(IN) :: DTTURBL
    REAL(KIND=r8)   ,INTENT(IN) :: USTAR
    !
    REAL(KIND=r8),INTENT(IN)   :: GH   (1:kMax-1)
    REAL(KIND=r8),INTENT(IN)   :: GM   (1:kMax-1)
    REAL(KIND=r8),INTENT(INOUT)   :: EL   (1:kMax-1)
    !
    REAL(KIND=r8),INTENT(INOUT)   :: Q2   (1:kMax)  
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: K
    !
    REAL(KIND=r8) :: ADEN,AEQU,ANUM,ARHS,BDEN,BEQU,BNUM,BRHS,CDEN,CRHS        &
         &       ,DLOQ1,ELOQ11,ELOQ12,ELOQ13,ELOQ21,ELOQ22,ELOQ31,ELOQ32   &
         &       ,ELOQ41,ELOQ42,ELOQ51,ELOQ52,ELOQN,EQOL2,GHL,GML          &
         &       ,RDEN1,RDEN2,RHS2,RHSP1,RHSP2,RHST2
    !
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    !
    main_integration: DO K=1,LMH-1
       GML=GM(K)
       GHL=GH(K)
       !
       !----------------------------------------------------------------------
       !***  COEFFICIENTS OF THE EQUILIBRIUM EQUATION
       !----------------------------------------------------------------------
       !
       AEQU=(AEQM*GML+AEQH*GHL)*GHL
       BEQU= BEQM*GML+BEQH*GHL
       !
       !----------------------------------------------------------------------
       !***  EQUILIBRIUM SOLUTION FOR L/Q
       !----------------------------------------------------------------------
       !
       EQOL2=-0.5_r8*BEQU+SQRT(BEQU*BEQU*0.25_r8-AEQU)
       !
       !----------------------------------------------------------------------
       !***  IS THERE PRODUCTION/DISSIPATION ?
       !----------------------------------------------------------------------
       !
       IF((GML+GHL*GHL<=EPSTRB)                                       &
            &   .OR.(GHL>=EPSGH.AND.GML/GHL<=REQU)                            &
            &   .OR.(EQOL2<=EPS2))THEN
          !
          !----------------------------------------------------------------------
          !***  NO TURBULENCE
          !----------------------------------------------------------------------
          !
          Q2(K)=EPSQ2
          EL(K)=EPSL
          !----------------------------------------------------------------------
          !
       ELSE
          !
          !----------------------------------------------------------------------
          !***  TURBULENCE
          !----------------------------------------------------------------------
          !----------------------------------------------------------------------
          !***  COEFFICIENTS OF THE TERMS IN THE NUMERATOR
          !----------------------------------------------------------------------
          !
          ANUM=(ANMM*GML+ANMH*GHL)*GHL
          BNUM= BNMM*GML+BNMH*GHL
          !
          !----------------------------------------------------------------------
          !***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
          !----------------------------------------------------------------------
          !
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          CDEN= 1.0_r8
          !
          !----------------------------------------------------------------------
          !***  COEFFICIENTS OF THE NUMERATOR OF THE LINEARIZED EQ.
          !----------------------------------------------------------------------
          !
          ARHS=-(ANUM*BDEN-BNUM*ADEN)*2.0_r8
          BRHS=- ANUM*4.0_r8
          CRHS=- BNUM*2.0_r8
          !
          !----------------------------------------------------------------------
          !***  INITIAL VALUE OF L/Q
          !----------------------------------------------------------------------
          !
          DLOQ1=EL(K)/SQRT(Q2(K))
          !
          !----------------------------------------------------------------------
          !***  FIRST ITERATION FOR L/Q, RHS=0
          !----------------------------------------------------------------------
          !
          ELOQ21=1.0_r8/EQOL2
          ELOQ11=SQRT(ELOQ21)
          ELOQ31=ELOQ21*ELOQ11
          ELOQ41=ELOQ21*ELOQ21
          ELOQ51=ELOQ21*ELOQ31
          !
          !----------------------------------------------------------------------
          !***  1./DENOMINATOR
          !----------------------------------------------------------------------
          !
          RDEN1=1.0_r8/(ADEN*ELOQ41+BDEN*ELOQ21+CDEN)
          !
          !----------------------------------------------------------------------
          !***  D(RHS)/D(L/Q)
          !----------------------------------------------------------------------
          !
          RHSP1=(ARHS*ELOQ51+BRHS*ELOQ31+CRHS*ELOQ11)*RDEN1*RDEN1
          !
          !----------------------------------------------------------------------
          !***  FIRST-GUESS SOLUTION
          !----------------------------------------------------------------------
          !
          ELOQ12=ELOQ11+(DLOQ1-ELOQ11)*EXP(RHSP1*DTTURBL)
          ELOQ12=MAX(ELOQ12,EPS1)
          !
          !----------------------------------------------------------------------
          !***  SECOND ITERATION FOR L/Q
          !----------------------------------------------------------------------
          !
          ELOQ22=ELOQ12*ELOQ12
          ELOQ32=ELOQ22*ELOQ12
          ELOQ42=ELOQ22*ELOQ22
          ELOQ52=ELOQ22*ELOQ32
          !
          !----------------------------------------------------------------------
          !***  1./DENOMINATOR
          !----------------------------------------------------------------------
          !
          RDEN2=1.0_r8/(ADEN*ELOQ42+BDEN*ELOQ22+CDEN)
          RHS2 =-(ANUM*ELOQ42+BNUM*ELOQ22)*RDEN2+RB1
          RHSP2= (ARHS*ELOQ52+BRHS*ELOQ32+CRHS*ELOQ12)*RDEN2*RDEN2
          RHST2=RHS2/RHSP2
          !
          !----------------------------------------------------------------------
          !***  CORRECTED SOLUTION
          !----------------------------------------------------------------------
          !
          ELOQ13=ELOQ12-RHST2+(RHST2+DLOQ1-ELOQ12)*EXP(RHSP2*DTTURBL)
          ELOQ13=MAX(ELOQ13,EPS1)
          !
          !----------------------------------------------------------------------
          !***  TWO ITERATIONS IS ENOUGH IN MOST CASES ...
          !----------------------------------------------------------------------
          !
          ELOQN=ELOQ13
          !
          IF(ELOQN>EPS1)THEN
             Q2(K)=EL(K)*EL(K)/(ELOQN*ELOQN)
             Q2(K)=MAX(Q2(K),EPSQ2)
             !
             IF(Q2(K)==EPSQ2)THEN
                EL(K)=EPSL
             ENDIF
             !
          ELSE
             Q2(K)=EPSQ2
             EL(K)=EPSL
          ENDIF
          !
          !----------------------------------------------------------------------
          !***  END OF TURBULENT BRANCH
          !----------------------------------------------------------------------
          !
       ENDIF
       !----------------------------------------------------------------------
       !***  END OF PRODUCTION/DISSIPATION LOOP
       !----------------------------------------------------------------------
       !
    ENDDO main_integration
    !
    !----------------------------------------------------------------------
    !***  LOWER BOUNDARY CONDITION FOR Q2
    !----------------------------------------------------------------------
    !
    Q2(LMH)=MAX(B1**(2.0_r8/3.0_r8)*USTAR*USTAR,EPSQ2)
    !----------------------------------------------------------------------
    !
  END SUBROUTINE PRODQ2
  !
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE DIFCOF                           &
       !   ******************************************************************
       !   *                                                                *
       !   *                LEVEL 2.5 DIFFUSION COEFFICIENTS                *
       !   *                                                                *
       !   ******************************************************************
    &(LMH,GM,GH,EL,Q2,Z,AKM,AKH                                &
         &,kMax)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    INTEGER,INTENT(IN) :: kMax
    ! 
    INTEGER,INTENT(IN) :: LMH
    !
    REAL(KIND=r8),INTENT(IN) :: Q2  (1:kMax)  
    REAL(KIND=r8),INTENT(IN) :: EL  (1:kMax-1)
    REAL(KIND=r8),INTENT(IN) :: GH  (1:kMax-1)
    REAL(KIND=r8),INTENT(IN) :: GM  (1:kMax-1)
    REAL(KIND=r8),INTENT(IN) :: Z   (1:kMax+1)
    !
    REAL(KIND=r8),INTENT(OUT) :: AKH(1:kMax-1)
    REAL(KIND=r8),INTENT(OUT) :: AKM(1:kMax-1)
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: K
    !
    REAL(KIND=r8) :: ADEN,BDEN,BESH,BESM,CDEN,ELL,ELOQ2,ELOQ4,ELQDZ &
         &       ,ESH,ESM,GHL,GML,Q1L,RDEN,RDZ
    !
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    !
    DO K=1,LMH-1
       ELL=EL(K)
       !
       ELOQ2=ELL*ELL/Q2(K)
       ELOQ4=ELOQ2*ELOQ2
       !
       GML=GM(K)
       GHL=GH(K)
       !
       !----------------------------------------------------------------------
       !***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
       !----------------------------------------------------------------------
       !
       ADEN=(ADNM*GML+ADNH*GHL)*GHL
       BDEN= BDNM*GML+BDNH*GHL
       CDEN= 1.0_r8
       !
       !----------------------------------------------------------------------
       !***  COEFFICIENTS FOR THE SM DETERMINANT
       !----------------------------------------------------------------------
       !
       BESM=BSMH*GHL
       !
       !----------------------------------------------------------------------
       !***  COEFFICIENTS FOR THE SH DETERMINANT
       !----------------------------------------------------------------------
       !
       BESH=BSHM*GML+BSHH*GHL
       !
       !----------------------------------------------------------------------
       !***  1./DENOMINATOR
       !----------------------------------------------------------------------
       !
       RDEN=1.0_r8/(ADEN*ELOQ4+BDEN*ELOQ2+CDEN)
       !
       !----------------------------------------------------------------------
       !***  SM AND SH
       !----------------------------------------------------------------------
       !
       ESM=(BESM*ELOQ2+CESM)*RDEN
       ESH=(BESH*ELOQ2+CESH)*RDEN
       !
       !----------------------------------------------------------------------
       !***  DIFFUSION COEFFICIENTS
       !----------------------------------------------------------------------
       !
       RDZ=2.0_r8/(Z(K)-Z(K+2))
       Q1L=SQRT(Q2(K))
       ELQDZ=ELL*Q1L*RDZ
       AKM(K)=ELQDZ*ESM
       AKH(K)=ELQDZ*ESH
       !----------------------------------------------------------------------
    ENDDO
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    !***  INVERSIONS
    !----------------------------------------------------------------------
    !
    !     IF(LMXL==LMH-1)THEN
    !
    !       KINV=LMH
    !       DO K=LMH/2,LMH-1
    !         D2T=T(K-1)-2.0_r8*T(K)+T(K+1)
    !         IF(D2T<0.0_r8)KINV=K
    !       ENDDO
    !
    !       IF(KINV<LMH)THEN
    !         DO K=KINV-1,LMH-1
    !           RDZ=2.0_r8/(Z(K)-Z(K+2))
    !           AKMIN=0.5_r8*RDZ
    !           AKM(K)=MAX(AKM(K),AKMIN)
    !           AKH(K)=MAX(AKH(K),AKMIN)
    !         ENDDO
    !       ENDIF
    !
    !     ENDIF
    !----------------------------------------------------------------------
    !
  END SUBROUTINE DIFCOF
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE VDIFQ                            &
       !   ******************************************************************
       !   *                                                                *
       !   *               VERTICAL DIFFUSION OF Q2 (TKE)                   *
       !   *                                                                *
       !   ******************************************************************
    &(LMH,DTDIF,Q2,EL,Z,kMax)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    INTEGER,INTENT(IN) :: kMax
    !
    INTEGER,INTENT(IN) :: LMH
    !
    REAL(KIND=r8)   ,INTENT(IN) :: DTDIF
    !
    REAL(KIND=r8)   ,INTENT(IN) :: EL  (1:kMax-1)
    REAL(KIND=r8)   ,INTENT(IN) :: Z  (1:kMax+1)
    !
    REAL(KIND=r8)   ,INTENT(INOUT) :: Q2  (1:kMax)  
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: K
    !
    REAL(KIND=r8) :: AKQS,CF,DTOZS   &
         &       ,ESQHF
    !
    REAL(KIND=r8),DIMENSION(1:kMax-2) :: AKQ,CM,CR,DTOZ,RSQ2
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    !***
    !***  VERTICAL TURBULENT DIFFUSION
    !***
    !----------------------------------------------------------------------
    ESQHF=0.5_r8*ESQ
    !
    DO K=1,LMH-2
       DTOZ(K)=(DTDIF+DTDIF)/(Z(K)-Z(K+2))
       AKQ(K)=SQRT((Q2(K)+Q2(K+1))*0.5_r8)*(EL(K)+EL(K+1))*ESQHF         &
            &        /(Z(K+1)-Z(K+2))
       CR(K)=-DTOZ(K)*AKQ(K)
    ENDDO
    !
    CM(1)=DTOZ(1)*AKQ(1)+1.0_r8
    RSQ2(1)=Q2(1)
    !
    DO K=1+1,LMH-2
       CF=-DTOZ(K)*AKQ(K-1)/CM(K-1)
       CM(K)=-CR(K-1)*CF+(AKQ(K-1)+AKQ(K))*DTOZ(K)+1.0_r8
       RSQ2(K)=-RSQ2(K-1)*CF+Q2(K)
    ENDDO
    !
    DTOZS=(DTDIF+DTDIF)/(Z(LMH-1)-Z(LMH+1))
    AKQS=SQRT((Q2(LMH-1)+Q2(LMH))*0.5_r8)*(EL(LMH-1)+ELZ0)*ESQHF        &
         &    /(Z(LMH)-Z(LMH+1))
    !
    CF=-DTOZS*AKQ(LMH-2)/CM(LMH-2)
    !
    Q2(LMH-1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH-2)*CF+Q2(LMH-1))          &
         &        /((AKQ(LMH-2)+AKQS)*DTOZS-CR(LMH-2)*CF+1.0_r8)
    !
    DO K=LMH-2,1,-1
       Q2(K)=(-CR(K)*Q2(K+1)+RSQ2(K))/CM(K)
    ENDDO
    !----------------------------------------------------------------------
    !
  END SUBROUTINE VDIFQ
  !
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------
  SUBROUTINE VDIFH(DTDIF,LMH,THZ0,QZ0,AKHS,CHKLOWQ,CT             &
       &                ,THE,Q,CWM,AKH,Z                                &
       &                ,kMax)
    !     ***************************************************************
    !     *                                                             *
    !     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
    !     *                                                             *
    !     ***************************************************************
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !---------------------------------------------------------------------
    INTEGER,INTENT(IN)    :: kMax
    INTEGER,INTENT(IN)    :: LMH
    REAL(KIND=r8)   ,INTENT(IN)    :: AKHS,CHKLOWQ,CT,DTDIF,QZ0,THZ0
    !
    REAL(KIND=r8)   ,INTENT(IN)    :: AKH   (1:kMax-1)
    REAL(KIND=r8)   ,INTENT(IN)    :: Z     (1:kMax+1)
    REAL(KIND=r8)   ,INTENT(INOUT) :: CWM   (1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: Q     (1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: THE  (1:kMax)
    ! 
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: K
    !
    REAL(KIND=r8)    :: AKHH,AKQS,CF,CMB,CMCB,CMQB,CMTB,CTHF,DTOZL,DTOZS         &
         &           ,RCML,RSCB,RSQB,RSTB
    !
    REAL(KIND=r8),DIMENSION(1:kMax-1) :: AKCT,CM,CR,DTOZ,RSC,RSQ,RST
    !
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    CTHF=0.5_r8*CT
    !
    DO K=1,LMH-1
       DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
       CR(K)=-DTOZ(K)*AKH(K)
       AKCT(K)=AKH(K)*(Z(K)-Z(K+2))*CTHF
    ENDDO
    !
    CM(1)=DTOZ(1)*AKH(1)+1.0_r8
    !----------------------------------------------------------------------
    RST(1)=-AKCT(1)*DTOZ(1)+THE(1)
    RSQ(1)=Q(1)
    RSC(1)=CWM(1)
    !----------------------------------------------------------------------
    DO K=1+1,LMH-1
       DTOZL=DTOZ(K)
       CF=-DTOZL*AKH(K-1)/CM(K-1)
       CM(K)=-CR(K-1)*CF+(AKH(K-1)+AKH(K))*DTOZL+1.0_r8
       RST(K)=-RST(K-1)*CF+(AKCT(K-1)-AKCT(K))*DTOZL+THE(K)
       RSQ(K)=-RSQ(K-1)*CF+Q(K)
       RSC(K)=-RSC(K-1)*CF+CWM(K)
    ENDDO
    !
    DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
    AKHH=AKH(LMH-1)
    !
    CF=-DTOZS*AKHH/CM(LMH-1)
    AKQS=AKHS*CHKLOWQ
    !
    CMB=CR(LMH-1)*CF
    CMTB=-CMB+(AKHH+AKHS)*DTOZS+1.0_r8
    CMQB=-CMB+(AKHH+AKQS)*DTOZS+1.0_r8
    CMCB=-CMB+(AKHH     )*DTOZS+1.0_r8
    !
    RSTB=-RST(LMH-1)*CF+(AKCT(LMH-1)-AKHS*CT)*DTOZS+THE(LMH)
    RSQB=-RSQ(LMH-1)*CF+Q(LMH)
    RSCB=-RSC(LMH-1)*CF+CWM(LMH)
    !----------------------------------------------------------------------
    THE(LMH)=(DTOZS*AKHS*THZ0+RSTB)/CMTB
    Q(LMH)  =(DTOZS*AKQS*QZ0 +RSQB)/CMQB
    CWM(LMH)=(                RSCB)/CMCB
    !----------------------------------------------------------------------
    DO K=LMH-1,1,-1
       RCML=1.0_r8/CM(K)
       THE(K)=(-CR(K)*THE(K+1)+RST(K))*RCML
       Q(K)  =(-CR(K)*  Q(K+1)+RSQ(K))*RCML
       CWM(K)=(-CR(K)*CWM(K+1)+RSC(K))*RCML
    ENDDO
    !----------------------------------------------------------------------
    !
  END SUBROUTINE VDIFH

  !---------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------
  SUBROUTINE VDIFV(LMH,DTDIF,UZ0,VZ0,AKMS,U,V,AKM,Z               &
       ,kMax)
    !     ***************************************************************
    !     *                                                             *
    !     *        VERTICAL DIFFUSION OF VELOCITY COMPONENTS            *
    !     *                                                             *
    !     ***************************************************************
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !---------------------------------------------------------------------
    INTEGER,INTENT(IN) :: kMax
    !
    INTEGER,INTENT(IN) :: LMH
    !
    REAL(KIND=r8)   ,INTENT(IN) :: AKMS,DTDIF,UZ0,VZ0
    !
    REAL(KIND=r8)   ,INTENT(IN) :: AKM (1:kMax-1)
    REAL(KIND=r8)   ,INTENT(IN) :: Z   (1:kMax+1)
    !
    REAL(KIND=r8)   ,INTENT(INOUT) :: U(1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: V(1:kMax)
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: K
    !
    REAL(KIND=r8) :: AKMH,CF,DTOZAK,DTOZL,DTOZS,RCML,RCMVB
    !
    REAL(KIND=r8),DIMENSION(1:kMax-1) :: CM,CR,DTOZ,RSU,RSV
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    DO K=1,LMH-1
       DTOZ(K)= DTDIF/(Z(K)-Z(K+1))
       CR  (K)=-DTOZ(K)*AKM(K)
    ENDDO
    !
    CM(1)=DTOZ(1)*AKM(1)+1.0_r8
    RSU(1)=U(1)
    RSV(1)=V(1)
    !----------------------------------------------------------------------
    DO K=2,LMH-1
       DTOZL=DTOZ(K)
       CF=-DTOZL*AKM(K-1)/CM(K-1)
       CM(K)=-CR(K-1)*CF+(AKM(K-1)+AKM(K))*DTOZL+1.0_r8
       RSU(K)=-RSU(K-1)*CF+U(K)
       RSV(K)=-RSV(K-1)*CF+V(K)
    ENDDO
    !----------------------------------------------------------------------
    DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
    AKMH=AKM(LMH-1)
    !
    CF=-DTOZS*AKMH/CM(LMH-1)
    RCMVB=1.0_r8/((AKMH+AKMS)*DTOZS-CR(LMH-1)*CF+1.0_r8)
    DTOZAK=DTOZS*AKMS
    !----------------------------------------------------------------------
    U(LMH)=(DTOZAK*UZ0-RSU(LMH-1)*CF+U(LMH))*RCMVB
    V(LMH)=(DTOZAK*VZ0-RSV(LMH-1)*CF+V(LMH))*RCMVB
    !----------------------------------------------------------------------
    DO K=LMH-1,1,-1
       RCML=1.0_r8/CM(K)
       U(K)=(-CR(K)*U(K+1)+RSU(K))*RCML
       V(K)=(-CR(K)*V(K+1)+RSV(K))*RCML
    ENDDO
    !----------------------------------------------------------------------
    !
  END SUBROUTINE VDIFV
  !
  !----------------------------------------------------------------------

  !
  !----------------------------------------------------------------------

  SUBROUTINE InitPbl_MellorYamada1()

    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !LOGICAL :: RESTART!INTENT(IN )
    !INTEGER :: ibMax!INTENT(IN )
    !INTEGER :: jbMax!INTENT(IN )
    !INTEGER :: kMax!INTENT(IN )
    !INTEGER :: I,J,K
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !ALLOCATE(RUBLTEN (ibMax,kMax,jbMax))
    !ALLOCATE(RVBLTEN (ibMax,kMax,jbMax))
    !ALLOCATE(RTHBLTEN(ibMax,kMax,jbMax))
    !ALLOCATE(RQVBLTEN(ibMax,kMax,jbMax))
    !ALLOCATE(TKE_MYJ (ibMax,kMax,jbMax))
    !ALLOCATE(EXCH_H  (ibMax,kMax,jbMax))


    !IF(.NOT.RESTART)THEN
    !   DO J=1,jbMax
    !      DO K=1,kMax
    !         DO I=1,ibMax
    !            TKE_MYJ (I,K,J)=EPSQ2
    !            RUBLTEN (I,K,J)=0.0_r8
    !            RVBLTEN (I,K,J)=0.0_r8
    !            RTHBLTEN(I,K,J)=0.0_r8
    !            RQVBLTEN(I,K,J)=0.0_r8
    !            EXCH_H  (I,K,J)=0.0_r8
    !         ENDDO
    !      ENDDO
    !   ENDDO
    !ENDIF

  END SUBROUTINE InitPbl_MellorYamada1
END MODULE Pbl_MellorYamada1

