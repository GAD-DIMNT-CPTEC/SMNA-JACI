!
!  $Author: pkubota $
!  $Date: 2007/03/23 20:23:38 $
!  $Revision: 1.11 $
!
!

MODULE Sfc_MellorYamada1
  USE Constants, ONLY :     &
       r8,gasr,grav

  USE Options, ONLY :  &
       nferr, nfprt,microphys,nClass,nAeros

    IMPLICIT NONE
  SAVE

  PRIVATE
  !  2. Following are constants for use in defining real number bounds.

  !  A really small number.

  REAL(KIND=r8)     , PARAMETER :: epsilon     = 1.E-15_r8
  !  4. Following is information related to the physical constants.

  !  These are the physical constants used within the model.

  ! JM NOTE -- can we name this grav instead?
  REAL(KIND=r8)    , PARAMETER :: g = 9.81_r8  ! acceleration due to gravity (m {s}^-2)

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

  REAL(KIND=r8) , PARAMETER ::  pq0   =379.90516_r8
  !      REAL(KIND=r8) , PARAMETER ::  epsq2=0.2_r8
  REAL(KIND=r8) , PARAMETER ::  epsq2 =  0.02_r8
  REAL(KIND=r8) , PARAMETER ::  a2    = 17.2693882_r8
  REAL(KIND=r8) , PARAMETER ::  a3    =273.16_r8
  REAL(KIND=r8) , PARAMETER ::  a4    = 35.86_r8
  REAL(KIND=r8) , PARAMETER ::  epsq  =  1.e-12_r8
  REAL(KIND=r8) , PARAMETER ::  p608  =rvovrd-1.0_r8 !461.6_r8/287.04_r8 -1!0.608
  REAL(KIND=r8) , PARAMETER ::  climit=   1.e-20_r8
  REAL(KIND=r8) , PARAMETER ::  cm1   =2937.4_r8
  REAL(KIND=r8) , PARAMETER ::  cm2   =   4.9283_r8
  REAL(KIND=r8) , PARAMETER ::  cm3   =  23.5518_r8
  REAL(KIND=r8) , PARAMETER ::  defc  =0.0_r8
  REAL(KIND=r8) , PARAMETER ::  defm  =99999.0_r8
  REAL(KIND=r8) , PARAMETER ::  epsfc =1.0_r8/1.05_r8
  REAL(KIND=r8) , PARAMETER ::  epswet=0.0_r8
  REAL(KIND=r8) , PARAMETER ::  fcdif =1.0_r8/3.0_r8
  REAL(KIND=r8) , PARAMETER ::  fcm   = 0.0_r8
  REAL(KIND=r8) , PARAMETER ::  gma   =-r_d*(1.0_r8-rcp)*0.5_r8
  REAL(KIND=r8) , PARAMETER ::  p400  =40000.0_r8
  REAL(KIND=r8) , PARAMETER ::  phitp =15000.0_r8
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
  !----------------------------------------------------------------------
  !
  ! REFERENCES:  Janjic (2001), NCEP Office Note 437
  !              Mellor and Yamada (1982), Rev. Geophys. Space Phys.
  !              Mellor and Yamada (1974), J. Atmos. Sci.
  !
  ! ABSTRACT:
  !     SfcPbl_MYJ1 GENERATES THE SURFACE EXCHANGE COEFFICIENTS FOR VERTICAL
  !     TURBULENT EXCHANGE BASED UPON MONIN_OBUKHOV THEORY WITH
  !     VARIOUS REFINEMENTS.
  !
  !----------------------------------------------------------------------
  !
  INTEGER :: ITRMX=5 ! Iteration count for sfc layer computations
  !
  REAL(KIND=r8),PARAMETER :: VKARMAN=0.4_r8
  REAL(KIND=r8),PARAMETER :: CAPA=R_D/CP
  REAL(KIND=r8),PARAMETER :: ELOCP=2.72E6_r8/CP !K
  REAL(KIND=r8),PARAMETER :: RCAP=1.0_r8/CAPA
  REAL(KIND=r8),PARAMETER :: GOCP02=G/CP*2.0_r8,GOCP10=G/CP*10.0_r8
  REAL(KIND=r8),PARAMETER :: EPSL=0.10_r8,EPSRU=1.E-7_r8,EPSRS=1.E-7_r8 
  REAL(KIND=r8),PARAMETER :: EPSU2=1.E-4_r8,EPSUST=0.07_r8,EPSZT=1.E-5_r8
  REAL(KIND=r8),PARAMETER :: A1=0.659888514560862645_r8                        &
       &                 ,A2x=0.6574209922667784586_r8                      &
       &                 ,B1=11.87799326209552761_r8                        &
       &                 ,B2=7.226971804046074028_r8                        &
       &                 ,C1=0.000830955950095854396_r8
  REAL(KIND=r8),PARAMETER :: A2S=17.2693882_r8,A3S=273.16_r8,A4S=35.86_r8
  REAL(KIND=r8),PARAMETER :: SEAFC=0.98_r8,PQ0SEA=PQ0*SEAFC
  REAL(KIND=r8),PARAMETER :: BETA=1.0_r8/273.0_r8,CZIL=0.1_r8
  REAL(KIND=r8),PARAMETER :: EXCML=0.001_r8
  REAL(KIND=r8),PARAMETER :: EXCMS=0.001_r8
  REAL(KIND=r8),PARAMETER :: GLKBR=10.0_r8,GLKBS=30.0_r8,PI=3.1415926_r8               &
       &                 ,QVISC=2.1E-5_r8,RIC=0.505_r8,SMALL=0.35_r8              &
       &                 ,SQPR=0.84_r8,SQSC=0.84_r8,SQVISC=258.2_r8,TVISC=2.1E-5_r8  &
       &                 ,USTC=0.007_r8,USTR=0.00225_r8,VISC=1.5E-5_r8                &
       &                 ,WWST=1.2_r8,ZTFC=1.0_r8
  !
  REAL(KIND=r8),PARAMETER :: BTG=BETA*G,CZIV=SMALL*GLKBS                    &
  &                 ,GRRS=GLKBR/GLKBS               &
       &                 ,RB1=1.0_r8/B1,RTVISC=1.0_r8/TVISC,RVISC=1.0_r8/VISC        &
       &                 ,ZQRZT=SQSC/SQPR
  !
  REAL(KIND=r8),PARAMETER :: ADNH= 9.0_r8*A1*A2x*A2x*(12.0_r8*A1+3.0_r8*B2)*BTG*BTG     &                  
       &                    ,ADNM=18.0_r8*A1*A1*A2x*(B2-3.0_r8*A2x)*BTG             & 
       &                    ,ANMH=-9.0_r8*A1*A2x*A2x*BTG*BTG                    &
       &                    ,ANMM=-3.0_r8*A1*A2x*(3.0_r8*A2x+3.0_r8*B2*C1+18.0_r8*A1*C1-B2) &
       &                                 *BTG                             &   
       &                    ,BDNH= 3.0_r8*A2x*(7.0_r8*A1+B2)*BTG                    &
       &                    ,BDNM= 6.0_r8*A1*A1                                 &
       &                    ,BEQH= A2x*B1*BTG+3.0_r8*A2x*(7.0_r8*A1+B2)*BTG         &
       &                    ,BEQM=-A1*B1*(1.0_r8-3.0_r8*C1)+6.0_r8*A1*A1                &
       &                    ,BNMH=-A2x*BTG                                  &     
       &                    ,BNMM=A1*(1.0_r8-3.0_r8*C1)                             &
       &                    ,BSHH=9.0_r8*A1*A2x*A2x*BTG                         &
       &                    ,BSHM=18.0_r8*A1*A1*A2x*C1                          &
       &                    ,BSMH=-3.0_r8*A1*A2x*(3.0_r8*A2x+3.0_r8*B2*C1+12.0_r8*A1*C1-B2) &
       &                                *BTG                             &
       &                 ,CESH=A2x                                       &
       &                 ,CESM=A1*(1.0_r8-3.0_r8*C1)                             &
       &                 ,CNV=EP_1*G/BTG                                 &
       &                 ,ELFCS=VKARMAN*BTG                              &
       &                 ,FZQ1=RTVISC*QVISC*ZQRZT                        &
       &                 ,FZQ2=RTVISC*QVISC*ZQRZT                        &
       &                 ,FZT1=RVISC *TVISC*SQPR                         &
       &                 ,FZT2=CZIV*GRRS*TVISC*SQPR                      &
       &                 ,FZU1=CZIV*VISC                                 &
       &                 ,PIHF=0.5_r8*PI                                    &
       &                 ,RQVISC=1.0_r8/QVISC                                &
       &                 ,RRIC=1.0_r8/RIC                                    &
       &                 ,USTFC=0.018_r8/G                                  &
       &                 ,WWST2=WWST*WWST                                &
       &                 ,ZILFC=-CZIL*VKARMAN*SQVISC
  !
  !
  !----------------------------------------------------------------------
  !***  FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2
  !----------------------------------------------------------------------
  !
  REAL(KIND=r8),PARAMETER :: AEQH=9.0_r8*A1*A2x*A2x*B1*BTG*BTG                  &
       &                      +9.0_r8*A1*A2x*A2x*(12.0_r8*A1+3.0_r8*B2)*BTG*BTG      &
       &                 ,AEQM=3.0_r8*A1*A2x*B1*(3.0_r8*A2x+3.0_r8*B2*C1+18.0_r8*A1*C1-B2) &
       &                      *BTG+18.0_r8*A1*A1*A2x*(B2-3.0_r8*A2x)*BTG
  !
  !----------------------------------------------------------------------
  !***  FORBIDDEN TURBULENCE AREA
  !----------------------------------------------------------------------
  !
  REAL(KIND=r8),PARAMETER :: REQU=-AEQH/AEQM                                &
       &                    ,EPSGH=1.E-9_r8
  REAL(KIND=r8),PARAMETER :: EPSGM=REQU*EPSGH! =1/s**2
  !----------------------------------------------------------------------
  !***  NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT
  !----------------------------------------------------------------------
  !
  REAL(KIND=r8),PARAMETER :: UBRYL=(18.0_r8*REQU*A1*A1*A2x*B2*C1*BTG            &
       &                         +9.0_r8*A1*A2x*A2x*B2*BTG*BTG)              &
       &                        /(REQU*ADNM+ADNH)                        &
       &                 ,UBRY=(1.0_r8+EPSRS)*UBRYL,UBRY3=3.0_r8*UBRY
  !
  REAL(KIND=r8),PARAMETER :: AUBH=27.0_r8*A1*A2x*A2x*B2*BTG*BTG-ADNH*UBRY3      &
       &                 ,AUBM=54.0_r8*A1*A1*A2x*B2*C1*BTG -ADNM*UBRY3       &
       &                 ,BUBH=(9.0_r8*A1*A2x+3.0_r8*A2x*B2)*BTG-BDNH*UBRY3      &
       &                 ,BUBM=18.0_r8*A1*A1*C1           -BDNM*UBRY3        &
       &                 ,CUBR=1.0_r8                     -     UBRY3        &
       &                 ,RCUBR=1.0_r8/CUBR
  !----------------------------------------------------------------------

  INTEGER, PARAMETER :: KZTM=10001,KZTM2=KZTM-2
  REAL(KIND=r8) :: DZETA1,DZETA2,FH01,FH02,ZTMAX1,ZTMAX2,ZTMIN1,ZTMIN2
  REAL(KIND=r8),DIMENSION(KZTM) :: PSIH1,PSIH2,PSIM1,PSIM2

  REAL(KIND=r8),    PARAMETER :: my_eps   =   0.608_r8
  REAL(KIND=r8),    PARAMETER :: my_facl  =   0.1_r8
  INTEGER, PARAMETER :: my_nitr  =   2
  REAL(KIND=r8),    PARAMETER :: my_gkm0  =   1.00_r8
  REAL(KIND=r8),    PARAMETER :: my_gkh0  =   0.10_r8
  REAL(KIND=r8),    PARAMETER :: my_gkm1  =  300.0_r8
  REAL(KIND=r8),    PARAMETER :: my_gkh1  =  300.0_r8
  REAL(KIND=r8),    PARAMETER :: my_vk0   =   0.4_r8
  INTEGER, PARAMETER :: my_kmean =   1
  REAL(KIND=r8),    PARAMETER :: my_a1    =   0.92_r8
  REAL(KIND=r8),    PARAMETER :: my_a2    =   0.74_r8
  REAL(KIND=r8),    PARAMETER :: my_b1    =  16.6_r8
  REAL(KIND=r8),    PARAMETER :: my_b2    =  10.1_r8
  REAL(KIND=r8),    PARAMETER :: my_c1    =   0.08_r8
  REAL(KIND=r8),    PARAMETER :: my_deltx =   0.0_r8
  REAL(KIND=r8), PARAMETER :: my_x        = 1.0_r8
  REAL(KIND=r8), PARAMETER :: my_xx       = my_x*my_x
  REAL(KIND=r8), PARAMETER :: my_g        = 1.0_r8
  REAL(KIND=r8), PARAMETER :: my_gravi    = 9.81_r8
  REAL(KIND=r8), PARAMETER :: my_agrav    = 1.0_r8/my_gravi
  REAL(KIND=r8), PARAMETER :: my_rgas     = 287.0_r8
  REAL(KIND=r8), PARAMETER :: my_akwnmb   = 2.5e-05_r8
  REAL(KIND=r8), PARAMETER :: my_lstar    = 1.0_r8/my_akwnmb
  REAL(KIND=r8), PARAMETER :: my_gocp     = my_gravi/1005.0_r8

  INTEGER         :: my_nbase
  INTEGER         :: my_nthin,my_nthinp

  REAL(KIND=r8)               :: my_alfa
  REAL(KIND=r8)               :: my_beta
  REAL(KIND=r8)               :: my_gama
  REAL(KIND=r8)               :: my_dela
  REAL(KIND=r8)               :: my_r1
  REAL(KIND=r8)               :: my_r2
  REAL(KIND=r8)               :: my_r3
  REAL(KIND=r8)               :: my_r4
  REAL(KIND=r8)               :: my_s1
  REAL(KIND=r8)               :: my_s2
  REAL(KIND=r8)               :: my_rfc
  REAL(KIND=r8), ALLOCATABLE  ::  RMOL_g (:,:)!MONIN-OBUKHOV LENGTH
    REAL(KIND=r8)     :: gbyr
    REAL(KIND=r8)     :: akappa


  PUBLIC :: InitSfc_MellorYamada1
  PUBLIC :: SfcPbl_MYJ1
CONTAINS
  SUBROUTINE InitSfc_MellorYamada1( &
       RESTART, &!(IN   )
       ibMax  , &!(IN   )
       jbMax  , &!(IN   )
       USTAR  , &!(INOUT)
       LOWLYR   )!(INOUT)
    !----------------------------------------------------------------------
    IMPLICIT NONE
    !----------------------------------------------------------------------
    LOGICAL      ,INTENT(IN   ) :: RESTART
    INTEGER      ,INTENT(IN   ) :: ibMax
    INTEGER      ,INTENT(IN   ) :: jbMax
    REAL(KIND=r8),INTENT(INOUT) :: USTAR  (1:ibMax,1:jbMax)
    INTEGER      ,INTENT(INOUT) :: LOWLYR (1:ibMax,1:jbMax)
    !
    REAL(KIND=r8) :: VZ0TBL   (0:30)
    REAL(KIND=r8) :: VZ0TBL_24(0:30)
    !
    INTEGER       :: I
    INTEGER       :: J
    INTEGER       :: K

    REAL(KIND=r8) :: X
    REAL(KIND=r8) :: ZETA1
    REAL(KIND=r8) :: ZETA2
    REAL(KIND=r8) :: ZRNG1
    REAL(KIND=r8) :: ZRNG2
    REAL(KIND=r8) :: PIHF  =3.1415926_r8/2.0_r8
    REAL(KIND=r8) :: EPS   =1.E-6_r8
    REAL(KIND=r8)     :: gam1
    REAL(KIND=r8)     :: gam2
    ALLOCATE(RMOL_g  (ibMax,jbMax))
    gam1=1.0_r8/3.0_r8-2.0_r8*my_a1/my_b1
    
    gam2=(my_b2+6.0_r8*my_a1)/my_b1
    
    my_alfa=my_b1*(gam1-my_c1)+3.0_r8*(my_a2+2.0_r8*my_a1)
    
    my_beta=my_b1*(gam1-my_c1)
    
    my_gama=my_a2/my_a1*(my_b1*(gam1+gam2)-3.0_r8*my_a1)
    
    my_dela=my_a2/my_a1* my_b1* gam1
    my_r1  =0.5_r8*my_gama/my_alfa
    my_r2  =    my_beta/my_gama
    my_r3  =2.0_r8*(2.0_r8*my_alfa*my_dela-my_gama*my_beta)/(my_gama*my_gama)
    my_r4  =my_r2*my_r2
    my_s1  =3.0_r8*my_a2* gam1
    my_s2  =3.0_r8*my_a2*(gam1+gam2)
    !     
    !     critical flux richardson number
    !     
    my_rfc =my_s1/my_s2
    akappa=gasr/cp
    !
    gbyr        =(my_gravi/gasr)**2!(m/sec**2)/(J/(Kg*K))=(m/sec**2)/((Kg*(m/sec**2)*m)/(Kg*K))
    !----------------------------------------------------------------------
    VZ0TBL=                                                          &
         &  (/0.000_r8,                                                          &
         &    2.653_r8,0.826_r8,0.563_r8,1.089_r8,0.854_r8,0.856_r8,0.035_r8,0.238_r8,0.065_r8,0.076_r8  &
         &   ,0.011_r8,0.035_r8,0.011_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8  &
         &   ,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8,0.000_r8/)

    VZ0TBL_24= (/0.000_r8,                                                 &
         &            1.000_r8, 0.070_r8, 0.070_r8, 0.070_r8, 0.070_r8, 0.150_r8,   &
         &            0.080_r8, 0.030_r8, 0.050_r8, 0.860_r8, 0.800_r8, 0.850_r8,   &
         &            2.650_r8, 1.090_r8, 0.800_r8, 0.001_r8, 0.040_r8, 0.050_r8,   &
         &            0.010_r8, 0.040_r8, 0.060_r8, 0.050_r8, 0.030_r8, 0.001_r8,   &
         &            0.000_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.000_r8/)

    !----------------------------------------------------------------------
    !
    !
    !
    !***  FOR NOW, ASSUME SIGMA MODE FOR LOWEST MODEL LAYER
    !
    DO J=1,jbMax
       DO I=1,ibMax
          LOWLYR(I,J)=1
          !       USTAR(I,J)=EPSUST
       ENDDO
    ENDDO
    !----------------------------------------------------------------------
    IF(.NOT.RESTART)THEN
       DO J=1,jbMax
          DO I=1,ibMax
             USTAR(I,J)=0.1_r8
          ENDDO
       ENDDO
    ENDIF
    !----------------------------------------------------------------------
    !
    !***  COMPUTE SURFACE LAYER INTEGRAL FUNCTIONS
    !
    !----------------------------------------------------------------------
    FH01=1.0_r8
    FH02=1.0_r8
    !
          ZTMIN1=-10.0_r8
          ZTMAX1=2.0_r8
          ZTMIN2=-10.0_r8
          ZTMAX2=2.0_r8
    !ZTMIN1=-5.0_r8
    !ZTMAX1=1.0_r8
    !ZTMIN2=-5.0_r8
    !ZTMAX2=1.0_r8
    !
    ZRNG1=ZTMAX1-ZTMIN1
    ZRNG2=ZTMAX2-ZTMIN2
    !
    DZETA1=ZRNG1/(KZTM-1)
    DZETA2=ZRNG2/(KZTM-1)
    !
    !----------------------------------------------------------------------
    !***  FUNCTION DEFINITION LOOP
    !----------------------------------------------------------------------
    !
    ZETA1=ZTMIN1
    ZETA2=ZTMIN2
    !
    DO K=1,KZTM
       !
       !----------------------------------------------------------------------
       !***  UNSTABLE RANGE
       !----------------------------------------------------------------------
       !
       IF(ZETA1<0.0_r8)THEN
          !
          !----------------------------------------------------------------------
          !***  PAULSON 1970 FUNCTIONS
          !----------------------------------------------------------------------
          X=SQRT(SQRT(1.0_r8-16.0_r8*ZETA1))
          !
          PSIM1(K)=-2.0_r8*LOG((X+1.0_r8)/2.0_r8)-LOG((X*X+1.0_r8)/2.0_r8)+2.0_r8*ATAN(X)-PIHF
          PSIH1(K)=-2.0_r8*LOG((X*X+1.0_r8)/2.0_r8)
          !
          !----------------------------------------------------------------------
          !***  STABLE RANGE
          !----------------------------------------------------------------------
          !
       ELSE
          !
          !----------------------------------------------------------------------
          !***  PAULSON 1970 FUNCTIONS
          !----------------------------------------------------------------------
          !
          !         PSIM1(K)=5.0_r8*ZETA1
          !         PSIH1(K)=5.0_r8*ZETA1
          !----------------------------------------------------------------------
          !***   HOLTSLAG AND DE BRUIN 1988
          !----------------------------------------------------------------------
          !
          PSIM1(K)=0.7_r8*ZETA1+0.75_r8*ZETA1*(6.0_r8-0.35_r8*ZETA1)*EXP(-0.35_r8*ZETA1)
          PSIH1(K)=0.7_r8*ZETA1+0.75_r8*ZETA1*(6.0_r8-0.35_r8*ZETA1)*EXP(-0.35_r8*ZETA1)
          !----------------------------------------------------------------------
          !
       ENDIF
       !
       !----------------------------------------------------------------------
       !***  UNSTABLE RANGE
       !----------------------------------------------------------------------
       !
       IF(ZETA2<0.0_r8)THEN
          !
          !----------------------------------------------------------------------
          !***  PAULSON 1970 FUNCTIONS
          !----------------------------------------------------------------------
          !
          X=SQRT(SQRT(1.0_r8-16.0_r8*ZETA2))
          !
          PSIM2(K)=-2.0_r8*LOG((X+1.0_r8)/2.0_r8)-LOG((X*X+1.0_r8)/2.0_r8)+2.0_r8*ATAN(X)-PIHF
          PSIH2(K)=-2.0_r8*LOG((X*X+1.0_r8)/2.0_r8)
          !----------------------------------------------------------------------
          !***  STABLE RANGE
          !----------------------------------------------------------------------
          !
       ELSE
          !
          !----------------------------------------------------------------------
          !***  PAULSON 1970 FUNCTIONS
          !----------------------------------------------------------------------
          !
          !         PSIM2(K)=5.0_r8*ZETA2
          !         PSIH2(K)=5.0_r8*ZETA2
          !
          !----------------------------------------------------------------------
          !***  HOLTSLAG AND DE BRUIN 1988
          !----------------------------------------------------------------------
          !
          PSIM2(K)=0.7_r8*ZETA2+0.75_r8*ZETA2*(6.0_r8-0.35_r8*ZETA2)*EXP(-0.35_r8*ZETA2)
          PSIH2(K)=0.7_r8*ZETA2+0.75_r8*ZETA2*(6.0_r8-0.35_r8*ZETA2)*EXP(-0.35_r8*ZETA2)
          !----------------------------------------------------------------------
          !
       ENDIF
       !
       !----------------------------------------------------------------------
       IF(K==KZTM)THEN
          ZTMAX1=ZETA1
          ZTMAX2=ZETA2
       ENDIF
       !
       ZETA1=ZETA1+DZETA1
       ZETA2=ZETA2+DZETA2
       !----------------------------------------------------------------------
    ENDDO
    !----------------------------------------------------------------------
    ZTMAX1=ZTMAX1-EPS
    ZTMAX2=ZTMAX2-EPS
    !----------------------------------------------------------------------
    !
  END SUBROUTINE InitSfc_MellorYamada1



  !----------------------------------------------------------------------
  SUBROUTINE SfcPbl_MYJ1(prsi,prsl,phii,phil,ITIMESTEP,HT,DZ       & 
       &            ,PMID,PINT,TH,T,QV,QC,U,V,Q2                    &
       &            ,TSK,QSFC,THZ0,QZ0,UZ0,VZ0                      &
       &            ,LOWLYR,XLAND                                    &
       &            ,USTAR,ZNT,PBLH,ELM,MAVAIL            &
       &            ,AKHS,AKMS                                      &
       &            ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC            &
       &            ,QGH,CPM,CT                                     &
       &            ,U10,V10,TSHLTR,TH10,QSHLTR,Q10,PSHLTR            &
       &            ,tmsfc,qmsfc,umsfc,gt,gq,gu,gv,bps,PBL_CoefKm   &
       &            ,PBL_CoefKh,latitu,dt &
       &            ,nCols,kMax,nClass,nAeros)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    INTEGER,INTENT(IN) :: nCols
    INTEGER,INTENT(IN) :: kMax
    INTEGER,INTENT(IN) :: nClass
    INTEGER,INTENT(IN) :: nAeros
    !
    INTEGER,INTENT(IN) :: ITIMESTEP!-- itimestep     number of time steps
    REAL(KIND=r8),    INTENT(IN   ) :: prsi  (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8),    INTENT(IN   ) :: prsl  (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8),    INTENT(IN   ) :: phii  (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8),    INTENT(IN   ) :: phil  (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)

    !
    INTEGER,INTENT(IN) :: LOWLYR(1:nCols)!index of lowest model layer above ground
    !
    REAL(KIND=r8)   ,INTENT(IN) :: HT    (1:nCols)!"HGT" "Terrain Height"   "m"
    REAL(KIND=r8)   ,INTENT(IN) :: XLAND (1:nCols)!land mask (1 for land, 2 for water)
    REAL(KIND=r8)   ,INTENT(IN) :: TSK   (1:nCols)!surface temperature (K)
    REAL(KIND=r8)   ,INTENT(IN) :: MAVAIL(1:nCols)!surface moisture availability (between 0 and 1)
    !
    REAL(KIND=r8)   ,INTENT(IN) :: DZ    (1:nCols,1:kMAx)! dz between full levels (m)
    REAL(KIND=r8)   ,INTENT(IN) :: PMID  (1:nCols,1:kMAx)! p_phy         pressure (Pa)
    REAL(KIND=r8)   ,INTENT(IN) :: PINT  (1:nCols,1:kMAx)! p8w           pressure at full levels (Pa)
    REAL(KIND=r8)   ,INTENT(IN) :: Q2    (1:nCols,1:kMAx)! tke_myj            turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
    REAL(KIND=r8)   ,INTENT(IN) :: QC    (1:nCols,1:kMAx)! "Cloud water mixing ratio"      "kg kg-1"
    REAL(KIND=r8)   ,INTENT(IN) :: QV    (1:nCols,1:kMAx)! "Water vapor mixing ratio"      "kg kg-1"
    REAL(KIND=r8)   ,INTENT(IN) :: T     (1:nCols,1:kMAx)!t_phy            temperature (K)
    REAL(KIND=r8)   ,INTENT(IN) :: TH    (1:nCols,1:kMAx)!th_phy potential temperature (K)
    REAL(KIND=r8)   ,INTENT(IN) :: U     (1:nCols,1:kMAx)!u_phy            u-velocity interpolated to theta points (m/s)
    REAL(KIND=r8)   ,INTENT(IN) :: V     (1:nCols,1:kMAx)!v_phy            v-velocity interpolated to theta points (m/s)
    !
    REAL(KIND=r8)   ,INTENT(OUT) :: FLX_LH(1:nCols)!-- LH            net upward latent heat flux at surface (W/m^2)
    REAL(KIND=r8)   ,INTENT(OUT) :: HFX   (1:nCols)!-- HFX           net upward heat flux at the surface (W/m^2)
    REAL(KIND=r8)   ,INTENT(OUT) :: PSHLTR(1:nCols)!-- pshltr        diagnostic shelter (2m) pressure from MYJ (Pa)
    REAL(KIND=r8)   ,INTENT(OUT) :: QFX   (1:nCols)!-- QFX           net upward moisture flux at the surface (kg/m^2/s)
    REAL(KIND=r8)   ,INTENT(OUT) :: Q10   (1:nCols)!-- q10           diagnostic 10-m specific humidity from MYJ
    REAL(KIND=r8)   ,INTENT(OUT) :: QSHLTR(1:nCols)!-- qshltr        diagnostic 2-m specific humidity from MYJ 
    REAL(KIND=r8)   ,INTENT(OUT) :: TH10  (1:nCols)!-- th10          diagnostic 10-m theta from MYJ
    REAL(KIND=r8)   ,INTENT(OUT) :: TSHLTR(1:nCols)!-- tshltr        diagnostic 2-m theta from MYJ
    REAL(KIND=r8)   ,INTENT(OUT) :: U10   (1:nCols)! u10 diagnostic 10-m u component from surface layer
    REAL(KIND=r8)   ,INTENT(OUT) :: V10   (1:nCols)! v10 diagnostic 10-m v component from surface layer
    !
    REAL(KIND=r8)  ,INTENT(INOUT) :: AKHS (1:nCols)! akhs sfc exchange coefficient of heat/moisture from MYJ
    REAL(KIND=r8)  ,INTENT(INOUT) :: AKMS (1:nCols)! akms sfc exchange coefficient of momentum from MYJ
    REAL(KIND=r8)  ,INTENT(INOUT) :: PBLH (1:nCols)! PBLH PBL height (m)
    REAL(KIND=r8)  ,INTENT(INOUT) :: ELM  (1:nCols,1:kMAx)!
    REAL(KIND=r8)  ,INTENT(INOUT) :: QSFC (1:nCols)! qsfc specific humidity at lower boundary (kg/kg)
    !
    REAL(KIND=r8)  ,INTENT(INOUT) :: QZ0  (1:nCols)
    REAL(KIND=r8)  ,INTENT(INOUT) :: THZ0 (1:nCols)! thz0 potential temperature at roughness length (K)
    REAL(KIND=r8)  ,INTENT(INOUT) :: USTAR(1:nCols)! UST  u* in similarity theory (m/s)
    REAL(KIND=r8)  ,INTENT(INOUT) :: UZ0  (1:nCols)! uz0  u wind component at roughness length (m/s)
    REAL(KIND=r8)  ,INTENT(INOUT) :: VZ0  (1:nCols)! vz0  v wind component at roughness length (m/s)
    REAL(KIND=r8)  ,INTENT(INOUT) :: ZNT  (1:nCols)! ZNT  time-varying roughness length (m)
    !
    REAL(KIND=r8)  ,INTENT(OUT )  :: CHS  (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: CHS2 (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: CQS2 (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: CPM  (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: CT   (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: FLHC (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: FLQC (1:nCols)
    REAL(KIND=r8)  ,INTENT(OUT )  :: QGH  (1:nCols)
    REAL(KIND=r8)  ,INTENT(inout) :: tmsfc(ncols,kmax,3)
    REAL(KIND=r8)  ,INTENT(inout) :: qmsfc(ncols,kmax,5+nClass+nAeros)
    REAL(KIND=r8)  ,INTENT(inout) :: umsfc(ncols,kmax,4)
    REAL(KIND=r8)  ,INTENT(IN)    :: dt
    REAL(KIND=r8)  ,INTENT(IN)    :: gt  (1:nCols,1:kMAx)
    REAL(KIND=r8)  ,INTENT(IN)    :: gq  (1:nCols,1:kMAx)
    REAL(KIND=r8)  ,INTENT(IN)    :: gu  (1:nCols,1:kMAx)
    REAL(KIND=r8)  ,INTENT(IN)    :: gv  (1:nCols,1:kMAx)   
    REAL(KIND=r8)  ,INTENT(IN)    :: bps (1:nCols,1:kMAx)   
    REAL(KIND=r8), INTENT(INOUT) :: PBL_CoefKm(ncols, kmax)
    REAL(KIND=r8), INTENT(INOUT) :: PBL_CoefKh(ncols, kmax)
    REAL(KIND=r8)    :: aa0(ncols,kmax)
    REAL(KIND=r8)    :: bb0(ncols,kmax)
    REAL(KIND=r8)    :: tt0(ncols,kmax)
    REAL(KIND=r8)    :: tt1(ncols,kmax)
    REAL(KIND=r8)    :: sigr  (ncols,kmax)
    REAL(KIND=r8)    :: sigriv(ncols,kmax)
    REAL(KIND=r8)    :: con0  (ncols,kmax)
    REAL(KIND=r8)    :: con1  (ncols,kmax)
    REAL(KIND=r8)    :: con2  (ncols,kmax)

    INTEGER, INTENT(in        ) :: latitu

    REAL(KIND=r8)  :: Pbl_ATemp(1:nCols,1:kMAx)
    REAL(KIND=r8)  :: Pbl_ITemp(1:nCols,1:kMAx)
    REAL(KIND=r8)  :: WSTAR(1:nCols)
    REAL(KIND=r8)  :: RIB(1:nCols,kmax)
    REAL(KIND=r8)  :: twodti,twodt
    REAL(KIND=r8)  :: my_a(ncols,kMax)
    REAL(KIND=r8)  :: my_b(ncols,kMax)
    REAL(KIND=r8)  :: RMOL (1:nCols)!MONIN-OBUKHOV LENGTH
    
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: I,K,KFLIP,LMH,LPBL,NTSD
    !
    REAL(KIND=r8) :: A,ADEN,APESFC,AUBR,B,BDEN,BUBR,CWMLOW               &
         &       ,ELOQ2X,FIS,GHK,GMK                                       &
         &       ,P02P,P10P,PLOW,PSFC,PTOP,QLOW,QOL2ST,QOL2UN,QS02,QS10    &
         &       ,RAPA,RAPA02,RAPA10,RATIOMX,RDZ,SEAMASK,SM                &
         &       ,T02P,T10P,TEM,TH02P,TH10P,THLOW,THELOW,THM               &
         &       ,TLOW,TZ0,ULOW,VLOW,FAC
    !
    REAL(KIND=r8),DIMENSION(1:kMax) :: CWMK,PK,Q2K,QK,THEK,THK,TK,UK,VK
    !
    REAL(KIND=r8),DIMENSION(1:kMax-1) :: GH,GM
    !
    REAL(KIND=r8),DIMENSION(1:kMax+1) :: ZHK
    !
    REAL(KIND=r8),DIMENSION(1:nCols) :: THSK,ZSL
    !
    REAL(KIND=r8),DIMENSION(1:nCols,1:kMax+1) :: ZINT
    !
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    !
    !***  MAKE PREPARATIONS
    !
    !----------------------------------------------------------------------
    !        setup_integration:  DO J=jMax0,jMax

    DO K=1,kMax+1
       DO I=1,nCols
          ZINT(I,K)=0.0_r8
       ENDDO
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ENDDO
    !
    DO I=1,nCols
       RMOL (i)=RMOL_g (i,latitu)!MONIN-OBUKHOV LENGTH
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      DO J=jMax0,jMax
    DO I=1,nCols
       ZINT(I,kMax+1)=HT(I)     ! Z at bottom of lowest sigma layer
       PBLH(I)=-1.0_r8
       !
!!!!!!
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!
!!!!!!  ZINT(I,kMax+1)=1.E-4_r8         ! Z of bottom of lowest eta layer
!!!!!!  ZHK(kMax+1)=1.E-4_r8            ! Z of bottom of lowest eta layer
       !
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ENDDO
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     DO J=jMax0,jMax
    DO K=kMax,1,-1
       KFLIP=kMax+1-K
       DO I=1,nCols
          ZINT(I,K)=ZINT(I,K+1)+DZ(I,KFLIP)
       ENDDO
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ENDDO
    !
    NTSD=ITIMESTEP
    !
    IF(NTSD==1)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        DO J=jMax0,jMax
       DO I=1,nCols
          USTAR(I)=0.1_r8
          FIS=HT(I)*G
          SM=XLAND(I)-1.0_r8
!!!       Z0 (I)=SM*Z0SEA+(1.0_r8-SM)*(Z0 (I)*Z0MAX+FIS*FCM+Z0LAND)
!!!       ZNT(I)=SM*Z0SEA+(1.0_r8-SM)*(ZNT(I)*Z0MAX+FIS*FCM+Z0LAND)
       ENDDO
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ENDDO

    !
!!!!  IF(NTSD==1)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        DO J=jMax0,jMax
    DO I=1,nCols
       CT(I)=0.0_r8
    ENDDO
    !    ENDDO
!!!!  ENDIF
    !
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !
    DO I=1,nCols
       !
       !***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
       !
       LMH=kMax-LOWLYR(I)+1
       !
       PTOP=PINT(I,kMax)      ! kMax+1=kMAx
       PSFC=PINT(I,LOWLYR(I))
       ! Define THSK here (for first timestep mostly)
       THSK(I)=TSK(I)*bps(i,1)!/(PSFC*1.E-5_r8)**CAPA
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
          THK (K)=TH(I,KFLIP) ! K
          TK  (K)=T(I,KFLIP)! K
          RATIOMX=QV(I,KFLIP)! kg/kg
          QK  (K)=RATIOMX/(1.0_r8 + RATIOMX)! kg/kg
          PK  (K)=PMID(I,KFLIP)!Pa
          CWMK(K)=QC(I,KFLIP)! kg/kg
          THEK(K)=(CWMK(K)*(-ELOCP/TK(K))+1.0_r8)*THK(K)! K
          Q2K (K)=2.0_r8*Q2(I,KFLIP)
          !
          !
          !***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
          !
          ZHK(K)=ZINT(I,K)
          !
       ENDDO
       ZHK(kMax+1)=HT(I) !m         ! Z at bottom of lowest sigma layer
       !
       DO K=kMax,1,-1
          KFLIP=kMax+1-K
          UK(K)=U(I,KFLIP)!m/s
          VK(K)=V(I,KFLIP)!m/s
       ENDDO
       !
       !----------------------------------------------------------------------
       !***  COMPUTE THE HEIGHT OF THE BOUNDARY LAYER
       !----------------------------------------------------------------------
       !
       DO K=1,LMH-1
          !
          !          2.0_r8
          ! RDZ = ------------------
          !        ZHK(K)-ZHK(K+2)
          !
          RDZ=2.0_r8/(ZHK(K)-ZHK(K+2))
          !
          !
          !    4*[ ((UK(K)-UK(K+1))**2  + (VK(K) - VK(K+1))**2)]     1
          !GMK=--------------------------------------------------=-------
          !            (ZHK(K) - ZHK(K+2))**2                        s*s
          !
          !
          GMK=((UK(K)-UK(K+1))**2+(VK(K)-VK(K+1))**2)*RDZ*RDZ
          !
          !               1
          ! GM(K) = -------
          !              s*s
          !
          GM(K)=MAX(GMK,EPSGM)
          !
          !             TK(K) + TK(K+1)
          ! TEM = -------------------- = K
          !                  2
          !
          TEM=(TK  (K) + TK  (K+1))*0.5_r8
          !
          !             THEK(K) + THEK(K+1)
          ! THM = ----------------------- = K
          !                  2
          !
          THM=(THEK(K)+THEK(K+1))*0.5_r8
          !
          !
          !
          !
          A=THM*P608! =T*(Rv/Rd - 1)
          !
          !
          !
          !
          B=(ELOCP/TEM - 1.0_r8-P608)*THM
          !
          GHK=((THEK(K)-THEK(K+1)+CT(I))                                                &
               &          *((QK(K)+QK(K+1)+CWMK(K)+CWMK(K+1))*(0.5_r8*P608)+1.0_r8)     &
               &          +( QK(K)-QK(K+1)+CWMK(K)-CWMK(K+1))*A                          &
               &          +(CWMK(K)-CWMK(K+1))*B)*RDZ
          !
          !
          !     
          IF(ABS(GHK)<=EPSGH)GHK=EPSGH
          !
          GH(K)=GHK ! K/m
       ENDDO
       !
       !***  FIND MAXIMUM MIXING LENGTHS AND THE LEVEL OF THE PBL TOP
       !
       LPBL=LMH
       !
       DO K=1,LMH-1
          GMK=GM(K)!1/s*s
          GHK=GH(K)!K/m
          !
          IF(GHK>=EPSGH)THEN
             !
             IF(GMK/GHK<=REQU)THEN
                ELM(I,K)=EPSL
                LPBL=K
             ELSE
                AUBR=(AUBM*GMK+AUBH*GHK)*GHK
                BUBR= BUBM*GMK+BUBH*GHK
                QOL2ST=(-0.5_r8*BUBR+SQRT(BUBR*BUBR*0.25_r8-AUBR*CUBR))*RCUBR
                ELOQ2X=1.0_r8/QOL2ST
                ELM(I,K)=MAX(SQRT(ELOQ2X*Q2K(K)),EPSL)
             ENDIF
             !
          ELSE
             ADEN=(ADNM*GMK+ADNH*GHK)*GHK
             BDEN= BDNM*GMK+BDNH*GHK
             QOL2UN=-0.5_r8*BDEN+SQRT(BDEN*BDEN*0.25_r8-ADEN)
             ELOQ2X=1.0_r8/(QOL2UN+EPSRU)  !  repsr1/qol2un
             ELM(I,K)=MAX(SQRT(ELOQ2X*Q2K(K)),EPSL)
          ENDIF
          !
       ENDDO
       !
       !
       IF(ELM(I,LMH-1)==EPSL)LPBL=LMH
       !
       !***  THE HEIGHT OF THE PBL
       !
       PBLH(I)=ZHK(LPBL)-ZHK(LMH+1)
       !
       !----------------------------------------------------------------------
       !***
       !***  FIND THE SURFACE EXCHANGE COEFFICIENTS
       !***
       !----------------------------------------------------------------------
       PLOW=PK(LMH)
       TLOW=TK(LMH)
       THLOW=THK(LMH)
       THELOW=THEK(LMH)
       QLOW=QK(LMH)
       CWMLOW=CWMK(LMH)
       ULOW=UK(LMH)
       VLOW=VK(LMH)
       ZSL(I)=(ZHK(LMH)-ZHK(LMH+1))*0.5_r8
       APESFC=(PSFC*1.E-5_r8)**CAPA
       TZ0=THZ0(I)/bps(i,1)!*APESFC
       !
       CALL SFCDIF(NTSD,SEAMASK,THSK(I),QSFC(I),PSFC            &
            &      ,UZ0(I),VZ0(I),TZ0,THZ0(I),QZ0(I)             &
            &      ,USTAR(I),ZNT(I),CT(I),RMOL(I)             &
            &      ,AKMS(I),AKHS(I),PBLH(I),MAVAIL(I)             &
            &      ,CHS(I),CHS2(I),CQS2(I)                       &
            &      ,HFX(I),QFX(I),FLX_LH(I)                       &
            &      ,FLHC(I),FLQC(I),QGH(I),CPM(I)             &
            &      ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW             &
            &      ,ZSL(I),PLOW                                             &
            &      ,U10(I),V10(I),TSHLTR(I),TH10(I)             &
            &      ,QSHLTR(I),Q10(I),PSHLTR(I),WSTAR(I),RIB(I,k))
       !
       !***  REMOVE SUPERATURATION AT 2M AND 10M
       !
       RAPA=APESFC
       TH02P=TSHLTR(I)
       TH10P=TH10(I)
       !
       RAPA02=RAPA-GOCP02/TH02P
       RAPA10=RAPA-GOCP10/TH10P
       !
       T02P=TH02P*RAPA02
       T10P=TH10P*RAPA10
       !
       P02P=(RAPA02**RCAP)*1.E5_r8
       P10P=(RAPA10**RCAP)*1.E5_r8
       !
       QS02=PQ0/P02P*EXP(A2*(T02P-A3)/(T02P-A4))
       QS10=PQ0/P10P*EXP(A2*(T10P-A3)/(T10P-A4))
       !
       IF(QSHLTR(I)>QS02)QSHLTR(I)=QS02
       IF(Q10   (I)>QS10)Q10   (I)=QS10
       !----------------------------------------------------------------------
       !
    ENDDO
    !          Pbl_ATemp (7)   temperature at the interface of two adjacent layers
    !
    !                  k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(2)
    !                  k=3/2----si,ric,rf,km,kh,b,l -----------
    !                  k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(1)
    !                  k=1/2----si ----------------------------
    !
    !                           sl(k+1)-si(k+1)
    !            t0    (   k)=-------------------
    !                           sl(k+1)-sl(k  )
    !
    !                           si(k+1)-sl(k  )
    !            t1    (   k)=-------------------
    !                           sl(k+1)-sl(k  )
    !
    !
    !                           sl(k+1)-si(k+1)                si(k+1)-sl(k  )
    !        Pbl_ATemp(i,k)  =-------------------*gt(i,k) + ------------------*gt(i,k+1)
    !                           sl(k+1)-sl(k  )                sl(k+1)-sl(k  )
    !
    !
    !
    !
    DO k = 1, kmax-1
       DO i = 1, ncols
          Pbl_CoefKm(i,k) = -AKMS(I)*(PBLH(I)*WSTAR(I))
          Pbl_CoefKh(i,k) = -AKHS(I)*(PBLH(I)*WSTAR(I))
       END DO
    END DO
    fac=0.25_r8
    IF (kmax >= 4) THEN
       DO k = 2, kmax-2
              DO i = 1, ncols
                 !
                 !             k=2  ****Km(k),sl*** } -----------
                 !             k=3/2----si,ric,rf,km,kh,b,l -----------
                 !             k=1  ****Km(k),sl*** } -----------
                 !             k=1/2----si ----------------------------
                 !
                 !       Km(k-1) + 2*Km(k) + Km(k+1)
                 ! Km = -------------------------------
                 !                  4
                 !
                 Pbl_CoefKm(i,k)=fac*(Pbl_CoefKm(i,k-1)+2.0_r8*Pbl_CoefKm(i,k)+Pbl_CoefKm(i,k+1))
                 !
                 !       Kh(k-1) + 2*Kh(k) + Kh(k+1)
                 ! Kh = -------------------------------
                 !                  4
                 !
                 Pbl_CoefKh(i,k)=fac*(Pbl_CoefKh(i,k-1)+2.0_r8*Pbl_CoefKh(i,k)+Pbl_CoefKh(i,k+1))
              END DO
       END DO
    END IF
    DO i = 1, ncols
       !
       !       Km(1) + Km(2)
       ! Km = ---------------
       !            2
       !
       !
       !       Kh(1) + Kh(2)
       ! Kh = ---------------
       !            2
       !
       Pbl_CoefKm(i,         1)=0.5_r8*(Pbl_CoefKm(i,     1) + Pbl_CoefKm(i,     2))
       Pbl_CoefKh(i,         1)=0.5_r8*(Pbl_CoefKh(i,     1) + Pbl_CoefKh(i,     2))
       Pbl_CoefKm(i,kmax-1)=0.5_r8*(Pbl_CoefKm(i,kmax-1) + Pbl_CoefKm(i,kmax-2))
       Pbl_CoefKh(i,kmax-1)=0.5_r8*(Pbl_CoefKh(i,kmax-1) + Pbl_CoefKh(i,kmax-2))
    END DO

    DO k = 1, kmax-1
       DO i = 1, ncols
          Pbl_CoefKm(i,k) = MIN(my_gkm1,MAX(my_gkm0,Pbl_CoefKm(i,k)))
          Pbl_CoefKh(i,k) = MIN(my_gkh1,MAX(my_gkh0,Pbl_CoefKh(i,k)))
       END DO
    END DO
    DO i = 1, ncols
       Pbl_CoefKm(i,kmax) = MIN(my_gkm1,MAX(my_gkm0,Pbl_CoefKm(i,kmax-1)))
       Pbl_CoefKh(i,kmax) = MIN(my_gkh1,MAX(my_gkh0,Pbl_CoefKh(i,kmax-1)))
    END DO

    twodt=(2.0_r8*dt)
    twodti=1.0_r8/twodt

    DO i=1,nCols
       aa0    (i,kmax)=0.0_r8
       tt0   (i,kmax)=0.0_r8
       tt1   (i,kmax)=0.0_r8
       bb0    (i,   1)=0.0_r8
       sigr  (i,kmax)=0.0_r8
       sigriv(i,   1)=0.0_r8

    END DO
    DO k = 1, kmax-1
       DO i=1,nCols
       ! con0  (k)=gasr*delsig(k)/(grav*sig(k))

        con0  (i,k)=gasr*((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))/(grav*(prsl(i,k)/prsi(i,1)))

        !con1  (   k)=grav*sigml(k+1)/(gasr*(sig(k)-sig(k+1)))

        con1  (i,k) = grav*(prsi(i,k+1)/prsi(i,1))/(gasr*(((prsl(i,k)/prsi(i,1)) - (prsl(i,k+1)/prsi(i,1)))))

        !con2  (   k)=grav*con1(k)
        con2  (i,k)=grav*con1(i,k)

        !con1  (   k)=con1(k)*con1(k)
        con1  (i,k)=con1(i,k)*con1(i,k)

!       t0    (   k)=(sig(k+1)-sigml(k+1))/(sig(k+1)-sig(k))

        tt0    (i,   k)=((prsl(i,k+1)/prsi(i,1)) - (prsi  (i,k+1)/prsi  (i,1))) / &
                        ((prsl(i,k+1)/prsi(i,1)) - (prsl(i,k)/prsi(i,1)))


 !       t1    (   k)=(sigml(k+1)-sig(k  ))/(sig(k+1)-sig(k))
        tt1    (i,   k)=((prsi  (i,k+1)/prsi  (i,1)) - (prsl(i,k)/prsi(i,1)))/&
                        ((prsl  (i,k+1)/prsi  (i,1)) - (prsl(i,k)/prsi(i,1)))

       !sigr  (   k)=sigk(k)*sigkiv(k+1)
        sigr  (i,k)= ((prsl(i,k)/prsi(i,1))**akappa)*  (1.0_r8/((prsl(i,k+1)/prsi(i,1))**akappa))
       !sigriv( k+1)=sigk(k+1)*sigkiv(k)
        sigriv(i,k+1)= ((prsl(i,k+1)/prsi(i,1))**akappa )  *  (1.0_r8/((prsl(i,k)/prsi(i,1))**akappa))

       !a0    (k   )=gbyr*sigml(k+1)**2/(delsig(k  )*(sig(k)-sig(k+1)))
       !REAL(KIND=r8),    INTENT(IN   ) :: prsi  (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
       !REAL(KIND=r8),    INTENT(IN   ) :: prsl  (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa

       aa0    (i,k   )=gbyr*((prsi  (i,k+1)/prsi  (i,1))**2)/(((prsi(i,k)/prsi(i,1)) - (prsi(i,k+1)/prsi(i,1)))*((prsl(i,k)/prsi(i,1)) - (prsl(i,k+1)/prsi(i,1))))

       !bb0    (i,k+1 )=gbyr*sigml(k+1)**2/(delsig(k+1)*(sig(k)-sig(k+1)))

       bb0    (i,k+1 )=gbyr*((prsi  (i,k+1)/prsi  (i,1))**2)/(((prsi(i,k+1)/prsi(i,1)) - (prsi(i,k+2)/prsi(i,1)))*((prsl(i,k)/prsi(i,1)) - (prsl(i,k+1)/prsi(i,1))))

       END DO
    END DO

    DO k = 1, kmax
       DO i=1,nCols

       !                 --                     --
       !    1           1        | g          si(k+1)      |     con1(k)
       !  ------ =--- * |--- * ----------------| = ----------
       !    DZ     T        | R        sl(k) -sl(k+1) |       T
       !                 --                     --
       !                -- --
       !    DA      d  |     |
       !  ------ =---- | W'A'|
       !    DT      dZ |     |
       !                -- --
       !
       !                          ----
       !    DA      d              d  |    |
       !  ------ =---- * K * --- | A  |
       !    DT      dZ        dZ |    |
       !                          ----
       !                     
       !    DA         d      dA 
       !  ------ =K * ---- * ---- 
       !    DT         dZ     dZ 
       !                     
       !
       !
       !a0(k)  =gbyr*sigml(k+1)**2/(delsig(k  )*(sig(k)-sig(k+1)))
       !
       !                grav **2
       !              --------    * si(k+1)**2
       !                gasr**2 
       !a0(k)  = ---------------------------------
       !          ((si(k)-si(k+1))) * (sl(k)-sl(k+1)))
       !
       !         --  -- 2                              -- -- 2
       !        |   g  |       si(k+1)**2             |  1  | 
       !a0(k)  =| -----| * --------------------    = | --- |
       !        |   R  |     ((si(k)-si(k+1))  )     |  dZ | 
       !         --  --                               -- --
       !                   --        --        -- -- 2            --                        -- 2
       !                  |          |    |  1  |            |         m         kg * K   |  
       my_a(i,k)=twodt*aa0(i,k)!  |  2*Dt |  * | --- |    ==> s * |   ------- * --------  | 
       !                  |          |    |  dZ |            |        s**2           J          |  
       !                   --        --        -- --                   --                        --
       !J = F*DX = kg m/s**2 *m = kg * m**2/s**2
       !                   --        --        -- -- 2            --                            -- 2
       !                  |          |    |  1  |            |    m         kg * K *s**2 |    K**2 * s
       !a(k)=twodt*a0(k)! |  2*Dt |  * | --- |    ==> s * | ------- * --------------- | = -----------
       !                  |          |    |  dZ |            |   s**2     kg * m**2      |    m**2  
       !                   --        --        -- --                   --                            --
       !                   --        --        -- -- 2  
       !                  |          |    |  1  |         K**2 * s
       !a(k)=twodt*a0(k)! |  2*Dt |  * | --- |    ==> -----------
       !                  |          |    |  dZ |           m**2  
       !                   --        --        -- --         

       !
       !                   gbyr*sigml(k+1)**2
       !    b0(k) = -----------------------------------------------
       !               (((si(k+1)-si(k+1+1))  ) *(sig(k)-sig(k+1)))

       my_b(i,k)=twodt*bb0(i,k)! s * K**2/m**2 
       END DO

    END DO
    DO k = 1, kmax-1
       DO i = 1, ncols
          Pbl_ATemp (i,k) = tt0(i,k)*T(i,k) + tt1(i,k)*T(i,k+1)
          Pbl_ITemp (i,k) = 1.0_r8/(Pbl_ATemp(i,k) * Pbl_ATemp(i,k))
       END DO
    END DO

    !
    !----------------------------------------------------------------------
    !      ENDDO setup_integration
    !----------------------------------------------------------------------
    !                                 --              --
    !           d -(U'W')           d        |         d U        | 
    !        ------------- = --------|KM *  -------- |
    !           dt                   dt        |         d z        |
    !                                 --              --
    !          momentum diffusion
    ! 
    CALL VDIFV(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKm,my_a,my_b,gu,gv,umsfc)
    !
    !                                 --              --
    !           d -(Q'W')           d        |         d Q        | 
    !        ------------- = --------|KH *  -------- |
    !           dt                   dt        |         d z        |
    !                                 --              --
    !          
    !          water vapour diffusion
    !
    !CALL VDIFH(nClass,nAeros,kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,my_a,my_b,gq,qmsfc)
 
    CALL VDIFH(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,my_a,my_b,gq,qmsfc)       
    !IF (microphys)THEN
    !   IF (microphys.and. (nClass+nAeros)>0 .and. PRESENT(gvar))THEN
    !      CALL VDIFHM(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,my_a,my_b,gice,gliq,qmsfc,gvar)
    !   ELSE 
    !      CALL VDIFHM(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,my_a,my_b,gice,gliq,qmsfc)
    !   END IF
    !END IF

    !          
    !                                 --              --
    !           d -(T'W')           d        |         d T        | 
    !        ------------- = --------|KH *  -------- |
    !           dt                   dt        |         d z        |
    !                                 --              --
    !         
    !          sensible heat diffusion
    ! 
    !CALL VDIFT(kMax,nCols,twodti,sigr,sigriv,Pbl_ITemp,Pbl_CoefKh,my_a,my_b,gt,tmsfc)
    CALL VDIFT(kMax,nCols,twodti,sigr,sigriv,Pbl_ITemp,Pbl_CoefKh,my_a,my_b,gt,tmsfc)

    ! 
    DO I=1,nCols
      RMOL_g (i,latitu)= RMOL (i)!MONIN-OBUKHOV LENGTH
    END DO

  END SUBROUTINE SfcPbl_MYJ1
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  !---------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------
  SUBROUTINE VDIFV(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKm,a,b,gu,gv,gmu)
    !     ***************************************************************
    !     *                                                             *
    !     *        VERTICAL DIFFUSION OF VELOCITY COMPONENTS            *
    !     *                                                             *
    !     ***************************************************************
    !---------------------------------------------------------------------
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)      !s * K**2/m**2 
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)      !s * K**2/m**2 
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)!1/K**2
    REAL(KIND=r8),    INTENT(inout) :: Pbl_CoefKm(nCols,kMax)!(m/sec**2)
    REAL(KIND=r8),    INTENT(in   ) :: gu        (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv        (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmu       (nCols,kMax,4)

    REAL(KIND=r8) :: Pbl_DifVzn(nCols,kMax) 
    REAL(KIND=r8) :: Pbl_DifVmd(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KHbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_KMbyDZ(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_TendU(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_TendV(nCols,kMax) !1/(sec * K)

    INTEGER       :: i
    INTEGER       :: k

    !     
    !     momentum diffusion
    !     
    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKm(i,k  )=Pbl_CoefKm(i,k)*Pbl_ITemp(i,k)!(m/sec**2) * (1/K**2) = m/(sec**2 * K**2)
          !
          !                   --   --      -- -- 2   
          !                  |       |    |  1  |         K**2 * s
          !a(k)=             |  2*Dt |  * | --- |     ==>-----------
          !                  |       |    |  dZ |           m**2  
          !                   --   --      -- --     
          !                                         K**2 * s         m              1
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKm(i,k)!----------- * -------------  = --------
          !                                           m**2          s**2*K**2     m * s
          !
          !                                         K**2 * s        m**2             1
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKm(i,k)!----------- * -------------  = --------
          !                                           m**2          s**2*K**2         s
          !
          !     gwrk(1)   difference of pseudo v wind ( km is destroyed )
          !     gwrk(5)   difference of pseudo u wind ( b  is destroyed )
          !     
          Pbl_DifVzn(i,k)=gu(i,k)-gu(i,k+1)
          Pbl_DifVmd(i,k)=gv(i,k)-gv(i,k+1)
       END DO
    END DO
    DO i = 1, nCols
       Pbl_KMbyDZ2   (i,1)=0.0_r8
       Pbl_KMbyDZ1   (i,1)=1.0_r8 + Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2   (i,1)=       - Pbl_KHbyDZ(i,1)
       !                          --             --
       ! DU       d(w'u')      d |          d U    |       m
       !------ = ------- =    ---| - Km *  ------  |  = --------
       ! Dt       dz           dz|          dz     |     s * s
       !                          --             --
       Pbl_TendU(i,1)=-twodti*Pbl_KHbyDZ(i,1)*Pbl_DifVzn(i,1)!(1/m)*(m/s)
       !
       !                     1          1           m           m
       !Pbl_TendU(i,1) = - ------ * -------- *  -------- =   --------
       !                     s          s           s         s * s
       !
       Pbl_TendV(i,1)=-twodti*Pbl_KHbyDZ(i,1)*Pbl_DifVmd(i,1)!m/s**2

       Pbl_KMbyDZ2   (i,kMax)=       - Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1   (i,kMax)=1.0_r8 + Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2   (i,kMax)=0.0_r8
       Pbl_TendU(i,kMax)=twodti*Pbl_KMbyDZ(i,kMax)*Pbl_DifVzn(i,kMax-1)
       Pbl_TendV(i,kMax)=twodti*Pbl_KMbyDZ(i,kMax)*Pbl_DifVmd(i,kMax-1)
    END DO
    DO k = 2, kMax-1
       DO i = 1, nCols
          Pbl_KMbyDZ2   (i,k)=      -Pbl_KMbyDZ(i,k)                   !1/(sec * K)
          Pbl_KMbyDZ1   (i,k)=1.0_r8+Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)   !1/(sec * K)
          Pbl_KHbyDZ2   (i,k)=      -Pbl_KHbyDZ(i,k)                   !1/(sec * K)
          Pbl_TendU(i,k)=(Pbl_KMbyDZ(i,k)*Pbl_DifVzn(i,k-1)  - Pbl_KHbyDZ(i,k)*Pbl_DifVzn(i,k  )) * twodti
          Pbl_TendV(i,k)=(Pbl_KMbyDZ(i,k)*Pbl_DifVmd(i,k-1)  - Pbl_KHbyDZ(i,k)*Pbl_DifVmd(i,k  )) * twodti
       END DO
    END DO
    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          !
          !                                - Pbl_KMbyDZ_1(i,k) 
          !Pbl_KHbyDZ2   (i,k)=-------------------------------------------------
          !                     1.0 + Pbl_KMbyDZ_1(i,k+1) + Pbl_KMbyDZ_2(i,k+1)
          !
          Pbl_KHbyDZ2   (i,k)=Pbl_KHbyDZ2 (i,k)/Pbl_KMbyDZ1(i,k+1)
          !
          !                                                                         Pbl_KMbyDZ_1(i,k)*Pbl_KMbyDZ_2(i,k+1)   
          !Pbl_KMbyDZ1   (i,k)=1.0_r8 + Pbl_KMbyDZ_1(i,k) + Pbl_KMbyDZ_2(i,k) - -------------------------------------------------
          !                                                                      1.0 + Pbl_KMbyDZ_1(i,k+1) + Pbl_KMbyDZ_2(i,k+1)
          !
          Pbl_KMbyDZ1   (i,k)=Pbl_KMbyDZ1 (i,k) - Pbl_KHbyDZ2(i,k)*Pbl_KMbyDZ2(i,k+1)
          !
          !
          !
          !
          !
          Pbl_TendU(i,k)=Pbl_TendU(i,k) - Pbl_KHbyDZ2(i,k)*Pbl_TendU(i,k+1)
          Pbl_TendV(i,k)=Pbl_TendV(i,k) - Pbl_KHbyDZ2(i,k)*Pbl_TendV(i,k+1)
       END DO
    END DO
    DO i = 1, ncols
       gmu(i,1,1)=Pbl_KMbyDZ2 (i,1)
       gmu(i,1,2)=Pbl_KMbyDZ1 (i,1)
       gmu(i,1,3)= 1.0_r8*Pbl_TendU   (i,1)
       gmu(i,1,4)= 1.0_r8*Pbl_TendV   (i,1)
    END DO
  END SUBROUTINE VDIFV
  !---------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------
  SUBROUTINE VDIFH(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,a,b,gq,gmq)
    !     ***************************************************************
    !     *                                                             *
    !     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
    !     *                                                             *
    !     ***************************************************************
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_CoefKh(nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gq        (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmq       (nCols,kMax,5+nClass+nAeros)

    REAL(KIND=r8) :: Pbl_DifQms(nCols,kMax) 
    REAL(KIND=r8) :: Pbl_CoefKh2(nCols,kMax)

    REAL(KIND=r8) :: Pbl_KHbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_KMbyDZ(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_TendQ(nCols,kMax) !1/(sec * K)

    INTEGER       :: i
    INTEGER       :: k

    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKh2(i,k  )=Pbl_CoefKh(i,k)*Pbl_ITemp(i,k)
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKh2(i,k)
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKh2(i,k)
          !     
          !     Pbl_DifQms(1)   difference of specific humidity
          !     
          Pbl_DifQms(i,k)=gq(i,k)-gq(i,k+1)  
       END DO
    END DO
    DO i = 1, ncols
       Pbl_KMbyDZ2(i,1)=0.0_r8
       Pbl_KMbyDZ1(i,1)=1.0_r8  + Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2(i,1)=        - Pbl_KHbyDZ(i,1)
       Pbl_TendQ  (i,1)=-twodti * Pbl_KHbyDZ(i,1) * Pbl_DifQms(i,1)

       Pbl_KMbyDZ2(i,kMax)=        - Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1(i,kMax)=1.0_r8  + Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2(i,kMax)=0.0_r8
       Pbl_TendQ  (i,kMax)= twodti * Pbl_KMbyDZ(i,kMax) * Pbl_DifQms(i,kMax-1)
    END DO
    DO k = 2, kmax-1
       DO i = 1, ncols
          Pbl_KMbyDZ2(i,k)=       - Pbl_KMbyDZ(i,k)
          Pbl_KMbyDZ1(i,k)=1.0_r8 + Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)
          Pbl_KHbyDZ2(i,k)=       - Pbl_KHbyDZ(i,k)
          Pbl_TendQ(i,k)=(Pbl_KMbyDZ(i,k) * Pbl_DifQms(i,k-1)- &
               Pbl_KHbyDZ(i,k) * Pbl_DifQms(i,k  ))*twodti
       END DO
    END DO
    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          Pbl_KHbyDZ2(i,k)=Pbl_KHbyDZ2(i,k) / Pbl_KMbyDZ1(i,k+1)
          Pbl_KMbyDZ1(i,k)=Pbl_KMbyDZ1(i,k) - Pbl_KHbyDZ2(i,k)*Pbl_KMbyDZ2(i,k+1)
          Pbl_TendQ  (i,k)=Pbl_TendQ  (i,k) - Pbl_KHbyDZ2(i,k)*Pbl_TendQ  (i,k+1)
       END DO
    END DO
    DO i = 1, ncols
       gmq(i,1,1)=Pbl_KMbyDZ2(i,1)
       gmq(i,1,2)=Pbl_KMbyDZ1(i,1)
       gmq(i,1,3)=Pbl_TendQ  (i,1)
    END DO
  END SUBROUTINE VDIFH
  SUBROUTINE VDIFHM(kMax,nCols,twodti,Pbl_ITemp,Pbl_CoefKh,a,b,gice ,gliq,gmq,gvar)
    !     ***************************************************************
    !     *                                                             *
    !     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
    !     *                                                             *
    !     ***************************************************************
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_CoefKh(nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gice      (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gliq      (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmq       (nCols,kMax,5+nClass+nAeros)
    REAL(KIND=r8),OPTIONAL, INTENT(in   ) :: gvar(nCols,kMax,nClass+nAeros)

    REAL(KIND=r8) :: Pbl_DifQms (nCols,kMax,3+nClass+nAeros) 
    REAL(KIND=r8) :: Pbl_TendQ  (nCols,kMax,3+nClass+nAeros) !1/(sec * K)

    REAL(KIND=r8) :: Pbl_KHbyDZ (nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_KMbyDZ (nCols,kMax) !1/(sec * K)
    REAL(KIND=r8) :: Pbl_CoefKh2(nCols,kMax)
    INTEGER , PARAMETER      :: ice=2
    INTEGER , PARAMETER      :: iliq=3

    INTEGER       :: i
    INTEGER       :: k
    INTEGER       :: kk

    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKh2(i,k  )=Pbl_CoefKh(i,k)*Pbl_ITemp  (i,k)
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKh2(i,k)
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKh2(i,k)
          !     
          !     Pbl_DifQms(1)   difference of specific humidity
          !     
          Pbl_DifQms(i,k,ice ) = gice(i,k)-gice(i,k+1)  
          Pbl_DifQms(i,k,iliq) = gliq(i,k)-gliq(i,k+1)  
          !t1        (i,k)=g1(i,k)-g1(i,k+1)
       END DO
    END DO
    IF (PRESENT(gvar)) THEN
       DO kk=1,nClass+nAeros
          DO k = 1, kMax-1
             DO i = 1, nCols
                !     
                !     Pbl_DifQms(1)   difference of specific humidity
                !     
                Pbl_DifQms(i,k,3+kk) = gvar(i,k,kk)-gvar(i,k+1,kk)  
                !t1        (i,k)=g1(i,k)-g1(i,k+1)
             END DO
          END DO
       END DO
    END IF
!-------------------------------------------------------------------------------------------------
    DO i = 1, ncols
       Pbl_KMbyDZ2(i,1)=0.0_r8
       Pbl_KMbyDZ1(i,1)=1.0_r8  + Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2(i,1)=        - Pbl_KHbyDZ(i,1)
       Pbl_TendQ  (i,1,ice )=-twodti * Pbl_KHbyDZ(i,1) * Pbl_DifQms(i,1,ice )
       Pbl_TendQ  (i,1,iliq)=-twodti * Pbl_KHbyDZ(i,1) * Pbl_DifQms(i,1,iliq)
       !t1        (i,1)=-twodti * Pbl_KHbyDZ(i,1) * t1(i,1)

       Pbl_KMbyDZ2(i,kMax)=        - Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1(i,kMax)=1.0_r8  + Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2(i,kMax)=0.0_r8
       Pbl_TendQ  (i,kMax,ice )= twodti * Pbl_KMbyDZ(i,kMax) * Pbl_DifQms(i,kMax-1,ice )
       Pbl_TendQ  (i,kMax,iliq)= twodti * Pbl_KMbyDZ(i,kMax) * Pbl_DifQms(i,kMax-1,iliq)

       !t2(i,kMax)        = twodti * Pbl_KMbyDZ(i,kMax) * t2(i,kMax-1)

    END DO
    IF (PRESENT(gvar)) THEN
       DO kk=1,nClass+nAeros
          DO i = 1, ncols
             Pbl_TendQ  (i,   1,3+kk)=-twodti * Pbl_KHbyDZ(i,1) * Pbl_DifQms(i,1,3+kk)
             !t1        (i,1)=-twodti * Pbl_KHbyDZ(i,1) * t1(i,1) 

             Pbl_TendQ  (i,kMax,3+kk)= twodti * Pbl_KMbyDZ(i,kMax) * Pbl_DifQms(i,kMax-1,3+kk)
 
             !t2(i,kMax)        = twodti * Pbl_KMbyDZ(i,kMax) * t2(i,kMax-1)
          END DO
       END DO
    END IF
!-------------------------------------------------------------------------------------------------

    DO k = 2, kmax-1
       DO i = 1, ncols
          Pbl_KMbyDZ2(i,k)=       - Pbl_KMbyDZ(i,k)
          Pbl_KMbyDZ1(i,k)=1.0_r8 + Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)
          Pbl_KHbyDZ2(i,k)=       - Pbl_KHbyDZ(i,k)
          Pbl_TendQ(i,k,ice )= (Pbl_KMbyDZ(i,k) * Pbl_DifQms(i,k-1,ice )- Pbl_KHbyDZ(i,k) * Pbl_DifQms(i,k,ice   ))*twodti
          Pbl_TendQ(i,k,iliq)= (Pbl_KMbyDZ(i,k) * Pbl_DifQms(i,k-1,iliq)- Pbl_KHbyDZ(i,k) * Pbl_DifQms(i,k,iliq  ))*twodti

         ! t1(i,k)       = (Pbl_KMbyDZ(i,k) * t1(i,k-1)        - Pbl_KHbyDZ(i,k) * t1(i,k  )        )*twodti

       END DO
    END DO
    IF (PRESENT(gvar)) THEN
       DO kk=1,nClass+nAeros
          DO k = 2, kmax-1
             DO i = 1, ncols
                Pbl_TendQ(i,k,3+kk )= (Pbl_KMbyDZ(i,k) * Pbl_DifQms(i,k-1,3+kk )- Pbl_KHbyDZ(i,k) * Pbl_DifQms(i,k,3+kk   ))*twodti

                ! t1(i,k)       = (Pbl_KMbyDZ(i,k) * t1(i,k-1)        - Pbl_KHbyDZ(i,k) * t1(i,k  )        )*twodti

             END DO
          END DO
       END DO
    END IF
!-------------------------------------------------------------------------------------------------

    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          Pbl_KHbyDZ2(i,k)=Pbl_KHbyDZ2(i,k) / Pbl_KMbyDZ1(i,k+1)
          Pbl_KMbyDZ1(i,k)=Pbl_KMbyDZ1(i,k) - Pbl_KHbyDZ2(i,k)*Pbl_KMbyDZ2(i,k+1)
          Pbl_TendQ  (i,k,ice )=Pbl_TendQ  (i,k,ice ) - Pbl_KHbyDZ2(i,k)*Pbl_TendQ  (i,k+1,ice )
          Pbl_TendQ  (i,k,iliq)=Pbl_TendQ  (i,k,iliq) - Pbl_KHbyDZ2(i,k)*Pbl_TendQ  (i,k+1,iliq)

         !t1         (i,k)=t1         (i,k) - Pbl_KHbyDZ2(i,k)*t1         (i,k+1)

       END DO
    END DO
    IF (PRESENT(gvar)) THEN
       DO kk=1,nClass+nAeros
          DO k = kmax-1, 1, -1
             DO i = 1, ncols
                Pbl_TendQ  (i,k,3+kk)=Pbl_TendQ  (i,k,3+kk) - Pbl_KHbyDZ2(i,k)*Pbl_TendQ  (i,k+1,3+kk)
                !t1         (i,k)=t1         (i,k) - Pbl_KHbyDZ2(i,k)*t1         (i,k+1)
             END DO
          END DO
       END DO
    END IF
!-------------------------------------------------------------------------------------------------

    DO k = 1, kmax
       DO i = 1, ncols
          gmq(i,k,1)=Pbl_KMbyDZ2(i,k)
          gmq(i,k,2)=Pbl_KMbyDZ1(i,k)
          gmq(i,k,3+ice -1)=Pbl_TendQ  (i,k,ice ) 
          gmq(i,k,3+iliq-1)=Pbl_TendQ  (i,k,iliq)
       END DO
    END DO
    IF (PRESENT(gvar)) THEN
       DO kk=1,nClass+nAeros
          DO k = 1, kmax
             DO i = 1, ncols
                gmq(i,k,6+kk -1)=Pbl_TendQ  (i,k,3+kk ) 
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE VDIFHM
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE VDIFT(kMax,nCols,twodti,sigr,sigriv,Pbl_ITemp,Pbl_CoefKh,a,b,gt,gmt)
    !     ***************************************************************
    !     *                                                             *
    !     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
    !     *                                                             *
    !     ***************************************************************
    INTEGER      ,    INTENT(in   ) :: nCols
    INTEGER      ,    INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: twodti
    REAL(KIND=r8),    INTENT(in   ) :: a         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: b         (nCols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: sigr      (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: sigriv    (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_ITemp (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: Pbl_CoefKh(nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gt        (nCols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: gmt       (nCols,kMax,3)

    REAL(KIND=r8) :: Pbl_CoefKh2(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KHbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ(nCols,kMax) !1/(sec * K) 
    REAL(KIND=r8) :: Pbl_KMbyDZ2(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KMbyDZ1(nCols,kMax)
    REAL(KIND=r8) :: Pbl_KHbyDZ2(nCols,kMax)
    REAL(KIND=r8) :: Pbl_TendT(nCols,kMax)
    INTEGER       :: i
    INTEGER       :: k
    !     
    !     sensible heat diffusion
    !     
    !   sigk  (   k)=       sl(k)**akappa
    !   sigkiv(   k)=1.0_r8/sl(k)**akappa
    !   sigr  (   k)=sigk(k  )*sigkiv(k+1) =(sl(k  )**akappa)*(1.0_r8/sl(k+1)**akappa)
    !   sigriv( k+1)=sigk(k+1)*sigkiv(k  ) =(sl(k+1)**akappa)*(1.0_r8/sl(k  )**akappa)
    !
    DO k = 1, kMax-1
       DO i = 1, nCols
          Pbl_CoefKh2(i,k  )=Pbl_CoefKh(i,k)*Pbl_ITemp(i,k)
          Pbl_KHbyDZ(i,k  )=a(i,k  )*Pbl_CoefKh2(i,k)
          Pbl_KMbyDZ(i,k+1)=b(i,k+1)*Pbl_CoefKh2(i,k)
       END DO
    END DO
    DO i = 1, nCols
       Pbl_KMbyDZ2(i,1)=  0.0_r8
       Pbl_KMbyDZ1(i,1)=  1.0_r8+Pbl_KHbyDZ(i,1)
       Pbl_KHbyDZ2(i,1)=-sigr(i,1)*Pbl_KHbyDZ(i,1)
       Pbl_TendT  (i,1)=-Pbl_KHbyDZ(i,1)*(gt(i,1)-sigr(i,1)*gt(i,1+1))*twodti

       Pbl_KMbyDZ2(i,kMax)=-sigriv(i,kMax)*Pbl_KMbyDZ(i,kMax)
       Pbl_KMbyDZ1(i,kMax)=    1.0_r8+Pbl_KMbyDZ(i,kMax)
       Pbl_KHbyDZ2(i,kMax)=0.0_r8
       Pbl_TendT  (i,kMax)=twodti*Pbl_KMbyDZ(i,kMax)*(sigriv(i,kMax)*gt(i,kMax-1)-gt(i,kMax))
    END DO
    DO k = 2, kMax-1
       DO i = 1, nCols
          Pbl_KMbyDZ2(i,k)=-sigriv(i,k)*Pbl_KMbyDZ(i,k)
          Pbl_KMbyDZ1(i,k)=1.0_r8+Pbl_KHbyDZ(i,k)+Pbl_KMbyDZ(i,k)
          Pbl_KHbyDZ2(i,k)=-sigr  (i,k)*Pbl_KHbyDZ(i,k)
          Pbl_TendT  (i,k)=( Pbl_KMbyDZ(i,k)*(sigriv(i,k)*gt(i,k-1) - gt(i,k))&
               -Pbl_KHbyDZ(i,k)*(gt(i,k)- sigr(i,k)*gt(i,k+1))  )*twodti
       END DO
    END DO

    DO k = kmax-1, 1, -1
       DO i = 1, ncols
          Pbl_KHbyDZ2(i,k)=Pbl_KHbyDZ2(i,k)/Pbl_KMbyDZ1(i,k+1)
          Pbl_KMbyDZ1(i,k)=Pbl_KMbyDZ1(i,k)-Pbl_KHbyDZ2(i,k  )*Pbl_KMbyDZ2(i,k+1)
          Pbl_TendT  (i,k)=Pbl_TendT  (i,k)-Pbl_KHbyDZ2(i,k  )*Pbl_TendT  (i,k+1)
       END DO
    END DO

    DO i = 1, ncols
       gmt(i,1,1)=Pbl_KMbyDZ2(i,1)
       gmt(i,1,2)=Pbl_KMbyDZ1(i,1)
       gmt(i,1,3)=Pbl_TendT  (i,1)
    END DO
  END SUBROUTINE VDIFT


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  SUBROUTINE SFCDIF(NTSD,SEAMASK,THS,QS,PSFC                       &
       &                 ,UZ0,VZ0,TZ0,THZ0,QZ0                           &
       &                 ,USTAR,Z0,CT,RLMO,AKMS,AKHS,PBLH,WETM           &
       &                 ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC,QGH,CPM &
       &                 ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW        &
       &                 ,ZSL,PLOW                                       &
       &                 ,U10,V10,TH02,TH10,Q02,Q10,PSHLTR,WSTAR,RIB)
    !     ****************************************************************
    !     *                                                              *
    !     *                       SURFACE LAYER                          *
    !     *                                                              *
    !     ****************************************************************
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------------
    !
    INTEGER,INTENT(IN) :: NTSD
    !
    REAL(KIND=r8),INTENT(IN) :: CWMLOW,PBLH,PLOW,QLOW,PSFC,SEAMASK            &
         &                  ,THELOW,THLOW,THS,TLOW,TZ0,ULOW,VLOW,WETM,ZSL
    !
    REAL(KIND=r8),INTENT(OUT) :: CHS,CHS2,CPM,CQS2,CT,FLHC,FLQC,FLX_LH,HFX    &
         &                   ,PSHLTR,Q02,Q10,QFX,QGH,RLMO,TH02,TH10,U10,V10 &
         &                   ,WSTAR,RIB
    !
    REAL(KIND=r8),INTENT(INOUT) :: AKHS,AKMS,QZ0,THZ0,USTAR,UZ0,VZ0,Z0,QS
    !----------------------------------------------------------------------
    !***
    !***  LOCAL VARIABLES
    !***
    INTEGER :: ITR,K
    !
    REAL(KIND=r8) :: A,B,BTGH,BTGX,CXCHL,CXCHS,DTHV,DU2,ELFC,FCT              &
         &       ,HLFLX,HSFLX,HV,PSH02,PSH10,PSHZ,PSHZL,PSM10,PSMZ,PSMZL   &
         &       ,RDZ,RDZT                              &
         &       ,RLOGT,RLOGU,RWGH,RZ,RZST,RZSU,SIMH,SIMM,TEM,THM          &
         &       ,UMFLX,USTARK,VMFLX,WGHT,WGHTT,WGHTQ,WSTAR2               &
         &       ,ZETALT,ZETALU          &
         &       ,ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU
    !
    !***  DIAGNOSTICS
    !
    REAL(KIND=r8) :: AKHS02,AKHS10,AKMS10,EKMS10          &
         &       ,RLNT02,RLNT10,RLNU10,SIMH02,SIMH10,SIMM10        &
         &       ,TERM1,RLOW,U10E,V10E            &
         &       ,ZT02,ZT10,ZTAT02,ZTAT10   &
         &       ,ZTAU10,ZU10,ZUUZ
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    RDZ=1.0_r8/ZSL
    !
    !        m**2      1
    !CXCHL= ------- * -----
    !        s**2       m
    !
    CXCHL=EXCML*RDZ ! m/s**2
    !
    !
    !
    CXCHS=EXCMS*RDZ ! m/s**2
    !
    !
    BTGX=G/THLOW
    ELFC=VKARMAN*BTGX
    ! 
    IF(PBLH>1000.0_r8)THEN
       BTGH=BTGX*PBLH
    ELSE
       BTGH=BTGX*1000.0_r8
    ENDIF 
    !
    !----------------------------------------------------------------------
    !
    !***  SEA POINTS
    !
    !----------------------------------------------------------------------
    !
    IF(SEAMASK>0.5_r8)THEN 
       !
       ZT=Z0
       !----------------------------------------------------------------------
       DO ITR=1,ITRMX
          !----------------------------------------------------------------------
          
          !Z0=MAX(USTFC*USTAR*USTAR,1.59E-5_r8)
          !
          !***  VISCOUS SUBLAYER, JANJIC MWR 1994
          !
          !----------------------------------------------------------------------
          IF(USTAR<USTC)THEN
             !----------------------------------------------------------------------
             !
             IF(USTAR<USTR)THEN
                !
                IF(NTSD==1)THEN
                   AKMS=CXCHL
                   AKHS=CXCHS
                   QS=QLOW
                ENDIF
                !
                ZU=FZU1*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
                WGHT=AKMS*ZU*RVISC
                RWGH=WGHT/(WGHT+1.0_r8)
                UZ0=(ULOW*RWGH+UZ0)*0.5_r8
                VZ0=(VLOW*RWGH+VZ0)*0.5_r8
                !
                ZT=FZT1*ZU
                ZQ=FZQ1*ZT
                WGHTT=AKHS*ZT*RTVISC
                WGHTQ=AKHS*ZQ*RQVISC
                !
                IF(NTSD>1)THEN
                   THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.0_r8)+THZ0)*0.5_r8
                   QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.0_r8)+QZ0)*0.5_r8
                ELSE
                   THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.0_r8)
                   QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.0_r8)
                ENDIF
                !
             ENDIF
             !
             IF(USTAR>=USTR.AND.USTAR<USTC)THEN
                ZU=Z0
                UZ0=0.0_r8
                VZ0=0.0_r8
                !
                ZT=FZT2*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
                ZQ=FZQ2*ZT
                WGHTT=AKHS*ZT*RTVISC
                WGHTQ=AKHS*ZQ*RQVISC
                !
                IF(NTSD>1)THEN
                   THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.0_r8)+THZ0)*0.5_r8
                   QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.0_r8)+QZ0)*0.5_r8
                ELSE
                   THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.0_r8)
                   QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.0_r8)
                ENDIF
                !
             ENDIF
             !----------------------------------------------------------------------
          ELSE
             !----------------------------------------------------------------------
             ZU=Z0
             UZ0=0.0_r8
             VZ0=0.0_r8
             !
             ZT=MAX(EXP(ZILFC*SQRT(USTAR*ZU))*ZU,EPSZT)
             THZ0=THS
             !
             ZQ=Z0
             QZ0=QS
             !----------------------------------------------------------------------
          ENDIF
          !----------------------------------------------------------------------
          TEM=(TLOW+TZ0)*0.5_r8
          THM=(THELOW+THZ0)*0.5_r8
          !
          A=THM*P608
          B=(ELOCP/TEM-1.0_r8-P608)*THM
          !
          DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5_r8*P608)+1.0_r8)        &
               &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B)
          !
          DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
          RIB=BTGX*DTHV*ZSL/DU2
          !----------------------------------------------------------------------
          !         IF(RIB>=RIC)THEN
          !----------------------------------------------------------------------
          !           AKMS=MAX( VISC*RDZ,CXCHL)
          !           AKHS=MAX(TVISC*RDZ,CXCHS)
          !----------------------------------------------------------------------
          !         ELSE  !  turbulent branch
          !----------------------------------------------------------------------
          ZSLU=ZSL+ZU
          ZSLT=ZSL+ZT
          !
          RZSU=ZSLU/ZU
          RZST=ZSLT/ZT
          !
          RLOGU=LOG(RZSU)
          RLOGT=LOG(RZST)
          !
          !----------------------------------------------------------------------
          !***  1.0_r8/MONIN-OBUKHOV LENGTH
          !----------------------------------------------------------------------
          !
          RLMO=ELFC*AKHS*DTHV/USTAR**3
          !
          ZETALU=ZSLU*RLMO
          ZETALT=ZSLT*RLMO
          ZETAU=ZU*RLMO
          ZETAT=ZT*RLMO
          !
          ZETALU=MIN(MAX(ZETALU,ZTMIN1),ZTMAX1)
          ZETALT=MIN(MAX(ZETALT,ZTMIN1),ZTMAX1)
          ZETAU=MIN(MAX(ZETAU,ZTMIN1/RZSU),ZTMAX1/RZSU)
          ZETAT=MIN(MAX(ZETAT,ZTMIN1/RZST),ZTMAX1/RZST)
          !
          !----------------------------------------------------------------------
          !***   WATER FUNCTIONS
          !----------------------------------------------------------------------
          !
          RZ=(ZETAU-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSMZ=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
          !
          RZ=(ZETALU-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSMZL=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
          !
          SIMM=PSMZL-PSMZ+RLOGU
          !
          RZ=(ZETAT-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSHZ=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
          !
          RZ=(ZETALT-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSHZL=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
          !
          SIMH=(PSHZL-PSHZ+RLOGT)*FH01
          !----------------------------------------------------------------------
          USTARK=USTAR*VKARMAN
          AKMS=MAX(USTARK/SIMM,CXCHL)
          AKHS=MAX(USTARK/SIMH,CXCHS)
          !
          !----------------------------------------------------------------------
          !***  BELJAARS CORRECTION FOR USTAR
          !----------------------------------------------------------------------
          !
          WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2.0_r8/3.0_r8)
          USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
          !----------------------------------------------------------------------
          !         ENDIF  !  End of turbulent branch
          !----------------------------------------------------------------------
          !
       ENDDO  !  End of the iteration loop over sea points
       !
       !----------------------------------------------------------------------
       !
       !***  LAND POINTS
       !
       !----------------------------------------------------------------------
       !
    ELSE  
       !
       !----------------------------------------------------------------------
       !
       IF(NTSD==1)THEN
          QS=QLOW
       ENDIF
       !
       ZU=Z0
       UZ0=0.0_r8
       VZ0=0.0_r8
       !
       ZT=ZU*ZTFC
       THZ0=THS
       QZ0=QS
       !----------------------------------------------------------------------
       TEM=(TLOW+TZ0)*0.5_r8
       THM=(THELOW+THZ0)*0.5_r8
       !
       A=THM*P608
       B=(ELOCP/TEM-1.0_r8-P608)*THM
       !
       DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5_r8*P608)+1.0_r8)          &
            &       +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B)
       !
       DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
       RIB=BTGX*DTHV*ZSL/DU2
       !----------------------------------------------------------------------
       !       IF(RIB>=RIC)THEN
       !         AKMS=MAX( VISC*RDZ,CXCHL)
       !         AKHS=MAX(TVISC*RDZ,CXCHS)
       !----------------------------------------------------------------------
       !       ELSE  !  Turbulent branch
       !----------------------------------------------------------------------
       ZSLU=ZSL+ZU
       RZSU=ZSLU/ZU
       RLOGU=LOG(RZSU)
       ZSLT=ZSL+ZU
       !----------------------------------------------------------------------
       !
       DO ITR=1,ITRMX
          !
          !----------------------------------------------------------------------
          !***  ZILITINKEVITCH FIX FOR ZT
          !----------------------------------------------------------------------
          !
          ZT=MAX(EXP(ZILFC*SQRT(USTAR*ZU))*ZU,EPSZT)
          !zj         ZT=EXP(ZILFC*SQRT(USTAR*ZU))*ZU
          RZST=ZSLT/ZT
          RLOGT=LOG(RZST)
          !
          !----------------------------------------------------------------------
          !***  1.0_r8/MONIN-OBUKHOV LENGTH-SCALE
          !----------------------------------------------------------------------
          !
          RLMO=ELFC*AKHS*DTHV/USTAR**3
          ZETALU=ZSLU*RLMO
          ZETALT=ZSLT*RLMO
          ZETAU=ZU*RLMO
          ZETAT=ZT*RLMO
          !
          ZETALU=MIN(MAX(ZETALU,ZTMIN2),ZTMAX2)
          ZETALT=MIN(MAX(ZETALT,ZTMIN2),ZTMAX2)
          ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
          ZETAT=MIN(MAX(ZETAT,ZTMIN2/RZST),ZTMAX2/RZST)
          !
          !----------------------------------------------------------------------
          !***  LAND FUNCTIONS
          !----------------------------------------------------------------------
          !
          RZ=(ZETAU-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSMZ=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
          !
          RZ=(ZETALU-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSMZL=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
          !
          SIMM=PSMZL-PSMZ+RLOGU
          !
          RZ=(ZETAT-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSHZ=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
          !
          RZ=(ZETALT-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K,KIND=r8)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSHZL=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
          !
          SIMH=(PSHZL-PSHZ+RLOGT)*FH02
          !----------------------------------------------------------------------
          USTARK=USTAR*VKARMAN
          AKMS=MAX(USTARK/SIMM,CXCHL)
          AKHS=MAX(USTARK/SIMH,CXCHS)
          !
          !----------------------------------------------------------------------
          !***  BELJAARS CORRECTION FOR USTAR
          !----------------------------------------------------------------------
          !
          WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2.0_r8/3.0_r8)
          USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
          !
          !----------------------------------------------------------------------
       ENDDO  !  End of iMaxration for land points
       !----------------------------------------------------------------------
       !
       !       ENDIF  !  End of turbulant branch over land
       !
       !----------------------------------------------------------------------
       !
    ENDIF  !  End of land/sea branch
    !
    !----------------------------------------------------------------------
    !***  COUNTERGRADIENT FIX
    !----------------------------------------------------------------------
    HV=-AKHS*DTHV
    IF(HV>0.0_r8)THEN
       FCT=-10.0_r8*(BTGX)**(-1.0_r8/3.0_r8)
       CT=FCT*(HV/(PBLH*PBLH))**(2.0_r8/3.0_r8)
    ELSE
       CT=0.0_r8
    ENDIF
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !***  THE FOLLOWING DIAGNOSTIC BLOCK PRODUCES 2-m and 10-m VALUES
    !***  FOR TEMPERATURE, MOISTURE, AND WINDS.  IT IS DONE HERE SINCE
    !***  THE VARIOUS QUANTITIES NEEDED FOR THE COMPUTATION ARE LOST
    !***  UPON EXIT FROM THE ROTUINE.
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !
    WSTAR=SQRT(WSTAR2)/WWST
    !
    !                   --          --
    !             |       d T    |
    !        -(T'W') = |KH * -------- |
    !                  |       d z    |
    !              --          --
    !                                

    !                                 --              --
    !           d -(T'W')           d        |         d T        | 
    !        ------------- = --------|KH *  -------- |
    !           dt                   dt        |         d z        |
    !                                 --              --

    UMFLX = AKMS * (ULOW -UZ0 )
    VMFLX = AKMS * (VLOW -VZ0 )
    HSFLX = AKHS * (THLOW-THZ0)
    HLFLX = AKHS * (QLOW -QZ0 )
    !----------------------------------------------------------------------
    !     IF(RIB>=RIC)THEN
    !----------------------------------------------------------------------
    !       IF(SEAMASK>0.5_r8)THEN
    !         AKMS10=MAX( VISC/10.0_r8,CXCHS)
    !         AKHS02=MAX(TVISC/02.0_r8,CXCHS)
    !         AKHS10=MAX(TVISC/10.0_r8,CXCHS)
    !       ELSE
    !         AKMS10=MAX( VISC/10.0_r8,CXCHL)
    !         AKHS02=MAX(TVISC/02.0_r8,CXCHS)
    !         AKHS10=MAX(TVISC/10.0_r8,CXCHS)
    !       ENDIF
    !----------------------------------------------------------------------
    !     ELSE
    !----------------------------------------------------------------------
    ZU10=ZU+10.0_r8
    ZT02=ZT+02.0_r8
    ZT10=ZT+10.0_r8
    !
    RLNU10=LOG(ZU10/ZU)
    RLNT02=LOG(ZT02/ZT)
    RLNT10=LOG(ZT10/ZT)
    !
    ZTAU10=ZU10*RLMO
    ZTAT02=ZT02*RLMO
    ZTAT10=ZT10*RLMO
    !
    !----------------------------------------------------------------------
    !***  SEA
    !----------------------------------------------------------------------
    !
    IF(SEAMASK>0.5_r8)THEN
       !
       !----------------------------------------------------------------------
       ZTAU10=MIN(MAX(ZTAU10,ZTMIN1),ZTMAX1)
       ZTAT02=MIN(MAX(ZTAT02,ZTMIN1),ZTMAX1)
       ZTAT10=MIN(MAX(ZTAT10,ZTMIN1),ZTMAX1)
       !----------------------------------------------------------------------
       RZ=(ZTAU10-ZTMIN1)/DZETA1
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSM10=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
       !
       SIMM10=PSM10-PSMZ+RLNU10
       !
       RZ=(ZTAT02-ZTMIN1)/DZETA1
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSH02=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
       !
       SIMH02=(PSH02-PSHZ+RLNT02)*FH01
       !
       RZ=(ZTAT10-ZTMIN1)/DZETA1
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSH10=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
       !
       SIMH10=(PSH10-PSHZ+RLNT10)*FH01
       !
       AKMS10=MAX(USTARK/SIMM10,CXCHL)
       AKHS02=MAX(USTARK/SIMH02,CXCHS)
       AKHS10=MAX(USTARK/SIMH10,CXCHS)
       !
       !----------------------------------------------------------------------
       !***  LAND
       !----------------------------------------------------------------------
       !
    ELSE
       !
       !----------------------------------------------------------------------
       ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
       ZTAT02=MIN(MAX(ZTAT02,ZTMIN2),ZTMAX2)
       ZTAT10=MIN(MAX(ZTAT10,ZTMIN2),ZTMAX2)
       !----------------------------------------------------------------------
       RZ=(ZTAU10-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
       !
       SIMM10=PSM10-PSMZ+RLNU10
       !
       RZ=(ZTAT02-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSH02=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
       !
       SIMH02=(PSH02-PSHZ+RLNT02)*FH02
       !
       RZ=(ZTAT10-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSH10=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
       !
       SIMH10=(PSH10-PSHZ+RLNT10)*FH02
       !
       AKMS10=MAX(USTARK/SIMM10,CXCHL)
       AKHS02=MAX(USTARK/SIMH02,CXCHS)
       AKHS10=MAX(USTARK/SIMH10,CXCHS)
       !----------------------------------------------------------------------
    ENDIF
    !----------------------------------------------------------------------
    !     ENDIF
    !----------------------------------------------------------------------
    U10 =UMFLX/AKMS10+UZ0
    V10 =VMFLX/AKMS10+VZ0
    TH02=HSFLX/AKHS02+THZ0
    TH10=HSFLX/AKHS10+THZ0
    Q02 =HLFLX/AKHS02+QZ0
    Q10 =HLFLX/AKHS10+QZ0
    TERM1=-0.068283_r8/TLOW
    PSHLTR=PSFC*EXP(TERM1)
    !
    !----------------------------------------------------------------------
    !***  COMPUTE "EQUIVALENT" Z0 TO APPROXIMATE LOCAL SHELTER READINGS.
    !----------------------------------------------------------------------
    !
    U10E=U10
    V10E=V10
    !
    IF(SEAMASK<0.5_r8)THEN
       !
       ZUUZ=MIN(0.5_r8*ZU,0.1_r8)
       ZU=MAX(0.1_r8*ZU,ZUUZ)
       !
       ZU10=ZU+10.0_r8
       RZSU=ZU10/ZU
       RLNU10=LOG(RZSU)
       !
       ZETAU=ZU*RLMO
       ZTAU10=ZU10*RLMO
       !
       ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
       ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
       !
       RZ=(ZTAU10-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K,KIND=r8)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
       SIMM10=PSM10-PSMZ+RLNU10
       EKMS10=MAX(USTARK/SIMM10,CXCHL)
       !
       U10E=UMFLX/EKMS10+UZ0
       V10E=VMFLX/EKMS10+VZ0
       !
    ENDIF
    !
    U10=U10E
    V10=V10E
    !
    !----------------------------------------------------------------------
    !***  SET OTHER WRF DRIVER ARRAYS
    !----------------------------------------------------------------------
    !
    RLOW=PLOW/(R_D*TLOW)
    CHS=AKHS
    CHS2=AKHS02
    CQS2=AKHS02
    HFX=-RLOW*CP*HSFLX
    QFX=-RLOW*HLFLX*WETM
    FLX_LH=XLV*QFX
    FLHC=RLOW*CP*AKHS
    FLQC=RLOW*AKHS*WETM
!!!   QGH=PQ0/PSHLTR*EXP(A2S*(TSK-A3S)/(TSK-A4S))
    QGH=((1.0_r8-SEAMASK)*PQ0+SEAMASK*PQ0SEA)                            &
         &     /PLOW*EXP(A2S*(TLOW-A3S)/(TLOW-A4S))
    QGH=QGH/(1.0_r8-QGH)    !Convert to mixing ratio
    CPM=CP*(1.0_r8+0.8_r8*QLOW)
    !
    !***  DO NOT COMPUTE QS OVER LAND POINTS HERE SINCE IT IS
    !***  A PROGNOSTIC VARIABLE THERE.  IT IS OKAY TO USE IT 
    !***  AS A DIAGNOSTIC OVER WATER SINCE IT WILL CAUSE NO
    !***  INTERFERENCE BEFORE BEING RECOMPUTED IN MYJPBL.
    !
    IF(SEAMASK>0.5_r8)THEN
       QS=QLOW + QFX /(RLOW*AKHS)
       QS=QS/(1.0_r8-QS)
    ENDIF
    !----------------------------------------------------------------------
    !
  END SUBROUTINE SFCDIF
END MODULE Sfc_MellorYamada1
