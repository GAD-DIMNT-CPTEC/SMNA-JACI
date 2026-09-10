MODULE Cu_RAS3PHASE
  USE Utils, ONLY: &
       IJtoIBJB ,&
       NearestIJtoIBJB,&
       linearijtoibjb
  
  USE Options, ONLY: &
       reducedGrid,tbase,sthick
  
  USE Parallelism, ONLY: &
     MsgOne, &
     FatalError

  USE PhysicalFunctions,Only : fpvs2es5

  USE Mod_GET_PRS, ONLY: GET_PRS,GET_PHI,sig2press

  IMPLICIT NONE
SAVE

  PRIVATE
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(13,60) ! the '60' maps to 64-bit real
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers

  !  --- ...  Geophysics/Astronomy constants

  REAL(kind=r8),PARAMETER:: con_g      =9.80665e+0_r8     ! gravity           (m/s2)
  REAL(kind=r8),PARAMETER:: con_pi     =3.1415926535897931 ! pi
  REAL(kind=r8),PARAMETER:: con_rerth  =6.3712e+6      ! radius of earth   (m)
  REAL(kind=r8),PARAMETER::  ki=1
  !  --- ...  Thermodynamics constants


  REAL(kind=r8),PARAMETER:: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_rv     =4.6150e+2_r8      ! gas constant H2O    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_hvap   =2.5000e+6_r8      ! lat heat H2O cond   (J/kg)
  REAL(kind=r8),PARAMETER:: con_hfus   =3.3358e+5_r8      ! lat heat H2O fusion (J/kg)
  REAL(kind=r8),PARAMETER:: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cvap   =1.8460e+3_r8      ! spec heat H2O gas   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cliq   =4.1855e+3_r8      ! spec heat H2O liq   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_csol   =2.1060e+3_r8      ! spec heat H2O ice   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_ttp    =2.7316e+2_r8      ! temp at H2O 3pt     (K)
  REAL(kind=r8),PARAMETER:: con_psat   =6.1078e+2_r8      ! pres at H2O 3pt     (Pa)  
  REAL (KIND=r8), PARAMETER   :: rmwmd  =                0.622e0_r8! fracao molar entre a agua e o ar 

  REAL (KIND=r8), PARAMETER   :: wgt  = 1.0e0_r8! 

  !  Secondary constants

  REAL(kind=r8),PARAMETER:: con_rocp   =con_rd/con_cp
  REAL(kind=r8),PARAMETER:: con_fvirt  =con_rv/con_rd-1.0_r8
  REAL(kind=r8),PARAMETER:: con_eps    =con_rd/con_rv
  REAL(kind=r8),PARAMETER:: con_epsm1  =con_rd/con_rv-1.0_r8

  !     module module_ras

  !     real, parameter :: con_FVirt=0.0
  !
  !INTEGER, PARAMETER :: nrcmax= 32 ! Maximum # of random clouds per 1200s
  !INTEGER, PARAMETER :: nrcmax= 32 ! Maximum # of random clouds per 1200s
  !INTEGER, PARAMETER :: nrcmax=12 ! Maximum # of random clouds per 1200s
  INTEGER, PARAMETER :: nrcmax=15 ! Maximum # of random clouds per 1200s
  !     integer, parameter :: nrcmax=15 ! Maximum # of random clouds per 1200s
  !     integer, parameter :: nrcmax=20
  INTEGER     ::  nrcm!     nrcm     - integer, number of random clouds                  1    !
  REAL (kind=r8), ALLOCATABLE :: xkt2(:,:,:)
  real (kind=r8)  :: delt_c=1800.0_r8/3600.0_r8
  !REAL (kind=r8), PARAMETER :: delt_c=1800.0_r8
  !     Adjustment time scales in hrs for deep and shallow clouds
!    &,                                   adjts_d=3.0, adjts_s=0.5
!    &,                                   adjts_d=2.5, adjts_s=0.5
!PK  real (kind=r8), parameter :: adjts_d=2.0, adjts_s=0.5
  real (kind=r8), parameter :: adjts_d=1.0_r8      ,adjts_s=0.5_r8
  real (kind=r8), parameter :: adjts_d_ocean=9.0_r8,adjts_s_ocean=0.5_r8

  LOGICAL       , PARAMETER :: fix_ncld_hr=.TRUE.
  !
  REAL(kind=r8), PARAMETER :: ZERO=0.0_r8
  REAL(kind=r8), PARAMETER :: HALF=0.5_r8
  REAL(kind=r8), PARAMETER :: ONE=1.0_r8
  REAL(kind=r8), PARAMETER :: TWO=2.0_r8
  REAL(kind=r8), PARAMETER :: FOUR_P2=4.E2_r8
  REAL(kind=r8), PARAMETER :: FOUR=4.0_r8
  REAL(kind=r8), PARAMETER :: ONE_M1=1.E-1_r8
  REAL(kind=r8), PARAMETER :: ONE_M2=1.E-2_r8
  REAL(kind=r8), PARAMETER :: ONE_M5=1.E-5_r8
  REAL(kind=r8), PARAMETER :: ONE_M6=1.E-6_r8
  REAL(kind=r8), PARAMETER :: ONE_M10=1.E-10_r8
  !
  REAL(kind=r8), PARAMETER :: cmb2pa = 100.0_r8  ! Conversion from MB to PA
  REAL(kind=r8), PARAMETER :: onebg       = ONE / con_g
  REAL(kind=r8), PARAMETER :: gravcon     = cmb2pa * ONEBG
  REAL(kind=r8), PARAMETER :: gravfac     = con_g / CMB2PA
  REAL(kind=r8), PARAMETER :: elocp       = con_hvap  / con_cp
  REAL(kind=r8), PARAMETER :: elfocp      = (con_hvap +con_hfus) / con_cp
  REAL(kind=r8), PARAMETER :: rkapi       = ONE / con_rocp
  REAL(kind=r8), PARAMETER :: rkpp1i      = ONE / (ONE+con_rocp)
  REAL(kind=r8), PARAMETER :: zfac        = 0.28888889E-4_r8 * ONEBG
  !
  !
  !     logical, parameter :: advcld=.false. advups=.true.
  !     logical, parameter :: advcld=.true., advups=.true., advtvd=.false.
  LOGICAL, PARAMETER :: advcld=.FALSE.
  LOGICAL, PARAMETER :: advups=.FALSE.
  LOGICAL, PARAMETER :: advtvd=.TRUE.
  LOGICAL, PARAMETER :: flipv = .TRUE.!     flipv    - logical, flag for vertical direction              1    !

  !     logical, parameter :: advcld=.false., advups=.false.
  !
  REAL(kind=r8), ALLOCATABLE  ::  RASAL(:)
  REAL(kind=r8) :: AFC
  !REAL(kind=r8) :: facdt

  !     PARAMETER (DD_DP=1000.0, RKNOB=1.0, EKNOB=1.0)   ! No downdraft!
  !     PARAMETER (DD_DP=100.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=200.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=250.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=300.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=450.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=500.0,  RKNOB=0.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=500.0,  RKNOB=0.70, EKNOB=1.0)
  !     PARAMETER (DD_DP=500.0,  RKNOB=0.75, EKNOB=1.0)
  REAL(kind=r8), PARAMETER :: DD_DP=500.0_r8,  RKNOB=1.0_r8, EKNOB=1.0_r8
  !     PARAMETER (DD_DP=500.0,  RKNOB=1.5, EKNOB=1.0)
!!!!! PARAMETER (DD_DP=450.0,  RKNOB=1.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=450.0,  RKNOB=2.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=450.0,  RKNOB=0.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=0.5, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=1.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=2.0, EKNOB=1.0)
  !     PARAMETER (DD_DP=350.0,  RKNOB=3.0, EKNOB=1.0)
  !
  REAL(kind=r8), PARAMETER :: RHMAX=1.0_r8     !  MAX RELATIVE HUMIDITY
  REAL(kind=r8), PARAMETER :: QUAD_LAM=1.0_r8  !  MASK FOR QUADRATIC LAMBDA
  !     PARAMETER (RHRAM=0.15)    !  PBL RELATIVE HUMIDITY RAMP
  REAL(kind=r8), PARAMETER :: RHRAM=0.05_r8    !  PBL RELATIVE HUMIDITY RAMP
  REAL(kind=r8), PARAMETER :: HCRITD=4000.0   ! Critical Moist Static Energy
  !     PARAMETER (RHRAM=0.10)    !  PBL RELATIVE HUMIDITY RAMP
  REAL(kind=r8), PARAMETER :: HCRIT=4000.0_r8  !  Critical Moist Static Energy
  REAL(kind=r8), PARAMETER :: qudfac=QUAD_LAM*half
  !     parameter (qudfac=QUAD_LAM*0.25)    ! Yogesh's
  REAL(kind=r8), PARAMETER :: shalfac=3.0

  REAL(kind=r8), PARAMETER :: testmb=0.1_r8
  REAL(kind=r8), PARAMETER :: tstmbi=one/testmb
  !
  !
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=0.00E-5, ALMAX=1.0E-1)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=0.00E-5, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=1.00E-6, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=5.00E-6, ALMIN2=2.50E-5, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=2.0E-2)
!!!   PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-2)
  !!    PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-3)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=1.00E-5, ALMAX=1.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.00E-5, ALMAX=1.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-2)
!  REAL(kind=r8), PARAMETER :: ALMIN1=0.00E-6_r8, ALMIN2=0.00E-5_r8, ALMAX=1.0E-2_r8
  REAL(kind=r8), PARAMETER :: ALMIN1=0.00E-6_r8, ALMIN2=4.00E-5_r8, ALMAX=1.0E-2_r8
  !cnt  PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=5.0E-3)
  !LL   PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=4.0E-3)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=1.00E-5, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.00E-6, ALMIN2=5.00E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.10E-4, ALMIN2=0.15E-4, ALMAX=1.0E-1)
  !     PARAMETER (ALMIN1=0.00E-4, ALMIN2=0.40E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.20E-4, ALMIN2=0.40E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.25E-4, ALMIN2=0.50E-4, ALMAX=2.0E-2)
  !     PARAMETER (ALMIN1=0.40E-4, ALMIN2=0.50E-4, ALMAX=2.0E-2)
  !
  !REAL(kind=r8), PARAMETER :: BLDMAX = 200.0_r8 !mb
  REAL(kind=r8), PARAMETER :: BLDMAX = 300.0_r8 !mb
  !
  !INTEGER  :: KBLMX

  !     PARAMETER (QI0=1.0E-4, QW0=1.0E-4)
  !     PARAMETER (QI0=0.0E-5, QW0=0.0E-0)
  REAL(kind=r8), PARAMETER :: QI0=1.0E-5_r8, QW0=1.0E-5_r8
  !     PARAMETER (QI0=1.0E-4, QW0=1.0E-5) ! 20050509
  !     PARAMETER (QI0=1.0E-5, QW0=1.0E-6)
  !     PARAMETER (QI0=0.0E-5, QW0=0.0E-5)
!!!   PARAMETER (QI0=5.0E-4, QW0=1.0E-5)
  !     PARAMETER (QI0=5.0E-4, QW0=5.0E-4)
  !     PARAMETER (QI0=2.0E-4, QW0=2.0E-5)
  !     PARAMETER (QI0=2.0E-5, QW0=2.0E-5)
  !     PARAMETER (QI0=2.0E-4, QW0=1.0E-4)
  !     PARAMETER (QI0=2.0E-4, QW0=1.0E-5)
  !     PARAMETER (QI0=1.0E-3, QW0=2.0E-5)
  !     PARAMETER (QI0=1.0E-3, QW0=7.0E-4)
  REAL(kind=r8), PARAMETER :: C00I=1.0E-3
  REAL(kind=r8), PARAMETER :: c00=2.0e-3

  !     PARAMETER (C0I=5.0E-4)
  !     PARAMETER (C0I=4.0E-4)
  REAL(kind=r8), PARAMETER :: C0I=1.0E-3_r8

  !     parameter (c0=1.0e-3)
  !     parameter (c0=1.5e-3)
  REAL(kind=r8), PARAMETER ::  c0=2.0e-3_r8
  !     parameter (c0=1.0e-3, KBLMX=10, ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
  !     parameter (c0=2.0e-3, KBLMX=10, ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
  !
  !     parameter (TF=130.16, TCR=160.16, TCRF=1.0/(TCR-TF),TCL=2.0)
  !     parameter (TF=230.16, TCR=260.16, TCRF=1.0/(TCR-TF))
  REAL(kind=r8), PARAMETER ::  TF=233.16_r8, TCR=263.16_r8, TCRF=1.0_r8/(TCR-TF),TCL=2.0_r8
  !
  !     For Tilting Angle Specification
  !
  REAL(kind=r8) :: TLBPL(7) 
  REAL(kind=r8) :: drdp(5)
  REAL(kind=r8) :: VTP
  !
  REAL(KIND=r8), PARAMETER :: PLAC(1:8)=(/100.0_r8, 200.0_r8, 300.0_r8, 400.0_r8, 500.0_r8, 600.0_r8, 700.0_r8, 800.0_r8/)
  !     DATA TLAC/ 37.0,  25.0,  17.0,  12.0,  10.0,  8.0,  6.0,  5.0/
  !     DATA TLAC/ 35.0,  24.0,  17.0,  12.0,  10.0,  8.0,  6.0,  5.0/
  !     DATA TLAC/ 35.0,  25.0,  20.0,  17.5,  15.0,  12.5,  10.0,  5.0/
  REAL(KIND=r8), PARAMETER :: TLAC(1:8)=(/ 35.0_r8,  25.0_r8,  20.0_r8,  17.5_r8,  15.0_r8,  12.5_r8,  10.0_r8,  7.5_r8/)
  !     DATA TLAC/ 37.0,  26.0,  18.0,  14.0,  10.0,  8.0,  6.0,  5.0/
  !     DATA TLAC/ 25.0,  22.5,  20.0,  17.5,  15.0,  12.5,  10.0,  10.0/
  REAL(KIND=r8), PARAMETER :: REFP(1:6)=(/500.0_r8, 300.0_r8, 250.0_r8, 200.0_r8, 150.0_r8, 100.0_r8/)
  !     DATA REFR/ 0.25,   0.5,  0.75,   1.0,   1.5,   2.0/
  !     DATA REFR/ 0.5,   1.0,  1.5,   2.0,   3.0,   4.0/
  REAL(KIND=r8), PARAMETER :: REFR(1:6)=(/ 1.0_r8,   2.0_r8,  3.0_r8,   4.0_r8,   6.0_r8,   8.0_r8/)
  !
  REAL(kind=r8) :: AC(16)
  REAL(kind=r8) :: AD(16)
  !
  INTEGER, PARAMETER :: NQRP=500001

  REAL(kind=r8) C1XQRP, C2XQRP, TBQRP(NQRP), TBQRA(NQRP)     &
       &,                    TBQRB(NQRP)
  !
  REAL(kind=r8) :: rasalf

  INTEGER, PARAMETER :: NVTP=10001
  REAL(kind=r8)      :: C1XVTP, C2XVTP, TBVTP(NVTP)
  ! 
  !  END    module module_ras
  !
  !
  !      module module_rascnv
  !
  LOGICAL      , PARAMETER ::  vsmooth=.TRUE.!PK.false.

  REAL(KIND=r8), PARAMETER :: frac=0.5_r8, crtmsf=0.0_r8
  !     PARAMETER (MAX_NEG_BOUY=0.25, REVAP=.TRUE., CUMFRC=.true.)
  !     PARAMETER (MAX_NEG_BOUY=0.20, REVAP=.TRUE., CUMFRC=.true.)
  !     PARAMETER (MAX_NEG_BOUY=0.15, REVAP=.true., CUMFRC=.true.)
  !LL3  PARAMETER (MAX_NEG_BOUY=0.10, REVAP=.TRUE., CUMFRC=.true.)
  !     PARAMETER (MAX_NEG_BOUY=0.15, REVAP=.true., CUMFRC=.false.)
  REAL(KIND=r8), PARAMETER :: MAX_NEG_BOUY=0.15_r8
  LOGICAL      , PARAMETER :: REVAP=.TRUE., CUMFRC=.TRUE.
  !     PARAMETER (MAX_NEG_BOUY=0.15, REVAP=.true., CUMFRC=.false.)
  !     PARAMETER (MAX_NEG_BOUY=0.05, REVAP=.true., CUMFRC=.true.)


  LOGICAL, PARAMETER :: WRKFUN = .FALSE.
  !LOGICAL, PARAMETER :: CRTFUN = .TRUE.,   CALKBL = .TRUE., BOTOP=.TRUE.
   LOGICAL, PARAMETER :: CRTFUN = .TRUE.,    BOTOP=.TRUE.
  !
  !     parameter (rhfacs=0.70, rhfacl=0.70)
  !     parameter (rhfacs=0.75, rhfacl=0.75)
  !REAL(KIND=r8), PARAMETER ::  rhfacs=0.80_r8, rhfacl=0.80_r8
  REAL(KIND=r8), PARAMETER ::  rhfacs=0.70_r8, rhfacl=0.70_r8
  !     parameter (rhfacs=0.80, rhfacl=0.85)
  REAL(KIND=r8), PARAMETER :: FACE=5.0_r8, DELX=10000.0_r8, DDFAC=FACE*DELX*0.001_r8
  !
  !     real (kind=r8), parameter :: pgftop=0.7, pgfbot=0.3        &
  !     real (kind=r8), parameter :: pgftop=0.75, pgfbot=0.35      &
  !REAL (kind=r8), PARAMETER :: pgftop=0.0_r8, pgfbot=0.0_r8      
  REAL (kind=r8), PARAMETER :: pgftop=0.80_r8, pgfbot=0.30_r8      
  REAL (kind=r8), PARAMETER :: pgfgrad=(pgfbot-pgftop)*0.001_r8
  !
  !
  INTEGER       , ALLOCATABLE :: nlons(:)
  !REAL (kind=r8), ALLOCATABLE ::si(:)
  !REAL (kind=r8), ALLOCATABLE ::sl(:)
  INTEGER      :: latg

  !      end module module_rascnv
  !
!  Parameters
        INTEGER,PARAMETER:: n=624
        INTEGER,PARAMETER:: m=397
        INTEGER(KIND=i8),PARAMETER:: mata=-1727483681_i8 ! constant vector a
        INTEGER(KIND=i8),PARAMETER:: umask=-2147483648_i8 ! most significant w-r bits
        INTEGER(KIND=i8),PARAMETER:: lmask =2147483647_i8 ! least significant r bits
        INTEGER(KIND=i8),PARAMETER:: tmaskb=-1658038656_i8 ! tempering parameter
        INTEGER(KIND=i8),PARAMETER:: tmaskc=-272236544_i8 ! tempering parameter
        INTEGER(KIND=i8),PARAMETER:: mag01(0:1)=(/0_i8,mata/)
        INTEGER,PARAMETER:: iseed=4357

!  Defined types
        TYPE random_stat
          PRIVATE
          INTEGER:: mti=n+1
          INTEGER:: mt(0:n-1)
          INTEGER:: iset
          REAL(KIND=r8)   :: gset
        END TYPE
!  Saved data
        TYPE(random_stat),SAVE:: sstat
!  Overloaded interfaces
!        INTERFACE random_setseed
!          MODULE PROCEDURE random_setseed_s
!          MODULE PROCEDURE random_setseed_t
!        END INTERFACE
!        INTERFACE random_number_ras
!          MODULE PROCEDURE random_number_i
!          MODULE PROCEDURE random_number_s
!          MODULE PROCEDURE random_number_t
!        END INTERFACE


  PUBLIC :: Init_Cu_RAS3PHASE
  PUBLIC :: Run_Cu_RAS3PHASE
  PUBLIC :: Finalize_Cu_RAS3PHASE
  
  
CONTAINS


!  Subprogram random_setseed_s
!  Sets seed in saved mode.
        SUBROUTINE random_setseed_s(inseed)
          IMPLICIT NONE
          INTEGER,INTENT(in):: inseed
          INTEGER :: rec
          CALL random_setseed_t(inseed,rec)
        END SUBROUTINE random_setseed_s
!  Subprogram random_setseed_t
!  Sets seed in thread-safe mode.
        SUBROUTINE random_setseed_t(inseed,rec)
          IMPLICIT NONE
          INTEGER,INTENT(in):: inseed
          INTEGER           :: rec
          INTEGER ii,mti
          ii=inseed
          IF(ii.EQ.0) ii=iseed
          sstat%mti=n
          sstat%mt(0)=IAND(ii,-1)
          DO mti=1,n-1
            sstat%mt(mti)=IAND(69069*sstat%mt(mti-1),-1)
          ENDDO
          sstat%iset=0
          sstat%gset=0.
        END SUBROUTINE random_setseed_t
!  Subprogram random_number_i
!  Generates random numbers in interactive mode.
        SUBROUTINE random_number_i(harvest,inseed)
          IMPLICIT NONE
          REAL(KIND=r8),INTENT(out):: harvest(:)
          INTEGER,INTENT(in):: inseed
          INTEGER :: rec
         ! TYPE(random_stat) stat
          CALL random_setseed_t(inseed,rec)
          CALL random_number_t(harvest,rec)
        END SUBROUTINE random_number_i
!  Generates random numbers in saved mode; overloads Fortran 90 standard.
        SUBROUTINE random_number_s(harvest)
          IMPLICIT NONE
          REAL(KIND=r8),INTENT(out):: harvest(:)
          INTEGER :: rec
          IF(sstat%mti.EQ.n+1) CALL random_setseed_t(iseed,rec)
          CALL random_number_t(harvest,rec)
        END SUBROUTINE random_number_s
!  Subprogram random_number_t
!  Generates random numbers in thread-safe mode.
        SUBROUTINE random_number_t(harvest,stat)
          IMPLICIT NONE
          REAL(KIND=r8), PARAMETER  :: twop32=2.0_r8**32
          REAL(KIND=r8), PARAMETER  :: twop32m1i=1.0_r8/(twop32-1.0_r8)
          REAL(KIND=r8),INTENT(out):: harvest(:)
          INTEGER,INTENT(inout):: stat
          INTEGER j,kk,y
          INTEGER tshftu,tshfts,tshftt,tshftl
          y=0
          tshftu=ISHFT(y,-11)
          tshfts=ISHFT(y,7)
          tshftt=ISHFT(y,15)
          tshftl=ISHFT(y,-18)
!          DO j=1,SIZE(harvest)
!             IF(sstat%mti.GE.n) THEN
!                DO kk=0,n-m-1
!                   y=INT(IOR(IAND(sstat%mt(kk),umask),IAND(sstat%mt(kk+1),lmask)),KIND=i4)
!                   sstat%mt(kk)=INT(IEOR(IEOR(sstat%mt(kk+m),ISHFT(y,-1)), &
!                        &                           mag01(IAND(y,1))),kind=i4)
!                ENDDO
!                DO kk=n-m,n-2
!                   y=INT(IOR(IAND(sstat%mt(kk),umask),IAND(sstat%mt(kk+1),lmask)),KIND=i4)
!                   sstat%mt(kk)=INT(IEOR(IEOR(sstat%mt(kk+(m-n)),ISHFT(y,-1)), &
!                        &                           mag01(IAND(y,1))),KIND=i4)
!                ENDDO
!                y=INT(IOR(IAND(sstat%mt(n-1),umask),IAND(sstat%mt(0),lmask)),KIND=i4)
!                sstat%mt(n-1)=INT(IEOR(IEOR(sstat%mt(m-1),ISHFT(y,-1)), &
!                     &                      mag01(IAND(y,1))),KIND=i4)
!                sstat%mti=0
!             ENDIF
!             y=sstat%mt(sstat%mti)
!             y=IEOR(y,tshftu)
!             y=INT(IEOR(y,IAND(tshfts,tmaskb)),KIND=i4)
!             y=INT(IEOR(y,IAND(tshftt,tmaskc)),KIND=i4)
!             y=IEOR(y,tshftl)
!             IF(y.LT.0) THEN
!                harvest(j)=(REAL(y,KIND=r8)+twop32)*twop32m1i
!             ELSE
!                harvest(j)=REAL(y,KIND=r8)*twop32m1i
!             ENDIF
!             sstat%mti=sstat%mti+1
!          ENDDO
        END SUBROUTINE random_number_t

  !-----------------------------------------------------------------------------------------
  SUBROUTINE Init_Cu_RAS3PHASE(kMax,dt,fhour,idate,iMax,jMax,ibMax,jbMax)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(kind=r8), INTENT(IN   ) :: dt
    REAL(kind=r8), INTENT(IN   ) ::  fhour
    INTEGER      , INTENT(IN   ) :: idate(4)!=(/00,01,29,2015/)
    INTEGER      , INTENT(IN   ) :: iMax,jMax,ibMax,jbMax
    !REAL(kind=r8), INTENT(IN   ) ::  si_in(kMax+1)
    !REAL(kind=r8), INTENT(IN   ) ::  sl_in(kMax)
    !real (kind=r8) :: rannum(lonr*latr*nrcm)

    INTEGER       :: nrc,i,j,ij,ib,jb
    INTEGER       :: iseed
    INTEGER       :: seed0
    REAL(kind=r8),ALLOCATABLE ::  brf(:,:)
    REAL(kind=r8) ::  wrk(1)
    REAL(kind=r8), PARAMETER :: cons_0=0.0_r8,   cons_24=24.0_r8
    REAL(kind=r8), PARAMETER :: cons_99=99.0_r8, cons_1p0d9=1.0E9_r8
    REAL (kind=r8), ALLOCATABLE :: rannum(:)

    PRINT*,'Init_Cu_RAS3PHASE'
    !ALLOCATE (si(kMax+1));si=si_in
    !ALLOCATE (sl(kMax));sl=sl_in

    delt_c=dt/3600.0_r8
    nrcm = MAX(nrcmax, NINT((nrcmax*dt)/600.0_r8))

    ALLOCATE(rannum(iMax*jMax));rannum=0.0_r8
    ALLOCATE(brf(iMax,jMax));brf=0.0_r8


    ALLOCATE(xkt2  (ibMax,jbMax,nrcm));xkt2=0.0_r8
    ALLOCATE(nlons(ibMax));nlons=iMax
    !
    latg=jMax

    seed0 = idate(1) + idate(2) + idate(3) + idate(4)
    CALL random_setseed_s(seed0)
    wrk=0.0_r8
    CALL random_number_s(wrk)

    DO nrc=1,nrcm
       seed0 = seed0 + NINT(wrk(1)*1000.0_r8)
       iseed = INT( MOD(100.0_r8*SQRT(fhour*3600.0_r8),cons_1p0d9) + 1 + seed0 ,KIND=i4)
       CALL random_setseed_s(iseed)
       CALL random_number_s(rannum)
       ij=0
       DO j=1,jmax
          DO i=1,iMax
             ij=ij+1
              brf(i,j) = rannum(ij)
          END DO
       END DO

       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(brf(:,:),xkt2(:,:,nrc))
       ELSE
          CALL IJtoIBJB      (brf(:,:) ,xkt2(:,:,nrc) )
       END IF
    END DO


    
    DEALLOCATE(rannum)
    DEALLOCATE(brf)
    CALL ras_init(kMax)

  END SUBROUTINE Init_Cu_RAS3PHASE
  !
  SUBROUTINE set_ras_afc(dt)
    IMPLICIT NONE
    REAL(kind=r8) :: DT
    !     AFC = -(1.04E-4*DT)*(3600./DT)**0.578
    AFC = -(1.01097E-4_r8*DT)*(3600.0_r8/DT)**0.57777778_r8
  END SUBROUTINE set_ras_afc

  SUBROUTINE ras_init(levs)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: levs

    !
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: actop
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(KIND=r8) :: A(15)
    INTEGER       :: i, l
    !     PARAMETER (ACTP=1.7,   FACM=1.20)
    REAL(kind=r8), PARAMETER :: ACTP=1.7_r8
    REAL(kind=r8), PARAMETER :: FACM=1.00_r8
    !     PARAMETER (ACTP=1.7,   FACM=0.90)
    !     PARAMETER (ACTP=1.7,   FACM=0.75)
    !     PARAMETER (ACTP=1.7,   FACM=0.60)
    !     PARAMETER (ACTP=1.7,   FACM=0.5)   ! cnt
    !     PARAMETER (ACTP=1.7,   FACM=0.4)
    !     PARAMETER (ACTP=1.7,   FACM=0.0)
    !
    !
    REAL(kind=r8), PARAMETER :: PH(1:15)=(/150.0_r8, 200.0_r8, 250.0_r8, 300.0_r8, 350.0_r8, 400.0_r8, 450.0_r8, 500.0_r8,    &
         550.0_r8, 600.0_r8, 650.0_r8, 700.0_r8, 750.0_r8, 800.0_r8, 850.0_r8/)
    !
    REAL(kind=r8), PARAMETER ::  AUX(1:15)=(/ 1.6851_r8, 1.1686_r8, 0.7663_r8, 0.5255_r8, 0.4100_r8, 0.3677_r8,           &
         &       0.3151_r8, 0.2216_r8, 0.1521_r8, 0.1082_r8, 0.0750_r8, 0.0664_r8,            &
         &       0.0553_r8, 0.0445_r8, 0.0633_r8/)
    !
    !LOGICAL :: first=.TRUE.
    !
    !IF (first) THEN
    !
    ALLOCATE (rasal(levs))
    !                                   set critical workfunction arrays
    ACTOP = ACTP*FACM
    DO L=1,15
       A(L) = AUX(L)*FACM
    ENDDO
    DO L=2,15
       TEM   = 1.0_r8 / (PH(L) - PH(L-1))
       AC(L) = (PH(L)*A(L-1) - PH(L-1)*A(L)) * TEM
       AD(L) = (A(L) - A(L-1)) * TEM
    ENDDO
    AC(1)  = ACTOP
    AC(16) = A(15)
    AD(1)  = 0.0_r8
    AD(16) = 0.0_r8
    !
    !       CALL SETES
    CALL SETQRP()
    CALL SETVTP()
    !
    !kblmx = levs / 2
    !
    !       RASALF  = 0.10
    !       RASALF  = 0.20
    RASALF  = 0.30_r8
    !       RASALF  = 0.35
    !
    DO L=1,LEVS
       RASAL(L) = RASALF
    ENDDO
    !
    !
    DO i=1,7
       tlbpl(i) = (tlac(i)-tlac(i+1)) / (plac(i)-plac(i+1))
    ENDDO
    DO i=1,5
       drdp(i)  = (REFR(i+1)-REFR(i)) / (REFP(i+1)-REFP(i))
    ENDDO
    !
    VTP    = 36.34_r8*SQRT(1.2_r8)* (0.001_r8)**0.1364_r8
    !
    !IF (me .EQ. 0) PRINT *,' NO DOWNDRAFT FOR CLOUD TYPES'          &
    !     &,        ' DETRAINING WITHIN THE BOTTOM ',DD_DP,' hPa LAYERS'
    !
    !first = .FALSE.
    !ENDIF
    !
  END SUBROUTINE ras_init
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  SUBROUTINE Run_Cu_RAS3PHASE(microphys,nClass,nAeros,nCols,latco,ntrac, kMax,dt,kdt,terr,&
                              t2,t3,q2,q3,ql2,ql3,qi2,qi3,u2,v2,prsi_i,prsl_i,phii_i,phil_i,&
                              dtdt ,dqdt ,dqldt,dqidt ,dudt,dvdt,colrad,ustar,pblh,tke,mask2,&!,KPBL,ustar,&
                              kbot ,ktop ,kuo  ,raincv,gvarm,gvarp,daerdt)
    IMPLICIT NONE
    LOGICAL      , INTENT(IN   ) :: microphys
    INTEGER      , INTENT(IN   ) :: nClass
    INTEGER      , INTENT(IN   ) :: nAeros
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: latco
    INTEGER      , INTENT(IN   ) :: ntrac
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: dt
    INTEGER      , INTENT(IN   ) :: kdt
    REAL(KIND=r8), INTENT(IN   ) :: terr(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: t2 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: t3 (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: q2 (nCols,kMax)!     tgrs     - real, layer mean specific humidty ( kg/kg )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: q3 (nCols,kMax)!     tgrs     - real, layer mean specific humidty ( kg/kg )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: ql2 (nCols,kMax)!     tgrs     - real, layer mean cloud liquid water ( kg/kg  )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: ql3 (nCols,kMax)!     tgrs     - real, layer mean cloud liquid water ( kg/kg  )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: qi2 (nCols,kMax)!     tgrs     - real, layer mean cloud ice water  ( kg/kg  )      ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: qi3 (nCols,kMax)!     tgrs     - real, layer mean cloud ice water ( kg/kg  )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: u2 (nCols,kMax)!     tgrs     - real, layer mean zonal wind ( m/s )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: v2 (nCols,kMax)!     tgrs     - real, layer mean meridional wind ( m/s )      ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: phii_i (1:nCols,1:kMax+1) 
    REAL(KIND=r8), INTENT(IN   ) :: phil_i (1:nCols,1:kMax  )   
    REAL(KIND=r8), INTENT(IN   ) :: prsi_i (1:nCols,1:kMax+1)   !interface level pressure Pa
    REAL(KIND=r8), INTENT(IN   ) :: prsl_i (1:nCols,1:kMax  )   !mean  level pressure Pa
    REAL(KINd=r8), INTENT(OUT  ) :: dtdt (nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dqdt (nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dqldt(nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dqidt(nCols,kMax)
    REAL(KINd=r8), INTENT(OUT  ) :: dudt(nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dvdt(nCols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: colrad (nCols)     ! colatitudes 
    REAL(KIND=r8), INTENT(IN   ) :: ustar(nCols)        ! scale velocity
    REAL(KIND=r8), INTENT(IN   ) :: pblh(nCols)      ! pbl higth  
    REAL(KIND=r8), INTENT(IN   ) :: tke   (1:nCols,1:kMax) ! cinetic energy
    INTEGER      , INTENT(IN   ) :: mask2(nCols)! sea -land mask 
    REAL(KIND=r8), OPTIONAL,   INTENT(INOUT) :: gvarm (nCols,kmax,nClass+nAeros) ! scalar
    REAL(KIND=r8), OPTIONAL,   INTENT(INOUT) :: gvarp (nCols,kmax,nClass+nAeros) ! scalar
    REAL(KIND=r8), OPTIONAL,   INTENT(OUT  ) :: daerdt(nCols,kMax,nClass+nAeros) 

    
!    INTEGER      , INTENT(IN   ) ::  KPBL(nCols)
!    REAL(KIND=r8), INTENT(IN   ) :: ustar(nCols)
    INTEGER      , INTENT(OUT  ) ::  kbot (nCols)
    INTEGER      , INTENT(OUT  ) ::  ktop (nCols)
    INTEGER      , INTENT(OUT  ) ::  kuo  (nCols)
    REAL(kind=r8) , INTENT(OUT  ) :: raincv(nCols)
    INTEGER       :: KPBL(nCols)



    REAL(kind=r8) :: rain1(nCols)

    REAL(KIND=r8) :: tgrs (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )ix,levs !
    REAL(KIND=r8) :: qgrs (nCols,kMax,3)! qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
    REAL(KIND=r8) :: qliq (nCols,kMax)!     qgrs     - real, layer mean tracer concentrationix,levs,ntrac!
    REAL(KIND=r8) :: qice (nCols,kMax)!     qgrs     - real, layer mean tracer concentrationix,levs,ntrac!

    REAL(KIND=r8) :: ugrs (nCols,kMax)!     ugrs,vgrs- real, u/v component of layer wind ix,levs !
    REAL(KIND=r8) :: vgrs (nCols,kMax)
!    REAL(KIND=r8) :: ps(nCols)        !cb

    REAL (kind=r8), PARAMETER ::                               &
!    & dxmax=-8.8818363,   dxmin=-5.2574954,   dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-9.5468126,   dxmin=-5.2574954,   dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-18.40047804, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-15.8949521,  dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-15.31958795, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-14.95494484, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-14.50865774, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
     & dxmax=-16.118095651,dxmin=-9.800790154, dxinv=1.0/(dxmin-dxmax)
      REAL (kind=r8), PARAMETER :: rhc_max = 0.9999      ! 20060512
!     real (kind=r8), parameter :: rhc_max = 0.999       ! for pry
      REAL (kind=r8), PARAMETER :: wg = 0.5_r8 !wg - gustiness factor m/s
      REAL(kind=r8), PARAMETER :: cb2mb   = 10.0
!     ccwf     - real, multiplication factor for critical cloud         !
!                      workfunction for RAS                        2    !
      REAL (kind=r8), PARAMETER :: ccwf(2)     = (/1.0,1.0/)  !         print *,' RAS Convection scheme used with ccwf=',ccwf
!     dlqf     - real, factor for cloud condensate detrainment from     !
!                      cloud edges (RAS)                           2    !
     REAL (kind=r8), PARAMETER :: dlqf(2)     =  (/0.0,0.0/)

    INTEGER, PARAMETER ::   num_p3d=3
    !     lonf,latg- integer, number of lon/lat points                 1    !
    REAL(KIND=r8) :: prsi(nCols,kMax+1)!     prsi     - real, pressure at layer interfaces             ix,levs+1
    REAL(KIND=r8) :: prsl(nCols,kMax)!     prsl     - real, mean layer presure                       ix,levs !
    REAL(KIND=r8) :: prsik(nCols,kMax+1)
    REAL(KIND=r8) :: sik(nCols,kMax+1)
    REAL(KIND=r8) :: sikp1(nCols,kMax+1)
    REAL(KIND=r8) :: prslk(nCols,kMax)
    REAL(KIND=r8) :: DeltaP(nCols,kMax)
    REAL(KIND=r8) :: slk(nCols,kMax)
    REAL(KIND=r8) :: phii(nCols,kMax+1)
    REAL(KIND=r8) :: phil(nCols,kMax)
    REAL(KIND=r8) :: del(nCols,kMax)
    !LOCAL 
    REAL(KIND=r8) :: coslat(nCols)!     coslat   - real, cos of latitude                             nCols   !
    REAL(KIND=r8) :: work1(nCols)
    REAL(KIND=r8) :: work2(nCols)
    REAL(KIND=r8) :: rhc(nCols,kMax)
    REAL(KIND=r8) :: clw(nCols,kMax,2+nClass+nAeros)
    LOGICAL       :: CALKBLMsk(nCols)
    REAL(kind=r8) :: rannum(nCols,nrcm)
    INTEGER       :: i,k,nn
    INTEGER       :: itc
    LOGICAL       :: trans_trac       = .TRUE.    ! This is effective only for RAS
    INTEGER       :: tottracer=0
    INTEGER       :: ntcw=3
    INTEGER       :: ntoz=0
    INTEGER       :: ncld=1
    INTEGER       :: itrc
    INTEGER       :: trc_shft=0
    INTEGER       :: tracers
    REAL(kind=r8) :: crtrh(3)
    REAL(kind=r8) :: rhbbot 
    REAL(kind=r8) :: rhpbl  
    REAL(kind=r8) :: theta
    REAL(KIND=r8) :: tem
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(kind=r8) :: rhbtop 
    REAL(kind=r8) :: DDVEL(nCols)
    REAL(kind=r8) :: garea(nCols)
    REAL(kind=r8) :: ua(nCols)
    REAL(kind=r8) :: cd(nCols)
    REAL(kind=r8) :: pgrk(nCols)
    REAL(kind=r8) :: SS
    INTEGER       :: lmh(nCols) 
    REAL(kind=r8) :: ccwfac(nCols)
    REAL(kind=r8) :: dlqfac(nCols)
    REAL(kind=r8) :: ud_mf(nCols,kMax)
    REAL(kind=r8) :: dd_mf(nCols,kMax)
    REAL(kind=r8) :: det_mf(nCols,kMax)
    REAL(kind=r8) :: fscav(2+nClass+nAeros)
    REAL(kind=r8), PARAMETER :: PT01=0.01_r8 !1/100  cb
    REAL(kind=r8), PARAMETER :: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
    REAL(kind=r8), PARAMETER :: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
    REAL(kind=r8), PARAMETER :: con_rocp   =con_rd/con_cp
    REAL(kind=r8), PARAMETER :: rkap = con_rocp
    REAL(kind=r8), PARAMETER :: rk = con_rocp
    REAL(kind=r8)  ::  RKAPI 
    REAL(kind=r8)  ::  RKAPP1
    LOGICAL       :: lprnt
    !IF(nClass>0 .and. PRESENT(gvarm))THEN
       tottracer=nClass+nAeros
    !ELSE
    !   tottracer=ntrac
    !END IF
    !tottracer=ntrac
    DO i=1, nCols
      ! ps(i)=ps2(i) !cb
       lmh(i) = kMax
       DDVEL(i)=0.0_r8
       rain1(i)=0.0_r8
       prsi(i,kMax+1) =prsi_i(i,kMax+1)/1000.0_r8
       phii(i,kMax+1) =(phii_i(i,kMax+1)+terr(i))*con_g
    END DO
    DO k=1,kMax
       DO i=1, nCols
           tgrs (i,k)= t3 (i,k) ! layer mean temperature ( k )        K
           ugrs (i,k)= u2 (i,k)
           vgrs (i,k)= v2 (i,k)
           qliq (i,k)=ql3(i,k)
           qice (i,k)=qi3(i,k)
           prsi (i,k)=prsi_i(i,k)/1000.0_r8
           prsl (i,k)=prsl_i(i,k)/1000.0_r8
           phii (i,k)=(phii_i(i,k)+terr(i))*con_g
           phil (i,k)=(phil_i(i,k)+terr(i))*con_g

           dtdt (i,k)=0.0_r8
           dqdt (i,k)=0.0_r8
           dqldt(i,k)=0.0_r8
           dqidt(i,k)=0.0_r8

       END DO      
    END DO
    fscav=0.0_r8!1.0_r8/(2000.0_r8)!(5.0e-4_r8)!m
    fscav=0.0_r8!1.0_r8/(2.0_r8)!(5.0e-1_r8)!km
    DO itrc=1,3
       DO k=1,kMax 
          DO i=1, nCols
             IF(itrc ==1)THEN
                qgrs (i,k,itrc)= q3(i,k)!qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
             END IF 

             IF(itrc ==2)THEN
                qgrs (i,k,itrc)= qi3(i,k)!qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
             END IF 
             IF(itrc ==3)THEN
                qgrs (i,k,itrc)= ql3(i,k)!qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
             END IF 
          END DO      
       END DO
    END DO
    DO k=1,kMax
      DO i=1,nCols
          DeltaP(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO

    RKAPI  = 1.0_r8 / RKAP
    RKAPP1 = 1.0_r8 + RKAP
    DO k=1,kMax+1
       DO i=1,nCols 
          !sik(k)   = si(k) ** rkap
          sik(i,k)   = (prsi(i,k)/prsi(i,1)) ** rkap
          !sikp1(k) = si(k) ** rkapp1
          sikp1(i,k) = (prsi(i,k)/prsi(i,1)) ** rkapp1
       END DO
    END DO
    DO k=1,kMax
       DO i=1,nCols 
          tem        = rkapp1 * DeltaP(i,k)
          slk(i,k)   = (sikp1(i,k)-sikp1(i,k+1))/tem
          !sl(k)    = slk(k) ** rkapi
       END DO
    END DO

    DO i=1,nCols
         pgrk(i)         = (prsi(i,1)*pt01) ** rk
         prsik(i,kMax+1) = sik(i,kMax+1) * pgrk(i)
    END DO
      DO k=1,kMax
        DO i=1,nCols
           !prsi(i,k)  = si(k)*pgr(i)               ! prsi are now pressures cb
           !prsl(i,k)  = sl(k)*pgr(i)               !cb 
          prsik(i,k) = sik(i,k) * pgrk(i)
          prslk(i,k) = slk(i,k) * pgrk(i)
        END DO
      END DO



!-------------------------------------------
! calculate PBL TOP level
!-----------------------------
       ! grell mask
       !      mask2(i)=0 ! land
       !      mask2(i)=1 ! water/ocean

    DO i=1,nCols
       IF(mask2(i) == 0)THEN
          !      mask2(i)=0 ! land
          CALKBLMsk(i)=.TRUE.
       ELSE
         ! CALKBLMsk(i)=.FALSE.
          CALKBLMsk(i)=.FALSE.
       END IF
    END DO


    kpbl=1
    DO k=1,kMax
       DO i = 1, nCols
         !IF(tke(i,k) >= 0.17_r8 .and. phii(i,k)/con_g <  pblh(i))THEN   !!
          IF(tke(i,k) >= 0.028999999999_r8 .and. phii(i,k)/con_g <  (terr(i)+pblh(i)))THEN   !!
            kpbl(i)=k
          END IF
       END DO
    END DO
!---------------------------------------

    lprnt=.FALSE.
    crtrh(1)         = 0.95_r8
    crtrh(2)         = 0.90_r8    
    crtrh(3)         = 0.85_r8
    !  rhbbot = crtrh(1)
    !  rhpbl  = crtrh(2)
   !   rhbtop = crtrh(3)

    DO itc=1,nrcm
       DO i=1, nCols
          rannum(i,itc) =xkt2  (i,latco,itc)
       END DO
    END DO
      DO i = 1, nCols
         ua(i)=SQRT((ugrs(i,1)**2)+(vgrs (i,1)**2))
         SS=(ua(i) + wg)           !velocity incl. gustiness param.

         CD(i)=(ustar(i)/SS)**2

       ! colrad.....colatitude  colrad=0-3.14 from np to sp in radians
       ! colrad.....colatitude  colrad=0-pi   from np to sp in radians

       !theta = 90.0_r8   -(180.0_r8/pai)*colrad(i) ! colatitude -> latitude
       theta = (con_pi/2)-colrad(i) ! colatitude -> latitude
       coslat(i)   = ABS(cos(((colrad(i)))-(3.1415926e0_r8/2.0_r8)))

       !coslat(i) = COS(theta)
       ! the 180 degrees are divided into 37 bands with 5deg each
       ! except for the first and last, which have 2.5 deg
       ! The centers of the bands are located at:
       !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
 
        tem1       = con_rerth * (con_pi + con_pi)*coslat(i) / nlons(i)
        tem2       = con_rerth * con_pi/latg
        garea(i)   = tem1 * tem2
        !dlength(i) = sqrt( tem1*tem1 + tem2*tem2 )
      ENDDO
!
!  --- ...  figure out how many extra tracers are there
!
!      IF ( trans_trac ) THEN
!        IF ( ntcw > 0 ) THEN
!          IF ( ntoz < ntcw ) THEN
!            trc_shft = ntcw + ncld
!          ELSE
!            trc_shft = ntoz
!          ENDIF
!        ELSEIF ( ntoz > 0 ) THEN
!          trc_shft = ntoz
!        ELSE
!          trc_shft = 1
!        ENDIF
!
!        tracers   = ntrac - trc_shft
!        IF ( ntoz > 0 ) tottracer = tracers + 1  ! ozone is added separately
!      ELSE
!        tottracer = 0                            ! no convective transport of tracers
!      ENDIF

!     allocate ( clw(nCols,kMax,tottracer+2) )
      DO k = 1, kMax
        DO i = 1, nCols
          clw(i,k,1) = 0.0
          clw(i,k,2) = -999.9
        ENDDO
      ENDDO
!  --- ...  for convective tracer transport (while using ras)

!      if ( ras ) then
!        IF ( tottracer > 0 ) THEN
!
!          IF ( ntoz > 0 ) THEN
!            clw(:,:,3) = qgrs(:,:,ntoz)
!
!            IF ( tracers > 0 ) THEN
!              DO nn = 1, tracers
!                clw(:,:,3+nn) = qgrs(:,:,nn+trc_shft)
!              ENDDO
!            ENDIF
!          ELSE
!            DO nn = 1, tracers
!              clw(:,:,2+nn) = qgrs(:,:,nn+trc_shft)
!            ENDDO
!          ENDIF
!
!        ENDIF
!      endif   ! end if_ras
!  --- ...  calling precipitation processes

      DO i = 1, nCols
        work1(i) = abs((LOG(coslat(i) / (nlons(i)*latg)) - dxmin) * dxinv)
        work1(i) = MAX( 0.0_r8, MIN( 1.0_r8, work1(i) ) )
        work2(i) = 1.0_r8 - work1(i)
      ENDDO

!  --- ...  calling convective parameterization

      IF ( ntcw > 0 ) THEN
      rhbbot = crtrh(1)
      rhpbl  = crtrh(2)
      rhbtop = crtrh(3)

        DO k = 1, kMax
          DO i = 1, nCols
            rhc(i,k) = rhbbot - (rhbbot - rhbtop) * (1.0_r8 - prslk(i,k))
            rhc(i,k) = rhc_max*work1(i) + rhc(i,k)*work2(i)
            rhc(i,k) = MAX( 0.0_r8, MIN( 1.0_r8, rhc(i,k) ) )
          ENDDO
        ENDDO

        IF ( num_p3d == 3 ) THEN    ! call brad ferrier's microphysics
!  --- ...  algorithm to separate different hydrometeor species

          DO k = 1, kMax
            DO i = 1, nCols
              clw(i,k,1)  = qice(i,k)
              clw(i,k,2)  = qliq(i,k)
!  --- ...  array to track fraction of "cloud" in the form of ice
            ENDDO
          ENDDO
          IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
             DO itrc=1,nClass+nAeros
                DO k = 1, kMax
                   DO i = 1, nCols
                      clw(i,k,2+itrc)  = gvarp (i,k,itrc)
                   END DO
                END DO
            END DO
          END IF
        ELSE   ! if_num_p3d
        

          DO k = 1, kMax
            DO i = 1, nCols
              clw(i,k,1) = qgrs(i,k,1)
            ENDDO
          ENDDO

        ENDIF  ! end if_num_p3d

      ELSE    ! if_ntcw

        rhc(:,:) = 1.0

      ENDIF   ! end if_ntcw

         if (ccwf(1) >= 0.0 .or. ccwf(2) >= 0 ) then
          do i=1,nCols
            ccwfac(i) = ccwf(1)*work1(i) + ccwf(2)*work2(i)
            dlqfac(i) = dlqf(1)*work1(i) + dlqf(2)*work2(i)
          enddo
        else
          ccwfac = -999.0_r8
          dlqfac = 0.0_r8
        endif

    !    prepare input, erase output
    !
    DO i=1,nCols
       kuo (i)=0
       kbot(i)=1
       ktop(i)=1
    END DO


        call rascnv(   nCols,    kMax,   2*dt, dt, kdt,mask2,rannum,         &
                 tgrs ,    qgrs(:,:,1:1),   ugrs,    vgrs, qice,qliq,clw,tottracer ,fscav,&
                 prsi ,   prsl,   prsik,  prslk, phil,  phii,&
                 kpbl ,   cd,     rain1,  kbot,  ktop,  kuo,&
                 DDVEL, flipv, cb2mb,garea, lmh, ccwfac,dlqfac, &
                 nrcm, rhc,CALKBLMsk,  ud_mf, dd_mf, det_mf, lprnt)

    DO i=1,nCols
       kbot(i)=MAX(kbot(i),1)
       ktop(i)=MAX(ktop(i),1)
    END DO

    DO i = 1, nCols
        raincv(i) =rain1(i)*0.5_r8
       IF(RAINCV(i) > 0.0_r8)kuo(i)=1
    ENDDO

    DO k=1,kMax
       DO i=1, nCols
          IF(rain1(i) > 0.0_r8)THEN

              dtdt (i,k)=(tgrs (i,k  )-t3 (i,k))/(2*dt)

              dqdt (i,k)=(qgrs (i,k,1)-q3 (i,k))/(2*dt)

              dqldt(i,k)=(qliq (i,k  )-ql3(i,k))/(2*dt)

              dqidt(i,k)=(qice (i,k  )-qi3(i,k))/(2*dt)

              t3(i,k) = tgrs (i,k  ) !t3 (i,k) + (tgrs (i,k  ) - t2  (i,k))/(2*dt)! layer mean temperature ( k )K
              q3 (i,k)= qgrs (i,k,1)!q3 (i,k) + (qgrs (i,k,1) - q2  (i,k))/(2*dt)
              ql3(i,k)= qliq (i,k)!ql3(i,k) + (clw  (i,k,2) - qliq(i,k))/(2*dt)
              qi3(i,k)= qice (i,k)!qi3(i,k) + (clw  (i,k,1) - qice(i,k))/(2*dt)
              dudt(i,k) = wgt*((ugrs(i,k) - u2(i,k))/(2*dt))
              dvdt(i,k) = wgt*((vgrs(i,k) - v2(i,k))/(2*dt))
           END IF   
       END DO      
    END DO
          IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
!             DO itrc=nClass+1,nClass+nAeros
             DO itrc=1,nClass+nAeros

                DO k = 1, kMax
                   DO i = 1, nCols
                       daerdt(i,k,itrc) =(clw(i,k,2+itrc)-gvarp (i,k,itrc) )/(2*dt)
                       gvarp (i,k,itrc) =clw(i,k,2+itrc)
                      !PRINT*,clw(i,k,2+itrc)
                   END DO
                END DO
            END DO
          END IF

!      deallocate ( clw )

  END SUBROUTINE Run_Cu_RAS3PHASE

  !
  !  !-----------------------------------------------------------------------------------------

  SUBROUTINE rascnv(   nCols,     kMax,      dt,    dtf,  kdt,mask2,rannum      &
       &,                 tin,   qin,    uin,    vin,  qice,qliq,ccin,  trac,fscav &
       &,                 prsi,  prsl,   prsik,  prslk, phil,  phii       &
       &,                 KPBL,  CDRAG,  RAINC,  kbot,  ktop,  kuo        &
       &,                 DDVEL, FLIPV,  facmb,    garea, lmh, ccwfac,dlqfac     &
       &,                 nrcm,  rhc,  CALKBLMsk,  ud_mf, dd_mf,  det_mf,lprnt)
    !
    !*********************************************************************
    !*********************************************************************
    !************         Relaxed Arakawa-Schubert      ******************
    !************             Parameterization          ******************
    !************          Plug Compatible Driver       ******************
    !************               23 May 2002             ******************
    !************                                       ******************
    !************               Developed By            ******************
    !************                                       ******************
    !************             Shrinivas Moorthi         ******************
    !************                                       ******************
    !************                  EMC/NCEP             ******************
    !*********************************************************************
    !*********************************************************************
    !
    !
    !      use module_ras, DPD => DD_DP
    !      use module_rascnv
    IMPLICIT NONE
    !
    !
    !      input
    !
    INTEGER      , INTENT(IN   ) ::  nrcm!     nrcm     - integer, number of random clouds                  1    !
    INTEGER      , INTENT(IN   ) ::  trac!     ntrac    - integer, number of tracers                        1    !
    INTEGER      , INTENT(IN   ) :: nCols!     nCols, IX   - integer, horiz dimention and num of used pts      1    !
    INTEGER      , INTENT(IN   ) :: kMax !     levs     - integer, vertical layer dimension                 1    !
    REAL(kind=r8), INTENT(IN   ) :: DT !     dtp,dtf  - real, time interval (second)                      1    !
    REAL(kind=r8), INTENT(IN   ) :: dtf!     dtp,dtf  - real, time interval (second)                      1    !
    INTEGER      , INTENT(IN   ) :: kdt
    REAL(kind=r8), INTENT(IN   ) :: rannum(nCols,nrcm)
    INTEGER      , INTENT(IN   ) :: mask2(nCols)! sea -land mask 
    REAL(kind=r8), INTENT(INOUT) :: tin(nCols,kMax)
    REAL(kind=r8), INTENT(INOUT) :: qin(nCols,kMax)
    REAL(kind=r8), INTENT(INOUT) :: uin(nCols,kMax)
    REAL(kind=r8), INTENT(INOUT) :: vin(nCols,kMax)
    REAL(kind=r8), INTENT(INOUT) :: qice(nCols,kMax)
    REAL(kind=r8), INTENT(INOUT) :: qliq(nCols,kMax)
    REAL(kind=r8), INTENT(INOUT) :: ccin(nCols,kMax,trac+2)
    REAL(kind=r8), INTENT(IN   ) :: prsi(nCols,kMax+1) 
    REAL(kind=r8), INTENT(IN   ) :: prsik(nCols,kMax+1)
    REAL(kind=r8), INTENT(IN   ) :: prsl(nCols,kMax)
    REAL(kind=r8), INTENT(IN   ) :: prslk(nCols,kMax+1)
    REAL(kind=r8), INTENT(IN   ) :: phil(nCols,kMax)
    REAL(kind=r8), INTENT(IN   ) :: phii(nCols,kMax+1)
    INTEGER      , INTENT(IN   ) :: KPBL(nCols)
    REAL(kind=r8), INTENT(IN   ) :: CDRAG(nCols)

    REAL(kind=r8), INTENT(OUT  ) :: RAINC(nCols) 
    INTEGER      , INTENT(INOUT) ::  kbot(nCols)
    INTEGER      , INTENT(INOUT) ::  ktop(nCols)
    INTEGER      , INTENT(INOUT) ::  kuo(nCols)
    REAL(kind=r8), INTENT(INOUT) :: DDVEL(nCols)
    LOGICAL      , INTENT(IN   ) :: FLIPV
    REAL(kind=r8), INTENT(IN   ) :: facmb
    REAL(kind=r8), INTENT(IN   ) :: garea(nCols)
    INTEGER      , INTENT(IN   ) :: lmh(nCols)
    REAL(kind=r8), INTENT(IN   ) :: ccwfac(nCols)  
    REAL(kind=r8), INTENT(IN   ) :: dlqfac(nCols)  
    REAL(kind=r8), INTENT(IN   ) :: rhc(nCols,kMax)
    REAL(kind=r8), INTENT(OUT  ) :: ud_mf(nCols,kMax)
    REAL(kind=r8), INTENT(OUT  ) :: dd_mf(nCols,kMax)
    REAL(kind=r8), INTENT(OUT  ) :: det_mf(nCols,kMax)
    LOGICAL      , INTENT(IN   ) :: lprnt
    LOGICAL      , INTENT(IN   ) :: CALKBLMsk(nCols)

    real(kind=r8), intent(in) :: fscav(trac+2)



    !
    !     locals
    !
    INTEGER ::  NCRND(nCols)

    REAL(kind=r8) :: RAIN(nCols)
    REAL(kind=r8) :: toi(nCols,kMax)
    REAL(kind=r8) :: qoi(nCols,kMax)
    REAL(kind=r8) :: uvi(nCols,kMax,trac+2)   
    REAL(kind=r8) :: TCU(nCols,kMax)
    REAL(kind=r8) :: QCU(nCols,kMax)
    REAL(kind=r8) :: PCU(nCols,kMax)
    REAL(kind=r8) :: clw(nCols,kMax)
    REAL(kind=r8) :: cli(nCols,kMax)  
    REAL(kind=r8) :: QII(nCols,kMax)
    REAL(kind=r8) :: QLI(nCols,kMax)
    REAL(kind=r8) :: PRS(nCols,kMax+1)
    REAL(kind=r8) :: PSJ(nCols,kMax+1)     
    REAL(kind=r8) :: phi_l(nCols,kMax)
    REAL(kind=r8) :: phi_h(nCols,kMax+1)              
    REAL(kind=r8) :: RCU(nCols,kMax,trac+2)
    REAL(kind=r8) :: wfnc
    REAL(kind=r8) :: flx(nCols,kMax+1)
    REAL(kind=r8) :: FLXD(nCols,kMax+1)

    REAL(kind=r8) :: tla(nCols)
    !REAL(kind=r8) :: pl
    INTEGER       :: irnd,ib

    INTEGER      , PARAMETER :: ICM=100
    !REAL(KIND=r8), PARAMETER :: DAYLEN=86400.0_r8
    !REAL(KIND=r8), PARAMETER :: PFAC=1.0_r8/450.0_r8
    REAL(KIND=r8), PARAMETER :: clwmin=1.0e-10_r8
    INTEGER       :: IC(nCols,ICM)
    !
    INTEGER , PARAMETER :: ptrac=2
    REAL(kind=r8)  ::  ALFINT(nCols,kMax,trac+ptrac+4)
    !REAL(kind=r8) :: ALFINQ(kMax)
    REAL(kind=r8) :: PRSM(nCols,kMax)
    REAL(kind=r8) :: trcfac(nCols,trac+2,kMax)
    REAL(kind=r8) :: alfind(nCols,kMax)
    REAL(kind=r8) :: rhc_l(nCols,kMax)
    REAL(kind=r8) :: dtvd(nCols,2,4)
    !REAL(kind=r8) :: CFAC
    REAL(kind=r8) :: TEM
    !REAL(kind=r8) :: dpi
    REAL(kind=r8) :: ccwf(nCols)
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(KIND=r8) :: facdt
    !
    INTEGER       :: KCR(nCols)
    INTEGER       :: KFX(nCols)
    INTEGER       :: NCMX(nCols)
    INTEGER       :: NC
    INTEGER       :: KTEM(nCols)
    INTEGER       :: I
    INTEGER       :: L
    INTEGER       :: lm1
    INTEGER       :: ntrc
    !INTEGER       :: ia
    INTEGER       :: ll
    INTEGER       :: km1
    INTEGER       :: kp1
    INTEGER       :: ipt
    !INTEGER       :: lv
    INTEGER       :: KBL(nCols)
    INTEGER       :: n 
    INTEGER       :: KRMIN(nCols)
    INTEGER       :: KRMAX(nCols)
    INTEGER       :: KFMAX(nCols)
    INTEGER       :: kblmx(nCols)
    !
    LOGICAL       :: DNDRFT, lprint
      real(kind=r8) :: sgcs(nCols,kMax),C0(nCols), C0I(nCols)
!
!  Scavenging related parameters
!
      real(kind=r8) ::fscav_(trac+2)  ! Fraction scavenged per km

    !
    !     locals
    !
    ncrnd=0

    RAIN=0.0_r8
    toi=0.0_r8
    qoi=0.0_r8
    uvi=0.0_r8
    TCU=0.0_r8
    QCU=0.0_r8
    PCU=0.0_r8
    clw=0.0_r8
    cli=0.0_r8
    QII=0.0_r8
    QLI=0.0_r8
    PRS=0.0_r8
    PSJ=0.0_r8  
    phi_l=0.0_r8
    phi_h=0.0_r8
    RCU=0.0_r8
    wfnc=0.0_r8
    flx=0.0_r8
    FLXD=0.0_r8

    tla=0.0_r8
    IC=0
    !
    PRSM=0.0_r8
    trcfac=0.0_r8
    alfind=0.0_r8
    rhc_l=0.0_r8
    dtvd=0.0_r8
    TEM=0.0_r8
    ccwf=0.0_r8
    tem1=0.0_r8
    tem2=0.0_r8
    !
    KCR=0
    KFX=0
    NCMX=0
    NC=0
    KTEM=0
    I=0
    L=0
    lm1=0
    ntrc=0
    ll=0
    km1=0
    kp1=0
    ipt=0
    KBL=0
    n =0
    KRMIN=0
    KRMAX=0
    KFMAX=0
      fscav_ = 0.0                        ! By default no scavenging
      if (trac > 0) then
        do i=1,trac
          fscav_(i) = fscav(i)
        enddo
      endif

    !
    km1    = kMax - 1
    kp1    = kMax + 1
    !

    ntrc = trac
    trcfac(:,:,:) = 1.0_r8             !  For other tracers
    IF (CUMFRC) THEN
       ntrc = ntrc + ptrac
       !       trcfac(trac+1) = 0.45_r8       !  For press grad correction c=0.55_r8
       !       trcfac(trac+2) = 0.45_r8       !  in momentum mixing calculations
    ENDIF
    !
    !IF (.NOT. ALLOCATED(alfint))THEN
    !   ALLOCATE(alfint(kMax,ntrc+2+4))
       alfint=0.0_r8
    !END IF
    !
    CALL set_ras_afc(dt)
    !
    ccwf = 0.5_r8
    DO IPT=1,nCols
      tem     = 1.0 + dlqfac(ipt)
      c0(ipt)      = c00  * tem
      c0i(ipt)     = c00i * tem

       !CALKBL=CALKBLMsk(ipt)
       !
       ! Resolution dependent press grad correction momentum mixing
       !
       IF (CUMFRC) THEN
          IF (ccwfac(ipt) >= 0.0_r8) ccwf(ipt) = ccwfac(ipt)
       ENDIF
       DO l=1,kMax
          ud_mf (ipt,l)  = 0.0_r8
          dd_mf (ipt,l)  = 0.0_r8
          det_mf(ipt,l)  = 0.0_r8
       ENDDO
    END DO


!
!     Compute NCRND  : here LMH is the number of layers above the
!                      bottom surface.  For sigma coordinate LMH=K.
!                      if flipv is true, then input variables are from bottom
!                      to top while RAS goes top to bottom
!
    DO IPT=1,nCols

       IF (flipv) THEN
          ll  = kp1 - LMH(ipt)
          tem = 1.0_r8 / prsi(ipt,ll)
       ELSE
          ll  = LMH(ipt)
          tem = 1.0_r8 / prsi(ipt,ll+1)
       ENDIF
       KRMIN(ipt) = 1
       KRMAX(ipt) = km1
       KFMAX(ipt) = KRMAX(ipt)
       kblmx(ipt) = 1
       DO L=1,LMH(ipt)-1
          ll = l
          IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
          sgcs(ipt,l)  = prsl(ipt,ll) * tem
          IF (sgcs(ipt,l) .LE. 0.050) KRMIN(ipt) = L

!         IF (sgcs(ipt,l) .LE. 0.700) KRMAX(ipt) = L
!         IF (sgcs(ipt,l) .LE. 0.800) KRMAX(ipt) = L
          IF (sgcs(ipt,l) .LE. 0.760) KRMAX(ipt) = L

!         IF (sgcs(ipt,l) .LE. 0.930) KFMAX(ipt) = L
          IF(mask2(ipt) == 0 )THEN
          !      mask2(i)=0 ! land
             IF (sgcs(ipt,l) .LE. 1.0_r8) KFMAX(ipt) = L    ! Commented on 20060202
          ELSE
       !      mask2(i)=1 ! water/ocean

             IF (sgcs(ipt,l) .LE. 0.970_r8) KFMAX(ipt) = L    ! Commented on 20060202
          END IF
!         IF (sgcs(l,ipt) .LE. 0.700) kblmx(ipt) = L    ! Commented on 20101015
!         IF (sgcs(ipt,l) .LE. 0.650) kblmx(ipt) = L    ! Commented on 20060202
          IF (sgcs(ipt,l) .LE. 0.600) kblmx(ipt) = L    ! 

       ENDDO
        krmin(ipt) = max(krmin(ipt),2)
    END DO
    
    
    DO IPT=1,nCols

       !     if (lprnt .and. ipt .eq. ipr) print *,' krmin=',krmin(ipt),' krmax=',
       !    &krmax(ipt),' kfmax=',kfmax(ipt),' LMH(ipt)=',LMH(ipt),' tem=',tem
       !
       IF (fix_ncld_hr) THEN
!!!       NCRND(ipt) = min(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)) * (DTF/1200) + 0.50001
          !NCRND(ipt) = INT(min(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)) * (DTF/1200) + 0.10001,KIND=i4)
          NCRND(ipt) = INT(min(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)) * (DTF/DT) + 0.10001,KIND=i4)
!!        NCRND(ipt) = min(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)) * (DTF/600) + 0.50001
!         NCRND(ipt) = min(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)) * (DTF/360) + 0.50001
!    &                                         + 0.50001
!         NCRND(ipt) = INT(min(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)) * min(1.0,DTF/360) + 0.1,KIND=i4)
          facdt = delt_c / dt
       ELSE
          NCRND(ipt) = INT(MIN(nrcmax, (KRMAX(ipt)-KRMIN(ipt)+1)),KIND=i4)
          facdt = 1.0 / 3600.0
       ENDIF
       IF (DT .GT. DTF) NCRND(ipt) = (5*NCRND(ipt)) / 4
       NCRND(ipt)   = MAX(NCRND(ipt), 1)
       !
       KCR(ipt)    = MIN(LMH(ipt),KRMAX(ipt))
       KTEM(ipt)   = MIN(LMH(ipt),KFMAX(ipt))
       KFX (ipt)   = KTEM(ipt) - KCR(ipt)
       !     if(lprnt)print*,' enter RASCNV k=',k,' ktem=',ktem(ipt),' LMH(ipt)='
       !    &,                 LMH(ipt)
       !    &,               ' krmax=',krmax(ipt),' kfmax=',kfmax(ipt)
       !    &,               ' kcr=',kcr(ipt), ' cdrag=',cdrag(ipr)

       IF (KFX(ipt) .GT. 0) THEN
          IF (BOTOP) THEN
             DO NC=1,KFX(ipt)
                IC(ipt,NC) = KTEM(ipt) + 1 - NC
             ENDDO
          ELSE
             DO NC=KFX(ipt),1,-1
                IC(ipt,NC) = KTEM(ipt) + 1 - NC
             ENDDO
          ENDIF
       ENDIF
       !
       NCMX(ipt)  = KFX(ipt) + NCRND(ipt)
       IF (NCRND(ipt) .GT. 0) THEN
          DO I=1,NCRND(ipt)
             !IRND = INT((RANNUM(ipt,I)-0.0005_r8)*(KCR(ipt)-KRMIN(ipt)+1),KIND=i4)
             IRND = INT((RANNUM(ipt,I))*(KCR(ipt)-KRMIN(ipt)+1),KIND=i4)
             IC(ipt,KFX(ipt)+I) = IRND + KRMIN(ipt)
          ENDDO
       ENDIF
!       PRINT*,MAXVAL(IC),MINVAL(IC),NCRND(ipt),KCR(ipt),KFX(ipt),KRMIN(ipt),IRND
!       STOP
    END DO
    
    
    DO IPT=1,nCols

       !
       !     ia = 1
       !
       !     print *,' in rascnv: k=',k,'lat=',lat,' lprnt=',lprnt
       !     if (lprnt) then
       !        if (me .eq. 0) then
       !        print *,' tin',(tin(ia,l),l=k,1,-1)
       !        print *,' qin',(qin(ia,l),l=k,1,-1)
       !     endif
       !
       !
       !lprint = lprnt .AND. ipt .EQ. ipr
       lprint = lprnt
       !       kuo(ipt)  = 0
       DO l=1,kMax
          ll = l
          IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
          CLW(ipt,l)     = 0.0_r8       ! Assumes initial value of Cloud water
          CLI(ipt,l)     = 0.0_r8       ! Assumes initial value of Cloud ice
          ! to be zero i.e. no environmental condensate!!!
          !         CLT(ipt,l) = 0.0_r8
          QII(ipt,l)     = 0.0_r8
          QLI(ipt,l)     = 0.0_r8
          !                          Initialize heating, drying, cloudiness etc.
          tcu(ipt,l)     = 0.0_r8
          qcu(ipt,l)     = 0.0_r8
          pcu(ipt,l)     = 0.0_r8
          flx(ipt,l)     = 0.0_r8
          flxd(ipt,l)    = 0.0_r8
          rcu(ipt,l,1)   = 0.0_r8
          rcu(ipt,l,2)   = 0.0_r8
          !                          Transfer input prognostic data into local variable
          toi(ipt,l)     = tin(ipt,ll)
          qoi(ipt,l)     = qin(ipt,ll)
          uvi(ipt,l,trac+1) = uin(ipt,ll)
          uvi(ipt,l,trac+2) = vin(ipt,ll)
          !
          DO n=1,trac
             uvi(ipt,l,n) = ccin(ipt,ll,n+2)
          ENDDO
          !
       ENDDO
       flx(ipt,kMax+1)  = 0.0_r8
       flxd(ipt,kMax+1) = 0.0_r8
       !
       IF (ccin(ipt,1,2) .LE. -999.0_r8) THEN
          DO l=1,kMax
             ll = l
             IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
             !PK tem = ccin(ipt,ll,1)                                      &
             !PK     &            * MAX(ZERO, MIN(ONE, (TCR-toi(ipt,L))*TCRF))
             !PK ccin(ipt,ll,2) = ccin(ipt,ll,1) - tem
             !PK ccin(ipt,ll,1) = tem
             tem = qice(ipt,ll)                                      &
                  &            * MAX(ZERO, MIN(ONE, (TCR-toi(ipt,L))*TCRF))
             qliq(ipt,ll) = qice(ipt,ll) - tem
             qice(ipt,ll) = tem

          ENDDO
       ENDIF

       IF (advcld) THEN
          DO l=1,kMax
             ll = l
             IF (flipv) ll = kp1 -l ! Input variables are bottom to top!
             !PK QII(ipt,L) = ccin(ipt,ll,1)
             !PK QLI(ipt,L) = ccin(ipt,ll,2)
             QII(ipt,L) = qice(ipt,ll)
             QLI(ipt,L) = qliq(ipt,ll)

          ENDDO
       ENDIF
       !
       KBL(ipt)  = KPBL(ipt)
       ! IF (flipv) KBL(ipt)  = MAX(MIN(k, kp1-KPBL(ipt)), k/2)
       ! IF (flipv) KBL(ipt)  = MAX(MIN(kMax, kp1-KPBL(ipt)), kMax)
       IF (flipv) KBL(ipt)  = MIN(kMax, kp1-KPBL(ipt))

       rain(ipt) = 0.0_r8
       !
       DO L=1,kp1
          ll = l
          IF (flipv) ll = kp1 + 1 - l      ! Input variables are bottom to top!
          PRS(ipt,LL)   = prsi(ipt, L) * facmb ! facmb is for conversion to MB
          PSJ(ipt,LL)   = prsik(ipt,L)
          phi_h(ipt,LL) = phii(ipt,L)
       ENDDO
       !
       DO L=1,kMax
          ll = l
          IF (flipv) ll = kp1 - l          ! Input variables are bottom to top!
          PRSM(ipt,LL)  = prsl(ipt, L) * facmb ! facmb is for conversion to MB
          phi_l(ipt,LL) = phil(ipt,L)
          rhc_l(ipt,LL) = rhc(ipt,L)
       ENDDO
       !
       !     if(lprint) print *,' PRS=',PRS
       !     if(lprint) print *,' PRSM=',PRSMipt,
       !     if (lprint) then
       !        print *,' qns=',qns(ia),' qoi=',qn0(ipt,ia,kMax),'qin=',qin(ia,1)
       !        if (me .eq. 0) then
       !        print *,' toi',(tn0(ia,l),l=1,kMax)
       !        print *,' qoi',(qn0(ia,l),l=1,kMax),' kbl=',kbl(ipt)
       !     endif
       !
       !
!       do l=k,kctop(1),-1
!!        DPI(L)  = 1.0 / (PRS(L+1) - PRS(L))
!       enddo
!
       !     print *,' ipt=',ipt
    END DO
    
    
    DO IPT=1,nCols

       IF (advups) THEN               ! For first order upstream for updraft
          alfint(ipt,:,:) = 1.0_r8
       ELSEIF (advtvd) THEN           ! TVD flux limiter scheme for updraft
          alfint(ipt,:,:) = 1.0_r8
          l   = krmin(ipt)
          lm1 = l - 1
          dtvd(ipt,1,1) = con_cp*(toi(ipt,l)-toi(ipt,lm1)) + phi_l(ipt,l)-phi_l(ipt,lm1)        &
               &              + con_hvap *(qoi(ipt,l)-qoi(ipt,lm1))
          dtvd(ipt,1,2) = qoi(ipt,l) - qoi(ipt,lm1)
          dtvd(ipt,1,3) = qli(ipt,l) - qli(ipt,lm1)
          dtvd(ipt,1,4) = qii(ipt,l) - qii(ipt,lm1)
          DO l=krmin(ipt)+1,kMax
             lm1 = l - 1

             !     print *,' toi=',toi(ipt,l),toi(lm1),' phi_l=',phi_l(ipt,l),phi_l(ipt,lm1)
             !    &,' qoi=',qoi(ipt,l),qoi(ipt,lm1),' con_cp=',con_cp,' con_hvap =',con_hvap 

             dtvd(ipt,2,1)   = con_cp*(toi(ipt,l)-toi(ipt,lm1)) + phi_l(ipt,l)-phi_l(ipt,lm1)    &
                  &                  + con_hvap *(qoi(ipt,l)-qoi(ipt,lm1))
 
             !     print *,' l=',l,' dtvd=',dtvd(ipt,:,1)
 
             IF (ABS(dtvd(ipt,2,1)) > 1.0e-10_r8) THEN
                tem1        = dtvd(ipt,1,1) / dtvd(ipt,2,1)
                tem2        = ABS(tem1)
                alfint(ipt,l,1) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for h
             ENDIF

             !     print *,' alfint=',alfint(ipt,l,1),' l=',l,' ipt=',ipt

             dtvd(ipt,1,1)   = dtvd(ipt,2,1)
             !
             dtvd(ipt,2,2)   = qoi(ipt,l) - qoi(ipt,lm1)
 
             !     print *,' l=',l,' dtvd2=',dtvd(ipt,:,2)
 
             IF (ABS(dtvd(ipt,2,2)) > 1.0e-10_r8) THEN
                tem1        = dtvd(ipt,1,2) / dtvd(ipt,2,2)
                tem2        = ABS(tem1)
                alfint(ipt,l,2) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for q
             ENDIF
             dtvd(ipt,1,2)   = dtvd(ipt,2,2)
             !
             dtvd(ipt,2,3)   = qli(ipt,l) - qli(ipt,lm1)

             !     print *,' l=',l,' dtvd3=',dtvd(ipt,:,3)

             IF (ABS(dtvd(ipt,2,3)) > 1.0e-10_r8) THEN
                tem1        = dtvd(ipt,1,3) / dtvd(ipt,2,3)
                tem2        = ABS(tem1)
                alfint(ipt,l,3) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for ql
             ENDIF
             dtvd(ipt,1,3)   = dtvd(ipt,2,3)
             !
             dtvd(ipt,2,4)   = qii(ipt,l) - qii(ipt,lm1)

             !     print *,' l=',l,' dtvd4=',dtvd(ipt,:,4)

             IF (ABS(dtvd(ipt,2,4)) > 1.0e-10_r8) THEN
                tem1        = dtvd(ipt,1,4) / dtvd(ipt,2,4)
                tem2        = ABS(tem1)
                alfint(ipt,l,4) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2)   ! for qi
             ENDIF
             dtvd(ipt,1,4)   = dtvd(ipt,2,4)
          ENDDO
          !
          IF (ntrc > 0) THEN
             DO n=1,ntrc
                l = krmin(ipt)
                dtvd(ipt,1,1)   = uvi(ipt,l,n) - uvi(ipt,l-1,n)
                DO l=krmin(ipt)+1,kMax
                   dtvd(ipt,2,1)     = uvi(ipt,l,n) - uvi(ipt,l-1,n)

                   !     print *,' l=',l,' dtvdn=',dtvd(ipt,:,1),' n=',n,' l=',l

                   IF (ABS(dtvd(ipt,2,1)) > 1.0e-10_r8) THEN
                      tem1          = dtvd(ipt,1,1) / dtvd(ipt,2,1)
                      tem2          = ABS(tem1)
                      alfint(ipt,l,n+4) = 1.0_r8 - 0.5_r8*(tem1 + tem2)/(1.0_r8 + tem2) ! for tracers
                   ENDIF
                   dtvd(ipt,1,1)     = dtvd(ipt,2,1)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          alfint(ipt,:,:) = 0.5              ! For second order scheme
       ENDIF
       alfind(ipt,:)   = 0.5
    END DO
    
    
    DO IPT=1,nCols

       !
       !     print *,' after alfint for ipt=',ipt
       IF (CUMFRC) THEN

          DO l=krmin(ipt),kMax
             tem = 1.0_r8 - MAX(pgfbot, MIN(pgftop, pgftop+pgfgrad*prsm(ipt,l)))
             trcfac(ipt,trac+1,l) = tem
             trcfac(ipt,trac+2,l) = tem
          ENDDO
       ENDIF
       !
       lprint = lprnt !.AND. ipt .EQ. ipr
       !     if (lprint) then
       !       print *,' trcfac=',trcfac(1+trac,krmin(ipt):kMax)
       !       print *,' alfint=',alfint(krmin(ipt):kMax,1)
       !       print *,' alfinq=',alfint(krmin(ipt):kMax,2)
       !       print *,' alfini=',alfint(krmin(ipt):kMax,4)
       !       print *,' alfinu=',alfint(krmin(ipt):kMax,5)
       !     endif
       !
       !IF (CALKBL) kbl(ipt) = kMax
       IF(CALKBLMsk(ipt)) kbl(ipt) = kMax
    END DO
    
    
          !       lprint = lprnt .and. ipt .eq. ipr
          !       lprint = lprnt .and. ipt .eq. ipr .and. ib .eq. 41
          !
          DNDRFT = DD_DP .GT. 0.0_r8

    DO NC=1,MAXVAL(NCMX)
       DO IPT=1,nCols
          IF(NC <= NCMX(IPT) .AND. IC(ipt,NC) <= kbl(ipt))THEN
          IB = IC(ipt,NC)
          !
          !IF (ib .GT. kbl(ipt)) CYCLE
          !
          !
          WFNC = 0.0_r8
          DO L=IB,kMax+1
             FLX (ipt,L)    = 0.0_r8
             FLXD(ipt,L)   = 0.0_r8
          ENDDO

          !
          !
          !     if (me .eq. 0) then
          !     if(lprint)then
          !     print *,' CALLING CLOUD TYPE IB= ', IB,' DT=',DT,' K=',K
          !    &, 'ipt=',ipt
          !     print *,' TOI=',(TOI(ipt,L),L=IB,K)
          !     print *,' QOI=',(QOI(ipt,L),L=IB,K)
          !     endif
          !     print *,' alft=',alfint
          !
          TLA(ipt) = -10.0_r8
          END IF
       END DO
          !
          !     if (lprint) print *,' qliin=',qli
          !     if (lprint) print *,' qiiin=',qii
       DO IPT=1,nCols
          IF(NC <= NCMX(IPT) .AND. IC(ipt,NC) <= kbl(ipt))THEN
             IB = IC(ipt,NC)

          CALL CLOUD(&
                     LMH(ipt)            , &!INTEGER      , INTENT(IN     ) :: K
                     IB                  , &!INTEGER      , INTENT(IN     ) :: KD
                     ntrc                , &!INTEGER      , INTENT(IN     ) :: M     !trac +2
                     kdt                 , &!INTEGER      , INTENT(IN     ) :: kdt
                     mask2(ipt)          , &!
                     kblmx(ipt)          , &!INTEGER      , INTENT(IN     ) :: kblmx
                     FRAC                , &!REAL(kind=r8), INTENT(IN     ) :: FRACBL
                     MAX_NEG_BOUY        , &!REAL(kind=r8), INTENT(IN     ) :: MAX_NEG_BOUY
                     vsmooth             , &!LOGICAL      , INTENT(IN     ) :: vsmooth
                     facdt               , &!REAL(KIND=r8), INTENT(IN     ) :: facdt
                     ALFINT(ipt,1:kMax,1:ntrc+4), &!REAL(kind=r8), INTENT(IN     ) :: ALFINT (K,M+4)   ALFINT(kMax,trac+ptrac+4)
                     rhfacs              , &!REAL(kind=r8), INTENT(IN     ) :: RHFACS 
                     garea(ipt)          , &!REAL(kind=r8), INTENT(IN     ) :: garea
                     alfind(ipt,1:kMax)      , &!REAL(kind=r8), INTENT(IN     ) :: alfind (k)
                     rhc_l (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(IN     ) :: rhc_ls (k)
                     TOI   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: TOI    (K) 
                     QOI   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: QOI    (K )
                     UVI   (ipt,1:kMax,1:ntrc)  , &!REAL(kind=r8), INTENT(INOUT  ) :: ROI    (K,M)
                     PRS   (ipt,1:kMax+1)    , &!REAL(kind=r8), INTENT(IN     ) :: PRS    (K+1)
                     PRSM  (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(IN     ) :: PRSM   (K) 
                     phi_l (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(IN     ) :: PHIL   (K)
                     phi_h (ipt,1:kMax+1)    , &!REAL(kind=r8), INTENT(IN     ) :: PHIH   (K+1) 
                     QLI   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: QLI    (ipt,K)
                     QII   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: QII    (ipt,K)
                     KBL   (ipt)                , &!INTEGER      , INTENT(INOUT  ) :: KPBL
                     DDVEL(ipt)          , &!REAL(kind=r8), INTENT(INOUT  ) :: DSFC
                     CDRAG(ipt)          , &!REAL(kind=r8), INTENT(IN     ) :: CD
                     lprint              , &!LOGICAL      , INTENT(IN     ) :: lprnt
                     trcfac(ipt,1:ntrc,1:kMax)  , &!REAL(kind=r8), INTENT(IN     ) :: trcfac (M,k)
                     ccwf  (ipt)         , &!REAL(kind=r8), INTENT(IN     ) :: ccwf
                     TCU   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: TCU    (K)
                     QCU   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: QCU    (K)
                     RCU   (ipt,1:kMax,1:ntrc)  , &!REAL(kind=r8), INTENT(INOUT  ) :: RCU    (K,M)
                     PCU   (ipt,1:kMax)      , &!REAL(kind=r8), INTENT(INOUT  ) :: PCU    (K) 
                     FLX   (ipt,1:kMax+1)    , &!REAL(kind=r8), INTENT(INOUT  ) :: FLX    (K+1)
                     FLXD  (ipt,1:kMax+1)    , &!REAL(kind=r8), INTENT(INOUT  ) :: FLXD   (K+1)
                     RAIN (ipt)               , &!REAL(kind=r8), INTENT(INOUT  ) :: CUP
                     REVAP               , &!LOGICAL       , INTENT(IN    ) :: REVAP
                     DT                  , &!REAL(kind=r8), INTENT(IN     ) :: DT
                     WFNC                , &!REAL(kind=r8), INTENT(INOUT  ) :: WFNC
                     WRKFUN              , &!LOGICAL       , INTENT(IN    ) :: WRKFUN
                     CALKBLMsk(ipt)              , &!LOGICAL       , INTENT(IN    ) :: CALKBL
                     CRTFUN              , &!LOGICAL       , INTENT(IN    ) :: CRTFUN
                     TLA(ipt)            , &!REAL(kind=r8), INTENT(INOUT  ) :: TLA
                     DNDRFT              , &!LOGICAL       , INTENT(IN    ) ::  DNDRFT
                     DD_DP               , &!REAL(kind=r8), INTENT(IN     ) :: DPD
                     fscav_(1:ntrc)      , &!REAL(kind=r8), INTENT(IN     ) :: fscav_  (M)
                     dlqfac(ipt)         , &!REAL(kind=r8), INTENT(IN     ) :: dlq_fac
                     C0    (ipt)         , &!REAL(kind=r8), INTENT(IN     ) :: C0
                     C0I    (ipt)           )!REAL(kind=r8), INTENT(IN     ) :: C0I 

          END IF
       END DO

       DO IPT=1,nCols
          IF(NC <= NCMX(IPT) .AND. IC(ipt,NC) <= kbl(ipt))THEN
             IB = IC(ipt,NC)

             IF (lprint) THEN
                PRINT *,' after calling CLOUD TYPE IB= ', IB                      &
                     &,' rain=',rain(ipt),' prskd=',prs(ipt,ib),' qli=',qli(ipt,ib),' qii=',qii(ipt,ib)
             ENDIF
             !     if (lprint) print *,' qliou=',qli
             !     if (lprint) print *,' qiiou=',qii
             !
             DO L=IB,kMax
                ll = l
                IF (flipv) ll  = kp1 -l    ! Input variables are bottom to top!
                ud_mf(ipt,ll)  = ud_mf(ipt,ll)  + flx(ipt,l+1)
                dd_mf(ipt,ll)  = dd_mf(ipt,ll)  + flxd(ipt,l+1)
             ENDDO
             ll = ib
             IF (flipv) ll  = kp1 - ib
             det_mf(ipt,ll) = det_mf(ipt,ll) + flx(ipt,ib)
             ! 
             !     Compute cloud amounts for the Goddard radiation
             !
             !         IF (FLX(ipt,KBL(ipt)) .GT. 0.0_r8) THEN
             !           PL   = 0.5_r8 * (PRS(ipt,IB) + PRS(ipt,IB+1))
             !           CFAC = MIN(1.0_r8, MAX(0.0_r8, (850.0_r8-PL)*PFAC))
             !         ELSE
             !           CFAC = 0.0_r8
             !         ENDIF
              !
             !   Warining!!!!
              !   ------------
             !   By doing the following, CLOUD does not contain environmental
             !   condensate!
             !
             IF (.NOT. advcld) THEN
                DO l=1,kMax
                   !             clw(ipt,l ) = clw(ipt,l) + QLI(ipt,L) + QII(ipt,L)
                   clw(ipt,l ) = clw(ipt,l) + QLI(ipt,L)
                   cli(ipt,l ) = cli(ipt,l) + QII(ipt,L)
                   QLI(ipt,L)  = 0.0_r8
                   QII(ipt,L)  = 0.0_r8
                ENDDO
             ENDIF
             !
             IF (lprint)  PRINT*,IB,NC,kbl(ipt)  ,rain(ipt)
          END IF  
       ENDDO

       ENDDO                      ! End of the NC loop!



      DO IPT=1,nCols
          !
          RAINC(ipt) = rain(ipt) * 0.001_r8    ! Output rain is in meters
      ENDDO                            ! End of the IPT Loop!

       !     if(lprint)print*,' convective precip=',rain*86400/dt,' mm/day'
       !    1,               ' ipt=',ipt
       !
       !     if (lprint) then
       !        print *,' toi',(tn0(imax,l),l=1,k)
       !        print *,' qoi',(qn0(imax,l),l=1,kMax)
       !     endif
       !
       DO l=1,kMax
          ll = l
          IF (flipv) ll  = kp1 - l
             DO IPT=1,nCols
                  tin(ipt,ll)    = toi(ipt,l)                   ! Temperature
                  qin(ipt,ll)    = qoi(ipt,l)                   ! Specific humidity
                  uin(ipt,ll)    = uvi(ipt,l,trac+1)            ! U momentum
                  vin(ipt,ll)    = uvi(ipt,l,trac+2)            ! V momentum
                  !         clw(ipt,l)         = clw(ipt,l) + qli(ipt,l) + qii(ipt,l) ! Cloud condensate
                  !         ccin(ipt,ll,1) = ccin(ipt,ll,1) + clw(ipt,l)
                  DO n=1,trac
                     ccin(ipt,ll,n+2) = uvi(ipt,l,n)             ! Tracers
                  ENDDO
            ENDDO                            ! End of the IPT Loop!
       ENDDO

       IF (advcld) THEN
          DO l=1,kMax
             ll = l
             IF (flipv) ll  = kp1 - l
             !           ccin(ipt,ll,1) = qli(ipt,l) + qii(ipt,l) ! Cloud condensate
             DO IPT=1,nCols

                qice(ipt,ll) = qii(ipt,l)          ! Cloud ice
                qliq(ipt,ll) = qli(ipt,l)          ! Cloud water
             END DO
          ENDDO
       ELSE
          DO l=1,kMax
             ll = l
             IF (flipv) ll  = kp1 - l
             !           ccin(ipt,ll,1) = ccin(ipt,ll,1) + clw(ipt,l)
             DO IPT=1,nCols

                qice(ipt,ll) = qice(ipt,ll) + cli(ipt,l)
                qliq(ipt,ll) = qliq(ipt,ll) + clw(ipt,l)
             END DO
          ENDDO
       ENDIF
       !
       DO IPT=1,nCols

          kuo(ipt)  = 0
          !
          ktop(ipt) = kp1
          kbot(ipt) = kp1
       END DO

    DO IPT=1,nCols

       DO l=LMH(ipt)-1,1,-1

          IF (sgcs(ipt,l) < 0.93 .and. tcu(ipt,l) .ne. 0.0) THEN
             kuo(ipt) = 1
          ENDIF
             !  New test for convective clouds ! added in 08/21/96
             IF (clw(ipt,l)+cli(ipt,l) .GT. 0.0_r8 .OR. qli(ipt,l)+qii(ipt,l) .GT. clwmin) THEN
                 ktop(ipt) = l
             END IF 
       ENDDO
       DO l=1,km1
          IF ( tcu(ipt,l) .NE. 0.0_r8) THEN ! for r1 &rf
             IF (clw(ipt,l)+cli(ipt,l) .GT. 0.0_r8 .OR.qli(ipt,l)+qii(ipt,l) .GT. clwmin)THEN
                 kbot(ipt) = l
             END IF      
          END IF      
       ENDDO
       IF (flipv) THEN
          ktop(ipt) = kp1 - ktop(ipt)
          kbot(ipt) = kp1 - kbot(ipt)
       ENDIF
 
       !
       DDVEL(ipt) = DDVEL(ipt) * DDFAC * con_g / (prs(ipt,kMax+1)-prs(ipt,kMax))
       !
    ENDDO                            ! End of the IPT Loop!
    !
    RETURN
  END  SUBROUTINE rascnv

  !-----------------------------------------------------------------------------------------

  SUBROUTINE CRTWRK(PL, CCWF, ACR)
    !use module_ras , only : ac, ad
    IMPLICIT NONE
    !
    REAL(kind=r8), INTENT(IN   ) :: PL, CCWF
    REAL(kind=r8), INTENT(OUT  ) :: ACR
    INTEGER :: IWK
    !
    ACR=0.0_r8
    IWK = INT(PL * 0.02_r8 - 0.999999999_r8,KIND=i4)
    IWK = MAX(1, MIN(IWK,16))
    ACR = (AC(IWK) + PL * AD(IWK)) * CCWF
    !
    RETURN
  END SUBROUTINE CRTWRK


  !-----------------------------------------------------------------------------------------

  SUBROUTINE CLOUD(&
                   K           ,&!INTEGER      , INTENT(IN     ) :: K
                   KD          ,&!INTEGER      , INTENT(IN     ) :: KD
                   M           ,&!INTEGER      , INTENT(IN     ) :: M	  !trac +2
                   kdt         ,&!INTEGER      , INTENT(IN     ) :: kdt
                   mask2       ,&!I
                   kblmx       ,&!INTEGER      , INTENT(IN     ) :: kblmx
                   FRACBL      ,&!REAL(kind=r8), INTENT(IN     ) :: FRACBL
                   MAX_NEG_BOUY,&!REAL(kind=r8), INTENT(IN     ) :: MAX_NEG_BOUY
                   vsmooth     ,&!LOGICAL      , INTENT(IN     ) :: vsmooth
                   facdt       ,&!REAL(KIND=r8), INTENT(IN     ) :: facdt
                   ALFINT      ,&!REAL(kind=r8), INTENT(IN     ) :: ALFINT (K,M+4)
                   RHFACS      ,&!REAL(kind=r8), INTENT(IN     ) :: RHFACS 
                   garea       ,&!REAL(kind=r8), INTENT(IN     ) :: garea
                   alfind      ,&!REAL(kind=r8), INTENT(IN     ) :: alfind (k)
                   rhc_ls      ,&!REAL(kind=r8), INTENT(IN     ) :: rhc_ls (k)
                   TOI         ,&!REAL(kind=r8), INTENT(INOUT  ) :: TOI    (K) 
                   QOI         ,&!REAL(kind=r8), INTENT(INOUT  ) :: QOI    (K )
                   ROI         ,&!REAL(kind=r8), INTENT(INOUT  ) :: ROI    (K,M)
                   PRS         ,&!REAL(kind=r8), INTENT(IN     ) :: PRS    (K+1)
                   PRSM        ,&!REAL(kind=r8), INTENT(IN     ) :: PRSM   (K) 
                   phil        ,&!REAL(kind=r8), INTENT(IN     ) :: PHIL   (K)
                   phih        ,&!REAL(kind=r8), INTENT(IN     ) :: PHIH   (K+1) 
                   QLI         ,&!REAL(kind=r8), INTENT(INOUT  ) :: QLI    (K)
                   QII         ,&!REAL(kind=r8), INTENT(INOUT  ) :: QII    (K)
                   KPBL        ,&!INTEGER      , INTENT(INOUT  ) :: KPBL
                   DSFC        ,&!REAL(kind=r8), INTENT(INOUT  ) :: DSFC
                   CD          ,&!REAL(kind=r8), INTENT(IN     ) :: CD
                   lprnt       ,&!LOGICAL      , INTENT(IN     ) :: lprnt
                   trcfac      ,&!REAL(kind=r8), INTENT(IN     ) :: trcfac (M,k)
                   ccwf        ,&!REAL(kind=r8), INTENT(IN     ) :: ccwf
                   TCU         ,&!REAL(kind=r8), INTENT(INOUT  ) :: TCU    (K)
                   QCU         ,&!REAL(kind=r8), INTENT(INOUT  ) :: QCU    (K)
                   RCU         ,&!REAL(kind=r8), INTENT(INOUT  ) :: RCU    (K,M)
                   PCU         ,&!REAL(kind=r8), INTENT(INOUT  ) :: PCU    (K) 
                   FLX         ,&!REAL(kind=r8), INTENT(INOUT  ) :: FLX    (K+1)
                   FLXD        ,&!REAL(kind=r8), INTENT(INOUT  ) :: FLXD   (K+1)
                   CUP         ,&!REAL(kind=r8), INTENT(INOUT  ) :: CUP
                   REVAP       ,&!LOGICAL	, INTENT(IN    ) :: REVAP
                   DT          ,&!REAL(kind=r8), INTENT(IN     ) :: DT
                   WFNC        ,&!REAL(kind=r8), INTENT(INOUT  ) :: WFNC
                   WRKFUN      ,&!LOGICAL	, INTENT(IN    ) :: WRKFUN
                   CALKBL      ,&!LOGICAL	, INTENT(IN    ) :: CALKBL
                   CRTFUN      ,&!LOGICAL	, INTENT(IN    ) :: CRTFUN
                   TLA         ,&!REAL(kind=r8), INTENT(INOUT  ) :: TLA
                   DNDRFT      ,&!LOGICAL	, INTENT(IN    ) ::  DNDRFT
                   DPD         ,&!REAL(kind=r8), INTENT(IN     ) :: DPD
                   fscav_      ,&!REAL(kind=r8), INTENT(IN     ) :: fscav_  (M)
                   dlq_fac     ,&!REAL(kind=r8), INTENT(IN     ) :: dlq_fac
                   C0          ,&!REAL(kind=r8), INTENT(IN     ) :: C0
                   C0I          )!REAL(kind=r8), INTENT(IN     ) :: C0I 


    !
    !***********************************************************************
    !******************** Relaxed  Arakawa-Schubert ************************
    !****************** Plug Compatible Scalar Version *********************
    !************************ SUBROUTINE CLOUD  ****************************
    !************************  October 2004     ****************************
    !********************  VERSION 2.0  (modified) *************************
    !************* Shrinivas.Moorthi@noaa.gov (301) 763 8000(X7233) ********
    !***********************************************************************
    !*Reference:
    !-----------
    !     NOAA Technical Report NWS/NCEP 99-01:
    !     Documentation of Version 2 of Relaxed-Arakawa-Schubert
    !     Cumulus Parameterization with Convective Downdrafts, June 1999.
    !     by S. Moorthi and M. J. Suarez.
    !
    !***********************************************************************
    !
    !===>    UPDATES CLOUD TENDENCIES DUE TO A SINGLE CLOUD
    !===>    DETRAINING AT LEVEL KD.
    !
    !***********************************************************************
    !
    !===>  TOI(K)     INOUT   TEMPERATURE             KELVIN
    !===>  QOI(K)     INOUT   SPECIFIC HUMIDITY       NON-DIMENSIONAL
    !===>  ROI(K,M)   INOUT   TRACER                  ARBITRARY
    !===>  QLI(K)     INOUT   LIQUID WATER            NON-DIMENSIONAL
    !===>  QII(K)     INOUT   ICE                     NON-DIMENSIONAL

    !===>  PRS(K+1)   INPUT   PRESSURE @ EDGES        MB
    !===>  PRSM(K)    INPUT   PRESSURE @ LAYERS       MB
    !===>  PHIH(K+1)  INPUT   GEOPOTENTIAL @ EDGES  IN MKS units
    !===>  PHIL(K)    INPUT   GEOPOTENTIAL @ LAYERS IN MKS units
    !===>  PRJ(K+1)   INPUT   (P/P0)^KAPPA  @ EDGES   NON-DIMENSIONAL
    !===>  PRJM(K)    INPUT   (P/P0)^KAPPA  @ LAYERS  NON-DIMENSIONAL

    !===>  K      INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
    !===>  KD     INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
    !===>  M      INPUT   NUMBER OF TRACERS. MAY BE ZERO.
    !===>  DNDRFT INPUT   LOGICAL .TRUE. OR .FALSE.
    !===>  DPD    INPUT   Minumum Cloud Depth for DOWNDRFAT Computation hPa
    !
    !===>  TCU(K  )   UPDATE  TEMPERATURE TENDENCY       DEG
    !===>  QCU(K  )   UPDATE  WATER VAPOR TENDENCY       (G/G)
    !===>  RCU(K,M)   UPDATE  TRACER TENDENCIES          ND
    !===>  PCU(K-1)   UPDATE  PRECIP @ BASE OF LAYER     KG/M^2
    !===>  FLX(K  )   UPDATE  MASS FLUX @ TOP OF LAYER   KG/M^2
    !===>  CUP        UPDATE  PRECIPITATION AT THE SURFACE KG/M^2
    !
    !      use module_ras
    IMPLICIT NONE
    !
    !  INPUT ARGUMENTS
    INTEGER      , INTENT(IN   ) :: K
    INTEGER      , INTENT(IN   ) :: KD
    INTEGER      , INTENT(IN   ) :: M   !trac +2
    INTEGER      , INTENT(IN   ) :: kdt
    INTEGER      , INTENT(IN   ) :: mask2! sea -land mask 
    INTEGER      , INTENT(IN   ) :: kblmx
    REAL(kind=r8), INTENT(IN   ) :: FRACBL
    REAL(kind=r8), INTENT(IN   ) :: MAX_NEG_BOUY
    LOGICAL      , INTENT(IN   ) :: vsmooth
    REAL(KIND=r8), INTENT(IN   ) :: facdt
    REAL(kind=r8), INTENT(IN   ) :: ALFINT(K,M+4)
    REAL(kind=r8), INTENT(IN   ) :: RHFACS 
    REAL(kind=r8), INTENT(IN   ) :: garea
    REAL(kind=r8), INTENT(IN   ) :: alfind(k)
    REAL(kind=r8), INTENT(IN   ) :: rhc_ls(k)
    REAL(kind=r8), INTENT(INOUT) :: TOI(K) 
    REAL(kind=r8), INTENT(INOUT) :: QOI(K )
    REAL(kind=r8), INTENT(INOUT) :: ROI(K,M)
    REAL(kind=r8), INTENT(IN   ) :: PRS(K+1)
    REAL(kind=r8), INTENT(IN   ) :: PRSM(K) 
    REAL(kind=r8), INTENT(IN   ) :: PHIL(K)
    REAL(kind=r8), INTENT(IN   ) :: PHIH(K+1) 
    REAL(kind=r8), INTENT(INOUT) :: QLI(K)
    REAL(kind=r8), INTENT(INOUT) :: QII(K)
    INTEGER      , INTENT(INOUT) :: KPBL
    REAL(kind=r8), INTENT(INOUT) :: DSFC
    REAL(kind=r8), INTENT(IN   ) :: CD
    LOGICAL      , INTENT(IN   ) :: lprnt
    REAL(kind=r8), INTENT(IN   ) :: trcfac(M,k)
    REAL(kind=r8), INTENT(IN   ) :: ccwf
    REAL(kind=r8), INTENT(INOUT) :: TCU(K)
    REAL(kind=r8), INTENT(INOUT) :: QCU(K)
    REAL(kind=r8), INTENT(INOUT) :: RCU(K,M)
    REAL(kind=r8), INTENT(INOUT) :: PCU(K) 
    REAL(kind=r8), INTENT(INOUT) :: FLX(K+1)
    REAL(kind=r8), INTENT(INOUT) :: FLXD(K+1)
    REAL(kind=r8), INTENT(INOUT) :: CUP
    LOGICAL      , INTENT(IN   ) :: REVAP
    REAL(kind=r8), INTENT(IN   ) :: DT
    REAL(kind=r8), INTENT(INOUT) :: WFNC
    LOGICAL      , INTENT(IN   ) :: WRKFUN
    LOGICAL      , INTENT(IN   ) :: CALKBL
    LOGICAL      , INTENT(IN   ) :: CRTFUN
    REAL(kind=r8), INTENT(INOUT) :: TLA
    LOGICAL      , INTENT(IN   ) ::  DNDRFT
    REAL(kind=r8), INTENT(IN   ) :: DPD
    REAL(kind=r8), INTENT(IN   ) :: fscav_(M)
    REAL(kind=r8), INTENT(IN   ) :: dlq_fac
    REAL(kind=r8), INTENT(IN   ) :: C0
    REAL(kind=r8), INTENT(IN   ) :: C0I 

    LOGICAL       :: CALCUP


    !REAL(kind=r8) :: UFN
    INTEGER       :: KBL
    INTEGER       :: KB1

    !     real(kind=r8) RASALF, FRACBL, MAX_NEG_BOUY, ALFINT(K),     &
    !     real(kind=r8) ALFINQ(K), DPD, alfind(k), rhc_ls(k)

    !  UPDATE ARGUMENTS

    REAL(kind=r8)  ::    TCD(K),   QCD(K)

    !  TEMPORARY WORK SPACE

    REAL(kind=r8) :: HOL(1:K)
    REAL(kind=r8) :: QOL(1:K)
    REAL(kind=r8) :: GAF(1:K+1)
    REAL(kind=r8) :: HST(1:K)
    REAL(kind=r8) :: QST(1:K)
    REAL(kind=r8) :: TOL(1:K)
    REAL(kind=r8) :: GMH(1:K)
    REAL(kind=r8) :: GMS(1:K+1)
    REAL(kind=r8) :: GAM(1:K+1)
    REAL(kind=r8) :: AKT(1:K)
    REAL(kind=r8) :: AKC(1:K)
    REAL(kind=r8) :: BKC(1:K)
    REAL(kind=r8) :: LTL(1:K)
    REAL(kind=r8) :: RNN(1:K)
    REAL(kind=r8) :: FCO(1:K)
    REAL(kind=r8) :: PRI(1:K)
    !REAL(kind=r8) :: PRH(1:K)
    REAL(kind=r8) :: QIL(1:K)
    REAL(kind=r8) :: QLL(1:K)
    REAL(kind=r8) :: ZET(1:K)
    REAL(kind=r8) :: XI(1:K)
    REAL(kind=r8) :: RNS(1:K)
    REAL(kind=r8) :: Q0U(1:K)
    REAL(kind=r8) :: Q0D(1:K)
    REAL(kind=r8) :: vtf(1:K)
    REAL(kind=r8) :: DLB(1:K+1)
    REAL(kind=r8) :: DLT(1:K+1)
    REAL(kind=r8) :: ETA(1:K+1)
    REAL(kind=r8) :: PRL(1:K+1)
    REAL(kind=r8) :: CIL(1:K)
    REAL(kind=r8) :: CLL(1:K)
    REAL(kind=r8) :: ETAI(1:K)

    REAL(kind=r8) :: ALM
    REAL(kind=r8) :: DET
    REAL(kind=r8) :: HCC
    REAL(kind=r8) :: CLP
    REAL(kind=r8) :: HSU
    REAL(kind=r8) :: HSD
    REAL(kind=r8) :: QTL
    REAL(kind=r8) :: QTV
    REAL(kind=r8) :: AKM
    REAL(kind=r8) :: WFN
    REAL(kind=r8) :: HOS
    REAL(kind=r8) :: QOS
    REAL(kind=r8) :: AMB
    REAL(kind=r8) :: TX1
    REAL(kind=r8) :: TX2
    REAL(kind=r8) :: TX3
    REAL(kind=r8) :: TX4
    REAL(kind=r8) :: TX5
    REAL(kind=r8) :: QIS
    REAL(kind=r8) :: QLS
    REAL(kind=r8) :: HBL
    REAL(kind=r8) :: QBL
    REAL(kind=r8) :: RBL(M)
    REAL(kind=r8) :: QLB
    REAL(kind=r8) :: QIB
    REAL(kind=r8) :: PRIS
    !REAL(kind=r8) :: TX6
    REAL(kind=r8) :: ACR
    !REAL(kind=r8) :: TX7
    !EAL(kind=r8) :: TX8
    !REAL(kind=r8) :: TX9
    REAL(kind=r8) :: RHC
    REAL(kind=r8) :: hstkd
    REAL(kind=r8) :: qstkd
    REAL(kind=r8) :: ltlkd
    REAL(kind=r8) :: q0ukd
    REAL(kind=r8) :: q0dkd
    REAL(kind=r8) :: dlbkd
    REAL(kind=r8) :: qtp
    REAL(kind=r8) :: qw00
    REAL(kind=r8) :: qi00
    REAL(kind=r8) :: qrbkd
    REAL(kind=r8) :: hstold
    REAL(kind=r8) :: rel_fac
    REAL(kind=r8) :: prism

    REAL(kind=r8) :: wrk1(1:k), wrk2(1:k)

    !     INTEGER IA,  I1,  I2, ID1, ID2
    !     INTEGER IB,  I3

    LOGICAL :: UNSAT
    LOGICAL ::  ep_wfn

    LOGICAL :: LOWEST
    !LOGICAL ::  SKPDD

    REAL(kind=r8) :: TL
    REAL(kind=r8) :: PL
    REAL(kind=r8) :: QL
    REAL(kind=r8) :: QS
    REAL(kind=r8) :: DQS
    REAL(kind=r8) :: ST1
    !REAL(kind=r8) :: SGN
    REAL(kind=r8) :: TAU
    REAL(kind=r8) :: QTVP
    REAL(kind=r8) :: HB
    REAL(kind=r8) :: QB
    !REAL(kind=r8) :: TB
    !REAL(kind=r8) :: QQQ
    REAL(kind=r8) :: HCCP
    REAL(kind=r8) :: DS
    REAL(kind=r8) :: DH
    REAL(kind=r8) :: AMBMAX
    REAL(kind=r8) :: X00
    REAL(kind=r8) :: EPP
    REAL(kind=r8) :: QTLP
    REAL(kind=r8) :: DPI
    REAL(kind=r8) :: DPHIB
    REAL(kind=r8) :: DPHIT
    REAL(kind=r8) :: DEL_ETA
    REAL(kind=r8) :: DETP
    REAL(kind=r8) :: TEM
    REAL(kind=r8) :: TEM1
    REAL(kind=r8) :: TEM2
    REAL(kind=r8) :: TEM3
    REAL(kind=r8) :: TEM4
    REAL(kind=r8) :: ST2
    REAL(kind=r8) :: ST3
    REAL(kind=r8) :: ST4
    REAL(kind=r8) :: ST5
    !REAL(kind=r8) :: ERRH
    !REAL(kind=r8) :: ERRW
    !REAL(kind=r8) :: ERRE
    REAL(kind=r8) :: TEM5
    REAL(kind=r8) :: TEM6
    REAL(kind=r8) :: HBD
    REAL(kind=r8) :: QBD
    REAL(kind=r8) :: st1s
    !REAL(kind=r8), PARAMETER :: ERRMIN=0.0001_r8!, ERRMI2=0.1_r8*ERRMIN

    !     parameter (c0=1.0e-3_r8, KBLMX=20, ERRMIN=0.0001_r8, ERRMI2=0.1_r8*ERRMIN)

    !INTEGER    :: I
    INTEGER    :: L
    INTEGER    :: N
    INTEGER    :: KD1
    INTEGER    :: II 
    INTEGER    :: KP1
    !INTEGER    :: IT
    INTEGER    :: KM1
    INTEGER    :: KTEM
    !INTEGER    :: KK
    !INTEGER    :: KK1
    !INTEGER    :: LM1
    !INTEGER    :: LL
    !INTEGER    :: LP1
    INTEGER    :: kbls
    INTEGER    :: kmxh
    INTEGER    :: kblh, kblm, kblpmn, kmax, kmaxm1, kmaxp1, klcl, kmin, kmxb

    REAL(kind=r8) ::  avt
    REAL(kind=r8) ::  avq
    REAL(kind=r8) ::  avr
    REAL(kind=r8) ::  avh
    !
    !     REEVAPORATION
    !
    !     real(kind=r8), parameter ::
    !    &                   clfa = -0.452550814376093547E-03_r8
    !    &,                  clfb =  0.161398573159240791E-01_r8
    !    &,                  clfc = -0.163676268676807096_r8
    !    &,                  clfd =  0.447988962175259131_r8
    !    &,                  point3 = 0.3, point01=0.01_r8

    !     real(kind=r8), parameter :: rainmin=1.0e-9_r8
    REAL(kind=r8), PARAMETER :: rainmin=1.0e-8_r8
    !REAL(kind=r8), PARAMETER :: oneopt9=1.0_r8/0.09_r8
    !REAL(kind=r8), PARAMETER :: oneopt4=1.0_r8/0.04_r8

    REAL(kind=r8) :: CLFRAC, clvfr

    REAL(kind=r8) :: ACTEVAP
    !REAL(kind=r8) :: AREARAT
    REAL(kind=r8) :: DELTAQ
    !REAL(kind=r8) :: MASS
    !REAL(kind=r8) :: MASSINV
    REAL(kind=r8) :: POTEVAP  
    REAL(kind=r8) :: TEQ
    REAL(kind=r8) :: QSTEQ
    REAL(kind=r8) :: DQDT
    REAL(kind=r8) :: QEQ
    !
    !     Temporary workspace and parameters needed for downdraft
    !
    REAL(kind=r8) ::  GMF
    !
    REAL(kind=r8) :: BUY(1:K+1)
    REAL(kind=r8) :: QRB(1:K)
    REAL(kind=r8) :: QRT(1:K) 
    REAL(kind=r8) :: ETD(1:K+1)
    REAL(kind=r8) :: HOD(1:K+1)
    REAL(kind=r8) :: QOD(1:K+1) 
    REAL(kind=r8) :: GHD(1:K)
    REAL(kind=r8) :: GSD(1:K)
    REAL(kind=r8) :: EVP(1:K)    
    REAL(kind=r8) :: CLDFR(1:K)
    REAL(kind=r8) :: TRAIN
    REAL(kind=r8) :: DOF
    REAL(kind=r8) :: CLDFRD
    !REAL(kind=r8) :: FAC
    !REAL(kind=r8) :: RSUM1
    !REAL(kind=r8) :: RSUM2
    !REAL(kind=r8) :: RSUM3
    REAL(kind=r8) :: dpneg
    REAL(kind=r8) :: hcrit
    INTEGER       :: IDH
    LOGICAL       :: DDFT
    REAL(kind=r8) :: dhdpmn, dhdp(1:k)
    REAL(kind=r8) :: dlq(1:k),ETZI(1:K-1)
    REAL(kind=r8) :: ETZ(1:K)

!
!  Scavenging related parameters
!
      real(kind=r8)                delzkm          ! layer thickness in km
      real(kind=r8)                fnoscav         ! fraction of tracer *not* scavenged
    REAL(kind=r8) :: hmax,hmin,shal_fac
    INTEGER :: lcon

    !     real(kind=r8) eps, epsm1, rvi, facw, faci, hsub, tmix, DEN
    !     real(kind=r8) eps, epsm1, rv, rd, depth
    !     real(kind=r8) eps, epsm1, rv, rd, fpvs, depth
    !
    !
    !***********************************************************************
    !
    !  UPDATE ARGUMENTS

    TCD=0.0_r8;QCD=0.0_r8

    !  TEMPORARY WORK SPACE

    HOL=0.0_r8
    QOL=0.0_r8
    GAF=0.0_r8
    HST=0.0_r8
    QST=0.0_r8
    TOL=0.0_r8
    GMH=0.0_r8
    GMS=0.0_r8
    GAM=0.0_r8
    AKT=0.0_r8
    AKC=0.0_r8
    BKC=0.0_r8
    LTL=0.0_r8
    RNN=0.0_r8
    FCO=0.0_r8
    PRI=0.0_r8
    !REAL(kind=r8) :: PRH(KD:K)
    QIL=0.0_r8
    QLL=0.0_r8
    ZET=0.0_r8
    XI =0.0_r8
    RNS=0.0_r8
    Q0U=0.0_r8
    Q0D=0.0_r8
    vtf=0.0_r8
    DLB=0.0_r8
    DLT=0.0_r8
    ETA=0.0_r8
    PRL=0.0_r8
    CIL=0.0_r8
    CLL=0.0_r8
    ETAI=0.0_r8

    ALM=0.0_r8
    DET=0.0_r8
    HCC=0.0_r8
    CLP=0.0_r8
    HSU=0.0_r8
    HSD=0.0_r8
    QTL=0.0_r8
    QTV=0.0_r8
    AKM=0.0_r8
    WFN=0.0_r8
    HOS=0.0_r8
    QOS=0.0_r8
    AMB=0.0_r8
    TX1=0.0_r8
    TX2=0.0_r8
    TX3=0.0_r8
    TX4=0.0_r8
    TX5=0.0_r8
    QIS=0.0_r8
    QLS=0.0_r8
    HBL=0.0_r8
    QBL=0.0_r8
    RBL=0.0_r8
    QLB=0.0_r8
    QIB=0.0_r8
    PRIS=0.0_r8
    ACR=0.0_r8
    RHC=0.0_r8
    hstkd=0.0_r8
    qstkd=0.0_r8
    ltlkd=0.0_r8
    q0ukd=0.0_r8
    q0dkd=0.0_r8
    dlbkd=0.0_r8
    qtp=0.0_r8
    qw00=0.0_r8
    qi00=0.0_r8
    qrbkd=0.0_r8
    hstold=0.0_r8
    rel_fac=0.0_r8

    !     INTEGER IA,  I1,  I2, ID1, ID2
    !     INTEGER IB,  I3


    !LOGICAL ::  SKPDD

    TL=0.0_r8
    PL=0.0_r8
    QL=0.0_r8
    QS=0.0_r8
    DQS=0.0_r8
    ST1=0.0_r8
    TAU=0.0_r8
    QTVP=0.0_r8
    HB=0.0_r8
    QB=0.0_r8
    HCCP=0.0_r8
    DS=0.0_r8
    DH=0.0_r8
    AMBMAX=0.0_r8
    X00=0.0_r8
    EPP=0.0_r8
    QTLP=0.0_r8
    DPI=0.0_r8
    DPHIB=0.0_r8
    DPHIT=0.0_r8
    DEL_ETA=0.0_r8
    DETP=0.0_r8
    TEM=0.0_r8
    TEM1=0.0_r8
    TEM2=0.0_r8
    TEM3=0.0_r8
    TEM4=0.0_r8
    ST2=0.0_r8
    ST3=0.0_r8
    ST4=0.0_r8
    ST5=0.0_r8
    TEM5=0.0_r8
    TEM6=0.0_r8
    HBD=0.0_r8
    QBD=0.0_r8
    st1s=0.0_r8

    L=0
    N=0
    KD1=0
    II =0
    KP1=0
    KM1=0
    KTEM=0
    kbls=0
    kmxh=0

    avt=0.0_r8
    avq=0.0_r8
    avr=0.0_r8
    avh=0.0_r8
    CLFRAC=0.0_r8
    ACTEVAP=0.0_r8
    DELTAQ=0.0_r8
    POTEVAP =0.0_r8 
    TEQ=0.0_r8
    QSTEQ=0.0_r8
    DQDT=0.0_r8
    QEQ=0.0_r8
    !
    !     Temporary workspace and parameters needed for downdraft
    !
    !
    BUY=0.0_r8
    QRB=0.0_r8
    QRT =0.0_r8
    ETD=0.0_r8
    HOD=0.0_r8
    QOD =0.0_r8
    GHD=0.0_r8
    GSD=0.0_r8
    EVP =0.0_r8
    ETZ=0.0_r8
    CLDFR=0.0_r8
    TRAIN=0.0_r8
    DOF=0.0_r8
    CLDFRD=0.0_r8
    dpneg=0.0_r8
    IDH=0

    !
    DO l=1,K
       tcd(L) = 0.0_r8
       qcd(L) = 0.0_r8
    ENDDO
    !
    KP1     = K  + 1
    KM1     = K  - 1
    KD1     = KD + 1
    !kblmx   = k
    !
    !     if (lprnt) print *,' IN CLOUD for KD=',kd
    !     if (lprnt) print *,' prs=',prs(Kd:K+1)
    !     if (lprnt) print *,' phil=',phil(KD:K)
    !     if (lprnt) print *,' phih=',phih(KD:K+1)
    !     if (lprnt) print *,' toi=',toi
    !     if (lprnt) print *,' qoi=',qoi
    !
    !
    CLDFRD   = 0.0_r8
    DOF      = 0.0_r8
    PRL(KP1) = PRS(KP1)
    !
    DO L=KD,K
       RNN(L) = 0.0_r8
       ZET(L) = 0.0_r8
       XI(L)  = 0.0_r8
       !
       TOL(L) = TOI(L)
       QOL(L) = QOI(L)
       PRL(L) = PRS(L)
       BUY(L) = 0.0_r8
       CLL(L) = QLI(L)
       CIL(L) = QII(L)
    ENDDO

    if (vsmooth) then
       do l=kd,k
         wrk1(l) = tol(l)
         wrk2(l) = qol(l)
       enddo
       do l=kd1,km1
          tol(l) = 0.25*wrk1(l-1) + 0.5*wrk1(l) + 0.25*wrk1(l+1)
          qol(l) = 0.25*wrk2(l-1) + 0.5*wrk2(l) + 0.25*wrk2(l+1)
       enddo
    endif

    !
    DO L=KD, K
       DPI    = ONE / (PRL(L+1) - PRL(L))
       PRI(L) = GRAVFAC * DPI
       !
       PL     = PRSM(L)
       TL     = TOL(L)

       AKT(L) = (PRL(L+1) - PL) * DPI
       !
       CALL QSATCN(TL, PL, QS, DQS)

       QST(L) = QS
       GAM(L) = DQS * ELOCP
       ST1    = ONE + GAM(L)
       GAF(L) = (ONE/con_hvap ) * (GAM(L)/(ONE + GAM(L)))

       QL     = MAX(MIN(QS*RHMAX,QOL(L)), ONE_M10)
       QOL(L) = QL

       TEM    = con_cp * TL
       LTL(L) = TEM * ST1 / (ONE+con_FVirt*(QST(L)+TL*DQS))
       vtf(L) = 1.0_r8 + con_FVirt * QL
       ETA(L) = ONE / (LTL(L) * VTF(L))

       HOL(L) = TEM + QL * con_hvap 
       HST(L) = TEM + QS * con_hvap 
       !
    ENDDO
    !
    ETA(K+1) = ZERO
    GMS(K)   = ZERO
    !
    AKT(KD)  = HALF
    GMS(KD)  = ZERO
    !
    CLP      = ZERO
    !
    GAM(K+1) = GAM(K)
    GAF(K+1) = GAF(K)
    !
    DO L=K,KD1,-1
       DPHIB  = PHIL(L) - PHIH(L+1)
       DPHIT  = PHIH(L) - PHIL(L)
       !
       DLB(L) = DPHIB * ETA(L)
       DLT(L) = DPHIT * ETA(L)
       !
       QRB(L) = DPHIB
       QRT(L) = DPHIT
       !
       ETA(L) = ETA(L+1) + DPHIB

       !
       HOL(L) = HOL(L) + ETA(L)
       hstold = hst(l)
       HST(L) = HST(L) + ETA(L)
       !
       ETA(L) = ETA(L) + DPHIT
    ENDDO
    !
    !     For the cloud top layer
    !
    L = KD

    DPHIB  = PHIL(L) - PHIH(L+1)
    !
    DLB(L) = DPHIB * ETA(L)
    !
    QRB(L) = DPHIB
    QRT(L) = DPHIB
    !
    ETA(L) = ETA(L+1) + DPHIB

    HOL(L) = HOL(L) + ETA(L)
    HST(L) = HST(L) + ETA(L)
    !
    !     if (kd .eq. 12) then
    !     if (lprnt) print *,' IN CLOUD for KD=',KD,' K=',K
    !     if (lprnt) print *,' l=',l,' hol=',hol(l),' hst=',hst(l)
    !     if (lprnt) print *,' TOL=',tol
    !     if (lprnt) print *,' qol=',qol
    !     if (lprnt) print *,' hol=',hol
    !     if (lprnt) print *,' hst=',hst
    !     endif
    !
    !     To determine KBL internally -- If KBL is defined externally
    !     the following two loop should be skipped
    !
    !     if (lprnt) print *,' calkbl=',calkbl

    hcrit = hcritd
    IF (CALKBL) THEN
         KTEM = MAX(KD+1, KBLMX)
         hmin = hol(k)
         kmin = k
         do l=km1,kd,-1
           if (hmin > hol(l)) then
             hmin = hol(l)
             kmin = l
           endif
         enddo
         if (kmin == k) return
         hmax = hol(k)
         kmax = k
         do l=km1,ktem,-1
           if (hmax < hol(l)) then
             hmax = hol(l)
             kmax = l
           endif
         enddo
         kmxb = kmax
         if (kmax < kmin) then
           kmax = k
           kmxb = k
           hmax = hol(kmax)
         elseif (kmax < k) then
           do l=kmax+1,k
             if (abs(hol(kmax)-hol(l)) > 0.5 * hcrit) then
               kmxb = l - 1
               exit
             endif
           enddo
         endif
         kmaxm1 = kmax - 1
         kmaxp1 = kmax + 1
         kblpmn = kmax
!
         dhdp(kmax:k) = 0.0
         dhdpmn = dhdp(kmax)
         do l=kmaxm1,ktem,-1
           dhdp(l) = (HOL(L)-HOL(L+1)) / (PRL(L+2)-PRL(L))
           if (dhdp(l) < dhdpmn) then
             dhdpmn = dhdp(l)
             kblpmn = l + 1
           elseif (dhdp(l) > 0.0 .and. l <= kmin) then
             exit
           endif
         enddo
         kbl = kmax
         if (kblpmn < kmax) then
           do l=kblpmn,kmaxm1
             if (hmax-hol(l) < 0.5*hcrit) then
               kbl = l
               exit
             endif
           enddo
         endif
       
!     if(lprnt) print *,' kbl=',kbl,' kbls=',kbls,' kmax=',kmax
!
         klcl = kd1
         if (kmax > kd1) then
           do l=kmaxm1,kd1,-1
             if (hmax > hst(l)) then
               klcl = l+1
               exit
             endif
           enddo
         endif
!        if(lprnt) print *,' klcl=',klcl,' ii=',ii
!        if (klcl == kd .or. klcl < ktem) return

!        This is to handle mid-level convection from quasi-uniform h

         if (kmax < kmxb) then
           kmax   = max(kd1, min(kmxb,k))
           kmaxm1 = kmax - 1
           kmaxp1 = kmax + 1
         endif


!        if (prl(Kmaxp1) - prl(klcl) > 250.0 ) return

         ii  = max(kbl,kd1)
         kbl = max(klcl,kd1)
         tem = min(50.0,max(10.0,(prl(kmaxp1)-prl(kd))*0.10))
         if (prl(kmaxp1) - prl(ii) > tem .and. ii > kbl) kbl = ii

!        if(lprnt) print *,' kbl2=',kbl,' ii=',ii
       IF (kbl .NE. ii) THEN
          IF (PRL(K+1)-PRL(KBL) .GT. bldmax) kbl = MAX(kbl,ii)
       ENDIF
       if (kbl < ii) then
         if (hol(ii)-hol(ii-1) > 0.5*hcrit) kbl = ii
       endif
       !
       !PK TESTE if (prl(kbl) - prl(klcl) > 250.0 ) return   
       if (prl(kbl) - prl(klcl) > 250.0 ) return

       KBL  = MIN(k, MAX(KBL,K-KBLMX))
       !KBL  = min(KPBL,KBL)!PK
       !KBL  = max(KPBL,KBL)!PK TESTE
       ! 
!       tem1 = MAX(prl(k+1)-prl(k),                                    &
!            &                     MIN((prl(kbl) - prl(kd))*0.05_r8, 20.0_r8))
!       !    &                     min((prl(kbl) - prl(kd))*0.05_r8, 30.0_r8))
!       IF (prl(k+1)-prl(kbl) .LT. tem1) THEN
!          KTEM = MAX(KD+1, K-KBLMX)
!          DO l=k,KTEM,-1
!             tem = prl(k+1) - prl(l)
!             IF (tem .GT. tem1) THEN
!                kbl = MIN(kbl,l)
!                EXIT
!!             ENDIF
!          ENDDO
!       ENDIF
!        if (kbl == kblmx .and. kmax >= k-1) kbl = k - 1

       KPBL = KBL
       !     if(lprnt)print*,' 1st kbl=',kbl,' kblmx=',kblmx,' kd=',kd
       !     if(lprnt)print*,' tx3=',tx3,' tx1=',tx1,' tem=',tem
       !    1,               ' hcrit=',hcrit
    ELSE
       KBL  = KPBL
       Kmaxp1=k+1
       Kmax=k
       kmaxm1=k-1
       !     if(lprnt)print*,' 2nd kbl=',kbl
    ENDIF

    !     if(lprnt)print*,' after CALKBL l=',l,' hol=',hol(l)
    !    1,               ' hst=',hst(l)
    !
    KBL      = min(kmax,MAX(KBL,KD+2))
    KB1      = KBL - 1
!!
!     if (lprnt) print *,' kbl=',kbl,' prlkbl=',prl(kbl),prl(k+1)

      if(PRL(Kmaxp1)-PRL(KBL) .gt. bldmax .or. kb1 .le. kd) then
        return
      endif
    !
    !     if (lprnt) print *,' kbl=',kbl
!     write(0,*)' kbl=',kbl,' kmax=',kmax,' kmaxp1=',kmaxp1,' k=',k
    !
    PRIS     = ONE / (PRL(K+1)-PRL(KBL))
    PRISM    = ONE / (PRL(Kmaxp1)-PRL(KBL))
    TX1      = ETA(KBL)
    !
    GMS(KBL) = 0.0_r8
    XI(KBL)  = 0.0_r8
    ZET(KBL) = 0.0_r8
    !
    shal_fac = 1.0
!   if (prl(kbl)-prl(kd) < 300.0 .and. kmax == k) shal_fac = shalfac
    if (prl(kbl)-prl(kd) < 350.0 .and. kmax == k) shal_fac = shalfac
    DO L=K,KD,-1
       IF (L .GE. KBL) THEN
          ETA(L) = (PRL(K+1)-PRL(L)) * PRIS
       ELSE
          ZET(L) = (ETA(L) - TX1) * ONEBG
          !XI(L)  =  ZET(L) * ZET(L) * QUDFAC
          XI(L)  =  ZET(L) * ZET(L) * (QUDFAC*shal_fac)
          ETA(L) =  ZET(L) - ZET(L+1)
          GMS(L) =  XI(L)  - XI(L+1)
       ENDIF
!       if (lprnt) print *,' l=',l,' eta=',eta(l),' kbl=',kbl
    ENDDO
      if (kmax < k) then
        do l=kmaxp1,kp1
          eta(l) = 0.0
        enddo
      endif

    !
    HBL = HOL(kmax) * ETA(kmax)
    QBL = QOL(kmax) * ETA(kmax)
    QLB = CLL(kmax) * ETA(kmax)
    QIB = CIL(kmax) * ETA(kmax)
    TX1 = QST(kmax) * ETA(kmax)
    !
    DO L=Kmaxm1,KBL,-1
       TEM = ETA(L) - ETA(L+1)
       HBL = HBL + HOL(L) * TEM
       QBL = QBL + QOL(L) * TEM
       QLB = QLB + CLL(L) * TEM
       QIB = QIB + CIL(L) * TEM
       TX1 = TX1 + QST(L) * TEM
    ENDDO

!     if (ctei .and. sgcs(ipt,l) > 0.65) then
!        hbl = hbl * hpert_fac
!        qbl = qbl * hpert_fac
!     endif

    !     if (lprnt) print *,' hbl=',hbl,' qbl=',qbl
    !                                   Find Min value of HOL in TX2
    TX2 = HOL(KD)
    IDH = KD1
    DO L=KD1,KB1
       IF (HOL(L) .LT. TX2) THEN
          TX2 = HOL(L)
          IDH = L             ! Level of minimum moist static energy!
       ENDIF
    ENDDO
    IDH = 1
    IDH = MAX(KD1, IDH)
    !
    TEM1 = HBL - HOL(KD)
    TEM  = HBL - HST(KD1) - LTL(KD1) *( con_FVirt *(QOL(KD1)-QST(KD1)))
    LOWEST = KD .EQ. KB1
    !
      lcon = kd
      do l=kb1,kd1,-1
        if (hbl >= hst(l)) then
          lcon = l
          exit
        endif
      enddo
!
      if (lcon == kd .or. kbl <= kd .or. prl(kbl)-prsm(lcon) > 150.0)   &
     &                                    return

    TX1   = RHFACS - QBL / TX1       !     Average RH

    UNSAT = (TEM .GT. ZERO .OR. (LOWEST .AND. TEM1 .GE. ZERO)).AND. (TX1 .LT. RHRAM)

    !     if(lprnt) print *,' unsat=',unsat,' tem=',tem,' tem1=',tem1
    !    &,' tx1=',tx1,' rhram=',rhram,' kbl=',kbl,' kd=',kd,' lowest='
    !    &,lowest,' rhfacs=',rhfacs,' ltl=',ltl(kd1),' qol=',qol(kd1)
    !    &,' qst=',qst(kd1),' hst=',hst(kd1),' con_FVirt=',con_FVirt

    !
    !===>  IF NO SOUNDING MEETS FIRST CONDITION, RETURN
    !     if(lprnt .and. (.not. unsat)) print *,' tx1=',tx1,' rhfacs='
    !    &,rhfacs, ' tem=',tem,' hst=',hst(kd1)

    IF (.NOT. UNSAT) RETURN
    !
    RHC    = MAX(ZERO, MIN(ONE, EXP(-20.0_r8*TX1) ))
    !
    if (m > 0) then
       DO N=1,M
          RBL(N) = ROI(Kmax,N) * ETA(Kmax)
       ENDDO
       DO N=1,M
          DO L=KmaxM1,KBL,-1
             RBL(N) = RBL(N) + ROI(L,N)*(ETA(L)-ETA(L+1))
          ENDDO
       ENDDO
    endif
    !
    TX4    = 0.0_r8
    TX5    = 0.0_r8
    !
    TX3      = QST(KBL) - GAF(KBL) * HST(KBL)
    QIL(KBL) = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KBL))*TCRF))
    !
    DO L=KB1,KD1,-1
       TEM      = QST(L) - GAF(L) * HST(L)
       TEM1     = (TX3 + TEM) * 0.5_r8
       ST2      = (GAF(L)+GAF(L+1)) * 0.5_r8
       !
       FCO(L+1) =            TEM1 + ST2 * HBL

       !     if(lprnt) print *,' fco=',fco(l+1),' tem1=',tem1,' st2=',st2
       !    &,' hbl=',hbl,' tx3=',tx3,' tem=',tem,' gaf=',gaf(l),' l=',l

       RNN(L+1) = ZET(L+1) * TEM1 + ST2 * TX4
       GMH(L+1) = XI(L+1)  * TEM1 + ST2 * TX5
       !
       TX3      = TEM
       TX4      = TX4 + ETA(L) * HOL(L)
       TX5      = TX5 + GMS(L) * HOL(L)
       !
       QIL(L)   = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(L))*TCRF))
       QLL(L+1) = (0.5_r8*con_hfus) * ST2 * (QIL(L)+QIL(L+1)) + ONE
    ENDDO
    !
    !     FOR THE CLOUD TOP -- L=KD
    !
    L = KD
    !
    TEM      = QST(L) - GAF(L) * HST(L)
    TEM1     = (TX3 + TEM) * 0.5_r8
    ST2      = (GAF(L)+GAF(L+1)) * 0.5_r8
    !
    FCO(L+1) =            TEM1 + ST2 * HBL
    RNN(L+1) = ZET(L+1) * TEM1 + ST2 * TX4
    GMH(L+1) = XI(L+1)  * TEM1 + ST2 * TX5
    !
    FCO(L)   = TEM + GAF(L) * HBL
    RNN(L)   = TEM * ZET(L) + (TX4 + ETA(L)*HOL(L)) * GAF(L)
    GMH(L)   = TEM * XI(L)  + (TX5 + GMS(L)*HOL(L)) * GAF(L)
    !
    !   Replace FCO for the Bottom
    !
    FCO(KBL) = QBL
    RNN(KBL) = 0.0_r8
    GMH(KBL) = 0.0_r8
    !
    QIL(KD)  =  MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KD))*TCRF))
    QLL(KD1) = (0.5_r8*con_hfus) * ST2 * (QIL(KD) + QIL(KD1)) + ONE
    QLL(KD ) = con_hfus * GAF(KD) * QIL(KD) + ONE
    !
    !     if (lprnt) print *,' fco=',fco(kd:kbl)
    !     if (lprnt) print *,' qil=',qil(kd:kbl)
    !     if (lprnt) print *,' qll=',qll(kd:kbl)
    !
    st1  = qil(kd)
    st2  = c0i * st1
    tem  = c0  * (1.0_r8-st1)
    tem2 = st2*qi0 + tem*qw0
    !
    DO L=KD,KB1
       tx2    = akt(l) * eta(l)
       tx1    = tx2 * tem2
       q0u(l) = tx1
       FCO(L) = FCO(L+1) - FCO(L) + tx1
       RNN(L) = RNN(L+1) - RNN(L)                                     &
            &          + ETA(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*zet(l)
       GMH(L) = GMH(L+1) - GMH(L)                                     &
            &          + GMS(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*xi(l)
       !
       tem1   = (1.0_r8-akt(l)) * eta(l)

       !     if(lprnt) print *,' qll=',qll(l),' st2=',st2,' tem=',tem
       !    &,' tx2=',tx2,' akt=',akt(l),' eta=',eta(l)

       AKT(L) = QLL(L)   + (st2 + tem) * tx2

       !     if(lprnt) print *,' akt==',akt(l),' l==',l

       AKC(L) = 1.0_r8 / AKT(L)
       !
       st1    = 0.5_r8 * (qil(l)+qil(l+1))
       st2    = c0i * st1
       tem    = c0  * (1.0_r8-st1)
       tem2   = st2*qi0 + tem*qw0
       !
       BKC(L) = QLL(L+1) - (st2 + tem) * tem1
       !
       tx1    = tem1*tem2
       q0d(l) = tx1
       FCO(L) = FCO(L) + tx1
       RNN(L) = RNN(L) + tx1*zet(l+1)
       GMH(L) = GMH(L) + tx1*xi(l+1)
    ENDDO

    !     if(lprnt) print *,' akt=',akt(kd:kb1)
    !     if(lprnt) print *,' akc=',akc(kd:kb1)

    qw00 = qw0
    qi00 = qi0
    ii = 0
777 CONTINUE
    !
    !     if (lprnt) print *,' after 777 ii=',ii,' ep_wfn=',ep_wfn
    !
    ep_wfn = .FALSE.
    RNN(KBL) = 0.0_r8
    TX3      = bkc(kb1) * (QIB + QLB)
    TX4      = 0.0_r8
    TX5      = 0.0_r8
    DO L=KB1,KD1,-1
       TEM    = BKC(L-1)       * AKC(L)
       !     if (lprnt) print *,' tx3=',tx3,' fco=',fco(l),' akc=',akc(l)
       !    &,' bkc=',bkc(l-1), ' l=',l
       TX3    = (TX3 + FCO(L)) * TEM
       TX4    = (TX4 + RNN(L)) * TEM
       TX5    = (TX5 + GMH(L)) * TEM
    ENDDO
    IF (KD .LT. KB1) THEN
       HSD   = HST(KD1) + LTL(KD1) *  con_FVirt *(QOL(KD1)-QST(KD1))
    ELSE
       HSD   = HBL
    ENDIF
    !
    !     if (lprnt) print *,' tx3=',tx3,' fco=',fco(kd),' akc=',akc(kd)
    TX3 = (TX3 + FCO(KD)) * AKC(KD)
    TX4 = (TX4 + RNN(KD)) * AKC(KD)
    TX5 = (TX5 + GMH(KD)) * AKC(KD)
    ALM = con_hfus*QIL(KD) - LTL(KD) * VTF(KD)
    !
    HSU = HST(KD) + LTL(KD) * con_FVirt * (QOL(KD)-QST(KD))

    !     if (lprnt) print *,' hsu=',hsu,' hst=',hst(kd),
    !    &' ltl=',ltl(kd),' qol=',qol(kd),' qst=',qst(kd)
    !
    !===> VERTICAL INTEGRALS NEEDED TO COMPUTE THE ENTRAINMENT PARAMETER
    !
    TX1 = ALM * TX4
    TX2 = ALM * TX5

    DO L=KD,KB1
       TAU = HOL(L) - HSU
       TX1 = TX1 + TAU * ETA(L)
       TX2 = TX2 + TAU * GMS(L)
    ENDDO
    !
    !     MODIFY HSU TO INCLUDE CLOUD LIQUID WATER AND ICE TERMS
    !
    !     if (lprnt) print *,' hsu=',hsu,' alm=',alm,' tx3=',tx3

    HSU   = HSU - ALM * TX3
    !
    CLP   = ZERO
    ALM   = -100.0_r8
    HOS   = HOL(KD)
    QOS   = QOL(KD)
    QIS   = CIL(KD)
    QLS   = CLL(KD)
    UNSAT = HBL .GT. HSU .AND. ABS(tx1) .GT. 1.0e-4_r8

    !     if (lprnt) print *,' ii=',ii,' unsat=',unsat,' hsu=',hsu
    !    &,' hbl=',hbl,' tx1=',tx1,' hsd=',hsd


    !***********************************************************************


    ST1  = HALF*(HSU + HSD)
    IF (UNSAT) THEN
       !
       !  STANDARD CASE:
       !   CLOUD CAN BE NEUTRALLY BOUYANT AT MIDDLE OF LEVEL KD W/ +VE LAMBDA.
       !   EPP < .25 IS REQUIRED TO HAVE REAL ROOTS.
       !
       clp = 1.0_r8
       st2 = hbl - hsu

       !     if(lprnt) print *,' tx2=',tx2,' tx1=',tx1,' st2=',st2
       !
       IF (tx2 .EQ. 0.0_r8) THEN
          alm = - st2 / tx1
          IF (alm .GT. almax) alm = -100.0_r8
       ELSE
          x00 = tx2 + tx2
          epp = tx1 * tx1 - (x00+x00)*st2
          IF (epp .GT. 0.0_r8) THEN
             x00  = 1.0_r8 / x00
             tem  = SQRT(epp)
             tem1 = (-tx1-tem)*x00
             tem2 = (-tx1+tem)*x00
             IF (tem1 .GT. almax) tem1 = -100.0_r8
             IF (tem2 .GT. almax) tem2 = -100.0_r8
             alm  = MAX(tem1,tem2)

             !     if (lprnt) print *,' tem1=',tem1,' tem2=',tem2,' alm=',alm
             !    &,' tx1=',tx1,' tem=',tem,' epp=',epp,' x00=',x00,' st2=',st2

          ENDIF
       ENDIF

       !     if (lprnt) print *,' almF=',alm,' ii=',ii,' qw00=',qw00
       !    &,' qi00=',qi00
       !
       !  CLIP CASE:
       !   NON-ENTRAINIG CLOUD DETRAINS IN LOWER HALF OF TOP LAYER.
       !   NO CLOUDS ARE ALLOWED TO DETRAIN BELOW THE TOP LAYER.
       !
    ELSEIF ( (HBL .LE. HSU) .AND.                                    &
         &          (HBL .GT. ST1   )     ) THEN
       ALM = ZERO
!        CLP = (HBL-ST1) / (HSU-ST1)    ! commented on Jan 16, 2010
    ENDIF
    !
    UNSAT = .TRUE.
    IF (ALMIN1 .GT. 0.0_r8) THEN
       IF (ALM .GE. ALMIN1) UNSAT = .FALSE.
    ELSE
       LOWEST   = KD .EQ. KB1
       IF ( (ALM .GT. ZERO) .OR.                                       &
            &      (.NOT. LOWEST .AND. ALM .EQ. ZERO) ) UNSAT = .FALSE.
    ENDIF
    !
    !===>  IF NO SOUNDING MEETS SECOND CONDITION, RETURN
    !
    IF (UNSAT) THEN
       IF (ii .GT. 0 .OR. (qw00 .EQ. 0.0_r8 .AND. qi00 .EQ. 0.0_r8)) RETURN
       CLP = 1.0_r8
       ep_wfn = .TRUE.
       GO TO 888
    ENDIF
    !
    !     if (lprnt) print *,' hstkd=',hst(kd),' qstkd=',qst(kd)
    !    &,' ii=',ii,' clp=',clp

    st1s = ONE
    IF(CLP.GT.ZERO .AND. CLP.LT.ONE) THEN
       ST1     = HALF*(ONE+CLP)
       ST2     = ONE - ST1
       st1s    = st1
       hstkd   = hst(kd)
       qstkd   = qst(kd)
       ltlkd   = ltl(kd)
       q0ukd   = q0u(kd)
       q0dkd   = q0d(kd)
       dlbkd   = dlb(kd)
       qrbkd   = qrb(kd)
       !
       HST(KD) = HST(KD)*ST1 + HST(KD1)*ST2
       HOS     = HOL(KD)*ST1 + HOL(KD1)*ST2
       QST(KD) = QST(KD)*ST1 + QST(KD1)*ST2
       QOS     = QOL(KD)*ST1 + QOL(KD1)*ST2
       QLS     = CLL(KD)*ST1 + CLL(KD1)*ST2
       QIS     = CIL(KD)*ST1 + CIL(KD1)*ST2
       LTL(KD) = LTL(KD)*ST1 + LTL(KD1)*ST2
       !
       DLB(KD) = DLB(KD)*CLP
       qrb(KD) = qrb(KD)*CLP
       ETA(KD) = ETA(KD)*CLP
       GMS(KD) = GMS(KD)*CLP
       Q0U(KD) = Q0U(KD)*CLP
       Q0D(KD) = Q0D(KD)*CLP
    ENDIF
    !
    !
    !***********************************************************************
    !
    !    Critical workfunction is included in this version
    !
    ACR = 0.0_r8
    TEM = PRL(KD1) - (PRL(KD1)-PRL(KD)) * CLP * HALF
    tx1 = PRL(KBL) - TEM
    tx2 = MIN(900.0_r8,MAX(tx1,100.0_r8))
    tem1    = LOG(tx2*0.01_r8) / LOG(10.0_r8)
    !rel_fac = (dt * facdt)  / (3600.0_r8 * (tem1*3.0_r8 + (1-tem1)*1.0_r8))
!-----------------------------
       ! grell mask
       !      mask2(i)=0 ! land

    IF(mask2 == 0)THEN
          !      mask2(i)=0 ! land
       if ( kdt == 1 ) then
           rel_fac = (dt * facdt)  / (tem1*12.0_r8 + (1-tem1)*3.0_r8)
       else
           rel_fac = (dt * facdt) / (tem1*adjts_d + (1-tem1)*adjts_s)
       endif
    ELSE
       !      mask2(i)=1 ! water/ocean
       if ( kdt == 1 ) then
           rel_fac = (dt * facdt)  / (tem1*12.0_r8 + (1-tem1)*3.0_r8)
       else
           rel_fac = (dt * facdt) / (tem1*adjts_d_ocean + (1-tem1)*adjts_s_ocean)
       endif
    END IF

    !rel_fac = MAX(zero, MIN(one,rel_fac))
    rel_fac = max(zero, min(half,rel_fac))

    IF (CRTFUN) THEN
       CALL CRTWRK(TEM, CCWF, ST1)
       ACR = TX1 * ST1
    ENDIF
    !
    !===>  NORMALIZED MASSFLUX
    !
    !  ETA IS THE THICKNESS COMING IN AND THE MASS FLUX GOING OUT.
    !  GMS IS THE THICKNESS OF THE SQUARE; IT IS LATER REUSED FOR GAMMA_S
    !
    !     ETA(K) = ONE

    DO L=KB1,KD,-1
       ETA(L)  = ETA(L+1) + ALM * (ETA(L) + ALM * GMS(L))
    ENDDO
    DO L=KD,KBL
       ETAI(L) = 1.0_r8 / ETA(L)
    ENDDO

    !     if (lprnt) print *,' eta=',eta,' ii=',ii,' alm=',alm
    !
    !===>  CLOUD WORKFUNCTION
    !
    WFN   = ZERO
    AKM   = ZERO
    DET   = ZERO
    HCC   = HBL
    UNSAT = .FALSE.
    QTL   = QST(KB1) - GAF(KB1)*HST(KB1)
    TX1   = HBL
    !!
    qtv   = qbl
    det   = qlb + qib
    !
    tx2   = 0.0_r8
    dpneg = 0.0_r8
    !
    DO L=KB1,KD1,-1
       DEL_ETA = ETA(L) - ETA(L+1)
       HCCP = HCC + DEL_ETA*HOL(L)
       !
       QTLP = QST(L-1) - GAF(L-1)*HST(L-1)
       QTVP = 0.5_r8 * ((QTLP+QTL)*ETA(L)                                &
            &              + (GAF(L)+GAF(L-1))*HCCP)
       ST1  = ETA(L)*Q0U(L) + ETA(L+1)*Q0D(L)
       DETP = (BKC(L)*DET - (QTVP-QTV)                                &
            &        + DEL_ETA*(QOL(L)+CLL(L)+CIL(L)) + ST1)  * AKC(L)

       !     if(lprnt) print *,' detp=',detp,' bkc=',bkc(l),' det=',det
       !     if (lprnt .and. kd .eq. 15) 
       !    &          print *,' detp=',detp,' bkc=',bkc(l),' det=',det
       !    &,' qtvp=',qtvp,' qtv=',qtv,' del_eta=',del_eta,' qol='
       !    &,qol(l),' st1=',st1,' akc=',akc(l)
       !
       TEM1   = AKT(L)   - QLL(L)
       TEM2   = QLL(L+1) - BKC(L)
       RNS(L) = TEM1*DETP  + TEM2*DET - ST1

       qtp    = 0.5_r8 * (qil(L)+qil(L-1))
       tem2   = MIN(qtp*(detp-eta(l)*qw00),                           &
            &               (1.0_r8-qtp)*(detp-eta(l)*qi00))
       st1    = MIN(tx2,tem2)
       tx2    = tem2
       !
       IF (rns(l) .LT. zero .OR. st1 .LT. zero) ep_wfn = .TRUE.
       IF (DETP .LE. ZERO) UNSAT = .TRUE.

       ST1  = HST(L) - LTL(L)*con_FVirt*(QST(L)-QOL(L))


       TEM2 = HCCP   + DETP   * QTP * con_hfus
       !
       !     if(lprnt) print *,' hst=',hst(l),' ltl=',ltl(l),' con_FVirt=',con_FVirt
       !     if (lprnt .and. kd .eq. 15) 
       !    &          print *,' hst=',hst(l),' ltl=',ltl(l),' con_FVirt=',con_FVirt
       !    &,' qst=',qst(l),' qol=',qol(l),' hccp=',hccp,' detp=',detp
       !    *,' qtp=',qtp,' con_hfus=',con_hfus,' vtf=',vtf(l)

       ST2  = LTL(L) * VTF(L)
       TEM5 = CLL(L) + CIL(L)
       TEM3 = (TX1  - ETA(L+1)*ST1 - ST2*(DET-TEM5*eta(l+1))) * DLB(L)
       TEM4 = (TEM2 - ETA(L  )*ST1 - ST2*(DETP-TEM5*eta(l)))  * DLT(L)
       !
       !     if (lprnt) then
       !     if (lprnt .and. kd .eq. 12) then 
       !       print *,' tem3=',tem3,' tx1=',tx1,' st1=',st1,' eta1=',eta(l+1)
       !    &, ' st2=',st2,' det=',det,' tem5=',tem5,' dlb=',dlb(l)
       !       print *,' tem4=',tem4,' tem2=',tem2,' detp=',detp
       !    &, ' eta=',eta(l),' dlt=',dlt(l),' rns=',rns(l),' l=',l
       !       print *,' bt1=',tem3/(eta(l+1)*qrb(l))
       !    &,         ' bt2=',tem4/(eta(l)*qrt(l))
       !      endif

       ST1  = TEM3 + TEM4

       !     if (lprnt) print *,' wfn=',wfn,' st1=',st1,' l=',l,' ep_wfn=',
       !    &ep_wfn,' akm=',akm

       WFN = WFN + ST1       
       AKM = AKM - MIN(ST1,ZERO)

       !     if (lprnt) print *,' wfn=',wfn,' akm=',akm

       IF (st1 .LT. zero .AND. wfn .LT. zero) THEN
          dpneg = dpneg + prl(l+1) - prl(l)
       ENDIF

       BUY(L) = 0.5_r8 * (tem3/(eta(l+1)*qrb(l)) + tem4/(eta(l)*qrt(l)))
       !
       HCC = HCCP
       DET = DETP
       QTL = QTLP
       QTV = QTVP
       TX1 = TEM2

    ENDDO

    DEL_ETA = ETA(KD) - ETA(KD1)
    HCCP    = HCC + DEL_ETA*HOS
    !
    QTLP = QST(KD) - GAF(KD)*HST(KD)
    QTVP = QTLP*ETA(KD) + GAF(KD)*HCCP
    ST1  = ETA(KD)*Q0U(KD) + ETA(KD1)*Q0D(KD)
    DETP = (BKC(KD)*DET - (QTVP-QTV)                                  &
         &     + DEL_ETA*(QOS+QLS+QIS) + ST1) * AKC(KD)
    !
    TEM1    = AKT(KD)  - QLL(KD)
    TEM2    = QLL(KD1) - BKC(KD)
    RNS(KD) = TEM1*DETP  + TEM2*DET - ST1
    !
    IF (rns(kd) .LT. zero) ep_wfn = .TRUE.
    IF (DETP.LE.ZERO) UNSAT = .TRUE.
    !
888 CONTINUE

    !     if (lprnt) print *,' ep_wfn=',ep_wfn,' ii=',ii,' rns=',rns(kd)
    !    &,' clp=',clp,' hst(kd)=',hst(kd)

    IF (ep_wfn) THEN
       IF ((qw00 .EQ. 0.0_r8 .AND. qi00 .EQ. 0.0_r8)) RETURN
       IF (ii .EQ. 0) THEN
          ii  = 1
          IF (clp .GT. 0.0_r8 .AND. clp .LT. 1.0_r8) THEN
             hst(kd) = hstkd
             qst(kd) = qstkd
             ltl(kd) = ltlkd
             q0u(kd) = q0ukd
             q0d(kd) = q0dkd
             dlb(kd) = dlbkd
             qrb(kd) = qrbkd
          ENDIF
          DO l=kd,kb1
             FCO(L) = FCO(L) - q0u(l) - q0d(l)
             RNN(L) = RNN(L) - q0u(l)*zet(l) - q0d(l)*zet(l+1)
             GMH(L) = GMH(L) - q0u(l)*xi(l)  - q0d(l)*zet(l+1)
             ETA(L) = ZET(L) - ZET(L+1)
             GMS(L) = XI(L)  - XI(L+1)
             Q0U(L) = 0.0_r8
             Q0D(L) = 0.0_r8
          ENDDO
          qw00 = 0.0_r8
          qi00 = 0.0_r8

          !     if (lprnt) print *,' returning to 777 : ii=',ii,' qw00=',qw00,qi00
          !    &,' clp=',clp,' hst(kd)=',hst(kd)

          go to 777
       ELSE
          unsat = .TRUE.
       ENDIF
    ENDIF
    !
    !
    !     ST1 = 0.5 * (HST(KD)  - LTL(KD)*con_FVirt*(QST(KD)-QOS)
    !    &          +  HST(KD1) - LTL(KD1)*con_FVirt*(QST(KD1)-QOL(KD1)))
    !
    ST1 = HST(KD)  - LTL(KD)*con_FVirt*(QST(KD)-QOS)
    ST2 = LTL(KD)  * VTF(KD)
    TEM5 = (QLS + QIS) * eta(kd1)
    ST1  = HALF * (TX1-ETA(KD1)*ST1-ST2*(DET-TEM5))*DLB(KD)
    !
    !     if (lprnt) print *,' st1=',st1,' st2=',st2,' ltl=',ltl(kd)
    !    *,ltl(kd1),' qos=',qos,qol(kd1)

    WFN = WFN + ST1
    AKM = AKM - MIN(ST1,ZERO)   ! Commented on 08/26/02 - does not include top
    !

    BUY(KD) = ST1 / (ETA(KD1)*qrb(kd))
    !
    !     if (lprnt) print *,' wfn=',wfn,' akm=',akm,' st1=',st1
    !    &,' dpneg=',dpneg

    DET = DETP
    HCC = HCCP
    AKM = AKM / WFN


    !***********************************************************************
    !
    !     If only to calculate workfunction save it and return
    !
    IF (WRKFUN) THEN
       IF (WFN .GE. 0.0_r8) WFNC = WFN
       RETURN
    ELSEIF (.NOT. CRTFUN) THEN
       ACR = WFNC
    ENDIF
    !
    !===>  THIRD CHECK BASED ON CLOUD WORKFUNCTION
    !
    CALCUP = .FALSE.

    TEM  =  max(0.005, MIN(CD*200.0, MAX_NEG_BOUY))

    IF (WFN .GT. ACR .AND.  (.NOT. UNSAT)                             &
       !    & .and. dpneg .lt. 100.0_r8  .AND. AKM .LE. TEM) THEN
         & .AND. dpneg .LT. 150.0_r8  .AND. AKM .LE. TEM) THEN
       !    & .and. dpneg .lt. 200.0_r8  .AND. AKM .LE. TEM) THEN
       !
       CALCUP = .TRUE.
    ENDIF

    !     if (lprnt) print *,' calcup=',calcup,' akm=',akm,' tem=',tem
    !    *,' unsat=',unsat,' clp=',clp,' rhc=',rhc,' cd=',cd,' acr=',acr
    !
    !===>  IF NO SOUNDING MEETS THIRD CONDITION, RETURN
    !
    !     if (lprnt .and. kd .eq. 15) stop
    IF (.NOT. CALCUP) RETURN
    !
    ! This is for not LL - 20050601
    IF (ALMIN2 .NE. 0.0_r8) THEN
       IF (ALMIN1 .NE. ALMIN2) ST1 = 1.0_r8 / MAX(ONE_M10,(ALMIN2-ALMIN1))
       IF (ALM .LT. ALMIN2) THEN
          CLP = CLP * MAX(0.0_r8, MIN(1.0_r8,(0.3_r8 + 0.7_r8*(ALM-ALMIN1)*ST1)))
          !          CLP = CLP * max(0.0_r8, min(1.0_r8,(0.2_r8 + 0.8_r8*(ALM-ALMIN1)*ST1)))
          !          CLP = CLP * max(0.0_r8, min(1.0_r8,(0.1_r8 + 0.9_r8*(ALM-ALMIN1)*ST1)))
       ENDIF
    ENDIF
    !
    !     if (lprnt) print *,' clp=',clp
    !
    CLP = CLP * RHC
    dlq = 0.0
    tem = 1.0 / (1.0 + dlq_fac)
    DO l=kd,kb1
       rnn(l) = rns(l)
       dlq(l) = rns(l) * tem * dlq_fac
    ENDDO
    DO L=KBL,K 
       RNN(L) = 0.0_r8 
    ENDDO
    !     if (lprnt) print *,' rnn=',rnn
    !
    !     If downdraft is to be invoked, do preliminary check to see
    !     if enough rain is available and then call DDRFT.
    !
    DDFT = .FALSE.
    IF (DNDRFT) THEN
       !
       TRAIN = 0.0_r8
       IF (CLP .GT. 0.0_r8) THEN
          DO L=KD,KB1
             TRAIN = TRAIN + RNN(L)
          ENDDO
       ENDIF

       PL = (PRL(KD1) + PRL(KD))*HALF
       TEM = PRL(K+1)*(1.0_r8-DPD*0.001_r8)
       IF (TRAIN .GT. 1.0E-4_r8 .AND. PL .LE. TEM) DDFT  = .TRUE.
    ENDIF
    !
    !     if (lprnt) print *,' BEFORE CALLING DDRFT KD=',kd,' DDFT=',DDFT
    !    &,                  ' PL=',PL,' TRAIN=',TRAIN
    !     if (lprnt) print *,' buy=',(buy(l),l=kd,kb1)

    IF (DDFT) THEN
       !
       !     Call Downdraft scheme based on (Cheng and Arakawa, 1997)
       !
       CALL DDRFT( &
                   K              , &!INTEGER      , INTENT(IN   ) :: K
                   KD             , &!INTEGER      , INTENT(IN   ) :: KD
                   TLA            , &!REAL(kind=r8), INTENT(INOUT) :: TLA
                   ALFIND(1:K)    , &!REAL(kind=r8), INTENT(IN   ) :: ALFIND(K)
                   TOL   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: TOL(KD:K)
                   QOL   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: QOL(KD:K)
                   HOL   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: HOL(KD:K)
                   PRL   (1:K+1) , &!REAL(kind=r8), INTENT(IN   ) :: PRL(KD:K+1)
                   QST   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: QST(KD:K)  
                   HST   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: HST(KD:K) 
                   GAM   (1:K+1) , &!REAL(kind=r8), INTENT(IN   ) :: GAM(KD:K+1)
                   GAF   (1:K+1) , &!REAL(kind=r8), INTENT(IN   ) :: GAF(KD:K+1)
                   QRB   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: QRB(KD:K)
                   QRT   (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: QRT(KD:K) 
                   BUY   (1:K+1) , &!REAL(kind=r8), INTENT(INOUT) :: BUY(KD:K+1)
                   KBL            , &!INTEGER,       INTENT(IN   ) :: KBL
                   IDH            , &!INTEGER,       INTENT(IN   ) :: IDH
                   ETA   (1:K+1) , &!REAL(kind=r8), INTENT(IN   ) :: ETA(KD:K+1)
                   RNN   (1:K)   , &!REAL(kind=r8), INTENT(INOUT) :: RNN(KD:K)
                   ETAI  (1:K)   , &!REAL(kind=r8), INTENT(IN   ) :: ETAI(KD:K)
                   ALM            , &!REAL(kind=r8), INTENT(IN   ) :: ALM
                   TRAIN          , &!REAL(kind=r8), INTENT(IN   ) :: TRAIN
                   DDFT           , &!LOGICAL,       INTENT(INOUT) :: DDFT
                   ETD   (1:K+1) , &!REAL(kind=r8), INTENT(INOUT) :: ETD(KD:K+1)
                   HOD   (1:K+1) , &!REAL(kind=r8), INTENT(INOUT) :: HOD(KD:K+1)
                   QOD   (1:K+1) , &!REAL(kind=r8), INTENT(INOUT) :: QOD(KD:K+1)
                   EVP   (1:K)   , &!REAL(kind=r8), INTENT(INOUT) :: EVP(KD:K)
                   DOF            , &!REAL(kind=r8), INTENT(OUT  ) :: DOF
                   CLDFR (1:K)   , &!REAL(kind=r8), INTENT(OUT  ) :: CLDFRD(KD:K)  
                   ETZ   (1:K)   , &!REAL(kind=r8), INTENT(OUT  ) :: WCB(KD:K)
                   GMS   (1:K+1) , &!REAL(kind=r8), INTENT(OUT  ) :: GMS(KD:K+1)
                   GSD   (1:K)   , &!REAL(kind=r8), INTENT(OUT  ) :: GSD(KD:K)
                   GHD   (1:K)     )!REAL(kind=r8), INTENT(OUT  ) :: GHD(KD:K)

    ENDIF
    !
    !  No Downdraft case (including case with no downdraft soln)
    !  ---------------------------------------------------------
    !
    IF (.NOT. DDFT) THEN
       DO L=KD,K+1
          ETD(L) = 0.0_r8
          HOD(L) = 0.0_r8
          QOD(L) = 0.0_r8
       ENDDO
       DO L=KD,K
          EVP(L) = 0.0_r8
          ETZ(L) = 0.0_r8
       ENDDO

    ENDIF
    !     if (lprnt) print *,' hod=',hod
    !     if (lprnt) print *,' etd=',etd
    !
    !
    !===> CALCULATE GAMMAS  i.e. TENDENCIES PER UNIT CLOUD BASE MASSFLUX
    !           Includes downdraft terms!

    avh = 0.0_r8

    !
    !     Fraction of detrained condensate evaporated
    !
    !     tem1 = max(ZERO, min(HALF, (prl(kd)-FOUR_P2)*ONE_M2))
    !     tem1 = max(ZERO, min(HALF, (prl(kd)-300.0_r8)*0.005_r8))
    tem1 = 0.0_r8
    !     tem1 = 1.0_r8
    !     if (kd1 .eq. kbl) tem1 = 0.0_r8
    !
    tem2    = 1.0_r8 - tem1
    TEM = DET * QIL(KD)


    st1 = (HCC+con_hfus*TEM-ETA(KD)*HST(KD)) / (1.0_r8+gam(KD))
    DS  = ETA(KD1) * (HOS- HOL(KD)) - con_hvap *(QOS - QOL(KD))
    DH  = ETA(KD1) * (HOS- HOL(KD))


    GMS(KD) = (DS + st1 - tem1*det*con_hvap -tem*con_hfus) * PRI(KD)
    GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOS + DH)


    !     if (lprnt) print *,' gmhkd=',gmh(kd),' gmskd=',gms(kd)
    !    &,' det=',det,' tem=',tem,' tem1=',tem1,' tem2=',tem2
    !
    !      TENDENCY FOR SUSPENDED ENVIRONMENTAL ICE AND/OR LIQUID WATER
    !
      QLL(KD) = (tem2*(DET-TEM) + ETA(KD1)*(QLS-CLL(KD))                &
     &        + (1.0-QIL(KD))*dlq(kd) - ETA(KD)*QLS ) * PRI(KD)

      QIL(KD) =     (tem2*TEM + ETA(KD1)*(QIS-CIL(KD))                  &
     &        + QIL(KD)*dlq(kd) - ETA(KD)*QIS ) * PRI(KD)
    !
    GHD(KD) = 0.0_r8
    GSD(KD) = 0.0_r8
    !
    DO L=KD1,K
       ST1 = ONE - ALFINT(L,1)
       ST2 = ONE - ALFINT(L,2)
       ST3 = ONE - ALFINT(L,3)
       ST4 = ONE - ALFINT(L,4)
       ST5 = ONE - ALFIND(L)
       HB       = ALFINT(L,1)*HOL(L-1) + ST1*HOL(L)
       QB       = ALFINT(L,2)*QOL(L-1) + ST2*QOL(L)

       TEM      = ALFINT(L,4)*CIL(L-1) + ST4*CIL(L)
       TEM2     = ALFINT(L,3)*CLL(L-1) + ST3*CLL(L)

       TEM1     = ETA(L) * (TEM - CIL(L))
       TEM3     = ETA(L) * (TEM2 - CLL(L))

       HBD      = ALFIND(L)*HOL(L-1) + ST5*HOL(L)
       QBD      = ALFIND(L)*QOL(L-1) + ST5*QOL(L)

       TEM5     = ETD(L) * (HOD(L) - HBD)
       TEM6     = ETD(L) * (QOD(L) - QBD)
       !
       DH       = ETA(L) * (HB - HOL(L)) + TEM5
       DS       = DH - con_hvap  * (ETA(L) * (QB - QOL(L)) + TEM6)

       GMH(L)   = DH * PRI(L)
       GMS(L)   = DS * PRI(L)

       !     if (lprnt) print *,' gmh=',gmh(l),' gms=',gms(l)
       !    &,' dh=',dh,' ds=',ds,' qb=',qb,' qol=',qol(l),' eta=',eta(l)
       !    &,' hb=',hb,' hol=',hol(l),' l=',l,' hod=',hod(l)
       !    &,' etd=',etd(l),' qod=',qod(l),' tem5=',tem5,' tem6=',tem6
       !
       GHD(L)   = TEM5 * PRI(L)
       GSD(L)   = (TEM5 - con_hvap  * TEM6) * PRI(L)
       !
       QIL(L)   = TEM1 * PRI(L)
       QLL(L)   = TEM3 * PRI(L)

       TEM1     = ETA(L) * (CIL(L-1) - TEM)
       TEM3     = ETA(L) * (CLL(L-1) - TEM2)

       DH       = ETA(L) * (HOL(L-1) - HB) - TEM5
       DS       = DH - con_hvap  * ETA(L) * (QOL(L-1) - QB)                &
            &                 + con_hvap  * (TEM6 - EVP(L-1))

       GMH(L-1) = GMH(L-1) + DH * PRI(L-1)
       GMS(L-1) = GMS(L-1) + DS * PRI(L-1)
       !
       !     if (lprnt) print *,' gmh1=',gmh(l-1),' gms1=',gms(l-1)
       !    &,' dh=',dh,' ds=',ds,' qb=',qb,' qol=',qol(l-1)
       !    &,' hb=',hb,' hol=',hol(l-1),' evp=',evp(l-1)
       !
       GHD(L-1) = GHD(L-1) - TEM5 * PRI(L-1)
       GSD(L-1) = GSD(L-1) - (TEM5-con_hvap *(TEM6-EVP(L-1))) * PRI(L-1)

       QIL(L-1) = QIL(L-1) + TEM1 * PRI(L-1)
       QLL(L-1) = QLL(L-1) + TEM3 * PRI(L-1)

       !     if (lprnt) print *,' gmh=',gmh(l),' gms=',gms(l)
       !    &,' dh=',dh,' ds=',ds,' qb=',qb,' qol=',qol(l),' eta=',eta(l)
       !    &,' hb=',hb,' hol=',hol(l),' l=',l
       !
       avh = avh + gmh(l-1)*(prs(l)-prs(l-1))

    ENDDO

    HBD  = HOL(K)
    QBD  = QOL(K)
    TEM5 =  ETD(K+1) * (HOD(K+1) - HBD)
    TEM6 =  ETD(K+1) * (QOD(K+1) - QBD)
    DH   = - TEM5
    DS   = DH  + con_hvap  * TEM6
    TEM1 = DH * PRI(K)
    TEM2 = (DS - con_hvap  * EVP(K)) * PRI(K)
    !!    TEM2 = - con_hvap  * EVP(K) * PRI(K)
    GMH(K) = GMH(K) + TEM1
    GMS(K) = GMS(K) + TEM2
    GHD(K) = GHD(K) + TEM1
    GSD(K) = GSD(K) + TEM2

    !     if (lprnt) print *,' gmhk=',gmh(k),' gmsk=',gms(k)
    !    &,' tem1=',tem1,' tem2=',tem2,' dh=',dh,' ds=',ds
    !
    avh = avh + gmh(K)*(prs(KP1)-prs(K))
    !
    tem4   = - GRAVFAC * pris
    TX1    = DH * tem4
    TX2    = DS * tem4
    !
    DO L=KBL,K
       GMH(L) = GMH(L) + TX1
       GMS(L) = GMS(L) + TX2
       GHD(L) = GHD(L) + TX1
       GSD(L) = GSD(L) + TX2
       !
       avh = avh + tx1*(prs(l+1)-prs(l))
    ENDDO

    !
    !     if (lprnt) then
    !        print *,' gmh=',gmh
    !        print *,' gms=',gms(KD:K)
    !     endif
    !
    !***********************************************************************
    !***********************************************************************

    !===>  KERNEL (AKM) CALCULATION BEGINS

    !===>  MODIFY SOUNDING WITH UNIT MASS FLUX
    !
    !     TESTMB = 0.01_r8

    DO L=KD,K

       TEM1   = GMH(L)
       TEM2   = GMS(L)
       HOL(L) = HOL(L) +  TEM1*TESTMB
       QOL(L) = QOL(L) + (TEM1-TEM2)  * (TESTMB/con_hvap )
       HST(L) = HST(L) +  TEM2*(ONE+GAM(L))*TESTMB
       QST(L) = QST(L) +  TEM2*GAM(L)*(TESTMB/con_hvap )
       CLL(L) = CLL(L) + QLL(L) * TESTMB
       CIL(L) = CIL(L) + QIL(L) * TESTMB
    ENDDO
    !

    IF (alm .GT. 0.0_r8) THEN
       HOS = HOS + GMH(KD)  * TESTMB
       QOS = QOS + (GMH(KD)-GMS(KD)) * (TESTMB/con_hvap )
       QLS     = QLS + QLL(KD) * TESTMB
       QIS     = QIS + QIL(KD) * TESTMB
    ELSE
       st2 = 1.0_r8 - st1s
       HOS = HOS + (st1s*GMH(KD)+st2*GMH(KD1))  * TESTMB
       QOS = QOS + (st1s * (GMH(KD)-GMS(KD))                           &
            &            +  st2  * (GMH(KD1)-GMS(KD1))) * (TESTMB/con_hvap )
       HST(kd) = HST(kd) + (st1s*GMS(kd)*(ONE+GAM(kd))                 &
            &                    +  st2*gms(kd1)*(ONE+GAM(kd1))) * TESTMB
       QST(kd) = QST(kd) + (st1s*GMS(kd)*GAM(kd)                       &
            &                    +  st2*gms(kd1)*gam(kd1)) * (TESTMB/con_hvap )

       QLS     = QLS + (st1s*QLL(KD)+st2*QLL(KD1)) * TESTMB
       QIS     = QIS + (st1s*QIL(KD)+st2*QIL(KD1)) * TESTMB
    ENDIF

    !
      TEM = PRL(Kmaxp1) - PRL(Kmax)
      HBL = HOL(Kmax) * TEM
      QBL = QOL(Kmax) * TEM
      QLB = CLL(Kmax) * TEM
      QIB = CIL(Kmax) * TEM
      DO L=KmaxM1,KBL,-1
        TEM = PRL(L+1) - PRL(L)
        HBL = HBL + HOL(L) * TEM
        QBL = QBL + QOL(L) * TEM
        QLB = QLB + CLL(L) * TEM
        QIB = QIB + CIL(L) * TEM
      ENDDO
    HBL = HBL * PRIS
    QBL = QBL * PRIS
    QLB = QLB * PRIS
    QIB = QIB * PRIS

!     if (ctei .and. sgcs(ipt,l) > 0.65) then
!        hbl = hbl * hpert_fac
!        qbl = qbl * hpert_fac
!     endif

    !     if (lprnt) print *,' hbla=',hbl,' qbla=',qbl

    !***********************************************************************

    !===>  CLOUD WORKFUNCTION FOR MODIFIED SOUNDING, THEN KERNEL (AKM)
    !
    AKM = ZERO
    TX1 = ZERO
    QTL = QST(KB1) - GAF(KB1)*HST(KB1)
    QTV = QBL
    HCC = HBL
    TX2 = HCC
    TX4 = (con_hfus*0.5_r8)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(KB1))*TCRF))
    !
    qtv   = qbl
    tx1   = qib + qlb
    !

    DO L=KB1,KD1,-1
       DEL_ETA = ETA(L) - ETA(L+1)
       HCCP = HCC + DEL_ETA*HOL(L)
       !
       QTLP = QST(L-1) - GAF(L-1)*HST(L-1)
       QTVP = 0.5_r8 * ((QTLP+QTL)*ETA(L) +(GAF(L)+GAF(L-1))*HCCP)

       DETP = (BKC(L)*TX1 - (QTVP-QTV)                                &
            &        +  DEL_ETA*(QOL(L)+CLL(L)+CIL(L))                         &
            &        +  ETA(L)*Q0U(L) + ETA(L+1)*Q0D(L)) * AKC(L)
       IF (DETP .LE. ZERO) UNSAT = .TRUE.

       ST1  = HST(L) - LTL(L)*con_FVirt*(QST(L)-QOL(L))

       TEM2 = (con_hfus*0.5_r8)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(L-1))*TCRF))
       TEM1 = HCCP + DETP * (TEM2+TX4)

       ST2  = LTL(L) * VTF(L)
       TEM5 = CLL(L) + CIL(L)
       AKM  = AKM +                                                   &
            &     (  (TX2  -ETA(L+1)*ST1-ST2*(TX1-TEM5*eta(l+1))) * DLB(L)     &
            &      + (TEM1 -ETA(L  )*ST1-ST2*(DETP-TEM5*eta(l)))  * DLT(L) )
       !
       HCC  = HCCP
       TX1  = DETP
       TX2  = TEM1
       QTL  = QTLP
       QTV  = QTVP
       TX4  = TEM2
    ENDDO
    !
    IF (unsat) RETURN
    !
    !  Eventhough we ignore the change in lambda, we still assume
    !  that the cLoud-top contribution is zero; as though we still
    !  had non-bouyancy there.
    !
    !
    ST1 = HST(KD)  - LTL(KD)*con_FVirt*(QST(KD)-QOS)
    ST2 = LTL(KD)  * VTF(KD)
    TEM5 = (QLS + QIS) * eta(kd1)
    AKM  = AKM + HALF * (TX2-ETA(KD1)*ST1-ST2*(TX1-TEM5)) * DLB(KD)
    !
    AKM = (AKM - WFN) * (ONE/TESTMB)


    !***********************************************************************

    !===>   MASS FLUX

    tem2 = rel_fac
    !
    AMB = - (WFN-ACR) / AKM
    !
    !IF(lprnt) PRINT *,' wfn=',wfn,' acr=',acr,' akm=',akm             &
    !     &,' amb=',amb,' KD=',kd,' cldfrd=',cldfrd,' tem2=',tem2            &
    !     &,' rel_fac=',rel_fac,' prskd=',prs(kd)

    !===>   RELAXATION AND CLIPPING FACTORS
    !
    AMB = AMB * CLP * tem2

!!!   if (DDFT) AMB = MIN(AMB, ONE/CLDFRD)

    !===>   SUB-CLOUD LAYER DEPTH LIMIT ON MASS FLUX

      AMBMAX = (PRL(KMAXP1)-PRL(KBL))*(FRACBL*GRAVCON)
    AMB    = MAX(MIN(AMB, AMBMAX),ZERO)


    !     if(lprnt) print *,' AMB=',amb,' clp=',clp,' ambmax=',ambmax
    !***********************************************************************
    !*************************RESULTS***************************************
    !***********************************************************************

    !===>  PRECIPITATION AND CLW DETRAINMENT
    !
    avt = 0.0_r8
    avq = 0.0_r8
    avr = dof

    !
    DSFC = DSFC + AMB * ETD(K) * (1.0_r8/DT)
    !
    !     DO L=KBL,KD,-1
    DO L=K,KD,-1
       PCU(L) = PCU(L) + AMB*RNN(L)      !  (A40)
       avr = avr + rnn(l)
       !     if(lprnt) print *,' avr=',avr,' rnn=',rnn(l),' l=',l
    ENDDO
    !
    !===> TEMPARATURE AND Q CHANGE AND CLOUD MASS FLUX DUE TO CLOUD TYPE KD
    !
    TX1  = AMB * (ONE/con_cp)
    TX2  = AMB * (ONE/con_hvap )
    DO L=KD,K
       ST1    = GMS(L)*TX1
       TOI(L) = TOI(L) + ST1
       TCU(L) = TCU(L) + ST1
       TCD(L) = TCD(L) + GSD(L) * TX1
       !
       !       st1 = st1 - (con_hfus/con_cp) * QIL(L) * AMB
       st1 = st1 - (con_hvap /con_cp) * (QIL(L) + QLL(L)) * AMB

       avt = avt + st1 * (prs(l+1)-prs(l))

       FLX(L)  = FLX(L)  + ETA(L)*AMB
       FLXD(L) = FLXD(L) + ETD(L)*AMB
       !
       QII(L)  = QII(L) + QIL(L) * AMB
       TEM     = 0.0_r8

       QLI(L)  = QLI(L) + QLL(L) * AMB + TEM

       ST1     = (GMH(L)-GMS(L)) * TX2

       QOI(L)  = QOI(L) + ST1
       QCU(L)  = QCU(L) + ST1
       QCD(L)  = QCD(L) + (GHD(L)-GSD(L)) * TX2
       !
       avq = avq + (st1+(QLL(L)+QIL(L))*amb) * (prs(l+1)-prs(l))
       !       avq = avq + st1 * (prs(l+1)-prs(l))
       !       avr = avr + (QLL(L) + QIL(L)*(1+con_hfus/con_hvap ))
       !       avr = avr + (QLL(L) + QIL(L))
       !    *                  * (prs(l+1)-prs(l)) * gravcon

       !     if(lprnt) print *,' avr=',avr,' qll=',qll(l),' l=',l
       !    &,' qil=',qil(l)

    ENDDO
    avr = avr * amb
    !
    !      Correction for negative condensate!
    !     if (advcld) then
    !       do l=kd,k
    !         if (qli(l) .lt. 0.0_r8) then
    !           qoi(l) = qoi(l) + qli(l)
    !           toi(l) = toi(l) - (con_hvap /con_cp) * qli(l)
    !           qli(l) = 0.0_r8
    !         endif
    !         if (qii(l) .lt. 0.0_r8) then
    !           qoi(l) = qoi(l) + qii(l)
    !           toi(l) = toi(l) - ((con_hvap +con_hfus)/con_cp) * qii(l)
    !           qii(l) = 0.0_r8
    !         endif
    !       enddo
    !     endif

    !
    !
    !     if (lprnt) then
    !       print *,' For KD=',KD
    !       avt = avt * con_cp * 100.0_r8*86400.0_r8 / (con_hvap *DT*con_g)
    !       avq = avq *  100.0_r8*86400.0_r8 / (DT*con_g)
    !       avr = avr * 86400.0_r8 / DT
    !       print *,' avt=',avt,' avq=',avq,' avr=',avr,' avh='
    !    *   ,avh,' alm=',alm,' DDFT=',DDFT,' KD=',KD
    !    &,' TOIK-',toi(k),' TOIK-1=',toi(k-1),' TOIK-2=',toi(k-2)
    !        if (kd .eq. 12 .and. .not. ddft) stop
    !       if (avh .gt. 0.1_r8 .or. abs(avt+avq) .gt. 1.0e-5_r8 .or.
    !    &      abs(avt-avr) .gt. 1.0e-5_r8 .or. abs(avr+avq) .gt. 1.0e-5_r8) stop
    !
    !     if (lprnt) then
    !       print *,' For KD=',KD
    !       print *,' TCU=',(tcu(l),l=kd,k)
    !       print *,' QCU=',(Qcu(l),l=kd,k)
    !     endif
    !
    TX1 = 0.0_r8
    TX2 = 0.0_r8
    !
    !     REEVAPORATION OF FALLING CONVECTIVE RAIN
    !
    IF (REVAP) THEN

       !
       tem = 0.0_r8
       DO l=kd,kbl
          !        tem = tem + pcu(l)
          IF (L .LT. IDH .OR. (.NOT. DDFT)) THEN
             tem = tem + amb * rnn(l)
          ENDIF
       ENDDO
       tem = tem + amb * dof
       tem = tem * (3600.0_r8/dt)
!!!!   tem1 = max(1.0, min(100.0,sqrt((5.0E10/max(garea,one)))))
!      tem1 = max(1.0, min(100.0,(7.5E10/max(garea,one))))
!      tem1 = max(1.0, min(100.0,(5.0E10/max(garea,one))))
!      tem1 = max(1.0, min(100.0,(4.0E10/max(garea,one))))
!!     tem1 = sqrt(max(1.0, min(100.0,(4.0E10/max(garea,one))))) ! 20100902
!       tem1 = MAX(1.0_r8, MIN(100.0_r8,SQRT((5.0E10_r8/MAX(garea,one)))))
       tem1 = sqrt(max(1.0_r8, min(100.0_r8,(6.25E10_r8/max(garea,one))))) ! 20110530

       !      if (lprnt) print *,' clfr0=',clf(tem),' tem=',tem,' tem1=',tem1

       !clfrac = MAX(ZERO, MIN(ONE, rknob*clf(tem)*tem1))
!      clfrac = max(ZERO, min(0.25, rknob*clf(tem)*tem1))
       clfrac = max(ZERO, min(half, rknob*clf(tem)*tem1))

       !      if (lprnt) print *,' cldfrd=',cldfrd,' amb=',amb
       !    &,' clfrac=',clfrac

       !      TX3    = AMB*ETA(KD)*PRI(KD)
       !
       !cnt   DO L=KD,K
       DO L=KD,KBL         ! Testing on 20070926
          IF (L .GE. IDH .AND. DDFT) THEN
             TX2    = TX2 + AMB * RNN(L)
             CLDFRD = MIN(AMB*CLDFR(L), clfrac)
          ELSE
             TX1 = TX1 + AMB * RNN(L)
          ENDIF
          tx4 = zfac * phil(l)
          tx4 = (one - tx4 * (one - half*tx4)) * afc
          !
          IF (TX1 .GT. 0.0_r8 .OR. TX2 .GT. 0.0_r8) THEN
             TEQ     = TOI(L)
             QEQ     = QOI(L)
             PL      = 0.5_r8 * (PRL(L+1)+PRL(L))

             ST1     = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
             ST2     = ST1*ELFOCP + (1.0_r8-ST1)*ELOCP

             CALL QSATCN ( TEQ,PL,QSTEQ,DQDT)
!            CALL QSATCN ( TEQ,PL,QSTEQ,DQDT,.false.)
             !
             DELTAQ = 0.5_r8 * (QSTEQ*rhc_ls(l)-QEQ) / (1.0_r8+ST2*DQDT)
             !
             QEQ    = QEQ + DELTAQ
             TEQ    = TEQ - DELTAQ*ST2
             !
             TEM1   = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
             TEM2   = TEM1*ELFOCP + (1.0_r8-TEM1)*ELOCP

             CALL QSATCN ( TEQ,PL,QSTEQ,DQDT)
!            CALL QSATCN ( TEQ,PL,QSTEQ,DQDT,.false.)             !

             DELTAQ = (QSTEQ*rhc_ls(l)-QEQ) / (1.0_r8+TEM2*DQDT)
             !
             QEQ    = QEQ + DELTAQ
             TEQ    = TEQ - DELTAQ*TEM2

             IF (QEQ .GT. QOI(L)) THEN
                POTEVAP = (QEQ-QOI(L))*(PRL(L+1)-PRL(L))*GRAVCON

                !           TEM3    = SQRT(PL*0.001_r8)
                tem4    = 0.0_r8
                IF (tx1 .GT. 0.0_r8)                                           &
                     &      TEM4    = POTEVAP * (1.0_r8 - EXP( tx4*TX1**0.57777778_r8 ) )
                !    &      TEM4    = POTEVAP * (1.0_r8 - EXP( AFC*tx4*SQRT(TX1) ) )
                ACTEVAP = MIN(TX1, TEM4*CLFRAC)

                !     if(lprnt) print *,' L=',L,' actevap=',actevap,' tem4=',tem4,
                !    &' clfrac='
                !    &,clfrac,' potevap=',potevap,'efac=',AFC*SQRT(TX1*TEM3)
                !    &,' tx1=',tx1

                IF (tx1 .LT. rainmin*dt) actevap = MIN(tx1, potevap)
                !
                tem4    = 0.0_r8
                IF (tx2 .GT. 0.0_r8)                                           &
                     &      TEM4    = POTEVAP * (1.0_r8 - EXP( tx4*TX2**0.57777778_r8 ) )
!    &      TEM4    = POTEVAP * (1. - EXP( AFC*tx4*SQRT(TX2) ) )
                TEM4    = MIN(MIN(TX2, TEM4*CLDFRD), potevap-actevap)
                IF (tx2 .LT. rainmin*dt) tem4 = MIN(tx2, potevap-actevap)
                !
                TX1     = TX1 - ACTEVAP
                TX2     = TX2 - TEM4
                ST1     = (ACTEVAP+TEM4) * PRI(L)
                QOI(L)  = QOI(L) + ST1
                QCU(L)  = QCU(L) + ST1
                !

                ST1     = ST1 * ELOCP
                TOI(L)  = TOI(L) - ST1 
                TCU(L)  = TCU(L) - ST1
             ENDIF
          ENDIF
       ENDDO
       !
       CUP = CUP + TX1 + TX2 + DOF * AMB
    ELSE
       DO L=KD,K
          TX1 = TX1 + AMB * RNN(L)
       ENDDO
       CUP = CUP + TX1 + DOF * AMB
    ENDIF

    !     if (lprnt) print *,' tx1=',tx1,' tx2=',tx2,' dof=',dof
    !    &,' cup=',cup*86400/dt,' amb=',amb
    !    &,' amb=',amb,' cup=',cup,' clfrac=',clfrac,' cldfrd=',cldfrd
    !    &,' ddft=',ddft,' kd=',kd,' kbl=',kbl,' k=',k
!
!    Convective transport (mixing) of passive tracers
!
    !
    !    MIXING OF PASSIVE TRACERS
    !
      if (m > 0) then
        do l=kd,k-1
          if (etz(l) .ne. zero) etzi(l) = one / etz(l)
        enddo

    DO N=1,M

       DO L=KD,K
          HOL(L) = ROI(L,N)
       ENDDO
       !
       HCC     = RBL(N)
       HOD(KD) = HOL(KD)
       !      Compute downdraft properties for the tracer
       DO L=KD1,K
          ST1 = ONE - ALFIND(L)
          HB  = ALFIND(L)  * HOL(L-1) + ST1 * HOL(L)
          IF (ETZ(L-1) .NE. 0.0_r8) THEN
             DEL_ETA = ETD(L) - ETD(L-1)
            ! TEM     = 1.0_r8 / ETZ(L-1)
              TEM = ETZI(L-1)
             IF (DEL_ETA .GT. 0.0_r8) THEN
                HOD(L) = (ETD(L-1)*(HOD(L-1)-HOL(L-1))                     &
                       +  ETD(L)  *(HOL(L-1)-HB) +  ETZ(L-1)*HB) * TEM
             ELSE
                HOD(L) = (ETD(L-1)*(HOD(L-1)-HB) + ETZ(L-1)*HB) * TEM
             ENDIF
          ELSE
             HOD(L) = HB
          ENDIF
       ENDDO

       DO L=KB1,KD,-1
          HCC = HCC + (ETA(L)-ETA(L+1))*HOL(L)
       ENDDO
!
!         Scavenging -- fscav   - fraction scavenged [km-1]
!                       delz    - distance from the entrainment to detrainment layer [km]
!                       fnoscav - the fraction not scavenged
!                                 following Liu et al. [JGR,2001] Eq 1

       if (FSCAV_(N) > 0.0_r8) then
         DELZKM = ( PHIL(KD) - PHIH(KD1) ) *(onebg*0.001_r8)
         FNOSCAV = exp(- FSCAV_(N) * DELZKM)
       else
         FNOSCAV = 1.0_r8
       endif
       !
       GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOL(KD)) * trcfac(n,k)      &
                                                    * FNOSCAV
       DO L=KD1,K
         if (FSCAV_(N) > 0.0) then
             DELZKM = ( PHIL(KD) - PHIH(L+1) ) *(onebg*0.001_r8)
             FNOSCAV = exp(- FSCAV_(N) * DELZKM)
          end if
          ST1 = ONE - ALFINT(L,N+4)
          ST2 = ONE - ALFIND(L)
          HB       = ALFINT(L,N+4) * HOL(L-1) + ST1 * HOL(L)
          HBD      = ALFIND(L) * HOL(L-1) + ST2 * HOL(L)
          TEM5     = ETD(L)    * (HOD(L) - HBD)
          DH       = ETA(L)    * (HB - HOL(L))   * FNOSCAV + TEM5
          GMH(L  ) = DH * PRI(L) * trcfac(n,l)
          DH       = ETA(L)    * (HOL(L-1) - HB) * FNOSCAV - TEM5
          GMH(L-1) = GMH(L-1)  + DH * PRI(L-1) * trcfac(n,l)
       ENDDO
       !
       DO L=KD,K
          ST1      = GMH(L)*AMB
          ROI(L,N) = HOL(L)   + ST1
          RCU(L,N) = RCU(L,N) + ST1
       ENDDO
    ENDDO                             ! Tracer loop M
      endif

    !     if (lprnt) print *,' toio=',toi
    !     if (lprnt) print *,' qoio=',qoi
    RETURN
  END SUBROUTINE CLOUD


  !-----------------------------------------------------------------------------------------

!  SUBROUTINE DDRFT(                                                 &
!       &                  K, KD                                           &
!       &,                 TLA, ALFIND                                     &
!       &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF&
!!       &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF, HBL, QBL&
!       &,                 QRB, QRT, BUY, KBL, IDH, ETA, RNN, ETAI         &
!       &,                 ALM, TRAIN, DDFT                           &
!!       &,                 ALM, WFN, TRAIN, DDFT                           &
!       &,                 ETD, HOD, QOD, EVP, DOF, CLDFRD, WCB            &
!       &,                 GMS, GSD, GHD)                   


  SUBROUTINE DDRFT( &
                   K            , &!INTEGER, INTENT(IN   ) :: K
                   KD           , &!INTEGER, INTENT(IN   ) :: KD
                   TLA          , &!REAL(kind=r8), INTENT(INOUT) :: TLA
                   ALFIND       , &!REAL(kind=r8), INTENT(IN   ) :: ALFIND(K)
                   TOL          , &!REAL(kind=r8), INTENT(IN   ) :: TOL(KD:K)
                   QOL          , &!REAL(kind=r8), INTENT(IN   ) :: QOL(KD:K)
                   HOL          , &!REAL(kind=r8), INTENT(IN   ) :: HOL(KD:K)
                   PRL          , &!REAL(kind=r8), INTENT(IN   ) :: PRL(KD:K+1)
                   QST          , &!REAL(kind=r8), INTENT(IN   ) :: QST(KD:K)  
                   HST          , &!REAL(kind=r8), INTENT(IN   ) :: HST(KD:K) 
                   GAM          , &!REAL(kind=r8), INTENT(IN   ) :: GAM(KD:K+1)
                   GAF          , &!REAL(kind=r8), INTENT(IN   ) :: GAF(KD:K+1)
                   QRB          , &!REAL(kind=r8), INTENT(IN   ) :: QRB(KD:K)
                   QRT          , &!REAL(kind=r8), INTENT(IN   ) :: QRT(KD:K) 
                   BUY          , &!REAL(kind=r8), INTENT(INOUT) :: BUY(KD:K+1)
                   KBL          , &!INTEGER,INTENT(IN   ) :: KBL
                   IDH          , &!INTEGER,INTENT(IN   ) :: IDH
                   ETA          , &!REAL(kind=r8), INTENT(IN   ) :: ETA(KD:K+1)
                   RNN          , &!REAL(kind=r8), INTENT(INOUT) :: RNN(KD:K)
                   ETAI         , &!REAL(kind=r8), INTENT(IN   ) :: ETAI(KD:K)
                   ALM          , &!REAL(kind=r8), INTENT(IN   ) :: ALM
                   TRAIN        , &!REAL(kind=r8), INTENT(IN   ) :: TRAIN
                   DDFT         , &!LOGICAL,INTENT(INOUT) :: DDFT
                   ETD          , &!REAL(kind=r8), INTENT(INOUT) :: ETD(KD:K+1)
                   HOD          , &!REAL(kind=r8), INTENT(INOUT) :: HOD(KD:K+1)
                   QOD          , &!REAL(kind=r8), INTENT(INOUT) :: QOD(KD:K+1)
                   EVP          , &!REAL(kind=r8), INTENT(INOUT) :: EVP(KD:K)
                   DOF          , &!REAL(kind=r8), INTENT(OUT  ) :: DOF
                   CLDFRD       , &!REAL(kind=r8), INTENT(OUT  ) :: CLDFRD(KD:K)  
                   WCB          , &!REAL(kind=r8), INTENT(OUT  ) :: WCB(KD:K)
                   GMS          , &!REAL(kind=r8), INTENT(OUT  ) :: GMS(KD:K+1)
                   GSD          , &!REAL(kind=r8), INTENT(OUT  ) :: GSD(KD:K)
                   GHD            )!REAL(kind=r8), INTENT(OUT  ) :: GHD(KD:K)



    !
    !***********************************************************************
    !******************** Cumulus Downdraft Subroutine *********************
    !****************** Based on Cheng and Arakawa (1997)  ****** **********
    !************************ SUBROUTINE DDRFT  ****************************
    !*************************  October 2004  ******************************
    !***********************************************************************
    !***********************************************************************
    !************* Shrinivas.Moorthi@noaa.gov (301) 763 8000(X7233) ********
    !***********************************************************************
    !***********************************************************************
    !23456789012345678901234567890123456789012345678901234567890123456789012
    !
    !===>  TOL(K)     INPUT   TEMPERATURE            KELVIN
    !===>  QOL(K)     INPUT   SPECIFIC HUMIDITY      NON-DIMENSIONAL

    !===>  PRL(K+1)   INPUT   PRESSURE @ EDGES       MB

    !===>  K     INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
    !===>  KD    INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
    !     
    !      use module_ras
    IMPLICIT NONE
    !
    !  INPUT ARGUMENTS
    !
    INTEGER      , INTENT(IN   ) :: K
    INTEGER      , INTENT(IN   ) :: KD
    REAL(kind=r8), INTENT(INOUT) :: TLA
    REAL(kind=r8), INTENT(IN   ) :: ALFIND(K)
    REAL(kind=r8), INTENT(IN   ) :: TOL(1:K)
    REAL(kind=r8), INTENT(IN   ) :: QOL(1:K)
    REAL(kind=r8), INTENT(IN   ) :: HOL(1:K)
    REAL(kind=r8), INTENT(IN   ) :: PRL(1:K+1)
    REAL(kind=r8), INTENT(IN   ) :: QST(1:K)  
    REAL(kind=r8), INTENT(IN   ) :: HST(1:K) 
    REAL(kind=r8), INTENT(IN   ) :: GAM(1:K+1)
    REAL(kind=r8), INTENT(IN   ) :: GAF(1:K+1)
!    REAL(kind=r8), INTENT(IN   ) :: HBL  ! not used
!    REAL(kind=r8), INTENT(IN   ) :: QBL ! not used
    REAL(kind=r8), INTENT(IN   ) :: QRB(1:K)
    REAL(kind=r8), INTENT(IN   ) :: QRT(1:K) 
    REAL(kind=r8), INTENT(INOUT) :: BUY(1:K+1)
    INTEGER,       INTENT(IN   ) :: KBL
    INTEGER,       INTENT(IN   ) :: IDH
    REAL(kind=r8), INTENT(IN   ) :: ETA(1:K+1)
    REAL(kind=r8), INTENT(INOUT) :: RNN(1:K)
    REAL(kind=r8), INTENT(IN   ) :: ETAI(1:K)
    REAL(kind=r8), INTENT(IN   ) :: ALM
    REAL(kind=r8), INTENT(IN   ) :: TRAIN
    LOGICAL,       INTENT(INOUT) :: DDFT
    REAL(kind=r8), INTENT(INOUT) :: ETD(1:K+1)
    REAL(kind=r8), INTENT(INOUT) :: HOD(1:K+1)
    REAL(kind=r8), INTENT(INOUT) :: QOD(1:K+1)
    REAL(kind=r8), INTENT(INOUT) :: EVP(1:K)
    REAL(kind=r8), INTENT(OUT  ) :: DOF
    REAL(kind=r8), INTENT(OUT  ) :: CLDFRD(1:K)  
    REAL(kind=r8), INTENT(OUT  ) :: WCB(1:K)
    REAL(kind=r8), INTENT(OUT  ) :: GMS(1:K+1)
    REAL(kind=r8), INTENT(OUT  ) :: GSD(1:K)
    REAL(kind=r8), INTENT(OUT  ) :: GHD(1:K)


    LOGICAL :: SKPDD, SKPUP
    INTEGER :: KB1

    REAL(kind=r8) :: RNS(1:K) 
    !
    !REAL(kind=r8) :: PRIS 
    !
    !     TEMPORARY WORK SPACE
    !
    REAL(kind=r8) :: TX1
    REAL(kind=r8) :: TX2
    REAL(kind=r8) :: TX3
    REAL(kind=r8) :: TX4
    REAL(kind=r8) :: TX5
    REAL(kind=r8) :: TX6
    !REAL(kind=r8) :: TX7
    REAL(kind=r8) :: TX8
    REAL(kind=r8) :: TX9
    LOGICAL       :: UNSAT

    !REAL(kind=r8) :: TL
    !REAL(kind=r8) :: PL
    !REAL(kind=r8) :: QL
    !REAL(kind=r8) :: QS
    !REAL(kind=r8) :: DQS
    REAL(kind=r8) :: ST1
    !REAL(kind=r8) :: HB
    !REAL(kind=r8) :: QB
    !REAL(kind=r8) :: TB  
    REAL(kind=r8) :: QQQ
    REAL(kind=r8) :: PICON
    REAL(kind=r8) :: DEL_ETA     
    REAL(kind=r8) :: TEM
    REAL(kind=r8) :: TEM1
    REAL(kind=r8) :: TEM2
    REAL(kind=r8) :: TEM3
    REAL(kind=r8) :: TEM4
    REAL(kind=r8) :: ST2   
    REAL(kind=r8) :: ERRH
    REAL(kind=r8) :: ERRW
    REAL(kind=r8) :: ERRE
    !REAL(kind=r8) :: TEM5  
    REAL(kind=r8) :: TEM6
    !REAL(kind=r8) :: HBD
    !REAL(kind=r8) :: QBD

    INTEGER  L,  N,  KD1, II                                     &
         &,       KP1,  KM1, KTEM, KK, KK1, LM1, LL, LP1                 &
         &,        ntla

    !
    INTEGER, PARAMETER :: NUMTLA=2
    !     integer, parameter :: NUMTLA=4
    REAL(KINd=r8),PARAMETER :: ERRMIN=0.0001_r8, ERRMI2=0.1_r8*ERRMIN
    !     parameter (ERRMIN=0.00001_r8, ERRMI2=0.1_r8*ERRMIN)
    !
    REAL(kind=r8) :: STLA
    REAL(kind=r8) :: CTL2
    REAL(kind=r8) :: CTL3
    !REAL(kind=r8) :: CTLA
    !REAL(kind=r8) :: VTRM
    REAL(kind=r8) :: VTPEXP 
    REAL(kind=r8) :: WCMIN
    REAL(kind=r8) :: WCBASE
    REAL(kind=r8) :: F2
    REAL(kind=r8) :: QRAF
    REAL(kind=r8) :: QRBF
    REAL(kind=r8) :: CMPOR 
    REAL(kind=r8) :: del_tla
    !REAL(kind=r8) :: sialf
    !
    REAL(kind=r8),PARAMETER :: ONPG=1.0_r8+0.5_r8, GMF=1.0_r8/ONPG, RPART=0.0_r8
    !     parameter (ONPG=1.0_r8+0.5_r8, GMF=1.0_r8/ONPG, RPART=1.0_r8)
    !     parameter (ONPG=1.0_r8+0.5_r8, GMF=1.0_r8/ONPG, RPART=0.5_r8)
    !     PARAMETER (AA1=1.0_r8, BB1=1.5_r8, CC1=1.1_r8, DD1=0.85_r8, F3=CC1, F5=2.5_r8)
    !     PARAMETER (AA1=2.0_r8, BB1=1.5_r8, CC1=1.1_r8, DD1=0.85_r8, F3=CC1, F5=2.5_r8)
    REAL(kind=r8),PARAMETER ::  AA1=1.0_r8, BB1=1.0_r8, CC1=1.0_r8, DD1=1.0_r8, F3=CC1,  F5=1.0_r8
    REAL(kind=r8),PARAMETER ::  QRMIN=1.0E-6_r8, WC2MIN=0.01_r8, GMF1=GMF/AA1, GMF5=GMF/F5
    !     parameter (QRMIN=1.0E-6_r8, WC2MIN=1.00_r8, GMF1=GMF/AA1, GMF5=GMF/F5)
    !     parameter (sialf=0.5_r8)
    !
    REAL(kind=r8),PARAMETER :: PI=3.1415926535897931_r8, PIINV=1.0_r8/PI
    INTEGER :: ITR
    !     PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=7)
!PK     INTEGER,PARAMETER :: ITRMU=25, ITRMD=25, ITRMIN=12, ITRMND=12
    INTEGER,PARAMETER :: ITRMU=25, ITRMD=25, ITRMIN=12, ITRMND=12

    !     PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=12)
    !     PARAMETER (ITRMU=14, ITRMD=18, ITRMIN=7)
    !     PARAMETER (ITRMU=10, ITRMD=10, ITRMIN=5)
    REAL(kind=r8) :: QRP(1:K+1)
    REAL(kind=r8) :: WVL(1:K+1)
    REAL(kind=r8) :: AL2
    REAL(kind=r8) :: WVLO(1:K+1)
    !
    REAL(kind=r8) :: RNF(1:K)
    REAL(kind=r8) :: ROR(1:K+1)
    REAL(kind=r8) :: STLT(1:K)
    REAL(kind=r8) :: RNT
    REAL(kind=r8) :: RNB
    REAL(kind=r8) :: ERRQ
    REAL(kind=r8) :: RNTP
    INTEGER       :: IDW 
    INTEGER       :: IDN(K)
    INTEGER       :: idnm
    !REAL(kind=r8) :: ELM(K)
    !     real(kind=r8) EM(K*K), ELM(K)
    REAL(kind=r8) :: EDZ
    REAL(kind=r8) :: DDZ
    REAL(kind=r8) :: CE
    REAL(kind=r8) :: QHS
    REAL(kind=r8) :: FAC
    REAL(kind=r8) :: FACG
   ! REAL(kind=r8) :: ASIN
    REAL(kind=r8) :: RSUM1
    REAL(kind=r8) :: RSUM2
!    REAL(kind=r8) :: RSUM3
!    REAL(kind=r8) :: CEE

    LOGICAL       :: DDLGK
    !
    REAL(kind=r8) :: AA(1:K,1:K+1)
    REAL(kind=r8) :: QW(1:K,1:K)          
    REAL(kind=r8) :: BUD(1:K)
    REAL(kind=r8) :: VT(2)
    REAL(kind=r8) :: VRW(2)
    REAL(kind=r8) :: TRW(2)    
    REAL(kind=r8) :: GQW(1:K)         
    REAL(kind=r8) :: QA(3)
    REAL(kind=r8) :: WA(3)
    REAL(kind=r8) :: DOFW         
    REAL(kind=r8) :: QRPI(1:K)
    !REAL(kind=r8) :: QRPS(1:K)


    !    &,                    GQW(KD:K), WCB(KD:K)

    !***********************************************************************
    DOF    =0.0_r8
    CLDFRD(1:K) =0.0_r8
    WCB(1:K) =0.0_r8
    GMS(1:K+1) =0.0_r8
    GSD(1:K) =0.0_r8
    GHD(1:K) =0.0_r8

   ! REAL(kind=r8) :: QRPF
    DOF=0.0_r8
    CLDFRD(1:K)  =0.0_r8
    WCB(1:K)=0.0_r8
    GMS(1:K+1)=0.0_r8
    GSD(1:K)=0.0_r8
    GHD(1:K)=0.0_r8


    KB1=0

    RNS(1:K) =0.0_r8
    TX1=0.0_r8
    TX2=0.0_r8
    TX3=0.0_r8
    TX4=0.0_r8
    TX5=0.0_r8
    TX6=0.0_r8
    TX8=0.0_r8
    TX9=0.0_r8

    ST1=0.0_r8
    QQQ=0.0_r8
    PICON=0.0_r8
    DEL_ETA=0.0_r8
    TEM=0.0_r8
    TEM1=0.0_r8
    TEM2=0.0_r8
    TEM3=0.0_r8
    TEM4=0.0_r8
    ST2=0.0_r8
    ERRH=0.0_r8
    ERRW=0.0_r8
    ERRE=0.0_r8
    TEM6=0.0_r8

    L=0  
    N=0  
    KD1=0 
    II=0      
    KP1=0  
    KM1=0 
    KTEM=0 
    KK=0 
    KK1=0 
    LM1=0 
    LL=0 
    LP1=0
    ntla=0
    STLA=0.0_r8
    CTL2=0.0_r8
    CTL3=0.0_r8
    VTPEXP =0.0_r8
    WCMIN=0.0_r8
    WCBASE=0.0_r8
    F2=0.0_r8
    QRAF=0.0_r8
    QRBF=0.0_r8
    CMPOR =0.0_r8
    del_tla=0.0_r8
    !
    ITR=0
    QRP=0.0_r8
    WVL=0.0_r8
    AL2=0.0_r8
    WVLO=0.0_r8
    !
    RNF=0.0_r8
    ROR=0.0_r8
    STLT=0.0_r8
    RNT=0.0_r8
    RNB=0.0_r8
    ERRQ=0.0_r8
    RNTP=0.0_r8
    IDW =0
    IDN(K)=0
    idnm=0
    EDZ=0.0_r8
    DDZ=0.0_r8
    CE=0.0_r8
    QHS=0.0_r8
    FAC=0.0_r8
    FACG=0.0_r8
    RSUM1=0.0_r8
    RSUM2=0.0_r8
    AA=0.0_r8
    QW=0.0_r8
    BUD=0.0_r8
    VT=0.0_r8
    VRW=0.0_r8
    TRW=0.0_r8
    GQW=0.0_r8
    QA=0.0_r8
    WA=0.0_r8
    DOFW=0.0_r8
    QRPI=0.0_r8
    !

    !     if(lprnt) print *,' K=',K,' KD=',KD,' In Downdrft'

    KD1    = KD + 1
    KP1    = K  + 1
    KM1    = K  - 1
    KB1    = KBL - 1
    !
    CMPOR  = CMB2PA / con_rd
    !
    !     VTP    = 36.34_r8*SQRT(1.2_r8)* (0.001_r8)**0.1364_r8
    VTPEXP = -0.3636_r8
    !     PIINV  = 1.0_r8 / PI
    PICON  = PI * ONEBG * 0.5_r8
    !
    !
    !     Compute Rain Water Budget of the Updraft (Cheng and Arakawa, 1997)
    !
    CLDFRD = 0.0_r8
    RNTP   = 0.0_r8
    DOF    = 0.0_r8
    ERRQ   = 10.0_r8
    RNB    = 0.0_r8
    RNT    = 0.0_r8
    TX2    = PRL(KBL)
    !
    TX1      = (PRL(KD) + PRL(KD1)) * 0.5_r8
    ROR(KD)  = CMPOR*TX1 / (TOL(KD)*(1.0_r8+con_FVirt*QOL(KD)))
    !     GMS(KD)  = VTP * ROR(KD) ** VTPEXP
    GMS(KD)  = VTP * VTPF(ROR(KD))
    !
    QRP(KD)  = QRMIN
    !
    TEM      = TOL(K) * (1.0_r8 + con_FVirt * QOL(K))
    ROR(K+1) = 0.5_r8 * CMPOR * (PRL(K+1)+PRL(K)) / TEM
    GMS(K+1) = VTP * VTPF(ROR(K+1))
    QRP(K+1) = QRMIN
    !
    kk = kbl
    DO L=KD1,K
       TEM = 0.5_r8 * (TOL(L)+TOL(L-1))                                   &
            &      * (1.0_r8 + (0.5_r8*con_FVirt) * (QOL(L)+QOL(L-1)))
       ROR(L) = CMPOR * PRL(L) / TEM
       !       GMS(L) = VTP * ROR(L) ** VTPEXP
       GMS(L) = VTP * VTPF(ROR(L))
       QRP(L) = QRMIN
       IF (buy(l) .LE. 0.0_r8 .AND. kk .EQ. KBL) THEN
          kk = l
       ENDIF
    ENDDO
    IF (kk .NE. kbl) THEN
       DO l=kk,kbl
          buy(l) = 0.9_r8 * buy(l-1)
       ENDDO
    ENDIF
    !
    DO l=kd,k
       qrpi(l) = buy(l)
    ENDDO
    DO l=kd1,kb1
       buy(l) = 0.25_r8 * (qrpi(l-1)+qrpi(l)+qrpi(l)+qrpi(l+1))
    ENDDO

    !
    !     CALL ANGRAD(TX1, ALM, STLA, CTL2, AL2, PI, TLA, TX2, WFN, TX3)
    tx1 = 1000.0_r8 + tx1 - prl(k+1)
    CALL ANGRAD(TX1, ALM,  AL2, TLA)
    !
    !    Following Ucla approach for rain profile
    !
    F2      = 2.0_r8*BB1*ONEBG/(PI*0.2_r8)
    WCMIN   = SQRT(WC2MIN)
    WCBASE  = WCMIN
    !
    !     del_tla = TLA * 0.2_r8
    !     del_tla = TLA * 0.25_r8
    del_tla = TLA * 0.3_r8
    TLA     = TLA - DEL_TLA
    !
    DO L=KD,K
       RNF(L)   = 0.0_r8
       RNS(L)   = 0.0_r8
       WVL(L)   = 0.0_r8
       STLT(L)  = 0.0_r8
       GQW(L)   = 0.0_r8
       QRP(L)   = QRMIN
       DO N=KD,K
          QW(N,L) = 0.0_r8
       ENDDO
    ENDDO
    !
    !-----QW(N,L) = D(W(N)*W(N))/DQR(L)
    !
    KK = KBL
    QW(KD,KD)  = -QRB(KD)  * GMF1
    GHD(KD)    = ETA(KD)   * ETA(KD)
    GQW(KD)    = QW(KD,KD) * GHD(KD)
    GSD(KD)    = ETAI(KD)  * ETAI(KD)
    !
    GQW(KK)    = -  QRB(KK-1) * (GMF1+GMF1)
    !
    WCB(KK)    = WCBASE * WCBASE

    TX1        = WCB(KK)
    GSD(KK)    = 1.0_r8
    GHD(KK)    = 1.0_r8
    !
    TEM        = GMF1 + GMF1
    DO L=KB1,KD1,-1
       GHD(L)  = ETA(L)  * ETA(L)
       GSD(L)  = ETAI(L) * ETAI(L)
       GQW(L)  = - GHD(L) * (QRB(L-1)+QRT(L)) * TEM
       QW(L,L) = - QRT(L) * TEM
       !
       st1     = 0.5_r8 * (eta(l) + eta(l+1))
       TX1     = TX1 + BUY(L) * TEM * (qrb(l)+qrt(l)) * st1 * st1
       WCB(L)  = TX1 * GSD(L)
    ENDDO
    !
    TEM1        = (QRB(KD) + QRT(KD1) + QRT(KD1)) * GMF1
    GQW(KD1)    = - GHD(KD1) * TEM1
    QW(KD1,KD1) = - QRT(KD1) * TEM
    st1     = 0.5_r8 * (eta(kd) + eta(kd1))
    WCB(KD)     = (TX1 + BUY(KD)*TEM*qrb(kd)*st1*st1) * GSD(KD)
    !
    DO L=KD1,KBL
       DO N=KD,L-1
          QW(N,L) = GQW(L) * GSD(N)
       ENDDO
    ENDDO
    QW(KBL,KBL) = 0.0_r8
    !
    DO ntla=1,numtla
       !
       !     if (errq .lt. 1.0_r8 .or. tla .gt. 45.0_r8) cycle
       IF (errq .LT. 0.1_r8 .OR. tla .GT. 45.0_r8) CYCLE
       !
       tla = tla + del_tla
       STLA = SIN(TLA*PI/180.0_r8)
       CTL2 = 1.0_r8 - STLA * STLA
       !
       !     if (lprnt) print *,' tla=',tla,' al2=',al2,' ptop='
       !    &,0.5_r8*(prl(kd)+prl(kd1)),' ntla=',ntla,' f2=',f2,' stla=',stla
       !     if (lprnt) print *,' buy=',(buy(l),l=kd,kbl)
       !
       STLA = F2     * STLA * AL2
       CTL2 = DD1    * CTL2
       CTL3 = 0.1364_r8 * CTL2
       !
       DO L=KD,K
          RNF(L)   = 0.0_r8
          WVL(L)   = 0.0_r8
          STLT(L)  = 0.0_r8
          QRP(L)   = QRMIN
       ENDDO
       WVL(KBL)    = WCBASE
       STLT(KBL)   = 1.0_r8 / WCBASE
       !
       DO L=KD,K+1
          DO N=KD,K
             AA(N,L) = 0.0_r8
          ENDDO
       ENDDO
       !
       SKPUP = .FALSE.
       !
       DO ITR=1,ITRMU               ! Rain Profile Iteration starts!
          IF (.NOT. SKPUP) THEN
             wvlo = wvl
             !
             !-----CALCULATING THE VERTICAL VELOCITY
             !
             TX1      = 0.0_r8
             QRPI(KBL) = 1.0_r8 / QRP(KBL)
             DO L=KB1,KD,-1
                TX1     = TX1    + QRP(L+1) * GQW(L+1)
                ST1     = WCB(L) + QW(L,L)  * QRP(L)                        &
                     &                       + TX1      * GSD(L)
                IF (st1 .GT. wc2min) THEN
                   !             WVL(L)  = SQRT(ST1)
                   WVL(L)  = 0.5_r8 * (SQRT(ST1) + WVL(L))
                   !             if (itr .eq. 1) wvl(l) = wvl(l) * 0.25_r8
                ELSE
                   !     if (lprnt)  print *,' l=',l,' st1=',st1,' wcb=',wcb(l),' qw='
                   !    &,qw(l,l),' qrp=',qrp(l),' tx1=',tx1,' gsd=',gsd(l),' ite=',itr
                   !             wvl(l) = 0.5_r8*(wcmin+wvl(l))
                   wvl(l) = 0.5_r8 * (wvl(l) + wvl(l+1))
                   qrp(l) = 0.5_r8 * ((wvl(l)*wvl(l)-wcb(l)-tx1*gsd(l))/qw(l,l) &
                        &                     + qrp(l))
                   !!            wvl(l) = 0.5_r8 * (wvl(l) + wvl(l+1))
                ENDIF
                !           wvl(l)  = 0.5_r8 * (wvl(l) + wvlo(l))
                !           WVL(L)  = SQRT(MAX(ST1,WC2MIN))
                wvl(l)  = MAX(wvl(l), wcbase)
                STLT(L) = 1.0_r8 / WVL(L)
                QRPI(L) = 1.0_r8 / QRP(L)
             ENDDO
             !
             !     if (lprnt) then
             !     print *,' ITR=',ITR,' ITRMU=',ITRMU
             !     print *,' WVL=',(WVL(L),L=KD,KBL)
             !     print *,' qrp=',(qrp(L),L=KD,KBL)
             !     print *,' qrpi=',(qrpi(L),L=KD,KBL)
             !     print *,' rnf=',(rnf(L),L=KD,KBL)
             !     endif
             !
             !-----CALCULATING TRW, VRW AND OF
             !
             !         VT(1)   = GMS(KD) * QRP(KD)**0.1364_r8
             VT(1)   = GMS(KD) * QRPF(QRP(KD))
             TRW(1)  = ETA(KD) * QRP(KD) * STLT(KD)
             TX6     = TRW(1) * VT(1)
             VRW(1)  = F3*WVL(KD) - CTL2*VT(1)
             BUD(KD) = STLA * TX6 * QRB(KD) * 0.5_r8
             RNF(KD) = BUD(KD)
             DOF     = 1.1364_r8 * BUD(KD) * QRPI(KD)
             DOFW    = -BUD(KD) * STLT(KD)
             !
             RNT     = TRW(1) * VRW(1)
             TX2     = 0.0_r8
             TX4     = 0.0_r8
             RNB     = RNT
             TX1     = 0.5_r8
             TX8     = 0.0_r8
             !
             IF (RNT .GE. 0.0_r8) THEN
                TX3 = (RNT-CTL3*TX6) * QRPI(KD)
                TX5 = CTL2 * TX6 * STLT(KD)
             ELSE
                TX3 = 0.0_r8
                TX5 = 0.0_r8
                RNT = 0.0_r8
                RNB = 0.0_r8
             ENDIF
             !
             DO L=KD1,KB1
                KTEM    = MAX(L-2, KD)
                LL      = L - 1
                !
                !           VT(2)   = GMS(L) * QRP(L)**0.1364_r8
                VT(2)   = GMS(L) * QRPF(QRP(L))
                TRW(2)  = ETA(L) * QRP(L) * STLT(L)
                VRW(2)  = F3*WVL(L) - CTL2*VT(2)
                QQQ     = STLA * TRW(2) * VT(2)
                ST1     = TX1  * QRB(LL)
                BUD(L)  = QQQ * (ST1 + QRT(L))
                !
                QA(2)   = DOF
                WA(2)   = DOFW
                DOF     = 1.1364_r8 * BUD(L) * QRPI(L)
                DOFW    = -BUD(L) * STLT(L)
                !
                RNF(LL) = RNF(LL) + QQQ * ST1
                RNF(L)  =           QQQ * QRT(L)
                !
                TEM3    = VRW(1) + VRW(2)
                TEM4    = TRW(1) + TRW(2)
                !
                TX6     = .25_r8 * TEM3 * TEM4
                TEM4    = TEM4 * CTL3
                !
                !-----BY QR ABOVE
                !
                !           TEM1    = .25_r8*(TRW(1)*TEM3 - TEM4*VT(1))*TX7
                TEM1    = .25_r8*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
                ST1     = .25_r8*(TRW(1)*(CTL2*VT(1)-VRW(2))                   &
                     &                  * STLT(LL) + F3*TRW(2))
                !-----BY QR BELOW
                TEM2    = .25_r8*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
                ST2     = .25_r8*(TRW(2)*(CTL2*VT(2)-VRW(1))                   &
                     &                 * STLT(L)  + F3*TRW(1))
                !
                !      From top to  the KBL-2 layer
                !
                QA(1)   = TX2
                QA(2)   = QA(2) + TX3 - TEM1
                QA(3)   = -TEM2
                !
                WA(1)   = TX4
                WA(2)   = WA(2) + TX5 - ST1
                WA(3)   = -ST2
                !
                TX2     = TEM1
                TX3     = TEM2
                TX4     = ST1
                TX5     = ST2
                !
                VT(1)   = VT(2)
                TRW(1)  = TRW(2)
                VRW(1)  = VRW(2)
                !
                IF (WVL(KTEM) .EQ. WCMIN) WA(1) = 0.0_r8
                IF (WVL(LL)   .EQ. WCMIN) WA(2) = 0.0_r8
                IF (WVL(L)    .EQ. WCMIN) WA(3) = 0.0_r8
                DO N=KTEM,KBL
                   AA(LL,N) = (WA(1)*QW(KTEM,N) * STLT(KTEM)                 &
                        &                 +  WA(2)*QW(LL,N)   * STLT(LL)                   &
                        &                 +  WA(3)*QW(L,N)    * STLT(L) ) * 0.5_r8
                ENDDO
                AA(LL,KTEM) = AA(LL,KTEM) + QA(1)
                AA(LL,LL)   = AA(LL,LL)   + QA(2)
                AA(LL,L)    = AA(LL,L)    + QA(3)
                BUD(LL)     = (TX8 + RNN(LL)) * 0.5_r8                         &
                     &                    - RNB + TX6 - BUD(LL)
                AA(LL,KBL+1) = BUD(LL)
                RNB = TX6
                TX1 = 1.0_r8
                TX8 = RNN(LL)
             ENDDO
             L  = KBL
             LL = L - 1
             !         VT(2)   = GMS(L) * QRP(L)**0.1364_r8
             VT(2)   = GMS(L) * QRPF(QRP(L))
             TRW(2)  = ETA(L) * QRP(L) * STLT(L)
             VRW(2)  = F3*WVL(L) - CTL2*VT(2)
             ST1     = STLA * TRW(2) * VT(2) * QRB(LL)
             BUD(L)  = ST1

             QA(2)   = DOF
             WA(2)   = DOFW
             DOF     = 1.1364_r8 * BUD(L) * QRPI(L)
             DOFW    = -BUD(L) * STLT(L)
             !
             RNF(LL) = RNF(LL) + ST1
             !
             TEM3    = VRW(1) + VRW(2)
             TEM4    = TRW(1) + TRW(2)
             !
             TX6     = .25_r8 * TEM3 * TEM4
             TEM4    = TEM4 * CTL3
             !
             !-----BY QR ABOVE
             !
             TEM1    = .25_r8*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
             ST1     = .25_r8*(TRW(1)*(CTL2*VT(1)-VRW(2))                     &
                  &                * STLT(LL) + F3*TRW(2))
             !-----BY QR BELOW
             TEM2    = .25_r8*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
             ST2     = .25_r8*(TRW(2)*(CTL2*VT(2)-VRW(1))                     &
                  &                 * STLT(L)  + F3*TRW(1))
             !
             !      For the layer next to the top of the boundary layer
             !
             QA(1)   = TX2
             QA(2)   = QA(2) + TX3 - TEM1
             QA(3)   = -TEM2
             !
             WA(1)   = TX4
             WA(2)   = WA(2) + TX5 - ST1
             WA(3)   = -ST2
             !
             TX2     = TEM1
             TX3     = TEM2
             TX4     = ST1
             TX5     = ST2
             !
             IDW     = MAX(L-2, KD)
             !
             IF (WVL(IDW) .EQ. WCMIN) WA(1) = 0.0_r8
             IF (WVL(LL)  .EQ. WCMIN) WA(2) = 0.0_r8
             IF (WVL(L)   .EQ. WCMIN) WA(3) = 0.0_r8
             !
             KK = IDW
             DO N=KK,L
                AA(LL,N) = (WA(1)*QW(KK,N) * STLT(KK)                       &
                     &               +  WA(2)*QW(LL,N) * STLT(LL)                       &
                     &               +  WA(3)*QW(L,N)  * STLT(L) ) * 0.5_r8

             ENDDO
             !
             AA(LL,IDW) = AA(LL,IDW) + QA(1)
             AA(LL,LL)  = AA(LL,LL)  + QA(2)
             AA(LL,L)   = AA(LL,L)   + QA(3)
             BUD(LL)    = (TX8+RNN(LL)) * 0.5_r8 - RNB + TX6 - BUD(LL)
             !
             AA(LL,L+1) = BUD(LL)
             !
             RNB        = TRW(2) * VRW(2)
             !
             !      For the top of the boundary layer
             !
             IF (RNB .LT. 0.0_r8) THEN
                KK    = KBL
                TEM   = VT(2) * TRW(2)
                QA(2) = (RNB - CTL3*TEM) * QRPI(KK)
                WA(2) = CTL2 * TEM * STLT(KK)
             ELSE
                RNB   = 0.0_r8
                QA(2) = 0.0_r8
                WA(2) = 0.0_r8
             ENDIF
             !
             QA(1) = TX2
             QA(2) = DOF + TX3 - QA(2)
             QA(3) = 0.0_r8
             !
             WA(1) = TX4
             WA(2) = DOFW + TX5 - WA(2)
             WA(3) = 0.0_r8
             !
             KK = KBL
             IF (WVL(KK-1) .EQ. WCMIN) WA(1) = 0.0_r8
             IF (WVL(KK)   .EQ. WCMIN) WA(2) = 0.0_r8
             !
             DO II=1,2
                N = KK + II - 2
                AA(KK,N) = (WA(1)*QW(KK-1,N) * STLT(KK-1)                  &
                     &                +  WA(2)*QW(KK,N)   * STLT(KK)) * 0.5_r8
             ENDDO
             FAC = 0.5_r8
             LL  = KBL
             L   = LL + 1
             LM1 = LL - 1
             AA(LL,LM1)  = AA(LL,LM1) + QA(1)
             AA(LL,LL)   = AA(LL,LL)  + QA(2)
             BUD(LL)     = 0.5_r8*RNN(LM1) - TX6 + RNB - BUD(LL)
             AA(LL,LL+1) = BUD(LL)
             !
             !-----SOLVING THE BUDGET EQUATIONS FOR DQR
             !
             DO L=KD1,KBL
                LM1  = L - 1
                UNSAT = ABS(AA(LM1,LM1)) .LT. ABS(AA(L,LM1))
                DO  N=LM1,KBL+1
                   IF (UNSAT) THEN
                      TX1       = AA(LM1,N)
                      AA(LM1,N) = AA(L,N)
                      AA(L,N)   = TX1
                   ENDIF
                ENDDO
                TX1 = AA(L,LM1) / AA(LM1,LM1)
                DO  N=L,KBL+1
                   AA(L,N) = AA(L,N) - TX1 * AA(LM1,N)
                ENDDO
             ENDDO
             !
             !-----BACK SUBSTITUTION AND CHECK IF THE SOLUTION CONVERGES
             !
             KK = KBL
             KK1 = KK + 1
             AA(KK,KK1) = AA(KK,KK1) / AA(KK,KK)      !   Qr correction !
             TX2        = ABS(AA(KK,KK1)) * QRPI(KK)  !   Error Measure !
             !     if (lprnt) print *,' tx2a=',tx2,' aa1=',aa(kk,kk1)
             !    &,' qrpi=',qrpi(kk)
             !
             KK = KBL + 1
             DO L=KB1,KD,-1
                LP1   = L + 1
                TX1  = 0.0_r8
                DO N=LP1,KBL
                   TX1  = TX1 + AA(L,N) * AA(N,KK)
                ENDDO
                AA(L,KK) = (AA(L,KK) - TX1) / AA(L,L)       ! Qr correction !
                TX2      = MAX(TX2, ABS(AA(L,KK))*QRPI(L))  ! Error Measure !
                !     if (lprnt) print *,' tx2b=',tx2,' aa1=',aa(l,kk)
                !    &,' qrpi=',qrpi(l),' L=',L
             ENDDO
             !
             !         tem = 0.5_r8
             IF (tx2 .GT. 1.0_r8 .AND. ABS(errq-tx2) .GT. 0.1_r8) THEN
                tem = 0.5_r8
                !!        elseif (tx2 .lt. 0.1_r8) then
                !!          tem = 1.2_r8
             ELSE
                tem = 1.0_r8
             ENDIF
             !
             DO L=KD,KBL
                !            QRP(L) = MAX(QRP(L)+AA(L,KBL+1), QRMIN)
                QRP(L) = MAX(QRP(L)+AA(L,KBL+1)*tem, QRMIN)
             ENDDO
             !
             !     if (lprnt) print *,' itr=',itr,' tx2=',tx2
             IF (ITR .LT. ITRMIN) THEN
                TEM = ABS(ERRQ-TX2) 
                IF (TEM .GE. ERRMI2 .AND. TX2 .GE. ERRMIN) THEN 
                   ERRQ  = TX2                              ! Further iteration !
                ELSE 
                   SKPUP = .TRUE.                           ! Converges      !
                   ERRQ  = 0.0_r8                              ! Rain profile exists!
                   !     print *,' here1',' tem=',tem,' tx2=',tx2,' errmi2=',
                   !    *errmi2,' errmin=',errmin
                ENDIF
             ELSE
                TEM = ERRQ - TX2
                !            IF (TEM .LT. ZERO .AND. ERRQ .GT. 0.1_r8) THEN
                IF (TEM .LT. ZERO .AND. ERRQ .GT. 0.5_r8) THEN
                   !            IF (TEM .LT. ZERO .and.                                    &
                   !    &          (ntla .lt. numtla .or. ERRQ .gt. 0.5_r8)) THEN
                   !     if (lprnt) print *,' tx2=',tx2,' errq=',errq,' tem=',tem
                   SKPUP = .TRUE.                           ! No convergence !
                   ERRQ = 10.0_r8                              ! No rain profile!
!!!!         ELSEIF (ABS(TEM).LT.ERRMI2 .OR. TX2.LT.ERRMIN) THEN
                ELSEIF (TX2.LT.ERRMIN) THEN
                   SKPUP = .TRUE.                           ! Converges      !
                   ERRQ = 0.0_r8                               ! Rain profile exists!
                   !     print *,' here2'
                ELSEIF (tem .LT. zero .AND. errq .LT. 0.1_r8) THEN
                   skpup = .TRUE.
                   !              if (ntla .eq. numtla .or. tem .gt. -0.003) then
                   errq  = 0.0_r8
                   !              else
                   !                errq = 10.0
                   !              endif
                ELSE
                   ERRQ = TX2                               ! Further iteration !
                   !     if (lprnt) print *,' itr=',itr,' errq=',errq
                   !              if (itr .eq. itrmu .and. ERRQ .GT. ERRMIN*10             &
                   !    &            .and. ntla .eq. 1) ERRQ = 10.0 
                ENDIF
             ENDIF
             !
             !         if (lprnt) print *,' ERRQ=',ERRQ

          ENDIF                                           ! SKPUP  ENDIF!
          !
       ENDDO                                          ! End of the ITR Loop!!
       !     enddo                                          ! End of ntla loop
       !
       !     if(lprnt) then
       !       print *,' QRP=',(QRP(L),L=KD,KBL)
       !       print *,'RNF=',(RNF(L),L=KD,KBL),' RNT=',RNT,' RNB=',RNB
       !    &,' errq=',errq
       !     endif
       !
       IF (ERRQ .LT. 0.1_r8) THEN
          DDFT = .TRUE.
          RNB  = - RNB
          !    do l=kd1,kb1-1
          !      if (wvl(l)-wcbase .lt. 1.0E-9) ddft = .false.
          !    enddo
       ELSE
          DDFT = .FALSE.
       ENDIF
       !
       !     Caution !! Below is an adjustment to rain flux to maintain
       !                conservation of precip!
       !
       IF (DDFT) THEN
          TX1 = 0.0_r8
          DO L=KD,KB1
             TX1 = TX1 + RNF(L)
          ENDDO
          !     if (lprnt) print *,' tx1+rnt+rnb=',tx1+rnt+rnb, ' train=',train
          TX1 = TRAIN / (TX1+RNT+RNB)
          IF (ABS(TX1-1.0_r8) .LT. 0.2_r8) THEN
             RNT = MAX(RNT*TX1,ZERO)
             RNB = RNB * TX1
          ELSE
             DDFT = .FALSE.
             ERRQ = 10.0_r8
          ENDIF
       ENDIF
    ENDDO                                          ! End of ntla loop
    !
    DOF = 0.0_r8
    IF (.NOT. DDFT) RETURN     ! Rain profile did not converge!
    !

    DO L=KD,KB1
       RNF(L) = RNF(L) * TX1

    ENDDO
    !     if (lprnt) print *,' TRAIN=',TRAIN
    !     if (lprnt) print *,' RNF=',RNF
    !
    !     Adjustment is over
    !
    !     Downdraft
    !
    DO L=KD,K
       WCB(L) = 0.0_r8
    ENDDO
    !
    SKPDD = .NOT. DDFT
    !
    ERRQ  = 10.0_r8
    IF (.NOT. SKPDD) THEN
       !
       !     Calculate Downdraft Properties
       !

       KK = MAX(KB1,KD1)
       DO L=KK,K
          STLT(L) = STLT(L-1)
       ENDDO
       TEM1 = 1.0_r8 / BB1
       !
       DO L=KD,K
          IF (L .LE. KBL) THEN
             TEM     = STLA * TEM1
             STLT(L) = ETA(L) * STLT(L) * TEM / ROR(L)
          ELSE
             STLT(L) = 0.0_r8
          ENDIF
       ENDDO
       !       if (lprnt) print *,' STLT=',stlt

       rsum1 = 0.0_r8
       rsum2 = 0.0_r8

       !
       IDN      = 99
       DO L=KD,K+1
          ETD(L)  = 0.0_r8
          WVL(L)  = 0.0_r8
          !         QRP(L)  = 0.0_r8
       ENDDO
       DO L=KD,K
          EVP(L)   = 0.0_r8
          BUY(L)   = 0.0_r8
          QRP(L+1) = 0.0_r8
       ENDDO
       HOD(KD)  = HOL(KD)
       QOD(KD)  = QOL(KD)
       TX1      = 0.0_r8                               ! sigma at the top
!!!     TX1      = STLT(KD)*QRB(KD)*ONE              ! sigma at the top
       !       TX1      = MIN(STLT(KD)*QRB(KD)*ONE, ONE)    ! sigma at the top
       !       TX1      = MIN(STLT(KD)*QRB(KD)*0.5_r8, ONE)    ! sigma at the top
       RNTP     = 0.0_r8
       TX5      = TX1
       QA(1)    = 0.0_r8
       !     if(lprnt) print *,' stlt=',stlt(kd),' qrb=',qrb(kd)
       !    *,' tx1=',tx1,' ror=',ror(kd),' gms=',gms(kd),' rpart=',rpart
       !    *,' rnt=',rnt
       !
       !       Here we assume RPART of detrained rain RNT goes to Pd
       !
       IF (RNT .GT. 0.0_r8) THEN
          IF (TX1 .GT. 0.0_r8) THEN
             QRP(KD) = (RPART*RNT / (ROR(KD)*TX1*GMS(KD)))               &
                  &                                          ** (1.0_r8/1.1364_r8)
          ELSE
             tx1 = RPART*RNT / (ROR(KD)*GMS(KD)*QRP(KD)**1.1364_r8)
          ENDIF
          RNTP    = (1.0_r8 - RPART) * RNT
          BUY(KD) = - ROR(KD) * TX1 * QRP(KD)
       ELSE
          QRP(KD) = 0.0_r8
       ENDIF
       !
       !     L-loop for the downdraft iteration from KD1 to K+1 (bottom surface)
       !
       !     BUD(KD) = ROR(KD)
       idnm = 1
       DO L=KD1,K+1

          QA(1) = 0.0_r8
          ddlgk = idn(idnm) .EQ. 99
          IF (.NOT. ddlgk) CYCLE
          IF (L .LE. K) THEN
             ST1   = 1.0_r8 - ALFIND(L)
             WA(1) = ALFIND(L)*HOL(L-1) + ST1*HOL(L)
             WA(2) = ALFIND(L)*QOL(L-1) + ST1*QOL(L)
             WA(3) = ALFIND(L)*TOL(L-1) + ST1*TOL(L)
             QA(2) = ALFIND(L)*HST(L-1) + ST1*HST(L)
             QA(3) = ALFIND(L)*QST(L-1) + ST1*QST(L)
          ELSE
             WA(1) = HOL(K)
             WA(2) = QOL(K)
             WA(3) = TOL(K)
             QA(2) = HST(K)
             QA(3) = QST(K)
          ENDIF
          !
          FAC = 2.0_r8
          IF (L .EQ. KD1) FAC = 1.0_r8

          FACG    = FAC * 0.5_r8 * GMF5     !  12/17/97
          !
          !         DDLGK   =  IDN(idnm) .EQ. 99
          BUD(KD) = ROR(L)

          !         IF (DDLGK) THEN
          TX1    = TX5
          WVL(L) = MAX(WVL(L-1),ONE_M1)

          QRP(L) = MAX(QRP(L-1),QRP(L))
          !
          !           VT(1)  = GMS(L-1) * QRP(L-1) ** 0.1364_r8
          VT(1)  = GMS(L-1) * QRPF(QRP(L-1))
          RNT    = ROR(L-1) * (WVL(L-1)+VT(1))*QRP(L-1)
          !     if(lprnt) print *,' l=',l,' qa=',qa(1), ' tx1RNT=',RNT*tx1,
          !    *' wvl=',wvl(l-1)
          !    *,' qrp=',qrp(l-1),' tx5=',tx5,' tx1=',tx1,' rnt=',rnt

          !

          !           TEM    = MAX(ALM, 2.5E-4) * MAX(ETA(L), 1.0)
          TEM    = MAX(ALM,ONE_M6) * MAX(ETA(L), ONE)
          !           TEM    = MAX(ALM, 1.0E-5) * MAX(ETA(L), 1.0)
          TRW(1) = PICON*TEM*(QRB(L-1)+QRT(L-1))
          TRW(2) = 1.0_r8 / TRW(1)
          !
          VRW(1) = 0.5_r8 * (GAM(L-1) + GAM(L))
          VRW(2) = 1.0_r8 / (VRW(1) + VRW(1))
          !
          TX4    =  (QRT(L-1)+QRB(L-1))*(ONEBG*FAC*500.00_r8*EKNOB)
          !
          DOFW   = 1.0_r8 / (WA(3) * (1.0_r8 + con_FVirt*WA(2)))      !  1.0_r8 / TVbar!
          !
          ETD(L) = ETD(L-1)
          HOD(L) = HOD(L-1)
          QOD(L) = QOD(L-1)
          !
          ERRQ   = 10.0_r8

          !
          IF (L .LE. KBL) THEN
             TX3 = STLT(L-1) * QRT(L-1) * (0.5_r8*FAC)
             TX8 = STLT(L)   * QRB(L-1) * (0.5_r8*FAC)
             TX9 = TX8 + TX3
          ELSE
             TX3 = 0.0_r8
             TX8 = 0.0_r8
             TX9 = 0.0_r8
          ENDIF
          !
          TEM  = WVL(L-1) + VT(1)
          IF (TEM .GT. 0.0_r8) THEN
             TEM1 = 1.0_r8 / (TEM*ROR(L-1))
             TX3 = VT(1) * TEM1 * ROR(L-1) * TX3
             TX6 = TX1 * TEM1
          ELSE
             TX6 = 1.0_r8
          ENDIF
          !         ENDIF
          !
          IF (L .EQ. KD1) THEN
             IF (RNT .GT. 0.0_r8) THEN
                TEM    = MAX(QRP(L-1),QRP(L))
                WVL(L) = TX1 * TEM * QRB(L-1)*(FACG*5.0_r8)
             ENDIF
             WVL(L) = MAX(ONE_M2, WVL(L))
             TRW(1) = TRW(1) * 0.5_r8
             TRW(2) = TRW(2) + TRW(2)
          ELSE
             IF (DDLGK) EVP(L-1) = EVP(L-2)
          ENDIF
          !
          !       No downdraft above level IDH
          !

          IF (L .LT. IDH) THEN

             ETD(L)   = 0.0_r8
             HOD(L)   = WA(1)
             QOD(L)   = WA(2)
             EVP(L-1) = 0.0_r8
             WVL(L)   = 0.0_r8
             QRP(L)   = 0.0_r8
             BUY(L)   = 0.0_r8
             TX5      = TX9
             ERRQ     = 0.0_r8
             RNTP     = RNTP + RNT * TX1
             RNT      = 0.0_r8
             WCB(L-1) = 0.0_r8
          ENDIF
          !         BUD(KD) = ROR(L)
          !
          !       Iteration loop for a given level L begins
          !
          !         if (lprnt) print *,' tx8=',tx8,' tx9=',tx9,' tx5=',tx5
          !    &,                      ' tx1=',tx1
          DO ITR=1,ITRMD
             !
             !           UNSAT =  DDLGK .AND. (ERRQ .GT. ERRMIN)
             UNSAT =  ERRQ .GT. ERRMIN
             IF (UNSAT) THEN
                !
                !             VT(1)  = GMS(L) * QRP(L) ** 0.1364
                VT(1)  = GMS(L) * QRPF(QRP(L))
                TEM    =  WVL(L) + VT(1)
                !
                IF (TEM .GT. 0.0_r8) THEN
                   ST1    = ROR(L) * TEM * QRP(L) + RNT
                   IF (ST1 .NE. 0.0_r8) ST1 = 2.0_r8 * EVP(L-1) / ST1
                   TEM1   = 1.0_r8 / (TEM*ROR(L))
                   TEM2   = VT(1) * TEM1 * ROR(L) * TX8
                ELSE
                   TEM1   = 0.0_r8
                   TEM2   = TX8
                   ST1    = 0.0_r8
                ENDIF
                !     if (lprnt) print *,' st1=',st1,' tem=',tem,' ror=',ror(l)
                !    &,' qrp=',qrp(l),' rnt=',rnt,' ror1=',ror(l-1),' wvl=',wvl(l)
                !    &,' wvl1=',wvl(l-1),' tem2=',tem2,' vt=',vt(1),' tx3=',tx3
                !
                st2 = tx5
                TEM = ROR(L)*WVL(L) - ROR(L-1)*WVL(L-1)
                IF (tem .GT. 0.0_r8) THEN
                   TX5 = (TX1 - ST1 + TEM2 + TX3)/(1.0_r8+tem*tem1)
                ELSE
                   TX5 = TX1 - tem*tx6 - ST1 + TEM2 + TX3
                ENDIF
                TX5   = MAX(TX5,ZERO)
                tx5 = 0.5_r8 * (tx5 + st2)
                !
                !             qqq = 1.0_r8 + tem * tem1 * (1.0_r8 - sialf)
                !
                !             if (qqq .gt. 0.0_r8) then
                !               TX5   = (TX1 - sialf*tem*tx6 - ST1 + TEM2 + TX3) / qqq
                !             else
                !               TX5   = (TX1 - tem*tx6 - ST1 + TEM2 + TX3)
                !             endif
                !
                !     if(lprnt) print *,' tx51=',tx5,' tx1=',tx1,' st1=',st1,' tem2='
                !     if(tx5 .le. 0.0_r8 .and. l .gt. kd+2)
                !    * print *,' tx51=',tx5,' tx1=',tx1,' st1=',st1,' tem2='
                !    *,tem2,' tx3=',tx3,' tem=',tem,' tem1=',tem1,' wvl=',wvl(l-1),
                !    &wvl(l),' l=',l,' itr=',itr,' evp=',evp(l-1),' vt=',vt(1)
                !    *,' qrp=',qrp(l),' rnt=',rnt,' kd=',kd
                !     if (lprnt) print *,' etd=',etd(l),' wvl=',wvl(l)
                !    &,' trw=',trw(1),trw(2),' ror=',ror(l),' wa=',wa


                !
                TEM1   = ETD(L)
                ETD(L) = ROR(L) * TX5 * MAX(WVL(L),ZERO)
                !
                IF (etd(l) .GT. 0.0_r8) etd(l) = 0.5_r8 * (etd(l) + tem1)
                !

                DEL_ETA = ETD(L) - ETD(L-1)

                !               TEM       = DEL_ETA * TRW(2)
                !               TEM2      = MAX(MIN(TEM, 1.0_r8), -1.0_r8)
                !               IF (ABS(TEM) .GT. 1.0_r8 .AND. ETD(L) .GT. 0.0_r8 ) THEN
                !                 DEL_ETA = TEM2 * TRW(1)
                !                 ETD(L)  = ETD(L-1) + DEL_ETA
                !               ENDIF
                !               IF (WVL(L) .GT. 0.0_r8) TX5 = ETD(L) / (ROR(L)*WVL(L))
                !
                ERRE  = ETD(L) - TEM1
                !
                tem  = MAX(ABS(del_eta), trw(1))
                tem2 = del_eta / tem
                TEM1 = SQRT(MAX((tem+DEL_ETA)*(tem-DEL_ETA),ZERO))
                !               TEM1 = SQRT(MAX((TRW(1)+DEL_ETA)*(TRW(1)-DEL_ETA),0.0_r8))

                EDZ  = (0.5_r8 + ASIN(TEM2)*PIINV)*DEL_ETA + TEM1*PIINV

                DDZ   = EDZ - DEL_ETA
                WCB(L-1) = ETD(L) + DDZ
                !
                TEM1  = HOD(L)
                IF (DEL_ETA .GT. 0.0_r8) THEN
                   QQQ    = 1.0_r8 / (ETD(L) + DDZ)
                   HOD(L) = (ETD(L-1)*HOD(L-1) + DEL_ETA*HOL(L-1)          &
                        &                                            + DDZ*WA(1)) * QQQ
                   QOD(L) = (ETD(L-1)*QOD(L-1) + DEL_ETA*QOL(L-1)          &
                        &                                            + DDZ*WA(2)) * QQQ
                ELSEIF((ETD(L-1) + EDZ) .GT. 0.0_r8) THEN
                   QQQ    = 1.0_r8 / (ETD(L-1) + EDZ)
                   HOD(L) = (ETD(L-1)*HOD(L-1) + EDZ*WA(1)) * QQQ
                   QOD(L) = (ETD(L-1)*QOD(L-1) + EDZ*WA(2)) * QQQ
                ENDIF
                ERRH  = HOD(L) - TEM1
                ERRQ  = ABS(ERRH/HOD(L))  + ABS(ERRE/MAX(ETD(L),ONE_M5))
                !     if (lprnt) print *,' ERRQP=',errq,' errh=',errh,' hod=',hod(l)
                !    &,' erre=',erre,' etd=',etd(l),' del_eta=',del_eta
                DOF   = DDZ
                VT(2) = QQQ

                !
                DDZ  = DOF
                TEM4 = QOD(L)
                TEM1 = VRW(1)
                !
                QHS  = QA(3) + 0.5_r8 * (GAF(L-1)+GAF(L))                    &
                     &                           * (HOD(L)-QA(2))
                !
                !                                           First iteration       !
                !
                ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                TEM2 = ROR(L) * QRP(L)
                CALL QRABF(TEM2,QRAF,QRBF)
                TEM6 = TX5 * (1.6_r8 + 124.9_r8 * QRAF) * QRBF * TX4
                !
                CE   = TEM6 * ST2 / ((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                !
                TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*QOD(L))
                TEM3   = (1.0_r8 + TEM1) * QHS * (QOD(L)+CE)
                TEM    = MAX(TEM2*TEM2 - 4.0_r8*TEM1*TEM3,ZERO)
                QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
                !

                !
                !                                            second iteration   !
                !
                ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                CE   = TEM6 * ST2 / ((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                !             CEE  = CE * (ETD(L)+DDZ)
                !


                TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*tem4)
                TEM3   = (1.0_r8 + TEM1) * QHS * (tem4+CE)
                TEM    = MAX(TEM2*TEM2 - 4.0_r8*TEM1*TEM3,ZERO)
                QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
                !                                              Evaporation in Layer L-1
                !

                EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)
                !                                              Calculate Pd (L+1/2)
                QA(1)    = TX1*RNT + RNF(L-1) - EVP(L-1)
                !
                !     if(lprnt) print *,' etd=',etd(l),' tx5=',tx5,' rnt=',rnt
                !    *,' rnf=',rnf(l-1),' evp=',evp(l-1),' itr=',itr,' L=',L

                !
                IF (qa(1) .GT. 0.0_r8) THEN
                   IF (ETD(L) .GT. 0.0_r8) THEN
                      TEM    = QA(1) / (ETD(L)+ROR(L)*TX5*VT(1))
                      QRP(L) = MAX(TEM,ZERO)
                   ELSEIF (TX5 .GT. 0.0_r8) THEN
                      QRP(L) = (MAX(ZERO,QA(1)/(ROR(L)*TX5*GMS(L))))           &
                           &                                          ** (1.0_r8/1.1364_r8)
                   ELSE
                      QRP(L) = 0.0_r8
                   ENDIF
                ELSE
                   qrp(l) = 0.5_r8 * qrp(l)
                ENDIF
                !                                              Compute Buoyancy
                TEM1   = WA(3)+(HOD(L)-WA(1)-con_hvap*(QOD(L)-WA(2)))         &
                     &                                                  * (1.0_r8/con_cp)
                !             if (lprnt) print *,' tem1=',tem1,' wa3=',wa(3),' hod='
                !    &,hod(l),' wa1=',wa(1),' qod=',qod(l),' wa2=',wa(2),' con_hvap=',con_hvap
                !    &,' cmpor=',cmpor,' dofw=',dofw,' prl=',prl(l),' qrp=',qrp(l)
                TEM1   = TEM1 * (1.0_r8 + con_FVirt*QOD(L))
                ROR(L) = CMPOR * PRL(L) / TEM1
                TEM1   = TEM1 * DOFW
!!!           TEM1   = TEM1 * (1.0_r8 + con_FVirt*QOD(L)) * DOFW

                BUY(L) = (TEM1 - 1.0_r8 - QRP(L)) * ROR(L) * TX5
                !                                              Compute W (L+1/2)

                TEM1   = WVL(L)
                !             IF (ETD(L) .GT. 0.0_r8) THEN
                WVL(L) = VT(2) * (ETD(L-1)*WVL(L-1) - FACG                &
                     &                 * (BUY(L-1)*QRT(L-1)+BUY(L)*QRB(L-1)))
                !
                !             if (lprnt) print *,' wvl=',wvl(l),'vt2=',vt(2),' buy1='
                !    &,buy(l-1),' buy=',buy(l),' qrt1=',qrt(l-1),' qrb1=',qrb(l-1)
                !    &,' etd1=',etd(l-1),' wvl1=',wvl(l-1)
                !             ENDIF
                !
                IF (wvl(l) .LT. 0.0_r8) THEN
                   !               WVL(L) = max(wvl(l), 0.1_r8*tem1)
                   !               WVL(L) = 0.5_r8*tem1
                   !               WVL(L) = 0.1_r8*tem1
                   !               WVL(L) = 0.0_r8
                   WVL(L) = 1.0e-10_r8
                ELSE
                   WVL(L) = 0.5_r8*(WVL(L)+TEM1)
                ENDIF

                !
                !             WVL(L) = max(0.5_r8*(WVL(L)+TEM1), 0.0_r8)

                ERRW   = WVL(L) - TEM1
                !
                ERRQ   = ERRQ + ABS(ERRW/MAX(WVL(L),ONE_M5))

                !     if (lprnt) print *,' errw=',errw,' wvl=',wvl(l)
                !     if(lprnt .or. tx5 .eq. 0.0_r8) then
                !     if(tx5 .eq. 0.0_r8 .and. l .gt. kbl) then
                !        print *,' errq=',errq,' itr=',itr,' l=',l,' wvl=',wvl(l)
                !    &,' tx5=',tx5,' idnm=',idnm,' etd1=',etd(l-1),' etd=',etd(l)
                !    &,' kbl=',kbl
                !     endif
                !
                !     if(lprnt) print *,' itr=',itr,' itrmnd=',itrmnd,' itrmd=',itrmd
                !             IF (ITR .GE. MIN(ITRMIN,ITRMD/2)) THEN
                IF (ITR .GE. MIN(ITRMND,ITRMD/2)) THEN
                   !     if(lprnt) print *,' itr=',itr,' etd1=',etd(l-1),' errq=',errq
                   IF (ETD(L-1) .EQ. 0.0_r8 .AND. ERRQ .GT. 0.2_r8) THEN
                      !     if(lprnt) print *,' bud=',bud(kd),' wa=',wa(1),wa(2)
                      ROR(L)   = BUD(KD)
                      ETD(L)   = 0.0_r8
                      WVL(L)   = 0.0_r8
                      ERRQ     = 0.0_r8
                      HOD(L)   = WA(1)
                      QOD(L)   = WA(2)
                      !                 TX5      = TX1 + TX9
                      IF (L .LE. KBL) THEN
                         TX5      = TX9
                      ELSE
                         TX5 = (STLT(KB1) * QRT(KB1)                         &
                              &                  +  STLT(KBL) * QRB(KB1)) * (0.5_r8*FAC)
                      ENDIF

                      !     if(lprnt) print *,' tx1=',tx1,' rnt=',rnt,' rnf=',rnf(l-1)
                      !    *,' evp=',evp(l-1),' l=',l
                      EVP(L-1) = 0.0_r8
                      TEM      = MAX(TX1*RNT+RNF(L-1),ZERO)
                      QA(1)    = TEM - EVP(L-1)
                      !                 IF (QA(1) .GT. 0.0_r8) THEN
                      !     if(lprnt) print *,' ror=',ror(l),' tx5=',tx5,' tx1=',tx1
                      !    *,' tx9=',tx9,' gms=',gms(l),' qa=',qa(1)
                      !     if(lprnt) call mpi_quit(13)
                      !     if (tx5 .eq. 0.0_r8 .or. gms(l) .eq. 0.0_r8)
                      !     if (lprnt) 
                      !    *  print *,' Atx5=',tx5,' gms=',gms(l),' ror=',ror(l)
                      !    *,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9
                      !    *,' kbl=',kbl,' etd1=',etd(l-1),' idnm=',idnm,' idn=',idn(idnm)
                      !    *,' errq=',errq

                      QRP(L)   = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
                           &                                            ** (1.0_r8/1.1364_r8)
                      !                 endif
                      BUY(L)   = - ROR(L) * TX5 * QRP(L)
                      WCB(L-1) = 0.0_r8
                   ENDIF
                   !
                   DEL_ETA = ETD(L) - ETD(L-1)
                   IF(DEL_ETA .LT. 0.0_r8 .AND. ERRQ .GT. 0.1_r8) THEN
                      ROR(L)   = BUD(KD)
                      ETD(L)   = 0.0_r8
                      WVL(L)   = 0.0_r8
!!!!!             TX5      = TX1 + TX9
                      CLDFRD(L-1) = TX5
                      !
                      DEL_ETA  = - ETD(L-1)
                      EDZ      = 0.0_r8
                      DDZ      = -DEL_ETA
                      WCB(L-1) = DDZ

                      !
                      HOD(L)   = HOD(L-1)
                      QOD(L)   = QOD(L-1)

                      !
                      TEM4     = QOD(L)
                      TEM1     = VRW(1)
                      !
                      QHS      = QA(3) + 0.5_r8 * (GAF(L-1)+GAF(L))            &
                           &                                   * (HOD(L)-QA(2))

                      !
                      !                                           First iteration       !
                      !
                      ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                      TEM2 = ROR(L) * QRP(L-1)
                      CALL QRABF(TEM2,QRAF,QRBF)
                      TEM6 = TX5 * (1.6_r8 + 124.9_r8 * QRAF) * QRBF * TX4
                      !
                      CE   = TEM6*ST2/((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                      !

                      TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*QOD(L))
                      TEM3   = (1.0_r8 + TEM1) * QHS * (QOD(L)+CE)
                      TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                      QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
                      !
                      !                                            second iteration   !
                      !
                      ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                      CE   = TEM6*ST2/((5.4E5_r8*ST2 + 2.55E6_r8)*(ETD(L)+DDZ))
                      !                 CEE  = CE * (ETD(L)+DDZ)
                      !


                      TEM2   = - ((1.0_r8+TEM1)*(QHS+CE) + TEM1*tem4)
                      TEM3   = (1.0_r8 + TEM1) * QHS * (tem4+CE)
                      TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                      QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))

                      !                                              Evaporation in Layer L-1
                      !
                      EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)

                      !                                               Calculate Pd (L+1/2)
                      !                 RNN(L-1) = TX1*RNT + RNF(L-1) - EVP(L-1)
                      QA(1)    = TX1*RNT + RNF(L-1)
                      EVP(L-1) = MIN(EVP(L-1), QA(1))
                      QA(1)    = QA(1) - EVP(L-1)
                      qrp(l)   = 0.0_r8
                      !
                      !     if (tx5 .eq. 0.0_r8 .or. gms(l) .eq. 0.0_r8)
                      !     if (lprnt)
                      !    *  print *,' Btx5=',tx5,' gms=',gms(l),' ror=',ror(l)
                      !    *,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9
                      !    *,' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
                      !    &,' evp=',evp(l-1)
                      !
                      !                 IF (QA(1) .GT. 0.0_r8) THEN
                      !!                  RNS(L-1) = QA(1)
!!!                 tx5      = tx9
                      !                   QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
                      !    &                                         ** (1.0_r8/1.1364_r8)
                      !                 endif
                      !                 ERRQ   = 0.0_r8
                      !                                              Compute Buoyancy
                      !                 TEM1   = WA(3)+(HOD(L)-WA(1)-con_hvap*(QOD(L)-WA(2)))     &
                      !    &                                                  * (1.0_r8/con_cp)
                      !                 TEM1   = TEM1 * (1.0_r8 + con_FVirt*QOD(L)) * DOFW

                      !                 BUY(L) = (TEM1 - 1.0_r8 - QRP(L)) * ROR(L) * TX5
                      !
                      !                 IF (QA(1) .GT. 0.0_r8) RNS(L) = QA(1)
                      IF (L .LE. K) THEN
                         RNS(L) = QA(1)
                         QA(1)  = 0.0_r8
                      ENDIF
                      tx5      = tx9
                      ERRQ     = 0.0_r8
                      QRP(L)   = 0.0_r8
                      BUY(L)   = 0.0_r8

                      !
                   ENDIF
                ENDIF
             ENDIF
             !

          ENDDO                ! End of the iteration loop  for a given L!
          IF (L .LE. K) THEN
             IF (ETD(L-1) .EQ. 0.0_r8                                       &
                  &         .AND. ERRQ .GT. 0.1_r8 .AND. l .LE. kbl) THEN
!!!  &         .AND. ERRQ .GT. ERRMIN*10.0 .and. l .le. kbl) THEN
                !    &         .AND. ERRQ .GT. ERRMIN*10.0) THEN
                ROR(L)   = BUD(KD)
                HOD(L)   = WA(1)
                QOD(L)   = WA(2)
                TX5      =       TX9     ! Does not make too much difference!
                !              TX5      = TX1 + TX9
                EVP(L-1) = 0.0_r8
                !              EVP(L-1) = CEE * (1.0_r8 - qod(l)/qa(3))
                QA(1)    = TX1*RNT + RNF(L-1)
                EVP(L-1) = MIN(EVP(L-1), QA(1))
                QA(1)    = QA(1) - EVP(L-1)

                !              QRP(L)   = 0.0_r8
!                IF (tx5 .EQ. 0.0_r8 .OR. gms(l) .EQ. 0.0_r8) THEN
!                   PRINT *,' Ctx5=',tx5,' gms=',gms(l),' ror=',ror(l)              &
!                        &,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9                     &
!                        &,' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
!                ENDIF
                !              IF (QA(1) .GT. 0.0_r8) THEN
                QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))                 &
                     &                                         ** (1.0_r8/1.1364_r8)
                !              ENDIF
                ETD(L)   = 0.0_r8
                WVL(L)   = 0.0_r8
                ST1      = 1.0_r8 - ALFIND(L)

                ERRQ     = 0.0_r8
                BUY(L)   = - ROR(L) * TX5 * QRP(L)
                WCB(L-1) = 0.0_r8
             ENDIF
          ENDIF
          !
          LL = MIN(IDN(idnm), K+1)
          IF (ERRQ .LT. 1.0_r8 .AND. L .LE. LL) THEN
             IF (ETD(L-1) .GT. 0.0_r8 .AND. ETD(L) .EQ. 0.0_r8) THEN
                IDN(idnm) = L
                wvl(l)    = 0.0_r8
                IF (L .LT. KBL .OR. tx5 .GT. 0.0_r8) idnm  = idnm + 1
                errq      = 0.0_r8
             ENDIF
             IF (etd(l) .EQ. 0.0_r8 .AND. l .GT. kbl) THEN
                idn(idnm) = l
                IF (tx5 .GT. 0.0_r8) idnm  = idnm + 1
             ENDIF
          ENDIF

          !       if (lprnt) then
          !       print *,' ERRQ=',ERRQ,' IDN=',IDN(idnm),' idnm=',idnm
          !       print *,' L=',L,' QRP=',QRP(L),' ETD=',ETD(L),' QA=',QA(1)
          !    *,' evp=',evp(l-1),' rnf=',rnf(l-1)
          !       endif

          ! 
          !     If downdraft properties are not obtainable, (i.e.solution does
          !      not converge) , no downdraft is assumed
          !
          !          IF (ERRQ .GT. ERRMIN*100.0_r8 .AND. IDN(idnm) .EQ. 99)          &
          IF (ERRQ .GT. 0.1_r8 .AND. IDN(idnm) .EQ. 99)                   &
               &                          DDFT = .FALSE.
          !
          !
          DOF = 0.0_r8
          IF (.NOT. DDFT) RETURN
          !
          !     if (ddlgk .or. l .le. idn(idnm)) then
          !     rsum2 = rsum2 + evp(l-1)
          !     print *,' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' qa=',qa(1)
          !    *,' evp=',evp(l-1)
          !     else
          !     rsum1 = rsum1 + rnf(l-1)
          !     print *,' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' rnf=',rnf(l-1)
          !     endif

       ENDDO                      ! End of the L Loop of downdraft !

       TX1 = 0.0_r8

       DOF = QA(1)
       !
       !     print *,' dof=',dof,' rntp=',rntp,' rnb=',rnb
       !     print *,' total=',(rsum1+dof+rntp+rnb)

    ENDIF                       ! SKPDD endif
    !

    RNN(KD) = RNTP
    TX1     = EVP(KD)
    TX2     = RNTP + RNB + DOF

    !     if (lprnt) print *,' tx2=',tx2
    II = IDH
    IF (II .GE. KD1+1) THEN
       RNN(KD)   = RNN(KD) + RNF(KD)
       TX2       = TX2 + RNF(KD)
       RNN(II-1) = 0.0_r8
       TX1       = EVP(II-1)
    ENDIF
    !     if (lprnt) print *,' tx2=',tx2,' idnm=',idnm,' idn=',idn(idnm)
    DO L=KD,K
       II = IDH
       IF (L .GT. KD1 .AND. L .LT. II) THEN
          RNN(L-1) = RNF(L-1)
          TX2      = TX2 + RNN(L-1)
       ELSEIF (L .GE. II .AND. L .LT. IDN(idnm)) THEN
          rnn(l) = rns(l)
          tx2    = tx2 + rnn(l)
          TX1    = TX1 + EVP(L)
       ELSEIF (L .GE. IDN(idnm)) THEN
          ETD(L+1) = 0.0_r8
          HOD(L+1) = 0.0_r8
          QOD(L+1) = 0.0_r8
          EVP(L)   = 0.0_r8
          RNN(L)   = RNF(L) + RNS(L)
          TX2      = TX2    + RNN(L)
       ENDIF
       !     if (lprnt) print *,' tx2=',tx2,' L=',L,' rnn=',rnn(l)
    ENDDO
    !
    !      For Downdraft case the rain is that falls thru the bottom

    L = KBL

    RNN(L)    = RNN(L) + RNB
    CLDFRD(L) = TX5

    !
    !     Caution !! Below is an adjustment to rain flux to maintain
    !                conservation of precip!

    !
    !     if (lprnt) print *,' train=',train,' tx2=',tx2,' tx1=',tx1

    IF (TX1 .GT. 0.0_r8) THEN
       TX1 = (TRAIN - TX2) / TX1
    ELSE
       TX1 = 0.0_r8
    ENDIF

    DO L=KD,K
       EVP(L) = EVP(L) * TX1
    ENDDO
    !
    !***********************************************************************
    !***********************************************************************

    RETURN
  END   SUBROUTINE DDRFT

  !-----------------------------------------------------------------------------------------

  SUBROUTINE QSATCN(TT,P,Q,DQDT)

    !      USE FUNCPHYS , ONLY : fpvs
    IMPLICIT NONE
    !     include 'constant.h'
    !
    REAL(kind=r8), INTENT(IN   ) :: TT
    REAL(KIND=r8), INTENT(IN   ) :: P
    REAL(KIND=r8), INTENT(OUT  ) :: Q 
    REAL(KIND=r8), INTENT(OUT  ) :: DQDT
    !
    REAL(KIND=r8),PARAMETER :: RVI=1.0_r8/con_rv
    REAL(KIND=r8),PARAMETER :: FACW=con_cvap-con_cliq
    REAL(KIND=r8),PARAMETER :: FACI=con_cvap-con_CSOL

    REAL(KIND=r8), PARAMETER :: HSUB=con_hvap+con_hfus
    REAL(KIND=r8), PARAMETER :: tmix=con_ttp-20.0_r8 
    REAL(KIND=r8), PARAMETER :: DEN=1.0_r8/(con_ttp-TMIX)
    REAL(KIND=r8), PARAMETER :: ZERO=0.0_r8
    REAL(KIND=r8), PARAMETER :: ONE=1.0_r8
    REAL(KIND=r8), PARAMETER :: ONE_M10=1.E-10_r8
    !
    !CFPP$ NOCONCUR R
    REAL(kind=r8) es, d, hlorv, W
    DQDT=0.0_r8;es=0.0_r8;d=0.0_r8;hlorv=0.0_r8;w=0.0_r8
    q=0.0_r8 
    !
    !     es    = 10.0_r8 * fpvs(tt)                ! fpvs is in centibars!
    es    = 0.01_r8 * fpvs2es5(tt)                ! fpvs is in Pascals!
    D     = 1.0_r8 / MAX(p+con_epsm1*es,ONE_M10)
    !
    q     = MIN(con_eps*es*D, ONE)
    !
    W     = MAX(ZERO, MIN(ONE, (TT - TMIX)*DEN))
    hlorv = ( W      * (con_hvap + FACW * (tt-con_ttp))                       &
         &       + (1.0_r8-W) * (HSUB + FACI * (tt-con_ttp)) ) * RVI
    dqdt  = p * q * hlorv *  D / (tt*tt)
    !
    RETURN
  END SUBROUTINE QSATCN
  !-----------------------------------------------------------------------------------------

  SUBROUTINE ANGRAD( PRES, ALM,  AL2, TLA)
    !     SUBROUTINE ANGRAD( PRES, ALM, STLA, CTL2, AL2                     &
    !    &,                  PI, TLA, PRB, WFN, UFN)
    !      use module_ras , only : refp, refr, tlac, plac, tlbpl, drdp, almax
    IMPLICIT NONE

    !     real(kind=r8) PRES, STLA, CTL2, pi,  pifac                 &
    REAL(kind=r8), INTENT(IN   ) :: PRES
    REAL(kind=r8), INTENT(IN   ) :: ALM
    REAL(kind=r8), INTENT(OUT  ) :: AL2
    REAL(kind=r8), INTENT(INOUT) :: TLA
    REAL(kind=r8) :: TEM 
    !
    !
    !     pifac = pi / 180.0_r8
    !     print *,' pres=',pres
    IF (TLA .LT. 0.0_r8) THEN
       IF (PRES .LE. PLAC(1)) THEN
          TLA = TLAC(1)
       ELSEIF (PRES .LE. PLAC(2)) THEN
          TLA = TLAC(2) + (PRES-PLAC(2))*tlbpl(1)
       ELSEIF (PRES .LE. PLAC(3)) THEN
          TLA = TLAC(3) + (PRES-PLAC(3))*tlbpl(2)
       ELSEIF (PRES .LE. PLAC(4)) THEN
          TLA = TLAC(4) + (PRES-PLAC(4))*tlbpl(3)
       ELSEIF (PRES .LE. PLAC(5)) THEN
          TLA = TLAC(5) + (PRES-PLAC(5))*tlbpl(4)
       ELSEIF (PRES .LE. PLAC(6)) THEN
          TLA = TLAC(6) + (PRES-PLAC(6))*tlbpl(5)
       ELSEIF (PRES .LE. PLAC(7)) THEN
          TLA = TLAC(7) + (PRES-PLAC(7))*tlbpl(6)
       ELSEIF (PRES .LE. PLAC(8)) THEN
          TLA = TLAC(8) + (PRES-PLAC(8))*tlbpl(7)
       ELSE
          TLA = TLAC(8)
       ENDIF
       !         tla = tla * 1.5

       !         STLA = SIN(TLA*PIFAC)
       !         TEM1 = COS(TLA*PIFAC)
       !         CTL2 = TEM1 * TEM1

    ELSE
       !         STLA = SIN(TLA*PIFAC)
       !         TEM1 = COS(TLA*PIFAC)
       !         CTL2 = TEM1 * TEM1

    ENDIF
    IF (PRES .GE. REFP(1)) THEN
       TEM = REFR(1)
    ELSEIF (PRES .GE. REFP(2)) THEN
       TEM = REFR(1) + (PRES-REFP(1)) * drdp(1)
    ELSEIF (PRES .GE. REFP(3)) THEN
       TEM = REFR(2) + (PRES-REFP(2)) * drdp(2)
    ELSEIF (PRES .GE. REFP(4)) THEN
       TEM = REFR(3) + (PRES-REFP(3)) * drdp(3)
    ELSEIF (PRES .GE. REFP(5)) THEN
       TEM = REFR(4) + (PRES-REFP(4)) * drdp(4)
    ELSEIF (PRES .GE. REFP(6)) THEN
       TEM = REFR(5) + (PRES-REFP(5)) * drdp(5)
    ELSE
       TEM = REFR(6)
    ENDIF
    !!      AL2 = min(ALMAX, MAX(ALM, 2.0E-4/TEM))
    !       AL2 = min(2.0E-3, MAX(ALM, 2.0E-4/TEM))
    !
    tem = 2.0E-4_r8 / tem
    al2 = MIN(4.0_r8*tem, MAX(alm, tem))
    !
    RETURN
  END SUBROUTINE ANGRAD
  !-----------------------------------------------------------------------------------------

  SUBROUTINE SETQRP()
    ! use module_ras , only : NQRP,C1XQRP,C2XQRP,TBQRP,TBQRA,TBQRB
    IMPLICIT NONE

    REAL(kind=r8) :: tem2,tem1,x,xinc,xmax,xmin
    INTEGER :: jx
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !CFPP$ NOCONCUR R
    !     XMIN=1.0E-6
    XMIN=0.0_r8
    XMAX=5.0_r8
    XINC=(XMAX-XMIN)/(NQRP-1)
    C2XQRP=1.0_r8/XINC
    C1XQRP=1.0_r8-XMIN*C2XQRP
    TEM1 = 0.001_r8 ** 0.2046_r8
    TEM2 = 0.001_r8 ** 0.525_r8
    DO JX=1,NQRP
       X         = XMIN + (JX-1)*XINC
       TBQRP(JX) =        X ** 0.1364_r8
       TBQRA(JX) = TEM1 * X ** 0.2046_r8
       TBQRB(JX) = TEM2 * X ** 0.525_r8
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END  SUBROUTINE SETQRP

  !-----------------------------------------------------------------------------------------

  REAL(kind=r8) FUNCTION QRPF(QRP)
    !
    !      use module_ras , only : NQRP,C1XQRP,C2XQRP,TBQRP,TBQRA,TBQRB
    IMPLICIT NONE

    REAL(kind=r8) QRP, XJ, REAL_NQRP
    REAL(kind=r8), PARAMETER :: ONE=1.0_r8
    INTEGER JX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    REAL_NQRP=REAL(NQRP)
    XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),REAL_NQRP)
    !     XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),FLOAT(NQRP))
    JX   = INT(MIN(XJ,NQRP-ONE),KIND=i4)
    QRPF = TBQRP(JX)  + (XJ-JX) * (TBQRP(JX+1)-TBQRP(JX))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END FUNCTION QRPF
  !-----------------------------------------------------------------------------------------

  SUBROUTINE QRABF(QRP,QRAF,QRBF)
    !      use module_ras , only : NQRP,C1XQRP,C2XQRP,TBQRP,TBQRA,TBQRB
    IMPLICIT NONE
    !
    REAL(kind=r8) QRP, QRAF, QRBF, XJ, REAL_NQRP
    REAL(kind=r8), PARAMETER :: ONE=1.0_r8
    INTEGER JX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    REAL_NQRP=REAL(NQRP)
    XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),REAL_NQRP)
    JX   = INT(MIN(XJ,NQRP-ONE),KIND=i4)
    XJ   = XJ - JX
    QRAF = TBQRA(JX)  + XJ * (TBQRA(JX+1)-TBQRA(JX))
    QRBF = TBQRB(JX)  + XJ * (TBQRB(JX+1)-TBQRB(JX))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END SUBROUTINE QRABF

  !-----------------------------------------------------------------------------------------

  SUBROUTINE SETVTP()
    !use module_ras , only : NVTP,C1XVTP,C2XVTP,TBVTP
    IMPLICIT NONE

    REAL(kind=r8) :: xinc,x,xmax,xmin
    INTEGER :: jx
    REAL(kind=r8), PARAMETER :: VTPEXP=-0.3636_r8
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !CFPP$ NOCONCUR R
    XMIN=0.05_r8
    XMAX=1.5_r8
    XINC=(XMAX-XMIN)/(NVTP-1)
    C2XVTP=1.0_r8/XINC
    C1XVTP=1.0_r8-XMIN*C2XVTP
    DO JX=1,NVTP
       X         = XMIN + (JX-1)*XINC
       TBVTP(JX) =        X ** VTPEXP
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END  SUBROUTINE SETVTP


  !-----------------------------------------------------------------------------------------

 REAL(kind=r8)  FUNCTION VTPF(ROR)
    !
    !use module_ras , only : NVTP,C1XVTP,C2XVTP,TBVTP
    IMPLICIT NONE
    REAL(kind=r8) :: ROR, XJ, REAL_NVTP
    REAL(kind=r8), PARAMETER :: ONE=1.0_r8
    INTEGER :: JX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    REAL_NVTP=REAL(NVTP)
    XJ   = MIN(MAX(C1XVTP+C2XVTP*ROR,ONE),REAL_NVTP)
    JX   = INT(MIN(XJ,NVTP-ONE),KIND=i4)
    VTPF = TBVTP(JX)  + (XJ-JX) * (TBVTP(JX+1)-TBVTP(JX))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RETURN
  END  FUNCTION VTPF
  !-----------------------------------------------------------------------------------------


    REAL(kind=r8) FUNCTION CLF(PRATE)
    !
    IMPLICIT NONE
    REAL(kind=r8) PRATE
    !
    REAL (kind=r8), PARAMETER :: ccf1=0.30_r8, ccf2=0.09_r8                 &
         &,                                   ccf3=0.04_r8, ccf4=0.01_r8          &
         &,                                   pr1=1.0_r8,   pr2=5.0_r8            &
         &,                                   pr3=20.0_r8
    !
    IF (prate .LT. pr1) THEN
       clf = ccf1
    ELSEIF (prate .LT. pr2) THEN
       clf = ccf2
    ELSEIF (prate .LT. pr3) THEN
       clf = ccf3
    ELSE
       clf = ccf4
    ENDIF
    !
    RETURN
  END FUNCTION CLF



  !-----------------------------------------------------------------------------------------
  SUBROUTINE Finalize_Cu_RAS3PHASE()
    IMPLICIT NONE
    DEALLOCATE (rasal)

  END SUBROUTINE Finalize_Cu_RAS3PHASE
  !-----------------------------------------------------------------------------------------

END MODULE Cu_RAS3PHASE


!PROGRAM Main
! USE Cu_RAS3PHASE, ONLY: Init_Cu_RAS3PHASE
! IMPLICIT NONE
! INTEGER :: kMax=28
! REAL(KIND=8) :: dt=1200.0
! REAL(KIND=8) :: fhour=0.0
! INTEGER :: idate(1:4)=(/00,01,29,2015/)
! INTEGER :: iMax=2
! INTEGER :: jMax=1
!! INTEGER :: ibMax=2
! INTEGER :: jbMax=1
! REAL(kind=8) ::  si_in(kMax+1)
! REAL(kind=8) ::  sl_in(kMax)
!
! CALL Init_Cu_RAS3PHASE(kMax,dt,fhour,idate,iMax,jMax,ibMax,jbMax,si_in,sl_in)
! PRINT*,' --  '
! CALL Init_Cu_RAS3PHASE(kMax,dt,fhour,idate,iMax,jMax,ibMax,jbMax,si_in,sl_in)
!
!END PROGRAM Main
