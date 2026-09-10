MODULE Pbl_HostlagBoville
  USE Constants, ONLY :     &
       cp,            &
       grav,          &
       gasr,          &
       r8,i4,i8

 USE PBL_Entrain, ONLY :  &
       PBLNASA

 USE Parallelism, ONLY: &
       MsgOne, FatalError

  USE Options, ONLY : &
      jdt        ,&
      nferr      , &!
      nfprt      , &!
      PBLEntrain , &!
      Wgh1       , &!1.0_r8/3.0_r8 ! pbl Hostlag Boville
      Wgh2       , &!1.0_r8/3.0_r8 ! pbl Mellor Yamada 2.0
      Wgh3         !1.0_r8/3.0_r8 ! pbl Mellor Yamada 2.5 

  IMPLICIT NONE
  SAVE
  PRIVATE
  ! Mellor & Yamada 1982, eq (45)
  ! Previous values from Mellor (1973) are shown as comments
  REAL(KIND=r8), PARAMETER :: epsq2  = 0.2_r8
  REAL(KIND=r8), PARAMETER :: FH     = 1.01_r8
  REAL(KIND=r8), PARAMETER :: onet   = 1.0_r8/3.0_r8   ! 1/3 power in wind gradient expression  common /compbl/
  REAL(KIND=r8), PARAMETER :: fak    = 8.5_r8     ! Constant in surface temperature excess      common /compbl/    
  REAL(KIND=r8), PARAMETER :: sffrac = 0.01_r8     ! Surface layer fraction of boundary layer common /compbl/
  REAL(KIND=r8), PARAMETER :: vk     = 0.40_r8    ! Von Karman's constant common /compbl/
  REAL(KIND=r8), PARAMETER :: fakn   = 7.2_r8     ! Constant in turbulent prandtl number
  REAL(KIND=r8), PARAMETER :: ccon   = fak*sffrac*vk  ! fak * sffrac * vk  common /compbl/

  REAL(KIND=r8), PARAMETER :: betam  = 15.0_r8 ! Constant in wind gradient expression common /compbl/
  REAL(KIND=r8), PARAMETER :: betas  =  5.0_r8 ! Constant in surface layer gradient expression common /compbl/
  REAL(KIND=r8), PARAMETER :: betah  = 15.0_r8 ! Constant in temperature gradient expression  common /compbl/
  !REAL(KIND=r8), PARAMETER :: ricr   =  0.1_r8  ! Critical richardson number
  REAL(KIND=r8) :: ricr=  0.1_r8 
    !
    ! sffrac = 0.1_r8     ! Surface layer fraction of boundary layer common /compbl/
    ! 
    ! betam  = 15.0_r8 ! Constant in wind gradient expression common /compbl/
    !
  REAL(KIND=r8), PARAMETER :: binm  = betam*sffrac  ! betam * sffrac  common /compbl/	binm = betam*sffrac
  REAL(KIND=r8), PARAMETER :: binh  = betah*sffrac  ! betah * sffrac common /compbl/   binh = betah*sffrac

  REAL(KIND=r8), PARAMETER :: cpair  = 1004.64_r8 !heat capacity dry air at const pres (j/kg/kelvin)
  REAL(KIND=r8), PARAMETER :: rair   = 287.04_r8  !gas constant for dry air (j/kg/kelvin)
  REAL(KIND=r8), PARAMETER :: gravit = 9.80616_r8 !Acceleration due to gravity common/comvd/
  REAL(KIND=r8), PARAMETER :: g= gravit      ! Gravitational acceleration  common /compbl/
    !
    !-----------------------------------------------------------------------
    !
    ! Hard-wired numbers.
    ! zkmin = minimum k = kneutral*f(ri)
    !
  REAL(KIND=r8), PARAMETER    :: zkmin   = 0.01_r8 ! Minimum kneutral*f(ri) common/comvd/
  INTEGER      , PARAMETER    :: ntopfl=1  ! Top level to which vertical diffusion is applied. common/comvd/
  INTEGER                     :: npbl      ! Maximum number of levels in pbl from surface common/comvd/

  REAL(KIND=r8), ALLOCATABLE  :: ml2    (:)    ! Mixing lengths squared common/comvd/
  REAL(KIND=r8), ALLOCATABLE  :: hypm   (:)    ! reference pressures at midpoints
  REAL(KIND=r8), ALLOCATABLE  :: qmin   (:)    ! Global minimum constituent concentration
  REAL(KIND=r8), ALLOCATABLE  :: qmincg (:)    ! Min. constituent concentration counter-gradient term

  REAL(KIND=r8), ALLOCATABLE :: TKEMYJ_MY20(:,:,:) 
  REAL(KIND=r8), ALLOCATABLE :: TKEMYJ_MY25(:,:,:) 
  REAL(KIND=r8), ALLOCATABLE :: TKEMYJ_HSBO(:,:,:)  

  REAL(KIND=r8),    PARAMETER :: CPWV   =  1.810e3_r8 ! specific heat of water vap
  REAL(KIND=r8),    PARAMETER :: cpwvx  = CPWV      ! spec. heat of water vapor at const. pressure
  REAL(KIND=r8),    PARAMETER :: cpairx = cpair      ! specific heat of dry air
  REAL(KIND=r8),    PARAMETER :: cpvir  = cpwvx/cpairx - 1.0_r8 ! Derived constant for cp moist air common/comvd/
  REAL(KIND=r8),    PARAMETER :: gkm0   = 1.00_r8
  REAL(KIND=r8),    PARAMETER :: gkh0   = 0.10_r8
  REAL(KIND=r8),    PARAMETER :: gkm1   = 300.0_r8
  REAL(KIND=r8),    PARAMETER :: gkh1   = 300.0_r8
  INTEGER      ,    PARAMETER :: kmean  =   1
  REAL(KIND=r8),    PARAMETER :: facl   =   0.05_r8
  REAL(KIND=r8),    PARAMETER :: a1     =   0.92_r8
  REAL(KIND=r8),    PARAMETER :: a2     =   0.74_r8
  REAL(KIND=r8),    PARAMETER :: b1     =  16.6_r8
  REAL(KIND=r8),    PARAMETER :: b2     =  10.1_r8
  REAL(KIND=r8),    PARAMETER :: c1     =   0.08_r8
  REAL(KIND=r8),    PARAMETER :: deltx  =   0.0_r8
  REAL(KIND=r8),    PARAMETER :: eps    =   0.608_r8
  !  gbyr        =(grav/gasr)**2!(m/sec**2)/(J/(Kg*K))=(m/sec**2)/((Kg*(m/sec**2)*m)/(Kg*K))
    !(m/sec**2)/((Kg*(m**2/sec**2))/(Kg*K))
    !(m/sec**2)/(m**2/sec**2*K)=K**2/m**2

  REAL(KIND=r8),    PARAMETER :: gbyr= (gravit/gasr)**2!(m/sec**2)/(J/(Kg*K))=(m/sec**2)/((Kg*(m/sec**2)*m)/(Kg*K))
    !(m/sec**2)/((Kg*(m**2/sec**2))/(Kg*K))
    !(m/sec**2)/(m**2/sec**2*K)=K**2/m**2
  INTEGER     , PARAMETER     :: nitr  =   2
  REAL(KIND=r8)               :: alfa
  REAL(KIND=r8)               :: beta
  REAL(KIND=r8)               :: gama
  REAL(KIND=r8)               :: dela
  REAL(KIND=r8)               :: r1
  REAL(KIND=r8)               :: r2
  REAL(KIND=r8)               :: r3
  REAL(KIND=r8)               :: r4
  REAL(KIND=r8)               :: s1
  REAL(KIND=r8)               :: s2
  REAL(KIND=r8)               :: rfc

  REAL(KIND=r8),PARAMETER     :: VKARMAN=0.4_r8
  !---------------------------------------------------------------------
  ! --- CONSTANTS
  !---------------------------------------------------------------------

  REAL(KIND=r8) :: aa1,   aa2,   bb1,  bb2,  ccc
  REAL(KIND=r8) :: ckm1,  ckm2,  ckm3, ckm4, ckm5, ckm6, ckm7, ckm8
  REAL(KIND=r8) :: ckh1,  ckh2,  ckh3, ckh4
  REAL(KIND=r8) :: cvfq1, cvfq2, bcq
    
  REAL(KIND=r8), PARAMETER :: aa1_old =  0.78_r8
  REAL(KIND=r8), PARAMETER :: aa2_old =  0.79_r8
  REAL(KIND=r8), PARAMETER :: bb1_old = 15.0_r8
  REAL(KIND=r8), PARAMETER :: bb2_old =  8.0_r8
  REAL(KIND=r8), PARAMETER :: ccc_old =  0.056_r8
  
  REAL(KIND=r8), PARAMETER :: aa1_new =  0.92_r8
  REAL(KIND=r8), PARAMETER :: aa2_new =  0.74_r8
  REAL(KIND=r8), PARAMETER :: bb1_new = 16.0_r8
  REAL(KIND=r8), PARAMETER :: bb2_new = 10.0_r8
  REAL(KIND=r8), PARAMETER :: ccc_new =  0.08_r8
  
  REAL(KIND=r8), PARAMETER :: cc1     =  0.27_r8
  REAL(KIND=r8), PARAMETER :: t00     =  2.7248e2_r8 
  REAL(KIND=r8), PARAMETER :: small2   =  1.0e-10_r8
  REAL(KIND=r8), PARAMETER :: VONKARM = 0.40_r8     

  !---------------------------------------------------------------------
  ! --- NAMELIST
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER    :: TKEmax       =   5.0_r8
  REAL(KIND=r8), PARAMETER    :: TKEmin       =   1.0e-12_r8
  REAL(KIND=r8), PARAMETER    :: el0max1      =   1.0e6_r8
  REAL(KIND=r8), PARAMETER    :: el0min1      =   0.0_r8
  REAL(KIND=r8), PARAMETER    :: alpha_land   =   0.10_r8
  REAL(KIND=r8), PARAMETER    :: alpha_sea    =   0.10_r8
  REAL(KIND=r8), PARAMETER    :: akmax        = 300.0_r8
  REAL(KIND=r8), PARAMETER    :: akmin_land   =   5.0_r8
!  REAL(KIND=r8), PARAMETER    :: akmin_sea    =   0.5_r8
  REAL(KIND=r8), PARAMETER    :: akmin_sea    =   5.0_r8

  INTEGER, PARAMETER :: nk_lim       =  2
  INTEGER, PARAMETER :: init_iters   = 1! 20
  LOGICAL, PARAMETER :: do_thv_stab  = .TRUE.
  LOGICAL, PARAMETER :: use_old_cons = .FALSE.
  LOGICAL, PARAMETER :: ensamble     = .TRUE.

  REAL(KIND=r8), PARAMETER    :: kcrit        =  0.01_r8

    REAL(KIND=r8) :: rairx        ! gas constant for dry air
    REAL(KIND=r8) :: akappa


  PUBLIC :: vdinti
  PUBLIC :: vdintr
CONTAINS

  SUBROUTINE vdinti(ibMax,jbMax,plev,a_hybr,b_hybr,pnats,pcnst)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Initialization of time independent fields for vertical diffusion.
    ! 
    ! Method: 
    ! Call initialization routine for boundary layer scheme.
    ! 
    ! Author: J. Rosinski
    !
    !-----------------------------------------------------------------------
    !   use precision
    !   use pmgrid

    IMPLICIT NONE

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER,  INTENT(in) :: ibMax
    INTEGER,  INTENT(in) :: jbMax
    INTEGER,  INTENT(in) :: plev         ! number of vertical levels
    REAL(KIND=r8),    INTENT(IN) :: a_hybr(plev+1)
    REAL(KIND=r8),    INTENT(IN) :: b_hybr(plev+1)
    INTEGER,OPTIONAL, INTENT(IN) :: pnats      ! number of non-advected trace species
    INTEGER,OPTIONAL, INTENT(IN) :: pcnst      ! number of constituents (including water vapor)

    REAL(KIND=r8) :: gam1
    REAL(KIND=r8) :: gam2

    REAL(KIND=r8) :: ps0    
    !REAL(KIND=r8) :: sigk(plev)

    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER :: k         ! vertical loop index 
    INTEGER :: l         ! vertical loop index reversed 
    INTEGER :: m
    CHARACTER(LEN=*), PARAMETER :: tab="    "
    CHARACTER(LEN=500)          :: line
    CHARACTER(LEN=*), PARAMETER :: h="**(Pbl_HostlagBoville)**"

    !
    ALLOCATE(hypm   (plev))               ;hypm  =0.0_r8
    ALLOCATE(ml2    (plev + 1))           ;ml2   =0.0_r8
    ALLOCATE(qmin   (pcnst+pnats))        ;qmin  =0.0_r8
    ALLOCATE(qmincg (pcnst+pnats))        ;qmincg=0.0_r8 
    ALLOCATE(TKEMYJ_MY20(ibMax,plev+1,jbMax) ); TKEMYJ_MY20=0.0_r8
    ALLOCATE(TKEMYJ_MY25(ibMax,plev+1,jbMax) ); TKEMYJ_MY25=0.0_r8
    ALLOCATE(TKEMYJ_HSBO(ibMax,plev+1,jbMax) ); TKEMYJ_HSBO=0.0_r8

    gam1=1.0_r8/3.0_r8-2.0_r8*a1/b1
    gam2=(b2+6.0_r8*a1)/b1
    alfa=b1*(gam1-c1)+3.0_r8*(a2+2.0_r8*a1) 
    beta=b1*(gam1-c1)
    gama=a2/a1*(b1*(gam1+gam2)-3.0_r8*a1)
    dela=a2/a1* b1* gam1
    r1  =0.5_r8*gama/alfa
    r2  =    beta/gama
    r3  =2.0_r8*(2.0_r8*alfa*dela-gama*beta)/(gama*gama)
    r4  =r2*r2
    s1  =3.0_r8*a2* gam1
    s2  =3.0_r8*a2*(gam1+gam2)
    !     
    !     critical flux richardson number
    !     
    rfc =s1/s2
    ricr=rfc
    ! (Use ricr = 0.3 in this formulation)
    akappa=gasr/cp

    rairx  =rair
      

   ! hypm     reference state midpoint pressures
    ps0    = 1.0e5_r8            ! Base state surface pressure (pascals)
    DO k=plev,1,-1
!  SB
!      hypm(k) =  ps0*sig(k)
       l = plev+1-k  ! inversion from top to bottom to bottom to top
       hypm(k) =  0.5_r8 * (ps0*(b_hybr(l)+b_hybr(l+1)) + a_hybr(l) + &
                            a_hybr(l+1) )
    END DO
    qmin(1) = 1.0e-21_r8        ! Minimum mixing ratio for moisture
    DO m=2,pcnst
       qmin(m) = 1.02-21_r8
    END DO
    qmincg(1) = 1.e-21_r8
    DO m=2,pcnst
       qmincg(m) = qmin(m)
    END DO
    !
    !-----------------------------------------------------------------------
    !
    ! Hard-wired numbers.
    ! zkmin = minimum k = kneutral*f(ri)
    !
    ! zkmin = 0.01_r8
    !
    ! Set physical constants for vertical diffusion and pbl:
    !
    !
    ! Derived constants
    ! ntopfl = top level to which v-diff is applied
    ! npbl = max number of levels (from bottom) in pbl
    !
    !
    ! Limit pbl height to regions below 400 mb
    !
    npbl = 0
    DO k=plev,1,-1
       IF (hypm(k) >= 4.e4_r8) THEN!40000
          npbl = npbl + 1
       END IF
    END DO
    npbl = MAX(npbl,1)

    WRITE(line,'(A44,I5,A22,F12.5,A8,F12.5)')'VDINTI: PBL height will be limited to bottom',npbl, &
         ' model levels. Top is ',hypm(plev + 1-npbl),' pascals',ricr
    WRITE(line,'(A44)')'VDINTI: PBL height will be limited to bottom'
    CALL MsgOne(h,TRIM(line))

    !ntopfl = 1
    !IF (plev.EQ.1) ntopfl = 0
    !
    ! Set the square of the mixing lengths.
    !
    ml2(1) = 0.0_r8
    DO k=2,plev
       ml2(k) = 10.0_r8**2
    END DO
    ml2(plev + 1) = 0.0_r8
    !
    ! Initialize pbl variables
    !
    !CALL pbinti()
    !
    !---------------------------------------------------------------------
    ! --- Initialize constants
    !---------------------------------------------------------------------

    IF( use_old_cons ) THEN
       aa1 = aa1_old 
       aa2 = aa2_old  
       bb1 = bb1_old 
       bb2 = bb2_old  
       ccc = ccc_old  
    ELSE
       aa1 = aa1_new 
       aa2 = aa2_new  
       bb1 = bb1_new 
       bb2 = bb2_new  
       ccc = ccc_new  
    END IF

    ckm1 = ( 1.0_r8 - 3.0_r8*ccc )*aa1
    ckm3 =  3.0_r8 * aa1*aa2*    ( bb2 - 3.0_r8*aa2 )
    ckm4 =  9.0_r8 * aa1*aa2*ccc*( bb2 + 4.0_r8*aa1 )
    ckm5 =  6.0_r8 * aa1*aa1
    ckm6 = 18.0_r8 * aa1*aa1*aa2*( bb2 - 3.0_r8*aa2 )
    ckm7 =  3.0_r8 * aa2*        ( bb2 + 7.0_r8*aa1 )
    ckm8 = 27.0_r8 * aa1*aa2*aa2*( bb2 + 4.0_r8*aa1 )
    ckm2 =  ckm3 - ckm4
    ckh1 =  aa2
    ckh2 =  6.0_r8 * aa1*aa2
    ckh3 =  3.0_r8 * aa2*( bb2 + 4.0_r8*aa1 )
    ckh4 =  2.0e-6_r8 * aa2
    cvfq1 = 5.0_r8 * cc1 / 3.0_r8
    cvfq2 = 1.0_r8 / bb1
    bcq   = 0.5_r8 * ( bb1**(2.0_r8/3.0_r8) )

    !---------------------------------------------------------------------


    RETURN
  END SUBROUTINE vdinti



  SUBROUTINE pbinti()
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Initialize time independent variables of pbl package.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: B. Boville
    ! 
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    !
    ! Basic constants
    !
    !
    ! Derived constants
    !
    !ccon = fak*sffrac*vk
    !
    ! sffrac = 0.1_r8     ! Surface layer fraction of boundary layer common /compbl/
    ! 
    ! betam  = 15.0_r8 ! Constant in wind gradient expression common /compbl/
    !
    !binm = betam*sffrac
    !binh = betah*sffrac

    !
    RETURN
  END SUBROUTINE pbinti

  SUBROUTINE vdintr(&
                    latco       , &! INTENT(IN        ) plon ! number of longitudes
                    plon       , &! INTENT(IN        ) plon ! number of longitudes
                    plond      , &! INTENT(IN        ) plond ! slt extended domain longitude
                    plev       , &! INTENT(IN        ) plev ! number of vertical levels
                    pcnst      , &! INTENT(IN        ) pcnst ! number of constituents (including water vapor)
                    ztodt      , &! INTENT(IN        ) ztodt ! 2 delta-t
                    colrad     , &! INTENT(IN        ) Cosino the colatitude [radian]
                    gl0        , &! INTENT(IN        ) Maximum mixing length l0 in blackerdar's formula [m]
                    bstar      , &! INTENT(IN        ) surface_bouyancy_scale m s-2
                    FRLAND     , &! INTENT(IN        ) Fraction Land [%]
                    z0         , &! INTENT(IN        ) Rougosiness [m]
                    LwCoolRate , &! INTENT(IN        ) air_temperature_tendency_due_to_longwave [K s-1]
                    LwCoolRateC, &! INTENT(IN        ) clear sky_air_temperature_tendency_lw [K s-1]
                    cldtot     , &! INTENT(IN        ) cloud fraction [%]
                    qliq       , &! INTENT(IN        ) cloud fraction [kg/kg]
                    pmidm1     , &! INTENT(IN        ) pmidm1(plond,plev) ! midpoint pressures
                    pintm1     , &! INTENT(IN        ) pintm1(plond,plev + 1)    ! interface pressures
                    psomc      , &! INTENT(IN        ) psomc(plond,plev)         ! (psm1/pmidm1)**cappa
                    thm        , &! INTENT(IN        ) thm(plond,plev) ! potential temperature midpoints
                    zm         , &! INTENT(IN        ) zm(plond,plev) ! midpoint geopotential height above sfc
                    zhalf      , &! INTENT(IN        ) zhalf(plond,plev) ! interface pressures geopotential height
                    psfcpa     , &! INTENT(IN        ) surface pressure [mb]
                    USTAR      , &! INTENT(IN        ) scale velocity turbulent [m/s]
                    TSK        , &! INTENT(IN        ) surface temperature
                    QSFC       , &! INTENT(IN        ) surface specific temperature
                    rpdel      , &! INTENT(IN        ) rpdel(plond,plev)! 1./pdel (thickness between interfaces)
                    rpdeli     , &! INTENT(IN        ) rpdeli(plond,plev)! 1./pdeli (thickness between midpoints)
                    um1        , &! INTENT(IN        ) um1(plond,plev)         ! u-wind input
                    vm1        , &! INTENT(IN        ) vm1(plond,plev)         ! v-wind input
                    tm1        , &! INTENT(IN        ) tm1(plond,plev)         ! temperature input
                    taux       , &! INTENT(IN        ) taux(plond)                   ! x surface stress ![N/m**2]
                    tauy       , &! INTENT(IN        ) tauy(plond)                   ! y surface stress ![N/m**2]
                    shflx      , &! INTENT(IN        ) shflx(plond)! surface sensible heat flux (w/m2)
                    cflx       , &! INTENT(IN        ) cflx(plond,pcnst)! surface constituent flux (kg/m2/s)
                    qm1        , &! INTENT(IN    ) qm1(plond,plev,pcnst)  ! initial/final constituent field
                    dtv        , &! INTENT(OUT  ) dtv(plond,plev)     ! temperature tendency (heating)
                    dqv        , &! INTENT(OUT  ) dqv(plond,plev,pcnst)
                    duv        , &! INTENT(OUT  ) duv(plond,plev)     ! u-wind tendency
                    dvv        , &! INTENT(OUT  ) dvv(plond,plev)     ! v-wind tendency
                    up1        , &! INTENT(OUT  ) up1(plond,plev)     ! u-wind after vertical diffusion
                    vp1        , &! INTENT(OUT  ) vp1(plond,plev)     ! v-wind after vertical diffusion
                    pblh       , &! INTENT(OUT  ) pblh(plond)! planetary boundary layer height
                    rino       , &! INTENT(INOUT) bulk Richardson no. from level to ref lev
                    tpert      , &! INTENT(OUT  ) tpert(plond)! convective temperature excess
                    qpert      , &! INTENT(OUT  ) qpert(plond,pcnst)! convective humidity and constituent excess
                    TKE        , &! INTENT(INOUT) Turbulent kinetic Energy [m/s]^2
                    kvh        , &! INTENT(OUT  ) Heat Coeficient Difusivity 
                    kvm        , &! INTENT(OUT  ) Momentun Coeficient Difusivity 
                    obklen     , &! INTENT(OUT  ) Heat Coeficient Difusivity 
                    phiminv    , &! INTENT(OUT  ) Momentum Stability Function 
                    phihinv    , &! INTENT(OUT  ) Heat Stability Function 
                    tstar      , &
                    wstar        )
    !-----------------------------------------------------------------------
    !
    ! interface routine for vertical diffusion and pbl scheme
    !
    ! calling sequence:
    !
    !    vdinti        initializes vertical diffustion constants
    !    pbinti        initializes pbl constants
    !     .
    !     .
    !    vdintr        interface for vertical diffusion and pbl scheme
    !      vdiff       performs vert diff and pbl
    !        pbldif    boundary layer scheme
    !        mvdiff    diffuse momentum
    !        qvdiff    diffuse constituents
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  J. Rosinski
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          P. Rasch, B. Boville, August 1992
    ! Reviewed:          P. Rasch, April 1996
    ! Reviewed:          B. Boville, April 1996
    !
    !-----------------------------------------------------------------------
    !
    ! $Id: vdintr.F,v 1.1.1.1 2001/03/09 00:29:29 mirin Exp $
    !
    !-----------------------------------------------------------------------
    !!      use precision
    !!      use pmgrid
    !-----------------------------------------------------------------------
    !!#include <implicit.h>
    !------------------------------Commons----------------------------------
    !!#include <comtrcnm.h>
    !-----------------------------------------------------------------------
    !!#include <comvd.h>
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) ::  latco
    INTEGER, INTENT(IN   ) ::  plon       ! number of longitudes
    INTEGER, INTENT(IN   ) ::  plev       ! number of vertical levels
    INTEGER, INTENT(IN   ) ::  plond      ! slt extended domain longitude
    INTEGER, INTENT(IN   ) ::  pcnst      ! number of constituents (including water vapor)
    REAL(kind=r8), INTENT(INOUT) ::  gl0        (plond)     !
    REAL(kind=r8), INTENT(IN   ) ::  bstar      (plond)     !
    REAL(kind=r8), INTENT(IN   ) ::  FRLAND     (plond)     !
    REAL(kind=r8), INTENT(IN   ) ::  LwCoolRate (plond,plev)     !
    REAL(kind=r8), INTENT(IN   ) ::  LwCoolRateC(plond,plev)     !
    REAL(kind=r8), INTENT(IN   ) ::  cldtot     (plond,plev)     !
    REAL(kind=r8), INTENT(IN   ) ::  qliq       (plond,plev)     !
    REAL(kind=r8), INTENT(IN   ) ::  pmidm1     (plond,plev)     ! midpoint pressures
    REAL(kind=r8), INTENT(IN   ) ::  pintm1     (plond,plev + 1) ! interface pressures
    REAL(kind=r8), INTENT(IN   ) ::  psomc      (plond,plev)     ! (psm1/pmidm1)**cappa
    REAL(kind=r8), INTENT(IN   ) ::  thm        (plond,plev)     ! potential temperature midpoints
    REAL(kind=r8), INTENT(IN   ) ::  zm         (plond,plev)     ! midpoint geopotential height above sfc
    REAL(kind=r8), INTENT(IN   ) ::  zhalf      (plond,plev+1) 
    REAL(kind=r8), INTENT(IN   ) ::  psfcpa     (plond)
    REAL(kind=r8), INTENT(INOUT) ::  USTAR      (plond)
    REAL(kind=r8), INTENT(IN   ) ::  TSK        (plond)
    REAL(kind=r8), INTENT(IN   ) ::  QSFC       (plond)
    REAL(kind=r8), INTENT(IN   ) ::  z0         (plond)
    REAL(kind=r8), INTENT(IN   ) ::  rpdel      (plond,plev)     ! 1./pdel (thickness between interfaces)
    REAL(kind=r8), INTENT(IN   ) ::  rpdeli     (plond,plev)     ! 1./pdeli (thickness between midpoints)
    REAL(kind=r8), INTENT(IN   ) ::  um1        (plond,plev)     ! u-wind input
    REAL(kind=r8), INTENT(IN   ) ::  vm1        (plond,plev)     ! v-wind input
    REAL(kind=r8), INTENT(IN   ) ::  tm1        (plond,plev)     ! temperature input
    REAL(kind=r8), INTENT(IN   ) ::  taux       (plond)          ! x surface stress ![N/m**2]
    REAL(kind=r8), INTENT(IN   ) ::  tauy       (plond)          ! y surface stress ![N/m**2]
    REAL(kind=r8), INTENT(IN   ) ::  shflx      (plond)          ! surface sensible heat flux (w/m2)
    REAL(kind=r8), INTENT(IN   ) ::  cflx       (plond,pcnst)    ! surface constituent flux (kg/m2/s)
    REAL(kind=r8), INTENT(IN   ) ::  ztodt                       ! 2 delta-t
    REAL(kind=r8), INTENT(IN   ) ::  colrad     (plond)    
    !
    ! Input/output arguments
    !
    REAL(kind=r8), INTENT(IN   ) :: qm1(plond,plev,pcnst)  ! initial/final constituent field
    !
    ! Output arguments
    !
    REAL(kind=r8), INTENT(OUT  ) :: dtv(plond,plev)        ! temperature tendency (heating)
    REAL(kind=r8), INTENT(OUT  ) :: dqv(plond,plev,pcnst)  ! constituent diffusion tendency

    REAL(kind=r8), INTENT(OUT  ) :: duv(plond,plev)        ! u-wind tendency
    REAL(kind=r8), INTENT(OUT  ) :: dvv(plond,plev)        ! v-wind tendency
    REAL(kind=r8), INTENT(OUT  ) :: up1(plond,plev)        ! u-wind after vertical diffusion
    REAL(kind=r8), INTENT(OUT  ) :: vp1(plond,plev)        ! v-wind after vertical diffusion
    REAL(kind=r8), INTENT(OUT  ) :: pblh(plond)            ! planetary boundary layer height
    REAL(kind=r8), INTENT(OUT  ) :: tpert(plond)           ! convective temperature excess
    REAL(kind=r8), INTENT(OUT  ) :: qpert(plond,pcnst)     ! convective humidity and constituent excess
    REAL(kind=r8), INTENT(INOUT) :: TKE(plond,plev+1)
    REAL(kind=r8), INTENT(INOUT) :: rino(plond,plev)        ! bulk Richardson no. from level to ref lev
    REAL(kind=r8), INTENT(OUT  ) :: obklen(plond)           ! Obukhov length
    REAL(kind=r8), INTENT(OUT  ) :: phiminv(plond)          ! inverse phi function for momentum
    REAL(kind=r8), INTENT(OUT  ) :: phihinv(plond)          ! inverse phi function for heat 
    REAL(kind=r8), INTENT(INOUT) :: kvh(plond,plev + 1)        ! diffusion coefficient for heat
    REAL(kind=r8), INTENT(INOUT) :: kvm(plond,plev + 1)        ! diffusion coefficient for momentum
    REAL(kind=r8), INTENT(OUT  ) ::tstar(plond)    
    REAL(kind=r8), INTENT(OUT  ) ::wstar(plond)    
    !
    !---------------------------Local storage-------------------------------
    !
    INTEGER :: i,k,m               ! longitude,level,constituent indices
    REAL(kind=r8) ::  denom                  ! denominator of expression
    REAL(kind=r8) ::  qp1   (plond,plev,pcnst)  ! constituents after vdiff
    REAL(kind=r8) ::  thp   (plond,plev)        ! potential temperature after vdiff
    REAL(kind=r8) ::  cgs   (plond,plev + 1)       ! counter-gradient star (cg/flux)
    REAL(kind=r8) ::  rztodt                 ! 1./ztodt
    REAL(kind=r8) ::  br    (plond)
    REAL(kind=r8) ::  THETA1(plond) 
    REAL(kind=r8) ::  TVS   (plond) 
    REAL(kind=r8) ::  THV1  (plond) 
    REAL(kind=r8) ::  DTV2  (plond) 
    REAL(kind=r8) ::  PS1   (plond) 
    REAL(kind=r8) ::  PS    (plond) 
    REAL(kind=r8) ::  Z1    (plond)
    REAL(kind=r8) ::  wspd  (plond)
    REAL(kind=r8) ::  TV1   (plond)
    REAL(KIND=r8),PARAMETER :: con_rv               =4.6150e+2_r8 ! gas constant H2O    (J/kg/K)
    REAL(KIND=r8),PARAMETER :: con_rd               =2.8705e+2_r8 ! gas constant air    (J/kg/K)
    REAL(KIND=r8),PARAMETER :: RVRDM1               =con_rv/con_rd-1.0_r8

    !REAL(kind=r8) :: tke(plond,plev + 1)
    REAL(KIND=r8) :: PBLH_MY20  (plond)
    REAL(KIND=r8) :: PBLH_MY25  (plond)
    REAL(KIND=r8) :: PBLH_HSBO  (plond)

    REAL(KIND=r8) :: KmMixl (plond,plev) 
    REAL(KIND=r8) :: KhMixl (plond,plev)
    REAL(KIND=r8) :: AKH2   (1:plond,1:plev+1)
    REAL(KIND=r8) :: AKM2   (1:plond,1:plev+1)
    REAL(KIND=r8) :: el2    (plond  ,plev + 1)         ! el       (1:nCols,1:kMax+1)!       
    REAL(KIND=r8) :: el0    (plond)
    REAL(KIND=r8) :: rrho   (1:plond)
    REAL(KIND=r8) :: ustr        
    REAL(KIND=r8) :: twodt

    !
    dtv  =0.0_r8! temperature tendency (heating)
    dqv  =0.0_r8! constituent diffusion tendency
    duv  =0.0_r8      ! u-wind tendency
    dvv  =0.0_r8      ! v-wind tendency
    obklen=0.0_r8           ! Obukhov length
    phiminv =0.0_r8        ! inverse phi function for momentum
    phihinv =0.0_r8        ! inverse phi function for heat 
    tstar   =0.0_r8  
    wstar   =0.0_r8  
    cgs =0.0_r8  
    br    =0.0_r8  
    THETA1=0.0_r8  
    TVS   =0.0_r8  
    THV1  =0.0_r8  
    DTV2  =0.0_r8  
    PS1   =0.0_r8  
    PS    =0.0_r8  
    Z1    =0.0_r8  
    wspd  =0.0_r8  
    TV1   =0.0_r8  
    KmMixl=0.0_r8  
    KhMixl=0.0_r8  
    AKH2  =0.0_r8  
    AKM2  =0.0_r8  
    el2   =0.0_r8  
    el0   =0.0_r8  
    rrho  =0.0_r8  
    ustr  =0.0_r8  
    twodt=0.0_r8  
    up1    =um1
    vp1    =vm1
    thp    =tm1
    qp1    =qm1
    pblh   =0.0_r8
    ustar  =0.0_r8
    kvh    =0.0_r8
    kvm    =0.0_r8
    tpert  =0.0_r8
    qpert  =0.0_r8
    cgs    =0.0_r8
    KmMixl =0.0_r8
    KhMixl =0.0_r8
    AKH2   =0.0_r8
    AKM2   =0.0_r8
    el2    =0.0_r8
    el0    =0.0_r8
    twodt =ztodt
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------    
    CALL MY20_TURB (       &
                   pcnst                                ,& ! INTENT(IN   ) :: pcnst
                   plond                                ,& ! INTENT(IN   ) :: plond
                   plev                                 ,& ! INTENT(IN   ) :: plev
                   twodt                                ,& ! INTENT(IN   ) :: ztodt
                   colrad     (1:plond)                 ,& ! INTENT(IN   ) :: colrad    (plond)
                   pintm1     (1:plond,1:plev+1) ,& ! phalf    (1:plond,1:plev+1)!!          phalf    -  Pressure at half levels, 
                   pmidm1     (1:plond,1:plev)   ,& ! pfull    (1:plond,1:plev)!!          pfull    -  Pressure at full levels! , zhalf, zfull
                   tm1        (1:plond,1:plev)          ,& ! INTENT(IN   ) :: tm1       (plond,plev)! temperature input
                   qp1        (1:plond,1:plev,1:pcnst)  ,& ! INTENT(IN   ) :: qp1       (plond,plev,pcnst)  ! constituents after vdiff
                   um1        (1:plond,1:plev)          ,& ! INTENT(IN   ) :: um1       (plond,plev)! u-wind input
                   vm1        (1:plond,1:plev)          ,& ! INTENT(IN   ) :: vm1       (plond,plev)! v-wind input
                   gl0        (1:plond)                 ,& ! INTENT(INOUT) :: gl0       (plond)!
                   PBLH_MY20  (1:plond)                 ,& ! INTENT(OUT  ) :: PBLH_MY20(plond,plev)   
                   TKEMYJ_MY20(1:plond,1:plev+1,latco)  ,& ! INTENT(OUT  ) :: TKEMYJ_MY20(plond,plev)   
                   KmMixl     (1:plond,1:plev)          ,& ! INTENT(OUT  ) :: Pbl_KmMixl(plond,plev)   
                   KhMixl     (1:plond,1:plev)           ) ! INTENT(OUT  ) :: Pbl_KhMixl(plond,plev)

    
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
    DO I=1,plond
       rrho(i)  = rair*tm1(i,plev)/pmidm1(i,plev)
       ustr     = SQRT(SQRT(taux(i)**2 + tauy(i)**2)*rrho(i))
       ustar(i) = MAX(ustr,0.01_r8) ! [m/s]
    END DO
    CALL MY25_TURB(&
       plond                         ,& ! plond
       plev                          ,& ! plev
       twodt                         ,& ! delt    !       delt     -  Time step in seconds
       um1        (1:plond,1:plev)   ,& ! um       (1:plond,1:plev)!Zonal wind components
       vm1        (1:plond,1:plev)   ,& ! vm       (1:plond,1:plev)!Meridional Wind components
       thm        (1:plond,1:plev)   ,& ! theta    (1:plond,1:plev)!theta     -  Potential temperature
       zhalf      (1:plond,1:plev+1) ,& ! zhalf    (1:plond,1:plev+1)!zfull         -  Height at full levels
       zm         (1:plond,1:plev)   ,& ! zfull    (1:plond,1:plev)!         zhalf    -  Height at half levels
       pintm1     (1:plond,1:plev+1) ,& ! phalf    (1:plond,1:plev+1)!!          phalf    -  Pressure at half levels, 
       pmidm1     (1:plond,1:plev)   ,& ! pfull    (1:plond,1:plev)!!          pfull    -  Pressure at full levels! , zhalf, zfull
       z0         (1:plond)          ,& ! z0       (1:plond)!   z0-  Roughness length
       FRLAND     (1:plond)          ,& ! fracland (1:plond)! fracland -  Fractional amount of land beneath a grid box
       ustar      (1:plond)          ,& ! ustar    (1:plond)! ustar         -  OPTIONAL:friction velocity (m/sec)
       el0        (1:plond)          ,& ! el0      (1:plond)!el0  -  characteristic length scale
       gl0        (1:plond)          ,& ! gl0      (plond)!
       PBLH_MY25  (1:plond)          ,& ! pblh     (1:plond)!el0  -  characteristic length scale
       akm2       (1:plond,1:plev+1) ,& ! akm      (1:plond,1:plev+1)!       akm  -  mixing coefficient for momentum
       akh2       (1:plond,1:plev+1) ,& ! akh      (1:plond,1:plev+1)!       akh  -  mixing coefficient for heat and moisture
       el2        (1:plond,1:plev+1) ,& ! el       (1:plond,1:plev+1)!       el   -  master length scale
       TKEMYJ_MY25(1:plond,1:plev+1,latco)  ) !TKE       (1:nCols,1:plev+1)!
!    !----------------------------------------------------------------------
!    !**********************************************************************
!    !----------------------------------------------------------------------

    DO i=1,plond
       wspd  (i) = SQRT((um1(i,plev)**2) + (vm1(i,plev)**2))
       wspd  (i) = MAX(wspd  (i),0.25_r8)
       THETA1(i) = tsk   (i) * psomc(i,plev)       
       TVS   (i) = tsk   (i) * (1.0_r8 + RVRDM1 * QSFC(i))
       THV1  (i) = THETA1(i) * (1.0_r8 + RVRDM1 * QSFC(i))
       TV1   (i) = tm1 (i,plev) * (1.0_r8 + RVRDM1 * QSFC(i))
       PS1   (i) = pmidm1(i,plev)
       PS    (i) = psfcpa(i)
       DTV2   (i) = THV1(i) - TVS(i)
       Z1    (i) = -con_rd * TV1(i) * LOG(PS1(i)/PS(i)) / grav
       br    (i) = grav * DTV2(i) * Z1(i) / (0.5_r8 * (THV1(i) + TVS(i))&
     &          * wspd(i) * wspd(i))
       br(i) = MAX(br(i),-5000.0_r8)
    END DO                  
    !
    !-----------------------------------------------------------------------
    !
    ! Call vertical diffusion code. No setup work is required.
    !
    CALL vdiff(bstar       (1:plond)           , &! INTENT(IN        ) :: 
               FRLAND      (1:plond)           , &! INTENT(IN        ) :: 
               LwCoolRate  (1:plond,1:plev)    , &! INTENT(IN        ) :: 
               LwCoolRateC (1:plond,1:plev)    , &! INTENT(IN        ) :: 
               cldtot      (1:plond,1:plev)    , &! INTENT(IN        ) :: 
               qliq        (1:plond,1:plev)    , &! INTENT(IN        ) :: 
               um1         (1:plond,1:plev)    , &! INTENT(IN        ) :: um1(plond,plev)            ! u wind input
               vm1         (1:plond,1:plev)    , &! INTENT(IN        ) :: vm1(plond,plev)            ! v wind input
               tm1         (1:plond,1:plev)    , &! INTENT(IN        ) :: tm1(plond,plev)            ! temperature input
               qm1         (1:plond,1:plev,1:pcnst), &! INTENT(IN        ) :: qm1(plond,plev,pcnst)  ! moisture and trace constituent input
               pmidm1      (1:plond,1:plev)    , &! INTENT(IN        ) :: pmidm1(plond,plev)     ! midpoint pressures
               pintm1      (1:plond,1:plev+1)  , &! INTENT(IN        ) :: pintm1(plond,plev + 1)    ! interface pressures
               rpdel       (1:plond,1:plev)    , &! INTENT(IN        ) :: rpdel(plond,plev)      ! 1./pdel  (thickness bet interfaces)
               rpdeli      (1:plond,1:plev)    , &! INTENT(IN        ) :: rpdeli(plond,plev)     ! 1./pdeli (thickness bet midpoints)
               ztodt                           , &! INTENT(IN        ) :: ztodt                    ! 2 delta-t
               thm         (1:plond,1:plev)    , &! INTENT(IN        ) :: thm(plond,plev)            ! potential temperature
               zm          (1:plond,1:plev)    , &! INTENT(IN        ) :: zm(plond,plev)            ! midpoint geoptl height above sfc
               taux        (1:plond)           , &! INTENT(IN        ) :: taux(plond)            ! x surface stress (n)
               tauy        (1:plond)           , &! INTENT(IN        ) :: tauy(plond)            ! y surface stress (n)
               shflx       (1:plond)           , &! INTENT(IN        ) :: shflx(plond)            ! surface sensible heat flux (w/m2)
               cflx        (1:plond,1:pcnst)   , &! INTENT(IN        ) :: cflx(plond,pcnst)      ! surface constituent flux (kg/m2/s)
               up1         (1:plond,1:plev)    , &! INTENT(OUT  ) :: up1(plond,plev)            ! u-wind after vertical diffusion
               vp1         (1:plond,1:plev)    , &! INTENT(OUT  ) :: vp1(plond,plev)            ! v-wind after vertical diffusion
               thp         (1:plond,1:plev)    , &! INTENT(OUT  ) :: thp(plond,plev)            ! pot temp after vert. diffusion
               qp1         (1:plond,1:plev,1:pcnst), &! INTENT(OUT  ) :: qp1(plond,plev,pcnst)  ! moist, tracers after vert. diff
               pblh        (1:plond)           , &! INTENT(OUT  ) :: pblh(plond)            ! planetary boundary layer height
               ustar       (1:plond)           , &! INTENT(OUT  ) :: ustar(plond)            ! surface friction velocity
               kvh         (1:plond,1:plev + 1), &! INTENT(OUT  ) :: kvh(plond,plev + 1)       ! coefficient for heat and tracers
               kvm         (1:plond,1:plev + 1), &! INTENT(OUT  ) :: kvm(plond,plev + 1)       ! coefficient for momentum
               tke         (1:plond,1:plev + 1), &! INTENT(OUT  ) :: kvm(plond,plev + 1)       ! coefficient for momentum
               tpert       (1:plond)           , &! INTENT(OUT  ) :: tpert(plond)            ! convective temperature excess
               qpert       (1:plond,1:pcnst)   , &! INTENT(OUT  ) :: qpert(plond)            ! convective humidity excess
               cgs         (1:plond,1:plev+1)  , &! INTENT(OUT  ) :: cgs(plond,plev + 1)       ! counter-grad star (cg/flux)
               plon                            , &! INTENT(IN        ) :: plon        ! number of longitudes
               plev                            , &! INTENT(IN        ) :: plev        ! number of vertical levels
               plond                           , &! INTENT(IN        ) :: plond        ! slt extended domain longitude
               pcnst                           , &! INTENT(IN        ) :: pcnst        ! number of constituents (including water vapor)
               fakn                            , &! INTENT(IN        ) :: fakn    ! Constant in turbulent prandtl number
               ricr                            , &! INTENT(IN        ) :: ricr    ! Critical richardson number
               tsk         (1:plond)           , &! (plond)
               qsfc        (1:plond)           , &
               KmMixl      (1:plond,1:plev)    , &
               KhMixl      (1:plond,1:plev)    , &
               AKH2        (1:plond,1:plev+1)  , &! akh      (1:nCols,1:plev+1)!  mixing coefficient for heat and moisture 
               AKM2        (1:plond,1:plev+1)  , &! akm      (1:nCols,1:plev+1)!  mixing coefficient for momentum
               rino        (1:plond,1:plev)    , &! REAL(KIND=r8), INTENT(INOUT) :: rino    (plond,plev)
               psomc       (1:plond,1:plev)    , &!REAL(KIND=r8), INTENT(IN   ) ::  psomc(plond,plev)      ! (psm1/pmidm1)**cappa
               obklen      (1:plond)           , &! REAL(KIND=r8), INTENT(OUT  ) :: obklen  (plond)           ! Obukhov length
               phiminv     (1:plond)           , &!REAL(KIND=r8), INTENT(OUT  ) :: phiminv (plond)          ! inverse phi function for momentum
               phihinv     (1:plond)           , &!REAL(KIND=r8), INTENT(OUT  ) :: phihinv (plond)          ! inverse phi function for heat 
               tstar       (1:plond)           , & !REAL(KIND=r8), INTENT(OUT  ) :: tstar(plond)
               wstar       (1:plond)           , &   !REAL(KIND=r8), INTENT(OUT  ) :: wstar(plond)
               PBLH_MY20   (1:plond)           , &
               PBLH_MY25   (1:plond)           , &
               PBLH_HSBO   (1:plond)           , &
               TKEMYJ_MY20 (1:plond,1:plev+1,latco),&
               TKEMYJ_MY25 (1:plond,1:plev+1,latco),&
               TKEMYJ_HSBO (1:plond,1:plev+1,latco))
    !----------------------------------------------------------------------
    !**********************************************************************
    !----------------------------------------------------------------------
 
    !
    ! Convert the diffused fields back to diffusion tendencies.
    ! Add the diffusion tendencies to the cummulative physics tendencies,
    ! except for constituents. The diffused values of the constituents
    ! replace the input values.
    !
    rztodt = 1.0_r8/ztodt
    DO k=ntopfl,plev
       DO i=1,plon
          duv(i,k) = (up1(i,k)*SIN( colrad(i)) - um1(i,k)*SIN( colrad(i)))*rztodt
          dvv(i,k) = (vp1(i,k)*SIN( colrad(i)) - vm1(i,k)*SIN( colrad(i)))*rztodt

          denom    = cpair*(1.0_r8 + cpvir*qm1(i,k,1))
          dtv(i,k) = (thp(i,k) - thm(i,k))*rztodt

       END DO
       DO m=1,pcnst
          DO i=1,plon
             dqv(i,k,m) = (qp1(i,k,m) - qm1(i,k,m))*rztodt
          END DO
       END DO
    END DO    
    DO i=1,plon
       IF(obklen  (i) <1.0e-12_r8 .and. obklen  (i) >-1.0e-12_r8)obklen  (i)=0.0_r8
       obklen  (i)=MIN(MAX(obklen (i),-1.e12_r8),1.e12_r8) 
       IF(phiminv  (i) <1.0e-12_r8 .and. phiminv  (i) >-1.0e-12_r8)phiminv  (i)=0.0_r8
       phiminv (i)=MIN(MAX(phiminv (i),-1.e12_r8),1.e12_r8) 
       IF(phihinv  (i) <1.0e-12_r8 .and. phihinv  (i) >-1.0e-12_r8)phihinv  (i)=0.0_r8
       phihinv (i)=MIN(MAX(phihinv (i),-1.e12_r8),1.e12_r8) 
    END DO    

    RETURN
  END SUBROUTINE vdintr

  !
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !---------------------------------------------------------------------

  SUBROUTINE vdiff(bstar , &! INTENT(IN   ) :: 
                   FRLAND, &! INTENT(IN   ) :: 
                   LwCoolRate  , &! INTENT(IN   ) :: 
                   LwCoolRateC , &! INTENT(IN   ) :: 
                   cldtot      , &! INTENT(IN   ) :: 
                   qliq    , &! INTENT(IN   ) :: 
                   um1     , &! INTENT(IN   ) :: um1(plond,plev)        ! u wind input
                   vm1     , &! INTENT(IN   ) :: vm1(plond,plev)        ! v wind input
                   tm1     , &! INTENT(IN   ) :: tm1(plond,plev)        ! temperature input
                   qm1     , &! INTENT(IN   ) :: qm1(plond,plev,pcnst)  ! moisture and trace constituent input
                   pmidm1  , &! INTENT(IN   ) :: pmidm1(plond,plev)     ! midpoint pressures
                   pintm1  , &! INTENT(IN   ) :: pintm1(plond,plev + 1)    ! interface pressures
                   rpdel   , &! INTENT(IN   ) :: rpdel(plond,plev)      ! 1./pdel  (thickness bet interfaces)
                   rpdeli  , &! INTENT(IN   ) :: rpdeli(plond,plev)     ! 1./pdeli (thickness bet midpoints)
                   ztodt   , &! INTENT(IN   ) :: ztodt                  ! 2 delta-t
                   thm     , &! INTENT(IN   ) :: thm(plond,plev)        ! potential temperature
                   zm      , &! INTENT(IN   ) :: zm(plond,plev)         ! midpoint geoptl height above sfc
                   taux           , &! INTENT(IN   ) :: taux(plond)            ! x surface stress (n)
                   tauy           , &! INTENT(IN   ) :: tauy(plond)            ! y surface stress (n)
                   shflx   , &! INTENT(IN   ) :: shflx(plond)           ! surface sensible heat flux (w/m2)
                   cflx    , &! INTENT(IN   ) :: cflx(plond,pcnst)      ! surface constituent flux (kg/m2/s)
                   up1     , &! INTENT(OUT  ) :: up1(plond,plev)        ! u-wind after vertical diffusion
                   vp1     , &! INTENT(OUT  ) :: vp1(plond,plev)        ! v-wind after vertical diffusion
                   thp     , &! INTENT(OUT  ) :: thp(plond,plev)        ! pot temp after vert. diffusion
                   qp1     , &! INTENT(OUT  ) :: qp1(plond,plev,pcnst)  ! moist, tracers after vert. diff
                   pblh    , &! INTENT(OUT  ) :: pblh(plond)            ! planetary boundary layer height
                   ustar   , &! INTENT(OUT  ) :: ustar(plond)           ! surface friction velocity
                   kvh     , &! INTENT(OUT  ) :: kvh(plond,plev + 1)       ! coefficient for heat and tracers
                   kvm           , &! INTENT(OUT  ) :: kvm(plond,plev + 1)       ! coefficient for momentum
                   tke           , &! INTENT(OUT  ) :: tke(plond,plev + 1)       ! coefficient for momentum
                   tpert   , &! INTENT(OUT  ) :: tpert(plond)           ! convective temperature excess
                   qpert   , &! INTENT(OUT  ) :: qpert(plond)           ! convective humidity excess
                   cgs     , &! INTENT(OUT  ) :: cgs(plond,plev + 1)       ! counter-grad star (cg/flux)
                   plon    , &! INTENT(IN   ) :: plon       ! number of longitudes
                   plev    , &! INTENT(IN   ) :: plev       ! number of vertical levels
                   plond   , &! INTENT(IN   ) :: plond      ! slt extended domain longitude
                   pcnst   , &! INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)
                   fakn    , &! INTENT(IN   ) :: fakn    ! Constant in turbulent prandtl number
                   ricr    , &! INTENT(IN   ) :: ricr    ! Critical richardson number
                   tsk     , &
                   qsfc    , &
                   KmMixl  , &
                   KhMixl  , &
                   AKH     , &
                   AKM     , &
                   rino    , &
                   psomc   , &
                   obklen  , &
                   phiminv , &
                   phihinv , &
                   tstar   , &
                   wstar   , &
                   PBLH_MY20 , &
                   PBLH_MY25 , &
                   PBLH_HSBO , &
                   TKEMYJ_MY20,&
                   TKEMYJ_MY25,&
                   TKEMYJ_HSBO)
    !-----------------------------------------------------------------------
    !
    ! Driver routine to compute vertical diffusion of momentum,
    ! moisture, trace constituents and potential temperature.
    !
    ! Free atmosphere diffusivities are computed first; then modified
    ! by the boundary layer scheme; then passed to individual
    ! parameterizations mvdiff, qvdiff
    !
    ! The free atmosphere diffusivities are based on standard mixing length 
    ! forms for the neutral diffusivity multiplied by functns of Richardson 
    ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for 
    ! momentum, potential temperature, and constitutents.
    ! The stable Richardson num function (Ri>0) is taken from Holtslag and 
    ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))
    ! The unstable Richardson number function (Ri<0) is taken from  CCM1.
    ! f = sqrt(1 - 18*Ri)
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  CCM1
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          P. Rasch, B. Boville, August 1992
    ! Reviewed:          P. Rasch, March 1996
    ! Reviewed:          B. Boville, April 1996
    !
    !-----------------------------------------------------------------------
    !
    ! $Id: vdiff.F,v 1.1.1.1 2001/03/09 00:29:28 mirin Exp $
    !
    !-----------------------------------------------------------------------
    !      use precision
    !      use pmgrid
    !#if ( ! defined CRAY )
    !      use srchutil, only: wheneq
    !#endif
    !-----------------------------------------------------------------------
    !#include <implicit.h>
    !------------------------------Commons----------------------------------
    !#include <comvd.h>
    !-----------------------------------------------------------------------
    !#include <comqmin.h>
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) ::  plon       ! number of longitudes
    INTEGER, INTENT(IN   ) :: plev       ! number of vertical levels
    INTEGER, INTENT(IN   ) :: plond      ! slt extended domain longitude
    INTEGER, INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)
    REAL(KIND=r8), INTENT(IN   ) :: fakn    ! Constant in turbulent prandtl number
    REAL(KIND=r8), INTENT(IN   ) :: ricr    ! Critical richardson number

    REAL(KIND=r8), INTENT(IN   ) :: bstar      (plond)      
    REAL(KIND=r8), INTENT(IN   ) :: FRLAND  (plond)      
    REAL(KIND=r8), INTENT(IN   ) :: LwCoolRate (plond,plev)        
    REAL(KIND=r8), INTENT(IN   ) :: LwCoolRateC(plond,plev)        
    REAL(KIND=r8), INTENT(IN   ) :: cldtot     (plond,plev)      
    REAL(KIND=r8), INTENT(IN   ) :: qliq       (plond,plev)      
    REAL(KIND=r8), INTENT(IN   ) :: um1(plond,plev)        ! u wind input
    REAL(KIND=r8), INTENT(IN   ) :: vm1(plond,plev)        ! v wind input
    REAL(KIND=r8), INTENT(IN   ) :: tm1(plond,plev)        ! temperature input
    REAL(KIND=r8), INTENT(IN   ) :: qm1(plond,plev,pcnst)  ! moisture and trace constituent input
    REAL(KIND=r8), INTENT(IN   ) :: pmidm1(plond,plev)     ! midpoint pressures
    REAL(KIND=r8), INTENT(IN   ) :: pintm1(plond,plev + 1)    ! interface pressures
    REAL(KIND=r8), INTENT(IN   ) :: rpdel(plond,plev)      ! 1./pdel  (thickness bet interfaces)
    REAL(KIND=r8), INTENT(IN   ) :: rpdeli(plond,plev)     ! 1./pdeli (thickness bet midpoints)
    REAL(KIND=r8), INTENT(IN   ) :: ztodt                  ! 2 delta-t
    REAL(KIND=r8), INTENT(IN   ) :: thm(plond,plev)        ! potential temperature
    REAL(KIND=r8), INTENT(IN   ) :: zm(plond,plev)         ! midpoint geoptl height above sfc
    REAL(KIND=r8), INTENT(IN   ) :: taux(plond)            ! x surface stress (n)
    REAL(KIND=r8), INTENT(IN   ) :: tauy(plond)            ! y surface stress (n)
    REAL(KIND=r8), INTENT(IN   ) :: shflx(plond)           ! surface sensible heat flux (w/m2)
    REAL(KIND=r8), INTENT(IN   ) :: cflx(plond,pcnst)      ! surface constituent flux (kg/m2/s)
    REAL(KIND=r8), INTENT(IN   ) :: tsk   (plond)     
    REAL(KIND=r8), INTENT(IN   ) :: qsfc   (plond)    
    REAL(KIND=r8), INTENT(IN   ) :: KmMixl(plond,plev) 
    REAL(KIND=r8), INTENT(IN   ) :: KhMixl(plond,plev) 
    REAL(KIND=r8), INTENT(IN   ) :: AKH(plond,plev+1) 
    REAL(KIND=r8), INTENT(IN   ) :: AKM(plond,plev+1) 
    REAL(KIND=r8), INTENT(IN   ) ::  psomc(plond,plev)      ! (psm1/pmidm1)**cappa
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(OUT  ) :: up1     (plond,plev)        ! u-wind after vertical diffusion
    REAL(KIND=r8), INTENT(OUT  ) :: vp1     (plond,plev)        ! v-wind after vertical diffusion
    REAL(KIND=r8), INTENT(OUT  ) :: thp     (plond,plev)        ! pot temp after vert. diffusion
    REAL(KIND=r8), INTENT(OUT  ) :: qp1     (plond,plev,pcnst)  ! moist, tracers after vert. diff
    REAL(KIND=r8), INTENT(INOUT) :: pblh    (plond)            ! planetary boundary layer height
    REAL(KIND=r8), INTENT(INOUT) :: rino    (plond,plev)        ! bulk Richardson no. from level to ref lev
    REAL(KIND=r8), INTENT(OUT  ) :: ustar   (plond)           ! surface friction velocity
    REAL(KIND=r8), INTENT(OUT  ) :: kvh     (plond,plev + 1)       ! coefficient for heat and tracers
    REAL(KIND=r8), INTENT(OUT  ) :: kvm     (plond,plev + 1)       ! coefficient for momentum
    REAL(KIND=r8), INTENT(OUT  ) :: tke     (plond,plev + 1)       ! coefficient for momentum
    REAL(KIND=r8), INTENT(OUT  ) :: tpert   (plond)           ! convective temperature excess
    REAL(KIND=r8), INTENT(OUT  ) :: qpert   (plond)           ! convective humidity excess
    REAL(KIND=r8), INTENT(OUT  ) :: cgs     (plond,plev + 1)       ! counter-grad star (cg/flux)
    REAL(KIND=r8), INTENT(OUT  ) :: obklen  (plond)           ! Obukhov length
    REAL(KIND=r8), INTENT(OUT  ) :: phiminv (plond)          ! inverse phi function for momentum
    REAL(KIND=r8), INTENT(OUT  ) :: phihinv (plond)          ! inverse phi function for heat 
    REAL(KIND=r8), INTENT(OUT  ) :: tstar(plond)
    REAL(KIND=r8), INTENT(OUT  ) :: wstar(plond)
    REAL(KIND=r8), INTENT(INOUT) :: TKEMYJ_MY20(plond,plev + 1)   
    REAL(KIND=r8), INTENT(INOUT) :: TKEMYJ_MY25(plond,plev + 1)   
    REAL(KIND=r8), INTENT(INOUT) :: TKEMYJ_HSBO(plond,plev + 1)   
    REAL(KIND=r8), INTENT(INOUT) :: PBLH_MY20(plond)
    REAL(KIND=r8), INTENT(INOUT) :: PBLH_MY25(plond)
    REAL(KIND=r8), INTENT(INOUT) :: PBLH_HSBO(plond)

    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(KIND=r8) cah   (plond,plev)        ! -upper diag for heat and constituts
    REAL(KIND=r8) cam   (plond,plev)        ! -upper diagonal for momentum
    REAL(KIND=r8) cch   (plond,plev)        ! -lower diag for heat and constits
    REAL(KIND=r8) ccm   (plond,plev)        ! -lower diagonal for momentum
    REAL(KIND=r8) cgh   (plond,plev + 1)       ! countergradient term for heat
    REAL(KIND=r8) cgq   (plond,plev + 1,pcnst) ! countergrad term for constituent
    REAL(KIND=r8) dvdz2                  ! (du/dz)**2 + (dv/dz)**2
    REAL(KIND=r8) dz                     ! delta z between midpoints
    REAL(KIND=r8) fstab                  ! stable f(ri)
    REAL(KIND=r8) funst                  ! unstable f(ri)
    REAL(KIND=r8) kvf   (plond,plev + 1)       ! free atmosphere kv at interfaces
    REAL(KIND=r8) rinub                  ! richardson no=(g/theta)(dtheta/dz)/
    !                                                 (du/dz**2+dv/dz**2)
    REAL(KIND=r8) sstab                  ! static stability = g/th  * dth/dz
    REAL(KIND=r8) potbar(plond,plev + 1)    ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
    REAL(KIND=r8) tmp1  (plond)            ! temporary storage
    REAL(KIND=r8) tmp2                   ! temporary storage
    REAL(KIND=r8) rcpair                 ! 1./cpair
    REAL(KIND=r8) ztodtgor               ! ztodt*gravit/rair
    REAL(KIND=r8) gorsq                  ! (gravit/rair)**2
    REAL(KIND=r8) dubot (plond)           ! lowest layer u change from stress
    REAL(KIND=r8) dvbot (plond)           ! lowest layer v change from stress
    REAL(KIND=r8) dtbot (plond)           ! lowest layer t change from heat flx
    REAL(KIND=r8) dqbot (plond,pcnst)     ! lowest layer q change from const flx
    REAL(KIND=r8) thx   (plond,plev)        ! temperature input + counter gradient
    REAL(KIND=r8) thv   (plond,plev)        ! virtual potential temperature
    REAL(KIND=r8) qmx   (plond,plev,pcnst)  ! constituents input + counter grad
    REAL(KIND=r8) zeh   (plond,plev)        ! term in tri-diag. matrix system (T & Q)
    REAL(KIND=r8) zem   (plond,plev)        ! term in tri-diag. matrix system (momentum)
    REAL(KIND=r8) termh (plond,plev)      ! 1./(1. + cah(k) + cch(k) - cch(k)*zeh(k-1))
    REAL(KIND=r8) termm (plond,plev)      ! 1./(1. + cam(k) + ccm(k) - ccm(k)*zem(k-1))
    REAL(KIND=r8) kvn                    ! neutral Kv
    REAL(KIND=r8) ksx   (plond)             ! effective surface drag factor (x)
    REAL(KIND=r8) ksy   (plond)             ! effective surface drag factor (y)
    REAL(KIND=r8) sufac (plond)           ! lowest layer u implicit stress factor
    REAL(KIND=r8) svfac (plond)           ! lowest layer v implicit stress factor
    INTEGER indx(plond)          ! array of indices of potential q<0
    INTEGER ilogic(plond)        ! 1 => adjust vertical profile
    INTEGER nval                ! num of values which meet criteria
    INTEGER ii                  ! longitude index of found points
    INTEGER i                   ! longitude index
    INTEGER k                   ! vertical index
    INTEGER m                   ! constituent index
    INTEGER ktopbl(plond)       ! index of first midpoint inside pbl
    INTEGER ktopblmn(plond)            ! min value of ktopbl
    up1=0.0_r8;vp1  =0.0_r8; thp  =0.0_r8;  qp1=0.0_r8;  pblh =0.0_r8;  rino =0.0_r8;    ustar=0.0_r8
    kvh  =0.0_r8;    kvm  =0.0_r8;    tpert=0.0_r8;    qpert=0.0_r8;    cgs  =0.0_r8;    obklen=0.0_r8;   phiminv=0.0_r8;  phihinv=0.0_r8
    cah =0.0_r8;   cam =0.0_r8;   cch =0.0_r8;   ccm =0.0_r8;   cgh =0.0_r8;   cgq =0.0_r8;   dvdz2=0.0_r8;  dz   =0.0_r8
    fstab=0.0_r8;  funst=0.0_r8;  kvf  =0.0_r8;  rinub=0.0_r8;  sstab =0.0_r8; potbar=0.0_r8; tmp1  =0.0_r8; tmp2 =0.0_r8
    rcpair  =0.0_r8;   ztodtgor =0.0_r8;  gorsq =0.0_r8;  dubot =0.0_r8;  dvbot =0.0_r8;  dtbot =0.0_r8
    dqbot =0.0_r8;  thx   =0.0_r8;  thv   =0.0_r8;  qmx   =0.0_r8;  zeh   =0.0_r8;  zem   =0.0_r8;  termh =0.0_r8;  termm =0.0_r8
    kvn   =0.0_r8;  ksx   =0.0_r8;  ksy   =0.0_r8;  sufac =0.0_r8;  svfac =0.0_r8;nval=0;ktopbl=0;ktopblmn=0
    !
    !-----------------------------------------------------------------------
    !
    ! Convert the surface fluxes to lowest level tendencies.
    ! Stresses are converted to effective drag coefficients if these are >0
    !
    rcpair = 1.0_r8/cpair
    DO i=1,plond
       tmp1(i) = ztodt*gravit*rpdel(i,plev)

       ksx(i) = -taux(i) / um1(i,plev)
       IF (ksx(i) .GT. 0.0_r8) THEN
          sufac(i) = tmp1(i) * ksx(i)
          dubot(i) = 0.0_r8
       ELSE
          ksx(i)   = 0.0_r8
          sufac(i) = 0.0_r8
          dubot(i) = tmp1(i) * taux(i)
       END IF

       ksy(i) = -tauy(i) / vm1(i,plev)
       IF (ksy(i) .GT. 0.0_r8) THEN
          svfac(i) = tmp1(i) * ksy(i)
          dvbot(i) = 0.0_r8
       ELSE
          ksy(i)   = 0.0_r8
          svfac(i) = 0.0_r8
          dvbot(i) = tmp1(i) * tauy(i)
       END IF

       dqbot(i,1)   = cflx(i,1)*tmp1(i)
       dtbot(i)     = shflx(i)*tmp1(i)*rcpair
       kvf(i,plev + 1) = 0.0_r8
    END DO
    DO m=2,pcnst
       DO i=1,plond
          dqbot(i,m) = cflx(i,m)*tmp1(i)
       END DO
    END DO
    !
    ! Set the vertical diffusion coefficient above the top diffusion level
    !
    DO k=1,ntopfl
       DO i=1,plond
          kvf(i,k) = 0.0_r8
       END DO
    END DO
    !
    ! Compute virtual potential temperature for use in static stability 
    ! calculation.  0.61 is 1. - R(water vapor)/R(dry air).  Use 0.61 instead
    ! of a computed variable in order to obtain an identical simulation to
    ! Case 414.
    !
    CALL virtem(thm     (1:plond,1:plev), &! INTENT(IN   )! REAL(KIND=r8), INTENT(IN   ) ::  t(plond,plev)! temperature
                qm1     (1:plond,1:plev,1), &! INTENT(IN   )! REAL(KIND=r8), INTENT(IN   ) ::  q(plond,plev)! specific humidity
                0.61_r8                 , &! INTENT(IN   )! REAL(KIND=r8), INTENT(IN   ) ::  zvir! virtual temperature constant
                thv     (1:plond,1:plev), &! INTENT(OUT  )! REAL(KIND=r8), INTENT(OUT  ) ::  tv(plond,plev)! virtual temperature
                plond                   , &! INTENT(IN   )! INTEGER, INTENT(IN   ) ::  plond! slt extended domain longitude
                plev                    , &! INTENT(IN   )! INTEGER, INTENT(IN   ) ::  plev! number of vertical levels
                plon                      )! INTENT(IN   )! INTEGER, INTENT(IN   ) ::  plon! number of longitudes
    !
    ! Compute the free atmosphere vertical diffusion coefficients
    ! kvh = kvq = kvm. 
    !
    DO k=ntopfl,plev-1
       DO i=1,plon
          !
          ! Vertical shear squared, min value of (delta v)**2 prevents zero shear.
          !
          dvdz2 = (um1(i,k)-um1(i,k+1))**2 + (vm1(i,k)-vm1(i,k+1))**2
          !
          dvdz2 = MAX(dvdz2,1.e-36_r8)
          dz    = zm(i,k) - zm(i,k+1)
          !
          !            _             _ 2        _              _ 2
          !           |               |        |               |
          !           |U(z) - U(z+1)) |     +  |V(z) - V(z+1)) |  
          !           |_             _|        |_             _|  
          !  dV^2      
          ! ------ = -------------------------------------------
          !  dz^2                _             _ 2   
          !                     |               |
          !                     |Z(z) - Z(z+1)) |  
          !                     |_             _|  
          !
          !
          dvdz2 = dvdz2/(dz**2)

          !
          ! Static stability (use virtual potential temperature)
          !
          !                _                                              _
          !               |                                                |
          !               |         THETA(z) - THETA(z+1)                  |
          !  sstab = 2*g* |----------------------------------------------- |
          !               |  _                    _     _              _   |
          !               | |                      |   |                |  |
          !               | |THETA(z) + THETA(z+1) | * |Z(z) - Z(z+1))  |  |
          !               | |_                    _|   |_              _|  |
          !               |_                                              _|
          !
          sstab = gravit*2.0_r8*( thv(i,k) - thv(i,k+1))/ &
                                ((thv(i,k) + thv(i,k+1))*dz)
          !
          ! Richardson number, stable and unstable modifying functions
          !
          !                 _                                                                          _
          !                |                                         _               _                  |    
          !                |                                        |                 |                 |   
          !                |             THETA(z) - THETA(z+1)  *   | Z(z) - Z(z+1))  |                 |   
          !                !                                        |_               _|                 !   
          !rinub =    2*g* |--------------------------------------------------------------------------- |  
          !                |                             _                                           _  |
          !                |  _                    _    |  _              _ 2      _             _ 2  | |
          !                | |                      |   | |                |      |               |   | |
          !                | |THETA(z) + THETA(z+1) | * | | U(z) - U(z+1)) |   +  |V(z) - V(z+1)) |   | |
          !                | |_                    _|   | |_              _|      |_             _|   | |
          !                |                            |_                                           _| |
          !                |_                                                                          _|
          !
          rinub = sstab/dvdz2

          !
          !
          !                              1 
          ! fstab = ------------------------------------------------
          !             _                                 _
          !            |                 _           _    | 
          !            |                |             |   | 
          !            | 1 + 10*rinub * | 1 + 8*rinub |   | 
          !            |                |_           _|   | 
          !            |_                                _|
          !
          fstab = 1.0_r8/(1.0_r8 + 10.0_r8*rinub*(1.0_r8 + 8.0_r8*rinub))
          !
          !
          !
          funst = MAX(1.0_r8 - 18.0_r8*rinub,0.0_r8)
          !
          ! Select the appropriate function of the richardson number
          !
          !         -                                              - 0.5
          !        |                                                |
          !        |                                                |
          !        |                         1                      |
          ! fstab =|------------------------------------------------|
          !        |    _                                _          |
          !        |   |                 _           _    |         | 
          !        |   |                |             |   |         | 
          !        |   | 1 + 10*rinub * | 1 + 8*rinub |   |         | 
          !        |   |                |_           _|   |         | 
          !        |   |_                                _|         | 
          !        |_                                              _|
          !
          IF (rinub < 0.0_r8) fstab = SQRT(funst)
          !
          !            Km 
          !Km = l^2 * -----
          !            l^2
          !
          !
          ! Neutral diffusion coefficient
          ! compute mixing length (z), where z is the interface height estimated
          ! with an 8 km scale height.
          !
          !      0 > l < ml2
          !
          !
          ! Set the square of the mixing lengths.
          !
          !IF(k==1)THEN
          !   MixLgh(i,k) = 0.0_r8
          !ELSE IF(k==plev + 1)THEN
          !   MixLgh(i,k) = 0.0_r8
          !ELSE
          !   MixLgh(i,k) = MIN(MAX(100.0_r8,MixLgh(i,k)),ml2(k))
          !END IF 
          !IF(XLAND(i) >1.0_r8)THEN
          ! OCEAN
          !   kvn = MixLgh(i,k)*SQRT(dvdz2)
          !ELSE
          !   ! LAND
          kvn = ml2(k)*SQRT(dvdz2)
          !END IF
          !
          ! Full diffusion coefficient (modified by f(ri)),
          !
          kvf(i,k+1) = MAX(zkmin,kvn*fstab)

          !kvf(i,k+1)=MIN(gkh1,MAX(gkh0,kvn*fstab)) !(m/sec**2)
       END DO
    END DO
    !
    ! Determine the boundary layer kvh (=kvq), kvm, 
    ! counter gradient terms (cgh, cgq, cgs)
    ! boundary layer height (pblh) and 
    ! the perturbation temperature and moisture (tpert and qpert)
    ! The free atmosphere kv is returned above the boundary layer top.
    !
    CALL pbldif(FRLAND     (1:plond), &
                thm       (1:plond,1:plev)   , &! INTENT(IN   ) ::  th(plond,plev)               ! potential temperature [K]
                qm1        (1:plond,1:plev,1:pcnst), &!INTENT(IN) ::  q(plond,plev,pcnst)     ! specific humidity [kg/kg]
                zm         (1:plond,1:plev), &! INTENT(IN   ) ::  z(plond,plev)               ! height above surface [m]
                um1        (1:plond,1:plev), &! INTENT(IN   ) ::  u(plond,plev)               ! windspeed x-direction [m/s]
                vm1        (1:plond,1:plev), &! INTENT(IN   ) ::  v(plond,plev)               ! windspeed y-direction [m/s]
                tm1        (1:plond,1:plev), &! INTENT(IN   ) ::  t(plond,plev)               ! temperature (used for density)
                pmidm1     (1:plond,1:plev), &! INTENT(IN   ) ::  pmid(plond,plev)        ! midpoint pressures
                kvf        (1:plond,1:plev+1), &! INTENT(IN   ) ::  kvf(plond,plev + 1)     ! free atmospheric eddy diffsvty [m2/s]
                cflx       (1:plond,1:pcnst), &! INTENT(IN   ) ::  cflx(plond,pcnst)       ! surface constituent flux (kg/m2/s)
                shflx      (1:plond), &! INTENT(IN   ) ::  shflx(plond)               ! surface heat flux (W/m2)
                taux       (1:plond), &! INTENT(IN   ) ::  taux(plond)               ! surface u stress (N)
                tauy       (1:plond), &! INTENT(IN   ) ::  tauy(plond)               ! surface v stress (N)
                ustar      (1:plond), &! INTENT(OUT  ) ::  ustar(plond)               ! surface friction velocity [m/s]
                kvm        (1:plond,1:plev+1), &! INTENT(OUT  ) ::  kvm(plond,plev + 1)     ! eddy diffusivity for momentum [m2/s]
                kvh        (1:plond,1:plev+1), &! INTENT(OUT  ) ::  kvh(plond,plev + 1)     ! eddy diffusivity for heat [m2/s]
                TKEMYJ_HSBO(1:plond,1:plev+1), &! INTENT(OUT  ) ::  tke(plond,plev + 1)     ! eddy diffusivity 
                cgh        (1:plond,1:plev+1), &! INTENT(OUT  ) ::  cgh(plond,plev + 1)     ! counter-gradient term for heat [K/m]
                cgq        (1:plond,1:plev+1,1:pcnst), &! INTENT(OUT  ) ::  cgq(plond,plev + 1,pcnst)! counter-gradient term for constituents
                cgs        (1:plond,1:plev+1), &! INTENT(OUT  ) ::  cgs(plond,plev + 1)     ! counter-gradient star (cg/flux)
                PBLH_HSBO  (1:plond), &! INTENT(OUT  ) ::  pblh(plond)               ! boundary-layer height [m]
                tpert      (1:plond), &! INTENT(OUT  ) ::  tpert(plond)               ! convective temperature excess
                qpert      (1:plond), &! INTENT(OUT  ) ::  qpert(plond)               ! convective humidity excess
                ktopbl     (1:plond), &! INTENT(OUT  ) ::  ktopbl(plond)               ! index of first midpoint inside pbl
                ktopblmn   (1:plond), &! INTENT(OUT  ) ::  ktopblmn                ! min value of ktopbl
                plond      , &! INTENT(IN   ) :: plond         ! slt extended domain longitude
                pcnst      , &! INTENT(IN   ) :: pcnst         ! number of constituents (including water vapor)
                plon       , &! INTENT(IN   ) :: plon         ! number of longitudes
                plev       , &! INTENT(IN   ) ::  plev          ! number of vertical levels
                fakn       , &! INTENT(IN   ) :: fakn ! Constant in turbulent prandtl number
                ricr       , &! INTENT(IN   ) :: ricr ! Critical richardson number
                tsk        (1:plond), &! REAL(KIND=r8), INTENT(IN   ) ::  tsk   (plond)    
                qsfc       (1:plond), &! REAL(KIND=r8), INTENT(IN   ) ::  qsfc  (plond)    
                rino       (1:plond,1:plev), &! REAL(KIND=r8), INTENT(INOUT ) ::  rino   (plond,plev)! bulk Richardson no. from level to ref lev
                psomc      (1:plond,1:plev), &! REAL(KIND=r8), INTENT(IN   ) ::  psomc(plond,plev)     ! (psm1/pmidm1)**cappa
                obklen     (1:plond), &! REAL(KIND=r8), INTENT(OUT   ) :: obklen  (plond)! Obukhov length
                phiminv    (1:plond), &! REAL(KIND=r8), INTENT(OUT   ) :: phiminv (plond)! inverse phi function for momentum
                phihinv    (1:plond), &! REAL(KIND=r8), INTENT(OUT   ) :: phihinv (plond)! inverse phi function for heat 
                tstar      (1:plond), &! REAL(KIND=r8), INTENT(OUT   ) :: tstar   (plond)
                wstar      (1:plond)  )! REAL(KIND=r8), INTENT(OUT   ) :: wstar   (plond)


    !

    !
    !     Wgh1=0.333!1.0_r8/3.0_r8 ! pbl Hostlag Boville
    !     Wgh2=0.333!1.0_r8/3.0_r8 ! pbl Mellor Yamada 2.0
    !     Wgh3=0.333!1.0_r8/3.0_r8 ! pbl Mellor Yamada 2.5
    ! 
    !          where kvm(plond,plev + 1) and kvh(plond,plev + 1) are the 
    !     diffusion coefficients for momentum 
    !          and heat, respectively, l is the master turbulence length scale, 
    !          q2 is the turbulent kinetic energy (so q is the magnitude of 
    !          the turbulent wind velocity), and SM and SH are momentum 
    !          flux and heat flux stability parameters, respectively    
    ! 
    DO k=1,plev + 1        
       DO i=1,plon
          TKE(i,k) = TKEMYJ_HSBO(i,k)*Wgh1  +  TKEMYJ_MY20(i,k)*Wgh2  +  TKEMYJ_MY25(i,k)*Wgh3
       END DO
    END DO
    DO i=1,plon
       pblh(i)     = PBLH_HSBO(i)*Wgh1  + PBLH_MY20(i)*Wgh2 + PBLH_MY25(i)*Wgh3
    END DO
    IF (kmean == 1) THEN
       IF (plev >= 4) THEN
              DO k = 2, plev  ! plev + 1
                 DO i = 1, plon
                    !
                    !              k=2  ****Km(k),sl*** } -----------
                    !              k=3/2----si,ric,rf,km,kh,b,l -----------
                    !              k=1  ****Km(k),sl*** } -----------
                    !              k=1/2----si ----------------------------
                    !
                    !        Km(k-1) + 2*Km(k) + Km(k+1)
                    ! Km = -------------------------------
                    !                   4
                    !
                    kvm(i,k)=(0.25_r8*(kvm(i,k-1)+2.0_r8*kvm(i,k)+kvm(i,k+1)))
                    IF(ensamble) THEN
                       IF(k<=plev-1) kvm(i,k)=( kvm(i,k)*Wgh1 + KmMixl(i,k)*Wgh2 + AKM(i,k)*Wgh3 )
                    ELSE
                       IF(k<=plev-1) kvm(i,k)=MAX( kvm(i,k), KmMixl(i,k), AKM(i,k))
                    END IF
                    !AKH
                    !        Kh(k-1) + 2*Kh(k) + Kh(k+1)
                    ! Kh = -------------------------------
                    !                   4
                    !
                    kvh(i,k)=(0.25_r8*(kvh(i,k-1)+2.0_r8*kvh(i,k)+kvh(i,k+1)))
                    IF(ensamble) THEN
                       IF(k<=plev-1) kvh(i,k)=( kvh(i,k)*Wgh1 + KhMixl(i,k)*Wgh2 + AKH(i,k)*Wgh3 )
                    ELSE
                       IF(k<=plev-1) kvh(i,k)=MAX( kvh(i,k),KhMixl(i,k), AKH(i,k))
                    END IF
                 END DO
              END DO
       END IF
       DO i = 1, plon
              !
              !          Km(1) + Km(2)
              ! Km = ---------------
              !               2
              !
              !
              !          Kh(1) + Kh(2)
              ! Kh = ---------------
              !               2
              !
              IF(ensamble) THEN
                 kvm(i,     1  )=       ( kvm(i,     1)*Wgh1 + KmMixl(i,     1)*Wgh2 + AKM(i,     1)*Wgh3 )
                 kvh(i,     1  )=       ( kvh(i,     1)*Wgh1 + KhMixl(i,     1)*Wgh2 + AKH(i,     1)*Wgh3 )
              ELSE
                 kvm(i,     1  )=       MAX( kvm(i,     1),KmMixl(i,     1), AKM(i,     1))
                 kvh(i,     1  )=       MAX( kvh(i,     1),KhMixl(i,     1), AKH(i,     1))
              END IF
              !
              !          Km(k-1) + 2*Km(k) + Km(k+1)
              ! Km = -------------------------------
              !                     4
              !
              !kvm(i,plev  )=(0.25_r8*(kvm(i,plev+1)+2.0_r8*kvm(i,plev  )+kvm(i,plev-1)))
              kvm(i,plev+1)=(0.25_r8*(kvm(i,plev-1)+2.0_r8*kvm(i,plev+1)+kvm(i,plev  )))
              !
              !          Kh(k-1) + 2*Kh(k) + Kh(k+1)
              ! Kh = -------------------------------
              !                     4
              !
              !kvh(i,plev  )=(0.25_r8*(kvh(i,plev+1)+2.0_r8*kvh(i,plev  )+kvh(i,plev-1)))
              kvh(i,plev+1)=(0.25_r8*(kvh(i,plev-1)+2.0_r8*kvh(i,plev+1)+kvh(i,plev  )))
       END DO

    END IF

    IF(PBLEntrain) THEN
    CALL PBLNASA(&
                 plond                            , &           !INTEGER         , INTENT(IN   ) :: IM
                 plev                             , &           !INTEGER         , INTENT(IN   ) :: LM
                 pcnst                            , &           !INTEGER         , INTENT(IN   ) :: pcnst
                 ztodt                            , &           !REAL(KIND=r8), INTENT(IN   ) :: DT
                 kvm        (1:plond,1:plev + 1)  , &           !REAL(KIND=r8), INTENT(INOUT) :: kvm        (IM,LM + 1) 
                 kvh        (1:plond,1:plev + 1)  , &           !REAL(KIND=r8), INTENT(INOUT) :: kvh        (IM,LM + 1 )
                 psomc      (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: psomc      (IM,LM)        !(a
                 zm         (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: zm        (IM,LM)                  ! u wind input 
                 um1        (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: um1       (IM,LM)                  ! u wind input
                 vm1        (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: vm1       (IM,LM)                  ! v wind input
                 tm1        (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: tm1       (IM,LM)                  ! temperature input
                 qm1        (1:plond,1:plev,1:pcnst), &         !REAL(KIND=r8), INTENT(IN   ) :: qm1       (IM,LM,pcnst)  ! moisture and trace constituent input
                 pmidm1     (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: pmidm1    (IM,LM)               ! midpoint pressures
                 pintm1     (1:plond,1:plev + 1)  , &           !REAL(KIND=r8), INTENT(IN   ) :: pintm1    (IM,LM + 1)    ! interface pressures
                 LwCoolRateC(1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: LwCoolRateC     (IM,LM) !clearsky_air_temperature_tendency_lw K s-1
                 LwCoolRate (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: LwCoolRate      (IM,LM)!                 air_temperature_tendency_due_to_longwave K s-1
                 cldtot     (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: cldtot    (IM,LM)             
                 qliq       (1:plond,1:plev)      , &           !REAL(KIND=r8), INTENT(IN   ) :: qliq      (IM,LM)             
                 bstar      (1:plond)             , &           !REAL(KIND=r8), INTENT(IN   ) :: bstar     (IM)   !surface_bouyancy_scale m s-2
                 ustar      (1:plond)             , &           !REAL(KIND=r8), INTENT(IN   ) :: USTAR      (IM) !surface_velocity_scale m s-1
                 FRLAND     (1:plond)                )          !REAL(KIND=r8), INTENT(IN   ) :: FRLAND     (IM)
        END IF         
    !
    ! Add the counter grad terms to potential temp, specific humidity
    ! and other constituents in the bdry layer. Note, ktopblmn gives the 
    ! minimum vertical index of the first midpoint within the boundary layer.
    !
    ! first set values above boundary layer
    !
    DO k=1,plev
       DO i=1,plon
          IF(k <= ktopblmn(i)-2)THEN
             thx(i,k)   = thm(i,k)
             qmx(i,k,1) = qm1(i,k,1)
          END IF
       END DO
       DO m=2,pcnst
          DO i=1,plon
             IF(k <= ktopblmn(i)-2)THEN
                qmx(i,k,m) = qm1(i,k,m)
             END IF
          END DO
       END DO
    END DO
    
    DO k=2,plev
       DO i=1,plon
          potbar(i,k) = pintm1(i,k)/(0.5_r8*(tm1(i,k) + tm1(i,k-1)))
       END DO
    END DO
    DO i=1,plon
       potbar(i,plev + 1) = pintm1(i,plev + 1)/tm1(i,plev)
    END DO
    !
    ! now focus on the boundary layer
    !
    ztodtgor = ztodt*gravit/rair
    DO k=1,plev
       DO i=1,plon
          IF(k >= ktopblmn(i)-1 ) THEN
             tmp1(i   ) = ztodtgor*rpdel(i,k)
             thx(i,k  ) = thm(i,k  ) + tmp1(i) * (potbar(i,k+1)*kvh(i,k+1)*cgh(i,k+1  ) - potbar(i,k  )*kvh(i,k  )*cgh(i,k    ))
             qmx(i,k,1) = qm1(i,k,1) + tmp1(i) * (potbar(i,k+1)*kvh(i,k+1)*cgq(i,k+1,1) - potbar(i,k  )*kvh(i,k  )*cgq(i,k  ,1))
          END IF 
       END DO
       DO m=2,pcnst
          DO i=1,plon
             IF(k >= ktopblmn(i)-1 ) THEN
                qmx(i,k,m) = qm1(i,k,m) + tmp1(i) * (potbar(i,k+1)*kvh(i,k+1)*cgq(i,k+1,m) - potbar(i,k  )*kvh(i,k  )*cgq(i,k  ,m))
             END IF
          END DO
       END DO
    END DO
    !
    ! Check for neg q's in each constituent and put the original vertical
    ! profile back if a neg value is found. A neg value implies that the
    ! quasi-equilibrium conditions assumed for the countergradient term are
    ! strongly violated.
    ! Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.
    !
    DO m=1,pcnst
       DO i=1,plon
          ilogic(i) = 0
       END DO
       DO k=1,plev
          DO i=1,plon
             IF(k>=ktopblmn(i)-1)THEN
                IF (qmx(i,k,m).LT.qmincg(m)) ilogic(i) = 1
             END IF
          END DO
       END DO
       !
       ! Find long indices of those columns for which negatives were found
       !
       CALL wheneq(plon  , &! INTENT(in)
                   ilogic, &! INTENT(in)
                   1     , &! INTENT(in)
                   1     , &! INTENT(in)
                   indx  , &! INTENT(out)
                   nval    )! INTENT(out)
       !
       ! Replace those columns with original values
       !
       IF (nval.GT.0) THEN
          DO k=1,plev
             DO ii=1,nval
                IF(k >= ktopblmn(indx(ii))-1)THEN
                   i=indx(ii)
                   qmx(i,k,m) = qm1(i,k,m)
                END IF
             END DO
          END DO
       END IF
    END DO
    !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
    ! tridiagonal diffusion matrix. the diagonal elements are a combination of 
    ! ca and cc; they are not explicitly provided to the solver
    !
    gorsq = (gravit/rair)**2
    DO k=ntopfl,plev-1
       DO i=1,plon
          tmp2 = ztodt*gorsq*rpdeli(i,k)*(potbar(i,k+1)**2)
          cah(i,k  ) = kvh(i,k+1)*tmp2*rpdel(i,k  )
          cam(i,k  ) = kvm(i,k+1)*tmp2*rpdel(i,k  )
          cch(i,k+1) = kvh(i,k+1)*tmp2*rpdel(i,k+1)
          ccm(i,k+1) = kvm(i,k+1)*tmp2*rpdel(i,k+1)
       END DO
    END DO
    !
    ! The last element of the upper diagonal is zero.
    !
    DO i=1,plon
       cah(i,plev) = 0.0_r8
       cam(i,plev) = 0.0_r8
    END DO
    !
    ! Calculate e(k) for heat & momentum vertical diffusion.  This term is 
    ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.
    !
    DO i=1,plon
       termh(i,ntopfl) = 1.0_r8/(1.0_r8 + cah(i,ntopfl))
       termm(i,ntopfl) = 1.0_r8/(1.0_r8 + cam(i,ntopfl))
       zeh  (i,ntopfl) = cah(i,ntopfl)*termh(i,ntopfl)
       zem  (i,ntopfl) = cam(i,ntopfl)*termm(i,ntopfl)
    END DO
    DO k=ntopfl+1,plev-1
       DO i=1,plon
          termh(i,k) = 1.0_r8/ &
               (1.0_r8 + cah(i,k) + cch(i,k) - cch(i,k)*zeh(i,k-1))
          termm(i,k) = 1.0_r8/ &
               (1.0_r8 + cam(i,k) + ccm(i,k) - ccm(i,k)*zem(i,k-1))
          zeh(i,k) = cah(i,k)*termh(i,k)
          zem(i,k) = cam(i,k)*termm(i,k)
       END DO
    END DO
    !
    ! Diffuse momentum
    !

    CALL mvdiff(um1     , & ! INTENT(IN  ) :: um1(plond,plev)             ! u wind input
                vm1     , & ! INTENT(IN  ) :: vm1(plond,plev)             ! v wind input
                dubot   , & ! INTENT(IN  ) :: uflx(plond)    ! sfc u flux into lowest model level
                dvbot   , & ! INTENT(IN  ) :: vflx(plond)    ! sfc v flux into lowest model level
                sufac   , & ! INTENT(IN  ) :: sufac(plond)   ! lowest layer u implicit stress factor
                svfac   , & ! INTENT(IN  ) :: svfac(plond)   ! lowest layer v implicit stress factor
                ccm     , & ! INTENT(IN  ) :: cc(plond,plev) ! -lower diag coeff. of tri-diag matrix
                zem     , & ! INTENT(IN  ) :: ze(plond,plev) ! term in tri-diag. matrix system
                termm   , & ! INTENT(IN  ) :: term(plond,plev)  ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                up1     , & ! INTENT(OUT ) :: up1(plond,plev)        ! u wind after diffusion
                vp1     , & ! INTENT(OUT ) :: vp1(plond,plev)        ! v wind after diffusion
                plond   , & ! INTENT(IN  ) :: plond         ! slt extended domain longitude
                plev    , & ! INTENT(IN  ) :: plev    ! number of vertical levels
                plon      ) ! INTENT(IN  ) :: plon    ! number of longitudes
    !+
    ! Determine the difference between the implicit stress and the
    ! externally specified stress. Apply over boundary layer depth.
    !-
   
    DO i = 1, plon
       IF (ksx(i) .GT. 0.0_r8) sufac(i) = (taux(i) + ksx(i)*up1(i,plev)) &
            * ztodt * gravit / (pintm1(i,plev + 1) - pintm1(i,ktopbl(i)))
       IF (ksy(i) .GT. 0.0_r8) svfac(i) = (tauy(i) + ksy(i)*vp1(i,plev)) &
            * ztodt * gravit / (pintm1(i,plev + 1) - pintm1(i,ktopbl(i)))
       !$$$         print 1, i, lat
       !$$$     $        , taux(i), -ksx(i)*up1(i,plev)
       !$$$     $        , sufac(i)*86400.0_r8/ztodt, ksx(i), um1(i,plev),up1(i,plev)
       !$$$    1 format (1x,i4,i4, 1P, 4E20.10_r8, 2f10.4)
    END DO
    DO k = plev, 1, -1
       DO i = 1, plon
          IF (k .GE. ktopbl(i) .and. k >= ktopblmn(i)) THEN
             up1(i,k) = up1(i,k) + sufac(i)
             vp1(i,k) = vp1(i,k) + svfac(i)
          END IF
       END DO
    END DO
    !
    ! Diffuse constituents
    !

    CALL qvdiff(pcnst   , &! INTENT(IN   ) :: ncnst                     ! number of constituents being diffused
                qmx     , &! INTENT(IN   ) :: qm1   (plond,plev,ncnst)  ! initial constituent
                dqbot   , &! INTENT(IN   ) :: qflx  (plond,ncnst)      ! sfc q flux into lowest model level
                cch     , &! INTENT(IN   ) :: cc    (plond,plev)         ! -lower diag coeff.of tri-diag matrix
                zeh     , &! INTENT(IN   ) :: ze    (plond,plev)         ! term in tri-diag. matrix system
                termh   , &! INTENT(IN   ) :: term  (plond,plev)       ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                qp1     , &! INTENT(OUT  ) :: qp1   (plond,plev,ncnst)  ! final constituent
                plon    , &! INTENT(IN   ) :: plon   ! number of longitudes
                plond   , &! INTENT(IN   ) :: plond  ! slt extended domain longitude
                plev    , &! INTENT(IN   ) :: plev   ! number of vertical levels
                pcnst   , &! INTENT(IN   ) :: pcnst  ! number of constituents (including water vapor)
                ntopfl   ) ! INTENT(IN   ) :: ntopfl ! Top level to which vertical diffusion is applied.
    !
    ! Identify and correct constituents exceeding user defined bounds
    ! 
    CALL qneg3(plond     , & ! INTENT(IN   ) :: plond           ! slt extended domain longitude
               plev      , & ! INTENT(IN   ) :: plev           ! number of vertical levels
               pcnst     , & ! INTENT(IN   ) :: pcnst           ! number of constituents (including water vapor)
               plon      , & ! INTENT(IN   ) :: plon           ! number of longitudes
               'VDIFF   ', & ! INTENT(IN   ) :: subnam       ! name of calling routine
               qp1(1,1,1)  ) ! INTENT(INOUT) :: q(plond,plev,pcnst) ! moisture/tracer field
    !
    ! Diffuse potential temperature
    !
    CALL qvdiff(1       , & ! INTENT(IN   ) :: ncnst                      ! number of constituents being diffused
                thx     , & ! INTENT(IN   ) :: qm1   (plond,plev,ncnst)  ! initial constituent
                dtbot   , & ! INTENT(IN   ) :: qflx  (plond,ncnst)        ! sfc q flux into lowest model level
                cch     , & ! INTENT(IN   ) :: cc    (plond,plev)          ! -lower diag coeff.of tri-diag matrix
                zeh     , & ! INTENT(IN   ) :: ze    (plond,plev)          ! term in tri-diag. matrix system
                termh   , & ! INTENT(IN   ) :: term  (plond,plev)        ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                thp     , & ! INTENT(OUT  ) :: qp1   (plond,plev,ncnst)  ! final constituent
                plon    , & ! INTENT(IN   ) :: plon   ! number of longitudes
                plond   , & ! INTENT(IN   ) :: plond  ! slt extended domain longitude
                plev    , & ! INTENT(IN   ) :: plev   ! number of vertical levels
                pcnst   , & ! INTENT(IN   ) :: pcnst  ! number of constituents (including water vapor)
                ntopfl   )  ! INTENT(IN   ) :: ntopfl ! Top level to which vertical diffusion is applied.
    !
    RETURN
  END SUBROUTINE vdiff



  !===============================================================================

  SUBROUTINE wheneq (n, array, inc, TARGET, index, nval)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Determine indices of "array" which equal "target"
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE
    !
    ! Arguments
    !
    INTEGER, INTENT(in) :: n

    INTEGER, INTENT(in) :: array(n)    ! array to be searched
    INTEGER, INTENT(in) :: TARGET      ! value to compare against
    INTEGER, INTENT(in) :: inc         ! increment to move through array

    INTEGER, INTENT(out) :: nval       ! number of values meeting criteria
    INTEGER, INTENT(out) :: INDEX(n)   ! output array of indices
    !
    ! Local workspace
    !
    INTEGER :: i
    INTEGER :: ina

    ina=1
    nval=0
    INDEX=0
    IF (inc .LT. 0) ina=(-inc)*(n-1)+1
    DO i=1,n
       IF(array(ina) .EQ. TARGET) THEN
          nval=nval+1
          INDEX(nval)=i
       END IF
       ina=ina+inc
    ENDDO
    RETURN
  END SUBROUTINE wheneq

  SUBROUTINE virtem(&
                   t       , &! REAL(KIND=r8), INTENT(IN   ) ::  t(plond,plev)       ! temperature
                   q       , &! REAL(KIND=r8), INTENT(IN   ) ::  q(plond,plev)       ! specific humidity
                   zvir    , &! REAL(KIND=r8), INTENT(IN   ) ::  zvir                ! virtual temperature constant
                   tv      , &! REAL(KIND=r8), INTENT(OUT  ) ::  tv(plond,plev)      ! virtual temperature
                   plond   , &! INTEGER, INTENT(IN   ) ::  plond      ! slt extended domain longitude
                   plev    , &! INTEGER, INTENT(IN   ) ::  plev       ! number of vertical levels
                   plon      )! INTEGER, INTENT(IN   ) ::  plon       ! number of longitudes
    !-----------------------------------------------------------------------
    !
    ! Compute the virtual temperature.
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  B. Boville
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, J. Hack, August 1992
    ! Reviewed:          D. Williamson, March 1996
    !
    !-----------------------------------------------------------------------
    !
    ! $Id: virtem.F,v 1.1.1.1 2001/03/09 00:29:29 mirin Exp $
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) ::  plond      ! slt extended domain longitude
    INTEGER, INTENT(IN   ) ::  plev       ! number of vertical levels
    INTEGER, INTENT(IN   ) ::  plon       ! number of longitudes

    REAL(KIND=r8), INTENT(IN   ) ::  t(plond,plev)       ! temperature
    REAL(KIND=r8), INTENT(IN   ) ::  q(plond,plev)       ! specific humidity
    REAL(KIND=r8), INTENT(IN   ) ::  zvir                ! virtual temperature constant
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(OUT  ) ::  tv(plond,plev)      ! virtual temperature
    !
    !---------------------------Local storage-------------------------------
    !
    INTEGER i,k              ! longitude and level indexes
    !
    tv=0.0_r8
    DO k=1,plev
       DO i=1,plon
          tv(i,k) = t(i,k)*(1.0_r8 + zvir*q(i,k))
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE virtem



  SUBROUTINE mvdiff(um1     , &! INTENT(IN  ) :: um1(plond,plev)        ! u wind input
                    vm1     , &! INTENT(IN  ) :: vm1(plond,plev)        ! v wind input
                    uflx    , &! INTENT(IN  ) :: uflx(plond)        ! sfc u flux into lowest model level
                    vflx    , &! INTENT(IN  ) :: vflx(plond)        ! sfc v flux into lowest model level
                    sufac   , &! INTENT(IN  ) :: sufac(plond)        ! lowest layer u implicit stress factor
                    svfac   , &! INTENT(IN  ) :: svfac(plond)        ! lowest layer v implicit stress factor
                    cc      , &! INTENT(IN  ) :: cc(plond,plev)        ! -lower diag coeff. of tri-diag matrix
                    ze      , &! INTENT(IN  ) :: ze(plond,plev)        ! term in tri-diag. matrix system
                    term    , &! INTENT(IN  ) :: term(plond,plev)  ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                    up1     , &! INTENT(OUT ) :: up1(plond,plev)   ! u wind after diffusion
                    vp1     , &! INTENT(OUT ) :: vp1(plond,plev)   ! v wind after diffusion
                    plond   , &! INTENT(IN  ) :: plond      ! slt extended domain longitude
                    plev    , &! INTENT(IN  ) :: plev         ! number of vertical levels
                    plon      )! INTENT(IN  ) :: plon         ! number of longitudes
    !-----------------------------------------------------------------------
    !
    ! Vertical momentum diffusion with explicit surface flux.
    ! Solve the vertical diffusion equation for momentum.
    ! Procedure for solution of the implicit equation follows
    ! Richtmyer and Morton (1967,pp 198-199)
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  CCM1
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          P. Rasch, B. Boville, August 1992
    ! Reviewed:          P. Rasch, April 1996
    ! Reviewed:          B. Boville, May 1996
    !
    !-----------------------------------------------------------------------
    !
    ! $Id: mvdiff.F,v 1.1.1.1 2001/03/09 00:29:27 mirin Exp $
    ! $Author: mirin $
    !
    !-----------------------------------------------------------------------
    !      use precision
    !      use pmgrid
    !-----------------------------------------------------------------------
    !!#include <implicit.h>
    !------------------------------Commons----------------------------------
    !!#include <comvd.h>
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: plond      ! slt extended domain longitude
    INTEGER, INTENT(IN   ) :: plev         ! number of vertical levels
    INTEGER, INTENT(IN   ) :: plon         ! number of longitudes

    REAL(KIND=r8), INTENT(IN  ) :: um1(plond,plev)        ! u wind input
    REAL(KIND=r8), INTENT(IN  ) :: vm1(plond,plev)        ! v wind input
    REAL(KIND=r8), INTENT(IN  ) :: uflx(plond)        ! sfc u flux into lowest model level
    REAL(KIND=r8), INTENT(IN  ) :: vflx(plond)        ! sfc v flux into lowest model level
    REAL(KIND=r8), INTENT(IN  ) :: sufac(plond)        ! lowest layer u implicit stress factor
    REAL(KIND=r8), INTENT(IN  ) :: svfac(plond)        ! lowest layer v implicit stress factor
    REAL(KIND=r8), INTENT(IN  ) :: cc(plond,plev)        ! -lower diag coeff. of tri-diag matrix
    REAL(KIND=r8), INTENT(IN  ) :: ze(plond,plev)        ! term in tri-diag. matrix system
    REAL(KIND=r8), INTENT(IN  ) :: term(plond,plev)  ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(OUT ) :: up1(plond,plev)   ! u wind after diffusion
    REAL(KIND=r8), INTENT(OUT ) :: vp1(plond,plev)   ! v wind after diffusion
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(KIND=r8) zfu(plond,plev)   ! terms appearing in soln of tridiag system
    REAL(KIND=r8) zfv(plond,plev)   ! terms appearing in soln of tridiag system
    INTEGER i,k            ! longitude,vertical indices
    !
    !-----------------------------------------------------------------------
    !
    ! Calc fu(k) and fv(k). These are terms required in solution of 
    ! tridiagonal matrix defined by implicit diffusion eqn.  Note that only 
    ! levels ntopfl through plev need be solved for. No diffusion is 
    ! applied above this level.
    !
    up1=0.0_r8
    vp1=0.0_r8
    zfu=0.0_r8
    zfv=0.0_r8
    DO i=1,plon
       zfu(i,ntopfl) = um1(i,ntopfl)*term(i,ntopfl)
       zfv(i,ntopfl) = vm1(i,ntopfl)*term(i,ntopfl)
    END DO
    DO k=ntopfl+1,plev-1
       DO i=1,plon
          zfu(i,k) = (um1(i,k) + cc(i,k)*zfu(i,k-1))*term(i,k)
          zfv(i,k) = (vm1(i,k) + cc(i,k)*zfv(i,k-1))*term(i,k)
       END DO
    END DO
    !
    ! Bottom level: (includes  surface fluxes as either an explicit RHS or
    ! as an implicit stress)
    !
    DO i=1,plon
       zfu(i,plev) = (um1(i,plev) + uflx(i) + cc(i,plev)*zfu(i,plev-1)) &
            / (1.0_r8 + cc(i,plev) + sufac(i) - cc(i,plev)*ze(i,plev-1))
       zfv(i,plev) = (vm1(i,plev) + vflx(i) + cc(i,plev)*zfv(i,plev-1)) &
            / (1.0_r8 + cc(i,plev) + svfac(i) - cc(i,plev)*ze(i,plev-1))
    END DO
    !
    ! Perform back substitution
    !
    DO i=1,plon
       up1(i,plev) = zfu(i,plev)
       vp1(i,plev) = zfv(i,plev)
    END DO
    DO k=plev-1,ntopfl,-1
       DO i=1,plon
          up1(i,k) = zfu(i,k) + ze(i,k)*up1(i,k+1)
          vp1(i,k) = zfv(i,k) + ze(i,k)*vp1(i,k+1)
       END DO
    END DO
    RETURN
  END SUBROUTINE mvdiff


  SUBROUTINE qneg3(plond  , &! INTENT(IN   ) :: plond      ! slt extended domain longitude
                   plev   , &! INTENT(IN   ) :: plev       ! number of vertical levels
                   pcnst  , &! INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)
                   plon   , &! INTENT(IN   ) :: plon       ! number of longitudes
                   subnam , &! INTENT(IN   ) :: subnam       ! name of calling routine
                   q       ) ! INTENT(INOUT) :: q(plond,plev,pcnst) ! moisture/tracer field
    !-----------------------------------------------------------------------
    !
    ! Check moisture and tracers for minimum value, reset any below
    ! minimum value to minimum value and return information to allow
    ! warning message to be printed. The global average is NOT preserved.
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  J. Rosinski
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, August 1992, March 1996
    !
    !-----------------------------------------------------------------------
    !
    ! $Id: qneg3.F,v 1.1.1.1 2001/03/09 00:29:27 mirin Exp $
    !
    !-----------------------------------------------------------------------
    !      use precision
    !      use pmgrid
    !-----------------------------------------------------------------------
    !#include <implicit.h>
    !------------------------------Commons----------------------------------
    !#include <comqmin.h>
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: plon       ! number of longitudes
    INTEGER, INTENT(IN   ) :: plond      ! slt extended domain longitude
    INTEGER, INTENT(IN   ) :: plev       ! number of vertical levels
    INTEGER, INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)

    CHARACTER*8, INTENT(IN   ) :: subnam       ! name of calling routine
    !
    ! Input/Output arguments
    !
    REAL(KIND=r8), INTENT(INOUT) :: q(plond,plev,pcnst) ! moisture/tracer field
    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER indx(plond)      ! array of indices of points < qmin
    INTEGER nval             ! number of points < qmin for 1 level
    INTEGER nvals            ! number of values found < qmin
    INTEGER i,ii,k           ! longitude, level indices
    INTEGER m                ! constituent index
    INTEGER iw,kw            ! i,k indices of worst violator

    LOGICAL found            ! true => at least 1 minimum violator found

    REAL(KIND=r8) worst               ! biggest violator
    !
    !-----------------------------------------------------------------------
    !
    DO m=1,pcnst
       DO k=1,plev
          DO i=1,plon
             q(i,k,m)=max(q(i,k,m),qmin(m))
          END DO   
       END DO          
    END DO
    RETURN
    
    DO m=1,pcnst
       nvals = 0
       found = .FALSE.
       worst = 0.0_r8
       !
       ! Test all field values for being less than minimum value. Set q = qmin
       ! for all such points. Trace offenders and identify worst one.
       !
       DO k=1,plev
          nval = 0
          DO i=1,plon
             IF (q(i,k,m) < qmin(m)) THEN
                nval = nval + 1
                indx(nval) = i
             END IF
          END DO

          IF (nval.GT.0) THEN
             found = .TRUE.
             nvals = nvals + nval
             DO ii=1,nval
                i = indx(ii)
                IF (q(i,k,m).LT.worst) THEN
                   !worst = q(i,k,m)
                   kw = k
                   iw = i
                END IF
                q(i,k,m) = qmin(m)
             END DO
          END IF
       END DO
       IF (found) THEN
          WRITE(6,9000)subnam,m,nvals,qmin(m),worst,iw,kw
       END IF
    END DO
    !
    RETURN
9000 FORMAT(' QNEG3 from ',a8,':m=',i3, &
         ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
         1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
  END SUBROUTINE qneg3



  SUBROUTINE pbldif(FRLAND  , &
                    th      , &! INTENT(IN   ) ::  th(plond,plev)          ! potential temperature [K]
                    q       , &! INTENT(IN   ) ::  q(plond,plev,pcnst)     ! specific humidity [kg/kg]
                    z       , &! INTENT(IN   ) ::  z(plond,plev)           ! height above surface [m]
                    u       , &! INTENT(IN   ) ::  u(plond,plev)           ! windspeed x-direction [m/s]
                    v       , &! INTENT(IN   ) ::  v(plond,plev)           ! windspeed y-direction [m/s]
                    t       , &! INTENT(IN   ) ::  t(plond,plev)           ! temperature (used for density)
                    pmid    , &! INTENT(IN   ) ::  pmid(plond,plev)        ! midpoint pressures
                    kvf     , &! INTENT(IN   ) ::  kvf(plond,plev + 1)     ! free atmospheric eddy diffsvty [m2/s]
                    cflx    , &! INTENT(IN   ) ::  cflx(plond,pcnst)           ! surface constituent flux (kg/m2/s)
                    shflx   , &! INTENT(IN   ) ::  shflx(plond)            ! surface heat flux (W/m2)
                    taux    , &! INTENT(IN   ) ::  taux(plond)             ! surface u stress (N)
                    tauy    , &! INTENT(IN   ) ::  tauy(plond)             ! surface v stress (N)
                    ustar   , &! INTENT(OUT  ) ::  ustar(plond)            ! surface friction velocity [m/s]
                    kvm     , &! INTENT(OUT  ) ::  kvm(plond,plev + 1)     ! eddy diffusivity for momentum [m2/s]
                    kvh     , &! INTENT(OUT  ) ::  kvh(plond,plev + 1)     ! eddy diffusivity for heat [m2/s]
                    tke     , &! INTENT(OUT  ) ::  tke(plond,plev + 1)     ! eddy diffusivity 
                    cgh     , &! INTENT(OUT  ) ::  cgh(plond,plev + 1)     ! counter-gradient term for heat [K/m]
                    cgq     , &! INTENT(OUT  ) ::  cgq(plond,plev + 1,pcnst)! counter-gradient term for constituents
                    cgs     , &! INTENT(OUT  ) ::  cgs(plond,plev + 1)     ! counter-gradient star (cg/flux)
                    pblh    , &! INTENT(OUT  ) ::  pblh(plond)             ! boundary-layer height [m]
                    tpert   , &! INTENT(OUT  ) ::  tpert(plond)            ! convective temperature excess
                    qpert   , &! INTENT(OUT  ) ::  qpert(plond)            ! convective humidity excess
                    ktopbl  , &! INTENT(OUT  ) ::  ktopbl(plond)           ! index of first midpoint inside pbl
                    ktopblmn, &! INTENT(OUT  ) ::  ktopblmn                ! min value of ktopbl
                    plond   , &! INTENT(IN   ) :: plond      ! slt extended domain longitude
                    pcnst   , &! INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)
                    plon    , &! INTENT(IN   ) :: plon       ! number of longitudes
                    plev    , &! INTENT(IN   ) ::  plev       ! number of vertical levels
                    fakn    , &! INTENT(IN   ) :: fakn    ! Constant in turbulent prandtl number
                    ricr    , &! INTENT(IN   ) :: ricr    ! Critical richardson number
                    tsk     , &! REAL(KIND=r8), INTENT(IN   ) ::  tsk   (plond)    
                    qsfc    , &! REAL(KIND=r8), INTENT(IN   ) ::  qsfc  (plond)    
                    rino    , &! REAL(KIND=r8), INTENT(INOUT ) ::  rino   (plond,plev)        ! bulk Richardson no. from level to ref lev
                    psomc   , &! REAL(KIND=r8), INTENT(IN   ) ::  psomc(plond,plev)! (psm1/pmidm1)**cappa
                    obklen  , &! REAL(KIND=r8), INTENT(OUT   ) :: obklen  (plond)           ! Obukhov length
                    phiminv , &! REAL(KIND=r8), INTENT(OUT   ) :: phiminv (plond)  ! inverse phi function for momentum
                    phihinv , &! REAL(KIND=r8), INTENT(OUT   ) :: phihinv (plond)  ! inverse phi function for heat 
                    tstar   , &! REAL(KIND=r8), INTENT(OUT   ) :: tstar   (plond)
                    wstar      )! REAL(KIND=r8), INTENT(OUT   ) :: wstar   (plond)
    !------------------------------------------------------------------------
    ! 
    ! Atmospheric boundary layer computation.
    !
    ! Nonlocal scheme that determines eddy diffusivities based on a
    ! diagnosed boundary layer height and a turbulent velocity scale;
    ! also, countergradient effects for heat and moisture, and constituents
    ! are included, along with temperature and humidity perturbations which 
    ! measure the strength of convective thermals in the lower part of the 
    ! atmospheric boundary layer.
    !
    ! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
    ! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
    ! Model. J. Clim., vol. 6., p. 1825--1842.
    !
    ! Updated by Holtslag and Hack to exclude the surface layer from the
    ! definition of the boundary layer Richardson number. Ri is now defined
    ! across the outer layer of the pbl (between the top of the surface
    ! layer and the pbl top) instead of the full pbl (between the surface and
    ! the pbl top). For simiplicity, the surface layer is assumed to be the
    ! region below the first model level (otherwise the boundary layer depth 
    ! determination would require iteration).
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  B. Boville
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          B. Boville, P. Rasch, August 1992
    ! Reviewed:          B. Boville, P. Rasch, April 1996
    !
    ! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
    ! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
    !
    !-----------------------------------------------------------------------
    !
    ! $Id: pbldif.F,v 1.1.1.1 2001/03/09 00:29:27 mirin Exp $
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: plond      ! slt extended domain longitude
    INTEGER, INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)
    INTEGER, INTENT(IN   ) :: plon       ! number of longitudes
    INTEGER, INTENT(IN   ) ::  plev       ! number of vertical levels
    REAL(KIND=r8), INTENT(IN   ) :: fakn    ! Constant in turbulent prandtl number
    REAL(KIND=r8), INTENT(IN   ) :: ricr    ! Critical richardson number
    REAL(KIND=r8), INTENT(IN   ) ::  th(plond,plev)          ! potential temperature [K]
    REAL(KIND=r8), INTENT(IN   ) ::  q(plond,plev,pcnst)     ! specific humidity [kg/kg]
    REAL(KIND=r8), INTENT(IN   ) ::  z(plond,plev)           ! height above surface [m]
    REAL(KIND=r8), INTENT(IN   ) ::  u(plond,plev)           ! windspeed x-direction [m/s]
    REAL(KIND=r8), INTENT(IN   ) ::  v(plond,plev)           ! windspeed y-direction [m/s]
    REAL(KIND=r8), INTENT(IN   ) ::  t(plond,plev)           ! temperature (used for density)
    REAL(KIND=r8), INTENT(IN   ) ::  pmid(plond,plev)        ! midpoint pressures
    REAL(KIND=r8), INTENT(IN   ) ::  kvf(plond,plev + 1)        ! free atmospheric eddy diffsvty [m2/s]
    REAL(KIND=r8), INTENT(IN   ) ::  cflx(plond,pcnst)       ! surface constituent flux (kg/m2/s)
    REAL(KIND=r8), INTENT(IN   ) ::  shflx(plond)            ! surface heat flux (W/m2)
    REAL(KIND=r8), INTENT(IN   ) ::  taux(plond)             ! surface u stress (N)
    REAL(KIND=r8), INTENT(IN   ) ::  tauy(plond)             ! surface v stress (N)
    REAL(KIND=r8), INTENT(IN   ) ::  tsk   (plond)    
    REAL(KIND=r8), INTENT(IN   ) ::  qsfc  (plond)    
    REAL(KIND=r8), INTENT(IN   ) ::  psomc(plond,plev)      ! (psm1/pmidm1)**cappa
    REAL(KIND=r8), INTENT(IN   ) ::  FRLAND (plond)    
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(OUT   ) ::  ustar  (plond)            ! surface friction velocity [m/s]
    REAL(KIND=r8), INTENT(OUT   ) ::  kvm    (plond,plev + 1)        ! eddy diffusivity for momentum [m2/s]
    REAL(KIND=r8), INTENT(OUT   ) ::  kvh    (plond,plev + 1)        ! eddy diffusivity for heat [m2/s]
    REAL(KIND=r8), INTENT(OUT   ) ::  cgh    (plond,plev + 1)        ! counter-gradient term for heat [K/m]
    REAL(KIND=r8), INTENT(OUT   ) ::  cgq    (plond,plev + 1,pcnst)  ! counter-gradient term for constituents
    REAL(KIND=r8), INTENT(OUT   ) ::  cgs    (plond,plev + 1)        ! counter-gradient star (cg/flux)
    REAL(KIND=r8), INTENT(INOUT ) ::  pblh   (plond)             ! boundary-layer height [m]
    REAL(KIND=r8), INTENT(INOUT ) ::  rino   (plond,plev)        ! bulk Richardson no. from level to ref lev
    REAL(KIND=r8), INTENT(OUT   ) ::  tke    (plond,plev + 1)
    REAL(KIND=r8), INTENT(OUT   ) ::  tpert  (plond)            ! convective temperature excess
    REAL(KIND=r8), INTENT(OUT   ) ::  qpert  (plond)            ! convective humidity excess

    INTEGER      , INTENT(OUT   ) :: ktopbl  (plond)        ! index of first midpoint inside pbl
    INTEGER      , INTENT(OUT   ) :: ktopblmn(plond)             ! min value of ktopbl
    REAL(KIND=r8), INTENT(OUT   ) :: obklen  (plond)           ! Obukhov length
    REAL(KIND=r8), INTENT(OUT   ) :: phiminv (plond)          ! inverse phi function for momentum
    REAL(KIND=r8), INTENT(OUT   ) :: phihinv (plond)          ! inverse phi function for heat 
    REAL(KIND=r8), INTENT(OUT   ) :: tstar   (plond)
    REAL(KIND=r8), INTENT(OUT   ) :: wstar   (plond)

    !
    !---------------------------Local parameters----------------------------
    !
    REAL(KIND=r8) tiny                    ! lower bound for wind magnitude
    PARAMETER (tiny=1.e-36_r8)
    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER i                    ! longitude index
    INTEGER k                    ! level index
    INTEGER m                    ! constituent index
    REAL(KIND=r8) heatv  (plond)            ! surface virtual heat flux
    REAL(KIND=r8) thvsrf (plond)           ! sfc (bottom) level virtual temperature
    REAL(KIND=r8) thvref (plond)           ! reference level virtual temperature
    REAL(KIND=r8) tkv                     ! model level potential temperature
    REAL(KIND=r8) therm  (plond)            ! thermal virtual temperature excess
    REAL(KIND=r8) wm     (plond)               ! turbulent velocity scale for momentum
    REAL(KIND=r8) vvk                     ! velocity magnitude squared
    REAL(KIND=r8) zm     (plond)               ! current level height
    REAL(KIND=r8) zp     (plond)               ! current level height + one level up
    REAL(KIND=r8) khfs   (plond)             ! surface kinematic heat flux [mK/s]
    REAL(KIND=r8) kqfs   (plond,pcnst)       ! sfc kinematic constituent flux [m/s]
    REAL(KIND=r8) zmzp                    ! level height halfway between zm and zp
    REAL(KIND=r8) tlv    (plond)              ! ref. level pot tmp + tmp excess
    REAL(KIND=r8) fak1   (plond)             ! k*ustar*pblh
    REAL(KIND=r8) fak2   (plond)             ! k*wm*pblh
    REAL(KIND=r8) fak3   (plond)             ! fakn*wstr/wm 
    REAL(KIND=r8) pblk   (plond)             ! level eddy diffusivity for momentum
    REAL(KIND=r8) pr     (plond)               ! Prandtl number for eddy diffusivities
    REAL(KIND=r8) zl     (plond)               ! zmzp / Obukhov length
    REAL(KIND=r8) zh     (plond)               ! zmzp / pblh
    REAL(KIND=r8) zzh    (plond)              ! (1-(zmzp/pblh))**2
    REAL(KIND=r8) wstr   (plond)             ! w*, convective velocity scale
    REAL(KIND=r8) rrho   (plond)             ! 1./bottom level density (temporary)
    REAL(KIND=r8) ustr                    ! unbounded ustar
    REAL(KIND=r8) term                    ! intermediate calculation
    REAL(KIND=r8) fac                     ! interpolation factor
    REAL(KIND=r8) FacLand                     ! interpolation factor
    REAL(KIND=r8) pblmin                  ! min pbl height due to mechanical mixing

    LOGICAL unstbl (plond)        ! pts w/unstbl pbl (positive virtual ht flx)
    LOGICAL stblev (plond)        ! stable pbl with levels within pbl
    LOGICAL unslev (plond)        ! unstbl pbl with levels within pbl
    LOGICAL unssrf (plond)        ! unstb pbl w/lvls within srf pbl lyr
    LOGICAL unsout (plond)        ! unstb pbl w/lvls in outer pbl lyr
    LOGICAL check  (plond)         ! True=>chk if Richardson no.>critcal
    ustar=0.0_r8;kvm    =0.0_r8;kvh    =0.0_r8;cgh    =0.0_r8  
    cgq    =0.0_r8;cgs    =0.0_r8;TKE=0.0_r8;tpert=0.0_r8;qpert=0.0_r8  
    obklen  =0.0_r8;phiminv =0.0_r8;    phihinv =0.0_r8  
    heatv =0.0_r8;    thvsrf=0.0_r8;    thvref=0.0_r8;    tkv   =0.0_r8;
    therm =0.0_r8;    wm    =0.0_r8;    vvk   =0.0_r8;    zm    =0.0_r8;
    zp    =0.0_r8;    khfs  =0.0_r8;    kqfs  =0.0_r8;    zmzp  =0.0_r8;
    tlv   =0.0_r8;    fak1  =0.0_r8;    fak2  =0.0_r8;    fak3  =0.0_r8;
    pblk  =0.0_r8;    pr    =0.0_r8;    zl    =0.0_r8;    zh    =0.0_r8;
    zzh   =0.0_r8;    wstr  =0.0_r8;    rrho  =0.0_r8;    ustr  =0.0_r8;
    term  =0.0_r8;    fac   =0.0_r8;    pblmin=0.0_r8;
    unstbl=.true.
    stblev=.true.
    unslev=.true.
    unssrf=.true.
    unsout=.true.
    check =.true.
    !
    ! Compute kinematic surface fluxes
    !
    DO i=1,plon
       !j/kg/kelvin
       !
       ! P = rho * R * T
       !
       !            P
       ! rho  = -------
       !          R * T
       !
       !           1             R * T
       !  rrho = ----- =   -------
       !          rho              P
       !
       !
       rrho(i) = rair*t(i,plev)/pmid(i,plev)

       !
       !
       ! tau = rho * ustar*2
       !
       !
       ustr = SQRT(SQRT(taux(i)**2 + tauy(i)**2)*rrho(i))
       !
       ! surface friction velocity [m/s]
       !
       ustar(i) = MAX(ustr,0.01_r8)

       !
       ! 
       ! surface kinematic heat flux [mK/s] ustar*tstar
       !
       ! H=Ho=-rho*cp*ustar*tstar  
       !
       ! surface heat flux (W/m2)  (j/kg/kelvin)
       !
       !                  H           j       kg K   m*m*m           m K
       ! ustar*tstar = ---------   = ------------------------- =   ------------
       !                -rho*cp      m*m*s   j         kg              s
       !
       khfs(i) = shflx(i)*rrho(i)/cpair
       !
       !
       !
       ! surface kinematic constituent flux [m/s]
       !
       !  J            N *m         Kg * m * m             m*m
       !------ = --------- =------------------- = ---------= hvap = 2.5104e+6! latent heat of vaporization of water (J kg-1)
       !  kg             kg          s*s Kg              s*s
       !
       ! (kg/m2/s)
       !   W            J              N *m         Kg * m * m          s*s              kg 
       ! ------ = ------- = --------- = ------------  * -------  = -------
       !  m*m           m*m*s      m*m*s          s*s*s*m*m          m*m            M*m*s
       !
       !         W           kg        
       !   =   ------ *  ------  =  LFlux / hvap
       !        m*m           J  
       !
       ! E=Eo=-rho*hl*ustar*qstar
       !
       !                  E               kg          m * m * m          m 
       ! ustar*qstar = ------------ = ----------- * ------------- = -----------
       !                -rho *hl       m * m * s       kg                s
       !
       !
       ! surface constituent flux (kg/m2/s)
       !
       kqfs(i,1) = cflx(i,1)*rrho(i)
       !
    END DO
    
    DO m=2,pcnst
       DO i=1,plon
          kqfs(i,m)= cflx(i,m)*rrho(i)
       END DO
    END DO

    !
    ! Initialize output arrays with free atmosphere values
    !
    DO k=1,plev + 1
       DO i=1,plon
          kvm(i,k)   = kvf(i,k)
          kvh(i,k)   = kvf(i,k)
          cgh(i,k)   = 0.0_r8
          cgq(i,k,1) = 0.0_r8
          cgs(i,k)   = 0.0_r8
       END DO
    END DO
    DO m=2,pcnst
       DO k=1,plev + 1
          DO i=1,plon
             cgq(i,k,m) = 0.0_r8
          END DO
       END DO
    END DO
    !
    ! Compute various arrays for use later:
    !
    DO i=1,plon
       !thvsrf(i) = th(i,plev)*(1.0_r8 + 0.61_r8*q(i,plev,1))

       thvsrf(i) = (tsk(i)*psomc(i,plev))*(1.0_r8 + 0.61_r8*qsfc(i))
       !
       !
       !  ustar*tstar = ustar*tstar + ustar*qstar*(0.61*Th)
       !
       !
       !heatv(i)  = khfs(i) + 0.61_r8*th(i,plev)*kqfs(i,1)
       heatv(i)  = khfs(i) *(1.0_r8 + 0.61_r8*qsfc(i))!+ 0.61_r8*(tsk(i)*psomc(i,plev))*kqfs(i,1)
       wm(i)     = 0.0_r8
       therm(i)  = 0.0_r8
       qpert(i)  = 0.0_r8
       tpert(i)  = 0.0_r8
       fak3(i)   = 0.0_r8  
       zh(i)     = 0.0_r8  
       !
       !! Obukhov length
       !
       !
       !                 ustar^3
       ! L = ----------------------------------
       !       -                           -
       !      |         g          Ho       |
       !      | k * -------- * ------------ |
       !      |         To       rho * cp   |
       !       -                           -
       ! L > 0 indicates stable conditions
       ! L < 0 indicates unstable conditions
       ! L --> infinito  applies to neutral conditions
       ! 
       obklen(i) = -thvsrf(i)*ustar(i)**3/ &
            (g*vk*(heatv(i) + SIGN(1.e-10_r8,heatv(i))))
       
       IF(obklen(i) == 0.0_r8)obklen(i) = 1.0e-12_r8
    END DO
    !
    ! 
    ! >>>> Define first a new factor fac=100 for use in Richarson number
    !      Calculate virtual potential temperature first level
    !      and initialize pbl height to z1
    !
    fac = 1.00_r8
    FacLand=1.00_r8
    !
    DO i=1,plon
       thvref(i) = th(i,plev)*(1.0_r8 + 0.61_r8*q(i,plev,1))
       pblh(i)   = z(i,plev)
       check(i)  = .TRUE.
       !
       ! Initialization of lowest level Ri number 
       ! (neglected in initial Holtslag implementation)
       !
       !      vv2(i) = fac*ustar(i)**2
       !      vv2(i) = max(vv2(i),tiny)
       !      rino(i,plev) = g*(thvsrf(i) - thvref(i))*z(i,plev)/ &
       !                    (thvref(i)*vv2(i))
       rino(i,plev) = 0.0_r8
    END DO
    !
    ! PBL height calculation:
    ! Search for level of pbl. Scan upward until the Richardson number between
    ! the first level and the current level exceeds the "critical" value.
    !
    DO k=plev-1,plev-npbl+1,-1
!    DO k=plev-1,1,-1

       DO i=1,plon
          IF (check(i)) THEN
             ! y = a + b * x
             IF(FRLAND(i) >=0.5_r8)THEN
                vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2  + FacLand*ustar(i)**2
             ELSE
                vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2  + fac*ustar(i)**2
             END IF
             !
             !     PARAMETER (tiny=1.e-36_r8)
             !
             vvk = MAX(vvk,tiny)
             tkv = th(i,k)*(1.0_r8 + 0.61_r8*q(i,k,1))
             !
             !           (Tpv(i) - Tpvo) * ( Z (i) - Zo) 
             !  Ri = g ------------------------------------------------
             !          Tpvo * ( (DU(i)^2  +  DV(i))^2  + fac*ustar^2 )
             !
             rino(i,k) = g*(tkv - thvref(i))*(z(i,k)-z(i,plev)) &
                  /(thvref(i)*vvk)                  
             IF(rino(i,k).GE.ricr) THEN
                !               -                -
                !              |                  |
                !         Ricr | u(h)^2 +  v(h)^2 |
                !              |                  |
                !               -                -
                ! h = ----------------------------------------
                !                -    -           -                    -
                !               |   g  |   |                     |
                !               | ---- | * | Thv(h) +  Ths    |
                !               |  Ths |   |                     |
                !               -    -           -                    -
                !
                !
                !
                !                Ri(k+1)       Ric/Ri - 1
                ! h = Z(k+1)   +          -------------------- *(Z(k) - Z(k+1))
                !                Ri(k+1)      Ri(k)/Ri - 1
                !
                !

                !
                !
                !                   Ric - Ri(k+1)
                ! h = Z(k+1)   +  -------------------- *(Z(k) - Z(k+1))
                !                   Ri(k) - Ri(k+1)
                !
                !
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/ &
                     (rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
                check(i) = .FALSE.
             END IF
          END IF
       END DO
    END DO
    !
    ! Set pbl height to maximum value where computation exceeds number of
    ! layers allowed
    !
    DO i=1,plon
       IF (check(i)) pblh(i) = z(i,plev + 1-npbl)
    END DO
    !
    ! Improve estimate of pbl height for the unstable points.
    ! Find unstable points (virtual heat flux is positive):
    !
    DO i=1,plon
          !
          ! heatv = ustar*tstar 
          !
       IF (heatv(i) .GT. 0.0_r8) THEN
          unstbl(i) = .TRUE.
          check(i) = .TRUE.
       ELSE
          unstbl(i) = .FALSE.
          check(i) = .FALSE.
       END IF
    END DO
    !
    ! For the unstable case, compute velocity scale and the
    ! convective temperature excess:
    !
    DO i=1,plon
       IF (check(i)) THEN
          !
          ! sffrac = 0.1_r8     ! Surface layer fraction of boundary layer common /compbl/
          ! 
          ! betam  = 15.0_r8 ! Constant in wind gradient expression common /compbl/
          !
          ! binm = betam*sffrac
          ! onet = 1.0_r8/3.0_r8
          !
          !              -                                        -  1/3
          !             |                                         |
          !             |                         HPBL                 |  
          !  phiminv =  |1.0_r8 - 0.1 * 15 * ------------------- |
          !             |                           L                 |  
          !             |                                         |
          !              -                                        -
          !
          phiminv(i) = (1.0_r8 - binm*pblh(i)/obklen(i))**onet
          !IF(phiminv(i) > 10.0_r8)WRITE(nfprt,*)'1',phiminv(i), pblh(i),obklen(i)

          !
          !            -                                                      -  1/3
          !           |                                                       |
          !           |                                         HPBL               |  
          ! um =   |ustar^3 - ustar^3 * 0.1 * 15 * ------------------- |
          !           |                                           L           |  
          !           |                                                       |
          !            -                                                      -
          !
          wm(i) = ustar(i)*phiminv(i)
          therm(i) = heatv(i)*fak/wm(i)       
          !--         rino(i,plev) = -g*therm(i)*z(i,plev)/(thvref(i)*vv2(i))
          rino(i,plev) = 0.0_r8
          tlv(i) = thvref(i) + therm(i)
       END IF
    END DO
    !
    ! Improve pblh estimate for unstable conditions using the
    ! convective temperature excess:
    !
    DO k=plev-1,plev-npbl+1,-1
!    DO k=plev-1,1,-1
       DO i=1,plon
          IF (check(i)) THEN
             IF(FRLAND(i) >=0.5_r8)THEN
                vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2 + FacLand*ustar(i)**2
             ELSE
                vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2 + fac*ustar(i)**2
             END IF
             vvk = MAX(vvk,tiny)
             tkv = th(i,k)*(1.0_r8 + 0.61_r8*q(i,k,1))
             !
             !
             !             
             ! Ri = g * -----------------------
             !
             !
             !
             rino(i,k) = g*(tkv - tlv(i))*(z(i,k)-z(i,plev)) &
                  /(thvref(i)*vvk)
             IF(rino(i,k).GE.ricr) THEN
                !
                !
                !                   Ric - Ri(k+1)
                ! h = Z(k+1)   +  -------------------- *(Z(k) - Z(k+1))
                !                   Ri(k) - Ri(k+1)
                !
                !
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/ &
                     (rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
                check(i) = .FALSE.
             END IF
          END IF
       END DO
    END DO
    !
    ! Points for which pblh exceeds number of pbl layers allowed;
    ! set to maximum
    !
    !
    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*.  We 
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often 
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter.  Here we make use of the experimental results of Koracin 
    ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700.
    !
    DO i=1,plon
!       pblmin  = 700.0_r8*ustar(i)
       IF(FRLAND(i) >=0.5_r8)THEN
          pblmin  = 400.0_r8*ustar(i)                       ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'. 
       ELSE
          pblmin  = 400.0_r8*ustar(i)                       ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'. 
       END IF
!       pblmin  = 200.0_r8*ustar(i)                       ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'. 

       pblh(i) = MAX(MAX(pblh(i),pblmin),1.0_r8)
    END DO
    !
    ! pblh is now available; do preparation for diffusivity calculation:
    !
    DO i=1,plon
       pblk(i) = 0.0_r8
       fak1(i) = ustar(i)*pblh(i)*vk
       !
       ! Do additional preparation for unstable cases only, set temperature
       ! and moisture perturbations depending on stability.
       !
       IF (unstbl(i)) THEN
          !phiminv(i) = (1.0_r8 - binm*pblh(i)/obklen(i))**onet
          phiminv(i) = (1.0_r8 - binm*pblh(i)/obklen(i))**onet
          !IF(phiminv(i) > 10.00_r8)WRITE(nfprt,*)'2',phiminv(i), pblh(i),obklen(i)
          
          phihinv(i) = SQRT(1.0_r8 - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*g*pblh(i)/thvref(i))**onet 
          fak3(i)    = fakn*wstr(i)/wm(i)
          tpert(i)   = MAX(khfs(i)*fak/wm(i),0.0_r8)   
          qpert(i)   = MAX(kqfs(i,1)*fak/wm(i),0.0_r8)    
       ELSE
          tpert(i)   = MAX(khfs(i)*fak/ustar(i),0.0_r8) 
          qpert(i)   = MAX(kqfs(i,1)*fak/ustar(i),0.0_r8) 
       END IF
    END DO
    !
    ! temperature scale
    !
    DO i=1,plon
       !j/kg/kelvin
       !
       ! P = rho * R * T
       !
       !            P
       ! rho  = -------
       !          R * T
       !
       !           1        R * T
       !  rrho = ----- =   -------
       !          rho         P
       !
       !
       rrho(i) = gasr*t(i,plev)/pmid(i,plev)
       !
       !          _       _  2/3             _      _ -1/3
       !         |    Ho    |               | g*h    |
       ! Tstar = | -------- |           *   |--------|          =   K
       !         |   rho*cp |               |  To    |
       !          -       -                  -      -
       !
       tstar(i) = (((MAX(shflx(i),0.0001_r8)*rrho(i)/cp)**(2.0_r8))**(1.0_r8/3.0_r8))   * &
                    (1.0_r8/(((grav*MAX(pblh(i),0.1_r8))/t(i,plev))**(1.0_r8/3.0_r8)))

       !
       !
       !          _                  _ 1/3
       !         |    Ho        g*h   |
       ! wstar = | -------- * --------|   =   K
       !         |  rho*cp      To    |
       !          -                  -
       ! 
       wstar(i) =  ((MAX(shflx(i),0.0001_r8)*rrho(i)/cp) * ((grav*MAX(pblh(i),0.1_r8))/t(i,plev)) )**(1.0_r8/3.0_r8)

    END DO
    !
    ! Main level loop to compute the diffusivities and 
    ! counter-gradient terms:
    !
    ktopbl=1
    DO 1000 k=plev,plev-npbl+2,-1
    !DO 1000 k=plev,2,-1
       !
       ! Find levels within boundary layer:
       !
       DO i=1,plon
          unslev(i) = .FALSE.! unstbl pbl with levels within pbl
          stblev(i) = .FALSE.
          zm(i) = z(i,k)
          zp(i) = z(i,k-1)
          IF (zkmin.EQ.0.0_r8 .AND. zp(i).GT.pblh(i)) zp(i) = pblh(i)
          IF (zm(i) .LT. pblh(i)) THEN
             ktopbl(i) = k
             zmzp = 0.5_r8*(zm(i) + zp(i))
             !
             !        0.5 *( Z(k) - Z(k-1) )
             !  zh = ------------------------
             !                 h
             !
             !
             !             
             zh(i) = zmzp/pblh(i)
             !
             !        0.5 *( Z(k) - Z(k-1) )
             !  zl = ------------------------
             !                 L
             !
             !
             !
             !
             !
             !                 ustar^3
             ! L = ----------------------------------
             !       -                           -
             !      |         g          Ho       |
             !      | k * -------- * ------------ |
             !      |         To       rho * cp   |
             !       -                           -
             !

             zl(i) = zmzp/obklen(i)
             !
             zzh(i) = 0.0_r8
             !              -           - 2
             !             |          z  |
             !zzh =        |1   -   -----|
             !             |          h  |
             !              -           -
             !
             IF (zh(i).LE.1.0_r8) zzh(i) = (1.0_r8 - zh(i))**2
             !
             ! stblev for points zm < plbh and stable and neutral
             ! unslev for points zm < plbh and unstable
             !
             IF (unstbl(i)) THEN
                unslev(i) = .TRUE. ! unstbl pbl with levels within pbl
             ELSE
                stblev(i) = .TRUE.
             END IF
          END IF
       END DO
       !
       ! Stable and neutral points; set diffusivities; counter-gradient
       ! terms zero for stable case:
       !
       DO i=1,plon
          IF (stblev(i)) THEN
             IF (zl(i).LE.1.0_r8) THEN
                !
                !               -           - 2
                !              |          z  |
                !  Kc = k*wi*z*|1   -   -----|
                !              |          h  |
                !               -           -
                !
                !  fak1(i) = ustar(i)*pblh(i)*vk
                !
                pblk(i) = fak1(i)*zh(i)*zzh(i)/(1.0_r8 + betas*zl(i))
             ELSE
                !
                !               -           - 2
                !              |          z  |
                !  Kc = k*wi*z*|1   -   -----|
                !              |          h  |
                !               -           -
                !
                pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas + zl(i))
             END IF
             kvm(i,k) = MAX(pblk(i),kvf(i,k))
             kvh(i,k) = kvm(i,k)
          END IF
       END DO
       !
       ! unssrf, unstable within surface layer of pbl
       ! unsout, unstable within outer   layer of pbl
       !
       DO i=1,plon
          unssrf(i) = .FALSE.
          unsout(i) = .FALSE.
          IF (unslev(i)) THEN
             ! unstbl pbl with levels within pbl
             IF (zh(i).LT.sffrac) THEN
                unssrf(i) = .TRUE.
             ELSE
                unsout(i) = .TRUE.
             END IF
          END IF
       END DO
       !
       ! Unstable for surface layer; counter-gradient terms zero
       !
       DO i=1,plon
          IF (unssrf(i)) THEN
             !
             !        0.5 *( Z(k) - Z(k-1) )
             !  zl = ------------------------
             !                 L
             !
             !
             !          -                                       - 1/3
             !         |                 0.5 *( Z(k) - Z(k-1) )  |
             !  term = |1   -  betam *  ------------------------ |
             !         |                          L              |
             !          -                                        -
             !
             term = (1.0_r8 - betam*zl(i))**onet
             !
             !               -           - 2
             !              |               z  |
             !  Kc = k*wi*z*|1   -   -----|
             !              |               h  |
             !               -           -
             !
             ! fak1(i) = ustar(i)*pblh(i)*vk
             !
             !                                                         -           - 2   -                                       - 1/3
             !                               0.5 *( Z(k) - Z(k-1) )   |          z  |   |                  0.5 *( Z(k) - Z(k-1) ) |
             !  Kc =  ustar(i)* h(i) * k *  ------------------------ *|1   -   -----| * |1   -  betam *  ------------------------ |
             !                                        h               |          h  |   |                           L             |
             !                                                         -           -     -                                       -
             !
             !                                               -           - 2   -                                       - 1/3
             !                                              |          z  |   |                  0.5 *( Z(k) - Z(k-1) )  |
             !  Kc =  ustar(i)* k * 0.5 *( Z(k) - Z(k-1) ) *|1   -   -----| * |1   -  betam *  ------------------------  |
             !                                              |          h  |   |                        L                 |
             !                                               -           -     -                                        -

             pblk(i) = fak1(i)*zh(i)*zzh(i)*term

             pr(i) = term/SQRT(1.0_r8 - betah*zl(i))
          END IF
       END DO
       !
       ! Unstable for outer layer; counter-gradient terms non-zero:
       !
       DO i=1,plon
          IF (unsout(i)) THEN
             !        -           - 2
             !       |          z  |
             !zzh =  |1   -   -----|
             !       |          h  |
             !        -           -
             !
             !
             !        0.5 *( Z(k) - Z(k-1) )
             !  zh = ------------------------
             !                 h
             !
             !            ____
             !  DCo     d w*Co
             ! ----- = -------
             !  Dt      dz
             !
             !               -         - 
             !  DCo     d   |      dCo  |   
             ! ----- = ---* |Kc * ----- | 
             !  Dt      dz  |       dz  |   
             !               -         - 
             ! fakn   = 7.2_r8     ! Constant in turbulent prandtl number

             !        ____
             ! fak3 = w*Co = fakn*wstr(i)/wm(i)
             !
             ! fak2(i)    = wm(i)*h(i)*vk
             !
             !                                                  -               - 2
             !                        0.5 *( Z(k) - Z(k-1) )   |              z  |
             ! Kc = wm  * h  * vk *  ------------------------ *|1       -   -----|
             !                                   h             |              h  |
             !                                                  -               -
             !
             !                                               -                - 2
             !                                              |               z  |
             ! Kc = wm  *  vk *  0.5 *( Z(k) - Z(k-1) )   * |1        -   -----|
             !                                              |               h  |
             !                                               -                -
             !
             !               -           - 2
             !              |          z  |
             !  Kc = k*wi*z*|1   -   -----|
             !              |          h  |
             !               -           -
             !
             pblk(i) = fak2(i)*zh(i)*zzh(i)
             !         ____
             !         w*Co
             ! Yc = d*---------
             !         Wstar*h
             !
             cgs (i,k) = fak3(i)/(pblh(i)*wm(i))
             !
             !  khfs = surface kinematic heat flux [mK/s]
             !
             !
             ! 
             ! surface kinematic heat flux [mK/s] ustar*tstar
             !
             ! H=Ho=-rho*cp*ustar*tstar  
             !
             ! surface heat flux (W/m2)  (j/kg/kelvin)
             !
             !                  H           j       kg K   m*m*m           m K
             ! ustar*tstar = ---------   = ------------------------- =   ------------
             !                -rho*cp      m*m*s   j         kg              s
             !
             ! khfs(i) = shflx(i)*rrho(i)/cpair
             !
             !
             !cgh (i,k)  = counter-gradient term for heat [K/m]
             !
             cgh (i,k) = khfs(i)*cgs(i,k)
             !
             !       Qh             z     Wstar
             ! Pr = ---- + a * k * --- * -------
             !       Qm             h     Wm 
             !
             !
             !       phiminv              fak3 
             ! Pr = ----------  + ccon * ----- 
             !       phihinv              fak  
             !
             !    ccon = fak*sffrac*vk

             ! fakn   = 7.2_r8     ! Constant in turbulent prandtl number

             !        ____
             ! fak3 = w*Co = fakn*wstr(i)/wm(i)
             !
             !       phiminv                        fakn*wstr(i)
             ! Pr = ----------  + fak*sffrac*vk * -------------------- 
             !       phihinv                         fak * wm(i)
             !
             !
             !       phiminv                         wstr(i)
             ! Pr = ----------  + fakn*vk*sffrac * ----------------- 
             !       phihinv                         wm(i)
             !

             pr  (i) = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
             
             !
             ! kqfs = sfc kinematic constituent flux [m/s]
             !
             !
             ! surface kinematic constituent flux [m/s]
             !
             ! E=Eo=-rho*cp*ustar*qstar
             !
             !                  E            kg          m * m * m          m 
             ! ustar*qstar = --------- = ----------- * ------------- = -----------
             !                -rho        m * m * s       kg                s
             !
             ! 
             ! surface constituent flux (kg/m2/s)
             !
             ! kqfs(i,1) = cflx(i,1)*rrho(i)
             !
             ! gq (i,k,1) = counter-gradient term for constituents
             !
             cgq (i,k,1) = kqfs(i,1)*cgs(i,k)
          END IF
       END DO
       DO m=2,pcnst
          DO i=1,plon
             IF (unsout(i)) cgq(i,k,m) = kqfs(i,m)*cgs(i,k)
          END DO
       END DO
       !
       ! For all unstable layers, set diffusivities
       !
       DO i=1,plon
          IF (unslev(i)) THEN
             !
             ! unstbl pbl with levels within pbl
             !
             kvm(i,k) = MAX(pblk(i)      ,kvf(i,k))
             kvh(i,k) = MAX(pblk(i)/pr(i),kvf(i,k))
          END IF
       END DO
1000   CONTINUE           ! end of level loop
       !+
       ! Check whether last allowed midpoint is within pbl, determine ktopblmn
       !-
       ktopblmn = plev
       k = plev-npbl+1
       DO i = 1, plon
          IF (z(i,k) .LT. pblh(i)) ktopbl(i) = k
          ktopblmn(i) = MIN(ktopblmn(i), ktopbl(i))
       END DO
       DO i=1,plon
          ktopblmn(i)=MAX(ktopblmn(i),2)
       END DO
       !
       !
       ! Compute various arrays for use later:
       !
       TKE=0.0_r8
       DO  k=plev,1,-1
          DO i=1,plon
             IF(z(i,k) <= pblh(i))THEN
                IF(ABS(obklen(i)) > pblh(i))obklen(i) = (obklen(i)/ABS(obklen(i)))*pblh(i)
                IF(z(i,k)/obklen(i) >= 0.0_r8)THEN
                   IF(z(i,k) <= 0.1_r8*pblh(i))THEN
                      !  In the surface layer, the TKE e and EDR   are given by Hogstrom, 1996; Rao and
                      !  Nappo, 1998
                      tke(i,k+1) = 6*(ustar(i)**2)
                   ELSE IF(z(i,k) > 0.1_r8*pblh(i))THEN 
                      !   In the neutral and moderately stable boundary layer, the TKE and EDR are given by
                      ! Hogstrom, 1996; Rao and Nappo, 1998
                      tke(i,k+1) = (6.0_r8*(ustar(i)**2))*(1.0_r8 - (z(i,k)/pblh(i)))**1.75
                   ELSE
                      tke(i,k+1) = 6*(ustar(i)**2)
                   END IF
                ELSE IF(z(i,k)/obklen(i) < 0.0_r8)THEN
                   IF(z(i,k) <= 0.1_r8*pblh(i))THEN
                      !In the surface layer, the TKE and EDR are given by Arya, 2000
                      tke(i,k+1) = (0.36_r8*(wstar(i)**2))   + &
                         (0.85_r8*(ustar(i)**2)*(1.0_r8 - 3.0_r8*(z(i,k)/pblh(i)))**2/3)
                   ELSE IF(z(i,k) > 0.1_r8*pblh(i))THEN 
                      !In the mixed layer, the TKE is given by Arya, 2000
                      tke(i,k+1) = (0.36_r8 + 0.9_r8*((z(i,k)/pblh(i))**2/3) * &
                             (1.0_r8 - 0.8_r8*((z(i,k)/pblh(i))))**2)*(wstar(i)**2)
                   ELSE
                      tke(i,k+1) = 6*(wstar(i)**2)
                   END IF
                END IF  
             END IF
          END DO
       END DO
     !




       RETURN
     END SUBROUTINE pbldif





     SUBROUTINE qvdiff(ncnst   , &! INTENT(IN   ) :: ncnst                  ! number of constituents being diffused
                       qm1     , &! INTENT(IN   ) :: qm1   (plond,plev,ncnst)  ! initial constituent
                       qflx    , &! INTENT(IN   ) :: qflx  (plond,ncnst)      ! sfc q flux into lowest model level
                       cc      , &! INTENT(IN   ) :: cc    (plond,plev)         ! -lower diag coeff.of tri-diag matrix
                       ze      , &! INTENT(IN   ) :: ze    (plond,plev)         ! term in tri-diag. matrix system
                       term    , &! INTENT(IN   ) :: term  (plond,plev)       ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
                       qp1     , &! INTENT(OUT  ) :: qp1   (plond,plev,ncnst)  ! final constituent
                       plon    , &! INTENT(IN   ) :: plon   ! number of longitudes
                       plond   , &! INTENT(IN   ) :: plond  ! slt extended domain longitude
                       plev    , &! INTENT(IN   ) :: plev   ! number of vertical levels
                       pcnst   , &! INTENT(IN   ) :: pcnst  ! number of constituents (including water vapor)
                       ntopfl   ) ! INTENT(IN   ) :: ntopfl ! Top level to which vertical diffusion is applied.
       !-----------------------------------------------------------------------
       !
       ! Solve vertical diffusion eqtn for constituent with explicit srfc flux.
       ! Procedure for solution of the implicit equation follows Richtmyer and 
       ! Morton (1967,pp 198-199).
       !
       !---------------------------Code history--------------------------------
       !
       ! Original version:  CCM1
       ! Standardized:      J. Rosinski, June 1992
       ! Reviewed:          P. Rasch, B. Boville, August 1992
       ! Reviewed:          P. Rasch, April 1996
       ! Reviewed:          B. Boville, May 1996
       !
       !-----------------------------------------------------------------------
       !
       ! $Id: qvdiff.F,v 1.1.1.1 2001/03/09 00:29:27 mirin Exp $
       !
       !-----------------------------------------------------------------------
       !-----------------------------------------------------------------------
       !------------------------------Commons----------------------------------
       !------------------------------Arguments--------------------------------
       !
       ! Input arguments
       !
       !      integer,parameter  :: plon   = 1 ! number of longitudes
       !      integer,parameter  :: plond  = 1 ! slt extended domain longitude
       !      integer,parameter  :: plev   =28 ! number of vertical levels
       !      integer,parameter  :: pcnst  = 1 ! number of constituents (including water vapor)
       !      integer,parameter  :: ntopfl =1 ! Top level to which vertical diffusion is applied.
       ! ntopfl = top level to which v-diff is applied
       IMPLICIT NONE
       INTEGER, INTENT(IN   )  :: plon   ! number of longitudes
       INTEGER, INTENT(IN   )  :: plond  ! slt extended domain longitude
       INTEGER, INTENT(IN   )  :: plev   ! number of vertical levels
       INTEGER, INTENT(IN   )  :: pcnst  ! number of constituents (including water vapor)
       INTEGER, INTENT(IN   )  :: ntopfl ! Top level to which vertical diffusion is applied.
       ! ntopfl = top level to which v-diff is applied

       INTEGER, INTENT(IN   )  :: ncnst                  ! number of constituents being diffused

       REAL(KIND=r8), INTENT(IN   ) :: qm1   (plond,plev,ncnst)  ! initial constituent
       REAL(KIND=r8), INTENT(IN   ) :: qflx  (plond,ncnst)      ! sfc q flux into lowest model level
       REAL(KIND=r8), INTENT(IN   ) :: cc    (plond,plev)         ! -lower diag coeff.of tri-diag matrix
       REAL(KIND=r8), INTENT(INOUT) :: ze    (plond,plev)         ! term in tri-diag. matrix system
       REAL(KIND=r8), INTENT(IN   ) :: term  (plond,plev)       ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
       !
       ! Output arguments
       !
       REAL(KIND=r8), INTENT(OUT  ) :: qp1   (plond,plev,ncnst)  ! final constituent
       !
       !---------------------------Local workspace-----------------------------
       !
       REAL(KIND=r8) :: zfq   (plond,plev,pcnst)  ! terms appear in soln of tri-diag sys
       REAL(KIND=r8) :: tmp1d (plond)           ! temporary workspace (1d array)

       INTEGER i,k                 ! longitude, vertical indices
       INTEGER m                   ! constituent index


       qp1=0.0_r8
       zfq=0.0_r8   
       tmp1d =0.0_r8
       !
       !-----------------------------------------------------------------------
       !
       ! Calculate fq(k).  Terms fq(k) and e(k) are required in solution of 
       ! tridiagonal matrix defined by implicit diffusion eqn.
       ! Note that only levels ntopfl through plev need be solved for.
       ! No vertical diffusion is applied above this level
       !
       DO m=1,ncnst
          DO i=1,plon
             zfq(i,ntopfl,m) = qm1(i,ntopfl,m)*term(i,ntopfl)
          END DO
          DO k=ntopfl+1,plev-1
             DO i=1,plon
                zfq(i,k,m) = (qm1(i,k,m) + cc(i,k)*zfq(i,k-1,m))*term(i,k)
             END DO
          END DO
       END DO
       !
       ! Bottom level: (includes  surface fluxes)
       !
       DO i=1,plon
          tmp1d(i) = 1.0_r8/(1.0_r8 + cc(i,plev) - cc(i,plev)*ze(i,plev-1))
          ze(i,plev) = 0.0_r8
       END DO
       DO m=1,ncnst
          DO i=1,plon
             zfq(i,plev,m) = (qm1(i,plev,m) + qflx(i,m) + &
                  cc(i,plev)*zfq(i,plev-1,m))*tmp1d(i)
          END DO
       END DO
       !
       ! Perform back substitution
       !
       DO m=1,ncnst
          DO i=1,plon
             qp1(i,plev,m) = zfq(i,plev,m)
          END DO
          DO k=plev-1,ntopfl,-1
             DO i=1,plon
                qp1(i,k,m) = zfq(i,k,m) + ze(i,k)*qp1(i,k+1,m)
             END DO
          END DO
       END DO
       RETURN
     END SUBROUTINE qvdiff
 
     SUBROUTINE MY20_TURB (       &
                     pcnst      ,& ! INTENT(IN   ) :: pcnst
                     plond      ,& ! INTENT(IN   ) :: plond
                     plev       ,& ! INTENT(IN   ) :: plev
                     ztodt      ,& ! INTENT(IN   ) :: ztodt
                     colrad     ,& ! INTENT(IN   ) :: colrad        (plond)    
                     prsi       ,& !(1:plond,1:plev+1) ,& ! phalf    (1:plond,1:plev+1)!!          phalf    -  Pressure at half levels, 
                     prsl       ,& !(1:plond,1:plev)   ,& ! pfull    (1:plond,1:plev)!!          pfull    -  Pressure at full levels! , zhalf, zfull
                     tm1        ,& ! INTENT(IN   ) :: tm1        (plond,plev)            ! temperature input
                     qp1        ,& ! INTENT(IN   ) :: qp1        (plond,plev,pcnst)  ! constituents after vdiff
                     um1        ,& ! INTENT(IN   ) :: um1        (plond,plev)            ! u-wind input
                     vm1        ,& ! INTENT(IN   ) :: vm1        (plond,plev)            ! v-wind input
                     gl0        ,& ! INTENT(INOUT) :: gl0        (plond)     !
                     pblh       ,& ! INTENT(IN   ) :: pblh        (plond)    
                     tkemyj     ,& ! INTENT(IN   ) :: tkemyj  (plond,plev)      
                     KmMixl     ,& ! INTENT(OUT  ) :: KmMixl(plond,plev)   
                     KhMixl      ) ! INTENT(OUT  ) :: KhMixl(plond,plev)
      IMPLICIT NONE
      INTEGER       , INTENT(IN   ) :: pcnst
      INTEGER       , INTENT(IN   ) :: plond
      INTEGER       , INTENT(IN   ) :: plev
      REAL(kind=r8), INTENT(IN   ) :: ztodt
      REAL(kind=r8), INTENT(IN   ) :: colrad    (plond)    
      REAL(kind=r8), INTENT(IN   ) :: prsi      (plond,plev+1)  ! phalf    (1:plond,1:plev+1)!!          phalf    -  Pressure at half levels, 
      REAL(kind=r8), INTENT(IN   ) :: prsl      (plond,plev)    ! pfull    (1:plond,1:plev)!!          pfull    -  Pressure at full levels! , zhalf, zfull
      REAL(kind=r8), INTENT(IN   ) :: tm1       (plond,plev)        ! temperature input
      REAL(kind=r8), INTENT(IN   ) :: qp1       (plond,plev,pcnst)  ! constituents after vdiff
      REAL(kind=r8), INTENT(IN   ) :: um1       (plond,plev)        ! u-wind input
      REAL(kind=r8), INTENT(IN   ) :: vm1       (plond,plev)        ! v-wind input
      REAL(kind=r8), INTENT(INOUT) :: gl0       (plond)             !
      REAL(KIND=r8), INTENT(inOUT) :: pblh      (plond)  
      REAL(KIND=r8), INTENT(OUT  ) :: tkemyj    (plond,plev+1)  
      REAL(KIND=r8), INTENT(OUT  ) :: KmMixl    (plond,plev) 
      REAL(KIND=r8), INTENT(OUT  ) :: KhMixl    (plond,plev)
  
      REAL(kind=r8) :: twodt 
      REAL(kind=r8) :: twodti
      REAL(KIND=r8) :: a   (plond,plev)
      REAL(KIND=r8) :: b   (plond,plev)
      REAL(KIND=r8) :: x,y
      REAL(KIND=r8) :: Pbl_HgtLyI(plond,plev)
      REAL(KIND=r8) :: Pbl_ATemp (plond,plev) 
      REAL(KIND=r8) :: Pbl_Stabil(plond,plev)
      REAL(KIND=r8) :: Pbl_ITemp (plond,plev)
      REAL(KIND=r8) :: Pbl_Shear (plond,plev)
      REAL(KIND=r8) :: Pbl_KM    (plond,plev) 
      REAL(KIND=r8) :: Pbl_KH    (plond,plev) 
      REAL(KIND=r8) :: Pbl_BRich (plond,plev)
      REAL(KIND=r8) :: Pbl_NRich (plond,plev)
      REAL(KIND=r8) :: Pbl_ShBar (plond,plev)
      REAL(KIND=r8) :: Pbl_SmBar (plond,plev) 
      REAL(KIND=r8) :: Pbl_Sqrtw (plond,plev)
      REAL(KIND=r8) :: Pbl_EddEner(plond,plev) 
      REAL(KIND=r8) :: Pbl_MixLgh (plond,plev)
      REAL(KIND=r8) :: MixLgh     (plond,plev+1) 
      REAL(KIND=r8) :: Pbl_PotTep (plond,plev)
 
      REAL(KIND=r8) :: Pbl_KmMixl(plond,plev)
      REAL(KIND=r8) :: Pbl_KhMixl(plond,plev)
      REAL(KIND=r8) :: gln(plond)
      REAL(KIND=r8) :: gld(plond)
      REAL(KIND=r8) :: csqiv(plond)
      REAL(KIND=r8) :: s1ms2g 
      REAL(KIND=r8) :: fac
      REAL(KIND=r8) :: PK_a0    (plond,plev)
      REAL(KIND=r8) :: PK_b0    (plond,plev)
      REAL(KIND=r8) :: PK_t0    (plond,plev)
      REAL(KIND=r8) :: PK_t1    (plond,plev)
      REAL(KIND=r8) :: PK_sigr  (plond,plev)
      REAL(KIND=r8) :: PK_sigriv(plond,plev)
      REAL(KIND=r8) :: PK_con0  (plond,plev)
      REAL(KIND=r8) :: PK_con1  (plond,plev)
      REAL(KIND=r8) :: PK_con2  (plond,plev)
      REAL(KIND=r8) ::  phalf    (1:plond,1:plev+1)!!	  phalf    -  Pressure at half levels, 
       REAL(KIND=r8) :: pfull    (1:plond,1:plev  )!!  	pfull	 -  Pressure at full levels! , zhalf, zfull
      INTEGER       :: i
      INTEGER       :: k,l
      INTEGER       :: itr
      INTEGER       :: nsurf
      KmMixl=0.0_r8
      KhMixl=0.0_r8 

       DO k = 1, plev+1
          DO i=1,plond
             phalf(i,plev+1-k+1) = prsi(i,k)
         END DO
      END DO
       DO k = 1, plev
          DO i=1,plond
             pfull(i,plev-k+1) = prsl(i,k)
         END DO
      END DO

      DO i=1,plond,1 
          csqiv(i)        = 1.0e0_r8/SIN(colrad(i))**2
      END DO

      nsurf=plev
      twodt =ztodt
      twodti=1.0_r8/twodt
 
      DO k = 1, plev
          DO i=1,plond
           ! con0  (k)=gasr*delsig(k)/(grav*sig(k)) 
           PK_con0  (i,k)=gasr*((phalf(i,k)/phalf(i,1)) - (phalf(i,k+1)/phalf(i,1)))/(grav*(pfull(i,k)/phalf(i,1)))
         END DO
      END DO 
 
      DO i=1,plond
         PK_a0    (i,plev)=0.0_r8
         PK_t0    (i,plev)=0.0_r8
         PK_t1    (i,plev)=0.0_r8
         PK_b0    (i,   1)=0.0_r8
         PK_sigr  (i,plev)=0.0_r8
         PK_sigriv(i,   1)=0.0_r8
      END DO

       DO k = 1, plev-1
          DO i=1,plond
           ! con0  (k)=gasr*delsig(k)/(grav*sig(k)) 

           !PK_con0  (i,k)=gasr*((phalf(i,k)/phalf(i,1)) - (phalf(i,k+1)/phalf(i,1)))/(grav*(pfull(i,k)/phalf(i,1)))

           !con1  (   k)=grav*sigml(k+1)/(gasr*(sig(k)-sig(k+1)))

           PK_con1  (i,k) = grav*(phalf(i,k+1)/phalf(i,1))/(gasr*(((pfull(i,k)/phalf(i,1)) - (pfull(i,k+1)/phalf(i,1)))))

           !con2  (   k)=grav*con1(k)
           PK_con2  (i,k)=grav*PK_con1(i,k)

           !con1  (   k)=con1(k)*con1(k)
           PK_con1  (i,k)=PK_con1(i,k)*PK_con1(i,k)

!          t0    (   k)=(sig(k+1)-sigml(k+1))/(sig(k+1)-sig(k))

           PK_t0    (i,   k)=((pfull(i,k+1)/phalf(i,1)) - (phalf  (i,k+1)/phalf  (i,1))) / &
                             ((pfull(i,k+1)/phalf(i,1)) - (pfull(i,k)/phalf(i,1)))


 !          t1    (   k)=(sigml(k+1)-sig(k  ))/(sig(k+1)-sig(k))
           PK_t1    (i,   k)=((phalf  (i,k+1)/phalf  (i,1)) - (pfull(i,k)/phalf(i,1)))/&
                             ((pfull  (i,k+1)/phalf  (i,1)) - (pfull(i,k)/phalf(i,1)))

          !sigr  (   k)=sigk(k)*sigkiv(k+1)
           PK_sigr  (i,k)= ((pfull(i,k)/phalf(i,1))**akappa)*  (1.0_r8/((pfull(i,k+1)/phalf(i,1))**akappa))
          !sigriv( k+1)=sigk(k+1)*sigkiv(k)
           PK_sigriv(i,k+1)= ((pfull(i,k+1)/phalf(i,1))**akappa )  *  (1.0_r8/((pfull(i,k)/phalf(i,1))**akappa))

          !a0    (k   )=gbyr*sigml(k+1)**2/(delsig(k  )*(sig(k)-sig(k+1)))
          !REAL(KIND=r8),    INTENT(IN   ) :: phalf  (nCols,kMax+1)  !     phalf     - real, pressure at layer interfaces             ix,levs+1  Pa
          !REAL(KIND=r8),    INTENT(IN   ) :: pfull  (nCols,kMax)    !     pfull     - real, mean layer presure                       ix,levs   Pa

          PK_a0    (i,k   )=gbyr*((phalf  (i,k+1)/phalf  (i,1))**2)/(((phalf(i,k)/phalf(i,1)) - (phalf(i,k+1)/phalf(i,1)))*((pfull(i,k)/phalf(i,1)) - (pfull(i,k+1)/phalf(i,1))))

          !bb0    (i,k+1 )=gbyr*sigml(k+1)**2/(delsig(k+1)*(sig(k)-sig(k+1)))

          PK_b0    (i,k+1 )=gbyr*((phalf  (i,k+1)/phalf  (i,1))**2)/(((phalf(i,k+1)/phalf(i,1)) - (phalf(i,k+2)/phalf(i,1)))*((pfull(i,k)/phalf(i,1)) - (pfull(i,k+1)/phalf(i,1))))

          END DO
       END DO

      DO k = 1, plev
          DO i=1,plond

          !                 --                  --
          !    1      1    | g       si(k+1)      |     con1(k)
          !  ------ =--- * |--- * ----------------| = ----------
          !    DZ     T    | R     sl(k) -sl(k+1) |       T
          !                 --                  --
          !                -- --
          !    DA      d  |     |
          !  ------ =---- | W'A'|
          !    DT      dZ |     |
          !                -- --
          !
          !                          ----
          !    DA      d         d  |    |
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
          !         --  -- 2                             -- -- 2
          !        |   g  |       si(k+1)**2            |  1  | 
          !a0(k)  =| -----| * --------------------    = | --- |
          !        |   R  |     ((si(k)-si(k+1))  )     |  dZ | 
          !         --  --                               -- --
          !                   --   --      -- -- 2            --                   -- 2
          !                  |       |    |  1  |            |      m       kg * K   |  
          a(i,k)=twodt*PK_a0(i,k)!  |  2*Dt |  * | --- |    ==> s * |   ------- * --------  | 
          !                  |       |    |  dZ |            |     s**2       J      |  
          !                   --   --      -- --              --                   --
          !J = F*DX = kg m/s**2 *m = kg * m**2/s**2
          !                   --   --      -- -- 2            --                       -- 2
          !                  |       |    |  1  |            |    m         kg * K *s**2 |    K**2 * s
          !a(k)=twodt*a0(k)! |  2*Dt |  * | --- |    ==> s * | ------- * --------------- | = -----------
          !                  |       |    |  dZ |            |   s**2     kg * m**2      |    m**2  
          !                   --   --      -- --              --                       --
          !                   --   --      -- -- 2  
          !                  |       |    |  1  |         K**2 * s
          !a(k)=twodt*a0(k)! |  2*Dt |  * | --- |    ==> -----------
          !                  |       |    |  dZ |           m**2  
          !                   --   --      -- --    

          !
          !                   gbyr*sigml(k+1)**2
          !    b0(k) = -----------------------------------------------
          !               (((si(k+1)-si(k+1+1))  ) *(sig(k)-sig(k+1)))

          b(i,k)=twodt*PK_b0(i,k)! s * K**2/m**2 
         END DO
      END DO
      !c=twodt*c0
       !
       !                          1
       !      sigkiv(k   )=--------------------
       !                      (sl(k)**akappa)
       !
       !                       --                       --
       !                      | grav        si(k+1)       |
       !   con2  (   k)=grav* |------- * ---------------- |
       !                      | gasr      sl(k)-sl(k+1))  |
       !                       --                       --
       !
       !                --                       --  2
       !               |   grav        si(k+1)     |
       !   con1  (k) = |------- * ---------------- |
       !               |  gasr      sl(k)-sl(k+1)) |
       !                --                       --
       !
       ! sl.........sigma value at midpoint of                  gasr/cp
       !                                         each layer : (k=287/1005)=R/cp
       !
       !                                                                     1
       !                                             +-                   + ---
       !                                             !     k+1         k+1!  k
       !                                             !si(l)   - si(l+1)   !
       !                                     sl(l) = !--------------------!
       !                                             !(k+1) (si(l)-si(l+1)!
       !  --          --      --           --
       ! |  P(k)-P(k+1) |    |               |
       ! |--------------| =  |sl(k)-sl(k+1)  |
       ! |      P0      |    |               |
       !  --          --      --           --


       ! THERMODYNAMIC TEMPERATURE (K) at the top of of k-th layer (interface)
       !
       ! Do a linear interpolation in pressure to find the value
       ! between each two layers. In this way, interface k is
       ! above layer k, and below layer k+1
       !
       !                      sl(k+1)-si(k+1)
       !       t0    (   k)=-------------------
       !                      sl(k+1)-sl(k  )
       !
       !                      si(k+1)-sl(k  )
       !       t1    (   k)=-------------------
       !                      sl(k+1)-sl(k  )
       !
       !
       !                     sl(k+1)-si(k+1)              si(k+1)-sl(k  )
       !  Pbl_ATemp(i,k) = -------------------*gt(i,k) + ------------------*gt(i,k+1)
       !                     sl(k+1)-sl(k  )              sl(k+1)-sl(k  )
       !
       DO k = 1, plev-1
          DO i = 1, plond
             Pbl_ATemp(i,k)  = PK_t0(i,k)*tm1 (i,(plev+1-k)) + PK_t1(i,k)*tm1(i,(plev+1-k)-1)
          END DO
       END DO
 
       ! VIRTUAL POTENTIAL TEMPERATURE (K) 
       !   Pbl_PotTep (k) => in the middle of k-th layer
       !   Pbl_PotTintf(k) => in the inteface of k and k+1 layers
       !
       ! gt(k) thermodynamic temperature (K) in the middle of k-th layer
       ! gq(k) specific humidity (kg/kg) in the middle of k-th layer
       ! sigkiv(k) poisson factor (unitless) in the middle of k-th layer
       !     = (P(k)/P0)^-Rm/Cpm = sl(k)^-Rm/Cpm =~ sl(k)^Rd/Cp,d
       !
       !     A better approximation would be 
       !
       !        Rm/Cp,m = Rd/Cp,d * (1-0.251*q)
       !
       ! virtual temperature   = T * (1+eps*q)
       ! potential temperature = T * (P/P0)^-Rm/Cpm 
       !
       !
       !                      si(k+1)-sl(k  )
       !       t1    (   k)=-------------------
       !                      sl(k+1)-sl(k  )
       !
       !                          1
       !      sigkiv(k   )=--------------------
       !                      (sl(k)**akappa)
       !
       !                       --                       --  
       !                      | grav        si(k+1)       |
       !   con2  (   k)=grav* |------- * ---------------- |  
       !                      | gasr      sl(k)-sl(k+1))  |
       !                       --                       -- 
       !
       !                --                       --  2
       !               |   grav        si(k+1)     |
       !   con1  (k) = |------- * ---------------- |         
       !               |  gasr      sl(k)-sl(k+1)) |
       !                --                       -- 
       !
       ! sl.........sigma value at midpoint of                  gasr/cp
       !                                         each layer : (k=287/1005)=R/cp
       !
       !                                                                     1
       !                                             +-                   + ---
       !                                             !     k+1         k+1!  k
       !                                             !si(l)   - si(l+1)   !
       !                                     sl(l) = !--------------------!
       !                                             !(k+1) (si(l)-si(l+1)!
       !             --   --  -(R/Cp)                +-                  -+     
       !            |  P    |
       !Thetav = Tv |-------|
       !            |  P0   | 
       !             --   --
       !
       !             --   --  -(R/Cp)        --   --  -(R/Cp)
       !            |  P    |               |       |
       !sigkiv(k) = |-------|            == |sl(k)  |
       !            |  P0   |               |       | 
       !             --   --                 --   --
       !  --   --      --   -- 
       ! |  P    |    |       |
       ! |-------| =  |sl(k)  |
       ! |  P0   |    |       |
       !  --   --      --   --
       !  --          --      --           -- 
       ! |  P(k)-P(k+1) |    |               |
       ! |--------------| =  |sl(k)-sl(k+1)  |
       ! |      P0      |    |               |
       !  --          --      --           --




       DO k = 1, plev
          DO i = 1, plond
             !Pbl_PotTep(i,k)=(1.0_r8/((prsl(i,(plev+1-k))/prsi(i,nsurf+1))**akappa))*tm1 (i,(plev+1-k))*(1.0_r8+eps*qp1(i,(plev+1-k),1))
             Pbl_PotTep(i,k)=(1.0_r8/((pfull(i,k)/phalf(i,1))**akappa))*tm1 (i,(plev+1-k))*(1.0_r8+eps*qp1(i,(plev+1-k),1))
          END DO
       END DO
       !     
       !     Pbl_Stabil(2)   stability
       !     Pbl_Shear (6)   square of vertical wind shear
       !  
       gln =1.0e-5_r8 
       DO k = 1, plev-1
          DO i = 1, plond
             !
             !
             !                 g         D (Theta)
             !Pbl_Stabil(2) =------- * -------------
             !                Theta      D (Z)
             !
             !
             ! P =rho*R*T and P = rho*g*Z
             !
             !                           P
             ! DP = rho*g*DZ and rho = ----
             !                          R*T
             !        P
             ! DP = ----*g*DZ
             !       R*T
             !
             !        R*T
             ! DZ = ------*DP
             !        g*P
             !
             !    1       g*P       1   
             !  ------ = ------ * ------
             !    DZ      R*T       DP  
             !
             !    T       g         P   
             !  ------ = ------ * ------
             !    DZ      R         DP  
             !
             !  --   --      --   -- 
             ! |  P    |    |       |
             ! |-------| =  |sl(k)  |
             ! |  P0   |    |       |
             !  --   --      --   --
             !
             !    T       g            si(k+1)   
             !  ------ = ------ * ----------------- =con1(k)
             !    DZ      R         sl(k) -sl(k+1)  
             !
             !                 --                  --
             !    1      1    | g       si(k+1)      |     con1(k)
             !  ------ =--- * |--- * ----------------| = ----------
             !    DZ     T    | R     sl(k) -sl(k+1) |       T
             !                 --                  --
             !                                  --                  --
             !                    g       1    | g       si(k+1)      |   
             ! Pbl_Stabil(2) = ------- * --- * |--- * ----------------|* D(Theta)
             !                  Theta     T    | R     sl(k) -sl(k+1) | 
             !                                  --                  --
             !
             !                   con2(k)*(Pbl_PotTep(i,k+1)-Pbl_PotTep(i,k))
             !Pbl_Stabil(i,k)=-----------------------------------------------------------------
             !                  ((t0(k)*Pbl_PotTep(i,k)+t1(k)*Pbl_PotTep(i,k+1))*Pbl_ATemp(i,k))


             Pbl_Stabil(i,k)=PK_con2(i,k)*(Pbl_PotTep(i,k+1)-Pbl_PotTep(i,k)) &
                  /((PK_t0(i,k)*Pbl_PotTep(i,k)+PK_t1(i,k)*Pbl_PotTep(i,k+1))*Pbl_ATemp(i,k))
                  
             Pbl_ITemp(i,k)=1.0_r8/(Pbl_ATemp(i,k)*Pbl_ATemp(i,k))
             ! 
             !                  --  -- 2    --  -- 2
             !                 |  dU  |    |  dV  |
             !Pbl_Shear(i,k) = | -----|  + | -----| 
             !                 |  DZ  |    |  DZ  |
             !                  --  --      --  --
             !
             !                  --  -- 2    --     -- 2
             !                 |   1  |    |         |
             !Pbl_Shear(i,k) = | -----|  * | DU + DV | 
             !                 |  DZ  |    |         |
             !                  --  --      --     --
             !
             !                         --                  -- 2    --     -- 2
             !                 1      | g       si(k+1)      |    |         |          1.0e0_r8
             !Pbl_Shear(i,k) =----- * |--- * ----------------|  * | Du + Dv |   * ---------------------
             !                 T*T    | R     sl(k) -sl(k+1) |    |         |       SIN(colrad(i))**2
             !                         --                  --      --     --
             !
             ! --  -- 2              --                       -- 2
             !|   1  |     1          | g            si(k+1)         | 
             !| -----|  = ----- * |--- * ----------------| 
             !|  DZ  |     T*T    | R          sl(k) -sl(k+1) | 
             ! --  --                --                       --  
             !
             ! --  --               --                       -- 
             !|   1  |     1          | g            si(k+1)         | 
             !| -----|  = ----- * |--- * ----------------| = SQRT((Pbl_ITemp(i,k) * con1(k)))
             !|  DZ  |     T      | R          sl(k) -sl(k+1) | 
             ! --  --                --                       --  
             !
             !                    --                        -- 
             !       1           | g           si(k+1)          |                  1
             !DZ = ----- *       |--- * ----------------| = ---------------------------------
             !       T           | R         sl(k) -sl(k+1)   |   SQRT((Pbl_ITemp(i,k) * con1(k)))
             !                    --                        --  
             !
             !                  1.0e0_r8
             ! csqiv (i) = ---------------------
             !               SIN(colrad(i))**2
             !
             Pbl_Shear(i,k)=(PK_con1(i,k))*Pbl_ITemp(i,k)*((um1(i,(plev+1-k))-um1(i,(plev+1-k)-1))**2+&
             (vm1(i,(plev+1-k))-vm1(i,(plev+1-k)-1))**2)
             Pbl_Shear(i,k)=MAX(gln(i),Pbl_Shear(i,k))          
          END DO
       END DO
       !     
       !     Pbl_BRich(4)        richardson number
       !     Pbl_NRich(5)   flux richardson number
       !     Pbl_NRich(8)   flux richardson number
       !     
       Pbl_KmMixl=0.0_r8 ;Pbl_Km =0.0_r8
       Pbl_KhMixl=0.0_r8 ;Pbl_Kh =0.0_r8
       DO k = 1,(plev-1)
          DO i = 1, plond
             !                g         D (Theta)
             !              ------- * -------------
             !               Theta      D (Z)
             ! Pbl_BRich = -----------------------------
             !                 --  -- 2    --  -- 2
             !                |  dU  |    |  dV  |
             !                | -----|  + | -----| 
             !                |  DZ  |    |  DZ  |
             !                 --  --      --  --
             !
             Pbl_BRich(i,k)=Pbl_Stabil(i,k)/Pbl_Shear(i,k)
             !r1 = 0.5_r8*gama/alfa
             Pbl_NRich(i,k)= r1*(Pbl_BRich(i,k)+r2 &
                  -SQRT(Pbl_BRich(i,k)*(Pbl_BRich(i,k)-r3)+r4))
             Pbl_NRich(i,k)=MIN(rfc,Pbl_NRich(i,k))
             Pbl_NRich(i,k)=Pbl_NRich(i,k)
             !
             !    Pbl_SmBar and Pbl_ShBar are momentum flux and heat flux 
             !    stability parameters, respectively
             !     
             !     Pbl_ShBar(3)   shbar
             !     Pbl_SmBar(4)   smbar
             !     
             !     eliminate negative value for s1-s2*gwrk(i,1,5):
             !     gwrk(i,1,5) is s1/s2 under some circumstances
             !     which makes this expression zero.  machine roundoff
             !     can produce an unphysical negative value in this case.
             !     it is used as sqrt argument in later loop.
             !     s1  =3.0_r8*a2* gam1
             !     s2  =3.0_r8*a2*(gam1+gam2)
             !
             !     s1-s2 = 3.0*a2* gam1 - 3.0*a2*(gam1+gam2)
             !     s1-s2 = 3.0*a2* (gam1 - (gam1+gam2))
             !     
             s1ms2g=s1-s2*Pbl_NRich(i,k)
             IF (ABS(s1ms2g) < 1.0e-10_r8) s1ms2g=0.0_r8
             !
             !                     s1ms2g
             !Pbl_ShBar(i,k)=---------------------------
             !                 (1.0_r8-Pbl_NRich(i,k))
             !
             ! a1    =   0.92_r8
             ! a2    =   0.74_r8
             ! b1    =  16.6_r8
             ! b2    =  10.1_r8
             !
             !         1.0               a1
             !gam1 = -------  - 2.0  * ------
             !         3.0               b1
             !
             !       b2          a1
             !gam2=(----) + (6*------)
             !       b1          b1
             !
             !gam2=(b2  + 6.0*a1)
             !     --------------
             !          b1
             !                 --                   --
             !                | 1.0               a1  |
             !s1 = 3.0 * a2 * |------  - 2.0  * ------|
             !                | 3.0               b1  |
             !                 --                   --
             !s2  =3.0_r8*a2*(gam1+gam2)
             !
             !                   s1-s2*Pbl_NRich(i,k)
             !Pbl_ShBar(i,k)=---------------------------
             !                 (1.0_r8-Pbl_NRich(i,k))
             !
             !                        ((gam1 - (gam1+gam2)))*Pbl_NRich(i,k)
             !Pbl_ShBar(i,k)=3.0*a2* ---------------------------------------- =Sm
             !                           (1.0_r8-Pbl_NRich(i,k))

             Pbl_ShBar(i,k)=s1ms2g/(1.0_r8-Pbl_NRich(i,k))
             !     
             !     gwrk(i,1,3)=(s1-s2*gwrk(i,1,5))/(1.0_r8-gwrk(i,1,5))
             !     end of  negative sqrt argument trap
             ! c1    =   0.08_r8
             ! alfa=b1*(gam1-c1)+3.0_r8*(a2+2.0_r8*a1)
             ! beta=b1*(gam1-c1)
             !        --                              --
             !       | a2                               |
             ! gama=-|---- *(b1*(gam1+gam2)-3.0_r8*a1)  |
             !       | a1                               |
             !        --                              --
             !       --          --
             !      | a2           |
             ! dela=|---- b1* gam1 |
             !      | a1           |
             !       --          --
             !                                     ((b1*(gam1-c1))-(b1*(gam1-c1)+3.0*(a2+2.0*a1))*Pbl_NRich(i,k))
             ! Pbl_SmBar(i,k)=Pbl_ShBar(i,k) * ------------------------------------------------------------------
             !                                     --          --     --                           --
             !                                    | a2           |   | a2                            |
             !                                    |---- b1* gam1 | - |---- *(b1*(gam1+gam2)-3.0*a1)  |*Pbl_NRich(i,k))
             !                                    | a1           |   | a1                            |
             !                                     --          --     --                           --
             !
             !                                     ((b1*(gam1-c1)) - (b1*(gam1-c1) + 3.0*a2 + 6.0*a1))*Pbl_NRich(i,k))
             ! Pbl_SmBar(i,k)=Pbl_ShBar(i,k) * --------------------------------------------------------------
             !                                     --  --     --                             --
             !                                    | a2   |   |                                 |
             !                                    |----  | * |b1* gam1 -(b1*(gam1+gam2)-3.0*a1)|*Pbl_NRich(i,k))
             !                                    | a1   |   |                                 |
             !                                     --  --     --                              --

             !   
             Pbl_SmBar(i,k)=Pbl_ShBar(i,k)*(beta-alfa*Pbl_NRich(i,k))/ &
                  (dela-gama*Pbl_NRich(i,k))
             !     
             !     
             !     Pbl_Sqrtw(5)   sqrt(w) or b/l
             !     Pbl_KmMixl(4)   km/l**2
             !     Pbl_KhMixl(3)   kh/l**2
             ! 
             !     The ratio of SH to SM is equal to the ratio of the turbulent
             !     flux Richardson number to the bulk (large scale) Richardson
             !     number. 
             !
             !     u^2 =  (1-2*gam1)q^2
             !     
             !     SQRT(B1*SM*(1 - Ri)*Shear)
             !             GH    SM
             !     Ri = - ---- =-----Rf
             !             GM    SH
             !             GH
             !     GM = - ----
             !             Ri   
             !                --                  --  
             !               | --  -- 2    --  -- 2 |
             !          l^2  ||  dU  |    |  dV  |  |
             !     GM =-----*|| -----|  + | -----|  |
             !          q^2  ||  DZ  |    |  DZ  |  |
             !               | --  --               |
             !                --                  --  
             ! 
             !                          --      -- 2 
             !             l^2         |  dTHETA  |  
             !     GH = - -----*beta*g*| ---------|  
             !             q^2         |  DZ      |  
             !                          --      --   

             !
             !          --                       --  1/2
             !         |      --  -- 2    --  -- 2 |
             !     q   | 1   |  dU  |    |  dV  |  |
             !    ----=|----*| -----|  + | -----|  |
             !     l   | GM  |  DZ  |    |  DZ  |  |
             !         |      --  --               |
             !          --                       --   
             !                                               
             !  1  
             ! ---- = B1*Sm*(1 - Rf)
             !  GM 
             !
             Pbl_Sqrtw (i,k)=SQRT(b1*Pbl_SmBar(i,k)*(1.0_r8-Pbl_NRich(i,k))*Pbl_Shear(i,k))
             !
             ! Km = l*q*Sm
             !  
             ! Km     q 
             !---- = ---*Sm 
             ! l^2    l
             !
             ! Kh = l*q*Sh
             !
             !  
             ! Kh     q 
             !---- = ---*Sh 
             ! l^2    l
             !
             Pbl_KmMixl(i,k)=Pbl_Sqrtw(i,k)*Pbl_SmBar(i,k)
             Pbl_KhMixl(i,k)=Pbl_Sqrtw(i,k)*Pbl_ShBar(i,k)
          END DO
       END DO
      !          
      !          Pbl_HgtLyI(1)   height at the layer interface 
      !
      !             R*T
      ! DZ = ------ * DP
      !             g*P
      !
      !             R*T      DP
      ! DZ = ------ * -----
      !              g       P
      !
      !              R      DP
      ! DZ = ------ * ----- * T
      !              g       P
      !
      !           
      ! DZ = con0(k)  * T
      !
      !                                                   kg*m*m    
      !                                                  --------  
      !                   R          DP           (j/kg/k)           kg*K*s*s  
      ! con0(k) = ------ * ----- =  --------   =   ------------         = m/K
      !                   g           P            m/s*s              m       
      !                                                    -----    
      !                                                     s*s          
      !                 gasr             si(k) - si(k+1)
      ! con0(k)= -------- * ------------------------
      !                 grav                  sl(k)
      !
      ! Pbl_HgtLyI = con0(k)  * T
      !
      k=1
      DO i = 1, plond
         Pbl_HgtLyI(i,k)=PK_con0(i,k)* tm1 (i,plev)
      END DO
    
      DO k = 2, plev
         DO i = 1, plond
             Pbl_HgtLyI(i,k)=Pbl_HgtLyI(i,k-1)+PK_con0(i,k)*tm1 (i,(plev+1-k))
         END DO
      END DO



       DO itr = 1, nitr
          !     
          !     Pbl_EddEner(2)   mixing length
          !     Pbl_EddEner(2)   b     :b**2 is eddy enegy
          !     
          !..gl0    maximum mixing length l0 in blackerdar's formula
          !         this is retained as a first guess for next time step
          !                                      k0*z
          !                                l = --------
          !                                    (1 + k0*z/l0)

          gln = 0.0_r8
          gld = 0.0_r8
          Pbl_EddEner(1:plond,1) = 0.0_r8
          DO k = 1, plev-1
             DO i = 1, plond
                !
                !                 vk0*gl0(i)*Z(i,k) 
                ! Pbl_EddEner = ------------------------------
                !                (gl0(i) + vk0*Z(i,k))
                !
                Pbl_EddEner(i,k+1)=VKARMAN*gl0(i)*Pbl_HgtLyI(i,k) / (gl0(i)+VKARMAN*Pbl_HgtLyI(i,k))
                !
                !                 vk0*gl0(i)*Z(i,k)           q          q
                ! Pbl_EddEner = ------------------------ * --------- = ---------
                !                (gl0(i) + vk0*Z(i,k))        l^2        l
                !
                Pbl_EddEner(i,k+1)=Pbl_EddEner(i,k+1)*Pbl_Sqrtw(i,k)
             END DO
          END DO
          k=1
          DO i = 1, plond
             Pbl_EddEner(i,k)= 1.0e-3_r8
          END DO
          k=1
          DO i = 1, plond
              !x=0.5_r8*delsig2(k)*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
             !x=0.5_r8*((prsi(i,(plev+1-k)+1)/prsi(i,nsurf+1)) - (prsi(i,(plev+1-k))/prsi(i,nsurf+1)))*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
             x=0.5_r8*((phalf(i,k)/phalf(i,1)) - (phalf(i,k+1)/phalf(i,1)))*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))

             y=x*0.5_r8*Pbl_HgtLyI(i,k)
             gld(i)=gld(i)+x
             gln(i)=gln(i)+y
          END DO
          k=plev
          DO i = 1, plond
             !x=0.5_r8*delsig2(k)*Pbl_EddEner(i,k)
             !x=0.5_r8*((prsi(i,(plev+1-k)+1)/prsi(i,nsurf+1)) - (prsi(i,(plev+1-k))/prsi(i,nsurf+1)))*Pbl_EddEner(i,k)
             x=0.5_r8*((phalf(i,k)/phalf(i,1)) - (phalf(i,k+1)/phalf(i,1)))*Pbl_EddEner(i,k)
             y=x*0.5_r8*(Pbl_HgtLyI(i,k)+Pbl_HgtLyI(i,k-1))
             gln(i)=gln(i)+y
             gld(i)=gld(i)+x
          END DO
          IF (plev > 2) THEN
             DO k = 2, plev-1
                DO i = 1, plond
                   !x=0.5_r8*delsig2(k)*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
                   !x=0.5_r8*((prsi(i,(plev+1-k)+1)/prsi(i,nsurf+1)) - (prsi(i,(plev+1-k))/prsi(i,nsurf+1)))*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
                   x=0.5_r8*((phalf(i,k)/phalf(i,1)) - (phalf(i,k+1)/phalf(i,1)))*(Pbl_EddEner(i,k)+Pbl_EddEner(i,k+1))
                   y=x*0.5_r8*(Pbl_HgtLyI(i,k)+Pbl_HgtLyI(i,k-1))
                   gln(i)=gln(i)+y
                   gld(i)=gld(i)+x
                END DO
             END DO
          END IF
          DO i = 1, plond
             !                0.5 * del(k) * (q(i,k)+q(i,k+1)) * 0.5*(Z(i,k)+Z(i,k-1))
             !gl0(i)=facl * -----------------------------------------------------------
             !                        0.5 * del(k) * (q(i,k)+q(i,k+1))
             gl0(i)=facl*gln(i)/gld(i)
          END DO
          !     
          !     iteration that determines mixing length
          !     
       END DO
       !
       !  Pbl_MixLgh(5)   mixing length
       !
       Pbl_MixLgh=0.0_r8
       DO k = 1, plev-1
          DO i=1,plond
                 !           vk0*gl0(i)*Z(i,k)             m^2
                 ! l  = -------------------------------- = ------
                 !           gl0(i) + vk0*Z(i,k)       m
                 Pbl_MixLgh(i,k)=VKARMAN*gl0(i)*Pbl_HgtLyI(i,k)/(gl0(i)+VKARMAN*Pbl_HgtLyI(i,k))
          END DO
       END DO
       !
       !     PBLH(5)   HEIGHT OF THE PBL
       !
       PBLH=gl0
       tkemyj=0.0_r8
       DO k = 2, plev
          DO i=1,plond
             !
!             tkemyj(i,k)=Pbl_EddEner(i,K)**2 
             tkemyj(i,k)=MIN(Pbl_EddEner(i,K)**2,12.0_r8)
             !                         0.02_r8
             IF(tkemyj(i,k) > EPSQ2*FH)THEN
                PBLH(i)=MIN(MAX(Pbl_HgtLyI(i,k),gl0(i)),3500.0_r8)
             END IF
             tkemyj(i,k)=tkemyj(i,k)/2.0_r8
          END DO
       END DO
       DO k = 1, plev
          DO i=1,plond
             tkemyj(i,k)=MIN(Pbl_EddEner(i,K)**2,12.0_r8)
             tkemyj(i,k)=tkemyj(i,k)/2.0_r8
          END DO
       END DO       
       
       DO i=1,plond
          tkemyj(i,plev+1)=tkemyj(i,plev)
       END DO
       DO k = 1, plev+1
          DO i=1,plond
             MixLgh(i,plev+2-k)=tkemyj(i,k)
          END DO
       END DO
       tkemyj=MixLgh
       !     
       !     Pbl_CoefKm(1)   km = l*q*Sm
       !     Pbl_CoefKh(2)   kh = l*q*Sh
       ! 
       !     where KM and KH are the diffusion coefficients for momentum 
       !     and heat, respectively, l is the master turbulence length scale, 
       !     q2 is the turbulent kinetic energy (so q is the magnitude of 
       !     the turbulent wind velocity), and SM and SH are momentum 
       !     flux and heat flux stability parameters, respectively    
       !    
       IF (kmean == 0) THEN
          DO k = 1, plev-1
             DO i = 1, plond
                !
                !             Km 
                !Km = l^2 * -----
                !            l^2
                !
                KmMixl(i,plev+1-k)=MIN(gkm1,MAX(gkm0,Pbl_MixLgh(i,k)**2*Pbl_KmMixl(i,k)))
                !
                !             Kh 
                !Kh = l^2 * -----
                !            l^2
                !
                KhMixl(i,plev+1-k)=MIN(gkh1,MAX(gkh0,Pbl_MixLgh(i,k)**2*Pbl_KhMixl(i,k)))
             END DO
          END DO
       ELSE
          DO k = 1, plev-1
             DO i = 1, plond
                !
                !             Km 
                !Km = l^2 * -----
                !            l^2
                !
                Pbl_KM(i,plev+1-k)=Pbl_MixLgh(i,k)**2*Pbl_KmMixl(i,k)
                !
                !             Kh 
                !Kh = l^2 * -----
                !            l^2
                !
                Pbl_KH(i,plev+1-k)=Pbl_MixLgh(i,k)**2*Pbl_KhMixl(i,k)
             END DO
          END DO
          fac=0.25_r8
          IF (plev >= 4) THEN
             DO k = 2, plev-2
                DO i = 1, plond
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
                   KmMixl(i,k)=fac*(Pbl_KM(i,k-1)+2.0_r8*Pbl_KM(i,k)+Pbl_KM(i,k+1))
                   !
                   !       Kh(k-1) + 2*Kh(k) + Kh(k+1)
                   ! Kh = -------------------------------
                   !                  4
                   !
                   KhMixl(i,k)=fac*(Pbl_KH(i,k-1)+2.0_r8*Pbl_KH(i,k)+Pbl_KH(i,k+1))
                END DO
             END DO
          END IF
          DO i = 1, plond
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
             KmMixl(i,     1)=0.5_r8*(Pbl_KM(i,     1)+Pbl_KM(i,     2))
             KhMixl(i,     1)=0.5_r8*(Pbl_KH(i,     1)+Pbl_KH(i,     2))
             KmMixl(i,plev-1)=0.5_r8*(Pbl_KM(i,plev-1)+Pbl_KM(i,plev-2))
             KhMixl(i,plev-1)=0.5_r8*(Pbl_KH(i,plev-1)+Pbl_KH(i,plev-2))
          END DO
          DO k = 1, plev-1
             DO i = 1, plond
                KmMixl(i,k)=MIN(gkm1,MAX(gkm0,KmMixl(i,k))) !(m/sec**2)
                KhMixl(i,k)=MIN(gkh1,MAX(gkh0,KhMixl(i,k))) !(m/sec**2)
             END DO
          END DO 
       END IF

 END SUBROUTINE  MY20_TURB    
     
     
  SUBROUTINE MY25_TURB(&
       nCols    ,& ! nCols
       kMax     ,& ! kMax
       delt     ,& ! delt    !       delt     -  Time step in seconds
       um       ,& ! um       (1:nCols,1:kMax)!Zonal wind components
       vm       ,& ! vm       (1:nCols,1:kMax)!Meridional Wind components
       theta    ,& ! theta    (1:nCols,1:kMax)!theta         -  Potential temperature
       zhalf    ,& ! zhalf    (1:nCols,1:kMax+1)!zfull         -  Height at full levels
       zfull    ,& ! zfull    (1:nCols,1:kMax)!         zhalf    -  Height at half levels
       phalf    ,& ! phalf    (1:nCols,1:kMax+1)!!          phalf    -  Pressure at half levels, 
       pfull    ,& ! pfull    (1:nCols,1:kMax)!!          pfull    -  Pressure at full levels! , zhalf, zfull
       z0       ,& ! z0       (1:nCols)!          z0           -  Roughness length
       fracland ,& ! fracland (1:nCols)!        fracland -  Fractional amount of land beneath a grid box
       ustar    ,& ! ustar    (1:nCols)!        ustar         -  OPTIONAL:friction velocity (m/sec)
       el0      ,& ! el0      (1:nCols)!       el0  -  characteristic length scale
       gl0      ,& ! gl0      (1:nCols)!       el0  -  characteristic length scale
       PBLH     ,& ! el0      (1:nCols)!       el0  -  characteristic length scale
       akm      ,& ! akm      (1:nCols,1:kMax+1)!       akm  -  mixing coefficient for momentum
       akh      ,& ! akh      (1:nCols,1:kMax+1)!       akh  -  mixing coefficient for heat and moisture
       el       ,& ! el       (1:nCols,1:kMax+1)!       el   -  master length scale
       TKE      ) ! tke       (1:nCols,1:kMax+1)!       tke  -  master length scale
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: delt    !       delt     -  Time step in seconds
    REAL(KIND=r8)   , INTENT(in   ) :: um       (1:nCols,1:kMax)!Zonal wind components
    REAL(KIND=r8)   , INTENT(in   ) :: vm       (1:nCols,1:kMax)!Meridional Wind components
    REAL(KIND=r8)   , INTENT(in   ) :: theta    (1:nCols,1:kMax)!theta    -  Potential temperature
    REAL(KIND=r8)   , INTENT(in   ) :: zfull    (1:nCols,1:kMax)!zfull    -  Height at full levels
    REAL(KIND=r8)   , INTENT(in   ) :: zhalf    (1:nCols,1:kMax+1)!       zhalf    -  Height at half levels
    REAL(KIND=r8)   , INTENT(in   ) :: phalf    (1:nCols,1:kMax+1)!!       phalf    -  Pressure at half levels, 
    REAL(KIND=r8)   , INTENT(in   ) :: pfull    (1:nCols,1:kMax)!!       pfull    -  Pressure at full levels! , zhalf, zfull
    REAL(KIND=r8)   , INTENT(in   ) :: z0       (1:nCols)   !          z0          -  Roughness length
    REAL(KIND=r8)   , INTENT(in   ) :: fracland (1:nCols)!       fracland -  Fractional amount of land beneath a grid box
    REAL(KIND=r8)   , INTENT(in   ) :: ustar    (1:nCols)!       ustar    -  OPTIONAL:friction velocity (m/sec)
    REAL(KIND=r8)   , INTENT(in   ) :: gl0      (1:nCols)!       
    REAL(KIND=r8)   , INTENT(out  ) :: el0      (1:nCols)
    REAL(KIND=r8)   , INTENT(out  ) :: PBLH     (1:nCols)
    REAL(KIND=r8)   , INTENT(out  ) :: akm      (1:nCols,1:kMax+1)
    REAL(KIND=r8)   , INTENT(out  ) :: akh      (1:nCols,1:kMax+1)
    REAL(KIND=r8)   , INTENT(out  ) :: el       (1:nCols,1:kMax+1)
    REAL(KIND=r8)   , INTENT(inout) :: TKE      (1:nCols,1:kMax+1)

    !       mask     -  OPTIONAL; floating point mask (0. or 1.) designating
    !                   where data is present

    REAL(KIND=r8)    :: zsfc  (1:nCols)
    REAL(KIND=r8)    :: el2   (1:nCols,1:kMax-1)

    REAL(KIND=r8)    :: dsdzh (1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: qm2   (1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: qm3   (1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: qm4   (1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: buoync(1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: dsdz  (1:nCols,1:kMax) 
    REAL(KIND=r8)    :: xx1   (1:nCols,1:kMax) 
    REAL(KIND=r8)    :: xx2   (1:nCols,1:kMax) 
    REAL(KIND=r8)    :: qm    (1:nCols,1:kMax) 
    REAL(KIND=r8)    :: x1    (1:nCols)
    REAL(KIND=r8)    :: x2    (1:nCols)
    REAL(KIND=r8)    :: akmin (1:nCols)
    REAL(KIND=r8)    :: aaa   (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: bbb   (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: ddd   (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: ccc   (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: xxm1  (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: xxm2  (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: xxm3  (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: xxm4  (1:nCols,1:kMax-1) 
    REAL(KIND=r8)    :: xxm5  (1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: shear (1:nCols,1:kMax-1)
    REAL(KIND=r8)    :: cvfqdt
    REAL(KIND=r8)    :: dvfqdt
    REAL(KIND=r8),    PARAMETER :: epsq2=0.2_r8
    REAL(KIND=r8),    PARAMETER :: FH=1.01_r8

    INTEGER  :: i
    INTEGER  :: k
    INTEGER  :: it
    INTEGER  :: klim
    el0= 0.0_r8;akm= 0.0_r8;akh= 0.0_r8;el =0.0_r8 
    zsfc= 0.0_r8;el2= 0.0_r8;dsdzh= 0.0_r8;qm2= 0.0_r8;qm3= 0.0_r8;qm4= 0.0_r8   
    buoync= 0.0_r8;dsdz= 0.0_r8; xx1 = 0.0_r8; xx2 = 0.0_r8;qm = 0.0_r8;x1  = 0.0_r8;  
    x2  = 0.0_r8;akmin = 0.0_r8;aaa  = 0.0_r8;bbb  = 0.0_r8;ddd  = 0.0_r8;ccc  = 0.0_r8; 
    xxm1 = 0.0_r8;xxm2 = 0.0_r8;xxm3 = 0.0_r8;xxm4  = 0.0_r8;xxm5  = 0.0_r8;shear = 0.0_r8;
    cvfqdt= 0.0_r8;dvfqdt= 0.0_r8
    !====================================================================
    ! --- SURFACE HEIGHT     
    !====================================================================

    !if( PRESENT( kbot ) ) then
    !  do i = 1,nCols
    !      k = kbot(i) + 1
    !     zsfc(i) = zhalf(i,k)
    !  end do
    !else
    DO i = 1,nCols
       zsfc(i) = zhalf(i,kMax+1)
    END DO
    !endif
    !====================================================================
    ! --- D( )/DZ OPERATORS: AT FULL LEVELS & AT HALF LEVELS          
    !====================================================================
    DO k = 1,kMax
       DO i = 1,nCols
          dsdz (i,k)  = 1.0_r8 / ( zhalf(i,k+1) - zhalf(i,k) )
       END DO
    END DO
    DO k = 1,kMax-1
       DO i = 1,nCols
          dsdzh(i,k) = 1.0_r8 / ( zfull(i,k+1)  - zfull(i,k) )
       END DO
    END DO
    !====================================================================
    ! --- WIND SHEAR                 
    !====================================================================
    DO k = 1,kMax-1
       DO i = 1,nCols
          xxm1 (i,k) = dsdzh(i,k) * ( um(i,k+1) - um(i,k) )
          xxm2 (i,k) = dsdzh(i,k) * ( vm(i,k+1) - vm(i,k) )
          shear(i,k) = xxm1 (i,k) * xxm1(i,k) + xxm2(i,k) * xxm2(i,k)
       END DO
    END DO

    !====================================================================
    ! --- BUOYANCY                 
    !====================================================================
    DO k = 1,kMax-1
       DO i = 1,nCols
          xxm1(i,k) = theta(i,k+1) - theta(i,k) 
          IF( do_thv_stab ) THEN
             xxm2(i,k) = 0.5_r8*( theta(i,k+1) + theta(i,k) )
          ELSE
             xxm2(i,k) = t00
          END IF

          buoync(i,k) = grav * dsdzh(i,k) * xxm1(i,k) / xxm2(i,k)

       END DO
    END DO

    !====================================================================
    ! --- MASK OUT UNDERGROUND VALUES FOR ETA COORDINATE
    !====================================================================
    !IF(jdt<=2)
     CALL TKE_SURF ( nCols                            , &
                     kMax                             , &
                     ustar(1:nCols)                   , &
                     TKE  (1:nCols,1:kMax+1)    )
    !  if( PRESENT( mask ) ) then
    !    do k=2,kMax
    !       do i=1,nCols
    !         if(mask(i,k) < 0.1_r8) then
    !            TKE(ism+i,jsm+j,k+1) = 0.0_r8
    !            dsdz(i,j,k  ) = 0.0_r8
    !            dsdzh(i,j,k-1) = 0.0_r8
    !            shear(i,j,k-1) = 0.0_r8
    !            buoync(i,j,k-1) = 0.0_r8
    !         endif
    !       enddo
    !      enddo
    !    enddo
    !  endif

    !====================================================================
    ! --- SET ITERATION LOOP IF INITALIZING TKE
    !====================================================================


    ! $$$$$$$$$$$$$$$$$
    DO it = 1,init_iters
       ! $$$$$$$$$$$$$$$$$

       !====================================================================
       ! --- SOME TKE STUFF
       !====================================================================

       DO k=1,kMax
          DO i=1,nCols
             xx1(i,k) = MAX(2*TKE(i,k+1),1.0e-12_r8)
             IF(xx1(i,k) > 0.0_r8) THEN
                qm(i,k) = SQRT(xx1(i,k))
             ELSE
                qm(i,k) = 0.0001_r8
             END IF
          END DO
       END DO

       DO k=1,kMax-1
          DO i=1,nCols
             qm2(i,k)  = xx1(i,k) 
             qm3(i,k)  = qm (i,k) * qm2(i,k) 
             qm4(i,k)  = qm2(i,k) * qm2(i,k) 
          END DO
       END DO

       !====================================================================
       ! --- CHARACTERISTIC LENGTH SCALE                         
       !====================================================================
       DO k=1,kMax-1
          DO i=1,nCols
!             xx1(i,k) = qm(i,k)*( pfull(i,k+1) - pfull(i,k) )
             xx1(i,k) = qm(i,k) * ( MAX(pfull(i,k+1) - pfull(i,k),0.0000001_r8 ))
          END DO
       END DO

       DO k = 1, kMax-1
          DO i=1,nCols
             xx2(i,k) = xx1(i,k)  * ( zhalf(i,k+1) - zsfc(i) )
          END DO
       END DO

       !if( PRESENT( kbot ) ) then
       !  xx1(:,:,kMax) = 0.0_r8
       !  xx2(:,:,kMax) = 0.0_r8
       !   do i = 1,nCols
       !     k = kbot(i) 
       !     xx1(i,k)  =  qm(i,k) * ( phalf(i,k+1) - pfull(i,k) )
       !     xx2(i,k)  = xx1(i,k) * z0(i)
       !   end do
       !else
       DO i = 1,nCols
!          xx1(i,kMax) =  qm(i,kMax) * ( MAX(phalf(i,kMax+1) - phalf(i,kMax),1.0_r8) )
          xx1(i,kMax) =  qm(i,kMax) * ( MAX(phalf(i,kMax+1) - phalf(i,kMax),0.00000000000001_r8) )
          xx2(i,kMax) = xx1(i,kMax) * z0(i)
       END DO
       !end if

       !  if (PRESENT(mask)) then
       !    x1 = SUM( xx1, 3, mask=mask.gt.0.1_r8 )
       !    x2 = SUM( xx2, 3, mask=mask.gt.0.1_r8 )
       !  else
       DO i = 1,nCols
          x1(i) = SUM( xx1(i,1:kMax))
          x2(i) = SUM( xx2(i,1:kMax))
          IF(x1(i) <= 0.0_r8)THEN
            x1(i)=1e-12_r8
          END IF
       END DO

       !  endif

       !---- should never be equal to zero ----

       IF (COUNT(x1 <= 0.0_r8) > 0) THEN
          WRITE(*,*)' MY25_TURB','divid by zero, x1 <= 0.0_r8'
          STOP
       END IF

       DO i = 1,nCols
          el0(i) = x2(i) / x1(i)
          el0(i) = el0(i) * (alpha_land*fracland (i)+ alpha_sea*(1.0_r8-fracland(i)))

          el0(i) = MIN( el0(i), el0max1 )
          el0(i) = MAX( el0(i), el0min1 )
       END DO
       !====================================================================
       ! --- MASTER LENGTH SCALE 
       !====================================================================

       DO k = 1, kmax-1
          DO i = 1,nCols
             xx1(i,k)  = vonkarm * ( zhalf(i,k+1) - zsfc(i) )
          END DO
       END DO
       DO i = 1,nCols
          x1(i) = vonkarm * z0(i) 
       END DO

       !if( PRESENT( kbot ) ) then
       !    do i = 1,nCols
       !      do k = kbot(i), kMax
       !        xx1(i,j,k) = x1(i,j)
       !      end do
       !    end do
       !else
       DO i = 1,nCols
          xx1(i,kMax) = x1(i)
       END DO
       !endif 

       DO k = 1,kMax
          DO i = 1,nCols
             el(i,k+1) = xx1(i,k) / ( 1.0_r8 + xx1(i,k) / el0(i) )
          END DO
       END DO
       DO i = 1,nCols
          el(i,1)   = el0(i)
       END DO
       DO k = 1,kMax-1
          DO i = 1,nCols
             el2(i,k) = el(i,k+1) * el(i,k+1) 
          END DO
       END DO
       !====================================================================
       ! --- MIXING COEFFICIENTS                     
       !====================================================================
       DO k = 1,kMax-1    
          DO i = 1,nCols
             xxm3(i,k) = el2(i,k)* buoync(i,k)
             xxm4(i,k) = el2(i,k)* shear(i,k)
             xxm5(i,k) = el (i,k+1)*qm3(i,k)
          END DO
       END DO
       !-------------------------------------------------------------------
       ! --- MOMENTUM 
       !-------------------------------------------------------------------
       DO k = 1,kMax-1
          DO i = 1,nCols
             xxm1(i,k) = xxm5(i,k)*( ckm1*qm2(i,k) + ckm2*xxm3(i,k) )
             xxm2(i,k) = qm4(i,k) + ckm5*qm2(i,k)*xxm4(i,k) + xxm3(i,k)*( ckm6*xxm4(i,k) + ckm7*qm2(i,k) + ckm8*xxm3(i,k) )

             xxm2(i,k) = MAX( xxm2(i,k), 0.2_r8*qm4(i,k) )
             xxm2(i,k) = MAX( xxm2(i,k), small2  )

             akm(i,1)    = 0.0_r8
             akm(i,k+1)  = xxm1(i,k) / xxm2(i,k)

             akm(i,k+1) = MAX( akm(i,k+1), 0.0_r8 )

             !-------------------------------------------------------------------
             ! --- HEAT AND MOISTURE 
             !-------------------------------------------------------------------

             xxm1(i,k) = ckh1*xxm5(i,k) - ckh2*xxm4(i,k)*akm(i,k+1)
             xxm2(i,k) = qm2      (i,k) + ckh3*xxm3(i,k)

             xxm1(i,k) = MAX( xxm1(i,k), ckh4*xxm5(i,k) )
             xxm2(i,k) = MAX( xxm2(i,k),  0.4_r8*qm2(i,k)   )
             xxm2(i,k) = MAX( xxm2(i,k), small2     )

             akh(i,1)    = 0.0_r8
             akh(i,k+1) = xxm1(i,k) / xxm2(i,k)

             akh(i,k+1) = MAX( akh(i,k+1), 0.0_r8 )
             !-------------------------------------------------------------------
             ! --- BOUNDS 
             !-------------------------------------------------------------------

             ! --- UPPER BOUND
             !akm(i,k) = MIN( akm(i,k), akmax )
             !akh(i,k) = MIN( akh(i,k), akmax )
          END DO
       END DO
       DO k = 1,kMax
          DO i = 1,nCols
             !-------------------------------------------------------------------
             ! --- BOUNDS 
             !-------------------------------------------------------------------

             ! --- UPPER BOUND
             akm(i,k) = MIN( akm(i,k), akmax )
             akh(i,k) = MIN( akh(i,k), akmax )
          END DO
       END DO

       ! --- LOWER BOUND 
       !klim = kMax - nk_lim + 1
       !DO  k = klim,kMax+1
       !   DO i = 1,nCols
       !      IF(akm(i,k) < small2 )  akm(i,k) = 0.0_r8
       !      IF(akh(i,k) < small2 )  akh(i,k) = 0.0_r8
       !   END DO
       !END DO

       ! --- LOWER BOUND NEAR SURFACE
       DO i = 1,nCols
          akmin(i) = akmin_land*fracland(i) + akmin_sea*(1.0_r8-fracland(i))
       END DO
       !DO i = 1,nCols
       !      !
       !      !       Km(1) + Km(2)
       !      ! Km = ---------------
       !      !            2
       !      !
       !      !
       !      !       Kh(1) + Kh(2)
       !      ! Kh = ---------------
       !      !            2
       !      !
       !      akm(i,     1)=0.5_r8*(akm(i,     1)+akm(i,     2))
       !      akh(i,     1)=0.5_r8*(akh(i,     1)+akh(i,     2))
       !      akm(i,kMax-1)=0.5_r8*(akm(i,kMax-1)+akm(i,kMax-2))
       !      akh(i,kMax-1)=0.5_r8*(akh(i,kMax-1)+akh(i,kMax-2))
       !      akm(i,kMax  )=0.5_r8*(akm(i,kMax  )+akm(i,kMax-1))
       !      akh(i,kMax  )=0.5_r8*(akh(i,kMax  )+akh(i,kMax-1))
       !      akm(i,kMax+1)=0.5_r8*(akm(i,kMax+1)+akm(i,kMax))
       !      akh(i,kMax+1)=0.5_r8*(akh(i,kMax+1)+akh(i,kMax))
       !END DO

       !if( PRESENT( kbot ) ) then
       !      do i = 1,nCols
       !           klim = kbot(i) - nk_lim + 1
       !         do  k = klim,kbot(i)
       !             akm(i,k) = MAX( akm(i,k), akmin(i) )
       !             akh(i,k) = MAX( akh(i,k), akmin(i) )
       !         end do
       !      end do
       !else
       klim = kMax - nk_lim + 1
       DO  k = klim,kMax
          DO i = 1,nCols
             akm(i,k) = MAX( akm(i,k), akmin(i) )
             akh(i,k) = MAX( akh(i,k), akmin(i) )
          END DO
       END DO
       !endif

       !-------------------------------------------------------------------
       ! --- MASK OUT UNDERGROUND VALUES FOR ETA COORDINATE
       !-------------------------------------------------------------------

       !if( PRESENT( mask ) ) then
       !    akm(:,:,1:kx) = akm(:,:,1:kx) * mask(:,:,1:kx) 
       !    akh(:,:,1:kx) = akh(:,:,1:kx) * mask(:,:,1:kx) 
       !endif

       !====================================================================
       ! --- PROGNOSTICATE TURBULENT KE 
       !====================================================================

       cvfqdt = cvfq1 * delt
       dvfqdt = cvfq2 * delt * 2.0_r8

       !-------------------------------------------------------------------
       ! --- PART OF LINEARIZED ENERGY DISIIPATION TERM 
       !-------------------------------------------------------------------
       DO k = 1,kMax-1
          DO i = 1,nCols
             xxm1(i,k) = dvfqdt * qm(i,k) / el(i,k+1)
          END DO
       END DO

       !-------------------------------------------------------------------
       ! --- PART OF LINEARIZED VERTICAL DIFFUSION TERM
       !-------------------------------------------------------------------
       DO k = 1,kMax
          DO i = 1,nCols

             xx1(i,k ) = el(i,k+1) * qm(i,k)
          END DO
       END DO
       DO i = 1,nCols
          xx2(i,1)    = 0.5_r8*  xx1(i,1)
       END DO
       DO k = 1,kMax-1
          DO i = 1,nCols
             xx2(i,k+1 ) = 0.5_r8*( xx1(i,k+1 ) + xx1(i,k) )
          END DO
       END DO
       DO k = 1,kMax
          DO i = 1,nCols
             xx1 (i,k )= xx2(i,k ) * dsdz(i,k )
          END DO
       END DO

       !-------------------------------------------------------------------
       ! --- IMPLICIT TIME DIFFERENCING FOR VERTICAL DIFFUSION 
       ! --- AND ENERGY DISSIPATION TERM 
       !-------------------------------------------------------------------

       DO k=1,kMax-1
          DO i=1,nCols
             aaa(i,k) = -cvfqdt * xx1(i,k+1) * dsdzh(i,k)
             ccc(i,k) = -cvfqdt * xx1(i,k  ) * dsdzh(i,k)
             bbb(i,k) =     1.0_r8 - aaa(i,k  ) -   ccc(i,k) 
             bbb(i,k) =             bbb(i,k  ) +  xxm1(i,k)
             ddd(i,k) =             TKE(i,k+1)
          END DO
       END DO

       ! correction for vertical diffusion of TKE surface boundary condition

       !if (present(kbot)) then
       !   do i = 1,nCols
       !        k = kbot(i)
       !        ddd(i,k-1) = ddd(i,k-1) - aaa(i,k-1) * TKE(i,latco,k+1)
       !   end do
       !else
       DO i = 1,nCols
          ddd(i,kMax-1) = ddd(i,kMax-1) - aaa(i,kMax-1) * TKE(i,kMax+1)
       END DO
       !endif

       ! mask out terms below ground

       !  if (present(mask)) then
       !     where (mask(:,:,2:kx) < 0.1_r8) ddd(:,:,1:kxm) = 0.0_r8
       !  endif


       CALL TRI_INVERT( &
            nCols,  & 
            kmax-1, & 
            xxm1  (1:nCols,1:kmax-1)     , & ! intent(out)  :: x (:,:)
            ddd   (1:nCols,1:kmax-1)     , & ! intent(in)   :: d (:,:)
            aaa   (1:nCols,1:kmax-1)     , & ! optional     :: a (:,:)
            bbb   (1:nCols,1:kmax-1)     , & ! optional     :: b (:,:)
            ccc   (1:nCols,1:kmax-1)       ) ! optional     :: c (:,:)


       !-------------------------------------------------------------------
       ! --- MASK OUT UNDERGROUND VALUES FOR ETA COORDINATE
       !-------------------------------------------------------------------

       ! if( PRESENT( mask ) ) then
       !   do k=1,kxm
       !   do j=1,jx
       !   do i=1,ix
       !     if(mask(i,j,k+1) < 0.1_r8) xxm1(i,j,k) = TKE(ism+i,jsm+j,k+1)
       !   enddo
       !   enddo
       !   enddo
       ! endif

       !-------------------------------------------------------------------
       ! --- SHEAR AND BUOYANCY TERMS
       !-------------------------------------------------------------------
       DO k=1,kMax-1
          DO i=1,nCols

             xxm2(i,k) =  delt*( akm(i, k+1 )* shear (i,k)    &
                  - akh(i,k+1  )* buoync(i,k ) )
          END DO
       END DO

       !-------------------------------------------------------------------
       ! --- UPDATE TURBULENT KINETIC ENERGY
       !-------------------------------------------------------------------

       DO i=1,nCols
          TKE(i,1) = 0.0_r8
       ENDDO

       DO k=2,kMax
          DO i=1,nCols
             TKE(i,k) = xxm1(i,k-1) + xxm2(i,k-1)
          ENDDO
       ENDDO

       !====================================================================
       ! --- BOUND TURBULENT KINETIC ENERGY
       !====================================================================
       DO k=1,kMax+1
          DO i=1,nCols
             TKE(i,k) = MIN( TKE(i,k), TKEmax )
             TKE(i,k) = MAX( TKE(i,k), TKEmin )
          ENDDO
       ENDDO     

       !DO k=1,kMax+1
       !  DO i=1,nCols
             !IF(TKE(i,k) == TKEmin) TKE(i,k)=0.0_r8
       !   ENDDO
       !ENDDO     

       !  if( PRESENT( mask ) ) then
       !    do k=1,kx
       !    do j=1,jx
       !    do i=1,ix
       !      if(mask(i,j,k) < 0.1_r8) TKE(ism+i,jsm+j,k+1) = 0.0_r8
       !    enddo
       !    enddo
       !    enddo
       !  endif

       !====================================================================
       ! --- COMPUTE PBL DEPTH IF DESIRED
       !====================================================================

       !  if (present(h)) then

       !      if (.not.present(ustar)) then
       !          CALL ERROR_MESG( ' MY25_TURB',     &
       !              'cannot request pbl depth diagnostic if ustar'// &
       !              ' and bstar are not also supplied', FATAL )
       !      end if
       !if (present(kbot)) then
       !    call k_pbl_depth(ustar (1:nCols)                , &
       !                           akm   (1:nCols,1:kMax+1)       , &
       !                           akh   (1:nCols,1:kMax+1)       , &
       !                           zsfc  (1:nCols)                , &
       !                           zfull (1:nCols,1:kMax)         , &
       !                           zhalf (1:nCols,1:kMax+1)       , &
       !                           h     (1:nCols)                , &
       !                           kbot(1:nCols)                    )
       !      else
       !          call k_pbl_depth(ustar (1:nCols)                  , &
       !                           akm   (1:nCols,1:kMax+1)          , &
       !                           akh   (1:nCols,1:kMax+1)          , &
       !                           zsfc  (1:nCols)                  , &
       !                           zfull (1:nCols,1:kMax)          , &
       !                           zhalf (1:nCols,1:kMax+1)          , &
       !                           h     (1:nCols)                    )
       !      end if

       !  end if
       !====================================================================
       !
       !     PBLH(5)   HEIGHT OF THE PBL
       !
       PBLH=gl0
       DO k =  kmax+1,1,-1
          DO i=1,ncols
             !
             !                         0.02_r8
             IF(TKE(i,k)> EPSQ2*FH)THEN
                PBLH(i)=MIN(MAX(zhalf(i,k),gl0(i)),3500.0_r8)
             END IF
          END DO
       END DO

       ! $$$$$$$$$$$$$$$$$
    END DO
    ! $$$$$$$$$$$$$$$$$

  !#######################################################################
  END SUBROUTINE MY25_TURB

  SUBROUTINE TKE_SURF ( nCols,kMax,u_star,TKE )

    !=======================================================================
    !---------------------------------------------------------------------
    ! Arguments (Intent in)
    !       u_star -  surface friction velocity (m/s)
    !       kbot   -  OPTIONAL;lowest model level index (integer);
    !                 at levels > Kbot, Mask = 0.
    !---------------------------------------------------------------------
    INTEGER , INTENT(in)   :: nCols
    INTEGER , INTENT(in)   :: kMax

    REAL(KIND=r8), INTENT(in)   :: u_star(nCols) 
    REAL(KIND=r8), INTENT(inout)  :: TKE   (nCols,kMax+1) 
    !---------------------------------------------------------------------
    !  (Intent local)
    !---------------------------------------------------------------------
    REAL(KIND=r8) :: x1(nCols)
    INTEGER  :: ix, kxp, i

    !=======================================================================

    ix  = nCols
    kxp = kMax+1
    !TKE=0.0_r8
    !---------------------------------------------------------------------
    DO i = 1,ix
       ! bb1 = 16.0_r8
       ! bcq   = 0.5 * ( 16.0_r8**(2.0/3.0) )
       ! tk2(i) = ( 0.5 * ( 16.0_r8**(2.0/3.0) )) * u_star(i) * u_star(i)
       x1(i) = MAX(bcq * u_star(i) * u_star(i),1.0e-12_r8 )
       TKE(i,kxp) = MAX(x1(i),TKE(i,kxp)) 
    END DO

    !=======================================================================
  END SUBROUTINE TKE_SURF


  SUBROUTINE tri_invert(&
       nCols   , &
       kMax    , &
       x       , &  ! intent(out)  :: x (:,:)
       d       , &  ! intent(in)   :: d (:,:)
       a       , &  ! optional          :: a (:,:)
       b       , &  ! optional          :: b (:,:)
       c         )  ! optional          :: c (:,:)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(out)  :: x  (nCols,kMax)
    REAL(KIND=r8), INTENT(in)   :: d  (nCols,kMax)
    REAL(KIND=r8), OPTIONAL     :: a  (nCols,kMax)
    REAL(KIND=r8), OPTIONAL     :: b  (nCols,kMax)
    REAL(KIND=r8), OPTIONAL     :: c  (nCols,kMax)

    REAL(KIND=r8)  :: f (nCols,kMax)

    REAL(KIND=r8)  :: e (nCols,kMax)
    REAL(KIND=r8)  :: g (nCols,kMax)
    REAL(KIND=r8)  :: bb(nCols)
    REAL(KIND=r8)  :: cc(nCols,kMax)

    INTEGER :: k
    INTEGER :: i
    
    f =0.0_r8;e=0.0_r8;g=0.0_r8;bb=0.0_r8;cc=0.0_r8;x=0.0_r8
    
    DO i=1,nCols
       e(i,1) = - a(i,1)/b(i,1)
       a(i,kMax) = 0.0_r8
    END DO

    DO  k= 2,kMax
       DO i=1,nCols
          g(i,k) = 1/(b(i,k)+c(i,k)*e(i,k-1))
          e(i,k) = - a(i,k)*g(i,k)
       END DO
    END DO
    DO  k= 1,kMax
       DO i=1,nCols
          cc(i,k) = c(i,k)
       END DO
    END DO
    DO i=1,nCols
       bb(i) = 1.0_r8/b(i,1)
    END DO

    ! if(.not.init_tridiagonal) error
    DO i=1,nCols
       f(i,1) =  d(i,1)*bb(i)  
    END DO
    DO k= 2, kMax
       DO i=1,nCols
          f(i,k) = (d(i,k) - cc(i,k)*f(i,k-1))*g(i,k)
       END DO
    END DO
    DO i=1,nCols
       x(i,kMax) = f(i,kMax)
    END DO

    DO k = kMax-1,1,-1
       DO i=1,nCols
          !IF(ABS(e(i,k)) < 1.0e-12_r8)PRINT*,'e(i,k)=',e(i,k)
          !IF(ABS(x(i,k)) < 1.0e-12_r8)PRINT*,'x(i,k)=',x(i,k+1)
          !IF(ABS(f(i,k)) < 1.0e-12_r8)PRINT*,'f(i,k)=',f(i,k)
          x(i,k) = e(i,k)*x(i,k+1)+f(i,k)
       END DO
    END DO

    RETURN
  END SUBROUTINE tri_invert



   END MODULE Pbl_HostlagBoville
