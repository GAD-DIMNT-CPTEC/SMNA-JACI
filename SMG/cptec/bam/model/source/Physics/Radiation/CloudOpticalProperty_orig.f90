MODULE CloudOpticalProperty
 USE PhysicalFunctions,Only : fpvs2es5

    IMPLICIT NONE
  SAVE

  PRIVATE
  INTEGER      , PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER      , PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER      , PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER      , PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER      , PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers
  REAL(kind=r8), PARAMETER :: con_csol   =2.1060e+3_r8      ! spec heat H2O ice   (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_psat   =6.1078e+2_r8      ! pres at H2O 3pt     (Pa)  
  REAL(kind=r8), PARAMETER :: con_pi     =3.1415926535897931_r8 ! pi
  REAL(kind=r8), PARAMETER :: con_cvap   =1.8460e+3_r8      ! spec heat H2O gas   (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_cliq   =4.1855e+3_r8      ! spec heat H2O liq   (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_rv     =4.6150e+2_r8      ! gas constant H2O    (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_fvirt  =con_rv/con_rd-1.0_r8
  REAL(kind=r8), PARAMETER :: con_g      =9.80665e+0_r8     ! gravity           (m/s2)
  REAL(kind=r8), PARAMETER :: con_t0c    =2.7315e+2_r8      ! temp at 0C          (K)
  REAL(kind=r8), PARAMETER :: con_hfus   =3.3358e+5_r8      ! lat heat H2O fusion (J/kg)
  REAL(kind=r8), PARAMETER :: con_hvap   =2.5000e+6_r8      ! lat heat H2O cond   (J/kg)
  REAL(kind=r8), PARAMETER :: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8), PARAMETER :: con_ttp    =2.7316e+2_r8      ! temp at H2O 3pt     (K)
  REAL(kind=r8), PARAMETER :: con_eps    =con_rd/con_rv
  REAL(kind=r8), PARAMETER :: con_epsm1  =con_rd/con_rv-1.0_r8
  REAL(KIND=r8), PARAMETER :: pptop = 0.005_r8       ! Model-top presure                                 
  REAL(KIND=r8), PARAMETER :: tmelt =273.16_r8

  INTEGER      , PARAMETER :: iflip=1
  !   iflip           : control flag for direction of vertical index      !
  !                     =0: index from toa to surface                     !
  !                     =1: index from surface to toa                     !
  !  ---  set constant parameters

  REAL (kind=r8), PARAMETER :: gfac=1.0e5_r8/con_g
  REAL (kind=r8), PARAMETER :: gord=con_g/con_rd
  INTEGER       , PARAMETER :: NF_CLDS = 9   ! number of fields in cloud array

  !  ---  pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
  REAL (kind=r8), PARAMETER :: ptopc(1:4,1:2)=  RESHAPE( (/ 1050.0_r8, 650.0_r8, 400.0_r8, 0.0_r8,  &
       1050.0_r8, 750.0_r8, 500.0_r8, 0.0_r8/),(/4,2/) )

  !org  data ptopc / 1050., 642., 350., 0.0,  1050., 750., 500., 0.0 /
  !DATA ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /(4,2)

  !     real (kind=r8), parameter :: climit = 0.01
  REAL (kind=r8), PARAMETER :: climit_cld = 0.001_r8
  REAL (kind=r8), PARAMETER :: climit2=0.05_r8
  REAL (kind=r8), PARAMETER :: ovcst  = 1.0_r8 - 1.0e-8_r8
  !  ---  set default quantities as parameters (for prognostic cloud)

  REAL (kind=r8), PARAMETER :: reliq_def = 10.0_r8    ! default liq radius to 10 micron
  REAL (kind=r8), PARAMETER :: reice_def = 50.0_r8    ! default ice radius to 50 micron
  REAL (kind=r8), PARAMETER :: rrain_def = 1000.0_r8  ! default rain radius to 1000 micron
  REAL (kind=r8), PARAMETER :: rsnow_def = 250.0_r8   ! default snow radius to 250 micron

  REAL (KIND=r8), PARAMETER   :: rmwmdi =                 1.61e0_r8!
  REAL (KIND=r8), PARAMETER   :: e0c  =                   6.11e0_r8!

  !
  !--- Common block of constants used in column microphysics
  !
  REAL(KIND=r8) ::  ABFR
  REAL(KIND=r8) ::  CBFR
  REAL(KIND=r8) ::  CIACW
  REAL(KIND=r8) ::  CIACR
  REAL(KIND=r8) ::  C_N0r0
  REAL(KIND=r8) ::  CN0r0
  REAL(KIND=r8) ::  CN0r_DMRmin
  REAL(KIND=r8) ::  CN0r_DMRmax
  REAL(KIND=r8) ::  CRACW
  REAL(KIND=r8) ::  CRAUT
  REAL(KIND=r8) ::  ESW0
  REAL(KIND=r8) ::  QAUTx
  REAL(KIND=r8) ::  RFmax
  REAL(KIND=r8) ::  RQR_DR1
  REAL(KIND=r8) ::  RQR_DR2
  REAL(KIND=r8) ::  RQR_DR3
  REAL(KIND=r8) ::  RQR_DRmin
  REAL(KIND=r8) ::  RQR_DRmax
  REAL(KIND=r8) ::  RR_DRmin
  REAL(KIND=r8) ::  RR_DR1
  REAL(KIND=r8) ::  RR_DR2
  REAL(KIND=r8) ::  RR_DR3
  REAL(KIND=r8) ::  RR_DRmax
  !
  INTEGER :: mic_step
  !
  !--- Common block for lookup table used in calculating growth rates of
  !    nucleated ice crystals growing in water saturated conditions
  !--- Discretized growth rates of small ice crystals after their nucleation
  !     at 1 C intervals from -1 C to -35 C, based on calculations by Miller
  !     and Young (1979, JAS) after 600 s of growth.  Resultant growth rates
  !     are multiplied by physics time step in GSMCONST.
  !
  INTEGER, PARAMETER :: MY_T1=1
  INTEGER, PARAMETER :: MY_T2=35
  REAL(KIND=r8)               :: MY_GROWTH(MY_T1:MY_T2)
  !
  !--- Parameters for ice lookup tables, which establish the range of mean ice
  !    particle diameters; from a minimum mean diameter of 0.05 mm (DMImin) to a
  !    maximum mean diameter of 1.00 mm (DMImax).  The tables store solutions
  !    at 1 micron intervals (DelDMI) of mean ice particle diameter.
  !

  REAL(KIND=r8)   , PARAMETER :: DMImin=.05e-3_r8
  REAL(KIND=r8)   , PARAMETER :: DMImax=1.e-3_r8
  REAL(KIND=r8)   , PARAMETER :: XMImin=1.e6_r8*DMImin
  REAL(KIND=r8)   , PARAMETER :: XMImax=1.e6_r8*DMImax
  REAL(KIND=r8)   , PARAMETER :: DelDMI=1.e-6_r8
  INTEGER         , PARAMETER :: MDImin=INT(XMImin)
  INTEGER         , PARAMETER :: MDImax=INT(XMImax)
  !
  !!
  !--- Various ice lookup tables
  !
  REAL(KIND=r8)  :: ACCRI(MDImin:MDImax)
  REAL(KIND=r8)  :: MASSI(MDImin:MDImax)
  REAL(KIND=r8)  :: SDENS(MDImin:MDImax)
  REAL(KIND=r8)  :: VSNOWI(MDImin:MDImax)
  REAL(KIND=r8)  :: VENTI1(MDImin:MDImax)
  REAL(KIND=r8)  :: VENTI2(MDImin:MDImax)

  !--- Mean rain drop diameters varying from 50 microns (0.05 mm) to 450 microns
  !      (0.45 mm), assuming an exponential size distribution.
  !
  REAL(KIND=r8), PRIVATE,PARAMETER :: DMRmin=0.05e-3_r8
  REAL(KIND=r8), PRIVATE,PARAMETER :: DMRmax=0.45e-3_r8
  REAL(KIND=r8), PRIVATE,PARAMETER :: XMRmin=1.e6*DMRmin
  REAL(KIND=r8), PRIVATE,PARAMETER :: XMRmax=1.e6*DMRmax
  REAL(KIND=r8), PRIVATE,PARAMETER :: DelDMR=1.e-6_r8
  REAL(KIND=r8), PRIVATE,PARAMETER :: NLImin=100.0_r8
  REAL(KIND=r8), PRIVATE,PARAMETER :: NLImax=20.E3_r8
  INTEGER      , PRIVATE,PARAMETER :: MDRmin=INT(XMRmin)
  INTEGER      , PRIVATE,PARAMETER :: MDRmax=INT(XMRmax)

  !
  !--- Factor of 1.5 for RECImin, RESNOWmin, & RERAINmin accounts for
  !    integrating exponential distributions for effective radius
  !    (i.e., the r**3/r**2 moments).
  !
  !     INTEGER, PRIVATE, PARAMETER :: INDEXSmin=300
  !!    INTEGER, PRIVATE, PARAMETER :: INDEXSmin=200
  INTEGER      , PRIVATE, PARAMETER :: INDEXSmin=100
  REAL(KIND=r8), PRIVATE, PARAMETER :: RERAINmin=1.5_r8*XMRmin
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=8.0
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=7.5
  REAL(KIND=r8)   , PRIVATE, PARAMETER :: RECImin=1.5_r8*XMImin
  REAL(KIND=r8)   , PRIVATE, PARAMETER :: RESNOWmin=1.5_r8*INDEXSmin
  REAL(KIND=r8)   , PRIVATE, PARAMETER :: RECWmin=10.0_r8

  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=15.
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=5.


  !
  !--- Various rain lookup tables
  !--- Rain lookup tables for mean rain drop diameters from DMRmin to DMRmax,
  !      assuming exponential size distributions for the rain drops
  !
  REAL(KIND=r8), PRIVATE :: ACCRR(MDRmin:MDRmax)
  REAL(KIND=r8), PRIVATE :: MASSR(MDRmin:MDRmax)
  REAL(KIND=r8), PRIVATE :: RRATE(MDRmin:MDRmax)
  REAL(KIND=r8), PRIVATE :: VRAIN(MDRmin:MDRmax)
  REAL(KIND=r8), PRIVATE :: VENTR1(MDRmin:MDRmax)
  REAL(KIND=r8), PRIVATE :: VENTR2(MDRmin:MDRmax)
  !(MDRmin:MDRmax)
  !--- Common block for riming tables
  !--- VEL_RF - velocity increase of rimed particles as functions of crude
  !      particle size categories (at 0.1 mm intervals of mean ice particle
  !      sizes) and rime factor (different values of Rime Factor of 1.1**N,
  !      where N=0 to Nrime).
  !
  INTEGER, PRIVATE,PARAMETER :: Nrime=40
  REAL(KIND=r8)   , PRIVATE :: VEL_RF(2:9,0:Nrime)
  !
  !--- The following variables are for microphysical statistics
  !
  INTEGER, PARAMETER :: ITLO=-60
  INTEGER, PARAMETER :: ITHI=40
  INTEGER  NSTATS(ITLO:ITHI,4)
  REAL(KIND=r8)     QMAX(ITLO:ITHI,5)
  REAL(KIND=r8)     QTOT(ITLO:ITHI,22)
  !
  !    &  T_ICE=-10., T_ICE_init=-5.      !- Ver1
!!!  &, T_ICE=-20.                      !- Ver2
  REAL(KIND=r8), PRIVATE,  PARAMETER :: T_ICE=-40.0_r8
  REAL(KIND=r8), PRIVATE,  PARAMETER :: T_ICE_init=-15.0_r8     !- Ver2
  !    &  T_ICE=-30., T_ICE_init=-5.      !- Ver2

  !     Some other miscellaneous parameters
  !
  REAL(KIND=r8), PRIVATE, PARAMETER :: Thom=T_ICE
  REAL(KIND=r8), PRIVATE, PARAMETER :: TNW=50.0_r8
  REAL(KIND=r8), PRIVATE, PARAMETER :: TOLER=1.0E-20_r8
  !     REAL(KIND=r8), PRIVATE, PARAMETER :: Thom=T_ICE, TNW=50., TOLER=5.E-7
  !     REAL(KIND=r8), PRIVATE, PARAMETER :: Thom=-35., TNW=50., TOLER=5.E-7

  ! Assume fixed cloud ice effective radius
  REAL(KIND=r8), PRIVATE, PARAMETER ::      RECICE=RECImin
  REAL(KIND=r8), PRIVATE, PARAMETER ::      EPSQ=1.0E-20_r8
  !      REAL(KIND=r8), PRIVATE, PARAMETER ::     &, EPSQ=1.E-12
  REAL(KIND=r8), PRIVATE, PARAMETER ::      FLG0P1=0.1_r8
  REAL(KIND=r8), PRIVATE, PARAMETER ::        FLG0P2=0.2_r8
  REAL(KIND=r8), PRIVATE, PARAMETER ::        FLG1P0=1.0_r8


  !
  !------------------------------------------------------------------------- 
  !------- Key parameters, local variables, & important comments ---------
  !-----------------------------------------------------------------------
  !
  !--- KEY Parameters:
  !
  !---- Comments on 14 March 2002
  !    * Set EPSQ to the universal value of 1.e-12 throughout the code
  !      condensate.  The value of EPSQ will need to be changed in the other 
  !      subroutines in order to make it consistent throughout the Eta code.  
  !    * Set CLIMIT=10.*EPSQ as the lower limit for the total mass of 
  !      condensate in the current layer and the input flux of condensate
  !      from above (TOT_ICE, TOT_ICEnew, TOT_RAIN, and TOT_RAINnew).
  !
  !-- NLImax - maximum number concentration of large ice crystals (20,000 /m**3, 20 per liter)
  !-- NLImin - minimum number concentration of large ice crystals (100 /m**3, 0.1 per liter)
  !
  REAL(KIND=r8), PARAMETER ::   RHOL=1000.0_r8
  REAL(KIND=r8), PARAMETER ::   XLS=con_hvap+con_hfus

  !    &, T_ICE=-10.          !- Ver1
  !    &, T_ICE_init=-5.      !- Ver1
!!!  &, T_ICE=-20.          !- Ver2
  !    &, T_ICE=-40.          !- Ver2
  !    &, T_ICE_init=-15.,    !- Ver2
  !
  !    & CLIMIT=10.*EPSQ, EPS1=con_rv/RD-1., RCP=1./CP,

  REAL(KIND=r8), PARAMETER :: RCP=1._r8/con_cp
  REAL(KIND=r8), PARAMETER :: RCPRV=RCP/con_rv
  REAL(KIND=r8), PARAMETER :: RRHOL=1._r8/RHOL
  REAL(KIND=r8), PARAMETER :: XLS1=XLS*RCP
  REAL(KIND=r8), PARAMETER :: XLS2=XLS*XLS*RCPRV
  REAL(KIND=r8), PARAMETER :: XLS3=XLS*XLS/con_rv
  REAL(KIND=r8), PARAMETER :: C1=1._r8/3._r8
  REAL(KIND=r8), PARAMETER :: C2=1._r8/6._r8
  REAL(KIND=r8), PARAMETER :: C3=3.31_r8/6._r8
  REAL(KIND=r8), PARAMETER :: DMR1=.1E-3_r8
  REAL(KIND=r8), PARAMETER :: DMR2=.2E-3_r8
  REAL(KIND=r8), PARAMETER :: DMR3=.32E-3_r8
  REAL(KIND=r8), PARAMETER :: N0r0=8.E6_r8
  REAL(KIND=r8), PARAMETER :: N0rmin=1.e4_r8
  REAL(KIND=r8), PARAMETER :: N0s0=4.E6_r8
  REAL(KIND=r8), PARAMETER :: RHO0=1.194_r8
  REAL(KIND=r8), PARAMETER :: XMR1=1.e6_r8*DMR1
  REAL(KIND=r8), PARAMETER :: XMR2=1.e6_r8*DMR2
  REAL(KIND=r8), PARAMETER :: XMR3=1.e6_r8*DMR3
  REAL(KIND=r8), PARAMETER :: Xratio=.025_r8
  INTEGER      , PARAMETER :: MDR1=INT(XMR1)
  INTEGER      , PARAMETER :: MDR2=INT(XMR2)
  INTEGER      , PARAMETER :: MDR3=INT(XMR3)

  INTEGER  :: llyr  !       llyr : upper limit of boundary layer clouds



  INTEGER       , PARAMETER   :: ntrac=3
  INTEGER,PARAMETER :: nxpvsl=7501
  REAL(KIND=r8)          :: tbpvsl(nxpvsl)
  REAL(KIND=r8)          :: c1xpvsl
  REAL(KIND=r8)          :: c2xpvsl
  INTEGER,PARAMETER:: nxpvsi=7501
  REAL(KIND=r8)          :: tbpvsi(nxpvsi)
  REAL(KIND=r8)          :: c1xpvsi
  REAL(KIND=r8)          :: c2xpvsi
  INTEGER,PARAMETER:: nxpvs=7501
  REAL(KIND=r8) c1xpvs,c2xpvs,tbpvs(nxpvs)
  REAL(KIND=r8), ALLOCATABLE :: si(:)
  REAL(KIND=r8), ALLOCATABLE :: sl(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  INTEGER :: nswbands = 14
  INTEGER :: nbndsw = 14
  ! number of lw bands
  INTEGER :: nlwbands = 16
  INTEGER :: nbndlw = 16
  INTEGER, PARAMETER, PUBLIC :: iulog =0
  REAL   , PUBLIC, PARAMETER :: scalefactor = 1._r8 !500._r8/917._r8

  ! Wavenumbers of band boundaries
  !
  ! Note: Currently rad_solar_var extends the lowest band down to
  ! 100 cm^-1 if it is too high to cover the far-IR. Any changes meant
  ! to affect IR solar variability should take note of this.
  INTEGER ,PARAMETER:: nbndsw_aux = 14
  INTEGER,PARAMETER :: nbndlw_aux = 16

  REAL(R8),PARAMETER :: wavenum_low(nbndsw_aux) = & ! in cm^-1
       (/2600._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
       8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,  820._r8/)
  REAL(R8),PARAMETER :: wavenum_high(nbndsw_aux) = & ! in cm^-1
       (/3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, 8050._r8, &
       12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8, 2600._r8/)

  ! Solar irradiance at 1 A.U. in W/m^2 assumed by radiation code
  ! Rescaled so that sum is precisely 1368.22 and fractional amounts sum to 1.0
  REAL(R8), PARAMETER :: solar_ref_band_irradiance(nbndsw_aux) = & 
       (/ &
       12.11_r8,  20.3600000000001_r8, 23.73_r8, &
       22.43_r8,  55.63_r8, 102.93_r8, 24.29_r8, &
       345.74_r8, 218.19_r8, 347.20_r8, &
       129.49_r8,  50.15_r8,   3.08_r8, 12.89_r8 &
       /)

  INTEGER , PARAMETER :: n_g_d = 300
  ! d_eff:long_name = "effective diameter" ;
  ! d_eff:units = "microns" ;
  REAL(r8), ALLOCATABLE :: g_d_eff(:)      !;d_eff :long_name = "effective diameter"
  !;d_eff:units = "microns" ;
  REAL(r8), ALLOCATABLE :: ext_sw_ice(:,:) !;lw_abs:long_name = "Longwave mass specific absorption for in-cloud ice water path" 
  !;lw_abs:units = "meter^2 kilogram^-1" ;
  REAL(r8), ALLOCATABLE :: ssa_sw_ice(:,:) !;sw_ext:long_name = "Shortwave extinction"
  !;sw_ext:units = "meter^2 kilogram^-1" ;
  REAL(r8), ALLOCATABLE :: asm_sw_ice(:,:) !;sw_ssa:long_name = "Shortwave single scattering albedo"
  !;sw_ssa:units = "fraction" ;
  REAL(r8), ALLOCATABLE :: abs_lw_ice(:,:) !;sw_asm:long_name = "Shortwave asymmetry parameter"
  !;sw_asm:units = "fraction" ;

  INTEGER, PARAMETER :: nmu = 20
  INTEGER, PARAMETER :: nlambda = 50

  REAL(r8), ALLOCATABLE :: g_mu(:)           ! mu:units       = "unitless" ;
  REAL(r8), ALLOCATABLE :: g_lambda(:,:)     ! lambda:units   = "meter^-1" ;
  REAL(r8), ALLOCATABLE :: ext_sw_liq(:,:,:) ! k_ext_sw:units = "meters^2/kg" ;
  REAL(r8), ALLOCATABLE :: ssa_sw_liq(:,:,:) ! ssa_sw  :units = "0 to 1 unitless" ;
  REAL(r8), ALLOCATABLE :: asm_sw_liq(:,:,:) ! asm_sw:units   = "-1 to 1 unitless" ;
  REAL(r8), ALLOCATABLE :: abs_lw_liq(:,:,:) ! k_abs_lw:units = "meters^2/kg" ;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Init_Optical_Properties
  PUBLIC :: RunCloudOpticalProperty
  PUBLIC :: RunCloudOpticalProperty2
  PUBLIC :: Cloud_Micro_WRF
CONTAINS

  SUBROUTINE Init_Optical_Properties(dtpg,NLAY,si_in,sl_in,nbndsw_in,nbndlw_in,fNameCldOptSW , &
             fNameCldOptLW )
    IMPLICIT NONE
    REAL(KIND=r8)   , INTENT(in) :: dtpg
    INTEGER         , INTENT(in) :: NLAY          !   NLAY            : vertical layer number                             !
    REAL(KIND=r8)   , INTENT(in) :: si_in(NLAY+1)
    REAL(KIND=r8)   , INTENT(in) :: sl_in(NLAY)
    INTEGER         , INTENT(in) :: nbndsw_in
    INTEGER         , INTENT(in) :: nbndlw_in

    CHARACTER(LEN=*), INTENT(IN) :: fNameCldOptSW 
    CHARACTER(LEN=*), INTENT(IN) :: fNameCldOptLW 
    nbndsw   = nbndsw_in
    nbndlw   = nbndlw_in
    nswbands = nbndsw_in
    nlwbands = nbndlw_in
    !...................................
    CALL cldinit(dtpg,NLAY,si_in,sl_in)
    !-----------------------------------
    !
    !--- Parameters & data statement for local calculations
    !
    CALL cloud_rad_props_init(fNameCldOptSW ,fNameCldOptLW )

  END SUBROUTINE Init_Optical_Properties


  !-----------------------------------
  SUBROUTINE RunCloudOpticalProperty2(       &
       IM            ,&
       LM            ,&
       ILCON         ,&
       delsig        ,&
       imask         ,&
       colrad        ,&
       ps            ,&!mb
       tgrs          ,&
       qgrs          ,&
       FlipPbot      ,&
       fcice         ,&
       fcliq         ,&
       frain         ,&
       tskn          ,&
       tsea          ,&
       cld           ,&
       clu           ,&
       EFFCS         ,&
       EFFIS         ,&
       clwp          ,&
       lmixr         ,&
       fice          ,&
       rei           ,&
       rel           ,&
       taud          ,&
       cicewp        ,&
       cliqwp        ,&
       c_cld_tau     ,&
       c_cld_tau_w   ,&
       c_cld_tau_w_g ,&
       c_cld_tau_w_f ,&
       c_cld_lw_abs  ,&
       cldfprime      )

    INTEGER         , INTENT(IN   ) :: IM
    INTEGER         , INTENT(IN   ) :: LM
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    REAL(KIND=r8)   , INTENT(IN   ) :: delsig(LM) 
    INTEGER(KIND=i8), INTENT(IN   ) :: imask (IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad(IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: ps    (IM)!mb
    REAL(KIND=r8)   , INTENT(IN   ) :: tgrs  (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: qgrs  (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: tskn  (IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: tsea  (IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: FlipPbot     (IM,LM)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8)   , INTENT(IN   ) :: fcice (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: fcliq (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: frain (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: cld   (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: clu   (IM,LM) 

    REAL(KIND=r8)   , INTENT(INOUT) :: EFFCS (IM,LM) 
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFIS (IM,LM) 
    ! Cloud properties
    REAL(KIND=r8), INTENT(OUT) :: clwp (IM,LM) ! Cloud Liquid Water Path
    REAL(KIND=r8), INTENT(OUT) :: lmixr(IM,LM) ! Ice/Water mixing ratio
    REAL(KIND=r8), INTENT(OUT) :: fice (IM,LM) ! Fractional amount of cloud that is ice
    REAL(KIND=r8), INTENT(OUT) :: rei  (IM,LM) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: rel  (IM,LM) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: taud (IM,LM) ! Shortwave cloud optical depth
    REAL(KIND=r8), INTENT(OUT) :: cicewp (IM,LM)
    REAL(KIND=r8), INTENT(OUT) :: cliqwp (IM,LM)
    ! combined cloud radiative parameters are "in cloud" not "in cell"
    REAL(r8)     , INTENT(OUT) :: c_cld_tau    (1:nbndsw,1:IM,1:LM) ! cloud extinction optical depth
    REAL(r8)     , INTENT(OUT) :: c_cld_tau_w  (1:nbndsw,1:IM,1:LM) ! cloud single scattering albedo * tau
    REAL(r8)     , INTENT(OUT) :: c_cld_tau_w_g(1:nbndsw,1:IM,1:LM) ! cloud assymetry parameter * w * tau
    REAL(r8)     , INTENT(OUT) :: c_cld_tau_w_f(1:nbndsw,1:IM,1:LM) ! cloud forward scattered fraction * w * tau
    REAL(r8)     , INTENT(OUT) :: c_cld_lw_abs (1:nbndlw,1:IM,1:LM) ! cloud absorption optics depth (LW)
    REAL(r8)     , INTENT(OUT) :: cldfprime    (1:IM,1:LM)          ! combined cloud fraction (snow plus regular)


    INTEGER            :: np3d   !=4  zhao/moorthi's prognostic cloud scheme
    !=3  ferrier's microphysics
    INTEGER, PARAMETER :: ncld=1 !      ncld            : only used when ntcw .gt. 0                     !
    INTEGER, PARAMETER :: NTRAC=1
    INTEGER, PARAMETER :: ntcw=1
    !LOGICAL, PARAMETER :: norad_precip    !      norad_precip    : logical flag for not using precip in radiation !
    LOGICAL, PARAMETER :: norad_precip = .FALSE.   ! This is effective only for Ferrier/Moorthi
    LOGICAL, PARAMETER :: crick_proof  = .TRUE.
    LOGICAL, PARAMETER :: ccnorm       = .FALSE.
    LOGICAL, PARAMETER :: sashal       = .TRUE.
    INTEGER, PARAMETER :: iovrsw =0 !      iovrsw/iovrlw   : control flag for cloud overlap (sw/lw rad)     !
    !                        =0 random overlapping clouds                   !
    !                        =1 max/ran overlapping clouds                  !
    INTEGER, PARAMETER ::iflip= 1!      iflip           : control flag for in/out vertical indexing      !
    !                        =0 index from toa to surface                   !
    !                        =1 index from surface to toa                   !
    REAL (kind=r8)  :: flgmin(IM)!      flgmin          : minimim large ice fraction                     !
    REAL (kind=r8)  :: clw(IM,LM) 
    REAL (kind=r8)  :: clouds(IM,LM,NF_CLDS)

    !     2. cloud profiles:      (defined in 'module_radiation_clouds')    !
    !                ---  for  prognostic cloud  ---                        !
    !          clouds(:,:,1)  -  layer total cloud fraction                 !
    !          clouds(:,:,2)  -  layer cloud liq water path                 !
    !          clouds(:,:,3)  -  mean effective radius for liquid cloud     !
    !          clouds(:,:,4)  -  layer cloud ice water path                 !
    !          clouds(:,:,5)  -  mean effective radius for ice cloud        !
    !          clouds(:,:,6)  -  layer rain drop water path                 !
    !          clouds(:,:,7)  -  mean effective radius for rain drop        !
    !          clouds(:,:,8)  -  layer snow flake water path                !
    !          clouds(:,:,9)  -  mean effective radius for snow flake       !
    !                ---  for  diagnostic cloud  ---                        !
    !          clouds(:,:,1)  -  layer total cloud fraction                 !
    !          clouds(:,:,2)  -  layer cloud optical depth                  !
    !          clouds(:,:,3)  -  layer cloud single scattering albedo       !
    !          clouds(:,:,4)  -  layer cloud asymmetry factor               !
    !                                                                       !
    REAL (kind=r8)   :: cldsa (IM,5)
    INTEGER          :: mbota (IM,3)
    INTEGER          :: mtopa (IM,3)
    REAL (kind=r8)   :: tlyr  (IM,LM)
    REAL (kind=r8)   :: rhly  (IM,LM)
    REAL (kind=r8)   :: qstl  (IM,LM)
    REAL (kind=r8)   :: qlyr  (IM,LM)
    REAL (kind=r8)   :: plyr  (IM,LM)
    REAL(KIND=r8)    :: prsl  (IM,LM)    !      prsi  (IX,LM+1) : model level pressure in cb      (kPa) !
    REAL(KIND=r8)    :: prsi  (IM,LM+1)  !      prsl  (IX,LM)   : model layer mean pressure in cb (kPa)          !
    REAL(KIND=r8)    :: slmsk (IM)       !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
    REAL(KIND=r8)    :: rrime (IM,LM) 
    REAL(KIND=r8)    :: cldfsnow (IM,LM) 
    REAL(KIND=r8)    :: lamc     (IM,LM) 
    REAL(KIND=r8)    :: pgam    (IM,LM) 
    LOGICAL          :: dosw
    LOGICAL          :: dolw
    LOGICAL          :: oldcldoptics
    CHARACTER(LEN=200) :: liqcldoptics
    CHARACTER(LEN=200) :: icecldoptics


    REAL (kind=r8)   :: xlat  (IM)!xlat (IM)  : grid longitude/latitude in radians             !
    REAL (kind=r8) :: plvl(IM,LM+1)
    REAL(KIND=r8) :: Zibot  (IM,LM+1) ! Height at middle of layer (m)
    REAL(KIND=r8) :: emziohl(IM,LM+1) ! exponential of Minus zi Over hl (no dim)
    REAL(KIND=r8) :: pdel   (IM,LM) ! Moist pressure difference across layer Pressure thickness [Pa] > 0
!    REAL(KIND=r8) :: gicewp (IM,LM)! grid-box cloud ice water path
!    REAL(KIND=r8) :: gliqwp (IM,LM)! grid-box cloud liquid water path
!    REAL(KIND=r8) :: cicewp (IM,LM)! in-cloud cloud ice water path
!    REAL(KIND=r8) :: cliqwp (IM,LM)! in-cloud cloud liquid water path
!    REAL(KIND=r8) :: cwp    (IM,LM)! in-cloud cloud (total) water path

    REAL(KIND=r8) :: hl     (IM)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: rhl    (IM)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: pw     (IM)        ! precipitable water (kg/m2)
    REAL(KIND=r8) :: kabs                   ! longwave absorption coeff (m**2/g)
    REAL(KIND=r8) :: kabsi                  ! ice absorption coefficient
    REAL(KIND=r8) :: emis(IM,LM)       ! cloud emissivity (fraction)
    REAL(KIND=r8) :: ocnfrac     (1:IM)   ! Ocean fraction
    REAL(KIND=r8) :: rrlrv
    REAL(KIND=r8) :: const

    !  ---  outputs:
    REAL(KIND=r8)    :: es
    REAL(KIND=r8)    :: qs
    INTEGER :: k
    INTEGER :: i,j
    INTEGER :: LP1
    REAL(KIND=r8), PARAMETER :: clwc0   = 0.21_r8 ! Reference liquid water concentration (g/m3)        
    REAL(KIND=r8), PARAMETER :: clwc0_Emirical = 1.0_r8       ! Model-top presure                                 
    REAL(KIND=r8), PARAMETER :: kabsl=0.090361_r8               ! longwave liquid absorption coeff (m**2/g)
    ! --- abssnow is the snow flake absorption coefficient (micron)
    REAL (kind=r8), PARAMETER :: abssnow0=1.5_r8             ! fu   coeff

    ! --- absrain is the rain drop absorption coefficient (m2/g)
    !     data absrain / 3.07e-3 /          ! chou coeff
    REAL (kind=r8), PARAMETER :: absrain=0.33e-3_r8          ! ncar coeff

    ! --- abssnow is the snow flake absorption coefficient (m2/g)
    !      data kabsl / 2.34e-3 /         ! ncar coeff

    !     absliqn is the liquid water absorption coefficient (m2/g).
    ! === for iflagliq = 1,
    !      data kabsl  absliq1 / 0.0602410 /
    REAL (kind=r8), PARAMETER ::  absliq1=0.0602410_r8
    !  ---  constant values
    REAL (kind=r8), PARAMETER  :: QMIN=1.0e-10_r8
    REAL (kind=r8), PARAMETER  :: QME5=1.0e-7_r8
    REAL (kind=r8), PARAMETER  :: QME6=1.0e-7_r8
    REAL (kind=r8), PARAMETER  :: EPSQ=1.0e-12_r8
    INTEGER                    :: cldfsnow_idx
    clwp=0.0_r8;lmixr=0.0_r8;fice=0.0_r8;rei=0.0_r8;rel=0.0_r8;taud=0.0_r8
    flgmin=0.0_r8;clw=0.0_r8;clouds=0.0_r8;cldsa=0.0_r8;mbota=-1;mtopa=-1;
    tlyr=0.0_r8
    cicewp=0.0_r8
    cliqwp=0.0_r8
    c_cld_tau=0.0_r8
    c_cld_tau_w=0.0_r8
    c_cld_tau_w_g=0.0_r8
    c_cld_tau_w_f=0.0_r8
    c_cld_lw_abs=0.0_r8

    IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
       np3d=4 !=4  zhao/moorthi's prognostic cloud scheme
       !=3  ferrier's microphysics
    ELSE IF ( TRIM(ILCON) == 'MIC'.OR. TRIM(ILCON).EQ.'HWRF' .OR. TRIM(ILCON).EQ.'HGFS'.OR.&
         TRIM(ILCON).EQ.'UKMO' .OR. TRIM(ILCON).EQ.'MORR'.OR.TRIM(ILCON).EQ.'HUMO') THEN
       np3d=3 !=4  zhao/moorthi's prognostic cloud scheme
       !=3  ferrier's microphysics
    END IF

    ! taud(i,k) = absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)

    !if (crick_proof) print *,' CRICK-Proof cloud water used in radiation '
    !if (ccnorm) print *,' Cloud condensate normalized by cloud cover for radiation'
    !if (sashal) print *,' New Massflux based shallow convection used'
    flgmin(:)        = 0.20_r8
    pw=0.0_r8
    !
    !===> ...  begin here
    !
    LP1 = LM + 1
    !  --- ...  prepare atmospheric profiles for radiation input
    !           convert pressure unit from mb to cb

    DO i=1,IM
       IF(iMask(i) > 0_i8)THEN
          ! land
          slmsk (i) = 1.0_r8  !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
          ocnfrac(i)=0.0_r8
       ELSE
          ! water/ocean
          slmsk (i) = 0.0_r8  !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
          ocnfrac  (i)=1.0_r8
          IF(ocnfrac(i).GT.0.01_r8.AND.ABS(tsea(i)).LT.260.0_r8) THEN
             slmsk (i) = 2.0_r8  !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
             ocnfrac(i) = 1.0_r8
          ENDIF
       END IF
    END DO

    DO i = 1, IM
       prsi(i,LM + 1)=MAX(ps(i)*si(LM + 1)/10.0_r8,1.0e-12_r8) !mb  -- > cb
    END DO

    DO k = 1, LM
       DO i = 1, IM
          rrime(i,k)= 1.0_r8
          prsl(i,k)=ps(i)*si(k)/10.0_r8
          prsi(i,k)=ps(i)*si(k)/10.0_r8
       ENDDO
    ENDDO

    !      tgrs(i,k) 
    !      qgrs(i,k)
    !      tskn(i) 
    !
    !  --- ...  prepare atmospheric profiles for radiation input
    !           convert pressure unit from cb to mb
    DO k = 1, LM
       DO i = 1, IM
          es  = MIN( prsl(i,k), 0.001_r8 * fpvs( tgrs(i,k) ) )   ! fpvs in pa
          qs  = MAX( QMIN, con_eps * es / (prsl(i,k) + con_epsm1*es) )
          rhly(i,k) = MAX( 0.0_r8, MIN( 1.0_r8, MAX(QMIN, qgrs(i,k))/qs ) )
          qstl(i,k) = qs
       ENDDO
    ENDDO

    DO k = 1, LM
       DO i = 1, IM
          qlyr(i,k) = MAX( QME6    , qgrs(i,k) )
       ENDDO
    ENDDO
    !  --- ...  prepare atmospheric profiles for radiation input
    !           convert pressure unit from cb to mb


    DO k = 1, LM
       DO i = 1, IM
          plvl(i,k) = 10.0_r8 * prsi(i,k)!cb -- >mb
          plyr(i,k) = 10.0_r8 * prsl(i,k)!cb -- >mb
          tlyr(i,k) = tgrs(i,k)
       ENDDO
    ENDDO

    DO i = 1, IM
       xlat=(((colrad(i)))-(3.1415926e0_r8/2.0_r8))
       plvl(i,LM + 1) = 10.0_r8 * prsi(i,LM + 1)!cb -- >mb
    ENDDO

    ! Heights corresponding to sigma at middle of layer: sig(k)
    ! Assuming isothermal atmosphere within each layer
    DO i=1,IM
       Zibot(i,1) = 0.0_r8
       DO k=2,LM
          Zibot(i,k) = Zibot(i,k-1) + (con_rd/con_g)*tgrs(i,k-1)* &
               !               LOG(sigbot(k-1)/sigbot(k))
               LOG(FlipPbot(i,LM+2-k)/FlipPbot(i,LM+1-k))
       END DO
    END DO

    DO i=1,IM
       Zibot(i,LM+1)=Zibot(i,LM)+(con_rd/con_g)*tgrs(i,LM)* &
            LOG(FlipPbot(i,1)/pptop)
    END DO

    ! precitable water, pw = sum_k { delsig(k) . Qe(k) } . Ps . 100 / g
    !                   pw = sum_k { Dp(k) . Qe(k) } / g
    !
    ! 100 is to change from mbar to pascal
    ! Dp(k) is the difference of pressure (N/m2) between bottom and top of layer
    ! Qe(k) is specific humidity in (g/g)
    ! gravity is m/s2 => so pw is in Kg/m2
    DO k=1,LM
       DO i = 1,IM
          pw(i) = pw(i) + delsig(k)*(qgrs(i,k))
       END DO
    END DO
    DO i = 1,IM
       pw(i)=100.0_r8*pw(i)*ps(i)/con_g
    END DO
    !
    ! diagnose liquid water scale height from precipitable water
    DO i=1,IM
       hl(i)  = 700.0_r8*LOG(MAX(pw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,LM
       DO i=1,IM
           pdel(i,k)   = (prsi(i,k)-prsi(i,k+1))*1000.0_r8! cb -- > Pa
       END DO
    END DO

    DO k=1,LM+1
       DO i=1,IM
          !          emziohl(i,k) = EXP(-zibot(i,k)/hl(i))
          emziohl(i,k) = EXP(-zibot(i,k)*rhl(i))
       END DO
    END DO
    !    DO i=1,ncols
    !       emziohl(i,kmax+1) = 0.0_r8
    !    END DO
    !
    ! Liquid water path is a mesure of total amount of liquid water present
    ! between two points  int he atmosphere
    ! Typical values of liquid water path in marine stratocumulus can be 
    ! of the order of 20-80 [g/m*m].
    !        --Po
    !       \    ql
    ! clwp = \ ----- * dp
    !        /   g
    !       /
    !        --0
    ! The units are g/m2.
    DO k=1,LM
       DO i=1,IM
          clwp(i,k) = clwc0_Emirical*clwc0*hl(i)*(emziohl(i,k) - emziohl(i,k+1))
       END DO
    END DO
    ! If we want to calculate the 'droplets/cristals' mixing ratio, we need
    ! to find the amount of dry air in each layer. 
    !
    !             dry_air_path = int rho_air dz  
    !
    ! This can be simply done using the hydrostatic equation:
    !
    !              dp/dz = -rho grav
    !              dp/grav = -rho dz
    !
    !
    ! The units are g/m2. The factor 1e5 accounts for the change
    !  mbar to Pa and kg/m2 to g/m2. 
    DO k=1,LM
       DO i=1,IM
          lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/Ps(i)
       END DO
    END DO



    IF (ntcw > 0) THEN                   ! prognostic cloud scheme

       DO k = 1, LM
          DO i = 1, IM
             clw(i,k) = 0.0_r8
          END DO
          IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
             DO j = 1, ncld
                DO i = 1, IM
                   clw(i,k) =  fcice(i,k) + fcliq(i,k) + frain(i,k)   ! cloud condensate amount
                ENDDO
             ENDDO
          ELSE IF ( TRIM(ILCON) == 'MIC'.OR. TRIM(ILCON).EQ.'HWRF' .OR. TRIM(ILCON).EQ.'HGFS'.OR.&
               TRIM(ILCON).EQ.'UKMO' .OR. TRIM(ILCON).EQ.'MORR'.OR.TRIM(ILCON).EQ.'HUMO') THEN
             DO j = 1, ncld
                DO i = 1, IM
                   clw(i,k) =  fcice(i,k) + fcliq(i,k) + frain(i,k)   ! cloud condensate amount
                ENDDO
             ENDDO
          END IF

       ENDDO

       WHERE (clw < EPSQ)
          clw = 0.0_r8
       END WHERE

       IF (np3d == 4) THEN              ! zhao/moorthi's prognostic cloud scheme

          CALL progcld1 (&
                                !  ---  inputs:
               plyr   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: plyr (:,:) : model layer mean pressure in mb (100Pa)  
               plvl   (1:IM,1:LM+1)   , &!real (kind=r8), intent(in) :: plvl  (:,:): model level pressure in mb (100Pa)       
               tlyr   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: tlyr (:,:): model layer mean temperature in k           
               qlyr   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: qlyr (:,:): layer specific humidity in gm/gm               
               qstl   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: qstl (:,:): layer saturate humidity in gm/gm                      !
               rhly   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: rhly (:,:): layer relative humidity (=qlyr/qstl) 
               clw    (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: clw  (:,:): layer cloud condensate amount      
               xlat   (1:IM)          , &!real (kind=r8), intent(in) :: xlat(:): grid latitude in radians                
               slmsk  (1:IM)          , &!real (kind=r8), intent(in) :: slmsk(:): sea/land mask array (sea:0,land:1,sea-ice:2)
               IM                     , &!integer,  intent(in) :: IX : horizontal dimention                     
               LM                     , &!integer,  intent(in) :: NLAY : vertical layer/level dimensions           
               iflip                  , &!integer,  intent(in) :: iflip: control flag for in/out vertical indexing 
               iovrsw                 , &!integer,  intent(in) :: iovr : control flag for cloud overlap            
               sashal                 , &!logical, intent(in) :: sashal
               crick_proof            , &!logical, intent(in) :: crick_proof
               ccnorm                 , &!logical, intent(in) :: ccnorm
                                !  ---  outputs:
               clouds (1:IM,1:LM,1:NF_CLDS)     , &!real (kind=r8),intent(out) :: clouds(:,:,:) : cloud profiles        
               cldsa  (1:IM,1:5)      , &!real (kind=r8),   intent(out) :: clds(:,:): fraction of clouds for low, mid, hi, tot, bl        !
               mtopa  (1:IM,1:3)      , &!integer,,   intent(out) :: mtop(:,:) : vertical indices for low, mid, hi cloud tops         !
               mbota  (1:IM,1:3)        )!integer,,   intent(out) :: mbot(:,:) : vertical indices for low, mid, hi cloud bases         !

       ELSEIF (np3d == 3) THEN          ! ferrier's microphysics

          !     print *,' in RunCloudOpticalProperty : calling progcld2'
          CALL progcld2 ( &
                                !  ---  inputs:
               plyr    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: plyr(:,:)
               plvl    (1:IM,1:LM+1)   ,&!REAL (kind=r8), INTENT(in) :: plvl(:,:)
               tlyr    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: tlyr(:,:)
               qlyr    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: qlyr(:,:)
               qstl    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: qstl(:,:)
               rhly    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: rhly(:,:)
               clw     (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: clw(:,:)
               xlat    (1:IM)          ,&!REAL (kind=r8), INTENT(in) :: xlat(:)
               fcice   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: f_ice(:,:)
               frain   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: f_rain(:,:)
               fcliq   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: f_rain(:,:)
               rrime   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: r_rime(:,:)
               flgmin  (1:IM)          ,&!REAL (kind=r8), INTENT(in) :: flgmin(:)
               IM                      ,&!INTEGER,  INTENT(in) :: IX
               LM                      ,&!INTEGER,  INTENT(in) :: NLAY
               iflip                   ,&!INTEGER,  INTENT(in) :: iflip
               iovrsw                  ,&!INTEGER,  INTENT(in) :: iovr
               sashal                  ,&!LOGICAL       , INTENT(in) :: sashal
               norad_precip            ,&!LOGICAL       , INTENT(in) :: norad_precip
               crick_proof             ,&!LOGICAL       , INTENT(in) :: crick_proof
               ccnorm                  ,&!LOGICAL       , INTENT(in) :: ccnorm
                                !  ---  outputs   
               clouds   (1:IM,1:LM,1:NF_CLDS)    ,&!REAL (kind=r8), INTENT(out) :: clouds(:,:,:)
               cldsa    (1:IM,1:5)    ,&!REAL (kind=r8), INTENT(out) :: clds(:,:)
               mtopa    (1:IM,1:3)    ,&!INTEGER      , INTENT(out) :: mtop(:,:)
               mbota    (1:IM,1:3)     )!INTEGER      , INTENT(out) :: mbot(:,:)

       ENDIF                            ! end if_np3d
    END IF
    ! define fractional amount of cloud that is ice
    ! if warmer than -10 degrees c then water phase
    ! docs CCM3, eq 4.a.16.1     
    ! allcld_liq = state%q(:,:,ixcldliq)
    ! allcld_ice = state%q(:,:,ixcldice)
    IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
       DO k=1,LM
          DO i=1,IM
             fice(i,k)=MAX(MIN((263.16_r8-tgrs(i,k))*0.05_r8,1.0_r8),0.0_r8)
          END DO
       END DO
    ELSE IF ( TRIM(ILCON) == 'MIC'.OR. TRIM(ILCON).EQ.'HWRF' .OR. TRIM(ILCON).EQ.'HGFS'.OR.&
         TRIM(ILCON).EQ.'UKMO' .OR. TRIM(ILCON).EQ.'MORR'.OR.TRIM(ILCON).EQ.'HUMO') THEN
       DO k=1,LM
          DO i=1,IM
             !fice(i,k) = MIN(MAX( fcice(i,k) /max(1.e-10_r8,(fcice(i,k) + fcliq(i,k) + frain(i,k))),0.000_r8),1.0_r8)
             fice(i,k) = MIN(MAX( fcice(i,k) /MAX(1.e-10_r8,(fcice(i,k) + fcliq(i,k) )),0.000_r8),1.0_r8)
          END DO
       END DO
    END IF
 
   DO k=1,LM
      DO i=1,IM
         clouds(i,k,1) =MAX(cld(i,LM-k+1),clu(i,LM-k+1))
         clouds(i,k,1) =MAX(clouds(i,k,1) , 2.0e-80_r8)
      END DO
   END DO

    cldfsnow_idx=0
    cldfsnow=0.0_r8
    lamc=0.0_r8
    pgam=0.0_r8
    oldcldoptics=.FALSE.
    liqcldoptics='slingo'
    icecldoptics='ebertcurry'!  ('mitchell')
    dosw=.TRUE.  
    dolw=.TRUE. 
    CALL Run_Optical_Properties( &
       IM                        , &
       IM                        , &
       LM                        , &
       cldfsnow_idx              , &
       cldfsnow    (1:IM,1:LM)   , &
       clouds      (1:IM,1:LM,1) , &
       fcliq       (1:IM,1:LM)   , &
       fcice       (1:IM,1:LM)   , &
       clouds      (1:IM,1:LM,3) , &!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
       clouds      (1:IM,1:LM,5) , &
       clouds      (1:IM,1:LM,2) , &
       clouds      (1:IM,1:LM,4) , &
       clwp        (1:IM,1:LM)   , &
       pdel        (1:IM,1:LM)   , &
       lamc        (1:IM,1:LM)   , &
       pgam        (1:IM,1:LM)   , &
    ! combined cloud radiative parameters are "in cloud" not "in cell"
       c_cld_tau    (1:nbndsw,1:IM,1:LM), & ! cloud extinction optical depth
       c_cld_tau_w  (1:nbndsw,1:IM,1:LM), & ! cloud single scattering albedo * tau
       c_cld_tau_w_g(1:nbndsw,1:IM,1:LM), & ! cloud assymetry parameter * w * tau
       c_cld_tau_w_f(1:nbndsw,1:IM,1:LM), & ! cloud forward scattered fraction * w * tau
       c_cld_lw_abs (1:nbndlw,1:IM,1:LM) ,& ! cloud absorption optics depth (LW)
       cldfprime    (1:IM,1:LM)  , &
       dosw                      , &
       dolw                      , &
       oldcldoptics              , &
       icecldoptics              , &
       liqcldoptics              )

    ! Compute optical depth from liquid water

    DO k=1,LM
       DO i=1,IM

          !      clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
          !      clouds(:,:,1) - layer total cloud fraction                       !

          !      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
          cliqwp(i,k)=clouds(i,k,2) 
          
          !IF ( TRIM(ILCON).EQ.'HUMO') THEN
          !   rel(i,k)=EFFCS (i,k)  ! EFFCS -  DROPLET EFFECTIVE RADIUS   (MICRON)
          !   rei(i,k)=EFFIS (i,k)  ! EFFIS -  CLOUD ICE EFFECTIVE RADIUS (MICRON)
          !ELSE
          rel(i,k)=clouds(i,k,3)!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
          !      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
          cicewp(i,k)=clouds(i,k,4)
          rei(i,k)=clouds(i,k,5)!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
          !END IF
          !      clouds(:,:,6) - layer rain drop water path         not assigned  !
          !      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !

          !  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
          !      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !

          !  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !

          !   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !

          !   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !

          !   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !


          !WRITE(*,'(A5,3e12.5)'),'CAM5',clwp(i,k),&
          !                  clouds(i,k,2)+clouds(i,k,4)+clouds(i,k,6)+clouds(i,k,8),& 
          !        clwp(i,k)+clouds(i,k,2)+clouds(i,k,4)+clouds(i,k,6)+clouds(i,k,8)
          IF ( TRIM(ILCON) .EQ. 'MORR' .OR.TRIM(ILCON).EQ.'HUMO') THEN

             !pdel(i,k)   = (prsi(i,k)-prsi(i,k+1))*1000.0_r8! cb -- > Pa
             !               g/g       Pa       /m/s2
             !
             !        N        kg * m         kg
             !  Pa =-----  =  ----------- = --------
             !        m*m       m*m**s2       m s**2
             !
             !             kg     s**2        kg
             !  P*g= = -------- * -----  = --------
             !           m s**2    m*m       m*m
             !
             !                     g
             !  P*g *1000 =    --------
             !                     m*m
             !
             !gicewp(i,k) = fcice(i,k)*pdel(i,k)/con_g*1000.0_r8                     ! Grid box ice water path.[g/m2]
             !gliqwp(i,k) = fcliq(i,k)*pdel(i,k)/con_g*1000.0_r8                     ! Grid box liquid water path.[g/m2]

             !cicewp(i,k) = gicewp(i,k) / max(0.01_r8,clouds(i,k,1))                 ! In-cloud ice water path.[g/m2]
             !cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,clouds(i,k,1))                 ! In-cloud liquid water path.[g/m2]

             !cwp   (i,k) = cicewp(i,k) + cliqwp(i,k) ! in-cloud cloud (total) water path[g/m2]


             ! lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/ps(i)! ! Ice/Water mixing ratio
             !lmixr(i,k)=fcice(i,k)/max(fcliq(i,k),1.0e-12_r8) ! ! Ice/Water mixing ratio

             !note that optical properties for ice valid only
             !in range of 13 > rei > 130 micron (Ebert and Curry 92)
             !if ( TRIM(ILCON) .eq. 'MORR' ) then

             !kabsi = 0.005_r8 + 1.0_r8/min(max(13.0_r8,rei(i,k)),130.0_r8)
             !else if ( microp_scheme .eq. 'RK' ) then
             !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
             !END IF
             !     (m**2/g)
             !kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 

             ! cloud emissivity (fraction)
             !emis(i,k) = 1.0_r8 - exp(-1.66_r8*kabs*cwp(i,k))  ! In-cloud liquid water path.
             ! cloud optical depth
             !taud (i,k) = kabs*cwp(i,k)                        ! In-cloud liquid water path.
             !taud(i,k) =   absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)

             clwp(i,k)=MAX(MIN(MAX(clwp(i,k), clouds(i,k,2)),500.0_r8),1.0e-8_r8)

             lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/ps(i) ! ! Ice/Water mixing ratio

             !note that optical properties for ice valid only
             !in range of 13 > rei > 130 micron (Ebert and Curry 92)
             !if ( TRIM(ILCON) .eq. 'MORR' ) then

             kabsi = 0.005_r8 + 1.0_r8/MIN(MAX(13.0_r8,rei(i,k)),130.0_r8)
             !else if ( microp_scheme .eq. 'RK' ) then
             !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
             !END IF
             !     (m**2/g)
             kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 
             ! cloud emissivity (fraction)
             emis(i,k) = 1.0_r8 - EXP(-1.66_r8*kabs*clwp(i,k))
             ! cloud optical depth
             !taud(i,k) = kabs*clwp(i,k)! g/m2
             !taud(i,k) =  MAX(MIN(MAX(kabs*clwp(i,k), absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)),500.0_r8),1.0e-8_r8)
             !WRITE(*,'(6E12.4)')SUM(c_cld_tau(1:nbndsw,i,k)),SUM(c_cld_tau_w (1:nbndsw,i,k)),SUM(c_cld_tau_w_g(1:nbndsw,i,k)),SUM(c_cld_tau_w_f(1:nbndsw,i,k)), &
             !SUM(c_cld_lw_abs (1:nbndlw,i,k)),taud(i,k)
             !taud(i,k) = SUM(c_cld_lw_abs (1:nbndlw,i,k))
             taud(i,k) =  MAX(MIN(MAX(kabs*clwp(i,k), absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)),500.0_r8),1.0e-8_r8)

          ELSE
             !clwp(i,k)=MAX(MIN((clwp(i,k),clouds(i,k,2)),500.0_r8),1.0e-8_r8)
             clwp(i,k)=MAX(MIN(MAX(clwp(i,k), clouds(i,k,2)),500.0_r8),1.0e-8_r8)

             lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/ps(i) ! ! Ice/Water mixing ratio

             !note that optical properties for ice valid only
             !in range of 13 > rei > 130 micron (Ebert and Curry 92)
             !if ( TRIM(ILCON) .eq. 'MORR' ) then

             kabsi = 0.005_r8 + 1.0_r8/MIN(MAX(13.0_r8,rei(i,k)),130.0_r8)
             !else if ( microp_scheme .eq. 'RK' ) then
             !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
             !END IF
             !     (m**2/g)
             kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 
             ! cloud emissivity (fraction)
             emis(i,k) = 1.0_r8 - EXP(-1.66_r8*kabs*clwp(i,k))
             ! cloud optical depth
             !taud(i,k) = kabs*clwp(i,k)! g/m2
             !taud(i,k) = MAX(MIN( kabs*clwp(i,k) + absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8),500.0_r8),1.0e-8_r8)

             !WRITE(*,'(6E12.4)')SUM(c_cld_tau(1:nbndsw,i,k)),SUM(c_cld_tau_w (1:nbndsw,i,k)),SUM(c_cld_tau_w_g(1:nbndsw,i,k)),SUM(c_cld_tau_w_f(1:nbndsw,i,k)), &
             !SUM(c_cld_lw_abs (1:nbndlw,i,k)),taud(i,k)
             !taud(i,k) = SUM(c_cld_lw_abs (1:nbndlw,i,k))
             taud(i,k) = MAX(MIN( kabs*clwp(i,k) + absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8),500.0_r8),1.0e-8_r8)
             !taud(i,k) =  MAX(MIN(MAX(kabs*clwp(i,k), absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)),500.0_r8),1.0e-8_r8)

          END IF
       END DO
    END DO


  END SUBROUTINE RunCloudOpticalProperty2


  !-----------------------------------
  SUBROUTINE RunCloudOpticalProperty(       &
       IM    ,&
       LM    ,&
       ILCON ,&
       delsig,&
       imask ,&
       colrad,&
       ps    ,&!mb
       tgrs  ,&
       qgrs  ,&
       FlipPbot ,&
       fcice ,&
       fcliq ,&
       frain ,&
       tskn  ,&
       tsea  ,&
       EFFCS ,&
       EFFIS ,&
       clwp  ,&
       lmixr ,&
       fice  ,&
       rei   ,&
       rel   ,&
       taud   )

    INTEGER         , INTENT(IN   ) :: IM
    INTEGER         , INTENT(IN   ) :: LM
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    REAL(KIND=r8)   , INTENT(IN   ) :: delsig(LM) 
    INTEGER(KIND=i8), INTENT(IN   ) :: imask (IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad(IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: ps    (IM)!mb
    REAL(KIND=r8)   , INTENT(IN   ) :: tgrs  (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: qgrs  (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: tskn  (IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: tsea  (IM)
    REAL(KIND=r8)   , INTENT(IN   ) :: FlipPbot     (IM,LM)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8)   , INTENT(IN   ) :: fcice (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: fcliq (IM,LM) 
    REAL(KIND=r8)   , INTENT(IN   ) :: frain (IM,LM) 
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFCS (IM,LM) 
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFIS (IM,LM) 
    ! Cloud properties
    REAL(KIND=r8), INTENT(OUT) :: clwp (IM,LM) ! Cloud Liquid Water Path
    REAL(KIND=r8), INTENT(OUT) :: lmixr(IM,LM) ! Ice/Water mixing ratio
    REAL(KIND=r8), INTENT(OUT) :: fice (IM,LM) ! Fractional amount of cloud that is ice
    REAL(KIND=r8), INTENT(OUT) :: rei  (IM,LM) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: rel  (IM,LM) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: taud (IM,LM) ! Shortwave cloud optical depth


    INTEGER            :: np3d   !=4  zhao/moorthi's prognostic cloud scheme
    !=3  ferrier's microphysics
    INTEGER, PARAMETER :: ncld=1 !      ncld            : only used when ntcw .gt. 0                     !
    INTEGER, PARAMETER :: NTRAC=1
    INTEGER, PARAMETER :: ntcw=1
    !LOGICAL, PARAMETER :: norad_precip    !      norad_precip    : logical flag for not using precip in radiation !
    LOGICAL, PARAMETER :: norad_precip = .FALSE.   ! This is effective only for Ferrier/Moorthi
    LOGICAL, PARAMETER :: crick_proof  = .TRUE.
    LOGICAL, PARAMETER :: ccnorm       = .FALSE.
    LOGICAL, PARAMETER :: sashal       = .TRUE.
    INTEGER, PARAMETER :: iovrsw =0 !      iovrsw/iovrlw   : control flag for cloud overlap (sw/lw rad)     !
    !                        =0 random overlapping clouds                   !
    !                        =1 max/ran overlapping clouds                  !
    INTEGER, PARAMETER ::iflip= 1!      iflip           : control flag for in/out vertical indexing      !
    !                        =0 index from toa to surface                   !
    !                        =1 index from surface to toa                   !
    REAL (kind=r8)  :: flgmin(IM)!      flgmin          : minimim large ice fraction                     !
    REAL (kind=r8)  :: clw(IM,LM) 
    REAL (kind=r8)  :: clouds(IM,LM,NF_CLDS)

    !     2. cloud profiles:      (defined in 'module_radiation_clouds')    !
    !                ---  for  prognostic cloud  ---                        !
    !          clouds(:,:,1)  -  layer total cloud fraction                 !
    !          clouds(:,:,2)  -  layer cloud liq water path                 !
    !          clouds(:,:,3)  -  mean effective radius for liquid cloud     !
    !          clouds(:,:,4)  -  layer cloud ice water path                 !
    !          clouds(:,:,5)  -  mean effective radius for ice cloud        !
    !          clouds(:,:,6)  -  layer rain drop water path                 !
    !          clouds(:,:,7)  -  mean effective radius for rain drop        !
    !          clouds(:,:,8)  -  layer snow flake water path                !
    !          clouds(:,:,9)  -  mean effective radius for snow flake       !
    !                ---  for  diagnostic cloud  ---                        !
    !          clouds(:,:,1)  -  layer total cloud fraction                 !
    !          clouds(:,:,2)  -  layer cloud optical depth                  !
    !          clouds(:,:,3)  -  layer cloud single scattering albedo       !
    !          clouds(:,:,4)  -  layer cloud asymmetry factor               !
    !                                                                       !
    REAL (kind=r8)   :: cldsa (IM,5)
    INTEGER          :: mbota (IM,3)
    INTEGER          :: mtopa (IM,3)
    REAL (kind=r8)   :: tlyr  (IM,LM)
    REAL (kind=r8)   :: rhly  (IM,LM)
    REAL (kind=r8)   :: qstl  (IM,LM)
    REAL (kind=r8)   :: qlyr  (IM,LM)
    REAL (kind=r8)   :: plyr  (IM,LM)
    REAL(KIND=r8)    :: prsl  (IM,LM)    !      prsi  (IX,LM+1) : model level pressure in cb      (kPa) !
    REAL(KIND=r8)    :: prsi  (IM,LM+1)  !      prsl  (IX,LM)   : model layer mean pressure in cb (kPa)          !
    REAL(KIND=r8)    :: slmsk (IM)       !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
    REAL(KIND=r8)    :: rrime (IM,LM) 


    REAL (kind=r8)   :: xlat  (IM)!xlat (IM)  : grid longitude/latitude in radians             !
    REAL (kind=r8) :: plvl(IM,LM+1)
    REAL(KIND=r8) :: Zibot  (IM,LM+1) ! Height at middle of layer (m)
    REAL(KIND=r8) :: emziohl(IM,LM+1) ! exponential of Minus zi Over hl (no dim)
    REAL(KIND=r8) :: pdel   (IM,LM) ! Moist pressure difference across layer Pressure thickness [Pa] > 0
!    REAL(KIND=r8) :: gicewp (IM,LM)! grid-box cloud ice water path
!    REAL(KIND=r8) :: gliqwp (IM,LM)! grid-box cloud liquid water path
!    REAL(KIND=r8) :: cicewp (IM,LM)! in-cloud cloud ice water path
!    REAL(KIND=r8) :: cliqwp (IM,LM)! in-cloud cloud liquid water path
!    REAL(KIND=r8) :: cwp    (IM,LM)! in-cloud cloud (total) water path

    REAL(KIND=r8) :: hl     (IM)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: rhl    (IM)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: pw     (IM)        ! precipitable water (kg/m2)
    REAL(KIND=r8) :: kabs                   ! longwave absorption coeff (m**2/g)
    REAL(KIND=r8) :: kabsi                  ! ice absorption coefficient
    REAL(KIND=r8) :: emis(IM,LM)       ! cloud emissivity (fraction)
    REAL(KIND=r8) :: ocnfrac     (1:IM)   ! Ocean fraction
    REAL(KIND=r8) ::  rrlrv
    REAL(KIND=r8) ::  const

    !  ---  outputs:
    REAL(KIND=r8)    :: es
    REAL(KIND=r8)    :: qs
    INTEGER :: k
    INTEGER :: i,j
    INTEGER :: LP1
    REAL(KIND=r8), PARAMETER :: pptop = 0.5       ! Model-top presure                                 
    REAL(KIND=r8), PARAMETER :: clwc0   = 0.21_r8 ! Reference liquid water concentration (g/m3)        
    REAL(KIND=r8), PARAMETER :: clwc0_Emirical = 1.0_r8       ! Model-top presure                                 

    REAL(KIND=r8), PARAMETER :: kabsl=0.090361_r8               ! longwave liquid absorption coeff (m**2/g)
    ! --- abssnow is the snow flake absorption coefficient (micron)
    REAL (kind=r8), PARAMETER :: abssnow0=1.5_r8             ! fu   coeff

    ! --- absrain is the rain drop absorption coefficient (m2/g)
    !     data absrain / 3.07e-3 /          ! chou coeff
    REAL (kind=r8), PARAMETER :: absrain=0.33e-3_r8          ! ncar coeff

    ! --- abssnow is the snow flake absorption coefficient (m2/g)
    !      data kabsl / 2.34e-3 /         ! ncar coeff

    !     absliqn is the liquid water absorption coefficient (m2/g).
    ! === for iflagliq = 1,
    !      data kabsl  absliq1 / 0.0602410 /
    REAL (kind=r8), PARAMETER ::  absliq1=0.0602410_r8
    !  ---  constant values
    REAL (kind=r8), PARAMETER  :: QMIN=1.0e-10_r8
    REAL (kind=r8), PARAMETER  :: QME5=1.0e-7_r8
    REAL (kind=r8), PARAMETER  :: QME6=1.0e-7_r8
    REAL (kind=r8), PARAMETER  :: EPSQ=1.0e-12_r8
    clwp=0.0_r8;lmixr=0.0_r8;fice=0.0_r8;rei=0.0_r8;rel=0.0_r8;taud=0.0_r8
    flgmin=0.0_r8;clw=0.0_r8;clouds=0.0_r8;cldsa=0.0_r8;mbota=-1;mtopa=-1;
    tlyr=0.0_r8
    IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
       np3d=4 !=4  zhao/moorthi's prognostic cloud scheme
       !=3  ferrier's microphysics
    ELSE IF ( TRIM(ILCON) == 'MIC'.OR. TRIM(ILCON).EQ.'HWRF' .OR. TRIM(ILCON).EQ.'HGFS'.OR.&
         TRIM(ILCON).EQ.'UKMO' .OR. TRIM(ILCON).EQ.'MORR'.OR.TRIM(ILCON).EQ.'HUMO') THEN
       np3d=3 !=4  zhao/moorthi's prognostic cloud scheme
       !=3  ferrier's microphysics
    END IF

    ! taud(i,k) = absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)

    !if (crick_proof) print *,' CRICK-Proof cloud water used in radiation '
    !if (ccnorm) print *,' Cloud condensate normalized by cloud cover for radiation'
    !if (sashal) print *,' New Massflux based shallow convection used'
    flgmin(:)        = 0.20_r8
    pw=0.0_r8
    !
    !===> ...  begin here
    !
    LP1 = LM + 1
    !  --- ...  prepare atmospheric profiles for radiation input
    !           convert pressure unit from mb to cb

    DO i=1,IM
       IF(iMask(i) > 0_i8)THEN
          ! land
          slmsk (i) = 1.0_r8  !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
          ocnfrac(i)=0.0_r8
       ELSE
          ! water/ocean
          slmsk (i) = 0.0_r8  !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
          ocnfrac  (i)=1.0_r8
          IF(ocnfrac(i).GT.0.01_r8.AND.ABS(tsea(i)).LT.260.0_r8) THEN
             slmsk (i) = 2.0_r8  !      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
             ocnfrac(i) = 1.0_r8
          ENDIF
       END IF
    END DO

    DO i = 1, IM
       prsi(i,LM + 1)=MAX(ps(i)*si(LM + 1)/10.0_r8,1.0e-12_r8) !mb  -- > cb
    END DO

    DO k = 1, LM
       DO i = 1, IM
          rrime(i,k)= 1.0_r8
          prsl(i,k)=ps(i)*si(k)/10.0_r8
          prsi(i,k)=ps(i)*si(k)/10.0_r8
       ENDDO
    ENDDO

    !      tgrs(i,k) 
    !      qgrs(i,k)
    !      tskn(i) 
    !
    !  --- ...  prepare atmospheric profiles for radiation input
    !           convert pressure unit from cb to mb
    !
    ! Find saturated moisture
    !
    rrlrv = -con_hvap/(rmwmdi*con_rd)
    const = e0c*EXP(-rrlrv/tmelt)

    DO k = 1, LM
       DO i = 1, IM
          !IF(tgrs(i,k)>270.0_r8)THEN
          !   es = const*EXP(rrlrv/tgrs(i,k))*0.1!mb convert to cb
          !   qs  = MAX( QMIN, con_eps * es /MAX(( (prsl(i,k)         + con_epsm1*es)    ),1.0e-12_r8))
          !ELSE
             es = 0.001_r8 * fpvs2es5(tgrs(i,k)) !Pa ->cb
             qs  = MAX( QMIN, con_eps * es /MAX(( (prsl(i,k)         + con_epsm1*es)    ),1.0e-12_r8))
          !END IF
          !es  = MIN( prsl(i,k), 0.001_r8 * fpvs( tgrs(i,k) ) )   ! fpvs in pa convert to cb
          !qs  = MAX( QMIN, con_eps * es / (prsl(i,k) + con_epsm1*es) )
          rhly(i,k) = MAX( 0.0_r8, MIN( 1.0_r8, MAX(QMIN, qgrs(i,k))/qs ) )
          qstl(i,k) = qs
       ENDDO
    ENDDO

    DO k = 1, LM
       DO i = 1, IM
          qlyr(i,k) = MAX( QME6    , qgrs(i,k) )
       ENDDO
    ENDDO
    !  --- ...  prepare atmospheric profiles for radiation input
    !           convert pressure unit from cb to mb


    DO k = 1, LM
       DO i = 1, IM
          plvl(i,k) = 10.0_r8 * prsi(i,k)!cb -- >mb
          plyr(i,k) = 10.0_r8 * prsl(i,k)!cb -- >mb
          tlyr(i,k) = tgrs(i,k)
       ENDDO
    ENDDO

    DO i = 1, IM
       xlat=(((colrad(i)))-(3.1415926e0_r8/2.0_r8))
       plvl(i,LM + 1) = 10.0_r8 * prsi(i,LM + 1)!cb -- >mb
    ENDDO

    ! Heights corresponding to sigma at middle of layer: sig(k)
    ! Assuming isothermal atmosphere within each layer
    DO i=1,IM
       Zibot(i,1) = 0.0_r8
       DO k=2,LM
          Zibot(i,k) = Zibot(i,k-1) + (con_rd/con_g)*tgrs(i,k-1)* &
               !               LOG(sigbot(k-1)/sigbot(k))
               LOG(FlipPbot(i,LM+2-k)/FlipPbot(i,LM+1-k))
       END DO
    END DO

    DO i=1,IM
       Zibot(i,LM+1)=Zibot(i,LM)+(con_rd/con_g)*tgrs(i,LM)* &
            LOG(FlipPbot(i,1)/pptop)
    END DO

    ! precitable water, pw = sum_k { delsig(k) . Qe(k) } . Ps . 100 / g
    !                   pw = sum_k { Dp(k) . Qe(k) } / g
    !
    ! 100 is to change from mbar to pascal
    ! Dp(k) is the difference of pressure (N/m2) between bottom and top of layer
    ! Qe(k) is specific humidity in (g/g)
    ! gravity is m/s2 => so pw is in Kg/m2
    DO k=1,LM
       DO i = 1,IM
          pw(i) = pw(i) + delsig(k)*(qgrs(i,k))
       END DO
    END DO
    DO i = 1,IM
       pw(i)=100.0_r8*pw(i)*ps(i)/con_g
    END DO
    !
    ! diagnose liquid water scale height from precipitable water
    DO i=1,IM
       hl(i)  = 700.0_r8*LOG(MAX(pw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,LM+1
       DO i=1,IM
          !          emziohl(i,k) = EXP(-zibot(i,k)/hl(i))
          emziohl(i,k) = EXP(-zibot(i,k)*rhl(i))
       END DO
    END DO
    !    DO i=1,ncols
    !       emziohl(i,kmax+1) = 0.0_r8
    !    END DO
    !
    ! Liquid water path is a mesure of total amount of liquid water present
    ! between two points  int he atmosphere
    ! Typical values of liquid water path in marine stratocumulus can be 
    ! of the order of 20-80 [g/m*m].
    !        --Po
    !       \    ql
    ! clwp = \ ----- * dp
    !        /   g
    !       /
    !        --0
    ! The units are g/m2.
    DO k=1,LM
       DO i=1,IM
          clwp(i,k) = clwc0_Emirical*clwc0*hl(i)*(emziohl(i,k) - emziohl(i,k+1))
       END DO
    END DO

    ! If we want to calculate the 'droplets/cristals' mixing ratio, we need
    ! to find the amount of dry air in each layer. 
    !
    !             dry_air_path = int rho_air dz  
    !
    ! This can be simply done using the hydrostatic equation:
    !
    !              dp/dz = -rho grav
    !              dp/grav = -rho dz
    !
    !
    ! The units are g/m2. The factor 1e5 accounts for the change
    !  mbar to Pa and kg/m2 to g/m2. 
    DO k=1,LM
       DO i=1,IM
          lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/Ps(i)
       END DO
    END DO



    IF (ntcw > 0) THEN                   ! prognostic cloud scheme

       DO k = 1, LM
          DO i = 1, IM
             clw(i,k) = 0.0_r8
          END DO
          IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
             DO j = 1, ncld
                DO i = 1, IM
                   clw(i,k) =  fcice(i,k) + fcliq(i,k) + frain(i,k)   ! cloud condensate amount
                ENDDO
             ENDDO
          ELSE IF ( TRIM(ILCON) == 'MIC'.OR. TRIM(ILCON).EQ.'HWRF' .OR. TRIM(ILCON).EQ.'HGFS'.OR.&
               TRIM(ILCON).EQ.'UKMO' .OR. TRIM(ILCON).EQ.'MORR'.OR.TRIM(ILCON).EQ.'HUMO') THEN
             DO j = 1, ncld
                DO i = 1, IM
                   clw(i,k) =  fcice(i,k) + fcliq(i,k) + frain(i,k)   ! cloud condensate amount
                ENDDO
             ENDDO
          END IF

       ENDDO

       WHERE (clw < EPSQ)
          clw = 0.0_r8
       END WHERE

       IF (np3d == 4) THEN              ! zhao/moorthi's prognostic cloud scheme

          CALL progcld1 (&
                                !  ---  inputs:
               plyr   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: plyr (:,:) : model layer mean pressure in mb (100Pa)  
               plvl   (1:IM,1:LM+1)   , &!real (kind=r8), intent(in) :: plvl  (:,:): model level pressure in mb (100Pa)       
               tlyr   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: tlyr (:,:): model layer mean temperature in k           
               qlyr   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: qlyr (:,:): layer specific humidity in gm/gm               
               qstl   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: qstl (:,:): layer saturate humidity in gm/gm                      !
               rhly   (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: rhly (:,:): layer relative humidity (=qlyr/qstl) 
               clw    (1:IM,1:LM)     , &!real (kind=r8), intent(in) :: clw  (:,:): layer cloud condensate amount      
               xlat   (1:IM)          , &!real (kind=r8), intent(in) :: xlat(:): grid latitude in radians                
               slmsk  (1:IM)          , &!real (kind=r8), intent(in) :: slmsk(:): sea/land mask array (sea:0,land:1,sea-ice:2)
               IM                     , &!integer,  intent(in) :: IX : horizontal dimention                     
               LM                     , &!integer,  intent(in) :: NLAY : vertical layer/level dimensions           
               iflip                  , &!integer,  intent(in) :: iflip: control flag for in/out vertical indexing 
               iovrsw                 , &!integer,  intent(in) :: iovr : control flag for cloud overlap            
               sashal                 , &!logical, intent(in) :: sashal
               crick_proof            , &!logical, intent(in) :: crick_proof
               ccnorm                 , &!logical, intent(in) :: ccnorm
                                !  ---  outputs:
               clouds (1:IM,1:LM,1:NF_CLDS)     , &!real (kind=r8),intent(out) :: clouds(:,:,:) : cloud profiles        
               cldsa  (1:IM,1:5)      , &!real (kind=r8),   intent(out) :: clds(:,:): fraction of clouds for low, mid, hi, tot, bl        !
               mtopa  (1:IM,1:3)      , &!integer,,   intent(out) :: mtop(:,:) : vertical indices for low, mid, hi cloud tops         !
               mbota  (1:IM,1:3)        )!integer,,   intent(out) :: mbot(:,:) : vertical indices for low, mid, hi cloud bases         !

       ELSEIF (np3d == 3) THEN          ! ferrier's microphysics

          !     print *,' in RunCloudOpticalProperty : calling progcld2'
          CALL progcld2 ( &
                                !  ---  inputs:
               plyr    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: plyr(:,:)
               plvl    (1:IM,1:LM+1)   ,&!REAL (kind=r8), INTENT(in) :: plvl(:,:)
               tlyr    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: tlyr(:,:)
               qlyr    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: qlyr(:,:)
               qstl    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: qstl(:,:)
               rhly    (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: rhly(:,:)
               clw     (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: clw(:,:)
               xlat    (1:IM)          ,&!REAL (kind=r8), INTENT(in) :: xlat(:)
               fcice   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: f_ice(:,:)
               frain   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: f_rain(:,:)
               fcliq   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: f_rain(:,:)
               rrime   (1:IM,1:LM)     ,&!REAL (kind=r8), INTENT(in) :: r_rime(:,:)
               flgmin  (1:IM)          ,&!REAL (kind=r8), INTENT(in) :: flgmin(:)
               IM                      ,&!INTEGER,  INTENT(in) :: IX
               LM                      ,&!INTEGER,  INTENT(in) :: NLAY
               iflip                   ,&!INTEGER,  INTENT(in) :: iflip
               iovrsw                  ,&!INTEGER,  INTENT(in) :: iovr
               sashal                  ,&!LOGICAL       , INTENT(in) :: sashal
               norad_precip            ,&!LOGICAL       , INTENT(in) :: norad_precip
               crick_proof             ,&!LOGICAL       , INTENT(in) :: crick_proof
               ccnorm                  ,&!LOGICAL       , INTENT(in) :: ccnorm
                                !  ---  outputs   
               clouds   (1:IM,1:LM,1:NF_CLDS)    ,&!REAL (kind=r8), INTENT(out) :: clouds(:,:,:)
               cldsa    (1:IM,1:5)    ,&!REAL (kind=r8), INTENT(out) :: clds(:,:)
               mtopa    (1:IM,1:3)    ,&!INTEGER      , INTENT(out) :: mtop(:,:)
               mbota    (1:IM,1:3)     )!INTEGER      , INTENT(out) :: mbot(:,:)

       ENDIF                            ! end if_np3d
    END IF
    ! define fractional amount of cloud that is ice
    ! if warmer than -10 degrees c then water phase
    ! docs CCM3, eq 4.a.16.1     
    ! allcld_liq = state%q(:,:,ixcldliq)
    ! allcld_ice = state%q(:,:,ixcldice)
    IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
       DO k=1,LM
          DO i=1,IM
             fice(i,k)=MAX(MIN((263.16_r8-tgrs(i,k))*0.05_r8,1.0_r8),0.0_r8)
          END DO
       END DO
    ELSE IF ( TRIM(ILCON) == 'MIC'.OR. TRIM(ILCON).EQ.'HWRF' .OR. TRIM(ILCON).EQ.'HGFS'.OR.&
         TRIM(ILCON).EQ.'UKMO' .OR. TRIM(ILCON).EQ.'MORR'.OR.TRIM(ILCON).EQ.'HUMO') THEN
       DO k=1,LM
          DO i=1,IM
             !fice(i,k) = MIN(MAX( fcice(i,k) /max(1.e-10_r8,(fcice(i,k) + fcliq(i,k) + frain(i,k))),0.000_r8),1.0_r8)
             fice(i,k) = MIN(MAX( fcice(i,k) /MAX(1.e-10_r8,(fcice(i,k) + fcliq(i,k) )),0.000_r8),1.0_r8)
          END DO
       END DO
    END IF


    ! Compute optical depth from liquid water

    DO k=1,LM
       DO i=1,IM

          !      clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
          !      clouds(:,:,1) - layer total cloud fraction                       !

          !      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
          !IF ( TRIM(ILCON).EQ.'HUMO') THEN
          !   rel(i,k)=EFFCS (i,k)  ! EFFCS -  DROPLET EFFECTIVE RADIUS   (MICRON)
          !   rei(i,k)=EFFIS (i,k)  ! EFFIS -  CLOUD ICE EFFECTIVE RADIUS (MICRON)
          !ELSE
          rel(i,k)=clouds(i,k,3)!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
          !      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
          rei(i,k)=clouds(i,k,5)!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
          !END IF
          !      clouds(:,:,6) - layer rain drop water path         not assigned  !
          !      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !

          !  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
          !      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !

          !  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !

          !   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !

          !   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !

          !   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !

          !WRITE(*,'(A5,3e12.5)'),'CAM5',clwp(i,k),&
          !                  clouds(i,k,2)+clouds(i,k,4)+clouds(i,k,6)+clouds(i,k,8),& 
          !        clwp(i,k)+clouds(i,k,2)+clouds(i,k,4)+clouds(i,k,6)+clouds(i,k,8)
          IF ( TRIM(ILCON) .EQ. 'MORR' .OR.TRIM(ILCON).EQ.'HUMO') THEN

             !pdel(i,k)   = (prsi(i,k)-prsi(i,k+1))*1000.0_r8! cb -- > Pa
             !               g/g       Pa       /m/s2
             !
             !        N        kg * m         kg
             !  Pa =-----  =  ----------- = --------
             !        m*m       m*m**s2       m s**2
             !
             !             kg     s**2        kg
             !  P*g= = -------- * -----  = --------
             !           m s**2    m*m       m*m
             !
             !                     g
             !  P*g *1000 =    --------
             !                     m*m
             !
             !gicewp(i,k) = fcice(i,k)*pdel(i,k)/con_g*1000.0_r8                     ! Grid box ice water path.[g/m2]
             !gliqwp(i,k) = fcliq(i,k)*pdel(i,k)/con_g*1000.0_r8                     ! Grid box liquid water path.[g/m2]

             !cicewp(i,k) = gicewp(i,k) / max(0.01_r8,clouds(i,k,1))                 ! In-cloud ice water path.[g/m2]
             !cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,clouds(i,k,1))                 ! In-cloud liquid water path.[g/m2]

             !cwp   (i,k) = cicewp(i,k) + cliqwp(i,k) ! in-cloud cloud (total) water path[g/m2]


             ! lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/ps(i)! ! Ice/Water mixing ratio
             !lmixr(i,k)=fcice(i,k)/max(fcliq(i,k),1.0e-12_r8) ! ! Ice/Water mixing ratio

             !note that optical properties for ice valid only
             !in range of 13 > rei > 130 micron (Ebert and Curry 92)
             !if ( TRIM(ILCON) .eq. 'MORR' ) then

             !kabsi = 0.005_r8 + 1.0_r8/min(max(13.0_r8,rei(i,k)),130.0_r8)
             !else if ( microp_scheme .eq. 'RK' ) then
             !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
             !END IF
             !     (m**2/g)
             !kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 

             ! cloud emissivity (fraction)
             !emis(i,k) = 1.0_r8 - exp(-1.66_r8*kabs*cwp(i,k))  ! In-cloud liquid water path.
             ! cloud optical depth
             !taud (i,k) = kabs*cwp(i,k)                        ! In-cloud liquid water path.
             !taud(i,k) =   absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)

             clwp(i,k)=MAX(MIN(MAX(clwp(i,k), clouds(i,k,2)),500.0_r8),1.0e-8_r8)

             lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/ps(i) ! ! Ice/Water mixing ratio

             !note that optical properties for ice valid only
             !in range of 13 > rei > 130 micron (Ebert and Curry 92)
             !if ( TRIM(ILCON) .eq. 'MORR' ) then

             kabsi = 0.005_r8 + 1.0_r8/MIN(MAX(13.0_r8,rei(i,k)),130.0_r8)
             !else if ( microp_scheme .eq. 'RK' ) then
             !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
             !END IF
             !     (m**2/g)
             kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 
             ! cloud emissivity (fraction)
             emis(i,k) = 1.0_r8 - EXP(-1.66_r8*kabs*clwp(i,k))
             ! cloud optical depth
             !taud(i,k) = kabs*clwp(i,k)! g/m2
             taud(i,k) =  MAX(MIN(MAX(kabs*clwp(i,k), absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)),500.0_r8),1.0e-8_r8)

          ELSE
             !clwp(i,k)=MAX(MIN((clwp(i,k),clouds(i,k,2)),500.0_r8),1.0e-8_r8)
             clwp(i,k)=MAX(MIN(MAX(clwp(i,k), clouds(i,k,2)),500.0_r8),1.0e-8_r8)

             lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/ps(i) ! ! Ice/Water mixing ratio

             !note that optical properties for ice valid only
             !in range of 13 > rei > 130 micron (Ebert and Curry 92)
             !if ( TRIM(ILCON) .eq. 'MORR' ) then

             kabsi = 0.005_r8 + 1.0_r8/MIN(MAX(13.0_r8,rei(i,k)),130.0_r8)
             !else if ( microp_scheme .eq. 'RK' ) then
             !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
             !END IF
             !     (m**2/g)
             kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 
             ! cloud emissivity (fraction)
             emis(i,k) = 1.0_r8 - EXP(-1.66_r8*kabs*clwp(i,k))
             ! cloud optical depth
             !taud(i,k) = kabs*clwp(i,k)! g/m2
             ! output variables:                                                     !
             !   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
             !      clouds(:,:,1) - layer total cloud fraction                       !
             !      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
             !      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
             !      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
             !      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
             !      clouds(:,:,6) - layer rain drop water path         (g/m**2)      !
             !      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
             !  *** clouds(:,:,8) - layer snow flake water path        (g/m**2)      !
             !      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
             taud(i,k) = MAX(MIN(MAX(kabs*clwp(i,k) , absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8)),500.0_r8),1.0e-8_r8)
             !taud(i,k) = MAX(MIN(    kabs*clwp(i,k) + absliq1* clouds(i,k,2) + absrain * clouds(i,k,6) + abssnow0 * clouds(i,k,8),500.0_r8),1.0e-8_r8)
          END IF
       END DO
    END DO


  END SUBROUTINE RunCloudOpticalProperty



  SUBROUTINE Run_Optical_Properties(ncol,pcols,pver,cldfsnow_idx,cldfsnow,cldn,ql,qi,rel,rei,&
       cliqwp,iciwp,clwp,pdel,lamc,pgam,c_cld_tau ,c_cld_tau_w , c_cld_tau_w_g,c_cld_tau_w_f, &
       c_cld_lw_abs,cldfprime , dosw,dolw,&
       oldcldoptics,icecldoptics,liqcldoptics)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: ncol
    INTEGER      , INTENT(IN   ) :: pcols
    INTEGER      , INTENT(IN   ) :: pver
    INTEGER      , INTENT(IN   ) :: cldfsnow_idx
    REAL(KIND=r8), INTENT(IN   ) :: cldfsnow (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: cldn     (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: ql       (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: qi       (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: rel      (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: rei      (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: cliqwp    (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: iciwp    (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: clwp     (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: pdel     (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: lamc     (pcols,pver)
    REAL(KIND=r8), INTENT(IN   ) :: pgam     (pcols,pver)
    LOGICAL      , INTENT(IN   ) :: dosw
    LOGICAL      , INTENT(IN   ) :: dolw
    LOGICAL      , INTENT(IN   ) :: oldcldoptics
    CHARACTER(LEN=*), INTENT(IN   ) :: liqcldoptics
    CHARACTER(LEN=*), INTENT(IN   ) :: icecldoptics
    ! combined cloud radiative parameters are "in cloud" not "in cell"
    REAL(r8)     , INTENT(OUT   ) :: c_cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
    REAL(r8)     , INTENT(OUT   ) :: c_cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
    REAL(r8)     , INTENT(OUT   ) :: c_cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
    REAL(r8)     , INTENT(OUT   ) :: c_cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
    REAL(r8)     , INTENT(OUT   ) :: c_cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)
    REAL(r8)     , INTENT(OUT   ) :: cldfprime    (pcols,pver)             ! combined cloud fraction (snow plus regular)
    REAL(KIND=r8) :: iclwp     (pcols,pver)
    REAL(KIND=r8) :: dei      (pcols,pver)

    ! cloud radiative parameters are "in cloud" not "in cell"
    REAL(r8) :: cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
    REAL(r8) :: cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
    REAL(r8) :: cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
    REAL(r8) :: cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
    REAL(r8) :: cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    REAL(r8) :: ice_tau    (nbndsw,pcols,pver) ! ice extinction optical depth
    REAL(r8) :: ice_tau_w  (nbndsw,pcols,pver) ! ice single scattering albedo * tau
    REAL(r8) :: ice_tau_w_g(nbndsw,pcols,pver) ! ice assymetry parameter * tau * w
    REAL(r8) :: ice_tau_w_f(nbndsw,pcols,pver) ! ice forward scattered fraction * tau * w
    REAL(r8) :: ice_lw_abs (nbndlw,pcols,pver) ! ice absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    REAL(r8) :: snow_tau    (nbndsw,pcols,pver) ! snow extinction optical depth
    REAL(r8) :: snow_tau_w  (nbndsw,pcols,pver) ! snow single scattering albedo * tau
    REAL(r8) :: snow_tau_w_g(nbndsw,pcols,pver) ! snow assymetry parameter * tau * w
    REAL(r8) :: snow_tau_w_f(nbndsw,pcols,pver) ! snow forward scattered fraction * tau * w
    REAL(r8) :: snow_lw_abs (nbndlw,pcols,pver)   ! snow absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    REAL(r8) :: liq_tau    (nbndsw,pcols,pver) ! liquid extinction optical depth
    REAL(r8) :: liq_tau_w  (nbndsw,pcols,pver) ! liquid single scattering albedo * tau
    REAL(r8) :: liq_tau_w_g(nbndsw,pcols,pver) ! liquid assymetry parameter * tau * w
    REAL(r8) :: liq_tau_w_f(nbndsw,pcols,pver) ! liquid forward scattered fraction * tau * w
    REAL(r8) :: liq_lw_abs (nbndlw,pcols,pver) ! liquid absorption optics depth (LW)

    INTEGER  :: i
    INTEGER  :: k
    c_cld_tau     =0.0_r8
    c_cld_tau_w   =0.0_r8
    c_cld_tau_w_g =0.0_r8
    c_cld_tau_w_f =0.0_r8
    c_cld_lw_abs  =0.0_r8
    cldfprime     =0.0_r8

    cld_tau     =0.0_r8
    cld_tau_w   =0.0_r8
    cld_tau_w_g =0.0_r8
    cld_tau_w_f =0.0_r8
    cld_lw_abs  =0.0_r8
    ice_tau     =0.0_r8
    ice_tau_w   =0.0_r8
    ice_tau_w_g =0.0_r8
    ice_tau_w_f =0.0_r8
    ice_lw_abs  =0.0_r8
    snow_tau    =0.0_r8
    snow_tau_w  =0.0_r8
    snow_tau_w_g=0.0_r8
    snow_tau_w_f=0.0_r8
    snow_lw_abs =0.0_r8
    liq_tau     =0.0_r8
    liq_tau_w   =0.0_r8
    liq_tau_w_g =0.0_r8
    liq_tau_w_f =0.0_r8
    liq_lw_abs  =0.0_r8
    dei=rei * 2._r8
    DO k=1,pver
       DO i=1,pcols
          iclwp(i,k)=MAX(MIN(MAX(clwp(i,k), cliqwp(i,k)),500.0_r8),1.0e-8_r8)
       END DO
    END DO
    
    IF (dosw) THEN
       IF(oldcldoptics) THEN
          CALL ec_ice_optics_sw(ncol,pcols,pver,cldn,rei,iciwp,qi,pdel,&
               ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.FALSE.)
          CALL slingo_liq_optics_sw(ncol,pcols,pver,cldn,rel,iclwp,ql,pdel,&
               liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.FALSE.)

       ELSE
          SELECT CASE (TRIM(icecldoptics))
          CASE ('ebertcurry')
             CALL  ec_ice_optics_sw(ncol,pcols,pver,cldn,rei,iciwp,qi,pdel,&
                  ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.TRUE.)
          CASE ('mitchell')
             CALL get_ice_optics_sw(ncol,pcols,pver,dei,iciwp,&
                  ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)
          CASE default
             !call endrun('iccldoptics must be one either ebertcurry or mitchell')
             STOP 'iccldoptics must be one either ebertcurry or mitchell'
          END SELECT

          SELECT CASE (TRIM(liqcldoptics))
          CASE ('slingo')
             CALL slingo_liq_optics_sw(ncol,pcols,pver,cldn,rel,iclwp,ql,pdel,&
                  liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.TRUE.)
          CASE ('gammadist')
             CALL get_liquid_optics_sw(ncol,pcols,pver,lamc,pgam,iclwp, &
                  liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
          CASE default
             !call endrun('liqcldoptics must be either slingo or gammadist')
             STOP 'liqcldoptics must be either slingo or gammadist'
          END SELECT
       ENDIF
       cld_tau    (:,1:ncol,:) =  liq_tau    (:,1:ncol,:) + ice_tau    (:,1:ncol,:)
       cld_tau_w  (:,1:ncol,:) =  liq_tau_w  (:,1:ncol,:) + ice_tau_w  (:,1:ncol,:)
       cld_tau_w_g(:,1:ncol,:) =  liq_tau_w_g(:,1:ncol,:) + ice_tau_w_g(:,1:ncol,:)
       cld_tau_w_f(:,1:ncol,:) =  liq_tau_w_f(:,1:ncol,:) + ice_tau_w_f(:,1:ncol,:)
       IF (cldfsnow_idx > 0) THEN
          ! add in snow
          CALL get_snow_optics_sw(ncol,pcols,pver,dei, iciwp, snow_tau, snow_tau_w, snow_tau_w_g, snow_tau_w_f)
          DO i=1,ncol
             DO k=1,pver
                cldfprime(i,k)=MAX(cldn(i,k),cldfsnow(i,k))
                IF(cldfprime(i,k) > 0.)THEN
                   c_cld_tau    (1:nbndsw,i,k)= &
                        (cldfsnow(i,k)*snow_tau    (1:nbndsw,i,k) + cldn(i,k)*cld_tau    (1:nbndsw,i,k))/cldfprime(i,k)
                   c_cld_tau_w  (1:nbndsw,i,k)= &
                        (cldfsnow(i,k)*snow_tau_w  (1:nbndsw,i,k) + cldn(i,k)*cld_tau_w  (1:nbndsw,i,k))/cldfprime(i,k)
                   c_cld_tau_w_g(1:nbndsw,i,k)= &
                        (cldfsnow(i,k)*snow_tau_w_g(1:nbndsw,i,k) + cldn(i,k)*cld_tau_w_g(1:nbndsw,i,k))/cldfprime(i,k)
                   c_cld_tau_w_f(1:nbndsw,i,k)= &
                        (cldfsnow(i,k)*snow_tau_w_f(1:nbndsw,i,k) + cldn(i,k)*cld_tau_w_f(1:nbndsw,i,k))/cldfprime(i,k)
                ELSE
                   c_cld_tau    (1:nbndsw,i,k)= 0._r8
                   c_cld_tau_w  (1:nbndsw,i,k)= 0._r8
                   c_cld_tau_w_g(1:nbndsw,i,k)= 0._r8
                   c_cld_tau_w_f(1:nbndsw,i,k)= 0._r8
                ENDIF
             ENDDO
          ENDDO
       ELSE
          c_cld_tau    (1:nbndsw,1:ncol,:)= cld_tau    (:,1:ncol,:)
          c_cld_tau_w  (1:nbndsw,1:ncol,:)= cld_tau_w  (:,1:ncol,:)
          c_cld_tau_w_g(1:nbndsw,1:ncol,:)= cld_tau_w_g(:,1:ncol,:)
          c_cld_tau_w_f(1:nbndsw,1:ncol,:)= cld_tau_w_f(:,1:ncol,:)
          IF (.NOT.(cldfsnow_idx > 0)) THEN
             cldfprime(1:ncol,:)=cldn(1:ncol,:)
          ENDIF
       ENDIF

    END IF

    IF (dolw) THEN
       IF(oldcldoptics) THEN
          CALL cloud_rad_props_get_lw(ncol,pcols,pver,pdel,qi,ql,&     
               rei,dei,cldn,lamc,pgam,iclwp,iciwp,&
               cld_lw_abs, oldliq=.FALSE., oldice=.FALSE., oldcloud=.TRUE.)
       ELSE
          SELECT CASE (TRIM(icecldoptics))
          CASE ('ebertcurry')
             CALL ec_ice_get_rad_props_lw(ncol,pcols,pver,cldn,rei,qi,ql,pdel,iclwp,&
                  iciwp,&
                  ice_lw_abs, oldicewp=.TRUE.)
          CASE ('mitchell')
             CALL ice_cloud_get_rad_props_lw(ncol,pcols,pver,dei,iciwp,&
                  ice_lw_abs)
          CASE default
             STOP 'iccldoptics must be one either ebertcurry or mitchell'
          END SELECT
          SELECT CASE (TRIM(liqcldoptics))
          CASE ('slingo')
             CALL slingo_liq_get_rad_props_lw(ncol,pcols,pver,iciwp,qi,ql,pdel,iclwp,iciwp, &
                  liq_lw_abs, oldliqwp=.TRUE.)

          CASE ('gammadist')
             CALL liquid_cloud_get_rad_props_lw( ncol,pcols,pver,lamc,pgam,iclwp, &
                  liq_lw_abs)
          CASE default
             STOP 'liqcldoptics must be either slingo or gammadist'
          END SELECT
          cld_lw_abs(:,1:ncol,:) = liq_lw_abs(:,1:ncol,:) + ice_lw_abs(:,1:ncol,:)
       ENDIF
       !call cloud_rad_props_get_lw(state,  pbuf, cld_lw_abs, oldliq=.true., oldice=.true.)
       !call cloud_rad_props_get_lw(state,  pbuf, cld_lw_abs, oldcloud=.true.)
       !call cloud_rad_props_get_lw(state,  pbuf, cld_lw_abs, oldliq=.true., oldice=.true.)

       IF (cldfsnow_idx > 0) THEN
          ! add in snow
          CALL snow_cloud_get_rad_props_lw( ncol,pcols,pver,dei,iciwp,&
               snow_lw_abs)
          DO i=1,ncol
             DO k=1,pver
                cldfprime(i,k)=MAX(cldn(i,k),cldfsnow(i,k))
                IF(cldfprime(i,k) > 0.)THEN
                   c_cld_lw_abs(1:nbndlw,i,k)= &
                        (cldfsnow(i,k)*snow_lw_abs(1:nbndlw,i,k) + cldn(i,k)*cld_lw_abs(1:nbndlw,i,k))/cldfprime(i,k)
                ELSE
                   c_cld_lw_abs(1:nbndlw,i,k)= 0._r8
                ENDIF
             ENDDO
          ENDDO
       ELSE
          c_cld_lw_abs(1:nbndlw,1:ncol,:)=cld_lw_abs(:,1:ncol,:)
       ENDIF

       IF (.NOT.(cldfsnow_idx > 0)) THEN
          cldfprime(1:ncol,:)=cldn(1:ncol,:)
       ENDIF


    ENDIF

  END SUBROUTINE Run_Optical_Properties

  !==============================================================================

  SUBROUTINE snow_cloud_get_rad_props_lw( ncol,pcols,pver,dei,iciwpth,&
       abs_od)

    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver   
    REAL(r8), INTENT(in  ) :: dei    (pcols,pver)
    REAL(r8), INTENT(in  ) :: iciwpth(pcols,pver)

    REAL(r8), INTENT(out) :: abs_od(nlwbands,pcols,pver)

    REAL(r8) :: dlimited ! d limited to range dmin,dmax

    INTEGER :: i,k,i_d_grid,k_d_eff,i_lwband
    REAL(r8) :: wd, onemwd, absor

    abs_od = 0._r8


    ! note that this code makes the "ice path" point to the "snow path from CAM"

    ! call pbuf_get_field(pbuf, i_icswp, iciwpth)
    ! call pbuf_get_field(pbuf, i_des,   dei)

    ! note that this code makes the "ice path" point to the "snow path from CAM"

    DO i = 1,ncol
       DO k = 1,pver
          dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
          ! if ice water path is too small, OD := 0
          IF( iciwpth(i,k) < 1.e-80_r8 .OR. dlimited .EQ. 0._r8) THEN
             abs_od (:,i,k) = 0._r8
             !else if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
             !   write(iulog,*) 'dlimited prognostic cldwat2m',dlimited
             !   write(iulog,*) 'grid values of deff ice from optics file',g_d_eff(1),' -> ',g_d_eff(n_g_d)
             !   !call endrun ('deff of ice exceeds limits')
          ELSE
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             IF (dlimited <= g_d_eff(1)) THEN
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             ELSEIF (dlimited >= g_d_eff(n_g_d)) THEN
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             ELSE
                DO i_d_grid = 2, n_g_d
                   k_d_eff = i_d_grid
                   IF(g_d_eff(i_d_grid) > dlimited) EXIT
                ENDDO
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             ENDIF
             ! interpolate into grid and extract radiative properties
             DO i_lwband = 1, nlwbands
                absor = wd*abs_lw_ice(k_d_eff-1,i_lwband) + &
                     onemwd*abs_lw_ice(k_d_eff  ,i_lwband)
                abs_od (i_lwband,i,k)=  iciwpth(i,k) * absor 
             ENDDO
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE snow_cloud_get_rad_props_lw


  SUBROUTINE slingo_liq_get_rad_props_lw(ncol,pcols,pver,cldn,qi,ql,pdel,iclwpth,iciwpth,&
       abs_od, oldliqwp)

    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: cldn    (pcols,pver)
    REAL(r8),INTENT(IN   ) :: qi      (pcols,pver)
    REAL(r8),INTENT(IN   ) :: ql      (pcols,pver)
    REAL(r8),INTENT(IN   ) :: pdel    (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iclwpth (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iciwpth (pcols,pver)

    REAL(r8), INTENT(out) :: abs_od(nlwbands,pcols,pver)
    LOGICAL, INTENT(in) :: oldliqwp

    REAL(r8) :: gicewp(pcols,pver)
    REAL(r8) :: gliqwp(pcols,pver)
    REAL(r8) :: cicewp(pcols,pver)
    REAL(r8) :: cliqwp(pcols,pver)
    REAL(r8) :: ficemr(pcols,pver)
    REAL(r8) :: cwp   (pcols,pver)
    REAL(r8) :: cldtau(pcols,pver)

    INTEGER ::  lwband, i, k 

    REAL(r8) :: kabs
    REAL(r8), PARAMETER :: kabsl = 0.090361_r8          ! longwave liquid absorption coeff (m**2/g)
    REAL(R8),PARAMETER :: SHR_CONST_G = 9.80616_R8      ! acceleration of gravity ~ m/s^2


    !itim  =  pbuf_old_tim_idx()
    !call pbuf_get_field(pbuf, rei_idx,   rei)
    !call pbuf_get_field(pbuf, cld_idx,   cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/))

    IF (oldliqwp) THEN
       DO k=1,pver
          DO i = 1,ncol
             gicewp(i,k) = qi(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = ql(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
             ficemr(i,k) = qi(i,k) /MAX(1.e-10_r8,(qi(i,k)+ql(i,k)))
          END DO
       END DO
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
    ELSE
       !if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
       !   call endrun('slingo_liq_get_rad_props_lw: oldliqwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
       !endif
       !call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
       !call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
       DO k=1,pver
          DO i = 1,ncol
             cwp   (i,k) = 1000.0_r8 * iclwpth(i,k) + 1000.0_r8 * iciwpth(i, k)
             ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(MAX(1.e-18_r8, cwp(i,k)))
          END DO
       END DO
    ENDIF


    DO k=1,pver
       DO i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use liquid water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabs = kabsl*(1._r8-ficemr(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       END DO
    END DO
    !
    DO lwband = 1,nlwbands
       abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
    ENDDO


  END SUBROUTINE slingo_liq_get_rad_props_lw


  !==============================================================================

  SUBROUTINE ec_ice_get_rad_props_lw(ncol,pcols,pver,cldn,rei,qi,ql,pdel,iclwpth,&
       iciwpth,&
       abs_od, oldicewp)
    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: cldn    (pcols,pver)
    REAL(r8),INTENT(IN   ) :: rei     (pcols,pver)
    REAL(r8),INTENT(IN   ) :: qi      (pcols,pver)
    REAL(r8),INTENT(IN   ) :: ql      (pcols,pver)
    REAL(r8),INTENT(IN   ) :: pdel    (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iclwpth (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iciwpth (pcols,pver)
    REAL(r8), INTENT(out) :: abs_od(nlwbands,pcols,pver)
    LOGICAL, INTENT(in) :: oldicewp

    REAL(r8) :: gicewp(pcols,pver)
    REAL(r8) :: gliqwp(pcols,pver)
    REAL(r8) :: cicewp(pcols,pver)
    REAL(r8) :: cliqwp(pcols,pver)
    REAL(r8) :: ficemr(pcols,pver)
    REAL(r8) :: cwp   (pcols,pver)
    REAL(r8) :: cldtau(pcols,pver)

    INTEGER :: lwband, i, k

    REAL(r8) :: kabs, kabsi

    REAL(r8), PARAMETER :: kabsl = 0.090361_r8         ! longwave liquid absorption coeff (m**2/g)
    REAL(R8),PARAMETER :: SHR_CONST_G = 9.80616_R8      ! acceleration of gravity ~ m/s^2




    !   itim  =  pbuf_old_tim_idx()
    !   call pbuf_get_field(pbuf, rei_idx,   rei)
    !   call pbuf_get_field(pbuf, cld_idx,   cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/))


    IF(oldicewp) THEN
       DO k=1,pver
          DO i = 1,ncol
             gicewp(i,k) = qi(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = ql(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
             ficemr(i,k) = qi(i,k) /MAX(1.e-10_r8,(qi(i,k)+ql(i,k)))
          END DO
       END DO
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
    ELSE
       !if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
       !   call endrun('ec_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
       !endif
       !call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
       !call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
       DO k=1,pver
          DO i = 1,ncol
             cwp(i,k)    = 1000.0_r8 *iciwpth(i,k) + 1000.0_r8 *iclwpth(i,k)
             ficemr(i,k) = 1000.0_r8 *iciwpth(i,k)/(MAX(1.e-18_r8,cwp(i,k)))
          END DO
       END DO
    ENDIF

    DO k=1,pver
       DO i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use ice water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/MIN(MAX(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs =  kabsi*ficemr(i,k) ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       END DO
    END DO
    !
    DO lwband = 1,nlwbands
       abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
    ENDDO

    !if(oldicewp) then
    !  call outfld('CIWPTH_OLD',cicewp(:,:)/1000,pcols,lchnk)
    !else
    !  call outfld('CIWPTH_OLD',iciwpth(:,:),pcols,lchnk)
    !endif
    !call outfld('CI_OD_LW_OLD',cldtau(:,:),pcols,lchnk)

  END SUBROUTINE ec_ice_get_rad_props_lw
  !==============================================================================


  !==============================================================================

  SUBROUTINE cloud_rad_props_get_lw(ncol,pcols,pver,pdel,qi,ql,&     
       rei,dei,cldn,lamc,pgam,iclwpth,iciwpth,&
       cld_abs_od, oldliq, oldice, oldcloud)

    ! Purpose: Compute cloud longwave absorption optical depth
    !    cloud_rad_props_get_lw() is called by radlw() 

    ! Arguments
    INTEGER ,INTENT(IN   ) :: ncol
    INTEGER ,INTENT(IN   ) :: pcols
    INTEGER ,INTENT(IN   ) :: pver
    REAL(r8),INTENT(in   ) :: pdel   (pcols,pver)
    REAL(r8),INTENT(in   ) :: qi     (pcols,pver)
    REAL(r8),INTENT(in   ) :: ql     (pcols,pver)
    REAL(r8),INTENT(in   ) :: rei    (pcols,pver)
    REAL(r8),INTENT(in   ) :: dei    (pcols,pver)
    REAL(r8),INTENT(in   ) :: cldn   (pcols,pver)
    REAL(r8),INTENT(IN   ) :: lamc   (pcols,pver)
    REAL(r8),INTENT(IN   ) :: pgam   (pcols,pver)
    REAL(r8),INTENT(in   ) :: iclwpth(pcols,pver)
    REAL(r8),INTENT(in   ) :: iciwpth(pcols,pver)

    REAL(r8),            INTENT(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
    LOGICAL, OPTIONAL,   INTENT(in)  :: oldliq  ! use old liquid optics
    LOGICAL, OPTIONAL,   INTENT(in)  :: oldice  ! use old ice optics
    LOGICAL, OPTIONAL,   INTENT(in)  :: oldcloud  ! use old optics for both (b4b)

    ! Local variables

!    INTEGER :: bnd_idx     ! LW band index
!    INTEGER :: i           ! column index
!    INTEGER :: k           ! lev index

    ! rad properties for liquid clouds
    REAL(r8) :: liq_tau_abs_od(nlwbands,pcols,pver) ! liquid cloud absorption optical depth

    ! rad properties for ice clouds
    REAL(r8) :: ice_tau_abs_od(nlwbands,pcols,pver) ! ice cloud absorption optical depth

    !-----------------------------------------------------------------------------

    ! compute optical depths cld_absod 
    cld_abs_od = 0._r8

    IF(PRESENT(oldcloud))THEN
       IF(oldcloud) THEN
          ! make diagnostic calls to these first to output ice and liq OD's
          ! call old_liq_get_rad_props_lw(state, pbuf, liq_tau_abs_od, oldliqwp=.false.)
          ! call old_ice_get_rad_props_lw(state, pbuf, ice_tau_abs_od, oldicewp=.false.)
          ! This affects climate (cld_abs_od)
          CALL oldcloud_lw(ncol,pcols,pver,rei,cldn,qi,ql,pdel,iclwpth,iciwpth, &
               cld_abs_od,oldwp=.FALSE.)

          RETURN
       ENDIF
    ENDIF

    IF(PRESENT(oldliq))THEN
       IF(oldliq) THEN
          CALL old_liq_get_rad_props_lw(ncol,pcols,pver,rei,cldn,qi,ql,pdel,iclwpth,iciwpth,&
               liq_tau_abs_od, oldliqwp=.FALSE.)
       ELSE
          CALL liquid_cloud_get_rad_props_lw( ncol,pcols,pver,lamc,pgam,iclwpth, &
               liq_tau_abs_od)
       ENDIF
    ELSE
       CALL liquid_cloud_get_rad_props_lw( ncol,pcols,pver,lamc,pgam,iclwpth, &
            liq_tau_abs_od)
    ENDIF

    IF(PRESENT(oldice))THEN
       IF(oldice) THEN
          CALL old_ice_get_rad_props_lw(ncol,pcols,pver,cldn,rei,iclwpth,iciwpth,& 
               qi,ql,pdel,&
               ice_tau_abs_od, oldicewp=.FALSE.)

       ELSE
          CALL ice_cloud_get_rad_props_lw(ncol,pcols,pver,dei,iciwpth,&
               ice_tau_abs_od)
       ENDIF
    ELSE
       CALL ice_cloud_get_rad_props_lw(ncol,pcols,pver,dei,iciwpth,&
            ice_tau_abs_od)
    ENDIF
    cld_abs_od(:,1:ncol,:) = liq_tau_abs_od(:,1:ncol,:) + ice_tau_abs_od(:,1:ncol,:) 

  END SUBROUTINE cloud_rad_props_get_lw


  !==============================================================================
  ! Private methods
  !==============================================================================
  !==============================================================================
  !==============================================================================

  SUBROUTINE liquid_cloud_get_rad_props_lw( ncol,pcols,pver,lamc,pgam,iclwpth, &
       abs_od)
    INTEGER ,INTENT(IN   ) :: ncol
    INTEGER ,INTENT(IN   ) :: pcols
    INTEGER ,INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: lamc   (pcols,pver)
    REAL(r8),INTENT(IN   ) :: pgam   (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iclwpth(pcols,pver)
    REAL(r8),INTENT(out  ) :: abs_od (nlwbands,pcols,pver)

    INTEGER :: i, k

    abs_od = 0._r8

    DO k = 1,pver
       DO i = 1,ncol
          IF(lamc(i,k) > 0._r8) THEN ! This seems to be the clue for no cloud from microphysics formulation
             CALL gam_liquid_lw(iclwpth(i,k), lamc(i,k), pgam(i,k), abs_od(1:nlwbands,i,k))
          ELSE
             abs_od(1:nlwbands,i,k) = 0._r8
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE liquid_cloud_get_rad_props_lw
  !==============================================================================
  ! Private methods
  !==============================================================================
  !==============================================================================

  SUBROUTINE old_ice_get_rad_props_lw(ncol,pcols,pver,cldn,rei,iclwpth,iciwpth,& 
       qi,ql,pdel,abs_od, oldicewp)

    INTEGER , INTENT(IN   ) :: ncol
    INTEGER , INTENT(IN   ) :: pcols
    INTEGER , INTENT(IN   ) :: pver
    REAL(r8), INTENT(in   ) :: cldn    (pcols,pver)
    REAL(r8), INTENT(in   ) :: rei     (pcols,pver)
    REAL(r8), INTENT(in   ) :: iclwpth (pcols,pver)
    REAL(r8), INTENT(in   ) :: iciwpth (pcols,pver)
    REAL(r8), INTENT(in   ) :: qi      (pcols,pver)
    REAL(r8), INTENT(in   ) :: ql      (pcols,pver)
    REAL(r8), INTENT(in   ) :: pdel    (pcols,pver)

    REAL(r8), INTENT(out  ) :: abs_od  (nlwbands,pcols,pver)
    LOGICAL , INTENT(in   ) :: oldicewp

    REAL(r8) :: gicewp(pcols,pver)
    REAL(r8) :: gliqwp(pcols,pver)
    REAL(r8) :: cicewp(pcols,pver)
    REAL(r8) :: cliqwp(pcols,pver)
    REAL(r8) :: ficemr(pcols,pver)
    REAL(r8) :: cwp   (pcols,pver)
    REAL(r8) :: cldtau(pcols,pver)

    INTEGER  ::  lwband, i, k

    REAL(r8) :: kabs, kabsi

    REAL(r8),PARAMETER :: kabsl       = 0.090361_r8     ! longwave liquid absorption coeff (m**2/g)
    REAL(R8),PARAMETER :: SHR_CONST_G = 9.80616_R8      ! acceleration of gravity ~ m/s^2


    IF(oldicewp) THEN
       DO k=1,pver
          DO i = 1,ncol
             gicewp(i,k) = qi(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = ql(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
             ficemr(i,k) = qi(i,k) /MAX(1.e-10_r8,(qi(i,k)+ql(i,k)))
          END DO
       END DO
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
    ELSE
       !if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
       !   call endrun('old_ice_get_rad_props_lw: oldicewp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
       !endif
       !call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
       !call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
       DO k=1,pver
          DO i = 1,ncol
             cwp(i,k)    = 1000.0_r8 *iciwpth(i,k) + 1000.0_r8 *iclwpth(i,k)
             ficemr(i,k) = 1000.0_r8 *iciwpth(i,k)/(MAX(1.e-18_r8,cwp(i,k)))
          END DO
       END DO
    ENDIF

    DO k=1,pver
       DO i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use ice water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/MIN(MAX(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs =  kabsi*ficemr(i,k) ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       END DO
    END DO
    !
    DO lwband = 1,nlwbands
       abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
    ENDDO

    !if(oldicewp) then
    !  call outfld('CIWPTH_OLD',cicewp(:,:)/1000,pcols,lchnk)
    !else
    !  call outfld('CIWPTH_OLD',iciwpth(:,:),pcols,lchnk)
    !endif
    !call outfld('CI_OD_LW_OLD',cldtau(:,:),pcols,lchnk)

  END SUBROUTINE old_ice_get_rad_props_lw
  !==============================================================================

  !==============================================================================
  ! Private methods
  !==============================================================================
  !==============================================================================
  !==============================================================================


  !==============================================================================

  SUBROUTINE ice_cloud_get_rad_props_lw(ncol,pcols,pver,dei,iciwpth,abs_od)

    INTEGER ,INTENT(IN   ) :: ncol
    INTEGER ,INTENT(IN   ) :: pcols
    INTEGER ,INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: dei     (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iciwpth (pcols,pver)
    REAL(r8),INTENT(out  ) :: abs_od  (nlwbands,pcols,pver)

    REAL(r8) :: dlimited ! d limited to range dmin,dmax

    INTEGER :: i,k,i_d_grid,k_d_eff,i_lwband
    REAL(r8) :: wd, onemwd, absor

    abs_od = 0._r8

    DO i = 1,ncol
       DO k = 1,pver
          dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
          ! if ice water path is too small, OD := 0
          IF( iciwpth(i,k) < 1.e-80_r8 .OR. dlimited .EQ. 0._r8) THEN
             abs_od (:,i,k) = 0._r8
             !else if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
             !   write(iulog,*) 'dlimited prognostic cldwat2m',dlimited
             !   write(iulog,*) 'grid values of deff ice from optics file',g_d_eff(1),' -> ',g_d_eff(n_g_d)
             !   !call endrun ('deff of ice exceeds limits')
          ELSE
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             IF (dlimited <= g_d_eff(1)) THEN
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             ELSEIF (dlimited >= g_d_eff(n_g_d)) THEN
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             ELSE
                DO i_d_grid = 2, n_g_d
                   k_d_eff = i_d_grid
                   IF(g_d_eff(i_d_grid) > dlimited) EXIT
                ENDDO
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             ENDIF
             ! interpolate into grid and extract radiative properties
             DO i_lwband = 1, nlwbands
                absor = wd*abs_lw_ice(k_d_eff-1,i_lwband) + &
                     onemwd*abs_lw_ice(k_d_eff  ,i_lwband)
                abs_od (i_lwband,i,k)=  iciwpth(i,k) * absor 
             ENDDO
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE ice_cloud_get_rad_props_lw

  !==============================================================================


  !==============================================================================
  ! Private methods
  !==============================================================================
  !==============================================================================

  SUBROUTINE gam_liquid_lw(clwptn, lamc, pgam, abs_od)
    REAL(r8), INTENT(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
    REAL(r8), INTENT(in) :: lamc   ! prognosed value of lambda for cloud
    REAL(r8), INTENT(in) :: pgam   ! prognosed value of mu for cloud
    REAL(r8), INTENT(out) :: abs_od(1:nlwbands)
    ! for interpolating into mu/lambda
    INTEGER :: imu, kmu
    INTEGER :: ilambda, klambda
    INTEGER :: lwband ! sw band index
    REAL(r8) :: absc, wmu, onemwmu ,wlambda, onemwlambda, lambdaplus, lambdaminus

    IF (clwptn < 1.e-80_r8) THEN
       abs_od = 0._r8
       RETURN
    ENDIF

    IF (pgam < g_mu(1) .OR. pgam > g_mu(nmu)) THEN
       WRITE(iulog,*)'pgam from prognostic cldwat2m',pgam
       WRITE(iulog,*)'g_mu from file',g_mu
       STOP'pgam exceeds limits' !call endrun ('pgam exceeds limits')
    ENDIF
    DO imu = 1, nmu
       kmu = imu
       IF (g_mu(kmu) > pgam) EXIT
    ENDDO
    wmu = (g_mu(kmu) - pgam)/(g_mu(kmu) - g_mu(kmu-1))
    onemwmu = 1._r8 - wmu

    DO ilambda = 1, nlambda
       klambda = ilambda
       IF (wmu*g_lambda(kmu-1,ilambda) + onemwmu*g_lambda(kmu,ilambda) < lamc) EXIT
    ENDDO
    IF (klambda <= 1 .OR. klambda > nlambda) STOP'lamc  exceeds limits'! call endrun('lamc  exceeds limits')
    lambdaplus = wmu*g_lambda(kmu-1,klambda  ) + onemwmu*g_lambda(kmu,klambda  )
    lambdaminus= wmu*g_lambda(kmu-1,klambda-1) + onemwmu*g_lambda(kmu,klambda-1)
    wlambda = (lambdaplus - lamc) / (lambdaplus - lambdaminus)
    onemwlambda = 1._r8 - wlambda

    DO lwband = 1, nlwbands
       absc=     wlambda*    wmu*abs_lw_liq(kmu-1,klambda-1,lwband) + &
            onemwlambda*    wmu*abs_lw_liq(kmu-1,klambda  ,lwband) + &
            wlambda*onemwmu*abs_lw_liq(kmu  ,klambda-1,lwband) + &
            onemwlambda*onemwmu*abs_lw_liq(kmu  ,klambda  ,lwband)

       abs_od(lwband) = clwptn * absc
    ENDDO

    RETURN
  END SUBROUTINE gam_liquid_lw


  !==============================================================================
  ! Private methods
  !==============================================================================

  !==============================================================================
  SUBROUTINE old_liq_get_rad_props_lw(ncol,pcols,pver,rei,cldn,qi,ql,pdel,iclwpth,iciwpth,&
       abs_od, oldliqwp)

    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver
    REAL(r8), INTENT(in ) :: rei    (pcols,pver)
    REAL(r8), INTENT(in ) :: pdel   (pcols,pver)
    REAL(r8), INTENT(in ) :: cldn   (pcols,pver)
    REAL(r8), INTENT(in ) :: qi     (pcols,pver)
    REAL(r8), INTENT(in ) :: ql     (pcols,pver)
    REAL(r8), INTENT(in ) :: iclwpth(pcols,pver)
    REAL(r8), INTENT(in ) :: iciwpth(pcols,pver)
    REAL(r8), INTENT(out) :: abs_od(nlwbands,pcols,pver)
    LOGICAL, INTENT(in) :: oldliqwp

    REAL(r8) :: gicewp(pcols,pver)
    REAL(r8) :: gliqwp(pcols,pver)
    REAL(r8) :: cicewp(pcols,pver)
    REAL(r8) :: cliqwp(pcols,pver)
    REAL(r8) :: ficemr(pcols,pver)
    REAL(r8) :: cwp(pcols,pver)
    REAL(r8) :: cldtau(pcols,pver)

    INTEGER :: lwband, i, k 

    REAL(r8) :: kabs, kabsi
    REAL(r8),PARAMETER :: kabsl = 0.090361_r8                 ! longwave liquid absorption coeff (m**2/g)
    REAL(R8),PARAMETER :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2



    ! call pbuf_get_field(pbuf, rei_idx,   rei)
    ! call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim/), kount=(/pcols,pver,1/))

    IF (oldliqwp) THEN
       DO k=1,pver
          DO i = 1,ncol
             gicewp(i,k) = qi  (i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = ql  (i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
             ficemr(i,k) = qi(i,k) /MAX(1.e-10_r8,(qi(i,k)+ql(i,k)))
          END DO
       END DO
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
    ELSE
       ! if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
       !    call endrun('old_liq_get_rad_props_lw: oldliqwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
       ! endif
       ! call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
       ! call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
       DO k=1,pver
          DO i = 1,ncol
             cwp(i,k)    = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
             ficemr(i,k) = 1000.0 * iciwpth(i,k)/(MAX(1.e-18_r8,cwp(i,k)))
          END DO
       END DO
    ENDIF


    DO k=1,pver
       DO i=1,ncol

          ! Note from Andrew Conley:
          !  Optics for RK no longer supported, This is constructed to get
          !  close to bit for bit.  Otherwise we could simply use liquid water path
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/MIN(MAX(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs = kabsl*(1._r8-ficemr(i,k)) ! + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       END DO
    END DO
    !
    DO lwband = 1,nlwbands
       abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
    ENDDO


  END SUBROUTINE old_liq_get_rad_props_lw
  !==============================================================================


  SUBROUTINE cloud_rad_props_init(fNameCldOptSW ,fNameCldOptLW )
    CHARACTER(LEN=*), INTENT(IN) :: fNameCldOptSW 
    CHARACTER(LEN=*), INTENT(IN) :: fNameCldOptLW 
    INTEGER :: ios
    ALLOCATE(g_mu      (nmu))
    ALLOCATE(g_lambda  (nmu,nlambda))
    ALLOCATE(ext_sw_liq(nmu,nlambda,nswbands) )
    ALLOCATE(ssa_sw_liq(nmu,nlambda,nswbands))
    ALLOCATE(asm_sw_liq(nmu,nlambda,nswbands))
    ALLOCATE(abs_lw_liq(nmu,nlambda,nlwbands))
    OPEN(1,FILE=TRIM(fNameCldOptSW),&
         ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED',&  
         STATUS='OLD',&
         ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
       WRITE(0, '(/a)') TRIM(fNameCldOptSW)
       WRITE(0, '(/a)') &
           '*** Error: Spectral file could not be opened.'
        STOP
    END IF

    READ(1)g_mu   
    READ(1)g_lambda  
    READ(1)ext_sw_liq
    READ(1)ssa_sw_liq
    READ(1)asm_sw_liq
    READ(1)abs_lw_liq
    CLOSE(1,STATUS='KEEP')
    ! I forgot to convert kext from m^2/Volume to m^2/Kg
    ext_sw_liq = ext_sw_liq / 0.9970449e3_r8 
    abs_lw_liq = abs_lw_liq / 0.9970449e3_r8 

    ALLOCATE(g_d_eff(n_g_d))
    ALLOCATE(abs_lw_ice(n_g_d,nlwbands))
    ALLOCATE(ext_sw_ice(n_g_d,nswbands))
    ALLOCATE(ssa_sw_ice(n_g_d,nswbands))
    ALLOCATE(asm_sw_ice(n_g_d,nswbands))
    OPEN(1,FILE=TRIM(fNameCldOptLW),&
         ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED',&  
         STATUS='OLD',&
         ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
       WRITE(0, '(/a)') TRIM(fNameCldOptLW)
       WRITE(0, '(/a)') &
           '*** Error: Spectral file could not be opened.'
        STOP 
    END IF

    READ(1)g_d_eff
    READ(1)abs_lw_ice
    READ(1)ext_sw_ice
    READ(1)ssa_sw_ice
    READ(1)asm_sw_ice
    CLOSE(1,STATUS='KEEP')

  END SUBROUTINE cloud_rad_props_init
  !==============================================================================
  ! Private methods
  !==============================================================================
  !==============================================================================

  SUBROUTINE oldcloud_lw(ncol,pcols,pver,rei,cldn,qi,ql,pdel,iclwpth,iciwpth, &
       cld_abs_od,oldwp)
    INTEGER , INTENT(IN   ) :: ncol
    INTEGER , INTENT(IN   ) :: pcols
    INTEGER , INTENT(IN   ) :: pver
    REAL(r8),INTENT(in ) :: pdel(pcols,pver)
    REAL(r8),INTENT(in ) :: qi  (pcols,pver)
    REAL(r8),INTENT(in ) :: ql  (pcols,pver)
    REAL(r8),INTENT(in ) :: rei (pcols,pver)
    REAL(r8),INTENT(in ) :: cldn(pcols,pver)
    REAL(r8),INTENT(in ) :: iclwpth(pcols,pver)
    REAL(r8),INTENT(in ) :: iciwpth(pcols,pver)
    REAL(r8),INTENT(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer

    LOGICAL ,INTENT(in ) :: oldwp                           ! use old definition of waterpath


    REAL(r8) :: gicewp (pcols,pver)
    REAL(r8) :: gliqwp (pcols,pver)
    REAL(r8) :: cicewp (pcols,pver)
    REAL(r8) :: cliqwp (pcols,pver)
    REAL(r8) :: ficemr (pcols,pver)
    REAL(r8) :: cwp    (pcols,pver)
    REAL(r8) :: cldtau (pcols,pver)

    INTEGER ::  lwband, i, k
    REAL(r8) :: kabs, kabsi
    REAL(r8),PARAMETER :: kabsl = 0.090361_r8  ! longwave liquid absorption coeff (m**2/g)
    REAL(R8),PARAMETER :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2




    !   itim  =  pbuf_old_tim_idx()
    !   call pbuf_get_field(pbuf, rei_idx,   rei)
    !   call pbuf_get_field(pbuf, cld_idx,   cldn,   start=(/1,1,itim/), kount=(/pcols,pver,1/))

    IF (oldwp) THEN
       DO k=1,pver
          DO i = 1,ncol
             gicewp(i,k) = qi(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = ql(i,k)*pdel(i,k)/SHR_CONST_G*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / MAX(0.01_r8,cldn(i,k))                 ! In-cloud liquid water path.
             ficemr(i,k) = qi(i,k)/MAX(1.e-10_r8,(qi(i,k)+ql(i,k)))
          END DO
       END DO
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)
    ELSE
       !if (iclwp_idx<=0 .or. iciwp_idx<=0) then 
       !
       !         call endrun('oldcloud_lw: oldwp must be set to true since ICIWP and/or ICLWP were not found in pbuf')
       !      endif
       !      call pbuf_get_field(pbuf, iclwp_idx, iclwpth)
       !      call pbuf_get_field(pbuf, iciwp_idx, iciwpth)
       DO k=1,pver
          DO i = 1,ncol
             cwp(i,k) = 1000.0_r8 *iclwpth(i,k) + 1000.0_r8 *iciwpth(i, k)
             ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(MAX(1.e-18_r8,cwp(i,k)))
          END DO
       END DO
    ENDIF

    DO k=1,pver
       DO i=1,ncol

          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabsi = 0.005_r8 + 1._r8/MIN(MAX(13._r8,scalefactor*rei(i,k)),130._r8)
          kabs = kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          !emis(i,k) = 1._r8 - exp(-1.66_r8*kabs*clwp(i,k))
          cldtau(i,k) = kabs*cwp(i,k)
       END DO
    END DO
    !
    DO lwband = 1,nlwbands
       cld_abs_od(lwband,1:ncol,1:pver)=cldtau(1:ncol,1:pver)
    ENDDO

  END SUBROUTINE oldcloud_lw


  !==============================================================================
  ! Private methods
  !==============================================================================

  SUBROUTINE get_ice_optics_sw   (ncol,pcols,pver,dei,iciwpth,&
       tau, tau_w, tau_w_g, tau_w_f)

    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: dei    (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iciwpth(pcols,pver)
    REAL(r8),INTENT(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
    REAL(r8),INTENT(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
    REAL(r8),INTENT(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
    REAL(r8),INTENT(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

    REAL(r8) :: dlimited ! d limited to dmin,dmax range

    INTEGER :: i,k,i_d_grid,k_d_eff,i_swband
    REAL(r8) :: wd, onemwd, ext, ssa, asm

    tau    =0.0_r8;tau_w  =0.0_r8;tau_w_g=0.0_r8;tau_w_f=0.0_r8
    ! note that this code makes the "ice path" 

    !call pbuf_get_field(pbuf, i_iciwp, iciwpth)
    !call pbuf_get_field(pbuf, i_dei,   dei)

    DO i = 1,ncol
       DO k = 1,pver
          dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
          IF( iciwpth(i,k) < 1.e-80_r8 .OR. dlimited .EQ. 0._r8) THEN
             ! if ice water path is too small, OD := 0
             tau    (:,i,k) = 0._r8
             tau_w  (:,i,k) = 0._r8
             tau_w_g(:,i,k) = 0._r8
             tau_w_f(:,i,k) = 0._r8
          ELSE 
             IF (dlimited < g_d_eff(1) .OR. dlimited > g_d_eff(n_g_d)) THEN
                WRITE(iulog,*) 'dei from prognostic cldwat2m',dei(i,k)
                WRITE(iulog,*) 'grid values of deff ice from optics file',g_d_eff
                STOP 'deff of ice exceeds limits'
                !call endrun ('deff of ice exceeds limits')
             ENDIF
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             IF (dlimited <= g_d_eff(1)) THEN
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             ELSEIF (dlimited >= g_d_eff(n_g_d)) THEN
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             ELSE
                DO i_d_grid = 1, n_g_d
                   k_d_eff = i_d_grid
                   IF(g_d_eff(i_d_grid) > dlimited) EXIT
                ENDDO
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             ENDIF
             ! interpolate into grid and extract radiative properties
             DO i_swband = 1, nswbands
                ext = wd*ext_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*ext_sw_ice(k_d_eff  ,i_swband) 
                ssa = wd*ssa_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*ssa_sw_ice(k_d_eff  ,i_swband) 
                asm = wd*asm_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*asm_sw_ice(k_d_eff  ,i_swband) 
                tau    (i_swband,i,k)=iciwpth(i,k) * ext
                tau_w  (i_swband,i,k)=iciwpth(i,k) * ext * ssa
                tau_w_g(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm
                tau_w_f(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm * asm
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_ice_optics_sw

  !==============================================================================


  !==============================================================================

  SUBROUTINE cloud_rad_props_get_sw(ncol,pcols,pver,cldn,rei,iciwp,qi,pdel, &
       tau, tau_w, tau_w_g, tau_w_f)

    ! return totaled (across all species) layer tau, omega, g, f 
    ! for all spectral interval for aerosols affecting the climate

    ! Arguments
    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: cldn(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: rei(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: iciwp(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: qi(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: pdel(pcols,pver) 

    REAL(r8), INTENT(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
    REAL(r8), INTENT(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
    REAL(r8), INTENT(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
    REAL(r8), INTENT(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

!    LOGICAL, OPTIONAL, INTENT(in) :: oldliq
!    LOGICAL, OPTIONAL, INTENT(in) :: oldice

    ! Local variables

!    INTEGER :: k, i    ! lev and daycolumn indices
!    INTEGER :: iswband ! sw band indices

    !-----------------------------------------------------------------------------


    ! initialize to conditions that would cause failure
    tau     (:,:,:) = -100._r8
    tau_w   (:,:,:) = -100._r8
    tau_w_g (:,:,:) = -100._r8
    tau_w_f (:,:,:) = -100._r8

    ! initialize layers to accumulate od's
    tau    (:,1:ncol,:) = 0._r8
    tau_w  (:,1:ncol,:) = 0._r8
    tau_w_g(:,1:ncol,:) = 0._r8
    tau_w_f(:,1:ncol,:) = 0._r8


    CALL ec_ice_optics_sw   (ncol,pcols,pver,cldn,rei,iciwp,qi,pdel, tau, tau_w,&
         tau_w_g, tau_w_f, oldicewp=.TRUE.)


  END SUBROUTINE cloud_rad_props_get_sw
  !==============================================================================

  !==============================================================================

  SUBROUTINE slingo_cloud_rad_props_get_sw(ncol,pcols,pver,cldn,rel,iclwp,ql,pdel,&
       tau, tau_w, tau_w_g, tau_w_f)

    ! return totaled (across all species) layer tau, omega, g, f 
    ! for all spectral interval for aerosols affecting the climate

    ! Arguments
    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pver
    INTEGER, INTENT(IN   ) :: pcols
    REAL(r8),INTENT(IN   ) :: cldn(pcols,pver)
    REAL(r8),INTENT(IN   ) :: rel (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iclwp(pcols,pver)
    REAL(r8),INTENT(IN   ) :: ql   (pcols,pver)
    REAL(r8),INTENT(IN   ) :: pdel  (pcols,pver)


    REAL(r8), INTENT(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
    REAL(r8), INTENT(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
    REAL(r8), INTENT(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
    REAL(r8), INTENT(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

    ! Local variables

!    INTEGER :: i    ! lev and daycolumn indices
!    INTEGER :: iswband ! sw band indices

!    REAL(r8) :: liq_tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
!    REAL(r8) :: liq_tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
!    REAL(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
!    REAL(r8) :: liq_tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w


    !-----------------------------------------------------------------------------


    CALL slingo_liq_optics_sw( ncol,pcols,pver,cldn,rel,iclwp,ql,pdel,&
         tau, tau_w, tau_w_g, tau_w_f, oldliqwp=.TRUE. )

  END SUBROUTINE slingo_cloud_rad_props_get_sw
  !==============================================================================


  !==============================================================================
  ! Private methods
  !==============================================================================

  SUBROUTINE ec_ice_optics_sw   (ncol,pcols,pver,cldn,rei,iciwp,qi,pdel,ice_tau, ice_tau_w,&
       ice_tau_w_g, ice_tau_w_f, oldicewp)

    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver
    REAL(r8),INTENT(IN   ) :: cldn(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: rei(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: iciwp(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: qi(pcols,pver) 
    REAL(r8),INTENT(IN   ) :: pdel(pcols,pver) 
    REAL(r8),INTENT(out) :: ice_tau    (nswbands,pcols,pver) ! extinction optical depth
    REAL(r8),INTENT(out) :: ice_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
    REAL(r8),INTENT(out) :: ice_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
    REAL(r8),INTENT(out) :: ice_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w
    LOGICAL, INTENT(in ) :: oldicewp

    REAL(r8), DIMENSION(pcols,pver) :: cicewp
    REAL(r8), DIMENSION(nswbands) :: wavmin
    REAL(r8), DIMENSION(nswbands) :: wavmax
    REAL(R8),PARAMETER :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2

    !
    ! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
    REAL(r8) :: abari(4) = &     ! a coefficient for extinction optical depth
         (/ 3.448e-03_r8, 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8/)
    REAL(r8) :: bbari(4) = &     ! b coefficient for extinction optical depth
         (/ 2.431_r8    , 2.431_r8    ,2.431_r8    ,2.431_r8    /)
    REAL(r8) :: cbari(4) = &     ! c coefficient for single scat albedo
         (/ 1.00e-05_r8 , 1.10e-04_r8 ,1.861e-02_r8,.46658_r8   /)
    REAL(r8) :: dbari(4) = &     ! d coefficient for single scat albedo
         (/ 0.0_r8      , 1.405e-05_r8,8.328e-04_r8,2.05e-05_r8 /)
    REAL(r8) :: ebari(4) = &     ! e coefficient for asymmetry parameter
         (/ 0.7661_r8   , 0.7730_r8   ,0.794_r8    ,0.9595_r8   /)
    REAL(r8) :: fbari(4) = &     ! f coefficient for asymmetry parameter
         (/ 5.851e-04_r8, 5.665e-04_r8,7.267e-04_r8,1.076e-04_r8/)

    REAL(r8) :: abarii           ! A coefficient for current spectral band
    REAL(r8) :: bbarii           ! B coefficient for current spectral band
    REAL(r8) :: cbarii           ! C coefficient for current spectral band
    REAL(r8) :: dbarii           ! D coefficient for current spectral band
    REAL(r8) :: ebarii           ! E coefficient for current spectral band
    REAL(r8) :: fbarii           ! F coefficient for current spectral band

    ! Minimum cloud amount (as a fraction of the grid-box area) to 
    ! distinguish from clear sky
    REAL(r8), PARAMETER :: cldmin = 1.0e-80_r8

    ! Decimal precision of cloud amount (0 -> preserve full resolution;
    ! 10^-n -> preserve n digits of cloud amount)
    REAL(r8), PARAMETER :: cldeps = 0.0_r8

    INTEGER :: ns, i, k, indxsl,  Nday
    REAL(r8) :: tmp1i, tmp2i, tmp3i, g

    Nday = ncol

    !itim = pbuf_old_tim_idx()
    !call pbuf_get_field(pbuf, cld_idx,cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/))
    !call pbuf_get_field(pbuf, rei_idx,rei)

    IF(oldicewp) THEN
       DO k=1,pver
          DO i = 1,Nday
             cicewp(i,k) = 1000.0_r8*qi(i,k)*pdel(i,k) /(SHR_CONST_G* MAX(0.01_r8,cldn(i,k)))
          END DO
       END DO
    ELSE
       !     if (iciwp_idx<=0) then 
       !        call endrun('ec_ice_optics_sw: oldicewp must be set to true since ICIWP was not found in pbuf')
       !     endif
       !call pbuf_get_field(pbuf, iciwp_idx, iciwp)
       cicewp(1:pcols,1:pver) =  1000.0_r8*iciwp(1:pcols,1:pver)
    ENDIF

    CALL get_sw_spectral_boundaries(wavmin,wavmax,'microns')

    DO ns = 1, nswbands

       IF(wavmax(ns) <= 0.7_r8) THEN
          indxsl = 1
       ELSE IF(wavmax(ns) <= 1.25_r8) THEN
          indxsl = 2
       ELSE IF(wavmax(ns) <= 2.38_r8) THEN
          indxsl = 3
       ELSE IF(wavmax(ns) > 2.38_r8) THEN
          indxsl = 4
       END IF

       abarii = abari(indxsl)
       bbarii = bbari(indxsl)
       cbarii = cbari(indxsl)
       dbarii = dbari(indxsl)
       ebarii = ebari(indxsl)
       fbarii = fbari(indxsl)

       DO k=1,pver
          DO i=1,Nday

             ! note that optical properties for ice valid only
             ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
             IF (cldn(i,k) >= cldmin .AND. cldn(i,k) >= cldeps) THEN
                tmp1i = abarii + bbarii/MAX(13._r8,MIN(scalefactor*rei(i,k),130._r8))
                ice_tau(ns,i,k) = cicewp(i,k)*tmp1i
             ELSE
                ice_tau(ns,i,k) = 0.0_r8
             ENDIF

             tmp2i = 1._r8 - cbarii - dbarii*MIN(MAX(13._r8,scalefactor*rei(i,k)),130._r8)
             tmp3i = fbarii*MIN(MAX(13._r8,scalefactor*rei(i,k)),130._r8)
             ! Do not let single scatter albedo be 1.  Delta-eddington solution
             ! for non-conservative case has different analytic form from solution
             ! for conservative case, and raddedmx is written for non-conservative case.
             ice_tau_w(ns,i,k) = ice_tau(ns,i,k) * MIN(tmp2i,.999999_r8)
             g = ebarii + tmp3i
             ice_tau_w_g(ns,i,k) = ice_tau_w(ns,i,k) * g
             ice_tau_w_f(ns,i,k) = ice_tau_w(ns,i,k) * g * g

          END DO ! End do i=1,Nday
       END DO    ! End do k=1,pver
    END DO ! nswbands

  END SUBROUTINE ec_ice_optics_sw


  !==============================================================================
  ! Private methods
  !==============================================================================


  SUBROUTINE slingo_liq_optics_sw(ncol,pcols,pver,cldn,rel,iclwp,ql,pdel,&
       liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp)

    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pver
    INTEGER, INTENT(IN   ) :: pcols
    REAL(r8),INTENT(IN   ) :: cldn(pcols,pver)
    REAL(r8),INTENT(IN   ) :: rel (pcols,pver)
    REAL(r8),INTENT(IN   ) :: iclwp(pcols,pver)
    REAL(r8),INTENT(IN   ) :: ql   (pcols,pver)
    REAL(r8),INTENT(IN   ) :: pdel  (pcols,pver)
    REAL(r8),INTENT(out) :: liq_tau    (nswbands,pcols,pver) ! extinction optical depth
    REAL(r8),INTENT(out) :: liq_tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
    REAL(r8),INTENT(out) :: liq_tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
    REAL(r8),INTENT(out) :: liq_tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w
    LOGICAL, INTENT(in ) :: oldliqwp

    REAL(r8), DIMENSION(pcols,pver) :: cliqwp
    REAL(r8), DIMENSION(nswbands)   :: wavmin
    REAL(r8), DIMENSION(nswbands)   :: wavmax
    REAL(R8), PARAMETER :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2

    ! Minimum cloud amount (as a fraction of the grid-box area) to 
    ! distinguish from clear sky
    REAL(r8), PARAMETER :: cldmin = 1.0e-80_r8

    ! Decimal precision of cloud amount (0 -> preserve full resolution;
    ! 10^-n -> preserve n digits of cloud amount)
    REAL(r8), PARAMETER :: cldeps = 0.0_r8

    ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
    ! Parameterization for the Shortwave Properties of Water Clouds' JAS
    ! vol. 46 may 1989 pp 1419-1427)
    REAL(r8) :: abarl(4) = &  ! A coefficient for extinction optical depth
         (/ 2.817e-02_r8, 2.682e-02_r8,2.264e-02_r8,1.281e-02_r8/)
    REAL(r8) :: bbarl(4) = &  ! B coefficient for extinction optical depth
         (/ 1.305_r8    , 1.346_r8    ,1.454_r8    ,1.641_r8    /)
    REAL(r8) :: cbarl(4) = &  ! C coefficient for single scat albedo
         (/-5.62e-08_r8 ,-6.94e-06_r8 ,4.64e-04_r8 ,0.201_r8    /)
    REAL(r8) :: dbarl(4) = &  ! D coefficient for single  scat albedo
         (/ 1.63e-07_r8 , 2.35e-05_r8 ,1.24e-03_r8 ,7.56e-03_r8 /)
    REAL(r8) :: ebarl(4) = &  ! E coefficient for asymmetry parameter
         (/ 0.829_r8    , 0.794_r8    ,0.754_r8    ,0.826_r8    /)
    REAL(r8) :: fbarl(4) = &  ! F coefficient for asymmetry parameter
         (/ 2.482e-03_r8, 4.226e-03_r8,6.560e-03_r8,4.353e-03_r8/)

    REAL(r8) :: abarli        ! A coefficient for current spectral band
    REAL(r8) :: bbarli        ! B coefficient for current spectral band
    REAL(r8) :: cbarli        ! C coefficient for current spectral band
    REAL(r8) :: dbarli        ! D coefficient for current spectral band
    REAL(r8) :: ebarli        ! E coefficient for current spectral band
    REAL(r8) :: fbarli        ! F coefficient for current spectral band

    ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
    ! greater than 20 micro-meters

    INTEGER :: ns, i, k, indxsl, Nday
    REAL(r8) :: tmp1l, tmp2l, tmp3l, g
!    REAL(r8) :: kext(pcols,pver)
    liq_tau=0.0_r8;liq_tau_w  =0.0_r8;liq_tau_w_g=0.0_r8;liq_tau_w_f=0.0_r8
    cliqwp=0.0_r8;wavmin=0.0_r8;wavmax=0.0_r8;
    Nday = ncol
    !itim = pbuf_old_tim_idx()
    !call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/))
    !call pbuf_get_field(pbuf, rel_idx, rel)

    IF (oldliqwp) THEN
       DO k=1,pver
          DO i = 1,Nday
             cliqwp(i,k) = ql(i,k)*pdel(i,k)/(SHR_CONST_G*MAX(0.01_r8,cldn(i,k)))
          END DO
       END DO
    ELSE
       !if (iclwp_idx<=0) then 
       !  call endrun('slingo_liq_optics_sw: oldliqwp must be set to true since ICLWP was not found in pbuf')
       !endif
       ! The following is the eventual target specification for in cloud liquid water path.
       !call pbuf_get_field(pbuf, iclwp_idx, tmpptr)
       cliqwp = iclwp
    ENDIF

    CALL get_sw_spectral_boundaries(wavmin,wavmax,'microns')

    DO ns = 1, nswbands
       ! Set index for cloud particle properties based on the wavelength,
       ! according to A. Slingo (1989) equations 1-3:
       ! Use index 1 (0.25 to 0.69 micrometers) for visible
       ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
       ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
       ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
       IF(wavmax(ns) <= 0.7_r8) THEN
          indxsl = 1
       ELSE IF(wavmax(ns) <= 1.25_r8) THEN
          indxsl = 2
       ELSE IF(wavmax(ns) <= 2.38_r8) THEN
          indxsl = 3
       ELSE IF(wavmax(ns) > 2.38_r8) THEN
          indxsl = 4
       END IF

       ! Set cloud extinction optical depth, single scatter albedo,
       ! asymmetry parameter, and forward scattered fraction:
       abarli = abarl(indxsl)
       bbarli = bbarl(indxsl)
       cbarli = cbarl(indxsl)
       dbarli = dbarl(indxsl)
       ebarli = ebarl(indxsl)
       fbarli = fbarl(indxsl)

       DO k=1,pver
          DO i=1,Nday

             ! note that optical properties for liquid valid only
             ! in range of 4.2 > rel > 16 micron (Slingo 89)
             IF (cldn(i,k) >= cldmin .AND. cldn(i,k) >= cldeps) THEN
                tmp1l = abarli + bbarli/MIN(MAX(4.2_r8,rel(i,k)),16._r8)
                liq_tau(ns,i,k) = 1000._r8*cliqwp(i,k)*tmp1l
             ELSE
                liq_tau(ns,i,k) = 0.0_r8
             ENDIF

             tmp2l = 1._r8 - cbarli - dbarli*MIN(MAX(4.2_r8,rel(i,k)),16._r8)
             tmp3l = fbarli*MIN(MAX(4.2_r8,rel(i,k)),16._r8)
             ! Do not let single scatter albedo be 1.  Delta-eddington solution
             ! for non-conservative case has different analytic form from solution
             ! for conservative case, and raddedmx is written for non-conservative case.
             liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * MIN(tmp2l,.999999_r8)
             g = ebarli + tmp3l
             liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g
             liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g

          END DO ! End do i=1,Nday
       END DO    ! End do k=1,pver
    END DO ! nswbands

    !call outfld('CL_OD_SW_OLD',liq_tau(idx_sw_diag,:,:), pcols, lchnk)
    !call outfld('REL_OLD',rel(:,:), pcols, lchnk)
    !call outfld('CLWPTH_OLD',cliqwp(:,:), pcols, lchnk)
    !call outfld('KEXT_OLD',kext(:,:), pcols, lchnk)


  END SUBROUTINE slingo_liq_optics_sw

  !==============================================================================

  SUBROUTINE get_liquid_optics_sw(ncol,pcols,pver,lamc,pgam,iclwpth, &
       tau, tau_w, tau_w_g, tau_w_f)
    !   type(physics_state), intent(in)   :: state
    !   type(physics_buffer_desc),pointer :: pbuf(:)
    INTEGER, INTENT(IN ) :: ncol
    INTEGER, INTENT(IN ) :: pcols
    INTEGER, INTENT(IN ) :: pver
    REAL(r8),INTENT(IN ) :: lamc   (pcols,pver) 
    REAL(r8),INTENT(IN ) :: pgam   (pcols,pver) 
    REAL(r8),INTENT(IN ) :: iclwpth(pcols,pver) 
    REAL(r8),INTENT(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
    REAL(r8),INTENT(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
    REAL(r8),INTENT(out) :: tau_w_g(nswbands,pcols,pver) ! asymetry parameter * tau * w
    REAL(r8),INTENT(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

!    REAL(r8), DIMENSION(pcols,pver) :: kext
    INTEGER i,k


    !   call pbuf_get_field(pbuf, i_lambda,  lamc   )
    !   call pbuf_get_field(pbuf, i_mu,      pgam   )
    !   call pbuf_get_field(pbuf, i_iclwp,   iclwpth)

    DO k = 1,pver
       DO i = 1,ncol
          IF(lamc(i,k) > 0._r8) THEN ! This seems to be clue from microphysics of no cloud
             CALL gam_liquid_sw(iclwpth(i,k), lamc(i,k), pgam(i,k), &
                  tau(1:nswbands,i,k), tau_w(1:nswbands,i,k), tau_w_g(1:nswbands,i,k), tau_w_f(1:nswbands,i,k))
          ELSE
             tau(1:nswbands,i,k)     = 0._r8
             tau_w(1:nswbands,i,k)   = 0._r8
             tau_w_g(1:nswbands,i,k) = 0._r8
             tau_w_f(1:nswbands,i,k) = 0._r8
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE get_liquid_optics_sw


  !==============================================================================

  SUBROUTINE gam_liquid_sw(clwptn, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f)
    REAL(r8), INTENT(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
    REAL(r8), INTENT(in) :: lamc   ! prognosed value of lambda for cloud
    REAL(r8), INTENT(in) :: pgam   ! prognosed value of mu for cloud
    REAL(r8), INTENT(out) :: tau(1:nswbands), tau_w(1:nswbands), tau_w_f(1:nswbands), tau_w_g(1:nswbands)
    ! for interpolating into mu/lambda
    INTEGER :: imu, kmu
    INTEGER :: ilambda, klambda
    INTEGER :: swband ! sw band index
    REAL(r8) :: ext, ssa, asm,wmu,  onemwmu,wlambda,onemwlambda,   lambdaplus, lambdaminus

    IF (clwptn < 1.e-80_r8) THEN
       tau     = 0._r8
       tau_w   = 0._r8
       tau_w_g = 0._r8
       tau_w_f = 0._r8
       RETURN
    ENDIF

    IF (pgam < g_mu(1) .OR. pgam > g_mu(nmu)) THEN
       WRITE(iulog,*)'pgam from prognostic cldwat2m',pgam
       WRITE(iulog,*)'g_mu from file',g_mu
       STOP 'pgam exceeds limits'
       !call endrun ('pgam exceeds limits')
    ENDIF
    DO imu = 1, nmu
       kmu = imu
       IF (g_mu(kmu) > pgam) EXIT
    ENDDO
    wmu = (g_mu(kmu) - pgam)/(g_mu(kmu) - g_mu(kmu-1))
    onemwmu = 1._r8 - wmu

    DO ilambda = 1, nlambda
       klambda = ilambda
       IF (wmu*g_lambda(kmu-1,ilambda) + onemwmu*g_lambda(kmu,ilambda) < lamc) EXIT
    ENDDO
    IF (klambda <= 1 .OR. klambda > nlambda)STOP 'lamc  exceeds limits'! call endrun('lamc  exceeds limits')
    lambdaplus = wmu*g_lambda(kmu-1,klambda  ) + onemwmu*g_lambda(kmu,klambda  )
    lambdaminus= wmu*g_lambda(kmu-1,klambda-1) + onemwmu*g_lambda(kmu,klambda-1)
    wlambda = (lambdaplus - lamc) / (lambdaplus - lambdaminus)
    onemwlambda = 1._r8 - wlambda

    DO swband = 1, nswbands
       ext =     wlambda*    wmu*ext_sw_liq(kmu-1,klambda-1,swband) + &
            onemwlambda*    wmu*ext_sw_liq(kmu-1,klambda  ,swband) + &
            wlambda*onemwmu*ext_sw_liq(kmu  ,klambda-1,swband) + &
            onemwlambda*onemwmu*ext_sw_liq(kmu  ,klambda  ,swband)
       ! probably should interpolate ext*ssa
       ssa =     wlambda*    wmu*ssa_sw_liq(kmu-1,klambda-1,swband) + &
            onemwlambda*    wmu*ssa_sw_liq(kmu-1,klambda  ,swband) + &
            wlambda*onemwmu*ssa_sw_liq(kmu  ,klambda-1,swband) + &
            onemwlambda*onemwmu*ssa_sw_liq(kmu  ,klambda  ,swband)
       ! probably should interpolate ext*ssa*asm
       asm =     wlambda*    wmu*asm_sw_liq(kmu-1,klambda-1,swband) + &
            onemwlambda*    wmu*asm_sw_liq(kmu-1,klambda  ,swband) + &
            wlambda*onemwmu*asm_sw_liq(kmu  ,klambda-1,swband) + &
            onemwlambda*onemwmu*asm_sw_liq(kmu  ,klambda  ,swband)
       ! compute radiative properties
       tau(swband) = clwptn * ext
       tau_w(swband) = clwptn * ext * ssa
       tau_w_g(swband) = clwptn * ext * ssa * asm
       tau_w_f(swband) = clwptn * ext * ssa * asm * asm
    ENDDO

    RETURN
  END SUBROUTINE gam_liquid_sw

  !==============================================================================


  !==============================================================================

  SUBROUTINE get_snow_optics_sw   (ncol,pcols,pver,dei, iciwpth, tau, tau_w, tau_w_g, tau_w_f)
    INTEGER, INTENT(IN   ) :: ncol
    INTEGER, INTENT(IN   ) :: pcols
    INTEGER, INTENT(IN   ) :: pver

    REAL(r8),INTENT(in ) :: dei     (pcols,pver) 
    REAL(r8),INTENT(IN ) :: iciwpth(pcols,pver) 
    REAL(r8),INTENT(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
    REAL(r8),INTENT(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
    REAL(r8),INTENT(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
    REAL(r8),INTENT(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

    REAL(r8) :: dlimited ! d limited to dmin,dmax range

    INTEGER  :: i,k,i_d_grid,k_d_eff,i_swband
    REAL(r8) :: wd, onemwd, ext, ssa, asm


    ! temporary code to support diagnostics of snow radiation
    !   call pbuf_get_field(pbuf, i_icswp, iciwpth)
    !   call pbuf_get_field(pbuf, i_des,   dei)
    ! temporary code to support diagnostics of snow radiation

    DO i = 1,ncol
       DO k = 1,pver
          dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
          IF( iciwpth(i,k) < 1.e-80_r8 .OR. dlimited .EQ. 0._r8) THEN
             ! if ice water path is too small, OD := 0
             tau    (:,i,k) = 0._r8
             tau_w  (:,i,k) = 0._r8
             tau_w_g(:,i,k) = 0._r8
             tau_w_f(:,i,k) = 0._r8
          ELSE 
             !if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
             !write(iulog,*) 'dei from prognostic cldwat2m',dei(i,k)
             !write(iulog,*) 'grid values of deff ice from optics file',g_d_eff
             !call endrun ('deff of ice exceeds limits')
             !endif
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             IF (dlimited <= g_d_eff(1)) THEN
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             ELSEIF (dlimited >= g_d_eff(n_g_d)) THEN
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             ELSE
                DO i_d_grid = 1, n_g_d
                   k_d_eff = i_d_grid
                   IF(g_d_eff(i_d_grid) > dlimited) EXIT
                ENDDO
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             ENDIF
             ! interpolate into grid and extract radiative properties
             DO i_swband = 1, nswbands
                ext = wd*ext_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*ext_sw_ice(k_d_eff  ,i_swband) 
                ssa = wd*ssa_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*ssa_sw_ice(k_d_eff  ,i_swband) 
                asm = wd*asm_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*asm_sw_ice(k_d_eff  ,i_swband) 
                tau    (i_swband,i,k)=iciwpth(i,k) * ext
                tau_w  (i_swband,i,k)=iciwpth(i,k) * ext * ssa
                tau_w_g(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm
                tau_w_f(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm * asm
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_snow_optics_sw


  !------------------------------------------------------------------------------
  SUBROUTINE get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)
    ! provide spectral boundaries of each shortwave band

    REAL(r8), INTENT(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
    CHARACTER(*), INTENT(in) :: units ! requested units
    low_boundaries=0.0_r8;high_boundaries=0.0_r8
 
    SELECT CASE (units)
    CASE ('inv_cm','cm^-1','cm-1')
       low_boundaries = wavenum_low
       high_boundaries = wavenum_high
    CASE('m','meter','meters')
       low_boundaries = 1.e-2_r8/wavenum_high
       high_boundaries = 1.e-2_r8/wavenum_low
    CASE('nm','nanometer','nanometers')
       low_boundaries = 1.e7_r8/wavenum_high
       high_boundaries = 1.e7_r8/wavenum_low
    CASE('um','micrometer','micrometers','micron','microns')
       low_boundaries = 1.e4_r8/wavenum_high
       high_boundaries = 1.e4_r8/wavenum_low
    CASE('cm','centimeter','centimeters')
       low_boundaries  = 1._r8/wavenum_high
       high_boundaries = 1._r8/wavenum_low
    CASE default
       ! call endrun('rad_constants.F90: spectral units not acceptable'//units)
       STOP 'rad_constants.F90: spectral units not acceptable'
    END SELECT

  END SUBROUTINE get_sw_spectral_boundaries


  !...................................
  SUBROUTINE cldinit(dtpg,NLAY,si_in,sl_in)
    !-----------------------------------
    REAL(KIND=r8)   , INTENT(in) ::  dtpg
    INTEGER, INTENT(in) :: NLAY !   NLAY            : vertical layer number                             !
    REAL(KIND=r8) , INTENT(in) :: si_in(NLAY+1)
    REAL(KIND=r8) , INTENT(in) :: sl_in(NLAY)
    !
    !--- Parameters & data statement for local calculations
    !
    REAL(KIND=r8)    :: dtph
    REAL(KIND=r8)    :: bbfr
    INTEGER :: kl
    INTEGER :: k,i
    ALLOCATE(si(NLAY+1))
    ALLOCATE(sl(NLAY))
    sl=sl_in
    si=si_in
    !  ---  compute llyr - the top of bl cld and is topmost non cld(low) layer
    !       for stratiform (at or above lowest 0.1 of the atmosphere)

    IF (iflip == 0) THEN      ! data from toa to sfc
       kl = NLAY
       DO k = NLAY+1, 2, -1
          kl = k
          IF (si(k) < 0.9e0_r8) EXIT
       ENDDO
       llyr = kl + 1
    ELSE                      ! data from sfc to top
       kl = 2
       DO k = 1, NLAY
          kl = k
          IF (si(k) < 0.9e0_r8) EXIT
       ENDDO
       llyr = kl - 1
    ENDIF                     ! end_if_iflip

    CALL gpvsl()
    CALL gpvsi()
    CALL gpvs()
    !-------------------------------------------------------------------------------
    ! ABSTRACT:
    !   * Reads various microphysical lookup tables used in COLUMN_MICRO
    !   * Lookup tables were created "offline" and are read in during execution
    !   * Creates lookup tables for saturation vapor pressure w/r/t water & ice
    !-------------------------------------------------------------------------------
    !   SUBROUTINES:
    !     MY_GROWTH_RATES - lookup table for growth of nucleated ice
    !
    !------------------------------------------------------------------------
    !  *************  Parameters used in ETA model -- Not used in Global Model *****
    !
    !--- DPHD, DLMD are delta latitude and longitude at the model (NOT geodetic) equator
    !    => "DX" is the hypotenuse of the model zonal & meridional grid increments.
    !
    !     DX=111.*(DPHD**2+DLMD**2)**.5         ! Resolution at MODEL equator (km)
    !     DX=MIN(100., MAX(5., DX) )
    !
    !--- Assume the following functional relationship for key constants that
    !    depend on grid resolution from DXmin (5 km) to DXmax (100 km) resolution:
    !
    !     DXmin=5.
    !     DXmax=100.
    !     DX=MIN(DXmax, MAX(DXmin, DX) )
    !
    !--- EXtune determines the degree to which the coefficients change with resolution.
    !    The larger EXtune is, the more sensitive the parameter.
    !
    !     EXtune=1.

    !
    !--- FXtune ==> F(DX) is the grid-resolution tuning parameter (from 0 to 1)
    !
    !     FXtune=((DXmax-DX)/(DXmax-DXmin))**EXtune
    !     FXtune=MAX(0., MIN(1., FXtune))
    !
    !--- Calculate grid-averaged RH for the onset of condensation (RHgrd) based on
    !    simple ***ASSUMED*** (user-specified) values at DXmax and at DXmin.
    !
    !     RH_DXmax=.90              !-- 90% RH at DXmax=100 km
    !     RH_DXmin=.98              !-- 98% RH at DXmin=5 km
    !
    !--- Note that RHgrd is right now fixed throughout the domain!!
    !
    !     RHgrd=RH_DXmax+(RH_DXmin-RH_DXmax)*FXtune
    !   ********************************************************************************
    !
    CALL ICE_LOOKUP  ()                 ! Lookup tables for ice
    CALL RAIN_LOOKUP ()                 ! Lookup tables for rain

    ABFR=-0.66_r8
    BBFR=100.0_r8
    CBFR=20.0_r8*con_pi*con_pi*BBFR*RHOL*1.E-21_r8
    !
    !--- QAUT0 is the threshold cloud content for autoconversion to rain
    !      needed for droplets to reach a diameter of 20 microns (following
    !      Manton and Cotton, 1977; Banta and Hanson, 1987, JCAM).  It is
    !      **STRONGLY** affected by the assumed droplet number concentrations
    !     XNCW!  For example, QAUT0=1.2567, 0.8378, or 0.4189 g/m**3 for
    !     droplet number concentrations of 300, 200, and 100 cm**-3, respectively.
    !
    !--- Calculate grid-averaged XNCW based on simple ***ASSUMED*** (user-specified)
    !    values at DXmax and at DXmin.
    !
    !     XNCW_DXmax=50.E6          !--  50 /cm**3 at DXmax=100 km
    !     XNCW_DXmin=200.E6         !-- 200 /cm**3 at DXmin=5 km
    !
    !--- Note that XNCW is right now fixed throughout the domain!!
    !
    !     XNCW=XNCW_DXmax+(XNCW_DXmin-XNCW_DXmax)*FXtune
    !
    !     QAUT0=con_pi*RHOL*XNCW*(20.E-6)**3/6.
    QAUTx=con_pi*RHOL*1.0E6_r8*(20.E-6_r8)**3/6.0_r8
    !
    !--- Based on rain lookup tables for mean diameters from 0.05 to 0.45 mm
    !    * Four different functional relationships of mean drop diameter as
    !      a function of rain rate (RR), derived based on simple fits to
    !      mass-weighted fall speeds of rain as functions of mean diameter
    !      from the lookup tables.
    !
    RR_DRmin=N0r0*RRATE(MDRmin)     ! RR for mean drop diameter of .05 mm
    RR_DR1=N0r0*RRATE(MDR1)         ! RR for mean drop diameter of .10 mm
    RR_DR2=N0r0*RRATE(MDR2)         ! RR for mean drop diameter of .20 mm
    RR_DR3=N0r0*RRATE(MDR3)         ! RR for mean drop diameter of .32 mm
    RR_DRmax=N0r0*RRATE(MDRmax)     ! RR for mean drop diameter of .45 mm
    !
    RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
    RQR_DR1=N0r0*MASSR(MDR1)        ! Rain content for mean drop diameter of .10 mm
    RQR_DR2=N0r0*MASSR(MDR2)        ! Rain content for mean drop diameter of .20 mm
    RQR_DR3=N0r0*MASSR(MDR3)        ! Rain content for mean drop diameter of .32 mm
    RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
    C_N0r0=con_pi*RHOL*N0r0
    CN0r0=1.E6_r8/C_N0r0**0.25_r8
    CN0r_DMRmin=1.0_r8/(con_pi*RHOL*DMRmin**4)
    CN0r_DMRmax=1.0_r8/(con_pi*RHOL*DMRmax**4)
    !
    !
    !     Find out what microphysics time step should be
    !
    mic_step = MAX(1, INT(dtpg/600.0_r8+0.5_r8))
    !     mic_step = max(1, int(dtpg/300.0_r8+0.5_r8))
    dtph     = dtpg / mic_step
    !      if (mype == 0) print *,' DTPG=',DTPG,' mic_step=',mic_step        &
    !     &,                ' dtph=',dtph
    !
    !--- Calculates coefficients for growth rates of ice nucleated in water
    !    saturated conditions, scaled by physics time step (lookup table)
    !
    CALL MY_GROWTH_RATES (DTPH)
    !
    !--- CIACW is used in calculating riming rates
    !      The assumed effective collection efficiency of cloud water rimed onto
    !      ice is =0.5_r8 below:
    !
    !Moor CIACW=DTPH*0.25_r8*con_pi*0.5_r8*(1.E5_r8)**C1   ! commented on 20050422
    !      ice is =0.1_r8 below:
    CIACW=DTPH*0.25_r8*con_pi*0.1_r8*(1.E5_r8)**C1
    !     CIACW = 0.0_r8      ! Brad's suggestion 20040614
    !
    !--- CIACR is used in calculating freezing of rain colliding with large ice
    !      The assumed collection efficiency is 1.0_r8
    !
    CIACR=con_pi*DTPH
    !
    !--- CRACW is used in calculating collection of cloud water by rain (an
    !      assumed collection efficiency of 1.0_r8)
    !
    !Moor CRACW=DTPH*0.25_r8*con_pi*1.0_r8                 ! commented on 20050422
    !
    !      assumed collection efficiency of 0.1_r8)
    CRACW=DTPH*0.25_r8*con_pi*0.1_r8
    !     CRACW = 0.0_r8      ! Brad's suggestion 20040614
    !
    ESW0=FPVSL(con_t0c)           ! Saturation vapor pressure at 0C
    RFmax=1.1_r8**Nrime          ! Maximum rime factor allowed
    !
    !------------------------------------------------------------------------
    !--------------- Constants passed through argument list -----------------
    !------------------------------------------------------------------------
    !
    !--- Important parameters for self collection (autoconversion) of
    !    cloud water to rain.
    !
    !--- CRAUT is proportional to the rate that cloud water is converted by
    !      self collection to rain (autoconversion rate)
    !
    CRAUT=1.0_r8-EXP(-1.E-3_r8*DTPH)
    !
    !     IF (MYPE == 0)
    !    & WRITE(6,"(A, A,F6.2,A, A,F5.4, A,F7.3,A, A,F6.2,A, A,F5.3,A)")
    !    &   'KEY MICROPHYSICAL PARAMETERS FOR '
    !    &  ,'DX=',DX,' KM:'
    !    &  ,'   FXtune=',FXtune
    !    &  ,'   RHgrd=',100.*RHgrd,' %'
    !    &  ,'   NCW=',1.E-6*XNCW,' /cm**3'
    !    &  ,'   QAUT0=',1.E3*QAUT0,' g/kg'
    !
    !--- For calculating snow optical depths by considering bulk density of
    !      snow based on emails from Q. Fu (6/27-28/01), where optical
    !      depth (T) = 1.5*SWP/(Reff*DENS), SWP is snow water path, Reff
    !      is effective radius, and DENS is the bulk density of snow.
    !
    !    SWP (kg/m**2)=(1.E-3 kg/g)*SWPrad, SWPrad in g/m**2 used in radiation
    !    T = 1.5*1.E3*SWPrad/(Reff*DENS)
    !
    !    See derivation for MASSI(INDEXS), note equal to RHO*QSNOW/NSNOW
    !
    !    SDENS=1.5e3/DENS, DENS=MASSI(INDEXS)/[con_pi*(1.E-6*INDEXS)**3]
    !
    DO I=MDImin,MDImax
       !MoorthiSDENS(I)=con_pi*1.5E-15*FLOAT(I*I*I)/MASSI(I)
       SDENS(I)=con_pi*1.0E-15_r8*FLOAT(I*I*I)/MASSI(I)
    ENDDO
    !
    !-----------------------------------------------------------------------
    !

    !...................................
  END SUBROUTINE cldinit

  !-----------------------------------
  SUBROUTINE progcld1(&
                                !...................................
       
                                !  ---  inputs:
       plyr       , &!real (kind=r8), intent(in) :: plyr (:,:): model layer mean pressure in mb (100Pa)  
       plvl       , &!real (kind=r8), intent(in) :: plvl  (:,:): model level pressure in mb (100Pa)       
       tlyr       , &!real (kind=r8), intent(in) :: tlyr (:,:): model layer mean temperature in k      
       qlyr       , &! real (kind=r8), intent(in) :: qlyr (:,:): layer specific humidity in gm/gm          
       qstl       , &!real (kind=r8), intent(in) :: qstl (:,:): layer saturate humidity in gm/gm                  !
       rhly       , &!real (kind=r8), intent(in) :: rhly (:,:): layer relative humidity (=qlyr/qstl) 
       clw        , &!real (kind=r8), intent(in) :: clw  (:,:): layer cloud condensate amount      
       xlat       , &!real (kind=r8), intent(in) :: xlat(:): grid latitude in radians               
       slmsk      , &!real (kind=r8), intent(in) :: slmsk(:): sea/land mask array (sea:0,land:1,sea-ice:2)
       IX         , &!integer,  intent(in) :: IX : horizontal dimention                 
       NLAY       , &!integer,  intent(in) :: NLAY : vertical layer/level dimensions          
       iflip      , &!integer,  intent(in) :: iflip: control flag for in/out vertical indexing 
       iovr       , &!integer,  intent(in) :: iovr : control flag for cloud overlap            
       sashal     , &!logical, intent(in) :: sashal
       crick_proof, &!logical, intent(in) :: crick_proof
       ccnorm     , &!logical, intent(in) :: ccnorm
                                !  ---  outputs:
       clouds     , &!real (kind=r8),intent(out) :: clouds(:,:,:) : cloud profiles        
       clds       , &!real (kind=r8),   intent(out) :: clds(:,:): fraction of clouds for low, mid, hi, tot, bl      !
       mtop       , &!integer,             ,   intent(out) :: mtop(:,:) : vertical indices for low, mid, hi cloud tops      !
       mbot         )!integer,             ,   intent(out) :: mbot(:,:) : vertical indices for low, mid, hi cloud bases     !

    ! =================   subprogram documentation block   ================ !
    !                                                                       !
    ! subprogram:    progcld1    computes cloud related quantities using    !
    !   zhao/moorthi's prognostic cloud microphysics scheme.                !
    !                                                                       !
    ! abstract:  this program computes cloud fractions from cloud           !
    !   condensates, calculates liquid/ice cloud droplet effective radius,  !
    !   and computes the low, mid, high, total and boundary layer cloud     !
    !   fractions and the vertical indices of low, mid, and high cloud      !
    !   top and base.  the three vertical cloud domains are set up in the   !
    !   initial subroutine "cldinit".                                       !
    !                                                                       !
    ! program history log:                                                  !
    !      11-xx-1992   y.h., k.a.c, a.k. - cloud parameterization          !
    !         'cldjms' patterned after slingo and slingo's work (jgr,       !
    !         1992), stratiform clouds are allowed in any layer except      !
    !         the surface and upper stratosphere. the relative humidity     !
    !         criterion may cery in different model layers.                 !
    !      10-25-1995   kenneth campana   - tuned cloud rh curves           !
    !         rh-cld relation from tables created using mitchell-hahn       !
    !         tuning technique on airforce rtneph observations.             !
    !      11-02-1995   kenneth campana   - the bl relationships used       !
    !         below llyr, except in marine stratus regions.                 !
    !      04-11-1996   kenneth campana   - save bl cld amt in cld(,5)      !
    !      12-29-1998   s. moorthi        - prognostic cloud method         !
    !      04-15-2003   yu-tai hou        - rewritten in frotran 90         !
    !         modulized form, seperate prognostic and diagnostic methods    !
    !         into two packages.                                            !
    !                                                                       !
    ! usage:         call progcld1                                          !
    !                                                                       !
    ! subprograms called:   gethml                                          !
    !                                                                       !
    ! attributes:                                                           !
    !   language:   fortran 90                                              !
    !   machine:    ibm-sp, sgi                                             !
    !                                                                       !
    !                                                                       !
    !  ====================  defination of variables  ====================  !
    !                                                                       !
    ! input variables:                                                      !
    !   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
    !   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
    !   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
    !   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
    !   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
    !   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
    !   clw   (IX,NLAY) : layer cloud condensate amount                     !
    !   xlat  (IX)      : grid latitude in radians                          !
    !   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
    !   IX              : horizontal dimention                              !
    !   NLAY,NLP1       : vertical layer/level dimensions                   !
    !   iflip           : control flag for in/out vertical indexing         !
    !                     =0: index from toa to surface                     !
    !                     =1: index from surface to toa                     !
    !   iovr            : control flag for cloud overlap                    !
    !                     =0 random overlapping clouds                      !
    !                     =1 max/ran overlapping clouds                     !
    !                                                                       !
    ! output variables:                                                     !
    !   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
    !      clouds(:,:,1) - layer total cloud fraction                       !
    !      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
    !      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
    !      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
    !      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
    !      clouds(:,:,6) - layer rain drop water path         not assigned  !
    !      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
    !  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
    !      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
    !  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
    !   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
    !   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
    !   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
    !                                                                       !
    !  ====================    end of description    =====================  !
    !
    IMPLICIT NONE

    !  ---  inputs
    INTEGER,  INTENT(in) :: IX
    INTEGER,  INTENT(in) :: NLAY
    INTEGER,  INTENT(in) :: iflip
    INTEGER,  INTENT(in) :: iovr

    REAL (kind=r8), INTENT(in) :: plvl (1:IX,1:NLAY+1) 
    REAL (kind=r8), INTENT(in) :: plyr (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: tlyr (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: qlyr (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: qstl (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: rhly (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: clw  (1:IX,1:NLAY) 

    REAL (kind=r8), INTENT(in) :: xlat  (1:IX)  
    REAL (kind=r8), INTENT(in) :: slmsk (1:IX)  
    LOGICAL, INTENT(in) :: sashal
    LOGICAL, INTENT(in) :: crick_proof
    LOGICAL, INTENT(in) :: ccnorm

    !  ---  outputs
    REAL (kind=r8),INTENT(out) :: clouds(1:IX,1:NLAY,1:NF_CLDS)

    REAL (kind=r8),   INTENT(out) :: clds(1:IX,1:5)

    INTEGER              ,   INTENT(out) :: mtop(1:IX,1:3)
    INTEGER              ,   INTENT(out) :: mbot(1:IX,1:3)

    !  ---  local variables:
    REAL (kind=r8) :: cldtot(IX,NLAY)
    REAL (kind=r8) :: cldcnv(IX,NLAY)
    REAL (kind=r8) :: cwp(IX,NLAY)
    REAL (kind=r8) :: cip(IX,NLAY)
    REAL (kind=r8) :: crp(IX,NLAY)
    REAL (kind=r8) :: csp(IX,NLAY)
    REAL (kind=r8) :: rew(IX,NLAY)
    REAL (kind=r8) :: rei(IX,NLAY)
    REAL (kind=r8) :: res(IX,NLAY)
    REAL (kind=r8) :: rer(IX,NLAY)
    REAL (kind=r8) :: delp(IX,NLAY)
    REAL (kind=r8) :: tem2d(IX,NLAY)
    REAL (kind=r8) :: clwf(IX,NLAY)

    REAL (kind=r8) :: ptop1(IX,4)

    REAL (kind=r8) :: clwmin
    REAL (kind=r8) :: clwm
    REAL (kind=r8) :: clwt
    REAL (kind=r8) :: onemrh
    REAL (kind=r8) :: value
    REAL (kind=r8) ::  tem1
    REAL (kind=r8) :: tem2
    REAL (kind=r8) :: tem3


    INTEGER :: i, k, id, nf


    !  ---  local variables:
    clds = 0.0_r8
    mtop = 0
    mbot = 0
    cldtot = 0.0_r8
    cldcnv = 0.0_r8
    cwp = 0.0_r8
    cip = 0.0_r8
    crp = 0.0_r8
    csp = 0.0_r8
    rew = 0.0_r8
    rei = 0.0_r8
    res = 0.0_r8
    rer = 0.0_r8
    delp = 0.0_r8
    tem2d = 0.0_r8
    clwf = 0.0_r8
    ptop1 = 0.0_r8

    !
    !===> ... begin here
    !
    DO nf=1,nf_clds
       DO k=1,nlay
          DO i=1,ix
             clouds(i,k,nf) = 0.0_r8
          ENDDO
       ENDDO
    ENDDO
    !     clouds(:,:,:) = 0.0

    DO k = 1, NLAY
       DO i = 1, IX
          cldtot(i,k) = 0.0_r8
          cldcnv(i,k) = 0.0_r8
          cwp   (i,k) = 0.0_r8
          cip   (i,k) = 0.0_r8
          crp   (i,k) = 0.0_r8
          csp   (i,k) = 0.0_r8
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = MIN( 1.0_r8, MAX( 0.0_r8, (con_ttp-tlyr(i,k))*0.05_r8 ) )
          clwf(i,k)   = 0.0_r8
       ENDDO
    ENDDO
    !
    IF (crick_proof) THEN
       DO i = 1, IX
          clwf(i,1)    = 0.75_r8*clw(i,1)    + 0.25_r8*clw(i,2)
          clwf(i,nlay) = 0.75_r8*clw(i,nlay) + 0.25_r8*clw(i,nlay-1)
       ENDDO
       DO k = 2, NLAY-1
          DO i = 1, IX
             clwf(i,K) = 0.25_r8*clw(i,k-1) + 0.5_r8*clw(i,k) + 0.25_r8*clw(i,k+1)
          ENDDO
       ENDDO
    ELSE
       DO k = 1, NLAY
          DO i = 1, IX
             clwf(i,k) = clw(i,k)
          ENDDO
       ENDDO
    ENDIF

    !  ---  find top pressure for each cloud domain for given latitude
    !       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
    !  ---  i=1,2 are low-lat (<45 degree) and pole regions)

    DO id = 1, 4
       tem1 = ptopc(id,2) - ptopc(id,1)

       DO i =1, IX
          ptop1(i,id) = ptopc(id,1) +                                   &
               tem1 * MAX( 0.0_r8, 4.0_r8*ABS(xlat(i))/con_pi-1.0_r8 )
       ENDDO
    ENDDO

    !  ---  compute liquid/ice condensate path in g/m**2


    IF (iflip == 0) THEN             ! input data from toa to sfc
       DO k = 1, NLAY
          DO i = 1, IX
             delp(i,k) = plvl(i,k+1) - plvl(i,k)
             clwt     = MAX(0.0_r8, clwf(i,k)) * gfac * delp(i,k)
             cip(i,k) = clwt * tem2d(i,k)
             cwp(i,k) = clwt - cip(i,k)
          ENDDO
       ENDDO
    ELSE                             ! input data from sfc to toa
       DO k = 1, NLAY
          DO i = 1, IX
             delp(i,k) = plvl(i,k) - plvl(i,k+1)
             clwt     = MAX(0.0_r8, clwf(i,k)) * gfac * delp(i,k)
             cip(i,k) = clwt * tem2d(i,k)
             cwp(i,k) = clwt - cip(i,k)
          ENDDO
       ENDDO
    ENDIF                            ! end_if_iflip

    !  ---  effective liquid cloud droplet radius over land

    DO k = 1, NLAY
       DO i = 1, IX
          IF (NINT(slmsk(i)) == 1) THEN
             rew(i,k) = 5.0_r8 + 5.0_r8 * tem2d(i,k)
          ENDIF
       ENDDO
    ENDDO

    !  ---  layer cloud fraction

    IF (iflip == 0) THEN                 ! input data from toa to sfc

       clwmin = 0.0_r8
       IF (.NOT. sashal) THEN
          DO k = NLAY, 1, -1
             DO i = 1, IX
                clwt = 1.0e-6_r8 * (plyr(i,k)*0.001_r8)
                !           clwt = 2.0e-6 * (plyr(i,k)*0.001_r8)

                IF (clwf(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   tem1  = MIN(MAX(SQRT(SQRT(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   tem1  = 2000.0_r8 / tem1
                   !             tem1  = 1000.0_r8 / tem1

                   value = MAX( MIN( tem1*(clwf(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )

                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO k = NLAY, 1, -1
             DO i = 1, IX
                clwt = 1.0e-6_r8 * (plyr(i,k)*0.001_r8)
                !           clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)

                IF (clwf(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   !             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   !             tem1  = 2000.0_r8 / tem1

                   tem1  = MIN(MAX((onemrh*qstl(i,k))**0.49_r8,0.0001_r8),1.0_r8)  !jhan
                   tem1  = 100.0_r8 / tem1
                   !
                   !             tem1  = 2000.0_r8 / tem1
                   !             tem1  = 1000.0_r8 / tem1
                   !

                   value = MAX( MIN( tem1*(clwf(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )
                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ELSE                                 ! input data from sfc to toa

       clwmin = 0.0_r8
       IF (.NOT. sashal) THEN
          DO k = 1, NLAY
             DO i = 1, IX
                clwt = 1.0e-6_r8 * (plyr(i,k)*0.001_r8)
                !           clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)

                IF (clwf(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   tem1  = MIN(MAX(SQRT(SQRT(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   tem1  = 2000.0_r8 / tem1

                   !             tem1  = 1000.0_r8 / tem1

                   value = MAX( MIN( tem1*(clwf(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )

                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO k = 1, NLAY
             DO i = 1, IX
                clwt = 1.0e-6_r8 * (plyr(i,k)*0.001_r8)
                !           clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)

                IF (clwf(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   !             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   !             tem1  = 2000.0_r8 / tem1

                   tem1  = MIN(MAX((onemrh*qstl(i,k))**0.49_r8,0.0001_r8),1.0_r8)  !jhan
                   tem1  = 100.0_r8 / tem1
                   !
                   !             tem1  = 2000.0_r8 / tem1
                   !             tem1  = 1000.0_r8 / tem1
                   !
                   value = MAX( MIN( tem1*(clwf(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )
                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ENDIF                                ! end_if_flip

    DO k = 1, NLAY
       DO i = 1, IX
          IF (cldtot(i,k) < climit_cld) THEN
             cldtot(i,k) = 0.0_r8
             cwp(i,k)    = 0.0_r8
             cip(i,k)    = 0.0_r8
             crp(i,k)    = 0.0_r8
             csp(i,k)    = 0.0_r8
          ENDIF
       ENDDO
    ENDDO
    !     where (cldtot < climit_cld)
    !       cldtot = 0.0_r8
    !       cwp    = 0.0_r8
    !       cip    = 0.0_r8
    !       crp    = 0.0_r8
    !       csp    = 0.0_r8
    !     endwhere
    !
    IF (ccnorm) THEN
       DO k = 1, NLAY
          DO i = 1, IX
             IF (cldtot(i,k) >= climit_cld) THEN
                tem1 = 1.0_r8 / MAX(climit2, cldtot(i,k))
                cwp(i,k) = cwp(i,k) * tem1
                cip(i,k) = cip(i,k) * tem1
                crp(i,k) = crp(i,k) * tem1
                csp(i,k) = csp(i,k) * tem1
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    !  ---  effective ice cloud droplet radius

    DO k = 1, NLAY
       DO i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          IF (cip(i,k) > 0.0_r8) THEN
             tem3 = gord * cip(i,k) * ( plyr(i,k) / delp(i,k) )          &
                  / (tlyr(i,k) * (1.0_r8 + con_fvirt * qlyr(i,k)))

             IF (tem2 < -50.0_r8) THEN
                rei(i,k) = (1250.0_r8/9.917_r8) * tem3 ** 0.109_r8
             ELSEIF (tem2 < -40.0_r8) THEN
                rei(i,k) = (1250.0_r8/9.337_r8) * tem3 ** 0.08_r8
             ELSEIF (tem2 < -30.0_r8) THEN
                rei(i,k) = (1250.0_r8/9.208_r8) * tem3 ** 0.055_r8
             ELSE
                rei(i,k) = (1250.0_r8/9.387_r8) * tem3 ** 0.031_r8
             ENDIF
             rei(i,k)   = MAX(20.0_r8, MIN(rei(i,k), 300.0_r8))
             !           rei(i,k)   = max(10.0_r8, min(rei(i,k), 100.0_r8))
          ENDIF
       ENDDO
    ENDDO

    !
    DO k = 1, NLAY
       DO i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          !         clouds(i,k,6) = 0.0_r8
          clouds(i,k,7) = rer(i,k)
          !         clouds(i,k,8) = 0.0_r8
          clouds(i,k,9) = rei(i,k)
       ENDDO
    ENDDO


    !  ---  compute low, mid, high, total, and boundary layer cloud fractions
    !       and clouds top/bottom layer indices for low, mid, and high clouds.
    !       The three cloud domain boundaries are defined by ptopc.  The cloud
    !       overlapping method is defined by control flag 'iovr', which is
    !  ---  also used by the lw and sw radiation programs.

    CALL gethml (&
                                !  ---  inputs:
         plyr    (1:IX,1:NLAY) , &!real   , intent(in) :: plyr(:,:) model layer mean pressure in mb (100Pa)
         ptop1   (1:IX,1:4)    , &!real   , intent(in) :: ptop1(:,:)pressure limits of cloud domain interfaces  in mb (100Pa)
         cldtot  (1:IX,1:NLAY) , &!real   , intent(in) :: cldtot(:,:)total or straiform cloud profile in fraction
         cldcnv  (1:IX,1:NLAY) , &!real   , intent(in) :: cldcnv(:,:)convective cloud (for diagnostic scheme only)
         IX                    , &!integer, intent(in) :: IXhorizontal dimention
         NLAY                  , &!integer, intent(in) :: NLAY    vertical layer dimensions
         iflip                 , &!integer, intent(in) :: iflip  control flag for in/out vertical indexing
         iovr                  , &!integer, intent(in) :: iovr  control flag for cloud overlap
                                !  ---  outputs:
         clds    (1:IX,1:5)    , &!real   , intent(out) :: clds(:,:)fraction of clouds for low, mid, hi, tot, bl   
         mtop    (1:IX,1:3)    , &!integer, intent(out) :: mtop(:,:)vertical indices for low, mid, hi cloud tops
         mbot    (1:IX,1:3)      )!integer, intent(out) :: mbot(:,:)vertical indices for low, mid, hi cloud bases


    !
    RETURN
    !...................................
  END SUBROUTINE progcld1
  !-----------------------------------


  !-----------------------------------
  SUBROUTINE progcld2( &
                                !...................................
       
                                !  ---  inputs:
       plyr         , &!REAL (kind=r8), INTENT(in) :: plyr(:,:)
       plvl         , &!REAL (kind=r8), INTENT(in) :: plvl(:,:)
       tlyr         , &!REAL (kind=r8), INTENT(in) :: tlyr(:,:)
       qlyr         , &!REAL (kind=r8), INTENT(in) :: qlyr(:,:)
       qstl         , &!REAL (kind=r8), INTENT(in) :: qstl(:,:)
       rhly         , &!REAL (kind=r8), INTENT(in) :: rhly(:,:)
       clw          , &!REAL (kind=r8), INTENT(in) :: clw(:,:)
       xlat         , &!REAL (kind=r8), INTENT(in) :: xlat(:)
       f_ice        , &!REAL (kind=r8), INTENT(in) :: f_ice(:,:)
       f_rain       , &!REAL (kind=r8), INTENT(in) :: f_rain(:,:)
       f_liq       , &!REAL (kind=r8), INTENT(in) :: f_liq(:,:)
       r_rime       , &!REAL (kind=r8), INTENT(in) :: r_rime(:,:)
       flgmin       , &!REAL (kind=r8), INTENT(in) :: flgmin(:)
       IX           , &!INTEGER,  INTENT(in) :: IX
       NLAY         , &!INTEGER,  INTENT(in) :: NLAY
       iflip        , &!INTEGER,  INTENT(in) :: iflip
       iovr         , &!INTEGER,  INTENT(in) :: iovr
       sashal       , &!LOGICAL       , INTENT(in) :: sashal
       norad_precip , &!LOGICAL       , INTENT(in) :: norad_precip
       crick_proof  , &!LOGICAL       , INTENT(in) :: crick_proof
       ccnorm       , &!LOGICAL       , INTENT(in) :: ccnorm
                                !  ---  outputs:   
       clouds       , &!REAL (kind=r8), INTENT(out) :: clouds(:,:,:)
       clds         , &!REAL (kind=r8), INTENT(out) :: clds(:,:)
       mtop         , &!INTEGER              , INTENT(out) :: mtop(:,:)
       mbot           )!INTEGER              , INTENT(out) :: mbot(:,:)

    ! =================   subprogram documentation block   ================ !
    !                                                                       !
    ! subprogram:    progcld2    computes cloud related quantities using    !
    !   ferrier's prognostic cloud microphysics scheme.                     !
    !                                                                       !
    ! abstract:  this program computes cloud fractions from cloud           !
    !   condensates, calculates liquid/ice cloud droplet effective radius,  !
    !   and computes the low, mid, high, total and boundary layer cloud     !
    !   fractions and the vertical indices of low, mid, and high cloud      !
    !   top and base.  the three vertical cloud domains are set up in the   !
    !   initial subroutine "cldinit".                                       !
    !                                                                       !
    ! program history log:                                                  !
    !        -  -       brad ferrier      - original development            !
    !        -  -2003   s. moorthi        - adapted to ncep gfs model       !
    !      05-05-2004   yu-tai hou        - rewritten as a separated        !
    !                   program in the cloud module.                        !
    !                                                                       !
    ! usage:         call progcld2                                          !
    !                                                                       !
    ! subprograms called:   gethml                                          !
    !                                                                       !
    ! attributes:                                                           !
    !   language:   fortran 90                                              !
    !   machine:    ibm-sp, sgi                                             !
    !                                                                       !
    !                                                                       !
    !  ====================  defination of variables  ====================  !
    !                                                                       !
    ! input variables:                                                      !
    !   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
    !   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
    !   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
    !   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
    !   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
    !   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
    !   clw   (IX,NLAY) : layer cloud condensate amount                     !
    !   f_ice (IX,NLAY) : fraction of layer cloud ice  (ferrier micro-phys) !
    !   f_rain(IX,NLAY) : fraction of layer rain water (ferrier micro-phys) !
    !   r_rime(IX,NLAY) : mass ratio of total ice to unrimed ice (>=1)      !
    !   xlat  (IX)      : grid latitude in radians                          !
    !   IX              : horizontal dimention                              !
    !   NLAY,NLP1       : vertical layer/level dimensions                   !
    !   iflip           : control flag for in/out vertical indexing         !
    !                     =0: index from toa to surface                     !
    !                     =1: index from surface to toa                     !
    !   iovr            : control flag for cloud overlap                    !
    !                     =0 random overlapping clouds                      !
    !                     =1 max/ran overlapping clouds                     !
    !                                                                       !
    ! output variables:                                                     !
    !   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
    !      clouds(:,:,1) - layer total cloud fraction                       !
    !      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
    !      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
    !      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
    !      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
    !      clouds(:,:,6) - layer rain drop water path         (g/m**2)      !
    !      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
    !  *** clouds(:,:,8) - layer snow flake water path        (g/m**2)      !
    !      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
    !  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
    !   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
    !   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
    !   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
    !                                                                       !
    !  ====================    end of description    =====================  !
    !
    IMPLICIT NONE
    !  ---  inputs
    INTEGER,  INTENT(in) :: IX
    INTEGER,  INTENT(in) :: NLAY
    INTEGER,  INTENT(in) :: iflip
    INTEGER,  INTENT(in) :: iovr

    REAL (kind=r8), INTENT(in) :: plvl  (1:IX,1:NLAY+1) 
    REAL (kind=r8), INTENT(in) :: plyr  (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: tlyr  (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: qlyr  (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: qstl  (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: rhly  (1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: clw   (1:IX,1:NLAY) 
    REAL (kind=r8), INTENT(in) :: f_ice (1:IX,1:NLAY) 
    REAL (kind=r8), INTENT(in) :: f_rain(1:IX,1:NLAY)
    REAL (kind=r8), INTENT(in) :: f_liq (1:IX,1:NLAY) 
    REAL (kind=r8), INTENT(in) :: r_rime(1:IX,1:NLAY)

    REAL (kind=r8), INTENT(in) :: xlat  (1:IX)  
    LOGICAL       , INTENT(in) :: sashal
    LOGICAL       , INTENT(in) :: norad_precip
    LOGICAL       , INTENT(in) :: crick_proof
    LOGICAL       , INTENT(in) :: ccnorm
    REAL (kind=r8), INTENT(in) :: flgmin (1:IX)

    !  ---  outputs
    REAL (kind=r8), INTENT(out) :: clouds(1:IX,1:NLAY,1:NF_CLDS)

    REAL (kind=r8), INTENT(out) :: clds  (1:IX,1:5)

    INTEGER              , INTENT(out) :: mtop(1:IX,1:3) 
    INTEGER              , INTENT(out) :: mbot(1:IX,1:3)

    !  ---  local variables:
    REAL (kind=r8) :: cldtot(IX,NLAY)
    REAL (kind=r8) :: cldcnv(IX,NLAY)
    REAL (kind=r8) :: cwp(IX,NLAY)
    REAL (kind=r8) :: cip(IX,NLAY)
    REAL (kind=r8) :: crp(IX,NLAY)
    REAL (kind=r8) :: csp(IX,NLAY)
    REAL (kind=r8) :: rew(IX,NLAY)
    REAL (kind=r8) :: rei(IX,NLAY)
    REAL (kind=r8) :: res(IX,NLAY)
    REAL (kind=r8) :: rer(IX,NLAY)
    REAL (kind=r8) :: tem2d(IX,NLAY)
    REAL (kind=r8) :: clw2(IX,NLAY)
    REAL (kind=r8) :: qcwat(IX,NLAY)
    REAL (kind=r8) :: qcice(IX,NLAY)
    REAL (kind=r8) :: qrain(IX,NLAY)
    REAL (kind=r8) :: fcice(IX,NLAY)
    REAL (kind=r8) :: frain(IX,NLAY)
    REAL (kind=r8) :: fcliq(IX,NLAY)

    REAL (kind=r8) :: rrime(IX,NLAY)
    REAL (kind=r8) :: rsden(IX,NLAY)
    REAL (kind=r8) :: clwf(IX,NLAY)

    REAL (kind=r8) :: ptop1(IX,4)

    REAL (kind=r8) :: clwmin
    REAL (kind=r8) :: clwm
    REAL (kind=r8) :: clwt
    REAL (kind=r8) :: onemrh
    REAL (kind=r8) :: value
    REAL (kind=r8) :: tem1
    REAL (kind=r8) :: tem2
    REAL (kind=r8) :: tem3


    INTEGER :: i, k, id

    !  ---  constants
    REAL (kind=r8), PARAMETER :: EPSQ = 1.0e-12_r8


    clouds= 0.0_r8

    clds= 0.0_r8

    mtop= 0
    mbot= 0

    !  ---  local variables:
    tem2d= 0.0_r8
    clw2= 0.0_r8
    qcwat= 0.0_r8
    qcice= 0.0_r8
    qrain= 0.0_r8
    fcice= 0.0_r8
    frain= 0.0_r8
    fcliq= 0.0_r8

    rrime= 0.0_r8
    rsden= 0.0_r8
    clwf= 0.0_r8

    ptop1= 0.0_r8

    clwmin= 0.0_r8
    clwm= 0.0_r8
    clwt= 0.0_r8
    onemrh= 0.0_r8
    value= 0.0_r8
    tem1= 0.0_r8
    tem2= 0.0_r8
    tem3= 0.0_r8

    !
    !===> ... begin here
    !
    !     clouds(:,:,:) = 0.0_r8

    DO k = 1, NLAY
       DO i = 1, IX
          cldtot(i,k) = 0.0_r8
          cldcnv(i,k) = 0.0_r8
          cwp   (i,k) = 0.0_r8
          cip   (i,k) = 0.0_r8
          crp   (i,k) = 0.0_r8
          csp   (i,k) = 0.0_r8
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          fcice (i,k) = MAX(0.0_r8, MIN(1.0_r8, f_ice(i,k)))
          frain (i,k) = MAX(0.0_r8, MIN(1.0_r8, f_rain(i,k)))
          fcliq (i,k) = MAX(0.0_r8, MIN(1.0_r8, f_liq(i,k)))
          rrime (i,k) = MAX(1.0_r8, r_rime(i,k))
          tem2d (i,k) = tlyr(i,k) - con_t0c
       ENDDO
    ENDDO
    !
    IF (crick_proof) THEN
       DO i = 1, IX
          clwf(i,1)    = 0.75_r8*clw(i,1)    + 0.25_r8*clw(i,2)
          clwf(i,nlay) = 0.75_r8*clw(i,nlay) + 0.25_r8*clw(i,nlay-1)
       ENDDO
       DO k = 2, NLAY-1
          DO i = 1, IX
             clwf(i,K) = 0.25_r8*clw(i,k-1) + 0.5_r8*clw(i,k) + 0.25_r8*clw(i,k+1)
          ENDDO
       ENDDO
    ELSE
       DO k = 1, NLAY
          DO i = 1, IX
             clwf(i,k) = clw(i,k)
          ENDDO
       ENDDO
    ENDIF

    !  ---  find top pressure for each cloud domain for given latitude
    !       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
    !  ---  i=1,2 are low-lat (<45 degree) and pole regions)

    DO id = 1, 4
       tem1 = ptopc(id,2) - ptopc(id,1)

       DO i =1, IX
          tem2 = MAX( 0.0_r8, 4.0_r8*ABS(xlat(i))/con_pi-1.0_r8 )
          ptop1(i,id) = ptopc(id,1) + tem1*tem2
       ENDDO
    ENDDO

    !  ---  separate cloud condensate into liquid, ice, and rain types, and
    !       save the liquid+ice condensate in array clw2 for later
    !       calculation of cloud fraction

    DO k = 1, NLAY
       DO i = 1, IX
          IF (tem2d(i,k) > -40.0_r8) THEN
             !qcice(i,k) = clwf(i,k) * fcice(i,k)
             qcice(i,k) =  fcice(i,k)
             !tem1       = clwf(i,k) - qcice(i,k)
             !qrain(i,k) = tem1 * frain(i,k)
             qrain(i,k) = frain(i,k)

             !qcwat(i,k) = tem1 - qrain(i,k)
             qcwat(i,k) = fcliq(i,k)
             clw2 (i,k) = qcwat(i,k) + qcice(i,k)
          ELSE
             qcice(i,k) = clwf(i,k)
             qrain(i,k) = 0.0_r8
             qcwat(i,k) = 0.0_r8
             clw2 (i,k) = clwf(i,k)
          ENDIF
       ENDDO
    ENDDO

    CALL  rsipath2 ( &
                                !  ---  inputs:
         plyr  (1:IX,1:NLAY)   , & !real   , intent(in) :: plyr(:,:)model layer mean pressure in mb (100Pa)
         plvl  (1:IX,1:NLAY+1) , & !real   , intent(in) :: plvl(:,:)model level pressure in mb (100Pa)     
         tlyr  (1:IX,1:NLAY)   , & !real   , intent(in) :: tlyr(:,:)model layer mean temperature in k    
         qlyr  (1:IX,1:NLAY)   , & !real   , intent(in) :: qlyr(:,:)layer specific humidity in gm/gm         
         qcwat (1:IX,1:NLAY)   , & !real   , intent(in) :: qcwat(:,:)layer cloud liquid water condensate amount    
         qcice (1:IX,1:NLAY)   , & !real   , intent(in) :: qcice(:,:)layer cloud ice water condensate amount              
         qrain (1:IX,1:NLAY)   , & !real   , intent(in) :: qrain(:,:)layer rain drop water amount 
         rrime (1:IX,1:NLAY)   , & !real   , intent(in) :: rrime(:,:)mass ratio of total to unrimed ice ( >= 1 )    
         IX                    , & !integer, intent(in) :: IM         horizontal dimention                 
         NLAY                  , & !integer, intent(in) :: LEVS       vertical layer dimensions     
         iflip                 , & !integer, intent(in) :: iflip      control flag for in/out vertical indexing     
         flgmin(1:IX)          , & !real   , intent(in) :: flgmin(:)  Minimum large ice fraction                
                                !  ---  outputs:
         cwp   (1:IX,1:NLAY)   , & !real   , intent(out) :: cwatp(:,:) layer cloud liquid water path 
         cip   (1:IX,1:NLAY)   , & !real   , intent(out) :: cicep(:,:)layer cloud ice water path           
         crp   (1:IX,1:NLAY)   , & !real   , intent(out) :: rainp(:,:)layer rain water path                         
         csp   (1:IX,1:NLAY)   , & !real   , intent(out) :: snowp(:,:)layer snow water path             
         rew   (1:IX,1:NLAY)   , & !real   , intent(out) :: recwat(:,:)layer cloud eff radius for liqid water (micron) 
         rer   (1:IX,1:NLAY)   , & !real   , intent(out) :: rerain(:,:)layer rain water effective radius      (micron)
         res   (1:IX,1:NLAY)   , & !real   , intent(out) :: resnow(:,:)layer snow flake effective radius      (micron)
         rsden (1:IX,1:NLAY)    ) !real   , intent(out) :: snden(:,:)1/snow density                 


    IF (iflip == 0) THEN             ! input data from toa to sfc
       DO k = 1, NLAY
          DO i = 1, IX
             tem2d(i,k) = (con_g * plyr(i,k))                            &
                  &                 / (con_rd* (plvl(i,k+1) - plvl(i,k)))
          ENDDO
       ENDDO
    ELSE                             ! input data from sfc to toa
       DO k = 1, NLAY
          DO i = 1, IX
             tem2d(i,k) = (con_g * plyr(i,k))                            &
                  &                 / (con_rd* (plvl(i,k) - plvl(i,k+1)))
          ENDDO
       ENDDO
    ENDIF                            ! end_if_iflip

    !  ---  layer cloud fraction

    IF (iflip == 0) THEN                 ! input data from toa to sfc

       clwmin = 0.0_r8
       IF (.NOT. sashal) THEN
          DO k = NLAY, 1, -1
             DO i = 1, IX
                !           clwt = 1.0e-7 * (plyr(i,k)*0.001)
                !           clwt = 1.0e-6 * (plyr(i,k)*0.001)
                clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)
                !           clwt = 5.0e-6 * (plyr(i,k)*0.001)
                !           clwt = 5.0e-6

                IF (clw2(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   !             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001_r8),1.0_r8)
                   !             tem1  = 100.0_r8 / tem1

                   tem1  = MIN(MAX(SQRT(SQRT(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   tem1  = 2000.0_r8 / tem1
                   !             tem1  = 2400.0_r8 / tem1
                   !cnt          tem1  = 2500.0_r8 / tem1
                   !             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001_r8),1.0_r8)
                   !             tem1  = 2000.0_r8 / tem1
                   !             tem1  = 1000.0_r8 / tem1
                   !             tem1  = 100.0_r8 / tem1

                   value = MAX( MIN( tem1*(clw2(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )

                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO k = NLAY, 1, -1
             DO i = 1, IX
                !           clwt = 1.0e-6_r8 * (plyr(i,k)*0.001_r8)
                clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)

                IF (clw2(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   tem1  = MIN(MAX((onemrh*qstl(i,k))**0.49_r8,0.0001_r8),1.0_r8)    !jhan
                   tem1  = 100.0_r8 / tem1

                   !             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
                   !             tem1  = 2000.0 / tem1
                   !
                   !             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
                   !             tem1  = 2200.0 / tem1
                   !             tem1  = 2400.0 / tem1
                   !             tem1  = 2500.0 / tem1
                   !             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
                   !             tem1  = 2000.0 / tem1
                   !             tem1  = 1000.0 / tem1
                   !             tem1  = 100.0 / tem1
                   !

                   value = MAX( MIN( tem1*(clw2(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )

                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ELSE                                 ! input data from sfc to toa

       clwmin = 0.0e-6_r8
       IF (.NOT. sashal) THEN
          DO k = 1, NLAY
             DO i = 1, IX
                !           clwt = 1.0e-7 * (plyr(i,k)*0.001)
                !           clwt = 1.0e-6 * (plyr(i,k)*0.001)
                clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)
                !           clwt = 5.0e-6 * (plyr(i,k)*0.001)
                !           clwt = 5.0e-6

                IF (clw2(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   !             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
                   !             tem1  = 100.0 / tem1

                   tem1  = MIN(MAX(SQRT(SQRT(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   tem1  = 2000.0_r8 / tem1
                   !             tem1  = 2400.0 / tem1
                   !cnt          tem1  = 2500.0 / tem1
                   !             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
                   !             tem1  = 2000.0 / tem1
                   !             tem1  = 1000.0 / tem1
                   !             tem1  = 100.0 / tem1

                   value = MAX( MIN( tem1*(clw2(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )

                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO k = 1, NLAY
             DO i = 1, IX
                !           clwt = 1.0e-6 * (plyr(i,k)*0.001)
                clwt = 2.0e-6_r8 * (plyr(i,k)*0.001_r8)

                IF (clw2(i,k) > clwt) THEN

                   onemrh= MAX( 1.e-10_r8, 1.0_r8-rhly(i,k) )
                   clwm  = clwmin / MAX( 0.01_r8, plyr(i,k)*0.001_r8 )

                   tem1  = MIN(MAX((onemrh*qstl(i,k))**0.49_r8,0.0001_r8),1.0_r8)   !jhan
                   tem1  = 100.0_r8 / tem1

                   !             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001_r8),1.0_r8)
                   !             tem1  = 2000.0_r8 / tem1
                   !
                   !             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
                   !             tem1  = 2200.0 / tem1
                   !             tem1  = 2400.0 / tem1
                   !             tem1  = 2500.0 / tem1
                   !             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
                   !             tem1  = 2000.0 / tem1
                   !             tem1  = 1000.0 / tem1
                   !             tem1  = 100.0 / tem1

                   value = MAX( MIN( tem1*(clw2(i,k)-clwm), 50.0_r8 ), 0.0_r8 )
                   tem2  = SQRT( SQRT(rhly(i,k)) )

                   cldtot(i,k) = MAX( tem2*(1.0_r8-EXP(-value)), 0.0_r8 )
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ENDIF                                ! end_if_flip

    DO k = 1, NLAY
       DO i = 1, IX
          IF (cldtot(i,k) < climit_cld) THEN
             cldtot(i,k) = 0.0_r8
             cwp(i,k)    = 0.0_r8
             cip(i,k)    = 0.0_r8
             crp(i,k)    = 0.0_r8
             csp(i,k)    = 0.0_r8
          ENDIF
       ENDDO
    ENDDO
    !     where (cldtot < climit_cld)
    !       cldtot = 0.0
    !       cwp    = 0.0
    !       cip    = 0.0
    !       crp    = 0.0
    !       csp    = 0.0
    !     endwhere 

    !     When norad_precip = .true. snow/rain has no impact on radiation
    IF (norad_precip) THEN
       DO k = 1, NLAY
          DO i = 1, IX
             crp(i,k) = 0.0_r8
             csp(i,k) = 0.0_r8
          ENDDO
       ENDDO
    ENDIF
    !
    IF (ccnorm) THEN
       DO k = 1, NLAY
          DO i = 1, IX
             IF (cldtot(i,k) >= climit_cld) THEN
                tem1 = 1.0_r8 / MAX(climit2, cldtot(i,k))
                cwp(i,k) = cwp(i,k) * tem1
                cip(i,k) = cip(i,k) * tem1
                crp(i,k) = crp(i,k) * tem1
                csp(i,k) = csp(i,k) * tem1
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    !  ---  effective ice cloud droplet radius

    DO k = 1, NLAY
       DO i = 1, IX
          tem1 = tlyr(i,k) - con_ttp
          tem2 = cip(i,k)

          IF (tem2 > 0.0_r8) THEN
             tem3 = tem2d(i,k) * tem2                                    &
                  &           / (tlyr(i,k) * (1.0_r8 + con_fvirt * qlyr(i,k)))

             IF (tem1 < -50.0_r8) THEN
                rei(i,k) = (1250.0_r8/9.917_r8) * tem3 ** 0.109_r8
             ELSEIF (tem1 < -40.0_r8) THEN
                rei(i,k) = (1250.0_r8/9.337_r8) * tem3 ** 0.08_r8
             ELSEIF (tem1 < -30.0_r8) THEN
                rei(i,k) = (1250.0_r8/9.208_r8) * tem3 ** 0.055_r8
             ELSE
                rei(i,k) = (1250.0_r8/9.387_r8) * tem3 ** 0.031_r8
             ENDIF

             !           if (lprnt .and. k == l) print *,' reiL=',rei(i,k),' icec=', &
             !    &        icec,' cip=',cip(i,k),' tem=',tem,' delt=',delt

             rei(i,k)   = MAX(10.0_r8, MIN(rei(i,k), 300.0_r8))
             !           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!!!!        rei(i,k)   = max(30.0, min(rei(i,k), 300.0))
             !           rei(i,k)   = max(50.0, min(rei(i,k), 300.0))
             !           rei(i,k)   = max(100.0, min(rei(i,k), 300.0))
          ENDIF
       ENDDO
    ENDDO
    !
    DO k = 1, NLAY
       DO i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
          !         clouds(i,k,8) = csp(i,k)     !ncar scheme
          clouds(i,k,8) = csp(i,k) * rsden(i,k)  !fu's scheme
          clouds(i,k,9) = rei(i,k)
       ENDDO
    ENDDO


    !  ---  compute low, mid, high, total, and boundary layer cloud fractions
    !       and clouds top/bottom layer indices for low, mid, and high clouds.
    !       The three cloud domain boundaries are defined by ptopc.  The cloud
    !       overlapping method is defined by control flag 'iovr', which is
    !  ---  also used by the lw and sw radiation programs.

    CALL gethml ( &
                                !  ---  inputs:
         plyr   (1:IX,1:NLAY) , &!real   , intent(in) :: plyr(:,:) model layer mean pressure in mb (100Pa)
         ptop1  (1:IX,1:4)    , &!real   , intent(in) :: ptop1(:,:)pressure limits of cloud domain interfaces  in mb (100Pa)
         cldtot (1:IX,1:NLAY) , &!real   , intent(in) :: cldtot(:,:)total or straiform cloud profile in fraction
         cldcnv (1:IX,1:NLAY) , &!real   , intent(in) :: cldcnv(:,:)convective cloud (for diagnostic scheme only)
         IX                   , &!integer, intent(in) :: IX      horizontal dimention               
         NLAY                 , &!integer, intent(in) :: NLAY    vertical layer dimensions            
         iflip                , &!integer, intent(in) :: iflip  control flag for in/out vertical indexing
         iovr                 , &!integer, intent(in) :: iovr  control flag for cloud overlap          
                                !  ---  outputs:
         clds   (1:IX,1:5)    , &!real   , intent(out) :: clds(:,:)fraction of clouds for low, mid, hi, tot, bl   
         mtop   (1:IX,1:3)    , &!integer, intent(out) :: mtop(:,:)vertical indices for low, mid, hi cloud tops        
         mbot   (1:IX,1:3)      )!integer, intent(out) :: mbot(:,:)vertical indices for low, mid, hi cloud bases


    !
    RETURN
    !...................................
  END SUBROUTINE progcld2
  !-----------------------------------

  !-----------------------------------                                    !
  SUBROUTINE gethml (&
                                !...................................                                    !
       
                                !  ---  inputs:
       plyr      , &!real   , intent(in) :: plyr(:,:) model layer mean pressure in mb (100Pa)
       ptop1     , &!real   , intent(in) :: ptop1(:,:)pressure limits of cloud domain interfaces  in mb (100Pa)
       cldtot    , &!real   , intent(in) :: cldtot(:,:)total or straiform cloud profile in fraction
       cldcnv    , &!real   , intent(in) :: cldcnv(:,:)convective cloud (for diagnostic scheme only)
       IX        , &!integer, intent(in) :: IX      horizontal dimention              
       NLAY      , &!integer, intent(in) :: NLAY    vertical layer dimensions            
       iflip     , &!integer, intent(in) :: iflip  control flag for in/out vertical indexing
       iovr      , &!integer, intent(in) :: iovr  control flag for cloud overlap         
                                !  ---  outputs:
       clds      , &!real   , intent(out) :: clds(:,:)fraction of clouds for low, mid, hi, tot, bl   
       mtop      , &!integer, intent(out) :: mtop(:,:)vertical indices for low, mid, hi cloud tops     
       mbot        )!integer, intent(out) :: mbot(:,:)vertical indices for low, mid, hi cloud bases

    !  ===================================================================  !
    !                                                                       !
    ! abstract: compute high, mid, low, total, and boundary cloud fractions !
    !   and cloud top/bottom layer indices for model diagnostic output.     !
    !   the three cloud domain boundaries are defined by ptopc.  the cloud  !
    !   overlapping method is defined by control flag 'iovr', which is also !
    !   used by lw and sw radiation programs.                               !
    !                                                                       !
    ! program history log:                                                  !
    !      04-29-2004   yu-tai hou        - separated to become individule  !
    !         subprogram to calculate averaged h,m,l,bl cloud amounts.      !
    !                                                                       !
    ! usage:         call gethml                                            !
    !                                                                       !
    ! subprograms called:  none                                             !
    !                                                                       !
    ! attributes:                                                           !
    !   language:   fortran 90                                              !
    !   machine:    ibm-sp, sgi                                             !
    !                                                                       !
    !                                                                       !
    !  ====================  defination of variables  ====================  !
    !                                                                       !
    ! input variables:                                                      !
    !   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
    !   ptop1 (IX,4)    : pressure limits of cloud domain interfaces        !
    !                     (sfc,low,mid,high) in mb (100Pa)                  !
    !   cldtot(IX,NLAY) : total or straiform cloud profile in fraction      !
    !   cldcnv(IX,NLAY) : convective cloud (for diagnostic scheme only)     !
    !   IX              : horizontal dimention                              !
    !   NLAY            : vertical layer dimensions                         !
    !   iflip           : control flag for in/out vertical indexing         !
    !                     =0: index from toa to surface                     !
    !                     =1: index from surface to toa                     !
    !   iovr            : control flag for cloud overlap                    !
    !                     =0 random overlapping clouds                      !
    !                     =1 max/ran overlapping clouds                     !
    !                                                                       !
    ! output variables:                                                     !
    !   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
    !   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
    !   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
    !                                                                       !
    !  ====================    end of description    =====================  !
    !
    IMPLICIT NONE!

    !  ---  inputs:
    INTEGER, INTENT(in) :: IX
    INTEGER, INTENT(in) :: NLAY
    INTEGER, INTENT(in) :: iflip
    INTEGER, INTENT(in) :: iovr

    REAL (kind=r8), INTENT(in) :: plyr   (1:IX,1:NLAY) !real   , intent(in) :: plyr(:,:) model layer mean pressure in mb (100Pa)
    REAL (kind=r8), INTENT(in) :: ptop1  (1:IX,1:4)    !real   , intent(in) :: ptop1(:,:)pressure limits of cloud domain interfaces  in mb (100Pa)
    REAL (kind=r8), INTENT(in) :: cldtot (1:IX,1:NLAY) !real   , intent(in) :: cldtot(:,:)total or straiform cloud profile in fraction
    REAL (kind=r8), INTENT(in) :: cldcnv (1:IX,1:NLAY) !real   , intent(in) :: cldcnv(:,:)convective cloud (for diagnostic scheme only)

    !  ---  outputs
    REAL (kind=r8), INTENT(out) :: clds  (1:IX,1:5)   !real   , intent(out) :: clds(:,:)fraction of clouds for low, mid, hi, tot, bl   

    INTEGER       , INTENT(out) :: mtop  (1:IX,1:3)!integer, intent(out) :: mtop(:,:)vertical indices for low, mid, hi cloud tops   
    INTEGER       , INTENT(out) :: mbot  (1:IX,1:3)!integer, intent(out) :: mbot(:,:)vertical indices for low, mid, hi cloud bases

    !  ---  local variables:
    REAL (kind=r8) :: cl1(IX)
    REAL (kind=r8) :: cl2(IX)

    REAL (kind=r8) :: pcur
    REAL (kind=r8) :: pnxt
    REAL (kind=r8) :: ccur
    REAL (kind=r8) :: cnxt

    INTEGER :: idom(IX)
    INTEGER :: kbt1(IX)
    INTEGER :: kth1(IX)
    INTEGER :: kbt2(IX)
    INTEGER :: kth2(IX)

    INTEGER :: i, k, id, id1, kstr, kend, kinc

    !
    !===> ... begin here
    !
    DO i = 1, IX
       clds(i,1) = 0.0_r8
       clds(i,2) = 0.0_r8
       clds(i,3) = 0.0_r8
       clds(i,4) = 0.0_r8
       clds(i,5) = 0.0_r8
       mtop(i,1) = 1
       mtop(i,2) = 1
       mtop(i,3) = 1
       mbot(i,1) = 1
       mbot(i,2) = 1
       mbot(i,3) = 1
       cl1 (i) = 1.0_r8
       cl2 (i) = 1.0_r8
    ENDDO

    !  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
    !       layer processed from surface and up

    IF (iflip == 0) THEN                      ! input data from toa to sfc
       kstr = NLAY
       kend = 1
       kinc = -1
    ELSE                                      ! input data from sfc to toa
       kstr = 1
       kend = NLAY
       kinc = 1
    ENDIF                                     ! end_if_iflip

    IF (iovr == 0) THEN                       ! random overlap

       DO k = kstr, kend, kinc
          DO i = 1, IX
             ccur = MIN( ovcst, MAX( cldtot(i,k), cldcnv(i,k) ))
             IF (ccur >= climit_cld) cl1(i) = cl1(i) * (1.0_r8 - ccur)
          ENDDO

          IF (k == llyr) THEN
             DO i = 1, IX
                clds(i,5) = 1.0_r8 - cl1(i)          ! save bl cloud
             ENDDO
          ENDIF
       ENDDO

       DO i = 1, IX
          clds(i,4) = 1.0_r8 - cl1(i)              ! save total cloud
       ENDDO

    ELSE                                      ! max/ran overlap

       DO k = kstr, kend, kinc
          DO i = 1, IX
             ccur = MIN( ovcst, MAX( cldtot(i,k), cldcnv(i,k) ))
             IF (ccur >= climit_cld) THEN             ! cloudy layer
                cl2(i) = MIN( cl2(i), (1.0_r8 - ccur) )
             ELSE                                ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0_r8
             ENDIF
          ENDDO

          IF (k == llyr) THEN
             DO i = 1, IX
                clds(i,5) = 1.0_r8 - cl1(i) * cl2(i) ! save bl cloud
             ENDDO
          ENDIF
       ENDDO

       DO i = 1, IX
          clds(i,4) = 1.0_r8 - cl1(i) * cl2(i)     ! save total cloud
       ENDDO

    ENDIF                                     ! end_if_iovr

    !  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
    !       layer processed from one layer below llyr and up
    !  ---  change! layer processed from surface to top, so low clouds will
    !       contains both bl and low clouds.

    IF (iflip == 0) THEN                      ! input data from toa to sfc

       DO i = 1, IX
          cl1 (i) = 0.0_r8
          cl2 (i) = 0.0_r8
          kbt1(i) = NLAY
          kbt2(i) = NLAY
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
       ENDDO

       !org    do k = llyr-1, 1, -1
       !PK DO k = NLAY, 1, -1
       DO k = llyr-1, 1, -1
          DO i = 1, IX
             id = idom(i)
             id1= id + 1

             pcur = plyr(i,k)
             ccur = MIN( ovcst, MAX( cldtot(i,k), cldcnv(i,k) ))

             IF (k > 1) THEN
                pnxt = plyr(i,k-1)
                cnxt = MIN( ovcst, MAX( cldtot(i,k-1), cldcnv(i,k-1) ))
             ELSE
                pnxt = -1.0_r8
                cnxt = 0.0_r8
             ENDIF

             IF (pcur < ptop1(i,id1)) THEN
                id = id + 1
                id1= id1 + 1
                idom(i) = id
             ENDIF

             IF (ccur >= climit_cld) THEN
                IF (kth2(i) == 0) kbt2(i) = k
                kth2(i) = kth2(i) + 1

                IF (iovr == 0) THEN
                   cl2(i) = cl2(i) + ccur - cl2(i)*ccur
                ELSE
                   cl2(i) = MAX( cl2(i), ccur )
                ENDIF

                IF (cnxt < climit_cld .OR. pnxt < ptop1(i,id1)) THEN
                   kbt1(i) = NINT( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
                        / (cl1(i) + cl2(i)) )
                   kth1(i) = NINT( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
                        / (cl1(i) + cl2(i)) )
                   cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                   kbt2(i) = 1
                   kth2(i) = 0
                   cl2 (i) = 0.0_r8
                ENDIF   ! end_if_cnxt_or_pnxt
             ENDIF     ! end_if_ccur

             IF (pnxt < ptop1(i,id1)) THEN
                clds(i,id) = cl1(i)
                mtop(i,id) = MIN( kbt1(i), kbt1(i)-kth1(i)+1 )
                mbot(i,id) = kbt1(i)

                cl1 (i) = 0.0_r8
                kbt1(i) = 1
                kth1(i) = 0
             ENDIF     ! end_if_pnxt

          ENDDO       ! end_do_i_loop
       ENDDO         ! end_do_k_loop

    ELSE                                      ! input data from sfc to toa

       DO i = 1, IX
          cl1 (i) = 0.0_r8
          cl2 (i) = 0.0_r8
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
       ENDDO

       !org    do k = llyr+1, NLAY
       !PK DO k = 1, NLAY
       DO k = llyr+1, NLAY
          DO i = 1, IX
             id = idom(i)
             id1= id + 1

             pcur = plyr(i,k)
             ccur = MIN( ovcst, MAX( cldtot(i,k), cldcnv(i,k) ))

             IF (k < NLAY) THEN
                pnxt = plyr(i,k+1)
                cnxt = MIN( ovcst, MAX( cldtot(i,k+1), cldcnv(i,k+1) ))
             ELSE
                pnxt = -1.0_r8
                cnxt = 0.0_r8
             ENDIF

             IF (pcur < ptop1(i,id1)) THEN
                id = id + 1
                id1= id1 + 1
                idom(i) = id
             ENDIF

             IF (ccur >= climit_cld) THEN
                IF (kth2(i) == 0) kbt2(i) = k
                kth2(i) = kth2(i) + 1

                IF (iovr == 0) THEN
                   cl2(i) = cl2(i) + ccur - cl2(i)*ccur
                ELSE
                   cl2(i) = MAX( cl2(i), ccur )
                ENDIF

                IF (cnxt < climit_cld .OR. pnxt < ptop1(i,id1)) THEN
                   kbt1(i) = NINT( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
                        / (cl1(i) + cl2(i)) )
                   kth1(i) = NINT( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
                        / (cl1(i) + cl2(i)) )
                   cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                   kbt2(i) = 1
                   kth2(i) = 0
                   cl2 (i) = 0.0_r8
                ENDIF     ! end_if_cnxt_or_pnxt
             ENDIF       ! end_if_ccur

             IF (pnxt < ptop1(i,id1)) THEN
                clds(i,id) = cl1(i)
                mtop(i,id) = MAX( kbt1(i), kbt1(i)+kth1(i)-1 )
                mbot(i,id) = kbt1(i)

                cl1 (i) = 0.0_r8
                kbt1(i) = 1
                kth1(i) = 0
             ENDIF     ! end_if_pnxt

          ENDDO       ! end_do_i_loop
       ENDDO         ! end_do_k_loop

    ENDIF                                     ! end_if_iflip

    !
    RETURN
    !...................................
  END SUBROUTINE gethml




  !-----------------------------------
  SUBROUTINE rsipath2( &
                                !...................................
       
                                !  ---  inputs:  
       plyr   , &!real   , intent(in) :: plyr(:,:)model layer mean pressure in mb (100Pa)
       plvl   , &!real   , intent(in) :: plvl(:,:)model level pressure in mb (100Pa)     
       tlyr   , &!real   , intent(in) :: tlyr(:,:)model layer mean temperature in k    
       qlyr   , &!real   , intent(in) :: qlyr(:,:)layer specific humidity in gm/gm    
       qcwat2 , &!real   , intent(in) :: qcwat(:,:)layer cloud liquid water condensate amount    
       qcice2 , &!real   , intent(in) :: qcice(:,:)layer cloud ice water condensate amount         
       qrain2 , &!real   , intent(in) :: qrain(:,:)layer rain drop water amount 
       rrime  , &!real   , intent(in) :: rrime(:,:)mass ratio of total to unrimed ice ( >= 1 )    
       IM     , &!integer, intent(in) :: IM         horizontal dimention              
       LEVS   , &!integer, intent(in) :: LEVS       vertical layer dimensions     
       iflip  , &!integer, intent(in) :: iflip      control flag for in/out vertical indexing     
       flgmin , &!real   , intent(in) :: flgmin(:)  Minimum large ice fraction               
                                !  ---  outputs:
       cwatp  , &!real   , intent(out) :: cwatp(:,:) layer cloud liquid water path 
       cicep  , &!real   , intent(out) :: cicep(:,:)layer cloud ice water path          
       rainp  , &!real   , intent(out) :: rainp(:,:)layer rain water path                     
       snowp  , &!real   , intent(out) :: snowp(:,:)layer snow water path         
       recwat , &!real   , intent(out) :: recwat(:,:)layer cloud eff radius for liqid water (micron) 
       rerain , &!real   , intent(out) :: rerain(:,:)layer rain water effective radius      (micron)
       resnow , &!real   , intent(out) :: resnow(:,:)layer snow flake effective radius      (micron)
       snden    )!real   , intent(out) :: snden(:,:)1/snow density            

    ! =================   subprogram documentation block   ================ !
    !                                                                       !
    ! abstract:  this program is a modified version of ferrier's original   !
    !   "rsipath" subprogram.  it computes layer's cloud liquid, ice, rain, !
    !   and snow water condensate path and the partical effective radius    !
    !   for liquid droplet, rain drop, and snow flake.                      !
    !                                                                       !
    !  ====================  defination of variables  ====================  !
    !                                                                       !
    ! input variables:                                                      !
    !   plyr  (IM,LEVS) : model layer mean pressure in mb (100Pa)           !
    !   plvl  (IM,LEVS+1):model level pressure in mb (100Pa)                !
    !   tlyr  (IM,LEVS) : model layer mean temperature in k                 !
    !   qlyr  (IM,LEVS) : layer specific humidity in gm/gm                  !
    !   qcwat (IM,LEVS) : layer cloud liquid water condensate amount        !
    !   qcice (IM,LEVS) : layer cloud ice water condensate amount           !
    !   qrain (IM,LEVS) : layer rain drop water amount                      !
    !   rrime (IM,LEVS) : mass ratio of total to unrimed ice ( >= 1 )       !
    !   IM              : horizontal dimention                              !
    !   LEVS            : vertical layer dimensions                         !
    !   iflip           : control flag for in/out vertical indexing         !
    !                     =0: index from toa to surface                     !
    !                     =1: index from surface to toa                     !
    !   flgmin          : Minimum large ice fraction                        !
    !   lprnt           : logical check print control flag                  !
    !                                                                       !
    ! output variables:                                                     !
    !   cwatp (IM,LEVS) : layer cloud liquid water path                     !
    !   cicep (IM,LEVS) : layer cloud ice water path                        !
    !   rainp (IM,LEVS) : layer rain water path                             !
    !   snowp (IM,LEVS) : layer snow water path                             !
    !   recwat(IM,LEVS) : layer cloud eff radius for liqid water (micron)   !
    !   rerain(IM,LEVS) : layer rain water effective radius      (micron)   !
    !   resnow(IM,LEVS) : layer snow flake effective radius      (micron)   !
    !   snden (IM,LEVS) : 1/snow density                                    !
    !                                                                       !
    !                                                                       !
    ! usage:     call rsipath2                                              !
    !                                                                       !
    ! subroutines called:  none                                             !
    !                                                                       !
    ! program history log:                                                  !
    !      xx-xx-2001   b. ferrier     - original program                   !
    !      xx-xx-2004   s. moorthi     - modified for use in gfs model      !
    !      05-20-2004   y. hou         - modified, added vertical index flag!
    !                     to reduce data flipping, and rearrange code to    !
    !                     be comformable with radiation part programs.      !
    !                                                                       !
    !  ====================    end of description    =====================  !
    !

    IMPLICIT NONE


    !  ---  inputs:
    INTEGER, INTENT(in) :: IM
    INTEGER, INTENT(in) :: LEVS

    REAL(KIND=r8) , INTENT(in) ::       plyr  (1:IM,1:LEVS)    !real   , intent(in) :: plyr(:,:)model layer mean pressure in mb (100Pa)
    REAL(KIND=r8) , INTENT(in) ::       plvl  (1:IM,1:LEVS+1)  !real   , intent(in) :: plvl(:,:)model level pressure in mb (100Pa)     
    REAL(KIND=r8) , INTENT(in) ::       tlyr  (1:IM,1:LEVS)    !real   , intent(in) :: tlyr(:,:)model layer mean temperature in k    
    REAL(KIND=r8) , INTENT(in) ::       qlyr  (1:IM,1:LEVS)    !real   , intent(in) :: qlyr(:,:)layer specific humidity in gm/gm 
    REAL(KIND=r8) , INTENT(in) ::       qcwat2(1:IM,1:LEVS)    !real   , intent(in) :: qcwat(:,:)layer cloud liquid water condensate amount    
    REAL(KIND=r8) , INTENT(in) ::       qcice2(1:IM,1:LEVS)    !real   , intent(in) :: qcice(:,:)layer cloud ice water condensate amount       
    REAL(KIND=r8) , INTENT(in) ::       qrain2(1:IM,1:LEVS)    !real   , intent(in) :: qrain(:,:)layer rain drop water amount 
    REAL(KIND=r8) , INTENT(in) ::       rrime (1:IM,1:LEVS)    !real   , intent(in) :: rrime(:,:)mass ratio of total to unrimed ice ( >= 1 )   

    INTEGER, INTENT(in) :: iflip
    REAL(KIND=r8) ,   INTENT(in) :: flgmin(1:IM)           !real   , intent(in) :: flgmin(:)  Minimum large ice fraction                

    !  ---  output:
    REAL(KIND=r8),  INTENT(out) :: cwatp (1:IM,1:LEVS)    !real   , intent(out) :: cwatp(:,:) layer cloud liquid water path 
    REAL(KIND=r8),  INTENT(out) :: cicep (1:IM,1:LEVS)    !real   , intent(out) :: cicep(:,:)layer cloud ice water path
    REAL(KIND=r8),  INTENT(out) :: rainp (1:IM,1:LEVS)    !real   , intent(out) :: rainp(:,:)layer rain water path
    REAL(KIND=r8),  INTENT(out) :: snowp (1:IM,1:LEVS)    !real   , intent(out) :: snowp(:,:)layer snow water path
    REAL(KIND=r8),  INTENT(out) :: recwat(1:IM,1:LEVS)    !real   , intent(out) :: recwat(:,:)layer cloud eff radius for liqid water (micron) 
    REAL(KIND=r8),  INTENT(out) :: rerain(1:IM,1:LEVS)    !real   , intent(out) :: rerain(:,:)layer rain water effective radius      (micron)
    REAL(KIND=r8),  INTENT(out) :: resnow(1:IM,1:LEVS)    !real   , intent(out) :: resnow(:,:)layer snow flake effective radius      (micron)
    REAL(KIND=r8),  INTENT(out) :: snden (1:IM,1:LEVS)    !real   , intent(out) :: snden(:,:)1/snow density

    !  ---  locals:
    REAL(KIND=r8)    :: qcwat(im,LEVS)
    REAL(KIND=r8)    :: qcice(im,LEVS)
    REAL(KIND=r8)    :: qrain(im,LEVS)
    REAL(KIND=r8)    :: recw1
    REAL(KIND=r8)    :: dsnow
    REAL(KIND=r8)    :: qsnow
    REAL(KIND=r8)    :: qqcice
    REAL(KIND=r8)    :: flarge
    REAL(KIND=r8)    :: xsimass
    REAL(KIND=r8)    :: pfac
    REAL(KIND=r8)    :: nlice
    REAL(KIND=r8)    :: xli
    REAL(KIND=r8)    :: nlimax
    REAL(KIND=r8)    :: dum
    REAL(KIND=r8)    :: tem
    REAL(KIND=r8)    :: rho
    REAL(KIND=r8)    :: cpath
    REAL(KIND=r8)    :: totcnd
    REAL(KIND=r8)    :: tc

    INTEGER :: i
    INTEGER :: k
    INTEGER :: indexs
    INTEGER :: ksfc
    INTEGER :: k1
    !  ---  constant parameter:
    REAL(KIND=r8), PARAMETER :: CEXP= 1.0_r8/3.0_r8

    cwatp = 0.0_r8
    cicep = 0.0_r8
    rainp = 0.0_r8
    snowp = 0.0_r8
    recwat = 0.0_r8
    rerain = 0.0_r8
    resnow = 0.0_r8
    snden = 0.0_r8

    !  ---  locals:
    qcwat = 0.0_r8
    qcice = 0.0_r8
    qrain = 0.0_r8
    recw1 = 0.0_r8
    dsnow = 0.0_r8
    qsnow = 0.0_r8
    qqcice = 0.0_r8
    flarge = 0.0_r8
    xsimass = 0.0_r8
    pfac = 0.0_r8
    nlice = 0.0_r8
    xli = 0.0_r8
    nlimax = 0.0_r8
    dum = 0.0_r8
    tem = 0.0_r8
    rho = 0.0_r8
    cpath = 0.0_r8
    totcnd = 0.0_r8
    tc = 0.0_r8

    !
    !===>  ...  begin here
    !
    recw1 = 620.3505_r8 / TNW**CEXP         ! cloud droplet effective radius

    DO k = 1, LEVS
       DO i = 1, IM
          totcnd = qcwat2(i,k) + qcice2(i,k) + qrain2(i,k)
          qcwat(i,k)= qcwat2(i,k)/MAX(totcnd,EPSQ)
          qcice(i,k)= qcice2(i,k)/MAX(totcnd,EPSQ)
          qrain(i,k)= qrain2(i,k)/MAX(totcnd,EPSQ)

          !--- hydrometeor's optical path
          cwatp(i,k) = 0.0_r8
          cicep(i,k) = 0.0_r8
          rainp(i,k) = 0.0_r8
          snowp(i,k) = 0.0_r8
          snden(i,k) = 0.0_r8
          !--- hydrometeor's effective radius
          recwat(i,k) = RECWmin
          rerain(i,k) = RERAINmin
          resnow(i,k) = RESNOWmin
       ENDDO
    ENDDO

    !  ---  set up pressure elated arrays, convert unit from mb to cb (10Pa)
    !       cause the rest part uses cb in computation

    IF (iflip == 0) THEN        ! data from toa to sfc
       ksfc = levs + 1
       k1   = 0
    ELSE                        ! data from sfc to top
       ksfc = 1
       k1   = 1
    ENDIF                       ! end_if_iflip
    !
    DO k = 1, LEVS
       DO i = 1, IM
          totcnd = qcwat(i,k) + qcice(i,k) + qrain(i,k)
          qsnow = 0.0_r8
          IF(totcnd > EPSQ) THEN

             !  ---  air density (rho), model mass thickness (cpath), temperature in c (tc)

             rho   = 0.1_r8 * plyr(i,k)                                     &
                  / (con_rd* tlyr(i,k) * (1.0_r8 + con_fvirt*qlyr(i,k)))
             cpath = ABS(plvl(i,k+1) - plvl(i,k)) * (100000.0_r8 / con_g)!g/m2
             tc    = tlyr(i,k) - con_t0c

             !! cloud water
             !
             !  ---  effective radius (recwat) & total water path (cwatp):
             !       assume monodisperse distribution of droplets (no factor of 1.5)

             IF (qcwat(i,k) > 0.0_r8) THEN
                recwat(i,k) = MAX(RECWmin,recw1*(rho*(qcwat(i,k)))**CEXP)
                cwatp (i,k) = cpath * qcwat(i,k)           ! cloud water path
                !WRITE(*,*)recwat(i,k),qcwat(i,k),recw1,rho,recw1*(rho*(qcwat(i,k)/MAX(totcnd,EPSQ)))**CEXP
                !             tem         = 5.0_r8*(1.0_r8 + max(0.0_r8, min(1.0_r8,-0.05_r8*tc)))
                !             recwat(i,k) = max(recwat(i,k), tem)
             ENDIF

             !! rain
             !
             !  ---  effective radius (rerain) & total water path (rainp):
             !       factor of 1.5_r8 accounts for r**3/r**2 moments for exponentially
             !       distributed drops in effective radius calculations
             !       (from m.d. chou's code provided to y.-t. hou)

             IF (qrain(i,k) > 0.0_r8) THEN
                tem         = CN0r0 * SQRT(SQRT(rho*(qrain(i,k))))
                rerain(i,k) = 1.5_r8 * MAX(XMRmin, MIN(XMRmax, tem))
                rainp (i,k) = cpath * qrain(i,k)           ! rain water path
             ENDIF

             !! snow (large ice) & cloud ice
             !
             !  ---  effective radius (resnow) & total ice path (snowp) for snow, and
             !       total ice path (cicep) for cloud ice:
             !       factor of 1.5 accounts for r**3/r**2 moments for exponentially
             !       distributed ice particles in effective radius calculations
             !       separation of cloud ice & "snow" uses algorithm from subroutine gsmcolumn

             !           pfac = max(0.5, sqrt(sqrt(min(1.0, pp1(i,k)*0.00004))))
             !go         pfac = max(0.5, (sqrt(min(1.0, pp1(i,k)*0.000025))))
             pfac = 1.0_r8

             IF (qcice(i,k) > 0.0_r8) THEN

                !  ---  mean particle size following houze et al. (jas, 1979, p. 160),
                !       converted from fig. 5 plot of lamdas.  an analogous set of
                !       relationships also shown by fig. 8 of ryan (bams, 1996, p. 66),
                !       but with a variety of different relationships that parallel
                !       the houze curves.

                !             dum = max(0.05, min(1.0, exp(0.0536*tc) ))
                dum = MAX(0.05_r8, MIN(1.0_r8, EXP(0.0564_r8*tc) ))
                indexs = MIN(MDImax, MAX(MDImin, INT(XMImax*dum) ))
                DUM=MAX(FLGmin(i)*pfac, DUM)

                !  ---  assumed number fraction of large ice to total (large & small) ice
                !       particles, which is based on a general impression of the literature.
                !       small ice are assumed to have a mean diameter of 50 microns.

                IF (tc >= 0.0_r8) THEN
                   flarge = FLG1P0
                ELSE
                   flarge = dum
                   !               flarge = max(FLGmin*pfac, dum)
                ENDIF
                !------------------------commented by moorthi -----------------------------
                !             elseif (tc >= -25.0_r8) then
                !
                !  ---  note that absence of cloud water (qcwat) is used as a quick
                !       substitute for calculating water subsaturation as in gsmcolumn
                !
                !               if (qcwat(i,k) <= 0.0_r8 .or. tc < -8.0_r8                 &
                !    &                                .or. tc > -3.0_r8) then
                !                 flarge = FLG0P2
                !               else
                !
                !  ---  parameterize effects of rime splintering by increasing
                !       number of small ice particles
                !
                !                 flarge = FLG0P1
                !               endif
                !             elseif (tc <= -50.0_r8) then
                !               flarge = 0.01_r8
                !             else
                !               flarge = 0.2_r8 * exp(0.1198_r8*(tc+25.0_r8))
                !             endif
                !____________________________________________________________________________

                xsimass = MASSI(MDImin) * (1.0_r8 - flarge) / flarge
                !             nlimax = 20.0e3_r8                                      !- ver3
                !             NLImax=50.E3_r8                 !- Ver3 => comment this line out
                NLImax=10.E3_r8/SQRT(DUM)       !- Ver3
                !             NLImax=5.E3_r8/sqrt(DUM)        !- Ver3
                !             NLImax=6.E3_r8/sqrt(DUM)        !- Ver3
                !             NLImax=7.5E3_r8/sqrt(DUM)       !- Ver3

                !             indexs = min(MDImax, max(MDImin, int(XMImax*dum) ))
                !moorthi      dsnow  = XMImax * exp(0.0536_r8*tc)
                !moorthi      indexs = max(INDEXSmin, min(MDImax, int(dsnow)))

                !             if (lprnt) print *,' rrime=',rrime,' xsimass=',xsimass,   &
                !    &       ' indexs=',indexs,' massi=',massi(indexs),' flarge=',flarge

                tem = rho * qcice(i,k)
                nlice = tem / (xsimass +rrime(i,k)*MASSI(indexs))

                !  ---  from subroutine gsmcolumn:
                !       minimum number concentration for large ice of NLImin=10/m**3
                !       at t>=0c.  done in order to prevent unrealistically small
                !       melting rates and tiny amounts of snow from falling to
                !       unrealistically warm temperatures.

                IF (tc >= 0.0_r8) THEN

                   nlice = MAX(NLImin, nlice)

                ELSEIF (nlice > nlimax) THEN

                   !  ---  ferrier 6/13/01:  prevent excess accumulation of ice

                   xli = (tem/nlimax - xsimass) / rrime(i,k)

                   IF (xli <= MASSI(450) ) THEN
                      dsnow = 9.5885e5_r8 * xli**0.42066_r8
                   ELSE
                      dsnow = 3.9751e6_r8 * xli** 0.49870_r8
                   ENDIF

                   indexs = MIN(MDImax, MAX(indexs, INT(dsnow)))
                   nlice = tem / (xsimass + rrime(i,k)*MASSI(indexs))

                ENDIF                               ! end if_tc block

                !             if (abs(plvl(i,ksfc)-plvl(i,k+k1)) < 300.0                &
                !             if (abs(plvl(i,ksfc)-plvl(i,k+k1)) < 400.0                &
                !             if (plvl(i,k+k1) > 600.0                                  &
                !    &                            .and. indexs >= INDEXSmin) then
                !             if (tc > -20.0 .and. indexs >= indexsmin) then
                IF (plvl(i,ksfc) > 850.0_r8 .AND.                            &
                                !    &            plvl(i,k+k1) > 600.0 .and. indexs >= indexsmin) then
                     &            plvl(i,k+k1) > 700.0_r8 .AND. indexs >= indexsmin) THEN ! 20060516
                   !!            if (plvl(i,ksfc) > 800.0 .and.                            &
                   !!   &            plvl(i,k+k1) > 700.0 .and. indexs >= indexsmin) then
                   !             if (plvl(i,ksfc) > 700.0 .and.                            &
                   !    &            plvl(i,k+k1) > 600.0 .and. indexs >= indexsmin) then
                   qsnow = MIN( qcice(i,k),                                &
                        &                       nlice*rrime(i,k)*MASSI(indexs)/rho )
                ENDIF

                qqcice      = MAX(0.0_r8, qcice(i,k)-qsnow)
                cicep (i,k) = cpath * qqcice          ! cloud ice path
                resnow(i,k) = 1.5_r8 * float(indexs)
                snden (i,k) = SDENS(indexs) / rrime(i,k)   ! 1/snow density
                snowp (i,k) = cpath*qsnow             ! snow path
                !             snowp (i,k) = cpath*qsnow*snden(i,k)  ! snow path / snow density

                !             if (lprnt .and. i == ipr) then
                !             if (i == 2) then
                !               print *,' L=',k,' snowp=',snowp(i,k),' cpath=',cpath,   &
                !    &         ' qsnow=',qsnow,' sden=',snden(i,k),' rrime=',rrime(i,k),&
                !    &         ' indexs=',indexs,' sdens=',sdens(indexs),' resnow=',    &
                !    &           resnow(i,k),' qcice=',qqcice,' cicep=',cicep(i,k)
                !           endif

             ENDIF                                 ! end if_qcice block
          ENDIF                                   ! end if_totcnd block

       ENDDO
    ENDDO
    !
    !...................................
  END SUBROUTINE rsipath2
  !-----------------------------------


  !#######################################################################
  !--------------- Creates lookup tables for ice processes ---------------
  !#######################################################################
  !
  SUBROUTINE ice_lookup()
    !
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------
    !
    !---- Key diameter values in mm
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Key concepts:
    !       - Actual physical diameter of particles (D)
    !       - Ratio of actual particle diameters to mean diameter (x=D/MD)
    !       - Mean diameter of exponentially distributed particles, which is the
    !         same as 1./LAMDA of the distribution (MD)
    !       - All quantitative relationships relating ice particle characteristics as
    !         functions of their diameter (e.g., ventilation coefficients, normalized
    !         accretion rates, ice content, and mass-weighted fall speeds) are a result
    !         of using composite relationships for ice crystals smaller than 1.5 mm
    !         diameter merged with analogous relationships for larger sized aggregates.
    !         Relationships are derived as functions of mean ice particle sizes assuming
    !         exponential size spectra and assuming the properties of ice crystals at
    !         sizes smaller than 1.5 mm and aggregates at larger sizes.  
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Actual ice particle diameters for which integrated distributions are derived
    !       - DminI - minimum diameter for integration (.02 mm, 20 microns)
    !       - DmaxI - maximum diameter for integration (2 cm)
    !       - DdelI - interval for integration (1 micron)
    !
    REAL(KIND=r8), PARAMETER :: DminI=.02e-3_r8
    REAL(KIND=r8), PARAMETER :: DmaxI=20.e-3_r8
    REAL(KIND=r8), PARAMETER :: DdelI=1.e-6_r8
    REAL(KIND=r8), PARAMETER :: XImin=1.e6_r8*DminI
    REAL(KIND=r8), PARAMETER :: XImax=1.e6_r8*DmaxI
    INTEGER      , PARAMETER :: IDImin=INT(XImin)
    INTEGER      , PARAMETER :: IDImax=INT(XImax)
    !
    !---- Meaning of the following arrays:
    !        - diam - ice particle diameter (m)
    !        - mass - ice particle mass (kg)
    !        - vel  - ice particle fall speeds (m/s)
    !        - vent1 - 1st term in ice particle ventilation factor
    !        - vent2 - 2nd term in ice particle ventilation factor
    !
    REAL(KIND=r8)   :: diam(IDImin:IDImax)
    REAL(KIND=r8)   :: mass(IDImin:IDImax)
    REAL(KIND=r8)   :: vel(IDImin:IDImax)

    REAL(KIND=r8)   :: vent1(IDImin:IDImax)
    REAL(KIND=r8)   :: vent2(IDImin:IDImax)
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Found from trial & error that the m(D) & V(D) mass & velocity relationships
    !       between the ice crystals and aggregates overlapped & merged near a particle
    !       diameter sizes of 1.5 mm.  Thus, ice crystal relationships are used for
    !       sizes smaller than 1.5 mm and aggregate relationships for larger sizes.
    !
    REAL(KIND=r8), PARAMETER :: d_crystal_max=1.5_r8
    !
    !---- The quantity xmax represents the maximum value of "x" in which the
    !       integrated values are calculated.  For xmax=20., this means that
    !       integrated ventilation, accretion, mass, and precipitation rates are
    !       calculated for ice particle sizes less than 20.*mdiam, the mean particle diameter.
    !
    REAL(KIND=r8), PARAMETER :: xmax=20.0_r8
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Meaning of the following arrays:
    !        - mdiam - mean diameter (m)
    !        - VENTI1 - integrated quantity associated w/ ventilation effects
    !                   (capacitance only) for calculating vapor deposition onto ice
    !        - VENTI2 - integrated quantity associated w/ ventilation effects
    !                   (with fall speed) for calculating vapor deposition onto ice
    !        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
    !        - MASSI  - integrated quantity associated w/ ice mass 
    !        - VSNOWI - mass-weighted fall speed of snow, used to calculate precip rates
    !
    !--- Mean ice-particle diameters varying from 50 microns to 1000 microns (1 mm), 
    !      assuming an exponential size distribution.  
    !
    REAL(KIND=r8) :: mdiam
    !
    !-----------------------------------------------------------------------------------
    !------------- Constants & parameters for ventilation factors of ice ---------------
    !-----------------------------------------------------------------------------------
    !
    !---- These parameters are used for calculating the ventilation factors for ice
    !       crystals between 0.2 and 1.5 mm diameter (Hall and Pruppacher, JAS, 1976).  
    !       From trial & error calculations, it was determined that the ventilation
    !       factors of smaller ice crystals could be approximated by a simple linear
    !       increase in the ventilation coefficient from 1.0 at 50 microns (.05 mm) to 
    !       1.1 at 200 microns (0.2 mm), rather than using the more complex function of
    !       1.0 + .14*(Sc**.33*Re**.5)**2 recommended by Hall & Pruppacher.
    !
    REAL(KIND=r8), PARAMETER :: cvent1i=.86_r8
    REAL(KIND=r8), PARAMETER :: cvent2i=.28_r8
    !
    !---- These parameters are used for calculating the ventilation factors for larger
    !       aggregates, where D>=1.5 mm (see Rutledge and Hobbs, JAS, 1983; 
    !       Thorpe and Mason, 1966).
    !
    REAL(KIND=r8), PARAMETER :: cvent1a=.65_r8
    REAL(KIND=r8), PARAMETER :: cvent2a=.44_r8
    !
    REAL(KIND=r8)    :: m_agg
    REAL(KIND=r8)    :: m_bullet
    REAL(KIND=r8)    :: m_column
    REAL(KIND=r8)    :: m_ice
    REAL(KIND=r8)    :: m_plate
    !
    !---- Various constants
    !
    REAL(KIND=r8), PARAMETER :: c1=2.0_r8/3.0_r8
    REAL(KIND=r8), PARAMETER :: cexp=1.0_r8/3.0_r8
    !
    !      logical :: iprint
    !      logical, parameter :: print_diag=.false.
    !
    !-----------------------------------------------------------------------------------
    !- Constants & parameters for calculating the increase in fall speed of rimed ice --
    !-----------------------------------------------------------------------------------
    !
    !---- Constants & arrays for estimating increasing fall speeds of rimed ice.
    !     Based largely on theory and results from Bohm (JAS, 1989, 2419-2427).
    !
    !-------------------- Standard atmosphere conditions at 1000 mb --------------------
    !
    REAL(KIND=r8), PARAMETER :: t_std=288.0_r8
    REAL(KIND=r8), PARAMETER :: dens_std=1000.e2_r8/(287.04_r8*288.0_r8)
    !
    !---- These "bulk densities" are the actual ice densities in the ice portion of the 
    !     lattice.  They are based on text associated w/ (12) on p. 2425 of Bohm (JAS, 
    !     1989).  Columns, plates, & bullets are assumed to have an average bulk density 
    !     of 850 kg/m**3.  Aggregates were assumed to have a slightly larger bulk density 
    !     of 600 kg/m**3 compared with dendrites (i.e., the least dense, most "lacy" & 
    !     tenous ice crystal, which was assumed to be ~500 kg/m**3 in Bohm).  
    !
    REAL(KIND=r8), PARAMETER :: dens_crystal=850.0_r8
    REAL(KIND=r8), PARAMETER :: dens_agg=600.0_r8
    !
    !--- A value of Nrime=40 for a logarithmic ratio of 1.1 yields a maximum rime factor
    !      of 1.1**40 = 45.26 that is resolved in these tables.  This allows the largest
    !      ice particles with a mean diameter of MDImax=1000 microns to achieve bulk 
    !      densities of 900 kg/m**3 for rimed ice.  
    !
    INTEGER, PARAMETER :: Nrime=40
    REAL(KIND=r8) :: m_rime
    REAL(KIND=r8) :: rime_factor(0:Nrime)
    REAL(KIND=r8) :: rime_vel   (0:Nrime)
    REAL(KIND=r8) :: vel_rime   (IDImin:IDImax,Nrime)
    REAL(KIND=r8) :: ivel_rime  (MDImin:MDImax,Nrime)
    !
    INTEGER :: i
    INTEGER :: j
    INTEGER :: jj
    INTEGER :: k
    INTEGER :: icount
    REAL(KIND=r8)    :: c2
    REAL(KIND=r8)    :: cbulk
    REAL(KIND=r8)    :: cbulk_ice
    REAL(KIND=r8)    :: px
    REAL(KIND=r8)    :: dynvis_std
    REAL(KIND=r8)    :: crime1
    REAL(KIND=r8)    :: crime2
    REAL(KIND=r8)    :: crime3
    REAL(KIND=r8)    :: crime4
    REAL(KIND=r8)    :: crime5
    REAL(KIND=r8)    :: d
    REAL(KIND=r8)    :: c_avg
    REAL(KIND=r8)    :: c_agg
    REAL(KIND=r8)    :: c_bullet
    REAL(KIND=r8)    :: c_column
    REAL(KIND=r8)    :: c_plate
    REAL(KIND=r8)    :: cl_agg
    REAL(KIND=r8)    :: cl_bullet
    REAL(KIND=r8)    :: cl_column
    REAL(KIND=r8)    :: cl_plate
    REAL(KIND=r8)    :: v_agg
    REAL(KIND=r8)    :: v_bullet
    REAL(KIND=r8)    :: v_column
    REAL(KIND=r8)    :: v_plate
    REAL(KIND=r8)    :: wd
    REAL(KIND=r8)    :: ecc_column
    REAL(KIND=r8)    :: cvent1
    REAL(KIND=r8)    :: cvent2
    REAL(KIND=r8)    :: crime_best
    REAL(KIND=r8)    :: rime_m1
    REAL(KIND=r8)    :: rime_m2
    REAL(KIND=r8)    :: x_rime
    REAL(KIND=r8)    :: re_rime
    REAL(KIND=r8)    :: smom3
    REAL(KIND=r8)    :: pratei
    REAL(KIND=r8)    :: expf
    REAL(KIND=r8)    :: bulk_dens
    REAL(KIND=r8)    :: xmass
    REAL(KIND=r8)    :: ecc_plate
    REAL(KIND=r8)    :: dx
    !
    !-----------------------------------------------------------------------------------
    !----------------------------- BEGIN EXECUTION -------------------------------------
    !-----------------------------------------------------------------------------------
    !
    !
    c2=1.0_r8/SQRT(3.0_r8)
    !     pi=acos(-1.0_r8)
    cbulk=6.0_r8/con_pi
    cbulk_ice=900.0_r8*con_pi/6.0_r8    ! Maximum bulk ice density allowed of 900 kg/m**3
    px=.4_r8**cexp             ! Convert fall speeds from 400 mb (Starr & Cox) to 1000 mb
    !
    !--------------------- Dynamic viscosity (1000 mb, 288 K) --------------------------
    !
    dynvis_std=1.496e-6_r8*t_std**1.5_r8/(t_std+120.0_r8)
    crime1=con_pi/24.0_r8
    crime2=8.0_r8*9.81_r8*dens_std/(con_pi*dynvis_std**2)
    crime3=crime1*dens_crystal
    crime4=crime1*dens_agg
    crime5=dynvis_std/dens_std

    DO i=0,Nrime
       rime_factor(i)=1.1_r8**i
    ENDDO
    !
    !#######################################################################
    !      Characteristics as functions of actual ice particle diameter 
    !#######################################################################
    !
    !----   M(D) & V(D) for 3 categories of ice crystals described by Starr 
    !----   & Cox (1985). 
    !
    !----   Capacitance & characteristic lengths for Reynolds Number calculations
    !----   are based on Young (1993; p. 144 & p. 150).  c-axis & a-axis 
    !----   relationships are from Heymsfield (JAS, 1972; Table 1, p. 1351).
    !
    icount=60
    !
    !      if (print_diag)                                                   & 
    !     &  write(7,"(2a)") '---- Increase in fall speeds of rimed ice',    &
    !     &    ' particles as function of ice particle diameter ----'
    DO i=IDImin,IDImax
       !        if (icount == 60 .and. print_diag) then
       !          write(6,"(/2a/3a)") 'Particle masses (mg), fall speeds ',     &
       !     &      '(m/s), and ventilation factors',                           &
       !     &      '  D(mm)  CR_mass   Mass_bull   Mass_col  Mass_plat ',      &
       !     &      '  Mass_agg   CR_vel  V_bul CR_col CR_pla Aggreg',          &
       !     &      '    Vent1      Vent2 '                               
       !          write(7,"(3a)") '        <----------------------------------',&
       !     &      '---------------  Rime Factor  --------------------------', &
       !     &      '--------------------------->'
       !          write(7,"(a,23f5.2)") '  D(mm)',(rime_factor(k), k=1,5),      &
       !     &       (rime_factor(k), k=6,40,2)
       !          icount=0
       !        endif
       d=(float(i)+0.5_r8)*1.e-3_r8    ! in mm
       c_avg=0.0_r8
       c_agg=0.0_r8
       c_bullet=0.0_r8
       c_column=0.0_r8
       c_plate=0.0_r8
       cl_agg=0.0_r8
       cl_bullet=0.0_r8
       cl_column=0.0_r8
       cl_plate=0.0_r8
       m_agg=0.0_r8
       m_bullet=0.0_r8
       m_column=0.0_r8
       m_plate=0.0_r8
       v_agg=0.0_r8
       v_bullet=0.0_r8
       v_column=0.0_r8
       v_plate=0.0_r8
       IF (d < d_crystal_max) THEN
          !
          !---- This block of code calculates bulk characteristics based on average
          !     characteristics of bullets, plates, & column ice crystals <1.5 mm size
          !
          !---- Mass-diameter relationships from Heymsfield (1972) & used
          !       in Starr & Cox (1985), units in mg
          !---- "d" is maximum dimension size of crystal in mm, 
          !
          ! Mass of pure ice for spherical particles, used as an upper limit for the
          !   mass of small columns (<~ 80 microns) & plates (<~ 35 microns)
          !
          m_ice=0.48_r8*d**3   ! Mass of pure ice for spherical particle
          !
          m_bullet=MIN(0.044_r8*d**3, m_ice)
          m_column=MIN(0.017_r8*d**1.7_r8, m_ice)
          m_plate=MIN(0.026_r8*d**2.5_r8, m_ice)
          !
          mass(i)=m_bullet+m_column+m_plate
          !
          !---- These relationships are from Starr & Cox (1985), applicable at 400 mb
          !---- "d" is maximum dimension size of crystal in mm, dx in microns
          !
          dx=1000.0_r8*d            ! Convert from mm to microns
          IF (dx <= 200.0_r8) THEN
             v_column=8.114e-5_r8*dx**1.585_r8
             v_bullet=5.666e-5_r8*dx**1.663_r8
             v_plate=1.e-3_r8*dx
          ELSE IF (dx <= 400.0_r8) THEN
             v_column=4.995e-3_r8*dx**0.807_r8
             v_bullet=3.197e-3_r8*dx**0.902_r8
             v_plate=1.48e-3_r8*dx**0.926_r8
          ELSE IF (dx <= 600.0_r8) THEN
             v_column=2.223e-2_r8*dx**0.558_r8
             v_bullet=2.977e-2_r8*dx**0.529_r8
             v_plate=9.5e-4_r8*dx
          ELSE IF (dx <= 800.0_r8) THEN
             v_column=4.352e-2_r8*dx**0.453_r8
             v_bullet=2.144e-2_r8*dx**0.581_r8
             v_plate=3.161e-3_r8*dx**0.812_r8
          ELSE 
             v_column=3.833e-2_r8*dx**0.472_r8
             v_bullet=3.948e-2_r8*dx**0.489_r8
             v_plate=7.109e-3_r8*dx**0.691_r8
          ENDIF
          !
          !---- Reduce fall speeds from 400 mb to 1000 mb
          !
          v_column=px*v_column
          v_bullet=px*v_bullet
          v_plate=px*v_plate
          !
          !---- DIFFERENT VERSION!  CALCULATES MASS-WEIGHTED CRYSTAL FALL SPEEDS
          !
          vel(i)=(m_bullet*v_bullet+m_column*v_column+m_plate*v_plate)/ &
               mass(i)
          mass(i)=mass(i)/3.0_r8
          !
          !---- Shape factor and characteristic length of various ice habits,
          !     capacitance is equal to 4*con_pi*(Shape factor)
          !       See Young (1993, pp. 143-152 for guidance)
          !
          !---- Bullets:
          !
          !---- Shape factor for bullets (Heymsfield, 1975)
          c_bullet=0.5_r8*d
          !---- Length-width functions for bullets from Heymsfield (JAS, 1972)
          IF (d > 0.3_r8) THEN
             wd=0.25_r8*d**0.7856_r8     ! Width (mm); a-axis
          ELSE
             wd=0.185_r8*d**0.552_r8
          ENDIF
          !---- Characteristic length for bullets (see first multiplicative term on right
          !       side of eq. 7 multiplied by crystal width on p. 821 of Heymsfield, 1975)
          cl_bullet=0.5_r8*con_pi*wd*(0.25_r8*wd+d)/(d+wd)
          !
          !---- Plates:
          !
          !---- Length-width function for plates from Heymsfield (JAS, 1972)
          wd=0.0449_r8*d**0.449_r8      ! Width or thickness (mm); c-axis
          !---- Eccentricity & shape factor for thick plates following Young (1993, p. 144)
          ecc_plate=SQRT(1.0_r8-wd*wd/(d*d))         ! Eccentricity
          c_plate=d*ecc_plate/ASIN(ecc_plate)    ! Shape factor
          !---- Characteristic length for plates following Young (1993, p. 150, eq. 6.6)
          cl_plate=d+2.0_r8*wd      ! Characteristic lengths for plates
          !
          !---- Columns:
          !
          !---- Length-width function for columns from Heymsfield (JAS, 1972)
          IF (d > 0.2_r8) THEN
             wd=0.1973_r8*d**0.414_r8    ! Width (mm); a-axis
          ELSE
             wd=0.5_r8*d             ! Width (mm); a-axis
          ENDIF
          !---- Eccentricity & shape factor for columns following Young (1993, p. 144)
          ecc_column=SQRT(1.0_r8-wd*wd/(d*d))                     ! Eccentricity
          c_column=ecc_column*d/LOG((1.0_r8+ecc_column)*d/wd)     ! Shape factor
          !---- Characteristic length for columns following Young (1993, p. 150, eq. 6.7)
          cl_column=(wd+2.0_r8*d)/(c1+c2*d/wd)       ! Characteristic lengths for columns
          !
          !---- Convert shape factor & characteristic lengths from mm to m for 
          !       ventilation calculations
          !
          c_bullet=0.001_r8*c_bullet
          c_plate=0.001_r8*c_plate
          c_column=0.001_r8*c_column
          cl_bullet=0.001_r8*cl_bullet
          cl_plate=0.001_r8*cl_plate
          cl_column=0.001_r8*cl_column
          !
          !---- Make a smooth transition between a ventilation coefficient of 1.0 at 50 microns
          !       to 1.1 at 200 microns
          !
          IF (d > 0.2_r8) THEN
             cvent1=cvent1i
             cvent2=cvent2i/3.0_r8
          ELSE
             cvent1=1.0_r8+0.1_r8*MAX(0.0_r8, d-0.05_r8)/0.15_r8
             cvent2=0.0_r8
          ENDIF
          !
          !---- Ventilation factors for ice crystals:
          !
          vent1(i)=cvent1*(c_bullet+c_plate+c_column)/3.0_r8
          vent2(i)=cvent2*(c_bullet*SQRT(v_bullet*cl_bullet)            &
               +c_plate*SQRT(v_plate*cl_plate)               &
               +c_column*SQRT(v_column*cl_column) )
          crime_best=crime3     ! For calculating Best No. of rimed ice crystals
       ELSE
          !
          !---- This block of code calculates bulk characteristics based on average
          !     characteristics of unrimed aggregates >= 1.5 mm using Locatelli & 
          !     Hobbs (JGR, 1974, 2185-2197) data.
          !
          !----- This category is a composite of aggregates of unrimed radiating 
          !-----   assemblages of dendrites or dendrites; aggregates of unrimed
          !-----   radiating assemblages of plates, side planes, bullets, & columns;
          !-----   aggregates of unrimed side planes (mass in mg, velocity in m/s)
          !
          m_agg=(0.073_r8*d**1.4_r8+0.037_r8*d**1.9_r8+0.04_r8*d**1.4_r8)/3.0_r8
          v_agg=(0.8_r8*d**0.16_r8+0.69_r8*d**0.41_r8+0.82_r8*d**0.12_r8)/3.0_r8
          mass(i)=m_agg
          vel(i)=v_agg
          !
          !---- Assume spherical aggregates
          !
          !---- Shape factor is the same as for bullets, = D/2
          c_agg=0.001_r8*0.5_r8*d         ! Units of m
          !---- Characteristic length is surface area divided by perimeter
          !       (.25*con_pi*D**2)/(con_pi*D**2) = D/4
          cl_agg=0.5_r8*c_agg         ! Units of m
          !
          !---- Ventilation factors for aggregates:
          !
          vent1(i)=cvent1a*c_agg
          vent2(i)=cvent2a*c_agg*SQRT(v_agg*cl_agg)
          crime_best=crime4     ! For calculating Best No. of rimed aggregates
       ENDIF
       !
       !---- Convert from shape factor to capacitance for ventilation factors
       !
       vent1(i)=4.0_r8*con_pi*vent1(i)
       vent2(i)=4.0_r8*con_pi*vent2(i)
       diam(i)=1.e-3_r8*d             ! Convert from mm to m
       mass(i)=1.e-6_r8*mass(i)       ! Convert from mg to kg
       !
       !---- Calculate increase in fall speeds of individual rimed ice particles
       !
       DO k=0,Nrime
          !---- Mass of rimed ice particle associated with rime_factor(k)
          rime_m1=rime_factor(k)*mass(i)
          rime_m2=cbulk_ice*diam(i)**3
          m_rime=MIN(rime_m1, rime_m2)
          !---- Best Number (X) of rimed ice particle combining eqs. (8) & (12) in Bohm
          x_rime=crime2*m_rime*(crime_best/m_rime)**0.25_r8
          !---- Reynolds Number for rimed ice particle using eq. (11) in Bohm
          re_rime=8.5_r8*(SQRT(1.0_r8+0.1519_r8*SQRT(x_rime))-1.0_r8)**2
          rime_vel(k)=crime5*re_rime/diam(i)
       ENDDO
       DO k=1,Nrime
          vel_rime(i,k)=rime_vel(k)/rime_vel(0)
       ENDDO
       !        if (print_diag) then
       !
       !---- Determine if statistics should be printed out.
       !
       !          iprint=.false.
       !          if (d <= 1.) then
       !            if (mod(i,10) == 0) iprint=.true.
       !          else
       !            if (mod(i,100) == 0) iprint=.true.
       !          endif
       !          if (iprint) then
       !            write(6,"(f7.4,5e11.4,1x,5f7.4,1x,2e11.4)")                 & 
       !     &        d,1.e6*mass(i),m_bullet,m_column,m_plate,m_agg,           &
       !     &        vel(i),v_bullet,v_column,v_plate,v_agg,                   &
       !     &        vent1(i),vent2(i)
       !            write(7,"(f7.4,23f5.2)") d,(vel_rime(i,k), k=1,5),          &
       !     &        (vel_rime(i,k), k=6,40,2)
       !            icount=icount+1
       !          endif
       !        endif
    ENDDO
    !
    !#######################################################################
    !      Characteristics as functions of mean particle diameter
    !#######################################################################
    !
    VENTI1=0.0_r8
    VENTI2=0.0_r8
    ACCRI=0.0_r8
    MASSI=0.0_r8
    VSNOWI=0.0_r8
    VEL_RF=0.0_r8
    ivel_rime=0.0_r8
    icount=0
    !      if (print_diag) then
    !        icount=60
    !        write(6,"(/2a)") '------------- Statistics as functions of ',   &
    !     &    'mean particle diameter -------------'
    !        write(7,"(/2a)") '------ Increase in fall speeds of rimed ice', &
    !     &    ' particles as functions of mean particle diameter -----'
    !      endif
    DO j=MDImin,MDImax
       !        if (icount == 60 .AND. print_diag) then
       !          write(6,"(/2a)") 'D(mm)    Vent1      Vent2    ',             &
       !     &       'Accrete       Mass     Vel  Dens  '
       !          write(7,"(3a)") '      <----------------------------------',  &
       !     &      '---------------  Rime Factor  --------------------------', &
       !     &      '--------------------------->'
       !          write(7,"(a,23f5.2)") 'D(mm)',(rime_factor(k), k=1,5),        &
       !     &       (rime_factor(k), k=6,40,2)
       !          icount=0
       !        endif
       mdiam=DelDMI*float(j)       ! in m
       smom3=0.0_r8
       pratei=0.0_r8
       rime_vel=0.0_r8                 ! Note that this array is being reused!
       DO i=IDImin,IDImax
          dx=diam(i)/mdiam
          IF (dx <= xmax) THEN      ! To prevent arithmetic underflows
             expf=EXP(-dx)*DdelI
             VENTI1(J)=VENTI1(J)+vent1(i)*expf
             VENTI2(J)=VENTI2(J)+vent2(i)*expf
             ACCRI(J)=ACCRI(J)+diam(i)*diam(i)*vel(i)*expf
             xmass=mass(i)*expf
             DO k=1,Nrime
                rime_vel(k)=rime_vel(k)+xmass*vel_rime(i,k)
             ENDDO
             MASSI(J)=MASSI(J)+xmass
             pratei=pratei+xmass*vel(i)
             smom3=smom3+diam(i)**3*expf
          ELSE
             EXIT
          ENDIF
       ENDDO
       !
       !--- Increased fall velocities functions of mean diameter (j),
       !      normalized by ice content, and rime factor (k) 
       !
       DO k=1,Nrime
          ivel_rime(j,k)=rime_vel(k)/MASSI(J)
       ENDDO
       !
       !--- Increased fall velocities functions of ice content at 0.1 mm
       !      intervals (j_100) and rime factor (k); accumulations here
       !
       jj=j/100
       IF (jj >= 2 .AND. jj <= 9) THEN
          DO k=1,Nrime
             VEL_RF(jj,k)=VEL_RF(jj,k)+ivel_rime(j,k)
          ENDDO
       ENDIF
       bulk_dens=cbulk*MASSI(J)/smom3
       VENTI1(J)=VENTI1(J)/mdiam
       VENTI2(J)=VENTI2(J)/mdiam
       ACCRI(J)=ACCRI(J)/mdiam
       VSNOWI(J)=pratei/MASSI(J)
       MASSI(J)=MASSI(J)/mdiam
       !        if (mod(j,10) == 0 .AND. print_diag) then
       !          xmdiam=1.e3*mdiam
       !          write(6,"(f5.3,4e11.4,f6.3,f8.3)") xmdiam,VENTI1(j),VENTI2(j),&
       !     &      ACCRI(j),MASSI(j),VSNOWI(j),bulk_dens
       !          write(7,"(f5.3,23f5.2)") xmdiam,(ivel_rime(j,k), k=1,5),      &
       !     &       (ivel_rime(j,k), k=6,40,2)
       !          icount=icount+1
       !        endif
    ENDDO
    !
    !--- Average increase in fall velocities rimed ice as functions of mean
    !      particle diameter (j, only need 0.1 mm intervals) and rime factor (k)
    !
    !      if (print_diag) then
    !        write(7,"(/2a)") ' ------- Increase in fall speeds of rimed ',  &
    !     &    'ice particles at reduced, 0.1-mm intervals  --------'
    !        write(7,"(3a)") '        <----------------------------------',  &
    !     &    '---------------  Rime Factor  --------------------------',   &
    !     &    '--------------------------->'
    !        write(7,"(a,23f5.2)") 'D(mm)',(rime_factor(k), k=1,5),          &
    !     &    (rime_factor(k), k=6,40,2)
    !      endif
    DO j=2,9
       VEL_RF(j,0)=1.0_r8
       DO k=1,Nrime
          VEL_RF(j,k)=0.01_r8*VEL_RF(j,k)
       ENDDO
       !        if (print_diag) write(7,"(f3.1,2x,23f5.2)") 0.1*j,              &
       !     &    (VEL_RF(j,k), k=1,5),(VEL_RF(j,k), k=6,40,2)
    ENDDO
    !
    !-----------------------------------------------------------------------------------
    !
  END SUBROUTINE ice_lookup
  !
  !#######################################################################
  !-------------- Creates lookup tables for rain processes ---------------
  !#######################################################################
  !
  SUBROUTINE rain_lookup()
    IMPLICIT NONE
    !
    !--- Parameters & arrays for fall speeds of rain as a function of rain drop
    !      diameter.  These quantities are integrated over exponential size
    !      distributions of rain drops at 1 micron intervals (DdelR) from minimum 
    !      drop sizes of .05 mm (50 microns, DminR) to maximum drop sizes of 10 mm 
    !      (DmaxR). 
    !
    REAL(KIND=r8)   , PARAMETER :: DminR=.05e-3_r8
    REAL(KIND=r8)   , PARAMETER :: DmaxR=10.e-3_r8
    REAL(KIND=r8)   , PARAMETER :: DdelR=1.e-6_r8
    REAL(KIND=r8)   , PARAMETER :: XRmin=1.e6_r8*DminR
    REAL(KIND=r8)   , PARAMETER :: XRmax=1.e6_r8*DmaxR
    INTEGER, PARAMETER :: IDRmin=INT(XRmin)
    INTEGER, PARAMETER :: IDRmax=INT(XRmax)
    REAL(KIND=r8)               :: diam(IDRmin:IDRmax)
    REAL(KIND=r8)               :: vel(IDRmin:IDRmax)
    !
    !--- Parameters rain lookup tables, which establish the range of mean drop
    !      diameters; from a minimum mean diameter of 0.05 mm (DMRmin) to a 
    !      maximum mean diameter of 0.45 mm (DMRmax).  The tables store solutions
    !      at 1 micron intervals (DelDMR) of mean drop diameter.  
    !
    REAL(KIND=r8)  :: mdiam
    !
    LOGICAL, PARAMETER :: print_diag=.FALSE.
    !
    REAL(KIND=r8)    :: d, cmass, pi2, expf
    INTEGER :: i, j, i1, i2
    !
    !-----------------------------------------------------------------------
    !------- Fall speeds of rain as function of rain drop diameter ---------
    !-----------------------------------------------------------------------
    !
    DO i=IDRmin,IDRmax
       diam(i)=float(i)*DdelR
       d=100.0_r8*diam(i)         ! Diameter in cm
       IF (d <= .42_r8) THEN
          !
          !--- Rutledge & Hobbs (1983); vel (m/s), d (cm)
          !
          !          vel(i)=max(0.0_r8, -0.267_r8+51.5_r8*d-102.25_r8*d*d+75.7_r8*d**3)
          vel(i)=MAX(0.0_r8, -0.267_r8+51.5_r8*d-102.25_r8*d*d+75.7_r8*d**3)

       ELSE IF (d > 0.42_r8 .AND. d <= .58_r8) THEN
          !
          !--- Linear interpolation of Gunn & Kinzer (1949) data
          !
          vel(i)=8.92_r8+.25_r8/(.58_r8-.42_r8)*(d-.42_r8)
       ELSE
          vel(i)=9.17_r8
       ENDIF
    ENDDO
    DO i=1,100
       i1=(i-1)*100+IDRmin
       i2=i1+90
       !
       !--- Print out rain fall speeds only for D<=5.8 mm (.58 cm)
       !
       IF (diam(i1) > .58e-2_r8) EXIT
       !        if (print_diag) then
       !          write(6,"(1x)")
       !          write(6,"('D(mm)->  ',10f7.3)") (1000.*diam(j), j=i1,i2,10)
       !          write(6,"('V(m/s)-> ',10f7.3)") (vel(j), j=i1,i2,10)
       !        endif
    ENDDO
    !
    !-----------------------------------------------------------------------
    !------------------- Lookup tables for rain processes ------------------
    !-----------------------------------------------------------------------
    !
    !     pi=acos(-1.)
    pi2=2.0_r8*con_pi
    cmass=1000.0_r8*con_pi/6.0_r8
    !      if (print_diag) then
    !        write(6,"(/'Diam - Mean diameter (mm)'                          &
    !     &          /'VENTR1 - 1st ventilation coefficient (m**2)'          &
    !     &          /'VENTR2 - 2nd ventilation coefficient (m**3/s**.5)'    &
    !     &          /'ACCRR - accretion moment (m**4/s)'                    &
    !     &          /'RHO*QR - mass content (kg/m**3) for N0r=8e6'          &
    !     &          /'RRATE - rain rate moment (m**5/s)'                    &
    !     &          /'VR - mass-weighted rain fall speed (m/s)'             &
    !     &    /' Diam      VENTR1      VENTR2       ACCRR      ',           &
    !     &    'RHO*QR       RRATE    VR')")
    !      endif
    DO j=MDRmin,MDRmax
       mdiam=float(j)*DelDMR
       VENTR2(J)=0.0_r8
       ACCRR(J)=0.0_r8
       MASSR(J)=0.0_r8
       RRATE(J)=0.0_r8
       DO i=IDRmin,IDRmax
          expf=EXP(-diam(i)/mdiam)*DdelR
          VENTR2(J)=VENTR2(J)+diam(i)**1.5_r8*vel(i)**0.5_r8*expf
          ACCRR(J)=ACCRR(J)+diam(i)*diam(i)*vel(i)*expf
          MASSR(J)=MASSR(J)+diam(i)**3*expf
          RRATE(J)=RRATE(J)+diam(i)**3*vel(i)*expf
       ENDDO
       !
       !--- Derived based on ventilation, F(D)=0.78+.31*Schmidt**(1/3)*Reynold**.5,
       !      where Reynold=(V*D*rho/dyn_vis), V is velocity, D is particle diameter,
       !      rho is air density, & dyn_vis is dynamic viscosity.  Only terms 
       !      containing velocity & diameter are retained in these tables.  
       !
       VENTR1(J)=.78_r8*pi2*mdiam**2
       VENTR2(J)=.31_r8*pi2*VENTR2(J)
       !
       MASSR(J)=cmass*MASSR(J)
       RRATE(J)=cmass*RRATE(J)
       VRAIN(J)=RRATE(J)/MASSR(J)
       !        if (print_diag) write(6,"(f5.3,5g12.5,f6.3)") 1000.*mdiam,      &
       !     &    ventr1(j),ventr2(j),accrr(j),8.e6*massr(j),rrate(j),vrain(j)
    ENDDO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE rain_lookup
  !
  !#######################################################################
  !--- Sets up lookup table for calculating initial ice crystal growth ---
  !#######################################################################
  !
  SUBROUTINE MY_GROWTH_RATES (DTPH)
    !
    IMPLICIT NONE
    !
    !--- Below are tabulated values for the predicted mass of ice crystals
    !    after 600 s of growth in water saturated conditions, based on 
    !    calculations from Miller and Young (JAS, 1979).  These values are
    !    crudely estimated from tabulated curves at 600 s from Fig. 6.9 of
    !    Young (1993).  Values at temperatures colder than -27C were 
    !    assumed to be invariant with temperature.  
    !
    !--- Used to normalize Miller & Young (1979) calculations of ice growth
    !    over large time steps using their tabulated values at 600 s.
    !    Assumes 3D growth with time**1.5 following eq. (6.3) in Young (1993).
    !
    REAL(KIND=r8) :: dtph, dt_ice
    REAL(KIND=r8) :: MY_600(MY_T1:MY_T2)
    !
    !-- 20090714: These values are in g and need to be converted to kg below
    DATA MY_600 /                                                     &
         & 5.5e-8_r8,  1.4E-7_r8,  2.8E-7_r8, 6.E-7_r8,   3.3E-6_r8,                       & !  -1 to  -5 deg C
         & 2.E-6_r8,   9.E-7_r8,   8.8E-7_r8, 8.2E-7_r8,  9.4e-7_r8,                       & !  -6 to -10 deg C
         & 1.2E-6_r8,  1.85E-6_r8, 5.5E-6_r8, 1.5E-5_r8,  1.7E-5_r8,                       & ! -11 to -15 deg C
         & 1.5E-5_r8,  1.E-5_r8,   3.4E-6_r8, 1.85E-6_r8, 1.35E-6_r8,                      & ! -16 to -20 deg C
         & 1.05E-6_r8, 1.E-6_r8,   9.5E-7_r8, 9.0E-7_r8, 9.5E-7_r8,                       &  ! -21 to -25 deg C
         & 9.5E-7_r8,  9.E-7_r8,   9.E-7_r8,  9.E-7_r8,   9.E-7_r8,                        &  ! -26 to -30 deg C
         & 9.E-7_r8,   9.E-7_r8,   9.E-7_r8,  9.E-7_r8,   9.E-7_r8/                         ! -31 to -35 deg C
    !
    !-----------------------------------------------------------------------
    !
    DT_ICE=(DTPH/600.0_r8)**1.5_r8
    !     MY_GROWTH=DT_ICE*MY_600          ! original version
    MY_GROWTH=DT_ICE*MY_600*1.E-3_r8    !-- 20090714: Convert from g to kg
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE MY_GROWTH_RATES

  !-------------------------------------------------------------------------------
  SUBROUTINE gpvsl
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvsl        Compute saturation vapor pressure table over liquid
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvsl.
    !   Exact saturation vapor pressures are calculated in subprogram fpvslx.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gpvsl
    !
    ! Subprograms called:
    !   (fpvslx)   inlinable function to compute saturation vapor pressure
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER  :: jx
    REAL(KIND=r8) :: xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvsl-1)
    c2xpvsl=1.0_r8/xinc
    c1xpvsl=1.0_r8-xmin*c2xpvsl
    DO jx=1,nxpvsl
       x=xmin+(jx-1)*xinc
       t=x
       tbpvsl(jx)=fpvslx(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvsl
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvslx(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvslx       Compute saturation vapor pressure over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas and liquid, and neglects the volume of the liquid.
    !   The model does account for the variation of the latent heat
    !   of condensation with temperature.  The ice option is not included.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   pvsl=fpvslx(t)
    !
    !   Input argument list:
    !     t          REAL(KIND=r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvslx     REAL(KIND=r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8) :: fpvslx
    REAL(KIND=r8),INTENT(in):: t
    REAL(KIND=r8),PARAMETER:: dldt=con_cvap-con_cliq
    REAL(KIND=r8),PARAMETER:: heat=con_hvap
    REAL(KIND=r8),PARAMETER:: xpona=-dldt/con_rv
    REAL(KIND=r8),PARAMETER:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(KIND=r8) tr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    fpvslx=con_psat*(tr**xpona)*EXP(xponb*(1.0_r8-tr))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvslx

  !
  !   fpvsl           Elementally compute saturation vapor pressure over liquid
  !     function result REAL(KIND=r8) saturation vapor pressure in Pascals
  !     t               REAL(KIND=r8) temperature in Kelvin
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsl(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsl        Compute saturation vapor pressure over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvsl. See documentation for fpvslx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvsl is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   pvsl=fpvsl(t)
    !
    !   Input argument list:
    !     t          REAL(KIND=r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsl      REAL(KIND=r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8)           :: fpvsl
    REAL(KIND=r8),INTENT(in):: t
    INTEGER  :: jx
    REAL(KIND=r8) :: xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsl+c2xpvsl*t,1.0_r8),REAL(nxpvsl,r8))
    jx=INT(MIN(xj,nxpvsl-1.0_r8))
    fpvsl=tbpvsl(jx)+(xj-jx)*(tbpvsl(jx+1)-tbpvsl(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsl

  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsi(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsi        Compute saturation vapor pressure over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvsi. See documentation for fpvsix for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvsi is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvsi=fpvsi(t)
    !
    !   Input argument list:
    !     t          REAL(KIND=r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsi      REAL(KIND=r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8) ::  fpvsi
    REAL(KIND=r8),INTENT(in):: t
    INTEGER ::jx
    REAL(KIND=r8):: xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsi+c2xpvsi*t,1.0_r8),REAL(nxpvsi,r8))
    jx=INT(MIN(xj,nxpvsi-1.0_r8))
    fpvsi=tbpvsi(jx)+(xj-jx)*(tbpvsi(jx+1)-tbpvsi(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsi


  !-------------------------------------------------------------------------------
  SUBROUTINE gpvsi()
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvsi        Compute saturation vapor pressure table over ice
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvsi.
    !   Exact saturation vapor pressures are calculated in subprogram fpvsix.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gpvsi
    !
    ! Subprograms called:
    !   (fpvsix)   inlinable function to compute saturation vapor pressure
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER  :: jx
    REAL(KIND=r8) :: xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvsi-1)
    c2xpvsi=1.0_r8/xinc
    c1xpvsi=1.0_r8-xmin*c2xpvsi
    DO jx=1,nxpvsi
       x=xmin+(jx-1)*xinc
       t=x
       tbpvsi(jx)=fpvsix(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvsi


  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsix(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsix       Compute saturation vapor pressure over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas and ice, and neglects the volume of the ice.
    !   The model does account for the variation of the latent heat
    !   of condensation with temperature.  The liquid option is not included.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsi=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvsi=fpvsix(t)
    !
    !   Input argument list:
    !     t          REAL(KIND=r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsix     REAL(KIND=r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8) fpvsix
    REAL(KIND=r8),INTENT(in):: t
    REAL(KIND=r8),PARAMETER:: dldt=con_cvap-con_csol
    REAL(KIND=r8),PARAMETER:: heat=con_hvap+con_hfus
    REAL(KIND=r8),PARAMETER:: xpona=-dldt/con_rv
    REAL(KIND=r8),PARAMETER:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(KIND=r8) tr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    fpvsix=con_psat*(tr**xpona)*EXP(xponb*(1.0_r8-tr))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsix

  !-------------------------------------------------------------------------------
  SUBROUTINE gpvs()
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvs         Compute saturation vapor pressure table
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvs.
    !   Exact saturation vapor pressures are calculated in subprogram fpvsx.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gpvs
    !
    ! Subprograms called:
    !   (fpvsx)    inlinable function to compute saturation vapor pressure
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(KIND=r8) xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvs-1)
    !   c1xpvs=1.-xmin/xinc
    c2xpvs=1./xinc
    c1xpvs=1.-xmin*c2xpvs
    DO jx=1,nxpvs
       x=xmin+(jx-1)*xinc
       t=x
       tbpvs(jx)=fpvsx(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvs
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsx(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsx        Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The saturation vapor pressure over either liquid and ice is computed
    !   over liquid for temperatures above the triple point,
    !   over ice for temperatures 20 degress below the triple point,
    !   and a linear combination of the two for temperatures in between.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas, liquid and ice, and neglects the volume of the condensate.
    !   The model does account for the variation of the latent heat
    !   of condensation and sublimation with temperature.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   The reference for this computation is Emanuel(1994), pages 116-117.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvsx(t)
    !
    !   Input argument list:
    !     t          REAL(KIND=r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsx      REAL(KIND=r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8) fpvsx
    REAL(KIND=r8),INTENT(in):: t
    REAL(KIND=r8),PARAMETER:: tliq=con_ttp
    REAL(KIND=r8),PARAMETER:: tice=con_ttp-20.0_r8
    REAL(KIND=r8),PARAMETER:: dldtl=con_cvap-con_cliq
    REAL(KIND=r8),PARAMETER:: heatl=con_hvap
    REAL(KIND=r8),PARAMETER:: xponal=-dldtl/con_rv
    REAL(KIND=r8),PARAMETER:: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
    REAL(KIND=r8),PARAMETER:: dldti=con_cvap-con_csol
    REAL(KIND=r8),PARAMETER:: heati=con_hvap+con_hfus
    REAL(KIND=r8),PARAMETER:: xponai=-dldti/con_rv
    REAL(KIND=r8),PARAMETER:: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
    REAL(KIND=r8) tr,w,pvl,pvi
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    IF(t.GE.tliq) THEN
       fpvsx=con_psat*(tr**xponal)*EXP(xponbl*(1._r8-tr))
    ELSEIF(t.LT.tice) THEN
       fpvsx=con_psat*(tr**xponai)*EXP(xponbi*(1._r8-tr))
    ELSE
       w=(t-tice)/(tliq-tice)
       pvl=con_psat*(tr**xponal)*EXP(xponbl*(1._r8-tr))
       pvi=con_psat*(tr**xponai)*EXP(xponbi*(1._r8-tr))
       fpvsx=w*pvl+(1._r8-w)*pvi
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsx
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvs(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvs         Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvs. See documentation for fpvsx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvs is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvs(t)
    !
    !   Input argument list:
    !     t          REAL(KIND=r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvs       REAL(KIND=r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(KIND=r8)           :: fpvs
    REAL(KIND=r8),INTENT(in):: t
    INTEGER :: jx
    REAL(KIND=r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvs+c2xpvs*t,1.0_r8),REAL(nxpvs,r8))
    jx=INT(MIN(xj,nxpvs-1._r8))
    fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvs

  
  SUBROUTINE Cloud_Micro_WRF(&
       ! Model info
       ncols, kmax ,  delsig, imask   , ILCON,schemes,&
       ! Atmospheric Fields
       Ps   , Te   , Qe    , tsea  , FlipPbot ,        &
       QCF  ,QCL   ,  cld ,clu ,                                 &
       ! Cloud properties
       clwp , lmixr, fice  , rei   , rel   , taud ,cicewp ,      &
       cliqwp       ,      &
       c_cld_tau    ,      &
       c_cld_tau_w  ,      &
       c_cld_tau_w_g,      &
       c_cld_tau_w_f,      &
       c_cld_lw_abs ,      &
       cldfprime          )
    IMPLICIT NONE

    ! As in the CCM2, cloud optical properties in the CCM3 are accounted for using
    ! the Slingo (1989) parameterization for liquid water droplet clouds. This
    ! scheme relates the extinction optical depth, the single-scattering albedo,
    ! and the asymmetry parameter to the cloud liquid water path and cloud drop
    ! effective radius. The latter two microphysical cloud properties were
    ! statically specified in the CCM2. In particular, in-cloud liquid water paths
    ! were evaluated from a prescribed, meridionally and height varying, but
    ! time independent, cloud liquid water density profile, rho_l(z), which
    ! was analytically determined on the basis of a meridionally specified
    ! liquid water scale height (e.g. see Kiehl et al., 1994; Kiehl, 1991).
    ! The cloud drop effective radius was simplly specified to be 10microns
    ! for all clouds. The CCM3 continues to diagnose cloud optical properties,
    ! but relaxes the rigid CCM2 framework. CCM3 employs the same exponentially
    ! decaying vertical profile for in-cloud water concentration
    !
    !             rho_l=rho_l^0*exp(-z/h_l)               eq 4.a.11
    !
    ! , where rho_l^0=0.21g/m3. Instead of specifying a zonally symmetric meridional
    ! dependence for the cloud water scale heigh, h_l, it is locally diagnosed
    ! as a function of the vertically integrated water vapor (precipitable water) 
    !
    !          h_l=700 ln [1+\frac{1}{g} \int_pT^ps q dp]  eq 4.a.12
    !
    ! hmjb> It is not explained, but the units of h_l must be meters, the same 
    ! hmjb> of the height, z.
    !
    ! The cloud water path (CWP) is determined by integrating the liquid
    ! water concentration using
    !
    !                 cwp = int rho_l dz     eq. 4.a.13
    ! 
    ! Which can be analytically evaluated for an arbitrary layer k as
    !
    !  rho_l^0 h_l [exp(-z_bot(k)/h_l) - exp(-z_top(k)/h_l)]   eq. 4.a.14
    !
    ! Where z_bot and z_top are the heights of the k'th layer interfaces.
    !
    ! hmjb> It is not explained, but the units of clwp must be g/m2
    ! hmjb> since it is the integral of rho_l*dz (eq.4.a.13)
    !
    ! CCM3 Documentation, pg 50
    ! Observational studies have shown a distinct difference between
    ! maritime and continental effective cloud drop size, r_e, for warm
    ! clouds. For this reason, the CCM3 differentiates between the cloud
    ! drop effective radius for clouds diagnosed over maritime and
    ! continental regimes (Kiehl, 1994). Over the ocean, the cloud drop
    ! effective radius for liquid water clouds, r_el, is specified to be
    ! 10microns, as in the CCM3. Over land masses r_el is determinedusing
    !
    ! r_el = 5 microns             T > -10oC
    !      = 5-5(t+10)/20 microns  -30oC <= T <= -10oC     eq. 4.a.14.1
    !      = r_ei                  T < -30oC
    !
    ! An ice particle effective radius, r_ei, is also diagnosed by CCM3,
    ! which at the moment amounts to a specification of ice radius as a
    ! function of normalized pressure
    !
    ! r_ei = 10 microns                                 p/ps > p_I^high   
    !      = r_ei^max - (r_ei^max - r_ei^min)           p/ps <= p_I^high       eq. 4.a.15.1
    !            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    ! where r_ei^max=30microns, r_ei^min=10microns, p_I^high=0.4 and p_I^low=0.0
    !
    ! hmjb>> I think there is a typo in the equation, otherwise the 
    ! hmjb>> expression for r_ei is not a continuous funcion of p/ps.
    ! hmjb>> For p/ps=p_I^high, r_ei should be r_ei^min and not r_ei^max.
    ! hmjb>> The correct equation is:
    ! hmjb>> r_ei = 10 microns                                 p/ps > p_I^high   
    ! hmjb>>      = R_EI^MIN - (r_ei^max - r_ei^min)           p/ps <= p_I^high
    ! hmjb>>            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    !--------------------------------------------------------------------------------- 
    ! Input/Output Variables
    !--------------------------------------------------------------------------------- 

    ! Model info
    INTEGER         , INTENT(IN   ) :: ncols  
    INTEGER         , INTENT(IN   ) :: kmax   
    REAL(KIND=r8)   , INTENT(IN   ) :: delsig   (kmax)  ! Layer thickness (sigma)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask    (ncols) ! Ocean/Land mask
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    INTEGER         , INTENT(IN   ) :: schemes
    ! Atmospheric Fields
    REAL(KIND=r8)   , INTENT(IN   ) :: Ps       (ncols)      ! Surface pressure (mb)
    REAL(KIND=r8)   , INTENT(IN   ) :: Te       (ncols,kmax) ! Temperature (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: Qe       (ncols,kmax) ! Specific Humidity (g/g)
    REAL(KIND=r8)   , INTENT(IN   ) :: tsea     (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: FlipPbot (ncols,kmax)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCF      (ncols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCL      (ncols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cld      (ncols,kMax) 
    REAL(KIND=r8)   , INTENT(IN   ) :: clu      (ncols,kMax) 

    ! Cloud properties
    REAL(KIND=r8)  , INTENT(OUT   ) :: clwp         (ncols,kmax) ! Cloud Liquid Water Path
    REAL(KIND=r8)  , INTENT(OUT   ) :: lmixr        (ncols,kmax) ! Ice/Water mixing ratio
    REAL(KIND=r8)  , INTENT(OUT   ) :: fice         (ncols,kmax) ! Fractional amount of cloud that is ice
    REAL(KIND=r8)  , INTENT(OUT   ) :: rei          (ncols,kmax) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8)  , INTENT(OUT   ) :: rel          (ncols,kmax) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8)  , INTENT(OUT   ) :: taud         (ncols,kmax) ! Shortwave cloud optical depth
    REAL(KIND=r8)  , INTENT(OUT   ) :: cicewp       (ncols,kmax)
    REAL(KIND=r8)  , INTENT(OUT   ) :: cliqwp       (ncols,kmax)
    REAL(r8)       , INTENT(OUT   ) :: c_cld_tau    (1:nbndsw,1:ncols,1:kmax) ! cloud extinction optical depth
    REAL(r8)       , INTENT(OUT   ) :: c_cld_tau_w  (1:nbndsw,1:ncols,1:kmax) ! cloud single scattering albedo * tau
    REAL(r8)       , INTENT(OUT   ) :: c_cld_tau_w_g(1:nbndsw,1:ncols,1:kmax) ! cloud assymetry parameter * w * tau
    REAL(r8)       , INTENT(OUT   ) :: c_cld_tau_w_f(1:nbndsw,1:ncols,1:kmax) ! cloud forward scattered fraction * w * tau
    REAL(r8)       , INTENT(OUT   ) :: c_cld_lw_abs (1:nbndlw,1:ncols,1:kmax) ! cloud absorption optics depth (LW)
    REAL(r8)       , INTENT(OUT   ) :: cldfprime    (1:ncols,1:kmax)          ! combined cloud fraction (snow plus regular)


    !--------------------------------------------------------------------------------- 
    ! Parameters
    !--------------------------------------------------------------------------------- 

    REAL(KIND=r8), PARAMETER :: abarl=2.261e-2_r8
    REAL(KIND=r8), PARAMETER :: bbarl=1.4365_r8
    REAL(KIND=r8), PARAMETER :: abari=3.448e-3_r8
    REAL(KIND=r8), PARAMETER :: bbari=2.431_r8
    REAL(KIND=r8), PARAMETER :: kabsl=0.090361_r8               ! longwave liquid absorption coeff (m**2/g)

    REAL(KIND=r8), PARAMETER :: clwc0   = 0.21_r8 ! Reference liquid water concentration (g/m3)        
    REAL(KIND=r8), PARAMETER :: reimin  = 10.0_r8 ! Minimum of Ice particle efective radius (microns)  
    REAL(KIND=r8), PARAMETER :: reirnge = 20.0_r8 ! Range of Ice particle efective radius (microns)   
    REAL(KIND=r8), PARAMETER :: sigrnge = 0.4_r8  ! Normalized pressure range                         
    REAL(KIND=r8), PARAMETER :: sigmax  = 0.4_r8  ! Normalized pressure maximum                       

    !REAL(KIND=r8), PARAMETER :: pptop = 0.005_r8       ! Model-top presure                                 

    !--------------------------------------------------------------------------------- 
    ! Local Variables
    !--------------------------------------------------------------------------------- 
    REAL(KIND=r8)      :: landm    (1:ncols)   ! Land fraction ramped
    REAL(KIND=r8)      :: icefrac  (1:ncols)   ! Ice fraction
    REAL(KIND=r8)      :: snowh    (1:ncols)   ! Snow depth over land, water equivalent (m)
    REAL(KIND=r8)      :: ocnfrac  (1:nCols)   ! Ocean fraction
    REAL(KIND=r8)      :: landfrac (1:ncols)   ! Land fraction

    REAL(KIND=r8)      :: hl       (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8)      :: rhl      (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8)      :: pw       (ncols)        ! precipitable water (kg/m2)
    REAL(KIND=r8)      :: Zibot    (ncols,kmax+1) ! Height at middle of layer (m)
    REAL(KIND=r8)      :: emziohl  (ncols,kmax+1) ! exponential of Minus zi Over hl (no dim)
    REAL(KIND=r8)      :: pdel     (ncols,kmax)   ! Moist pressure difference across layer Pressure thickness [Pa] > 0
    REAL(KIND=r8)      :: prsl     (ncols,kmax)   !      prsi  (IX,LM+1) : model level pressure in cb      (kPa) !
    REAL(KIND=r8)      :: prsi     (ncols,kmax+1) !      prsl  (IX,LM)   : model layer mean pressure in cb (kPa)          !
    REAL(KIND=r8)      :: cldfrac  (ncols,kmax)   ! 
    REAL(KIND=r8)      :: cldfsnow (ncols,kmax) 
    REAL(KIND=r8)      :: lamc     (ncols,kmax) 
    REAL(KIND=r8)      :: pgam     (ncols,kmax) 
    REAL(KIND=r8)      :: emis     (ncols,kMax)       ! cloud emissivity (fraction)
    LOGICAL            :: dosw
    LOGICAL            :: dolw
    LOGICAL            :: oldcldoptics
    CHARACTER(LEN=200) :: liqcldoptics
    CHARACTER(LEN=200) :: icecldoptics
    INTEGER            :: cldfsnow_idx

    REAL(KIND=r8)      ::gicewp
    REAL(KIND=r8)      ::gliqwp
    !-- Aux variables

    INTEGER :: i,k
!    REAL(KIND=r8) :: weight
    REAL(KIND=r8) :: kabs                   ! longwave absorption coeff (m**2/g)
    REAL(KIND=r8) :: kabsi                  ! ice absorption coefficient
    

    !--------------------------------------------------------------------------------- 
    !--------------------------------------------------------------------------------- 

    clwp=0.0_r8
    lmixr=0.0_r8
    fice=0.0_r8
    rei=0.0_r8
    rel=0.0_r8
    taud=0.0_r8
    cicewp=0.0_r8
    cliqwp=0.0_r8
    c_cld_tau=0.0_r8
    c_cld_tau_w=0.0_r8
    c_cld_tau_w_g=0.0_r8
    c_cld_tau_w_f=0.0_r8
    c_cld_lw_abs =0.0_r8
    cldfprime=0.0_r8
    landm=0.0_r8
    snowh=0.0_r8
    icefrac=0.0_r8
    landfrac=0.0_r8
    ocnfrac=0.0_r8

    hl   =0.0_r8
    rhl =0.0_r8
    pw =0.0_r8
    Zibot =0.0_r8
    emziohl  =0.0_r8
    pdel =0.0_r8
    prsl =0.0_r8
    prsi =0.0_r8
    cldfrac  =0.0_r8
    cldfsnow =0.0_r8
    lamc =0.0_r8
    pgam =0.0_r8
    emis =0.0_r8

    DO i=1,nCols
       IF(schemes == 1 .and. imask(i) == 13_i8) snowh(i)=5.0_r8
       IF(schemes == 2 .and. imask(i) == 13_i8) snowh(i)=5.0_r8
       IF(schemes == 3 .and. imask(i) == 15_i8) snowh(i)=5.0_r8
       IF(imask(i) >   0_i8)THEN
          ! land
          icefrac  (i)=0.0_r8
          landfrac (i)=1.0_r8
          ocnfrac  (i)=0.0_r8
       ELSE
          ! water/ocean
          landfrac  (i) =0.0_r8
          ocnfrac   (i) =1.0_r8
          IF(ocnfrac(i).GT.0.01_r8.AND.ABS(tsea(i)).LT.260.0_r8) THEN
             icefrac(i) = 1.0_r8
             ocnfrac(i) = 1.0_r8
          ENDIF
       END IF
    END DO

    DO i = 1, nCols
       prsi(i,kMax + 1)=MAX(ps(i)*si(kMax + 1)/10.0_r8,1.0e-12_r8) !mb  -- > cb
    END DO

    DO k = 1, kMax
       DO i = 1, nCols
          prsl(i,k)=ps(i)*si(k)/10.0_r8 ! cb
          prsi(i,k)=ps(i)*si(k)/10.0_r8 ! cb
       ENDDO
    ENDDO

    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,kMax
       DO i=1,nCols
           pdel(i,k)   = (prsi(i,k)-prsi(i,k+1))*1000.0_r8! cb -- > Pa
       END DO
    END DO

    ! Heights corresponding to sigma at middle of layer: sig(k)
    ! Assuming isothermal atmosphere within each layer

    DO i=1,nCols
       Zibot(i,1) = 0.0_r8
       DO k=2,kMax
          Zibot(i,k) = Zibot(i,k-1) + (con_rd/con_g)*Te(i,k-1)* &
               !               LOG(sigbot(k-1)/sigbot(k))
               LOG(FlipPbot(i,kMax+2-k)/FlipPbot(i,kMax+1-k))
       END DO
    END DO

    DO i=1,nCols
       Zibot(i,kMax+1)=Zibot(i,kMax)+(con_rd/con_g)*Te(i,kMax)* &
            LOG(FlipPbot(i,1)/pptop)
    END DO
    

    ! precitable water, pw = sum_k { delsig(k) . Qe(k) } . Ps . 100 / g
    !                   pw = sum_k { Dp(k) . Qe(k) } / g
    !
    ! 100 is to change from mbar to pascal
    ! Dp(k) is the difference of pressure (N/m2) between bottom and top of layer
    ! Qe(k) is specific humidity in (g/g)
    ! gravity is m/s2 => so pw is in Kg/m2
    DO k=1,kmax
       DO i = 1,ncols
          pw(i) = pw(i) + delsig(k)*Qe(i,k)
       END DO
    END DO
    DO i = 1,ncols
       pw(i)=100.0_r8*pw(i)*Ps(i)/con_g
    END DO
    !
    ! diagnose liquid water scale height from precipitable water
    DO i=1,ncols
       hl(i)  = 700.0_r8*LOG(MAX(pw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,kmax+1
       DO i=1,ncols
          emziohl(i,k) = EXP(-zibot(i,k)*rhl(i))
       END DO
    END DO
    !    DO i=1,ncols
    !       emziohl(i,kmax+1) = 0.0_r8
    !    END DO

    ! The units are g/m2.
    DO k=1,kmax
       DO i=1,ncols
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k) - emziohl(i,k+1))
       END DO
    END DO

    ! If we want to calculate the 'droplets/cristals' mixing ratio, we need
    ! to find the amount of dry air in each layer. 
    !
    !             dry_air_path = int rho_air dz  
    !
    ! This can be simply done using the hydrostatic equation:
    !
    !              dp/dz   = -rho grav
    !              dp/grav = -rho dz
    !
    !
    ! The units are g/m2. The factor 1e5 accounts for the change
    !  mbar to Pa and kg/m2 to g/m2. 
    DO k=1,kmax
       DO i=1,ncols
          lmixr(i,k)=clwp(i,k)*con_g*1.0e-5_r8/delsig(k)/Ps(i)
       END DO
    END DO

    !
    ! Cloud water and ice particle sizes, saved in physics buffer for radiation
    ! Author: Byron Boville  Sept 06, 2002, assembled from existing subroutines
    !
    CALL cldefr(&
         ncols                             , & !INTEGER , INTENT(in ) :: pcols                ! number of atmospheric columns
         kmax                              , & !INTEGER , INTENT(in ) :: kMax                 ! number of vertical levels
         Te         (1:ncols,1:kMax)       , & !REAL(r8), INTENT(in ) :: t       (pcols,kMax) ! Temperature
         rel        (1:ncols,1:kMax)       , & !REAL(r8), INTENT(out) :: rel     (pcols,kMax) ! Liquid effective drop size (microns)
         rei        (1:ncols,1:kMax)       , & !REAL(r8), INTENT(out) :: rei     (pcols,kMax) ! Ice effective drop size (microns)
         landm      (1:ncols)              , & !REAL(r8), INTENT(in ) :: landm   (pcols)      !
         icefrac    (1:ncols)              , & !REAL(r8), INTENT(in ) :: icefrac (pcols)      ! Ice fraction
         snowh      (1:ncols)                ) !REAL(r8), INTENT(in ) :: snowh   (pcols)      ! Snow depth over land, water equivalent (m)
 
 
       DO k=1,kmax
          DO i=1,ncols
             cldfrac(i,k) =MAX(cld(i,kmax-k+1),clu(i,kmax-k+1))
             cldfrac(i,k) =MAX(cldfrac(i,k) , 2.0e-80_r8)
          END DO
       END DO

! From module_ra_cam: Convert liquid and ice mixing ratios to water paths;
! pdel is in mb here; convert back to Pa (*100.)
! Water paths are in units of g/m2
! snow added as ice cloud (JD 091022)
    do k = 1,kmax
       DO i=1,ncols

          ! gicewp = (QCF(i,k)+qs1d(k)) * pdel(ncol,k)*100.0 / gravmks * 1000.0     ! Grid box ice water path.
          gicewp = (QCF(i,k)) * pdel(i,k) / con_g * 1000.0_r8     ! Grid box ice water path.
          gliqwp = (QCL(i,k)) * pdel(i,k) / con_g * 1000.0_r8       ! Grid box liquid water path.
          cicewp(i,k) = gicewp / max(0.01_r8,cldfrac(i,k))                    ! In-cloud ice water path.
          cliqwp(i,k) = gliqwp / max(0.01_r8,cldfrac(i,k))                    ! In-cloud liquid water path.
       END DO
    end do


    ! define fractional amount of cloud that is ice
    ! if warmer than -10 degrees c then water phase
    ! docs CCM3, eq 4.a.16.1     
    ! allcld_liq = state%q(:,:,ixcldliq)
    ! allcld_ice = state%q(:,:,ixcldice)
    IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
       DO k=1,kmax
          DO i=1,ncols
             fice(i,k)=MAX(MIN((263.16_r8-Te(i,k))*0.05_r8,1.0_r8),0.0_r8)
             !fice(i,k) = allcld_ice(i,k) /max(1.e-10_r8,(allcld_ice(i,k) + allcld_liq(i,k)))  
          END DO
       END DO    
    ELSE IF ( TRIM(ILCON) == 'MIC'.or. TRIM(ILCON).EQ.'HWRF' .or. TRIM(ILCON).EQ.'HGFS'.or.&
              TRIM(ILCON).EQ.'UKMO' .or. TRIM(ILCON).EQ.'MORR' .or.TRIM(ILCON).EQ.'HUMO') THEN
       DO k=1,kmax
          DO i=1,ncols
             fice(i,k) = MIN(MAX( QCF(i,k) /max(1.e-10_r8,(QCF(i,k) + QCL(i,k))),0.000_r8),1.0_r8)
             !fice(i,k) = allcld_ice(i,k) /max(1.e-10_r8,(allcld_ice(i,k) + allcld_liq(i,k)))
          END DO
       END DO
    END IF
    cldfsnow_idx=0
    cldfsnow=0.0_r8
    lamc=0.0_r8
    pgam=0.0_r8
    oldcldoptics=.FALSE.
    liqcldoptics='slingo'
    icecldoptics='ebertcurry'!  ('mitchell')
    dosw=.TRUE.  
    dolw=.TRUE. 
    CALL Run_Optical_Properties( &
       ncols                         , &
       ncols                         , &
       kmax                          , &
       cldfsnow_idx                  , &
       cldfsnow     (1:nCols,1:kMax) , &
       cldfrac      (1:nCols,1:kMax) , &
       QCL          (1:nCols,1:kMax) , &
       QCF          (1:nCols,1:kMax) , &
       rel          (1:nCols,1:kMax) , &!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
       rei          (1:nCols,1:kMax) , &
       cliqwp       (1:nCols,1:kMax) , &
       cicewp       (1:nCols,1:kMax) , &
       clwp         (1:nCols,1:kMax) , &
       pdel         (1:nCols,1:kMax) , &
       lamc         (1:nCols,1:kMax) , &
       pgam         (1:nCols,1:kMax) , &
    ! combined cloud radiative parameters are "in cloud" not "in cell"
       c_cld_tau    (1:nbndsw,1:nCols,1:kMax), & ! cloud extinction optical depth
       c_cld_tau_w  (1:nbndsw,1:nCols,1:kMax), & ! cloud single scattering albedo * tau
       c_cld_tau_w_g(1:nbndsw,1:nCols,1:kMax), & ! cloud assymetry parameter * w * tau
       c_cld_tau_w_f(1:nbndsw,1:nCols,1:kMax), & ! cloud forward scattered fraction * w * tau
       c_cld_lw_abs (1:nbndlw,1:nCols,1:kMax) ,& ! cloud absorption optics depth (LW)
       cldfprime    (1:nCols,1:kMax)  , &
       dosw                      , &
       dolw                      , &
       oldcldoptics              , &
       icecldoptics              , &
       liqcldoptics              )
    ! Compute optical depth from liquid water
    DO k=1,kMax
       DO i=1,nCols
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          !if ( microp_scheme .eq. 'MG' ) then
          kabsi = 0.005_r8 + 1.0_r8/min(max(13._r8,rei(i,k)),130._r8)
          !else if ( microp_scheme .eq. 'RK' ) then
          !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
          !END IF
          !     (m**2/g)
          kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 
          ! cloud emissivity (fraction)
          emis(i,k) = 1.0_r8 - exp(-1.66_r8*kabs*clwp(i,k))
          ! cloud optical depth
          taud(i,k) = kabs*clwp(i,k)! g/m2
       END DO
    END DO


  END SUBROUTINE Cloud_Micro_WRF

 
  !===============================================================================
  SUBROUTINE cldefr( &
       ncols    , &!INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
       kMax    , &!INTEGER, INTENT(in) :: kMax                  ! number of vertical levels
       t       , &!REAL(r8), INTENT(in) :: t       (pcols,kMax)        ! Temperature
       rel     , &!REAL(r8), INTENT(out) :: rel(pcols,kMax)      ! Liquid effective drop size (microns)
       rei     , &!REAL(r8), INTENT(out) :: rei(pcols,kMax)      ! Ice effective drop size (microns)
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
    INTEGER, INTENT(in) :: ncols                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: kMax                  ! number of vertical levels

    REAL(r8), INTENT(in) :: icefrac (ncols)       ! Ice fraction
    REAL(r8), INTENT(in) :: t       (ncols,kMax)  ! Temperature
    REAL(r8), INTENT(in) :: landm   (ncols)
    REAL(r8), INTENT(in) :: snowh   (ncols)       ! Snow depth over land, water equivalent (m)
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: rel(ncols,kMax)      ! Liquid effective drop size (microns)
    REAL(r8), INTENT(out) :: rei(ncols,kMax)      ! Ice effective drop size (microns)
    !

    !++pjr
    ! following Kiehl
    CALL reltab(ncols,kMax,  t,  landm, icefrac, rel, snowh)

    ! following Kristjansson and Mitchell
    CALL reitab(ncols,kMax, t, rei)
    !--pjr
    !
    !
    RETURN
  END SUBROUTINE cldefr


  !===============================================================================
  SUBROUTINE reltab(nCols,kMax, t, landm, icefrac, rel, snowh)
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
    INTEGER  , INTENT(in) :: nCols
    INTEGER, INTENT(in) :: kMax                  ! number of vertical levels

    REAL(r8), INTENT(in) :: icefrac(nCols)       ! Ice fraction
    REAL(r8), INTENT(in) :: snowh(nCols)         ! Snow depth over land, water equivalent (m)
    REAL(r8), INTENT(in) :: landm(nCols)         ! Land fraction ramping to zero over ocean
    REAL(r8), INTENT(in) :: t(nCols,kMax)        ! Temperature

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: rel(nCols,kMax)      ! Liquid effective drop size (microns)
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
    rliqice   = 14.0_r8
    !rliqocean = 14.0_r8
    rliqocean = 10.0_r8
    rliqland  = 8.0_r8
    DO k=1,kMax
       DO i=1,nCols
          ! jrm Reworked effective radius algorithm
          ! Start with temperature-dependent value appropriate for continental air
          ! Note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * MIN(1.0_r8,MAX(0.0_r8,(tmelt-t(i,k))*0.05_r8))
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
  SUBROUTINE reitab(nCols,kMax, t, re)
    !

    INTEGER  , INTENT(in) :: nCols
    INTEGER  , INTENT(in) :: kMax
    REAL(r8), INTENT(in ) :: t(nCols,kMax)
    REAL(r8), INTENT(out) :: re(nCols,kMax)
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
    DO k=1,kMax
       DO i=1,nCols
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


END MODULE CloudOpticalProperty
!PROGRAM Main
!   USE CloudOpticalProperty
!END PROGRAM MAin
