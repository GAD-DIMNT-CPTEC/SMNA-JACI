!
!  $Author: pkubota $
!  $Date: 2011/05/27 23:15:12 $
!  $Revision: 1.3 $
!
MODULE Rad_Clirad

  USE Constants, ONLY :  &
       pai,              & ! constant pi=3.1415926e0
       tf,               & ! Temperatura de congelamento (K)=273.16e0
       r8,i8

  USE Options, ONLY : &
       nfprt

  !
  ! cliradintf---
  !             |
  !             |
  !             |
  !             |
  !             |
  !             |
  !             |
  !             |
  !             | cloudy*--| soradcld ----| soluvcld ---| cldscale
  !                                       |             | deledd
  !                                       |             | cldflx
  !                                       |
  !                                       | solircld ---| cldscale
  !                                       |             | deledd
  !                                       |             | cldflx
  !                                       | rflx

  USE Options, ONLY :       &
       asolc, &
       asolm

  IMPLICIT NONE

  PRIVATE
 ! parameters for co2 transmission tables
 INTEGER, PARAMETER :: nu=43
 INTEGER, PARAMETER :: nw=37
 INTEGER, PARAMETER :: nx=62
 INTEGER, PARAMETER :: ny=101
 ! parameters
 INTEGER, PARAMETER :: nm=11
 INTEGER, PARAMETER :: nt=9
 INTEGER, PARAMETER :: na=11
 ! include the pre-computed table of mcai for scaling the cloud optical
 ! thickness under the assumption that clouds are maximally overlapped
 ! caib is for scaling the cloud optical thickness for direct radiation
 REAL(KIND=r8), ALLOCATABLE :: caib(:,:,:)
 ! caif is for scaling the cloud optical thickness for diffuse radiation
 REAL(KIND=r8), ALLOCATABLE :: caif(:,:)

 ! cah is the co2 absorptance in band 10
 REAL(KIND=r8), ALLOCATABLE  :: coa  (:,:)
 ! coa is the co2 absorptance in strong absorption regions of band 11
 REAL(KIND=r8), ALLOCATABLE :: cah (:,:)

 ! parameters
 INTEGER, PARAMETER :: nband=8

 ! hk is the fractional extra-terrestrial solar flux in each
 ! of the 8 bands. the sum of hk is 0.47074. (table 3)
 REAL(KIND=r8), PARAMETER, DIMENSION(nband) :: hk3 = (/ &
      0.00057_r8, 0.00367_r8, 0.00083_r8, 0.00417_r8, 0.00600_r8, &
      0.00556_r8, 0.05913_r8, 0.39081_r8 /)

 ! zk is the ozone absorption coefficient. unit: /(cm-atm)stp (table 3)
 REAL(KIND=r8), PARAMETER, DIMENSION(nband) :: zk3 = (/ &
      30.47_r8, 187.2_r8, 301.9_r8, 42.83_r8, 7.09_r8, 1.25_r8, 0.0345_r8, 0.0572_r8 /)

 ! wk is the water vapor absorption coefficient. unit: cm**2/g (table 3)
 REAL(KIND=r8), PARAMETER, DIMENSION(nband) :: wk3 = (/ 0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,&
                                                        0.0_r8,0.00075_r8/)

 ! ry is the extinction coefficient for rayleigh scattering. unit: /mb. (table 3)
 REAL(KIND=r8), PARAMETER, DIMENSION(nband) :: ry3 = (/ &
      0.00604_r8, 0.00170_r8, 0.00222_r8, 0.00132_r8, 0.00107_r8, 0.00091_r8, &
      0.00055_r8, 0.00012_r8 /)

 ! coefficients for computing the extinction coefficients of ice,
 ! water, and rain particles, independent of spectral band. (table 4)
 REAL(KIND=r8), PARAMETER, DIMENSION(2) :: aib3 = (/  3.33e-4_r8, 2.52_r8 /)
 REAL(KIND=r8), PARAMETER, DIMENSION(2) :: awb3 = (/ -6.59e-3_r8, 1.65_r8 /)
 REAL(KIND=r8), PARAMETER, DIMENSION(2) :: arb3 = (/  3.07e-3_r8, 0.00_r8 /)

 ! coefficients for computing the asymmetry factor of ice, water,
 ! and rain particles, independent of spectral band. (table 6)

 REAL(KIND=r8), PARAMETER, DIMENSION(3) :: aig3 = (/ 0.74625_r8, 0.0010541_r8, -0.00000264_r8 /)
 REAL(KIND=r8), PARAMETER, DIMENSION(3) :: awg3 = (/ 0.82562_r8, 0.0052900_r8, -0.00014866_r8 /)
 REAL(KIND=r8), PARAMETER, DIMENSION(3) :: arg3 = (/ 0.883_r8,   0.0_r8,      0.0_r8/)

  INTEGER ::    nbndsw=14
! may need to rename these - from v2.6
      integer, parameter :: jpband   = 29
      integer, parameter :: jpb1     = 16   !istart
      integer, parameter :: jpb2     = 29   !iend
      real(kind=r8) :: wavenum2(jpb1:jpb2)


!------------------------------------------------------------------
! rrtmg_sw cloud property coefficients
!
! Initial: J.-J. Morcrette, ECMWF, oct1999
! Revised: J. Delamere/MJIacono, AER, aug2005
! Revised: MJIacono, AER, nov2005
! Revised: MJIacono, AER, jul2006
!------------------------------------------------------------------
!
!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! xxxliq1 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from
!                    Hu & Stamnes, j. clim., 6, 728-742, 1993.  
! xxxice2 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from streamer v3.0,
!                    Key, streamer user's guide, cooperative institude 
!                    for meteorological studies, 95 pp., 2001.
! xxxice3 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from
!                    Fu, j. clim., 9, 1996.
! xbari   : real   : optical property coefficients for five spectral 
!                    intervals (2857-4000, 4000-5263, 5263-7692, 7692-14285,
!                    and 14285-40000 wavenumbers) following 
!                    Ebert and Curry, jgr, 97, 3831-3836, 1992.
!------------------------------------------------------------------

      real(kind=r8) :: extliq1(58,16:29), ssaliq1(58,16:29), asyliq1(58,16:29)
      real(kind=r8) :: extice2(43,16:29), ssaice2(43,16:29), asyice2(43,16:29)
      real(kind=r8) :: extice3(46,16:29), ssaice3(46,16:29), asyice3(46,16:29)
      real(kind=r8) :: fdlice3(46,16:29)
      real(kind=r8) :: abari(5),bbari(5),cbari(5),dbari(5),ebari(5),fbari(5)
      
 PUBLIC :: InitCliradSW
 PUBLIC :: cliradsw

CONTAINS
   SUBROUTINE InitCliradSW()
      IMPLICIT NONE 
      ALLOCATE(caib(nm,nt,na))
      ALLOCATE(caif(nt,na))
      ALLOCATE(coa (nx,ny))
      ALLOCATE(cah (nu,nw))

      CALL read_table()
      CALL swcldpr()
   END SUBROUTINE InitCliradSW


  !-----------------------------------------------------------------------
  !
  ! Subroutine: CLIRADSW (Clirad Interface)
  !
  ! This subroutine replaces SWRAD in the global model. It uses nearly the 
  ! same input and output and calls the main routines of Chou and Suarez
  ! radiation code.
  !
  ! ACRONYMS: 
  !   CDLGP...CLOUDY DAYTIME LATITUDE GRID POINTS
  !   DLGP....DAYTIME LATITUDE GRID POINTS
  !   LGP.....LATITUDE GRID POINTS
  !   NSOL....NUMBER OF DAYTIME LATITUDE GRID POINTS
  !   NCLD....NUMBER OF CLOUDY DLGP
  !
  ! Authors: Tatiana Tarasova & Henrique Barbosa
  !
  !-----------------------------------------------------------------------

  SUBROUTINE cliradsw(&
       ! Model Info and flags
       ncols , kmax  , nls   , noz   , &
       icld  , inalb , s0    , cosz  , &
       schemes     , &
       ! Atmospheric fields
       pl20  , dpl   , tl    , ql    , &
       o3l   , co2l  ,gps   , imask , &
       ! SURFACE:  albedo
       alvdf , alndf , alvdr , alndr , &
       ! SW Radiation fields 
       swinc ,                         &
       radvbc, radvdc, radnbc, radndc, &
       radvbl, radvdl, radnbl, radndl, &
       dswclr, dswtop, ssclr , ss    , &
       aslclr, asl   ,                 &
       ! Cloud field
       cld   , clu   , fice  , &
       rei   , rel   , taud  ,  &
       cicewp_in,cliqwp_in , &
       E_cld_tau_in,    E_cld_tau_w_in  ,&
       E_cld_tau_w_g_in,E_cld_tau_w_f_in ,cldfprime_in  )

    IMPLICIT NONE
  !==========================================================================
  !
  ! ___________________________
  ! MODEL INFORMATION AND FLAGS
  !
  ! ncols.....Number of grid points on a gaussian latitude circle
  ! kmax......Number of sigma levels
  ! nls.......Number of layers in the stratosphere.
  ! noz.......
  ! icld......Select three types of cloud emissivity/optical depth
  !           =1 : old cloud emissivity (optical depth) setting
  !                ccu :  0.05 *dp
  !                css :  0.025*dp       for ice cloud t<253.0
  !                       0.05 *dp       for ice cloud t>253.0
  !           =2 : new cloud emissivity (optical depth) setting
  !                ccu : (0.16)*dp
  !                css :  0.0                              t<-82.5
  !                      (2.0e-6*(t-tcrit)**2)*dp    -82.5<t<-10.0
  !                      (6.949e-3*(t-273)+.08)*dp   -10.0<t<  0.0
  !                      (0.08)*dp                   -10.0<t<  0.0
  !           =3    : ccm3 based cloud emissivity
  ! inalb.....Select two types of surface albedo
  !           =1 : input two  types surface albedo (2 diffused)
  !                direct beam albedos are calculated by the subr.
  !           =2 : input four types surfc albedo (2 diff,2 direct)
  ! s0........Solar constant  at proper sun-earth distance
  ! cosz......Cosine of solar zenith angle
  !
  ! __________________
  ! ATMOSPHERIC FIELDS
  !
  !    flip arrays (k=1 means top of atmosphere)
  !
  ! pl20......Flip array of pressure at bottom of layers (mb)
  ! dpl.......Flip array of pressure difference bettween levels
  ! tl........Flip array of temperature (K)
  ! ql........Flip array of specific humidity (g/g)
  ! o3l.......Flip array of ozone mixing ratio (g/g)
  ! gps.......Surface pressure (mb)
  ! imask.....Sea/Land mask
  !
  ! ________________
  ! SURFACE:  albedo
  !
  ! alvdf.....Visible diffuse surface albedo
  ! alndf.....Near-ir diffuse surface albedo
  ! alvdr.....Visible beam surface albedo
  ! alndr.....Near-ir beam surface albedo
  !
  ! ________________________________________
  ! SW Radiation fields 
  !
  ! swinc.....Incident SW at top 
  
  ! radvbc....Down Sfc SW flux visible beam    (all-sky)
  ! radnbc....Down Sfc SW flux Near-IR beam    (all-sky)
  ! radvdc....Down Sfc SW flux visible diffuse (all-sky)
  ! radndc....Down Sfc SW flux Near-IR diffuse (all-sky)
  
  ! radvbl....Down Sfc SW flux visible beam    (clear)  
  ! radnbl....Down Sfc SW flux Near-IR beam    (clear)  
  ! radvdl....Down Sfc SW flux visible diffuse (clear)  
  ! radndl....Down Sfc SW flux Near-IR diffuse (clear)  
  
  ! dswclr....Net SW flux at TOA (clear)   = Abs by Atmos + Sfc
  ! dswtop....Net SW flux at TOA (all-sky) = Abs by Atmos + Sfc
  ! ssclr.....Net SW flux at SFC (clear)   = Abs by Sfc
  ! ss........Net SW flux at SFC (all-sky) = Abs by Sfc
  
  ! asl.......Heating rate due to shortwave         (K/s)
  ! aslclr....Heating rate due to shortwave (clear) (K/s)
  !
  ! ___________
  ! Cloud field
  !
  ! cld.......Supersaturation cloud fraction
  ! clu.......Convective cloud fraction     
  ! fice......Fraction of cloud water in the form of ice
  ! rei.......Ice particle Effective Radius (microns)
  ! rel.......Liquid particle Effective Radius (microns)
  ! taud......Shortwave cloud optical depth
  !
  !==========================================================================

    ! Model Info and flags
    INTEGER      ,    INTENT(in   ) :: ncols
    INTEGER      ,    INTENT(in   ) :: kmax
    INTEGER      ,    INTENT(in   ) :: nls
    LOGICAL      ,    INTENT(in   ) :: noz
    INTEGER      ,    INTENT(in   ) :: icld
    INTEGER      ,    INTENT(in   ) :: inalb
    REAL(KIND=r8),    INTENT(in   ) :: s0
    REAL(KIND=r8),    INTENT(in   ) :: cosz   (ncols)
    INTEGER      ,    INTENT(IN   ) :: schemes
    ! Atmospheric fields
    REAL(KIND=r8),    INTENT(in   ) :: pl20   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: dpl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tl     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: ql     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: o3l    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: co2l   (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gps    (ncols)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask  (ncols) 

    ! SURFACE:  albedo
    REAL(KIND=r8),    INTENT(in   ) :: alvdf  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: alndf  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: alvdr  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: alndr  (ncols)
    ! SW Radiation fields 
    REAL(KIND=r8),    INTENT(out) :: swinc  (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvbc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvdc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radnbc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radndc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvbl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvdl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radnbl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radndl (ncols)
    REAL(KIND=r8),    INTENT(out) :: aslclr (ncols,kmax)
    REAL(KIND=r8),    INTENT(out) :: asl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(out) :: ss     (ncols)
    REAL(KIND=r8),    INTENT(out) :: ssclr  (ncols)
    REAL(KIND=r8),    INTENT(out) :: dswtop (ncols)
    REAL(KIND=r8),    INTENT(out) :: dswclr (ncols)

    ! Cloud field and Microphysics
    REAL(KIND=r8),    INTENT(in   ) :: cld    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: clu    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: fice   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: rei    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: rel    (ncols,kmax) 
    REAL(KIND=r8),    INTENT(inout) :: taud   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) ::  cicewp_in(ncols,kmax) 
    REAL(KIND=r8),    INTENT(in   ) ::  cliqwp_in(ncols,kmax) 
    REAL(KIND=r8),    INTENT(in   ) ::  cldfprime_in (ncols,kmax)

    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_in    (nbndsw, nCols, kMax)      ! cloud optical depth
    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_w_in  (nbndsw, nCols, kMax)      ! cloud optical 
    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_w_g_in(nbndsw, nCols, kMax)      ! cloud optical 
    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_w_f_in(nbndsw, nCols, kMax)      ! cloud optical 

    ! local  variables
    INTEGER i,k,ns
    INTEGER :: nsol ! NUMBER OF SOLAR LATITUDE GRID POINTS (COSZ>0.01)

    REAL(KIND=r8):: scosz (ncols) ! DOWNWARD FLUX AT TOP IN DLGP
    REAL(KIND=r8):: cmu   (ncols) ! COSINE OF SOLAR ZENITH ANGLE IN DLGP
    INTEGER      :: dmask (ncols) ! sib-mask in DLGP
    !
    !     Downward ground fluxes in dlgp
    !
    REAL(KIND=r8):: rvbl  (ncols) ! VISIBLE, BEAM, CLEAR
    REAL(KIND=r8):: rvdl  (ncols) ! VISIBLE, DIFFUSE, CLEAR
    REAL(KIND=r8):: rnbl  (ncols) ! NearIR, BEAM, CLEAR
    REAL(KIND=r8):: rndl  (ncols) ! NearIR, DIFFUSE, CLEAR
    
    REAL(KIND=r8):: rvbc  (ncols) ! VISIBLE, BEAM, CLOUDY
    REAL(KIND=r8):: rvdc  (ncols) ! VISIBLE, DIFFUSE, CLOUDY
    REAL(KIND=r8):: rnbc  (ncols) ! NearIR, BEAM, CLOUDY
    REAL(KIND=r8):: rndc  (ncols) ! NearIR, DIFFUSE, CLOUDY
    
    REAL(KIND=r8):: sl    (ncols) ! NET GROUND CLEAR FLUX IN DLGP
    REAL(KIND=r8):: sc    (ncols) ! NET GROUND CLOUDY FLUX IN DLGP
    !
    !     Other fluxes
    !
    REAL(KIND=r8):: dsclr (ncols)        ! ABSORPTION OF CLEAR ATMOSPHERE AND GROUND
    REAL(KIND=r8):: dscld (ncols)        ! DOWNWARD FLUX AT TOP IN DLGP
    REAL(KIND=r8):: aclr  (ncols,kmax+1) ! ABSORPTION IN CLEAR ATM.
    REAL(KIND=r8):: acld  (ncols,kmax+1) ! HEATING RATE (CLOUDY)
    !
    !     Albedo in DLGP
    !
    REAL(KIND=r8):: agv   (ncols)        ! GROUND VISIBLE DIFFUSE ALBEDO IN DLGP
    REAL(KIND=r8):: agn   (ncols)        ! GROUND NEAR IR DIFFUSE ALBEDO IN DLGP
    REAL(KIND=r8):: rsurfv(ncols)        ! GROUND VISIBLE BEAM ALBEDO IN DLGP
    REAL(KIND=r8):: rsurfn(ncols)        ! GROUND NEAR IR BEAM ALBEDO IN DLGP
    !
    !     Atmospheric state in DLGP
    !
    REAL(KIND=r8):: ps (ncols)        ! surface pressure in DLGP (mb)
    REAL(KIND=r8):: ta (ncols,kmax+1) ! LAYER TEMPERATURE IN DLGP
    REAL(KIND=r8):: wa (ncols,kmax+1) ! LAYER SPECIFIC HUMIDITY IN DLGP
    REAL(KIND=r8):: oa (ncols,kmax+1) ! LAYER OZONE MIXING RATIO IN DLGP
    REAL(KIND=r8):: pu (ncols,kmax+2) ! PRESSURE AT BOTTOM OF LAYER IN DLGP
    REAL(KIND=r8):: dp (ncols,kmax+1) ! PRESSURE DIFFERENCE IN DLGP
    !
    !     Clouds in DLGP
    !
    REAL(KIND=r8):: css    (ncols,kmax+1) ! LARGE SCALE CLOUD AMOUNT IN DLGP
    REAL(KIND=r8):: ccu    (ncols,kmax+1) ! CUMULUS CLOUD AMOUNT IN DLGP
    REAL(KIND=r8):: e0     (ncols,kmax+1) ! Cloud optical depth
    REAL(KIND=r8):: frcice (ncols,kmax+1) ! Fraction of cloud water in the form of ice
    REAL(KIND=r8):: radliq (ncols,kmax+1) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8):: radice (ncols,kmax+1) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8):: CO2    (nCols,kMax+1) ! mol/mol                
    REAL(KIND=r8):: cldfrac   (nCols,1:kMax)
    REAL(KIND=r8):: tauc_sw   (1:nbndsw,nCols,1:kMax) 
    REAL(KIND=r8):: ssac_sw   (1:nbndsw,nCols,1:kMax) 
    REAL(KIND=r8):: asmc_sw   (1:nbndsw,nCols,1:kMax)  
    REAL(KIND=r8):: fsfc_sw   (1:nbndsw,nCols,1:kMax)  
    REAL(KIND=r8):: cicewp      (nCols,1:kMax)
    REAL(KIND=r8):: cliqwp      (nCols,1:kMax)
    REAL(KIND=r8):: reicmcl     (nCols,1:kMax)
    REAL(KIND=r8):: relqmcl     (nCols,1:kMax)
    REAL(KIND=r8):: taucldorig  (1:kMax,1:jpband) 
    REAL(KIND=r8):: taucloud    (1:kMax,1:jpband)
    REAL(KIND=r8):: TauCloudTot (ncols,kmax) !
    REAL(KIND=r8):: tauk        (ncols,kmax+1) ! total optical depth with k=1 at top
    !
    !     Logical and working vectors
    !
    INTEGER  :: litx  (ncols) ! NUMBERS OF DLGP IN ALL LAYERS
    INTEGER  :: listim(ncols) ! =1,2,3...NCOLS*(KMAX+1)
    LOGICAL  :: bitx  (ncols) ! TRUE IN SOLAR LATITUDE GRID POINTS
    LOGICAL  :: bitc  (ncols,kmax+1) ! Working logical

    ! Initialize local vectors and output variables
    scosz  =0.0_r8
    cmu    =0.0_r8
    dmask  =0
    rvbl   =0.0_r8
    rvdl   =0.0_r8
    rnbl   =0.0_r8
    rndl   =0.0_r8
    rvbc   =0.0_r8
    rvdc   =0.0_r8
    rnbc   =0.0_r8
    rndc   =0.0_r8
    sl     =0.0_r8
    sc     =0.0_r8
    dsclr =0.0_r8
    dscld =0.0_r8
    aclr  =0.0_r8
    acld  =0.0_r8
    cldfrac    =0.0_r8
    tauc_sw    =0.0_r8
    ssac_sw    =0.0_r8
    asmc_sw    =0.0_r8
    fsfc_sw    =0.0_r8
    cicewp     =0.0_r8
    cliqwp     =0.0_r8
    reicmcl    =0.0_r8
    relqmcl    =0.0_r8
    taucldorig =0.0_r8
    taucloud   =0.0_r8
    agv    =0.0_r8
    agn    =0.0_r8
    rsurfv =0.0_r8
    rsurfn =0.0_r8
    ps    =0.0_r8
    ta    =0.0_r8
    wa    =0.0_r8
    oa    =0.0_r8
    pu    =0.0_r8
    dp    =0.0_r8
    css   =0.0_r8
    ccu   =0.0_r8
    e0    =0.0_r8
    frcice=0.0_r8
    radliq=0.0_r8
    radice=0.0_r8
    CO2   =0.0_r8
    nsol  =0
    tauk  =0.0_r8
    litx  =0
    listim=0
    bitx  =.FALSE.
    bitc  =.FALSE.
    swinc =0.0_r8
    radvbc=0.0_r8
    radvdc=0.0_r8
    radnbc=0.0_r8
    radndc=0.0_r8
    radvbl=0.0_r8
    radvdl=0.0_r8
    radnbl=0.0_r8
    radndl=0.0_r8
    aslclr=0.0_r8
    asl   =0.0_r8
    ss    =0.0_r8
    ssclr =0.0_r8
    dswtop=0.0_r8
    dswclr=0.0_r8

    scosz=0.0_r8
    cmu=0.0_r8
    listim=0
    bitx=.FALSE.

    ps=0.0_r8
    ta=0.0_r8
    wa=0.0_r8
    oa=0.0_r8
    pu=0.0_r8

    rsurfv=0.0_r8
    rsurfn=0.0_r8
    agv=0.0_r8
    agn=0.0_r8

    ccu=0.0_r8
    css=0.0_r8
    dp=0.0_r8
    TauCloudTot=0.0_r8
    litx=0
    dmask=0

    ! Subroutine starts here

    !
    ! Set array listim = i, WHEN I=1,ncols*(kmax+2)
    !
    !FORALL (I=1:ncols) listim(I)=I
    DO i=1,ncols
        listim(i)=i
    END DO
    !
    ! set bits for daytime grid points
    ! BITX=.TRUE. IF COSZ>DUM(1)....0.01
    !
    !bitx(1:ncols)=cosz(1:ncols).ge.0.01e0_r8
    nsol=0
    DO i=1,ncols
       IF(cosz(i).ge.0.01e0_r8)THEN
           nsol=nsol+1
           bitx(i)=.TRUE.
       END IF
    END DO
    !
    ! Calculate nsol = number of daytime latitude grid points
    !
    !nsol=COUNT(bitx(1:ncols))
    !
    ! Set zero to all latitude grids surface fluxes
    !
    swinc  = 0.0_r8
    ss     = 0.0_r8
    ssclr  = 0.0_r8
    dswtop = 0.0_r8
    dswclr = 0.0_r8
    radvbl = 0.0_r8
    radvdl = 0.0_r8
    radnbl = 0.0_r8
    radndl = 0.0_r8
    radvbc = 0.0_r8
    radvdc = 0.0_r8
    radnbc = 0.0_r8
    radndc = 0.0_r8
    asl    = 0.0_r8
    aslclr = 0.0_r8
    !
    ! If there are no daytime points then
    !
    IF(nsol.eq.0) RETURN
    !
    ! Set integer array litx (nsol*(kmax+1))
    ! numbers of latitude daytime grid points at first level
    !
    !litx(1:nsol) = PACK(listim(1:ncols), bitx(1:ncols))
    nsol=0
    DO  i=1,ncols
       IF(bitx(i))THEN
          nsol=nsol+1
          litx(nsol) = listim(i)
       END IF  
    END DO
    !
    !  Transform  two-size  input arrays:
    !        pl20(ncols,kmax),dpl,tl,ql,cld,clu
    !  in two-size arrays:
    !        pu(nsol,kmax+2),dp(nsol,kmax+1),ta,wa,css,ccu
    !  in daytime latitude grid points at all levels
    !
    
    DO k=2,kMax+1
       DO i=1,nsol
          IF(litx(I).le.ncols)THEN
          pu (i,k+1)=pl20(litx(i),k-1)
          dp (i,k)  =dpl (litx(i),k-1)
          ta (i,k)  =tl  (litx(i),k-1)
          co2(i,k)  =co2l(litx(i),k-1)
          wa (i,k)  =max (0.1e-22_r8,ql(litx(i),k-1))
          oa (I,k)  =max (0.1e-9_r8,o3l(litx(i),k-1))
          css(i,k)  =cld (litx(i),k-1)
          ccu(i,k)  =clu (litx(i),k-1)
          ps (i)    =gps (litx(i))
          END IF
       END DO
    END DO
    !
    ! Set some parameters at first, second
    !
    !#TO Correcao de hmjb
    DO i=1,nsol
       pu(i,2) = PU(I,3)/2.0_r8  ! 0.5_r8
       pu(i,1) = PU(I,2)/10.0_r8 ! 0.05_r8
       dp(i,1) = pu(i,2)-pu(i,1) ! pressure differense
       ta(i,1) = ta(i,2)       ! temperature
       wa(i,1) = wa(i,2)       ! specific humidity
       oa(i,1) = oa(i,2)       ! ozone
       co2(i,1) = co2(i,2)
    ENDDO

    ! if no ozone
    IF (noz) oa=0.0_r8

    !
    ! Set some parameters at stratospheric levels
    !
    css(1:nsol,1:nls+1)=0.0_r8
    ccu(1:nsol,1:nls+1)=0.0_r8
    wa (1:nsol,1:nls  )=3.0e-6_r8

    !
    ! The same transformation as mentioned above for:
    !   visible surface albedo....alvdf to agv
    !   nearir surface albedo.....alndf to agn
    !   cosine of solar zenith angle..cosz to cmu
    !
    FORALL (I=1:nsol,litx(I).le.ncols)
       agv  (I) = alvdf (litx(I))
       agn  (I) = alndf (litx(I))
       cmu  (I) = cosz  (litx(I))
       dmask(I) = INT(imask (litx(i)))
    ENDFORALL
    !
    ! If direct beam albedos are given then
    ! alvdr transform to rsurfv(nsol) and alndr to rsurfn(nsol)
    ! in daytime grid points
    !
    !hmjb inalb=2 is hardcoded in physics!!!
    IF (inalb .eq. 2) THEN
       FORALL (I=1:nsol,litx(I).le.ncols)
          rsurfv(I)=alvdr(litx(I))
          rsurfn(I)=alndr(litx(I))
       ENDFORALL
    ELSE
       !
       ! If direct beam albedos are not given then do the reverse
       ! calculate direct beam surface albedo
       !
       rvbl(1:nsol)=acos(cmu(1:nsol)) ! rvbl... solar zenith angle
       DO i=1,nsol
          rvdc(i)  =  -18.0_r8 * (0.5_r8 * pai - rvbl(i)) / pai
          rvbc(i)  =  exp(rvdc(i))
       ENDDO
       DO i=1,nsol
          rvdc(i)  = (agv(i) - 0.054313_r8) / 0.945687_r8
          rndc(i)  = (agn(i) - 0.054313_r8) / 0.945687_r8
          rsurfv(i) = rvdc(i)+(1.0-rvdc(i))*rvbc(i)
          rsurfn(i) = rndc(i)+(1.0-rndc(i))*rvbc(i)
       ENDDO
       DO i=1,ncols
          alvdr(i) = 0.0_r8
          alndr(i) = 0.0_r8
       ENDDO
       FORALL (I=1:nsol,litx(I).le.ncols) alvdr(litx(I))=rsurfv(I)
       FORALL (I=1:nsol,litx(I).le.ncols) alndr(litx(I))=rsurfn(I)
    ENDIF
    !
    ! CMU.......COSINE OF SOLAR ZENITH ANGLE AT DLGP
    !
    DO i=1,nsol
       scosz(i) = s0 * cmu(i)  ! DOWNWARD SOLAR FLUX AT TOP
    ENDDO
    !
    ! Transform scosz(nsol) to swinc(ncols) at all lgp
    !
    FORALL(I=1:nsol,litx(I).le.ncols) swinc(litx(I))=scosz(I)
    !
    ! Calculate solar fluxes
    ! Calls cloudy to calculate:
    !   dsclr,sl,aclr,rvbl,rvdl,rnbl,rndl
    !   dscld,sc,acld,rvbc,rvdc,rnbc,rndc
    ! The values are packed at the begining of the arrays.
    ! Instead of occupying 1..ncols, they cover only the range 1..nsol
    !

    ! Set cloud amount as Maximum cloudiness
    css=max(ccu,css)

    IF (icld.eq.1) THEN

       e0  (1:nsol,1:kmax+1) = 0.05_r8
       bitc(1:nsol,1:kmax+1) = (ta(1:nsol,1:kmax+1).lt.253.0_r8).and.(ccu(1:nsol,1:kmax+1).eq.0.0_r8)

       ! IF BITC=.TRUE. EO=0.025_r8
       WHERE (bitc(1:nsol,1:kmax+1)) e0(1:nsol,1:kmax+1)=0.025_r8

       ! the extra cloud fraction idea from ncar is not used with clirad
       ! because clirad properly acounts for combination between different
       ! layers of clouds
       !WHERE(css.gt.0.0_r8) e0=e0*dp*css
       WHERE(css(1:nsol,1:kmax+1).gt.0.0_r8) e0(1:nsol,1:kmax+1)=e0(1:nsol,1:kmax+1)*dp(1:nsol,1:kmax+1)
       ! Prepare tau, fice, rel e rei
       DO k=1,kmax
          DO i=1,nsol
             IF(litx(I).le.ncols)THEN
                taud(litx(i),k)= e0  (i,k+1) 
             END IF
          END DO
       END DO

    ELSE IF(icld.eq.4 .or. icld == 5.or. icld == 6 .or. icld == 7) THEN
       IF(icld == 6)THEN

       DO k=1,kmax
          DO i=1,nsol
             IF(litx(I).le.ncols)THEN
                reicmcl(i,k)   = rei (litx(i),k)
                relqmcl(i,k)   = rel (litx(i),k)
                cicewp (i,k)   = cicewp_in(litx(i),k)
                cliqwp (i,k)   = cliqwp_in(litx(i),k)
                cldfrac(i,k)   = cldfprime_in(litx(i),k)
             END IF
          END DO
       END DO
       do ns = 1, nbndsw
          DO k=1,kmax
             DO i=1,nsol
                IF(litx(I).le.ncols)THEN
                   if (E_cld_tau_w_in(ns,litx(i),k) > 0._r8) then
                      fsfc_sw(ns,i,k)=E_cld_tau_w_f_in(ns,litx(i),k)/E_cld_tau_w_in(ns,litx(i),k)
                      asmc_sw(ns,i,k)=E_cld_tau_w_g_in(ns,litx(i),k)/E_cld_tau_w_in(ns,litx(i),k)
                   else
                      fsfc_sw(ns,i,k) = 0._r8
                      asmc_sw(ns,i,k) = 0._r8
                   endif
                   tauc_sw(ns,i,k)=E_cld_tau_in(ns,litx(i),k)
                   if (tauc_sw(ns,i,k) > 0._r8) then
                     ssac_sw(ns,i,k)=E_cld_tau_w_in(ns,litx(i),k)/tauc_sw(ns,i,k)
                   else
                      tauc_sw(ns,i,k) = 0._r8
                      fsfc_sw(ns,i,k) = 0._r8
                      asmc_sw(ns,i,k) = 0._r8
                      ssac_sw(ns,i,k) = 1._r8
                   endif
                END IF
             END DO
          END DO
       END DO
       DO i=1,nsol
          call cldprop_sw( &             !             cldprop_sw(
                            kmax                              , &!          kmax   , &
                            cldfrac (           i,1:kmax)     , &!          cldfrac   (1:kmax)                 , &
                            tauc_sw   (1:nbndsw,i,1:kmax)     , &!          tauc_sw   (1:nbndsw,iplon,1:kmax)  , &
                            ssac_sw   (1:nbndsw,i,1:kmax)     , &!          ssac_sw   (1:nbndsw,iplon,1:kmax)  , &
                            asmc_sw   (1:nbndsw,i,1:kmax)     , &!          asmc_sw   (1:nbndsw,iplon,1:kmax)  , &
                            cicewp    (         i,1:kmax)     , &!          cicewp    (iplon,1:kmax)               , &
                            cliqwp    (         i,1:kmax)     , &!          cliqwp    (iplon,1:kmax)               , &
                            reicmcl   (         i,1:kmax)     , &!          reicmcl   (iplon,1:kmax)               , &
                            2*reicmcl (         i,1:kmax)     , &!          2*reicmcl (iplon,1:kmax)               , &
                            relqmcl   (         i,1:kmax)     , &!          relqmcl   (iplon,1:kmax)               , &
                            taucldorig(1:kmax,1:jpband)       , &!          taucldorig(1:kmax,1:jpband)          , &
                            taucloud  (1:kmax,1:jpband)         )!          taucloud  (1:kmax,1:jpband)          , &
          DO k=1,kmax
             DO ns = 1, jpband
                 TauCloudTot(i,k)=TauCloudTot(i,k)+taucloud (k,ns)
             END DO
          ENDDO
       END DO
       END IF
       ! Prepare tau, fice, rel e rei
       DO k=1,kmax
          DO i=1,nsol
             IF(litx(I).le.ncols)THEN
             IF(icld == 6)THEN
                tauk  (i,k+1) = TauCloudTot(i,k)
!                tauk  (i,k+1) = MAX(MIN((FlipTauCloudTotSW(litx(i),k)+taud(litx(i),k))/2.0_r8,500.0_r8),1.0e-8_r8)
             ELSE
                tauk  (i,k+1) = taud(litx(i),k)
             END IF
             frcice(i,k+1) = fice(litx(i),k)
             radice(i,k+1) = rei (litx(i),k)
             radliq(i,k+1) = rel (litx(i),k)
             END IF
          END DO
       END DO
       tauk  (1:nsol,1)=0.0_r8
       frcice(1:nsol,1)=frcice(1:nsol,2)
       radice(1:nsol,1)=radice(1:nsol,2)
       radliq(1:nsol,1)=radliq(1:nsol,2)

       ! the extra cloud fraction idea from ncar is not used with clirad
       !e0(1:nsol,1:kmax+1) = sqrt(css(1:nsol,1:kmax+1))*tauk(1:nsol,1:kmax+1)
       e0(1:nsol,1:kmax+1) = tauk(1:nsol,1:kmax+1)
    ELSE
       WRITE(nfprt,*) 'error! icld must be 1 ,4 ,6 , 5 ,6 ,7 with Clirad-sw-m '
       STOP
    ENDIF

    !
    ! Call subroutine cloudy to calculate all-sky fluxes
    !
    CALL cloudy( &
         schemes                ,& !INTEGER      , INTENT(IN   ) :: schemes
         s0                     ,& !REAL(KIND=r8), INTENT(IN   ) :: s0
         CO2 (1:nsol,1:kmax+1)  ,& !REAL(KIND=r8), INTENT(IN   ) :: rco2(m,np)
         nsol                   ,& !INTEGER      , INTENT(IN   ) :: m
         kmax+1                 ,& !INTEGER      , INTENT(IN   ) :: np
         pu(1:nsol,1:kmax+2)    ,& !REAL(KIND=r8), INTENT(IN   ) :: pl      (m,np+1)
         ta(1:nsol,1:kmax+1)    ,& !REAL(KIND=r8), INTENT(IN   ) :: ta      (m,np  )
         wa(1:nsol,1:kmax+1)    ,& !REAL(KIND=r8), INTENT(IN   ) :: wa      (m,np  )
         oa(1:nsol,1:kmax+1)    ,& !REAL(KIND=r8), INTENT(IN   ) :: oa      (m,np  )
         cmu(1:nsol)            ,& !REAL(KIND=r8), INTENT(IN   ) :: cosz    (m)  
         rsurfv(1:nsol)         ,& !REAL(KIND=r8), INTENT(IN   ) :: rsuvbm  (m)       
         agv(1:nsol)            ,& !REAL(KIND=r8), INTENT(IN   ) :: rsuvdf  (m)       
         rsurfn(1:nsol)         ,& !REAL(KIND=r8), INTENT(IN   ) :: rsirbm  (m)       
         agn(1:nsol)            ,& !REAL(KIND=r8), INTENT(IN   ) :: rsirdf  (m)       
         dscld(1:nsol)          ,& !REAL(KIND=r8), INTENT(OUT  ) :: dscld1  (m)   
         sc(1:nsol)             ,& !REAL(KIND=r8), INTENT(OUT  ) :: sc1     (m)   
         acld(1:nsol,1:kmax+1)  ,& !REAL(KIND=r8), INTENT(OUT  ) :: acld1   (m,np)  
         rvbc(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rvbc1   (m)       
         rvdc(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rvdc1   (m)       
         rnbc(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rnbc1   (m)       
         rndc(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rndc1   (m)       
         dsclr(1:nsol)          ,& !REAL(KIND=r8), INTENT(OUT  ) :: dsclr1  (m)     
         sl(1:nsol)             ,& !REAL(KIND=r8), INTENT(OUT  ) :: sl1     (m)       
         aclr(1:nsol,1:kmax+1)  ,& !REAL(KIND=r8), INTENT(OUT  ) :: aclr1   (m,np)  
         rvbl(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rvbl1   (m)  
         rvdl(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rvdl1   (m)  
         rnbl(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rnbl1   (m)  
         rndl(1:nsol)           ,& !REAL(KIND=r8), INTENT(OUT  ) :: rndl1   (m)  
         e0(1:nsol,1:kmax+1)    ,& !REAL(KIND=r8), INTENT(IN   ) :: tauc    (m,np  )
         css(1:nsol,1:kmax+1)   ,& !REAL(KIND=r8), INTENT(IN   ) :: csscgp  (m,np  )
         ps(1:nsol)             ,& !REAL(KIND=r8), INTENT(IN   ) :: psc     (m)     
         dmask(1:nsol)          ,& !INTEGER      , INTENT(IN   ) :: dmask   (m) !sib-mask in DLGP
         frcice(1:nsol,1:kmax+1),& !REAL(KIND=r8), intent(in   ) :: fice    (m,np) 
         radliq(1:nsol,1:kmax+1),& !REAL(KIND=r8), intent(in   ) :: rel     (m,np)  
         radice(1:nsol,1:kmax+1),& !REAL(KIND=r8), intent(in   ) :: rei     (m,np)  
         icld                    ) !INTEGER      , INTENT(IN   ) :: icld ! new cloud microphysics
 

    !
    ! SET SOLAR FLUXES IN ALL GRID POINTS
    ! All values are nsol-packed and need to be unpacked
    ! This is done by copying values from positions (1:nsol) to
    ! positions litx(1:nsol).
    !
    FORALL(I=1:nsol,litx(I).le.ncols)
       ! clear
       ssclr (litx(I))=sl   (I)
       dswclr(litx(I))=dsclr(I)
       radvbl(litx(I))=rvbl (I)
       radvdl(litx(I))=rvdl (I)
       radnbl(litx(I))=rnbl (I)
       radndl(litx(I))=rndl (I)

       ! cloudy
       ss    (litx(I))=sc   (I)
       dswtop(litx(I))=dscld(I)
       radvbc(litx(I))=rvbc (I)
       radvdc(litx(I))=rvdc (I)
       radnbc(litx(I))=rnbc (I)
       radndc(litx(I))=rndc (I)
    ENDFORALL

    DO k=1,kmax
       DO i=1,nsol
          IF(litx(I).le.ncols)THEN
            aslclr(litx(i),k)=aclr(i,k+1)
            asl   (litx(i),k)=acld(i,k+1)
          END IF
       END DO
    END DO


    !
    ! Calculation of solar heating rate in k/s
    !
    DO k=1,kmax
       DO i=1,ncols
          IF(aslclr(i,k).lt.1.0e-22_r8) aslclr(i,k) = 0.0_r8
          IF(asl   (i,k).lt.1.0e-22_r8) asl   (i,k) = 0.0_r8

          aslclr   (i,k) = aslclr(i,k) * 1.1574e-5_r8
          asl      (i,k) = asl   (i,k) * 1.1574e-5_r8
       ENDDO
    ENDDO

  END SUBROUTINE cliradsw

  !
  ! Subroutine: CLOUDY
  !
  ! $Author: pkubota $
  ! Modifications: H.Barbosa 2005
  !
  ! Description:
  !
  !NEW! continental aerosol model is included
  !NEW! the k-distributions of Tarasova and Fomin (2000)
  !NEW! 28 layers
  !
  !this  is  the source  program  for  computing  solar fluxes  due  to
  !absorption  by water  vapor, ozone,  co2, o2,  clouds,  and aerosols
  !anddue to scattering by clouds, aerosols, and gases.
  !
  !this is a vectorized code.   it computes fluxes simultaneously for m
  !soundings.
  !
  !the meaning, units and DIMENSION  of the input and output parameters
  !are given in the SUBROUTINE sorad.
  !
  ! Inputs:
  !setsw (global)         Clirad
  !
  !ncld       ncld
  !ncols
  !kmax
  !
  !puu(ncols*(kmax+2))....level pressures mb       pl
  !taa(ncols*(kmax+1))....layer temp in K       ta
  !waa(ncols*(kmax+1))....layer specific humidity g/g      wa
  !oaa(ncols*(kmax+1))....layer ozone concentration g/g    oa
  !tauc(ncols*(kmax+1))...cloud optical depth
  !css(ncols*(kmax+1))....cloud amount       fcld
  !cmu(ncols).............cosine solar zenith angle        cosz
  !rsurfv(ncols)..........Vis Beam Albedo       rsuvbm
  !agv(ncols).............Vis Diffuse Albedo       rsuvdf
  !rsurfn(ncols)..........Nir Beam Albedo       rsirbm
  !agn(ncols).............Nir Diffuse Albedo       rsirdf
  ! psc   =  surface pressure   (mb)
  ! Outputs:
  !  dscld(ncols)    ABSORPTION IN THE CLOUDY ATMOSPHERE AND AT THE GROUND
  !  sc(ncols)    ABSORPTION AT THE GROUND IN CLOUDY CASE
  !acld(ncols*(kmax+1)) HEATING RATE (CLOUDY) in K/sec
  !rvbc(ncols)    VISIBLE BEAM  FLUXES (CLOUDY)
  !rvdc(ncols)    VISIBLE DIFFUSE FLUXES (CLOUDY)
  !rnbc(ncols)    NEAR-IR BEAM  FLUXES (CLOUDY)
  !rndc(ncols)    NEAR-IR DIFFUSE FLUXES  (CLOUDY)
  !
  !
  !
  !
  SUBROUTINE cloudy( &
       schemes ,&!INTEGER      , INTENT(IN   ) :: schemes
       s0      ,&!REAL(KIND=r8), INTENT(IN   ) :: s0
       rco2    ,&!REAL(KIND=r8), INTENT(IN   ) :: rco2(m,np)
       m       ,&!INTEGER      , INTENT(IN   ) :: m
       np      ,&!INTEGER      , INTENT(IN   ) :: np
       pl      ,&!REAL(KIND=r8), INTENT(IN   ) :: pl      (m,np+1)
       ta      ,&!REAL(KIND=r8), INTENT(IN   ) :: ta      (m,np  )
       wa      ,&!REAL(KIND=r8), INTENT(IN   ) :: wa      (m,np  )
       oa      ,&!REAL(KIND=r8), INTENT(IN   ) :: oa      (m,np  )
       cosz    ,&!REAL(KIND=r8), INTENT(IN   ) :: cosz    (m)  
       rsuvbm  ,&!REAL(KIND=r8), INTENT(IN   ) :: rsuvbm  (m)       
       rsuvdf  ,&!REAL(KIND=r8), INTENT(IN   ) :: rsuvdf  (m)       
       rsirbm  ,&!REAL(KIND=r8), INTENT(IN   ) :: rsirbm  (m)       
       rsirdf  ,&!REAL(KIND=r8), INTENT(IN   ) :: rsirdf  (m)       
       dscld1  ,&!REAL(KIND=r8), INTENT(OUT  ) :: dscld1  (m)   
       sc1     ,&!REAL(KIND=r8), INTENT(OUT  ) :: sc1     (m)   
       acld1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: acld1   (m,np)  
       rvbc1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rvbc1   (m)       
       rvdc1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rvdc1   (m)       
       rnbc1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rnbc1   (m)       
       rndc1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rndc1   (m)       
       dsclr1  ,&!REAL(KIND=r8), INTENT(OUT  ) :: dsclr1  (m)     
       sl1     ,&!REAL(KIND=r8), INTENT(OUT  ) :: sl1     (m)       
       aclr1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: aclr1   (m,np)  
       rvbl1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rvbl1   (m)  
       rvdl1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rvdl1   (m)  
       rnbl1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rnbl1   (m)  
       rndl1   ,&!REAL(KIND=r8), INTENT(OUT  ) :: rndl1   (m)  
       tauc    ,&!REAL(KIND=r8), INTENT(IN   ) :: tauc    (m,np  )
       csscgp  ,&!REAL(KIND=r8), INTENT(IN   ) :: csscgp  (m,np  )
       psc     ,&!REAL(KIND=r8), INTENT(IN   ) :: psc     (m)     
       dmask   ,&!INTEGER      , INTENT(IN   ) :: dmask   (m) !sib-mask in DLGP
       fice    ,&!REAL(KIND=r8), intent(in   ) :: fice    (m,np) 
       rel     ,&!REAL(KIND=r8), intent(in   ) :: rel     (m,np)  
       rei     ,&!REAL(KIND=r8), intent(in   ) :: rei     (m,np)  
       icld     )!INTEGER      , INTENT(IN   ) :: icld ! new cloud microphysics

    IMPLICIT NONE

    ! input variables
    INTEGER      , INTENT(IN   ) :: schemes
    INTEGER      , INTENT(IN   ) :: m
    INTEGER      , INTENT(IN   ) :: np
    INTEGER      , INTENT(IN   ) :: dmask (m) !sib-mask in DLGP
    INTEGER      , INTENT(IN   ) :: icld ! new cloud microphysics
    REAL(KIND=r8), INTENT(IN   ) :: s0
    REAL(KIND=r8), INTENT(IN   ) :: rco2(m,np)
    ! pl,ta,wa and oa are, respectively, the level pressure (mb), layer
    ! temperature (k), layer specific humidity (g/g), and layer ozone
    ! concentration (g/g)
    REAL(KIND=r8), INTENT(IN   ) :: pl      (m,np+1)
    REAL(KIND=r8), INTENT(IN   ) :: ta      (m,np  )
    REAL(KIND=r8), INTENT(IN   ) :: wa      (m,np  )
    REAL(KIND=r8), INTENT(IN   ) :: oa      (m,np  )
    REAL(KIND=r8), INTENT(IN   ) :: tauc    (m,np  )
    REAL(KIND=r8), INTENT(IN   ) :: cosz    (m)     
    REAL(KIND=r8), INTENT(IN   ) :: rsuvbm  (m)     
    REAL(KIND=r8), INTENT(IN   ) :: rsuvdf  (m)     
    REAL(KIND=r8), INTENT(IN   ) :: rsirbm  (m)     
    REAL(KIND=r8), INTENT(IN   ) :: rsirdf  (m)     
    REAL(KIND=r8), INTENT(IN   ) :: psc     (m)     

    REAL(KIND=r8), intent(in   ) :: fice    (m,np)  
    REAL(KIND=r8), intent(in   ) :: rel     (m,np)  
    REAL(KIND=r8), intent(in   ) :: rei     (m,np)  
    REAL(KIND=r8), INTENT(IN   ) :: csscgp  (m,np  )

    ! output variables
    REAL(KIND=r8), INTENT(OUT  ) :: dscld1  (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: sc1     (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: acld1   (m,np)  
    REAL(KIND=r8), INTENT(OUT  ) :: rvbc1   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: rvdc1   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: rnbc1   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: rndc1   (m)     

    REAL(KIND=r8), INTENT(OUT  ) :: dsclr1  (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: sl1     (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: aclr1   (m,np)  
    REAL(KIND=r8), INTENT(OUT  ) :: rvbl1   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: rvdl1   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: rnbl1   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: rndl1   (m)     

    ! Parameters

    ! specify aerosol properties of continental (c) aerosols
    REAL(KIND=r8), PARAMETER, DIMENSION(11) :: tau_c = (/ &
         2.432_r8,2.1_r8,2.1_r8,1.901_r8,1.818_r8,1.76_r8,1.562_r8, &
         1.0_r8,0.568_r8, 0.281_r8,0.129_r8 /)

    REAL(KIND=r8), PARAMETER, DIMENSION(11) :: ssa_c = (/ &
         0.653_r8,0.78_r8,0.78_r8,0.858_r8,0.886_r8,0.892_r8,0.903_r8, &
         0.891_r8,0.836_r8,0.765_r8,0.701_r8 /)

    REAL(KIND=r8), PARAMETER, DIMENSION(11) :: asym_c = (/ &
         0.726_r8,0.686_r8,0.686_r8,0.666_r8,0.658_r8,0.656_r8,0.65_r8, &
         0.637_r8,0.632_r8,0.66_r8,0.776_r8 /)

    ! specify aerosol properties of maritime (m) aerosols

    REAL(KIND=r8), PARAMETER, DIMENSION(11) :: tau_m = (/ &
         1.33_r8,1.241_r8,1.241_r8,1.192_r8,1.173_r8,1.161_r8,1.117_r8, &
         1.00_r8,0.906_r8,0.799_r8,0.603_r8 /)

    REAL(KIND=r8), PARAMETER, DIMENSION(11) :: ssa_m = (/ &
         0.84_r8,0.921_r8,0.921_r8,0.962_r8,0.977_r8,0.98_r8,0.987_r8, &
         0.989_r8,0.987_r8,0.986_r8,0.861_r8 /)

    REAL(KIND=r8), PARAMETER, DIMENSION(11) :: asym_m = (/ &
         0.774_r8,0.757_r8,0.757_r8,0.746_r8,0.745_r8,0.744_r8,0.743_r8, &
         0.743_r8,0.756_r8,0.772_r8,0.802_r8 /)


    ! local  variables

    REAL(KIND=r8), DIMENSION(m,np,11) :: taual,ssaal,asyal
    REAL(KIND=r8), DIMENSION(m,np,3)  :: taucld,reff
    REAL(KIND=r8), DIMENSION(m,np)    :: fcld,dzint,aotb_c,aotb_m
    REAL(KIND=r8), DIMENSION(m,np+1)  :: flx,flc,flx_d,flx_u,flc_d,flc_u
    REAL(KIND=r8), DIMENSION(m)       :: fdiruv,fdifuv,fdirpar,fdifpar,fdirir,fdifir

    REAL(KIND=r8), DIMENSION(m)       :: fdiruv_c,fdifuv_c,fdirpar_c,fdifpar_c,fdirir_c,fdifir_c



    !hmjb new indexes for high/mid/low clouds
    INTEGER :: ict(m),icb(m)

    REAL(KIND=r8) :: topa(m),hzmask(m,np),heat
    INTEGER :: i,k,ib

    ! Initialize local vectors and output variables
    dscld1=0.0_r8;sc1=0.0_r8;rvbc1=0.0_r8;rvdc1=0.0_r8;rnbc1=0.0_r8;rndc1=0.0_r8;acld1=0.0_r8
    dsclr1=0.0_r8;sl1=0.0_r8;rvbl1=0.0_r8;rvdl1=0.0_r8;rnbl1=0.0_r8;rndl1=0.0_r8;aclr1=0.0_r8
    taual=0.0_r8;ssaal=0.0_r8;asyal=0.0_r8
    taucld=0.0_r8;reff=0.0_r8
    fcld=0.0_r8;dzint=0.0_r8;aotb_c=0.0_r8;aotb_m=0.0_r8
    flx=0.0_r8;flc=0.0_r8;flx_d=0.0_r8;flx_u=0.0_r8;flc_d=0.0_r8;flc_u=0.0_r8
    fdiruv=0.0_r8;fdifuv=0.0_r8;fdirpar=0.0_r8;fdifpar=0.0_r8;fdirir=0.0_r8;fdifir=0.0_r8
    fdiruv_c=0.0_r8;fdifuv_c=0.0_r8;fdirpar_c=0.0_r8;fdifpar_c=0.0_r8;fdirir_c=0.0_r8;fdifir_c=0.0_r8
    topa=0.0_r8;hzmask=0.0_r8;heat=0.0_r8
    ! subroutine starts here
    ict=1
    icb=1
    dscld1=0.0_r8
    sc1=0.0_r8
    rvbc1=0.0_r8
    rvdc1=0.0_r8
    rnbc1=0.0_r8
    rndc1=0.0_r8
    acld1=0.0_r8

    dsclr1=0.0_r8
    sl1=0.0_r8
    rvbl1=0.0_r8
    rvdl1=0.0_r8
    rnbl1=0.0_r8
    rndl1=0.0_r8
    aclr1=0.0_r8
    
    ! specify level indices separating high clouds from middle clouds
    ! (ict), and middle clouds from low clouds (icb).  this levels
    ! correspond to 400mb and 700 mb roughly.

    ! CPTEC-GCM works in sigma levels, hence in all columns the same
    ! layer will correspond to 0.4 and 0.7.Therefore, search is
    ! done only in the 1st column
    DO k=1,np
       DO i=1,m
          IF (pl(i,k)/psc(i).le.0.4_r8.and.pl(i,k+1)/psc(i).gt.0.4_r8) ict(i)=k
          IF (pl(i,k)/psc(i).le.0.7_r8.and.pl(i,k+1)/psc(i).gt.0.7_r8) icb(i)=k
       END DO
    ENDDO

    ! specify cloud optical thickness (taucld), amount (fcld), effective
    ! particle size (reff).cloud ice (index 1), liquid (index 2), and
    ! rain (index 3) are allowed to co-exit in a layer.
    ! cwc is the cloud ice/water concentration. if cldwater=.true.,
    ! taucld is computed from cwc and reff.  if cldwater=.false.,
    ! taucld is an input parameter, and cwc is irrelevent
    !m,np,3
    IF (icld.eq.1) THEN
       DO k=1,np
          DO i=1,m
             taucld(i,k,1)=0.0_r8
             taucld(i,k,2)=tauc(i,k)
             reff  (i,k,1) = 80.0_r8   ! ice particles
             reff  (i,k,2) = 5.25_r8   ! water particles
          END DO
       END DO
    ELSEIF (icld.eq.4 .or. icld.eq.5 .or. icld.eq.6 .or. icld.eq.7) THEN
       DO k=1,np
          DO i=1,m
             taucld(i,k,1)=tauc(i,k)*fice(i,k)
             taucld(i,k,2)=tauc(i,k)*(1.0_r8-fice(i,k))
             reff  (i,k,1) = rei(i,k)   ! ice particles
             reff  (i,k,2) = rel(i,k)   ! water particles
          END DO
       END DO
    ELSE
       WRITE(nfprt,*) 'error! icld must be 1 or 4 with Clirad-sw-m '
       STOP
    ENDIF
    DO k=1,np
       DO i=1,m
          taucld(i,k,3)  = 0.0_r8 ! no droplets
          reff  (i,k,3)  = 0.0_r8 ! no droplets
          fcld  (i,k)    = csscgp(i,k)  ! cloud field
       END DO
    END DO

    ! calculation of background aerosol optical depth profile: aotb(m,np)

    ! calculation of layer depth in km
    DO i=1,m
       DO k=1,np
          dzint(i,k)=0.0660339_r8*(log10(pl(i,k+1))-log10(pl(i,k)))*ta(i,k)
       ENDDO
    ENDDO

    ! calculation of number of layers in 2 km depth (nta)
    !hmjb Now we save, for each (i,k), a number hzmask.
    !     This number is 1, if bottom of layer below 2km
    !     This number is 0, if bottom of layer above 2km
    !  I did this because it is faster to multiply by this
    !  matrix than to do a loop in (i) with two loops in
    !  k inside (1:2km and 2km:ground).

    hzmask=0.0_r8
    DO i=1,m
       k=np+1
       topa(i)=0.0_r8
       DO while (topa(i).lt.2.0_r8)
          k=k-1
          topa(i)=topa(i)+dzint(i,k)
       ENDDO
       hzmask(i,k:np)=1.0_r8
    ENDDO

    ! background aerosol profile with optical depth
    ! extinction coefficient 0.1 or 0.5 km-1 in each layer

    !hmjb The total column aerosol (0.22 or 0.14) should be distributed
    !     over the first 2km. However, we see from the calculation above
    !     that we distribute the aersol from the first layer above 2km
    !     down to the ground. Therefore, we must consider that probably
    !     the height of this region will be more than 2km.
    ! I changed the loop above so that we keep track of the total
    !     height, in each column, of the levels where we will spread the
    !     aerosol. Now the distribution of the aerosol inside each layer.
    !     More than that, if we sum the total aerosol in the full column,
    !     it will add to the stipulated value.
    DO k=1,np
       DO i=1,m
          aotb_c(i,k)=asolc*dzint(i,k)*hzmask(i,k)/2.0_r8
          aotb_m(i,k)=asolm*dzint(i,k)*hzmask(i,k)/2.0_r8
          !need more testing
          !new   aotb_c(i,:)=asolc*dzint(i,:)*hzmask(i,:)/topa(i)
          !new   aotb_m(i,:)=asolm*dzint(i,:)*hzmask(i,:)/topa(i)
       END DO
    END DO

    ! specify aerosol optical thickness (taual), single-scattering
    ! albedo (ssaal), and asymmetry factor (asyal)
    ! nta is top level of aerosol layer over the ground

    DO k=1,np

       DO ib=1,11
          DO i=1,m
             IF(schemes == 3)THEN
                IF (dmask(i).gt.0.and.dmask(i).le.14) THEN
                   ! sibmask=1..12 is land with different vegetation types
                   taual(i,k,ib) = tau_c(ib)*aotb_c(i,k)*hzmask(i,k)
                   ssaal(i,k,ib) = ssa_c(ib)
                   asyal(i,k,ib) = asym_c(ib)
                ELSE
                   ! sibmask -1 or 0 means ice/water and 13 is permanent ice (greenland and antartic)
                   taual(i,k,ib) = tau_m(ib)*aotb_m(i,k)*hzmask(i,k)
                   ssaal(i,k,ib) = ssa_m(ib)
                   asyal(i,k,ib) = asym_m(ib)
                END IF
             ELSE
                IF (dmask(i).gt.0.and.dmask(i).le.12) THEN
                   ! sibmask=1..12 is land with different vegetation types
                   taual(i,k,ib) = tau_c(ib)*aotb_c(i,k)*hzmask(i,k)
                   ssaal(i,k,ib) = ssa_c(ib)
                   asyal(i,k,ib) = asym_c(ib)
                ELSE
                   ! sibmask -1 or 0 means ice/water and 13 is permanent ice (greenland and antartic)
                   taual(i,k,ib) = tau_m(ib)*aotb_m(i,k)*hzmask(i,k)
                   ssaal(i,k,ib) = ssa_m(ib)
                   asyal(i,k,ib) = asym_m(ib)
                END IF
             END IF
          END DO
       END DO
    END DO

    ! compute solar fluxes

    CALL soradcld ( &
         m                       ,&! INTEGER      , INTENT(IN   ) :: m
         np                      ,&! INTEGER      , INTENT(IN   ) :: np
         pl       (1:m,1:np+1 )  ,&! REAL(KIND=r8), INTENT(IN   ) :: pl       (1:m,1:np+1 )
         ta       (1:m,1:np   )  ,&! REAL(KIND=r8), INTENT(IN   ) :: ta       (1:m,1:np   )
         wa       (1:m,1:np   )  ,&! REAL(KIND=r8), INTENT(IN   ) :: wa       (1:m,1:np   )
         oa       (1:m,1:np   )  ,&! REAL(KIND=r8), INTENT(IN   ) :: oa       (1:m,1:np   )
         rco2     (1:m,1:np)     ,&! REAL(KIND=r8), INTENT(IN   ) :: co2      (1:m,1:np)
         taucld   (1:m,1:np,1:3) ,&! REAL(KIND=r8), INTENT(IN   ) :: taucld   (1:m,1:np,1:3)
         reff     (1:m,1:np,1:3 ),&! REAL(KIND=r8), INTENT(IN   ) :: reff     (1:m,1:np,1:3 )
         fcld     (1:m,1:np     ),&! REAL(KIND=r8), INTENT(IN   ) :: fcld     (1:m,1:np     )
         ict      (1:m)          ,&! INTEGER      , INTENT(IN   ) :: ict      (1:m)
         icb      (1:m)          ,&! INTEGER      , INTENT(IN   ) :: icb      (1:m)
         taual    (1:m,1:np,1:11),&! REAL(KIND=r8), INTENT(IN   ) :: taual    (1:m,1:np,1:11)
         ssaal    (1:m,1:np,1:11),&! REAL(KIND=r8), INTENT(IN   ) :: ssaal    (1:m,1:np,1:11)
         asyal    (1:m,1:np,1:11),&! REAL(KIND=r8), INTENT(IN   ) :: asyal    (1:m,1:np,1:11)
         cosz     (1:m)          ,&! REAL(KIND=r8), INTENT(IN   ) :: cosz     (1:m)    
         rsuvbm   (1:m)          ,&! REAL(KIND=r8), INTENT(IN   ) :: rsuvbm   (1:m)    
         rsuvdf   (1:m)          ,&! REAL(KIND=r8), INTENT(IN   ) :: rsuvdf   (1:m)    
         rsirbm   (1:m)          ,&! REAL(KIND=r8), INTENT(IN   ) :: rsirbm   (1:m)    
         rsirdf   (1:m)          ,&! REAL(KIND=r8), INTENT(IN   ) :: rsirdf   (1:m)    
         flx      (1:m,1:np+1)   ,&! REAL(KIND=r8), INTENT(INOUT) :: flx      (1:m,1:np+1)
         flc      (1:m,1:np+1)   ,&! REAL(KIND=r8), INTENT(INOUT) :: flc      (1:m,1:np+1)
         fdiruv   (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdiruv   (1:m)    
         fdifuv   (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifuv   (1:m)    
         fdirpar  (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirpar  (1:m)    
         fdifpar  (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifpar  (1:m)    
         fdirir   (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirir   (1:m)    
         fdifir   (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifir   (1:m)    
         fdiruv_c (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdiruv_c (1:m)    
         fdifuv_c (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifuv_c (1:m)    
         fdirpar_c(1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirpar_c(1:m)    
         fdifpar_c(1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifpar_c(1:m)    
         fdirir_c (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirir_c (1:m)    
         fdifir_c (1:m)          ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifir_c (1:m)    
         flx_d    (1:m,1:np+1)   ,&! REAL(KIND=r8), INTENT(INOUT) :: flx_d    (1:m,1:np+1)
         flx_u    (1:m,1:np+1)   ,&! REAL(KIND=r8), INTENT(INOUT) :: flx_u    (1:m,1:np+1)
         flc_d    (1:m,1:np+1)   ,&! REAL(KIND=r8), INTENT(INOUT) :: flc_d    (1:m,1:np+1)
         flc_u    (1:m,1:np+1)    )! REAL(KIND=r8), INTENT(INOUT) :: flc_u    (1:m,1:np+1)

    ! convert the units of flx and flc from fraction to w/m^2
    ! transfer to the global model fluxes

    DO i=1,m
       dscld1(i)=flx(i,1)*s0*cosz(i)
       sc1   (i)=flx(i,np+1)*s0*cosz(i)
       rvbc1 (i)=(fdiruv(i)+fdirpar(i))*s0*cosz(i)
       rvdc1 (i)=(fdifuv(i)+fdifpar(i))*s0*cosz(i)
       rnbc1 (i)=fdirir(i)*s0*cosz(i)
       rndc1 (i)=fdifir(i)*s0*cosz(i)
       !
       dsclr1(i)=flc(i,1)*s0*cosz(i)
       sl1   (i)=flc(i,np+1)*s0*cosz(i)
       rvbl1 (i)=(fdiruv_c(i)+fdirpar_c(i))*s0*cosz(i)
       rvdl1 (i)=(fdifuv_c(i)+fdifpar_c(i))*s0*cosz(i)
       rnbl1 (i)=fdirir_c(i)*s0*cosz(i)
       rndl1 (i)=fdifir_c(i)*s0*cosz(i)
    ENDDO


    ! compute heating rates, c/day

    DO k=1,np
       DO i=1,m
          heat=8.4410_r8*s0*cosz(i)/(pl(i,k+1)-pl(i,k))
          aclr1(i,k)=(flc(i,k)-flc(i,k+1))*heat
          acld1(i,k)=(flx(i,k)-flx(i,k+1))*heat
       ENDDO
    ENDDO

  END SUBROUTINE cloudy


  ! ==============
  !
  !  clirad-sw
  !
  ! ==============
  !
  ! Subroutine: Soradcld
  !
  ! $Author: pkubota $
  ! Modifications: T. Tarasova, 2005
  ! Modifications: H. Barbosa, 2005
  !
  ! Description:
  !
  !following  the nasa  technical  memorandum (nasa/tm-1999-104606,vol.
  !15) of chou and suarez (1999), this routine computes solarfluxes due
  !to absorption by water vapor, ozone, co2, o2, clouds,andaerosols and
  !due to scattering by clouds, aerosols, and gases.
  !
  !this code computes fluxes simultaneously for m soundings.
  !
  !cloud ice, liquid,  and rain particles are allowed  to co-exist in a
  !layer.
  !
  !there is an option of  providing either cloud ice/water mixing ratio
  !(cwc) or optical thickness (taucld).  if the former is provided, set
  !cldwater=.true.,  and taucld  is computed  from  cwc and  reff as  a
  !function  of  spectra band.   otherwise,  set cldwater=.false.,  and
  !specify taucld, independent of spectral band.
  !
  !if  no information  is available  for the  effective  particle size,
  !reff, default values of 10 micron for liquid water and 75 micron for
  !ice may be  used.  the size of raindrops,  reff(3), is irrelevant in
  !this code. it can be set to any values.  for a clear layer, reff can
  !be set to any values except zero.
  !
  !the  maximum-random   assumption  is  appliedfor  treating  cloud
  !overlapping. clouds  are grouped into  high, middle, and  low clouds
  !separated  by  the level  indices  ict  and  icb.  for  detail,  see
  !SUBROUTINE "cldscale".
  !
  !in   a  high   spatial-resolution   atmospheric  model,   fractional
  !cloudcover might be  computed to be either 0 or 1.   in such a case,
  !scaling  of the  cloud optical  thickness isnot  necessary,  and the
  !computation   can  bemade  faster   by   setting  overcast=.true.
  !otherwise, set the option overcast=.false.
  !
  !aerosol optical  thickness, single-scattering albedo,  and asymmetry
  !factor can be specified as functions of height and spectral band.
  !
  ! Inputs: units        size
  !
  !      m: number of soundings n/d      1
  !     np: number of atmospheric layers n/d      1
  !     pl: level pressure  mb      m*(np+1)
  !     ta: layer temperature k      m*np
  !     wa: layer specific humidity gm/gm        m*np
  !     oa: layer ozone concentration gm/gm        m*np
  !    co2: co2 mixing ratio by volume pppv      1
  !  overcast: option for scaling cloud optical thickness n/d      1
  !   "true"  = scaling is not required
  !   "fasle" = scaling is required
  !  cldwater: input option for cloud optical thickness n/d      1
  !   "true"  = taucld is provided
  !   "false" = cwp is provided
  !    cwc: cloud water mixing ratio gm/gm        m*np*3
  !   index 1 for ice particles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  ! taucld: cloud optical thickness n/d      m*np*3
  !   index 1 for ice particles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  !   reff: effective cloud-particle size   micrometer   m*np*3
  !   index 1 for ice paticles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  !   fcld: cloud amount fraction     m*np
  !    ict: level index separating high and middle clouds   n/d      m
  !    icb: level indiex separating middle and low clouds   n/d      m
  !  taual: aerosol optical thickness n/d      m*np*11
  !  ssaal: aerosol single-scattering albedo n/d      m*np*11
  !  asyal: aerosol asymmetry factor n/d      m*np*11
  ! in the uv region :
  !    index  1 for the 0.175-0.225 micron band
  !    index  2 for the 0.225-0.245; 0.260-0.280 micron band
  !    index  3 for the 0.245-0.260 micron band
  !    index  4 for the 0.280-0.295 micron band
  !    index  5 for the 0.295-0.310 micron band
  !    index  6 for the 0.310-0.320 micron band
  !    index  7 for the 0.325-0.400 micron band
  ! in the par region :
  !    index  8 for the 0.400-0.700 micron band
  ! in the infrared region :
  !    index  9 for the 0.700-1.220 micron band
  !    index 10 for the 1.220-2.270 micron band
  !    index 11 for the 2.270-10.00 micron band
  !   cosz: cosine of solar zenith angle       n/d    m
  ! rsuvbm: uv+vis sfc albedo for beam rad for wl<0.7 micron      fraction     m
  ! rsuvdf: uv+vis sfc albedo for diffuse rad  for wl<0.7 micron  fraction     m
  ! rsirbm: ir sfc albedo for beam rad for wl>0.7 micron       fraction     m
  ! rsirdf: ir sfc albedo for diffuse rad         fraction     m
  !
  ! Outputs: (updated parameters)
  !
  !    flx: all-sky   net downward flux       fraction     m*(np+1)
  !    flc: clear-sky net downward flux       fraction     m*(np+1)
  ! fdiruv: all-sky direct  downward uv (.175-.4 micron) flux sfc fraction     m
  ! fdifuv: all-sky diffuse downward uv flux at the surface       fraction     m
  !fdirpar: all-sky direct  downward par (.4-.7 micron) flux sfc  fraction     m
  !fdifpar: all-sky diffuse downward par flux at the surface      fraction     m
  ! fdirir: all-sky direct  downward ir (.7-10 micron) flux sfc   fraction     m
  ! fdifir: all-sky diffuse downward ir flux at the surface ()    fraction     m
  !
  !
  !
  ! NOTES
  !
  ! (1) the unit of output fluxes (flx,flc,etc.) is fraction of the
  !     insolation at the top of the atmosphere.  therefore, fluxes
  !     are the output fluxes multiplied by the extra-terrestrial solar
  !     flux and the cosine of the solar zenith angle.
  ! (2) pl( ,1) is the pressure at the top of the model, and
  !     pl( ,np+1) is the surface pressure.
  ! (3) the pressure levels ict and icb correspond approximately
  !     to 400 and 700 mb.
  !
  !  if coding errors are found, please notify ming-dah chou at
  !  chou@climate.gsfc.nasa.gov
  !
  !
  !
  !
  SUBROUTINE Soradcld (&
       m         ,&! INTEGER      , INTENT(IN   ) :: m
       np        ,&! INTEGER      , INTENT(IN   ) :: np
       pl        ,&! REAL(KIND=r8), INTENT(IN   ) :: pl     (m,np+1 )
       ta        ,&! REAL(KIND=r8), INTENT(IN   ) :: ta     (m,np   )
       wa        ,&! REAL(KIND=r8), INTENT(IN   ) :: wa     (m,np   )
       oa        ,&! REAL(KIND=r8), INTENT(IN   ) :: oa     (m,np   )
       co2       ,&! REAL(KIND=r8), INTENT(IN   ) :: co2      (m,np)
       taucld    ,&! REAL(KIND=r8), INTENT(IN   ) :: taucld   (m,np,3 )
       reff      ,&! REAL(KIND=r8), INTENT(IN   ) :: reff     (m,np,3 )
       fcld      ,&! REAL(KIND=r8), INTENT(IN   ) :: fcld     (m,np   )
       ict       ,&! INTEGER      , INTENT(IN   ) :: ict      (m)
       icb       ,&! INTEGER      , INTENT(IN   ) :: icb      (m)
       taual     ,&! REAL(KIND=r8), INTENT(IN   ) :: taual    (m,np,11)
       ssaal     ,&! REAL(KIND=r8), INTENT(IN   ) :: ssaal    (m,np,11)
       asyal     ,&! REAL(KIND=r8), INTENT(IN   ) :: asyal    (m,np,11)
       cosz      ,&! REAL(KIND=r8), INTENT(IN   ) :: cosz     (m)    
       rsuvbm    ,&! REAL(KIND=r8), INTENT(IN   ) :: rsuvbm   (m)    
       rsuvdf    ,&! REAL(KIND=r8), INTENT(IN   ) :: rsuvdf   (m)    
       rsirbm    ,&! REAL(KIND=r8), INTENT(IN   ) :: rsirbm   (m)    
       rsirdf    ,&! REAL(KIND=r8), INTENT(IN   ) :: rsirdf   (m)    
       flx       ,&! REAL(KIND=r8), INTENT(INOUT) :: flx      (m,np+1)
       flc       ,&! REAL(KIND=r8), INTENT(INOUT) :: flc      (m,np+1)
       fdiruv    ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdiruv   (m)    
       fdifuv    ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifuv   (m)    
       fdirpar   ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirpar  (m)    
       fdifpar   ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifpar  (m)    
       fdirir    ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirir   (m)    
       fdifir    ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifir   (m)    
       fdiruv_c  ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdiruv_c (m)    
       fdifuv_c  ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifuv_c (m)    
       fdirpar_c ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirpar_c(m)    
       fdifpar_c ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifpar_c(m)    
       fdirir_c  ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdirir_c (m)    
       fdifir_c  ,&! REAL(KIND=r8), INTENT(OUT  ) :: fdifir_c (m)    
       flx_d     ,&! REAL(KIND=r8), INTENT(INOUT) :: flx_d    (m,np+1)
       flx_u     ,&! REAL(KIND=r8), INTENT(INOUT) :: flx_u    (m,np+1)
       flc_d     ,&! REAL(KIND=r8), INTENT(INOUT) :: flc_d    (m,np+1)
       flc_u      )! REAL(KIND=r8), INTENT(INOUT) :: flc_u    (m,np+1)

    IMPLICIT NONE
    INTEGER i


    ! parameters for co2 transmission tables
    !INTEGER, PARAMETER :: nu=43
    !INTEGER, PARAMETER :: nw=37
    !INTEGER, PARAMETER :: nx=62
    !INTEGER, PARAMETER :: ny=101

    ! cah is the co2 absorptance in band 10
    !REAL(KIND=r8), DIMENSION(nx,ny) :: coa
    !INCLUDE "coa.data90"

    ! coa is the co2 absorptance in strong absorption regions of band 11
    !REAL(KIND=r8), DIMENSION(nu,nw) :: cah
    !INCLUDE "cah.data90"

    ! input variables
    INTEGER      , INTENT(IN   ) :: m
    INTEGER      , INTENT(IN   ) :: np
    INTEGER      , INTENT(IN   ) :: ict    (m)
    INTEGER      , INTENT(IN   ) :: icb    (m)
    REAL(KIND=r8), INTENT(IN   ) :: taual  (m,np,11)
    REAL(KIND=r8), INTENT(IN   ) :: ssaal  (m,np,11)
    REAL(KIND=r8), INTENT(IN   ) :: asyal  (m,np,11)
    REAL(KIND=r8), INTENT(IN   ) :: taucld (m,np,3 )
    REAL(KIND=r8), INTENT(IN   ) :: reff   (m,np,3 )
    REAL(KIND=r8), INTENT(IN   ) :: pl     (m,np+1 )
    REAL(KIND=r8), INTENT(IN   ) :: ta     (m,np   )
    REAL(KIND=r8), INTENT(IN   ) :: wa     (m,np   )
    REAL(KIND=r8), INTENT(IN   ) :: oa     (m,np   )
    REAL(KIND=r8), INTENT(IN   ) :: fcld   (m,np   )
    REAL(KIND=r8), INTENT(IN   ) :: cosz   (m)    
    REAL(KIND=r8), INTENT(IN   ) :: rsuvbm (m)    
    REAL(KIND=r8), INTENT(IN   ) :: rsuvdf (m)    
    REAL(KIND=r8), INTENT(IN   ) :: rsirbm (m)    
    REAL(KIND=r8), INTENT(IN   ) :: rsirdf (m)    
    REAL(KIND=r8), INTENT(IN   ) :: co2    (m,np)

    ! output variables
    REAL(KIND=r8), INTENT(INOUT) :: flx      (m,np+1)
    REAL(KIND=r8), INTENT(INOUT) :: flc      (m,np+1)
    REAL(KIND=r8), INTENT(INOUT) :: flx_d    (m,np+1)
    REAL(KIND=r8), INTENT(INOUT) :: flx_u    (m,np+1)
    REAL(KIND=r8), INTENT(INOUT) :: flc_d    (m,np+1)
    REAL(KIND=r8), INTENT(INOUT) :: flc_u    (m,np+1)
    REAL(KIND=r8), INTENT(OUT  ) :: fdiruv   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdifuv   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdirpar  (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdifpar  (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdirir   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdifir   (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdiruv_c (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdifuv_c (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdirpar_c(m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdifpar_c(m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdirir_c (m)     
    REAL(KIND=r8), INTENT(OUT  ) :: fdifir_c (m)     

    ! local  variables
    INTEGER :: k,ntop
    INTEGER, DIMENSION(m) :: nctop

    REAL(KIND=r8)                    :: x, w1,dw,u1,du
    REAL(KIND=r8), DIMENSION(m)      :: snt!, cnt
    REAL(KIND=r8), DIMENSION(m,np)   :: dp, wh, oh, scal,cnt
    REAL(KIND=r8), DIMENSION(m,np+1) :: swu, swh, so2, df

    ! initialize local vectors and output variables
    nctop=0
    x=0.0_r8; w1=0.0_r8;dw=0.0_r8;u1=0.0_r8;du=0.0_r8
    snt=0.0_r8; cnt=0.0_r8
    dp=0.0_r8; wh=0.0_r8; oh=0.0_r8; scal=0.0_r8
    swu=0.0_r8; swh=0.0_r8; so2=0.0_r8; df=0.0_r8
   
    !flx=0.0_r8;flc=0.0_r8;flx_d=0.0_r8;flx_u=0.0_r8; flc_d=0.0_r8;flc_u=0.0_r8
    fdiruv=0.0_r8;fdifuv=0.0_r8;fdirpar=0.0_r8;fdifpar=0.0_r8;fdirir=0.0_r8;fdifir=0.0_r8
    fdiruv_c=0.0_r8;fdifuv_c=0.0_r8;fdirpar_c=0.0_r8;fdifpar_c=0.0_r8;fdirir_c=0.0_r8;fdifir_c=0.0_r8
    ! subroutine starts here
    DO i=1,m
       swh(i,1)=0.0_r8
       so2(i,1)=0.0_r8
       snt(i)=1.0_r8/cosz(i) ! snt is the secant of the solar zenith angle
    ENDDO

    DO k=1,np
       DO i=1,m

          ! compute layer thickness. indices for the surface level and
          ! surface layer are np+1 and np, respectively.
          dp(i,k)=pl(i,k+1)-pl(i,k)
          !
          ! compute scaled water vapor amount following eqs. (3.3) and (3.5)
          ! unit is g/cm**2
          !
          scal(i,k)=dp(i,k)*(0.5_r8*(pl(i,k)+pl(i,k+1))/300.0_r8)**0.8_r8
          wh(i,k)=1.02_r8*wa(i,k)*scal(i,k)*(1.0_r8+0.00135_r8*(ta(i,k)-240.0_r8))+1.e-11_r8
          swh(i,k+1)=swh(i,k)+wh(i,k)
          !
          ! compute ozone amount, unit is (cm-atm)stp
          ! the number 466.7 is the unit conversion factor
          ! from g/cm**2 to (cm-atm)stp
          !
          oh(i,k)=1.02_r8*oa(i,k)*dp(i,k)*466.7_r8 +1.0e-11_r8
          !
          ! compute layer cloud water amount (gm/m**2)
          ! the index is 1 for ice crystals, 2 for liquid drops, and
          ! 3 for rain drops
          !
          !     x=1.02*10000.*dp(i,k)
          !     cwp(i,k,1)=x*cwc(i,k,1)
          !     cwp(i,k,2)=x*cwc(i,k,2)
          !     cwp(i,k,3)=x*cwc(i,k,3)
       ENDDO
    ENDDO

    ! initialize fluxes for all-sky (flx), clear-sky (flc), and
    ! flux reduction (df)

    DO k=1,np+1
       DO i=1,m
          flx   (i,k)=0.0_r8
          flx_d (i,k)=0.0_r8    ! new
          flx_u (i,k)=0.0_r8    ! new
          flc   (i,k)=0.0_r8
          flc_d (i,k)=0.0_r8    ! new
          flc_u (i,k)=0.0_r8    ! new
          df    (i,k)=0.0_r8
       ENDDO
    ENDDO

    ! compute solar uv and par fluxes

    CALL soluvcld (&
         m                        , &!    INTEGER, INTENT(IN) :: m
         np                       , &!    INTEGER, INTENT(IN) :: np
         wh        (1:m,1:np)     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: wh
         oh        (1:m,1:np)     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: oh
         dp        (1:m,1:np)     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: dp
         taucld    (1:m,1:np,1:3) , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: taucld
         reff      (1:m,1:np,1:3) , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: reff
         ict       (1:m)          , &!    INTEGER, INTENT(IN) :: ict(m)
         icb       (1:m)          , &!    INTEGER, INTENT(IN) :: icb(m)
         fcld      (1:m,1:np)     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np):: fcld
         cosz      (1:m)          , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m):: cosz
         taual     (1:m,1:np,1:11), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: taual
         ssaal     (1:m,1:np,1:11), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: ssaal
         asyal     (1:m,1:np,1:11), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: asyal
         rsuvbm    (1:m)          , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: rsuvbm
         rsuvdf    (1:m)          , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: rsuvdf
         flx       (1:m,1:np+1)   , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx
         flc       (1:m,1:np+1)   , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc
         fdiruv    (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdiruv
         fdifuv    (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifuv
         fdirpar   (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdirpar
         fdifpar   (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifpar
         fdiruv_c  (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdiruv_c
         fdifuv_c  (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifuv_c
         fdirpar_c (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdirpar_c
         fdifpar_c (1:m)          , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifpar_c
         flx_d     (1:m,1:np+1)   , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_d
         flx_u     (1:m,1:np+1)   , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_u
         flc_d     (1:m,1:np+1), &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc_d
         flc_u     (1:m,1:np+1)  )!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc_u


    ! compute and update solar ir fluxes
    CALL solircld (m,np,wh,dp,  &
         taucld,reff,ict,icb,fcld,cosz,  &
         taual,ssaal,asyal,rsirbm,rsirdf,  &
         flx,flc,fdirir,fdifir,fdirir_c,fdifir_c,  &
         flx_d,flx_u,flc_d,flc_u)  ! new


    ! compute pressure-scaled o2 amount following eq. (3.5) with
    !     f=1. unit is (cm-atm)stp.
    !     the constant 165.22 equals (1000/980)*23.14%*(22400/32)

    DO k=1,np
       DO i=1,m
         cnt(i,k)=165.22_r8*snt(i)
       ENDDO
    END DO

    DO k=1,np
       DO i=1,m
          so2(i,k+1)=so2(i,k)+scal(i,k)*cnt(i,k)
       ENDDO
    ENDDO

    ! compute flux reduction due to oxygen following eq. (3.18)
    !     the constant 0.0633 is the fraction of insolation contained
    !     in the oxygen bands

    DO k=2,np+1
       DO i=1,m
          x=so2(i,k)
          df(i,k)=0.0633_r8*(1.0_r8-exp(-0.000145_r8*sqrt(x)))
       ENDDO
    ENDDO

    ! for solar heating due to co2
    DO k=1,np
       DO i=1,m
          cnt(i,k)=co2(i,k)*snt(i)
       END DO
    END DO

    ! scale co2 amounts following eq. (3.5) with f=1.
    !     unit is (cm-atm)stp.
    !     the constant 789 equals (1000/980)*(44/28.97)*(22400/44)

    DO k=1,np
       DO i=1,m
          x=789.0_r8*cnt(i,k)
          so2(i,k+1)=so2(i,k)+x*scal(i,k)+1.0e-11_r8
       ENDDO
    ENDDO
    ! for co2 absorption in band 10 where absorption due to
    !     water vapor and co2 are both moderate

    u1=-3.0_r8
    du=0.15_r8
    w1=-4.0_r8
    dw=0.15_r8

    ! so2 and swh are the co2 and water vapor amounts integrated
    !     from the top of the atmosphere

    DO k=2,np+1
       DO i=1,m
          swu(i,k)=log10(so2(i,k))
          swh(i,k)=log10(swh(i,k)*snt(i))
       ENDDO
    ENDDO

    ! df is the updated flux reduction given by the second term on the
    !     right-hand-side of eq. (3.24) divided by so

    CALL rflx(m,np,swu,u1,du,nu,swh,w1,dw,nw,cah,df)

    ! for co2 absorption in band 11 where the co2 absorption has
    !     a large impact on the heating of middle atmosphere.

    u1=0.000250_r8
    du=0.000050_r8
    w1=-2.0_r8
    dw=0.05_r8

    DO k=2,np+1
       DO i=1,m
          swu(i,k)=co2(i,k-1)*snt(i)
       ENDDO
    END DO

    ! co2 mixing ratio is independent of space

!    DO k=2,np+1
       DO i=1,m
          swu(i,1)=swu(i,2)
       ENDDO
!    ENDDO

    ! swh is the logarithm of pressure

    DO k=2,np+1
       DO i=1,m
          swh(i,k)=log10(pl(i,k))
       ENDDO
    ENDDO

    ! df is the updated flux reduction derived from the table given by
    !     eq. (3.19)

    CALL rflx(m,np,swu,u1,du,nx,swh,w1,dw,ny,coa,df)

    ! adjustment for the effect of o2 and co2 on clear-sky fluxes.
    !     both flc and df are positive quantities

    DO k=1,np+1
       DO i=1,m
          flc(i,k)=flc(i,k)-df(i,k)
          flc_d(i,k)=flc_d(i,k)-df(i,k)   ! new
       ENDDO
    ENDDO

    ! adjustment for the direct downward flux (CLEAR)

    DO i=1,m
       fdirir_c(i)=fdirir_c(i)-df(i,np+1)
       IF (fdirir_c(i) .lt. 0.0_r8) fdirir_c(i)=0.0_r8
    ENDDO

    ! identify top cloud-layer

    DO i=1,m
       nctop(i)=np+1
    ENDDO

    DO k=1,np
       DO i=1,m
          IF (fcld(i,k).gt.0.01_r8 .and. nctop(i).eq.np+1) THEN
             nctop(i)=k
          ENDIF
       ENDDO
    ENDDO

    DO i=1,m
       !hmjb ERROR
       !     ntop=nctop(m)
       ntop=nctop(i)
       !
       ! adjust fluxes above clouds following eq. (6.17)
       DO k=1,ntop
          flx(i,k)=flx(i,k)-df(i,k)
          flx_d(i,k)=flx_d(i,k)-df(i,k)   ! new
       ENDDO
       ! adjust fluxes below cloud top following eq. (6.18)
       IF (ntop.lt.np+1) THEN
          DO k=ntop+1,np+1
             df(i,k)=df(i,k)*(flx(i,k)/flc(i,k))
             flx(i,k)=flx(i,k)-df(i,k)
             flx_d(i,k)=flx_d(i,k)-df(i,k)! new
          ENDDO
       ENDIF
    ENDDO
    ! adjustment for the direct downward flux
    DO i=1,m
       fdirir(i)=fdirir(i)-df(i,np+1)
       IF (fdirir(i) .lt. 0.0_r8) fdirir(i)=0.0_r8
    ENDDO

  END SUBROUTINE Soradcld

  !
  ! Subroutine: SOLIRCLD
  !
  ! $Author: pkubota $
  ! Modifications: T. Tarasova, 2005
  ! Modifications: H. Barbosa, 2005
  !
  ! Description:
  !  compute solar flux in the infrared region. the spectrum is divided
  !into three bands:
  !
  !       band   wavenumber(/cm)  wavelength (micron)
  !       1( 9)14280-8200   0.70-1.22
  !       2(10) 8200-4400   1.22-2.27
  !       3(11) 4400-1000   2.27-10.0
  !
  ! Inputs: units        size
  !
  !      m: number of soundings n/d      1
  !     np: number of atmospheric layers n/d      1
  !     wh: layer scaled-water vapor content gm/cm^2      m*np
  !  overcast: option for scaling cloud optical thickness n/d      1
  !   "true"  = scaling is not required
  !   "fasle" = scaling is required
  !  cldwater: input option for cloud optical thickness n/d      1
  !   "true"  = taucld is provided
  !   "false" = cwp is provided
  !    cwp: cloud water amount gm/m**2      m*np*3
  !   index 1 for ice particles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  ! taucld: cloud optical thickness n/d      m*np*3
  !   index 1 for ice particles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  !   reff: effective cloud-particle size   micrometer   m*np*3
  !   index 1 for ice paticles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  !    ict: level index separating high and middle clouds   n/d      m
  !    icb: level indiex separating middle and low clouds   n/d      m
  !   fcld: cloud amount fraction     m*np
  !  taual: aerosol optical thickness n/d      m*np*11
  !  ssaal: aerosol single-scattering albedo n/d      m*np*11
  !  asyal: aerosol asymmetry factor n/d      m*np*11
  ! rsirbm: near ir surface albedo for beam radiation fraction     m
  ! rsirdf: near ir surface albedo for diffuse radiation fraction     m
  !
  ! Outputs: (updated parameters)
  !
  !    flx: all-sky   net downward flux    fraction  m*(np+1)
  !    flc: clear-sky net downward flux    fraction  m*(np+1)
  ! fdirir: all-sky direct  downward ir flux at the surface    fraction  m
  ! fdifir: all-sky diffuse downward ir flux at the surface    fraction  m
  !
  ! Local Variables
  !
  !  tauclb: scaled cloud optical thickness for beam radiation    n/d   m*np
  !  tauclf: scaled cloud optical thickness for diffuse radiation n/d   m*np
  !
  !
  !
  SUBROUTINE Solircld (m,np,wh,dp,  &
       taucld,reff,ict,icb,fcld,cosz,  &
       taual,ssaal,asyal,  &
       rsirbm,rsirdf,flx,flc,fdirir,fdifir,fdirir_c,fdifir_c,  &
       flx_d,flx_u,flc_d,flc_u)   !  new

    IMPLICIT NONE

    ! parameters
    INTEGER, PARAMETER :: nk=10
    INTEGER, PARAMETER :: nband=3

    ! water vapor absorption coefficient for 10 k-intervals. unit: cm^2/gm (table 2)
    REAL(KIND=r8), PARAMETER, DIMENSION(nk) :: xk2 = (/    &
         0.0010_r8, 0.0133_r8, 0.0422_r8, 0.1334_r8, 0.4217_r8, &
         1.334_r8,  5.623_r8,  31.62_r8,  177.8_r8,  1000.0_r8 /)

    ! water vapor k-distribution function,
    ! the sum of hk is 0.52926. unit: fraction (table 2)
    ! --- new coefficients (tarasova and fomin, 2000)
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,nk) :: hk2 = RESHAPE( &
         SOURCE = (/ &
         0.19310_r8, 0.06924_r8, 0.00310_r8, 0.05716_r8, 0.01960_r8, 0.00637_r8, &
         0.02088_r8, 0.00795_r8, 0.00526_r8, 0.02407_r8, 0.01716_r8, 0.00641_r8, &
         0.01403_r8, 0.01118_r8, 0.00542_r8, 0.00582_r8, 0.01377_r8, 0.00312_r8, &
         0.00246_r8, 0.02008_r8, 0.00368_r8, 0.00163_r8, 0.00265_r8, 0.00346_r8, &
         0.00101_r8, 0.00282_r8, 0.00555_r8, 0.00041_r8, 0.00092_r8, 0.00098_r8 /) , &
         SHAPE = (/nband , nk/) )

    ! ry is the extinction coefficient for rayleigh scattering. unit: /mb (table 3)
    REAL(KIND=r8), PARAMETER, DIMENSION(nband) :: ry2 = (/ 0.0000156_r8, 0.0000018_r8, 0.000000_r8 /)

    ! coefficients for computing the extinction coefficients of
    ! ice, water, and rain particles (table 4)
    !REAL(KIND=r8), PARAMETER, DIMENSION(nband,2) :: aib = RESHAPE( &
    !     SHAPE = (/ nband, 2 /), SOURCE = (/ &
    !     0.000333_r8, 0.000333_r8, 0.000333_r8, 2.52_r8,    2.52_r8,2.52_r8 /) )
    !REAL(KIND=r8), PARAMETER, DIMENSION(nband,2) :: awb = RESHAPE( &
    !     SHAPE = (/ nband, 2 /), SOURCE = (/ &
    !     -0.0101_r8, -0.0166_r8, -0.0339_r8, 1.72_r8,    1.85_r8,2.16_r8 /) )
    !REAL(KIND=r8), PARAMETER, DIMENSION(nband,2) :: arb = RESHAPE( &
    !     SHAPE = (/ nband, 2 /), SOURCE = (/ &
    !     0.00307_r8, 0.00307_r8, 0.00307_r8, 0.0_r8    , 0.0_r8    , 0.0_r8  /) )

    ! coefficients for computing the single-scattering co-albedo of
    !     ice, water, and rain particles (table 5)
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,3) :: aia  = RESHAPE( &
         SHAPE = (/ nband, 3 /), SOURCE = (/ &
         -0.00000260_r8,  0.00215346_r8,  0.08938331_r8, &
         0.00000746_r8,  0.00073709_r8,  0.00299387_r8, &
         0.00000000_r8, -0.00000134_r8, -0.00001038_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,3) :: awa = RESHAPE( &
         SHAPE = (/ nband, 3 /), SOURCE = (/ &
         0.00000007_r8,-0.00019934_r8, 0.01209318_r8, &
         0.00000845_r8, 0.00088757_r8, 0.01784739_r8, &
         -0.00000004_r8,-0.00000650_r8,-0.00036910_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,3) :: ara = RESHAPE( &
         SHAPE = (/ nband, 3 /), SOURCE = (/ &
         0.029_r8,  0.342_r8,    0.466_r8, &
         0.0000_r8,  0.000_r8,    0.000_r8, &
         0.0000_r8,  0.000_r8,    0.000_r8 /) )

    ! coefficients for computing the asymmetry factor of
    !     ice, water, and rain particles (table 6)

    REAL(KIND=r8), PARAMETER, DIMENSION(nband,3) :: aig = RESHAPE(  &
         SHAPE = (/ nband, 3 /), SOURCE = (/ &
         0.74935228_r8, 0.76098937_r8, 0.84090400_r8, &
         0.00119715_r8, 0.00141864_r8, 0.00126222_r8, &
         -0.00000367_r8,-0.00000396_r8,-0.00000385_r8 /) )

    REAL(KIND=r8), PARAMETER, DIMENSION(nband,3) :: awg = RESHAPE( &
         SHAPE = (/ nband, 3 /), SOURCE = (/ &
         0.79375035_r8, 0.74513197_r8, 0.83530748_r8, &
         0.00832441_r8, 0.01370071_r8, 0.00257181_r8, &
         -0.00023263_r8,-0.00038203_r8, 0.00005519_r8 /) )

    REAL(KIND=r8), PARAMETER, DIMENSION(nband,3) :: arg = RESHAPE( &
         SHAPE = (/ nband, 3 /), SOURCE = (/ &
         0.891_r8,  0.948_r8,    0.971_r8, &
         0.0000_r8,  0.000_r8,    0.000_r8, &
         0.0000_r8,  0.000_r8,    0.000_r8 /) )

    ! input variables
    INTEGER, INTENT(IN) :: m,np
    INTEGER, INTENT(IN) :: ict(m),icb(m)
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: taual, ssaal, asyal
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3) :: taucld, reff
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: fcld, wh, dp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: rsirbm, rsirdf, cosz

    ! output variables
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx,flc,flx_d,flx_u, flc_d,flc_u
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m) :: fdirir,fdifir,fdirir_c,fdifir_c

    ! local  variables
    INTEGER :: i,k,ib,ik,iv
    INTEGER :: ih1,ih2,im1,im2,is1,is2
    REAL(KIND=r8)    :: taurs,tauwv
    REAL(KIND=r8)    :: taux,reff1,reff2,w1,w2,w3,g1,g2,g3

    REAL(KIND=r8), DIMENSION(m)        :: dsm, fsdir, fsdif, fsdir_c, fsdif_c,  ssaclt
    REAL(KIND=r8), DIMENSION(m,3)      :: cc
    REAL(KIND=r8), DIMENSION(m,np)     :: tauclb,tauclf,asycl,tautof,ssatof,asytof,ssacl
    REAL(KIND=r8), DIMENSION(m,np)     :: tausto,ssatau,asysto,tautob,ssatob,asytob
    REAL(KIND=r8), DIMENSION(m,np)     :: dum,rrt,ttt,tdt,rst,tst
    REAL(KIND=r8), DIMENSION(m,np+1)   :: fall,fclr,fall_d,fall_u,fclr_d,fclr_u
    REAL(KIND=r8), DIMENSION(m,np+1,2) :: rr,tt,td,rs,ts

    ! Initialize local vectors
     taurs=0.0_r8;tauwv=0.0_r8
     taux=0.0_r8;reff1=0.0_r8;reff2=0.0_r8;w1=0.0_r8;w2=0.0_r8;w3=0.0_r8;g1=0.0_r8;g2=0.0_r8;g3=0.0_r8
     dsm=0.0_r8; fsdir=0.0_r8; fsdif=0.0_r8; fsdir_c=0.0_r8; fsdif_c=0.0_r8;  ssaclt=0.0_r8
     cc=0.0_r8
     tauclb=0.0_r8;tauclf=0.0_r8;asycl=0.0_r8;tautof=0.0_r8;ssatof=0.0_r8;asytof=0.0_r8;ssacl=0.0_r8
     tausto=0.0_r8;ssatau=0.0_r8;asysto=0.0_r8;tautob=0.0_r8;ssatob=0.0_r8;asytob=0.0_r8
     dum=0.0_r8;rrt=0.0_r8;ttt=0.0_r8;tdt=0.0_r8;rst=0.0_r8;tst=0.0_r8
     fall=0.0_r8;fclr=0.0_r8;fall_d=0.0_r8;fall_u=0.0_r8;fclr_d=0.0_r8;fclr_u=0.0_r8
     rr=0.0_r8;tt=0.0_r8;td=0.0_r8;rs=0.0_r8;ts=0.0_r8
    ! Subroutine starts here
    !
    ! initialize surface fluxes, reflectances, and transmittances.
    ! the reflectance and transmittance of the clear and cloudy portions
    ! of a layer are denoted by 1 and 2, respectively.
    ! cc is the maximum cloud cover in each of the high, middle, and low
    ! cloud groups.
    ! 1/dsm=1/cos(53)=1.66

    dsm=0.602_r8
    fdirir=0.0_r8
    fdifir=0.0_r8
    fdirir_c=0.0_r8
    fdifir_c=0.0_r8
    rr(1:m,np+1,1)=rsirbm(1:m)
    rr(1:m,np+1,2)=rsirbm(1:m)
    rs(1:m,np+1,1)=rsirdf(1:m)
    rs(1:m,np+1,2)=rsirdf(1:m)
    td=0.0_r8
    tt=0.0_r8
    ts=0.0_r8
    cc=0.0_r8

    ! integration over spectral bands
    DO ib=1,nband
       iv=ib+8
       ! scale cloud optical thickness in each layer from taucld (with
       !     cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
       !     tauclb is the scaled optical thickness for beam radiation and
       !     tauclf is for diffuse radiation.
       CALL cldscale( &
                    m           , & !INTEGER, INTENT(IN) :: m
                    np          , & !INTEGER, INTENT(IN) :: np
                    cosz        , & !REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: cosz
                    fcld        , & !REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: fcld
                    taucld      , & !REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3) :: taucld
                    ict         , & !INTEGER, INTENT(IN) :: ict
                    icb         , & !INTEGER, INTENT(IN) :: icb
                    cc          , & !REAL(KIND=r8), INTENT(OUT), DIMENSION(m,3)  :: cc
                    tauclb      , & !REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclb
                    tauclf        ) !REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclf
       ! compute cloud single scattering albedo and asymmetry factor
       !     for a mixture of ice and liquid particles.
       !     eqs.(4.6)-(4.8), (6.2)-(6.4)
       DO k=1,np
          DO i=1,m
             ssaclt(i)=0.99999_r8
             asycl(i,k)=1.0_r8

             taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
             IF (taux.gt.0.02_r8 .and. fcld(i,k).gt.0.01_r8) THEN
                reff1=min(reff(i,k,1),130.0_r8)
                reff2=min(reff(i,k,2),20.0_r8)

                w1=(1.0_r8-(aia(ib,1)+(aia(ib,2)+ &
                     aia(ib,3)*reff1)*reff1))*taucld(i,k,1)
                w2=(1.0_r8-(awa(ib,1)+(awa(ib,2)+ &
                     awa(ib,3)*reff2)*reff2))*taucld(i,k,2)
                w3=(1.0_r8- ara(ib,1))*taucld(i,k,3)
                ssaclt(i)=(w1+w2+w3)/taux

                g1=(aig(ib,1)+(aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
                g2=(awg(ib,1)+(awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
                g3= arg(ib,1)*w3
                asycl(i,k)=(g1+g2+g3)/(w1+w2+w3)
             ENDIF
          ENDDO

          DO i=1,m
             ssacl(i,k)=ssaclt(i)
          ENDDO
!          DO i=1,m
!             asycl(i,k)=asyclt(i)
!          ENDDO
       ENDDO

       ! integration over the k-distribution function
       DO ik=1,nk
          DO k=1,np
             DO i=1,m
                taurs=ry2(ib)*dp(i,k)
                tauwv=xk2(ik)*wh(i,k)
                ! compute clear-sky optical thickness, single scattering albedo,
                !and asymmetry factor. eqs.(6.2)-(6.4)
                tausto(i,k)=taurs+tauwv+taual(i,k,iv)+1.0e-8_r8
                ssatau(i,k)=ssaal(i,k,iv)*taual(i,k,iv)+taurs
                asysto(i,k)=asyal(i,k,iv)*ssaal(i,k,iv)*taual(i,k,iv)

                ! compute reflectance and transmittance of the clear portion of a layer
                tautob(i,k)=tausto(i,k)
                ssatob(i,k)=ssatau(i,k)/tautob(i,k)+1.0e-8_r8
                ssatob(i,k)=min(ssatob(i,k),0.999999_r8)
                asytob(i,k)=asysto(i,k)/(ssatob(i,k)*tautob(i,k))
             ENDDO
          ENDDO

          ! for direct incident radiation

          CALL deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

          ! diffuse incident radiation is approximated by beam radiation with
          !  an incident angle of 53 degrees, eqs. (6.5) and (6.6)

          CALL deledd (m,np,tautob,ssatob,asytob,dsm,rst,tst,dum)

          DO k=1,np
             DO i=1,m
                rr(i,k,1)=rrt(i,k)
                tt(i,k,1)=ttt(i,k)
                td(i,k,1)=tdt(i,k)
                rs(i,k,1)=rst(i,k)
                ts(i,k,1)=tst(i,k)
             ENDDO
          ENDDO

          ! compute reflectance and transmittance of the cloudy portion of a layer
          DO k=1,np
             DO i=1,m
                ! for direct incident radiation. eqs.(6.2)-(6.4)
                tautob(i,k)=tausto(i,k)+tauclb(i,k)
                ssatob(i,k)=(ssatau(i,k)+ssacl(i,k)*tauclb(i,k)) &
                     /tautob(i,k)+1.0e-8_r8
                ssatob(i,k)=min(ssatob(i,k),0.999999_r8)
                asytob(i,k)=(asysto(i,k)+asycl(i,k)*ssacl(i,k)*tauclb(i,k)) &
                     /(ssatob(i,k)*tautob(i,k))

                ! for diffuse incident radiation
                tautof(i,k)=tausto(i,k)+tauclf(i,k)
                ssatof(i,k)=(ssatau(i,k)+ssacl(i,k)*tauclf(i,k)) &
                     /tautof(i,k)+1.0e-8_r8
                ssatof(i,k)=min(ssatof(i,k),0.999999_r8)
                asytof(i,k)=(asysto(i,k)+asycl(i,k)*ssacl(i,k)*tauclf(i,k)) &
                     /(ssatof(i,k)*tautof(i,k))
             ENDDO
          ENDDO

          ! for direct incident radiation

          CALL deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

          ! diffuse incident radiation is approximated by beam radiation with
          !  an incident angle of 53 degrees, eqs.(6.5) and (6.6)

          CALL deledd (m,np,tautof,ssatof,asytof,dsm,rst,tst,dum)

          DO k=1,np
             DO i=1,m
                rr(i,k,2)=rrt(i,k)
                tt(i,k,2)=ttt(i,k)
                td(i,k,2)=tdt(i,k)
                rs(i,k,2)=rst(i,k)
                ts(i,k,2)=tst(i,k)
             ENDDO
          ENDDO
          ! flux calculations

          ! initialize clear-sky flux (fclr), all-sky flux (fall),
          !  and surface downward fluxes (fsdir and fsdif)
          !hmjb they are initialized inside cldfx()
          !     fclr   = 0.0_r8
          !     fall   = 0.0_r8
          !     fclr_d = 0.0_r8
          !     fclr_u = 0.0_r8
          !     fall_d = 0.0_r8
          !     fall_u = 0.0_r8
          !
          !     fsdir   = 0.0_r8
          !     fsdif   = 0.0_r8
          !     fsdir_c = 0.0_r8
          !     fsdif_c = 0.0_r8

          ! for clear- and all-sky fluxes
          !  the all-sky flux, fall is the summation inside the brackets
          !  of eq. (7.11)

          ih1=1
          ih2=2
          im1=1
          im2=2
          is1=1
          is2=2

          CALL cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,  &
               cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif, fsdir_c,fsdif_c, &
               fclr_d,fclr_u,fall_d,fall_u)! new

          ! flux integration following eq. (6.1)

          DO k=1,np+1
             DO i=1,m
                flx_d(i,k)=flx_d(i,k)+fall_d(i,k)*hk2(ib,ik) !new
                flx_u(i,k)=flx_u(i,k)+fall_u(i,k)*hk2(ib,ik) !new
                flx(i,k) = flx(i,k)+fall(i,k)*hk2(ib,ik)
             ENDDO

             DO i=1,m
                flc(i,k) = flc(i,k)+fclr(i,k)*hk2(ib,ik)
                flc_d(i,k)=flc_d(i,k)+fclr_d(i,k)*hk2(ib,ik) !new
                flc_u(i,k)=flc_u(i,k)+fclr_u(i,k)*hk2(ib,ik) !new
             ENDDO
          ENDDO

          ! compute downward surface fluxes in the ir region

          DO i=1,m
             fdirir(i) = fdirir(i)+fsdir(i)*hk2(ib,ik)
             fdifir(i) = fdifir(i)+fsdif(i)*hk2(ib,ik)
             fdirir_c(i) = fdirir_c(i)+fsdir_c(i)*hk2(ib,ik)
             fdifir_c(i) = fdifir_c(i)+fsdif_c(i)*hk2(ib,ik)
          ENDDO

       ENDDO ! integration over the k-distribution function
    ENDDO ! integration over spectral bands

  END SUBROUTINE solircld

  !
  ! Subroutine: Soluvcld
  !
  ! $Author: pkubota $
  ! Modifications: T. Tarasova, 2005
  ! Modifications: H. Barbosa
  !
  ! Description:
  !  compute solar fluxes in the uv+par region. the spectrum is
  !  grouped into 8 bands:
  !
  !   band     micrometer
  !
  !    uv-c    1.     .175 - .225
  !    2.     .225 - .245
  !   .260 - .280
  !    3.     .245 - .260
  !
  !    uv-b    4.     .280 - .295
  !    5.     .295 - .310
  !    6.     .310 - .320
  !
  !    uv-a    7.     .320 - .400
  !
  !    par     8.     .400 - .700
  !
  ! Inputs: units        size
  !
  !      m: number of soundings n/d      1
  !     np: number of atmospheric layers n/d      1
  !     wh: layer scaled-water vapor content gm/cm^2      m*np
  !     oh: layer ozone content (cm-atm)stp  m*np
  !     dp: layer pressure thickness mb      m*np
  !  overcast: option for scaling cloud optical thickness n/d      1
  !   "true"  = scaling is not required
  !   "fasle" = scaling is required
  !  cldwater: input option for cloud optical thickness n/d      1
  !   "true"  = taucld is provided
  !   "false" = cwp is provided
  !    cwp: cloud water amount gm/m**2      m*np*3
  !   index 1 for ice particles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  ! taucld: cloud optical thickness n/d      m*np*3
  !   index 1 for ice particles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  !   reff: effective cloud-particle size   micrometer   m*np*3
  !   index 1 for ice paticles
  !   index 2 for liquid drops
  !   index 3 for rain drops
  !    ict: level index separating high and middle clouds   n/d      m
  !    icb: level indiex separating middle and low clouds   n/d      m
  !   fcld: cloud amount fraction     m*np
  !   cosz: cosine of solar zenith angle n/d      m
  !  taual: aerosol optical thickness n/d      m*np*11
  !  ssaal: aerosol single-scattering albedo n/d      m*np*11
  !  asyal: aerosol asymmetry factor n/d      m*np*11
  ! rsuvbm: uv+par surface albedo for beam radiation fraction     m
  ! rsuvdf: uv+par surface albedo for diffuse radiation fraction     m
  !
  ! Outputs: (updated parameters)
  !
  !    flx: all-sky   net downward flux    fraction  m*(np+1)
  !    flc: clear-sky net downward flux    fraction  m*(np+1)
  ! fdiruv: all-sky direct  downward uv  flux at the surface   fraction  m
  ! fdifuv: all-sky diffuse downward uv  flux at the surface   fraction  m
  !fdirpar: all-sky direct  downward par flux at the surface   fraction  m
  !fdifpar: all-sky diffuse downward par flux at the surface   fraction  m
  !
  ! Local Variables
  !
  !  tauclb: scaled cloud optical thickness for beam radiation    n/d   m*np
  !  tauclf: scaled cloud optical thickness for diffuse radiation n/d   m*np
  !
  !
  !
  !
  SUBROUTINE Soluvcld ( &
       m        , &!    INTEGER, INTENT(IN) :: m
       np       , &!    INTEGER, INTENT(IN) :: np
       wh       , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: wh
       oh       , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: oh
       dp       , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: dp
       taucld   , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: taucld
       reff     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: reff
       ict      , &!    INTEGER, INTENT(IN) :: ict(m)
       icb      , &!    INTEGER, INTENT(IN) :: icb(m)
       fcld     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: fcld
       cosz     , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m)       :: cosz
       taual    , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: taual
       ssaal    , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: ssaal
       asyal    , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: asyal
       rsuvbm   , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m)       :: rsuvbm
       rsuvdf   , &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m)       :: rsuvdf
       flx      , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx
       flc      , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc
       fdiruv   , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdiruv
       fdifuv   , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdifuv
       fdirpar  , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdirpar
       fdifpar  , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdifpar
       fdiruv_c , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdiruv_c
       fdifuv_c , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdifuv_c
       fdirpar_c, &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdirpar_c
       fdifpar_c, &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)        :: fdifpar_c
       flx_d    , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_d
       flx_u    , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_u
       flc_d    , &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc_d
       flc_u      )!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc_u


    IMPLICIT NONE

    ! input variables
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: np
    INTEGER, INTENT(IN) :: ict(m)
    INTEGER, INTENT(IN) :: icb(m)
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: taual
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: ssaal
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: asyal
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: taucld
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: reff
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: fcld
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: wh
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: oh
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)    :: dp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m)       :: rsuvbm
    REAL(KIND=r8), INTENT(IN), DIMENSION(m)       :: rsuvdf
    REAL(KIND=r8), INTENT(IN), DIMENSION(m)       :: cosz

    ! output variables
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_d
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_u
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc_d
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc_u
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdiruv
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdifuv
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdirpar
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdifpar
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdiruv_c
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdifuv_c
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdirpar_c
    REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)      :: fdifpar_c

    ! local  variables
    INTEGER :: k,ib,i
    INTEGER :: ih1,ih2,im1,im2,is1,is2
    REAL(KIND=r8)    :: taurs,tauoz,tauwv
    REAL(KIND=r8)    :: taux,reff1,reff2,g1,g2,g3

    REAL(KIND=r8), DIMENSION(m)        :: dsm, fsdir, fsdif, fsdir_c, fsdif_c
    REAL(KIND=r8), DIMENSION(m,3)      :: cc
    REAL(KIND=r8), DIMENSION(m,np)     :: tauclb,tauclf,asycl,tautof,ssatof,asytof
    REAL(KIND=r8), DIMENSION(m,np)     :: tausto,ssatau,asysto,tautob,ssatob,asytob
    REAL(KIND=r8), DIMENSION(m,np)     :: dum,rrt,ttt,tdt,rst,tst
    REAL(KIND=r8), DIMENSION(m,np+1)   :: fall,fclr,fall_d,fall_u,fclr_d,fclr_u
    REAL(KIND=r8), DIMENSION(m,np+1,2) :: rr,tt,td,rs,ts

   ! INTEGER :: i


    ! initialize local vectors
    ! subroutine starts here

    !  initialize fdiruv, fdifuv, surface reflectances and transmittances.
    !  the reflectance and transmittance of the clear and cloudy portions
    !  of a layer are denoted by 1 and 2, respectively.
    !  cc is the maximum cloud cover in each of the high, middle, and low
    !  cloud groups.
    ! 1/dsm=1/cos(53) = 1.66
    ih1=0;ih2=0;im1=0;im2=0;is1=0;is2=0;
    dsm=0.0_r8; fsdir=0.0_r8; fsdif=0.0_r8; fsdir_c=0.0_r8; fsdif_c=0.0_r8
    cc=0.0_r8
    tauclb=0.0_r8;tauclf=0.0_r8;asycl=0.0_r8;tautof=0.0_r8;ssatof=0.0_r8;asytof=0.0_r8
    tausto=0.0_r8;ssatau=0.0_r8;asysto=0.0_r8;tautob=0.0_r8;ssatob=0.0_r8;asytob=0.0_r8
    dum=0.0_r8;rrt=0.0_r8;ttt=0.0_r8;tdt=0.0_r8;rst=0.0_r8;tst=0.0_r8
    fall=0.0_r8;fclr=0.0_r8;fall_d=0.0_r8;fall_u=0.0_r8;fclr_d=0.0_r8;fclr_u=0.0_r8
    rr=0.0_r8;tt=0.0_r8;td=0.0_r8;rs=0.0_r8;ts=0.0_r8
    
    dsm=0.602_r8
    fdiruv=0.0_r8
    fdifuv=0.0_r8
    fdirpar=0.0_r8
    fdifpar=0.0_r8
    fdiruv_c=0.0_r8
    fdifuv_c=0.0_r8
    fdirpar_c=0.0_r8
    fdifpar_c=0.0_r8
    DO i=1,m
       rr(i,np+1,1)=rsuvbm(i)
       rr(i,np+1,2)=rsuvbm(i)
       rs(i,np+1,1)=rsuvdf(i)
       rs(i,np+1,2)=rsuvdf(i)
    ENDDO
    td=0.0_r8
    tt=0.0_r8
    ts=0.0_r8
    cc=0.0_r8


    ! scale cloud optical thickness in each layer from taucld (with
    ! cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
    ! tauclb is the scaled optical thickness for beam radiation and
    ! tauclf is for diffuse radiation (see section 7).

    CALL cldscale( &
                 m          , & !INTEGER, INTENT(IN) :: m
                 np         , & !INTEGER, INTENT(IN) :: np
                 cosz       , & !REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: cosz
                 fcld       , & !REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: fcld
                 taucld     , & !REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3) :: taucld
                 ict        , & !INTEGER, INTENT(IN) :: ict
                 icb        , & !INTEGER, INTENT(IN) :: icb
                 cc         , & !REAL(KIND=r8), INTENT(OUT), DIMENSION(m,3)  :: cc
                 tauclb     , & !REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclb
                 tauclf       ) !REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclf

    ! cloud asymmetry factor for a mixture of liquid and ice particles.
    ! unit of reff is micrometers. eqs. (4.8) and (6.4)

    DO k=1,np
       DO i=1,m
          asycl(i,k)=1.0_r8
          taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
          IF (taux.gt.0.02_r8 .and. fcld(i,k).gt.0.01_r8) THEN
             reff1=min(reff(i,k,1),130.0_r8)
             reff2=min(reff(i,k,2),20.0_r8)
             g1=(aig3(1)+(aig3(2)+aig3(3)*reff1)*reff1)*taucld(i,k,1)
             g2=(awg3(1)+(awg3(2)+awg3(3)*reff2)*reff2)*taucld(i,k,2)
             g3= arg3(1)*taucld(i,k,3)
             asycl(i,k)=(g1+g2+g3)/taux
          ENDIF
       ENDDO

!       DO i=1,m
!          asycl(i,k)=asyclt(i)
!       ENDDO

    ENDDO

    ! integration over spectral bands
    DO ib=1,nband
       DO k=1,np
          DO i=1,m
             ! compute rayleigh, ozone and water vapor optical thicknesses
             taurs=ry3(ib)*dp(i,k)
             tauoz=zk3(ib)*oh(i,k)
             tauwv=wk3(ib)*wh(i,k)

             ! compute clear-sky optical thickness, single scattering albedo,
             ! and asymmetry factor (eqs. 6.2-6.4)
             tausto(i,k)=taurs+tauoz+tauwv+taual(i,k,ib)+1.0e-8_r8
             ssatau(i,k)=ssaal(i,k,ib)*taual(i,k,ib)+taurs
             asysto(i,k)=asyal(i,k,ib)*ssaal(i,k,ib)*taual(i,k,ib)

             ! compute reflectance and transmittance of the clear portion of a layer
             tautob(i,k)=tausto(i,k)
             ssatob(i,k)=ssatau(i,k)/tautob(i,k)+1.0e-8_r8
             ssatob(i,k)=min(ssatob(i,k),0.999999_r8)
             asytob(i,k)=asysto(i,k)/(ssatob(i,k)*tautob(i,k))
          ENDDO
       ENDDO

       ! for direct incident radiation
       CALL deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

       ! diffuse incident radiation is approximated by beam radiation with
       ! an incident angle of 53 degrees, eqs. (6.5) and (6.6)
       CALL deledd (m,np,tautob,ssatob,asytob,dsm,rst,tst,dum)

       DO k=1,np
          DO i=1,m
             rr(i,k,1)=rrt(i,k)
             tt(i,k,1)=ttt(i,k)
             td(i,k,1)=tdt(i,k)
             rs(i,k,1)=rst(i,k)
             ts(i,k,1)=tst(i,k)
          ENDDO
       ENDDO

       ! compute reflectance and transmittance of the cloudy portion of a layer
       DO k=1,np
          DO i=1,m

             ! for direct incident radiation
             ! the effective layer optical properties. eqs. (6.2)-(6.4)
             tautob(i,k)=tausto(i,k)+tauclb(i,k)
             ssatob(i,k)=(ssatau(i,k)+tauclb(i,k))/tautob(i,k)+1.0e-8_r8
             ssatob(i,k)=min(ssatob(i,k),0.999999_r8)
             asytob(i,k)=(asysto(i,k)+asycl(i,k)*tauclb(i,k)) &
                  /(ssatob(i,k)*tautob(i,k))

             ! for diffuse incident radiation
             tautof(i,k)=tausto(i,k)+tauclf(i,k)
             ssatof(i,k)=(ssatau(i,k)+tauclf(i,k))/tautof(i,k)+1.0e-8_r8
             ssatof(i,k)=min(ssatof(i,k),0.999999_r8)
             asytof(i,k)=(asysto(i,k)+asycl(i,k)*tauclf(i,k))  &
                  /(ssatof(i,k)*tautof(i,k))
          ENDDO
       ENDDO

       ! for direct incident radiation
       ! note that the cloud optical thickness is scaled differently for direct
       ! and diffuse insolation, eqs. (7.3) and (7.4).

       CALL deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

       ! diffuse incident radiation is approximated by beam radiation with
       ! an incident angle of 53 degrees, eqs. (6.5) and (6.6)

       CALL deledd (m,np,tautof,ssatof,asytof,dsm,rst,tst,dum)

       DO k=1,np
          DO i=1,m
             rr(i,k,2)=rrt(i,k)
             tt(i,k,2)=ttt(i,k)
             td(i,k,2)=tdt(i,k)
             rs(i,k,2)=rst(i,k)
             ts(i,k,2)=tst(i,k)
          ENDDO
       ENDDO

       ! flux calculations

       ! initialize clear-sky flux (fclr), all-sky flux (fall),
       ! and surface downward fluxes (fsdir and fsdif)

       !hmjb they are initialized inside cldfx()
       !     fclr   = 0.0_r8
       !     fall   = 0.0_r8
       !     fclr_d = 0.0_r8
       !     fclr_u = 0.0_r8
       !     fall_d = 0.0_r8
       !     fall_u = 0.0_r8
       !
       !     fsdir   = 0.0_r8
       !     fsdif   = 0.0_r8
       !     fsdir_c = 0.0_r8
       !     fsdif_c = 0.0_r8

       ! for clear- and all-sky fluxes
       ! the all-sky flux, fall is the summation inside the brackets
       ! of eq. (7.11)
       ih1=1
       ih2=2
       im1=1
       im2=2
       is1=1
       is2=2

       CALL cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2, &
            cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,fsdir_c,fsdif_c, &
            fclr_d,fclr_u,fall_d,fall_u)   ! new

       ! flux integration, eq. (6.1)

       DO k=1,np+1
          DO i=1,m
             flx  (i,k)=flx  (i,k)+fall  (i,k)*hk3(ib)
             flx_d(i,k)=flx_d(i,k)+fall_d(i,k)*hk3(ib) !new
             flx_u(i,k)=flx_u(i,k)+fall_u(i,k)*hk3(ib) !new

             flc  (i,k)=flc  (i,k)+fclr  (i,k)*hk3(ib)
             flc_d(i,k)=flc_d(i,k)+fclr_d(i,k)*hk3(ib) !new
             flc_u(i,k)=flc_u(i,k)+fclr_u(i,k)*hk3(ib) !new
          ENDDO
       ENDDO

       ! compute direct and diffuse downward surface fluxes in the uv
       ! and par regions
       IF(ib.lt.8) THEN
          DO i=1,m
             fdiruv(i)=fdiruv(i)+fsdir(i)*hk3(ib)
             fdifuv(i)=fdifuv(i)+fsdif(i)*hk3(ib)
             fdiruv_c(i) = fdiruv_c(i)+fsdir_c(i)*hk3(ib)
             fdifuv_c(i) = fdifuv_c(i)+fsdif_c(i)*hk3(ib)
          ENDDO
       ELSE
          DO i=1,m
             fdirpar(i)=fsdir(i)*hk3(ib)
             fdifpar(i)=fsdif(i)*hk3(ib)
             fdirpar_c(i) = fsdir_c(i)*hk3(ib)
             fdifpar_c(i) = fsdif_c(i)*hk3(ib)
          ENDDO
       ENDIF

    ENDDO ! integration over spectral bands

  end SUBROUTINE soluvcld

  !
  ! Subroutine: CLDSCALE
  !
  ! $Author: pkubota $
  !
  ! Description:
  !this SUBROUTINE computes the high, middle, and low cloud
  ! amounts and scales the cloud optical thickness (section 7)
  !
  !to simplify calculations in a cloudy atmosphere, clouds are
  ! grouped into high, middle and low clouds separated by the levels
  ! ict and icb (level 1 is the top of the model atmosphere).
  !
  !within each of the three groups, clouds are assumed maximally
  ! overlapped, and the cloud cover (cc) of a group is the maximum
  ! cloud cover of all the layers in the group.  the optical thickness
  ! (taucld) of a given layer is then scaled to new values (tauclb and
  ! tauclf) so that the layer reflectance corresponding to the cloud
  ! cover cc is the same as the original reflectance with optical
  ! thickness taucld and cloud cover fcld.
  !
  ! Inputs:
  !   m:  number of atmospheric soundings
  !  np:  number of atmospheric layers
  !cosz:  cosine of the solar zenith angle
  !fcld:  fractional cloud cover
  ! taucld:  cloud optical thickness
  ! ict:  index separating high and middle clouds
  ! icb:  index separating middle and low clouds
  !
  ! Outputs:
  !   cc:  fractional cover of high, middle, and low cloud groups
  !  tauclb:  scaled cloud optical thickness for direct  radiation
  !  tauclf:  scaled cloud optical thickness for diffuse radiation
  !
  !
  !
  SUBROUTINE cldscale (&
         m             , &!INTEGER, INTENT(IN) :: m
         np            , &!INTEGER, INTENT(IN) :: np
         cosz          , &!REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: cosz
         fcld          , &!REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: fcld
         taucld        , &!REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3) :: taucld
         ict           , &!INTEGER, INTENT(IN) :: ict
         icb           , &!INTEGER, INTENT(IN) :: icb
         cc            , &!REAL(KIND=r8), INTENT(OUT), DIMENSION(m,3)  :: cc
         tauclb        , &!REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclb
         tauclf          )!REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclf

    IMPLICIT NONE

    ! parameters
    !INTEGER, PARAMETER :: nm=11
    !INTEGER, PARAMETER :: nt=9
    !INTEGER, PARAMETER :: na=11

    REAL(KIND=r8), PARAMETER :: dm=0.1_r8     ! size of cosz-interval
    REAL(KIND=r8), PARAMETER :: dt=0.30103_r8  ! size of taucld-interval
    REAL(KIND=r8), PARAMETER :: da=0.1_r8      ! size of cloud amount-interval
    REAL(KIND=r8), PARAMETER :: t1=-0.9031_r8

    INTEGER i
    ! include the pre-computed table of mcai for scaling the cloud optical
    ! thickness under the assumption that clouds are maximally overlapped

    ! caib is for scaling the cloud optical thickness for direct radiation
    !REAL(KIND=r8) :: caib(nm,nt,na)

    ! caif is for scaling the cloud optical thickness for diffuse radiation
    !REAL(KIND=r8) :: caif(nt,na)

    !INCLUDE "mcai.data90"

    ! input variables
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: np
    INTEGER, INTENT(IN) :: ict(m)
    INTEGER, INTENT(IN) :: icb(m)

    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: cosz
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: fcld
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3) :: taucld

    ! output variables
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,3)  :: cc
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclb
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: tauclf

    !  local  variables
    INTEGER :: k,im,it,ia,kk
    REAL(KIND=r8) :: fm,ft,fa,xai,taux

    ! initialize local vectors and output variables
    fm=0.0_r8;ft=0.0_r8;fa=0.0_r8;xai=0.0_r8;taux=0.0_r8
    ! subroutine starts here 
    cc=0.0_r8
    tauclb=0.0_r8
    tauclf=0.0_r8
    ! clouds within each of the high, middle, and low clouds are assumed
    !     to be maximally overlapped, and the cloud cover (cc) for a group
    !     (high, middle, or low) is the maximum cloud cover of all the layers
    !     within a group

    DO i=1,m
       cc(i,1)=0.0_r8
       cc(i,2)=0.0_r8
       cc(i,3)=0.0_r8
    ENDDO

!    DO k=1,MAXVAL(ict)-1
    DO k=1,np
       DO i=1,m
          IF(k <= ict(i)-1)THEN
             cc(i,1)=max(cc(i,1),fcld(i,k))
          END IF
       ENDDO
    ENDDO

!    DO k=MINVAL(ict),MAXVAL(icb)-1
    DO k=1,np
       DO i=1,m
         IF( k >= ict(i) .and. k <= icb(i)-1)THEN
          cc(i,2)=max(cc(i,2),fcld(i,k))
         END IF
       ENDDO
    ENDDO

!    DO k=MINVAL(icb),np
    DO k=1,np
       DO i=1,m
          IF(k>=icb(i))THEN
             cc(i,3)=max(cc(i,3),fcld(i,k))
          END IF
       ENDDO
    ENDDO

    ! scale the cloud optical thickness.
    !     taucld(i,k,1) is the optical thickness for ice particles
    !     taucld(i,k,2) is the optical thickness for liquid particles
    !     taucld(i,k,3) is the optical thickness for rain drops

    DO k=1,np
       DO i=1,m

       IF(k.lt.ict(i)) THEN
          kk=1
       ELSEIF(k.ge.ict(i) .and. k.lt.icb(i)) THEN
          kk=2
       ELSE
          kk=3
       ENDIF

!       DO i=1,m
          !hmjb ta feito no comeco
          tauclb(i,k) = 0.0_r8
          tauclf(i,k) = 0.0_r8
          taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)

          IF (taux.gt.0.02_r8 .and. fcld(i,k).gt.0.01_r8) THEN

             ! normalize cloud cover following eq. (7.8)

             fa=fcld(i,k)/cc(i,kk)

             ! table look-up

             taux=min(taux,32.0_r8)

             fm=cosz(i)/dm
             ft=(log10(taux)-t1)/dt
             fa=fa/da

             im=int(fm+1.5_r8)
             it=int(ft+1.5_r8)
             ia=int(fa+1.5_r8)

             im=max(im,2)
             it=max(it,2)
             ia=max(ia,2)

             im=min(im,nm-1)
             it=min(it,nt-1)
             ia=min(ia,na-1)

             fm=fm-real(im-1,kind=r8)
             ft=ft-real(it-1,kind=r8)
             fa=fa-real(ia-1,kind=r8)

             ! scale cloud optical thickness for beam radiation following eq. (7.3)
             !     the scaling factor, xai, is a function of the solar zenith
             !     angle, optical thickness, and cloud cover.

             xai=    (-caib(im-1,it,ia)*(1.0_r8-fm)+  &
                  caib(im+1,it,ia)*(1.0_r8+fm))*fm*0.5_r8+caib(im,it,ia)*(1.0_r8-fm*fm)

             xai=xai+(-caib(im,it-1,ia)*(1.0_r8-ft)+ &
                  caib(im,it+1,ia)*(1.0_r8+ft))*ft*0.5_r8+caib(im,it,ia)*(1.0_r8-ft*ft)

             xai=xai+(-caib(im,it,ia-1)*(1.0_r8-fa)+ &
                  caib(im,it,ia+1)*(1.0_r8+fa))*fa*0.5_r8+caib(im,it,ia)*(1.0_r8-fa*fa)

             xai= xai-2.0_r8*caib(im,it,ia)
             xai=max(xai,0.0_r8)

             tauclb(i,k) = taux*xai

             ! scale cloud optical thickness for diffuse radiation following eq. (7.4)
             !     the scaling factor, xai, is a function of the cloud optical
             !     thickness and cover but not the solar zenith angle.

             xai=    (-caif(it-1,ia)*(1.0_r8-ft)+  &
                  caif(it+1,ia)*(1.0_r8+ft))*ft*0.5_r8+caif(it,ia)*(1.0_r8-ft*ft)

             xai=xai+(-caif(it,ia-1)*(1.0_r8-fa)+  &
                  caif(it,ia+1)*(1.0_r8+fa))*fa*0.5_r8+caif(it,ia)*(1.0_r8-fa*fa)

             xai= xai-caif(it,ia)
             xai=max(xai,0.0_r8)

             tauclf(i,k) = taux*xai

          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE cldscale

  !
  ! Subroutine: DELEDD
  !
  ! $Author: pkubota $
  !
  ! Description:
  !  uses the delta-eddington approximation to compute the
  !  bulk scattering properties of a single layer
  !  coded following king and harshvardhan (jas, 1986)
  !
  ! Inputs:
  !   m:  number of soundings
  !  np:  number of atmospheric layers
  ! tau:  optical thickness
  ! ssc:  single scattering albedo
  !  g0:  asymmetry factor
  ! cza:  cosine of the zenith angle
  !
  ! Outputs:
  !   rr:  reflection of the direct beam
  !   tt:  total diffuse transmission of the direct beam
  !   td:  direct transmission of the direct beam
  !
  !
  !
  SUBROUTINE deledd(m,np,tau,ssc,g0,cza,rr,tt,td)

    IMPLICIT NONE

    ! parameters
    REAL(KIND=r8), PARAMETER :: zero=0.0_r8
    REAL(KIND=r8), PARAMETER :: one=1.0_r8
    REAL(KIND=r8), PARAMETER :: two=2.0_r8
    REAL(KIND=r8), PARAMETER :: three=3.0_r8
    REAL(KIND=r8), PARAMETER :: four=4.0_r8
    REAL(KIND=r8), PARAMETER :: fourth=0.25_r8
    REAL(KIND=r8), PARAMETER :: seven=7.0_r8
    REAL(KIND=r8), PARAMETER :: thresh=1.0e-8_r8

    ! input variables
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: np

    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,np) :: tau
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,np) :: ssc
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,np) :: g0
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m) :: cza

    ! output variables
    REAL(KIND=r8),  INTENT(OUT), DIMENSION(m,np) :: rr
    REAL(KIND=r8),  INTENT(OUT), DIMENSION(m,np) :: tt
    REAL(KIND=r8),  INTENT(OUT), DIMENSION(m,np) :: td

    ! local  variables
    INTEGER :: i,k
    REAL(KIND=r8) :: zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2, &
         temp_all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4

     rr=0.0_r8
     tt=0.0_r8
     td=0.0_r8

     zth=0.0_r8;ff=0.0_r8;xx=0.0_r8;taup=0.0_r8;sscp=0.0_r8;gp=0.0_r8;gm1=0.0_r8;gm2=0.0_r8;gm3=0.0_r8;akk=0.0_r8;alf1=0.0_r8;alf2=0.0_r8
     temp_all=0.0_r8;bll=0.0_r8;st7=0.0_r8;st8=0.0_r8;cll=0.0_r8;dll=0.0_r8;fll=0.0_r8;ell=0.0_r8;st1=0.0_r8;st2=0.0_r8;st3=0.0_r8;st4=0.0_r8 
    ! initialize local vectors and output variables
    ! subroutine starts here
    DO k=1,np
       DO i=1,m
          zth = cza(i)
          !  delta-eddington scaling of single scattering albedo,
          !  optical thickness, and asymmetry factor,
          !  k & h eqs(27-29)

          ff  = g0(i,k)*g0(i,k)
          xx  = one-ff *ssc(i,k)
          taup= tau(i,k)*xx
          sscp= ssc(i,k)*(one-ff)/xx
          gp  = g0(i,k) /(one+g0(i,k))

          !  gamma1, gamma2, and gamma3. see table 2 and eq(26) k & h
          !  ssc and gp are the d-s single scattering
          !  albedo and asymmetry factor.

          xx  =  three*gp
          gm1 =  (seven - sscp*(four+xx))*fourth
          gm2 = -(one   - sscp*(four-xx))*fourth

          !  akk is k as defined in eq(25) of k & h

          akk = sqrt((gm1+gm2)*(gm1-gm2))

          xx  = akk * zth
          st7 = one - xx
          st8 = one + xx
          st3 = st7 * st8

          IF (abs(st3) .lt. thresh) THEN
             zth = zth + 0.001_r8
             xx  = akk * zth
             st7 = one - xx
             st8 = one + xx
             st3 = st7 * st8
          ENDIF

          !  extinction of the direct beam transmission
          td(i,k)  = exp(-taup/zth)

          !  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of k & h
          gm3  = (two - zth*three*gp)*fourth
          xx   = gm1 - gm2
          alf1 = gm1 - gm3 * xx
          alf2 = gm2 + gm3 * xx

          ! all is last term in eq(21) of k & h
          ! bll is last term in eq(22) of k & h

          xx  = akk * two
          temp_all = (gm3 - alf2 * zth    )*xx*td(i,k)
          bll = (one - gm3 + alf1*zth)*xx

          xx  = akk * gm3
          cll = (alf2 + xx) * st7
          dll = (alf2 - xx) * st8

          xx  = akk * (one-gm3)
          fll = (alf1 + xx) * st8
          ell = (alf1 - xx) * st7

          st2 = exp(-akk*taup)
          st4 = st2 * st2

          st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)

          ! rr is r-hat of eq(21) of k & h
          ! tt is diffuse part of t-hat of eq(22) of k & h

          rr(i,k) =( cll-dll*st4      -temp_all*st2)*st1
          tt(i,k) = - ((fll-ell*st4)*td(i,k)-bll*st2)*st1

          rr(i,k) = max(rr(i,k),zero)
          tt(i,k) = max(tt(i,k),zero)

          tt(i,k) = tt(i,k)+td(i,k)
       ENDDO
    ENDDO

  END SUBROUTINE deledd

  !
  ! Subroutine: RFLX
  !
  ! $Author: pkubota $
  !
  ! Description:
  !Computes the reduction of clear-sky downward solar flux
  !due to co2 absorption.
  !
  !
  !
  !
  SUBROUTINE rflx(m,np,swc,u1,du,nu,swh,w1,dw,nw,tbl,df)

    IMPLICIT NONE

    ! input variables
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: np
    INTEGER, INTENT(IN) :: nu
    INTEGER, INTENT(IN) :: nw

    REAL(KIND=r8),  INTENT(IN) :: u1
    REAL(KIND=r8),  INTENT(IN) :: du
    REAL(KIND=r8),  INTENT(IN) :: w1
    REAL(KIND=r8),  INTENT(IN) :: dw
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,np+1) :: swc
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,np+1) :: swh
    REAL(KIND=r8),  INTENT(IN), DIMENSION(nu,nw)  :: tbl

    ! output variables (updated)
    REAL(KIND=r8),  INTENT(INOUT), DIMENSION(m,np+1) :: df

    ! local  variables
    INTEGER :: i,k,ic,iw
    REAL(KIND=r8) :: temp_clog,wlog,dc,dd,x1,x2,y1,y2

    ! subroutine starts here
    temp_clog=0.0_r8;wlog=0.0_r8;dc=0.0_r8;dd=0.0_r8;x1=0.0_r8;x2=0.0_r8;y1=0.0_r8;y2=0.0_r8;
    ! table look-up for the reduction of clear-sky solar
    x1=u1-0.5_r8*du
    y1=w1-0.5_r8*dw

    DO k= 2, np+1
       DO i= 1, m
          temp_clog=swc(i,k)
          wlog=swh(i,k)
          ic=int( (temp_clog-x1)/du+1.0_r8)
          iw=int( (wlog-y1)/dw+1.0_r8)
          IF(ic.lt.2)ic=2
          IF(iw.lt.2)iw=2
          IF(ic.gt.nu)ic=nu
          IF(iw.gt.nw)iw=nw
          dc=temp_clog-real(ic-2,kind=r8)*du-u1
          dd=wlog-real(iw-2,kind=r8)*dw-w1
          x2=tbl(ic-1,iw-1)+(tbl(ic-1,iw)-tbl(ic-1,iw-1))/dw*dd
          y2=x2+(tbl(ic,iw-1)-tbl(ic-1,iw-1))/du*dc
          df(i,k)=df(i,k)+y2
       ENDDO
    ENDDO

  END SUBROUTINE rflx


  !
  ! Subroutine: CLDFLX
  !
  ! $Author: pkubota $
  !
  ! Description:
  !  compute upward and downward fluxes using a two-stream adding method
  !  following equations (6.9)-(6.16).
  !
  !  clouds are grouped into high, middle, and low clouds which are assumed
  !  randomly overlapped. it involves a maximum of 8 sets of calculations.
  !  in each set of calculations, each atmospheric layer is homogeneous,
  !  either totally filled with clouds or without clouds.
  !
  ! Inputs:
  !    m:  number of soundings
  !   np:  number of atmospheric layers
  !  ict:  the level separating high and middle clouds
  !  icb:  the level separating middle and low clouds
  ! ih1,ih2,
  ! im1,im2,
  ! is1,is2: indices for three group of clouds
  !   cc:  effective cloud covers for high, middle and low clouds
  !   rr:  reflection of a layer illuminated by beam radiation
  !   tt:  total diffuse transmission of a layer illuminated by beam radiation
  !   td:  direct beam transmission
  !   rs:  reflection of a layer illuminated by diffuse radiation
  !   ts:  transmission of a layer illuminated by diffuse radiation
  !
  ! Outputs:
  !  fclr:  clear-sky flux (downward minus upward)
  !  fall:  all-sky flux (downward minus upward)
  !  fsdir: surface direct downward flux
  !  fsdif: surface diffuse downward flux
  !
  !
  !
  !
  SUBROUTINE cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2, &
       cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,fsdir_c,fsdif_c, &
       fclr_d,fclr_u,fall_d,fall_u)   ! new

    IMPLICIT NONE

    ! input variables
    INTEGER,  INTENT(IN) :: m,np,ih1,ih2,im1,im2,is1,is2
    INTEGER,  INTENT(IN) :: ict(m),icb(m)
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,np+1,2) :: rr,tt,td,rs,ts
    REAL(KIND=r8),  INTENT(IN), DIMENSION(m,3) :: cc

    ! output variables
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np+1) :: fclr, fclr_u, fclr_d
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np+1) :: fall, fall_u, fall_d
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m)      :: fsdir, fsdif, fsdir_c, fsdif_c

    ! local  variables
    REAL(KIND=r8), DIMENSION(m,np+1,2,2) :: tta ! composite total transmittance illuminated by beam radiation
    REAL(KIND=r8), DIMENSION(m,np+1,2,2) :: tda ! composite transmittance illuminated by beam radiation
    REAL(KIND=r8), DIMENSION(m,np+1,2,2) :: rsa ! composite reflectance illuminated from below by diffuse radiation
    REAL(KIND=r8), DIMENSION(m,np+1,2,2) :: rra ! composite reflectance illuminated by beam radiation
    REAL(KIND=r8), DIMENSION(m,np+1,2,2) :: rxa ! the composite reflectance illuminated from above by diffuse radiation

    REAL(KIND=r8), DIMENSION(m,np+1) :: fdn   ! downward fluxes
    REAL(KIND=r8), DIMENSION(m,np+1) :: fup   ! upward fluxes
    REAL(KIND=r8), DIMENSION(m,np+1) :: flxdn ! net flux

    REAL(KIND=r8), DIMENSION(m) :: fdndir    ! direct  downward flux
    REAL(KIND=r8), DIMENSION(m) :: fdndif    ! diffuse downward flux
    REAL(KIND=r8), DIMENSION(m) :: ch, cm, ct

    REAL(KIND=r8) :: fupdif      ! diffuse upward flux
    REAL(KIND=r8) :: denm,xx,yy
    INTEGER :: i,k,ih,im,is
     tta  = 0.0_r8
     tda  = 0.0_r8
     rsa  = 0.0_r8
     rra  = 0.0_r8
     rxa  = 0.0_r8
     fdn    = 0.0_r8
     fup    = 0.0_r8
     flxdn  = 0.0_r8
     fdndir     = 0.0_r8
     fdndif     = 0.0_r8
     ch = 0.0_r8; cm = 0.0_r8; ct = 0.0_r8
     fupdif = 0.0_r8
     denm = 0.0_r8;xx = 0.0_r8;yy = 0.0_r8
     
    ! subroutine starts here
    fclr   = 0.0_r8
    fall   = 0.0_r8
    fclr_d = 0.0_r8
    fclr_u = 0.0_r8
    fall_d = 0.0_r8
    fall_u = 0.0_r8

    fsdir   = 0.0_r8
    fsdif   = 0.0_r8
    fsdir_c = 0.0_r8
    fsdif_c = 0.0_r8

    ! compute transmittances and reflectances for a composite of
    !     layers. layers are added one at a time, going down from the top.
    !     tta and rsa are computed from eqs. (6.10) and (6.12)

    ! for high clouds
    !     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition
    DO ih=ih1,ih2

       DO i=1,m
          tda(i,1,ih,1)=td(i,1,ih)
          tta(i,1,ih,1)=tt(i,1,ih)
          rsa(i,1,ih,1)=rs(i,1,ih)
          tda(i,1,ih,2)=td(i,1,ih)
          tta(i,1,ih,2)=tt(i,1,ih)
          rsa(i,1,ih,2)=rs(i,1,ih)
       ENDDO

!       DO k=2,MAXVAL(ict)-1
       DO k=2,np
          DO i=1,m
             IF(k<=ict(i)-1) then
             denm = ts(i,k,ih)/( 1.0_r8-rsa(i,k-1,ih,1)*rs(i,k,ih))
             tda(i,k,ih,1)= tda(i,k-1,ih,1)*td(i,k,ih)
             tta(i,k,ih,1)= tda(i,k-1,ih,1)*tt(i,k,ih) &
                  +(tda(i,k-1,ih,1)*rsa(i,k-1,ih,1)*rr(i,k,ih) &
                  +tta(i,k-1,ih,1)-tda(i,k-1,ih,1))*denm
             rsa(i,k,ih,1)= rs(i,k,ih)+ts(i,k,ih) &
                  *rsa(i,k-1,ih,1)*denm
             tda(i,k,ih,2)= tda(i,k,ih,1)
             tta(i,k,ih,2)= tta(i,k,ih,1)
             rsa(i,k,ih,2)= rsa(i,k,ih,1)
             END IF
          ENDDO

       ENDDO

       ! for middle clouds
       !     im=1 for clear-sky condition, im=2 for cloudy-sky condition
       DO im=im1,im2

!          DO k=MINVAL(ict),MAXVAL(icb)-1
           DO k=1,np
             DO i=1,m
                IF(k>=ict(i) .and. k <= icb(i)-1)THEN
                denm = ts(i,k,im)/( 1.0_r8-rsa(i,k-1,ih,im)*rs(i,k,im))
                tda(i,k,ih,im)= tda(i,k-1,ih,im)*td(i,k,im)
                tta(i,k,ih,im)= tda(i,k-1,ih,im)*tt(i,k,im)  &
                     +(tda(i,k-1,ih,im)*rsa(i,k-1,ih,im)*rr(i,k,im)  &
                     +tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
                rsa(i,k,ih,im)= rs(i,k,im)+ts(i,k,im)  &
                     *rsa(i,k-1,ih,im)*denm
                END IF
             ENDDO
          ENDDO

       ENDDO         ! end im loop
    ENDDO    ! end ih loop

    ! layers are added one at a time, going up from the surface.
    !     rra and rxa are computed from eqs. (6.9) and (6.11)

    ! for the low clouds
    !     is=1 for clear-sky condition, is=2 for cloudy-sky condition
    DO is=is1,is2

       DO i=1,m
          rra(i,np+1,1,is)=rr(i,np+1,is)
          rxa(i,np+1,1,is)=rs(i,np+1,is)
          rra(i,np+1,2,is)=rr(i,np+1,is)
          rxa(i,np+1,2,is)=rs(i,np+1,is)
       ENDDO

!       DO k=np,MinVAL(icb),-1
       DO k=np,1,-1
          DO i=1,m
             IF(k >=icb(i) )THEN  
             denm=ts(i,k,is)/( 1.0_r8-rs(i,k,is)*rxa(i,k+1,1,is) )
             rra(i,k,1,is)=rr(i,k,is)+(td(i,k,is)*rra(i,k+1,1,is) &
                  +(tt(i,k,is)-td(i,k,is))*rxa(i,k+1,1,is))*denm
             rxa(i,k,1,is)= rs(i,k,is)+ts(i,k,is)  &
                  *rxa(i,k+1,1,is)*denm
             rra(i,k,2,is)=rra(i,k,1,is)
             rxa(i,k,2,is)=rxa(i,k,1,is)
            END IF
          ENDDO
       ENDDO

       ! for middle clouds
       DO im=im1,im2
 !         DO k=icb-1,MAXVAL(ict),-1
          DO k=np,1,-1
             DO i=1,m
                IF(k<= icb(i)-1 .and. k >= ict(i))then
                denm=ts(i,k,im)/( 1.0_r8-rs(i,k,im)*rxa(i,k+1,im,is) )
                rra(i,k,im,is)= rr(i,k,im)+(td(i,k,im)*rra(i,k+1,im,is) &
                     +(tt(i,k,im)-td(i,k,im))*rxa(i,k+1,im,is))*denm
                rxa(i,k,im,is)= rs(i,k,im)+ts(i,k,im) &
                     *rxa(i,k+1,im,is)*denm
                END IF
             ENDDO
          ENDDO

       ENDDO       ! end im loop
    ENDDO  ! end is loop

    ! integration over eight sky situations.
    !     ih, im, is denotes high, middle and low cloud groups.
    DO ih=ih1,ih2

       ! clear portion
       IF(ih.eq.1) THEN
          DO i=1,m
             ch(i)=1.0_r8-cc(i,1)
          ENDDO
       ELSE
          ! cloudy portion
          DO i=1,m
             ch(i)=cc(i,1)
          ENDDO
       ENDIF

       DO im=im1,im2
          ! clear portion
          IF(im.eq.1) THEN
             DO i=1,m
                cm(i)=ch(i)*(1.0_r8-cc(i,2))
             ENDDO
          ELSE
             ! cloudy portion
             DO i=1,m
                cm(i)=ch(i)*cc(i,2)
             ENDDO
          ENDIF

          DO is=is1,is2
             ! clear portion
             IF(is.eq.1) THEN
                DO i=1,m
                   ct(i)=cm(i)*(1.0_r8-cc(i,3))
                ENDDO
             ELSE
                ! cloudy portion
                DO i=1,m
                   ct(i)=cm(i)*cc(i,3)
                ENDDO
             ENDIF

             ! add one layer at a time, going down.
       !      DO k=MINVAL(icb),np
             DO k=1,np
                DO i=1,m
                   IF(k >= icb(i))then
                   denm = ts(i,k,is)/( 1.0_r8-rsa(i,k-1,ih,im)*rs(i,k,is) )
                   tda(i,k,ih,im)= tda(i,k-1,ih,im)*td(i,k,is)
                   tta(i,k,ih,im)=  tda(i,k-1,ih,im)*tt(i,k,is)  &
                        +(tda(i,k-1,ih,im)*rr(i,k,is)  &
                        *rsa(i,k-1,ih,im)+tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
                   rsa(i,k,ih,im)= rs(i,k,is)+ts(i,k,is)  &
                        *rsa(i,k-1,ih,im)*denm
                   endif
                ENDDO
             ENDDO

             ! add one layer at a time, going up.
          !   DO k=MINVAL(ict)-1,1,-1
             DO k=np,1,-1
                DO i=1,m
                   IF(k<=ict(i)-1)THEN 
                   denm =ts(i,k,ih)/(1.0_r8-rs(i,k,ih)*rxa(i,k+1,im,is))
                   rra(i,k,im,is)= rr(i,k,ih)+(td(i,k,ih)*rra(i,k+1,im,is) &
                        +(tt(i,k,ih)-td(i,k,ih))*rxa(i,k+1,im,is))*denm
                   rxa(i,k,im,is)= rs(i,k,ih)+ts(i,k,ih) &
                        *rxa(i,k+1,im,is)*denm
                   END IF
                ENDDO
             ENDDO

             ! compute fluxes following eq. (6.15) for fupdif and
             !     eq. (6.16) for (fdndir+fdndif)
             DO k=2,np+1
                DO i=1,m
                   denm= 1.0_r8/(1.0_r8-rsa(i,k-1,ih,im)*rxa(i,k,im,is))
                   fdndir(i)= tda(i,k-1,ih,im)
                   xx= tda(i,k-1,ih,im)*rra(i,k,im,is)
                   yy= tta(i,k-1,ih,im)-tda(i,k-1,ih,im)
                   fdndif(i)= (xx*rsa(i,k-1,ih,im)+yy)*denm

                   ! calculation of downward (fdn(i,k)) and upward (fupdif(i,k)) fluxes
                   ! as well as net flux
                   fupdif= (xx+yy*rxa(i,k,im,is))*denm
                   flxdn(i,k)= fdndir(i)+fdndif(i)-fupdif  ! net flux
                   fup(i,k)= (xx+yy*rxa(i,k,im,is))*denm   !new
                   fdn(i,k)= fdndir(i)+fdndif(i)  !new
                ENDDO
             ENDDO

             DO i=1,m
                flxdn(i,1)=1.0_r8-rra(i,1,im,is)
                fdn(i,1)=1.0_r8
                fup(i,1)=rra(i,1,im,is)
             ENDDO

             ! summation of fluxes over all sky situations;
             !     the term in the brackets of eq. (7.11)
             DO k=1,np+1
                DO i=1,m
                   IF(ih.eq.1 .and. im.eq.1 .and. is.eq.1) THEN
                      fclr(i,k)=flxdn(i,k)
                      fclr_d(i,k)=fdn(i,k)
                      fclr_u(i,k)=fup(i,k)
                   ENDIF
                   fall(i,k)=fall(i,k)+flxdn(i,k)*ct(i)
                   fall_d(i,k)=fall_d(i,k)+fdn(i,k)*ct(i)
                   fall_u(i,k)=fall_u(i,k)+fup(i,k)*ct(i)
                ENDDO
             ENDDO

             DO i=1,m
                fsdir(i)=fsdir(i)+fdndir(i)*ct(i)
                fsdif(i)=fsdif(i)+fdndif(i)*ct(i)
                !hmjb
                IF(ih.eq.1 .and. im.eq.1 .and. is.eq.1) THEN
                   fsdir_c(i)=fdndir(i)
                   fsdif_c(i)=fdndif(i)
                ENDIF
             ENDDO

          ENDDO ! end is loop
       ENDDO        ! end im loop
    ENDDO   ! end ih loop

  END SUBROUTINE cldflx






! ----------------------------------------------------------------------------
      subroutine cldprop_sw( &             !             cldprop_sw(
                            nlayers     , &!          nlayers   , &
                            cldfrac     , &!          cldfrac   (1:nlayers)                 , &
                            tauc        , &!          tauc_sw   (1:nbndsw,iplon,1:nlayers)  , &
                            ssac        , &!          ssac_sw   (1:nbndsw,iplon,1:nlayers)  , &
                            asmc        , &!          asmc_sw   (1:nbndsw,iplon,1:nlayers)  , &
                            ciwp        , &!          cicewp    (iplon,1:nlayers)               , &
                            clwp        , &!          cliqwp    (iplon,1:nlayers)               , &
                            rei         , &!          reicmcl   (iplon,1:nlayers)               , &
                            dge         , &!          2*reicmcl (iplon,1:nlayers)               , &
                            rel         , &!          relqmcl   (iplon,1:nlayers)               , &
                            taucldorig  , &!          taucldorig(1:nlayers,1:jpband)          , &
                            taucloud      )!          taucloud  (1:nlayers,1:jpband)          , &
! ----------------------------------------------------------------------------

! Purpose: Compute the cloud optical properties for each cloudy layer.
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=2,3 are available;
! (Hu & Stamnes, Key, and Fu) are implemented.

! ------- Input -------

      integer, intent(in) :: nlayers         ! total number of layers

      real(kind=r8), intent(in) :: cldfrac(1:nlayers)           ! cloud fraction
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: ciwp(1:nlayers)            ! cloud ice water path
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: clwp(1:nlayers)            ! cloud liquid water path
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: rei(1:nlayers)             ! cloud ice particle effective radius (microns)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: dge(1:nlayers)             ! cloud ice particle generalized effective size (microns)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: rel(1:nlayers)             ! cloud liquid particle effective radius (microns)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: tauc(1:nbndsw,1:nlayers)          ! cloud optical depth
                                                        !    Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: ssac(1:nbndsw,1:nlayers)          ! single scattering albedo
                                                        !    Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: asmc(1:nbndsw,1:nlayers)          ! asymmetry parameter
                                                        !    Dimensions: (nbndsw,nlayers)

! ------- Output -------

      real(kind=r8), intent(out) :: taucloud(1:nlayers,1:jpband)      ! cloud optical depth (delta scaled)
                                                        !    Dimensions: (nlayers,jpband)
      real(kind=r8), intent(out) :: taucldorig(1:nlayers,1:jpband)    ! cloud optical depth (non-delta scaled)
                                                        !    Dimensions: (nlayers,jpband)
      real(kind=r8) :: ssacloud(1:nlayers,1:jpband)     ! single scattering albedo (delta scaled)
                                                        !    Dimensions: (nlayers,jpband)
      real(kind=r8) :: asmcloud(1:nlayers,1:jpband)      ! asymmetry parameter (delta scaled)
                                                        !    Dimensions: (nlayers,jpband)

! ------- Local -------
      integer, PARAMETER :: inflag =0         ! see definitions
      integer, PARAMETER :: iceflag  =0       ! see definitions
      integer, PARAMETER :: liqflag    =0     ! see definitions

!      integer :: ncbands
      integer :: ib, ib1, ib2, lay, istr, index, icx

      real(kind=r8), parameter :: eps = 1.e-06_r8     ! epsilon
      real(kind=r8), parameter :: cldmin = 1.e-80_r8  ! minimum value for cloud quantities
      real(kind=r8) :: cwp                            ! total cloud water path
      real(kind=r8) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=r8) :: radice                         ! cloud ice effective radius (microns)
      real(kind=r8) :: dgeice                         ! cloud ice generalized effective size (microns)
      real(kind=r8) :: factor
      real(kind=r8) :: fint
      real(kind=r8) :: tauctot(nlayers)               ! band integrated cloud optical depth

      real(kind=r8) :: taucldorig_a, ssacloud_a, taucloud_a, ffp, ffp1, ffpssa
      real(kind=r8) :: tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq

      real(kind=r8) :: fdelta(jpb1:jpb2)
      real(kind=r8) :: extcoice(jpb1:jpb2), gice(jpb1:jpb2)
      real(kind=r8) :: ssacoice(jpb1:jpb2), forwice(jpb1:jpb2)
      real(kind=r8) :: extcoliq(jpb1:jpb2), gliq(jpb1:jpb2)
      real(kind=r8) :: ssacoliq(jpb1:jpb2), forwliq(jpb1:jpb2)

! Initialize

!      hvrcld = '$Revision: 1.4 $'
      taucloud=0.0_r8
      taucldorig=0.0_r8
      ssacloud=0.0_r8
      asmcloud=0.0_r8
!      ncbands = 29
      ib1 = jpb1
      ib2 = jpb2
      tauctot(:) = 0._r8

      do lay = 1, nlayers
         do ib = ib1 , ib2
            taucloud(lay,ib) = 0.0_r8
            ssacloud(lay,ib) = 1.0_r8
            asmcloud(lay,ib) = 0.0_r8
            tauctot(lay) = tauctot(lay) + tauc(ib-15,lay)
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers

         cwp = ciwp(lay) + clwp(lay)
         if (cldfrac(lay) .ge. cldmin .and. &
            (cwp .ge. cldmin .or. tauctot(lay) .ge. cldmin)) then

! (inflag=0): Cloud optical properties input directly
! Cloud optical properties already defined in tauc, ssac, asmc are unscaled;
! Apply delta-M scaling here
            if (inflag .eq. 0) then

               do ib = ib1 , ib2
                  taucldorig_a = tauc(ib-15,lay)
                  ffp = asmc(ib-15,lay)**2
                  ffp1 = 1.0_r8 - ffp
                  ffpssa = 1.0_r8 - ffp * ssac(ib-15,lay)
                  ssacloud_a = ffp1 * ssac(ib-15,lay) / ffpssa
                  taucloud_a = ffpssa * taucldorig_a

                  taucldorig(lay,ib) = taucldorig_a
                  ssacloud(lay,ib) = ssacloud_a
                  taucloud(lay,ib) = taucloud_a
                  asmcloud(lay,ib) = (asmc(ib-15,lay) - ffp) / (ffp1)
               enddo

! (inflag=2): Separate treatement of ice clouds and water clouds.
            elseif (inflag .eq. 2) then       
               radice = rei(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwp(lay) .eq. 0.0_r8) then
                  do ib = ib1 , ib2
                     extcoice(ib) = 0.0_r8
                     ssacoice(ib) = 0.0_r8
                     gice(ib)     = 0.0_r8
                     forwice(ib)  = 0.0_r8
                  enddo

! (iceflag = 1): 
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat ineffective for large ice particles
               elseif (iceflag .eq. 1) then
                  do ib = ib1, ib2
                     if (wavenum2(ib) .gt. 1.43e04_r8) then
                        icx = 1
                     elseif (wavenum2(ib) .gt. 7.7e03_r8) then
                        icx = 2
                     elseif (wavenum2(ib) .gt. 5.3e03_r8) then
                        icx = 3
                     elseif (wavenum2(ib) .gt. 4.0e03_r8) then
                        icx = 4
                     elseif (wavenum2(ib) .ge. 2.5e03_r8) then
                        icx = 5
                     endif
                     extcoice(ib) = abari(icx) + bbari(icx)/radice
                     ssacoice(ib) = 1._r8 - cbari(icx) - dbari(icx) * radice
                     gice(ib) = ebari(icx) + fbari(icx) * radice

! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice(ib) .ge. 1.0_r8) gice(ib) = 1.0_r8 - eps
                     forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0_r8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0_r8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0_r8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0_r8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0_r8) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

! For iceflag=2 option, combine with iceflag=0 option to handle large particle sizes.
! Use iceflag=2 option for ice particle effective radii from 5.0 to 131.0 microns
! and use iceflag=0 option for ice particles greater than 131.0 microns.
! *** NOTE: Transition between two methods has not been smoothed.

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0_r8) stop 'ICE RADIUS OUT OF BOUNDS'
                  if (radice .ge. 5.0_r8 .and. radice .le. 131._r8) then
                     factor = (radice - 2._r8)/3._r8
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - float(index)
                     do ib = ib1, ib2
                        extcoice(ib) = extice2(index,ib) + fint * &
                                      (extice2(index+1,ib) -  extice2(index,ib))
                        ssacoice(ib) = ssaice2(index,ib) + fint * &
                                      (ssaice2(index+1,ib) -  ssaice2(index,ib))
                        gice(ib) = asyice2(index,ib) + fint * &
                                      (asyice2(index+1,ib) -  asyice2(index,ib))
                        forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ib) .lt. 0.0_r8) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ib) .gt. 1.0_r8) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ib) .lt. 0.0_r8) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ib) .gt. 1.0_r8) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ib) .lt. 0.0_r8) stop 'ICE ASYM LESS THAN 0.0'
                     enddo
                  elseif (radice .gt. 131._r8) then
                     do ib = ib1, ib2
                        if (wavenum2(ib) .gt. 1.43e04_r8) then
                           icx = 1
                        elseif (wavenum2(ib) .gt. 7.7e03_r8) then
                           icx = 2
                        elseif (wavenum2(ib) .gt. 5.3e03_r8) then
                           icx = 3
                        elseif (wavenum2(ib) .gt. 4.0e03_r8) then
                           icx = 4
                        elseif (wavenum2(ib) .ge. 2.5e03_r8) then
                           icx = 5
                        endif
                        extcoice(ib) = abari(icx) + bbari(icx) / radice
                        ssacoice(ib) = 1._r8 - cbari(icx) - dbari(icx) * radice
                        gice(ib) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                        if (gice(ib) .ge. 1.0_r8) gice(ib) = 1.0_r8 - eps
                        forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ib) .lt. 0.0_r8) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ib) .gt. 1.0_r8) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ib) .lt. 0.0_r8) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ib) .gt. 1.0_r8) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ib) .lt. 0.0_r8) stop 'ICE ASYM LESS THAN 0.0'
                     enddo
                  endif

! For iceflag=3 option, combine with iceflag=0 option to handle large particle sizes
! Use iceflag=3 option for ice particle effective radii from 3.2 to 91.0 microns
! (generalized effective size, dge, from 5 to 140 microns), and use iceflag=0 option
! for ice particle effective radii greater than 91.0 microns (dge = 140 microns).
! *** NOTE: Fu parameterization requires particle size in generalized effective size.
! *** NOTE: Transition between two methods has not been smoothed. 

               elseif (iceflag .eq. 3) then
                  dgeice = dge(lay)
                  if (dgeice .lt. 5.0_r8) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                  if (dgeice .ge. 5.0_r8 .and. dgeice .le. 140._r8) then
                     factor = (dgeice - 2._r8)/3._r8
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - float(index)

                     do ib = ib1 , ib2
                        extcoice(ib) = extice3(index,ib) + fint * &
                                      (extice3(index+1,ib) - extice3(index,ib))
                        ssacoice(ib) = ssaice3(index,ib) + fint * &
                                      (ssaice3(index+1,ib) - ssaice3(index,ib))
                        gice(ib) = asyice3(index,ib) + fint * &
                                  (asyice3(index+1,ib) - asyice3(index,ib))
                        fdelta(ib) = fdlice3(index,ib) + fint * &
                                    (fdlice3(index+1,ib) - fdlice3(index,ib))
                        if (fdelta(ib) .lt. 0.0_r8) stop 'FDELTA LESS THAN 0.0'
                        if (fdelta(ib) .gt. 1.0_r8) stop 'FDELTA GT THAN 1.0'                     
                        forwice(ib) = fdelta(ib) + 0.5_r8 / ssacoice(ib)
! See Fu 1996 p. 2067 
                        if (forwice(ib) .gt. gice(ib)) forwice(ib) = gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ib) .lt. 0.0_r8) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ib) .gt. 1.0_r8) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ib) .lt. 0.0_r8) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ib) .gt. 1.0_r8) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ib) .lt. 0.0_r8) stop 'ICE ASYM LESS THAN 0.0'
                     enddo
                  elseif (dgeice .gt. 140._r8) then
                     do ib = ib1, ib2
                        if (wavenum2(ib) .gt. 1.43e04_r8) then
                           icx = 1
                        elseif (wavenum2(ib) .gt. 7.7e03_r8) then
                           icx = 2
                        elseif (wavenum2(ib) .gt. 5.3e03_r8) then
                           icx = 3
                        elseif (wavenum2(ib) .gt. 4.0e03_r8) then
                           icx = 4
                        elseif (wavenum2(ib) .ge. 2.5e03_r8) then
                           icx = 5
                        endif
                        extcoice(ib) = abari(icx) + bbari(icx)/radice
                        ssacoice(ib) = 1._r8 - cbari(icx) - dbari(icx) * radice
                        gice(ib) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                        if (gice(ib).ge.1.0_r8) gice(ib) = 1.0_r8-eps
                        forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ib) .lt. 0.0_r8) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ib) .gt. 1.0_r8) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ib) .lt. 0.0_r8) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ib) .gt. 1.0_r8) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ib) .lt. 0.0_r8) stop 'ICE ASYM LESS THAN 0.0'
                     enddo
                  endif
                endif
                  
! Calculation of absorption coefficients due to water clouds.
                if (clwp(lay) .eq. 0.0_r8) then
                   do ib = ib1 , ib2
                      extcoliq(ib) = 0.0_r8
                      ssacoliq(ib) = 0.0_r8
                      gliq(ib) = 0.0_r8
                      forwliq(ib) = 0.0_r8
                   enddo

                elseif (liqflag .eq. 1) then
                   radliq = rel(lay)
                   if (radliq .lt. 1.5_r8 .or. radliq .gt. 60._r8) stop &
                      'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                   index = int(radliq - 1.5_r8)
                   if (index .eq. 0) index = 1
                   if (index .eq. 58) index = 57
                   fint = radliq - 1.5_r8 - float(index)

                   do ib = ib1 , ib2
                      extcoliq(ib) = extliq1(index,ib) + fint * &
                                    (extliq1(index+1,ib) - extliq1(index,ib))
                      ssacoliq(ib) = ssaliq1(index,ib) + fint * &
                                    (ssaliq1(index+1,ib) - ssaliq1(index,ib))
                      if (fint .lt. 0._r8 .and. ssacoliq(ib) .gt. 1._r8) &
                                     ssacoliq(ib) = ssaliq1(index,ib)
                      gliq(ib) = asyliq1(index,ib) + fint * &
                                (asyliq1(index+1,ib) - asyliq1(index,ib))
                      forwliq(ib) = gliq(ib)*gliq(ib)
! Check to ensure all calculated quantities are within physical limits.
                      if (extcoliq(ib) .lt. 0.0_r8) stop 'LIQUID EXTINCTION LESS THAN 0.0'
                      if (ssacoliq(ib) .gt. 1.0_r8) stop 'LIQUID SSA GRTR THAN 1.0'
                      if (ssacoliq(ib) .lt. 0.0_r8) stop 'LIQUID SSA LESS THAN 0.0'
                      if (gliq(ib) .gt. 1.0_r8) stop 'LIQUID ASYM GRTR THAN 1.0'
                      if (gliq(ib) .lt. 0.0_r8) stop 'LIQUID ASYM LESS THAN 0.0'
                   enddo
                endif

                do ib = ib1 , ib2
                   tauliqorig = clwp(lay) * extcoliq(ib)
                   tauiceorig = ciwp(lay) * extcoice(ib)
                   taucldorig(lay,ib) = tauliqorig + tauiceorig

                   ssaliq = ssacoliq(ib) * (1.0_r8 - forwliq(ib)) / &
                           (1.0_r8 - forwliq(ib) * ssacoliq(ib))
                   tauliq = (1.0_r8 - forwliq(ib) * ssacoliq(ib)) * tauliqorig
                   ssaice = ssacoice(ib) * (1.0_r8 - forwice(ib)) / &
                           (1.0_r8 - forwice(ib) * ssacoice(ib))
                   tauice = (1.0_r8 - forwice(ib) * ssacoice(ib)) * tauiceorig

                   scatliq = ssaliq * tauliq
                   scatice = ssaice * tauice

                   taucloud(lay,ib) = tauliq + tauice

! Ensure non-zero taucmc and scatice
                   if (taucloud(lay,ib).eq.0.0_r8) taucloud(lay,ib) = cldmin
                   if (scatice.eq.0.0_r8) scatice = cldmin

                   ssacloud(lay,ib) = (scatliq + scatice) / taucloud(lay,ib)

                   if (iceflag .eq. 3) then
! In accordance with the 1996 Fu paper, equation A.3, 
! the moments for ice were calculated depending on whether using spheres
! or hexagonal ice crystals.
                      istr = 1
                      asmcloud(lay,ib) = (1.0_r8/(scatliq+scatice)) * &
                         (scatliq*(gliq(ib)**istr - forwliq(ib)) / &
                         (1.0_r8 - forwliq(ib)) + scatice * ((gice(ib)-forwice(ib)) / &
                         (1.0_r8 - forwice(ib)))**istr)
                   else 
! This code is the standard method for delta-m scaling. 
                      istr = 1
                      asmcloud(lay,ib) = (scatliq *  &
                         (gliq(ib)**istr - forwliq(ib)) / &
                         (1.0_r8 - forwliq(ib)) + scatice * (gice(ib)**istr - forwice(ib)) / &
                         (1.0_r8 - forwice(ib)))/(scatliq + scatice)
                   endif 

                enddo

            endif

         endif

! End layer loop
      enddo

      end subroutine cldprop_sw

 
   SUBROUTINE read_table()
    IMPLICIT NONE
    REAL(KIND=r8) :: data2(nm*nt*na)
    REAL(KIND=r8) :: data1(nt*na)
    INTEGER :: i,j,k,it 

    data2(1:nm*nt*na)=RESHAPE(SOURCE=(/&
     0.000_r8,0.068_r8,0.140_r8,0.216_r8,0.298_r8,0.385_r8,0.481_r8,0.586_r8,0.705_r8,0.840_r8,1.000_r8, &
     0.000_r8,0.052_r8,0.106_r8,0.166_r8,0.230_r8,0.302_r8,0.383_r8,0.478_r8,0.595_r8,0.752_r8,1.000_r8, &
     0.000_r8,0.038_r8,0.078_r8,0.120_r8,0.166_r8,0.218_r8,0.276_r8,0.346_r8,0.438_r8,0.582_r8,1.000_r8, &
     0.000_r8,0.030_r8,0.060_r8,0.092_r8,0.126_r8,0.164_r8,0.206_r8,0.255_r8,0.322_r8,0.442_r8,1.000_r8, &
     0.000_r8,0.025_r8,0.051_r8,0.078_r8,0.106_r8,0.136_r8,0.170_r8,0.209_r8,0.266_r8,0.462_r8,1.000_r8, &
     0.000_r8,0.023_r8,0.046_r8,0.070_r8,0.095_r8,0.122_r8,0.150_r8,0.187_r8,0.278_r8,0.577_r8,1.000_r8, &
     0.000_r8,0.022_r8,0.043_r8,0.066_r8,0.089_r8,0.114_r8,0.141_r8,0.187_r8,0.354_r8,0.603_r8,1.000_r8, &
     0.000_r8,0.021_r8,0.042_r8,0.063_r8,0.086_r8,0.108_r8,0.135_r8,0.214_r8,0.349_r8,0.565_r8,1.000_r8, &
     0.000_r8,0.021_r8,0.041_r8,0.062_r8,0.083_r8,0.105_r8,0.134_r8,0.202_r8,0.302_r8,0.479_r8,1.000_r8, &
     0.000_r8,0.088_r8,0.179_r8,0.272_r8,0.367_r8,0.465_r8,0.566_r8,0.669_r8,0.776_r8,0.886_r8,1.000_r8, &
     0.000_r8,0.079_r8,0.161_r8,0.247_r8,0.337_r8,0.431_r8,0.531_r8,0.637_r8,0.749_r8,0.870_r8,1.000_r8, &
     0.000_r8,0.065_r8,0.134_r8,0.207_r8,0.286_r8,0.372_r8,0.466_r8,0.572_r8,0.692_r8,0.831_r8,1.000_r8, &
     0.000_r8,0.049_r8,0.102_r8,0.158_r8,0.221_r8,0.290_r8,0.370_r8,0.465_r8,0.583_r8,0.745_r8,1.000_r8, &
     0.000_r8,0.037_r8,0.076_r8,0.118_r8,0.165_r8,0.217_r8,0.278_r8,0.354_r8,0.459_r8,0.638_r8,1.000_r8, &
     0.000_r8,0.030_r8,0.061_r8,0.094_r8,0.130_r8,0.171_r8,0.221_r8,0.286_r8,0.398_r8,0.631_r8,1.000_r8, &
     0.000_r8,0.026_r8,0.052_r8,0.081_r8,0.111_r8,0.146_r8,0.189_r8,0.259_r8,0.407_r8,0.643_r8,1.000_r8, &
     0.000_r8,0.023_r8,0.047_r8,0.072_r8,0.098_r8,0.129_r8,0.170_r8,0.250_r8,0.387_r8,0.598_r8,1.000_r8, &
     0.000_r8,0.022_r8,0.044_r8,0.066_r8,0.090_r8,0.118_r8,0.156_r8,0.224_r8,0.328_r8,0.508_r8,1.000_r8, &
     0.000_r8,0.094_r8,0.189_r8,0.285_r8,0.383_r8,0.482_r8,0.582_r8,0.685_r8,0.788_r8,0.894_r8,1.000_r8, &
     0.000_r8,0.088_r8,0.178_r8,0.271_r8,0.366_r8,0.465_r8,0.565_r8,0.669_r8,0.776_r8,0.886_r8,1.000_r8, &
     0.000_r8,0.079_r8,0.161_r8,0.247_r8,0.337_r8,0.431_r8,0.531_r8,0.637_r8,0.750_r8,0.870_r8,1.000_r8, &
     0.000_r8,0.066_r8,0.134_r8,0.209_r8,0.289_r8,0.375_r8,0.470_r8,0.577_r8,0.697_r8,0.835_r8,1.000_r8, &
     0.000_r8,0.050_r8,0.104_r8,0.163_r8,0.227_r8,0.300_r8,0.383_r8,0.483_r8,0.606_r8,0.770_r8,1.000_r8, &
     0.000_r8,0.038_r8,0.080_r8,0.125_r8,0.175_r8,0.233_r8,0.302_r8,0.391_r8,0.518_r8,0.710_r8,1.000_r8, &
     0.000_r8,0.031_r8,0.064_r8,0.100_r8,0.141_r8,0.188_r8,0.249_r8,0.336_r8,0.476_r8,0.689_r8,1.000_r8, &
     0.000_r8,0.026_r8,0.054_r8,0.084_r8,0.118_r8,0.158_r8,0.213_r8,0.298_r8,0.433_r8,0.638_r8,1.000_r8, &
     0.000_r8,0.023_r8,0.048_r8,0.074_r8,0.102_r8,0.136_r8,0.182_r8,0.254_r8,0.360_r8,0.542_r8,1.000_r8, &
     0.000_r8,0.096_r8,0.193_r8,0.290_r8,0.389_r8,0.488_r8,0.589_r8,0.690_r8,0.792_r8,0.896_r8,1.000_r8, &
     0.000_r8,0.092_r8,0.186_r8,0.281_r8,0.378_r8,0.477_r8,0.578_r8,0.680_r8,0.785_r8,0.891_r8,1.000_r8, &
     0.000_r8,0.086_r8,0.174_r8,0.264_r8,0.358_r8,0.455_r8,0.556_r8,0.660_r8,0.769_r8,0.882_r8,1.000_r8, &
     0.000_r8,0.074_r8,0.153_r8,0.235_r8,0.323_r8,0.416_r8,0.514_r8,0.622_r8,0.737_r8,0.862_r8,1.000_r8, &
     0.000_r8,0.061_r8,0.126_r8,0.195_r8,0.271_r8,0.355_r8,0.449_r8,0.555_r8,0.678_r8,0.823_r8,1.000_r8, &
     0.000_r8,0.047_r8,0.098_r8,0.153_r8,0.215_r8,0.286_r8,0.370_r8,0.471_r8,0.600_r8,0.770_r8,1.000_r8, &
     0.000_r8,0.037_r8,0.077_r8,0.120_r8,0.170_r8,0.230_r8,0.303_r8,0.401_r8,0.537_r8,0.729_r8,1.000_r8, &
     0.000_r8,0.030_r8,0.062_r8,0.098_r8,0.138_r8,0.187_r8,0.252_r8,0.343_r8,0.476_r8,0.673_r8,1.000_r8, &
     0.000_r8,0.026_r8,0.053_r8,0.082_r8,0.114_r8,0.154_r8,0.207_r8,0.282_r8,0.391_r8,0.574_r8,1.000_r8, &
     0.000_r8,0.097_r8,0.194_r8,0.293_r8,0.392_r8,0.492_r8,0.592_r8,0.693_r8,0.794_r8,0.897_r8,1.000_r8, &
     0.000_r8,0.094_r8,0.190_r8,0.286_r8,0.384_r8,0.483_r8,0.584_r8,0.686_r8,0.789_r8,0.894_r8,1.000_r8, &
     0.000_r8,0.090_r8,0.181_r8,0.274_r8,0.370_r8,0.468_r8,0.569_r8,0.672_r8,0.778_r8,0.887_r8,1.000_r8, &
     0.000_r8,0.081_r8,0.165_r8,0.252_r8,0.343_r8,0.439_r8,0.539_r8,0.645_r8,0.757_r8,0.874_r8,1.000_r8, &
     0.000_r8,0.069_r8,0.142_r8,0.218_r8,0.302_r8,0.392_r8,0.490_r8,0.598_r8,0.717_r8,0.850_r8,1.000_r8, &
     0.000_r8,0.054_r8,0.114_r8,0.178_r8,0.250_r8,0.330_r8,0.422_r8,0.529_r8,0.656_r8,0.810_r8,1.000_r8, &
     0.000_r8,0.042_r8,0.090_r8,0.141_r8,0.200_r8,0.269_r8,0.351_r8,0.455_r8,0.589_r8,0.764_r8,1.000_r8, &
     0.000_r8,0.034_r8,0.070_r8,0.112_r8,0.159_r8,0.217_r8,0.289_r8,0.384_r8,0.515_r8,0.703_r8,1.000_r8, &
     0.000_r8,0.028_r8,0.058_r8,0.090_r8,0.128_r8,0.174_r8,0.231_r8,0.309_r8,0.420_r8,0.602_r8,1.000_r8, &
     0.000_r8,0.098_r8,0.196_r8,0.295_r8,0.394_r8,0.494_r8,0.594_r8,0.695_r8,0.796_r8,0.898_r8,1.000_r8, &
     0.000_r8,0.096_r8,0.193_r8,0.290_r8,0.389_r8,0.488_r8,0.588_r8,0.690_r8,0.792_r8,0.895_r8,1.000_r8, &
     0.000_r8,0.092_r8,0.186_r8,0.281_r8,0.378_r8,0.477_r8,0.577_r8,0.680_r8,0.784_r8,0.891_r8,1.000_r8, &
     0.000_r8,0.086_r8,0.174_r8,0.264_r8,0.358_r8,0.455_r8,0.556_r8,0.661_r8,0.769_r8,0.882_r8,1.000_r8, &
     0.000_r8,0.075_r8,0.154_r8,0.237_r8,0.325_r8,0.419_r8,0.518_r8,0.626_r8,0.741_r8,0.865_r8,1.000_r8, &
     0.000_r8,0.062_r8,0.129_r8,0.201_r8,0.279_r8,0.366_r8,0.462_r8,0.571_r8,0.694_r8,0.836_r8,1.000_r8, &
     0.000_r8,0.049_r8,0.102_r8,0.162_r8,0.229_r8,0.305_r8,0.394_r8,0.501_r8,0.631_r8,0.793_r8,1.000_r8, &
     0.000_r8,0.038_r8,0.080_r8,0.127_r8,0.182_r8,0.245_r8,0.323_r8,0.422_r8,0.550_r8,0.730_r8,1.000_r8, &
     0.000_r8,0.030_r8,0.064_r8,0.100_r8,0.142_r8,0.192_r8,0.254_r8,0.334_r8,0.448_r8,0.627_r8,1.000_r8, &
     0.000_r8,0.098_r8,0.198_r8,0.296_r8,0.396_r8,0.496_r8,0.596_r8,0.696_r8,0.797_r8,0.898_r8,1.000_r8, &
     0.000_r8,0.097_r8,0.194_r8,0.293_r8,0.392_r8,0.491_r8,0.591_r8,0.693_r8,0.794_r8,0.897_r8,1.000_r8, &
     0.000_r8,0.094_r8,0.190_r8,0.286_r8,0.384_r8,0.483_r8,0.583_r8,0.686_r8,0.789_r8,0.894_r8,1.000_r8, &
     0.000_r8,0.089_r8,0.180_r8,0.274_r8,0.369_r8,0.467_r8,0.568_r8,0.672_r8,0.778_r8,0.887_r8,1.000_r8, &
     0.000_r8,0.081_r8,0.165_r8,0.252_r8,0.344_r8,0.440_r8,0.541_r8,0.646_r8,0.758_r8,0.875_r8,1.000_r8, &
     0.000_r8,0.069_r8,0.142_r8,0.221_r8,0.306_r8,0.397_r8,0.496_r8,0.604_r8,0.722_r8,0.854_r8,1.000_r8, &
     0.000_r8,0.056_r8,0.116_r8,0.182_r8,0.256_r8,0.338_r8,0.432_r8,0.540_r8,0.666_r8,0.816_r8,1.000_r8, &
     0.000_r8,0.043_r8,0.090_r8,0.143_r8,0.203_r8,0.273_r8,0.355_r8,0.455_r8,0.583_r8,0.754_r8,1.000_r8, &
     0.000_r8,0.034_r8,0.070_r8,0.111_r8,0.157_r8,0.210_r8,0.276_r8,0.359_r8,0.474_r8,0.650_r8,1.000_r8, &
     0.000_r8,0.099_r8,0.198_r8,0.298_r8,0.398_r8,0.497_r8,0.598_r8,0.698_r8,0.798_r8,0.899_r8,1.000_r8, &
     0.000_r8,0.098_r8,0.196_r8,0.295_r8,0.394_r8,0.494_r8,0.594_r8,0.695_r8,0.796_r8,0.898_r8,1.000_r8, &
     0.000_r8,0.096_r8,0.193_r8,0.290_r8,0.390_r8,0.489_r8,0.589_r8,0.690_r8,0.793_r8,0.896_r8,1.000_r8, &
     0.000_r8,0.093_r8,0.186_r8,0.282_r8,0.379_r8,0.478_r8,0.578_r8,0.681_r8,0.786_r8,0.892_r8,1.000_r8, &
     0.000_r8,0.086_r8,0.175_r8,0.266_r8,0.361_r8,0.458_r8,0.558_r8,0.663_r8,0.771_r8,0.883_r8,1.000_r8, &
     0.000_r8,0.076_r8,0.156_r8,0.240_r8,0.330_r8,0.423_r8,0.523_r8,0.630_r8,0.744_r8,0.867_r8,1.000_r8, &
     0.000_r8,0.063_r8,0.130_r8,0.203_r8,0.282_r8,0.369_r8,0.465_r8,0.572_r8,0.694_r8,0.834_r8,1.000_r8, &
     0.000_r8,0.049_r8,0.102_r8,0.161_r8,0.226_r8,0.299_r8,0.385_r8,0.486_r8,0.611_r8,0.774_r8,1.000_r8, &
     0.000_r8,0.038_r8,0.078_r8,0.122_r8,0.172_r8,0.229_r8,0.297_r8,0.382_r8,0.498_r8,0.672_r8,1.000_r8, &
     0.000_r8,0.099_r8,0.199_r8,0.298_r8,0.398_r8,0.498_r8,0.598_r8,0.699_r8,0.799_r8,0.899_r8,1.000_r8, &
     0.000_r8,0.099_r8,0.198_r8,0.298_r8,0.398_r8,0.497_r8,0.598_r8,0.698_r8,0.798_r8,0.899_r8,1.000_r8, &
     0.000_r8,0.098_r8,0.196_r8,0.295_r8,0.394_r8,0.494_r8,0.594_r8,0.695_r8,0.796_r8,0.898_r8,1.000_r8, &
     0.000_r8,0.096_r8,0.193_r8,0.290_r8,0.389_r8,0.488_r8,0.588_r8,0.690_r8,0.792_r8,0.895_r8,1.000_r8, &
     0.000_r8,0.092_r8,0.185_r8,0.280_r8,0.376_r8,0.474_r8,0.575_r8,0.678_r8,0.782_r8,0.890_r8,1.000_r8, &
     0.000_r8,0.084_r8,0.170_r8,0.259_r8,0.351_r8,0.447_r8,0.547_r8,0.652_r8,0.762_r8,0.878_r8,1.000_r8, &
     0.000_r8,0.071_r8,0.146_r8,0.224_r8,0.308_r8,0.398_r8,0.494_r8,0.601_r8,0.718_r8,0.850_r8,1.000_r8, &
     0.000_r8,0.056_r8,0.114_r8,0.178_r8,0.248_r8,0.325_r8,0.412_r8,0.514_r8,0.638_r8,0.793_r8,1.000_r8, &
     0.000_r8,0.042_r8,0.086_r8,0.134_r8,0.186_r8,0.246_r8,0.318_r8,0.405_r8,0.521_r8,0.691_r8,1.000_r8, &
     0.000_r8,0.100_r8,0.200_r8,0.300_r8,0.400_r8,0.500_r8,0.600_r8,0.700_r8,0.800_r8,0.900_r8,1.000_r8, &
     0.000_r8,0.100_r8,0.200_r8,0.300_r8,0.400_r8,0.500_r8,0.600_r8,0.700_r8,0.800_r8,0.900_r8,1.000_r8, &
     0.000_r8,0.100_r8,0.200_r8,0.300_r8,0.400_r8,0.500_r8,0.600_r8,0.700_r8,0.800_r8,0.900_r8,1.000_r8, &
     0.000_r8,0.100_r8,0.199_r8,0.298_r8,0.398_r8,0.498_r8,0.598_r8,0.698_r8,0.798_r8,0.899_r8,1.000_r8, &
     0.000_r8,0.098_r8,0.196_r8,0.294_r8,0.392_r8,0.491_r8,0.590_r8,0.691_r8,0.793_r8,0.896_r8,1.000_r8, &
     0.000_r8,0.092_r8,0.185_r8,0.278_r8,0.374_r8,0.470_r8,0.570_r8,0.671_r8,0.777_r8,0.886_r8,1.000_r8, &
     0.000_r8,0.081_r8,0.162_r8,0.246_r8,0.333_r8,0.424_r8,0.521_r8,0.625_r8,0.738_r8,0.862_r8,1.000_r8, &
     0.000_r8,0.063_r8,0.128_r8,0.196_r8,0.270_r8,0.349_r8,0.438_r8,0.540_r8,0.661_r8,0.809_r8,1.000_r8, &
     0.000_r8,0.046_r8,0.094_r8,0.146_r8,0.202_r8,0.264_r8,0.337_r8,0.426_r8,0.542_r8,0.710_r8,1.000_r8, &
     0.000_r8,0.101_r8,0.202_r8,0.302_r8,0.402_r8,0.502_r8,0.602_r8,0.702_r8,0.802_r8,0.901_r8,1.000_r8, &
     0.000_r8,0.102_r8,0.202_r8,0.303_r8,0.404_r8,0.504_r8,0.604_r8,0.703_r8,0.802_r8,0.902_r8,1.000_r8, &
     0.000_r8,0.102_r8,0.205_r8,0.306_r8,0.406_r8,0.506_r8,0.606_r8,0.706_r8,0.804_r8,0.902_r8,1.000_r8, &
     0.000_r8,0.104_r8,0.207_r8,0.309_r8,0.410_r8,0.510_r8,0.609_r8,0.707_r8,0.805_r8,0.902_r8,1.000_r8, &
     0.000_r8,0.106_r8,0.208_r8,0.309_r8,0.409_r8,0.508_r8,0.606_r8,0.705_r8,0.803_r8,0.902_r8,1.000_r8, &
     0.000_r8,0.102_r8,0.202_r8,0.298_r8,0.395_r8,0.493_r8,0.590_r8,0.690_r8,0.790_r8,0.894_r8,1.000_r8, &
     0.000_r8,0.091_r8,0.179_r8,0.267_r8,0.357_r8,0.449_r8,0.545_r8,0.647_r8,0.755_r8,0.872_r8,1.000_r8, &
     0.000_r8,0.073_r8,0.142_r8,0.214_r8,0.290_r8,0.372_r8,0.462_r8,0.563_r8,0.681_r8,0.822_r8,1.000_r8, &
     0.000_r8,0.053_r8,0.104_r8,0.158_r8,0.217_r8,0.281_r8,0.356_r8,0.446_r8,0.562_r8,0.726_r8,1.000_r8/),SHAPE=(/nm*nt*na/))
    it=0
    DO k=1,nm
       DO i=1,nt
          DO j=1,na
             it=it+1
             caib(k,i,j) =data2(it)
          END DO
       END DO
    END DO

    data1(1:nt*na)=RESHAPE(SOURCE=(/&
    0.000_r8,0.099_r8,0.198_r8,0.297_r8,0.397_r8,0.496_r8,0.597_r8,0.697_r8,0.798_r8,0.899_r8,1.000_r8, &
    0.000_r8,0.098_r8,0.196_r8,0.294_r8,0.394_r8,0.494_r8,0.594_r8,0.694_r8,0.796_r8,0.898_r8,1.000_r8, &
    0.000_r8,0.096_r8,0.192_r8,0.290_r8,0.388_r8,0.487_r8,0.587_r8,0.689_r8,0.792_r8,0.895_r8,1.000_r8, &
    0.000_r8,0.092_r8,0.185_r8,0.280_r8,0.376_r8,0.476_r8,0.576_r8,0.678_r8,0.783_r8,0.890_r8,1.000_r8, &
    0.000_r8,0.085_r8,0.173_r8,0.263_r8,0.357_r8,0.454_r8,0.555_r8,0.659_r8,0.768_r8,0.881_r8,1.000_r8, &
    0.000_r8,0.076_r8,0.154_r8,0.237_r8,0.324_r8,0.418_r8,0.517_r8,0.624_r8,0.738_r8,0.864_r8,1.000_r8, &
    0.000_r8,0.063_r8,0.131_r8,0.203_r8,0.281_r8,0.366_r8,0.461_r8,0.567_r8,0.688_r8,0.830_r8,1.000_r8, &
    0.000_r8,0.052_r8,0.107_r8,0.166_r8,0.232_r8,0.305_r8,0.389_r8,0.488_r8,0.610_r8,0.770_r8,1.000_r8, &
    0.000_r8,0.043_r8,0.088_r8,0.136_r8,0.189_r8,0.248_r8,0.317_r8,0.400_r8,0.510_r8,0.675_r8,1.000_r8/),SHAPE=(/nt*na/))
    it=0
    DO i=1,nt
       DO j=1,na
          it=it+1
          caif(i,j)=data1(it)
       END DO
    END DO

    coa=RESHAPE(SOURCE=(/&
       0.0000080_r8,  0.0000089_r8,  0.0000098_r8,  0.0000106_r8,  0.0000114_r8,  &
       0.0000121_r8,  0.0000128_r8,  0.0000134_r8,  0.0000140_r8,  0.0000146_r8,  &
       0.0000152_r8,  0.0000158_r8,  0.0000163_r8,  0.0000168_r8,  0.0000173_r8,  &
       0.0000178_r8,  0.0000182_r8,  0.0000186_r8,  0.0000191_r8,  0.0000195_r8,  &
       0.0000199_r8,  0.0000202_r8,  0.0000206_r8,  0.0000210_r8,  0.0000213_r8,  &
       0.0000217_r8,  0.0000220_r8,  0.0000223_r8,  0.0000226_r8,  0.0000229_r8,  &
       0.0000232_r8,  0.0000235_r8,  0.0000238_r8,  0.0000241_r8,  0.0000244_r8,  &
       0.0000246_r8,  0.0000249_r8,  0.0000252_r8,  0.0000254_r8,  0.0000257_r8,  &
       0.0000259_r8,  0.0000261_r8,  0.0000264_r8,  0.0000266_r8,  0.0000268_r8,  &
       0.0000271_r8,  0.0000273_r8,  0.0000275_r8,  0.0000277_r8,  0.0000279_r8,  &
       0.0000281_r8,  0.0000283_r8,  0.0000285_r8,  0.0000287_r8,  0.0000289_r8,  &
       0.0000291_r8,  0.0000293_r8,  0.0000295_r8,  0.0000297_r8,  0.0000298_r8,  &
       0.0000300_r8,  0.0000302_r8, &
       0.0000085_r8,  0.0000095_r8,  0.0000104_r8,  0.0000113_r8,  0.0000121_r8,  &
       0.0000128_r8,  0.0000136_r8,  0.0000143_r8,  0.0000149_r8,  0.0000155_r8,  &
       0.0000161_r8,  0.0000167_r8,  0.0000172_r8,  0.0000178_r8,  0.0000183_r8,  &
       0.0000187_r8,  0.0000192_r8,  0.0000196_r8,  0.0000201_r8,  0.0000205_r8,  &
       0.0000209_r8,  0.0000213_r8,  0.0000217_r8,  0.0000220_r8,  0.0000224_r8,  &
       0.0000227_r8,  0.0000231_r8,  0.0000234_r8,  0.0000237_r8,  0.0000240_r8,  &
       0.0000243_r8,  0.0000246_r8,  0.0000249_r8,  0.0000252_r8,  0.0000255_r8,  &
       0.0000258_r8,  0.0000260_r8,  0.0000263_r8,  0.0000266_r8,  0.0000268_r8,  &
       0.0000271_r8,  0.0000273_r8,  0.0000275_r8,  0.0000278_r8,  0.0000280_r8,  &
       0.0000282_r8,  0.0000285_r8,  0.0000287_r8,  0.0000289_r8,  0.0000291_r8,  &
       0.0000293_r8,  0.0000295_r8,  0.0000297_r8,  0.0000299_r8,  0.0000301_r8,  &
       0.0000303_r8,  0.0000305_r8,  0.0000307_r8,  0.0000309_r8,  0.0000311_r8,  &
       0.0000313_r8,  0.0000314_r8, &
       0.0000095_r8,  0.0000106_r8,  0.0000116_r8,  0.0000125_r8,  0.0000134_r8,  &
       0.0000143_r8,  0.0000150_r8,  0.0000158_r8,  0.0000165_r8,  0.0000171_r8,  &
       0.0000178_r8,  0.0000184_r8,  0.0000189_r8,  0.0000195_r8,  0.0000200_r8,  &
       0.0000205_r8,  0.0000210_r8,  0.0000215_r8,  0.0000219_r8,  0.0000223_r8,  &
       0.0000228_r8,  0.0000232_r8,  0.0000235_r8,  0.0000239_r8,  0.0000243_r8,  &
       0.0000247_r8,  0.0000250_r8,  0.0000253_r8,  0.0000257_r8,  0.0000260_r8,  &
       0.0000263_r8,  0.0000266_r8,  0.0000269_r8,  0.0000272_r8,  0.0000275_r8,  &
       0.0000278_r8,  0.0000281_r8,  0.0000283_r8,  0.0000286_r8,  0.0000289_r8,  &
       0.0000291_r8,  0.0000294_r8,  0.0000296_r8,  0.0000299_r8,  0.0000301_r8,  &
       0.0000303_r8,  0.0000306_r8,  0.0000308_r8,  0.0000310_r8,  0.0000312_r8,  &
       0.0000315_r8,  0.0000317_r8,  0.0000319_r8,  0.0000321_r8,  0.0000323_r8,  &
       0.0000325_r8,  0.0000327_r8,  0.0000329_r8,  0.0000331_r8,  0.0000333_r8,  &
       0.0000335_r8,  0.0000329_r8, &
       0.0000100_r8,  0.0000111_r8,  0.0000122_r8,  0.0000131_r8,  0.0000141_r8,  &
       0.0000149_r8,  0.0000157_r8,  0.0000165_r8,  0.0000172_r8,  0.0000179_r8,  &
       0.0000185_r8,  0.0000191_r8,  0.0000197_r8,  0.0000203_r8,  0.0000208_r8,  &
       0.0000213_r8,  0.0000218_r8,  0.0000223_r8,  0.0000227_r8,  0.0000232_r8,  &
       0.0000236_r8,  0.0000240_r8,  0.0000244_r8,  0.0000248_r8,  0.0000252_r8,  &
       0.0000255_r8,  0.0000259_r8,  0.0000262_r8,  0.0000266_r8,  0.0000269_r8,  &
       0.0000272_r8,  0.0000275_r8,  0.0000278_r8,  0.0000281_r8,  0.0000284_r8,  &
       0.0000287_r8,  0.0000290_r8,  0.0000293_r8,  0.0000295_r8,  0.0000298_r8,  &
       0.0000300_r8,  0.0000303_r8,  0.0000306_r8,  0.0000308_r8,  0.0000310_r8,  &
       0.0000313_r8,  0.0000315_r8,  0.0000317_r8,  0.0000320_r8,  0.0000322_r8,  &
       0.0000324_r8,  0.0000326_r8,  0.0000328_r8,  0.0000331_r8,  0.0000333_r8,  &
       0.0000335_r8,  0.0000330_r8,  0.0000339_r8,  0.0000341_r8,  0.0000343_r8,  &
       0.0000345_r8,  0.0000346_r8, &
       0.0000109_r8,  0.0000121_r8,  0.0000132_r8,  0.0000143_r8,  0.0000152_r8,  &
       0.0000161_r8,  0.0000170_r8,  0.0000178_r8,  0.0000185_r8,  0.0000192_r8,  &
       0.0000199_r8,  0.0000205_r8,  0.0000211_r8,  0.0000217_r8,  0.0000222_r8,  &
       0.0000228_r8,  0.0000233_r8,  0.0000238_r8,  0.0000242_r8,  0.0000247_r8,  &
       0.0000251_r8,  0.0000255_r8,  0.0000259_r8,  0.0000263_r8,  0.0000267_r8,  &
       0.0000271_r8,  0.0000275_r8,  0.0000278_r8,  0.0000282_r8,  0.0000285_r8,  &
       0.0000288_r8,  0.0000291_r8,  0.0000295_r8,  0.0000298_r8,  0.0000301_r8,  &
       0.0000304_r8,  0.0000307_r8,  0.0000309_r8,  0.0000312_r8,  0.0000315_r8,  &
       0.0000318_r8,  0.0000320_r8,  0.0000323_r8,  0.0000325_r8,  0.0000328_r8,  &
       0.0000330_r8,  0.0000333_r8,  0.0000335_r8,  0.0000330_r8,  0.0000340_r8,  &
       0.0000342_r8,  0.0000344_r8,  0.0000346_r8,  0.0000348_r8,  0.0000351_r8,  &
       0.0000353_r8,  0.0000355_r8,  0.0000357_r8,  0.0000359_r8,  0.0000361_r8,  &
       0.0000363_r8,  0.0000365_r8, &
       0.0000117_r8,  0.0000130_r8,  0.0000142_r8,  0.0000153_r8,  0.0000163_r8,  &
       0.0000173_r8,  0.0000181_r8,  0.0000190_r8,  0.0000197_r8,  0.0000204_r8,  &
       0.0000211_r8,  0.0000218_r8,  0.0000224_r8,  0.0000230_r8,  0.0000235_r8,  &
       0.0000241_r8,  0.0000246_r8,  0.0000251_r8,  0.0000256_r8,  0.0000260_r8,  &
       0.0000265_r8,  0.0000269_r8,  0.0000273_r8,  0.0000277_r8,  0.0000281_r8,  &
       0.0000285_r8,  0.0000289_r8,  0.0000293_r8,  0.0000296_r8,  0.0000299_r8,  &
       0.0000303_r8,  0.0000306_r8,  0.0000309_r8,  0.0000313_r8,  0.0000316_r8,  &
       0.0000319_r8,  0.0000322_r8,  0.0000324_r8,  0.0000327_r8,  0.0000330_r8,  &
       0.0000333_r8,  0.0000336_r8,  0.0000331_r8,  0.0000341_r8,  0.0000343_r8,  &
       0.0000346_r8,  0.0000348_r8,  0.0000351_r8,  0.0000353_r8,  0.0000355_r8,  &
       0.0000358_r8,  0.0000360_r8,  0.0000362_r8,  0.0000365_r8,  0.0000367_r8,  &
       0.0000369_r8,  0.0000371_r8,  0.0000373_r8,  0.0000375_r8,  0.0000377_r8,  &
       0.0000379_r8,  0.0000381_r8, &
       0.0000125_r8,  0.0000139_r8,  0.0000151_r8,  0.0000163_r8,  0.0000173_r8,  &
       0.0000183_r8,  0.0000192_r8,  0.0000200_r8,  0.0000208_r8,  0.0000216_r8,  &
       0.0000223_r8,  0.0000229_r8,  0.0000236_r8,  0.0000242_r8,  0.0000247_r8,  &
       0.0000253_r8,  0.0000258_r8,  0.0000263_r8,  0.0000268_r8,  0.0000273_r8,  &
       0.0000277_r8,  0.0000282_r8,  0.0000286_r8,  0.0000290_r8,  0.0000294_r8,  &       
       0.0000298_r8,  0.0000302_r8,  0.0000306_r8,  0.0000309_r8,  0.0000313_r8,  &
       0.0000316_r8,  0.0000320_r8,  0.0000323_r8,  0.0000326_r8,  0.0000329_r8,  &
       0.0000332_r8,  0.0000335_r8,  0.0000331_r8,  0.0000341_r8,  0.0000344_r8,  &
       0.0000347_r8,  0.0000350_r8,  0.0000352_r8,  0.0000355_r8,  0.0000358_r8,  &
       0.0000360_r8,  0.0000363_r8,  0.0000365_r8,  0.0000368_r8,  0.0000370_r8,  &
       0.0000372_r8,  0.0000375_r8,  0.0000377_r8,  0.0000379_r8,  0.0000382_r8,  &
       0.0000384_r8,  0.0000386_r8,  0.0000388_r8,  0.0000390_r8,  0.0000392_r8,  &
       0.0000394_r8,  0.0000396_r8, &
       0.0000132_r8,  0.0000147_r8,  0.0000160_r8,  0.0000172_r8,  0.0000183_r8,  &
       0.0000193_r8,  0.0000202_r8,  0.0000210_r8,  0.0000218_r8,  0.0000226_r8,  &
       0.0000233_r8,  0.0000240_r8,  0.0000246_r8,  0.0000252_r8,  0.0000258_r8,  &
       0.0000264_r8,  0.0000269_r8,  0.0000274_r8,  0.0000279_r8,  0.0000284_r8,  &
       0.0000289_r8,  0.0000293_r8,  0.0000298_r8,  0.0000302_r8,  0.0000306_r8,  &
       0.0000310_r8,  0.0000314_r8,  0.0000318_r8,  0.0000321_r8,  0.0000325_r8,  &
       0.0000328_r8,  0.0000332_r8,  0.0000335_r8,  0.0000331_r8,  0.0000342_r8,  &
       0.0000345_r8,  0.0000348_r8,  0.0000351_r8,  0.0000354_r8,  0.0000357_r8,  &
       0.0000360_r8,  0.0000363_r8,  0.0000365_r8,  0.0000368_r8,  0.0000371_r8,  &
       0.0000373_r8,  0.0000376_r8,  0.0000378_r8,  0.0000381_r8,  0.0000383_r8,  &
       0.0000386_r8,  0.0000388_r8,  0.0000391_r8,  0.0000393_r8,  0.0000395_r8,  &
       0.0000397_r8,  0.0000400_r8,  0.0000402_r8,  0.0000404_r8,  0.0000406_r8,  &
       0.0000408_r8,  0.0000411_r8, &
       0.0000143_r8,  0.0000158_r8,  0.0000172_r8,  0.0000184_r8,  0.0000195_r8,  &
       0.0000206_r8,  0.0000215_r8,  0.0000224_r8,  0.0000232_r8,  0.0000240_r8,  &
       0.0000247_r8,  0.0000254_r8,  0.0000261_r8,  0.0000267_r8,  0.0000273_r8,  &
       0.0000279_r8,  0.0000284_r8,  0.0000290_r8,  0.0000295_r8,  0.0000300_r8,  &
       0.0000305_r8,  0.0000309_r8,  0.0000314_r8,  0.0000318_r8,  0.0000322_r8,  &
       0.0000326_r8,  0.0000330_r8,  0.0000334_r8,  0.0000331_r8,  0.0000342_r8,  &
       0.0000345_r8,  0.0000349_r8,  0.0000352_r8,  0.0000356_r8,  0.0000359_r8,  &
       0.0000362_r8,  0.0000365_r8,  0.0000368_r8,  0.0000371_r8,  0.0000374_r8,  &
       0.0000377_r8,  0.0000380_r8,  0.0000383_r8,  0.0000386_r8,  0.0000389_r8,  &
       0.0000391_r8,  0.0000394_r8,  0.0000397_r8,  0.0000399_r8,  0.0000402_r8,  &
       0.0000404_r8,  0.0000407_r8,  0.0000409_r8,  0.0000412_r8,  0.0000414_r8,  &
       0.0000416_r8,  0.0000419_r8,  0.0000421_r8,  0.0000423_r8,  0.0000426_r8,  &
       0.0000428_r8,  0.0000430_r8, &
       0.0000153_r8,  0.0000169_r8,  0.0000183_r8,  0.0000196_r8,  0.0000207_r8,   &
       0.0000218_r8,  0.0000227_r8,  0.0000236_r8,  0.0000245_r8,  0.0000253_r8,   &
       0.0000260_r8,  0.0000267_r8,  0.0000274_r8,  0.0000281_r8,  0.0000287_r8,   &
       0.0000293_r8,  0.0000298_r8,  0.0000304_r8,  0.0000309_r8,  0.0000314_r8,   &
       0.0000319_r8,  0.0000324_r8,  0.0000328_r8,  0.0000333_r8,  0.0000330_r8,   &
       0.0000341_r8,  0.0000345_r8,  0.0000349_r8,  0.0000353_r8,  0.0000357_r8,   &
       0.0000361_r8,  0.0000364_r8,  0.0000368_r8,  0.0000371_r8,  0.0000375_r8,   &
       0.0000378_r8,  0.0000381_r8,  0.0000384_r8,  0.0000387_r8,  0.0000391_r8,   &
       0.0000394_r8,  0.0000397_r8,  0.0000399_r8,  0.0000402_r8,  0.0000405_r8,   &
       0.0000408_r8,  0.0000411_r8,  0.0000413_r8,  0.0000416_r8,  0.0000419_r8,   &
       0.0000421_r8,  0.0000424_r8,  0.0000426_r8,  0.0000429_r8,  0.0000431_r8,   &
       0.0000434_r8,  0.0000436_r8,  0.0000439_r8,  0.0000441_r8,  0.0000443_r8,   &
       0.0000446_r8,  0.0000448_r8, &
       0.0000165_r8,  0.0000182_r8,  0.0000196_r8,  0.0000209_r8,  0.0000221_r8,  &
       0.0000232_r8,  0.0000242_r8,  0.0000251_r8,  0.0000260_r8,  0.0000268_r8,  &
       0.0000276_r8,  0.0000283_r8,  0.0000290_r8,  0.0000297_r8,  0.0000303_r8,  &
       0.0000309_r8,  0.0000315_r8,  0.0000321_r8,  0.0000326_r8,  0.0000331_r8,  &
       0.0000336_r8,  0.0000341_r8,  0.0000346_r8,  0.0000350_r8,  0.0000355_r8,  &
       0.0000359_r8,  0.0000363_r8,  0.0000367_r8,  0.0000371_r8,  0.0000375_r8,  &
       0.0000379_r8,  0.0000383_r8,  0.0000386_r8,  0.0000390_r8,  0.0000394_r8,  &
       0.0000397_r8,  0.0000400_r8,  0.0000404_r8,  0.0000407_r8,  0.0000410_r8,  &
       0.0000413_r8,  0.0000416_r8,  0.0000419_r8,  0.0000422_r8,  0.0000425_r8,  &
       0.0000428_r8,  0.0000431_r8,  0.0000434_r8,  0.0000437_r8,  0.0000439_r8,  &
       0.0000442_r8,  0.0000445_r8,  0.0000447_r8,  0.0000450_r8,  0.0000453_r8,  &
       0.0000455_r8,  0.0000458_r8,  0.0000460_r8,  0.0000463_r8,  0.0000465_r8,  &
       0.0000468_r8,  0.0000470_r8, &
       0.0000173_r8,  0.0000190_r8,  0.0000205_r8,  0.0000219_r8,  0.0000231_r8,  &
       0.0000242_r8,  0.0000252_r8,  0.0000262_r8,  0.0000271_r8,  0.0000279_r8,  &
       0.0000287_r8,  0.0000294_r8,  0.0000301_r8,  0.0000308_r8,  0.0000314_r8,  &
       0.0000320_r8,  0.0000326_r8,  0.0000332_r8,  0.0000330_r8,  0.0000343_r8,  &
       0.0000348_r8,  0.0000353_r8,  0.0000358_r8,  0.0000362_r8,  0.0000367_r8,  &
       0.0000371_r8,  0.0000376_r8,  0.0000380_r8,  0.0000384_r8,  0.0000388_r8,  &
       0.0000392_r8,  0.0000396_r8,  0.0000399_r8,  0.0000403_r8,  0.0000407_r8,  &
       0.0000410_r8,  0.0000414_r8,  0.0000417_r8,  0.0000420_r8,  0.0000424_r8,  &
       0.0000427_r8,  0.0000430_r8,  0.0000433_r8,  0.0000436_r8,  0.0000439_r8,  &
       0.0000442_r8,  0.0000445_r8,  0.0000448_r8,  0.0000451_r8,  0.0000454_r8,  &
       0.0000457_r8,  0.0000459_r8,  0.0000462_r8,  0.0000465_r8,  0.0000468_r8,  &
       0.0000470_r8,  0.0000473_r8,  0.0000475_r8,  0.0000478_r8,  0.0000481_r8,  &
       0.0000483_r8,  0.0000486_r8, &
       0.0000186_r8,  0.0000204_r8,  0.0000219_r8,  0.0000233_r8,  0.0000246_r8,  &
       0.0000257_r8,  0.0000268_r8,  0.0000277_r8,  0.0000286_r8,  0.0000295_r8,  &
       0.0000303_r8,  0.0000311_r8,  0.0000318_r8,  0.0000325_r8,  0.0000331_r8,  &
       0.0000331_r8,  0.0000344_r8,  0.0000350_r8,  0.0000355_r8,  0.0000361_r8,  &
       0.0000366_r8,  0.0000371_r8,  0.0000376_r8,  0.0000381_r8,  0.0000386_r8,  &
       0.0000390_r8,  0.0000395_r8,  0.0000399_r8,  0.0000403_r8,  0.0000407_r8,  &
       0.0000412_r8,  0.0000416_r8,  0.0000419_r8,  0.0000423_r8,  0.0000427_r8,  &
       0.0000431_r8,  0.0000434_r8,  0.0000438_r8,  0.0000441_r8,  0.0000445_r8,  &
       0.0000448_r8,  0.0000451_r8,  0.0000455_r8,  0.0000458_r8,  0.0000461_r8,  &
       0.0000464_r8,  0.0000467_r8,  0.0000470_r8,  0.0000473_r8,  0.0000476_r8,  &
       0.0000479_r8,  0.0000482_r8,  0.0000485_r8,  0.0000488_r8,  0.0000491_r8,  &
       0.0000494_r8,  0.0000497_r8,  0.0000499_r8,  0.0000502_r8,  0.0000505_r8,  &
       0.0000507_r8,  0.0000510_r8, &
       0.0000198_r8,  0.0000216_r8,  0.0000232_r8,  0.0000246_r8,  0.0000259_r8,  &
       0.0000271_r8,  0.0000281_r8,  0.0000291_r8,  0.0000301_r8,  0.0000310_r8,  &
       0.0000318_r8,  0.0000326_r8,  0.0000333_r8,  0.0000340_r8,  0.0000347_r8,  &
       0.0000354_r8,  0.0000360_r8,  0.0000366_r8,  0.0000372_r8,  0.0000377_r8,  &
       0.0000383_r8,  0.0000388_r8,  0.0000393_r8,  0.0000398_r8,  0.0000403_r8,  &
       0.0000408_r8,  0.0000412_r8,  0.0000417_r8,  0.0000421_r8,  0.0000425_r8,  &
       0.0000430_r8,  0.0000434_r8,  0.0000438_r8,  0.0000442_r8,  0.0000446_r8,  &
       0.0000449_r8,  0.0000453_r8,  0.0000457_r8,  0.0000461_r8,  0.0000464_r8,  &
       0.0000468_r8,  0.0000471_r8,  0.0000475_r8,  0.0000478_r8,  0.0000481_r8,  &
       0.0000485_r8,  0.0000488_r8,  0.0000491_r8,  0.0000494_r8,  0.0000498_r8,  &
       0.0000501_r8,  0.0000504_r8,  0.0000507_r8,  0.0000510_r8,  0.0000513_r8,  &
       0.0000516_r8,  0.0000519_r8,  0.0000522_r8,  0.0000524_r8,  0.0000527_r8,  &
       0.0000530_r8,  0.0000533_r8, &
       0.0000209_r8,  0.0000228_r8,  0.0000244_r8,  0.0000258_r8,  0.0000271_r8,  &
       0.0000283_r8,  0.0000294_r8,  0.0000305_r8,  0.0000314_r8,  0.0000323_r8,  &
       0.0000332_r8,  0.0000340_r8,  0.0000347_r8,  0.0000354_r8,  0.0000361_r8,  &
       0.0000368_r8,  0.0000375_r8,  0.0000381_r8,  0.0000387_r8,  0.0000392_r8,  &
       0.0000398_r8,  0.0000404_r8,  0.0000409_r8,  0.0000414_r8,  0.0000419_r8,  &
       0.0000424_r8,  0.0000429_r8,  0.0000433_r8,  0.0000438_r8,  0.0000442_r8,  &
       0.0000447_r8,  0.0000451_r8,  0.0000455_r8,  0.0000459_r8,  0.0000463_r8,  &
       0.0000467_r8,  0.0000471_r8,  0.0000475_r8,  0.0000479_r8,  0.0000483_r8,  &
       0.0000486_r8,  0.0000490_r8,  0.0000493_r8,  0.0000497_r8,  0.0000501_r8,  &
       0.0000504_r8,  0.0000507_r8,  0.0000511_r8,  0.0000514_r8,  0.0000518_r8,  &
       0.0000521_r8,  0.0000524_r8,  0.0000527_r8,  0.0000530_r8,  0.0000534_r8,  &
       0.0000537_r8,  0.0000540_r8,  0.0000543_r8,  0.0000546_r8,  0.0000549_r8,  &
       0.0000552_r8,  0.0000555_r8, &
       0.0000221_r8,  0.0000240_r8,  0.0000257_r8,  0.0000272_r8,  0.0000285_r8,  &
       0.0000297_r8,  0.0000308_r8,  0.0000319_r8,  0.0000329_r8,  0.0000331_r8,  &
       0.0000347_r8,  0.0000355_r8,  0.0000363_r8,  0.0000370_r8,  0.0000377_r8,  &
       0.0000384_r8,  0.0000391_r8,  0.0000397_r8,  0.0000404_r8,  0.0000409_r8,  &
       0.0000415_r8,  0.0000421_r8,  0.0000426_r8,  0.0000432_r8,  0.0000437_r8,  &
       0.0000442_r8,  0.0000447_r8,  0.0000452_r8,  0.0000456_r8,  0.0000461_r8,  &
       0.0000466_r8,  0.0000470_r8,  0.0000475_r8,  0.0000479_r8,  0.0000483_r8,  &
       0.0000487_r8,  0.0000491_r8,  0.0000496_r8,  0.0000500_r8,  0.0000503_r8,  &
       0.0000507_r8,  0.0000511_r8,  0.0000515_r8,  0.0000519_r8,  0.0000523_r8,  &
       0.0000526_r8,  0.0000530_r8,  0.0000533_r8,  0.0000537_r8,  0.0000540_r8,  &
       0.0000544_r8,  0.0000547_r8,  0.0000551_r8,  0.0000554_r8,  0.0000558_r8,  &
       0.0000561_r8,  0.0000564_r8,  0.0000567_r8,  0.0000571_r8,  0.0000574_r8,  &
       0.0000577_r8,  0.0000580_r8, &
       0.0000234_r8,  0.0000254_r8,  0.0000271_r8,  0.0000286_r8,  0.0000300_r8,  &
       0.0000312_r8,  0.0000324_r8,  0.0000335_r8,  0.0000345_r8,  0.0000354_r8,  &
       0.0000363_r8,  0.0000372_r8,  0.0000380_r8,  0.0000387_r8,  0.0000395_r8,  &
       0.0000402_r8,  0.0000409_r8,  0.0000415_r8,  0.0000422_r8,  0.0000428_r8,  &
       0.0000434_r8,  0.0000440_r8,  0.0000446_r8,  0.0000451_r8,  0.0000457_r8,  &
       0.0000462_r8,  0.0000467_r8,  0.0000472_r8,  0.0000477_r8,  0.0000482_r8,  &
       0.0000487_r8,  0.0000492_r8,  0.0000496_r8,  0.0000501_r8,  0.0000505_r8,  &
       0.0000510_r8,  0.0000514_r8,  0.0000518_r8,  0.0000523_r8,  0.0000527_r8,  &
       0.0000531_r8,  0.0000535_r8,  0.0000539_r8,  0.0000543_r8,  0.0000547_r8,  &
       0.0000551_r8,  0.0000555_r8,  0.0000559_r8,  0.0000562_r8,  0.0000566_r8,  &
       0.0000570_r8,  0.0000573_r8,  0.0000577_r8,  0.0000581_r8,  0.0000584_r8,  &
       0.0000588_r8,  0.0000591_r8,  0.0000595_r8,  0.0000598_r8,  0.0000602_r8,  &
       0.0000605_r8,  0.0000608_r8, &
       0.0000248_r8,  0.0000268_r8,  0.0000285_r8,  0.0000301_r8,  0.0000315_r8,  &
       0.0000328_r8,  0.0000340_r8,  0.0000351_r8,  0.0000362_r8,  0.0000371_r8,  &
       0.0000381_r8,  0.0000389_r8,  0.0000398_r8,  0.0000406_r8,  0.0000413_r8,  &
       0.0000421_r8,  0.0000428_r8,  0.0000435_r8,  0.0000442_r8,  0.0000448_r8,  &
       0.0000454_r8,  0.0000460_r8,  0.0000466_r8,  0.0000472_r8,  0.0000478_r8,  &
       0.0000484_r8,  0.0000489_r8,  0.0000494_r8,  0.0000500_r8,  0.0000505_r8,  &
       0.0000510_r8,  0.0000515_r8,  0.0000520_r8,  0.0000525_r8,  0.0000530_r8,  &
       0.0000534_r8,  0.0000539_r8,  0.0000544_r8,  0.0000548_r8,  0.0000553_r8,  &
       0.0000557_r8,  0.0000561_r8,  0.0000566_r8,  0.0000570_r8,  0.0000574_r8,  &
       0.0000578_r8,  0.0000582_r8,  0.0000586_r8,  0.0000590_r8,  0.0000594_r8,  &
       0.0000598_r8,  0.0000602_r8,  0.0000606_r8,  0.0000610_r8,  0.0000614_r8,  &
       0.0000618_r8,  0.0000621_r8,  0.0000625_r8,  0.0000629_r8,  0.0000633_r8,  &
       0.0000636_r8,  0.0000640_r8, &
       0.0000260_r8,  0.0000281_r8,  0.0000299_r8,  0.0000315_r8,  0.0000330_r8,  &
       0.0000343_r8,  0.0000355_r8,  0.0000367_r8,  0.0000377_r8,  0.0000388_r8,  &
       0.0000397_r8,  0.0000406_r8,  0.0000415_r8,  0.0000423_r8,  0.0000431_r8,  &
       0.0000439_r8,  0.0000446_r8,  0.0000453_r8,  0.0000460_r8,  0.0000467_r8,  &
       0.0000474_r8,  0.0000480_r8,  0.0000487_r8,  0.0000493_r8,  0.0000499_r8,  &
       0.0000505_r8,  0.0000510_r8,  0.0000516_r8,  0.0000522_r8,  0.0000527_r8,  &
       0.0000533_r8,  0.0000538_r8,  0.0000543_r8,  0.0000548_r8,  0.0000553_r8,  &
       0.0000558_r8,  0.0000563_r8,  0.0000568_r8,  0.0000573_r8,  0.0000578_r8,  &
       0.0000582_r8,  0.0000587_r8,  0.0000591_r8,  0.0000596_r8,  0.0000601_r8,  &
       0.0000605_r8,  0.0000609_r8,  0.0000614_r8,  0.0000618_r8,  0.0000622_r8,  &
       0.0000626_r8,  0.0000631_r8,  0.0000635_r8,  0.0000639_r8,  0.0000643_r8,  &
       0.0000647_r8,  0.0000651_r8,  0.0000655_r8,  0.0000659_r8,  0.0000663_r8,  &
       0.0000667_r8,  0.0000670_r8, &
       0.0000275_r8,  0.0000296_r8,  0.0000315_r8,  0.0000332_r8,  0.0000347_r8,  &
       0.0000360_r8,  0.0000373_r8,  0.0000385_r8,  0.0000396_r8,  0.0000407_r8,  &
       0.0000417_r8,  0.0000426_r8,  0.0000435_r8,  0.0000444_r8,  0.0000452_r8,  &
       0.0000460_r8,  0.0000468_r8,  0.0000476_r8,  0.0000483_r8,  0.0000490_r8,  &
       0.0000497_r8,  0.0000504_r8,  0.0000511_r8,  0.0000517_r8,  0.0000524_r8,  &
       0.0000530_r8,  0.0000536_r8,  0.0000542_r8,  0.0000548_r8,  0.0000554_r8,  &
       0.0000560_r8,  0.0000566_r8,  0.0000571_r8,  0.0000577_r8,  0.0000582_r8,  &
       0.0000587_r8,  0.0000593_r8,  0.0000598_r8,  0.0000603_r8,  0.0000608_r8,  &
       0.0000613_r8,  0.0000618_r8,  0.0000623_r8,  0.0000628_r8,  0.0000633_r8,  &
       0.0000638_r8,  0.0000642_r8,  0.0000647_r8,  0.0000652_r8,  0.0000656_r8,  &
       0.0000661_r8,  0.0000665_r8,  0.0000670_r8,  0.0000674_r8,  0.0000678_r8,  &
       0.0000683_r8,  0.0000687_r8,  0.0000691_r8,  0.0000695_r8,  0.0000700_r8,  &
       0.0000704_r8,  0.0000708_r8, &
       0.0000290_r8,  0.0000312_r8,  0.0000331_r8,  0.0000349_r8,  0.0000364_r8,  &
       0.0000379_r8,  0.0000392_r8,  0.0000404_r8,  0.0000416_r8,  0.0000427_r8,  &
       0.0000437_r8,  0.0000447_r8,  0.0000457_r8,  0.0000466_r8,  0.0000475_r8,  &
       0.0000483_r8,  0.0000492_r8,  0.0000500_r8,  0.0000507_r8,  0.0000515_r8,  &
       0.0000523_r8,  0.0000530_r8,  0.0000537_r8,  0.0000544_r8,  0.0000551_r8,  &
       0.0000558_r8,  0.0000564_r8,  0.0000571_r8,  0.0000577_r8,  0.0000583_r8,  &
       0.0000589_r8,  0.0000596_r8,  0.0000602_r8,  0.0000607_r8,  0.0000613_r8,  &
       0.0000619_r8,  0.0000625_r8,  0.0000630_r8,  0.0000636_r8,  0.0000641_r8,  &
       0.0000647_r8,  0.0000652_r8,  0.0000657_r8,  0.0000663_r8,  0.0000668_r8,  &
       0.0000673_r8,  0.0000678_r8,  0.0000683_r8,  0.0000688_r8,  0.0000693_r8,  &
       0.0000698_r8,  0.0000702_r8,  0.0000707_r8,  0.0000712_r8,  0.0000716_r8,  &
       0.0000721_r8,  0.0000726_r8,  0.0000730_r8,  0.0000735_r8,  0.0000739_r8,  &
       0.0000744_r8,  0.0000748_r8, &
       0.0000306_r8,  0.0000329_r8,  0.0000349_r8,  0.0000366_r8,  0.0000383_r8,  &
       0.0000398_r8,  0.0000411_r8,  0.0000424_r8,  0.0000436_r8,  0.0000448_r8,  &
       0.0000459_r8,  0.0000469_r8,  0.0000479_r8,  0.0000489_r8,  0.0000499_r8,  &
       0.0000508_r8,  0.0000516_r8,  0.0000525_r8,  0.0000533_r8,  0.0000542_r8,  &
       0.0000549_r8,  0.0000557_r8,  0.0000565_r8,  0.0000572_r8,  0.0000580_r8,  &
       0.0000587_r8,  0.0000594_r8,  0.0000601_r8,  0.0000608_r8,  0.0000615_r8,  &
       0.0000621_r8,  0.0000628_r8,  0.0000634_r8,  0.0000640_r8,  0.0000647_r8,  &
       0.0000653_r8,  0.0000659_r8,  0.0000665_r8,  0.0000671_r8,  0.0000677_r8,  &
       0.0000683_r8,  0.0000688_r8,  0.0000694_r8,  0.0000700_r8,  0.0000705_r8,  &
       0.0000711_r8,  0.0000716_r8,  0.0000721_r8,  0.0000727_r8,  0.0000732_r8,  &
       0.0000737_r8,  0.0000742_r8,  0.0000747_r8,  0.0000752_r8,  0.0000757_r8,  &
       0.0000762_r8,  0.0000767_r8,  0.0000772_r8,  0.0000777_r8,  0.0000782_r8,  &
       0.0000786_r8,  0.0000791_r8, &
       0.0000323_r8,  0.0000347_r8,  0.0000368_r8,  0.0000386_r8,  0.0000403_r8,  &
       0.0000419_r8,  0.0000433_r8,  0.0000447_r8,  0.0000459_r8,  0.0000472_r8,  &
       0.0000483_r8,  0.0000494_r8,  0.0000505_r8,  0.0000516_r8,  0.0000526_r8,  &
       0.0000535_r8,  0.0000545_r8,  0.0000554_r8,  0.0000563_r8,  0.0000572_r8,  &
       0.0000580_r8,  0.0000589_r8,  0.0000597_r8,  0.0000605_r8,  0.0000613_r8,  &
       0.0000621_r8,  0.0000628_r8,  0.0000636_r8,  0.0000643_r8,  0.0000650_r8,  &
       0.0000657_r8,  0.0000664_r8,  0.0000671_r8,  0.0000678_r8,  0.0000685_r8,  &
       0.0000692_r8,  0.0000698_r8,  0.0000705_r8,  0.0000711_r8,  0.0000717_r8,  &
       0.0000724_r8,  0.0000730_r8,  0.0000736_r8,  0.0000742_r8,  0.0000748_r8,  &
       0.0000754_r8,  0.0000760_r8,  0.0000765_r8,  0.0000771_r8,  0.0000777_r8,  &
       0.0000782_r8,  0.0000788_r8,  0.0000793_r8,  0.0000799_r8,  0.0000804_r8,  &
       0.0000809_r8,  0.0000815_r8,  0.0000820_r8,  0.0000825_r8,  0.0000830_r8,  &
       0.0000835_r8,  0.0000840_r8, &
       0.0000341_r8,  0.0000365_r8,  0.0000387_r8,  0.0000406_r8,  0.0000424_r8,  &
       0.0000440_r8,  0.0000456_r8,  0.0000470_r8,  0.0000483_r8,  0.0000496_r8,  &
       0.0000509_r8,  0.0000521_r8,  0.0000532_r8,  0.0000543_r8,  0.0000554_r8,  &
       0.0000564_r8,  0.0000574_r8,  0.0000584_r8,  0.0000594_r8,  0.0000603_r8,  &
       0.0000613_r8,  0.0000622_r8,  0.0000630_r8,  0.0000639_r8,  0.0000648_r8,  &
       0.0000656_r8,  0.0000664_r8,  0.0000672_r8,  0.0000680_r8,  0.0000688_r8,  &
       0.0000696_r8,  0.0000703_r8,  0.0000711_r8,  0.0000718_r8,  0.0000725_r8,  &
       0.0000732_r8,  0.0000739_r8,  0.0000746_r8,  0.0000753_r8,  0.0000760_r8,  &
       0.0000767_r8,  0.0000773_r8,  0.0000780_r8,  0.0000786_r8,  0.0000793_r8,  &
       0.0000799_r8,  0.0000805_r8,  0.0000811_r8,  0.0000817_r8,  0.0000823_r8,  &
       0.0000829_r8,  0.0000835_r8,  0.0000841_r8,  0.0000847_r8,  0.0000853_r8,  &
       0.0000858_r8,  0.0000864_r8,  0.0000870_r8,  0.0000875_r8,  0.0000881_r8,  &
       0.0000886_r8,  0.0000892_r8, &
       0.0000359_r8,  0.0000385_r8,  0.0000408_r8,  0.0000428_r8,  0.0000447_r8, &
       0.0000464_r8,  0.0000480_r8,  0.0000495_r8,  0.0000510_r8,  0.0000524_r8, &
       0.0000537_r8,  0.0000550_r8,  0.0000562_r8,  0.0000574_r8,  0.0000585_r8, &
       0.0000597_r8,  0.0000608_r8,  0.0000618_r8,  0.0000629_r8,  0.0000639_r8, &
       0.0000649_r8,  0.0000658_r8,  0.0000668_r8,  0.0000677_r8,  0.0000686_r8, &
       0.0000695_r8,  0.0000704_r8,  0.0000713_r8,  0.0000721_r8,  0.0000730_r8, &
       0.0000738_r8,  0.0000746_r8,  0.0000754_r8,  0.0000762_r8,  0.0000770_r8, &
       0.0000777_r8,  0.0000785_r8,  0.0000792_r8,  0.0000800_r8,  0.0000807_r8, &
       0.0000814_r8,  0.0000821_r8,  0.0000828_r8,  0.0000835_r8,  0.0000842_r8, &
       0.0000849_r8,  0.0000856_r8,  0.0000862_r8,  0.0000869_r8,  0.0000875_r8, &
       0.0000882_r8,  0.0000888_r8,  0.0000894_r8,  0.0000900_r8,  0.0000907_r8, &
       0.0000913_r8,  0.0000919_r8,  0.0000925_r8,  0.0000931_r8,  0.0000936_r8, &
       0.0000942_r8,  0.0000948_r8, &
       0.0000380_r8,  0.0000407_r8,  0.0000431_r8,  0.0000453_r8,  0.0000473_r8,  &
       0.0000491_r8,  0.0000508_r8,  0.0000525_r8,  0.0000540_r8,  0.0000555_r8,  &
       0.0000569_r8,  0.0000583_r8,  0.0000596_r8,  0.0000609_r8,  0.0000622_r8,  &
       0.0000634_r8,  0.0000646_r8,  0.0000657_r8,  0.0000668_r8,  0.0000679_r8,  &
       0.0000690_r8,  0.0000700_r8,  0.0000711_r8,  0.0000721_r8,  0.0000731_r8,  &
       0.0000740_r8,  0.0000750_r8,  0.0000759_r8,  0.0000769_r8,  0.0000778_r8,  &
       0.0000786_r8,  0.0000795_r8,  0.0000804_r8,  0.0000812_r8,  0.0000821_r8,  &
       0.0000829_r8,  0.0000837_r8,  0.0000845_r8,  0.0000853_r8,  0.0000861_r8,  &
       0.0000869_r8,  0.0000876_r8,  0.0000884_r8,  0.0000891_r8,  0.0000899_r8,  &
       0.0000906_r8,  0.0000913_r8,  0.0000920_r8,  0.0000927_r8,  0.0000934_r8,  &
       0.0000941_r8,  0.0000948_r8,  0.0000955_r8,  0.0000961_r8,  0.0000968_r8,  &
       0.0000974_r8,  0.0000981_r8,  0.0000987_r8,  0.0000994_r8,  0.0001000_r8,  &
       0.0001006_r8,  0.0001012_r8, &
       0.0000403_r8,  0.0000431_r8,  0.0000456_r8,  0.0000479_r8,  0.0000500_r8,  &
       0.0000520_r8,  0.0000538_r8,  0.0000556_r8,  0.0000573_r8,  0.0000589_r8,  &  
       0.0000604_r8,  0.0000619_r8,  0.0000633_r8,  0.0000647_r8,  0.0000661_r8,  &
       0.0000674_r8,  0.0000686_r8,  0.0000699_r8,  0.0000711_r8,  0.0000723_r8,  &
       0.0000734_r8,  0.0000746_r8,  0.0000757_r8,  0.0000768_r8,  0.0000778_r8,  &
       0.0000789_r8,  0.0000799_r8,  0.0000809_r8,  0.0000819_r8,  0.0000829_r8,  &
       0.0000838_r8,  0.0000848_r8,  0.0000857_r8,  0.0000866_r8,  0.0000875_r8,  &
       0.0000884_r8,  0.0000893_r8,  0.0000902_r8,  0.0000910_r8,  0.0000919_r8,  &
       0.0000927_r8,  0.0000935_r8,  0.0000943_r8,  0.0000951_r8,  0.0000959_r8,  &
       0.0000967_r8,  0.0000974_r8,  0.0000982_r8,  0.0000990_r8,  0.0000997_r8,  &
       0.0001004_r8,  0.0001012_r8,  0.0001019_r8,  0.0001026_r8,  0.0001033_r8,  &
       0.0001040_r8,  0.0001047_r8,  0.0001054_r8,  0.0001061_r8,  0.0001067_r8,  &
       0.0001074_r8,  0.0001080_r8, &
       0.0000426_r8,  0.0000456_r8,  0.0000482_r8,  0.0000507_r8,  0.0000529_r8,  &
       0.0000550_r8,  0.0000570_r8,  0.0000589_r8,  0.0000607_r8,  0.0000624_r8,  &
       0.0000641_r8,  0.0000657_r8,  0.0000672_r8,  0.0000687_r8,  0.0000702_r8,  &
       0.0000716_r8,  0.0000730_r8,  0.0000743_r8,  0.0000756_r8,  0.0000769_r8,  &
       0.0000781_r8,  0.0000794_r8,  0.0000806_r8,  0.0000817_r8,  0.0000829_r8,  &
       0.0000840_r8,  0.0000851_r8,  0.0000862_r8,  0.0000873_r8,  0.0000883_r8,  &
       0.0000893_r8,  0.0000904_r8,  0.0000913_r8,  0.0000923_r8,  0.0000933_r8,  &
       0.0000943_r8,  0.0000952_r8,  0.0000961_r8,  0.0000970_r8,  0.0000979_r8,  &
       0.0000988_r8,  0.0000997_r8,  0.0001006_r8,  0.0001014_r8,  0.0001023_r8,  &
       0.0001031_r8,  0.0001039_r8,  0.0001047_r8,  0.0001055_r8,  0.0001063_r8,  &
       0.0001071_r8,  0.0001079_r8,  0.0001087_r8,  0.0001094_r8,  0.0001102_r8,  &
       0.0001109_r8,  0.0001116_r8,  0.0001124_r8,  0.0001131_r8,  0.0001138_r8,  &
       0.0001145_r8,  0.0001152_r8, &
       0.0000451_r8,  0.0000482_r8,  0.0000511_r8,  0.0000537_r8,  0.0000561_r8,  &
       0.0000584_r8,  0.0000605_r8,  0.0000626_r8,  0.0000645_r8,  0.0000664_r8,  &
       0.0000682_r8,  0.0000699_r8,  0.0000715_r8,  0.0000732_r8,  0.0000747_r8,  &
       0.0000763_r8,  0.0000777_r8,  0.0000792_r8,  0.0000806_r8,  0.0000820_r8,  &
       0.0000833_r8,  0.0000846_r8,  0.0000859_r8,  0.0000872_r8,  0.0000884_r8,  &
       0.0000896_r8,  0.0000908_r8,  0.0000920_r8,  0.0000931_r8,  0.0000942_r8,  &
       0.0000953_r8,  0.0000964_r8,  0.0000975_r8,  0.0000986_r8,  0.0000996_r8,  &
       0.0001006_r8,  0.0001016_r8,  0.0001026_r8,  0.0001036_r8,  0.0001046_r8,  &
       0.0001055_r8,  0.0001064_r8,  0.0001074_r8,  0.0001083_r8,  0.0001092_r8,  &
       0.0001101_r8,  0.0001110_r8,  0.0001118_r8,  0.0001127_r8,  0.0001135_r8,  &
       0.0001144_r8,  0.0001152_r8,  0.0001160_r8,  0.0001168_r8,  0.0001176_r8,  &
       0.0001184_r8,  0.0001192_r8,  0.0001200_r8,  0.0001207_r8,  0.0001215_r8,  &
       0.0001222_r8,  0.0001230_r8, &
       0.0000478_r8,  0.0000512_r8,  0.0000543_r8,  0.0000571_r8,  0.0000597_r8,  &
       0.0000621_r8,  0.0000644_r8,  0.0000666_r8,  0.0000687_r8,  0.0000708_r8,  &
       0.0000727_r8,  0.0000746_r8,  0.0000764_r8,  0.0000781_r8,  0.0000798_r8,  &
       0.0000814_r8,  0.0000830_r8,  0.0000846_r8,  0.0000861_r8,  0.0000876_r8,  &
       0.0000891_r8,  0.0000905_r8,  0.0000919_r8,  0.0000932_r8,  0.0000945_r8,  &
       0.0000958_r8,  0.0000971_r8,  0.0000984_r8,  0.0000996_r8,  0.0001008_r8,  &
       0.0001020_r8,  0.0001032_r8,  0.0001043_r8,  0.0001055_r8,  0.0001066_r8,  &
       0.0001077_r8,  0.0001088_r8,  0.0001098_r8,  0.0001109_r8,  0.0001119_r8,  &
       0.0001129_r8,  0.0001139_r8,  0.0001149_r8,  0.0001159_r8,  0.0001168_r8,  &
       0.0001178_r8,  0.0001187_r8,  0.0001197_r8,  0.0001206_r8,  0.0001215_r8,  &
       0.0001224_r8,  0.0001233_r8,  0.0001241_r8,  0.0001250_r8,  0.0001258_r8,  &
       0.0001267_r8,  0.0001275_r8,  0.0001283_r8,  0.0001292_r8,  0.0001300_r8,  &
       0.0001308_r8,  0.0001316_r8, &
       0.0000508_r8,  0.0000544_r8,  0.0000577_r8,  0.0000607_r8,  0.0000635_r8,  &
       0.0000661_r8,  0.0000686_r8,  0.0000710_r8,  0.0000733_r8,  0.0000754_r8,  &
       0.0000775_r8,  0.0000795_r8,  0.0000815_r8,  0.0000834_r8,  0.0000852_r8,  &
       0.0000870_r8,  0.0000887_r8,  0.0000904_r8,  0.0000920_r8,  0.0000936_r8,  &
       0.0000952_r8,  0.0000967_r8,  0.0000982_r8,  0.0000996_r8,  0.0001011_r8,  &
       0.0001025_r8,  0.0001038_r8,  0.0001052_r8,  0.0001065_r8,  0.0001078_r8,  &
       0.0001091_r8,  0.0001103_r8,  0.0001116_r8,  0.0001128_r8,  0.0001140_r8,  &
       0.0001151_r8,  0.0001163_r8,  0.0001174_r8,  0.0001186_r8,  0.0001197_r8,  &
       0.0001207_r8,  0.0001218_r8,  0.0001229_r8,  0.0001239_r8,  0.0001249_r8,  &
       0.0001260_r8,  0.0001270_r8,  0.0001279_r8,  0.0001289_r8,  0.0001299_r8,  &
       0.0001308_r8,  0.0001318_r8,  0.0001327_r8,  0.0001336_r8,  0.0001317_r8,  &
       0.0001325_r8,  0.0001363_r8,  0.0001372_r8,  0.0001380_r8,  0.0001389_r8,  &
       0.0001397_r8,  0.0001406_r8, &
       0.0000540_r8,  0.0000579_r8,  0.0000615_r8,  0.0000647_r8,  0.0000677_r8,  &
       0.0000706_r8,  0.0000733_r8,  0.0000758_r8,  0.0000783_r8,  0.0000806_r8,  &
       0.0000829_r8,  0.0000851_r8,  0.0000872_r8,  0.0000892_r8,  0.0000912_r8,  &
       0.0000931_r8,  0.0000950_r8,  0.0000968_r8,  0.0000985_r8,  0.0001003_r8,  &
       0.0001020_r8,  0.0001036_r8,  0.0001052_r8,  0.0001068_r8,  0.0001083_r8,  &
       0.0001098_r8,  0.0001113_r8,  0.0001127_r8,  0.0001142_r8,  0.0001156_r8,  &
       0.0001169_r8,  0.0001183_r8,  0.0001196_r8,  0.0001209_r8,  0.0001222_r8,  &
       0.0001234_r8,  0.0001246_r8,  0.0001259_r8,  0.0001270_r8,  0.0001282_r8,  &
       0.0001294_r8,  0.0001305_r8,  0.0001317_r8,  0.0001328_r8,  0.0001339_r8,  &
       0.0001321_r8,  0.0001360_r8,  0.0001371_r8,  0.0001381_r8,  0.0001391_r8,  &
       0.0001401_r8,  0.0001411_r8,  0.0001421_r8,  0.0001431_r8,  0.0001440_r8,  &
       0.0001450_r8,  0.0001459_r8,  0.0001469_r8,  0.0001478_r8,  0.0001487_r8,  &
       0.0001496_r8,  0.0001505_r8, &
       0.0000575_r8,  0.0000617_r8,  0.0000655_r8,  0.0000690_r8,  0.0000723_r8,  &
       0.0000754_r8,  0.0000783_r8,  0.0000810_r8,  0.0000837_r8,  0.0000862_r8,  &
       0.0000887_r8,  0.0000910_r8,  0.0000933_r8,  0.0000955_r8,  0.0000976_r8,  &
       0.0000997_r8,  0.0001017_r8,  0.0001036_r8,  0.0001055_r8,  0.0001074_r8,  &
       0.0001092_r8,  0.0001110_r8,  0.0001127_r8,  0.0001144_r8,  0.0001160_r8,  &
       0.0001176_r8,  0.0001192_r8,  0.0001208_r8,  0.0001223_r8,  0.0001238_r8,  &
       0.0001252_r8,  0.0001267_r8,  0.0001281_r8,  0.0001295_r8,  0.0001308_r8,  &
       0.0001322_r8,  0.0001335_r8,  0.0001319_r8,  0.0001360_r8,  0.0001373_r8,  &
       0.0001385_r8,  0.0001397_r8,  0.0001409_r8,  0.0001421_r8,  0.0001433_r8,  &
       0.0001444_r8,  0.0001456_r8,  0.0001467_r8,  0.0001478_r8,  0.0001489_r8,  &
       0.0001499_r8,  0.0001510_r8,  0.0001520_r8,  0.0001531_r8,  0.0001541_r8,  &
       0.0001551_r8,  0.0001561_r8,  0.0001571_r8,  0.0001581_r8,  0.0001590_r8,  &
       0.0001600_r8,  0.0001609_r8, &
       0.0000613_r8,  0.0000659_r8,  0.0000700_r8,  0.0000738_r8,  0.0000773_r8,  &
       0.0000806_r8,  0.0000838_r8,  0.0000868_r8,  0.0000896_r8,  0.0000924_r8,  &
       0.0000950_r8,  0.0000976_r8,  0.0001000_r8,  0.0001024_r8,  0.0001047_r8,  &
       0.0001069_r8,  0.0001091_r8,  0.0001112_r8,  0.0001132_r8,  0.0001152_r8,  &
       0.0001172_r8,  0.0001191_r8,  0.0001209_r8,  0.0001227_r8,  0.0001245_r8,  &
       0.0001262_r8,  0.0001279_r8,  0.0001296_r8,  0.0001312_r8,  0.0001328_r8,  &
       0.0001344_r8,  0.0001359_r8,  0.0001374_r8,  0.0001389_r8,  0.0001403_r8,  &
       0.0001417_r8,  0.0001432_r8,  0.0001445_r8,  0.0001459_r8,  0.0001472_r8,  &
       0.0001485_r8,  0.0001498_r8,  0.0001511_r8,  0.0001524_r8,  0.0001536_r8,  &
       0.0001548_r8,  0.0001560_r8,  0.0001572_r8,  0.0001584_r8,  0.0001595_r8,  &
       0.0001607_r8,  0.0001618_r8,  0.0001629_r8,  0.0001640_r8,  0.0001651_r8,  &
       0.0001661_r8,  0.0001672_r8,  0.0001682_r8,  0.0001693_r8,  0.0001703_r8,  &
       0.0001713_r8,  0.0001723_r8, &
       0.0000654_r8,  0.0000703_r8,  0.0000747_r8,  0.0000789_r8,  0.0000827_r8,  &
       0.0000863_r8,  0.0000897_r8,  0.0000929_r8,  0.0000960_r8,  0.0000990_r8,  &
       0.0001018_r8,  0.0001046_r8,  0.0001072_r8,  0.0001098_r8,  0.0001123_r8,  &
       0.0001147_r8,  0.0001170_r8,  0.0001193_r8,  0.0001214_r8,  0.0001236_r8,  &
       0.0001257_r8,  0.0001277_r8,  0.0001297_r8,  0.0001316_r8,  0.0001335_r8,  &
       0.0001325_r8,  0.0001372_r8,  0.0001389_r8,  0.0001407_r8,  0.0001424_r8,  &
       0.0001440_r8,  0.0001457_r8,  0.0001473_r8,  0.0001488_r8,  0.0001504_r8,  &
       0.0001519_r8,  0.0001534_r8,  0.0001548_r8,  0.0001563_r8,  0.0001577_r8,  &
       0.0001591_r8,  0.0001605_r8,  0.0001618_r8,  0.0001631_r8,  0.0001645_r8,  &
       0.0001658_r8,  0.0001670_r8,  0.0001683_r8,  0.0001695_r8,  0.0001707_r8,  &
       0.0001720_r8,  0.0001732_r8,  0.0001743_r8,  0.0001755_r8,  0.0001767_r8,  &
       0.0001778_r8,  0.0001789_r8,  0.0001800_r8,  0.0001811_r8,  0.0001822_r8,  &
       0.0001833_r8,  0.0001844_r8, &
       0.0000699_r8,  0.0000752_r8,  0.0000800_r8,  0.0000844_r8,  0.0000886_r8,  &
       0.0000925_r8,  0.0000962_r8,  0.0000997_r8,  0.0001030_r8,  0.0001062_r8,  &
       0.0001093_r8,  0.0001123_r8,  0.0001151_r8,  0.0001179_r8,  0.0001205_r8,  &
       0.0001231_r8,  0.0001256_r8,  0.0001280_r8,  0.0001304_r8,  0.0001327_r8,  &
       0.0001321_r8,  0.0001371_r8,  0.0001392_r8,  0.0001413_r8,  0.0001433_r8,  &
       0.0001453_r8,  0.0001472_r8,  0.0001491_r8,  0.0001509_r8,  0.0001527_r8,  &
       0.0001545_r8,  0.0001562_r8,  0.0001579_r8,  0.0001596_r8,  0.0001612_r8,  &
       0.0001629_r8,  0.0001644_r8,  0.0001660_r8,  0.0001675_r8,  0.0001690_r8,  &
       0.0001705_r8,  0.0001720_r8,  0.0001734_r8,  0.0001749_r8,  0.0001762_r8,  &
       0.0001776_r8,  0.0001790_r8,  0.0001803_r8,  0.0001817_r8,  0.0001830_r8,  &
       0.0001842_r8,  0.0001855_r8,  0.0001868_r8,  0.0001880_r8,  0.0001892_r8,  &
       0.0001905_r8,  0.0001917_r8,  0.0001928_r8,  0.0001940_r8,  0.0001952_r8,  &
       0.0001963_r8,  0.0001975_r8, &
       0.0000748_r8,  0.0000805_r8,  0.0000858_r8,  0.0000906_r8,  0.0000951_r8, &
       0.0000993_r8,  0.0001033_r8,  0.0001071_r8,  0.0001107_r8,  0.0001142_r8, &
       0.0001175_r8,  0.0001207_r8,  0.0001238_r8,  0.0001267_r8,  0.0001296_r8, &
       0.0001323_r8,  0.0001322_r8,  0.0001376_r8,  0.0001401_r8,  0.0001426_r8, &
       0.0001450_r8,  0.0001473_r8,  0.0001496_r8,  0.0001518_r8,  0.0001539_r8, &
       0.0001560_r8,  0.0001581_r8,  0.0001601_r8,  0.0001620_r8,  0.0001640_r8, &
       0.0001659_r8,  0.0001677_r8,  0.0001695_r8,  0.0001713_r8,  0.0001731_r8, &
       0.0001748_r8,  0.0001765_r8,  0.0001781_r8,  0.0001798_r8,  0.0001814_r8, &
       0.0001830_r8,  0.0001845_r8,  0.0001861_r8,  0.0001876_r8,  0.0001891_r8, &
       0.0001905_r8,  0.0001920_r8,  0.0001934_r8,  0.0001948_r8,  0.0001962_r8, &
       0.0001976_r8,  0.0001990_r8,  0.0002003_r8,  0.0002017_r8,  0.0002030_r8, &
       0.0002043_r8,  0.0002056_r8,  0.0002068_r8,  0.0002081_r8,  0.0002093_r8, &
       0.0002106_r8,  0.0002118_r8, &
       0.0000802_r8,  0.0000863_r8,  0.0000920_r8,  0.0000972_r8,  0.0001021_r8, &
       0.0001067_r8,  0.0001110_r8,  0.0001151_r8,  0.0001190_r8,  0.0001227_r8, &
       0.0001263_r8,  0.0001297_r8,  0.0001330_r8,  0.0001362_r8,  0.0001393_r8, &
       0.0001422_r8,  0.0001451_r8,  0.0001479_r8,  0.0001506_r8,  0.0001532_r8, &
       0.0001557_r8,  0.0001582_r8,  0.0001606_r8,  0.0001630_r8,  0.0001653_r8, &
       0.0001675_r8,  0.0001697_r8,  0.0001719_r8,  0.0001740_r8,  0.0001760_r8, &
       0.0001780_r8,  0.0001800_r8,  0.0001819_r8,  0.0001839_r8,  0.0001857_r8, &
       0.0001876_r8,  0.0001894_r8,  0.0001911_r8,  0.0001929_r8,  0.0001946_r8, &
       0.0001963_r8,  0.0001980_r8,  0.0001996_r8,  0.0002012_r8,  0.0002028_r8, &
       0.0002044_r8,  0.0002060_r8,  0.0002075_r8,  0.0002090_r8,  0.0002105_r8, &
       0.0002120_r8,  0.0002135_r8,  0.0002149_r8,  0.0002164_r8,  0.0002178_r8, &
       0.0002192_r8,  0.0002205_r8,  0.0002219_r8,  0.0002233_r8,  0.0002246_r8, &
       0.0002259_r8,  0.0002273_r8, &
       0.0000859_r8,  0.0000926_r8,  0.0000987_r8,  0.0001044_r8,  0.0001097_r8, &
       0.0001146_r8,  0.0001193_r8,  0.0001237_r8,  0.0001279_r8,  0.0001319_r8, &
       0.0001358_r8,  0.0001395_r8,  0.0001430_r8,  0.0001464_r8,  0.0001497_r8, &
       0.0001528_r8,  0.0001559_r8,  0.0001589_r8,  0.0001617_r8,  0.0001645_r8, &
       0.0001673_r8,  0.0001699_r8,  0.0001725_r8,  0.0001750_r8,  0.0001774_r8, &
       0.0001798_r8,  0.0001822_r8,  0.0001845_r8,  0.0001867_r8,  0.0001889_r8, &
       0.0001911_r8,  0.0001932_r8,  0.0001953_r8,  0.0001973_r8,  0.0001993_r8, &
       0.0002013_r8,  0.0002032_r8,  0.0002051_r8,  0.0002070_r8,  0.0002088_r8, &
       0.0002107_r8,  0.0002124_r8,  0.0002142_r8,  0.0002160_r8,  0.0002177_r8, &
       0.0002194_r8,  0.0002211_r8,  0.0002227_r8,  0.0002243_r8,  0.0002260_r8, &
       0.0002276_r8,  0.0002291_r8,  0.0002307_r8,  0.0002322_r8,  0.0002338_r8, &
       0.0002353_r8,  0.0002368_r8,  0.0002382_r8,  0.0002397_r8,  0.0002412_r8, &
       0.0002426_r8,  0.0002440_r8, &
       0.0000922_r8,  0.0000995_r8,  0.0001061_r8,  0.0001122_r8,  0.0001179_r8, &
       0.0001233_r8,  0.0001283_r8,  0.0001331_r8,  0.0001376_r8,  0.0001419_r8, &
       0.0001460_r8,  0.0001500_r8,  0.0001538_r8,  0.0001574_r8,  0.0001609_r8, &
       0.0001643_r8,  0.0001676_r8,  0.0001707_r8,  0.0001738_r8,  0.0001768_r8, &
       0.0001797_r8,  0.0001825_r8,  0.0001853_r8,  0.0001880_r8,  0.0001906_r8, &
       0.0001932_r8,  0.0001957_r8,  0.0001981_r8,  0.0002006_r8,  0.0002029_r8, &
       0.0002052_r8,  0.0002075_r8,  0.0002097_r8,  0.0002119_r8,  0.0002141_r8, &
       0.0002162_r8,  0.0002183_r8,  0.0002203_r8,  0.0002223_r8,  0.0002243_r8, &
       0.0002263_r8,  0.0002282_r8,  0.0002301_r8,  0.0002320_r8,  0.0002339_r8, &
       0.0002357_r8,  0.0002375_r8,  0.0002393_r8,  0.0002411_r8,  0.0002428_r8, &
       0.0002446_r8,  0.0002463_r8,  0.0002480_r8,  0.0002496_r8,  0.0002513_r8, &
       0.0002529_r8,  0.0002546_r8,  0.0002562_r8,  0.0002578_r8,  0.0002593_r8, &
       0.0002609_r8,  0.0002625_r8, &
       0.0000990_r8,  0.0001069_r8,  0.0001141_r8,  0.0001207_r8,  0.0001268_r8,  &
       0.0001326_r8,  0.0001380_r8,  0.0001431_r8,  0.0001480_r8,  0.0001526_r8,  &
       0.0001570_r8,  0.0001612_r8,  0.0001653_r8,  0.0001692_r8,  0.0001729_r8,  &
       0.0001766_r8,  0.0001801_r8,  0.0001835_r8,  0.0001868_r8,  0.0001900_r8,  &
       0.0001931_r8,  0.0001961_r8,  0.0001991_r8,  0.0002019_r8,  0.0002048_r8,  &
       0.0002075_r8,  0.0002102_r8,  0.0002129_r8,  0.0002154_r8,  0.0002180_r8,  &
       0.0002205_r8,  0.0002229_r8,  0.0002253_r8,  0.0002277_r8,  0.0002300_r8,  &
       0.0002323_r8,  0.0002346_r8,  0.0002368_r8,  0.0002390_r8,  0.0002411_r8,  &
       0.0002432_r8,  0.0002453_r8,  0.0002474_r8,  0.0002494_r8,  0.0002515_r8,  &
       0.0002535_r8,  0.0002554_r8,  0.0002574_r8,  0.0002593_r8,  0.0002612_r8,  &
       0.0002631_r8,  0.0002649_r8,  0.0002668_r8,  0.0002686_r8,  0.0002704_r8,  &
       0.0002722_r8,  0.0002740_r8,  0.0002757_r8,  0.0002775_r8,  0.0002792_r8,  &
       0.0002809_r8,  0.0002826_r8, &
       0.0001063_r8,  0.0001148_r8,  0.0001226_r8,  0.0001297_r8,  0.0001363_r8,  &
       0.0001425_r8,  0.0001483_r8,  0.0001538_r8,  0.0001590_r8,  0.0001639_r8,  &
       0.0001687_r8,  0.0001732_r8,  0.0001775_r8,  0.0001817_r8,  0.0001857_r8,  &
       0.0001896_r8,  0.0001933_r8,  0.0001970_r8,  0.0002005_r8,  0.0002039_r8,  &
       0.0002073_r8,  0.0002105_r8,  0.0002137_r8,  0.0002168_r8,  0.0002198_r8,  &
       0.0002228_r8,  0.0002257_r8,  0.0002286_r8,  0.0002314_r8,  0.0002341_r8,  &
       0.0002368_r8,  0.0002394_r8,  0.0002420_r8,  0.0002446_r8,  0.0002471_r8,  &
       0.0002496_r8,  0.0002520_r8,  0.0002544_r8,  0.0002568_r8,  0.0002591_r8,  &
       0.0002615_r8,  0.0002637_r8,  0.0002660_r8,  0.0002682_r8,  0.0002704_r8,  &
       0.0002726_r8,  0.0002747_r8,  0.0002768_r8,  0.0002789_r8,  0.0002810_r8,  &
       0.0002831_r8,  0.0002851_r8,  0.0002871_r8,  0.0002891_r8,  0.0002911_r8,  &
       0.0002930_r8,  0.0002950_r8,  0.0002969_r8,  0.0002988_r8,  0.0003007_r8,  &
       0.0003025_r8,  0.0003044_r8, &
       0.0001141_r8,  0.0001233_r8,  0.0001316_r8,  0.0001393_r8,  0.0001464_r8, &
       0.0001531_r8,  0.0001593_r8,  0.0001652_r8,  0.0001707_r8,  0.0001760_r8, &
       0.0001811_r8,  0.0001859_r8,  0.0001905_r8,  0.0001950_r8,  0.0001993_r8, &
       0.0002035_r8,  0.0002075_r8,  0.0002114_r8,  0.0002152_r8,  0.0002189_r8, &
       0.0002225_r8,  0.0002260_r8,  0.0002294_r8,  0.0002328_r8,  0.0002360_r8, &
       0.0002393_r8,  0.0002424_r8,  0.0002455_r8,  0.0002485_r8,  0.0002515_r8, &
       0.0002544_r8,  0.0002573_r8,  0.0002601_r8,  0.0002629_r8,  0.0002656_r8, &
       0.0002683_r8,  0.0002709_r8,  0.0002736_r8,  0.0002762_r8,  0.0002787_r8, &
       0.0002812_r8,  0.0002837_r8,  0.0002862_r8,  0.0002886_r8,  0.0002910_r8, &
       0.0002934_r8,  0.0002957_r8,  0.0002980_r8,  0.0003003_r8,  0.0003026_r8, &
       0.0003048_r8,  0.0003071_r8,  0.0003093_r8,  0.0003114_r8,  0.0003136_r8, &
       0.0003157_r8,  0.0003179_r8,  0.0003200_r8,  0.0003221_r8,  0.0003241_r8, &
       0.0003262_r8,  0.0003282_r8, &
       0.0001224_r8,  0.0001323_r8,  0.0001413_r8,  0.0001496_r8,  0.0001572_r8, &
       0.0001643_r8,  0.0001709_r8,  0.0001772_r8,  0.0001832_r8,  0.0001888_r8, &
       0.0001943_r8,  0.0001994_r8,  0.0002044_r8,  0.0002092_r8,  0.0002138_r8, &
       0.0002183_r8,  0.0002226_r8,  0.0002269_r8,  0.0002309_r8,  0.0002349_r8, &
       0.0002388_r8,  0.0002426_r8,  0.0002463_r8,  0.0002499_r8,  0.0002535_r8, &
       0.0002570_r8,  0.0002604_r8,  0.0002637_r8,  0.0002670_r8,  0.0002702_r8, &
       0.0002734_r8,  0.0002765_r8,  0.0002796_r8,  0.0002826_r8,  0.0002856_r8, &
       0.0002886_r8,  0.0002915_r8,  0.0002943_r8,  0.0002972_r8,  0.0002999_r8, &
       0.0003027_r8,  0.0003054_r8,  0.0003081_r8,  0.0003108_r8,  0.0003134_r8, &
       0.0003160_r8,  0.0003185_r8,  0.0003211_r8,  0.0003236_r8,  0.0003261_r8, &
       0.0003286_r8,  0.0003310_r8,  0.0003334_r8,  0.0003358_r8,  0.0003382_r8, &
       0.0003405_r8,  0.0003428_r8,  0.0003451_r8,  0.0003474_r8,  0.0003497_r8, &
       0.0003519_r8,  0.0003542_r8, &
       0.0001312_r8,  0.0001419_r8,  0.0001515_r8,  0.0001603_r8,  0.0001685_r8, &
       0.0001761_r8,  0.0001832_r8,  0.0001899_r8,  0.0001963_r8,  0.0002024_r8, &
       0.0002082_r8,  0.0002138_r8,  0.0002191_r8,  0.0002243_r8,  0.0002292_r8, &
       0.0002341_r8,  0.0002387_r8,  0.0002433_r8,  0.0002477_r8,  0.0002520_r8, &
       0.0002562_r8,  0.0002603_r8,  0.0002644_r8,  0.0002683_r8,  0.0002722_r8, &
       0.0002759_r8,  0.0002796_r8,  0.0002833_r8,  0.0002869_r8,  0.0002904_r8, &
       0.0002939_r8,  0.0002973_r8,  0.0003006_r8,  0.0003039_r8,  0.0003072_r8, &
       0.0003104_r8,  0.0003136_r8,  0.0003167_r8,  0.0003198_r8,  0.0003228_r8, &
       0.0003259_r8,  0.0003288_r8,  0.0003318_r8,  0.0003347_r8,  0.0003376_r8, &
       0.0003404_r8,  0.0003432_r8,  0.0003460_r8,  0.0003487_r8,  0.0003515_r8, &
       0.0003542_r8,  0.0003568_r8,  0.0003595_r8,  0.0003621_r8,  0.0003647_r8, &
       0.0003673_r8,  0.0003698_r8,  0.0003724_r8,  0.0003749_r8,  0.0003773_r8, &
       0.0003798_r8,  0.0003822_r8, &
       0.0001406_r8,  0.0001520_r8,  0.0001623_r8,  0.0001718_r8,  0.0001805_r8, &
       0.0001886_r8,  0.0001963_r8,  0.0002035_r8,  0.0002103_r8,  0.0002168_r8, &
       0.0002231_r8,  0.0002291_r8,  0.0002348_r8,  0.0002404_r8,  0.0002458_r8, &
       0.0002510_r8,  0.0002561_r8,  0.0002610_r8,  0.0002658_r8,  0.0002705_r8, &
       0.0002750_r8,  0.0002795_r8,  0.0002839_r8,  0.0002882_r8,  0.0002924_r8, &
       0.0002965_r8,  0.0003005_r8,  0.0003045_r8,  0.0003084_r8,  0.0003123_r8, &
       0.0003161_r8,  0.0003198_r8,  0.0003235_r8,  0.0003271_r8,  0.0003307_r8, &
       0.0003342_r8,  0.0003376_r8,  0.0003411_r8,  0.0003445_r8,  0.0003478_r8, &
       0.0003511_r8,  0.0003544_r8,  0.0003576_r8,  0.0003608_r8,  0.0003639_r8, &
       0.0003670_r8,  0.0003701_r8,  0.0003731_r8,  0.0003762_r8,  0.0003791_r8, &
       0.0003821_r8,  0.0003850_r8,  0.0003879_r8,  0.0003908_r8,  0.0003936_r8, &
       0.0003965_r8,  0.0003992_r8,  0.0004020_r8,  0.0004047_r8,  0.0004075_r8, &
       0.0004102_r8,  0.0004128_r8, &
       0.0001506_r8,  0.0001628_r8,  0.0001739_r8,  0.0001840_r8,  0.0001934_r8, &
       0.0002021_r8,  0.0002103_r8,  0.0002180_r8,  0.0002254_r8,  0.0002324_r8, &
       0.0002391_r8,  0.0002456_r8,  0.0002518_r8,  0.0002579_r8,  0.0002637_r8, &
       0.0002694_r8,  0.0002749_r8,  0.0002802_r8,  0.0002854_r8,  0.0002905_r8, &
       0.0002955_r8,  0.0003004_r8,  0.0003052_r8,  0.0003099_r8,  0.0003145_r8, &
       0.0003190_r8,  0.0003234_r8,  0.0003278_r8,  0.0003320_r8,  0.0003362_r8, &
       0.0003404_r8,  0.0003445_r8,  0.0003485_r8,  0.0003524_r8,  0.0003564_r8, &
       0.0003602_r8,  0.0003640_r8,  0.0003678_r8,  0.0003715_r8,  0.0003751_r8, &
       0.0003787_r8,  0.0003823_r8,  0.0003858_r8,  0.0003893_r8,  0.0003928_r8, &
       0.0003962_r8,  0.0003995_r8,  0.0004029_r8,  0.0004062_r8,  0.0004094_r8, &
       0.0004127_r8,  0.0004159_r8,  0.0004190_r8,  0.0004222_r8,  0.0004253_r8, &
       0.0004283_r8,  0.0004314_r8,  0.0004344_r8,  0.0004374_r8,  0.0004404_r8, &
       0.0004433_r8,  0.0004462_r8, &
       0.0001613_r8,  0.0001744_r8,  0.0001863_r8,  0.0001971_r8,  0.0002072_r8, &
       0.0002165_r8,  0.0002254_r8,  0.0002337_r8,  0.0002417_r8,  0.0002493_r8, &
       0.0002565_r8,  0.0002636_r8,  0.0002703_r8,  0.0002769_r8,  0.0002832_r8, &
       0.0002894_r8,  0.0002954_r8,  0.0003013_r8,  0.0003070_r8,  0.0003125_r8, &
       0.0003180_r8,  0.0003233_r8,  0.0003286_r8,  0.0003337_r8,  0.0003387_r8, &
       0.0003436_r8,  0.0003485_r8,  0.0003533_r8,  0.0003579_r8,  0.0003625_r8, &
       0.0003671_r8,  0.0003716_r8,  0.0003760_r8,  0.0003803_r8,  0.0003846_r8, &
       0.0003888_r8,  0.0003929_r8,  0.0003971_r8,  0.0004011_r8,  0.0004051_r8, &
       0.0004091_r8,  0.0004130_r8,  0.0004168_r8,  0.0004206_r8,  0.0004244_r8, &
       0.0004281_r8,  0.0004318_r8,  0.0004354_r8,  0.0004390_r8,  0.0004426_r8, &
       0.0004461_r8,  0.0004496_r8,  0.0004530_r8,  0.0004565_r8,  0.0004598_r8, &
       0.0004632_r8,  0.0004665_r8,  0.0004698_r8,  0.0004730_r8,  0.0004763_r8, &
       0.0004795_r8,  0.0004826_r8, &
       0.0001728_r8,  0.0001868_r8,  0.0001996_r8,  0.0002112_r8,  0.0002220_r8, &
       0.0002321_r8,  0.0002417_r8,  0.0002507_r8,  0.0002593_r8,  0.0002676_r8, &
       0.0002755_r8,  0.0002831_r8,  0.0002905_r8,  0.0002977_r8,  0.0003046_r8, &
       0.0003113_r8,  0.0003179_r8,  0.0003243_r8,  0.0003305_r8,  0.0003366_r8, &
       0.0003426_r8,  0.0003484_r8,  0.0003542_r8,  0.0003598_r8,  0.0003653_r8, &
       0.0003707_r8,  0.0003760_r8,  0.0003812_r8,  0.0003863_r8,  0.0003914_r8, &
       0.0003963_r8,  0.0004012_r8,  0.0004060_r8,  0.0004108_r8,  0.0004154_r8, &
       0.0004201_r8,  0.0004246_r8,  0.0004291_r8,  0.0004335_r8,  0.0004379_r8, &
       0.0004422_r8,  0.0004464_r8,  0.0004506_r8,  0.0004548_r8,  0.0004589_r8, &
       0.0004629_r8,  0.0004669_r8,  0.0004709_r8,  0.0004748_r8,  0.0004787_r8, &
       0.0004825_r8,  0.0004863_r8,  0.0004900_r8,  0.0004937_r8,  0.0004974_r8, &
       0.0005010_r8,  0.0005046_r8,  0.0005082_r8,  0.0005117_r8,  0.0005152_r8, &
       0.0005187_r8,  0.0005221_r8, &
       0.0001851_r8,  0.0002003_r8,  0.0002139_r8,  0.0002265_r8,  0.0002382_r8, &
       0.0002491_r8,  0.0002595_r8,  0.0002693_r8,  0.0002787_r8,  0.0002877_r8, &
       0.0002963_r8,  0.0003047_r8,  0.0003127_r8,  0.0003205_r8,  0.0003281_r8, &
       0.0003355_r8,  0.0003427_r8,  0.0003497_r8,  0.0003565_r8,  0.0003632_r8, &
       0.0003697_r8,  0.0003761_r8,  0.0003824_r8,  0.0003885_r8,  0.0003945_r8, &
       0.0004004_r8,  0.0004062_r8,  0.0004119_r8,  0.0004175_r8,  0.0004230_r8, &
       0.0004285_r8,  0.0004338_r8,  0.0004390_r8,  0.0004442_r8,  0.0004493_r8, &
       0.0004543_r8,  0.0004593_r8,  0.0004641_r8,  0.0004689_r8,  0.0004737_r8, &
       0.0004784_r8,  0.0004830_r8,  0.0004875_r8,  0.0004920_r8,  0.0004965_r8, &
       0.0005009_r8,  0.0005052_r8,  0.0005095_r8,  0.0005138_r8,  0.0005179_r8, &
       0.0005221_r8,  0.0005262_r8,  0.0005302_r8,  0.0005342_r8,  0.0005268_r8, &
       0.0005421_r8,  0.0005460_r8,  0.0005499_r8,  0.0005537_r8,  0.0005574_r8, &
       0.0005612_r8,  0.0005649_r8, &
       0.0001985_r8,  0.0002149_r8,  0.0002297_r8,  0.0002433_r8,  0.0002559_r8, &
       0.0002679_r8,  0.0002791_r8,  0.0002898_r8,  0.0003001_r8,  0.0003099_r8, &
       0.0003193_r8,  0.0003285_r8,  0.0003373_r8,  0.0003459_r8,  0.0003542_r8, &
       0.0003622_r8,  0.0003701_r8,  0.0003778_r8,  0.0003853_r8,  0.0003926_r8, &
       0.0003997_r8,  0.0004067_r8,  0.0004135_r8,  0.0004202_r8,  0.0004268_r8, &
       0.0004333_r8,  0.0004396_r8,  0.0004458_r8,  0.0004519_r8,  0.0004579_r8, &
       0.0004638_r8,  0.0004696_r8,  0.0004753_r8,  0.0004809_r8,  0.0004864_r8, &
       0.0004919_r8,  0.0004972_r8,  0.0005025_r8,  0.0005077_r8,  0.0005129_r8, &
       0.0005179_r8,  0.0005229_r8,  0.0005279_r8,  0.0005327_r8,  0.0005375_r8, &
       0.0005423_r8,  0.0005470_r8,  0.0005516_r8,  0.0005562_r8,  0.0005607_r8, &
       0.0005652_r8,  0.0005696_r8,  0.0005739_r8,  0.0005782_r8,  0.0005825_r8, &
       0.0005867_r8,  0.0005909_r8,  0.0005951_r8,  0.0005991_r8,  0.0006032_r8, &
       0.0006072_r8,  0.0006112_r8, &
       0.0002132_r8,  0.0002309_r8,  0.0002469_r8,  0.0002617_r8,  0.0002755_r8, &
       0.0002885_r8,  0.0003008_r8,  0.0003125_r8,  0.0003237_r8,  0.0003345_r8, &
       0.0003449_r8,  0.0003549_r8,  0.0003645_r8,  0.0003739_r8,  0.0003830_r8, &
       0.0003918_r8,  0.0004004_r8,  0.0004088_r8,  0.0004170_r8,  0.0004250_r8, &
       0.0004328_r8,  0.0004404_r8,  0.0004478_r8,  0.0004551_r8,  0.0004623_r8, &
       0.0004693_r8,  0.0004762_r8,  0.0004829_r8,  0.0004895_r8,  0.0004960_r8, &
       0.0005024_r8,  0.0005087_r8,  0.0005149_r8,  0.0005210_r8,  0.0005269_r8, &
       0.0005328_r8,  0.0005272_r8,  0.0005443_r8,  0.0005500_r8,  0.0005555_r8, &
       0.0005609_r8,  0.0005663_r8,  0.0005717_r8,  0.0005769_r8,  0.0005821_r8, &
       0.0005872_r8,  0.0005922_r8,  0.0005972_r8,  0.0006021_r8,  0.0006070_r8, &
       0.0006118_r8,  0.0006165_r8,  0.0006212_r8,  0.0006258_r8,  0.0006304_r8, &
       0.0006349_r8,  0.0006394_r8,  0.0006438_r8,  0.0006482_r8,  0.0006525_r8, &
       0.0006568_r8,  0.0006611_r8, &
       0.0002293_r8,  0.0002485_r8,  0.0002660_r8,  0.0002821_r8,  0.0002972_r8, &
       0.0003114_r8,  0.0003249_r8,  0.0003377_r8,  0.0003500_r8,  0.0003618_r8, &
       0.0003732_r8,  0.0003841_r8,  0.0003947_r8,  0.0004049_r8,  0.0004149_r8, &
       0.0004245_r8,  0.0004339_r8,  0.0004430_r8,  0.0004520_r8,  0.0004606_r8, &
       0.0004691_r8,  0.0004774_r8,  0.0004855_r8,  0.0004934_r8,  0.0005012_r8, &
       0.0005087_r8,  0.0005162_r8,  0.0005235_r8,  0.0005306_r8,  0.0005377_r8, &
       0.0005446_r8,  0.0005513_r8,  0.0005580_r8,  0.0005646_r8,  0.0005710_r8, &
       0.0005773_r8,  0.0005836_r8,  0.0005897_r8,  0.0005957_r8,  0.0006017_r8, &
       0.0006076_r8,  0.0006134_r8,  0.0006191_r8,  0.0006247_r8,  0.0006302_r8, &
       0.0006357_r8,  0.0006411_r8,  0.0006464_r8,  0.0006517_r8,  0.0006569_r8, &
       0.0006620_r8,  0.0006671_r8,  0.0006721_r8,  0.0006771_r8,  0.0006820_r8, &
       0.0006868_r8,  0.0006916_r8,  0.0006963_r8,  0.0007010_r8,  0.0007056_r8, &
       0.0007102_r8,  0.0007147_r8, &
       0.0002471_r8,  0.0002680_r8,  0.0002871_r8,  0.0003048_r8,  0.0003214_r8, &
       0.0003369_r8,  0.0003517_r8,  0.0003658_r8,  0.0003792_r8,  0.0003921_r8, &
       0.0004045_r8,  0.0004165_r8,  0.0004281_r8,  0.0004392_r8,  0.0004501_r8, &
       0.0004606_r8,  0.0004708_r8,  0.0004807_r8,  0.0004903_r8,  0.0004998_r8, &
       0.0005089_r8,  0.0005179_r8,  0.0005267_r8,  0.0005352_r8,  0.0005436_r8, &
       0.0005518_r8,  0.0005598_r8,  0.0005677_r8,  0.0005754_r8,  0.0005829_r8, &
       0.0005903_r8,  0.0005976_r8,  0.0006048_r8,  0.0006118_r8,  0.0006187_r8, &
       0.0006255_r8,  0.0006322_r8,  0.0006388_r8,  0.0006453_r8,  0.0006516_r8, &
       0.0006579_r8,  0.0006641_r8,  0.0006702_r8,  0.0006762_r8,  0.0006822_r8, &
       0.0006880_r8,  0.0006938_r8,  0.0006995_r8,  0.0007051_r8,  0.0007106_r8, &
       0.0007161_r8,  0.0007215_r8,  0.0007269_r8,  0.0007321_r8,  0.0007374_r8, &
       0.0007425_r8,  0.0007476_r8,  0.0007527_r8,  0.0007576_r8,  0.0007626_r8, &
       0.0007674_r8,  0.0007723_r8, &
       0.0002669_r8,  0.0002898_r8,  0.0003107_r8,  0.0003300_r8,  0.0003482_r8, &
       0.0003653_r8,  0.0003815_r8,  0.0003969_r8,  0.0004116_r8,  0.0004257_r8, &
       0.0004392_r8,  0.0004522_r8,  0.0004648_r8,  0.0004769_r8,  0.0004887_r8, &
       0.0005001_r8,  0.0005111_r8,  0.0005218_r8,  0.0005323_r8,  0.0005425_r8, &
       0.0005524_r8,  0.0005620_r8,  0.0005714_r8,  0.0005807_r8,  0.0005897_r8, &
       0.0005985_r8,  0.0006071_r8,  0.0006155_r8,  0.0006238_r8,  0.0006319_r8, &
       0.0006398_r8,  0.0006476_r8,  0.0006553_r8,  0.0006628_r8,  0.0006702_r8, &
       0.0006775_r8,  0.0006846_r8,  0.0006917_r8,  0.0006986_r8,  0.0007054_r8, &
       0.0007121_r8,  0.0007187_r8,  0.0007252_r8,  0.0007316_r8,  0.0007379_r8, &
       0.0007442_r8,  0.0007503_r8,  0.0007564_r8,  0.0007624_r8,  0.0007683_r8, &
       0.0007741_r8,  0.0007798_r8,  0.0007855_r8,  0.0007911_r8,  0.0007967_r8, &
       0.0008022_r8,  0.0008076_r8,  0.0008129_r8,  0.0008182_r8,  0.0008235_r8, &
       0.0008286_r8,  0.0008337_r8, &
       0.0002889_r8,  0.0003140_r8,  0.0003369_r8,  0.0003582_r8,  0.0003780_r8, &
       0.0003967_r8,  0.0004144_r8,  0.0004312_r8,  0.0004473_r8,  0.0004626_r8, &
       0.0004773_r8,  0.0004914_r8,  0.0005050_r8,  0.0005182_r8,  0.0005309_r8, &
       0.0005432_r8,  0.0005551_r8,  0.0005666_r8,  0.0005779_r8,  0.0005888_r8, &
       0.0005995_r8,  0.0006098_r8,  0.0006200_r8,  0.0006298_r8,  0.0006395_r8, &
       0.0006489_r8,  0.0006581_r8,  0.0006672_r8,  0.0006760_r8,  0.0006847_r8, &
       0.0006932_r8,  0.0007015_r8,  0.0007097_r8,  0.0007177_r8,  0.0007256_r8, &
       0.0007333_r8,  0.0007409_r8,  0.0007484_r8,  0.0007558_r8,  0.0007630_r8, &
       0.0007702_r8,  0.0007772_r8,  0.0007841_r8,  0.0007909_r8,  0.0007976_r8, &
       0.0008043_r8,  0.0008108_r8,  0.0008172_r8,  0.0008236_r8,  0.0008298_r8, &
       0.0008360_r8,  0.0008421_r8,  0.0008481_r8,  0.0008541_r8,  0.0008600_r8, &
       0.0008658_r8,  0.0008715_r8,  0.0008772_r8,  0.0008828_r8,  0.0008883_r8, &
       0.0008938_r8,  0.0008992_r8, &
       0.0003135_r8,  0.0003410_r8,  0.0003662_r8,  0.0003895_r8,  0.0004112_r8, &
       0.0004316_r8,  0.0004509_r8,  0.0004692_r8,  0.0004866_r8,  0.0005032_r8, &
       0.0005191_r8,  0.0005344_r8,  0.0005491_r8,  0.0005632_r8,  0.0005769_r8, &
       0.0005901_r8,  0.0006029_r8,  0.0006153_r8,  0.0006274_r8,  0.0006391_r8, &
       0.0006505_r8,  0.0006616_r8,  0.0006725_r8,  0.0006830_r8,  0.0006933_r8, &
       0.0007034_r8,  0.0007132_r8,  0.0007229_r8,  0.0007323_r8,  0.0007415_r8, &
       0.0007506_r8,  0.0007595_r8,  0.0007682_r8,  0.0007767_r8,  0.0007851_r8, &
       0.0007933_r8,  0.0008014_r8,  0.0008093_r8,  0.0008172_r8,  0.0008249_r8, &
       0.0008324_r8,  0.0008399_r8,  0.0008472_r8,  0.0008544_r8,  0.0008615_r8, &
       0.0008685_r8,  0.0008755_r8,  0.0008823_r8,  0.0008890_r8,  0.0008956_r8, &
       0.0009021_r8,  0.0009086_r8,  0.0009149_r8,  0.0009212_r8,  0.0009274_r8, &
       0.0009335_r8,  0.0009396_r8,  0.0009455_r8,  0.0009514_r8,  0.0009573_r8, &
       0.0009630_r8,  0.0009687_r8, &
       0.0003409_r8,  0.0003711_r8,  0.0003987_r8,  0.0004241_r8,  0.0004478_r8, &
       0.0004700_r8,  0.0004909_r8,  0.0005107_r8,  0.0005295_r8,  0.0005474_r8, &
       0.0005645_r8,  0.0005810_r8,  0.0005968_r8,  0.0006120_r8,  0.0006267_r8, &
       0.0006408_r8,  0.0006545_r8,  0.0006678_r8,  0.0006807_r8,  0.0006932_r8, &
       0.0007054_r8,  0.0007173_r8,  0.0007288_r8,  0.0007401_r8,  0.0007510_r8, &
       0.0007618_r8,  0.0007722_r8,  0.0007825_r8,  0.0007925_r8,  0.0008023_r8, &
       0.0008119_r8,  0.0008213_r8,  0.0008306_r8,  0.0008396_r8,  0.0008485_r8, &
       0.0008572_r8,  0.0008658_r8,  0.0008742_r8,  0.0008824_r8,  0.0008906_r8, &
       0.0008986_r8,  0.0009064_r8,  0.0009142_r8,  0.0009218_r8,  0.0009293_r8, &
       0.0009367_r8,  0.0009439_r8,  0.0009511_r8,  0.0009582_r8,  0.0009652_r8, &
       0.0009720_r8,  0.0009788_r8,  0.0009855_r8,  0.0009921_r8,  0.0009986_r8, &
       0.0010050_r8,  0.0010113_r8,  0.0010176_r8,  0.0010238_r8,  0.0010299_r8, &
       0.0010359_r8,  0.0010419_r8, &
       0.0003715_r8,  0.0004046_r8,  0.0004346_r8,  0.0004623_r8,  0.0004880_r8, &
       0.0005120_r8,  0.0005346_r8,  0.0005560_r8,  0.0005762_r8,  0.0005955_r8, &
       0.0006139_r8,  0.0006315_r8,  0.0006485_r8,  0.0006648_r8,  0.0006804_r8, &
       0.0006956_r8,  0.0007102_r8,  0.0007244_r8,  0.0007382_r8,  0.0007515_r8, &
       0.0007645_r8,  0.0007771_r8,  0.0007893_r8,  0.0008012_r8,  0.0008129_r8, &
       0.0008243_r8,  0.0008354_r8,  0.0008463_r8,  0.0008569_r8,  0.0008673_r8, &
       0.0008774_r8,  0.0008874_r8,  0.0008971_r8,  0.0009067_r8,  0.0009160_r8, &
       0.0009252_r8,  0.0009342_r8,  0.0009431_r8,  0.0009518_r8,  0.0009604_r8, &
       0.0009688_r8,  0.0009770_r8,  0.0009851_r8,  0.0009931_r8,  0.0010010_r8, &
       0.0010088_r8,  0.0010164_r8,  0.0010239_r8,  0.0010313_r8,  0.0010386_r8, &
       0.0010458_r8,  0.0010529_r8,  0.0010598_r8,  0.0010667_r8,  0.0010735_r8, &
       0.0010802_r8,  0.0010869_r8,  0.0010934_r8,  0.0010998_r8,  0.0011062_r8, &
       0.0011125_r8,  0.0011187_r8, &
       0.0004055_r8,  0.0004415_r8,  0.0004742_r8,  0.0005042_r8,  0.0005320_r8, &
       0.0005579_r8,  0.0005822_r8,  0.0006052_r8,  0.0006269_r8,  0.0006476_r8, &
       0.0006673_r8,  0.0006862_r8,  0.0007042_r8,  0.0007216_r8,  0.0007383_r8, &
       0.0007545_r8,  0.0007701_r8,  0.0007851_r8,  0.0007997_r8,  0.0008139_r8, &
       0.0008276_r8,  0.0008410_r8,  0.0008540_r8,  0.0008666_r8,  0.0008789_r8, &
       0.0008910_r8,  0.0009027_r8,  0.0009141_r8,  0.0009253_r8,  0.0009362_r8, &
       0.0009469_r8,  0.0009574_r8,  0.0009676_r8,  0.0009777_r8,  0.0009875_r8, &
       0.0009972_r8,  0.0010066_r8,  0.0010159_r8,  0.0010250_r8,  0.0010339_r8, &
       0.0010427_r8,  0.0010514_r8,  0.0010599_r8,  0.0010682_r8,  0.0010764_r8, &
       0.0010845_r8,  0.0010925_r8,  0.0011003_r8,  0.0011080_r8,  0.0011156_r8, &
       0.0011231_r8,  0.0011304_r8,  0.0011377_r8,  0.0011449_r8,  0.0011519_r8, &
       0.0011589_r8,  0.0011658_r8,  0.0011726_r8,  0.0011793_r8,  0.0011859_r8, &
       0.0011924_r8,  0.0011989_r8, &
       0.0004429_r8,  0.0004821_r8,  0.0005175_r8,  0.0005499_r8,  0.0005798_r8, &
       0.0006076_r8,  0.0006337_r8,  0.0006583_r8,  0.0006816_r8,  0.0007037_r8, &
       0.0007247_r8,  0.0007448_r8,  0.0007640_r8,  0.0007825_r8,  0.0008003_r8, &
       0.0008174_r8,  0.0008339_r8,  0.0008499_r8,  0.0008653_r8,  0.0008803_r8, &
       0.0008948_r8,  0.0009089_r8,  0.0009226_r8,  0.0009359_r8,  0.0009488_r8, &
       0.0009615_r8,  0.0009738_r8,  0.0009858_r8,  0.0009975_r8,  0.0010089_r8, &
       0.0010201_r8,  0.0010311_r8,  0.0010418_r8,  0.0010523_r8,  0.0010625_r8, &
       0.0010726_r8,  0.0010825_r8,  0.0010921_r8,  0.0011016_r8,  0.0011109_r8, &
       0.0011201_r8,  0.0011291_r8,  0.0011379_r8,  0.0011466_r8,  0.0011551_r8, &
       0.0011635_r8,  0.0011718_r8,  0.0011799_r8,  0.0011879_r8,  0.0011958_r8, &
       0.0012035_r8,  0.0012112_r8,  0.0012187_r8,  0.0012261_r8,  0.0012335_r8, &
       0.0012407_r8,  0.0012478_r8,  0.0012548_r8,  0.0012618_r8,  0.0012686_r8, &
       0.0012754_r8,  0.0012821_r8, &
       0.0004840_r8,  0.0005264_r8,  0.0005646_r8,  0.0005994_r8,  0.0006316_r8, &
       0.0006614_r8,  0.0006893_r8,  0.0007155_r8,  0.0007403_r8,  0.0007638_r8, &
       0.0007862_r8,  0.0008075_r8,  0.0008279_r8,  0.0008475_r8,  0.0008663_r8, &
       0.0008844_r8,  0.0009018_r8,  0.0009186_r8,  0.0009349_r8,  0.0009506_r8, &
       0.0009658_r8,  0.0009806_r8,  0.0009950_r8,  0.0010089_r8,  0.0010225_r8, &
       0.0010356_r8,  0.0010485_r8,  0.0010610_r8,  0.0010733_r8,  0.0010852_r8, &
       0.0010968_r8,  0.0011082_r8,  0.0011194_r8,  0.0011303_r8,  0.0011409_r8, &
       0.0011514_r8,  0.0011616_r8,  0.0011717_r8,  0.0011815_r8,  0.0011912_r8, &
       0.0012007_r8,  0.0012100_r8,  0.0012191_r8,  0.0012281_r8,  0.0012370_r8, &
       0.0012457_r8,  0.0012542_r8,  0.0012627_r8,  0.0012709_r8,  0.0012791_r8, &
       0.0012871_r8,  0.0012951_r8,  0.0013029_r8,  0.0013106_r8,  0.0013181_r8, &
       0.0013256_r8,  0.0013330_r8,  0.0013403_r8,  0.0013475_r8,  0.0013546_r8, &
       0.0013616_r8,  0.0013685_r8, &
       0.0005290_r8,  0.0005747_r8,  0.0006157_r8,  0.0006530_r8,  0.0006874_r8, &
       0.0007192_r8,  0.0007490_r8,  0.0007769_r8,  0.0008032_r8,  0.0008281_r8, &
       0.0008518_r8,  0.0008743_r8,  0.0008959_r8,  0.0009165_r8,  0.0009362_r8, &
       0.0009552_r8,  0.0009735_r8,  0.0009911_r8,  0.0010082_r8,  0.0010246_r8, &
       0.0010405_r8,  0.0010559_r8,  0.0010709_r8,  0.0010854_r8,  0.0010995_r8, &
       0.0011132_r8,  0.0011266_r8,  0.0011396_r8,  0.0011523_r8,  0.0011647_r8, &
       0.0011768_r8,  0.0011886_r8,  0.0012002_r8,  0.0012115_r8,  0.0012225_r8, &
       0.0012334_r8,  0.0012440_r8,  0.0012544_r8,  0.0012646_r8,  0.0012746_r8, &
       0.0012844_r8,  0.0012941_r8,  0.0013036_r8,  0.0013129_r8,  0.0013220_r8, &
       0.0013310_r8,  0.0013399_r8,  0.0013486_r8,  0.0013572_r8,  0.0013657_r8, &
       0.0013740_r8,  0.0013823_r8,  0.0013904_r8,  0.0013983_r8,  0.0014062_r8, &
       0.0014140_r8,  0.0014217_r8,  0.0014292_r8,  0.0014367_r8,  0.0014441_r8, &
       0.0014514_r8,  0.0014586_r8, &
       0.0005778_r8,  0.0006269_r8,  0.0006708_r8,  0.0007107_r8,  0.0007473_r8, &
       0.0007812_r8,  0.0008127_r8,  0.0008423_r8,  0.0008701_r8,  0.0008964_r8, &
       0.0009213_r8,  0.0009450_r8,  0.0009676_r8,  0.0009892_r8,  0.0010099_r8, &
       0.0010297_r8,  0.0010488_r8,  0.0010671_r8,  0.0010848_r8,  0.0011020_r8, &
       0.0011185_r8,  0.0011345_r8,  0.0011500_r8,  0.0011651_r8,  0.0011798_r8, &
       0.0011940_r8,  0.0012078_r8,  0.0012213_r8,  0.0012345_r8,  0.0012473_r8, &
       0.0012599_r8,  0.0012721_r8,  0.0012841_r8,  0.0012958_r8,  0.0013073_r8, &
       0.0013185_r8,  0.0013295_r8,  0.0013403_r8,  0.0013509_r8,  0.0013612_r8, &
       0.0013714_r8,  0.0013815_r8,  0.0013913_r8,  0.0014010_r8,  0.0014105_r8, &
       0.0014199_r8,  0.0014291_r8,  0.0014382_r8,  0.0014471_r8,  0.0014559_r8, &
       0.0014646_r8,  0.0014732_r8,  0.0014816_r8,  0.0014900_r8,  0.0014982_r8, &
       0.0015063_r8,  0.0015143_r8,  0.0015222_r8,  0.0015300_r8,  0.0015377_r8, &
       0.0015453_r8,  0.0015528_r8, &
       0.0006307_r8,  0.0006832_r8,  0.0007301_r8,  0.0007725_r8,  0.0008114_r8, &
       0.0008472_r8,  0.0008805_r8,  0.0009116_r8,  0.0009409_r8,  0.0009684_r8, &
       0.0009945_r8,  0.0010193_r8,  0.0010428_r8,  0.0010653_r8,  0.0010868_r8, &
       0.0011075_r8,  0.0011273_r8,  0.0011464_r8,  0.0011647_r8,  0.0011825_r8, &
       0.0011996_r8,  0.0012162_r8,  0.0012323_r8,  0.0012480_r8,  0.0012631_r8, &
       0.0012779_r8,  0.0012922_r8,  0.0013062_r8,  0.0013199_r8,  0.0013332_r8, &
       0.0013462_r8,  0.0013589_r8,  0.0013713_r8,  0.0013835_r8,  0.0013954_r8, &
       0.0014071_r8,  0.0014186_r8,  0.0014298_r8,  0.0014408_r8,  0.0014516_r8, &
       0.0014623_r8,  0.0014727_r8,  0.0014830_r8,  0.0014931_r8,  0.0015030_r8, &
       0.0015128_r8,  0.0015225_r8,  0.0015319_r8,  0.0015413_r8,  0.0015505_r8, &
       0.0015596_r8,  0.0015686_r8,  0.0015774_r8,  0.0015862_r8,  0.0015948_r8, &
       0.0016033_r8,  0.0016117_r8,  0.0016200_r8,  0.0016282_r8,  0.0016363_r8, &
       0.0016443_r8,  0.0016522_r8, &
       0.0006876_r8,  0.0007436_r8,  0.0007934_r8,  0.0008383_r8,  0.0008793_r8, &
       0.0009170_r8,  0.0009520_r8,  0.0009846_r8,  0.0010150_r8,  0.0010439_r8, &
       0.0010710_r8,  0.0010968_r8,  0.0011213_r8,  0.0011446_r8,  0.0011669_r8, &
       0.0011883_r8,  0.0012089_r8,  0.0012287_r8,  0.0012477_r8,  0.0012661_r8, &
       0.0012839_r8,  0.0013011_r8,  0.0013178_r8,  0.0013340_r8,  0.0013498_r8, &
       0.0013651_r8,  0.0013800_r8,  0.0013946_r8,  0.0014088_r8,  0.0014227_r8, &
       0.0014362_r8,  0.0014495_r8,  0.0014624_r8,  0.0014751_r8,  0.0014876_r8, &
       0.0014998_r8,  0.0015117_r8,  0.0015235_r8,  0.0015350_r8,  0.0015464_r8, &
       0.0015575_r8,  0.0015685_r8,  0.0015792_r8,  0.0015899_r8,  0.0016003_r8, &
       0.0016106_r8,  0.0016207_r8,  0.0016307_r8,  0.0016405_r8,  0.0016502_r8, &
       0.0016598_r8,  0.0016693_r8,  0.0016786_r8,  0.0016878_r8,  0.0016969_r8, &
       0.0017059_r8,  0.0017147_r8,  0.0017235_r8,  0.0017322_r8,  0.0017407_r8, &
       0.0017492_r8,  0.0017575_r8, &
       0.0007485_r8,  0.0008080_r8,  0.0008606_r8,  0.0009079_r8,  0.0009509_r8, &
       0.0009904_r8,  0.0010269_r8,  0.0010608_r8,  0.0010926_r8,  0.0011225_r8, &
       0.0011507_r8,  0.0011774_r8,  0.0012028_r8,  0.0012270_r8,  0.0012501_r8, &
       0.0012723_r8,  0.0012936_r8,  0.0013142_r8,  0.0013339_r8,  0.0013531_r8, &
       0.0013716_r8,  0.0013895_r8,  0.0014069_r8,  0.0014238_r8,  0.0014402_r8, &
       0.0014562_r8,  0.0014718_r8,  0.0014870_r8,  0.0015019_r8,  0.0015164_r8, &
       0.0015306_r8,  0.0015445_r8,  0.0015581_r8,  0.0015714_r8,  0.0015845_r8, &
       0.0015973_r8,  0.0016099_r8,  0.0016223_r8,  0.0016344_r8,  0.0016464_r8, &
       0.0016581_r8,  0.0016697_r8,  0.0016811_r8,  0.0016922_r8,  0.0017033_r8, &
       0.0017141_r8,  0.0017249_r8,  0.0017354_r8,  0.0017458_r8,  0.0017561_r8, &
       0.0017662_r8,  0.0017762_r8,  0.0017861_r8,  0.0017959_r8,  0.0018055_r8, &
       0.0018150_r8,  0.0018244_r8,  0.0018337_r8,  0.0018429_r8,  0.0018520_r8, &
       0.0018610_r8,  0.0018698_r8, &
       0.0008135_r8,  0.0008762_r8,  0.0009315_r8,  0.0009811_r8,  0.0010259_r8, &
       0.0010670_r8,  0.0011050_r8,  0.0011402_r8,  0.0011732_r8,  0.0012042_r8, &
       0.0012334_r8,  0.0012611_r8,  0.0012874_r8,  0.0013126_r8,  0.0013366_r8, &
       0.0013597_r8,  0.0013819_r8,  0.0014033_r8,  0.0014239_r8,  0.0014439_r8, &
       0.0014632_r8,  0.0014820_r8,  0.0015002_r8,  0.0015179_r8,  0.0015351_r8, &
       0.0015519_r8,  0.0015683_r8,  0.0015843_r8,  0.0016000_r8,  0.0016153_r8, &
       0.0016302_r8,  0.0016449_r8,  0.0016592_r8,  0.0016733_r8,  0.0016871_r8, &
       0.0017007_r8,  0.0017140_r8,  0.0017271_r8,  0.0017400_r8,  0.0017526_r8, &
       0.0017651_r8,  0.0017773_r8,  0.0017894_r8,  0.0018012_r8,  0.0018129_r8, &
       0.0018245_r8,  0.0018358_r8,  0.0018471_r8,  0.0018581_r8,  0.0018690_r8, &
       0.0018798_r8,  0.0018904_r8,  0.0019009_r8,  0.0019113_r8,  0.0019215_r8, &
       0.0019316_r8,  0.0019416_r8,  0.0019515_r8,  0.0019613_r8,  0.0019709_r8, &
       0.0019805_r8,  0.0019899_r8, &
       0.0008823_r8,  0.0009481_r8,  0.0010059_r8,  0.0010575_r8,  0.0011042_r8, &
       0.0011468_r8,  0.0011862_r8,  0.0012227_r8,  0.0012569_r8,  0.0012891_r8, &
       0.0013195_r8,  0.0013483_r8,  0.0013757_r8,  0.0014020_r8,  0.0014271_r8, &
       0.0014512_r8,  0.0014744_r8,  0.0014968_r8,  0.0015185_r8,  0.0015395_r8, &
       0.0015598_r8,  0.0015795_r8,  0.0015987_r8,  0.0016174_r8,  0.0016356_r8, &
       0.0016533_r8,  0.0016707_r8,  0.0016876_r8,  0.0017041_r8,  0.0017203_r8, &
       0.0017362_r8,  0.0017517_r8,  0.0017669_r8,  0.0017819_r8,  0.0017965_r8, &
       0.0018109_r8,  0.0018251_r8,  0.0018390_r8,  0.0018527_r8,  0.0018661_r8, &
       0.0018793_r8,  0.0018924_r8,  0.0019052_r8,  0.0019178_r8,  0.0019303_r8, &
       0.0019425_r8,  0.0019546_r8,  0.0019665_r8,  0.0019783_r8,  0.0019899_r8, &
       0.0020014_r8,  0.0020127_r8,  0.0020238_r8,  0.0020348_r8,  0.0020457_r8, &
       0.0020565_r8,  0.0020671_r8,  0.0020776_r8,  0.0020880_r8,  0.0020983_r8, &
       0.0021084_r8,  0.0021185_r8, &
       0.0009546_r8,  0.0010233_r8,  0.0010834_r8,  0.0011370_r8,  0.0011854_r8, &
       0.0012297_r8,  0.0012705_r8,  0.0013085_r8,  0.0013441_r8,  0.0013776_r8, &
       0.0014093_r8,  0.0014395_r8,  0.0014682_r8,  0.0014957_r8,  0.0015221_r8, &
       0.0015475_r8,  0.0015720_r8,  0.0015956_r8,  0.0016185_r8,  0.0016406_r8, &
       0.0016621_r8,  0.0016830_r8,  0.0017033_r8,  0.0017231_r8,  0.0017424_r8, &
       0.0017613_r8,  0.0017797_r8,  0.0017976_r8,  0.0018152_r8,  0.0018324_r8, &
       0.0018493_r8,  0.0018658_r8,  0.0018820_r8,  0.0018979_r8,  0.0019135_r8, &
       0.0019288_r8,  0.0019439_r8,  0.0019587_r8,  0.0019732_r8,  0.0019875_r8, &
       0.0020016_r8,  0.0020155_r8,  0.0020291_r8,  0.0020425_r8,  0.0020558_r8, &
       0.0020688_r8,  0.0020817_r8,  0.0020944_r8,  0.0021069_r8,  0.0021192_r8, &
       0.0021314_r8,  0.0021434_r8,  0.0021095_r8,  0.0021670_r8,  0.0021786_r8, &
       0.0021900_r8,  0.0022013_r8,  0.0022125_r8,  0.0022235_r8,  0.0022344_r8, &
       0.0022452_r8,  0.0022558_r8, &
       0.0010302_r8,  0.0011016_r8,  0.0011640_r8,  0.0012195_r8,  0.0012698_r8, &
       0.0013159_r8,  0.0013584_r8,  0.0013981_r8,  0.0014354_r8,  0.0014705_r8, &
       0.0015038_r8,  0.0015355_r8,  0.0015658_r8,  0.0015948_r8,  0.0016227_r8, &
       0.0016496_r8,  0.0016755_r8,  0.0017005_r8,  0.0017248_r8,  0.0017483_r8, &
       0.0017712_r8,  0.0017934_r8,  0.0018150_r8,  0.0018360_r8,  0.0018566_r8, &
       0.0018766_r8,  0.0018962_r8,  0.0019153_r8,  0.0019340_r8,  0.0019524_r8, &
       0.0019703_r8,  0.0019879_r8,  0.0020052_r8,  0.0020221_r8,  0.0020387_r8, &
       0.0020550_r8,  0.0020710_r8,  0.0020867_r8,  0.0021022_r8,  0.0021174_r8, &
       0.0021324_r8,  0.0021471_r8,  0.0021159_r8,  0.0021759_r8,  0.0021900_r8, &
       0.0022039_r8,  0.0022175_r8,  0.0022310_r8,  0.0022443_r8,  0.0022574_r8, &
       0.0022703_r8,  0.0022830_r8,  0.0022956_r8,  0.0023080_r8,  0.0023203_r8, &
       0.0023324_r8,  0.0023443_r8,  0.0023562_r8,  0.0023678_r8,  0.0023794_r8, &
       0.0023908_r8,  0.0024020_r8, &
       0.0011087_r8,  0.0011828_r8,  0.0012476_r8,  0.0013054_r8,  0.0013579_r8, &
       0.0014060_r8,  0.0014506_r8,  0.0014923_r8,  0.0015315_r8,  0.0015685_r8, &
       0.0016038_r8,  0.0016373_r8,  0.0016694_r8,  0.0017002_r8,  0.0017298_r8, &
       0.0017583_r8,  0.0017859_r8,  0.0018125_r8,  0.0018384_r8,  0.0018634_r8, &
       0.0018877_r8,  0.0019114_r8,  0.0019344_r8,  0.0019568_r8,  0.0019787_r8, &
       0.0020000_r8,  0.0020209_r8,  0.0020412_r8,  0.0020612_r8,  0.0020807_r8, &
       0.0020998_r8,  0.0021185_r8,  0.0021368_r8,  0.0021090_r8,  0.0021725_r8, &
       0.0021898_r8,  0.0022068_r8,  0.0022235_r8,  0.0022399_r8,  0.0022561_r8, &
       0.0022720_r8,  0.0022876_r8,  0.0023029_r8,  0.0023181_r8,  0.0023330_r8, &
       0.0023477_r8,  0.0023621_r8,  0.0023764_r8,  0.0023904_r8,  0.0024042_r8, &
       0.0024179_r8,  0.0024314_r8,  0.0024446_r8,  0.0024577_r8,  0.0024707_r8, &
       0.0024834_r8,  0.0024960_r8,  0.0025085_r8,  0.0025208_r8,  0.0025329_r8, &
       0.0025449_r8,  0.0025568_r8, &
       0.0011902_r8,  0.0012672_r8,  0.0013347_r8,  0.0013952_r8,  0.0014502_r8,  &
       0.0015008_r8,  0.0015478_r8,  0.0015919_r8,  0.0016334_r8,  0.0016727_r8,  &
       0.0017101_r8,  0.0017457_r8,  0.0017799_r8,  0.0018126_r8,  0.0018442_r8,  &
       0.0018746_r8,  0.0019039_r8,  0.0019323_r8,  0.0019598_r8,  0.0019865_r8,  &
       0.0020124_r8,  0.0020376_r8,  0.0020621_r8,  0.0020859_r8,  0.0021092_r8,  &
       0.0021319_r8,  0.0021083_r8,  0.0021757_r8,  0.0021969_r8,  0.0022176_r8,  &
       0.0022379_r8,  0.0022577_r8,  0.0022772_r8,  0.0022962_r8,  0.0023149_r8,  &
       0.0023333_r8,  0.0023513_r8,  0.0023690_r8,  0.0023863_r8,  0.0024034_r8,  &
       0.0024202_r8,  0.0024366_r8,  0.0024529_r8,  0.0024688_r8,  0.0024845_r8,  &
       0.0025000_r8,  0.0025152_r8,  0.0025302_r8,  0.0025450_r8,  0.0025595_r8,  &
       0.0025739_r8,  0.0025880_r8,  0.0026019_r8,  0.0026157_r8,  0.0026293_r8,  &
       0.0026426_r8,  0.0026557_r8,  0.0026689_r8,  0.0026817_r8,  0.0026944_r8,  &
       0.0027070_r8,  0.0027194_r8, &
       0.0012749_r8,  0.0013552_r8,  0.0014259_r8,  0.0014895_r8,  0.0015475_r8, &
       0.0016010_r8,  0.0016509_r8,  0.0016977_r8,  0.0017418_r8,  0.0017836_r8, &
       0.0018234_r8,  0.0018614_r8,  0.0018978_r8,  0.0019328_r8,  0.0019664_r8, &
       0.0019987_r8,  0.0020300_r8,  0.0020602_r8,  0.0020895_r8,  0.0021179_r8, &
       0.0021454_r8,  0.0021722_r8,  0.0021982_r8,  0.0022235_r8,  0.0022482_r8, &
       0.0022723_r8,  0.0022958_r8,  0.0023187_r8,  0.0023411_r8,  0.0023630_r8, &
       0.0023844_r8,  0.0024054_r8,  0.0024259_r8,  0.0024460_r8,  0.0024658_r8, &
       0.0024851_r8,  0.0025040_r8,  0.0025226_r8,  0.0025409_r8,  0.0025588_r8, &
       0.0025764_r8,  0.0025937_r8,  0.0026107_r8,  0.0026275_r8,  0.0026439_r8, &
       0.0026601_r8,  0.0026760_r8,  0.0026917_r8,  0.0027071_r8,  0.0027223_r8, &
       0.0027373_r8,  0.0027520_r8,  0.0027665_r8,  0.0027808_r8,  0.0027950_r8, &
       0.0028089_r8,  0.0028226_r8,  0.0028361_r8,  0.0028495_r8,  0.0028627_r8, &
       0.0028757_r8,  0.0028885_r8, &
       0.0013631_r8,  0.0014474_r8,  0.0015220_r8,  0.0015892_r8,  0.0016507_r8, &
       0.0017076_r8,  0.0017607_r8,  0.0018105_r8,  0.0018575_r8,  0.0019021_r8, &
       0.0019445_r8,  0.0019850_r8,  0.0020238_r8,  0.0020610_r8,  0.0020967_r8, &
       0.0021312_r8,  0.0021186_r8,  0.0021965_r8,  0.0022276_r8,  0.0022577_r8, &
       0.0022868_r8,  0.0023152_r8,  0.0023427_r8,  0.0023695_r8,  0.0023956_r8, &
       0.0024210_r8,  0.0024457_r8,  0.0024699_r8,  0.0024935_r8,  0.0025165_r8, &
       0.0025390_r8,  0.0025610_r8,  0.0025826_r8,  0.0026037_r8,  0.0026243_r8, &
       0.0026445_r8,  0.0026643_r8,  0.0026838_r8,  0.0027028_r8,  0.0027215_r8, &
       0.0027399_r8,  0.0027579_r8,  0.0027756_r8,  0.0027930_r8,  0.0028101_r8, &
       0.0028269_r8,  0.0028434_r8,  0.0028596_r8,  0.0028756_r8,  0.0028913_r8, &
       0.0029068_r8,  0.0029220_r8,  0.0029370_r8,  0.0029518_r8,  0.0029664_r8, &
       0.0029807_r8,  0.0029949_r8,  0.0030088_r8,  0.0030225_r8,  0.0030361_r8, &
       0.0030494_r8,  0.0030626_r8, &
       0.0014557_r8,  0.0015446_r8,  0.0016236_r8,  0.0016950_r8,  0.0017605_r8, &
       0.0018211_r8,  0.0018777_r8,  0.0019308_r8,  0.0019809_r8,  0.0020284_r8, &
       0.0020736_r8,  0.0021167_r8,  0.0021121_r8,  0.0021973_r8,  0.0022353_r8, &
       0.0022718_r8,  0.0023069_r8,  0.0023409_r8,  0.0023737_r8,  0.0024055_r8, &
       0.0024362_r8,  0.0024661_r8,  0.0024950_r8,  0.0025232_r8,  0.0025506_r8, &
       0.0025772_r8,  0.0026031_r8,  0.0026284_r8,  0.0026531_r8,  0.0026771_r8, &
       0.0027006_r8,  0.0027235_r8,  0.0027460_r8,  0.0027679_r8,  0.0027893_r8, &
       0.0028103_r8,  0.0028309_r8,  0.0028510_r8,  0.0028707_r8,  0.0028900_r8, &
       0.0029090_r8,  0.0029276_r8,  0.0029458_r8,  0.0029637_r8,  0.0029813_r8, &
       0.0029986_r8,  0.0030156_r8,  0.0030322_r8,  0.0030486_r8,  0.0030647_r8, &
       0.0030806_r8,  0.0030962_r8,  0.0031115_r8,  0.0031266_r8,  0.0031414_r8, &
       0.0031560_r8,  0.0031704_r8,  0.0031846_r8,  0.0031986_r8,  0.0032123_r8, &
       0.0032259_r8,  0.0032392_r8, &
       0.0015532_r8,  0.0016476_r8,  0.0017317_r8,  0.0018078_r8,  0.0018775_r8, &
       0.0019422_r8,  0.0020024_r8,  0.0020590_r8,  0.0021123_r8,  0.0021169_r8, &
       0.0022106_r8,  0.0022563_r8,  0.0022999_r8,  0.0023416_r8,  0.0023817_r8, &
       0.0024202_r8,  0.0024572_r8,  0.0024929_r8,  0.0025273_r8,  0.0025607_r8, &
       0.0025929_r8,  0.0026241_r8,  0.0026543_r8,  0.0026837_r8,  0.0027122_r8, &
       0.0027399_r8,  0.0027668_r8,  0.0027931_r8,  0.0028186_r8,  0.0028435_r8, &
       0.0028678_r8,  0.0028915_r8,  0.0029146_r8,  0.0029372_r8,  0.0029592_r8, &
       0.0029808_r8,  0.0030019_r8,  0.0030225_r8,  0.0030427_r8,  0.0030625_r8, &
       0.0030819_r8,  0.0031009_r8,  0.0031195_r8,  0.0031377_r8,  0.0031556_r8, &
       0.0031732_r8,  0.0031904_r8,  0.0032073_r8,  0.0032239_r8,  0.0032402_r8, &
       0.0032563_r8,  0.0032720_r8,  0.0032875_r8,  0.0033027_r8,  0.0033177_r8, &
       0.0033324_r8,  0.0033468_r8,  0.0033611_r8,  0.0033751_r8,  0.0033889_r8, &
       0.0034025_r8,  0.0034159_r8, &
       0.0016566_r8,  0.0017571_r8,  0.0018467_r8,  0.0019278_r8,  0.0020021_r8, &
       0.0020708_r8,  0.0021349_r8,  0.0021949_r8,  0.0022514_r8,  0.0023047_r8, &
       0.0023554_r8,  0.0024035_r8,  0.0024494_r8,  0.0024932_r8,  0.0025352_r8, &
       0.0025755_r8,  0.0026142_r8,  0.0026515_r8,  0.0026874_r8,  0.0027220_r8, &
       0.0027555_r8,  0.0027878_r8,  0.0028191_r8,  0.0028495_r8,  0.0028789_r8, &
       0.0029075_r8,  0.0029352_r8,  0.0029621_r8,  0.0029883_r8,  0.0030138_r8, &
       0.0030387_r8,  0.0030629_r8,  0.0030864_r8,  0.0031094_r8,  0.0031319_r8, &
       0.0031538_r8,  0.0031752_r8,  0.0031961_r8,  0.0032166_r8,  0.0032365_r8, &
       0.0032561_r8,  0.0032752_r8,  0.0032940_r8,  0.0033123_r8,  0.0033303_r8, &
       0.0033479_r8,  0.0033652_r8,  0.0033821_r8,  0.0033988_r8,  0.0034151_r8, &
       0.0034310_r8,  0.0034467_r8,  0.0034622_r8,  0.0034773_r8,  0.0034922_r8, &
       0.0035068_r8,  0.0035211_r8,  0.0035353_r8,  0.0035491_r8,  0.0035628_r8, &
       0.0035762_r8,  0.0035894_r8, &
       0.0017664_r8,  0.0018736_r8,  0.0019689_r8,  0.0020552_r8,  0.0021342_r8, &
       0.0022071_r8,  0.0022749_r8,  0.0023383_r8,  0.0023978_r8,  0.0024539_r8, &
       0.0025070_r8,  0.0025574_r8,  0.0026054_r8,  0.0026511_r8,  0.0026948_r8, &
       0.0027366_r8,  0.0027767_r8,  0.0028153_r8,  0.0028523_r8,  0.0028880_r8, &
       0.0029224_r8,  0.0029556_r8,  0.0029877_r8,  0.0030187_r8,  0.0030487_r8, &
       0.0030778_r8,  0.0031060_r8,  0.0031334_r8,  0.0031599_r8,  0.0031857_r8, &
       0.0032108_r8,  0.0032352_r8,  0.0032590_r8,  0.0032821_r8,  0.0033047_r8, &
       0.0033267_r8,  0.0033481_r8,  0.0033690_r8,  0.0033894_r8,  0.0034094_r8, &
       0.0034289_r8,  0.0034479_r8,  0.0034665_r8,  0.0034847_r8,  0.0035025_r8, &
       0.0035200_r8,  0.0035371_r8,  0.0035538_r8,  0.0035702_r8,  0.0035862_r8, &
       0.0036020_r8,  0.0036174_r8,  0.0036325_r8,  0.0036474_r8,  0.0036620_r8, &
       0.0036763_r8,  0.0036903_r8,  0.0037041_r8,  0.0037177_r8,  0.0037310_r8, &
       0.0037441_r8,  0.0037569_r8, &
       0.0018832_r8,  0.0019973_r8,  0.0020987_r8,  0.0021902_r8,  0.0022737_r8, &
       0.0023505_r8,  0.0024220_r8,  0.0024885_r8,  0.0025508_r8,  0.0026093_r8, &
       0.0026646_r8,  0.0027169_r8,  0.0027666_r8,  0.0028138_r8,  0.0028589_r8, &
       0.0029018_r8,  0.0029430_r8,  0.0029824_r8,  0.0030202_r8,  0.0030565_r8, &
       0.0030915_r8,  0.0031252_r8,  0.0031576_r8,  0.0031890_r8,  0.0032192_r8, &
       0.0032485_r8,  0.0032768_r8,  0.0033042_r8,  0.0033308_r8,  0.0033566_r8, &
       0.0033816_r8,  0.0034059_r8,  0.0034295_r8,  0.0034525_r8,  0.0034748_r8, &
       0.0034966_r8,  0.0035177_r8,  0.0035384_r8,  0.0035585_r8,  0.0035781_r8, &
       0.0035972_r8,  0.0036159_r8,  0.0036341_r8,  0.0036520_r8,  0.0036694_r8, &
       0.0036864_r8,  0.0037031_r8,  0.0037193_r8,  0.0037353_r8,  0.0037509_r8, &
       0.0037662_r8,  0.0037811_r8,  0.0037958_r8,  0.0038102_r8,  0.0038243_r8, &
       0.0038381_r8,  0.0038517_r8,  0.0038650_r8,  0.0038780_r8,  0.0038908_r8, &
       0.0039034_r8,  0.0039158_r8, &
       0.0020072_r8,  0.0021282_r8,  0.0022356_r8,  0.0023321_r8,  0.0024200_r8, &
       0.0025006_r8,  0.0025751_r8,  0.0026443_r8,  0.0027089_r8,  0.0027695_r8, &
       0.0028265_r8,  0.0028803_r8,  0.0029311_r8,  0.0029794_r8,  0.0030252_r8, &
       0.0030689_r8,  0.0031106_r8,  0.0031504_r8,  0.0031886_r8,  0.0032251_r8, &
       0.0032602_r8,  0.0032939_r8,  0.0033263_r8,  0.0033575_r8,  0.0033876_r8, &
       0.0034166_r8,  0.0034447_r8,  0.0034718_r8,  0.0034980_r8,  0.0035234_r8, &
       0.0035480_r8,  0.0035719_r8,  0.0035950_r8,  0.0036175_r8,  0.0036393_r8, &
       0.0036605_r8,  0.0036811_r8,  0.0037012_r8,  0.0037207_r8,  0.0037398_r8, &
       0.0037583_r8,  0.0037764_r8,  0.0037940_r8,  0.0038112_r8,  0.0038280_r8, &
       0.0038444_r8,  0.0038604_r8,  0.0038761_r8,  0.0038914_r8,  0.0039063_r8, &
       0.0039210_r8,  0.0039353_r8,  0.0039494_r8,  0.0039631_r8,  0.0039766_r8, &
       0.0039898_r8,  0.0040027_r8,  0.0040154_r8,  0.0040278_r8,  0.0040400_r8, &
       0.0040520_r8,  0.0040638_r8, &
       0.0021381_r8,  0.0022661_r8,  0.0023791_r8,  0.0024803_r8,  0.0025719_r8, &
       0.0026557_r8,  0.0027328_r8,  0.0028041_r8,  0.0028705_r8,  0.0029325_r8, &
       0.0029905_r8,  0.0030451_r8,  0.0030966_r8,  0.0031453_r8,  0.0031914_r8, &
       0.0032353_r8,  0.0032769_r8,  0.0033166_r8,  0.0033545_r8,  0.0033908_r8, &
       0.0034255_r8,  0.0034588_r8,  0.0034907_r8,  0.0035214_r8,  0.0035509_r8, &
       0.0035793_r8,  0.0036067_r8,  0.0036331_r8,  0.0036586_r8,  0.0036833_r8, &
       0.0037071_r8,  0.0037302_r8,  0.0037526_r8,  0.0037743_r8,  0.0037953_r8, &
       0.0038157_r8,  0.0038356_r8,  0.0038548_r8,  0.0038736_r8,  0.0038916_r8, &
       0.0039095_r8,  0.0039268_r8,  0.0039436_r8,  0.0039600_r8,  0.0039761_r8, &
       0.0039917_r8,  0.0040069_r8,  0.0040218_r8,  0.0040364_r8,  0.0040506_r8, &
       0.0040645_r8,  0.0040781_r8,  0.0040914_r8,  0.0041044_r8,  0.0041172_r8, &
       0.0041296_r8,  0.0041419_r8,  0.0041539_r8,  0.0041656_r8,  0.0041772_r8, &
       0.0041885_r8,  0.0041996_r8, &
       0.0022756_r8,  0.0024100_r8,  0.0025280_r8,  0.0026332_r8,  0.0027280_r8, &
       0.0028142_r8,  0.0028931_r8,  0.0029658_r8,  0.0030331_r8,  0.0030957_r8, &
       0.0031541_r8,  0.0032089_r8,  0.0032603_r8,  0.0033088_r8,  0.0033545_r8, &
       0.0033978_r8,  0.0034389_r8,  0.0034780_r8,  0.0035151_r8,  0.0035506_r8, &
       0.0035844_r8,  0.0036168_r8,  0.0036478_r8,  0.0036775_r8,  0.0037061_r8, &
       0.0037335_r8,  0.0037599_r8,  0.0037853_r8,  0.0038098_r8,  0.0038335_r8, &
       0.0038563_r8,  0.0038784_r8,  0.0038998_r8,  0.0039205_r8,  0.0039405_r8, &
       0.0039599_r8,  0.0039788_r8,  0.0039971_r8,  0.0040149_r8,  0.0040322_r8, &
       0.0040490_r8,  0.0040653_r8,  0.0040813_r8,  0.0040968_r8,  0.0041119_r8, &
       0.0041267_r8,  0.0041411_r8,  0.0041552_r8,  0.0041689_r8,  0.0041823_r8, &
       0.0041954_r8,  0.0042082_r8,  0.0042208_r8,  0.0042331_r8,  0.0042451_r8, &
       0.0042569_r8,  0.0042684_r8,  0.0042797_r8,  0.0042908_r8,  0.0043017_r8, &
       0.0043123_r8,  0.0043228_r8, &
       0.0024185_r8,  0.0025586_r8,  0.0026808_r8,  0.0027890_r8,  0.0028859_r8, &
       0.0029735_r8,  0.0030533_r8,  0.0031264_r8,  0.0031938_r8,  0.0032562_r8, &
       0.0033142_r8,  0.0033683_r8,  0.0034190_r8,  0.0034665_r8,  0.0035113_r8, &
       0.0035535_r8,  0.0035934_r8,  0.0036313_r8,  0.0036672_r8,  0.0037014_r8, &
       0.0037340_r8,  0.0037651_r8,  0.0037948_r8,  0.0038232_r8,  0.0038505_r8, &
       0.0038767_r8,  0.0039018_r8,  0.0039260_r8,  0.0039493_r8,  0.0039717_r8, &
       0.0039934_r8,  0.0040143_r8,  0.0040345_r8,  0.0040540_r8,  0.0040730_r8, &
       0.0040913_r8,  0.0041091_r8,  0.0041264_r8,  0.0041432_r8,  0.0041595_r8, &
       0.0041753_r8,  0.0041907_r8,  0.0042057_r8,  0.0042204_r8,  0.0042346_r8, &
       0.0042485_r8,  0.0042621_r8,  0.0042754_r8,  0.0042883_r8,  0.0043009_r8, &
       0.0043133_r8,  0.0043254_r8,  0.0043372_r8,  0.0043488_r8,  0.0043601_r8, &
       0.0043712_r8,  0.0043821_r8,  0.0043928_r8,  0.0044032_r8,  0.0044135_r8, &
       0.0044236_r8,  0.0044335_r8, &
       0.0025653_r8,  0.0027099_r8,  0.0028350_r8,  0.0029450_r8,  0.0030428_r8, &
       0.0031307_r8,  0.0032102_r8,  0.0032828_r8,  0.0033493_r8,  0.0034106_r8, &
       0.0034673_r8,  0.0035200_r8,  0.0035692_r8,  0.0036152_r8,  0.0036583_r8, &
       0.0036990_r8,  0.0037373_r8,  0.0037735_r8,  0.0038078_r8,  0.0038404_r8, &
       0.0038714_r8,  0.0039009_r8,  0.0039291_r8,  0.0039560_r8,  0.0039818_r8, &
       0.0040065_r8,  0.0040303_r8,  0.0040531_r8,  0.0040751_r8,  0.0040962_r8, &
       0.0041166_r8,  0.0041363_r8,  0.0041553_r8,  0.0041738_r8,  0.0041916_r8, &
       0.0042089_r8,  0.0042256_r8,  0.0042419_r8,  0.0042577_r8,  0.0042730_r8, &
       0.0042879_r8,  0.0043025_r8,  0.0043166_r8,  0.0043304_r8,  0.0043439_r8, &
       0.0043570_r8,  0.0043698_r8,  0.0043823_r8,  0.0043945_r8,  0.0044064_r8, &
       0.0044181_r8,  0.0044296_r8,  0.0044407_r8,  0.0044517_r8,  0.0044624_r8, &
       0.0044730_r8,  0.0044833_r8,  0.0044934_r8,  0.0045033_r8,  0.0045130_r8, &
       0.0045226_r8,  0.0045320_r8, &
       0.0027124_r8,  0.0028597_r8,  0.0029861_r8,  0.0030963_r8,  0.0031936_r8, &
       0.0032805_r8,  0.0033586_r8,  0.0034295_r8,  0.0034941_r8,  0.0035534_r8, &
       0.0036081_r8,  0.0036587_r8,  0.0037058_r8,  0.0037497_r8,  0.0037908_r8, &
       0.0038294_r8,  0.0038657_r8,  0.0039000_r8,  0.0039324_r8,  0.0039631_r8, &
       0.0039923_r8,  0.0040201_r8,  0.0040467_r8,  0.0040720_r8,  0.0040963_r8, &
       0.0041195_r8,  0.0041418_r8,  0.0041633_r8,  0.0041839_r8,  0.0042038_r8, &
       0.0042230_r8,  0.0042415_r8,  0.0042594_r8,  0.0042768_r8,  0.0042935_r8, &
       0.0043098_r8,  0.0043256_r8,  0.0043409_r8,  0.0043558_r8,  0.0043703_r8, &
       0.0043844_r8,  0.0043982_r8,  0.0044115_r8,  0.0044246_r8,  0.0044373_r8, &
       0.0044497_r8,  0.0044619_r8,  0.0044737_r8,  0.0044853_r8,  0.0044966_r8, &
       0.0045077_r8,  0.0045186_r8,  0.0045292_r8,  0.0045397_r8,  0.0045499_r8, &
       0.0045599_r8,  0.0045697_r8,  0.0045793_r8,  0.0045888_r8,  0.0045981_r8, &
       0.0046072_r8,  0.0046162_r8, &
       0.0028481_r8,  0.0029956_r8,  0.0031209_r8,  0.0032292_r8,  0.0033241_r8, &
       0.0034083_r8,  0.0034836_r8,  0.0035515_r8,  0.0036132_r8,  0.0036696_r8, &
       0.0037214_r8,  0.0037692_r8,  0.0038136_r8,  0.0038548_r8,  0.0038934_r8, &
       0.0039296_r8,  0.0039636_r8,  0.0039956_r8,  0.0040260_r8,  0.0040547_r8, &
       0.0040820_r8,  0.0041080_r8,  0.0041328_r8,  0.0041565_r8,  0.0041792_r8, &
       0.0042010_r8,  0.0042219_r8,  0.0042419_r8,  0.0042613_r8,  0.0042800_r8, &
       0.0042980_r8,  0.0043154_r8,  0.0043322_r8,  0.0043485_r8,  0.0043643_r8, &
       0.0043796_r8,  0.0043945_r8,  0.0044090_r8,  0.0044231_r8,  0.0044367_r8, &
       0.0044501_r8,  0.0044631_r8,  0.0044757_r8,  0.0044881_r8,  0.0045001_r8, &
       0.0045119_r8,  0.0045234_r8,  0.0045347_r8,  0.0045457_r8,  0.0045564_r8, &
       0.0045670_r8,  0.0045773_r8,  0.0045874_r8,  0.0045973_r8,  0.0046070_r8, &
       0.0046166_r8,  0.0046259_r8,  0.0046351_r8,  0.0046441_r8,  0.0046530_r8, &
       0.0046616_r8,  0.0046702_r8, &
       0.0029341_r8,  0.0030768_r8,  0.0031968_r8,  0.0032997_r8,  0.0033892_r8, &
       0.0034681_r8,  0.0035383_r8,  0.0036014_r8,  0.0036584_r8,  0.0037104_r8, &
       0.0037581_r8,  0.0038020_r8,  0.0038427_r8,  0.0038805_r8,  0.0039159_r8, &
       0.0039490_r8,  0.0039801_r8,  0.0040095_r8,  0.0040373_r8,  0.0040637_r8, &
       0.0040888_r8,  0.0041127_r8,  0.0041355_r8,  0.0041573_r8,  0.0041782_r8, &
       0.0041983_r8,  0.0042175_r8,  0.0042361_r8,  0.0042540_r8,  0.0042713_r8, &
       0.0042880_r8,  0.0043041_r8,  0.0043197_r8,  0.0043349_r8,  0.0043495_r8, &
       0.0043638_r8,  0.0043777_r8,  0.0043911_r8,  0.0044042_r8,  0.0044170_r8, &
       0.0044294_r8,  0.0044416_r8,  0.0044534_r8,  0.0044649_r8,  0.0044762_r8, &
       0.0044872_r8,  0.0044980_r8,  0.0045085_r8,  0.0045188_r8,  0.0045289_r8, &
       0.0045387_r8,  0.0045484_r8,  0.0045579_r8,  0.0045672_r8,  0.0045763_r8, &
       0.0045852_r8,  0.0045940_r8,  0.0046026_r8,  0.0046110_r8,  0.0046193_r8, &
       0.0046275_r8,  0.0046355_r8, &
       0.0029122_r8,  0.0030427_r8,  0.0031513_r8,  0.0032438_r8,  0.0033237_r8, &
       0.0033938_r8,  0.0034559_r8,  0.0035116_r8,  0.0035619_r8,  0.0036076_r8, &
       0.0036495_r8,  0.0036882_r8,  0.0037238_r8,  0.0037573_r8,  0.0037884_r8, &
       0.0038176_r8,  0.0038451_r8,  0.0038711_r8,  0.0038957_r8,  0.0039191_r8, &
       0.0039413_r8,  0.0039626_r8,  0.0039829_r8,  0.0040023_r8,  0.0040209_r8, &
       0.0040389_r8,  0.0040561_r8,  0.0040727_r8,  0.0040887_r8,  0.0041042_r8, &
       0.0041192_r8,  0.0041337_r8,  0.0041477_r8,  0.0041613_r8,  0.0041745_r8, &
       0.0041874_r8,  0.0041998_r8,  0.0042120_r8,  0.0042238_r8,  0.0042353_r8, &
       0.0042465_r8,  0.0042574_r8,  0.0042681_r8,  0.0042785_r8,  0.0042887_r8, &
       0.0042986_r8,  0.0043084_r8,  0.0043179_r8,  0.0043272_r8,  0.0043363_r8, &
       0.0043452_r8,  0.0043539_r8,  0.0043625_r8,  0.0043709_r8,  0.0043791_r8, &
       0.0043872_r8,  0.0043951_r8,  0.0044029_r8,  0.0044105_r8,  0.0044180_r8, &
       0.0044254_r8,  0.0044326_r8, &
       0.0027405_r8,  0.0028512_r8,  0.0029426_r8,  0.0030199_r8,  0.0030864_r8, &
       0.0031447_r8,  0.0031962_r8,  0.0032424_r8,  0.0032841_r8,  0.0033221_r8, &
       0.0033569_r8,  0.0033891_r8,  0.0034190_r8,  0.0034468_r8,  0.0034729_r8, &
       0.0034974_r8,  0.0035206_r8,  0.0035424_r8,  0.0035632_r8,  0.0035830_r8, &
       0.0036018_r8,  0.0036198_r8,  0.0036370_r8,  0.0036535_r8,  0.0036694_r8, &
       0.0036847_r8,  0.0036993_r8,  0.0037135_r8,  0.0037272_r8,  0.0037404_r8, &
       0.0037532_r8,  0.0037656_r8,  0.0037776_r8,  0.0037892_r8,  0.0038005_r8, &
       0.0038115_r8,  0.0038222_r8,  0.0038326_r8,  0.0038427_r8,  0.0038526_r8, &
       0.0038622_r8,  0.0038716_r8,  0.0038808_r8,  0.0038897_r8,  0.0038984_r8, &
       0.0039070_r8,  0.0039153_r8,  0.0039235_r8,  0.0039314_r8,  0.0039393_r8, &
       0.0039469_r8,  0.0039544_r8,  0.0039618_r8,  0.0039690_r8,  0.0039761_r8, &
       0.0039830_r8,  0.0039898_r8,  0.0039965_r8,  0.0040030_r8,  0.0040095_r8, &
       0.0040158_r8,  0.0040220_r8, &
       0.0024633_r8,  0.0025514_r8,  0.0026239_r8,  0.0026851_r8,  0.0027377_r8, &
       0.0027838_r8,  0.0028247_r8,  0.0028613_r8,  0.0028946_r8,  0.0029249_r8, &
       0.0029529_r8,  0.0029787_r8,  0.0030028_r8,  0.0030253_r8,  0.0030464_r8, &
       0.0030663_r8,  0.0030851_r8,  0.0031029_r8,  0.0031199_r8,  0.0031360_r8, &
       0.0031514_r8,  0.0031661_r8,  0.0031803_r8,  0.0031938_r8,  0.0032068_r8, &
       0.0032194_r8,  0.0032314_r8,  0.0032431_r8,  0.0032544_r8,  0.0032652_r8, &
       0.0032758_r8,  0.0032860_r8,  0.0032959_r8,  0.0033055_r8,  0.0033148_r8, &
       0.0033239_r8,  0.0033327_r8,  0.0033413_r8,  0.0033497_r8,  0.0033578_r8, &
       0.0033658_r8,  0.0033735_r8,  0.0033811_r8,  0.0033885_r8,  0.0033957_r8, &
       0.0034028_r8,  0.0034097_r8,  0.0034164_r8,  0.0034230_r8,  0.0034295_r8, &
       0.0034359_r8,  0.0034421_r8,  0.0034482_r8,  0.0034541_r8,  0.0034600_r8, &
       0.0034657_r8,  0.0034714_r8,  0.0034769_r8,  0.0034823_r8,  0.0034877_r8, &
       0.0034929_r8,  0.0034981_r8, &
       0.0021142_r8,  0.0022278_r8,  0.0022837_r8,  0.0023309_r8,  0.0023717_r8, &
       0.0024075_r8,  0.0024394_r8,  0.0024681_r8,  0.0024943_r8,  0.0025182_r8, &
       0.0025404_r8,  0.0025609_r8,  0.0025801_r8,  0.0025980_r8,  0.0026149_r8, &
       0.0026308_r8,  0.0026459_r8,  0.0026602_r8,  0.0026738_r8,  0.0026868_r8, &
       0.0026992_r8,  0.0027111_r8,  0.0027225_r8,  0.0027334_r8,  0.0027439_r8, &
       0.0027541_r8,  0.0027638_r8,  0.0027733_r8,  0.0027824_r8,  0.0027912_r8, &
       0.0027997_r8,  0.0028080_r8,  0.0028160_r8,  0.0028238_r8,  0.0028314_r8, &
       0.0028387_r8,  0.0028459_r8,  0.0028529_r8,  0.0028597_r8,  0.0028663_r8, &
       0.0028727_r8,  0.0028791_r8,  0.0028852_r8,  0.0028912_r8,  0.0028971_r8, &
       0.0029028_r8,  0.0029084_r8,  0.0029139_r8,  0.0029193_r8,  0.0029246_r8, &
       0.0029297_r8,  0.0029348_r8,  0.0029398_r8,  0.0029446_r8,  0.0029494_r8, &
       0.0029541_r8,  0.0029587_r8,  0.0029632_r8,  0.0029676_r8,  0.0029720_r8, &
       0.0029761_r8,  0.0029805_r8, &
       0.0018726_r8,  0.0019238_r8,  0.0019660_r8,  0.0020019_r8,  0.0020331_r8, &
       0.0020606_r8,  0.0020852_r8,  0.0021074_r8,  0.0021278_r8,  0.0021464_r8, &
       0.0021179_r8,  0.0021798_r8,  0.0021948_r8,  0.0022089_r8,  0.0022221_r8, &
       0.0022347_r8,  0.0022466_r8,  0.0022578_r8,  0.0022686_r8,  0.0022789_r8, &
       0.0022887_r8,  0.0022981_r8,  0.0023071_r8,  0.0023158_r8,  0.0023241_r8, &
       0.0023321_r8,  0.0023399_r8,  0.0023474_r8,  0.0023546_r8,  0.0023616_r8, &
       0.0023684_r8,  0.0023750_r8,  0.0023814_r8,  0.0023876_r8,  0.0023936_r8, &
       0.0023995_r8,  0.0024052_r8,  0.0024108_r8,  0.0024162_r8,  0.0024215_r8, &
       0.0024266_r8,  0.0024317_r8,  0.0024366_r8,  0.0024414_r8,  0.0024461_r8, &
       0.0024507_r8,  0.0024552_r8,  0.0024596_r8,  0.0024639_r8,  0.0024681_r8, &
       0.0024722_r8,  0.0024763_r8,  0.0024802_r8,  0.0024841_r8,  0.0024880_r8, &
       0.0024917_r8,  0.0024954_r8,  0.0024990_r8,  0.0025026_r8,  0.0025060_r8, &
       0.0025095_r8,  0.0025129_r8, &
       0.0016337_r8,  0.0016718_r8,  0.0017033_r8,  0.0017303_r8,  0.0017537_r8, &
       0.0017745_r8,  0.0017931_r8,  0.0018100_r8,  0.0018254_r8,  0.0018397_r8, &
       0.0018529_r8,  0.0018651_r8,  0.0018766_r8,  0.0018874_r8,  0.0018976_r8, &
       0.0019073_r8,  0.0019164_r8,  0.0019251_r8,  0.0019334_r8,  0.0019413_r8, &
       0.0019489_r8,  0.0019562_r8,  0.0019631_r8,  0.0019698_r8,  0.0019763_r8, &
       0.0019825_r8,  0.0019885_r8,  0.0019944_r8,  0.0020000_r8,  0.0020054_r8, &
       0.0020107_r8,  0.0020158_r8,  0.0020208_r8,  0.0020256_r8,  0.0020303_r8, &
       0.0020349_r8,  0.0020394_r8,  0.0020437_r8,  0.0020479_r8,  0.0020521_r8, &
       0.0020561_r8,  0.0020600_r8,  0.0020639_r8,  0.0020676_r8,  0.0020713_r8, &
       0.0020749_r8,  0.0020784_r8,  0.0020819_r8,  0.0020852_r8,  0.0020886_r8, &
       0.0020918_r8,  0.0020950_r8,  0.0020981_r8,  0.0021011_r8,  0.0021041_r8, &
       0.0021071_r8,  0.0021100_r8,  0.0021128_r8,  0.0021156_r8,  0.0021183_r8, &
       0.0021210_r8,  0.0021236_r8, &
       0.0014740_r8,  0.0015024_r8,  0.0015259_r8,  0.0015460_r8,  0.0015636_r8, &
       0.0015791_r8,  0.0015931_r8,  0.0016058_r8,  0.0016174_r8,  0.0016282_r8, &
       0.0016381_r8,  0.0016474_r8,  0.0016561_r8,  0.0016643_r8,  0.0016720_r8, &
       0.0016793_r8,  0.0016863_r8,  0.0016929_r8,  0.0016992_r8,  0.0017052_r8, &
       0.0017110_r8,  0.0017165_r8,  0.0017219_r8,  0.0017270_r8,  0.0017319_r8, &
       0.0017367_r8,  0.0017413_r8,  0.0017458_r8,  0.0017501_r8,  0.0017543_r8, &
       0.0017584_r8,  0.0017623_r8,  0.0017661_r8,  0.0017699_r8,  0.0017735_r8, &
       0.0017770_r8,  0.0017804_r8,  0.0017838_r8,  0.0017870_r8,  0.0017902_r8, &
       0.0017933_r8,  0.0017964_r8,  0.0017993_r8,  0.0018023_r8,  0.0018051_r8, &
       0.0018079_r8,  0.0018106_r8,  0.0018132_r8,  0.0018158_r8,  0.0018184_r8, &
       0.0018209_r8,  0.0018233_r8,  0.0018258_r8,  0.0018281_r8,  0.0018304_r8, &
       0.0018327_r8,  0.0018349_r8,  0.0018371_r8,  0.0018393_r8,  0.0018414_r8, &
       0.0018434_r8,  0.0018455_r8, &
       0.0013895_r8,  0.0014110_r8,  0.0014289_r8,  0.0014441_r8,  0.0014574_r8, &
       0.0014692_r8,  0.0014798_r8,  0.0014894_r8,  0.0014982_r8,  0.0015064_r8, &
       0.0015139_r8,  0.0015210_r8,  0.0015277_r8,  0.0015338_r8,  0.0015398_r8, &
       0.0015454_r8,  0.0015508_r8,  0.0015558_r8,  0.0015607_r8,  0.0015653_r8, &
       0.0015698_r8,  0.0015740_r8,  0.0015782_r8,  0.0015821_r8,  0.0015859_r8, &
       0.0015896_r8,  0.0015932_r8,  0.0015966_r8,  0.0016000_r8,  0.0016032_r8, &
       0.0016064_r8,  0.0016094_r8,  0.0016124_r8,  0.0016153_r8,  0.0016181_r8, &
       0.0016208_r8,  0.0016235_r8,  0.0016261_r8,  0.0016286_r8,  0.0016311_r8, &
       0.0016335_r8,  0.0016358_r8,  0.0016381_r8,  0.0016404_r8,  0.0016426_r8, &
       0.0016447_r8,  0.0016468_r8,  0.0016489_r8,  0.0016509_r8,  0.0016529_r8, &
       0.0016548_r8,  0.0016567_r8,  0.0016586_r8,  0.0016604_r8,  0.0016622_r8, &
       0.0016639_r8,  0.0016657_r8,  0.0016673_r8,  0.0016690_r8,  0.0016706_r8, &
       0.0016722_r8,  0.0016738_r8, &
       0.0013502_r8,  0.0013669_r8,  0.0013807_r8,  0.0013924_r8,  0.0014027_r8, &
       0.0014118_r8,  0.0014200_r8,  0.0014274_r8,  0.0014343_r8,  0.0014406_r8, &
       0.0014465_r8,  0.0014520_r8,  0.0014571_r8,  0.0014620_r8,  0.0014666_r8, &
       0.0014710_r8,  0.0014751_r8,  0.0014791_r8,  0.0014829_r8,  0.0014865_r8, &
       0.0014900_r8,  0.0014933_r8,  0.0014966_r8,  0.0014997_r8,  0.0015027_r8, &
       0.0015055_r8,  0.0015083_r8,  0.0015109_r8,  0.0015136_r8,  0.0015162_r8, &
       0.0015186_r8,  0.0015210_r8,  0.0015234_r8,  0.0015256_r8,  0.0015278_r8, &
       0.0015299_r8,  0.0015320_r8,  0.0015340_r8,  0.0015360_r8,  0.0015380_r8, &
       0.0015398_r8,  0.0015417_r8,  0.0015435_r8,  0.0015452_r8,  0.0015469_r8, &
       0.0015486_r8,  0.0015503_r8,  0.0015519_r8,  0.0015534_r8,  0.0015550_r8, &
       0.0015565_r8,  0.0015580_r8,  0.0015594_r8,  0.0015608_r8,  0.0015622_r8, &
       0.0015636_r8,  0.0015649_r8,  0.0015663_r8,  0.0015676_r8,  0.0015688_r8, &
       0.0015701_r8,  0.0015713_r8, &
       0.0013341_r8,  0.0013476_r8,  0.0013588_r8,  0.0013683_r8,  0.0013766_r8, &
       0.0013840_r8,  0.0013907_r8,  0.0013967_r8,  0.0014023_r8,  0.0014074_r8, &
       0.0014122_r8,  0.0014167_r8,  0.0014209_r8,  0.0014248_r8,  0.0014286_r8, &
       0.0014321_r8,  0.0014355_r8,  0.0014387_r8,  0.0014418_r8,  0.0014447_r8, &
       0.0014476_r8,  0.0014503_r8,  0.0014529_r8,  0.0014554_r8,  0.0014578_r8, &
       0.0014602_r8,  0.0014624_r8,  0.0014646_r8,  0.0014667_r8,  0.0014688_r8, &
       0.0014708_r8,  0.0014727_r8,  0.0014746_r8,  0.0014764_r8,  0.0014782_r8, &
       0.0014799_r8,  0.0014816_r8,  0.0014833_r8,  0.0014849_r8,  0.0014864_r8, &
       0.0014879_r8,  0.0014894_r8,  0.0014909_r8,  0.0014923_r8,  0.0014937_r8, &
       0.0014950_r8,  0.0014964_r8,  0.0014977_r8,  0.0014989_r8,  0.0015002_r8, &
       0.0015014_r8,  0.0015026_r8,  0.0015038_r8,  0.0015049_r8,  0.0015060_r8, &
       0.0015071_r8,  0.0015082_r8,  0.0015093_r8,  0.0015103_r8,  0.0015114_r8, &
       0.0015124_r8,  0.0015134_r8, &
       0.0013255_r8,  0.0013373_r8,  0.0013470_r8,  0.0013554_r8,  0.0013626_r8, &
       0.0013691_r8,  0.0013749_r8,  0.0013803_r8,  0.0013851_r8,  0.0013896_r8, &
       0.0013938_r8,  0.0013977_r8,  0.0014014_r8,  0.0014049_r8,  0.0014082_r8, &
       0.0014113_r8,  0.0014142_r8,  0.0014171_r8,  0.0014197_r8,  0.0014223_r8, &
       0.0014248_r8,  0.0014272_r8,  0.0014294_r8,  0.0014316_r8,  0.0014337_r8, &
       0.0014358_r8,  0.0014378_r8,  0.0014397_r8,  0.0014415_r8,  0.0014433_r8, &
       0.0014450_r8,  0.0014467_r8,  0.0014483_r8,  0.0014499_r8,  0.0014515_r8, &
       0.0014530_r8,  0.0014544_r8,  0.0014559_r8,  0.0014573_r8,  0.0014586_r8, &
       0.0014599_r8,  0.0014612_r8,  0.0014625_r8,  0.0014637_r8,  0.0014649_r8, &
       0.0014661_r8,  0.0014672_r8,  0.0014684_r8,  0.0014695_r8,  0.0014705_r8, &
       0.0014716_r8,  0.0014726_r8,  0.0014736_r8,  0.0014746_r8,  0.0014756_r8, &
       0.0014766_r8,  0.0014775_r8,  0.0014784_r8,  0.0014793_r8,  0.0014802_r8, &
       0.0014811_r8,  0.0014820_r8, &
       0.0013126_r8,  0.0013234_r8,  0.0013324_r8,  0.0013401_r8,  0.0013469_r8, &
       0.0013529_r8,  0.0013583_r8,  0.0013632_r8,  0.0013677_r8,  0.0013719_r8, &
       0.0013758_r8,  0.0013795_r8,  0.0013829_r8,  0.0013861_r8,  0.0013891_r8, &
       0.0013920_r8,  0.0013947_r8,  0.0013974_r8,  0.0013998_r8,  0.0014022_r8, &
       0.0014045_r8,  0.0014067_r8,  0.0014088_r8,  0.0014108_r8,  0.0014127_r8, &
       0.0014146_r8,  0.0014164_r8,  0.0014181_r8,  0.0014198_r8,  0.0014215_r8, &
       0.0014230_r8,  0.0014246_r8,  0.0014261_r8,  0.0014275_r8,  0.0014289_r8, &
       0.0014303_r8,  0.0014316_r8,  0.0014329_r8,  0.0014341_r8,  0.0014354_r8, &
       0.0014366_r8,  0.0014377_r8,  0.0014389_r8,  0.0014400_r8,  0.0014411_r8, &
       0.0014421_r8,  0.0014432_r8,  0.0014442_r8,  0.0014452_r8,  0.0014462_r8, &
       0.0014471_r8,  0.0014480_r8,  0.0014490_r8,  0.0014499_r8,  0.0014507_r8, &
       0.0014516_r8,  0.0014525_r8,  0.0014533_r8,  0.0014541_r8,  0.0014549_r8, &
       0.0014557_r8,  0.0014565_r8, &
       0.0012882_r8,  0.0012983_r8,  0.0013066_r8,  0.0013138_r8,  0.0013202_r8,  &
       0.0013258_r8,  0.0013309_r8,  0.0013355_r8,  0.0013398_r8,  0.0013437_r8,  &
       0.0013473_r8,  0.0013507_r8,  0.0013539_r8,  0.0013569_r8,  0.0013598_r8,  &
       0.0013625_r8,  0.0013650_r8,  0.0013674_r8,  0.0013697_r8,  0.0013719_r8,  &
       0.0013740_r8,  0.0013760_r8,  0.0013780_r8,  0.0013798_r8,  0.0013816_r8,  &
       0.0013833_r8,  0.0013850_r8,  0.0013865_r8,  0.0013881_r8,  0.0013896_r8,  &
       0.0013910_r8,  0.0013924_r8,  0.0013938_r8,  0.0013951_r8,  0.0013963_r8,  &
       0.0013976_r8,  0.0013988_r8,  0.0013999_r8,  0.0014011_r8,  0.0014022_r8,  &
       0.0014033_r8,  0.0014043_r8,  0.0014054_r8,  0.0014064_r8,  0.0014073_r8,  &
       0.0014083_r8,  0.0014092_r8,  0.0014101_r8,  0.0014110_r8,  0.0014119_r8,  &
       0.0014128_r8,  0.0014136_r8,  0.0014145_r8,  0.0014153_r8,  0.0014161_r8,  &
       0.0014168_r8,  0.0014176_r8,  0.0014183_r8,  0.0014191_r8,  0.0014198_r8,  &
       0.0014205_r8,  0.0014212_r8/),SHAPE=(/nx,ny/))







       cah=RESHAPE(SOURCE=(/&
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8, &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8, &
       0.0000559_r8,  0.0000772_r8,  0.0001059_r8,  0.0001439_r8,  0.0001936_r8, &
       0.0002575_r8,  0.0003384_r8,  0.0004400_r8,  0.0005662_r8,  0.0007219_r8, &
       0.0009131_r8,  0.0011470_r8,  0.0014327_r8,  0.0017806_r8,  0.0022021_r8, &
       0.0027093_r8,  0.0033141_r8,  0.0040280_r8,  0.0048609_r8,  0.0058217_r8, &
       0.0069177_r8,  0.0081559_r8,  0.0095430_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8, &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8, &
       0.0000559_r8,  0.0000772_r8,  0.0001059_r8,  0.0001439_r8,  0.0001936_r8, &
       0.0002575_r8,  0.0003384_r8,  0.0004400_r8,  0.0005662_r8,  0.0007219_r8, &
       0.0009130_r8,  0.0011470_r8,  0.0014326_r8,  0.0017805_r8,  0.0022020_r8, &
       0.0027091_r8,  0.0033139_r8,  0.0040276_r8,  0.0048605_r8,  0.0058211_r8, &
       0.0069170_r8,  0.0081551_r8,  0.0095420_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8, &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8, &
       0.0000559_r8,  0.0000772_r8,  0.0001059_r8,  0.0001439_r8,  0.0001936_r8, &
       0.0002574_r8,  0.0003384_r8,  0.0004399_r8,  0.0005661_r8,  0.0007218_r8, &
       0.0009129_r8,  0.0011468_r8,  0.0014325_r8,  0.0017803_r8,  0.0022017_r8, &
       0.0027088_r8,  0.0033135_r8,  0.0040271_r8,  0.0048599_r8,  0.0058204_r8, &
       0.0069161_r8,  0.0081539_r8,  0.0095406_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8, &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8, &
       0.0000559_r8,  0.0000772_r8,  0.0001059_r8,  0.0001439_r8,  0.0001936_r8, &
       0.0002574_r8,  0.0003384_r8,  0.0004399_r8,  0.0005661_r8,  0.0007217_r8, &
       0.0009128_r8,  0.0011467_r8,  0.0014323_r8,  0.0017800_r8,  0.0022014_r8, &
       0.0027084_r8,  0.0033130_r8,  0.0040265_r8,  0.0048591_r8,  0.0058194_r8, &
       0.0069148_r8,  0.0081524_r8,  0.0095387_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8, &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8, &
       0.0000559_r8,  0.0000772_r8,  0.0001059_r8,  0.0001439_r8,  0.0001935_r8, &
       0.0002574_r8,  0.0003383_r8,  0.0004398_r8,  0.0005660_r8,  0.0007216_r8, &
       0.0009127_r8,  0.0011465_r8,  0.0014320_r8,  0.0017797_r8,  0.0022010_r8, &
       0.0027078_r8,  0.0033123_r8,  0.0040256_r8,  0.0048580_r8,  0.0058180_r8, &
       0.0069132_r8,  0.0081503_r8,  0.0095361_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8, &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8, &
       0.0000559_r8,  0.0000772_r8,  0.0001059_r8,  0.0001439_r8,  0.0001935_r8, &
       0.0002573_r8,  0.0003383_r8,  0.0004398_r8,  0.0005659_r8,  0.0007215_r8, &
       0.0009125_r8,  0.0011462_r8,  0.0014317_r8,  0.0017792_r8,  0.0022004_r8, &
       0.0027071_r8,  0.0033113_r8,  0.0040244_r8,  0.0048565_r8,  0.0058162_r8, &
       0.0069109_r8,  0.0081476_r8,  0.0095328_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000053_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000147_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8,  &
       0.0000559_r8,  0.0000772_r8,  0.0001058_r8,  0.0001438_r8,  0.0001935_r8,  &
       0.0002573_r8,  0.0003382_r8,  0.0004396_r8,  0.0005657_r8,  0.0007213_r8,  &
       0.0009122_r8,  0.0011459_r8,  0.0014312_r8,  0.0017786_r8,  0.0021996_r8,  &
       0.0027061_r8,  0.0033100_r8,  0.0040228_r8,  0.0048545_r8,  0.0058137_r8,  &
       0.0069079_r8,  0.0081439_r8,  0.0095283_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8,  &
       0.0000558_r8,  0.0000772_r8,  0.0001058_r8,  0.0001438_r8,  0.0001934_r8,  &
       0.0002572_r8,  0.0003381_r8,  0.0004395_r8,  0.0005655_r8,  0.0007210_r8,  &
       0.0009119_r8,  0.0011454_r8,  0.0014306_r8,  0.0017778_r8,  0.0021985_r8,  &
       0.0027047_r8,  0.0033084_r8,  0.0040207_r8,  0.0048519_r8,  0.0058105_r8,  &
       0.0069040_r8,  0.0081391_r8,  0.0095225_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000206_r8,  0.0000288_r8,  0.0000402_r8,  &
       0.0000558_r8,  0.0000771_r8,  0.0001058_r8,  0.0001437_r8,  0.0001933_r8,  &
       0.0002571_r8,  0.0003379_r8,  0.0004393_r8,  0.0005652_r8,  0.0007206_r8,  &
       0.0009114_r8,  0.0011447_r8,  0.0014297_r8,  0.0017767_r8,  0.0021971_r8,  &
       0.0027030_r8,  0.0033061_r8,  0.0040180_r8,  0.0048485_r8,  0.0058064_r8,  &
       0.0068989_r8,  0.0081329_r8,  0.0095149_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000205_r8,  0.0000288_r8,  0.0000402_r8,  &
       0.0000558_r8,  0.0000771_r8,  0.0001057_r8,  0.0001437_r8,  0.0001932_r8,  &
       0.0002569_r8,  0.0003377_r8,  0.0004390_r8,  0.0005649_r8,  0.0007201_r8,  &
       0.0009107_r8,  0.0011439_r8,  0.0014286_r8,  0.0017753_r8,  0.0021953_r8,  &
       0.0027006_r8,  0.0033032_r8,  0.0040144_r8,  0.0048441_r8,  0.0058009_r8,  &
       0.0068922_r8,  0.0081248_r8,  0.0095051_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000205_r8,  0.0000287_r8,  0.0000401_r8,  &
       0.0000558_r8,  0.0000770_r8,  0.0001056_r8,  0.0001436_r8,  0.0001931_r8,  &
       0.0002567_r8,  0.0003375_r8,  0.0004387_r8,  0.0005644_r8,  0.0007195_r8,  &
       0.0009098_r8,  0.0011428_r8,  0.0014271_r8,  0.0017734_r8,  0.0021929_r8,  &
       0.0026976_r8,  0.0032995_r8,  0.0040097_r8,  0.0048384_r8,  0.0057939_r8,  &
       0.0068837_r8,  0.0081145_r8,  0.0094926_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000205_r8,  0.0000287_r8,  0.0000401_r8,  &
       0.0000557_r8,  0.0000770_r8,  0.0001055_r8,  0.0001434_r8,  0.0001929_r8,  &
       0.0002565_r8,  0.0003371_r8,  0.0004382_r8,  0.0005637_r8,  0.0007186_r8,  &
       0.0009087_r8,  0.0011413_r8,  0.0014252_r8,  0.0017709_r8,  0.0021898_r8,  &
       0.0026937_r8,  0.0032946_r8,  0.0040038_r8,  0.0048311_r8,  0.0057850_r8,  &
       0.0068729_r8,  0.0081013_r8,  0.0094768_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000205_r8,  0.0000287_r8,  0.0000400_r8,  &
       0.0000556_r8,  0.0000769_r8,  0.0001054_r8,  0.0001432_r8,  0.0001926_r8,  &
       0.0002561_r8,  0.0003366_r8,  0.0004376_r8,  0.0005629_r8,  0.0007175_r8,  &
       0.0009073_r8,  0.0011394_r8,  0.0014228_r8,  0.0017678_r8,  0.0021859_r8,  &
       0.0026888_r8,  0.0032885_r8,  0.0039963_r8,  0.0048218_r8,  0.0057738_r8,  &
       0.0068592_r8,  0.0080849_r8,  0.0094570_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000074_r8,  &
       0.0000104_r8,  0.0000146_r8,  0.0000205_r8,  0.0000286_r8,  0.0000400_r8,  &
       0.0000556_r8,  0.0000767_r8,  0.0001052_r8,  0.0001430_r8,  0.0001923_r8,  &
       0.0002557_r8,  0.0003361_r8,  0.0004368_r8,  0.0005619_r8,  0.0007161_r8,  &
       0.0009054_r8,  0.0011370_r8,  0.0014197_r8,  0.0017639_r8,  0.0021809_r8,  &
       0.0026826_r8,  0.0032809_r8,  0.0039869_r8,  0.0048103_r8,  0.0057597_r8,  &
       0.0068422_r8,  0.0080643_r8,  0.0094323_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000073_r8,  &
       0.0000103_r8,  0.0000145_r8,  0.0000204_r8,  0.0000286_r8,  0.0000399_r8,  &
       0.0000554_r8,  0.0000766_r8,  0.0001050_r8,  0.0001427_r8,  0.0001919_r8,  &
       0.0002552_r8,  0.0003353_r8,  0.0004358_r8,  0.0005605_r8,  0.0007144_r8,  &
       0.0009032_r8,  0.0011340_r8,  0.0014159_r8,  0.0017590_r8,  0.0021748_r8,  &
       0.0026750_r8,  0.0032715_r8,  0.0039752_r8,  0.0047961_r8,  0.0057424_r8,  &
       0.0068212_r8,  0.0080389_r8,  0.0094019_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000019_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000073_r8,  &
       0.0000103_r8,  0.0000145_r8,  0.0000204_r8,  0.0000285_r8,  0.0000398_r8,  &
       0.0000553_r8,  0.0000764_r8,  0.0001047_r8,  0.0001423_r8,  0.0001914_r8,  &
       0.0002545_r8,  0.0003344_r8,  0.0004345_r8,  0.0005589_r8,  0.0007122_r8,  &
       0.0009003_r8,  0.0011304_r8,  0.0014112_r8,  0.0017531_r8,  0.0021673_r8,  &
       0.0026656_r8,  0.0032598_r8,  0.0039609_r8,  0.0047786_r8,  0.0057211_r8,  &
       0.0067954_r8,  0.0080078_r8,  0.0093646_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000018_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000073_r8,  &
       0.0000103_r8,  0.0000145_r8,  0.0000203_r8,  0.0000284_r8,  0.0000397_r8,  &
       0.0000551_r8,  0.0000761_r8,  0.0001044_r8,  0.0001419_r8,  0.0001908_r8,  &
       0.0002536_r8,  0.0003332_r8,  0.0004330_r8,  0.0005568_r8,  0.0007095_r8,  &
       0.0008968_r8,  0.0011259_r8,  0.0014054_r8,  0.0017458_r8,  0.0021123_r8,  &
       0.0026542_r8,  0.0032457_r8,  0.0039435_r8,  0.0047573_r8,  0.0056951_r8,  &
       0.0067640_r8,  0.0079700_r8,  0.0093194_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000018_r8,  0.0000026_r8,  0.0000037_r8,  0.0000052_r8,  0.0000073_r8,  &
       0.0000102_r8,  0.0000144_r8,  0.0000202_r8,  0.0000283_r8,  0.0000395_r8,  &
       0.0000549_r8,  0.0000758_r8,  0.0001040_r8,  0.0001413_r8,  0.0001900_r8,  &
       0.0002525_r8,  0.0003318_r8,  0.0004311_r8,  0.0005543_r8,  0.0007063_r8,  &
       0.0008926_r8,  0.0011204_r8,  0.0013985_r8,  0.0017370_r8,  0.0021470_r8,  &
       0.0026404_r8,  0.0032285_r8,  0.0039224_r8,  0.0047315_r8,  0.0056637_r8,  &
       0.0067260_r8,  0.0079245_r8,  0.0092651_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000006_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000018_r8,  0.0000026_r8,  0.0000036_r8,  0.0000051_r8,  0.0000072_r8,  &
       0.0000102_r8,  0.0000143_r8,  0.0000201_r8,  0.0000282_r8,  0.0000393_r8,  &
       0.0000546_r8,  0.0000754_r8,  0.0001034_r8,  0.0001406_r8,  0.0001890_r8,  &
       0.0002512_r8,  0.0003300_r8,  0.0004287_r8,  0.0005513_r8,  0.0007023_r8,  &
       0.0008875_r8,  0.0011139_r8,  0.0013901_r8,  0.0017264_r8,  0.0021337_r8,  &
       0.0026238_r8,  0.0032080_r8,  0.0038971_r8,  0.0047005_r8,  0.0056261_r8,  &
       0.0066806_r8,  0.0078701_r8,  0.0092003_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000005_r8,  0.0000006_r8,  0.0000009_r8,  0.0000013_r8,  &
       0.0000018_r8,  0.0000026_r8,  0.0000036_r8,  0.0000051_r8,  0.0000072_r8,  &
       0.0000101_r8,  0.0000142_r8,  0.0000200_r8,  0.0000280_r8,  0.0000391_r8,  &
       0.0000543_r8,  0.0000750_r8,  0.0001028_r8,  0.0001397_r8,  0.0001878_r8,  &
       0.0002496_r8,  0.0003279_r8,  0.0004259_r8,  0.0005476_r8,  0.0006975_r8,  &
       0.0008813_r8,  0.0011060_r8,  0.0013802_r8,  0.0017138_r8,  0.0021179_r8,  &
       0.0026040_r8,  0.0031835_r8,  0.0038670_r8,  0.0046637_r8,  0.0055814_r8,  &
       0.0066267_r8,  0.0078055_r8,  0.0091235_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000005_r8,  0.0000006_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000018_r8,  0.0000025_r8,  0.0000036_r8,  0.0000051_r8,  0.0000071_r8, &
       0.0000100_r8,  0.0000141_r8,  0.0000198_r8,  0.0000278_r8,  0.0000388_r8, &
       0.0000539_r8,  0.0000744_r8,  0.0001020_r8,  0.0001386_r8,  0.0001863_r8, &
       0.0002477_r8,  0.0003253_r8,  0.0004226_r8,  0.0005432_r8,  0.0006918_r8, &
       0.0008740_r8,  0.0010966_r8,  0.0013683_r8,  0.0016988_r8,  0.0020991_r8, &
       0.0025806_r8,  0.0031545_r8,  0.0038313_r8,  0.0046201_r8,  0.0055285_r8, &
       0.0065630_r8,  0.0077294_r8,  0.0090332_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000018_r8,  0.0000025_r8,  0.0000036_r8,  0.0000050_r8,  0.0000071_r8, &
       0.0000100_r8,  0.0000140_r8,  0.0000197_r8,  0.0000275_r8,  0.0000384_r8, &
       0.0000534_r8,  0.0000737_r8,  0.0001011_r8,  0.0001373_r8,  0.0001846_r8, &
       0.0002453_r8,  0.0003222_r8,  0.0004185_r8,  0.0005265_r8,  0.0006850_r8, &
       0.0008652_r8,  0.0010855_r8,  0.0013541_r8,  0.0016809_r8,  0.0020768_r8, &
       0.0025528_r8,  0.0031202_r8,  0.0037892_r8,  0.0045688_r8,  0.0054664_r8, &
       0.0064883_r8,  0.0076402_r8,  0.0089277_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000009_r8,  0.0000013_r8, &
       0.0000018_r8,  0.0000025_r8,  0.0000035_r8,  0.0000050_r8,  0.0000070_r8, &
       0.0000098_r8,  0.0000138_r8,  0.0000194_r8,  0.0000272_r8,  0.0000380_r8, &
       0.0000528_r8,  0.0000729_r8,  0.0000999_r8,  0.0001357_r8,  0.0001825_r8, &
       0.0002425_r8,  0.0003185_r8,  0.0004137_r8,  0.0005316_r8,  0.0006769_r8, &
       0.0008548_r8,  0.0010722_r8,  0.0013373_r8,  0.0016599_r8,  0.0020504_r8, &
       0.0025201_r8,  0.0030799_r8,  0.0037398_r8,  0.0045087_r8,  0.0053938_r8, &
       0.0064013_r8,  0.0075366_r8,  0.0088053_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8, &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000009_r8,  0.0000012_r8, &
       0.0000017_r8,  0.0000025_r8,  0.0000035_r8,  0.0000049_r8,  0.0000069_r8, &
       0.0000097_r8,  0.0000137_r8,  0.0000192_r8,  0.0000268_r8,  0.0000375_r8, &
       0.0000520_r8,  0.0000719_r8,  0.0000986_r8,  0.0001339_r8,  0.0001800_r8, &
       0.0002392_r8,  0.0003142_r8,  0.0004079_r8,  0.0005242_r8,  0.0006673_r8, &
       0.0008426_r8,  0.0010567_r8,  0.0013177_r8,  0.0016352_r8,  0.0020196_r8, &
       0.0024820_r8,  0.0030330_r8,  0.0036825_r8,  0.0044391_r8,  0.0053098_r8, &
       0.0063007_r8,  0.0074172_r8,  0.0084815_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000009_r8,  0.0000012_r8,  &
       0.0000017_r8,  0.0000024_r8,  0.0000034_r8,  0.0000048_r8,  0.0000068_r8,  &
       0.0000096_r8,  0.0000134_r8,  0.0000189_r8,  0.0000264_r8,  0.0000369_r8,  &
       0.0000512_r8,  0.0000708_r8,  0.0000970_r8,  0.0001318_r8,  0.0001772_r8,  &
       0.0002354_r8,  0.0003091_r8,  0.0004013_r8,  0.0005156_r8,  0.0006562_r8,  &
       0.0008284_r8,  0.0010386_r8,  0.0012949_r8,  0.0016066_r8,  0.0019840_r8,  &
       0.0024379_r8,  0.0029788_r8,  0.0036164_r8,  0.0043590_r8,  0.0052135_r8,  &
       0.0061857_r8,  0.0072808_r8,  0.0085042_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000008_r8,  0.0000012_r8,  &
       0.0000017_r8,  0.0000024_r8,  0.0000034_r8,  0.0000047_r8,  0.0000067_r8,  &
       0.0000094_r8,  0.0000132_r8,  0.0000185_r8,  0.0000259_r8,  0.0000362_r8,  &
       0.0000503_r8,  0.0000695_r8,  0.0000952_r8,  0.0001294_r8,  0.0001739_r8,  &
       0.0002310_r8,  0.0003033_r8,  0.0003937_r8,  0.0005057_r8,  0.0006435_r8,  &
       0.0008121_r8,  0.0010180_r8,  0.0012688_r8,  0.0015739_r8,  0.0019434_r8,  &
       0.0023877_r8,  0.0029172_r8,  0.0035413_r8,  0.0042681_r8,  0.0051043_r8,  &
       0.0060554_r8,  0.0071267_r8,  0.0083234_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000008_r8,  0.0000012_r8,  &
       0.0000016_r8,  0.0000023_r8,  0.0000033_r8,  0.0000046_r8,  0.0000065_r8,  &
       0.0000092_r8,  0.0000129_r8,  0.0000181_r8,  0.0000254_r8,  0.0000355_r8,  &
       0.0000493_r8,  0.0000680_r8,  0.0000933_r8,  0.0001267_r8,  0.0001702_r8,  &
       0.0002261_r8,  0.0002968_r8,  0.0003852_r8,  0.0004946_r8,  0.0006291_r8,  &
       0.0007937_r8,  0.0009946_r8,  0.0012394_r8,  0.0015370_r8,  0.0018975_r8,  &
       0.0023310_r8,  0.0028478_r8,  0.0034568_r8,  0.0041660_r8,  0.0049818_r8,  &
       0.0059096_r8,  0.0069544_r8,  0.0081215_r8, &
       0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000008_r8,  0.0000011_r8,  &
       0.0000016_r8,  0.0000023_r8,  0.0000032_r8,  0.0000045_r8,  0.0000064_r8,  &
       0.0000090_r8,  0.0000126_r8,  0.0000177_r8,  0.0000248_r8,  0.0000346_r8,  &
       0.0000481_r8,  0.0000664_r8,  0.0000910_r8,  0.0001236_r8,  0.0001661_r8,  &
       0.0002206_r8,  0.0002895_r8,  0.0003755_r8,  0.0004821_r8,  0.0006130_r8,  &
       0.0007731_r8,  0.0009685_r8,  0.0012065_r8,  0.0014959_r8,  0.0018463_r8,  &
       0.0022680_r8,  0.0027705_r8,  0.0033629_r8,  0.0040526_r8,  0.0048459_r8,  &
       0.0057480_r8,  0.0067639_r8,  0.0078987_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000008_r8,  0.0000011_r8,  &
       0.0000016_r8,  0.0000022_r8,  0.0000031_r8,  0.0000044_r8,  0.0000062_r8,  &
       0.0000087_r8,  0.0000123_r8,  0.0000173_r8,  0.0000242_r8,  0.0000330_r8,  &
       0.0000468_r8,  0.0000646_r8,  0.0000886_r8,  0.0001203_r8,  0.0001616_r8,  &
       0.0002145_r8,  0.0002814_r8,  0.0003649_r8,  0.0004682_r8,  0.0005951_r8,  &
       0.0007503_r8,  0.0009396_r8,  0.0011701_r8,  0.0014505_r8,  0.0017900_r8,  &
       0.0021986_r8,  0.0026857_r8,  0.0032598_r8,  0.0039283_r8,  0.0046971_r8,  &
       0.0055713_r8,  0.0065558_r8,  0.0076557_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000005_r8,  0.0000008_r8,  0.0000011_r8,  &
       0.0000015_r8,  0.0000021_r8,  0.0000030_r8,  0.0000043_r8,  0.0000060_r8,  &
       0.0000085_r8,  0.0000119_r8,  0.0000167_r8,  0.0000234_r8,  0.0000327_r8,  &
       0.0000454_r8,  0.0000627_r8,  0.0000859_r8,  0.0001166_r8,  0.0001566_r8,  &
       0.0002078_r8,  0.0002724_r8,  0.0003531_r8,  0.0004529_r8,  0.0005755_r8,  &
       0.0007253_r8,  0.0009079_r8,  0.0011304_r8,  0.0014010_r8,  0.0017287_r8,  &
       0.0021232_r8,  0.0025935_r8,  0.0031480_r8,  0.0037936_r8,  0.0045361_r8,  &
       0.0053805_r8,  0.0063314_r8,  0.0073941_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000005_r8,  0.0000007_r8,  0.0000010_r8,  &
       0.0000015_r8,  0.0000021_r8,  0.0000029_r8,  0.0000041_r8,  0.0000058_r8,  &
       0.0000082_r8,  0.0000115_r8,  0.0000162_r8,  0.0000226_r8,  0.0000316_r8,  &
       0.0000438_r8,  0.0000605_r8,  0.0000829_r8,  0.0001125_r8,  0.0001510_r8,  &
       0.0002004_r8,  0.0002626_r8,  0.0003402_r8,  0.0004362_r8,  0.0005540_r8,  &
       0.0006980_r8,  0.0008736_r8,  0.0010874_r8,  0.0013476_r8,  0.0016627_r8,  &
       0.0020421_r8,  0.0024947_r8,  0.0030283_r8,  0.0036497_r8,  0.0043644_r8,  &
       0.0051772_r8,  0.0060928_r8,  0.0071164_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000003_r8,  0.0000004_r8,  0.0000005_r8,  0.0000007_r8,  0.0000010_r8,  &
       0.0000014_r8,  0.0000020_r8,  0.0000028_r8,  0.0000040_r8,  0.0000056_r8,  &
       0.0000079_r8,  0.0000111_r8,  0.0000155_r8,  0.0000218_r8,  0.0000303_r8,  &
       0.0000421_r8,  0.0000582_r8,  0.0000797_r8,  0.0001081_r8,  0.0001450_r8,  &
       0.0001923_r8,  0.0002519_r8,  0.0003262_r8,  0.0004180_r8,  0.0005308_r8,  &
       0.0006686_r8,  0.0008367_r8,  0.0010414_r8,  0.0012905_r8,  0.0015925_r8,  &
       0.0019561_r8,  0.0023900_r8,  0.0029017_r8,  0.0034978_r8,  0.0041836_r8,  &
       0.0049638_r8,  0.0058430_r8,  0.0068264_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8, &
       0.0000002_r8,  0.0000003_r8,  0.0000005_r8,  0.0000007_r8,  0.0000010_r8, &
       0.0000014_r8,  0.0000019_r8,  0.0000027_r8,  0.0000038_r8,  0.0000053_r8, &
       0.0000075_r8,  0.0000106_r8,  0.0000149_r8,  0.0000208_r8,  0.0000290_r8, &
       0.0000403_r8,  0.0000556_r8,  0.0000761_r8,  0.0001032_r8,  0.0001384_r8, &
       0.0001834_r8,  0.0002402_r8,  0.0003110_r8,  0.0003985_r8,  0.0005059_r8, &
       0.0006372_r8,  0.0007974_r8,  0.0009926_r8,  0.0012302_r8,  0.0015185_r8, &
       0.0018657_r8,  0.0022803_r8,  0.0027696_r8,  0.0033398_r8,  0.0039960_r8, &
       0.0047430_r8,  0.0055851_r8,  0.0065278_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000002_r8,  0.0000003_r8,  0.0000005_r8,  0.0000006_r8,  0.0000009_r8,  &
       0.0000013_r8,  0.0000018_r8,  0.0000026_r8,  0.0000036_r8,  0.0000051_r8,  &
       0.0000071_r8,  0.0000100_r8,  0.0000141_r8,  0.0000197_r8,  0.0000275_r8,  &
       0.0000382_r8,  0.0000527_r8,  0.0000722_r8,  0.0000979_r8,  0.0001312_r8,  &
       0.0001739_r8,  0.0002277_r8,  0.0002947_r8,  0.0003775_r8,  0.0004793_r8,  &
       0.0006038_r8,  0.0007558_r8,  0.0009412_r8,  0.0011671_r8,  0.0014412_r8,  &
       0.0017717_r8,  0.0021208_r8,  0.0026329_r8,  0.0031768_r8,  0.0038033_r8,  &
       0.0045168_r8,  0.0053220_r8,  0.0062240_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000002_r8,  &
       0.0000002_r8,  0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000009_r8,  &
       0.0000012_r8,  0.0000017_r8,  0.0000024_r8,  0.0000034_r8,  0.0000048_r8,  &
       0.0000067_r8,  0.0000095_r8,  0.0000133_r8,  0.0000186_r8,  0.0000259_r8,  &
       0.0000360_r8,  0.0000496_r8,  0.0000679_r8,  0.0000921_r8,  0.0001235_r8,  &
       0.0001637_r8,  0.0002143_r8,  0.0002773_r8,  0.0003554_r8,  0.0004513_r8,  &
       0.0005688_r8,  0.0007124_r8,  0.0008876_r8,  0.0011014_r8,  0.0013610_r8,  &
       0.0016745_r8,  0.0020493_r8,  0.0024925_r8,  0.0030099_r8,  0.0036066_r8,  &
       0.0042868_r8,  0.0050553_r8,  0.0059171_r8, &
       0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  &
       0.0000002_r8,  0.0000003_r8,  0.0000004_r8,  0.0000006_r8,  0.0000008_r8,  &
       0.0000011_r8,  0.0000016_r8,  0.0000022_r8,  0.0000032_r8,  0.0000045_r8,  &
       0.0000063_r8,  0.0000088_r8,  0.0000124_r8,  0.0000173_r8,  0.0000242_r8,  &
       0.0000336_r8,  0.0000463_r8,  0.0000634_r8,  0.0000860_r8,  0.0001153_r8,  &
       0.0001528_r8,  0.0002001_r8,  0.0002591_r8,  0.0003322_r8,  0.0004221_r8,  &
       0.0005323_r8,  0.0006672_r8,  0.0008322_r8,  0.0010335_r8,  0.0012785_r8,  &
       0.0015746_r8,  0.0019293_r8,  0.0023491_r8,  0.0028399_r8,  0.0034067_r8,  &
       0.0040539_r8,  0.0047860_r8,  0.0056083_r8, &
       0.0000000_r8,  0.0000000_r8,  0.0000001_r8,  0.0000001_r8,  0.0000001_r8,  &
       0.0000002_r8,  0.0000003_r8,  0.0000004_r8,  0.0000005_r8,  0.0000007_r8,  &
       0.0000010_r8,  0.0000015_r8,  0.0000021_r8,  0.0000029_r8,  0.0000041_r8,  &
       0.0000058_r8,  0.0000082_r8,  0.0000114_r8,  0.0000160_r8,  0.0000223_r8,  &
       0.0000310_r8,  0.0000428_r8,  0.0000586_r8,  0.0000795_r8,  0.0001067_r8,  &
       0.0001414_r8,  0.0001853_r8,  0.0002401_r8,  0.0003081_r8,  0.0003918_r8,  &
       0.0004947_r8,  0.0006208_r8,  0.0007751_r8,  0.0009639_r8,  0.0011940_r8,  &
       0.0014726_r8,  0.0018069_r8,  0.0022032_r8,  0.0026674_r8,  0.0032043_r8,  &
       0.0038186_r8,  0.0045147_r8,  0.0052979_r8/),SHAPE=(/nu,nw/))




   END SUBROUTINE read_table

 
!***********************************************************************
      subroutine swcldpr()
!***********************************************************************

! Purpose: Define cloud extinction coefficient, single scattering albedo
!          and asymmetry parameter data.
!

! ------- Modules -------

!      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
!                           extice2, ssaice2, asyice2, &
!                           extice3, ssaice3, asyice3, fdlice3, &
!                           abari, bbari, cbari, dbari, ebari, fbari
!
!      save

!-----------------------------------------------------------------------
!
! Explanation of the method for each value of INFLAG.  A value of
!  0 for INFLAG do not distingish being liquid and ice clouds.
!  INFLAG = 2 does distinguish between liquid and ice clouds, and
!    requires further user input to specify the method to be used to 
!    compute the aborption due to each.
!  INFLAG = 0:  For each cloudy layer, the cloud fraction, the cloud optical
!    depth, the cloud single-scattering albedo, and the
!    moments of the phase function (0:NSTREAM).  Note
!    that these values are delta-m scaled within this
!    subroutine.

!  INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
!    water path (g/m2), and cloud ice fraction are input.
!  ICEFLAG = 2:  The ice effective radius (microns) is input and the
!    optical properties due to ice clouds are computed from
!    the optical properties stored in the RT code, STREAMER v3.0 
!    (Reference: Key. J., Streamer User's Guide, Cooperative 
!    Institute for Meteorological Satellite Studies, 2001, 96 pp.).
!    Valid range of values for re are between 5.0 and
!    131.0 micron.
!    This version uses Ebert and Curry, JGR, (1992) method for 
!    ice particles larger than 131.0 microns. 
!  ICEFLAG = 3:  The ice generalized effective size (dge) is input
!    and the optical depths, single-scattering albedo,
!    and phase function moments are calculated as in
!    Q. Fu, J. Climate, (1996). Q. Fu provided high resolution
!    tables which were appropriately averaged for the
!    bands in RRTM_SW.  Linear interpolation is used to
!    get the coefficients from the stored tables.
!    Valid range of values for dge are between 5.0 and
!    140.0 micron. 
!    This version uses Ebert and Curry, JGR, (1992) method for 
!    ice particles larger than 140.0 microns. 
!  LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!    and the optical depths due to water clouds are computed 
!    as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!    The values for absorption coefficients appropriate for
!    the spectral bands in RRTM have been obtained for a 
!    range of effective radii by an averaging procedure 
!    based on the work of J. Pinto (private communication).
!    Linear interpolation is used to get the absorption 
!    coefficients for the input effective radius.
!
!     ------------------------------------------------------------------

! Everything below is for INFLAG = 2.

! Coefficients for Ebert and Curry method
      abari(:) = (/ &
        & 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8,3.448e-03_r8,3.448e-03_r8 /)
      bbari(:) = (/ &
        & 2.431e+00_r8,2.431e+00_r8,2.431e+00_r8,2.431e+00_r8,2.431e+00_r8 /)
      cbari(:) = (/ &
        & 1.000e-05_r8,1.100e-04_r8,1.240e-02_r8,3.779e-02_r8,4.666e-01_r8 /)
      dbari(:) = (/ &
        & 0.000e+00_r8,1.405e-05_r8,6.867e-04_r8,1.284e-03_r8,2.050e-05_r8 /)
      ebari(:) = (/ &
        & 7.661e-01_r8,7.730e-01_r8,7.865e-01_r8,8.172e-01_r8,9.595e-01_r8 /)
      fbari(:) = (/ &
        & 5.851e-04_r8,5.665e-04_r8,7.204e-04_r8,7.463e-04_r8,1.076e-04_r8 /)

! Extinction coefficient
      extliq1(:, 16) = (/ &
        & 8.981463e-01_r8,6.317895e-01_r8,4.557508e-01_r8,3.481624e-01_r8,2.797950e-01_r8,&
        & 2.342753e-01_r8,2.026934e-01_r8,1.800102e-01_r8,1.632408e-01_r8,1.505384e-01_r8,&
        & 1.354524e-01_r8,1.246520e-01_r8,1.154342e-01_r8,1.074756e-01_r8,1.005353e-01_r8,&
        & 9.442987e-02_r8,8.901760e-02_r8,8.418693e-02_r8,7.984904e-02_r8,7.593229e-02_r8,&
        & 7.237827e-02_r8,6.913887e-02_r8,6.617415e-02_r8,6.345061e-02_r8,6.094001e-02_r8,&
        & 5.861834e-02_r8,5.646506e-02_r8,5.446250e-02_r8,5.249596e-02_r8,5.081114e-02_r8,&
        & 4.922243e-02_r8,4.772189e-02_r8,4.630243e-02_r8,4.495766e-02_r8,4.368189e-02_r8,&
        & 4.246995e-02_r8,4.131720e-02_r8,4.021941e-02_r8,3.917276e-02_r8,3.817376e-02_r8,&
        & 3.721926e-02_r8,3.630635e-02_r8,3.543237e-02_r8,3.459491e-02_r8,3.379171e-02_r8,&
        & 3.302073e-02_r8,3.228007e-02_r8,3.156798e-02_r8,3.088284e-02_r8,3.022315e-02_r8,&
        & 2.958753e-02_r8,2.897468e-02_r8,2.838340e-02_r8,2.781258e-02_r8,2.726117e-02_r8,&
        & 2.672821e-02_r8,2.621278e-02_r8,2.5714e-02_r8 /)
      extliq1(:, 17) = (/ &
        & 8.293797e-01_r8,6.048371e-01_r8,4.465706e-01_r8,3.460387e-01_r8,2.800064e-01_r8,&
        & 2.346584e-01_r8,2.022399e-01_r8,1.782626e-01_r8,1.600153e-01_r8,1.457903e-01_r8,&
        & 1.334061e-01_r8,1.228548e-01_r8,1.138396e-01_r8,1.060486e-01_r8,9.924856e-02_r8,&
        & 9.326208e-02_r8,8.795158e-02_r8,8.320883e-02_r8,7.894750e-02_r8,7.509792e-02_r8,&
        & 7.160323e-02_r8,6.841653e-02_r8,6.549889e-02_r8,6.281763e-02_r8,6.034516e-02_r8,&
        & 5.805802e-02_r8,5.593615e-02_r8,5.396226e-02_r8,5.202302e-02_r8,5.036246e-02_r8,&
        & 4.879606e-02_r8,4.731610e-02_r8,4.591565e-02_r8,4.458852e-02_r8,4.332912e-02_r8,&
        & 4.213243e-02_r8,4.099390e-02_r8,3.990941e-02_r8,3.887522e-02_r8,3.788792e-02_r8,&
        & 3.694440e-02_r8,3.604183e-02_r8,3.517760e-02_r8,3.434934e-02_r8,3.355485e-02_r8,&
        & 3.279211e-02_r8,3.205925e-02_r8,3.135458e-02_r8,3.067648e-02_r8,3.002349e-02_r8,&
        & 2.939425e-02_r8,2.878748e-02_r8,2.820200e-02_r8,2.763673e-02_r8,2.709062e-02_r8,&
        & 2.656272e-02_r8,2.605214e-02_r8,2.5558e-02_r8 /)
      extliq1(:, 18) = (/ &
        & 9.193685e-01_r8,6.128292e-01_r8,4.344150e-01_r8,3.303048e-01_r8,2.659500e-01_r8,&
        & 2.239727e-01_r8,1.953457e-01_r8,1.751012e-01_r8,1.603515e-01_r8,1.493360e-01_r8,&
        & 1.323791e-01_r8,1.219335e-01_r8,1.130076e-01_r8,1.052926e-01_r8,9.855839e-02_r8,&
        & 9.262925e-02_r8,8.736918e-02_r8,8.267112e-02_r8,7.844965e-02_r8,7.463585e-02_r8,&
        & 7.117343e-02_r8,6.801601e-02_r8,6.512503e-02_r8,6.246815e-02_r8,6.001806e-02_r8,&
        & 5.775154e-02_r8,5.564872e-02_r8,5.369250e-02_r8,5.176284e-02_r8,5.011536e-02_r8,&
        & 4.856099e-02_r8,4.709211e-02_r8,4.570193e-02_r8,4.438430e-02_r8,4.313375e-02_r8,&
        & 4.194529e-02_r8,4.081443e-02_r8,3.973712e-02_r8,3.870966e-02_r8,3.772866e-02_r8,&
        & 3.679108e-02_r8,3.589409e-02_r8,3.503514e-02_r8,3.421185e-02_r8,3.342206e-02_r8,&
        & 3.266377e-02_r8,3.193513e-02_r8,3.123447e-02_r8,3.056018e-02_r8,2.991081e-02_r8,&
        & 2.928502e-02_r8,2.868154e-02_r8,2.809920e-02_r8,2.753692e-02_r8,2.699367e-02_r8,&
        & 2.646852e-02_r8,2.596057e-02_r8,2.5469e-02_r8 /)
      extliq1(:, 19) = (/ &
        & 9.136931e-01_r8,5.743244e-01_r8,4.080708e-01_r8,3.150572e-01_r8,2.577261e-01_r8,&
        & 2.197900e-01_r8,1.933037e-01_r8,1.740212e-01_r8,1.595056e-01_r8,1.482756e-01_r8,&
        & 1.312164e-01_r8,1.209246e-01_r8,1.121227e-01_r8,1.045095e-01_r8,9.785967e-02_r8,&
        & 9.200149e-02_r8,8.680170e-02_r8,8.215531e-02_r8,7.797850e-02_r8,7.420361e-02_r8,&
        & 7.077530e-02_r8,6.764798e-02_r8,6.478369e-02_r8,6.215063e-02_r8,5.972189e-02_r8,&
        & 5.747458e-02_r8,5.538913e-02_r8,5.344866e-02_r8,5.153216e-02_r8,4.989745e-02_r8,&
        & 4.835476e-02_r8,4.689661e-02_r8,4.551629e-02_r8,4.420777e-02_r8,4.296563e-02_r8,&
        & 4.178497e-02_r8,4.066137e-02_r8,3.959081e-02_r8,3.856963e-02_r8,3.759452e-02_r8,&
        & 3.666244e-02_r8,3.577061e-02_r8,3.491650e-02_r8,3.409777e-02_r8,3.331227e-02_r8,&
        & 3.255803e-02_r8,3.183322e-02_r8,3.113617e-02_r8,3.046530e-02_r8,2.981918e-02_r8,&
        & 2.919646e-02_r8,2.859591e-02_r8,2.801635e-02_r8,2.745671e-02_r8,2.691599e-02_r8,&
        & 2.639324e-02_r8,2.588759e-02_r8,2.5398e-02_r8 /)
      extliq1(:, 20) = (/ &
        & 8.447548e-01_r8,5.326840e-01_r8,3.921523e-01_r8,3.119082e-01_r8,2.597055e-01_r8,&
        & 2.228737e-01_r8,1.954157e-01_r8,1.741155e-01_r8,1.570881e-01_r8,1.431520e-01_r8,&
        & 1.302034e-01_r8,1.200491e-01_r8,1.113571e-01_r8,1.038330e-01_r8,9.725657e-02_r8,&
        & 9.145949e-02_r8,8.631112e-02_r8,8.170840e-02_r8,7.756901e-02_r8,7.382641e-02_r8,&
        & 7.042616e-02_r8,6.732338e-02_r8,6.448069e-02_r8,6.186672e-02_r8,5.945494e-02_r8,&
        & 5.722277e-02_r8,5.515089e-02_r8,5.322262e-02_r8,5.132153e-02_r8,4.969799e-02_r8,&
        & 4.816556e-02_r8,4.671686e-02_r8,4.534525e-02_r8,4.404480e-02_r8,4.281014e-02_r8,&
        & 4.163643e-02_r8,4.051930e-02_r8,3.945479e-02_r8,3.843927e-02_r8,3.746945e-02_r8,&
        & 3.654234e-02_r8,3.565518e-02_r8,3.480547e-02_r8,3.399088e-02_r8,3.320930e-02_r8,&
        & 3.245876e-02_r8,3.173745e-02_r8,3.104371e-02_r8,3.037600e-02_r8,2.973287e-02_r8,&
        & 2.911300e-02_r8,2.851516e-02_r8,2.793818e-02_r8,2.738101e-02_r8,2.684264e-02_r8,&
        & 2.632214e-02_r8,2.581863e-02_r8,2.5331e-02_r8 /)
      extliq1(:, 21) = (/ &
        & 7.727642e-01_r8,5.034865e-01_r8,3.808673e-01_r8,3.080333e-01_r8,2.586453e-01_r8,&
        & 2.224989e-01_r8,1.947060e-01_r8,1.725821e-01_r8,1.545096e-01_r8,1.394456e-01_r8,&
        & 1.288683e-01_r8,1.188852e-01_r8,1.103317e-01_r8,1.029214e-01_r8,9.643967e-02_r8,&
        & 9.072239e-02_r8,8.564194e-02_r8,8.109758e-02_r8,7.700875e-02_r8,7.331026e-02_r8,&
        & 6.994879e-02_r8,6.688028e-02_r8,6.406807e-02_r8,6.148133e-02_r8,5.909400e-02_r8,&
        & 5.688388e-02_r8,5.483197e-02_r8,5.292185e-02_r8,5.103763e-02_r8,4.942905e-02_r8,&
        & 4.791039e-02_r8,4.647438e-02_r8,4.511453e-02_r8,4.382497e-02_r8,4.260043e-02_r8,&
        & 4.143616e-02_r8,4.032784e-02_r8,3.927155e-02_r8,3.826375e-02_r8,3.730117e-02_r8,&
        & 3.638087e-02_r8,3.550013e-02_r8,3.465646e-02_r8,3.384759e-02_r8,3.307141e-02_r8,&
        & 3.232598e-02_r8,3.160953e-02_r8,3.092040e-02_r8,3.025706e-02_r8,2.961810e-02_r8,&
        & 2.900220e-02_r8,2.840814e-02_r8,2.783478e-02_r8,2.728106e-02_r8,2.674599e-02_r8,&
        & 2.622864e-02_r8,2.572816e-02_r8,2.5244e-02_r8 /)
      extliq1(:, 22) = (/ &
        & 7.416833e-01_r8,4.959591e-01_r8,3.775057e-01_r8,3.056353e-01_r8,2.565943e-01_r8,&
        & 2.206935e-01_r8,1.931479e-01_r8,1.712860e-01_r8,1.534837e-01_r8,1.386906e-01_r8,&
        & 1.281198e-01_r8,1.182344e-01_r8,1.097595e-01_r8,1.024137e-01_r8,9.598552e-02_r8,&
        & 9.031320e-02_r8,8.527093e-02_r8,8.075927e-02_r8,7.669869e-02_r8,7.302481e-02_r8,&
        & 6.968491e-02_r8,6.663542e-02_r8,6.384008e-02_r8,6.126838e-02_r8,5.889452e-02_r8,&
        & 5.669654e-02_r8,5.465558e-02_r8,5.275540e-02_r8,5.087937e-02_r8,4.927904e-02_r8,&
        & 4.776796e-02_r8,4.633895e-02_r8,4.498557e-02_r8,4.370202e-02_r8,4.248306e-02_r8,&
        & 4.132399e-02_r8,4.022052e-02_r8,3.916878e-02_r8,3.816523e-02_r8,3.720665e-02_r8,&
        & 3.629011e-02_r8,3.541290e-02_r8,3.457257e-02_r8,3.376685e-02_r8,3.299365e-02_r8,&
        & 3.225105e-02_r8,3.153728e-02_r8,3.085069e-02_r8,3.018977e-02_r8,2.955310e-02_r8,&
        & 2.893940e-02_r8,2.834742e-02_r8,2.777606e-02_r8,2.722424e-02_r8,2.669099e-02_r8,&
        & 2.617539e-02_r8,2.567658e-02_r8,2.5194e-02_r8 /)
      extliq1(:, 23) = (/ &
        & 7.058580e-01_r8,4.866573e-01_r8,3.712238e-01_r8,2.998638e-01_r8,2.513441e-01_r8,&
        & 2.161972e-01_r8,1.895576e-01_r8,1.686669e-01_r8,1.518437e-01_r8,1.380046e-01_r8,&
        & 1.267564e-01_r8,1.170399e-01_r8,1.087026e-01_r8,1.014704e-01_r8,9.513729e-02_r8,&
        & 8.954555e-02_r8,8.457221e-02_r8,8.012009e-02_r8,7.611136e-02_r8,7.248294e-02_r8,&
        & 6.918317e-02_r8,6.616934e-02_r8,6.340584e-02_r8,6.086273e-02_r8,5.851465e-02_r8,&
        & 5.634001e-02_r8,5.432027e-02_r8,5.243946e-02_r8,5.058070e-02_r8,4.899628e-02_r8,&
        & 4.749975e-02_r8,4.608411e-02_r8,4.474303e-02_r8,4.347082e-02_r8,4.226237e-02_r8,&
        & 4.111303e-02_r8,4.001861e-02_r8,3.897528e-02_r8,3.797959e-02_r8,3.702835e-02_r8,&
        & 3.611867e-02_r8,3.524791e-02_r8,3.441364e-02_r8,3.361360e-02_r8,3.284577e-02_r8,&
        & 3.210823e-02_r8,3.139923e-02_r8,3.071716e-02_r8,3.006052e-02_r8,2.942791e-02_r8,&
        & 2.881806e-02_r8,2.822974e-02_r8,2.766185e-02_r8,2.711335e-02_r8,2.658326e-02_r8,&
        & 2.607066e-02_r8,2.557473e-02_r8,2.5095e-02_r8 /)
      extliq1(:, 24) = (/ &
        & 6.822779e-01_r8,4.750373e-01_r8,3.634834e-01_r8,2.940726e-01_r8,2.468060e-01_r8,&
        & 2.125768e-01_r8,1.866586e-01_r8,1.663588e-01_r8,1.500326e-01_r8,1.366192e-01_r8,&
        & 1.253472e-01_r8,1.158052e-01_r8,1.076101e-01_r8,1.004954e-01_r8,9.426089e-02_r8,&
        & 8.875268e-02_r8,8.385090e-02_r8,7.946063e-02_r8,7.550578e-02_r8,7.192466e-02_r8,&
        & 6.866669e-02_r8,6.569001e-02_r8,6.295971e-02_r8,6.044642e-02_r8,5.812526e-02_r8,&
        & 5.597500e-02_r8,5.397746e-02_r8,5.211690e-02_r8,5.027505e-02_r8,4.870703e-02_r8,&
        & 4.722555e-02_r8,4.582373e-02_r8,4.449540e-02_r8,4.323497e-02_r8,4.203742e-02_r8,&
        & 4.089821e-02_r8,3.981321e-02_r8,3.877867e-02_r8,3.779118e-02_r8,3.684762e-02_r8,&
        & 3.594514e-02_r8,3.508114e-02_r8,3.425322e-02_r8,3.345917e-02_r8,3.269698e-02_r8,&
        & 3.196477e-02_r8,3.126082e-02_r8,3.058352e-02_r8,2.993141e-02_r8,2.930310e-02_r8,&
        & 2.869732e-02_r8,2.811289e-02_r8,2.754869e-02_r8,2.700371e-02_r8,2.647698e-02_r8,&
        & 2.596760e-02_r8,2.547473e-02_r8,2.4998e-02_r8 /)
      extliq1(:, 25) = (/ &
        & 6.666233e-01_r8,4.662044e-01_r8,3.579517e-01_r8,2.902984e-01_r8,2.440475e-01_r8,&
        & 2.104431e-01_r8,1.849277e-01_r8,1.648970e-01_r8,1.487555e-01_r8,1.354714e-01_r8,&
        & 1.244173e-01_r8,1.149913e-01_r8,1.068903e-01_r8,9.985323e-02_r8,9.368351e-02_r8,&
        & 8.823009e-02_r8,8.337507e-02_r8,7.902511e-02_r8,7.510529e-02_r8,7.155482e-02_r8,&
        & 6.832386e-02_r8,6.537113e-02_r8,6.266218e-02_r8,6.016802e-02_r8,5.786408e-02_r8,&
        & 5.572939e-02_r8,5.374598e-02_r8,5.189830e-02_r8,5.006825e-02_r8,4.851081e-02_r8,&
        & 4.703906e-02_r8,4.564623e-02_r8,4.432621e-02_r8,4.307349e-02_r8,4.188312e-02_r8,&
        & 4.075060e-02_r8,3.967183e-02_r8,3.864313e-02_r8,3.766111e-02_r8,3.672269e-02_r8,&
        & 3.582505e-02_r8,3.496559e-02_r8,3.414196e-02_r8,3.335198e-02_r8,3.259362e-02_r8,&
        & 3.186505e-02_r8,3.116454e-02_r8,3.049052e-02_r8,2.984152e-02_r8,2.921617e-02_r8,&
        & 2.861322e-02_r8,2.803148e-02_r8,2.746986e-02_r8,2.692733e-02_r8,2.640295e-02_r8,&
        & 2.589582e-02_r8,2.540510e-02_r8,2.4930e-02_r8 /)
      extliq1(:, 26) = (/ &
        & 6.535669e-01_r8,4.585865e-01_r8,3.529226e-01_r8,2.867245e-01_r8,2.413848e-01_r8,&
        & 2.083956e-01_r8,1.833191e-01_r8,1.636150e-01_r8,1.477247e-01_r8,1.346392e-01_r8,&
        & 1.236449e-01_r8,1.143095e-01_r8,1.062828e-01_r8,9.930773e-02_r8,9.319029e-02_r8,&
        & 8.778150e-02_r8,8.296497e-02_r8,7.864847e-02_r8,7.475799e-02_r8,7.123343e-02_r8,&
        & 6.802549e-02_r8,6.509332e-02_r8,6.240285e-02_r8,5.992538e-02_r8,5.763657e-02_r8,&
        & 5.551566e-02_r8,5.354483e-02_r8,5.170870e-02_r8,4.988866e-02_r8,4.834061e-02_r8,&
        & 4.687751e-02_r8,4.549264e-02_r8,4.417999e-02_r8,4.293410e-02_r8,4.175006e-02_r8,&
        & 4.062344e-02_r8,3.955019e-02_r8,3.852663e-02_r8,3.754943e-02_r8,3.661553e-02_r8,&
        & 3.572214e-02_r8,3.486669e-02_r8,3.404683e-02_r8,3.326040e-02_r8,3.250542e-02_r8,&
        & 3.178003e-02_r8,3.108254e-02_r8,3.041139e-02_r8,2.976511e-02_r8,2.914235e-02_r8,&
        & 2.854187e-02_r8,2.796247e-02_r8,2.740309e-02_r8,2.686271e-02_r8,2.634038e-02_r8,&
        & 2.583520e-02_r8,2.534636e-02_r8,2.4873e-02_r8 /)
      extliq1(:, 27) = (/ &
        & 6.448790e-01_r8,4.541425e-01_r8,3.503348e-01_r8,2.850494e-01_r8,2.401966e-01_r8,&
        & 2.074811e-01_r8,1.825631e-01_r8,1.629515e-01_r8,1.471142e-01_r8,1.340574e-01_r8,&
        & 1.231462e-01_r8,1.138628e-01_r8,1.058802e-01_r8,9.894286e-02_r8,9.285818e-02_r8,&
        & 8.747802e-02_r8,8.268676e-02_r8,7.839271e-02_r8,7.452230e-02_r8,7.101580e-02_r8,&
        & 6.782418e-02_r8,6.490685e-02_r8,6.222991e-02_r8,5.976484e-02_r8,5.748742e-02_r8,&
        & 5.537703e-02_r8,5.341593e-02_r8,5.158883e-02_r8,4.977355e-02_r8,4.823172e-02_r8,&
        & 4.677430e-02_r8,4.539465e-02_r8,4.408680e-02_r8,4.284533e-02_r8,4.166539e-02_r8,&
        & 4.054257e-02_r8,3.947283e-02_r8,3.845256e-02_r8,3.747842e-02_r8,3.654737e-02_r8,&
        & 3.565665e-02_r8,3.480370e-02_r8,3.398620e-02_r8,3.320198e-02_r8,3.244908e-02_r8,&
        & 3.172566e-02_r8,3.103002e-02_r8,3.036062e-02_r8,2.971600e-02_r8,2.909482e-02_r8,&
        & 2.849582e-02_r8,2.791785e-02_r8,2.735982e-02_r8,2.682072e-02_r8,2.629960e-02_r8,&
        & 2.579559e-02_r8,2.530786e-02_r8,2.4836e-02_r8 /)
      extliq1(:, 28) = (/ &
        & 6.422688e-01_r8,4.528453e-01_r8,3.497232e-01_r8,2.847724e-01_r8,2.400815e-01_r8,&
        & 2.074403e-01_r8,1.825502e-01_r8,1.629415e-01_r8,1.470934e-01_r8,1.340183e-01_r8,&
        & 1.230935e-01_r8,1.138049e-01_r8,1.058201e-01_r8,9.888245e-02_r8,9.279878e-02_r8,&
        & 8.742053e-02_r8,8.263175e-02_r8,7.834058e-02_r8,7.447327e-02_r8,7.097000e-02_r8,&
        & 6.778167e-02_r8,6.486765e-02_r8,6.219400e-02_r8,5.973215e-02_r8,5.745790e-02_r8,&
        & 5.535059e-02_r8,5.339250e-02_r8,5.156831e-02_r8,4.975308e-02_r8,4.821235e-02_r8,&
        & 4.675596e-02_r8,4.537727e-02_r8,4.407030e-02_r8,4.282968e-02_r8,4.165053e-02_r8,&
        & 4.052845e-02_r8,3.945941e-02_r8,3.843980e-02_r8,3.746628e-02_r8,3.653583e-02_r8,&
        & 3.564567e-02_r8,3.479326e-02_r8,3.397626e-02_r8,3.319253e-02_r8,3.244008e-02_r8,&
        & 3.171711e-02_r8,3.102189e-02_r8,3.035289e-02_r8,2.970866e-02_r8,2.908784e-02_r8,&
        & 2.848920e-02_r8,2.791156e-02_r8,2.735385e-02_r8,2.681507e-02_r8,2.629425e-02_r8,&
        & 2.579053e-02_r8,2.530308e-02_r8,2.4831e-02_r8 /)
      extliq1(:, 29) = (/ &
        & 4.614710e-01_r8,4.556116e-01_r8,4.056568e-01_r8,3.529833e-01_r8,3.060334e-01_r8,&
        & 2.658127e-01_r8,2.316095e-01_r8,2.024325e-01_r8,1.773749e-01_r8,1.556867e-01_r8,&
        & 1.455558e-01_r8,1.332882e-01_r8,1.229052e-01_r8,1.140067e-01_r8,1.062981e-01_r8,&
        & 9.955703e-02_r8,9.361333e-02_r8,8.833420e-02_r8,8.361467e-02_r8,7.937071e-02_r8,&
        & 7.553420e-02_r8,7.204942e-02_r8,6.887031e-02_r8,6.595851e-02_r8,6.328178e-02_r8,&
        & 6.081286e-02_r8,5.852854e-02_r8,5.640892e-02_r8,5.431269e-02_r8,5.252561e-02_r8,&
        & 5.084345e-02_r8,4.925727e-02_r8,4.775910e-02_r8,4.634182e-02_r8,4.499907e-02_r8,&
        & 4.372512e-02_r8,4.251484e-02_r8,4.136357e-02_r8,4.026710e-02_r8,3.922162e-02_r8,&
        & 3.822365e-02_r8,3.727004e-02_r8,3.635790e-02_r8,3.548457e-02_r8,3.464764e-02_r8,&
        & 3.384488e-02_r8,3.307424e-02_r8,3.233384e-02_r8,3.162192e-02_r8,3.093688e-02_r8,&
        & 3.027723e-02_r8,2.964158e-02_r8,2.902864e-02_r8,2.843722e-02_r8,2.786621e-02_r8,&
        & 2.731457e-02_r8,2.678133e-02_r8,2.6266e-02_r8 /)

! Single scattering albedo     
      ssaliq1(:, 16) = (/ &
        & 8.143821e-01_r8,7.836739e-01_r8,7.550722e-01_r8,7.306269e-01_r8,7.105612e-01_r8,&
        & 6.946649e-01_r8,6.825556e-01_r8,6.737762e-01_r8,6.678448e-01_r8,6.642830e-01_r8,&
        & 6.679741e-01_r8,6.584607e-01_r8,6.505598e-01_r8,6.440951e-01_r8,6.388901e-01_r8,&
        & 6.347689e-01_r8,6.315549e-01_r8,6.290718e-01_r8,6.271432e-01_r8,6.255928e-01_r8,&
        & 6.242441e-01_r8,6.229207e-01_r8,6.214464e-01_r8,6.196445e-01_r8,6.173388e-01_r8,&
        & 6.143527e-01_r8,6.105099e-01_r8,6.056339e-01_r8,6.108290e-01_r8,6.073939e-01_r8,&
        & 6.043073e-01_r8,6.015473e-01_r8,5.990913e-01_r8,5.969173e-01_r8,5.950028e-01_r8,&
        & 5.933257e-01_r8,5.918636e-01_r8,5.905944e-01_r8,5.894957e-01_r8,5.885453e-01_r8,&
        & 5.877209e-01_r8,5.870003e-01_r8,5.863611e-01_r8,5.857811e-01_r8,5.852381e-01_r8,&
        & 5.847098e-01_r8,5.841738e-01_r8,5.836081e-01_r8,5.829901e-01_r8,5.822979e-01_r8,&
        & 5.815089e-01_r8,5.806011e-01_r8,5.795521e-01_r8,5.783396e-01_r8,5.769413e-01_r8,&
        & 5.753351e-01_r8,5.734986e-01_r8,5.7141e-01_r8 /)
      ssaliq1(:, 17) = (/ &
        & 8.165821e-01_r8,8.002015e-01_r8,7.816921e-01_r8,7.634131e-01_r8,7.463721e-01_r8,&
        & 7.312469e-01_r8,7.185883e-01_r8,7.088975e-01_r8,7.026671e-01_r8,7.004020e-01_r8,&
        & 7.042138e-01_r8,6.960930e-01_r8,6.894243e-01_r8,6.840459e-01_r8,6.797957e-01_r8,&
        & 6.765119e-01_r8,6.740325e-01_r8,6.721955e-01_r8,6.708391e-01_r8,6.698013e-01_r8,&
        & 6.689201e-01_r8,6.680339e-01_r8,6.669805e-01_r8,6.655982e-01_r8,6.637250e-01_r8,&
        & 6.611992e-01_r8,6.578588e-01_r8,6.535420e-01_r8,6.584449e-01_r8,6.553992e-01_r8,&
        & 6.526547e-01_r8,6.501917e-01_r8,6.479905e-01_r8,6.460313e-01_r8,6.442945e-01_r8,&
        & 6.427605e-01_r8,6.414094e-01_r8,6.402217e-01_r8,6.391775e-01_r8,6.382573e-01_r8,&
        & 6.374413e-01_r8,6.367099e-01_r8,6.360433e-01_r8,6.354218e-01_r8,6.348257e-01_r8,&
        & 6.342355e-01_r8,6.336313e-01_r8,6.329935e-01_r8,6.323023e-01_r8,6.315383e-01_r8,&
        & 6.306814e-01_r8,6.297122e-01_r8,6.286110e-01_r8,6.273579e-01_r8,6.259333e-01_r8,&
        & 6.243176e-01_r8,6.224910e-01_r8,6.2043e-01_r8 /)
      ssaliq1(:, 18) = (/ &
        & 9.900163e-01_r8,9.854307e-01_r8,9.797730e-01_r8,9.733113e-01_r8,9.664245e-01_r8,&
        & 9.594976e-01_r8,9.529055e-01_r8,9.470112e-01_r8,9.421695e-01_r8,9.387304e-01_r8,&
        & 9.344918e-01_r8,9.305302e-01_r8,9.267048e-01_r8,9.230072e-01_r8,9.194289e-01_r8,&
        & 9.159616e-01_r8,9.125968e-01_r8,9.093260e-01_r8,9.061409e-01_r8,9.030330e-01_r8,&
        & 8.999940e-01_r8,8.970154e-01_r8,8.940888e-01_r8,8.912058e-01_r8,8.883579e-01_r8,&
        & 8.855368e-01_r8,8.827341e-01_r8,8.799413e-01_r8,8.777423e-01_r8,8.749566e-01_r8,&
        & 8.722298e-01_r8,8.695605e-01_r8,8.669469e-01_r8,8.643875e-01_r8,8.618806e-01_r8,&
        & 8.594246e-01_r8,8.570179e-01_r8,8.546589e-01_r8,8.523459e-01_r8,8.500773e-01_r8,&
        & 8.478516e-01_r8,8.456670e-01_r8,8.435219e-01_r8,8.414148e-01_r8,8.393439e-01_r8,&
        & 8.373078e-01_r8,8.353047e-01_r8,8.333330e-01_r8,8.313911e-01_r8,8.294774e-01_r8,&
        & 8.275904e-01_r8,8.257282e-01_r8,8.238893e-01_r8,8.220721e-01_r8,8.202751e-01_r8,&
        & 8.184965e-01_r8,8.167346e-01_r8,8.1499e-01_r8 /)
      ssaliq1(:, 19) = (/ &
        & 9.999916e-01_r8,9.987396e-01_r8,9.966900e-01_r8,9.950738e-01_r8,9.937531e-01_r8,&
        & 9.925912e-01_r8,9.914525e-01_r8,9.902018e-01_r8,9.887046e-01_r8,9.868263e-01_r8,&
        & 9.849039e-01_r8,9.832372e-01_r8,9.815265e-01_r8,9.797770e-01_r8,9.779940e-01_r8,&
        & 9.761827e-01_r8,9.743481e-01_r8,9.724955e-01_r8,9.706303e-01_r8,9.687575e-01_r8,&
        & 9.668823e-01_r8,9.650100e-01_r8,9.631457e-01_r8,9.612947e-01_r8,9.594622e-01_r8,&
        & 9.576534e-01_r8,9.558734e-01_r8,9.541275e-01_r8,9.522059e-01_r8,9.504258e-01_r8,&
        & 9.486459e-01_r8,9.468676e-01_r8,9.450921e-01_r8,9.433208e-01_r8,9.415548e-01_r8,&
        & 9.397955e-01_r8,9.380441e-01_r8,9.363022e-01_r8,9.345706e-01_r8,9.328510e-01_r8,&
        & 9.311445e-01_r8,9.294524e-01_r8,9.277761e-01_r8,9.261167e-01_r8,9.244755e-01_r8,&
        & 9.228540e-01_r8,9.212534e-01_r8,9.196748e-01_r8,9.181197e-01_r8,9.165894e-01_r8,&
        & 9.150851e-01_r8,9.136080e-01_r8,9.121596e-01_r8,9.107410e-01_r8,9.093536e-01_r8,&
        & 9.079987e-01_r8,9.066775e-01_r8,9.0539e-01_r8 /)
      ssaliq1(:, 20) = (/ &
        & 9.979493e-01_r8,9.964113e-01_r8,9.950014e-01_r8,9.937045e-01_r8,9.924964e-01_r8,&
        & 9.913546e-01_r8,9.902575e-01_r8,9.891843e-01_r8,9.881136e-01_r8,9.870238e-01_r8,&
        & 9.859934e-01_r8,9.849372e-01_r8,9.838873e-01_r8,9.828434e-01_r8,9.818052e-01_r8,&
        & 9.807725e-01_r8,9.797450e-01_r8,9.787225e-01_r8,9.777047e-01_r8,9.766914e-01_r8,&
        & 9.756823e-01_r8,9.746771e-01_r8,9.736756e-01_r8,9.726775e-01_r8,9.716827e-01_r8,&
        & 9.706907e-01_r8,9.697014e-01_r8,9.687145e-01_r8,9.678060e-01_r8,9.668108e-01_r8,&
        & 9.658218e-01_r8,9.648391e-01_r8,9.638629e-01_r8,9.628936e-01_r8,9.619313e-01_r8,&
        & 9.609763e-01_r8,9.600287e-01_r8,9.590888e-01_r8,9.581569e-01_r8,9.572330e-01_r8,&
        & 9.563176e-01_r8,9.554108e-01_r8,9.545128e-01_r8,9.536239e-01_r8,9.527443e-01_r8,&
        & 9.518741e-01_r8,9.510137e-01_r8,9.501633e-01_r8,9.493230e-01_r8,9.484931e-01_r8,&
        & 9.476740e-01_r8,9.468656e-01_r8,9.460683e-01_r8,9.452824e-01_r8,9.445080e-01_r8,&
        & 9.437454e-01_r8,9.429948e-01_r8,9.4226e-01_r8 /)
      ssaliq1(:, 21) = (/ &
        & 9.988742e-01_r8,9.982668e-01_r8,9.976935e-01_r8,9.971497e-01_r8,9.966314e-01_r8,&
        & 9.961344e-01_r8,9.956545e-01_r8,9.951873e-01_r8,9.947286e-01_r8,9.942741e-01_r8,&
        & 9.938457e-01_r8,9.933947e-01_r8,9.929473e-01_r8,9.925032e-01_r8,9.920621e-01_r8,&
        & 9.916237e-01_r8,9.911875e-01_r8,9.907534e-01_r8,9.903209e-01_r8,9.898898e-01_r8,&
        & 9.894597e-01_r8,9.890304e-01_r8,9.886015e-01_r8,9.881726e-01_r8,9.877435e-01_r8,&
        & 9.873138e-01_r8,9.868833e-01_r8,9.864516e-01_r8,9.860698e-01_r8,9.856317e-01_r8,&
        & 9.851957e-01_r8,9.847618e-01_r8,9.843302e-01_r8,9.839008e-01_r8,9.834739e-01_r8,&
        & 9.830494e-01_r8,9.826275e-01_r8,9.822083e-01_r8,9.817918e-01_r8,9.813782e-01_r8,&
        & 9.809675e-01_r8,9.805598e-01_r8,9.801552e-01_r8,9.797538e-01_r8,9.793556e-01_r8,&
        & 9.789608e-01_r8,9.785695e-01_r8,9.781817e-01_r8,9.777975e-01_r8,9.774171e-01_r8,&
        & 9.770404e-01_r8,9.766676e-01_r8,9.762988e-01_r8,9.759340e-01_r8,9.755733e-01_r8,&
        & 9.752169e-01_r8,9.748649e-01_r8,9.7452e-01_r8 /)
      ssaliq1(:, 22) = (/ &
        & 9.994441e-01_r8,9.991608e-01_r8,9.988949e-01_r8,9.986439e-01_r8,9.984054e-01_r8,&
        & 9.981768e-01_r8,9.979557e-01_r8,9.977396e-01_r8,9.975258e-01_r8,9.973120e-01_r8,&
        & 9.971011e-01_r8,9.968852e-01_r8,9.966708e-01_r8,9.964578e-01_r8,9.962462e-01_r8,&
        & 9.960357e-01_r8,9.958264e-01_r8,9.956181e-01_r8,9.954108e-01_r8,9.952043e-01_r8,&
        & 9.949987e-01_r8,9.947937e-01_r8,9.945892e-01_r8,9.943853e-01_r8,9.941818e-01_r8,&
        & 9.939786e-01_r8,9.937757e-01_r8,9.935728e-01_r8,9.933922e-01_r8,9.931825e-01_r8,&
        & 9.929739e-01_r8,9.927661e-01_r8,9.925592e-01_r8,9.923534e-01_r8,9.921485e-01_r8,&
        & 9.919447e-01_r8,9.917421e-01_r8,9.915406e-01_r8,9.913403e-01_r8,9.911412e-01_r8,&
        & 9.909435e-01_r8,9.907470e-01_r8,9.905519e-01_r8,9.903581e-01_r8,9.901659e-01_r8,&
        & 9.899751e-01_r8,9.897858e-01_r8,9.895981e-01_r8,9.894120e-01_r8,9.892276e-01_r8,&
        & 9.890447e-01_r8,9.888637e-01_r8,9.886845e-01_r8,9.885070e-01_r8,9.883314e-01_r8,&
        & 9.881576e-01_r8,9.879859e-01_r8,9.8782e-01_r8 /)
      ssaliq1(:, 23) = (/ &
        & 9.999138e-01_r8,9.998730e-01_r8,9.998338e-01_r8,9.997965e-01_r8,9.997609e-01_r8,&
        & 9.997270e-01_r8,9.996944e-01_r8,9.996629e-01_r8,9.996321e-01_r8,9.996016e-01_r8,&
        & 9.995690e-01_r8,9.995372e-01_r8,9.995057e-01_r8,9.994744e-01_r8,9.994433e-01_r8,&
        & 9.994124e-01_r8,9.993817e-01_r8,9.993510e-01_r8,9.993206e-01_r8,9.992903e-01_r8,&
        & 9.992600e-01_r8,9.992299e-01_r8,9.991998e-01_r8,9.991698e-01_r8,9.991398e-01_r8,&
        & 9.991098e-01_r8,9.990799e-01_r8,9.990499e-01_r8,9.990231e-01_r8,9.989920e-01_r8,&
        & 9.989611e-01_r8,9.989302e-01_r8,9.988996e-01_r8,9.988690e-01_r8,9.988386e-01_r8,&
        & 9.988084e-01_r8,9.987783e-01_r8,9.987485e-01_r8,9.987187e-01_r8,9.986891e-01_r8,&
        & 9.986598e-01_r8,9.986306e-01_r8,9.986017e-01_r8,9.985729e-01_r8,9.985443e-01_r8,&
        & 9.985160e-01_r8,9.984879e-01_r8,9.984600e-01_r8,9.984324e-01_r8,9.984050e-01_r8,&
        & 9.983778e-01_r8,9.983509e-01_r8,9.983243e-01_r8,9.982980e-01_r8,9.982719e-01_r8,&
        & 9.982461e-01_r8,9.982206e-01_r8,9.9820e-01_r8 /)
      ssaliq1(:, 24) = (/ &
        & 9.999985e-01_r8,9.999979e-01_r8,9.999972e-01_r8,9.999966e-01_r8,9.999961e-01_r8,&
        & 9.999955e-01_r8,9.999950e-01_r8,9.999944e-01_r8,9.999938e-01_r8,9.999933e-01_r8,&
        & 9.999927e-01_r8,9.999921e-01_r8,9.999915e-01_r8,9.999910e-01_r8,9.999904e-01_r8,&
        & 9.999899e-01_r8,9.999893e-01_r8,9.999888e-01_r8,9.999882e-01_r8,9.999877e-01_r8,&
        & 9.999871e-01_r8,9.999866e-01_r8,9.999861e-01_r8,9.999855e-01_r8,9.999850e-01_r8,&
        & 9.999844e-01_r8,9.999839e-01_r8,9.999833e-01_r8,9.999828e-01_r8,9.999823e-01_r8,&
        & 9.999817e-01_r8,9.999812e-01_r8,9.999807e-01_r8,9.999801e-01_r8,9.999796e-01_r8,&
        & 9.999791e-01_r8,9.999786e-01_r8,9.999781e-01_r8,9.999776e-01_r8,9.999770e-01_r8,&
        & 9.999765e-01_r8,9.999761e-01_r8,9.999756e-01_r8,9.999751e-01_r8,9.999746e-01_r8,&
        & 9.999741e-01_r8,9.999736e-01_r8,9.999732e-01_r8,9.999727e-01_r8,9.999722e-01_r8,&
        & 9.999718e-01_r8,9.999713e-01_r8,9.999709e-01_r8,9.999705e-01_r8,9.999701e-01_r8,&
        & 9.999697e-01_r8,9.999692e-01_r8,9.9997e-01_r8 /)
      ssaliq1(:, 25) = (/ &
        & 9.999999e-01_r8,9.999998e-01_r8,9.999997e-01_r8,9.999997e-01_r8,9.999997e-01_r8,&
        & 9.999996e-01_r8,9.999996e-01_r8,9.999995e-01_r8,9.999995e-01_r8,9.999994e-01_r8,&
        & 9.999994e-01_r8,9.999993e-01_r8,9.999993e-01_r8,9.999992e-01_r8,9.999992e-01_r8,&
        & 9.999991e-01_r8,9.999991e-01_r8,9.999991e-01_r8,9.999990e-01_r8,9.999989e-01_r8,&
        & 9.999989e-01_r8,9.999989e-01_r8,9.999988e-01_r8,9.999988e-01_r8,9.999987e-01_r8,&
        & 9.999987e-01_r8,9.999986e-01_r8,9.999986e-01_r8,9.999985e-01_r8,9.999985e-01_r8,&
        & 9.999984e-01_r8,9.999984e-01_r8,9.999984e-01_r8,9.999983e-01_r8,9.999983e-01_r8,&
        & 9.999982e-01_r8,9.999982e-01_r8,9.999982e-01_r8,9.999981e-01_r8,9.999980e-01_r8,&
        & 9.999980e-01_r8,9.999980e-01_r8,9.999979e-01_r8,9.999979e-01_r8,9.999978e-01_r8,&
        & 9.999978e-01_r8,9.999977e-01_r8,9.999977e-01_r8,9.999977e-01_r8,9.999976e-01_r8,&
        & 9.999976e-01_r8,9.999975e-01_r8,9.999975e-01_r8,9.999974e-01_r8,9.999974e-01_r8,&
        & 9.999974e-01_r8,9.999973e-01_r8,1.0000e+00_r8 /)
      ssaliq1(:, 26) = (/ &
        & 9.999997e-01_r8,9.999995e-01_r8,9.999993e-01_r8,9.999992e-01_r8,9.999990e-01_r8,&
        & 9.999989e-01_r8,9.999988e-01_r8,9.999987e-01_r8,9.999986e-01_r8,9.999985e-01_r8,&
        & 9.999984e-01_r8,9.999983e-01_r8,9.999982e-01_r8,9.999981e-01_r8,9.999980e-01_r8,&
        & 9.999978e-01_r8,9.999977e-01_r8,9.999976e-01_r8,9.999975e-01_r8,9.999974e-01_r8,&
        & 9.999973e-01_r8,9.999972e-01_r8,9.999970e-01_r8,9.999969e-01_r8,9.999968e-01_r8,&
        & 9.999967e-01_r8,9.999966e-01_r8,9.999965e-01_r8,9.999964e-01_r8,9.999963e-01_r8,&
        & 9.999962e-01_r8,9.999961e-01_r8,9.999959e-01_r8,9.999958e-01_r8,9.999957e-01_r8,&
        & 9.999956e-01_r8,9.999955e-01_r8,9.999954e-01_r8,9.999953e-01_r8,9.999952e-01_r8,&
        & 9.999951e-01_r8,9.999949e-01_r8,9.999949e-01_r8,9.999947e-01_r8,9.999946e-01_r8,&
        & 9.999945e-01_r8,9.999944e-01_r8,9.999943e-01_r8,9.999942e-01_r8,9.999941e-01_r8,&
        & 9.999940e-01_r8,9.999939e-01_r8,9.999938e-01_r8,9.999937e-01_r8,9.999936e-01_r8,&
        & 9.999935e-01_r8,9.999934e-01_r8,9.9999e-01_r8 /)
      ssaliq1(:, 27) = (/ &
        & 9.999984e-01_r8,9.999976e-01_r8,9.999969e-01_r8,9.999962e-01_r8,9.999956e-01_r8,&
        & 9.999950e-01_r8,9.999945e-01_r8,9.999940e-01_r8,9.999935e-01_r8,9.999931e-01_r8,&
        & 9.999926e-01_r8,9.999920e-01_r8,9.999914e-01_r8,9.999908e-01_r8,9.999903e-01_r8,&
        & 9.999897e-01_r8,9.999891e-01_r8,9.999886e-01_r8,9.999880e-01_r8,9.999874e-01_r8,&
        & 9.999868e-01_r8,9.999863e-01_r8,9.999857e-01_r8,9.999851e-01_r8,9.999846e-01_r8,&
        & 9.999840e-01_r8,9.999835e-01_r8,9.999829e-01_r8,9.999824e-01_r8,9.999818e-01_r8,&
        & 9.999812e-01_r8,9.999806e-01_r8,9.999800e-01_r8,9.999795e-01_r8,9.999789e-01_r8,&
        & 9.999783e-01_r8,9.999778e-01_r8,9.999773e-01_r8,9.999767e-01_r8,9.999761e-01_r8,&
        & 9.999756e-01_r8,9.999750e-01_r8,9.999745e-01_r8,9.999739e-01_r8,9.999734e-01_r8,&
        & 9.999729e-01_r8,9.999723e-01_r8,9.999718e-01_r8,9.999713e-01_r8,9.999708e-01_r8,&
        & 9.999703e-01_r8,9.999697e-01_r8,9.999692e-01_r8,9.999687e-01_r8,9.999683e-01_r8,&
        & 9.999678e-01_r8,9.999673e-01_r8,9.9997e-01_r8 /)
      ssaliq1(:, 28) = (/ &
        & 9.999981e-01_r8,9.999973e-01_r8,9.999965e-01_r8,9.999958e-01_r8,9.999951e-01_r8,&
        & 9.999943e-01_r8,9.999937e-01_r8,9.999930e-01_r8,9.999924e-01_r8,9.999918e-01_r8,&
        & 9.999912e-01_r8,9.999905e-01_r8,9.999897e-01_r8,9.999890e-01_r8,9.999883e-01_r8,&
        & 9.999876e-01_r8,9.999869e-01_r8,9.999862e-01_r8,9.999855e-01_r8,9.999847e-01_r8,&
        & 9.999840e-01_r8,9.999834e-01_r8,9.999827e-01_r8,9.999819e-01_r8,9.999812e-01_r8,&
        & 9.999805e-01_r8,9.999799e-01_r8,9.999791e-01_r8,9.999785e-01_r8,9.999778e-01_r8,&
        & 9.999771e-01_r8,9.999764e-01_r8,9.999757e-01_r8,9.999750e-01_r8,9.999743e-01_r8,&
        & 9.999736e-01_r8,9.999729e-01_r8,9.999722e-01_r8,9.999715e-01_r8,9.999709e-01_r8,&
        & 9.999701e-01_r8,9.999695e-01_r8,9.999688e-01_r8,9.999682e-01_r8,9.999675e-01_r8,&
        & 9.999669e-01_r8,9.999662e-01_r8,9.999655e-01_r8,9.999649e-01_r8,9.999642e-01_r8,&
        & 9.999636e-01_r8,9.999630e-01_r8,9.999624e-01_r8,9.999618e-01_r8,9.999612e-01_r8,&
        & 9.999606e-01_r8,9.999600e-01_r8,9.9996e-01_r8 /)
      ssaliq1(:, 29) = (/ &
        & 8.505737e-01_r8,8.465102e-01_r8,8.394829e-01_r8,8.279508e-01_r8,8.110806e-01_r8,&
        & 7.900397e-01_r8,7.669615e-01_r8,7.444422e-01_r8,7.253055e-01_r8,7.124831e-01_r8,&
        & 7.016434e-01_r8,6.885485e-01_r8,6.767340e-01_r8,6.661029e-01_r8,6.565577e-01_r8,&
        & 6.480013e-01_r8,6.403373e-01_r8,6.334697e-01_r8,6.273034e-01_r8,6.217440e-01_r8,&
        & 6.166983e-01_r8,6.120740e-01_r8,6.077796e-01_r8,6.037249e-01_r8,5.998207e-01_r8,&
        & 5.959788e-01_r8,5.921123e-01_r8,5.881354e-01_r8,5.891285e-01_r8,5.851143e-01_r8,&
        & 5.814653e-01_r8,5.781606e-01_r8,5.751792e-01_r8,5.724998e-01_r8,5.701016e-01_r8,&
        & 5.679634e-01_r8,5.660642e-01_r8,5.643829e-01_r8,5.628984e-01_r8,5.615898e-01_r8,&
        & 5.604359e-01_r8,5.594158e-01_r8,5.585083e-01_r8,5.576924e-01_r8,5.569470e-01_r8,&
        & 5.562512e-01_r8,5.555838e-01_r8,5.549239e-01_r8,5.542503e-01_r8,5.535420e-01_r8,&
        & 5.527781e-01_r8,5.519374e-01_r8,5.509989e-01_r8,5.499417e-01_r8,5.487445e-01_r8,&
        & 5.473865e-01_r8,5.458466e-01_r8,5.4410e-01_r8 /)

! asymmetry parameter
      asyliq1(:, 16) = (/ &
        & 8.133297e-01_r8,8.133528e-01_r8,8.173865e-01_r8,8.243205e-01_r8,8.333063e-01_r8,&
        & 8.436317e-01_r8,8.546611e-01_r8,8.657934e-01_r8,8.764345e-01_r8,8.859837e-01_r8,&
        & 8.627394e-01_r8,8.824569e-01_r8,8.976887e-01_r8,9.089541e-01_r8,9.167699e-01_r8,&
        & 9.216517e-01_r8,9.241147e-01_r8,9.246743e-01_r8,9.238469e-01_r8,9.221504e-01_r8,&
        & 9.201045e-01_r8,9.182299e-01_r8,9.170491e-01_r8,9.170862e-01_r8,9.188653e-01_r8,&
        & 9.229111e-01_r8,9.297468e-01_r8,9.398950e-01_r8,9.203269e-01_r8,9.260693e-01_r8,&
        & 9.309373e-01_r8,9.349918e-01_r8,9.382935e-01_r8,9.409030e-01_r8,9.428809e-01_r8,&
        & 9.442881e-01_r8,9.451851e-01_r8,9.456331e-01_r8,9.456926e-01_r8,9.454247e-01_r8,&
        & 9.448902e-01_r8,9.441503e-01_r8,9.432661e-01_r8,9.422987e-01_r8,9.413094e-01_r8,&
        & 9.403594e-01_r8,9.395102e-01_r8,9.388230e-01_r8,9.383594e-01_r8,9.381810e-01_r8,&
        & 9.383489e-01_r8,9.389251e-01_r8,9.399707e-01_r8,9.415475e-01_r8,9.437167e-01_r8,&
        & 9.465399e-01_r8,9.500786e-01_r8,9.5439e-01_r8 /)
      asyliq1(:, 17) = (/ &
        & 8.794448e-01_r8,8.819306e-01_r8,8.837667e-01_r8,8.853832e-01_r8,8.871010e-01_r8,&
        & 8.892675e-01_r8,8.922584e-01_r8,8.964666e-01_r8,9.022940e-01_r8,9.101456e-01_r8,&
        & 8.839999e-01_r8,9.035610e-01_r8,9.184568e-01_r8,9.292315e-01_r8,9.364282e-01_r8,&
        & 9.405887e-01_r8,9.422554e-01_r8,9.419703e-01_r8,9.402759e-01_r8,9.377159e-01_r8,&
        & 9.348345e-01_r8,9.321769e-01_r8,9.302888e-01_r8,9.297166e-01_r8,9.310075e-01_r8,&
        & 9.347080e-01_r8,9.413643e-01_r8,9.515216e-01_r8,9.306286e-01_r8,9.361781e-01_r8,&
        & 9.408374e-01_r8,9.446692e-01_r8,9.477363e-01_r8,9.501013e-01_r8,9.518268e-01_r8,&
        & 9.529756e-01_r8,9.536105e-01_r8,9.537938e-01_r8,9.535886e-01_r8,9.530574e-01_r8,&
        & 9.522633e-01_r8,9.512688e-01_r8,9.501370e-01_r8,9.489306e-01_r8,9.477126e-01_r8,&
        & 9.465459e-01_r8,9.454934e-01_r8,9.446183e-01_r8,9.439833e-01_r8,9.436519e-01_r8,&
        & 9.436866e-01_r8,9.441508e-01_r8,9.451073e-01_r8,9.466195e-01_r8,9.487501e-01_r8,&
        & 9.515621e-01_r8,9.551185e-01_r8,9.5948e-01_r8 /)
      asyliq1(:, 18) = (/ &
        & 8.478817e-01_r8,8.269312e-01_r8,8.161352e-01_r8,8.135960e-01_r8,8.173586e-01_r8,&
        & 8.254167e-01_r8,8.357072e-01_r8,8.461167e-01_r8,8.544952e-01_r8,8.586776e-01_r8,&
        & 8.335562e-01_r8,8.524273e-01_r8,8.669052e-01_r8,8.775014e-01_r8,8.847277e-01_r8,&
        & 8.890958e-01_r8,8.911173e-01_r8,8.913038e-01_r8,8.901669e-01_r8,8.882182e-01_r8,&
        & 8.859692e-01_r8,8.839315e-01_r8,8.826164e-01_r8,8.825356e-01_r8,8.842004e-01_r8,&
        & 8.881223e-01_r8,8.948131e-01_r8,9.047837e-01_r8,8.855951e-01_r8,8.911796e-01_r8,&
        & 8.959229e-01_r8,8.998837e-01_r8,9.031209e-01_r8,9.056939e-01_r8,9.076609e-01_r8,&
        & 9.090812e-01_r8,9.100134e-01_r8,9.105167e-01_r8,9.106496e-01_r8,9.104712e-01_r8,&
        & 9.100404e-01_r8,9.094159e-01_r8,9.086568e-01_r8,9.078218e-01_r8,9.069697e-01_r8,&
        & 9.061595e-01_r8,9.054499e-01_r8,9.048999e-01_r8,9.045683e-01_r8,9.045142e-01_r8,&
        & 9.047962e-01_r8,9.054730e-01_r8,9.066037e-01_r8,9.082472e-01_r8,9.104623e-01_r8,&
        & 9.133079e-01_r8,9.168427e-01_r8,9.2113e-01_r8 /)
      asyliq1(:, 19) = (/ &
        & 8.216697e-01_r8,7.982871e-01_r8,7.891147e-01_r8,7.909083e-01_r8,8.003833e-01_r8,&
        & 8.142516e-01_r8,8.292290e-01_r8,8.420356e-01_r8,8.493945e-01_r8,8.480316e-01_r8,&
        & 8.212381e-01_r8,8.394984e-01_r8,8.534095e-01_r8,8.634813e-01_r8,8.702242e-01_r8,&
        & 8.741483e-01_r8,8.757638e-01_r8,8.755808e-01_r8,8.741095e-01_r8,8.718604e-01_r8,&
        & 8.693433e-01_r8,8.670686e-01_r8,8.655464e-01_r8,8.652872e-01_r8,8.668006e-01_r8,&
        & 8.705973e-01_r8,8.771874e-01_r8,8.870809e-01_r8,8.678284e-01_r8,8.732315e-01_r8,&
        & 8.778084e-01_r8,8.816166e-01_r8,8.847146e-01_r8,8.871603e-01_r8,8.890116e-01_r8,&
        & 8.903266e-01_r8,8.911632e-01_r8,8.915796e-01_r8,8.916337e-01_r8,8.913834e-01_r8,&
        & 8.908869e-01_r8,8.902022e-01_r8,8.893873e-01_r8,8.885001e-01_r8,8.875986e-01_r8,&
        & 8.867411e-01_r8,8.859852e-01_r8,8.853891e-01_r8,8.850111e-01_r8,8.849089e-01_r8,&
        & 8.851405e-01_r8,8.857639e-01_r8,8.868372e-01_r8,8.884185e-01_r8,8.905656e-01_r8,&
        & 8.933368e-01_r8,8.967899e-01_r8,9.0098e-01_r8 /)
      asyliq1(:, 20) = (/ &
        & 8.063610e-01_r8,7.938147e-01_r8,7.921304e-01_r8,7.985092e-01_r8,8.101339e-01_r8,&
        & 8.242175e-01_r8,8.379913e-01_r8,8.486920e-01_r8,8.535547e-01_r8,8.498083e-01_r8,&
        & 8.224849e-01_r8,8.405509e-01_r8,8.542436e-01_r8,8.640770e-01_r8,8.705653e-01_r8,&
        & 8.742227e-01_r8,8.755630e-01_r8,8.751004e-01_r8,8.733491e-01_r8,8.708231e-01_r8,&
        & 8.680365e-01_r8,8.655035e-01_r8,8.637381e-01_r8,8.632544e-01_r8,8.645665e-01_r8,&
        & 8.681885e-01_r8,8.746346e-01_r8,8.844188e-01_r8,8.648180e-01_r8,8.700563e-01_r8,&
        & 8.744672e-01_r8,8.781087e-01_r8,8.810393e-01_r8,8.833174e-01_r8,8.850011e-01_r8,&
        & 8.861485e-01_r8,8.868183e-01_r8,8.870687e-01_r8,8.869579e-01_r8,8.865441e-01_r8,&
        & 8.858857e-01_r8,8.850412e-01_r8,8.840686e-01_r8,8.830263e-01_r8,8.819726e-01_r8,&
        & 8.809658e-01_r8,8.800642e-01_r8,8.793260e-01_r8,8.788099e-01_r8,8.785737e-01_r8,&
        & 8.786758e-01_r8,8.791746e-01_r8,8.801283e-01_r8,8.815955e-01_r8,8.836340e-01_r8,&
        & 8.863024e-01_r8,8.896592e-01_r8,8.9376e-01_r8 /)
      asyliq1(:, 21) = (/ &
        & 7.885899e-01_r8,7.937172e-01_r8,8.020658e-01_r8,8.123971e-01_r8,8.235502e-01_r8,&
        & 8.343776e-01_r8,8.437336e-01_r8,8.504711e-01_r8,8.534421e-01_r8,8.514978e-01_r8,&
        & 8.238888e-01_r8,8.417463e-01_r8,8.552057e-01_r8,8.647853e-01_r8,8.710038e-01_r8,&
        & 8.743798e-01_r8,8.754319e-01_r8,8.746786e-01_r8,8.726386e-01_r8,8.698303e-01_r8,&
        & 8.667724e-01_r8,8.639836e-01_r8,8.619823e-01_r8,8.612870e-01_r8,8.624165e-01_r8,&
        & 8.658893e-01_r8,8.722241e-01_r8,8.819394e-01_r8,8.620216e-01_r8,8.671239e-01_r8,&
        & 8.713983e-01_r8,8.749032e-01_r8,8.776970e-01_r8,8.798385e-01_r8,8.813860e-01_r8,&
        & 8.823980e-01_r8,8.829332e-01_r8,8.830500e-01_r8,8.828068e-01_r8,8.822623e-01_r8,&
        & 8.814750e-01_r8,8.805031e-01_r8,8.794056e-01_r8,8.782407e-01_r8,8.770672e-01_r8,&
        & 8.759432e-01_r8,8.749275e-01_r8,8.740784e-01_r8,8.734547e-01_r8,8.731146e-01_r8,&
        & 8.731170e-01_r8,8.735199e-01_r8,8.743823e-01_r8,8.757625e-01_r8,8.777191e-01_r8,&
        & 8.803105e-01_r8,8.835953e-01_r8,8.8763e-01_r8 /)
      asyliq1(:, 22) = (/ &
        & 7.811516e-01_r8,7.962229e-01_r8,8.096199e-01_r8,8.212996e-01_r8,8.312212e-01_r8,&
        & 8.393430e-01_r8,8.456236e-01_r8,8.500214e-01_r8,8.524950e-01_r8,8.530031e-01_r8,&
        & 8.251485e-01_r8,8.429043e-01_r8,8.562461e-01_r8,8.656954e-01_r8,8.717737e-01_r8,&
        & 8.750020e-01_r8,8.759022e-01_r8,8.749953e-01_r8,8.728027e-01_r8,8.698461e-01_r8,&
        & 8.666466e-01_r8,8.637257e-01_r8,8.616047e-01_r8,8.608051e-01_r8,8.618483e-01_r8,&
        & 8.652557e-01_r8,8.715487e-01_r8,8.812485e-01_r8,8.611645e-01_r8,8.662052e-01_r8,&
        & 8.704173e-01_r8,8.738594e-01_r8,8.765901e-01_r8,8.786678e-01_r8,8.801517e-01_r8,&
        & 8.810999e-01_r8,8.815713e-01_r8,8.816246e-01_r8,8.813185e-01_r8,8.807114e-01_r8,&
        & 8.798621e-01_r8,8.788290e-01_r8,8.776713e-01_r8,8.764470e-01_r8,8.752152e-01_r8,&
        & 8.740343e-01_r8,8.729631e-01_r8,8.720602e-01_r8,8.713842e-01_r8,8.709936e-01_r8,&
        & 8.709475e-01_r8,8.713041e-01_r8,8.721221e-01_r8,8.734602e-01_r8,8.753774e-01_r8,&
        & 8.779319e-01_r8,8.811825e-01_r8,8.8519e-01_r8 /)
      asyliq1(:, 23) = (/ &
        & 7.865744e-01_r8,8.093340e-01_r8,8.257596e-01_r8,8.369940e-01_r8,8.441574e-01_r8,&
        & 8.483602e-01_r8,8.507096e-01_r8,8.523139e-01_r8,8.542834e-01_r8,8.577321e-01_r8,&
        & 8.288960e-01_r8,8.465308e-01_r8,8.597175e-01_r8,8.689830e-01_r8,8.748542e-01_r8,&
        & 8.778584e-01_r8,8.785222e-01_r8,8.773728e-01_r8,8.749370e-01_r8,8.717419e-01_r8,&
        & 8.683145e-01_r8,8.651816e-01_r8,8.628704e-01_r8,8.619077e-01_r8,8.628205e-01_r8,&
        & 8.661356e-01_r8,8.723803e-01_r8,8.820815e-01_r8,8.616715e-01_r8,8.666389e-01_r8,&
        & 8.707753e-01_r8,8.741398e-01_r8,8.767912e-01_r8,8.787885e-01_r8,8.801908e-01_r8,&
        & 8.810570e-01_r8,8.814460e-01_r8,8.814167e-01_r8,8.810283e-01_r8,8.803395e-01_r8,&
        & 8.794095e-01_r8,8.782971e-01_r8,8.770613e-01_r8,8.757610e-01_r8,8.744553e-01_r8,&
        & 8.732031e-01_r8,8.720634e-01_r8,8.710951e-01_r8,8.703572e-01_r8,8.699086e-01_r8,&
        & 8.698084e-01_r8,8.701155e-01_r8,8.708887e-01_r8,8.721872e-01_r8,8.740698e-01_r8,&
        & 8.765957e-01_r8,8.798235e-01_r8,8.8381e-01_r8 /)
      asyliq1(:, 24) = (/ &
        & 8.069513e-01_r8,8.262939e-01_r8,8.398241e-01_r8,8.486352e-01_r8,8.538213e-01_r8,&
        & 8.564743e-01_r8,8.576854e-01_r8,8.585455e-01_r8,8.601452e-01_r8,8.635755e-01_r8,&
        & 8.337383e-01_r8,8.512655e-01_r8,8.643049e-01_r8,8.733896e-01_r8,8.790535e-01_r8,&
        & 8.818295e-01_r8,8.822518e-01_r8,8.808533e-01_r8,8.781676e-01_r8,8.747284e-01_r8,&
        & 8.710690e-01_r8,8.677229e-01_r8,8.652236e-01_r8,8.641047e-01_r8,8.648993e-01_r8,&
        & 8.681413e-01_r8,8.743640e-01_r8,8.841007e-01_r8,8.633558e-01_r8,8.682719e-01_r8,&
        & 8.723543e-01_r8,8.756621e-01_r8,8.782547e-01_r8,8.801915e-01_r8,8.815318e-01_r8,&
        & 8.823347e-01_r8,8.826598e-01_r8,8.825663e-01_r8,8.821135e-01_r8,8.813608e-01_r8,&
        & 8.803674e-01_r8,8.791928e-01_r8,8.778960e-01_r8,8.765366e-01_r8,8.751738e-01_r8,&
        & 8.738670e-01_r8,8.726755e-01_r8,8.716585e-01_r8,8.708755e-01_r8,8.703856e-01_r8,&
        & 8.702483e-01_r8,8.705229e-01_r8,8.712687e-01_r8,8.725448e-01_r8,8.744109e-01_r8,&
        & 8.769260e-01_r8,8.801496e-01_r8,8.8414e-01_r8 /)
      asyliq1(:, 25) = (/ &
        & 8.252182e-01_r8,8.379244e-01_r8,8.471709e-01_r8,8.535760e-01_r8,8.577540e-01_r8,&
        & 8.603183e-01_r8,8.618820e-01_r8,8.630578e-01_r8,8.644587e-01_r8,8.666970e-01_r8,&
        & 8.362159e-01_r8,8.536817e-01_r8,8.666387e-01_r8,8.756240e-01_r8,8.811746e-01_r8,&
        & 8.838273e-01_r8,8.841191e-01_r8,8.825871e-01_r8,8.797681e-01_r8,8.761992e-01_r8,&
        & 8.724174e-01_r8,8.689593e-01_r8,8.663623e-01_r8,8.651632e-01_r8,8.658988e-01_r8,&
        & 8.691064e-01_r8,8.753226e-01_r8,8.850847e-01_r8,8.641620e-01_r8,8.690500e-01_r8,&
        & 8.731026e-01_r8,8.763795e-01_r8,8.789400e-01_r8,8.808438e-01_r8,8.821503e-01_r8,&
        & 8.829191e-01_r8,8.832095e-01_r8,8.830813e-01_r8,8.825938e-01_r8,8.818064e-01_r8,&
        & 8.807787e-01_r8,8.795704e-01_r8,8.782408e-01_r8,8.768493e-01_r8,8.754557e-01_r8,&
        & 8.741193e-01_r8,8.728995e-01_r8,8.718561e-01_r8,8.710484e-01_r8,8.705360e-01_r8,&
        & 8.703782e-01_r8,8.706347e-01_r8,8.713650e-01_r8,8.726285e-01_r8,8.744849e-01_r8,&
        & 8.769933e-01_r8,8.802136e-01_r8,8.8421e-01_r8 /)
      asyliq1(:, 26) = (/ &
        & 8.370583e-01_r8,8.467920e-01_r8,8.537769e-01_r8,8.585136e-01_r8,8.615034e-01_r8,&
        & 8.632474e-01_r8,8.642468e-01_r8,8.650026e-01_r8,8.660161e-01_r8,8.677882e-01_r8,&
        & 8.369760e-01_r8,8.543821e-01_r8,8.672699e-01_r8,8.761782e-01_r8,8.816454e-01_r8,&
        & 8.842103e-01_r8,8.844114e-01_r8,8.827872e-01_r8,8.798766e-01_r8,8.762179e-01_r8,&
        & 8.723500e-01_r8,8.688112e-01_r8,8.661403e-01_r8,8.648758e-01_r8,8.655563e-01_r8,&
        & 8.687206e-01_r8,8.749072e-01_r8,8.846546e-01_r8,8.636289e-01_r8,8.684849e-01_r8,&
        & 8.725054e-01_r8,8.757501e-01_r8,8.782785e-01_r8,8.801503e-01_r8,8.814249e-01_r8,&
        & 8.821620e-01_r8,8.824211e-01_r8,8.822620e-01_r8,8.817440e-01_r8,8.809268e-01_r8,&
        & 8.798699e-01_r8,8.786330e-01_r8,8.772756e-01_r8,8.758572e-01_r8,8.744374e-01_r8,&
        & 8.730760e-01_r8,8.718323e-01_r8,8.707660e-01_r8,8.699366e-01_r8,8.694039e-01_r8,&
        & 8.692271e-01_r8,8.694661e-01_r8,8.701803e-01_r8,8.714293e-01_r8,8.732727e-01_r8,&
        & 8.757702e-01_r8,8.789811e-01_r8,8.8297e-01_r8 /)
      asyliq1(:, 27) = (/ &
        & 8.430819e-01_r8,8.510060e-01_r8,8.567270e-01_r8,8.606533e-01_r8,8.631934e-01_r8,&
        & 8.647554e-01_r8,8.657471e-01_r8,8.665760e-01_r8,8.676496e-01_r8,8.693754e-01_r8,&
        & 8.384298e-01_r8,8.557913e-01_r8,8.686214e-01_r8,8.774605e-01_r8,8.828495e-01_r8,&
        & 8.853287e-01_r8,8.854393e-01_r8,8.837215e-01_r8,8.807161e-01_r8,8.769639e-01_r8,&
        & 8.730053e-01_r8,8.693812e-01_r8,8.666321e-01_r8,8.652988e-01_r8,8.659219e-01_r8,&
        & 8.690419e-01_r8,8.751999e-01_r8,8.849360e-01_r8,8.638013e-01_r8,8.686371e-01_r8,&
        & 8.726369e-01_r8,8.758605e-01_r8,8.783674e-01_r8,8.802176e-01_r8,8.814705e-01_r8,&
        & 8.821859e-01_r8,8.824234e-01_r8,8.822429e-01_r8,8.817038e-01_r8,8.808658e-01_r8,&
        & 8.797887e-01_r8,8.785323e-01_r8,8.771560e-01_r8,8.757196e-01_r8,8.742828e-01_r8,&
        & 8.729052e-01_r8,8.716467e-01_r8,8.705666e-01_r8,8.697250e-01_r8,8.691812e-01_r8,&
        & 8.689950e-01_r8,8.692264e-01_r8,8.699346e-01_r8,8.711795e-01_r8,8.730209e-01_r8,&
        & 8.755181e-01_r8,8.787312e-01_r8,8.8272e-01_r8 /)
      asyliq1(:, 28) = (/ &
        & 8.452284e-01_r8,8.522700e-01_r8,8.572973e-01_r8,8.607031e-01_r8,8.628802e-01_r8,&
        & 8.642215e-01_r8,8.651198e-01_r8,8.659679e-01_r8,8.671588e-01_r8,8.690853e-01_r8,&
        & 8.383803e-01_r8,8.557485e-01_r8,8.685851e-01_r8,8.774303e-01_r8,8.828245e-01_r8,&
        & 8.853077e-01_r8,8.854207e-01_r8,8.837034e-01_r8,8.806962e-01_r8,8.769398e-01_r8,&
        & 8.729740e-01_r8,8.693393e-01_r8,8.665761e-01_r8,8.652247e-01_r8,8.658253e-01_r8,&
        & 8.689182e-01_r8,8.750438e-01_r8,8.847424e-01_r8,8.636140e-01_r8,8.684449e-01_r8,&
        & 8.724400e-01_r8,8.756589e-01_r8,8.781613e-01_r8,8.800072e-01_r8,8.812559e-01_r8,&
        & 8.819671e-01_r8,8.822007e-01_r8,8.820165e-01_r8,8.814737e-01_r8,8.806322e-01_r8,&
        & 8.795518e-01_r8,8.782923e-01_r8,8.769129e-01_r8,8.754737e-01_r8,8.740342e-01_r8,&
        & 8.726542e-01_r8,8.713934e-01_r8,8.703111e-01_r8,8.694677e-01_r8,8.689222e-01_r8,&
        & 8.687344e-01_r8,8.689646e-01_r8,8.696715e-01_r8,8.709156e-01_r8,8.727563e-01_r8,&
        & 8.752531e-01_r8,8.784659e-01_r8,8.8245e-01_r8 /)
      asyliq1(:, 29) = (/ &
        & 7.800869e-01_r8,8.091120e-01_r8,8.325369e-01_r8,8.466266e-01_r8,8.515495e-01_r8,&
        & 8.499371e-01_r8,8.456203e-01_r8,8.430521e-01_r8,8.470286e-01_r8,8.625431e-01_r8,&
        & 8.402261e-01_r8,8.610822e-01_r8,8.776608e-01_r8,8.904485e-01_r8,8.999294e-01_r8,&
        & 9.065860e-01_r8,9.108995e-01_r8,9.133503e-01_r8,9.144187e-01_r8,9.145855e-01_r8,&
        & 9.143320e-01_r8,9.141402e-01_r8,9.144933e-01_r8,9.158754e-01_r8,9.187716e-01_r8,&
        & 9.236677e-01_r8,9.310503e-01_r8,9.414058e-01_r8,9.239108e-01_r8,9.300719e-01_r8,&
        & 9.353612e-01_r8,9.398378e-01_r8,9.435609e-01_r8,9.465895e-01_r8,9.489829e-01_r8,&
        & 9.508000e-01_r8,9.521002e-01_r8,9.529424e-01_r8,9.533860e-01_r8,9.534902e-01_r8,&
        & 9.533143e-01_r8,9.529177e-01_r8,9.523596e-01_r8,9.516997e-01_r8,9.509973e-01_r8,&
        & 9.503121e-01_r8,9.497037e-01_r8,9.492317e-01_r8,9.489558e-01_r8,9.489356e-01_r8,&
        & 9.492311e-01_r8,9.499019e-01_r8,9.510077e-01_r8,9.526084e-01_r8,9.547636e-01_r8,&
        & 9.575331e-01_r8,9.609766e-01_r8,9.6515e-01_r8 /)

! Spherical Ice Particle Parameterization
! extinction units (ext coef/iwc): [(m^-1)/(g m^-3)]
      extice2(:, 16) = (/ &
! band 16
        & 4.101824e-01_r8,2.435514e-01_r8,1.713697e-01_r8,1.314865e-01_r8,1.063406e-01_r8,&
        & 8.910701e-02_r8,7.659480e-02_r8,6.711784e-02_r8,5.970353e-02_r8,5.375249e-02_r8,&
        & 4.887577e-02_r8,4.481025e-02_r8,4.137171e-02_r8,3.842744e-02_r8,3.587948e-02_r8,&
        & 3.365396e-02_r8,3.169419e-02_r8,2.995593e-02_r8,2.840419e-02_r8,2.701091e-02_r8,&
        & 2.575336e-02_r8,2.461293e-02_r8,2.357423e-02_r8,2.262443e-02_r8,2.175276e-02_r8,&
        & 2.095012e-02_r8,2.020875e-02_r8,1.952199e-02_r8,1.888412e-02_r8,1.829018e-02_r8,&
        & 1.773586e-02_r8,1.721738e-02_r8,1.673144e-02_r8,1.627510e-02_r8,1.584579e-02_r8,&
        & 1.544122e-02_r8,1.505934e-02_r8,1.469833e-02_r8,1.435654e-02_r8,1.403251e-02_r8,&
        & 1.372492e-02_r8,1.343255e-02_r8,1.315433e-02_r8 /)
      extice2(:, 17) = (/ &
! band 17
        & 3.836650e-01_r8,2.304055e-01_r8,1.637265e-01_r8,1.266681e-01_r8,1.031602e-01_r8,&
        & 8.695191e-02_r8,7.511544e-02_r8,6.610009e-02_r8,5.900909e-02_r8,5.328833e-02_r8,&
        & 4.857728e-02_r8,4.463133e-02_r8,4.127880e-02_r8,3.839567e-02_r8,3.589013e-02_r8,&
        & 3.369280e-02_r8,3.175027e-02_r8,3.002079e-02_r8,2.847121e-02_r8,2.707493e-02_r8,&
        & 2.581031e-02_r8,2.465962e-02_r8,2.360815e-02_r8,2.264363e-02_r8,2.175571e-02_r8,&
        & 2.093563e-02_r8,2.017592e-02_r8,1.947015e-02_r8,1.881278e-02_r8,1.819901e-02_r8,&
        & 1.762463e-02_r8,1.708598e-02_r8,1.657982e-02_r8,1.610330e-02_r8,1.565390e-02_r8,&
        & 1.522937e-02_r8,1.482768e-02_r8,1.444706e-02_r8,1.408588e-02_r8,1.374270e-02_r8,&
        & 1.341619e-02_r8,1.310517e-02_r8,1.280857e-02_r8 /)
      extice2(:, 18) = (/ &
! band 18
        & 4.152673e-01_r8,2.436816e-01_r8,1.702243e-01_r8,1.299704e-01_r8,1.047528e-01_r8,&
        & 8.756039e-02_r8,7.513327e-02_r8,6.575690e-02_r8,5.844616e-02_r8,5.259609e-02_r8,&
        & 4.781531e-02_r8,4.383980e-02_r8,4.048517e-02_r8,3.761891e-02_r8,3.514342e-02_r8,&
        & 3.298525e-02_r8,3.108814e-02_r8,2.940825e-02_r8,2.791096e-02_r8,2.656858e-02_r8,&
        & 2.535869e-02_r8,2.426297e-02_r8,2.326627e-02_r8,2.235602e-02_r8,2.152164e-02_r8,&
        & 2.075420e-02_r8,2.004613e-02_r8,1.939091e-02_r8,1.878296e-02_r8,1.821744e-02_r8,&
        & 1.769015e-02_r8,1.719741e-02_r8,1.673600e-02_r8,1.630308e-02_r8,1.589615e-02_r8,&
        & 1.551298e-02_r8,1.515159e-02_r8,1.481021e-02_r8,1.448726e-02_r8,1.418131e-02_r8,&
        & 1.389109e-02_r8,1.361544e-02_r8,1.335330e-02_r8 /)
      extice2(:, 19) = (/ &
! band 19
        & 3.873250e-01_r8,2.331609e-01_r8,1.655002e-01_r8,1.277753e-01_r8,1.038247e-01_r8,&
        & 8.731780e-02_r8,7.527638e-02_r8,6.611873e-02_r8,5.892850e-02_r8,5.313885e-02_r8,&
        & 4.838068e-02_r8,4.440356e-02_r8,4.103167e-02_r8,3.813804e-02_r8,3.562870e-02_r8,&
        & 3.343269e-02_r8,3.149539e-02_r8,2.977414e-02_r8,2.823510e-02_r8,2.685112e-02_r8,&
        & 2.560015e-02_r8,2.446411e-02_r8,2.342805e-02_r8,2.247948e-02_r8,2.160789e-02_r8,&
        & 2.080438e-02_r8,2.006139e-02_r8,1.937238e-02_r8,1.873177e-02_r8,1.813469e-02_r8,&
        & 1.757689e-02_r8,1.705468e-02_r8,1.656479e-02_r8,1.610435e-02_r8,1.567081e-02_r8,&
        & 1.526192e-02_r8,1.487565e-02_r8,1.451020e-02_r8,1.416396e-02_r8,1.383546e-02_r8,&
        & 1.352339e-02_r8,1.322657e-02_r8,1.294392e-02_r8 /)
      extice2(:, 20) = (/ &
! band 20
        & 3.784280e-01_r8,2.291396e-01_r8,1.632551e-01_r8,1.263775e-01_r8,1.028944e-01_r8,&
        & 8.666975e-02_r8,7.480952e-02_r8,6.577335e-02_r8,5.866714e-02_r8,5.293694e-02_r8,&
        & 4.822153e-02_r8,4.427547e-02_r8,4.092626e-02_r8,3.804918e-02_r8,3.555184e-02_r8,&
        & 3.336440e-02_r8,3.143307e-02_r8,2.971577e-02_r8,2.817912e-02_r8,2.679632e-02_r8,&
        & 2.554558e-02_r8,2.440903e-02_r8,2.337187e-02_r8,2.242173e-02_r8,2.154821e-02_r8,&
        & 2.074249e-02_r8,1.999706e-02_r8,1.930546e-02_r8,1.866212e-02_r8,1.806221e-02_r8,&
        & 1.750152e-02_r8,1.697637e-02_r8,1.648352e-02_r8,1.602010e-02_r8,1.558358e-02_r8,&
        & 1.517172e-02_r8,1.478250e-02_r8,1.441413e-02_r8,1.406498e-02_r8,1.373362e-02_r8,&
        & 1.341872e-02_r8,1.311911e-02_r8,1.283371e-02_r8 /)
      extice2(:, 21) = (/ &
! band 21
        & 3.719909e-01_r8,2.259490e-01_r8,1.613144e-01_r8,1.250648e-01_r8,1.019462e-01_r8,&
        & 8.595358e-02_r8,7.425064e-02_r8,6.532618e-02_r8,5.830218e-02_r8,5.263421e-02_r8,&
        & 4.796697e-02_r8,4.405891e-02_r8,4.074013e-02_r8,3.788776e-02_r8,3.541071e-02_r8,&
        & 3.324008e-02_r8,3.132280e-02_r8,2.961733e-02_r8,2.809071e-02_r8,2.671645e-02_r8,&
        & 2.547302e-02_r8,2.434276e-02_r8,2.331102e-02_r8,2.236558e-02_r8,2.149614e-02_r8,&
        & 2.069397e-02_r8,1.995163e-02_r8,1.926272e-02_r8,1.862174e-02_r8,1.802389e-02_r8,&
        & 1.746500e-02_r8,1.694142e-02_r8,1.644994e-02_r8,1.598772e-02_r8,1.555225e-02_r8,&
        & 1.514129e-02_r8,1.475286e-02_r8,1.438515e-02_r8,1.403659e-02_r8,1.370572e-02_r8,&
        & 1.339124e-02_r8,1.309197e-02_r8,1.280685e-02_r8 /)
      extice2(:, 22) = (/ &
! band 22
        & 3.713158e-01_r8,2.253816e-01_r8,1.608461e-01_r8,1.246718e-01_r8,1.016109e-01_r8,&
        & 8.566332e-02_r8,7.399666e-02_r8,6.510199e-02_r8,5.810290e-02_r8,5.245608e-02_r8,&
        & 4.780702e-02_r8,4.391478e-02_r8,4.060989e-02_r8,3.776982e-02_r8,3.530374e-02_r8,&
        & 3.314296e-02_r8,3.123458e-02_r8,2.953719e-02_r8,2.801794e-02_r8,2.665043e-02_r8,&
        & 2.541321e-02_r8,2.428868e-02_r8,2.326224e-02_r8,2.232173e-02_r8,2.145688e-02_r8,&
        & 2.065899e-02_r8,1.992067e-02_r8,1.923552e-02_r8,1.859808e-02_r8,1.800356e-02_r8,&
        & 1.744782e-02_r8,1.692721e-02_r8,1.643855e-02_r8,1.597900e-02_r8,1.554606e-02_r8,&
        & 1.513751e-02_r8,1.475137e-02_r8,1.438586e-02_r8,1.403938e-02_r8,1.371050e-02_r8,&
        & 1.339793e-02_r8,1.310050e-02_r8,1.281713e-02_r8 /)
      extice2(:, 23) = (/ &
! band 23
        & 3.605883e-01_r8,2.204388e-01_r8,1.580431e-01_r8,1.229033e-01_r8,1.004203e-01_r8,&
        & 8.482616e-02_r8,7.338941e-02_r8,6.465105e-02_r8,5.776176e-02_r8,5.219398e-02_r8,&
        & 4.760288e-02_r8,4.375369e-02_r8,4.048111e-02_r8,3.766539e-02_r8,3.521771e-02_r8,&
        & 3.307079e-02_r8,3.117277e-02_r8,2.948303e-02_r8,2.796929e-02_r8,2.660560e-02_r8,&
        & 2.537086e-02_r8,2.424772e-02_r8,2.322182e-02_r8,2.228114e-02_r8,2.141556e-02_r8,&
        & 2.061649e-02_r8,1.987661e-02_r8,1.918962e-02_r8,1.855009e-02_r8,1.795330e-02_r8,&
        & 1.739514e-02_r8,1.687199e-02_r8,1.638069e-02_r8,1.591845e-02_r8,1.548276e-02_r8,&
        & 1.507143e-02_r8,1.468249e-02_r8,1.431416e-02_r8,1.396486e-02_r8,1.363318e-02_r8,&
        & 1.331781e-02_r8,1.301759e-02_r8,1.273147e-02_r8 /)
      extice2(:, 24) = (/ &
! band 24
        & 3.527890e-01_r8,2.168469e-01_r8,1.560090e-01_r8,1.216216e-01_r8,9.955787e-02_r8,&
        & 8.421942e-02_r8,7.294827e-02_r8,6.432192e-02_r8,5.751081e-02_r8,5.199888e-02_r8,&
        & 4.744835e-02_r8,4.362899e-02_r8,4.037847e-02_r8,3.757910e-02_r8,3.514351e-02_r8,&
        & 3.300546e-02_r8,3.111382e-02_r8,2.942853e-02_r8,2.791775e-02_r8,2.655584e-02_r8,&
        & 2.532195e-02_r8,2.419892e-02_r8,2.317255e-02_r8,2.223092e-02_r8,2.136402e-02_r8,&
        & 2.056334e-02_r8,1.982160e-02_r8,1.913258e-02_r8,1.849087e-02_r8,1.789178e-02_r8,&
        & 1.733124e-02_r8,1.680565e-02_r8,1.631187e-02_r8,1.584711e-02_r8,1.540889e-02_r8,&
        & 1.499502e-02_r8,1.460354e-02_r8,1.423269e-02_r8,1.388088e-02_r8,1.354670e-02_r8,&
        & 1.322887e-02_r8,1.292620e-02_r8,1.263767e-02_r8 /)
      extice2(:, 25) = (/ &
! band 25
        & 3.477874e-01_r8,2.143515e-01_r8,1.544887e-01_r8,1.205942e-01_r8,9.881779e-02_r8,&
        & 8.366261e-02_r8,7.251586e-02_r8,6.397790e-02_r8,5.723183e-02_r8,5.176908e-02_r8,&
        & 4.725658e-02_r8,4.346715e-02_r8,4.024055e-02_r8,3.746055e-02_r8,3.504080e-02_r8,&
        & 3.291583e-02_r8,3.103507e-02_r8,2.935891e-02_r8,2.785582e-02_r8,2.650042e-02_r8,&
        & 2.527206e-02_r8,2.415376e-02_r8,2.313142e-02_r8,2.219326e-02_r8,2.132934e-02_r8,&
        & 2.053122e-02_r8,1.979169e-02_r8,1.910456e-02_r8,1.846448e-02_r8,1.786680e-02_r8,&
        & 1.730745e-02_r8,1.678289e-02_r8,1.628998e-02_r8,1.582595e-02_r8,1.538835e-02_r8,&
        & 1.497499e-02_r8,1.458393e-02_r8,1.421341e-02_r8,1.386187e-02_r8,1.352788e-02_r8,&
        & 1.321019e-02_r8,1.290762e-02_r8,1.261913e-02_r8 /)
      extice2(:, 26) = (/ &
! band 26
        & 3.453721e-01_r8,2.130744e-01_r8,1.536698e-01_r8,1.200140e-01_r8,9.838078e-02_r8,&
        & 8.331940e-02_r8,7.223803e-02_r8,6.374775e-02_r8,5.703770e-02_r8,5.160290e-02_r8,&
        & 4.711259e-02_r8,4.334110e-02_r8,4.012923e-02_r8,3.736150e-02_r8,3.495208e-02_r8,&
        & 3.283589e-02_r8,3.096267e-02_r8,2.929302e-02_r8,2.779560e-02_r8,2.644517e-02_r8,&
        & 2.522119e-02_r8,2.410677e-02_r8,2.308788e-02_r8,2.215281e-02_r8,2.129165e-02_r8,&
        & 2.049602e-02_r8,1.975874e-02_r8,1.907365e-02_r8,1.843542e-02_r8,1.783943e-02_r8,&
        & 1.728162e-02_r8,1.675847e-02_r8,1.626685e-02_r8,1.580401e-02_r8,1.536750e-02_r8,&
        & 1.495515e-02_r8,1.456502e-02_r8,1.419537e-02_r8,1.384463e-02_r8,1.351139e-02_r8,&
        & 1.319438e-02_r8,1.289246e-02_r8,1.260456e-02_r8 /)
      extice2(:, 27) = (/ &
! band 27
        & 3.417883e-01_r8,2.113379e-01_r8,1.526395e-01_r8,1.193347e-01_r8,9.790253e-02_r8,&
        & 8.296715e-02_r8,7.196979e-02_r8,6.353806e-02_r8,5.687024e-02_r8,5.146670e-02_r8,&
        & 4.700001e-02_r8,4.324667e-02_r8,4.004894e-02_r8,3.729233e-02_r8,3.489172e-02_r8,&
        & 3.278257e-02_r8,3.091499e-02_r8,2.924987e-02_r8,2.775609e-02_r8,2.640859e-02_r8,&
        & 2.518695e-02_r8,2.407439e-02_r8,2.305697e-02_r8,2.212303e-02_r8,2.126273e-02_r8,&
        & 2.046774e-02_r8,1.973090e-02_r8,1.904610e-02_r8,1.840801e-02_r8,1.781204e-02_r8,&
        & 1.725417e-02_r8,1.673086e-02_r8,1.623902e-02_r8,1.577590e-02_r8,1.533906e-02_r8,&
        & 1.492634e-02_r8,1.453580e-02_r8,1.416571e-02_r8,1.381450e-02_r8,1.348078e-02_r8,&
        & 1.316327e-02_r8,1.286082e-02_r8,1.257240e-02_r8 /)
      extice2(:, 28) = (/ &
! band 28
        & 3.416111e-01_r8,2.114124e-01_r8,1.527734e-01_r8,1.194809e-01_r8,9.804612e-02_r8,&
        & 8.310287e-02_r8,7.209595e-02_r8,6.365442e-02_r8,5.697710e-02_r8,5.156460e-02_r8,&
        & 4.708957e-02_r8,4.332850e-02_r8,4.012361e-02_r8,3.736037e-02_r8,3.495364e-02_r8,&
        & 3.283879e-02_r8,3.096593e-02_r8,2.929589e-02_r8,2.779751e-02_r8,2.644571e-02_r8,&
        & 2.522004e-02_r8,2.410369e-02_r8,2.308271e-02_r8,2.214542e-02_r8,2.128195e-02_r8,&
        & 2.048396e-02_r8,1.974429e-02_r8,1.905679e-02_r8,1.841614e-02_r8,1.781774e-02_r8,&
        & 1.725754e-02_r8,1.673203e-02_r8,1.623807e-02_r8,1.577293e-02_r8,1.533416e-02_r8,&
        & 1.491958e-02_r8,1.452727e-02_r8,1.415547e-02_r8,1.380262e-02_r8,1.346732e-02_r8,&
        & 1.314830e-02_r8,1.284439e-02_r8,1.255456e-02_r8 /)
      extice2(:, 29) = (/ &
! band 29
        & 4.196611e-01_r8,2.493642e-01_r8,1.761261e-01_r8,1.357197e-01_r8,1.102161e-01_r8,&
        & 9.269376e-02_r8,7.992985e-02_r8,7.022538e-02_r8,6.260168e-02_r8,5.645603e-02_r8,&
        & 5.139732e-02_r8,4.716088e-02_r8,4.356133e-02_r8,4.046498e-02_r8,3.777303e-02_r8,&
        & 3.541094e-02_r8,3.332137e-02_r8,3.145954e-02_r8,2.978998e-02_r8,2.828419e-02_r8,&
        & 2.691905e-02_r8,2.567559e-02_r8,2.453811e-02_r8,2.349350e-02_r8,2.253072e-02_r8,&
        & 2.164042e-02_r8,2.081464e-02_r8,2.004652e-02_r8,1.933015e-02_r8,1.866041e-02_r8,&
        & 1.803283e-02_r8,1.744348e-02_r8,1.688894e-02_r8,1.636616e-02_r8,1.587244e-02_r8,&
        & 1.540539e-02_r8,1.496287e-02_r8,1.454295e-02_r8,1.414392e-02_r8,1.376423e-02_r8,&
        & 1.340247e-02_r8,1.305739e-02_r8,1.272784e-02_r8 /)

! single-scattering albedo: unitless
      ssaice2(:, 16) = (/ &
! band 16
        & 6.630615e-01_r8,6.451169e-01_r8,6.333696e-01_r8,6.246927e-01_r8,6.178420e-01_r8,&
        & 6.121976e-01_r8,6.074069e-01_r8,6.032505e-01_r8,5.995830e-01_r8,5.963030e-01_r8,&
        & 5.933372e-01_r8,5.906311e-01_r8,5.881427e-01_r8,5.858395e-01_r8,5.836955e-01_r8,&
        & 5.816896e-01_r8,5.798046e-01_r8,5.780264e-01_r8,5.763429e-01_r8,5.747441e-01_r8,&
        & 5.732213e-01_r8,5.717672e-01_r8,5.703754e-01_r8,5.690403e-01_r8,5.677571e-01_r8,&
        & 5.665215e-01_r8,5.653297e-01_r8,5.641782e-01_r8,5.630643e-01_r8,5.619850e-01_r8,&
        & 5.609381e-01_r8,5.599214e-01_r8,5.589328e-01_r8,5.579707e-01_r8,5.570333e-01_r8,&
        & 5.561193e-01_r8,5.552272e-01_r8,5.543558e-01_r8,5.535041e-01_r8,5.526708e-01_r8,&
        & 5.518551e-01_r8,5.510561e-01_r8,5.502729e-01_r8 /)
      ssaice2(:, 17) = (/ &
! band 17
        & 7.689749e-01_r8,7.398171e-01_r8,7.205819e-01_r8,7.065690e-01_r8,6.956928e-01_r8,&
        & 6.868989e-01_r8,6.795813e-01_r8,6.733606e-01_r8,6.679838e-01_r8,6.632742e-01_r8,&
        & 6.591036e-01_r8,6.553766e-01_r8,6.520197e-01_r8,6.489757e-01_r8,6.461991e-01_r8,&
        & 6.436531e-01_r8,6.413075e-01_r8,6.391375e-01_r8,6.371221e-01_r8,6.352438e-01_r8,&
        & 6.334876e-01_r8,6.318406e-01_r8,6.302918e-01_r8,6.288315e-01_r8,6.274512e-01_r8,&
        & 6.261436e-01_r8,6.249022e-01_r8,6.237211e-01_r8,6.225953e-01_r8,6.215201e-01_r8,&
        & 6.204914e-01_r8,6.195055e-01_r8,6.185592e-01_r8,6.176492e-01_r8,6.167730e-01_r8,&
        & 6.159280e-01_r8,6.151120e-01_r8,6.143228e-01_r8,6.135587e-01_r8,6.128177e-01_r8,&
        & 6.120984e-01_r8,6.113993e-01_r8,6.107189e-01_r8 /)
      ssaice2(:, 18) = (/ &
! band 18
        & 9.956167e-01_r8,9.814770e-01_r8,9.716104e-01_r8,9.639746e-01_r8,9.577179e-01_r8,&
        & 9.524010e-01_r8,9.477672e-01_r8,9.436527e-01_r8,9.399467e-01_r8,9.365708e-01_r8,&
        & 9.334672e-01_r8,9.305921e-01_r8,9.279118e-01_r8,9.253993e-01_r8,9.230330e-01_r8,&
        & 9.207954e-01_r8,9.186719e-01_r8,9.166501e-01_r8,9.147199e-01_r8,9.128722e-01_r8,&
        & 9.110997e-01_r8,9.093956e-01_r8,9.077544e-01_r8,9.061708e-01_r8,9.046406e-01_r8,&
        & 9.031598e-01_r8,9.017248e-01_r8,9.003326e-01_r8,8.989804e-01_r8,8.976655e-01_r8,&
        & 8.963857e-01_r8,8.951389e-01_r8,8.939233e-01_r8,8.927370e-01_r8,8.915785e-01_r8,&
        & 8.904464e-01_r8,8.893392e-01_r8,8.882559e-01_r8,8.871951e-01_r8,8.861559e-01_r8,&
        & 8.851373e-01_r8,8.841383e-01_r8,8.831581e-01_r8 /)
      ssaice2(:, 19) = (/ &
! band 19
        & 9.723177e-01_r8,9.452119e-01_r8,9.267592e-01_r8,9.127393e-01_r8,9.014238e-01_r8,&
        & 8.919334e-01_r8,8.837584e-01_r8,8.765773e-01_r8,8.701736e-01_r8,8.643950e-01_r8,&
        & 8.591299e-01_r8,8.542942e-01_r8,8.498230e-01_r8,8.456651e-01_r8,8.417794e-01_r8,&
        & 8.381324e-01_r8,8.346964e-01_r8,8.314484e-01_r8,8.283687e-01_r8,8.254408e-01_r8,&
        & 8.226505e-01_r8,8.199854e-01_r8,8.174348e-01_r8,8.149891e-01_r8,8.126403e-01_r8,&
        & 8.103808e-01_r8,8.082041e-01_r8,8.061044e-01_r8,8.040765e-01_r8,8.021156e-01_r8,&
        & 8.002174e-01_r8,7.983781e-01_r8,7.965941e-01_r8,7.948622e-01_r8,7.931795e-01_r8,&
        & 7.915432e-01_r8,7.899508e-01_r8,7.884002e-01_r8,7.868891e-01_r8,7.854156e-01_r8,&
        & 7.839779e-01_r8,7.825742e-01_r8,7.812031e-01_r8 /)
      ssaice2(:, 20) = (/ &
! band 20
        & 9.933294e-01_r8,9.860917e-01_r8,9.811564e-01_r8,9.774008e-01_r8,9.743652e-01_r8,&
        & 9.718155e-01_r8,9.696159e-01_r8,9.676810e-01_r8,9.659531e-01_r8,9.643915e-01_r8,&
        & 9.629667e-01_r8,9.616561e-01_r8,9.604426e-01_r8,9.593125e-01_r8,9.582548e-01_r8,&
        & 9.572607e-01_r8,9.563227e-01_r8,9.554347e-01_r8,9.545915e-01_r8,9.537888e-01_r8,&
        & 9.530226e-01_r8,9.522898e-01_r8,9.515874e-01_r8,9.509130e-01_r8,9.502643e-01_r8,&
        & 9.496394e-01_r8,9.490366e-01_r8,9.484542e-01_r8,9.478910e-01_r8,9.473456e-01_r8,&
        & 9.468169e-01_r8,9.463039e-01_r8,9.458056e-01_r8,9.453212e-01_r8,9.448499e-01_r8,&
        & 9.443910e-01_r8,9.439438e-01_r8,9.435077e-01_r8,9.430821e-01_r8,9.426666e-01_r8,&
        & 9.422607e-01_r8,9.418638e-01_r8,9.414756e-01_r8 /)
      ssaice2(:, 21) = (/ &
! band 21
        & 9.900787e-01_r8,9.828880e-01_r8,9.779258e-01_r8,9.741173e-01_r8,9.710184e-01_r8,&
        & 9.684012e-01_r8,9.661332e-01_r8,9.641301e-01_r8,9.623352e-01_r8,9.607083e-01_r8,&
        & 9.592198e-01_r8,9.578474e-01_r8,9.565739e-01_r8,9.553856e-01_r8,9.542715e-01_r8,&
        & 9.532226e-01_r8,9.522314e-01_r8,9.512919e-01_r8,9.503986e-01_r8,9.495472e-01_r8,&
        & 9.487337e-01_r8,9.479549e-01_r8,9.472077e-01_r8,9.464897e-01_r8,9.457985e-01_r8,&
        & 9.451322e-01_r8,9.444890e-01_r8,9.438673e-01_r8,9.432656e-01_r8,9.426826e-01_r8,&
        & 9.421173e-01_r8,9.415684e-01_r8,9.410351e-01_r8,9.405164e-01_r8,9.400115e-01_r8,&
        & 9.395198e-01_r8,9.390404e-01_r8,9.385728e-01_r8,9.381164e-01_r8,9.376707e-01_r8,&
        & 9.372350e-01_r8,9.368091e-01_r8,9.363923e-01_r8 /)
      ssaice2(:, 22) = (/ &
! band 22
        & 9.986793e-01_r8,9.985239e-01_r8,9.983911e-01_r8,9.982715e-01_r8,9.981606e-01_r8,&
        & 9.980562e-01_r8,9.979567e-01_r8,9.978613e-01_r8,9.977691e-01_r8,9.976798e-01_r8,&
        & 9.975929e-01_r8,9.975081e-01_r8,9.974251e-01_r8,9.973438e-01_r8,9.972640e-01_r8,&
        & 9.971855e-01_r8,9.971083e-01_r8,9.970322e-01_r8,9.969571e-01_r8,9.968830e-01_r8,&
        & 9.968099e-01_r8,9.967375e-01_r8,9.966660e-01_r8,9.965951e-01_r8,9.965250e-01_r8,&
        & 9.964555e-01_r8,9.963867e-01_r8,9.963185e-01_r8,9.962508e-01_r8,9.961836e-01_r8,&
        & 9.961170e-01_r8,9.960508e-01_r8,9.959851e-01_r8,9.959198e-01_r8,9.958550e-01_r8,&
        & 9.957906e-01_r8,9.957266e-01_r8,9.956629e-01_r8,9.955997e-01_r8,9.955367e-01_r8,&
        & 9.954742e-01_r8,9.954119e-01_r8,9.953500e-01_r8 /)
      ssaice2(:, 23) = (/ &
! band 23
        & 9.997944e-01_r8,9.997791e-01_r8,9.997664e-01_r8,9.997547e-01_r8,9.997436e-01_r8,&
        & 9.997327e-01_r8,9.997219e-01_r8,9.997110e-01_r8,9.996999e-01_r8,9.996886e-01_r8,&
        & 9.996771e-01_r8,9.996653e-01_r8,9.996533e-01_r8,9.996409e-01_r8,9.996282e-01_r8,&
        & 9.996152e-01_r8,9.996019e-01_r8,9.995883e-01_r8,9.995743e-01_r8,9.995599e-01_r8,&
        & 9.995453e-01_r8,9.995302e-01_r8,9.995149e-01_r8,9.994992e-01_r8,9.994831e-01_r8,&
        & 9.994667e-01_r8,9.994500e-01_r8,9.994329e-01_r8,9.994154e-01_r8,9.993976e-01_r8,&
        & 9.993795e-01_r8,9.993610e-01_r8,9.993422e-01_r8,9.993230e-01_r8,9.993035e-01_r8,&
        & 9.992837e-01_r8,9.992635e-01_r8,9.992429e-01_r8,9.992221e-01_r8,9.992008e-01_r8,&
        & 9.991793e-01_r8,9.991574e-01_r8,9.991352e-01_r8 /)
      ssaice2(:, 24) = (/ &
! band 24
        & 9.999949e-01_r8,9.999947e-01_r8,9.999943e-01_r8,9.999939e-01_r8,9.999934e-01_r8,&
        & 9.999927e-01_r8,9.999920e-01_r8,9.999913e-01_r8,9.999904e-01_r8,9.999895e-01_r8,&
        & 9.999885e-01_r8,9.999874e-01_r8,9.999863e-01_r8,9.999851e-01_r8,9.999838e-01_r8,&
        & 9.999824e-01_r8,9.999810e-01_r8,9.999795e-01_r8,9.999780e-01_r8,9.999764e-01_r8,&
        & 9.999747e-01_r8,9.999729e-01_r8,9.999711e-01_r8,9.999692e-01_r8,9.999673e-01_r8,&
        & 9.999653e-01_r8,9.999632e-01_r8,9.999611e-01_r8,9.999589e-01_r8,9.999566e-01_r8,&
        & 9.999543e-01_r8,9.999519e-01_r8,9.999495e-01_r8,9.999470e-01_r8,9.999444e-01_r8,&
        & 9.999418e-01_r8,9.999392e-01_r8,9.999364e-01_r8,9.999336e-01_r8,9.999308e-01_r8,&
        & 9.999279e-01_r8,9.999249e-01_r8,9.999219e-01_r8 /)
      ssaice2(:, 25) = (/ &
! band 25
        & 9.999997e-01_r8,9.999997e-01_r8,9.999997e-01_r8,9.999996e-01_r8,9.999996e-01_r8,&
        & 9.999995e-01_r8,9.999994e-01_r8,9.999993e-01_r8,9.999993e-01_r8,9.999992e-01_r8,&
        & 9.999991e-01_r8,9.999989e-01_r8,9.999988e-01_r8,9.999987e-01_r8,9.999986e-01_r8,&
        & 9.999984e-01_r8,9.999983e-01_r8,9.999981e-01_r8,9.999980e-01_r8,9.999978e-01_r8,&
        & 9.999976e-01_r8,9.999974e-01_r8,9.999972e-01_r8,9.999971e-01_r8,9.999969e-01_r8,&
        & 9.999966e-01_r8,9.999964e-01_r8,9.999962e-01_r8,9.999960e-01_r8,9.999957e-01_r8,&
        & 9.999955e-01_r8,9.999953e-01_r8,9.999950e-01_r8,9.999947e-01_r8,9.999945e-01_r8,&
        & 9.999942e-01_r8,9.999939e-01_r8,9.999936e-01_r8,9.999934e-01_r8,9.999931e-01_r8,&
        & 9.999928e-01_r8,9.999925e-01_r8,9.999921e-01_r8 /)
      ssaice2(:, 26) = (/ &
! band 26
        & 9.999997e-01_r8,9.999996e-01_r8,9.999996e-01_r8,9.999995e-01_r8,9.999994e-01_r8,&
        & 9.999993e-01_r8,9.999992e-01_r8,9.999991e-01_r8,9.999990e-01_r8,9.999989e-01_r8,&
        & 9.999987e-01_r8,9.999986e-01_r8,9.999984e-01_r8,9.999982e-01_r8,9.999980e-01_r8,&
        & 9.999978e-01_r8,9.999976e-01_r8,9.999974e-01_r8,9.999972e-01_r8,9.999970e-01_r8,&
        & 9.999967e-01_r8,9.999965e-01_r8,9.999962e-01_r8,9.999959e-01_r8,9.999956e-01_r8,&
        & 9.999954e-01_r8,9.999951e-01_r8,9.999947e-01_r8,9.999944e-01_r8,9.999941e-01_r8,&
        & 9.999938e-01_r8,9.999934e-01_r8,9.999931e-01_r8,9.999927e-01_r8,9.999923e-01_r8,&
        & 9.999920e-01_r8,9.999916e-01_r8,9.999912e-01_r8,9.999908e-01_r8,9.999904e-01_r8,&
        & 9.999899e-01_r8,9.999895e-01_r8,9.999891e-01_r8 /)
      ssaice2(:, 27) = (/ &
! band 27
        & 9.999987e-01_r8,9.999987e-01_r8,9.999985e-01_r8,9.999984e-01_r8,9.999982e-01_r8,&
        & 9.999980e-01_r8,9.999978e-01_r8,9.999976e-01_r8,9.999973e-01_r8,9.999970e-01_r8,&
        & 9.999967e-01_r8,9.999964e-01_r8,9.999960e-01_r8,9.999956e-01_r8,9.999952e-01_r8,&
        & 9.999948e-01_r8,9.999944e-01_r8,9.999939e-01_r8,9.999934e-01_r8,9.999929e-01_r8,&
        & 9.999924e-01_r8,9.999918e-01_r8,9.999913e-01_r8,9.999907e-01_r8,9.999901e-01_r8,&
        & 9.999894e-01_r8,9.999888e-01_r8,9.999881e-01_r8,9.999874e-01_r8,9.999867e-01_r8,&
        & 9.999860e-01_r8,9.999853e-01_r8,9.999845e-01_r8,9.999837e-01_r8,9.999829e-01_r8,&
        & 9.999821e-01_r8,9.999813e-01_r8,9.999804e-01_r8,9.999796e-01_r8,9.999787e-01_r8,&
        & 9.999778e-01_r8,9.999768e-01_r8,9.999759e-01_r8 /)
      ssaice2(:, 28) = (/ &
! band 28
        & 9.999989e-01_r8,9.999989e-01_r8,9.999987e-01_r8,9.999986e-01_r8,9.999984e-01_r8,&
        & 9.999982e-01_r8,9.999980e-01_r8,9.999978e-01_r8,9.999975e-01_r8,9.999972e-01_r8,&
        & 9.999969e-01_r8,9.999966e-01_r8,9.999962e-01_r8,9.999958e-01_r8,9.999954e-01_r8,&
        & 9.999950e-01_r8,9.999945e-01_r8,9.999941e-01_r8,9.999936e-01_r8,9.999931e-01_r8,&
        & 9.999925e-01_r8,9.999920e-01_r8,9.999914e-01_r8,9.999908e-01_r8,9.999902e-01_r8,&
        & 9.999896e-01_r8,9.999889e-01_r8,9.999883e-01_r8,9.999876e-01_r8,9.999869e-01_r8,&
        & 9.999861e-01_r8,9.999854e-01_r8,9.999846e-01_r8,9.999838e-01_r8,9.999830e-01_r8,&
        & 9.999822e-01_r8,9.999814e-01_r8,9.999805e-01_r8,9.999796e-01_r8,9.999787e-01_r8,&
        & 9.999778e-01_r8,9.999769e-01_r8,9.999759e-01_r8 /)
      ssaice2(:, 29) = (/ &
! band 29
        & 7.042143e-01_r8,6.691161e-01_r8,6.463240e-01_r8,6.296590e-01_r8,6.166381e-01_r8,&
        & 6.060183e-01_r8,5.970908e-01_r8,5.894144e-01_r8,5.826968e-01_r8,5.767343e-01_r8,&
        & 5.713804e-01_r8,5.665256e-01_r8,5.620867e-01_r8,5.579987e-01_r8,5.542101e-01_r8,&
        & 5.506794e-01_r8,5.473727e-01_r8,5.442620e-01_r8,5.413239e-01_r8,5.385389e-01_r8,&
        & 5.358901e-01_r8,5.333633e-01_r8,5.309460e-01_r8,5.286277e-01_r8,5.263988e-01_r8,&
        & 5.242512e-01_r8,5.221777e-01_r8,5.201719e-01_r8,5.182280e-01_r8,5.163410e-01_r8,&
        & 5.145062e-01_r8,5.127197e-01_r8,5.109776e-01_r8,5.092766e-01_r8,5.076137e-01_r8,&
        & 5.059860e-01_r8,5.043911e-01_r8,5.028266e-01_r8,5.012904e-01_r8,4.997805e-01_r8,&
        & 4.982951e-01_r8,4.968326e-01_r8,4.953913e-01_r8 /)

! asymmetry factor: unitless
      asyice2(:, 16) = (/ &
! band 16
        & 7.946655e-01_r8,8.547685e-01_r8,8.806016e-01_r8,8.949880e-01_r8,9.041676e-01_r8,&
        & 9.105399e-01_r8,9.152249e-01_r8,9.188160e-01_r8,9.216573e-01_r8,9.239620e-01_r8,&
        & 9.258695e-01_r8,9.274745e-01_r8,9.288441e-01_r8,9.300267e-01_r8,9.310584e-01_r8,&
        & 9.319665e-01_r8,9.327721e-01_r8,9.334918e-01_r8,9.341387e-01_r8,9.347236e-01_r8,&
        & 9.352551e-01_r8,9.357402e-01_r8,9.361850e-01_r8,9.365942e-01_r8,9.369722e-01_r8,&
        & 9.373225e-01_r8,9.376481e-01_r8,9.379516e-01_r8,9.382352e-01_r8,9.385010e-01_r8,&
        & 9.387505e-01_r8,9.389854e-01_r8,9.392070e-01_r8,9.394163e-01_r8,9.396145e-01_r8,&
        & 9.398024e-01_r8,9.399809e-01_r8,9.401508e-01_r8,9.403126e-01_r8,9.404670e-01_r8,&
        & 9.406144e-01_r8,9.407555e-01_r8,9.408906e-01_r8 /)
      asyice2(:, 17) = (/ &
! band 17
        & 9.078091e-01_r8,9.195850e-01_r8,9.267250e-01_r8,9.317083e-01_r8,9.354632e-01_r8,&
        & 9.384323e-01_r8,9.408597e-01_r8,9.428935e-01_r8,9.446301e-01_r8,9.461351e-01_r8,&
        & 9.474555e-01_r8,9.486259e-01_r8,9.496722e-01_r8,9.506146e-01_r8,9.514688e-01_r8,&
        & 9.522476e-01_r8,9.529612e-01_r8,9.536181e-01_r8,9.542251e-01_r8,9.547883e-01_r8,&
        & 9.553124e-01_r8,9.558019e-01_r8,9.562601e-01_r8,9.566904e-01_r8,9.570953e-01_r8,&
        & 9.574773e-01_r8,9.578385e-01_r8,9.581806e-01_r8,9.585054e-01_r8,9.588142e-01_r8,&
        & 9.591083e-01_r8,9.593888e-01_r8,9.596569e-01_r8,9.599135e-01_r8,9.601593e-01_r8,&
        & 9.603952e-01_r8,9.606219e-01_r8,9.608399e-01_r8,9.610499e-01_r8,9.612523e-01_r8,&
        & 9.614477e-01_r8,9.616365e-01_r8,9.618192e-01_r8 /)
      asyice2(:, 18) = (/ &
! band 18
        & 8.322045e-01_r8,8.528693e-01_r8,8.648167e-01_r8,8.729163e-01_r8,8.789054e-01_r8,&
        & 8.835845e-01_r8,8.873819e-01_r8,8.905511e-01_r8,8.932532e-01_r8,8.955965e-01_r8,&
        & 8.976567e-01_r8,8.994887e-01_r8,9.011334e-01_r8,9.026221e-01_r8,9.039791e-01_r8,&
        & 9.052237e-01_r8,9.063715e-01_r8,9.074349e-01_r8,9.084245e-01_r8,9.093489e-01_r8,&
        & 9.102154e-01_r8,9.110303e-01_r8,9.117987e-01_r8,9.125253e-01_r8,9.132140e-01_r8,&
        & 9.138682e-01_r8,9.144910e-01_r8,9.150850e-01_r8,9.156524e-01_r8,9.161955e-01_r8,&
        & 9.167160e-01_r8,9.172157e-01_r8,9.176959e-01_r8,9.181581e-01_r8,9.186034e-01_r8,&
        & 9.190330e-01_r8,9.194478e-01_r8,9.198488e-01_r8,9.202368e-01_r8,9.206126e-01_r8,&
        & 9.209768e-01_r8,9.213301e-01_r8,9.216731e-01_r8 /)
      asyice2(:, 19) = (/ &
! band 19
        & 8.116560e-01_r8,8.488278e-01_r8,8.674331e-01_r8,8.788148e-01_r8,8.865810e-01_r8,&
        & 8.922595e-01_r8,8.966149e-01_r8,9.000747e-01_r8,9.028980e-01_r8,9.052513e-01_r8,&
        & 9.072468e-01_r8,9.089632e-01_r8,9.104574e-01_r8,9.117713e-01_r8,9.129371e-01_r8,&
        & 9.139793e-01_r8,9.149174e-01_r8,9.157668e-01_r8,9.165400e-01_r8,9.172473e-01_r8,&
        & 9.178970e-01_r8,9.184962e-01_r8,9.190508e-01_r8,9.195658e-01_r8,9.200455e-01_r8,&
        & 9.204935e-01_r8,9.209130e-01_r8,9.213067e-01_r8,9.216771e-01_r8,9.220262e-01_r8,&
        & 9.223560e-01_r8,9.226680e-01_r8,9.229636e-01_r8,9.232443e-01_r8,9.235112e-01_r8,&
        & 9.237652e-01_r8,9.240074e-01_r8,9.242385e-01_r8,9.244594e-01_r8,9.246708e-01_r8,&
        & 9.248733e-01_r8,9.250674e-01_r8,9.252536e-01_r8 /)
      asyice2(:, 20) = (/ &
! band 20
        & 8.047113e-01_r8,8.402864e-01_r8,8.570332e-01_r8,8.668455e-01_r8,8.733206e-01_r8,&
        & 8.779272e-01_r8,8.813796e-01_r8,8.840676e-01_r8,8.862225e-01_r8,8.879904e-01_r8,&
        & 8.894682e-01_r8,8.907228e-01_r8,8.918019e-01_r8,8.927404e-01_r8,8.935645e-01_r8,&
        & 8.942943e-01_r8,8.949452e-01_r8,8.955296e-01_r8,8.960574e-01_r8,8.965366e-01_r8,&
        & 8.969736e-01_r8,8.973740e-01_r8,8.977422e-01_r8,8.980820e-01_r8,8.983966e-01_r8,&
        & 8.986889e-01_r8,8.989611e-01_r8,8.992153e-01_r8,8.994533e-01_r8,8.996766e-01_r8,&
        & 8.998865e-01_r8,9.000843e-01_r8,9.002709e-01_r8,9.004474e-01_r8,9.006146e-01_r8,&
        & 9.007731e-01_r8,9.009237e-01_r8,9.010670e-01_r8,9.012034e-01_r8,9.013336e-01_r8,&
        & 9.014579e-01_r8,9.015767e-01_r8,9.016904e-01_r8 /)
      asyice2(:, 21) = (/ &
! band 21
        & 8.179122e-01_r8,8.480726e-01_r8,8.621945e-01_r8,8.704354e-01_r8,8.758555e-01_r8,&
        & 8.797007e-01_r8,8.825750e-01_r8,8.848078e-01_r8,8.865939e-01_r8,8.880564e-01_r8,&
        & 8.892765e-01_r8,8.903105e-01_r8,8.911982e-01_r8,8.919689e-01_r8,8.926446e-01_r8,&
        & 8.932419e-01_r8,8.937738e-01_r8,8.942506e-01_r8,8.946806e-01_r8,8.950702e-01_r8,&
        & 8.954251e-01_r8,8.957497e-01_r8,8.960477e-01_r8,8.963223e-01_r8,8.965762e-01_r8,&
        & 8.968116e-01_r8,8.970306e-01_r8,8.972347e-01_r8,8.974255e-01_r8,8.976042e-01_r8,&
        & 8.977720e-01_r8,8.979298e-01_r8,8.980784e-01_r8,8.982188e-01_r8,8.983515e-01_r8,&
        & 8.984771e-01_r8,8.985963e-01_r8,8.987095e-01_r8,8.988171e-01_r8,8.989195e-01_r8,&
        & 8.990172e-01_r8,8.991104e-01_r8,8.991994e-01_r8 /)
      asyice2(:, 22) = (/ &
! band 22
        & 8.169789e-01_r8,8.455024e-01_r8,8.586925e-01_r8,8.663283e-01_r8,8.713217e-01_r8,&
        & 8.748488e-01_r8,8.774765e-01_r8,8.795122e-01_r8,8.811370e-01_r8,8.824649e-01_r8,&
        & 8.835711e-01_r8,8.845073e-01_r8,8.853103e-01_r8,8.860068e-01_r8,8.866170e-01_r8,&
        & 8.871560e-01_r8,8.876358e-01_r8,8.880658e-01_r8,8.884533e-01_r8,8.888044e-01_r8,&
        & 8.891242e-01_r8,8.894166e-01_r8,8.896851e-01_r8,8.899324e-01_r8,8.901612e-01_r8,&
        & 8.903733e-01_r8,8.905706e-01_r8,8.907545e-01_r8,8.909265e-01_r8,8.910876e-01_r8,&
        & 8.912388e-01_r8,8.913812e-01_r8,8.915153e-01_r8,8.916419e-01_r8,8.917617e-01_r8,&
        & 8.918752e-01_r8,8.919829e-01_r8,8.920851e-01_r8,8.921824e-01_r8,8.922751e-01_r8,&
        & 8.923635e-01_r8,8.924478e-01_r8,8.925284e-01_r8 /)
      asyice2(:, 23) = (/ &
! band 23
        & 8.387642e-01_r8,8.569979e-01_r8,8.658630e-01_r8,8.711825e-01_r8,8.747605e-01_r8,&
        & 8.773472e-01_r8,8.793129e-01_r8,8.808621e-01_r8,8.821179e-01_r8,8.831583e-01_r8,&
        & 8.840361e-01_r8,8.847875e-01_r8,8.854388e-01_r8,8.860094e-01_r8,8.865138e-01_r8,&
        & 8.869634e-01_r8,8.873668e-01_r8,8.877310e-01_r8,8.880617e-01_r8,8.883635e-01_r8,&
        & 8.886401e-01_r8,8.888947e-01_r8,8.891298e-01_r8,8.893477e-01_r8,8.895504e-01_r8,&
        & 8.897393e-01_r8,8.899159e-01_r8,8.900815e-01_r8,8.902370e-01_r8,8.903833e-01_r8,&
        & 8.905214e-01_r8,8.906518e-01_r8,8.907753e-01_r8,8.908924e-01_r8,8.910036e-01_r8,&
        & 8.911094e-01_r8,8.912101e-01_r8,8.913062e-01_r8,8.913979e-01_r8,8.914856e-01_r8,&
        & 8.915695e-01_r8,8.916498e-01_r8,8.917269e-01_r8 /)
      asyice2(:, 24) = (/ &
! band 24
        & 8.522208e-01_r8,8.648132e-01_r8,8.711224e-01_r8,8.749901e-01_r8,8.776354e-01_r8,&
        & 8.795743e-01_r8,8.810649e-01_r8,8.822518e-01_r8,8.832225e-01_r8,8.840333e-01_r8,&
        & 8.847224e-01_r8,8.853162e-01_r8,8.858342e-01_r8,8.862906e-01_r8,8.866962e-01_r8,&
        & 8.870595e-01_r8,8.873871e-01_r8,8.876842e-01_r8,8.879551e-01_r8,8.882032e-01_r8,&
        & 8.884316e-01_r8,8.886425e-01_r8,8.888380e-01_r8,8.890199e-01_r8,8.891895e-01_r8,&
        & 8.893481e-01_r8,8.894968e-01_r8,8.896366e-01_r8,8.897683e-01_r8,8.898926e-01_r8,&
        & 8.900102e-01_r8,8.901215e-01_r8,8.902272e-01_r8,8.903276e-01_r8,8.904232e-01_r8,&
        & 8.905144e-01_r8,8.906014e-01_r8,8.906845e-01_r8,8.907640e-01_r8,8.908402e-01_r8,&
        & 8.909132e-01_r8,8.909834e-01_r8,8.910507e-01_r8 /)
      asyice2(:, 25) = (/ &
! band 25
        & 8.578202e-01_r8,8.683033e-01_r8,8.735431e-01_r8,8.767488e-01_r8,8.789378e-01_r8,&
        & 8.805399e-01_r8,8.817701e-01_r8,8.827485e-01_r8,8.835480e-01_r8,8.842152e-01_r8,&
        & 8.847817e-01_r8,8.852696e-01_r8,8.856949e-01_r8,8.860694e-01_r8,8.864020e-01_r8,&
        & 8.866997e-01_r8,8.869681e-01_r8,8.872113e-01_r8,8.874330e-01_r8,8.876360e-01_r8,&
        & 8.878227e-01_r8,8.879951e-01_r8,8.881548e-01_r8,8.883033e-01_r8,8.884418e-01_r8,&
        & 8.885712e-01_r8,8.886926e-01_r8,8.888066e-01_r8,8.889139e-01_r8,8.890152e-01_r8,&
        & 8.891110e-01_r8,8.892017e-01_r8,8.892877e-01_r8,8.893695e-01_r8,8.894473e-01_r8,&
        & 8.895214e-01_r8,8.895921e-01_r8,8.896597e-01_r8,8.897243e-01_r8,8.897862e-01_r8,&
        & 8.898456e-01_r8,8.899025e-01_r8,8.899572e-01_r8 /)
      asyice2(:, 26) = (/ &
! band 26
        & 8.625615e-01_r8,8.713831e-01_r8,8.755799e-01_r8,8.780560e-01_r8,8.796983e-01_r8,&
        & 8.808714e-01_r8,8.817534e-01_r8,8.824420e-01_r8,8.829953e-01_r8,8.834501e-01_r8,&
        & 8.838310e-01_r8,8.841549e-01_r8,8.844338e-01_r8,8.846767e-01_r8,8.848902e-01_r8,&
        & 8.850795e-01_r8,8.852484e-01_r8,8.854002e-01_r8,8.855374e-01_r8,8.856620e-01_r8,&
        & 8.857758e-01_r8,8.858800e-01_r8,8.859759e-01_r8,8.860644e-01_r8,8.861464e-01_r8,&
        & 8.862225e-01_r8,8.862935e-01_r8,8.863598e-01_r8,8.864218e-01_r8,8.864800e-01_r8,&
        & 8.865347e-01_r8,8.865863e-01_r8,8.866349e-01_r8,8.866809e-01_r8,8.867245e-01_r8,&
        & 8.867658e-01_r8,8.868050e-01_r8,8.868423e-01_r8,8.868778e-01_r8,8.869117e-01_r8,&
        & 8.869440e-01_r8,8.869749e-01_r8,8.870044e-01_r8 /)
      asyice2(:, 27) = (/ &
! band 27
        & 8.587495e-01_r8,8.684764e-01_r8,8.728189e-01_r8,8.752872e-01_r8,8.768846e-01_r8,&
        & 8.780060e-01_r8,8.788386e-01_r8,8.794824e-01_r8,8.799960e-01_r8,8.804159e-01_r8,&
        & 8.807660e-01_r8,8.810626e-01_r8,8.813175e-01_r8,8.815390e-01_r8,8.817335e-01_r8,&
        & 8.819057e-01_r8,8.820593e-01_r8,8.821973e-01_r8,8.823220e-01_r8,8.824353e-01_r8,&
        & 8.825387e-01_r8,8.826336e-01_r8,8.827209e-01_r8,8.828016e-01_r8,8.828764e-01_r8,&
        & 8.829459e-01_r8,8.830108e-01_r8,8.830715e-01_r8,8.831283e-01_r8,8.831817e-01_r8,&
        & 8.832320e-01_r8,8.832795e-01_r8,8.833244e-01_r8,8.833668e-01_r8,8.834071e-01_r8,&
        & 8.834454e-01_r8,8.834817e-01_r8,8.835164e-01_r8,8.835495e-01_r8,8.835811e-01_r8,&
        & 8.836113e-01_r8,8.836402e-01_r8,8.836679e-01_r8 /)
      asyice2(:, 28) = (/ &
! band 28
        & 8.561110e-01_r8,8.678583e-01_r8,8.727554e-01_r8,8.753892e-01_r8,8.770154e-01_r8,&
        & 8.781109e-01_r8,8.788949e-01_r8,8.794812e-01_r8,8.799348e-01_r8,8.802952e-01_r8,&
        & 8.805880e-01_r8,8.808300e-01_r8,8.810331e-01_r8,8.812058e-01_r8,8.813543e-01_r8,&
        & 8.814832e-01_r8,8.815960e-01_r8,8.816956e-01_r8,8.817839e-01_r8,8.818629e-01_r8,&
        & 8.819339e-01_r8,8.819979e-01_r8,8.820560e-01_r8,8.821089e-01_r8,8.821573e-01_r8,&
        & 8.822016e-01_r8,8.822425e-01_r8,8.822801e-01_r8,8.823150e-01_r8,8.823474e-01_r8,&
        & 8.823775e-01_r8,8.824056e-01_r8,8.824318e-01_r8,8.824564e-01_r8,8.824795e-01_r8,&
        & 8.825011e-01_r8,8.825215e-01_r8,8.825408e-01_r8,8.825589e-01_r8,8.825761e-01_r8,&
        & 8.825924e-01_r8,8.826078e-01_r8,8.826224e-01_r8 /)
      asyice2(:, 29) = (/ &
! band 29
        & 8.311124e-01_r8,8.688197e-01_r8,8.900274e-01_r8,9.040696e-01_r8,9.142334e-01_r8,&
        & 9.220181e-01_r8,9.282195e-01_r8,9.333048e-01_r8,9.375689e-01_r8,9.412085e-01_r8,&
        & 9.443604e-01_r8,9.471230e-01_r8,9.495694e-01_r8,9.517549e-01_r8,9.537224e-01_r8,&
        & 9.555057e-01_r8,9.571316e-01_r8,9.586222e-01_r8,9.599952e-01_r8,9.612656e-01_r8,&
        & 9.624458e-01_r8,9.635461e-01_r8,9.645756e-01_r8,9.655418e-01_r8,9.664513e-01_r8,&
        & 9.673098e-01_r8,9.681222e-01_r8,9.688928e-01_r8,9.696256e-01_r8,9.703237e-01_r8,&
        & 9.709903e-01_r8,9.716280e-01_r8,9.722391e-01_r8,9.728258e-01_r8,9.733901e-01_r8,&
        & 9.739336e-01_r8,9.744579e-01_r8,9.749645e-01_r8,9.754546e-01_r8,9.759294e-01_r8,&
        & 9.763901e-01_r8,9.768376e-01_r8,9.772727e-01_r8 /)

! Hexagonal Ice Particle Parameterization
! extinction units (ext coef/iwc): [(m^-1)/(g m^-3)]
      extice3(:, 16) = (/ &
! band 16
        & 5.194013e-01_r8,3.215089e-01_r8,2.327917e-01_r8,1.824424e-01_r8,1.499977e-01_r8,&
        & 1.273492e-01_r8,1.106421e-01_r8,9.780982e-02_r8,8.764435e-02_r8,7.939266e-02_r8,&
        & 7.256081e-02_r8,6.681137e-02_r8,6.190600e-02_r8,5.767154e-02_r8,5.397915e-02_r8,&
        & 5.073102e-02_r8,4.785151e-02_r8,4.528125e-02_r8,4.297296e-02_r8,4.088853e-02_r8,&
        & 3.899690e-02_r8,3.727251e-02_r8,3.569411e-02_r8,3.424393e-02_r8,3.290694e-02_r8,&
        & 3.167040e-02_r8,3.052340e-02_r8,2.945654e-02_r8,2.846172e-02_r8,2.753188e-02_r8,&
        & 2.666085e-02_r8,2.584322e-02_r8,2.507423e-02_r8,2.434967e-02_r8,2.366579e-02_r8,&
        & 2.301926e-02_r8,2.240711e-02_r8,2.182666e-02_r8,2.127551e-02_r8,2.075150e-02_r8,&
        & 2.025267e-02_r8,1.977725e-02_r8,1.932364e-02_r8,1.889035e-02_r8,1.847607e-02_r8,&
        & 1.807956e-02_r8 /)
      extice3(:, 17) = (/ &
! band 17
        & 4.901155e-01_r8,3.065286e-01_r8,2.230800e-01_r8,1.753951e-01_r8,1.445402e-01_r8,&
        & 1.229417e-01_r8,1.069777e-01_r8,9.469760e-02_r8,8.495824e-02_r8,7.704501e-02_r8,&
        & 7.048834e-02_r8,6.496693e-02_r8,6.025353e-02_r8,5.618286e-02_r8,5.263186e-02_r8,&
        & 4.950698e-02_r8,4.673585e-02_r8,4.426164e-02_r8,4.203904e-02_r8,4.003153e-02_r8,&
        & 3.820932e-02_r8,3.654790e-02_r8,3.502688e-02_r8,3.362919e-02_r8,3.234041e-02_r8,&
        & 3.114829e-02_r8,3.004234e-02_r8,2.901356e-02_r8,2.805413e-02_r8,2.715727e-02_r8,&
        & 2.631705e-02_r8,2.552828e-02_r8,2.478637e-02_r8,2.408725e-02_r8,2.342734e-02_r8,&
        & 2.280343e-02_r8,2.221264e-02_r8,2.165242e-02_r8,2.112043e-02_r8,2.061461e-02_r8,&
        & 2.013308e-02_r8,1.967411e-02_r8,1.923616e-02_r8,1.881783e-02_r8,1.841781e-02_r8,&
        & 1.803494e-02_r8 /)
      extice3(:, 18) = (/ &
! band 18
        & 5.056264e-01_r8,3.160261e-01_r8,2.298442e-01_r8,1.805973e-01_r8,1.487318e-01_r8,&
        & 1.264258e-01_r8,1.099389e-01_r8,9.725656e-02_r8,8.719819e-02_r8,7.902576e-02_r8,&
        & 7.225433e-02_r8,6.655206e-02_r8,6.168427e-02_r8,5.748028e-02_r8,5.381296e-02_r8,&
        & 5.058572e-02_r8,4.772383e-02_r8,4.516857e-02_r8,4.287317e-02_r8,4.079990e-02_r8,&
        & 3.891801e-02_r8,3.720217e-02_r8,3.563133e-02_r8,3.418786e-02_r8,3.285686e-02_r8,&
        & 3.162569e-02_r8,3.048352e-02_r8,2.942104e-02_r8,2.843018e-02_r8,2.750395e-02_r8,&
        & 2.663621e-02_r8,2.582160e-02_r8,2.505539e-02_r8,2.433337e-02_r8,2.365185e-02_r8,&
        & 2.300750e-02_r8,2.239736e-02_r8,2.181878e-02_r8,2.126937e-02_r8,2.074699e-02_r8,&
        & 2.024968e-02_r8,1.977567e-02_r8,1.932338e-02_r8,1.889134e-02_r8,1.847823e-02_r8,&
        & 1.808281e-02_r8 /)
      extice3(:, 19) = (/ &
! band 19
        & 4.881605e-01_r8,3.055237e-01_r8,2.225070e-01_r8,1.750688e-01_r8,1.443736e-01_r8,&
        & 1.228869e-01_r8,1.070054e-01_r8,9.478893e-02_r8,8.509997e-02_r8,7.722769e-02_r8,&
        & 7.070495e-02_r8,6.521211e-02_r8,6.052311e-02_r8,5.647351e-02_r8,5.294088e-02_r8,&
        & 4.983217e-02_r8,4.707539e-02_r8,4.461398e-02_r8,4.240288e-02_r8,4.040575e-02_r8,&
        & 3.859298e-02_r8,3.694016e-02_r8,3.542701e-02_r8,3.403655e-02_r8,3.275444e-02_r8,&
        & 3.156849e-02_r8,3.046827e-02_r8,2.944481e-02_r8,2.849034e-02_r8,2.759812e-02_r8,&
        & 2.676226e-02_r8,2.597757e-02_r8,2.523949e-02_r8,2.454400e-02_r8,2.388750e-02_r8,&
        & 2.326682e-02_r8,2.267909e-02_r8,2.212176e-02_r8,2.159253e-02_r8,2.108933e-02_r8,&
        & 2.061028e-02_r8,2.015369e-02_r8,1.971801e-02_r8,1.930184e-02_r8,1.890389e-02_r8,&
        & 1.852300e-02_r8 /)
      extice3(:, 20) = (/ &
! band 20
        & 5.103703e-01_r8,3.188144e-01_r8,2.317435e-01_r8,1.819887e-01_r8,1.497944e-01_r8,&
        & 1.272584e-01_r8,1.106013e-01_r8,9.778822e-02_r8,8.762610e-02_r8,7.936938e-02_r8,&
        & 7.252809e-02_r8,6.676701e-02_r8,6.184901e-02_r8,5.760165e-02_r8,5.389651e-02_r8,&
        & 5.063598e-02_r8,4.774457e-02_r8,4.516295e-02_r8,4.284387e-02_r8,4.074922e-02_r8,&
        & 3.884792e-02_r8,3.711438e-02_r8,3.552734e-02_r8,3.406898e-02_r8,3.272425e-02_r8,&
        & 3.148038e-02_r8,3.032643e-02_r8,2.925299e-02_r8,2.825191e-02_r8,2.731612e-02_r8,&
        & 2.643943e-02_r8,2.561642e-02_r8,2.484230e-02_r8,2.411284e-02_r8,2.342429e-02_r8,&
        & 2.277329e-02_r8,2.215686e-02_r8,2.157231e-02_r8,2.101724e-02_r8,2.048946e-02_r8,&
        & 1.998702e-02_r8,1.950813e-02_r8,1.905118e-02_r8,1.861468e-02_r8,1.819730e-02_r8,&
        & 1.779781e-02_r8 /)
      extice3(:, 21) = (/ &
! band 21
        & 5.031161e-01_r8,3.144511e-01_r8,2.286942e-01_r8,1.796903e-01_r8,1.479819e-01_r8,&
        & 1.257860e-01_r8,1.093803e-01_r8,9.676059e-02_r8,8.675183e-02_r8,7.861971e-02_r8,&
        & 7.188168e-02_r8,6.620754e-02_r8,6.136376e-02_r8,5.718050e-02_r8,5.353127e-02_r8,&
        & 5.031995e-02_r8,4.747218e-02_r8,4.492952e-02_r8,4.264544e-02_r8,4.058240e-02_r8,&
        & 3.870979e-02_r8,3.700242e-02_r8,3.543933e-02_r8,3.400297e-02_r8,3.267854e-02_r8,&
        & 3.145345e-02_r8,3.031691e-02_r8,2.925967e-02_r8,2.827370e-02_r8,2.735203e-02_r8,&
        & 2.648858e-02_r8,2.567798e-02_r8,2.491555e-02_r8,2.419710e-02_r8,2.351893e-02_r8,&
        & 2.287776e-02_r8,2.227063e-02_r8,2.169491e-02_r8,2.114821e-02_r8,2.062840e-02_r8,&
        & 2.013354e-02_r8,1.966188e-02_r8,1.921182e-02_r8,1.878191e-02_r8,1.837083e-02_r8,&
        & 1.797737e-02_r8 /)
      extice3(:, 22) = (/ &
! band 22
        & 4.949453e-01_r8,3.095918e-01_r8,2.253402e-01_r8,1.771964e-01_r8,1.460446e-01_r8,&
        & 1.242383e-01_r8,1.081206e-01_r8,9.572235e-02_r8,8.588928e-02_r8,7.789990e-02_r8,&
        & 7.128013e-02_r8,6.570559e-02_r8,6.094684e-02_r8,5.683701e-02_r8,5.325183e-02_r8,&
        & 5.009688e-02_r8,4.729909e-02_r8,4.480106e-02_r8,4.255708e-02_r8,4.053025e-02_r8,&
        & 3.869051e-02_r8,3.701310e-02_r8,3.547745e-02_r8,3.406631e-02_r8,3.276512e-02_r8,&
        & 3.156153e-02_r8,3.044494e-02_r8,2.940626e-02_r8,2.843759e-02_r8,2.753211e-02_r8,&
        & 2.668381e-02_r8,2.588744e-02_r8,2.513839e-02_r8,2.443255e-02_r8,2.376629e-02_r8,&
        & 2.313637e-02_r8,2.253990e-02_r8,2.197428e-02_r8,2.143718e-02_r8,2.092649e-02_r8,&
        & 2.044032e-02_r8,1.997694e-02_r8,1.953478e-02_r8,1.911241e-02_r8,1.870855e-02_r8,&
        & 1.832199e-02_r8 /)
      extice3(:, 23) = (/ &
! band 23
        & 5.052816e-01_r8,3.157665e-01_r8,2.296233e-01_r8,1.803986e-01_r8,1.485473e-01_r8,&
        & 1.262514e-01_r8,1.097718e-01_r8,9.709524e-02_r8,8.704139e-02_r8,7.887264e-02_r8,&
        & 7.210424e-02_r8,6.640454e-02_r8,6.153894e-02_r8,5.733683e-02_r8,5.367116e-02_r8,&
        & 5.044537e-02_r8,4.758477e-02_r8,4.503066e-02_r8,4.273629e-02_r8,4.066395e-02_r8,&
        & 3.878291e-02_r8,3.706784e-02_r8,3.549771e-02_r8,3.405488e-02_r8,3.272448e-02_r8,&
        & 3.149387e-02_r8,3.035221e-02_r8,2.929020e-02_r8,2.829979e-02_r8,2.737397e-02_r8,&
        & 2.650663e-02_r8,2.569238e-02_r8,2.492651e-02_r8,2.420482e-02_r8,2.352361e-02_r8,&
        & 2.287954e-02_r8,2.226968e-02_r8,2.169136e-02_r8,2.114220e-02_r8,2.062005e-02_r8,&
        & 2.012296e-02_r8,1.964917e-02_r8,1.919709e-02_r8,1.876524e-02_r8,1.835231e-02_r8,&
        & 1.795707e-02_r8 /)
      extice3(:, 24) = (/ &
! band 24
        & 5.042067e-01_r8,3.151195e-01_r8,2.291708e-01_r8,1.800573e-01_r8,1.482779e-01_r8,&
        & 1.260324e-01_r8,1.095900e-01_r8,9.694202e-02_r8,8.691087e-02_r8,7.876056e-02_r8,&
        & 7.200745e-02_r8,6.632062e-02_r8,6.146600e-02_r8,5.727338e-02_r8,5.361599e-02_r8,&
        & 5.039749e-02_r8,4.754334e-02_r8,4.499500e-02_r8,4.270580e-02_r8,4.063815e-02_r8,&
        & 3.876135e-02_r8,3.705016e-02_r8,3.548357e-02_r8,3.404400e-02_r8,3.271661e-02_r8,&
        & 3.148877e-02_r8,3.034969e-02_r8,2.929008e-02_r8,2.830191e-02_r8,2.737818e-02_r8,&
        & 2.651279e-02_r8,2.570039e-02_r8,2.493624e-02_r8,2.421618e-02_r8,2.353650e-02_r8,&
        & 2.289390e-02_r8,2.228541e-02_r8,2.170840e-02_r8,2.116048e-02_r8,2.063950e-02_r8,&
        & 2.014354e-02_r8,1.967082e-02_r8,1.921975e-02_r8,1.878888e-02_r8,1.837688e-02_r8,&
        & 1.798254e-02_r8 /)
      extice3(:, 25) = (/ &
! band 25
        & 5.022507e-01_r8,3.139246e-01_r8,2.283218e-01_r8,1.794059e-01_r8,1.477544e-01_r8,&
        & 1.255984e-01_r8,1.092222e-01_r8,9.662516e-02_r8,8.663439e-02_r8,7.851688e-02_r8,&
        & 7.179095e-02_r8,6.612700e-02_r8,6.129193e-02_r8,5.711618e-02_r8,5.347351e-02_r8,&
        & 5.026796e-02_r8,4.742530e-02_r8,4.488721e-02_r8,4.260724e-02_r8,4.054790e-02_r8,&
        & 3.867866e-02_r8,3.697435e-02_r8,3.541407e-02_r8,3.398029e-02_r8,3.265824e-02_r8,&
        & 3.143535e-02_r8,3.030085e-02_r8,2.924551e-02_r8,2.826131e-02_r8,2.734130e-02_r8,&
        & 2.647939e-02_r8,2.567026e-02_r8,2.490919e-02_r8,2.419203e-02_r8,2.351509e-02_r8,&
        & 2.287507e-02_r8,2.226903e-02_r8,2.169434e-02_r8,2.114862e-02_r8,2.062975e-02_r8,&
        & 2.013578e-02_r8,1.966496e-02_r8,1.921571e-02_r8,1.878658e-02_r8,1.837623e-02_r8,&
        & 1.798348e-02_r8 /)
      extice3(:, 26) = (/ &
! band 26
        & 5.068316e-01_r8,3.166869e-01_r8,2.302576e-01_r8,1.808693e-01_r8,1.489122e-01_r8,&
        & 1.265423e-01_r8,1.100080e-01_r8,9.728926e-02_r8,8.720201e-02_r8,7.900612e-02_r8,&
        & 7.221524e-02_r8,6.649660e-02_r8,6.161484e-02_r8,5.739877e-02_r8,5.372093e-02_r8,&
        & 5.048442e-02_r8,4.761431e-02_r8,4.505172e-02_r8,4.274972e-02_r8,4.067050e-02_r8,&
        & 3.878321e-02_r8,3.706244e-02_r8,3.548710e-02_r8,3.403948e-02_r8,3.270466e-02_r8,&
        & 3.146995e-02_r8,3.032450e-02_r8,2.925897e-02_r8,2.826527e-02_r8,2.733638e-02_r8,&
        & 2.646615e-02_r8,2.564920e-02_r8,2.488078e-02_r8,2.415670e-02_r8,2.347322e-02_r8,&
        & 2.282702e-02_r8,2.221513e-02_r8,2.163489e-02_r8,2.108390e-02_r8,2.056002e-02_r8,&
        & 2.006128e-02_r8,1.958591e-02_r8,1.913232e-02_r8,1.869904e-02_r8,1.828474e-02_r8,&
        & 1.788819e-02_r8 /)
      extice3(:, 27) = (/ &
! band 27
        & 5.077707e-01_r8,3.172636e-01_r8,2.306695e-01_r8,1.811871e-01_r8,1.491691e-01_r8,&
        & 1.267565e-01_r8,1.101907e-01_r8,9.744773e-02_r8,8.734125e-02_r8,7.912973e-02_r8,&
        & 7.232591e-02_r8,6.659637e-02_r8,6.170530e-02_r8,5.748120e-02_r8,5.379634e-02_r8,&
        & 5.055367e-02_r8,4.767809e-02_r8,4.511061e-02_r8,4.280423e-02_r8,4.072104e-02_r8,&
        & 3.883015e-02_r8,3.710611e-02_r8,3.552776e-02_r8,3.407738e-02_r8,3.274002e-02_r8,&
        & 3.150296e-02_r8,3.035532e-02_r8,2.928776e-02_r8,2.829216e-02_r8,2.736150e-02_r8,&
        & 2.648961e-02_r8,2.567111e-02_r8,2.490123e-02_r8,2.417576e-02_r8,2.349098e-02_r8,&
        & 2.284354e-02_r8,2.223049e-02_r8,2.164914e-02_r8,2.109711e-02_r8,2.057222e-02_r8,&
        & 2.007253e-02_r8,1.959626e-02_r8,1.914181e-02_r8,1.870770e-02_r8,1.829261e-02_r8,&
        & 1.789531e-02_r8 /)
      extice3(:, 28) = (/ &
! band 28
        & 5.062281e-01_r8,3.163402e-01_r8,2.300275e-01_r8,1.807060e-01_r8,1.487921e-01_r8,&
        & 1.264523e-01_r8,1.099403e-01_r8,9.723879e-02_r8,8.716516e-02_r8,7.898034e-02_r8,&
        & 7.219863e-02_r8,6.648771e-02_r8,6.161254e-02_r8,5.740217e-02_r8,5.372929e-02_r8,&
        & 5.049716e-02_r8,4.763092e-02_r8,4.507179e-02_r8,4.277290e-02_r8,4.069649e-02_r8,&
        & 3.881175e-02_r8,3.709331e-02_r8,3.552008e-02_r8,3.407442e-02_r8,3.274141e-02_r8,&
        & 3.150837e-02_r8,3.036447e-02_r8,2.930037e-02_r8,2.830801e-02_r8,2.738037e-02_r8,&
        & 2.651132e-02_r8,2.569547e-02_r8,2.492810e-02_r8,2.420499e-02_r8,2.352243e-02_r8,&
        & 2.287710e-02_r8,2.226604e-02_r8,2.168658e-02_r8,2.113634e-02_r8,2.061316e-02_r8,&
        & 2.011510e-02_r8,1.964038e-02_r8,1.918740e-02_r8,1.875471e-02_r8,1.834096e-02_r8,&
        & 1.794495e-02_r8 /)
      extice3(:, 29) = (/ &
! band 29
        & 1.338834e-01_r8,1.924912e-01_r8,1.755523e-01_r8,1.534793e-01_r8,1.343937e-01_r8,&
        & 1.187883e-01_r8,1.060654e-01_r8,9.559106e-02_r8,8.685880e-02_r8,7.948698e-02_r8,&
        & 7.319086e-02_r8,6.775669e-02_r8,6.302215e-02_r8,5.886236e-02_r8,5.517996e-02_r8,&
        & 5.189810e-02_r8,4.895539e-02_r8,4.630225e-02_r8,4.389823e-02_r8,4.171002e-02_r8,&
        & 3.970998e-02_r8,3.787493e-02_r8,3.618537e-02_r8,3.462471e-02_r8,3.317880e-02_r8,&
        & 3.183547e-02_r8,3.058421e-02_r8,2.941590e-02_r8,2.832256e-02_r8,2.729724e-02_r8,&
        & 2.633377e-02_r8,2.542675e-02_r8,2.457136e-02_r8,2.376332e-02_r8,2.299882e-02_r8,&
        & 2.227443e-02_r8,2.158707e-02_r8,2.093400e-02_r8,2.031270e-02_r8,1.972091e-02_r8,&
        & 1.915659e-02_r8,1.861787e-02_r8,1.810304e-02_r8,1.761055e-02_r8,1.713899e-02_r8,&
        & 1.668704e-02_r8 /)

! single-scattering albedo: unitless
      ssaice3(:, 16) = (/ &
! band 16
        & 6.749442e-01_r8,6.649947e-01_r8,6.565828e-01_r8,6.489928e-01_r8,6.420046e-01_r8,&
        & 6.355231e-01_r8,6.294964e-01_r8,6.238901e-01_r8,6.186783e-01_r8,6.138395e-01_r8,&
        & 6.093543e-01_r8,6.052049e-01_r8,6.013742e-01_r8,5.978457e-01_r8,5.946030e-01_r8,&
        & 5.916302e-01_r8,5.889115e-01_r8,5.864310e-01_r8,5.841731e-01_r8,5.821221e-01_r8,&
        & 5.802624e-01_r8,5.785785e-01_r8,5.770549e-01_r8,5.756759e-01_r8,5.744262e-01_r8,&
        & 5.732901e-01_r8,5.722524e-01_r8,5.712974e-01_r8,5.704097e-01_r8,5.695739e-01_r8,&
        & 5.687747e-01_r8,5.679964e-01_r8,5.672238e-01_r8,5.664415e-01_r8,5.656340e-01_r8,&
        & 5.647860e-01_r8,5.638821e-01_r8,5.629070e-01_r8,5.618452e-01_r8,5.606815e-01_r8,&
        & 5.594006e-01_r8,5.579870e-01_r8,5.564255e-01_r8,5.547008e-01_r8,5.527976e-01_r8,&
        & 5.507005e-01_r8 /)
      ssaice3(:, 17) = (/ &
! band 17
        & 7.628550e-01_r8,7.567297e-01_r8,7.508463e-01_r8,7.451972e-01_r8,7.397745e-01_r8,&
        & 7.345705e-01_r8,7.295775e-01_r8,7.247881e-01_r8,7.201945e-01_r8,7.157894e-01_r8,&
        & 7.115652e-01_r8,7.075145e-01_r8,7.036300e-01_r8,6.999044e-01_r8,6.963304e-01_r8,&
        & 6.929007e-01_r8,6.896083e-01_r8,6.864460e-01_r8,6.834067e-01_r8,6.804833e-01_r8,&
        & 6.776690e-01_r8,6.749567e-01_r8,6.723397e-01_r8,6.698109e-01_r8,6.673637e-01_r8,&
        & 6.649913e-01_r8,6.626870e-01_r8,6.604441e-01_r8,6.582561e-01_r8,6.561163e-01_r8,&
        & 6.540182e-01_r8,6.519554e-01_r8,6.499215e-01_r8,6.479099e-01_r8,6.459145e-01_r8,&
        & 6.439289e-01_r8,6.419468e-01_r8,6.399621e-01_r8,6.379686e-01_r8,6.359601e-01_r8,&
        & 6.339306e-01_r8,6.318740e-01_r8,6.297845e-01_r8,6.276559e-01_r8,6.254825e-01_r8,&
        & 6.232583e-01_r8 /)
      ssaice3(:, 18) = (/ &
! band 18
        & 9.924147e-01_r8,9.882792e-01_r8,9.842257e-01_r8,9.802522e-01_r8,9.763566e-01_r8,&
        & 9.725367e-01_r8,9.687905e-01_r8,9.651157e-01_r8,9.615104e-01_r8,9.579725e-01_r8,&
        & 9.544997e-01_r8,9.510901e-01_r8,9.477416e-01_r8,9.444520e-01_r8,9.412194e-01_r8,&
        & 9.380415e-01_r8,9.349165e-01_r8,9.318421e-01_r8,9.288164e-01_r8,9.258373e-01_r8,&
        & 9.229027e-01_r8,9.200106e-01_r8,9.171589e-01_r8,9.143457e-01_r8,9.115688e-01_r8,&
        & 9.088263e-01_r8,9.061161e-01_r8,9.034362e-01_r8,9.007846e-01_r8,8.981592e-01_r8,&
        & 8.955581e-01_r8,8.929792e-01_r8,8.904206e-01_r8,8.878803e-01_r8,8.853562e-01_r8,&
        & 8.828464e-01_r8,8.803488e-01_r8,8.778616e-01_r8,8.753827e-01_r8,8.729102e-01_r8,&
        & 8.704421e-01_r8,8.679764e-01_r8,8.655112e-01_r8,8.630445e-01_r8,8.605744e-01_r8,&
        & 8.580989e-01_r8 /)
      ssaice3(:, 19) = (/ &
! band 19
        & 9.629413e-01_r8,9.517182e-01_r8,9.409209e-01_r8,9.305366e-01_r8,9.205529e-01_r8,&
        & 9.109569e-01_r8,9.017362e-01_r8,8.928780e-01_r8,8.843699e-01_r8,8.761992e-01_r8,&
        & 8.683536e-01_r8,8.608204e-01_r8,8.535873e-01_r8,8.466417e-01_r8,8.399712e-01_r8,&
        & 8.335635e-01_r8,8.274062e-01_r8,8.214868e-01_r8,8.157932e-01_r8,8.103129e-01_r8,&
        & 8.050336e-01_r8,7.999432e-01_r8,7.950294e-01_r8,7.902798e-01_r8,7.856825e-01_r8,&
        & 7.812250e-01_r8,7.768954e-01_r8,7.726815e-01_r8,7.685711e-01_r8,7.645522e-01_r8,&
        & 7.606126e-01_r8,7.567404e-01_r8,7.529234e-01_r8,7.491498e-01_r8,7.454074e-01_r8,&
        & 7.416844e-01_r8,7.379688e-01_r8,7.342485e-01_r8,7.305118e-01_r8,7.267468e-01_r8,&
        & 7.229415e-01_r8,7.190841e-01_r8,7.151628e-01_r8,7.111657e-01_r8,7.070811e-01_r8,&
        & 7.028972e-01_r8 /)
      ssaice3(:, 20) = (/ &
! band 20
        & 9.942270e-01_r8,9.909206e-01_r8,9.876775e-01_r8,9.844960e-01_r8,9.813746e-01_r8,&
        & 9.783114e-01_r8,9.753049e-01_r8,9.723535e-01_r8,9.694553e-01_r8,9.666088e-01_r8,&
        & 9.638123e-01_r8,9.610641e-01_r8,9.583626e-01_r8,9.557060e-01_r8,9.530928e-01_r8,&
        & 9.505211e-01_r8,9.479895e-01_r8,9.454961e-01_r8,9.430393e-01_r8,9.406174e-01_r8,&
        & 9.382288e-01_r8,9.358717e-01_r8,9.335446e-01_r8,9.312456e-01_r8,9.289731e-01_r8,&
        & 9.267255e-01_r8,9.245010e-01_r8,9.222980e-01_r8,9.201147e-01_r8,9.179496e-01_r8,&
        & 9.158008e-01_r8,9.136667e-01_r8,9.115457e-01_r8,9.094359e-01_r8,9.073358e-01_r8,&
        & 9.052436e-01_r8,9.031577e-01_r8,9.010763e-01_r8,8.989977e-01_r8,8.969203e-01_r8,&
        & 8.948423e-01_r8,8.927620e-01_r8,8.906778e-01_r8,8.885879e-01_r8,8.864907e-01_r8,&
        & 8.843843e-01_r8 /)
      ssaice3(:, 21) = (/ &
! band 21
        & 9.934014e-01_r8,9.899331e-01_r8,9.865537e-01_r8,9.832610e-01_r8,9.800523e-01_r8,&
        & 9.769254e-01_r8,9.738777e-01_r8,9.709069e-01_r8,9.680106e-01_r8,9.651862e-01_r8,&
        & 9.624315e-01_r8,9.597439e-01_r8,9.571212e-01_r8,9.545608e-01_r8,9.520605e-01_r8,&
        & 9.496177e-01_r8,9.472301e-01_r8,9.448954e-01_r8,9.426111e-01_r8,9.403749e-01_r8,&
        & 9.381843e-01_r8,9.360370e-01_r8,9.339307e-01_r8,9.318629e-01_r8,9.298313e-01_r8,&
        & 9.278336e-01_r8,9.258673e-01_r8,9.239302e-01_r8,9.220198e-01_r8,9.201338e-01_r8,&
        & 9.182700e-01_r8,9.164258e-01_r8,9.145991e-01_r8,9.127874e-01_r8,9.109884e-01_r8,&
        & 9.091999e-01_r8,9.074194e-01_r8,9.056447e-01_r8,9.038735e-01_r8,9.021033e-01_r8,&
        & 9.003320e-01_r8,8.985572e-01_r8,8.967766e-01_r8,8.949879e-01_r8,8.931888e-01_r8,&
        & 8.913770e-01_r8 /)
      ssaice3(:, 22) = (/ &
! band 22
        & 9.994833e-01_r8,9.992055e-01_r8,9.989278e-01_r8,9.986500e-01_r8,9.983724e-01_r8,&
        & 9.980947e-01_r8,9.978172e-01_r8,9.975397e-01_r8,9.972623e-01_r8,9.969849e-01_r8,&
        & 9.967077e-01_r8,9.964305e-01_r8,9.961535e-01_r8,9.958765e-01_r8,9.955997e-01_r8,&
        & 9.953230e-01_r8,9.950464e-01_r8,9.947699e-01_r8,9.944936e-01_r8,9.942174e-01_r8,&
        & 9.939414e-01_r8,9.936656e-01_r8,9.933899e-01_r8,9.931144e-01_r8,9.928390e-01_r8,&
        & 9.925639e-01_r8,9.922889e-01_r8,9.920141e-01_r8,9.917396e-01_r8,9.914652e-01_r8,&
        & 9.911911e-01_r8,9.909171e-01_r8,9.906434e-01_r8,9.903700e-01_r8,9.900967e-01_r8,&
        & 9.898237e-01_r8,9.895510e-01_r8,9.892784e-01_r8,9.890062e-01_r8,9.887342e-01_r8,&
        & 9.884625e-01_r8,9.881911e-01_r8,9.879199e-01_r8,9.876490e-01_r8,9.873784e-01_r8,&
        & 9.871081e-01_r8 /)
      ssaice3(:, 23) = (/ &
! band 23
        & 9.999343e-01_r8,9.998917e-01_r8,9.998492e-01_r8,9.998067e-01_r8,9.997642e-01_r8,&
        & 9.997218e-01_r8,9.996795e-01_r8,9.996372e-01_r8,9.995949e-01_r8,9.995528e-01_r8,&
        & 9.995106e-01_r8,9.994686e-01_r8,9.994265e-01_r8,9.993845e-01_r8,9.993426e-01_r8,&
        & 9.993007e-01_r8,9.992589e-01_r8,9.992171e-01_r8,9.991754e-01_r8,9.991337e-01_r8,&
        & 9.990921e-01_r8,9.990505e-01_r8,9.990089e-01_r8,9.989674e-01_r8,9.989260e-01_r8,&
        & 9.988846e-01_r8,9.988432e-01_r8,9.988019e-01_r8,9.987606e-01_r8,9.987194e-01_r8,&
        & 9.986782e-01_r8,9.986370e-01_r8,9.985959e-01_r8,9.985549e-01_r8,9.985139e-01_r8,&
        & 9.984729e-01_r8,9.984319e-01_r8,9.983910e-01_r8,9.983502e-01_r8,9.983094e-01_r8,&
        & 9.982686e-01_r8,9.982279e-01_r8,9.981872e-01_r8,9.981465e-01_r8,9.981059e-01_r8,&
        & 9.980653e-01_r8 /)
      ssaice3(:, 24) = (/ &
! band 24
        & 9.999978e-01_r8,9.999965e-01_r8,9.999952e-01_r8,9.999939e-01_r8,9.999926e-01_r8,&
        & 9.999913e-01_r8,9.999900e-01_r8,9.999887e-01_r8,9.999873e-01_r8,9.999860e-01_r8,&
        & 9.999847e-01_r8,9.999834e-01_r8,9.999821e-01_r8,9.999808e-01_r8,9.999795e-01_r8,&
        & 9.999782e-01_r8,9.999769e-01_r8,9.999756e-01_r8,9.999743e-01_r8,9.999730e-01_r8,&
        & 9.999717e-01_r8,9.999704e-01_r8,9.999691e-01_r8,9.999678e-01_r8,9.999665e-01_r8,&
        & 9.999652e-01_r8,9.999639e-01_r8,9.999626e-01_r8,9.999613e-01_r8,9.999600e-01_r8,&
        & 9.999587e-01_r8,9.999574e-01_r8,9.999561e-01_r8,9.999548e-01_r8,9.999535e-01_r8,&
        & 9.999522e-01_r8,9.999509e-01_r8,9.999496e-01_r8,9.999483e-01_r8,9.999470e-01_r8,&
        & 9.999457e-01_r8,9.999444e-01_r8,9.999431e-01_r8,9.999418e-01_r8,9.999405e-01_r8,&
        & 9.999392e-01_r8 /)
      ssaice3(:, 25) = (/ &
! band 25
        & 9.999994e-01_r8,9.999993e-01_r8,9.999991e-01_r8,9.999990e-01_r8,9.999989e-01_r8,&
        & 9.999987e-01_r8,9.999986e-01_r8,9.999984e-01_r8,9.999983e-01_r8,9.999982e-01_r8,&
        & 9.999980e-01_r8,9.999979e-01_r8,9.999977e-01_r8,9.999976e-01_r8,9.999975e-01_r8,&
        & 9.999973e-01_r8,9.999972e-01_r8,9.999970e-01_r8,9.999969e-01_r8,9.999967e-01_r8,&
        & 9.999966e-01_r8,9.999965e-01_r8,9.999963e-01_r8,9.999962e-01_r8,9.999960e-01_r8,&
        & 9.999959e-01_r8,9.999957e-01_r8,9.999956e-01_r8,9.999954e-01_r8,9.999953e-01_r8,&
        & 9.999952e-01_r8,9.999950e-01_r8,9.999949e-01_r8,9.999947e-01_r8,9.999946e-01_r8,&
        & 9.999944e-01_r8,9.999943e-01_r8,9.999941e-01_r8,9.999940e-01_r8,9.999939e-01_r8,&
        & 9.999937e-01_r8,9.999936e-01_r8,9.999934e-01_r8,9.999933e-01_r8,9.999931e-01_r8,&
        & 9.999930e-01_r8 /)
      ssaice3(:, 26) = (/ &
! band 26
        & 9.999997e-01_r8,9.999995e-01_r8,9.999992e-01_r8,9.999990e-01_r8,9.999987e-01_r8,&
        & 9.999985e-01_r8,9.999983e-01_r8,9.999980e-01_r8,9.999978e-01_r8,9.999976e-01_r8,&
        & 9.999973e-01_r8,9.999971e-01_r8,9.999969e-01_r8,9.999967e-01_r8,9.999965e-01_r8,&
        & 9.999963e-01_r8,9.999960e-01_r8,9.999958e-01_r8,9.999956e-01_r8,9.999954e-01_r8,&
        & 9.999952e-01_r8,9.999950e-01_r8,9.999948e-01_r8,9.999946e-01_r8,9.999944e-01_r8,&
        & 9.999942e-01_r8,9.999939e-01_r8,9.999937e-01_r8,9.999935e-01_r8,9.999933e-01_r8,&
        & 9.999931e-01_r8,9.999929e-01_r8,9.999927e-01_r8,9.999925e-01_r8,9.999923e-01_r8,&
        & 9.999920e-01_r8,9.999918e-01_r8,9.999916e-01_r8,9.999914e-01_r8,9.999911e-01_r8,&
        & 9.999909e-01_r8,9.999907e-01_r8,9.999905e-01_r8,9.999902e-01_r8,9.999900e-01_r8,&
        & 9.999897e-01_r8 /)
      ssaice3(:, 27) = (/ &
! band 27
        & 9.999991e-01_r8,9.999985e-01_r8,9.999980e-01_r8,9.999974e-01_r8,9.999968e-01_r8,&
        & 9.999963e-01_r8,9.999957e-01_r8,9.999951e-01_r8,9.999946e-01_r8,9.999940e-01_r8,&
        & 9.999934e-01_r8,9.999929e-01_r8,9.999923e-01_r8,9.999918e-01_r8,9.999912e-01_r8,&
        & 9.999907e-01_r8,9.999901e-01_r8,9.999896e-01_r8,9.999891e-01_r8,9.999885e-01_r8,&
        & 9.999880e-01_r8,9.999874e-01_r8,9.999869e-01_r8,9.999863e-01_r8,9.999858e-01_r8,&
        & 9.999853e-01_r8,9.999847e-01_r8,9.999842e-01_r8,9.999836e-01_r8,9.999831e-01_r8,&
        & 9.999826e-01_r8,9.999820e-01_r8,9.999815e-01_r8,9.999809e-01_r8,9.999804e-01_r8,&
        & 9.999798e-01_r8,9.999793e-01_r8,9.999787e-01_r8,9.999782e-01_r8,9.999776e-01_r8,&
        & 9.999770e-01_r8,9.999765e-01_r8,9.999759e-01_r8,9.999754e-01_r8,9.999748e-01_r8,&
        & 9.999742e-01_r8 /)
      ssaice3(:, 28) = (/ &
! band 28
        & 9.999975e-01_r8,9.999961e-01_r8,9.999946e-01_r8,9.999931e-01_r8,9.999917e-01_r8,&
        & 9.999903e-01_r8,9.999888e-01_r8,9.999874e-01_r8,9.999859e-01_r8,9.999845e-01_r8,&
        & 9.999831e-01_r8,9.999816e-01_r8,9.999802e-01_r8,9.999788e-01_r8,9.999774e-01_r8,&
        & 9.999759e-01_r8,9.999745e-01_r8,9.999731e-01_r8,9.999717e-01_r8,9.999702e-01_r8,&
        & 9.999688e-01_r8,9.999674e-01_r8,9.999660e-01_r8,9.999646e-01_r8,9.999631e-01_r8,&
        & 9.999617e-01_r8,9.999603e-01_r8,9.999589e-01_r8,9.999574e-01_r8,9.999560e-01_r8,&
        & 9.999546e-01_r8,9.999532e-01_r8,9.999517e-01_r8,9.999503e-01_r8,9.999489e-01_r8,&
        & 9.999474e-01_r8,9.999460e-01_r8,9.999446e-01_r8,9.999431e-01_r8,9.999417e-01_r8,&
        & 9.999403e-01_r8,9.999388e-01_r8,9.999374e-01_r8,9.999359e-01_r8,9.999345e-01_r8,&
        & 9.999330e-01_r8 /)
      ssaice3(:, 29) = (/ &
! band 29
        & 4.526500e-01_r8,5.287890e-01_r8,5.410487e-01_r8,5.459865e-01_r8,5.485149e-01_r8,&
        & 5.498914e-01_r8,5.505895e-01_r8,5.508310e-01_r8,5.507364e-01_r8,5.503793e-01_r8,&
        & 5.498090e-01_r8,5.490612e-01_r8,5.481637e-01_r8,5.471395e-01_r8,5.460083e-01_r8,&
        & 5.447878e-01_r8,5.434946e-01_r8,5.421442e-01_r8,5.407514e-01_r8,5.393309e-01_r8,&
        & 5.378970e-01_r8,5.364641e-01_r8,5.350464e-01_r8,5.336582e-01_r8,5.323140e-01_r8,&
        & 5.310283e-01_r8,5.298158e-01_r8,5.286914e-01_r8,5.276704e-01_r8,5.267680e-01_r8,&
        & 5.260000e-01_r8,5.253823e-01_r8,5.249311e-01_r8,5.246629e-01_r8,5.245946e-01_r8,&
        & 5.247434e-01_r8,5.251268e-01_r8,5.257626e-01_r8,5.266693e-01_r8,5.278653e-01_r8,&
        & 5.293698e-01_r8,5.312022e-01_r8,5.333823e-01_r8,5.359305e-01_r8,5.388676e-01_r8,&
        & 5.422146e-01_r8 /)

! asymmetry factor: unitless
      asyice3(:, 16) = (/ &
! band 16
        & 8.340752e-01_r8,8.435170e-01_r8,8.517487e-01_r8,8.592064e-01_r8,8.660387e-01_r8,&
        & 8.723204e-01_r8,8.780997e-01_r8,8.834137e-01_r8,8.882934e-01_r8,8.927662e-01_r8,&
        & 8.968577e-01_r8,9.005914e-01_r8,9.039899e-01_r8,9.070745e-01_r8,9.098659e-01_r8,&
        & 9.123836e-01_r8,9.146466e-01_r8,9.166734e-01_r8,9.184817e-01_r8,9.200886e-01_r8,&
        & 9.215109e-01_r8,9.227648e-01_r8,9.238661e-01_r8,9.248304e-01_r8,9.256727e-01_r8,&
        & 9.264078e-01_r8,9.270505e-01_r8,9.276150e-01_r8,9.281156e-01_r8,9.285662e-01_r8,&
        & 9.289806e-01_r8,9.293726e-01_r8,9.297557e-01_r8,9.301435e-01_r8,9.305491e-01_r8,&
        & 9.309859e-01_r8,9.314671e-01_r8,9.320055e-01_r8,9.326140e-01_r8,9.333053e-01_r8,&
        & 9.340919e-01_r8,9.349861e-01_r8,9.360000e-01_r8,9.371451e-01_r8,9.384329e-01_r8,&
        & 9.398744e-01_r8 /)
      asyice3(:, 17) = (/ &
! band 17
        & 8.728160e-01_r8,8.777333e-01_r8,8.823754e-01_r8,8.867535e-01_r8,8.908785e-01_r8,&
        & 8.947611e-01_r8,8.984118e-01_r8,9.018408e-01_r8,9.050582e-01_r8,9.080739e-01_r8,&
        & 9.108976e-01_r8,9.135388e-01_r8,9.160068e-01_r8,9.183106e-01_r8,9.204595e-01_r8,&
        & 9.224620e-01_r8,9.243271e-01_r8,9.260632e-01_r8,9.276788e-01_r8,9.291822e-01_r8,&
        & 9.305817e-01_r8,9.318853e-01_r8,9.331012e-01_r8,9.342372e-01_r8,9.353013e-01_r8,&
        & 9.363013e-01_r8,9.372450e-01_r8,9.381400e-01_r8,9.389939e-01_r8,9.398145e-01_r8,&
        & 9.406092e-01_r8,9.413856e-01_r8,9.421511e-01_r8,9.429131e-01_r8,9.436790e-01_r8,&
        & 9.444561e-01_r8,9.452517e-01_r8,9.460729e-01_r8,9.469270e-01_r8,9.478209e-01_r8,&
        & 9.487617e-01_r8,9.497562e-01_r8,9.508112e-01_r8,9.519335e-01_r8,9.531294e-01_r8,&
        & 9.544055e-01_r8 /)
      asyice3(:, 18) = (/ &
! band 18
        & 7.897566e-01_r8,7.948704e-01_r8,7.998041e-01_r8,8.045623e-01_r8,8.091495e-01_r8,&
        & 8.135702e-01_r8,8.178290e-01_r8,8.219305e-01_r8,8.258790e-01_r8,8.296792e-01_r8,&
        & 8.333355e-01_r8,8.368524e-01_r8,8.402343e-01_r8,8.434856e-01_r8,8.466108e-01_r8,&
        & 8.496143e-01_r8,8.525004e-01_r8,8.552737e-01_r8,8.579384e-01_r8,8.604990e-01_r8,&
        & 8.629597e-01_r8,8.653250e-01_r8,8.675992e-01_r8,8.697867e-01_r8,8.718916e-01_r8,&
        & 8.739185e-01_r8,8.758715e-01_r8,8.777551e-01_r8,8.795734e-01_r8,8.813308e-01_r8,&
        & 8.830315e-01_r8,8.846799e-01_r8,8.862802e-01_r8,8.878366e-01_r8,8.893534e-01_r8,&
        & 8.908350e-01_r8,8.922854e-01_r8,8.937090e-01_r8,8.951099e-01_r8,8.964925e-01_r8,&
        & 8.978609e-01_r8,8.992192e-01_r8,9.005718e-01_r8,9.019229e-01_r8,9.032765e-01_r8,&
        & 9.046369e-01_r8 /)
      asyice3(:, 19) = (/ &
! band 19
        & 7.812615e-01_r8,7.887764e-01_r8,7.959664e-01_r8,8.028413e-01_r8,8.094109e-01_r8,&
        & 8.156849e-01_r8,8.216730e-01_r8,8.273846e-01_r8,8.328294e-01_r8,8.380166e-01_r8,&
        & 8.429556e-01_r8,8.476556e-01_r8,8.521258e-01_r8,8.563753e-01_r8,8.604131e-01_r8,&
        & 8.642481e-01_r8,8.678893e-01_r8,8.713455e-01_r8,8.746254e-01_r8,8.777378e-01_r8,&
        & 8.806914e-01_r8,8.834948e-01_r8,8.861566e-01_r8,8.886854e-01_r8,8.910897e-01_r8,&
        & 8.933779e-01_r8,8.955586e-01_r8,8.976402e-01_r8,8.996311e-01_r8,9.015398e-01_r8,&
        & 9.033745e-01_r8,9.051436e-01_r8,9.068555e-01_r8,9.085185e-01_r8,9.101410e-01_r8,&
        & 9.117311e-01_r8,9.132972e-01_r8,9.148476e-01_r8,9.163905e-01_r8,9.179340e-01_r8,&
        & 9.194864e-01_r8,9.210559e-01_r8,9.226505e-01_r8,9.242784e-01_r8,9.259476e-01_r8,&
        & 9.276661e-01_r8 /)
      asyice3(:, 20) = (/ &
! band 20
        & 7.640720e-01_r8,7.691119e-01_r8,7.739941e-01_r8,7.787222e-01_r8,7.832998e-01_r8,&
        & 7.877304e-01_r8,7.920177e-01_r8,7.961652e-01_r8,8.001765e-01_r8,8.040551e-01_r8,&
        & 8.078044e-01_r8,8.114280e-01_r8,8.149294e-01_r8,8.183119e-01_r8,8.215791e-01_r8,&
        & 8.247344e-01_r8,8.277812e-01_r8,8.307229e-01_r8,8.335629e-01_r8,8.363046e-01_r8,&
        & 8.389514e-01_r8,8.415067e-01_r8,8.439738e-01_r8,8.463560e-01_r8,8.486568e-01_r8,&
        & 8.508795e-01_r8,8.530274e-01_r8,8.551039e-01_r8,8.571122e-01_r8,8.590558e-01_r8,&
        & 8.609378e-01_r8,8.627618e-01_r8,8.645309e-01_r8,8.662485e-01_r8,8.679178e-01_r8,&
        & 8.695423e-01_r8,8.711251e-01_r8,8.726697e-01_r8,8.741792e-01_r8,8.756571e-01_r8,&
        & 8.771065e-01_r8,8.785307e-01_r8,8.799331e-01_r8,8.813169e-01_r8,8.826854e-01_r8,&
        & 8.840419e-01_r8 /)
      asyice3(:, 21) = (/ &
! band 21
        & 7.602598e-01_r8,7.651572e-01_r8,7.699014e-01_r8,7.744962e-01_r8,7.789452e-01_r8,&
        & 7.832522e-01_r8,7.874205e-01_r8,7.914538e-01_r8,7.953555e-01_r8,7.991290e-01_r8,&
        & 8.027777e-01_r8,8.063049e-01_r8,8.097140e-01_r8,8.130081e-01_r8,8.161906e-01_r8,&
        & 8.192645e-01_r8,8.222331e-01_r8,8.250993e-01_r8,8.278664e-01_r8,8.305374e-01_r8,&
        & 8.331153e-01_r8,8.356030e-01_r8,8.380037e-01_r8,8.403201e-01_r8,8.425553e-01_r8,&
        & 8.447121e-01_r8,8.467935e-01_r8,8.488022e-01_r8,8.507412e-01_r8,8.526132e-01_r8,&
        & 8.544210e-01_r8,8.561675e-01_r8,8.578554e-01_r8,8.594875e-01_r8,8.610665e-01_r8,&
        & 8.625951e-01_r8,8.640760e-01_r8,8.655119e-01_r8,8.669055e-01_r8,8.682594e-01_r8,&
        & 8.695763e-01_r8,8.708587e-01_r8,8.721094e-01_r8,8.733308e-01_r8,8.745255e-01_r8,&
        & 8.756961e-01_r8 /)
      asyice3(:, 22) = (/ &
! band 22
        & 7.568957e-01_r8,7.606995e-01_r8,7.644072e-01_r8,7.680204e-01_r8,7.715402e-01_r8,&
        & 7.749682e-01_r8,7.783057e-01_r8,7.815541e-01_r8,7.847148e-01_r8,7.877892e-01_r8,&
        & 7.907786e-01_r8,7.936846e-01_r8,7.965084e-01_r8,7.992515e-01_r8,8.019153e-01_r8,&
        & 8.045011e-01_r8,8.070103e-01_r8,8.094444e-01_r8,8.118048e-01_r8,8.140927e-01_r8,&
        & 8.163097e-01_r8,8.184571e-01_r8,8.205364e-01_r8,8.225488e-01_r8,8.244958e-01_r8,&
        & 8.263789e-01_r8,8.281993e-01_r8,8.299586e-01_r8,8.316580e-01_r8,8.332991e-01_r8,&
        & 8.348831e-01_r8,8.364115e-01_r8,8.378857e-01_r8,8.393071e-01_r8,8.406770e-01_r8,&
        & 8.419969e-01_r8,8.432682e-01_r8,8.444923e-01_r8,8.456706e-01_r8,8.468044e-01_r8,&
        & 8.478952e-01_r8,8.489444e-01_r8,8.499533e-01_r8,8.509234e-01_r8,8.518561e-01_r8,&
        & 8.527528e-01_r8 /)
      asyice3(:, 23) = (/ &
! band 23
        & 7.575066e-01_r8,7.606912e-01_r8,7.638236e-01_r8,7.669035e-01_r8,7.699306e-01_r8,&
        & 7.729046e-01_r8,7.758254e-01_r8,7.786926e-01_r8,7.815060e-01_r8,7.842654e-01_r8,&
        & 7.869705e-01_r8,7.896211e-01_r8,7.922168e-01_r8,7.947574e-01_r8,7.972428e-01_r8,&
        & 7.996726e-01_r8,8.020466e-01_r8,8.043646e-01_r8,8.066262e-01_r8,8.088313e-01_r8,&
        & 8.109796e-01_r8,8.130709e-01_r8,8.151049e-01_r8,8.170814e-01_r8,8.190001e-01_r8,&
        & 8.208608e-01_r8,8.226632e-01_r8,8.244071e-01_r8,8.260924e-01_r8,8.277186e-01_r8,&
        & 8.292856e-01_r8,8.307932e-01_r8,8.322411e-01_r8,8.336291e-01_r8,8.349570e-01_r8,&
        & 8.362244e-01_r8,8.374312e-01_r8,8.385772e-01_r8,8.396621e-01_r8,8.406856e-01_r8,&
        & 8.416476e-01_r8,8.425479e-01_r8,8.433861e-01_r8,8.441620e-01_r8,8.448755e-01_r8,&
        & 8.455263e-01_r8 /)
      asyice3(:, 24) = (/ &
! band 24
        & 7.568829e-01_r8,7.597947e-01_r8,7.626745e-01_r8,7.655212e-01_r8,7.683337e-01_r8,&
        & 7.711111e-01_r8,7.738523e-01_r8,7.765565e-01_r8,7.792225e-01_r8,7.818494e-01_r8,&
        & 7.844362e-01_r8,7.869819e-01_r8,7.894854e-01_r8,7.919459e-01_r8,7.943623e-01_r8,&
        & 7.967337e-01_r8,7.990590e-01_r8,8.013373e-01_r8,8.035676e-01_r8,8.057488e-01_r8,&
        & 8.078802e-01_r8,8.099605e-01_r8,8.119890e-01_r8,8.139645e-01_r8,8.158862e-01_r8,&
        & 8.177530e-01_r8,8.195641e-01_r8,8.213183e-01_r8,8.230149e-01_r8,8.246527e-01_r8,&
        & 8.262308e-01_r8,8.277483e-01_r8,8.292042e-01_r8,8.305976e-01_r8,8.319275e-01_r8,&
        & 8.331929e-01_r8,8.343929e-01_r8,8.355265e-01_r8,8.365928e-01_r8,8.375909e-01_r8,&
        & 8.385197e-01_r8,8.393784e-01_r8,8.401659e-01_r8,8.408815e-01_r8,8.415240e-01_r8,&
        & 8.420926e-01_r8 /)
      asyice3(:, 25) = (/ &
! band 25
        & 7.548616e-01_r8,7.575454e-01_r8,7.602153e-01_r8,7.628696e-01_r8,7.655067e-01_r8,&
        & 7.681249e-01_r8,7.707225e-01_r8,7.732978e-01_r8,7.758492e-01_r8,7.783750e-01_r8,&
        & 7.808735e-01_r8,7.833430e-01_r8,7.857819e-01_r8,7.881886e-01_r8,7.905612e-01_r8,&
        & 7.928983e-01_r8,7.951980e-01_r8,7.974588e-01_r8,7.996789e-01_r8,8.018567e-01_r8,&
        & 8.039905e-01_r8,8.060787e-01_r8,8.081196e-01_r8,8.101115e-01_r8,8.120527e-01_r8,&
        & 8.139416e-01_r8,8.157764e-01_r8,8.175557e-01_r8,8.192776e-01_r8,8.209405e-01_r8,&
        & 8.225427e-01_r8,8.240826e-01_r8,8.255585e-01_r8,8.269688e-01_r8,8.283117e-01_r8,&
        & 8.295856e-01_r8,8.307889e-01_r8,8.319198e-01_r8,8.329767e-01_r8,8.339579e-01_r8,&
        & 8.348619e-01_r8,8.356868e-01_r8,8.364311e-01_r8,8.370930e-01_r8,8.376710e-01_r8,&
        & 8.381633e-01_r8 /)
      asyice3(:, 26) = (/ &
! band 26
        & 7.491854e-01_r8,7.518523e-01_r8,7.545089e-01_r8,7.571534e-01_r8,7.597839e-01_r8,&
        & 7.623987e-01_r8,7.649959e-01_r8,7.675737e-01_r8,7.701303e-01_r8,7.726639e-01_r8,&
        & 7.751727e-01_r8,7.776548e-01_r8,7.801084e-01_r8,7.825318e-01_r8,7.849230e-01_r8,&
        & 7.872804e-01_r8,7.896020e-01_r8,7.918862e-01_r8,7.941309e-01_r8,7.963345e-01_r8,&
        & 7.984951e-01_r8,8.006109e-01_r8,8.026802e-01_r8,8.047009e-01_r8,8.066715e-01_r8,&
        & 8.085900e-01_r8,8.104546e-01_r8,8.122636e-01_r8,8.140150e-01_r8,8.157072e-01_r8,&
        & 8.173382e-01_r8,8.189063e-01_r8,8.204096e-01_r8,8.218464e-01_r8,8.232148e-01_r8,&
        & 8.245130e-01_r8,8.257391e-01_r8,8.268915e-01_r8,8.279682e-01_r8,8.289675e-01_r8,&
        & 8.298875e-01_r8,8.307264e-01_r8,8.314824e-01_r8,8.321537e-01_r8,8.327385e-01_r8,&
        & 8.332350e-01_r8 /)
      asyice3(:, 27) = (/ &
! band 27
        & 7.397086e-01_r8,7.424069e-01_r8,7.450955e-01_r8,7.477725e-01_r8,7.504362e-01_r8,&
        & 7.530846e-01_r8,7.557159e-01_r8,7.583283e-01_r8,7.609199e-01_r8,7.634888e-01_r8,&
        & 7.660332e-01_r8,7.685512e-01_r8,7.710411e-01_r8,7.735009e-01_r8,7.759288e-01_r8,&
        & 7.783229e-01_r8,7.806814e-01_r8,7.830024e-01_r8,7.852841e-01_r8,7.875246e-01_r8,&
        & 7.897221e-01_r8,7.918748e-01_r8,7.939807e-01_r8,7.960380e-01_r8,7.980449e-01_r8,&
        & 7.999995e-01_r8,8.019000e-01_r8,8.037445e-01_r8,8.055311e-01_r8,8.072581e-01_r8,&
        & 8.089235e-01_r8,8.105255e-01_r8,8.120623e-01_r8,8.135319e-01_r8,8.149326e-01_r8,&
        & 8.162626e-01_r8,8.175198e-01_r8,8.187025e-01_r8,8.198089e-01_r8,8.208371e-01_r8,&
        & 8.217852e-01_r8,8.226514e-01_r8,8.234338e-01_r8,8.241306e-01_r8,8.247399e-01_r8,&
        & 8.252599e-01_r8 /)
      asyice3(:, 28) = (/ &
! band 28
        & 7.224533e-01_r8,7.251681e-01_r8,7.278728e-01_r8,7.305654e-01_r8,7.332444e-01_r8,&
        & 7.359078e-01_r8,7.385539e-01_r8,7.411808e-01_r8,7.437869e-01_r8,7.463702e-01_r8,&
        & 7.489291e-01_r8,7.514616e-01_r8,7.539661e-01_r8,7.564408e-01_r8,7.588837e-01_r8,&
        & 7.612933e-01_r8,7.636676e-01_r8,7.660049e-01_r8,7.683034e-01_r8,7.705612e-01_r8,&
        & 7.727767e-01_r8,7.749480e-01_r8,7.770733e-01_r8,7.791509e-01_r8,7.811789e-01_r8,&
        & 7.831556e-01_r8,7.850791e-01_r8,7.869478e-01_r8,7.887597e-01_r8,7.905131e-01_r8,&
        & 7.922062e-01_r8,7.938372e-01_r8,7.954044e-01_r8,7.969059e-01_r8,7.983399e-01_r8,&
        & 7.997047e-01_r8,8.009985e-01_r8,8.022195e-01_r8,8.033658e-01_r8,8.044357e-01_r8,&
        & 8.054275e-01_r8,8.063392e-01_r8,8.071692e-01_r8,8.079157e-01_r8,8.085768e-01_r8,&
        & 8.091507e-01_r8 /)
      asyice3(:, 29) = (/ &
! band 29
        & 8.850026e-01_r8,9.005489e-01_r8,9.069242e-01_r8,9.121799e-01_r8,9.168987e-01_r8,&
        & 9.212259e-01_r8,9.252176e-01_r8,9.289028e-01_r8,9.323000e-01_r8,9.354235e-01_r8,&
        & 9.382858e-01_r8,9.408985e-01_r8,9.432734e-01_r8,9.454218e-01_r8,9.473557e-01_r8,&
        & 9.490871e-01_r8,9.506282e-01_r8,9.519917e-01_r8,9.531904e-01_r8,9.542374e-01_r8,&
        & 9.551461e-01_r8,9.559298e-01_r8,9.566023e-01_r8,9.571775e-01_r8,9.576692e-01_r8,&
        & 9.580916e-01_r8,9.584589e-01_r8,9.587853e-01_r8,9.590851e-01_r8,9.593729e-01_r8,&
        & 9.596632e-01_r8,9.599705e-01_r8,9.603096e-01_r8,9.606954e-01_r8,9.611427e-01_r8,&
        & 9.616667e-01_r8,9.622826e-01_r8,9.630060e-01_r8,9.638524e-01_r8,9.648379e-01_r8,&
        & 9.659788e-01_r8,9.672916e-01_r8,9.687933e-01_r8,9.705014e-01_r8,9.724337e-01_r8,&
        & 9.746084e-01_r8 /)

! fdelta: unitless
      fdlice3(:, 16) = (/ &
! band 16
        & 4.959277e-02_r8,4.685292e-02_r8,4.426104e-02_r8,4.181231e-02_r8,3.950191e-02_r8,&
        & 3.732500e-02_r8,3.527675e-02_r8,3.335235e-02_r8,3.154697e-02_r8,2.985578e-02_r8,&
        & 2.827395e-02_r8,2.679666e-02_r8,2.541909e-02_r8,2.413640e-02_r8,2.294378e-02_r8,&
        & 2.183639e-02_r8,2.080940e-02_r8,1.985801e-02_r8,1.897736e-02_r8,1.816265e-02_r8,&
        & 1.740905e-02_r8,1.671172e-02_r8,1.606585e-02_r8,1.546661e-02_r8,1.490917e-02_r8,&
        & 1.438870e-02_r8,1.390038e-02_r8,1.343939e-02_r8,1.300089e-02_r8,1.258006e-02_r8,&
        & 1.217208e-02_r8,1.177212e-02_r8,1.137536e-02_r8,1.097696e-02_r8,1.057210e-02_r8,&
        & 1.015596e-02_r8,9.723704e-03_r8,9.270516e-03_r8,8.791565e-03_r8,8.282026e-03_r8,&
        & 7.737072e-03_r8,7.151879e-03_r8,6.521619e-03_r8,5.841467e-03_r8,5.106597e-03_r8,&
        & 4.312183e-03_r8 /)
      fdlice3(:, 17) = (/ &
! band 17
        & 5.071224e-02_r8,5.000217e-02_r8,4.933872e-02_r8,4.871992e-02_r8,4.814380e-02_r8,&
        & 4.760839e-02_r8,4.711170e-02_r8,4.665177e-02_r8,4.622662e-02_r8,4.583426e-02_r8,&
        & 4.547274e-02_r8,4.514007e-02_r8,4.483428e-02_r8,4.455340e-02_r8,4.429544e-02_r8,&
        & 4.405844e-02_r8,4.384041e-02_r8,4.363939e-02_r8,4.345340e-02_r8,4.328047e-02_r8,&
        & 4.311861e-02_r8,4.296586e-02_r8,4.282024e-02_r8,4.267977e-02_r8,4.254248e-02_r8,&
        & 4.240640e-02_r8,4.226955e-02_r8,4.212995e-02_r8,4.198564e-02_r8,4.183462e-02_r8,&
        & 4.167494e-02_r8,4.150462e-02_r8,4.132167e-02_r8,4.112413e-02_r8,4.091003e-02_r8,&
        & 4.067737e-02_r8,4.042420e-02_r8,4.014854e-02_r8,3.984840e-02_r8,3.952183e-02_r8,&
        & 3.916683e-02_r8,3.878144e-02_r8,3.836368e-02_r8,3.791158e-02_r8,3.742316e-02_r8,&
        & 3.689645e-02_r8 /)
      fdlice3(:, 18) = (/ &
! band 18
        & 1.062938e-01_r8,1.065234e-01_r8,1.067822e-01_r8,1.070682e-01_r8,1.073793e-01_r8,&
        & 1.077137e-01_r8,1.080693e-01_r8,1.084442e-01_r8,1.088364e-01_r8,1.092439e-01_r8,&
        & 1.096647e-01_r8,1.100970e-01_r8,1.105387e-01_r8,1.109878e-01_r8,1.114423e-01_r8,&
        & 1.119004e-01_r8,1.123599e-01_r8,1.128190e-01_r8,1.132757e-01_r8,1.137279e-01_r8,&
        & 1.141738e-01_r8,1.146113e-01_r8,1.150385e-01_r8,1.154534e-01_r8,1.158540e-01_r8,&
        & 1.162383e-01_r8,1.166045e-01_r8,1.169504e-01_r8,1.172741e-01_r8,1.175738e-01_r8,&
        & 1.178472e-01_r8,1.180926e-01_r8,1.183080e-01_r8,1.184913e-01_r8,1.186405e-01_r8,&
        & 1.187538e-01_r8,1.188291e-01_r8,1.188645e-01_r8,1.188580e-01_r8,1.188076e-01_r8,&
        & 1.187113e-01_r8,1.185672e-01_r8,1.183733e-01_r8,1.181277e-01_r8,1.178282e-01_r8,&
        & 1.174731e-01_r8 /)
      fdlice3(:, 19) = (/ &
! band 19
        & 1.076195e-01_r8,1.065195e-01_r8,1.054696e-01_r8,1.044673e-01_r8,1.035099e-01_r8,&
        & 1.025951e-01_r8,1.017203e-01_r8,1.008831e-01_r8,1.000808e-01_r8,9.931116e-02_r8,&
        & 9.857151e-02_r8,9.785939e-02_r8,9.717230e-02_r8,9.650774e-02_r8,9.586322e-02_r8,&
        & 9.523623e-02_r8,9.462427e-02_r8,9.402484e-02_r8,9.343544e-02_r8,9.285358e-02_r8,&
        & 9.227675e-02_r8,9.170245e-02_r8,9.112818e-02_r8,9.055144e-02_r8,8.996974e-02_r8,&
        & 8.938056e-02_r8,8.878142e-02_r8,8.816981e-02_r8,8.754323e-02_r8,8.689919e-02_r8,&
        & 8.623517e-02_r8,8.554869e-02_r8,8.483724e-02_r8,8.409832e-02_r8,8.332943e-02_r8,&
        & 8.252807e-02_r8,8.169175e-02_r8,8.081795e-02_r8,7.990419e-02_r8,7.894796e-02_r8,&
        & 7.794676e-02_r8,7.689809e-02_r8,7.579945e-02_r8,7.464834e-02_r8,7.344227e-02_r8,&
        & 7.217872e-02_r8 /)
      fdlice3(:, 20) = (/ &
! band 20
        & 1.119014e-01_r8,1.122706e-01_r8,1.126690e-01_r8,1.130947e-01_r8,1.135456e-01_r8,&
        & 1.140199e-01_r8,1.145154e-01_r8,1.150302e-01_r8,1.155623e-01_r8,1.161096e-01_r8,&
        & 1.166703e-01_r8,1.172422e-01_r8,1.178233e-01_r8,1.184118e-01_r8,1.190055e-01_r8,&
        & 1.196025e-01_r8,1.202008e-01_r8,1.207983e-01_r8,1.213931e-01_r8,1.219832e-01_r8,&
        & 1.225665e-01_r8,1.231411e-01_r8,1.237050e-01_r8,1.242561e-01_r8,1.247926e-01_r8,&
        & 1.253122e-01_r8,1.258132e-01_r8,1.262934e-01_r8,1.267509e-01_r8,1.271836e-01_r8,&
        & 1.275896e-01_r8,1.279669e-01_r8,1.283134e-01_r8,1.286272e-01_r8,1.289063e-01_r8,&
        & 1.291486e-01_r8,1.293522e-01_r8,1.295150e-01_r8,1.296351e-01_r8,1.297104e-01_r8,&
        & 1.297390e-01_r8,1.297189e-01_r8,1.296480e-01_r8,1.295244e-01_r8,1.293460e-01_r8,&
        & 1.291109e-01_r8 /)
      fdlice3(:, 21) = (/ &
! band 21
        & 1.133298e-01_r8,1.136777e-01_r8,1.140556e-01_r8,1.144615e-01_r8,1.148934e-01_r8,&
        & 1.153492e-01_r8,1.158269e-01_r8,1.163243e-01_r8,1.168396e-01_r8,1.173706e-01_r8,&
        & 1.179152e-01_r8,1.184715e-01_r8,1.190374e-01_r8,1.196108e-01_r8,1.201897e-01_r8,&
        & 1.207720e-01_r8,1.213558e-01_r8,1.219389e-01_r8,1.225194e-01_r8,1.230951e-01_r8,&
        & 1.236640e-01_r8,1.242241e-01_r8,1.247733e-01_r8,1.253096e-01_r8,1.258309e-01_r8,&
        & 1.263352e-01_r8,1.268205e-01_r8,1.272847e-01_r8,1.277257e-01_r8,1.281415e-01_r8,&
        & 1.285300e-01_r8,1.288893e-01_r8,1.292173e-01_r8,1.295118e-01_r8,1.297710e-01_r8,&
        & 1.299927e-01_r8,1.301748e-01_r8,1.303154e-01_r8,1.304124e-01_r8,1.304637e-01_r8,&
        & 1.304673e-01_r8,1.304212e-01_r8,1.303233e-01_r8,1.301715e-01_r8,1.299638e-01_r8,&
        & 1.296983e-01_r8 /)
      fdlice3(:, 22) = (/ &
! band 22
        & 1.145360e-01_r8,1.153256e-01_r8,1.161453e-01_r8,1.169929e-01_r8,1.178666e-01_r8,&
        & 1.187641e-01_r8,1.196835e-01_r8,1.206227e-01_r8,1.215796e-01_r8,1.225522e-01_r8,&
        & 1.235383e-01_r8,1.245361e-01_r8,1.255433e-01_r8,1.265579e-01_r8,1.275779e-01_r8,&
        & 1.286011e-01_r8,1.296257e-01_r8,1.306494e-01_r8,1.316703e-01_r8,1.326862e-01_r8,&
        & 1.336951e-01_r8,1.346950e-01_r8,1.356838e-01_r8,1.366594e-01_r8,1.376198e-01_r8,&
        & 1.385629e-01_r8,1.394866e-01_r8,1.403889e-01_r8,1.412678e-01_r8,1.421212e-01_r8,&
        & 1.429469e-01_r8,1.437430e-01_r8,1.445074e-01_r8,1.452381e-01_r8,1.459329e-01_r8,&
        & 1.465899e-01_r8,1.472069e-01_r8,1.477819e-01_r8,1.483128e-01_r8,1.487976e-01_r8,&
        & 1.492343e-01_r8,1.496207e-01_r8,1.499548e-01_r8,1.502346e-01_r8,1.504579e-01_r8,&
        & 1.506227e-01_r8 /)
      fdlice3(:, 23) = (/ &
! band 23
        & 1.153263e-01_r8,1.161445e-01_r8,1.169932e-01_r8,1.178703e-01_r8,1.187738e-01_r8,&
        & 1.197016e-01_r8,1.206516e-01_r8,1.216217e-01_r8,1.226099e-01_r8,1.236141e-01_r8,&
        & 1.246322e-01_r8,1.256621e-01_r8,1.267017e-01_r8,1.277491e-01_r8,1.288020e-01_r8,&
        & 1.298584e-01_r8,1.309163e-01_r8,1.319736e-01_r8,1.330281e-01_r8,1.340778e-01_r8,&
        & 1.351207e-01_r8,1.361546e-01_r8,1.371775e-01_r8,1.381873e-01_r8,1.391820e-01_r8,&
        & 1.401593e-01_r8,1.411174e-01_r8,1.420540e-01_r8,1.429671e-01_r8,1.438547e-01_r8,&
        & 1.447146e-01_r8,1.455449e-01_r8,1.463433e-01_r8,1.471078e-01_r8,1.478364e-01_r8,&
        & 1.485270e-01_r8,1.491774e-01_r8,1.497857e-01_r8,1.503497e-01_r8,1.508674e-01_r8,&
        & 1.513367e-01_r8,1.517554e-01_r8,1.521216e-01_r8,1.524332e-01_r8,1.526880e-01_r8,&
        & 1.528840e-01_r8 /)
      fdlice3(:, 24) = (/ &
! band 24
        & 1.160842e-01_r8,1.169118e-01_r8,1.177697e-01_r8,1.186556e-01_r8,1.195676e-01_r8,&
        & 1.205036e-01_r8,1.214616e-01_r8,1.224394e-01_r8,1.234349e-01_r8,1.244463e-01_r8,&
        & 1.254712e-01_r8,1.265078e-01_r8,1.275539e-01_r8,1.286075e-01_r8,1.296664e-01_r8,&
        & 1.307287e-01_r8,1.317923e-01_r8,1.328550e-01_r8,1.339149e-01_r8,1.349699e-01_r8,&
        & 1.360179e-01_r8,1.370567e-01_r8,1.380845e-01_r8,1.390991e-01_r8,1.400984e-01_r8,&
        & 1.410803e-01_r8,1.420429e-01_r8,1.429840e-01_r8,1.439016e-01_r8,1.447936e-01_r8,&
        & 1.456579e-01_r8,1.464925e-01_r8,1.472953e-01_r8,1.480642e-01_r8,1.487972e-01_r8,&
        & 1.494923e-01_r8,1.501472e-01_r8,1.507601e-01_r8,1.513287e-01_r8,1.518511e-01_r8,&
        & 1.523252e-01_r8,1.527489e-01_r8,1.531201e-01_r8,1.534368e-01_r8,1.536969e-01_r8,&
        & 1.538984e-01_r8 /)
      fdlice3(:, 25) = (/ &
! band 25
        & 1.168725e-01_r8,1.177088e-01_r8,1.185747e-01_r8,1.194680e-01_r8,1.203867e-01_r8,&
        & 1.213288e-01_r8,1.222923e-01_r8,1.232750e-01_r8,1.242750e-01_r8,1.252903e-01_r8,&
        & 1.263187e-01_r8,1.273583e-01_r8,1.284069e-01_r8,1.294626e-01_r8,1.305233e-01_r8,&
        & 1.315870e-01_r8,1.326517e-01_r8,1.337152e-01_r8,1.347756e-01_r8,1.358308e-01_r8,&
        & 1.368788e-01_r8,1.379175e-01_r8,1.389449e-01_r8,1.399590e-01_r8,1.409577e-01_r8,&
        & 1.419389e-01_r8,1.429007e-01_r8,1.438410e-01_r8,1.447577e-01_r8,1.456488e-01_r8,&
        & 1.465123e-01_r8,1.473461e-01_r8,1.481483e-01_r8,1.489166e-01_r8,1.496492e-01_r8,&
        & 1.503439e-01_r8,1.509988e-01_r8,1.516118e-01_r8,1.521808e-01_r8,1.527038e-01_r8,&
        & 1.531788e-01_r8,1.536037e-01_r8,1.539764e-01_r8,1.542951e-01_r8,1.545575e-01_r8,&
        & 1.547617e-01_r8 /)
      fdlice3(:, 26) = (/ &
!band 26
        & 1.180509e-01_r8,1.189025e-01_r8,1.197820e-01_r8,1.206875e-01_r8,1.216171e-01_r8,&
        & 1.225687e-01_r8,1.235404e-01_r8,1.245303e-01_r8,1.255363e-01_r8,1.265564e-01_r8,&
        & 1.275888e-01_r8,1.286313e-01_r8,1.296821e-01_r8,1.307392e-01_r8,1.318006e-01_r8,&
        & 1.328643e-01_r8,1.339284e-01_r8,1.349908e-01_r8,1.360497e-01_r8,1.371029e-01_r8,&
        & 1.381486e-01_r8,1.391848e-01_r8,1.402095e-01_r8,1.412208e-01_r8,1.422165e-01_r8,&
        & 1.431949e-01_r8,1.441539e-01_r8,1.450915e-01_r8,1.460058e-01_r8,1.468947e-01_r8,&
        & 1.477564e-01_r8,1.485888e-01_r8,1.493900e-01_r8,1.501580e-01_r8,1.508907e-01_r8,&
        & 1.515864e-01_r8,1.522428e-01_r8,1.528582e-01_r8,1.534305e-01_r8,1.539578e-01_r8,&
        & 1.544380e-01_r8,1.548692e-01_r8,1.552494e-01_r8,1.555767e-01_r8,1.558490e-01_r8,&
        & 1.560645e-01_r8 /)
      fdlice3(:, 27) = (/ &
! band 27
        & 1.200480e-01_r8,1.209267e-01_r8,1.218304e-01_r8,1.227575e-01_r8,1.237059e-01_r8,&
        & 1.246739e-01_r8,1.256595e-01_r8,1.266610e-01_r8,1.276765e-01_r8,1.287041e-01_r8,&
        & 1.297420e-01_r8,1.307883e-01_r8,1.318412e-01_r8,1.328988e-01_r8,1.339593e-01_r8,&
        & 1.350207e-01_r8,1.360813e-01_r8,1.371393e-01_r8,1.381926e-01_r8,1.392396e-01_r8,&
        & 1.402783e-01_r8,1.413069e-01_r8,1.423235e-01_r8,1.433263e-01_r8,1.443134e-01_r8,&
        & 1.452830e-01_r8,1.462332e-01_r8,1.471622e-01_r8,1.480681e-01_r8,1.489490e-01_r8,&
        & 1.498032e-01_r8,1.506286e-01_r8,1.514236e-01_r8,1.521863e-01_r8,1.529147e-01_r8,&
        & 1.536070e-01_r8,1.542614e-01_r8,1.548761e-01_r8,1.554491e-01_r8,1.559787e-01_r8,&
        & 1.564629e-01_r8,1.568999e-01_r8,1.572879e-01_r8,1.576249e-01_r8,1.579093e-01_r8,&
        & 1.581390e-01_r8 /)
      fdlice3(:, 28) = (/ &
! band 28
        & 1.247813e-01_r8,1.256496e-01_r8,1.265417e-01_r8,1.274560e-01_r8,1.283905e-01_r8,&
        & 1.293436e-01_r8,1.303135e-01_r8,1.312983e-01_r8,1.322964e-01_r8,1.333060e-01_r8,&
        & 1.343252e-01_r8,1.353523e-01_r8,1.363855e-01_r8,1.374231e-01_r8,1.384632e-01_r8,&
        & 1.395042e-01_r8,1.405441e-01_r8,1.415813e-01_r8,1.426140e-01_r8,1.436404e-01_r8,&
        & 1.446587e-01_r8,1.456672e-01_r8,1.466640e-01_r8,1.476475e-01_r8,1.486157e-01_r8,&
        & 1.495671e-01_r8,1.504997e-01_r8,1.514117e-01_r8,1.523016e-01_r8,1.531673e-01_r8,&
        & 1.540073e-01_r8,1.548197e-01_r8,1.556026e-01_r8,1.563545e-01_r8,1.570734e-01_r8,&
        & 1.577576e-01_r8,1.584054e-01_r8,1.590149e-01_r8,1.595843e-01_r8,1.601120e-01_r8,&
        & 1.605962e-01_r8,1.610349e-01_r8,1.614266e-01_r8,1.617693e-01_r8,1.620614e-01_r8,&
        & 1.623011e-01_r8 /)
      fdlice3(:, 29) = (/ &
! band 29
        & 1.006055e-01_r8,9.549582e-02_r8,9.063960e-02_r8,8.602900e-02_r8,8.165612e-02_r8,&
        & 7.751308e-02_r8,7.359199e-02_r8,6.988496e-02_r8,6.638412e-02_r8,6.308156e-02_r8,&
        & 5.996942e-02_r8,5.703979e-02_r8,5.428481e-02_r8,5.169657e-02_r8,4.926719e-02_r8,&
        & 4.698880e-02_r8,4.485349e-02_r8,4.285339e-02_r8,4.098061e-02_r8,3.922727e-02_r8,&
        & 3.758547e-02_r8,3.604733e-02_r8,3.460497e-02_r8,3.325051e-02_r8,3.197604e-02_r8,&
        & 3.077369e-02_r8,2.963558e-02_r8,2.855381e-02_r8,2.752050e-02_r8,2.652776e-02_r8,&
        & 2.556772e-02_r8,2.463247e-02_r8,2.371415e-02_r8,2.280485e-02_r8,2.189670e-02_r8,&
        & 2.098180e-02_r8,2.005228e-02_r8,1.910024e-02_r8,1.811781e-02_r8,1.709709e-02_r8,&
        & 1.603020e-02_r8,1.490925e-02_r8,1.372635e-02_r8,1.247363e-02_r8,1.114319e-02_r8,&
        & 9.727157e-03_r8 /)


      wavenum2(:) = (/3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, 8050._r8, &
                     12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8, 2600._r8/)

      end subroutine swcldpr

END MODULE Rad_Clirad
