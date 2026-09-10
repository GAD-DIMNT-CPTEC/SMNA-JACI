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
       rei   , rel   , taud              )

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

    ! local  variables
    INTEGER i,k
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

    REAL(KIND=r8):: tauk    (ncols,kmax+1) ! total optical depth with k=1 at top
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
    nsol  =0.0_r8
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
    DO i=1,nsol
      IF(litx(i) .le. nCols) THEN
         co2(i,1) = co2l(litx(i),1)
      END IF
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

       ! Prepare tau, fice, rel e rei
       DO k=1,kmax
          DO i=1,nsol
             IF(litx(I).le.ncols)THEN
             tauk  (i,k+1) = taud(litx(i),k)
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
    !CO2=co2val*1.0E-6_r8

    CALL cloudy( &
         s0                 , CO2                ,   nsol               ,   kmax+1,   &
         pu(1:nsol,1:kmax+2), ta(1:nsol,1:kmax+1), wa(1:nsol,1:kmax+1)  , oa(1:nsol,1:kmax+1), &
         cmu(1:nsol),  rsurfv(1:nsol), agv(1:nsol), rsurfn(1:nsol), agn(1:nsol),   &
         dscld(1:nsol), sc(1:nsol),  acld(1:nsol,1:kmax+1), rvbc(1:nsol), rvdc(1:nsol), rnbc(1:nsol), rndc(1:nsol),  &
         dsclr(1:nsol), sl(1:nsol),  aclr(1:nsol,1:kmax+1), rvbl(1:nsol), rvdl(1:nsol), rnbl(1:nsol), rndl(1:nsol),  &
         e0(1:nsol,1:kmax+1), css(1:nsol,1:kmax+1),  ps(1:nsol),   &
         dmask(1:nsol), frcice(1:nsol,1:kmax+1), radliq(1:nsol,1:kmax+1), radice(1:nsol,1:kmax+1),icld)
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
  SUBROUTINE cloudy(s0,rco2,m,np,&
       pl,ta,wa,oa, &
       cosz,rsuvbm, rsuvdf,rsirbm,rsirdf,&
       dscld1,sc1,acld1,rvbc1,rvdc1,rnbc1,rndc1, &
       dsclr1,sl1,aclr1,rvbl1,rvdl1,rnbl1,rndl1, &
       tauc,csscgp, psc, &
       dmask,fice, rel, rei, icld)

    IMPLICIT NONE

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

    ! input variables
    INTEGER, INTENT(IN) :: m, np
    INTEGER, INTENT(IN) :: dmask (m) !sib-mask in DLGP
    INTEGER, INTENT(IN) :: icld ! new cloud microphysics

    REAL(KIND=r8), INTENT(IN) :: s0, rco2(m,np)
    ! pl,ta,wa and oa are, respectively, the level pressure (mb), layer
    ! temperature (k), layer specific humidity (g/g), and layer ozone
    ! concentration (g/g)
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np+1) :: pl
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: ta, wa, oa, tauc, csscgp
    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: cosz,rsuvbm,rsuvdf,rsirbm,rsirdf
    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: psc

    REAL(KIND=r8), intent(in), DIMENSION(m,np) :: fice, rel, rei
    ! output variables
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m) :: dscld1,sc1, rvbc1,rvdc1,rnbc1,rndc1
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: acld1

    REAL(KIND=r8), INTENT(OUT), DIMENSION(m) :: dsclr1,sl1, rvbl1,rvdl1,rnbl1,rndl1
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m,np) :: aclr1

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
          END DO
       END DO
    END DO

    ! compute solar fluxes

    CALL soradcld (m,np,pl,ta,wa,oa,rco2,  &
         taucld,reff,fcld,ict,icb,  &
         taual,ssaal,asyal,  &
         cosz,rsuvbm,rsuvdf,rsirbm,rsirdf,  &
         flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,fdirir,fdifir,  &
         fdiruv_c,fdifuv_c,fdirpar_c,fdifpar_c,fdirir_c,fdifir_c,  &
         flx_d,flx_u,flc_d,flc_u)    ! new

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
  SUBROUTINE Soradcld (m,np,pl,ta,wa,oa,co2,  &
       taucld,reff,fcld,ict,icb,  &
       taual,ssaal,asyal,  &
       cosz,rsuvbm,rsuvdf,rsirbm,rsirdf,  &
       flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,fdirir,fdifir,  &
       fdiruv_c,fdifuv_c,fdirpar_c,fdifpar_c,fdirir_c,fdifir_c,  &
       flx_d,flx_u,flc_d,flc_u)

    IMPLICIT NONE
    INTEGER i,j


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
    INTEGER, INTENT(IN) :: m,np
    INTEGER, INTENT(IN) :: ict(m),icb(m)
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: taual, ssaal, asyal
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3) :: taucld, reff
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np+1) :: pl
    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np) :: ta, wa, oa, fcld
    REAL(KIND=r8), INTENT(IN), DIMENSION(m) :: cosz, rsuvbm, rsuvdf, rsirbm, rsirdf
    REAL(KIND=r8), INTENT(IN) :: co2(m,np)

    ! output variables
    REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx,flc,flx_d,flx_u, flc_d,flc_u
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m) :: fdiruv,fdifuv,fdirpar,fdifpar,fdirir,fdifir
    REAL(KIND=r8), INTENT(OUT), DIMENSION(m) :: fdiruv_c,fdifuv_c,fdirpar_c,fdifpar_c,fdirir_c,fdifir_c

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
         m         , &!    INTEGER, INTENT(IN) :: m
         np        , &!    INTEGER, INTENT(IN) :: np
         wh        (1:m,1:np), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)	 :: wh
         oh        (1:m,1:np), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)	 :: oh
         dp        (1:m,1:np), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)	 :: dp
         taucld    (1:m,1:np,1:3), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: taucld
         reff      (1:m,1:np,1:3), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,3)  :: reff
         ict       (1:m), &!    INTEGER, INTENT(IN) :: ict(m)
         icb       (1:m), &!    INTEGER, INTENT(IN) :: icb(m)
         fcld      (1:m,1:np), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np)	 :: fcld
         cosz      (1:m), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m)	 :: cosz
         taual     (1:m,1:np,1:11), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: taual
         ssaal     (1:m,1:np,1:11), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: ssaal
         asyal     (1:m,1:np,1:11), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m,np,11) :: asyal
         rsuvbm    (1:m), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m)	 :: rsuvbm
         rsuvdf    (1:m), &!    REAL(KIND=r8), INTENT(IN), DIMENSION(m)	 :: rsuvdf
         flx       (1:m,1:np+1), &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx
         flc       (1:m,1:np+1), &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flc
         fdiruv    (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdiruv
         fdifuv    (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifuv
         fdirpar   (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdirpar
         fdifpar   (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifpar
         fdiruv_c  (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdiruv_c
         fdifuv_c  (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifuv_c
         fdirpar_c (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdirpar_c
         fdifpar_c (1:m), &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)  :: fdifpar_c
         flx_d     (1:m,1:np+1), &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_d
         flx_u     (1:m,1:np+1), &!REAL(KIND=r8), INTENT(INOUT), DIMENSION(m,np+1) :: flx_u
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

    DO k=2,np+1
       DO i=1,m
          swu(i,k)=swu(i,k)
       ENDDO
    ENDDO

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
          IF (fcld(i,k).gt.0.02_r8 .and. nctop(i).eq.np+1) THEN
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
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,2) :: aib = RESHAPE( &
         SHAPE = (/ nband, 2 /), SOURCE = (/ &
         0.000333_r8, 0.000333_r8, 0.000333_r8, 2.52_r8,    2.52_r8,2.52_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,2) :: awb = RESHAPE( &
         SHAPE = (/ nband, 2 /), SOURCE = (/ &
         -0.0101_r8, -0.0166_r8, -0.0339_r8, 1.72_r8,    1.85_r8,2.16_r8 /) )
    REAL(KIND=r8), PARAMETER, DIMENSION(nband,2) :: arb = RESHAPE( &
         SHAPE = (/ nband, 2 /), SOURCE = (/ &
         0.00307_r8, 0.00307_r8, 0.00307_r8, 0.0_r8    , 0.0_r8    , 0.0_r8  /) )

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
       fdiruv   , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdiruv
       fdifuv   , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdifuv
       fdirpar  , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdirpar
       fdifpar  , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdifpar
       fdiruv_c , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdiruv_c
       fdifuv_c , &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdifuv_c
       fdirpar_c, &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdirpar_c
       fdifpar_c, &!REAL(KIND=r8), INTENT(OUT  ), DIMENSION(m)	:: fdifpar_c
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

    INTEGER i,j
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


END MODULE Rad_Clirad
