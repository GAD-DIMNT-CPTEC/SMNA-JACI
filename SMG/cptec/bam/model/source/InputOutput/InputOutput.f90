!
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.15 $
!
!  Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
!  Modifications (10)are marked by 
!  Fine mode 2000 is also included
!  !tar begin and  !tar end
!  For using with RRTMG and CRD, CRDTF
!
MODULE InputOutput

  USE Parallelism, ONLY: &
       myId,             &
       myId_four,        &
       maxNodes      ,   &
       FatalError

  USE Constants, ONLY: &
       r4, i4, r8, i8, root2, ngrmx, numx, ncf, ncf2

  USE IOLowLevel, ONLY: &
       ReadHead     , &
       GReadHead    , &
       ReadField    , &
       GReadField   , &
       WriteHead    , &
       GWriteHead   , &
       WriteField   , &
       GWriteField  , &
       ReadGetALB   , &
       ReadGetSST   , &
       ReadGetSLM   , &
       ReadGetSLM3D , &
       ReadGetSNW   , &
       ReadGetSST2  , &
       ReadOzone    , &
       ReadTracer   , &
!tar begin 
!Climate aerosol parameter reading
       Read_Aeros  
!tar end


  USE Options, ONLY: &
       nfprt, &
       nfctrl, &
       nfsst,nfco2fx, &
       nfndvi, &
       nfSoilMostSib2,&
       nfsnw, &
       nfalb, &
       nfslm, &
       nfauntbl, &
       nfcnftbl, &
       nfcnf2tb, &
       nflooktb, &
       reducedGrid,&
       labelsi,&
       labelsj,&
       labelsi_soilm,&
       labelsj_soilm,&
       labelsi_flxco2, &
       labelsj_flxco2, &
       ifco2, ifozone, & !hmjb
       nfco2, nfozone, & !hmjb
       nftrc,iftracer, &
       co2val,StrFormat, &            !hmjb for new co2 values
       fNameSnow  , &
       fNameSSTAOI, &
       fNameCO2FLX, &
       fNameNDVIAOI,&
       fNameSoilms, &
       fNameSoilmsWkl, &
       fNameSoilMoistSib2,&
       fNameAlbedo, &
       fNameCO2   , &
       fNameOzone , &
       fNametracer, &
       schemes    , &
       ifndvi     , &
       ifslmSib2  , &
       intndvi    , &
       ifalb      , &
       ifsst      , &
       ifco2flx   , &
       isimco2    , &
       ifslm      , &
       ifsnw      , &
       ifozone    , &

!tar begin
!Climate aerosol file names and file numbers of coarse mode
      
       ifaeros    , &
       fNameClimAodRRTM, &
       fNameClimAsyRRTM, &       
       fNameClimSsaRRTM, & 
       fNameClimAodVrt, &
       nfaod, &
       nfasy, &
       nfssa, &
       nfaodvrt, &             
!tar end
!
!tar begin
!Climate aerosol file names and file numbers of fine mode 
      
       fNameFineAodRRTM, &
       fNameFineAsyRRTM, &       
       fNameFineSsaRRTM, & 
       fNameFineAodVrt, &
       nfaodF, &
       nfasyF, &
       nfssaF, &
       nfaodvrtF, &             
!tar end
       SetBCCte   , &
       Flxco2lag  , &
       sstlag     , &
       soilmlag   , &
       intsst     , &
       intsoilm   , &
       intflxco2  , &
       fint       , &
       yrl        , &
       monl       , &
       jull

  USE Utils, ONLY: &
       IJtoIBJB      ,&
       AveBoxIJtoIBJB,&
       NearestIJtoIBJB

  USE Sizes, ONLY: &
       mymnmax,    &
       myjmax_d,   &
       imax,       &
       jmax,       &
       kmaxloc,    &
       myfirstlev, &
       mylastlev,  &
       HaveM1,     &
       mymmax,     &
       msinproc,   &
       mnmap,      &
       mymnmap


  USE Communications, ONLY: &
       Collect_Spec,        &
       Collect_Grid_Sur,    &
       Collect_Grid_d

   IMPLICIT NONE
  SAVE       


  PRIVATE
  PUBLIC :: InitInputOutput
  PUBLIC :: UpDateGetsbc
  PUBLIC :: cnvray
  PUBLIC :: scloutsp
  PUBLIC :: scloutgr
  PUBLIC :: WillGetSbc
  PUBLIC :: getsbc
  PUBLIC :: gread
  PUBLIC :: gread4
  PUBLIC :: gwrite
  PUBLIC :: fsbc
  PUBLIC :: aunits

  INTEGER              :: mMax
  INTEGER              :: nMax
  INTEGER              :: mnMax
  INTEGER              :: kMax

  LOGICAL              :: fsbc
  CHARACTER(LEN=100)   :: path

  CHARACTER(LEN=16), ALLOCATABLE :: aunits(:)
  INTEGER,           ALLOCATABLE :: looku (:,:,:)
  REAL(KIND=r8),              ALLOCATABLE :: cnfac (:)
  REAL(KIND=r8),              ALLOCATABLE :: cnfac2(:)

  REAL(KIND=r8),    PARAMETER   :: undef =1.0e53_r8


CONTAINS



  ! InitInputOutput: Initializes module

  SUBROUTINE InitInputOutput ( &
       mMax_in, nMax_in, mnMax_in, kmax_in, &
       path_in, fNameCnfTbl, &
       fNameCnf2Tb, fNameLookTb, fNameUnitTb)

    INTEGER,          INTENT(IN) :: mMax_in
    INTEGER,          INTENT(IN) :: nMax_in
    INTEGER,          INTENT(IN) :: mnMax_in
    INTEGER,          INTENT(IN) :: kmax_in
    CHARACTER(LEN=*), INTENT(IN) :: path_in
    CHARACTER(LEN=*), INTENT(IN) :: fNameCnfTbl
    CHARACTER(LEN=*), INTENT(IN) :: fNameCnf2Tb
    CHARACTER(LEN=*), INTENT(IN) :: fNameLookTb
    CHARACTER(LEN=*), INTENT(IN) :: fNameUnitTb
    INTEGER :: ierr


    OPEN(UNIT=nfauntbl, FILE=TRIM(fNameUnitTb), FORM='formatted',ACCESS='sequential',&
         ACTION='read', STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameUnitTb), ierr
       STOP "**(ERROR)**"
    END IF

    OPEN(UNIT=nfcnftbl, FILE=TRIM(fNameCnfTbl), FORM='formatted',ACCESS='sequential',&
         ACTION='read', STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameCnfTbl), ierr
       STOP "**(ERROR)**"
    END IF

    OPEN(UNIT=nfcnf2tb, FILE=TRIM(fNameCnf2Tb), FORM='formatted',ACCESS='sequential',&
         ACTION='read', STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameCnf2Tb), ierr
       STOP "**(ERROR)**"
    END IF

    OPEN(UNIT=nflooktb, FILE=TRIM(fNameLookTb), FORM='formatted',ACCESS='sequential',&
         ACTION='read', STATUS='old', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameLookTb), ierr
       STOP "**(ERROR)**"
    END IF

    path  = path_in
    mMax  = mMax_in
    nMax  = nMax_in
    mnMax = mnMax_in
    kmax  = kmax_in

    ALLOCATE(aunits(-1:numx))
    REWIND nfauntbl
    READ(UNIT=nfauntbl,FMT="(A16)") aunits
    REWIND nfauntbl
    CLOSE(UNIT=nfauntbl,status='KEEP')

    ALLOCATE(cnfac(ncf))
    REWIND nfcnftbl
    READ(UNIT=nfcnftbl,FMT="(5e16.8)") cnfac
    REWIND nfcnftbl
    CLOSE(UNIT=nfcnftbl,status='KEEP')

    ALLOCATE(cnfac2(ncf2))
    REWIND nfcnf2tb
    READ(UNIT=nfcnf2tb,FMT="(5e16.8)") cnfac2
    REWIND nfcnf2tb
    CLOSE(UNIT=nfcnf2tb,status='KEEP')

    ALLOCATE(looku(0:9,0:9,0:ngrmx))
    REWIND nflooktb
    READ(UNIT=nflooktb,FMT="(20i4)") looku
    REWIND nflooktb
    CLOSE(UNIT=nflooktb,status='KEEP')

  END SUBROUTINE InitInputOutput



  ! cnvray: convert array



  SUBROUTINE cnvray (array, idim, ifr, ito)
    INTEGER, INTENT(IN   ) :: idim
    REAL(KIND=r8),    INTENT(INOUT) :: array(idim)
    INTEGER, INTENT(IN   ) :: ifr
    INTEGER, INTENT(IN   ) :: ito

    CHARACTER(LEN=20) :: c0
    CHARACTER(LEN=20) :: c1
    INTEGER           :: i
    INTEGER           :: icf
    INTEGER           :: igpf
    INTEGER           :: iuf
    INTEGER           :: igpt
    INTEGER           :: iut
    REAL(KIND=r8)              :: cf
    REAL(KIND=r8)              :: cf2

    CHARACTER(LEN=*), PARAMETER :: h="**(cnvray)**"

    ! consistency

    IF (idim.eq.0) RETURN

    IF (ifr <= -1) THEN
       WRITE(c0,"(i20)") ifr
       WRITE(UNIT=nfprt,FMT="(a)") h//" ERROR: ifr ("//TRIM(ADJUSTL(c0))//") <= -1 "
       STOP h
    ELSE IF (ito <= -1) THEN
       WRITE(c0,"(i20)") ito
       WRITE(UNIT=nfprt,FMT="(a)") h//" ERROR: ito ("//TRIM(ADJUSTL(c0))//") <= -1 "
       STOP h
    ELSE IF (idim < 0) THEN
       WRITE(c0,"(i20)") idim
       WRITE(UNIT=nfprt,FMT="(a)") h//" ERROR: idim ("//TRIM(ADJUSTL(c0))//") < 0 "
       STOP h
    ELSE IF (ito /= ifr) THEN
       igpf=ifr/10
       igpt=ito/10
       IF (igpf /= igpt) THEN
          WRITE(c0,"(i20)") igpf
          WRITE(c1,"(i20)") igpt
          WRITE(UNIT=nfprt,FMT="(a)") h//" ERROR: igpf ("//TRIM(ADJUSTL(c0))//&
               &") /= igpt ("//TRIM(ADJUSTL(c1))//")"
          STOP h
       ELSE IF (igpf > ngrmx) THEN
          WRITE(c0,"(i20)") igpf
          WRITE(c1,"(i20)") ngrmx
          WRITE(UNIT=nfprt,FMT="(a)") h//" ERROR: igpf ("//TRIM(ADJUSTL(c0))//&
               &") > ngrmx ("//TRIM(ADJUSTL(c1))//")"
          STOP h
       ELSE

          ! table look-up

          iuf=MOD(ifr,10)
          iut=MOD(ito,10)
          icf=looku(iuf,iut,igpf)

          ! consistency, again

          IF (icf < 1 .OR. icf > ncf) THEN
             WRITE(c0,"(i20)") icf
             WRITE(c1,"(i20)") ncf
             WRITE(UNIT=nfprt,FMT="(a)") h//" ERROR: icf ("//TRIM(ADJUSTL(c0))//&
                  &") < 1 or > ncf ("//TRIM(ADJUSTL(c1))//")"
             STOP h
          END IF

          ! get coeficients

          cf=cnfac(icf)
          IF (icf <= ncf2) THEN
             cf2=cnfac2(icf)
          ELSE
             cf2=0.0_r8
          END IF

          ! convert array

          DO i = 1, idim
             IF (array(i) /= undef) THEN
                array(i)=cf*array(i)+cf2
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE cnvray


  !
  ! scale, convert to 32 bits and output field
  !
  SUBROUTINE scloutsp(unit, field, levs, levsg, fact1, nufr, nuto)
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: levs
    INTEGER, INTENT(IN) :: levsg
    REAL(KIND=r8),    INTENT(IN) :: field(2*mymnmax,levs)
    REAL(KIND=r8),    INTENT(IN) :: fact1
    INTEGER, INTENT(IN) :: nufr
    INTEGER, INTENT(IN) :: nuto
    REAL(KIND=r8) :: fldaux(2*mymnmax,levs), fout(2*mnmax,levsg)

    fldaux = fact1 * field

    CALL cnvray(fldaux,2*levs*mymnmax,nufr,nuto)

    IF (Maxnodes.eq.1) THEN

       ! inversion from top to bottom to bottom to top
       CALL WriteField(unit, fldaux(:,levsg:1:-1))

    ELSE
       CALL Collect_Spec(fldaux, fout, levs, levsg, 0)
       ! inversion from top to bottom to bottom to top
       IF (myid.eq.0) CALL WriteField(unit, fout(:,levsg:1:-1))
    ENDIF



  END SUBROUTINE scloutsp




  SUBROUTINE scloutgr(unit, field, levs, fact1, nufr, nuto)
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: levs
    REAL(KIND=r8),    INTENT(IN) :: field(imax*myjmax_d,levs)
    REAL(KIND=r8),    INTENT(IN) :: fact1
    INTEGER, INTENT(IN) :: nufr
    INTEGER, INTENT(IN) :: nuto
    REAL(KIND=r8) :: fldaux(imax*myjmax_d,levs), fout(imax*jmax,levs)
    fldaux=0.0_r8
    fout=0.0_r8
    fldaux = fact1 * field

    CALL cnvray(fldaux,imax*myjmax_d*levs,nufr,nuto)

    IF (Maxnodes.eq.1) THEN

       CALL WriteField(unit, fldaux)

    ELSE
       IF(levs.eq.1) THEN
          CALL Collect_Grid_Sur(fldaux, fout, 0)
       ELSE
          CALL Collect_Grid_d(fldaux, fout, levs, 0)
       ENDIF
       IF(myid.eq.0) CALL WriteField(unit, fout)

    ENDIF



  END SUBROUTINE scloutgr




  LOGICAL FUNCTION WillGetSbc(idate, tod, fint)
    INTEGER, INTENT(IN) :: idate(4)
    REAL(KIND=r8),    INTENT(IN) :: tod
    REAL(KIND=r8),    INTENT(IN) :: fint
    REAL(KIND=r8)                :: fhr

    WillGetSbc = .TRUE.
    IF (fint > 0.0_r8) THEN
       fhr=REAL(idate(1),r8)+tod/3600.0_r8+1.0e-3_r8
       WillGetSbc = fsbc .OR. ABS( MOD(fhr,fint)) <= 1.0e-2_r8
    END IF
  END FUNCTION WillGetSbc
 
   SUBROUTINE UpDateGetsbc(ifday ,tod   ,idate ,idatec,ibMax,jbMax,kMax,ibMaxPerJB,&
                        AlbVisDiff,gtsea,gco2flx,gndvi,soilm,wsib3d,sheleg,o3mix,tracermix,&
!tar begin
!Climate aerosol parameters of coarse and fine modes
                        aod,asy,ssa,z_aer,aodF,asyF,ssaF,z_aerF)
!tar end

    !
    ! getsbc :read surface boundary conditions.
    !
    INTEGER      , INTENT(in   ) :: ifday
    REAL(KIND=r8), INTENT(in   ) :: tod
    INTEGER      , INTENT(in   ) :: idate     (4)
    INTEGER      , INTENT(in   ) :: idatec    (4)
    INTEGER      , INTENT(in   ) :: ibMax
    INTEGER      , INTENT(in   ) :: jbMax
    INTEGER      , INTENT(in   ) :: kMax
    INTEGER      , INTENT(in   ) :: ibMaxPerJB(jbMax)
    REAL(KIND=r8), INTENT(inout) :: AlbVisDiff(ibMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: gtsea     (ibMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: gco2flx   (ibMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: soilm     (ibMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: gndvi     (ibMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: sheleg    (ibMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: o3mix     (ibMax,kMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: tracermix (ibMax,kMax,jbMax)
    REAL(KIND=r8), INTENT(inout) :: wsib3d    (ibMax,kMax,8)
!tar begin  
!climate aerosol parameters of coarse mode
    REAL(KIND=r8), INTENT(inout) :: z_aer(ibMax,jbMax,40)
    REAL(KIND=r8), INTENT(inout) :: aod(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(inout) :: asy(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(inout) :: ssa(ibMax,jbMax,14)
!tar end 

!tar begin  
!climate aerosol parameters of fine mode
    REAL(KIND=r8), INTENT(inout) :: z_aerF(ibMax,jbMax,40)
    REAL(KIND=r8), INTENT(inout) :: aodF(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(inout) :: asyF(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(inout) :: ssaF(ibMax,jbMax,14)
!tar end  

    REAL(KIND=r8)                   :: tice  =271.16e0_r8

    IF (WillGetSbc(idate, tod, fint)) THEN

      CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,&
                   sheleg,o3mix,tracermix,wsib3d,&
!tar begin 
!climate aerosol parameters of coarse mode
               aod, asy, ssa, z_aer, ifaeros, & 
!tar end
!
!tar begin 
!climate aerosol parameters of fine mode
               aodF, asyF, ssaF, z_aerF, & 
!tar end
               ifday , tod  ,idate ,idatec, &
               ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
               sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
               yrl  ,monl,ibMax,jbMax,ibMaxPerJB)

    END IF

  END SUBROUTINE UpDateGetsbc


  !
  ! getsbc :read surface/atmosphere boundary conditions.
  !
  SUBROUTINE getsbc (imax,jmax,kmax,galb ,gsst ,gco2flx,gndvi,gslm,&
       gsnw,gozo,tracermix,wsib3d,&
!tar begin 
!climate aerosol parameters of coarse mode     
       aod, asy, ssa, z_aer, ifaeros, &  
!tar end 
!
!tar begin 
!climate aerosol parameters of fine mode     
       aodF, asyF, ssaF, z_aerF, &  
!tar end                    
       ifday,tod,idate,idatec,&
       ifalb,ifsst,ifco2flx,ifndvi,ifslm,ifslmSib2,ifsnw,ifozone,iftracer,&
       sstlag,intsst,intndvi,intsoilm,fint,tice,&
       yrl ,monl,ibMax,jbMax,ibMaxPerJB)
    IMPLICIT NONE
    !
    ! INPUT/OUTPUT VARIABLES
    !
    ! Real size of the grid
    INTEGER, INTENT(in   ) :: imax
    INTEGER, INTENT(in   ) :: jmax
    INTEGER, INTENT(in   ) :: kmax
    ! Size of block divided grid
    INTEGER, INTENT(in   ) :: ibMax
    INTEGER, INTENT(in   ) :: jbMax
    INTEGER, INTENT(in   ) :: ibMaxPerJB(jbMax)

    ! Boundary fields output
    REAL(KIND=r8), INTENT(out  ) :: galb(ibMax,jbMax) ! albedo
    REAL(KIND=r8), INTENT(out  ) :: gndvi(ibMax,jbMax) ! ndvi    
    REAL(KIND=r8), INTENT(out  ) :: gsst(ibMax,jbMax) ! sst
    REAL(KIND=r8), INTENT(out  ) :: gco2flx(ibMax,jbMax) ! sst
    REAL(KIND=r8), INTENT(out  ) :: gslm(ibMax,jbMax) ! soil moisture
    REAL(KIND=r8), INTENT(out  ) :: gsnw(ibMax,jbMax) ! snow
    REAL(KIND=r8), INTENT(out  ) :: wsib3d(ibMax,jbMax,8) ! moisture
!tar begin
!climate aerosol parameters of coarse mode
    REAL(KIND=r8), INTENT(out  ) :: z_aer(ibMax,jbMax,40)
    REAL(KIND=r8), INTENT(out  ) :: aod(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(out  ) :: asy(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(out  ) :: ssa(ibMax,jbMax,14)
!tar end 
!
!tar begin
!climate aerosol parameters of fine mode
    REAL(KIND=r8), INTENT(out  ) :: z_aerF(ibMax,jbMax,40)
    REAL(KIND=r8), INTENT(out  ) :: aodF(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(out  ) :: asyF(ibMax,jbMax,14)
    REAL(KIND=r8), INTENT(out  ) :: ssaF(ibMax,jbMax,14)
!tar end 

    !hmjb o ozonio nao pode ser apenas 'out' pois, no caso de usar a antiga
    !  getoz(), ele sairia daqui com valores indefinidos... Com inout,
    !  ele entra e,  se nao for alterado, sai como entrou
    REAL(KIND=r8), INTENT(inout) :: gozo(ibMax,kMax,jbMax) ! ozone
    REAL(KIND=r8), INTENT(inout) :: tracermix(ibMax,kMax,jbMax) ! tracer

    ! Options for reading boundary fields
    INTEGER, INTENT(inout) :: ifalb
    INTEGER, INTENT(inout) :: ifsst
    INTEGER, INTENT(inout) :: ifco2flx
    INTEGER, INTENT(inout) :: ifndvi
    INTEGER, INTENT(inout) :: ifslm
    INTEGER, INTENT(inout) :: ifsnw
    INTEGER, INTENT(inout) :: ifslmSib2
    INTEGER, INTENT(inout) :: ifozone
!tar begin   
    INTEGER, INTENT(inout) :: ifaeros
!tar end
    INTEGER, INTENT(inout) :: iftracer
    ! Time
    INTEGER, INTENT(in   ) :: ifday
    REAL(KIND=r8), INTENT(in   ) :: tod
    INTEGER, INTENT(in   ) :: idate(4)
    INTEGER, INTENT(in   ) :: idatec(4)
    REAL(KIND=r8), INTENT(in   ) :: sstlag
    INTEGER, INTENT(in   ) :: intsst
    INTEGER, INTENT(in   ) :: intndvi
    INTEGER, INTENT(IN   ) :: intsoilm
    REAL(KIND=r8), INTENT(in   ) :: fint
    REAL(KIND=r8), INTENT(in   ) :: tice
    REAL(KIND=r8), INTENT(in   ) :: yrl
    INTEGER, INTENT(in   ) :: monl(12)
    !
    ! LOCAL VARIABLES
    !
    REAL(KIND=r8)                :: xfco2   (ibMax,jbMax)
    REAL(KIND=r8)                :: xndvi   (ibMax,jbMax)
    REAL(KIND=r8)                :: xsoilm   (ibMax,jbMax)
    REAL(KIND=r8)                :: xsst    (ibMax,jbMax)
    REAL(KIND=r8)                :: bfr_in  (imax,jmax)
    REAL(KIND=r8)                :: bfrw_in  (imax,jmax,3)
    REAL(KIND=r8)                :: bfrw_out  (ibmax,jbmax,3)
    REAL(KIND=r4)                :: rbrfw3d    (iMax,jMax,3)
!tar begin 
!for reading climate aerosol parameters
    Real(KIND=r8)                :: abfr14     (iMax,jMax,14) 
    Real(KIND=r8)                :: abfr40     (iMax,jMax,40)
    Real(KIND=r8)                :: aers14_wrk (ibMax,jbMax,14)
    Real(KIND=r8)                :: aers40_wrk (ibMax,jbMax,40)
!tar end 

    REAL(KIND=r8)                :: bfr_in3 (imax,kmax,jmax)
    REAL(KIND=r8)                :: bfr_out (ibMax,jbMax)
    REAL(KIND=r8)                :: bfr_out3(ibMax,kmax,jbMax)
    REAL(KIND=r4)                :: rbrf    (iMax,jMax)
    REAL(KIND=r4)                :: rbrf3   (iMax,kmax,jMax)

    !
    !
    INTEGER                :: lrecl,LRecIn
    REAL(KIND=r8)          :: fhr
    INTEGER                :: mf
    INTEGER                :: mn
    INTEGER                :: mf_ndvi
    INTEGER                :: mn_ndvi
    INTEGER                :: mf_co2flx
    INTEGER                :: mn_co2flx
    INTEGER                :: mf_aer
    INTEGER                :: mn_aer
    
    INTEGER                :: irec_co2flx
    INTEGER                :: irec_aer
    INTEGER                :: irec_soilm
    INTEGER                :: irec2_soilm
    REAL(KIND=r8)          :: f1_soilm
    REAL(KIND=r8)          :: f2_soilm
    INTEGER                :: mf_soilm
    INTEGER                :: mn_soilm
    INTEGER                :: month
    INTEGER                :: mm
    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: k
    INTEGER                :: km
    INTEGER                :: irec
    INTEGER                :: irec2
    INTEGER                :: irec_ndvi

    REAL(KIND=r8)                :: f1
    REAL(KIND=r8)                :: f2
    REAL(KIND=r8)                :: f1_ndvi
    REAL(KIND=r8)                :: f2_ndvi
    REAL(KIND=r8)                :: f1_co2flx
    REAL(KIND=r8)                :: f2_co2flx
    REAL(KIND=r8)                :: f1_aer
    REAL(KIND=r8)                :: f2_aer

    REAL(KIND=r8)                :: gmax
    REAL(KIND=r8)                :: gmin
    REAL(KIND=r8)                :: fsst
    REAL(KIND=r8)                :: fsoilm
    REAL(KIND=r8)                :: fxco2
    REAL(KIND=r8)                :: fndvi
    REAL(KIND=r8)                :: fisst
    REAL(KIND=r8)                :: fico2
    REAL(KIND=r8)                :: fisoilm
    REAL(KIND=r8)                :: findvi
    REAL(KIND=r8)                :: xx1
    REAL(KIND=r8)                :: xx2
    REAL(KIND=r8)                :: xday
    INTEGER :: ierr
    !
    !   ifxxx=0    xxx is not processed
    !   ifxxx=1    xxx is set to month=idate(2) in the first call,
    !              but not processed from the subsequent calls.
    !              ifxxx is set to zero after interpolation
    !   ifxxx=2    xxx is interpolated to current day and time every fint
    !              hours synchronized to 00z regardless of initial time.
    !              interpolation is continuous (every time step) if fint<0.
    !   ifxxx=3    xxx is interpolated to current day and time when ifday=0
    !              and tod=0.0 but not processed otherwise
    !              ( appropriate only when xxx is predicted )
    !
    !              the following are for sst only (fint applies as in
    !              ifxxx=2):
    !   ifsst=4    sst is linearly interpolated from continuous direct
    !              access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days or every
    !              calendar month is intsst < 0.
    !   ifsst=5    sst is expanded from piecewise cubic coefficients in
    !              direct access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days.
    !   note:      for ifsst=4 or 5 sstlag must be set.  sstlag is the
    !              number of days plus any fraction prior to the initial
    !              condition date and time the data set begins if intsst > 0
    !              sstlag is the number of whole months prior to the initial
    !              condition date the data set begins if intsst < 0.
    !
    !     ifsst=-1 for numerical weather forecasting using mean weekly sst:
    !              sst is read in the first call as the second record of
    !              the archieve but not processed from the subsequent calls.
    !
    rbrf=0.0_r4
    rbrfw3d=0.0_r4
    rbrf3=0.0_r4
    IF (fint > 0.0_r8) THEN
       fhr=REAL(idate(1),r8)+tod/3600.0_r8+1.0e-3_r8
       IF (.NOT. fsbc .AND. ABS( MOD(fhr,fint)) > 1.0e-2_r8) THEN
          RETURN
       END IF
    END IF
    IF (ifsst == 4 .AND. intsst <= 0) THEN
       CALL GetRecWgtMonthlySST &
            (idate, idatec, tod, labelsi, labelsj, &
            irec, f1, f2, mf, mn,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    ELSE IF (ifsst == 5 .AND. intsst <= 0) THEN
       CALL GetRecWgtDailySST &
            (idate, idatec, tod, labelsi, labelsj, &
            irec,irec2, f1, f2, mf, mn,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2

    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1, f2,mf)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ifslm == 4 .AND. intsoilm <= 0) THEN
       CALL GetRecWgtMonthlySST &
            (idate, idatec, tod, labelsi_soilm, labelsj_soilm, &
            irec_soilm, f1_soilm, f2_soilm, mf_soilm, mn_soilm,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    ELSE IF (ifslm == 5 .AND. intsoilm <= 0) THEN
       CALL GetRecWgtDailySST &
            (idate, idatec, tod, labelsi_soilm, labelsj_soilm, &
            irec_soilm,irec2_soilm, f1_soilm, f2_soilm, mf_soilm, mn_soilm,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2

    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1_soilm, f2_soilm,mf_soilm)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ifco2flx == 4 .AND. intflxco2 <= 0) THEN
       CALL GetRecWgtMonthlySST &
            (idate, idatec, tod, labelsi_flxco2, labelsj_flxco2, &
            irec_co2flx, f1_co2flx, f2_co2flx, mf_co2flx, mn_co2flx,monl)

!!$        WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlyCO2'
!!$        WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$             ' reci = ', irec_co2flx, ' recf = ', irec_co2flx+1, &
!!$             ' mra = ', mf_co2flx, ' mrb = ', mf_co2flx+1
!!$        WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1_co2flx, ' fb  (*mrb) = ', f2_co2flx
    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1_co2flx, f2_co2flx,mf_co2flx)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ifndvi == 4 .AND. intndvi <= 0) THEN
       CALL GetRecWgtMonthlySST &
            (idate, idatec, tod, labelsi, labelsj, &
            irec_ndvi, f1_ndvi, f2_ndvi, mf_ndvi, mn_ndvi,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1_ndvi, f2_ndvi,mf_ndvi)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    IF (ifaeros == 4 .AND. ifaeros <= 0) THEN
!       CALL GetRecWgtMonthlySST &
!            (idate, idatec, tod, labelsi, labelsj, &
!            irec_caer, f1_aer, f2_aer, mf_aer, mn_aer,monl)

!!$       WRITE (UNIT=nfprt, FMT='(A)') ' GetRecWgtMonthlySST'
!!$       WRITE (UNIT=nfprt, FMT='(/,4(A,I5),/)') &
!!$            ' reci = ', irec, ' recf = ', irec+1, &
!!$            ' mra = ', mf, ' mrb = ', mf+1
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
!    ELSE
       CALL GetWeightsOld(yrl,monl,idatec, tod, f1_aer, f2_aer,mf_aer)
!!$       WRITE (UNIT=nfprt, FMT=*) ' fa  (*mra) = ', f1, ' fb  (*mrb) = ', f2
!    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! process albedo file
    !
    IF (ifalb /= 0) THEN
       IF (ifalb == 1) THEN
          month=idate(2)
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfalb,FILE=TRIM(fNameAlbedo),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameAlbedo), ierr
             STOP "**(ERROR)**"
          END IF
          irec=month
          CALL ReadGetALB(nfalb,irec,bfr_in)

          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,galb)
          ELSE
             CALL IJtoIBJB(bfr_in ,galb)
          END IF
          CLOSE(UNIT=nfalb)
          ifalb=0
       ELSE IF (&
            (ifalb == 2) .OR. &
            (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0)) THEN
          rbrf=0.0_r4
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfalb,FILE=TRIM(fNameAlbedo),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameAlbedo), ierr
             STOP "**(ERROR)**"
          END IF
          irec=mf
          CALL ReadGetALB(nfalb,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,galb)
          ELSE
             CALL IJtoIBJB(bfr_in ,galb)
          END IF
          IF (irec == 12) THEN
             irec=1
          ELSE   
             irec=irec+1
          END IF
          CALL ReadGetALB(nfalb,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF
          CLOSE(UNIT=nfalb)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                galb(i,j)=f2*galb(i,j)+f1*bfr_out(i,j)
                gmax=MAX(gmax,galb(i,j))
                gmin=MIN(gmin,galb(i,j))
             END DO
          END DO

          IF (ifalb == 3 .AND. tod == 0.0_r8 .AND. ifday == 0) THEN
             ifalb=0
          END IF

          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=888) mf,f1,f2,gmax,gmin
          END IF

       ELSE
          WRITE(UNIT=nfprt,FMT=999)
          STOP
       END IF
    END IF

    !
    ! process CO2 Flux file
    !
    IF (ifco2flx /= 0 .and. isimco2 /= 0 ) THEN
       !   ifxxx=0    xxx is not processed
       IF (ifco2flx == -1) THEN
         !     ifsst=-1 for numerical weather forecasting using mean weekly sst:
         !               sst is read in the first call as the second record of
         !              the archieve but not processed from the subsequent calls.
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfco2fx, FILE=TRIM(fNameCO2FLX),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameCO2FLX), ierr
             STOP "**(ERROR)**"
          END IF
          irec_co2flx=1
          CALL ReadGetSST(nfco2fx,irec_co2flx,bfr_in)
          irec_co2flx=2
          CALL ReadGetSST(nfco2fx,irec_co2flx,bfr_in)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,gco2flx)
          ELSE
             CALL IJtoIBJB(bfr_in ,gco2flx)
          END IF
          CLOSE(UNIT=nfco2fx)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (abs(gco2flx(i,j)) > 10.0_r8) THEN
                   gco2flx(i,j)=0.0_r8
                ELSE IF (ABS(gco2flx(i,j)) < 10.0_r8) THEN
                   gco2flx(i,j)=gco2flx(i,j)
                ELSE
                   PRINT *, " OPTION ifco2flx=-1 INCORRECT VALUE OF CO2 "
                   STOP "**(ERROR)**"
                END IF
                gmax=MAX(gmax,gco2flx(i,j))
                gmin=MIN(gmin,gco2flx(i,j))
             END DO
          END DO

          WRITE(UNIT=nfprt,FMT=667) ifco2flx,gmax,gmin
          ifco2flx=0

       ELSE IF (ifco2flx == 1) THEN
    !
    !   ifxxx=1    xxx is set to month=idate(2) in the first call,
    !              but not processed from the subsequent calls.
    !              ifxxx is set to zero after interpolation
          OPEN(UNIT=nfco2fx, FILE=TRIM(fNameCO2FLX), FORM='unformatted', ACCESS='sequential',&
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameCO2FLX), ierr
             STOP "**(ERROR)**"
          END IF
          READ(UNIT=nfco2fx)
          month=idate(2)
          DO mm=1,month
          irec_co2flx=mm
             CALL ReadGetSST(nfco2fx,irec_co2flx,bfr_in)
          END DO
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,gco2flx)
          ELSE
             CALL IJtoIBJB(bfr_in ,gco2flx)
          END IF
          CLOSE(UNIT=nfco2fx)
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (abs(gco2flx(i,j)) > 10.0_r8) THEN
                   gco2flx(i,j)=0.0_r8
                ELSE IF (ABS(gco2flx(i,j)) < 10.0_r8) THEN
                   gco2flx(i,j)=gco2flx(i,j)
                ELSE
                   PRINT *, " OPTION ifco2flx=-1 INCORRECT VALUE OF CO2 "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO
          ifco2flx=0
       ELSE IF (ifco2flx == 2.OR.(ifco2flx == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
    !   ifxxx=2    xxx is interpolated to current day and time every fint
    !              hours synchronized to 00z regardless of initial time.
    !              interpolation is continuous (every time step) if fint<0.
    !   ifxxx=3    xxx is interpolated to current day and time when ifday=0
    !              and tod=0.0 but not processed otherwise
    !              ( appropriate only when xxx is predicted )
    !
    !              the following are for sst only (fint applies as in
    !              ifxxx=2):
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfco2fx, FILE=TRIM(fNameCO2FLX),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
             ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameCO2FLX), ierr
             STOP "**(ERROR)**"
          END IF

          irec_co2flx = mf_co2flx+1
          CALL ReadGetSST(nfco2fx,irec_co2flx,bfr_in)

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,xsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,xsst)
          END IF
          
          IF (irec_co2flx == 13) THEN
             irec_co2flx=2
          ELSE
             irec_co2flx=irec_co2flx+1
          END IF
          CALL ReadGetSST(nfco2fx,irec_co2flx,bfr_in)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF

          CLOSE(UNIT=nfco2fx)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)

                fsst=f2_co2flx*xsst(i,j)+f1_co2flx*bfr_out(i,j)
                IF (fsst > gmax) THEN
                   gmax=fsst
                END IF
                IF (fsst < gmin) THEN
                   gmin=fsst
                END IF
 
                IF (abs(fsst) > 10.0_r8) THEN
                   xsst(i,j)=0.0
                ELSE IF (ABS(fsst) < 10.0_r8) THEN
                   xsst(i,j)=fsst
                ELSE
                   PRINT *, " OPTION ifco2flx=-1 INCORRECT VALUE OF CO2 "
                   STOP "**(ERROR)**"
                END IF

             END DO
          END DO
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
               ! IF (tod == 0.0_r8.AND.ifday == 0) THEN
                   gco2flx(i,j)=xsst(i,j)
               ! END IF
             END DO
          END DO

          IF (ifco2flx == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifco2flx=0
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=666) mf,f1,f2,gmax,gmin
          END IF
          xsst=0.0_r8;bfr_out=0.0_r8;fsst=0.0_r8;irec_co2flx=0
       ELSE IF (ifco2flx == 4) THEN
    !   ifsst=4    sst is linearly interpolated from continuous direct
    !              access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days or every
    !              calendar month is intsst < 0.
          IF (intflxco2 > 0) THEN
             fico2=REAL(intflxco2,r8)
             xday=ifday+tod/86400.0_r8+Flxco2lag
             irec_co2flx=INT(xday/fico2+1.0e-3_r8+1.0_r8)
             xx1= MOD(xday,fico2)/fico2
             xx2=1.0_r8-xx1
          ELSE
             xx1=f1_co2flx
             xx2=f2_co2flx
          END IF
          INQUIRE (IOLENGTH=lrecl) bfr_in
          lrecl=lrecl/2
          OPEN(UNIT=nfco2fx,FILE=TRIM(fNameCO2FLX),FORM='unformatted',ACCESS='direct',&
                  RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameCO2FLX), ierr
                STOP "**(ERROR)**"
          END IF

          CALL ReadGetSST2(nfco2fx,bfr_in,irec_co2flx)

          IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,xfco2)
          ELSE
                CALL IJtoIBJB(bfr_in ,xfco2)
          END IF

          CALL ReadGetSST2(nfco2fx,bfr_in,irec_co2flx+1)

          IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,bfr_out)
          ELSE
                CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF
          CLOSE(UNIT=nfco2fx)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)

                fxco2=xx2*xfco2(i,j)+xx1*bfr_out(i,j)
                IF (fxco2 > gmax) THEN
                   gmax=fxco2
                END IF
                IF (fxco2 < gmin) THEN
                      gmin=fxco2
                END IF

                IF (abs(fxco2) > 10.0_r8) THEN
                   xfco2(i,j)=0.0_r8
                ELSE IF (abs(fxco2) <= 10.0_r8) THEN
                   xfco2(i,j)=fxco2
                ELSE
                   PRINT *, " OPTION ifco2flx=4 INCORRECT VALUE OF CO2 "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO

          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=666) irec_co2flx,xx1,xx2,gmax,gmin
          END IF

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (tod == 0.0_r8.AND.ifday == 0) THEN
                   gco2flx(i,j)=xfco2(i,j)
                ELSE 
                   gco2flx(i,j)=xfco2(i,j)
                END IF
             END DO
          END DO

!          WRITE(UNIT=nfprt,FMT=1999)
!          STOP

       ELSE IF (ifco2flx == 5) THEN
    !   ifsst=5    sst is expanded from piecewise cubic coefficients in
    !              direct access data set to current day and time.  data set
    !              is assumed to be spaced every intsst days.
          WRITE(UNIT=nfprt,FMT=1999)
          STOP

       ELSE
    !   note:      for ifsst=4 or 5 sstlag must be set.  sstlag is the
    !              number of days plus any fraction prior to the initial
    !              condition date and time the data set begins if intsst > 0
    !              sstlag is the number of whole months prior to the initial
    !              condition date the data set begins if intsst < 0.
          WRITE(UNIT=nfprt,FMT=1999)
          STOP
       END IF

    END IF

    !
    ! Aerosol  file
    !

!tar begin
!Reading of climate aerosol parameters (Kinne, 2013) (coarse mode)
!----------------------------------------------------------------   
    IF (ifaeros /= 0) THEN
       IF (ifaeros == 2) THEN
       abfr14=0.0_r8
!
       INQUIRE (IOLENGTH=LRecIn) abfr14
!
!     AOD  reading    
       OPEN (UNIT=nfaod,FILE=TRIM(fNameClimAodRRTM),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='read', STATUS='old', IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                 TRIM(fNameClimAodRRTM), ierr
          STOP "**(ERROR)**"
       END IF
!
!       IF(myid.EQ.0) then
!     write(*,*) 'nfaod=',nfaod, 'fNameClimAodRRTM=', fNameClimAodRRTM
!   END IF    
              
       irec_aer=mf_aer
       CALL Read_Aeros(nfaod,irec_aer,abfr14)
       DO k=1,14
          IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),aod(:,:,k))
             CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aod(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aod(1:ibMax,1:jbMax,k))
          END IF
       END DO
!
       IF (irec_aer == 12) THEN
          irec_aer=1
       ELSE   
          irec_aer=irec_aer+1
       END IF
       CALL Read_Aeros(nfaod,irec_aer,abfr14)
!
       DO k=1,14
          IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),aers14_wrk(:,:,k))
             CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
          END IF
       END DO
       CLOSE(UNIT=nfaod)
!
       DO k=1,14
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                aod(i,j,k)=f2_aer*aod(i,j,k)+f1_aer*aers14_wrk(i,j,k)
             END DO
          END DO
       END DO
!
! ASY reading
       OPEN (UNIT=nfasy,FILE=TRIM(fNameClimAsyRRTM),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='read', STATUS='old', IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameClimAsyRRTM), ierr
          STOP "**(ERROR)**"
       END IF
!
       irec_aer=mf_aer
       CALL Read_Aeros(nfasy,irec_aer,abfr14)
       DO k=1,14
          IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),asy(:,:,k))
             CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),asy(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),asy(1:ibMax,1:jbMax,k))
          END IF
       END DO
!
       IF (irec_aer == 12) THEN
          irec_aer=1
       ELSE   
          irec_aer=irec_aer+1
       END IF
       CALL Read_Aeros(nfasy,irec_aer,abfr14)
!
       DO k=1,14
          IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),aers14_wrk(:,:,k))
             CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
          END IF
       END DO
       CLOSE(UNIT=nfasy)
!
!tar print 1 in InputOutput(getsbc)
!!     IF(myid.EQ.95) then
!!         IF(idatec(1).EQ.6.AND.idatec(3).EQ.1) then
!    
!!     OPEN(unit=75,file='/scratchin/grupos/mcga/home/t.tarassova/OUTPUT_T/Tar.txt', &
!!     ACCESS='APPEND', STATUS='OLD')
!
!!     WRITE(75,*) 'in InputOutput,  getsbc' 
!!      WRITE(75,*) 'myid=', myid   
!     WRITE(75,*) 'asy(1:192,1,3(0.512))=', asy(:,1,3)
!!     WRITE(75,*) 'f2_aer=', f2_aer, 'f1_aer=', f1_aer 
!!     WRITE(75,*) 'idatec(4)=h,m,d,y', idatec(1),idatec(2),idatec(3),idatec(4)         
!!     CLOSE(75) 
!!          ENDIF 
!!    ENDIF
!       
 
       DO k=1,14
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                asy(i,j,k)=f2_aer*asy(i,j,k)+f1_aer*aers14_wrk(i,j,k)
             END DO
          END DO
       END DO


!  SSA reading
       OPEN (UNIT=nfssa,FILE=TRIM(fNameClimSsaRRTM),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
                ACTION='read', STATUS='old', IOSTAT=ierr)
       IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                   TRIM(fNameClimSsaRRTM), ierr
           STOP "**(ERROR)**"
       END IF
!
       irec_aer=mf_aer
       CALL Read_Aeros(nfssa,irec_aer,abfr14)
       DO k=1,14
          IF (reducedGrid) THEN
!            CALL AveBoxIJtoIBJB(abfr14(:,:,k),ssa(:,:,k))
             CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),ssa(1:ibMax,1:jbMax,k))
          ELSE
               CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),ssa(1:ibMax,1:jbMax,k))
            END IF
       END DO
!
       IF (irec_aer == 12) THEN
           irec_aer=1
       ELSE   
           irec_aer=irec_aer+1
       END IF
       CALL Read_Aeros(nfssa,irec_aer,abfr14)
!
       DO k=1,14
          IF (reducedGrid) THEN
!                CALL AveBoxIJtoIBJB(abfr14(:,:,k),aers14_wrk(:,:,k))
             CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
          END IF
       END DO
       CLOSE(UNIT=nfssa)
!
       DO k=1,14
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                ssa(i,j,k)=f2_aer*ssa(i,j,k)+f1_aer*aers14_wrk(i,j,k)
             END DO
          END DO
       END DO
!
!  reading of vertical AOD distribution (AODVRT)
       abfr40 =0.0_r8
       INQUIRE (IOLENGTH=LRecIn) abfr40
       OPEN (UNIT=nfaodvrt,FILE=TRIM(fNameClimAodVrt),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='read', STATUS='old', IOSTAT=ierr)
       IF (ierr /= 0) THEN
            WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameClimAodVrt), ierr
            STOP "**(ERROR)**"
       END IF
!
       irec_aer=mf_aer
       CALL Read_Aeros(nfaodvrt,irec_aer,abfr40)
       DO k=1,40
          IF (reducedGrid) THEN
!            CALL AveBoxIJtoIBJB(abfr40(:,:,k),z_aer(:,:,k))
             CALL NearestIJtoIBJB(abfr40(1:iMax,1:jMax,k),z_aer(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr40(1:iMax,1:jMax,k),z_aer(1:ibMax,1:jbMax,k))
          END IF
       END DO
! 
       IF (irec_aer == 12) THEN
             irec_aer=1
       ELSE   
             irec_aer=irec_aer+1
       END IF
       CALL Read_Aeros(nfaodvrt,irec_aer,abfr40)
! 
       DO k=1,40
          IF (reducedGrid) THEN
!            CALL AveBoxIJtoIBJB(abfr40(:,:,k),aers40_wrk(:,:,k))
             CALL NearestIJtoIBJB(abfr40(1:iMax,1:jMax,k),aers40_wrk(1:ibMax,1:jbMax,k))
          ELSE
             CALL IJtoIBJB(abfr40(1:iMax,1:jMax,k),aers40_wrk(1:ibMax,1:jbMax,k))
          END IF
       END DO
       CLOSE(UNIT=nfaodvrt)
!
       DO k=1,40
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                z_aer (i,j,k)=f2_aer*z_aer(i,j,k)+f1_aer*aers40_wrk(i,j,k)
             END DO
          END DO
       END DO

       END IF
    END IF

!tar print 2
!
!!       IF(myid.EQ.95) then
!!           IF(idatec(1).EQ.6.AND.idatec(3).EQ.1) then
!
!!     OPEN(unit=75,file='/scratchin/grupos/mcga/home/t.tarassova/OUTPUT_T/Tar.txt', &
!!     ACCESS='APPEND', STATUS='OLD')
!
!!     WRITE(75,*) 'in InputOutput,  getsbc' 
!!     WRITE(75,*) 'myid=', myid 
!!     WRITE(75,*)  'ibMaxPerJB(1)=', ibMaxPerJB(1) 
!             
!!     WRITE(75,*) 'time mean aod(1:192,1,3(0.512))=', aod(:,1,3)
!!     WRITE(75,*) 'time mean asy(1:192,1,3(0.512))=', asy(:,1,3)
!!     WRITE(75,*) 'time mean ssa(1:192,1,3(0.512))=', ssa(:,1,3)
     
!     WRITE(75,*) 'time mean aod(1:192,1,1(0.252))=', aod(:,1,1)
!     WRITE(75,*) 'time mean asy(1:192,1,1(0.252))=', asy(:,1,1)
!     WRITE(75,*) 'time mean ssa(1:192,1,1(0.252))=', ssa(:,1,1)

!     WRITE(75,*) 'time mean aod(1:192,1,8(6.135))=', aod(:,1,8)
!     WRITE(75,*) 'time mean asy(1:192,1,8(6.135))=', asy(:,1,8)
!     WRITE(75,*) 'time mean ssa(1:192,1,8(6.135))=', ssa(:,1,8)
!     

!!     WRITE(75,*) 'z_aer(1,1,1:40)=', z_aer(1,1,:)
!!     WRITE(75,*) 'z_aer(70,1,1:40)=', z_aer(70,1,:)
!!     WRITE(75,*) 'z_aer(ibMaxPerJB(1),1,1:40)=', z_aer(ibMaxPerJB(1),1,:)
            
!!     CLOSE(75)
!!           ENDIF
!!       ENDIF


!----------------------------------------------------------------------
!tar end
!
!tar begin
!Reading of climate aerosol parameters (Kinne, 2013) (fine mode)
!----------------------------------------------------------------
    IF (ifaeros /= 0) THEN
       IF (ifaeros == 2) THEN
          abfr14=0.0_r8
!
          INQUIRE (IOLENGTH=LRecIn) abfr14
!
!   AOD  reading
!
          OPEN (UNIT=nfaodF,FILE=TRIM(fNameFineAodRRTM),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                 TRIM(fNameFineAodRRTM), ierr
             STOP "**(ERROR)**"
          END IF
!
!         IF(myid.EQ.0) then
!         write(*,*) 'nfaodF=',nfaodF, 'fNameFineAodRRTM=', fNameFineAodRRTM
!         END IF

          irec_aer=mf_aer
          CALL Read_Aeros(nfaodF,irec_aer,abfr14)
          DO k=1,14
             IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),aodF(:,:,k))
                CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aodF(1:ibMax,1:jbMax,k))
             ELSE
                CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aodF(1:ibMax,1:jbMax,k))
             END IF
          END DO
!
          IF (irec_aer == 12) THEN
             irec_aer=1
          ELSE   
             irec_aer=irec_aer+1
          END IF
          CALL Read_Aeros(nfaodF,irec_aer,abfr14)
!
          DO k=1,14
             IF (reducedGrid) THEN
!                CALL AveBoxIJtoIBJB(abfr14(:,:,k),aers14_wrk(:,:,k))
                CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
             ELSE
                CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
             END IF
          END DO
          CLOSE(UNIT=nfaodF)
!
          DO k=1,14
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   aodF(i,j,k)=f2_aer*aodF(i,j,k)+f1_aer*aers14_wrk(i,j,k)
                END DO
             END DO
          END DO
!
! ASY reading
!
          OPEN (UNIT=nfasyF,FILE=TRIM(fNameFineAsyRRTM),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
                ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
              WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameFineAsyRRTM), ierr
              STOP "**(ERROR)**"
          END IF
!
          irec_aer=mf_aer
          CALL Read_Aeros(nfasyF,irec_aer,abfr14)
          DO k=1,14
             IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),asyF(:,:,k))
                CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),asyF(1:ibMax,1:jbMax,k))
             ELSE
                CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),asyF(1:ibMax,1:jbMax,k))
             END IF
          END DO
!
          IF (irec_aer == 12) THEN
              irec_aer=1
          ELSE   
              irec_aer=irec_aer+1
          END IF
          CALL Read_Aeros(nfasyF,irec_aer,abfr14)
!
          DO k=1,14
             IF (reducedGrid) THEN
!                CALL AveBoxIJtoIBJB(abfr14(:,:,k),aers14_wrk(:,:,k))
                 CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
             ELSE
                 CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
             END IF
          END DO
          CLOSE(UNIT=nfasyF)
!
!tar print 1 in InputOutput(getsbc)
!!     IF(myid.EQ.95) then
!!         IF(idatec(1).EQ.6.AND.idatec(3).EQ.1) then
!    
!!     OPEN(unit=75,file='/scratchin/grupos/mcga/home/t.tarassova/OUTPUT_T/Tar.txt', &
!!     ACCESS='APPEND', STATUS='OLD')
!
!!     WRITE(75,*) 'in InputOutput,  getsbc' 
!!      WRITE(75,*) 'myid=', myid   
!     WRITE(75,*) 'asy(1:192,1,3(0.512))=', asy(:,1,3)
!!     WRITE(75,*) 'f2_aer=', f2_aer, 'f1_aer=', f1_aer 
!!     WRITE(75,*) 'idatec(4)=h,m,d,y', idatec(1),idatec(2),idatec(3),idatec(4)         
!!     CLOSE(75) 
!!          ENDIF 
!!    ENDIF
!

          DO k=1,14
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   asyF(i,j,k)=f2_aer*asyF(i,j,k)+f1_aer*aers14_wrk(i,j,k)
                END DO
             END DO
          END DO


!  SSA reading
          OPEN (UNIT=nfssaF,FILE=TRIM(fNameFineSsaRRTM),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
              WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameFineSsaRRTM), ierr
              STOP "**(ERROR)**"
          END IF
!
         irec_aer=mf_aer
         CALL Read_Aeros(nfssaF,irec_aer,abfr14)
         DO k=1,14
            IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),ssaF(:,:,k))
                CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),ssaF(1:ibMax,1:jbMax,k))
            ELSE
                CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),ssaF(1:ibMax,1:jbMax,k))
            END IF
         END DO
!
         IF (irec_aer == 12) THEN
             irec_aer=1
         ELSE   
             irec_aer=irec_aer+1
         END IF
         CALL Read_Aeros(nfssaF,irec_aer,abfr14)
!
         DO k=1,14
            IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr14(:,:,k),aers14_wrk(:,:,k))
               CALL NearestIJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
            ELSE
               CALL IJtoIBJB(abfr14(1:iMax,1:jMax,k),aers14_wrk(1:ibMax,1:jbMax,k))
           END IF
         ENDDO
         CLOSE(UNIT=nfssaF)
!
         DO k=1,14
            DO j=1,jbMax
               DO i=1,ibMaxPerJB(j)
                  ssaF(i,j,k)=f2_aer*ssaF(i,j,k)+f1_aer*aers14_wrk(i,j,k)
               END DO
            END DO
         END DO
!
!  reading of vertical AOD distribution (AODVRT)
!
         abfr40 =0.0_r8
         INQUIRE (IOLENGTH=LRecIn) abfr40
         OPEN (UNIT=nfaodvrtF,FILE=TRIM(fNameFineAodVrt),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='read', STATUS='old', IOSTAT=ierr)
         IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameFineAodVrt), ierr
             STOP "**(ERROR)**"
         END IF
!
         irec_aer=mf_aer
         CALL Read_Aeros(nfaodvrtF,irec_aer,abfr40)
         DO k=1,40
            IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr40(:,:,k),z_aerF(:,:,k))
               CALL NearestIJtoIBJB(abfr40(1:iMax,1:jMax,k),z_aerF(1:ibMax,1:jbMax,k))
            ELSE
               CALL IJtoIBJB(abfr40(1:iMax,1:jMax,k),z_aerF(1:ibMax,1:jbMax,k))
           END IF
         END DO
!
         IF (irec_aer == 12) THEN
             irec_aer=1
         ELSE   
             irec_aer=irec_aer+1
         END IF
         CALL Read_Aeros(nfaodvrtF,irec_aer,abfr40)
!
         DO k=1,40
            IF (reducedGrid) THEN
!               CALL AveBoxIJtoIBJB(abfr40(:,:,k),aers40_wrk(:,:,k))
                CALL NearestIJtoIBJB(abfr40(1:iMax,1:jMax,k),aers40_wrk(1:ibMax,1:jbMax,k))
            ELSE
               CALL IJtoIBJB(abfr40(1:iMax,1:jMax,k),aers40_wrk(1:ibMax,1:jbMax,k))
            END IF
         END DO
         CLOSE(UNIT=nfaodvrtF)
!
         DO k=1,40
            DO j=1,jbMax
               DO i=1,ibMaxPerJB(j)
                  z_aerF (i,j,k)=f2_aer*z_aerF(i,j,k)+f1_aer*aers40_wrk(i,j,k)
               END DO
            END DO
         END DO
      END IF
   END IF

!tar print 2
!
!!       IF(myid.EQ.95) then
!!           IF(idatec(1).EQ.6.AND.idatec(3).EQ.1) then
!   
!!     OPEN(unit=75,file='/scratchin/grupos/mcga/home/t.tarassova/OUTPUT_T/Tar.txt', &
!!     ACCESS='APPEND', STATUS='OLD')
!
!!     WRITE(75,*) 'in InputOutput,  getsbc' 
!!     WRITE(75,*) 'myid=', myid 
!!     WRITE(75,*)  'ibMaxPerJB(1)=', ibMaxPerJB(1) 
!             
!!     WRITE(75,*) 'time mean aod(1:192,1,3(0.512))=', aod(:,1,3)
!!     WRITE(75,*) 'time mean asy(1:192,1,3(0.512))=', asy(:,1,3)
!!     WRITE(75,*) 'time mean ssa(1:192,1,3(0.512))=', ssa(:,1,3)
     
!     WRITE(75,*) 'time mean aod(1:192,1,1(0.252))=', aod(:,1,1)
!     WRITE(75,*) 'time mean asy(1:192,1,1(0.252))=', asy(:,1,1)
!     WRITE(75,*) 'time mean ssa(1:192,1,1(0.252))=', ssa(:,1,1)

!     WRITE(75,*) 'time mean aod(1:192,1,8(6.135))=', aod(:,1,8)
!     WRITE(75,*) 'time mean asy(1:192,1,8(6.135))=', asy(:,1,8)
!     WRITE(75,*) 'time mean ssa(1:192,1,8(6.135))=', ssa(:,1,8)
!     

!!     WRITE(75,*) 'z_aer(1,1,1:40)=', z_aer(1,1,:)
!!     WRITE(75,*) 'z_aer(70,1,1:40)=', z_aer(70,1,:)    
!!     WRITE(75,*) 'z_aer(ibMaxPerJB(1),1,1:40)=', z_aer(ibMaxPerJB(1),1,:)
            
!!     CLOSE(75) 
!!           ENDIF 
!!       ENDIF


!----------------------------------------------------------------------
!tar end
!

    !
    ! process sst file
    !

    IF (ifsst /= 0) THEN
       IF (ifsst == -1) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsst, FILE=TRIM(fNameSSTAOI),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
               ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF
          irec=1
          CALL ReadGetSST(nfsst,irec,bfr_in)
          irec=2
          CALL ReadGetSST(nfsst,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,gsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,gsst)
          END IF
          CLOSE(UNIT=nfsst)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (gsst(i,j) > 10.0_r8) THEN
                   gsst(i,j)=-gsst(i,j)
                ELSE IF (gsst(i,j) < 0.0_r8) THEN
                   gsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
                gmax=MAX(gmax,gsst(i,j))
                gmin=MIN(gmin,gsst(i,j))
             END DO
          END DO

          WRITE(UNIT=nfprt,FMT=667) ifsst,gmax,gmin
          IF(SetBCCte == 1)THEN
             ifsst=-1 !PKUBOTA ASSIMILATION
          ELSE
             ifsst=0 
          ENDIF
       ELSE IF (ifsst == 1) THEN
          OPEN(UNIT=nfsst, FILE=TRIM(fNameSSTAOI), FORM='unformatted', ACCESS='sequential',&
               ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF
          READ(UNIT=nfsst)
          month=idate(2)
          DO mm=1,month
             CALL ReadGetSST(nfsst,irec,bfr_in)
          END DO
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,gsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,gsst)
          END IF
          CLOSE(UNIT=nfsst)
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (gsst(i,j) > 10.0_r8) THEN
                   gsst(i,j)=gsst(i,j)
                ELSE IF (gsst(i,j) < 0.0_r8) THEN
                   gsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO
          ifsst=0
       ELSE IF (ifsst == 2.OR. &
            (ifsst == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsst, FILE=TRIM(fNameSSTAOI),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
             ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF

          irec = mf+1
          CALL ReadGetSST(nfsst,irec,bfr_in)

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,xsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,xsst)
          END IF
          
          IF (irec == 13) THEN
             irec=2
          ELSE
             irec=irec+1
          END IF
          CALL ReadGetSST(nfsst,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF

          CLOSE(UNIT=nfsst)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)

                fsst=f2*xsst(i,j)+f1*bfr_out(i,j)
                IF (fsst > gmax) THEN
                   gmax=fsst
                END IF
                IF (fsst < gmin) THEN
                   gmin=fsst
                END IF
                IF (fsst > 10.0_r8) THEN
                   xsst(i,j)=-fsst
                ELSE IF (fsst < 0.0_r8) THEN
                   xsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (tod == 0.0_r8.AND.ifday == 0) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(xsst(i,j)) >= tice) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(gsst(i,j)) >= tice) THEN
                   gsst(i,j)=-tice+1.0e-2_r8
                END IF
             END DO
          END DO

          IF (ifsst == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsst=0
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=666) mf,f1,f2,gmax,gmin
          END IF
       ELSE IF (ifsst == 4) THEN
          IF (intsst > 0) THEN
             fisst=REAL(intsst,r8)
             xday=ifday+tod/86400.0_r8+sstlag
             irec=INT(xday/fisst+1.0e-3_r8+1.0_r8)
             xx1= MOD(xday,fisst)/fisst
             xx2=1.0_r8-xx1
          ELSE
             xx1=f1
             xx2=f2
          END IF
          INQUIRE (IOLENGTH=lrecl) rbrf
          OPEN(UNIT=nfsst,FILE=TRIM(fNameSSTAOI),FORM='unformatted',ACCESS='direct',&
               RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSSTAOI), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadGetSST2(nfsst,bfr_in,irec)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,xsst)
          ELSE
             CALL IJtoIBJB(bfr_in ,xsst)
          END IF
          CALL ReadGetSST2(nfsst,bfr_in,irec+1)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in ,bfr_out)
          ELSE
             CALL IJtoIBJB(bfr_in ,bfr_out)
          END IF
          CLOSE(UNIT=nfsst)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                fsst=xx2*xsst(i,j)+xx1*bfr_out(i,j)
                IF (fsst > gmax) THEN
                   gmax=fsst
                END IF
                IF (fsst < gmin) THEN
                   gmin=fsst
                END IF
                IF (fsst > 10.0_r8) THEN
                   xsst(i,j)=-fsst
                ELSE IF (fsst < 0.0_r8) THEN
                   xsst(i,j)=290.0_r8
                ELSE
                   PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                   STOP "**(ERROR)**"
                END IF
             END DO
          END DO

          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=666) irec,xx1,xx2,gmax,gmin
          END IF

          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                IF (tod == 0.0_r8.AND.ifday == 0) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(xsst(i,j)) >= tice) THEN
                   gsst(i,j)=xsst(i,j)
                ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(gsst(i,j)) >= tice) THEN
                   gsst(i,j)=-tice+1.0e-2_r8
                END IF
             END DO
          END DO

       ELSE IF (ifsst == 5) THEN
         IF (intsst > 0) THEN
            fisst=REAL(intsst,r8)
            xday=ifday+tod/86400.0_r8+sstlag
            irec=INT(xday/fisst+1.0e-3_r8+1.0_r8)
            xx1= MOD(xday,fisst)/fisst
            xx2=1.0_r8-xx1
         ELSE
            xx1=f1
            xx2=f2
         END IF
         INQUIRE (IOLENGTH=lrecl) rbrf
         OPEN(UNIT=nfsst,FILE=TRIM(fNameSSTAOI),FORM='unformatted',ACCESS='direct',&
              RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                 TRIM(fNameSSTAOI), ierr
            STOP "**(ERROR)**"
         END IF
         CALL ReadGetSST2(nfsst,bfr_in,irec2)
         IF (reducedGrid) THEN
            CALL NearestIJtoIBJB(bfr_in ,xsst)
         ELSE
            CALL IJtoIBJB(bfr_in ,xsst)
         END IF
         CALL ReadGetSST2(nfsst,bfr_in,irec2+1)
         IF (reducedGrid) THEN
            CALL NearestIJtoIBJB(bfr_in ,bfr_out)
         ELSE
            CALL IJtoIBJB(bfr_in ,bfr_out)
         END IF
         CLOSE(UNIT=nfsst)
         gmax=-1.0e10_r8
         gmin=+1.0e10_r8
         DO j=1,jbMax
            DO i=1,ibMaxPerJB(j)
               fsst=xx2*xsst(i,j)+xx1*bfr_out(i,j)
               IF (fsst > gmax) THEN
                  gmax=fsst
               END IF
               IF (fsst < gmin) THEN
                  gmin=fsst
               END IF
               IF (fsst > 10.0_r8) THEN
                  xsst(i,j)=-fsst
               ELSE IF (fsst < 0.0_r8) THEN
                  xsst(i,j)=290.0_r8
               ELSE
                  PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SST "
                  STOP "**(ERROR)**"
               END IF
            END DO
         END DO

         IF (nfctrl(23) >= 1) THEN
            WRITE(UNIT=nfprt,FMT=666) irec,xx1,xx2,gmax,gmin
         END IF

         DO j=1,jbMax
            DO i=1,ibMaxPerJB(j)
               IF (tod == 0.0_r8.AND.ifday == 0) THEN
                  gsst(i,j)=xsst(i,j)
               ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(xsst(i,j)) >= tice) THEN
                  gsst(i,j)=xsst(i,j)
               ELSE IF (xsst(i,j) < 0.0_r8.AND.ABS(xsst(i,j)) >= tice) THEN
                  gsst(i,j)=-tice+1.0e-2_r8
               ELSE
                  gsst(i,j)=xsst(i,j)
               END IF
            END DO
         END DO

       ELSE
          WRITE(UNIT=nfprt,FMT=1999)
          STOP
       END IF
    END IF


    IF (schemes ==2) THEN

       !
       ! process ndvi file
       !

       IF (ifndvi /= 0) THEN
          IF (ifndvi == -1) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfndvi, FILE=TRIM(fNameNDVIAOI),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
             ACTION='READ',STATUS='OLD', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF
             irec_ndvi=1
             CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,gndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,gndvi)
             END IF
             CLOSE(UNIT=nfndvi)
             gmax=0.0e0_r8
             gmin=+1.0e0_r8

             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (gndvi(i,j) > 0.0_r8) THEN
                      gndvi(i,j)=gndvi(i,j)
                   ELSE IF (gndvi(i,j) <=0.0_r8) THEN
                      gndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                   gmax=MAX(gmax,gndvi(i,j))
                   gmin=MIN(gmin,gndvi(i,j))
                END DO
             END DO

             WRITE(UNIT=nfprt,FMT=667) ifndvi,gmax,gmin
             ifndvi=0
          ELSE IF (ifndvi == 1) THEN
             OPEN(UNIT=nfndvi, FILE=TRIM(fNameNDVIAOI), FORM='unformatted', ACCESS='sequential',&
                  ACTION='read', STATUS='old', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF
             month=idate(2)
             DO mm=1,month
                CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)
             END DO
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,gndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,gndvi)
             END IF
             CLOSE(UNIT=nfndvi)
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (gndvi(i,j) > 0.0_r8) THEN
                         gndvi(i,j)=gndvi(i,j)
                   ELSE IF (gndvi(i,j) <=0.0_r8) THEN
                      gndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO
             ifndvi=0
          ELSE IF (ifndvi == 2.OR. &
               (ifndvi == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfndvi, FILE=TRIM(fNameNDVIAOI),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn,&
                   ACTION='READ',STATUS='OLD', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF

             irec_ndvi = mf
             CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)

             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,xndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,xndvi)
             END IF
          
             IF (irec_ndvi == 12) THEN
                irec_ndvi=1
             ELSE
                irec_ndvi=irec_ndvi+1
             END IF
             CALL ReadGetSST(nfndvi,irec_ndvi,bfr_in)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in ,bfr_out)
             END IF

             CLOSE(UNIT=nfndvi)
             gmax=0.0e0_r8
             gmin=+1.0e0_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)

                   fndvi=f2_ndvi*xndvi(i,j)+f1_ndvi*bfr_out(i,j)
                   IF (fndvi > gmax) THEN
                      gmax=fndvi
                   END IF
                   IF (fndvi < gmin) THEN
                      gmin=fndvi
                   END IF
                   IF (fndvi > 0.0_r8) THEN
                      xndvi(i,j)=fndvi
                   ELSE IF (fndvi <= 0.0_r8) THEN
                      xndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (tod == 0.0_r8.AND.ifday == 0) THEN
                      gndvi(i,j)=xndvi(i,j)
                   ELSE
                      gndvi(i,j)=xndvi(i,j)
                   END IF
                END DO
             END DO

             IF (ifndvi == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
                ifndvi=0
             END IF
             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=666) mf,f1_ndvi,f2_ndvi,gmax,gmin
             END IF
          ELSE IF (ifndvi == 4) THEN
             IF (intndvi > 0) THEN
                findvi=REAL(intndvi,r8)
                xday=ifday+tod/86400.0_r8+sstlag
                irec_ndvi=INT(xday/findvi+1.0e-3_r8+1.0_r8)
                xx1= MOD(xday,findvi)/findvi
                xx2=1.0_r8-xx1
             ELSE
                xx1=f1_ndvi
                xx2=f2_ndvi
             END IF
             INQUIRE (IOLENGTH=lrecl) bfr_in
             lrecl=lrecl/2
             OPEN(UNIT=nfndvi,FILE=TRIM(fNameNDVIAOI),FORM='unformatted',ACCESS='direct',&
                  RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameNDVIAOI), ierr
                STOP "**(ERROR)**"
             END IF
             CALL ReadGetSST2(nfndvi,bfr_in,irec_ndvi)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,xndvi)
             ELSE
                CALL IJtoIBJB(bfr_in ,xndvi)
             END IF
             CALL ReadGetSST2(nfndvi,bfr_in,irec_ndvi+1)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in ,bfr_out)
             END IF
             CLOSE(UNIT=nfndvi)
             gmax=0.0e0_r8
             gmin=+1.0e0_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   fndvi=xx2*xndvi(i,j)+xx1*bfr_out(i,j)
                   IF (fndvi > gmax) THEN
                      gmax=fndvi
                   END IF
                   IF (fndvi < gmin) THEN
                         gmin=fndvi
                   END IF
                   IF (fndvi > 0.0_r8) THEN
                      xndvi(i,j)=fndvi
                   ELSE IF (fndvi <= 0.0_r8) THEN
                      xndvi(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifndvi=-1 INCORRECT VALUE OF SST "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO

             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=666) irec_ndvi,xx1,xx2,gmax,gmin
             END IF

             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (tod == 0.0_r8.AND.ifday == 0) THEN
                      gndvi(i,j)=xndvi(i,j)
                   ELSE
                      gndvi(i,j)=xndvi(i,j)
                   END IF
                END DO
             END DO

          ELSE IF (ifndvi == 5.AND.intndvi > 0) THEN

             !*(JP)* Eliminei este caso pelas obs do Bonatti e minhas

             PRINT *, " OPTION ifndvi=5 NOT CORRECTLY IMPLEMENTED "
             STOP "**(ERROR)**"

          ELSE
             WRITE(UNIT=nfprt,FMT=1999)
             STOP
          END IF
       END IF
    END IF
    !
    ! process snow file
    !

    IF (ifsnw /= 0) THEN
       IF (ifsnw == 1) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsnw,FILE=TRIM(fNameSnow), FORM='UNFORMATTED', ACCESS='DIRECT', &
                RECL=LRecIn, ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSnow), ierr
             STOP "**(ERROR)**"
          END IF
          irec=1
          CALL ReadGetSNW(nfsnw,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,gsnw)
          ELSE
             CALL IJtoIBJB(bfr_in,gsnw)
          END IF
          CLOSE(UNIT=nfsnw)
          ifsnw=0
       ELSE IF (ifsnw == 2.OR. &
            (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf
          OPEN (UNIT=nfsnw,FILE=TRIM(fNameSnow), FORM='UNFORMATTED', ACCESS='DIRECT', &
                RECL=LRecIn, ACTION='READ',STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSnow), ierr
             STOP "**(ERROR)**"
          END IF
          irec=1
          CALL ReadGetSNW(nfsnw,irec,bfr_in)
          IF (reducedGrid) THEN
             CALL AveBoxIJtoIBJB(bfr_in,gsnw)
          ELSE
             CALL IJtoIBJB(bfr_in,gsnw)
          END IF
          CLOSE(UNIT=nfsnw)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                gmax=MAX(gmax,gsnw(i,j))
                gmin=MIN(gmin,gsnw(i,j))
             END DO
          END DO

          IF (ifsnw == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
             ifsnw=0
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=444) gmax,gmin
          END IF
       ELSE
          WRITE(UNIT=nfprt,FMT=555)
          STOP
       END IF
    END IF

    IF (schemes ==5) THEN
       !
       ! process soil moisture file
       !

       IF (ifslmSib2 /= 0) THEN
          IF (ifslmSib2 == 1) THEN
             rbrfw3d=0.0_r4
             INQUIRE (IOLENGTH=LRecIn) rbrfw3d
             OPEN (UNIT=nfSoilMostSib2,FILE=TRIM(fNameSoilMoistSib2),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilMoistSib2), ierr
                STOP "**(ERROR)**"
             END IF
          
             irec=idate(2)
             CALL ReadGetSLM3D(nfSoilMostSib2,irec,bfrw_in)
     
             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(1:iMax,1:jMax,k),wsib3d(1:ibMax,1:jbMax,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(1:iMax,1:jMax,k),wsib3d(1:ibMax,1:jbMax,k))
                END IF
             END DO
     
             CLOSE(UNIT=nfSoilMostSib2)
             ifslmSib2=0
          ELSE IF (ifslmSib2 == 2.OR. &
               (ifslmSib2 == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
             rbrfw3d=0.0_r4
             INQUIRE (IOLENGTH=LRecIn) rbrfw3d
             OPEN (UNIT=nfSoilMostSib2,FILE=TRIM(fNameSoilMoistSib2),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilMoistSib2), ierr
                STOP "**(ERROR)**"
             END IF
             irec=mf
             CALL ReadGetSLM3D(nfSoilMostSib2,irec,bfrw_in)
 
             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(1:iMax,1:jMax,k),wsib3d(1:ibMax,1:jbMax,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(1:iMax,1:jMax,k),wsib3d(1:ibMax,1:jbMax,k))
                END IF
             END DO

             IF (irec == 12) THEN
                irec=1
             ELSE
                irec=irec+1    
             END IF
  
             CALL ReadGetSLM3D(nfSoilMostSib2,irec,bfrw_in)

             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(1:iMax,1:jMax,k),bfrw_out(1:ibMax,1:jbMax,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(1:iMax,1:jMax,k),bfrw_out(1:ibMax,1:jbMax,k))
                END IF
             END DO
             CLOSE(UNIT=nfSoilMostSib2)
             gmax=-1.0e0_r8
             gmin=+1.0e0_r8
             DO k=1,3
                DO j=1,jbMax
                   DO i=1,ibMaxPerJB(j)
                      wsib3d(i,j,k)=f2*wsib3d(i,j,k)+f1*bfrw_out(i,j,k)
                      gmax=MAX(gmax,wsib3d(i,j,k))
                      gmin=MIN(gmin,wsib3d(i,j,k))
                   END DO
                END DO
             END DO
             IF (ifslmSib2 == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
                ifslmSib2=0
             END IF
             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=222) mf,f1,f2,gmax,gmin
             END IF
          ELSE
             WRITE(UNIT=nfprt,FMT=333)
             STOP
          END IF
       END IF
       
    ELSE
       !
       ! process soil moisture file
       !
       IF (ifslm /= 0) THEN
          IF (ifslm == -1) THEN
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfslm,FILE=TRIM(fNameSoilmsWkl),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                   WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSoilmsWkl), ierr
                STOP "**(ERROR)**"
             END IF
             irec=1
             CALL ReadGetSST(nfslm,irec,bfr_in)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,gslm)
             ELSE
                CALL IJtoIBJB(bfr_in ,gslm)
             END IF
             CLOSE(UNIT=nfslm)
             gmax=-1.0e10_r8
             gmin=+1.0e10_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (gslm(i,j) > 0.0_r8) THEN
                      gslm(i,j)=gslm(i,j)
                   ELSE IF (gslm(i,j) < 1.0_r8) THEN
                      gslm(i,j)=gslm(i,j)
                   ELSE
                      PRINT *, " OPTION ifslm=-1 INCORRECT VALUE OF SOIL MOISTURE "
                      STOP "**(ERROR)**"
                   END IF
                   gmax=MAX(gmax,gslm(i,j))
                   gmin=MIN(gmin,gslm(i,j))
                END DO
             END DO
             WRITE(UNIT=nfprt,FMT=667) ifslm,gmax,gmin
             ifslm=0
          ELSE IF (ifslm == 1) THEN
             !     ifxxx=1    xxx is set to month=idatec(2) in the first call,
             !                but not processed from the subsequent calls.
             !                ifxxx is set to zero after interpolation
             INQUIRE (IOLENGTH=LRecIn) rbrf
             !---------------------------
             OPEN (UNIT=nfslm,FILE=TRIM(fNameSoilmsWkl),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                   WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSoilmsWkl), ierr
                STOP "**(ERROR)**"
             END IF
             ! irec correspond the record by soil layer
             DO irec=1,8
                CALL ReadGetSLM(nfslm,irec,bfr_in)
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfr_in,gslm)
                ELSE
                   CALL IJtoIBJB(bfr_in,gslm)
                END IF
                DO j=1,jbMax
                   DO i=1,ibMaxPerJB(j)
                      wsib3d(i,j,8+1-irec)=gslm(i,j)
                   END DO
                END DO
             END DO
             CLOSE(UNIT=nfslm)
             
             !---------------------------

             OPEN (UNIT=nfslm,FILE=TRIM(fNameSoilms),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                   WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameSoilms), ierr
                STOP "**(ERROR)**"
             END IF
             
             irec=idate(2)
             CALL ReadGetSLM(nfslm,irec,bfr_in)
             IF (reducedGrid) THEN
                CALL AveBoxIJtoIBJB(bfr_in,gslm)
             ELSE
                CALL IJtoIBJB(bfr_in,gslm)
             END IF
             CLOSE(UNIT=nfslm)

             ifslm=0
          ELSE IF (ifslm == 2.OR. &
               (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
               !     ifxxx=2    xxx is interpolated to current day and time every fint
               !                hours synchronized to 00z regardless of initial time.
               !                interpolation is continuous (every time step) if fint<0.
             INQUIRE (IOLENGTH=LRecIn) rbrf
             OPEN (UNIT=nfslm,FILE=TRIM(fNameSoilms),FORM='UNFORMATTED', ACCESS='DIRECT', &
                  ACTION='read', RECL=LRecIn, STATUS='OLD', IOSTAT=ierr) 
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilms), ierr
                STOP "**(ERROR)**"
             END IF

             irec=mf
             CALL ReadGetSLM(nfslm,irec,bfr_in)
 
             IF (reducedGrid) THEN
                CALL AveBoxIJtoIBJB(bfr_in,gslm)
             ELSE
                CALL IJtoIBJB(bfr_in,gslm)
             END IF

             IF (irec == 12) THEN
                irec=1
             ELSE
                irec=irec+1    
             END IF
  
             CALL ReadGetSLM(nfslm,irec,bfr_in)
          
             IF (reducedGrid) THEN
                CALL AveBoxIJtoIBJB(bfr_in,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in,bfr_out)
             END IF

             CLOSE(UNIT=nfslm)
             gmax=-1.0e10_r8
             gmin=+1.0e10_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   gslm(i,j)=f2*gslm(i,j)+f1*bfr_out(i,j)
                   gmax=MAX(gmax,gslm(i,j))
                   gmin=MIN(gmin,gslm(i,j))
                END DO
             END DO
             IF (ifslm == 3.AND.tod == 0.0_r8.AND.ifday == 0) THEN
                ifslm=0
             END IF
             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=222) mf,f1,f2,gmax,gmin
             END IF
          ELSE IF (ifslm == 4) THEN
              !   isoilm=4    soil moisture is linearly interpolated from continuous direct
              !              access data set to current day and time.  data set
              !              is assumed to be spaced every intsoim days or every
              !              calendar month is intsoim < 0. 
             IF (intsoilm > 0) THEN
                fisoilm=REAL(intsoilm,r8)
                xday=ifday+tod/86400.0_r8+soilmlag
                irec=INT(xday/fisoilm+1.0e-3_r8+1.0_r8)
                xx1= MOD(xday,fisoilm)/fisoilm
                xx2=1.0_r8-xx1
             ELSE
                xx1=f1
                xx2=f2
             END IF
             INQUIRE (IOLENGTH=lrecl) rbrf
             OPEN(UNIT=nfslm,FILE=TRIM(fNameSoilms),FORM='unformatted',ACCESS='direct',&
                  RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilms), ierr
                STOP "**(ERROR)**"
             END IF
             CALL ReadGetSST2(nfslm,bfr_in,irec_soilm)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,xsoilm)
             ELSE
                CALL IJtoIBJB(bfr_in ,xsoilm)
             END IF
             CALL ReadGetSST2(nfslm,bfr_in,irec_soilm+1)
             IF (reducedGrid) THEN
                CALL NearestIJtoIBJB(bfr_in ,bfr_out)
             ELSE
                CALL IJtoIBJB(bfr_in ,bfr_out)
             END IF
             CLOSE(UNIT=nfslm)
             gmax=-1.0e10_r8
             gmin=+1.0e10_r8
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   fsoilm=xx2*xsoilm(i,j)+xx1*bfr_out(i,j)
                   IF (fsoilm > gmax) THEN
                      gmax=fsoilm
                   END IF
                   IF (fsoilm < gmin) THEN
                      gmin=fsoilm
                   END IF
                   IF (fsoilm > 0.0_r8) THEN
                      xsoilm(i,j)=fsoilm
                   ELSE IF (fsoilm <= 0.0_r8) THEN
                      xsoilm(i,j)=0.0_r8
                   ELSE
                      PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SOIL MOISTURE "
                      STOP "**(ERROR)**"
                   END IF
                END DO
             END DO

             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=666) irec,xx1,xx2,gmax,gmin
             END IF

             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (tod == 0.0_r8.AND.ifday == 0) THEN
                      gslm(i,j)=xsoilm(i,j)
                   ELSE 
                      gslm(i,j)=xsoilm(i,j)
                   END IF
                END DO
             END DO

          ELSE IF (ifslm == 5) THEN
              !   ifslm=5    soil moisture is linearly interpolated from continuous direct
              !              access data set to current day and time.  data set
              !              is assumed to be spaced every intsoim days or every
              !              calendar month is intsoim < 0. 
             IF (intsoilm > 0) THEN
                fisoilm=REAL(intsoilm,r8)
                xday=ifday+tod/86400.0_r8+soilmlag
                irec=INT(xday/fisoilm+1.0e-3_r8+1.0_r8)
                xx1= MOD(xday,fisoilm)/fisoilm
                xx2=1.0_r8-xx1
             ELSE
                xx1=f1
                xx2=f2
             END IF

             INQUIRE (IOLENGTH=lrecl) rbrfw3d
             OPEN(UNIT=nfslm,FILE=TRIM(fNameSoilms),FORM='unformatted',ACCESS='direct',&
                  RECL=lrecl,ACTION='read', STATUS='old', IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                     TRIM(fNameSoilms), ierr
                STOP "**(ERROR)**"
             END IF

             !    REAL(KIND=r8)                :: bfrw_in  (imax,jmax,3)
             !   REAL(KIND=r8)                :: bfrw_out  (ibmax,jbmax,3)
             !   REAL(KIND=r4)                :: rbrfw3d    (iMax,jMax,3)
             CALL ReadGetSLM3D(nfslm,irec2_soilm,bfrw_in)
             DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(1:iMax,1:jMax,k),wsib3d(1:ibMax,1:jbMax,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(1:iMax,1:jMax,k),wsib3d(1:ibMax,1:jbMax,k))
                END IF
             END DO
             CALL ReadGetSLM3D(nfslm,irec2_soilm+1,bfrw_in)
              DO k=1,3
                IF (reducedGrid) THEN
                   CALL AveBoxIJtoIBJB(bfrw_in(1:iMax,1:jMax,k),bfrw_out(1:ibMax,1:jbMax,k))
                ELSE
                   CALL IJtoIBJB(bfrw_in(1:iMax,1:jMax,k),bfrw_out(1:ibMax,1:jbMax,k))
                END IF
             END DO
             CLOSE(UNIT=nfslm)
             gmax=-1.0e10_r8
             gmin=+1.0e10_r8
             DO k=1,3
                DO j=1,jbMax
                   DO i=1,ibMaxPerJB(j)
                      fsoilm=xx2*wsib3d(i,j,k)+xx1*bfrw_out(i,j,k)
                      wsib3d(i,j,k) =fsoilm
                      IF (fsoilm > gmax) THEN
                         gmax=fsoilm
                      END IF
                      IF (fsoilm < gmin) THEN
                         gmin=fsoilm
                      END IF
                      IF (fsoilm > 0.0_r8) THEN
                         wsib3d(i,j,k)=fsoilm
                      ELSE IF (fsoilm <= 0.0_r8) THEN
                         wsib3d(i,j,k)=0.0_r8
                      ELSE
                         PRINT *, " OPTION ifsst=-1 INCORRECT VALUE OF SOIL MOISTURE "
                         STOP "**(ERROR)**"
                      END IF
                   END DO
                END DO
             END DO

             IF (nfctrl(23) >= 1) THEN
                WRITE(UNIT=nfprt,FMT=666) irec,xx1,xx2,gmax,gmin
             END IF

             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   IF (tod == 0.0_r8.AND.ifday == 0) THEN
                      gslm(i,j)=150*wsib3d(i,j,3)
                   ELSE 
                      gslm(i,j)=150*wsib3d(i,j,3)
                   END IF
                END DO
             END DO
                
          ELSE
             WRITE(UNIT=nfprt,FMT=333)
             STOP
          END IF
       END IF
    END IF

    !
    ! Process CO2 file/field/value
    !

    IF(ifco2.EQ.-1) THEN
       CALL getco2(idatec,co2val)
    ELSEIF(ifco2.EQ.1) THEN
       !CALL READ_MONTH_CO2
    ELSEIF(ifco2.EQ.2) THEN
    ELSEIF(ifco2.EQ.3) THEN
    ELSEIF(ifco2.EQ.4) THEN
    ENDIF

    !
    ! Process ozone file
    !

    IF (ifozone /= 0) THEN
       !   =1    read field from single month file (first call only)
       IF (ifozone == 1) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf3
          OPEN (UNIT=nfozone, FILE=TRIM(fNameOzone), FORM='UNFORMATTED', &
          ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNameOzone), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadOzone(nfozone,bfr_in3,1)
          CLOSE(UNIT=nfozone)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,bfr_out3)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,bfr_out3)
          END IF
          DO k=1,kMax
          !  invert to top to bottom
             km = kMax+1-k
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   gozo(i,k,j)=bfr_out3(i,k,j)
                END DO
             END DO
          END DO
          ifozone=-1
          !   =2    interpolated to current day and time from 12 month clim
          !   =3    interpolated to current day and time from 12 month predicted field
       ELSE IF (ifozone == 2.OR. &
            (ifozone == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=LRecIn) rbrf3
          OPEN (UNIT=nfozone, FILE=TRIM(fNameOzone), FORM='UNFORMATTED', &
          ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i12)") &
                  TRIM(fNameOzone), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadOzone(nfozone,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,gozo)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,gozo)
          END IF

          mf=mf+1
          IF (mf == 13) mf=1
          CALL ReadOzone(nfozone,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,bfr_out3)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,bfr_out3)
          END IF

          CLOSE (UNIT=nfozone)
          gmax=-1.0e18_r8
          gmin=+1.0e18_r8
          DO j=1,jbMax
             DO k=1,kMax
                DO i=1,ibMaxPerJB(j)
                   bfr_out3(i,k,j)=f2*gozo(i,k,j)+f1*bfr_out3(i,k,j)
                END DO
             END DO
             DO k=1,kMax
             !  invert to top to bottom
                km = kMax+1-k
                DO i=1,ibMaxPerJB(j)
                   gozo(i,k,j)=bfr_out3(i,km,j)
                   gmax=MAX(gmax,gozo(i,k,j))
                   gmin=MIN(gmin,gozo(i,k,j))
                END DO
             END DO
          END DO
          IF (ifozone == 3) THEN
             ifozone=-3
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=223) mf,f1,f2,gmax,gmin
          END IF
          !   =4    interpolated from continuous direct access data set to current day and time
       ELSE IF (ifozone == 4) THEN
          WRITE(UNIT=nfprt,FMT=*) 'ERROR: DIRECT ACCESS OZONE FILE NOT IMPLEMENTED! ABORTING...'
          STOP
       END IF
    END IF


    !
    ! Process tracer file
    !

    IF (iftracer /= 0) THEN
       !   =1    read field from single month file (first call only)
       IF (iftracer == 1) THEN
          !INQUIRE (IOLENGTH=LRecIn) rbrf
          !OPEN (UNIT=nftrc, FILE=TRIM(fNametracer), FORM='UNFORMATTED', &
          !ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          OPEN (UNIT=nftrc, FILE=TRIM(fNametracer), FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNametracer), ierr
             STOP "**(ERROR)**"
          END IF
          CALL ReadTracer(nftrc,bfr_in3)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,bfr_out3)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,bfr_out3)
          END IF
          CLOSE(UNIT=nftrc)
          DO k=1,kMax
          !  invert to top to bottom
             km = kMax+1-k
             DO j=1,jbMax
                DO i=1,ibMaxPerJB(j)
                   tracermix(i,k,j)=bfr_out3(i,km,j)
                END DO
             END DO
          END DO
          iftracer=-1
          !   =2    interpolated to current day and time from 12 month clim
          !   =3    interpolated to current day and time from 12 month predicted field
       ELSE IF (iftracer == 2.OR. &
            (iftracer == 3.AND.tod == 0.0_r8.AND.ifday == 0)) THEN
          INQUIRE (IOLENGTH=lrecl) bfr_in3
          lrecl=lrecl/2
          OPEN(UNIT=nftrc,file=TRIM(fNametracer),ACCESS='direct',&
               FORM='unformatted',RECL=lrecl,STATUS='old')
          IF (ierr /= 0) THEN
             WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                  TRIM(fNametracer), ierr
             STOP "**(ERROR)**"
          END IF

          CALL ReadTracer(nftrc,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,tracermix)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,tracermix)
          END IF

          mf=mf+1
          IF (mf == 13) mf=1
          CALL ReadTracer(nftrc,bfr_in3,mf)
          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(bfr_in3 ,bfr_out3)
          ELSE
             CALL IJtoIBJB(bfr_in3 ,bfr_out3)
          END IF

          CLOSE (UNIT=nftrc)
          gmax=-1.0e10_r8
          gmin=+1.0e10_r8
          DO j=1,jbMax
             DO k=1,kMax
                DO i=1,ibMaxPerJB(j)
                   bfr_out3(i,k,j)=f2*tracermix(i,k,j)+f1*bfr_out3(i,k,j)
                   gmax=MAX(gmax,tracermix(i,k,j))
                   gmin=MIN(gmin,tracermix(i,k,j))
                END DO
             END DO
             DO k=1,kMax
             !  invert to top to bottom
                km = kMax+1-k
                DO i=1,ibMaxPerJB(j)
                   tracermix(i,k,j)=bfr_out3(i,km,j)
                   gmax=MAX(gmax,tracermix(i,k,j))
                   gmin=MIN(gmin,tracermix(i,k,j))
                END DO
             END DO
          END DO
          IF (iftracer == 3) THEN
             iftracer=-3
          END IF
          IF (nfctrl(23) >= 1) THEN
             WRITE(UNIT=nfprt,FMT=223) mf,f1,f2,gmax,gmin
          END IF
          !   =4    interpolated from continuous direct access data set to current day and time
       ELSE IF (iftracer == 4) THEN
          WRITE(UNIT=nfprt,FMT=*) 'ERROR: DIRECT ACCESS OZONE FILE NOT IMPLEMENTED! ABORTING...'
          STOP
       END IF
    END IF


222 FORMAT(' SOILM   START MONTH=',i2,'  F1,F2=',2f6.3,'  MAX,MIN=',2e12.5)
223 FORMAT(' OZONE   START MONTH=',i2,'  F1,F2=',2f6.3,'  MAX,MIN=',2e12.5)
333 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SOILM  INTERPOLATION')
444 FORMAT(' SNOW HAS ONLY ONE FILE','  MAX,MIN=',2E12.5)
555 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT SNOW   INTERPOLATION')
666 FORMAT(' SST START REC (MONTH+2) =',I5, &
         '  F1,F2=',2G13.6,'  MAX,MIN=',2G12.5)
667 FORMAT(' SST:  IFSST=',I2,'  MAX,MIN=',2G12.5)
888 FORMAT(' ALBEDO  START MONTH=',I2, &
         '  F1,F2=',2F6.3,'  MAX,MIN=',2E12.5)
999 FORMAT(' ABNORMAL END IN SUBR.GETSBC AT ALBEDO INTERPOLATION')
1999 FORMAT('ABNORMAL END IN SUBR.GETSBC AT SST   INTERPOLATION')
  END SUBROUTINE getsbc


  SUBROUTINE GetRecWgtMonthlySST &
       (idate, idatec, tod, labelsi, labelsj, &
       irec, f1, f2, mra, mrb,monl)

    IMPLICIT NONE

    ! Computes the Corresponding Records to do Linear
    ! Time Interpolation and the Respectives Weights.

    INTEGER            , INTENT (IN) :: idate(4), idatec(4),monl(12)
    REAL      (KIND=r8), INTENT (IN) :: tod
    CHARACTER (LEN=10) , INTENT (IN) :: labelsi, labelsj

    INTEGER, INTENT (OUT) :: irec, mra, mrb
    REAL (KIND=r8), INTENT (OUT) :: f1, f2

    ! Local Constants
    INTEGER :: ysi, msi, dsi,si, ysj, msj, dsj, ndij, nd, &
         tmca, tmcb, tmcf,yi,mi,di,hi,LenYearbyDay,nday2y
    REAL (KIND=r8) :: xday2
    REAL (KIND=r8) :: xday, zdayf, zdaya, zdayb, tc

    ! Get Year, Month and Day of the Initial and Second Medium Date
    ! for SST Direct Access File Data

    READ (labelsi(1:4), '(I4)') ysi
    READ (labelsi(5:6), '(I2)') msi
    READ (labelsi(7:8), '(I2)') dsi
    READ (labelsj(1:4), '(I4)') ysj
    READ (labelsj(5:6), '(I2)') msj
    READ (labelsj(7:8), '(I2)') dsj

    ! Lag of Days for SST Data:
    ! Just for Checking if the Scale is a Month
    ndij=0
    IF (msi+1 <= msj-1) THEN
       DO nd=msi+1,msj-1
          ndij=ndij+monl(nd)
       END DO
    ELSE
       DO nd=msi+1,12
          ndij=ndij+monl(nd)
       END DO
       DO nd=1,msj-1
          ndij=ndij+monl(nd)
       END DO
    END IF
    ndij=ndij+monl(msi)-dsi+dsj+365*(ysj-ysi-1)

    ! Check for Monthly Scale SST Data
    IF (ABS(ndij) <= 27 .OR. ABS(ndij) >= 32) THEN
       WRITE (UNIT=0, FMT='(/,A)') ' *** Error: The SST Data Is Not On Monthly Scale   ***'
       WRITE (UNIT=0, FMT='(/,A,I8,12X,A,/)') ' *** Lag Of Days For SST Data: ', ndij, '***'
       WRITE (UNIT=0, FMT='(A,/)') ' *** Program STOP: SUBROUTINE GetRecWgtMonthlySST  ***'
       STOP
    END IF

    ! Length in Days of the Date of Forecasting
    yi=idatec(4)
    mi=idatec(2)
    di=idatec(3)
    hi=idatec(1)
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
    tmcf=monl(idatec(2))
    IF (idatec(2) == 2 .AND. LenYearbyDay == 366) tmcf=29
    ! Medium Day of the Month of Forecasting
    zdayf=0.5_r8*REAL(tmcf,r8)+1.0_r8
    ! Fractional Day of Forecasting
    tc=REAL(idate(1),r8)/24.0_r8+tod/86400.0_r8
    ! Correcting Factor if Necessary (tc is in Days)
    IF (tc >= 1.0_r8) tc=tc-1.0_r8
    xday=REAL(idatec(3),r8)+tc
    ! Getting the Corresponding Record in SST Data
    irec=12-msi+idatec(2)+12*(idatec(4)-ysi-1)+2
    IF (xday >= zdayf) irec=irec+1

    ! Months for the Linear Time Interpolation Related to the Records
    mra=MOD(irec-3+msi,12)
    IF (mra == 0) mra=12
    mrb=mra+1
    IF (mrb > 12) mrb=1

    ! Length in Days for the First Month of Interpolation
    tmca=monl(mra)
    yi=ysi
    mi=mra
    di=dsi
    hi=0
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
    IF (mra == 2 .AND. LenYearbyDay == 366) tmca=29
    ! Medium Fracitonal Day for the First Month of Interpolation
    zdaya=0.5_r8*REAL(tmca,r8)+1.0_r8-REAL(tmca,r8)
    ! Length in Days for the Second Month of Interpolation
    tmcb=monl(mrb)
    yi=ysj
    mi=mrb
    di=dsj
    hi=0
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
    IF (mrb == 2 .AND. LenYearbyDay == 366) tmcb=29
    ! Medium Fracitonal Day for the Second Month of Interpolation
    zdayb=0.5_r8*REAL(tmcb,r8)+1.0_r8
    ! Scaling Fractional Day of Forecasting, if Necessary
    IF (xday >= zdayf) xday=xday-REAL(tmca,r8)
    ! Interpolation Factors
    f1=(xday-zdaya)/(zdayb-zdaya)
    f2=1.0_r8-f1

  END SUBROUTINE GetRecWgtMonthlySST




  SUBROUTINE GetRecWgtDailySST &
       (idate, idatec, tod, labelsi, labelsj, &
       irec,irec2, f1, f2, mra, mrb,monl)

    IMPLICIT NONE

    ! Computes the Corresponding Records to do Linear
    ! Time Interpolation and the Respectives Weights.

    INTEGER, INTENT (IN) :: idate(4), idatec(4),monl(12)
    REAL (KIND=r8), INTENT (IN) :: tod
    CHARACTER (LEN=10), INTENT (IN) :: labelsi, labelsj

    INTEGER, INTENT (OUT) :: irec,irec2, mra, mrb
    REAL (KIND=r8), INTENT (OUT) :: f1, f2

    ! Local Constants
    INTEGER :: ysi, msi, dsi, ysj, msj, dsj, ndij, nd, &
         tmca, tmcb, tmcf
    INTEGER :: yi,mi,di,hi,LenYearbyDay,nday2y
    REAL (KIND=r8) :: xday2
    REAL (KIND=r8) :: xday, zdayf, zdaya, zdayb, tc
    INTEGER :: iday
    INTEGER :: firstday
    INTEGER :: lastday

    INTEGER :: imonth
    INTEGER :: firstmonth
    INTEGER :: lastmonth

    INTEGER :: iyear
    INTEGER :: firstyear
    INTEGER :: lastyear
    INTEGER :: ndays
    ! Get Year, Month and Day of the Initial and Second Medium Date
    ! for SST Direct Access File Data

    READ (labelsi(1:4), '(I4)') ysi
    READ (labelsi(5:6), '(I2)') msi
    READ (labelsi(7:8), '(I2)') dsi
    READ (labelsj(1:4), '(I4)') ysj
    READ (labelsj(5:6), '(I2)') msj
    READ (labelsj(7:8), '(I2)') dsj

    ! Lag of Days for SST Data:
    ! Just for Checking if the Scale is a Month
    ndij=0
    IF (msi+1 <= msj-1) THEN
       DO nd=msi+1,msj-1
          ndij=ndij+monl(nd)
       END DO
    ELSE
       DO nd=msi+1,12
          ndij=ndij+monl(nd)
       END DO
       DO nd=1,msj-1
          ndij=ndij+monl(nd)
       END DO
    END IF
    ndij=ndij+monl(msi)-dsi+dsj+365*(ysj-ysi-1)

    ! Check for Monthly Scale SST Data
    IF (ABS(ndij) <= 27 .OR. ABS(ndij) >= 32) THEN
       WRITE (UNIT=0, FMT='(/,A)') ' *** Error: The SST Data Is Not On Monthly Scale   ***'
       WRITE (UNIT=0, FMT='(/,A,I8,12X,A,/)') ' *** Lag Of Days For SST Data: ', ndij, '***'
       WRITE (UNIT=0, FMT='(A,/)') ' *** Program STOP: SUBROUTINE GetRecWgtDailySST  ***'
       STOP
    END IF

    ! Length in Days of the Date of Forecasting
    tmcf=monl(idatec(2))
    yi=idatec(4)
    mi=idatec(2)
    di=idatec(3)
    hi=idatec(1)
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
    IF (idatec(2) == 2 .AND. LenYearbyDay == 366) tmcf=29
    ! Medium Day of the Month of Forecasting
    zdayf=0.5_r8*REAL(tmcf,r8)+1.0_r8
    ! Fractional Day of Forecasting
    tc=REAL(idate(1),r8)/24.0_r8+tod/86400.0_r8
    ! Correcting Factor if Necessary (tc is in Days)
    IF (tc >= 1.0_r8) tc=tc-1.0_r8
    xday=REAL(idatec(3),r8)+tc
    ! Getting the Corresponding Record in SST Data
    irec=12-msi+idatec(2)+12*(idatec(4)-ysi-1)+2
    FirstDay=dsi 
    FirstMonth=msi
    FirstYear=ysi
    LastYear=idatec(4)
    irec2=2
    DO iyear=FirstYear,LastYear
       IF(iyear < LastYear)THEN
          LastMonth=12
          DO imonth=FirstMonth,LastMonth
             ndays=monl(imonth)
             yi=iyear
             mi=imonth
             di=1
             hi=0
             CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
             IF(LenYearbyDay == 366 .and. imonth == 2)ndays=29
             LastDay=ndays
             DO iday=FirstDay,LastDay
                irec2=irec2+1
             END DO
             FirstDay=1
          END DO
       ELSE IF(iyear == LastYear)THEN
          LastMonth=idatec(2)
          DO imonth=FirstMonth,LastMonth
             IF(imonth<LastMonth)THEN
                ndays=monl(imonth)
               yi=iyear
               mi=imonth
               di=1
               hi=0
               CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
              IF(LenYearbyDay == 366 .and. imonth == 2)ndays=29
               LastDay=ndays
             ELSE
                LastDay=idatec(3)
             END IF
             DO iday=FirstDay,LastDay
                irec2=irec2+1
             END DO
             FirstDay=1
          END DO
       END IF
       FirstMonth=1
    END DO

    IF (xday >= zdayf) irec=irec+1

    ! Months for the Linear Time Interpolation Related to the Records
    mra=MOD(irec-3+msi,12)
    IF (mra == 0) mra=12
    mrb=mra+1
    IF (mrb > 12) mrb=1

    ! Length in Days for the First Month of Interpolation
    tmca=monl(mra)
    yi=ysi
    mi=mra
    di=dsi
    hi=0
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay) 
    IF (mra == 2 .AND. LenYearbyDay == 366 ) tmca=29
    ! Medium Fracitonal Day for the First Month of Interpolation
    zdaya=0.5_r8*REAL(tmca,r8)+1.0_r8-REAL(tmca,r8)
    ! Length in Days for the Second Month of Interpolation
    tmcb=monl(mrb)
    yi=ysj
    mi=mrb
    di=dsj
    hi=0
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay) 
    IF (mrb == 2 .AND.  LenYearbyDay == 366) tmcb=29
    ! Medium Fracitonal Day for the Second Month of Interpolation
    zdayb=0.5_r8*REAL(tmcb,r8)+1.0_r8
    ! Scaling Fractional Day of Forecasting, if Necessary
    IF (xday >= zdayf) xday=xday-REAL(tmca,r8)
    ! Interpolation Factors
    f1=(xday-zdaya)/(zdayb-zdaya)
    f2=1.0_r8-f1
    ! Interpolation Factors
    f1=(86400.0_r8-tod)/86400.0_r8 
    f2=(tod)/86400.0_r8 

  END SUBROUTINE GetRecWgtDailySST


  SUBROUTINE GetWeightsOld (yrl,monl,idatec, tod, f1, f2,mf)

    IMPLICIT NONE

    ! Computes Weights as in getsbc:

    INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)
    INTEGER, INTENT (IN) :: idatec(4)
    INTEGER, INTENT (IN) :: monl(12)
    REAL (KIND=r8), INTENT (IN) :: tod
    REAL (KIND=r8), INTENT (IN) :: yrl
    REAL (KIND=r8), INTENT (OUT) :: f1, f2
    INTEGER,  INTENT (OUT):: mf
    INTEGER :: mon, mnl, mn, mnlf, mnln
    INTEGER :: yi,mi,di,hi,LenYearbyDay,nday2y
    REAL (KIND=r8) :: xday2
    REAL (KIND=r8) :: yday, add
    LOGICAL :: ly

    mon=idatec(2)
    yday=REAL(idatec(3),r8)+REAL(idatec(1),r8)/24.0_r8+MOD(tod,3600.0_r8)/86400.0_r8
    mf=mon-1
    yi=idatec(4)
    mi=idatec(2)
    di=idatec(3)
    hi=idatec(1)
    CALL jull(yi,mi,di,hi,tod,xday2,nday2y,LenYearbyDay)
    ly= yrl == 365.25_r8 .AND. LenYearbyDay == 366
    mnl=monl(mon)
    IF (ly .AND. mon == 2) mnl=29
    ! Em getsbc seria apenas >
    ! As consideracoes de interpolacao leva a >=
    IF (yday >= 1.0_r8+0.5_r8*REAL(mnl,r8)) mf=mon
    mn=mf+1
    IF (mf < 1) mf=12
    IF (mn > 12) mn=1
    mnlf=monl(mf)
    IF (ly .AND. mf == 2) mnlf=29
    add=0.5_r8*REAL(mnlf,r8)-1.0_r8
    IF (mf == mon) add=-add-2.0_r8
    mnln=monl(mn)
    IF (ly .AND. mn == 2) mnln=29
    f1=2.0_r8*(yday+add)/REAL(mnlf+mnln,r8)
    f2=1.0_r8-f1

  END SUBROUTINE GetWeightsOld

  ! gread : reads in history carrying variables for one time step,
  !         surface geopotential, and sigma coordinate levels.
  !         checks sigma coordinate levels for consistency.



  SUBROUTINE gread(n, ifday, tod, idate, idatec, &
       qgzs, qlnp, qtmp, qdiv, qrot, qq, a, b)
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(OUT) :: ifday
    REAL(KIND=r8),    INTENT(OUT) :: tod
    INTEGER, INTENT(OUT) :: idate(4)
    INTEGER, INTENT(OUT) :: idatec(4)
    REAL(KIND=r8),    INTENT(OUT) :: qgzs(2*mymnMax)
    REAL(KIND=r8),    INTENT(OUT) :: qlnp(2*mymnMax)
    REAL(KIND=r8),    INTENT(OUT) :: qtmp(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qdiv(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qrot(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qq  (2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN ) :: a(kMax+1)
    REAL(KIND=r8),    INTENT(IN ) :: b(kMax+1)

    INTEGER              :: k, i1, i2, m, nn, mm
    REAL(KIND=r8)                 :: aux(2*mnMax,kMax)
    REAL(KIND=r8)                 :: dphi(kMax+1)
    REAL(KIND=r8)                 :: dlam(kMax+1)
    !
    !     spectral data file format
    !     hour,idate(4),si( kMax+1 ),sl( kMax )
    !     zln qlnp qtmp qdiv qrot
    !
    !     Spectral fields are reversed here (from top to bottom of atmosphere)
    !     since in the files used in pre and pos processing they range from
    !     bottom to top
    !
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=999) n,kMax,kMax+1
    CALL GReadHead(n, ifday, tod, idate, idatec, dphi, dlam, kMax)

    IF (maxnodes.eq.1) THEN
       !
       ! Spectral Coefficients of Orography (m)
       !
       CALL GReadField(n, qgzs)
       !
       ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
       !
       CALL GReadField(n, qlnp)
       !
       ! Spectral Coefficients of Virtual Temp (K)
       !
       DO k = kMax,1,-1
          CALL GReadField(n, qtmp(:,k))
       END DO
       !
       ! Spectral Coefficients of Divergence and Vorticity (1/seg)
       !
       IF(TRIM(StrFormat) == 'old')THEN
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Divergence
             CALL GReadField(n, qdiv(:,k))
             !Vorticity
             CALL GReadField(n, qrot(:,k))
          END DO
       else if(TRIM(StrFormat) == 'new') then
          !
          !Spectral Coefficients of Divergence
          !
          DO k = kMax,1,-1
             !Divergence
             CALL GReadField(n, qdiv(:,k))
          END DO
          DO k = kMax,1,-1
             !Vorticity
             CALL GReadField(n, qrot(:,k))
          END DO
       ELSE
         CALL FatalError('Invalid StrFormat. see at the namelist, Aborting Model!')
       END IF
       !
       ! Spectral Coefficients of Specific Humidity (kg/kg)
       !
       DO k = kMax,1,-1
          CALL GReadField(n, qq(:,k))
       END DO
       !
    ELSE
       !
       ! Spectral Coefficients of Orography (m)
       !
       CALL GReadField(n, aux(:,1))
       DO mm=1,mymmax
          m = msinproc(mm,myid_four)
          i1 = 2*mnmap(m,m)-1
          i2 = 2*mymnmap(mm,m)-1
          DO nn=0,2*(mmax-m)+1
             qgzs(i2+nn) = aux(i1+nn,1)
          ENDDO
       ENDDO
       !
       ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
       !
       CALL GReadField(n, aux(:,1))
       DO mm=1,mymmax
          m = msinproc(mm,myid_four)
          i1 = 2*mnmap(m,m)-1
          i2 = 2*mymnmap(mm,m)-1
          DO nn=0,2*(mmax-m)+1
             qlnp(i2+nn) = aux(i1+nn,1)
          ENDDO
       ENDDO
       !
       ! Spectral Coefficients of Virtual Temp (K)
       !
       CALL GReadField(n, aux)
       DO k=myfirstlev,mylastlev
          DO mm=1,mymmax
             m = msinproc(mm,myid_four)
             i1 = 2*mnmap(m,m)-1
             i2 = 2*mymnmap(mm,m)-1
             DO nn=0,2*(mmax-m)+1
                qtmp(i2+nn,k+1-myfirstlev) = aux(i1+nn,kmax+1-k)
             ENDDO
          ENDDO
       ENDDO
       !
       ! Spectral Coefficients of Divergence and Vorticity (1/seg)
       !
       IF(TRIM(StrFormat) == 'old')THEN
          DO k = kMax,1,-1
             !Divergence
             CALL GReadField(n, aux(:,1))
             !Vorticity
             CALL GReadField(n, aux(:,2))
             IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                DO mm=1,mymmax
                   m = msinproc(mm,myid_four)
                   i1 = 2*mnmap(m,m)-1
                   i2 = 2*mymnmap(mm,m)-1
                   DO nn=0,2*(mmax-m)+1
                     qdiv(i2+nn,k+1-myfirstlev) = aux(i1+nn,1)
                     qrot(i2+nn,k+1-myfirstlev) = aux(i1+nn,2)
                   ENDDO
                ENDDO
             END IF
          ENDDO
       else if(TRIM(StrFormat) == 'new') then
          !
          !Spectral Coefficients of Divergence
          !
          DO k = kMax,1,-1
             !Divergence
             CALL GReadField(n, aux(:,1))
             IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                DO mm=1,mymmax
                   m = msinproc(mm,myid_four)
                   i1 = 2*mnmap(m,m)-1
                   i2 = 2*mymnmap(mm,m)-1
                   DO nn=0,2*(mmax-m)+1
                     qdiv(i2+nn,k+1-myfirstlev) = aux(i1+nn,1)
                   ENDDO
                ENDDO
             END IF
          ENDDO
          !
          !Spectral Coefficients of Vorticity
          !
          DO k = kMax,1,-1
             !Vorticity
             CALL GReadField(n, aux(:,2))
             IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                DO mm=1,mymmax
                   m = msinproc(mm,myid_four)
                   i1 = 2*mnmap(m,m)-1
                   i2 = 2*mymnmap(mm,m)-1
                   DO nn=0,2*(mmax-m)+1
                     qrot(i2+nn,k+1-myfirstlev) = aux(i1+nn,2)
                   ENDDO
                ENDDO
             END IF
          ENDDO

       ELSE
         CALL FatalError('Invalid StrFormat. see at the namelist, Aborting Model!')
       END IF
       !
       ! Spectral Coefficients of Specific Humidity (kg/kg)
       !
       CALL GReadField(n, aux)
       DO k=myfirstlev,mylastlev
          DO mm=1,mymmax
             m = msinproc(mm,myid_four)
             i1 = 2*mnmap(m,m)-1
             i2 = 2*mymnmap(mm,m)-1
             DO nn=0,2*(mmax-m)+1
                qq(i2+nn,k+1-myfirstlev) = aux(i1+nn,kmax+1-k)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

    CLOSE(UNIT=n)

    !cdir novector
    DO k=1, kMax+1
       dphi(k)=dphi(k)-a(k)
       dlam(k)=dlam(k)-b(k)
    END DO
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=100) (dphi(k),k=1, kMax+1 )
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=100) (dlam(k),k=1, kMax+1 )
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=101) n,ifday,tod,idate,idatec
100 FORMAT(' ', 13(E9.3))
101 FORMAT (' ', 'IF ABOVE TWO ROWS NOT ZERO, ', &
         'INCONSISTENCY IN HYBRID COORDINATE  DEFINITION ON N=',I2/' AT DAY=',I8, &
         ' TIME=',F8.1,' STARTING',3I3,I5,' CURRENT',3I3,I5)
999 FORMAT(' N,KMAX,KMAXP=',3I4)
  END SUBROUTINE gread


  SUBROUTINE gread4 (n, ifday, tod, idate, idatec, &
       qgzs, qlnp, qtmp, qdiv, qrot, qq, a, b, dodyn, nfdyn)
    INTEGER, INTENT(IN ) :: n
    INTEGER, INTENT(OUT) :: ifday
    REAL(KIND=r8),    INTENT(OUT) :: tod
    INTEGER, INTENT(OUT) :: idate(4)
    INTEGER, INTENT(OUT) :: idatec(4)
    REAL(KIND=r8),    INTENT(OUT) :: qgzs(2*mymnMax)
    REAL(KIND=r8),    INTENT(OUT) :: qlnp(2*mymnMax)
    REAL(KIND=r8),    INTENT(OUT) :: qtmp(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qdiv(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qrot(2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(OUT) :: qq  (2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN)  :: a(kMax+1)
    REAL(KIND=r8),    INTENT(IN)  :: b(kMax+1)
    LOGICAL, INTENT(IN ) :: dodyn
    INTEGER, INTENT(IN ) :: nfdyn

    INTEGER              :: k, i1, i2, m, nn, mm
    REAL(KIND=r8)                 :: aux(2*mnMax,kMax)
    REAL(KIND=r8)                 :: aux1(2*mnMax)
    REAL(KIND=r8)                 :: dphi(kMax+1)
    REAL(KIND=r8)                 :: dlam(kMax+1)


    INTEGER(KIND=i4)        :: ifday4
    INTEGER(KIND=i4)        :: idat4(4)
    INTEGER(KIND=i4)        :: idat4c(4)
    REAL(KIND=r4)           :: tod4
    REAL(KIND=r4)           :: dph4(kmax+1)
    REAL(KIND=r4)           :: dla4(kmax+1)
    INTEGER, SAVE        :: ifdyn = 0
    !
    !     spectral data file format
    !     hour,idate(4),a( kmax+1 ),b( kmax+1 )
    !     zln qlnp qtmp qdiv qrot
    !
    !     Spectral fields are reversed here (from top to bottom of atmosphere)
    !     since in the files used in pre and pos processing they range from
    !     bottom to top
    !
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=999) n,kmax,kmax+1
    CALL ReadHead(n, ifday4, tod4, idat4, idat4c, dph4, dla4, kMax)
    ifday=ifday4
    tod=tod4
    !dph4=si
    !dla4=sl
    !si=dph4
    !sl=dla4
    DO k=1,4
       idate(k)=idat4(k)
       idatec(k)=idat4c(k)
    ENDDO
    DO k=1,kmax+1
       dphi(k)=dph4(k)
       dlam(k)=dla4(k)
    ENDDO


    IF (maxnodes.eq.1) THEN
       !
       ! Spectral Coefficients of Orography (m)
       !
       CALL ReadField(n, qgzs)

       IF (ifdyn .EQ. 0) THEN
          ifdyn=1
          IF (dodyn) THEN
             WRITE(UNIT=nfprt,FMT='(A,I5,A,F15.2,A)') ' ifday=',ifday4,' tod=',tod4,' dyn'
             WRITE(UNIT=nfdyn) ifday4,tod4
             WRITE(UNIT=nfdyn) qgzs
          END IF
       END IF
       !
       ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
       !
       CALL ReadField(n, qlnp)
       !      transform surface pressure from cbar to Pascal
       !
       IF (HaveM1) qlnp(1) = qlnp(1) + log(1000._r8) * root2
       !
       ! Spectral Coefficients of Virtual Temp (K)
       !
       DO k = kMax,1,-1
          CALL ReadField(n, qtmp(:,k))
       END DO
       !
       ! Spectral Coefficients of Divergence and Vorticity (1/seg)
       !
       IF(TRIM(StrFormat) == 'old')THEN
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Divergence
             CALL ReadField(n, qdiv(:,k))
             !Vorticity
             CALL ReadField(n, qrot(:,k))
          END DO
       else if(TRIM(StrFormat) == 'new') then
          DO k = kMax,1,-1
             !Divergence
             CALL ReadField(n, qdiv(:,k))
          END DO
          DO k = kMax,1,-1
             !Vorticity
             CALL ReadField(n, qrot(:,k))
          END DO
       ELSE
         CALL FatalError('Invalid StrFormat. see at the namelist, Aborting Model!')
       END IF
       !
       ! Spectral Coefficients of Specific Humidity (kg/kg)
       !
       DO k = kMax,1,-1
          CALL ReadField(n, qq(:,k))
       END DO
    ELSE
       !
       ! Spectral Coefficients of Orography (m)
       !
       CALL ReadField(n, aux1)
       IF(myid.eq.0) THEN
          IF (ifdyn .EQ. 0) THEN
             ifdyn=1
             IF (dodyn) THEN
                WRITE (UNIT=nfprt,FMT='(A,I5,A,F15.2,A)') ' ifday=',ifday4,' tod=',tod4,' dyn'
                WRITE (UNIT=nfdyn) ifday4,tod4
                WRITE (UNIT=nfdyn) aux1
             END IF
          END IF
       END IF

       DO mm=1,mymmax
          m = msinproc(mm,myid_four)
          i1 = 2*mnmap(m,m)-1
          i2 = 2*mymnmap(mm,m)-1
          DO nn=0,2*(mmax-m)+1
             qgzs(i2+nn) = aux1(i1+nn)
          ENDDO
       ENDDO
       !
       ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
       !
       CALL ReadField(n, aux1)
       DO mm=1,mymmax
          m = msinproc(mm,myid_four)
          i1 = 2*mnmap(m,m)-1
          i2 = 2*mymnmap(mm,m)-1
          DO nn=0,2*(mmax-m)+1
             qlnp(i2+nn) = aux1(i1+nn)
          ENDDO
       ENDDO
!      transform surface pressure from cbar to Pascal
!
       IF (HaveM1) qlnp(1) = qlnp(1) + log(1000._r8) * root2
       !
       ! Spectral Coefficients of Virtual Temp (K)
       !
       CALL ReadField(n, aux)
       DO k=myfirstlev,mylastlev
          DO mm=1,mymmax
             m = msinproc(mm,myid_four)
             i1 = 2*mnmap(m,m)-1
             i2 = 2*mymnmap(mm,m)-1
             DO nn=0,2*(mmax-m)+1
                qtmp(i2+nn,k+1-myfirstlev) = aux(i1+nn,kmax+1-k)
             ENDDO
          ENDDO
       ENDDO
       !
       ! Spectral Coefficients of Divergence and Vorticity (1/seg)
       !
       IF(TRIM(StrFormat) == 'old')THEN
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Divergence
             CALL ReadField(n, aux(:,1))
             !Vorticity
             CALL ReadField(n, aux(:,2))
             IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                DO mm=1,mymmax
                   m = msinproc(mm,myid_four)
                   i1 = 2*mnmap(m,m)-1
                   i2 = 2*mymnmap(mm,m)-1
                   DO nn=0,2*(mmax-m)+1
                      qdiv(i2+nn,k+1-myfirstlev) = aux(i1+nn,1)
                      qrot(i2+nn,k+1-myfirstlev) = aux(i1+nn,2)
                   ENDDO
                ENDDO
             END IF
          ENDDO
       else if(TRIM(StrFormat) == 'new') then
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Divergence
             CALL ReadField(n, aux(:,1))
             IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                DO mm=1,mymmax
                   m = msinproc(mm,myid_four)
                   i1 = 2*mnmap(m,m)-1
                   i2 = 2*mymnmap(mm,m)-1
                   DO nn=0,2*(mmax-m)+1
                      qdiv(i2+nn,k+1-myfirstlev) = aux(i1+nn,1)
                   ENDDO
                ENDDO
             END IF
          ENDDO
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Vorticity
             CALL ReadField(n, aux(:,2))
             IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                DO mm=1,mymmax
                   m = msinproc(mm,myid_four)
                   i1 = 2*mnmap(m,m)-1
                   i2 = 2*mymnmap(mm,m)-1
                   DO nn=0,2*(mmax-m)+1
                      qrot(i2+nn,k+1-myfirstlev) = aux(i1+nn,2)
                   ENDDO
                ENDDO
             END IF
          ENDDO
       ELSE
         CALL FatalError('Invalid StrFormat. see at the namelist, Aborting Model!')
       END IF
       !
       ! Spectral Coefficients of Specific Humidity (kg/kg)
       !
       CALL ReadField(n, aux)
       DO k=myfirstlev,mylastlev
          DO mm=1,mymmax
             m = msinproc(mm,myid_four)
             i1 = 2*mnmap(m,m)-1
             i2 = 2*mymnmap(mm,m)-1
             DO nn=0,2*(mmax-m)+1
                qq(i2+nn,k+1-myfirstlev) = aux(i1+nn,kmax+1-k)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

    CLOSE(UNIT=n)

    DO k=1, kmax+1
       dphi(k)=dphi(k)-a(k)
       dlam(k)=dlam(k)-b(k)
    END DO

    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=100) (dphi(k),k=1, kmax+1 )
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=100) (dlam(k),k=1, kmax+1 )
    IF(nfctrl(35).GE.1)WRITE(UNIT=nfprt,FMT=101) n,ifday,tod,idate,idatec
100 FORMAT(' ', 13(E9.3))
101 FORMAT (' IF ABOVE TWO ROWS NOT ZERO, ', &
         'INCONSISTENCY IN HYBRID COORDINATE  DEFINITION ON N=',I2/' AT DAY=',I8, &
         ' TIME=',F8.1,' STARTING',3I3,I5,' CURRENT',3I3,I5)
999 FORMAT(' N,KMAX,KMAXP=',3I4)
  END SUBROUTINE gread4

  !     gwrite : writes out the surface geopotential and history carrying
  !              fields of the spectral model after first inverting the
  !              laplacian to recapture the surface geopotential field.
  !       Fields are written from bottom to top of the atmosphere, in the 
  !       way they are used in pre and pos processing.

  SUBROUTINE gwrite(n, ifday, tod, idate, idatec, &
       qlnp, qtmp, qdiv, qrot, qq, a, b, qgzs)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: ifday
    REAL(KIND=r8),    INTENT(IN) :: tod
    INTEGER, INTENT(IN) :: idate(4)
    INTEGER, INTENT(IN) :: idatec(4)
    REAL(KIND=r8),    INTENT(IN) :: qgzs (2*mymnMax)
    REAL(KIND=r8),    INTENT(INOUT) :: qlnp (2*mymnMax)
    REAL(KIND=r8),    INTENT(IN) :: qtmp (2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN) :: qdiv (2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN) :: qrot (2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN) :: qq   (2*mymnMax,kMaxloc)
    REAL(KIND=r8),    INTENT(IN) :: a(kMax+1)
    REAL(KIND=r8),    INTENT(IN) :: b(kMax+1)

    INTEGER             :: k
    REAL(KIND=r8)                :: aux (2*mnMax,kMax)
    REAL(KIND=r8)                :: aux1(2*mnMax,kMax)

    IF(myid.eq.0)  CALL GWriteHead(n, ifday, tod, idate, idatec, a, b)

    IF(maxnodes.eq.1) THEN
       !
       ! Spectral Coefficients of Orography (m)
       !
       CALL GWriteField(n, qgzs)
       !
       ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
       !
       !      transform surface pressure from Pascal to cbar
       !
       qlnp(1) = qlnp(1) - log(1000._r8) * root2

       CALL GWriteField(n, qlnp)
       !      transform surface pressure back from cbar to Pascal
       !
       qlnp(1) = qlnp(1) + log(1000._r8) * root2
       !
       ! Spectral Coefficients of Virtual Temp (K)
       !
       DO k = kMax,1,-1
          CALL GWriteField(n, qtmp(:,k))
       END DO
       !
       ! Spectral Coefficients of Divergence and Vorticity (1/seg)
       !
       IF(TRIM(StrFormat) == 'old')THEN
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Divergence
             CALL GWriteField(n, qdiv(:,k))
             !Vorticity
             CALL GWriteField(n, qrot(:,k))
          END DO
       else if(TRIM(StrFormat) == 'new') then
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          DO k = kMax,1,-1
             !Divergence
             CALL GWriteField(n, qdiv(:,k))
          END DO

          DO k = kMax,1,-1
             !Vorticity
             CALL GWriteField(n, qrot(:,k))
          END DO
       ELSE
         CALL FatalError('Invalid StrFormat. see at the namelist, Aborting Model!')
       END IF
       !
       ! Spectral Coefficients of Specific Humidity (kg/kg)
       !
       DO k = kMax,1,-1
          CALL GWriteField(n, qq(:,k))
       END DO
       !
    ELSE
       !
       ! Spectral Coefficients of Orography (m)
       !
       CALL Collect_Spec(qgzs, aux(:,1), 1, 1, 0)
       IF(myid.eq.0) CALL GWriteField(n, aux(:,1))
       !
       ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
       !
       CALL Collect_Spec(qlnp, aux(:,1), 1, 1, 0)
       !      transform surface pressure from Pascal to cbar
       !
       IF(myid.eq.0) aux(1,1) = aux(1,1) - log(1000._r8) * root2
       IF(myid.eq.0) CALL GWriteField(n, aux(:,1))
       !
       ! Spectral Coefficients of Virtual Temp (K)
       !
       CALL Collect_Spec(qtmp, aux, kmaxloc, kmax, 0)
       IF(myid.eq.0) THEN
          DO k = kMax,1,-1
             CALL GWriteField(n, aux(:,k))
          END DO
       ENDIF
       !
       ! Spectral Coefficients of Divergence and Vorticity (1/seg)
       !
      IF(TRIM(StrFormat) == 'old')THEN
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          !Divergence
          CALL Collect_Spec(qdiv, aux, kmaxloc, kmax, 0)
          !Vorticity
          CALL Collect_Spec(qrot, aux1, kmaxloc, kmax, 0)
          IF(myid.eq.0) THEN
             DO k = kMax,1,-1
                !Divergence
                CALL GWriteField(n, aux(:,k))
                !Vorticity
                CALL GWriteField(n, aux1(:,k))
             END DO
          ENDIF
       else if(TRIM(StrFormat) == 'new') then
          !
          !Spectral Coefficients of Divergence and Vorticity
          !
          !Divergence
          CALL Collect_Spec(qdiv, aux, kmaxloc, kmax, 0)
          !Vorticity
          CALL Collect_Spec(qrot, aux1, kmaxloc, kmax, 0)
          IF(myid.eq.0) THEN
             DO k = kMax,1,-1
                !Divergence
                CALL GWriteField(n, aux(:,k))
             END DO
             DO k = kMax,1,-1
                !Vorticity
                CALL GWriteField(n, aux1(:,k))
             END DO
          ENDIF
       ELSE
         CALL FatalError('Invalid StrFormat. see at the namelist, Aborting Model!')
       END IF
       !
       ! Spectral Coefficients of Specific Humidity (kg/kg)
       !
       CALL Collect_Spec(qq, aux, kmaxloc, kmax, 0)
       IF(myid.eq.0) THEN
          DO k = kMax,1,-1
             CALL GWriteField(n, aux(:,k))
          END DO
       ENDIF
    ENDIF

    IF(nfctrl(43).GE.1)WRITE(UNIT=nfprt,FMT=3001)ifday,tod,idate,idatec,n
3001 FORMAT(' GWRITE IFDAY=',I8,' TOD=',F8.1,2(2X,3I3,I5), 2X,'N=',I2)
  END SUBROUTINE gwrite


  !hmjb
  SUBROUTINE getco2(time,co2val)
    !==========================================================================
    ! getco2: Interpolates Mauna Loa data for a given time
    !
    ! *** Atmospheric CO2 concentrations (ppmv) derived from in situ  ***
    ! *** air samples collected at Mauna Loa Observatory, Hawaii      ***
    !
    ! Data:
    !
    !   http://cdiac.ornl.gov/trends/co2/contents.htm
    !   http://cdiac.ornl.gov/ftp/trends/co2/maunaloa.co2
    !
    ! Parabolic fitting by hbarbosa@cptec.inpe.br, 17 Jan 2007:
    !
    !   co2val = a*(time-2000)^2 + b*(time-2000) + c
    !
    !       a  = 0.0116696   +/- 0.0005706    (4.89%)
    !       b  = 1.79984     +/- 0.022        (1.222%)
    !       c  = 369         +/- 0.1794       (0.04863%)
    !
    !==========================================================================
    !     time.......date of current data
    !     time(1)....hour(00/12)
    !     time(2)....month
    !     time(3)....day of month
    !     time(4)....year
    !
    !    co2val....co2val is wgne standard value in ppm "co2val = /345.0/
    !==========================================================================

    IMPLICIT NONE
    REAL(KIND=r8), PARAMETER :: A = 0.0116696
    REAL(KIND=r8), PARAMETER :: B = 1.79984
    REAL(KIND=r8), PARAMETER :: C = 369.0

    INTEGER,       INTENT(IN ) :: time(4)
    REAL(KIND=r8), INTENT(OUT) :: co2val

    REAL(KIND=r8) :: TDIF

    tdif=time(4) + (time(2)-1.)/12. + (time(3)-1.+ time(1)/24.)/365. - 2000.

    co2val = A*tdif**2 + B*tdif + C

    !    WRITE(*,123) time,tdif+2000.,co2val
    !123 format('hmjb co2val date=',3(I2,1x),I4,' fyear=',F10.5,' val=',F7.3)

    RETURN
  END SUBROUTINE getco2
  !hmjb
END MODULE InputOutput
