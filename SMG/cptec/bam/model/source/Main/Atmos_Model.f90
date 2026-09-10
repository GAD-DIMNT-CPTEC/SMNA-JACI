
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.28 $
!  $Modificacao Solange 27-01-2012: Inclue a saida de dois arquivos para o G3DVAR/ Grupo Assim_Dados (GSI)
!                                   nameg  -->>  pnew*    - GSI
!                                   namer  -->>  sfcfg*   - GSI
!       P.S.: Manteve-se a saida de namef  -->>  fct padrao CPTEC
!
! Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
! Fine aerosol mode 2000 is included
!  Modifications (2) are marked by 
!  !tar begin and  !tar end
!

MODULE AtmosModelMod
  USE Constants, ONLY:    r8,i8

  USE Parallelism, ONLY:   &
       CreateParallelism,  &
       DestroyParallelism, &
       MsgOne,             &
       FatalError,         &
       unitDump,           &
       myId

  USE Watches, ONLY:  &
       CreateWatches, &
       NameWatch,     &
       ChangeWatch,   &
       DumpWatches,   &
       DestroyWatches

!  USE Dumpgraph, ONLY: &
!       dumpgra, &
!       writectl


   IMPLICIT NONE
  SAVE       

  PRIVATE
  TYPE land_ice_atmos_boundary_type
      !variables of this type are declared by coupler_main, allocated by flux_exchange_init
      !quantities going from land+ice to atmos
      REAL, DIMENSION(:,:), POINTER :: albedo => NULL()
  END TYPE land_ice_atmos_boundary_type

  TYPE surf_diff_type
      REAL, POINTER, DIMENSION(:,:) :: dtmass  =>NULL()
  END TYPE surf_diff_type

  TYPE atmos_data_type
      REAL, POINTER, DIMENSION(:,:) :: t_bot =>NULL() 
  END TYPE atmos_data_type

  CHARACTER(LEN=200)   :: roperm
  CHARACTER(LEN=  9)   :: namee
  CHARACTER(LEN=  9)   :: nameh
  CHARACTER(LEN=  9)   :: namef
  CHARACTER(LEN=  9)   :: nameg
  CHARACTER(LEN=  9)   :: namer
  CHARACTER(LEN= 10)   :: labeli
  CHARACTER(LEN= 10)   :: labelc
  CHARACTER(LEN= 10)   :: labelf
  CHARACTER(LEN=  4)   :: PRC='    '  
  LOGICAL              :: lreststep
  LOGICAL              :: restart
  LOGICAL              :: bckhum
  LOGICAL              :: dotrac
  LOGICAL              :: dohum

  REAL(KIND=r8)   , ALLOCATABLE :: qgzs_orig(:)
  REAL(KIND=r8)   , ALLOCATABLE :: lsmk(:)
  REAL(KIND=r8)   , ALLOCATABLE :: a_in(:)
  REAL(KIND=r8)   , ALLOCATABLE :: b_in(:)

!  REAL(KIND=r8)   , ALLOCATABLE :: del_in(:)

  REAL(KIND=r8)   , ALLOCATABLE :: ct_in(:)
  REAL(KIND=r8)   , ALLOCATABLE :: cq_in(:)

!  REAL(KIND=r8)   , ALLOCATABLE :: tequi(:,:,:)

  LOGICAL         , ALLOCATABLE :: cehl(:)
  REAL(KIND=r8)   , ALLOCATABLE :: cehr(:)

  REAL(KIND=r8)   , ALLOCATABLE :: rlsm(:,:)  
  REAL(KIND=r8)   , ALLOCATABLE :: IVGTYP(:,:)

  REAL(KIND=r8)      :: fdh
  REAL(KIND=r8)      :: dth
  REAL(KIND=r8)      :: delth
  REAL(KIND=r8)      :: fdayh
  REAL(KIND=r8)      :: cthw
  INTEGER             :: maxt0
  REAL(KIND=r8)      :: fa
  REAL(KIND=r8)      :: fb
  REAL(KIND=r8)      :: fb1

  INTEGER              :: ids(4)
  INTEGER              :: idc(4)

  INTEGER              :: ifday
  REAL(KIND=r8)        :: tod
  LOGICAL              :: enhdifl
  INTEGER              :: iovmax
  INTEGER              :: nsca_save
  INTEGER              :: ifdy=0
  REAL(KIND=r8)       :: todcld=0.0_r8
  REAL(KIND=r8)       :: todsib
  INTEGER              :: limlow
  REAL(KIND=r8)       :: delta2

  INTEGER              :: ngra
  LOGICAL, PARAMETER :: instrument=.TRUE.

  PUBLIC :: atmos_model_init
  PUBLIC :: atmos_model_run
  PUBLIC :: atmos_model_finalize

CONTAINS
  SUBROUTINE atmos_model_init(myid_in,maxnodes_in)

  USE  Options, ONLY: &
       SetOutPut,SetOutPutDHN,ReadNameList, DumpOptions,CreateFileName,FNameRestInput2,&
       FNameRestOutput2,FNameSibPrgInp0,FNameConvClInp0,&
       FNameGDHN ,FNameGDYN,FNameGPRC,FNamenDrGH,FNameTopGH,FNameOutGH,&
       maxtim, trunc, vert,dt,initlz,idate,idatec,idatef,delt,&
       maxtid,nstep,istrt,ifilt,filta,filtb,tk,dk, &
       jdt,ddelt,ifsst,ifco2flx,intsst,intflxco2,sstlag,ktm,kt,yrl,ktp,&
       ifslm,intsoilm,soilmlag,Flxco2lag,&
       maxtfm,dctd,mdxtfm,dcte,cteh0,mextfm,nfctrl,cdhl,ctdh0,monl,&
       isimp,doprec,dogwd,enhdif,igwd,grhflg,allghf,dodyn,&
       start,reducedGrid,linearGrid,slagr,nlnminit,SL_twotime_scheme,slhum,nfsibd,&
       grid_difus,eigeninit,diabatic,&
       nfprt, nfin0, nfin1, nfout0, nfout1, nfsibo, nfsibi, nfdrct, &
       nfin1,nfprt,nfsibi,nfsibo,nfin0,nfdrct,nfout1,&
       nfcnv0,nfdiag,nffcst,nftmp,nedrct,ndhndrct,neprog,nefcst,ndhnfcst,nfcnv1,&
       nhdhn,nfdhn,nfdyn,nfprc,nfghdr,nfghds,nfghloc,nfghtop,nfghou,&
       fNameInput1,mgiven,gaussgiven,fNameInput0,fNameNmi,&
       fNameIBISDeltaTemp,fNameSandMask,fNameClayMask,fNameClimaTemp,fNameSSTAOI,&
       fNameCnfTbl,fNameCnf2Tb,fNameSibVeg,fNameIBISMask,fNameLookTb,fNameUnitTb, &
       fNameDTable, fNameSibAlb,fNameSibmsk ,fNameMicro,fNameSoilms,fNameCO2FLX ,&
       iqdif,path_in,dirfNameOutput,PREFX ,EXDW,EXTW,UNIFIED,&
       cthl,schemes,nfsibt,LV,TRCG,nscalars,microphys ,fNameSlabOcen, nfsoiltp, &
       record_type,fNameSoilType,GenAssFiles,nClass,nAeros,co2val,isimco2,ifozone,typechem,indexchem

  USE InputOutput, ONLY:       &
       InitInputOutput       , &
       gread4                , &
       gwrite                , &
       fsbc

  USE IOLowLevel, ONLY:    &
       InitReadWriteSpec         , &
       ReadHead                  , &
       GReadHead                 , &
       !LandSoilmMask             , &
       LandFlxCO2Mask            ,& 
       ReadLandSeaMask2          , &
       LandSeaMask

  USE Utils, ONLY:     &
       rcl,            &
       colrad,         &
       InitTimeStamp,  &
       total_mass,&
       total_flux,nTtimes ,aTfluxco2,totflux

  USE Diagnostics, ONLY:      &
       InitDiagnostics       , &
       StartStorDiag         , &
       rsdiag                , &
       Prec_Diag             , &
       accpf                 , &
       wridia                , &
       wdhnprog              , &
       weprog                , &
       wrprog

  USE FieldsDynamics, ONLY: &
       fgq, fgqm, fgqmm, fgum, fgvm, fgtmpm, fgdivm, &
       fgtlamm, fgtphim, fglnpm, fgplamm, fgpphim, &
       qgzs, qlnpp, qlnp_nabla4,qtmpp, qdivp, qrotp, qqp,fgtmp,fgzs, &
       fgpass_scalars, adr_scalars, fgice, fgicem, fgliq, fgliqm,omg, fgvar, fgvarm,&
       fgprsl,fgprsi,fgphil,fgphii
       
  USE Sizes, ONLY:            &
       myMNMax               , &
       imaxperj              , &
       ibMax                 , &
       jbMax                 , &
       mnMax                 , &
       mMax                  , &
       nMax                  , &
       mnMap                 , &
       kMax                  , &
       imax                  , &
       jmax                  , &
       ijMax                 , &
       ijMaxGauQua           , &
       ibPerIJ               , &
       jbPerIJ               , &
       ibMaxPerJB            , &
       a_hybr                , &
       b_hybr                , &
! add solange 21-01-2012
       havesurf              , &
       kMaxloc
! fim add

  USE FieldsPhysics, ONLY:     &
       InitFieldsPhyscs,  &
       capac0,w0 , prct, prcc,geshem ,gtsea,td0,&
       iMask,SoilMask,lowlyr,ustar,z0,temp2m,umes2m,tracermix,&
! add solange 27-01-2012
       o3mix,ustar,&
       sheleg,tg0,tc0,convc,convt,convb,&
       sm0,mlsi,topoi,   &
       AlbNirBeam, AlbNirDiff, AlbVisBeam, AlbVisDiff,sfc,std,statec

       !avisb,avisd,anirb,anird
! fim add
      
  USE SemiLagrangian, ONLY: &
       InitSL , &
       ulonm  , &
       ulatm  , &
       uetam  , &
       ulonm2D, &
       ulatm2D

  USE Griddynamics, ONLY:      &
       SetJablo,               &
       do_globconserv,         &
       do_globfluxconserv,     &
       init_globconserv,       &
       init_globfluxconserv

  USE ModTimeStep, ONLY:      &
       SfcGeoTrans,           &
       InitTrans,             &
       TimeStep   ,           &
       InitBoundSimpPhys,     &
       GetSfcTemp

  USE GridHistory, ONLY:      &
       InitGridHistory       , &
       WriteGridHistory      , &
       TurnOnGridHistory     , &
       IsGridHistoryOn       , &
       WriteGridHistoryTopo

  USE NonLinearNMI, ONLY:  &
       Nlnmi                     , &
       Diaten, &
       Getmod

  USE Surface, ONLY:          &
       InitSurface

  USE SFC_SSiB, ONLY:          &
       InitCheckSSiBFile,&
       InitSurfTemp,&
       ReStartSSiB

  USE SFC_SiB2, ONLY:          &
       InitCheckSiB2File,&
       InitSurfTempSiB2,&
       ReStartSiB2

  USE Sfc_Ibis_Fiels, ONLY:          &
       ReStartIBIS,&
       vegtype0
       
  USE GwddDriver, ONLY:   &
       InitGWDDDriver

  USE Init, ONLY :            &
       InitAll, nls

  USE PblDriver, ONLY:   &
       InitPBLDriver

  USE SfcPBLDriver, ONLY: &
       InitSfcPBL_Driver

  USE Convection, ONLY:       &
       InitConvection,&
       InitCheckFileConvec,&
       ReStartConvec

  USE ModRadiationDriver, ONLY:        &
       InitRadiationDriver         

  USE PhysicsDriver , ONLY: &
       InitSimpPhys

  USE SpecDynamics, ONLY:      &
       bmcm

  USE PhysicalFunctions, ONLY: &
      InitPhysicalFunctions     

! add solange 27-01-2012
  USE SpecDump , ONLY:InitSpecDump, write_sigma_file
  USE GridDump , ONLY:InitGridDump, write_GridSigma_file
!fim add  
      


    INCLUDE 'mpif.h'

    INTEGER,    INTENT(IN), OPTIONAL :: myid_in
    INTEGER,    INTENT(IN), OPTIONAL :: maxnodes_in

    CHARACTER (LEN=10)   :: DateInit_s
    REAL(KIND=r8)                 :: SumDel

    CHARACTER(LEN=*), PARAMETER :: h="**(atmos_model_init)**" 
    INTEGER              :: nPtland
    REAL(r8), PARAMETER :: amdc  = 0.658114_r8! Molecular weight of dry air / carbon dioxide
    REAL(KIND=r8)                 :: zero=0.0_r8
    REAL(KIND=r8)                 :: ahour
    INTEGER          :: ierr
    INTEGER          :: i
    INTEGER          :: j 
    INTEGER          :: k, km
    INTEGER          :: l,l1,l2,m
    INTEGER          :: ij,itr
    IF (PRESENT(myid_in) .AND. PRESENT(maxnodes_in) ) THEN
       CALL MsgOne(h," atmos_model_init")
    END IF

  ! engage MPI 

  CALL CreateParallelism()

  ! execution time instrumentation

  IF (instrument) THEN
     CALL CreateWatches(12, 1)
     CALL NameWatch(1,"Initialize     ")
     CALL NameWatch(2,"Integrate      ")
     CALL NameWatch(3,"Outros         ")
     CALL NameWatch(4,"Backtrans      ")
     CALL NameWatch(5,"GRPcomp        ")
     CALL NameWatch(6,"Slagr          ")
     CALL NameWatch(7,"Dirtrans       ")
     CALL NameWatch(8,"Semiimplicit   ")
     CALL NameWatch(9,"Humidphys      ")
     CALL NameWatch(10,"HumidDirtrans  ")
     CALL NameWatch(11,"Humidbacktrans  ")
     CALL NameWatch(12,"ShortWave  ")
  END IF

  ! read name list and fill all options

  CALL ReadNameList()

  fsbc    =.TRUE.

  ALLOCATE(cehl(0:maxtid))
  ALLOCATE(cehr(1:maxtid))
  !ALLOCATE(cdhr(1:maxtid))
  !
  !     Get Initial Labels and Vertical Discretization
  !
  ALLOCATE (a_in (vert+1))
  ALLOCATE (b_in (vert+1))
  !ALLOCATE (del_in(vert  ))
  ALLOCATE (ct_in (vert  ))
  ALLOCATE (cq_in (vert  ))

  CALL MsgOne(h," ALLocation done ")
  IF (myid < 100) WRITE(PRC,'(a1,i3.3)')'P',myid
  IF (myid >= 100 .AND. myid < 1000) WRITE(PRC,'(a1,i3)')'P',myid

  CALL CreateFileName()

  IF (TRIM(start) == "warm" ) THEN
     OPEN(UNIT=nfin1, FILE=TRIM(fNameInput1)//TRIM(PRC), FORM='unformatted', ACCESS='sequential', &
          ACTION='read', STATUS='old',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(fNameInput1), ierr
        STOP "**(ERROR)**"
     END IF
     if(myid.eq.0) write(*,*) 'will call GReadHead ',nfin1
     CALL GReadHead (nfin1, ifday, tod, idate, idatec, a_in, b_in, vert)
     CALL InitTimeStamp (DateInit_s, idate)
     labeli  = DateInit_s
     CALL InitTimeStamp (DateInit_s, idatec)
     labelc  = DateInit_s
     WRITE(labelf,'(I4.4, 3I2.2)' ) (idatef(i),i=4,1,-1)
     nfcnv0 = nfcnv1
     nfsibi = nfsibo
     CLOSE(UNIT=nfin1)
  ELSE

     OPEN(UNIT=nfin1, FILE=TRIM(fNameInput1), FORM='unformatted',  ACCESS='sequential',&
          ACTION='read', STATUS='old', IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(fNameInput1), ierr
        STOP "**(ERROR)**"
     END IF

     CALL ReadHead (nfin1, ifday, tod, idate, idatec, a_in, b_in, vert)
     CALL InitTimeStamp (DateInit_s, idate)
     labeli  = DateInit_s
     WRITE(labelf,'(I4.4, 3I2.2)' ) (idatef(i),i=4,1,-1)
     labelc  = labelf
     idatec  = idate
     IF ((initlz < 0))THEN
        nfcnv0 = nfcnv1
        nfsibi = nfsibo
     END IF
     REWIND (nfin1)
  END IF
  !
  ! initialize modules
  !
  IF(TRIM(isimp).EQ.'YES'.OR.slhum) iqdif='NO'
  if (myid.eq.0) write(*,*) 'calling initall '
  CALL InitAll(trunc, vert, reducedGrid, linearGrid, mgiven, gaussgiven, &
       ct_in, cq_in, a_in, b_in, dk, tk ,dt)
  CALL MsgOne(h," After initall   ")

  IF( TRIM(start) == "warm" )THEN
     fNameInput0=TRIM(fNameInput0)//TRIM(PRC)
     fNameInput1=TRIM(fNameInput1)//TRIM(PRC)
  END IF






  CALL InitPhysicalFunctions()
     CALL MsgOne(h," After InitPhysicalFunctions  ")

  CALL InitInputOutput (mMax, nMax , mnMax, kmax, &
       path_in, fNameCnfTbl, &
       fNameCnf2Tb, fNameLookTb, fNameUnitTb)

     CALL MsgOne(h," After InitInputOutput    ")
  CALL InitDiagnostics (doprec, dodyn  ,    colrad , &
       mMax   ,nMax   ,mnMax  , iMax   ,   jMax, &
       kMax   ,ibMax  ,jbMax  , ibMaxPerJB, fNameDTable)

     CALL MsgOne(h," After InitDiagnostics    ")
  CALL InitReadWriteSpec(&
       mMax        ,mnMax   ,kMax      ,ijMaxGauQua ,iMax    , &
       jMax        ,ibMax   ,jbMax     )

     CALL MsgOne(h," After InitReadWriteSpec  ")
  CALL InitFieldsPhyscs(ibMax, kMax, jbMax,nClass,nAeros)

     CALL MsgOne(h," After InitFieldsPhyscs   ")
  DO jdt=0,maxtim
     cehl(jdt)=.FALSE.
     cdhl(jdt)=.FALSE.
     cthl(jdt)=.FALSE.
  END DO
  !
  dogwd=1
  delt     = ddelt
!PK  IF (nstep.EQ.1) nstep=7

  IF(TRIM(enhdif).EQ.'YES') THEN
     enhdifl = .TRUE.
  ELSE
     enhdifl = .FALSE.
  ENDIF

     CALL MsgOne(h," After TRIM(enhdif)       ")
  IF(TRIM(isimp).EQ.'YES') THEN
     WRITE(UNIT=nfprt,FMT=7)
     initlz=0
     nfcnv0=0
     microphys=.FALSE.
     igwd  ='NO'
  END IF
  !
  ! Initialize modules
  !
  IF (slagr.or.slhum) THEN
     CALL MsgOne(h," calling InitSL           ")
     CALL InitSL
     CALL MsgOne(h," after   InitSL           ")
     iovmax = 2
     IF (nstep.EQ.1) nstep=7
   ELSE
     iovmax = 0
     IF (nstep.EQ.1) nstep=7
  END IF
  !
  ! prepare output files
  !
  roperm  = dirfNameOutput
  namee   = "GPRG"//TRIM(PREFX)
  nameh   = "GDHN"//TRIM(PREFX)
  namef   = "GFCT"//TRIM(PREFX)
! add solange 27-01-2012
  nameg   = "GASS"//TRIM(PREFX)
  namer   = "GSFC"//TRIM(PREFX)
! fim add

  !
  ! Criate restart file name   
  !  
  IF(initlz < 0)THEN 
     CALL CreateFileName(TRIM(PRC),1,'RESTAT_SIB')
  ELSE
     CALL CreateFileName(TRIM(PRC),1)
  END IF
  IF( TRIM(start) == "warm" )THEN
     nfcnv0 = nfcnv1
     nfsibi = nfsibo
     OPEN (UNIT=nfsibi, FILE=TRIM(FNameSibPrgInp0),FORM='unformatted',ACCESS='sequential',&
          ACTION='read',STATUS='old',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(FNameSibPrgInp0), ierr
        STOP "**(ERROR)**"
     END IF
     OPEN (UNIT=nfcnv0, FILE=TRIM(FNameConvClInp0),FORM='unformatted',ACCESS='sequential',&
          ACTION='read',STATUS='old',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(FNameConvClInp0), ierr
        STOP "**(ERROR)**"
     END IF
  ELSE IF ( TRIM(start) == "cold" .AND. (initlz < 0))THEN
      IF(initlz == -3)THEN 
        nfcnv0 = nfcnv1
        nfsibi = nfsibo
      ELSE IF(initlz == -2)THEN 
        nfsibi = nfsibo
      ELSE
        STOP "**(ERROR)**"      
      END IF 
     OPEN (UNIT=nfsibi, FILE=TRIM(FNameSibPrgInp0),FORM='unformatted',ACCESS='sequential',&
          ACTION='read',STATUS='old',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(FNameSibPrgInp0), ierr
        STOP "**(ERROR)**"
     END IF
     IF(nfcnv0 /= 0) THEN
        OPEN (UNIT=nfcnv0, FILE=TRIM(FNameConvClInp0),FORM='unformatted',ACCESS='sequential',&
             ACTION='read',STATUS='old',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                TRIM(FNameConvClInp0), ierr
           STOP "**(ERROR)**"
        END IF
     END IF
  ENDIF
  RESTART = nfin0 .NE. nfin1
  if (myid.eq.0) write (*,*) ' defining restart  ', RESTART
  !
  !     start reading initial values (t=t   )
  !
  IF (nfin0 == nfin1) THEN
      !
      ! cold initialization
      !
      CALL gread4 (nfin0, ifday, tod, idate, idatec, qgzs, qlnpp, &
                   qtmpp, qdivp, qrotp, qqp, a_hybr, b_hybr, dodyn, &
                   nfdyn)
      REWIND nfin0
      !$OMP PARALLEL
      CALL GetSfcTemp ()
      !$OMP END PARALLEL
  END IF

  ALLOCATE(lsmk(ijMaxGauQua))
  ALLOCATE(rlsm(imax,jmax))
  ALLOCATE(IVGTYP(iBMax,jBMax))
  IF (ifco2flx .GE. 4) THEN
    CALL LandFlxCO2Mask(ifco2flx,labeli,intflxco2,Flxco2lag,fNameCO2FLX,rlsm)
  END IF
  IF (ifsst .GE. 4) THEN
     CALL LandSeaMask(ifsst,labeli,intsst,sstlag,fNameSSTAOI,rlsm)
     ij=0
     DO j=1,jmax
        DO i=1,imax
           ij=ij+1
           lsmk(ij)=rlsm(i,j)
        END DO
     END DO
  ELSE
     CALL ReadLandSeaMask2 (TRIM(fNameSSTAOI), lsmk)
  END IF
!  IF (ifslm .GE. 4) THEN
!     CALL LandSoilmMask(ifslm,labeli,intsoilm,soilmlag,fNameSoilms,rlsm)
!  END IF
! IF(GenAssFiles)THEN
!     CALL InitSpecDump(si,sl,idatec,idate)
!     CALL InitGridDump(si,sl,idatec,reducedGrid,iMax,&
!                  jMax,kMax,ibMax,jbMax,&
!                  record_type,nfsoiltp,fNameSoilType,nfprt,ibMaxPerJB)
! END IF
!  IF(GenAssFiles)THEN
!      CALL InitSpecDump(si,sl,idatec,idate)
!      CALL InitGridDump(si,sl,idatec,reducedGrid,iMax,&
!                   jMax,kMax,ibMax,jbMax,&
!                   record_type,nfsoiltp,fNameSoilType,nfprt,ibMaxPerJB)
!  END IF
  CALL InitPBLDriver(ibMax,jbMax,kmax,a_hybr,b_hybr,RESTART)
  call msgone(h, " after   InitPBLDriver ")
  CALL InitGWDDDriver(ibMax,jbMax,iMax,jMax,kmax,ibMaxPerJB)
  call msgone(h, " after  InitGWDDDrivey ")
  CALL InitConvection(std,a_hybr ,b_hybr     , &
                     kmax    ,iMax       ,jMax       ,ibMax,jbMax,&
                     trunc   ,ifdy       ,todcld     ,ids        , &
                     idc     ,ifday      ,tod        ,fNameMicro,path_in,  idate     )
  call msgone(h, " after  InitCOnvection ")
  CALL InitSurface(ibMax             ,jbMax             ,iMax              ,jMax          , &
                   kMax              ,path_in           ,fNameSibVeg       , &
                   fNameSibAlb       ,idate             ,idatec            ,dt            , &
                   nfsibd            ,nfprt             ,nfsibt            ,fNameSibmsk   , &
                   ifday             ,ibMaxPerJB        ,tod               ,ids           , &
                   idc               ,ifdy              ,todsib            ,fNameIBISMask , &
                   fNameIBISDeltaTemp,fNameSandMask     ,fNameClayMask     ,fNameClimaTemp, &
                   RESTART           ,imask             ,gtsea             ,fgtmp(:,kmax:1:-1,:), &
                   fgq(:,kmax:1:-1,:),topoi             ,fNameSlabOcen     )
  call msgone(h, " after  Initsurface    ")
  CALL InitRadiationDriver(monl,yrl,kmax,a_hybr,b_hybr,dt,nls)
  call msgone(h, " after  Initradiation  ")
  !
  ! Write problem options to stdio
  !
  CALL DumpOptions()
!  IF (nhdhn /= 0) THEN
!      !
!      ! Criate DHN file name   
!      !  
!     CALL CreateFileName(TRIM(PRC),"GDHN")
!     OPEN (UNIT=nfdhn, FILE=TRIM(FNameGDHN), FORM='unformatted', ACCESS='sequential',&
!          ACTION='write', STATUS='unknown',IOSTAT=ierr)
!     IF (ierr /= 0) THEN
!        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
!             TRIM(FNameGDHN), ierr
!        STOP "**(ERROR)**"
!     END IF
!  END IF

  IF (dodyn) THEN
      !
      ! Criate Dynamics file name   
      !  
     CALL CreateFileName(TRIM(PRC),"GDYN")
     OPEN (UNIT=nfdyn, FILE=TRIM(FNameGDYN),FORM='unformatted', ACCESS='sequential',&
          ACTION='write', STATUS='replace',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(FNameGDYN), ierr
        STOP "**(ERROR)**"
     END IF
  END IF

  IF (doprec) THEN
     CALL CreateFileName(TRIM(PRC),"GPRC")
     OPEN (UNIT=nfprc, FILE=TRIM(FNameGPRC), FORM='unformatted', ACCESS='sequential',&
          ACTION='write', STATUS='replace',IOSTAT=ierr)
     IF (ierr /= 0) THEN
        WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
             TRIM(FNameGPRC), ierr
        STOP "**(ERROR)**"
     END IF
  END IF

  IF (grhflg) THEN
     CALL CreateFileName(TRIM(PRC),"GFGH")
     IF(myid ==0 ) THEN
        OPEN (UNIT=nfghdr, FILE=TRIM(FNamenDrGH), FORM='formatted', ACCESS='sequential',&
             ACTION='write',STATUS='replace',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                TRIM(FNamenDrGH), ierr
           STOP "**(ERROR)**"
        END IF
     END IF
  END IF
  !
  !     read cloud dataset - logic assumes that initialization not
  !     performed for warm start
  !
  jdt=0
  CALL InitGridHistory (idate,idatec ,iovmax, allghf,grhflg,nfghds  , &
       nfghloc   ,nfghdr ,iMax  ,jMax  ,ibMax ,jbMax ,ibPerIJ, &
       jbPerIJ,kMax,ibMaxPerJB,iMaxPerJ,start)
























  IF (nscalars.gt.0) THEN
        DO k=1,kMax
           DO j=1,jbmax
              DO i=1,ibmax
                 fgpass_scalars(i,k,j,1,1) = tracermix(i,k,j)
                 fgpass_scalars(i,k,j,1,2) = tracermix(i,k,j)
                 fgpass_scalars(i,k,j,2,1) = tracermix(i,k,j)*0.000000000008_r8
                 fgpass_scalars(i,k,j,2,2) = tracermix(i,k,j)*0.000000000008_r8
             END DO
          END DO
       END DO
  END IF
  IF (microphys) THEN
        DO k=1,kMax
           DO j=1,jbmax
              DO i=1,ibmax
                 fgicem(i,k,j) = tracermix(i,k,j)*0.000000000008_r8
                 fgice (i,k,j) = tracermix(i,k,j)*0.000000000008_r8 
                 fgliqm(i,k,j) = tracermix(i,k,j)
                 fgliq (i,k,j) = tracermix(i,k,j) 
             END DO
          END DO
       END DO

       IF((nClass+nAeros)>0)THEN
          DO m=nClass+1,nClass+nAeros
                IF(TRIM(typechem(m))=='CO2')THEN
                   DO k=1,kmax
                      DO j=1,jbmax
                         DO i=1,ibmax
                            fgvarm  (i,k,j,indexchem(m))=co2val*1.0e-6_r8 / amdc  !convert kg/kg to mol/mol
                         END DO
                      END DO
                   END DO
                END IF
          END DO
       END IF

       IF(ifozone /= 0)THEN
          DO m=nClass+1,nClass+nAeros
             IF(TRIM(typechem(m))=='O3')THEN
                DO k=1,kmax
                   DO j=1,jbmax
                      DO i=1,ibmax
                         fgvarm  (i,k,j,indexchem(m))=o3mix(i,k,j)   !convert kg/kg to mol/mol
                         !fgvar   (i,k,j,indexchem(m))=o3mix(i,k,j)   !convert kg/kg to mol/mol
                      END DO
                   END DO
                END DO
             END IF
          END DO
       END IF
  END IF

  LOWLYR=1
  IVGTYP=REAL(imask,kind=r8)
  !
  CALL InitSfcPBL_Driver(RESTART,ibMax  ,jbMax  ,USTAR  , &!(INOUT)
                        LOWLYR)


  !
  tod=0.0_r8
  ! passive scalars are not used in normal mode initialization
  !
  do_globconserv = .FALSE.
  do_globfluxconserv =.FALSE.
  init_globconserv = .FALSE.
  init_globfluxconserv= .FALSE.
  !
  IF (instrument) THEN
     CALL ChangeWatch(3)
  END IF
  IF(iabs(initlz) >=  1 .AND. initlz >= -3) THEN
     CALL MsgOne(h," Init: Diabatic Tendencies")
     CALL rsdiag
     CALL Diaten(slagr,fNameInput0,ifday, tod, idate, idatec)
     ktm=0
     ktp=0
     kt =0
     !
     !     snow reinitialization after surface temperature initialization
     !
     IF(initlz >= 1)THEN
        IF(schemes == 1) THEN
           CALL InitSurfTemp (jbMax,ibMaxPerJB )
        ELSE IF(schemes == 2) THEN
           CALL     InitSurfTempSiB2 (jbMax,ibMaxPerJB )
        ELSE IF(schemes == 3) THEN   
          !STOP 'vazio'
        END IF  
     END IF
     !
     !   reset old time step values on the grid to zero
     !   ----------------------------------------------
     fgqm = 0.0_r8
     fgqmm = 0.0_r8
     fgum = 0.0_r8
     fgvm = 0.0_r8
     fgtmpm = 0.0_r8
     fgdivm = 0.0_r8
     fgtlamm = 0.0_r8
     fgtphim = 0.0_r8
     fglnpm = 0.0_r8
     fgplamm = 0.0_r8
     fgpphim = 0.0_r8
     IF (microphys) THEN
        fgicem = 0.0_r8
        !fgice  = 0.0_r8
        fgliqm = 0.0_r8
        !fgliq  = 0.0_r8
        IF((nClass+nAeros)>0)THEN
           fgvarm= 0.0_r8
           DO itr=nClass+1,nClass+nAeros
              IF(TRIM(typechem(itr))=='CO2')THEN
                 DO j=1,jbMax
                    DO k=1,kMax
                       DO i=1,ibMaxPerJB(j)
                          ! convert mol/mol   to   kg/kg 
                          !kg/kg  =            mol/mol
                          fgvarm(i,k,j,indexchem(itr)) = co2val*1.0e-6_r8 / amdc
                       END DO
                    END DO   
                 END DO
              END IF

              IF(ifozone /= 0)THEN
                 IF(TRIM(typechem(itr))=='O3')THEN
                    DO j=1,jbmax
                       DO k=1,kmax
                          DO i=1,ibMaxPerJB(j)
                             fgvarm  (i,k,j,indexchem(itr))=o3mix(i,k,j)   !convert kg/kg to mol/mol
                          END DO
                      END DO
                   END DO
                END IF
             END IF


           END DO
        END IF
     END IF
  END IF

  ALLOCATE(qgzs_orig(2*myMNMax))

  CLOSE(UNIT=nfin1)
  OPEN(UNIT=nfin0, FILE=TRIM(fNameInput0), FORM='unformatted',ACCESS='sequential',&
       ACTION='read', STATUS='old', IOSTAT=ierr)
  IF (ierr /= 0) THEN
     WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
          TRIM(fNameInput0), ierr
     STOP "**(ERROR)**"
  END IF
  OPEN(UNIT=nfin1, FILE=TRIM(fNameInput1), FORM='unformatted',ACCESS='sequential', &
       ACTION='read', STATUS='old',  IOSTAT=ierr)
  IF (ierr /= 0) THEN
     WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
          TRIM(fNameInput1), ierr
     STOP "**(ERROR)**"
  END IF
  !
  !     start reading initial values (t=t   )
  !
  IF (nfin0 .EQ. nfin1) THEN
     !
     ! cold inicialization
     !
     restart = .FALSE.
     CALL gread4 (nfin0, ifday, tod  , idate, idatec, qgzs  , qlnpp, &
                  qtmpp, qdivp, qrotp, qqp  , a_hybr, b_hybr, dodyn, &
                  nfdyn)
!!   Jablo
!    qlnpp = 0.
!    qtmpp = 0.
!    qdivp = 0.
!    qrotp = 0.
!    qgzs  = 0.
!    qqp   = 0.
!!   Jablo
     qgzs_orig = qgzs
     REWIND nfin0
  ELSE
     !
     ! warm inicialization
     !
     restart = .TRUE.
     READ(UNIT=nfin0) fgqm, fgum, fgvm, fgtmpm, fgdivm, fglnpm, fgtlamm, &
          fgtphim, fgplamm, fgpphim
     IF (slhum) READ(UNIT=nfin0) fgq
     IF(UNIFIED)READ(UNIT=nfin0) fgprsl,fgprsi,fgphil,fgphii
     IF (microphys) THEN
        IF((nClass+nAeros)>0)THEN
           READ(UNIT=nfin0) fgicem, fgice, fgliqm, fgliq,omg, fgvarm,fgvar
        ELSE
           READ(UNIT=nfin0) fgicem, fgice, fgliqm, fgliq,omg
        END IF
     END IF  
     IF (slagr.or.slhum) THEN
        READ(UNIT=nfin0) ulonm ,ulatm ,uetam,ulonm2D,ulatm2D
     END IF
     if (nscalars>0) then
        read(UNIT=nfin0) fgpass_scalars, adr_scalars
     endif




     READ(UNIT=nfin1) ifday,tod,idate,idatec,a_hybr,b_hybr
     READ(UNIT=nfin1) qgzs,qlnpp,qlnp_nabla4,qtmpp,qdivp,qrotp, &
                      qqp,total_mass,total_flux,nTtimes ,aTfluxco2,totflux
     qgzs_orig = qgzs
     REWIND nfin0
     REWIND nfin1
  END IF
  !
  IF (TRIM(isimp) == 'YES') THEN
     IF (TRIM(start) == "cold") THEN
        !$OMP PARALLEL
        CALL InitBoundSimpPhys ()
        !$OMP END PARALLEL
     END IF
    !CALL InitSimpPhys(fgtmp,tequi,sl,dt)
  END IF
  !
  IF(ifsst.GT.3.AND.sstlag.LT.zero)THEN
     WRITE(UNIT=nfprt,FMT=336)ifsst,sstlag
     STOP 336
  END IF
  !
  !
  !     write diagnostics/prognostics directory
  !
  !
  !     write uninitialized initial condition prognostics on tape
  !     the use of swrk in wrprog destroys wsib
  !
  IF(ifday.EQ.0.AND.tod.EQ.zero) THEN
     
      IF(schemes==3)THEN
         IVGTYP=vegtype0
      ELSE
         DO j=1,jbMax
            nPtland=0
            DO i=1,ibMaxPerJB(j)
               IF(imask(i,j) >=1_i8)THEN
                  nPtland=nPtland+1
                  IVGTYP(nPtland,j)=imask(i,j)
               END IF  
            END DO
         END DO   
      END IF

     CALL wrprog ( &
          nfdrct     ,nfdiag   ,ifday  ,tod    ,idate     ,idatec   , &
          qrotp      ,qdivp    ,qqp    ,qlnpp  ,qtmpp     ,gtsea    , &
          td0        ,SoilMask ,capac0 ,w0     ,imask     ,IVGTYP   , &
          temp2m     ,umes2m   ,nffcst ,nftmp  ,a_hybr    , b_hybr  ,qgzs_orig, &
          lsmk       ,tg0      ,sheleg ,mlsi   , &
          roperm     ,namef    ,labeli ,labelf ,extw      ,exdw     , &
          TRIM(TRCG) ,TRIM(LV)  ,.FALSE. &
          )

! add solange 27-01-2012
!    IF(GenAssFiles)THEN
!       CALL write_sigma_file ( &
!           ifday      ,tod       ,idate     ,idatec    , &
!           qrotp      ,qdivp     ,qqp       ,qlnpp     , &
!           qtmpp      ,qgzs_orig ,o3mix     ,del_in    , &
!           ijMaxGauQua,imax      ,jmax      ,ibMax     , &
!           jbMax      ,kMax      ,roperm    ,nameg     , &
!           labeli     ,labelf    ,extw      ,exdw      , &
!           TRIM(TRCG) ,TRIM(LV)  ,si        ,sl        , &
!           havesurf   ,kmaxloc)
!       CALL write_GridSigma_file ( &
!           ifday      ,tod       ,idate     ,idatec    ,&
!           td0        ,tg0       ,tc0       ,z0        ,&
!           convc      ,convt     ,convb     ,gtsea     ,&
!           AlbVisBeam ,AlbVisDiff,AlbNirBeam,AlbNirDiff,&
!           ustar      ,sm0       ,sheleg    ,lsmk      ,&
!           imask      ,mlsi      ,del_in    ,roperm    ,&
!           namer      ,labeli    ,labelf    ,extw      ,&
!           exdw       ,TRIM(TRCG),TRIM(LV)  ,nfprt)
!    END IF
! fim add

     IF (nhdhn .GT. 0) &
          CALL wdhnprog ( &
          ndhndrct  ,nfdhn     ,ndhnfcst  ,ifday     ,tod       ,idate    ,&
          idatec    ,qgzs_orig ,lsmk      ,qlnpp     ,qdivp    ,&
          qrotp     ,qqp       ,qtmpp     ,gtsea     ,td0      ,&
          capac0    ,w0        ,imask     ,IVGTYP    ,temp2m    ,umes2m   ,&
          roperm    ,nameh     ,labeli    ,labelf    ,a_hybr    , b_hybr  , &
          extw      ,exdw      ,TRIM(TRCG),TRIM(LV) )

     IF (mextfm .GT. 0) &
          CALL weprog ( &
          nedrct    ,neprog    ,nefcst    ,ifday     ,tod       ,idate     ,&
          idatec    ,a_hybr    ,b_hybr    ,qgzs_orig ,lsmk      ,qlnpp     ,qdivp    ,&
          qrotp     ,qqp       ,qtmpp     , gtsea     ,td0      ,SoilMask  ,&
          capac0    ,w0        ,imask     ,IVGTYP    ,temp2m    ,umes2m    ,&
          roperm    ,namee     ,labeli    ,labelf    , &
          extw      ,exdw      ,TRIM(TRCG),TRIM(LV) )
  END IF
  !
  ! compute the spectral coefficients of Laplacian of topography
  !
  CALL SfcGeoTrans(slagr)

  CALL WriteGridHistoryTopo (fgzs,FNameTopGH,nfghtop)

  !
  ! compute normal modes if necessary
  !
  IF (eigeninit) THEN
     IF (myid.eq.0) CALL Getmod(fNameNmi)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  END IF
  !
  ! Non-linear normal mode initialization
  !
  IF (nlnminit) THEN
      !
      !     do machenhauer's non-linear normal mode initialization
      !
      IF(iabs(initlz) == 2 .or. initlz == -3.or. initlz == 0) THEN
!     IF(iabs(initlz) == 2 .or. initlz == -3) THEN
        CALL MsgOne(h," Init: Non-linear Normal Modes")
        CALL rsdiag
        CALL MsgOne(h," After rsdiag   ")

        CALL Nlnmi(nlnminit,diabatic,.FALSE.,fNameNmi, ifday, tod, idatec, ktm)
        IF(schemes==3)THEN
           IVGTYP=vegtype0
        ELSE
           DO j=1,jbMax
              nPtland=0
              DO i=1,ibMaxPerJB(j)
                 IF(imask(i,j) >=1_i8)THEN
                    nPtland=nPtland+1
                    IVGTYP(nPtland,j)=imask(i,j)
                 END IF  
              END DO
           END DO   
        END IF

        CALL wrprog ( &
             nfdrct     ,nfdiag   ,ifday  ,tod    ,idate     ,idatec    , &
             qrotp      ,qdivp    ,qqp    ,qlnpp  ,qtmpp     ,gtsea     , &
             td0        ,SoilMask ,capac0 ,w0     ,imask     ,IVGTYP    , &
             temp2m     ,umes2m   ,nffcst ,nftmp  ,a_hybr    ,b_hybr    ,qgzs_orig, &
             lsmk       ,tg0      ,sheleg ,mlsi   , &
             roperm     ,namef    ,labeli ,labelf ,extw      ,exdw      , &
             TRIM(TRCG) ,TRIM(LV)  ,.FALSE.)

! add solange 27-01-2012
!       IF(GenAssFiles)THEN
!          CALL write_sigma_file ( &
!               ifday      ,tod       ,idate     ,idatec     , &
!               qrotp      ,qdivp     ,qqp       ,qlnpp      , &
!               qtmpp      ,qgzs_orig ,o3mix     ,del_in     , &
!               ijMaxGauQua,imax      ,jmax      ,ibMax      , & 
!               jbMax      ,kMax      ,roperm    ,nameg      , &
!               labeli     ,labelf    ,extw      ,exdw       , &
!               TRIM(TRCG) ,TRIM(LV)  ,si        ,sl         , &
!               havesurf   ,kmaxloc)
!
!          CALL write_GridSigma_file (&
!               ifday     ,tod       ,idate     ,idatec     ,&
!               td0       ,tg0       ,tc0       ,z0         ,&
!               convc     ,convt     ,convb     ,gtsea      ,&
!               AlbVisBeam,AlbVisDiff,AlbNirBeam,AlbNirDiff ,&
!               ustar     ,sm0       ,sheleg    ,lsmk       ,&
!               imask     ,mlsi      ,del_in    ,roperm     ,&
!               namer     ,labeli    ,labelf    ,extw       ,&
!               exdw      ,TRIM(TRCG),TRIM(LV)  ,nfprt)
!       END IF
! fim add
        IF (nhdhn .GT. 0) &
             CALL wdhnprog (&
             ndhndrct   ,nfdhn     ,ndhnfcst  ,ifday    ,tod       ,idate     ,&
             idatec     ,qgzs_orig ,lsmk      ,qlnpp    ,qdivp     ,&
             qrotp      ,qqp       ,qtmpp     ,gtsea    ,td0       ,&
             capac0     ,w0        ,imask     ,IVGTYP   ,temp2m    ,umes2m    ,&
             roperm     ,nameh     ,labeli    ,labelf   ,a_hybr    , b_hybr   ,&
             extw       ,exdw      ,TRIM(TRCG),TRIM(LV) )

        IF (mextfm .GT. 0) &
             CALL weprog (&
             nedrct     ,neprog    ,nefcst    ,ifday    ,tod       ,idate     ,&
             idatec     ,a_hybr    , b_hybr   ,qgzs_orig ,lsmk     ,qlnpp     ,qdivp     ,&
             qrotp      ,qqp       ,qtmpp     ,gtsea     ,td0      ,SoilMask  ,&
             capac0     ,w0        ,imask     ,IVGTYP   ,temp2m    ,umes2m    ,&
             roperm     ,namee     ,labeli    ,labelf   ,&
             extw       ,exdw      ,TRIM(TRCG),TRIM(LV) )

        ! Writes the Initialized Fields as an Initial Condition for Further Run

        IF (maxtim <= 0) THEN
           IF(myid.eq.0) THEN
              OPEN (UNIT=nfout1, FILE=TRIM(FNameRestInput2),FORM='unformatted', &
                   ACCESS='sequential',ACTION='write', STATUS='replace', IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                      TRIM(FNameRestInput2), ierr
                 STOP "**(ERROR)**"
              END IF
           END IF
           CALL gwrite(nfout1,ifday ,tod   ,idate ,idatec,qlnpp ,qtmpp , &
                qdivp ,qrotp ,qqp , a_hybr, b_hybr, qgzs_orig)
           IF(myid.eq.0) CLOSE(UNIT=nfout1)
        END IF
     END IF
  END IF
  nlnminit = .FALSE.

  ! Force stop when labeli=labec=labelf
  IF (maxtim <= 0) STOP ' Model Ended Without Time Integration'

  fdh=24.0_r8
  dth=86400.0_r8/fdh
  delth=delt/dth
  fdayh=REAL(ifday,r8)*fdh+tod/dth
  cthw=fdayh

  CALL SetOutPut (tod,idatec)
  CALL SetOutPutDHN (tod,idatec)
!  IF (dctd .GT. 0.0_r8) THEN
!     cdhr(1)=ctdh0
!     DO l=2,mdxtfm
!        cdhr(l)=cdhr(l-1)+dctd
!     END DO
!     cdhr(mdxtfm:maxtid) = -100.0_r8
!  ELSE
!     cdhr = -100.0_r8
!  END IF

  IF (dcte .GT. 0.0_r8) THEN
     !hmjb: This is the same as doing:
     !    cehr(l)=cteh0+dcte*(l-1)
     ! but avoids the multiplications
     cehr(1)=cteh0
     DO l=2,mextfm
        cehr(l)=cehr(l-1)+dcte
     END DO
  END IF

!!$  WRITE(UNIT=nfprt,FMT=*)' '
!!$  WRITE(UNIT=nfprt,FMT=*)' Time Step List Output, cehl:'

!hmjb
! Searching cehr() between 1:mextfm may take too long if 
! the user set NHEXT>50years. This might be the case, for
! instance, when running a long climatic run and asking 
! the model to save PRG's every 6hr until the end of the
! run.
! 
! The alternative code below uses the fact that 
!    cehr(l)=cteh0+dcte*(l-1)
! and 
!    -0.00001_r8 <= cehr(l)-cthw <= 0.00001_r8
! hence
!    (chtw-cteh0-0.00001_r8)/dcte <= l <= (chtw-cteh0+0.00001_r8)/dcte
!
! Therefore, cehl(jdt)=T only when there is an integer
! value between these two real limits.
!

  cthw=fdayh
  DO jdt=1,maxtim
     cthw=cthw+delth
    IF(dcte == 0.0) CYCLE
     l1 = 1+CEILING((cthw-cteh0-0.00001_r8)/dcte)
     l2 = 1+FLOOR  ((cthw-cteh0+0.00001_r8)/dcte)
     IF (l1.eq.l2) THEN
        cehl(jdt)=.TRUE.
     ENDIF
  END DO

!  cthw=fdayh
!  DO jdt=1,maxtim
!     cthw=fdayh + (jdt-1)*delth
!     DO l=1,maxtfm
!        IF (ABS(cdhr(l)-cthw) .LE. 0.00001_r8) THEN
!           cdhl(jdt)=.TRUE.
!        END IF
!     END DO
!  END DO
  !
  maxt0=NINT((REAL(ifday,r8)*fdh+tod/dth)/delth)
  !
  !     this is to remove accumulations from nlnmi
  !
  geshem=0.0_r8
  !
  !     clear all diagnostic accumulators
  !
  CALL rsdiag
  !
  !     check files
  !     if nfin0=nfin1   then  cold start
  !
  limlow=1
  IF(nfin0.EQ.nfin1)THEN
     CALL MsgOne(h," Init: Cold Start")
     !
     !     read cloud dataset for cold start
     !
     IF(schemes==1)THEN
        CALL InitCheckSSiBFile(iMax,jMax,ibMax,&
             jbMax  ,kMax, ifdy  ,ids   ,idc   ,ifday , &
             tod   ,idate ,idatec,todsib  ,ibMaxPerJB)
     ELSE IF (schemes==2) THEN
       CALL InitCheckSiB2File(iMax,jMax,ibMax,&
             jbMax  ,kMax, ifdy ,ids   ,idc   ,ifday , &
             tod   ,idate ,idatec   ,todsib,ibMaxPerJB )
     ELSE IF (schemes==3) THEN
        ! STOP "vazio"
     END IF          
     CALL InitCheckFileConvec(ifdy  ,todcld,ids   ,idc  )
     !
     !     cold start (at first delt/4 ,then delt/2 )
     !
     limlow =2
     dt= delt /4.0_r8
     !
     ! filter arguments for first time step
     !
!!   fgqm = 0.0_r8
!!   fgqmm = 0.0_r8
!!   fgum = 0.0_r8
!!   fgvm = 0.0_r8
!!   fgtmpm = 0.0_r8
!!   fgdivm = 0.0_r8
!!   fgtlamm = 0.0_r8
!!   fgtphim = 0.0_r8
!!   fglnpm = 0.0_r8
!!   fgplamm = 0.0_r8
!!   fgpphim = 0.0_r8
!     ngra = 0
!     OPEN(unit=199,file='graphout',access='sequential',form='unformatted')
!    CALL SetJablo
!    CALL Inittrans(rcl)
!    CALL SfcGeoTrans(slagr)
     fa  = 0.0_r8
     fb  = 1.0_r8
     fb1 = 1.0_r8
     nlnminit = .FALSE.
     init_globconserv = .TRUE.
!    init_globfluxconserv= .TRUE.
     bckhum = .TRUE.
     dotrac=  .not.slhum
     dohum = .not.slhum
     IF (nscalars.gt.0) THEN
        !do_globconserv = .FALSE.
        do_globconserv = .TRUE.
     ENDIF
     IF (nAeros.gt.0) THEN
        !do_globfluxconserv = .FALSE.
        do_globfluxconserv = .TRUE.
     ENDIF

     DO jdt=1,2
        !
        !     snow reinitialization after surface temperature initialization
        !
        IF( initlz >=1 .AND. jdt == 2 ) THEN
           CALL InitSurfTemp(jbMax ,ibMaxPerJB)
        END IF

        IF (jdt == 2) THEN
           CALL rsdiag()
           CALL TurnOnGridHistory()
        END IF

        istrt=jdt
        IF(nfctrl(7).GE.1)WRITE(UNIT=nfprt,FMT=104) jdt
        !
        !     calculate matrices for semi-implicit integration
        !
        delta2 = dt
        IF(slagr.and.SL_twotime_scheme) delta2 = dt/2._r8
        CALL bmcm(delta2)
        !
        ! perform time step
        !
        ifilt=0
        !$OMP PARALLEL
        CALL TimeStep(fb1,fa,fb,dotrac,slagr,slhum,dohum,bckhum,grid_difus,&
	              nlnminit,.FALSE.,.FALSE.,dt,jdt,ifday,tod,idatec,jdt)
        !$OMP END PARALLEL
!        if (jdt.eq.1) then
!        ngra = ngra + 1
!        call dumpgra
!        endif
        IF(jdt.EQ.2) THEN
           !
           !     accumulate time mean prognostic fields if requested
           !
           CALL accpf (ifday, tod, qtmpp, qrotp, qdivp, qqp, qlnpp, nfdyn)
           !
           !     diagnostic of preciptation if requested
           !
           CALL Prec_Diag (ifday, tod, prct, prcc, nfprc)

        END IF


        ! prepare next time step, including filter arguments

        dt=dt*2.0_r8
        IF(IsGridHistoryOn())THEN
           IF(myid.eq.0)THEN
              OPEN (UNIT=nfghou, FILE=TRIM(FNameOutGH),FORM='unformatted',&
                   ACCESS='sequential',ACTION='write',STATUS='replace',IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                      TRIM(FNameOutGH), ierr
                 STOP "**(ERROR)**"
              END IF
           END IF
           CALL WriteGridHistory (nfghou, ifday, tod, idate)
        END IF
        ktm=kt
        fsbc=.FALSE.
        CALL MsgOne(h,"Time integration starts")
        StartStorDiag=.TRUE.
        fb1 = 0.0_r8
        init_globconserv = .FALSE.
        init_globfluxconserv= .FALSE.
        IF (slhum) bckhum = .FALSE.
     END DO
     tod=dt !PYK
  ELSE
     CALL MsgOne(h,"Time integration starts")
     StartStorDiag=.TRUE.
     CALL TurnOnGridHistory()  
    IF(IsGridHistoryOn())THEN
       IF(myid.eq.0)THEN
              OPEN (UNIT=nfghou, FILE=TRIM(FNameOutGH),FORM='unformatted',&
                   ACCESS='sequential',ACTION='write',STATUS='replace',IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                      TRIM(FNameOutGH), ierr
                 STOP "**(ERROR)**"
              END IF
       END IF
       !CALL WriteGridHistory (nfghou, ifday, tod, idatec)
    END IF
  END IF
  !
  !     smooth start
  !
  !IF(TRIM(igwd).EQ.'YES')dogwd=0
  IF(TRIM(igwd).EQ.'YES' .or. TRIM(igwd).EQ.'CAM' .or. TRIM(igwd).EQ.'USS'.or. TRIM(igwd).EQ.'GMB')dogwd=0
  ahour=ifday*24.0_r8+tod/3600.0_r8
  istrt=0
  !
  ! semi-implicit matrices
  !
  delta2 = dt
  IF(slagr.and.SL_twotime_scheme) delta2 = dt/2._r8
  CALL bmcm(delta2)
  !
  ! filter arguments for all remaining time steps
  !
  fa     = filta
  fb     = filtb
  fb1    = 0.0_r8
  dohum = .not.slhum
  bckhum = .not.slhum
  dotrac = .not.slhum
  init_globconserv = .FALSE.
  init_globfluxconserv= .FALSE.
  IF(slagr.and.SL_twotime_scheme) THEN
     fa     = 1.0_r8
     fb     = 0.0_r8
  END IF
  IF( TRIM(start) == "warm" )THEN
     CLOSE(UNIT=nfsibi,STATUS='KEEP')
     CLOSE(UNIT=nfcnv0,STATUS='KEEP')
     CLOSE(UNIT=nfin0,STATUS='KEEP')
     CLOSE(UNIT=nfin1,STATUS='KEEP')
     fb1    = fb
     IF(slagr.and.SL_twotime_scheme)   fb1 = 0.0_r8
  ENDIF

    RETURN 

7   FORMAT( ' SIMPLIFIED PHYSICS OPTION IN EFFECT'/ &
         '  ALL OTHER OPTIONS OVERRIDDEN'/ &
         '  INITIAL CONDITIONS ISOTHERMAL RESTING ATMOSPHERE'/ &
         '  FLAT LOWER BOUNDARY'/ &
         ' PHYSICS:'/ &
         '  NEWTONIAN COOLING'/ &
         '  RALEIGH DAMPING'/ &
         '  NO WATER VAPOR OR CONDENSATION EFFECTS')
104 FORMAT(' ITERATION COUNT FOR THE COLD START=',I2)
336 FORMAT(' FOR IFSST=',I5,' SSTLAG MUST BE SET NONNEGATIVE.  NOT ',G12.5)

  END SUBROUTINE atmos_model_init




  SUBROUTINE atmos_model_run()

  USE FieldsPhysics, ONLY:     &
       capac0,geshem,gtsea, prct, prcc,td0,w0,gndvi, &
       z0,temp2m,umes2m,sfc,statec,&
       AlbVisDiff,o3mix,wsib3d,iMask,SoilMask,sheleg,soilm,tracermix,gco2flx,&
!tar begin
!climate aerosol parameters of coarse mode
       aod,asy,ssa,z_aer, &
!climate aerosol parameters of fine mode 2000
       aodF,asyF,ssaF,z_aerF, &      
!tar end       
! add solange 27-01-2012
       ustar,&
       tg0,tc0,convc,convt,convb,&
       sm0,mlsi,   &
       AlbNirBeam, AlbNirDiff, AlbVisBeam,sfc

! fim add

  USE Sizes, ONLY:             &
       myMNMax               , &
       imaxperj              , &
       ibMax                 , &
       jbMax                 , &
       mnMax                 , &
       mMax                  , &
       nMax                  , &
       mnMap                 , &
       kMax                  , &
       imax                  , &
       jmax                  , &
       a_hybr                , &
       b_hybr                , &
       ijMax                 , &
       ijMaxGauQua           , &
       ibPerIJ               , &
       jbPerIJ               , &
       ibMaxPerJB            , &
! add solange 21-01-2012
       havesurf              , &
       kMaxloc
! fim add



  USE GridHistory, ONLY:      &
       WriteGridHistory      , &
       TurnOnGridHistory     , &
       IsGridHistoryOn       , &
       WriteGridHistoryTopo

  USE InputOutput, ONLY:       &
       UpDateGetsbc          , &       
       gread4                , &
       gwrite                , &
       fsbc

  USE Diagnostics, ONLY:      &
       StartStorDiag         , &
       rsdiag                , &
       Prec_Diag             , &
       accpf                 , &
       wridia                , &
       wdhnprog              , &
       weprog                , &
       wrprog

  USE FieldsDynamics, ONLY: &
       fgq, fgqm, fgqmm, fgum, fgvm, fgtmpm, fgdivm, &
       fgtlamm, fgtphim, fglnpm, fgplamm, fgpphim, &
       qgzs, qlnpp,qlnp_nabla4, qtmpp, qdivp, qrotp, qqp,fgtmp,fgzs, &
       fgpass_scalars, adr_scalars, fgice, fgicem, fgliq, fgliqm,omg, fgvar, fgvarm,&
       fgprsl,fgprsi,fgphil,fgphii
       
  USE SFC_SSiB, ONLY:          &
       ReStartSSiB

  USE SFC_SiB2, ONLY:          &
       ReStartSiB2

  USE Sfc_Ibis_Fiels, ONLY:          &
       ReStartIBIS,&
       vegtype0

  USE Convection, ONLY:       &
       ReStartConvec

  USE Utils, ONLY:     &
       tmstmp2,        &
       colrad,         &
       colrad2D,       &
       InitTimeStamp,  &
       TimeStamp,      &
       total_mass,    &
       total_flux,nTtimes ,aTfluxco2,totflux

  USE ModTimeStep, ONLY:      &
       SfcGeoTrans,           &
       TimeStep  
       
  USE SemiLagrangian, ONLY: &
       InitSL , &
       ulonm  , &
       ulatm  , &
       uetam  , &
       ulonm2D, &
       ulatm2D


  USE  Options, ONLY: &
       CreateFileName,&
       FNameRestOutput2,FNameRestInput2,&
       FNameSibPrgOut1,FNameConvClInp0,FNameConvClOut1,FNameRestInput1,FNameRestOutput1,FNameSibPrgInp0,&
       idate, idatec, idatef,maxtim,&
       ifilt,&
       jdt,kt,ktp,ktm,dt,&
       mextfm,nfctrl,&
       isimp,&
       nlnminit,SL_twotime_scheme,slhum,grid_difus,reststep,slagr, initlz,&
       rmRestFiles,GenRestFiles,&
       nfdrct,nfout0,nfin1,nfprt,nfsibo,nfout1,nfdhn,nhdhn,& 
       nefcst,ndhnfcst,nfdiag, nfcnv1,nedrct,ndhndrct,ndhndrct,neprog,nffcst,nftmp, &
       nfghou,nfdyn,nfprc ,&
       EXDW ,EXTW ,UNIFIED,&
       cthl,cdhl,&
       LV,TRCG,schemes,nscalars,microphys,nClass,nAeros,GenAssFiles

! add solange 27-01-2012
  USE SpecDump , ONLY: write_sigma_file
  USE GridDump , ONLY: write_GridSigma_file
!fim add

    IMPLICIT NONE

    INCLUDE 'mpif.h'
    INTEGER              :: jhr
    INTEGER              :: jmon
    INTEGER              :: jday
    INTEGER              :: jyr
    REAL(KIND=r8)       :: ahour
    INTEGER              :: ierr
    INTEGER              :: maxstp
    INTEGER              :: j
    INTEGER              :: i
    INTEGER              :: nPtland
    LOGICAL              :: ExistFile


    !$ INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
    CHARACTER(LEN=*), PARAMETER :: h="**(atmos_model_run)**" 
    !
    ! time step loop
    !
  ifilt=1
  IF (instrument) THEN
     CALL ChangeWatch(3)
  END IF
  DO jdt=limlow,maxtim
     !
     !     step loop starts
     !
     IF(TRIM(isimp).NE.'YES') THEN
        CALL UpDateGetsbc(ifday, tod, idate, idatec,ibMax,jbMax,kMax,ibMaxPerJB,&
                        AlbVisDiff,gtsea,gco2flx,gndvi,soilm,wsib3d,sheleg,o3mix,tracermix, &
!tar begin
!climate aerosol parameters
               aod,asy,ssa,z_aer,aodF,asyF,ssaF,z_aerF)
!tar end
     END IF
     tod=tod+dt
     IF(ABS( MOD(tod+0.03125_r8,86400.0_r8)-0.03125_r8).LT.0.0625_r8)THEN
        tod=0.0_r8
        ifday=ifday+1
     END IF
     CALL tmstmp2(idate,ifday,tod,jhr,jday,jmon,jyr)
     idatec(1)=jhr
     idatec(2)=jmon
     idatec(3)=jday
     idatec(4)=jyr
     ahour=(ifday*24.0e0_r8)+(tod/3.6e3_r8)

     kt   =INT(ahour-(1.0e-2_r8))
     ktp  =INT(ahour+(dt/3.6e3_r8)-(1.0e-2_r8))
     IF(jdt.EQ.maxtim) THEN
        ktm=kt
     END IF
     IF (slhum) dohum = cthl(jdt)
     !
     ! perform time step
     !     
     !$OMP PARALLEL
     CALL TimeStep(fb1,fa,fb,dotrac,slagr,slhum,dohum,bckhum,grid_difus, &
                   nlnminit,.FALSE.,enhdifl,dt,jdt,ifday,tod,idatec,jdt)
     !$OMP END PARALLEL

     !
     !     accumulate time mean prognostic fields if requested
     !
     CALL accpf (ifday, tod, qtmpp, qrotp, qdivp, qqp, qlnpp, nfdyn)

     !
     !     diagnostic of preciptation if requested
     !
     CALL Prec_Diag (ifday, tod, prct, prcc, nfprc)

     IF(IsGridHistoryOn())THEN
        CALL WriteGridHistory (nfghou, ifday, tod, idate)
     END IF
     fsbc=.FALSE.
     ktm=kt

     !
     ! output, if desired
     !
     !*JPB IF(MOD(jdt,reststep)==0.OR.jdt.EQ.maxtim) THEN
!     IF (reststep == 0) THEN
!        lreststep=(jdt == maxtim .AND. GenRestFiles)
!     ELSE IF (reststep > 0) THEN
!        lreststep=((MOD(jdt,reststep) == 0 .OR. jdt == maxtim) .AND. GenRestFiles)
!     ELSE IF (reststep < 0 .and. reststep /= -365) THEN
!        lreststep=((idatec(3) == 1 .AND. idatec(1) ==  idate(1) .AND. &
!                  tod == 0.0_r8 .OR. jdt == maxtim) .AND. GenRestFiles)
!     ELSE IF (reststep == -365) THEN
!        lreststep=((idatec(3) == 1 .AND. idatec(1) ==  idate(1) .AND. idatec(2) == 1  .AND.&
!                  tod == 0.0_r8 .OR. jdt == maxtim) .AND. GenRestFiles)
!     END IF

     IF (reststep == 0) THEN
        lreststep=((jdt == maxtim .and. .not. GenAssFiles) .AND. GenRestFiles)
     ELSE IF (reststep > 0) THEN
        lreststep=((MOD(jdt,reststep) == 0 .OR. (jdt == maxtim .and. .not. GenAssFiles)) .AND. GenRestFiles)
     ELSE IF (reststep < 0 .and. reststep /= -365) THEN
        lreststep=((idatec(3) == 1 .AND. idatec(1) ==  idate(1) .AND. &
                  tod == 0.0_r8 .OR. (jdt == maxtim .and. .not. GenAssFiles)) .AND. GenRestFiles)
     ELSE IF (reststep == -365) THEN
        lreststep=((idatec(3) == 1 .AND. idatec(1) ==  idate(1) .AND. idatec(2) == 1  .AND.&
                  tod == 0.0_r8 .OR. (jdt == maxtim .and. .not. GenAssFiles)) .AND. GenRestFiles)
     END IF


     IF (lreststep) THEN
        !
        !     write history wave-data
        !
        WRITE(labeli,'(I4.4,3I2.2)')idate (4),idate (2),idate (3),idate (1)
        WRITE(labelc,'(I4.4,3I2.2)')idatec(4),idatec(2),idatec(3),idatec(1)

        IF (rmRestFiles) THEN

           INQUIRE(File=TRIM(FNameSibPrgInp0), Exist=ExistFile)
           IF(ExistFile)THEN
              OPEN(UNIT=nfsibo,FILE=TRIM(FNameSibPrgInp0),FORM='unformatted',ACCESS='sequential',&
                ACTION='write',STATUS='unknown',IOSTAT=ierr)
           ELSE 
              OPEN(UNIT=nfsibo,FILE=TRIM(FNameSibPrgInp0),FORM='unformatted',ACCESS='sequential',&
                ACTION='write',STATUS='unknown',IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                      TRIM(FNameSibPrgInp0), ierr
                 STOP "**(ERROR)**"
              END IF
           END IF 
           INQUIRE(File=TRIM(FNameConvClInp0), Exist=ExistFile)
           IF(ExistFile)THEN
              OPEN(UNIT=nfcnv1, FILE=TRIM(FNameConvClInp0 ),FORM='unformatted',ACCESS='sequential',&
                  ACTION='write',STATUS='unknown',IOSTAT=ierr)
           ELSE 
              OPEN(UNIT=nfcnv1, FILE=TRIM(FNameConvClInp0 ),FORM='unformatted',ACCESS='sequential',&
                  ACTION='write',STATUS='unknown',IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                   TRIM(FNameConvClInp0 ), ierr
                 STOP "**(ERROR)**"
              END IF
           END IF
           INQUIRE(File=TRIM(FNameRestInput2), Exist=ExistFile)
           IF(ExistFile)THEN
              OPEN(UNIT=nfout1, FILE=TRIM(FNameRestInput2),FORM='unformatted',ACCESS='sequential',&
                   ACTION='write',STATUS='unknown',IOSTAT=ierr)
           ELSE
              OPEN(UNIT=nfout1, FILE=TRIM(FNameRestInput2),FORM='unformatted',ACCESS='sequential',&
                   ACTION='write',STATUS='unknown',IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                      TRIM(FNameRestInput2), ierr
                 STOP "**(ERROR)**"
              END IF
           END IF
           INQUIRE(File=TRIM(FNameRestInput1), Exist=ExistFile)
           IF(ExistFile)THEN
              OPEN(UNIT=nfout0, FILE=TRIM(FNameRestInput1),FORM='unformatted',ACCESS='sequential',&
                ACTION='write', STATUS='unknown',IOSTAT=ierr)
           ELSE
              OPEN(UNIT=nfout0, FILE=TRIM(FNameRestInput1),FORM='unformatted',ACCESS='sequential',&
                ACTION='write', STATUS='unknown',IOSTAT=ierr)
              IF (ierr /= 0) THEN
                 WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                      TRIM(FNameRestInput1), ierr
                STOP "**(ERROR)**"
              END IF
           END IF
           CLOSE(UNIT=nfout0,STATUS='DELETE')
           CLOSE(UNIT=nfsibo,STATUS='DELETE')
           CLOSE(UNIT=nfcnv1,STATUS='DELETE')
           CLOSE(UNIT=nfout1,STATUS='DELETE')

        END IF

        IF(initlz < 0)THEN 
           CALL CreateFileName(TRIM(PRC),2,'RESTAT_SIB')
        ELSE
           CALL CreateFileName(TRIM(PRC),2)
        END IF

        OPEN(UNIT=nfsibo, FILE=TRIM(FNameSibPrgOut1), FORM='unformatted',ACCESS='sequential',&
             ACTION='write',STATUS='replace',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                TRIM(FNameSibPrgOut1), ierr
           STOP "**(ERROR)**"
        END IF
        OPEN(UNIT=nfcnv1, FILE=TRIM(FNameConvClOut1), FORM='unformatted',ACCESS='sequential',&
             ACTION='write',STATUS='replace',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                TRIM(FNameConvClOut1), ierr
           STOP "**(ERROR)**"
        END IF
        OPEN(UNIT=nfout1, FILE=TRIM(FNameRestOutput2), FORM='unformatted',ACCESS='sequential',&
             ACTION='write',STATUS='replace',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                TRIM(FNameRestOutput2), ierr
           STOP "**(ERROR)**"
        END IF
        OPEN(UNIT=nfout0, FILE=TRIM(FNameRestOutput1), FORM='unformatted',ACCESS='sequential',&
             ACTION='write',STATUS='replace',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
                TRIM(FNameRestOutput1), ierr
           STOP "**(ERROR)**"
        END IF

        WRITE(UNIT=nfout0) fgqm, fgum, fgvm, fgtmpm, fgdivm, fglnpm, fgtlamm, &
             fgtphim, fgplamm, fgpphim
        IF (slhum) WRITE(UNIT=nfout0) fgq
        IF(UNIFIED)WRITE(UNIT=nfout0) fgprsl,fgprsi,fgphil,fgphii
        IF (microphys) THEN
           IF((nClass+nAeros) >0 )THEN
              WRITE(UNIT=nfout0) fgicem, fgice, fgliqm, fgliq,omg, fgvarm,fgvar
           ELSE
              WRITE(UNIT=nfout0) fgicem, fgice, fgliqm, fgliq,omg
           END IF
        END IF 
        IF (slagr.or.slhum) THEN
           WRITE(UNIT=nfout0)ulonm ,ulatm ,uetam,ulonm2D,ulatm2D
        END IF
        IF (nscalars>0) THEN
           WRITE(UNIT=nfout0) fgpass_scalars, adr_scalars
        END IF

        WRITE(UNIT=nfout1) ifday,tod,idate,idatec,a_hybr,b_hybr
        WRITE(UNIT=nfout1) qgzs_orig,qlnpp,qlnp_nabla4,qtmpp,qdivp,qrotp, &
	                    qqp,total_mass,total_flux,nTtimes ,aTfluxco2,totflux


        CALL ReStartConvec (ifday,tod ,idate ,idatec, &
             nfcnv1)
        IF(schemes==1) THEN
           CALL ReStartSSiB   (jbMax,ifday,tod ,idate ,idatec, &
                nfsibo,ibMaxPerJB)
        ELSE IF(schemes==2) THEN
           CALL ReStartSiB2   (jbMax,ifday,tod ,idate ,idatec, &
                nfsibo,ibMaxPerJB)                
        ELSE IF(schemes==3) THEN
           CALL ReStartIBIS(nfsibo)
           !STOP 'vazio'
        END IF
        FNameRestInput2 = FNameRestOutput2
        FNameRestInput1 = FNameRestOutput1
        FNameConvClInp0 = FNameConvClOut1
        FNameSibPrgInp0 = FNameSibPrgOut1
        CLOSE(UNIT=nfsibo)
        CLOSE(UNIT=nfcnv1)
        CLOSE(UNIT=nfout1)
        CLOSE(UNIT=nfout0)

     ENDIF

     IF(cthl(jdt)) THEN
        maxstp=NINT((REAL(ifday,r8)*fdh+tod/dth)/delth)-maxt0

!!$        WRITE(UNIT=nfprt,FMT="('Write file at timestep ',i5)") jdt
        !
        !     reset precip. every maxstp time steps
        !
        IF(schemes==3)THEN
           IVGTYP=vegtype0
        ELSE
           DO j=1,jbMax
              nPtland=0
              DO i=1,ibMaxPerJB(j)
                 IF(imask(i,j) >=1_i8)THEN
                    nPtland=nPtland+1
                    IVGTYP(nPtland,j)=imask(i,j)
                 END IF  
              END DO
           END DO   
        END IF
!      ngra = ngra + 1
!      CALL dumpgra
        CALL wrprog (&
             nfdrct   ,nfdiag   ,ifday     ,tod       ,idate  ,idatec , &
             qrotp    ,qdivp    ,qqp       ,qlnpp     ,qtmpp  ,gtsea  , &
             td0      ,SoilMask ,capac0    ,w0        ,imask  ,IVGTYP , &
             temp2m   ,umes2m   ,nffcst    ,nftmp     ,a_hybr ,b_hybr , qgzs_orig, &
             lsmk     ,tg0      ,sheleg    ,mlsi      , &
             roperm   ,namef    ,labeli    ,labelf    ,extw   ,exdw   , &
             TRIM(TRCG),TRIM(LV)  ,.TRUE. &
             )
! add solange 27-01-2012
!       IF(GenAssFiles)THEN
!          CALL write_sigma_file (&
!               ifday      ,tod        ,idate     ,idatec     , &
!               qrotp      ,qdivp      ,qqp       ,qlnpp      , &
!               qtmpp      ,qgzs_orig  ,o3mix     ,del_in     , &
!               ijMaxGauQua,imax       ,jmax      ,ibMax      , &
!               jbMax      ,kMax       ,roperm    ,nameg      , &
!               labeli     ,labelf     ,extw      ,exdw       , &
!               TRIM(TRCG) ,TRIM(LV)   ,si        ,sl         , &
!               havesurf   ,kmaxloc)

!          CALL write_GridSigma_file (&
!               ifday     ,tod        ,idate     ,idatec     ,&
!               td0       ,tg0        ,tc0       ,z0         ,&
!               convc     ,convt      ,convb     ,gtsea      ,&
!               AlbVisBeam,AlbVisDiff ,AlbNirBeam,AlbNirDiff ,&
!               ustar     ,sm0        ,sheleg    ,lsmk       ,&
!               imask     ,mlsi       ,del_in    ,roperm     ,&
!               namer     ,labeli     ,labelf    ,extw       ,&
!               exdw      ,TRIM(TRCG) ,TRIM(LV)  ,nfprt)
!       END IF
! fim add

        CALL wridia(nfdiag, maxstp, idatec)
        !
        !     zero reset diagnostic fields
        !
        CALL rsdiag
        geshem=0.0_r8
        limlow=1

        IF(jdt.NE.maxtim)THEN
           maxt0=NINT((REAL(ifday,r8)*fdh+tod/dth)/delth)
        END IF
     END IF

     IF (cehl(jdt)) THEN
        maxstp=NINT((REAL(ifday,r8)*fdh+tod/dth)/delth)-maxt0

        IF(nfctrl(8).GE.1) THEN
!!$           WRITE(UNIT=nfprt,FMT=102)dt,ifday,tod, &
!!$                REAL(ifday,r8)*fdh+tod/dth,jdt,maxstp
        END IF
        IF(schemes==3)THEN
           IVGTYP=vegtype0
        ELSE
           DO j=1,jbMax
              nPtland=0
              DO i=1,ibMaxPerJB(j)
                 IF(imask(i,j) >=1_i8)THEN
                    nPtland=nPtland+1
                    IVGTYP(nPtland,j)=imask(i,j)
                 END IF  
              END DO
           END DO   
        END IF
        
        IF (mextfm .GT. 0) &
             CALL weprog (&
             nedrct    ,neprog    ,nefcst    ,ifday     ,tod       ,idate    , &
             idatec    ,a_hybr    ,b_hybr    ,qgzs_orig ,lsmk      ,qlnpp     ,qdivp    , &
             qrotp     ,qqp       ,qtmpp     ,gtsea     ,td0       ,SoilMask , &
             capac0    ,w0        ,imask     ,IVGTYP    ,temp2m    ,umes2m   , &
             roperm    ,namee     ,labeli    ,labelf    , &
             extw      ,exdw      ,TRIM(TRCG),TRIM(LV) )
     END IF

     IF (cdhl(jdt)) THEN
        maxstp=NINT((REAL(ifday,r8)*fdh+tod/dth)/delth)-maxt0

        IF(nfctrl(8).GE.1) THEN
!!$           WRITE(UNIT=nfprt,FMT=102)dt,ifday,tod, &
!!$                REAL(ifday,r8)*fdh+tod/dth,jdt,maxstp
        END IF
        IF(schemes==3)THEN
           IVGTYP=vegtype0
        ELSE
           DO j=1,jbMax
              nPtland=0
              DO i=1,ibMaxPerJB(j)
                 IF(imask(i,j) >=1_i8)THEN
                    nPtland=nPtland+1
                    IVGTYP(nPtland,j)=imask(i,j)
                 END IF  
              END DO
           END DO   
        END IF
        IF (nhdhn .GT. 0) &
             CALL wdhnprog (&
             ndhndrct  ,nfdhn     ,ndhnfcst  ,ifday     ,tod       ,idate    , &
             idatec    ,qgzs_orig ,lsmk      ,qlnpp     ,qdivp    , &
             qrotp     ,qqp       ,qtmpp     ,gtsea     ,td0      , &
             capac0    ,w0        ,imask     ,IVGTYP    ,temp2m    ,umes2m   , &
             roperm    ,nameh     ,labeli    ,labelf    ,a_hybr    ,b_hybr   , &
             extw      ,exdw      ,TRIM(TRCG),TRIM(LV) )
        
     END IF

     fb1 = fb
     IF(slagr.and.SL_twotime_scheme)   fb1 = 0.0_r8
  ENDDO
!  if (myid.eq.0) CALL writectl(ngra)
 
  RETURN
  END SUBROUTINE atmos_model_run


  SUBROUTINE atmos_model_finalize(myid_in,maxnodes_in)
    USE  Options, ONLY: nfout0,schemes


    USE SFC_SSiB, ONLY:          &
       Finalize_SSiB

    USE SFC_SiB2, ONLY:          &
       Finalize_SiB2

  USE Sfc_Ibis_Fiels, ONLY:          &
       Finalize_IBIS

   USE Watches, ONLY:  &
       DumpWatches,   &
       DestroyWatches

  USE Parallelism, ONLY:   &
       DestroyParallelism, &
       MsgOne
      
    INCLUDE 'mpif.h'
    INTEGER,    INTENT(IN), OPTIONAL :: myid_in
    INTEGER,    INTENT(IN), OPTIONAL :: maxnodes_in
    CHARACTER(LEN=*), PARAMETER :: h="**(atmos_model_finalize)**" 

    REWIND nfout0
   DEALLOCATE(cehl  )
   DEALLOCATE(cehr  )  
   !DEALLOCATE(cdhr  )
   DEALLOCATE(ct_in )
   DEALLOCATE(cq_in )
   DEALLOCATE(lsmk)
   DEALLOCATE(rlsm)
   DEALLOCATE(qgzs_orig)





   IF(schemes ==1)THEN
      CALL Finalize_SSiB()
   ELSE IF(schemes ==2) THEN
      CALL Finalize_SiB2()
   ELSE IF(schemes ==3) THEN
     CALL Finalize_IBIS()
   END IF    
    !
    !  finish MPI
    !
    IF (instrument) THEN
       CALL DumpWatches(unitDump)
       CALL DestroyWatches()
    END IF
    IF (PRESENT(myid_in) .AND. PRESENT(maxnodes_in) ) THEN
       CALL MsgOne(h," atmos_model_finalize:MPI_finalize is executed at other place")
    ELSE
       CALL DestroyParallelism("***MODEL EXECUTION ENDS NORMALY***")
    END IF

    RETURN

  END SUBROUTINE atmos_model_finalize

END MODULE AtmosModelMod
