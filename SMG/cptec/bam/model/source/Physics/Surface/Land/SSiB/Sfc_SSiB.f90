!
!  $Author: pkubota $
!  $Date: 2008/04/09 12:42:57 $
!  $Revision: 1.12 $
!
! Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
!  Fine aerosol mode 2000 included
!  Modifications (5) are marked by 
!  !tar begin   and  !tar end
!
MODULE SFC_SSiB

  ! InitSSiB
  !
  !fysiks -----| pbl ------| root
  !            |           |
  !            |           | raduse
  !            |           |
  !            |           | stomat
  !            |           |
  !            |           | interc
  !            |           |
  !            |           | sflxes ------| vntlax
  !            |                          |
  !            |                          | rbrd
  !            |                          |
  !            |                          | cut
  !            |                          |
  !            |                          | stres2
  !            |                          |
  !            |                          | temres
  !            |                          |
  !            |                          | update
  !            |                          |
  !            |                          | airmod
  !            |
  !            | snowm
  !            |
  !            | runoff
  !            |
  !            | seasfc ------| vntlt1
  !            |
  !            | sextrp
  !            |
  !            | sibwet ------| extrak
  !            |
  !            | sibwet_GLSM ------| extrak
  !            | 
  !            | Albedo --- radalb
  !            |
  !            | radalb
  !            |
  !            | vegin
  !            |
  !            |re_assign_sib_soil_prop
  !            |
  !            | wheat

  USE Constants, ONLY :     &
       ityp, imon, icg, iwv, idp, ibd, tice,&
       gasr,          &
       pie,           &
       cp,            &
       hl,            &
       grav,          &
       stefan,        &
       snomel,        &
       tf,            &
       epsfac,        &
       clai,          &
       athird,        &
       cw,            &
       z0ice,         &
       oceald   ,     &
       icealn   ,     &
       icealv   ,     &
       r8,i8,r4,i4

  USE Options, ONLY: &
       nfsibt,fNameSibmsk,&
       nftgz0,fNameTg3zrl,&
       nfzol,fNameRouLen,&
       reducedGrid,initlz,&
       ifalb  , ifsst ,ifco2flx, &
       ifslm  , ifsnw , &
       ifozone, sstlag, &
!tar begin
!climate aerosol selection parameter
       ifaeros, &
!tar end
       intsst , fint  , &
       iglsm_w, mxiter, &
       record_type, fNameSoilType ,&
       fNameVegType,fNameSoilMoist,&
       nfsoiltp,nfvegtp,nfslmtp,nfsibi,isimp,&
       nfprt  , nfctrl, nfsibd, nfalb,filta,epsflt,istrt,Model1D,yrl,monl,schemes,&
       ifndvi,intsoilm,&
       ifslmSib2,&
       intndvi,iftracer,OCFLUX,omlmodel,oml_hml0,SLABOCEAN,atmpbl,ICEMODEL
 
  USE Parallelism, ONLY: &
       MsgOne, FatalError

  USE Utils, ONLY: &
       IJtoIBJB, &
       NearestIJtoIBJB, &
       LinearIJtoIBJB,  &
       SplineIJtoIBJB, &
       AveBoxIJtoIBJB, &
       FreqBoxIJtoIBJB, &
       vfirec

  USE InputOutput, ONLY: &
       WillGetSbc, &
       getsbc

  USE IOLowLevel, ONLY: &
       ReadVar      , &
       ReadGetNFTGZ

  USE Sfc_SeaFlux_Interface   , Only :  seasfc

  USE SlabOceanModel  , Only : GetOceanAlb

  USE Sfc_SeaIceFlux_WRF_Model  , Only : GetIceOceanAlb,&
      TC_SeaIce     ,&
      TGS_SeaIce ,&
      TD_SeaIce  ,&
      TA_SeaIce  ,&
      SNOA_SeaIce,&
      SNOB_SeaIce

  USE FieldsPhysics, ONLY: &
      npatches     , &
      npatches_actual, &
      nzg            , &
      sheleg         , &
      imask          , &
      SoilMask       , &
      AlbVisDiff     , &
      gtsea          , &
      gco2flx        , &
      gndvi          , &
      soilm          , &
      o3mix          , &
!tar begin  
!climate aerosol optical parameters of coarse mode   
      aod            , &
      asy            , &
      ssa            , &
      z_aer          , &    
!tar end
!
!tar begin  
!climate aerosol optical parameters of coarse mode   
      aodF            , &
      asyF            , &
      ssaF            , &
      z_aerF          , &    
!tar end 
      tg1            , &
      tg2            , &
      tg3            , &
      rVisDiff       , &
      ssib           , &
      wsib3d         , &
      ppli           , &
      ppci           , &
      capac0         , &
      gl0            , &
      Mmlen          , &
      tseam          , &
      w0             , &
      td0            , &
      tg0            , &
      tc0            , &
      tm0            , &
      qm0            , &
      tmm            , &
      qmm            , &
      tcm            , &
      tgm            , &
      tdm            , &
      wm             , &
      capacm         , &
      qsfc0          , &
      tsfc0          , &
      qsfcm          , &
      tsfcm          , &
      z0             , &
      zorl           , &
      tkemyj         , &
      MskAnt         , &
      tracermix      , &
      HML            , &
      HUML           , &
      HVML           , &
      TSK            , &
      z0sea          , &
      mlsi           , &  ! add solange 13-11-2012
      sm0            , &  ! add solange 13-11-2012
      laymld,       hbath,     tdeep,sdeep,&
      sfc,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 

  IMPLICIT NONE
SAVE


  PRIVATE
  PUBLIC :: Init_SSiB
  PUBLIC :: Finalize_SSiB  
  PUBLIC :: InitCheckSSiBFile
  PUBLIC :: InitSurfTemp
  PUBLIC :: ReStartSSiB
  PUBLIC :: Albedo
  PUBLIC :: Phenology
  PUBLIC :: SSiB_Driver

  PUBLIC :: x0x
  PUBLIC :: xd
  PUBLIC :: xdc
  PUBLIC :: xbc
  REAL(KIND=r8)   :: expcut

  REAL(KIND=r8)   , ALLOCATABLE :: cedfu (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cedir (:,:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cedfu1(:,:,:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cedir1(:,:,:,:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cedfu2(:,:,:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cedir2(:,:,:,:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cledfu(:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cledir(:,:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xmiu  (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: cether(:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xmiw  (:,:)

  REAL(KIND=r8)   , ALLOCATABLE :: ystpar(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: yopt  (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yll   (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yu    (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yefac (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yh1   (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yh2   (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yootd (:)
  REAL(KIND=r8)   , ALLOCATABLE :: yreen (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ycover(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ylt   (:,:)


  REAL(KIND=r8)   , ALLOCATABLE :: rstpar_fixed    (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: chil_fixed      (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: topt_fixed      (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: tll_fixed          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: tu_fixed          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: defac_fixed     (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ph1_fixed          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ph2_fixed          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: rootd     (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: bee          (:)
  REAL(KIND=r8)   , ALLOCATABLE :: phsat     (:)
  REAL(KIND=r8)   , ALLOCATABLE :: satco     (:)
  REAL(KIND=r8)   , ALLOCATABLE :: poros     (:)
  REAL(KIND=r8)   , ALLOCATABLE :: zdepth    (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: green_fixed     (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xcover_fixed    (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: zlt_fixed          (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: x0x          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xd          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: z2          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: z1          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xdc          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: xbc          (:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: zlt(:, :, :)
  REAL(KIND=r8)   , ALLOCATABLE :: xcover  (:, :, :)
  REAL(KIND=r8)   , ALLOCATABLE :: ph2(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ph1(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: green(:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: defac(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: tu(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: tll(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: topt(:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: rstpar(:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: chil(:,:)

  REAL(KIND=r8), TARGET   , ALLOCATABLE :: vcover_gbl (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: zlt_gbl    (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: green_gbl  (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: chil_gbl   (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: topt_gbl   (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: tll_gbl    (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: tu_gbl     (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: defac_gbl  (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ph2_gbl    (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: ph1_gbl    (:,:,:)
  REAL(KIND=r8)   , ALLOCATABLE :: rstpar_gbl (:,:,:,:)
!  REAL(KIND=r8)   , ALLOCATABLE :: snow_gbl   (:,:,:)


  REAL(KIND=r8), ALLOCATABLE :: satcap3d  (:,:,:)    
  REAL(KIND=r8), ALLOCATABLE :: extk_gbl  (:,:,:,:,:)
  REAL(KIND=r8), ALLOCATABLE :: radfac_gbl(:,:,:,:,:)
  REAL(KIND=r8), ALLOCATABLE :: closs_gbl (:,:)      
  REAL(KIND=r8), ALLOCATABLE :: gloss_gbl (:,:)      
  REAL(KIND=r8), ALLOCATABLE :: thermk_gbl(:,:)      
  REAL(KIND=r8), ALLOCATABLE :: p1f_gbl   (:,:)      
  REAL(KIND=r8), ALLOCATABLE :: p2f_gbl   (:,:)      
  REAL(KIND=r8), ALLOCATABLE :: tgeff_gbl (:,:)      
  REAL(KIND=r8), PUBLIC, ALLOCATABLE :: zlwup_SSiB (:,:)      
  REAL(KIND=r8), ALLOCATABLE :: AlbGblSSiB(:,:,:,:)  
  INTEGER(KIND=i8), ALLOCATABLE :: MskAntSSiB (:,:)
  INTEGER(KIND=i8), ALLOCATABLE :: iMaskSSiB  (:,:)
  REAL(KIND=r8), ALLOCATABLE :: glsm_w  (:,:,:) ! initial soil wetness data at soil model
  REAL(KIND=r8), ALLOCATABLE:: veg_type (:,:,:)! SIB veg type
  REAL(KIND=r8), PUBLIC, ALLOCATABLE :: frac_occ(:,:,:) ! fractional area
  ! coverage
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:  ) :: soil_type! FAO/USDA soil texture

  CHARACTER(LEN=200) :: path_in
  CHARACTER(LEN=200) :: fNameSibVeg
  CHARACTER(LEN=200) :: fNameSibAlb

  real, parameter, public :: MAPL_AIRMW  = 28.97                  ! kg/Kmole
  real, parameter, public :: MAPL_H2OMW  = 18.01                  ! kg/Kmole

  real, parameter, public :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0   ! --

CONTAINS


  SUBROUTINE Init_SSiB(ibMax         ,jbMax         ,iMax          ,jMax          , &
                       kMax          ,path          ,fNameSibVeg_in, &
                       fNameSibAlb_in,ifdy          ,ids           , &
                       idc           ,ifday         ,tod           ,todsib        , &
                       idate         ,idatec        , &
                       ibMaxPerJB  )
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: jbMax
    INTEGER, INTENT(IN) :: iMax
    INTEGER, INTENT(IN) :: jMax
    INTEGER, INTENT(IN) :: kMax
    CHARACTER(LEN=*), INTENT(IN   ) :: path
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibVeg_in
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibAlb_in    
    INTEGER         , INTENT(OUT  ) :: ifdy
    INTEGER         , INTENT(OUT  ) :: ids(4)
    INTEGER         , INTENT(OUT  ) :: idc(4)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    REAL(KIND=r8)   , INTENT(OUT  ) :: todsib
    INTEGER         , INTENT(IN   ) :: idate(4)
    INTEGER         , INTENT(IN   ) :: idatec(4)
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(jbMax)

    expcut=- LOG(1.0e53_r8)
    path_in=path
    fNameSibVeg = fNameSibVeg_in
    fNameSibAlb = fNameSibAlb_in
    ALLOCATE(vcover_gbl (ibMax,jbMax,icg) )     
    vcover_gbl=0.0_r8
    sfc%vcover=>vcover_gbl(1:ibMax,1:jbMax,1) 
    ALLOCATE(zlt_gbl    (ibMax,jbMax,icg))      
    zlt_gbl   =0.0_r8
    ALLOCATE(green_gbl  (ibMax,jbMax,icg))      
    green_gbl =0.0_r8
    ALLOCATE(chil_gbl   (ibMax,jbMax,icg))      
    chil_gbl  =0.0_r8
    ALLOCATE(topt_gbl   (ibMax,jbMax,icg))      
    topt_gbl  =0.0_r8
    ALLOCATE(tll_gbl    (ibMax,jbMax,icg))      
    tll_gbl =0.0_r8  
    ALLOCATE(tu_gbl     (ibMax,jbMax,icg))      
    tu_gbl=0.0_r8 
    ALLOCATE(defac_gbl  (ibMax,jbMax,icg))      
    defac_gbl =0.0_r8
    ALLOCATE(ph2_gbl    (ibMax,jbMax,icg))      
    ph2_gbl   =0.0_r8
    ALLOCATE(ph1_gbl    (ibMax,jbMax,icg))      
    ph1_gbl  =0.0_r8 
    ALLOCATE(rstpar_gbl (ibMax,jbMax,icg,iwv)) 
    rstpar_gbl=0.0_r8
!    ALLOCATE(snow_gbl   (ibMax,jbMax,icg)) 
!    snow_gbl =0.0_r8 
    ALLOCATE( satcap3d  (ibMax,icg,jbMax))
    satcap3d=0.0_r8
    ALLOCATE( extk_gbl  (ibMax,icg,iwv,ibd,jbMax))
    extk_gbl=0.0_r8
    ALLOCATE( radfac_gbl(ibMax,icg,iwv,ibd,jbMax))
    radfac_gbl=0.0_r8
    ALLOCATE( closs_gbl (ibMax,jbMax))
    closs_gbl=0.0_r8
    ALLOCATE( gloss_gbl (ibMax,jbMax))
    gloss_gbl=0.0_r8
    ALLOCATE( thermk_gbl(ibMax,jbMax))
    thermk_gbl=0.0_r8
    ALLOCATE( p1f_gbl   (ibMax,jbMax))
    p1f_gbl=0.0_r8
    ALLOCATE( p2f_gbl   (ibMax,jbMax))
    p2f_gbl=0.0_r8
    ALLOCATE( tgeff_gbl (ibMax,jbMax))
    tgeff_gbl=0.0_r8
    ALLOCATE(  zlwup_SSiB(ibMax,jbMax))
    zlwup_SSiB=0.0_r8
    ALLOCATE( AlbGblSSiB(ibMax,2,2,jbMax))
    AlbGblSSiB=0.0_r8
    ALLOCATE(MskAntSSiB(ibMax,jbMax))   
    MskAntSSiB=0_i8
    ALLOCATE(iMaskSSiB(ibMax,jbMax)) ;iMaskSSiB=0_i8
    ALLOCATE(glsm_w   (ibMax,jbMax,nzg     ))
    glsm_w=0.0_r8
    ALLOCATE(veg_type (ibMax,jbMax,npatches))
    veg_type=0.0_r8
    ALLOCATE(frac_occ (ibMax,jbMax,npatches))
    frac_occ=0.0_r8
    ALLOCATE(soil_type(ibMax,jbMax         ))
    soil_type=0.0_r8
    IF(TRIM(isimp).NE.'YES') THEN
    
       CALL InitSSiBBoundCond(iMax      ,jMax      ,ibMax,jbMax,kMax      ,ifday     ,&
                              tod       ,idate     ,&
                              idatec    ,ibMaxPerJB,ifdy      ,&
                              ids       ,idc       ,todsib    )
    END IF
   
  END SUBROUTINE Init_SSiB

  SUBROUTINE InitSSiBBoundCond(iMax      ,&
                               jMax      ,&
                               ibMax     ,& 
                               jbMax     ,&
                               kMax      ,&
                               ifday     ,&
                               tod       ,&
                               idate     ,&
                               idatec    ,&
                               ibMaxPerJB,&
                               ifdy      ,&
                               ids       ,&
                               idc       ,& 
                               todsib    &
                                )
   IMPLICIT NONE
   INTEGER      , INTENT(IN   ) :: iMax
   INTEGER      , INTENT(IN   ) :: jMax
   INTEGER      , INTENT(IN   ) :: ibMax
   INTEGER      , INTENT(IN   ) :: jbMax
   INTEGER      , INTENT(IN   ) :: kMax
   INTEGER      , INTENT(IN   ) :: ifday
   REAL(KIND=r8), INTENT(IN   ) :: tod
   INTEGER      , INTENT(IN   ) :: idate (4)
   INTEGER      , INTENT(IN   ) :: idatec(4)
   INTEGER      , INTENT(IN   ) :: ibMaxPerJB(jbmax)
   INTEGER      , INTENT(OUT  ) :: ifdy
   INTEGER      , INTENT(OUT  ) :: ids       (4)
   INTEGER      , INTENT(OUT  ) :: idc       (4)
   REAL(KIND=r8), INTENT(OUT  ) :: todsib

   INTEGER                      :: ncount
   INTEGER                      :: LRecIN
   INTEGER                      :: irec
   INTEGER                      :: ierr
   INTEGER(KIND=i4)             :: ibuf     (iMax,jMax)
   REAL(KIND=r4)                :: brf      (iMax,jMax)
   REAL(KIND=r8)                :: buf      (iMax,jMax,4)
   INTEGER(KIND=i8)             :: imask_in (iMax,jMax)
   INTEGER(KIND=i8)             :: mskant_in(iMax,jMax)
   REAL(KIND=r8)                :: VegType  (iMax,jMax,npatches) ! SIB veg type

   REAL(KIND=r8), ALLOCATABLE   :: wsib     (:,:)   
   REAL(KIND=r8), PARAMETER     :: t0 =271.17_r8
   REAL(KIND=r8)                :: sinmax
   REAL(KIND=r8), PARAMETER     :: xl0   =10.0_r8
   REAL(KIND=r8), PARAMETER     :: zero  =0.0e3_r8
   REAL(KIND=r8), PARAMETER     :: thousd=1.0e3_r8
   CHARACTER(LEN=*), PARAMETER :: h='**(InitSSiBBoundCond)**'
   INTEGER :: ier(iMax,jMax)
!   INTEGER :: ibMax
   !INTEGER :: jbMax
   INTEGER :: i
   INTEGER :: j
   ibuf=0
   brf=0.0_r4
   buf=0.0_r8
   imask_in =0_i8
   mskant_in=0_i8
   sheleg=0.0_r8
   buf=0.0_r8
   ier=0
  ! ibMax=size(imask,1)
   !jbMax=size(imask,2)
   ALLOCATE(wsib(ibMax,jbMax));wsib=0.0_r8

   CALL vegin ()
   
   INQUIRE (IOLENGTH=LRecIN) ibuf
   OPEN (UNIT=nfsibt, FILE=TRIM(fNameSibmsk),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN,&
      ACTION='READ',STATUS='OLD', IOSTAT=ierr)
   IF (ierr /= 0) THEN
      WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
              TRIM(fNameSibmsk), ierr
      STOP "**(ERROR)**"
   END IF
   brf=0.0_r4
   INQUIRE (IOLENGTH=LRecIN) brf
   OPEN (UNIT=nftgz0,FILE=TRIM(fNameTg3zrl), FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN, &
      ACTION='read', STATUS='OLD', IOSTAT=ierr) 
   IF (ierr /= 0) THEN
      WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
              TRIM(fNameTg3zrl), ierr
      STOP "**(ERROR)**"
   END IF
   
   INQUIRE (IOLENGTH=LRecIN) brf
   OPEN (UNIT=nfzol,FILE=TRIM(fNameRouLen),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIN, &
           ACTION='READ', STATUS='OLD', IOSTAT=ierr)
   IF (ierr /= 0) THEN
      WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
              TRIM(fNameRouLen), ierr
      STOP "**(ERROR)**"
   END IF

   READ(UNIT=nfsibt, REC=1) ibuf
   imask_in=ibuf
   IF (reducedGrid) THEN
      CALL FreqBoxIJtoIBJB(imask_in,iMaskSSiB)
   ELSE
      CALL IJtoIBJB( imask_in,iMaskSSiB)
   END IF
   imask=iMaskSSiB
   DO j=1,jMax
      DO i=1,iMax        
            IF (imask_in(i,j) >= 1_i8) THEN
                ier(i,j) = 0
             ELSE
                 ier(i,j) = 1
             END IF 
      END DO
      IF (ANY( ier(1:iMax,j) /= 0)) THEN
            DO i=1,iMax              
                mskant_in(i,j) = 1_i8
            END DO     
      ELSE
            DO i=1,iMax              
                mskant_in(i,j) = 0_i8
            END DO                  
      END IF              
   END DO
   IF (reducedGrid) THEN
      CALL FreqBoxIJtoIBJB(mskant_in,MskAntSSiB)
   ELSE
      CALL IJtoIBJB( mskant_in,MskAntSSiB)
   END IF
   MskAnt=MskAntSSiB       


   DO j=1,jbMax
      ncount=0
      DO i=1,ibMaxPerJB(j)
         IF(imask(i,j) > 0_i8) THEN
            ncount=ncount+1
            SoilMask(ncount,j)=textclass(bee(imask(i,j)),poros(imask(i,j) ))
         END IF
      END DO
   END DO    
   !
   !
   !         initialize sib variables
   !

   IF( ifday == 0.AND. tod == 0.0_r8 .AND. initlz >= 0 ) THEN
      call MsgOne(h,'Cold start SSib variables')
      CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin 
!climate aerosol parameters of coarse mode
       aod,asy,ssa,z_aer,ifaeros,&
!tar end
!
!tar begin 
!climate aerosol parameters of fine mode
       aodF,asyF,ssaF,z_aerF, &
!tar end  
           ifday , tod  ,idate ,idatec,&
           ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
           sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
           yrl  ,monl,ibMax,jbMax,ibMaxPerJB)
      irec=1
      CALL ReadGetNFTGZ(nftgz0,irec,buf(1:iMax,1:jMax,1),buf(1:iMax,1:jMax,2),buf(1:iMax,1:jMax,3))
      READ (UNIT=nfzol, REC=1) brf
      
      buf(1:iMax,1:jMax,4)=brf(1:iMax,1:jMax)
      
      IF (reducedGrid) THEN
         CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,1),tg1)
            !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,1),tg1)
      ELSE
         CALL IJtoIBJB(buf(1:iMax,1:jMax,1) ,tg1 )
      END IF

      IF (reducedGrid) THEN
         CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,2),tg2)
         !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,2),tg2)
      ELSE
         CALL IJtoIBJB(buf(1:iMax,1:jMax,2) ,tg2 )
      END IF

      IF (reducedGrid) THEN
         CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,3),tg3)
         !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,3),tg3)
      ELSE
         CALL IJtoIBJB(buf(1:iMax,1:jMax,3) ,tg3 )
      END IF
 
      IF (reducedGrid) THEN
         CALL LinearIJtoIBJB(buf(1:iMax,1:jMax,4),zorl)
         !CALL AveBoxIJtoIBJB(buf(1:iMax,1:jMax,4),zorl)
      ELSE
         CALL IJtoIBJB(buf(1:iMax,1:jMax,4),zorl )
      END IF
      z0    =zorl
      sinmax=150.0_r8
      !
      !     use rVisDiff as temporary for abs(soilm)
      !
      DO j=1,jbMax
         DO i=1,ibMaxPerJB(j)
            rVisDiff(i,j)=ABS(soilm(i,j))
         END DO
      END DO

      !-srf--------------------------------
      IF(iglsm_w == 0) THEN
         IF(ifslm == 5 .AND. intsoilm <=0)THEN
            DO j=1,jbMax
               ncount=0    
               DO i=1,ibMaxPerJB(j)
                  wsib(i,j)=wsib3d(i,j,3)
                  ssib(i,j)=wsib3d(i,j,3)
               END DO
            END DO
         ELSE
            CALL sibwet(ibMax,jbMax,rVisDiff,sinmax,iMaskSSiB,wsib,ssib, &
                 mxiter,ibMaxPerJB)
         END IF
      ELSE
         !
         !- rotina para chamar leitura da umidade do solo
         !
         CALL read_gl_sm_bc(imax           , & 
                            jmax           , & 
                            jbMax           , & 
                            nzg            , &
                            npatches       , &
                            npatches_actual, &
                            ibMaxPerJB     , & 
                            record_type    , & 
                            fNameSoilType  , & 
                            fNameVegType   , & 
                            fNameSoilMoist , & 
                            glsm_w         , &
                            soil_type      , &
                            veg_type       , &
                            frac_occ       , &
                            VegType             ) 
                                        
                                        
         CALL re_assign_sib_soil_prop(imax              , & ! IN
                                      jmax              , & ! IN
                                      npatches        , & ! IN
                                      imask_in        , & ! INOUT
                                      VegType                ) ! IN

         IF (reducedGrid) THEN
            CALL FreqBoxIJtoIBJB(imask_in,iMaskSSiB)
         ELSE
            CALL IJtoIBJB( imask_in,iMaskSSiB)
         END IF

         !
         !- for output isurf, use rlsm array
         !
         DO j=1,jbMax
            DO i=1,ibMaxPerJB(j)
               !  rlsm(i,j)=REAL(imask(i,j),r8)
            END DO
         END DO
         
         CALL sibwet_GLSM (ibMax              , & ! IN
                           jbMax              , & ! IN
                           iMaskSSiB          , & ! IN
                           wsib               , & ! OUT
                           ssib               , & ! IN
                           mxiter             , & ! IN
                           ibMaxPerJB         , & ! IN
                           soilm              , & ! OUT
                           nzg                , & ! in
                           wsib3d             , & ! OUT
                           glsm_w)                 ! IN

      END IF
      !-srf--------------------------------
      ppli  =0.0_r8
      ppci  =0.0_r8
      capac0=0.0_r8
      capacm=0.0_r8
      !
      !     td0 (deep soil temp) is temporarily defined as tg3
      !
      !$OMP PARALLEL DO PRIVATE(ncount,i)
      DO j=1,jbMax
         ncount=0
         DO i=1,ibMaxPerJB(j)
            gl0(i,j)=xl0
            Mmlen(i,j)=xl0
            IF(iMaskSSiB(i,j) >= 1_i8)gtsea(i,j)=290.0_r8
            IF(iMaskSSiB(i,j) >= 1_i8)gco2flx(i,j)=0.0_r8
            tseam(i,j)=gtsea(i,j)
            TSK (I,J)=ABS(gtsea(i,j))
            IF (omlmodel) THEN
                HML  (i,j) = oml_hml0 - 13.5_r8*log(MAX(ABS(TSK(i,j))-tice+0.01_r8,1.0_r8))
                HUML (I,J)=0.0_r8
                HVML (I,J)=0.0_r8
            END IF
            IF(iMaskSSiB(i,j) == 0_i8) THEN
               IF(-gtsea(i,j).LT.t0) THEN
                  iMaskSSiB(i,j)=-1_i8
                  imask(i,j)=-1_i8
               END IF
            ELSE
               ncount=ncount+1
               IF(iglsm_w == 0) THEN
                  w0        (ncount,1,j)=wsib(i,j)
                  w0        (ncount,2,j)=wsib(i,j)
                  w0        (ncount,3,j)=wsib(i,j)
               ELSE
                  !-srf--------------------------------
                  w0        (ncount,1,j)=wsib3d(i,j,1)
                  w0        (ncount,2,j)=wsib3d(i,j,2)
                  w0        (ncount,3,j)=wsib3d(i,j,3)
                  !-srf--------------------------------
               END IF
               sm0(ncount,1,j)=w0(ncount,1,j)*poros(imask(i,j))
               sm0(ncount,2,j)=w0(ncount,2,j)*poros(imask(i,j))
               sm0(ncount,3,j)=w0(ncount,3,j)*poros(imask(i,j))

               td0   (ncount,  j)=tg3 (i,j)

               IF(iglsm_w == 0) THEN
                  wm        (ncount,1,j)=wsib(i,j)
                  wm        (ncount,2,j)=wsib(i,j)
                  wm        (ncount,3,j)=wsib(i,j)
               ELSE
                  !-srf--------------------------------
                  wm        (ncount,1,j)=wsib3d(i,j,1)
                  wm        (ncount,2,j)=wsib3d(i,j,2)
                  wm        (ncount,3,j)=wsib3d(i,j,3)
                  !-srf--------------------------------
               END IF

               tdm   (ncount,  j)=tg3 (i,j)
               tgm   (ncount,  j)=tg3 (i,j)
               tcm   (ncount,  j)=tg3 (i,j)
               ssib  (ncount,j  )=0.0_r8
               IF(soilm(i,j).LT.0.0_r8)ssib(ncount,j)=wsib(i,j)

               IF(sheleg(i,j).GT.zero) THEN
                  capac0(ncount,2,j)=sheleg(i,j)/thousd
                  capacm(ncount,2,j)=sheleg(i,j)/thousd
               END IF

            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

 
 

   ELSE

      call MsgOne(h,'Warm start SSib variables')


      READ(UNIT=nfsibi) ifdy,todsib,ids,idc
      READ(UNIT=nfsibi) tm0   ,tmm
      READ(UNIT=nfsibi) qm0   ,qmm
      READ(UNIT=nfsibi) td0   ,tdm
      READ(UNIT=nfsibi) tg0   ,tgm
      READ(UNIT=nfsibi) tc0   ,tcm
      READ(UNIT=nfsibi) w0    ,wm
      READ(UNIT=nfsibi) capac0,capacm
      READ(UNIT=nfsibi) ppci  ,ppli,tkemyj
      READ(UNIT=nfsibi) gl0   ,zorl  ,gtsea ,gco2flx,tseam,qsfc0,tsfc0,qsfcm,tsfcm,HML,HUML,HVML,TSK,z0sea,&
      TC_SeaIce,TGS_SeaIce,TD_SeaIce,TA_SeaIce,SNOA_SeaIce,SNOB_SeaIce,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 
      READ(UNIT=nfsibi)  laymld,       hbath,     tdeep,sdeep
      Mmlen=gl0
      REWIND nfsibi


      IF (initlz < 0 .AND. initlz >= -3 )THEN

         IF(initlz == -2 .or. initlz == -3 )ifsst=-1
         IF(ifco2flx == -2 .or. ifco2flx == -3 )ifco2flx=-1

         CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin 
!climate aerosol parameters of coarse mode
           aod,asy,ssa,z_aer,ifaeros,&
!tar end 
!
!tar begin 
!climate aerosol parameters of fine mode
           aodF,asyF,ssaF,z_aerF, &
!tar end   
           ifday , tod  ,idate ,idatec,&
           ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
           sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
           yrl  ,monl,ibMax,jbMax,ibMaxPerJB)

         IF( initlz == -2  .or. initlz == -3 ) THEN
            !$OMP PARALLEL DO PRIVATE(i)
            DO j=1,jbMax
               DO i=1,ibMaxPerJB(j)
                  IF(iMaskSSiB(i,j) >= 1_i8) gtsea(i,j)=290.0_r8
                  IF(iMaskSSiB(i,j) >= 1_i8) gco2flx(i,j)=0.0_r8
                  !tseam(i,j) = gtsea(i,j)
                  TSK  (I,J) = ABS(gtsea(i,j))
                  IF (omlmodel) THEN
                     HML  (i,j) = oml_hml0 - 13.5_r8*log(MAX(ABS(TSK(i,j))-tice+0.01_r8,1.0_r8))
                     HUML (I,J)=0.0_r8
                     HVML (I,J)=0.0_r8
                  END IF
                  IF(iMaskSSiB(i,j) == 0_i8) THEN
                     IF(-gtsea(i,j).LT.t0) THEN
                        iMaskSSiB(i,j)=-1_i8
                        imask(i,j)=-1_i8
                     END IF
                  END IF
               END DO
            END DO
            !$OMP END PARALLEL DO
         END IF


      END IF

      !$OMP PARALLEL DO PRIVATE(ncount,i)
      DO j=1,jbMax
         ncount=0
         DO i=1,ibMaxPerJB(j)
            IF(iMaskSSiB(i,j) >= 1_i8)THEN
               ncount=ncount+1
               ssib(ncount,j)=0.0_r8
               IF(w0(ncount,1,j).LT.0.0_r8)THEN
                  ssib(ncount,j)=ABS(w0(ncount,1,j))
                  w0(ncount,1,j)=ABS(w0(ncount,1,j))
                  w0(ncount,2,j)=ABS(w0(ncount,2,j))
                  w0(ncount,3,j)=ABS(w0(ncount,3,j))
                  wm(ncount,1,j)=ABS(wm(ncount,1,j))
                  wm(ncount,2,j)=ABS(wm(ncount,2,j))
                  wm(ncount,3,j)=ABS(wm(ncount,3,j))
               END IF
               sm0(ncount,1,j)=w0(ncount,1,j)*poros(imask(i,j))
               sm0(ncount,2,j)=w0(ncount,2,j)*poros(imask(i,j))
               sm0(ncount,3,j)=w0(ncount,3,j)*poros(imask(i,j))
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

      IF(nfctrl(5).GE.1)WRITE(UNIT=nfprt,FMT=444) ifdy,todsib,ids,idc
   END IF 
   DEALLOCATE(wsib)

444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
     
  END SUBROUTINE InitSSiBBoundCond


  SUBROUTINE SSiB_Driver(&
       jdt                ,latitu               ,dtc3x            ,ncols            ,&
       nmax               ,kMax                 ,ktm              ,initlz           ,&
       kt                 ,nsx                  ,iswrad           ,ilwrad           ,&
       gt                 ,gq                   ,gu               ,gv               ,&
       prsi               ,prsl                 ,phii             ,phil             ,&
       gps                ,tmtx                 ,qmtx             ,umtx             ,&
       zenith             ,colrad               ,cos2             ,mon              ,&
       cosz               ,beam_visb            ,beam_visd        ,beam_nirb        ,&
       beam_nird          ,dlwbot               ,xvisb            ,xvisd            ,&
       xnirb              ,xnird                ,slrad            ,ppli             ,&
       ppci               ,tsea                 ,ssib             ,intg             ,&
       tseam              ,tsurf                ,qsurf            ,&
       imask              ,itype                ,tg               ,&
       ra                 ,rb                   ,rd               ,rc               ,&
       rg                 ,ta                   ,ea               ,etc              ,&
       etg                ,radt                 ,rst              ,rsoil            ,&
       ect                ,eci                  ,egt              ,egi              ,&
       egs                ,ec                   ,eg               ,hc               ,&
       hg                 ,egmass               ,etmass           ,hflux            ,&
       chf                ,shf                  ,fluxef           ,roff             ,&
       drag               ,cu                   ,ustar            ,hr               ,&
       sens               ,evap                 ,umom             ,vmom             ,&
       zorl               ,rmi                  ,rhi              ,cond             ,&
       stor               ,z0x                  ,speedm           ,Ustarm           ,&
       z0sea              ,rho                  ,d                ,qsfc             ,&
       tsfc               ,mskant               ,bstar            , &
       HML                ,HUML                 ,HVML             , &
       TSK                ,cldtot               ,ySwSfcNet        ,LwSfcNet         ,&
       pblh               ,QCF                  ,QCL              , &
       sm0                ,mlsi                 ,LwSfcDown        ,month            ,&
       Mmlen              ,idatec              ,dump)! add solange: sm0, mlsi 13-11-2012


    IMPLICIT NONE
    INTEGER      , INTENT(in   ) :: jdt
    INTEGER      , INTENT(in   ) :: latitu
    REAL(KIND=r8), INTENT(in   ) :: dtc3x
    INTEGER      , INTENT(in   ) :: nCols
    INTEGER      , INTENT(in   ) :: nmax
    INTEGER      , INTENT(in   ) :: kMax
    INTEGER      , INTENT(IN   ) :: ktm
    INTEGER      , INTENT(IN   ) :: initlz
    INTEGER      , INTENT(IN   ) :: kt
    INTEGER      , INTENT(IN   ) :: nsx(nCols)
    CHARACTER(len=*),INTENT(IN ) :: iswrad
    CHARACTER(len=*),INTENT(IN ) :: ilwrad
    INTEGER      , INTENT(IN   ) :: idatec(1:4) 
    REAL(KIND=r8), INTENT(INOUT) :: gt       (nCols,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: gq       (nCols,kmax)
    REAL(KIND=r8), INTENT(IN   ) :: gu       (nCols,kmax)
    REAL(KIND=r8), INTENT(IN   ) :: gv       (nCols,kmax)
    REAL(KIND=r8), INTENT(IN   ) :: prsi     (nCols,kMax+1)  !     prsi     - real, pressure at layer interfaces             ix,levs+1  Pa
    REAL(KIND=r8), INTENT(IN   ) :: prsl     (nCols,kMax)    !     prsl     - real, mean layer presure                       ix,levs   Pa
    REAL(KIND=r8), INTENT(IN   ) :: phii     (nCols,kMax+1)  !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8), INTENT(IN   ) :: phil     (nCols,kMax)    !===>  PHIL(K)    INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8), INTENT(IN   ) :: gps      (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: tmtx     (nCols,kmax,3)
    REAL(KIND=r8), INTENT(INOUT) :: qmtx     (nCols,kmax,3)
    REAL(KIND=r8), INTENT(INOUT) :: umtx     (nCols,kmax,4)

    REAL(KIND=r8), INTENT(IN   ) :: zenith   (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: colrad   (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: cos2     (nCols)
    INTEGER      , INTENT(inout) :: mon      (ncols)
    REAL(KIND=r8), INTENT(inout) :: cosz     (ncols)

    REAL(KIND=r8), INTENT(IN   ) :: beam_visb(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: beam_visd(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: beam_nirb(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: beam_nird(nCols)
    REAL(KIND=r8), INTENT(IN   ) :: dlwbot   (nCols)

    REAL(KIND=r8), INTENT(IN   ) :: xvisb    (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: xvisd    (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: xnirb    (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: xnird    (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: slrad    (nCols)

    REAL(KIND=r8), INTENT(IN   ) :: ppli     (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: ppci     (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: tsea     (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: ssib     (ncols)
    REAL(KIND=r8), INTENT(INOUT) :: sm0      (ncols,3)
    INTEGER      , INTENT(IN   ) :: intg

    REAL(KIND=r8), INTENT(INOUT) :: tseam    (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: tsurf    (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: qsurf    (nCols)
    
    INTEGER(KIND=i8),INTENT(IN ) :: imask    (nCols)    
    INTEGER      , INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8), INTENT(inout) :: tg   (ncols)
    REAL(KIND=r8)    ,INTENT(IN OUT) ::dump(1:nCols,1:kMax )

    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8), INTENT(inout) :: ra    (ncols)
    REAL(KIND=r8), INTENT(inout) :: rb    (ncols)
    REAL(KIND=r8), INTENT(inout) :: rd    (ncols)
    REAL(KIND=r8), INTENT(inout) :: rc    (ncols)
    REAL(KIND=r8), INTENT(inout) :: rg    (ncols)
    REAL(KIND=r8), INTENT(inout) :: ta    (ncols)
    REAL(KIND=r8), INTENT(inout) :: ea    (ncols)
    REAL(KIND=r8), INTENT(inout) :: etc   (ncols)
    REAL(KIND=r8), INTENT(inout) :: etg   (ncols)
    REAL(KIND=r8), INTENT(inout) :: radt  (ncols,icg)
    REAL(KIND=r8), INTENT(inout) :: rst   (ncols,icg)
    REAL(KIND=r8), INTENT(inout) :: rsoil (ncols)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8), INTENT(inout) :: ect   (ncols)
    REAL(KIND=r8), INTENT(inout) :: eci   (ncols)
    REAL(KIND=r8), INTENT(inout) :: egt   (ncols)
    REAL(KIND=r8), INTENT(inout) :: egi   (ncols)
    REAL(KIND=r8), INTENT(inout) :: egs   (ncols)
    REAL(KIND=r8), INTENT(inout) :: ec    (ncols)
    REAL(KIND=r8), INTENT(inout) :: eg    (ncols)
    REAL(KIND=r8), INTENT(inout) :: hc    (ncols)
    REAL(KIND=r8), INTENT(inout) :: hg    (ncols)
    REAL(KIND=r8), INTENT(inout) :: egmass(ncols)
    REAL(KIND=r8), INTENT(inout) :: etmass(ncols)
    REAL(KIND=r8), INTENT(inout) :: hflux (ncols)
    REAL(KIND=r8), INTENT(inout) :: chf   (ncols)
    REAL(KIND=r8), INTENT(inout) :: shf   (ncols)
    REAL(KIND=r8), INTENT(inout) :: fluxef(ncols)
    REAL(KIND=r8), INTENT(inout) :: roff  (ncols)
    REAL(KIND=r8), INTENT(inout) :: drag  (ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8), INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8), INTENT(inout) :: ustar (ncols)
    REAL(KIND=r8), INTENT(inout) :: hr    (ncols)




    REAL(KIND=r8), INTENT(INOUT) :: sens     (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: evap     (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: umom     (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: vmom     (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: zorl     (nCols)

    REAL(KIND=r8), INTENT(INOUT) :: rmi      (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: rhi      (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: cond     (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: stor     (nCols)

    REAL(KIND=r8), INTENT(INOUT) :: z0x      (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: speedm   (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: Ustarm   (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: z0sea    (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: rho      (nCols)
    REAL(KIND=r8), INTENT(INOUT) :: d        (ncols)
    REAL(KIND=r8), INTENT(INOUT) :: qsfc     (ncols)
    REAL(KIND=r8), INTENT(INOUT) :: tsfc     (ncols)
    INTEGER(KIND=i8),INTENT(IN ) :: mskant   (ncols)
    INTEGER(KIND=i8), INTENT(INOUT) :: mlsi  (ncols)
    REAL(KIND=r8),INTENT(out   ) :: bstar  (ncols)

    REAL(KIND=r8),    INTENT(INOUT) :: HML  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HUML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HVML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: TSK  (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cldtot (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ySwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: LwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCF(nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCL(nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: LwSfcDown(1:nCols )
    INTEGER         , INTENT(IN   ) :: month(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: Mmlen (1:nCols)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8) :: qm  (ncols)
    REAL(KIND=r8) :: tm  (ncols)
    REAL(KIND=r8) :: um  (ncols)
    REAL(KIND=r8) :: vm  (ncols)
    REAL(KIND=r8) :: psur(ncols)
    REAL(KIND=r8) :: ppc (ncols)
    REAL(KIND=r8) :: ppl (ncols)
    REAL(KIND=r8) :: radn(ncols,3,2)
    !
    !     prognostic variables
    !
    REAL(KIND=r8) :: tc   (ncols)
    REAL(KIND=r8) :: td   (ncols)
    REAL(KIND=r8) :: capac(ncols,2)
    REAL(KIND=r8) :: w    (ncols,3)
    REAL(KIND=r8) :: tcta  (ncols)
    REAL(KIND=r8) :: tgta  (ncols)
    REAL(KIND=r8) :: btc   (ncols)
    REAL(KIND=r8) :: btg   (ncols)
    REAL(KIND=r8) :: u2    (ncols)
    REAL(KIND=r8) :: par   (ncols,icg)
    REAL(KIND=r8) :: pd    (ncols,icg)

    REAL(KIND=r8) :: phroot(ncols,icg)
    REAL(KIND=r8) :: hrr   (ncols)
    REAL(KIND=r8) :: phsoil(ncols,idp)
    REAL(KIND=r8) :: cc    (ncols)
    REAL(KIND=r8) :: cg    (ncols)
    REAL(KIND=r8) :: satcap(ncols,icg)
    REAL(KIND=r8) :: snow  (ncols,icg)
    REAL(KIND=r8) :: dtc   (ncols)
    REAL(KIND=r8) :: dtg   (ncols)
    REAL(KIND=r8) :: dtm   (ncols)
    REAL(KIND=r8) :: dqm   (ncols)
    REAL(KIND=r8) :: stm   (ncols,icg)
    REAL(KIND=r8) :: extk  (ncols,icg,iwv,ibd)
    REAL(KIND=r8) :: radfac(ncols,icg,iwv,ibd)
    REAL(KIND=r8) :: closs (ncols)
    REAL(KIND=r8) :: gloss (ncols)
    REAL(KIND=r8) :: thermk(ncols)
    REAL(KIND=r8) :: p1f   (ncols)
    REAL(KIND=r8) :: p2f   (ncols)

    REAL(KIND=r8) :: ecidif(ncols)
    REAL(KIND=r8) :: egidif(ncols)
    REAL(KIND=r8) :: ecmass(ncols)

    REAL(KIND=r8) :: bps   (ncols)
    REAL(KIND=r8) :: psb   (ncols)
    REAL(KIND=r8) :: dzm   (ncols)
    REAL(KIND=r8) :: em    (ncols)
    REAL(KIND=r8) :: gmt   (ncols,3)
    REAL(KIND=r8) :: gmq   (ncols,3)
    REAL(KIND=r8) :: gmu   (ncols,4)

    REAL(KIND=r8) :: cuni  (ncols)
    REAL(KIND=r8) :: ctni  (ncols)

    REAL(KIND=r8) :: sinclt(ncols)
    REAL(KIND=r8) :: rhoair(ncols)
    REAL(KIND=r8) :: psy   (ncols)
    REAL(KIND=r8) :: rcp   (ncols)
    REAL(KIND=r8) :: wc    (ncols)
    REAL(KIND=r8) :: wg    (ncols)
    REAL(KIND=r8) :: fc    (ncols)
    REAL(KIND=r8) :: fg    (ncols)

    REAL(KIND=r8) :: zlwup_local    (ncols)
    REAL(KIND=r8) :: salb     (ncols,2,2)
    REAL(KIND=r8) :: tgeff    (ncols)
    
    REAL(KIND=r8) :: rstpar2 (ncols,icg,iwv)
    REAL(KIND=r8) :: zlt2    (ncols,icg)
    REAL(KIND=r8) :: green2  (ncols,icg)
    REAL(KIND=r8) :: chil2   (ncols,icg)
    REAL(KIND=r8) :: vcover  (ncols,icg)
    REAL(KIND=r8) :: rdc     (ncols)
    REAL(KIND=r8) :: rbc     (ncols)
    REAL(KIND=r8) :: z0      (ncols)
    REAL(KIND=r8) :: topt2   (ncols,icg)
    REAL(KIND=r8) :: tll2    (ncols,icg)
    REAL(KIND=r8) :: tu2     (ncols,icg)
    REAL(KIND=r8) :: defac2  (ncols,icg)
    REAL(KIND=r8) :: xsea    (nCols)
    REAL(KIND=r8) :: tmin    (ncols)
    REAL(KIND=r8) :: tmax    (ncols)
    REAL(KIND=r8) :: z0xloc  (ncols)
    REAL(KIND=r8) :: ph22    (ncols,icg)
    REAL(KIND=r8) :: ph12    (ncols,icg)
    REAL(KIND=r8) :: bstar1  (ncols)
    REAL(KIND=r8) :: cpsy
    REAL(KIND=r8) :: r100
    REAL(KIND=r8) :: dzcut
    LOGICAL       :: InitMod
    REAL(KIND=r8)   :: GSW (nCols)
    REAL(KIND=r8)   :: GLW (nCols)

    INTEGER :: ncount,i,n,m,k,ll,j,itr,ind,nint,IntSib
    pd=0.0_r8
    phsoil=0.0_r8
    PHROOT=0.0_r8
    PAR=0.0_r8
    SNOW=0.0_r8
    STM=0.0_r8
    !
    !            cp
    ! cpsy = ------------
    !          L*epsfac
    !
    cpsy=cp/(hl*epsfac)
    r100=100.0e0_r8 /gasr
    DO i=1,nCols
       GSW(I) = xvisb  (i)+xvisd  (i)+xnirb  (i)+xnird  (i)
       GLW(I) = dlwbot  (i)
    END DO

    DO j=1,icg
       DO i=1,nmax
          vcover     (i,j)   =  vcover_gbl (i,latitu,j)
          zlt2       (i,j)   =  zlt_gbl    (i,latitu,j)
          green2     (i,j)   =  green_gbl  (i,latitu,j)
          chil2      (i,j)   =  chil_gbl   (i,latitu,j)
          topt2      (i,j)   =  topt_gbl   (i,latitu,j)
          tll2       (i,j)   =  tll_gbl    (i,latitu,j)
          tu2        (i,j)   =  tu_gbl     (i,latitu,j)
          defac2     (i,j)   =  defac_gbl  (i,latitu,j)
          ph12       (i,j)   =  ph1_gbl    (i,latitu,j)
          ph22       (i,j)   =  ph2_gbl    (i,latitu,j)
          rstpar2    (i,j,1) =  rstpar_gbl (i,latitu,j,1)
          rstpar2    (i,j,2) =  rstpar_gbl (i,latitu,j,2)
          rstpar2    (i,j,3) =  rstpar_gbl (i,latitu,j,3)
          satcap     (i,j)   =  satcap3d   (i,j,latitu)
          !snow       (i,j)   =  snow_gbl   (i,latitu,j)
       END DO
    END DO    
    
    DO i=1,nmax  
       closs(i)    = closs_gbl(i,latitu)
       gloss(i)    = gloss_gbl(i,latitu)
       thermk(i)   = thermk_gbl(i,latitu)
       p1f  (i)    = p1f_gbl(i,latitu)  
       p2f  (i)    = p2f_gbl(i,latitu)  
       zlwup_local(i)    = zlwup_SSiB(i,latitu)
       salb  (i,1,1) = AlbGblSSiB(i,1,1,latitu)
       salb  (i,1,2) = AlbGblSSiB(i,1,2,latitu)
       salb  (i,2,1) = AlbGblSSiB(i,2,1,latitu)
       salb  (i,2,2) = AlbGblSSiB(i,2,2,latitu)
       tgeff (i)     = tgeff_gbl(i,latitu)
    END DO

    DO m=1,ibd
       DO n=1,iwv
          DO ll=1,icg
             DO i=1,nmax
                extk  (i,ll,n,m) = extk_gbl  (i,ll,n,m,latitu)
                radfac(i,ll,n,m) = radfac_gbl(i,ll,n,m,latitu)
             END DO
          END DO
       END DO
    END DO    
    !
    ! Forcing the soil moisture with observation data [specific case]
    !
    IF(ifslm == 5) THEN
      DO k=1,3
         ncount=0
         DO i=1,nCols
            IF(imask(i) >= 1_i8)THEN
               ncount=ncount+1
               wm(ncount,k,latitu)=wsib3d(i,j,k)
               w0(ncount,k,latitu)=wsib3d(i,j,k)
            END IF  
         END DO
       END DO
    END IF
    !
    !  Due implicit scheme of the ssib. Initializate the prognostics variable with time step t-1
    ! 
    DO i=1,nmax
          td   (i)    = tdm   (i,latitu)
          tg   (i)    = tgm   (i,latitu)
          tc   (i)    = tcm   (i,latitu)
          capac(i,1)  = capacm(i,1,latitu)
          capac(i,2)  = capacm(i,2,latitu)
          w    (i,1)  = wm    (i,1,latitu)
          w    (i,2)  = wm    (i,2,latitu)
          w    (i,3)  = wm    (i,3,latitu)
    END DO

    ncount=0
    DO i=1,nCols
       IF(imask(i).GE.1_i8) THEN
          ncount=ncount+1
          mlsi(i) = 1_i8   !add solange 13-11-2012
          sinclt(ncount) = SIN(colrad(i))
          psur(ncount)=gps(i)
          tm  (ncount)=gt (i,1)
          qm  (ncount)=gq (i,1)
          um  (ncount)=gu (i,1)/SIN( colrad(i))
          vm  (ncount)=gv (i,1)/SIN( colrad(i))
          psy (ncount)=cpsy*psur(ncount)
          !
          ! P =rho*R*T   and     P = rho*g*Z
          !
          !                           P
          ! DP = rho*g*DZ and rho = ----
          !                          R*T
          !        P
          ! DP = ----*g*DZ
          !       R*T
          !
          !       RT      DP
          ! DZ = ---- * -----
          !       g       P
          !
          !              RT     P2 - P1
          ! Z2 - Z1 =  ----- *---------
          !              g        P2
          !
          !
          !              R      si(k) - si(k+1)
          ! dzm   (i) = --- * (----------------) * tm(i)
          !              g            2
          !
          !            (J/kg/K)
          ! dzm   (i) =--------*K
          !              m/s^2
          !
          !                J*s^2
          ! dzm   (i) = -----------*K
          !               kg*K * m
          !
          !               (kg*m*s^-2*m)*s^2
          ! dzm   (i) = ---------------------*K
          !                  kg*K * m
          !
          !               (m*m)
          ! dzm   (i) = --------- = m
          !                 m
          !dzm   (ncount)=rbyg*tm(ncount)           
          dzm   (ncount)=0.5_r8*MAX((phii(i,2) - phii(i,1)),0.5_r8)
          !
          ! presure of vapor at atmosphere
          !
          em    (ncount)=qm(ncount)*psur(ncount)/(epsfac+qm(ncount))
          !
          ! Factor conversion to potention temperature
          !
          !bps   (ncount)=sigki(1)
          bps   (ncount)= (prsi(i,1)/(prsl(i,1)))**(gasr/cp)
          !
          ! Difference of pressure
          !
          psb   (ncount)=psur(ncount)*(prsi(i,1)/prsi(i,1)) -(prsi(i,2)/prsi(i,1) )
          !
          !Density of air
          !
          !
          ! P =rho*R*T
          !
          !        P
          !rho =-------
          !       R*T
          !
          !               1       100*psur         1             Pa
          !rhoair(i) = ------- * ---------- = ----------- *-----------------
          !               R         Tm         (J/kg/K)           K
          !
          !                kg*K            N/m^2
          !rhoair(i) =  ----------- *-----------------
          !                  J              K
          !
          !                kg*K               kg*m*s^-2* m^-2
          !rhoair(i) =  ---------------- * -----------------
          !               (kg*m*s^-2*m)          K
          !
          !                 kg
          !rhoair(i) =  --------
          !                 m^3
          !
          rhoair(ncount)=r100*psur(ncount)/tm(ncount)
          !
          !         J          kg          J
          ! rcp = -------- * ------- = ----------
          !        kg * K      m^3       K * m^3
          !
          rcp   (ncount)=cp*rhoair(ncount)
          !
          gmt             (ncount,1)=tmtx(i,1,1)
          gmt             (ncount,2)=tmtx(i,1,2)
          gmt             (ncount,3)=tmtx(i,1,3)
          gmq             (ncount,1)=qmtx(i,1,1)
          gmq             (ncount,2)=qmtx(i,1,2)
          gmq             (ncount,3)=qmtx(i,1,3)
          gmu             (ncount,1)=umtx(i,1,1)
          gmu             (ncount,2)=umtx(i,1,2)
          gmu             (ncount,3)=umtx(i,1,3)
          gmu             (ncount,4)=umtx(i,1,4)
          rbc        (ncount)     =  xbc   (itype(ncount),mon(ncount))
          rdc        (ncount)     =  xdc   (itype(ncount),mon(ncount))
!          z0x        (ncount)     =  x0x   (itype(ncount),mon(ncount))
          z0xloc     (ncount)     =  x0x   (itype(ncount),mon(ncount))    ! tvsgm - Global Mean Surface Virtual Temperature
          IF(itype(ncount) == 1) THEN
             ! dz - mean height of the first model layer
             !tvsgm=288.16_r8
             !dz=(gasr*tvsgm/grav)*LOG(si1/sl1)
             ! Forest
             !dzcut=0.75_r8*dz
             dzcut=0.6_r8*dzm(ncount)
             d (ncount)     = MIN(xd    (itype(ncount),mon(ncount)),dzcut)
             IF((dzm(ncount) - d (ncount))  < 1.0_r8)dzm(ncount) = d (ncount) + 1.0_r8
             
             !xd(1,1:imon)=MIN(xd(1,1:imon),dzcut)
          ELSE
            ! Other
            ! SiB calibration values
            ! 45 m - height of the first tower level of measurements
            ! 27 m - maximum calibrated displacement height
            dzcut=(27.0_r8/45.0_r8)*dzm(ncount)
            d (ncount)     = MIN(xd    (itype(ncount),mon(ncount)),dzcut)
            IF((dzm(ncount) - d (ncount))  < 1.0_r8)dzm(ncount) = d (ncount) + 1.0_r8
           !xd(2:ityp,1:imon)=MIN(xd(2:ityp,1:imon),dzcut)
          END IF
          !d          (ncount)     =  xd    (itype(ncount),mon(ncount))
       END IF
    END DO


    InitMod = (initlz >= 0 .AND. ktm == -1 .AND. kt == 0 .AND. nmax >= 1)

    IF(InitMod)THEN
       nint=2
       IntSib=5
    ELSE
       nint=1
       IntSib=1
    END IF

    IF(TRIM(iswrad).NE.'NON'.AND.TRIM(ilwrad).NE.'NON') THEN
       IF(InitMod)THEN

          DO ind=1,nint
             ncount=0
             DO i=1,nCols
                IF(imask(i).GE.1_i8) THEN
                   ncount=ncount+1
                   IF(ind.EQ.1) THEN
                      !
                      !     night
                      !
                      radn(ncount,1,1)=0.0e0_r8
                      radn(ncount,1,2)=0.0e0_r8
                      radn(ncount,2,1)=0.0e0_r8
                      radn(ncount,2,2)=0.0e0_r8
                      cosz(ncount)    =0.0e0_r8
                   ELSE
                      !
                      !     noon
                      !
                      radn(ncount,1,1)=beam_visb (i)  ! radn(1,1)=!Downward Surface shortwave fluxe visible beam (cloudy)
                      radn(ncount,1,2)=beam_visd (i)  ! radn(1,2)=!Downward Surface shortwave fluxe visible diffuse (cloudy)
                      radn(ncount,2,1)=beam_nirb (i)  ! radn(2,1)=!Downward Surface shortwave fluxe Near-IR beam (cloudy)
                      radn(ncount,2,2)=beam_nird (i)  ! radn(2,2)=!Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
                      cosz(ncount)    =cos2(i)
                   END IF
                   radn(ncount,3,1)=0.0e0_r8
                   radn(ncount,3,2)=dlwbot(i)
                   !
                   !     precipitation
                   !
                   ppl (ncount)    =0.0e0_r8
                   ppc (ncount)    =0.0e0_r8
                END IF
             END DO
             DO itr=1,IntSib
                CALL radalb( &
                     nmax              ,mon(1:nmax)         ,nmax                ,itype(1:nmax)       , &
                     tc(1:nmax)        ,tg(1:nmax)          ,capac(1:nmax,:)     , &
                     satcap(1:nmax,:)  ,extk(1:nmax,:,:,:)  ,radfac(1:nmax,:,:,:),closs(1:nmax)       , &
                     gloss(1:nmax)     ,thermk(1:nmax)      ,p1f(1:nmax)         ,p2f(1:nmax)         , &
                     zlwup_local(1:nmax),salb(1:nmax,:,:)    ,tgeff(1:nmax)       ,cosz(1:nmax)        , &
                     nsx  (1:nmax)      ,latitu  )

                CALL fysiks(&
                     vcover(1:nmax,:)  ,z0xloc(1:nmax)      ,d(1:nmax)       ,rdc(1:nmax)     ,&
                     rbc(1:nmax)       ,z0(1:nmax)          ,jdt             ,latitu          ,&
                     bps(1:nmax)       ,psb(1:nmax)         ,dzm(1:nmax)     ,em(1:nmax)      ,&
                     gmt(1:nmax,:)     ,gmq(1:nmax,:)       ,gmu(1:nmax,:)   ,cu(1:nmax)      ,&
                     cuni(1:nmax)      ,ctni(1:nmax)        ,ustar(1:nmax)   ,cosz(1:nmax)    ,&
                     sinclt(1:nmax)    ,rhoair(1:nmax)      ,psy(1:nmax)     ,rcp(1:nmax)     ,&
                     wc(1:nmax)        ,wg(1:nmax)          ,fc(1:nmax)      ,fg(1:nmax)      ,&
                     hr(1:nmax)        ,ect(1:nmax)         ,eci(1:nmax)     ,egt(1:nmax)     ,&
                     egi(1:nmax)       ,egs(1:nmax)         ,ec(1:nmax)      ,eg(1:nmax)      ,&
                     hc(1:nmax)        ,hg(1:nmax)          ,ecidif(1:nmax)  ,egidif(1:nmax)  ,&
                     ecmass(1:nmax)    ,egmass(1:nmax)      ,etmass(1:nmax)  ,hflux(1:nmax)   ,&
                     chf(1:nmax)       ,shf(1:nmax)         ,fluxef(1:nmax)  ,roff(1:nmax)    ,&
                     drag(1:nmax)      ,ra(1:nmax)          ,rb(1:nmax)      ,rd(1:nmax)      ,&
                     rc(1:nmax)        ,rg(1:nmax)          ,tcta(1:nmax)    ,tgta(1:nmax)    ,&
                     ta(1:nmax)        ,ea(1:nmax)          ,etc(1:nmax)     ,etg(1:nmax)     ,&
                     btc(1:nmax)       ,btg(1:nmax)         ,u2(1:nmax)      ,radt(1:nmax,:)  ,&
                     par(1:nmax,:)     ,pd(1:nmax,:)        ,rst(1:nmax,:)   ,rsoil(1:nmax)   ,&
                     phroot(1:nmax,:)  ,hrr(1:nmax)         ,phsoil(1:nmax,:),cc(1:nmax)      ,&
                     cg(1:nmax)        ,satcap(1:nmax,:)    ,snow(1:nmax,:)  ,dtc(1:nmax)     ,&
                     dtg(1:nmax)       ,dtm(1:nmax)         ,dqm(1:nmax)     ,stm(1:nmax,:)   ,&
                     extk(1:nmax,:,:,:),radfac(1:nmax,:,:,:),closs(1:nmax)   ,gloss(1:nmax)   ,&
                     thermk(1:nmax)    ,p1f(1:nmax)         ,p2f(1:nmax)     ,tc(1:nmax)      ,&
                     tg(1:nmax)        ,td(1:nmax)          ,capac(1:nmax,:) ,w(1:nmax,:)     ,&
                     qm(1:nmax)        ,tm(1:nmax)          ,um(1:nmax)      ,vm(1:nmax)      ,&
                     psur(1:nmax)      ,ppc(1:nmax)         ,ppl(1:nmax)     ,radn(1:nmax,:,:),&
                     itype(1:nmax)     ,dtc3x               ,mon (1:nmax)    ,nmax            ,&
                     nmax              ,zlt2(1:nmax,:)      ,green2(1:nmax,:),chil2(1:nmax,:) ,&
                     rstpar2(1:nmax,:,:),topt2(1:nmax,:)    ,tll2(1:nmax,:)  ,tu2(1:nmax,:)   ,&
                     defac2(1:nmax,:)  ,ph12(1:nmax,:)      ,ph22(1:nmax,:)  ,bstar1(1:nmax))
                ncount=0
                DO i=1,nCols
                   IF(imask(i).GE.1_i8) THEN
                      ncount=ncount+1
                      tm (ncount  )=gt  (i,1)
                      qm (ncount  )=gq  (i,1)
                      gmt(ncount,1)=tmtx(i,1,1)
                      gmt(ncount,2)=tmtx(i,1,2)
                      gmt(ncount,3)=tmtx(i,1,3)
                      gmq(ncount,1)=qmtx(i,1,1)
                      gmq(ncount,2)=qmtx(i,1,2)
                      gmq(ncount,3)=qmtx(i,1,3)
                      gmu(ncount,1)=umtx(i,1,1)
                      gmu(ncount,2)=umtx(i,1,2)
                      gmu(ncount,3)=umtx(i,1,3)
                      gmu(ncount,4)=umtx(i,1,4)
                   END IF
                END DO
             END DO
             DO i=1,nmax
                capac(i,1)=capacm(i,1,latitu)
                capac(i,2)=capacm(i,2,latitu)
                w    (i,1)=wm    (i,1,latitu)
                w    (i,2)=wm    (i,2,latitu)
                w    (i,3)=wm    (i,3,latitu)
                td   (i)  =tdm   (i,latitu)
                tc   (i)  =tcm   (i,latitu)
                IF(ind.EQ.1) THEN
                   tmin (i) =tg (i)
                ELSE
                   tmax (i) =tg (i)
                END IF
                tg   (i) =tgm(i,latitu)
             END DO
          END DO
          DO i=1,nmax
             td   (i) =0.9_r8*0.5_r8*(tmax(i)+tmin(i))+0.1_r8*tdm(i,latitu)
             tdm  (i,latitu) =td(i)
             td0  (i,latitu) =td(i)
          END DO
          !
          !     this is a start of equilibrium tg,tc comp.
          !
          ncount=0
          DO i=1,nCols
             IF(imask(i).GE.1_i8) THEN
                ncount=ncount+1
                cosz(ncount)    =zenith(i)
             END IF
          END DO
          DO i=1,nmax
             IF(cosz(i).LT.0.0e0_r8) THEN
                tgm  (i,latitu)  =tmin(i)
                tg0  (i,latitu)  =tmin(i)
             END IF
          END DO
          CALL radalb ( &
               nmax              ,mon(1:nmax)         ,nmax                ,itype(1:nmax)       , &
               tc(1:nmax)        ,tg(1:nmax)          ,capac(1:nmax,:)     , &
               satcap(1:nmax,:)  ,extk(1:nmax,:,:,:)  ,radfac(1:nmax,:,:,:),closs(1:nmax)       , &
               gloss(1:nmax)     ,thermk(1:nmax)      ,p1f(1:nmax)         ,p2f(1:nmax)         , &
               zlwup_local(1:nmax),salb(1:nmax,:,:)    ,tgeff(1:nmax)       ,cosz(1:nmax)        , &
               nsx (1:nmax)              ,latitu    )
       END IF
    END IF
    IF(nmax.GE.1) THEN
       ncount=0
       DO i=1,nCols
          IF(imask(i).GE.1_i8) THEN
             ncount=ncount+1
             !
             !     this is for radiation interpolation
             !
             IF(cosz(ncount).GE.0.01746e0_r8 ) THEN
                radn(ncount,1,1)=xvisb (i)
                radn(ncount,1,2)=xvisd (i)
                radn(ncount,2,1)=xnirb (i)
                radn(ncount,2,2)=xnird (i)
             ELSE
                radn(ncount,1,1)=0.0e0_r8
                radn(ncount,1,2)=0.0e0_r8
                radn(ncount,2,1)=0.0e0_r8
                radn(ncount,2,2)=0.0e0_r8
             END IF
             radn(ncount,3,1)=0.0e0_r8
             radn(ncount,3,2)=dlwbot(i)
             !
             !     precipitation
             !
             ppl (ncount)    =ppli  (i)
             ppc (ncount)    =ppci  (i)
          END IF
       END DO
       CALL fysiks(&
            vcover(1:nmax,:)  ,z0xloc(1:nmax)      ,d(1:nmax)       ,rdc(1:nmax)     ,&
            rbc(1:nmax)       ,z0(1:nmax)          ,jdt             ,latitu          ,&
            bps(1:nmax)       ,psb(1:nmax)         ,dzm(1:nmax)     ,em(1:nmax)      ,&
            gmt(1:nmax,:)     ,gmq(1:nmax,:)       ,gmu(1:nmax,:)   ,cu(1:nmax)      ,&
            cuni(1:nmax)      ,ctni(1:nmax)        ,ustar(1:nmax)   ,cosz(1:nmax)    ,&
            sinclt(1:nmax)    ,rhoair(1:nmax)      ,psy(1:nmax)     ,rcp(1:nmax)     ,&
            wc(1:nmax)        ,wg(1:nmax)          ,fc(1:nmax)      ,fg(1:nmax)      ,&
            hr(1:nmax)        ,ect(1:nmax)         ,eci(1:nmax)     ,egt(1:nmax)     ,&
            egi(1:nmax)       ,egs(1:nmax)         ,ec(1:nmax)      ,eg(1:nmax)      ,&
            hc(1:nmax)        ,hg(1:nmax)          ,ecidif(1:nmax)  ,egidif(1:nmax)  ,&
            ecmass(1:nmax)    ,egmass(1:nmax)      ,etmass(1:nmax)  ,hflux(1:nmax)   ,&
            chf(1:nmax)       ,shf(1:nmax)         ,fluxef(1:nmax)  ,roff(1:nmax)    ,&
            drag(1:nmax)      ,ra(1:nmax)          ,rb(1:nmax)      ,rd(1:nmax)      ,&
            rc(1:nmax)        ,rg(1:nmax)          ,tcta(1:nmax)    ,tgta(1:nmax)    ,&
            ta(1:nmax)        ,ea(1:nmax)          ,etc(1:nmax)     ,etg(1:nmax)     ,&
            btc(1:nmax)       ,btg(1:nmax)         ,u2(1:nmax)      ,radt(1:nmax,:)  ,&
            par(1:nmax,:)     ,pd(1:nmax,:)        ,rst(1:nmax,:)   ,rsoil(1:nmax)   ,&
            phroot(1:nmax,:)  ,hrr(1:nmax)         ,phsoil(1:nmax,:),cc(1:nmax)      ,&
            cg(1:nmax)        ,satcap(1:nmax,:)    ,snow(1:nmax,:)  ,dtc(1:nmax)     ,&
            dtg(1:nmax)       ,dtm(1:nmax)         ,dqm(1:nmax)     ,stm(1:nmax,:)   ,&
            extk(1:nmax,:,:,:),radfac(1:nmax,:,:,:),closs(1:nmax)   ,gloss(1:nmax)   ,&
            thermk(1:nmax)    ,p1f(1:nmax)         ,p2f(1:nmax)     ,tc(1:nmax)      ,&
            tg(1:nmax)        ,td(1:nmax)          ,capac(1:nmax,:) ,w(1:nmax,:)     ,&
            qm(1:nmax)        ,tm(1:nmax)          ,um(1:nmax)      ,vm(1:nmax)      ,&
            psur(1:nmax)      ,ppc(1:nmax)         ,ppl(1:nmax)     ,radn(1:nmax,:,:),&
            itype(1:nmax)     ,dtc3x               ,mon (1:nmax)    ,nmax            ,&
            nmax              ,zlt2(1:nmax,:)      ,green2(1:nmax,:),chil2(1:nmax,:) ,&
            rstpar2(1:nmax,:,:),topt2(1:nmax,:)    ,tll2(1:nmax,:)  ,tu2(1:nmax,:)   ,&
            defac2(1:nmax,:)  ,ph12(1:nmax,:)      ,ph22(1:nmax,:)  ,bstar1(1:nmax))
    END IF
    !     
    !     temperature and snow depths in Antarctica and Groenland
    !     
    ncount=0
    DO i=1,nCols
       IF(imask(i).GE.1_i8 ) THEN
          ncount=ncount+1
          IF ( imask(i).EQ.13_i8 ) THEN
            w (ncount,1)  = 1.0_r8
            w (ncount,2)  = 1.0_r8
            w (ncount,3)  = 1.0_r8
            TD(ncount  )  = MAX(MIN(TD(ncount) ,273.15_r8),218.15_r8)
            TC(ncount  )  = MAX(MIN(TC(ncount) ,273.15_r8),218.15_r8)
            tg(ncount  )  = MAX(MIN(tg(ncount) ,273.15_r8),218.15_r8)
          END IF
       END IF
    END DO

    !
    !     sib time integaration and time filter
    !
    DO i=1,nmax
       !tm(i)=ABS(ta(i))/bps(i)
       !qm(i)=0.622e0_r8*EXP(21.65605e0_r8 -5418.0e0_r8 /tm(i))/gps(i)
       qm(i)=MAX(1.0e-12_r8,qm(i))
    END DO
    CALL sextrp ( &
         td(1:nCols)               ,tg(1:nCols)               ,tc(1:nmax)                ,w(1:nCols,1:3)        , &
         capac(1:nCols,1:2)        ,td0(1:nCols,latitu)       ,tg0(1:nCols,latitu)       ,tc0 (1:nCols,latitu)  , &
         w0(1:nCols,1:3,latitu)    ,capac0(1:nCols,1:2,latitu),tdm(1:nCols,latitu)       ,tgm(1:nCols,latitu)   , &
         tcm(1:nCols,latitu)       ,wm(1:nCols,1:3,latitu)    ,capacm(1:nCols,1:2,latitu),istrt                 , &
         ncols                     ,nmax                      ,epsflt                    ,intg                  , &
         latitu                    ,tm0(1:nCols,latitu)       ,qm0(1:nCols,latitu)       ,tm(1:nCols)           , &
         qm(1:nCols)               ,tmm(1:nCols,latitu)       ,qmm(1:nCols,latitu)    )
    !
    !     fix soil moisture at selected locations
    !
    DO i=1,nmax
       IF(ssib(i).GT.0.0_r8)THEN
          qm(i)=MAX(1.0e-12_r8,qm(i))
          w0(i,1,latitu)=ssib(i)
          w0(i,2,latitu)=ssib(i)
          w0(i,3,latitu)=ssib(i)
          wm(i,1,latitu)=ssib(i)
          wm(i,2,latitu)=ssib(i)
          wm(i,3,latitu)=ssib(i)
       END IF
    END DO
    ncount=0
    DO i=1,nCols
       IF(imask(i).GE.1_i8) THEN
          ncount=ncount+1
          tmtx(i,1,3)=gmt(ncount,3)
          qmtx(i,1,3)=gmq(ncount,3)
          umtx(i,1,3)=gmu(ncount,3)
          umtx(i,1,4)=gmu(ncount,4)
          tsea(i)    =tgeff(ncount)
          TSK (i)=tgeff(ncount)
          z0x(ncount)=z0(ncount)
       END IF
    END DO
    ncount=0
    DO i=1,nCols
       IF(imask(i).GE.1_i8 ) THEN
          ncount=ncount+1
          IF ( imask(i).EQ.13_i8 ) THEN
             sm0 (ncount,1)    = poros (imask(i))
          ELSE
             sm0 (ncount,1)    = w0(ncount,1,latitu)* poros (imask(i))
             sm0 (ncount,2)    = w0(ncount,2,latitu)* poros (imask(i))
             sm0 (ncount,3)    = w0(ncount,3,latitu)* poros (imask(i))
          END IF
       END IF
    END DO
    !
    !     sea or sea ice
    ! gu gv gps colrad sigki delsig sens evap umom vmom rmi rhi cond stor zorl rnet ztn2 THETA_2M VELC_2m MIXQ_2M
    ! THETA_10M VELC_10M MIXQ_10M
    ! including case 1D physics
    DO i=1,nCols
       IF(MskAnt(i) == 1_i8)THEN
          xsea (i) = tseam(i)
          tsfc (i) = tsfcm(i,latitu)
          qsfc (i) = qsfcm(i,latitu)
       END IF   
    END DO

   CALL seasfc( &
           tmtx  (1:nCols,1:kMax,1:3)  ,umtx  (1:nCols,1:kMax,1:4),qmtx  (1:nCols,1:kMax,1:3)  ,&
           kmax                        ,kmax                      ,slrad (1:nCols)             ,&
           tsurf (1:nCols)             ,qsurf (1:nCols)           ,gu    (1:nCols,1:kMax)      ,&
           gv    (1:nCols,1:kMax)      ,gt    (1:nCols,1:kMax)    ,gq    (1:nCols,1:kMax)      ,&
           prsi(1:ncols,1:kMax+1),prsl(1:ncols,1:kMax),phii(1:nCols,1:kMax+1),phil(1:nCols,1:kMax),&
           xsea  (1:nCols)           ,dtc3x                       ,SIN(colrad(1:nCols))        ,&
           sens  (1:nCols)             ,evap  (1:nCols)           ,umom  (1:nCols)             ,&
           vmom  (1:nCols)             ,rmi   (1:nCols)           ,rhi   (1:nCols)             ,&
           cond  (1:nCols)             ,stor  (1:nCols)           ,zorl  (1:nCols)             ,&
           nCols                       ,speedm(1:nCols)           ,bstar (1:nCols)             ,&
           Ustarm(1:nCols)             ,z0sea (1:nCols)           ,rho   (1:nCols)             ,&
           qsfc  (1:nCols)             ,tsfc  (1:nCols)           ,MskAnt(1:nCols)             ,&
           iMask (1:nCols)             ,zenith (1:nCols)          ,ppli  (1:nCols)             ,&
           ppci  (1:nCols)             ,LwSfcDown(1:nCols)        ,xvisb (1:nCols)             ,&
           xvisd (1:nCols)             ,xnirb(1:nCols)            ,xnird (1:nCols)             ,&
           HML   (1:nCols)             ,HUML (1:nCols)            ,HVML (1:nCols)              ,&
           TSK   (1:nCols)             ,GSW(1:nCols)              ,GLW(1:nCols)                ,&
           cldtot(1:nCols,1:kMax)      ,ySwSfcNet(1:nCols)        ,month(1:nCols)             ,& 
           LwSfcNet(1:nCols)           ,pblh  (1:nCols)           ,QCF (1:nCols,1:kMax)        ,&
           QCL  (1:nCols,1:kMax)       ,mlsi  (1:nCols)           ,latitu                      ,&
           Mmlen (1:nCols)             ,colrad(1:nCols)           ,idatec ,dump(1:nCols,1:kMax ))


    DO i=1,nCols
       IF(MskAnt(i) == 1_i8 .and. tsea(i).LE.0.0e0_r8.AND.tsurf(i).LT.tice+0.01e0_r8 ) THEN
              IF(intg.EQ.2) THEN
                 IF(istrt.EQ.0) THEN
                    tseam(i)=filta*tsea (i) + epsflt*(tseam(i)+xsea(i))
                    qsfc (i)=MAX(1.0e-12_r8,qsfc(i))
                    tsfcm(i,latitu)=filta*tsfc0 (i,latitu) + epsflt*(tsfcm(i,latitu)+tsfc(i))
                    qsfcm(i,latitu)=filta*qsfc0 (i,latitu) + epsflt*(qsfcm(i,latitu)+qsfc(i))
                 END IF
                 tsea (i) = xsea(i)
                 qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
                 tsfc0(i,latitu) = tsfc(i)
                 qsfc0(i,latitu) = qsfc(i)
              ELSE
                 tsea (i) = xsea(i)
                 tseam(i) = xsea(i)
                 qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
                 tsfc0(i,latitu) = tsfc(i)
                 qsfc0(i,latitu) = qsfc(i)
                 tsfcm(i,latitu) = tsfc(i)
                 qsfcm(i,latitu) = qsfc(i)
              END IF
       END IF
       IF(MskAnt(i) == 1_i8 .and. tsea(i).LT.0.0e0_r8.AND.tsurf(i).GE.tice+0.01e0_r8) THEN
              tseam(i) = tsea (i)
              tsfcm(i,latitu) = tsfc0(i,latitu)
              qsfcm(i,latitu) = qsfc0(i,latitu)
       END IF
    END DO
    
    ncount=0
    DO i=1,nCols
       IF(imask(i).GE.1_i8) THEN
          ncount=ncount+1
          bstar (i)       = bstar1(ncount)
          sens  (i)       = (hc   (ncount) + hg(ncount))*(1.0_r8/dtc3x)
          evap  (i)       = (ec   (ncount) + eg(ncount))*(1.0_r8/dtc3x)
          QSfc0(i,latitu)=MAX(1.0e-12_r8,qm0 (ncount,latitu))
          QSfcm(i,latitu)=MAX(1.0e-12_r8,qmm (ncount,latitu))
          TSfc0(i,latitu)=tm0 (ncount,latitu)
          TSfcm(i,latitu)=tmm (ncount,latitu)
          zlwup_SSiB(ncount,latitu)  =  zlwup_local(ncount)  
          z0x          (ncount) =  z0    (ncount) 
       END IF
    END DO       
!    DO j=1,icg
!       ncount=0
!       DO i=1,nCols
!          IF(imask(i).GE.1_i8) THEN
!             ncount=ncount+1
!             snow_gbl     (ncount,latitu,j) =snow       (ncount,j)
!          END IF
!       END DO
!    END DO

  END SUBROUTINE SSiB_Driver

  SUBROUTINE CopySurfaceData(itype,mon,colrad2,xday,idatec,nsx,nCols,nmax,latitu)
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: nmax
    INTEGER      , INTENT(IN   ) :: latitu
    INTEGER      , INTENT(in )   :: itype   (nCols)
    INTEGER      , INTENT(in )   :: mon     (nCols)
    REAL(KIND=r8), INTENT(in )   :: colrad2 (nCols)
    REAL(KIND=r8), INTENT(in )   :: xday
    INTEGER      , INTENT(in )   :: idatec(4)
    INTEGER      , INTENT(inout) :: nsx (nCols)
    INTEGER :: i,j
    xcover = xcover_fixed
    zlt    = zlt_fixed
    green  = green_fixed
    ph2    = ph2_fixed
    ph1    = ph1_fixed
    defac  = defac_fixed
    tu     = tu_fixed
    tll    = tll_fixed
    topt   = topt_fixed
    rstpar = rstpar_fixed
    chil   = chil_fixed
    DO j=1,icg
       DO i=1,nmax
          vcover_gbl (i,latitu,j) =  xcover_fixed(itype(i),mon(i),j)
          zlt_gbl    (i,latitu,j) =  zlt_fixed   (itype(i),mon(i),j)
          green_gbl  (i,latitu,j) =  green_fixed (itype(i),mon(i),j)
          chil_gbl   (i,latitu,j) =  chil_fixed  (itype(i),j)
          topt_gbl   (i,latitu,j) =  topt_fixed  (itype(i),j)
          tll_gbl    (i,latitu,j) =  tll_fixed   (itype(i),j)
          tu_gbl     (i,latitu,j) =  tu_fixed    (itype(i),j)
          defac_gbl  (i,latitu,j) =  defac_fixed (itype(i),j)
          ph1_gbl    (i,latitu,j) =  ph1_fixed   (itype(i),j)
          ph2_gbl    (i,latitu,j) =  ph2_fixed   (itype(i),j)
          rstpar_gbl (i,latitu,j,1)= rstpar_fixed(itype(i),j,1)
          rstpar_gbl (i,latitu,j,2)= rstpar_fixed(itype(i),j,2)
          rstpar_gbl (i,latitu,j,3)= rstpar_fixed(itype(i),j,3)
       END DO
    END DO
    CALL wheat (latitu,itype ,nmax  ,colrad2 ,mon ,xday   ,yrl   , &
         idatec,monl  ,nsx    )

  END SUBROUTINE CopySurfaceData
  ! airmod :alteration of aerodynamic transfer properties in case of snow
  !         accumulation.
  !



  SUBROUTINE airmod (tg, capac, z0x, d, rdc, rbc, itype, &
       mon, nmax, ncols)
    !
    !
    !-----------------------------------------------------------------------
    !       input parameters
    !-----------------------------------------------------------------------
    !   tg............ground temperature
    !   tf............freezing point
    !   z2............height of canopy top
    !   capac(cg).....liquid water stored on canopy/ground cover foliage
    !                                                            (m)
    !   d.............displacement height                        (m)
    !   z0x...........roughness length                           (m)
    !   rdc...........constant related to aerodynamic resistance
    !                 between ground and canopy air space
    !   rbc...........constant related to bulk boundary layer
    !                 resistance
    !-----------------------------------------------------------------------
    !      output parameters
    !-----------------------------------------------------------------------
    !   d.............displacement height                        (m)
    !   z0x...........roughness length                           (m)
    !   rdc...........constant related to aerodynamic resistance
    !                 between ground and canopy air space
    !   rbc...........constant related to bulk boundary layer
    !                 resistance
    !-----------------------------------------------------------------------
    !=======================================================================
    !   ncols.........Numero de ponto por faixa de latitude
    !   ityp..........Numero do tipo de solo      13
    !   imon..........Numero maximo de meses no ano (12)
    !   mon...........Numero do mes do ano (1-12)
    !   nmax
    !   xd............Deslocamento do plano zero (m)
    !   itype.........Classe de textura do solo
    !=======================================================================
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: mon(ncols)
    INTEGER, INTENT(in   ) :: nmax
    !
    !     vegetation and soil parameters
    !
    INTEGER, INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0x   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: d     (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rdc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rbc   (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(in   ) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: capac(ncols,2)
    !
    REAL(KIND=r8)    :: sdep(ncols)
    REAL(KIND=r8)    :: xz  (ncols)
    !
    INTEGER :: i
    INTEGER :: ntyp

    DO i = 1, nmax
       IF( (tg(i) <= tf) .AND. (capac(i,2) >= 0.001_r8) )THEN
          ntyp=itype(i)
          xz  (i)=z2(ntyp,mon(i))
          sdep(i)=capac(i,2)*5.0_r8
          sdep(i)=MIN( sdep(i) , xz(i)*0.95_r8 )
          d  (i)=xz (i)-( xz(i)- d(i) )/xz(i)*(xz(i)-sdep(i))
          z0x(i)=z0x(i)/( xz(i)-xd(ntyp,mon(i)))*(xz(i)-d   (i))
          rdc(i)=rdc(i)*( xz(i)-sdep(i) )/xz(i)
          rbc(i)=rbc(i)*xz(i)/( xz(i)-sdep(i) )
       END IF
    END DO
  END SUBROUTINE airmod





  SUBROUTINE temres(&
       bps   ,psb   ,em    ,gmt   ,gmq   ,psy   ,rcp   ,wc    ,wg    , &
       fc    ,fg    ,hr    ,hgdtg ,hgdtc ,hgdtm ,hcdtg ,hcdtc ,hcdtm , &
       egdtg ,egdtc ,egdqm ,ecdtg ,ecdtc ,ecdqm ,deadtg,deadtc,deadqm, &
       ect   ,eci   ,egt   ,egi   ,egs   ,ec    ,eg    ,hc    ,hg    , &
       ecidif,egidif,ra    ,rb    ,rd    ,rc    ,rg    ,ta    ,ea    , &
       etc   ,etg   ,btc   ,btg   ,radt  ,rst   ,rsoil ,hrr   ,cc    , &
       cg    ,satcap,dtc   ,dtg   ,dtm   ,dqm   ,thermk,tc    ,tg    , &
       td    ,capac ,qm    ,tm    ,psur  ,dtc3x , &
       nmax  ,vcover,ncols )
    !
    !-----------------------------------------------------------------------
    ! temres :performs temperature tendency equations with interception loss.
    !-----------------------------------------------------------------------
    !     ncols.......Numero de ponto por faixa de latitude
    !     ityp........numero das classes de solo 13
    !     imon........Numero maximo de meses no ano (12)
    !     icg.........Parametros da vegetacao (icg=1 topo e icg=2 base)
    !     pie.........Constante Pi=3.1415926e0
    !     stefan .....Constante de Stefan Boltzmann
    !     cp..........specific heat of air (j/kg/k)
    !     hl..........heat of evaporation of water   (j/kg)
    !     grav........gravity constant      (m/s**2)
    !     tf..........Temperatura de congelamento (K)
    !     epsfac......Constante 0.622 Razao entre as massas moleculares do vapor
    !                 de agua e do ar seco
    !     dtc3x.......time increment dt
    !     nmax........
    !     xcover......Fracao de cobertura vegetal icg=1 topo
    !     xcover......Fracao de cobertura vegetal icg=2 base
    !     vcover......Fracao de cobertura vegetal icg=1 topo
    !     vcover......Fracao de cobertura vegetal icg=2 topo
    !     qm..........specific humidity of reference (fourier)
    !     tm..........Temperature of reference (fourier)
    !     psur........surface pressure in mb
    !     tc..........Temperatura da copa "dossel" canopy leaf temperature(K)
    !     tg..........Temperatura da superficie do solo ground temperature (K)
    !     td .........Temperatura do solo profundo (K)
    !     capac(iv)...Agua interceptada iv=1 no dossel "water store capacity of leaves"(m)
    !     capac(iv)...Agua interceptada iv=2 na cobertura do solo (m)
    !     ra..........Resistencia Aerodinamica (s/m)
    !     rb..........bulk boundary layer resistance             (s/m)
    !     rd..........aerodynamic resistance between ground
    !                 and canopy air space                       (s/m)
    !     rc..........Resistencia do topo da copa (s/m)
    !     rg..........Resistencia da base da copa (s/m)
    !     ta..........Temperatura no nivel de fonte de calor do dossel (K)
    !     ea..........Pressao de vapor
    !     etc.........Pressure of vapor at top of the copa
    !     etg.........Pressao de vapor no base da copa
    !     btc.........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i)).
    !     btg.........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))
    !     radt........net heat received by canopy/ground vegetation
    !     rst.........Resisttencia Estomatica "Stomatal resistence" (s/m)
    !     rsoil ......Resistencia do solo (s/m)
    !     hrr.........rel. humidity in top layer
    !     cc..........heat capacity of the canopy
    !     cg..........heat capacity of the ground
    !     satcap......saturation liquid water capacity         (m)
    !     dtc.........dtc(i)=pblsib(i,2,5)*dtc3x
    !     dtg.........dtg(i)=pblsib(i,1,5)*dtc3x
    !     dtm.........dtm(i)=pblsib(i,3,5)*dtc3x
    !     dqm.........dqm(i)=pblsib(i,4,5)*dtc3x
    !     thermk......canopy emissivity
    !     ect.........Transpiracao(J/m*m)
    !     eci.........Evaporacao da agua interceptada (J/m*m)
    !     egt.........Transpiracao na base da copa (J/m*m)
    !     egi.........Evaporacao da neve (J/m*m)
    !     egs.........Evaporacao do solo arido (J/m*m)
    !     ec..........Soma da Transpiracao e Evaporacao da agua interceptada pelo
    !                 topo da copa   ec   (i)=eci(i)+ect(i)
    !     eg..........Soma da transpiracao na base da copa +  Evaporacao do solo arido
    !                 +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"
    !     hc..........total sensible heat lost of top from the veggies.
    !     hg..........total sensible heat lost of base from the veggies.
    !     ecidif......check if interception loss term has exceeded canopy storage
    !                 ecidif(i)=MAX(0.0   , eci(i)-capac(i,1)*hl3 )
    !     egidif......check if interception loss term has exceeded canopy storage
    !                 ecidif(i)=MAX(0.0   , egi(i)-capac(i,1)*hl3 )
    !     hgdtg ......n.b. fluxes expressed in joules m-2
    !     hgdtc.......n.b. fluxes expressed in joules m-2
    !     hgdtm.......n.b. fluxes expressed in joules m-2
    !     hcdtg.......n.b. fluxes expressed in joules m-2
    !     hcdtc.......n.b. fluxes expressed in joules m-2
    !     hcdtm.......n.b. fluxes expressed in joules m-2
    !     egdtg.......partial derivative calculation for latent heat
    !     egdtc.......partial derivative calculation for latent heat
    !     egdqm.......partial derivative calculation for latent heat
    !     ecdtg ......partial derivative calculation for latent heat
    !     ecdtc ......partial derivative calculation for latent heat
    !     ecdqm.......partial derivative calculation for latent heat
    !     deadtg......
    !     deadtc......
    !     deadqm......
    !     bps.........
    !     psb.........
    !     em..........Pressao de vapor da agua
    !     gmt.........
    !     gmq.........specific humidity of reference (fourier)
    !     psy.........(cp/(hl*epsfac))*psur(i)
    !     rcp.........densidade do ar vezes o calor especifico do ar
    !     wc..........Minimo entre 1 e a razao entre a agua interceptada pelo
    !                 indice de area foliar no topo da copa
    !     wg..........Minimo entre 1 e a razao entre a agua interceptada pelo
    !                 indice de area foliar na base da copa
    !     fc..........Condicao de oravalho 0 ou 1 na topo da copa
    !     fg..........Condicao de oravalho 0 ou 1 na base da copa
    !     hr..........rel. humidity in top layer
    !-----------------------------------------------------------------------

    INTEGER, INTENT(in   ) :: ncols

    REAL(KIND=r8),    INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: nmax
    !
    !     vegetation and soil parameters
    !
    REAL(KIND=r8),    INTENT(in) :: vcover(ncols,icg)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8), INTENT(in   ) :: qm  (ncols)
    REAL(KIND=r8), INTENT(in   ) :: tm  (ncols)
    REAL(KIND=r8), INTENT(in   ) :: psur(ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8), INTENT(in   ) :: tc   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: tg   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: td   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: capac(ncols,2)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8), INTENT(in   ) :: ra    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: rb    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: rd    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: rc    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: rg    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ta    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ea    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: etc   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: etg   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: btc   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: btg   (ncols)
    REAL(KIND=r8), INTENT(inout) :: radt  (ncols,icg)
    REAL(KIND=r8), INTENT(inout) :: rst   (ncols,icg)
    REAL(KIND=r8), INTENT(in   ) :: rsoil (ncols)
    REAL(KIND=r8), INTENT(in   ) :: hrr   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: cc    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: cg    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: satcap(ncols,icg)
    REAL(KIND=r8), INTENT(inout  ) :: dtc   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: dtg   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: dtm   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: dqm   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: thermk(ncols)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8), INTENT(inout  ) :: ect   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: eci   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egt   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egi   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egs   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ec    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: eg    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hc    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hg    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ecidif(ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egidif(ncols)
    !
    !     derivatives
    !
    REAL(KIND=r8), INTENT(inout  ) :: hgdtg (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hgdtc (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hgdtm (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hcdtg (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hcdtc (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hcdtm (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egdtg (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egdtc (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: egdqm (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ecdtg (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ecdtc (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: ecdqm (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: deadtg(ncols)
    REAL(KIND=r8), INTENT(inout  ) :: deadtc(ncols)
    REAL(KIND=r8), INTENT(inout  ) :: deadqm(ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8), INTENT(in   ) :: bps   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: psb   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: em    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: gmt   (ncols,3)
    REAL(KIND=r8), INTENT(in   ) :: gmq   (ncols,3)
    REAL(KIND=r8), INTENT(in   ) :: psy   (ncols)
    REAL(KIND=r8), INTENT(in   ) :: rcp   (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: wc    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: wg    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: fc    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: fg    (ncols)
    REAL(KIND=r8), INTENT(inout  ) :: hr    (ncols)


    REAL(KIND=r8) :: vcover2(ncols,icg)

    REAL(KIND=r8)    :: pblsib(ncols,4,5)
    REAL(KIND=r8)    :: coc
    REAL(KIND=r8)    :: rsurf
    REAL(KIND=r8)    :: cog1
    REAL(KIND=r8)    :: cog2
    REAL(KIND=r8)    :: d1
    REAL(KIND=r8)    :: d2
    REAL(KIND=r8)    :: d1i
    REAL(KIND=r8)    :: top
    REAL(KIND=r8)    :: ak    (ncols)
    REAL(KIND=r8)    :: ah    (ncols)
    REAL(KIND=r8)    :: cci   (ncols)
    REAL(KIND=r8)    :: cgi   (ncols)
    REAL(KIND=r8)    :: ecpot (ncols)
    REAL(KIND=r8)    :: egpot (ncols)
    REAL(KIND=r8)    :: ecf
    REAL(KIND=r8)    :: egf
    REAL(KIND=r8)    :: coct
    REAL(KIND=r8)    :: cogt
    REAL(KIND=r8)    :: cogs1
    REAL(KIND=r8)    :: cogs2
    REAL(KIND=r8)    :: psyi  (ncols)
    REAL(KIND=r8)    :: fac1
    REAL(KIND=r8)    :: rcdtc (ncols)
    REAL(KIND=r8)    :: rcdtg (ncols)
    REAL(KIND=r8)    :: rgdtc (ncols)
    REAL(KIND=r8)    :: rgdtg (ncols)
    REAL(KIND=r8), PARAMETER :: capi  =1.0_r8/4.0e-3_r8
    REAL(KIND=r8)    :: timcon
    REAL(KIND=r8)    :: timcn2
    REAL(KIND=r8)    :: tim
    REAL(KIND=r8)    :: dtc3xi
    REAL(KIND=r8)    :: fak
    REAL(KIND=r8)    :: fah
    INTEGER :: i
    REAL(KIND=r8)    :: stb4
    REAL(KIND=r8)    :: stb8
    REAL(KIND=r8)    :: hlat3

    timcon = pie/86400.0_r8
    timcn2 = 2.0_r8   * timcon
    tim    = 1.0_r8   + timcn2*dtc3x
    dtc3xi = 1.0_r8   / dtc3x
    fak    = 0.01_r8  * grav/cp
    fah    = 0.01_r8  * grav/hl
    vcover2=vcover

    DO i = 1, nmax
       !
       !                  --                                                   --
       !                 |      razao entre a agua interceptada no topo da copa  |
       !     wc = Minimo*| 1 , --------------------------------------------------|
       !                 |     indice de area foliar no topo da copa             |
       !                  --                                                   --
       !
       wc  (i)=MIN( 1.0_r8   , capac(i,1)/satcap(i,1))
       !                  --                                                    --
       !                 |       razao entre a agua interceptada na base da copa  |
       !     wg = Minimo*| 1  ,---------------------------------------------------|
       !                 |      indice de area foliar na base da copa             |
       !                  --                                                    --
       !
       wg  (i)=MIN( 1.0_r8   , capac(i,2)/satcap(i,2))
       !
       ! Temperatura de congelamento (K)
       !
       IF (tg(i) <= tf) THEN
          vcover2(i,2)=1.0_r8
          wg    (i)  =MIN(1.0_r8   ,capac(i,2)*capi)
          !
          !     rsoil ......Resistencia do solo (s/m)
          !
          rst   (i,2)=rsoil(i)
       END IF
       !
       !
       !  DT     d    d[w'T']
       ! ---- = --   ----------
       !  Dt     dt     dz
       !
       !H =rho*cp*w'T'
       !
       !          H
       !w'T' = -------
       !        rho*cp
       !
       !  DT     d     1          dH
       ! ---- = ---- --------- * -----
       !  Dt     dt    rho*cp      dz
       !
       ! P =rho*R*T and P = rho*g*Z
       !
       !                           P
       ! DP = rho*g*DZ and rho = ----
       !                          R*T
       !
       ! 1             1
       !----= rho*g* -----
       ! DZ            DP
       !
       ! 1     1         1
       !--- * ---- = g* -----
       !rho    DZ        DP
       !
       !
       !  DT     d      g      dH
       ! ---- = ---- *----- * -----
       !  Dt     dt     cp     dP
       !
       !           g      d
       ! ak(i) = ----- * -----
       !           cp     dP

       !                grav            1
       !ak(i) = 0.01 * ------ * ------------------
       !                 cp      (psb(i)*bps(i))
       !
       !                                        -(R/Cp)
       !                  g               sl(k)
       ! ak(i) = 0.01 * ------- * -------------------------------
       !                  cp        (P * (si(k) - si(k+1)))
       !
       ak  (i) =fak/(psb(i)*bps(i))
       !
       !L =rho*hl*w'Q'
       !
       !
       !                 g                1
       !ah(i) = 0.01 * ------ * --------------------------
       !                 hl       (P * (si(k) - si(k+1)))
       !
       ah  (i) =fah/ psb(i)
       !
       !     cc..........heat capacity of the canopy
       !     cg..........heat capacity of the ground
       !
       cgi (i) =1.0_r8   / cg(i)
       cci (i) =1.0_r8   / cc(i)
       !
       ! rcp ---- densidade do ar vezes o calor especifico do ar
       !
       !(cp/(hl*epsfac))*psur(i)
       !
       psyi(i) =rcp(i)/psy(i)
    END DO
    !
    !     partial derivative calculations for sensible heat
    !
    DO i = 1, nmax
       !
       !           1          1          1
       ! d1     =------- + -------- + --------
       !          ra(i)      rb(i)      rd(i)
       !
       !            rb(i)*rd(i) + ra(i)*rd(i) + ra(i)*rb(i)
       ! d1     = --------------------------------------------
       !                     ra(i)*rb(i)*rd(i)
       !
       d1     =1.0_r8/ra(i) + 1.0_r8/rb(i) + 1.0_r8/rd(i)
       !
       !          rcp(i)     rcp(i)     rcp(i)
       ! d1i =  --------- + -------- + --------
       !          ra(i)      rb(i)      rd(i)
       !
       !
       !          rb(i)*rd(i)*rcp(i)  +   ra(i)*rd(i)*rcp(i)  +   ra(i)*rb(i)*rcp(i)
       ! d1i =  ----------------------------------------------------------------------
       !                              ra(i)*rb(i)*rd(i)
       !
       d1i    =rcp(i)/d1
       !
       !       --                                --
       !      | tg(i)     tc(i)      tm(i)*bps(i)  |  /
       !ta(i)=|------- + -------- + -------------- | /d1
       !      | rd(i)      rb(i)         ra(i)     |/
       !       --                                --
       !
       ta(i)=( tg(i)/rd(i) + tc(i)/rb(i) + tm(i)*bps(i)/ra(i) )/d1
       !
       !dtc3x = time increment dt
       !rcp----densidade do ar vezes o calor especifico do ar
       !
       !
       !   total sensible heat lost of top from the veggies.
       !                  (tc(i)-ta(i))
       !hc(i) = rcp(i) * ----------------*dt
       !                      rb(i)
       !
       hc(i)=rcp(i) * ( tc(i) - ta(i) ) / rb(i) * dtc3x
       !
       !   total sensible heat lost of base from the veggies.
       !
       !                  (tg(i)-ta(i))
       !hg(i) = rcp(i) * ---------------*dt
       !                      rd(i)
       !
       hg(i)=rcp(i) * ( tg(i) - ta(i) ) / rd(i) * dtc3x
       !                                              J
       !     n.b. fluxes expressed in joules m-2  = ------
       !                                             m^2
       !          rcp(i)     rcp(i)     rcp(i)
       ! d1i =  --------- + -------- + --------
       !          ra(i)      rb(i)      rd(i)
       !
       !          rb(i)*rd(i)*rcp(i)  +   ra(i)*rd(i)*rcp(i)  +   ra(i)*rb(i)*rcp(i)
       ! d1i =  ----------------------------------------------------------------------
       !                              ra(i)*rb(i)*rd(i)
       !
       !                       --               --
       !              d1i     |  1.0       1.0    |
       !hcdtc(i) =  ------- * | ------ + -------- |
       !             rb(i)    |  ra(i)     rd(i)  |
       !                       --               --
       !
       hcdtc(i)= d1i   / rb(i)*( 1.0_r8/ra(i) + 1.0_r8/rd(i) )
       !
       !                -d1i
       !hcdtg(i) =  ---------------
       !             rb(i) * rd(i)
       !
       hcdtg(i)=-d1i   / ( rb(i)*rd(i) )
       !
       !                     ra(i)*rb(i)*rd(i)*rcp(i)
       ! d1i     = --------------------------------------------
       !              rb(i)*rd(i) + ra(i)*rd(i) + ra(i)*rb(i)
       !
       !               -d1i
       ! hcdtm(i)= ----------------- * bps(i)
       !            ( rb(i)*ra(i) )
       !
       !                          - rd(i)*rcp(i)
       ! hcdtm(i)= ---------------------------------------------- * bps(i)
       !              rb(i)*rd(i) + ra(i)*rd(i) + ra(i)*rb(i)
       !

       hcdtm(i)=-d1i   / ( rb(i)*ra(i) ) *bps(i)
       !
       !                        --               --
       !               d1i     |  1.0       1.0    |
       ! hgdtg(i) =  ------- * | ------ + -------- |
       !              rd(i)    |  ra(i)     rb(i)  |
       !                        --               --
       !
       hgdtg(i)= d1i   / rd(i)*( 1.0_r8/ra(i) + 1.0_r8/rb(i))
       !
       !                -d1i
       !hgdtc(i) = -----------------
       !            ( rd(i)*rb(i) )
       !
       !
       hgdtc(i)=-d1i   / ( rd(i)*rb(i) )
       !
       !                    -d1i                (R/Cp)
       ! hgdtm(i) = ----------------- *   sl(k)
       !             ( rd(i)*ra(i) )
       !
       !
       !                          - rb(i)*rcp(i)                          (R/Cp)
       ! hgdtm(i)= ---------------------------------------------- * sl(k)
       !              rb(i)*rd(i) + ra(i)*rd(i) + ra(i)*rb(i)

       hgdtm(i)=-d1i   / ( rd(i)*ra(i) ) *bps(i)
       !
    END DO
    !
    !     partial derivative calculations for longwave radiation flux
    !
    stb4  = 4.0_r8 * stefan
    stb8  = 8.0_r8 * stefan
    !
    DO i = 1, nmax
       fac1     = vcover2(i,1)*(1.0_r8  - thermk(i))
       rcdtc(i) = fac1 * stb8 * tc(i)*tc(i)*tc(i)
       rcdtg(i) =-fac1 * stb4 * tg(i)*tg(i)*tg(i)
       rgdtc(i) =-fac1 * stb4 * tc(i)*tc(i)*tc(i)
       rgdtg(i) =        stb4 * tg(i)*tg(i)*tg(i)
    END DO
    DO i = 1, nmax
       !
       !     partial derivative calculation for latent heat
       !     modification for soil dryness : hr=rel. humidity in top layer
       !
       hr   (i)  = hrr(i)   * fg(i) + 1.0_r8 - fg(i)
       !
       !     fc = Condicao de oravalho 0 ou 1 na topo da copa
       !
       rc   (i)  = rst(i,1) * fc(i) + 2.0_r8 * rb(i)
       !
       !        ( 1.0_r8 - wc(i) )         wc(i)
       ! coc = -------------------- + ------------------
       !              rc(i)            (2.0_r8 * rb(i))
       !
       coc       = ( 1.0_r8 - wc(i) ) / rc(i) + wc(i)/(2.0_r8 * rb(i))
       !
       ! fg = Condicao de oravalho 0 ou 1 na base da copa
       !
       rg   (i)  = rst(i,2)*fg(i)
       !
       rsurf     = rsoil(i)*fg(i)
       !
       !     hr..........rel. humidity in top layer
       !     vcover......Fracao de cobertura vegetal icg=1 topo
       !     vcover......Fracao de cobertura vegetal icg=2 topo
       !
       !                      (1 - wg(i))            (1 - vcover(i,2))              vcover(i,2)
       ! cog1 = vcover(i,2)*--------------- + hr(i)*------------------- + hr(i)*----------------------
       !                     (rg(i)+rd(i))            (rsurf + rd(i))            (rsurf + rd(i) + 44)
       !
       cog1      =   vcover2(i,2)*(1.0_r8 - wg(i))/(rg(i)+rd(i)) &
            + (1.0_r8 - vcover2(i,2))/(rsurf + rd(i)) * hr(i) &
            + vcover2(i,2) / (rsurf + rd(i) + 44.0_r8) * hr(i)
       !
       !                      (1 - wg(i))      (1 - vcover(i,2))        vcover(i,2)
       ! cog2 = vcover(i,2)*--------------- + ------------------- + ----------------------
       !                     (rg(i)+rd(i))      (rsurf + rd(i))      (rsurf + rd(i) + 44)
       !
       cog2      = vcover2(i,2)*(1.0_r8 - wg(i))/(rg(i)+rd(i)) &
            +     (1.0_r8 - vcover2(i,2))/(rsurf + rd(i)) &
            +      vcover2(i,2)/(rsurf   +rd(i)+44.0_r8)

       !                       (1 - wg(i))      hr(i)*(1 - vcover(i,2))      hr(i)*vcover(i,2)      wg(i)* vcover(i,2)
       ! cog1 = vcover(i,2) * -------------- + ------------------------- + --------------------- + --------------------
       !                      (rg(i)+rd(i))     (rsurf + rd(i))            (rsurf + rd(i) + 44)           rd(i)
       !
       !
       cog1      = cog1 + wg(i) / rd(i)*vcover2(i,2)
       !
       !                      wg(i)
       ! cog2      = cog2 + -------- * vcover(i,2)
       !                      rd(i)
       !
       cog2      = cog2 + wg(i)/rd(i)*vcover2(i,2)
       !
       !        1.0       ( 1.0_r8 - wc(i) )         wc(i)
       !d2 = --------- + -------------------- + ------------------ + cog2
       !       ra(i)            rc(i)            (2.0_r8 * rb(i))
       !
       d2        = 1.0_r8/ra(i) + coc + cog2
       !
       !                                      em(i)
       !top = coc * etc(i) + cog1 * etg(i) + -------
       !                                      ra(i)
       !
       top       = coc * etc(i) + cog1 * etg(i) + em(i)/ra(i)
       !
       ea (i)    = top/d2
       !
       !       psyi(i) =rcp(i)/psy(i)
       !
       ! The rate of evaporation from the wetted portions of the vegetation
       !
       !        ( 1 - wc(i) )         wc(i)
       ! coc = ---------------- + --------------
       !            rc(i)          ( 2 * rb(i) )
       !
       !The latent heat fluxes from the canopy is defined by:
       !
       !
       !            --        --                      --             --
       !           |            |     rho(i) * cp    |  wc     1 - wc  |
       !ec = LEc = | e[Tc] - ea | * -------------- * | ---- + ---------|
       !           |            |       psy(i)       |  rb     rb + rc |
       !            --        --                      --             --
       !
       ec (i)    = ( etc(i)-ea(i) )  *  coc * psyi(i) * dtc3x
       !
       !The latent heat fluxes from the ground is defined by:
       !
       !             --           --                       --        --
       !            |               |      rho(i) * cp    |     1      |
       ! eg = LEgs =|fh*e[Tgs] - ea | *  -------------- * |------------|
       !            |               |        psy(i)       | rsurf + rd |
       !             --           --                       --        --
       !
       eg (i)    = (etg(i)*cog1   - ea(i)*cog2   )*psyi(i)*dtc3x
       !
       deadtc(i) = btc(i) *  coc / d2
       !
       deadtg(i) = btg(i) * cog1 / d2
       !
       !                             psur(i)
       ! deadqm(i) = epsfac * ------------------------------------------
       !                       ( ( epsfac + qm(i) )**2  *  ra(i)*d2    )
       !
       deadqm(i) = epsfac * psur(i)/( (epsfac+qm(i))**2 * ra(i)*d2    )
       !
       ecdtc(i)  = (btc(i) - deadtc(i) ) * coc * psyi(i)
       !
       ecdtg(i)  = -deadtg(i) * coc   * psyi(i)
       !
       ecdqm(i)  = -deadqm(i) * coc   * psyi(i)
       !
       egdtg(i)  = ( btg(i) * cog1 - deadtg(i) * cog2 )*psyi(i)
       !
       egdtc(i)  = -deadtc(i) * cog2   * psyi(i)
       !
       egdqm(i)  = -deadqm(i) * cog2   * psyi(i)
       !
    END DO
    !
    !     solve for time changes of pbl and sib variables,
    !     using a semi-implicit scheme.
    !
    DO i = 1, nmax
       !
       !     tg equation
       !
       !     cc..........heat capacity of the canopy
       !     cg..........heat capacity of the ground
       !
       !            1.0                  1.0
       ! cgi (i) = ----- = ---------------------------------
       !           cg(i)     heat capacity of the ground
       !
       !           1.0                         1.0
       ! cci (i) =------   =  --------------------------------------
       !           cc(i)          heat capacity of the canopy
       !
       !             2 * pi * dt     s
       !tim = 1.0 + ------------- = ---
       !               86400.0       s
       !
       pblsib(i,1,1) = tim + dtc3x * cgi(i) * (hgdtg(i) + egdtg(i) + rgdtg(i))
       pblsib(i,1,2) =       dtc3x * cgi(i) * (hgdtc(i) + egdtc(i) + rgdtc(i))
       pblsib(i,1,3) =       dtc3x * cgi(i) * hgdtm(i)
       pblsib(i,1,4) =       dtc3x * cgi(i) * egdqm(i)
       !
       !     tc equation
       !
       pblsib(i,2,1) =          dtc3x * cci(i) * ( hcdtg(i) + ecdtg(i) + rcdtg(i) )
       !
       pblsib(i,2,2) = 1.0_r8 + dtc3x * cci(i) * ( hcdtc(i) + ecdtc(i) + rcdtc(i) )
       !
       !
       !               -d1i                  (R/Cp)
       ! hcdtm(i)= ----------------- *  sl(k)
       !            Cc(i)*( rb(i)*ra(i) )
       !
       pblsib(i,2,3) = dtc3x * cci(i) * hcdtm(i)
       !
       pblsib(i,2,4) = dtc3x * cci(i) * ecdqm(i)
       !
       !     tm equation
       !
       !                                        -(R/Cp)
       !                  g               sl(k)
       ! ak(i) = 0.01 * ------- * -------------------------------
       !                  cp        (P * (si(k) - si(k+1)))
       !
       !
       pblsib(i,3,1) = -dtc3x * ak(i) * ( hgdtg(i) + hcdtg(i) )
       !
       pblsib(i,3,2) = -dtc3x * ak(i) * ( hgdtc(i) + hcdtc(i) )
       !
       !
       !      --   --  -(R/Cp)    --   --  -(R/Cp)
       !     |  P    |           |       |
       !bps  |-------|        == |sl(k)  |
       !     |  P0   |           |       |
       !      --   --             --   --
       !
       !             --   --  -(R/Cp)
       !            |  P    |
       !Tpot =  T * |-------|
       !            |  P0   |
       !             --   --
       !
       !
       ! P =rho*R*T and P = rho*g*Z
       !
       !                           P
       ! DP = rho*g*DZ and rho = ----
       !                          R*T
       !
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
       !                                        -(R/Cp)
       !                  g               sl(k)
       ! ak(i) = 0.01 * ------- * -------------------------------
       !                  cp        (P * (si(k) - si(k+1)))
       !
       !
       ! pblsib(i,3,3) = gmt(i,2) - 2*dt*ak(i)*(hgdtm(i) + hcdtm(i))
       !
       !                             T      1
       !Pbl_KMbyDZ_1  =   2*Dt * Km*---- * ----
       !                             dZ     dZ
       !
       !                             T      1
       !Pbl_KMbyDZ_1  =   2*Dt * Km*---- * ----
       !                             dZ     dZ
       !
       !           --                                       --     --                                           --
       !          |                                           |   |    Pbl_KMbyDZ_1(i,k)*Pbl_KMbyDZ_2(i,k+1)      |
       !gmt(i,2) =|1.0 + Pbl_KMbyDZ_1(i,k) + Pbl_KMbyDZ_2(i,k)| - |-----------------------------------------------|
       !          |                                           |   |1.0 + Pbl_KMbyDZ_1(i,k+1) + Pbl_KMbyDZ_2(i,k+1)|
       !           --                                       --     --                                           --

       pblsib(i,3,3) = gmt(i,2) - dtc3x * ak(i) * ( hgdtm(i) + hcdtm(i) )
       pblsib(i,3,4) = 0.0_r8
       !
       !     qm equation
       !
       !
       !                 g                1
       !ah(i) = 0.01 * ------ * --------------------------
       !                 hl       (P * (si(k) - si(k+1)))
       !
       pblsib(i,4,1) = - dtc3x * ah(i) * ( egdtg(i) + ecdtg(i) )
       pblsib(i,4,2) = - dtc3x * ah(i) * ( egdtc(i) + ecdtc(i) )
       pblsib(i,4,3) =   0.0_r8
       !
       pblsib(i,4,4) =   gmq(i,2) - dtc3x * ah(i) * ( egdqm(i) + ecdqm(i) )
       !
       !                                                           Rngs
       ! radt = net heat received by canopy/ground vegetation  = --------
       !                                                            dt
       !
       !      dTgs                          2*PI*Cgs
       ! Cgs*------ = Rngs - Hgs - LHgs - ------------ * (Tgs - Td)
       !      dt                            dayleg
       !
       !      dTgs     Rngs      Hgs     LHgs         2*PI
       !     ------ = ------ - ------ - ------- - -------------- * (Tgs - Td)
       !      dt       Cgs       Cgs      Cgs         dayleg
       !
       !               --                     --
       !      dTgs    |  Rngs      Hgs     LHgs |        2*PI
       !     ------ = | ------ - ------ - ------| - -------------- * (Tgs - Td)
       !      dt      |  Cgs       Cgs      Cgs |        dayleg
       !               --                     --
       !
       !               --                   --
       !      dTgs    |                       |    1           2*PI
       !     ------ = |  Rngs - ( Hgs + LHgs )| *------ - -------------- * (Tgs - Td)
       !      dt      |                       |    Cgs        dayleg
       !               --                   --
       !
       !                                                          cgi(i)         2*pi
       !pblsib(i,1,5) = (radt(i,2)* cgi(i) - ( hg(i) + eg(i) ) * -------- )  - --------- * ( tg(i) - td(i) )
       !                                                             dt          86400.0
       !            2*pi
       ! timcn2 = ---------
       !           86400.0
       pblsib(i,1,5) = (radt(i,2) - ( hg(i) + eg(i) ) * dtc3xi ) * cgi(i) - timcn2 * ( tg(i) - td(i) )
       !
       !      dTc
       ! Cc*------ = Rnc - Hc - LHc
       !      dt
       !
       !             --             --
       !      dTc   |                 |     1
       !    ------ =|Rnc - Hc - LHc   | * -----
       !      dt    |                 |     Cc
       !             --             --
       !                 --                                     --
       !                |                                   1     |     1
       !pblsib(i,2,5) = |radt(i,1) - ( hc(i) + ec(i) ) * -------- | * -------
       !                |                                   dt    |    cc(i)
       !                 --                                     --
       !                 --                                     --
       !                |                                   1     |
       !pblsib(i,2,5) = |radt(i,1) - ( hc(i) + ec(i) ) * -------- | * cci(i)
       !                |                                   dt    |
       !                 --                                     --
       pblsib(i,2,5) = (radt(i,1) - ( hc(i) + ec(i) ) * dtc3xi ) * cci(i)
       !
       !                                         -(R/Cp)              --           --
       !   dTm              g               sl(k)                    | hg(i) + hc(i) |
       ! ------ =  0.01 * ------- * ------------------------------- *| --------------|
       !   dt               cp        (P * (si(k) - si(k+1)))        |       dt      |
       !                                                              --           --
       !                                                                  --     --
       !   dTm              m*kg*K                    1                  |   J     |
       ! ------ =  0.01 * ----------- * ------------------------------- *| ------  |
       !   dt              s^2*J                      Pa                 | m^2*s   |
       !                                                                  --     --
       !   dTm              m*kg*K                  m^2                  |   N*m   |
       ! ------ =  0.01 * ----------- * ------------------------------- *| --------|
       !   dt              s^2*N*m                   N                   | m^2*s   |

       !   dTm              m*Kg*K*s^2         m^3
       ! ------ =  0.01 * ------------- * ---------------
       !   dt              s^2*kg*m*m          m^2*s

       !   dTm                K            m
       ! ------ =  0.01 * ------------- *-------
       !   dt                 m            s
       !
       !   dTm             K
       ! ------ =  0.01 * ----
       !   dt              s

       !
       !                                               -(R/Cp)              --           --
       !                          g               sl(k)                    | hg(i) + hc(i) |
       ! pblsib(i,3,5) = 0.01 * ------- * ------------------------------- *| --------------|
       !                          cp        (P * (si(k) - si(k+1)))        |       dt      |
       !                                                                    --           --
       !
       pblsib(i,3,5) =  gmt(i,3) + ak(i) * ( hg(i) + hc(i) ) * dtc3xi
       !
       !                                                             --           --
       !                         g                1                 | eg(i) + ec(i) |
       !pblsib(i,4,5) = 0.01 * ------ * ------------------------- * | --------------|
       !                         hl       (P * (si(k) - si(k+1)))   |       dt      |
       !                                                             --           --
       !
       pblsib(i,4,5) =  gmq(i,3) + ah(i) * ( eg(i) + ec(i) ) * dtc3xi
    END DO
    !
    !     solve 4 x 5 matrix equation
    !
    DO i = 1, nmax
       pblsib(i,2,2) =  pblsib(i,2,2) - pblsib(i,2,1) * ( pblsib(i,1,2) / pblsib(i,1,1) )
       pblsib(i,2,3) =  pblsib(i,2,3) - pblsib(i,2,1) * ( pblsib(i,1,3) / pblsib(i,1,1) )
       pblsib(i,2,4) =  pblsib(i,2,4) - pblsib(i,2,1) * ( pblsib(i,1,4) / pblsib(i,1,1) )
       pblsib(i,2,5) =  pblsib(i,2,5) - pblsib(i,2,1) * ( pblsib(i,1,5) / pblsib(i,1,1) )
       pblsib(i,3,2) =  pblsib(i,3,2) - pblsib(i,3,1) * ( pblsib(i,1,2) / pblsib(i,1,1) )
       pblsib(i,3,3) =  pblsib(i,3,3) - pblsib(i,3,1) * ( pblsib(i,1,3) / pblsib(i,1,1) )
       pblsib(i,3,4) =  pblsib(i,3,4) - pblsib(i,3,1) * ( pblsib(i,1,4) / pblsib(i,1,1) )
       pblsib(i,3,5) =  pblsib(i,3,5) - pblsib(i,3,1) * ( pblsib(i,1,5) / pblsib(i,1,1) )
       pblsib(i,4,2) =  pblsib(i,4,2) - pblsib(i,4,1) * ( pblsib(i,1,2) / pblsib(i,1,1) )
       pblsib(i,4,3) =  pblsib(i,4,3) - pblsib(i,4,1) * ( pblsib(i,1,3) / pblsib(i,1,1) )
       pblsib(i,4,4) =  pblsib(i,4,4) - pblsib(i,4,1) * ( pblsib(i,1,4) / pblsib(i,1,1) )
       pblsib(i,4,5) =  pblsib(i,4,5) - pblsib(i,4,1) * ( pblsib(i,1,5) / pblsib(i,1,1) )
       pblsib(i,3,3) =  pblsib(i,3,3) - pblsib(i,3,2) * ( pblsib(i,2,3) / pblsib(i,2,2) )
       pblsib(i,3,4) =  pblsib(i,3,4) - pblsib(i,3,2) * ( pblsib(i,2,4) / pblsib(i,2,2) )
       pblsib(i,3,5) =  pblsib(i,3,5) - pblsib(i,3,2) * ( pblsib(i,2,5) / pblsib(i,2,2) )
       pblsib(i,4,3) =  pblsib(i,4,3) - pblsib(i,4,2) * ( pblsib(i,2,3) / pblsib(i,2,2) )
       pblsib(i,4,4) =  pblsib(i,4,4) - pblsib(i,4,2) * ( pblsib(i,2,4) / pblsib(i,2,2) )
       pblsib(i,4,5) =  pblsib(i,4,5) - pblsib(i,4,2) * ( pblsib(i,2,5) / pblsib(i,2,2) )
       pblsib(i,4,4) =  pblsib(i,4,4) - pblsib(i,4,3) * ( pblsib(i,3,4) / pblsib(i,3,3) )
       pblsib(i,4,5) =  pblsib(i,4,5) - pblsib(i,4,3) * ( pblsib(i,3,5) / pblsib(i,3,3) )

       pblsib(i,4,5) =     pblsib(i,4,5) / pblsib(i,4,4)

       pblsib(i,3,5) =   ( pblsib(i,3,5) / pblsib(i,3,3) ) &
            - ( pblsib(i,3,4) / pblsib(i,3,3) ) * pblsib(i,4,5)

       pblsib(i,2,5) =   ( pblsib(i,2,5) / pblsib(i,2,2) ) &
            - ( pblsib(i,2,4) / pblsib(i,2,2) ) * pblsib(i,4,5) &
            - ( pblsib(i,2,3) / pblsib(i,2,2) ) * pblsib(i,3,5)

       pblsib(i,1,5) =   ( pblsib(i,1,5) / pblsib(i,1,1) ) &
            - ( pblsib(i,1,4) / pblsib(i,1,1) ) * pblsib(i,4,5) &
            - ( pblsib(i,1,3) / pblsib(i,1,1) ) * pblsib(i,3,5) &
            - ( pblsib(i,1,2) / pblsib(i,1,1) ) * pblsib(i,2,5)
    END DO
    DO i = 1, nmax
       dtg(i) = pblsib(i,1,5) * dtc3x
       dtc(i) = pblsib(i,2,5) * dtc3x
       dtm(i) = pblsib(i,3,5) * dtc3x
       dqm(i) = pblsib(i,4,5) * dtc3x
       hc (i) = hc(i) + dtc3x * ( hcdtc(i) * dtc(i) + hcdtg(i) * dtg(i) + hcdtm(i) * dtm(i) )
       hg (i) = hg(i) + dtc3x * ( hgdtc(i) * dtc(i) + hgdtg(i) * dtg(i) + hgdtm(i) * dtm(i) )
       !
       !     check if interception loss term has exceeded canopy storage
       !
       ecpot(i)=( etc(i) - ea(i) ) + ( btc(i) - deadtc(i) ) * dtc(i) &
            -deadtg(i) * dtg(i) - deadqm(i) * dqm(i)
       egpot(i)=( etg(i) - ea(i) ) + ( btg(i) - deadtg(i) ) * dtg(i) &
            -deadtc(i) * dtc(i) - deadqm(i) * dqm(i)
    END DO
    !----------------------------------------------------------------------
    !     EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY
    !     ( hl*1000.) LOSS IS IN M M-2 (hl(J/kg))(1 J/kg ==> 1000J/m-3)
    !     MASS TERMS ARE IN KG M-2 DT-1
    !----------------------------------------------------------------------
    hlat3=1.0e+03_r8*hl
    DO i = 1, nmax
       eci   (i) = ecpot(i) * wc(i) * psyi(i) / ( 2.0_r8 * rb(i) ) * dtc3x
       ecidif(i) = MAX( 0.0_r8   , eci(i) - capac(i,1) * hlat3 )
       hc    (i) = hc(i) + ecidif(i)
       eci   (i) = MIN( eci(i) , capac(i,1) * hlat3 )
       egi   (i) = egpot(i) * vcover2(i,2) * wg(i) * psyi(i) / rd(i)*dtc3x
       egidif(i) = MAX( 0.0_r8 , egi(i) - capac(i,2) * hlat3 )
       hg    (i) = hg(i) + egidif(i)
       egi   (i) = MIN( egi(i) , capac(i,2) * hlat3 )
       !
       !     evaporation is given in j m-2, calculated from gradients
       !
       rsurf     = rsoil(i) * fg(i)
       coct      = ( 1.0_r8 - wc(i) )/rc(i)
       cogt      = vcover2(i,2) * ( 1.0_r8 - wg(i) ) / ( rg(i) + rd(i) )
       cogs1     = ( 1.0_r8 - vcover2(i,2) ) * hr(i)/( rd(i) + rsurf ) &
            + vcover2(i,2) / ( rd(i) + rsurf + 44.0_r8 ) * hr(i)
       cogs2     = cogs1/hr(i)
       ect  (i)  = ecpot(i)*coct*psyi(i)*dtc3x
       ec   (i)  = eci(i)+ect(i)
       egt  (i)  = egpot(i)*cogt*psyi(i)*dtc3x
       egs  (i)  = (etg(i)+btg(i)*dtg(i))*cogs1    &
            -(ea(i)+deadtg(i)*dtg(i)+deadtc(i)*dtc(i)+deadqm(i)*dqm(i) &
            )   *cogs2
       egs  (i)  = egs(i)*psyi(i)*dtc3x
       eg   (i)  = egt(i)+egs(i)+egi(i)
       !vcover2(i,2)=xcover(itype(i),mon(i),2)
    END DO
    !
    !     test of dew condition. recalculation ensues if necessary.
    !
    DO i = 1, nmax
       radt(i,1) = radt(i,1) - rcdtc(i) * dtc(i) - rcdtg(i) * dtg(i)
       radt(i,2) = radt(i,2) - rgdtc(i) * dtc(i) - rgdtg(i) * dtg(i)
       ecf    = SIGN(1.0_r8   ,ecpot(i)) * ( fc(i) * 2.0_r8 - 1.0_r8 )
       egf    = SIGN(1.0_r8   ,egpot(i)) * ( fg(i) * 2.0_r8 - 1.0_r8 )
       IF ( ecf <= 0.0_r8 ) THEN
          hc (i) = hc(i) + eci(i) + ect(i)
          eci(i) = 0.0_r8
          ect(i) = 0.0_r8
          ec (i) = 0.0_r8
       END IF
       IF (egf    <= 0.0_r8) THEN
          hg (i) = hg(i)+egi(i)+egt(i)+egs(i)
          egi(i) = 0.0_r8
          egt(i) = 0.0_r8
          egs(i) = 0.0_r8
          eg (i) = 0.0_r8
       END IF
    END DO
  END SUBROUTINE temres



  ! cut    :performs vapor pressure calculation at level "a".



  SUBROUTINE cut( &
       icheck,em    ,rhoair,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    , &
       ra    ,rb    ,rd    ,rc    ,rg    ,ea    ,etc   ,etg   ,rst   , &
       rsoil ,vcover,nmax  ,ncols )
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! input parameters
    !   fc      fg      hr      wc      wg      rhoair  cp
    !   rst     ra      rb      rg      rd      rsurf   vcover
    !   etc     etg     em
    !-----------------------------------------------------------------------
    ! output parameters
    !   ea
    !-----------------------------------------------------------------------
    ! ncols......Numero de ponto por faixa de latitude
    ! icg........Parametros da vegetacao (icg=1 topo e icg=2 base)
    ! cp.........specific heat of air (j/kg/k)
    ! nmax.......
    ! vcover(iv).Fracao de cobertura da vegetacao iv=1 topo ()
    ! vcover(iv).Fracao de cobertura da vegetacao iv=2 bottom ()
    ! ra.........Resistencia Aerodinamica (s/m)
    ! rb.........bulk boundary layer resistance             (s/m)
    ! rd.........aerodynamic resistance between ground
    !            and canopy air space                       (s/m)
    ! rc.........Resistencia do topo da copa
    ! rg.........Resistencia da base da copa
    ! ea.........Pressao de vapor
    ! etc........Pressao de vapor no topo da copa
    ! etg........Pressao de vapor no base da copa
    ! rst........Resistencia stomatal (s/m)
    ! rsoil......Resistencia do solo (s/m)
    ! em.........Pressao de vapor da agua
    ! rhoair.....Desnsidade do ar
    ! rcp........densidade do ar vezes o calor especifico do ar
    ! wc.........Minimo entre 1 e a razao entre a agua interceptada pelo
    !            indice de area foliar no topo da copa
    ! wg.........Minimo entre 1 e a razao entre a agua interceptada pelo
    !             indice de area foliar na parte inferior da copa
    ! fc.........Condicao de oravalho 0 ou 1 no topo da copa
    ! fg.........Condicao de oravalho 0 ou 1 na base da copa
    ! hr.........Rel. humidity in top layer
    ! icheck
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: nmax
    REAL(KIND=r8),    INTENT(in   ) :: vcover(ncols,icg)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(in   ) :: ra    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rb    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rd    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ea    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: etc   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: etg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rst   (ncols,icg)
    REAL(KIND=r8),    INTENT(in   ) :: rsoil (ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: em    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rhoair(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rcp   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: wc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: wg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: fc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: fg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: hr    (ncols)
    INTEGER, INTENT(in   ) :: icheck(ncols)

    REAL(KIND=r8) :: coc
    REAL(KIND=r8) :: rsurf
    REAL(KIND=r8) :: cog1
    REAL(KIND=r8) :: cog2
    REAL(KIND=r8) :: d2
    REAL(KIND=r8) :: top
    REAL(KIND=r8) :: xnum
    REAL(KIND=r8) :: tem
    INTEGER :: i

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          rcp  (i) = rhoair(i)*cp
          rc   (i) = rst(i,1)*fc(i)+rb(i)+rb(i)*fc(i)
          coc      = (1.0_r8 -wc(i))/rc(i)+wc(i)/(2.0_r8 *rb(i))
          rg   (i) = rst(i,2)*fg(i)
          rsurf    = rsoil(i)*fg(i)
          tem      = vcover(i,2)*(1.0_r8-wg(i))/(rg(i)+rd(i))
          cog2     = tem    &
               + (1.0_r8 -vcover(i,2))/(rsurf   +rd(i)) &
               + vcover(i,2)/(rsurf   +rd(i)+44.0_r8)
          cog1     = (cog2   -tem   )*hr(i)+tem
          xnum     = wg(i)/rd(i)*vcover(i,2)
          cog1     = cog1   +xnum
          cog2     = cog2   +xnum
          d2       = 1.0_r8 /ra(i)+coc+cog2
          top      = coc*etc(i)+em(i)/ra(i)+cog1   *etg(i)
          !
          !     vapor pressure at level "a"
          !
          ea (i)  = top   /d2
       END IF
    END DO
  END SUBROUTINE cut



  ! rbrd   :calculates bulk boundary layer resistance and aerodynamic
  !         resistence betweenground and canopi air space4 as functions
  !         of wind speed at top of canopy and temperatures.



  SUBROUTINE rbrd(rb    ,rd    ,tcta  ,tgta  ,u2    ,tg    ,rdc   ,rbc   ,itype , &
       z2    ,mon   ,nmax  ,ncols , zlt2)
    !
    !
    !         rb and rd as functions of u2 and temperatures. simplified( xue et
    !         al. 1991)
    !
    !-----------------------------------------------------------------------
    !       input parameters
    !-----------------------------------------------------------------------
    !   tcta..........diferenca entre tc-ta                      (k)
    !   tgta..........diferenca entre tg-ta                      (k)
    !   tg............ground temperature                         (k)
    !   u2............wind speed at top of canopy                (m/s)
    !   z2............height of canopy top                       (m)
    !   zlt(cg).......canopy/ground cover leaf and stem area density
    !                                                            (m**2/m**3)
    !   rbc...........constant related to bulk boundary layer
    !                 resistance
    !   rdc...........constant related to aerodynamic resistance
    !                 between ground and canopy air space
    !-----------------------------------------------------------------------
    !      output parameters
    !-----------------------------------------------------------------------
    !   rb............bulk boundary layer resistance             (s/m)
    !   rd............aerodynamic resistance between ground      (s/m)
    !                 and canopy air space
    !-----------------------------------------------------------------------
    !   ncols.........Numero de ponto por faixa de latitude
    !   ityp..........numero das classes de solo 13
    !   imon...........Numero maximo de meses no ano (12)
    !   icg...........Parametros da vegetacao (icg=1 topo e icg=2 base)
    !   mon...........Numero do mes do ano (1-12)
    !   nmax .........
    !   itype.........Classe de textura do solo
    !=======================================================================
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: mon (ncols)
    INTEGER, INTENT(in   ) :: nmax
    !
    !     vegetation and soil parameters
    !
    REAL(KIND=r8),    INTENT(in   ) :: z2    (ityp,imon)
    INTEGER, INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rdc   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rbc   (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(in   ) :: tg   (ncols)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout  ) :: rb    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rd    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tcta  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tgta  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: u2    (ncols)
    REAL(KIND=r8),    INTENT(in  ) :: zlt2(ncols,icg)

    REAL(KIND=r8) :: temdif(ncols)
    REAL(KIND=r8) :: fih   (ncols)

    REAL(KIND=r8), PARAMETER :: factg=88.29_r8
    INTEGER :: i
    INTEGER :: ntyp

    DO i = 1, nmax
       ntyp=itype(i)
       IF (tcta(i) > 0.0_r8 ) THEN
          temdif(i)=tcta(i)+0.1_r8
       ELSE
          temdif(i)=        0.1_r8
       END IF
       rb (i)=1.0_r8  /(SQRT(u2(i))/rbc(i)+zlt2(i,1)*0.004_r8 )
       IF (tgta(i) > 0) THEN
          temdif(i)=tgta(i)+0.1_r8
       ELSE
          temdif(i)=        0.1_r8
       END IF
       fih(i)=sqrt &
            (1.0_r8 +factg*temdif(i)*z2(ntyp,mon(i))/(tg(i)*u2(i)*u2(i)))
       rd(i) =rdc(i)/(u2(i)*fih(i))
    END DO
  END SUBROUTINE rbrd



  ! vntlax :performs ventilation mass flux, based on deardorff, mwr, 1972?.


  SUBROUTINE vntlax(ustarn, &
       icheck,bps   ,dzm   ,cu    ,ct,cuni  ,ctni  ,ustar ,ra    ,ta    , &
       u2    ,tm    ,um    ,vm    ,d     ,z0    ,itype ,z2    , &
       mon   ,nmax  ,jstneu,ncols )
    !
    !
    !-----------------------------------------------------------------------
    !         input parameters
    !-----------------------------------------------------------------------
    !
    !   ea..........Pressao de vapor
    !   ta..........Temperatura no nivel de fonte de calor do dossel (K)
    !   um..........Razao entre zonal pseudo-wind (fourier) e seno da
    !               colatitude
    !   vm..........Razao entre meridional pseudo-wind (fourier) e seno da
    !               colatitude
    !   qm..........specific humidity of reference (fourier)
    !   tm..........Temperature of reference (fourier)
    !   dzm  .......Altura media de referencia  para o vento para o calculo
    !               da estabilidade do escoamento
    !   grav........gravity constant      (m/s**2)
    !   cpair.......specific heat of air (j/kg/k)
    !   gasr........gas constant of dry air      (j/kg/k)
    !   bps ........
    !   z2..........height of canopy top
    !   d...........displacement height                        (m)
    !   epsfac......parametro para o gas 0.622
    !
    !
    !
    !-----------------------------------------------------------------------
    !        output parameters
    !-----------------------------------------------------------------------
    !
    !   ustar.........surface friction velocity  (m/s)
    !   ra............Resistencia Aerodinamica (s/m)
    !   u2............wind speed at top of canopy                (m/s)
    !   ventmf........ventilation mass flux
    !-----------------------------------------------------------------------
    !=======================================================================
    !   ncols........Numero de ponto por faixa de latitude
    !   ityp.........Numero do tipo de solo
    !   imon.........Numero maximo de meses no ano (12)
    !   jstneu.......The first call to vntlat just gets the neutral values
    !                of ustar and ventmf para jstneu=.TRUE..
    !   mon..........Numero do mes do ano (1-12)
    !   nmax.........
    !   itype........Classe de textura do solo
    !   z0...........roughness length
    !   bps..........bps   (i)=sigki(1)=1.0e0/EXP(akappa*LOG(sig(k)))
    !   cu...........friction  transfer coefficients.
    !   ct...........heat transfer coefficients.
    !   cuni.........neutral friction transfer  coefficients.
    !   ctni.........neutral heat transfer coefficients.
    !   icheck.......this version assumes dew-free conditions "icheck=1" to
    !                estimate ea for buoyancy term in vntmf or ra.
    !=======================================================================
    INTEGER, INTENT(in   ) :: ncols

    LOGICAL, INTENT(in   ) :: jstneu
    INTEGER, INTENT(in   ) :: mon(ncols)
    INTEGER, INTENT(in   ) :: nmax
    !
    !     vegetation and soil parameters
    !
    REAL(KIND=r8),    INTENT(in   ) :: z2    (ityp,imon)
    INTEGER, INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: d     (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: z0    (ncols)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(in   ) :: tm  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: um  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vm  (ncols)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout) :: ra    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: u2    (ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: bps   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    INTEGER, INTENT(in   ) :: icheck(ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustarn(ncols)


    !REAL(KIND=r8) :: thm(ncols)      !**(JP)** scalar
    REAL(KIND=r8) :: thm
    !REAL(KIND=r8) :: ros(ncols)      !**(JP)** unused
    REAL(KIND=r8) :: speedm(ncols)
    !REAL(KIND=r8) :: thvgm(ncols)    !**(JP)** scalar
    REAL(KIND=r8) :: thvgm
    !REAL(KIND=r8) :: rib(ncols)      !**(JP)** scalar
    REAL(KIND=r8) :: rib
    !REAL(KIND=r8) :: cui(ncols)      !**(JP)** scalar
    REAL(KIND=r8) :: cui
    !REAL(KIND=r8) :: ran(ncols)      !**(JP)** unused
    !REAL(KIND=r8) :: cti(ncols)      !**(JP)** scalar
    REAL(KIND=r8) :: cti
    !REAL(KIND=r8) :: ct      (ncols)      !**(JP)** unused

    REAL(KIND=r8), PARAMETER ::  vkrmn=0.40_r8
    REAL(KIND=r8), PARAMETER ::  fsc=66.85_r8
    REAL(KIND=r8), PARAMETER ::  ftc=0.904_r8
    REAL(KIND=r8), PARAMETER ::  fvc=0.315_r8
    REAL(KIND=r8) :: rfac
    REAL(KIND=r8) :: vkrmni
    REAL(KIND=r8) :: g2
    REAL(KIND=r8) :: zl
    REAL(KIND=r8) :: xct1
    REAL(KIND=r8) :: xct2
    REAL(KIND=r8) :: xctu1
    REAL(KIND=r8) :: xctu2
    REAL(KIND=r8) :: grib
    REAL(KIND=r8) :: grzl
    REAL(KIND=r8) :: grz2
    REAL(KIND=r8) :: fvv
    REAL(KIND=r8) :: ftt
    REAL(KIND=r8) :: rzl
    REAL(KIND=r8) :: rz2
    INTEGER :: i
    INTEGER :: ntyp


    rfac  =1.0e2_r8 /gasr

    vkrmni=1.0_r8  /vkrmn
    g2 = 0.75_r8

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          speedm(i)=SQRT(um(i)**2+vm(i)**2)
          speedm(i)=MAX(2.0_r8  ,speedm(i))
       END IF
    END DO
    !
    !     cu and ct are the friction and heat transfer coefficients.
    !     cun and ctn are the neutral friction and heat transfer
    !     coefficients.
    !
    IF (jstneu) THEN
       DO i = 1, nmax
          ntyp=itype(i)
          zl = z2(ntyp,mon(i)) + 11.785_r8  * z0(i)
!ANNE          cuni(i)=LOG((dzm(i)-d(i))/z0(i))*vkrmni
          cuni(i)=LOG(MAX((dzm(i)-d(i))/z0(i),0.0001_r8))*vkrmni
          ustarn(i)=speedm(i)/cuni(i)
          IF (zl < dzm(i)) THEN
!ANNE              xct1 = LOG((dzm(i)-d(i))/(zl-d(i)))
!ANNE              xct2 = LOG((zl-d(i))/z0(i))
             xct1 = LOG(MAX((dzm(i)-d(i))/(zl-d(i)),0.0001_r8))
             xct2 = LOG(MAX((zl-d(i))/z0(i),0.0001_r8))
             xctu1 = xct1
!ANNE              xctu2 = LOG((zl-d(i))/(z2(ntyp,mon(i))-d(i)))
             xctu2 = LOG(MAX((zl-d(i))/(z2(ntyp,mon(i))-d(i)),0.0001_r8))
             ctni(i) = (xct1 + g2 * xct2) *vkrmni
          ELSE
!ANNE                 xct2 =  LOG((dzm(i)-d(i))/z0(i))
             xct2 =  LOG(MAX((dzm(i)-d(i))/z0(i),0.0001_r8))
             xctu1 =  0.0_r8
!ANNE             xctu2 =  LOG((dzm(i)-d(i))/(z2(ntyp,mon(i))-d(i))) 
             xctu2 =  LOG(MAX((dzm(i)-d(i))/(z2(ntyp,mon(i))-d(i)),0.0001_r8))
             ctni(i) = g2 * xct2 *vkrmni
          END IF
          !
          !     neutral values of ustar and ventmf
          !
!ANNE            u2(i) = speedm(i) - ustarn(i)*vkrmni*(xctu1 + g2*xctu2)
          u2(i) = max(speedm(i) - ustarn(i)*vkrmni*(xctu1 + g2*xctu2),0.01_r8)

       END DO
       RETURN
    END IF
    !
    !     stability branch based on bulk richardson number.
    !
    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          !
          !     freelm(i)=.false.
          !
          thm= tm(i)*bps(i)
          ntyp=itype(i)
          zl = z2(ntyp,mon(i)) + 11.785_r8  * z0(i)
          thvgm   = ta(i)-thm
!ANNE           rib     =-thvgm   *grav*(dzm(i)-d(i)) &
          rib     =-thvgm   *grav*(MAX(dzm(i)-d(i),0.0001_r8)) &
               /(thm*(speedm(i)-u2(i))**2)
          ! Manzi Suggestion:
          ! rib   (i)=max(-10.0_r8  ,rib(i))
          rib      =MAX(-1.5_r8  ,rib   )
          rib      =MIN( 0.165_r8  ,rib   )
          IF (rib    < 0.0_r8) THEN
             grib = -rib
!ANNE             grzl = -rib   * (zl-d(i))/(dzm(i)-d(i))
!ANNE             grz2 = -rib   * z0(i)/(dzm(i)-d(i))
             grzl = -rib   * (zl-d(i))/(MAX(dzm(i)-d(i),0.0001_r8))
             grz2 = -rib   * z0(i)/(MAX(dzm(i)-d(i),0.0001_r8))
             fvv = fvc*grib
             IF (zl < dzm(i)) THEN
                ftt = (ftc*grib) + (g2-1.0_r8) * (ftc*grzl) - g2 * (ftc*grz2)
             ELSE
                ftt = g2*((ftc*grib) - (ftc*grz2))
             END IF
             cui    = cuni(i) - fvv
             cti    = ctni(i) - ftt
          ELSE
!ANNE             rzl = rib   /(dzm(i)-d(i))*(zl-d(i))
!ANNE             rz2 = rib   /(dzm(i)-d(i))*z0(i)
             rzl = rib   /(MAX(dzm(i)-d(i),0.0001_r8))*(zl-d(i))
             rz2 = rib   /(MAX(dzm(i)-d(i),0.0001_r8))*z0(i)
             fvv = fsc*rib
             IF (zl < dzm(i)) THEN
                ftt = (fsc*rib) + (g2-1) * (fsc*rzl) - g2 * (fsc*rz2)
             ELSE
                ftt = g2 * ((fsc*rib) - (fsc*rz2))
             END IF
             cui    = cuni(i) + fvv
             cti    = ctni(i) + ftt
          ENDIF
          cu    (i)=1.0_r8/cui
          !**(JP)** ct is not used anywhere else
          ct    (i)=1.0_r8/cti
          !
          !
          !     surface friction velocity and ventilation mass flux
          !
          ustar (i)=speedm(i)*cu(i)
          ra(i) = cti    / ustar(i)
          !**(JP)** ran is not used anywhere else
          !ran(i) = ctni(i) / ustarn(i)
          !ran(i) = MAX(ran(i), 0.8_r8 )
          ra(i) = MAX(ra(i), 0.8_r8 )
       END IF
    END DO
  END SUBROUTINE vntlax




  ! runoff :performs inter-layer moisture exchanges.



  SUBROUTINE runoff( &
       roff  ,tg    ,td    ,capac ,w     ,itype ,dtc3x ,nmax  ,ncols    )
    !
    !-----------------------------------------------------------------------
    ! input parameters
    !-----------------------------------------------------------------------
    !   w(3)     roff     slope    bee      satco     zdepth
    !   phsat    poros    pie      dtc3x    snomel
    !   w(3)
    !
    !-----------------------------------------------------------------------
    ! output parameters
    !-----------------------------------------------------------------------
    !   w(3)     roff
    !-----------------------------------------------------------------------
    !
    ! roff.......Runoff (escoamente superficial e drenagem)(m)
    ! slope......Inclinacao de perda hidraulica na camada profunda do solo
    ! bee........Fator de retencao da umidade no solo (expoente da umidade do
    !            solo)
    ! satco......Condutividade hidraulica do solo saturado(m/s)
    ! zdepth(id).Profundidade das camadas de solo id=1 superficial
    ! zdepth(id).Profundidade das camadas de solo id=2 camada de raizes
    ! zdepth(id).Profundidade das camadas de solo id=3 camada de drenagem
    ! phsat......Potencial matricial do solo saturado(m) (tensao do solo em
    !            saturacao)
    ! poros......Porosidade do solo
    ! pie........pi = 3.1415926e0
    ! dtc3x......time increment dt
    ! snomel.....Calor latente de fusao(J/kg)
    ! w(id)......Grau de saturacao de umidade do solo id=1 na camada superficial
    ! w(id)......Grau de saturacao de umidade do solo id=2 na camada de raizes
    ! w(id)......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    ! capac(iv)..Agua interceptada iv=1 no dossel (m)
    ! capac(iv)..Agua interceptada iv=2 na cobertura do solo (m)
    ! tg.........Temperatura da superficie do solo  (K)
    ! td.........Temperatura do solo profundo (K)
    ! itype......Classe de textura do solo
    ! tf.........Temperatura de congelamento (K)
    ! idp........Parametro para as camadas de solo idp=1->3
    ! nmax.......
    ! ncols......Numero de ponto por faixa de latitude
    ! ityp.......13
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: ncols

    REAL(KIND=r8),    INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: nmax

    INTEGER, INTENT(in   ) :: itype (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(in   ) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: td   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: capac(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: w    (ncols,3)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8),    INTENT(inout) :: roff (ncols)

    REAL(KIND=r8)    :: q3g   (ncols)
    REAL(KIND=r8)    :: div   (ncols)
    REAL(KIND=r8)    :: twi   (ncols,3)
    REAL(KIND=r8)    :: twip  (ncols,3)
    REAL(KIND=r8)    :: twipp (ncols,3)
    REAL(KIND=r8)    :: avk   (ncols)
    REAL(KIND=r8)    :: aaa_1, aaa_2
    REAL(KIND=r8)    :: bbb_1, bbb_2
    REAL(KIND=r8)    :: ccc_1, ccc_2
    REAL(KIND=r8)    :: qqq_1, qqq_2

    REAL(KIND=r8)    :: subdt
    REAL(KIND=r8)    :: subdti
    REAL(KIND=r8)    :: slop
    REAL(KIND=r8)    :: pows
    REAL(KIND=r8)    :: wmax
    REAL(KIND=r8)    :: wmin
    REAL(KIND=r8)    :: pmax
    REAL(KIND=r8)    :: pmin
    REAL(KIND=r8)    :: dpdw
    REAL(KIND=r8)    :: rsame
    REAL(KIND=r8)    :: tsnow
    REAL(KIND=r8)    :: areas
    REAL(KIND=r8)    :: tgs
    REAL(KIND=r8)    :: ts
    REAL(KIND=r8)    :: props
    REAL(KIND=r8)    :: dpdwdz
    REAL(KIND=r8)    :: denom
    REAL(KIND=r8)    :: rdenom
    REAL(KIND=r8)    :: qmax
    REAL(KIND=r8)    :: qmin
    REAL(KIND=r8)    :: excess
    REAL(KIND=r8)    :: deficit
    INTEGER :: n
    INTEGER :: i
    INTEGER :: ntyp
    REAL(KIND=r8),    PARAMETER     :: smal2 = 1.0e-3_r8

    subdt =dtc3x
    subdti=1.0_r8 /dtc3x
    q3g=0.0_r8
    !
    !     eliinate negative soil moisture
    !
    DO n = 1, nmax
       IF (w(n,1) < 0.0_r8) w(n,1)=smal2
       IF (w(n,2) < 0.0_r8) w(n,2)=smal2
       IF (w(n,3) < 0.0_r8) w(n,3)=smal2
    END DO

    DO i = 1, 3
       DO n = 1, nmax
          ntyp      =itype(n)
          twi(n,i)=MIN(1.0_r8, MAX(0.03_r8,w(n,i)))
          twip(n,i) =EXP(-bee(ntyp)*LOG(twi(n,i)))
          twipp(n,i)=EXP((2.0_r8*bee(ntyp)+3.0_r8)*LOG(MIN(1.0_r8,twi(n,i))))
       END DO
    END DO

    DO n = 1, nmax
       ntyp = itype(n)
       slop = 0.1736_r8
       IF (poros(ntyp) == 0.4352_r8) slop = 0.0872_r8
       IF (poros(ntyp) == 0.4577_r8) slop = 0.3420_r8
       !
       !     calculation of gravitationally driven drainage from w(3) : taken
       !     as an integral of time varying conductivity.addition of liston
       !     baseflow term to original q3g to insure flow in
       !     dry season. modified liston baseflow constant scaled
       !     by available water.
       !
       !     q3g (q3) : equation (62) , se-86
       !
       pows    = 2.0_r8 *bee(ntyp)+2.0_r8
       q3g (n) = EXP(-pows*LOG(twi(n,3))) &
            +satco(ntyp)/(zdepth(ntyp,3)*poros(ntyp))* &
            slop*pows*subdt
       q3g (n) = EXP(LOG(q3g(n))/pows)
       q3g (n) =-(1.0_r8 /q3g(n)-w(n,3)) &
            *poros(ntyp)*zdepth(ntyp,3)*subdti
       q3g (n) = MAX(0.0_r8 ,q3g(n))
       q3g (n) = MIN(q3g(n), w(n,3)*poros(ntyp)*zdepth(ntyp,3) &
            *subdti)
       q3g (n) = q3g(n)+0.002_r8*poros(ntyp)*zdepth(ntyp,3)*0.5_r8 &
            /86400.0_r8*w(n,3)
    END DO
    !
    !     calculation of inter-layer exchanges of water due to gravitation
    !     and hydraulic gradient. the values of w(x) + dw(x) are used to
    !     calculate the potential gradients between layers.
    !     modified calculation of mean conductivities follows milly and
    !     eagleson (1982 ), reduces recharge flux to top layer.
    !
    !      dpdw           : estimated derivative of soil moisture potential
    !                       with respect to soil wetness. assumption of
    !                       gravitational drainage used to estimate likely
    !                       minimum wetness over the time step.
    !
    !      qqq  (q     )  : equation (61) , s-86
    !             i,i+1
    !            -
    !      avk  (k     )  : equation (4.14) , milly and eagleson (1982)
    !             i,i+1
    !
    DO n = 1, nmax
       ntyp=itype(n)
       wmax = MAX( w(n,1), w(n,2), w(n,3), 0.05_r8 )
       wmax = MIN( wmax, 1.0_r8 )
       pmax = EXP(-bee(ntyp)*LOG(wmax))
       wmin = EXP(-1.0_r8/bee(ntyp)*LOG(pmax-2.0_r8/(phsat(ntyp) &
            *(zdepth(ntyp,1)+2.0_r8*zdepth(ntyp,2)+zdepth(ntyp,3)))))
       wmin = MIN( w(n,1), w(n,2), w(n,3), wmin )
       wmin = MAX( wmin, 0.02_r8 )
       pmin = EXP(-bee(ntyp)*LOG(wmin))
       dpdw = phsat(ntyp)*( pmax-pmin )/( wmax-wmin )

       ! hand unrolling of next do loop, first iteration

       rsame = 0.0_r8
       avk(n)    =twip(n,1)*twipp(n,1)-twip(n,1+1)*twipp(n,1+1)
       div(n)    =twip(n,1+1) - twip(n,1)
       IF(ABS(div(n)).LE.1.0e-7_r8) rsame = 1.0_r8
       avk(n)=satco(ntyp)*avk(n)/ &
            ((1.0_r8 +3.0_r8 /bee(ntyp))*div(n)+ rsame)
       avk(n)=MAX(avk(n),satco(ntyp)*MIN(twipp(n,1),twipp(n,1+1)))
       avk(n)=MIN(avk(n),1.01_r8*(satco(ntyp) &
            *MAX(twipp(n,1),twipp(n,1+1))))
       !
       !     conductivities and base flow reduced when temperature drops below
       !     freezing
       !
       tsnow = MIN (tf-0.01_r8, tg(n))
       areas = MIN (0.999_r8,13.2_r8*capac(n,2))
       tgs = tsnow*areas + tg(n)*(1.0_r8-areas)
       ts  = tgs*(2-1) + td(n)*(1-1)
       props = (ts-(tf-10.0_r8))/10.0_r8
       props = MAX(0.05_r8,MIN(1.0_r8, props))
       avk(n) = avk(n) * props
       q3g(n) = q3g(n) * props
       !
       !     backward implicit calculation of flows between soil layers
       !
       dpdwdz= dpdw * 2.0_r8/( zdepth(ntyp,1) + zdepth(ntyp,1+1) )
       aaa_1=1.0_r8+avk(n)*dpdwdz* &
            (1.0_r8/zdepth(ntyp,1)+1.0_r8/zdepth(ntyp,1+1))      &
            *subdt/poros(ntyp)
       bbb_1 =-avk(n)* dpdwdz*1.0_r8/zdepth(ntyp,2)*subdt/poros(ntyp)
       ccc_1 = avk(n) * (dpdwdz * ( w(n,1)-w(n,1+1) )+1.0_r8+(1-1) &
            *dpdwdz*q3g(n)*1.0_r8/zdepth(ntyp,3)*subdt/poros(ntyp))

       ! hand unrolling of next do loop, second iteration

       rsame = 0.0_r8
       avk(n)    =twip(n,2)*twipp(n,2)-twip(n,2+1)*twipp(n,2+1)
       div(n)    =twip(n,2+1) - twip(n,2)
       IF(ABS(div(n)).LE.1.0e-7_r8) rsame = 1.0_r8
       avk(n)=satco(ntyp)*avk(n)/ &
            ((1.0_r8 +3.0_r8 /bee(ntyp))*div(n)+ rsame)
       avk(n)=MAX(avk(n),satco(ntyp)*MIN(twipp(n,2),twipp(n,2+1)))
       avk(n)=MIN(avk(n),1.01_r8*(satco(ntyp) &
            *MAX(twipp(n,2),twipp(n,2+1))))
       !
       !     conductivities and base flow reduced when temperature drops below
       !     freezing
       !
       tsnow = MIN (tf-0.01_r8, tg(n))
       areas = MIN (0.999_r8,13.2_r8*capac(n,2))
       tgs = tsnow*areas + tg(n)*(1.0_r8-areas)
       ts  = tgs*(2-2) + td(n)*(2-1)
       props = (ts-(tf-10.0_r8))/10.0_r8
       props = MAX(0.05_r8,MIN(1.0_r8, props))
       avk(n) = avk(n) * props
       q3g(n) = q3g(n) * props
       !
       !     backward implicit calculation of flows between soil layers
       !
       dpdwdz= dpdw * 2.0_r8/( zdepth(ntyp,2) + zdepth(ntyp,2+1) )
       aaa_2=1.0_r8+avk(n)*dpdwdz* &
            (1.0_r8/zdepth(ntyp,2)+1.0_r8/zdepth(ntyp,2+1))      &
            *subdt/poros(ntyp)
       bbb_2 =-avk(n)* dpdwdz*1.0_r8/zdepth(ntyp,2)*subdt/poros(ntyp)
       ccc_2 = avk(n) * (dpdwdz * ( w(n,2)-w(n,2+1) )+1.0_r8+(2-1) &
            *dpdwdz*q3g(n)*1.0_r8/zdepth(ntyp,3)*subdt/poros(ntyp))


       !       DO i = 1, 2
       !          rsame = 0.0_r8
       !          avk(n)    =twip(n,i)*twipp(n,i)-twip(n,i+1)*twipp(n,i+1)
       !          div(n)    =twip(n,i+1) - twip(n,i)
       !          IF(ABS(div(n)).LE.1.0e-7_r8) rsame = 1.0_r8
       !          avk(n)=satco(ntyp)*avk(n)/ &
       !               ((1.0_r8 +3.0_r8 /bee(ntyp))*div(n)+ rsame)
       !          avk(n)=MAX(avk(n),satco(ntyp)*MIN(twipp(n,i),twipp(n,i+1)))
       !          avk(n)=MIN(avk(n),1.01_r8*(satco(ntyp) &
       !               *MAX(twipp(n,i),twipp(n,i+1))))
       !          !
       !          !     conductivities and base flow reduced when temperature drops below
       !          !     freezing
       !          !
       !          tsnow = MIN (tf-0.01_r8, tg(n))
       !          areas = MIN (0.999_r8,13.2_r8*capac(n,2))
       !          tgs = tsnow*areas + tg(n)*(1.0_r8-areas)
       !          ts  = tgs*(2-i) + td(n)*(i-1)
       !          props = (ts-(tf-10.0_r8))/10.0_r8
       !          props = MAX(0.05_r8,MIN(1.0_r8, props))
       !          avk(n) = avk(n) * props
       !          q3g(n) = q3g(n) * props
       !          !
       !          !     backward implicit calculation of flows between soil layers
       !          !
       !          dpdwdz= dpdw * 2.0_r8/( zdepth(ntyp,i) + zdepth(ntyp,i+1) )
       !          aaa(i)=1.0_r8+avk(n)*dpdwdz* &
       !               (1.0_r8/zdepth(ntyp,i)+1.0_r8/zdepth(ntyp,i+1))      &
       !               *subdt/poros(ntyp)
       !          bbb(i) =-avk(n)* dpdwdz*1.0_r8/zdepth(ntyp,2)*subdt/poros(ntyp)
       !          ccc(i) = avk(n) * (dpdwdz * ( w(n,i)-w(n,i+1) )+1.0_r8+(i-1) &
       !               *dpdwdz*q3g(n)*1.0_r8/zdepth(ntyp,3)*subdt/poros(ntyp))
       !       END DO
       denom    = ( aaa_1*aaa_2 - bbb_1*bbb_2 )
       rdenom   = 0.0_r8
       IF (ABS(denom) < 1.e-6_r8 ) rdenom = 1.0_r8
       rdenom   = ( 1.0_r8-rdenom)/( denom + rdenom )
       qqq_1   = ( aaa_2*ccc_1 - bbb_1*ccc_2 ) * rdenom
       qqq_2   = ( aaa_1*ccc_2 - bbb_2*ccc_1 ) * rdenom
       !
       !     update wetness of each soil moisture layer due to layer interflow
       !     and base flow.
       !
       w(n,3)  = w(n,3) - q3g(n)*subdt/(poros(ntyp)*zdepth(ntyp,3))
       roff(n) = roff(n) + q3g(n) * subdt

       ! hand unrolling of next do loop, first iteration

       qmax     =  w(n,1)   * (poros(ntyp)*zdepth(ntyp,1)  /subdt)
       qmin     = -w(n,1+1) * (poros(ntyp)*zdepth(ntyp,1+1)/subdt)
       qqq_1   =  MIN( qqq_1,qmax)
       qqq_1   =  MAX( qqq_1,qmin)
       w(n,1)   =  w(n,1)  -qqq_1/(poros(ntyp)*zdepth(ntyp,1) /subdt)
       w(n,1+1) =  w(n,1+1)+ &
            qqq_1/(poros(ntyp)*zdepth(ntyp,1+1)/subdt)

       ! hand unrolling of next do loop, second iteration

       qmax     =  w(n,2)   * (poros(ntyp)*zdepth(ntyp,2)  /subdt)
       qmin     = -w(n,2+1) * (poros(ntyp)*zdepth(ntyp,2+1)/subdt)
       qqq_2   =  MIN( qqq_2,qmax)
       qqq_2   =  MAX( qqq_2,qmin)
       w(n,2)   =  w(n,2)  -qqq_2/(poros(ntyp)*zdepth(ntyp,2) /subdt)
       w(n,2+1) =  w(n,2+1)+ &
            qqq_2/(poros(ntyp)*zdepth(ntyp,2+1)/subdt)
       !     DO i = 1, 2
       !        qmax     =  w(n,i)   * (poros(ntyp)*zdepth(ntyp,i)  /subdt)
       !        qmin     = -w(n,i+1) * (poros(ntyp)*zdepth(ntyp,i+1)/subdt)
       !        qqq(i)   =  MIN( qqq(i),qmax)
       !        qqq(i)   =  MAX( qqq(i),qmin)
       !        w(n,i)   =  w(n,i)  -qqq(i)/(poros(ntyp)*zdepth(ntyp,i) /subdt)
       !        w(n,i+1) =  w(n,i+1)+ &
       !             qqq(i)/(poros(ntyp)*zdepth(ntyp,i+1)/subdt)
       !     END DO

       ! hand unrolling of next do loop, first iteration

       excess   = MAX(0.0_r8,(w(n,1) - 1.0_r8))
       w(n,1)   = w(n,1) - excess
       roff(n)  = roff(n) + excess * poros(ntyp)*zdepth(ntyp,1)

       ! hand unrolling of next do loop, second iteration

       excess   = MAX(0.0_r8,(w(n,2) - 1.0_r8))
       w(n,2)   = w(n,2) - excess
       roff(n)  = roff(n) + excess * poros(ntyp)*zdepth(ntyp,2)

       ! hand unrolling of next do loop, third iteration

       excess   = MAX(0.0_r8,(w(n,3) - 1.0_r8))
       w(n,3)   = w(n,3) - excess
       roff(n)  = roff(n) + excess * poros(ntyp)*zdepth(ntyp,3)

       !     DO i = 1, 3
       !        excess   = MAX(0.0_r8,(w(n,i) - 1.0_r8))
       !        w(n,i)   = w(n,i) - excess
       !        roff(n)  = roff(n) + excess * poros(ntyp)*zdepth(ntyp,i)
       !     END DO

       ! hand unrolling of next do loop, first iteration

       deficit   = MAX (0.0_r8,(1.e-12_r8 - w(n,1)))
       w(n,1)    = w(n,1) + deficit
       w(n,1+1)  = w(n,1+1)-deficit*zdepth(ntyp,1)/zdepth(ntyp,1+1)

       ! hand unrolling of next do loop, second iteration

       deficit   = MAX (0.0_r8,(1.e-12_r8 - w(n,2)))
       w(n,2)    = w(n,2) + deficit
       w(n,2+1)  = w(n,2+1)-deficit*zdepth(ntyp,2)/zdepth(ntyp,2+1)

       !
       !     prevent negative values of www(i)
       !

       !       DO i = 1,2
       !          deficit   = MAX (0.0_r8,(1.e-12_r8 - w(n,i)))
       !          w(n,i)    = w(n,i) + deficit
       !          w(n,i+1)  = w(n,i+1)-deficit*zdepth(ntyp,i)/zdepth(ntyp,i+1)
       !       END DO

       w(n,3)      = MAX (w(n,3),1.0e-12_r8)
    END DO
  END SUBROUTINE runoff


  ! stres2 :calculates the adjustment to light dependent stomatal resistance
  !         by temperature, humidity and stress factors (simplified).



  SUBROUTINE stres2( &
       icount,ft1   ,fp1   ,icheck,ta    ,ea    ,rst   ,phsoil,stm   , &
       tc    ,tg    ,w     ,vcover,itype , &
       rootd ,zdepth,nmax  ,ncols ,topt2 ,tll2  ,tu2   , &
       defac2,ph12  ,ph22)
    !
    !
    !-----------------------------------------------------------------------
    ! ityp........numero das classes de solo 13
    ! icg.........Parametros da vegetacao (icg=1 topo e icg=2 base)
    ! idp.........Parametro para as camadas de solo idp=1->3
    ! icount......
    ! ft1.........temperature  factor   simplified
    ! fp1.........soil water potential factor simplified
    ! hl........heat of evaporation of water   (j/kg)
    ! nmax........
    ! topt........Temperatura ideal de funcionamento estomatico
    ! tll.........Temperatura minima de funcionamento estomatico
    ! tu..........Temperatura maxima de funcionamento estomatico
    ! defac.......Parametro de deficit de pressao de vapor d'agua
    ! ph1.........Coeficiente para o efeito da agua no solo
    ! ph2 ........Potencial de agua no solo para ponto de Wilting
    ! rootd.......Profundidade das raizes
    ! zdepth......Profundidade para as tres camadas de solo
    ! itype.......Classe de textura do solo
    ! vcover(iv)..Fracao de cobertura de vegetacao iv=1 Top
    ! vcover(iv)..Fracao de cobertura de vegetacao iv=2 Bottom
    ! tc..........Temperatura da copa "dossel"(K)
    ! tg .........Temperatura da superficie do solo (K)
    ! w(id).......Grau de saturacao de umidade do solo id=1 na camada superficial
    ! w(id).......Grau de saturacao de umidade do solo id=2 na camada de raizes
    ! w(id).......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    ! ta..........Temperatura no nivel de fonte de calor do dossel (K)
    ! ea..........Pressao de vapor
    ! rst ........Resisttencia Estomatica "Stomatal resistence" (s/m)
    ! phsoil......soil moisture potential of the i-th soil layer
    ! stm.........Resisttencia Estomatica "Stomatal resistence" (s/m)
    ! icheck......this version assumes dew-free conditions "icheck=1" to
    !             estimate ea for buoyancy term in vntmf or ra.
    !-----------------------------------------------------------------------
    INTEGER, INTENT(IN   ) :: ncols
    INTEGER, INTENT(IN   ) :: icount
    REAL(KIND=r8),    INTENT(INOUT) :: ft1   (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: fp1   (ncols)
    !
    INTEGER, INTENT(in   ) :: nmax
    !
    !     vegetation and soil parameters
    !
    REAL(KIND=r8),    INTENT(in   ) :: rootd (ityp,icg)
    REAL(KIND=r8),    INTENT(in   ) :: zdepth(ityp,idp)
    INTEGER, INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vcover(ncols,icg)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(in   ) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(in   ) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ea    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rst   (ncols,icg)
    REAL(KIND=r8),    INTENT(in   ) :: phsoil(ncols,idp)
    REAL(KIND=r8),    INTENT(in   ) :: stm   (ncols,icg)
    INTEGER, INTENT(in   ) :: icheck(ncols)
    REAL(KIND=r8)   , INTENT(in   ) :: topt2 (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tll2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tu2   (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: defac2(ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph12  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph22  (ncols,icg)
    !
    REAL(KIND=r8)    :: tv  (ncols)
    REAL(KIND=r8)    :: d2  (ncols)
    REAL(KIND=r8)    :: ft  (ncols)
    REAL(KIND=r8)    :: drop(ncols)
    REAL(KIND=r8)    :: fd  (ncols)
    REAL(KIND=r8)    :: fp  (ncols)
    REAL(KIND=r8)    :: ftpd(ncols)
    REAL(KIND=r8)    :: dep(3)
    REAL(KIND=r8)    :: hl3i
    REAL(KIND=r8)    :: xrot
    REAL(KIND=r8)    :: xdr
    REAL(KIND=r8)    :: arg
    INTEGER :: iveg
    INTEGER :: i
    INTEGER :: ntyp
    !
    !     humidity, temperature and transpiration factors
    !
    tv=0.0_r8   !  CALL reset(tv,ncols*13)
    hl3i=1.0_r8   /(hl*1000.0_r8  )
    iveg=1

    IF (icount == 1) THEN
       !cdir novector
       DO i = 1, nmax
          IF (icheck(i) == 1) THEN
             ntyp=itype(i)
             IF ((ntyp == 11) .OR. (ntyp == 13)) THEN
                CONTINUE
             ELSE
                IF (iveg == 1) THEN
                   tv  (i)=tc(i)
                ELSE
                   tv  (i)=tg(i)
                END IF
                tv(i)=MIN((tu2 (i,iveg)-0.1_r8   ),tv(i))
                tv(i)=MAX((tll2(i,iveg)+0.1_r8   ),tv(i))
                d2(i)=(tu2  (i,iveg)-topt2(i,iveg)) &
                     /(topt2(i,iveg)-tll2 (i,iveg))
                ft(i)=(tv(i)-tll2(i,iveg))/ &
                     (topt2(i,iveg)-tll2(i,iveg)) &
                     *EXP(d2(i)*LOG( &
                     (tu2 (i,iveg)-tv(i))/ &
                     (tu2(i,iveg)-topt2(i,iveg)) ) )
                ft(i) = MIN(ft(i), 1.e0_r8)
                ft(i) = MAX(ft(i), 1.e-5_r8)
                ft1(i) = ft(i)
                !
                !  simplified calculation of soil water potential factor, fp
                !
                xrot = rootd(ntyp,iveg)
                dep(1) = 0.0e0_r8
                dep(2) = 0.0e0_r8
                dep(3) = 0.0e0_r8
                dep(1) = MIN(zdepth(ntyp,1), xrot)
                xrot = xrot - zdepth(ntyp,1)
                IF (xrot > 0.0e0_r8) THEN
                   dep(2) = MIN(zdepth(ntyp,2), xrot)
                   xrot = xrot - zdepth(ntyp,2)
                ENDIF
                IF (xrot > 0.0e0_r8) THEN
                   dep(3) = MIN(zdepth(ntyp,3), xrot)
                   xrot = xrot - zdepth(ntyp,3)
                ENDIF
                xdr = (phsoil(i,1) * dep(1) + phsoil(i,2) * dep(2) &
                     +phsoil(i,3) * dep(3)) / rootd(ntyp,iveg)
                xdr = - xdr
                IF (xdr <= 1.0e-5_r8) xdr = 1.0e-5_r8
                xdr = LOG (xdr)
                arg = -ph12(i,1)*(ph22(i,1)-xdr)
                arg = MIN(arg,0.0_r8)
                fp(i) = 1.e0_r8 - EXP(arg)
                IF ((w(i,2) > 0.15e0_r8) .AND. (fp(i) < 0.05e0_r8)) fp(i)=0.05e0_r8
                fp(i) = MIN(fp(i), 1.e0_r8)
                fp(i) = MAX(fp(i), 1.e-5_r8)
                fp1(i) = fp(i)
             END IF
          END IF
       END DO
    END IF

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          ntyp=itype(i)
          drop(i)=EXP(21.65605_r8   -5418.0_r8   /ta(i))      -ea(i)
          fd(i) = MAX( 1.0e-5_r8,  1.0_r8/(1.0_r8+ defac2(i,iveg)*drop(i)))
          fd(i) = MIN(fd(i), 1.e0_r8)
       END IF
    END DO

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          ntyp=itype(i)
          rst(i,2) = 1.e5_r8
          IF ((ntyp == 11) .OR. (ntyp == 13)) THEN
             rst(i,1) = 1.0e5_r8
             CYCLE
          END IF
          ftpd(i)    =  fd(i)* ft1(i) * fp1(i)
          rst(i,iveg)=stm(i,iveg)/(ftpd(i)*vcover(i,iveg))
          rst(i,iveg)=MIN(rst(i,iveg),1.0e5_r8)
       END IF
    END DO

  END SUBROUTINE stres2



  ! update :performs the updating of soil moisture stores
  !         and interception capacity.


  SUBROUTINE update( &
       bps   ,deadtg,deadtc,deadqm,ect   ,eci   ,egt   ,egi   ,egs   , &
       eg    ,hc    ,hg    ,ecmass,egmass,etmass,hflux ,chf   ,shf   , &
       ra    ,rb    ,rd    ,ea    ,etc   ,etg   ,btc   ,btg   ,cc    , &
       cg    ,dtc   ,dtg   ,dtm   ,dqm   ,tc    ,tg    ,td    ,capac , &
       tm    ,nmax  ,dtc3x ,ncols)
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !   ncols.......Numero de ponto por faixa de latitude
    !   pie.........Constante Pi=3.1415926e0
    !   hl..........heat of evaporation of water   (j/kg)
    !   snomel......heat of melting
    !   tf..........Temperatura de congelamento (K)
    !   dtc3x.......time increment dt
    !   nmax........
    !   tm..........Temperature of reference (fourier)
    !   tc..........Temperatura da copa "dossel" canopy leaf temperature(K)
    !   tg..........Temperatura da superficie do solo ground temperature (K)
    !   td..........Temperatura do solo profundo (K)
    !   capac ......Agua interceptada iv=1 no dossel "water store capacity of leaves"(m)
    !   capac.......Agua interceptada iv=2 na cobertura do solo (m)
    !   ra..........Resistencia Aerodinamica (s/m)
    !   rb .........bulk boundary layer resistance             (s/m)
    !   rd..........aerodynamic resistance between ground
    !               and canopy air space
    !   ea..........Pressao de vapor
    !   etc.........Pressure of vapor at top of the copa
    !   etg.........Pressao de vapor no base da copa
    !   btc.........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i))
    !   btg.........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))
    !   cc..........heat capacity of the canopy
    !   cg..........heat capacity of the ground
    !   dtc.........dtc(i)=pblsib(i,2,5)*dtc3x
    !   dtg.........dtg(i)=pblsib(i,1,5)*dtc3x
    !   dtm.........dtm(i)=pblsib(i,3,5)*dtc3x
    !   dqm.........dqm(i)=pblsib(i,4,5)*dtc3x
    !   ect.........Transpiracao(J/m*m)
    !   eci.........Evaporacao da interceptacao da agua (J/m*m)
    !   egt ........Transpiracao na base da copa (J/m*m)  .
    !   egi.........Evaporacao da neve (J/m*m)
    !   egs.........Evaporacao do solo arido (J/m*m)
    !   eg..........Soma da transpiracao na base da copa +  Evaporacao do solo arido
    !              +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"
    !   hc..........total sensible heat lost of top from the veggies.
    !   hg..........total sensible heat lost of base from the veggies.
    !   ecmass......Mass of water lost of top from the veggies.
    !   egmass......Mass of water lost of base from the veggies.
    !   etmass......total mass of water lost from the veggies.
    !   hflux.......total sensible heat lost from the veggies.
    !   chf.........heat fluxes into the canopy  in w/m**2
    !   shf.........heat fluxes into the ground, in w/m**2
    !   deadtg......
    !   deadtc......
    !   deadqm......
    !   bps.........
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: ncols

    REAL(KIND=r8),    INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: nmax
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(in   ) :: tm  (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(in   ) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: td   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: capac(ncols,2)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(in   ) :: ra  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rb  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rd  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: etc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: etg (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: btc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: btg (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cc  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cg  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dtc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dtg (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dtm (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dqm (ncols)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8),    INTENT(inout) :: ect   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: eci   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egt   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egi   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egs   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: eg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: hc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: hg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hflux (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: chf   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: shf   (ncols)
    !
    !     derivatives
    !
    REAL(KIND=r8),    INTENT(in   ) :: deadtg(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: deadtc(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: deadqm(ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: bps  (ncols)




    REAL(KIND=r8) :: tgen  (ncols)
    REAL(KIND=r8) :: tcen  (ncols)
    REAL(KIND=r8) :: tmen  (ncols)
    REAL(KIND=r8) :: taen  (ncols)
    REAL(KIND=r8) :: eaen  (ncols)
    REAL(KIND=r8) :: d1    (ncols)
    REAL(KIND=r8) :: estarc(ncols)
    REAL(KIND=r8) :: estarg(ncols)
    REAL(KIND=r8) :: facks (ncols)
    INTEGER :: i
    REAL(KIND=r8) :: timcon
    REAL(KIND=r8) :: dtc3xi
    REAL(KIND=r8) :: hlati
    REAL(KIND=r8) :: hlat3i
    REAL(KIND=r8) :: snofac
    !
    !     adjustment of temperatures and vapor pressure ,
    !     sensible heat fluxes. n.b. latent heat fluxes cannot be derived
    !     from estarc, estarg, ea due to linear result of implicit method
    !
    !
    !
    !
    DO i = 1, nmax
       tgen(i)=tg(i)+dtg(i)
       tcen(i)=tc(i)+dtc(i)
       tmen(i)=tm(i)+dtm(i)
       d1(i)=1.0_r8   /ra(i)+1.0_r8   /rb(i)+1.0_r8   /rd(i)
       !
       !     compute the fluxes consistent with the differencing scheme.
       !
       taen(i)=(tgen(i)/rd(i)+tcen(i)/ &
            rb(i)+tmen(i)*bps(i)/ra(i))/d1(i)
       eaen(i)=ea(i)+deadtc(i)*dtc(i)+deadtg(i)* &
            dtg(i)+deadqm(i)*dqm(i)
       !
       !     vapor pressures within the canopy and the moss.
       !
       estarc(i)=etc(i)+btc(i)*dtc(i)
       estarg(i)=etg(i)+btg(i)*dtg(i)
    END DO
    DO i = 1, nmax
       IF (tgen(i) <= tf) THEN
          egs(i)=eg(i)-egi(i)
          egt(i)=0.0_r8
       END IF
    END DO
    !
    !     heat fluxes into the canopy and the ground, in w/m**2
    !
    timcon=2.0_r8   *pie/86400.0_r8
    dtc3xi=1.0_r8   /dtc3x
    hlati =1.0_r8   /        hl
    hlat3i=1.0_r8   /(1.0e3_r8*hl)
    DO i = 1, nmax
       chf(i)=dtc3xi*cc(i)*dtc(i)
       shf(i)=dtc3xi*cg(i)*dtg(i) + timcon*cg(i)*(tg(i)+dtg(i)-td(i))
    END DO
    !
    !     evaporation losses are expressed in j m-2 : when divided by
    !     ( hl*1000.0_r8) loss is in m m-2
    !
    snofac=1.0_r8   /( 1.0_r8   +snomel*hlat3i)
    DO i = 1, nmax
       facks(i)=1.0_r8
       IF (tcen(i) <= tf) facks(i)=snofac
       IF ((ect(i)+eci(i)) <= 0.0_r8) THEN
          eci(i)  =ect(i)+eci(i)
          ect(i)  =0.0_r8
          facks(i)=1.0_r8   /facks(i)
       END IF
    END DO
    DO i = 1, nmax
       capac(i,1)=capac(i,1)-eci(i)*facks(i)*hlat3i
       !
       !     mass terms are in kg m-2 dt-1
       !
       ecmass(i)=(ect(i)+eci(i)*facks(i))*hlati
    END DO
    DO i = 1, nmax
       facks(i)=1.0_r8
       IF (tgen(i) <= tf) facks(i)=snofac
       IF ((egt(i)+egi(i)) <= 0.0_r8) THEN
          egi(i)  =egt(i)+egi(i)
          egt(i)  =0.0_r8
          facks(i)=1.0_r8  /facks(i)
       END IF
    END DO
    DO i = 1, nmax
       capac(i,2)=capac(i,2)-egi(i)*facks(i)*hlat3i
       egmass(i)=(egt(i)+egs(i)+egi(i)*facks(i))*hlati
       !
       !     total mass of water and total sensible heat lost from the veggies.
       !
       etmass(i)=ecmass(i)+egmass(i)
       hflux (i)=hc(i)+hg(i)
    END DO
  END SUBROUTINE update





  SUBROUTINE sflxes(&
       hgdtg ,hgdtc ,hgdtm ,hcdtg ,hcdtc ,hcdtm ,egdtg ,egdtc ,egdqm , &
       ecdtg ,ecdtc ,ecdqm ,deadtg,deadtc,deadqm,icheck,bps   ,psb   , &
       dzm   ,em    ,gmt   ,gmq   ,cu    ,cuni  ,ctni  ,ustar ,rhoair, &
       psy   ,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    ,ect   ,eci   , &
       egt   ,egi   ,egs   ,ec    ,eg    ,hc    ,hg    ,ecidif,egidif, &
       ecmass,egmass,etmass,hflux ,chf   ,shf   ,ra    ,rb    ,rd    , &
       rc    ,rg    ,tcta  ,tgta  ,ta    ,ea    ,etc   ,etg   ,btc   , &
       btg   ,u2    ,radt  ,rst   ,rsoil ,hrr   ,phsoil,cc    ,cg    , &
       satcap,dtc   ,dtg   ,dtm   ,dqm   ,stm   ,thermk,tc    ,tg    , &
       td    ,capac ,w     ,qm    ,tm    ,um    ,vm    ,psur  ,vcover, &
       z0x   ,d     ,rdc   ,rbc   ,z0    ,itype ,dtc3x ,mon   ,nmax  , &
       jstneu,ncols ,zlt2 ,topt2  ,tll2  ,tu2   ,defac2,ph12  ,ph22  , &
       ct)

    !-----------------------------------------------------------------------
    ! sflxes :performs surface flux parameterization.
    !-----------------------------------------------------------------------
    !
    !  ncols........Numero de ponto por faixa de latitude
    !  ityp........numero das classes de solo 13
    !  imon........Numero maximo de meses no ano (12)
    !  icg.........Parametros da vegetacao (icg=1 topo e icg=2 base)
    !  idp.........Camadas de solo (1 a 3)
    !  pie.........Constante Pi=3.1415926e0
    !  stefan......Constante de Stefan Boltzmann
    !  cp..........specific heat of air (j/kg/k)
    !  hl..........heat of evaporation of water   (j/kg)
    !  grav........gravity constant      (m/s**2)
    !  snomel......heat of melting (j m-1)
    !  tf..........Temperatura de congelamento (K)
    !  gasr........Constant of dry air      (j/kg/k)
    !  epsfac......Constante 0.622 Razao entre as massas moleculares do vapor
    !              de agua e do ar seco
    !  jstneu......The first call to vntlat just gets the neutral values of ustar
    !              and ventmf para jstneu=.TRUE..
    !  dtc3x.......time increment dt
    !  mon.........Number of month at year (1-12)
    !  nmax........
    !  topt........Temperatura ideal de funcionamento estomatico
    !  tll.........Temperatura minima de funcionamento estomatico
    !  tu..........Temperatura maxima de funcionamento estomatico
    !  defac.......Parametro de deficit de pressao de vapor d'agua
    !  ph1.........Coeficiente para o efeito da agua no solo
    !  ph2.........Potencial de agua no solo para ponto de Wilting
    !  rootd.......Profundidade das raizes
    !  bee.........Expoente da curva de retencao "expoente para o solo umido"
    !  phsat.......Tensao do solo saturado " Potencial de agua no solo saturado"
    !  zdepth......Profundidade para as tres camadas de solo
    !  zlt(icg)....Indice de area foliar "LEAF AREA INDEX" icg=1 topo da copa
    !  zlt(icg)....Indice de area foliar "LEAF AREA INDEX" icg=2 base da copa
    !  x0x.........Comprimento de rugosidade
    !  xd..........Deslocamento do plano zero
    !  z2..........Altura do topo do dossel
    !  xdc.........Constant related to aerodynamic resistance
    !              between ground and canopy air space
    !  xbc.........Constant related to bulk boundary layer resistance
    !  itype.......Classe de textura do solo
    !  vcover(iv)..Fracao de cobertura de vegetacao iv=1 Top
    !  vcover(iv)..Fracao de cobertura de vegetacao iv=2 Botto
    !  z0x.........roughness length                           (m)
    !  d...........Displacement height                        (m)
    !  rdc.........Constant related to aerodynamic resistance
    !              between ground and canopy air space
    !  rbc.........Constant related to bulk boundary layer resistance
    !  z0..........Roughness length
    !  qm..........reference specific humidity (fourier)
    !  tm .........reference temperature    (fourier)                (k)
    !  um..........Razao entre zonal pseudo-wind (fourier) e seno da
    !              colatitude
    !  vm..........Razao entre meridional pseudo-wind (fourier) e seno da
    !              colatitude
    !  psur........surface pressure in mb
    !  tc..........Temperatura da copa "dossel"(K)
    !  tg..........Temperatura da superficie do solo (K)
    !  td..........Temperatura do solo profundo (K)
    !  capac(iv)...Agua interceptada iv=1 no dossel "water store capacity of leaves"(m)
    !  capac(iv)...Agua interceptada iv=2 na cobertura do solo (m)
    !  w(id).......Grau de saturacao de umidade do solo id=1 na camada superficial
    !  w(id).......Grau de saturacao de umidade do solo id=2 na camada de raizes
    !  w(id).......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    !  ra..........Resistencia Aerodinamica (s/m)
    !  rb..........bulk boundary layer resistance             (s/m)
    !  rd..........aerodynamic resistance between ground      (s/m)
    !              and canopy air space
    !  rc..........Resistencia do topo da copa
    !  rg......... Resistencia da base da copa
    !  tcta........Diferenca entre tc-ta                      (k)
    !  tgta........Diferenca entre tg-ta                      (k)
    !  ta..........Temperatura no nivel de fonte de calor do dossel (K)
    !  ea..........Pressure of vapor
    !  etc.........Pressure of vapor at top of the copa
    !  etg.........Pressao de vapor no base da copa
    !  btc.........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i)).
    !  btg.........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))
    !  u2..........wind speed at top of canopy                (m/s)
    !  radt........net heat received by canopy/ground vegetation
    !  rst ........Resisttencia Estomatica "Stomatal resistence" (s/m)
    !  rsoil ......Resistencia do solo (s/m)
    !  hrr.........rel. humidity in top layer
    !  phsoil......soil moisture potential of the i-th soil layer
    !  cc..........heat capacity of the canopy
    !  cg..........heat capacity of the ground
    !  satcap......saturation liquid water capacity         (m)
    !  dtc.........dtc(i)=pblsib(i,2,5)*dtc3x
    !  dtg.........dtg(i)=pblsib(i,1,5)*dtc3x
    !  dtm.........dtm(i)=pblsib(i,3,5)*dtc3x
    !  dqm.........dqm(i)=pblsib(i,4,5)*dtc3x
    !  stm ........Variavel utilizada mo cal. da Resistencia
    !  thermk......canopy emissivity
    !  ect.........Transpiracao no topo da copa (J/m*m)
    !  eci.........Evaporacao da agua interceptada no topo da copa (J/m*m)
    !  egt.........Transpiracao na base da copa (J/m*m)
    !  egi.........Evaporacao da neve (J/m*m)
    !  egs.........Evaporacao do solo arido (J/m*m)
    !  ec..........Soma da Transpiracao e Evaporacao da agua interceptada pelo
    !              topo da copa   ec   (i)=eci(i)+ect(i)
    !  eg..........Soma da transpiracao na base da copa +  Evaporacao do solo arido
    !              +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"
    !  hc..........total sensible heat lost of top from the veggies.
    !  hg..........total sensible heat lost of base from the veggies.
    !  ecidif......check if interception loss term has exceeded canopy storage
    !              ecidif(i)=MAX(0.0   , eci(i)-capac(i,1)*hlat3 )
    !  egidif......check if interception loss term has exceeded canopy storage
    !              ecidif(i)=MAX(0.0   , egi(i)-capac(i,1)*hlat3 )
    !  ecmass......Mass of water lost of top from the veggies.
    !  egmass......Mass of water lost of base from the veggies.
    !  etmass......total mass of water lost from the veggies.
    !  hflux.......total sensible heat lost from the veggies.
    !  chf.........heat fluxes into the canopy  in w/m**2
    !  shf.........heat fluxes into the ground, in w/m**2
    !  bps.........
    !  psb.........
    !  dzm.........Altura media de referencia  para o vento para o calculo
    !               da estabilidade do escoamento
    !  em..........Pressao de vapor da agua
    !  gmt.........
    !  gmq.........specific humidity of reference (fourier)
    !  cu..........Friction  transfer coefficients.
    !  cuni........neutral friction transfer  coefficients.
    !  ctni........neutral heat transfer coefficients.
    !  ustar.......surface friction velocity  (m/s)
    !  rhoair......Desnsidade do ar
    !  psy.........(cp/(hl*epsfac))*psur(i)
    !  rcp.........densidade do ar vezes o calor especifico do ar
    !  wc..........Minimo entre 1 e a razao entre a agua interceptada pelo
    !              indice de area foliar no topo da copa
    !  wg..........Minimo entre 1 e a razao entre a agua interceptada pelo
    !              indice de area foliar na base da copa
    !  fc..........Condicao de oravalho 0 ou 1 na topo da copa
    !  fg..........Condicao de oravalho 0 ou 1 na base da copa
    !  hr..........rel. humidity in top layer
    !  icheck......this version assumes dew-free conditions "icheck=1" to
    !              estimate ea for buoyancy term in vntmf or ra.
    !  hgdtg.......n.b. fluxes expressed in joules m-2
    !  hgdtc.......n.b. fluxes expressed in joules m-2
    !  hgdtm.......n.b. fluxes expressed in joules m-2
    !  hcdtg.......n.b. fluxes expressed in joules m-2
    !  hcdtc.......n.b. fluxes expressed in joules m-2
    !  hcdtm.......n.b. fluxes expressed in joules m-2
    !  egdtg.......partial derivative calculation for latent heat
    !  egdtc.......partial derivative calculation for latent heat
    !  egdqm.......partial derivative calculation for latent heat
    !  ecdtg.......partial derivative calculation for latent heat
    !  ecdtc.......partial derivative calculation for latent heat
    !  ecdqm.......partial derivative calculation for latent heat
    !  deadtg......
    !  deadtc......
    !  deadqm......
    !
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: ncols

    LOGICAL, INTENT(inout  ) :: jstneu

    REAL(KIND=r8)   , INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: mon(ncols)
    INTEGER, INTENT(in   ) :: nmax
    !
    INTEGER, INTENT(in   ) :: itype (ncols)
    !
    REAL(KIND=r8),    INTENT(in) :: vcover(ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: z0x   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: d     (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rdc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rbc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0    (ncols)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(inout) :: qm  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tm  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: um  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vm  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: psur(ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(inout) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: td   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: capac(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout) :: ra    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rb    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rd    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: tcta  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: tgta  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ea    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etc   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etg   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: btc   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: btg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: u2    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: radt  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: rst   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: rsoil (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hrr   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: phsoil(ncols,idp)
    REAL(KIND=r8),    INTENT(inout) :: cc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: satcap(ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: dtc   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dtg   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dtm   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dqm   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: stm   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: thermk(ncols)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8),    INTENT(inout) :: ect   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: eci   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egt   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egi   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egs   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ec    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: eg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecidif(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egidif(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hflux (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: chf   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: shf   (ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: bps   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: psb   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: em    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: gmt   (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: gmq   (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rhoair(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: psy   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rcp   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: fc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: fg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: hr    (ncols)
    INTEGER, INTENT(inout  ) :: icheck(ncols)
    !
    !     derivatives
    !
    REAL(KIND=r8),    INTENT(inout  ) :: hgdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hgdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hgdtm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hcdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hcdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hcdtm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egdqm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecdqm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: deadtg(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: deadtc(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: deadqm(ncols)
    !
    REAL(KIND=r8),    INTENT(in   ) :: zlt2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: topt2   (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tll2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tu2     (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: defac2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph12    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph22    (ncols,icg)
    !

    REAL(KIND=r8) :: ustarn(ncols)

    REAL(KIND=r8) :: psit
    REAL(KIND=r8) :: fac
    REAL(KIND=r8) :: y1
    REAL(KIND=r8) :: y2
    REAL(KIND=r8) :: ecf (ncols)
    REAL(KIND=r8) :: egf (ncols)
    REAL(KIND=r8) :: dewc(ncols)
    REAL(KIND=r8) :: dewg(ncols)
    !
    REAL(KIND=r8) :: tcsav (ncols)
    REAL(KIND=r8) :: tgsav (ncols)
    REAL(KIND=r8) :: tmsav (ncols)
    REAL(KIND=r8) :: qmsav (ncols)
    REAL(KIND=r8) :: tsav  (ncols)
    REAL(KIND=r8) :: esav  (ncols)
    REAL(KIND=r8) :: rdsav (ncols,2)
    REAL(KIND=r8) :: wt
    REAL(KIND=r8) :: ft1   (ncols)
    REAL(KIND=r8) :: fp1   (ncols)
    INTEGER :: idewco(ncols)
    !
    INTEGER, PARAMETER :: icmax = 10
    REAL(KIND=r8),    PARAMETER :: small = 1.0e-3_r8
    REAL(KIND=r8)    :: gxx
    REAL(KIND=r8)    :: capaci
    REAL(KIND=r8)    :: eee
    REAL(KIND=r8)    :: dtmdt
    REAL(KIND=r8)    :: dqmdt
    REAL(KIND=r8)    :: vcover2(ncols,icg)

    INTEGER :: i
    INTEGER :: ntyp
    INTEGER :: ncount
    INTEGER :: icount
    !
    vcover2=vcover
    DO i = 1, nmax
       tcsav(i)=tc(i)
       tgsav(i)=tg(i)
       tmsav(i)=tm(i)
       qmsav(i)=qm(i)
       rdsav(i,1)=radt(i,1)
       rdsav(i,2)=radt(i,2)
       stm(i,1)=rst(i,1)
       stm(i,2)=rst(i,2)
    END DO
    !
    !     airmod checks for the effects of snow
    !
    CALL airmod( &
         tg    ,capac ,z0x   ,d     ,rdc   ,rbc   ,itype , &
         mon   ,nmax  ,ncols )
    !
    !     sib roughness length
    !
    DO i = 1, nmax
       z0    (i)=z0x(i)
    END DO

    gxx   =grav/461.5_r8
    capaci= 1.0_r8 /0.004_r8

    DO i = 1, nmax
       ntyp   =itype(i)
       wc  (i)=MIN(1.0_r8 ,capac(i,1)/satcap(i,1))
       wg  (i)=MIN(1.0_r8 ,capac(i,2)/satcap(i,2))
       !
       !     rsoil function from fit to camillo and gurney (1985) data.
       !     wetness of upper 0.5 cm of soil calculated from approximation
       !     to milly flow equation with reduced (1/50) conductivity in
       !     top layer.
       !
       wt = MAX(small,w(i,1))
       wt=wt+(0.75_r8*zdepth(ntyp,1)/(zdepth(ntyp,1)+ &
            zdepth(ntyp,2)))*(wt-(w(i,2)*w(i,2)/wt))*0.5_r8*50.0_r8
       fac   =MIN(wt,0.99_r8)
       fac   =MAX(fac   ,small)
       rsoil(i)=101840.0_r8*(1.0_r8 - EXP(0.0027_r8 * LOG(fac   )))
       !
       !phsat =  " Potencial de agua no solo saturado"
       !
       psit = phsat(ntyp) * EXP(-bee(ntyp) * LOG(fac   ))
       !
       !        --     --
       !       |  PSI*g  |
       ! eee = |---------|
       !       |  Tg*R   |
       !        --     --
       !
       eee = psit * gxx/tg(i)
       !
       !The relative humidity of air at the soil surface
       !
       !           --     --
       !          |  PSI*g  |
       ! fh = exp*|---------|
       !          |  Tg*R   |
       !           --     --
       !
       hrr  (i)=MAX (small,EXP(eee))
       !
       hr   (i)=hrr(i)
       !
       IF (tg(i) <= tf) THEN
          vcover2(i,2)=1.0_r8
          wg    (i)  =MIN(1.0_r8 ,capac(i,2)*capaci)
          rst   (i,2)=rsoil(i)
          stm   (i,2)=rsoil(i)
       END IF
       !
       fc(i)=1.0_r8
       fg(i)=1.0_r8
    END DO
    !
    !     this is the start of iteration of time integration
    !     to avoid oscillation
    !
    ncount=0
7000 CONTINUE
    ncount=ncount+1
    DO i = 1, nmax
       icheck(i)=1
       !
       !  etc.........Pressure of vapor at top of the copa
       !  etg.........Pressao de vapor no base da copa
       !
       etc(i)=EXP(21.65605_r8  -5418.0_r8  / tc(i))
       etg(i)=EXP(21.65605_r8  -5418.0_r8  / tg(i))
    END DO
    !
    !     first guesses for ta and ea
    !
    IF (ncount == 1) THEN
       DO i = 1, nmax
          ta (i)=tc(i)
          !
          !  ea..........Pressure of vapor
          !
          ea (i)=qm(i)*psur(i)/(epsfac+qm(i))
       END DO
    END IF
    !
    !     the first call to vntlat just gets the neutral values of ustar
    !     and ventmf.
    !
    jstneu=.TRUE.

    CALL vntlax(ustarn, &
         icheck,bps   ,dzm   ,cu    ,ct    ,cuni  ,ctni  ,ustar ,ra    ,ta    , &
         u2    ,tm    ,um    ,vm    ,d     ,z0    ,itype ,z2    , &
         mon   ,nmax  ,jstneu,ncols )

    jstneu=.FALSE.

    CALL vntlax(ustarn, &
         icheck ,bps  ,dzm   ,cu    ,ct    ,cuni  ,ctni  ,ustar ,ra    ,ta    , &
         u2     ,tm   ,um    ,vm    ,d     ,z0    ,itype ,z2    , &
         mon    ,nmax ,jstneu,ncols )

    DO i = 1, nmax
       tcta(i)=tc(i)/bps(i)-tm(i)
       tgta(i)=tg(i)/bps(i)-tm(i)
    END DO
    CALL rbrd( &
         rb    ,rd    ,tcta  ,tgta  ,u2    ,tg    ,rdc   ,rbc   ,itype , &
         z2    ,mon   ,nmax  ,ncols ,zlt2)
    !
    !     iterate for air temperature and ventilation mass flux
    !     n.b. this version assumes dew-free conditions to estimate ea
    !     for buoyancy term in vntmf or ra.
    !
    icount = 0
2000 icount = icount + 1
    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          tsav(i) = ta (i)
          esav(i) = ea (i)
       END IF
    END DO
    CALL vntlax(ustarn, &
         icheck,bps   ,dzm   ,cu    ,ct    ,cuni  ,ctni  ,ustar ,ra    ,ta    , &
         u2    ,tm    ,um    ,vm    ,d     ,z0    ,itype ,z2    , &
         mon   ,nmax  ,jstneu,ncols )

    CALL cut( &
         icheck,em    ,rhoair,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    , &
         ra    ,rb    ,rd    ,rc    ,rg    ,ea    ,etc   ,etg   ,rst   , &
         rsoil ,vcover2,nmax  ,ncols )

    CALL stres2( &
         icount,ft1   ,fp1   ,icheck,ta    ,ea    ,rst   ,phsoil,stm   , &
         tc    ,tg    ,w     ,vcover2,itype ,&
         rootd ,zdepth,nmax  ,ncols ,topt2 ,tll2  ,tu2   , &
         defac2,ph12  ,ph22)

    CALL cut(  &
         icheck,em    ,rhoair,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    , &
         ra    ,rb    ,rd    ,rc    ,rg    ,ea    ,etc   ,etg   ,rst   , &
         rsoil ,vcover2,nmax  ,ncols )

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          ta(i)= (tg(i)/rd(i)+tc(i)/rb(i)+tm(i)/ra(i)*bps(i)) &
               /(1.0_r8 /rd(i)+1.0_r8 /rb(i)+1.0_r8 /ra(i))
       END IF
    END DO

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          y1   =ABS(ta(i)-tsav(i))
          y2   =ABS(ea(i)-esav(i))
          IF((y1 <= 1.0e-2_r8 .AND. y2    <= 5.0e-3_r8) &
               .OR. icount > icmax) THEN
             icheck(i)=0
          END IF
       END IF
    END DO

    DO i = 1, nmax
       IF (icheck(i) == 1) GOTO 2000
    END DO

    DO i = 1, nmax
       fc    (i)=1.0_r8
       fg    (i)=1.0_r8
       idewco(i)=0
       icheck(i)=1
    END DO

    DO i = 1, nmax
       tc(i)    =tcsav(i)
       tg(i)    =tgsav(i)
       tm(i)    =tmsav(i)
       qm(i)    =qmsav(i)
       radt(i,1)=rdsav(i,1)
       radt(i,2)=rdsav(i,2)
       etc(i)=EXP(21.65605_r8  -5418.0_r8  /tc(i))
       etg(i)=EXP(21.65605_r8  -5418.0_r8  /tg(i))
       btc(i)=EXP(30.25353_r8  -5418.0_r8  /tc(i))/(tc(i)*tc(i))
       btg(i)=EXP(30.25353_r8  -5418.0_r8  /tg(i))/(tg(i)*tg(i))
    END DO

3000 CONTINUE

    CALL cut( &
         icheck,em    ,rhoair,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    , &
         ra    ,rb    ,rd    ,rc    ,rg    ,ea    ,etc   ,etg   ,rst   , &
         rsoil ,vcover2,nmax  ,ncols )

    DO i = 1, nmax
       IF (icheck(i) == 1) THEN
          ecf (i)=SIGN(1.0_r8  ,etc(i)-ea(i))
          egf (i)=SIGN(1.0_r8  ,etg(i)-ea(i))
          dewc(i)=fc(i)*2.0_r8  -1.0_r8
          dewg(i)=fg(i)*2.0_r8  -1.0_r8
          ecf (i)=ecf(i)*dewc(i)
          egf (i)=egf(i)*dewg(i)
       END IF
    END DO

    DO i = 1, nmax
       IF ( (ecf(i) > 0.0_r8  .AND. egf(i) > 0.0_r8 ).OR. &
            idewco(i) == 3) THEN
          icheck(i)=0
       ELSE
          idewco(i)=idewco(i)+1
          IF (idewco(i) == 1) THEN
             fc(i)=0.0_r8
             fg(i)=1.0_r8
          ELSE IF (idewco(i) == 2) THEN
             fc(i)=1.0_r8
             fg(i)=0.0_r8
          ELSE IF (idewco(i) == 3) THEN
             fc(i)=0.0_r8
             fg(i)=0.0_r8
          END IF
       END IF
    END DO

    DO i=1,nmax
       IF (icheck(i) == 1) go to 3000
    END DO

    CALL temres(&
         bps   ,psb   ,em    ,gmt   ,gmq   ,psy   ,rcp   ,wc    ,wg    , &
         fc    ,fg    ,hr    ,hgdtg ,hgdtc ,hgdtm ,hcdtg ,hcdtc ,hcdtm , &
         egdtg ,egdtc ,egdqm ,ecdtg ,ecdtc ,ecdqm ,deadtg,deadtc,deadqm, &
         ect   ,eci   ,egt   ,egi   ,egs   ,ec    ,eg    ,hc    ,hg    , &
         ecidif,egidif,ra    ,rb    ,rd    ,rc    ,rg    ,ta    ,ea    , &
         etc   ,etg   ,btc   ,btg   ,radt  ,rst   ,rsoil ,hrr   ,cc    , &
         cg    ,satcap,dtc   ,dtg   ,dtm   ,dqm   ,thermk,tc    ,tg    , &
         td    ,capac ,qm    ,tm    ,psur  ,dtc3x , &
         nmax  ,vcover2,ncols)

    IF (ncount <= 1) THEN
       DO i = 1, nmax
          tc(i)=tc(i)+dtc(i)
          tg(i)=tg(i)+dtg(i)
          tm(i)=tm(i)+dtm(i)
          qm(i)=qm(i)+dqm(i)
       END DO
       go to 7000
    END IF

    CALL update( &
         bps   ,deadtg,deadtc,deadqm,ect   ,eci   ,egt   ,egi   ,egs   , &
         eg    ,hc    ,hg    ,ecmass,egmass,etmass,hflux ,chf   ,shf   , &
         ra    ,rb    ,rd    ,ea    ,etc   ,etg   ,btc   ,btg   ,cc    , &
         cg    ,dtc   ,dtg   ,dtm   ,dqm   ,tc    ,tg    ,td    ,capac , &
         tm    ,nmax  ,dtc3x ,ncols)

    DO i = 1, nmax
       fac     =grav/(100.0_r8 *psb(i)*dtc3x)
       dtmdt   =(gmt(i,3) + hflux (i) * fac   /(cp*bps(i)))/gmt(i,2)
       dqmdt   =(gmq(i,3) + etmass(i) * fac)  / gmq(i,2)
       dtm  (i)=dtmdt   *   dtc3x
       dqm  (i)=dqmdt   *   dtc3x
       gmt(i,3)=dtmdt
       gmq(i,3)=dqmdt
       tm   (i)=tm(i)+dtm(i)
       qm   (i)=qm(i)+dqm(i)
    END DO
    DO i = 1, nmax
       ntyp=itype(i)
       !vcover(i,2)=xcover(ntyp,mon(i),2)
       d     (i)=xd (ntyp,mon(i))
       z0x   (i)=x0x(ntyp,mon(i))
       rdc   (i)=xdc(ntyp,mon(i))
       rbc   (i)=xbc(ntyp,mon(i))
    END DO
  END SUBROUTINE sflxes
  !
  !
  !
  ! interc :calculation of (1) interception and drainage of rainfall and snow
  !                        (2) specific heat terms fixed for time step
  !                        (3) modifications for 4-th order model may not
  !                            conserve energy;
  !         modification: non-uniform precipitation convective ppn
  !                       is described by area-intensity relationship :-
  !
  !                       f(x)=a*exp(-b*x)+c
  !
  !                       throughfall, interception and infiltration
  !                       excess are functional on this relationship
  !                       and proportion of large-scale ppn.



  SUBROUTINE interc( &
       roff  ,cc    ,cg    ,satcap,snow  ,extk  ,tc    ,tg    ,td    , &
       capac ,w     ,tm    ,ppc   ,ppl   ,vcover,itype ,dtc3x , &
       nmax  ,ncols ,zlt2 )
    !
    !
    !        input parameters
    !-----------------------------------------------------------------------
    !   ppc.............precipitation rate ( cumulus )           (mm/s)
    !   ppl.............precipitation rate ( large scale )       (mm/s)
    !   w(1)............soil wetnessof ground surface
    !   poros...........porosity
    !   pie.............pai=3.14159..
    !   cw..............liquid water heat capacity               (j/m**3)
    !   clai............heat capacity of foliage
    !   capac(cg).......canopy/ground cover liquid water capacity(m)
    !   satcap(cg)......saturation liquid water capacity         (m)
    !   extk(cg,  ,  )..extinction coefficient
    !   zlt(1)..........canopy leaf and stem area density        (m**2/m**3)
    !   zlt(2)..........ground cover leaf and stem area index    (m**2/m**2)
    !   vcover(cg)......vegetation cover
    !   tm..............reference temperature                    (k)
    !   tc..............canopy temperature                       (k)
    !   tg..............ground temperature                       (k)
    !   tf..............freezing point                           (k)
    !   satco............mean soil hydraulic conductivity in the root zone
    !                                                            (m/s)
    !   dtc3x...........time increment dt
    !   snomel..........heat of melting                          (j/kg)
    !-----------------------------------------------------------------------
    !     in subr. parameters
    !-----------------------------------------------------------------------
    !   chisl...........soil conductivity
    !   difsl...........soil diffusivity
    !-----------------------------------------------------------------------
    !       output parameters
    !-----------------------------------------------------------------------
    !   roff............runoff
    !   snow............snow amount
    !   capac(cg).......canopy/ground cover liquid water capacity(m)
    !   cc..............heat capacity of the canopy
    !   cg..............heat capacity of the ground
    !   w(1)............soil wetnessof ground surface
    !-----------------------------------------------------------------------
    !   ncols...........Numero de ponto por faixa de latitude
    !   ityp............numero das classes de solo 13
    !   imon............Numero maximo de meses no ano (12)
    !   icg.............Parametros da vegetacao (icg=1 topo e icg=2 base)
    !   iwv.............Compriment de onda iwv=1=visivel, iwv=2=infravermelho
    !                   proximo, iwv=3 infravermelho termal
    !   idp.............Camadas de solo (1 a 3)
    !   ibd.............Estado da vegetacao  ibd=1 verde / ibd=2 seco
    !   mon.............Numero do mes do ano (1-12)
    !   nmax
    !   zdepth..........Profundidade para as tres camadas de solo
    !   itype...........Classe de textura do solo
    !   td..............Temperatura do solo profundo (K)
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: ncols

    REAL(KIND=r8),    INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: nmax

    INTEGER, INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vcover(ncols,icg)
    !
    !     the size of working area is ncols*187
    !
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(in   ) :: tm  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ppc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ppl (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(inout) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: td   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: capac(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout  ) :: cc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: cg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: satcap(ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: snow  (ncols,icg)
    REAL(KIND=r8),    INTENT(in   ) :: extk  (ncols,icg,iwv,ibd)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8),   INTENT(inout) :: roff  (ncols)
    REAL(KIND=r8)   , INTENT(in   ) :: zlt2    (ncols,icg)

    !
    REAL(KIND=r8)    :: ap         (ncols)
    REAL(KIND=r8)    :: cp         (ncols)
    REAL(KIND=r8)    :: totalp(ncols)
    REAL(KIND=r8)    :: thru  (ncols)
    REAL(KIND=r8)    :: fpi   (ncols)
    REAL(KIND=r8)    :: chisl (ncols)
    REAL(KIND=r8)    :: csoil (ncols)
    REAL(KIND=r8)    :: p0         (ncols)
    REAL(KIND=r8)    :: ts         (ncols)
    REAL(KIND=r8)    :: specht(ncols)
    REAL(KIND=r8)    :: spwet1(ncols)
    REAL(KIND=r8)    :: zload (ncols)
    REAL(KIND=r8)    :: ccp   (ncols)
    REAL(KIND=r8)    :: cct   (ncols)
    REAL(KIND=r8)    :: zmelt (ncols)
    REAL(KIND=r8)    :: xsc   (ncols)
    REAL(KIND=r8)    :: tti   (ncols)
    REAL(KIND=r8)    :: xs         (ncols)
    REAL(KIND=r8)    :: arg   (ncols)
    REAL(KIND=r8)    :: tex   (ncols)
    REAL(KIND=r8)    :: tsd   (ncols)
    REAL(KIND=r8)    :: pinf  (ncols)
    REAL(KIND=r8)    :: equdep(ncols)
    REAL(KIND=r8)    :: roffo (ncols)
    REAL(KIND=r8)    :: tsf   (ncols)
    REAL(KIND=r8)    :: diff  (ncols)
    REAL(KIND=r8)    :: freeze(ncols)
    REAL(KIND=r8)    :: ccc   (ncols)
    REAL(KIND=r8)    :: spwet (ncols)
    REAL(KIND=r8)    :: snowp (ncols,2)
    REAL(KIND=r8)    :: capacp(ncols,2)

    REAL(KIND=r8), PARAMETER :: pcoefs(2,2) = RESHAPE ( &
         (/20.0_r8    ,0.0001_r8  ,0.206e-8_r8,0.9999_r8  /), &
         (/2,2/))
    REAL(KIND=r8), PARAMETER :: bp = 20.0_r8
    REAL(KIND=r8), PARAMETER :: difsl = 5.0e-7_r8
    REAL(KIND=r8)    :: d1x
    REAL(KIND=r8)    :: theta
    INTEGER :: i
    INTEGER :: iveg
    INTEGER :: ntyp
    !
    !     diffusivity of the soil
    !            --          --
    !           |    86400.0   |
    !d1x   =SQRT|--------------|*0.5
    !           |  (pie*difsl  |
    !            --          --
    d1x   =SQRT(86400.0_r8 /(pie*difsl))*0.5_r8
    !
    ap    = 0.0_r8 !  CALL reset(ap,ncols*33)
    DO i = 1, nmax
       ap(i)=pcoefs(2,1)
       cp(i)=pcoefs(2,2)
       totalp(i) = ppc(i) + ppl(i)
       IF (totalp(i) >= 1.0e-8_r8) THEN
          !
          !       (ppc(i)*pcoefs(1,1) + ppl(i)*pcoefs(2,1))
          !ap(i)=---------------------------------------------
          !                     totalp(i)
          !
          ap(i)=(ppc(i)*pcoefs(1,1) + ppl(i)*pcoefs(2,1))/totalp(i)
          !
          !       (ppc(i)*pcoefs(1,1) + ppl(i)*pcoefs(2,1))
          !ap(i)=---------------------------------------------
          !                     totalp(i)
          !
          cp(i)=(ppc(i)*pcoefs(1,2) + ppl(i)*pcoefs(2,2))/totalp(i)
          !
       END IF
       roff(i)=0.0_r8
       thru(i)=0.0_r8
       fpi (i)=0.0_r8
       !
       !     conductivity of the soil, taking into account porosity
       !
       ntyp    = itype(i)
       !
       theta   = w(i,1)*poros(ntyp)
       !
       !            ( 9.8e-4 + 1.2e-3 * theta )
       !chisl(i) = -----------------------------
       !             ( 1.1 - 0.4 * theta )
       !
       chisl(i) = ( 9.8e-4_r8 + 1.2e-3_r8 *theta ) / ( 1.1_r8 - 0.4_r8 *theta )
       !
       chisl(i) = chisl(i)*4.186e2_r8
       !
       !     heat capacity of the soil
       !
       !            --          --
       !           |    86400.0   |
       !d1x   =SQRT|--------------|*0.5
       !           |  (pie*difsl) |
       !            --          --
       csoil(i)=chisl(i)*d1x
       !
       !     precipitation is given in mm
       !
       p0(i)=totalp(i)*0.001_r8
    END DO
    !
    !
    !
    DO iveg = 1, 2

       IF (iveg == 1) THEN
          DO i = 1, nmax
             ntyp     =itype(i)
             ts    (i)=tc (i)
             !  zlt(icg) = Indice de area foliar "LEAF AREA INDEX" icg=1 topo da copa /icg=2 base da copa
             !  clai     = heat capacity of foliage
             specht(i)=zlt2(i,1)*clai
          END DO
       ELSE
          DO i = 1, nmax
             ts    (i)=tg (i)
             specht(i)=csoil(i) !  heat capacity of the soil
          END DO
       END IF

       DO i = 1, nmax
          IF (iveg == 1 .OR. ts(i) > tf) THEN
             !
             ! capac(1/2) = canopy/ground cover liquid water capacity(m)
             ! satcap(cg) = saturation liquid water capacity         (m)
             !
             xsc(i) = MAX(0.0_r8  , capac(i,iveg) - satcap(i,iveg))
             !
             capac(i,iveg) = capac(i,iveg) - xsc(i)
             !
             roff(i) = roff(i) + xsc(i)
          END IF
       END DO

       DO i = 1, nmax
          ntyp=itype(i)
          !
          !   cw  = liquid water heat capacity (j/m**3)
          !
          spwet1(i)=MIN(0.05_r8 ,capac(i,iveg))*cw
          !
          capacp(i,iveg)=0.0_r8
          !
          snowp (i,iveg)=0.0_r8
          !
          IF (ts(i) > tf) THEN
             capacp(i,iveg)=capac (i,iveg)
          ELSE
             snowp (i,iveg)=capac (i,iveg)
          END IF
          !
          capac (i,iveg)=capacp(i,iveg)
          !
          snow  (i,iveg)=snowp (i,iveg)
          !
          zload (i)     =capac (i,iveg) + snow(i,iveg)
          !
          !                --                                            --
          !               |             --                              -- |
          !               |            | -extk(i,iveg,3,1) * zlt2(i,iveg) ||
          !fpi   (i)     =| 1.0  -  EXP|----------------------------------|| * vcover(i,iveg)
          !               |            |         vcover(i,iveg)           ||
          !               |              --                              --|
          !                --                                            --
          !
          fpi   (i)     =( 1.0_r8 -EXP(-extk(i,iveg,3,1)*zlt2(i,iveg) &
               /vcover(i,iveg))) *vcover(i,iveg)
          !
          tti(i)=p0(i)*( 1.0_r8 -fpi(i) )
          !
          IF (iveg.EQ.2) tti(i) = p0(i)
       END DO
       !
       !     proportional saturated area (xs) and leaf drainage(tex)
       !
       DO i = 1, nmax
          xs(i)=1.0_r8
          IF (p0(i) >= 1.0e-9_r8) THEN
             !
             !        (satcap(i,iveg) - zload(i))      cp(i)
             !arg(i)=----------------------------- - ---------
             !          (p0(i)*fpi(i)*ap(i))           ap(i)
             !
             arg(i)=(satcap(i,iveg)-zload(i))/ &
                  (p0(i)*fpi(i)*ap(i)) - cp(i)/ap(i)
             IF (arg(i) >= 1.0e-9_r8) THEN
                !
                !         -1.0
                !xs(i) = ------ * LOG(arg(i))
                !          bp
                !
                xs(i)=-1.0_r8/bp * LOG( arg(i) )
                xs(i)= MIN ( xs(i) , 1.0_r8 )
                xs(i)= MAX ( xs(i) , 0.0_r8 )
             END IF
          END IF
       END DO

       DO i = 1, nmax
          !                     --                                        --
          !                    | ap(i)                                      |
          !tex(i)=p0(i)*fpi(i)*|-------*(1.0 - EXP(-bp*xs(i))) + cp(i)*xs(i)|-(satcap(i,iveg) - zload(i))*xs(i)
          !                    |  bp                                        |
          !                     --                                        --
          tex(i)=p0(i)*fpi(i)*(ap(i)/bp*(1.0_r8 -EXP(-bp*xs(i)))+cp(i)*xs(i)) &
               -(satcap(i,iveg)-zload(i))*xs(i)

          tex(i)= MAX ( tex(i), 0.0_r8 )
          !
          IF (iveg == 2) tex(i) = 0.0_r8
          !
          !     total throughfall (thru) and store augmentation
          !
          thru(i)=tti(i)+tex(i)
          IF (iveg == 2 .AND. tg(i) <= tf) THEN
             thru(i)=0.0_r8
          END IF

          pinf(i)=p0(i) - thru(i)

          IF (tm(i) > tf) THEN
             capac(i,iveg) = capac(i,iveg) + pinf(i)
          ELSE
             snow (i,iveg) = snow (i,iveg) + pinf(i)
          END IF
       END DO

       IF (iveg == 2)   THEN
          DO i = 1, nmax
             ntyp=itype(i)
             IF (tm(i) <= tf) THEN
                snow  (i,iveg) = snowp(i,iveg) + p0(i)
                thru  (i)=0.0_r8
             ELSE
                !
                !     instantaneous overland flow contribution ( roff )
                !
                equdep(i)=satco(ntyp)*dtc3x
                xs(i)=1.0_r8
                IF (thru(i) >= 1.0e-9_r8) THEN
                   arg(i)=equdep(i)/( thru(i)*ap(i) ) -cp(i)/ap(i)
                   IF (arg(i) >= 1.0e-9_r8) THEN
                      xs(i)=-1.0_r8 /bp* LOG( arg(i) )
                      xs(i)= MIN ( xs(i), 1.0_r8 )
                      xs(i)= MAX ( xs(i), 0.0_r8 )
                   END IF
                END IF
                roffo(i)=thru(i)* &
                     (ap(i)/bp*(1.0_r8 -EXP(-bp*xs(i)))+cp(i)*xs(i)) &
                     -equdep(i)*xs(i)
                roffo(i)= MAX ( roffo(i), 0.0_r8 )
                roff (i)= roff (i)+roffo(i)
                w(i,1)=w(i,1)+(thru(i)-roffo(i))/ &
                     ( poros(ntyp)*zdepth(ntyp,1))
             END IF
          END DO
       END IF
       !
       !     temperature change due to addition of precipitation
       !
       DO i = 1, nmax
          diff(i)=(capac (i,iveg)+snow (i,iveg) &
               -capacp(i,iveg)-snowp(i,iveg))*cw
          ccp(i)=specht(i)+spwet1(i)
          cct(i)=specht(i)+spwet1(i)+diff(i)
          tsd(i)=( ts(i)*ccp(i)+tm(i)*diff(i) )/cct(i)
          tsf(i)=( ts(i)-tf)*( tm(i)-tf)
       END DO
       DO i = 1, nmax
          IF (tsf(i) < 0.0_r8) THEN
             IF (tsd(i) <= tf) THEN
                !
                !     freezing of water on canopy or ground
                !
                ccc(i)=capacp(i,iveg)*snomel
                IF (ts(i) < tm(i)) ccc(i)=diff(i)*snomel/cw
                tsd   (i)=( ts(i)*ccp(i)+tm(i)*diff(i)+ccc(i) )/cct(i)
                freeze(i)= tf*cct(i)-( ts(i)*ccp(i)+tm(i)*diff(i) )
                freeze(i)=( MIN ( ccc(i), freeze(i) ))/snomel
                IF (tsd(i) > tf) tsd(i) = tf - 0.1_r8
                snow (i,iveg)=snow (i,iveg)+freeze(i)
                capac(i,iveg)=capac(i,iveg)-freeze(i)
             ELSE
                !
                !     melting of water on canopy or ground
                !
                ccc(i)=- snow(i,iveg)*snomel
                IF (ts(i) > tm(i)) ccc(i)=- diff(i)*snomel/cw
                tsd   (i)=( ts(i)*ccp(i)+tm(i)*diff(i)+ccc(i) )/cct(i)
                freeze(i)=( tf*cct(i)-( ts(i)*ccp(i)+tm(i)*diff(i) ))
                freeze(i)= MAX ( ccc(i), freeze(i) ) /snomel
                IF (tsd(i) <= tf) tsd(i) = tf - 0.1_r8
                snow (i,iveg)=snow (i,iveg)+freeze(i)
                capac(i,iveg)=capac(i,iveg)-freeze(i)
             END IF
          END IF
       END DO
       DO i = 1, nmax
          IF (iveg == 1) THEN
             tc(i)=tsd(i)
          ELSE
             tg(i)=tsd(i)
          END IF
       END DO
       DO i = 1, nmax
          IF (snow(i,iveg) >= 0.0000001_r8 .OR. iveg == 2) THEN
             zmelt(i) = 0.0_r8
             IF (td(i) > tf) THEN
                zmelt(i)=capac(i,iveg)
             ELSE
                roff (i)=roff(i)+capac(i,iveg)
             END IF
             capac(i,iveg)=0.0_r8
             !
             !     if tg is less than tf water accumulates as snowpack in capac(2)
             !
             ntyp=itype(i)
             w(i,1)=w(i,1)+zmelt(i)/( poros(ntyp)*zdepth(ntyp,1))
          END IF
       END DO
       DO i = 1, nmax
          !
          !     these lines exist to eliminate a cray compiler error
          !
          IF (iveg == 2) THEN
             IF (snow(i,2) > 0.0_r8 .AND. tg(i) > 273.16_r8) THEN
             END IF
             IF (capac(i,2) > 0.0_r8 .AND. tg(i) > 273.16_r8) THEN
             END IF
          END IF
          capac(i,iveg)=capac(i,iveg)+snow(i,iveg)
          snow (i,iveg)=0.0_r8
          p0(i)=thru(i)
       END DO
    END DO
    !
    !     calculation of canopy and ground heat capacities.
    !
    DO i = 1, nmax
       ntyp=itype(i)
       cc(i)=zlt2(i,1)*clai+capac(i,1)*cw
       spwet(i)=MIN( 0.05_r8 , capac(i,2))*cw
       cg(i)=csoil(i)+spwet(i)
    END DO
  END SUBROUTINE interc



  ! stomat :performs stomatal resistance.



  SUBROUTINE stomat( &
       cosz  ,par   ,pd    ,rst   ,extk  ,vcover,itype ,nmax  ,ncols ,&
       zlt2  ,green2,chil2 ,rstpar2)
    !
    !
    !-----------------------------------------------------------------------
    !      input parameters
    !-----------------------------------------------------------------------
    !   cosz.............cosine of zenith angle
    !   extk(cg,vnt,bd)..extinction coefficient
    !   zlt   (cg).......leaf area index
    !   vcover(cg).......fraction of vegetation cover
    !   green (cg).......fraction of grenn leaves
    !   chil  (cg).......leaf orientation pameter
    !   rstpar(cg,3).....coefficints related to par influence on
    !                    stomatal resistance
    !   radn   (vnt,bd)..downward sw/lw radiation at the surface
    !   par   (cg).......par( photo-synthetic active radiation)
    !   pd    (cg).......ratio of par(beam) to par(beam+diffuse)
    !-----------------------------------------------------------------------
    !     output parameters
    !-----------------------------------------------------------------------
    !   rst(cg)..........stomatal reistance
    !-----------------------------------------------------------------------
    !   itype............Classe de textura do solo
    !   nmax
    !   pie..............Constante Pi=3.1415926e0
    !   athird...........Constante athird=1.0e0 /3.0e0
    !   ncols............Numero de ponto por faixa de latitude
    !   ityp.............numero das classes de solo 13
    !   imon.............Numero maximo de meses no ano (12)
    !   icg..............Parametros da vegetacao (icg=1 topo e icg=2 base)
    !   iwv..............Compriment de onda iwv=1=visivel, iwv=2=infravermelho
    !                    proximo, iwv=3 infravermelho termal
    !   ibd..............Estado da vegetacao ibd=1 verde / ibd=2 seco
    !-----------------------------------------------------------------------
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: nmax
    REAL(KIND=r8)   , INTENT(in   ) :: rstpar2(ncols,icg,iwv)
    INTEGER, INTENT(in   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vcover(ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: zlt2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: green2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: chil2   (ncols,icg)

    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(in   ) :: par   (ncols,icg)
    REAL(KIND=r8),    INTENT(in   ) :: pd    (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: rst   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: extk  (ncols,icg,iwv,ibd)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: cosz  (ncols)
    !
    REAL(KIND=r8)    :: f     (ncols)
    REAL(KIND=r8)    :: gamma (ncols)
    REAL(KIND=r8)    :: at    (ncols)
    REAL(KIND=r8)    :: power1(ncols)
    REAL(KIND=r8)    :: power2(ncols)
    REAL(KIND=r8)    :: aa    (ncols)
    REAL(KIND=r8)    :: bb    (ncols)
    REAL(KIND=r8)    :: zat   (ncols)
    REAL(KIND=r8)    :: zk    (ncols)
    REAL(KIND=r8)    :: ekat  (ncols)
    REAL(KIND=r8)    :: rho4  (ncols)
    REAL(KIND=r8)    :: avflux(ncols)
    !
    INTEGER :: i
    INTEGER :: iveg
    INTEGER :: irad
    REAL(KIND=r8)    :: fcon
    REAL(KIND=r8)    :: xabc
    REAL(KIND=r8)    :: xabd
    REAL(KIND=r8)    :: ftemp
    !
    !
    !     bounding of product of extinction coefficient and local l.a.i.
    !
    DO i = 1, nmax
       f(i) = MAX( cosz(i), 0.01746_r8 )
    END DO
    !
    DO iveg = 1, 2
       DO irad = 1, 2 !Estado da vegetacao irad=1 verde / irad=2 seco
          DO i = 1, nmax
             !             --                   --
             !            |             150       |
             ! extk = MIN |  extk, -------------- |
             !            |         zlt2 * vcover |
             !             --                   --
             !
             extk(i,iveg,1,irad)=min(extk(i,iveg,1,irad),150.0_r8 / &
                  zlt2(i,iveg)*vcover(i,iveg))
          END DO
       END DO
    END DO
    !
    fcon  =0.25_r8*pie + athird
    iveg=1
    !
    DO i = 1, nmax
       IF (itype(i) == 13 .OR. itype(i) == 11) THEN
          rst(i,iveg) = 1.0e5_r8
       ELSE
          !
          !        zlt2           leaf area index
          ! at = -------- = ------------------------------
          !       vcover     fraction of vegetation cover
          !
          !
          at(i) = zlt2(i,iveg)/vcover(i,iveg)
          !
          IF (par(i,iveg) <= 0.00101_r8) THEN
             !
             ! iwv........Compriment de onda iwv=1=visivel, iwv=2=infravermelho
             !             proximo, iwv=3 infravermelho termal
             !
             !            rstpar(visivel)
             ! xabc = ------------------------- + rstpar(infravermelho termal)
             !          rstpar(infravermelho)
             !
             xabc = rstpar2(i,iveg,1) / rstpar2(i,iveg,2) + rstpar2(i,iveg,3)
             !
             !        0.5
             ! xabd =------ *  at(i)
             !        xabc
             !
             xabd = 0.5_r8  / xabc * at(i)
             !
             !          1
             ! rst  = ------
             !         xabd
             !
             rst(i,iveg) = 1.0_r8 / xabd
          ELSE
             !
             !         (rstpar2(visivel) + rstpar2(infravermelho)* rstpar2(infravermelho termal))
             ! gamma =---------------------------------------------------------------------------
             !                          rstpar2(infravermelho termal)
             !
             gamma(i)  = (rstpar2(i,iveg,1) + rstpar2(i,iveg,2) &
                  * rstpar2(i,iveg,3))/ rstpar2(i,iveg,3)
             !
             !     single extinction coefficient using weighted
             !     values of direct and diffus contributions to p.a.r.
             !
             !
             !        zlt            leaf area index
             ! at = -------- = ------------------------------
             !       vcover     fraction of vegetation cover
             !
             !
             at(i)     = zlt2(i,iveg)/vcover(i,iveg)
             !
             !           zlt          150
             !power1 = -------- * --------------
             !          vcover     zlt2 * vcover
             !
             power1(i) = at(i)*extk(i,iveg,1,1)!Estado da vegetacao irad=1 verde
             power2(i) = at(i)*extk(i,iveg,1,2)!Estado da vegetacao irad=2 seco
             !
             ! chil2   Leaf orientation parameter
             !  icg    Parameters of vagetation (icg=1 top e icg=2 bottom)
             !
             ! aa(i)   = 0.5 - (0.633 + 0.33 * chil2(i,icg)) * chil2(i,icg)
             !
             aa(i)     = 0.5_r8 -(0.633_r8 + 0.33_r8 * chil2(i,iveg)) * chil2(i,iveg)
             !
             bb(i)     = 0.877_r8 -1.754_r8 *aa(i)
             !
             !
             !        LOG(( EXP(-power1(i)) + 1 ) * 0.5 ) * pd(i,iveg)
             !zat = ------------------------------------------------------
             !                       extk(i,iveg,1,1)
             !
             !
             zat(i)    = LOG(( EXP(-power1(i))+1.0_r8 )*0.5_r8 ) * pd(i,iveg) / extk(i,iveg,1,1)
             !
             zat(i)    = zat(i) + LOG((EXP(-power2(i)) + 1.0_r8 )*0.5_r8 )*( 1.0_r8 -pd(i,iveg))/extk(i,iveg,1,2)
             !
             zk(i)     = 1.0_r8 /zat(i) * LOG(pd(i,iveg) *EXP( power1(i)*zat(i)/at(i) ) &
                  + (1.0_r8 -pd(i,iveg))*EXP( power2(i)*zat(i)/at(i) ))
             !
             !     canopy and ground cover bulk resistances using
             !     ross-goudriaan leaf function , total par flux (avflux) and
             !     mean extinction coefficient (zk)
             !
             ftemp       = MIN( zk(i)*at(i),20.0_r8 )
             ekat (i)    = EXP( ftemp )
             !
             !                                    --          --
             !                                   |  aa(i)       |
             !avflux = par(i,iveg) * (pd(i,iveg)*|------ + bb(i)| + (1 - pd(i,iveg))*(bb(i) * fcon + aa(i)*1.5))
             !                                   |  f(i)        |
             !                                    --          --
             !
             avflux(i)   = par(i,iveg)*( pd(i,iveg)*( aa(i)/f(i)+bb(i)) &
                  + ( 1.0_r8 -pd(i,iveg))*( bb(i)*fcon+aa(i)*1.5_r8 ))
             !
             !                 gamma(i)
             !rho4(i)     = ----------------
             !                 avflux(i)
             !
             rho4(i)     = gamma(i)/avflux(i)
             !
             !                         rstpar2(i,iveg,2)
             !rst(i,iveg) = ----------------------------------------------
             !                            --                        --
             !                           |  (rho4(i) * ekat(i) + 1.0) |
             !               gamma(i)*LOG|----------------------------|
             !                           |   (rho4(i) + 1.0 )         |
             !                            --                        --
             !
             rst(i,iveg) = rstpar2(i,iveg,2) / gamma(i)*LOG((rho4(i)*ekat(i)+1.0_r8 )/(rho4(i)+1.0_r8 ))
             !
             !
             !
             !
             !
             rst(i,iveg)=rst(i,iveg) - LOG((rho4(i)+1.0_r8 /ekat(i))/(rho4(i)+1.0_r8 ))
             !
             !             rst(i,iveg)
             ! rst =----------------------------
             !         zk(i) * rstpar2(i,iveg,3)
             !
             !
             rst(i,iveg)=rst(i,iveg)/(zk(i)*rstpar2(i,iveg,3))
             !
             !                      1
             ! rst = --------------------------------
             !          rst(i,iveg) * green2(i,iveg)
             !
             rst(i,iveg)=1.0_r8 /( rst(i,iveg)*green2(i,iveg))
          END IF
       END IF
    END DO
    !
    DO i = 1, nmax
       rst(i,2) = 1.0e5_r8
    END DO
    !
  END SUBROUTINE stomat



  ! raduse :performs the absorption of radiation by surface.



  SUBROUTINE raduse(radt  ,par   ,pd    ,radfac,closs ,gloss ,thermk,p1f   , &
       p2f   ,radn  ,vcover,nmax  ,ncols )
    !
    !-----------------------------------------------------------------------
    ! input parameters
    !-----------------------------------------------------------------------
    !   tf...............freezing temperature
    !   tg...............ground   temperature
    !   polar............
    !   radsav...........passesd from subr.radalb
    !   radfac(cg,vn,bd).fractions of downward solar radiation at surface
    !                    passed from subr.radalb
    !   radn(vnt,bd).....downward sw/lw radiation at the surface
    !   vcover(cg).......vegetation cover
    !-----------------------------------------------------------------------
    ! output parameters
    !-----------------------------------------------------------------------
    !   radt(cg).........net heat received by canopy/ground vegetation
    !                    by radiation & conduction
    !   par(cg)..........par incident on canopy
    !   pd(cg)...........ratio of par beam to total par
    !-----------------------------------------------------------------------
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: nmax
    REAL(KIND=r8),    INTENT(in   ) :: vcover(ncols,icg)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(in   ) :: radn  (ncols,3,2)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout  ) :: radt  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: par   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: pd    (ncols,icg)
    REAL(KIND=r8),    INTENT(in   ) :: radfac(ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(in   ) :: closs (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gloss (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: thermk(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: p1f   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: p2f   (ncols)


    REAL(KIND=r8) :: p1 (ncols)
    REAL(KIND=r8) :: p2 (ncols)


    INTEGER :: i
    INTEGER :: iveg
    INTEGER :: iwave
    INTEGER :: irad
    DO i = 1, nmax
       radt(i,1)=0.0_r8
       radt(i,2)=0.0_r8
    END DO
    !
    ! radn(1,1)=!Downward Surface shortwave fluxe visible beam (cloudy)
    ! radn(1,2)=!Downward Surface shortwave fluxe visible diffuse (cloudy)
    ! radn(2,1)=!Downward Surface shortwave fluxe Near-IR beam (cloudy)
    ! radn(2,2)=!Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
    ! radfac(cg,vn,bd).fractions of downward solar radiation at surface
    !                    passed from subr.radalb
    !
    !     summation of radiation fractions for canopy and ground
    !
    DO iveg = 1, 2
       DO iwave = 1, 2
          DO irad = 1, 2
             DO i = 1, nmax
                radt(i,iveg)=radt(i,iveg)+radfac(i,iveg,iwave,irad)*radn(i,iwave,irad)
             END DO
          END DO
       END DO
    END DO
    !
    !     total long wave ( and polar ice conduction ) adjustments to
    !     canopy and ground net radiation terms
    !     thermk = canopy emissivity
    DO i = 1, nmax
       radt(i,1) = radt(i,1) + radn(i,3,2) * vcover(i,1)*(1.0_r8 -thermk(i)) - closs(i)
       radt(i,2) = radt(i,2) + radn(i,3,2) * (1.0_r8 - vcover(i,1)*(1.0_r8 -thermk(i))) - gloss(i)
       par(i,1)  = radn(i,1,1) + radn(i,1,2) + 0.001_r8! total par incident on canopy
       pd (i,1)  = (radn(i,1,1) + 0.001_r8 ) / par(i,1)! ratio of par beam on topo of the canopy to total par
       p1(i)     = p1f(i)*radn(i,1,1) + 0.001_r8 ! net par beam on topo of the canopy
       p2(i)     = p2f(i)*radn(i,1,2)            ! net par beam on base of the canopy
       par(i,2)  = p1(i)+p2(i)! net par incident on canopy and ground
       IF (par(i,1) <= 0.000001_r8) par(i,1) = 0.000001_r8
       IF (par(i,2) <= 0.000001_r8) par(i,2) = 0.000001_r8
       pd (i,2)  = p1(i)/par(i,2) !ratio of net par beam to net par incident on canopy and ground
    END DO
  END SUBROUTINE raduse



  ! root   :performs soil moisture potentials in root zone of each
  !         vegetation layer and summed soil+root resistance.



  SUBROUTINE root(phroot,phsoil,w     ,itype ,nmax  , ncols )
    !
    ! input parameters
    !-----------------------------------------------------------------------
    !   w(1).............wetness of surface store
    !   w(2).............wetness of root zone
    !   w(3).............wetness of recharge zone
    !   phsat............soil moisture potential at saturation   (m)
    !   bee..............empirical constant
    !   zdepth(3)........depth of the i-th soil layer            (m)
    !   rootd (cg).......rooting depth                           (m)
    !   satco............mean soil hydraulic conductivity in the root zone
    !                                                            (m/s)
    !   rootl(cg)........root density                            (m/m**3)
    !   rootca(cg).......root cross section                      (m**2)
    !   rdres(cg)........resistance per unit root length         (s/m)
    !   rplant(cg).......area averaged resistance imposed by the plant
    !                    vascular system                         (s)
    !-----------------------------------------------------------------------
    ! output parameters
    !-----------------------------------------------------------------------
    !   vroot............root volume density                     (m**3/m**3)
    !-----------------------------------------------------------------------
    ! output parameters
    !-----------------------------------------------------------------------
    !   phsoil(3)........soil moisture potential of the i-th soil layer
    !                                                            (m)
    !   rootr(cg)........root resistance                         (s)
    !-----------------------------------------------------------------------
    !
    !   imax.............Numero de ponto por faixa de latitude
    !   ityp.............numero das classes de solo 13
    !   icg..............Parametros da vegetacao (icg=1 topo e icg=2 base)
    !   idp..............Camadas de solo (1 a 3)
    !   nmax.............
    !   itype............Classe de textura do solo
    !   phroot...........Soil moisture potentials in root zone of each
    !                    vegetation layer and summed soil+root resistance.
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: nmax

    INTEGER, INTENT(in   ) :: itype (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(in   ) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout  ) :: phroot(ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: phsoil(ncols,idp)


    REAL(KIND=r8)     :: www   (ncols,3)
    INTEGER  :: i
    INTEGER  :: n

    DO i = 1, 3
       DO n = 1, nmax
          !      0
          !w = -----
          !      0s
          !
          www   (n,i) = MAX(0.10_r8 ,w(n,i))
          !                            --          --
          !                           |         0    |
          ! phsoil(n,i) = phsat * EXP |-b*LOG(-----) |
          !                           |         0s   |
          !                            --          --
          phsoil(n,i) = phsat(itype(n)) * EXP(-bee( itype(n))*LOG(www(n,i)))
          !
       END DO
    END DO
    !
    !   Soil moisture potentials in root zone of each
    !   vegetation layer and summed soil+root resistance.
    !
    DO n = 1, nmax
       phroot(n,1)   = phsoil(n,1) - 0.01_r8
       DO i = 2, 3
          phroot(n,1) = MAX( phroot(n,1), phsoil(n,i))
       END DO
       phroot(n,2)   = phroot(n,1)
    END DO

  END SUBROUTINE root



  ! pbl    :performs planetary boundary layer parameterization.



  SUBROUTINE pbl(jstneu, hgdtg , hgdtc , hgdtm , hcdtg , hcdtc , hcdtm , &
       egdtg , egdtc , egdqm , ecdtg , ecdtc , ecdqm , deadtg, &
       deadtc, deadqm,icheck , ect   , eci   , egt   , egi   , &
       egs   , ec    , eg    , hc    , hg    , ecidif, egidif, &
       ecmass, egmass, etmass, hflux , chf   , shf   , roff  , &
       bps   , psb   , dzm   , em    , gmt   , gmq   ,  cu   , &
       cuni  , ctni  , ustar , cosz  , rhoair, psy   , rcp   , &
       wc    , wg    , fc    , fg    , hr    , vcover, z0x   , &
       d     , rdc   , rbc   , z0    , qm    , tm    , um    , &
       vm    , psur  , ppc   , ppl   , radn  , ra    , rb    , &
       rd    , rc    , rg    , tcta  , tgta  , ta    , ea    , &
       etc   , etg   , btc   , btg   , u2    , radt  , par   , &
       pd    , rst   , rsoil , phroot,  hrr  , phsoil, cc    , &
       cg    , satcap, snow  , dtc   , dtg   , dtm   , dqm   , &
       stm   , extk  , radfac, closs , gloss , thermk, p1f   , &
       p2f   , tc    , tg    , td    , capac , w     , itype , &
       dtc3x , mon   , nmax  , ncols ,zlt2  ,green2,chil2,rstpar2,&
       topt2,tll2    ,tu2    , defac2,ph12  ,ph22    ,ct )
    !
    ! jstneu......The first call to vntlat just gets the neutral values of ustar
    !              and ventmf para jstneu=.TRUE..
    ! hgdtg.......n.b. fluxes expressed in joules m-2
    ! hgdtc.......n.b. fluxes expressed in joules m-2
    ! hgdtm.......n.b. fluxes expressed in joules m-2
    ! hcdtg.......n.b. fluxes expressed in joules m-2
    ! hcdtc.......n.b. fluxes expressed in joules m-2
    ! hcdtm.......n.b. fluxes expressed in joules m-2
    ! egdtg.......partial derivative calculation for latent heat
    ! egdtc.......partial derivative calculation for latent heat
    ! egdqm.......partial derivative calculation for latent heat
    ! ecdtg.......partial derivative calculation for latent heat
    ! ecdtc.......partial derivative calculation for latent heat
    ! ecdqm.......partial derivative calculation for latent heat
    ! deadtg
    ! deadtc
    ! deadqm
    ! icheck......this version assumes dew-free conditions "icheck=1" to
    !              estimate ea for buoyancy term in vntmf or ra.
    ! ect.........Transpiracao no topo da copa (J/m*m)
    ! eci.........Evaporacao da agua interceptada no topo da copa (J/m*m)
    ! egt.........Transpiracao na base da copa (J/m*m)
    ! egi.........Evaporacao da neve (J/m*m)
    ! egs.........Evaporacao do solo arido (J/m*m)
    ! ec..........Soma da Transpiracao e Evaporacao da agua interceptada pelo
    !              topo da copa   ec   (i)=eci(i)+ect(i)
    ! eg..........Soma da transpiracao na base da copa +  Evaporacao do solo arido
    !              +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"
    ! hc..........total sensible heat lost of top from the veggies.
    ! hg..........total sensible heat lost of base from the veggies.
    ! ecidif......check if interception loss term has exceeded canopy storage
    !              ecidif(i)=MAX(0.0   , eci(i)-capac(i,1)*hlat3 )
    ! egidif......check if interception loss term has exceeded canopy storage
    !              ecidif(i)=MAX(0.0   , egi(i)-capac(i,1)*hlat3 )
    ! ecmass......Mass of water lost of top from the veggies.
    ! egmass......Mass of water lost of base from the veggies.
    ! etmass......total mass of water lost from the veggies.
    ! hflux.......total sensible heat lost from the veggies.
    ! chf.........heat fluxes into the canopy  in w/m**2
    ! shf.........heat fluxes into the ground, in w/m**2
    ! roff........runoff
    ! pie.........Constante Pi=3.1415926e0
    ! stefan......Constante de Stefan Boltzmann
    ! cpair.......specific heat of air (j/kg/k)
    ! hlat........heat of evaporation of water   (j/kg)
    ! grav........gravity constant      (m/s**2)
    ! snomel......heat of melting (j m-1)
    ! tf..........Temperatura de congelamento (K)
    ! clai........heat capacity of foliage
    ! cw..........liquid water heat capacity               (j/m**3)
    ! gasr........Constant of dry air      (j/kg/k)
    ! epsfac......Constante 0.622 Razao entre as massas moleculares do vapor
    !              de agua e do ar seco
    ! athird......Constante athird=1.0e0 /3.0e0
    ! bps
    ! psb
    ! dzm.........Altura media de referencia  para o vento para o calculo
    !               da estabilidade do escoamento
    ! em..........Pressao de vapor da agua
    ! gmt(i,k,3)..virtual temperature tendency due to vertical diffusion
    ! gmq.........specific humidity of reference (fourier)
    ! cu..........Friction  transfer coefficients.
    ! cuni........neutral friction transfer  coefficients.
    ! ctni........neutral heat transfer coefficients.
    ! ustar.......surface friction velocity  (m/s)
    ! cosz........cosine of zenith angle
    ! rhoair......Desnsidade do ar
    ! psy ........(cp/(hl*epsfac))*psur(i)
    ! rcp.........densidade do ar vezes o calor especifico do ar
    ! wc..........Minimo entre 1 e a razao entre a agua interceptada pelo
    !              indice de area foliar no topo da copa
    ! wg..........Minimo entre 1 e a razao entre a agua interceptada pelo
    !              indice de area foliar na base da copa
    ! fc..........Condicao de oravalho 0 ou 1 na topo da copa
    ! fg..........Condicao de oravalho 0 ou 1 na base da copa
    ! hr..........rel. humidity in top layer
    ! vcover(iv)..Fracao de cobertura de vegetacao iv=1 Top
    ! vcover(iv)..Fracao de cobertura de vegetacao iv=2 Botto
    ! z0x.........roughness length
    ! d...........Displacement height
    ! rdc.........Constant related to aerodynamic resistance
    !              between ground and canopy air space
    ! rbc.........Constant related to bulk boundary layer resistance
    ! z0..........Roughness length
    ! qm..........reference specific humidity (fourier)
    ! tm .........reference temperature    (fourier)                (k)
    ! um..........Razao entre zonal pseudo-wind (fourier) e seno da
    !              colatitude
    ! vm..........Razao entre meridional pseudo-wind (fourier) e seno da
    !              colatitude
    ! psur........surface pressure in mb
    ! ppc.........precipitation rate ( cumulus )           (mm/s)
    ! ppl.........precipitation rate ( large scale )       (mm/s)
    ! radn........downward sw/lw radiation at the surface
    ! ra..........Resistencia Aerodinamica (s/m)
    ! rb..........bulk boundary layer resistance             (s/m)
    ! rd..........aerodynamic resistance between ground      (s/m)
    !              and canopy air space
    ! rc..........Resistencia do topo da copa
    ! rg......... Resistencia da base da copa
    ! tcta........Diferenca entre tc-ta                      (k)
    ! tgta........Diferenca entre tg-ta                      (k)
    ! ta..........Temperatura no nivel de fonte de calor do dossel (K)
    ! ea..........Pressure of vapor
    ! etc.........Pressure of vapor at top of the copa
    ! etg.........Pressao de vapor no base da copa
    ! btc.........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i)).
    ! btg.........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))
    ! u2..........wind speed at top of canopy                (m/s)
    ! radt........net heat received by canopy/ground vegetation
    ! par.........par incident on canopy
    ! pd..........ratio of par beam to total par
    ! rst ........Resisttencia Estomatica "Stomatal resistence" (s/m)
    ! rsoil ......Resistencia do solo (s/m)
    ! phroot......Soil moisture potentials in root zone of each
    !                    vegetation layer and summed soil+root resistance.
    ! hrr.........rel. humidity in top layer
    ! phsoil......soil moisture potential of the i-th soil layer
    ! cc..........heat capacity of the canopy
    ! cg..........heat capacity of the ground
    ! satcap......saturation liquid water capacity         (m)
    ! snow........snow amount
    ! dtc.........dtc(i)=pblsib(i,2,5)*dtc3x
    ! dtg.........dtg(i)=pblsib(i,1,5)*dtc3x
    ! dtm.........dtm(i)=pblsib(i,3,5)*dtc3x
    ! dqm ........dqm(i)=pblsib(i,4,5)*dtc3x
    ! stm ........Variavel utilizada mo cal. da Resisttencia
    ! extk........extinction coefficient
    ! radfac......fractions of downward solar radiation at surface
    !             passed from subr.radalb
    ! closs.......radiation loss from canopy
    ! gloss.......radiation loss from ground
    ! thermk......canopy emissivity
    ! p1f
    ! p2f
    ! tc..........Temperatura da copa "dossel"(K)
    ! tg..........Temperatura da superficie do solo (K)
    ! td..........Temperatura do solo profundo (K)
    ! capac(iv)...Agua interceptada iv=1 no dossel "water store capacity
    !             of leaves"(m)
    ! capac(iv)...Agua interceptada iv=2 na cobertura do solo (m)
    ! w(id).......Grau de saturacao de umidade do solo id=1 na camada superficial
    ! w(id).......Grau de saturacao de umidade do solo id=2 na camada de raizes
    ! w(id).......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    ! itype ......Classe de textura do solo
    ! rstpar(cg,3).coefficints related to par influence on
    !                    stomatal resistance
    ! chil........leaf orientation pameter
    ! topt........Temperatura ideal de funcionamento estomatico
    ! tll.........Temperatura minima de funcionamento estomatico
    ! tu..........Temperatura maxima de funcionamento estomatico
    ! defac.......Parametro de deficit de pressao de vapor d'agua
    ! ph1.........Coeficiente para o efeito da agua no solo
    ! ph2.........Potencial de agua no solo para ponto de Wilting
    ! rootd.......Profundidade das raizes
    ! bee.........Expoente da curva de retencao "expoente para o solo umido"
    ! phsat.......Tensao do solo saturado " Potencial de agua no solo saturado"
    ! satco.......mean soil hydraulic conductivity in the root zone
    ! poros.......porosity
    ! zdepth......Profundidade para as tres camadas de solo
    ! green.......fraction of grenn leaves
    ! xcover(iv)..Fracao de cobertura de vegetacao iv=2 Bottom
    ! zlt(icg)....Indice de area foliar "LEAF AREA INDEX" icg=1 topo da copa
    ! zlt(icg)....Indice de area foliar "LEAF AREA INDEX" icg=2 base da copa
    ! x0x.........Comprimento de rugosidade
    ! xd..........Deslocamento do plano zero
    ! z2..........Altura do topo do dossel
    ! xdc.........Constant related to aerodynamic resistance
    !             between ground and canopy air space
    ! xbc.........Constant related to bulk boundary layer resistance
    ! dtc3x.......time increment dt
    ! mon.........Number of month at year (1-12)
    ! nmax
    ! ityp........numero das classes de solo 13
    ! imon........Numero maximo de meses no ano (12)
    ! icg.........Parametros da vegetacao (icg=1 topo e icg=2 base)
    ! iwv.........Compriment de onda iwv=1=visivel, iwv=2=infravermelho
    !             proximo, iwv=3 infravermelho termal
    ! idp.........Camadas de solo (1 a 3)
    ! ibd.........Estado da vegetacao ibd=1 verde / ibd=2 seco
    ! ncols.......Numero de ponto por faixa de latitude
    !
    !
    INTEGER, INTENT(in   ) :: ncols

    REAL(KIND=r8)   , INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: mon(ncols)
    INTEGER, INTENT(in   ) :: nmax

    INTEGER, INTENT(in   ) :: itype(ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(inout) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: td   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: capac(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout) :: ra    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rb    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rd    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tcta  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tgta  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ea    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: etc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: etg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: btc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: btg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: u2    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: radt  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: par   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: pd    (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: rst   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: rsoil (ncols)
    REAL(KIND=r8),    INTENT(inout) :: phroot(ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: hrr   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: phsoil(ncols,idp)
    REAL(KIND=r8),    INTENT(inout) :: cc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: satcap(ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: snow  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: dtc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dtg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dtm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dqm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: stm   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: extk  (ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(in   ) :: radfac(ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(in   ) :: closs (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gloss (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thermk(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: p1f   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: p2f   (ncols)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(inout) :: qm  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tm  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: um  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vm  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: psur(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ppc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ppl (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: radn(ncols,3,2)

    REAL(KIND=r8)   , INTENT(in   ) :: zlt2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: green2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: chil2   (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: rstpar2(ncols,icg,iwv)
    REAL(KIND=r8),    INTENT(inout) :: vcover  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: z0x(ncols)
    REAL(KIND=r8),    INTENT(inout) :: d  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rdc(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rbc(ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0 (ncols)
    REAL(KIND=r8)   , INTENT(in   ) :: topt2   (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tll2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tu2     (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: defac2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph12    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph22    (ncols,icg)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: bps   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: psb   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: em    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: gmt   (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: gmq   (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cosz  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rhoair(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: psy   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rcp   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: fc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: fg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: hr    (ncols)

    REAL(KIND=r8),    INTENT(inout) :: ect   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: eci   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egt   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egi   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egs   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ec    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: eg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecidif(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egidif(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hflux (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: chf   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: shf   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: roff  (ncols)

    INTEGER, INTENT(inout  ) :: icheck(ncols)
    !
    !     derivatives
    !
    REAL(KIND=r8),    INTENT(inout  ) :: hgdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hgdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hgdtm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hcdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hcdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hcdtm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egdqm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecdtg (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecdtc (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecdqm (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: deadtg(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: deadtc(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: deadqm(ncols)
    LOGICAL, INTENT(inout  ) :: jstneu
    REAL(KIND=r8),    INTENT(inout  ) :: ct(ncols)

    !
    CALL root(phroot,phsoil,w     ,itype ,nmax  ,ncols  )

    CALL raduse(radt  ,par   ,pd    ,radfac,closs ,gloss ,thermk,p1f   , &
         p2f   ,radn  ,vcover,nmax  ,ncols   )

    CALL stomat(cosz  ,par   ,pd    ,rst   ,extk  ,vcover,itype , &
         nmax  ,ncols ,zlt2  ,green2,chil2 ,rstpar2)

    CALL interc( &
         roff  ,cc    ,cg    ,satcap,snow  ,extk  ,tc    ,tg    ,td    , &
         capac ,w     ,tm    ,ppc   ,ppl   ,vcover,itype ,dtc3x , &
         nmax  ,ncols ,zlt2 )
    !
    !     surface flux
    !
    CALL sflxes( &
         hgdtg ,hgdtc ,hgdtm ,hcdtg ,hcdtc ,hcdtm ,egdtg ,egdtc ,egdqm , &
         ecdtg ,ecdtc ,ecdqm ,deadtg,deadtc,deadqm,icheck,bps   ,psb   , &
         dzm   ,em    ,gmt   ,gmq   ,cu    ,cuni  ,ctni  ,ustar ,rhoair, &
         psy   ,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    ,ect   ,eci   , &
         egt   ,egi   ,egs   ,ec    ,eg    ,hc    ,hg    ,ecidif,egidif, &
         ecmass,egmass,etmass,hflux ,chf   ,shf   ,ra    ,rb    ,rd    , &
         rc    ,rg    ,tcta  ,tgta  ,ta    ,ea    ,etc   ,etg   ,btc   , &
         btg   ,u2    ,radt  ,rst   ,rsoil ,hrr   ,phsoil,cc    ,cg    , &
         satcap,dtc   ,dtg   ,dtm   ,dqm   ,stm   ,thermk,tc    ,tg    , &
         td    ,capac ,w     ,qm    ,tm    ,um    ,vm    ,psur  ,vcover, &
         z0x   ,d     ,rdc   ,rbc   ,z0    ,itype ,dtc3x ,mon   ,nmax  , &
         jstneu,ncols ,zlt2  ,topt2 ,tll2  ,tu2   , defac2,ph12  ,ph22 , &
         ct)
  END SUBROUTINE pbl






  SUBROUTINE snowm(&
       chf   ,shf   ,fluxef,roff  ,cc    ,cg    ,snow  ,dtc   ,dtg   , &
       tc    ,tg    ,td    ,capac ,w     ,itype ,dtc3x ,nmax  ,ncols  )
    !
    ! snowm  :calculates snowmelt and modification of temperatures;
    !         this version deals with refreezing of water;
    !         version modified to use force-restore heat fluxes.
    !
    !-----------------------------------------------------------------------
    ! chf.........Fluxo de calor na copa (J/m*m)
    ! shf.........Fluxo de calor no solo (J/m*m)
    ! fluxef......modified to use force-restore heat fluxes
    !             fluxef(i) = shf(i) - cg(i)*dtg(i)*dtc3xi " Garrat pg. 227"
    ! roff........runoff (escoamente superficial e drenagem)(m)
    ! cc..........heat capacity of the canopy
    ! cg..........heat capacity of the ground
    ! snow........snow amount
    ! dtc ........dtc(i)=pblsib(i,2,5)*dtc3x
    ! dtg ........dtg(i)=pblsib(i,1,5)*dtc3x
    ! tc..........Temperatura da copa "dossel"(K)
    ! tg..........Temperatura da superficie do solo (K)
    ! td..........Temperatura do solo profundo (K)
    ! capac(iv)...Agua interceptada iv=1 no dossel (m)
    ! capac(iv)...Agua interceptada iv=2 na cobertura do solo (m)
    ! w(id).......Grau de saturacao de umidade do solo id=1 na camada superficial
    ! w(id).......Grau de saturacao de umidade do solo id=2 na camada de raizes
    ! w(id).......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    ! poros.......Porosidade do solo (m"3/m"3)
    ! zdepth(id)..Profundidade das camadas de solo id=1 superficial
    ! zdepth(id)..Profundidade das camadas de solo id=2 camada de raizes
    ! zdepth(id)..Profundidade das camadas de solo id=3 camada de drenagem
    ! itype.......Classe de textura do solo
    ! ncols.......Numero de ponto por faixa de latitude
    ! ityp........13
    ! icg.........Parametros da vegetacao (icg=1 topo e icg=2 base)
    ! idp.........Camadas de solo (1 a 3)
    ! snomel......Calor latente de fusao(J/kg)
    ! tf..........Temperatura de congelamento (K)
    ! dtc3x.......time increment dt
    ! nmax........
    !-----------------------------------------------------------------------
    INTEGER, INTENT(in   ) :: ncols
    REAL(KIND=r8),    INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: nmax

    INTEGER, INTENT(in   ) :: itype (ncols)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(inout) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: td   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: capac(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(in   ) :: cc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: snow  (ncols,icg)
    REAL(KIND=r8),    INTENT(in   ) :: dtc   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dtg   (ncols)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8),    INTENT(in   ) :: chf   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: shf   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: fluxef(ncols)
    REAL(KIND=r8),    INTENT(inout) :: roff  (ncols)

    REAL(KIND=r8)    :: cct   (ncols)
    REAL(KIND=r8)    :: ts         (ncols)
    REAL(KIND=r8)    :: dts   (ncols)
    REAL(KIND=r8)    :: flux  (ncols)
    REAL(KIND=r8)    :: tta   (ncols)
    REAL(KIND=r8)    :: ttb   (ncols)
    REAL(KIND=r8)    :: dtf   (ncols)
    REAL(KIND=r8)    :: work  (ncols)
    REAL(KIND=r8)    :: hf         (ncols)
    REAL(KIND=r8)    :: fcap  (ncols)
    REAL(KIND=r8)    :: spwet (ncols)
    REAL(KIND=r8)    :: dtf2  (ncols)
    REAL(KIND=r8)    :: tn         (ncols)
    REAL(KIND=r8)    :: change(ncols)
    REAL(KIND=r8)    :: dtime1(ncols)
    REAL(KIND=r8)    :: dtime2(ncols)

    INTEGER :: i
    INTEGER :: iveg
    INTEGER :: ntyp
    REAL(KIND=r8)    :: dtc3xi


    cct=0.0_r8
    dtc3xi=1.0_r8 /dtc3x

    DO iveg = 1, 2
       IF (iveg == 1) THEN

          DO i = 1, nmax
             cct (i)=cc (i)
             ts  (i)=tc (i)
             dts (i)=dtc(i)
             flux(i)=chf(i)
          END DO

       ELSE

          DO i = 1, nmax
             cct (i)=cg (i)
             ts  (i)=tg (i)
             dts (i)=dtg(i)
             flux(i)=cct(i)*dtg(i)*dtc3xi
          END DO

       END IF

       DO i = 1, nmax
          tta(i) = ts(i) - dts(i)
          ttb(i) = ts(i)
       END DO

       DO i = 1, nmax
          IF (tta(i) <= tf) THEN
             snow (i,iveg) = capac(i,iveg)
             capac(i,iveg) = 0.0_r8
          ELSE
             snow (i,iveg) = 0.0_r8
          END IF
       END DO

       DO i = 1, nmax
          work(i)=(tta(i)-tf)*(ttb(i)-tf)
       END DO

       DO i = 1, nmax
          IF (work(i) < 0.0_r8) THEN
             ntyp=itype(i)
             dtf   (i)= tf - tta(i)
             dtime1(i)= cct (i)* dtf(i)/ flux(i)
             hf    (i)= flux(i)*(dtc3x-dtime1(i))
             spwet (i)=  MIN ( 5.0_r8 , snow(i,iveg) )
             IF (dts(i) <= 0.0_r8) THEN
                fcap (i) =-capac(i,iveg)* snomel
             ELSE
                fcap (i) = spwet(i)     * snomel
             END IF
             dtime2(i)= fcap(i) / flux(i)
             dtf2  (i)= flux(i) * (dtc3x-dtime1(i)-dtime2(i))/cct(i)
             tn(i)    = tf + dtf2(i)
             IF (ABS(hf(i)) < ABS(fcap(i))) THEN
                ts(i)    = tf -0.1_r8
             ELSE
                ts(i)    = tn(i)
             END IF
             IF (ABS(hf(i)) < ABS(fcap(i))) THEN
                change(i) = hf  (i)
             ELSE
                change(i) = fcap(i)
             END IF
             change(i)     =change(i)      / snomel
             snow  (i,iveg)=snow  (i,iveg) - change(i)
             capac (i,iveg)=capac (i,iveg) + change(i)
             IF (snow(i,iveg) < 1.e-10_r8) snow(i,iveg)=0.0e0_r8
             IF (iveg == 1)THEN
                tc(i)=ts(i)
             ELSE
                tg(i)=ts(i)
             END IF
             IF (iveg == 2) THEN
                IF (td(i) > tf) THEN
                   w (i,1)=w (i,1)+capac(i,iveg) &
                        /(poros(ntyp)*zdepth(ntyp,1))
                ELSE
                   roff(i)=roff(i)+capac(i,iveg)
                END IF
                capac(i,iveg) = 0.0_r8
             END IF
          END IF
       END DO
       DO i = 1, nmax
          capac(i,iveg) =  capac(i,iveg) + snow(i,iveg)
       END DO
    END DO

    !   modified to use force-restore heat fluxes

    DO i = 1, nmax
       fluxef(i) = shf(i) - cg(i)*dtg(i)*dtc3xi
    END DO

  END SUBROUTINE snowm




  ! fysiks :it is a physics driver; performs the following:
  !         a) soil water budget prior to calling pbl
  !         b) planetary boundary layer (pbl) parameterization
  !         c) update sib variables
  !         d) dumping of small capac values onto soil surface store
  !         e) snowmelt/refreeze calculation
  !         f) update deep soil temperature using effective soil heat flux
  !         g) bare soil evaporation loss
  !         h) extraction of transpiration loss from root zone
  !         i) interflow, infiltration excess and loss to groundwater
  !         j) increment prognostic variables and
  !            adjust theta and sh to be consistent with dew formation
  !         k) calculates soil water budget after calling pbl
  !            and compares with previous budget.

  SUBROUTINE fysiks(vcover, z0x  , d    , rdc  , rbc  , z0   ,ndt   , &
       latitu, bps  ,psb   ,dzm   ,em    ,gmt   ,gmq   , &
       gmu   ,cu    , cuni ,ctni  ,ustar ,cosz  ,sinclt,rhoair, &
       psy   ,rcp   , wc   ,wg    ,fc    ,fg    ,hr    , ect  , &
       eci   , egt  , egi  , egs  , ec   , eg   , hc   , hg   , &
       ecidif,egidif,ecmass,egmass,etmass,hflux , chf  , shf  , &
       fluxef, roff , drag ,ra    , rb   , rd   , rc   , rg   , &
       tcta  , tgta , ta   , ea   , etc  , etg  , btc  , btg  , &
       u2    , radt , par  , pd   , rst  ,rsoil ,phroot, hrr  , &
       phsoil, cc   , cg   ,satcap, snow , dtc  , dtg  , dtm  , &
       dqm   , stm  , extk ,radfac, closs,gloss ,thermk, p1f  , &
       p2f   , tc   , tg   , td   , capac, w    ,  qm  , tm   , &
       um    , vm   , psur , ppc  , ppl  , radn ,itype ,dtc3x , &
       mon   , nmax , ncols,zlt2  ,green2,chil2 ,rstpar2,topt2, &
       tll2  ,tu2   , defac2,ph12  ,ph22 ,bstar)
    !
    !
    !-----------------------------------------------------------------------
    !
    !  roff.......Runoff (escoamente superficial e drenagem)(m)
    !  slope......Inclinacao de perda hidraulica na camada profunda do solo
    !  bee........Fator de retencao da umidade no solo (expoente da umidade do
    !             solo)
    !  satco......Condutividade hidraulica do solo saturado(m/s)
    !  zdepth(id).Profundidade das camadas de solo id=1 superficial
    !  zdepth(id).Profundidade das camadas de solo id=2 camada de raizes
    !  zdepth(id).Profundidade das camadas de solo id=3 camada de drenagem
    !  phsat......Potencial matricial do solo saturado(m) (tensao do solo em
    !             saturacao)
    !  poros......Porosidade do solo
    !  dtc3x......time increment dt
    !  snomel.....Calor latente de fusao(J/kg)
    !  w(id)......Grau de saturacao de umidade do solo id=1 na camada superficial
    !  w(id)......Grau de saturacao de umidade do solo id=2 na camada de raizes
    !  w(id)......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    !  capac(iv)..Agua interceptada iv=1 no dossel (m)
    !  capac(iv)..Agua interceptada iv=2 na cobertura do solo (m)
    !  tg.........Temperatura da superficie do solo  (K)
    !  td.........Temperatura do solo profundo (K)
    !  itype......Classe de textura do solo
    !  tf.........Temperatura de congelamento (K)
    !  idp........Parametro para as camadas de solo idp=1->3
    !  nmax.......
    !  ncols......Number of grid points on a gaussian latitude circle
    !  ityp.......Numero das classes de solo 13
    !  imon.......Numero maximo de meses no ano (12)
    !  icg........Parametros da vegetacao (icg=1 topo e icg=2 base)
    !  iwv........Compriment de onda iwv=1=visivel, iwv=2=infravermelho
    !             proximo, iwv=3 infravermelho termal
    !  idp........Camadas de solo (1 a 3)
    !  ibd........Estado da vegetacao ibd=1 verde / ibd=2 seco
    !  pie........Constante Pi=3.1415926e0
    !  stefan.....Constante de Stefan Boltzmann
    !  cp.........specific heat of air (j/kg/k)
    !  hl ........heat of evaporation of water   (j/kg)
    !  grav.......gravity constant      (m/s**2)
    !  snomel.....heat of melting (j m-1)
    !  tf.........Temperatura de congelamento (K)
    !  clai.......heat capacity of foliage
    !  cw.........liquid water heat capacity               (j/m**3)
    !  gasr.......Constant of dry air      (j/kg/k)
    !  epsfac.....Constante 0.622 Razao entre as massas moleculares do vapor
    !             de agua e do ar seco
    !  athird.....Constante athird=1.0e0 /3.0e0
    !  dtc3x......time increment dt
    !  mon........Number of month at year (1-12)
    !  nmax
    !  rstpar.....Coefficints related to par influence on
    !             stomatal resistance
    !  chil.......Leaf orientation parameter
    !  topt.......Temperatura ideal de funcionamento estomatico
    !  tll........Temperatura minima de funcionamento estomatico
    !  tu.........Temperatura maxima de funcionamento estomatico
    !  defac......Parametro de deficit de pressao de vapor d'agua
    !  ph1........Coeficiente para o efeito da agua no solo
    !  ph2........Potencial de agua no solo para ponto de Wilting
    !  rootd......Profundidade das raizes
    !  bee........Expoente da curva de retencao "expoente para o solo umido"
    !  phsat......Tensao do solo saturado " Potencial de agua no solo saturado"
    !  satco......mean soil hydraulic conductivity in the root zone
    !  poros......Porosity
    !  zdepth.....Profundidade para as tres camadas de solo
    !  green......Fraction of grenn leaves
    !  xcover(iv).Fracao de cobertura de vegetacao iv=1 Top
    !  xcover(iv).Fracao de cobertura de vegetacao iv=2 Bottom
    !  zlt(icg)...Indice de area foliar "LEAF AREA INDEX" icg=1 topo da copa
    !  zlt(icg)...Indice de area foliar "LEAF AREA INDEX" icg=2 base da copa
    !  x0x........Comprimento de rugosidade
    !  xd.........Deslocamento do plano zero
    !  z2.........Altura do topo do dossel
    !  xdc........Constant related to aerodynamic resistance
    !             between ground and canopy air space
    !  xbc........Constant related to bulk boundary layer resistance
    !  itype......Classe de textura do solo
    !  qm.........Reference specific humidity (fourier)
    !  tm.........Reference temperature    (fourier)                (k)
    !  um.........Razao entre zonal pseudo-wind (fourier) e seno da
    !             colatitude
    !  vm.........Razao entre meridional pseudo-wind (fourier) e seno da
    !             colatitude
    !  psur.......Surface pressure in mb
    !  ppc........Precipitation rate ( cumulus )           (mm/s)
    !  ppl........Precipitation rate ( large scale )       (mm/s)
    !  radn.......Downward sw/lw radiation at the surface
    !  tc.........Temperatura da copa "dossel"(K)
    !  tg.........Temperatura da superficie do solo (K)
    !  td.........Temperatura do solo profundo (K)
    !  capac(iv)..Agua interceptada iv=1 no dossel "water store capacity
    !             of leaves"(m)
    !  capac(iv)..Agua interceptada iv=2 na cobertura do solo (m)
    !  w(id)......Grau de saturacao de umidade do solo id=1 na camada superficial
    !  w(id)......Grau de saturacao de umidade do solo id=2 na camada de raizes
    !  w(id)......Grau de saturacao de umidade do solo id=3 na camada de drenagem
    !  ra.........Resistencia Aerodinamica (s/m)
    !  rb.........bulk boundary layer resistance
    !  rd.........Aerodynamic resistance between ground      (s/m)
    !             and canopy air space
    !  rc.........Resistencia do topo da copa
    !  rg.........Resistencia da base da copa
    !  tcta.......Diferenca entre tc-ta                      (k)
    !  tgta.......Diferenca entre tg-ta                      (k)
    !  ta.........Temperatura no nivel de fonte de calor do dossel (K)
    !  ea.........Pressure of vapor
    !  etc........Pressure of vapor at top of the copa
    !  etg........Pressao de vapor no base da copa
    !  btc........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i)).
    !  btg........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))
    !  u2.........wind speed at top of canopy
    !  radt.......net heat received by canopy/ground vegetation
    !  par........par incident on canopy
    !  pd.........ratio of par beam to total par
    !  rst .......Resisttencia Estomatica "Stomatal resistence" (s/m)
    !  rsoil......Resistencia do solo (s/m)
    !  phroot.....Soil moisture potentials in root zone of each
    !             vegetation layer and summed soil+root resistance.
    !  hrr........rel. humidity in top layer
    !  phsoil.....soil moisture potential of the i-th soil layer
    !  cc.........heat capacity of the canopy
    !  cg.........heat capacity of the ground
    !  satcap.....saturation liquid water capacity         (m)
    !  snow.......snow amount
    !  dtc........dtc(i)=pblsib(i,2,5)*dtc3x
    !  dtg........dtg(i)=pblsib(i,1,5)*dtc3x
    !  dtm........dtm(i)=pblsib(i,3,5)*dtc3x
    !  dqm .......dqm(i)=pblsib(i,4,5)*dtc3x
    !  stm .......Variavel utilizada mo cal. da Resisttencia
    !  extk.......extinction coefficient
    !  radfac.....Fractions of downward solar radiation at surface
    !             passed from subr.radalb
    !  closs......Radiation loss from canopy
    !  gloss......Radiation loss from ground
    !  thermk.....Canopy emissivity
    !  p1f
    !  p2f
    !  ect........Transpiracao no topo da copa (J/m*m)
    !  eci........Evaporacao da agua interceptada no topo da copa (J/m*m)
    !  egt........Transpiracao na base da copa (J/m*m)
    !  egi........Evaporacao da neve (J/m*m)
    !  egs........Evaporacao do solo arido (J/m*m)
    !  ec.........Soma da Transpiracao e Evaporacao da agua interceptada pelo
    !             topo da copa   ec   (i)=eci(i)+ect(i)
    !  eg.........Soma da transpiracao na base da copa +  Evaporacao do solo arido
    !             +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"
    !  hc.........Total sensible heat lost of top from the veggies.
    !  hg.........Total sensible heat lost of base from the veggies.
    !  ecidif.....check if interception loss term has exceeded canopy storage
    !             ecidif(i)=MAX(0.0   , eci(i)-capac(i,1)*hlat3 )
    !  egidif.....check if interception loss term has exceeded canopy storage
    !             ecidif(i)=MAX(0.0   , egi(i)-capac(i,1)*hlat3 )
    !  ecmass.....Mass of water lost of top from the veggies.
    !  egmass.....Mass of water lost of base from the veggies.
    !  etmass.....Total mass of water lost from the veggies.
    !  hflux......Total sensible heat lost from the veggies
    !  chf........Heat fluxes into the canopy  in w/m**2
    !  shf........Heat fluxes into the ground, in w/m**2
    !  fluxef.....Modified to use force-restore heat fluxes
    !             fluxef(i) = shf(i) - cg(i)*dtg(i)*dtc3xi " Garrat pg. 227"
    !  roff.......runoff (escoamente superficial e drenagem)(m)
    !  drag.......tensao superficial
    !  bps
    !  psb
    !  dzm........Altura media de referencia  para o vento para o calculo
    !             da estabilidade do escoamento
    !  em.........Pressao de vapor da agua
    !  gmt(i,k,3).temperature related matrix virtual temperature tendency
    !             due to vertical diffusion
    !  gmq........specific humidity related matrix specific humidity of
    !             reference (fourier)
    !  gmu........wind related matrix
    !  cu.........Friction  transfer coefficients.
    !  cuni.......Neutral friction transfer  coefficients.
    !  ctni.......Neutral heat transfer coefficients.
    !  ustar......Surface friction velocity  (m/s)
    !  cosz.......Cosine of zenith angle
    !  sinclt.....sinclt=SIN(colrad(latitu))"seno da colatitude"
    !  rhoair.....Desnsidade do ar
    !  psy........(cp/(hl*epsfac))*psur(i)
    !  rcp........densidade do ar vezes o calor especifico do ar
    !  wc.........Minimo entre 1 e a razao entre a agua interceptada pelo
    !             indice de area foliar no topo da copa
    !  wg.........Minimo entre 1 e a razao entre a agua interceptada pelo
    !             indice de area foliar na base da copa
    !  fc.........Condicao de oravalho 0 ou 1 na topo da copa
    !  fg.........Condicao de oravalho 0 ou 1 na base da copa
    !  hr.........rel. humidity in top layer
    !  ndt
    !  latitu
    !  jstneu.....The first call to vntlat just gets the neutral values of ustar
    !             and ventmf para jstneu=.TRUE..
    !  hgdtg.......n.b. fluxes expressed in joules m-2
    !  hgdtc.......n.b. fluxes expressed in joules m-2
    !  hgdtm.......n.b. fluxes expressed in joules m-2
    !  hcdtg.......n.b. fluxes expressed in joules m-2
    !  hcdtc.......n.b. fluxes expressed in joules m-2
    !  hcdtm.......n.b. fluxes expressed in joules m-2
    !  egdtg.......partial derivative calculation for latent heat
    !  egdtc.......partial derivative calculation for latent heat
    !  egdqm.......partial derivative calculation for latent heat
    !  ecdtg.......partial derivative calculation for latent heat
    !  ecdtc.......partial derivative calculation for latent heat
    !  ecdqm.......partial derivative calculation for latent heat
    !  deadtg
    !  deadtc
    !  deadqm
    !  icheck......this version assumes dew-free conditions "icheck=1" to
    !              estimate ea for buoyancy term in vntmf or ra.
    !  vcover(iv)..Fracao de cobertura de vegetacao iv=1 Top
    !  vcover(iv)..Fracao de cobertura de vegetacao iv=2 Botto
    !  z0x.........roughness length
    !  d...........Displacement height
    !  rdc.........Constant related to aerodynamic resistance
    !              between ground and canopy air space
    !  rbc.........Constant related to bulk boundary layer resistance
    !  z0..........Roughness length
    !-----------------------------------------------------------------------
    !
    INTEGER, INTENT(in   ) :: ncols

    REAL(KIND=r8)   , INTENT(in   ) :: dtc3x
    INTEGER, INTENT(in   ) :: mon(ncols)
    INTEGER, INTENT(in   ) :: nmax

    INTEGER, INTENT(in   ) :: itype (ncols)
    !
    !     the size of working area is ncols*187
    !     atmospheric parameters as boudary values for sib
    !
    REAL(KIND=r8),    INTENT(inout) :: qm  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tm  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: um  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: vm  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: psur(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ppc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ppl (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: radn(ncols,3,2)
    !
    !     prognostic variables
    !
    REAL(KIND=r8),    INTENT(inout) :: tc   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: td   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: capac(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: w    (ncols,3)
    !
    !     variables calculated from above and ambient conditions
    !
    REAL(KIND=r8),    INTENT(inout) :: ra    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rb    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rd    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: tcta  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: tgta  (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ta    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ea    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etc   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etg   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: btc   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: btg   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: u2    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: radt  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: par   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: pd    (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: rst   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: rsoil (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: phroot(ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: hrr   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: phsoil(ncols,idp)
    REAL(KIND=r8),    INTENT(inout  ) :: cc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: cg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: satcap(ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: snow  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout  ) :: dtc   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dtg   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dtm   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: dqm   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: stm   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: extk  (ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(in   ) :: radfac(ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(in   ) :: closs (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gloss (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thermk(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: p1f   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: p2f   (ncols)
    !
    !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
    !
    REAL(KIND=r8),    INTENT(inout) :: ect   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: eci   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egt   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egi   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: egs   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ec    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: eg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hc    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hg    (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecidif(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egidif(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: ecmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: egmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: etmass(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: hflux (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: chf   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: shf   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: fluxef(ncols)
    REAL(KIND=r8),    INTENT(inout) :: roff  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: drag  (ncols)
    !
    !     this is for coupling with closure turbulence model
    !
    REAL(KIND=r8),    INTENT(in   ) :: bps   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: psb   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: em    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: gmt   (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: gmq   (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: gmu   (ncols,4)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: cosz  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: rhoair(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: psy   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rcp   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: fc    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: fg    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: hr    (ncols)

    INTEGER, INTENT(in   ) :: ndt
    INTEGER, INTENT(in   ) :: latitu

    REAL(KIND=r8)   , INTENT(in   ) :: rstpar2 (ncols,icg,iwv)
    REAL(KIND=r8)   , INTENT(in   ) :: zlt2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: green2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: chil2   (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: vcover  (ncols,icg)
    REAL(KIND=r8),    INTENT(inout) :: z0x(ncols)
    REAL(KIND=r8),    INTENT(inout) :: d  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rdc(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rbc(ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0 (ncols)
    REAL(KIND=r8)   , INTENT(in   ) :: topt2   (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tll2    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: tu2     (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: defac2  (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph12    (ncols,icg)
    REAL(KIND=r8)   , INTENT(in   ) :: ph22    (ncols,icg)
    REAL(KIND=r8)   , INTENT(out   ) :: bstar(ncols)
    LOGICAL :: jstneu
    INTEGER :: icheck(ncols)

    !
    !     derivatives
    !
    REAL(KIND=r8) :: hgdtg (ncols)
    REAL(KIND=r8) :: hgdtc (ncols)
    REAL(KIND=r8) :: hgdtm (ncols)
    REAL(KIND=r8) :: hcdtg (ncols)
    REAL(KIND=r8) :: hcdtc (ncols)
    REAL(KIND=r8) :: hcdtm (ncols)
    REAL(KIND=r8) :: egdtg (ncols)
    REAL(KIND=r8) :: egdtc (ncols)
    REAL(KIND=r8) :: egdqm (ncols)
    REAL(KIND=r8) :: ecdtg (ncols)
    REAL(KIND=r8) :: ecdtc (ncols)
    REAL(KIND=r8) :: ecdqm (ncols)
    REAL(KIND=r8) :: deadtg(ncols)
    REAL(KIND=r8) :: deadtc(ncols)
    REAL(KIND=r8) :: deadqm(ncols)

    REAL(KIND=r8)    :: ef    (ncols,3)
    REAL(KIND=r8)    :: absoil(ncols)
    REAL(KIND=r8)    :: totdep(ncols)
    REAL(KIND=r8)    :: div   (ncols)
    REAL(KIND=r8)    :: eft   (ncols)
    REAL(KIND=r8)    :: aaa   (ncols)
    REAL(KIND=r8)    :: dep   (ncols)


    INTEGER :: i
    INTEGER :: il
    INTEGER :: ntyp
    INTEGER :: iveg
    REAL(KIND=r8)    :: hlat3i
    REAL(KIND=r8)    :: gby100
    REAL(KIND=r8)    :: timcon
    REAL(KIND=r8)    :: totwb(ncols)
    REAL(KIND=r8)    :: endwb(ncols)
    REAL(KIND=r8)    :: cbal (ncols)
    REAL(KIND=r8)    :: gbal (ncols)
    REAL(KIND=r8)    :: ct(ncols)
    REAL(KIND=r8)    :: qh(ncols)
    REAL(KIND=r8)    :: speed(ncols)
    REAL(KIND=r8)    :: d1
    !
    !     calculates soil water budget prior to calling pbl
    !
    DO i = 1, nmax
       !
       !  capac(1)..Agua interceptada no dossel (m)
       !  capac(2)..Agua interceptada na cobertura do solo (m)
       !
       totwb(i)=w(i,1)*poros(itype(i))*zdepth(itype(i),1) &
            +w(i,2)*poros(itype(i))*zdepth(itype(i),2) &
            +w(i,3)*poros(itype(i))*zdepth(itype(i),3) &
            +capac(i,1) + capac(i,2)
    END DO
    !
    !     planetary boundary layer parameterization
    !
    CALL pbl(jstneu, hgdtg , hgdtc , hgdtm , hcdtg , hcdtc , hcdtm , &
         egdtg , egdtc , egdqm , ecdtg , ecdtc , ecdqm , deadtg, &
         deadtc, deadqm,icheck , ect   , eci   , egt   , egi   , &
         egs   , ec    , eg    , hc    , hg    , ecidif, egidif, &
         ecmass, egmass, etmass, hflux , chf   , shf   , roff  , &
         bps   , psb   , dzm   , em    , gmt   , gmq   ,  cu   , &
         cuni  , ctni  , ustar , cosz  , rhoair, psy   , rcp   , &
         wc   , wg       , fc       , fg       , hr       , vcover, z0x   , &
         d       , rdc   , rbc   , z0       , qm       , tm       , um       , &
         vm       , psur  , ppc   , ppl   , radn  , ra       , rb       , &
         rd       , rc       , rg       , tcta  , tgta  , ta       , ea       , &
         etc   , etg   , btc   , btg   , u2       , radt  , par   , &
         pd       , rst   , rsoil , phroot,  hrr  , phsoil, cc       , &
         cg       , satcap, snow  , dtc   , dtg   , dtm   , dqm   , &
         stm   , extk  , radfac, closs , gloss , thermk, p1f   , &
         p2f   , tc       , tg       , td       , capac , w       , itype , &
         dtc3x , mon   , nmax  , ncols ,zlt2  ,green2,chil2 ,rstpar2,&
         topt2 ,tll2  ,tu2   , defac2,ph12  ,ph22 ,ct)

    !
    !     continue to update sib variables
    !
    DO i = 1, nmax
       tc(i) = tc(i) + dtc(i)
       tg(i) = tg(i) + dtg(i)
    END DO
    !
    !     dumping of small capac values onto soil surface store
    !
    DO iveg = 1, 2
       DO i = 1, nmax
          ntyp  =itype(i)
          IF (capac(i,iveg) <= 1.e-6_r8)THEN
             w(i,1)=w(i,1)+capac(i,iveg)/(poros(ntyp)*zdepth(ntyp,1))
             capac(i,iveg)=0.0_r8
          END IF
       END DO
    END DO
    !
    !     snowmelt/refreeze calculation
    !
    CALL snowm(&
         chf   ,shf   ,fluxef,roff  ,cc    ,cg    ,snow  ,dtc   ,dtg   , &
         tc    ,tg    ,td    ,capac ,w     ,itype ,dtc3x ,nmax  ,ncols   )
    !
    !     update deep soil temperature using effective soil heat flux
    !
    timcon=dtc3x/(2.0_r8 *SQRT(pie*365.0_r8 ))

    DO i = 1, nmax
       td(i)=td(i)+fluxef(i)/cg(i)*timcon
    END DO
    !
    !     bare soil evaporation loss
    !
    hlat3i=1.0_r8/(hl*1000.0_r8 )
    DO i = 1, nmax
       ntyp=itype(i)
       w(i,1)=w(i,1)-egs(i)*hlat3i/(poros(ntyp)*zdepth(ntyp,1))
    END DO
    !
    !        extraction of transpiration loss from root zone
    !
    DO iveg = 1, 2
       IF (iveg == 1) THEN
          DO i = 1, nmax
             absoil(i)=ect(i)*hlat3i
          END DO
       ELSE
          DO i = 1, nmax
             absoil(i)=egt(i)*hlat3i
          END DO
       END IF
       DO i = 1, nmax
          ntyp=itype(i)
          ef(i,2)=0.0_r8
          ef(i,3)=0.0_r8
          totdep(i)=zdepth(ntyp,1)
       END DO
       DO il = 2, 3
          DO i = 1, nmax
             ntyp=itype(i)
             totdep(i)=totdep(i)+zdepth(ntyp,il)
             div(i)=rootd(ntyp,iveg)
             dep(i)=MAX(0.0_r8  ,rootd(ntyp,iveg)-totdep(i)+ &
                  zdepth(ntyp,il))
             dep(i)=MIN(dep(i),zdepth(ntyp,il))
             ef(i,il)=dep(i)/div(i)
          END DO
       END DO
       DO i = 1, nmax
          eft(i  )=ef(i,2)+ef (i,3)
          eft(i) = MAX(eft(i),0.1e-5_r8)
          ef (i,2)=ef(i,2)/eft(i)
          ef (i,3)=ef(i,3)/eft(i)
       END DO
       DO il = 2, 3
          DO i = 1, nmax
             ntyp=itype(i)
             w(i,il)=w(i,il)-absoil(i)*ef(i,il)/ &
                  (poros(ntyp)*zdepth(ntyp,il))
          END DO
       END DO
    END DO
    !
    !     interflow, infiltration excess and loss to
    !     groundwater .  all losses are assigned to variable 'roff' .
    !
    DO il = 1, 2
       DO i = 1, nmax
          IF (w(i,il) <= 0.0_r8) THEN
             ntyp=itype(i)
             w(i,il+1)=w(i,il+1)+w(i,il)* &
                  zdepth(ntyp,il)/zdepth(ntyp,il+1)
             w(i,il  )=0.0_r8
          END IF
       END DO
    END DO

    CALL runoff(&
         roff  ,tg    ,td    ,capac ,w     ,itype ,dtc3x ,nmax  ,ncols )

    DO i = 1, nmax
       ntyp = itype(i)
       IF (w(i,1) > 1.0_r8) THEN
          w(i,2)=w(i,2)+(w(i,1)-1.0_r8 )*zdepth(ntyp,1)/zdepth(ntyp,2)
          w(i,1)=1.0_r8
       ENDIF
       IF (w(i,2) > 1.0_r8) THEN
          w(i,3)= w(i,3)+(w(i,2)-1.0_r8 )*zdepth(ntyp,2)/zdepth(ntyp,3)
          w(i,2)=1.0_r8
       ENDIF
       IF (w(i,3) > 1.0_r8) THEN
          roff(i)=roff(i)+(w(i,3)-1.0_r8 )*poros(ntyp)*zdepth(ntyp,3)
          w(i,3)=1.0_r8
       END IF
    END DO
    !
    !     increment prognostic variables
    !
    !     adjust theta and sh to be consistent with dew formation
    !
    gby100 = 0.01_r8  * grav

    DO i = 1, nmax
       !
       !     solve implicit system for winds
       !
       ! psb(i) = psur(i) * ( si(k) - si(k+1) )
       !
       drag(i)  =rhoair(i)*cu(i)*ustar(i)
       !
       ! P=rho*G*Z ===> DP=rho*G*DZ
       !
       ! D                D
       !---- = rho * g * ----
       ! DZ               DP
       !                                D
       ! aaa = cu * ustar *  rho * g * ----
       !                                DP
       !
       !                                                   g
       ! aaa (i)  = rhoair(i)*cu(i)*ustar(i) * -------------------------------
       !                                       100*psur(i) * ( si(k) - si(k+1) )
       !
       aaa (i)  =drag  (i)*gby100/psb(i)

       gmu (i,2) =  gmu(i,2) + dtc3x*aaa(i)
       gmu (i,3) = (gmu(i,3) - aaa(i) * um(i)*sinclt(i) ) / gmu(i,2)
       gmu (i,4) = (gmu(i,4) - aaa(i) * vm(i)*sinclt(i) ) / gmu(i,2)

       d1     =1.0_r8/ra(i) + 1.0_r8/rb(i) + 1.0_r8/rd(i)

       ta(i)  =( tg(i)/rd(i) + tc(i)/rb(i) + tm(i)*bps(i)/ra(i) )/d1
       speed (i) = SQRT(um(i)**2 + vm(i)**2)
       speed (i) = MAX(2.0_r8  ,speed(i))
       qh(i)=0.622e0_r8*EXP(21.65605e0_r8 -5418.0e0_r8 /ta(i))/psur(i)

       !THAT  =effective_surface_skin_temperature [K]
       !BSTAR = (grav/(rhos*sqrt(CM*max(UU,1.e-30)/RHOS))) *  &
       !(CT*(TH-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ*(QH-QA))
       !
       !                            grav
       !BSTAR = -------------------------------------------------------------- *
       !         (rho(i)*sqrt(CU(i)*max(speed(i),1.e-30_r8)/RHOS(i)))
       !
       !
       !
       !
       !         (CT(i)*(tsfc(i)-gt(i)-(GRAV/CP)*dzm(i))/gt(i) + MAPL_VIREPS*CT(i)*(qsfc-gq(i)))
       !
       !
       !          cuni(i)=LOG((dzm(i)-d(i))/z0(i))*vkrmni
       ! 
       !
       !            0.4                    T0 -T
       !b_star = ------------------ * G *--------
       !         log(z(i)/z0(i))            T0
       !
       !bstar(i)=cu(i)*grav*(ct(i)*(ta(i)-tm(i)*BPS(I)-(grav/cp)*dzm(i))/tm(i)*BPS(I))!+mapl_vireps*ct(i)*(qh(i)-qm(i)))
        bstar(i)=cu(i)*grav*(ct(i)*(ta(i)/BPS(I)-tm(i))/tm(i)*BPS(I))!+mapl_vireps*ct(i)*(qh(i)-qm(i)))

       !bstar(i) = (grav/(rhoair(i)*sqrt(cu(i)*max(speed(i),1.e-30_r8)/rhoair(i)))) + &
       !             (ct(i)*(tg(i)-tm(i)-(grav/cp)*dzm(i))/tm(i) + mapl_vireps*ct(i)*(qh(i)-qm(i)))


    END DO
    !
    !     calculates soil water budget after calling pbl
    !     and compares with previous budget
    !
    DO i = 1, nmax
       ntyp=itype(i)
       endwb(i)=w(i,1)*poros(ntyp)*zdepth(ntyp,1) &
            +w(i,2)*poros(ntyp)*zdepth(ntyp,2) &
            +w(i,3)*poros(ntyp)*zdepth(ntyp,3) &
            +capac(i,1)+capac(i,2) &
            -(ppl(i)+ppc(i))/1000.0_r8 + etmass(i)/1000.0_r8 + roff(i)
       !IF (ABS(totwb(i)-endwb(i)) > 0.0001_r8) THEN
       !  WRITE(UNIT=nfprt,FMT=998) latitu,i,ntyp,ndt, &
       !       totwb(i),endwb(i),(totwb(i)-endwb(i)),w(i,1),w(i,2), &
       !       w(i,3),capac(i,1),capac(i,2),ppl(i),ppc(i),etmass(i), &
       !       roff(i),zlt(ntyp,12,1),zlt(ntyp,12,2), &
       !       tc(i),tg(i),td(i),tm(i)
       !END IF
       !
       !     calculates and compares energy budgets
       !
       cbal(i)=radt(i,1)-chf(i)-(ect(i)+hc(i)+eci(i))/dtc3x
       gbal(i)=radt(i,2)-shf(i)-(egt(i)+egi(i)+hg(i)+egs(i))/dtc3x
       !IF (ABS(cbal(i)-gbal(i)) > 5.0_r8) &
       !    WRITE(UNIT=nfprt,FMT=999)latitu,i,ntyp,ndt, &
       !    radt(i,1),radt(i,2),chf(i),shf(i),hflux(i), &
       !    ect(i),eci(i),egt(i),egi(i),egs(i)
    END DO
    !cdir critical
    DO i=1,nmax
       ntyp=itype(i)
       ! if(abs(totwb(i)-endwb(i)).gt.0.0001_r8) then
       IF(ABS(totwb(i)-endwb(i)).GT.0.0005_r8) THEN
          WRITE(UNIT=nfprt,FMT=998) latitu,i,ntyp,ndt, &
               totwb(i),endwb(i),(totwb(i)-endwb(i)),w(i,1),w(i,2), &
               w(i,3),capac(i,1),capac(i,2),ppl(i),ppc(i),etmass(i), &
               roff(i),zlt(ntyp,12,1),zlt(ntyp,12,2), &
               tc(i),tg(i),td(i),tm(i)
       END IF
       IF(ABS(cbal(i)-gbal(i)).GT.5.0_r8) &
            WRITE(UNIT=nfprt,FMT=999)latitu,i,ntyp,ndt, &
            radt(i,1),radt(i,2),chf(i),shf(i),hflux(i), &
            ect(i),eci(i),egt(i),egi(i),egs(i)
    END DO
    !cdir end critical

998    FORMAT(3I4,1X,'WATER BAL.',I8,/3E12.4/3E12.4/2E12.4/4E12.4/2E12.4/4E12.4)
999    FORMAT(3I4,1X,'ENERGY BAL.',I8/4E12.3/6E12.3)
  END SUBROUTINE fysiks


  SUBROUTINE sextrp &
       (td    ,tg    ,tc    ,w     ,capac ,td0   ,tg0   ,tc0   ,w0    , &
       capac0,tdm   ,tgm   ,tcm   ,wm    ,capacm,istrt ,ncols ,nmax  , &
       epsflt,intg  ,latitu,tm0   ,qm0   ,tm    ,qm    ,tmm    ,qmm     )
    INTEGER, INTENT(in   ) :: istrt
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: nmax
    REAL(KIND=r8)   , INTENT(in   ) :: epsflt
    INTEGER, INTENT(in   ) :: intg
    INTEGER, INTENT(in   ) :: latitu
    REAL(KIND=r8),    INTENT(in   ) :: tm    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: qm    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: td    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tg    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tc    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: w     (ncols,3)
    REAL(KIND=r8),    INTENT(in   ) :: capac (ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: td0   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tg0   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tc0   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: w0    (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: capac0(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: tdm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tgm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tcm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: wm    (ncols,3)
    REAL(KIND=r8),    INTENT(inout) :: capacm(ncols,2)
    REAL(KIND=r8),    INTENT(inout) :: tm0 (ncols)
    REAL(KIND=r8),    INTENT(inout) :: qm0 (ncols)
    REAL(KIND=r8),    INTENT(inout) :: tmm (ncols)
    REAL(KIND=r8),    INTENT(inout) :: qmm (ncols)
    INTEGER :: i, nc, ii
    ii=0
    IF (intg == 2) THEN
       IF (istrt >= 1) THEN
          DO i = 1, nmax
             tm0   (i)  =tm   (i)
             qm0   (i)  =qm   (i)
             td0   (i)  =td   (i)
             tg0   (i)  =tg   (i)
             tc0   (i)  =tc   (i)
             w0    (i,1)=w    (i,1)
             w0    (i,2)=w    (i,2)
             w0    (i,3)=w    (i,3)
             capac0(i,1)=capac(i,1)
             capac0(i,2)=capac(i,2)
             IF (capac0(i,2) > 0.0_r8 .AND. tg0(i) > 273.16_r8) THEN
                nc=0
                ii=0
                !DO ii = 1, ncols
                !   IF (imask(ii) >= 1) nc=nc+1
                !   IF (nc == i) EXIT
                !END DO
                WRITE(UNIT=nfprt,FMT=200)ii,latitu,i,capac0(i,2),tg0(i)
             END IF
          END DO
       ELSE
          DO i = 1, nmax
             td0(i)=td0(i)+epsflt*(td(i)+tdm(i)-2.0_r8  *td0(i))
             tg0(i)=tg0(i)+epsflt*(tg(i)+tgm(i)-2.0_r8  *tg0(i))
             tc0(i)=tc0(i)+epsflt*(tc(i)+tcm(i)-2.0_r8  *tc0(i))

             tm0(i)=tm0(i)+epsflt*(tm(i)+tmm(i)-2.0_r8  *tm0(i))
             qm0(i)=qm0(i)+epsflt*(qm(i)+qmm(i)-2.0_r8  *qm0(i))

             IF(w0    (i,1) > 0.0_r8 ) THEN
                w0(i,1)=w0(i,1)+epsflt*(w(i,1)+wm(i,1)-2.0_r8  *w0(i,1))
             END IF
             IF(w0    (i,2) > 0.0_r8 ) THEN
                w0(i,2)=w0(i,2)+epsflt*(w(i,2)+wm(i,2)-2.0_r8  *w0(i,2))
             END IF
             IF(w0    (i,3) > 0.0_r8 ) THEN
                w0(i,3)=w0(i,3)+epsflt*(w(i,3)+wm(i,3)-2.0_r8  *w0(i,3))
             END IF
             IF(capac0(i,1) > 0.0_r8 ) THEN
                capac0(i,1)=capac0(i,1) &
                     +epsflt*(capac(i,1)+capacm(i,1)-2.0_r8*capac0(i,1))
             END IF
             IF(capac0(i,2) > 0.0_r8 ) THEN
                capac0(i,2)=capac0(i,2) &
                     +epsflt*(capac(i,2)+capacm(i,2)-2.0_r8*capac0(i,2))
             END IF
          END DO
          DO i = 1, nmax
             tdm   (i)  =td0   (i)
             tgm   (i)  =tg0   (i)
             tcm   (i)  =tc0   (i)
             tmm   (i)  =tm0   (i)
             qmm   (i)  =qm0   (i)
             wm    (i,1)=w0    (i,1)
             wm    (i,2)=w0    (i,2)
             wm    (i,3)=w0    (i,3)
             capacm(i,1)=capac0(i,1)
             capacm(i,2)=capac0(i,2)
             IF (capacm(i,2) > 0.0_r8) tgm(i)=MIN(tgm(i),273.06_r8)
          END DO
          DO i = 1, nmax
             td0   (i)  =td    (i)
             tg0   (i)  =tg    (i)
             tc0   (i)  =tc    (i)
             tm0   (i)  =tm    (i)
             qm0   (i)  =qm    (i)
             w0    (i,1)=w     (i,1)
             w0    (i,2)=w     (i,2)
             w0    (i,3)=w     (i,3)
             capac0(i,1)=capac (i,1)
             capac0(i,2)=capac (i,2)
             IF (capac0(i,2) > 0.0_r8 .AND. tg0(i) > 273.16_r8) THEN
                nc=0
                !DO ii = 1, ncols
                !   IF (imask(ii) >= 1) nc=nc+1
                !   IF (nc == i) EXIT
                !END DO
                WRITE(UNIT=nfprt,FMT=200)ii,latitu,i,capac0(i,2),tg0(i)
             END IF
          END DO
       END IF
    ELSE
       DO i = 1, nmax
          tdm   (i)  =td   (i)
          tgm   (i)  =tg   (i)
          tcm   (i)  =tc   (i)
          tmm   (i)  =tm   (i)
          qmm   (i)  =qm   (i)
          wm    (i,1)=w    (i,1)
          wm    (i,2)=w    (i,2)
          wm    (i,3)=w    (i,3)
          capacm(i,1)=capac(i,1)
          capacm(i,2)=capac(i,2)
          IF (capacm(i,2) > 0.0_r8 .AND. tgm(i) > 273.16_r8) THEN
             nc=0
             !DO ii = 1, ncols
             !   IF (imask(ii) >= 1) nc=nc+1
             !   IF (nc == i) EXIT
             !END DO
             WRITE(UNIT=nfprt,FMT=650)ii,latitu,i,capacm(i,2),tgm(i)
          END IF
       END DO
       DO i = 1, nmax
          td0   (i)  =td   (i)
          tg0   (i)  =tg   (i)
          tc0   (i)  =tc   (i)
          tm0   (i)  =tm   (i)
          qm0   (i)  =qm   (i)
          w0    (i,1)=w    (i,1)
          w0    (i,2)=w    (i,2)
          w0    (i,3)=w    (i,3)
          capac0(i,1)=capac(i,1)
          capac0(i,2)=capac(i,2)
       END DO
    END IF
200 FORMAT(' CAPAC0 AND TG0 NOT CONSISTENT AT I,J,IS=',3I4, &
         ' CAPAC=',G16.8,' TG=',G16.8)
650 FORMAT(' CAPACM AND TGM NOT CONSISTENT AT I,J,IS=',3I4, &
         ' CAPAC=',G16.8,' TG=',G16.8)
  END SUBROUTINE sextrp


  SUBROUTINE Albedo( &
           ! Model information
            ncols     ,kMax      ,latco     ,&
            nmax      ,nsx       ,itype     ,&
            imask     , &
           ! Model Geometry
            cosz      ,zenith    , &
           ! Time info
            month2    ,month     , &
           ! Microphysics
            taud      , &
           ! Atmospheric fields
            wind      ,tsea      , &
           ! LW Radiation fields at last integer hour
            LwSfcDown ,&
           ! Radiation field (Interpolated) at time = tod
            xVisBeam  ,xVisDiff  ,xNirBeam   , &
            xNirDiff  , &
           ! Surface Albedo
            avisb     ,avisd     ,anirb      , &
            anird     )

   IMPLICIT NONE
   ! Model information
   INTEGER         , INTENT(IN   ) :: ncols
   INTEGER         , INTENT(IN   ) :: kmax
   INTEGER         , INTENT(IN   ) :: latco
   INTEGER         , INTENT(IN   ) :: nmax
   INTEGER         , INTENT(IN   ) :: nsx   (ncols)
   INTEGER         , INTENT(IN   ) :: itype (ncols)
   INTEGER(KIND=i8), INTENT(IN   ) :: imask (ncols)
   ! Model Geometry
   REAL   (KIND=r8), INTENT(IN   ) :: cosz  (ncols)
   REAL   (KIND=r8), INTENT(IN   ) :: zenith(ncols)
   ! Time info
   INTEGER         , INTENT(INOUT) :: month2(ncols)
   INTEGER         , INTENT(IN   ) :: month (ncols)
   ! Atmospheric fields
   REAL(KIND=r8) ,   INTENT(IN   ) :: wind  (ncols)!wind speed in m/s
   REAL(KIND=r8) ,   INTENT(IN   ) :: tsea  (ncols)
   ! Microphysics
   REAL(KIND=r8) ,   INTENT(IN   ) :: taud  (ncols,kMax)
   ! LW Radiation fields at last integer hour
   REAL(KIND=r8) ,   INTENT(IN   ) :: LwSfcDown(1:nCols)
   ! Radiation field (Interpolated) at time = tod
   REAL(KIND=r8) ,   INTENT(IN   ) :: xVisBeam (1:nCols)
   REAL(KIND=r8) ,   INTENT(IN   ) :: xVisDiff (1:nCols)
   REAL(KIND=r8) ,   INTENT(IN   ) :: xNirBeam (1:nCols)
   REAL(KIND=r8) ,   INTENT(IN   ) :: xNirDiff (1:nCols)
   ! Surface Albedo
   REAL(KIND=r8) ,   INTENT(OUT  ) :: avisb (ncols)
   REAL(KIND=r8) ,   INTENT(OUT  ) :: avisd (ncols)
   REAL(KIND=r8) ,   INTENT(OUT  ) :: anirb (ncols)
   REAL(KIND=r8) ,   INTENT(OUT  ) :: anird (ncols)

   REAL(KIND=r8) :: tc        (ncols)
   REAL(KIND=r8) :: tg    (ncols)
   REAL(KIND=r8) :: tm    (ncols)
   REAL(KIND=r8) :: qm    (ncols)
   REAL(KIND=r8) :: td    (ncols)
   REAL(KIND=r8) :: capac (ncols,2)
   REAL(KIND=r8) :: w          (ncols,3)
   REAL(KIND=r8) :: satcap(ncols,icg)
   REAL(KIND=r8) :: extk  (ncols,icg,iwv,ibd)
   REAL(KIND=r8) :: radfac(ncols,icg,iwv,ibd)
   REAL(KIND=r8) :: closs (ncols)
   REAL(KIND=r8) :: gloss (ncols)
   REAL(KIND=r8) :: thermk(ncols)
   REAL(KIND=r8) :: p1f        (ncols)
   REAL(KIND=r8) :: p2f        (ncols)
   REAL(KIND=r8) :: zlwup_local (ncols)
   REAL(KIND=r8) :: salb  (ncols,2,2)
   REAL(KIND=r8) :: tgeff (ncols)
   INTEGER       :: i
   INTEGER       :: k,m,n,ll
   REAL(KIND=r8) :: ocealb    
   REAL(KIND=r8) :: IceOceanAlb (nCols,2,2)  

!     --------------------------- INPUT ---------------------------------------
!   
!    specify the parameters for albedo here:

    REAL(KIND=r8) ::         tau(ncols)             !aerosol/cloud optical depth
    REAL(KIND=r8) ::         chl (ncols)            !chlorophyll concentration in mg/m3
    

   REAL(KIND=r8) :: f
   INTEGER       :: ncount

   IF(nmax.GE.1) THEN

       DO i=1,nmax
          tm      (i)    = tmm   (i,latco)
          qm      (i)    = qmm   (i,latco)
          td      (i)    = tdm   (i,latco)
          tg      (i)    = tgm   (i,latco)
          tc      (i)    = tcm   (i,latco)
          capac   (i,1)  = capacm(i,1,latco)
          capac   (i,2)  = capacm(i,2,latco)
          w       (i,1)  = wm    (i,1,latco)
          w       (i,2)  = wm    (i,2,latco)
          w       (i,3)  = wm    (i,3,latco)
          closs   (i)    = closs_gbl(i,latco)
          gloss   (i)    = gloss_gbl(i,latco)
          thermk  (i)    = thermk_gbl(i,latco)
          p1f     (i)    = p1f_gbl   (i,latco)
          p2f     (i)    = p2f_gbl   (i,latco)
          zlwup_local   (i)    = zlwup_SSiB(i,latco)
          salb  (i,1,1) = AlbGblSSiB(i,1,1,latco)
          salb  (i,1,2) = AlbGblSSiB(i,1,2,latco)
          salb  (i,2,1) = AlbGblSSiB(i,2,1,latco)
          salb  (i,2,2) = AlbGblSSiB(i,2,2,latco)
          tgeff   (i)    = tgeff_gbl (i,latco)
       END DO
       DO k=1,icg
          DO i=1,nmax
             satcap(i,k) = satcap3d(i,k,latco)
          END DO
       END DO
       DO m=1,ibd
          DO n=1,iwv
             DO ll=1,icg
                DO i=1,nmax
                   extk  (i,ll,n,m) = extk_gbl  (i,ll,n,m,latco)
                   radfac(i,ll,n,m) = radfac_gbl(i,ll,n,m,latco)
                END DO
             END DO
          END DO
       END DO    

       CALL radalb ( &
            nmax              ,month2(1:nmax)      ,nmax                ,itype(1:nmax)       , &
            tc(1:nmax)        ,tg(1:nmax)          ,capac(1:nmax,:)     ,satcap(1:nmax,:)    , &
            extk(1:nmax,:,:,:),radfac(1:nmax,:,:,:),closs(1:nmax)       ,gloss(1:nmax)       , &
            thermk(1:nmax)    ,p1f(1:nmax)         ,p2f(1:nmax)         ,zlwup_local(1:nmax) , &
            salb(1:nmax,:,:)  ,tgeff(1:nmax)       ,cosz(1:nmax)        ,nsx(1:nmax)         , &
            latco     )
    DO k=1,icg
       DO i=1,nmax
          satcap3d(i,k,latco) =satcap(i,k) 
       END DO
    END DO
    DO i=1,nmax
       closs_gbl (i,latco) = closs(i) 
       gloss_gbl (i,latco) = gloss(i)
       thermk_gbl(i,latco) = thermk(i)
       p1f_gbl   (i,latco) = p1f  (i)
       p2f_gbl   (i,latco) = p2f  (i)
       zlwup_SSiB(i,latco)  =  zlwup_local(i)  
       AlbGblSSiB(i,1,1,latco)  = salb  (i,1,1)
       AlbGblSSiB(i,1,2,latco)  = salb  (i,1,2)
       AlbGblSSiB(i,2,1,latco)  = salb  (i,2,1)
       AlbGblSSiB(i,2,2,latco)  = salb  (i,2,2)
       tgeff_gbl(i,latco)  = tgeff (i)   
    END DO
    
    
    DO m=1,ibd
       DO n=1,iwv
          DO ll=1,icg
             DO i=1,nmax
                extk_gbl  (i,ll,n,m,latco) = extk  (i,ll,n,m)
                radfac_gbl(i,ll,n,m,latco) = radfac(i,ll,n,m) 
             END DO
          END DO
       END DO
    END DO    

   END IF
!     --------------------------- INPUT ---------------------------------------
!   
!    specify the parameters for albedo here:
    IceOceanAlb=0.0_r8 
    chl = 0.10_r8              !chlorophyll concentration in mg/m3
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !
    tau=0.0_r8
    DO k=1,kMax
       DO i=1,ncols
          tau(i)=tau(i)+ taud(i,k) ! tau = SUM(taud(i,1:kMax))              !aerosol/cloud optical depth
       END DO
    END DO   
    DO i=1,ncols
        tau(i)=MIN(MAX(tau(i),0.0_r8),25.0_r8)      !aerosol/cloud optical depth
    END DO

    
    
    IF(TRIM(SLABOCEAN) == 'SLAB')THEN
       ncount=0
       DO i=1,ncols
          IF(imask(i).GE.1_i8) THEN
                 ncount=ncount+1
             avisb(i)=salb(ncount,1,1)
             avisd(i)=salb(ncount,1,2)
             anirb(i)=salb(ncount,2,1)
             anird(i)=salb(ncount,2,2)
          ELSE IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
             f=MAX(zenith(i),0.0e0_r8 )
             avisb(i)=GetOceanAlb(i,tau(i),f,wind(i),chl(i),0.2_r8,0.7_r8) !   s (shortwave)  = 0.2-0.7 micro-meters
             avisd(i)=oceald
             anirb(i)=GetOceanAlb(i,tau(i),f,wind(i),chl(i),0.7_r8,5.0_r8) !   l (longwave)   = 0.7-5.0 micro-meters
             anird(i)=oceald
          ELSE
             IF (TRIM(ICEMODEL)=='SSIB')THEN
                f=MAX(zenith(i),0.0e0_r8 )
                IceOceanAlb(i,:,:)=GetIceOceanAlb(i,latco,month(i),xVisBeam(i),xVisDiff(i),&
                                       xNirBeam(i),xNirDiff(i),f,LwSfcDown(i))
                avisb(i)=IceOceanAlb(i,1,1)
                avisd(i)=IceOceanAlb(i,1,2)
                anirb(i)=IceOceanAlb(i,2,1)
                anird(i)=IceOceanAlb(i,2,2)
             ELSE IF (TRIM(ICEMODEL)=='COLA')THEN
                avisb(i)=icealv
                avisd(i)=icealv
                anirb(i)=icealn
                anird(i)=icealn
             ELSE
                STOP "ICEMODEL ->OPTIONS"
             END IF      
          END IF
       END DO
    ELSE IF(TRIM(SLABOCEAN) == 'COLA')THEN
       ncount=0
       DO i=1,ncols
          IF(imask(i).GE.1_i8) THEN
             ncount=ncount+1
             avisb(i)=salb(ncount,1,1)
             avisd(i)=salb(ncount,1,2)
             anirb(i)=salb(ncount,2,1)
             anird(i)=salb(ncount,2,2)
          ELSE IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
             f=MAX(zenith(i),0.0e0_r8 )
             ocealb=0.12347e0_r8 +f*(0.34667e0_r8+f*(-1.7485e0_r8 + &
                  f*(2.04630e0_r8 -0.74839e0_r8 *f)))
             avisb(i)=ocealb
             avisd(i)=oceald
             anirb(i)=ocealb
             anird(i)=oceald
          ELSE
             IF (TRIM(ICEMODEL)=='SSIB')THEN
                f=MAX(zenith(i),0.0e0_r8 )
                IceOceanAlb(i,:,:)=GetIceOceanAlb(i,latco,month(i),xVisBeam(i),xVisDiff(i),&
                                       xNirBeam(i),xNirDiff(i),f,LwSfcDown(i))
                avisb(i)=IceOceanAlb(i,1,1)
                avisd(i)=IceOceanAlb(i,1,2)
                anirb(i)=IceOceanAlb(i,2,1)
                anird(i)=IceOceanAlb(i,2,2)
             ELSE IF (TRIM(ICEMODEL)=='COLA')THEN
                avisb(i)=icealv
                avisd(i)=icealv
                anirb(i)=icealn
                anird(i)=icealn
             ELSE
                STOP "ICEMODEL ->OPTIONS"
             END IF      
          END IF
       END DO
    ELSE
       WRITE(0,*)"ERRO SLABOCEAN",TRIM(SLABOCEAN)
       STOP 
    END IF
  END SUBROUTINE Albedo
  
  ! radalb :surface albedos via two stream approximation (direct and diffuse).



  SUBROUTINE radalb ( &
            ncols             ,mon                 ,nmax                ,itype               , &
            tc                ,tg                  ,capac               ,satcap              , &
            extk              ,radfac              ,closs               ,gloss               , &
            thermk            ,p1f                 ,p2f                 ,zlwup               , &
            salb              ,tgeff               ,cosz                ,nsx                 , &
            latitu )
    !
    !
    ! reference  : a simple biosphere model (xue et al 1991)
    !-----------------------------------------------------------------------
    !     *** indices ***
    !   cg =1...canopy
    !   cg =2...ground cover
    !   vn =1...visible      (0.0-0.7 micron)
    !   vn =2...near-infrared(0.7-3.0 micron)
    !   bd =1...beam
    !   bd =2...diffuse
    !   ld =1...live leaves
    !   ld =2...dead leaves
    !   vnt=1...visible      (0.0-0.7 micron)
    !   vnt=2...near-infrared(0.7-3.0 micron)
    !   vnt=3...thermal
    !-----------------------------------------------------------------------
    !        input parameters
    !-----------------------------------------------------------------------
    !   zlt(cg)..........leaf area index
    !   z1...............bottom height of canopy
    !   z2...............top    height of canopy
    !   ref (cg,vnt,ld)..reflectance   of vegetation
    !   tran(cg,vnt,ld)..transmittance of vegetation
    !   green (cg).......fraction of green leaf area
    !   chil  (cg).......leaf orientation factor
    !   vcover(cg).......fraction of vegetation cover
    !   soref (vnt)......ground albedo
    !   chil  (cg).......leaf orientation factor
    !   cosz.............cosine of solar zenith angle
    !   tf...............water freezing temperature
    !   tg...............ground temperature
    !   tc...............canopy leaf temperature
    !   capac(cg)........water store capacity of leaves
    !   stefan...........stefan-boltsman constant
    !-----------------------------------------------------------------------
    !     in-subr. parameters
    !-----------------------------------------------------------------------
    !   albedo(cg,vnt,bd)
    !-----------------------------------------------------------------------
    !       output parameters
    !-----------------------------------------------------------------------
    !   extk(cg,vnt,bd)..extinction coefficient
    !                    passed to subr.raduse through radsave
    !   radfac(cg,vn,bd).fractions of downward solar radiation at surface
    !                    passed to subr.raduse
    !   salb(vn,bd)......surface albedo
    !                    passed to subr.spmrad
    !   tgeff............effective ground temperature
    !                    passed to subr.spmrad
    !   thermk...........canopy emissivity
    !   radsav(1)........beam    extinction coefficient (par)
    !   radsav(2)........diffuse extinction coefficient (par)
    !   closs............radiation loss from canopy
    !   gloss............radiation loss from ground
    !-----------------------------------------------------------------------
    !
    !   ityp.......Numero das classes de solo 13
    !   imon.......Numero maximo de meses no ano (12)
    !   icg........Parametros da vegetacao (icg=1 topo e icg=2 base)
    !   iwv........Compriment de onda iwv=1=visivel, iwv=2=infravermelho
    !              proximo, iwv=3 infravermelho termal
    !   ibd........Estado da vegetacao ibd=1 verde / ibd=2 seco
    !   ncols......Number of grid points on a gaussian latitude circle
    !   mon........Number of month at year (1-12)
    !   nmax
    !   itype......Classe de textura do solo
    !   satcap.....saturation liquid water capacity         (m)
    !   p1f........
    !   p2f........
    !   zlwup......zlwup(i)= stefan*( fac1(i)*tc4(i)+ &
    !              (1.0  -vcover(i,1)*(1.0  -thermk(i)))*fac2(i)*tg4(i))
    !   nsx........
    !

    INTEGER, INTENT(IN   ) :: ncols
    INTEGER, INTENT(INOUT) :: mon(ncols)
    INTEGER, INTENT(IN   ) :: nmax
    INTEGER, INTENT(IN   ) :: itype (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: tc    (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: tg    (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: capac (ncols,2)
    REAL(KIND=r8),    INTENT(INOUT  ) :: satcap(ncols,icg)
    REAL(KIND=r8),    INTENT(INOUT  ) :: extk  (ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(INOUT  ) :: radfac(ncols,icg,iwv,ibd)
    REAL(KIND=r8),    INTENT(INOUT  ) :: closs (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: gloss (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: thermk(ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: p1f   (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: p2f   (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: zlwup (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: salb  (ncols,2,2)
    REAL(KIND=r8),    INTENT(INOUT  ) :: tgeff (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: cosz  (ncols)
    INTEGER, INTENT(IN   ) :: nsx (ncols)
    INTEGER, INTENT(IN   ) :: latitu

    REAL(KIND=r8)   :: zlt2    (ncols,icg)
    REAL(KIND=r8)   :: vcover  (ncols,icg)


    REAL(KIND=r8) :: f     (ncols)
    REAL(KIND=r8) :: deltg (ncols)
    REAL(KIND=r8) :: fmelt (ncols)
    REAL(KIND=r8) :: depcov(ncols)
    REAL(KIND=r8) :: scov  (ncols)
    REAL(KIND=r8) :: scov2 (ncols)
    REAL(KIND=r8) :: tc4   (ncols)
    REAL(KIND=r8) :: tg4   (ncols)
    REAL(KIND=r8) :: fac1  (ncols)
    REAL(KIND=r8) :: fac2  (ncols)
    REAL(KIND=r8) :: zkat  (ncols)

    INTEGER, PARAMETER :: nk=3
    REAL(KIND=r8)    :: temp(nmax,18)
    REAL(KIND=r8)    :: xmi1(nmax,nk)
    INTEGER :: i
    INTEGER :: ntyp(ncols)
    INTEGER :: monx(ncols)
    INTEGER :: jj
    INTEGER :: i1
    INTEGER :: ml(nmax)
    INTEGER :: k1
    INTEGER :: k2
    INTEGER :: ik
    REAL(KIND=r8)    :: capaci
    !    REAL(KIND=r8)    :: xf
    !    REAL(KIND=r8)    :: xf2
    !    REAL(KIND=r8)    :: sc1
    !    REAL(KIND=r8)    :: sc2
    REAL(KIND=r8)    :: xm1
    !    REAL(KIND=r8)    :: xm2
    REAL(KIND=r8)    :: xtm1
    REAL(KIND=r8)    :: xtm2
    REAL(KIND=r8)    :: stbi
    LOGICAL :: flagtyp(nmax)
    LOGICAL :: flagscov(nmax)

    DO i = 1, nmax

       zlt2       (i,1)   =  zlt_gbl    (i,latitu,1) !zlt   (itype(i),mon(i),1)
       zlt2       (i,2)   =  zlt_gbl    (i,latitu,2) !zlt   (itype(i),mon(i),2)
       vcover     (i,1)   =  vcover_gbl (i,latitu,1) !xcover(itype(i),mon(i),1)
       vcover     (i,2)   =  vcover_gbl (i,latitu,2) !xcover(itype(i),mon(i),2)
       f(i)= MAX ( cosz(i), 0.01746_r8  )
    END DO
    !
    !     maximum water storage values.
    !
    DO i = 1, nmax
       deltg(i)=tf-tg(i)
       fmelt(i)=1.0_r8
       IF (ABS(deltg(i)) < 0.5_r8 .AND. deltg(i) > 0.0_r8) THEN
          fmelt(i)=0.6_r8
       END IF
    END DO
    ntyp=itype
    DO i = 1, nmax
       !ntyp=itype(i)
       satcap(i,1)=zlt2(i,1)*1.0e-4_r8
       satcap(i,2)=zlt2(i,2)*1.0e-4_r8
       depcov(i  )=MAX(0.0_r8  ,capac(i,2)*5.0_r8  -z1(ntyp(i),mon(i)))
       depcov(i  )=MIN(depcov(i),(z2(ntyp(i),mon(i))-z1(ntyp(i),mon(i)))*0.95_r8  )
       satcap(i,1)=satcap(i,1) &
            *(1.0_r8  -depcov(i)/(z2(ntyp(i),mon(i))-z1(ntyp(i),mon(i))))
    END DO

    DO i = 1, nmax
       scov(i)=0.0_r8
       IF (tc(i) <= tf) THEN
          scov(i)= MIN( 0.5_r8  , capac(i,1)/satcap(i,1))
       END IF
    END DO
    capaci=1.0_r8  /0.004_r8
    DO i = 1, nmax
       IF (tg(i) > tf) THEN
          scov2(i)=0.0_r8
       ELSE
          scov2(i)=MIN( 1.0_r8  , capac(i,2)*capaci)
       END IF
    END DO
    !
    !     terms which multiply incoming short wave fluxes
    !     to give absorption of radiation by canopy and ground
    !
    monx = mon
    DO i = 1, nmax
       IF (fmelt(i) == 1.0_r8) THEN
          ml(i) = 1
       ELSE
          ml(i) = 2
       END IF
    END DO
    ntyp=itype

    DO i = 1, nmax
       mon(i) = monx(i)
       flagtyp(i) = .TRUE.
       IF (ntyp(i) == 13) ntyp(i) = 11
       IF (ntyp(i) == 12 .AND. nsx(i) > 0) THEN
          ntyp(i) = 13
          mon(i) = nsx(i)
          IF (nsx(i) == 1 .AND. (monx(i) >= 9 .AND. monx(i) <= 11)) mon(i) = 7
          flagtyp(i) = .FALSE.
       END IF
    END DO
    DO jj = 1, nk
       DO i=1, nmax
          xmi1(i,jj) = xmiu(mon(i),jj)
       END DO
    END DO
    DO jj = 1, nk
       DO i=1, nmax
          IF (.NOT.flagtyp(i))xmi1(i,jj) = xmiw(mon(i),jj)
       END DO
    END DO
    !
    !        snow free case
    !
    DO i = 1, nmax
       flagscov(i) = scov(i) < 0.025_r8 .AND. scov2(i) < 0.025_r8
    END DO

    DO i1 = 1, 9
       DO i = 1, nmax
          IF (flagscov(i)) THEN
             temp(i,i1) = cledir(ntyp(i),mon(i),i1,1) + cledir(ntyp(i),mon(i),i1,2) &
                  * f(i) + cledir(ntyp(i),mon(i),i1,3) * (f(i)*f(i))
             temp(i,i1+9) = cledfu(ntyp(i),mon(i),i1)
          END IF
       END DO
    END DO
    flagscov = .NOT. flagscov
    DO i1 = 1, 9
       DO i = 1, nmax
          IF (flagscov(i)) THEN
             !
             !     with snow cover
             !
             temp(i,i1) = cedir(ntyp(i),mon(i),i1,1) + f(i) * &
                  cedir(ntyp(i),mon(i),i1,2) + cedir(ntyp(i),mon(i),i1,3) * (f(i)*f(i))
             temp(i,i1+9) = cedfu(ntyp(i),mon(i),i1)
          END IF
       END DO
    END DO
    DO i1 = 1, 6
       DO i = 1, nmax
          IF (flagscov(i) .AND. ntyp(i) == 11) THEN
             !sc2 = scov2(i) * scov2(i)
             !sc1 = scov2(i)
             temp(i,i1)=cedir2(ml(i),ntyp(i),mon(i),i1,nk,1)+ &
                  cedir2(ml(i),ntyp(i),mon(i),i1,nk,2) &
                  *scov2(i) + cedir2(ml(i),ntyp(i),mon(i),i1,nk,3) *(scov2(i) * scov2(i)) + temp(i,i1)
             temp(i,i1+9) = cedfu2(ml(i),ntyp(i),mon(i),i1,1) +  &
                  cedfu2(ml(i),ntyp(i),mon(i),i1,2) &
                  * scov2(i) + cedfu2(ml(i),ntyp(i),mon(i),i1,3) * (scov2(i) * scov2(i)) + temp(i,i1+9)
          END IF
       END DO
    END DO
    DO i = 1, nmax
       IF (flagscov(i) .AND. ntyp(i) /= 11) THEN
          k2 = 1
          k1 = 2
          DO ik = nk, 1, -1
             IF (f(i) >= xmi1(i,ik)) THEN
                CONTINUE
             ELSE
                k1 = ik + 1
                k2 = ik
                EXIT
             END IF
          END DO
          !xm2 = xmi1(mon(i),k2)
          IF (k1 <= nk) xm1 = xmi1(i,k1)
          !
          !     snow cover at 1st layer
          !
          IF (scov(i) > 0.025_r8) THEN
             !sc2 = scov(i) * scov(i)
             !sc1 = scov(i)
             IF (k2 >= nk .OR. k2 <= 1) THEN
                DO i1 = 1, 6
                   temp(i,i1)=cedir1(ml(i),ntyp(i),mon(i),i1,k2,1)+ &
                        cedir1(ml(i),ntyp(i),mon(i),i1,k2,2)*scov(i) + &
                        cedir1(ml(i),ntyp(i),mon(i),i1,k2,3)*(scov(i) * scov(i)) + temp(i,i1)
                END DO
             ELSE
                DO i1 = 1, 6
                   xtm1=cedir1(ml(i),ntyp(i),mon(i),i1,k1,1)+ &
                        cedir1(ml(i),ntyp(i),mon(i),i1,k1,2)*scov(i) + &
                        cedir1(ml(i),ntyp(i),mon(i),i1,k1,3)*(scov(i) * scov(i))
                   xtm2=cedir1(ml(i),ntyp(i),mon(i),i1,k2,1)+ &
                        cedir1(ml(i),ntyp(i),mon(i),i1,k2,2)*scov(i)+ &
                        cedir1(ml(i),ntyp(i),mon(i),i1,k2,3) *(scov(i) * scov(i))
                   temp(i,i1) = (xtm1*((xmi1(i,k2))-f(i))+xtm2*(f(i)-xm1))/((xmi1(i,k2))-xm1) &
                        + temp(i,i1)
                END DO
             END IF
             DO i1 = 1, 6
                temp(i,i1+9) = cedfu1(ml(i),ntyp(i),mon(i),i1,1) +  &
                     cedfu1(ml(i),ntyp(i),mon(i),i1,2)*scov(i) + &
                     cedfu1(ml(i),ntyp(i),mon(i),i1,3) * (scov(i) * scov(i)) + temp(i,i1+9)
             END DO
          END IF
          !
          !     snow cover on ground
          !
          IF (scov2(i) > 0.025_r8) THEN
             !sc2 = scov2(i) * scov2(i)
             !sc1 = scov2(i)
             IF (k2 >= nk .OR. k2 <= 1) THEN
                DO i1 = 1, 6
                   temp(i,i1)=cedir2(ml(i),ntyp(i),mon(i),i1,k2,1)+ &
                        cedir2(ml(i),ntyp(i),mon(i),i1,k2,2)*scov2(i) +  &
                        cedir2(ml(i),ntyp(i),mon(i),i1,k2,3)*(scov2(i) * scov2(i)) + temp(i,i1)
                END DO
             ELSE
                DO i1 = 1, 6
                   xtm1=cedir2(ml(i),ntyp(i),mon(i),i1,k1,1)+ &
                        cedir2(ml(i),ntyp(i),mon(i),i1,k1,2)*scov2(i) + &
                        cedir2(ml(i),ntyp(i),mon(i),i1,k1,3) *(scov2(i) * scov2(i))
                   xtm2=cedir2(ml(i),ntyp(i),mon(i),i1,k2,1)+ &
                        cedir2(ml(i),ntyp(i),mon(i),i1,k2,2)*scov2(i)+ &
                        cedir2(ml(i),ntyp(i),mon(i),i1,k2,3) *(scov2(i) * scov2(i))
                   temp(i,i1) = (xtm1*((xmi1(i,k2))-f(i))+xtm2*(f(i)-xm1))/((xmi1(i,k2))-xm1) &
                        + temp(i,i1)
                END DO
             END IF
             DO i1 = 1, 6
                temp(i,i1+9) = cedfu2(ml(i),ntyp(i),mon(i),i1,1) +  &
                     cedfu2(ml(i),ntyp(i),mon(i),i1,2)* scov2(i) +  &
                     cedfu2(ml(i),ntyp(i),mon(i),i1,3) * (scov2(i) * scov2(i)) + temp(i,i1+9)
             END DO
          END IF
       END IF
    END DO
    !500    CONTINUE
    DO i = 1, nmax
       radfac(i,1,1,2) = temp(i,10)
       radfac(i,1,2,2) = temp(i,11)
       radfac(i,2,1,2) = temp(i,12)
       radfac(i,2,2,2) = temp(i,13)
       salb(i,1,2) = temp(i,14)
       salb(i,2,2) = temp(i,15)
       p2f(i) =  temp(i,16)
       extk(i,1,1,2) = temp(i,17)
       extk(i,2,1,2) = temp(i,18)
       radfac(i,1,1,1) = temp(i,1)
       radfac(i,1,2,1) = temp(i,2)
       radfac(i,2,1,1) = temp(i,3)
       radfac(i,2,2,1) = temp(i,4)
       salb(i,1,1) = temp(i,5)
       salb(i,2,1) = temp(i,6)
       p1f(i) =  temp(i,7)
       extk(i,1,1,1) = temp(i,8) / f(i)
       extk(i,2,1,1) = temp(i,9) / f(i)
       extk(i,1,3,1) = cether(ntyp(i),mon(i),1)
       extk(i,1,3,2) = cether(ntyp(i),mon(i),2)
       extk(i,2,3,1) = cether(ntyp(i),mon(i),1)
       extk(i,2,3,2) = cether(ntyp(i),mon(i),2)
    END DO
    mon = monx
    !
    !     long-wave flux terms from canopy and ground
    !
    stbi=1.0_r8  /stefan
    DO  i = 1, nmax
       tc4(i)=tc(i)*tc(i)*tc(i)*tc(i)
       tg4(i)=tg(i)*tg(i)*tg(i)*tg(i)
       !ntyp=itype(i)
       zkat(i)=extk(i,1,3,2)*zlt2(i,1)/vcover(i,1)
       zkat(i)=MAX(expcut  ,-zkat(i) )
       zkat(i)=MIN(-10.0e-5_r8, zkat(i) )
       thermk(i)=EXP(zkat(i))
       fac1 (i)=vcover(i,1)*( 1.0_r8  -thermk(i) )
       fac2 (i)=1.0_r8
       closs(i)=2.0_r8  *fac1(i)*stefan*tc4(i)
       closs(i)=closs(i)-fac2(i)*fac1(i)*stefan*tg4(i)
       gloss(i)= fac2(i)*stefan*tg4(i)
       gloss(i)= gloss(i)-fac1(i)*fac2(i)*stefan*tc4(i)
       !
       !     effective surface radiative temperature ( tgeff )
       !
       zlwup(i) = stefan*( fac1(i)*tc4(i) + (1.0_r8  - vcover(i,1) * (1.0_r8  -thermk(i)))*fac2(i)*tg4(i))
       tgeff(i)=SQRT ( SQRT (( zlwup(i)*stbi )))
    END DO
  END SUBROUTINE radalb


  !
  !
  !
  SUBROUTINE Phenology(latco,nCols,nmax,itype,colrad2, month2, xday, idatec,nsx)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: latco
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: nmax
    INTEGER      , INTENT(in )   :: itype   (ncols)
    REAL(KIND=r8), INTENT(in )   :: colrad2 (ncols)
    INTEGER      , INTENT(in )   :: month2  (ncols)
    REAL(KIND=r8), INTENT(in )   :: xday
    INTEGER      , INTENT(in )   :: idatec(4)
    INTEGER      , INTENT(inout) :: nsx(ncols)
    IF(schemes==1)CALL CopySurfaceData(itype,month2,colrad2,xday,idatec,nsx,nCols,nmax,latco)
    IF(schemes==2) PRINT*,'ERROR schemes 2'
  END SUBROUTINE Phenology



  ! vegin  :reads vegetation morphoLOGICAL and physioLOGICAL data.




  SUBROUTINE vegin ( )

    REAL(KIND=r8), PARAMETER :: si1=1
    REAL(KIND=r8), PARAMETER :: sl1=0.8888
    INTEGER, PARAMETER ::  njj=6,nj=9, nk=3,ild=2

   ! Vegetation and Soil Parameters

   REAL (KIND=r4) rstpar_r4(ityp,icg,iwv), &
                  chil_r4(ityp,icg), &
                  topt_r4(ityp,icg), &
                  tll_r4(ityp,icg), &
                  tu_r4(ityp,icg), &
                  defac_r4(ityp,icg), &
                  ph1_r4(ityp,icg), &
                  ph2_r4(ityp,icg), &
                  rootd_r4(ityp,icg), &
                  bee_r4(ityp), &
                  phsat_r4(ityp), &
                  satco_r4(ityp), &
                  poros_r4(ityp), &
                  zdepth_r4(ityp,idp), &
                  green_r4(ityp,imon,icg), &
                  xcover_r4(ityp,imon,icg), &
                  zlt_r4(ityp,imon,icg), &
                  x0x_r4(ityp,imon),&
                  xd_r4(ityp,imon), &
                  z2_r4   (ityp,imon), &
                  z1_r4   (ityp,imon), &
                  xdc_r4  (ityp,imon), &
                  xbc_r4  (ityp,imon)
                  
    REAL(KIND=r4) :: cedfu_r4 (ityp,imon,nj), &
                     cedir_r4 (ityp,imon,nj,3), &
                     cedfu1_r4(2,ityp,imon,njj,3), &
                     cedir1_r4(2,ityp,imon,njj,nk,3), &
                     cedfu2_r4(2,ityp,imon,njj,3), &
                     cedir2_r4(2,ityp,imon,njj,nk,3), &
                     cledfu_r4(ityp,imon,nj), &
                     cledir_r4(ityp,imon,nj,3), &
                     cether_r4(ityp,imon,2), &
                     xmiu_r4  (imon,nk), &
                     xmiw_r4  (imon,nk)

    INTEGER :: jcg
    INTEGER :: jmon
    INTEGER :: jtyp
    INTEGER :: iv
    INTEGER :: im
    INTEGER :: i
    REAL(KIND=r8)    :: f0001
    REAL(KIND=r8)    :: yhil (2)
    REAL(KIND=r8)    :: dz
    REAL(KIND=r8)    :: dzcut
    REAL(KIND=r8)    :: tvsgm
    INTEGER :: ierr
    !
    ALLOCATE(cedfu (13,12, 9)           )
    ALLOCATE(cedir (13,12, 9,3)      )
    ALLOCATE(cedfu1( 2,13,12,6,3)    )
    ALLOCATE(cedir1( 2,13,12,6,3,3)  )
    ALLOCATE(cedfu2( 2,13,12,6,3)    )
    ALLOCATE(cedir2( 2,13,12,6,3,3)  )
    ALLOCATE(cledfu(13,12, 9)           )
    ALLOCATE(cledir(13,12, 9,3)      )
    ALLOCATE(xmiu  (12, 3)           )
    ALLOCATE(cether(13,12, 2)           )
    ALLOCATE(xmiw  (12, 3)           )
    !
    ALLOCATE(ystpar(2,3)             )
    ALLOCATE(yopt  (2)               )
    ALLOCATE(yll   (2)               )
    ALLOCATE(yu    (2)               )
    ALLOCATE(yefac (2)               )
    ALLOCATE(yh1   (2)               )
    ALLOCATE(yh2   (2)               )
    ALLOCATE(yootd (2)               )
    ALLOCATE(yreen (12,2)            )
    ALLOCATE(ycover(12,2)            )
    ALLOCATE(ylt   (12,2)            )
    !
    !     vegetation and soil parameters
    !
    ALLOCATE(rstpar_fixed(ityp,icg,iwv)   )
    ALLOCATE(chil_fixed  (ityp,icg)          )
    ALLOCATE(topt_fixed  (ityp,icg)          )
    ALLOCATE(tll_fixed   (ityp,icg)          )
    ALLOCATE(tu_fixed    (ityp,icg)          )
    ALLOCATE(defac_fixed (ityp,icg)          )
    ALLOCATE(ph1_fixed   (ityp,icg)          )
    ALLOCATE(ph2_fixed   (ityp,icg)          )
    ALLOCATE(rootd (ityp,icg)          )
    ALLOCATE(bee   (ityp)          )
    ALLOCATE(phsat (ityp)          )
    ALLOCATE(satco (ityp)          )
    ALLOCATE(poros (ityp)          )
    ALLOCATE(zdepth(ityp,idp)          )
    ALLOCATE(green_fixed (ityp,imon,icg)  )
    ALLOCATE(xcover_fixed(ityp,imon,icg)  )
    ALLOCATE(zlt_fixed   (ityp,imon,icg)  )
    ALLOCATE(x0x   (ityp,imon)      )
    ALLOCATE(xd    (ityp,imon)      )
    ALLOCATE(z2    (ityp,imon)      )
    ALLOCATE(z1    (ityp,imon)      )
    ALLOCATE(xdc   (ityp,imon)      )
    ALLOCATE(xbc   (ityp,imon)      )
    ALLOCATE(zlt   (ityp,imon,icg)  )
    ALLOCATE(xcover  (ityp, imon, icg))
    ALLOCATE(ph2    (ityp,icg))
    ALLOCATE(ph1    (ityp,icg))
    ALLOCATE(green(ityp,imon,icg))
    ALLOCATE(defac(ityp,icg))
    ALLOCATE(tu   (ityp,icg))
    ALLOCATE(tll  (ityp,icg))
    ALLOCATE(topt (ityp,icg))
    ALLOCATE(rstpar(ityp,icg,iwv))
    ALLOCATE(chil  (ityp,icg))

    OPEN(UNIT=nfsibd, FILE=TRIM(fNameSibVeg),FORM='UNFORMATTED', ACCESS='SEQUENTIAL',&
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSibVeg), ierr
       STOP "**(ERROR)**"
    END IF

    OPEN (UNIT=nfalb, FILE=TRIM(fNameSibAlb),FORM='UNFORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ',STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameSibAlb), ierr
       STOP "**(ERROR)**"
    END IF
   
    READ (UNIT=nfsibd) rstpar_r4, chil_r4, topt_r4, tll_r4, tu_r4, defac_r4, ph1_r4, ph2_r4, &
                       rootd_r4, bee_r4, phsat_r4, satco_r4, poros_r4, zdepth_r4
    READ (UNIT=nfsibd) green_r4, xcover_r4, zlt_r4, x0x_r4, xd_r4, z2_r4, z1_r4, xdc_r4, xbc_r4

    rstpar_fixed = rstpar_r4
    chil_fixed   = chil_r4
    topt_fixed   = topt_r4
    tll_fixed    = tll_r4
    tu_fixed     = tu_r4
    defac_fixed  = defac_r4
    ph1_fixed    = ph1_r4
    ph2_fixed    = ph2_r4
    rootd        = rootd_r4
    bee          = bee_r4
    phsat        = phsat_r4
    satco        = satco_r4
    poros        = poros_r4
    zdepth       = zdepth_r4
    green_fixed  = green_r4
    xcover_fixed = xcover_r4
    zlt_fixed    = zlt_r4
    x0x          = x0x_r4
    xd           = xd_r4
    z2           = z2_r4
    z1           = z1_r4
    xdc          = xdc_r4
    xbc          = xbc_r4

    READ(UNIT=nfalb) cedfu_r4, cedir_r4, cedfu1_r4, cedir1_r4, cedfu2_r4, cedir2_r4, &
         cledfu_r4, cledir_r4, xmiu_r4, cether_r4, xmiw_r4
    cedfu  = REAL(cedfu_r4 ,KIND=r8) 
    cedir  = REAL(cedir_r4 ,KIND=r8) 
    cedfu1 = REAL(cedfu1_r4,KIND=r8) 
    cedir1 = REAL(cedir1_r4,KIND=r8) 
    cedfu2 = REAL(cedfu2_r4,KIND=r8) 
    cedir2 = REAL(cedir2_r4,KIND=r8) 
    cledfu = REAL(cledfu_r4,KIND=r8) 
    cledir = REAL(cledir_r4,KIND=r8) 
    cether = REAL(cether_r4,KIND=r8) 
    xmiu   = REAL(xmiu_r4  ,KIND=r8) 
    xmiw   = REAL(xmiw_r4  ,KIND=r8) 
    REWIND nfsibd

    REWIND nfalb

    f0001=0.0001_r8

    DO jcg =1, 2
       DO jmon=1,12
          DO jtyp=1,ityp
             green_fixed(jtyp,jmon,jcg)=MAX(f0001,green_fixed(jtyp,jmon,jcg))
          END DO
       END DO
    END DO

    DO iv =1, 2
       jtyp = 12
       IF (iv.EQ.2) jtyp = 13

       DO  im = 1,3
          ystpar(iv,im)=rstpar_fixed(jtyp,1,im)
       END DO

       yhil  (iv)=chil_fixed  (jtyp,1)
       yopt  (iv)=topt_fixed  (jtyp,1)
       yll   (iv)=tll_fixed   (jtyp,1)
       yu    (iv)=tu_fixed    (jtyp,1)
       yefac (iv)=defac_fixed (jtyp,1)
       yootd (iv)=rootd (jtyp,1)
       yh1   (iv)=ph1_fixed   (jtyp,1)
       yh2   (iv)=ph2_fixed   (jtyp,1)

    END DO

    DO jmon=1,12
       DO iv = 1,2
          jtyp = 12
          IF (iv.EQ.2) jtyp = 13
          ylt   (jmon,iv)=zlt_fixed(jtyp,jmon,1)
          yreen (jmon,iv)=green_fixed (jtyp,jmon,1)
          ycover(jmon,iv)=xcover_fixed(jtyp,jmon,1)
       END DO
    END DO

    DO iv = 1,2
       DO im = 1,3
          rstpar_fixed(13,iv,im) = 1000.0_r8
       END DO
       chil_fixed  (13,iv) = 0.01_r8
       topt_fixed  (13,iv) = 310.0_r8
       tll_fixed   (13,iv) = 300.0_r8
       tu_fixed    (13,iv) = 320.0_r8
       defac_fixed (13,iv) = 0.0_r8
       ph1_fixed   (13,iv) = 3.0_r8
       ph2_fixed   (13,iv) = 6.0_r8
       rootd       (13,iv) = 2.1_r8
    END DO

    bee(13) = 4.8_r8
    phsat(13) = -0.167_r8
    satco(13) = 0.762e-4_r8
    poros(13) = 0.4352_r8

    DO i = 1, imon
       zlt_fixed(13,i,1) = 0.0001_r8
       zlt_fixed(13,i,2) = 0.0001_r8
       z2(13,i) = 0.1_r8
       z1(13,i) = 0.0001_r8
       xcover_fixed(13,i,1) = 0.0001_r8
       xcover_fixed(13,i,2) = 0.0001_r8
       x0x(13,i) = 0.01_r8
       xd(13,i) = 0.0004_r8
       xbc(13,i) = 35461.0_r8
       xdc(13,i) = 28.5_r8
    END DO

    zdepth(13,1) = 1.0_r8
    zdepth(13,2) = 1.0_r8
    zdepth(13,3) = 1.0_r8

    ! tvsgm - Global Mean Surface Virtual Temperature
    ! dz - mean height of the first model layer
    tvsgm=288.16_r8
    dz=(gasr*tvsgm/grav)*LOG(si1/sl1)
    ! Forest
    !dzcut=0.75_r8*dz
    dzcut=0.6_r8*dz
    xd(1,1:imon)=MIN(xd(1,1:imon),dzcut)
    ! Other
    ! SiB calibration values
    ! 45 m - height of the first tower level of measurements
    ! 27 m - maximum calibrated displacement height
    dzcut=(27.0_r8/45.0_r8)*dz
    xd(2:ityp,1:imon)=MIN(xd(2:ityp,1:imon),dzcut)
  END SUBROUTINE vegin

  !
  !------------------------------------------------------------
  !
  SUBROUTINE re_assign_sib_soil_prop(iMax            , & ! IN
       jMax            , & ! IN
       npatches        , & ! IN
       imask           , & ! INOUT
       veg_type          ) ! IN
    INTEGER, INTENT(IN   )  :: iMax
    INTEGER, INTENT(IN   )  :: jMax
    INTEGER, INTENT(IN   )  :: npatches

    INTEGER(KIND=i8), INTENT(INOUT)  :: imask    (imax,jmax         )
    REAL(KIND=r8)   , INTENT(IN   )  :: veg_type (imax,jmax,npatches)
    REAL(KIND=r8)    :: GSWP_soil_input_data(10,12  )
    INTEGER :: nnn
    INTEGER :: i
    INTEGER :: j
    !
    !-------------------------------Soil data from GSWP-2 -------------------------------------
    !
    DATA GSWP_soil_input_data/  &
                                !1     2    3     4        5       6      7      8      9     10
                                !SAND(%) SILT CLAY QUARTZ  Wfc    Wwilt  Wsat    b    PHIsat  Ksat
         92.0_r8, 5.0_r8, 3.0_r8,0.92_r8,0.132_r8,0.033_r8,0.373_r8, 3.30_r8,-0.05_r8,2.45E-05_r8,&!1  Sand
         82.0_r8,12.0_r8, 6.0_r8,0.82_r8,0.156_r8,0.051_r8,0.386_r8, 3.80_r8,-0.07_r8,1.75E-05_r8,&!2  Loamy Sand
         58.0_r8,32.0_r8,10.0_r8,0.60_r8,0.196_r8,0.086_r8,0.419_r8, 4.34_r8,-0.16_r8,8.35E-06_r8,&!3  Sandy Loam
         10.0_r8,85.0_r8, 5.0_r8,0.25_r8,0.361_r8,0.045_r8,0.471_r8, 3.63_r8,-0.84_r8,1.10E-06_r8,&!4  Silt Loam
         17.0_r8,70.0_r8,13.0_r8,0.40_r8,0.270_r8,0.169_r8,0.476_r8, 5.25_r8,-0.65_r8,2.36E-06_r8,&!5  Loam
         58.0_r8,15.0_r8,27.0_r8,0.60_r8,0.253_r8,0.156_r8,0.412_r8, 7.32_r8,-0.12_r8,6.31E-06_r8,&!6  Sandy Clay Loam
         32.0_r8,34.0_r8,34.0_r8,0.10_r8,0.301_r8,0.211_r8,0.447_r8, 8.34_r8,-0.28_r8,2.72E-06_r8,&!7  Silty Clay Loam
         10.0_r8,56.0_r8,34.0_r8,0.35_r8,0.334_r8,0.249_r8,0.478_r8, 8.41_r8,-0.63_r8,1.44E-06_r8,&!8  Clay Loam
         52.0_r8, 6.0_r8,42.0_r8,0.52_r8,0.288_r8,0.199_r8,0.415_r8, 9.70_r8,-0.12_r8,4.25E-06_r8,&!9  Sandy Clay
         6.0_r8,47.0_r8,47.0_r8,0.10_r8,0.363_r8,0.286_r8,0.478_r8,10.78_r8,-0.58_r8,1.02E-06_r8,&!10 Silty Clay
         22.0_r8,20.0_r8,58.0_r8,0.25_r8,0.353_r8,0.276_r8,0.450_r8,12.93_r8,-0.27_r8,1.33E-06_r8,&!11 Clay
         43.0_r8,39.0_r8,18.0_r8,0.10_r8,0.250_r8,0.148_r8,0.437_r8, 5.96_r8,-0.24_r8,4.66E-06_r8 /!12 Silt
    !
    !-srf: avoid this for now, only use it when all arrays above are used like:
    ! bee(int(soil_type(lon,lat))) and not the usual way: bee(isurf(lon,lat))),
    ! where isurf is the vegetation index
    !
    !GO TO 332
    DO nnn = 1,12
       !
       !   sslfc(nnn)  = GSWP_soil_input_data(5,nnn)        !not in use
       !   sswlts(nnn) = GSWP_soil_input_data(6,nnn)        !not in use
       !   sswlts(nnn) = max(0.06_r8,GSWP_soil_input_data(6,n) !not in use nn)
       !
       ! print*,nnn,'poros bee phsat satco'
       ! print*,poros(nnn) , GSWP_soil_input_data(7,nnn)
       ! print*,bee(nnn)   ,GSWP_soil_input_data(8,nnn)
       ! print*,phsat(nnn) ,GSWP_soil_input_data(9,nnn)
       ! print*,satco(nnn) ,GSWP_soil_input_data(10,nnn)

       ! poros(nnn) = GSWP_soil_input_data(7,nnn)
       ! bee  (nnn) = GSWP_soil_input_data(8,nnn)
       ! phsat(nnn) = GSWP_soil_input_data(9,nnn)
       ! satco(nnn) = GSWP_soil_input_data(10,nnn)
    END DO

    !332 continue
    !
    !- for now, set isurf(:,:) as the veg data of the predominant biome:
    !
    DO j=1,jMax
       DO i= 1,iMax
          !imask(i,j) = 0 => ocean  / imask(i,j) = 13 => ice
          IF (imask(i,j) > 0_i8 .and. imask(i,j) < 13_i8) THEN
             !print*,'1',i,j,int(veg_type(i,j,2)),imask(i,j)
             imask(i,j) = int(veg_type(i,j,2))
          END IF
       END DO
    END DO
    !
    !stop 44433
    !srf- original SSIB from MCGA requires 13 soil classes, while USDA/GSWP2 has only 12
    !srf- the soil class 13 is not changed here (see vegin.f90)
    !  bee(13) = 4.8_r8
    !  phsat(13) = -0.167_r8
    !  satco(13) = 0.762e-4_r8
    !  poros(13) = 0.4352_r8
    !  zdepth(13,1) = 1.0_r8
    !  zdepth(13,2) = 1.0_r8
    !  zdepth(13,3) = 1.0_r8
    !
    RETURN
  END SUBROUTINE re_assign_sib_soil_prop

  ! wheat  :determine wheat phenology for latitude and julian day?.

  SUBROUTINE wheat (latitu,itype ,nmax  ,colrad ,month ,xday   ,yrl   , &
       idatec,monl  ,nsx    )
    !==========================================================================
    !==========================================================================
    !  ityp.......Numero das classes de solo vegetacao 13
    !  imon.......Number max of month at year (12)
    !  icg........Parameters of vagetation (icg=1 top e icg=2 bottom)
    !  iwv........Comprimento de onda iwv=1=visivel, iwv=2=infravermelho
    !             proximo, iwv=3 infravermelho termal
    !  nmax
    !  itype......Classe de textura do solo ou classe de vegetacao
    !  jmax.......Number of grid points on a gaussian longitude circle
    !  colrad.....colatitude
    !  month......Number of month at year (1-12)
    !  xday.......is julian day - 1 with fraction of day
    !  pie........Constante Pi=3.1415926e0
    !  yrl........length of year in days
    !  idatec.....idatec(1)=current hour of
    !            idatec(2)=current day of month.
    !            idatec(3)=current month of year.
    !            idatec(4)=current year.
    !  monl.......length of each month in days
    !  ystpar.....Coefficints related to par influence on
    !             stomatal resistance
    !  yopt.......Temperatura ideal de funcionamento estomatico
    !  yll........Temperatura minima de funcionamento estomatico
    !  yu.........Temperatura maxima de funcionamento estomatico
    !  yefac......Parametro de deficit de pressao de vapor d'agua
    !  yh1........Coeficiente para o efeito da agua no solo
    !  yh2........Potencial de agua no solo para ponto de Wilting
    !  rstpar.....Coefficints related to par influence on
    !             stomatal resistance
    !  chil.......Leaf orientation parameter
    !  topt.......Temperatura ideal de funcionamento estomatico
    !  tll........Temperatura minima de funcionamento estomatico
    !  tu.........Temperatura maxima de funcionamento estomatico
    !  defac......Parametro de deficit de pressao de vapor d'agua
    !  ph1........Coeficiente para o efeito da agua no solo
    !  ph2........Potencial de agua no solo para ponto de Wilting
    !  green......Fraction of grenn leaves
    !  xcover(iv).Fracao de cobertura de vegetacao iv=1 Top
    !  xcover(iv).Fracao de cobertura de vegetacao iv=2 Bottom
    !  nsx........phenology dates to fall within one year period
    !==========================================================================
    INTEGER , PARAMETER :: itveg = 13 ! Number of Vegetation Types
    INTEGER , PARAMETER :: isoil = 13 ! Number of Vegetation Types
    INTEGER , PARAMETER :: imon = 12 ! Number of Months with Defined Vegetation Types
    INTEGER , PARAMETER :: icg  = 2  ! Number of Vegetation Parameters
    INTEGER , PARAMETER :: iwv  = 3  ! Number of Radiation Wavelengths
    INTEGER , PARAMETER :: idp  = 3  ! Number of Soil Layer Parameters
    INTEGER , PARAMETER :: ibd  = 2  ! Number of Vegetation Stage

    INTEGER, INTENT(in ) :: nmax
    INTEGER, INTENT(in ) :: latitu
    INTEGER, INTENT(in ) :: itype (nmax)
    REAL(KIND=r8),    INTENT(in ) :: colrad(nmax)
    INTEGER, INTENT(in ) :: month (nmax)
    REAL(KIND=r8),    INTENT(in ) :: xday
    REAL(KIND=r8),    INTENT(in ) :: yrl
    INTEGER, INTENT(in ) :: idatec(4)
    INTEGER, INTENT(in ) :: monl  (12)
    INTEGER, INTENT(inout) :: nsx(nmax)
    REAL(KIND=r8)    :: rday
    REAL(KIND=r8)    :: thrsh
    REAL(KIND=r8)    :: phi(nmax)
    REAL(KIND=r8)    :: flip
    REAL(KIND=r8)    :: rootgc (nmax)
    REAL(KIND=r8)    :: chilw (nmax)
    REAL(KIND=r8)    :: tlai(nmax)
    REAL(KIND=r8)    :: xcover2(nmax)
    REAL(KIND=r8)    :: grlf (nmax)
    REAL(KIND=r8)    :: diff1 (nmax)
    REAL(KIND=r8)    :: diff2 (nmax)
    REAL(KIND=r8)    :: perc
    REAL(KIND=r8)    :: x1
    REAL(KIND=r8)    :: xdif1
    REAL(KIND=r8)    :: xdif2
    INTEGER :: i
    INTEGER :: kold
    INTEGER :: i1
    INTEGER :: ns
    INTEGER :: mind (nmax)
    INTEGER :: index (nmax)
    INTEGER :: icond
    INTEGER :: kk
    INTEGER :: mnl
    REAL(KIND=r8) :: pie=3.1415926e0_r8
    REAL(KIND=r8)    :: phenst(nmax,9)
    LOGICAL    :: test(nmax)
    INTEGER, PARAMETER :: iimon=12

    REAL(KIND=r8), PARAMETER :: wlai(9)=(/1.0_r8, 2.0_r8, 6.0_r8, 4.0_r8, 3.0_r8, 1.0_r8, 0.01_r8, 0.01_r8, 1.0_r8/)

    REAL(KIND=r8), PARAMETER :: xgren(iimon+1)=(/0.55_r8,0.68_r8,0.8_r8,0.9_r8,0.9_r8,0.9_r8,0.9_r8,0.81_r8,0.64_r8,&
         0.53_r8,0.49_r8,0.48_r8,0.55_r8/)

    REAL(KIND=r8), PARAMETER :: vlt(iimon+1)=(/1.0_r8,6.0_r8,6.0_r8,6.0_r8,6.0_r8,6.0_r8,6.0_r8,6.0_r8,6.0_r8,3.78_r8,&
         1.63_r8,1.0_r8,1.0_r8/)

    !    INTEGER, SAVE :: kmon(iimon+1)

    REAL(KIND=r8) :: xgreen(nmax,iimon+1)
    INTEGER :: kmon(imon+1)

    INTEGER, PARAMETER :: ihead = 3
    INTEGER, PARAMETER :: iwheat=12
    REAL(KIND=r8),    PARAMETER :: syr   =365.25e0_r8
    REAL(KIND=r8),    PARAMETER :: vcv   =0.569_r8
    !
    !     vlt and xgren are assumed to be correct at the beginning of the
    !     month
    !
    nsx = 0
    index= 0
    phenst=0.0_r8
    !
    !     xday is julian day - 1 with fraction of day
    !
    rday=xday
    !
    !     for standard length years, determine the offset for the year
    !     within the leap year period
    !
    thrsh=-MOD(idatec(4)+3,4)*0.25e0_r8
    test=.TRUE.
    DO i = 1, nmax
       !pi === 180
       !y  === x
       !
       ! X = (180 * Y)/pi
       !
       phi(i) = 90.0_r8-180.0e0_r8/pie * colrad(i)
       !
       !     constrain latitude range
       !
       !fixa o valor -55 ou +55 se o valor absoluto da latitude for maior que 55
       IF (ABS(phi(i)) > 55.0_r8) phi(i)=SIGN(55.0_r8,phi(i))
       !fixa o valor -20 ou +20 se a valor absoluto da latitude for menor que 20
       IF (ABS(phi(i)) < 20.0_r8) phi(i)=SIGN(20.0_r8,phi(i))

    ENDDO
    DO i1 = 1, iimon+1
       DO i = 1, nmax
          xgreen(i,i1)=xgren(i1)
       END DO
    END DO
    !
    !     search for any wheat vegetation points at this latitude
    !     if found, set sib parameters for latitude and time of year
    !
    kold=0
    DO i1 = 1, iimon
       kmon(i1)=kold
       !
       !     add extra day for leap years if using standard length year
       !
       IF (MOD(idatec(4),4) == 0 .AND. i1 == 2)kmon(i1)=kmon(i1)+1
       kold=kold+monl(i1)
    END DO
    DO i = 1, nmax
       IF (itype(i) /= iwheat) CYCLE
       flip =   0.0_r8
       IF (phi(i)< 0.0e0_r8) flip = yrl/2.0_r8
       !
       !     determine julian day - 1 for each wheat phenology for this
       !     latitude.  scale by length of year and adjust for south. hem.
       !
       phenst(i,2) = (4.50_r8 * ABS(phi(i)) - 65.0_r8) * (yrl/syr) + flip
       phenst(i,3) = (4.74_r8 * ABS(phi(i)) - 47.2_r8) * (yrl/syr) + flip
       phenst(i,4) = (4.86_r8 * ABS(phi(i)) - 31.8_r8) * (yrl/syr) + flip
       phenst(i,5) = (4.55_r8 * ABS(phi(i)) -  2.0_r8) * (yrl/syr) + flip
       phenst(i,6) = (4.35_r8 * ABS(phi(i)) + 10.5_r8) * (yrl/syr) + flip


       phenst(i,7) = phenst(i,6) + 3.0_r8 * (yrl/syr)
       phenst(i,1) = phenst(i,2) - ABS(5.21_r8 * ABS(phi(i)) - 0.3_r8)*(yrl/syr)
       phenst(i,9) = phenst(i,1)
       phenst(i,8) = phenst(i,9) - 5.0_r8*(yrl/syr)
    END DO
    DO ns = 1, 9
       DO i = 1, nmax
          IF (itype(i) /= iwheat) CYCLE
          !
          !     constrain phenology dates to fall within one year period
          !
          IF (phenst(i,ns) < 0.0e0_r8) phenst(i,ns) = phenst(i,ns) + yrl
          IF (phenst(i,ns) > yrl)      phenst(i,ns) = phenst(i,ns) - yrl
       END DO
    END DO

    DO i1 = 1, 12
       DO i = 1, nmax
          IF (itype(i) /= iwheat) CYCLE
          !
          !     find month of the head phenology stage for this latitude
          !
          IF (phenst(i,ihead) <= kmon(i1+1)) THEN
             mind(i) = i1
             IF (i1 <= 4) THEN
                xgreen(i,i1+1) = 0.9_r8
                xgreen(i,i1+2) = 0.9_r8
             END IF
          END IF
       END DO
    END DO

    DO ns = 1,8
       DO i = 1, nmax
          IF (itype(i) /= iwheat) CYCLE
          rootgc(i) = 1.0_r8
          chilw(i)  =-0.02_r8
          tlai(i)   = 0.5_r8
          grlf(i)   = 0.6_r8
          xcover2(i)=xcover(iwheat,month(i),1)
          !
          !     find growth stage given latitude and day
          !
          IF(test(i))THEN
             diff1(i) =  phenst(i,ns+1)- phenst(i,ns)
             diff2(i) = rday- phenst(i,ns)
             IF ( phenst(i,ns) >=  phenst(i,ns+1)) THEN
                IF ((rday <  phenst(i,ns)) .OR. (rday >  phenst(i,ns+1))) THEN
                   !
                   !     phenology stages overlap the end of year?
                   !
                   icond = 0
                   IF (rday >=  phenst(i,ns)   .AND. rday <= yrl  ) icond = 1
                   IF (rday >= thrsh .AND. rday <=  phenst(i,ns+1)) icond = 2
                   IF (icond /= 2) THEN
                      diff1(i) = yrl    -  phenst(i,ns) +  phenst(i,ns+1)
                      diff2(i) = rday   -  phenst(i,ns)
                   ELSE
                      diff1(i) = yrl   -  phenst(i,ns) + phenst(i,ns+1)
                      diff2(i) = yrl   -  phenst(i,ns) + rday
                   END IF
                END IF
                IF (icond /= 0) THEN
                   !
                   !     date found in phenology stage
                   !
                   perc =  diff2(i)/diff1(i)
                   !
                   !     kk is current month number
                   !
                   kk=idatec(2)
                   mnl=monl(kk)
                   IF (MOD(idatec(4),4) == 0 .AND. kk == 2)mnl=mnl+1
                   IF (rday > phenst(i,ihead)) THEN
                      IF (kk /= mind(i)) THEN
                         x1 = vlt(kk)
                         xdif1 = mnl
                         xdif2 = rday - kmon(kk)
                      ELSE
                         x1    = wlai(ihead)
                         xdif1 = kmon(kk+1) - phenst(i,ihead)
                         xdif2 = rday - phenst(i,ihead)
                      END IF
                      tlai(i) = x1 - (x1-vlt(kk+1)) / xdif1 * xdif2
                   ELSE
                      tlai(i) =  perc*(wlai(ns+1)-wlai(ns)) + wlai(ns)
                   END IF
                   IF (rday > phenst(i,ihead+1)) THEN
                      xcover2(i)=vcv + (0.9_r8 - vcv) * (yrl - rday)/(yrl - phenst(i,ihead+1))
                   ELSE
                      xcover2(i)=0.90_r8*(1.0_r8 - EXP(-tlai(i)))
                   END IF
                   grlf(i)   = xgreen(i,kk)-(xgreen(i,kk)-xgreen(i,kk+1))/mnl*(rday-kmon(kk))
                   rootgc(i) = 2910.0_r8 * (0.5_r8 + 0.5_r8 * tlai(i)/ wlai(ihead) * grlf(i))
                   IF (ns /= 1 .AND. ns /= 2) chilw(i)=-0.2_r8
                   test(i)=.FALSE.
                   index(i)=ns
                END IF
             END IF
          END IF
       END DO
    END DO

    DO i = 1, nmax
       IF (itype(i) /= iwheat) CYCLE
       nsx(i) = index(i)
       IF (nsx(i) == 9) nsx(i) = 1
       IF (nsx(i) >  6) nsx(i) = 6
       vcover_gbl (i,latitu,1) =   xcover2(i) !xcover(itype(i),month(i),1)
       zlt_gbl    (i,latitu,1) =   tlai(i)    !zlt   (itype(i),month(i),1)
       green_gbl  (i,latitu,1) =   grlf(i)    !green (itype(i),month(i),1)
       chil_gbl   (i,latitu,1) =   chilw(i)   !chil  (itype(i),1)
       topt_gbl   (i,latitu,1) =   yopt (2)   !topt  (itype(i),1)
       tll_gbl    (i,latitu,1) =   yll  (2)   !tll   (itype(i),1)
       tu_gbl     (i,latitu,1) =   yu(2)   !tu    (itype(i),1)
       defac_gbl  (i,latitu,1) =   yefac(2)   !defac (itype(i),1)
       ph1_gbl    (i,latitu,1) =   yh1  (2)   !ph1   (itype(i),1)
       ph2_gbl    (i,latitu,1) =   yh2  (2)   !ph2   (itype(i),1)
       rstpar_gbl (i,latitu,1,1)=  ystpar(2,1)!rstpar(itype(i),1,1)
       rstpar_gbl (i,latitu,1,2)=  ystpar(2,2)!rstpar(itype(i),1,2)
       rstpar_gbl (i,latitu,1,3)=  ystpar(2,3)!rstpar(itype(i),1,3)

    END DO
    RETURN
  END SUBROUTINE wheat



  ! sibwet :transform mintz-serafini and national meteoroLOGICAL center fields
  !         of soil moisture into sib compatible fields of soil moisture.




  SUBROUTINE sibwet &
       (ibmax,jbmax,sinp,sinmax,imask,wsib,ssib,mxiter,ibMaxPerJB)
    !
    !
    !     piers sellers : 29 april 1987
    !
    !
    !   input  :   sinp     = mintz-serafini or national meteoroLOGICAL
    !                         center soil moisture (mm)
    !              sinmax   = maximum value of sinp (mm)
    !              wsinp    = m-s or nmc fractional wetness
    !              ms       = 1, mintz-serafini
    !              nmc      = 1, national meteoroLOGICAL center
    !              bee      = sib : soil moisture potential factor
    !              phsat    = sib : soil potential at saturation (m)
    !              zdepth(3)= sib : depth of 3 soil layers (m)
    !              poros    = Porosidade do solo (m"3/m"3)
    !
    !   output :   wsibt    = sib : fractional wetness
    !              ssibt    = sib : soil moisture content (m)
    !              psit     = sib : soil moisture potential (m)
    !              factor   = sib : extraction factor
    !
    INTEGER, INTENT(in   ) :: ibmax
    INTEGER, INTENT(in   ) :: jbmax
    INTEGER, INTENT(in   ) :: mxiter
    REAL(KIND=r8)   , INTENT(in   ) :: sinp(ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(in   ) :: sinmax
    !

    INTEGER(KIND=i8), INTENT(in     ) :: imask (ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(inout  ) :: wsib  (ibmax,jbmax)
    REAL(KIND=r8)   , INTENT(inout  ) :: ssib  (ibmax,jbmax)
    INTEGER         , INTENT(in     ) :: ibMaxPerJB(jbmax)

    REAL(KIND=r8) :: sm   (ityp,mxiter)
    REAL(KIND=r8) :: time (ityp,mxiter)
    REAL(KIND=r8) :: fact (ityp,mxiter)

    REAL(KIND=r8), PARAMETER :: xph1(13,2) = RESHAPE( &
         (/-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8,-120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8, &
         -200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8,-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8, &
         -120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8,-200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8/), &
         (/13,2/))
    REAL(KIND=r8), PARAMETER :: xph2(13,2) = RESHAPE( &
         (/-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8, &
         -400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8,-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8, &
         -230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8,-400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8/) , &
         (/13,2/))

    REAL(KIND=r8)    :: tzdep  (3)
    REAL(KIND=r8)    :: tzltm  (2)
    REAL(KIND=r8)    :: sibmax(ityp)
    REAL(KIND=r8)    :: tphsat
    REAL(KIND=r8)    :: tbee
    REAL(KIND=r8)    :: tporos
    INTEGER :: imm1
    INTEGER :: imm2
    INTEGER :: is
    INTEGER :: im
    INTEGER :: imm
    INTEGER :: ivegm
    REAL(KIND=r8)    :: cover
    REAL(KIND=r8)    :: tph1
    REAL(KIND=r8)    :: tph2
    REAL(KIND=r8)    :: sref
    REAL(KIND=r8)    :: smin
    REAL(KIND=r8)    :: dssib
    REAL(KIND=r8)    :: dw
    REAL(KIND=r8)    :: times
    REAL(KIND=r8)    :: soilmo
    REAL(KIND=r8)    :: w
    REAL(KIND=r8)    :: rsoilm
    INTEGER :: iter
    INTEGER :: latmax
    INTEGER :: lonmax
    INTEGER :: lat
    INTEGER :: lon
    REAL(KIND=r8)    :: tsinp
    REAL(KIND=r8)    :: etp
    REAL(KIND=r8)    :: facmod
    REAL(KIND=r8)    :: ssibt
    REAL(KIND=r8)    :: psit
    REAL(KIND=r8)    :: factor
    REAL(KIND=r8)    :: dt
    INTEGER :: itsoil
    INTEGER :: itfac

    sm  =0.0_r8
    time=0.0_r8
    fact=0.0_r8
    ssib=0.0_r8
    wsib=0.0_r8

    lonmax=ibmax
    latmax=jbmax

    DO is = 1,ityp
       !zdepth(3)= sib : depth of 3 soil layers (m)
       tzdep (1)= zdepth(is,1)
       tzdep (2)= zdepth(is,2)
       tzdep (3)= zdepth(is,3)
       tphsat   = phsat (is)
       tbee     = bee   (is)
       tporos   = poros (is)
       imm1=1
       imm2=1
       tzltm(1)=zlt_fixed(is,1,1)
       tzltm(2)=zlt_fixed(is,1,2)
       DO im=2,12
          IF(tzltm(1).LE.zlt_fixed(is,im,1) ) THEN
             imm1=im
             tzltm(1)=zlt_fixed(is,im,1)
          END IF
          IF(tzltm(2).LE.zlt_fixed(is,im,2) )THEN
             imm2=im
             tzltm(2)=zlt_fixed(is,im,2)
          END IF
       END DO
       imm=imm1
       ivegm=1
       IF(tzltm(1).LE.tzltm(2)) THEN
          imm=imm2
          ivegm=2
       END IF
       !
       !     xcover......Fracao de cobertura vegetal icg=1 topo
       !     xcover......Fracao de cobertura vegetal icg=2 base
       !
       cover=xcover_fixed(is,imm,ivegm)
       tph1=xph1         (is,ivegm)
       tph2=xph2         (is,ivegm)
       !
       !                                                     m^3
       ! sibmax(is) =(Z1 + Z2 + Z3) * poros = [m + m + m] * ----- = m = Os
       !                                                     m^3
       !
       sibmax(is) = ( tzdep(1) + tzdep(2) + tzdep(3) ) * tporos
       !
       IF(nfctrl(83).GE.1)WRITE(UNIT=nfprt,FMT=999)is,sibmax(is),tzdep(1), &
            tzdep(2),tzdep(3),tporos
       !
       !            bee      = soil moisture potential factor
       !            phsat    = soil potential at saturation   (m)
       !
       !                   --              --
       !                  | log ( - tphsat/1)|
       !  O  = Os * EXP * | -----------------|
       !                  |        b         |
       !                   --              --
       !
       sref = sibmax(is) * EXP( LOG(tphsat /(-1.0e0_r8)) /tbee)
       !                   --                          --
       !                  | log ( - tphsat/(-1.0e10) )   |
       !Omin = Os * EXP * | -----------------------------|
       !                  |              b               |
       !                   --                          --
       !
       smin    = sibmax(is) * EXP( LOG(tphsat /(-1.0e10_r8)) / tbee)
       !
       !             O - Omin
       !dssib  = ------------------
       !              mxiter
       !
       dssib   = (sref - smin) / REAL(mxiter,r8)
       !
       !              O - Omin
       ! dw    =  ------------------
       !             mxiter*Os
       !
       dw      = dssib / sibmax(is)
       !
       times   = 0.0e0_r8
       soilmo  = sref
       !
       !       O
       ! w = -----
       !       Os
       !
       w = soilmo / sibmax(is)
       !
       !                      --             --
       !                     |       0.0027    |
       !rsoilm  = 101840.0 * |1.0 - w          |
       !                     |                 |
       !                      --             --
       !
       rsoilm  = 101840.0_r8 * (1.0_r8 - w**0.0027_r8)
       DO iter = 1, mxiter
          CALL extrak( w   ,dw  ,tbee,tphsat, rsoilm, cover, &
               tph1,tph2,psit,factor )
          !
          !       dssib
          !dt = ----------
          !       factor
          !
          dt            = dssib  / factor
          !
          soilmo        = soilmo - dssib
          !
          !       O
          ! w = -----
          !       Os
          !
          w             = soilmo / sibmax(is)
          times         = times  + dt
          sm  (is,iter) = soilmo
          time(is,iter) = times
          fact(is,iter) = factor
       END DO

    END DO
    !
    !     input soil moisture map is now transformed to sib fields.
    !
    DO lat = 1, latmax
       DO lon = 1, ibMaxPerJB(lat)
          is=INT(imask(lon,lat),kind=i4)
          IF(is.NE.0)THEN
             tsinp = sinp(lon,lat)
             tsinp = MAX (sinmax/100.0e3_r8 , tsinp )
             tsinp = MIN (sinmax,tsinp)
             IF (tsinp .GT. 0.75e0_r8*sinmax ) etp = sinmax - tsinp
             facmod=MIN(1.0e0_r8,tsinp/(0.75e0_r8*sinmax) )
             IF (tsinp .LE. 0.75e0_r8*sinmax ) THEN
                etp = 0.75e0_r8*sinmax*LOG(0.75e0_r8*sinmax/tsinp ) + 0.25e0_r8*sinmax
             END IF
             etp = etp / 1000.0e0_r8
             DO iter = 1, mxiter
                itsoil=iter
                IF ( time(is,iter) - etp .GT. 0.0e0_r8  ) EXIT
             END DO
             DO iter=1,mxiter
                itfac=iter
                IF( fact(is,iter)-facmod-0.01e0_r8.LT.0.0e0_r8)EXIT
             END DO
             ssibt=MIN(sm(is,itsoil),sm(is,itfac))
             DO iter=1,mxiter
                IF(ssibt.GT.sm(is,iter))EXIT
             END DO
             ssib(lon,lat) = sm(is,iter)
             !
             !          O
             ! wsib = -----
             !         Os
             !
             wsib(lon,lat) = sm(is,iter) / sibmax(is)
          END IF
       END DO
    END DO
999 FORMAT(' IS,MAX,D1,D2,D3,POROS=',I2,1X,5E12.5)
  END SUBROUTINE sibwet



  SUBROUTINE sibwet_GLSM (ibMax          , & ! IN
       jbMax          , & ! IN
       imask          , & ! IN
       wsib           , & ! IN
       ssib           , & ! IN
       mxiter         , & ! OUT
       ibMaxPerJB     , & ! OUT
       soilm          , & ! in
       nzg         , & ! in
       wsib3d         , & ! OUT
       glsm_w)            ! IN

    !
    ! $Author: pkubota $
    ! $Date: 2008/04/09 12:42:57 $
    ! $Revision: 1.12 $
    !
    ! sibwet :transform mintz-serafini and national meteoroLOGICAL center fields
    !         of soil moisture into sib compatible fields of soil moisture.
    !
    !     piers sellers : 29 april 1987
    !
    INTEGER, INTENT(IN   )            :: ibMax
    INTEGER, INTENT(IN   )            :: jbMax
    INTEGER, INTENT(IN   )            :: mxiter
    REAL(KIND=r8)   , INTENT(OUT  )            :: soilm          (ibMax,jbMax)
    INTEGER(KIND=i8), INTENT(IN   )            :: imask          (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  )            :: wsib           (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  )            :: ssib           (ibMax,jbMax)
    INTEGER, INTENT(in   )            :: ibMaxPerJB     (:)
    INTEGER, INTENT(in   )            :: nzg
    REAL(KIND=r8)   , INTENT(OUT  )            :: wsib3d    (ibMax,jbMax,8       )
    REAL(KIND=r8)   , INTENT(IN   )            :: glsm_w    (ibMax,jbMax,nzg     )

    REAL(KIND=r8)               :: sm  (ityp,mxiter)
    REAL(KIND=r8)               :: time(ityp,mxiter)
    REAL(KIND=r8)               :: fact(ityp,mxiter)
    !
    !-srf
    !
    INTEGER, PARAMETER :: nzgmax=20
    REAL(KIND=r8)               :: glsm_w1d  (0:nzgmax)     ! dummy 1d initial soil  wetness
    REAL(KIND=r8)               :: glsm_tzdep(0:3)          ! sib soil levels
    REAL(KIND=r8)               :: glsm_w_sib(0:3)          ! SIB dummy 1d initial and interpolated soil  wetness
    !
    !-srf
    !
    REAL(KIND=r8)               :: tzdep (3)
    REAL(KIND=r8)               :: tzltm (2)
    REAL(KIND=r8)               :: sibmax(ityp)
    INTEGER            :: k
    REAL(KIND=r8)               :: fx
    INTEGER            :: lonmax
    INTEGER            :: latmax
    INTEGER            :: is
    REAL(KIND=r8)               :: tphsat
    REAL(KIND=r8)               :: tbee
    REAL(KIND=r8)               :: tporos
    INTEGER            :: imm1
    INTEGER            :: imm2
    INTEGER            :: im
    INTEGER            :: imm
    INTEGER            :: ivegm
    REAL(KIND=r8)               :: cover
    REAL(KIND=r8)               :: tph1
    REAL(KIND=r8)               :: tph2
    REAL(KIND=r8)               :: sref
    REAL(KIND=r8)               :: smin
    REAL(KIND=r8)               :: dssib
    REAL(KIND=r8)               :: dw
    REAL(KIND=r8)               :: times
    REAL(KIND=r8)               :: soilmo
    REAL(KIND=r8)               :: w
    REAL(KIND=r8)               :: rsoilm
    INTEGER            :: iter
    REAL(KIND=r8)               :: psit
    REAL(KIND=r8)               :: factor
    REAL(KIND=r8)               :: dt
    INTEGER            :: lat
    INTEGER            :: lon

    !
    !              wsinp    = m-s or nmc fractional wetness
    !              ms       = 1, mintz-serafini
    !              nmc      = 1, national meteoroLOGICAL center
    !              bee      = sib : soil moisture potential factor
    !              phsat     = sib : soil potential at saturation (m)
    !              zdepth(3)= sib : depth of 3 soil layers (m)
    !              poros    = sib : soil porosity
    !              ph1      = sib : leaf potential, stress onset (m)
    !              ph2      = sib : leaf potential, no e-t (m)
    !
    !   output :   wsibt    = sib : fractional wetness
    !              ssibt    = sib : soil moisture content (m)
    !              psit     = sib : soil moisture potential (m)
    !              factor   = sib : extraction factor
    !
    REAL(KIND=r8), PARAMETER :: xph1(13,2) = RESHAPE( &
         (/-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8,-120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8, &
         -200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8,-100.0_r8,-190.0_r8,-200.0_r8,-200.0_r8,-200.0_r8, &
         -120.0_r8,-120.0_r8,-120.0_r8,-200.0_r8,-200.0_r8, -10.0_r8,-190.0_r8, -10.0_r8/), &
         (/13,2/))
    REAL(KIND=r8), PARAMETER :: xph2(13,2) = RESHAPE( &
         (/-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8, &
         -400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8,-500.0_r8,-250.0_r8,-250.0_r8,-250.0_r8,-250.0_r8, &
         -230.0_r8,-230.0_r8,-280.0_r8,-400.0_r8,-400.0_r8,-100.0_r8,-250.0_r8,-100.0_r8/) , &
         (/13,2/))

    !-srf
    !hmjb
    !    REAL, PARAMETER :: glsm_slz(0:nzgmax) = (/  0., 0.1, 0.25, 0.5, 1., 2., 3.,& !7  values
    !         0., 0.,  0.,   0.,  0., 0., 0., 0., 0., 0.,& !10 values
    !         0., 0.,  0.,   0.                         /) !4  values
    !versao para NZG=8 => 9 niveis no MCGA
    REAL(KIND=r8), PARAMETER :: glsm_slz(0:nzgmax) = (/  0.0_r8, 0.05_r8, 0.13_r8, 0.25_r8, 0.5_r8, 1.0_r8, 1.75_r8,& !9  values
         2.5_r8, 4.5_r8,  0.0_r8,   0.0_r8,  0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8,& !10 values
         0.0_r8, 0.0_r8,  0.0_r8,   0.0_r8                         /) !4  values
    !-srf
    sm =0.0_r8
    time  =0.0_r8
    fact  =0.0_r8
    ssib  =0.0_r8
    wsib  =0.0_r8

    lonmax=ibMax
    latmax=jbMax

    DO is = 1,ityp
       tzdep (1)= zdepth(is,1)
       tzdep (2)= zdepth(is,2)
       tzdep (3)= zdepth(is,3)
       tphsat   = phsat (is)
       tbee     = bee   (is)
       tporos   = poros (is)
       imm1=1
       imm2=1
       tzltm(1)=zlt_fixed(is,1,1)
       tzltm(2)=zlt_fixed(is,1,2)
       DO im=2,12
          IF (tzltm(1).le.zlt_fixed(is,im,1) ) THEN
             imm1=im
             tzltm(1)=zlt_fixed(is,im,1)
          END IF

          IF (tzltm(2).le.zlt_fixed(is,im,2) ) THEN
             imm2=im
             tzltm(2)=zlt_fixed(is,im,2)
          END IF
       END DO

       imm=imm1
       ivegm=1

       IF (tzltm(1).le.tzltm(2)) THEN
          imm=imm2
          ivegm=2
       END IF
       cover=xcover_fixed(is,imm,ivegm)
       tph1=xph1(is,ivegm)
       tph2=xph2(is,ivegm)
       !
       !srf- max water content
       !
       sibmax(is) = ( tzdep(1) + tzdep(2) + tzdep(3) ) * tporos
       IF (nfctrl(83).ge.1) WRITE(UNIT=nfprt,FMT=999)is,sibmax(is),tzdep(1), &
            tzdep(2),tzdep(3),tporos
       sref = sibmax(is) * exp( log(tphsat /(-1.0e0_r8)) /tbee )
       smin = sibmax(is) * exp( log(tphsat /(-1.0e10_r8)) /tbee )
       dssib= (sref - smin) / REAL(mxiter,r8)
       dw   = dssib / sibmax(is)
       times  = 0.0e0_r8
       soilmo = sref
       w      = soilmo / sibmax(is)
       rsoilm = 101840.0_r8 * (1.0_r8 - w**0.0027_r8)

       DO iter = 1, mxiter
          CALL extrak      ( w     , &  ! IN
               dw    , &  ! IN
               tbee  , &  ! IN
               tphsat, &  ! IN
               rsoilm, &  ! IN
               cover , &  ! IN
               tph1  , &  ! IN
               tph2  , &  ! IN
               psit  , &  ! OUT
               factor  )  ! OUT
          dt            = dssib / factor
          soilmo        = soilmo - dssib
          w             = soilmo / sibmax(is)
          times         = times + dt
          sm  (is,iter) = soilmo
          time(is,iter) = times
          fact(is,iter) = factor
       END DO
    END DO
    !
    !     input soil moisture map is now transformed to sib fields.
    !
    DO lat = 1, latmax
       DO lon = 1, ibMaxPerJB(lat)

          wsib3d(lon,lat,:) = 0.e0_r8

          is=INT(imask(lon,lat),kind=i4)
          IF (is.ne.0) THEN

             tzdep (1)= zdepth(is,1)
             tzdep (2)= zdepth(is,2)
             tzdep (3)= zdepth(is,3)
             tphsat   = phsat (is)
             tbee     = bee   (is)
             tporos   = poros (is)
             !
             !-sib soil levels
             !
             glsm_tzdep(0) = 0.e0_r8
             glsm_w_sib(0) = 0.e0_r8

             DO k=1,3
                glsm_tzdep (k) = zdepth(is,k) + glsm_tzdep (k-1)
                glsm_w_sib (k) = 0.e0_r8
             END DO
             !
             !- copy 3d soil moisture array to 1d column array
             !
             DO k=1,nzg
                glsm_w1d(k)=glsm_w(lon,lat,k)
             END DO
             !
             !- performs vertical interpolation from soil moisture
             !  levels to sib levels
             !
             CALL vert_interp(4               , &  ! IN
                  nzg+1           , &  ! IN
                  glsm_tzdep(0:3) , &  ! IN
                  glsm_slz(0:nzg) , &  ! IN
                  glsm_w1d(0:nzg) , &  ! IN
                  glsm_w_sib(0:3)  )   ! OUT


             !endif
             !
             !- stores 1d sib soil moisture at 3d array
             !
             DO k=1,3
                wsib3d(lon,lat,k) = glsm_w_sib(k)
             END DO
             !
             !------------------------- remove this later--------------------------------X
             !- for now fix zero soil moisture inside the land
             !- latter fix this at soil moisture original data
             !
             !IF (imask(lon,lat) > 0 ) THEN
             !   ssm=0.
             !   DO k=1,3
             !      ssm=ssm+wsib3d(lon,lat,k)
             !   END DO
             !
             !   IF (ssm < 0.15) THEN
             !      !
             !      !print*,'SM null inside land portion', imask(lon,lat)
             !      !print*,'1',lon,lat,wsib3d(lon,lat,:)
             !      !
             !      ssm1d(:) = 0.
             !      ncount = 0
             !      DO i=max(1,lon-4),min(lonmax,lon+4)
             !         DO j=max(1,lat-4),min(latmax,lat+4)
             !    IF (imask(i,j) > 0) THEN !only points inside the land
             !       ssm=0.
             !       DO k=1,3
             !  ssm=ssm+wsib3d(i,j,k)
             !       END DO
             !
             !       IF (ssm > 0.15) THEN
             !  ncount=ncount  + 1
             !  ssm1d(:) = ssm1d(:) + wsib3d(i,j,:)
             !       END IF
             !    END IF
             !         END DO
             !      END DO
             !
             !     IF (ncount > 1) THEN
             !         wsib3d(lon,lat,:)=ssm1d(:)/float(ncount)
             !      ELSE
             !         wsib3d(lon,lat,:)=0.5
             !      END IF
             !      !
             !      !print*,'2',lon,lat,wsib3d(lon,lat,:)
             !      !
             !   END IF
             !END IF
             !
             !-----------------------------------------------------------------------------X
             !
             ssib(lon,lat) = 0.0_r8
             wsib(lon,lat) = 0.0_r8

             DO k=1,3

                fx            = ( glsm_tzdep(k)-glsm_tzdep(k-1) ) / glsm_tzdep(3)
                wsib(lon,lat) = wsib(lon,lat) + glsm_w_sib(k) * fx
                ssib(lon,lat) = ssib(lon,lat) + glsm_w_sib(k) * fx * tporos

             END DO
             !
             ! total water in mm
             !
             soilm(lon,lat) = ( tzdep(1)*wsib3d(lon,lat,1) + &
                  tzdep(2)*wsib3d(lon,lat,2) + &
                  tzdep(3)*wsib3d(lon,lat,3) ) * tporos
             !
          END IF
       END DO
    END DO
999 FORMAT(' IS,MAX,D1,D2,D3,POROS=',I2,1X,5E12.5)
  END SUBROUTINE sibwet_GLSM

  SUBROUTINE extrak( w, dw, tbee, tphsat, rsoilm, cover, tph1, tph2, &
       psit, factor )
    REAL(KIND=r8), INTENT(in   ) :: w
    REAL(KIND=r8), INTENT(in   ) :: dw
    REAL(KIND=r8), INTENT(in   ) :: tbee
    REAL(KIND=r8), INTENT(in   ) :: tphsat
    REAL(KIND=r8), INTENT(in   ) :: rsoilm
    REAL(KIND=r8), INTENT(in   ) :: cover
    REAL(KIND=r8), INTENT(in   ) :: tph1
    REAL(KIND=r8), INTENT(in   ) :: tph2
    REAL(KIND=r8), INTENT(inout  ) :: psit
    REAL(KIND=r8), INTENT(inout  ) :: factor
    REAL(KIND=r8) :: rsoil
    REAL(KIND=r8) :: argg
    REAL(KIND=r8) :: hr
    REAL(KIND=r8) :: rplant
    !                --     -- (-b)
    !               |      dw |                  0
    ! psit = PHYs * | w - --- |      where w = -----
    !               |      2  |                  0s
    !                --     --
    psit   = tphsat * ( w-dw/2.0e0_r8 ) ** (-tbee)
    !
    !                      --                        --
    !                     |       --     -- (0.0027)   |
    !                     |      |      dw |           |
    !rsoil   = 101840.0 * |1.0 - | w - --- |           |
    !                     |      |      2  |           |
    !                     |       --     --            |
    !                      --                        --
    !
    rsoil  = 101840.0_r8 * (1.0_r8-( w-dw/2.0_r8) ** 0.0027_r8)
    !
    !                9.81       1
    !argg = psit * -------- * -------
    !               461.50     310.0
    !
    argg   = MAX ( -10.0e0_r8 , ((psit * 9.81e0_r8 / 461.5e0_r8) / 310.e0_r8))
    !
    !            --                       --
    !           |         9.81       1      |
    !hr   = EXP |psit * -------- * -------  |
    !           |        461.50     310.0   |
    !            --                       --
    !
    hr     = EXP ( argg )
    !
    !         rsoilm
    ! rsoil =--------- * hr
    !         rsoil
    !
    rsoil  = rsoilm /rsoil * hr
    !
    !          ( psit - tph2 - 50.0)
    !rplant = -------------------------
    !             ( tph1 - tph2 )
    !
    rplant = ( psit - tph2 -50.0_r8) / ( tph1 - tph2 )
    rplant = MAX ( 0.0e0_r8, MIN ( 1.0e0_r8, rplant ) )
    !                                                                     --                   --
    !                  --                 --                             |     --     -- (0.0027)|
    !                 |( psit - tph2 - 50)  |                            |    |      dw |        |
    !factor = cover * |---------------------| + (1 - cover) * 101840.0 * |1 - | w - --- |        |
    !                 |   ( tph1 - tph2 )   |                            |    |      2  |        |
    !                  --                 --                             |     --     --         |
    !                                                                     --                   --
    factor = cover * rplant + ( 1.0e0_r8 - cover ) * rsoil
    factor = MAX ( 1.e-6_r8, factor )
  END SUBROUTINE extrak
  !
  !------------------------------------------------------------
  !
  SUBROUTINE vert_interp(nsib      , & ! IN
       nzg       , & ! IN
       tzdep     , & ! IN
       glsm_slz  , & ! IN
       gl_sm     , & ! IN
       glsm_w_sib  ) ! OUT

    INTEGER, INTENT(IN   ) :: nsib
    INTEGER, INTENT(IN   ) :: nzg
    REAL(KIND=r8)   , INTENT(IN   ) :: tzdep     (nsib)
    REAL(KIND=r8)   , INTENT(IN   ) :: glsm_slz  (: )
    REAL(KIND=r8)   , INTENT(IN   ) :: gl_sm     (: )
    REAL(KIND=r8)   , INTENT(OUT  ) :: glsm_w_sib(nsib)

    REAL(KIND=r8)    :: zm        (nsib)
    REAL(KIND=r8)    :: wf        (nsib)
    REAL(KIND=r8)    :: zc        (nzg )
    REAL(KIND=r8)    :: wi        (nzg )
    REAL(KIND=r8)    :: dzlft
    INTEGER :: ZDM
    INTEGER :: k
    INTEGER :: kstart
    INTEGER :: L

    DO k=1,nzg
       zc(k)=glsm_slz(k)
    END DO

    DO k=1,nsib
       zm(k)=tzdep(k)
    END DO

    zdm=nsib
    KSTART=3
    !
    ! Transfere valores da grade de MAIOR resolucao (WI)
    !                     para a grade de MENOR resolucao (WF)
    !
    ! OS valores de WI devem estar definidos nos pontos de grade ZCS=zc/2
    ! OS valores de WF saem nos niveis ZMS = ZM/2
    !
    !
    !
    !    Dados da grade de maior resolucao
    !
    DO K=1,NZG
       WI(K) = gl_sm(k)
       !print*,'wi=',k,wi(k)
    END DO
    !
    !     Dado interpolado
    !
    wf(:)=0.0_r8
    !
    !     Valor de superficie:
    !
    WF(1)=WI(2)
    WF(2)=WI(2)
    !
    !
    DZLFT=0.0_r8
    L=2
    DO K=KSTART,ZDM
       !
       !    if(k==4) print*,'0',l,WF(K),WI(L),DZLFT
       !
       IF(DZLFT.NE.0.0_r8) THEN

          WF(K)=WF(K)+WI(L)*DZLFT
          !    if(k==4) print*,'1',l,WF(K),WI(L),DZLFT
          L=L+1

       END IF

70     CONTINUE

       IF(ZC(L).LE.ZM(K)) THEN

          WF(K)=WF(K)+WI(L)*(ZC(L)-ZC(L-1))

          !   if(k==4) print*,'2',l,WF(K),WI(L),ZC(L),zm(k)

          L=L+1
          DZLFT=0.0_r8
          IF (L>nzg) GO TO 1000
          GO TO 70
       ELSE

          WF(K)=WF(K)+WI(L)*(ZM(K)-ZC(L-1))
          DZLFT=ZC(L)-ZM(K)
       ENDIF
    ENDDO

1000 CONTINUE

    DO K=KSTART,ZDM
       !
       !   WF(K) =WF(K)/(ZM(K)-ZM(K-1))
       !         if(k==4)print*,zm(k),zc(nzg),ZM(K-1),WF(K)
       !
       IF (ZM(K) > ZC(nzg)) THEN
          WF(K) = WF(K)/(ZC(NZG)-ZM(K-1))
       ELSE
          WF(K)  = WF(K)/(ZM(K)-ZM(K-1))
       END IF
    END DO
    !
    !valores na grade do SIB
    !
    DO k=1,nsib
       glsm_w_sib(k)=WF(k)
       !print*,'SIB',k,glsm_w_sib(k)
    END DO
    !
    !check conservacao
    !srf - verifique se a integral de ambos calculos percorrem
    !srf - o mesmo intervalo
    !      print*,'        '
    !      sumf=0.0_r8
    !      DO K=2,ZDM
    !       sumf=sumf+wf(k)*(ZM(K)-ZM(K-1))
    !       print*,sumf,wf(k),zm(k),ZM(K)-ZM(K-1)
    !      ENDDO
    !      print*,'--------sumf-----',sumf
    !      sumi=0.0_r8
    !      DO K=2,nzg
    !       sumi=sumi+wi(k)*(glsm_slz(K)-glsm_slz(K-1))
    !       print*,k,sumi,wi(k),glsm_slz(K),(glsm_slz(K)-glsm_slz(K-1))
    !      ENDDO
    !      print*,'--------sumi-----',sumi, 100*(sumf-sumi)/sumi
    !
    RETURN
  END SUBROUTINE vert_interp


  !------------------------------------------------------------
  SUBROUTINE read_gl_sm_bc(iMax           , &!   IN
                           jMax           , &!   IN
                           jbMax          , &!   IN
                           nzg            , &!   IN
                           npatches       , &
                           npatches_actual, &
                           ibMaxPerJB     , &!   IN
                           record_type    , &! IN
                           fNameSoilType  , &! IN
                           fNameVegType   , &! IN
                           fNameSoilMoist , &
                           glsm_w         , &
                           soil_type      , &
                           veg_type       , &
                           frac_occ       , &
                           VegType         )! IN

    INTEGER, INTENT(IN   )            :: iMax
    INTEGER, INTENT(IN   )            :: jMax
    INTEGER, INTENT(IN   )            :: jbMax
    INTEGER, INTENT(IN   )            :: nzg
    INTEGER, INTENT(IN   )            :: npatches
    INTEGER, INTENT(IN   )            :: npatches_actual
    INTEGER, INTENT(IN   )            :: ibMaxPerJB(:)
    CHARACTER(LEN=*), INTENT(IN   )   :: record_type
    CHARACTER(LEN=*), INTENT(IN   )   :: fNameSoilType
    CHARACTER(LEN=*), INTENT(IN   )   :: fNameVegType
    CHARACTER(LEN=*), INTENT(IN   )   :: fNameSoilMoist
    REAL(KIND=r8) , INTENT(INOUT   )  :: glsm_w(:,:,:)
    REAL(KIND=r8) , INTENT(OUT   )  :: soil_type(:,:)
    REAL(KIND=r8) , INTENT(OUT   )  :: veg_type(:,:,:)      
    REAL(KIND=r8) , INTENT(OUT   )  :: frac_occ (:,:,:)     
    REAL(KIND=r8) , INTENT(INOUT   ):: VegType(iMax,jMax,npatches) ! SIB veg type
    REAL(KIND=r8) :: FracOcc(iMax,jMax,npatches) ! fractional area
    REAL(KIND=r8) :: glsm(iMax,jMax,nzg     )   ! initial soil wetness data
    REAL(KIND=r8) :: SoilType(iMax,jMax         )! FAO/USDA soil texture
    !
    ! Local
    !
    INTEGER            :: i
    INTEGER            :: j
    INTEGER            :: k
    INTEGER            :: ipatch
    INTEGER            :: ierr
    REAL(KIND=r8)      :: fractx
    !
    !------------------------- soil type initialization ------------
    !
    call MsgOne('**(read_gl_sm_bc)**','Opening GL soil file='//TRIM(fNameSoilType))

    FracOcc=0.0_r8
    glsm=0.0_r8
    SoilType=0.0_r8
    IF (record_type == 'seq') THEN      !sequential mode

       OPEN(UNIT=nfsoiltp,FILE=TRIM(fNameSoilType),FORM='unformatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilType), ierr
          STOP "**(ERROR)**"
       END IF
       READ(UNIT=nfsoiltp) ((SoilType(i,j),i=1,iMax),j=1,jMax)

       IF (reducedGrid) THEN
          CALL NearestIJtoIBJB(SoilType,soil_type)
       ELSE
          CALL IJtoIBJB( SoilType,soil_type)
       END IF

       CLOSE(UNIT=nfsoiltp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfsoiltp,FILE=TRIM(fNameSoilType),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilType), ierr
          STOP "**(ERROR)**"
       END IF

       CALL  vfirec(nfsoiltp,SoilType,imax*jmax,'LIN')

       IF (reducedGrid) THEN
          CALL NearestIJtoIBJB(SoilType,soil_type)
       ELSE
          CALL IJtoIBJB( SoilType,soil_type)
       END IF

       CLOSE(UNIT=nfsoiltp)

    END IF
    DO i=1,iMax
       DO j=1,jMax
          SoilType(i,j)=REAL(INT(SoilType(i,j)+0.1_r8),r8)
       END DO
    END DO
    DO j=1,jbMax
       DO i=1,ibMaxPerJB(j)
          soil_type(i,j)=REAL(INT(soil_type(i,j)+0.1_r8),r8)
       END DO
    END DO
    !
    !-------------------veg type and fractional area initialization ------------
    !
    call MsgOne('**(read_gl_sm_bc)**','Opening GL veg file='//TRIM(fNameVegType))

    IF (record_type == 'seq') THEN !sequential mode

       OPEN(UNIT=nfvegtp,FILE=TRIM(fNameVegType),FORM='unformatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameVegType), ierr
          STOP "**(ERROR)**"
       END IF

       DO ipatch=1,npatches_actual

          READ(UNIT=nfvegtp) ((VegType(i,j,ipatch),i=1,iMax),j=1,jMax) !veg dominante no patch
          READ(UNIT=nfvegtp) ((FracOcc(i,j,ipatch),i=1,iMax),j=1,jMax) !fracao ocupada pelo patch

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(VegType(:,:,ipatch) ,veg_type(:,:,ipatch) )
          ELSE
             CALL IJtoIBJB( VegType(:,:,ipatch) ,veg_type(:,:,ipatch) )
          END IF

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch) )
          ELSE
             CALL IJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch)  )
          END IF

       END DO

       CLOSE(UNIT=nfvegtp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfvegtp,FILE=TRIM(fNameVegType),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameVegType), ierr
          STOP "**(ERROR)**"
       END IF

       DO ipatch=1,npatches_actual
          !
          !print*,'=======================VEGET =======================',ipatch
          !
          CALL vfirec(nfvegtp,VegType(1,1,ipatch),iMax*jMax,'LIN') !veg dominante no patch
          !
          !print*,'=======================FRACA =======================',ipatch
          !
          CALL vfirec(nfvegtp,FracOcc(1,1,ipatch),iMax*jMax,'LIN') !fracao ocupada pelo patch

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(VegType(:,:,ipatch) ,veg_type(:,:,ipatch))
          ELSE
             CALL IJtoIBJB(VegType(:,:,ipatch) ,veg_type(:,:,ipatch) )
          END IF

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch) )
          ELSE
             CALL IJtoIBJB(FracOcc(:,:,ipatch) ,frac_occ(:,:,ipatch)  )
          END IF

       END DO

       CLOSE(UNIT=nfvegtp)

    END IF
    !
    ! fractional area normalization
    !
    DO j=1,jbMax
       DO i=1,ibMaxPerJB(j)
          IF(frac_occ(i,j,1) < 0.99999_r8) THEN
             fractx=0.0_r8

             DO ipatch=1,npatches_actual-1
                fractx=fractx+frac_occ(i,j,ipatch)
             END DO

             frac_occ(i,j,npatches_actual)= 1.0_r8 - fractx
          END IF

       END DO
    END DO

    !IF (reducedGrid) THEN
    !   CALL NearestIBJBtoIJ(frac_occ(:,:,npatches_actual),FracOcc(:,:,npatches_actual))
    !ELSE
    !   CALL IBJBtoIJ(frac_occ(:,:,npatches_actual),FracOcc(:,:,npatches_actual))
    !END IF
    !!
    !!-
    !!
    !DO ipatch=1,npatches_actual
    !   DO j=1,jbMax
    !      DO i=1,ibMaxPerJB(j)
    !         veg_type(i,j,ipatch)=REAL(INT(veg_type(i,j,ipatch)+0.1_r8),r8)
    !      END DO
    !   END DO
!
!       IF (reducedGrid) THEN
!          CALL NearestIBJBtoIJ(veg_type(:,:,npatches_actual),VegType(:,:,npatches_actual))
!       ELSE
!          CALL IBJBtoIJ(veg_type(:,:,npatches_actual),VegType(:,:,npatches_actual))
!       END IF
!
!    END DO
    !
    !------------------------- soil moisture initialization ------------
    !
    call MsgOne('**(read_gl_sm_bc)**','Opening GL_SM file='//TRIM(fNameSoilMoist))

    IF (record_type == 'seq') THEN !sequential mode

       OPEN(UNIT=nfslmtp,FILE=TRIM(fNameSoilMoist),FORM='unformatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilMoist), ierr
          STOP "**(ERROR)**"
       END IF

       ! do k=1,nzg   ! direct order

       DO k=nzg,1,-1 ! revert reading order
          READ(UNIT=nfslmtp) ((glsm(i,j,k),i=1,iMax),j=1,jMax) ! wetness

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          ELSE
             CALL IJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          END IF

       END DO

       CLOSE(UNIT=nfslmtp)

    ELSE IF (record_type == 'vfm') THEN !vformat model

       OPEN(UNIT=nfslmtp,FILE=TRIM(fNameSoilMoist),FORM='formatted',ACCESS='sequential',&
            ACTION='read',STATUS='old',IOSTAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
               TRIM(fNameSoilMoist), ierr
          STOP "**(ERROR)**"
       END IF


       ! do k=1,nzg   ! direct order

       DO k=nzg,1,-1 ! revert reading order
          !
          !print*,'================== GLSM for k====================',k
          !
          CALL vfirec(nfslmtp,glsm(1,1,k),iMax*jMax,'LIN')

          IF (reducedGrid) THEN
             CALL NearestIJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          ELSE
             CALL IJtoIBJB(glsm(:,:,k) ,glsm_w(:,:,k) )
          END IF

       END DO

       CLOSE(UNIT=nfslmtp)

    ELSE

       call FatalError('**(read_gl_sm_bc)** unknown record type')

    END IF

    call MsgOne('**(read_gl_sm_bc)**','DONE')

    RETURN
  END SUBROUTINE read_gl_sm_bc



  SUBROUTINE InitCheckSSiBFile(iMax,jMax,ibMax,&
       jbMax  ,kMax, ifdy  ,ids   ,idc   ,ifday , &
       tod   ,idate ,idatec   ,todsib,ibMaxPerJB )
    INTEGER      , INTENT(IN   ) :: iMax
    INTEGER      , INTENT(IN   ) :: jMax
    INTEGER      , INTENT(IN   ) :: ibMax
    INTEGER      , INTENT(IN   ) :: jbMax
    INTEGER      , INTENT(IN   ) :: kMax
    INTEGER      , INTENT(OUT  ) :: ifdy
    INTEGER      , INTENT(OUT  ) :: ids   (4)
    INTEGER      , INTENT(OUT  ) :: idc   (4)
    INTEGER      , INTENT(IN   ) :: ifday
    REAL(KIND=r8), INTENT(IN   ) :: tod
    INTEGER      , INTENT(IN   ) :: idate (4)
    INTEGER      , INTENT(IN   ) :: idatec(4)
    REAL(KIND=r8), INTENT(OUT  ) :: todsib
    INTEGER      , INTENT(IN   ) :: ibMaxPerJB(jbMax)

    INTEGER                :: j
    INTEGER                :: ncount
    INTEGER                :: i
    REAL(KIND=r8)                   :: tice  =271.16e0_r8
    CHARACTER(LEN=*), PARAMETER :: h="**(InitCheckSSiBFile)**"    
    !
    !     read ssib dataset for cold start
    !
    IF( initlz < 0 )THEN

       CALL MsgOne(h,'Read SSib variables from warm-start file')

       READ(UNIT=nfsibi) ifdy,todsib,ids,idc
       READ(UNIT=nfsibi) tm0   ,tmm
       READ(UNIT=nfsibi) qm0   ,qmm
       READ(UNIT=nfsibi) td0   ,tdm
       READ(UNIT=nfsibi) tg0   ,tgm
       READ(UNIT=nfsibi) tc0   ,tcm
       READ(UNIT=nfsibi) w0    ,wm
       READ(UNIT=nfsibi) capac0,capacm
       READ(UNIT=nfsibi) ppci  ,ppli,tkemyj
       READ(UNIT=nfsibi) gl0   ,zorl  ,gtsea,gco2flx,tseam,qsfc0,tsfc0,qsfcm,tsfcm,HML,HUML,HVML,TSK,z0sea,&
       TC_SeaIce,TGS_SeaIce,TD_SeaIce,TA_SeaIce,SNOA_SeaIce,SNOB_SeaIce,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 
       READ(UNIT=nfsibi) laymld,       hbath,     tdeep,sdeep
       Mmlen=gl0
       REWIND nfsibi

       CALL getsbc (iMax ,jMax  ,kMax, AlbVisDiff,gtsea,gco2flx,gndvi,soilm,sheleg,o3mix,tracermix,wsib3d,&
!tar begin 
!climate aerosol parameters of coarse mode
          aod,asy,ssa,z_aer,ifaeros,&
!tar end
!
!tar begin 
!climate aerosol parameters of fine mode
          aodF,asyF,ssaF,z_aerF, &
!tar end       
           ifday , tod  ,idate ,idatec,&
           ifalb,ifsst,ifco2flx,ifndvi,ifslm ,ifslmSib2,ifsnw,ifozone,iftracer, &
           sstlag,intsst,intndvi,intsoilm,fint ,tice  , &
           yrl  ,monl,ibMax,jbMax,ibMaxPerJB)

       DO j=1,jbMax
          ncount=0
          DO i=1,ibMaxPerJB(j)
             IF(imask(i,j) >= 1_i8)THEN
                ncount=ncount+1
                ssib(ncount,j)=0.0_r8
                IF(w0(ncount,1,j).LT.0.0_r8)THEN
                   ssib(ncount,  j)=ABS(w0(ncount,1,j))
                   w0  (ncount,1,j)=ABS(w0(ncount,1,j))
                   w0  (ncount,2,j)=ABS(w0(ncount,2,j))
                   w0  (ncount,3,j)=ABS(w0(ncount,3,j))
                   wm  (ncount,1,j)=ABS(wm(ncount,1,j))
                   wm  (ncount,2,j)=ABS(wm(ncount,2,j))
                   wm  (ncount,3,j)=ABS(wm(ncount,3,j))
                END IF
             END IF
          END DO
       END DO
       IF(nfctrl(5).GE.1)WRITE(UNIT=nfprt,FMT=444) ifdy,todsib,ids,idc
    END IF

444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
  END SUBROUTINE InitCheckSSiBFile


  SUBROUTINE InitSurfTemp (jbMax ,ibMaxPerJB)

    INTEGER, INTENT(IN   ) :: jbMax
    INTEGER, INTENT(IN   ) :: ibMaxPerJB(jbMax)
    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: ncount
    REAL(KIND=r8)                   :: zero  =0.0e3_r8
    REAL(KIND=r8)                   :: thousd=1.0e3_r8
    REAL(KIND=r8)                   :: tf    =273.16e0_r8
    capacm=0.0_r8
    capac0=0.0_r8
    DO j=1,jbMax
       ncount=0
       DO i=1,ibMaxPerJB(j)
          IF(imask(i,j) >= 1_i8) THEN
             ncount=ncount+1
             IF(sheleg(i,j).GT.zero) THEN
                capac0(ncount,2,j) = sheleg(i,j)/thousd
                capacm(ncount,2,j) = sheleg(i,j)/thousd
                tg0   (ncount,  j) = MIN(tg0(ncount,j),tf-0.01e0_r8)
                tgm   (ncount,  j) = MIN(tgm(ncount,j),tf-0.01e0_r8)
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE InitSurfTemp

  INTEGER FUNCTION textclass(BEE,Poros)
  IMPLICIT NONE
  REAL(KIND=r8), INTENT(IN   ) :: BEE
  REAL(KIND=r8), INTENT(IN   ) :: Poros
  INTEGER :: k
  REAL(KIND=r8) :: diff_BEE(12)
  REAL(KIND=r8) :: diff_POROS(12)

  REAL(KIND=r8) :: diff_BEE2(12)
  REAL(KIND=r8) :: diff_POROS2(12)

  REAL(KIND=r8) :: diff_TOTAL1(12)
  REAL(KIND=r8) :: diff_TOTAL2(12)

  REAL(KIND=r8) :: diffN_TOTAL1(12)
  REAL(KIND=r8) :: diffN_TOTAL2(12)

  INTEGER      :: minsoiltx1(0:1)
  INTEGER      :: minsoiltx2(0:1)
  INTEGER      :: maxsoiltx1(0:1)
  INTEGER      :: maxsoiltx2(0:1)

  REAL(KIND=r8), PARAMETER :: poros_def1(12)= RESHAPE( &
   (/0.373_r8, 0.386_r8, 0.407_r8, 0.461_r8, 0.480_r8, 0.436_r8,&
     0.416_r8, 0.423_r8, 0.449_r8, 0.476_r8, 0.480_r8, 0.465_r8/), (/12/))


  REAL(KIND=r8), PARAMETER :: poros_def2(12)= RESHAPE( &
  (/0.437_r8,0.437_r8,0.453_r8,0.501_r8,0.480_r8,0.463_r8,&
    0.398_r8,0.430_r8,0.464_r8,0.471_r8,0.479_r8,0.475_r8/), (/12/))


  REAL(KIND=r8), PARAMETER :: bee_def1(12)= RESHAPE( &
  (/  3.387_r8, 3.705_r8, 4.500_r8, 4.977_r8, 4.023_r8, 5.772_r8, &
      7.362_r8, 9.270_r8, 9.111_r8, 9.111_r8, 9.429_r8,13.245_r8/), (/12/))


  REAL(KIND=r8), PARAMETER :: bee_def2(12)= RESHAPE( &
  (/1.7_r8,2.1_r8,3.1_r8,4.7_r8,4.0_r8,4.5_r8,   &
    4.0_r8,6.0_r8,5.2_r8,6.6_r8,7.9_r8,7.6_r8/), (/12/))
    !Soil  BEE      Poros    Soil  Name            % clay  % sand
    ! 1   3.387     0.373     1    sand               3       92
    ! 2   3.705     0.386     2    loamy sand         5       82
    ! 3   4.500     0.407     3    sandy loam         10      65
    ! 4   4.977     0.461     4    silt loam          13      22
    ! 5   4.023     0.480     5    silt               7       7
    ! 6   5.772     0.436     6    loam               18      42
    ! 7   7.362     0.416     7    sandy clay loam    28      58
    ! 8   9.270     0.423     8    sandy clay         40      52
    ! 9   9.111     0.449     9    clay  loam         39      32
    !10   9.111     0.476    10    silty clay loam    39      10
    !11   9.429     0.480    11    silty clay         41      7
    !12  13.245     0.465    12    clay               65      19

    !------------------------------------------------------
    ! Soil properties data from Rawls et al. (1992)
    ! Organic properties data compiled by Mustapha El Maayar (2000)
    ! Organic FC and WP taken from Nijssen et al., 1997 (JGR; table 3 OBS-top)

    !------------------------------------------------------
    ! Variable column header definitions
    !------------------------------------------------------
    ! Sand     : sand fraction
    ! Silt     : silt fraction
    ! Clay     : clay fraction
    ! Porosity : porosity (volume fraction)
    ! FC       : field capacity (volume fraction)
    ! WP       : wilting point (volume fraction)
    ! bexp     : Campbell's 'b' exponent
    ! AEP      : air entry potential (m-H20)
    ! SHC      : saturated hydraulic conductivity (m s-1)
    !      ! dummyvarpk=0.0_r8
    !      dummyvarpk(1:108)=(/ &
    !------------------------------------------------------------------------------------------------------
    !      Sand   Silt   Clay  Porosity    FC      WP    bexp    AEP       SHC        Texture class
    !------------------------------------------------------------------------------------------------------
!  1!      0.92_r8,0.05_r8,0.03_r8,0.437_r8,0.091_r8,0.033_r8,1.7_r8,0.07_r8,5.83300E-05_r8,&  ! Sand
!  2!      0.81_r8,0.12_r8,0.07_r8,0.437_r8,0.125_r8,0.055_r8,2.1_r8,0.09_r8,1.69720E-05_r8,&  ! Loamy Sand
!  3!      0.65_r8,0.25_r8,0.10_r8,0.453_r8,0.207_r8,0.095_r8,3.1_r8,0.15_r8,7.19440E-06_r8,&  ! Sandy Loam
!  4!      0.20_r8,0.65_r8,0.15_r8,0.501_r8,0.330_r8,0.133_r8,4.7_r8,0.21_r8,1.88890E-06_r8,&  ! Silty Loam
! 5 !                              0.480_r8,                  4.0_r8                           ! Silt
!  6!      0.42_r8,0.40_r8,0.18_r8,0.463_r8,0.270_r8,0.117_r8,4.5_r8,0.11_r8,3.66670E-06_r8,&  ! Loam
!  7!      0.60_r8,0.13_r8,0.27_r8,0.398_r8,0.255_r8,0.148_r8,4.0_r8,0.28_r8,1.19440E-06_r8,&  ! Sandy Clay Loam
!  8!      0.53_r8,0.07_r8,0.40_r8,0.430_r8,0.339_r8,0.239_r8,6.0_r8,0.29_r8,3.33330E-07_r8,&  ! Sandy Clay
!  9!      0.32_r8,0.34_r8,0.34_r8,0.464_r8,0.318_r8,0.197_r8,5.2_r8,0.26_r8,6.38890E-07_r8,&  ! Clay Loam
! 10!      0.09_r8,0.58_r8,0.33_r8,0.471_r8,0.366_r8,0.208_r8,6.6_r8,0.33_r8,4.16670E-07_r8,&  ! Silty Clay Loam
! 11!      0.10_r8,0.45_r8,0.45_r8,0.479_r8,0.387_r8,0.250_r8,7.9_r8,0.34_r8,2.50000E-07_r8,&  ! Silty Clay
! 12!      0.20_r8,0.20_r8,0.60_r8,0.475_r8,0.396_r8,0.272_r8,7.6_r8,0.37_r8,1.66670E-07_r8,&  ! Clay
    !====================================================================

    !Call function that estimates texture and put into structure:

  DO k=1,12
     diff_BEE (k)  =ABS(BEE  -bee_def1  (k))
     diff_BEE2(k)  =ABS(BEE  -bee_def2  (k))
     diff_TOTAL1(k)= SQRT((diff_BEE(k)**2)+(diff_BEE2(k)**2))

     diff_POROS (k)=ABS(Poros-poros_def1(k))
     diff_POROS2(k)=ABS(Poros-poros_def2(k))
     diff_TOTAL2(k)= SQRT((diff_POROS(k)**2)+(diff_POROS2(k)**2))


  END DO
  minsoiltx1=MINLOC(diff_TOTAL1,dim=1)
  minsoiltx2=MINLOC(diff_TOTAL2,dim=1)

  maxsoiltx1=MAXLOC(diff_TOTAL1,dim=1)
  maxsoiltx2=MAXLOC(diff_TOTAL2,dim=1)

  IF(minsoiltx1(1) /= minsoiltx2(1))THEN
    DO k=1,12
       !Normalizacao segundo a amplitude
       !Justificativa: unidades diferentes ou dispersoes muito heterogeneas 
       diffN_TOTAL1(k)=(diff_TOTAL1(k)-diff_TOTAL1(minsoiltx1(1)))/(diff_TOTAL1(maxsoiltx1(1))-diff_TOTAL1(minsoiltx1(1)))
       diffN_TOTAL2(k)=(diff_TOTAL2(k)-diff_TOTAL2(minsoiltx2(1)))/(diff_TOTAL2(maxsoiltx2(1))-diff_TOTAL2(minsoiltx2(1)))

    END DO
    minsoiltx1=MINLOC(diffN_TOTAL1,dim=1)
    minsoiltx2=MINLOC(diffN_TOTAL2,dim=1)


    IF(diffN_TOTAL1(minsoiltx1(1))+ diffN_TOTAL2(minsoiltx1(1)) < diffN_TOTAL1(minsoiltx2(1))+ diffN_TOTAL2(minsoiltx2(1)) )THEN
      textclass=minsoiltx1(1)
    ELSE
      textclass=minsoiltx2(1)
    END IF
  ELSE
      textclass=minsoiltx1(1)
  END IF
  
  END FUNCTION textclass


  SUBROUTINE ReStartSSiB (jbMax,ifday,tod,idate ,idatec, &
       nfsibo,ibMaxPerJB)

    INTEGER           ,INTENT(IN   ) :: jbMax
    INTEGER           ,INTENT(IN   ) :: ifday
    REAL(KIND=r8)     ,INTENT(IN   ) :: tod
    INTEGER           ,INTENT(IN   ) :: idate(4)
    INTEGER           ,INTENT(IN   ) :: idatec(4)
    INTEGER           ,INTENT(IN   ) :: nfsibo
    INTEGER           ,INTENT(IN   ) :: ibMaxPerJB(jbMax)
    INTEGER                         :: i
    INTEGER                         :: j
    INTEGER                         :: ncount

    IF(TRIM(isimp).NE.'YES') THEN

       CALL MsgOne('**(restartphyscs)**','Saving physics state for restart')

       !$OMP DO PRIVATE(ncount, i)
       DO j=1,jbMax
          ncount=0
          DO i=1,ibMaxPerJB(j)
             IF(imask(i,j) >= 1_i8)THEN
                ncount=ncount+1
                IF(ssib(ncount,j).GT.0.0_r8)THEN
                   w0  (ncount,1,j)=-ssib(ncount,j)
                   w0  (ncount,2,j)=-ssib(ncount,j)
                   w0  (ncount,3,j)=-ssib(ncount,j)
                   wm  (ncount,1,j)=-ssib(ncount,j)
                   wm  (ncount,2,j)=-ssib(ncount,j)
                   wm  (ncount,3,j)=-ssib(ncount,j)
                END IF
             END IF
          END DO
       END DO


       !$OMP SINGLE
       WRITE(UNIT=nfsibo) ifday,tod,idate,idatec
       WRITE(UNIT=nfsibo) tm0,tmm
       WRITE(UNIT=nfsibo) qm0,qmm
       WRITE(UNIT=nfsibo) td0,tdm
       WRITE(UNIT=nfsibo) tg0,tgm
       WRITE(UNIT=nfsibo) tc0,tcm
       WRITE(UNIT=nfsibo) w0 ,wm
       WRITE(UNIT=nfsibo) capac0,capacm
       WRITE(UNIT=nfsibo)  ppci  ,ppli,tkemyj
       WRITE(UNIT=nfsibo) gl0 ,zorl,gtsea,gco2flx,tseam,qsfc0,tsfc0,qsfcm,tsfcm,HML,HUML,HVML,TSK,z0sea,&
       TC_SeaIce,TGS_SeaIce,TD_SeaIce,TA_SeaIce,SNOA_SeaIce,SNOB_SeaIce,PBL_CoefKm, PBL_CoefKh,tauresx,tauresy,poda,tmin2m   ,tmax2m 
       WRITE(UNIT=nfsibo) laymld,       hbath,     tdeep,sdeep
       !$OMP END SINGLE
    END IF

  END SUBROUTINE ReStartSSiB

  SUBROUTINE Finalize_SSiB()
    DEALLOCATE(vcover_gbl)
    DEALLOCATE(zlt_gbl   )
    DEALLOCATE(green_gbl )
    DEALLOCATE(chil_gbl  )
    DEALLOCATE(topt_gbl  )
    DEALLOCATE(tll_gbl   )
    DEALLOCATE(tu_gbl    )
    DEALLOCATE(defac_gbl )
    DEALLOCATE(ph2_gbl   )
    DEALLOCATE(ph1_gbl   )
    DEALLOCATE(rstpar_gbl)
    DEALLOCATE(satcap3d  )
    DEALLOCATE(extk_gbl  )
    DEALLOCATE(radfac_gbl)
    DEALLOCATE(closs_gbl )
    DEALLOCATE(gloss_gbl )
    DEALLOCATE(thermk_gbl)
    DEALLOCATE(p1f_gbl   )
    DEALLOCATE(p2f_gbl   )
    DEALLOCATE(tgeff_gbl )
    DEALLOCATE(zlwup_SSiB)
    DEALLOCATE(AlbGblSSiB)
    DEALLOCATE(MskAntSSiB)
    DEALLOCATE(glsm_w    )
    DEALLOCATE(veg_type  )
   
  END SUBROUTINE Finalize_SSiB

END MODULE SFC_SSiB
