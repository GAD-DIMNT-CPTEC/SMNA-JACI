MODULE MicroPhysics
  USE Constants, ONLY :  &
       delq             ,&
       r8,i8,qmin,grav,gasr

  USE Micro_Hack, ONLY:       &
      Init_Micro_Hack ,RunMicro_Hack   
 
  USE Micro_GTHOMPSON, ONLY:       &
      Init_Micro_thompson,RunMicro_thompson   
 
  USE Micro_Ferrier, ONLY:       &
      Init_Micro_Ferrier ,RunMicro_FERRIER

  USE Micro_UKME, ONLY:       &
      Init_Micro_UKME ,RunMicro_UKME

  USE Micro_MORR, ONLY:       &
      Init_Micro_MORR ,RunMicro_MORR

  USE Micro_HugMorr, ONLY:       &
      Init_Micro_HugMorr ,RunMicro_HugMorr

  USE Micro_LrgScl, ONLY:     &
      Init_Micro_LrgScl,RunMicro_LrgScl 

  USE FieldsPhysics, ONLY:  &
       LOWLYR             , &
       F_ICE_PHY          , &
       F_RAIN_PHY         , &
       F_RIMEF_PHY       


   IMPLICIT NONE
  SAVE

  PRIVATE

  INTEGER, PARAMETER :: ppcnst=3
  REAL(KIND=r8), ALLOCATABLE  :: ql2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: ql3(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qi3(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qr2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qr3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: qs2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qs3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: qg2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: qg3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NI2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NI3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NS2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NS3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: nifa2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: nifa3(:,:,:)  

  REAL(KIND=r8), ALLOCATABLE  :: nwfa2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: nwfa3(:,:,:)  
  
  REAL(KIND=r8), ALLOCATABLE  :: NR2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NR3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NG2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NG3(:,:,:)

  REAL(KIND=r8), ALLOCATABLE  :: NC2(:,:,:)
  REAL(KIND=r8), ALLOCATABLE  :: NC3(:,:,:)


  PUBLIC InitMicroPhysics
  PUBLIC RunMicroPhysics
  PUBLIC FinalizeMicroPhysics
CONTAINS 
  SUBROUTINE InitMicroPhysics( &
                          dt      ,a_hybr    ,b_hybr    ,restart,&
                          kMax    ,iMax      ,jMax      ,ibMax      , &
                          jbMax   ,fNameMicro,path_in   , ILCON    ,microphys  , &
                          nClass  ,nAeros    ,EFFCS     ,EFFIS)
    IMPLICIT NONE
    REAL(KIND=r8)   , INTENT(IN   ) :: dt
    INTEGER         , INTENT(IN   ) :: kMax
    INTEGER         , INTENT(IN   ) :: iMax
    INTEGER         , INTENT(IN   ) :: jMax
    INTEGER         , INTENT(IN   ) :: ibMax
    INTEGER         , INTENT(IN   ) :: jbMax
    REAL(KIND=r8)   , INTENT(IN   ) :: a_hybr (kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: b_hybr (kMax+1)
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameMicro
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    LOGICAL         , INTENT(IN   ) :: microphys
    CHARACTER(LEN=*), INTENT(IN   ) :: path_in
    INTEGER         , INTENT(IN   ) :: nClass
    INTEGER         , INTENT(IN   ) :: nAeros
    REAL(KIND=r8)   , INTENT(OUT  ) :: EFFCS(ibMax,kMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: EFFIS(ibMax,kMax,jbMax)
    CHARACTER(LEN=*), PARAMETER     :: h='**(InitMicroPhysics)**'
    LOGICAL         , INTENT(IN   )  :: restart

    IF(TRIM(ILCON).EQ.'YES'.or.TRIM(ILCON).EQ.'LSC') CALL Init_Micro_LrgScl()    
    IF(TRIM(ILCON).EQ.'MIC') CALL Init_Micro_Hack(kMax,jMax,ibMax,jbMax,ppcnst,a_hybr,b_hybr)

    IF(TRIM(ILCON).EQ.'HWRF') CALL Init_Micro_thompson(kMax,path_in,restart,a_hybr,b_hybr)

    IF(TRIM(ILCON).EQ.'HGFS') CALL Init_Micro_Ferrier (dt,iMax,jMax,ibmax,kMax,jbMax,F_ICE_PHY,F_RAIN_PHY ,F_RIMEF_PHY)
    IF(TRIM(ILCON).EQ.'UKMO') CALL Init_Micro_UKME(iMax)
    IF(TRIM(ILCON).EQ.'MORR') CALL Init_Micro_MORR ()
    IF(TRIM(ILCON).EQ.'HUMO') CALL Init_Micro_HugMorr (ibMax,kMax,jbMax,EFFCS,EFFIS)

    IF (microphys) THEN
       ALLOCATE(ql3(ibMax,kMax,jbMax));ql3=0.00001e-12_r8
       ALLOCATE(qi3(ibMax,kMax,jbMax));qi3=0.00001e-12_r8

       ALLOCATE(ql2(ibMax,kMax,jbMax));ql2=0.00001e-12_r8
       ALLOCATE(qi2(ibMax,kMax,jbMax));qi2=0.00001e-12_r8
       IF((nClass+nAeros)>0)THEN
          ALLOCATE(qr3(ibMax,kMax,jbMax));qr3=0.00001e-12_r8
          IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO' .or.TRIM(ILCON).EQ.'HWRF') THEN
             ALLOCATE(qs3(ibMax,kMax,jbMax));qs3=0.00001e-12_r8
             ALLOCATE(qg3(ibMax,kMax,jbMax));qg3=0.00001e-12_r8
             ALLOCATE(NI3(ibMax,kMax,jbMax));NI3=0.00001e-12_r8
             IF( TRIM(ILCON).EQ.'HWRF')THEN
                ALLOCATE(nifa3(ibMax,kMax,jbMax));nifa3=0.00001e-12_r8
             ELSE      
                ALLOCATE(NS3  (ibMax,kMax,jbMax));NS3=0.00001e-12_r8
             END IF
             ALLOCATE(NR3(ibMax,kMax,jbMax));NR3=0.00001e-12_r8

             IF( TRIM(ILCON).EQ.'HWRF')THEN
                ALLOCATE(nwfa3(ibMax,kMax,jbMax));nwfa3=0.00001e-12_r8
             ELSE      
                ALLOCATE(NG3(ibMax,kMax,jbMax));NG3=0.00001e-12_r8
             END IF

             IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF')THEN
                ALLOCATE(NC3(ibMax,kMax,jbMax));NC3=0.00001e-12_r8
             END IF
          END IF
          ALLOCATE(qr2(ibMax,kMax,jbMax));qr2=0.00001e-12_r8
          IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF') THEN
             ALLOCATE(qs2(ibMax,kMax,jbMax));qs2=0.00001e-12_r8
             ALLOCATE(qg2(ibMax,kMax,jbMax));qg2=0.00001e-12_r8
             ALLOCATE(NI2(ibMax,kMax,jbMax));NI2=0.00001e-12_r8
             IF( TRIM(ILCON).EQ.'HWRF')THEN
                ALLOCATE(nifa2(ibMax,kMax,jbMax));nifa2=0.00001e-12_r8
             ELSE
                ALLOCATE(NS2(ibMax,kMax,jbMax));NS2=0.00001e-12_r8
             END IF
             ALLOCATE(NR2(ibMax,kMax,jbMax));NR2=0.00001e-12_r8
             IF( TRIM(ILCON).EQ.'HWRF')THEN
                ALLOCATE(nwfa2(ibMax,kMax,jbMax));nwfa2=0.00001e-12_r8
             ELSE      
                ALLOCATE(NG2(ibMax,kMax,jbMax));NG2=0.00001e-12_r8
             END IF
             IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF') THEN 
                ALLOCATE(NC2(ibMax,kMax,jbMax));NC2=0.00001e-12_r8
             END IF
          END IF
       END IF 
    END IF

  END SUBROUTINE InitMicroPhysics
!**********************************************************************************
!**********************************************************************************


  SUBROUTINE RunMicroPhysics(&
      ! Run Flags
               mlrg        , & !INTEGER      , INTENT(in   ) :: mlrg
               ILCON       , & !CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
               microphys   , & !LOGICAL      , INTENT(in   ) :: microphys
      ! Time info
               jdt         , & !INTEGER      , INTENT(in   ) :: jdt
               dt          , & !REAL(KIND=r8), INTENT(in   ) :: dt
      ! Model Geometry
               colrad      , & !REAL(KIND=r8)   , INTENT(IN   ) :: colrad  (iMax)  
               LOWLYR      , & !INTEGER         , INTENT(IN   ) :: LOWLYR  (iMax)  
               terr        , & !REAL(KIND=r8)   , INTENT(IN   ) :: terr (iMax)  
      ! Model information
               mask        , & !INTEGER(KIND=i8), INTENT(IN   ) :: mask (iMax) 
               nClass      , & !INTEGER      , INTENT(in   ) :: nClass
               nAeros      , &
               iMax        , & !INTEGER      , INTENT(in   ) :: iMax
               kMax        , & !INTEGER      , INTENT(in   ) :: kMax
               latco       , & !INTEGER      , INTENT(in   ) :: latco
      ! Surface field
               tsfc        , & !REAL(KIND=r8)   , INTENT(IN   ) :: tsfc (iMax)
      ! PBL field
               PBL_CoefKh  , & !REAL(KIND=r8)   , INTENT(IN   ) :: PBL_CoefKh(1:iMax,1:kMax)         , &
               tke         , & !REAL(KIND=r8)   , INTENT(IN   ) :: tke(1:iMax,1:kMax)         , &
               pblh        , & !REAL(KIND=r8)   , INTENT(IN   ) :: pblh (iMax)  
               var         , &
      ! CONVECTION: Cloud field
               kuo         , & !INTEGER         , INTENT(in   ) :: kuo (iMax)  
               cmfmc       , & !REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc   (iMax,kMax+1)   ! convective mass flux--m sub c
               cmfmc2      , & !REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc2  (iMax,kMax+1)   ! shallow convective mass flux--m sub c
               dlf         , & !REAL(KIND=r8)   , INTENT(IN   ) :: dlf (iMax,kMax)! detrained water from ZM
               rliq        , & !REAL(KIND=r8)   , INTENT(IN   ) :: rliq (iMax)        ! vertical integral of liquid not yet in q(ixcldliq)
               concld      , & !REAL(KIND=r8)   , INTENT(INOUT) :: concld  (iMax,kMax)
               cld         , & !REAL(KIND=r8)   , INTENT(INOUT) :: cld     (iMax,kMax)
               cldtot      , &!(1:ibLim,:,latco)
               fdqn        , & !REAL(KIND=r8)   , INTENT(inout) :: fdqn (iMax,kMax)
               EFFCS       , & !REAL(KIND=r8)   , INTENT(INOUT) :: EFFCS(1:iMax,1:kMax)         , &
               EFFIS       , & !REAL(KIND=r8)   , INTENT(INOUT) :: EFFIS(1:iMax,1:kMax)         , &
               Total_Rain  , & !REAL(KIND=r8)   , INTENT(inout) :: Total_Rain  (iMax)
               Total_Snow  , & !REAL(KIND=r8)   , INTENT(inout) :: Total_Snow  (iMax)
               RAINNCV     , & !REAL(KIND=r8)   , INTENT(INOUT) :: RAINNCV     (iMax)
               SNOWNCV     , & !REAL(KIND=r8)   , INTENT(INOUT) :: SNOWNCV     (iMax)
               F_ICE_PHY   , & !REAL(KIND=r8)   , INTENT(INOUT) :: F_ICE_PHY   (1:iMax,1:kMax)
               F_RAIN_PHY  , & !REAL(KIND=r8)   , INTENT(INOUT) :: F_RAIN_PHY  (1:iMax,1:kMax)
               F_RIMEF_PHY , & !REAL(KIND=r8)   , INTENT(INOUT) :: F_RIMEF_PHY (1:iMax,1:kMax)
      ! Atmospheric fields
               prsi        , & 
               prsl        , & 
               t2          , & !REAL(KIND=r8)   , INTENT(IN   ) :: t2(iMax,kMax)
               t3          , & !REAL(KIND=r8)   , INTENT(inout) :: t3 (iMax,kMax) 
               q2          , & !REAL(KIND=r8)   , INTENT(IN   ) :: q2 (iMax,kMax)
               q3          , & !REAL(KIND=r8)   , INTENT(inout) :: q3 (iMax,kMax)
               ub          , & !REAL(KIND=r8)   , INTENT(IN   ) :: ub (iMax,kMax) ! (m/s) 
               vb          , & !REAL(KIND=r8)   , INTENT(IN   ) :: vb (iMax,kMax) ! (m/s)
               omgb        , & !REAL(KIND=r8)   , INTENT(IN   ) :: omgb (iMax,kMax) ! (Pa/s)
               dq          , & !REAL(KIND=r8)   , INTENT(inout) :: dq (iMax,kMax)
               tLrgs       , & !REAL(KIND=r8)   , INTENT(OUT  ) :: tLrgs       (1:iMax,1:kMax)
               qLrgs       , & !REAL(KIND=r8)   , INTENT(OUT  ) :: qLrgs       (1:iMax,1:kMax)
      ! Microphysics
               dtdt        , & 
               dqdt        , & 
               dqldt       , & 
               dqidt       , & 
               gicem       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gicem       (iMax,kmax)
               gicep       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gicep       (iMax,kmax)
               gicet       , & 
               gliqm       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gliqm       (iMax,kmax)
               gliqp       , & !REAL(KIND=r8)   , INTENT(INOUT) :: gliqp       (iMax,kmax)
               gliqt       , &
               gvarm       , & !REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)
               gvarp       , & !REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)
               gvart         ) 
    !************************************************************************
    !   The cumulus_driver subroutine calls deep and shallow cumulus
    !   parameterization schemes.
    !   more information nilo@cptec.inpe.br
    !   NOTE: This version is not official. You can use only for test.
    !************************************************************************
    !
    !  Definition/
    !---------------
    !             I N P U T  O U T P U T  F O R   G C M
    !             -------------------------------------
    ! INPUT
    !
    !** integer
    !    iMax                   ! end index for longitude domain
    !    kMax                   ! end index for u,v,t,p sigma levels
    !    jdt                    ! number of time step
    !   kuo                     ! convection yes(1) or not(0) for shallow convection
    !
    !** real
    !    dt                     ! time step (s)
    !    ta                     ! temperature (K) at time t-1
    !    tb                     ! temperature (K) at time t
    !    tc                     ! temperature (K) at time t+1
    !    qa                     ! water vapor mixing ratio (kg/kg) at time t-1
    !    qb                     ! water vapor mixing ratio (kg/kg) at time t
    !    qc                     ! water vapor mixing ratio (kg/kg) at time t+1
    !    psb                    ! surface pressure (cb)     at time t
    !    ub                     ! u-velocity (m/s) at time t
    !    vb                     ! v-velocity (m/s) at time t
    !    omgb                   ! vertical omega velocity (Pa/s) at time t
    !                           ! it is in half levels along with U,V,T,Q
    !    sl                     ! half sigma layers
    !    si                     ! full sigma layers
    !    del                    ! difference between full sigma layers
    !    xland                  ! land-sea mask (1 for land; 0 for water)
    !    zs                     ! topography (m)
    !    DX                     ! horizontal space interval (m)
    !    qrem,cldm              ! local variables for  RAS-Scheme
    !
    !    hrem,qrem              ! these arrays are needed for the heating 
    !                           ! and mostening from ras  scheme
    !
    !
    !    ktops, kbots           ! these arrays are needed for the new 
    !                           ! shallow convection scheme
    !    cldm                   ! needed for cloud fraction based on mass 
    !                           ! flux
    !    noshal1, kctop1, kcbot1! needed for cloud fraction based on mass 
    !                           ! flux new arrays needed for shallow 
    !                           ! convective clouds
    !     
    !
    !
    !   OUTPUT
    !**  integer
    !    kuo                    ! indix for shalow convection KUO,RAS,KUOG, GRELL
    !    ktop                   ! level of convective cloud top
    !    kbot                   ! level of convective cloud base
    !    plcl                   ! pressure at convection levl for shallow convection
    !                           ! in Kuo 
    !
    !** real
    !   RAINCV                  ! cumulus scheme precipitation (mm)
    !   tc                      ! new temperature (K) at time t+1  after CU precip
    !   qc                      ! new  water vapor mixing ratio (kg/kg) at time t+1.
    !
    !
    !*********************************************************************************
    IMPLICIT NONE
    !              I N P U T     O U T P U T    V A R I A B L E S
    !              ----------------------------------------------
    !              Xa at t-1   Xb at t    Xc at t+1

    INTEGER         , INTENT(in   ) :: jdt
    INTEGER         , INTENT(in   ) :: nClass
    INTEGER         , INTENT(IN   ) :: nAeros
    INTEGER         , INTENT(in   ) :: iMax
    INTEGER         , INTENT(in   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: dt
    INTEGER         , INTENT(in   ) :: mlrg
    INTEGER         , INTENT(in   ) :: latco
    LOGICAL         , INTENT(in   ) :: microphys
    CHARACTER(LEN=*), INTENT(IN   ) :: ILCON
    INTEGER(KIND=i8), INTENT(IN   ) :: mask    (iMax) 
    REAL(KIND=r8)   , INTENT(IN   ) :: tsfc    (iMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: terr    (iMax)  
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad  (iMax)  
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh    (iMax)  
    REAL(KIND=r8)   , INTENT(IN   ) :: var     (iMax)  

    INTEGER         , INTENT(in   ) :: kuo     (iMax)  
    INTEGER         , INTENT(IN   ) :: LOWLYR  (iMax)  
    REAL(KIND=r8)   , INTENT(in   ) :: prsi    (iMax,kMax+1)  
    REAL(KIND=r8)   , INTENT(in   ) :: prsl    (iMax,kMax  ) 
    REAL(KIND=r8)   , INTENT(inout) :: t2      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(inout) :: t3      (iMax,kMax) 
    REAL(KIND=r8)   , INTENT(inout) :: q2      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(inout) :: q3      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ub      (iMax,kMax) ! (m/s) 
    REAL(KIND=r8)   , INTENT(IN   ) :: vb      (iMax,kMax) ! (m/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: omgb    (iMax,kMax) ! (Pa/s)
    REAL(KIND=r8)   , INTENT(inout) :: dq      (iMax,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc   (iMax,kMax+1)   ! convective mass flux--m sub c
    REAL(KIND=r8)   , INTENT(IN   ) :: cmfmc2  (iMax,kMax+1)   ! shallow convective mass flux--m sub c
    REAL(KIND=r8)   , INTENT(IN   ) :: dlf     (iMax,kMax)    ! detrained water from ZM
    REAL(KIND=r8)   , INTENT(IN   ) :: rliq    (iMax)        ! vertical integral of liquid not yet in q(ixcldliq)
    REAL(KIND=r8)   , INTENT(INOUT) :: concld  (iMax,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: cld     (iMax,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cldtot  (iMax,kMax)
    REAL(KIND=r8)   , INTENT(inout) :: fdqn    (iMax,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFCS      (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: EFFIS      (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: PBL_CoefKh (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: tke(1:iMax,1:kMax) 
    REAL(KIND=r8)   , INTENT(inout) :: Total_Rain  (iMax)
    REAL(KIND=r8)   , INTENT(inout) :: Total_Snow  (iMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: RAINNCV     (iMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: SNOWNCV     (iMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: F_ICE_PHY   (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: F_RAIN_PHY  (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: F_RIMEF_PHY (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tLrgs       (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qLrgs       (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gicem       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gicep       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gicet       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gliqm       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gliqp       (iMax,kmax)
    REAL(KIND=r8)   , INTENT(INOUT) :: gliqt       (iMax,kmax)
    REAL(KINd=r8)   , INTENT(OUT  ) :: dtdt (iMax,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqdt (iMax,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqldt(iMax,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqidt(iMax,kMax)



    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gvart (iMax,kmax,nClass+nAeros)






    REAL(KIND=r8)    :: dqrdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dqsdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dqgdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dnidt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dnsdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dnrdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dNGdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dNCdt    (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dnifadt  (1:iMax, 1:kMax)
    REAL(KIND=r8)    :: dnwfadt  (1:iMax, 1:kMax)
    ! WIND COMPONENTS FOR GRELL ENSEMBLE
    REAL(KIND=r8) :: u2(iMax,kMax)
    REAL(KIND=r8) :: v2(iMax,kMax)
    REAL(KIND=r8) :: w2(iMax,kMax)
    REAL(KIND=r8) :: prec_str(iMax)  ! [Total] sfc flux of precip from stratiform (m/s) 
    REAL(KIND=r8) :: snow_str(iMax)  ! [Total] sfc flux of snow from stratiform   (m/s)
    REAL(KIND=r8) :: prec_sed(iMax)  ! surface flux of total cloud water from sedimentation
    REAL(KIND=r8) :: snow_sed(iMax)  ! surface flux of cloud ice from sedimentation
    REAL(KIND=r8) :: prec_pcw(iMax)  ! sfc flux of precip from microphysics(m/s)
    REAL(KIND=r8) :: snow_pcw(iMax)  ! sfc flux of snow from microphysics (m/s)
    REAL(KIND=r8) :: icefrac (iMax)
    REAL(KIND=r8) :: landfrac(iMax)
    REAL(KIND=r8) :: ocnfrac (iMax)
    REAL(KIND=r8) :: landm   (iMax)  ! land fraction ramped over water
    REAL(KIND=r8) :: snowh   (iMax)  ! Snow depth over land, water equivalent (m)
    REAL(KIND=r8) :: ts      (iMax)      ! surface temperature
    REAL(KIND=r8) :: sst     (iMax)       !sea surface temperature
    REAL(KIND=r8) :: tv   (iMax,kMax)
    REAL(KIND=r8) :: press(iMax,kMax)
    REAL(KIND=r8) :: RHO  (iMax,kMax)
    REAL(KIND=r8) :: delz (iMax,kMax)
    REAL(KIND=r8) :: DeltaP (iMax,kMax)
    REAL(KIND=r8) :: r1000
    REAL(KIND=r8) :: rbyg
    INTEGER       :: i
    INTEGER       :: k



    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
             ql3     (i,k,latco) = gliqp(i,k)
             qi3     (i,k,latco) = gicep(i,k)

             ql2     (i,k,latco) = gliqm(i,k) 
             qi2     (i,k,latco) = gicem(i,k)
             IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
                qr3     (i,k,latco) = gvarp(i,k,1)
                !gvart   (i,k,1)  =0.0_r8
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF') THEN
                   qs3     (i,k,latco) = gvarp(i,k,2)
                !   gvart   (i,k,2)  =0.0_r8
                   qg3     (i,k,latco) = gvarp(i,k,3)
                !   gvart   (i,k,3)  =0.0_r8
                   NI3     (i,k,latco) = gvarp(i,k,4)
                !   gvart   (i,k,4)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nifa3   (i,k,latco) = gvarp(i,k,5)
                !      gvart   (i,k,5)  =0.0_r8
                   ELSE
                      NS3     (i,k,latco) = gvarp(i,k,5)
                !      gvart   (i,k,5)  =0.0_r8
                   END IF
                   NR3     (i,k,latco) = gvarp(i,k,6)
                !   gvart   (i,k,6)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nwfa3   (i,k,latco) = gvarp(i,k,7)
                !     gvart   (i,k,7)  =0.0_r8
                   ELSE
                      NG3     (i,k,latco) = gvarp(i,k,7)
                !      gvart   (i,k,7)  =0.0_r8
                   END IF

                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF') THEN
                      NC3     (i,k,latco) = gvarp(i,k,8)
                !      gvart   (i,k,8)  =0.0_r8
                   END IF 
                END IF
                qr2     (i,k,latco) = gvarm(i,k,1)
                dqrdt   (i,k)  =0.0_r8
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF') THEN
                   qs2     (i,k,latco) = gvarm(i,k,2)
                   dqsdt   (i,k)  =0.0_r8
                   qg2     (i,k,latco) = gvarm(i,k,3)
                   dqgdt   (i,k)  =0.0_r8
                   NI2     (i,k,latco) = gvarm(i,k,4)
                   dnidt   (i,k)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nifa2   (i,k,latco) = gvarm(i,k,5)
                      dnifadt   (i,k)  =0.0_r8
                   ELSE
                      NS2     (i,k,latco) = gvarm(i,k,5)
                      dnsdt   (i,k)  =0.0_r8
                   END IF
                   NR2     (i,k,latco) = gvarm(i,k,6)
                   dnrdt   (i,k)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nwfa2   (i,k,latco) = gvarm(i,k,7)
                      dnwfadt   (i,k)  =0.0_r8
                   ELSE
                      NG2     (i,k,latco) = gvarm(i,k,7)
                      dngdt   (i,k)  =0.0_r8
                   END IF
                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF') THEN
                      NC2     (i,k,latco) = gvarm(i,k,8)
                      dncdt   (i,k)  =0.0_r8
                   END IF
                END IF
             END IF 
          END DO
       END DO
    ELSE
       DO k=1,kMax
          DO i=1,iMax
             !ql3     (i,k,latco) = 0.0_r8
             !qi3     (i,k,latco) = 0.0_r8

             !ql2     (i,k,latco) = 0.0_r8
             !qi2     (i,k,latco) = 0.0_r8
             IF((nClass+nAeros)>0 )THEN
                qr3     (i,k,latco) = 0.0_r8 
                dqrdt   (i,k)  =0.0_r8
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF') THEN
                   qs3     (i,k,latco) = 0.0_r8 
                   dqsdt   (i,k)  =0.0_r8
                   qg3     (i,k,latco) = 0.0_r8 
                   dqgdt   (i,k)  =0.0_r8
                   NI3     (i,k,latco) = 0.0_r8
                   dnidt   (i,k)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nifa3   (i,k,latco) = 0.0_r8
                      dnifadt   (i,k)  =0.0_r8
                   ELSE
                      NS3     (i,k,latco) = 0.0_r8
                      dnsdt   (i,k)  =0.0_r8
                   END IF
                   NR3     (i,k,latco) = 0.0_r8
                   dnrdt   (i,k)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nwfa3   (i,k,latco) = 0.0_r8
                      dnwfadt   (i,k)  =0.0_r8
                   ELSE
                      NG3     (i,k,latco) = 0.0_r8
                      dngdt   (i,k)  =0.0_r8
                   END IF

                   IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF')THEN
                      NC3     (i,k,latco) = 0.0_r8
                      dncdt   (i,k)  =0.0_r8
                   END IF
                END IF
                qr2     (i,k,latco) = 0.0_r8
                dqrdt   (i,k)  =0.0_r8
                IF( TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF') THEN
                   qs2     (i,k,latco) = 0.0_r8
                   dqsdt   (i,k)  =0.0_r8
                   qg2     (i,k,latco) = 0.0_r8
                   dqgdt   (i,k)  =0.0_r8
                   NI2     (i,k,latco) = 0.0_r8
                   dnidt   (i,k)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nifa2     (i,k,latco) = 0.0_r8
                      dnifadt   (i,k)  =0.0_r8
                   ELSE
                      NS2     (i,k,latco) = 0.0_r8
                      dnsdt   (i,k)  =0.0_r8
                   END IF 
                   NR2     (i,k,latco) = 0.0_r8
                   dnrdt   (i,k)  =0.0_r8
                   IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nwfa2   (i,k,latco) = 0.0_r8
                      dnwfadt   (i,k)  =0.0_r8
                   ELSE
                      NG2     (i,k,latco) = 0.0_r8
                      dngdt   (i,k)  =0.0_r8
                   END IF

                  IF( TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF') THEN
                     NC2     (i,k,latco) = 0.0_r8
                     dncdt   (i,k)  =0.0_r8
                  END IF
                END IF
             END IF
          END DO
       END DO
    END IF
    DO k=1,kMax
      DO i=1,iMax
          DeltaP(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO

    !-----------------------------------------------------------------
    ! Large Scale Precipitation
    !-----------------------------------------------------------------
    IF(TRIM(ILCON).EQ.'LSC' .OR. TRIM(ILCON).EQ.'YES' ) THEN
      CALL RunMicro_LrgScl(Total_Rain, t3, dq, q3,prsi,prsl,dtdt,dqdt, dt, &
                mlrg, latco, iMax, kMax)
      CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
    ELSE IF(TRIM(ILCON).EQ.'MIC') THEN
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
        snowh   =0.0_r8
        icefrac =0.0_r8
        landfrac=0.0_r8
        ocnfrac =0.0_r8
        DO i=1,iMax
           IF(mask(i) == 13_i8 .or. mask(i) == 15_i8 )snowh(i)=5.0_r8
           IF(mask(i).GT.0_i8)THEN
              ! land
              icefrac(i)=0.0_r8
              landfrac(i)=1.0_r8
              ocnfrac(i)=0.0_r8
           ELSE
              ! water/ocean
              landfrac(i)=0.0_r8
              ocnfrac(i) =1.0_r8
              IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                 icefrac(i)=1.0_r8
                 ocnfrac(i) =0.0_r8
              ENDIF
           END IF
        END DO
        DO k=1,kMax
           DO i=1,iMax
              u2(i,k)=ub  (i,k)
              v2(i,k)=vb  (i,k)
              w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
           END DO
        END DO
!----------------------------
        DO k=1,kMax
           DO i=1,iMax
              !press(i,k)=ps2(i)*1000_r8*sl(k)
              press(i,k)=prsl(i,k)
              tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
           END DO
        END DO    
        !
        !  Calculate the distance between the surface and the first layer of the model
        !
        r1000=1.0e0_r8 /gasr
        DO i=1,iMax
           rbyg=gasr/grav*DeltaP(i,1)*0.5e0_r8
           RHO     (i,1)=r1000*(prsl(i,1))/t2(i,1)
           delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
        END DO 

        DO k=2,kMax
           DO i=1,iMax
              RHO (i,k)=r1000*(prsl(i,k))/t2(i,k)
              delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                        LOG(press(i,k-1)/press(i,k))/grav
           END DO
        END DO
!--------------------------------------------
       landm          =0.0_r8  
       ts             =tsfc
       sst            =tsfc
       CALL RunMicro_Hack( &
       jdt                              , &! INTEGER , INTENT(in)  :: ibMax                     !number of columns (max)
       iMax                             , &! INTEGER , INTENT(in)  :: ibMax                     !number of columns (max)
       kMax                             , &! INTEGER , INTENT(in)  :: kMax                      !number of vertical levels
       kMax+1                           , &! INTEGER , INTENT(in)  :: kMax+1                    !number of vertical levels + 1
       ppcnst                           , &! INTEGER , INTENT(in)  :: ppcnst                    !number of constituent
       latco                            , &! INTEGER , INTENT(in)  :: latco                     !latitude
       dt                               , &! REAL(r8), INTENT(in)  :: dtime                     !timestep
       prsi                             , &
       prsl                             , &
       t2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t     (ibMax,kMax)  !(ibMax,kMax)! temperature (K)
       t3          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t     (ibMax,kMax)  !(ibMax,kMax)! temperature (K)
       q2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv    (ibMax,kMax)  !(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       q3          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv    (ibMax,kMax)  !(ibMax,kMax,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       ql2         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)  !(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       ql3         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)  !(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       qi2         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)  !(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       qi3         (1:iMax,1:kMax,latco), &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)  !(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
       w2          (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_omega (ibMax,kMax)  !(ibMax,kMax)! vertical pressure velocity (Pa/s) 
       icefrac     (1:iMax)             , &! REAL(r8), INTENT(in)  :: icefrac     (ibMax)       !sea ice fraction (fraction)
       landfrac    (1:iMax)             , &! REAL(r8), INTENT(in)  :: landfrac    (ibMax)       !land fraction (fraction)
       ocnfrac     (1:iMax)             , &! REAL(r8), INTENT(in)  :: ocnfrac     (ibMax)       !ocean fraction (fraction)
       landm       (1:iMax)             , &! REAL(r8), INTENT(in)  :: landm       (ibMax)       !land fraction ramped over water
       snowh       (1:iMax)             , &! REAL(r8), INTENT(in)  :: snowh       (ibMax)       !Snow depth over land, water equivalent (m)
       dlf         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_dlf   (ibMax,kMax)  !detrained water from ZM
       rliq        (1:iMax)             , &! REAL(r8), INTENT(in)  :: rliq        (ibMax)       !vertical integral of liquid not yet in q(ixcldliq)
       cmfmc       (1:iMax,1:kMax+1)    , &! REAL(r8), INTENT(in)  :: state_cmfmc (ibMax,kMax+1)!convective mass flux--m sub c
       cmfmc2      (1:iMax,1:kMax+1)    , &! REAL(r8), INTENT(in)  :: state_cmfmc2(ibMax,kMax+1)!shallow convective mass flux--m sub c
       concld      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_concld(ibMax,kMax)  !convective cloud cover
       cld         (1:iMax,1:kMax)      , &! REAL(r8), INTENT(out) :: state_cld   (ibMax,kMax)  !cloud fraction
       sst         (1:iMax)             , &! REAL(r8), INTENT(in)  :: sst         (ibMax)       !sea surface temperature
       !zdu         (1:iMax,1:kMax)     , &! REAL(r8), INTENT(in)  :: state_zdu  (ibMax,kMax)  !detrainment rate from deep convection
       prec_str    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_str   (ibMax)       ![Total] sfc flux of precip from stratiform (m/s) 
       snow_str    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_str   (ibMax)       ![Total] sfc flux of snow from stratiform   (m/s)
       prec_sed    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_sed   (ibMax)       !surface flux of total cloud water from sedimentation
       snow_sed    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_sed   (ibMax)       !surface flux of cloud ice from sedimentation
       prec_pcw    (1:iMax)             , &! REAL(r8), INTENT(out)  :: prec_pcw   (ibMax)       !sfc flux of precip from microphysics(m/s)
       snow_pcw    (1:iMax)             , &! REAL(r8), INTENT(out)  :: snow_pcw   (ibMax)       !sfc flux of snow from microphysics (m/s)
       dtdt        (1:iMax,1:kMax)      , &
       dqdt        (1:iMax,1:kMax)      , &
       dqldt       (1:iMax,1:kMax)      , &
       dqidt       (1:iMax,1:kMax)        )

       Total_Rain=Total_Rain+MAX((prec_sed*0.5_r8*dt),0.0_r8)+MAX((prec_pcw*0.5_r8*dt),0.0_r8)
       Total_Snow=Total_Snow+MAX((snow_sed*0.5_r8*dt),0.0_r8)+MAX((snow_pcw*0.5_r8*dt),0.0_r8)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'HGFS') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
!----------------------------
       DO k=1,kMax
          DO i=1,iMax
             press(i,k)=prsl(i,k)
             tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
          END DO
       END DO    
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       r1000=1.0e0_r8 /gasr
       !rbyg=gasr/grav*del(1)*0.5e0_r8
       DO i=1,iMax
          rbyg=gasr/grav*DeltaP(i,1)*0.5e0_r8
          RHO     (i,1)=r1000*(prsl(i,1))/t2(i,1)
          delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
       END DO

       DO k=2,kMax
          DO i=1,iMax
             RHO (i,k)=r1000*(prsl(i,k))/t2(i,k)
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                       LOG(press(i,k-1)/press(i,k))/grav
          END DO
       END DO
  !-----------------------------------------------------------------------
       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_FERRIER (&
       iMax                            , &!INTEGER      , INTENT(IN)     :: nCols
       kMax                            , &!INTEGER      , INTENT(IN)     :: kMax
       DT                              , &!REAL(KIND=r8), INTENT(IN)     :: DT
       RHO        (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(IN),    :: DEL       (ims:ime, kms:kme, jms:jme)
       t3         (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: gt (ims:ime, kms:kme, jms:jme)
       q3         (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: qv (ims:ime, kms:kme, jms:jme)
       ql3        (1:iMax,1:kMax,latco), &!REAL(KIND=r8), INTENT(INOUT), :: qc (ims:ime, kms:kme, jms:jme)
       qi3        (1:iMax,1:kMax,latco), &!REAL(KIND=r8), INTENT(INOUT), :: qi (ims:ime, kms:kme, jms:jme)
       qr3        (1:iMax,1:kMax,latco), &!REAL(KIND=r8), INTENT(INOUT), :: qr (ims:ime, kms:kme, jms:jme)
       F_ICE_PHY  (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)            , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)
       F_RAIN_PHY (1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)       , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)
       F_RIMEF_PHY(1:iMax,1:kMax)      , &!REAL(KIND=r8), INTENT(INOUT), :: F_RIMEF_PHY(ims:ime, kms:kme, jms:jme)
       RAINNCV    (1:iMax)             , &!REAL(KIND=r8), INTENT(INOUT), :: RAINNCV    (ims:ime,  jms:jme) 
       SNOWNCV    (1:iMax)             , &!REAL(KIND=r8), INTENT(INOUT), :: SNOWNCV    (ims:ime,  jms:jme) 
       prsi                            , &
       prsl                            , &
       colrad     (1:iMax)             , &
       dtdt       (1:iMax,1:kMax)      , &
       dqdt       (1:iMax,1:kMax)      , &
       dqldt      (1:iMax,1:kMax)      , &
       dqidt      (1:iMax,1:kMax)      , &
       dqrdt      (1:iMax,1:kMax)         )

!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
      ELSE IF(TRIM(ILCON).EQ.'HWRF') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
!----------------------------
       DO k=1,kMax
          DO i=1,iMax
             press(i,k)=prsl(i,k)
             tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
          END DO
       END DO    
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       r1000=1.0e0_r8 /gasr
       !rbyg=gasr/grav*del(1)*0.5e0_r8
       DO i=1,iMax
          rbyg=gasr/grav*DeltaP(i,1)*0.5e0_r8
          RHO     (i,1)=r1000*(prsl(i,1))/t2(i,1)
          delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
       END DO

       DO k=2,kMax
          DO i=1,iMax
             RHO (i,k)=r1000*(prsl(i,k))/t2(i,k)
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                       LOG(press(i,k-1)/press(i,k))/grav
          END DO
       END DO
  !-----------------------------------------------------------------------
       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8

  CALL RunMicro_thompson( &
       iMax                             , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                             , &!INTEGER      , INTENT(IN   ) :: kMax 
       prsi       (1:iMax,1:kMax+1)     , &
       prsl       (1:iMax,1:kMax)       , &
       t3         (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       q3         (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       ql3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       qr3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       qi3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       qs3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       qg3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NR3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       nifa3      (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       nwfa3      (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       nc3        (1:iMax,1:kMax,latco) , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)   
       dtdt       (1:iMax,1:kMax)        , &!
       dqdt       (1:iMax,1:kMax)        , &!
       dqldt      (1:iMax,1:kMax)        , &!
       dqrdt      (1:iMax,1:kMax)        , &!
       dqidt      (1:iMax,1:kMax)        , &!
       dqsdt      (1:iMax,1:kMax)        , &!
       dqgdt      (1:iMax,1:kMax)        , &!
       dnidt      (1:iMax,1:kMax)        , &!
       dnrdt      (1:iMax,1:kMax)        , &!
       dnifadt    (1:iMax,1:kMax)        , &!
       dnwfadt    (1:iMax,1:kMax)        , &!
       dncdt      (1:iMax,1:kMax)        , &!
       TKE        (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(IN   ) :: TKE (1:nCols, 1:kMax)   
       PBL_CoefKh (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax)   
       DT                               , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       w2         (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       EFFCS      (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       EFFIS      (1:iMax,1:kMax)       , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       RAINNCV    (1:iMax)              , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       SNOWNCV    (1:iMax)                )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)
       Total_Rain=Total_Rain + RAINNCV
       Total_Snow=Total_Snow + SNOWNCV
       snow_str=0.0_r8
       prec_sed=0.0_r8
       snow_sed=0.0_r8
       prec_pcw=0.0_r8
       snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1


!--------------------------------------------

      ELSE IF(TRIM(ILCON).EQ.'UKMO') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5000.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
!----------------------------
       DO k=1,kMax
          DO i=1,iMax
             press(i,k)=prsl(i,k)
             tv   (i,k)=t2 (i,k)*(1.0_r8+0.608_r8*q2(i,k))
          END DO
       END DO    
       !
       !  Calculate the distance between the surface and the first layer of the model
       !
       r1000=1.0e0_r8 /gasr
       !rbyg=gasr/grav*del(1)*0.5e0_r8
       DO i=1,iMax
          rbyg=gasr/grav*DeltaP(i,1)*0.5e0_r8
          RHO     (i,1)=r1000*(prsl(i,1))/t2(i,1)
          delz    (i,1)=MAX((rbyg * tv(i,1)),0.5_r8)*0.75_r8
       END DO

       DO k=2,kMax
          DO i=1,iMax
             RHO (i,k)=r1000*(prsl(i,k))/t2(i,k)
             delz(i,k)=0.5_r8*gasr*(tv(i,k-1)+tv(i,k))* &
                       LOG(press(i,k-1)/press(i,k))/grav
          END DO
       END DO
  !-----------------------------------------------------------------------

       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_UKME(&
                        kMax                            , &
                        iMax                            , &
                        prsi       (1:iMax,1:kMax+1)    , &
                        prsl       (1:iMax,1:kMax)      , &
                        dt                              , &
                        pblh       (1:iMax      )       , &
                        colrad     (1:iMax)             , &
                        kuo        (1:iMax)             , &
                        q3         (1:iMax,1:kMax)      , &!q              , &
                        qi3        (1:iMax,1:kMax,latco), &!QCF    , &
                        ql3        (1:iMax,1:kMax,latco), &!QCL    , &
                        qr3        (1:iMax,1:kMax,latco), &!qcf2   , &
                        t3         (1:iMax,1:kMax)      , &!T  , &
                        F_RIMEF_PHY(1:iMax,1:kMax      ), &! CF     , &
                        F_RAIN_PHY (1:iMax,1:kMax      ), &! CFL    , &
                        F_ICE_PHY  (1:iMax,1:kMax      ), &! CFF    , &
                        dtdt       (1:iMax,1:kMax      ), &
                        dqdt       (1:iMax,1:kMax      ), &
                        dqldt      (1:iMax,1:kMax      ), &
                        dqidt      (1:iMax,1:kMax      ), &
                        dqrdt      (1:iMax,1:kMax      ), &
                        snowh      (1:iMax)             , &
                        landfrac   (1:iMax)             , &
                        terr       (1:iMax)             , &
                        RAINNCV    (1:iMax)             , &
                        SNOWNCV    (1:iMax)               )


!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'MORR') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5000.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
  !-----------------------------------------------------------------------

       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_MORR( &
       iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax 
       prsi        (1:iMax,1:kMax+1)        , &
       prsl        (1:iMax,1:kMax)          , &
       t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       dtdt        (1:iMax,1:kMax)          , &
       dqdt        (1:iMax,1:kMax)          , &
       dqldt       (1:iMax,1:kMax)          , &
       dqrdt       (1:iMax,1:kMax)          , &
       dqidt       (1:iMax,1:kMax)          , &
       dqsdt       (1:iMax,1:kMax)          , &
       dqgdt       (1:iMax,1:kMax)          , &
       dnidt       (1:iMax,1:kMax)          , &
       dnsdt       (1:iMax,1:kMax)          , &
       dnrdt       (1:iMax,1:kMax)          , &
       dNGdt       (1:iMax,1:kMax)          , &
       DT                                   , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       w2         (1:iMax,1:kMax)           , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       RAINNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       SNOWNCV    (1:iMax)                    )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)

!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

      ELSE IF(TRIM(ILCON).EQ.'HUMO') THEN
       ! grell mask
       snowh   =0.0_r8
       icefrac =0.0_r8
       landfrac=0.0_r8
       ocnfrac =0.0_r8
       DO i=1,iMax
          IF(mask(i) == 13 .or. mask(i) == 15 )snowh(i)=5000.0_r8      
          IF(mask(i).GT.0_i8)THEN
             ! land
             icefrac (i)=0.0_r8
             landfrac(i)=1.0_r8
             ocnfrac (i)=0.0_r8
          ELSE
             ! water/ocean
             landfrac(i)=0.0_r8
             ocnfrac(i) =1.0_r8
             IF(ocnfrac(i).GT.0.01_r8.AND.tsfc(i).LT.260.0_r8) THEN
                icefrac(i)=1.0_r8
                ocnfrac(i) =0.0_r8
             ENDIF
          END IF
       END DO
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
  !-----------------------------------------------------------------------

       RAINNCV=0.0_r8
       SNOWNCV=0.0_r8
       CALL RunMicro_HugMorr( &
       iMax                                 , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                                 , &!INTEGER      , INTENT(IN   ) :: kMax 
       prsi        (1:iMax,1:kMax+1)        , &
       prsl        (1:iMax,1:kMax)          , &
       t3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: Tc (1:nCols, 1:kMax)
       q3          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(INOUT) :: qv (1:nCols, 1:kMax)
       ql3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qc (1:nCols, 1:kMax)
       qr3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qr (1:nCols, 1:kMax)
       qi3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qi (1:nCols, 1:kMax)
       QS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qs (1:nCols, 1:kMax)
       QG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: qg (1:nCols, 1:kMax)
       NI3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ni (1:nCols, 1:kMax)
       NS3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: ns (1:nCols, 1:kMax)
       NR3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: nr (1:nCols, 1:kMax)
       NG3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NG (1:nCols, 1:kMax)   
       NC3         (1:iMax,1:kMax,latco)    , &!REAL(KIND=r8), INTENT(INOUT) :: NC (1:nCols, 1:kMax)   
       dtdt        (1:iMax,1:kMax)          , &!
       dqdt        (1:iMax,1:kMax)          , &!
       dqldt       (1:iMax,1:kMax)          , &!
       dqrdt       (1:iMax,1:kMax)          , &!
       dqidt       (1:iMax,1:kMax)          , &!
       dqsdt       (1:iMax,1:kMax)          , &!
       dqgdt       (1:iMax,1:kMax)          , &!
       dnidt       (1:iMax,1:kMax)          , &!
       dnsdt       (1:iMax,1:kMax)          , &!
       dnrdt       (1:iMax,1:kMax)          , &!
       dNGdt       (1:iMax,1:kMax)          , &!
       dNCdt       (1:iMax,1:kMax)          , &!
       TKE         (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: TKE (1:nCols, 1:kMax) (m^2 s-2) 
       PBL_CoefKh  (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: KZH (1:nCols, 1:kMax) (M^2 S-1)  
       DT                                   , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       w2          (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(IN   ) :: omega  ! omega (Pa/s)
       EFFCS       (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFCS (1:nCols, 1:kMax)   ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       EFFIS       (1:iMax,1:kMax)          , &!REAL(KIND=r8), INTENT(OUT  ) :: EFFIS (1:nCols, 1:kMax)   ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
       RAINNCV    (1:iMax)                  , &!REAL(KIND=r8), INTENT(OUT) :: LSRAIN(1:nCols)
       SNOWNCV    (1:iMax)                    )!REAL(KIND=r8), INTENT(OUT) :: LSSNOW(1:nCols)

!--------------------------------------------
        Total_Rain=Total_Rain + RAINNCV
        Total_Snow=Total_Snow + SNOWNCV
        snow_str=0.0_r8
        prec_sed=0.0_r8
        snow_sed=0.0_r8
        prec_pcw=0.0_r8
        snow_pcw=0.0_r8
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

    END IF

!    DO i=1,iMax
!       PRINT*,TRIM(ILCON),'   ',RAINNCV(i),SNOWNCV(i)
!    END DO

    ! Save humd/temp after large scale convection
    DO k=1,kMax
      DO i=1,iMax
        ! PRINT*,t3(i,k),dtdt(i,k),q3(i,k),dqdt(i,k)

         q3  (i,k)       = MAX(q3  (i,k)      ,qmin)
         IF (microphys) THEN
            ql3 (i,k,latco) = MAX(ql3 (i,k,latco),qmin)
            qi3 (i,k,latco) = MAX(qi3 (i,k,latco),qmin)

            ql2 (i,k,latco) = MAX(ql2 (i,k,latco),qmin)
            qi2 (i,k,latco) = MAX(qi2 (i,k,latco),qmin)
            IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
               qr3 (i,k,latco) = MAX(qr3 (i,k,latco),1.0e-12_r8)
               IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF')THEN
                  qs3 (i,k,latco) = MAX(qs3 (i,k,latco),qmin)
                  qg3 (i,k,latco) = MAX(qg3 (i,k,latco),qmin)
                  NI3 (i,k,latco) = MAX(NI3 (i,k,latco),qmin)
                  IF( TRIM(ILCON).EQ.'HWRF')THEN
                      nifa3 (i,k,latco) = MAX(nifa3 (i,k,latco),qmin)
                  ELSE
                      NS3   (i,k,latco) = MAX(NS3 (i,k,latco),qmin)
                  END IF
                  NR3 (i,k,latco) = MAX(NR3 (i,k,latco),qmin)
                  IF( TRIM(ILCON).EQ.'HWRF')THEN
                     nwfa3 (i,k,latco) = MAX(nwfa3 (i,k,latco),qmin)
                  ELSE
                     NG3 (i,k,latco) = MAX(NG3 (i,k,latco),qmin)
                  END IF

                  IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF')NC3 (i,k,latco) = MAX(NC3 (i,k,latco),0.0e-12_r8)

               END IF
               qr2 (i,k,latco) = MAX(qr2 (i,k,latco),1.0e-12_r8)
               IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF')THEN
                  qs2 (i,k,latco) = MAX(qs2 (i,k,latco),qmin)
                  qg2 (i,k,latco) = MAX(qg2 (i,k,latco),qmin)
                  NI2 (i,k,latco) = MAX(NI2 (i,k,latco),qmin)
                  IF( TRIM(ILCON).EQ.'HWRF')THEN
                     nifa2 (i,k,latco) = MAX(nifa2 (i,k,latco),qmin)
                  ELSE
                     NS2 (i,k,latco) = MAX(NS2 (i,k,latco),qmin)
                  END IF
                  NR2 (i,k,latco) = MAX(NR2 (i,k,latco),qmin)
                  IF( TRIM(ILCON).EQ.'HWRF')THEN
                     nwfa2 (i,k,latco) = MAX(nwfa2 (i,k,latco),qmin)
                  ELSE
                     NG2 (i,k,latco) = MAX(NG2 (i,k,latco),qmin)
                  END IF
                  IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF')NC2 (i,k,latco) = MAX(NC2 (i,k,latco),qmin)
               END IF
            END IF
         END IF
         tLrgs(i,k)=t3(i,k)
         qLrgs(i,k)=q3(i,k)
      END DO
    END DO

    IF (microphys) THEN
       DO k=1,kMax
          DO i=1,iMax
              gliqp(i,k) =ql3     (i,k,latco) 
              gliqt(i,k) = gliqt(i,k) + dqldt    (i,k) 
      
              gicep(i,k) =qi3     (i,k,latco) 
              gicet(i,k) = gicet(i,k)+ dqidt    (i,k) 
 
              gliqm(i,k)   =ql2   (i,k,latco) 
              gicem(i,k)   =qi2   (i,k,latco) 
              IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
                 gvarp(i,k,1) =qr3   (i,k,latco) 
                 gvart(i,k,1) =  gvart(i,k,1)+dqrdt(i,k)
                 IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF')THEN
                    gvarp(i,k,2) = qs3  (i,k,latco) 
                    gvart(i,k,2) = gvart(i,k,2) + dqsdt(i,k)
                    gvarp(i,k,3) = qg3  (i,k,latco) 
                    gvart(i,k,3) =  gvart(i,k,3) + dqgdt(i,k)
                    gvarp(i,k,4) = NI3  (i,k,latco) 
                    gvart(i,k,4) =  gvart(i,k,4)  + dnidt(i,k)
                    IF( TRIM(ILCON).EQ.'HWRF')THEN 
                       gvarp(i,k,5) = nifa3  (i,k,latco) 
                       gvart(i,k,5) = gvart(i,k,5) + dnifadt(i,k)
                    ELSE
                       gvarp(i,k,5) = NS3  (i,k,latco) 
                       gvart(i,k,5) =  gvart(i,k,5) +dnsdt(i,k)
                    END IF
                    gvarp(i,k,6) = NR3   (i,k,latco) 
                    gvart(i,k,6) =  gvart(i,k,6)+dnrdt(i,k)
                    IF( TRIM(ILCON).EQ.'HWRF')THEN 
                       gvarp(i,k,7) =nwfa3 (i,k,latco) 
                       gvart(i,k,7) =  gvart(i,k,7) +dnwfadt(i,k)
                    ELSE
                       gvarp(i,k,7) =NG3   (i,k,latco) 
                       gvart(i,k,7) = gvart(i,k,7) + dngdt(i,k)
                    END IF
                    IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF') THEN
                       gvarp(i,k,8) =NC3   (i,k,latco) 
                       gvart(i,k,8) = gvart(i,k,8)+dncdt(i,k)
                    END IF
                 END IF
                 gvarm(i,k,1) =qr2   (i,k,latco) 
                 dqrdt(i,k)  =0.0_r8
                 IF(TRIM(ILCON).EQ.'MORR'.or.TRIM(ILCON).EQ.'HUMO'.OR.  TRIM(ILCON).EQ.'HWRF')THEN
                    gvarm(i,k,2) =qs2   (i,k,latco) 
                    dqsdt(i,k)  =0.0_r8
                    gvarm(i,k,3) =qg2   (i,k,latco) 
                    dqgdt(i,k)  =0.0_r8
                    gvarm(i,k,4) =NI2   (i,k,latco) 
                    dnidt(i,k)  =0.0_r8
                    IF( TRIM(ILCON).EQ.'HWRF')THEN
                       gvarm(i,k,5) =nifa2   (i,k,latco) 
                       dnifadt(i,k)  =0.0_r8
                    ELSE
                       gvarm(i,k,5) =NS2   (i,k,latco) 
                       dnsdt(i,k)  =0.0_r8
                    END IF
                    gvarm(i,k,6) =NR2   (i,k,latco) 
                    dnrdt(i,k)  =0.0_r8
                    IF( TRIM(ILCON).EQ.'HWRF')THEN
                       gvarm(i,k,7) =nwfa2   (i,k,latco) 
                       dnwfadt(i,k)  =0.0_r8
                    ELSE
                       gvarm(i,k,7) =NG2   (i,k,latco) 
                       dngdt(i,k)  =0.0_r8
                    END IF

                    IF(TRIM(ILCON).EQ.'HUMO'.or.TRIM(ILCON).EQ.'HWRF')THEN
                       gvarm(i,k,8) =NC2   (i,k,latco) 
                       dncdt(i,k)  =0.0_r8
                    END IF
                 END IF 
              END IF
          END DO
       END DO   
    END IF


  END SUBROUTINE RunMicroPhysics

!**********************************************************************************
!**********************************************************************************

  SUBROUTINE  qnegat2 (fq, fdq, rdt, prsi, iMax, kMax)
    !
    ! input: fq  specific humidity (dimensionless mixing ratio)
    !        fp  surface pressure (cb)
    ! ouput: fq  adjusted specific humidity
    !        fp  unchanged
    !        fdq distribution of moisture modification
    !
    ! iMax......Number of grid points on a gaussian latitude circle   
    ! kMax......Number of sigma levels  
    ! imx.......=iMax+1 or iMax+2   :this dimension instead of iMax
    !              is used in order to avoid bank conflict of memory
    !              access in fft computation and make it efficient. the
    !              choice of 1 or 2 depends on the number of banks and
    !              the declared type of grid variable (real*4,real*8)
    !              to be fourier transformed.
    !              cyber machine has the symptom.
    !              cray machine has no bank conflict, but the argument
    !              'imx' in subr. fft991 cannot be replaced by iMax    
    ! del.......sigma spacing for each layer computed in routine "setsig".  
    ! dfact.....del(k+1)/del(k)
    !
    INTEGER, INTENT(in   ) :: iMax  
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8)   , INTENT(in   ) :: rdt

    REAL(KIND=r8),    INTENT(inout) :: fq   (iMax,kMax)
    REAL(KIND=r8),    INTENT(inout) :: fdq  (iMax,kMax)  
    REAL(KIND=r8),    INTENT(in   ) :: prsi (iMax,kMax+1)

    REAL(KIND=r8)   :: dfact(iMax,kMax)
    REAL(KIND=r8)   :: DeltaP(iMax,kMax)  
    REAL(KIND=r8)   :: rdt2
    INTEGER :: klev
    INTEGER :: kblw
    INTEGER :: i
    INTEGER :: k  
    DO k=1,kMax
      DO i=1,iMax
          DeltaP(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO

    rdt2=rdt
    DO k=1,kMax-1
      DO i=1,iMax
         dfact(i,k+1) = DeltaP(i,k+1)/DeltaP(i,k)
      END DO
    END DO
    !     
    !     ecmwf vertical borrowing scheme
    !     fdq contains compensated borrowing above first level, uncompensated
    !     borrowing in first level
    !     
    DO k=1,kMax-1
       klev = kMax-k+1
       kblw = klev - 1
       DO i=1,iMax
          fdq(i,klev) = fq(i,klev)
          IF(fq(i,klev).LT.0.0e0_r8) fq(i,klev) = 1.0e-12_r8
          fdq(i,klev) = fq(i,klev) - fdq(i,klev)
          fq(i,kblw) = fq(i,kblw) - fdq(i,klev)*dfact(i,klev)
       END DO
    END DO

    DO i=1,iMax
       fdq(i,1) = fq(i,1)
       IF(fq(i,1).LT.0.0e0_r8) fq(i,1) = 1.0e-12_r8
       fdq(i,1) = fq(i,1) - fdq(i,1)
    END DO

  END SUBROUTINE qnegat2

!**********************************************************************************
!**********************************************************************************

  SUBROUTINE FinalizeMicroPhysics()
    IMPLICIT NONE
    DEALLOCATE(ql3)
    DEALLOCATE(qi3)
    DEALLOCATE(qr3)
    DEALLOCATE(qs3)
    DEALLOCATE(qg3)
    DEALLOCATE(NI3)
    DEALLOCATE(NS3)
    DEALLOCATE(nifa3)
    DEALLOCATE(NR3)
    DEALLOCATE(NG3)
    DEALLOCATE(nwfa3)
    DEALLOCATE(NC3)
    
    DEALLOCATE(ql2)
    DEALLOCATE(qi2)
    DEALLOCATE(qr2)
    DEALLOCATE(qs2)
    DEALLOCATE(qg2)
    DEALLOCATE(NI2)
    DEALLOCATE(NS2)
    DEALLOCATE(nifa2)
    DEALLOCATE(NR2)
    DEALLOCATE(NG2)
    DEALLOCATE(nwfa2)
    DEALLOCATE(NC2)

  END SUBROUTINE FinalizeMicroPhysics

END MODULE MicroPhysics
