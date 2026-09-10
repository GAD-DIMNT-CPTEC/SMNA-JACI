MODULE DeepConvection

  USE Constants, ONLY :  &
       delq             ,&
       r8,i8,qmin,grav,gasr

  USE Cu_Kuolcl, ONLY :       &
      Init_Cu_Kuolcl ,RunCu_Kuolcl   
      
  USE Cu_ZhangMcFarlane, ONLY:       &
      Init_Cu_ZhangMcFarlane ,RunCu_ZhangMcFarlane

  USE Cu_Grellens, ONLY: &
      Init_Cu_GrellEns,RunCu_GrellEns

  USE Cu_Grellens_CPTEC, ONLY:  &
      Init_Cu_GrellEnsCPTEC,RunCu_GrellEnsCPTEC
  
  USE Cu_GDM_BAM, Only:Init_Cu_GDM_BAM,RunCu_GDM_BAM

  USE ModConRas, ONLY: &
      Init_Cu_RelAraSch,RunCu_RelAraSch

  USe Cu_RAS3PHASE, Only:    &
      Init_Cu_RAS3PHASE,Run_Cu_RAS3PHASE,Finalize_Cu_RAS3PHASE
  
   IMPLICIT NONE
  SAVE

 PRIVATE


 PUBLIC :: InitDeepConvection
 PUBLIC :: RunDeepConvection
 PUBLIC :: FinalizeDeepConvection

CONTAINS

 SUBROUTINE InitDeepConvection( &
                          dt      ,a_hybr, b_hybr, trunc, &
                          kMax    ,iMax,jMax,ibMax,jbMax  , &
                          fhour   ,idate,iccon        )
  IMPLICIT NONE
    REAL(KIND=r8)   , INTENT(IN   ) :: dt
    INTEGER         , INTENT(IN   ) :: trunc
    INTEGER         , INTENT(IN   ) :: kMax
    INTEGER         , INTENT(IN   ) :: iMax
    INTEGER         , INTENT(IN   ) :: jMax
    INTEGER         , INTENT(IN   ) :: ibMax
    INTEGER         , INTENT(IN   ) :: jbMax
    REAL(KIND=r8)   , INTENT(IN   ) :: a_hybr (kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: b_hybr (kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: fhour
    INTEGER         , INTENT(IN   ) :: idate(4)
    CHARACTER(LEN=*), INTENT(IN   ) :: iccon

    CHARACTER(LEN=*), PARAMETER     :: h='**(InitDeepConvection)**'

    IF(TRIM(iccon).EQ.'ARA')  CALL Init_Cu_RelAraSch   (trunc)
    IF(TRIM(iccon).EQ.'KUO')  CALL Init_Cu_Kuolcl () 
    IF(TRIM(iccon).EQ.'GRE')  CALL Init_Cu_GrellEns  ()
    IF(TRIM(iccon).EQ.'GEC')  CALL Init_Cu_GrellEnsCPTEC()
    IF(TRIM(iccon).EQ.'GDN')  CALL Init_Cu_GDM_BAM()
    IF(TRIM(iccon).EQ.'ZMC')  CALL Init_Cu_ZhangMcFarlane(dt,a_hybr,b_hybr,ibMax,kMax,jbMax,jMax)

    IF(TRIM(iccon).EQ.'RAS')  CALL Init_Cu_RAS3PHASE(kMax,2*dt,fhour,idate,iMax,jMax,ibMax,jbMax)

 END SUBROUTINE InitDeepConvection
!**********************************************************************************
!**********************************************************************************


 SUBROUTINE RunDeepConvection( &
      ! Run Flags
            iccon     ,& !CHARACTER(LEN=*),    INTENT(IN   ) :: iccon
            cflric    ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: cflric
            microphys   , & !LOGICAL      , INTENT(in   ) :: microphys
      ! Time info
            tod       ,& 
            dt        ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: dt
      ! Model Geometry
      ! Model information
            jdt        ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            nClass    ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            nAeros    ,&
            iMax      ,& !INTEGER         ,    INTENT(in   ) :: iMax  
            kMax      ,& !INTEGER         ,    INTENT(in   ) :: kMax
            nls       ,& !INTEGER         ,    INTENT(in   ) :: nls 
            latco     ,& !INTEGER         ,    INTENT(in   ) :: latco  
            mask      ,& !INTEGER(KIND=i8),    INTENT(IN   ) :: mask      (iMax)
            terr      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: terr      (iMax)  
            colrad    ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: colrad    (iMax)
            lonrad    ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: lonrad    (iMax)  
      ! Surface field
            sens      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: sens      (iMax)
            tpert     ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tpert     (iMax)
            ustar     ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tpert     (iMax)
      ! PBL field
            pblh      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: pblh      (iMax)
            tke       ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: tke       (iMax,kMax)
      ! CONVECTION: Cloud field
            ktop      ,& !INTEGER         ,    INTENT(inout) :: ktop      (iMax)
            kbot      ,& !INTEGER         ,    INTENT(inout) :: kbot      (iMax)
            ktops     ,& !INTEGER         ,    INTENT(INOUT) :: ktops     (iMax)
            kbots     ,& !INTEGER         ,    INTENT(INOUT) :: kbots     (iMax)
            cape_old  ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cape_old  (iMax)
            cape      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cape      (iMax)
            cine      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: cine      (iMax)
            hrem      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: hrem      (iMax,kmax)
            qrem      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: qrem      (iMax,kmax)
            kuo       ,& !INTEGER         ,    INTENT(inout) :: kuo       (iMax)
            cldm      ,& !REAL(KIND=r8)   ,    INTENT(inout) :: cldm      (iMax)
            plcl      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: plcl      (iMax)
            ddql      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: ddql      (iMax,kMax)
            ddmu      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: ddmu      (iMax,kMax)
            fdqn      ,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: fdqn      (iMax,kMax)
            qliq      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: qliq      (iMax,kMax) !qrc liquid water content in
            deep_newcld,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: cld       (iMax,kMax)
            cmfmc     ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: cmfmc     (iMax,kMax+1)
            dlf       ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: dlf       (iMax,kMax)
            zdu       ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: zdu       (iMax,kMax)
            rliq      ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: rliq      (iMax)
            RAINCV    ,& !REAL(KIND=r8)   ,    INTENT(inout) :: RAINCV    (iMax)
            SNOWCV    ,& !REAL(KIND=r8)   ,    INTENT(inout) :: SNOWCV    (iMax)
            snow_zmc  ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: snow_zmc  (iMax)
            prec_zmc  ,& !REAL(KIND=r8)   ,    INTENT(OUT  ) :: prec_zmc  (iMax)
            Total_Rain,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Rain(iMax)
            Total_Snow,& !REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Snow(iMax)
            dudt      ,&
            dvdt      ,&
            dtdt      ,&
            dqdt      ,&
            dqldt     ,&
            dqidt     ,&
      ! Atmospheric fields
            tc        ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: tc        (iMax,kmax)
            qc        ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: qc        (iMax,kmax)
            phii      ,& 
            phil      ,& 
            prsi      ,& !
            prsl      ,& !
            ps2       ,& !REAL(KIND=r8)   ,    INTENT(in   ) :: ps2       (iMax)
            ub        ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: ub        (iMax,kMax) ! (m/s) 
            vb        ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: vb        (iMax,kMax) ! (m/s)
            omgb      ,& !REAL(KIND=r8)   ,    INTENT(IN   ) :: omgb      (iMax,kMax) ! (Pa/s)
            T2        ,& !REAL(KIND=r8)   ,    INTENT(inout) :: T2        (iMax,kmax)
            T3        ,& !REAL(KIND=r8)   ,    INTENT(inout) :: T3        (iMax,kmax)
            Q1        ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q1        (iMax,kmax)
            Q2        ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q2        (iMax,kmax)
            Q3        ,& !REAL(KIND=r8)   ,    INTENT(inout) :: Q3        (iMax,kmax)
            dq        ,& !REAL(KIND=r8)   ,    INTENT(inout) :: dq        (iMax,kmax)
      ! Microphysics
            ql2       , & !REAL(KIND=r8)  , INTENT(inout) :: ql2    (1:iMax,1:kmax)
            qi2       , & !REAL(KIND=r8)  , INTENT(inout) :: qi2     (1:iMax,1:kmax)
            ql3       , & !REAL(KIND=r8)   ,    INTENT(inout) :: ql3       (iMax,kmax)
            qi3       , & !REAL(KIND=r8)   ,    INTENT(inout) :: qi3       (iMax,kmax)
            daerdt    , &
            gvarm     , &
            gvarp       )
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
    !    iccon                  ! cu schemes ex. KUO, ARA, GRE ..
    !    kuo                     ! convection yes(1) or not(0) for shallow convection
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

    LOGICAL         ,    INTENT(in   ) :: microphys    !LOGICAL      , INTENT(in   ) :: microphys

    CHARACTER(LEN=*),    INTENT(IN   ) :: iccon
    REAL(KIND=r8)   ,    INTENT(in   ) :: tod
    REAL(KIND=r8)   ,    INTENT(in   ) :: dt
    REAL(KIND=r8)   ,    INTENT(in   ) :: cflric
    INTEGER         ,    INTENT(in   ) :: jdt                         !INTEGER         ,    INTENT(in   ) :: iMax  
    INTEGER         ,    INTENT(in   ) :: nClass                         !INTEGER         ,    INTENT(in   ) :: iMax  
    INTEGER         ,    INTENT(IN   ) :: nAeros
    INTEGER         ,    INTENT(in   ) :: iMax  
    INTEGER         ,    INTENT(in   ) :: kMax
    INTEGER         ,    INTENT(in   ) :: nls 
    INTEGER         ,    INTENT(in   ) :: latco  
    REAL(KIND=r8)   ,    INTENT(in   ) :: tc      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(in   ) :: qc      (iMax,kmax)
    INTEGER(KIND=i8),    INTENT(IN   ) :: mask    (iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: terr    (iMax)  
    REAL(KIND=r8)   ,    INTENT(IN   ) :: colrad  (iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: lonrad  (iMax)  
    REAL(KIND=r8)   ,    INTENT(IN   ) :: cape_old(1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: cape    (1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: cine    (1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: sens    (1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: pblh    (1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: tpert   (1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: ustar   (1:iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: phii     (1:iMax,1:kMax+1)  
    REAL(KIND=r8)   ,    INTENT(IN   ) :: phil     (1:iMax,1:kMax  )  
    REAL(KIND=r8)   ,    INTENT(IN   ) :: prsi    (1:iMax,1:kMax+1)   !
    REAL(KIND=r8)   ,    INTENT(IN   ) :: prsl    (1:iMax,1:kMax  )   !
    REAL(KIND=r8)   ,    INTENT(in   ) :: ps2     (iMax)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: ub      (iMax,kMax) ! (m/s) 
    REAL(KIND=r8)   ,    INTENT(IN   ) :: vb      (iMax,kMax) ! (m/s)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: omgb    (iMax,kMax) ! (Pa/s)
    REAL(KIND=r8)   ,    INTENT(IN   ) :: tke     (iMax,kMax)
    INTEGER         ,    INTENT(inout) :: ktop    (iMax)
    INTEGER         ,    INTENT(inout) :: kbot    (iMax)
    REAL(KIND=r8)   ,    INTENT(inout) :: RAINCV  (iMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: SNOWCV  (iMax)
    REAL(KIND=r8)   ,    INTENT(inout) :: hrem    (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: qrem    (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: T2      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: T3      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: Q1      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: Q2      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: Q3      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: dq      (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: ql2     (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: qi2     (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: ql3     (iMax,kmax)
    REAL(KIND=r8)   ,    INTENT(inout) :: qi3     (iMax,kmax)

    INTEGER         ,    INTENT(inout) :: kuo     (iMax)
    REAL(KIND=r8)   ,    INTENT(inout) :: cldm    (iMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: plcl    (iMax)
    INTEGER         ,    INTENT(INOUT) :: ktops   (iMax)
    INTEGER         ,    INTENT(INOUT) :: kbots   (iMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: ddql    (1:iMax,1:kMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: ddmu    (1:iMax,1:kMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: fdqn    (1:iMax,1:kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: qliq    (iMax,kMax) !qrc liquid water content in cloud after rainout
    REAL(KIND=r8)   ,    INTENT(INOUT) :: deep_newcld(1:iMax,1:kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: cmfmc   (1:iMax,1:kMax+1)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dlf     (1:iMax,1:kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: zdu     (1:iMax,1:kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: rliq    (1:iMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: snow_zmc(1:iMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: prec_zmc(1:iMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Rain(iMax)
    REAL(KIND=r8)   ,    INTENT(INOUT) :: Total_Snow(iMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dudt(iMax,kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dvdt(iMax,kMax)
    REAL(KINd=r8)   ,    INTENT(OUT  ) :: dtdt (iMax,kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dqdt (iMax,kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dqldt(iMax,kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dqidt(iMax,kMax)

    REAL(KIND=r8)   , OPTIONAL,   INTENT(OUT  ) :: daerdt(iMax,kMax,nClass+nAeros) 
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)
    REAL(KIND=r8)   , OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)

    ! WIND COMPONENTS FOR GRELL ENSEMBLE
    REAL(KIND=r8) :: u2(iMax,kMax)
    REAL(KIND=r8) :: v2(iMax,kMax)
    REAL(KIND=r8) :: w2(iMax,kMax)
    REAL(KIND=r8) :: rain1(iMax)
    INTEGER       :: mask2(iMax)
    INTEGER       :: i
    INTEGER       :: k,ntrac
    qliq    =0.0_r8
    deep_newcld     =0.0_r8
    cmfmc   =0.0_r8
    dlf     =0.0_r8
    zdu     =0.0_r8
    rliq    =0.0_r8
    snow_zmc=0.0_r8
    prec_zmc=0.0_r8
    u2=0.0_r8
    v2=0.0_r8
    w2=0.0_r8
    dudt=0.0_r8
    dvdt=0.0_r8
    dtdt =0.0_r8
    dqdt =0.0_r8
    dqldt=0.0_r8
    dqidt=0.0_r8
    IF((nClass+nAeros)>0 .and. PRESENT(daerdt))THEN
    daerdt=0.0_r8
    END IF
    IF(TRIM(iccon).EQ.'ARA') THEN
       ! grell mask
       DO i=1,iMax      
          IF(mask(i).GT.0_i8)THEN
             mask2(i)=0 ! land
          ELSE
             mask2(i)=1 ! water/ocean
          END IF
       END DO
       ! grell wind
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       ntrac=0

       CALL RunCu_RelAraSch(iMax, KMax,tod,dt,mask2,terr,&
       t2,t3,  q2,ql2,ql3,qi2,qi3, q3,    u2,    v2, w2, prsi,prsl,phii,phil, &
       dudt,dvdt,kbot,  ktop,  kuo,raincv,dtdt ,dqdt ,dqldt,dqidt )
       
       DO i=1,iMax
          plcl(i) = prsi(i,kbot(i))/1000.0_r8 !convert Pa to Cb
       END DO
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

    ! Copy deep convection rain to total raini

       Total_Rain=RAINCV
       Total_Snow=SNOWCV

    END IF

    IF(TRIM(iccon).EQ.'KUO')THEN
       CALL  RunCu_Kuolcl(dt,prsi,prsl, q1, q3, &
            T3, dq, RAINCV, kuo, plcl, ktop, &
            kbot, iMax, kMax,dtdt,dqdt)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
    ! Copy deep convection rain to total raini
       Total_Rain=RAINCV
       Total_Snow=SNOWCV
    END IF
!******************

    IF(TRIM(iccon).EQ.'RAS')THEN

       ! grell mask
       DO i=1,iMax      
          IF(mask(i).GT.0_i8)THEN
             mask2(i)=0 ! land
          ELSE
             mask2(i)=1 ! water/ocean
          END IF
       END DO

       ! grell wind
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       ntrac=0
       
       IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
          CALL Run_Cu_RAS3PHASE(microphys,nClass,nAeros,iMax,latco,ntrac, kMax,dt,jdt,terr,&
                             t2,t3,q2,q3,ql2,ql3,qi2,qi3,u2,v2,prsi,prsl,phii,phil,&
                             dtdt ,dqdt ,dqldt,dqidt ,dudt,dvdt,colrad,ustar,pblh,tke,mask2,&
                             kbot ,ktop ,kuo  ,RAINCV,gvarm,gvarp,daerdt)
       ELSE
          CALL Run_Cu_RAS3PHASE(microphys,nClass,nAeros,iMax,latco, ntrac,kMax,dt,jdt,terr,&
                             t2,t3,q2,q3,ql2,ql3,qi2,qi3,u2,v2,prsi,prsl,phii,phil,&
                             dtdt ,dqdt ,dqldt,dqidt ,dudt,dvdt,colrad,ustar,pblh,tke,mask2,&
                             kbot ,ktop ,kuo  ,RAINCV)
       END IF
       !CALL shllcl(dt    ,ps2   ,sl   ,si   ,q3    ,T3    ,kuo   , &
       !     plcl  ,ktops ,kbots ,iMax ,kMax  )
       DO i=1,iMax
          plcl(i) = prsl(i,kbot(i))/1000.0_r8 !convert Pa to cb
       END DO
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

    ! Copy deep convection rain to total raini

       Total_Rain=RAINCV
       Total_Snow=SNOWCV
    END IF

!******************
    IF(TRIM(iccon).EQ.'GRE')THEN

       ! grell mask
       DO i=1,iMax      
          IF(mask(i).GT.0_i8)THEN
             mask2(i)=0 ! land
          ELSE
             mask2(i)=1 ! water/ocean
          END IF
       END DO

       ! grell wind
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       
       CALL RunCu_GrellEns(&
              prsi  ,prsl  ,u2    ,v2    ,w2    ,t2    ,t3    ,q2    ,&
              q3    ,terr  ,mask2 ,dt    ,RAINCV,kuo   ,ktop  ,kbot  ,plcl  ,&
              qliq  ,iMax  ,kMax  ,tke   ,latco ,lonrad,colrad,ddql  ,ddmu  ,&
              cape  ,cine  ,sens  ,cape_old,dtdt,dqdt)

       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

    ! Copy deep convection rain to total raini

       Total_Rain=RAINCV
       Total_Snow=SNOWCV
    END IF

    IF(TRIM(iccon).EQ.'GEC')THEN

       ! grell mask
       DO i=1,iMax      
          IF(mask(i).GT.0_i8)THEN
             mask2(i)=0 ! land
          ELSE
             mask2(i)=1 ! water/ocean
          END IF
       END DO

       ! grell wind
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       
       CALL RunCu_GrellEnsCPTEC(&
       prsi   ,prsl   ,u2     ,v2     ,w2     ,T2     ,T3     , &
       q2     ,q3     ,terr   ,mask2  ,dt     ,RAINCV , &
       kuo    ,ktop   ,kbot   ,plcl   ,qliq   ,iMax   ,kMax,dtdt ,dqdt )

       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

    ! Copy deep convection rain to total raini

       Total_Rain=RAINCV
       Total_Snow=SNOWCV
    END IF

!******************
    IF(TRIM(iccon).EQ.'GDN')THEN

       ! grell mask
       DO i=1,iMax      
          IF(mask(i).GT.0_i8)THEN
             mask2(i)=0 ! land
          ELSE
             mask2(i)=1 ! water/ocean
          END IF
       END DO

       ! grell wind
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       
       CALL RunCu_GDM_BAM(&
              prsi  ,prsl  ,u2    ,v2    ,w2    ,t2    ,t3    ,q2    , &
              q3    ,ql3   ,terr  ,mask2 ,dt    ,jdt   ,RAINCV,kuo   ,ktop  , &
              kbot  ,plcl  ,qliq  ,iMax  ,kMax  ,tke   ,latco ,lonrad,colrad, &
              ddql  ,ddmu  ,deep_newcld,cape_old,cape  ,cine  ,sens ,dtdt  ,dqdt  ,dqldt)

       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1

    ! Copy deep convection rain to total raini

       Total_Rain=RAINCV
       Total_Snow=SNOWCV
    END IF


    IF(TRIM(iccon).EQ.'ZMC')THEN
       DO k=1,kMax
          DO i=1,iMax
             u2(i,k)=ub(i,k)
             v2(i,k)=vb(i,k)
             w2(i,k)=omgb(i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       CALL RunCu_ZhangMcFarlane( &
       iMax                         , &! INTEGER, INTENT(in) :: pcols   ! number of columns (max)
       kMax+1                       , &! INTEGER, INTENT(in) :: pverp   ! number of vertical levels + 1
       kMax                         , &! INTEGER, INTENT(in) :: pver    ! number of vertical levels
       latco                        , &! INTEGER, INTENT(in) :: latco    ! number of latitudes
       1                            , &! INTEGER, INTENT(in) :: pcnst=1 ! number of advected constituents (including water vapor)
       2                            , &! INTEGER, INTENT(in) :: pnats=2 ! number of non-advected constituents
       prsi    (1:iMax,1:kMax+1)    , &! 
       prsl    (1:iMax,1:kMax)      , &! 
       terr    (1:iMax)*grav        , &! REAL(r8), INTENT(in)  :: state_phis   (pcols)    !(pcols)  ! surface geopotential
       t3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t   (pcols,pver)!(pcols,pver)! temperature (K)
       q3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv  (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
       ql3     (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_ql (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
       qi3     (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qi (pcols,pver)!(pcols,pver,ppcnst)! ice  mixing ratio (kg/kg moist or dry air depending on type)
       w2      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
       deep_newcld(1:iMax,1:kMax)   , &! REAL(r8), INTENT(out) :: state_cld(ibMax,kMax)          !  cloud fraction
       prec_zmc(1:iMax)             , &! real(r8), intent(out) :: prec(pcols)   ! total precipitation(m/s)
       pblh    (1:iMax)             , &! real(r8), intent(in)  :: pblh(pcols)
       cmfmc   (1:iMax,1:kMax+1)    , &! real(r8), intent(out) :: cmfmc(pcols,pverp)
       tpert   (1:iMax)             , &! real(r8), intent(in)  :: tpert(pcols)
       dlf     (1:iMax,1:kMax)      , &! real(r8), intent(out) :: dlf(pcols,pver)! scattrd version of the detraining cld h2o tend
       zdu     (1:iMax,1:kMax)      , &! real(r8), intent(out) :: zdu(pcols,pver)
       rliq    (1:iMax)             , &! real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
       2.0_r8*dt                    , &! real(r8), intent(in)  :: ztodt       ! 2 delta t (model time increment)
       snow_zmc(1:iMax)             , &! real(r8), intent(out) :: snow(pcols)   ! snow from ZM convection 
       ktop    (1:iMax)             , &
       kbot    (1:iMax)             , &
       kuo     (1:iMax)             , &
       dtdt    (1:iMax,1:kMax)      , &
       dqdt    (1:iMax,1:kMax)      , &
       dqldt   (1:iMax,1:kMax)      , &
       dqidt   (1:iMax,1:kMax)        )

       RAINCV=MAX(prec_zmc*dt,0.0_r8)
       SNOWCV=MAX(snow_zmc*dt,0.0_r8)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
       Total_Rain=RAINCV
       Total_Snow=SNOWCV
    END IF


 END SUBROUTINE RunDeepConvection

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

 SUBROUTINE FinalizeDeepConvection()
  IMPLICIT NONE

 END SUBROUTINE FinalizeDeepConvection

END MODULE DeepConvection
