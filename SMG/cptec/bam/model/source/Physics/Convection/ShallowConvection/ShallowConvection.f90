MODULE ShallowConvection

  USE Constants, ONLY :  &
       delq             ,&
       r8,i8,qmin,grav,gasr

  USE Shall_Tied, ONLY:       &
      InitShall_Tied        , &
      shalv2
      
   USE Shall_JHack, ONLY:       &
      convect_shallow_init ,convect_shallow_tend

  USE Shall_Souza, ONLY: shallsouza
 
  USE Shall_MasFlux, ONLY: Init_Shall_MasFlux,run_shall_masflux


   IMPLICIT NONE
  SAVE

 PRIVATE



 PUBLIC  :: InitShallowConvection
 PUBLIC  :: RunShallowConvection
 PUBLIC  :: FinalizeShallowConvection

CONTAINS

 SUBROUTINE InitShallowConvection(a_hybr, b_hybr, &
                                   kMax    ,jMax,ibMax,jbMax,trunc, &
                                   ppcnst   ,ISCON        )
   IMPLICIT NONE
    INTEGER         , INTENT(IN   ) :: kMax
    INTEGER         , INTENT(IN   ) :: jMax
    INTEGER         , INTENT(IN   ) :: ibMax
    INTEGER         , INTENT(IN   ) :: jbMax
    INTEGER         , INTENT(IN   ) :: trunc
    INTEGER         , INTENT(IN   ) :: ppcnst
    REAL(KIND=r8)   , INTENT(IN   ) :: a_hybr (kMax+1)
    REAL(KIND=r8)   , INTENT(IN   ) :: b_hybr (kMax+1)
    CHARACTER(LEN=*), INTENT(IN   ) :: ISCON
    CHARACTER(LEN=*), PARAMETER     :: h='**(InitShallowConvection)**'


    IF(TRIM(ISCON).EQ.'TIED') CALL InitShall_Tied()
    IF(TRIM(ISCON).EQ.'JHK'.or. TRIM(ISCON).EQ.'UW')  CALL convect_shallow_init(ISCON,kMax,jMax,ibMax,jbMax,ppcnst,a_hybr,b_hybr)

    IF(TRIM(ISCON).EQ.'MFLX') CALL Init_Shall_MasFlux (kMax,trunc)

 END SUBROUTINE InitShallowConvection

 SUBROUTINE RunShallowConvection( &
      ! Run Flags
                      ISCON     , & !CHARACTER(LEN=*), INTENT(IN   ) :: ISCON
                      iccon     , & !CHARACTER(LEN=*), INTENT(IN   ) :: iccon
      ! Time info
                      tod                            ,& 
                      jdt       , & !INTEGER         , INTENT(in   ) :: jdt
                      dt        , & !REAL(KIND=r8)   , INTENT(in   ) :: dt
      ! Model information
                      iMax      , & !INTEGER         , INTENT(in   ) :: iMax
                      kmax      , & !INTEGER         , INTENT(in   ) :: kmax
                      latco     , & !INTEGER         , INTENT(in   ) :: latco
                      terr      , & !REAL(KIND=r8)   , INTENT(in   ) :: terr    (1:iMax)
      ! Model Geometry
      !                si        , & !REAL(KIND=r8)   , INTENT(in   ) :: si      (1:kMax+1)
      !                sl        , & !REAL(KIND=r8)   , INTENT(in   ) :: sl      (1:kmax)
      ! Surface field
                      mask      , & !INTEGER(KIND=i8),INTENT(IN ) :: mask (iMax) 
                      sens      , & !
                      evap      , & !
                      qpert     , & !REAL(KIND=r8)   , INTENT(in   ) :: qpert   (1:iMax)
      ! PBL field
                      pblh      , & !REAL(KIND=r8)   , INTENT(in   ) :: pblh    (1:iMax)
                      tke       , & !REAL(KIND=r8)   , INTENT(in   ) :: tke     (1:iMax,1:kMax)
      ! CONVECTION: Cloud field
                      dudt      , & !
                      dvdt      , & !
                      dtdt      , & !
                      dqdt      , & !
                      dqldt     , & !
                      dqidt     , & !
                      rliq      , & !REAL(KIND=r8)   , INTENT(in   ) :: rliq    (1:iMax)
                      ktop      , & !INTEGER         , INTENT(in   ) :: ktop    (1:iMax)
                      ktops     , & !INTEGER         , INTENT(in   ) :: ktops   (1:iMax)
                      kuo       , & !INTEGER         , INTENT(in   ) :: kuo     (1:iMax)
                      plcl      , & !REAL(KIND=r8)   , INTENT(inout) :: plcl    (1:iMax)
                      kcbot1    , & !INTEGER         , INTENT(inout) :: kcbot1  (1:iMax)
                      kctop1    , & !INTEGER         , INTENT(inout) :: kctop1  (1:iMax)
                      noshal    , & !INTEGER         , INTENT(inout) :: noshal  (1:iMax)
                      concld    , & !REAL(KIND=r8)   , INTENT(inout) :: concld  (1:iMax,1:kMax)
                      cld       , & !REAL(KIND=r8)   , INTENT(inout) :: cld     (1:iMax,1:kMax)
                      cmfmc     , & !REAL(KIND=r8)   , INTENT(inout) :: cmfmc   (1:iMax,1:kMax+1)
                      cmfmc2    , & !REAL(KIND=r8)   , INTENT(out  ) :: cmfmc2  (1:iMax,1:kMax+1)
                      dlf       , & !REAL(KIND=r8)   , INTENT(inout) :: dlf     (1:iMax,1:kMax)
                      fdqn      , & !REAL(KIND=r8)   , INTENT(inout) :: fdqn     (1:iMax,1:kMax)
                      rliq2     , & !REAL(KIND=r8)   , INTENT(out  ) :: rliq2   (1:iMax)
                      snow_cmf  , & !REAL(KIND=r8)   , INTENT(out  ) :: snow_cmf(1:iMax)
                      prec_cmf  , & !REAL(KIND=r8)   , INTENT(out  ) :: prec_cmf(1:iMax)
                      RAINCV    , & !REAL(KIND=r8)   , INTENT(out  ) :: RAINCV  (1:iMax)
                      SNOWCV    , & !REAL(KIND=r8)   , INTENT(out  ) :: SNOWCV  (1:iMax)
      ! Atmospheric fields
                      prsi      ,& !
                      prsl      ,& !
                      phii      ,& !
                      phil      ,& !
                      PS_work   , & !REAL(KIND=r8)   , INTENT(in   ) :: PS_work (1:iMax)
                      ps2       , & !REAL(KIND=r8)   , INTENT(in   ) :: ps2     (1:iMax)
                      ub        , & !REAL(KIND=r8)   , INTENT(in   ) :: ub      (1:iMax,1:kmax)
                      vb        , & !REAL(KIND=r8)   , INTENT(in   ) :: vb      (1:iMax,1:kmax)
                      omgb      , & !REAL(KIND=r8)   , INTENT(in   ) :: omgb    (1:iMax,1:kmax)
                      t3        , & !REAL(KIND=r8)   , INTENT(inout) :: t3      (1:iMax,1:kmax)
                      q3        , & !REAL(KIND=r8)   , INTENT(inout) :: q3      (1:iMax,1:kmax)
                      ql3       , & !REAL(KIND=r8)   , INTENT(inout) :: ql3     (1:iMax,1:kmax)
                      qi3         ) !REAL(KIND=r8)   , INTENT(inout) :: qi3     (1:iMax,1:kmax)

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
    !    iccon                  ! cu schemes ex. KUO, ARA, GRE, RAS, GDN..
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


    CHARACTER(LEN=*), INTENT(IN   ) :: ISCON
    CHARACTER(LEN=*), INTENT(IN   ) :: iccon
    INTEGER         , INTENT(in   ) :: jdt
    INTEGER         , INTENT(in   ) :: iMax
    INTEGER         , INTENT(in   ) :: kmax
    INTEGER         , INTENT(in   ) :: latco
    REAL(KIND=r8)   , INTENT(in   ) :: dt
    REAL(KIND=r8)   , INTENT(in   ) :: tod 
    REAL(KIND=r8)   , INTENT(in   ) :: prsi    (1:iMax,1:kMax+1)  
    REAL(KIND=r8)   , INTENT(in   ) :: prsl    (1:iMax,1:kMax  )  
    REAL(KIND=r8)   , INTENT(in   ) :: phii    (1:iMax,1:kMax+1)  
    REAL(KIND=r8)   , INTENT(in   ) :: phil    (1:iMax,1:kMax  )  
    REAL(KIND=r8)   , INTENT(in   ) :: PS_work (1:iMax)
    REAL(KIND=r8)   , INTENT(in   ) :: ps2     (1:iMax)
    REAL(KIND=r8)   , INTENT(in   ) :: terr    (1:iMax)
    INTEGER(KIND=i8), INTENT(IN   ) :: mask    (1:iMax) 
    REAL(KIND=r8)   , INTENT(in   ) :: sens    (1:iMax)
    REAL(KIND=r8)   , INTENT(in   ) :: evap    (1:iMax)
    REAL(KIND=r8)   , INTENT(in   ) :: qpert   (1:iMax)
    REAL(KIND=r8)   , INTENT(in   ) :: pblh    (1:iMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dudt    (1:iMax,1:kMax)   
    REAL(KIND=r8)   , INTENT(OUT  ) :: dvdt    (1:iMax,1:kMax)   
    REAL(KINd=r8)   , INTENT(OUT  ) :: dtdt    (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqdt    (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqldt   (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqidt   (1:iMax,1:kMax)

    REAL(KIND=r8)   , INTENT(inout) :: rliq    (1:iMax)
    REAL(KIND=r8)   , INTENT(in   ) :: tke     (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(in   ) :: ub      (1:iMax,1:kmax)
    REAL(KIND=r8)   , INTENT(in   ) :: vb      (1:iMax,1:kmax)
    REAL(KIND=r8)   , INTENT(in   ) :: omgb    (1:iMax,1:kmax)
    INTEGER         , INTENT(in   ) :: ktop    (1:iMax)
    INTEGER         , INTENT(in   ) :: ktops   (1:iMax)
    INTEGER         , INTENT(inout) :: kuo     (1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: plcl    (1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: t3      (1:iMax,1:kmax)
    REAL(KIND=r8)   , INTENT(inout) :: q3      (1:iMax,1:kmax)
    REAL(KIND=r8)   , INTENT(inout) :: ql3     (1:iMax,1:kmax)
    REAL(KIND=r8)   , INTENT(inout) :: qi3     (1:iMax,1:kmax)
    INTEGER         , INTENT(inout) :: kcbot1  (1:iMax)
    INTEGER         , INTENT(inout) :: kctop1  (1:iMax)
    INTEGER         , INTENT(inout) :: noshal  (1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: concld  (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(inout) :: cld     (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(inout) :: cmfmc   (1:iMax,1:kMax+1)
    REAL(KIND=r8)   , INTENT(out  ) :: cmfmc2  (1:iMax,1:kMax+1)
    REAL(KIND=r8)   , INTENT(inout) :: dlf     (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(inout) :: fdqn    (1:iMax,1:kMax)
    REAL(KIND=r8)   , INTENT(inout) :: rliq2   (1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: snow_cmf(1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: prec_cmf(1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: RAINCV  (1:iMax)
    REAL(KIND=r8)   , INTENT(inout) :: SNOWCV  (1:iMax)

    ! Wind components for grell ensemble
    REAL(KIND=r8) :: u2(iMax,kMax)
    REAL(KIND=r8) :: v2(iMax,kMax)
    REAL(KIND=r8) :: w2(iMax,kMax)
    LOGICAL       :: newr
    INTEGER       :: i
    INTEGER       :: k
    u2      =0.0_r8
    v2      =0.0_r8
    w2      =0.0_r8
    dudt=0.0_r8
    dvdt=0.0_r8
    dtdt =0.0_r8
    dqdt =0.0_r8
    dqldt=0.0_r8
    dqidt=0.0_r8
    !-----------------------------------------------------------------
    ! Shallow Convection
    !-----------------------------------------------------------------
    IF(TRIM(ISCON).EQ.'TIED')THEN
       IF(TRIM(iccon).EQ.'KUO'.OR.TRIM(iccon).EQ.'RAS'.OR.TRIM(iccon).EQ.'GRE' .OR. &
          TRIM(iccon).EQ.'ZMC'.OR.TRIM(iccon).EQ.'GEC'.OR.TRIM(iccon).EQ.'GDN' .OR. &
          TRIM(iccon).EQ.'ARA')THEN
          newr=.FALSE.
          CALL shalv2( t3, q3, prsi ,prsl , dt, &
               ktop, plcl, kuo, kMax+1, kctop1, kcbot1, noshal, &
               newr, iMax, kMax,dtdt ,dqdt)
       END IF

!       IF(TRIM(iccon).EQ.'ARA') THEN
!          newr=.TRUE.
!          CALL shalv2(si, sl, t3, q3, PS_work, dt, &
!               ktops, plcl, kuo, kMax+1, kctop1, kcbot1, noshal, &
!               newr, iMax, kMax)
!       END IF
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
    END IF
    
    IF(TRIM(ISCON).EQ.'MFLX')THEN
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub   (i,k)
             v2(i,k)=vb   (i,k)
             w2(i,k)=omgb (i,k)*1000.0_r8  ! (Pa/s)
          END DO 
       END DO
       CALL Run_Shall_MasFlux(iMax,kMax,dt,prsi ,prsl ,phii,phil,t3,q3,u2,v2,w2,ql3,qi3,&
                         dudt,dvdt,dtdt,dqdt,dqldt,dqidt,kcbot1,kctop1,kuo,noshal,mask,pblh,sens,evap)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
    END IF
     
    IF(TRIM(ISCON).EQ.'JHK' .OR. TRIM(ISCON).EQ.'UW' )THEN
       !IF (nscalars <= 0) THEN
       !    dlf   =0.0_r8
       !    rliq  =0.0_r8
       !    cmfmc =0.0_r8
       !END IF
       DO i=1,iMax
          DO k=1,kMax
             u2(i,k)=ub   (i,k)
             v2(i,k)=vb   (i,k)
             w2(i,k)=omgb (i,k)*1000.0_r8  ! (Pa/s)
          END DO
       END DO
       call convect_shallow_tend ( &
        latco      , &!INTEGER, INTENT(IN   )  :: latco
        iMax       , &!INTEGER, INTENT(IN   )  :: pcols
        kMax       , &!INTEGER, INTENT(IN   )  :: pver
        kMax+1     , &!INTEGER, INTENT(IN   )  :: pverp
        jdt        , &!INTEGER, INTENT(IN   )  :: nstep
        1          , &! INTEGER, INTENT(in):: pcnst=1 ! number of advected constituents (including water vapor)
        2          , &! INTEGER, INTENT(in):: pnats=2 ! number of non-advected constituents
        2*dt       , &!real(r8), intent(in) :: ztodt        ! 2 delta-t (seconds)
        prsi    (1:iMax,1:kMax+1)    , &
        prsl    (1:iMax,1:kMax)      , &
        terr    (1:iMax)*grav        , &! REAL(r8), INTENT(in)  :: state_phis   (pcols)  !(pcols)   ! surface geopotential
        t3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_t   (pcols,pver)!(pcols,pver)! temperature (K)
        q3      (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qv   (pcols,pver)!(pcols,pver,ppcnst)! vapor  mixing ratio (kg/kg moist or dry air depending on type)
        ql3     (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_ql    (pcols,pver)!(pcols,pver,ppcnst)! liquid mixing ratio (kg/kg moist or dry air depending on type)
        qi3     (1:iMax,1:kMax)      , &! REAL(r8), INTENT(in)  :: state_qi    (pcols,pver)!(pcols,pver,ppcnst)! ice    mixing ratio (kg/kg moist or dry air depending on type)
        u2      (1:iMax,1:kMax)      , &
        v2      (1:iMax,1:kMax)      , &
        w2      (1:iMax,1:kMax)      , &!REAL(r8), INTENT(in   )  :: state_omega  (pcols,pver)!(pcols,pver)! vertical pressure velocity (Pa/s) 
        concld  (1:iMax,1:kMax)      , &!REAL(r8), INTENT(in   )  ::
        cld     (1:iMax,1:kMax)      , &!REAL(r8), INTENT(in   )  ::
        qpert   (1:iMax)             , &!real(r8), intent(in   ) :: qpert(pcols)  ! PBL perturbation specific humidity
        pblh    (1:iMax)             , &!real(r8), intent(in   ) :: pblht(pcols)    ! PBL height (provided by PBL routine)
        cmfmc   (1:iMax,1:kMax+1)    , &!real(r8), intent(inout) :: cmfmc(pcols,pverp)  ! moist deep + shallow convection cloud mass flux
        cmfmc2  (1:iMax,1:kMax+1)    , &!real(r8), intent(out  ) :: cmfmc2(pcols,pverp)  ! moist shallow convection cloud mass flux
        prec_cmf(1:iMax)             , &!real(r8), intent(out  ) :: precc(pcols)     ! convective precipitation rate
        dlf     (1:iMax,1:kMax)      , &!real(r8), intent(inout) :: qc(pcols,pver)      ! dq/dt due to export of cloud water  ! detrained water 
        rliq    (1:iMax)             , &!real(r8), intent(inout) :: rliq(pcols) ! vertical integral of liquid not yet in q(ixcldliq)
        rliq2   (1:iMax)             , &!real(r8), intent(out  ) :: rliq2(pcols) 
        snow_cmf(1:iMax)             , &!real(r8), intent(out  ) :: snow(pcols)  ! snow from this convection     
        tke     (1:iMax,1:kMax)      , &
        ktop    (1:iMax)             , &
        KUO     (1:iMax)             , &
        dtdt    (1:iMax,1:kMax)      , &
        dqdt    (1:iMax,1:kMax)      , &
        dqldt   (1:iMax,1:kMax)      , &
        dqidt   (1:iMax,1:kMax)      , &
        kctop1  (1:iMax)             , &
        kcbot1  (1:iMax)             , &
        noshal  (1:iMax)               )
       prec_cmf=(prec_cmf*dt)
       snow_cmf=(snow_cmf*dt)
       CALL qnegat2 (q3, fdqn,(1.0_r8/dt), prsi, iMax, kMax)! time t+1
       
        RAINCV = RAINCV + MAX(prec_cmf,0.0_r8)
        SNOWCV = SNOWCV + MAX(snow_cmf,0.0_r8)
    END IF
    
    IF(TRIM(ISCON).EQ.'SOUZ')THEN
       WRITE(*,*)"it is not available yet. Problem in the water balance'"
       STOP "ERROR"
       !CALL Shallsouza(t3,q3,PS_work,sl,sens,evap,dt,iMax,kMax,kuo, &
       !     noshal, kctop1, kcbot1, 560.0_r8, 1.6_r8)
    END IF

 END SUBROUTINE RunShallowConvection
!-----------------------------------------------------------------
!-----------------------------------------------------------------
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

!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE FinalizeShallowConvection()
   IMPLICIT NONE
 END SUBROUTINE FinalizeShallowConvection

END MODULE ShallowConvection
