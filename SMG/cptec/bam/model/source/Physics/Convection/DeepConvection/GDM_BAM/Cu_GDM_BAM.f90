!  $Date: 09/2015
!--------------------------------------------------------------------
!
!
!  This convection scheme called Grell_Devenyi Modified Version (GDM) is a modified version of Grell Devenyi (GD) scheme.
!  The modifications are in many things: Closure, it is a single closure scheme; a new entrainment and detrainment
!  have been developed by using a SCM and Cloud Resolving Model (CRM), and modification in downdraft and updraft.
!  The documentation of this version will be  available soon. 
!  This scheme has been developed for the Brazilian Atmosphere Global Model (BAM).    
!  NOTE: It is not official version yet!! It can be used only for tests
!  more informations: nilo.figueroa@cptec.inpe.br
!  Reference: Silvio Nilo Figueroa, Paulo Kubota, Enver Ramirez and G. Grell (2015). Paper in preparation
!---------------------------------------------------------------------------------------------------------------------

MODULE Cu_GDM_BAM
  !
  ! GDM--|---GDM2-|
  !                         |cup_env
  !                         |  
  !                         |cup_env
  !                         |  
  !                         |cup_env_clev
  !                         |  
  !                         |cup_env_clev
  !                         |  
  !                         |cup_maximi
  !                         |  
  !                         |cup_kbcon
  !                         |  
  !                         |cup_minimi
  !                         |  
  !                         |cup_up_he
  !                         |  
  !                         |cup_up_he
  !                         |  
  !                         |cup_ktop
  !                         |  
  !                         |cup_minimi
  !                         |  
  !                         |cup_up_nms
  !                         |  
  !                         |cup_up_nms
  !                         |  
  !                         |cup_dd_nms
  !                         |  
  !                         |cup_dd_nms
  !                         |  
  !                         |cup_dd_he
  !                         |  
  !                         |cup_dd_he
  !                         |  
  !                         |cup_dd_moisture
  !                         |  
  !                         |cup_dd_moisture
  !                         |  
  !                         |cup_up_moisture
  !                         |  
  !                         |cup_up_moisture
  !                         |  
  !                         |cup_up_aa0
  !                         |  
  !                         |cup_up_aa0
  !                         |  
  !                         |cup_dd_edt
  !                         |  
  !                         |cup_dellabot
  !                         |  
  !                         |cup_dellabot
  !                         |  
  !                         |cup_dellas
  !                         |  
  !                         |cup_dellas
  !                         |  
  !                         |cup_env
  !                         |  
  !                         |cup_env_clev
  !                         |  
  !                         |cup_up_he
  !                         |  
  !                         |cup_up_nms
  !                         |  
  !                         |cup_dd_nms
  !                         |  
  !                         |cup_dd_he
  !                         |  
  !                         |cup_dd_moisture
  !                         |  
  !                         |cup_up_moisture
  !                         |  
  !                         |cup_up_aa0
  !                         |  
  !                         |cup_maximi
  !                         |  
  !                         |cup_kbcon
  !                         |  
  !                         |cup_forcing_ens_16
  !                         |  
  !                         |cup_output_ens
  ! 
  ! USE Diagnostics, ONLY: updia, dodia , &
  !        nDiag_deepdcape,&
  !        nDiag_deepwf, &
  !        nDiag_deepmu,  &
  !        nDiag_deepmd,  &
  !        nDiag_deepde, &
  !        nDiag_deepen, &
  !        nDiag_deepql ,  &
  !        nDiag_deepqc ,  &
  !        nDiag_deephe ,  &
  !        nDiag_deephes, &
  !        nDiag_deephc ,  &
  !        nDiag_deeprh 

  USE Constants, ONLY :  &
       grav          , &
       undef           , &
       r8
  !USE Options, ONLY :       &
  !     gdmpar1           , &!      
  !     gdmpar2           

  ! USE wv_saturation,Only : ,findsp_mask

  USE PhysicalFunctions,ONLY : fpvs2es5

  ! USE Parallelism, ONLY: &
  !      myId,&
  !      MsgOne, &
  !      FatalError

  IMPLICIT NONE
SAVE

  PRIVATE
  PUBLIC :: Init_Cu_GDM_BAM
  PUBLIC :: RunCu_GDM_BAM
  REAL(KIND=r8)   , PARAMETER :: cp   =1004.0_r8
  REAL(KIND=r8)   , PARAMETER :: xl   =2.5e06_r8
  REAL(KIND=r8)   , PARAMETER :: rv   =461.9_r8
  INTEGER, PARAMETER :: iens     =  1 
  INTEGER, PARAMETER :: iens_tmp =  1
  INTEGER, PARAMETER :: mjx      =  1 
  INTEGER, PARAMETER :: maxens   =  3 ! ensemble one on mbdt from PARAME
  INTEGER, PARAMETER :: maxens2  =  3 ! ensemble two on precip efficiency
  INTEGER, PARAMETER :: maxens3  = 1 ! 6! ensemble in cup_forcing  
  INTEGER, PARAMETER :: ensdim   = 1*maxens*maxens2*maxens3 !9 54 
  REAL(KIND=r8)   , PARAMETER :: tcrit    =   273.15_r8
  !
  ! workfunctions for downdraft
  !
  REAL(KIND=r8)            :: ae  (2)
  REAL(KIND=r8)            :: be  (2)
  REAL(KIND=r8)            :: ht  (2)

  REAL(KIND=r8)               :: radius
  REAL(KIND=r8)               :: entr_rate
  REAL(KIND=r8)               :: mentrd_rate
  REAL(KIND=r8)               :: mentr_rate
  REAL(KIND=r8)               :: edtmin
  REAL(KIND=r8)               :: edtmax
  REAL(KIND=r8)               :: edtmax1
  REAL(KIND=r8)               :: effmax
  REAL(KIND=r8)               :: depth_min
  REAL(KIND=r8)               :: cap_maxs 
  REAL(KIND=r8)               :: cap_maxs_land 
  REAL(KIND=r8)               :: cap_max_increment
  REAL(KIND=r8)               :: zkbmax
  REAL(KIND=r8)               :: zcutdown
  REAL(KIND=r8)               :: z_detr


CONTAINS
  SUBROUTINE Init_Cu_GDM_BAM()
    !
    ! specify entrainmentrate and detrainmentrate
    ! Larger radius will give less mass fluix and make cloud grow taller
    ! and shift heating. Recomend 12 km.
    ! snf        radius=12000.
    !
    radius=12000.0_r8
    !  gross entrainment rate
    entr_rate=0.2_r8/radius
    ! entrainment of mass
    !
    mentrd_rate=entr_rate
    mentr_rate =entr_rate
    !
    ! initial detrainmentrates
    !
    !
    ! max/min allowed value for epsilon (ratio downdraft base 
    ! mass flux/updraft base mass flux
    !
    edtmin=0.2_r8
    !
    ! snf   edtmax=0.60_r8
    ! snf   edtmax=0.75_r8
    !
    edtmax=0.99_r8    ! ok over land snf  
    edtmax1=0.99_r8   ! ok over ocean
    !---------------
    effmax=0.99_r8
    edtmax=effmax
    edtmax1=effmax 
    !
    !  minimum depth (m), clouds must have
    !
    depth_min=500.0_r8
    !
    ! maximum depth (mb) of capping
    !     if(iens == 3)cap_max2=50.0_r8
    !     if(iens == 2)cap_max2=75.0_r8 
    !     if(iens == 1)cap_max2=100.0_r8
    ! original    cap_maxs=125.0_r8 cap_max_increment=50.0_r8      !new
    !
    cap_maxs=120.0_r8  !!
    cap_maxs_land=120.0_r8  !!
    !--------------------------------
    !!    cap_max_increment=50.0_r8
    !!snf modified   
    cap_max_increment=30.0_r8 

    !
    ! max height(m) above ground where updraft air can originate
    !
    zkbmax=4000.0_r8
    !
    ! height(m) above which no downdrafts are allowed to originate
    !
    zcutdown=4000.0_r8 !!!!nilo3 3000.0_r8
    !
    ! depth(m) over which downdraft detrains all its mass
    !
    z_detr=1250.0_r8
    !
    ht(1)=xl/cp

    ht(2)=2.834e6_r8/cp

    be(1)=0.622_r8*ht(1)/0.286_r8

    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)

    be(2)=0.622_r8*ht(2)/0.286_r8

    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)

  END SUBROUTINE Init_Cu_GDM_BAM
  SUBROUTINE RunCu_GDM_BAM(&
       prsi    ,prsl   ,u2     ,v2     ,w2    ,t2     ,t3    ,q2     ,&
       q3      ,ql3    ,terr   ,xland  ,dtime  ,jdt   ,RAINCV ,kuo   ,ktop   ,&
       kbot    ,plcl   ,qliq   ,nCols  ,kMax   ,tke   ,latco  ,lon   ,lat    ,&
       ddql    ,ddmu   ,deep_newcld ,cape2   ,cape3   ,cine3  ,sens,dtdt  ,dqdt  ,dqldt )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime
    INTEGER, INTENT(IN   ) :: jdt
    REAL(KIND=r8)   , INTENT(IN   ) :: prsi  (1:nCols,1:kMax+1)   !
    REAL(KIND=r8)   , INTENT(IN   ) :: prsl  (1:nCols,1:kMax  )   !
    REAL(KIND=r8)   , INTENT(IN   ) :: u2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: v2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: w2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: t3    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: q3    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: ql3   (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: q2    (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: terr  (1:nCols       )
    INTEGER         , INTENT(IN   ) :: xland (1:nCols       )
    INTEGER         , INTENT(INOUT) :: kuo   (1:nCols       )
    INTEGER         , INTENT(INOUT) :: ktop  (1:nCols       )
    INTEGER         , INTENT(INOUT) :: kbot  (1:nCols       )
    REAL(KIND=r8)   , INTENT(OUT  ) :: qliq  (1:nCols,1:kMax)

    !
    ! output variables after cumulus parameterization
    !
    REAL(KIND=r8)   , INTENT(INOUT) :: RAINCV(1:nCols       )
    REAL(KIND=r8)   , INTENT(INOUT) :: plcl  (1:nCols       ) 
    !---------nilo
    INTEGER         , INTENT(IN   ) :: latco
    REAL(KIND=r8)   , INTENT(IN   ) :: lon(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: lat(nCols)
    REAL(KIND=r8)                   :: lat2(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: tke   (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: ddql (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: ddmu (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: deep_newcld (1:nCols,1:kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cape2(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cape3(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cine3(1:nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: sens (1:nCols)
    REAL(KINd=r8)   , INTENT(OUT  ) :: dtdt  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqdt  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dqldt (nCols,kMax)

    !------------------
    !
    ! local variables
    !
    ! at time t
    REAL(KIND=r8)                   :: t        (nCols,kMax)
    REAL(KIND=r8)                   :: q        (nCols,kMax)
    REAL(KIND=r8)                   :: p        (nCols,kMax)
    REAL(KIND=r8)                   :: us       (nCols,kMax)
    REAL(KIND=r8)                   :: vs       (nCols,kMax)
    ! at time t+1
    REAL(KIND=r8)                   :: tn       (nCols,kMax)
    REAL(KIND=r8)                   :: qn       (nCols,kMax)
    REAL(KIND=r8)                   :: pn       (nCols,kMax)
    REAL(KIND=r8)                   :: omeg     (nCols,kMax)
    REAL(KIND=r8)                   :: psur     (nCols     )
    ! tendencies
    REAL(KIND=r8)                   :: outt     (nCols,kMax) 
    REAL(KIND=r8)                   :: outq     (nCols,kMax) 
    REAL(KIND=r8)                   :: outqc    (nCols,kMax)
    ! Precipitation
    REAL(KIND=r8)                   :: pre1     (nCols     )

    ! Auxiliars

    INTEGER                         :: i,k, kk
    INTEGER                         :: ierr     (nCols)
    REAL(KIND=r8)                   :: bncy     (nCols,kMax)
    REAL(KIND=r8)                   :: massfln  (nCols,ensdim)
    REAL(KIND=r8)                   :: pw(nCols)
    REAL(KIND=r8)     :: dp1, dp2,aa,bb,cc,dd,ee

    !new 2011-dec-updated Sep2015
    REAL(KIND=r8)     :: xmb      (nCols     )
    REAL(KIND=r8)     :: out1_mass_flux_up(nCols,kMax)
    REAL(KIND=r8)     :: out2_mass_flux_d(nCols,kMax)
    REAL(KIND=r8)     :: out3_detr(nCols,kMax)
    REAL(KIND=r8)     :: out4_liq_water(nCols,kMax)
    REAL(KIND=r8)     :: out5_qc(nCols,kMax)
    REAL(KIND=r8)     :: out6_wf(nCols)
    REAL(KIND=r8)     :: out7_entr(nCols,kMax)
    REAL(KIND=r8)     :: out8_he(nCols,kMax)
    REAL(KIND=r8)     :: out9_hes(nCols,kMax)
    REAL(KIND=r8)     :: out10_hc(nCols,kMax)
    REAL(KIND=r8)     :: out11_rh(nCols,kMax)
    !-------------
    REAL(KIND=r8)     :: dcape(nCols)
    REAL(KIND=r8)     :: ddwf(nCols)
    !REAL(KIND=r8)     :: ddmu(nCols,kMax) !OUT
    REAL(KIND=r8)     :: ddmd(nCols,kMax)  
    REAL(KIND=r8)     :: ddde(nCols,kMax)
    REAL(KIND=r8)     :: dden(nCols,kMax)
    !REAL(KIND=r8)     :: ddql(nCols,kMax) !OUT
    REAL(KIND=r8)     :: ddqc(nCols,kMax)
    REAL(KIND=r8)     :: ddhe(nCols,kMax)
    REAL(KIND=r8)     :: ddhes(nCols,kMax)
    REAL(KIND=r8)     :: ddhc(nCols,kMax)
    REAL(KIND=r8)     :: ddrh(nCols,kMax)
    REAL(KIND=r8)     :: DeltaP(nCols,kMax)
    !---------------
    REAL(KIND=r8)     :: qlmedia(nCols)
    !------CLOUD
    INTEGER           :: k_newddcl
    aa=0.0_r8
    !
    !    prepare input, erase output
    !
    DO i=1,nCols
       kuo (i)=0
       kbot(i)=1
       ktop(i)=1
       ierr(i)=0
       plcl(i) = 1.0e0_r8
       lat2(i)=90.0_r8-lat(i)*180.0_r8/3.1415926_r8
    END DO
    DO k=1,kMax
      DO i=1,nCols
          DeltaP(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO

    !
    !move  variables from GCM to local variables
    !
    DO k = 1, kMax
       DO i = 1,nCols
          dtdt  (i,k)=0.0_r8
          dqdt  (i,k)=0.0_r8
          dqldt (i,k)=0.0_r8
          !p   (i,k) = ps2(i)*sl(k)*10.0_r8               ! pressure in mbar
          p   (i,k) = prsl(i,k)/100.0_r8
          !pn  (i,k) = ps2(i)*sl(k)*10.0_r8
          pn  (i,k) = prsl(i,k)/100.0_r8
          us  (i,k) = u2(i,k) 
          vs  (i,k) = v2(i,k) 
          omeg(i,k) = w2(i,k)   !test 
          !
          !
          q   (i,k) = q2(i,k)
          qn  (i,k) = q3(i,k) 
          !----------------------------------------------------------q to r
          ! q(i,k)=q2(i,k)/(1.0_r8-q2(i,k))
          ! qn(i,k)=q3(i,k)/(1.0_r8-q3(i,k))
          !-------------------------------------------------------------
          t   (i,k) = t2(i,k)
          tn  (i,k)=  t3(i,k)
          IF(TN(I,K) < 200.0_r8)    TN(I,K) = T(I,K)
          IF(QN(I,K) < 1.E-08_r8)  QN(I,K) = 1.E-08_r8
          psur (i) = prsi(i,1) /100.0_r8!ps2(1:nCols)*10.0_r8
       END DO
    END DO
    !
    ! call cumulus parameterization
    !
    pre1   = 0.0_r8
    outt   = 0.0_r8
    outq   = 0.0_r8
    outqc  = 0.0_r8
    xmb    =0.0_r8               !total base mass flux
    !-----
    out1_mass_flux_up=0.0_r8
    out2_mass_flux_d=0.0_r8
    out3_detr=0.0_r8
    out4_liq_water=0.0_r8
    out5_qc=0.0_r8
    out6_wf=0.0_r8
    out7_entr=0.0_r8
    out8_he=0.0_r8
    out9_hes=0.0_r8
    out10_hc=0.0_r8
    out11_rh=0.0_r8
    ddql=0.0_r8
    ddhc=0.0_r8
    ddmu=0.0_r8
    ddmd=0.0_r8
    dden=undef
    ddde=undef 
    ddwf=0.0_r8
    ddhe=0.0_r8
    ddhes=0.0_r8
    ddhc=0.0_r8
    ddrh=0.0_r8
    !---------------------------------------------
    ! calculate precipitable water and DCAPE
    pw=0.0_r8    
    dcape  =0.0_r8
    DO k=1,kmax
       DO i = 1,ncols
          pw(i) = pw(i) + DeltaP(i,k)*q3(i,k)
       END DO
    END DO
    DO i = 1,ncols
       pw(i)=100.0_r8*pw(i)*psur(i)/9.81_r8
       dd=cape3(i)-cape2(i)
       IF(cape3(i)>0.0_r8.AND.dd>0.0_r8)THEN
          dcape(i)=dd
       ELSE
          dcape(i)=0.0_r8
       ENDIF

    END DO
    !---------------------------
    CALL GDM2(&
         t      ,q      , p   ,us    ,vs     ,tn    ,qn     ,pn      ,psur      , &
         omeg   ,terr    ,dtime,jdt  ,nCols   ,kMax   ,xland  ,pw     ,lat2  ,lon     , &
         tke    ,dcape   ,cine3   ,sens                                            , &
         ierr   ,outt   ,outq   ,outqc  ,pre1   ,ktop    ,kbot    ,xmb           , &
         out1_mass_flux_up,    out2_mass_flux_d,  out3_detr,  out4_liq_water     , &
         out5_qc,    out6_wf ,out7_entr,out8_he, out9_hes, out10_hc, out11_rh       )

    !
    ! after cumulus parameterization
    ! out  tn1, qn1, prec, kuo,ktop, kbot
    !
    DO i = 1,nCols
       IF(ierr(i) == 0)THEN
          RAINCV(i)=dtime*pre1(i)             !in mm/sec(ditme),by 0.5_r8(if leap-frog or 2dt)
          RAINCV(i)=RAINCV(i)/1000.0_r8       !in m for gcm
!!!!!write(*,*)'sss',RAINCV(i)
       END IF
       IF(RAINCV(i) > 0.0_r8)kuo(i)=1
       kk=kbot(i)
       plcl(i)=p(i,kk)/10.0_r8  ! from mb to cb for Shallow convection
    END DO
    !------------------------------------------
    ! Mass FLUX UP (ddmu)=xmb(i)*out1_mass_flux_up(i,k) Fluxo de massa  (onde xmb=fluxo de massa na basae).
    ! Mass Flux DOWN (DDMD)=xmb(i)*out2_mass_flux_d(i,k)  fluxo de massa descendente (xmb= na base)
    ! ddqcc = cloud q (including liquid water) after entrainment
    !-----------------------------------------------


    DO k=1,kMax
       DO i=1,nCols   
          IF(ierr(i) == 0)THEN
             IF(RAINCV(i) > 0.0_r8)THEN
                dtdt  (i,k)=outt(i,k)
                dqdt  (i,k)=outq(i,k)
                dqldt (i,k)=outqc(i,k)

                t3(i,k) = t3(i,k)+ 2.0_r8*outt(i,k)*dtime
                !outq(i,k)=(outq(i,k)+outqc(i,k))/(1.0_r8+outq(i,k)+outqc(i,k))
                q3(i,k) = q3(i,k)+ 2.0_r8*outq(i,k)*dtime
                ql3(i,k)=ql3(i,k) +2.0_r8*outqc(i,k)*dtime
             END IF
          END IF
       END DO
    END DO
    !-------------------------------------------------up
    !2D
    !3D
    DO i=1,nCols
       IF(RAINCV(i) > 0.0_r8)THEN
          ddwf(i)=out6_wf(i) !! aa1, ou aa1-aa0, or xmb 
          DO k=kbot(i),ktop(i)
             !UPDRAFT
             ddmu(i,k)=xmb(i)*out1_mass_flux_up(i,k)
             IF(ddmu(i,k).LT.1.0e-8_r8)ddmu(i,k)=0.0_r8
             ddde(i,k)=MIN(MAX(0.0_r8,out3_detr(i,k)),1.0e-2_r8)
             dden(i,k)=MIN(MAX(0.0_r8,out7_entr(i,k)),1.0e-2_r8)
             ddql(i,k)=out4_liq_water(i,k)
             IF(ddql(i,k).LT.1.0e-8_r8)ddql(i,k)=1.0e-8_r8
             ddqc(i,k)=out5_qc(i,k)
             IF(ddqc(i,k).LT.1.0e-8_r8)ddqc(i,k)=1.0e-8_r8
             ddhe(i,k)=out8_he(i,k)
             ddhes(i,k)=out9_hes(i,k)
             ddhc(i,k)=out10_hc(i,k)
             ddrh(i,k)=out11_rh(i,k)
          ENDDO

          DO k=1,ktop(i)
             ddmd(i,k)=xmb(i)*out2_mass_flux_d(i,k)
          ENDDO
       END IF
    END DO
    !----------------------------
    ! NEW DEEP CLOUD
    !---------------------------- 
    deep_newcld=0.0_r8
    k_newddcl=2
    !-------------------------------------------scaled by Mass flux our Liq water
    IF(k_newddcl==1)THEN
       DO i=1, nCols
          IF(RAINCV(i) > 0.0_r8)THEN
             DO k=kbot(i),ktop(i)
                deep_newcld(i,k) =0.0_r8
                IF(ddrh(i,k)>=1.0_r8)THEN
                   deep_newcld(i,k)=1.0_r8
                ELSE
                   bb=(1000.0_r8*(q3(i,k)/ddrh(i,k)-q3(i,k)))**0.49_r8
                   cc=-100.0_r8*ddql(i,k)/bb
                   deep_newcld(i,k)=(aa**0.25)*(1.0_r8-EXP(cc))
                ENDIF
                deep_newcld(i,k)=MIN (MAX(deep_newcld(i,k), 0.00_r8), 0.6_r8)
             ENDDO
          ENDIF
       ENDDO
    ELSEIF(k_newddcl==2)THEN
       DO i=1, nCols
          IF(RAINCV(i) > 0.0_r8)THEN
             DO k=kbot(i),ktop(i)
                deep_newcld(i,k)=0.0_r8
                !!if(ddmu(i,k+1)>ssmu(i,k+1))then
                !!   bb=ddmu(i,k+1)-ssmu(i,k+1)
                !!   else
                bb=ddmu(i,k+1)
                !!endif
                dp1 = 0.14_r8              !!
                dp2 = 500.0_r8             !! 675.0_r8
                deep_newcld(i,k)=dp1*LOG(1.0_r8+dp2*bb)
                deep_newcld(i,k)=MIN (MAX(deep_newcld(i,k), 0.00_r8), 0.6_r8)
             ENDDO
          ENDIF
       ENDDO
    ENDIF

    !---------------------------------------------------
    !    (RAINCV(i)>=0.0_r8)
    !-------------------------------------------------
    !xmb     total base mass flux
    !------
    DO i=1,nCols
       IF(RAINCV(i)<=0.0_r8)THEN
          !2D
          dcape(i)=undef
          ddwf(i)=undef
          !3D
          ddmu(i,:)=undef
          ddmd(i,:)=undef
          ddde(i,:)=undef
          dden(i,:)=undef
          ddql(i,:)=undef  
          ddqc(i,:)=undef
          ddhe(i,:)=undef
          ddhes(i,:)=undef
          ddhc (i,:)=undef
          ddrh (i,:)=undef
       ENDIF
    ENDDO
    !2D
    !      if (dodia(nDiag_deepdcape))  call updia(dcape              ,nDiag_deepdcape   ,latco,jdt)
    !      if (dodia(nDiag_deepwf)) call updia(ddwf                   ,nDiag_deepwf ,latco,jdt)
    !3D
    !      if (dodia(nDiag_deepmu))  call updia(ddmu                  ,nDiag_deepmu  ,latco,jdt) 
    !      if (dodia(nDiag_deepmd))  call updia(ddmd                  ,nDiag_deepmd  ,latco,jdt)
    !      if (dodia(nDiag_deepde)) call updia(ddde                  ,nDiag_deepde ,latco,jdt)
    !      if (dodia(nDiag_deepen))  call updia(dden                 ,nDiag_deepen   ,latco,jdt)
    !      if (dodia(nDiag_deepql ))  call updia(ddql                  ,nDiag_deepql   ,latco,jdt)
    !      if (dodia(nDiag_deepqc ))  call updia(ddqc                  ,nDiag_deepqc   ,latco,jdt)
    !      if (dodia(nDiag_deephe ))  call updia(ddhe                  ,nDiag_deephe   ,latco,jdt)
    !      if (dodia(nDiag_deephes ))  call updia(ddhes                  ,nDiag_deephes   ,latco,jdt)
    !      if (dodia(nDiag_deephc ))  call updia(ddhc                  ,nDiag_deephc   ,latco,jdt)
    !      if (dodia(nDiag_deeprh ))  call updia(ddrh                  ,nDiag_deeprh   ,latco,jdt)

    RETURN
  END SUBROUTINE RunCu_GDM_BAM


  !*-------------------------

  SUBROUTINE GDM2(& 
       t      ,q      ,p    ,us    ,vs    ,tn     ,qo  ,po  ,psur        , &
       omeg   ,z1     ,dtime,jdt  ,nCols  ,kMax   ,mask   ,pwt    ,lat   ,lon      , &
       tke    ,dcape   ,cine3   ,sens                                            , &
       ierr   , outt   ,outq   ,outqc ,pre    ,ktop  , kbcon  ,xmb             , &
       out1_mass_flux_up,   out2_mass_flux_d, out3_detr,  out4_liq_water       , &
       out5_qc, out6_wf,out7_entr,out8_he, out9_hes,out10_hc,out11_rh       )

    IMPLICIT NONE
    !
    ! input variables
    ! IN
    INTEGER, INTENT(IN   )       :: nCols
    INTEGER, INTENT(IN   )       :: kMax
    INTEGER, INTENT(IN   )       :: jdt
    REAL(KIND=r8)   , INTENT(IN   )       :: dtime
    REAL(KIND=r8)   , INTENT(IN   )       :: pwt    (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: lat (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: lon (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: t    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: q    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: p   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: us   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: vs   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: tn   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: qo   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: po    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: omeg (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )       :: z1   (nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: tke  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT )       :: dcape(nCols)
    REAL(KIND=r8)   , INTENT(IN)       :: cine3(nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: sens(nCols)
    REAL(KIND=r8)   , INTENT(IN   )       :: psur (nCols)
    INTEGER, INTENT(IN   )       :: mask (nCols)
    REAL(KIND=r8)   , INTENT(OUT)       :: xmb(nCols)
    INTEGER, INTENT(INOUT)       :: ierr      (nCols)
    !----
    REAL(KIND=r8)                         :: massfln(nCols,ensdim) !new

    !
    ! output variables
    ! OUT


    REAL(KIND=r8)   , INTENT(INOUT)       :: outt (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: outq (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: outqc (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)       :: pre  (nCols) 
    INTEGER, INTENT(INOUT)       :: ktop (nCols)
    INTEGER, INTENT(INOUT)       :: kbcon(nCols)       

    !new-dec2011
    !----------------------------------------------------------------------------------
    REAL(KIND=r8)   , INTENT(OUT)      :: out1_mass_flux_up(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out2_mass_flux_d(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out3_detr(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out4_liq_water(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out5_qc(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out6_wf(nCols)
    REAL(KIND=r8)   , INTENT(OUT)      :: out7_entr(nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out8_he(nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out9_hes(nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out10_hc(nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT)      :: out11_rh(nCols, kMax)

    !
    ! LOCAL VARIABLES
    !
    INTEGER            :: i
    INTEGER            :: k
    INTEGER            :: j
    !not used    INTEGER            :: iedt
    INTEGER            :: istart 
    INTEGER            :: iend

    INTEGER            :: kdet1     (nCols)
    INTEGER            :: kdet      (nCols)
    REAL(KIND=r8)               :: mconv     (nCols)
    INTEGER            :: kzdown    (nCols)
    INTEGER            :: kbmax     (nCols)
    INTEGER            :: k22       (nCols)
    INTEGER            :: jmin      (nCols)
    INTEGER            :: kstabi    (nCols)
    INTEGER            :: kstabm    (nCols)
    INTEGER            :: KZI       (nCols)

    REAL(KIND=r8)               :: aaeq      (nCols)
    REAL(KIND=r8)               :: edt       (nCols)
    REAL(KIND=r8)               :: aa0       (nCols)  !at time t
    REAL(KIND=r8)               :: aa1       (nCols)  ! at time t+1
    REAL(KIND=r8)               :: hkb       (nCols)
    REAL(KIND=r8)               :: hkbo      (nCols)
    REAL(KIND=r8)               :: pwav      (nCols)
    REAL(KIND=r8)               :: pwev      (nCols)
    REAL(KIND=r8)               :: pwavo     (nCols)
    REAL(KIND=r8)               :: pwevo     (nCols)
    REAL(KIND=r8)               :: bu        (nCols)
    REAL(KIND=r8)               :: cap_max   (nCols)
    REAL(KIND=r8)               :: vshear    (nCols)
    REAL(KIND=r8)               :: sdp       (nCols)
    REAL(KIND=r8)               :: vws       (nCols)
    REAL(KIND=r8)               :: he        (nCols,kMax)
    REAL(KIND=r8)               :: hes       (nCols,kMax)
    REAL(KIND=r8)               :: qes       (nCols,kMax)
    REAL(KIND=r8)               :: z         (nCols,kMax)
    REAL(KIND=r8)               :: dby       (nCols,kMax)
    REAL(KIND=r8)               :: qc        (nCols,kMax)
    REAL(KIND=r8)               :: qrcd      (nCols,kMax)
    REAL(KIND=r8)               :: pwd       (nCols,kMax)
    REAL(KIND=r8)               :: pw        (nCols,kMax)
    REAL(KIND=r8)               :: heo       (nCols,kMax)
    REAL(KIND=r8)               :: heso      (nCols,kMax)
    REAL(KIND=r8)               :: qeso      (nCols,kMax)
    REAL(KIND=r8)               :: zo        (nCols,kMax)
    REAL(KIND=r8)               :: dbyo      (nCols,kMax)
    REAL(KIND=r8)               :: qco       (nCols,kMax)
    REAL(KIND=r8)               :: qrcdo     (nCols,kMax)
    REAL(KIND=r8)               :: pwdo      (nCols,kMax)
    REAL(KIND=r8)               :: pwo       (nCols,kMax)
    REAL(KIND=r8)               :: hcd       (nCols,kMax)
    REAL(KIND=r8)               :: hcdo      (nCols,kMax)
    REAL(KIND=r8)               :: qcd       (nCols,kMax)
    REAL(KIND=r8)               :: qcdo      (nCols,kMax)
    REAL(KIND=r8)               :: dbyd      (nCols,kMax)
    REAL(KIND=r8)               :: dbydo     (nCols,kMax)
    REAL(KIND=r8)               :: hc        (nCols,kMax)
    REAL(KIND=r8)               :: hco       (nCols,kMax)
    REAL(KIND=r8)               :: qrc       (nCols,kMax)
    REAL(KIND=r8)               :: qrco      (nCols,kMax)
    REAL(KIND=r8)               :: zu        (nCols,kMax)
    REAL(KIND=r8)               :: zuo       (nCols,kMax)
    REAL(KIND=r8)               :: zd        (nCols,kMax)
    REAL(KIND=r8)               :: zdo       (nCols,kMax)
    REAL(KIND=r8)               :: qes_cup   (nCols,kMax)
    REAL(KIND=r8)               :: q_cup     (nCols,kMax)
    REAL(KIND=r8)               :: he_cup    (nCols,kMax)
    REAL(KIND=r8)               :: hes_cup   (nCols,kMax)
    REAL(KIND=r8)               :: z_cup     (nCols,kMax)
    REAL(KIND=r8)               :: p_cup     (nCols,kMax)
    REAL(KIND=r8)               :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)               :: t_cup     (nCols,kMax)
    REAL(KIND=r8)               :: qeso_cup  (nCols,kMax)
    REAL(KIND=r8)               :: qo_cup    (nCols,kMax)
    REAL(KIND=r8)               :: heo_cup   (nCols,kMax)
    REAL(KIND=r8)               :: heso_cup  (nCols,kMax)
    REAL(KIND=r8)               :: zo_cup    (nCols,kMax)
    REAL(KIND=r8)               :: po_cup    (nCols,kMax)
    REAL(KIND=r8)               :: gammao_cup(nCols,kMax)
    REAL(KIND=r8)               :: tn_cup    (nCols,kMax)
    REAL(KIND=r8)               :: cd        (nCols,kMax)
    REAL(KIND=r8)               :: cdd       (nCols,kMax)

    REAL(KIND=r8)               :: dellat_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)               :: dellaq_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)               :: dellaqc_ens(nCols,kMax, maxens2)
    REAL(KIND=r8)               :: pwo_ens    (nCols,kMax, maxens2)
    REAL(KIND=r8)               :: xf_ens     (nCols,ensdim)
    REAL(KIND=r8)               :: outt_ens   (nCols,ensdim)
    REAL(KIND=r8)               :: pr_ens     (nCols,ensdim)
    REAL(KIND=r8)               :: edtc     (nCols,maxens2)
    REAL(KIND=r8)               :: dq 
    REAL(KIND=r8)               :: mbdt
    REAL(KIND=r8)               :: zktop
    !-------
    !new
    !---------
    REAL(KIND=r8)               :: dh2         (nCols) 
    REAL(KIND=r8)               :: xfac1       (nCols) 
    REAL(KIND=r8)               :: xfac_for_dn (nCols)
    INTEGER            :: left        (nCols)
    INTEGER            :: nLeft,ib   !not used ,nNewLeft
    INTEGER            :: maxens22
    INTEGER            :: kk1(nCols),kk
    REAL(KIND=r8)               ::tke1D(kMax)
    REAL(KIND=r8)               ::tkeMax(nCols)
    REAL(KIND=r8)               ::tkeMedia(nCols)
    REAL(KIND=r8)               :: den0 (nCols,kMax)
    REAL(KIND=r8)               :: den1 (nCols,kMax)
    REAL(KIND=r8)               :: da,dz
    REAL(KIND=r8)               :: hh(nCols)
    REAL(KIND=r8)               :: cine(nCols)

    REAL(KIND=r8)               :: rh(nCols,kMax)
    REAL(KIND=r8)               :: bb,cc,dd,ee,ff
    REAL(KIND=r8)               :: entr2D(nCols,kMax)
    REAL(KIND=r8)               :: detr2D(nCols,kMax)
    REAL(KIND=r8)               :: entr2,Hdeep
    REAL(KIND=r8)               :: factor
    INTEGER                     :: kentra,jdt24 !nilo1
    REAL(KIND=r8)               :: ucape0 (nCols), ucape1 (nCols),dcape0 (nCols), dcape1 (nCols)
    REAL(KIND=r8)               :: d_ucape (nCols),d_dcape(nCols)
    REAL(KIND=r8)               :: cine0(nCols),cine1(nCols)
    REAL(KIND=r8)               :: aaa(nCols)
    REAL(KIND=r8)               :: k_zhang(nCols),dcape_zhang(nCols)
    REAL(KIND=r8)               :: k_zhang1(nCols),cape_zhang1(nCols)

    !
    ! Compress Local Variable 
    !
    INTEGER :: nCols_gz
    REAL(KIND=r8)    :: edtc_gz       (nCols,maxens2)
    INTEGER :: ierr_gz       (nCols)
    REAL(KIND=r8)    :: dellat_ens_gz (nCols,kMax, maxens2)
    REAL(KIND=r8)    :: dellaq_ens_gz (nCols,kMax, maxens2)
    REAL(KIND=r8)    :: dellaqc_ens_gz(nCols,kMax, maxens2)
    REAL(KIND=r8)    :: pwo_ens_gz    (nCols,kMax, maxens2)
    REAL(KIND=r8)    :: heo_cup_gz    (nCols,kMax)
    REAL(KIND=r8)    :: zo_cup_gz     (nCols,kMax)
    REAL(KIND=r8)    :: po_cup_gz     (nCols,kMax)
    REAL(KIND=r8)    :: hcdo_gz       (nCols,kMax)
    REAL(KIND=r8)    :: zdo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: cdd_gz        (nCols,kMax)
    REAL(KIND=r8)    :: heo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: qo_cup_gz     (nCols,kMax)
    REAL(KIND=r8)    :: qrcdo_gz      (nCols,kMax)
    REAL(KIND=r8)    :: qo_gz         (nCols,kMax)
    REAL(KIND=r8)    :: zuo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: cd_gz         (nCols,kMax)
    REAL(KIND=r8)    :: hco_gz        (nCols,kMax)
    INTEGER :: ktop_gz       (nCols)
    INTEGER :: k22_gz        (nCols)
    INTEGER :: kbcon_gz      (nCols)
    INTEGER :: jmin_gz       (nCols)
    INTEGER :: kdet_gz       (nCols)
    REAL(KIND=r8)    :: qco_gz        (nCols,kMax)
    REAL(KIND=r8)    :: qrco_gz       (nCols,kMax)
    REAL(KIND=r8)    :: tn_gz         (nCols,kMax)
    REAL(KIND=r8)    :: po_gz         (nCols,kMax)
    REAL(KIND=r8)    :: z1_gz         (nCols)
    REAL(KIND=r8)    :: psur_gz       (nCols)
    REAL(KIND=r8)    :: gamma_cup_gz  (nCols,kMax)
    REAL(KIND=r8)    :: pr_ens_gz     (nCols,ensdim)
    REAL(KIND=r8)    :: pwo_gz        (nCols,kMax)
    REAL(KIND=r8)    :: pwdo_gz       (nCols,kMax)
    REAL(KIND=r8)    :: outt_ens_gz   (nCols,ensdim)    
    REAL(KIND=r8)    :: he_cup_gz     (nCols,kMax)
    INTEGER :: kbmax_gz      (nCols)
    REAL(KIND=r8)    :: heso_cup_gz   (nCols,kMax)
    REAL(KIND=r8)    :: cap_max_gz    (nCols)
    REAL(KIND=r8)    :: aa0_gz        (nCols)
    REAL(KIND=r8)    :: aa1_gz        (nCols)
    REAL(KIND=r8)    :: xmb_gz        (nCols)
    REAL(KIND=r8)    :: xf_ens_gz     (nCols,ensdim)
    INTEGER :: mask_gz       (nCols)
    REAL(KIND=r8)    :: mconv_gz      (nCols)
    REAL(KIND=r8)    :: omeg_gz       (nCols,kMax)
    REAL(KIND=r8)    :: massfln_gz    (nCols,ensdim)
    REAL(KIND=r8)    :: p_cup_gz      (nCols,kMax)
    INTEGER :: listim  (nCols)
    LOGICAL :: bitx    (nCols)
    INTEGER :: litx    (nCols)
    REAL(KIND=r8)    :: entr2D_gz   (nCols,kMax)
    REAL(KIND=r8)    :: detr2D_gz   (nCols,kMax)
    REAL(KIND=r8)    :: dcape_zhang_gz      (nCols)
    REAL(KIND=r8)    :: k_zhang_gz      (nCols)
    REAL(KIND=r8)    :: dcape_gz      (nCols)
    REAL(KIND=r8)    :: d_ucape_gz      (nCols)
    REAL(KIND=r8)    :: d_dcape_gz      (nCols)
    REAL(KIND=r8)    :: cine_gz      (nCols)
    REAL(KIND=r8)    :: sens_gz      (nCols)
    REAL(KIND=r8)    :: tkeMax_gz     (nCols)
    REAL(KIND=r8)    :: tkeMedia_gz     (nCols)
    REAL(KIND=r8)    :: den1_gz       (nCols,kMax)


    !
    !***************** the following are your basic environmental
    !                  variables. They carry a "_cup" if they are
    !                  on model cloud levels (staggered). They carry
    !                  an "o"-ending (z becomes zo), if they are the forced
    !                  variables. They are preceded by x (z becomes xz)
    !                  to indicate modification by some typ of cloud
    !
    ! z           = heights of model levels
    ! q           = environmental mixing ratio
    ! qes         = environmental saturation mixing ratio
    ! t           = environmental temp
    ! p           = environmental pressure
    ! he          = environmental moist static energy
    ! hes         = environmental saturation moist static energy
    ! z_cup       = heights of model cloud levels
    ! q_cup       = environmental q on model cloud levels
    ! qes_cup     = saturation q on model cloud levels
    ! t_cup       = temperature (Kelvin) on model cloud levels
    ! p_cup       = environmental pressure
    ! he_cup = moist static energy on model cloud levels
    ! hes_cup = saturation moist static energy on model cloud levels
    ! gamma_cup = gamma on model cloud levels
    !
    !
    ! hcd = moist static energy in downdraft
    ! zd normalized downdraft mass flux
    ! dby = buoancy term
    ! entr = entrainment rate
    ! zd   = downdraft normalized mass flux
    ! entr= entrainment rate
    ! hcd = h in model cloud
    ! bu = buoancy term
    ! zd = normalized downdraft mass flux
    ! gamma_cup = gamma on model cloud levels
    ! mentr_rate = entrainment rate
    ! qcd = cloud q (including liquid water) after entrainment
    ! qrch = saturation q in cloud
    ! pwd = evaporate at that level
    ! pwev = total normalized integrated evaoprate (I2)
    ! entr= entrainment rate
    ! z1 = terrain elevation
    ! entr = downdraft entrainment rate
    ! jmin = downdraft originating level
    ! kdet = level above ground where downdraft start detraining
    ! psur        = surface pressure
    ! z1          = terrain elevation
    ! pr_ens = precipitation ensemble
    ! xf_ens = mass flux ensembles
    ! massfln = downdraft mass flux ensembles used in next timestep
    ! omeg = omega from large scale model
    ! mconv = moisture convergence from large scale model
    ! zd      = downdraft normalized mass flux
    ! zu      = updraft normalized mass flux
    ! dir     = "storm motion"
    ! mbdt    = arbitrary numerical parameter
    ! dtime   = dt over which forcing is applied
    ! iact_gr_old = flag to tell where convection was active
    ! kbcon       = LFC of parcel from k22
    ! k22         = updraft originating level
    ! icoic       = flag if only want one closure (usually set to zero!)
    ! dby = buoancy term
    ! ktop = cloud top (output)
    ! xmb    = total base mass flux
    ! hc = cloud moist static energy
    ! hkb = moist static energy at originating level
    ! mentr_rate = entrainment rate

    !
    !snf parameter from namelist 
    !
    ! begin executable
    !
    outqc   =0.0_r8
    kzdown  =0
    kbmax   =0
    k22     =0
    jmin    =0
    kstabi  =0
    kstabm  =0
    KZI     =0
    left    =0
    aaeq    =0.0_r8
    edt     =0.0_r8
    aa1     =0.0_r8
    aa0     =0.0_r8
    hkb     =0.0_r8
    hkbo    =0.0_r8
    xmb     =0.0_r8
    pwav    =0.0_r8
    pwev    =0.0_r8
    pwavo   =0.0_r8
    pwevo   =0.0_r8
    bu      =0.0_r8
    cap_max =0.0_r8
    vshear  =0.0_r8
    sdp     =0.0_r8
    vws     =0.0_r8
    he      =0.0_r8
    hes     =0.0_r8
    qes     =0.0_r8
    z       =0.0_r8
    dby     =0.0_r8
    qc      =0.0_r8
    qrcd    =0.0_r8
    pwd     =0.0_r8
    pw      =0.0_r8
    heo     =0.0_r8
    heso    =0.0_r8
    qeso    =0.0_r8
    zo      =0.0_r8
    dbyo    =0.0_r8
    qco     =0.0_r8
    qrcdo  =0.0_r8
    pwdo  =0.0_r8
    pwo  =0.0_r8
    hcd  =0.0_r8
    hcdo  =0.0_r8
    qcd  =0.0_r8
    qcdo  =0.0_r8
    dbyd  =0.0_r8
    dbydo  =0.0_r8
    hc  =0.0_r8
    hco  =0.0_r8
    qrc  =0.0_r8
    qrco  =0.0_r8
    zu  =0.0_r8
    zuo  =0.0_r8
    zd  =0.0_r8
    zdo  =0.0_r8
    qes_cup   =0.0_r8
    q_cup  =0.0_r8
    he_cup    =0.0_r8
    hes_cup   =0.0_r8
    z_cup  =0.0_r8
    p_cup  =0.0_r8
    gamma_cup =0.0_r8
    t_cup  =0.0_r8
    qeso_cup  =0.0_r8
    qo_cup    =0.0_r8
    heo_cup   =0.0_r8
    heso_cup  =0.0_r8
    zo_cup    =0.0_r8
    po_cup    =0.0_r8
    gammao_cup=0.0_r8
    tn_cup    =0.0_r8
    cd  =0.0_r8
    cdd  =0.0_r8
    dellat_ens  =0.0_r8
    dellaq_ens  =0.0_r8
    dellaqc_ens =0.0_r8
    pwo_ens  =0.0_r8   
    xf_ens   =0.0_r8   
    outt_ens  =0.0_r8  
    pr_ens   =0.0_r8   
    edtc =0.0_r8
    dq =0.0_r8
    mbdt=0.0_r8
    zktop=0.0_r8



    dh2   =0.0_r8
    xfac1   =0.0_r8
    xfac_for_dn=0.0_r8
    dellat_ens =0.0_r8
    dellaq_ens =0.0_r8
    dellaqc_ens=0.0_r8
    pwo_ens    =0.0_r8
    hkbo=0
    istart=1
    iend=nCols
    maxens22=3 !!! 
    massfln=0.0_r8
    mconv=0.0_r8
    qrco=0.0_r8
    qrco_gz=0.0_r8
    DO i=istart,iend
       !
       ! prepare input, erase output
       !
       kdet  (i) =2
       kdet1 (i) =0
       pre   (I) =0.0_r8
    END DO

    !new
    aaa=0.0_r8 
    !
    ! calculate moisture convergence mconv
    !
    DO k=2,kMax-1
       DO i = istart,iend
          !dq      = 0.5_r8*(q(i,k+1)-q(i,k-1))
          dq      = 0.5_r8*(qo(i,k+1)-qo(i,k-1))
          mconv(i) = mconv(i) + omeg(i,k)*dq/9.81_r8
       END DO
    END DO

    DO I = istart,iend
       IF(mconv(I) < 0.0_r8)  mconv(I) = 0.0_r8
    END DO
    !
    ! initial detrainment rates
    !
    DO k=1,kMax
       DO i=istart,iend
          cd (i,k) = 0.1_r8*entr_rate            
          cdd(i,k) = 0.0_r8
       END DO
    END DO

    DO i=istart,iend
       aa0   (i)=0.0_r8                                   
       aa1   (i)=0.0_r8
       hh(i)=0.0_r8
       cine (i)=0.0_r8
       kstabm(i)=kMax-2
       aaeq  (i)=0.0_r8                              
       IF(aaeq(i) <  0.0_r8)THEN
          ierr(i)=20
       ELSE
          ierr (i)=0
       END IF
       tkeMedia(i)=0.0_r8
    END DO
    !
    DO i=istart,iend
       cap_max(i)=cap_maxs
       IF(mask(i).NE.1)cap_max(i)=cap_maxs_land
    END DO
    mbdt=(REAL(1,kind=r8)-3.0_r8)*dtime*1.e-3_r8 + dtime*5.e-03_r8  
    !
    ! environmental conditions, FIRST HEIGHTS
    !
    DO k=1,maxens*maxens22*maxens3
       DO i=istart,iend
          IF(ierr(i).NE.20)THEN
             xf_ens  (i,(iens-1)*maxens*maxens22*maxens3+k)= 0.0_r8
             pr_ens  (i,(iens-1)*maxens*maxens22*maxens3+k)= 0.0_r8
             outt_ens(i,(iens-1)*maxens*maxens22*maxens3+k)= 0.0_r8
          END IF
       END DO
    END DO
    !
    ! calculate moist static energy, heights, qes
    !
    ! at time t
    CALL cup_env(z      , & ! z      (out)
         qes    , & ! qes    (out)
         he     , & ! he     (inout)
         hes    , & ! hes    (out)
         t      , & ! t      (in)
         q      , & ! q      (inout)
         p      , & ! p      (in)
         z1     , & ! z1     (in)
         nCols  , & ! nCols  (in)
         kMax   , & ! kMax   (in)
         istart , & ! istart (in)
         iend   , & ! iend   (in)
         psur   , & ! psur   (in)
         ierr   , & ! ierr   (in)
         0      ,&  !tcrit  (in)
         den0     ) ! den    (out)
    !at time t+1
    CALL cup_env(zo     , & ! zo     (out)
         qeso   , & ! qeso   (out)
         heo    , & ! heo    (inout)
         heso   , & ! heso   (out)
         tn     , & ! tn     (in)
         qo     , & ! qo     (inout)
         po     , & ! po     (in)
         z1     , & ! z1     (in)
         nCols  , & ! nCols  (in)
         kMax   , & ! kMax   (in)
         istart , & ! istart (in)
         iend   , & ! iend   (in)
         psur   , & ! psur   (in)
         ierr   , & ! ierr   (in)
         0      , & ! tcrit  (in)
         den1      ) ! den    (out)
    !
    ! environmental values on cloud levels
    !
    !at time t
    CALL cup_env_clev(t        , & ! t         (in)
         qes      , & ! qes       (in)
         q        , & ! q         (in)
         he       , & ! he        (in)
         hes      , & ! hes       (in)
         z        , & ! z         (in)
         p        , & ! p         (in)
         qes_cup  , & ! qes_cup   (out)
         q_cup    , & ! q_cup     (out)
         he_cup   , & ! he_cup    (out)
         hes_cup  , & ! hes_cup   (out)
         z_cup    , & ! z_cup     (out)
         p_cup    , & ! p_cup     (out)
         gamma_cup, & ! gamma_cup (out)
         t_cup    , & ! t_cup     (out)
         psur     , & ! psur      (in)
         nCols    , & ! nCols     (in)
         kMax     , & ! kMax      (in)
         istart   , & ! istart    (in)
         iend     , & ! iend      (in)
         ierr     , & ! ierr      (in)
         z1         ) ! z1        (in)

    !at time t+1

    CALL cup_env_clev(tn        , &! tn        (in)
         qeso      , &! qeso      (in)
         qo        , &! qo        (in)
         heo       , &! heo       (in)
         heso      , &! heso      (in)
         zo        , &! zo        (in)
         po        , &! po        (in)
         qeso_cup  , &! qeso_cup  (out)
         qo_cup    , &! qo_cup    (out)
         heo_cup   , &! heo_cup   (out)
         heso_cup  , &! heso_cup  (out)
         zo_cup    , &! zo_cup    (out)
         po_cup    , &! po_cup    (out)
         gammao_cup, &! gammao_cup(out)
         tn_cup    , &! tn_cup    (out)
         psur      , &! psur      (in)
         nCols     , &! nCols     (in)
         kMax      , &! kMax      (in)
         istart    , &! istart    (in)
         iend      , &! iend      (in)
         ierr      , &! ierr      (in)
         z1          )! z1        (in)
    !
    !
    !
    kbmax=0
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. zo_cup(i,k) >  zkbmax+z1(i) .AND. kbmax(i) ==0)THEN
             kbmax(i)=k
          END IF
       END DO
    END DO
    !
    ! level where detrainment for downdraft starts
    !
    kdet1=0
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. zo_cup(i,k) >  z_detr+z1(i) .AND.kdet1(i) ==0)THEN
             kdet (i)=k
             kdet1(i)=k
          END IF
       END DO
    END DO
    !
    ! determine level with highest moist static energy content - k22
    ! kstart = 3 
    !
    CALL cup_maximi(heo_cup  , &  ! heo_cup (in)
         nCols    , &  ! nCols   (in)
         kMax     , &  ! kMax    (in) 
         3        , &  ! ks      (in) !era 3  
         kbmax    , &  ! kbmax   (in)
         k22      , &  ! k22     (out)
         istart   , &  ! istart  (in)
         iend     , &  ! iend    (in)
         ierr       )  ! ierr    (in)
    !                
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          kzi(i) = 1
          IF(k22(i) >= kbmax(i))ierr(i)=2
       END IF
    END DO
    !
    ! determine the level of convective cloud base  - kbcon
    !-------------------------------------
    !snf  call cup_kbcon for cap_max=cap_max-(1-1)*cap_max_increment
    !---------------
    !
    CALL cup_kbcon(&
         cap_max_increment, & ! cap_max_increment (in)
         1                , & ! iloop             (in)
         k22              , & ! k22               (inout)
         kbcon            , & ! kbcon             (out)
         heo_cup          , & ! heo_cup           (in)
         heso_cup         , & ! heso_cup          (in)
         nCols            , & ! nCols             (in)
         kMax             , & ! kMax              (in)
         istart           , & ! istart            (in)
         iend             , & ! iend              (in)
         ierr             , & ! ierr              (inout)
         kbmax            , & ! kbmax             (in)
         po_cup           , & ! po_cup            (in)
         cap_max            ) ! cap_max           (in)
    DO I=ISTART,IEND
       IF(ierr(I) == 0)THEN
          hkb(i)=hkbo(i)
       END IF
    END DO
    !----------------------------------new
    tkeMax=0.0_r8
    DO i=istart,iend
       tke1D=0.0_r8
       IF(ierr(I) == 0)THEN
          DO k=2,kbcon(i)
             tke1D(k)=tke(i,k)
          ENDDO
          tkeMax(i)=MAXVAL(tke1D)
       END IF
    ENDDO
    !-------------------------------------------
    ! calculate PBL TOP level
    !-----------------------------
    DO i=istart,iend
       kk1(i)=kbcon(i)
       DO k=2,kbcon(i)
          IF(tke(i,k)<=0.17_r8)THEN   !!
             kk1(i)=k
             EXIT
          ENDIF
       ENDDO
    ENDDO
    !---------------------------------------

    tkeMedia=0.0_r8
    DO I=ISTART,IEND
       IF(ierr(I) == 0)THEN
          DO k=1,kbcon(i)
             tkeMedia(i)=tkeMedia(i)+tke(i,k)
          ENDDO
          IF(kbcon(i)>=1)tkeMedia(i)=tkeMedia(i)/float(kbcon(i))
       END IF
    END DO
    !---------------------------------------------------------------
    !
    ! increase detrainment in stable layers
    !
    CALL cup_minimi( &
         heso_cup , &  ! heso_cup (in)
         nCols    , &  ! nCols    (in)
         kMax     , &  ! kMax     (in)
         kbcon    , &  ! kbcon    (in)
         kstabm   , &  ! kstabm   (in)
         kstabi   , &  ! kstabi   (out)
         istart   , &  ! istart   (in)
         iend     , &  ! iend     (in)
         ierr       )  ! ierr     (in)
    !--------------------------
    ! define entr and detr rate    nilo1
    !---------------------------
    rh=0.0_r8
    entr2D=0.0_r8
    detr2D=0.0_r8
    entr2=0.0_r8
    DO k=1,kMax
       DO i=istart,iend
          rh(i,k)=qo_cup(i,k)/qeso_cup(i,k)
          IF(rh(i,k)<=1.0e-8_r8)rh(i,k)=1.0e-8_r8 
          entr2D(i,k)=entr_rate
          detr2D(i,k)=mentr_rate
          entr2=entr_rate
       ENDDO
    ENDDO
    !--------------------------------------
    ! Below for NWP to reduce excess precipitation at first day.
    jdt24=MAX(INT(24.0_r8*3600.0_r8/dtime),2)   !integer
    factor=1.0_r8
    IF(jdt<jdt24)THEN
       !factor=2.0_r8-(jdt-1)*1.25_r8/(jdt24-1)
       ee=1.0_r8
       factor=3.0_r8+(jdt-1)*(ee-3.0_r8)/(jdt24-1)
    ENDIF
    !!factor=1.0_r8

    kentra=1
    IF(kentra==1)THEN   !! GRELL_DEVENYI 2002

       DO k=1,kMax
          DO i=istart,iend
             entr2D(i,k)=entr_rate
             detr2D(i,k)=mentr_rate
             entr2=entr_rate
          ENDDO
       ENDDO
    ELSEIF(kentra==2)THEN  !!!! ECMWF-BETCHTOLD-2008
       DO k=1,kMax
          DO i=istart,iend
             IF( ierr(i) == 0.AND.k>=kbcon(i)) THEN
                cc=qeso_cup(i,k)/qeso_cup(i,kbcon(i))
                bb=1.8e-3_r8 !3.0e-3_r8,2.0e-5_r8, 5.0e-5_r8 which depends from ee
                ee=3.0_r8    !(3,2,1, 0)
                !!entr2D(i,k)=bb*(1.3_r8-rh(i,k))*cc**3
                entr2D(i,k)=bb*(1.3_r8-rh(i,k))*cc**ee
                entr2D(i,k)=MAX(MIN(entr2D(i,k),1.0e-1_r8),0.0_r8)
             ENDIF
          END DO
       END DO
    ELSEIF(kentra==3)THEN  !NEW NILO 2015
       DO k=1,kMax
          DO i=istart,iend
             IF( ierr(i) == 0.AND.k>=kbcon(i)) THEN
                bb=zo_cup(i,k)-zo_cup(i,kbcon(i)-1)
                dd=(zo_cup(i,k)-zo_cup(i,1))/(radius*(rh(i,k)+0.1_r8))  
                cc=MIN(30.0_r8,dd)
                entr2D(i,k)=0.5_r8*2.0e-5_r8*(1.0e-4_r8*EXP(cc)+2.0e+3_r8/bb)
                IF(mask(i)==1)entr2D(i,k)=factor*entr2D(i,k)
                entr2D(i,k)=MAX(MIN(entr2D(i,k),1.0e-2_r8),0.0_r8)  
             ENDIF
          END DO
       END DO
    ENDIF
    !-------------------------------------------------------

    DO k=1,kMax
       DO i=istart,iend
          IF( ierr(i) == 0 .AND. kstabm(i)-1 >= kstabi(i) .AND. &
               k >= kstabi(i) .AND. k<= kstabm(i)-1 )THEN
             entr2=entr2D(i,k)
             cd(i,k)=cd(i,k-1)+1.5_r8*entr2  
!!!cd(i,k)=cd(i,k-1)+2.0_r8*entr2
             IF(iens >  4)THEN
                cd(i,k)=cd(i,k-1)+REAL(iens-4,kind=r8)*entr2&
                     /REAL(kstabm(i)-kstabi(i),kind=r8)
             ELSE
                cd(i,k)=cd(i,k)
             END IF
             IF(cd(i,k) >  10.0_r8*entr2) cd(i,k)=10.0_r8*entr2 !new
          END IF
       END DO
    END DO
    !
    ! calculate incloud moist static energy
    !
    !at time t
    CALL cup_up_he(  &
         k22       , & ! k22        (in)
         hkb       , & ! hkb        (out)
         z_cup     , & ! z_cup      (in)
         cd        , & ! cd         (in)
         mentr_rate, & ! mentr_rate (in)
         he_cup    , & ! he_cup     (in)
         hc        , & ! hc         (out)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         kbcon     , & ! kbcon      (in)
         ierr      , & ! ierr       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         dby       , & ! dby        (out)
         he        , & ! he         (in)
         hes_cup   , & !  ) ! hes_cup    (in)
         entr2D)
    !
    !at time t+1
    CALL cup_up_he(  &
         k22       , & ! k22        (in)
         hkbo      , & ! hkbo       (out)
         zo_cup    , & ! zo_cup     (in)
         cd        , & ! cd         (in)
         mentr_rate, & ! mentr_rate (in)
         heo_cup   , & ! heo_cup    (in)
         hco       , & ! hco        (out)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         kbcon     , & ! kbcon      (in)
         ierr      , & ! ierr       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         dbyo      , & ! dbyo       (out)
         heo       , & ! heo        (in)
         heso_cup  , & !
         entr2D)
    !

    ! determine cloud top - ktop
    !
    CALL cup_ktop(  &
         1        , & ! ilo    (in)
         dbyo     , & ! dbyo   (inout)
         kbcon    , & ! kbcon  (in)
         ktop     , & ! ktop   (out)
         nCols    , & ! nCols  (in)
         kMax     , & ! kMax   (in)
         istart   , & ! istart (in)
         iend     , & ! iend   (in)
         ierr       ) ! ierr   (inout)

    kzdown(istart:iend)=0
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             zktop=(zo_cup(i,ktop(i))-z1(i))*0.6_r8
             zktop=MIN(zktop+z1(i),zcutdown+z1(i))
             IF(zo_cup(i,k) >  zktop .AND. kzdown(i) == 0 )THEN
                kzdown(i) = k
             END IF
          END IF
       END DO
    END DO
    !
    ! downdraft originating level - jmin
    ! jmin output from cup_minimi
    !
    CALL cup_minimi( &
         heso_cup  , &! heso_cup (in)
         nCols     , &! nCols    (in)
         kMax      , &! kMax     (in)
         k22       , &! k22      (in)
         kzdown    , &! kzdown   (in)
         jmin      , &! jmin     (out)
         istart    , &! istart   (in)
         iend      , &! iend     (in)
         ierr        )! ierr     (in)
    !
    ! check whether it would have buoyancy, if there where
    ! no entrainment/detrainment
    !
    dh2(istart:iend)=0.0_r8    
    nLeft = 0
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          nLeft = nLeft + 1
          left(nLeft) = i
          IF (jmin(i)-1 <  kdet(i)) kdet(i)=jmin(i)-1
          IF (jmin(i) >= ktop(i)-1) jmin(i)=ktop(i)-2
       END IF
    END DO

    DO ib=1,nLeft
       i=left(ib)
101    CONTINUE
       DO k=jmin(i)-1,1,-1
          dh2(i)    = dh2(i) + (zo_cup  (i,k+1) - zo_cup  (i,k)) &
               * (heso_cup(i,jmin(i)) - heso_cup(i,k))
          IF(dh2(i) >  0.0_r8)THEN
             jmin(i)=jmin(i)-1
             IF(jmin(i) > 3)THEN
                IF (jmin(i)-1 <  kdet(i)  ) kdet(i)=jmin(i)-1
                IF (jmin(i)   >= ktop(i)-1) jmin(i)=ktop(i)-2
                dh2(i)=0.0_r8
                go to 101
             ELSE IF(jmin(i) <= 3 .AND. ierr(i) /= 9 )THEN
                ierr(i)=9
             END IF
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(jmin(i) <= 3 .AND. ierr(i) /= 4 .AND. dh2(i) <= 0.0_r8)THEN
             ierr(i)=4
          END IF
       END IF
    END DO
    !
    ! Must have at least depth_min m between cloud convective base
    ! and cloud top
    !
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)) <  depth_min)THEN
             ierr(i)=6
          END IF
       END IF
    END DO
    !
    ! normalized updraft mass flux profile
    !
    !
    CALL cup_up_nms( &
         zu        , & ! zu         (out)
         z_cup     , & ! z_cup      (in)
         mentr_rate, & ! mentr_rate (in)
         cd        , & ! cd         (in)
         kbcon     , & ! kbcon      (in)
         ktop      , & ! ktop       (in)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         ierr      , & ! ierr       (in)
         k22       , & ! ) ! k22        (in)
         entr2D)

    CALL cup_up_nms( &
         zuo       , & ! zuo        (out)
         zo_cup    , & ! zo_cup     (in)
         mentr_rate, & ! mentr_rate (in)
         cd        , & ! cd         (in)
         kbcon     , & ! kbcon      (in)
         ktop      , & ! ktop       (in)
         nCols     , & ! nCols      (in)
         kMax      , & ! kMax       (in)
         istart    , & ! istart     (in)
         iend      , & ! iend       (in)
         ierr      , & ! ierr       (in)
         k22       , & !  ) ! k22        (in)
         entr2D)
    !
    ! normalized downdraft mass flux profile,also work on bottom
    ! detrainment in this routin
    !
    CALL cup_dd_nms(  &
         zd         , & ! zd          (out)
         z_cup      , & ! z_cup       (in)
         cdd        , & ! cdd         (out)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (in)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         0          , & ! itest       (in)
         kdet       , & ! kdet        (in)
         z1         , & !  ) ! z1          (in)
         detr2D)

    CALL cup_dd_nms(  &
         zdo        , & ! zdo         (out)
         zo_cup     , & ! zo_cup      (in)
         cdd        , & ! cdd         (out)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (in)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         1          , & ! itest       (in)
         kdet       , & ! kdet        (in)
         z1         , & !  ) ! z1          (in)
         detr2D)
    !
    !  downdraft moist static energy
    !
    CALL cup_dd_he (  &
         hes_cup    , &! hes_cup     (in)
         hcd        , &! hcd         (out)
         z_cup      , &! z_cup       (in)
         cdd        , &! cdd         (in)
         mentrd_rate, &! mentrd_rate (in)
         jmin       , &! jmin        (in)
         ierr       , &! ierr        (in)
         nCols      , &! nCols       (in)
         kMax       , &! kMax        (in)
         istart     , &! istart      (in)
         iend       , &! iend        (in)
         he         , &! he          (in)
         dbyd       , & !  )! dbyd        (out)
         detr2D)
    CALL cup_dd_he (  &
         heso_cup   , &! heso_cup    (in)
         hcdo       , &! hcdo        (out)
         zo_cup     , &! zo_cup      (in)
         cdd        , &! cdd         (in)
         mentrd_rate, &! mentrd_rate (in)
         jmin       , &! jmin        (in)
         ierr       , &! ierr        (in)
         nCols      , &! nCols       (in)
         kMax       , &! kMax        (in)
         istart     , &! istart      (in)
         iend       , &! iend        (in)
         heo        , &! heo         (in)
         dbydo       , &! )! dbydo       (out)
         detr2D)

    !
    !  calculate moisture properties of downdraft
    !
    !
    !snf out  qcd = cloud q (including liquid water) after entrainment
    CALL cup_dd_moisture( &
         zd         , & ! zd          (in)
         hcd        , & ! hcd         (in)
         hes_cup    , & ! hes_cup     (in)
         qcd        , & ! qcd         (out)
         qes_cup    , & ! qes_cup     (in)
         pwd        , & ! pwd         (out)
         q_cup      , & ! q_cup       (in)
         z_cup      , & ! z_cup       (in)
         cdd        , & ! cdd         (in)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (inout)
         gamma_cup  , & ! gamma_cup   (in)
         pwev       , & ! pwev        (out)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         bu         , & ! bu          (out)
         qrcd       , & ! qrcd        (out)
         q          , & ! q           (in)
         2          , & !  ) ! iloop       (in)
         detr2D)
    CALL cup_dd_moisture( &
         zdo        , & ! zdo         (in)
         hcdo       , & ! hcdo        (in)
         heso_cup   , & ! heso_cup    (in)
         qcdo       , & ! qcdo        (out)
         qeso_cup   , & ! qeso_cup    (in)
         pwdo       , & ! pwdo        (out)
         qo_cup     , & ! qo_cup      (in)
         zo_cup     , & ! zo_cup      (in)
         cdd        , & ! cdd         (in)
         mentrd_rate, & ! mentrd_rate (in)
         jmin       , & ! jmin        (in)
         ierr       , & ! ierr        (inout)
         gammao_cup , & ! gammao_cup  (in)
         pwevo      , & ! pwevo       (out)
         nCols      , & ! nCols       (in)
         kMax       , & ! kMax        (in)
         istart     , & ! istart      (in)
         iend       , & ! iend        (in)
         bu         , & ! bu          (out)
         qrcdo      , & ! qrcdo       (out)
         qo         , & ! qo          (in)
         1          , & !  ) ! iloop       (in)
         detr2D)
    !
    ! calculate moisture properties of updraft
    !
    !snf
    !OUT
    ! qc = cloud q (including liquid water) after entrainment
    ! qrc = liquid water content in cloud after rainout
    ! pw = condensate that will fall out at that level


    CALL cup_up_moisture( &
         ierr       , & ! ierr       (in)
         z_cup      , & ! z_cup      (in)
         qc         , & ! qc         (out)
         qrc        , & ! qrc        (out)
         pw         , & ! pw         (out)
         pwav       , & ! pwav       (out)
         kbcon      , & ! kbcon      (in)
         ktop       , & ! ktop       (in)
         nCols      , & ! nCols      (in)
         kMax       , & ! kMax       (in)
         istart     , & ! istart     (in)
         iend       , & ! iend       (in)
         cd         , & ! cd         (in)
         dby        , & ! dby        (inout)
         mentr_rate , & ! mentr_rate (in)
         q          , & ! q          (in)
         gamma_cup  , & ! gamma_cup  (in)
         zu         , & ! zu         (in)
         qes_cup    , & ! qes_cup    (in)
         k22        , & ! k22        (in)
         q_cup      , & !  ) ! q_cup      (in)
         entr2D)

    CALL cup_up_moisture( &
         ierr       , & ! ierr       (in)
         zo_cup     , & ! zo_cup     (in)
         qco        , & ! qco        (out)
         qrco       , & ! qrco       (out)
         pwo        , & ! pwo        (out)
         pwavo      , & ! pwavo      (out)
         kbcon      , & ! kbcon      (in)
         ktop       , & ! ktop       (in)
         nCols      , & ! nCols      (in)
         kMax       , & ! kMax       (in)
         istart     , & ! istart     (in)
         iend       , & ! iend       (in)
         cd         , & ! cd         (in)
         dbyo       , & ! dbyo       (inout)
         mentr_rate , & ! mentr_rate (in)
         q          , & ! q          (in)
         gammao_cup , & ! gammao_cup (in)
         zuo        , & ! zuo        (in)
         qeso_cup   , & ! qeso_cup   (in)
         k22        , & ! k22        (in)
         qo_cup     , & !  ) ! qo_cup     (in)
         entr2D)


    !
    ! calculate workfunctions for updrafts
    !
    CALL cup_up_aa0(  &
         aa0        , & ! aa0       (inout)
         z          , & ! z         (in)
         zu         , & ! zu        (in)
         dby        , & ! dby       (in)
         gamma_cup  , & ! gamma_cup (in)
         t_cup      , & ! t_cup     (in)
         kbcon      , & ! kbcon     (in)
         ktop       , & ! ktop      (in)
         kMax       , & ! kMax      (in)
         nCols      , & ! nCols     (in)
         istart     , & ! istart    (in)
         iend       , & ! iend      (in)
         ierr         ) ! ierr      (inout)

    CALL cup_up_aa0(  &
         aa1        , & ! aa1       (inout)
         zo         , & ! z0        (in)
         zuo        , & ! zu0       (in)
         dbyo       , & ! dbyo      (in)
         gammao_cup , & ! gammao_cup(in)
         tn_cup     , & ! tn_cup    (in)
         kbcon      , & ! kbcon     (in)
         ktop       , & ! ktop      (in)
         kMax       , & ! kMax      (in)
         nCols      , & ! nCols     (in)
         istart     , & ! istart    (in)
         iend       , & ! iend      (in)
         ierr         ) ! ierr      (inout)
    d_ucape=0.0_r8
    d_dcape=0.0_r8
    !-------------------
    ! deleted cup_up_cine
    ! undiluted and diluted CAPE deleted
    !-------------------
    DO i=istart,iend
       cine(i)=cine3(i)
       dcape1(i)=0.0_r8
       dcape0(i)=0.0_r8
    ENDDO
    !-----------------
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(aa1(i) == 0.0_r8)THEN
             ierr(i)=17
          END IF
       END IF
    END DO
    !
    ! determine downdraft strength in terms of windshear
    !
    CALL cup_dd_edt( &
         ierr      , &! ierr    (in)
         us        , &! us      (in)
         vs        , &! vs      (in)
         zo        , &! zo      (in)
         ktop      , &! ktop    (in)
         kbcon     , &! kbcon   (in)
         edt       , &! edt     (out)
         po        , &! po      (in)
         pwavo     , &! pwavo   (in)
         pwevo     , &! pwevo   (in)
         nCols     , &! nCols   (in)
         kMax      , &! kMax    (in)
         istart    , &! istart  (in)
         iend      , &! iend    (in)
         edtmax    , &! edtmax  (in)
         edtmin    , &! edtmin  (in)
         maxens2   , &! maxens2 (in)
         edtc      , &! edtc    (out)
         vshear    , &! vshear  (out)
         sdp       , &! sdp     (out)
         vws       , &! vws     (out)
         mask      , &! mask    (in)
         edtmax1   , &! edtmax1 (in)
         maxens22    )! maxens22(in)

    !nilo new
    !------------------------------------------
    out1_mass_flux_up=0.0_r8
    out2_mass_flux_d=0.0_r8
    out3_detr=0.0_r8
    out4_liq_water=0.0_r8
    out5_qc=0.0_r8
    out6_wf=0.0_r8
    out7_entr=0.0_r8
    out8_he=0.0_r8
    out9_hes=0.0_r8
    out10_hc=0.0_r8
    out11_rh=0.0_r8

    DO i=istart,iend
       IF(ierr(i) ==0)THEN
          out6_wf(i)=aa1(i)-aa0(i)  ! aa1 (at t+1) and aao(at time t)
          DO k=1,kMax
             out1_mass_flux_up(i,k)=zuo(i,k)
             out2_mass_flux_d(i,k)=zdo(i,k)*edt(i)
             out3_detr(i,k)=cd(i,k)
             out4_liq_water(i,k)=qrc(i,k) 
             out5_qc(i,k)=qc(i,k)
             out7_entr(i,k)=entr2D(i,k)
             out8_he(i,k)=heo(i,k)
             out9_hes(i,k)=heso(i,k)
             out10_hc(i,k)=hco(i,k)
             out11_rh(i,k)=rh(i,k)
          ENDDO
       ENDIF
    ENDDO
    !-------------------------------------------

    DO i = 1, ncols
       listim(i)=i
    END DO
    nCols_gz=0
    DO i=1,ncols
       IF ((ierr(i) == 0 ) .OR. (ierr(i) > 995)) THEN
          nCols_gz=nCols_gz+1
          bitx(i)=.TRUE.
       END IF
    END DO

    IF (nCols_gz > 0 ) THEN

       nCols_gz=0
       DO  i=1,ncols
          IF(bitx(i))THEN
             nCols_gz=nCols_gz+1
             litx(nCols_gz) = listim(i)
          END IF
       END DO

       DO k=1,maxens2
          DO i=1,nCols_gz
             edtc_gz (i,k) = edtc  (litx(i),k)
          END DO
       END DO

       DO k=1,ensdim
          DO i=1,nCols_gz
             pr_ens_gz   (i,k) =pr_ens   (litx(i),k) 
             outt_ens_gz (i,k) =outt_ens (litx(i),k)
             xf_ens_gz   (i,k) =xf_ens   (litx(i),k) 
             massfln_gz  (i,k) =massfln  (litx(i),k) 
          END DO
       END DO


       DO j=1,kMax
          DO i=1,nCols_gz
             heo_cup_gz  (i,j)=heo_cup    (litx(i),j)
             zo_cup_gz   (i,j)=zo_cup     (litx(i),j)
             po_cup_gz   (i,j)=po_cup     (litx(i),j)
             hcdo_gz     (i,j)=hcdo       (litx(i),j)
             zdo_gz      (i,j)=zdo        (litx(i),j)
             cdd_gz      (i,j)=cdd        (litx(i),j)
             heo_gz      (i,j)=heo        (litx(i),j)
             qo_cup_gz   (i,j)=qo_cup     (litx(i),j)
             qrcdo_gz    (i,j)=qrcdo      (litx(i),j)
             qo_gz       (i,j)=qo         (litx(i),j)
             zuo_gz      (i,j)=zuo        (litx(i),j)
             cd_gz       (i,j)=cd         (litx(i),j)
             hco_gz      (i,j)=hco        (litx(i),j)
             qco_gz      (i,j)= qco       (litx(i),j)
             qrco_gz     (i,j)= qrco      (litx(i),j)
             tn_gz       (i,j)= tn        (litx(i),j)
             po_gz       (i,j)= po        (litx(i),j)
             gamma_cup_gz(i,j)= gamma_cup (litx(i),j)
             pwo_gz      (i,j)= pwo       (litx(i),j)
             pwdo_gz     (i,j)= pwdo      (litx(i),j)
             he_cup_gz   (i,j)= he_cup    (litx(i),j)
             heso_cup_gz (i,j)= heso_cup  (litx(i),j)
             omeg_gz     (i,j)= omeg      (litx(i),j)
             p_cup_gz    (i,j)= p_cup     (litx(i),j)
             entr2D_gz   (i,j)= entr2D    (litx(i),j)
             detr2D_gz   (i,j)= detr2D    (litx(i),j)
             den1_gz     (i,j)=den1       (litx(i),j)

          END DO
       END DO


       DO i=1,nCols_gz
          ierr_gz    (i)=ierr    (litx(i))
          ktop_gz    (i)=ktop    (litx(i))
          k22_gz     (i)=k22     (litx(i))
          kbcon_gz   (i)=kbcon   (litx(i))
          jmin_gz    (i)=jmin    (litx(i))
          kdet_gz    (i)=kdet    (litx(i))
          z1_gz      (i)=z1      (litx(i))
          psur_gz    (i)=psur    (litx(i))
          kbmax_gz   (i)=kbmax   (litx(i))
          cap_max_gz (i)=cap_max (litx(i))
          aa0_gz     (i)=aa0     (litx(i))
          aa1_gz     (i)=aa1     (litx(i))
          xmb_gz     (i)=xmb     (litx(i))
          mask_gz    (i)=mask    (litx(i))
          mconv_gz   (i)=mconv   (litx(i))
          dcape_zhang_gz    (i)=dcape_zhang    (litx(i))
          k_zhang_gz (i)=k_zhang               (litx(i))
          dcape_gz   (i)=dcape                 (litx(i))
          d_ucape_gz    (i)=d_ucape    (litx(i))
          d_dcape_gz    (i)=d_dcape    (litx(i))
          cine_gz    (i)=cine    (litx(i))
          sens_gz    (i)=cine    (litx(i))
          tkeMax_gz  (i)=tkeMax  (litx(i))
          tkeMedia_gz(i)=tkeMedia(litx(i))
       END DO

       CALL Ensemble(&
            istart                        , &     !INTEGER,(IN   )
            nCols_gz                      , &     !INTEGER,(IN   )
            nCols                         , &     !INTEGER,(IN   )
            kMax                          , &     !INTEGER,(IN   )
            maxens                        , &     !INTEGER,(IN   )
            maxens2                       , &     !INTEGER,(IN   )
            maxens22                      , &     !INTEGER,(IN   )
            maxens3                       , &     !INTEGER,(IN   )
            ensdim                        , &     !INTEGER,(IN   )
            mbdt                          , &     !REAL   ,(IN   )
            dtime                         , &     !REAL   ,(IN   )
            edtc_gz       , &          !REAL   ,(IN   )(nCols,     maxens2)        edtc            (1:nCols   ,1:maxens2        ), &         !REAL   ,(IN        )(nCols,     maxens2)
            ierr_gz       , &          !INTEGER,(INOUT)(nCols             )        ierr            (1:nCols                        ), &         !INTEGER,(INOUT)(nCols             )
            dellat_ens_gz , &          !REAL   ,(OUT  )(nCols,kMax,maxens2)     dellat_ens         (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            dellaq_ens_gz , &          !REAL   ,(OUT  )(nCols,kMax,maxens2)             dellaq_ens         (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            dellaqc_ens_gz, &          !REAL   ,(OUT  )(nCols,kMax,maxens2)             dellaqc_ens   (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            pwo_ens_gz    , &          !REAL   ,(OUT  )(nCols,kMax,maxens2)             pwo_ens         (1:nCols   ,1:kMax,1:maxens2), &     !REAL   ,(OUT  )(nCols,kMax,maxens2)
            heo_cup_gz    , &          !REAL   ,(IN   )(nCols,kMax             )             heo_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            zo_cup_gz     , &          !REAL   ,(IN   )(nCols,kMax             )             zo_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            po_cup_gz     , &          !REAL   ,(INOUT)(nCols,kMax             )             po_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(INOUT)(nCols,kMax         )
            hcdo_gz       , &          !REAL   ,(IN   )(nCols,kMax             )             hcdo          (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            zdo_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             zdo           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            cdd_gz        , &          !REAL   ,(INOUT)(nCols,kMax             )             cdd           (1:nCols   ,1:kMax  ), &     !REAL   ,(INOUT)(nCols,kMax         )
            heo_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             heo           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qo_cup_gz     , &          !REAL   ,(IN   )(nCols,kMax             )             qo_cup         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qrcdo_gz      , &          !REAL   ,(IN   )(nCols,kMax             )             qrcdo         (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qo_gz         , &          !REAL   ,(IN   )(nCols,kMax             )             qo                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            zuo_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             zuo           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            cd_gz         , &          !REAL   ,(IN   )(nCols,kMax             )             cd                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            hco_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             hco           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            ktop_gz       , &          !INTEGER,(IN   )(nCols             )        ktop            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            k22_gz        , &     !INTEGER,(IN   )(nCols              )       k22            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            kbcon_gz      , &          !INTEGER,(IN   )(nCols             )        kbcon            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            jmin_gz       , &          !INTEGER,(IN   )(nCols             )        jmin            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            kdet_gz       , &          !INTEGER,(IN   )(nCols             )        kdet            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            qco_gz        , &          !REAL   ,(IN   )(nCols,kMax             )             qco           (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            qrco_gz       , &          !REAL   ,(IN   )(nCols,kMax             )             qrco          (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            tn_gz         , &          !REAL   ,(IN   )(1:nCols_gz,kMax   )             tn                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            po_gz         , &          !REAL   ,(IN   )(nCols,kMax             )             po                 (1:nCols   ,1:kMax  ), &     !REAL   ,(IN   )(nCols,kMax         )
            z1_gz         , &          !REAL   ,(IN   )(nCols             )        z1            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            psur_gz       , &          !REAL   ,(IN   )(nCols             )        psur            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            gamma_cup_gz  , &            !REAL   ,(INOUT)(nCols,kMax        )             gamma_cup     (1:nCols   ,1:kMax  ), &        !REAL        ,(INOUT)(nCols,kMax           )
            pr_ens_gz     , &            !REAL   ,(INOUT)(nCols,ensdim      )             pr_ens           (1:nCols   ,1:ensdim), &        !REAL        ,(INOUT)(nCols,ensdim           )
            pwo_gz        , &            !REAL   ,(IN   )(nCols,kMax        )             pwo           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            pwdo_gz       , &            !REAL   ,(IN   )(nCols,kMax        )             pwdo           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            outt_ens_gz   , &            !REAL   ,(OUT  )(nCols,ensdim      )             outt_ens           (1:nCols   ,1:ensdim), &        !REAL        ,(OUT  )(nCols,ensdim           )        
            he_cup_gz     , &            !REAL   ,(IN   )(nCols,kMax        )             he_cup           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            kbmax_gz      , &          !INTEGER,(IN   )(nCols             )        kbmax            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            heso_cup_gz   , &          !REAL   ,(IN   )(nCols,kMax             )        heso_cup      (1:nCols   ,1:kMax  ), &         !REAL   ,(IN        )(nCols,kMax            )
            cap_max_gz    , &          !REAL   ,(IN   )(nCols             )        cap_max            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            aa0_gz        , &          !REAL   ,(INOUT)(nCols             )        aa0            (1:nCols                ), &         !REAL   ,(INOUT)(nCols             )
            aa1_gz        , &          !REAL   ,(IN   )(nCols             )        aa1            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            xmb_gz        , &          !REAL   ,(OUT  )(nCols             )        xmb            (1:nCols                ), &         !REAL   ,(OUT  )(nCols             )
            xf_ens_gz     , &          !REAL   ,(OUT  )(nCols,ensdim      )        xf_ens            (1:nCols   ,1:ensdim), &         !REAL   ,(OUT  )(nCols,ensdim      )
            mask_gz       , &          !INTEGER,(IN   )(nCols             )        mask            (1:nCols                ), &         !INTEGER,(IN        )(nCols             )
            mconv_gz      , &          !REAL   ,(IN   )(nCols             )        mconv            (1:nCols                ), &         !REAL   ,(IN        )(nCols             )
            omeg_gz       , &            !REAL   ,(IN   )(nCols,kMax        )             omeg           (1:nCols   ,1:kMax  ), &        !REAL        ,(IN   )(nCols,kMax           )
            massfln_gz    , &            !REAL   ,(OUT  )(nCols,ensdim      )             massfln           (1:nCols   ,1:ensdim), &        !REAL        ,(OUT  )(nCols,ensdim           )
            p_cup_gz      , & !!!nilo )  !REAL   ,(IN   )(nCols,kMax        )             p_cup           (1:nCols   ,1:kMax)  )        !REAL        ,(IN   )(nCols,kMax           )
            entr2D_gz     , &    !
            detr2D_gz     , &  !!)
            den1_gz       , &
            dcape_zhang_gz       , &
            k_zhang_gz        , &
            dcape_gz          , &
            d_ucape_gz       , &
            d_dcape_gz       , &
            cine_gz       , &
            sens_gz       , &
            tkeMax_gz     , &
            tkeMedia_gz   )


       DO i=1,nCols_gz
          ierr(litx(i)) = ierr_gz(i)
       END DO

       DO i=1,nCols_gz
          aa0  (litx(i)) = aa0_gz (i)
          xmb  (litx(i)) = xmb_gz( i)
       END DO

       DO k=1,kMax
          DO i=1,nCols_gz
             po_cup   (litx(i),k) = po_cup_gz   (i,k)
             cdd      (litx(i),k) = cdd_gz      (i,k) 
             gamma_cup(litx(i),k) = gamma_cup_gz(i,k)
             qrco     (litx(i),k) = qrco_gz     (i,k) 
          END DO
       END DO

       DO j=1,ensdim
          DO i=1,nCols_gz
             pr_ens   (litx(i),j) = pr_ens_gz   (i,j) 
             outt_ens (litx(i),j) = outt_ens_gz (i,j) 
             xf_ens   (litx(i),j) = xf_ens_gz   (i,j) 
             massfln  (litx(i),j) = massfln_gz  (i,j) 
          END DO
       END DO
       DO j=1,maxens2
          DO k=1,kMax
             DO i=1,nCols_gz
                dellat_ens  (litx(i),k,j)=dellat_ens_gz (i,k,j)
                dellaq_ens  (litx(i),k,j)=dellaq_ens_gz (i,k,j)
                dellaqc_ens (litx(i),k,j)=dellaqc_ens_gz(i,k,j)
                pwo_ens     (litx(i),k,j)=pwo_ens_gz    (i,k,j)
             END DO
          END DO
       END DO

       !
       !--- FEEDBACK
       !
       CALL cup_output_ens( &
            xf_ens     , & ! xf_ens      (in)
            ierr       , & ! ierr        (inout)
            dellat_ens , & ! dellat_ens  (in)
            dellaq_ens , & ! dellaq_ens  (in)
            dellaqc_ens, & ! dellaqc_ens (in)
            outt       , & ! outt        (inout) hmjb
            outq       , & ! outq        (inout) hmjb
            outqc      , & ! outqc       (out)
            pre        , & ! pre         (out)
            pwo_ens    , & ! pwo_ens     (in)
            xmb        , & ! xmb         (out)
            ktop       , & ! ktop        (in)
            nCols      , & ! nCols       (in)
            kMax       , & ! kMax        (in)
            istart     , & ! istart      (in)
            iend       , & ! iend        (in)
            maxens2    , & ! maxens2     (in)
            maxens     , & ! maxens      (in)
            iens       , & ! iens        (in)
            pr_ens     , & ! pr_ens      (inout)
            outt_ens   , & ! outt_ens    (inout)
            maxens3    , & ! maxens3     (in)
            ensdim     , & ! ensdim      (in)
            massfln    , & ! massfln     (inout)
            xfac1      , & ! xfac1       (out)
            xfac_for_dn, & ! xfac_for_dn (out) 
            maxens22    ) ! maxens22    (in)
       DO i=istart,iend
          pre(i)=MAX(pre(i),0.0_r8)
          !snf
       END DO
    ELSE
       outt =0.0_r8
       outq  =0.0_r8
       outqc =0.0_r8
       pre=0.0_r8
       xmb=0.0_r8
       pr_ens=0.0_r8
       outt_ens=0.0_r8
       massfln=0.0_r8
       xfac1=0.0_r8
       xfac_for_dn=0.0_r8
    END IF

    RETURN
  END SUBROUTINE GDM2
  !
  !END CUP
  !

  SUBROUTINE Ensemble(&
       istart     , &
       iend       , &
       nCols      , &
       kMax       , &
       maxens     , & 
       maxens2    , &
       maxens22   , &
       maxens3    , &
       ensdim     , &
       mbdt       , &
       dtime      , & 
       edtc       , &
       ierr       , &
       dellat_ens , &
       dellaq_ens , &
       dellaqc_ens, &
       pwo_ens    , &
       heo_cup    , &
       zo_cup     , &
       po_cup     , & 
       hcdo       , &
       zdo        , &
       cdd        , &
       heo        , &
       qo_cup     , &
       qrcdo      , &
       qo         , &
       zuo        , &
       cd         , &
       hco        , &
       ktop       , &
       k22        , &
       kbcon      , &
       jmin       , &
       kdet       , &
       qco        , &
       qrco       , &
       tn         , &
       po         , & 
       z1         , & 
       psur       , & 
       gamma_cup  , & 
       pr_ens     , &
       pwo        , &
       pwdo       , &
       outt_ens   , &
       he_cup     , &
       kbmax      , &
       heso_cup   , &    
       cap_max    , &
       aa0        , &
       aa1        , &
       xmb        , &
       xf_ens     , &
       mask       , &
       mconv      , &
       omeg       , &
       massfln    , &
       p_cup      , &  !!)
       entr2D     , &
       detr2D     , & !! )
       den1       , &
       dcape_zhang   , &
       k_zhang       , &
       dcape         , &
       d_ucape    , &
       d_dcape    , &
       cine       , &
       sens       , &
       tkeMax     , &
       tkeMedia   )


    INTEGER :: iedt
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: maxens
    INTEGER, INTENT(IN   ) :: maxens2,maxens22
    INTEGER, INTENT(IN   ) :: maxens3
    INTEGER, INTENT(IN   ) :: ensdim 
    REAL(KIND=r8)   , INTENT(IN   ) :: mbdt
    REAL(KIND=r8)   , INTENT(IN   ) :: dtime
    REAL(KIND=r8)   , INTENT(IN   ) :: edtc       (nCols,maxens2)
    INTEGER, INTENT(INOUT) :: ierr       (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dellat_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dellaq_ens (nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dellaqc_ens(nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwo_ens    (nCols,kMax, maxens2)
    REAL(KIND=r8)   , INTENT(IN   ) :: heo_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zo_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: po_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hcdo       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zdo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: cdd        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: heo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qo_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qrcdo      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qo         (nCols,kMax)! water vapor mixing ratio (kg/kg) at time t+1
    REAL(KIND=r8)   , INTENT(IN   ) :: zuo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cd         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hco        (nCols,kMax)
    INTEGER, INTENT(IN   ) :: ktop       (nCols)
    INTEGER, INTENT(IN   ) :: k22        (nCols)
    INTEGER, INTENT(IN   ) :: kbcon      (nCols)! level of convective cloud base       
    INTEGER, INTENT(IN   ) :: jmin       (nCols)
    INTEGER, INTENT(IN   ) :: kdet       (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: qco        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qrco       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: tn         (nCols,kMax)! temperature (K) at time t+1
    REAL(KIND=r8)   , INTENT(IN   ) :: po         (nCols,kMax)! pressao de superficie no tempo t mb
    REAL(KIND=r8)   , INTENT(IN   ) :: z1         (nCols)! topography (m)
    REAL(KIND=r8)   , INTENT(IN   ) :: psur       (nCols)! pressao de superficie no tempo t mb
    REAL(KIND=r8)   , INTENT(INOUT) :: gamma_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: pr_ens     (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(IN   ) :: pwo        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: pwdo       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: outt_ens   (nCols,ensdim)    
    REAL(KIND=r8)   , INTENT(IN   ) :: he_cup     (nCols,kMax)
    INTEGER, INTENT(IN   ) :: kbmax      (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: heso_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cap_max    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0        (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: aa1        (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: xmb        (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: xf_ens     (nCols,ensdim)
    INTEGER, INTENT(IN   ) :: mask       (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: mconv      (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: omeg       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: massfln    (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(IN   ) :: p_cup      (nCols,kMax)

    !
    ! LOCAL VARIABLE
    !
    REAL(KIND=r8)    :: dellat    (nCols,kMax)
    REAL(KIND=r8)    :: dellaq    (nCols,kMax)
    REAL(KIND=r8)    :: dellah    (nCols,kMax)
    REAL(KIND=r8)    :: dellaqc   (nCols,kMax)
    REAL(KIND=r8)    :: xhe       (nCols,kMax)
    REAL(KIND=r8)    :: edt       (nCols)
    REAL(KIND=r8)    :: bu        (nCols)    
    INTEGER          :: ierr2     (nCols)
    INTEGER          :: ierr3     (nCols)    
    REAL(KIND=r8)    :: dbyd      (nCols,kMax)
    REAL(KIND=r8)    :: xq        (nCols,kMax)
    REAL(KIND=r8)    :: xt        (nCols,kMax)
    REAL(KIND=r8)    :: xqes      (nCols,kMax)
    REAL(KIND=r8)    :: edto      (nCols)
    REAL(KIND=r8)    :: xhes      (nCols,kMax)
    REAL(KIND=r8)    :: xff_ens3  (nCols,maxens3)
    REAL(KIND=r8)    :: xk        (nCols,maxens) 
    REAL(KIND=r8)    :: xaa0_ens  (nCols,maxens)
    INTEGER          :: k22x      (nCols)    
    INTEGER          :: kbconx    (nCols)
    INTEGER          :: nallp     
    REAL(KIND=r8)    :: xaa0      (nCols)
    REAL(KIND=r8)    :: denx      (nCols,kMax)
    REAL(KIND=r8)    :: xt_cup    (nCols,kMax)
    REAL(KIND=r8)    :: xdby      (nCols,kMax)
    REAL(KIND=r8)    :: xzu       (nCols,kMax)
    REAL(KIND=r8)    :: xz        (nCols,kMax)
    REAL(KIND=r8)    :: xq_cup    (nCols,kMax)
    REAL(KIND=r8)    :: xqes_cup  (nCols,kMax)
    REAL(KIND=r8)    :: xpwav     (nCols)
    REAL(KIND=r8)    :: xpw       (nCols,kMax)
    REAL(KIND=r8)    :: xqrc      (nCols,kMax)
    REAL(KIND=r8)    :: xqc       (nCols,kMax)
    REAL(KIND=r8)    :: xz_cup    (nCols,kMax)
    REAL(KIND=r8)    :: xqrcd     (nCols,kMax)
    REAL(KIND=r8)    :: xpwev     (nCols)
    REAL(KIND=r8)    :: xpwd      (nCols,kMax)
    REAL(KIND=r8)    :: xqcd      (nCols,kMax)
    REAL(KIND=r8)    :: xhes_cup  (nCols,kMax)
    REAL(KIND=r8)    :: xhcd      (nCols,kMax)
    REAL(KIND=r8)    :: xzd       (nCols,kMax)
    REAL(KIND=r8)    :: xhc       (nCols,kMax)
    REAL(KIND=r8)    :: xhe_cup   (nCols,kMax)
    REAL(KIND=r8)    :: xhkb      (nCols)
    REAL(KIND=r8)    :: scr1      (nCols,kMax)
    REAL(KIND=r8)    :: dz
    REAL(KIND=r8)    :: massfld
    INTEGER :: i
    INTEGER :: k
    INTEGER :: nens
    INTEGER :: nens3
    REAL(KIND=r8)   , INTENT(IN   ) :: entr2D   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: den1     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: dcape_zhang (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: k_zhang (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: dcape(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: d_ucape (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: d_dcape (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cine (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: sens (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: tkeMax(nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: tkeMedia(nCols)

    dellat    =0.0_r8;dellaq    =0.0_r8
    dellah    =0.0_r8;dellaqc   =0.0_r8
    xhe       =0.0_r8;edt       =0.0_r8
    bu        =0.0_r8;ierr2     =0
    ierr3     =0;dbyd           =0.0_r8
    xq        =0.0_r8;xt        =0.0_r8
    xqes      =0.0_r8;edto      =0.0_r8
    xhes      =0.0_r8;xff_ens3  =0.0_r8
    xk        =0.0_r8;xaa0_ens  =0.0_r8
    k22x      =0;kbconx         =0
    nallp     =0;xaa0           =0.0_r8
    xt_cup    =0.0_r8;xdby      =0.0_r8
    xzu       =0.0_r8;xz        =0.0_r8
    xq_cup    =0.0_r8;xqes_cup  =0.0_r8
    xpwav     =0.0_r8;xpw       =0.0_r8
    xqrc      =0.0_r8;xqc       =0.0_r8
    xz_cup    =0.0_r8;xqrcd     =0.0_r8
    xpwev     =0.0_r8;xpwd      =0.0_r8
    xqcd      =0.0_r8;xhes_cup  =0.0_r8
    xhcd      =0.0_r8;xzd       =0.0_r8
    xhc       =0.0_r8;xhe_cup   =0.0_r8
    xhkb      =0.0_r8;scr1      =0.0_r8
    dz        =0.0_r8;massfld   =0.0_r8
    !
    ! LOOP FOR ENSEMBLE MAXENS2
    !
    DO 250 iedt=1,maxens22

       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             edt (i)=edtc(i,iedt)
             edto(i)=edtc(i,iedt)
          END IF
       END DO

       DO k=1,kMax
          DO i=istart,iend
             dellat_ens (i,k,iedt)=0.0_r8
             dellaq_ens (i,k,iedt)=0.0_r8
             dellaqc_ens(i,k,iedt)=0.0_r8
             pwo_ens    (i,k,iedt)=0.0_r8
          END DO
       END DO
       !
       !--- downdraft workfunctions
       !
       !
       !--- change per unit mass that a model cloud would modify the environment
       !
       !--- 1.0_r8 in bottom layer
       !
       CALL cup_dellabot( &
            heo_cup    , &  ! heo_cup     (in)
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            hcdo       , &  ! hcdo        (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            heo        , &  ! heo         (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellah     , &  ! dellah      (out)
            mentrd_rate, &  !  mentrd_rate (in)
            detr2D)

       CALL cup_dellabot(&
            qo_cup     , &  ! qo_cup      (in)
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            qrcdo      , &  ! qrcdo       (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            qo         , &  ! qo          (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellaq     , &  ! dellaq      (out)
            mentrd_rate, &  !  mentrd_rate (in)
            detr2D)
       !
       !--- 2. everywhere else
       !
       CALL cup_dellas(&
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            hcdo       , &  ! hcdo        (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            heo        , &  ! heo         (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellah     , &  ! dellah      (out)
            mentrd_rate, &  ! mentrd_rate (in)
            zuo        , &  ! zuo         (in)
            cd         , &  ! cd          (in)
            hco        , &  ! hco         (in)
            ktop       , &  ! ktop        (in)
            k22        , &  ! k22         (in)
            kbcon      , &  ! kbcon       (in)
            mentr_rate , &  ! mentr_rate  (in)
            jmin       , &  ! jmin        (in)
            heo_cup    , &  ! heo_cup     (in)
            kdet       , &  ! kdet        (in)
            k22        , &  !)    ! k22         (in)
            entr2D     , &  !  )  ! mentrd_rate (in)
            detr2D      )


       !-- take out cloud liquid water for detrainment
       DO k=1,kMax
          DO i=istart,iend
             scr1   (i,k)=0.0_r8
             dellaqc(i,k)=0.0_r8
             IF(ierr(i) == 0)THEN
                scr1(i,k)=qco(i,k)-qrco(i,k)
                IF(k == ktop(i)-0)dellaqc(i,k)=                 &
                     0.01_r8*zuo(i,ktop(i))*qrco(i,ktop(i))*        &
                     9.81_r8/(po_cup(i,k  )-po_cup(i,k+1))

                IF(k <  ktop(i)  .AND.k >  kbcon(i))THEN
                   dz=zo_cup(i,k+1)-zo_cup(i,k)
                   dellaqc(i,k)=0.01_r8*9.81_r8*cd(i,k)*dz*zuo(i,k)    &
                        *0.5_r8*(qrco(i,k)+qrco(i,k+1))/            &
                        (po_cup(i,k  )-po_cup(i,k+1))
                END IF
             END IF
          END DO
       END DO

       !
       CALL cup_dellas( &
            ierr       , &  ! ierr        (in)
            zo_cup     , &  ! zo_cup      (in)
            po_cup     , &  ! po_cup      (in)
            qrcdo      , &  ! qrcdo       (in)
            edto       , &  ! edto        (in)
            zdo        , &  ! zdo         (in)
            cdd        , &  ! cdd         (in)
            qo         , &  ! qo          (in)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            dellaq     , &  ! dellaq      (out)
            mentrd_rate, &  ! mentrd_rate (in)
            zuo        , &  ! zuo         (in)
            cd         , &  ! cd          (in)
            scr1       , &  ! scr1        (in)
            ktop       , &  ! ktop        (in)
            k22        , &  ! k22         (in)
            kbcon      , &  ! kbcon       (in)
            mentr_rate , &  ! mentr_rate  (in)
            jmin       , &  ! jmin        (in)
            qo_cup     , &  ! qo_cup      (in)
            kdet       , &  ! kdet        (in)
            k22        , &  !  )  ! k22         (in)
            entr2D     , & !  )  ! mentrd_rate (in)
            detr2D)
       !
       !--- using dellas, calculate changed environmental profiles
       !
       !................second loop..........................start 200
!!!old       do 200 nens=1,maxens
!!!old           mbdt=mbdt_ens(nens)
       !
       DO k=1,maxens
          DO i=istart,iend
             xaa0_ens(i,k)=0.0_r8
          END DO
       END DO
       !
       !-----------------------------
       !
       DO k=1,kMax-1
          DO i=istart,iend
             dellat(i,k)=0.0_r8
             IF(ierr(i) == 0)THEN
                xhe   (i,k)=dellah(i,k)*mbdt+heo(i,k)
                xq    (i,k)=dellaq(i,k)*mbdt+qo (i,k)
                dellat(i,k)=(1.0_r8/1004.0_r8)*(dellah(i,k)-2.5e06_r8*dellaq(i,k))
                xt    (i,k)= dellat(i,k)*mbdt+tn(i,k)
                IF(xq(i,k) <= 0.0_r8)xq(i,k)=1.e-08_r8
             END IF
          END DO
       END DO
       !
       !
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             xhe(i,kMax)=heo(i,kMax)
             xq (i,kMax)=qo (i,kMax)
             xt (i,kMax)=tn (i,kMax)
             IF(xq(i,kMax) <= 0.0_r8)xq(i,kMax)=1.e-08_r8
          END IF
       END DO
       !
       ! calculate moist static energy, heights, qes
       !
       CALL cup_env(&
            xz        , &  ! xz     (out)
            xqes      , &  ! xqes   (out)
            xhe       , &  ! xhe    (inout)
            xhes      , &  ! xhes   (out)
            xt        , &  ! xt     (in)
            xq        , &  ! xq     (inout)
            po        , &  ! po     (in)
            z1        , &  ! z1     (in)
            nCols     , &  ! nCols  (in)
            kMax      , &  ! kMax   (in)
            istart    , &  ! istart (in)
            iend      , &  ! iend   (in)
            psur      , &  ! psur   (in)
            ierr      , &  ! ierr   (in)
            2         , &  !  2     (in)
            denx        )  ! den    (out)
       !
       ! environmental values on cloud levels
       !
       CALL cup_env_clev( &
            xt        , &  ! xt        (in)
            xqes      , &  ! xqes      (in)
            xq        , &  ! xq        (in)
            xhe       , &  ! xhe       (in)
            xhes      , &  ! xhes      (in)
            xz        , &  ! xz        (in)
            po        , &  ! po        (in)
            xqes_cup  , &  ! xqes_cup  (out)
            xq_cup    , &  ! xq_cup    (out)
            xhe_cup   , &  ! xhe_cup   (out)
            xhes_cup  , &  ! xhes_cup  (out)
            xz_cup    , &  ! xz_cup    (out)
            po_cup    , &  ! po_cup    (out)
            gamma_cup , &  ! gamma_cup (out)
            xt_cup    , &  ! xt_cup    (out)
            psur      , &  ! psur      (in)
            nCols     , &  ! nCols     (in)
            kMax      , &  ! kMax      (in)
            istart    , &  ! istart    (in)
            iend      , &  ! iend      (in)
            ierr      , &  ! ierr      (in)
            z1          )  ! z1        (in)
       !
       !STATIC CONTROL
       !
       ! moist static energy inside cloud
       !
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             xhkb(i)=xhe(i,k22(i))
          END IF
       END DO

       CALL cup_up_he(&
            k22       , & ! k22        (in)
            xhkb      , & ! xhkb       (out)
            xz_cup    , & ! xz_cup     (in)
            cd        , & ! cd         (in)
            mentr_rate, & ! mentr_rate (in)
            xhe_cup   , & ! xhe_cup    (in)
            xhc       , & ! xhc        (out)
            nCols     , & ! nCols      (in)
            kMax      , & ! kMax       (in)
            kbcon     , & ! kbcon      (in)
            ierr      , & ! ierr       (in)
            istart    , & ! istart     (in)
            iend      , & ! iend       (in)
            xdby      , & ! xdby       (out)
            xhe       , & ! xhe        (in)
            xhes_cup  , & !  ) ! xhes_cup   (in)
            entr2D)
       !
       ! normalized mass flux profile
       !
       CALL cup_up_nms(&
            xzu       , & ! xzu        (out)
            xz_cup    , & ! xz_cup     (in)
            mentr_rate, & ! mentr_rate (in)
            cd        , & ! cd         (in)
            kbcon     , & ! kbcon      (in)
            ktop      , & ! ktop       (in)
            nCols     , & ! nCols      (in)
            kMax      , & ! kMax       (in)
            istart    , & ! istart     (in)
            iend      , & ! iend       (in)
            ierr      , & ! ierr       (in)
            k22       , & ! ) ! k22        (in)
            entr2D)

       CALL cup_dd_nms(xzd        , &! xzd        (out)
            xz_cup     , &! xz_cup     (in)
            cdd        , &! cdd        (out)
            mentrd_rate, &! mentrd_rate(in)
            jmin       , &! jmin       (in)
            ierr       , &! ierr       (in)
            nCols      , &! nCols      (in)
            kMax       , &! kMax       (in)
            istart     , &! istart     (in)
            iend       , &! iend       (in)
            1          , &! 1
            kdet       , &! kdet       (in)
            z1         , & !  )! z1         (in)
            detr2D)

       !
       ! moisture downdraft
       !
       CALL cup_dd_he(xhes_cup   , &! xhes_cup    (in)
            xhcd       , &! xhcd        (out)
            xz_cup     , &! xz_cup      (in)
            cdd        , &! cdd         (in)
            mentrd_rate, &! mentrd_rate (in)
            jmin       , &! jmin        (in)
            ierr       , &! ierr        (in)
            nCols      , &! nCols       (in)
            kMax       , &! kMax        (in)
            istart     , &! istart      (in)
            iend       , &! iend        (in)
            xhe        , &! xhe         (in)
            dbyd       , & !  )! dbyd        (out)
            detr2D)

       CALL cup_dd_moisture(xzd        , &  ! xzd         (in)
            xhcd       , &  ! xhcd        (in)
            xhes_cup   , &  ! xhes_cup    (in)
            xqcd       , &  ! xqcd        (out)
            xqes_cup   , &  ! xqes_cup    (in)
            xpwd       , &  ! xpwd        (out)
            xq_cup     , &  ! xq_cup      (in)
            xz_cup     , &  ! xz_cup      (in)
            cdd        , &  ! cdd         (in)
            mentrd_rate, &  ! mentrd_rate (in)
            jmin       , &  ! jmin        (in)
            ierr       , &  ! ierr        (inout)
            gamma_cup  , &  ! gamma_cup   (in)
            xpwev      , &  ! xpwev       (out)
            nCols      , &  ! nCols       (in)
            kMax       , &  ! kMax        (in)
            istart     , &  ! istart      (in)
            iend       , &  ! iend        (in)
            bu         , &  ! bu          (out)
            xqrcd      , &  ! xqrcd       (out)
            xq         , &  ! xq          (in)
            3          , & !  )  ! 3
            detr2D)

       !
       ! moisture updraft
       !

       CALL cup_up_moisture(ierr       , &  ! ierr       (in)
            xz_cup     , &  ! xz_cup     (in)
            xqc        , &  ! xqc        (out)
            xqrc       , &  ! xqrc       (out)
            xpw        , &  ! xpw        (out)
            xpwav      , &  ! xpwav      (out)
            kbcon      , &  ! kbcon      (in)
            ktop       , &  ! ktop       (in)
            nCols      , &  ! nCols      (in)
            kMax       , &  ! kMax       (in)
            istart     , &  ! istart     (in)
            iend       , &  ! iend       (in)
            cd         , &  ! cd         (in)
            xdby       , &  ! xdby       (inout)
            mentr_rate , &  ! mentr_rate (in)
            xq         , &  ! xq         (in)
            gamma_cup  , &  ! gamma_cup  (in)
            xzu        , &  ! xzu        (in)
            xqes_cup   , &  ! xqes_cup   (in)
            k22        , &  ! k22        (in)
            xq_cup     , &  !  )  ! xq_cup     (in)
            detr2D)

       !
       ! workfunctions for updraft
       !
       CALL cup_up_aa0(xaa0       , & ! xaa0      (inout)
            xz         , & ! xz        (in)
            xzu        , & ! xzu       (in)
            xdby       , & ! xdby      (in)
            gamma_cup  , & ! gamma_cup (in)
            xt_cup     , & ! xt_cup    (in)
            kbcon      , & ! kbcon     (in)
            ktop       , & ! ktop      (in)
            kMax       , & ! kMax      (in)
            nCols      , & ! nCols     (in)
            istart     , & ! istart    (in)
            iend       , & ! iend      (in)
            ierr         ) ! ierr      (in)

       !
       ! workfunctions for downdraft
       !---------0--------------
       ! 
       DO 200 nens=1,maxens
          DO i=istart,iend 
             IF(ierr(i) == 0)THEN
                xaa0_ens(i,nens)=xaa0(i)
             END IF
          END DO
          nallp=(iens-1)*maxens3*maxens*maxens22 &
               +(iedt-1)*maxens*maxens3 &
               +(nens-1)*maxens3
          DO nens3=1,maxens3
             DO k=1,kMax
                DO i=istart,iend
                   IF( k <= ktop(i) .AND. ierr(i) == 0 )THEN                
                      pr_ens(i,nallp+nens3)=pr_ens(i,nallp+nens3)+&
                           pwo(i,k)+1.0_r8*edto(i)*pwdo(i,k)
                   END IF
                END DO
             END DO
             DO i=istart,iend 
                IF(ierr(i) == 0)THEN
                   outt_ens (i,nallp+nens3)=dellat(i,1)
                   IF(pr_ens(i,nallp+nens3) < 0.0_r8)THEN
                      pr_ens(i,nallp+nens3)= 0.0_r8
                   END IF
                END IF
             END DO
          END DO
200    END DO
       !...............end 200
       !
       ! LARGE SCALE FORCING
       !
       CALL cup_maximi(he_cup    , & ! he_cup (in)
            nCols     , & ! nCols  (in)
            kMax      , & ! kMax   (in) 
            3         , & ! 3      (in)
            kbmax     , & ! kbmax  (in)
            k22x      , & ! k22x   (out)
            istart    , & ! istart (in)
            iend      , & ! iend   (in)
            ierr        ) ! ierr   (in)

       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             k22x (i)=k22(i)
          END IF
          ierr2(i)=ierr(i)
          ierr3(i)=ierr(i)
       END DO
       !
       ! --- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
       ! snf  call cup_kbcon for cap_max=cap_max-(2-1)*cap_max_increment
       !
       CALL cup_kbcon(&
            cap_max_increment, &  ! cap_max_increment (in)
            2         , &         ! 2
            k22x      , &         ! k22x              (inout)
            kbconx    , &         ! kbconx            (out)
            heo_cup   , &         ! heo_cup           (in)
            heso_cup  , &         ! heso_cup          (in)
            nCols     , &         ! nCols             (in)
            kMax      , &         ! kMax              (in)
            istart    , &         ! istart            (in)
            iend      , &         ! iend              (in)
            ierr2     , &         ! ierr2             (inout)
            kbmax     , &         ! kbmax             (in)
            po_cup    , &         ! po_cup            (in)
            cap_max     )         ! cap_max           (in)
       !
       ! snf  call cup_kbcon for cap_max=cap_max-(3-1)*cap_max_increment
       !
       CALL cup_kbcon(&
            cap_max_increment, &  ! cap_max_increment (in)
            3                , &  ! 3                 (in)
            k22x             , &  ! k22x              (inout)
            kbconx           , &  ! kbconx            (out)
            heo_cup          , &  ! heo_cup           (in)
            heso_cup         , &  ! heso_cup          (in)
            nCols            , &  ! nCols             (in)
            kMax             , &  ! kMax              (in)
            istart           , &  ! istart            (in)
            iend             , &  ! iend              (in)
            ierr3            , &  ! ierr3             (inout)
            kbmax            , &  ! kbmax             (in)
            po_cup           , &  ! po_cup            (in)
            cap_max            )  ! cap_max           (in)

       CALL cup_forcing_ens( &
            aa0       , & ! aa0      (inout)
            aa1       , & ! aa1      (in)
            xaa0_ens  , & ! xaa0_ens (in)
            mbdt      , & ! mbdt     (in)
            dtime     , & ! dtime    (in)
            xmb       , & ! xmb      (out)
            ierr      , & ! ierr     (inout)
            nCols     , & ! nCols    (in)
            kMax      , & ! kMax     (in)
            istart    , & ! istart   (in)
            iend      , & ! iend     (in)
            xf_ens    , & ! xf_ens   (out)
            'deeps'   , & ! 'deeps'  (in)
            mask      , & ! mask     (in)
            maxens    , & ! maxens   (in)
            iens      , & ! iens     (in)
            iedt      , & ! iedt     (in)
            maxens3   , & ! maxens3  (in)
            mconv     , & ! mconv    (in)
            omeg      , & ! omeg     (in)
            k22       , & ! k22      (in)
            pr_ens    , & ! pr_ens   (in)
            edto      , & ! edto     (in)
            kbcon     , & ! kbcon    (in)
            ensdim    , & ! ensdim   (in)
            massfln   , & ! massfln  (out)
            massfld   , & ! massfld  (inout)
            xff_ens3  , & ! xff_ens3 (out)
            xk        , & ! xk       (out)
            p_cup     , & ! p_cup    (in)
            ktop      , & ! ktop     (in)
            ierr2     , & ! ierr2    (in)
            ierr3     , & ! ierr3    (in)
            1         , & ! gdmpar  (in)
            maxens22  , & !  ) !maxens22  (in)
            den1      , &
            dcape_zhang  , &
            k_zhang      , &
            dcape        , &
            d_ucape   , &
            d_dcape   , &
            cine      , &
            sens      , &
            tkeMax    , &
            tkeMedia  , &
            zo_cup     )

       DO k=1,kMax
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                dellat_ens (i,k,iedt)=dellat (i,k)
                dellaq_ens (i,k,iedt)=dellaq (i,k)
                dellaqc_ens(i,k,iedt)=dellaqc(i,k)
                pwo_ens    (i,k,iedt)=pwo    (i,k)+edt(i)*pwdo(i,k)
             ELSE 
                dellat_ens (i,k,iedt)=0.0_r8
                dellaq_ens (i,k,iedt)=0.0_r8
                dellaqc_ens(i,k,iedt)=0.0_r8
                pwo_ens    (i,k,iedt)=0.0_r8
             END IF
          END DO
       END DO
250 END DO

  END SUBROUTINE Ensemble
  !
  !END CUP
  !-----------------------------------------------------------------------subroutines
  !*-----------
  SUBROUTINE cup_env( &
       z      ,qes    ,he     ,hes    ,t      ,q      , &
       p      ,z1     ,nCols  ,kMax   ,istart ,iend   ,psur   , &
       ierr   ,itest,den                                    )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: nCols
    INTEGER, INTENT(IN   )    :: kMax
    REAL(KIND=r8)   , INTENT(OUT  )    :: z   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )    :: qes (nCols,kMax)! pressure vapor
    REAL(KIND=r8)   , INTENT(OUT  )    :: den (nCols,kMax)!density
    REAL(KIND=r8)   , INTENT(INOUT)    :: he  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )    :: hes (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )    :: t   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT)    :: q   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )    :: p   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )    :: z1  (nCols)
    INTEGER, INTENT(IN   )    :: istart
    INTEGER, INTENT(IN   )    :: iend
    REAL(KIND=r8)   , INTENT(IN   )    :: psur(nCols)
    INTEGER, INTENT(IN   )    :: ierr(nCols)
    INTEGER, INTENT(IN   )    :: itest
    REAL(KIND=r8)             :: esft
    INTEGER                   :: kase
    !
    ! local variables
    !
    INTEGER                   :: i
    INTEGER                   :: k
    REAL(KIND=r8)                      :: tv  (nCols,kMax) ! virtual temperature
    den=0.0_r8
    kase=2   !nilo2
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             !
             ! sgb - IPH is for phase, dependent on TCRIT (water or ice)
             ! calculation of the pressure vapor
             !
             IF(kase==1)THEN
                esft=es5(t(i,k))
             ELSE
                esft=fpvs2es5(t(i,k))
             ENDIF
             !qes(i,k) = 0.622_r8*esft/(100.0_r8*p(i,k)-esft)
             qes(i,k) = 0.622_r8*esft/MAX((100.0_r8*p(i,k) - esft),1.0e-12_r8)
             IF(qes(i,k) <= 1.0e-08_r8  )      qes(i,k)=1.0e-08_r8
             IF(q(i,k)   >  qes(i,k))        q(i,k)=qes(i,k)
             !
             ! calculation of virtual temperature
             !
             tv(i,k) = t(i,k)+0.608_r8*q(i,k)*t(i,k)
             den(i,k)=100.0_r8*p(i,k)/(287.0_r8*tv(i,k))  !nilo

          END IF
       END DO
    END DO
    !
    ! z's are calculated with changed h's and q's and t's
    ! if itest=2
    !
    !
    ! calculate heights geopotential
    !
    IF(itest.NE.2)THEN
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             z(i,1) = MAX(0.0_r8,z1(i))-(LOG(p(i,1))-LOG(psur(i)) )*287.0_r8 &
                  * tv(i,1)/9.81_r8
             he  (i,1)=9.81_r8*z(i,1)+1004.0_r8*t(i,1)+2.5e06_r8*q  (i,1)
             hes (i,1)=9.81_r8*z(i,1)+1004.0_r8*t(i,1)+2.5e06_r8*qes(i,1)
             IF(he(i,1) >= hes(i,1))he(i,1)=hes(i,1)
          END IF
       END DO
       DO k=2,kMax
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                z(i,k) = z(i,k-1)     -(LOG(p(i,k))-LOG(p(i,k-1)))*287.0_r8 &
                     * (0.5_r8*tv(i,k)+0.5_r8*tv(i,k-1) )/9.81_r8
                he  (i,k)=9.81_r8*z(i,k)+1004.0_r8*t(i,k)+2.5e06_r8*q  (i,k)
                hes (i,k)=9.81_r8*z(i,k)+1004.0_r8*t(i,k)+2.5e06_r8*qes(i,k)
                IF(he(i,k) >= hes(i,k))he(i,k)=hes(i,k)
             END IF
          END DO
       END DO
    ELSE
       DO k=1,kMax
          DO i=istart,iend
             IF(ierr(i) == 0)THEN
                z(i,k)=(he(i,k)-1004.0_r8*t(i,k)-2.5e6_r8*q(i,k))/9.81_r8
                z(i,k)=MAX(1.0e-3_r8,z(i,k))
                hes(i,k)=9.81_r8*z(i,k)+1004.0_r8*t(i,k)+2.5e06_r8*qes(i,k)
                IF(he(i,k) >= hes(i,k))he(i,k)=hes(i,k)
             END IF
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE cup_env

  !*--------
  SUBROUTINE cup_env_clev( &
       t        ,qes      ,q      ,he       ,hes     ,z      , &
       p        ,qes_cup  ,q_cup  ,he_cup   ,hes_cup ,z_cup  ,p_cup  , &
       gamma_cup,t_cup    ,psur   ,nCols    ,kMax    ,istart ,iend   , &
       ierr     ,z1                                                 )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )                  :: nCols
    INTEGER, INTENT(IN   )                  :: kMax
    INTEGER, INTENT(IN   )                  :: istart
    INTEGER, INTENT(IN   )                  :: iend
    INTEGER, INTENT(IN   )                  :: ierr(nCols)
    REAL(KIND=r8)   , INTENT(IN   )                  :: t        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: qes      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: q        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: he       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: hes      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: z        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: p        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: qes_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: q_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: he_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: hes_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: z_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: p_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: gamma_cup(nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                  :: t_cup    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                  :: psur     (nCols)    
    REAL(KIND=r8)   , INTENT(IN   )                  :: z1       (nCols)    

    INTEGER                  :: i
    INTEGER                  :: k

    ! ierr error value, maybe modified in this routine
    ! q           = environmental mixing ratio
    ! q_cup       = environmental mixing ratio on cloud levels
    ! qes         = environmental saturation mixing ratio
    ! qes_cup     = environmental saturation mixing ratio on cloud levels
    ! t           = environmental temp
    ! t_cup       = environmental temp on cloud levels
    ! p           = environmental pressure
    ! p_cup       = environmental pressure on cloud levels
    ! z           = environmental heights
    ! z_cup       = environmental heights on cloud levels
    ! he          = environmental moist static energy
    ! he_cup      = environmental moist static energy on cloud levels
    ! hes         = environmental saturation moist static energy
    ! hes_cup     = environmental saturation moist static energy on cloud levels
    ! gamma_cup   = gamma on cloud levels
    ! psur        = surface pressure
    ! z1          = terrain elevation


    DO k=2,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             qes_cup(i,k) = 0.5_r8*(qes(i,k-1) + qes(i,k))
             q_cup  (i,k) = 0.5_r8*(  q(i,k-1) +   q(i,k))
             hes_cup(i,k) = 0.5_r8*(hes(i,k-1) + hes(i,k))
             he_cup (i,k) = 0.5_r8*( he(i,k-1) +  he(i,k))

             IF(he_cup(i,k)  >   hes_cup(i,k)) he_cup(i,k) = hes_cup(i,k)
             z_cup    (i,k) = 0.5_r8*(z(i,k-1) + z(i,k))
             p_cup    (i,k) = 0.5_r8*(p(i,k-1) + p(i,k))
             t_cup    (i,k) = 0.5_r8*(t(i,k-1) + t(i,k))
             gamma_cup(i,k) =(xl/cp)*(xl/(rv*t_cup(i,k)    &
                  *t_cup(i,k)))*qes_cup(i,k)
          END IF
       END DO
    END DO
    !
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          qes_cup  (i,1) =  qes(i,1)
          q_cup    (i,1) =    q(i,1)
          hes_cup  (i,1) =  hes(i,1)
          he_cup   (i,1) =   he(i,1)

          z_cup    (i,1) = 0.5_r8*( z(i,1) +   z1(i))
          p_cup    (i,1) = 0.5_r8*( p(i,1) + psur(i))
          t_cup    (i,1) =      t(i,1)
          gamma_cup(i,1) = xl/cp*(xl/(rv*t_cup(i,1)               &
               *t_cup(i,1)))*qes_cup(i,1)
       END IF
    END DO
  END SUBROUTINE cup_env_clev

  !*--------

  SUBROUTINE cup_maximi( &
       array    ,nCols    ,kMax      ,ks       ,ke       , &
       maxx     ,istart   ,iend     ,ierr)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols 
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ks
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    REAL(KIND=r8)   , INTENT(IN   ) :: array(nCols, kMax)
    INTEGER, INTENT(IN   ) :: ierr (nCols) 
    INTEGER, INTENT(OUT  ) :: maxx (nCols) 
    INTEGER, INTENT(IN   ) :: ke   (nCols)

    REAL(KIND=r8)                   :: x    (nCols)
    INTEGER                :: i
    INTEGER                :: k

    DO i=istart,iend
       maxx(i)=ks
       IF(ierr(i) == 0)THEN
          x(i)=array(i,ks)
       END IF
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= ks .AND.  k <= ke(i) ) THEN
             IF(array(i,k) >= x(i)) THEN
                x(i)=array(i,k)
                maxx(i)=k
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_maximi

  !*------

  SUBROUTINE cup_minimi(array     ,nCols    ,kMax       ,ks        ,kend      , &
       kt        ,istart   ,iend      ,ierr)

    IMPLICIT NONE
    INTEGER, INTENT(IN   )   :: nCols 
    INTEGER, INTENT(IN   )   :: kMax
    INTEGER, INTENT(IN   )   :: istart
    INTEGER, INTENT(IN   )   :: iend
    REAL(KIND=r8)   , INTENT(IN   )   :: array(nCols, kMax)
    INTEGER, INTENT(OUT  )   :: kt   (nCols)
    INTEGER, INTENT(IN   )   :: ks   (nCols)
    INTEGER, INTENT(IN   )   :: kend (nCols)
    INTEGER, INTENT(IN   )   :: ierr (nCols)

    REAL(KIND=r8)                     :: x    (nCols     )
    INTEGER                  :: kstop(nCols     )
    INTEGER                  :: i
    INTEGER                  :: k 

    DO i=istart,iend
       kt(i)=ks(i)
       IF(ierr(i) == 0)THEN
          x    (i)=array(i,ks(i))
          kstop(i)=MAX(ks(i)+1,kend(i))
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             IF (k >= ks(i)+1 .AND. k <= kstop(i)) THEN
                IF(array(i,k) <  x(i)) THEN
                   x(i)=array(i,k)
                   kt(i)=k
                END IF
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_minimi

  SUBROUTINE cup_kbcon(cap_inc   ,&
       iloop     ,k22       ,kbcon     ,he_cup    ,hes_cup   , &
       nCols     ,kMax      ,istart    ,iend      ,ierr      , &
       kbmax     ,p_cup     ,cap_max)
    IMPLICIT NONE
    INTEGER, INTENT(IN   )        :: nCols
    INTEGER, INTENT(IN   )        :: kMax
    INTEGER, INTENT(IN   )        :: istart
    INTEGER, INTENT(IN   )        :: iend
    INTEGER, INTENT(IN   )        :: iloop
    INTEGER, INTENT(OUT  )        :: kbcon   (nCols)
    INTEGER, INTENT(INOUT)        :: k22     (nCols)
    INTEGER, INTENT(INOUT)        :: ierr    (nCols)
    INTEGER, INTENT(IN   )        :: kbmax   (nCols)
    REAL(KIND=r8)   , INTENT(IN   )        :: cap_max (nCols) 
    REAL(KIND=r8)   , INTENT(IN   )        :: he_cup  (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )        :: hes_cup (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )        :: p_cup   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )        :: cap_inc
    !
    ! new
    !
    REAL(KIND=r8)                          :: plus
    REAL(KIND=r8)                          :: pbcdif
    INTEGER                       :: i
    INTEGER :: left(iend-istart+1)
    INTEGER :: nLeft
    INTEGER :: nNewLeft
    INTEGER :: toContinue(iend-istart+1)
    INTEGER :: nToContinue
    INTEGER :: cnt
    !
    ! determine the level of convective cloud base  - kbcon
    !
    kbcon=1
    nLeft = 0
    DO i=istart,iend
       IF(ierr(i) == 0 ) THEN
          kbcon(i)=k22(i)
          nLeft = nLeft + 1
          left(nLeft) = i
       ELSE
          kbcon(i)=1
       END IF
    END DO
    DO
       IF (nLeft == 0) THEN
          EXIT
       ELSE

          nNewLeft = 0
          nToContinue = 0
          !CDIR NODEP
          DO cnt = 1, nLeft
             i = left(cnt)
             IF(he_cup(i,k22(i)) <  hes_cup(i,kbcon(i))) THEN
                kbcon(i)=kbcon(i)+1
                IF(kbcon(i) >  kbmax(i)+2)THEN
                   IF (iloop <  4) THEN
                      ierr(i) =   3
                   ELSE IF (iloop == 4) THEN
                      ierr(i) = 997
                   END IF
                ELSE
                   nNewLeft = nNewLeft + 1
                   left(nNewLeft) = i
                END IF
             ELSE
                nToContinue = nToContinue + 1
                toContinue(nToContinue) = i
             END IF
          END DO

          !CDIR NODEP
          DO cnt = 1, nToContinue
             i = toContinue(cnt)
             IF(kbcon(i)-k22(i) /= 1) THEN
                !
                ! cloud base pressure and max moist static energy pressure
                !
                ! i.e., the depth (in mb) of the layer of negative buoyancy                  
                !
                pbcdif=-p_cup(i,kbcon(i))+p_cup(i,k22(i))
                plus  =MAX(25.0_r8, cap_max(i)-REAL(iloop-1,kind=r8)*cap_inc)   !new
                IF(pbcdif > plus)THEN
                   k22  (i)=k22(i)+1
                   kbcon(i)=k22(i)
                   nNewLeft = nNewLeft + 1
                   left(nNewLeft) = i
                END IF
             END IF
          END DO

          nLeft = nNewLeft

       END IF
    END DO
  END SUBROUTINE cup_kbcon


  SUBROUTINE cup_up_he(k22       ,hkb       ,z_cup     ,cd        ,entr      , &
       he_cup    ,hc        ,nCols     ,kMax      ,kbcon     , &
       ierr      ,istart    ,iend      ,dby       ,he        , &
       hes_cup,entr2D)

    IMPLICIT NONE
    INTEGER, INTENT(IN   )                       :: nCols
    INTEGER, INTENT(IN   )                       :: kMax
    INTEGER, INTENT(IN   )                       :: istart
    INTEGER, INTENT(IN   )                       :: iend
    REAL(KIND=r8)   , INTENT(IN   )                       :: entr
    INTEGER, INTENT(IN   )                       :: kbcon    (nCols)     
    INTEGER, INTENT(IN   )                       :: ierr     (nCols)     
    INTEGER, INTENT(IN   )                       :: k22      (nCols)     
    REAL(KIND=r8)   , INTENT(OUT  )                       :: hkb      (nCols)     
    REAL(KIND=r8)   , INTENT(IN   )                       :: he_cup   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                       :: hc       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: z_cup    (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: cd       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                       :: dby      (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: he       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: hes_cup  (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   )                       :: entr2D (nCols, kMax)
    REAL(KIND=r8):: entr2


    INTEGER                       :: i
    INTEGER                       :: k
    REAL(KIND=r8)                          :: dz
    !
    ! hc = cloud moist static energy
    ! hkb = moist static energy at originating level
    ! he = moist static energy on model levels
    ! he_cup = moist static energy on model cloud levels
    ! hes_cup = saturation moist static energy on model cloud levels
    ! dby = buoancy term
    ! cd= detrainment function
    ! z_cup = heights of model cloud levels
    ! entr = entrainment rate

    !
    ! moist static energy inside cloud
    !
    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          hkb(i)=he_cup(i,k22(i))
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF( ierr(i) == 0 .AND. k < k22(i)  .AND. k<= kbcon(i)-1)THEN
             hc(i,k)=he_cup(i,k)
!!!dby(i,k)=0.0_r8
             dby(i,k)=hc(i,k)-hes_cup(i,k)
          END IF
          IF(ierr(i) == 0 .AND. k >= k22(i)  .AND.k <= kbcon(i)-1)THEN
             hc(i,k)=hkb(i)
!!!!dby(i,k)=0.0_r8
             dby(i,k)=hc(i,k)-hes_cup(i,k)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          hc (i,kbcon(i))= hkb(i)
          dby(i,kbcon(i))= hkb(i)-hes_cup(i,kbcon(i))
       END IF
    END DO
    DO k=1,kMax-1
       DO i=istart,iend
          IF(k >= 2 .AND. k >= kbcon(i).AND.ierr(i) == 0)THEN
             dz=z_cup(i,k)-z_cup(i,k-1)
             entr2=entr2D(i,k)
             hc(i,k)=(hc(i,k-1)*(1.0_r8-0.5_r8*cd(i,k)*dz)+entr2*          &
                  dz*he(i,k-1))/(1.0_r8+entr2*dz-0.5_r8*cd(i,k)*dz)

             dby(i,k)=hc(i,k)-hes_cup(i,k)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_up_he


  SUBROUTINE cup_ktop( &
       ilo      ,dby      ,kbcon    ,ktop     ,nCols     ,kMax     , &
       istart   ,iend     ,ierr)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: ilo
    INTEGER, INTENT(INOUT) :: ierr   (nCols)
    INTEGER, INTENT(IN   ) :: kbcon  (nCols)
    INTEGER, INTENT(OUT  ) :: ktop   (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: dby    (nCols,kMax)

    INTEGER :: i
    INTEGER :: k

    ktop (istart:iend)=1
    DO k=1,kMax-2  
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= kbcon(i)+1 )THEN
             IF(dby(i,k) <= 0.0_r8 .AND. ktop(i) == 1) THEN
                ktop(i)=k-1
             END IF
          END IF
       END DO
    END DO

    DO i=istart,iend      
       IF(ierr(i) == 0 .AND. ktop(i) == 1)THEN    
          IF (ilo == 1) ierr(i)=5
          IF (ilo == 2) ierr(i)=998
       END IF
    END DO

    DO  k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND.k >= ktop(i)+1)THEN
             dby(i,k)=0.0_r8
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_ktop

  !*------

  SUBROUTINE cup_up_nms( &
       zu        ,z_cup     ,entr      ,cd        ,kbcon     , &
       ktop      ,nCols     ,kMax      ,istart    ,iend      , &
       ierr      ,k22,entr2D                                          )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    REAL(KIND=r8)   , INTENT(OUT  ) :: zu   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup(nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cd   (nCols, kMax)
    INTEGER, INTENT(IN   ) :: kbcon(nCols  )
    INTEGER, INTENT(IN   ) :: ktop (nCols  )
    INTEGER, INTENT(IN   ) :: k22  (nCols  )
    INTEGER, INTENT(IN   ) :: ierr (nCols  )
    REAL(KIND=r8)   , INTENT(IN   ) :: entr2D (nCols, kMax)
    REAL(KIND=r8)::entr2

    INTEGER                :: i
    INTEGER                :: k
    REAL(KIND=r8)                   :: dz

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             zu(i,k)=0.0_r8
          END IF
       END DO
    END DO
    DO k=1,kMax 
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= k22(i) .AND. k <= kbcon(i))THEN
             zu(i,k)=1.0_r8
          END IF
       END DO
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= kbcon(i)+1 .AND. k <= ktop(i))THEN
             dz=z_cup(i,k)-z_cup(i,k-1)
             entr2=entr2D(i,k)
             zu(i,k)=zu(i,k-1)*(1.0_r8+(entr2-cd(i,k))*dz)
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cup_up_nms

  !*--------

  SUBROUTINE cup_dd_nms(zd       ,z_cup     ,cdd       ,entr       ,jmin     , &
       ierr     ,nCols     ,kMax       ,istart     ,iend     , &
       itest    ,kdet      ,z1,detr2D)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: itest
    INTEGER, INTENT(IN   ) :: jmin   (nCols)
    INTEGER, INTENT(IN   ) :: ierr   (nCols)
    INTEGER, INTENT(IN   ) :: kdet   (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    REAL(KIND=r8)   , INTENT(OUT  ) :: zd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: cdd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z1     (nCols    )
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D(nCols,kMax)
    REAL(KIND=r8)  ::entr2
    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: ki
    REAL(KIND=r8)                   :: dz
    REAL(KIND=r8)                   :: a
    REAL(KIND=r8)                   :: perc
    !
    ! z_cup = height of cloud model level
    ! z1 = terrain elevation
    ! entr = downdraft entrainment rate
    ! jmin = downdraft originating level
    ! kdet = level above ground where downdraft start detraining
    ! itest = flag to whether to calculate cdd

    !
    ! perc is the percentage of mass left when hitting the ground
    !
    !perc=0.2_r8
    perc=0.03_r8       !it is ok - new.
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             zd(i,k)=0.0_r8
             IF(itest == 0)cdd(i,k)=0.0_r8
          END IF
       END DO
    END DO

    a=1.0_r8-perc

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          zd(i,jmin(i))=1.0_r8
       END IF
    END DO
    DO ki=kMax-1,1,-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. ki <= jmin(i)-1 .AND. ki >= 1)THEN
             !
             ! integrate downward, specify detrainment(cdd)!
             !
             dz=z_cup(i,ki+1)-z_cup(i,ki)
             entr2=detr2D(i,ki)
             IF(ki <= kdet(i).AND.itest == 0)THEN
                cdd(i,ki)=entr2+(1.0_r8- (a*(z_cup(i,ki)-z1(i))          &
                     +perc*(z_cup(i,kdet(i))-z1(i)) )           &
                     /(a*(z_cup(i,ki+1)-z1(i))                  &
                     +perc*(z_cup(i,kdet(i))-z1(i))))/dz
             END IF
             zd(i,ki)=zd(i,ki+1)*(1.0_r8+(entr2-cdd(i,ki))*dz)
             !zd(i,ki)=0.0_r8 !!zd(i,ki)  !!testing nilo10
             !
             !----------------------
             !
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_nms

  !*------

  SUBROUTINE cup_dd_he(hes_cup    ,hcd       ,z_cup     ,cdd        , &
       entr       ,jmin      ,ierr      ,nCols      ,kMax    , &
       istart     ,iend      ,he        ,dby, detr2D     )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cdd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: he      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: dby     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hcd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hes_cup (nCols,kMax)
    INTEGER, INTENT(IN   ) :: jmin    (nCols)
    INTEGER, INTENT(IN   ) :: ierr    (nCols)
    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: ki
    REAL(KIND=r8)                   :: dz
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D(nCols,kMax)
    REAL(KIND=r8):: entr2
    !-------
    ! hcd = downdraft moist static energy
    ! he = moist static energy on model levels
    ! he_cup = moist static energy on model cloud levels
    ! hes_cup = saturation moist static energy on model cloud levels
    ! dby = buoancy term
    ! cdd= detrainment function
    ! z_cup = heights of model cloud levels
    ! entr = entrainment rate
    ! zd   = downdraft normalized mass flux



    DO k=2,kMax
       DO i=istart,iend
          dby(i,k)=0.0_r8
          IF(ierr(I) == 0)THEN
             hcd(i,k)=hes_cup(i,k)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          hcd(i,jmin(i)) = hes_cup(i,jmin(i))
          dby(i,jmin(i)) = hcd    (i,jmin(i)) - hes_cup(i,jmin(i))
       END IF
    END DO

    DO ki=kMax-1,1,-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. ki <= jmin(i)-1 )THEN
             dz        = z_cup(i,ki+1)-z_cup(i,ki)
             entr2=detr2D(i,ki)
             hcd(i,ki) = (hcd(i,ki+1)*(1.0_r8-0.5_r8*cdd(i,ki)*dz)     &
                  + entr2*dz*he(i,ki)  )                   &
                  / (1.0_r8+entr2*dz-0.5_r8*cdd(i,ki)*dz)
             dby(i,ki) = hcd(i,ki)-hes_cup(i,ki)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_he

  !*------

  SUBROUTINE cup_dd_moisture( &
       zd        ,hcd       ,hes_cup   ,qcd       , &
       qes_cup    ,pwd       ,q_cup     ,z_cup     ,cdd       , &
       entr       ,jmin      ,ierr      ,gamma_cup ,pwev      , &
       nCols      ,kMax      ,istart    ,iend      ,bu        , &
       qrcd       ,q         , &
       iloop, detr2D                                                    )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: iloop
    REAL(KIND=r8)   , INTENT(IN   ) :: entr
    INTEGER, INTENT(IN   ) :: jmin      (nCols     )
    INTEGER, INTENT(INOUT) :: ierr      (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: bu        (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwev      (nCols     )
    REAL(KIND=r8)   , INTENT(IN   ) :: zd        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qcd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qrcd      (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hes_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hcd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qes_cup   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: q_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cdd       (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: q         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: detr2D(nCols,kMax)
    REAL(KIND=r8):: entr2

    INTEGER                :: i
    INTEGER                :: k
    INTEGER                :: ki
    REAL(KIND=r8)                   :: dz
    REAL(KIND=r8)                   :: dqeva
    REAL(KIND=r8)                   :: dh
    !------
    ! cdd= detrainment function
    ! q = environmental q on model levels
    ! q_cup = environmental q on model cloud levels
    ! qes_cup = saturation q on model cloud levels
    ! hes_cup = saturation h on model cloud levels
    ! hcd = h in model cloud
    ! bu = buoancy term
    ! zd = normalized downdraft mass flux
    ! gamma_cup = gamma on model cloud levels
    ! mentr_rate = entrainment rate
    ! qcd = cloud q (including liquid water) after entrainment
    ! qrch = saturation q in cloud
    ! pwd = evaporate at that level
    ! pwev = total normalized integrated evaoprate (I2)
    ! entr= entrainment rate


    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          bu  (i) = 0.0_r8
          pwev(i) = 0.0_r8
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             qcd (i,k) = 0.0_r8
             qrcd(i,k) = 0.0_r8
             pwd (i,k) = 0.0_r8
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          k              = jmin(i)
          dz             = z_cup(i,k+1)-z_cup(i,k)
          qcd (i,k)      = q_cup(i,k)
          qrcd(i,k)      = qes_cup(i,k)
          pwd (i,jmin(i))= MIN(0.0_r8,qcd(i,k)-qrcd(i,k))
          pwev(i)        = pwev(i)+pwd(i,jmin(i))
          qcd (i,k)      = qes_cup(i,k)
          dh             = hcd(i,k)-hes_cup(i,k)
          bu(i)          = dz*dh
       END IF
    END DO

    DO ki=kMax-1,1,-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. ki <= jmin(i)-1)THEN
             dz        = z_cup(i,ki+1)-z_cup(i,ki)
             entr2=detr2D(i,ki)
             qcd(i,ki) = (qcd(i,ki+1)*(1.0_r8-0.5_r8*cdd(i,ki)*dz)     &
                  + entr2*dz*q(i,ki)   )                   &
                  / (1.0_r8+entr2*dz-0.5_r8*cdd(i,ki)*dz)
             !
             ! to be negatively buoyant, hcd should be smaller than hes!
             !
             dh         = hcd(i,ki)-hes_cup(i,ki)
             bu  (i)    = bu(i)+dz*dh
             qrcd(i,ki) = qes_cup(i,ki)+(1.0_r8/xl)*(gamma_cup(i,ki)           &
                  / (1.0_r8+gamma_cup(i,ki)))*dh
             dqeva      = qcd(i,ki)-qrcd(i,ki)

             IF(dqeva >  0.0_r8) dqeva=0.0_r8

             pwd (i,ki) = zd  (i,ki)*dqeva
             qcd (i,ki) = qrcd(i,ki)
             pwev(i)    = pwev(i)+pwd(i,ki)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          IF(bu(i) >= 0.0_r8 .AND. iloop == 1)THEN
             ierr(i)=7
          END IF
       END IF
    END DO
    RETURN
  END SUBROUTINE cup_dd_moisture

  !*--------

  SUBROUTINE cup_up_moisture( &
       ierr       ,z_cup     ,qc        ,qrc        ,pw        , &
       pwav       ,kbcon     ,ktop      ,nCols      ,kMax      , &
       istart     ,iend      ,cd        ,dby        ,mentr_rate, &
       q          ,gamma_cup ,zu        ,qes_cup    ,k22       , &
       qe_cup, entr2D                                                    )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: mentr_rate
    INTEGER, INTENT(IN   ) :: kbcon      (nCols)
    INTEGER, INTENT(IN   ) :: ktop       (nCols)
    INTEGER, INTENT(IN   ) :: ierr       (nCols)
    INTEGER, INTENT(IN   ) :: k22        (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pwav       (nCols)     
    REAL(KIND=r8)   , INTENT(IN   ) :: q          (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zu         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup  (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qe_cup     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: dby        (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cd         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z_cup      (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: qes_cup    (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qc         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: qrc        (nCols, kMax)
    REAL(KIND=r8)   , INTENT(OUT  ) :: pw         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: entr2D(nCols, kMax)
    REAL(KIND=r8)  :: entr2


    INTEGER                    :: i
    INTEGER                    :: k
    INTEGER                    :: IALL
    REAL(KIND=r8)                       :: radius2
    REAL(KIND=r8)                       :: dz
    REAL(KIND=r8)                       :: qrch
    REAL(KIND=r8)                       :: c0
    !---------
    ! cd= detrainment function
    ! q = environmental q on model levels
    ! qe_cup = environmental q on model cloud levels
    ! qes_cup = saturation q on model cloud levels
    ! dby = buoancy term
    ! cd= detrainment function
    ! zu = normalized updraft mass flux
    ! gamma_cup = gamma on model cloud levels
    ! mentr_rate = entrainment rate
    !
    ! qc = cloud q (including liquid water) after entrainment
    ! qrch = saturation q in cloud
    ! qrc = liquid water content in cloud after rainout
    ! pw = condensate that will fall out at that level
    ! pwav = totan normalized integrated condensate (I1)
    ! c0 = conversion rate (cloud to rain)



    IALL=0
    c0=0.002_r8
    !
    ! no precip for small clouds
    !
    IF(mentr_rate >  0.0_r8)THEN
       radius2=0.2_r8/mentr_rate
       IF(radius2 <  900.0_r8)c0=0.0_r8
    ENDIF


    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          pwav(i)=0.0_r8
       END IF
    END DO

    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             pw(i,k) =0.0_r8
             !
             !snf        qc(i,k) =qes_cup(i,k)
             !
             qc(i,k) =qe_cup(i,k)   !new
             qrc(i,k)=0.0_r8
          END IF
       END DO
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k >= k22(i) .AND. k<= kbcon(i)-1)THEN
             qc(i,k)=qe_cup(i,k22(i))
          END IF
       END DO
    END DO

    !

    DO k=1,kMax
       DO  i=istart,iend
          IF(ierr(i) == 0 .AND. k >= kbcon(i) .AND. k <= ktop(i) ) THEN
             dz=z_cup(i,k)-z_cup(i,k-1)
             entr2=entr2D(i,k)
             !
             ! 1. steady state plume equation, for what could
             !    be in cloud without condensation
             !
             qc(i,k)=(qc(i,k-1)*(1.0_r8-0.5_r8*cd(i,k)*dz)+entr2*       &
                  dz*q(i,k-1))/(1.0_r8+entr2*dz-0.5_r8*cd(i,k)*dz)
             !
             !2. saturation  in cloud, this is what is allowed to be in it
             !
             qrch=qes_cup(i,k)+(1.0_r8/xl)*(gamma_cup(i,k)       &
                  /(1.0_r8+gamma_cup(i,k)))*dby(i,k)
             !
             ! liquid water content in cloud after rainout
             !
             qrc(i,k)=(qc(i,k)-qrch)/(1.0_r8+c0*dz)
             IF(qrc(i,k) <  0.0_r8)THEN
                qrc(i,k)=0.0_r8
             END IF
             !
             ! 3.Condensation
             !
             pw(i,k)=c0*dz*qrc(i,k)*zu(i,k)
             IF(IALL == 1)THEN
                qrc(i,k)=0.0_r8
                pw(i,k)=(qc(i,k)-qrch)*zu(i,k)
                IF(pw(i,k) <  0.0_r8)pw(i,k)=0.0_r8
             END IF
             !
             ! set next level
             !
             qc(i,k)=qrc(i,k)+qrch
             !
             ! integrated normalized ondensate
             !
             pwav(i)=pwav(i)+pw(i,k)
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cup_up_moisture

  !*---------- 

  SUBROUTINE cup_up_aa0( &
       aa0        ,z         ,zu        ,dby        ,gamma_cup , &
       t_cup      ,kbcon     ,ktop      ,kMax       ,nCols     , &
       istart     ,iend      ,ierr                              )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: kbcon     (nCols)     
    INTEGER, INTENT(IN   ) :: ktop      (nCols)     
    INTEGER, INTENT(IN   ) :: ierr      (nCols)     
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0       (nCols)     
    REAL(KIND=r8)   , INTENT(IN   ) :: z         (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zu        (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: dby       (nCols,kMax)

    REAL(KIND=r8)                         :: dz
    REAL(KIND=r8)                         :: da
    INTEGER                      :: i 
    INTEGER                      :: k

    !
    ! aa0 cloud work function
    ! gamma_cup = gamma on model cloud levels
    ! t_cup = temperature (Kelvin) on model cloud levels
    ! dby = buoancy term
    ! zu= normalized updraft mass flux
    ! z = heights of model levels
    ! ierr error value, maybe modified in this routine
    !


    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          aa0(i)=0.0_r8
       END IF
    END DO

    DO  k=2,kMax-1
       DO  i=istart,iend
          IF(ierr(i) == 0 .AND. k > kbcon(i) .AND. k <= ktop(i) ) THEN
             dz=z(i,k)-z(i,k-1)
             da=zu(i,k)*dz*(9.81_r8/(1004.0_r8*(   &
                  (t_cup(i,k)))))*dby(i,k-1)/ &
                  (1.0_r8+gamma_cup(i,k))
             IF (k == ktop(i) .AND. da <= 0.0_r8) THEN
                CYCLE
             ELSE
                aa0(i)=aa0(i)+da
                IF(aa0(i) <  0.0_r8)aa0(i)=0.0_r8
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_up_aa0

  !*-CINE DELETED---------

  SUBROUTINE cup_dd_edt(ierr      ,us        ,vs        ,z         ,ktop      , &
       kbcon     ,edt       ,p         ,pwav      ,pwev      , &
       nCols     ,kMax      ,istart    ,iend      ,edtmax    , &
       edtmin    ,maxens2   ,edtc      ,vshear    ,sdp       , &
       vws       ,mask      ,edtmax1,maxens22                        )


    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: maxens2,maxens22
    REAL(KIND=r8)   , INTENT(IN   ) :: edtmax
    REAL(KIND=r8)   , INTENT(IN   ) :: edtmin
    INTEGER, INTENT(IN   ) :: ktop   (nCols)
    INTEGER, INTENT(IN   ) :: kbcon  (nCols)
    INTEGER, INTENT(IN   ) :: ierr   (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: us     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: vs     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z      (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: p      (nCols, kMax)     
    REAL(KIND=r8)   , INTENT(OUT  ) :: edt    (nCols     )
    REAL(KIND=r8)   , INTENT(IN   ) :: pwav   (nCols     )
    REAL(KIND=r8)   , INTENT(IN   ) :: pwev   (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: vshear (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: sdp    (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: vws    (nCols     )
    REAL(KIND=r8)   , INTENT(OUT  ) :: edtc   (nCols, maxens2)                !new 
    REAL(KIND=r8)   , INTENT(IN   ) :: edtmax1
    INTEGER, INTENT(IN   ) :: mask   (nCols)    

    REAL(KIND=r8)                         :: pefb
    REAL(KIND=r8)                         :: prezk
    REAL(KIND=r8)                         :: zkbc
    REAL(KIND=r8)                         :: pef
    REAL(KIND=r8)                         :: einc
    REAL(KIND=r8)                         :: aa1
    REAL(KIND=r8)                         :: aa2
    INTEGER                      :: i
    INTEGER                      :: kk
    INTEGER                      :: k
    !
    ! determine downdraft strength in terms of windshear
    ! calculate an average wind shear over the depth of the cloud
    !
    DO i=istart,iend
       IF(ierr(i) == 0) THEN
          edt   (i)=0.0_r8
          vws   (i)=0.0_r8
          sdp   (i)=0.0_r8
          vshear(i)=0.0_r8
       END IF
    END DO

    DO kk = 1,kMax-1
       DO i=istart,iend
          IF(ierr(i) == 0) THEN
             IF(kk  <=  MIN(ktop(i),kMax-1) .AND. kk  >=  kbcon(i)) THEN
                aa1=ABS((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk)))
                aa2=ABS((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))
                vws(i) = vws(i)+(aa1+aa2)*(p(i,kk) - p(i,kk+1))
                sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
             END IF
             IF (kk  ==  kMax-1)vshear(i) = 1.0e3_r8 * vws(i) / sdp(i)
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          pef=(1.591_r8-0.639_r8*vshear(i)+0.0953_r8*(vshear(i)**2)        &
               -0.00496_r8*(vshear(i)**3))
          !-------------------------------------------
          !snf
          IF(mask(i) == 1)THEN
             IF(pef >  edtmax1)pef=edtmax1          
          ELSE
             IF(pef >  edtmax )pef=edtmax
          END IF
          !
          !------------------
          !
          IF(pef <  edtmin)pef=edtmin
          !
          ! cloud base precip efficiency
          !
          zkbc=z(i,kbcon(i))*3.281e-3_r8
          prezk=0.02_r8
          IF(zkbc >  3.0_r8)THEN
             prezk= 0.96729352_r8+zkbc*(-0.70034167_r8+zkbc*(0.162179896_r8+zkbc*(-  &
                  1.2569798e-2_r8+zkbc*(4.2772e-4_r8-zkbc*5.44e-6_r8))))
          END IF

          IF(zkbc >  25.0_r8)THEN
             prezk=2.4_r8
          END IF
          pefb=1.0_r8/(1.0_r8+prezk)
          !
          !          if(pefb >  edtmax)pefb=edtmax
          !-------------------------------------------
          !snf
          IF(mask(i) == 1)THEN
             IF(pefb >  edtmax1)pefb=edtmax1
          ELSE
             IF(pefb >  edtmax)pefb=edtmax
          END IF
          !------------------

          IF(pefb <  edtmin)pefb=edtmin
          edt(i)=1.0_r8-0.5_r8*(pefb+pef)
          !
          ! edt here is 1-precipeff
          !
       END IF
    END DO

    DO k=1,maxens22
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             einc      = edt(i) / REAL(maxens22+1  ,kind=r8  )
             edtc(i,k) = edt(i) - REAL(k,kind=r8)*einc
             !edtc(i,1) = edt(i)*0.75_r8
             !edtc(i,2) = edt(i)*0.50_r8
             !snf-new     !edtc(i,3) =edt(i)*0.25_r8
             ! forcando usar 0.25_r8 quando apenas 1 ensamble 
             !        if(maxens22.eq.1)edtc(i,k) = edt(i)*0.25_r8
          END IF
       END DO
    END DO

    DO k=1,maxens22
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             edtc(i,k)=-edtc(i,k)*pwav(i)/pwev(i)
             !
             !             if(edtc(i,k) >  edtmax)edtc(i,k)=edtmax
             !-------------------------------------------
             !snf
             IF(mask(i) == 1)THEN
                IF(edtc(i,k) >  edtmax1)edtc(i,k)=edtmax1
             ELSE
                IF(edtc(i,k) >  edtmax )edtc(i,k)=edtmax
             END IF
             !------------------
             IF(edtc(i,k) <  edtmin)edtc(i,k)=edtmin
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_edt

  !*--------

  SUBROUTINE cup_dd_aa0( &
       edt        ,ierr      ,aa0       ,jmin      ,gamma_cup , &
       t_cup      ,hcd       ,hes_cup   ,z         ,nCols     , &
       kMax       ,istart    ,iend      ,zd                     )

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: istart
    INTEGER, INTENT(IN   ) :: iend
    INTEGER, INTENT(IN   ) :: jmin      (nCols)     
    INTEGER, INTENT(IN   ) :: ierr      (nCols)     
    REAL(KIND=r8)   , INTENT(IN   ) :: gamma_cup (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: t_cup     (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: z         (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hes_cup   (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: zd        (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: hcd       (nCols, kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: edt       (nCols     )
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0       (nCols     )
    REAL(KIND=r8)    :: dz
    INTEGER :: i
    INTEGER :: k
    INTEGER :: kk


    DO k=1,kMax-1
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k <  jmin(i))THEN
             kk=jmin(i)-k
             !
             ! original
             !
             dz=(z(i,kk)-z(i,kk+1))
             aa0(i)=aa0(i)+zd(i,kk)*edt(i)*dz*(9.81_r8/(1004.0_r8*t_cup(i,kk)))  &
                  *((hcd(i,kk)-hes_cup(i,kk))/(1.0_r8+gamma_cup(i,kk)))
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_dd_aa0

  !*--------

  SUBROUTINE cup_dellabot( &
       he_cup    ,ierr       ,z_cup     ,p_cup     ,hcd         , &
       edt       ,zd         ,cdd       ,he        ,nCols       , &
       kMax       ,istart     ,iend      ,della    ,mentrd_rate,detr2D)

    IMPLICIT NONE
    INTEGER, INTENT(IN   )                   :: nCols
    INTEGER, INTENT(IN   )                   :: kMax
    INTEGER, INTENT(IN   )                   :: istart
    INTEGER, INTENT(IN   )                   :: iend
    REAL(KIND=r8)   , INTENT(IN   )                   :: z_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: p_cup  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: hcd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: zd     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: cdd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: he     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )                   :: della  (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: he_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )                   :: edt    (nCols)
    INTEGER, INTENT(IN   )                   :: ierr   (nCols)
    REAL(KIND=r8)   , INTENT(IN   )                   :: mentrd_rate
    REAL(KIND=r8)   , INTENT(IN   )                   :: detr2D(nCols,kMax)
    REAL(KIND=r8):: entr2


    REAL(KIND=r8)                      :: detdo1
    REAL(KIND=r8)                      :: detdo2
    REAL(KIND=r8)                      :: entdo
    REAL(KIND=r8)                      :: g
    REAL(KIND=r8)                      :: dp
    REAL(KIND=r8)                      :: dz
    REAL(KIND=r8)                      :: subin
    REAL(KIND=r8)                      :: detdo
    INTEGER                   :: i
    g=9.81_r8
    DO i=istart,iend
       della(i,1)=0.0_r8
       IF(ierr(i) == 0) THEN
          dz        =       z_cup(i,2)-z_cup(i,1)
          entr2     = detr2D(i,1)
          dp        = 100.0_r8*(p_cup(i,1)-p_cup(i,2))
          detdo1    = edt(i)*zd(i,2)*cdd(i,1)*dz
          detdo2    = edt(i)*zd(i,1)
          entdo     = edt(i)*zd(i,2)*entr2*dz
          !snf
          subin     =-edt(i)*zd(i,2)
          detdo     = detdo1+detdo2-entdo+subin
          della(i,1)= (  detdo1*0.5_r8*(hcd(i,1)+hcd(i,2))      &
               + detdo2*    hcd(i,1)                   &
               + subin *    he_cup(i,2)                &
               - entdo *    he(i,1)    )*g/dp
       END IF
    END DO
    RETURN
  END SUBROUTINE cup_dellabot

  !*--------

  SUBROUTINE cup_dellas( &
       ierr       ,z_cup      ,p_cup     ,hcd       ,edt       , &
       zd         ,cdd        ,he        ,nCols     ,kMax      , &
       istart     ,iend       ,della     , &
       mentrd_rate,zu         ,cd        ,hc        ,ktop      , &
       k22        ,kbcon      ,mentr_rate,jmin      ,he_cup    , &
       kdet       ,kpbl       ,entr2D,       detr2D              )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )           :: nCols
    INTEGER, INTENT(IN   )           :: kMax
    INTEGER, INTENT(IN   )           :: istart
    INTEGER, INTENT(IN   )           :: iend
    REAL(KIND=r8)   , INTENT(IN   )           :: z_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: p_cup (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: hcd   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: zd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: cdd   (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: he    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(OUT  )           :: della (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: hc    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: cd    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: zu    (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: he_cup(nCols,kMax)
    INTEGER, INTENT(IN   )           :: kbcon (nCols )
    INTEGER, INTENT(IN   )           :: ktop  (nCols )
    INTEGER, INTENT(IN   )           :: k22   (nCols )
    INTEGER, INTENT(IN   )           :: jmin  (nCols )
    INTEGER, INTENT(IN   )           :: ierr  (nCols )
    INTEGER, INTENT(IN   )           :: kdet  (nCols )
    INTEGER, INTENT(IN   )           :: kpbl  (nCols )
    REAL(KIND=r8)   , INTENT(IN   )           :: edt   (nCols )
    REAL(KIND=r8)   , INTENT(IN   )           :: mentrd_rate
    REAL(KIND=r8)   , INTENT(IN   )           :: mentr_rate
    REAL(KIND=r8)   , INTENT(IN   )           :: entr2D (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   )           :: detr2D (nCols,kMax)


    REAL(KIND=r8)                        :: entdo
    REAL(KIND=r8)                        :: g
    REAL(KIND=r8)                        :: dp
    REAL(KIND=r8)                        :: dz
    REAL(KIND=r8)                        :: subin
    REAL(KIND=r8)                        :: detdo
    REAL(KIND=r8)                        :: entup
    REAL(KIND=r8)                        :: detup
    REAL(KIND=r8)                        :: subdown
    REAL(KIND=r8)                        :: entdoj
    REAL(KIND=r8)                        :: entupk
    REAL(KIND=r8)                        :: detupk
    REAL(KIND=r8)                        :: totmas
    INTEGER                     :: ier  (nCols)
    INTEGER                     :: i
    INTEGER                     :: k

    g=9.81_r8
    ier=0
    DO k=2,kMax
       DO i=istart,iend
          della(i,k)=0.0_r8
       END DO
    END DO

    DO k=2,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k <= ktop(i) ) THEN 
             !
             ! specify detrainment of downdraft, has to be consistent
             ! with zd calculations in soundd
             !
             dz    = z_cup(i,k+1)-z_cup(i,k)
             detdo = edt(i)*cdd(i,k)   *dz*zd(i,k+1)
             !! entdo = edt(i)*mentrd_rate*dz*zd(i,k+1)
             entdo = edt(i)*detr2D(i,k)*dz*zd(i,k+1)
             subin = zu(i,k+1)-zd(i,k+1)*edt(i)
             entup = 0.0_r8
             detup = 0.0_r8
             subdown = (zu(i,k)-zd(i,k)*edt(i))                 !new
             entdoj  = 0.0_r8
             entupk  = 0.0_r8
             detupk  = 0.0_r8
             !
             !         if(k >= kbcon(i))entup=mentr_rate*dz*zu(i,k)     !old
             !         subdown=(zu(i,k)-zd(i,k)*edt(i))                 !old
             !
             IF(k >= kbcon(i) .AND. k <  ktop(i))THEN
                !!entup = mentr_rate*dz*zu(i,k)
                entup = entr2D(i,k)*dz*zu(i,k)
                detup = cd(i,k+1) *dz*zu(i,k)
             END IF

             IF(k == jmin(i))THEN
                entdoj  =edt(i)*zd(i,k)
             END IF
             !
             !         if(k == kpbl(i)-1)then                           !old
             !
             IF(k == k22(i)-1)THEN                            !new
                entupk  = zu(i,kpbl(i))
             END IF

             IF(k >  kdet(i))THEN
                detdo   = 0.0_r8
             END IF
             !
             !         if(k == ktop(i)-1)then                          !old
             !
             IF(k == ktop(i)-0)THEN                          !new
                detupk  = zu(i,ktop(i))
                subin   = 0.0_r8
             END IF
             IF(k <  kbcon(i))THEN
                detup   = 0.0_r8
             END IF
             !
             !changed due to subsidence and entrainment
             !
             totmas=subin-subdown+detup-entup-entdo+                       &
                  detdo-entupk-entdoj+detupk

             IF(ABS(totmas) >  1.e-6_r8)THEN                  !test new
!!! nilo3 ier(i)=1
             END IF
             !--
             dp =  100.0_r8*( p_cup(i,k)-p_cup(i,k+1) )
             della(i,k)=(                                                  &
                  subin  * he_cup(i,k+1)                                 &
                  - subdown* he_cup(i,k  )                                 &
                  + detup  * 0.5_r8*( hc(i,k+1)+ hc(i,k))                     &
                  + detdo  * 0.5_r8*(hcd(i,k+1)+hcd(i,k))                     &
                  - entup  * he    (i,k)                                       &
                  - entdo  * he    (i,k)                                       &
                  - entupk * he_cup(i,k22(i))                              &
                  - entdoj * he_cup(i,jmin(i))                             &
                  + detupk * hc(i,ktop(i))                                 &
                  )*g/dp
          END IF
       END DO
    END DO
    IF (ANY(ier /=0)) THEN
       WRITE (0, '( " some ier /= 0; will stop")')
       STOP "** ERROR AT cup_dellas **"
    END IF
    RETURN
  END SUBROUTINE cup_dellas

  !*------

  SUBROUTINE cup_forcing_ens( &
       aa0       ,aa1       ,xaa0      ,mbdt      ,dtime     , &
       xmb       ,ierr      ,nCols     ,kMax      ,istart    , &
       iend      ,xf        ,name      ,mask      ,maxens    , &
       iens      ,iedt      ,maxens3   ,mconv     , &
       omeg      ,k22       ,pr_ens    ,edt       ,kbcon     , &
       ensdim    ,massfln   ,massfld   ,xff_ens3  ,xk        , &
       p_cup     ,ktop      ,ierr2     ,ierr3     ,gdmpar    , &
       maxens22,den1,dcape_zhang,k_zhang,dcape,          &
       d_ucape,d_dcape,cine,sens,tkeMax,tkeMedia,zo_cup )

    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) ::  name
    INTEGER, INTENT(IN)    :: istart
    INTEGER, INTENT(IN)    :: iend
    INTEGER, INTENT(IN)    :: nCols
    INTEGER, INTENT(IN)    :: kMax
    INTEGER, INTENT(IN)    :: maxens
    INTEGER, INTENT(IN)    :: maxens3
    INTEGER, INTENT(IN)    :: ensdim
    INTEGER, INTENT(IN)    :: iens
    INTEGER, INTENT(IN)    :: iedt
    INTEGER, INTENT(IN)    :: maxens22
    INTEGER, INTENT(IN)    :: gdmpar
    INTEGER, INTENT(IN)    :: mask      (nCols)
    INTEGER, INTENT(IN)    :: k22       (nCols)
    INTEGER, INTENT(IN)    :: kbcon     (nCols)
    INTEGER, INTENT(INOUT) :: ierr      (nCols)
    INTEGER, INTENT(IN)    :: ierr2     (nCols)
    INTEGER, INTENT(IN)    :: ierr3     (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: aa0       (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: aa1       (nCols)       
    REAL(KIND=r8)   , INTENT(OUT)   :: xmb       (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: edt       (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: mconv     (nCols)       
    REAL(KIND=r8)   , INTENT(IN)    :: mbdt        
    REAL(KIND=r8)   , INTENT(OUT)   :: xk        (nCols,maxens )   
    REAL(KIND=r8)   , INTENT(OUT)   :: xff_ens3  (nCols,maxens3)   
    REAL(KIND=r8)   , INTENT(IN)    :: omeg      (nCols,kMax)   
    REAL(KIND=r8)   , INTENT(IN)    :: xaa0      (nCols,maxens)
    REAL(KIND=r8)   , INTENT(OUT)   :: xf        (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(IN)    :: pr_ens    (nCols,ensdim)
    REAL(KIND=r8)   , INTENT(OUT)   :: massfln   (nCols,ensdim)
    INTEGER, INTENT(IN)    :: ktop      (nCols)    
    REAL(KIND=r8)   , INTENT(IN)    :: p_cup     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(INOUT) :: massfld
    REAL(KIND=r8)   , INTENT(IN)    :: dtime
    !-new
    REAL(KIND=r8)   , INTENT(IN)    :: den1 (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN)    :: dcape_zhang (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: k_zhang (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: dcape (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: d_ucape (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: d_dcape (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: cine (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: sens (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: tkeMax (nCols)
    REAL(KIND=r8)   , INTENT(IN)    :: tkeMedia (nCols)
    REAL(KIND=r8)   , INTENT(IN)  :: zo_cup(nCols,kMax)

    REAL(KIND=r8)                     :: xff0(nCols) 
    INTEGER                  :: nens
    INTEGER                  :: n
    INTEGER                  :: nens3
    INTEGER                  :: iresult
    INTEGER                  :: iresultd
    INTEGER                  :: iresulte
    INTEGER                  :: i
    INTEGER                  :: k,kk
    INTEGER                  :: nall(nCols,maxens)
    REAL(KIND=r8)                     :: xff_max
    !
    REAL(KIND=r8)    :: a1,aa,bb,cc,ee,ff,gg   
    INTEGER :: kclim(nCols) 
    LOGICAL :: teste2(nCols,maxens3)
    LOGICAL :: teste3(nCols,maxens3)

    ! ierr error value, maybe modified in this routine
    ! pr_ens = precipitation ensemble
    ! xf_ens = mass flux ensembles
    ! massfln = downdraft mass flux ensembles used in next timestep
    ! omeg = omega from large scale model
    ! mconv = moisture convergence from large scale model
    ! zd      = downdraft normalized mass flux
    ! zu      = updraft normalized mass flux
    ! aa0     = cloud work function without forcing effects
    ! aa1     = cloud work function with forcing effects
    ! xaa0    = cloud work function with cloud effects (ensemble dependent)
    ! edt     = epsilon
    ! dir     = "storm motion"
    ! mbdt    = arbitrary numerical parameter
    ! dtime   = dt over which forcing is applied
    ! iact_gr_old = flag to tell where convection was active
    ! kbcon       = LFC of parcel from k22
    ! k22         = updraft originating level
    ! icoic       = flag if only want one closure (usually set to zero!)
    ! name        = deep or shallow convection flag



    DO i=istart,iend
       IF(ierr(i) ==  0 )THEN
          xff0    (i)  =       (aa1(i)-aa0(i))/dtime
          xff_ens3(i,1)=       (aa1(i)-aa0(i))/dtime
          DO kk=1,maxens3 
             xff_ens3(i,kk)=xff_ens3(i,1)
          ENDDO
       END IF
    END DO

    DO nens=1,maxens
       DO   i=istart,iend
          IF(ierr(i) == 0)THEN
             xk(i,nens)=(xaa0(i,nens)-aa1(i))/mbdt
             IF(xk(i,nens) <= 0.0_r8 .AND. xk(i,nens) > -1.0e-6_r8)xk(i,nens)=-1.0e-6_r8
             IF(xk(i,nens) >  0.0_r8 .AND. xk(i,nens) <  1.0e-6_r8)xk(i,nens)= 1.0e-6_r8
             nall(i,nens)=(iens-1)*maxens3*maxens*maxens22       &
                  +(iedt-1)*maxens3*maxens               &
                  +(nens-1)*maxens3
          END IF
       END DO
    END DO
    !
    !-----------------------------------------------
    ! observe the mass flux calculation:
    !-----------------------------------------------!
    ! ne   |     ierr     | mass flux               !
    ! 1    |     ierr =0  |  mf1 = xff_ens3/xk (ne) !
    ! 1    |     ierr >0  |  mf1 =  0               !
    ! 2    |     ierr2=0  |  mf2 = mf1              !
    ! 2    |     ierr2>0  |  mf2 =  0               !
    ! 3    |     ierr3=0  |  mf3 = mf1              !
    ! 3    |     ierr3>0  |  mf3 =  0               !
    ! 
    !
    ! xk(ne) is the same for any 'ne'.    
    !
    ! if ierr2 > 0 (convection was not permited for that cap_max)
    ! then equal to zero the mass flux for the second member of the ensemble (maxens)
    !
    teste2=.TRUE.
    DO nens3=1,maxens3
       DO   i=istart,iend
          IF(ierr(i) == 0 .AND. ierr2(i) > 0)THEN
             xf     (i,nall(i,2)+nens3)=0.0_r8
             massfln(i,nall(i,2)+nens3)=0.0_r8
             teste2(i,nens3)=.FALSE.
          END IF
       END DO
    END DO

    teste3=.TRUE.
    DO nens3=1,maxens3
       DO   i=istart,iend
          IF(ierr(i) == 0 .AND. ierr3(i) > 0)THEN
             xf(i,nall(i,3)+nens3)=0.0_r8
             massfln(i,nall(i,3)+nens3)=0.0_r8
             teste3(i,nens3)=.FALSE.
          END IF
       END DO
    END DO

    DO nens=1,maxens
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. teste2(i,nens) .AND. teste3(i,nens) ) THEN
             !
             !---------------------------------------------------
             !! for every xk, we have maxens3 xffs,
             !! iens is from outermost ensemble (most expensive!
             !! iedt (maxens2 belongs to it)
             !! is from second, next outermost, not so expensive
             !! so, for every outermost loop, we have maxens*maxens2*3
             !! ensembles!!! nall would be 0, if everything is on first
             !! loop index, then nens would start counting, then iedt, then iensi...
             !------------------------------------------------
             !
             iresultd=0
             iresulte=0
             !
             ! check for upwind convection
             !
             iresult=0
             massfld=0.0_r8

             IF(xk(i,nens) <  0.0_r8 .AND. xff0(i) >  0.0_r8)iresultd=1
             iresulte=MAX(iresult,iresultd)
             IF(iresulte == 1)THEN
                !  snf--------
                xff_max=0.05_r8 !!!0.05_r8
                !IF(mask(i) == 1)xff_max=0.10_r8  !!nilo!xfmax   
                !
                IF(xff0(i) >  xff_max)THEN 
                   DO k=1,maxens3 !snilo
                      xf(i,nall(i,nens)+k)=MAX(0.0_r8,-xff_ens3(i,k)/xk(i,nens))+massfld
                   ENDDO
                ELSE
                   DO k=1,maxens3
                      xf(i,nall(i,nens)+k)=massfld
                   ENDDO
                END IF
                !
                !------------------------------------
                gg=zo_cup(i,ktop(i))-zo_cup(i,kbcon(i))
                !IF(gg<5000.0_r8)then
                !ENDIF
                !!write(*,*)'ggggggg', aa,bb,gg,xf(i,nall(i,nens)+17)
                !-----------------------------------
             END IF
          END IF
       END DO
    END DO

    DO nens=1,maxens      
       DO nens3=1,maxens3
          DO i=istart,iend
             IF(ierr(i) == 0)THEN         
                IF(teste2(i,nens) .AND. teste3(i,nens) ) THEN
                   iresultd=0
                   iresulte=0
                   !
                   ! check for upwind convection
                   !
                   iresult=0
                   massfld=0.0_r8
                   IF(xk(i,nens) <  0.0_r8 .AND. xff0(i) >  0.0_r8)iresultd=1
                   iresulte=MAX(iresult,iresultd)
                   IF(iresulte == 1)THEN
                      !
                      !****************************************************************
                      !----- 1d closure ensemble -------------
                      !
                      IF(gdmpar>= 1 .AND. gdmpar <= maxens3)THEN
                         xf(i,nall(i,nens)+nens3)=xf(i,nall(i,nens)+gdmpar)
                      END IF
                      !
                      !-------------------------
                      ! store new for next time step
                      !-------
                      !
                      massfln(i,nall(i,nens)+nens3)=edt(i)*xf(i,nall(i,nens)+nens3)
                      massfln(i,nall(i,nens)+nens3)=MAX(0.0_r8,massfln(i,nall(i,nens)+nens3))
                   END IF
                END IF
             END IF
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE cup_forcing_ens

  !------------------------------------------------------------------------

  SUBROUTINE cup_output_ens( &
       xf_ens    ,ierr      ,dellat    ,dellaq    ,dellaqc    , &
       outt      ,outq      ,outqc     ,pre       ,pw         , &
       xmb       ,ktop      ,nCols     ,kMax       , &
       istart    ,iend      ,maxens2    , &
       maxens    ,iens      ,pr_ens    ,outt_ens  ,maxens3    , &
       ensdim    ,massfln   ,xfac1     ,xfac_for_dn,maxens22   )

    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: nCols
    INTEGER, INTENT(IN   )         :: kMax
    INTEGER, INTENT(IN   )         :: istart
    INTEGER, INTENT(IN   )         :: iend
    INTEGER, INTENT(IN   )         :: ensdim
    INTEGER, INTENT(IN   )         :: maxens2,maxens22
    INTEGER, INTENT(IN   )         :: maxens
    INTEGER, INTENT(IN   )         :: iens
    INTEGER, INTENT(IN   )         :: maxens3
    REAL(KIND=r8)   , INTENT(OUT  )         :: pre     (nCols)            
    REAL(KIND=r8)   , INTENT(OUT  )         :: xmb     (nCols)            
    REAL(KIND=r8)   , INTENT(OUT  )         :: xfac1   (nCols)            
    REAL(KIND=r8)   , INTENT(INOUT)         :: outt    (nCols,kMax)        
    REAL(KIND=r8)   , INTENT(INOUT)         :: outq    (nCols,kMax)   
    !hmjb    REAL(KIND=r8)   , INTENT(OUT  )         :: outt    (nCols,kMax)        
    !hmjb    REAL(KIND=r8)   , INTENT(OUT  )         :: outq    (nCols,kMax)                 
    REAL(KIND=r8)   , INTENT(OUT  )         :: outqc   (nCols,kMax)        
    REAL(KIND=r8)   , INTENT(IN   )         :: dellat  (nCols,kMax,maxens2)
    REAL(KIND=r8)   , INTENT(IN   )         :: dellaq  (nCols,kMax,maxens2)
    REAL(KIND=r8)   , INTENT(IN   )         :: pw      (nCols,kMax,maxens2)
    REAL(KIND=r8)   , INTENT(IN   )         :: xf_ens  (nCols,ensdim)     
    REAL(KIND=r8)   , INTENT(INOUT)         :: pr_ens  (nCols,ensdim)     
    REAL(KIND=r8)   , INTENT(INOUT)         :: massfln (nCols,ensdim)     
    REAL(KIND=r8)   , INTENT(INOUT)         :: outt_ens(nCols,ensdim)
    INTEGER, INTENT(IN   )         :: ktop    (nCols)
    INTEGER, INTENT(INOUT)         :: ierr    (nCols)
    REAL(KIND=r8)   , INTENT(OUT  )         :: xfac_for_dn(nCols)
    INTEGER          :: ncount  (nCols)
    !
    ! new
    !
    INTEGER                     :: i
    INTEGER                     :: k
    INTEGER                     :: n
    REAL(KIND=r8)                        :: outtes
    REAL(KIND=r8)                        :: ddtes
    REAL(KIND=r8)                        :: dtt     (nCols,kMax)
    REAL(KIND=r8)                        :: dtq     (nCols,kMax)
    REAL(KIND=r8)                        :: dtqc    (nCols,kMax)
    REAL(KIND=r8)                        :: dtpw    (nCols,kMax)
    REAL(KIND=r8)                        :: dellaqc (nCols,kMax,maxens2)
    !
    ! xf_ens = ensemble mass fluxes
    ! pr_ens = precipitation ensembles
    ! dellat = change of temperature per unit mass flux of cloud ensemble
    ! dellaq = change of q per unit mass flux of cloud ensemble
    ! dellaqc = change of qc per unit mass flux of cloud ensemble
    ! outtem = output temp tendency (per s)
    ! outq   = output q tendency (per s)
    ! outqc  = output qc tendency (per s)
    ! pre    = output precip
    ! xmb    = total base mass flux
    ! xfac1  = correction factor
    ! pw = pw -epsilon*pd (ensemble dependent)
    ! ierr error value, maybe modified in this routine


    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             outt (i,k) = 0.0_r8
             outq (i,k) = 0.0_r8
             outqc(i,k) = 0.0_r8
             dtt  (i,k) = 0.0_r8
             dtq  (i,k) = 0.0_r8
             dtqc (i,k) = 0.0_r8
             dtpw (i,k) = 0.0_r8
          END IF
       END DO
    END DO

    DO i=istart,iend
       IF(ierr(i) == 0)THEN
          pre  (i)      = 0.0_r8
          xmb  (i)      = 0.0_r8
          xfac1(i)      = 1.0_r8
          xfac_for_dn(i)= 1.0_r8
          ncount(i)     = 0
       END IF
    END DO
    !
    ! calculate mass fluxes
    !
    ! Simple average  (OLD)
    ! --------------
    !
    DO n=(iens-1)*maxens22*maxens*maxens3+1, iens*maxens22*maxens*maxens3 
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             pr_ens  (i,n) = pr_ens  (i,n)*xf_ens (i,n)
             outt_ens(i,n) = outt_ens(i,n)*xf_ens (i,n)
             IF(xf_ens(i,n) >= 0.0_r8)THEN
                xmb   (i) = xmb   (i) + xf_ens(i,n)
                ncount(i) = ncount(i) + 1
             END IF
          END IF
       END DO
    END DO
    DO i=istart,iend
       IF(ierr(i) == 0)THEN          
          IF(ncount(i) >  0)THEN
             xmb (i)=xmb(i)/REAL(ncount(i),kind=r8)
          ELSE
             xmb (i)=0.0_r8
             ierr(i)=13
          END IF
          xfac1(i)=xmb(i)!new1
       END IF
    END DO
    !
    !--------------------
    !! now do feedback 
    !----------------------
    !
    ddtes=250.0_r8   
    !
    DO n=1,maxens22
       DO k=1,kMax
          DO i=istart,iend
             IF(ierr(i) == 0 .AND. k <= ktop(i))THEN
                dtt (i,k)  = dtt (i,k) + dellat (i,k,n)
                dtq (i,k)  = dtq (i,k) + dellaq (i,k,n)
                dtqc(i,k)  = dtqc(i,k) + dellaqc(i,k,n)
                dtpw(i,k)  = dtpw(i,k) + pw     (i,k,n)
             END IF
          END DO
       END DO
    END DO
    DO k=1,kMax
       DO i=istart,iend
          IF(ierr(i) == 0 .AND. k <= ktop(i))THEN
             outtes = dtt(i,k)*xmb(i)*86400.0_r8/REAL(maxens22,kind=r8)
             IF(outtes  >   2.0_r8*ddtes .AND. k >  2)THEN
                xmb(i) = 2.0_r8*ddtes/outtes * xmb(i)
                outtes = 1.0_r8*ddtes
             END IF


             IF(outtes  <   -ddtes)THEN
                xmb(i) = -ddtes/outtes * xmb(i)
                outtes = -ddtes
             END IF

             IF(outtes  >   0.5_r8*ddtes .AND. k <= 2)THEN
                xmb(i) =    ddtes/outtes * xmb(i)
                outtes = 0.5_r8*ddtes
             END IF
             outt (i,k) = outt (i,k) + xmb(i)*dtt (i,k)/REAL(maxens22,kind=r8)
             outq (i,k) = outq (i,k) + xmb(i)*dtq (i,k)/REAL(maxens22,kind=r8)
             outqc(i,k) = outqc(i,k) + xmb(i)*dtqc(i,k)/REAL(maxens22,kind=r8)
             pre  (i)   = pre  (i)   + xmb(i)*dtpw(i,k)/REAL(maxens22,kind=r8)
          END IF
       END DO
    END DO
    !
    ! below is new  it is only for statistics?
    !
    DO k=(iens-1)*maxens22*maxens*maxens3+1,iens*maxens22*maxens*maxens3
       DO i=istart,iend
          IF(ierr(i) == 0)THEN
             xfac1   (i) = xmb(i)        / (xfac1(i)+1.e-16_r8)
             massfln (i,k) = massfln (i,k) * xfac1(i)*xfac_for_dn(i)
             pr_ens  (i,k) = pr_ens  (i,k) * xfac1(i)
             outt_ens(i,k) = outt_ens(i,k) * xfac1(i)
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE cup_output_ens
  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    REAL(KIND=r8), INTENT(IN) :: t

    IF (t <= tcrit) THEN
       es5 = EXP(ae(2)-be(2)/t)
    ELSE
       es5 = EXP(ae(1)-be(1)/t)
    END IF
  END FUNCTION es5
END MODULE Cu_GDM_BAM
