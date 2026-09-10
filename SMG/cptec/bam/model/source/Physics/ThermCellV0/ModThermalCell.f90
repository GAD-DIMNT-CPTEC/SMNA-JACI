!       thermcell_dv2
!       thermcell_flux2
!       thermcell_plume|______
!                             thermcell_qsat
!                             thermcell_dry
!                             thermcell_closure
!                             thermcell_condens
MODULE ThermalCell
    IMPLICIT NONE
  SAVE

  PRIVATE
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(13,60) ! the '60' maps to 64-bit real
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  REAL(KIND=r8), PARAMETER ::  RG=9.80665_r8
  REAL(KIND=r8), PARAMETER ::  RKBOL=1.380658E-23_r8
  REAL(KIND=r8), PARAMETER ::  RNAVO=6.0221367E+23_r8
  REAL(KIND=r8), PARAMETER ::  RMD=28.9644_r8
  REAL(KIND=r8), PARAMETER ::  RMV=18.0153_r8
  REAL(KIND=r8), PARAMETER ::  RTT=273.16_r8
  REAL(KIND=r8), PARAMETER ::  R4IES=7.66_r8
  REAL(KIND=r8), PARAMETER ::  R3IES=21.875_r8
  REAL(KIND=r8), PARAMETER ::  R3LES=17.269_r8
  REAL(KIND=r8), PARAMETER ::  R4LES=35.86_r8
  REAL(KIND=r8), PARAMETER ::  RLSTT=2.8345E+6_r8
  REAL(KIND=r8), PARAMETER ::  RLVTT=2.5008E+6_r8

  REAL(KIND=r8), PARAMETER ::  R=RNAVO*RKBOL

  REAL(KIND=r8), PARAMETER ::  RD=1000.0_r8*R/RMD
  REAL(KIND=r8), PARAMETER ::  RV=1000.0_r8*R/RMV
  REAL(KIND=r8), PARAMETER ::  RESTT=611.14_r8
  REAL(KIND=r8), PARAMETER ::  RCPV=4.0_r8 *RV

  REAL(KIND=r8), PARAMETER ::  R2ES=RESTT*RD/RV
  REAL(KIND=r8), PARAMETER ::  RCPD=3.5_r8*RD
  REAL(KIND=r8), PARAMETER ::  R5IES=R3IES*(RTT-R4IES)
  REAL(KIND=r8), PARAMETER ::  R5LES=R3LES*(RTT-R4LES)

  REAL(KIND=r8), PARAMETER ::  RETV=RV/RD-1.0_r8
  REAL(KIND=r8), PARAMETER ::  RVTMP2=RCPV/RCPD-1.0_r8

  REAL(KIND=r8), PARAMETER ::  RKAPPA=RD/RCPD

  REAL(KIND=r8), PARAMETER :: DDT0=0.01_r8

  REAL(KIND=r8), PARAMETER :: alp_bl_k = 1.0_r8
  INTEGER      , PARAMETER :: iflag_thermals_optflux=0 ! =0 orig
  INTEGER      , PARAMETER :: iflag_thermals_ed=10   !1 =orig
  INTEGER      , PARAMETER :: iflag_thermals=16
  REAL(KIND=r8), PARAMETER :: r_aspect_thermals=2.0_r8
!  REAL(KIND=r8), PARAMETER :: l_mix_thermals=30.0_r8
!  REAL(KIND=r8), PARAMETER :: tau_thermals = 7200.0_r8!1800.0_r8

  INTEGER      , PARAMETER :: iflag_coupl=3
  INTEGER      , PARAMETER :: prt_level=1

  LOGICAL      , PARAMETER :: plusqueun=.TRUE.
  LOGICAL      , PARAMETER :: centre=.FALSE.
  LOGICAL      , PARAMETER :: CFL_TEST=.TRUE.

  PUBLIC :: Diver_ThermCell
CONTAINS

  SUBROUTINE Diver_ThermCell(ncols    ,kMax    ,timestep,kt       ,&
                             jdt      ,initlz  ,nClass,nAeros,microphys,iswrad,&
                             prsi     ,prsl    ,phii    , phil   ,&
                             gt       ,gq      ,gu      ,gv      ,&
                             colrad   ,f0      ,fm0     ,entr0   ,&
                             detr0    ,pblh    ,imask   ,dump    ,&
                             dtv      ,dqv     ,duv     ,&
                             dvv     ,&
                             gice     ,gliq    ,gvar)
    IMPLICIT NONE
    INTEGER , INTENT(IN   ) :: ncols
    INTEGER , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: timestep
    INTEGER , INTENT(IN   ) :: kt
    INTEGER , INTENT(IN   ) :: jdt
    INTEGER , INTENT(IN   ) :: initlz
    INTEGER , INTENT(IN   ) :: nClass
    INTEGER , INTENT(IN   ) :: nAeros
    LOGICAL , INTENT(IN   ) :: microphys
    CHARACTER(LEN=*), INTENT(IN   ) :: iswrad
    REAL(KIND=r8), INTENT(IN   ) :: prsi   (1:ncols,1:kMax+1)
    REAL(KIND=r8), INTENT(IN   ) :: prsl   (1:ncols,1:kMax) 
    REAL(KIND=r8), INTENT(IN   ) :: phii   (1:nCols,1:kMax+1)     
    REAL(KIND=r8), INTENT(IN   ) :: phil   (1:nCols,1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: gt     (1:nCols,1:kMax) 
    REAL(KIND=r8), INTENT(IN   ) :: gq     (1:nCols,1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: gu     (1:nCols,1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: gv     (1:nCols,1:kMax)
    REAL(KIND=r8), INTENT(IN   ) :: colrad (1:nCols)
    REAL(KIND=r8), INTENT(INOUT) :: f0     (nCols)       
    REAL(KIND=r8), INTENT(INOUT) :: fm0    (nCols,kMax+1)
    REAL(KIND=r8), INTENT(INOUT) :: entr0  (nCols,kMax)
    REAL(KIND=r8), INTENT(INOUT) :: detr0  (nCols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: pblh   (1:nCols)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask(1:ncols)
    REAL(kind=r8), INTENT(INOUT) :: dump   (1:nCols,1:kMax)  
    REAL(kind=r8), INTENT(OUT  ) :: dtv    (nCols,kMax)        ! temperature tendency (heating)
    REAL(kind=r8), INTENT(OUT  ) :: dqv    (nCols,kMax,1:5+nClass+nAeros)! constituent diffusion tendency
    REAL(kind=r8), INTENT(OUT  ) :: duv    (nCols,kMax)        ! u-wind tendency
    REAL(kind=r8), INTENT(OUT  ) :: dvv    (nCols,kMax)        ! v-wind tendency
    REAL(KIND=r8), INTENT(in   ) :: gice   (ncols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: gliq   (ncols,kMax) 
    REAL(KIND=r8),OPTIONAL,  INTENT(inout) :: gvar (ncols,kMax,nClass+nAeros)


    ! LOCAL VARIABLE
    REAL(KIND=r8) :: fact(ncols)

    REAL(KIND=r8) :: pplay(ncols,kMax)  ! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
    REAL(KIND=r8) :: pplev(ncols,kMax+1)! paprs---input-R-pression pour chaque inter-couche (en Pa)
    REAL(KIND=r8) :: pphi(ncols,kMax)   ! pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
    REAL(KIND=r8) :: ppii(ncols,kMax+1)   ! pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
    REAL(KIND=r8) :: pu(ncols,kMax)     ! u-------input-R-vitesse dans la direction X (de O a E) en m/s
    REAL(KIND=r8) :: pv(ncols,kMax)     ! v-------input-R-vitesse Y (de S a N) en m/s
    REAL(KIND=r8) :: pt(ncols,kMax)     ! t-------input-R-temperature (K)
    REAL(KIND=r8) :: po(ncols,kMax)     ! qx------input-R-humidite specifique (kg/kg) et d'autres traceurs

    REAL(KIND=r8) :: pice(ncols,kMax)                   ! qx------input (kg/kg) et d'autres traceurs
    REAL(KIND=r8) :: pliq(ncols,kMax)                   ! qx------input (kg/kg) et d'autres traceurs
    REAL(KIND=r8) :: pvar(ncols,kMax,1:nClass+nAeros)     ! qx------input (kg/kg) et d'autres traceurs

    REAL(KIND=r8) :: pduadj(nCols,kMax)
    REAL(KIND=r8) :: pdvadj(nCols,kMax)
    REAL(KIND=r8) :: pdtadj(nCols,kMax)
    REAL(KIND=r8) :: pdoadj(nCols,kMax)
    REAL(KIND=r8) :: pdicedj(nCols,kMax)
    REAL(KIND=r8) :: pdliqdj(nCols,kMax)

    REAL(KIND=r8) :: pdaedj(nCols,kMax,1:5+nClass+nAeros)
    REAL(KIND=r8) :: zqta(nCols,kMax)
    REAL(KIND=r8) :: zqla(nCols,kMax)
    INTEGER       :: lmax(nCols)
    REAL(KIND=r8) :: ratqscth(nCols,kMax)
    REAL(KIND=r8) :: ratqsdiff(nCols,kMax)
    REAL(KIND=r8) :: zqsatth(nCols,kMax) 
    REAL(KIND=r8) :: Ale_bl(nCols)
    REAL(KIND=r8) :: Alp_bl(nCols)
    INTEGER       :: lalim_conv(nCols)
    REAL(KIND=r8) :: wght_th(nCols,kMax)
    REAL(KIND=r8) :: zmax0(nCols)!diagnostics
    REAL(KIND=r8) :: zw2(nCols,kMax+1)
    REAL(KIND=r8) :: fraca(nCols,kMax+1)
    REAL(KIND=r8) :: ztv(nCols,kMax)
    REAL(KIND=r8) :: zpspsk(nCols,kMax)
    REAL(KIND=r8) :: ztla(nCols,kMax)
    REAL(KIND=r8) :: zthl(nCols,kMax)
    REAL(KIND=r8) :: ptimestep
    LOGICAL       ::  logexpr1(nCols)
    REAL(KIND=r8)  :: tau_thermals
    LOGICAL       :: debut
    LOGICAL       :: flag_bidouille_stratocu
    INTEGER       :: nsplit_thermals=1
    INTEGER       :: i,k,kk
    REAL(KIND=r8) :: Factor=1.0_r8 
fact=0.0_r8
dtv=0.0_r8
dqv=0.0_r8
duv=0.0_r8
dvv=0.0_r8
pplay=0.0_r8  ! pplay---input-R-press
pplev=0.0_r8! paprs---input-R-press
pphi=0.0_r8   ! pphi----input-R-geopo
ppii=0.0_r8   ! pphi----input-R-geo
pu=0.0_r8     ! u-------input-R-vites
pv=0.0_r8     ! v-------input-R-vites
pt=0.0_r8     ! t-------input-R-tempe
po=0.0_r8     ! qx------input-R-humid
pduadj=0.0_r8
pdvadj=0.0_r8
pdtadj=0.0_r8
pdoadj=0.0_r8
pdicedj=0.0_r8
pdliqdj=0.0_r8
pdaedj=0.0_r8
zqta=0.0_r8
zqla=0.0_r8
lmax=0.0_r8
ratqscth=0.0_r8
ratqsdiff=0.0_r8
zqsatth=0.0_r8 
Ale_bl=0.0_r8
Alp_bl=0.0_r8
lalim_conv=0.0_r8
wght_th=0.0_r8
zmax0=0.0_r8!diagnostics
zw2=0.0_r8 
fraca=0.0_r8 
ztv=0.0_r8
zpspsk=0.0_r8
ztla=0.0_r8
zthl=0.0_r8
ptimestep=0.0_r8
pvar=0.0_r8
    nsplit_thermals=1
    ptimestep=2*timestep
    !tau_thermals=21600.0_r8!6*timestep
    tau_thermals=(2*timestep)/log(0.89483932_r8)

!          lambda=EXP(-ptimestep/tau_thermals)
!          ln(lambda)=-ptimestep/tau_thermals
!          tau_thermals=(2*timestep)/ln(lambda)
    IF(initlz >= 0 .AND. kt == 0 .AND. jdt == 1)THEN
      debut=.TRUE.
    ELSE
      debut=.FALSE.
    END IF  
    DO i=1,ncols
          pplev(i,kMax+1) = prsi(i,kMax+1)
          ppii (i,kMax+1) = phii(i,kMax+1)
    END DO
    IF (TRIM(iswrad).eq.'RRTMG') THEN
       Factor=0.1_r8 
    ELSE  IF (TRIM(iswrad).eq.'CRD') THEN
       Factor=1.0_r8 
    ELSE  IF (TRIM(iswrad).eq.'CRDTF') THEN
       Factor=1.0_r8 
    ELSE
       Factor=1.0_r8 
    END IF
    DO k=1,kMax
       DO i=1,ncols
          pplay(i,k) = prsl(i,k)  
          pplev(i,k) = prsi(i,k)
          ppii(i,k)  = phii(i,k)
          pphi(i,k)  = phil(i,k)*RG
          pu(i,k)    = gu  (i,k)/SIN( colrad(i))
          pv(i,k)    = gv  (i,k)/SIN( colrad(i))
          pt(i,k)    = gt  (i,k) 
          po(i,k)    = gq  (i,k)
          pice(i,k)  = gice(i,k)
          pliq(i,k)  = gliq(i,k)
       END DO
    END DO
    IF (microphys) THEN
       IF( (nClass+nAeros)>0 .and. PRESENT(gvar))THEN
          DO kk=1,nClass+nAeros
             DO k=1,kMax
                DO i=1,ncols
                    pvar(i,k,kk) =gvar (i,k,kk)
                END DO
             END DO
          END DO
      END IF
    END IF  
!fm0    =0.0_r8
!entr0  =0.0_r8
!detr0  =0.0_r8
    CALL thermcell_main( &
      nCols       ,&!INTEGER      , INTENT(IN   ) :: nCols
      kMax        ,&!INTEGER      , INTENT(IN   ) :: kMax
      nClass     ,&
      nAeros     ,&
      debut       ,&!LOGICAL      , INTENT(IN   ) :: debut
      ptimestep   ,&!REAL(KIND=r8), INTENT(IN   ) :: ptimestep
      tau_thermals  ,&
      pplay     (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(IN   ) :: pplay(nCols,kMax)
      pplev     (1:nCols,1:kMax+1),&!REAL(KIND=r8), INTENT(IN   ) :: pplev(nCols,kMax+1)
      pphi      (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(IN   ) :: pphi(nCols,kMax)
      pu        (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(IN   ) :: pu(nCols,kMax)
      pv        (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(IN   ) :: pv(nCols,kMax)
      pt        (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(IN   ) :: pt(nCols,kMax)
      po        (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(INOUT) :: po(nCols,kMax)
      pice      (1:nCols,1:kMax)  ,&!
      pliq      (1:nCols,1:kMax)  ,&!
      pvar      (1:ncols,1:kMax,1:nClass+nAeros), &
      pduadj    (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: pduadj(nCols,kMax)
      pdvadj    (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdvadj(nCols,kMax)
      pdtadj    (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdtadj(nCols,kMax)
      pdoadj    (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdoadj(nCols,kMax)
      pdicedj   (1:nCols,1:kMax)  ,&
      pdliqdj   (1:nCols,1:kMax)  ,&
      pdaedj    (1:nCols,1:kMax,1:5+nClass+nAeros)  ,&!
      fm0       (1:nCols,1:kMax+1),&!REAL(KIND=r8), INTENT(INOUT) :: fm0(nCols,kMax+1)
      entr0     (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(INOUT) :: entr0(nCols,kMax)
      detr0     (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(INOUT) :: detr0(nCols,kMax)
      pblh      (1:nCols       )  ,&!INTEGER      , INTENT(IN   ) :: pblh(nCols)
      imask     (1:nCols       )  ,&!INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
      dump      (1:nCols,1:kMax)  ,&!
      zqta      (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqta(nCols,kMax)
      zqla      (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqla(nCols,kMax)
      lmax      (1:nCols       )  ,&!INTEGER      , INTENT(OUT  ) :: lmax(nCols)
      ratqscth  (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ratqscth(nCols,kMax)
      ratqsdiff (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ratqsdiff(nCols,kMax)
      zqsatth   (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(nCols,kMax) 
      Ale_bl    (1:nCols       )  ,&!REAL(KIND=r8), INTENT(OUT  ) :: Ale_bl(nCols)
      Alp_bl    (1:nCols       )  ,&!REAL(KIND=r8), INTENT(OUT  ) :: Alp_bl(nCols)
      lalim_conv(1:nCols       )  ,&!INTEGER      , INTENT(OUT  ) :: lalim_conv(nCols)
      wght_th   (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: wght_th(nCols,kMax)
      zmax0     (1:nCols       )  ,&!REAL(KIND=r8), INTENT(INOUT) :: zmax0(nCols)           !diagnostics
      f0        (1:nCols       )  ,&!REAL(KIND=r8), INTENT(INOUT) :: f0(nCols)           !diagnostics
      zw2       (1:nCols,1:kMax+1),&!REAL(KIND=r8), INTENT(OUT  ) :: zw2(nCols,kMax+1)
      fraca     (1:nCols,1:kMax+1)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: fraca(nCols,kMax+1)
      ztv       (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztv(nCols,kMax)
      zpspsk    (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zpspsk(nCols,kMax)
      ztla      (1:nCols,1:kMax)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztla(nCols,kMax)
      zthl      (1:nCols,1:kMax)  )!REAL(KIND=r8), INTENT(OUT  ) :: zthl(nCols,kMax)

!      flag_bidouille_stratocu=iflag_thermals.eq.14.or.iflag_thermals.eq.16


!      fact(:)=0.0_r8

!      DO i=1,nCols
!       logexpr1(i)=flag_bidouille_stratocu
!       IF(logexpr1(i)) fact(i)=1./REAL(nsplit_thermals)
!      ENDDO
!    DO k=1,kMax
!!  transformation de la derivee en tendance
!            pduadj(:,k)=pduadj(:,k)*ptimestep*fact(:)
!            pdvadj(:,k)=pdvadj(:,k)*ptimestep*fact(:)
!            pdtadj(:,k)=pdtadj(:,k)*ptimestep*fact(:)
!            pdoadj(:,k)=pdoadj(:,k)*ptimestep*fact(:)

!            entr0  (:,k)=entr0  (:,k) +zfm_therm  (:,k)*fact(:)
!            entr_therm(:,k)=entr_therm(:,k) +zentr_therm(:,k)*fact(:)
!            detr_therm(:,k)=detr_therm(:,k) +zdetr_therm(:,k)*fact(:)
!    ENDDO
!    fm0(:,kMax+1)=0.0_r8
    DO k=1,kMax
       DO i=1,ncols
!          PRINT*, i,k,pduadj(i,k),pdvadj(i,k), pdtadj(i,k), pdoadj(i,k)

    !
    ! Convert the diffused fields back to diffusion tendencies.
    ! Add the diffusion tendencies to the cummulative physics tendencies,
    ! except for constituents. The diffused values of the constituents
    ! replace the input values.
    !
          duv(i,k) =  Factor*pduadj(i,k)*SIN( colrad(i)) ! (up1(i,k)*SIN( colrad(i)) - um1(i,k)*SIN( colrad(i)))*rztodt
          dvv(i,k) =  Factor*pdvadj(i,k)*SIN( colrad(i)) !(vp1(i,k)*SIN( colrad(i)) - vm1(i,k)*SIN( colrad(i)))*rztodt
          dtv(i,k) =  Factor*pdtadj(i,k)
       END DO


       DO i=1,ncols
             dqv(i,k,3) =  Factor*pdoadj(i,k)  !(qp1(i,k,m) - qm1(i,k,m))*rztodt
             dqv(i,k,4) =  Factor*pdicedj(i,k) !(qp1(i,k,m) - qm1(i,k,m))*rztodt
             dqv(i,k,5) =  Factor*pdliqdj(i,k) !(qp1(i,k,m) - qm1(i,k,m))*rztodt
       END DO

    END DO
    IF((nClass+nAeros)>0 .and. PRESENT(gvar))THEN
        DO kk=1,nClass+nAeros
           DO k=1,kmax
              DO i=1,ncols
                   dqv(i,k,6+kk-1)= Factor*pdaedj(i,k,6+kk-1)
              END DO
           END DO
        END DO
    END IF

    RETURN
  END SUBROUTINE Diver_ThermCell

  !
  ! $Id: thermcell_main.F90 1525 2011-05-25 10:55:27Z idelkadi $
  !
  SUBROUTINE thermcell_main( &
       ngrid      ,&!INTEGER      , INTENT(IN   ) :: ngrid
       nlay       ,&!INTEGER      , INTENT(IN   ) :: nlay
       nClass     ,&
       nAeros     ,&
       debut      ,&!LOGICAL      , INTENT(IN   ) :: debut
       ptimestep  ,&!REAL(KIND=r8), INTENT(IN   ) :: ptimestep
       tau_thermals,&
       pplay      ,&!REAL(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
       pplev      ,&!REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
       pphi       ,&!REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
       pu         ,&!REAL(KIND=r8), INTENT(IN   ) :: pu(ngrid,nlay)
       pv         ,&!REAL(KIND=r8), INTENT(IN   ) :: pv(ngrid,nlay)
       pt         ,&!REAL(KIND=r8), INTENT(IN   ) :: pt(ngrid,nlay)
       po         ,&!REAL(KIND=r8), INTENT(INOUT) :: po(ngrid,nlay)
       pice       ,&!
       pliq       ,&!
       pvar       , &
       pduadj     ,&!REAL(KIND=r8), INTENT(OUT  ) :: pduadj(ngrid,nlay)
       pdvadj     ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdvadj(ngrid,nlay)
       pdtadj     ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdtadj(ngrid,nlay)
       pdoadj     ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdoadj(ngrid,nlay)
       pdicedj    ,&
       pdliqdj    ,&
       pdaedj     ,&!REAL(KIND=r8), INTENT(OUT  ) :: pdaedj(ngrid,nlay)
       fm0        ,&!REAL(KIND=r8), INTENT(INOUT) :: fm0(ngrid,nlay+1)
       entr0      ,&!REAL(KIND=r8), INTENT(INOUT) :: entr0(ngrid,nlay)
       detr0      ,&!REAL(KIND=r8), INTENT(INOUT) :: detr0(ngrid,nlay)
       pblh       ,&!INTEGER      , INTENT(IN   ) :: pblh(nCols)
       imask      ,&!INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
       dump       ,&!INTEGER      , INTENT(IN   ) :: pblh(nCols)
       zqta       ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqta(ngrid,nlay)
       zqla       ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqla(ngrid,nlay)
       lmax       ,&!INTEGER      , INTENT(OUT  ) :: lmax(ngrid)
       ratqscth   ,&!REAL(KIND=r8), INTENT(OUT  ) :: ratqscth(ngrid,nlay)
       ratqsdiff  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ratqsdiff(ngrid,nlay)
       zqsatth    ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(ngrid,nlay) 
       Ale_bl     ,&!REAL(KIND=r8), INTENT(OUT  ) :: Ale_bl(ngrid)
       Alp_bl     ,&!REAL(KIND=r8), INTENT(OUT  ) :: Alp_bl(ngrid)
       lalim_conv ,&!INTEGER      , INTENT(OUT  ) :: lalim_conv(ngrid)
       wght_th    ,&!REAL(KIND=r8), INTENT(OUT  ) :: wght_th(ngrid,nlay)
       zmax0      ,&!REAL(KIND=r8), INTENT(INOUT) :: zmax0(ngrid)           !diagnostics
       f0         ,&!REAL(KIND=r8), INTENT(INOUT) :: f0(ngrid)           !diagnostics
       zw2        ,&!REAL(KIND=r8), INTENT(OUT  ) :: zw2(ngrid,nlay+1)
       fraca      ,&!REAL(KIND=r8), INTENT(OUT  ) :: fraca(ngrid,nlay+1)
       ztv        ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztv(ngrid,nlay)
       zpspsk     ,&!REAL(KIND=r8), INTENT(OUT  ) :: zpspsk(ngrid,nlay)
       ztla       ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztla(ngrid,nlay)
       zthl        )!REAL(KIND=r8), INTENT(OUT  ) :: zthl(ngrid,nlay)

    !      USE dimphy
    IMPLICIT NONE

    !=======================================================================
    !   Auteurs: Frederic Hourdin, Catherine Rio, Anne Mathieu
    !   Version du 09.02.07
    !   Calcul du transport vertical dans la couche limite en presence
    !   de "thermiques" explicitement representes avec processus nuageux
    !
    !   Reecriture a partir d'un listing papier a Habas, le 14/02/00
    !
    !   le thermique est suppose homogene et dissipe par melange avec
    !   son environnement. la longueur l_mix controle l'efficacite du
    !   melange
    !
    !   Le calcul du transport des differentes especes se fait en prenant
    !   en compte:
    !     1. un flux de masse montant
    !     2. un flux de masse descendant
    !     3. un entrainement
    !     4. un detrainement
    !
    !=======================================================================

    !-----------------------------------------------------------------------
    !   declarations:
    !   -------------

    !#include "dimensions.h"
    !#include "YOMCST.h"
    !#include "YOETHF.h"
    !#include "FCTTRE.h"
    !#include "iniprint.h"
    !#include "thermcell.h"

    !   arguments:
    !   ----------

    !IM 140508
    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    INTEGER      , INTENT(IN   ) :: nClass
    INTEGER      , INTENT(IN   ) :: nAeros
    LOGICAL      , INTENT(IN   ) :: debut              !debut---input-L-variable logique indiquant le premier passage
    REAL(KIND=r8), INTENT(IN   ) :: ptimestep
    REAL(KIND=r8), INTENT(IN   ) :: tau_thermals
    REAL(KIND=r8), INTENT(IN   ) :: pplay  (ngrid,nlay)  ! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
    REAL(KIND=r8), INTENT(IN   ) :: pplev  (ngrid,nlay+1)! paprs---input-R-pression pour chaque inter-couche (en Pa)
    REAL(KIND=r8), INTENT(IN   ) :: pphi   (ngrid,nlay)   ! pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
    REAL(KIND=r8), INTENT(IN   ) :: pu     (ngrid,nlay)     ! u-------input-R-vitesse dans la direction X (de O a E) en m/s
    REAL(KIND=r8), INTENT(IN   ) :: pv     (ngrid,nlay)     ! v-------input-R-vitesse Y (de S a N) en m/s
    REAL(KIND=r8), INTENT(IN   ) :: pt     (ngrid,nlay)     ! t-------input-R-temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: po     (ngrid,nlay)     ! qx------input-R-humidite specifique (kg/kg) et d'autres traceurs
    REAL(KIND=r8), INTENT(INOUT) :: pice   (ngrid,nlay)                   ! qx------input (kg/kg) et d'autres traceurs
    REAL(KIND=r8), INTENT(INOUT) :: pliq   (ngrid,nlay)                   ! qx------input (kg/kg) et d'autres traceurs
    REAL(KIND=r8), INTENT(INOUT) :: pvar   (1:ngrid,1:nlay,1:nClass+nAeros)     ! qx------input (kg/kg) et d'autres traceurs
    REAL(KIND=r8), INTENT(OUT  ) :: pduadj (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: pdvadj (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: pdtadj (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: pdoadj (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: pdicedj(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: pdliqdj(ngrid,nlay) 
    REAL(KIND=r8), INTENT(OUT  ) :: pdaedj (ngrid,nlay,1:5+nClass+nAeros)
    REAL(KIND=r8), INTENT(INOUT) :: fm0    (ngrid,nlay+1)
    REAL(KIND=r8), INTENT(INOUT) :: entr0  (ngrid,nlay)
    REAL(KIND=r8), INTENT(INOUT) :: detr0  (ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pblh   (ngrid)  !INTEGER      , INTENT(IN   ) :: pblh(nCols)
    INTEGER(KIND=i8), INTENT(IN) :: imask  (ngrid)
    REAL(KIND=r8), INTENT(INOUT) :: dump   (ngrid,nlay)  
    REAL(KIND=r8), INTENT(OUT  ) :: zqta   (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqla   (ngrid,nlay)
    INTEGER      , INTENT(OUT  ) :: lmax   (ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: ratqscth(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ratqsdiff(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(ngrid,nlay) 
    REAL(KIND=r8), INTENT(OUT  ) :: Ale_bl (ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: Alp_bl (ngrid)
    INTEGER      , INTENT(OUT  ) :: lalim_conv(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: wght_th(ngrid,nlay)
    REAL(KIND=r8), INTENT(INOUT) :: zmax0  (ngrid)
    REAL(KIND=r8), INTENT(INOUT) :: f0     (ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: zw2    (ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: fraca  (ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: ztv    (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zpspsk (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ztla   (ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zthl   (ngrid,nlay)

    !   local:
    !   ------

    !      integer icount
    !      data icount/0/
    !      save icount
!!!!$OMP THREADPRIVATE(icount)

    INTEGER           :: igout=1
!!!!!!$OMP THREADPRIVATE(igout)
    INTEGER,PARAMETER :: lunout1=6
!!!!!$OMP THREADPRIVATE(lunout1)
!!!    INTEGER,PARAMETER :: lev_out=10
!!!!!$OMP THREADPRIVATE(lev_out)

    INTEGER       :: ig
    INTEGER       :: k,kk
    INTEGER       :: l
!    INTEGER       :: ll
    INTEGER       :: ierr
!    REAL(KIND=r8) :: zsortie1d(ngrid)
    INTEGER       :: lmin(ngrid)
    INTEGER       :: lalim(ngrid)
    INTEGER       :: lmix(ngrid)
    INTEGER       :: lmix_bis(ngrid)
    REAL(KIND=r8) :: linter(ngrid)
    REAL(KIND=r8) :: zmix(ngrid)
    REAL(KIND=r8) :: zmax(ngrid)
    REAL(KIND=r8) :: ztva(ngrid,nlay)
    REAL(KIND=r8) :: zw_est(ngrid,nlay+1)
    REAL(KIND=r8) :: ztva_est(ngrid,nlay)
    !      REAL(KIND=r8) fraca(ngrid,nlay)

    REAL(KIND=r8) :: zmax_sec(ngrid)
    !on garde le zmax du pas de temps precedent
    !FH/IM     save zmax0

    REAL(KIND=r8) :: lambda

    REAL(KIND=r8) :: zlev(ngrid,nlay+1)
    REAL(KIND=r8) :: zlay(ngrid,nlay)
    REAL(KIND=r8) :: deltaz(ngrid,nlay)
    REAL(KIND=r8) :: zh(ngrid,nlay)
    REAL(KIND=r8) :: zdthladj(ngrid,nlay)
    REAL(KIND=r8) :: zu(ngrid,nlay)
    REAL(KIND=r8) :: zv(ngrid,nlay)
    REAL(KIND=r8) :: zo(ngrid,nlay)
    REAL(KIND=r8) :: zl(ngrid,nlay)
!    REAL(KIND=r8) :: zsortie(ngrid,nlay)
    REAL(KIND=r8) :: zva(ngrid,nlay)
    REAL(KIND=r8) :: zua(ngrid,nlay)
    REAL(KIND=r8) :: zoa(ngrid,nlay)
    REAL(KIND=r8) :: zoaaux(ngrid,nlay)

    REAL(KIND=r8) :: zta(ngrid,nlay)
    REAL(KIND=r8) :: zha(ngrid,nlay)
    REAL(KIND=r8) :: zf
    REAL(KIND=r8) :: zf2
    REAL(KIND=r8) :: thetath2(ngrid,nlay)
    REAL(KIND=r8) :: wth2(ngrid,nlay)
    REAL(KIND=r8) :: wth3(ngrid,nlay)
    REAL(KIND=r8) :: q2(ngrid,nlay)
    ! FH probleme de dimensionnement avec l'allocation dynamique
    !     common/comtherm/thetath2,wth2
    REAL(KIND=r8) :: wq(ngrid,nlay)
    REAL(KIND=r8) :: wthl(ngrid,nlay)
    REAL(KIND=r8) :: wthv(ngrid,nlay)

    REAL(KIND=r8) :: var
    REAL(KIND=r8) :: vardiff

    LOGICAL       :: sorties
    REAL(KIND=r8) :: rho(ngrid,nlay)
    REAL(KIND=r8) :: rhobarz(ngrid,nlay)
    REAL(KIND=r8) :: masse(ngrid,nlay)

    REAL(KIND=r8) :: wmax(ngrid)
    REAL(KIND=r8) :: wmax_tmp(ngrid)
    REAL(KIND=r8) :: wmax_sec(ngrid)
    REAL(KIND=r8) :: fm(ngrid,nlay+1)
    REAL(KIND=r8) :: entr(ngrid,nlay)
    REAL(KIND=r8) :: detr(ngrid,nlay)

    !niveau de condensation
    INTEGER       :: nivcon(ngrid)
    REAL(KIND=r8) :: zcon(ngrid)
    REAL(KIND=r8) :: CHI
    REAL(KIND=r8) :: zcon2(ngrid)
    REAL(KIND=r8) :: pcon(ngrid)
    REAL(KIND=r8) :: zqsat(ngrid,nlay)

    REAL(KIND=r8) :: f_star(ngrid,nlay+1)
    REAL(KIND=r8) :: entr_star(ngrid,nlay)
    REAL(KIND=r8) :: detr_star(ngrid,nlay)
    REAL(KIND=r8) :: alim_star_tot(ngrid)
    REAL(KIND=r8) :: alim_star(ngrid,nlay)
    REAL(KIND=r8) :: alim_star_clos(ngrid,nlay)
    REAL(KIND=r8) :: f(ngrid)
    !FH/IM     save f0
!    REAL(KIND=r8) :: zlevinter(ngrid)
    REAL(KIND=r8) :: seuil
    REAL(KIND=r8) :: csc(ngrid,nlay)

    !
    !nouvelles variables pour la convection
    REAL(KIND=r8) :: alp_int(ngrid)
    REAL(KIND=r8) :: dp_int(ngrid)
    REAL(KIND=r8) :: zdp
!    REAL(KIND=r8) :: ale_int(ngrid)
!    INTEGER       :: n_int(ngrid)
    REAL(KIND=r8) :: fm_tot(ngrid)
    !v1d     logical therm
    !v1d     save therm

!    CHARACTER*2 :: str2
!    CHARACTER*10 :: str10

!    CHARACTER (len=20) :: modname='thermcell_main'
    CHARACTER (len=80) :: abort_message

    !      EXTERNAL SCOPY
    !
 pduadj=0.0_r8; pdvadj=0.0_r8;
 pdtadj=0.0_r8; pdoadj=0.0_r8;pdicedj =0.0_r8;
pdliqdj =0.0_r8;;pdaedj=0.0_r8
 zqta=0.0_r8;zqla=0.0_r8;lmax=0
 ratqscth=0.0_r8;ratqsdiff=0.0_r8;
zqsatth=0.0_r8; Ale_bl=0.0_r8
Alp_bl=0.0_r8;lalim_conv=0.0_r8
wght_th=0.0_r8;zw2=0.0_r8;
fraca=0.0_r8;ztv=0.0_r8;
zpspsk=0.0_r8;ztla=0.0_r8;
zthl=0.0_r8;lmin=0.0_r8
lalim=0.0_r8;lmix=0.0_r8
lmix_bis=0.0_r8;linter=0.0_r8
zmix=0.0_r8;zmax=0.0_r8
ztva=0.0_r8;zw_est=0.0_r8;
ztva_est=0.0_r8;zmax_sec=0.0_r8; lambda=0.0_r8

 zlev=0.0_r8; zlay=0.0_r8;
 deltaz=0.0_r8; zh=0.0_r8;
 zdthladj=0.0_r8; zu=0.0_r8;
 zv=0.0_r8; zo=0.0_r8;
 zl=0.0_r8; zva=0.0_r8;
 zua=0.0_r8; zoa=0.0_r8;

 zta=0.0_r8; zha=0.0_r8;
 zf=0.0_r8; zf2=0.0_r8;
 thetath2=0.0_r8; wth2=0.0_r8;
 wth3=0.0_r8; q2=0.0_r8;
 wq=0.0_r8; wthl=0.0_r8;
 wthv=0.0_r8; var=0.0_r8;
 vardiff=0.0_r8; rho=0.0_r8;
 rhobarz=0.0_r8; masse=0.0_r8;

 wmax=0.0_r8; wmax_tmp=0.0_r8
 wmax_sec=0.0_r8; fm=0.0_r8;
 entr=0.0_r8; detr=0.0_r8;

 nivcon=0.0_r8; zcon=0.0_r8
 CHI=0.0_r8; zcon2=0.0_r8
 pcon=0.0_r8; zqsat=0.0_r8;

 f_star=0.0_r8; entr_star=0.0_r8;
 detr_star=0.0_r8; alim_star_tot=0.0_r8
 alim_star=0.0_r8; alim_star_clos=0.0_r8;
 f=0.0_r8; seuil=0.0_r8;
 csc=0.0_r8; alp_int=0.0_r8
 dp_int=0.0_r8; zdp=0.0_r8;
 fm_tot=0.0_r8
    !-----------------------------------------------------------------------
    !   initialisation:
    !   ---------------
    !

    seuil=0.25_r8

    IF (debut)  THEN
       fm0=0.0_r8
       entr0=0.0_r8
       detr0=0.0_r8


       !#undef wrgrads_thermcell
       !#ifdef wrgrads_thermcell
       !! Initialisation des sorties grads pour les thermiques.
       !! Pour l'instant en 1D sur le point igout.
       !! Utilise par thermcell_out3d.h
       !         str10='therm'
       !         call inigrads(1,1,rlond(igout),1.,-180.,180.,jjm, &
       !     &   rlatd(igout),-90.,90.,1.,llm,pplay(igout,:),1.,   &
       !     &   ptimestep,str10,'therm ')
       !#endif



    ENDIF

    fm=0.0_r8 ; entr=0.0_r8 ; detr=0.0_r8


    !      icount=icount+1

    !IM 090508 beg
    !print*,'====================================================================='
    !print*,'====================================================================='
    !print*,' PAS ',icount,' PAS ',icount,' PAS ',icount,' PAS ',icount
    !print*,'====================================================================='
    !print*,'====================================================================='
    !IM 090508 end

    !IF (prt_level.GE.1) PRINT*,'thermcell_main V4'

    sorties=.TRUE.
!    IF(ngrid.NE.ngrid) THEN
!       PRINT*
!       PRINT*,'STOP dans convadj'
!       PRINT*,'ngrid    =',ngrid
!       PRINT*,'ngrid  =',ngrid
!    ENDIF
    !
    !     write(lunout,*)'WARNING thermcell_main f0=max(f0,1.e-2)'
    DO ig=1,ngrid
       f0(ig)=MAX(f0(ig),1.e-2_r8)
       zmax0(ig)=MAX(zmax0(ig),40.0_r8)
       !IMmarche pas ?!       if (f0(ig)<1.e-2_r8) f0(ig)=1.e-2_r8
    ENDDO

!    IF (prt_level.GE.20) THEN
!       DO ig=1,ngrid
!          PRINT*,'th_main ig f0',ig,f0(ig)
!       ENDDO
!    ENDIF
    !-----------------------------------------------------------------------
    ! Calcul de T,q,ql a partir de Tl et qT dans l environnement
    !   --------------------------------------------------------------------
    !

    CALL thermcell_env( &
         ngrid                    ,&! INTEGER      , INTENT(IN   ) :: ngrid
         nlay                     ,&! INTEGER      , INTENT(IN   ) :: nlay
         po     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: po(1:ngrid,1:nlay)
         pt     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: pt(1:ngrid,1:nlay)
         pu     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: pu(1:ngrid,1:nlay)
         pv     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: pv(1:ngrid,1:nlay)
         pplay  (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: pplay (1:ngrid,1:nlay)
         pplev  (1:ngrid,1:nlay+1),&! REAL(KIND=r8), INTENT(IN   ) :: pplev (1:ngrid,1:nlay+1)
         zo     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zo(1:ngrid,1:nlay)
         zh     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zh(1:ngrid,1:nlay)
         zl     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zl(1:ngrid,1:nlay)
         ztv    (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztv(1:ngrid,1:nlay)
         zthl   (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zthl  (1:ngrid,1:nlay)
         zu     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zu(1:ngrid,1:nlay)
         zv     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zv(1:ngrid,1:nlay)
         zpspsk (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zpspsk(1:ngrid,1:nlay) 
         zqsat  (1:ngrid,1:nlay)   )! REAL(KIND=r8), INTENT(OUT  ) :: pqsat (1:ngrid,1:nlay)

!    DO ig=1,ngrid
!
!       dump(ig,1)= zqsat  (ig,3)
!       dump(ig,2)= entr0(ig,1)
!       dump(ig,3)= detr0(ig,1)
!       dump(ig,4)= masse(ig,1)
!       dump(ig,5)= zthl (ig,1)
!       dump(ig,6)= fraca(ig,1)
!       dump(ig,7)= entr (ig,1)
!       dump(ig,8)= detr (ig,1)
!
!    END DO

!    IF (prt_level.GE.1) PRINT*,'thermcell_main apres thermcell_env'

    !------------------------------------------------------------------------
    !                       --------------------
    !
    !
    !                       + + + + + + + + + + +
    !
    !
    !  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
    !  wh,wt,wo ...
    !
    !                       + + + + + + + + + + +  zh,zu,zv,zo,rho
    !
    !
    !                       --------------------   zlev(1)
    !                       \\\\\\\\\\\\\\\\\\\\
    !
    !

    !-----------------------------------------------------------------------
    !   Calcul des altitudes des couches
    !-----------------------------------------------------------------------

    DO l=2,nlay
       DO ig=1,ngrid
          zlev(ig,l)=0.5_r8*(pphi(ig,l)+pphi(ig,l-1))/RG
       END DO
    ENDDO
    DO ig=1,ngrid
       zlev(ig,1)=0.0_r8
       zlev(ig,nlay+1)=(2.0_r8*pphi(ig,nlay)-pphi(ig,nlay-1))/RG
    END DO
    DO l=1,nlay
       DO ig=1,ngrid
          zlay(ig,l)=pphi(ig,l)/RG
       END DO
    ENDDO
    !calcul de l epaisseur des couches
    DO l=1,nlay
       DO ig=1,ngrid
          deltaz(ig,l)=zlev(ig,l+1)-zlev(ig,l)
       END DO
    ENDDO

    !     print*,'2 OK convect8'
    !-----------------------------------------------------------------------
    !   Calcul des densites
    !-----------------------------------------------------------------------
    DO l=1,nlay
       DO ig=1,ngrid
          rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*ztv(ig,l))
       END DO
    END DO

!    IF (prt_level.GE.10)WRITE(lunout1,*)                                &
!         &    'WARNING thermcell_main rhobarz(:,1)=rho(:,1)'
    DO ig=1,ngrid
       rhobarz(ig,1)=rho(ig,1)
    END DO

    DO l=2,nlay
       DO ig=1,ngrid
          rhobarz(ig,l)=0.5_r8*(rho(ig,l)+rho(ig,l-1))
       END DO 
    ENDDO

    !calcul de la masse
    DO l=1,nlay
       DO ig=1,ngrid
          masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
       END DO
    ENDDO

!    IF (prt_level.GE.1) PRINT*,'thermcell_main apres initialisation'

    !------------------------------------------------------------------
    !
    !             /|\
    !    --------  |  F_k+1 -------   
    !                              ----> D_k
    !             /|\              <---- E_k , A_k
    !    --------  |  F_k --------- 
    !                              ----> D_k-1
    !                              <---- E_k-1 , A_k-1
    !
    !
    !
    !
    !
    !    ---------------------------
    !
    !    ----- F_lmax+1=0 ----------         \
    !            lmax     (zmax)              |
    !    ---------------------------          |
    !                                         |
    !    ---------------------------          |
    !                                         |
    !    ---------------------------          |
    !                                         |
    !    ---------------------------          |
    !                                         |
    !    ---------------------------          |
    !                                         |  E
    !    ---------------------------          |  D
    !                                         |
    !    ---------------------------          |
    !                                         |
    !    ---------------------------  \       |
    !            lalim                 |      |
    !    ---------------------------   |      |
    !                                  |      |
    !    ---------------------------   |      |
    !                                  | A    |
    !    ---------------------------   |      |
    !                                  |      |
    !    ---------------------------   |      |
    !    lmin  (=1 pour le moment)     |      |
    !    ----- F_lmin=0 ------------  /      /
    !
    !    ---------------------------
    !    //////////////////////////
    !
    !
    !=============================================================================
    !  Calculs initiaux ne faisant pas intervenir les changements de phase
    !=============================================================================

    !------------------------------------------------------------------
    !  1. alim_star est le profil vertical de l'alimentation a la base du
    !     panache thermique, calcule a partir de la flotabilite de l'air sec
    !  2. lmin et lalim sont les indices inferieurs et superieurs de alim_star
    !------------------------------------------------------------------
    !
    entr_star=0.0 ; detr_star=0.0 ; alim_star=0.0 ; alim_star_tot=0.0
    lmin=1

    !-----------------------------------------------------------------------------
    !  3. wmax_sec et zmax_sec sont les vitesses et altitudes maximum d'un
    !     panache sec conservatif (e=d=0) alimente selon alim_star 
    !     Il s'agit d'un calcul de type CAPE
    !     zmax_sec est utilise pour determiner la geometrie du thermique.
    !------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !calcul du melange et des variables dans le thermique
    !--------------------------------------------------------------------------------
    !
    !IF (prt_level.GE.1) PRINT*,'avant thermcell_plume ',lev_out
    !IM 140508   CALL thermcell_plume(ngrid,nlay,ptimestep,ztv,zthl,po,zl,rhobarz,  &

    ! Gestion temporaire de plusieurs appels à thermcell_plume au travers
    ! de la variable iflag_thermals

    !      print*,'THERM thermcell_main iflag_thermals_ed=',iflag_thermals_ed
    IF (iflag_thermals_ed<=9) THEN
       !         print*,'THERM NOUVELLE/NOUVELLE Arnaud'
       CALL thermcell_plume(&
            ngrid                          ,&! INTEGER      , INTENT(IN   ) :: ngrid
            nlay                           ,&! INTEGER      , INTENT(IN   ) :: nlay
!           ptimestep                      ,&! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
            ztv          (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: ztv            (1:ngrid,1:nlay)
            zthl         (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: zthl           (1:ngrid,1:nlay)
            po           (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: po             (1:ngrid,1:nlay)
            zl           (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: zl             (1:ngrid,1:nlay)
            rhobarz      (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: rhobarz        (1:ngrid,1:nlay)
            zlev         (1:ngrid,1:nlay+1),&! REAL(KIND=r8), INTENT(IN   ) :: zlev           (1:ngrid,1:nlay+1)
            pplev        (1:ngrid,1:nlay+1),&! REAL(KIND=r8), INTENT(IN   ) :: pplev          (1:ngrid,1:nlay+1)
            pphi         (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: pphi           (1:ngrid,1:nlay)
            zpspsk       (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(IN   ) :: zpspsk         (1:ngrid,1:nlay)
            alim_star    (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: alim_star      (1:ngrid,1:nlay)
            alim_star_tot(1:ngrid)         ,&! REAL(KIND=r8), INTENT(OUT  ) :: alim_star_tot  (1:ngrid)
            lalim        (1:ngrid)         ,&! INTEGER      , INTENT(OUT  ) :: lalim          (1:ngrid)
            f0           (1:ngrid)         ,&! REAL(KIND=r8), INTENT(IN   ) :: f0             (1:ngrid)
            detr_star    (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: detr_star      (1:ngrid,1:nlay)
            entr_star    (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: entr_star      (1:ngrid,1:nlay)
            f_star       (1:ngrid,1:nlay+1),&! REAL(KIND=r8), INTENT(OUT  ) :: f_star         (1:ngrid,1:nlay+1)
            csc          (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: csc            (1:ngrid,1:nlay)
            ztva         (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztva           (1:ngrid,1:nlay)
            ztla         (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztla           (1:ngrid,1:nlay)
            zqla         (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zqla           (1:ngrid,1:nlay)
            zqta         (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zqta           (1:ngrid,1:nlay)
            zha          (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zha            (1:ngrid,1:nlay)
            zw2          (1:ngrid,1:nlay+1),&! REAL(KIND=r8), INTENT(OUT  ) :: zw2            (1:ngrid,1:nlay+1)
            zw_est       (1:ngrid,1:nlay+1),&! REAL(KIND=r8), INTENT(OUT  ) :: w_est          (1:ngrid,1:nlay+1)
            ztva_est     (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztva_est       (1:ngrid,1:nlay)
            zqsatth      (1:ngrid,1:nlay)  ,&! REAL(KIND=r8), INTENT(OUT  ) :: zqsatth        (1:ngrid,1:nlay)
            lmix         (1:ngrid)         ,&! INTEGER      , INTENT(OUT  ) :: lmix           (1:ngrid)
            lmix_bis     (1:ngrid)         ,&! INTEGER      , INTENT(OUT  ) :: lmix_bis       (1:ngrid)
            linter       (1:ngrid)          )! REAL(KIND=r8), INTENT(OUT  ) :: linter         (1:ngrid)

    ELSEIF (iflag_thermals_ed>9) THEN
       !        print*,'THERM RIO et al 2010, version d Arnaud'
       CALL thermcellV1_plume( &
            ngrid                           ,&!INTEGER      , INTENT(IN   ) :: ngrid
            nlay                            ,&!INTEGER      , INTENT(IN   ) :: nlay
!           ptimestep                       ,&!REAL(KIND=r8), INTENT(IN   ) :: ptimestep
            ztv           (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: ztv           (1:ngrid,1:nlay)
            zthl          (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: zthl          (1:ngrid,1:nlay)
            po            (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: po            (1:ngrid,1:nlay)
            zl            (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: zl            (1:ngrid,1:nlay)
            rhobarz       (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: rhobarz       (1:ngrid,1:nlay)
            zlev          (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(IN   ) :: zlev          (1:ngrid,1:nlay+1)
            pplev         (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(IN   ) :: pplev         (1:ngrid,1:nlay+1)
            pphi          (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: pphi          (1:ngrid,1:nlay)
            zpspsk        (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(IN   ) :: zpspsk        (1:ngrid,1:nlay)
            alim_star     (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: alim_star     (1:ngrid,1:nlay)
            alim_star_tot (1:ngrid)         ,&!REAL(KIND=r8), INTENT(OUT  ) :: alim_star_tot (1:ngrid)
            lalim         (1:ngrid)         ,&!INTEGER      , INTENT(OUT  ) :: lalim         (1:ngrid)
            f0            (1:ngrid)         ,&!REAL(KIND=r8), INTENT(IN   ) :: f0            (1:ngrid)
            detr_star     (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: detr_star     (1:ngrid,1:nlay)
            entr_star     (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: entr_star     (1:ngrid,1:nlay)
            f_star        (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(OUT  ) :: f_star        (1:ngrid,1:nlay+1)
            csc           (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: csc           (1:ngrid,1:nlay)
            ztva          (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztva          (1:ngrid,1:nlay)
            ztla          (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztla          (1:ngrid,1:nlay)
            zqla          (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqla          (1:ngrid,1:nlay)
            zqta          (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqta          (1:ngrid,1:nlay)
            zha           (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zha           (1:ngrid,1:nlay)
            zw2           (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(OUT  ) :: zw2           (1:ngrid,1:nlay+1)
            zw_est        (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(OUT  ) :: w_est         (1:ngrid,1:nlay+1)
            ztva_est      (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztva_est      (1:ngrid,1:nlay)
            zqsatth       (1:ngrid,1:nlay)  ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqsatth       (1:ngrid,1:nlay)
            lmix          (1:ngrid)         ,&!INTEGER      , INTENT(OUT  ) :: lmix          (1:ngrid)
            lmix_bis      (1:ngrid)         ,&!INTEGER      , INTENT(OUT  ) :: lmix_bis      (1:ngrid)
            linter        (1:ngrid)          )!REAL(KIND=r8), INTENT(OUT  ) :: linter        (1:ngrid)

    ENDIF

!    DO ig=1,ngrid
!
!       dump(ig,1)= zqsat  (ig,3)
!       dump(ig,2)= detr_star(ig,3)
!       dump(ig,3)= entr_star(ig,3)
!       dump(ig,4)= f_star(ig,3)
!       dump(ig,5)= zthl (ig,1)
!       dump(ig,6)= fraca(ig,1)
!       dump(ig,7)= entr (ig,1)
!       dump(ig,8)= detr (ig,1)
!
!    END DO

!RETURN !pkubota

    !IF (prt_level.GE.1) PRINT*,'apres thermcell_plume ',lev_out

    CALL test_ltherm( &
         ngrid                        ,&! integer      , INTENT(IN   ) :: ngrid
         nlay                         ,&! integer      , INTENT(IN   ) :: nlay
         pplay      (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay  (1:ngrid,1:nlay)
         lalim      (1:ngrid)         ,&! integer      , INTENT(IN   ) :: long   (1:ngrid)
         ztv        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv    (1:ngrid,1:nlay)
         po         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po     (1:ngrid,1:nlay)
         ztva       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva   (1:ngrid,1:nlay)
         zqla       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla   (1:ngrid,1:nlay)
         f_star     (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star (1:ngrid,1:nlay)
         zw2        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2    (1:ngrid,1:nlay)
         'thermcell_plum lalim '       )! character(LEN=*), INTENT(IN) :: comment
    CALL test_ltherm(&
         ngrid                          ,&! integer      , INTENT(IN   ) :: ngrid
         nlay                           ,&! integer      , INTENT(IN   ) :: nlay
         pplay        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay  (1:ngrid,1:nlay)
         lmix         (1:ngrid)         ,&! integer      , INTENT(IN   ) :: long   (1:ngrid)
         ztv          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv    (1:ngrid,1:nlay)
         po           (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po     (1:ngrid,1:nlay)
         ztva         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva   (1:ngrid,1:nlay)
         zqla         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla   (1:ngrid,1:nlay)
         f_star       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star (1:ngrid,1:nlay)
         zw2          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2    (1:ngrid,1:nlay)
         'thermcell_plum lmix  '         )! character(LEN=*), INTENT(IN) :: comment

    !IF (prt_level.GE.1) PRINT*,'thermcell_main apres thermcell_plume'
 !   IF (prt_level.GE.10) THEN
 !      WRITE(lunout1,*) 'Dans thermcell_main 2'
 !      WRITE(lunout1,*) 'lmin ',lmin(igout)
 !      WRITE(lunout1,*) 'lalim ',lalim(igout)
 !      WRITE(lunout1,*) ' ig l alim_star entr_star detr_star f_star '
 !      WRITE(lunout1,'(i6,i4,4e15.5)') (igout,l,alim_star(igout,l),entr_star(igout,l),detr_star(igout,l) &
 !           &    ,f_star(igout,l+1),l=1,NINT(linter(igout))+5)
 !   ENDIF

    !-------------------------------------------------------------------------------
    ! Calcul des caracteristiques du thermique:zmax,zmix,wmax
    !-------------------------------------------------------------------------------
    !
    CALL thermcell_height( &
         ngrid                   ,&!INTEGER      , INTENT(IN   ) :: ngrid
         nlay                    ,&!INTEGER      , INTENT(IN   ) :: nlay
         lalim (1:ngrid)         ,&!INTEGER      , INTENT(INOUT) :: lalim  (1:ngrid)
         lmin  (1:ngrid)         ,&!INTEGER      , INTENT(INOUT) :: lmin   (1:ngrid)
         linter(1:ngrid)         ,&!REAL(KIND=r8), INTENT(IN   ) :: linter (1:ngrid)
         lmix  (1:ngrid)         ,&!INTEGER      , INTENT(INOUT) :: lmix   (1:ngrid)
         zw2   (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(INOUT) :: zw2    (1:ngrid,1:nlay+1)
         zlev  (1:ngrid,1:nlay+1),&!REAL(KIND=r8), INTENT(IN   ) :: zlev   (1:ngrid,1:nlay+1)
         lmax  (1:ngrid)         ,&!INTEGER      , INTENT(OUT  ) :: lmax   (1:ngrid)
         zmax  (1:ngrid)         ,&!REAL(KIND=r8), INTENT(OUT  ) :: zmax   (1:ngrid)
         zmax0 (1:ngrid)         ,&!REAL(KIND=r8), INTENT(INOUT) :: zmax0  (1:ngrid)
         zmix  (1:ngrid)         ,&!REAL(KIND=r8), INTENT(OUT  ) :: zmix   (1:ngrid)
         wmax  (1:ngrid)          )!REAL(KIND=r8), INTENT(OUT  ) :: wmax   (1:ngrid)
!         lev_out )!INTEGER      , INTENT(IN   ) :: lev_out
    ! Attention, w2 est transforme en sa racine carree dans cette routine
    ! Le probleme vient du fait que linter et lmix sont souvent égaux à 1.
!    DO ig=1,ngrid
!
!       dump(ig,1)= zqsat  (ig,3)
!       dump(ig,2)= detr_star(ig,3)
!       dump(ig,3)= entr_star(ig,3)
!       dump(ig,4)= f_star(ig,3)
!       dump(ig,5)= zmix (ig)
!       dump(ig,6)= zw2(ig,3)
!       dump(ig,7)= entr (ig,1)
!       dump(ig,8)= detr (ig,1)
!
!    END DO

    wmax_tmp=0.0_r8
    DO l=1,nlay
       DO ig=1,ngrid
          wmax_tmp(ig)=MAX(wmax_tmp(ig),zw2(ig,l))
       END DO
    ENDDO
    !     print*,"ZMAX ",lalim,lmin,linter,lmix,lmax,zmax,zmax0,zmix,wmax

!RETURN!pkubota2

    CALL test_ltherm(&
         ngrid                         ,&! integer      , INTENT(IN   ) :: ngrid
         nlay                          ,&! integer      , INTENT(IN   ) :: nlay
         pplay       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay  (ngrid,nlay)
         lalim       (1:ngrid)         ,&! integer      , INTENT(IN   ) :: long   (ngrid)
         ztv         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv    (ngrid,nlay)
         po          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po          (ngrid,nlay)
         ztva        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva   (ngrid,nlay)
         zqla        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla   (ngrid,nlay)
         f_star      (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star (ngrid,nlay)
         zw2         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2    (ngrid,nlay)
         'thermcell_heig lalim ' )! character(LEN=*), INTENT(IN) :: comment
    CALL test_ltherm(&
         ngrid                          ,&! integer      , INTENT(IN   ) :: ngrid
         nlay                           ,&! integer      , INTENT(IN   ) :: nlay
         pplay        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
         lmin         (1:ngrid)         ,&! integer      , INTENT(IN   ) :: long(ngrid)
         ztv          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
         po           (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
         ztva         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva(ngrid,nlay)
         zqla         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla(ngrid,nlay)
         f_star       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star(ngrid,nlay)
         zw2          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay)
         'thermcell_heig lmin  ' )! character(LEN=*), INTENT(IN) :: comment
    CALL test_ltherm(&
         ngrid                           ,&! integer      , INTENT(IN        ) :: ngrid
         nlay                            ,&! integer      , INTENT(IN        ) :: nlay
         pplay         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: pplay(ngrid,nlay)
         lmix          (1:ngrid)         ,&! integer      , INTENT(IN        ) :: long(ngrid)
         ztv           (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: ztv(ngrid,nlay)
         po            (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: po(ngrid,nlay)
         ztva          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: ztva(ngrid,nlay)
         zqla          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: zqla(ngrid,nlay)
         f_star        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: f_star(ngrid,nlay)
         zw2           (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN        ) :: zw2(ngrid,nlay)
         'thermcell_heig lmix  ' ) ! character(LEN=*), INTENT(IN) :: comment
    CALL test_ltherm( &
         ngrid                         ,&! integer      , INTENT(IN   ) :: ngrid
         nlay                          ,&! integer      , INTENT(IN   ) :: nlay
         pplay       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
         lmax        (1:ngrid)         ,&! integer      , INTENT(IN   ) :: long(ngrid)
         ztv         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
         po          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
         ztva        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva(ngrid,nlay)
         zqla        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla(ngrid,nlay)
         f_star      (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star(ngrid,nlay)
         zw2         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay)
         'thermcell_heig lmax  ' )! character(LEN=*), INTENT(IN) :: comment

    !IF (prt_level.GE.1) PRINT*,'thermcell_main apres thermcell_height'

    !-------------------------------------------------------------------------------
    ! Fermeture,determination de f
    !-------------------------------------------------------------------------------
    !
    !
    !!      write(lunout,*)'THERM NOUVEAU XXXXX'
    CALL thermcell_dry( &
         ngrid                       , &! INTEGER      , INTENT(IN   ) :: ngrid
         nlay                        , &! INTEGER      , INTENT(IN   ) :: nlay
         zlev      (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
         pphi      (1:ngrid,1:nlay  ), &! REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
         ztv       (1:ngrid,1:nlay  ), &! REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
         alim_star (1:ngrid,1:nlay  ), &! REAL(KIND=r8), INTENT(IN   ) :: alim_star(ngrid,nlay)
         lalim     (1:ngrid)         , &! INTEGER      , INTENT(IN   ) :: lalim(ngrid)
         lmin      (1:ngrid)         , &! INTEGER      , INTENT(IN   ) :: lmin(ngrid)
         zmax_sec  (1:ngrid)         , &! REAL(KIND=r8), INTENT(OUT  ) :: zmax(ngrid)
         wmax_sec  (1:ngrid)           )! REAL(KIND=r8), INTENT(OUT  ) :: wmax(ngrid)
!         lev_out     )! INTEGER      , INTENT(IN   ) :: lev_out      ! niveau pour les print

!    DO ig=1,ngrid
!
!       dump(ig,1)= zqsat  (ig,3)
!       dump(ig,2)= detr_star(ig,3)
!       dump(ig,3)= entr_star(ig,3)
!       dump(ig,4)= f_star(ig,3)
!       dump(ig,5)= zmix (ig,3)
!       dump(ig,6)= zw2(ig,3)
!       dump(ig,7)= zmax_sec (ig)
!       dump(ig,8)= zmax_sec (ig)
!
!    END DO

    CALL test_ltherm(&
         ngrid                         ,&! integer            , INTENT(IN        ) :: ngrid
         nlay                          ,&! integer            , INTENT(IN        ) :: nlay
         pplay       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: pplay(ngrid,nlay)
         lmin        (1:ngrid)         ,&! integer            , INTENT(IN        ) :: long(ngrid)
         ztv         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: ztv(ngrid,nlay)
         po          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: po(ngrid,nlay)
         ztva        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: ztva(ngrid,nlay)
         zqla        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: zqla(ngrid,nlay)
         f_star      (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: f_star(ngrid,nlay)
         zw2         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN           ) :: zw2(ngrid,nlay)
         'thermcell_dry  lmin  ' )! character(LEN=*), INTENT(IN) :: comment
    CALL test_ltherm( &
         ngrid                          ,&! integer      , INTENT(IN        ) :: ngrid
         nlay                           ,&! integer      , INTENT(IN        ) :: nlay
         pplay        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: pplay(ngrid,nlay)
         lalim        (1:ngrid)         ,&! integer             , INTENT(IN        ) :: long(ngrid)
         ztv          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: ztv(ngrid,nlay)
         po           (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: po(ngrid,nlay)
         ztva         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: ztva(ngrid,nlay)
         zqla         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: zqla(ngrid,nlay)
         f_star       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: f_star(ngrid,nlay)
         zw2          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN            ) :: zw2(ngrid,nlay)
         'thermcell_dry  lalim ' )! character(LEN=*), INTENT(IN) :: comment

    !IF (prt_level.GE.1) PRINT*,'thermcell_main apres thermcell_dry'
!    IF (prt_level.GE.10) THEN
!       WRITE(lunout1,*) 'Dans thermcell_main 1b'
!       WRITE(lunout1,*) 'lmin ',lmin(igout)
!       WRITE(lunout1,*) 'lalim ',lalim(igout)
!       WRITE(lunout1,*) ' ig l alim_star entr_star detr_star f_star '
!       WRITE(lunout1,'(i6,i4,e15.5)') (igout,l,alim_star(igout,l) &
!            &    ,l=1,lalim(igout)+4)
!    ENDIF




    ! Choix de la fonction d'alimentation utilisee pour la fermeture.
    ! Apparemment sans importance
    DO l=1,nlay
       DO ig=1,ngrid
          alim_star_clos(ig,l)=alim_star(ig,l)
          alim_star_clos(ig,l)=entr_star(ig,l)+alim_star(ig,l)
       END DO
    END DO 

    ! Appel avec la version seche
    CALL thermcell_closure( &
         ngrid                               ,&!  INTEGER      , INTENT(IN ) :: ngrid
         nlay                                ,&!  INTEGER      , INTENT(IN ) :: nlay
         r_aspect_thermals                   ,&!  REAL(KIND=r8), INTENT(IN ) :: r_aspect
         rho               (1:ngrid,1:nlay)  ,&!  REAL(KIND=r8), INTENT(IN ) :: rho(ngrid,nlay)
         zlev              (1:ngrid,1:nlay)  ,&!  REAL(KIND=r8), INTENT(IN ) :: zlev(ngrid,nlay)
         lalim             (1:ngrid)         ,&!  INTEGER      , INTENT(IN ) :: lalim(ngrid)
         alim_star_clos    (1:ngrid,1:nlay)  ,&!  REAL(KIND=r8), INTENT(IN ) :: alim_star(ngrid,nlay)
         zmax_sec          (1:ngrid)         ,&!  REAL(KIND=r8), INTENT(IN ) :: zmax(ngrid)
         wmax_sec          (1:ngrid)         ,&!  REAL(KIND=r8), INTENT(IN ) :: wmax(ngrid)
         f                 (1:ngrid)          )!  REAL(KIND=r8), INTENT(OUT) :: f(ngrid)

!    DO ig=1,ngrid
!
!       dump(ig,1)= zqsat  (ig,3)
!       dump(ig,2)= detr_star(ig,3)
!       dump(ig,3)= entr_star(ig,3)
!       dump(ig,4)= f_star(ig,3)
!       dump(ig,5)= zmix (ig,3)
!       dump(ig,6)= zw2(ig,3)
!       dump(ig,7)= zmax_sec (ig,1)
!       dump(ig,8)= zmax_sec (ig,1)
!       dump(ig,9)= f (ig)

!
!    END DO

!RETURN !pkubota3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Appel avec les zmax et wmax tenant compte de la condensation
    ! Semble moins bien marcher
    !     CALL thermcell_closure(ngrid,nlay,r_aspect_thermals,ptimestep,rho,  &
    !    &   zlev,lalim,alim_star,f_star,zmax,wmax,f,lev_out)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !IF(prt_level.GE.1)PRINT*,'thermcell_closure apres thermcell_closure'
    DO ig=1,ngrid
       IF (tau_thermals>1.0_r8) THEN
          lambda=EXP(-ptimestep/tau_thermals)
          lambda=0.7_r8
          f0(ig)=1.00_r8*((1.0_r8-lambda)*f(ig)+lambda*f0(ig))
       ELSE
          lambda=1.0_r8
          f0(ig)=f(ig)*lambda
       ENDIF
   END DO

    ! Test valable seulement en 1D mais pas genant
!    IF (.NOT. (f0(1).GE.0.0_r8) ) THEN
!       abort_message = '.not. (f0(1).ge.0.)'
!       STOP
!       !
!       !              CALL abort_gcm (modname,abort_message,1)
!    ENDIF

    !-------------------------------------------------------------------------------
    !deduction des flux
    !-------------------------------------------------------------------------------

    CALL thermcell_flux2( &
         ngrid                       , &! INTEGER        , INTENT(IN   ) :: ngrid
         nlay                        , &! INTEGER        , INTENT(IN   ) :: nlay
         ptimestep                   , &! REAL(KIND=r8)  , INTENT(IN   ) :: ptimestep
         masse     (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(IN   ) :: masse     (1:ngrid,1:nlay)
         lalim     (1:ngrid)         , &! INTEGER        , INTENT(IN   ) :: lalim     (1:ngrid)
         lmax      (1:ngrid)         , &! INTEGER        , INTENT(IN   ) :: lmax      (1:ngrid)
         alim_star (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(IN   ) :: alim_star (1:ngrid,1:nlay)
         entr_star (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(IN   ) :: entr_star (1:ngrid,1:nlay)
         detr_star (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(IN   ) :: detr_star (1:ngrid,1:nlay)
         f         (1:ngrid)         , &! REAL(KIND=r8)  , INTENT(IN   ) :: f         (1:ngrid)
         rhobarz   (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(IN   ) :: rhobarz   (1:ngrid,1:nlay)
         zw2       (1:ngrid,1:nlay+1), &! REAL(KIND=r8)  , INTENT(OUT  ) :: zw2       (1:ngrid,1:nlay+1)
         fm        (1:ngrid,1:nlay+1), &! REAL(KIND=r8)  , INTENT(OUT  ) :: fm        (1:ngrid,1:nlay+1)
         entr      (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(OUT  ) :: entr      (1:ngrid,1:nlay)
         detr      (1:ngrid,1:nlay)  , &! REAL(KIND=r8)  , INTENT(OUT  ) :: detra     (1:ngrid,1:nlay)
         lunout1                     , &! integer        , INTENT(IN   ) :: lunout1
         igout                         )! integer        , INTENT(INOUT) :: igout

    !IM 060508    &       detr,zqla,zmax,lev_out,lunout,igout)
!    DO ig=1,ngrid
!
!       dump(ig,1)= zqsat  (ig,3)
!       dump(ig,2)= detr_star(ig,3)
!       dump(ig,3)= entr_star(ig,3)
!       dump(ig,4)= f_star(ig,3)
!       dump(ig,5)= zmix (ig,3)
!       dump(ig,6)= zw2(ig,3)
!       dump(ig,7)= zmax_sec (ig,1)
!       dump(ig,8)= zmax_sec (ig,1)
!       dump(ig,9)= f (ig,3)
!       dump(ig,10)= zw2 (ig,3)
!       dump(ig,11)= fm (ig,3)
!       dump(ig,12)= entr (ig,3)
!       dump(ig,13)= detr (ig,3)

!
!    END DO

    !IF (prt_level.GE.1) PRINT*,'thermcell_main apres thermcell_flux'
    CALL test_ltherm(&
         ngrid                         ,&! integer        , INTENT(IN   ) :: ngrid
         nlay                          ,&! integer        , INTENT(IN   ) :: nlay
         pplay       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
         lalim       (1:ngrid)         ,&! integer        , INTENT(IN   ) :: long(ngrid)
         ztv         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
         po          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
         ztva        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva(ngrid,nlay)
         zqla        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla(ngrid,nlay)
         f_star      (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star(ngrid,nlay)
         zw2         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay)
         'thermcell_flux lalim ')! character(LEN=*), INTENT(IN) :: comment
    CALL test_ltherm(&
         ngrid                         ,&! integer        , INTENT(IN   ) :: ngrid
         nlay                          ,&! integer        , INTENT(IN   ) :: nlay
         pplay       (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
         lmax        (1:ngrid)         ,&! integer        , INTENT(IN   ) :: long(ngrid)
         ztv         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
         po          (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
         ztva        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: ztva(ngrid,nlay)
         zqla        (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zqla(ngrid,nlay)
         f_star      (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: f_star(ngrid,nlay)
         zw2         (1:ngrid,1:nlay)  ,&! real(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay)
         'thermcell_flux lmax  ' )! character(LEN=*), INTENT(IN) :: comment

    !------------------------------------------------------------------
    !   On ne prend pas directement les profils issus des calculs precedents
    !   mais on s'autorise genereusement une relaxation vers ceci avec
    !   une constante de temps tau_thermals (typiquement 1800s).
    !------------------------------------------------------------------
    DO ig=1,ngrid
       IF (tau_thermals>1.0_r8) THEN
          lambda=EXP(-ptimestep/tau_thermals)
          lambda=0.7_r8
          fm0 (ig,1:nlay+1) =1.00_r8*((1.0_r8-lambda)*fm(ig,1:nlay+1)  +lambda*fm0(ig,1:nlay+1)  )
       ELSE
          lambda=1.0_r8
          fm0(ig,1:nlay+1)=fm(ig,1:nlay+1)*lambda
       END IF
    END DO
    DO l=1,nlay
       DO ig=1,ngrid
          IF (tau_thermals>1.0_r8) THEN
             lambda=EXP(-ptimestep/tau_thermals)
             lambda=0.7_r8
             fm0 (ig,l) =1.00_r8*((1.0_r8-lambda)*  fm(ig,l)  +lambda*fm0(ig,l))
             entr0(ig,l)=1.00_r8*((1.0_r8-lambda)*entr(ig,l)+lambda*entr0(ig,l))
             detr0(ig,l)=1.00_r8*((1.0_r8-lambda)*detr(ig,l)+lambda*detr0(ig,l))
          ELSE
             lambda=1.0_r8
             fm0(ig,l)=fm(ig,l)*lambda
             entr0(ig,l)=entr(ig,l)*lambda
             detr0(ig,l)=detr(ig,l)*lambda
          ENDIF
       END DO
    END DO
!    DO ig=1,ngrid
!
!       dump(ig,14)= fm0  (ig,1)
!       dump(ig,15)= entr0(ig,1)
!       dump(ig,16)= detr0(ig,1)
!       dump(ig,17)= masse(ig,1)
!       dump(ig,18)= zthl (ig,1)
!       dump(ig,19)= fraca(ig,1)
!       dump(ig,20)= entr (ig,1)
!       dump(ig,21)= detr (ig,1)
!
!    END DO
    !c------------------------------------------------------------------
    !   calcul du transport vertical
    !------------------------------------------------------------------

    CALL thermcell_dq( &
         ngrid                         , &! integer	 , INTENT(IN   ) :: ngrid
         nlay                          , &! integer	 , INTENT(IN   ) :: nlay
         ptimestep                     , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
         fm0         (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
         entr0       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
         masse       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
         zthl        (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
         zdthladj    (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
         imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
         dump        (1:ngrid,1:nlay)  , &
         zta         (1:ngrid,1:nlay)    )! REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)
!         lev_out        )! integer      ,
    CALL thermcell_dq(&
         ngrid                         , &! integer	 , INTENT(IN   ) :: ngrid
         nlay                          , &! integer	 , INTENT(IN   ) :: nlay
         ptimestep                     , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
         fm0         (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
         entr0       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
         masse       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
         po          (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
         pdoadj      (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
         imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
         dump        (1:ngrid,1:nlay)  , &
         zoa         (1:ngrid,1:nlay)    )! REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)
!         lev_out        )! integer      ,

    CALL thermcell_dq(&
         ngrid                         , &! integer	 , INTENT(IN   ) :: ngrid
         nlay                          , &! integer	 , INTENT(IN   ) :: nlay
         ptimestep                     , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
         fm0         (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
         entr0       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
         masse       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
         pice        (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
         pdicedj     (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
         imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
         dump        (1:ngrid,1:nlay)  , &
         zoaaux      (1:ngrid,1:nlay)    )! REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)
!         lev_out        )! integer      ,

    CALL thermcell_dq(&
         ngrid                         , &! integer	 , INTENT(IN   ) :: ngrid
         nlay                          , &! integer	 , INTENT(IN   ) :: nlay
         ptimestep                     , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
         fm0         (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
         entr0       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
         masse       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
         pliq        (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
         pdliqdj     (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
         imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
         dump        (1:ngrid,1:nlay)  , &
         zoaaux      (1:ngrid,1:nlay)    )! REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)
!         lev_out       )! integer      , 

    DO kk=1,nClass+nAeros
       CALL thermcell_dq(&
            ngrid                              , &! integer	 , INTENT(IN   ) :: ngrid
            nlay                               , &! integer	 , INTENT(IN   ) :: nlay
            ptimestep                          , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
            fm0      (1:ngrid,1:nlay+1    )    , &! REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
            entr0    (1:ngrid,1:nlay      )    , &! REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
            masse    (1:ngrid,1:nlay      )    , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
            pvar     (1:ngrid,1:nlay,kk   )    , &! REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
            pdaedj   (1:ngrid,1:nlay,5+kk)     , &! REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
            imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
            dump        (1:ngrid,1:nlay)  , &
            zoaaux   (1:ngrid,1:nlay)          )  ! REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)
!         lev_out        )! integer      , INTENT(
    END DO

!    DO ig=1,ngrid
!
!       dump(ig,14)= fm0      (ig,3)
!       dump(ig,15)= entr0    (ig,3)
!       dump(ig,16)= detr0    (ig,3)
!       dump(ig,17)= masse    (ig,3)
!       dump(ig,18)= zdthladj (ig,3) 
!       dump(ig,19)= pdoadj   (ig,3) 
!       dump(ig,20)= pdicedj  (ig,3) 
!       dump(ig,21)= pdliqdj  (ig,3) 
!                                    
!    END DO
    !------------------------------------------------------------------
    ! Calcul de la fraction de l'ascendance
    !------------------------------------------------------------------
    DO ig=1,ngrid
       fraca(ig,1)=0.0_r8
       fraca(ig,nlay+1)=0.0_r8
    ENDDO
    DO l=2,nlay
       DO ig=1,ngrid
          IF (zw2(ig,l).GT.1.e-10_r8) THEN
             fraca(ig,l)=fm(ig,l)/(rhobarz(ig,l)*zw2(ig,l))
          ELSE
             fraca(ig,l)=0.0_r8
          ENDIF
       ENDDO
    ENDDO

    !------------------------------------------------------------------
    !  calcul du transport vertical du moment horizontal
    !------------------------------------------------------------------

    !IM 090508  
    IF (1.EQ.1) THEN
       !IM 070508 vers. _dq       
       !     if (1.eq.0) then


       ! Calcul du transport de V tenant compte d'echange par gradient
       ! de pression horizontal avec l'environnement

       CALL thermcell_dv2(&
            ngrid                         ,&!integer      , INTENT(IN    ) :: ngrid
            nlay                          ,&!integer      , INTENT(IN    ) :: nlay
            ptimestep                     ,&!real(KIND=r8), INTENT(IN    ) :: ptimestep
            zlay      (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(IN    ) :: zlay     (1:ngrid,1:nlay)
            pblh      (1:ngrid)           ,&!REAL(KIND=r8), INTENT(IN    ) :: pblh     (1:ngrid)
            imask     (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
            fm0       (1:ngrid,1:nlay+1)  ,&!real(KIND=r8), INTENT(IN    ) :: fm       (1:ngrid,1:nlay+1)
            entr0     (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(IN    ) :: entr     (1:ngrid,1:nlay)
            masse     (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(IN    ) :: masse    (1:ngrid,1:nlay)
            fraca     (1:ngrid,1:nlay+1)  ,&!real(KIND=r8), INTENT(IN    ) :: fraca    (1:ngrid,1:nlay+1)
            zmax      (1:ngrid)           ,&!real(KIND=r8), INTENT(IN    ) :: larga    (1:ngrid)
            zu        (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(IN    ) :: u        (1:ngrid,1:nlay)
            zv        (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(IN    ) :: v        (1:ngrid,1:nlay)
            pduadj    (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(OUT   ) :: du       (1:ngrid,1:nlay)
            pdvadj    (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(OUT   ) :: dv       (1:ngrid,1:nlay)
            zua       (1:ngrid,1:nlay)    ,&!real(KIND=r8), INTENT(OUT   ) :: ua       (1:ngrid,1:nlay)
            zva       (1:ngrid,1:nlay)     )!real(KIND=r8), INTENT(OUT   ) :: va       (1:ngrid,1:nlay)
!            lev_out     )!integer     , INTENT(IN     ) :: lev_out       ! niveau pour les print

    ELSE

       ! calcul purement conservatif pour le transport de V
       CALL thermcell_dq( &
            ngrid                       , &! integer      , INTENT(IN   ) :: ngrid
            nlay                        , &! integer      , INTENT(IN   ) :: nlay
            ptimestep                   , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
            fm0       (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: fm        (1:ngrid,1:nlay+1)
            entr0     (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(IN   ) :: entr      (1:ngrid,1:nlay)
            masse     (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: masse     (1:ngrid,1:nlay)
            zu        (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: q         (1:ngrid,1:nlay)
            pduadj    (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(OUT  ) :: dq        (1:ngrid,1:nlay)
            imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
            dump        (1:ngrid,1:nlay)  , &
            zua       (1:ngrid,1:nlay)    )! REAL(KIND=r8), INTENT(OUT  ) :: qa        (1:ngrid,1:nlay)
!            lev_out     )! integer      , INTENT(IN   ) :: lev_out
       CALL thermcell_dq( &
            ngrid                      , &! integer      , INTENT(IN   ) :: ngrid
            nlay                       , &! integer      , INTENT(IN   ) :: nlay
            ptimestep                  , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
            fm0      (1:ngrid,1:nlay+1), &! REAL(KIND=r8), INTENT(IN   ) :: fm       (1:ngrid,1:nlay+1)
            entr0    (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(IN   ) :: entr     (1:ngrid,1:nlay)
            masse    (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: masse    (1:ngrid,1:nlay)
            zv       (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(INOUT) :: q        (1:ngrid,1:nlay)
            pdvadj   (1:ngrid,1:nlay)  , &! REAL(KIND=r8), INTENT(OUT  ) :: dq       (1:ngrid,1:nlay)
            imask       (1:ngrid)         , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
            dump        (1:ngrid,1:nlay)  , &
            zva      (1:ngrid,1:nlay)    )! REAL(KIND=r8), INTENT(OUT  ) :: qa       (1:ngrid,1:nlay)
!            lev_out    )! integer      , INTENT(IN   ) :: lev_out
    ENDIF

    !     print*,'13 OK convect8'
    DO l=1,nlay
       DO ig=1,ngrid
          pdtadj(ig,l)=zdthladj(ig,l)*zpspsk(ig,l)  
       ENDDO
    ENDDO

    !IF (prt_level.GE.1) PRINT*,'14 OK convect8'
    !------------------------------------------------------------------
    !   Calculs de diagnostiques pour les sorties
    !------------------------------------------------------------------
    !calcul de fraca pour les sorties

    IF (sorties) THEN
       !IF (prt_level.GE.1) PRINT*,'14a OK convect8'
       ! calcul du niveau de condensation
       ! initialisation
       DO ig=1,ngrid
          nivcon(ig)=0
          zcon(ig)=0.0_r8
       ENDDO
       !nouveau calcul
       DO ig=1,ngrid
          CHI=zh(ig,1)/(1669.0_r8-122.0_r8*zo(ig,1)/zqsat(ig,1)-zh(ig,1))
          pcon(ig)=pplay(ig,1)*(zo(ig,1)/zqsat(ig,1))**CHI
       ENDDO
       !IM   do k=1,nlay
       DO k=1,nlay-1
          DO ig=1,ngrid
             IF ((pcon(ig).LE.pplay(ig,k))  &
                  &      .AND.(pcon(ig).GT.pplay(ig,k+1))) THEN
                zcon2(ig)=zlay(ig,k)-(pcon(ig)-pplay(ig,k))/(RG*rho(ig,k))/100.0_r8
             ENDIF
          ENDDO
       ENDDO
       !IM
       ierr=0
       DO ig=1,ngrid
          IF (pcon(ig).LE.pplay(ig,nlay)) THEN 
             zcon2(ig)=zlay(ig,nlay)-MAX(pcon(ig)-pplay(ig,nlay),0.01_r8)/(RG*rho(ig,nlay))/100.0_r8
             ierr=0!1
          ENDIF
       ENDDO
       IF (ierr==1) THEN
          abort_message = 'thermcellV0_main: les thermiques vont trop haut '
          STOP
          !           CALL abort_gcm (modname,abort_message,1)
       ENDIF

       !IF (prt_level.GE.1) PRINT*,'14b OK convect8'
       DO k=nlay,1,-1
          DO ig=1,ngrid
             IF (zqla(ig,k).GT.1e-10_r8) THEN
                nivcon(ig)=k
                zcon(ig)=zlev(ig,k)
             ENDIF
          ENDDO
       ENDDO
       !IF (prt_level.GE.1) PRINT*,'14c OK convect8'
       !calcul des moments
       !initialisation
       DO l=1,nlay
          DO ig=1,ngrid
             q2(ig,l)=0.0_r8
             wth2(ig,l)=0.0_r8
             wth3(ig,l)=0.0_r8
             ratqscth(ig,l)=0.0_r8
             ratqsdiff(ig,l)=0.0_r8
          ENDDO
       ENDDO
       !IF (prt_level.GE.1) PRINT*,'14d OK convect8'
       IF (prt_level.GE.10)WRITE(lunout1,*)                                &
            &     'WARNING thermcell_main wth2=0. si zw2 > 1.e-10'
       DO l=1,nlay
          DO ig=1,ngrid
             zf=fraca(ig,l)
             zf2=zf/(1.0_r8-zf)
             !
             thetath2(ig,l)=zf2*(ztla(ig,l)-zthl(ig,l))**2
             IF(zw2(ig,l).GT.1.e-10_r8) THEN
                wth2(ig,l)=zf2*(zw2(ig,l))**2
             ELSE
                wth2(ig,l)=0.0_r8
             ENDIF
             wth3(ig,l)=zf2*(1-2.0_r8*fraca(ig,l))/(1-fraca(ig,l))  &
                  &                *zw2(ig,l)*zw2(ig,l)*zw2(ig,l)
             q2(ig,l)=zf2*(zqta(ig,l)*1000.0_r8-po(ig,l)*1000.0_r8)**2
             !test: on calcul q2/po=ratqsc
             ratqscth(ig,l)=SQRT(MAX(q2(ig,l),1.e-6_r8)/(po(ig,l)*1000.0_r8))
          ENDDO
       ENDDO
       !calcul des flux: q, thetal et thetav
       DO l=1,nlay
          DO ig=1,ngrid
             wq(ig,l)=fraca(ig,l)*zw2(ig,l)*(zqta(ig,l)*1000.0_r8-po(ig,l)*1000.0_r8)
             wthl(ig,l)=fraca(ig,l)*zw2(ig,l)*(ztla(ig,l)-zthl(ig,l))
             wthv(ig,l)=fraca(ig,l)*zw2(ig,l)*(ztva(ig,l)-ztv(ig,l))
          ENDDO
       ENDDO
       !
       IF (prt_level.GE.10) THEN
          ig=igout
          DO l=1,nlay
             PRINT*,'14f OK convect8 ig,l,zha zh zpspsk ',ig,l,zha(ig,l),zh(ig,l),zpspsk(ig,l)
             PRINT*,'14g OK convect8 ig,l,po',ig,l,po(ig,l)
          ENDDO
       ENDIF

       !      print*,'avant calcul ale et alp' 
       !calcul de ALE et ALP pour la convection
       Alp_bl(:)=0.0_r8
       Ale_bl(:)=0.0_r8
       !          print*,'ALE,ALP ,l,zw2(ig,l),Ale_bl(ig),Alp_bl(ig)'
       DO l=1,nlay
          DO ig=1,ngrid
             Alp_bl(ig)=MAX(Alp_bl(ig),0.5_r8*rhobarz(ig,l)*wth3(ig,l) )
             Ale_bl(ig)=MAX(Ale_bl(ig),0.5_r8*zw2(ig,l)**2)
             !          print*,'ALE,ALP',l,zw2(ig,l),Ale_bl(ig),Alp_bl(ig)
          ENDDO
       ENDDO

       !test:calcul de la ponderation des couches pour KE
       !initialisations

       fm_tot(:)=0.0_r8
       wght_th(:,:)=1.0_r8
       lalim_conv(:)=lalim(:)

       DO k=1,nlay
          DO ig=1,ngrid
             IF (k<=lalim_conv(ig)) fm_tot(ig)=fm_tot(ig)+fm(ig,k)
          ENDDO
       ENDDO

       ! assez bizarre car, si on est dans la couche d'alim et que alim_star et
       ! plus petit que 1.e-10, on prend wght_th=1.
       DO k=1,nlay
          DO ig=1,ngrid
             IF (k<=lalim_conv(ig).AND.alim_star(ig,k)>1.e-10_r8) THEN
                wght_th(ig,k)=alim_star(ig,k)
             ENDIF
          ENDDO
       ENDDO

       !      print*,'apres wght_th'
       !test pour prolonger la convection
       DO ig=1,ngrid
          !v1d  if ((alim_star(ig,1).lt.1.e-10).and.(therm)) then
          IF ((alim_star(ig,1).LT.1.e-10_r8)) THEN
             lalim_conv(ig)=1
             wght_th(ig,1)=1.0_r8
             !      print*,'lalim_conv ok',lalim_conv(ig),wght_th(ig,1)
          ENDIF
       ENDDO

       !------------------------------------------------------------------------
       ! Modif CR/FH 20110310 : Alp integree sur la verticale.
       ! Integrale verticale de ALP.
       ! wth3 etant aux niveaux inter-couches, on utilise d play comme masse des
       ! couches
       !------------------------------------------------------------------------

       alp_int(:)=0.0_r8
       dp_int(:)=0.0_r8
       DO l=2,nlay
          DO ig=1,ngrid
             IF(l.LE.lmax(ig)) THEN
                zdp=pplay(ig,l-1)-pplay(ig,l)
                alp_int(ig)=alp_int(ig)+0.5_r8*rhobarz(ig,l)*wth3(ig,l)*zdp
                dp_int(ig)=dp_int(ig)+zdp
             ENDIF
          ENDDO
       ENDDO

       IF (iflag_coupl>=3 .AND. iflag_coupl<=5) THEN
          DO ig=1,ngrid
             !valeur integree de alp_bl * 0.5:
             IF (dp_int(ig)>0.0_r8) THEN
                Alp_bl(ig)=alp_int(ig)/dp_int(ig)
             ENDIF
          ENDDO!
       ENDIF


       ! Facteur multiplicatif sur Alp_bl
       Alp_bl(:)=alp_bl_k*Alp_bl(:)

       !------------------------------------------------------------------------


       !calcul du ratqscdiff
       !IF (prt_level.GE.1) PRINT*,'14e OK convect8'
       var=0.0_r8
       vardiff=0.0_r8
       ratqsdiff(:,:)=0.0_r8

       DO l=1,nlay
          DO ig=1,ngrid
             IF (l<=lalim(ig)) THEN
                var=var+alim_star(ig,l)*zqta(ig,l)*1000.0_r8
             ENDIF
          ENDDO
       ENDDO

       !IF (prt_level.GE.1) PRINT*,'14f OK convect8'

       DO l=1,nlay
          DO ig=1,ngrid
             IF (l<=lalim(ig)) THEN
                zf=fraca(ig,l)
                zf2=zf/(1.0_r8-zf)
                vardiff=vardiff+alim_star(ig,l)*(zqta(ig,l)*1000.0_r8-var)**2
             ENDIF
          ENDDO
       ENDDO

       !IF (prt_level.GE.1) PRINT*,'14g OK convect8'
       DO l=1,nlay
          DO ig=1,ngrid
             ratqsdiff(ig,l)=SQRT(vardiff)/(po(ig,l)*1000.0_r8)   
             !           write(11,*)'ratqsdiff=',ratqsdiff(ig,l)
          ENDDO
       ENDDO
       !--------------------------------------------------------------------    
       !
       !ecriture des fichiers sortie
       !     print*,'15 OK convect8 CCCCCCCCCCCCCCCCCCc'

       !#ifdef wrgrads_thermcell
       !      if (prt_level.ge.1) print*,'thermcell_main sorties 3D'
       !#include "thermcell_out3d.h"
       !#endif

    ENDIF

    !IF (prt_level.GE.1) PRINT*,'thermcell_main FIN  OK'

    RETURN
  END SUBROUTINE thermcell_main


  !-----------------------------------------------------------------------------

  SUBROUTINE test_ltherm( &
       ngrid    ,&! integer      , INTENT(IN   ) :: ngrid
       nlay    ,&! integer      , INTENT(IN   ) :: nlay
       pplay   ,&! real(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
       long    ,&! integer      , INTENT(IN   ) :: long(ngrid)
       ztv     ,&! real(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
       po      ,&! real(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
       ztva    ,&! real(KIND=r8), INTENT(IN   ) :: ztva(ngrid,nlay)
       zqla    ,&! real(KIND=r8), INTENT(IN   ) :: zqla(ngrid,nlay)
       f_star  ,&! real(KIND=r8), INTENT(IN   ) :: f_star(ngrid,nlay)
       zw2     ,&! real(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay)
       comment  )! character(LEN=*), INTENT(IN) :: comment
    IMPLICIT NONE
    !#include "iniprint.h"

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
!    REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
    INTEGER      , INTENT(IN   ) :: long(ngrid)
!    REAL(KIND=r8), INTENT(IN   ) :: seuil
    REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: ztva(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zqla(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: f_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay)
    CHARACTER(LEN=*), INTENT(IN) :: comment

    INTEGER :: i
    INTEGER :: k

!    IF (prt_level.GE.1) THEN
!       PRINT*,'WARNING !!! TEST ',comment
!    ENDIF
    RETURN

    !  test sur la hauteur des thermiques ...
    DO i=1,ngrid
       !IMtemp           if (pplay(i,long(i)).lt.seuil*pplev(i,1)) then
       IF (prt_level.GE.10) THEN
          PRINT*,'WARNING ',comment,' au point ',i,' K= ',long(i)
          PRINT*,'  K  P(MB)  THV(K)     Qenv(g/kg)THVA        QLA(g/kg)   F*        W2'
          DO k=1,nlay
             WRITE(6,'(i3,7f10.3)') k,pplay(i,k),ztv(i,k),1000*po(i,k),ztva(i,k),1000*zqla(i,k),f_star(i,k),zw2(i,k)
          ENDDO
       ENDIF
    ENDDO


    RETURN
  END SUBROUTINE test_ltherm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE thermcellV1_plume(        &
       ngrid         ,&!INTEGER      , INTENT(IN   ) :: ngrid
       nlay          ,&!INTEGER      , INTENT(IN   ) :: nlay
!       ptimestep     ,&!REAL(KIND=r8), INTENT(IN   ) :: ptimestep
       ztv           ,&!REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
       zthl          ,&!REAL(KIND=r8), INTENT(IN   ) :: zthl(ngrid,nlay)
       po            ,&!REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
       zl            ,&!REAL(KIND=r8), INTENT(IN   ) :: zl(ngrid,nlay)
       rhobarz       ,&!REAL(KIND=r8), INTENT(IN   ) :: rhobarz(ngrid,nlay)
       zlev          ,&!REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
       pplev         ,&!REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
       pphi          ,&!REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
       zpspsk        ,&!REAL(KIND=r8), INTENT(IN   ) :: zpspsk(ngrid,nlay)
       alim_star     ,&!REAL(KIND=r8), INTENT(OUT  ) :: alim_star(ngrid,nlay)
       alim_star_tot ,&!REAL(KIND=r8), INTENT(OUT  ) :: alim_star_tot(ngrid)
       lalim         ,&!INTEGER      , INTENT(OUT  ) :: lalim(ngrid)
       f0            ,&!REAL(KIND=r8), INTENT(IN   ) :: f0(ngrid)
       detr_star     ,&!REAL(KIND=r8), INTENT(OUT  ) :: detr_star(ngrid,nlay)
       entr_star     ,&!REAL(KIND=r8), INTENT(OUT  ) :: entr_star(ngrid,nlay)
       f_star        ,&!REAL(KIND=r8), INTENT(OUT  ) :: f_star(ngrid,nlay+1)
       csc           ,&!REAL(KIND=r8), INTENT(OUT  ) :: csc(ngrid,nlay)
       ztva          ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztva(ngrid,nlay)
       ztla          ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztla(ngrid,nlay)
       zqla          ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqla(ngrid,nlay)
       zqta          ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqta(ngrid,nlay)
       zha           ,&!REAL(KIND=r8), INTENT(OUT  ) :: zha(ngrid,nlay)
       zw2           ,&!REAL(KIND=r8), INTENT(OUT  ) :: zw2(ngrid,nlay+1)
       w_est         ,&!REAL(KIND=r8), INTENT(OUT  ) :: w_est(ngrid,nlay+1)
       ztva_est      ,&!REAL(KIND=r8), INTENT(OUT  ) :: ztva_est(ngrid,nlay)
       zqsatth       ,&!REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(ngrid,nlay)
       lmix          ,&!INTEGER      , INTENT(OUT  ) :: lmix(ngrid)
       lmix_bis      ,&!INTEGER      , INTENT(OUT  ) :: lmix_bis(ngrid)
       linter         )!REAL(KIND=r8), INTENT(OUT  ) :: linter(ngrid)

    !--------------------------------------------------------------------------
    !thermcell_plume: calcule les valeurs de qt, thetal et w dans l ascendance
    ! Version conforme a l'article de Rio et al. 2010.
    ! Code ecrit par Catherine Rio, Arnaud Jam et Frederic Hourdin
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    !#include "YOMCST.h"
    !#include "YOETHF.h"
    !#include "FCTTRE.h"
    !#include "iniprint.h"
    !#include "thermcell.h"

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
!    REAL(KIND=r8), INTENT(IN   ) :: ptimestep
    REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zthl(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zl(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: rhobarz(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zpspsk(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: alim_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: alim_star_tot(ngrid)
    INTEGER      , INTENT(OUT  ) :: lalim(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: f0(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: detr_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: entr_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: f_star(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: csc(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ztva(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ztla(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqla(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqta(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zha(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zw2(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: w_est(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: ztva_est(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(ngrid,nlay)
    INTEGER      , INTENT(OUT  ) :: lmix(ngrid)
    INTEGER      , INTENT(OUT  ) :: lmix_bis(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: linter(ngrid)


    INTEGER :: nbpb



!    REAL(KIND=r8) :: coefc
    REAL(KIND=r8) :: detr(ngrid,nlay)
    REAL(KIND=r8) :: entr(ngrid,nlay)


    REAL(KIND=r8) :: wa_moy(ngrid,nlay+1)

    REAL(KIND=r8) :: zqla_est(ngrid,nlay)
    REAL(KIND=r8) :: zta_est(ngrid,nlay)
    REAL(KIND=r8) :: ztemp(ngrid)
    REAL(KIND=r8) :: zqsat(ngrid)
    REAL(KIND=r8) :: zdw2
!    REAL(KIND=r8) :: zw2modif
    REAL(KIND=r8) :: zw2fact
    REAL(KIND=r8) :: zeps(ngrid,nlay)

    REAL(KIND=r8) :: wmaxa(ngrid)

    INTEGER ig,l

    REAL(KIND=r8) :: zdz
    REAL(KIND=r8) :: zbuoy(ngrid,nlay)
    REAL(KIND=r8) :: zalpha
    REAL(KIND=r8) :: gamma(ngrid,nlay)
    REAL(KIND=r8) :: zdqt(ngrid,nlay)
    REAL(KIND=r8) :: zw2m
    REAL(KIND=r8) :: zbuoybis
!    REAL(KIND=r8) :: zcor
!    REAL(KIND=r8) :: zdelta
!    REAL(KIND=r8) :: zcvm5
!    REAL(KIND=r8) :: qlbef
!    REAL(KIND=r8) :: zdz2
    REAL(KIND=r8) :: betalpha
    REAL(KIND=r8) :: zbetalpha
!    REAL(KIND=r8) :: eps
    REAL(KIND=r8) :: afact
!    REAL(KIND=r8) :: REPS
    REAL(KIND=r8) :: RLvCp
!    REAL(KIND=r8), PARAMETER :: DDT0=.01_r8
    LOGICAL :: Zsat
    LOGICAL :: active(ngrid)
    LOGICAL :: activetmp(ngrid)
!    REAL(KIND=r8) :: fact_gamma
    REAL(KIND=r8) :: fact_epsilon
!    REAL(KIND=r8) :: fact_gamma2
!    REAL(KIND=r8) :: fact_epsilon2
!    REAL(KIND=r8) :: c2(ngrid,nlay)





 alim_star=0.0_r8;
 alim_star_tot=0.0_r8
 lalim=0.0_r8
 detr_star=0.0_r8;
 entr_star=0.0_r8;
 f_star=0.0_r8
 csc=0.0_r8;
 ztva=0.0_r8;
 ztla=0.0_r8;
 zqla=0.0_r8;
 zqta=0.0_r8;
 zha=0.0_r8;
 zw2=0.0_r8
 w_est=0.0_r8
 ztva_est=0.0_r8;
 zqsatth=0.0_r8;
 lmix=0.0_r8
 lmix_bis=0.0_r8
 linter=0.0_r8
 detr=0.0_r8;
 entr=0.0_r8;
 wa_moy=0.0_r8
 zqla_est=0.0_r8;
 zta_est=0.0_r8;
 ztemp=0.0_r8
 zqsat=0.0_r8
 zdw2=0.0_r8;
 zw2fact=0.0_r8;
 zeps=0.0_r8;
 wmaxa=0.0_r8
 zdz=0.0_r8;
 zbuoy=0.0_r8;
 zalpha=0.0_r8;
 gamma=0.0_r8;
 zdqt=0.0_r8;
 zw2m=0.0_r8;
 zbuoybis=0.0_r8;
 betalpha=0.0_r8;
 zbetalpha=0.0_r8;
 afact=0.0_r8;
 RLvCp=0.0_r8;
 fact_epsilon=0.0_r8;

    Zsat=.FALSE.
    ! Initialisation

    RLvCp = RLVTT/RCPD
    fact_epsilon=0.002_r8
    betalpha=0.9_r8 
    afact=2.0_r8/3.0_r8            

    zbetalpha=betalpha/(1.0_r8+betalpha)


    ! Initialisations des variables reeles
    IF (1==0) THEN
       ztva(:,:)=ztv(:,:)
       ztva_est(:,:)=ztva(:,:)
       ztla(:,:)=zthl(:,:)
       zqta(:,:)=po(:,:)
       zha(:,:) = ztva(:,:)
    ELSE
       ztva(:,:)=0.0_r8
       ztva_est(:,:)=0.0_r8
       ztla(:,:)=0.0_r8
       zqta(:,:)=0.0_r8
       zha(:,:) =0.0_r8
    ENDIF

    zqla_est(:,:)=0.0_r8
    zqsatth(:,:)=0.0_r8
    zqla(:,:)=0.0_r8
    detr_star(:,:)=0.0_r8
    entr_star(:,:)=0.0_r8
    alim_star(:,:)=0.0_r8
    alim_star_tot(:)=0.0_r8
    csc(:,:)=0.0_r8
    detr(:,:)=0.0_r8
    entr(:,:)=0.0_r8
    zw2(:,:)=0.0_r8
    zbuoy(:,:)=0.0_r8
    gamma(:,:)=0.0_r8
    zeps(:,:)=0.0_r8
    w_est(:,:)=0.0_r8
    f_star(:,:)=0.0_r8
    wa_moy(:,:)=0.0_r8
    linter(:)=1.0_r8
    !     linter(:)=1.0_r8
    ! Initialisation des variables entieres
    lmix(:)=1
    lmix_bis(:)=2
    wmaxa(:)=0.0_r8
    lalim(:)=1


    !-------------------------------------------------------------------------
    ! On ne considere comme actif que les colonnes dont les deux premieres
    ! couches sont instables.
    !-------------------------------------------------------------------------
    active(:)=ztv(:,1)>ztv(:,2)

    !-------------------------------------------------------------------------
    ! Definition de l'alimentation a l'origine dans thermcell_init
    !-------------------------------------------------------------------------
    DO l=1,nlay-1
       DO ig=1,ngrid
          IF (ztv(ig,l)> ztv(ig,l+1) .AND. ztv(ig,1)>=ztv(ig,l) ) THEN
             alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.0_r8)  &
                  &                       *SQRT(zlev(ig,l+1)) 
             lalim(ig)=l+1
             alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          ENDIF
       ENDDO
    ENDDO
    DO l=1,nlay
       DO ig=1,ngrid 
          IF (alim_star_tot(ig) > 1.e-10_r8 ) THEN
             alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
          ENDIF
       ENDDO
    ENDDO
    alim_star_tot(:)=1.0_r8



    !------------------------------------------------------------------------------
    ! Calcul dans la premiere couche
    ! On decide dans cette version que le thermique n'est actif que si la premiere
    ! couche est instable.
    ! Pourrait etre change si on veut que le thermiques puisse se dÃ©clencher
    ! dans une couche l>1
    !------------------------------------------------------------------------------
    DO ig=1,ngrid
       ! Le panache va prendre au debut les caracteristiques de l'air contenu
       ! dans cette couche.
       IF (active(ig)) THEN
          ztla(ig,1)=zthl(ig,1) 
          zqta(ig,1)=po(ig,1)
          zqla(ig,1)=zl(ig,1)
          !cr: attention, prise en compte de f*(1)=1
          f_star(ig,2)=alim_star(ig,1)
          zw2(ig,2)=2.0_r8*RG*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
               &                     *(zlev(ig,2)-zlev(ig,1))  &
               &                     *0.4_r8*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1))
          w_est(ig,2)=zw2(ig,2)
       ENDIF
    ENDDO
    !

    !==============================================================================
    !boucle de calcul de la vitesse verticale dans le thermique
    !==============================================================================
    DO l=2,nlay-1
       !==============================================================================


       ! On decide si le thermique est encore actif ou non
       ! AFaire : Il faut sans doute ajouter entr_star a alim_star dans ce test
       DO ig=1,ngrid
          active(ig)=active(ig) &
               &                 .AND. zw2(ig,l)>1.e-10_r8 &
               &                 .AND. f_star(ig,l)+alim_star(ig,l)>1.e-10_r8
       ENDDO



       !---------------------------------------------------------------------------
       ! calcul des proprietes thermodynamiques et de la vitesse de la couche l
       ! sans tenir compte du detrainement et de l'entrainement dans cette
       ! couche
       ! C'est a dire qu'on suppose 
       ! ztla(l)=ztla(l-1) et zqta(l)=zqta(l-1)
       ! Ici encore, on doit pouvoir ajouter entr_star (qui peut etre calculer
       ! avant) a l'alimentation pour avoir un calcul plus propre
       !---------------------------------------------------------------------------

       DO ig=1,ngrid 
           ztemp(ig)=zpspsk(ig,l)*ztla(ig,l-1)
       END DO
       CALL thermcell_qsat(&
            ngrid                ,&! INTEGER      , INTENT(IN   ) :: ngrid
            active(1:ngrid)      ,&! LOGICAL      , INTENT(IN   ) :: active(ngrid)
            pplev (1:ngrid,l)    ,&! REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid)
            ztemp (1:ngrid)      ,&! REAL(KIND=r8), INTENT(IN   ) :: ztemp(ngrid)
            zqta  (1:ngrid,l-1)  ,&! REAL(KIND=r8), INTENT(IN   ) :: zqta(ngrid)
            zqsat (1:ngrid)       )! REAL(KIND=r8), INTENT(OUT  ) :: zqsat(ngrid)

       DO ig=1,ngrid 
          !       print*,'active',active(ig),ig,l
          IF(active(ig)) THEN 
             zqla_est(ig,l)=MAX(0.0_r8,zqta(ig,l-1)-zqsat(ig))
             ztva_est(ig,l) = ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*zqla_est(ig,l)
             zta_est(ig,l)=ztva_est(ig,l)
             ztva_est(ig,l) = ztva_est(ig,l)/zpspsk(ig,l)
             ztva_est(ig,l) = ztva_est(ig,l)*(1.0_r8+RETV*(zqta(ig,l-1)  &
                  &      -zqla_est(ig,l))-zqla_est(ig,l))

             !------------------------------------------------
             !AJAM:nouveau calcul de w²  
             !------------------------------------------------
             zdz=zlev(ig,l+1)-zlev(ig,l)
             zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)

             zw2fact=fact_epsilon*2.0_r8*zdz/(1.0_r8+betalpha)
             zdw2=(afact)*zbuoy(ig,l)/(fact_epsilon)
             w_est(ig,l+1)=MAX(0.0001_r8,EXP(-zw2fact)*(w_est(ig,l)-zdw2)+zdw2)


             IF (w_est(ig,l+1).LT.0.0_r8) THEN
                w_est(ig,l+1)=zw2(ig,l)
             ENDIF
          ENDIF
       ENDDO


       !-------------------------------------------------
       !calcul des taux d'entrainement et de detrainement
       !-------------------------------------------------

       DO ig=1,ngrid
          IF (active(ig)) THEN

             zw2m=MAX(0.5_r8*(w_est(ig,l)+w_est(ig,l+1)),0.1_r8)
             zw2m=w_est(ig,l+1)
             zdz=zlev(ig,l+1)-zlev(ig,l)
             zbuoy(ig,l)=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
             !          zbuoybis=zbuoy(ig,l)+RG*0.1/300.
             zbuoybis=zbuoy(ig,l)
             zalpha=f0(ig)*f_star(ig,l)/SQRT(w_est(ig,l+1))/rhobarz(ig,l)
             zdqt(ig,l)=MAX(zqta(ig,l-1)-po(ig,l),0.0_r8)/po(ig,l)


             entr_star(ig,l)=f_star(ig,l) * zdz *  zbetalpha*MAX(0.0_r8,afact*zbuoybis/zw2m - fact_epsilon )


             detr_star(ig,l)=f_star(ig,l)*zdz                        &
                  &     *MAX(1.e-3_r8, -afact*zbetalpha*zbuoy(ig,l)/zw2m          &
                  &     + 0.012_r8*(zdqt(ig,l)/zw2m)**0.5_r8 )

             ! En dessous de lalim, on prend le max de alim_star et entr_star pour
             ! alim_star et 0 sinon
             IF (l.LT.lalim(ig)) THEN
                alim_star(ig,l)=MAX(alim_star(ig,l),entr_star(ig,l))
                entr_star(ig,l)=0.0_r8
             ENDIF

             ! Calcul du flux montant normalise
             f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
                  &              -detr_star(ig,l)

          ENDIF
       ENDDO


       !----------------------------------------------------------------------------
       !calcul de la vitesse verticale en melangeant Tl et qt du thermique
       !---------------------------------------------------------------------------
       DO ig=1,ngrid
          activetmp(ig)=active(ig) .AND. f_star(ig,l+1)>1.e-10_r8
       END DO
       DO ig=1,ngrid
          IF (activetmp(ig)) THEN 
             Zsat=.FALSE.
             ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
                  &            (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))  &
                  &            /(f_star(ig,l+1)+detr_star(ig,l))
             zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+  &
                  &            (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))  &
                  &            /(f_star(ig,l+1)+detr_star(ig,l))

          ENDIF
       ENDDO
       DO ig=1,ngrid
           ztemp(ig)=zpspsk(ig,l)*ztla(ig,l)
       END DO
       CALL thermcell_qsat( &
            ngrid                  ,&! INTEGER      , INTENT(IN      ) :: ngrid
            activetmp(1:ngrid)     ,&! LOGICAL      , INTENT(IN      ) :: active(ngrid)
            pplev    (1:ngrid,l)   ,&! REAL(KIND=r8), INTENT(IN      ) :: pplev(ngrid)
            ztemp    (1:ngrid)     ,&! REAL(KIND=r8), INTENT(IN      ) :: ztemp(ngrid)
            zqta     (1:ngrid,l)   ,&! REAL(KIND=r8), INTENT(IN      ) :: zqta(ngrid)
            zqsatth  (1:ngrid,l)    )! REAL(KIND=r8), INTENT(OUT     ) :: zqsat(ngrid)

       DO ig=1,ngrid
          IF (activetmp(ig)) THEN
             ! on ecrit de maniere conservative (sat ou non)
             !          T = Tl +Lv/Cp ql
             zqla(ig,l)=MAX(0.0_r8,zqta(ig,l)-zqsatth(ig,l))
             ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
             ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
             !on rajoute le calcul de zha pour diagnostiques (temp potentielle)
             zha(ig,l) = ztva(ig,l)
             ztva(ig,l) = ztva(ig,l)*(1.0_r8+RETV*(zqta(ig,l)  &
                  &              -zqla(ig,l))-zqla(ig,l))
             zbuoy(ig,l)=RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
             zdz=zlev(ig,l+1)-zlev(ig,l)
             zeps(ig,l)=(entr_star(ig,l)+alim_star(ig,l))/(f_star(ig,l)*zdz)

             zw2fact=fact_epsilon*2.0_r8*zdz/(1.0_r8+betalpha)
             zdw2=afact*zbuoy(ig,l)/(fact_epsilon)
             zw2(ig,l+1)=MAX(0.0001_r8,EXP(-zw2fact)*(zw2(ig,l)-zdw2)+zdw2) 
          ENDIF
       ENDDO

       IF (prt_level.GE.20) PRINT*,'coucou calcul detr 460: ig, l',ig, l
       !
       !---------------------------------------------------------------------------
       !initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max 
       !---------------------------------------------------------------------------

       nbpb=0
       DO ig=1,ngrid
          IF (zw2(ig,l+1)>0.0_r8 .AND. zw2(ig,l+1).LT.1.e-10_r8) THEN
             !               stop'On tombe sur le cas particulier de thermcell_dry'
             !               print*,'On tombe sur le cas particulier de thermcell_plume'
             nbpb=nbpb+1
             zw2(ig,l+1)=0.0_r8
             linter(ig)=l+1
          ENDIF

          IF (zw2(ig,l+1).LT.0.0_r8) THEN 
             linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
                  &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
             zw2(ig,l+1)=0.0_r8
          ENDIF

          wa_moy(ig,l+1)=SQRT(zw2(ig,l+1)) 

          IF (wa_moy(ig,l+1).GT.wmaxa(ig)) THEN
             !   lmix est le niveau de la couche ou w (wa_moy) est maximum
             !on rajoute le calcul de lmix_bis
             IF (zqla(ig,l).LT.1.e-10_r8) THEN
                lmix_bis(ig)=l+1
             ENDIF
             lmix(ig)=l+1
             wmaxa(ig)=wa_moy(ig,l+1)
          ENDIF
       ENDDO

       IF (nbpb>0) THEN
          PRINT*,'WARNING on tombe ',nbpb,' x sur un pb pour l=',l,' dans thermcell_plume'
       ENDIF

       !=========================================================================
       ! FIN DE LA BOUCLE VERTICALE
    ENDDO
    !=========================================================================

    !on recalcule alim_star_tot
    DO ig=1,ngrid
       alim_star_tot(ig)=0.0_r8
    ENDDO
    DO ig=1,ngrid
       DO l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
       ENDDO
    ENDDO


    RETURN 
  END  SUBROUTINE thermcellV1_plume



  SUBROUTINE thermcell_height(&
       ngrid   ,&!INTEGER      , INTENT(IN   ) :: ngrid
       nlay    ,&!INTEGER      , INTENT(IN   ) :: nlay
       lalim   ,&!INTEGER      , INTENT(INOUT) :: lalim(ngrid)
       lmin    ,&!INTEGER      , INTENT(INOUT) :: lmin(ngrid)
       linter  ,&!REAL(KIND=r8), INTENT(IN   ) :: linter(ngrid)
       lmix    ,&!INTEGER      , INTENT(INOUT) :: lmix(ngrid)
       zw2     ,&!REAL(KIND=r8), INTENT(INOUT) :: zw2(ngrid,nlay+1)
       zlev    ,&!REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
       lmax    ,&!INTEGER      , INTENT(OUT  ) :: lmax(ngrid)
       zmax    ,&!REAL(KIND=r8), INTENT(OUT  ) :: zmax(ngrid)
       zmax0   ,&!REAL(KIND=r8), INTENT(INOUT) :: zmax0(ngrid)
       zmix    ,&!REAL(KIND=r8), INTENT(OUT  ) :: zmix(ngrid)
       wmax    )!REAL(KIND=r8), INTENT(OUT  ) :: wmax(ngrid)
!       lev_out  )!INTEGER      , INTENT(IN   ) :: lev_out

    !-----------------------------------------------------------------------------
    !thermcell_height: calcul des caracteristiques du thermique: zmax,wmax,zmix
    !-----------------------------------------------------------------------------
    IMPLICIT NONE
    !#include "iniprint.h"
    !#include "thermcell.h"

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    INTEGER      , INTENT(INOUT) :: lalim(ngrid)
    INTEGER      , INTENT(INOUT) :: lmin(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: linter(ngrid)
    INTEGER      , INTENT(INOUT) :: lmix(ngrid)
    REAL(KIND=r8), INTENT(INOUT) :: zw2(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
    INTEGER      , INTENT(OUT  ) :: lmax(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: zmax(ngrid)
    REAL(KIND=r8), INTENT(INOUT) :: zmax0(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: zmix(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: wmax(ngrid)
!    INTEGER      , INTENT(IN   ) :: lev_out                           ! niveau pour les print

    REAL(KIND=r8) :: num(ngrid)
    REAL(KIND=r8) :: denom(ngrid)
    INTEGER       :: ig
    INTEGER       :: l

    REAL(KIND=r8) :: zlevinter(ngrid)




    lmax=0.0_r8;
    zmax=0.0_r8;
    zmix=0.0_r8;
    wmax=0.0_r8;
    num=0.0_r8;
    denom=0.0_r8;
    zlevinter=0.0_r8;
    !calcul de la hauteur max du thermique
    DO ig=1,ngrid
       lmax(ig)=lalim(ig)
    ENDDO

    DO ig=1,ngrid
       DO l=nlay,lalim(ig)+1,-1
          IF (zw2(ig,l).LE.1.e-10_r8) THEN
             lmax(ig)=l-1
          ENDIF
       ENDDO
    ENDDO

    ! On traite le cas particulier qu'il faudrait éviter ou le thermique
    ! atteind le haut du modele ...
    DO ig=1,ngrid
       IF ( zw2(ig,nlay) > 1.e-10_r8 ) THEN
          PRINT*,'WARNING !!!!! W2 thermiques non nul derniere couche '
          lmax(ig)=nlay
       ENDIF
    ENDDO

    ! pas de thermique si couche 1 stable
    DO ig=1,ngrid
       IF (lmin(ig).GT.1) THEN
          lmax(ig)=1
          lmin(ig)=1
          lalim(ig)=1
       ENDIF
    ENDDO
    !    
    ! Determination de zw2 max
    DO ig=1,ngrid
       wmax(ig)=0.0_r8
    ENDDO

    DO l=1,nlay
       DO ig=1,ngrid
          IF (l.LE.lmax(ig)) THEN
             IF (zw2(ig,l).LT.0.0_r8)THEN
                zw2(ig,l)=0.0_r8
                !PRINT*,'pb2 zw2<0'
             ENDIF
             zw2(ig,l)=SQRT(zw2(ig,l))
             wmax(ig)=MAX(wmax(ig),zw2(ig,l))
          ELSE
             zw2(ig,l)=0.0_r8
          ENDIF
       ENDDO
    ENDDO

    !   Longueur caracteristique correspondant a la hauteur des thermiques.
    DO  ig=1,ngrid
       zmax(ig)=0.0_r8
       zlevinter(ig)=zlev(ig,1)
    ENDDO

    IF (iflag_thermals_ed.GE.1) THEN

       num(:)=0.0_r8
       denom(:)=0.0_r8
       DO ig=1,ngrid
          DO l=1,nlay
             num(ig)=num(ig)+zw2(ig,l)*zlev(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
             denom(ig)=denom(ig)+zw2(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
          ENDDO
       ENDDO
       DO ig=1,ngrid
          IF (denom(ig).GT.1.e-10_r8) THEN
             zmax(ig)=2.0_r8*num(ig)/denom(ig)
             zmax0(ig)=zmax(ig)
          ENDIF
       ENDDO

    ELSE

       DO  ig=1,ngrid
          ! calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*  &
               &    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)  &
               &    -zlev(ig,lmax(ig)))
          !pour le cas ou on prend tjs lmin=1
          !       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
          zmax(ig)=MAX(zmax(ig),zlevinter(ig)-zlev(ig,1))
          zmax0(ig)=zmax(ig)
       ENDDO


    ENDIF
    !endif iflag_thermals_ed
    !
    ! def de  zmix continu (profil parabolique des vitesses)
    DO ig=1,ngrid
       IF (lmix(ig).GT.1) THEN
          ! test 
          IF (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))  &
               &        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))  &
               &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))  &
               &        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))).GT.1e-10_r8)  &
               &        THEN
             !             
             zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))  &
                  &        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)  &
                  &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))  &
                  &        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))  &
                  &        /(2.0_r8*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))  &
                  &        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))  &
                  &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))  &
                  &        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
          ELSE
             zmix(ig)=zlev(ig,lmix(ig))
!             PRINT*,'pb zmix'
          ENDIF
       ELSE 
          zmix(ig)=0.0_r8
       ENDIF
       !test
       IF ((zmax(ig)-zmix(ig)).LE.0.0_r8) THEN
          zmix(ig)=0.9_r8*zmax(ig)
          !            print*,'pb zmix>zmax'
       ENDIF
    ENDDO
    !
    ! calcul du nouveau lmix correspondant
    DO ig=1,ngrid
       DO l=1,nlay
          IF (zmix(ig).GE.zlev(ig,l).AND.  &
               &          zmix(ig).LT.zlev(ig,l+1)) THEN
             lmix(ig)=l
          ENDIF
       ENDDO
    ENDDO
    !
    RETURN 
  END SUBROUTINE thermcell_height

  SUBROUTINE thermcell_dq( &
       ngrid    , &! integer      , INTENT(IN   ) :: ngrid
       nlay     , &! integer      , INTENT(IN   ) :: nlay
       ptimestep, &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
       fm       , &! REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
       entr     , &! REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
       masse    , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
       q        , &! REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
       dq       , &! REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
       imask    , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
       dump     , &! 
       qa         )! REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)


    IMPLICIT NONE

    !#include "iniprint.h"
    !=======================================================================
    !
    !   Calcul du transport verticale dans la couche limite en presence
    !   de "thermiques" explicitement representes
    !   calcul du dq/dt une fois qu'on connait les ascendances
    !
    !=======================================================================

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    REAL(KIND=r8), INTENT(IN   ) :: ptimestep
    REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
    REAL(KIND=r8), INTENT(INOUT) :: q(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: dq(ngrid,nlay)
    INTEGER(KIND=i8), INTENT(IN) :: imask (ngrid)
    REAL(KIND=r8), INTENT(INOUT) :: dump(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: qa(ngrid,nlay)
!    INTEGER      , INTENT(IN   ) :: lev_out                           ! niveau pour les print

    REAL(KIND=r8)    :: detr(ngrid,nlay)
    REAL(KIND=r8)    :: wqd(ngrid,nlay+1)

    REAL(KIND=r8)    :: zzm

    INTEGER          :: ig,k
    REAL(KIND=r8)    :: cfl

    REAL(KIND=r8)    :: qold(ngrid,nlay)
    REAL(KIND=r8)    :: ztimestep
    INTEGER          :: niter
    INTEGER          :: iter
!    CHARACTER (LEN=20) :: modname='thermcell_dq'
    CHARACTER (LEN=80) :: abort_message

!    RETURN
    
    dq= 0.0_r8
    qa= 0.0_r8
    detr= 0.0_r8
    wqd= 0.0_r8
    cfl= 0.0_r8
    zzm= 0.0_r8
    qold= 0.0_r8
    ztimestep= 0.0_r8
    ! Calcul du critere CFL pour l'advection dans la subsidence
    cfl = 0.0_r8
    DO k=1,nlay
       DO ig=1,ngrid
          zzm=masse(ig,k)/ptimestep
          cfl=MAX(cfl,fm(ig,k)/zzm)
          IF (entr(ig,k).GT.zzm) THEN
             PRINT*,'entr dt > m ',entr(ig,k)*ptimestep,masse(ig,k)
             entr(ig,k) = zzm
             !RETURN
             !               CALL abort_gcm (modname,abort_message,1)
          ENDIF
       ENDDO
    ENDDO

!   print*,'CFL CFL CFL CFL ',cfl

    IF(CFL_TEST)THEN
       ! On subdivise le calcul en niter pas de temps.
       niter=3          !INT(cfl)+1  not function 
    ELSE
       niter=1
    END IF

    ztimestep=ptimestep/niter
    qold=q


    DO iter=1,niter
!       IF (prt_level.GE.1) PRINT*,'Q2 THERMCEL_DQ 0'

       !   calcul du detrainement
       DO k=1,nlay
          DO ig=1,ngrid
             detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
             !           print*,'Q2 DQ ',detr(ig,k),fm(ig,k),entr(ig,k)
             !test
             IF (detr(ig,k).LT.0.0_r8) THEN
                entr(ig,k)=entr(ig,k)-detr(ig,k)
                detr(ig,k)=0.0_r8
                !               print*,'detr2<0!!!','ig=',ig,'k=',k,'f=',fm(ig,k),
                !     s         'f+1=',fm(ig,k+1),'e=',entr(ig,k),'d=',detr(ig,k)
             ENDIF
             IF (fm(ig,k+1).LT.0.0_r8) THEN
                !               print*,'fm2<0!!!'
             ENDIF
             IF (entr(ig,k).LT.0.0_r8) THEN
                !               print*,'entr2<0!!!'
             ENDIF
          ENDDO
       ENDDO

       !   calcul de la valeur dans les ascendances
       DO ig=1,ngrid
          qa(ig,1)=q(ig,1)
       ENDDO

       DO k=2,nlay
          DO ig=1,ngrid
             IF ((fm(ig,k+1)+detr(ig,k))*ztimestep.GT.  &
                  &         1.e-5_r8*masse(ig,k)) THEN
                qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))  &
                     &         /(fm(ig,k+1)+detr(ig,k))
             ELSE
                qa(ig,k)=q(ig,k)
             ENDIF
             IF (qa(ig,k).LT.0.0_r8) THEN
                !               print*,'qa<0!!!'
             ENDIF
             IF (q(ig,k).LT.0.0_r8) THEN
                !               print*,'q<0!!!'
             ENDIF
          ENDDO
       ENDDO

       ! Calcul du flux subsident

       DO k=2,nlay
          DO ig=1,ngrid
             IF(centre)THEN
                wqd(ig,k)=fm(ig,k)*0.5_r8*(q(ig,k-1)+q(ig,k))
             ELSE
                IF(plusqueun)THEN
                   ! Schema avec advection sur plus qu'une maille.
                   zzm=masse(ig,k)/ztimestep
                   IF (fm(ig,k)>zzm) THEN
                      wqd(ig,k)=zzm*q(ig,k)+(fm(ig,k)-zzm)*q(ig,k+1)
                   ELSE
                      wqd(ig,k)=fm(ig,k)*q(ig,k)
                   ENDIF
                ELSE
                   wqd(ig,k)=fm(ig,k)*q(ig,k)
                END IF
             END IF
             IF (wqd(ig,k).LT.0.0_r8) THEN
                !               print*,'wqd<0!!!'
             ENDIF
          ENDDO
       ENDDO
       DO ig=1,ngrid
          wqd(ig,1)=0.0_r8
          wqd(ig,nlay+1)=0.0_r8
       ENDDO


       ! Calcul des tendances
       DO k=1,nlay
          DO ig=1,ngrid
             q(ig,k)=q(ig,k)+(detr(ig,k)*qa(ig,k)-entr(ig,k)*q(ig,k)  &
                  &               -wqd(ig,k)+wqd(ig,k+1))  &
                  &               *ztimestep/masse(ig,k)
             !            if (dq(ig,k).lt.0.0_r8) then
             !               print*,'dq<0!!!'
             !            endif
          ENDDO
       ENDDO


    ENDDO


    ! Calcul des tendances
    DO k=1,nlay
       DO ig=1,ngrid
          IF(imask(ig) >=1_i8)THEN
             dq(ig,k)=0.1_r8*(q(ig,k)-qold(ig,k))/ptimestep
             q(ig,k)=qold(ig,k)
          ELSE 
             dq(ig,k)=(q(ig,k)-qold(ig,k))/ptimestep
             q(ig,k)=qold(ig,k)
          END IF 
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE thermcell_dq

  SUBROUTINE thermcell_dv2( &
       ngrid     ,&!integer      , INTENT(IN   ) :: ngrid
       nlay      ,&!integer      , INTENT(IN   ) :: nlay
       ptimestep ,&!real(KIND=r8), INTENT(IN   ) :: ptimestep
       zlay      ,&!(ngrid,nlay)
       pblh       ,&
       imask      , &!   INTEGER(KIND=i8), INTENT(INOUT) :: imask (ncols)
       fm        ,&!real(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
       entr      ,&!real(KIND=r8), INTENT(IN   ) :: entr(ngrid,nlay)
       masse     ,&!real(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
       fraca     ,&!real(KIND=r8), INTENT(IN   ) :: fraca(ngrid,nlay+1)
       larga     ,&!real(KIND=r8), INTENT(IN   ) :: larga(ngrid)
       u         ,&!real(KIND=r8), INTENT(IN   ) :: u(ngrid,nlay)
       v         ,&!real(KIND=r8), INTENT(IN   ) :: v(ngrid,nlay)
       du        ,&!real(KIND=r8), INTENT(OUT  ) :: du(ngrid,nlay)
       dv        ,&!real(KIND=r8), INTENT(OUT  ) :: dv(ngrid,nlay)
       ua        ,&!real(KIND=r8), INTENT(OUT  ) :: ua(ngrid,nlay)
       va        )!real(KIND=r8), INTENT(OUT  ) :: va(ngrid,nlay)
!       lev_out    )!integer      , INTENT(IN   ) :: lev_out      ! niveau pour les print
    IMPLICIT NONE

    !#include "iniprint.h"
    !=======================================================================
    !
    !   Calcul du transport verticale dans la couche limite en presence
    !   de "thermiques" explicitement representes
    !   calcul du dq/dt une fois qu'on connait les ascendances
    !
    ! Vectorisation, FH : 2010/03/08
    !
    !=======================================================================


    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    REAL(KIND=r8), INTENT(IN   ) :: ptimestep
    REAL(KIND=r8), INTENT(IN   ) :: zlay(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pblh (ngrid)
    INTEGER(KIND=i8), INTENT(IN) :: imask (ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: fm(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(INOUT) :: entr(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: fraca(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: larga(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: u(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: v(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: du(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: dv(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ua(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: va(ngrid,nlay)
!    INTEGER      , INTENT(IN   ) :: lev_out                           ! niveau pour les print




!    REAL(KIND=r8)    :: qa(ngrid,nlay)
    REAL(KIND=r8)    :: detr(ngrid,nlay)
    REAL(KIND=r8)    :: zf
    REAL(KIND=r8)    :: zf2
    REAL(KIND=r8)    :: wvd(ngrid,nlay+1)
    REAL(KIND=r8)    :: wud(ngrid,nlay+1)
    REAL(KIND=r8)    :: gamma0(ngrid,nlay+1)
    REAL(KIND=r8)    :: gamma(ngrid,nlay+1)
    REAL(KIND=r8)    :: ue(ngrid,nlay)
    REAL(KIND=r8)    :: ve(ngrid,nlay)
    LOGICAL          :: ltherm(ngrid,nlay)
    REAL(KIND=r8)    :: dua(ngrid,nlay)
    REAL(KIND=r8)    :: dva(ngrid,nlay)
    INTEGER          :: iter
    REAL(KIND=r8)    :: cfl,zzm
    INTEGER :: ig
    INTEGER :: k,niter
    INTEGER :: nlarga0
    INTEGER :: lunout=0

du=0.0_r8
dv=0.0_r8
ua=0.0_r8
va=0.0_r8
detr=0.0_r8
zf=0.0_r8
zf2=0.0_r8
wvd=0.0_r8
wud=0.0_r8
gamma0=0.0_r8
gamma=0.0_r8
ue=0.0_r8
ve=0.0_r8
dua=0.0_r8
dva=0.0_r8


!      RETURN!PKUBOTA4
    !-------------------------------------------------------------------------
    ! Calcul du critere CFL pour l'advection dans la subsidence
    cfl = 0.0_r8
    DO k=1,nlay
       DO ig=1,ngrid
          zzm=masse(ig,k)/ptimestep
          cfl=MAX(cfl,fm(ig,k)/zzm)
          IF (entr(ig,k).GT.zzm) THEN
             PRINT*,'entr dt > m ',entr(ig,k)*ptimestep,masse(ig,k)
             entr(ig,k) = zzm
             !RETURN!STOP
             !               CALL abort_gcm (modname,abort_message,1)
          ENDIF
       ENDDO
    ENDDO

    !IM 090508     print*,'CFL CFL CFL CFL ',cfl

    IF(CFL_TEST)THEN
       ! On subdivise le calcul en niter pas de temps.
       niter=INT(cfl)+1
    ELSE
       niter=1
    END IF

    !   calcul du detrainement
    !---------------------------

    !      print*,'THERMCELL DV2 OPTIMISE 3'

    nlarga0=0

    DO k=1,nlay
       DO ig=1,ngrid
          detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
       ENDDO
    ENDDO

    !   calcul de la valeur dans les ascendances
    DO ig=1,ngrid
       ua(ig,1)=u(ig,1)
       va(ig,1)=v(ig,1)
       ue(ig,1)=u(ig,1)
       ve(ig,1)=v(ig,1)
    ENDDO

    IF(prt_level>9)WRITE(lunout,*)                                    &
         &      'WARNING on initialise gamma(1:ngrid,1)=0.0_r8'
    gamma(1:ngrid,1)=0.0_r8
    DO k=2,nlay
       DO ig=1,ngrid
          ltherm(ig,k)=(fm(ig,k+1)+detr(ig,k))*ptimestep > 1.e-5_r8*masse(ig,k)
          IF(ltherm(ig,k).AND.larga(ig)>0.0_r8) THEN
             gamma0(ig,k)=masse(ig,k)  &
                  &         *SQRT( 0.5_r8*(fraca(ig,k+1)+fraca(ig,k)) )  &
                  &         *0.5_r8/larga(ig)  &
                  &         *1.0_r8
          ELSE
             gamma0(ig,k)=0.0_r8
          ENDIF
          IF (ltherm(ig,k).AND.larga(ig)<=0.0_r8) nlarga0=nlarga0+1
       ENDDO
    ENDDO

    gamma(:,:)=0.0_r8

    DO k=2,nlay

       DO ig=1,ngrid
          IF (ltherm(ig,k)) THEN
             dua(ig,k)=ua(ig,k-1)-u(ig,k-1)
             dva(ig,k)=va(ig,k-1)-v(ig,k-1)
          ELSE
             ua(ig,k)=u(ig,k)
             va(ig,k)=v(ig,k)
             ue(ig,k)=u(ig,k)
             ve(ig,k)=v(ig,k)
          ENDIF
       ENDDO


       ! Debut des iterations
       !----------------------
       DO iter=1,5
          DO ig=1,ngrid
             ! Pour memoire : calcul prenant en compte la fraction reelle
             zf=0.5_r8*(fraca(ig,k)+fraca(ig,k+1))
             zf2=1.0_r8/(1.0_r8-zf)
             ! Calcul avec fraction infiniement petite
             !zf=0.0_r8
             !zf2=1.0_r8

             !  la première fois on multiplie le coefficient de freinage
             !  par le module du vent dans la couche en dessous.
             !  Mais pourquoi donc ???
             IF (ltherm(ig,k)) THEN
                !   On choisit une relaxation lineaire.
                !                 gamma(ig,k)=gamma0(ig,k)
                !   On choisit une relaxation quadratique.
                gamma(ig,k)=gamma0(ig,k)*SQRT(dua(ig,k)**2+dva(ig,k)**2)
                ua(ig,k)=(fm(ig,k)*ua(ig,k-1)  &
                     &               +(zf2*entr(ig,k)+gamma(ig,k))*u(ig,k))  &
                     &               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2  &
                     &                 +gamma(ig,k))
                va(ig,k)=(fm(ig,k)*va(ig,k-1)  &
                     &               +(zf2*entr(ig,k)+gamma(ig,k))*v(ig,k))  &
                     &               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2  &
                     &                 +gamma(ig,k))
                !                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua(ig,k),dva(ig,k)
                dua(ig,k)=ua(ig,k)-u(ig,k)
                dva(ig,k)=va(ig,k)-v(ig,k)
                ue(ig,k)=(u(ig,k)-zf*ua(ig,k))*zf2
                ve(ig,k)=(v(ig,k)-zf*va(ig,k))*zf2
             ENDIF
          ENDDO
          ! Fin des iterations
          !--------------------
       ENDDO

    ENDDO ! k=2,nlay


    ! Calcul du flux vertical de moment dans l'environnement.
    !---------------------------------------------------------
    DO k=2,nlay
       DO ig=1,ngrid
          wud(ig,k)=fm(ig,k)*ue(ig,k)
          wvd(ig,k)=fm(ig,k)*ve(ig,k)
       ENDDO
    ENDDO
    DO ig=1,ngrid
       wud(ig,1)=0.0_r8
       wud(ig,nlay+1)=0.0_r8
       wvd(ig,1)=0.0_r8
       wvd(ig,nlay+1)=0.0_r8
    ENDDO

    ! calcul des tendances.
    !-----------------------
    DO k=1,nlay
       DO ig=1,ngrid
          IF(imask(ig) >=1_i8)THEN
!            IF(zlay(ig,k) <= 1000.0_r8)THEN
             IF(zlay(ig,k) <= MAX(pblh(ig),800.0_r8))THEN
                du(ig,k)=-0.10_r8*((detr(ig,k)+gamma(ig,k))*ua(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)  &
                  &               -wud(ig,k)+wud(ig,k+1))  &
                  &               /masse(ig,k)
                dv(ig,k)=-0.10_r8*((detr(ig,k)+gamma(ig,k))*va(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)  &
                  &               -wvd(ig,k)+wvd(ig,k+1))  &
                  &               /masse(ig,k)
             ELSE
                du(ig,k)=-0.10_r8*((detr(ig,k)+gamma(ig,k))*ua(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)  &
                  &               -wud(ig,k)+wud(ig,k+1))  &
                  &               /masse(ig,k)
                dv(ig,k)=-0.10_r8*((detr(ig,k)+gamma(ig,k))*va(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)  &
                  &               -wvd(ig,k)+wvd(ig,k+1))  &
                  &               /masse(ig,k)
             END IF
          ELSE
!            IF(zlay(ig,k) <= 1000.0_r8)THEN
             IF(zlay(ig,k) <= MAX(pblh(ig),800.0_r8))THEN
                du(ig,k)=-0.50_r8*((detr(ig,k)+gamma(ig,k))*ua(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)  &
                  &               -wud(ig,k)+wud(ig,k+1))  &
                  &               /masse(ig,k)
                dv(ig,k)=-0.50_r8*((detr(ig,k)+gamma(ig,k))*va(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)  &
                  &               -wvd(ig,k)+wvd(ig,k+1))  &
                  &               /masse(ig,k)
             ELSE
                du(ig,k)=-0.250_r8*((detr(ig,k)+gamma(ig,k))*ua(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)  &
                  &               -wud(ig,k)+wud(ig,k+1))  &
                  &               /masse(ig,k)
                dv(ig,k)=-0.250_r8*((detr(ig,k)+gamma(ig,k))*va(ig,k)  &
                  &               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)  &
                  &               -wvd(ig,k)+wvd(ig,k+1))  &
                  &               /masse(ig,k)
             END IF
          END IF
       ENDDO
    ENDDO


    ! Sorties eventuelles.
    !----------------------

    IF(prt_level.GE.10) THEN
       DO k=1,nlay
          DO ig=1,ngrid
             PRINT*,'th_dv2 ig k gamma entr detr ua ue va ve wud wvd masse',ig,k,gamma(ig,k), &
                  &   entr(ig,k),detr(ig,k),ua(ig,k),ue(ig,k),va(ig,k),ve(ig,k),&
                  &   wud(ig,k),wvd(ig,k),wud(ig,k+1),wvd(ig,k+1), &
                  &   masse(ig,k)
          ENDDO
       ENDDO
    ENDIF
    !
    !IF (nlarga0>0) THEN
    !   PRINT*,'WARNING !!!!!! DANS THERMCELL_DV2 '
    !   PRINT*,nlarga0,' points pour lesquels laraga=0. dans un thermique'
    !   PRINT*,'Il faudrait decortiquer ces points'
    !ENDIF

    RETURN
  END SUBROUTINE thermcell_dv2

  !
  ! $Id: thermcell_flux2.F90 1403 2010-07-01 09:02:53Z fairhead $
  !
  SUBROUTINE thermcell_flux2(&
       ngrid      , &! INTEGER      , INTENT(IN   ) :: ngrid
       nlay       , &! INTEGER      , INTENT(IN   ) :: nlay
       ptimestep  , &! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
       masse      , &! REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
       lalim      , &! INTEGER      , INTENT(IN   ) :: lalim(ngrid)
       lmax       , &! INTEGER      , INTENT(IN   ) :: lmax(ngrid)
       alim_star  , &! REAL(KIND=r8), INTENT(IN   ) :: alim_star(ngrid,nlay)
       entr_star  , &! REAL(KIND=r8), INTENT(IN   ) :: entr_star(ngrid,nlay)
       detr_star  , &! REAL(KIND=r8), INTENT(IN   ) :: detr_star(ngrid,nlay)
       f          , &! REAL(KIND=r8), INTENT(IN   ) :: f(ngrid)
       rhobarz    , &! REAL(KIND=r8), INTENT(IN   ) :: rhobarz(ngrid,nlay)
       zw2        , &! REAL(KIND=r8), INTENT(OUT  ) :: zw2(ngrid,nlay+1)
       fm         , &! REAL(KIND=r8), INTENT(OUT  ) :: fm(ngrid,nlay+1)
       entr       , &! REAL(KIND=r8), INTENT(OUT  ) :: entr(ngrid,nlay)
       detr       , &! REAL(KIND=r8), INTENT(OUT  ) :: detr(ngrid,nlay)
       lunout1    , &! integer      , INTENT(IN   ) :: lunout1
       igout        )! integer      , INTENT(INOUT) :: igout


    !IM 060508    &       detr,zqla,zmax,lev_out,lunout,igout)


    !---------------------------------------------------------------------------
    !thermcell_flux: deduction des flux
    !---------------------------------------------------------------------------

    IMPLICIT NONE
    !#include "iniprint.h"
    !#include "thermcell.h"

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    REAL(KIND=r8), INTENT(IN   ) :: ptimestep
    REAL(KIND=r8), INTENT(IN   ) :: masse(ngrid,nlay)
    INTEGER      , INTENT(IN   ) :: lalim(ngrid)
    INTEGER      , INTENT(IN   ) :: lmax(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: alim_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: entr_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: detr_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: f(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: rhobarz(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: fm(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: entr(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: detr(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zw2(ngrid,nlay+1)
    INTEGER      , INTENT(IN   ) :: lunout1
    INTEGER      , INTENT(INOUT) :: igout

    INTEGER :: ncorecfm1
    INTEGER :: ncorecfm2
    INTEGER :: ncorecfm3
    INTEGER :: ncorecalpha
    INTEGER :: ncorecfm4
    INTEGER :: ncorecfm5
    INTEGER :: ncorecfm6
    INTEGER :: ncorecfm7
    INTEGER :: ncorecfm8

!    REAL(KIND=r8) :: zmax(ngrid)

    REAL(KIND=r8)    :: zfm

    INTEGER :: lout

    REAL(KIND=r8)    :: f_old
    REAL(KIND=r8)    :: ddd0
    REAL(KIND=r8)    :: eee0
    REAL(KIND=r8)    :: ddd
    REAL(KIND=r8)    :: eee
    REAL(KIND=r8)    :: zzz
    INTEGER :: ig
    INTEGER :: l

    REAL(KIND=r8) , PARAMETER   ::    fomass_max=0.5_r8

    REAL(KIND=r8) , PARAMETER   ::    alphamax=0.7_r8

    LOGICAL :: check_debug
    LOGICAL :: labort_gcm

!    CHARACTER (len=20) :: modname='thermcell_flux2'
    CHARACTER (len=80) :: abort_message

 fm=0.0_r8
 entr=0.0_r8
 detr=0.0_r8
 zfm=0.0_r8
 f_old=0.0_r8
 ddd0=0.0_r8
 eee0=0.0_r8
 ddd=0.0_r8
 eee=0.0_r8
 zzz=0.0_r8
    ncorecfm1=0
    ncorecfm2=0
    ncorecfm3=0
    ncorecfm4=0
    ncorecfm5=0
    ncorecfm6=0
    ncorecfm7=0
    ncorecfm8=0
    ncorecalpha=0

    !initialisation
    fm(:,:)=0.0_r8

    IF (prt_level.GE.10) THEN
       WRITE(lunout1,*) 'Dans thermcell_flux 0'
       WRITE(lunout1,*) 'flux base ',f(igout)
       WRITE(lunout1,*) 'lmax ',lmax(igout)
       WRITE(lunout1,*) 'lalim ',lalim(igout)
       WRITE(lunout1,*) 'ig= ',igout
       WRITE(lunout1,*) ' l E*    A*     D*  '
       WRITE(lunout1,'(i4,3e15.5)') (l,entr_star(igout,l),alim_star(igout,l),detr_star(igout,l) &
            &    ,l=1,lmax(igout))
    ENDIF


    !-------------------------------------------------------------------------
    ! Verification de la nullite des entrainement et detrainement au dessus
    ! de lmax(ig)
    ! Active uniquement si check_debug=.true. ou prt_level>=10
    !-------------------------------------------------------------------------

    check_debug=.FALSE..OR.prt_level>=10

    IF (check_debug) THEN
       DO l=1,nlay
          DO ig=1,ngrid
             IF (l.LE.lmax(ig)) THEN
                IF (entr_star(ig,l).GT.1.0_r8) THEN
                   PRINT*,'WARNING thermcell_flux 1 ig,l,lmax(ig)',ig,l,lmax(ig)
                   PRINT*,'entr_star(ig,l)',entr_star(ig,l)
                   PRINT*,'alim_star(ig,l)',alim_star(ig,l)
                   PRINT*,'detr_star(ig,l)',detr_star(ig,l)
                ENDIF
             ELSE
                IF (ABS(entr_star(ig,l))+ABS(alim_star(ig,l))+ABS(detr_star(ig,l)).GT.0.0_r8) THEN
                   PRINT*,'cas 1 : ig,l,lmax(ig)',ig,l,lmax(ig)
                   PRINT*,'entr_star(ig,l)',entr_star(ig,l)
                   PRINT*,'alim_star(ig,l)',alim_star(ig,l)
                   PRINT*,'detr_star(ig,l)',detr_star(ig,l)
                   abort_message = ''
                   labort_gcm=.TRUE.
                   STOP
                   !                    CALL abort_gcm (modname,abort_message,1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    !-------------------------------------------------------------------------
    ! Multiplication par le flux de masse issu de la femreture
    !-------------------------------------------------------------------------

    DO l=1,nlay
       DO ig=1,ngrid
          entr(ig,l)=f(ig)*(entr_star(ig,l)+alim_star(ig,l))
          detr(ig,l)=f(ig)*detr_star(ig,l)
       END DO
    ENDDO

    IF (prt_level.GE.10) THEN
       WRITE(lunout1,*) 'Dans thermcell_flux 1'
       WRITE(lunout1,*) 'flux base ',f(igout)
       WRITE(lunout1,*) 'lmax ',lmax(igout)
       WRITE(lunout1,*) 'lalim ',lalim(igout)
       WRITE(lunout1,*) 'ig= ',igout
       WRITE(lunout1,*) ' l   E    D     W2'
       WRITE(lunout1,'(i4,3e15.5)') (l,entr(igout,l),detr(igout,l) &
            &    ,zw2(igout,l+1),l=1,lmax(igout))
    ENDIF

    fm(:,1)=0.0_r8
    DO l=1,nlay
       DO ig=1,ngrid
          IF (l.LT.lmax(ig)) THEN
             fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
          ELSEIF(l.EQ.lmax(ig)) THEN
             fm(ig,l+1)=0.0_r8
             detr(ig,l)=fm(ig,l)+entr(ig,l)
          ELSE
             fm(ig,l+1)=0.0_r8
          ENDIF
       ENDDO
    ENDDO



    ! Test provisoire : pour comprendre pourquoi on corrige plein de fois 
    ! le cas fm6, on commence par regarder une premiere fois avant les
    ! autres corrections.

    DO l=1,nlay
       DO ig=1,ngrid
          IF (detr(ig,l).GT.fm(ig,l)) THEN
             ncorecfm8=ncorecfm8+1
             !              igout=ig
          ENDIF
       ENDDO
    ENDDO

    !      if (prt_level.ge.10) &
    !    &    call printflux(ngrid,nlay,lunout1,igout,f,lmax,lalim, &
    !    &    ptimestep,masse,entr,detr,fm,'2  ')



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FH Version en cours de test;
    ! par rapport a thermcell_flux, on fait une grande boucle sur "l"
    ! et on modifie le flux avec tous les contrï¿½les appliques d'affilee
    ! pour la meme couche
    ! Momentanement, on duplique le calcule du flux pour pouvoir comparer
    ! les flux avant et apres modif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO l=1,nlay

       DO ig=1,ngrid
          IF (l.LT.lmax(ig)) THEN
             fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
          ELSEIF(l.EQ.lmax(ig)) THEN
             fm(ig,l+1)=0.0_r8
             detr(ig,l)=fm(ig,l)+entr(ig,l)
          ELSE
             fm(ig,l+1)=0.0_r8
          ENDIF
       ENDDO


       !-------------------------------------------------------------------------
       ! Verification de la positivite des flux de masse
       !-------------------------------------------------------------------------

       !     do l=1,nlay
       DO ig=1,ngrid
          IF (fm(ig,l+1).LT.0.0_r8) THEN
             !              print*,'fm1<0',l+1,lmax(ig),fm(ig,l+1)
             ncorecfm1=ncorecfm1+1
             fm(ig,l+1)=fm(ig,l)
             detr(ig,l)=entr(ig,l)
          ENDIF
       ENDDO
       !     enddo

       IF (prt_level.GE.10) &
            &   WRITE(lunout1,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
            &     entr(igout,l),detr(igout,l),fm(igout,l+1)

       !-------------------------------------------------------------------------
       !Test sur fraca croissant
       !-------------------------------------------------------------------------
       IF (iflag_thermals_optflux==0) THEN 
          !     do l=1,nlay
          DO ig=1,ngrid
             IF (l.GE.lalim(ig).AND.l.LE.lmax(ig) &
                  &    .AND.(zw2(ig,l+1).GT.1.e-10_r8).AND.(zw2(ig,l).GT.1.e-10_r8) ) THEN
                !  zzz est le flux en l+1 a frac constant
                zzz=fm(ig,l)*rhobarz(ig,l+1)*zw2(ig,l+1)  &
                     &                          /(rhobarz(ig,l)*zw2(ig,l))
                IF (fm(ig,l+1).GT.zzz) THEN
                   detr(ig,l)=detr(ig,l)+fm(ig,l+1)-zzz
                   fm(ig,l+1)=zzz
                   ncorecfm4=ncorecfm4+1
                ENDIF
             ENDIF
          ENDDO
          !     enddo
       ENDIF

       IF (prt_level.GE.10) &
            &   WRITE(lunout1,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
            &     entr(igout,l),detr(igout,l),fm(igout,l+1)


       !-------------------------------------------------------------------------
       !test sur flux de masse croissant
       !-------------------------------------------------------------------------
       IF (iflag_thermals_optflux==0) THEN
          !     do l=1,nlay
          DO ig=1,ngrid
             IF ((fm(ig,l+1).GT.fm(ig,l)).AND.(l.GT.lalim(ig))) THEN
                f_old=fm(ig,l+1)
                fm(ig,l+1)=fm(ig,l)
                detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
                ncorecfm5=ncorecfm5+1
             ENDIF
          ENDDO
          !     enddo
       ENDIF

       IF (prt_level.GE.10) &
            &   WRITE(lunout1,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
            &     entr(igout,l),detr(igout,l),fm(igout,l+1)

       !fin 1.eq.0
       !-------------------------------------------------------------------------
       !detr ne peut pas etre superieur a fm
       !-------------------------------------------------------------------------

       IF(1.EQ.1) THEN

          !     do l=1,nlay



          labort_gcm=.FALSE.
          DO ig=1,ngrid
             IF (entr(ig,l)<0.0_r8) THEN
                labort_gcm=.TRUE.
                igout=ig
                lout=l
             ENDIF
          ENDDO

          IF (labort_gcm) THEN
             PRINT*,'N1 ig,l,entr',igout,lout,entr(igout,lout)
             abort_message = 'entr negatif'
             STOP
             !            CALL abort_gcm (modname,abort_message,1)
          ENDIF

          DO ig=1,ngrid
             IF (detr(ig,l).GT.fm(ig,l)) THEN
                ncorecfm6=ncorecfm6+1
                detr(ig,l)=fm(ig,l)
                entr(ig,l)=fm(ig,l+1)

                ! Dans le cas ou on est au dessus de la couche d'alimentation et que le
                ! detrainement est plus fort que le flux de masse, on stope le thermique.
                !test:on commente
                !               if (l.gt.lalim(ig)) then
                !                  lmax(ig)=l
                !                  fm(ig,l+1)=0.
                !                  entr(ig,l)=0.
                !               else
                !                  ncorecfm7=ncorecfm7+1
                !               endif
             ENDIF

             IF(l.GT.lmax(ig)) THEN
                detr(ig,l)=0.0_r8
                fm(ig,l+1)=0.0_r8
                entr(ig,l)=0.0_r8
             ENDIF
          ENDDO

          labort_gcm=.FALSE.
          DO ig=1,ngrid
             IF (entr(ig,l).LT.0.0_r8) THEN
                labort_gcm=.TRUE.
                igout=ig
             ENDIF
          ENDDO
          IF (labort_gcm) THEN
             ig=igout
             PRINT*,'ig,l,lmax(ig)',ig,l,lmax(ig)
             PRINT*,'entr(ig,l)',entr(ig,l)
             PRINT*,'fm(ig,l)',fm(ig,l)
             abort_message = 'probleme dans thermcell flux'
             STOP
             !            CALL abort_gcm (modname,abort_message,1)
          ENDIF


          !     enddo
       ENDIF


       IF (prt_level.GE.10) &
            &   WRITE(lunout1,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
            &     entr(igout,l),detr(igout,l),fm(igout,l+1)

       !-------------------------------------------------------------------------
       !fm ne peut pas etre negatif
       !-------------------------------------------------------------------------

       !     do l=1,nlay
       DO ig=1,ngrid
          IF (fm(ig,l+1).LT.0.0_r8) THEN
             detr(ig,l)=detr(ig,l)+fm(ig,l+1)
             fm(ig,l+1)=0.0_r8
             ncorecfm2=ncorecfm2+1
          ENDIF
       ENDDO

       labort_gcm=.FALSE.
       DO ig=1,ngrid
          IF (detr(ig,l).LT.0.0_r8) THEN
             labort_gcm=.TRUE.
             igout=ig
          ENDIF
       ENDDO
       IF (labort_gcm) THEN
          ig=igout
          PRINT*,'cas 2 : ig,l,lmax(ig)',ig,l,lmax(ig)
          PRINT*,'detr(ig,l)',detr(ig,l)
          PRINT*,'fm(ig,l)',fm(ig,l)
          abort_message = 'probleme dans thermcell flux'
          STOP
          !               CALL abort_gcm (modname,abort_message,1)
       ENDIF
       !    enddo

       IF (prt_level.GE.10) &
            &   WRITE(lunout1,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
            &     entr(igout,l),detr(igout,l),fm(igout,l+1)

       !-----------------------------------------------------------------------
       !la fraction couverte ne peut pas etre superieure a 1            
       !-----------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! FH Partie a revisiter.
       ! Il semble qu'etaient codees ici deux optiques dans le cas
       ! F/ (rho *w) > 1
       ! soit limiter la hauteur du thermique en considerant que c'est 
       ! la derniere chouche, soit limiter F a rho w.
       ! Dans le second cas, il faut en fait limiter a un peu moins
       ! que ca parce qu'on a des 1 / ( 1 -alpha) un peu plus loin
       ! dans thermcell_main et qu'il semble de toutes facons deraisonable
       ! d'avoir des fractions de 1..
       ! Ci dessous, et dans l'etat actuel, le premier des  deux if est
       ! sans doute inutile.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !    do l=1,nlay
       DO ig=1,ngrid
          IF (zw2(ig,l+1).GT.1.e-10_r8) THEN
             zfm=rhobarz(ig,l+1)*zw2(ig,l+1)*alphamax
             IF ( fm(ig,l+1) .GT. zfm) THEN
                f_old=fm(ig,l+1)
                fm(ig,l+1)=zfm
                !             zw2(ig,l+1)=0.
                !             zqla(ig,l+1)=0.
                detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
                !             lmax(ig)=l+1
                !             zmax(ig)=zlev(ig,lmax(ig))
                !             print*,'alpha>1',l+1,lmax(ig)
                ncorecalpha=ncorecalpha+1
             ENDIF
          ENDIF
       ENDDO
       !    enddo
       !


       IF (prt_level.GE.10) &
            &   WRITE(lunout1,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
            &     entr(igout,l),detr(igout,l),fm(igout,l+1)

       ! Fin de la grande boucle sur les niveaux verticaux
    ENDDO

    !      if (prt_level.ge.10) &
    !    &    call printflux(ngrid,nlay,lunout1,igout,f,lmax,lalim, &
    !    &    ptimestep,masse,entr,detr,fm,'8  ')


    !-----------------------------------------------------------------------
    ! On fait en sorte que la quantite totale d'air entraine dans le 
    ! panache ne soit pas trop grande comparee a la masse de la maille
    !-----------------------------------------------------------------------

    IF (1.EQ.1) THEN
       labort_gcm=.FALSE.
       DO l=1,nlay-1
          DO ig=1,ngrid
             eee0=entr(ig,l)
             ddd0=detr(ig,l)
             eee=entr(ig,l)-masse(ig,l)*fomass_max/ptimestep
             ddd=detr(ig,l)-eee
             IF (eee.GT.0.0_r8) THEN
                ncorecfm3=ncorecfm3+1
                entr(ig,l)=entr(ig,l)-eee
                IF ( ddd.GT.0.0_r8) THEN
                   !   l'entrainement est trop fort mais l'exces peut etre compense par une
                   !   diminution du detrainement)
                   detr(ig,l)=ddd
                ELSE
                   !   l'entrainement est trop fort mais l'exces doit etre compense en partie
                   !   par un entrainement plus fort dans la couche superieure
                   IF(l.EQ.lmax(ig)) THEN
                      detr(ig,l)=fm(ig,l)+entr(ig,l)
                   ELSE
                      IF(l.GE.lmax(ig).AND.0.EQ.1) THEN
                         igout=ig
                         lout=l
                         labort_gcm=.TRUE.
                      ENDIF
                      entr(ig,l+1)=entr(ig,l+1)-ddd
                      detr(ig,l)=0.0_r8
                      fm(ig,l+1)=fm(ig,l)+entr(ig,l)
                      detr(ig,l)=0.0_r8
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF (labort_gcm) THEN
          ig=igout
          l=lout
          PRINT*,'ig,l',ig,l
          PRINT*,'eee0',eee0
          PRINT*,'ddd0',ddd0
          PRINT*,'eee',eee
          PRINT*,'ddd',ddd
          PRINT*,'entr',entr(ig,l)
          PRINT*,'detr',detr(ig,l)
          PRINT*,'masse',masse(ig,l)
          PRINT*,'fomass_max',fomass_max
          PRINT*,'masse(ig,l)*fomass_max/ptimestep',masse(ig,l)*fomass_max/ptimestep
          PRINT*,'ptimestep',ptimestep
          PRINT*,'lmax(ig)',lmax(ig)
          PRINT*,'fm(ig,l+1)',fm(ig,l+1)
          PRINT*,'fm(ig,l)',fm(ig,l)
          abort_message = 'probleme dans thermcell_flux'
          STOP
          !                         CALL abort_gcm (modname,abort_message,1)
       ENDIF
    ENDIF
    !                  
    !              ddd=detr(ig)-entre
    !on s assure que tout s annule bien en zmax
    DO ig=1,ngrid
       fm(ig,lmax(ig)+1)=0.0_r8
       entr(ig,lmax(ig))=0.0_r8
       detr(ig,lmax(ig))=fm(ig,lmax(ig))+entr(ig,lmax(ig))
    ENDDO

    !-----------------------------------------------------------------------
    ! Impression du nombre de bidouilles qui ont ete necessaires
    !-----------------------------------------------------------------------

    !IM 090508 beg
    !     if (ncorecfm1+ncorecfm2+ncorecfm3+ncorecfm4+ncorecfm5+ncorecalpha > 0 ) then
    !
    !         print*,'PB thermcell : on a du coriger ',ncorecfm1,'x fm1',&
    !   &     ncorecfm2,'x fm2',ncorecfm3,'x fm3 et', &
    !   &     ncorecfm4,'x fm4',ncorecfm5,'x fm5 et', &
    !   &     ncorecfm6,'x fm6', &
    !   &     ncorecfm7,'x fm7', &
    !   &     ncorecfm8,'x fm8', &
    !   &     ncorecalpha,'x alpha'
    !     endif
    !IM 090508 end

    !      if (prt_level.ge.10) &
    !    &    call printflux(ngrid,nlay,lunout1,igout,f,lmax,lalim, &
    !    &    ptimestep,masse,entr,detr,fm,'fin')


    RETURN
  END SUBROUTINE thermcell_flux2



  SUBROUTINE thermcell_env( &
       ngrid    ,&! INTEGER      , INTENT(IN   ) :: ngrid
       nlay     ,&! INTEGER      , INTENT(IN   ) :: nlay
       po       ,&! REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
       pt       ,&! REAL(KIND=r8), INTENT(IN   ) :: pt(ngrid,nlay)
       pu       ,&! REAL(KIND=r8), INTENT(IN   ) :: pu(ngrid,nlay)
       pv       ,&! REAL(KIND=r8), INTENT(IN   ) :: pv(ngrid,nlay)
       pplay    ,&! REAL(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
       pplev    ,&! REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
       zo       ,&! REAL(KIND=r8), INTENT(OUT  ) :: zo(ngrid,nlay)
       zh       ,&! REAL(KIND=r8), INTENT(OUT  ) :: zh(ngrid,nlay)
       zl       ,&! REAL(KIND=r8), INTENT(OUT  ) :: zl(ngrid,nlay)
       ztv      ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztv(ngrid,nlay)
       zthl     ,&! REAL(KIND=r8), INTENT(OUT  ) :: zthl(ngrid,nlay)
       zu       ,&! REAL(KIND=r8), INTENT(OUT  ) :: zu(ngrid,nlay)
       zv       ,&! REAL(KIND=r8), INTENT(OUT  ) :: zv(ngrid,nlay)
       zpspsk   ,&! REAL(KIND=r8), INTENT(OUT  ) :: zpspsk(ngrid,nlay) 
       pqsat      )! REAL(KIND=r8), INTENT(OUT  ) :: pqsat(ngrid,nlay)

    !--------------------------------------------------------------
    !thermcell_env: calcule les caracteristiques de l environnement
    !necessaires au calcul des proprietes dans le thermique
    !--------------------------------------------------------------

    IMPLICIT NONE

    !#include "YOMCST.h"
    !#include "YOETHF.h"
    !#include "FCTTRE.h"      
    !#include "iniprint.h"

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pt(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pu(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pv(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pplay(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: zo(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zh(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zl(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ztv(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zthl(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zu(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zv(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zpspsk(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: pqsat(ngrid,nlay)

    INTEGER       :: ig,ll

!    REAL(KIND=r8) :: dqsat_dT
    REAL(KIND=r8) :: RLvCp

    LOGICAL       :: mask(ngrid,nlay)


    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Initialisations :
    !------------------
    zo=0.0_r8
    zh=0.0_r8
    zl=0.0_r8
    ztv=0.0_r8
    zthl=0.0_r8
    zu=0.0_r8
    zv=0.0_r8
    zpspsk=0.0_r8
    pqsat=0.0_r8
    RLvCp=0.0_r8

    mask(:,:)=.TRUE.
    RLvCp = RLVTT/RCPD

    !
    ! calcul des caracteristiques de l environnement
    DO  ll=1,nlay
       DO ig=1,ngrid
          zo(ig,ll)=po(ig,ll)
          zl(ig,ll)=0.0_r8
          zh(ig,ll)=pt(ig,ll)
       ENDDO
    ENDDO
    !
    !
    ! Condensation :
    !---------------
    ! Calcul de l'humidite a saturation et de la condensation
    DO ll=1,nlay
       CALL thermcell_qsat( &
            ngrid                    ,&!  ngrid     ,&! INTEGER     , INTENT(IN   ) :: ngrid
            mask  (1:ngrid,ll)       ,&!  active   ,&! LOGICAL      , INTENT(IN   ) :: active(ngrid)
            pplev (1:ngrid,ll)       ,&!  pplev    ,&! REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid)
            pt    (1:ngrid,ll)       ,&!  ztemp    ,&! REAL(KIND=r8), INTENT(IN   ) :: ztemp(ngrid)
            po    (1:ngrid,ll)       ,&!  zqta     ,&! REAL(KIND=r8), INTENT(IN   ) :: zqta(ngrid)
            pqsat (1:ngrid,ll)        )!  zqsat    ) ! REAL(KIND=r8), INTENT(OUT  ) :: zqsat(ngrid)
    END DO
    DO ll=1,nlay
       DO ig=1,ngrid
          zl(ig,ll) = MAX(0.0_r8,po(ig,ll)-pqsat(ig,ll))
          zh(ig,ll) = pt(ig,ll)+RLvCp*zl(ig,ll)         !   T = Tl + Lv/Cp ql
          zo(ig,ll) = po(ig,ll)-zl(ig,ll)
       ENDDO
    ENDDO
    !
    !
    !-----------------------------------------------------------------------

   ! IF (prt_level.GE.1) PRINT*,'0 OK convect8'

    DO ll=1,nlay
       DO ig=1,ngrid
          zpspsk(ig,ll)=(pplay(ig,ll)/100000.0_r8)**RKAPPA
          zu(ig,ll)=pu(ig,ll)
          zv(ig,ll)=pv(ig,ll)
          !attention zh est maintenant le profil de T et plus le profil de theta !
          ! Quelle horreur ! A eviter.
          !
          !   T-> Theta
          ztv(ig,ll)=zh(ig,ll)/zpspsk(ig,ll)
          !Theta_v
          ztv(ig,ll)=ztv(ig,ll)*(1.0_r8+RETV*(zo(ig,ll))-zl(ig,ll))
          !Thetal
          zthl(ig,ll)=pt(ig,ll)/zpspsk(ig,ll)
          !            
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE thermcell_env

  !
  ! $Id: thermcell_plume.F90 1503 2011-03-23 11:57:52Z idelkadi $
  !
  SUBROUTINE thermcell_plume( &
       ngrid         ,&! INTEGER      , INTENT(IN   ) :: ngrid
       nlay          ,&! INTEGER      , INTENT(IN   ) :: nlay
!       ptimestep     ,&! REAL(KIND=r8), INTENT(IN   ) :: ptimestep
       ztv           ,&! REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
       zthl          ,&! REAL(KIND=r8), INTENT(IN   ) :: zthl(ngrid,nlay)
       po            ,&! REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
       zl            ,&! REAL(KIND=r8), INTENT(IN   ) :: zl(ngrid,nlay)
       rhobarz       ,&! REAL(KIND=r8), INTENT(IN   ) :: rhobarz(ngrid,nlay)
       zlev          ,&! REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
       pplev         ,&! REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
       pphi          ,&! REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
       zpspsk        ,&! REAL(KIND=r8), INTENT(IN   ) :: zpspsk(ngrid,nlay)
       alim_star     ,&! REAL(KIND=r8), INTENT(OUT  ) :: alim_star(ngrid,nlay)
       alim_star_tot ,&! REAL(KIND=r8), INTENT(OUT  ) :: alim_star_tot(ngrid)
       lalim         ,&! INTEGER      , INTENT(OUT  ) :: lalim(ngrid)
       f0            ,&! REAL(KIND=r8), INTENT(IN   ) :: f0(ngrid)
       detr_star     ,&! REAL(KIND=r8), INTENT(OUT  ) :: detr_star(ngrid,nlay)
       entr_star     ,&! REAL(KIND=r8), INTENT(OUT  ) :: entr_star(ngrid,nlay)
       f_star        ,&! REAL(KIND=r8), INTENT(OUT  ) :: f_star(ngrid,nlay+1)
       csc           ,&! REAL(KIND=r8), INTENT(OUT  ) :: csc(ngrid,nlay)
       ztva          ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztva(ngrid,nlay)
       ztla          ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztla(ngrid,nlay)
       zqla          ,&! REAL(KIND=r8), INTENT(OUT  ) :: zqla(ngrid,nlay)
       zqta          ,&! REAL(KIND=r8), INTENT(OUT  ) :: zqta(ngrid,nlay)
       zha           ,&! REAL(KIND=r8), INTENT(OUT  ) :: zha(ngrid,nlay)
       zw2           ,&! REAL(KIND=r8), INTENT(OUT  ) :: zw2(ngrid,nlay+1)
       w_est         ,&! REAL(KIND=r8), INTENT(OUT  ) :: w_est(ngrid,nlay+1)
       ztva_est      ,&! REAL(KIND=r8), INTENT(OUT  ) :: ztva_est(ngrid,nlay)
       zqsatth       ,&! REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(ngrid,nlay)
       lmix          ,&! INTEGER      , INTENT(OUT  ) :: lmix(ngrid)
       lmix_bis      ,&! INTEGER      , INTENT(OUT  ) :: lmix_bis(ngrid)
       linter         )! REAL(KIND=r8), INTENT(OUT  ) :: linter(ngrid)

    !--------------------------------------------------------------------------
    !thermcell_plume: calcule les valeurs de qt, thetal et w dans l ascendance
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    !#include "YOMCST.h"
    !#include "YOETHF.h"
    !#include "FCTTRE.h"
    !#include "iniprint.h"
    !#include "thermcell.h"

    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
!    REAL(KIND=r8), INTENT(IN   ) :: ptimestep
    REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zthl(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: po(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zl(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: rhobarz(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zpspsk(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: alim_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: alim_star_tot(ngrid)
    INTEGER      , INTENT(OUT  ) :: lalim(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: f0(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: detr_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: entr_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: f_star(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: csc(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ztva(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: ztla(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqla(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqta(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zha(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zw2(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: w_est(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(OUT  ) :: ztva_est(ngrid,nlay)
    REAL(KIND=r8), INTENT(OUT  ) :: zqsatth(ngrid,nlay)
    INTEGER      , INTENT(OUT  ) :: lmix(ngrid)
    INTEGER      , INTENT(OUT  ) :: lmix_bis(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: linter(ngrid)


    INTEGER :: nbpb
!    REAL(KIND=r8)    :: zcon2(ngrid)



!    REAL(KIND=r8)    :: coefc
    REAL(KIND=r8)    :: detr(ngrid,nlay)
    REAL(KIND=r8)    :: entr(ngrid,nlay)


    REAL(KIND=r8)    :: wa_moy(ngrid,nlay+1)

    REAL(KIND=r8)    :: zqla_est(ngrid,nlay)
    REAL(KIND=r8)    :: zta_est(ngrid,nlay)
!    REAL(KIND=r8)    :: zdw2
!    REAL(KIND=r8)    :: zw2modif
    REAL(KIND=r8)    :: zeps

    REAL(KIND=r8)    :: wmaxa(ngrid)

    INTEGER :: ig
    INTEGER :: l
!    INTEGER :: k

    REAL(KIND=r8)    :: zdz
    REAL(KIND=r8)    :: zfact
    REAL(KIND=r8)    :: zbuoy
    REAL(KIND=r8)    :: zalpha
    REAL(KIND=r8)    :: zdrag
!    REAL(KIND=r8)    :: zcor
!    REAL(KIND=r8)    :: zdelta
!    REAL(KIND=r8)    :: zcvm5
!    REAL(KIND=r8)    :: qlbef
!    REAL(KIND=r8)    :: Tbef
    REAL(KIND=r8)    :: qsatbef
!    REAL(KIND=r8)    :: dqsat_dT
!    REAL(KIND=r8)    :: DT
!    REAL(KIND=r8)    :: num
!    REAL(KIND=r8)    :: denom
!    REAL(KIND=r8)    :: REPS
    REAL(KIND=r8)    :: RLvCp
    LOGICAL :: Zsat
    LOGICAL :: active(ngrid)
    LOGICAL :: activetmp(ngrid)
    REAL(KIND=r8)    :: fact_gamma
    REAL(KIND=r8)    :: fact_epsilon
    REAL(KIND=r8)    :: fact_gamma2
!    REAL(KIND=r8)    :: c2(ngrid,nlay)
    REAL(KIND=r8)    :: a1,m

    REAL(KIND=r8)    :: zw2fact
    REAL(KIND=r8)    :: expa



    ! Initialisation

    alim_star=0.0_r8;
    alim_star_tot=0.0_r8;
    lalim=0;
    detr_star=0.0_r8;
    entr_star=0.0_r8;
    f_star=0.0_r8
    csc=0.0_r8;
    ztva=0.0_r8;
    ztla=0.0_r8;
    zqla=0.0_r8;
    zqta=0.0_r8;
    zha=0.0_r8;
    zw2=0.0_r8
    w_est=0.0_r8
    ztva_est=0.0_r8;
    zqsatth=0.0_r8;
    lmix=0;
    lmix_bis=0;
    linter=0.0_r8;
    nbpb=0;
  detr=0.0_r8;
  entr=0.0_r8;
  wa_moy=0.0_r8
  zqla_est=0.0_r8;
  zta_est=0.0_r8;
  zeps=0.0_r8;
  wmaxa=0.0_r8;

  ig=0
  l=0
  zdz=0.0_r8;
  zfact=0.0_r8;
  zbuoy=0.0_r8;
  zalpha=0.0_r8;
  zdrag=0.0_r8;
  qsatbef=0.0_r8;
  RLvCp=0.0_r8;
  fact_gamma=0.0_r8;
  fact_epsilon=0.0_r8;
  fact_gamma2=0.0_r8;
  a1=0.0_r8;m=0.0_r8;

  zw2fact=0.0_r8;
  expa=0.0_r8;



    Zsat=.FALSE.
    ! Initialisation
    RLvCp = RLVTT/RCPD


    fact_epsilon=0.002_r8
    a1=2.0_r8/3.0_r8
    fact_gamma=0.9_r8
    zfact=fact_gamma/(1+fact_gamma)
    fact_gamma2=zfact
    expa=0.0_r8


    ! Initialisations des variables reeles
    IF (1==1) THEN
       ztva(:,:)=ztv(:,:)
       ztva_est(:,:)=ztva(:,:)
       ztla(:,:)=zthl(:,:)
       zqta(:,:)=po(:,:)
       zha(:,:) = ztva(:,:)
    ELSE
       ztva(:,:)=0.0_r8
       ztva_est(:,:)=0.0_r8
       ztla(:,:)=0.0_r8
       zqta(:,:)=0.0_r8
       zha(:,:) =0.0_r8
    ENDIF

    zqla_est(:,:)=0.0_r8
    zqsatth(:,:)=0.0_r8
    zqla(:,:)=0.0_r8
    detr_star(:,:)=0.0_r8
    entr_star(:,:)=0.0_r8
    alim_star(:,:)=0.0_r8
    alim_star_tot(:)=0.0_r8
    csc(:,:)=0.0_r8
    detr(:,:)=0.0_r8
    entr(:,:)=0.0_r8
    zw2(:,:)=0.0_r8
    w_est(:,:)=0.0_r8
    f_star(:,:)=0.0_r8
    wa_moy(:,:)=0.0_r8
    linter(:)=1.0_r8
    linter(:)=1.0_r8

    ! Initialisation des variables entieres
    lmix(:)=1
    lmix_bis(:)=2
    wmaxa(:)=0.0_r8
    lalim(:)=1

    !-------------------------------------------------------------------------
    ! On ne considere comme actif que les colonnes dont les deux premieres
    ! couches sont instables.
    !-------------------------------------------------------------------------
    active(:)=ztv(:,1)>ztv(:,2)

    !-------------------------------------------------------------------------
    ! Definition de l'alimentation a l'origine dans thermcell_init
    !-------------------------------------------------------------------------
    DO l=1,nlay-1
       DO ig=1,ngrid
          IF (ztv(ig,l)> ztv(ig,l+1) .AND. ztv(ig,1)>=ztv(ig,l) ) THEN
             alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.0_r8)  &
                  &                       *SQRT(zlev(ig,l+1)) 
             lalim(ig)=l+1
             alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          ENDIF
       ENDDO
    ENDDO
    DO l=1,nlay
       DO ig=1,ngrid 
          IF (alim_star_tot(ig) > 1.e-10_r8 ) THEN
             alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
          ENDIF
       ENDDO
    ENDDO
    alim_star_tot(:)=1.0_r8


    !------------------------------------------------------------------------------
    ! Calcul dans la premiere couche
    ! On decide dans cette version que le thermique n'est actif que si la premiere
    ! couche est instable.
    ! Pourrait etre change si on veut que le thermiques puisse se dÃ©clencher
    ! dans une couche l>1
    !------------------------------------------------------------------------------
    DO ig=1,ngrid
       ! Le panache va prendre au debut les caracteristiques de l'air contenu
       ! dans cette couche.
       IF (active(ig)) THEN
          ztla(ig,1)=zthl(ig,1) 
          zqta(ig,1)=po(ig,1)
          zqla(ig,1)=zl(ig,1)
          !cr: attention, prise en compte de f*(1)=1
          f_star(ig,2)=alim_star(ig,1)
          zw2(ig,2)=2.0_r8*RG*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
               &                     *(zlev(ig,2)-zlev(ig,1))  &
               &                     *0.4_r8*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1))
          w_est(ig,2)=zw2(ig,2)
       ENDIF
    ENDDO
    !

    !==============================================================================
    !boucle de calcul de la vitesse verticale dans le thermique
    !==============================================================================
    DO l=2,nlay-1
       !==============================================================================


       ! On decide si le thermique est encore actif ou non
       ! AFaire : Il faut sans doute ajouter entr_star a alim_star dans ce test
       DO ig=1,ngrid
          active(ig)=active(ig) &
               &                 .AND. zw2(ig,l)>1.e-10_r8 &
               &                 .AND. f_star(ig,l)+alim_star(ig,l)>1.e-10_r8
       ENDDO



       ! Premier calcul de la vitesse verticale a partir de la temperature
       ! potentielle virtuelle
       !     if (1.eq.1) then
       !         w_est(ig,3)=zw2(ig,2)* &
       !    &      ((f_star(ig,2))**2) &
       !    &      /(f_star(ig,2)+alim_star(ig,2))**2+ &
       !    &      2.*RG*(ztva(ig,2)-ztv(ig,2))/ztv(ig,2) &
       !    &      *(zlev(ig,3)-zlev(ig,2))
       !      endif


       !---------------------------------------------------------------------------
       ! calcul des proprietes thermodynamiques et de la vitesse de la couche l
       ! sans tenir compte du detrainement et de l'entrainement dans cette
       ! couche
       ! Ici encore, on doit pouvoir ajouter entr_star (qui peut etre calculer
       ! avant) a l'alimentation pour avoir un calcul plus propre
       !---------------------------------------------------------------------------

       CALL thermcell_condens( &
            ngrid                  , & !    INTEGER      , INTENT(IN   )    :: ngrid
            active  (1:ngrid)      , & !    LOGICAL      , INTENT(IN   )    :: active(ngrid)
            zpspsk  (1:ngrid,l)    , & !    REAL(KIND=r8), INTENT(IN   )    :: zpspsk(ngrid)
            pplev   (1:ngrid,l)    , & !    REAL(KIND=r8), INTENT(IN   )    :: pplev (ngrid)
            ztla    (1:ngrid,l-1)  , & !    REAL(KIND=r8), INTENT(IN   )    :: ztla  (ngrid)
            zqta    (1:ngrid,l-1)  , & !    REAL(KIND=r8), INTENT(IN   )    :: zqta  (ngrid)
            zqla_est(1:ngrid,l)      ) !    REAL(KIND=r8), INTENT(INOUT)    :: zqla  (ngrid)

       DO ig=1,ngrid
          IF(active(ig)) THEN
             ztva_est(ig,l) = ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*zqla_est(ig,l)
             zta_est(ig,l)=ztva_est(ig,l)
             ztva_est(ig,l) = ztva_est(ig,l)/zpspsk(ig,l)
             ztva_est(ig,l) = ztva_est(ig,l)*(1.0_r8+RETV*(zqta(ig,l-1)  &
                  &      -zqla_est(ig,l))-zqla_est(ig,l))

             IF (1.EQ.0) THEN 
                !calcul de w_est sans prendre en compte le drag 
                w_est(ig,l+1)=zw2(ig,l)*  &
                     &                   ((f_star(ig,l))**2)  &
                     &                   /(f_star(ig,l)+alim_star(ig,l))**2+  &
                     &                   2.0_r8*RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)  &
                     &                   *(zlev(ig,l+1)-zlev(ig,l))
             ELSE

                zdz=zlev(ig,l+1)-zlev(ig,l)
                zalpha=f0(ig)*f_star(ig,l)/SQRT(w_est(ig,l))/rhobarz(ig,l)
                zbuoy=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
                zdrag=fact_epsilon/(zalpha**expa)
                zw2fact=zbuoy/zdrag*a1
                w_est(ig,l+1)=(w_est(ig,l)-zw2fact)*EXP(-2.0_r8*zdrag/(1+fact_gamma)*zdz) &
                     &    +zw2fact

             ENDIF

             IF (w_est(ig,l+1).LT.0.0_r8) THEN
                w_est(ig,l+1)=zw2(ig,l)
             ENDIF
          ENDIF
       ENDDO

       !-------------------------------------------------
       !calcul des taux d'entrainement et de detrainement
       !-------------------------------------------------

       DO ig=1,ngrid
          IF (active(ig)) THEN
             zdz=zlev(ig,l+1)-zlev(ig,l)
             zbuoy=RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)

             ! estimation de la fraction couverte par les thermiques
             zalpha=f0(ig)*f_star(ig,l)/SQRT(w_est(ig,l))/rhobarz(ig,l)

             !calcul de la soumission papier 
             ! Calcul  du taux d'entrainement entr_star (epsilon)
             entr_star(ig,l)=f_star(ig,l)*zdz * (  zfact * MAX(0.0_r8,  &     
                  &     a1*zbuoy/w_est(ig,l+1) &
                  &     - fact_epsilon/zalpha**expa  ) &
                  &     +0.0_r8 )

             !calcul du taux de detrainment (delta)
             !           detr_star(ig,l)=f_star(ig,l)*zdz * (                           &
             !     &      MAX(1.e-3, &
             !     &      -fact_gamma2*a1*zbuoy/w_est(ig,l+1)        &
             !     &      +0.01*(max(zqta(ig,l-1)-po(ig,l),0.)/(po(ig,l))/(w_est(ig,l+1)))**0.5    &    
             !     &     +0. ))

             m=0.5_r8

             detr_star(ig,l)=1.0_r8*f_star(ig,l)*zdz *                    &
                  &     MAX(5.e-4_r8,-fact_gamma2*a1*(1.0_r8/w_est(ig,l+1))*     &
                  ((1.0_r8-(1.0_r8-m)/(1.0_r8+70*zqta(ig,l-1)))*zbuoy        &
                  &     -40*(1.0_r8-m)*(MAX(zqta(ig,l-1)-po(ig,l),0.0_r8))/(1.0_r8+70*zqta(ig,l-1)) )   )

             !           detr_star(ig,l)=f_star(ig,l)*zdz * (                           &
             !     &      MAX(0.0_r8, &
             !     &      -fact_gamma2*a1*zbuoy/w_est(ig,l+1)        &
             !     &      +20*(max(zqta(ig,l-1)-po(ig,l),0.0_r8))**1*(zalpha/w_est(ig,l+1))**0.5_r8    &    
             !     &     +0.0_r8 ))


             ! En dessous de lalim, on prend le max de alim_star et entr_star pour
             ! alim_star et 0 sinon
             IF (l.LT.lalim(ig)) THEN
                alim_star(ig,l)=MAX(alim_star(ig,l),entr_star(ig,l))
                entr_star(ig,l)=0.0_r8
             ENDIF

             !attention test
             !        if (detr_star(ig,l).gt.(f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l))) then       
             !            detr_star(ig,l)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)
             !        endif
             ! Calcul du flux montant normalise
             f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
                  &              -detr_star(ig,l)

          ENDIF
       ENDDO

       !----------------------------------------------------------------------------
       !calcul de la vitesse verticale en melangeant Tl et qt du thermique
       !---------------------------------------------------------------------------
       DO ig=1,ngrid
          activetmp(ig)=active(ig) .AND. f_star(ig,l+1)>1.e-10_r8
       END DO
       DO ig=1,ngrid
          IF (activetmp(ig)) THEN 
             Zsat=.FALSE.
             ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
                  &            (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))  &
                  &            /(f_star(ig,l+1)+detr_star(ig,l))
             zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+  &
                  &            (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))  &
                  &            /(f_star(ig,l+1)+detr_star(ig,l))

          ENDIF
       ENDDO

       CALL thermcell_condens( &
            ngrid                   ,& ! INTEGER      , INTENT(IN   ) :: ngrid
            activetmp(1:ngrid)      ,& ! LOGICAL      , INTENT(IN   ) :: active(ngrid)
            zpspsk   (1:ngrid,l)    ,& ! REAL(KIND=r8), INTENT(IN   ) :: zpspsk(ngrid)
            pplev    (1:ngrid,l)    ,& ! REAL(KIND=r8), INTENT(IN   ) :: pplev (ngrid)
            ztla     (1:ngrid,l)    ,& ! REAL(KIND=r8), INTENT(IN   ) :: ztla  (ngrid)
            zqta     (1:ngrid,l)    ,& ! REAL(KIND=r8), INTENT(IN   ) :: zqta  (ngrid)
            zqla     (1:ngrid,l)     ) ! REAL(KIND=r8), INTENT(INOUT) :: zqla  (ngrid)


       DO ig=1,ngrid
          IF (activetmp(ig)) THEN
             ! on ecrit de maniere conservative (sat ou non)
             !          T = Tl +Lv/Cp ql
             ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
             ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
             !on rajoute le calcul de zha pour diagnostiques (temp potentielle)
             zha(ig,l) = ztva(ig,l)
             ztva(ig,l) = ztva(ig,l)*(1.0_r8+RETV*(zqta(ig,l)  &
                  &              -zqla(ig,l))-zqla(ig,l))

             !on ecrit zqsat 
             zqsatth(ig,l)=qsatbef  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !          zw2(ig,l+1)=&
             !     &                 zw2(ig,l)*(1-fact_epsilon/(1.0_r8+fact_gamma)*2.0_r8*(zlev(ig,l+1)-zlev(ig,l))) &
             !     &                 +2.0_r8*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
             !     &                 *1.0_r8/(1.0_r8+fact_gamma) &
             !     &                 *(zlev(ig,l+1)-zlev(ig,l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! La meme en plus modulaire :
             zbuoy=RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
             zdz=zlev(ig,l+1)-zlev(ig,l)


             zeps=(entr_star(ig,l)+alim_star(ig,l))/(f_star(ig,l)*zdz)

             !if (1==0) then
             !           zw2modif=zw2(ig,l)*(1-fact_epsilon/(1.0_r8+fact_gamma)*2.0_r8*zdz)
             !           zdw2=2.0_r8*zbuoy/(1.0_r8+fact_gamma)*zdz
             !           zw2(ig,l+1)=zw2modif+zdw2
             !else
             zdrag=fact_epsilon/(zalpha**expa)
             zw2fact=zbuoy/zdrag*a1
             zw2(ig,l+1)=(zw2(ig,l)-zw2fact)*EXP(-2.0_r8*zdrag/(1+fact_gamma)*zdz) &
                  &    +zw2fact
             !endif

          ENDIF
       ENDDO

       IF (prt_level.GE.20) PRINT*,'coucou calcul detr 460: ig, l',ig, l
       !
       !---------------------------------------------------------------------------
       !initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max 
       !---------------------------------------------------------------------------

       nbpb=0
       DO ig=1,ngrid
          IF (zw2(ig,l+1)>0.0_r8 .AND. zw2(ig,l+1).LT.1.e-10_r8) THEN
             !               stop'On tombe sur le cas particulier de thermcell_dry'
             !               print*,'On tombe sur le cas particulier de thermcell_plume'
             nbpb=nbpb+1
             zw2(ig,l+1)=0.0_r8
             linter(ig)=l+1
          ENDIF

          IF (zw2(ig,l+1).LT.0.0_r8) THEN 
             linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
                  &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
             zw2(ig,l+1)=0.0_r8
          ENDIF

          wa_moy(ig,l+1)=SQRT(zw2(ig,l+1)) 

          IF (wa_moy(ig,l+1).GT.wmaxa(ig)) THEN
             !   lmix est le niveau de la couche ou w (wa_moy) est maximum
             !on rajoute le calcul de lmix_bis
             IF (zqla(ig,l).LT.1.e-10_r8) THEN
                lmix_bis(ig)=l+1
             ENDIF
             lmix(ig)=l+1
             wmaxa(ig)=wa_moy(ig,l+1)
          ENDIF
       ENDDO

       IF (nbpb>0) THEN
          PRINT*,'WARNING on tombe ',nbpb,' x sur un pb pour l=',l,' dans thermcell_plume'
       ENDIF

       !=========================================================================
       ! FIN DE LA BOUCLE VERTICALE
    ENDDO
    !=========================================================================

    !on recalcule alim_star_tot
    DO ig=1,ngrid
       alim_star_tot(ig)=0.0_r8
    ENDDO
    DO ig=1,ngrid
       DO l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
       ENDDO
    ENDDO


    IF (prt_level.GE.20) PRINT*,'coucou calcul detr 470: ig, l', ig, l


    RETURN 
  END  SUBROUTINE thermcell_plume


  SUBROUTINE thermcell_qsat( &
       ngrid     ,&! INTEGER      , INTENT(IN   ) :: ngrid
       active   ,&! LOGICAL      , INTENT(IN   ) :: active(ngrid)
       pplev    ,&! REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid)
       ztemp    ,&! REAL(KIND=r8), INTENT(IN   ) :: ztemp(ngrid)
       zqta     ,&! REAL(KIND=r8), INTENT(IN   ) :: zqta(ngrid)
       zqsat    ) ! REAL(KIND=r8), INTENT(OUT  ) :: zqsat(ngrid)
    IMPLICIT NONE

    !#include "YOMCST.h"
    !#include "YOETHF.h"
    !#include "FCTTRE.h"


    !====================================================================
    ! DECLARATIONS
    !====================================================================

    ! Arguments
    INTEGER      , INTENT(IN   ) :: ngrid
    LOGICAL      , INTENT(IN   ) :: active(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: pplev(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: ztemp(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: zqta(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: zqsat(ngrid)

    !REAL(KIND=r8)    :: zpspsk(ngrid)

    ! Variables locales
    INTEGER :: ig
    INTEGER :: iter
    REAL(KIND=r8) :: Tbef(ngrid)
    REAL(KIND=r8) :: DT(ngrid)
!    REAL(KIND=r8) :: tdelta
    REAL(KIND=r8) :: qsatbef
    REAL(KIND=r8) :: zcor
    REAL(KIND=r8) :: qlbef
    REAL(KIND=r8) :: zdelta
    REAL(KIND=r8) :: zcvm5
!    REAL(KIND=r8) :: dqsat
    REAL(KIND=r8) :: num
    REAL(KIND=r8) :: denom
    REAL(KIND=r8) :: dqsat_dT
!    LOGICAL :: Zsat
    REAL(KIND=r8)    :: RLvCp
    !REAL(KIND=r8), SAVE :: DDT0=.01
    LOGICAL :: afaire(ngrid)
    LOGICAL :: tout_converge

    !====================================================================
    ! INITIALISATIONS
    !====================================================================
 zqsat=0.0_r8;
 Tbef=0.0_r8;
 DT=0.0_r8;
 qsatbef=0.0_r8;
 zcor=0.0_r8;
 qlbef=0.0_r8;
 zdelta=0.0_r8;
 zcvm5=0.0_r8;
 num=0.0_r8;
 denom=0.0_r8;
 dqsat_dT=0.0_r8;
 RLvCp=0.0_r8;



    RLvCp = RLVTT/RCPD
    tout_converge=.FALSE.
    afaire(:)=.FALSE.
    DT(:)=0.0_r8


    !====================================================================
    ! Routine a vectoriser en copiant active dans converge et en mettant
    ! la boucle sur les iterations a l'exterieur est en mettant
    ! converge= false des que la convergence est atteinte.
    !====================================================================

    DO ig=1,ngrid
       IF (active(ig)) THEN
          Tbef(ig)=ztemp(ig)
          zdelta=MAX(0.0_r8,SIGN(1.0_r8,RTT-Tbef(ig)))
          qsatbef= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig)
          qsatbef=MIN(0.5_r8,qsatbef)
          zcor=1.0_r8/(1.0_r8-retv*qsatbef)
          qsatbef=qsatbef*zcor
          qlbef=MAX(0.0_r8,zqta(ig)-qsatbef)
          DT(ig) = 0.5_r8*RLvCp*qlbef
          zqsat(ig)=qsatbef
       ENDIF
    ENDDO

    ! Traitement du cas ou il y a condensation mais faible
    ! On ne condense pas mais on dit que le qsat est le qta
    DO ig=1,ngrid
       IF (active(ig)) THEN
          IF (0.0_r8<ABS(DT(ig)).AND.ABS(DT(ig))<=DDT0) THEN
             zqsat(ig)=zqta(ig)
          ENDIF
       ENDIF
    ENDDO

    DO iter=1,10
       DO ig=1,ngrid
          afaire(ig)=ABS(DT(ig)).GT.DDT0
       ENDDO
       DO ig=1,ngrid
          IF (afaire(ig)) THEN
             Tbef(ig)=Tbef(ig)+DT(ig)
             zdelta=MAX(0.0_r8,SIGN(1.0_r8,RTT-Tbef(ig)))
             qsatbef= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig)
             qsatbef=MIN(0.5_r8,qsatbef)
             zcor=1.0_r8/(1.0_r8-retv*qsatbef)
             qsatbef=qsatbef*zcor
             qlbef=zqta(ig)-qsatbef
             zdelta=MAX(0.0_r8,SIGN(1.0_r8,RTT-Tbef(ig)))
             zcvm5=R5LES*(1.0_r8-zdelta) + R5IES*zdelta
             zcor=1.0_r8/(1.0_r8-retv*qsatbef)
             dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef,zcor)
             num=-Tbef(ig)+ztemp(ig)+RLvCp*qlbef
             denom=1.0_r8+RLvCp*dqsat_dT
             zqsat(ig) = qsatbef
             DT(ig)=num/denom
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE thermcell_qsat

  !
  ! $Id: thermcell_dry.F90 1403 2010-07-01 09:02:53Z fairhead $
  !
  SUBROUTINE thermcell_dry( &
       ngrid     , &! INTEGER , INTENT(IN   ) :: ngrid
       nlay      , &! INTEGER , INTENT(IN   ) :: nlay
       zlev      , &! REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
       pphi      , &! REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
       ztv       , &! REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
       alim_star , &! REAL(KIND=r8), INTENT(IN   ) :: alim_star(ngrid,nlay)
       lalim     , &! INTEGER , INTENT(IN   ) :: lalim(ngrid)
       lmin      , &! INTEGER , INTENT(IN   ) :: lmin(ngrid)
       zmax      , &! REAL(KIND=r8), INTENT(OUT  ) :: zmax(ngrid)
       wmax        )! REAL(KIND=r8), INTENT(OUT  ) :: wmax(ngrid)
!       lev_out     )! INTEGER , INTENT(IN   ) :: lev_out         ! niveau pour les print

    !--------------------------------------------------------------------------
    !thermcell_dry: calcul de zmax et wmax du thermique sec
    ! Calcul de la vitesse maximum et de la hauteur maximum pour un panache
    ! ascendant avec une fonction d'alimentation alim_star et sans changement 
    ! de phase.
    ! Le calcul pourrait etre sans doute simplifier.
    ! La temperature potentielle virtuelle dans la panache ascendant est
    ! la temperature potentielle virtuelle pondÃ©rÃ©e par alim_star.
    !--------------------------------------------------------------------------

    IMPLICIT NONE
    !#include "YOMCST.h"       
    !#include "iniprint.h"
    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay+1)
    REAL(KIND=r8), INTENT(IN   ) :: pphi(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: ztv(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: alim_star(ngrid,nlay)
    INTEGER      , INTENT(IN   ) :: lalim(ngrid)
    INTEGER      , INTENT(IN   ) :: lmin(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: zmax(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: wmax(ngrid)
!    INTEGER      , INTENT(IN   ) :: lev_out                           ! niveau pour les print

    !variables locales
    REAL(KIND=r8) :: zw2(ngrid,nlay+1)
    REAL(KIND=r8) :: f_star(ngrid,nlay+1)
    REAL(KIND=r8) :: ztva(ngrid,nlay+1)
    REAL(KIND=r8) :: wmaxa(ngrid)
    REAL(KIND=r8) :: wa_moy(ngrid,nlay+1)
    REAL(KIND=r8) :: linter(ngrid)
    REAL(KIND=r8) :: zlevinter(ngrid)
    INTEGER :: lmix(ngrid)
    INTEGER ::lmax(ngrid)
!    CHARACTER (LEN=20) :: modname='thermcell_dry'
!    CHARACTER (LEN=80) :: abort_message
    INTEGER :: l
    INTEGER :: ig


    !initialisations
 zmax=0.0_r8
 wmax=0.0_r8
 zw2=0.0_r8
 f_star=0.0_r8
 ztva=0.0_r8
 wmaxa=0.0_r8
 wa_moy=0.0_r8
 linter=0.0_r8
 zlevinter=0.0_r8
 lmix=0
 lmax=0
    !initialisations
    DO ig=1,ngrid
       DO l=1,nlay+1
          zw2(ig,l)=0.0_r8
          wa_moy(ig,l)=0.0_r8
       ENDDO
    ENDDO
    DO ig=1,ngrid
       DO l=1,nlay
          ztva(ig,l)=ztv(ig,l)
       ENDDO
    ENDDO
    DO ig=1,ngrid
       wmax(ig)=0.0_r8
       wmaxa(ig)=0.0_r8
    ENDDO
    !calcul de la vitesse a partir de la CAPE en melangeant thetav


    ! Calcul des F^*, integrale verticale de E^*
    DO ig=1,ngrid
       f_star(ig,1)=0.0_r8
    ENDDO

    DO l=1,nlay
       DO ig=1,ngrid
          f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)
       END DO
    ENDDO

    ! niveau (reel) auquel zw2 s'annule FH :n'etait pas initialise
    linter(:)=0.0_r8

    ! couche la plus haute concernee par le thermique. 
    lmax(:)=1

    ! Le niveau linter est une variable continue qui se trouve dans la couche
    ! lmax

    DO l=1,nlay-2
       DO ig=1,ngrid
          IF (l.EQ.lmin(ig).AND.lalim(ig).GT.1) THEN

             !------------------------------------------------------------------------
             !  Calcul de la vitesse en haut de la premiere couche instable.
             !  Premiere couche du panache thermique
             !------------------------------------------------------------------------

             zw2(ig,l+1)=2.0_r8*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)  &
                  &                     *(zlev(ig,l+1)-zlev(ig,l))  &
                  &                     *0.4_r8*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))

             !------------------------------------------------------------------------
             ! Tant que la vitesse en bas de la couche et la somme du flux de masse
             ! et de l'entrainement (c'est a dire le flux de masse en haut) sont
             ! positifs, on calcul
             ! 1. le flux de masse en haut  f_star(ig,l+1)
             ! 2. la temperature potentielle virtuelle dans la couche ztva(ig,l)
             ! 3. la vitesse au carré en haut zw2(ig,l+1)
             !------------------------------------------------------------------------

          ELSE IF (zw2(ig,l).GE.1e-10_r8) THEN

             ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+alim_star(ig,l)  &
                  &                    *ztv(ig,l))/f_star(ig,l+1)
             zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+  &
                  &                     2.0_r8*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
                  &                     *(zlev(ig,l+1)-zlev(ig,l))
          ENDIF
          ! determination de zmax continu par interpolation lineaire
          !------------------------------------------------------------------------

          IF (zw2(ig,l+1)>0.0_r8 .AND. zw2(ig,l+1).LT.1.e-10_r8) THEN
             !               stop'On tombe sur le cas particulier de thermcell_dry'
             !               print*,'On tombe sur le cas particulier de thermcell_dry'
             zw2(ig,l+1)=0.0_r8
             linter(ig)=l+1
             lmax(ig)=l
          ENDIF

          IF (zw2(ig,l+1).LT.0.0_r8) THEN
             linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
                  &           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
             zw2(ig,l+1)=0.0_r8
             lmax(ig)=l
          ENDIF

          wa_moy(ig,l+1)=SQRT(zw2(ig,l+1))

          IF (wa_moy(ig,l+1).GT.wmaxa(ig)) THEN
             !   lmix est le niveau de la couche ou w (wa_moy) est maximum
             lmix(ig)=l+1
             wmaxa(ig)=wa_moy(ig,l+1)
          ENDIF
       ENDDO
    ENDDO
    !IF (prt_level.GE.1) PRINT*,'fin calcul zw2'
    !
    ! Determination de zw2 max
    DO ig=1,ngrid
       wmax(ig)=0.0_r8
    ENDDO

    DO l=1,nlay
       DO ig=1,ngrid
          IF (l.LE.lmax(ig)) THEN
             zw2(ig,l)=SQRT(zw2(ig,l))
             wmax(ig)=MAX(wmax(ig),zw2(ig,l))
          ELSE
             zw2(ig,l)=0.0_r8
          ENDIF
       ENDDO
    ENDDO

    !   Longueur caracteristique correspondant a la hauteur des thermiques.
    DO  ig=1,ngrid
       zmax(ig)=0.0_r8
       zlevinter(ig)=zlev(ig,1)
    ENDDO
    DO  ig=1,ngrid
       ! calcul de zlevinter
       zlevinter(ig)=zlev(ig,lmax(ig)) + &
            &    (linter(ig)-lmax(ig))*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
       zmax(ig)=MAX(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
    ENDDO

    RETURN
  END SUBROUTINE thermcell_dry

  !
  ! $Header$
  !
  SUBROUTINE thermcell_closure( &
       ngrid      ,&!  INTEGER      , INTENT(IN) :: ngrid
       nlay       ,&!  INTEGER      , INTENT(IN) :: nlay
       r_aspect   ,&!  REAL(KIND=r8), INTENT(IN) :: r_aspect
       rho        ,&!  REAL(KIND=r8), INTENT(IN) :: rho(ngrid,nlay)
       zlev       ,&!  REAL(KIND=r8), INTENT(IN) :: zlev(ngrid,nlay)
       lalim      ,&!  INTEGER      , INTENT(IN) :: lalim(ngrid)
       alim_star  ,&!  REAL(KIND=r8), INTENT(IN) :: alim_star(ngrid,nlay)
       zmax       ,&!  REAL(KIND=r8), INTENT(IN) :: zmax(ngrid)
       wmax       ,&!  REAL(KIND=r8), INTENT(IN) :: wmax(ngrid)
       f          )!  REAL(KIND=r8), INTENT(OUT) :: f(ngrid)

    !-------------------------------------------------------------------------
    !thermcell_closure: fermeture, determination de f
    !
    ! Modification 7 septembre 2009
    ! 1. On enleve alim_star_tot des arguments pour le recalculer et etre ainis
    ! coherent avec l'integrale au numerateur.
    ! 2. On ne garde qu'une version des couples wmax,zmax et wmax_sec,zmax_sec
    ! l'idee etant que le choix se fasse a l'appel de thermcell_closure
    ! 3. Vectorisation en mettant les boucles en l l'exterieur avec des if
    !-------------------------------------------------------------------------
    IMPLICIT NONE

    !!#include "iniprint.h"
    !!#include "thermcell.h"
    INTEGER      , INTENT(IN   ) :: ngrid
    INTEGER      , INTENT(IN   ) :: nlay
    REAL(KIND=r8), INTENT(IN   ) :: r_aspect
    REAL(KIND=r8), INTENT(IN   ) :: rho(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zlev(ngrid,nlay)
    INTEGER      , INTENT(IN   ) :: lalim(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: alim_star(ngrid,nlay)
    REAL(KIND=r8), INTENT(IN   ) :: zmax(ngrid)
    REAL(KIND=r8), INTENT(IN   ) :: wmax(ngrid)
    REAL(KIND=r8), INTENT(OUT  ) :: f(ngrid)

    INTEGER :: ig
    INTEGER :: k       
!    REAL(KIND=r8)    ::  zdenom(ngrid)
    REAL(KIND=r8)    ::  alim_star2(ngrid)
    REAL(KIND=r8)    ::  alim_star_tot(ngrid)
    INTEGER :: llmax







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*,'THERMCELL CLOSURE 26E'

    alim_star2(:)=0.0_r8
    alim_star_tot(:)=0.0_r8
    f(:)=0.0_r8

    ! Indice vertical max (max de lalim) atteint par les thermiques sur le domaine
    llmax=1
    DO ig=1,ngrid
       IF (lalim(ig)>llmax) llmax=lalim(ig)
    ENDDO


    ! Calcul des integrales sur la verticale de alim_star et de
    !   alim_star^2/(rho dz)
    DO k=1,llmax-1
       DO ig=1,ngrid
          IF (k<lalim(ig)) THEN
             alim_star2(ig)=alim_star2(ig)+alim_star(ig,k)**2  &
                  /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,k)
          ENDIF
       ENDDO
    ENDDO


    DO ig=1,ngrid
       IF (alim_star2(ig)>1.e-10_r8) THEN
          f(ig)=wmax(ig)*alim_star_tot(ig)/  &
               (MAX(500.0_r8,zmax(ig))*r_aspect*alim_star2(ig))
       ENDIF
    ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TESTS POUR UNE NOUVELLE FERMETURE DANS LAQUELLE ALIM_STAR NE SERAIT
    ! PAS NORMALISE
    !           f(ig)=f(ig)*f_star(ig,2)/(f_star(ig,lalim(ig)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    RETURN
  END SUBROUTINE thermcell_closure


  SUBROUTINE thermcell_condens( &
       ngrid    ,&!    INTEGER      , INTENT(IN   )    :: ngrid
       active  ,&!    LOGICAL      , INTENT(IN   )    :: active(ngrid)
       zpspsk  ,&!    REAL(KIND=r8), INTENT(IN   )    :: zpspsk(ngrid)
       pplev   ,&!    REAL(KIND=r8), INTENT(IN   )    :: pplev (ngrid)
       ztla    ,&!    REAL(KIND=r8), INTENT(IN   )    :: ztla  (ngrid)
       zqta    ,&!    REAL(KIND=r8), INTENT(IN   )    :: zqta  (ngrid)
       zqla     )!    REAL(KIND=r8), INTENT(INOUT)    :: zqla  (ngrid)
    IMPLICIT NONE

    !#include "YOMCST.h"
    !#include "YOETHF.h"
    !#include "FCTTRE.h"


    !====================================================================
    ! DECLARATIONS
    !====================================================================

    ! Arguments
    INTEGER      , INTENT(IN   )    :: ngrid
    LOGICAL      , INTENT(IN   )    :: active(ngrid)
    REAL(KIND=r8), INTENT(IN   )    :: zpspsk(ngrid)
    REAL(KIND=r8), INTENT(IN   )    :: pplev (ngrid)
    REAL(KIND=r8), INTENT(IN   )    :: ztla  (ngrid)
    REAL(KIND=r8), INTENT(IN   )    :: zqta  (ngrid)
    REAL(KIND=r8), INTENT(INOUT)    :: zqla  (ngrid)

    ! Variables locales
    INTEGER :: ig
    INTEGER :: iter
    REAL(KIND=r8)    :: Tbef(ngrid)
    REAL(KIND=r8)    :: DT(ngrid)
!    REAL(KIND=r8)    :: tdelta! Unused variable 'dqsat' declared at (1)
    REAL(KIND=r8)    :: qsatbef
    REAL(KIND=r8)    :: zcor
    REAL(KIND=r8)    :: qlbef
    REAL(KIND=r8)    :: zdelta
    REAL(KIND=r8)    :: zcvm5
 !   REAL(KIND=r8)    :: dqsat ! Unused variable 'dqsat' declared at (1)
    REAL(KIND=r8)    :: num
    REAL(KIND=r8)    :: denom
    REAL(KIND=r8)    :: dqsat_dT
!    LOGICAL :: Zsat
    REAL(KIND=r8)    :: RLvCp
    !REAL(KIND=r8), SAVE :: DDT0=.01
    LOGICAL :: afaire(ngrid)
    LOGICAL :: tout_converge

    !====================================================================
    ! INITIALISATIONS
    !====================================================================
     Tbef=0.0_r8;
     DT=0.0_r8;
     qsatbef=0.0_r8;
     zcor=0.0_r8;
     qlbef=0.0_r8;
     zdelta=0.0_r8;
     zcvm5=0.0_r8;
     num=0.0_r8;
     denom=0.0_r8;
     dqsat_dT=0.0_r8;
     RLvCp=0.0_r8;
    
    
    RLvCp = RLVTT/RCPD
    tout_converge=.FALSE.
    afaire(:)=.FALSE.
    DT(:)=0.0_r8


    !====================================================================
    ! Routine a vectoriser en copiant active dans converge et en mettant
    ! la boucle sur les iterations a l'exterieur est en mettant
    ! converge= false des que la convergence est atteinte.
    !====================================================================

    DO ig=1,ngrid
       IF (active(ig)) THEN
          Tbef(ig)=ztla(ig)*zpspsk(ig)
          zdelta=MAX(0.0_r8,SIGN(1.0_r8,RTT-Tbef(ig)))
          qsatbef= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig)
          qsatbef=MIN(0.5_r8,qsatbef)
          zcor=1.0_r8/(1.0_r8-retv*qsatbef)
          qsatbef=qsatbef*zcor
          qlbef=MAX(0.0_r8,zqta(ig)-qsatbef)
          DT(ig) = 0.5_r8*RLvCp*qlbef
       ENDIF
    ENDDO

    DO iter=1,10
       afaire(:)=ABS(DT(:)).GT.DDT0
       DO ig=1,ngrid
          IF (afaire(ig)) THEN
             Tbef(ig)=Tbef(ig)+DT(ig)
             zdelta=MAX(0.0_r8,SIGN(1.0_r8,RTT-Tbef(ig)))
             qsatbef= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig)
             qsatbef=MIN(0.5_r8,qsatbef)
             zcor=1.0_r8/(1.0_r8-retv*qsatbef)
             qsatbef=qsatbef*zcor
             qlbef=zqta(ig)-qsatbef
             zdelta=MAX(0.0_r8,SIGN(1.0_r8,RTT-Tbef(ig)))
             zcvm5=R5LES*(1.0_r8-zdelta) + R5IES*zdelta
             zcor=1.0_r8/(1.0_r8-retv*qsatbef)
             dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef,zcor)
             num=-Tbef(ig)+ztla(ig)*zpspsk(ig)+RLvCp*qlbef
             denom=1.0_r8+RLvCp*dqsat_dT
             zqla(ig) = MAX(0.0_r8,zqta(ig)-qsatbef) 
             DT(ig)=num/denom
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE thermcell_condens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! $Header$
  !
  !
  !  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
  !                 veillez  n'utiliser que des ! pour les commentaires
  !                 et  bien positionner les & des lignes de continuation
  !                 (les placer en colonne 6 et en colonne 73)
  !
  !     ------------------------------------------------------------------
  !     This COMDECK includes the Thermodynamical functions for the cy39
  !       ECMWF Physics package.
  !       Consistent with YOMCST Basic physics constants, assuming the
  !       partial pressure of water vapour is given by a first order
  !       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
  !       in YOETHF
  !     ------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8) FUNCTION FOEEW( PTARG,PDELARG )
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: PTARG
    REAL(KIND=r8), INTENT(IN   ) :: PDELARG

    !    fonction psat(T)

    FOEEW  = EXP (                               &
         (R3LES*(1.0_r8-PDELARG)+R3IES*PDELARG) * (PTARG-RTT)     &
         / (PTARG-(R4LES*(1.0_r8-PDELARG)+R4IES*PDELARG)) )

  END FUNCTION FOEEW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8) FUNCTION FOEDE ( PTARG,PDELARG,P5ARG,PQSARG,PCOARG )
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: PTARG
    REAL(KIND=r8), INTENT(IN   ) :: PDELARG
    REAL(KIND=r8), INTENT(IN   ) :: P5ARG
    REAL(KIND=r8), INTENT(IN   ) :: PQSARG
    REAL(KIND=r8), INTENT(IN   ) :: PCOARG 

    FOEDE  = PQSARG*PCOARG*P5ARG                     &
         & / (PTARG-(R4LES*(1.0_r8-PDELARG)+R4IES*PDELARG))**2

  END FUNCTION FOEDE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8) FUNCTION qsats(ptarg)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: ptarg

    qsats = 100.0_r8 * 0.622_r8 * 10.0_r8                               &
         &           ** (2.07023_r8 - 0.00320991_r8 * ptarg                       &
         &           - 2484.896_r8 / ptarg + 3.56654_r8 * LOG10(ptarg))
  END FUNCTION qsats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8) FUNCTION qsatl(ptarg)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) ::  ptarg

    qsatl = 100.0_r8 * 0.622_r8 * 10.0_r8                                      &
         &           ** (23.8319_r8 - 2948.964_r8 / ptarg                         &
         &           - 5.028_r8 * LOG10(ptarg)                                 &
         &           - 29810.16_r8 * EXP( - 0.0699382_r8 * ptarg)                 &
         &           + 25.21935_r8 * EXP( - 2999.924_r8 / ptarg))

  END FUNCTION qsatl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8) FUNCTION dqsats(ptarg,pqsarg)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) ::  ptarg
    REAL(KIND=r8), INTENT(IN   ) ::  pqsarg

    dqsats = RLVTT/RCPD*pqsarg * (3.56654_r8/ptarg         &
         &                     +2484.896_r8*LOG(10.0_r8)/ptarg**2                  &
         &                     -0.00320991_r8*LOG(10.0_r8))

  END FUNCTION dqsats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8) FUNCTION dqsatl(ptarg,pqsarg)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) ::  ptarg
    REAL(KIND=r8), INTENT(IN   ) ::  pqsarg

    dqsatl = RLVTT/RCPD*pqsarg*LOG(10.0_r8)*                &
         &                (2948.964_r8/ptarg**2-5.028_r8/LOG(10.0_r8)/ptarg           &
         &                +25.21935_r8*2999.924_r8/ptarg**2*EXP(-2999.924_r8/ptarg)  &
         &                +29810.16_r8*0.0699382_r8*EXP(-0.0699382_r8*ptarg))
  END FUNCTION dqsatl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=r8) FUNCTION fsta(x)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: x
    fsta = 1.0_r8 / (1.0_r8+10.0_r8*x*(1+8.0_r8*x))
  END FUNCTION fsta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=r8) FUNCTION fins(x)
    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN   ) :: x

    fins = SQRT(1.0_r8-18.0_r8*x)
  END FUNCTION fins

  !
  !****************************************************************************************
  !
END MODULE ThermalCell


!PROGRAM Main
!  USE ThermalCell
!
!END PROGRAM Main
