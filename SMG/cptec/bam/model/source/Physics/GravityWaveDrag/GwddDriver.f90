!
!  $Author: pkubota $
!  $Date: 2008/08/19 16:57:04 $
!  $Revision: 1.3 $
!
MODULE GwddDriver

  USE Constants, ONLY :     &
       cp,            &
       grav,          &
       gasr,          &
       r8,           &
       r4,i8

  USE GwddSchemeAlpert, ONLY: &
      InitGwddSchAlpert  ,&
      GwddSchAlpert

  USE Options, ONLY: &
       reducedGrid,igwd,nfvar,fNameTopo,fNameOrgvar,fNameGtopog,fNameHPRIME,nfprt

  USE GwddSchemeCAM, ONLY: &
      gw_inti  ,&
      gw_intr

  USE Gwdd_ECMWF, ONLY: &
      InitGwdd_ECMWF  ,&
      Run_Gwdd_ECMWF

  USE GwddSchemeUSSP, ONLY: &
      Init_GwUSSP, &
      gw_ussp

  USE GwddSchemeCPTEC, ONLY: &
      Init_Gwave,&
      g_wave

  USE Diagnostics, ONLY: &
        updia,dodia , &
        StartStorDiag,&
        nDiag_txgwds, & ! gravity wave drag surface zonal stress
        nDiag_tygwds, & ! gravity wave drag surface meridional stress
        nDiag_gwduzc, & ! gravity wave drag zonal momentum change
        nDiag_gwdvmc, & ! gravity wave drag meridional momentum change
        nDiag_txgwdp, & ! gravity wave drag profile zonal stress 
        nDiag_tygwdp    ! gravity wave drag profile meridional stress

  USE Options, ONLY :       &
       igwd            

  USE Utils, ONLY: &
        IJtoIBJB, &
        LinearIJtoIBJB

  USE GridHistory, ONLY:       &
       IsGridHistoryOn, StoreGridHistory, StoreMaskedGridHistory, dogrh

  USE FieldsPhysics, ONLY: &
      var,topoi

  USE IOLowLevel, ONLY: &
       ReadVar      

    IMPLICIT NONE
  SAVE


  PRIVATE

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   G_dhdx !Zonal Gradient of Topography (m/m)
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   G_dhdy !Meridional Gradient of Topography (m/m)

  PUBLIC :: InitGWDDDriver
  PUBLIC :: Gwdd_Driver

CONTAINS
  SUBROUTINE InitGWDDDriver(ibMax,jbMax,iMax,jMax,kmax, ibMaxPerJB)
   IMPLICIT NONE
    INTEGER      , INTENT(IN) :: ibMax
    INTEGER      , INTENT(IN) :: jbMax
    INTEGER      , INTENT(IN) :: iMax
    INTEGER      , INTENT(IN) :: jMax
    INTEGER      , INTENT(IN) :: kmax

    INTEGER      , INTENT(IN) :: ibMaxPerJB(jbMax)

    ALLOCATE( G_dhdx (ibMax,jbMax))
    G_dhdx=0.0_r8
    ALLOCATE( G_dhdy (ibMax,jbMax))
    G_dhdy=0.0_r8

    CALL InitVariancia(iMax,jMax,igwd,nfvar,fNameOrgvar,fNameGtopog,fNameTopo)

    IF(TRIM(igwd) =='YES')CALL InitGwddSchAlpert()
    IF(TRIM(igwd) =='CAM')CALL gw_inti     (kMax)
    IF(TRIM(igwd) =='USS')CALL Init_GwUSSP ()
    IF(TRIM(igwd) =='USS')CALL Init_Gwave  ()
    IF(TRIM(igwd) =='GMB')THEN
      CALL InitGwddSchAlpert()
      CALL InitGwdd_ECMWF(iMax,jMax,kMax,ibMax,jbMax,nfprt,&
                                              nfvar,reducedGrid,fNameHPRIME,ibMaxPerJB)
      CALL Init_GwUSSP ()
      CALL Init_Gwave  ()
    END IF 
  END SUBROUTINE InitGWDDDriver

 SUBROUTINE InitVariancia(iMax,jMax,igwd,nfvar,fNameOrgvar,fNameGtopog,fNameTopo)
    INTEGER          , INTENT(IN   ) :: iMax
    INTEGER          , INTENT(IN   ) :: jMax
    CHARACTER(LEN=*) , INTENT(in   ) :: igwd
    INTEGER          , INTENT(in   ) :: nfvar
    CHARACTER(LEN=*) , INTENT(in   ) :: fNameOrgvar    
    CHARACTER(LEN=*) , INTENT(in   ) :: fNameGtopog
    CHARACTER(LEN=*) , INTENT(in   ) :: fNameTopo
    
    INTEGER       :: LRecIn,irec,ierr
    REAL(KIND=r8) ::   var_in (iMax,jMax)
    REAL(KIND=r4) ::   buffer (iMax,jMax)
    !
    ! VARIANCY TOPOGRAPHY
    !
    var_in  =0.0_r8
    buffer  =0.0_r4
    INQUIRE (IOLENGTH=LRecIn) buffer 
    OPEN (UNIT=nfvar,FILE=TRIM(fNameOrgvar),FORM='UNFORMATTED', ACCESS='DIRECT', &
         ACTION='READ',RECL=LRecIn,STATUS='OLD', IOSTAT=ierr)

    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameOrgvar), ierr
       STOP "**(ERROR)**"
    END IF


    IF(TRIM(igwd).EQ.'YES' .or. TRIM(igwd).EQ.'CAM' .or. TRIM(igwd).EQ.'USS' .or. TRIM(igwd) == 'GMB') THEN
       irec=1
       CALL ReadVar(nfvar,irec,var_in,0)
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(var_in,var)
       ELSE
          CALL IJtoIBJB(var_in,var)
       END IF
    END IF
    CLOSE(nfvar,STATUS='KEEP')
    !
    ! REAL TOPOGRAPHY
    !
    INQUIRE (IOLENGTH=LRecIn) buffer 
    OPEN (UNIT=nfvar,FILE=TRIM(fNameTopo),FORM='UNFORMATTED', ACCESS='DIRECT', &
         ACTION='READ',RECL=LRecIn,STATUS='OLD', IOSTAT=ierr)

    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameTopo), ierr
       STOP "**(ERROR)**"
    END IF


    IF(TRIM(igwd).EQ.'YES' .or. TRIM(igwd).EQ.'CAM' .or. TRIM(igwd).EQ.'USS' .or. TRIM(igwd) == 'GMB') THEN
       irec=1
       CALL ReadVar(nfvar,irec,var_in,0)
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(var_in,topoi)
       ELSE
          CALL IJtoIBJB(var_in,topoi)
       END IF
    END IF
    topoi=MAX(topoi,0.0_r8)
    CLOSE(nfvar,STATUS='KEEP')

    IF( TRIM(igwd) == 'USS'.or. TRIM(igwd) == 'GMB') THEN

       var_in  =0.0_r8
       buffer  =0.0_r4
       INQUIRE (IOLENGTH=LRecIn) buffer
       OPEN (UNIT=nfvar,FILE=TRIM(fNameGtopog),FORM='UNFORMATTED', ACCESS='DIRECT', &
         ACTION='READ',RECL=LRecIn,STATUS='OLD', IOSTAT=ierr)

       IF (ierr /= 0) THEN
          WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameGtopog), ierr
          STOP "**(ERROR)**"
       END IF


       irec=1
       CALL ReadVar(nfvar,irec,var_in,1)
       !IF (reducedGrid) THEN
       !   CALL LinearIJtoIBJB(var_in,topo)
       !ELSE
       !   CALL IJtoIBJB(var_in,topo)
       !END IF

       irec=2
       CALL ReadVar(nfvar,irec,var_in,1)
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(var_in,G_dhdx)
       ELSE
          CALL IJtoIBJB(var_in,G_dhdx)
       END IF

       irec=3
       CALL ReadVar(nfvar,irec,var_in,0)
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(var_in,G_dhdy)
       ELSE
          CALL IJtoIBJB(var_in,G_dhdy)
       END IF
       CLOSE(nfvar,STATUS='KEEP')

    END IF

  END SUBROUTINE InitVariancia
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE Gwdd_Driver (      prsi ,prsl  ,phii ,phil    ,&
            gu   ,gv   ,gt   , gq,chug, chvg, chtg,xdrag, ydrag, &
            var, varcut,  ncols,&
            kmax,latco,dt,imask,colrad,topog,pblh,cu_hr,cu_kbot,cu_ktop,cu_Kuo ,dump)
   IMPLICIT NONE
    INTEGER      , INTENT(in   ) :: ncols
    INTEGER      , INTENT(in   ) :: kmax
    INTEGER      , INTENT(in   ) :: latco
    REAL(KIND=r8), INTENT(in   ) :: dt
    REAL(KIND=r8), INTENT(in   ) :: prsi   (ncols,kMax+1)  !     prsi     - real, pressure at layer interfaces [Pa]
    REAL(KIND=r8), INTENT(in   ) :: prsl   (ncols,kMax)    !     prsl     - real, mean layer presure [Pa]
    REAL(KIND=r8), INTENT(in   ) :: phii   (nCols,kMax+1) !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8), INTENT(in   ) :: phil   (nCols,kMax)   !===>  PHIL(K) INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    !REAL(KIND=r8), INTENT(inout) :: ps     (ncols)       !ln(ps) cb   timestep t-1
    REAL(KIND=r8), INTENT(in   ) :: gu     (ncols,kmax)
    REAL(KIND=r8), INTENT(in   ) :: gv     (ncols,kmax)
    REAL(KIND=r8), INTENT(in   ) :: gt     (ncols,kmax)
    REAL(KIND=r8), INTENT(in   ) :: gq     (ncols,kmax) 
    REAL(KIND=r8), INTENT(inout) :: chug   (ncols,kmax)
    REAL(KIND=r8), INTENT(inout) :: chvg   (ncols,kmax)
    REAL(KIND=r8), INTENT(inout) :: chtg   (ncols,kmax)
    REAL(KIND=r8), INTENT(inout) :: xdrag  (ncols)
    REAL(KIND=r8), INTENT(inout) :: ydrag  (ncols)
    REAL(KIND=r8), INTENT(inout) :: var    (ncols)
    REAL(KIND=r8), INTENT(in   ) :: varcut
    INTEGER(KIND=i8)      , INTENT(in   ) :: imask(ncols)
    REAL(KIND=r8), INTENT(in   ) :: colrad(ncols)
    REAL(KIND=r8), INTENT(in   ) :: topog(ncols)
    REAL(KIND=r8), INTENT(in   ) :: pblh (ncols)
    REAL(KIND=r8), INTENT(in   ) :: cu_hr   (ncols,kmax)
    INTEGER, INTENT(in   ) :: cu_kbot (ncols)
    INTEGER, INTENT(in   ) :: cu_ktop (ncols)
    INTEGER, INTENT(in   ) :: cu_Kuo  (ncols)
    REAL(KIND=r8), INTENT(inout) :: dump(ncols,kmax)
    REAL(KIND=r8) :: dhdx(ncols)
    REAL(KIND=r8) :: dhdy(ncols)
    REAL(KIND=r8) :: dhdx_dhdy(ncols)
    REAL(KIND=r8) :: lnps(ncols) 

    REAL(KIND=r8) ::psi(ncols) 
    REAL(KIND=r8) ::tmp(ncols,kmax)
    REAL(KIND=r8) ::ums(ncols,kmax)
    REAL(KIND=r8) ::sgh(ncols)
    REAL(KIND=r8) :: up   (ncols,kmax)
    REAL(KIND=r8) :: vp   (ncols,kmax)
    REAL(KIND=r8) :: um   (ncols,kmax)
    REAL(KIND=r8) :: vm   (ncols,kmax)
    REAL(KIND=r8) :: du_dt (ncols,kmax)
    REAL(KIND=r8) :: dv_dt (ncols,kmax)
    REAL(KIND=r8) :: du2_dt (ncols,kmax)
    REAL(KIND=r8) :: dv2_dt (ncols,kmax)

    REAL(KIND=r8) :: xtens (ncols,kmax+1)
    REAL(KIND=r8) :: ytens (ncols,kmax+1)
    REAL(KIND=r8) :: chug2   (ncols,kmax)
    REAL(KIND=r8) :: chvg2   (ncols,kmax)

    INTEGER :: i,k,iret
    REAL(KIND=r8) ::landfrac(ncols)
    REAL(KIND=r8) ::rlat(ncols)
    chug=0.0_r8
    chvg=0.0_r8
    chtg=0.0_r8
    du_dt=0.0_r8
    dv_dt=0.0_r8
    du2_dt=0.0_r8
    dv2_dt=0.0_r8
    chug2=0.0_r8
    chvg2=0.0_r8
    IF(TRIM(igwd) =='YES')THEN
       DO i=1,nCols
          var(i) =MIN(varcut,var(i))
          var(i) =MAX(0.0_r8,var(i))
       END DO
       CALL GwddSchAlpert( prsi ,prsl  ,phii ,phil    ,&
                         gu    ,gv   ,gt   , chug  , chvg , xdrag, ydrag,xtens,&
                         ytens, var,varcut, ncols, kmax)
    ELSE IF(TRIM(igwd) =='CAM')THEN
        !PRINT*,TRIM(igwd)
       DO k=1,kMax
          DO i=1,nCols
             ums(i,k)=MAX(1.0e-12_r8,gq(i,k))
             tmp(i,k)=gt(i,k)/(1.0e0_r8+0.608e0_r8*ums(i,k))
             up (i,k)=gu(i,k)/SIN( colrad(i))
             vp (i,k)=gv(i,k)/SIN( colrad(i))
          END DO
       END DO   
       DO i=1,nCols
          psi(i)=topog(i)*grav
          !ps(i) = (EXP(ps(i))*10.0_r8)*100.0_r8   !ln(ps) cb   timestep t-1 ----> Pa
          var(i) =MIN(varcut,var(i))
          var(i) =MAX(0.0_r8,var(i))
          sgh(i) = sqrt(var(i))
          rlat(i) =  colrad(i)  -   (3.1415926e0_r8/2.0_r8)
          IF(imask(i).GE.1_i8) THEN
             landfrac(i)=1.0_r8
          ELSE
             landfrac(i)=0.0_r8
          END IF
       END DO
    
       CALL gw_intr ( &
         ncols    , &!INTEGER, INTENT(IN   ) :: pcols
         kMax     , &!INTEGER, INTENT(IN   ) :: pver
         kMax+1   , &!INTEGER, INTENT(IN   ) :: pverp
         prsi ,prsl  ,phii ,phil    ,&
         tmp      , &!REAL(r8), INTENT(in) :: gt (pcols,pver)  
         ums      , &!REAL(r8), INTENT(in) :: gq (pcols,pver)  
         up       , &!REAL(r8), INTENT(in) :: gu (pcols,pver)  
         vp       , &!REAL(r8), INTENT(in) :: gv (pcols,pver)  
         sgh      , &!REAL(r8), INTENT(in) :: sgh(pcols)                ! standard deviation of orography
         2*dt     , &!REAL(r8), INTENT(in) :: dt                        ! time step
         landfrac , &!REAL(r8), INTENT(in) :: landfrac(pcols)        ! Land fraction
         rlat     , &!REAL(r8), INTENT(in) :: rlat(pcols)             ! latitude in radians for columns
         psi      , &
         chug     , &
         chvg     , &
         chtg       )
       DO k=1,kMax
          DO i=1,nCols
             !PRINT*,chug(i,k),chvg(i,k)
             chug (i,k)=chug(i,k)*SIN( colrad(i))
             chvg (i,k)=chvg(i,k)*SIN( colrad(i))
             chtg (i,k)=chtg(i,k)*(1.0e0_r8+0.608e0_r8*ums(i,k))
          END DO
       END DO   


    ELSE IF(TRIM(igwd) == 'USS')THEN
       !PRINT*,TRIM(igwd)
       DO k=1,kMax
          DO i=1,nCols
             ums(i,k)=MAX(1.0e-12_r8,gq(i,k))
             tmp(i,k)=gt(i,k)/(1.0e0_r8+0.608e0_r8*ums(i,k))
             
             up (i,k)=gu(i,k)/SIN(colrad(i))
             vp (i,k)=gv(i,k)/SIN( colrad(i))
             
             um (i,k)=gu(i,k)/SIN( colrad(i))
             vm (i,k)=gv(i,k)/SIN( colrad(i))
          END DO
       END DO   
       DO i=1,nCols
          dhdx(i)     =G_dhdx(i,latco)
          dhdy(i)     =G_dhdy(i,latco)
          dhdx_dhdy(i)=dhdx(i)*dhdy(i)
          psi(i)      =topog(i)*grav
          !ps(i)       =(EXP(ps(i))*10.0_r8)*100.0_r8
          var(i) =MAX(0.0_r8,var(i))
          sgh(i)      =sqrt(var(i))
          rlat(i)     =colrad(i)  -   (3.1415926e0_r8/2.0_r8)
          IF(imask(i).GE.1_i8) THEN
             landfrac(i)=1.0_r8
          ELSE
             landfrac(i)=0.0_r8
          END IF
       END DO
       !PRINT*,'g_wave'
       CALL g_wave( &
             kMax         , & !INTEGER      , INTENT(IN   ) :: kMax
             nCols        , & !INTEGER      , INTENT(IN   ) :: nCols ! number of points per row
             prsi ,prsl  ,phii ,phil    ,&
             tmp          , & !REAL(KIND=r8), INTENT(IN   ) :: gt               (nCols,kMax)    !REAL(r8), INTENT(in) :: gt (nCols,kMax) 
             ums          , & !REAL(KIND=r8), INTENT(IN   ) :: gq               (nCols,kMax)    !REAL(r8), INTENT(in) :: gq (nCols,kMax) 
             up           , & !REAL(KIND=r8), INTENT(IN   ) :: gu               (nCols,kMax)    !REAL(r8), INTENT(in) :: gu (nCols,kMax) 
             vp           , & !REAL(KIND=r8), INTENT(IN   ) :: gv               (nCols,kMax)    !REAL(r8), INTENT(in) :: gv (nCols,kMax) 
             topog        , & !REAL(KIND=r8), INTENT(IN   ) :: topo           (nCols)
             !colrad       , & !REAL(KIND=r8), INTENT(IN   ) :: colrad           (nCols)    
             sgh          , & !REAL(KIND=r8) ,INTENT(in   ) :: sd_orog     (nCols)! standard deviation of orography (m)
             dhdx*dhdx    , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_xx(nCols)! (dh/dx)^2 grid box average of the
             dhdy*dhdy    , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_yy(nCols)! (dh/dy)^2 grid box average of the
             dhdx_dhdy    , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_xy(nCols)! (dh/dx)*(dh/dy)
             imask        , & !INTEGER      , INTENT(in   ) :: imask            (nCols)! index for land points
             2*dt         , & !REAL(KIND=r8),INTENT(in    ) :: timestep            ! timestep (s)
             du_dt        , & !REAL(KIND=r8),INTENT(out   ) :: du_dt(nCols,kMax)   ! total GWD du/dt on land/theta
             dv_dt        , & !REAL(KIND=r8),INTENT(out   ) :: dv_dt(nCols,kMax)   ! total GWD dv/dt on land/theta
             iret           ) !INTEGER      ,INTENT(OUT   ) :: iret                   ! return code : iret=0 normal exit


        CALL gw_ussp( &
               kMax       , &
               nCols      , &
               prsi ,prsl  ,phii ,phil    ,&
               tmp        , &
               ums        , &
               up         , &
               vp         , &
               topog      , &
               rlat       , &
               chug       , &
               chvg         &
                            )
       DO k=1,kMax
          DO i=1,nCols
             up   (i,k)= up(i,k) +  2.0_r8*dt*chug(i,k)
             vp   (i,k)= vp(i,k) +  2.0_r8*dt*chvg(i,k)
          END DO
       END DO   

       DO k=1,kMax
          DO i=1,nCols
             chug (i,k)=((du_dt(i,k) + chug(i,k)))*SIN( colrad(i))
             chvg (i,k)=((dv_dt(i,k) + chvg(i,k)))*SIN( colrad(i))
          END DO
       END DO   

    ELSE IF(TRIM(igwd) == 'GMB')THEN

       !PRINT*,TRIM(igwd)
       DO k=1,kMax
          DO i=1,nCols
             ums(i,k)=MAX(1.0e-12_r8,gq(i,k))
             tmp(i,k)=gt(i,k)/(1.0e0_r8+0.608e0_r8*ums(i,k))
             
             up (i,k)=gu(i,k)/SIN(colrad(i))
             vp (i,k)=gv(i,k)/SIN( colrad(i))
             
             um (i,k)=gu(i,k)/SIN( colrad(i))
             vm (i,k)=gv(i,k)/SIN( colrad(i))
          END DO
       END DO   
       DO i=1,nCols
          dhdx(i)     =G_dhdx(i,latco)
          dhdy(i)     =G_dhdy(i,latco)
          dhdx_dhdy(i)=dhdx(i)*dhdy(i)
          psi(i)      =topog(i)*grav
          !ps(i)       =(EXP(ps(i))*10.0_r8)*100.0_r8
          var(i) =MAX(0.0_r8,var(i))
          sgh(i)      =sqrt(var(i))
          rlat(i)     =colrad(i)  -   (3.1415926e0_r8/2.0_r8)
          IF(imask(i).GE.1_i8) THEN
             landfrac(i)=1.0_r8
          ELSE
             landfrac(i)=0.0_r8
          END IF
       END DO
       CALL g_wave( &
             kMax         , & !INTEGER      , INTENT(IN   ) :: kMax
             nCols        , & !INTEGER      , INTENT(IN   ) :: nCols ! number of points per row
             prsi ,prsl  ,phii ,phil    ,&
             tmp          , & !REAL(KIND=r8), INTENT(IN   ) :: gt               (nCols,kMax)    !REAL(r8), INTENT(in) :: gt (nCols,kMax) 
             ums          , & !REAL(KIND=r8), INTENT(IN   ) :: gq               (nCols,kMax)    !REAL(r8), INTENT(in) :: gq (nCols,kMax) 
             up           , & !REAL(KIND=r8), INTENT(IN   ) :: gu               (nCols,kMax)    !REAL(r8), INTENT(in) :: gu (nCols,kMax) 
             vp           , & !REAL(KIND=r8), INTENT(IN   ) :: gv               (nCols,kMax)    !REAL(r8), INTENT(in) :: gv (nCols,kMax) 
             topog        , & !REAL(KIND=r8), INTENT(IN   ) :: topo           (nCols)
             !colrad       , & !REAL(KIND=r8), INTENT(IN   ) :: colrad           (nCols)    
             sgh          , & !REAL(KIND=r8) ,INTENT(in   ) :: sd_orog     (nCols)! standard deviation of orography (m)
             dhdx*dhdx    , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_xx(nCols)! (dh/dx)^2 grid box average of the
             dhdy*dhdy    , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_yy(nCols)! (dh/dy)^2 grid box average of the
             dhdx_dhdy    , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_xy(nCols)! (dh/dx)*(dh/dy)
             imask        , & !INTEGER      , INTENT(in   ) :: imask            (nCols)! index for land points
             2*dt         , & !REAL(KIND=r8),INTENT(in    ) :: timestep            ! timestep (s)
             du2_dt        , & !REAL(KIND=r8),INTENT(out   ) :: du_dt(nCols,kMax)   ! total GWD du/dt on land/theta
             dv2_dt        , & !REAL(KIND=r8),INTENT(out   ) :: dv_dt(nCols,kMax)   ! total GWD dv/dt on land/theta
             iret           ) !INTEGER      ,INTENT(OUT   ) :: iret                   ! return code : iret=0 normal exit


        CALL gw_ussp( &
               kMax       , &
               nCols      , &
               prsi ,prsl  ,phii ,phil    ,&
               tmp        , &
               ums        , &
               up         , &
               vp         , &
               topog      , &
               rlat       , &
               chug2       , &
               chvg2         &
                            )
       DO k=1,kMax
          DO i=1,nCols
             chug2 (i,k)=((du2_dt(i,k) + chug2(i,k)))*SIN( colrad(i))
             chvg2 (i,k)=((dv2_dt(i,k) + chvg2(i,k)))*SIN( colrad(i))
          END DO
       END DO   

       DO i=1,nCols
          var(i) =MIN(varcut,var(i))
          var(i) =MAX(0.0_r8,var(i))
       END DO

       CALL GwddSchAlpert( prsi ,prsl  ,phii ,phil    ,&
                         gu    ,gv   ,gt   , chug  , chvg , xdrag, ydrag,xtens,&
                         ytens, var,varcut, ncols, kmax)

       DO k=1,kMax
          DO i=1,nCols
             ums(i,k)=MAX(1.0e-12_r8,gq(i,k))
             tmp(i,k)=gt(i,k)/(1.0e0_r8+0.608e0_r8*ums(i,k))
             
             up (i,k)=gu(i,k)/SIN(colrad(i))
             vp (i,k)=gv(i,k)/SIN( colrad(i))
             
          END DO
       END DO   


        CALL Run_Gwdd_ECMWF(      &
        nCols      , &
        kMax       , &
        latco      , &
        du_dt      , &
        dv_dt      , &
        prsi ,prsl  ,&
        up         , &
        vp         , &
        tmp        , &
        ums        , &
        pblh       , &
        cu_hr      , &
        cu_kbot    , &
        cu_ktop    , &
        cu_Kuo     , &
        colrad     , &
        2*dt        )

       DO k=1,kMax
          DO i=1,nCols
             du_dt (i,k)=du_dt(i,k)*SIN( colrad(i))
             dv_dt (i,k)=dv_dt(i,k)*SIN( colrad(i))
          END DO
       END DO
       DO k=1,kMax
          DO i=1,nCols
             IF(ABS(du_dt(i,k))>ABS(chug(i,k)))THEN 
                chug (i,k)=du_dt(i,k)
             ELSE
                chug (i,k)=chug(i,k)
             END IF
             IF(ABS(dv_dt(i,k))>ABS(chvg(i,k)))THEN
                chvg (i,k)=dv_dt(i,k)!/3
             ELSE
                chvg (i,k)=chvg(i,k)
             END IF
          END DO
       END DO



    END IF

    !-----------------
    ! Storage Diagnostic Fields
    !------------------
    IF( StartStorDiag)THEN
       CALL GwddDiagnStorage(latco,nCols,kMax,xdrag,ydrag,chug,chvg,xtens,ytens )
    END IF 
    !-----------------
    ! Storage GridHistory Fields
    !------------------
    IF(IsGridHistoryOn())THEN
       CALL GwddGridHistoryStorage()
    END IF

  END SUBROUTINE Gwdd_Driver
  
  SUBROUTINE GwddGridHistoryStorage()
   IMPLICIT NONE

   RETURN      
  END SUBROUTINE GwddGridHistoryStorage
  
  
  SUBROUTINE GwddDiagnStorage(latco,nCols,kMax,xdrag,ydrag,&
                              chug,chvg,xtens,ytens )
   IMPLICIT NONE
   INTEGER      , INTENT(IN   ) :: latco
   INTEGER      , INTENT(IN   ) :: nCols
   INTEGER      , INTENT(IN   ) :: kMax
   REAL(KIND=r8), INTENT(in   ) :: xdrag (ncols)
   REAL(KIND=r8), INTENT(in   ) :: ydrag (ncols)
   REAL(KIND=r8), INTENT(in   ) :: chug  (ncols,kmax)
   REAL(KIND=r8), INTENT(in   ) :: chvg  (ncols,kmax) 
   REAL(KIND=r8), INTENT(in   ) :: xtens (ncols,kmax+1)
   REAL(KIND=r8), INTENT(in   ) :: ytens (ncols,kmax+1)
   
   IF(dodia(nDiag_txgwdp))CALL updia(xtens (1:ncols,1:kmax),nDiag_txgwdp,latco)
   IF(dodia(nDiag_tygwdp))CALL updia(ytens (1:ncols,1:kmax),nDiag_tygwdp,latco)
   IF(dodia(nDiag_txgwds))CALL updia(xdrag,nDiag_txgwds,latco)
   IF(dodia(nDiag_tygwds))CALL updia(ydrag,nDiag_tygwds,latco)
   IF(dodia(nDiag_gwduzc))CALL updia(chug, nDiag_gwduzc,latco)
   IF(dodia(nDiag_gwdvmc))CALL updia(chvg, nDiag_gwdvmc,latco)

  END SUBROUTINE GwddDiagnStorage
  !----------------------------------------------------------------------
END MODULE GwddDriver
