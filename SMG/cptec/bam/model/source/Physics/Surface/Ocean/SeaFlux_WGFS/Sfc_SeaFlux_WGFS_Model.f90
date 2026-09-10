MODULE Sfc_SeaFlux_WGFS_Model
  IMPLICIT NONE
SAVE
  
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers
  !REAL(KIND=r8), PARAMETER   :: oml_hml0=200.0_r8
  !REAL(KIND=r8), PARAMETER   :: oml_hml0=50.0_r8


  INTEGER,PARAMETER:: NTYPE=9
  INTEGER,PARAMETER:: NGRID=22
  REAL(kind=r8)              DFK(NGRID,NTYPE),                            &
       &                     KTK(NGRID,NTYPE),                            &
       &                     DFKT(NGRID,NTYPE)
  !
  !  the nine soil types are:
  !    1  ... loamy sand (coarse) 
  !    2  ... silty clay loam (medium)
  !    3  ... light clay (fine) 
  !    4  ... sandy loam (coarse-medium) 
  !    5  ... sandy clay (coarse-fine)
  !    6  ... clay loam  (medium-fine)
  !    7  ... sandy clay loam (coarse-med-fine)
  !    8  ... loam  (organic)
  !    9  ... ice (use loamy sand property)
  !
  REAL(kind=r8), PARAMETER  :: b     (1:NTYPE)=(/4.26_r8,8.72_r8,11.55_r8,4.74_r8,10.73_r8,8.17_r8,6.77_r8,5.25_r8,4.26_r8/)
  REAL(kind=r8), PARAMETER  :: satpsi(1:NTYPE)=(/.04_r8,.62_r8,.47_r8,.14_r8,.10_r8,.26_r8,.14_r8,.36_r8,.04_r8/)
  REAL(kind=r8), PARAMETER  :: satkt (1:NTYPE)= (/1.41e-5_r8,.20e-5_r8,.10e-5_r8,.52e-5_r8,.72e-5_r8,&
       &           .25e-5_r8,.45e-5_r8,.34e-5_r8,1.41e-5_r8/)
  REAL(kind=r8), PARAMETER  :: TSAT  (1:NTYPE)=(/.421_r8,.464_r8,.468_r8,.434_r8,.406_r8,.465_r8,.404_r8,.439_r8,.421_r8/)
  !
  REAL(kind=r8),PUBLIC,PARAMETER::  CONKE_CH  =  1.0_r8    ! TUNABLE CONSTANT FOR EVAPORATION DUE  THERMAL  INSTABILITY
  REAL(kind=r8),PUBLIC,PARAMETER::  CONKE_CM  =  1.0_r8    ! TUNABLE CONSTANT FOR EVAPORATION DUE  DYNANICS INSTABILITY

CONTAINS
  SUBROUTINE Init_Sfc_SeaFlux_WGFS_Model()
    IMPLICIT NONE
    CALL GRDKT()
    CALL GRDDF()
  END SUBROUTINE Init_Sfc_SeaFlux_WGFS_Model


  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

  SUBROUTINE SeaFlux_WGFS_Model( &
       U3D       , &!   REAL, INTENT(IN   ) ::  U3D (1:nCols, 1:kMax) !-- U3D     3D u-velocity interpolated to theta points (m/s)
       V3D       , &!   REAL, INTENT(IN   ) ::  V3D (1:nCols, 1:kMax) !-- V3D     3D v-velocity interpolated to theta points (m/s)
       T3D       , &!   REAL, INTENT(IN   ) ::  T3D (1:nCols, 1:kMax) !-- T3D     temperature (K)
       QV3D      , &!   REAL, INTENT(IN   ) ::  QV3D(1:nCols, 1:kMax) !-- QV3D    3D water vapor mixing ratio (Kg/Kg)
       P3D       , &!   REAL, INTENT(IN   ) ::  P3D (1:nCols) !-- P3D     3D pressure (Pa)
       PSFC      , &!   REAL, INTENT(IN   ) ::  PSFC  (1:nCols)                    surface pressure (Pa)
       CHS       , &!   REAL, INTENT(OUT  ) ::  CHS  (1:nCols)
       CHS2      , &!   REAL, INTENT(OUT  ) ::  CHS2 (1:nCols)
       CQS2      , &!   REAL, INTENT(OUT  ) ::  CQS2 (1:nCols)
       CPM       , &!   REAL, INTENT(OUT  ) ::  CPM      REAL, DnColsNSION(1:nCols), INTENT(OUT) ::                 &
       ZNT       , &!   REAL, INTENT(INOUT) ::  ZNT  (1:nCols)
       UST       , &!   REAL, INTENT(INOUT) ::  UST  (1:nCols)
       PSIM      , &!   REAL, INTENT(OUT  ) ::  PSIM  (1:nCols)
       PSIH      , &!   REAL, INTENT(OUT  ) ::  PSIH  (1:nCols)
       XLAND     , &!   REAL, INTENT(IN   ) ::  XLAND (1:nCols)
       HFX       , &!   REAL, INTENT(OUT  ) ::  HFX   (1:nCols)
       QFX       , &!   REAL, INTENT(OUT  ) ::  QFX   (1:nCols)
       LH        , &!   REAL, INTENT(OUT  ) ::  LH    (1:nCols)
       TSK       , &!   REAL, INTENT(IN   ) ::  TSK   (1:nCols)
       FLHC      , &!   REAL, INTENT(OUT  ) ::  FLHC  (1:nCols)
       FLQC      , &!   REAL, INTENT(OUT  ) ::  FLQC  (1:nCols)
       QGH       , &!   REAL, INTENT(OUT  ) ::  QGH    (1:nCols),
       QSFC      , &!   REAL, INTENT(OUT  ) ::  QSFC    (1:nCols),
       U10       , &!   REAL, INTENT(OUT  ) ::  U10     (1:nCols),
       V10       , &!   REAL, INTENT(OUT  ) ::  V10     (1:nCols),
       GZ1OZ0    , &!   REAL, INTENT(OUT  ) ::  GZ1OZ0  (1:nCols),
       WSPD      , &!   REAL, INTENT(OUT  ) ::  WSPD    (1:nCols),
       BR        , &!   REAL, INTENT(OUT  ) ::  BR      (1:nCols),
       CHS_SEA   , &!   REAL, INTENT(OUT) ::, CHS_SEA  (1:nCols)
       CHS2_SEA  , &!   REAL, INTENT(OUT) ::, CHS2_SEA  (1:nCols)
       CPM_SEA   , &!   REAL, INTENT(OUT) ::, CPM_SEA  (1:nCols)
       CQS2_SEA  , &!   REAL, INTENT(OUT) ::, CQS2_SEA  (1:nCols)
       FLHC_SEA  , &!   REAL, INTENT(OUT) ::, FLHC_SEA  (1:nCols)
       FLQC_SEA  , &!   REAL, INTENT(OUT) ::, FLQC_SEA  (1:nCols)
       HFX_SEA   , &!   REAL, INTENT(OUT) ::, HFX_SEA  (1:nCols)
       LH_SEA    , &!   REAL, INTENT(OUT) ::, LH_SEA  (1:nCols)
       QFX_SEA   , &!   REAL, INTENT(OUT) ::, QFX_SEA  (1:nCols)
       QGH_SEA   , &!   REAL, INTENT(OUT) ::, QGH_SEA  (1:nCols)
       QSFC_SEA  , &!   REAL, INTENT(OUT) ::, QSFC_SEA  (1:nCols)
       UST_SEA   , &!   REAL, INTENT(OUT) ::, UST_SEA  (1:nCols)
       ZNT_SEA   , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
       CM_SEA    , &
       CH_SEA    , &
       WSPD_SEA  , &
       SST       , &!   REAL, INTENT(in ) ::, SST    (1:nCols)
       XICE      , &!   REAL, INTENT(in ) ::, XICE   (1:nCols)
       mskant    , &
       ztn       , &
       nCols       )
    IMPLICIT NONE  

    REAL(kind=r8), PARAMETER :: R       =2.8705e+2_r8 ! gas constant air    (J/kg/K)   
    REAL(kind=r8), PARAMETER :: r_v          = 461.6_r8! gas constant H2O    (J/kg/K)
    REAL(kind=r8), PARAMETER :: CP      =1.0046e+3_r8 ! spec heat air cp    (J/kg/K)

    REAL(kind=r8), PARAMETER :: ROVCP   =R/CP
    REAL(kind=r8), PARAMETER :: XLV          = 2.5E6_r8
    REAL(kind=r8), PARAMETER :: EP1=R_v/R-1.0_r8
    REAL(kind=r8), PARAMETER :: EP2=R/R_v
    REAL(kind=r8), PARAMETER ::  KARMAN=0.4_r8
    INTEGER      , PARAMETER ::  ISFFLX=1_r8   

    INTEGER         , INTENT(IN) :: nCols
    INTEGER(KIND=i8), INTENT(IN) :: mskant(nCols)
    REAL(KIND=r8),    INTENT(IN) :: ztn   (nCols)
    REAL(KIND=r8),    INTENT(IN) :: P3D   (1:nCols)
    REAL(KIND=r8),    INTENT(IN) :: QV3D  (1:nCols)
    REAL(KIND=r8),    INTENT(IN) :: T3D   (1:nCols)
    REAL(KIND=r8),    INTENT(IN) :: U3D   (1:nCols)
    REAL(KIND=r8),    INTENT(IN) :: V3D   (1:nCols)

    REAL(KIND=r8),    INTENT(IN   ) ::  TSK   (1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) ::  PSFC  (1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) ::  XLAND (1:nCols)
    REAL(KIND=r8),    INTENT(IN   ) ::  XICE  (1:nCols)

    REAL(KIND=r8),    INTENT(INOUT) :: UST  (1:nCols)
    REAL(KIND=r8),    INTENT(INOUT) :: ZNT  (1:nCols)
    REAL(KIND=r8),    INTENT(INOUT) :: SST  (1:nCols)

    REAL(KIND=r8), INTENT(OUT) ::  BR    (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  CHS   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  CHS2  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  CPM   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  CQS2  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  FLHC  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  FLQC  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  GZ1OZ0(1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  HFX   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  LH    (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  PSIM  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  PSIH  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  QFX   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  QGH   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  QSFC  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  U10   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  V10   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::  WSPD  (1:nCols)

    REAL(KIND=r8), INTENT(OUT) ::   CHS_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   CHS2_SEA (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   CPM_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   CQS2_SEA (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   FLHC_SEA (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   FLQC_SEA (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   HFX_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   LH_SEA   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   QFX_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   QGH_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   QSFC_SEA (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   UST_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   ZNT_SEA  (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   CM_SEA   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   CH_SEA   (1:nCols)
    REAL(KIND=r8), INTENT(OUT) ::   WSPD_SEA (1:nCols)



    !-------------------------------------------------------------------------
    !   Local
    !-------------------------------------------------------------------------
    INTEGER :: I
    REAL(KIND=r8) ::   BR_SEA      (1:nCols)
    REAL(KIND=r8) ::   GZ1OZ0_SEA  (1:nCols)
    REAL(KIND=r8) ::   PSIM_SEA    (1:nCols)
    REAL(KIND=r8) ::   PSIH_SEA    (1:nCols)
    REAL(KIND=r8) ::   U10_SEA     (1:nCols)
    REAL(KIND=r8) ::   V10_SEA     (1:nCols)
    REAL(KIND=r8) ::   XLAND_SEA   (1:nCols)
    REAL(KIND=r8) ::   TSK_SEA     (1:nCols)
    REAL(KIND=r8) ::   UST_HOLD    (1:nCols)
    REAL(KIND=r8) ::   ZNT_HOLD    (1:nCols)
    REAL(KIND=r8) ::   TSK_LOCAL   (1:nCols)
    REAL(KIND=r8) ::   CM          (1:nCols)
    REAL(KIND=r8) ::   CH          (1:nCols)
    REAL(KIND=r8) ::   XICE_THRESHOLD

    !
    ! Set up for frozen ocean call for sea ice points
    !
    BR          =0.0_r8;   CHS         =0.0_r8;   CHS2        =0.0_r8;   CPM         =0.0_r8
    CQS2        =0.0_r8;   FLHC        =0.0_r8;   FLQC        =0.0_r8;   GZ1OZ0      =0.0_r8
    HFX         =0.0_r8;   LH          =0.0_r8;   PSIM        =0.0_r8;   PSIH        =0.0_r8
    QFX         =0.0_r8;   QGH         =0.0_r8;   QSFC        =0.0_r8;   U10         =0.0_r8
    V10         =0.0_r8;   WSPD        =0.0_r8;   CHS_SEA     =0.0_r8;   CHS2_SEA    =0.0_r8
    CPM_SEA     =0.0_r8;   CQS2_SEA    =0.0_r8;   FLHC_SEA    =0.0_r8;   FLQC_SEA    =0.0_r8
    HFX_SEA     =0.0_r8;   LH_SEA      =0.0_r8;   QFX_SEA     =0.0_r8;   QGH_SEA     =0.0_r8
    QSFC_SEA    =0.0_r8;   UST_SEA     =0.0_r8;   ZNT_SEA     =0.0_r8;   CM_SEA      =0.0_r8
    CH_SEA      =0.0_r8;   WSPD_SEA    =0.0_r8;   BR_SEA      =0.0_r8;   GZ1OZ0_SEA  =0.0_r8
    PSIM_SEA    =0.0_r8;   PSIH_SEA    =0.0_r8;   U10_SEA     =0.0_r8;   V10_SEA     =0.0_r8
    XLAND_SEA   =0.0_r8;   TSK_SEA     =0.0_r8;   TSK_LOCAL   =0.0_r8;   CM          =0.0_r8
    CH          =0.0_r8
    ! Strictly INTENT(IN), Should be unchanged by SF_GFS:
    !     CP
    !     EP1
    !     EP2
    !     KARMAN
    !     R
    !     ROVCP
    !     XLV
    !     P3D
    !     QV3D
    !     T3D
    !     U3D
    !     V3D
    !     TSK
    !     PSFC
    !     XLAND
    !     ISFFLX
    !     ITIMESTEP


    ! Intent (INOUT), original value is used and changed by SF_GFS.
    !     UST
    !     ZNT

    ZNT_HOLD = ZNT
    UST_HOLD = UST
    !if ( fractional_seaice == 0 ) then
    xice_threshold = 0.5_r8
    !else if ( fractional_seaice == 1 ) then
    !  xice_threshold = 0.02
    !endif

    ! Strictly INTENT (OUT), set by SF_GFS:
    !     BR
    !     CHS     -- used by LSM routines
    !     CHS2    -- used by LSM routines
    !     CPM     -- used by LSM routines
    !     CQS2    -- used by LSM routines
    !     FLHC
    !     FLQC
    !     GZ1OZ0
    !     HFX     -- used by LSM routines
    !     LH      -- used by LSM routines
    !     PSIM
    !     PSIH
    !     QFX     -- used by LSM routines
    !     QGH     -- used by LSM routines
    !     QSFC    -- used by LSM routines
    !     U10
    !     V10
    !     WSPD

    !
    ! Frozen ocean / true land call.
    !
    CALL SF_GFS(U3D,V3D,T3D,QV3D,P3D,                  &
         CP,ROVCP,R,XLV,PSFC,CHS,CHS2,CQS2,CPM_SEA,    &
         ZNT,UST,PSIM,PSIH,                            &
         XLAND,HFX,QFX,LH,TSK_LOCAL,FLHC,FLQC,         &
         QGH,QSFC,U10,V10,                             &
         GZ1OZ0,WSPD,BR,ISFFLX,                        &
         EP1,EP2,KARMAN,CM,CH,mskant, ztn,                    &
         nCols)

    ! Set up for open-water call

    !     DO j = JTS , JTE
    DO i = 1 , nCols
       IF(mskant(i) == 1_i8)THEN
          IF ( ( XICE(I) .GE. XICE_THRESHOLD ) .AND. ( XICE(i) .LE. 1.0_r8 ) ) THEN
             ! Sets up things for open ocean fraction of sea-ice points
             XLAND_SEA(i)=2.0_r8
             ZNT_SEA(I) = 0.0001_r8
             IF ( SST(i) .LT. 271.4_r8 ) THEN
                SST(i) = 271.4_r8
             ENDIF
             TSK_SEA(i) = SST(i)
          ELSE
             ! Fully open ocean or true land points
             XLAND_SEA(i)=xland(i)
             ZNT_SEA(I) = ZNT_HOLD(I)
             UST_SEA(i) = UST_HOLD(i)
             TSK_SEA(i) = TSK(i)
          ENDIF
       END IF
    ENDDO
    !    ENDDO
    ! Open-water call
    ! _SEA variables are held for later use as the result of the open-water call.
    CALL SF_GFS(U3D,V3D,T3D,QV3D,P3D,                  &
         CP,ROVCP,R,XLV,PSFC,CHS_SEA,CHS2_SEA,CQS2_SEA,CPM,        &
         ZNT_SEA,UST_SEA,PSIM_SEA,PSIH_SEA,                        &
         XLAND,HFX_SEA,QFX_SEA,LH_SEA,TSK_SEA,FLHC_SEA,FLQC_SEA,   &
         QGH_SEA,QSFC_SEA,U10_SEA,V10_SEA,                         &
         GZ1OZ0_SEA,WSPD_SEA,BR_SEA,ISFFLX,                        &
         EP1,EP2,KARMAN,CM_SEA,CH_SEA,mskant ,ztn,                    &
         nCols)

    ! Weighting, after our two calls to SF_GFS

    !     DO j = JTS , JTE
    DO i = 1 , nCols       
       IF(mskant(i) == 1_i8)THEN
          ! Over sea-ice points, weight the results.  Otherwise, just take the results from the
          ! first call to SF_GFS_
          IF ( ( XICE(I) .GE. XICE_THRESHOLD ) .AND. ( XICE(i) .LE. 1.0_r8 ) ) THEN
             ! Weight a number of fields (between open-water results
             ! and full ice results) by sea-ice fraction.

             BR(i)     = ( BR(i)     * XICE(i) ) + ( (1.0_r8-XICE(i)) * BR_SEA(i)     )
             ! CHS, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! CHS2, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! CPM, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! CQS2, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! FLHC(i,j) = ( FLHC(i,j) * XICE(i,j) ) + ( (1.0-XICE(i,j)) * FLHC_SEA(i,j) )
             ! FLQC(i,j) = ( FLQC(i,j) * XICE(i,j) ) + ( (1.0-XICE(i,j)) * FLQC_SEA(i,j) )
             GZ1OZ0(i) = ( GZ1OZ0(i) * XICE(i) ) + ( (1.0_r8-XICE(i)) * GZ1OZ0_SEA(i) )
             ! HFX, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! LH, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             PSIM(i)   = ( PSIM(i)   * XICE(i) ) + ( (1.0_r8-XICE(i)) * PSIM_SEA(i)   )
             PSIH(i)   = ( PSIH(i)   * XICE(i) ) + ( (1.0_r8-XICE(i)) * PSIH_SEA(i)   )
             ! QFX, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! QGH, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! QSFC, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             U10(i)    = ( U10(i)    * XICE(i) ) + ( (1.0_r8-XICE(i)) * U10_SEA(i)    )
             V10(i)    = ( V10(i)    * XICE(i) ) + ( (1.0_r8-XICE(i)) * V10_SEA(i)    )
             WSPD(i)   = ( WSPD(i)   * XICE(i) ) + ( (1.0_r8-XICE(i)) * WSPD_SEA(i)   )
             ! UST, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
             ! ZNT, used by the LSM routines, is not updated yet.  Return results from both calls in separate variables
          END IF
       ENDIF
    ENDDO
    !     ENDDO

  END SUBROUTINE SeaFlux_WGFS_Model


  !-------------------------------------------------------------------
  SUBROUTINE SF_GFS(U3D,V3D,T3D,QV3D,P3D,                              &
       CP,ROVCP,R,XLV,PSFC,CHS,CHS2,CQS2,CPM,             &
       ZNT,UST,PSIM,PSIH,                                 &
       XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,                    &
       QGH,QSFC,U10,V10,                                  &
       GZ1OZ0,WSPD,BR,ISFFLX,                             &
       EP1,EP2,KARMAN,CM,CH,mskant ,  ztn,                       &
       nCols)
    !-------------------------------------------------------------------

    USE PhysicalFunctions , ONLY : fpvs
    !-------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    !-- U3D            3D u-velocity interpolated to theta points (m/s)
    !-- V3D            3D v-velocity interpolated to theta points (m/s)
    !-- T3D             temperature (K)
    !-- QV3D        3D water vapor mixing ratio (Kg/Kg)
    !-- P3D         3D pressure (Pa)
    !-- CP                heat capacity at constant pressure for dry air (J/kg/K)
    !-- ROVCP       R/CP
    !-- R           gas constant for dry air (J/kg/K)
    !-- XLV         latent heat of vaporization for water (J/kg)
    !-- PSFC        surface pressure (Pa)
    !-- ZNT                roughness length (m)
    !-- UST                u* in similarity theory (m/s)
    !-- PSIM        similarity stability function for momentum
    !-- PSIH        similarity stability function for heat
    !-- XLAND        land mask (1 for land, 2 for water)
    !-- HFX                upward heat flux at the surface (W/m^2)
    !-- QFX           upward moisture flux at the surface (kg/m^2/s)
    !-- LH          net upward latent heat flux at surface (W/m^2)
    !-- TSK                surface temperature (K)
    !-- FLHC        exchange coefficient for heat (m/s)
    !-- FLQC        exchange coefficient for moisture (m/s)
    !-- QGH         lowest-level saturated mixing ratio
    !-- U10         diagnostic 10m u wind
    !-- V10         diagnostic 10m v wind
    !-- GZ1OZ0      log(z/z0) where z0 is roughness length
    !-- WSPD        wind speed at lowest model level (m/s)
    !-- BR          bulk Richardson number in surface layer
    !-- ISFFLX      isfflx=1 for surface heat and moisture fluxes
    !-- EP1         constant for virtual temperature (R_v/R_d - 1) (dimensionless)
    !-- KARMAN      Von Karman constant
    !-- ids         start index for i in domain
    !--          end index for i in domain
    !-- jds         start index for j in domain
    !--         end index for j in domain
    !-- kds         start index for k in domain
    !-- ims         start index for i in memory
    !-- nCols         end index for i in memory
    !-- jms         start index for j in memory
    !--          end index for j in memory
    !-- kms         start index for k in memory
    !-- kMax         end index for k in memory
    !--          start index for i in tile
    !--          end index for i in tile
    !-- jts         start index for j in tile
    !-- jte         end index for j in tile
    !--          start index for k in tile
    !--          end index for k in tile
    !-------------------------------------------------------------------

    INTEGER, INTENT(IN)          ::  nCols
    INTEGER, INTENT(IN)          ::  ISFFLX
    REAL(KIND=r8),    INTENT(IN) ::  ztn(nCols)
    INTEGER(KIND=i8), INTENT(IN) ::  mskant(nCols)
    REAL(KIND=r8),    INTENT(IN) ::  CP
    REAL(KIND=r8),    INTENT(IN) ::  EP1
    REAL(KIND=r8),    INTENT(IN) ::  EP2
    REAL(KIND=r8),    INTENT(IN) ::  KARMAN
    REAL(KIND=r8),    INTENT(IN) ::  R 
    REAL(KIND=r8),    INTENT(IN) ::  ROVCP
    REAL(KIND=r8),    INTENT(IN) ::  XLV

    REAL(KIND=r8), INTENT(IN)   :: P3D  (1:nCols)
    REAL(KIND=r8), INTENT(IN)   :: QV3D (1:nCols)
    REAL(KIND=r8), INTENT(IN)   :: T3D  (1:nCols)
    REAL(KIND=r8), INTENT(IN)   :: U3D  (1:nCols)
    REAL(KIND=r8), INTENT(IN)   :: V3D  (1:nCols)

    REAL(KIND=r8), INTENT(IN   ) ::     TSK  (1:nCols)
    REAL(KIND=r8), INTENT(IN   ) ::     PSFC  (1:nCols)
    REAL(KIND=r8), INTENT(IN   ) ::     XLAND (1:nCols)

    REAL(KIND=r8), INTENT(INOUT) :: UST (1:nCols)
    REAL(KIND=r8), INTENT(INOUT) :: ZNT (1:nCols)

    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        BR    (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        CHS   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        CHS2  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        CPM   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        CQS2  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        FLHC  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        FLQC  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        GZ1OZ0 (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        HFX   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        LH    (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        PSIM  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        PSIH  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        QFX   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        QGH   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        QSFC  (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        U10   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        V10   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        WSPD   (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        CH     (1:nCols)
    REAL(KIND=r8), DIMENSION(1:nCols), INTENT(OUT) ::        CM     (1:nCols)


    !--------------------------- LOCAL VARS ------------------------------

    REAL     (KIND=r8) ::     ESAT
    REAL     (kind=r8) ::     DDVEL(1:nCols)
    REAL     (kind=r8) ::     DRAIN(1:nCols)
    REAL     (kind=r8) ::     EP(1:nCols)
    REAL     (kind=r8) ::     EVAP(1:nCols)
    REAL     (kind=r8) ::     FH(1:nCols)
    REAL     (kind=r8) ::     FH2(1:nCols)
    REAL     (kind=r8) ::     FM(1:nCols)
    REAL     (kind=r8) ::     PH(1:nCols)
    REAL     (kind=r8) ::     PM(1:nCols)
    REAL     (kind=r8) ::     PRSL1(1:nCols)
    REAL     (kind=r8) ::     PRSLKI(1:nCols)
    REAL     (kind=r8) ::     PS(1:nCols)
    REAL     (kind=r8) ::     Q1(1:nCols)
    REAL     (kind=r8) ::     Q2M(1:nCols)
    REAL     (kind=r8) ::     QSS(1:nCols)
    REAL     (kind=r8) ::     RB(1:nCols)
    REAL     (kind=r8) ::     RCL(1:nCols)
    REAL     (kind=r8) ::     RHO1(1:nCols)
    REAL     (kind=r8) ::     SLIMSK(1:nCols)
    REAL     (kind=r8) ::     STRESS(1:nCols)
    REAL     (kind=r8) ::     T1(1:nCols)
    REAL     (kind=r8) ::     T2M(1:nCols)
    REAL     (kind=r8) ::     THGB(1:nCols)
    REAL     (kind=r8) ::     THX(1:nCols)
    REAL     (kind=r8) ::     TSKIN(1:nCols)
    REAL     (kind=r8) ::     U1(1:nCols)
    REAL     (kind=r8) ::     U10M(1:nCols)
    REAL     (kind=r8) ::     USTAR(1:nCols)
    REAL     (kind=r8) ::     V1(1:nCols)
    REAL     (kind=r8) ::     V10M(1:nCols)
    REAL     (kind=r8) ::     WIND(1:nCols)
    REAL     (kind=r8) ::     Z0RL(1:nCols)
    REAL     (kind=r8) ::     Z1(1:nCols)


    INTEGER :: I
    INTEGER :: IM
    INTEGER :: KM


         ESAT=0.0_r8

         BR=0.0_r8
         CHS=0.0_r8   
         CHS2=0.0_r8  
         CPM=0.0_r8   
         CQS2=0.0_r8  
         FLHC=0.0_r8  
         FLQC=0.0_r8  
         GZ1OZ0=0.0_r8
         HFX=0.0_r8   
         LH=0.0_r8    
         PSIM=0.0_r8  
         PSIH=0.0_r8  
         QFX=0.0_r8   
         QGH=0.0_r8   
         QSFC=0.0_r8  
         U10=0.0_r8   
         V10=0.0_r8   
         WSPD=0.0_r8  
         CH=0.0_r8            
         CM=0.0_r8            
         DDVEL=0.0_r8         
         DRAIN=0.0_r8         
         EP=0.0_r8            
         EVAP=0.0_r8          
         FH=0.0_r8            
         FH2=0.0_r8           
         FM=0.0_r8            
         PH=0.0_r8            
         PM=0.0_r8            
         PRSL1=0.0_r8         
         PRSLKI=0.0_r8        
         PS=0.0_r8            
         Q1=0.0_r8            
         Q2M=0.0_r8           
         QSS=0.0_r8           
         RB=0.0_r8            
         RCL=0.0_r8           
         RHO1=0.0_r8          
         SLIMSK=0.0_r8        
         STRESS=0.0_r8        
         T1=0.0_r8            
         T2M=0.0_r8           
         THGB=0.0_r8          
         THX=0.0_r8           
         TSKIN=0.0_r8         
         U1=0.0_r8            
         U10M=0.0_r8          
         USTAR=0.0_r8         
         V1=0.0_r8            
         V10M=0.0_r8          
         WIND=0.0_r8          
         Z0RL=0.0_r8          
         Z1=0.0_r8          
    !   if(itimestep.eq.0) then
    !     CALL GFUNCPHYS
    !   endif

    IM=nCols-1+1
    KM=1

    !   DO J=jts,jte

    DO i=1,nCols
       IF(mskant(i) == 1.0_r8)THEN
          DDVEL(I)=0.0_r8
          RCL(i)=1.0_r8
          PRSL1(i)=P3D(i)*0.001_r8
          PS(i)=PSFC(i)*.001_r8
          Q1(I) = QV3D(i)
          !        QSURF(I)=QSFC(I)
          !-- XLAND        land mask (1 for land, 2 for water)

          !SLIMSK(i)=ABS(XLAND(i)-2.0_r8) !land mask -1.0 for land, 0.0 for water)
          SLIMSK(i)=0.0_r8!         ABS(XLAND(i)-2.0_r8) !land mask -1.0 for land, 0.0 for water)

          TSKIN(i )=TSK(i)
          T1(I)    = T3D(i)
          U1(I)    = U3D(i)
          USTAR(I) = UST(i)
          V1(I)    = V3D(i)
          Z0RL(I)  = ZNT(i)*100.0_r8
       END IF
    ENDDO

    DO i=1,nCols       
       IF(mskant(i) == 1.0_r8)THEN
          PRSLKI(i)=(PS(I)/PRSL1(I))**ROVCP
          THGB(I)=TSKIN(i)*(100.0_r8/PS(I))**ROVCP
          THX(I)=T1(i)*(100.0_r8/PRSL1(I))**ROVCP
          RHO1(I)=PRSL1(I)*1000.0_r8/(R*T1(I)*(1.0_r8+EP1*Q1(I)))
          Q1(I)=Q1(I)/(1.0_r8+Q1(I))
       END IF
    ENDDO


    CALL PROGTM(IM,KM,PS,U1,V1,T1,Q1,                                 &
         TSKIN,                                   &
         Z0RL,                                                 &
         U10M,V10M,T2M,Q2M,                                    &
         CM,CH,RB,                                             &
         RCL,PRSL1,PRSLKI,SLIMSK,                              &
         DRAIN,EVAP,STRESS,EP,                            &
         FM,FH,USTAR,WIND,DDVEL,                               &
         PM,PH,FH2,QSS,Z1 ,mskant    ,ztn                                 )


    DO i=1,nCols       
       IF(mskant(i)  == 1.0_r8)THEN
          U10(i)=U10M(i)
          V10(i)=V10M(i)
          BR(i)=RB(i)
          CHS(I)=CH(I)*WIND(I)
          CHS2(I)=USTAR(I)*KARMAN/FH2(I)
          CPM(I)=CP*(1.0_r8+0.8_r8*QV3D(i))
          esat = fpvs(t1(i))
          QGH(I)=ep2*esat/(1000.0_r8*ps(i)-esat)
          QSFC(I)=qss(i)
          PSIH(i)=PH(i)
          PSIM(i)=PM(i)
          UST(i)=ustar(i)
          WSPD(i)=WIND(i)
          ZNT(i)=Z0RL(i)*.01_r8
       END IF
    ENDDO

    DO i=1,nCols    
       IF(mskant(i) == 1.0_r8)THEN

          FLHC(i)=CPM(I)*RHO1(I)*CHS(I)
          FLQC(i)=       RHO1(I)*CHS(I)
          GZ1OZ0(i)=LOG(Z1(I)/(Z0RL(I)*.01_r8))
          CQS2(i)=CHS2(I)
       END IF
    ENDDO

    IF (ISFFLX.EQ.0) THEN
       DO i=1,nCols    
          IF(mskant(i) == 1.0_r8)THEN
             HFX(i)=0.0_r8
             LH(i)=0.0_r8
             QFX(i)=0.0_r8
          END IF
       ENDDO
    ELSE
       DO i=1,nCols         
          IF(mskant(i) == 1.0_r8)THEN
             IF(XLAND(I)-1.5_r8.GT.0._r8)THEN
                HFX(I)=FLHC(I)*(THGB(I)-THX(I))
             ELSEIF(XLAND(I)-1.5_r8.LT.0.0_r8)THEN
                HFX(I)=FLHC(I)*(THGB(I)-THX(I))
                HFX(I)=MAX(HFX(I),-250.0_r8)
             ENDIF
             QFX(I)=FLQC(I)*(QSFC(I)-Q1(I))
             QFX(I)=QFX(I)
             LH(I)=XLV*QFX(I)          
            ! PRINT*,'LH(I)',LH(I) ,HFX(i) , CH(i) , CM(I)
          END IF
       ENDDO
    ENDIF


    !   ENDDO


  END SUBROUTINE SF_GFS




  !-------------------------------------------------------------------

  SUBROUTINE PROGTM( &
       IM         , &!INTEGER  , INTENT(IN   ) :: IM, km
       KM         , &!INTEGER  , INTENT(IN   ) :: km
       PS         , &! REAL(kind=r8)   , INTENT(IN   ) :: PS(IM) 
       U1         , &! REAL(kind=r8)   , INTENT(IN   ) :: U1(IM)
       V1         , &! REAL(kind=r8)   , INTENT(IN   ) :: V1(IM) 
       T1         , &! REAL(kind=r8)   , INTENT(IN   ) :: T1(IM)
       Q1         , &! REAL(kind=r8)   , INTENT(IN   ) :: Q1(IM)
       TSKIN      , &! REAL(kind=r8)   , INTENT(IN   ) :: TSKIN(IM)
       Z0RL       , &! REAL(kind=r8)   , INTENT(INOUT) :: Z0RL(IM)
       U10M       , &! REAL(kind=r8)   , INTENT(OUT  ) :: U10M(IM)
       V10M       , &! REAL(kind=r8)   , INTENT(OUT  ) :: V10M(IM)
       T2M        , &! REAL(kind=r8)   , INTENT(OUT  ) :: T2M(IM)
       Q2M        , &! REAL(kind=r8)   , INTENT(OUT  ) :: Q2M(IM)
       CM         , &! REAL(kind=r8)   , INTENT(OUT  ) :: CM(IM)
       CH         , &! REAL(kind=r8)   , INTENT(OUT  ) :: CH(IM)
       RB         , &! REAL(kind=r8)   , INTENT(OUT  ) :: RB(IM)
       RCL        , &! REAL(kind=r8)   , INTENT(IN   ) ::RCL(IM)
       PRSL1      , &! REAL(kind=r8)   , INTENT(IN   ) ::PRSL1(IM)
       PRSLKI     , &! REAL(kind=r8)   , INTENT(IN   ) ::PRSLKI(IM)
       SLIMSK     , &! REAL(kind=r8)   , INTENT(IN   ) ::SLIMSK(IM)
       DRAIN      , &! REAL(kind=r8)   , INTENT(INOUT) ::DRAIN(IM)
       EVAP       , &! REAL(kind=r8)   , INTENT(OUT  ) ::EVAP(IM)
       STRESS     , &! REAL(kind=r8)   , INTENT(INOUT) ::STRESS(IM)
       EP         , &! REAL(kind=r8)   , INTENT(OUT  ) ::EP(IM)
       FM         , &! REAL(kind=r8)   , INTENT(OUT  ) ::FM(IM)
       FH         , &! REAL(kind=r8)   , INTENT(OUT  ) ::FH(IM)
       USTAR      , &! REAL(kind=r8)   , INTENT(INOUT) ::USTAR(IM)
       WIND       , &! REAL(kind=r8)   , INTENT(OUT  ) ::WIND(IM)
       DDVEL      , &! REAL(kind=r8)   , INTENT(INOUT) ::DDVEL(IM)
       PM         , &! REAL(kind=r8)   , INTENT(OUT  ) :: PM(IM)
       PH         , &! REAL(kind=r8)   , INTENT(OUT  ) :: PH(IM)
       FH2        , &! REAL(kind=r8)   , INTENT(OUT  ) :: FH2(IM)
       QSS        , &! REAL(kind=r8)   , INTENT(OUT  ) ::QSS(IM)    
       Z1         , &! REAL(kind=r8)   , INTENT(OUT  ) ::Z1(IM)    
       mskant      , &
       ztn                            )
    !
    USE PhysicalFunctions, ONLY : fpvs
    IMPLICIT NONE      
    REAL(kind=r8),PARAMETER :: con_rd      =2.8705e+2_r8 ! gas constant air    (J/kg/K)
    REAL(kind=r8),PARAMETER :: con_rv      =4.6150e+2_r8 ! gas constant H2O    (J/kg/K)

    REAL(kind=r8),PARAMETER:: RD = con_RD    ! gas constant air    (J/kg/K)

    REAL(kind=r8),PARAMETER:: grav     =9.80665e+0_r8! gravity             (m/s2)
    REAL(kind=r8),PARAMETER:: SBC      =5.6730e-8_r8 ! stefan-boltzmann    (W/m2/K4)
    REAL(kind=r8),PARAMETER:: HVAP     =2.5000e+6_r8 ! lat heat H2O cond   (J/kg)
    REAL(kind=r8),PARAMETER:: CP       =1.0046e+3_r8 ! spec heat air @p    (J/kg/K)
    REAL(kind=r8),PARAMETER:: HFUS     =3.3358e+5_r8 ! lat heat H2O fusion (J/kg)
    REAL(kind=r8),PARAMETER:: JCAL     =4.1855E+0_r8 ! JOULES PER CALORIE  ()
    REAL(kind=r8),PARAMETER:: EPS      =con_rd/con_rv
    REAL(kind=r8),PARAMETER:: EPSM1    =con_rd/con_rv-1.0_r8
    REAL(kind=r8),PARAMETER:: t0c      =2.7315e+2_r8 ! temp at 0C          (K)
    REAL(kind=r8),PARAMETER:: RVRDM1   =con_rv/con_rd-1.0_r8
    REAL(kind=r8), PARAMETER :: cpinv=1.0_r8/cp, HVAPI=1.0_r8/HVAP

    !
    !     include 'constant.h'
    !
    INTEGER  , INTENT(IN   ) :: IM, km
    !
    REAL(kind=r8)   , INTENT(IN   ) :: ztn(im)
    INTEGER(KINd=i8), INTENT(IN   ) :: mskant(Im)
    REAL(kind=r8)   , INTENT(IN   ) :: PS(IM) 
    REAL(kind=r8)   , INTENT(IN   ) :: U1(IM)
    REAL(kind=r8)   , INTENT(IN   ) :: V1(IM) 
    REAL(kind=r8)   , INTENT(IN   ) ::  T1(IM)
    REAL(kind=r8)   , INTENT(IN   ) ::  Q1(IM)
    REAL(kind=r8)   , INTENT(IN   ) ::  TSKIN(IM)
    REAL(kind=r8)   , INTENT(INOUT) :: Z0RL(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: U10M(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: V10M(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: T2M(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: Q2M(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: CM(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: CH(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: RB(IM)
    REAL(kind=r8)   , INTENT(IN   ) ::PRSL1(IM)
    REAL(kind=r8)   , INTENT(IN   ) ::PRSLKI(IM)
    REAL(kind=r8)   , INTENT(IN   ) ::SLIMSK(IM)
    REAL(kind=r8)   , INTENT(INOUT) ::DRAIN(IM)
    REAL(kind=r8)   , INTENT(OUT  ) ::EVAP(IM)
    REAL(kind=r8)   , INTENT(OUT  ) ::EP(IM)
    REAL(kind=r8)   , INTENT(OUT  ) ::FM(IM)
    REAL(kind=r8)   , INTENT(OUT  ) ::FH(IM)
    REAL(kind=r8)   , INTENT(INOUT) ::USTAR(IM)
    REAL(kind=r8)   , INTENT(OUT  ) ::WIND(IM)
    REAL(kind=r8)   , INTENT(INOUT) ::DDVEL(IM)
    REAL(kind=r8)   , INTENT(INOUT) ::STRESS(IM)
    REAL(kind=r8)   , INTENT(IN   ) :: RCL(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: PH(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: PM(IM)
    REAL(kind=r8)   , INTENT(OUT  ) :: FH2(IM)
    REAL(kind=r8)    ::PH2(IM)
    REAL(kind=r8)   , INTENT(OUT  ) ::QSS(IM)   
    REAL(kind=r8)   , INTENT(OUT  ) ::Z1(IM)    
    REAL(kind=r8)    :: DM(IM)
    REAL(kind=r8)    :: SNOWMT(IM)
    REAL(kind=r8)    :: GFLUX(IM)
    REAL(kind=r8)    :: F10M(IM)
    REAL(kind=r8)    :: RHSCNPY(IM)
    REAL(kind=r8)    :: RHSMC(IM,KM)
    REAL(kind=r8)   ::AIM(IM,KM)
    REAL(kind=r8)   ::BIM(IM, KM)
    REAL(kind=r8)   ::CIM(IM, KM)
    !
    !     Locals
    !
    INTEGER        ::      k,i
    !
    REAL(kind=r8) :: CANFAC(IM),                                  &
         &                     DTV(IM),     EC(IM),            &
         &                     EDIR(IM),    ETPFAC(IM),        &
         &                          FM10(IM),          &
         &                     FX(IM),                             &
         &                     HL1(IM),     HL12(IM),          &
         &                     HLINF(IM),   PM10(IM),          &
         &                     PSURF(IM),   Q0(IM),      QS1(IM),           &
         &                      RAT(IM),            &
         &                     RHO(IM),     RS(IM),            &
         &                     THETA1(IM),  THV1(IM),          &
         &                     TSURF(IM),   TV1(IM),           &
         &                     TVS(IM),                        &
         &                      XRCL(IM),                      &
         &                     Z0(IM),      Z0MAX(IM),           &
         &                     ZTMAX(IM),   PS1(IM)

    REAL(kind=r8) :: OLINF
    REAL(kind=r8) :: RESTAR
    REAL(kind=r8) :: HL0
    REAL(kind=r8) :: HL0INF
    REAL(kind=r8) :: HL110
    REAL(kind=r8) :: HLT
    REAL(kind=r8) :: SIG2K
    REAL(kind=r8) :: AA
    REAL(kind=r8) :: BB0
    REAL(kind=r8) :: CQ
    REAL(kind=r8) :: HLTINF
    REAL(kind=r8) :: FHS
    REAL(kind=r8) :: FMS
    REAL(kind=r8) :: BB
    REAL(kind=r8) :: AA0
    REAL(kind=r8) :: ADTV
    REAL(KIND=r8), PARAMETER ::  CHARNOCK=.014_r8
    REAL(KIND=r8), PARAMETER ::  CA=.4_r8          !C CA IS THE VON KARMAN CONSTANT
    REAL(KIND=r8), PARAMETER ::  G=grav
    REAL(KIND=r8), PARAMETER ::  sigma=sbc


    REAL(KIND=r8), PARAMETER ::  ALPHA=5._r8
    REAL(KIND=r8), PARAMETER ::  A0=-3.975_r8
    REAL(KIND=r8), PARAMETER ::  A1=12.32_r8
    REAL(KIND=r8), PARAMETER ::  B1=-7.755_r8
    REAL(KIND=r8), PARAMETER ::  B2=6.041_r8
    REAL(KIND=r8), PARAMETER ::  A0P=-7.941_r8
    REAL(KIND=r8), PARAMETER ::  A1P=24.75_r8
    REAL(KIND=r8), PARAMETER ::  B1P=-8.705_r8
    REAL(KIND=r8), PARAMETER ::  B2P=7.899_r8
    REAL(KIND=r8), PARAMETER ::  VIS=1.4E-5_r8
    REAL(KIND=r8), PARAMETER ::  AA1=-1.076_r8
    REAL(KIND=r8), PARAMETER ::  BB1=.7045_r8
    REAL(KIND=r8), PARAMETER ::  CC1=-.05808_r8
    REAL(KIND=r8), PARAMETER ::  BB2=-.1954_r8
    REAL(KIND=r8), PARAMETER ::  CC2=.009999_r8
    REAL(KIND=r8), PARAMETER ::  ELOCP=HVAP/CP
    REAL(KIND=r8), PARAMETER ::  DFSNOW=.31_r8
    REAL(KIND=r8), PARAMETER ::  CH2O=4.2E6_r8
    REAL(KIND=r8), PARAMETER ::  CSOIL=1.26E6_r8
    REAL(KIND=r8), PARAMETER ::  SCANOP=.5_r8
    REAL(KIND=r8), PARAMETER ::  CFACTR=.5_r8
    REAL(KIND=r8), PARAMETER ::  ZBOT=-3._r8
    REAL(KIND=r8), PARAMETER ::  TGICE=271.2_r8
    REAL(KIND=r8), PARAMETER ::  CICE=1880.*917._r8
    REAL(KIND=r8), PARAMETER ::  topt=298._r8
    REAL(KIND=r8), PARAMETER ::  RHOH2O=1000._r8
    REAL(KIND=r8), PARAMETER ::  CONVRAD=JCAL*1.E4/60._r8
    REAL(KIND=r8), PARAMETER ::  CTFIL1=.5_r8
    REAL(KIND=r8), PARAMETER ::  CTFIL2=1.-CTFIL1
    REAL(KIND=r8), PARAMETER ::  RNU=1.51E-5_r8
    REAL(KIND=r8), PARAMETER ::  ARNU=.135_r8*RNU
    REAL(KIND=r8), PARAMETER ::  snomin=1.0e-9_r8
    !

    !
    !  the 13 vegetation types are:
    !
    !  1  ...  broadleave-evergreen trees (tropical forest)
    !  2  ...  broadleave-deciduous trees
    !  3  ...  broadleave and needle leave trees (mixed forest)
    !  4  ...  needleleave-evergreen trees
    !  5  ...  needleleave-deciduous trees (larch)
    !  6  ...  broadleave trees with groundcover (savanna)
    !  7  ...  groundcover only (perenial)
    !  8  ...  broadleave shrubs with perenial groundcover
    !  9  ...  broadleave shrubs with bare soil
    ! 10  ...  dwarf trees and shrubs with ground cover (trunda)
    ! 11  ...  bare soil
    ! 12  ...  cultivations (use parameters from type 7)
    ! 13  ...  glacial
    !
    REAL(KIND=r8), PARAMETER :: rsmax(13)=(/5000._r8,5000._r8,5000._r8,5000._r8,5000._r8,5000._r8, &
         5000._r8,5000._r8,5000._r8,5000._r8,5000._r8,5000._r8,5000._r8/)
    REAL(KIND=r8), PARAMETER :: rsmin(13)=(/150._r8,100._r8,125._r8,150._r8,100._r8,70._r8,40._r8,                      &
         &           300._r8,400._r8,150._r8,999._r8,40._r8,999._r8/)
    REAL(KIND=r8), PARAMETER ::  rgl(13)=(/30._r8,30._r8,30._r8,30._r8,30._r8,65._r8,100._r8,100._r8,100._r8,100._r8,999._r8,100._r8,999._r8/)
    REAL(KIND=r8), PARAMETER ::  hs(13)=(/41.69_r8,54.53_r8,51.93_r8,47.35_r8,47.35_r8,54.53_r8,36.35_r8, &
         &                               42.00_r8,42.00_r8,42.00_r8,999._r8,36.35_r8,999._r8/)
    REAL(KIND=r8), PARAMETER ::  smmax(9)=(/.421_r8,.464_r8,.468_r8,.434_r8,.406_r8,.465_r8,.404_r8,.439_r8,.421_r8/)
    REAL(KIND=r8), PARAMETER ::  smdry(9)=(/.07_r8,.14_r8,.22_r8,.08_r8,.18_r8,.16_r8,.12_r8,.10_r8,.07_r8/)
    REAL(KIND=r8), PARAMETER ::  smref(9)=(/.283_r8,.387_r8,.412_r8,.312_r8,.338_r8,.382_r8,.315_r8,.329_r8,.283_r8/)
    REAL(KIND=r8), PARAMETER ::  smwlt(9)=(/.029_r8,.119_r8,.139_r8,.047_r8,.010_r8,.103_r8,.069_r8,.066_r8,.029_r8/)
   
    U10M=0.0_r8;    V10M=0.0_r8;    T2M=0.0_r8;    Q2M=0.0_r8
    CM=0.0_r8;    CH=0.0_r8;    RB=0.0_r8;    EVAP=0.0_r8
    EP=0.0_r8;    FM=0.0_r8;    FH=0.0_r8;    WIND=0.0_r8


!!!   save rsmax, rsmin, rgl, hs, smmax, smdry, smref, smwlt
    !
    DM=0.0_r8;    SNOWMT=0.0_r8;    GFLUX=0.0_r8;    F10M=0.0_r8
    RHSCNPY=0.0_r8;    RHSMC=0.0_r8;    AIM=0.0_r8;    BIM=0.0_r8
    CIM=0.0_r8;       DTV=0.0_r8; EC=0.0_r8
    EDIR=0.0_r8; ETPFAC=0.0_r8;    FH2=0.0_r8; FM10=0.0_r8
    FX=0.0_r8;    HL1=0.0_r8; HL12=0.0_r8; HLINF=0.0_r8;   PH=0.0_r8
    PH2=0.0_r8; PM=0.0_r8;      PM10=0.0_r8; PSURF=0.0_r8;   Q0=0.0_r8;      QS1=0.0_r8
    QSS=0.0_r8; RAT=0.0_r8;    RHO=0.0_r8; RS=0.0_r8
    THETA1=0.0_r8;  THV1=0.0_r8;    TSURF=0.0_r8;   TV1=0.0_r8
    TVS=0.0_r8;    XRCL=0.0_r8;Z0=0.0_r8;Z0MAX=0.0_r8;Z1=0.0_r8
    ZTMAX=0.0_r8;   PS1=0.0_r8;    OLINF=0.0_r8
    RESTAR=0.0_r8;    HL0=0.0_r8;    HL0INF=0.0_r8;    HL110=0.0_r8
    HLT=0.0_r8;    SIG2K=0.0_r8;    AA=0.0_r8;    BB0=0.0_r8
    CQ=0.0_r8;    HLTINF=0.0_r8;    FHS=0.0_r8;    FMS=0.0_r8
    BB=0.0_r8;    AA0=0.0_r8;    ADTV=0.0_r8

    !WRF      DELT2 = DELT * 2.
    !
    !     ESTIMATE SIGMA ** K AT 2 M
    !
    SIG2K = 1._r8 - 4._r8 * G * 2._r8 / (CP * 280._r8)
    !
    !  INITIALIZE VARIABLES. ALL UNITS ARE SUPPOSEDLY M.K.S. UNLESS SPECIFIE
    !  PSURF IS IN PASCALS
    !  WIND IS WIND SPEED, THETA1 IS ADIABATIC SURFACE TEMP FROM LEVEL 1
    !  RHO IS DENSITY, QS1 IS SAT. HUM. AT LEVEL1 AND QSS IS SAT. HUM. AT
    !  SURFACE
    !  CONVERT SLRAD TO THE CIVILIZED UNIT FROM LANGLEY MINUTE-1 K-4
    !  SURFACE ROUGHNESS LENGTH IS CONVERTED TO M FROM CM
    !
    !!
    !     qs1 = fpvs(t1)
    !     qss = fpvs(tskin)
    DO I=1,IM
       IF(mskant(i) == 1_r8) THEN
          XRCL(I)  = SQRT(RCL(I))
          PSURF(I) = 1000._r8 * PS(I)
          PS1(I)   = 1000._r8 * PRSL1(I)
          !       SLWD(I)  = SLRAD(I) * CONVRAD
          !WRF        SLWD(I)  = SLRAD(I)
          !
          !  DLWFLX has been given a negative sign for downward longwave
          !  snet is the net shortwave flux
          !
          !WRF        SNET(I) = -SLWD(I) - DLWFLX(I)
          WIND(I) = XRCL(I) * SQRT(U1(I) * U1(I) + V1(I) * V1(I)) + MAX(0.0_r8, MIN(DDVEL(I), 30.0_r8))
          WIND(I) = MAX(WIND(I),0.10_r8)
          Q0(I)   = MAX(Q1(I),1.E-8_r8)
          TSURF(I) = TSKIN(I)
          THETA1(I) = T1(I) * PRSLKI(I)
          TV1(I)  = T1(I)     * (1._r8 + RVRDM1 * Q0(I))
          THV1(I) = THETA1(I) * (1._r8 + RVRDM1 * Q0(I))
          TVS(I)  = TSURF(I) * (1._r8 + RVRDM1 * Q0(I))
          RHO(I) = PS1(I) / (RD * TV1(I))
          !jfe    QS1(I) = 1000._r8 * FPVS(T1(I))
          qs1(i) = fpvs(t1(i))
          QS1(I) = EPS * QS1(I) / (PS1(I) + EPSM1 * QS1(I))
          QS1(I) = MAX(QS1(I), 1.E-8_r8)
          Q0(I) = MIN(QS1(I),Q0(I))
          !jfe    QSS(I) = 1000. * FPVS(TSURF(I))
          qss(i) = fpvs(tskin(i))
          QSS(I) = EPS * QSS(I) / (PSURF(I) + EPSM1 * QSS(I))
          !       RS = PLANTR
          RS(I) = 0._r8
          !WRF        if(VEGTYPE(I).gt.0.) RS(I) = rsmin(VEGTYPE(I))
          Z0(I) = .01_r8 * Z0RL(i)
          !WRF        CANOPY(I)= MAX(CANOPY(I),0._r8)
          DM(I) = 1._r8
          !WRF
       END IF
    ENDDO
    !    rbyg  =rgas/grav*ztn*0.5_r8
    !
    !REAL(kind=r8), PARAMETER :: R       =2.8705e+2_r8 ! gas constant air    (J/kg/K)   
    !REAL(kind=r8), PARAMETER :: r_v          = 461.6_r8! gas constant H2O    (J/kg/K)
    !REAL(kind=r8), PARAMETER :: CP      =1.0046e+3_r8 ! spec heat air cp    (J/kg/K)
    !REAL(kind=r8),PARAMETER:: grav     =9.80665e+0_r8! gravity             (m/s2)


    DO I=1,IM       
       IF(mskant(i) == 1_r8) THEN
          !Z1(I) = -RD * TV1(I) * LOG(PS1(I)/PSURF(I)) / G          
          Z1(I) = ((RD/grav)*ztn(i)*0.5_r8)*TV1(I)
          DRAIN(I) = 0._r8
       END IF
    ENDDO

    !!
    DO K = 1, KM
       DO I=1,IM    
          IF(mskant(i) == 1_r8) THEN
             RHSMC(I,K) = 0._r8
             AIM(I,K) = 0._r8
             BIM(I,K) = 1._r8
             CIM(I,K) = 0._r8
          END IF
       ENDDO
    ENDDO

    DO I=1,IM       
       IF(mskant(i) == 1_r8) THEN
          EDIR(I) = 0._r8
          EC(I) = 0._r8
          EVAP(I) = 0._r8
          EP(I) = 0._r8
          SNOWMT(I) = 0._r8
          GFLUX(I) = 0._r8
          RHSCNPY(I) = 0._r8
          FX(I) = 0._r8
          ETPFAC(I) = 0._r8
          CANFAC(I) = 0._r8
       END IF
    ENDDO
    !
    !  COMPUTE STABILITY DEPENDENT EXCHANGE COEFFICIENTS
    !
    !  THIS PORTION OF THE CODE IS PRESENTLY SUPPRESSED
    !
    DO I=1,IM       
       IF(mskant(i) == 1_r8) THEN

          IF(SLIMSK(I).EQ.0._r8) THEN
             !OCEAN
             ! pkubota USTAR(I) = SQRT(G * Z0(I) / CHARNOCK)
             USTAR(I) = SQRT(G * Z0(I) / CHARNOCK)
          ENDIF
          !
          !  COMPUTE STABILITY INDICES (RB AND HLINF)
          !

          Z0MAX(I) = MIN(Z0(I),0.1_r8 * Z1(I))
          ZTMAX(I) = Z0MAX(I)
          IF(SLIMSK(I).EQ.0._r8) THEN
             !OCEAN
             RESTAR = USTAR(I) * Z0MAX(I) / VIS
             RESTAR = MAX(RESTAR,.000001_r8)
             !         RESTAR = ALOG(RESTAR)
             !         RESTAR = MIN(RESTAR,5.)
             !         RESTAR = MAX(RESTAR,-5.)
             !         RAT(I) = AA1 + BB1 * RESTAR + CC1 * RESTAR ** 2
             !         RAT(I) = RAT(I) / (1. + BB2 * RESTAR
             !    &                       + CC2 * RESTAR ** 2)
             !  Rat taken from Zeng, Zhao and Dickinson 1997
             RAT(I) = 2.67_r8 * restar ** .25_r8 - 2.57_r8
             RAT(I) = MIN(RAT(I),7._r8)
             ZTMAX(I) = Z0MAX(I) * EXP(-RAT(I))
          ENDIF
       END IF
    ENDDO

    DO I = 1, IM    
       IF(mskant(i) == 1_r8) THEN
          DTV(I) = THV1(I) - TVS(I)
          ADTV = ABS(DTV(I))
          ADTV = MAX(ADTV,.001_r8)
          DTV(I) = SIGN(1._r8,DTV(I)) * ADTV
          RB(I) = G * DTV(I) * Z1(I) / (.5_r8 * (THV1(I) + TVS(I))           &
               &          * WIND(I) * WIND(I))
          RB(I) = MAX(RB(I),-5000.0_r8)
          !        FM(I) = LOG((Z0MAX(I)+Z1(I)) / Z0MAX(I))
          !        FH(I) = LOG((ZTMAX(I)+Z1(I)) / ZTMAX(I))
          FM(I) = LOG((Z1(I)) / Z0MAX(I))
          FH(I) = LOG((Z1(I)) / ZTMAX(I))
          HLINF(I) = RB(I) * FM(I) * FM(I) / FH(I)
          FM10(I) = LOG((Z0MAX(I)+10._r8) / Z0MAX(I))
          FH2(I) = LOG((ZTMAX(I)+2._r8) / ZTMAX(I))
       END IF
    ENDDO
    !
    !  STABLE CASE
    !
    DO I = 1, IM    
       IF(mskant(i) == 1_r8) THEN
          IF(DTV(I).GE.0._r8) THEN
             HL1(I) = HLINF(I)
          ENDIF
          IF(DTV(I).GE.0._r8.AND.HLINF(I).GT..25_r8) THEN
             !PRINT*,   Z0MAX(I) ,ZTMAX(I), HLINF(I), Z1(I) , DTV(I)

             HL0INF = Z0MAX(I) * HLINF(I) / Z1(I)
             HLTINF = ZTMAX(I) * HLINF(I) / Z1(I)
             !IF(4.0_r8 * ALPHA * HLINF(I) < -1.0_r8)THEN
             !   AA     = SQRT(1._r8 + 4._r8 * ALPHA * ABS(HLINF(I))) 
             !ELSE
                AA     = SQRT(1._r8 + 4._r8 * ALPHA * HLINF(I))
             !END IF
             !IF( 4.0_r8 * ALPHA * HL0INF <-1.0_r8)THEN
             !   AA0    = SQRT(1._r8 + 4._r8 * ALPHA * ABS(HL0INF))
             !ELSE
                AA0    = SQRT(1._r8 + 4._r8 * ALPHA * HL0INF)
             !END IF 
             BB     = AA
             !IF(4.0_r8 * ALPHA * HLTINF < -1.0_r8)THEN 
             !   BB0 = SQRT(1.0_r8 + 4.0_r8 * ALPHA * ABS(HLTINF))
             !ELSE
                BB0 = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HLTINF)
             !END IF 
             PM(I)  = AA0 - AA + LOG((AA + 1._r8) / (AA0 + 1._r8))
             PH(I)  = BB0 - BB + LOG((BB + 1._r8) / (BB0 + 1._r8))
             FMS    = FM(I) - PM(I)
             FHS    = FH(I) - PH(I)
             HL1(I) = FMS * FMS * RB(I) / FHS
          END IF
       END IF
    ENDDO
    !
    !  SECOND ITERATION
    !
    DO I = 1, IM    
       IF(mskant(i) == 1_r8) THEN
          IF(DTV(I).GE.0._r8) THEN
             HL0 = Z0MAX(I) * HL1(I) / Z1(I)
             HLT = ZTMAX(I) * HL1(I) / Z1(I)
             !PRINT*,   HL1(I)
             !IF( 4.0_r8 * ALPHA * HL1(I) < -1.0_r8)THEN
             !    AA = SQRT(1.0_r8 + 4.0_r8 * ALPHA * ABS(HL1(I)))
             !ELSE          
                 AA = SQRT(1.0_r8 + 4.0_r8 * ALPHA * HL1(I))
             !END IF
             !IF(4._r8 * ALPHA * HL0 < -1.0_r8)THEN
             !   AA0 = SQRT(1._r8 + 4._r8 * ALPHA * ABS(HL0))
             !ELSE     
                AA0 = SQRT(1._r8 + 4._r8 * ALPHA * HL0)
             !END IF
             BB = AA
             !IF(4._r8 * ALPHA * HLT < -1.0_r8)THEN
             !   BB0 = SQRT(1._r8 + 4._r8 * ALPHA * ABS(HLT))
             !ELSE
                BB0 = SQRT(1._r8 + 4._r8 * ALPHA * HLT)
             !END IF
             PM(I) = AA0 - AA + LOG((AA + 1._r8) / (AA0 + 1._r8))
             PH(I) = BB0 - BB + LOG((BB + 1._r8) / (BB0 + 1._r8))
             HL110 = HL1(I) * 10._r8 / Z1(I)
             !IF(4._r8 * ALPHA * HL110 < -1.0_r8)THEN
             !   AA= SQRT(1._r8 + 4._r8 * ALPHA * ABS(HL110))
             !ELSE
                AA = SQRT(1._r8 + 4._r8 * ALPHA * HL110)
             !END IF
             PM10(I) = AA0 - AA + LOG((AA + 1.0_r8) / (AA0 + 1._r8))
             HL12(I) = HL1(I) * 2._r8 / Z1(I)
             !         AA = SQRT(1. + 4. * ALPHA * HL12(I))
             !IF(4._r8 * ALPHA * HL12(I) < -1.0_r8)THEN
             !   BB = SQRT(1._r8 + 4._r8 * ALPHA * ABS(HL12(I)))
             !ELSE
                BB = SQRT(1._r8 + 4._r8 * ALPHA * HL12(I))
             !END IF
             PH2(I) = BB0 - BB + LOG((BB + 1._r8) / (BB0 + 1._r8))
          ENDIF
       END IF
    ENDDO
    !
    !  UNSTABLE CASE
    !
    !
    !  CHECK FOR UNPHYSICAL OBUKHOV LENGTH
    !
    DO I=1,IM     
       IF(mskant(i) == 1_r8) THEN
          IF(DTV(I).LT.0._r8) THEN
             OLINF = Z1(I) / HLINF(I)
             IF(ABS(OLINF).LE.50._r8 * Z0MAX(I)) THEN
                HLINF(I) = -Z1(I) / (50._r8 * Z0MAX(I))
             ENDIF
          ENDIF
       END IF
    ENDDO
    !
    !  GET PM AND PH
    !
    DO I = 1, IM   
       IF(mskant(i) == 1_r8) THEN
          IF(DTV(I).LT.0._r8.AND.HLINF(I).GE.-.5_r8) THEN
             HL1(I) = HLINF(I)
             PM(I) = (A0 + A1 * HL1(I)) * HL1(I)                           &
                  &            / (1._r8 + B1 * HL1(I) + B2 * HL1(I) * HL1(I))
             PH(I) = (A0P + A1P * HL1(I)) * HL1(I)                         &
                  &            / (1._r8 + B1P * HL1(I) + B2P * HL1(I) * HL1(I))
             HL110 = HL1(I) * 10._r8 / Z1(I)
             PM10(I) = (A0 + A1 * HL110) * HL110                           &
                  &            / (1._r8 + B1 * HL110 + B2 * HL110 * HL110)
             HL12(I) = HL1(I) * 2._r8 / Z1(I)
             PH2(I) = (A0P + A1P * HL12(I)) * HL12(I)                      &
                  &            / (1._r8 + B1P * HL12(I) + B2P * HL12(I) * HL12(I))
          ENDIF
          IF(DTV(I).LT.0_r8.AND.HLINF(I).LT.-.5_r8) THEN
             HL1(I) = -HLINF(I)
             PM(I) = LOG(HL1(I)) + 2._r8 * HL1(I) ** (-.25_r8) - .8776_r8
             PH(I) = LOG(HL1(I)) + .5_r8 * HL1(I) ** (-.5_r8) + 1.386_r8
             HL110 = HL1(I) * 10._r8 / Z1(I)
             PM10(I) = LOG(HL110) + 2._r8 * HL110 ** (-.25_r8) - .8776_r8
             HL12(I) = HL1(I) * 2._r8 / Z1(I)
             PH2(I) = LOG(HL12(I)) + .5_r8 * HL12(I) ** (-.5_r8) + 1.386_r8
          ENDIF
       END IF
    ENDDO
    !
    !  FINISH THE EXCHANGE COEFFICIENT COMPUTATION TO PROVIDE FM AND FH
    !
    DO I = 1, IM
       IF(mskant(i) == 1_r8) THEN

          FM(I) = FM(I) - PM(I)
          FH(I) = FH(I) - PH(I)
          FM10(I) = FM10(I) - PM10(I)
          FH2(I) = FH2(I) - PH2(I)
          CM(I) = CONKE_CH*(CA * CA / (FM(I) * FM(I)))
          CH(I) = CONKE_CM*(CA * CA / (FM(I) * FH(I)))
          CQ = CH(I)
          STRESS(I) = CM(I) * WIND(I) * WIND(I)
          USTAR(I)  = SQRT(STRESS(I))
          !       USTAR(I) = SQRT(CM(I) * WIND(I) * WIND(I))
       END IF
    ENDDO
    !##DG  IF(LAT.EQ.LATD) THEN
    !##DG    PRINT *, ' FM, FH, CM, CH(I), USTAR =',
    !##DG &   FM, FH, CM, ch, USTAR
    !##DG  ENDIF
    !
    !  UPDATE Z0 OVER OCEAN
    !
    DO I = 1, IM      
       IF(mskant(i) == 1_r8) THEN
          IF(SLIMSK(I).EQ.0._r8) THEN
             !OCEAN
             Z0(I) = (CHARNOCK / G) * USTAR(I) ** 2
             !  NEW IMPLEMENTATION OF Z0
             !         CC = USTAR(I) * Z0 / RNU
             !         PP = CC / (1. + CC)
             !         FF = G * ARNU / (CHARNOCK * USTAR(I) ** 3)
             !         Z0 = ARNU / (USTAR(I) * FF ** PP)
             Z0(I) = MIN(Z0(I),.1_r8)
             Z0(I) = MAX(Z0(I),1.E-7_r8)
             Z0RL(I) = 100._r8 * Z0(I)
          ENDIF
       END IF
    ENDDO


    !
    !  CALCULATE SENSIBLE HEAT FLUX
    !
    !
    !  THE REST OF THE OUTPUT
    !
    !  CONVERT SNOW DEPTH BACK TO MM OF WATER EQUIVALENT
    !

    DO I = 1, IM       
       IF(mskant(i) == 1_r8) THEN
          F10M(I) = FM10(I) / FM(I)
          F10M(I) = MIN(F10M(I),1._r8)
          U10M(I) = F10M(I) * XRCL(I) * U1(I)
          V10M(I) = F10M(I) * XRCL(I) * V1(I)
          T2M(I) = 0
          Q2M(I)= 0
          !WRF         T2M(I) = TSKIN(I) * (1. - FH2(I) / FH(I))                      &
          !WRF     &          + THETA1(I) * FH2(I) / FH(I)
          !WRF         T2M(I) = T2M(I) * SIG2K
          !        Q2M(I) = QSURF(I) * (1. - FH2(I) / FH(I))                      &
          !    &         + Q1(I) * FH2(I) / FH(I)
          !       T2M(I) = T1
          !       Q2M(I) = Q1
          !WRF        IF(EVAP(I).GE.0.) THEN
          !
          !  IN CASE OF EVAPORATION, USE THE INFERRED QSURF TO DEDUCE Q2M
          !
          !WRF          Q2M(I) = QSURF(I) * (1. - FH2(I) / FH(I))                     &
          !WRF     &         + Q1(I) * FH2(I) / FH(I)
          !WRF        ELSE
          !
          !  FOR DEW FORMATION SITUATION, USE SATURATED Q AT TSKIN
          !
          !jfe      QSS(I) = 1000. * FPVS(TSKIN(I))
          !WRF          qss(I) = fpvs(tskin(I))
          !WRF          QSS(I) = EPS * QSS(I) / (PSURF(I) + EPSM1 * QSS(I))
          !WRF          Q2M(I) = QSS(I) * (1. - FH2(I) / FH(I))                       &
          !WRF     &         + Q1(I) * FH2(I) / FH(I)
          !WRF        ENDIF
          !jfe    QSS(I) = 1000. * FPVS(T2M(I))
          !WRF        QSS(I) = fpvs(t2m(I))
          !       QSS(I) = 1000. * T2MO(I)
          !WRF        QSS(I) = EPS * QSS(I) / (PSURF(I) + EPSM1 * QSS(I))
          !WRF        Q2M(I) = MIN(Q2M(I),QSS(I))
       END IF
    ENDDO

    RETURN
  END SUBROUTINE PROGTM


  FUNCTION KTSOIL(THETA,KTYPE)
    !

    IMPLICIT NONE

    INTEGER              ktype,kw
    REAL(kind=r8) ktsoil, theta, w
    !
    W = (THETA / TSAT(KTYPE)) * 20._r8 + 1._r8
    KW = W
    KW = MIN(KW,21)
    KW = MAX(KW,1)
    KTSOIL = DFKT(KW,KTYPE)                                           &
         &         + (W - KW) * (DFKT(KW+1,KTYPE) - DFKT(KW,KTYPE))
    RETURN
  END FUNCTION KTSOIL

  FUNCTION FUNCDF(THETA,KTYPE)
    !
    IMPLICIT NONE     
    INTEGER              ktype,kw
    REAL(kind=r8) funcdf,theta,w
    !
    W = (THETA / TSAT(KTYPE)) * 20._r8 + 1._r8
    KW = W
    KW = MIN(KW,21)
    KW = MAX(KW,1)
    FUNCDF = DFK(KW,KTYPE)                                            &
         &         + (W - KW) * (DFK(KW+1,KTYPE) - DFK(KW,KTYPE))
    RETURN
  END FUNCTION FUNCDF


  FUNCTION FUNCKT(THETA,KTYPE)
    !
    IMPLICIT NONE      
    INTEGER             ktype,kw
    REAL(kind=r8) funckt,theta,w
    !
    W = (THETA / TSAT(KTYPE)) * 20._r8 + 1._r8
    KW = W
    KW = MIN(KW,21)
    KW = MAX(KW,1)
    FUNCKT = KTK(KW,KTYPE)                                            &
         &         + (W - KW) * (KTK(KW+1,KTYPE) - KTK(KW,KTYPE))
    RETURN
  END FUNCTION FUNCKT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GRDDF
    IMPLICIT NONE
    INTEGER              i,    k
    REAL(kind=r8) dynw, f1, f2, theta
    !
    !  GRDDF SETS UP MOISTURE DIFFUSIVITY AND HYDROLIC CONDUCTIVITY
    !  FOR ALL SOIL TYPES
    !  GRDDFS SETS UP THERMAL DIFFUSIVITY FOR ALL SOIL TYPES
    !
    DO K = 1, NTYPE
       DYNW = TSAT(K) * .05_r8
       F1 = B(K) * SATKT(K) * SATPSI(K) / TSAT(K) ** (B(K) + 3._r8)
       F2 = SATKT(K) / TSAT(K) ** (B(K) * 2._r8 + 3._r8)
       !
       !  CONVERT FROM M/S TO KG M-2 S-1 UNIT
       !
       F1 = F1 * 1000._r8
       F2 = F2 * 1000._r8
       DO I = 1, NGRID
          THETA = real(I-1,kind=r8) * DYNW
          THETA = MIN(TSAT(K),THETA)
          DFK(I,K) = F1 * THETA ** (B(K) + 2._r8)
          KTK(I,K) = F2 * THETA ** (B(K) * 2._r8 + 3._r8)
       ENDDO
    ENDDO
  END SUBROUTINE GRDDF


  SUBROUTINE GRDKT()
    IMPLICIT NONE
    INTEGER              i,    k
    REAL(kind=r8) dynw, f1, theta, pf
    DO K = 1, NTYPE
       DYNW = TSAT(K) * .05_r8
       F1 = LOG10(SATPSI(K)) + B(K) * LOG10(TSAT(K)) + 2._r8
       DO I = 1, NGRID
          THETA = real(I-1,kind=r8) * DYNW
          THETA = MIN(TSAT(K),THETA)
          IF(THETA.GT.0._r8) THEN
             PF = F1 - B(K) * LOG10(THETA)
          ELSE
             PF = 5.2_r8
          ENDIF
          IF(PF.LE.5.1_r8) THEN
             DFKT(I,K) = EXP(-(2.7_r8+PF)) * 420._r8
          ELSE
             DFKT(I,K) = .1744_r8
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE GRDKT









 !-------------------------------------------------------------------

  !----------------------------------------------------------------
  SUBROUTINE OCEANML( &
       ust      , & !REAL,     INTENT(IN   ) :: UST  ( 1:nCols )
       u_phy    , & !REAL,     INTENT(IN   ) :: U_PHY( 1:nCols )
       v_phy    , & !REAL,     INTENT(IN   ) :: V_PHY( 1:nCols )
       mskant    , & !REAL,     INTENT(IN   ) :: mskant( 1:nCols )
       HFX      , & !REAL,     INTENT(IN   ) :: HFX  ( 1:nCols )
       LH       , & !REAL,     INTENT(IN   ) :: LH   ( 1:nCols )
       tsea     , & !REAL,     INTENT(IN   ) :: tsea ( 1:nCols )
       TSK      , & !REAL,     INTENT(INOUT) :: TSK  ( 1:nCols )
       HML      , & !REAL,     INTENT(INOUT) :: HML  ( 1:nCols )
       HUML     , & !REAL,     INTENT(INOUT) :: HUML ( 1:nCols )
       HVML     , & !REAL,     INTENT(INOUT) :: HVML ( 1:nCols )
       GSW      , & !REAL,     INTENT(IN   ) :: GSW  ( 1:nCols )
       GLW      , & !REAL,     INTENT(IN   ) :: GLW  ( 1:nCols )
       EMISS    , & !REAL,     INTENT(IN   ) :: EMISS( 1:nCols )
       DT       , & !REAL,     INTENT(IN   ) :: DT
       H0ML , &
       nCols      ) !INTEGER,  INTENT(IN   ) :: nCols

    !----------------------------------------------------------------
    IMPLICIT NONE
    !----------------------------------------------------------------
    !
    !  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK)
    !  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON
    !  (Pollard, Rhines and Thompson (1973).
    !
    !-- TML         ocean mixed layer temperature (K)
    !-- T0ML        ocean mixed layer temperature (K) at initial time
    !-- TMOML       top 200 m ocean mean temperature (K) at initial time
    !-- HML         ocean mixed layer depth (m)
    !-- H0ML        ocean mixed layer depth (m) at initial time
    !-- HUML        ocean mixed layer u component of wind
    !-- HVML        ocean mixed layer v component of wind
    !-- OML_GAMMA   deep water lapse rate (K m-1)
    !-- UAIR,VAIR   lowest model level wind component
    !-- UST         frictional velocity
    !-- HFX         upward heat flux at the surface (W/m^2)
    !-- LH          latent heat flux at the surface (W/m^2)
    !-- TSK         surface temperature (K)
    !-- GSW         downward short wave flux at ground surface (W/m^2)
    !-- GLW         downward long wave flux at ground surface (W/m^2)
    !-- EMISS       emissivity of the surface
    !-- XLAND       land mask (1 for land, 2 for water)
    !-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
    !-- F           Coriolis parameter      fcor   =  5.0e-5        ! Coriolis parameter (1/s)
    !-- DT          time step (second)
    !-- G           acceleration due to gravity
    INTEGER,  INTENT(IN   ) :: nCols
    REAL(KIND=r8),     INTENT(IN   ) :: H0ML( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: DT
    REAL(KIND=r8),     INTENT(IN   ) :: EMISS( 1:nCols )
    INTEGER(KIND=I8),  INTENT(IN   ) :: mskant( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: GSW  ( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: GLW  ( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: HFX  ( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: LH   ( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: tsea ( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: U_PHY( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: V_PHY( 1:nCols )
    REAL(KIND=r8),     INTENT(IN   ) :: UST  ( 1:nCols )
    REAL(KIND=r8),     INTENT(INOUT) :: TSK  ( 1:nCols )
    REAL(KIND=r8),     INTENT(INOUT) :: HML  ( 1:nCols )
    REAL(KIND=r8),     INTENT(INOUT) :: HUML ( 1:nCols )
    REAL(KIND=r8),     INTENT(INOUT) :: HVML ( 1:nCols )

    REAL(KIND=r8)  :: F    ( 1:nCols )
    REAL(KIND=r8),     PARAMETER     :: G        =   9.81_r8
    REAL(KIND=r8),     PARAMETER     :: OML_GAMMA=  0.14_r8 !-- oml_gamma     lapse rate below mixed layer in ocean (default 0.14 K m-1)
    REAL(KIND=r8)    , PARAMETER     ::  STBOLT=5.67051E-8_r8
    ! LOCAL VARS
    REAL(KIND=r8) :: T0ML (nCols)
    REAL(KIND=r8) :: TMOML(nCols)
    REAL(KIND=r8) :: TML  (nCols)

    INTEGER ::  I
    DO i=1,nCols
       TML  (i)=tsea(I)
       T0ML (I)=tsea(I)
       TMOML(I)=tsea(I)-3.0_r8
       F    ( i )= 5.0e-5_r8   
    END DO
    CALL OML1D(  &
         TML  ( 1:nCols ), & !REAL, INTENT(INOUT)  :: TML  ( 1:nCols )
         T0ML ( 1:nCols ), & !REAL, INTENT(IN   )  :: T0ML ( 1:nCols )
         HML  ( 1:nCols ), & !REAL, INTENT(INOUT)  :: H    ( 1:nCols )
         H0ML ( 1:nCols ), & !REAL, INTENT(IN   )  :: H0   ( 1:nCols )
         HUML ( 1:nCols ), & !REAL, INTENT(INOUT)  :: HUML ( 1:nCols )
         HVML ( 1:nCols ), & !REAL, INTENT(INOUT)  :: HVML ( 1:nCols )
         TSK  ( 1:nCols ), & !REAL, INTENT(INOUT)  :: TSK  ( 1:nCols )
         HFX  ( 1:nCols ), & !REAL, INTENT(IN   )  :: HFX  ( 1:nCols )
         LH   ( 1:nCols ), & !REAL, INTENT(IN   )  :: LH   ( 1:nCols )
         GSW  ( 1:nCols ), & !REAL, INTENT(IN   )  :: GSW  ( 1:nCols )
         GLW  ( 1:nCols ), & !REAL, INTENT(IN   )  :: GLW  ( 1:nCols )
         TMOML( 1:nCols ), & !REAL, INTENT(IN   )  :: TMOML( 1:nCols )
         U_PHY( 1:nCols ), & !REAL, INTENT(IN   )  :: UAIR ( 1:nCols )
         V_PHY( 1:nCols ), & !REAL, INTENT(IN   )  :: VAIR ( 1:nCols )
         UST  ( 1:nCols ), & !REAL, INTENT(IN   )  :: UST  ( 1:nCols )
         F    ( 1:nCols ), & !REAL, INTENT(IN   )  :: F    ( 1:nCols )
         EMISS( 1:nCols ), & !REAL, INTENT(IN   )  :: EMISS( 1:nCols )
         mskant( 1:nCols ), & !REAL, INTENT(IN   )  :: mskant( 1:nCols )
         STBOLT          , & !REAL, INTENT(IN   )  :: STBOLT
         G               , & !REAL, INTENT(IN   )  :: G
         DT              , & !REAL, INTENT(IN   )  :: DT
         OML_GAMMA       , & !REAL, INTENT(IN   )  :: OML_GAMMA
         nCols             ) !INTEGER, INTENT(IN   )  :: nCols

  END SUBROUTINE OCEANML

  !----------------------------------------------------------------
  SUBROUTINE OML1D(TML      , & !REAL,    INTENT(INOUT)    :: TML  ( 1:nCols )
       T0ML     , & !REAL,    INTENT(IN   )    :: T0ML ( 1:nCols )
       H        , & !REAL,    INTENT(INOUT)    :: H    ( 1:nCols )
       H0       , & !REAL,    INTENT(IN   )    :: H0   ( 1:nCols )
       HUML     , & !REAL,    INTENT(INOUT)    :: HUML ( 1:nCols )
       HVML     , & !REAL,    INTENT(INOUT)    :: HVML ( 1:nCols )
       TSK      , & !REAL,    INTENT(INOUT)    :: TSK  ( 1:nCols )
       HFX      , & !REAL,    INTENT(IN   )    :: HFX  ( 1:nCols )
       LH       , & !REAL,    INTENT(IN   )    :: LH   ( 1:nCols )
       GSW      , & !REAL,    INTENT(IN   )    :: GSW  ( 1:nCols )
       GLW      , & !REAL,    INTENT(IN   )    :: GLW  ( 1:nCols )
       TMOML    , & !REAL,    INTENT(IN   )    :: TMOML( 1:nCols )
       UAIR     , & !REAL,    INTENT(IN   )    :: UAIR ( 1:nCols )
       VAIR     , & !REAL,    INTENT(IN   )    :: VAIR ( 1:nCols )
       UST      , & !REAL,    INTENT(IN   )    :: UST  ( 1:nCols )
       F        , & !REAL,    INTENT(IN   )    :: F    ( 1:nCols )
       EMISS    , & !REAL,    INTENT(IN   )    :: EMISS( 1:nCols )
       mskant    , & !REAL,    INTENT(IN   )    :: mskant( 1:nCols )
       STBOLT   , & !REAL,    INTENT(IN   )    :: STBOLT
       G        , & !REAL,    INTENT(IN   )    :: G
       DT       , & !REAL,    INTENT(IN   )    :: DT
       OML_GAMMA, & !REAL,    INTENT(IN   )    :: OML_GAMMA
       nCols      ) !INTEGER, INTENT(IN   )    :: nCols

    !----------------------------------------------------------------
    IMPLICIT NONE
    !----------------------------------------------------------------
    !
    !  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK) 
    !  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON 
    !  (Pollard, Rhines and Thompson (1973).
    !
    !-- TML         ocean mixed layer temperature (K)
    !-- T0ML        ocean mixed layer temperature (K) at initial time
    !-- TMOML       top 200 m ocean mean temperature (K) at initial time
    !-- H           ocean mixed layer depth (m)
    !-- H0          ocean mixed layer depth (m) at initial time
    !-- HUML        ocean mixed layer u component of wind
    !-- HVML        ocean mixed layer v component of wind
    !-- OML_GAMMA   deep water lapse rate (K m-1)
    !-- OMLCALL     whether to call oml model
    !-- UAIR,VAIR   lowest model level wind component
    !-- UST         frictional velocity
    !-- HFX         upward heat flux at the surface (W/m^2)
    !-- LH          latent heat flux at the surface (W/m^2)
    !-- TSK         surface temperature (K)
    !-- GSW         downward short wave flux at ground surface (W/m^2)
    !-- GLW         downward long wave flux at ground surface (W/m^2)
    !-- EMISS       emissivity of the surface
    !-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
    !-- F           Coriolis parameter
    !-- DT          time step (second)
    !-- G           acceleration due to gravity
    !
    !----------------------------------------------------------------
    INTEGER, INTENT(IN   )    :: nCols
    REAL(KIND=r8),    INTENT(INOUT)    :: TML ( 1:nCols )
    REAL(KIND=r8),    INTENT(INOUT)    :: H   ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: H0  ( 1:nCols )
    REAL(KIND=r8),    INTENT(INOUT)    :: HUML( 1:nCols )
    REAL(KIND=r8),    INTENT(INOUT)    :: HVML( 1:nCols )
    REAL(KIND=r8),    INTENT(INOUT)    :: TSK ( 1:nCols )

    REAL(KIND=r8),    INTENT(IN   )    :: T0ML ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: TMOML( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: HFX  ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: LH   ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: GSW  ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: GLW  ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: UAIR ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: VAIR ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: UST  ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: F    ( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: EMISS( 1:nCols )
    INTEGER(KIND=I8),    INTENT(IN   )    :: mskant( 1:nCols )
    REAL(KIND=r8),    INTENT(IN   )    :: STBOLT, G, DT, OML_GAMMA

    ! Local
    REAL(KIND=r8) :: rhoair( 1:nCols ) 
    REAL(KIND=r8) :: rhowater( 1:nCols ) 
    REAL(KIND=r8) :: Gam( 1:nCols ) 
    REAL(KIND=r8) :: alp 
    REAL(KIND=r8) :: BV2 
    REAL(KIND=r8) :: A1 
    REAL(KIND=r8) :: A2 
    REAL(KIND=r8) ::  B2 
    REAL(KIND=r8) :: u 
    REAL(KIND=r8) :: v 
    REAL(KIND=r8) :: wspd
    REAL(KIND=r8) :: hu1 
    REAL(KIND=r8) :: hv1 
    REAL(KIND=r8) :: hu2 
    REAL(KIND=r8) :: hv2 
    REAL(KIND=r8) :: taux 
    REAL(KIND=r8) :: tauy 
    REAL(KIND=r8) :: tauxair 
    REAL(KIND=r8) :: tauyair 
    REAL(KIND=r8) :: q( 1:nCols ) 
    REAL(KIND=r8) :: hold( 1:nCols )
    REAL(KIND=r8) :: hsqrd 
    REAL(KIND=r8) :: thp 
    REAL(KIND=r8) ::  cwater( 1:nCols ) 
    REAL(KIND=r8) :: ust2
    !CHARACTER(LEN=120) :: time_series
    INTEGER :: i

    DO i=1,nCols
       IF (mskant(i) == 1_i8) THEN

          hu1=huml(i)
          hv1=hvml(i)
          rhoair(i)=1.0_r8
          rhowater(i)=1000.0_r8
          cwater(i)=4200.0_r8
          ! Deep ocean lapse rate (K/m) - from Rich
          Gam(i)=oml_gamma
          !     if(i.eq.1 .eq.1 .or. i.eq.105.eq.105) print *, 'gamma = ', gam
          !     Gam=0.14
          !     Gam=5.6/40.
          !     Gam=5./100.
          ! Thermal expansion coeff (/K)
          !     alp=.0002
          !     temp dependence (/K)
          alp=MAX((tml(i)-273.15_r8)*1.e-5_r8, 1.e-6_r8)
          BV2=alp*g*Gam(i)
          thp=t0ml(i)-Gam(i)*(h(i)-h0(i))
          A1=(tml(i)-thp)*h(i) - 0.5_r8*Gam(i)*h(i)*h(i)
          IF(h(i).NE.0.0_r8)THEN
             u=hu1/h(i)
             v=hv1/h(i)
          ELSE
             u=0.0_r8
             v=0.0_r8
          ENDIF

          !  time step

          q(i)=(-hfx(i)-lh(i)+gsw(i)+glw(i)-stbolt*emiss(i)*tml(i)*tml(i)*tml(i)*tml(i))/(rhowater(i)*cwater(i))
          !       wspd=max(sqrt(uair*uair+vair*vair),0.1)
          wspd=SQRT(uair(i)*uair(i)+vair(i)*vair(i))
          IF (wspd .LT. 1.e-10_r8 ) THEN
             !          print *, 'i,wspd are ', i,wspd
             wspd = 1.e-10_r8
          ENDIF
          ! limit ust to 1.6 to give a value of ust for water of 0.05
          !       ust2=min(ust, 1.6)
          ! new limit for ust: reduce atmospheric ust by half for ocean
          ust2=0.5_r8*ust(i)
          tauxair=ust2*ust2*uair(i)/wspd
          taux=rhoair(i)/rhowater(i)*tauxair
          tauyair=ust2*ust2*vair(i)/wspd
          tauy=rhoair(i)/rhowater(i)*tauyair
          ! note: forward-backward coriolis force for effective time-centering
          hu2=hu1+dt*( f(i)*hv1 + taux)
          hv2=hv1+dt*(-f(i)*hu2 + tauy)
          ! consider the flux effect
          A2=A1+q(i)*dt

          huml(i)=hu2
          hvml(i)=hv2

          hold(i)=h(i)
          B2=hu2*hu2+hv2*hv2
          hsqrd=-A2/Gam(i) + SQRT(A2*A2/(Gam(i)*Gam(i)) + 2.0_r8*B2/BV2)
          h(i)=SQRT(MAX(hsqrd,0.0_r8))
          ! limit to positive h change
          IF(h(i).LT.hold(i))h(i)=hold(i)

          ! no change unless tml is warmer than layer mean temp tmol or tsk-5 (see omlinit)
          IF(tml(i).GE.tmoml(i) .AND. h(i).NE.0.0_r8)THEN
             tml(i)=MAX(t0ml(i) - Gam(i)*(h(i)-h0(i)) + 0.5_r8*Gam(i)*h(i) + A2/h(i), tmoml(i))
             u=hu2/h(i)
             v=hv2/h(i)
          ELSE
             tml(i)=t0ml(i)
             u=0.0_r8
             v=0.0_r8
          ENDIF
          tsk(i)=tml(i)

       END IF
    ENDDO
  END SUBROUTINE OML1D




END MODULE Sfc_SeaFlux_WGFS_Model




