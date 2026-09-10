!
!  $Author: pkubota $
!  $Date: 2010/04/20 20:18:04 $
!  $Revision: 1.18 $
!
MODULE GridDynamics

  USE Constants,   ONLY : &
       tref,              & ! intent(in)
       psref,             & ! intent(in)
       gasr,              & ! intent(in)
       p000,              & ! intent(in)
       cp,                & ! intent(in)
       cpv,               & ! intent(in)
       qpsm1,             & ! intent(in)
       rk,                & ! intent(in)
       qmin,              & ! intent(in)
       i8,                & ! intent(in)
       r8,                & ! intent(in)
       grav 

  USE FieldsDynamics, ONLY : &
       fgtmp  ,  & ! intent(in)
!       fgtmpC ,  & ! intent(inout)
       fgdivm ,  & ! intent(inout)
       fgdiv  ,  & ! intent(in)
!       fgdivC ,  & ! intent(inout)
       fgum   ,  & ! intent(inout)
       fgu    ,  & ! intent(in)
!       fguC   ,  & ! intent(inout)
       fgvm   ,  & ! intent(inout)
       fgv    ,  & ! intent(in)
!       fgvC   ,  & ! intent(inout)
       fgw    ,  & ! intent(in)
       fgwm   ,  & ! intent(inout)
       fgqm   ,  & ! intent(inout)
       fgqmm  ,  & ! intent(out)
       fgq    ,  & ! intent(inout)
       fgqp   ,  & ! intent(inout)
       fgice  ,  & !
       fgicem ,  & !
       fgicep ,  & !
       fgliq  ,  & !
       fgliqm ,  & !
       fgliqp ,  & !
       fgvar  ,  & !
       fgvarm ,  & !
       fgvarp ,  & !
       fgtlamm,  & ! intent(inout)
       fgtlam ,  & ! intent(in)
!       fgtlamC,  & ! intent(inout)
       fgtphim,  & ! intent(inout)
       fgtphi ,  & ! intent(in)
!      fgtphiC,  & ! intent(inout)
       fglnpm ,  & ! intent(inout)
       fglnps ,  & ! intent(in)
       fglnps_nabla4,  & ! intent(in)
!       fglnpsC,  & ! intent(inout)
       fgplamm,  & ! intent(inout)
       fgplam ,  & ! intent(in)
!       fgplamC,  & ! intent(inout)
       fgpphim,  & ! intent(inout)
       fgpphi,   & ! intent(in)
!       fgpphiC,  & ! intent(inout)
       fgyum,    & ! intent(in)
       fgyvm,    & ! intent(in)
       fgtdm,    & ! intent(in)
       fgvdlnpm, & ! intent(in)
       fgtmpm,   & ! intent(inout)
       fgyu,     & ! intent(inout)
       fgyv,     & ! intent(inout)
       fgtd,     & ! intent(inout)
       fgqd,     & ! intent(inout)
       fgvdlnp,  & ! intent(inout)
       fgps,     & ! intent(in)
       fgzs,     & ! intent(in)
       fgpass_scalars,  &
       adr_scalars,     &
       fgpass_fluxscalars, &
       fgpsp

  USE Sizes, ONLY:   &
       ibMaxPerJB,   &
       iPerIJB,      &
       jPerIJB,      &
       ibPerIJ,      & 
       jbPerIJ,      &
       imaxperj,     &
       myfirstlat_diag, &
       mylastlat_diag,  &
       myfirstlon,   &
       mylastlon,    &
       myjmax_d,     &
       jmax,         &
       imax,         &
       nlatsinproc_d,&
       ibMax,        & ! intent(in)
       jbMax,        & ! intent(in)
       kMax,         &
       a_hybr,       &
       b_hybr,       &
       c_hybr,       &
       delb,         & ! intent(in)
       p_r,          &
       delp_r,       &
       alpha_r,      &
       rpi_r,        &
       del_eta,      & ! intent(in)
       eta 




  USE PhysicsDriver, ONLY : &
       DryPhysics,          &
       SimpPhys


  USE GridHistory, ONLY:  &
       StoreGridHistory , &
       IsGridHistoryOn  , &
       dogrh            , &
       nGHis_temper    , &
       nGHis_uzonal    , &
       nGHis_vmerid    , &
       nGHis_spchum    , &
       nGHis_tvirsf    , &
       nGHis_uzonsf    , &
       nGHis_vmersf    , &
       nGHis_sphusf    , &
       nGHis_snowdp    , &
       nGHis_rouglg    , &
       nGHis_tcanop    , &
       nGHis_tgfccv    , &
       nGHis_tgdeep    , &
       nGHis_swtsfz    , &
       nGHis_swtrtz    , &
       nGHis_swtrcz    , &
       nGHis_mostca    , &
       nGHis_mostgc    , &
       nGHis_vegtyp    , &
       nGHis_presfc    , &
       nGHis_tspres

  USE Diagnostics, ONLY:   &
       StartStorDiag     , &
       pwater            , &
       updia             , &
       dodia             , &
       nDiag_tmpsfc      , & ! time mean surface pressure
       nDiag_tmtsfc      , & ! time mean surface temperature
       nDiag_omegav      , & ! omega
       nDiag_sigdot      , & ! sigma dot
       nDiag_pwater      , & ! precipitable water
       nDiag_divgxq      , & ! divergence * specific humidity
       nDiag_vmoadv      , & ! vertical moisture advection
       nDiag_tmtdps      , & ! time mean deep soil temperature
       nDiag_tgfccv      , & ! ground/surface cover temperature
       nDiag_tcanop      , & ! canopy temperature
       nDiag_homtvu      , & ! Horizontal Momentum Transport
       nDiag_vzmtwu      , & ! Vertical Zonal Momentum Transport
       nDiag_vmmtwv      , & ! Vertical Meridional Momentum Transport
       nDiag_mshtvt      , & ! Meridional Sensible Heat Transport
       nDiag_zshtut      , & ! Zonal Sensible Heat Transport
       nDiag_vshtwt      , & ! Vertical Sensible Heat Transport
       nDiag_mshtuq      , & ! Meridional Specific Humidity Transport
       nDiag_zshtuq      , & ! Zonal Specific Humidity Transport
       nDiag_vshtwq      , & ! Vertical Specific Humidity Transport
       nDiag_dewptt      , & ! Dew Point Temperature K
       nDiag_tspres          ! TIME MEAN MAXIMUM TENDENCY SFC PRESSURE (Pa)

  USE FieldsPhysics, ONLY: &
       dump              , &
       zorl              , &
       sheleg            , &
       imask             , &
       gtsea             , &
       tcm               , &
       tgm               , &
       tdm               , &
       wm                , &
       capacm

  USE Options, ONLY:       &
       alfa              , &!Nilo new fil
       isimp             , &
       eps_sic           , &
       nscalars          , &
       nClass            , &
       nAeros            , &
       microphys         , &
       SL_twotime_scheme , &
       istrt

  USE Utils  , ONLY:       &
       coslat            , &
       sinlat            , &
       longit            , &
       cel_area          , &
       massconsrv        , &
       fconsrv           , &
       fconsrv_flux      , &
       totmas            , &
       totflux           , &
       total_mass        , &
       total_flux

  USE Parallelism, ONLY:       &
       MsgOne,                 &
       maxnodes

  USE Communications, ONLY:    &
       Collect_Gauss

   IMPLICIT NONE
  SAVE       


  INCLUDE 'mpif.h'

  PRIVATE
  PUBLIC :: AddTend
  PUBLIC :: GrpComp
  PUBLIC :: TimeFilterStep1
  PUBLIC :: TimeFilterStep2
  PUBLIC :: GlobConservation
  PUBLIC :: GlobFluxConservation
  PUBLIC :: UpdateConserv
  PUBLIC :: End_temp_difus
! PUBLIC :: Scalardiffusion
  PUBLIC :: SetJablo

  REAL(KIND=r8) :: totmass
  REAL(KIND=r8), ALLOCATABLE :: fg (:,:,:)
  REAL(KIND=r8), ALLOCATABLE :: fgs(:,:,:)
  REAL(KIND=r8), ALLOCATABLE :: fg_flux (:,:,:)
  REAL(KIND=r8), ALLOCATABLE :: fgs_flux(:,:,:)

  REAL(KIND=r8), ALLOCATABLE :: fps(:,:)
  INTEGER      , ALLOCATABLE :: displ(:)
  INTEGER      , ALLOCATABLE :: displ_flux(:)

  LOGICAL      , PUBLIC      :: do_globconserv
  LOGICAL      , PUBLIC      :: do_globfluxconserv
  LOGICAL      , PUBLIC      :: init_globconserv
  LOGICAL      , PUBLIC      :: init_globfluxconserv
  INTEGER :: ierr
  CHARACTER(LEN=*), PARAMETER :: h="**(grpcomp)**"

CONTAINS


  SUBROUTINE GrpComp(&
       gyu    , gyv    , gtd    , gqd    , &
       gvdlnp , gdiv   , gtmp   , grot   , &
       gu     , gv     , gw     , gq     , &
       gplam  , gpphi  , gum    , gzs    , &
       gvm    , gtm    , gqm    , omg    , &
       ps     , gtlam  , gtphi  , gqlam  , &
       gqphi  , gulam  , guphi  , gvlam  , &
       gvphi  , gtlamm , gtphim , gplamm , &
       gpphim , glnpm  , gdivm  , gzslam , &
       gzsphi , gyum   , gyvm   , gtdm   , &
       gqdm   , gvdlnpm, colrad , rcl    , &
       vmax   , ifday  , tod    ,          &
       ibMax  , kMax   , ibLim  , slagr  , &
       slhum  , jb     , lonrad , cos2d  , &
       intcosz, cos2lat, ercossin, fcor  , &
       cosiv  , initial, gprsl  ,  gprsi , &
       gphil  , gphii  , &
       gicem  , gice   ,gicet  , &
       gliqm  , gliq   ,gliqt  , &
       gvarm  ,gvar  , gvart)

    !
    ! grpcomp: grid-point computations (all tendencies are computed) 
    !
    !
    ! slagr is the option for eulerian (slagr=.false.) or
    ! semi-Lagrangian integration (slagr=.true.)
    !
    !
    !
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: initial
    REAL(KIND=r8),    INTENT(OUT  ) :: gyu    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gyv    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gtd    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gqd    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gvdlnp (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gdiv   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gtmp   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: grot   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gu     (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gv     (ibMax, kMax)
    REAL(KIND=r8),    INTENT(OUT  ) :: gw     (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gq     (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gplam  (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gpphi  (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gzs    (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gum    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gvm    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gtm    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gqm    (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: omg    (ibMax, kMax) 
    REAL(KIND=r8),    INTENT(IN   ) :: ps     (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gtlam  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gtphi  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gqlam  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gqphi  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gulam  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: guphi  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gvlam  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gvphi  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gtlamm (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gtphim (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gplamm (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gpphim (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: glnpm  (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gdivm  (ibMax, kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gzslam (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gzsphi (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gyum   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gyvm   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gtdm   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gqdm   (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gvdlnpm(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: colrad (ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: rcl    (ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: vmax   (kMax)

    REAL(KIND=r8),    INTENT(INOUT) :: gprsl (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gprsi (ibMax, kMax+1)
    REAL(KIND=r8),    INTENT(INOUT) :: gphil (ibMax, kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: gphii (ibMax, kMax+1)

    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gicem  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gice   (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gicet  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliqm  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliq   (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gliqt  (ibMax, kMax)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gvarm  (ibMax, kMax,nClass+nAeros)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gvar   (ibMax, kMax,nClass+nAeros)
    REAL(KIND=r8),    OPTIONAL,   INTENT(INOUT) :: gvart  (ibMax, kMax,nClass+nAeros)
     
    INTEGER, INTENT(IN   ) :: ifday
    REAL(KIND=r8),    INTENT(IN   ) :: tod
    INTEGER, INTENT(IN   ) :: ibLim
    LOGICAL, INTENT(IN   ) :: slagr
    LOGICAL, INTENT(IN   ) :: slhum
    INTEGER, INTENT(IN   ) :: jb
    REAL(KIND=r8)   , INTENT(IN   ) :: lonrad (ibMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cos2d  (ibMax)    
    REAL(KIND=r8)   , INTENT(IN   ) :: cos2lat(ibMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ercossin(ibMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: fcor   (ibMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: cosiv  (ibMax)
    LOGICAL, INTENT(IN   ) :: intcosz
    !
    !  local variables 
    !
    REAL(KIND=r8)   , DIMENSION(ibMax,max(kMax,9)) :: gdt
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: psint
    REAL(KIND=r8)   , DIMENSION(ibMax     ) :: zsint
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: adveps
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: divint
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: delp 
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: rpi  
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: alpha
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: beta
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax+1) :: eta_dp_eta
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: dewpwtr     ! used for dewpoint and pwater
    REAL(KIND=r8)   , DIMENSION(ibMax,kMax) :: dewpoint
    REAL(KIND=r8)   , DIMENSION(ibMax) :: tendpress
    INTEGER                        :: i
    INTEGER                        :: k
    INTEGER                        :: ki
    INTEGER                        :: ncount
    INTEGER                        :: latco
    tendpress=0.0_r8
    gyu=0.0_r8
    gyv=0.0_r8
    gtd=0.0_r8
    gqd=0.0_r8
    gvdlnp=0.0_r8
    gdt=0.0_r8
    IF (microphys) THEN
       gicet =0.0_r8
       gliqt=0.0_r8
       IF((nClass+nAeros) >0 .and. PRESENT(gvarm))THEN
           gvart  =0.0_r8
       END IF
    END IF

    IF (initial.eq.1.or..not.slagr.or..not.SL_twotime_scheme) THEN
       gyum=0.0_r8
       gyvm=0.0_r8
       gtdm=0.0_r8
       gqdm=0.0_r8
       gvdlnpm=0.0_r8
    END IF
    latco=jb
    ! 
    !
    ! enforce humidity to be above a certain level (avoid negative values...)
    ! ----------------------------------------------------------------------
    !
    IF(TRIM(isimp).ne.'YES')THEN
       gq=MAX(gq,qmin)
    ENDIF
    !
    IF(dodia(nDiag_tmtdps).or.dodia(nDiag_tgfccv).or.dodia(nDiag_tcanop).or. &
       dodia(nDiag_tmtsfc))THEN
       ncount=0
       DO i=1,ibLim
         IF(imask(i,jb).ge.1_i8)THEN
             ncount=ncount+1
             gdt(i,1)=tcm(ncount,jb)
             gdt(i,2)=tgm(ncount,jb)
             gdt(i,3)=tdm(ncount,jb)
         ELSE 
             gdt(i,1)=ABS(gtsea (i,jb))
             gdt(i,2)=ABS(gtsea (i,jb))
             gdt(i,3)=ABS(gtsea (i,jb))
         END IF
       END DO
       IF(dodia(nDiag_tmtdps))CALL updia(gdt(1:ibLim,3:3),nDiag_tmtdps,latco)
       IF(dodia(nDiag_tgfccv))CALL updia(gdt(1:ibLim,2:2),nDiag_tgfccv,latco)
       IF(dodia(nDiag_tcanop))CALL updia(gdt(1:ibLim,1:1),nDiag_tcanop,latco)
       IF(dodia(nDiag_tmtsfc))CALL updia(ABS(gtsea(1:ibLim,jb)),nDiag_tmtsfc,latco)
    END IF
    !
    !     obtain grid history fields if requested
    !
    IF(IsGridHistoryOn())THEN





       ! invert fields for Diagnostics (from top to bottom to bottom to top )
       IF(dogrh(nGHis_temper,latco)) CALL StoreGridHistory(gtmp(1:ibLim,kMax:1:-1),nGHis_temper,latco)
       IF(dogrh(nGHis_uzonal,latco)) CALL StoreGridHistory( gu(1:ibLim,kMax:1:-1),nGHis_uzonal,latco,sqrt(rcl(1:ibLim)))
       IF(dogrh(nGHis_vmerid,latco)) CALL StoreGridHistory( gv(1:ibLim,kMax:1:-1),nGHis_vmerid,latco,sqrt(rcl(1:ibLim)))
       IF(dogrh(nGHis_spchum,latco)) CALL StoreGridHistory( gq(1:ibLim,kMax:1:-1),nGHis_spchum,latco)

       !  Surface is now at kMax
       IF(dogrh(nGHis_tvirsf,latco)) CALL StoreGridHistory(gtmp(1:ibLim,kMax),nGHis_tvirsf,latco)
       IF(dogrh(nGHis_uzonsf,latco)) CALL StoreGridHistory( gu(1:ibLim,kMax),nGHis_uzonsf,latco,sqrt(rcl(1:ibLim)))
       IF(dogrh(nGHis_vmersf,latco)) CALL StoreGridHistory( gv(1:ibLim,kMax),nGHis_vmersf,latco,sqrt(rcl(1:ibLim)))
       IF(dogrh(nGHis_sphusf,latco)) CALL StoreGridHistory( gq(1:ibLim,kMax),nGHis_sphusf,latco)

       IF(dogrh(nGHis_snowdp,latco)) CALL StoreGridHistory(sheleg(1:ibLim,latco),nGHis_snowdp,latco)
       IF(dogrh(nGHis_rouglg,latco)) CALL StoreGridHistory(  zorl(1:ibLim,latco),nGHis_rouglg,latco)

       ncount=0
       DO i=1,ibLim
          gdt(i,1)=gtsea(i,latco)
          gdt(i,2)=gtsea(i,latco)
          gdt(i,3)=gtsea(i,latco)
          gdt(i,4)=1.0_r8
          gdt(i,5)=1.0_r8
          gdt(i,6)=1.0_r8
          gdt(i,7)=0.0e0_r8
          gdt(i,8)=0.0e0_r8
          gdt(i,9)=imask(i,latco)
          IF(imask(i,latco).ge.1_i8)THEN
             ncount=ncount+1
             gdt(i,1)=tcm(ncount,latco)
             gdt(i,2)=tgm(ncount,latco)
             gdt(i,3)=tdm(ncount,latco)
             gdt(i,4)=wm (ncount,1,latco)
             gdt(i,5)=wm (ncount,2,latco)
             gdt(i,6)=wm (ncount,3,latco)
             gdt(i,7)=capacm(ncount,1,latco)
             gdt(i,8)=capacm(ncount,2,latco)
          END IF
       END DO

       IF(dogrh(nGHis_tcanop,latco)) CALL StoreGridHistory(gdt(1:ibLim,1),nGHis_tcanop,latco)
       IF(dogrh(nGHis_tgfccv,latco)) CALL StoreGridHistory(gdt(1:ibLim,2),nGHis_tgfccv,latco)
       IF(dogrh(nGHis_tgdeep,latco)) CALL StoreGridHistory(gdt(1:ibLim,3),nGHis_tgdeep,latco)
       IF(dogrh(nGHis_swtsfz,latco)) CALL StoreGridHistory(gdt(1:ibLim,4),nGHis_swtsfz,latco)
       IF(dogrh(nGHis_swtrtz,latco)) CALL StoreGridHistory(gdt(1:ibLim,5),nGHis_swtrtz,latco)
       IF(dogrh(nGHis_swtrcz,latco)) CALL StoreGridHistory(gdt(1:ibLim,6),nGHis_swtrcz,latco)
       IF(dogrh(nGHis_mostca,latco)) CALL StoreGridHistory(gdt(1:ibLim,7),nGHis_mostca,latco,1000.0_r8)
       IF(dogrh(nGHis_mostgc,latco)) CALL StoreGridHistory(gdt(1:ibLim,8),nGHis_mostgc,latco,1000.0_r8)
       IF(dogrh(nGHis_vegtyp,latco)) CALL StoreGridHistory(gdt(1:ibLim,9),nGHis_vegtyp,latco)

    END IF
    gdt=0.0_r8
    !
    !     computation of maximum wind 
    !     ---------------------------
    !
    DO k=1,kMax
       DO i=1,ibLim
          vmax(k)=MAX(vmax(k),cosiv(i)*SQRT(gu(i,k)*gu(i,k)+gv(i,k)*gv(i,k)))
       ENDDO
    ENDDO
    !
    !     Computation of tendencies (part related to intermediate time-step)
    !     ------------------------------------------------------------------
    !
    !     wind derivatives with respect to phi 
    !     ------------------------------------
    !
    IF (.NOT. slagr) THEN
       CALL delwind(gulam,gvlam,grot,gdiv,guphi,gvphi,cos2lat,&
            ibMax, kMax, ibLim)
    END IF
    !
    !     computation of gw, eta_dp_eta and vertical integrals of div and lnps
    !     --------------------------------------------------------------------
    !
    CALL vertint(slagr,zsint,psint,adveps,divint,gw,ps, &
        a_hybr,b_hybr,delb,delp,gprsi,rpi,c_hybr,eta_dp_eta,beta, &
        alpha,gu,gv,gdiv,gplam,gpphi,gzslam,gzsphi,rcl,ibMax,kMax,ibLim)
    !
    IF (.NOT. slagr) THEN
       !
       !     eulerian horizontal advection of wind
       !     -------------------------------------
       !
       CALL hadvec(gu,gv,gulam,guphi,gyu,rcl,ibMax,kMax,ibLim)
       CALL hadvec(gu,gv,gvlam,gvphi,gyv,rcl,ibMax,kMax,ibLim)
       !
       !     eulerian vertical advection of wind
       !     -----------------------------------
       !
       CALL vadvec(gu, eta_dp_eta, delp, gyu, ibMax, kMax, ibLim)
       CALL vadvec(gv, eta_dp_eta, delp, gyv, ibMax, kMax, ibLim)
       !
       !     metric term
       !     -----------
       !
       CALL metric(gu,gv,gyv,ercossin,ibMax,kMax,ibLim)
       !
    ENDIF
    !
    !     coriolis terms
    !     --------------
    !
    CALL coriol(gu,gv,gyu,gyv,fcor,ibMax,kMax,ibLim)
    !
    !     non-linear part of pressure gradient
    !     ------------------------------------
    !
    CALL nlprgr(gplam,gpphi,ps,gtlam,gtphi,gtmp,gyu,gyv,gasr, & 
                alpha,alpha_r,rpi,rpi_r,delp,delb,b_hybr,beta, &
                tref,gprsi,gphil,gphii,gzs,gzslam,gzsphi,ibMax,kMax,ibLim)
    !
    IF (.NOT. slagr) THEN
    !
    !     eulerian horizontal advection of temperature 
    !     --------------------------------------------
    !
       CALL hadvec(gu,gv,gtlam,gtphi,gtd,rcl,ibMax,kMax,ibLim)
    !
    !     eulerian vertical advection of temperature 
    !     ------------------------------------------
    !
       CALL vadvec(gtmp, eta_dp_eta, delp, gtd, ibMax, kMax, ibLim)
    !
    END IF
    !
    !     complete non-linear part of temperature and log pressure tendencies
    !     -------------------------------------------------------------------
    !
    CALL tmp_lnp_tend(gtd,gvdlnp,gtmp,gdiv,gq,ps,tref,psint,adveps, &
               divint,rk,cp,cpv,qpsm1,rpi,delp,alpha,c_hybr,rpi_r,delp_r, &
               delb,zsint,slagr,psref,ibMax,kMax,ibLim,latco,omg)
    !
    !     computation of omega and pressure at full levels
    !     ------------------------------------------------
    !
    CALL press_and_omega (gprsl, gprsi, omg, ibMax, kMax, ibLim)
    !
    !     Computation of Diagnostic for Transportation Fluxes
    !     ---------------------------------------------------
    !
    !           invert fields for Diagnostics (from top to bottom to bottom to top )
    !
    !     Horizontal Momentum Transport
    IF (dodia(nDiag_homtvu)) CALL updia (gv(1:ibLim,kMax:1:-1)*gu(1:ibLim,kMax:1:-1), nDiag_homtvu, latco)
    !     Vertical Zonal Momentum Transport
    IF (dodia(nDiag_vzmtwu)) CALL updia (omg(1:ibLim,kMax:1:-1)*gu(1:ibLim,kMax:1:-1),  nDiag_vzmtwu, latco)
    !     Vertical Meridional Momentum Transport
    IF (dodia(nDiag_vmmtwv)) CALL updia (omg(1:ibLim,kMax:1:-1)*gv(1:ibLim,kMax:1:-1),  nDiag_vmmtwv, latco)
    IF (dodia(nDiag_mshtvt) .OR. dodia(nDiag_zshtut) .OR. dodia(nDiag_vshtwt)) THEN
       !  Get Dry Absolute Temperature
       DO k=1,kMax
          DO i=1,ibLim
             gdt(i,k)=gtmp(i,k)/(1.0_r8+0.608_r8*gq(i,k))
          END DO
       END DO
       !     Meridional Sensible Heat Transport
       IF (dodia(nDiag_mshtvt)) CALL updia (gv(1:ibLim,kMax:1:-1)*gdt(1:ibLim,kMax:1:-1),  nDiag_mshtvt, latco)
       !     Zonal Sensible Heat Transport
       IF (dodia(nDiag_zshtut)) CALL updia (gu(1:ibLim,kMax:1:-1)*gdt(1:ibLim,kMax:1:-1),  nDiag_zshtut, latco)
       !     Vertical Sensible Heat Transport
       IF (dodia(nDiag_vshtwt)) CALL updia (omg(1:ibLim,kMax:1:-1)*gdt(1:ibLim,kMax:1:-1), nDiag_vshtwt, latco)
    END IF
    !     Meridional Specific Humidity Transport
    IF (dodia(nDiag_mshtuq)) CALL updia (gv(1:ibLim,kMax:1:-1)*gq(1:ibLim,kMax:1:-1),   nDiag_mshtuq, latco)
    !     Zonal Specific Humidity Transport
    IF (dodia(nDiag_zshtuq)) CALL updia (gu(1:ibLim,kMax:1:-1)*gq(1:ibLim,kMax:1:-1),   nDiag_zshtuq, latco)
    !     Vertical Specific Humidity Transport
    IF (dodia(nDiag_vshtwq)) CALL updia (omg(1:ibLim,kMax:1:-1)*gq(1:ibLim,kMax:1:-1),  nDiag_vshtwq, latco)
    !
    !     Dew Point Temperature K
    IF (dodia(nDiag_dewptt)) THEN
       DO k=1,kMax
          ! invert to bottom to top for diagnostics
          ki = kMax + 1 -k
          DO i=1,ibLim
             dewpwtr(i,ki)=(log(gq(i,k)*(a_hybr(k)+ps(i)*b_hybr(k))/380.042_r8)&
                           * 35.86_r8 - (4717.4732_r8))&
                           /(log(gq(i,k)*(a_hybr(k)+ps(i)*b_hybr(k))&
                           /380.042_r8) - 17.27_r8)
          END DO
       END DO
       CALL updia (dewpwtr(1:ibLim,1:kMax),nDiag_dewptt,latco)
    END IF
    !
    !
    IF (.NOT. slagr .AND. .NOT. slhum) THEN
       !
       !     eulerian horizontal advection of humidity 
       !     -----------------------------------------
       !
       CALL hadvec(gu,gv,gqlam,gqphi,gqd,rcl,ibMax,kMax,ibLim)
       !
       !     eulerian vertical advection of humidity
       !     ---------------------------------------
       !
       CALL vadvec(gq, eta_dp_eta, delp, gqd, ibMax, kMax, ibLim)
       !
    ENDIF
    IF(dodia(nDiag_divgxq)) CALL updia(psint (1:ibLim,kMax:1:-1),nDiag_divgxq,latco)
    IF(dodia(nDiag_vmoadv)) CALL updia(divint(1:ibLim,kMax:1:-1),nDiag_vmoadv,latco)
    IF(dodia(nDiag_omegav)) CALL updia(omg   (1:ibLim,kMax:1:-1),nDiag_omegav,latco)
    !
    !
    !
    IF (IsGridHistoryOn()) THEN
       IF(dogrh(nGHis_presfc,latco)) CALL StoreGridHistory(ps(1:ibLim),nGHis_presfc,latco,10.0_r8)
    END IF
    !     
    !     sigma gke computed only at interior interfaces.
    !
    !
    IF(dodia(nDiag_sigdot))CALL updia(gw(1:ibLim,kMax:1:-1),nDiag_sigdot,latco)

    !
    !     tendency from old time-step
    !     ---------------------------
    IF (slagr.and.SL_twotime_scheme) THEN
       IF (initial.ne.2) THEN
          CALL tndtold(gyum,gyvm,gtdm,gvdlnpm,gplam,gpphi,gzslam,gzsphi,gdiv,  &
                       gtlam,gtphi,alpha_r,rpi_r,delp_r,tref,psref,gasr,rk, &
                       ibMax,kMax,ibLim)
       ENDIF
    ELSE
       CALL tndtold(gyum,gyvm,gtdm,gvdlnpm,gplamm,gpphim,gzslam,gzsphi,gdivm,  &
                    gtlamm,gtphim,alpha_r,rpi_r,delp_r,tref,psref,gasr,rk, &
                    ibMax,kMax,ibLim)
    ENDIF

    IF(TRIM(isimp).ne.'YES')THEN
       !
       !     gplam surface pressure in mb
       !
       IF (microphys) THEN
          IF((nClass+nAeros)>0 .and. PRESENT(gvarm))THEN
             CALL DryPhysics & 
               (gzs, gtm, gqm, gum, gvm, ps, gyu, gyv, gtd, gqd, colrad, &
                ifday , tod, gtmp, gq, omg, gprsl, gprsi, gphil, gphii, &
                jb, lonrad, glnpm, cos2d, intcosz, &
                gicem  ,gice   ,gicet , &
                gliqm  ,gliq   ,gliqt ,&
                gvarm  ,gvar   ,gvart)
          ELSE
             CALL DryPhysics & 
               (gzs, gtm, gqm, gum, gvm, ps, gyu, gyv, gtd, gqd, colrad, &
                ifday , tod, gtmp, gq, omg, gprsl, gprsi, gphil, gphii, &
                jb, lonrad, glnpm, cos2d, intcosz, &
                gicem  ,gice   ,gicet , &
                gliqm  ,gliq   , gliqt  )
          END IF
         ELSE
          CALL DryPhysics & 
               (gzs, gtm, gqm, gum, gvm, ps, gyu, gyv, gtd, gqd, colrad, &
                ifday , tod, gtmp, gq, omg, gprsl, gprsi, gphil, gphii, &
                jb, lonrad, glnpm, cos2d, intcosz)
       ENDIF

       !     
       !     diagnostic of precipitable water
       !
       CALL pwater(gq    ,dewpwtr   ,delp  ,ibMax,ibLim ,kMax  )
       IF(dodia(nDiag_pwater))CALL updia(dewpwtr(1:ibLim,kMax:1:-1) ,nDiag_pwater,latco)
       !
    END IF
    !
    !     diagnostic of time mean surface pressure
    !
    IF(dodia(nDiag_tmpsfc))CALL updia(ps(1:ibLim)/1000.0_r8  ,nDiag_tmpsfc,latco)
    !
    IF (initial.eq.1.and.slagr.and.SL_twotime_scheme) THEN
       gyum = gyum - 0.5_r8 * gyu 
       gyvm = gyvm - 0.5_r8 * gyv
       gtdm = gtdm - 0.5_r8 * gtd
       gqdm = gqdm - 0.5_r8 * gqd
       gvdlnpm = gvdlnpm - 0.5_r8 * gvdlnp
    END IF
  END SUBROUTINE GrpComp






  SUBROUTINE delwind(ulam, vlam, vor, div, uphi, vphi, cos2lat, &
       ibMax, kMax, ibLim)
    INTEGER, INTENT(IN ) :: ibMax
    INTEGER, INTENT(IN ) :: kMax
    REAL(KIND=r8),    INTENT(IN ) :: cos2lat(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: ulam(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN ) :: vlam(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN ) :: div (ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN ) :: vor (ibMax,kMax)
    REAL(KIND=r8),    INTENT(OUT) :: uphi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(OUT) :: vphi(ibMax,kMax)
    INTEGER, INTENT(IN ) :: ibLim
    INTEGER :: ib, k
    !      
    !      From the vorticity, divergence and the e-w derivatives of 
    !      U and V, computes the values of cos(phi) d/d phi F , where F = U,
    !      and F = V.
    !
    DO k=1,kMax   
       DO ib=1,ibLim
          uphi(ib,k) = vlam(ib,k)  - cos2lat(ib) * vor(ib,k)
          vphi(ib,k) = cos2lat(ib) * div(ib,k) - ulam(ib,k)
       ENDDO
    ENDDO
  END SUBROUTINE delwind



  SUBROUTINE press_and_omega (press, phalf, omega, ibMax, kMax, ibLim)

    INTEGER, INTENT(IN ) :: ibMax
    INTEGER, INTENT(IN ) :: kMax
    REAL(KIND=r8),    INTENT(IN )   :: phalf(ibMax,kMax+1)
    REAL(KIND=r8),    INTENT(OUT)   :: press(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: omega(ibMax,kMax)
    INTEGER, INTENT(IN ) :: ibLim
    INTEGER :: ib, k
    !      
    !      From pressure at interfaces compute the layer pressures and
    !      finish evaluation of w (w/p at entrance)
    !
    !       Pressures are in Pascal!  Omega is in Pa / s !
    !
    DO k=1,kMax   
       DO ib=1,ibLim
          press(ib,k) = 0.5_r8 * (phalf(ib,k) + phalf(ib,k+1))
          omega(ib,k) = omega(ib,k) * press(ib,k)
       ENDDO
    ENDDO
  END SUBROUTINE press_and_omega


  SUBROUTINE vertint(slagr,zsint,psint, adveps, divint, w, ps, &
     a_hybr, b_hybr, delb, delp, phalf, rpi, c_hybr, eta_dp_eta, beta, &
     alpha, u, v, div, plam, pphi, zlam, zphi, rcl, ibMax, kMax, ibLim)
    LOGICAL, INTENT(IN ) :: slagr
    INTEGER, INTENT(IN ) :: ibMax
    INTEGER, INTENT(IN ) :: kMax
    INTEGER, INTENT(IN ) :: ibLim
    REAL(KIND=r8),    INTENT(IN ) :: u(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN ) :: v(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN ) :: div(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN ) :: ps(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: plam(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: pphi(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: zlam(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: zphi(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: rcl(ibMax)
    REAL(KIND=r8),    INTENT(IN ) :: a_hybr(kMax+1)
    REAL(KIND=r8),    INTENT(IN ) :: b_hybr(kMax+1)
    REAL(KIND=r8),    INTENT(IN ) :: c_hybr(kMax)
    REAL(KIND=r8),    INTENT(IN ) :: delb(kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: delp(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: phalf(ibMax,kMax+1)
    REAL(KIND=r8),    INTENT(INOUT) :: beta(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: alpha(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: eta_dp_eta(ibMax,kMax+1)
    REAL(KIND=r8),    INTENT(INOUT) :: rpi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: zsint(ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: psint(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: w(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: divint(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: adveps(ibMax,kMax)
    !      
    !      Pressure and related variables are derived from surface pressure.
    !      Computation of the vertical integral (in finite differences)
    !      of the divergence field (to be stored in divint). The scalar product  
    !      of the wind field (in each level) with the gradient of the surface
    !      pressure is stored in adveps. Its vertical integral is stored in psint.
    !      The vertical velocity  is stored in w.
    !
    INTEGER :: i, k
    REAL(KIND=r8) :: dela
    !

    
    DO i=1,ibLim
       phalf(i,kmax+1) = ps(i)
    ENDDO
    DO k=kMax,1,-1
       DO i=1,ibLim
          phalf(i,k) = a_hybr(k) + b_hybr(k) * ps(i)
          delp(i,k) = phalf(i,k+1)-phalf(i,k)
          adveps(i,k) = rcl(i)*(u(i,k) * plam(i) + v(i,k) * pphi(i) )
       ENDDO
    ENDDO
    k=1
    DO i=1,ibLim
       eta_dp_eta(i,kmax+1) = 0.0_r8
       eta_dp_eta(i,1) = 0.0_r8
       alpha(i,k) = log(2.0_r8)
       psint(i,k) = delb(k) * adveps(i,k)
       divint(i,k) = delp(i,k) * div(i,k)
    ENDDO
    IF (slagr) THEN
       DO i=1,ibLim
          zsint(i) = delb(k) * rcl(i)*(u(i,k)*zlam(i)+v(i,k)*zphi(i))
       ENDDO
    ENDIF
    DO k=2,kMax
       DO i=1,ibLim
          psint(i,k) = psint(i,k-1) + delb(k) * adveps(i,k)
          divint(i,k) = divint(i,k-1) + delp(i,k) * div(i,k)
          beta(i,k) =  - c_hybr(k) / (phalf(i,k+1) * phalf(i,k) )
          rpi(i,k) = log (phalf(i,k+1)/phalf(i,k))
          alpha(i,k) = 1.0_r8 - phalf(i,k) * rpi(i,k) / delp(i,k) 
       ENDDO
       IF (slagr) THEN
          DO i=1,ibLim
             zsint(i) = zsint(i) + delb(k)*rcl(i)*(u(i,k)*zlam(i)+v(i,k)*zphi(i))
          ENDDO
       ENDIF
    ENDDO
    DO k=2,kMax
       DO i=1,ibLim
          eta_dp_eta(i,k) = b_hybr(k) * ( divint(i,kmax) + ps(i) * psint(i,kmax) )&
                            -  divint(i,k-1) - ps(i) * psint(i,k-1) 
       ENDDO
    ENDDO
!PK    IF (slagr) THEN
       DO k=1,kMax
          dela = a_hybr(k+1) - a_hybr(k)
          DO i=1,ibLim
             w(i,k) = ( eta_dp_eta(i,k) + eta_dp_eta(i,k+1) ) / 2.0
             w(i,k) = w(i,k) * ( dela / p000 + delb(k) ) / &
                        ( dela + ps(i) * delb(k) )
          ENDDO
       ENDDO
!PK    ENDIF
  END SUBROUTINE vertint

  SUBROUTINE hadvec(u, v, flam, fphi, tend, rcl, ibMax, kMax, ibLim)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    REAL(KIND=r8),    INTENT(IN   ) :: rcl(ibMax) ! 1.0_r8 / ( cos(lat)**2 )
    REAL(KIND=r8),    INTENT(IN   ) :: u(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: v(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: flam(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: fphi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tend(ibMax,kMax)
    !      
    ! Computes the horizontal advection term of field f (whose horizontal
    ! derivatives are given in flam and fphi) and add its contribution to current
    ! tendency (stored in tend).
    !
    INTEGER :: i, k
    DO k=1,kMax   
       DO i=1,ibLim
          tend(i,k) = tend(i,k) - rcl(i) * (u(i,k)*flam(i,k) + v(i,k)*fphi(i,k))
       END DO
    END DO
  END SUBROUTINE hadvec





  SUBROUTINE vadvec(f, eta_dp_eta, delp, tend, ibMax, kMax, ibLim)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    REAL(KIND=r8),    INTENT(IN   ) :: f(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: eta_dp_eta(ibMax,kMax+1)
    REAL(KIND=r8),    INTENT(IN   ) :: delp(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tend(ibMax,kMax)
    !      
    ! Computes the vertical advection of field f (in finite differences)
    ! and add its contribution to current tendency (in tend)
    !
    INTEGER :: i, k
    REAL(KIND=r8) :: dp2
    DO k=1,kMax-1
       DO i=1,ibLim
          dp2 = 0.5_r8 * eta_dp_eta(i,k+1)*(f(i,k+1)-f(i,k))
          tend(i,k  ) = tend(i,k  ) - dp2 / delp(i,k  )
          tend(i,k+1) = tend(i,k+1) - dp2 / delp(i,k+1)
       ENDDO
    ENDDO
  END SUBROUTINE vadvec






  SUBROUTINE metric(u, v, tend, ercossin, ibMax, kMax, ibLim)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    REAL(KIND=r8),    INTENT(IN   ) :: ercossin(ibMax) ! sin(lat) / ( er * cos(lat)**2 )
    REAL(KIND=r8),    INTENT(IN   ) :: u(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: v(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tend(ibMax,kMax)
    INTEGER :: i, k
    !      
    ! Computes the metric term and add its contribution  to current v-tendency
    !     ercossin(1:ibLim) = COS(colrad(1:ibLim)) * rcl(1:ibLim) / er
    !                                             1.0_r8/cos(latitude)**2.0_r8
    !     ercossin(1:ibLim) =   1.0_r8/cos(latitude)/ er
    !                                             
    !
    DO k=1,kMax   
       DO i=1,ibLim
          tend(i,k) = tend(i,k) - ercossin(i)*(u(i,k)*u(i,k) + v(i,k)*v(i,k))
       END DO
    END DO
  END SUBROUTINE metric






  SUBROUTINE coriol(u, v, tendu, tendv, fcor, ibMax, kMax, ibLim)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    REAL(KIND=r8),    INTENT(IN   ) :: fcor(ibMax) ! 2 * omega * sin(phi)
    REAL(KIND=r8),    INTENT(IN   ) :: u(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: v(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tendu(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tendv(ibMax,kMax)
    INTEGER :: i, k
    !      
    ! Computes the coriolis contributions  to current u- and v- tendencies
    !
    !
    DO k=1,kMax   
       DO i=1,ibLim
          tendu(i,k) = tendu(i,k) + fcor(i) * v(i,k)
          tendv(i,k) = tendv(i,k) - fcor(i) * u(i,k)
       ENDDO
    ENDDO
  END SUBROUTINE coriol






  SUBROUTINE nlprgr(plam, pphi, ps, tlam, tphi, tmp, tendu, tendv, gasr, & 
                    alpha, alpha_r, rpi, rpi_r, delp, delb, b_hybr, beta, &
                    tref, phalf, phil, phii, zs, zslam, zsphi, ibMax, kMax, ibLim)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    REAL(KIND=r8),    INTENT(IN   ) :: gasr
    REAL(KIND=r8),    INTENT(IN   ) :: tref
    REAL(KIND=r8),    INTENT(IN   ) :: alpha(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: beta(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: alpha_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: rpi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: rpi_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: delp(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: delb(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: b_hybr(kMax+1)
    REAL(KIND=r8),    INTENT(IN   ) :: phalf(ibMax,kMax+1)
    REAL(KIND=r8),    INTENT(IN   ) :: tmp(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: ps(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: plam(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: pphi(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: zs(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: zslam(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: zsphi(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: tlam(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: tphi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: phil(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: phii(ibMax,kMax+1)
    REAL(KIND=r8),    INTENT(INOUT) :: tendu(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tendv(ibMax,kMax)
    REAL(KIND=r8) :: apu(ibMax), apv(ibMax)
    REAL(KIND=r8) :: au, av, at, cg, c
    INTEGER :: i, j, k
    !      
    ! Computes the non-linear part of the pressure gradient contribution for the 
    ! current u-  and v- tendencies.
    !
    ! Also evaluate the heights of the atmospheric layers and interfaces
    !
    apu = 0.0_r8
    apv = 0.0_r8
!   phii(1:iblim,kMax+1) = zs(1:iblim) / grav
    phii(1:iblim,kMax+1) = 0.0_r8
    
    cg = gasr / grav
    c = eps_sic
    DO k=kMax,1,-1   
       DO i=1,ibLim
          IF (k.lt.kmax) THEN
             j = k+1
             apu(i) = apu(i) + (rpi(i,j)-c*rpi_r(j))*tlam(i,j) + &
                              tmp(i,j)*beta(i,j)*ps(i)*plam(i)
             apv(i) = apv(i) + (rpi(i,j)-c*rpi_r(j))*tphi(i,j) + &
                              tmp(i,j)*beta(i,j)*ps(i)*pphi(i)
          ENDIF
!         at = ( (rpi(i,k)*b_hybr(k)+alpha(i,k)*delb(k))/delp(i,k) - &
!                  c_grad_alpha(i,k) ) * tmp(i,k)*ps(i)
          at = b_hybr(k+1) * tmp(i,k) * ps(i) / phalf(i,k+1)
          au = (alpha(i,k)-c*alpha_r(k))*tlam(i,k) + (at - c*tref) * plam(i)
          av = (alpha(i,k)-c*alpha_r(k))*tphi(i,k) + (at - c*tref) * pphi(i)
!         tendu(i,k) = tendu(i,k) - gasr * (au +apu(i)) - zslam(i)
!         tendv(i,k) = tendv(i,k) - gasr * (av +apv(i)) - zsphi(i)
          tendu(i,k) = tendu(i,k) - gasr * (au +apu(i))
          tendv(i,k) = tendv(i,k) - gasr * (av +apv(i))
          phii(i,k) = phii(i,k+1) + cg * rpi(i,k) * tmp(i,k)
          phil(i,k) = phii(i,k+1) + cg * alpha(i,k) * tmp(i,k)
       ENDDO
    ENDDO
  END SUBROUTINE nlprgr


  SUBROUTINE tmp_lnp_tend(tend_t, tend_l, tmp, div, q, ps, tref, psint, adveps, &
       divint, rk, cp, cpv, qpsm1, rpi, delp, alpha, c_hybr, rpi_r, delp_r, &
       delb, zsint, slagr, psref, ibMax, kMax, ibLim,latco, omg)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim,latco
    LOGICAL, INTENT(IN   ) :: slagr
    REAL(KIND=r8),    INTENT(IN   ) :: tref, psref, rk, cp, cpv, qpsm1 
    REAL(KIND=r8),    INTENT(IN   ) :: psint(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: divint(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: adveps(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: tmp(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: div(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: q(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: rpi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: delp(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: alpha(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: zsint(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: ps(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: c_hybr(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: rpi_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: delp_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: delb(kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: omg(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tend_t(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tend_l(ibMax)
    INTEGER :: i, k
    REAL(KIND=r8) :: dm1, s(ibMax), s1, s2, s3, cq, a
    !      
    ! Computes the non-linear contributions to the temperature and
    ! log pressure tendencies
    !
    dm1 = cpv / cp - 1.0_r8
    s = 0.0_r8
    s2 = rk*tref*eps_sic
    DO k=2,kmax
       s1 = s2*rpi_r(k)/delp_r(k)
       s3 = s2 * alpha_r(k)
       DO i=1,ibLim
          cq = (1.0_r8 + qpsm1 * q(i,k)) / (1.0_r8 + dm1 * q(i,k))
          a = rk * tmp(i,k) * cq 
          omg(i,k) =  rpi(i,k)*(psint(i,k-1)*ps(i) + divint(i,k-1)) + &
                      alpha(i,k) * (delp(i,k)*div(i,k)+&
                      adveps(i,k)*ps(i)*delb(k) ) - ps(i) * (delb(k) + &
                      c_hybr(k) * rpi(i,k) / delp(i,k) ) * adveps(i,k) 
 
          omg(i,k) = - omg(i,k) / delp(i,k)
          s(i) = s(i) + delp_r(k-1) * div(i,k-1)
          tend_t(i,k) = tend_t(i,k) + a * omg(i,k) + s1 * s(i) + s3 * div(i,k)
       ENDDO
    ENDDO
    k = 1   
    s3 = s2 * alpha_r(k)
    DO i=1,ibLim
       cq = (1.0_r8 + qpsm1 * q(i,k)) / (1.0_r8 + dm1 * q(i,k))
       a = rk * tmp(i,k) * cq 
       omg(i,k) = alpha(i,k) * (delp(i,k)*div(i,k) + &
                  adveps(i,k)*ps(i)*delb(k) ) - ps(i) * delb(k) * adveps(i,k)
       omg(i,k) = - omg(i,k) / delp(i,k)
       tend_t(i,k) = tend_t(i,k) +  a * omg(i,k) + s3 * div(i,k)
       s(i) = s(i) + delp_r(kmax) * div(i,kmax)
    ENDDO
    IF (slagr) THEN
       DO i=1,ibLim
          tend_l(i) = tend_l(i) - (divint(i,kmax)/ps(i)-eps_sic*s(i)/psref) + &
                                   zsint(i) / (gasr*tref)
       ENDDO
     ELSE
       DO i=1,ibLim
          tend_l(i) = tend_l(i) - (divint(i,kmax)/ps(i)-eps_sic*s(i)/psref) - psint(i,kmax)
       ENDDO
    ENDIF
    
    
  END SUBROUTINE tmp_lnp_tend






  SUBROUTINE tndtold(tdu, tdv, tdt, tdlnp, plam, pphi, zslam, zsphi, div, tlam, &
       tphi, alpha_r, rpi_r, delp_r, tref, psref, rc, rk, ibMax, kMax, ibLim)
    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: ibLim
    REAL(KIND=r8),    INTENT(IN   ) :: rc
    REAL(KIND=r8),    INTENT(IN   ) :: rk
    REAL(KIND=r8),    INTENT(IN   ) :: tref
    REAL(KIND=r8),    INTENT(IN   ) :: psref
    REAL(KIND=r8),    INTENT(IN   ) :: alpha_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: rpi_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: delp_r(kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: tlam(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: tphi(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: div(ibMax,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: plam(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: pphi(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: zslam(ibMax)
    REAL(KIND=r8),    INTENT(IN   ) :: zsphi(ibMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tdu(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tdv(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tdt(ibMax,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: tdlnp(ibMax)
    INTEGER :: i, k
    REAL(KIND=r8) :: half=0.5e0_r8
    REAL(KIND=r8) :: s, sk
    REAL(KIND=r8) :: intu(ibMax), intv(ibMax), intl(ibMax)

    !      
    ! Computes the part of the tendencies relative to the old time-step
    !
    !
    !   pressure gradient terms
    !   -----------------------
    intu = 0.0_r8
    intv = 0.0_r8
    intl = 0.0_r8
    DO k=kMax,1,-1  
       DO i=1,ibLim
          tdu(i,k) = tdu(i,k) - half * ( zslam(i) + rc * eps_sic * & 
                     ( intu(i) + alpha_r(k) * tlam(i,k) + tref * plam(i) ) )
          tdv(i,k) = tdv(i,k) - half * ( zsphi(i) + rc * eps_sic * & 
                     ( intv(i) + alpha_r(k) * tphi(i,k) + tref * pphi(i) ) )
          intu(i) = intu(i) + rpi_r(k) * tlam(i,k)
          intv(i) = intv(i) + rpi_r(k) * tphi(i,k)
       ENDDO
    ENDDO
    !
    !   Temperature tendency
    !   --------------------
    intl = 0.0_r8
    s = half * eps_sic * rk * tref
    DO k=1,kMax   
       sk = rpi_r(k) / delp_r(k)
       DO i=1,ibLim
          tdt(i,k) = tdt(i,k) - s * ( sk * intl(i) + alpha_r(k) * div(i,k) )
          intl(i) = intl(i) + delp_r(k) * div(i,k)
       ENDDO
    ENDDO
    !
    !   log pressure tendency
    !   ---------------------
    s = half * eps_sic / psref
    DO i=1,ibLim
       tdlnp(i) = tdlnp(i) - s * intl(i)
    ENDDO
  END SUBROUTINE tndtold

  !
  ! addtend: finish tendency computations, adding contributions from
  !          old and current time step.

  SUBROUTINE AddTend(dt, nlnminit, jbFirst, jbLast, slhum)
    REAL(KIND=r8),    INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: nlnminit
    LOGICAL, INTENT(IN) :: slhum
    INTEGER, INTENT(IN) :: jbFirst
    INTEGER, INTENT(IN) :: jbLast
    INTEGER :: ib, jb, k
    REAL(KIND=r8) :: dt2
    dt2 = dt + dt
    IF (.NOT.nlnminit) THEN

       DO jb = jbFirst, jbLast
          DO k = 1, kMax
             DO ib = 1, ibMaxPerJB(jb)
                fgyu(ib,k,jb) = fgum(ib,k,jb) + &
                     dt2 * ( fgyu(ib,k,jb) + fgyum(ib,k,jb) ) 
                fgyv(ib,k,jb) = fgvm(ib,k,jb) + &
                     dt2 * ( fgyv(ib,k,jb) + fgyvm(ib,k,jb) ) 
                fgtd(ib,k,jb) = fgtmpm(ib,k,jb) + &
                     dt2 * ( fgtd(ib,k,jb) + fgtdm(ib,k,jb) ) 
             END DO
             IF (.not.slhum) THEN
                DO ib = 1, ibMaxPerJB(jb)
                   fgqd(ib,k,jb) = fgqm(ib,k,jb) + &
                        dt2 * fgqd(ib,k,jb) 
                END DO
             ENDIF
          END DO
          DO ib = 1, ibMaxPerJB(jb)
             fgvdlnp(ib,jb) = fglnpm(ib,jb) + &
                  dt2 * ( fgvdlnp(ib,jb) + fgvdlnpm(ib,jb) )
          END DO
       END DO

    ELSE

       DO jb = jbFirst, jbLast
          DO k = 1, kMax
             DO ib = 1, ibMaxPerJB(jb)
                fgyu(ib,k,jb) = fgyu(ib,k,jb) + 2.0_r8 * fgyum(ib,k,jb)
                fgyv(ib,k,jb) = fgyv(ib,k,jb) + 2.0_r8 * fgyvm(ib,k,jb)
                fgtd(ib,k,jb) = fgtd(ib,k,jb) + 2.0_r8 * fgtdm(ib,k,jb)
                fgqd(ib,k,jb) = fgqd(ib,k,jb)
             END DO
          END DO
          DO ib = 1, ibMaxPerJB(jb)
             fgvdlnp(ib,jb) = fgvdlnp(ib,jb) + 2.0_r8 * fgvdlnpm(ib,jb)
          END DO
       END DO

    ENDIF
  END SUBROUTINE AddTend

  ! SetJablo:  sets initial fields for Jablonowski test case

  SUBROUTINE SetJablo
    INTEGER :: ib, jb, k, i, j
    REAL(KIND=r8) :: tsig, phisig, sigv, cossigv, a, ag, b, c
    REAL(KIND=r8) :: zlon, zlat, szlat, czlat, up, br, pai, paihalf
    REAL(KIND=r8) :: slat, clat, s2lat, alon, uprime
    REAL(KIND=r8),   PARAMETER :: p0 = 1.0d05
    REAL(KIND=r8),   PARAMETER :: h  = 7.340d0
    REAL(KIND=r8),   PARAMETER :: u0 = 35.d0
    REAL(KIND=r8),   PARAMETER :: T0 = 288.d0
    REAL(KIND=r8),   PARAMETER :: st = 0.2d0
    REAL(KIND=r8),   PARAMETER :: delT0 = 4.8d05
    REAL(KIND=r8),   PARAMETER :: omg = 7.29212d-05
    REAL(KIND=r8),   PARAMETER :: rearth  = 6.371229d06
    REAL(KIND=r8),   PARAMETER :: r = 287.d0
    REAL(KIND=r8),   PARAMETER :: g = 9.80616d0
    REAL(KIND=r8),   PARAMETER :: eta0 = 0.252d0
    REAL(KIND=r8),   PARAMETER :: gama = 0.005d0
    REAL(KIND=r8),   PARAMETER :: rgamg = r*gama/g
    REAL(KIND=r8),   PARAMETER :: tggam = T0*g/gama
    REAL(KIND=r8),   PARAMETER :: rdt = r * delT0

    pai = acos(-1.d0)
    paihalf = pai / 2.d0
    DO k=1,kmax
       tsig = T0 * eta(k) ** rgamg - 300.d0
!      phisig = tggam * (1.d0-eta(k) ** rgamg)
       if (eta(k).lt.0.20d0) tsig = tsig + delT0 * (0.20d0-eta(k)) ** 5
!      if (eta(k).lt.0.20d0) &
!           phisig = phisig - rdt * ( (log(eta(k)/0.2d0) + &
!           137.d0 / 60.d0) * 0.2d0 ** 5 - 5.d0 * eta(k) * 0.2d0 ** 4 &
!          + 5.d0 * eta(k) ** 2 * 0.2d0 ** 3 - 10.d0 / 3.d0 * &
!          eta(k) ** 3 * 0.2d0 ** 2 + 1.25d0 * eta(k) ** 4 * 0.2d0 &
!          - eta(k) ** 5 * 0.2d0 )
       sigv = (eta(k)-eta0)*paihalf
       cossigv = cos(sigv)**1.50d0
       a = 0.750d0 * eta(k)*pai*u0/r*sin(sigv)*sqrt(cos(sigv))
!      ag = u0*cossigv
       DO j = 1, jbmax
          DO i = 1, ibMaxPerJB(j)
             slat = sinlat(i,j)
             clat = coslat(i,j)
             s2lat = 2.0d0 * slat * clat
             b = 2.0d0*u0*cossigv*(10.0d0/63.0d0- &
                 2.0d0*slat**6*(clat**2+1.d0/3.0d0))
             c = omg*rearth*(1.6d0*clat**3*(slat**2+2.d0/3.0d0)-pai/4.d0)
             fgu(i,k,j) = u0 * cossigv * s2lat * s2lat
             fgtmp(i,k,j) =  tsig + a * (b+c)
!            geop(i,k,j) = phisig + ag * (b/2.d0+c)
             fgv(i,k,j) = 0.0
          ENDDO
       ENDDO
    ENDDO
    sigv = (1.0d0-eta0)*paihalf
    cossigv = u0*cos(sigv)**1.5d0
    zlon = pai / 9.d0
    zlat = 2.d0 * pai / 9.d0
    szlat = sin(zlat)
    czlat = cos(zlat)
    up = 1.
    DO j = 1, jbmax
       DO i = 1, ibMaxPerJB(j)
          slat = sinlat(i,j)
          clat = coslat(i,j)
          b = (10.d0/63.d0-2.d0*slat**6*(clat**2+1.d0/3.d0))*cossigv
          c = omg*rearth*(1.6d0*clat**3*(slat**2+2.d0/3.d0)-pai/4.d0)
          alon = longit(i,j)
          br = szlat*slat+czlat*clat*cos(alon-zlon)
          br = 10.d0 * acos(br)
          uprime = up * exp(-br**2)
          DO k=1,kmax
             fgu(i,k,j) = (fgu(i,k,j) + uprime) * coslat(i,j)
!            fgu(i,k,j) = fgu(i,k,j) * coslat(i,j)
          ENDDO
          fgzs(i,j) = cossigv * (b+c)
          fglnps(i,j) = log(p0)
          fgps(i,j) = p0
       ENDDO
    ENDDO

  END SUBROUTINE SetJablo


  ! TimeFilterStep1: First part of the asselin/robert time-filter
  ! (computes a partially filtered value of fold)

  SUBROUTINE TimeFilterStep1(fa, fb,fb1, jbFirst, jbLast, slhum)
    REAL(KIND=r8),    INTENT(IN) :: fa
    REAL(KIND=r8),    INTENT(IN) :: fb
    REAL(KIND=r8),    INTENT(IN) :: fb1
    INTEGER, INTENT(IN) :: jbFirst
    INTEGER, INTENT(IN) :: jbLast
    LOGICAL, INTENT(IN) :: slhum
    INTEGER :: ib, jb, k,kk
    REAL(KIND=r8):: alfa1, sigma,fc
    !alfa=0.53_r8
    alfa1=alfa                          !filtro de Willians
    sigma=(1.0_r8-fa)/alfa1     !filtro de Asselin
    fc=(1.0_r8-alfa1)*sigma/2.0_r8

    DO jb = jbFirst, jbLast
       DO k = 1, kMax
          DO ib = 1, ibMaxPerJB(jb)
             fgtmpm  (ib,k,jb) = fa*fgtmp (ib,k,jb) + fb*fgtmpm (ib,k,jb)
             fgdivm  (ib,k,jb) = fa*fgdiv (ib,k,jb) + fb*fgdivm (ib,k,jb)
             fgum    (ib,k,jb) = fa*fgu   (ib,k,jb) + fb*fgum   (ib,k,jb)
             fgvm    (ib,k,jb) = fa*fgv   (ib,k,jb) + fb*fgvm   (ib,k,jb)
             fgtlamm (ib,k,jb) = fa*fgtlam(ib,k,jb) + fb*fgtlamm(ib,k,jb)
             fgtphim (ib,k,jb) = fa*fgtphi(ib,k,jb) + fb*fgtphim(ib,k,jb)
          END DO
       END DO
       IF(.not.slhum) THEN
          DO k = 1, kMax
             DO ib = 1, ibMaxPerJB(jb)
                fgqm(ib,k,jb) = fa*fgq(ib,k,jb) + fb*fgqm(ib,k,jb)
             END DO
          END DO
        ELSE 
          IF (fa.ne.0.0_r8) THEN 
             DO k = 1, kMax
                DO ib = 1, ibMaxPerJB(jb)
                   fgqm(ib,k,jb) = fa*fgq(ib,k,jb) + fb*(fgqm(ib,k,jb)+fgqp(ib,k,jb))
                END DO
             END DO
             IF (microphys) THEN
                DO k = 1, kMax
                   DO ib = 1, ibMaxPerJB(jb)
                      fgicem(ib,k,jb) = fa*fgice(ib,k,jb) + fb*(fgicem(ib,k,jb)+fgicep(ib,k,jb))
                      fgliqm(ib,k,jb) = fa*fgliq(ib,k,jb) + fb*(fgliqm(ib,k,jb)+fgliqp(ib,k,jb))
                   END DO
                END DO
                DO kk=1,nClass+nAeros
                   DO k = 1, kMax
                      DO ib = 1, ibMaxPerJB(jb)
                         fgvarm(ib,k,jb,kk) = fa*fgvar(ib,k,jb,kk) + fb*(fgvarm(ib,k,jb,kk)+fgvarp(ib,k,jb,kk))
                      END DO
                   END DO
                END DO
             ENDIF
          ENDIF
          DO k = 1, kMax
             DO ib = 1, ibMaxPerJB(jb)
                fgq(ib,k,jb) = fgqp(ib,k,jb)
             END DO
          END DO
          IF (microphys) THEN
             DO k = 1, kMax
                DO ib = 1, ibMaxPerJB(jb)
                   fgice(ib,k,jb) = fgicep(ib,k,jb)
                   fgliq(ib,k,jb) = fgliqp(ib,k,jb)
                END DO
             END DO
             DO kk=1,nClass+nAeros
                DO k = 1, kMax
                   DO ib = 1, ibMaxPerJB(jb)
                      fgvar(ib,k,jb,kk) = fgvarp(ib,k,jb,kk)
                   END DO
                END DO
             END DO
          ENDIF
       ENDIF
       DO ib = 1, ibMaxPerJB(jb)
          fglnpm  (ib,jb) = fa*fglnps(ib,jb) + fb*fglnpm (ib,jb)
          fgplamm (ib,jb) = fa*fgplam(ib,jb) + fb*fgplamm(ib,jb)
          fgpphim (ib,jb) = fa*fgpphi(ib,jb) + fb*fgpphim(ib,jb)
       END DO
    END DO
  END SUBROUTINE TimeFilterStep1


  ! TimeFilterStep2: Second part of the asselin/robert time-filter
  ! (the partially filtered value of fold is filtered completely)

  SUBROUTINE TimeFilterStep2(fb1, fa,jbFirst, jbLast, slhum, bckhum)
    REAL(KIND=r8),    INTENT(IN) :: fb1
    REAL(KIND=r8),    INTENT(IN) :: fa
    INTEGER, INTENT(IN) :: jbFirst
    INTEGER, INTENT(IN) :: jbLast
    LOGICAL, INTENT(IN) :: slhum
    LOGICAL, INTENT(IN) :: bckhum
    INTEGER :: ib, jb, k
    REAL(KIND=r8):: alfa1, sigma,fc
    !alfa=0.53_r8
    alfa1=alfa
    sigma=(1.0_r8-fa)/alfa1
    fc=(1.0_r8-alfa1)*sigma/2.0_r8

    DO jb = jbFirst, jbLast
       DO k = 1, kMax
          DO ib = 1, ibMaxPerJB(jb)
             fgtmpm (ib,k,jb) = fgtmpm (ib,k,jb) + fb1*fgtmp (ib,k,jb)
             fgdivm (ib,k,jb) = fgdivm (ib,k,jb) + fb1*fgdiv (ib,k,jb)
             fgum   (ib,k,jb) = fgum   (ib,k,jb) + fb1*fgu   (ib,k,jb)
             fgvm   (ib,k,jb) = fgvm   (ib,k,jb) + fb1*fgv   (ib,k,jb)
             fgtlamm(ib,k,jb) = fgtlamm(ib,k,jb) + fb1*fgtlam(ib,k,jb)
             fgtphim(ib,k,jb) = fgtphim(ib,k,jb) + fb1*fgtphi(ib,k,jb)
          END DO
       END DO
       IF(.not.slhum.or.bckhum) THEN
          DO k = 1, kMax
             DO ib = 1, ibMaxPerJB(jb)
                fgqm   (ib,k,jb) = fgqm   (ib,k,jb) + fb1*fgq   (ib,k,jb)
             END DO
          END DO
!          IF (microphys) THEN
!             DO k = 1, kMax
!                DO ib = 1, ibMaxPerJB(jb)
!                   fgicem(ib,k,jb) = fgicem(ib,k,jb) + fb1*(fgice(ib,k,jb))
!                   fgliqm(ib,k,jb) = fgliqm(ib,k,jb) + fb1*(fgliq(ib,k,jb))
!!!                END DO
!             END DO
!          ENDIF
       ENDIF
       DO ib = 1, ibMaxPerJB(jb)
          fglnpm (ib,jb) = fglnpm (ib,jb) + fb1*fglnps(ib,jb)
          fgplamm(ib,jb) = fgplamm(ib,jb) + fb1*fgplam(ib,jb)
          fgpphim(ib,jb) = fgpphim(ib,jb) + fb1*fgpphi(ib,jb)
       END DO
    END DO
  END SUBROUTINE TimeFilterStep2


  SUBROUTINE End_temp_difus(jbFirst, jbLast)
    INTEGER, INTENT(IN) :: jbFirst
    INTEGER, INTENT(IN) :: jbLast
    INTEGER :: ib, jb, k
    REAL(KIND=r8) :: kh
    REAL(KIND=r8) :: delp(ibMax,kMax)
    REAL(KIND=r8) :: ph(ibMax,kMax+1)

    DO jb = jbFirst, jbLast
       DO ib = 1, ibMaxPerJB(jb)
          ph(ib,kMax+1) = fgps(ib,jb)
       ENDDO
       DO k=kMax,1,-1
          DO ib = 1, ibMaxPerJB(jb)
             ph(ib,k) = a_hybr(k) + b_hybr(k) * fgps(ib,jb)
             delp(ib,k) = ph(ib,k+1)-ph(ib,k)
          ENDDO
       ENDDO
       k=1
       DO ib = 1, ibMaxPerJB(jb)
          ph(ib,k) =  fgps(ib,jb) * &
              ( b_hybr(k+1) * (fgtmp(ib,k+1,jb)-fgtmp(ib,k,jb)) ) / delp(ib,k)
       ENDDO
       DO k=2,kMax-1
          DO ib = 1, ibMaxPerJB(jb)
             ph(ib,k) = fgps(ib,jb) * &
                  ( b_hybr(k+1) * (fgtmp(ib,k+1,jb)-fgtmp(ib,k  ,jb)) + &
                    b_hybr(k  ) * (fgtmp(ib,k  ,jb)-fgtmp(ib,k-1,jb)) ) / delp(ib,k)
          ENDDO
       ENDDO
       k=kMax
       DO ib = 1, ibMaxPerJB(jb)
          ph(ib,k) = fgps(ib,jb) * &
                 ( b_hybr(k) * (fgtmp(ib,k,jb)-fgtmp(ib,k-1,jb)) ) / delp(ib,k)
       ENDDO
       DO k=1,kMax
          DO ib = 1, ibMaxPerJB(jb)
                IF(1.0_r8 + ((fgzs(ib,jb)/grav )/10000.0_r8) > 1.35_r8)THEN   
                   kh=1.0e-1_r8
                   fgtmp(ib,k,jb) = fgtmp(ib,k,jb) +  kh*ph(ib,k) * fglnps_nabla4(ib,jb)
                ELSE
                   kh=1.0_r8
                   fgtmp(ib,k,jb) = fgtmp(ib,k,jb) +  kh*ph(ib,k) * fglnps_nabla4(ib,jb)
                END IF 
          ENDDO 
       ENDDO
    ENDDO
  END SUBROUTINE End_temp_difus


  SUBROUTINE GlobConservation(jFirst, jLast, jbFirst, jbLast, &
                              jFirst_d, jLast_d)
    INTEGER, INTENT(IN) :: jFirst
    INTEGER, INTENT(IN) :: jLast
    INTEGER, INTENT(IN) :: jbFirst
    INTEGER, INTENT(IN) :: jbLast
    INTEGER, INTENT(IN) :: jFirst_d
    INTEGER, INTENT(IN) :: jLast_d
    INTEGER :: ib, jb, k
    INTEGER :: i, j, ns, j1, ins
    REAL(KIND=r8) :: s
    
    !
    !  Mass Conservation
    !  -----------------
    !
    ! fgpsp = ps global field, arbitray grid
    ! fps = ps global field, 1 latitude per processor
    !
    ! To guarantee numerical reproducibility we need to make sure the
    ! mass summation is done the same way (order) no matter how the
    ! grid is divided. The problem is that the division of grid-cells
    ! between processors is arbitrary, and hence the order and the
    ! number of terms in the summation done by each processor will be
    ! different! To avoid that we redistribute the grid-cells as
    ! one-full-latitude-circle per processor
    !
    !$OMP SINGLE
    IF (.NOT.ALLOCATED(fps)) THEN
       ALLOCATE (fps(imax,myjmax_d))
       IF (nscalars > 0) ALLOCATE (fgs(imax,myjmax_d,nscalars))
       IF (nscalars > 0) ALLOCATE (fg(ibmax,nscalars,jbmax))
       ALLOCATE (displ(0:maxnodes-1))
    ENDIF
    IF (init_globconserv) THEN
       CALL Collect_Gauss(fgps, fps, 1) 
      ELSE
       CALL Collect_Gauss(fgpsp, fps, 1) 
    ENDIF
    !$OMP END SINGLE
    DO j=jFirst_d, jLast_d
       s = 0.0      
       j1 = j - myfirstlat_diag + 1
       DO i=1,imaxperj(j)
          s = s + fps(i,j1)
       ENDDO
       massconsrv(j) = s * cel_area(j)
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE
    displ(0) = 0
    DO i=1,maxnodes-1
       displ(i) = displ(i-1) + nlatsinproc_d(i-1)
    ENDDO
    ! Up to now, this processor has filled only the j-th positions
    ! assigned to him. The remaining positions are now filled with the
    ! values calculated by the other processors. 
    CALL MPI_ALLGATHERV(massconsrv(myfirstlat_diag),myjmax_d,MPI_DOUBLE_PRECISION, &
                        massconsrv,nlatsinproc_d,displ,MPI_DOUBLE_PRECISION, & 
                        MPI_COMM_WORLD, ierr)
    totmass = 0.
    DO j=1,jmax
       totmass = totmass + massconsrv(j)
    ENDDO
    !$OMP END SINGLE

    ! If that is the first time the total mass is calculated then save
    ! the result for future reference, otherwise enforce mass
    ! conservation.
    IF (.NOT. init_globconserv) THEN
       s = total_mass(0) / MAX(totmass,1e-21_r8)
       DO j=jbFirst,jbLast
          DO i=1,ibMaxPerJB(j)
             fgpsp(i,j) = fgpsp(i,j)*s
          ENDDO
       ENDDO
     ELSE
       total_mass(0) = totmass
    ENDIF

    IF (.not. do_globconserv) RETURN
    !
    !  Passive scalars Conservation
    !  ----------------------------
    !$OMP BARRIER
    ins = adr_scalars
    DO ns=1,nscalars
       DO j=jFirst,jLast
          ! Before redistributing the grid points, each processor sum
          ! over the kmax vertical levels, which are independent of
          ! the grid division.
          DO i=myfirstlon(j),mylastlon(j)
             ib = ibPerIJ(i,j)
             jb = jbPerIJ(i,j)
             s = 0.0      
             DO k=1,kmax
                s = s + fgpass_scalars(ib,k,jb,ns,ins)*delb(k)   
             ENDDO
             IF (.not. init_globconserv) THEN
                fg(ib,ns,jb) = s*fgpsp(ib,jb)
               ELSE
                fg(ib,ns,jb) = s*fgps(ib,jb)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE
    !
    CALL Collect_Gauss(fg, fgs, nscalars) 
    !$OMP END SINGLE
    DO ns=1,nscalars
       DO j=jFirst_d, jLast_d
          s = 0.0      
          j1 = j-myfirstlat_diag+1
          DO i=1,imaxperj(j)
             s = s + fgs(i,j1,ns)
          ENDDO
          fconsrv(ns,j) = s * cel_area(j)
       ENDDO
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE

    displ(0) = 0
    DO i=1,maxnodes-1
       displ(i) = displ(i-1) + nlatsinproc_d(i-1)*nscalars
    ENDDO
    CALL MPI_ALLGATHERV(&
         fconsrv(1,myfirstlat_diag), myjmax_d*nscalars, MPI_DOUBLE_PRECISION, &
         fconsrv, nlatsinproc_d(:)*nscalars,displ, MPI_DOUBLE_PRECISION, & 
         MPI_COMM_WORLD, ierr)
    DO ns=1,nscalars
       totmas(ns) = 0.
       DO j=1,jmax
          totmas(ns) = totmas(ns) + fconsrv(ns,j)
       ENDDO
    ENDDO
    !$OMP END SINGLE
    IF (.not. init_globconserv) THEN
       DO ns=1,nscalars
          IF (totmas(ns).ne.0.) THEN
             !$OMP SINGLE
             totmas(ns) = total_mass(ns) / totmas(ns)
             !$OMP END SINGLE
             DO j=jbFirst,jbLast
                DO i=1,ibMaxPerJB(j)
                   fgpass_scalars(i,:,j,ns,ins)=fgpass_scalars(i,:,j,ns,ins) * &
                                                   totmas(ns)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
     ELSE
       !$OMP SINGLE
       DO ns=1,nscalars
          total_mass(ns) = totmas(ns)
       ENDDO
       !$OMP END SINGLE
    ENDIF
       
  END SUBROUTINE GlobConservation

  SUBROUTINE GlobFluxConservation(jFirst, jLast,  jbFirst, jbLast, &
                                 jFirst_d, jLast_d)
    INTEGER, INTENT(IN) :: jFirst
    INTEGER, INTENT(IN) :: jLast
    INTEGER, INTENT(IN) :: jbFirst
    INTEGER, INTENT(IN) :: jbLast
    INTEGER, INTENT(IN) :: jFirst_d
    INTEGER, INTENT(IN) :: jLast_d
    INTEGER :: ib, jb, k
    INTEGER :: i, j, ns, j1, ins
    REAL(KIND=r8) :: s

    IF (nAeros == 0) RETURN

    !$OMP SINGLE
    IF (.NOT.ALLOCATED(fgs_flux)) THEN
       IF (nAeros > 0) ALLOCATE (fgs_flux(imax,myjmax_d,nAeros))
       IF (nAeros > 0) ALLOCATE (fg_flux(ibmax,nAeros   ,jbmax))
       ALLOCATE (displ_flux(0:maxnodes-1))
    ENDIF
    !$OMP END SINGLE
    IF (.not. do_globfluxconserv) RETURN
    !
    !  Passive scalars Conservation
    !  ----------------------------
    DO ns=1,nAeros
       DO j=jFirst,jLast
          DO i=myfirstlon(j),mylastlon(j)
             ib = ibPerIJ(i,j)
             jb = jbPerIJ(i,j)
             fg_flux(ib,ns,jb) =  fgpass_fluxscalars(ib,jb,ns)
          ENDDO
       ENDDO
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE
    CALL Collect_Gauss(fg_flux, fgs_flux, nAeros)
    !$OMP END SINGLE
    DO ns=1,nAeros
       DO j=jFirst_d, jLast_d
          s = 0.0_r8
          j1 = j-myfirstlat_diag+1
          DO i=1,imaxperj(j)
             s = s + fgs_flux(i,j1,ns)
          ENDDO
          fconsrv_flux(ns,j) = s * cel_area(j)
       ENDDO
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE
    displ_flux(0) = 0
    DO i=1,maxnodes-1
       displ_flux(i) = displ_flux(i-1) + nlatsinproc_d(i-1)*nAeros
    ENDDO
    CALL MPI_ALLGATHERV(&
         fconsrv_flux(1,myfirstlat_diag),myjmax_d*nAeros,MPI_DOUBLE_PRECISION, &
         fconsrv_flux,nlatsinproc_d(:)*nAeros,displ_flux,MPI_DOUBLE_PRECISION, &
         MPI_COMM_WORLD, ierr)
    DO ns=1,nAeros
       totflux(ns) = 0.0_r8
       DO j=1,jmax
          totflux(ns) = totflux(ns) + fconsrv_flux(ns,j)
       ENDDO
    ENDDO
    !$OMP END SINGLE

    IF (.not. init_globfluxconserv) THEN
       DO ns=1,nAeros
          IF (totflux(ns).ne.0.) THEN
             !$OMP SINGLE
             totflux(ns) = (total_flux(ns) - totflux(ns))/REAL(imax*jMax)
             !$OMP END SINGLE
             DO j=jbFirst,jbLast
                DO i=1,ibMaxPerJB(j)
                    fgpass_fluxscalars(i,j,ns) =  totflux(ns)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
     ELSE
       !$OMP SINGLE
       DO ns=1,nAeros
          total_flux(ns) = totflux(ns)
       ENDDO
       !$OMP END SINGLE
    ENDIF

  END SUBROUTINE GlobFluxConservation

  !
  !  Vertical Diffusion of Scalar variables
  !  --------------------------------------
  !
  !  Solves diffusion equation implicitly:
  !
  !      dq/dt = - d(<w'q'>)/dz = d[ K dq/dz ]/dz
  !
  !  By finite differencing the equation as:
  !
  !      (q^n+1 - q^n-1)/dt = d[ K dq^n+1/dz ]/dz
  !
  !  This lead to kmax equations that can be written in matrix form
  !
  !                                               (n+1)                    (n-1)
  !      | a(1) c(1)                 0   |   |q(1)|        |q(1)-2 Dt Fs/D1|
  !      | b(2) a(2) c(2)                |   |    |        |               | 
  !      |      b(3) a(3)   c(3)         |   | .  |        | .             | 
  !      |                               | * | .  |      = | .             | 
  !      |                               |   | .  |        | .             | 
  !      |          b(m-1) a(m-1) c(m-1) |   |    |        |               | 
  !      | 0                b(m)   a(m)  |   |q(m)|        |q(m)           | 
  !
  !
! SUBROUTINE Scalardiffusion(iblim, jb, deltat, Kh, tv, gq)

    ! IN/OUT VARIBLES --------------------------------------------------
!   INTEGER, INTENT(IN) :: iblim
!   INTEGER, INTENT(IN) :: jb
!   REAL(KIND=r8), INTENT(IN) :: deltat

!   ! Diffusion coefficient on top each layer (m^2/s)
    ! Note that:
    !  - at k=kmax it should be zero (no flux across the model top)
    !  - there is no k=0 (surface) value because this is acounted 
    !    for in the tendency
!   REAL(KIND=r8), INTENT(IN) :: Kh(ibmax,kmax)

    ! Virtual temperature in the middle of layer (K)
!   REAL(KIND=r8), INTENT(IN) :: tv(ibmax,kmax)

    ! Specific humidity in the middle of layer (kg/kg)
!   REAL(KIND=r8), INTENT(IN) :: gq(ibmax,kmax)

    ! LOCAL VARIABLES --------------------------------------------------

!   INTEGER :: i, k, ns, ins ! counters
!   REAL(KIND=r8) :: s1, s2 ! aux variables

    ! matrix coefficients
!   REAL(KIND=r8) :: a(ibmax,kmax), b(ibmax,kmax), c(ibmax,kmax), dt2

    ! Thermodinamic temperature in the middle of layer (K)
!   REAL(KIND=r8) :: gt(ibmax,kmax)

    !
    !  Compute coefficients
    !  --------------------

!   dt2 = - deltat

    ! Calculate thermodinamic temperature
!   DO k=1,kmax
!      DO i=1,iblim
!         gt(i,k) = tv(i,k)/(1.0_r8+0.608_r8*gq(i,k))
!      ENDDO
!   ENDDO

    ! Calculate b_k coefficients
!   DO i=1,iblim
!      b(i,1) = 0.0_r8
!   ENDDO
!   DO k=2,kmax
!      s2 = dt2 * (grav/gasr)**2 * sl(k) * si(k) / (delcl(k-1)*del(k))
!      DO i=1,iblim
!         b(i,k) = s2 * Kh(i,k-1) * 2._r8 / (gt(i,k-1)+gt(i,k)) / gt(i,k)
!      ENDDO
!   ENDDO

    ! Calculate c_k coefficients 
!   DO k=1,kmax-1
!      s1 = dt2 * (grav/gasr)**2 * sl(k) * si(k+1) / (delcl(k)*del(k))
!      DO i=1,iblim
!         c(i,k) = s1 * Kh(i,k) * 2._r8 / (gt(i,k)+gt(i,k+1)) / gt(i,k)
!      ENDDO
!   ENDDO
!   DO i=1,iblim
!      c(i,kmax) = 0.0_r8
!   ENDDO

    ! Calculate a_k coefficients
    ! We can run 1:kmax because we set b(1)=c(kmax)=0.0
!   DO k=1,kmax
!      DO i=1,iblim
!         a(i,k) = 1.0_r8 - b(i,k) - c(i,k)
!      ENDDO
!   ENDDO

    !
    !  Solve implicit systems
    !  ----------------------
    ! 
    ! First pass, going down and eliminating all b(k)
!   ins = adr_scalars
!   DO k=2,kmax
!      DO i=1,iblim
!         a(i,k) = a(i,k) - c(i,k-1) * b(i,k) / a(i,k-1)
!      ENDDO
!      DO ns=1,nscalars
!         DO i=1,iblim
!            fgpass_scalars(i,k,jb,ns,ins) = fgpass_scalars(i,k,jb,ns,ins) &
!                 - fgpass_scalars(i,k-1,jb,ns,ins) * b(i,k) / a(i,k-1)
!         ENDDO
!      ENDDO
!   ENDDO
    ! The m-th equation is now trivial as only a(kmax) is non-zero
!   DO ns=1,nscalars
!      DO i=1,iblim
!         fgpass_scalars(i,kmax,jb,ns,ins) = fgpass_scalars(i,kmax,jb,ns,ins) &
!                                             /  a(i,kmax)
!      ENDDO
!   ENDDO
    ! Now that we have q(kmax), go upwards calculating q(kmax-1), ... q(1)
!   DO k=kmax-1,1,-1
!      DO ns=1,nscalars
!         DO i=1,iblim
!            fgpass_scalars(i,k,jb,ns,ins) = ( fgpass_scalars(i,k,jb,ns,ins) &
!                          - c(i,k) * fgpass_scalars(i,k+1,jb,ns,ins) ) / a(i,k)
!         ENDDO
!      ENDDO
!   ENDDO
       
! END SUBROUTINE Scalardiffusion


  SUBROUTINE UpdateConserv(jFirst, jLast, &
                          jFirst_d, jLast_d)
    INTEGER, INTENT(IN) :: jFirst
    INTEGER, INTENT(IN) :: jLast
    INTEGER, INTENT(IN) :: jFirst_d
    INTEGER, INTENT(IN) :: jLast_d
    INTEGER :: ib, jb, k
    INTEGER :: i, j, ns, j1, ins
    REAL(KIND=r8) :: s

    !
    !  Passive scalars Conservation
    !  ----------------------------
    ins = adr_scalars
    DO ns=1,nscalars
       DO j=jFirst,jLast
          DO i=myfirstlon(j),mylastlon(j)
             ib = ibPerIJ(i,j)
             jb = jbPerIJ(i,j)
             s = 0.0
             DO k=1,kmax
                s = s + fgpass_scalars(ib,k,jb,ns,ins)*delb(k)
             ENDDO
             fg(ib,ns,jb) = s*fgpsp(ib,jb)
          ENDDO
       ENDDO
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE
    CALL Collect_Gauss(fg, fgs, nscalars)
    !$OMP END SINGLE
    DO ns=1,nscalars
       DO j=jFirst_d, jLast_d
          s = 0.0
          j1 = j-myfirstlat_diag+1!hmjb
          DO i=1,imaxperj(j)
             s = s + fgs(i,j1,ns)
          ENDDO
          fconsrv(ns,j) = s * cel_area(j)
       ENDDO
    ENDDO
    !$OMP BARRIER
    !$OMP SINGLE
    displ(0) = 0
    DO i=1,maxnodes-1
       displ(i) = displ(i-1) + nlatsinproc_d(i-1)*nscalars
    ENDDO
    CALL MPI_ALLGATHERV(&
         fconsrv(1,myfirstlat_diag),myjmax_d*nscalars,MPI_DOUBLE_PRECISION, &
         fconsrv,nlatsinproc_d(:)*nscalars,displ,MPI_DOUBLE_PRECISION, &
         MPI_COMM_WORLD, ierr)
    DO ns=1,nscalars
       total_mass(ns) = 0.
       DO j=1,jmax
          total_mass(ns) = total_mass(ns) + fconsrv(ns,j)
       ENDDO
    ENDDO
    !$OMP END SINGLE

  END SUBROUTINE UpdateConserv

       
END MODULE GridDynamics
