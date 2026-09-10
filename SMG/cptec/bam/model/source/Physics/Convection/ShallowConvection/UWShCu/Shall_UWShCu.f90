
MODULE uwshcu

  !  use cam_history,    only: outfld, addfld, phys_decomp
  !  use error_function, only: erfc
  !  use cam_logfile,    only: iulog
  !  USE ppgrid,         ONLY: nco
  !  use abortutils,     only: endrun

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC init_uwshcu
  PUBLIC compute_uwshcu
  PUBLIC compute_uwshcu_inv
  public fqsatd   ! Function version of vqsatd

  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! 4 byte real
  INTEGER , PARAMETER :: r8 = SELECTED_REAL_KIND(12)    !  8 byte real  
  REAL(r8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number ~ molecules/kmole
  REAL(r8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant ~ J/K/molecule  
  REAL(r8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(r8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air ~ kg/kmole
  REAL(r8),PARAMETER :: SHR_CONST_LATVAP = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_LATICE = 3.337e5_r8      ! latent heat of fusion ~ J/kg
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(r8),PARAMETER :: SHR_CONST_TKFRZ  = 273.16_r8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)


  REAL(r8), PUBLIC, PARAMETER :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights 
  REAL(r8), PUBLIC, PARAMETER :: latvap = shr_const_latvap ! Latent heat of vaporization
  REAL(r8), PUBLIC, PARAMETER :: latice = shr_const_latice ! Latent heat of fusion
  REAL(r8), PUBLIC, PARAMETER :: rh2o    =SHR_CONST_RWV   !! Gas constant for water vapor
  REAL(r8), PUBLIC, PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: tmelt = shr_const_tkfrz   ! Freezing point of water

  REAL(r8)            :: xlv                            !  Latent heat of vaporization
  REAL(r8)            :: xlf                            !  Latent heat of fusion
  REAL(r8)            :: xls                            !  Latent heat of sublimation = xlv + xlf
  REAL(r8)            :: cp                             !  Specific heat of dry air
  REAL(r8)            :: zvir                           !  rh2o/rair - 1
  REAL(r8)            :: r                              !  Gas constant for dry air
  REAL(r8)            :: g                              !  Gravitational constant
  REAL(r8)            :: ep2                            !  mol wgt water vapor / mol wgt dry air 
  REAL(r8)            :: p00                            !  Reference pressure for exner function
  REAL(r8)            :: rovcp                          !  R/cp

  LOGICAL   :: WACCM_MOZART=.TRUE.
  INTERFACE erfc
     MODULE PROCEDURE erfc_r4
     MODULE PROCEDURE derfc
  END INTERFACE

  INTEGER , PARAMETER :: iulog =0 
  !
  ! Data
  !
  INTEGER    plenest  ! length of saturation vapor pressure table
  PARAMETER (plenest=250)
  !
  ! Table of saturation vapor pressure values es from tmin degrees
  ! to tmax+1 degrees k in one degree increments.  ttrice defines the
  ! transition region where es is a combination of ice & water values
  !
  REAL(r8) estbl(plenest)      ! table values of saturation vapor pressure
  REAL(r8) tmin       ! min temperature (K) for table
  REAL(r8) tmax       ! max temperature (K) for table
  REAL(r8) ttrice     ! transition range from es over H2O to es over ice
  REAL(r8) pcf(6)     ! polynomial coeffs -> es transition water to ice
  REAL(r8) epsqs      ! Ratio of h2o to dry air molecular weights 
  REAL(r8) rgasv      ! Gas constant for water vapor
  REAL(r8) hlatf      ! Latent heat of vaporization
  REAL(r8) hlatv      ! Latent heat of fusion
  !REAL(r8) cp         ! specific heat of dry air
  !REAL(r8) tmelt      ! Melting point of water (K)
  LOGICAL icephs  ! false => saturation vapor press over water only

  INTEGER, PARAMETER  :: ppcnst=3
  REAL(r8) :: qmin(3)
  CHARACTER*3, PUBLIC :: cnst_type(ppcnst)          ! wet or dry mixing ratio
  CHARACTER*3, PARAMETER :: mixtype(ppcnst)=(/'wet','wet', 'wet'/) ! mixing ratio type (dry, wet)
  INTEGER   , PARAMETER, PUBLIC :: pcnst  = 3!3PCNST      ! number of advected constituents (including water vapor)
  CHARACTER*3, PARAMETER ::cnst_get_type_byind (pcnst)=(/'wet','wet', 'wet'/)


CONTAINS

  REAL(r8) FUNCTION exnf(pressure)
    REAL(r8), INTENT(in)              :: pressure
    exnf = (pressure/p00)**rovcp
    RETURN
  END FUNCTION exnf

  SUBROUTINE init_uwshcu( kind, xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in, ep2_in )

    !------------------------------------------------------------- ! 
    ! Purpose:                                                     !
    ! Initialize key constants for the shallow convection package. !
    !------------------------------------------------------------- !

    !    use cam_history,   only: outfld, addfld, phys_decomp
    IMPLICIT NONE
    INTEGER , INTENT(in) :: kind       !  kind of reals being passed in
    REAL(r8), INTENT(in) :: xlv_in     !  Latent heat of vaporization
    REAL(r8), INTENT(in) :: xlf_in     !  Latent heat of fusion
    REAL(r8), INTENT(in) :: cp_in      !  Specific heat of dry air
    REAL(r8), INTENT(in) :: zvir_in    !  rh2o/rair - 1
    REAL(r8), INTENT(in) :: r_in       !  Gas constant for dry air
    REAL(r8), INTENT(in) :: g_in       !  Gravitational constant
    REAL(r8), INTENT(in) :: ep2_in     !  mol wgt water vapor / mol wgt dry air 
    REAL(r8),PARAMETER :: tmn  = 173.16_r8          ! Minimum temperature entry in table
    REAL(r8),PARAMETER :: tmx  = 375.16_r8          ! Maximum temperature entry in table
    REAL(r8),PARAMETER :: trice  =  20.00_r8         ! Trans range from es over h2o to es over ice
    LOGICAL ip           ! Ice phase (true or false)
    INTEGER :: ind
    ! ------------------------- !
    ! Internal Output Variables !
    ! ------------------------- !

    !call addfld( 'qtflx_Cu'       , 'kg/m2/s' , pverp , 'A' , 'Convective qt flux'                                  , phys_decomp )
    !call addfld( 'slflx_Cu'       , 'J/m2/s'  , pverp , 'A' , 'Convective sl flux'                                  , phys_decomp )
    !call addfld( 'uflx_Cu'        , 'kg/m/s2' , pverp , 'A' , 'Convective  u flux'                                  , phys_decomp )
    !call addfld( 'vflx_Cu'        , 'kg/m/s2' , pverp , 'A' , 'Convective  v flux'                                  , phys_decomp )

    !call addfld( 'qtten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qt tendency by convection'                           , phys_decomp )
    !call addfld( 'slten_Cu'       , 'J/kg/s'  , pver  , 'A' , 'sl tendency by convection'                           , phys_decomp )
    !call addfld( 'uten_Cu'        , 'm/s2'    , pver  , 'A' , ' u tendency by convection'                           , phys_decomp )
    !call addfld( 'vten_Cu'        , 'm/s2'    , pver  , 'A' , ' v tendency by convection'                           , phys_decomp )
    !call addfld( 'qvten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qv tendency by convection'                           , phys_decomp )
    !call addfld( 'qlten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'ql tendency by convection'                           , phys_decomp )
    !call addfld( 'qiten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qi tendency by convection'                           , phys_decomp )

    !call addfld( 'cbmf_Cu'        , 'kg/m2/s' , 1     , 'A' , 'Cumulus base mass flux'                              , phys_decomp )
    !call addfld( 'ufrcinvbase_Cu' , 'fraction', 1     , 'A' , 'Cumulus fraction at PBL top'                         , phys_decomp ) 
    !call addfld( 'ufrclcl_Cu'     , 'fraction', 1     , 'A' , 'Cumulus fraction at LCL'                             , phys_decomp )
    !call addfld( 'winvbase_Cu'    , 'm/s'     , 1     , 'A' , 'Cumulus vertical velocity at PBL top'                , phys_decomp )
    !call addfld( 'wlcl_Cu'        , 'm/s'     , 1     , 'A' , 'Cumulus vertical velocity at LCL'                    , phys_decomp )
    !call addfld( 'plcl_Cu'        , 'Pa'      , 1     , 'A' , 'LCL of source air'                                   , phys_decomp )
    !call addfld( 'pinv_Cu'        , 'Pa'      , 1     , 'A' , 'PBL top pressure'                                    , phys_decomp )
    !call addfld( 'plfc_Cu'        , 'Pa'      , 1     , 'A' , 'LFC of source air'                                   , phys_decomp )
    !call addfld( 'pbup_Cu'        , 'Pa'      , 1     , 'A' , 'Highest interface level of positive cumulus buoyancy', phys_decomp )
    !call addfld( 'ppen_Cu'        , 'Pa'      , 1     , 'A' , 'Highest level where cumulus w is 0'                  , phys_decomp )
    !call addfld( 'qtsrc_Cu'       , 'kg/kg'   , 1     , 'A' , 'Cumulus source air qt'                               , phys_decomp )
    !call addfld( 'thlsrc_Cu'      , 'K'       , 1     , 'A' , 'Cumulus source air thl'                              , phys_decomp )
    !call addfld( 'thvlsrc_Cu'     , 'K'       , 1     , 'A' , 'Cumulus source air thvl'                             , phys_decomp )
    !call addfld( 'emfkbup_Cu'     , 'kg/m2/s' , 1     , 'A' , 'Penetrative mass flux at kbup'                       , phys_decomp )
    !call addfld( 'cin_Cu'         , 'J/kg'    , 1     , 'A' , 'CIN upto LFC'                                        , phys_decomp )
    !call addfld( 'cinlcl_Cu'      , 'J/kg'    , 1     , 'A' , 'CIN upto LCL'                                        , phys_decomp )
    !call addfld( 'cbmflimit_Cu'   , 'kg/m2/s' , 1     , 'A' , 'cbmflimiter'                                         , phys_decomp ) 
    !call addfld( 'tkeavg_Cu'      , 'm2/s2'   , 1     , 'A' , 'Average tke within PBL for convection scheme'        , phys_decomp ) 
    !call addfld( 'zinv_Cu'        , 'm'       , 1     , 'A' , 'PBL top height'                                      , phys_decomp )
    !call addfld( 'rcwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus LWP+IWP'                                     , phys_decomp )
    !call addfld( 'rlwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus LWP'                                         , phys_decomp )
    !call addfld( 'riwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus IWP'                                         , phys_decomp )
    !call addfld( 'tophgt_Cu'      , 'm'       , 1     , 'A' , 'Cumulus top height'                                  , phys_decomp )

    !call addfld( 'wu_Cu'          , 'm/s'     , pverp , 'A' , 'Convective updraft vertical velocity'                , phys_decomp )
    !call addfld( 'ufrc_Cu'        , 'fraction', pverp , 'A' , 'Convective updraft fractional area'                  , phys_decomp )
    !call addfld( 'qtu_Cu'         , 'kg/kg'   , pverp , 'A' , 'Cumulus updraft qt'                                  , phys_decomp )
    !call addfld( 'thlu_Cu'        , 'K'       , pverp , 'A' , 'Cumulus updraft thl'                                 , phys_decomp )
    !call addfld( 'thvu_Cu'        , 'K'       , pverp , 'A' , 'Cumulus updraft thv'                                 , phys_decomp )
    !call addfld( 'uu_Cu'          , 'm/s'     , pverp , 'A' , 'Cumulus updraft uwnd'                                , phys_decomp )
    !call addfld( 'vu_Cu'          , 'm/s'     , pverp , 'A' , 'Cumulus updraft vwnd'                                , phys_decomp )
    !call addfld( 'qtu_emf_Cu'     , 'kg/kg'   , pverp , 'A' , 'qt of penatratively entrained air'                   , phys_decomp )
    !call addfld( 'thlu_emf_Cu'    , 'K'       , pverp , 'A' , 'thl of penatratively entrained air'                  , phys_decomp )
    !call addfld( 'uu_emf_Cu'      , 'm/s'     , pverp , 'A' , 'uwnd of penatratively entrained air'                 , phys_decomp )
    !call addfld( 'vu_emf_Cu'      , 'm/s'     , pverp , 'A' , 'vwnd of penatratively entrained air'                 , phys_decomp )
    !call addfld( 'umf_Cu'         , 'kg/m2/s' , pverp , 'A' , 'Cumulus updraft mass flux'                           , phys_decomp )
    !call addfld( 'uemf_Cu'        , 'kg/m2/s' , pverp , 'A' , 'Cumulus net ( updraft + entrainment ) mass flux'     , phys_decomp )
    !call addfld( 'qcu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft LWC+IWC'                             , phys_decomp )
    !call addfld( 'qlu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft LWC'                                 , phys_decomp )
    !call addfld( 'qiu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft IWC'                                 , phys_decomp )
    !call addfld( 'cufrc_Cu'       , 'fraction', pver  , 'A' , 'Cumulus cloud fraction'                              , phys_decomp )
    !call addfld( 'fer_Cu'         , '1/m'     , pver  , 'A' , 'Cumulus lateral fractional entrainment rate'         , phys_decomp )
    !call addfld( 'fdr_Cu'         , '1/m'     , pver  , 'A' , 'Cumulus lateral fractional detrainment Rate'         , phys_decomp )

    !call addfld( 'dwten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Expellsion rate of cumulus cloud water to env.'      , phys_decomp )
    !call addfld( 'diten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Expellsion rate of cumulus ice water to env.'        , phys_decomp )
    !call addfld( 'qrten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Production rate of rain by cumulus'                  , phys_decomp )
    !call addfld( 'qsten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Production rate of snow by cumulus'                  , phys_decomp )
    !call addfld( 'flxrain_Cu'     , 'kg/m2/s' , pverp , 'A' , 'Rain flux induced by Cumulus'                        , phys_decomp )
    !call addfld( 'flxsnow_Cu'     , 'kg/m2/s' , pverp , 'A' , 'Snow flux induced by Cumulus'                        , phys_decomp )
    !call addfld( 'ntraprd_Cu'     , 'kg/kg/s' , pver  , 'A' , 'Net production rate of rain by Cumulus'              , phys_decomp )
    !call addfld( 'ntsnprd_Cu'     , 'kg/kg/s' , pver  , 'A' , 'Net production rate of snow by Cumulus'              , phys_decomp )

    !call addfld( 'excessu_Cu'     , 'no'      , pver  , 'A' , 'Updraft saturation excess'                           , phys_decomp )
    !call addfld( 'excess0_Cu'     , 'no'      , pver  , 'A' , 'Environmental saturation excess'                     , phys_decomp )
    !call addfld( 'xc_Cu'          , 'no'      , pver  , 'A' , 'Critical mixing ratio'                               , phys_decomp )
    !call addfld( 'aquad_Cu'       , 'no'      , pver  , 'A' , 'aquad'                                               , phys_decomp )
    !call addfld( 'bquad_Cu'       , 'no'      , pver  , 'A' , 'bquad'                                               , phys_decomp )
    !call addfld( 'cquad_Cu'       , 'no'      , pver  , 'A' , 'cquad'                                               , phys_decomp )
    !call addfld( 'bogbot_Cu'      , 'no'      , pver  , 'A' , 'Cloud buoyancy at the bottom interface'              , phys_decomp )
    !call addfld( 'bogtop_Cu'      , 'no'      , pver  , 'A' , 'Cloud buoyancy at the top interface'                 , phys_decomp )

    !call addfld('exit_UWCu_Cu'    , 'no'      , 1     , 'A' , 'exit_UWCu'                                           , phys_decomp ) 
    !call addfld('exit_conden_Cu'  , 'no'      , 1     , 'A' , 'exit_conden'                                         , phys_decomp ) 
    !call addfld('exit_klclmkx_Cu' , 'no'      , 1     , 'A' , 'exit_klclmkx'                                        , phys_decomp ) 
    !call addfld('exit_klfcmkx_Cu' , 'no'      , 1     , 'A' , 'exit_klfcmkx'                                        , phys_decomp ) 
    !call addfld('exit_ufrc_Cu'    , 'no'      , 1     , 'A' , 'exit_ufrc'                                           , phys_decomp ) 
    !call addfld('exit_wtw_Cu'     , 'no'      , 1     , 'A' , 'exit_wtw'                                            , phys_decomp ) 
    !call addfld('exit_drycore_Cu' , 'no'      , 1     , 'A' , 'exit_drycore'                                        , phys_decomp ) 
    !call addfld('exit_wu_Cu'      , 'no'      , 1     , 'A' , 'exit_wu'                                             , phys_decomp ) 
    !call addfld('exit_cufilter_Cu', 'no'      , 1     , 'A' , 'exit_cufilter'                                       , phys_decomp ) 
    !call addfld('exit_kinv1_Cu'   , 'no'      , 1     , 'A' , 'exit_kinv1'                                          , phys_decomp ) 
    !call addfld('exit_rei_Cu'     , 'no'      , 1     , 'A' , 'exit_rei'                                            , phys_decomp ) 

    !call addfld('limit_shcu_Cu'   , 'no'      , 1     , 'A' , 'limit_shcu'                                          , phys_decomp ) 
    !call addfld('limit_negcon_Cu' , 'no'      , 1     , 'A' , 'limit_negcon'                                        , phys_decomp ) 
    !call addfld('limit_ufrc_Cu'   , 'no'      , 1     , 'A' , 'limit_ufrc'                                          , phys_decomp ) 
    !call addfld('limit_ppen_Cu'   , 'no'      , 1     , 'A' , 'limit_ppen'                                          , phys_decomp ) 
    !call addfld('limit_emf_Cu'    , 'no'      , 1     , 'A' , 'limit_emf'                                           , phys_decomp ) 
    !call addfld('limit_cinlcl_Cu' , 'no'      , 1     , 'A' , 'limit_cinlcl'                                        , phys_decomp ) 
    !call addfld('limit_cin_Cu'    , 'no'      , 1     , 'A' , 'limit_cin'                                           , phys_decomp ) 
    !call addfld('limit_cbmf_Cu'   , 'no'      , 1     , 'A' , 'limit_cbmf'                                          , phys_decomp ) 
    !call addfld('limit_rei_Cu'    , 'no'      , 1     , 'A' , 'limit_rei'                                           , phys_decomp ) 
    !call addfld('ind_delcin_Cu'   , 'no'      , 1     , 'A' , 'ind_delcin'                                          , phys_decomp ) 

    IF( kind .NE. r8 ) THEN
       WRITE(iulog,*) 'wrong KIND of reals passed to init_uwshcu -- exiting.'
       CALL endrun
    ENDIF
    xlv   = xlv_in
    xlf   = xlf_in
    xls   = xlv + xlf
    cp    = cp_in
    zvir  = zvir_in
    r     = r_in
    g     = g_in
    ep2   = ep2_in
    p00   = 1.e5_r8
    rovcp = r/cp

    qmin=1.0e-12_r8

    DO ind=1,ppcnst
       ! set constituent mixing ratio type
       !if ( present(mixtype) )then
       cnst_type(ind) = mixtype(ind) 
       !else
       !   cnst_type(ind) = 'wet'
       !end if
    END DO

    ip    = .TRUE.

    CALL gestbl(tmn     ,tmx     ,trice   ,ip      ,epsilo   , &
         latvap  ,latice  ,rh2o    ,cpair   ,tmelt   )


  END SUBROUTINE init_uwshcu

  SUBROUTINE compute_uwshcu_inv( &
       mix           , &
       mkx           , &
       iend          , &
       ncnst         , &
       dt            , & 
       ps0_inv       , &
       zs0_inv       , &
       p0_inv        , &
       z0_inv        , &
       dp0_inv       , &
       u0_inv        , &
       v0_inv        , &
       qv0_inv       , &
       ql0_inv       , &
       qi0_inv       , &
       t0_inv        , &
       s0_inv        , &
       tr0_inv       , &
       tke_inv       , &
       cldfrct_inv   , &
       concldfrct_inv, &
       pblh          , &
       cush          , & 
       umf_inv       , &
       slflx_inv     , &
       qtflx_inv     , & 
       qvten_inv     , &
       qlten_inv     , &
       qiten_inv     , &
       sten_inv      , &
       uten_inv      , &
       vten_inv      , &
       trten_inv     , &  
       qrten_inv     , &
       qsten_inv     , &
       precip        , &
       snow          , &
       evapc_inv     , &
       cufrc_inv     , &
       qcu_inv       , &
       qlu_inv       , &
       qiu_inv       , &   
       cbmf          , &
       qc_inv        , &
       rliq          , &
       cnt_inv       , &
       cnb_inv       , &
       qsat          , &
       dpdry0_inv       ) 

    IMPLICIT NONE
    INTEGER , INTENT(in)    :: mix
    INTEGER , INTENT(in)    :: mkx
    INTEGER , INTENT(in)    :: iend
    INTEGER , INTENT(in)    :: ncnst
    REAL(r8), INTENT(in)    :: dt                       !  Time step : 2*delta_t [ s ]
    REAL(r8), INTENT(in)    :: ps0_inv   (mix,mkx+1)       !  Environmental pressure at the interfaces [ Pa ]
    REAL(r8), INTENT(in)    :: zs0_inv   (mix,mkx+1)       !  Environmental height at the interfaces   [ m ]
    REAL(r8), INTENT(in)    :: p0_inv    (mix,mkx)          !  Environmental pressure at the layer mid-point [ Pa ]
    REAL(r8), INTENT(in)    :: z0_inv    (mix,mkx)          !  Environmental height at the layer mid-point [ m ]
    REAL(r8), INTENT(in)    :: dp0_inv   (mix,mkx)         !  Environmental layer pressure thickness [ Pa ] > 0.
    REAL(r8), INTENT(in)    :: dpdry0_inv(mix,mkx)      !  Environmental dry layer pressure thickness [ Pa ]
    REAL(r8), INTENT(in)    :: u0_inv    (mix,mkx)          !  Environmental zonal wind [ m/s ]
    REAL(r8), INTENT(in)    :: v0_inv    (mix,mkx)          !  Environmental meridional wind [ m/s ]
    REAL(r8), INTENT(in)    :: qv0_inv   (mix,mkx)         !  Environmental water vapor specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: ql0_inv   (mix,mkx)         !  Environmental liquid water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: qi0_inv   (mix,mkx)         !  Environmental ice specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: t0_inv    (mix,mkx)          !  Environmental temperature [ K ]
    REAL(r8), INTENT(in)    :: s0_inv    (mix,mkx)          !  Environmental dry static energy [ J/kg ]
    REAL(r8), INTENT(in)    :: tr0_inv   (mix,mkx,ncnst)   !  Environmental tracers [ #, kg/kg ]
    REAL(r8), INTENT(in)    :: tke_inv(mix,mkx+1)       !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    REAL(r8), INTENT(in)    :: cldfrct_inv(mix,mkx)     !  Total cloud fraction at the previous time step [ fraction ]
    REAL(r8), INTENT(in)    :: concldfrct_inv(mix,mkx)  !  Total convective ( shallow + deep ) cloud fraction at the previous time step [ fraction ]
    REAL(r8), INTENT(in)    :: pblh(mix)                !  Height of PBL [ m ]
    REAL(r8), INTENT(inout) :: cush(mix)                !  Convective scale height [ m ]
    REAL(r8), INTENT(out)   :: umf_inv(mix,mkx+1)       !  Updraft mass flux at the interfaces [ kg/m2/s ]
    REAL(r8), INTENT(out)   :: qvten_inv(mix,mkx)       !  Tendency of water vapor specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: qlten_inv(mix,mkx)       !  Tendency of liquid water specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: qiten_inv(mix,mkx)       !  Tendency of ice specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: sten_inv(mix,mkx)        !  Tendency of dry static energy [ J/kg/s ]
    REAL(r8), INTENT(out)   :: uten_inv(mix,mkx)        !  Tendency of zonal wind [ m/s2 ]
    REAL(r8), INTENT(out)   :: vten_inv(mix,mkx)        !  Tendency of meridional wind [ m/s2 ]
    REAL(r8), INTENT(out)   :: trten_inv(mix,mkx,ncnst) !  Tendency of tracers [ #/s, kg/kg/s ]
    REAL(r8), INTENT(out)   :: qrten_inv(mix,mkx)       !  Tendency of rain water specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: qsten_inv(mix,mkx)       !  Tendency of snow specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: precip(mix)              !  Precipitation ( rain + snow ) flux at the surface [ m/s ]
    REAL(r8), INTENT(out)   :: snow(mix)                !  Snow flux at the surface [ m/s ]
    REAL(r8), INTENT(out)   :: evapc_inv(mix,mkx)       !  Evaporation of precipitation [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: rliq(mix)                !  Vertical integral of tendency of detrained cloud condensate qc [ m/s ]
    REAL(r8), INTENT(out)   :: slflx_inv(mix,mkx+1)     !  Updraft liquid static energy flux [ J/kg * kg/m2/s ]
    REAL(r8), INTENT(out)   :: qtflx_inv(mix,mkx+1)     !  Updraft total water flux [ kg/kg * kg/m2/s ]
    REAL(r8), INTENT(out)   :: cufrc_inv(mix,mkx)       !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    REAL(r8), INTENT(out)   :: qcu_inv(mix,mkx)         !  Liquid+ice specific humidity within cumulus updraft [ kg/kg ]
    REAL(r8), INTENT(out)   :: qlu_inv(mix,mkx)         !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
    REAL(r8), INTENT(out)   :: qiu_inv(mix,mkx)         !  Ice specific humidity within cumulus updraft [ kg/kg ]
    REAL(r8), INTENT(out)   :: qc_inv(mix,mkx)          !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: cbmf(mix)                !  Cumulus base mass flux [ kg/m2/s ]
    REAL(r8), INTENT(out)   :: cnt_inv(mix)             !  Cumulus top  interface index, cnt = kpen [ no ]
    REAL(r8), INTENT(out)   :: cnb_inv(mix)             !  Cumulus base interface index, cnb = krel - 1 [ no ]
    INTEGER , EXTERNAL      :: qsat                     !  Function pointer to sat vap pressure function

    REAL(r8)                :: ps0(mix,0:mkx)           !  Environmental pressure at the interfaces [ Pa ]
    REAL(r8)                :: zs0(mix,0:mkx)           !  Environmental height at the interfaces   [ m ]
    REAL(r8)                :: p0(mix,mkx)              !  Environmental pressure at the layer mid-point [ Pa ]
    REAL(r8)                :: z0(mix,mkx)              !  Environmental height at the layer mid-point [ m ]
    REAL(r8)                :: dp0(mix,mkx)             !  Environmental layer pressure thickness [ Pa ] > 0.
    REAL(r8)                :: dpdry0(mix,mkx)          !  Environmental dry layer pressure thickness [ Pa ]
    REAL(r8)                :: u0(mix,mkx)              !  Environmental zonal wind [ m/s ]
    REAL(r8)                :: v0(mix,mkx)              !  Environmental meridional wind [ m/s ]
    REAL(r8)                :: tke(mix,0:mkx)           !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    REAL(r8)                :: cldfrct(mix,mkx)         !  Total cloud fraction at the previous time step [ fraction ]
    REAL(r8)                :: concldfrct(mix,mkx)      !  Total convective ( shallow + deep ) cloud fraction at the previous time step [ fraction ]
    REAL(r8)                :: qv0(mix,mkx)             !  Environmental water vapor specific humidity [ kg/kg ]
    REAL(r8)                :: ql0(mix,mkx)             !  Environmental liquid water specific humidity [ kg/kg ]
    REAL(r8)                :: qi0(mix,mkx)             !  Environmental ice specific humidity [ kg/kg ]
    REAL(r8)                :: t0(mix,mkx)              !  Environmental temperature [ K ]
    REAL(r8)                :: s0(mix,mkx)              !  Environmental dry static energy [ J/kg ]
    REAL(r8)                :: tr0(mix,mkx,ncnst)       !  Environmental tracers [ #, kg/kg ]
    REAL(r8)                :: umf(mix,0:mkx)           !  Updraft mass flux at the interfaces [ kg/m2/s ]
    REAL(r8)                :: qvten(mix,mkx)           !  Tendency of water vapor specific humidity [ kg/kg/s ]
    REAL(r8)                :: qlten(mix,mkx)           !  Tendency of liquid water specific humidity [ kg/kg/s ]
    REAL(r8)                :: qiten(mix,mkx)           !  tendency of ice specific humidity [ kg/kg/s ]
    REAL(r8)                :: sten(mix,mkx)            !  Tendency of static energy [ J/kg/s ]
    REAL(r8)                :: uten(mix,mkx)            !  Tendency of zonal wind [ m/s2 ]
    REAL(r8)                :: vten(mix,mkx)            !  Tendency of meridional wind [ m/s2 ]
    REAL(r8)                :: trten(mix,mkx,ncnst)     !  Tendency of tracers [ #/s, kg/kg/s ]
    REAL(r8)                :: qrten(mix,mkx)           !  Tendency of rain water specific humidity [ kg/kg/s ]
    REAL(r8)                :: qsten(mix,mkx)           !  Tendency of snow speficif humidity [ kg/kg/s ]
    REAL(r8)                :: evapc(mix,mkx)           !  Tendency of evaporation of precipitation [ kg/kg/s ]
    REAL(r8)                :: slflx(mix,0:mkx)         !  Updraft liquid static energy flux [ J/kg * kg/m2/s ]
    REAL(r8)                :: qtflx(mix,0:mkx)         !  Updraft total water flux [ kg/kg * kg/m2/s ]
    REAL(r8)                :: cufrc(mix,mkx)           !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    REAL(r8)                :: qcu(mix,mkx)             !  Condensate water specific humidity within cumulus updraft at the layer mid-point [ kg/kg ]
    REAL(r8)                :: qlu(mix,mkx)             !  Liquid water specific humidity within cumulus updraft at the layer mid-point [ kg/kg ]
    REAL(r8)                :: qiu(mix,mkx)             !  Ice specific humidity within cumulus updraft at the layer mid-point [ kg/kg ]
    REAL(r8)                :: qc(mix,mkx)              !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
    REAL(r8)                :: cnt(mix)                 !  Cumulus top  interface index, cnt = kpen [ no ]
    REAL(r8)                :: cnb(mix)                 !  Cumulus base interface index, cnb = krel - 1 [ no ] 
    INTEGER                 :: k                        !  Vertical index for local fields [ no ] 
    INTEGER                 :: k_inv                    !  Vertical index for incoming fields [ no ]
    INTEGER                 :: m                        !  Tracer index [ no ]

    DO k = 1, mkx
       k_inv               = mkx + 1 - k
       p0(:iend,k)         = p0_inv(:iend,k_inv)
       u0(:iend,k)         = u0_inv(:iend,k_inv)
       v0(:iend,k)         = v0_inv(:iend,k_inv)
       z0(:iend,k)         = z0_inv(:iend,k_inv)
       dp0(:iend,k)        = dp0_inv(:iend,k_inv)
       dpdry0(:iend,k)     = dpdry0_inv(:iend,k_inv)
       qv0(:iend,k)        = qv0_inv(:iend,k_inv)
       ql0(:iend,k)        = ql0_inv(:iend,k_inv)
       qi0(:iend,k)        = qi0_inv(:iend,k_inv)
       t0(:iend,k)         = t0_inv(:iend,k_inv)
       s0(:iend,k)         = s0_inv(:iend,k_inv)
       cldfrct(:iend,k)    = cldfrct_inv(:iend,k_inv)
       concldfrct(:iend,k) = concldfrct_inv(:iend,k_inv)
       DO m = 1, ncnst
          tr0(:iend,k,m)   = tr0_inv(:iend,k_inv,m)
       ENDDO
    ENDDO

    DO k = 0, mkx
       k_inv               = mkx + 1 - k
       ps0(:iend,k)        = ps0_inv(:iend,k_inv)
       zs0(:iend,k)        = zs0_inv(:iend,k_inv)
       tke(:iend,k)        = tke_inv(:iend,k_inv)
    END DO

    CALL compute_uwshcu( mix  , mkx    , iend      , ncnst , dt   , &
         ps0  , zs0    , p0        , z0    , dp0  , &
         u0   , v0     , qv0       , ql0   , qi0  , & 
         t0   , s0     , tr0       ,                & 
         tke  , cldfrct, concldfrct, pblh  , cush , & 
         umf  , slflx  , qtflx     ,                &  
         qvten, qlten  , qiten     ,                & 
         sten , uten   , vten      , trten ,        &
         qrten, qsten  , precip    , snow  , evapc, &
         cufrc, qcu    , qlu       , qiu   ,        &
         cbmf , qc     , rliq      ,                &
         cnt  , cnb    , qsat      , dpdry0 )

    ! Reverse cloud top/base interface indices

    cnt_inv(:iend) = mkx + 1 - cnt(:iend)
    cnb_inv(:iend) = mkx + 1 - cnb(:iend)

    DO k = 0, mkx
       k_inv                  = mkx + 1 - k
       umf_inv  (:iend,k_inv) = umf  (:iend,k)       
       slflx_inv(:iend,k_inv) = slflx(:iend,k)     
       qtflx_inv(:iend,k_inv) = qtflx(:iend,k)     
    END DO

    DO k = 1, mkx
       k_inv                         = mkx + 1 - k
       qvten_inv(:iend,k_inv)        = qvten(:iend,k)   
       qlten_inv(:iend,k_inv)        = qlten(:iend,k)   
       qiten_inv(:iend,k_inv)        = qiten(:iend,k)   
       sten_inv (:iend,k_inv)        = sten (:iend,k)    
       uten_inv (:iend,k_inv)        = uten (:iend,k)    
       vten_inv (:iend,k_inv)        = vten (:iend,k)    
       qrten_inv(:iend,k_inv)        = qrten(:iend,k)   
       qsten_inv(:iend,k_inv)        = qsten(:iend,k)   
       evapc_inv(:iend,k_inv)        = evapc(:iend,k)
       cufrc_inv(:iend,k_inv)        = cufrc(:iend,k)   
       qcu_inv  (:iend,k_inv)        = qcu  (:iend,k)     
       qlu_inv  (:iend,k_inv)        = qlu  (:iend,k)     
       qiu_inv  (:iend,k_inv)        = qiu  (:iend,k)     
       qc_inv   (:iend,k_inv)        = qc   (:iend,k)      
       DO m = 1, ncnst
          trten_inv(:iend,k_inv,m)   = trten(:iend,k,m) 
       ENDDO
    ENDDO

  END SUBROUTINE compute_uwshcu_inv


  SUBROUTINE compute_uwshcu( mix      , mkx       , iend         , ncnst    , dt        , &
       ps0_in   , zs0_in    , p0_in        , z0_in    , dp0_in    , &
       u0_in    , v0_in     , qv0_in       , ql0_in   , qi0_in    , &
       t0_in    , s0_in     , tr0_in       ,                        &
       tke_in   , cldfrct_in, concldfrct_in,  pblh_in , cush_inout, & 
       umf_out  , slflx_out , qtflx_out    ,                        &
       qvten_out, qlten_out , qiten_out    ,                        & 
       sten_out , uten_out  , vten_out     , trten_out,             &
       qrten_out, qsten_out , precip_out   , snow_out , evapc_out , &
       cufrc_out, qcu_out   , qlu_out      , qiu_out  ,             &
       cbmf_out , qc_out    , rliq_out     ,                        &
       cnt_out  , cnb_out   , qsat         ,  dpdry0_in )

    ! ------------------------------------------------------------ !
    !                                                              !  
    !  University of Washington Shallow Convection Scheme          !
    !                                                              !
    !  Described in Park and Bretherton. 2008. J. Climate :        !
    !                                                              !
    ! 'The University of Washington shallow convection and         !
    !  moist turbulent schemes and their impact on climate         !
    !  simulations with the Community Atmosphere Model'            !
    !                                                              !
    !  Coded by Sungsu Park. Oct.2005.                             ! 
    !                        May.2008.                             !
    !  For questions, send an email to sungsup@ucar.edu or         ! 
    !                                  sungsu@atmos.washington.edu !
    !                                                              !
    ! ------------------------------------------------------------ !

    !   use cam_history,     only : outfld, addfld, phys_decomp
    !    use constituents,    only : qmin, cnst_get_type_byind, cnst_get_ind
    !#ifdef MODAL_AERO
    !    use modal_aero_data, only : ntot_amode, numptr_amode
    !#endif

    IMPLICIT NONE

    ! ---------------------- !
    ! Input-Output Variables !
    ! ---------------------- !

    INTEGER , INTENT(in)    :: mix
    INTEGER , INTENT(in)    :: mkx
    INTEGER , INTENT(in)    :: iend
    INTEGER , INTENT(in)    :: ncnst
    REAL(r8), INTENT(in)    :: dt                             !  Time step : 2*delta_t [ s ]
    REAL(r8), INTENT(in)    :: ps0_in(mix,0:mkx)              !  Environmental pressure at the interfaces [ Pa ]
    REAL(r8), INTENT(in)    :: zs0_in(mix,0:mkx)              !  Environmental height at the interfaces [ m ]
    REAL(r8), INTENT(in)    :: p0_in(mix,mkx)                 !  Environmental pressure at the layer mid-point [ Pa ]
    REAL(r8), INTENT(in)    :: z0_in(mix,mkx)                 !  Environmental height at the layer mid-point [ m ]
    REAL(r8), INTENT(in)    :: dp0_in(mix,mkx)                !  Environmental layer pressure thickness [ Pa ] > 0.
    REAL(r8), INTENT(in)    :: dpdry0_in(mix,mkx)             !  Environmental dry layer pressure thickness [ Pa ]
    REAL(r8), INTENT(in)    :: u0_in(mix,mkx)                 !  Environmental zonal wind [ m/s ]
    REAL(r8), INTENT(in)    :: v0_in(mix,mkx)                 !  Environmental meridional wind [ m/s ]
    REAL(r8), INTENT(in)    :: qv0_in(mix,mkx)                !  Environmental water vapor specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: ql0_in(mix,mkx)                !  Environmental liquid water specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: qi0_in(mix,mkx)                !  Environmental ice specific humidity [ kg/kg ]
    REAL(r8), INTENT(in)    :: t0_in(mix,mkx)                 !  Environmental temperature [ K ]
    REAL(r8), INTENT(in)    :: s0_in(mix,mkx)                 !  Environmental dry static energy [ J/kg ]
    REAL(r8), INTENT(in)    :: tr0_in(mix,mkx,ncnst)          !  Environmental tracers [ #, kg/kg ]
    REAL(r8), INTENT(in)    :: tke_in(mix,0:mkx)              !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    REAL(r8), INTENT(in)    :: cldfrct_in(mix,mkx)            !  Total cloud fraction at the previous time step [ fraction ]
    REAL(r8), INTENT(in)    :: concldfrct_in(mix,mkx)         !  Total convective cloud fraction at the previous time step [ fraction ]
    REAL(r8), INTENT(in)    :: pblh_in(mix)                   !  Height of PBL [ m ]
    REAL(r8), INTENT(inout) :: cush_inout(mix)                !  Convective scale height [ m ]

    REAL(r8)                   tw0_in(mix,mkx)                !  Wet bulb temperature [ K ]
    REAL(r8)                   qw0_in(mix,mkx)                !  Wet-bulb specific humidity [ kg/kg ]

    REAL(r8), INTENT(out)   :: umf_out(mix,0:mkx)             !  Updraft mass flux at the interfaces [ kg/m2/s ]
    REAL(r8), INTENT(out)   :: qvten_out(mix,mkx)             !  Tendency of water vapor specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: qlten_out(mix,mkx)             !  Tendency of liquid water specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: qiten_out(mix,mkx)             !  Tendency of ice specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: sten_out(mix,mkx)              !  Tendency of dry static energy [ J/kg/s ]
    REAL(r8), INTENT(out)   :: uten_out(mix,mkx)              !  Tendency of zonal wind [ m/s2 ]
    REAL(r8), INTENT(out)   :: vten_out(mix,mkx)              !  Tendency of meridional wind [ m/s2 ]
    REAL(r8), INTENT(out)   :: trten_out(mix,mkx,ncnst)       !  Tendency of tracers [ #/s, kg/kg/s ]
    REAL(r8), INTENT(out)   :: qrten_out(mix,mkx)             !  Tendency of rain water specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: qsten_out(mix,mkx)             !  Tendency of snow specific humidity [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: precip_out(mix)                !  Precipitation ( rain + snow ) rate at surface [ m/s ]
    REAL(r8), INTENT(out)   :: snow_out(mix)                  !  Snow rate at surface [ m/s ]
    REAL(r8), INTENT(out)   :: evapc_out(mix,mkx)             !  Tendency of evaporation of precipitation [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: slflx_out(mix,0:mkx)           !  Updraft/pen.entrainment liquid static energy flux [ J/kg * kg/m2/s ]
    REAL(r8), INTENT(out)   :: qtflx_out(mix,0:mkx)           !  updraft/pen.entrainment total water flux [ kg/kg * kg/m2/s ]
    REAL(r8), INTENT(out)   :: cufrc_out(mix,mkx)             !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    REAL(r8), INTENT(out)   :: qcu_out(mix,mkx)               !  Condensate water specific humidity within cumulus updraft [ kg/kg ]
    REAL(r8), INTENT(out)   :: qlu_out(mix,mkx)               !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
    REAL(r8), INTENT(out)   :: qiu_out(mix,mkx)               !  Ice specific humidity within cumulus updraft [ kg/kg ]
    REAL(r8), INTENT(out)   :: cbmf_out(mix)                  !  Cloud base mass flux [ kg/m2/s ]
    REAL(r8), INTENT(out)   :: qc_out(mix,mkx)                !  Tendency of detrained cumulus condensate into the environment [ kg/kg/s ]
    REAL(r8), INTENT(out)   :: rliq_out(mix)                  !  Vertical integral of qc_out [ m/s ]
    REAL(r8), INTENT(out)   :: cnt_out(mix)                   !  Cumulus top  interface index, cnt = kpen [ no ]
    REAL(r8), INTENT(out)   :: cnb_out(mix)                   !  Cumulus base interface index, cnb = krel - 1 [ no ] 

    !
    ! Internal Output Variables
    !

    INTEGER , EXTERNAL      :: qsat 
    REAL(r8)                   qtten_out(mix,mkx)             !  Tendency of qt [ kg/kg/s ]
    REAL(r8)                   slten_out(mix,mkx)             !  Tendency of sl [ J/kg/s ]
    REAL(r8)                   ufrc_out(mix,0:mkx)            !  Updraft fractional area at the interfaces [ fraction ]
    REAL(r8)                   uflx_out(mix,0:mkx)            !  Updraft/pen.entrainment zonal momentum flux [ m/s/m2/s ]
    REAL(r8)                   vflx_out(mix,0:mkx)            !  Updraft/pen.entrainment meridional momentum flux [ m/s/m2/s ]
    REAL(r8)                   fer_out(mix,mkx)               !  Fractional lateral entrainment rate [ 1/Pa ]
    REAL(r8)                   fdr_out(mix,mkx)               !  Fractional lateral detrainment rate [ 1/Pa ]
    REAL(r8)                   cinh_out(mix)                  !  Convective INhibition upto LFC (CIN) [ J/kg ]
    REAL(r8)                   trflx_out(mix,0:mkx,ncnst)     !  Updraft/pen.entrainment tracer flux [ #/m2/s, kg/kg/m2/s ] 

    ! -------------------------------------------- !
    ! One-dimensional variables at each grid point !
    ! -------------------------------------------- !

    ! 1. Input variables

    REAL(r8)    ps0(0:mkx)                                    !  Environmental pressure at the interfaces [ Pa ]
    REAL(r8)    zs0(0:mkx)                                    !  Environmental height at the interfaces [ m ]
    REAL(r8)    p0(mkx)                                       !  Environmental pressure at the layer mid-point [ Pa ]
    REAL(r8)    z0(mkx)                                       !  Environmental height at the layer mid-point [ m ]
    REAL(r8)    dp0(mkx)                                      !  Environmental layer pressure thickness [ Pa ] > 0.
    REAL(r8)    dpdry0(mkx)                                   !  Environmental dry layer pressure thickness [ Pa ]
    REAL(r8)    u0(mkx)                                       !  Environmental zonal wind [ m/s ]
    REAL(r8)    v0(mkx)                                       !  Environmental meridional wind [ m/s ]
    REAL(r8)    tke(0:mkx)                                    !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
    REAL(r8)    cldfrct(mkx)                                  !  Total cloud fraction at the previous time step [ fraction ]
    REAL(r8)    concldfrct(mkx)                               !  Total convective cloud fraction at the previous time step [ fraction ]
    REAL(r8)    qv0(mkx)                                      !  Environmental water vapor specific humidity [ kg/kg ]
    REAL(r8)    ql0(mkx)                                      !  Environmental liquid water specific humidity [ kg/kg ]
    REAL(r8)    qi0(mkx)                                      !  Environmental ice specific humidity [ kg/kg ]
    REAL(r8)    t0(mkx)                                       !  Environmental temperature [ K ]
    REAL(r8)    s0(mkx)                                       !  Environmental dry static energy [ J/kg ]
    REAL(r8)    pblh                                          !  Height of PBL [ m ]
    REAL(r8)    cush                                          !  Convective scale height [ m ]
    REAL(r8)    tr0(mkx,ncnst)                                !  Environmental tracers [ #, kg/kg ]

    ! 2. Environmental variables directly derived from the input variables

    REAL(r8)    qt0(mkx)                                      !  Environmental total specific humidity [ kg/kg ]
    REAL(r8)    thl0(mkx)                                     !  Environmental liquid potential temperature [ K ]
    REAL(r8)    thvl0(mkx)                                    !  Environmental liquid virtual potential temperature [ K ]
    REAL(r8)    ssqt0(mkx)                                    !  Linear internal slope of environmental total specific humidity [ kg/kg/Pa ]
    REAL(r8)    ssthl0(mkx)                                   !  Linear internal slope of environmental liquid potential temperature [ K/Pa ]
    REAL(r8)    ssu0(mkx)                                     !  Linear internal slope of environmental zonal wind [ m/s/Pa ]
    REAL(r8)    ssv0(mkx)                                     !  Linear internal slope of environmental meridional wind [ m/s/Pa ]
    REAL(r8)    thv0bot(mkx)                                  !  Environmental virtual potential temperature at the bottom of each layer [ K ]
    REAL(r8)    thv0top(mkx)                                  !  Environmental virtual potential temperature at the top of each layer [ K ]
    REAL(r8)    thvl0bot(mkx)                                 !  Environmental liquid virtual potential temperature at the bottom of each layer [ K ]
    REAL(r8)    thvl0top(mkx)                                 !  Environmental liquid virtual potential temperature at the top of each layer [ K ]
    REAL(r8)    exn0(mkx)                                     !  Exner function at the layer mid points [ no ]
    REAL(r8)    exns0(0:mkx)                                  !  Exner function at the interfaces [ no ]
    REAL(r8)    sstr0(mkx,ncnst)                              !  Linear slope of environmental tracers [ #/Pa, kg/kg/Pa ]

    ! 2-1. For preventing negative condensate at the provisional time step

    REAL(r8)    qv0_star(mkx)                                 !  Environmental water vapor specific humidity [ kg/kg ]
    REAL(r8)    ql0_star(mkx)                                 !  Environmental liquid water specific humidity [ kg/kg ]
    REAL(r8)    qi0_star(mkx)                                 !  Environmental ice specific humidity [ kg/kg ]
    REAL(r8)    s0_star(mkx)                                  !  Environmental dry static energy [ J/kg ]

    ! 3. Variables associated with cumulus convection

    REAL(r8)    umf(0:mkx)                                    !  Updraft mass flux at the interfaces [ kg/m2/s ]
    REAL(r8)    emf(0:mkx)                                    !  Penetrative entrainment mass flux at the interfaces [ kg/m2/s ]
    REAL(r8)    qvten(mkx)                                    !  Tendency of water vapor specific humidity [ kg/kg/s ]
    REAL(r8)    qlten(mkx)                                    !  Tendency of liquid water specific humidity [ kg/kg/s ]
    REAL(r8)    qiten(mkx)                                    !  Tendency of ice specific humidity [ kg/kg/s ]
    REAL(r8)    sten(mkx)                                     !  Tendency of dry static energy [ J/kg ]
    REAL(r8)    uten(mkx)                                     !  Tendency of zonal wind [ m/s2 ]
    REAL(r8)    vten(mkx)                                     !  Tendency of meridional wind [ m/s2 ]
    REAL(r8)    qrten(mkx)                                    !  Tendency of rain water specific humidity [ kg/kg/s ]
    REAL(r8)    qsten(mkx)                                    !  Tendency of snow specific humidity [ kg/kg/s ]
    REAL(r8)    precip                                        !  Precipitation rate ( rain + snow) at the surface [ m/s ]
    REAL(r8)    snow                                          !  Snow rate at the surface [ m/s ]
    REAL(r8)    evapc(mkx)                                    !  Tendency of evaporation of precipitation [ kg/kg/s ]
    REAL(r8)    slflx(0:mkx)                                  !  Updraft/pen.entrainment liquid static energy flux [ J/kg * kg/m2/s ]
    REAL(r8)    qtflx(0:mkx)                                  !  Updraft/pen.entrainment total water flux [ kg/kg * kg/m2/s ]
    REAL(r8)    uflx(0:mkx)                                   !  Updraft/pen.entrainment flux of zonal momentum [ m/s/m2/s ]
    REAL(r8)    vflx(0:mkx)                                   !  Updraft/pen.entrainment flux of meridional momentum [ m/s/m2/s ]
    REAL(r8)    cufrc(mkx)                                    !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    REAL(r8)    qcu(mkx)                                      !  Condensate water specific humidity within convective updraft [ kg/kg ]
    REAL(r8)    qlu(mkx)                                      !  Liquid water specific humidity within convective updraft [ kg/kg ]
    REAL(r8)    qiu(mkx)                                      !  Ice specific humidity within convective updraft [ kg/kg ]
    REAL(r8)    dwten(mkx)                                    !  Detrained water tendency from cumulus updraft [ kg/kg/s ]
    REAL(r8)    diten(mkx)                                    !  Detrained ice   tendency from cumulus updraft [ kg/kg/s ]
    REAL(r8)    fer(mkx)                                      !  Fractional lateral entrainment rate [ 1/Pa ]
    REAL(r8)    fdr(mkx)                                      !  Fractional lateral detrainment rate [ 1/Pa ]
    REAL(r8)    uf(mkx)                                       !  Zonal wind at the provisional time step [ m/s ]
    REAL(r8)    vf(mkx)                                       !  Meridional wind at the provisional time step [ m/s ]
    REAL(r8)    qc(mkx)                                       !  Tendency due to detrained 'cloud water + cloud ice' (without rain-snow contribution) [ kg/kg/s ]
    REAL(r8)    qc_l(mkx)                                     !  Tendency due to detrained 'cloud water' (without rain-snow contribution) [ kg/kg/s ]
    REAL(r8)    qc_i(mkx)                                     !  Tendency due to detrained 'cloud ice' (without rain-snow contribution) [ kg/kg/s ]
    REAL(r8)    qc_lm
    REAL(r8)    qc_im
    REAL(r8)    nc_lm
    REAL(r8)    nc_im
    REAL(r8)    ql_emf_kbup
    REAL(r8)    qi_emf_kbup
    REAL(r8)    nl_emf_kbup
    REAL(r8)    ni_emf_kbup
    REAL(r8)    qlten_det
    REAL(r8)    qiten_det
    REAL(r8)    rliq                                          !  Vertical integral of qc [ m/s ] 
    REAL(r8)    cnt                                           !  Cumulus top  interface index, cnt = kpen [ no ]
    REAL(r8)    cnb                                           !  Cumulus base interface index, cnb = krel - 1 [ no ] 
    REAL(r8)    qtten(mkx)                                    !  Tendency of qt [ kg/kg/s ]
    REAL(r8)    slten(mkx)                                    !  Tendency of sl [ J/kg/s ]
    REAL(r8)    ufrc(0:mkx)                                   !  Updraft fractional area [ fraction ]
    REAL(r8)    trten(mkx,ncnst)                              !  Tendency of tracers [ #/s, kg/kg/s ]
    REAL(r8)    trflx(0:mkx,ncnst)                            !  Flux of tracers due to convection [ # * kg/m2/s, kg/kg * kg/m2/s ]
    REAL(r8)    trflx_d(0:mkx)                                !  Adjustive downward flux of tracers to prevent negative tracers
    REAL(r8)    trflx_u(0:mkx)                                !  Adjustive upward   flux of tracers to prevent negative tracers
    REAL(r8)    trmin                                         !  Minimum concentration of tracers allowed
    REAL(r8)    pdelx, dum 

    !----- Variables used for the calculation of condensation sink associated with compensating subsidence
    !      In the current code, this 'sink' tendency is simply set to be zero.

    REAL(r8)    uemf(0:mkx)                                   !  Net updraft mass flux at the interface ( emf + umf ) [ kg/m2/s ]
    REAL(r8)    comsub(mkx)                                   !  Compensating subsidence at the layer mid-point ( unit of mass flux, umf ) [ kg/m2/s ]
    REAL(r8)    qlten_sink(mkx)                               !  Liquid condensate tendency by compensating subsidence/upwelling [ kg/kg/s ]
    REAL(r8)    qiten_sink(mkx)                               !  Ice    condensate tendency by compensating subsidence/upwelling [ kg/kg/s ]
    REAL(r8)    nlten_sink(mkx)                               !  Liquid droplets # tendency by compensating subsidence/upwelling [ kg/kg/s ]
    REAL(r8)    niten_sink(mkx)                               !  Ice    droplets # tendency by compensating subsidence/upwelling [ kg/kg/s ]
    REAL(r8)    thlten_sub, qtten_sub                         !  Tendency of conservative scalars by compensating subsidence/upwelling
    REAL(r8)    qlten_sub, qiten_sub                          !  Tendency of ql0, qi0             by compensating subsidence/upwelling
    REAL(r8)    nlten_sub, niten_sub                          !  Tendency of nl0, ni0             by compensating subsidence/upwelling
    REAL(r8)    thl_prog, qt_prog                             !  Prognosed 'thl, qt' by compensating subsidence/upwelling 

    !----- Variables describing cumulus updraft

    REAL(r8)    wu(0:mkx)                                     !  Updraft vertical velocity at the interface [ m/s ]
    REAL(r8)    thlu(0:mkx)                                   !  Updraft liquid potential temperature at the interface [ K ]
    REAL(r8)    qtu(0:mkx)                                    !  Updraft total specific humidity at the interface [ kg/kg ]
    REAL(r8)    uu(0:mkx)                                     !  Updraft zonal wind at the interface [ m/s ]
    REAL(r8)    vu(0:mkx)                                     !  Updraft meridional wind at the interface [ m/s ]
    REAL(r8)    thvu(0:mkx)                                   !  Updraft virtual potential temperature at the interface [ m/s ]
    REAL(r8)    rei(mkx)                                      !  Updraft fractional mixing rate with the environment [ 1/Pa ]
    REAL(r8)    tru(0:mkx,ncnst)                              !  Updraft tracers [ #, kg/kg ]

    !----- Variables describing conservative scalars of entraining downdrafts  at the 
    !      entraining interfaces, i.e., 'kbup <= k < kpen-1'. At the other interfaces,
    !      belows are simply set to equal to those of updraft for simplicity - but it
    !      does not influence numerical calculation.

    REAL(r8)    thlu_emf(0:mkx)                               !  Penetrative downdraft liquid potential temperature at entraining interfaces [ K ]
    REAL(r8)    qtu_emf(0:mkx)                                !  Penetrative downdraft total water at entraining interfaces [ kg/kg ]
    REAL(r8)    uu_emf(0:mkx)                                 !  Penetrative downdraft zonal wind at entraining interfaces [ m/s ]
    REAL(r8)    vu_emf(0:mkx)                                 !  Penetrative downdraft meridional wind at entraining interfaces [ m/s ]
    REAL(r8)    tru_emf(0:mkx,ncnst)                          !  Penetrative Downdraft tracers at entraining interfaces [ #, kg/kg ]    

    !----- Variables associated with evaporations of convective 'rain' and 'snow'

    REAL(r8)    flxrain(0:mkx)                                !  Downward rain flux at each interface [ kg/m2/s ]
    REAL(r8)    flxsnow(0:mkx)                                !  Downward snow flux at each interface [ kg/m2/s ]
    REAL(r8)    ntraprd(mkx)                                  !  Net production ( production - evaporation +  melting ) rate of rain in each layer [ kg/kg/s ]
    REAL(r8)    ntsnprd(mkx)                                  !  Net production ( production - evaporation + freezing ) rate of snow in each layer [ kg/kg/s ]
    !REAL(r8)    flxsntm                                       !  Downward snow flux at the top of each layer after melting [ kg/m2/s ]
    REAL(r8)    snowmlt                                       !  Snow melting tendency [ kg/kg/s ]
    REAL(r8)    subsat                                        !  Sub-saturation ratio (1-qv/qs) [ no unit ]
    REAL(r8)    evprain                                       !  Evaporation rate of rain [ kg/kg/s ]
    REAL(r8)    evpsnow                                       !  Evaporation rate of snow [ kg/kg/s ]
    REAL(r8)    evplimit                                      !  Limiter of 'evprain + evpsnow' [ kg/kg/s ]
    REAL(r8)    evplimit_rain                                 !  Limiter of 'evprain' [ kg/kg/s ]
    REAL(r8)    evplimit_snow                                 !  Limiter of 'evpsnow' [ kg/kg/s ]
    REAL(r8)    evpint_rain                                   !  Vertically-integrated evaporative flux of rain [ kg/m2/s ]
    REAL(r8)    evpint_snow                                   !  Vertically-integrated evaporative flux of snow [ kg/m2/s ]
    REAL(r8)    kevp                                          !  Evaporative efficiency [ complex unit ]

    !----- Other internal variables

    INTEGER     kk,  k, i, m, kp1, km1
    INTEGER     iter_scaleh, iter_xc
    INTEGER     id_check, status
    INTEGER     klcl                                          !  Layer containing LCL of source air
    INTEGER     kinv                                          !  Inversion layer with PBL top interface as a lower interface
    INTEGER     krel                                          !  Release layer where buoyancy sorting mixing occurs for the first time
    INTEGER     klfc                                          !  LFC layer of cumulus source air
    INTEGER     kbup                                          !  Top layer in which cloud buoyancy is positive at the top interface
    INTEGER     kpen                                          !  Highest layer with positive updraft vertical velocity - top layer cumulus can reach
    LOGICAL     id_exit   
    LOGICAL     forcedCu                                      !  If 'true', cumulus updraft cannot overcome the buoyancy barrier just above the PBL top.
    REAL(r8)    thlsrc, qtsrc, usrc, vsrc, thvlsrc            !  Updraft source air properties
    REAL(r8)    PGFc, uplus, vplus
    REAL(r8)    trsrc(ncnst), tre(ncnst)
    REAL(r8)    plcl, plfc, prel, wrel
    REAL(r8)    frc_rasn
    REAL(r8)    ee2, ud2, wtw, wtwb
    REAL(r8)    xc                                       
    REAL(r8)    cldhgt, scaleh, tscaleh, cridis, rle, rkm
    REAL(r8)    rkfre, sigmaw, epsvarw, tkeavg, dpsum, dpi, thvlmin
    REAL(r8)    thlxsat, qtxsat, thvxsat, x_cu, x_en, thv_x0, thv_x1
    REAL(r8)    thj, qvj, qlj, qij, thvj, tj, thv0j, rho0j, rhos0j, qse 
    REAL(r8)    cin, cinlcl
    REAL(r8)    pe, dpe, exne, thvebot, thle, qte, ue, ve, thlue, qtue, wue
    REAL(r8)    mu, mumin0, mumin1, mumin2, mulcl, mulclstar
    REAL(r8)    cbmf, wcrit, winv, wlcl, ufrcinv, ufrclcl, rmaxfrac
    REAL(r8)    criqc, exql, exqi, rpen, ppen
    REAL(r8)    thl0top, thl0bot, qt0bot, qt0top, thvubot, thvutop
    REAL(r8)    thlu_top, qtu_top, qlu_top, qiu_top, qlu_mid, qiu_mid, exntop
    REAL(r8)    thl0lcl, qt0lcl, thv0lcl, thv0rel, rho0inv, autodet
    REAL(r8)    aquad, bquad, cquad, xc1, xc2, excessu, excess0, xsat, xs1, xs2
    REAL(r8)    bogbot, bogtop, delbog, drage, expfac, rbuoy, rdrag
    REAL(r8)    rcwp, rlwp, riwp, qcubelow, qlubelow, qiubelow
    REAL(r8)    rainflx, snowflx                     
    REAL(r8)    es(1)                               
    REAL(r8)    qs(1)                               
    REAL(r8)    gam(1)                                        !  (L/cp)*dqs/dT
    REAL(r8)    qsat_arg             
    REAL(r8)    xsrc, xmean, xtop, xbot, xflx(0:mkx)
    REAL(r8)    tmp1, tmp2

    !----- Some diagnostic internal output variables

    REAL(r8)  ufrcinvbase_out(mix)                            !  Cumulus updraft fraction at the PBL top [ fraction ]
    REAL(r8)  ufrclcl_out(mix)                                !  Cumulus updraft fraction at the LCL ( or PBL top when LCL is below PBL top ) [ fraction ]
    REAL(r8)  winvbase_out(mix)                               !  Cumulus updraft velocity at the PBL top [ m/s ]
    REAL(r8)  wlcl_out(mix)                                   !  Cumulus updraft velocity at the LCL ( or PBL top when LCL is below PBL top ) [ m/s ]
    REAL(r8)  plcl_out(mix)                                   !  LCL of source air [ Pa ]
    REAL(r8)  pinv_out(mix)                                   !  PBL top pressure [ Pa ]
    REAL(r8)  plfc_out(mix)                                   !  LFC of source air [ Pa ]
    REAL(r8)  pbup_out(mix)                                   !  Highest interface level of positive buoyancy [ Pa ]
    REAL(r8)  ppen_out(mix)                                   !  Highest interface evel where Cu w = 0 [ Pa ]
    REAL(r8)  qtsrc_out(mix)                                  !  Sourse air qt [ kg/kg ]
    REAL(r8)  thlsrc_out(mix)                                 !  Sourse air thl [ K ]
    REAL(r8)  thvlsrc_out(mix)                                !  Sourse air thvl [ K ]
    REAL(r8)  emfkbup_out(mix)                                !  Penetrative downward mass flux at 'kbup' interface [ kg/m2/s ]
    REAL(r8)  cinlclh_out(mix)                                !  Convective INhibition upto LCL (CIN) [ J/kg = m2/s2 ]
    REAL(r8)  tkeavg_out(mix)                                 !  Average tke over the PBL [ m2/s2 ]
    REAL(r8)  cbmflimit_out(mix)                              !  Cloud base mass flux limiter [ kg/m2/s ]
    REAL(r8)  zinv_out(mix)                                   !  PBL top height [ m ]
    REAL(r8)  rcwp_out(mix)                                   !  Layer mean Cumulus LWP+IWP [ kg/m2 ] 
    REAL(r8)  rlwp_out(mix)                                   !  Layer mean Cumulus LWP [ kg/m2 ] 
    REAL(r8)  riwp_out(mix)                                   !  Layer mean Cumulus IWP [ kg/m2 ] 
    REAL(r8)  wu_out(mix,0:mkx)                               !  Updraft vertical velocity ( defined from the release level to 'kpen-1' interface )
    REAL(r8)  qtu_out(mix,0:mkx)                              !  Updraft qt [ kg/kg ]
    REAL(r8)  thlu_out(mix,0:mkx)                             !  Updraft thl [ K ]
    REAL(r8)  thvu_out(mix,0:mkx)                             !  Updraft thv [ K ]
    REAL(r8)  uu_out(mix,0:mkx)                               !  Updraft zonal wind [ m/s ] 
    REAL(r8)  vu_out(mix,0:mkx)                               !  Updraft meridional wind [ m/s ]
    REAL(r8)  qtu_emf_out(mix,0:mkx)                          !  Penetratively entrained qt [ kg/kg ]   
    REAL(r8)  thlu_emf_out(mix,0:mkx)                         !  Penetratively entrained thl [ K ]
    REAL(r8)  uu_emf_out(mix,0:mkx)                           !  Penetratively entrained u [ m/s ]
    REAL(r8)  vu_emf_out(mix,0:mkx)                           !  Penetratively entrained v [ m/s ]
    REAL(r8)  uemf_out(mix,0:mkx)                             !  Net upward mass flux including penetrative entrainment (umf+emf) [ kg/m2/s ]
    REAL(r8)  tru_out(mix,0:mkx,ncnst)                        !  Updraft tracers [ #, kg/kg ]   
    REAL(r8)  tru_emf_out(mix,0:mkx,ncnst)                    !  Penetratively entrained tracers [ #, kg/kg ]

    REAL(r8)  wu_s(0:mkx)                                     !  Same as above but for implicit CIN
    REAL(r8)  qtu_s(0:mkx)
    REAL(r8)  thlu_s(0:mkx)
    REAL(r8)  thvu_s(0:mkx)
    REAL(r8)  uu_s(0:mkx)
    REAL(r8)  vu_s(0:mkx)
    REAL(r8)  qtu_emf_s(0:mkx) 
    REAL(r8)  thlu_emf_s(0:mkx)  
    REAL(r8)  uu_emf_s(0:mkx)   
    REAL(r8)  vu_emf_s(0:mkx)
    REAL(r8)  uemf_s(0:mkx)   
    REAL(r8)  tru_s(0:mkx,ncnst)
    REAL(r8)  tru_emf_s(0:mkx,ncnst)   

    REAL(r8)  dwten_out(mix,mkx)
    REAL(r8)  diten_out(mix,mkx)
    REAL(r8)  flxrain_out(mix,0:mkx)  
    REAL(r8)  flxsnow_out(mix,0:mkx)  
    REAL(r8)  ntraprd_out(mix,mkx)    
    REAL(r8)  ntsnprd_out(mix,mkx)    

    REAL(r8)  dwten_s(mkx)
    REAL(r8)  diten_s(mkx)
    REAL(r8)  flxrain_s(0:mkx)  
    REAL(r8)  flxsnow_s(0:mkx)  
    REAL(r8)  ntraprd_s(mkx)    
    REAL(r8)  ntsnprd_s(mkx)    

    REAL(r8)  excessu_arr_out(mix,mkx)
    REAL(r8)  excessu_arr(mkx) 
    REAL(r8)  excessu_arr_s(mkx)
    REAL(r8)  excess0_arr_out(mix,mkx)
    REAL(r8)  excess0_arr(mkx)
    REAL(r8)  excess0_arr_s(mkx)
    REAL(r8)  xc_arr_out(mix,mkx)
    REAL(r8)  xc_arr(mkx)
    REAL(r8)  xc_arr_s(mkx)
    REAL(r8)  aquad_arr_out(mix,mkx)
    REAL(r8)  aquad_arr(mkx)
    REAL(r8)  aquad_arr_s(mkx)
    REAL(r8)  bquad_arr_out(mix,mkx)
    REAL(r8)  bquad_arr(mkx)
    REAL(r8)  bquad_arr_s(mkx)
    REAL(r8)  cquad_arr_out(mix,mkx) 
    REAL(r8)  cquad_arr(mkx)
    REAL(r8)  cquad_arr_s(mkx)
    REAL(r8)  bogbot_arr_out(mix,mkx)
    REAL(r8)  bogbot_arr(mkx)
    REAL(r8)  bogbot_arr_s(mkx)
    REAL(r8)  bogtop_arr_out(mix,mkx)
    REAL(r8)  bogtop_arr(mkx)
    REAL(r8)  bogtop_arr_s(mkx)

    REAL(r8)  exit_UWCu(mix)
    REAL(r8)  exit_conden(mix)
    REAL(r8)  exit_klclmkx(mix)
    REAL(r8)  exit_klfcmkx(mix)
    REAL(r8)  exit_ufrc(mix)
    REAL(r8)  exit_wtw(mix)
    REAL(r8)  exit_drycore(mix)
    REAL(r8)  exit_wu(mix)
    REAL(r8)  exit_cufilter(mix)
    REAL(r8)  exit_kinv1(mix)
    REAL(r8)  exit_rei(mix)

    REAL(r8)  limit_shcu(mix)
    REAL(r8)  limit_negcon(mix)
    REAL(r8)  limit_ufrc(mix)
    REAL(r8)  limit_ppen(mix)
    REAL(r8)  limit_emf(mix)
    REAL(r8)  limit_cinlcl(mix)
    REAL(r8)  limit_cin(mix)
    REAL(r8)  limit_cbmf(mix)
    REAL(r8)  limit_rei(mix)
    REAL(r8)  ind_delcin(mix)

    REAL(r8) :: ufrcinvbase_s, ufrclcl_s, winvbase_s, wlcl_s, plcl_s, pinv_s, plfc_s, &
         qtsrc_s, thlsrc_s, thvlsrc_s, emfkbup_s, cinlcl_s, pbup_s, ppen_s, cbmflimit_s, &
         tkeavg_s, zinv_s, rcwp_s, rlwp_s, riwp_s 
    REAL(r8) :: ufrcinvbase, winvbase,   emfkbup, cbmflimit  

    !----- Variables for implicit CIN computation

    REAL(r8), DIMENSION(mkx)         :: qv0_s  , ql0_s   , qi0_s   , s0_s    , u0_s    ,           & 
         v0_s   , t0_s    , qt0_s   , qvten_s , &
         qlten_s, qiten_s , qrten_s , qsten_s , sten_s  , evapc_s , &
         uten_s , vten_s  , cufrc_s , qcu_s   , qlu_s   , qiu_s   , &
         fer_s  , fdr_s   , qc_s    , qtten_s , slten_s 
    REAL(r8), DIMENSION(0:mkx)       :: umf_s  , slflx_s , qtflx_s , ufrc_s  , uflx_s , vflx_s
    REAL(r8)                         :: cush_s , precip_s, snow_s  , cin_s   , rliq_s, cbmf_s, cnt_s, cnb_s
    REAL(r8)                         :: cin_i,cin_f,del_CIN,ke,alpha
    REAL(r8)                         :: cinlcl_i,cinlcl_f,del_cinlcl
    INTEGER                          :: iter

    REAL(r8), DIMENSION(mkx,ncnst)   :: tr0_s, trten_s
    REAL(r8), DIMENSION(0:mkx,ncnst) :: trflx_s

    !----- Variables for temporary storages

    REAL(r8), DIMENSION(mkx)         :: qv0_o, ql0_o, qi0_o, t0_o, s0_o, u0_o, v0_o
    REAL(r8), DIMENSION(mkx)         :: qt0_o    , thl0_o   , thvl0_o   ,                         &
         thv0bot_o, thv0top_o, thvl0bot_o, thvl0top_o,             &
         ssthl0_o , ssqt0_o  , ssu0_o    , ssv0_o      
    REAL(r8)                         :: tkeavg_o , thvlmin_o, qtsrc_o  , thvlsrc_o, thlsrc_o ,    &
         usrc_o   , vsrc_o   , plcl_o   , plfc_o   ,               &
         thv0lcl_o 
    INTEGER                          :: kinv_o   , klcl_o   , klfc_o  

    REAL(r8), DIMENSION(mkx,ncnst)   :: tr0_o
    REAL(r8), DIMENSION(mkx,ncnst)   :: sstr0_o  
    REAL(r8), DIMENSION(ncnst)       :: trsrc_o
    INTEGER                          :: ixnumliq, ixnumice

    ! ------------------ !
    !                    !
    ! Define Parameters  !
    !                    !
    ! ------------------ !

    ! ------------------------ !
    ! Iterative xc calculation !
    ! ------------------------ !

    INTEGER , PARAMETER              :: niter_xc = 2

    ! ----------------------------------------------------------- !
    ! Choice of 'CIN = cin' (.true.) or 'CIN = cinlcl' (.false.). !
    ! ----------------------------------------------------------- !

    LOGICAL , PARAMETER              :: use_CINcin = .TRUE.

    ! --------------------------------------------------------------- !
    ! Choice of 'explicit' ( 1 ) or 'implicit' ( 2 )  CIN.            !
    !                                                                 !
    ! When choose 'CIN = cinlcl' above,  it is recommended not to use ! 
    ! implicit CIN, i.e., do 'NOT' choose simultaneously :            !
    !            [ 'use_CINcin=.false. & 'iter_cin=2' ]               !
    ! since 'cinlcl' will be always set to zero whenever LCL is below !
    ! the PBL top interface in the current code. So, averaging cinlcl !
    ! of two iter_cin steps is likely not so good. Except that,   all !
    ! the other combinations of  'use_CINcin'  & 'iter_cin' are OK.   !
    !                                                                 !
    ! Feb 2007, Bundy: Note that use_CINcin = .false. will try to use !
    !           a variable (del_cinlcl) that is not currently set     !
    !                                                                 !
    ! --------------------------------------------------------------- !

    INTEGER , PARAMETER              :: iter_cin = 2

    ! ---------------------------------------------------------------- !
    ! Choice of 'self-detrainment' by negative buoyancy in calculating !
    ! cumulus updraft mass flux at the top interface in each layer.    !
    ! ---------------------------------------------------------------- !

    LOGICAL , PARAMETER              :: use_self_detrain = .FALSE.

    ! --------------------------------------------------------- !
    ! Cumulus momentum flux : turn-on (.true.) or off (.false.) !
    ! --------------------------------------------------------- !

    LOGICAL , PARAMETER              :: use_momenflx = .TRUE.

    ! ----------------------------------------------------------------------------------------- !
    ! Penetrative Entrainment : Cumulative ( .true. , original ) or Non-Cumulative ( .false. )  !
    ! This option ( .false. ) is designed to reduce the sensitivity to the vertical resolution. !
    ! ----------------------------------------------------------------------------------------- !

    LOGICAL , PARAMETER              :: use_cumpenent = .TRUE.

    ! --------------------------------------------------------------------------------------------------------------- !
    ! Computation of the grid-mean condensate tendency.                                                               !
    !     use_expconten = .true.  : explcitly compute tendency by condensate detrainment and compensating subsidence  !
    !     use_expconten = .false. : use the original proportional condensate tendency equation. ( original )          !
    ! --------------------------------------------------------------------------------------------------------------- !

    LOGICAL , PARAMETER              :: use_expconten = .TRUE.

    ! --------------------------------------------------------------------------------------------------------------- !
    ! Treatment of reserved condensate                                                                                !
    !     use_unicondet = .true.  : detrain condensate uniformly over the environment ( original )                    !
    !     use_unicondet = .false. : detrain condensate into the pre-existing stratus                                  !
    ! --------------------------------------------------------------------------------------------------------------- !

    LOGICAL , PARAMETER              :: use_unicondet = .FALSE.

    ! ----------------------- !
    ! For lateral entrainment !
    ! ----------------------- !

    PARAMETER (rle = 0.1_r8)         !  For critical stopping distance for lateral entrainment [no unit]
    !   parameter (rkm = 16.0_r8)        !  Determine the amount of air that is involved in buoyancy-sorting [no unit] 
    PARAMETER (rkm = 14.0_r8)        !  Determine the amount of air that is involved in buoyancy-sorting [no unit]

    PARAMETER (rpen = 10.0_r8)       !  For penetrative entrainment efficiency
    PARAMETER (rkfre = 1.0_r8)       !  Vertical velocity variance as fraction of  tke. 
    PARAMETER (rmaxfrac = 0.10_r8)   !  Maximum allowable 'core' updraft fraction
    PARAMETER (mumin1 = 0.906_r8)    !  Normalized CIN ('mu') corresponding to 'rmaxfrac' at the PBL top
    !  obtaind by inverting 'rmaxfrac = 0.5*erfc(mumin1)'.
    !  [ rmaxfrac:mumin1 ] = [ 0.05:1.163, 0.075:1.018, 0.1:0.906, 0.15:0.733, 0.2:0.595, 0.25:0.477 ] 
    PARAMETER (rbuoy = 1.0_r8)       !  For nonhydrostatic pressure effects on updraft [no unit]
    PARAMETER (rdrag = 1.0_r8)       !  Drag coefficient [no unit]

    PARAMETER (epsvarw = 5.e-4_r8)   !  Variance of w at PBL top by meso-scale component [m2/s2]          
    PARAMETER (PGFc = 0.7_r8)        !  This is used for calculating vertical variations cumulus  
    !  'u' & 'v' by horizontal PGF during upward motion [no unit]

    ! ---------------------------------------- !
    ! Bulk microphysics controlling parameters !
    ! --------------------------------------------------------------------------- ! 
    ! criqc    : Maximum condensate that can be hold by cumulus updraft [kg/kg]   !
    ! frc_rasn : Fraction of precipitable condensate in the expelled cloud water  !
    !            from cumulus updraft. The remaining fraction ('1-frc_rasn')  is  !
    !            'suspended condensate'.                                          !
    !                0 : all expelled condensate is 'suspended condensate'        ! 
    !                1 : all expelled condensate is 'precipitable condensate'     !
    ! kevp     : Evaporative efficiency                                           !
    ! noevap_krelkpen : No evaporation from 'krel' to 'kpen' layers               ! 
    ! --------------------------------------------------------------------------- !    

    PARAMETER ( criqc    = 0.7e-3_r8 ) 
    PARAMETER ( frc_rasn = 1.0_r8    )
    PARAMETER ( kevp     = 2.e-6_r8  )
    LOGICAL, PARAMETER :: noevap_krelkpen = .FALSE.

    !------------------------!
    !                        !
    ! Start Main Calculation !
    !                        !
    !------------------------!

    !call cnst_get_ind( 'NUMLIQ', ixnumliq )
    !call cnst_get_ind( 'NUMICE', ixnumice )
    ixnumliq=2
    ixnumice=3
    ! ------------------------------------------------------- !
    ! Initialize output variables defined for all grid points !
    ! ------------------------------------------------------- !
    del_cinlcl=0.0_r8
    umf_out(:iend,0:mkx)         = 0.0_r8
    slflx_out(:iend,0:mkx)       = 0.0_r8
    qtflx_out(:iend,0:mkx)       = 0.0_r8
    qvten_out(:iend,:mkx)        = 0.0_r8
    qlten_out(:iend,:mkx)        = 0.0_r8
    qiten_out(:iend,:mkx)        = 0.0_r8
    sten_out(:iend,:mkx)         = 0.0_r8
    uten_out(:iend,:mkx)         = 0.0_r8
    vten_out(:iend,:mkx)         = 0.0_r8
    qrten_out(:iend,:mkx)        = 0.0_r8
    qsten_out(:iend,:mkx)        = 0.0_r8
    precip_out(:iend)            = 0.0_r8
    snow_out(:iend)              = 0.0_r8
    evapc_out(:iend,:mkx)        = 0.0_r8
    cufrc_out(:iend,:mkx)        = 0.0_r8
    qcu_out(:iend,:mkx)          = 0.0_r8
    qlu_out(:iend,:mkx)          = 0.0_r8
    qiu_out(:iend,:mkx)          = 0.0_r8
    fer_out(:iend,:mkx)          = 0.0_r8
    fdr_out(:iend,:mkx)          = 0.0_r8
    cinh_out(:iend)              = -1.0_r8
    cinlclh_out(:iend)           = -1.0_r8
    cbmf_out(:iend)              = 0.0_r8
    qc_out(:iend,:mkx)           = 0.0_r8
    rliq_out(:iend)              = 0.0_r8
    cnt_out(:iend)               = REAL(mkx, r8)
    cnb_out(:iend)               = 0.0_r8
    qtten_out(:iend,:mkx)        = 0.0_r8
    slten_out(:iend,:mkx)        = 0.0_r8
    ufrc_out(:iend,0:mkx)        = 0.0_r8

    uflx_out(:iend,0:mkx)        = 0.0_r8
    vflx_out(:iend,0:mkx)        = 0.0_r8

    trten_out(:iend,:mkx,:ncnst) = 0.0_r8
    trflx_out(:iend,0:mkx,:ncnst)= 0.0_r8

    ufrcinvbase_out(:iend)       = 0.0_r8
    ufrclcl_out(:iend)           = 0.0_r8
    winvbase_out(:iend)          = 0.0_r8
    wlcl_out(:iend)              = 0.0_r8
    plcl_out(:iend)              = 0.0_r8
    pinv_out(:iend)              = 0.0_r8
    plfc_out(:iend)              = 0.0_r8
    pbup_out(:iend)              = 0.0_r8
    ppen_out(:iend)              = 0.0_r8
    qtsrc_out(:iend)             = 0.0_r8
    thlsrc_out(:iend)            = 0.0_r8
    thvlsrc_out(:iend)           = 0.0_r8
    emfkbup_out(:iend)           = 0.0_r8
    cbmflimit_out(:iend)         = 0.0_r8
    tkeavg_out(:iend)            = 0.0_r8
    zinv_out(:iend)              = 0.0_r8
    rcwp_out(:iend)              = 0.0_r8
    rlwp_out(:iend)              = 0.0_r8
    riwp_out(:iend)              = 0.0_r8

    wu_out(:iend,0:mkx)          = 0.0_r8
    qtu_out(:iend,0:mkx)         = 0.0_r8
    thlu_out(:iend,0:mkx)        = 0.0_r8
    thvu_out(:iend,0:mkx)        = 0.0_r8
    uu_out(:iend,0:mkx)          = 0.0_r8
    vu_out(:iend,0:mkx)          = 0.0_r8
    qtu_emf_out(:iend,0:mkx)     = 0.0_r8
    thlu_emf_out(:iend,0:mkx)    = 0.0_r8
    uu_emf_out(:iend,0:mkx)      = 0.0_r8
    vu_emf_out(:iend,0:mkx)      = 0.0_r8
    uemf_out(:iend,0:mkx)        = 0.0_r8

    tru_out(:iend,0:mkx,:ncnst)     = 0.0_r8
    tru_emf_out(:iend,0:mkx,:ncnst) = 0.0_r8

    dwten_out(:iend,:mkx)        = 0.0_r8
    diten_out(:iend,:mkx)        = 0.0_r8
    flxrain_out(:iend,0:mkx)     = 0.0_r8  
    flxsnow_out(:iend,0:mkx)     = 0.0_r8
    ntraprd_out(:iend,mkx)       = 0.0_r8
    ntsnprd_out(:iend,mkx)       = 0.0_r8

    excessu_arr_out(:iend,:mkx)  = 0.0_r8
    excess0_arr_out(:iend,:mkx)  = 0.0_r8
    xc_arr_out(:iend,:mkx)       = 0.0_r8
    aquad_arr_out(:iend,:mkx)    = 0.0_r8
    bquad_arr_out(:iend,:mkx)    = 0.0_r8
    cquad_arr_out(:iend,:mkx)    = 0.0_r8
    bogbot_arr_out(:iend,:mkx)   = 0.0_r8
    bogtop_arr_out(:iend,:mkx)   = 0.0_r8

    exit_UWCu(:iend)             = 0.0_r8 
    exit_conden(:iend)           = 0.0_r8 
    exit_klclmkx(:iend)          = 0.0_r8 
    exit_klfcmkx(:iend)          = 0.0_r8 
    exit_ufrc(:iend)             = 0.0_r8 
    exit_wtw(:iend)              = 0.0_r8 
    exit_drycore(:iend)          = 0.0_r8 
    exit_wu(:iend)               = 0.0_r8 
    exit_cufilter(:iend)         = 0.0_r8 
    exit_kinv1(:iend)            = 0.0_r8 
    exit_rei(:iend)              = 0.0_r8 

    limit_shcu(:iend)            = 0.0_r8 
    limit_negcon(:iend)          = 0.0_r8 
    limit_ufrc(:iend)            = 0.0_r8
    limit_ppen(:iend)            = 0.0_r8
    limit_emf(:iend)             = 0.0_r8
    limit_cinlcl(:iend)          = 0.0_r8
    limit_cin(:iend)             = 0.0_r8
    limit_cbmf(:iend)            = 0.0_r8
    limit_rei(:iend)             = 0.0_r8

    ind_delcin(:iend)            = 0.0_r8

    !--------------------------------------------------------------!
    !                                                              !
    ! Start the column i loop where i is a horozontal column index !
    !                                                              !
    !--------------------------------------------------------------!

    ! Compute wet-bulb temperature and specific humidity
    ! for treating evaporation of precipitation.

    CALL findsp( iend,mkx, qv0_in, t0_in, p0_in, tw0_in, qw0_in )

    DO i = 1, iend                                      

       id_exit = .FALSE.

       ! -------------------------------------------- !
       ! Define 1D input variables at each grid point !
       ! -------------------------------------------- !

       ps0(0:mkx)       = ps0_in(i,0:mkx)
       zs0(0:mkx)       = zs0_in(i,0:mkx)
       p0(:mkx)         = p0_in(i,:mkx)
       z0(:mkx)         = z0_in(i,:mkx)
       dp0(:mkx)        = dp0_in(i,:mkx)
       dpdry0(:mkx)     = dpdry0_in(i,:mkx)
       u0(:mkx)         = u0_in(i,:mkx)
       v0(:mkx)         = v0_in(i,:mkx)
       qv0(:mkx)        = qv0_in(i,:mkx)
       ql0(:mkx)        = ql0_in(i,:mkx)
       qi0(:mkx)        = qi0_in(i,:mkx)
       t0(:mkx)         = t0_in(i,:mkx)
       s0(:mkx)         = s0_in(i,:mkx)
       tke(0:mkx)       = tke_in(i,0:mkx)
       cldfrct(:mkx)    = cldfrct_in(i,:mkx)
       concldfrct(:mkx) = concldfrct_in(i,:mkx)
       pblh             = pblh_in(i)
       cush             = cush_inout(i)
       DO m = 1, ncnst
          tr0(:mkx,m)   = tr0_in(i,:mkx,m)
       ENDDO

       ! --------------------------------------------------------- !
       ! Compute other basic thermodynamic variables directly from ! 
       ! the input variables at each grid point                    !
       ! --------------------------------------------------------- !

       !----- 1. Compute internal environmental variables

       exn0(:mkx)   = (p0(:mkx)/p00)**rovcp
       exns0(0:mkx) = (ps0(0:mkx)/p00)**rovcp
       qt0(:mkx)    = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
       thl0(:mkx)   = (t0(:mkx) - xlv*ql0(:mkx)/cp - xls*qi0(:mkx)/cp)/exn0(:mkx)
       thvl0(:mkx)  = (1._r8 + zvir*qt0(:mkx))*thl0(:mkx)

       !----- 2. Compute slopes of environmental variables in each layer
       !         Dimension of ssthl0(:mkx) is implicit.

       ssthl0       = slope(mkx,thl0,p0) 
       ssqt0        = slope(mkx,qt0 ,p0)
       ssu0         = slope(mkx,u0  ,p0)
       ssv0         = slope(mkx,v0  ,p0)
       DO m = 1, ncnst
          sstr0(:mkx,m) = slope(mkx,tr0(:mkx,m),p0)
       ENDDO

       !----- 3. Compute "thv0" and "thvl0" at the top/bottom interfaces in each layer
       !         There are computed from the reconstructed thl, qt at the top/bottom.

       DO k = 1, mkx

          thl0bot = thl0(k) + ssthl0(k)*(ps0(k-1) - p0(k))
          qt0bot  = qt0(k)  + ssqt0(k) *(ps0(k-1) - p0(k))
          CALL conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check,qsat)
          IF( id_check .EQ. 1 ) THEN
             exit_conden(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          END IF
          thv0bot(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
          thvl0bot(k) = thl0bot*(1._r8 + zvir*qt0bot)

          thl0top = thl0(k) + ssthl0(k)*(ps0(k) - p0(k))
          qt0top  =  qt0(k) + ssqt0(k) *(ps0(k) - p0(k))
          CALL conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check,qsat)
          IF( id_check .EQ. 1 ) THEN
             exit_conden(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          END IF
          thv0top(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
          thvl0top(k) = thl0top*(1._r8 + zvir*qt0top)

       END DO

       ! ------------------------------------------------------------ !
       ! Save input and related environmental thermodynamic variables !
       ! for use at "iter_cin=2" when "del_CIN >= 0"                  !
       ! ------------------------------------------------------------ !

       qv0_o(:mkx)          = qv0(:mkx)
       ql0_o(:mkx)          = ql0(:mkx)
       qi0_o(:mkx)          = qi0(:mkx)
       t0_o(:mkx)           = t0(:mkx)
       s0_o(:mkx)           = s0(:mkx)
       u0_o(:mkx)           = u0(:mkx)
       v0_o(:mkx)           = v0(:mkx)
       qt0_o(:mkx)          = qt0(:mkx)
       thl0_o(:mkx)         = thl0(:mkx)
       thvl0_o(:mkx)        = thvl0(:mkx)
       ssthl0_o(:mkx)       = ssthl0(:mkx)
       ssqt0_o(:mkx)        = ssqt0(:mkx)
       thv0bot_o(:mkx)      = thv0bot(:mkx)
       thv0top_o(:mkx)      = thv0top(:mkx)
       thvl0bot_o(:mkx)     = thvl0bot(:mkx)
       thvl0top_o(:mkx)     = thvl0top(:mkx)
       ssu0_o(:mkx)         = ssu0(:mkx) 
       ssv0_o(:mkx)         = ssv0(:mkx) 
       DO m = 1, ncnst
          tr0_o(:mkx,m)     = tr0(:mkx,m)
          sstr0_o(:mkx,m)   = sstr0(:mkx,m)
       ENDDO

       ! ---------------------------------------------- !
       ! Initialize output variables at each grid point !
       ! ---------------------------------------------- !

       umf(0:mkx)          = 0.0_r8
       emf(0:mkx)          = 0.0_r8
       slflx(0:mkx)        = 0.0_r8
       qtflx(0:mkx)        = 0.0_r8
       uflx(0:mkx)         = 0.0_r8
       vflx(0:mkx)         = 0.0_r8
       qvten(:mkx)         = 0.0_r8
       qlten(:mkx)         = 0.0_r8
       qiten(:mkx)         = 0.0_r8
       sten(:mkx)          = 0.0_r8
       uten(:mkx)          = 0.0_r8
       vten(:mkx)          = 0.0_r8
       qrten(:mkx)         = 0.0_r8
       qsten(:mkx)         = 0.0_r8
       dwten(:mkx)         = 0.0_r8
       diten(:mkx)         = 0.0_r8
       precip              = 0.0_r8
       snow                = 0.0_r8
       evapc(:mkx)         = 0.0_r8
       cufrc(:mkx)         = 0.0_r8
       qcu(:mkx)           = 0.0_r8
       qlu(:mkx)           = 0.0_r8
       qiu(:mkx)           = 0.0_r8
       fer(:mkx)           = 0.0_r8
       fdr(:mkx)           = 0.0_r8
       cin                 = 0.0_r8
       cbmf                = 0.0_r8
       qc(:mkx)            = 0.0_r8
       qc_l(:mkx)          = 0.0_r8
       qc_i(:mkx)          = 0.0_r8
       rliq                = 0.0_r8
       cnt                 = REAL(mkx, r8)
       cnb                 = 0.0_r8
       qtten(:mkx)         = 0.0_r8
       slten(:mkx)         = 0.0_r8   
       ufrc(0:mkx)         = 0.0_r8  

       thlu(0:mkx)         = 0.0_r8
       qtu(0:mkx)          = 0.0_r8
       uu(0:mkx)           = 0.0_r8
       vu(0:mkx)           = 0.0_r8
       wu(0:mkx)           = 0.0_r8
       thvu(0:mkx)         = 0.0_r8
       thlu_emf(0:mkx)     = 0.0_r8
       qtu_emf(0:mkx)      = 0.0_r8
       uu_emf(0:mkx)       = 0.0_r8
       vu_emf(0:mkx)       = 0.0_r8

       ufrcinvbase         = 0.0_r8
       ufrclcl             = 0.0_r8
       winvbase            = 0.0_r8
       wlcl                = 0.0_r8
       emfkbup             = 0.0_r8 
       cbmflimit           = 0.0_r8
       excessu_arr(:mkx)   = 0.0_r8
       excess0_arr(:mkx)   = 0.0_r8
       xc_arr(:mkx)        = 0.0_r8
       aquad_arr(:mkx)     = 0.0_r8
       bquad_arr(:mkx)     = 0.0_r8
       cquad_arr(:mkx)     = 0.0_r8
       bogbot_arr(:mkx)    = 0.0_r8
       bogtop_arr(:mkx)    = 0.0_r8

       uemf(0:mkx)         = 0.0_r8
       comsub(:mkx)        = 0.0_r8
       qlten_sink(:mkx)    = 0.0_r8
       qiten_sink(:mkx)    = 0.0_r8 
       nlten_sink(:mkx)    = 0.0_r8
       niten_sink(:mkx)    = 0.0_r8 

       DO m = 1, ncnst
          trflx(0:mkx,m)   = 0.0_r8
          trten(:mkx,m)    = 0.0_r8
          tru(0:mkx,m)     = 0.0_r8
          tru_emf(0:mkx,m) = 0.0_r8
       ENDDO

       !-----------------------------------------------! 
       ! Below 'iter' loop is for implicit CIN closure !
       !-----------------------------------------------!

       ! ----------------------------------------------------------------------------- ! 
       ! It is important to note that this iterative cin loop is located at the outest !
       ! shell of the code. Thus, source air properties can also be changed during the !
       ! iterative cin calculation, because cumulus convection induces non-zero fluxes !
       ! even at interfaces below PBL top height through 'fluxbelowinv' subroutine.    !
       ! ----------------------------------------------------------------------------- !

       DO iter = 1, iter_cin

          ! ---------------------------------------------------------------------- ! 
          ! Cumulus scale height                                                   ! 
          ! In contrast to the premitive code, cumulus scale height is iteratively !
          ! calculated at each time step, and at each iterative cin step.          !
          ! It is not clear whether I should locate below two lines within or  out !
          ! of the iterative cin loop.                                             !
          ! ---------------------------------------------------------------------- !

          tscaleh = cush                        
          cush    = -1._r8

          ! ----------------------------------------------------------------------- !
          ! Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer !
          ! index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the !
          ! layer index having PBLH as a lower interface.                           !
          ! In the previous code, I set the lower limit of 'kinv' by 2  in order to !
          ! be consistent with the other parts of the code. However in the modified !
          ! code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program !
          ! without performing cumulus convection. This new approach seems to be    !
          ! more reasonable: if PBL height is within 'kinv=1' layer, surface is STL !
          ! interface (bflxs <= 0) and interface just above the surface should be   !
          ! either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this !
          ! interface is identified as a base external interface of upperlying CL.  !
          ! Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case !
          ! it is reasonable to assume that cumulus convection does not happen.     !
          ! When these is SBCL, PBL height from the PBL scheme is likely to be very !
          ! close at 'kinv-1' interface, but not exactly, since 'zi' information is !
          ! changed between two model time steps. In order to ensure correct identi !
          ! fication of 'kinv' for general case including SBCL, I imposed an offset !
          ! of 5 [m] in the below 'kinv' finding block.                             !
          ! ----------------------------------------------------------------------- !

          DO k = mkx - 1, 1, -1 
             IF( (pblh + 5._r8 - zs0(k))*(pblh + 5._r8 - zs0(k+1)) .LT. 0._r8 ) THEN
                kinv = k + 1 
                go to 15
             ENDIF
          END DO
          kinv = 1
15        CONTINUE    

          IF( kinv .LE. 1 ) THEN          
             exit_kinv1(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          ENDIF
          ! From here, it must be 'kinv >= 2'.

          ! -------------------------------------------------------------------------- !
          ! Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL !
          ! In the current code, 'tkeavg' is obtained by averaging all interfacial TKE !
          ! within the PBL. However, in order to be conceptually consistent with   PBL !
          ! scheme, 'tkeavg' should be calculated by considering surface buoyancy flux.!
          ! If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE !
          ! should be included in calculating 'tkeavg', while if bflxs <= 0,   surface !
          ! interfacial TKE should not be included in calculating 'tkeavg'.   I should !
          ! modify the code when 'bflxs' is available as an input of cumulus scheme.   !
          ! 'thvlmin' is a minimum 'thvl' within PBL obtained by comparing top &  base !
          ! interface values of 'thvl' in each layers within the PBL.                  !
          ! -------------------------------------------------------------------------- !

          dpsum    = 0._r8
          tkeavg   = 0._r8
          thvlmin  = 1000._r8
          DO k = 0, kinv - 1   ! Here, 'k' is an interfacial layer index.  
             IF( k .EQ. 0 ) THEN
                dpi = ps0(0) - p0(1)
             ELSEIF( k .EQ. (kinv-1) ) THEN 
                dpi = p0(kinv-1) - ps0(kinv-1)
             ELSE
                dpi = p0(k) - p0(k+1)
             ENDIF
             dpsum  = dpsum  + dpi  
             tkeavg = tkeavg + dpi*tke(k) 
             IF( k .NE. 0 ) thvlmin = MIN(thvlmin,MIN(thvl0bot(k),thvl0top(k)))
          END DO
          tkeavg  = tkeavg/dpsum

          ! ------------------------------------------------------------------ !
          ! Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc !
          ! Note that 'thlsrc' was con-cocked using 'thvlsrc' and 'qtsrc'.     !
          ! 'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc' !
          ! is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined !
          ! as the values just below the PBL top interface.                    !
          ! ------------------------------------------------------------------ !

          qtsrc   = qt0(1)                     
          thvlsrc = thvlmin 
          thlsrc  = thvlsrc / ( 1._r8 + zvir * qtsrc )  
          usrc    = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )             
          vsrc    = v0(kinv-1) + ssv0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )             
          DO m = 1, ncnst
             trsrc(m) = tr0(1,m)
          ENDDO

          ! ------------------------------------------------------------------ !
          ! Find LCL of the source air and a layer index containing LCL (klcl) !
          ! When the LCL is exactly at the interface, 'klcl' is a layer index  ! 
          ! having 'plcl' as the lower interface similar to the 'kinv' case.   !
          ! In the previous code, I assumed that if LCL is located within the  !
          ! lowest model layer ( 1 ) or the top model layer ( mkx ), then  no  !
          ! convective adjustment is performed and just exited.   However, in  !
          ! the revised code, I relaxed the first constraint and  even though  !
          ! LCL is at the lowest model layer, I allowed cumulus convection to  !
          ! be initiated. For this case, cumulus convection should be started  !
          ! from the PBL top height, as shown in the following code.           !
          ! When source air is already saturated even at the surface, klcl is  !
          ! set to 1.                                                          !
          ! ------------------------------------------------------------------ !

          plcl = qsinvert(qtsrc,thlsrc,ps0(0),qsat)
          DO k = 0, mkx
             IF( ps0(k) .LT. plcl ) THEN
                klcl = k
                go to 25
             ENDIF
          END DO
          klcl = mkx
25        CONTINUE
          klcl = MAX(1,klcl)

          IF( plcl .LT. 30000._r8 ) THEN               
             ! if( klcl .eq. mkx ) then          
             exit_klclmkx(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          ENDIF

          ! ------------------------------------------------------------- !
          ! Calculate environmental virtual potential temperature at LCL, !
          !'thv0lcl' which is solely used in the 'cin' calculation. Note  !
          ! that 'thv0lcl' is calculated first by calculating  'thl0lcl'  !
          ! and 'qt0lcl' at the LCL, and performing 'conden' afterward,   !
          ! in fully consistent with the other parts of the code.         !
          ! ------------------------------------------------------------- !

          thl0lcl = thl0(klcl) + ssthl0(klcl) * ( plcl - p0(klcl) )
          qt0lcl  = qt0(klcl)  + ssqt0(klcl)  * ( plcl - p0(klcl) )
          CALL conden(plcl,thl0lcl,qt0lcl,thj,qvj,qlj,qij,qse,id_check,qsat)
          IF( id_check .EQ. 1 ) THEN
             exit_conden(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          END IF
          thv0lcl = thj * ( 1._r8 + zvir * qvj - qlj - qij )

          ! ------------------------------------------------------------------------ !
          ! Compute Convective Inhibition, 'cin' & 'cinlcl' [J/kg]=[m2/s2] TKE unit. !
          !                                                                          !
          ! 'cin' (cinlcl) is computed from the PBL top interface to LFC (LCL) using ! 
          ! piecewisely reconstructed environmental profiles, assuming environmental !
          ! buoyancy profile within each layer ( or from LCL to upper interface in   !
          ! each layer ) is simply a linear profile. For the purpose of cin (cinlcl) !
          ! calculation, we simply assume that lateral entrainment does not occur in !
          ! updrafting cumulus plume, i.e., cumulus source air property is conserved.!
          ! Below explains some rules used in the calculations of cin (cinlcl).   In !
          ! general, both 'cin' and 'cinlcl' are calculated from a PBL top interface !
          ! to LCL and LFC, respectively :                                           !
          ! 1. If LCL is lower than the PBL height, cinlcl = 0 and cin is calculated !
          !    from PBL height to LFC.                                               !
          ! 2. If LCL is higher than PBL height,   'cinlcl' is calculated by summing !
          !    both positive and negative cloud buoyancy up to LCL using 'single_cin'!
          !    From the LCL to LFC, however, only negative cloud buoyancy is counted !
          !    to calculate final 'cin' upto LFC.                                    !
          ! 3. If either 'cin' or 'cinlcl' is negative, they are set to be zero.     !
          ! In the below code, 'klfc' is the layer index containing 'LFC' similar to !
          ! 'kinv' and 'klcl'.                                                       !
          ! ------------------------------------------------------------------------ !

          cin    = 0._r8
          cinlcl = 0._r8
          plfc   = 0._r8
          klfc   = mkx

          ! ------------------------------------------------------------------------- !
          ! Case 1. LCL height is higher than PBL interface ( 'pLCL <= ps0(kinv-1)' ) !
          ! ------------------------------------------------------------------------- !

          IF( klcl .GE. kinv ) THEN

             DO k = kinv, mkx - 1
                IF( k .LT. klcl ) THEN
                   thvubot = thvlsrc
                   thvutop = thvlsrc  
                   cin     = cin + single_cin(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop)
                ELSEIF( k .EQ. klcl ) THEN
                   !----- Bottom to LCL
                   thvubot = thvlsrc
                   thvutop = thvlsrc
                   cin     = cin + single_cin(ps0(k-1),thv0bot(k),plcl,thv0lcl,thvubot,thvutop)
                   IF( cin .LT. 0._r8 ) limit_cinlcl(i) = 1._r8
                   cinlcl  = MAX(cin,0._r8)
                   cin     = cinlcl
                   !----- LCL to Top
                   thvubot = thvlsrc
                   CALL conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   thvutop = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                   CALL getbuoy(plcl,thv0lcl,ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                   IF( plfc .GT. 0._r8 ) THEN 
                      klfc = k 
                      go to 35
                   END IF
                ELSE
                   thvubot = thvutop
                   CALL conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   thvutop = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                   CALL getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                   IF( plfc .GT. 0._r8 ) THEN 
                      klfc = k
                      go to 35
                   END IF
                ENDIF
             END DO

             ! ----------------------------------------------------------------------- !
             ! Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)' ) !
             ! ----------------------------------------------------------------------- !

          ELSE
             cinlcl = 0._r8 
             DO k = kinv, mkx - 1
                CALL conden(ps0(k-1),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                thvubot = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                CALL conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                thvutop = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                CALL getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                IF( plfc .GT. 0._r8 ) THEN 
                   klfc = k
                   go to 35
                END IF
             END DO
          ENDIF  ! End of CIN case selection

35        CONTINUE
          IF( cin .LT. 0._r8 ) limit_cin(i) = 1._r8
          cin = MAX(0._r8,cin)
          IF( klfc .GE. mkx ) THEN
             klfc = mkx
             ! write(iulog,*) 'klfc >= mkx'
             exit_klfcmkx(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          ENDIF

          ! ---------------------------------------------------------------------- !
          ! In order to calculate implicit 'cin' (or 'cinlcl'), save the initially !
          ! calculated 'cin' and 'cinlcl', and other related variables. These will !
          ! be restored after calculating implicit CIN.                            !
          ! ---------------------------------------------------------------------- !

          IF( iter .EQ. 1 ) THEN 
             cin_i       = cin
             cinlcl_i    = cinlcl
             ke          = rbuoy / ( rkfre * tkeavg + epsvarw ) 
             kinv_o      = kinv     
             klcl_o      = klcl     
             klfc_o      = klfc    
             plcl_o      = plcl    
             plfc_o      = plfc     
             tkeavg_o    = tkeavg   
             thvlmin_o   = thvlmin
             qtsrc_o     = qtsrc    
             thvlsrc_o   = thvlsrc  
             thlsrc_o    = thlsrc
             usrc_o      = usrc     
             vsrc_o      = vsrc     
             thv0lcl_o   = thv0lcl  
             DO m = 1, ncnst
                trsrc_o(m) = trsrc(m)
             ENDDO
          ENDIF

          ! Modification : If I impose w = max(0.1_r8, w) up to the top interface of
          !                klfc, I should only use cinlfc.  That is, if I want to
          !                use cinlcl, I should not impose w = max(0.1_r8, w).
          !                Using cinlcl is equivalent to treating only 'saturated'
          !                moist convection. Note that in this sense, I should keep
          !                the functionality of both cinlfc and cinlcl.
          !                However, the treatment of penetrative entrainment level becomes
          !                ambiguous if I choose 'cinlcl'. Thus, the best option is to use
          !                'cinlfc'.

          ! -------------------------------------------------------------------------- !
          ! Calculate implicit 'cin' by averaging initial and final cins.    Note that !
          ! implicit CIN is adopted only when cumulus convection stabilized the system,!
          ! i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note !
          ! also that since 'cinlcl' is set to zero whenever LCL is below the PBL top, !
          ! (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not !
          ! good. Thus, when using implicit CIN, always try to only use 'implicit CIN= !
          ! cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'  !
          ! are good when using explicit CIN.                                          !
          ! -------------------------------------------------------------------------- !

          IF( iter .NE. 1 ) THEN

             cin_f = cin
             cinlcl_f = cinlcl
             IF( use_CINcin ) THEN
                del_CIN = cin_f - cin_i
             ELSE
                del_CIN = cinlcl_f - cinlcl_i
             ENDIF

             IF( del_CIN .GT. 0._r8 ) THEN

                ! -------------------------------------------------------------- ! 
                ! Calculate implicit 'cin' and 'cinlcl'. Note that when we chose !
                ! to use 'implicit CIN = cin', choose 'cinlcl = cinlcl_i' below: !
                ! because iterative CIN only aims to obtain implicit CIN,  once  !
                ! we obtained 'implicit CIN=cin', it is good to use the original !
                ! profiles information for all the other variables after that.   !
                ! Note 'cinlcl' will be explicitly used in calculating  'wlcl' & !
                ! 'ufrclcl' after calculating 'winv' & 'ufrcinv'  at the PBL top !
                ! interface later, after calculating 'cbmf'.                     !
                ! -------------------------------------------------------------- !

                alpha = compute_alpha( del_CIN, ke ) 
                cin   = cin_i + alpha * del_CIN
                IF( use_CINcin ) THEN
                   cinlcl = cinlcl_i                 
                ELSE
                   cinlcl = cinlcl_i + alpha * del_cinlcl   
                ENDIF

                ! ----------------------------------------------------------------- !
                ! Restore the original values from the previous 'iter_cin' step (1) !
                ! to compute correct tendencies for (n+1) time step by implicit CIN !
                ! ----------------------------------------------------------------- !

                kinv      = kinv_o     
                klcl      = klcl_o     
                klfc      = klfc_o    
                plcl      = plcl_o    
                plfc      = plfc_o     
                tkeavg    = tkeavg_o   
                thvlmin   = thvlmin_o
                qtsrc     = qtsrc_o    
                thvlsrc   = thvlsrc_o  
                thlsrc    = thlsrc_o
                usrc      = usrc_o     
                vsrc      = vsrc_o     
                thv0lcl   = thv0lcl_o  
                DO m = 1, ncnst
                   trsrc(m) = trsrc_o(m)
                ENDDO

                qv0(:mkx)            = qv0_o(:mkx)
                ql0(:mkx)            = ql0_o(:mkx)
                qi0(:mkx)            = qi0_o(:mkx)
                t0(:mkx)             = t0_o(:mkx)
                s0(:mkx)             = s0_o(:mkx)
                u0(:mkx)             = u0_o(:mkx)
                v0(:mkx)             = v0_o(:mkx)
                qt0(:mkx)            = qt0_o(:mkx)
                thl0(:mkx)           = thl0_o(:mkx)
                thvl0(:mkx)          = thvl0_o(:mkx)
                ssthl0(:mkx)         = ssthl0_o(:mkx)
                ssqt0(:mkx)          = ssqt0_o(:mkx)
                thv0bot(:mkx)        = thv0bot_o(:mkx)
                thv0top(:mkx)        = thv0top_o(:mkx)
                thvl0bot(:mkx)       = thvl0bot_o(:mkx)
                thvl0top(:mkx)       = thvl0top_o(:mkx)
                ssu0(:mkx)           = ssu0_o(:mkx) 
                ssv0(:mkx)           = ssv0_o(:mkx) 
                DO m = 1, ncnst
                   tr0(:mkx,m)   = tr0_o(:mkx,m)
                   sstr0(:mkx,m) = sstr0_o(:mkx,m)
                ENDDO

                ! ------------------------------------------------------ !
                ! Initialize all fluxes, tendencies, and other variables ! 
                ! in association with cumulus convection.                !
                ! ------------------------------------------------------ ! 

                umf(0:mkx)          = 0.0_r8
                emf(0:mkx)          = 0.0_r8
                slflx(0:mkx)        = 0.0_r8
                qtflx(0:mkx)        = 0.0_r8
                uflx(0:mkx)         = 0.0_r8
                vflx(0:mkx)         = 0.0_r8
                qvten(:mkx)         = 0.0_r8
                qlten(:mkx)         = 0.0_r8
                qiten(:mkx)         = 0.0_r8
                sten(:mkx)          = 0.0_r8
                uten(:mkx)          = 0.0_r8
                vten(:mkx)          = 0.0_r8
                qrten(:mkx)         = 0.0_r8
                qsten(:mkx)         = 0.0_r8
                dwten(:mkx)         = 0.0_r8
                diten(:mkx)         = 0.0_r8
                precip              = 0.0_r8
                snow                = 0.0_r8
                evapc(:mkx)         = 0.0_r8
                cufrc(:mkx)         = 0.0_r8
                qcu(:mkx)           = 0.0_r8
                qlu(:mkx)           = 0.0_r8
                qiu(:mkx)           = 0.0_r8
                fer(:mkx)           = 0.0_r8
                fdr(:mkx)           = 0.0_r8
                qc(:mkx)            = 0.0_r8
                qc_l(:mkx)          = 0.0_r8
                qc_i(:mkx)          = 0.0_r8
                rliq                = 0.0_r8
                cbmf                = 0.0_r8
                cnt                 = REAL(mkx, r8)
                cnb                 = 0.0_r8
                qtten(:mkx)         = 0.0_r8
                slten(:mkx)         = 0.0_r8
                ufrc(0:mkx)         = 0.0_r8

                thlu(0:mkx)         = 0.0_r8
                qtu(0:mkx)          = 0.0_r8
                uu(0:mkx)           = 0.0_r8
                vu(0:mkx)           = 0.0_r8
                wu(0:mkx)           = 0.0_r8
                thvu(0:mkx)         = 0.0_r8
                thlu_emf(0:mkx)     = 0.0_r8
                qtu_emf(0:mkx)      = 0.0_r8
                uu_emf(0:mkx)       = 0.0_r8
                vu_emf(0:mkx)       = 0.0_r8

                DO m = 1, ncnst
                   trflx(0:mkx,m)   = 0.0_r8
                   trten(:mkx,m)    = 0.0_r8
                   tru(0:mkx,m)     = 0.0_r8
                   tru_emf(0:mkx,m) = 0.0_r8
                ENDDO

                ! -------------------------------------------------- !
                ! Below are diagnostic output variables for detailed !
                ! analysis of cumulus scheme.                        !
                ! -------------------------------------------------- ! 

                ufrcinvbase         = 0.0_r8
                ufrclcl             = 0.0_r8
                winvbase            = 0.0_r8
                wlcl                = 0.0_r8
                emfkbup             = 0.0_r8 
                cbmflimit           = 0.0_r8
                excessu_arr(:mkx)   = 0.0_r8
                excess0_arr(:mkx)   = 0.0_r8
                xc_arr(:mkx)        = 0.0_r8
                aquad_arr(:mkx)     = 0.0_r8
                bquad_arr(:mkx)     = 0.0_r8
                cquad_arr(:mkx)     = 0.0_r8
                bogbot_arr(:mkx)    = 0.0_r8
                bogtop_arr(:mkx)    = 0.0_r8

             ELSE ! When 'del_CIN < 0', use explicit CIN instead of implicit CIN.

                ! ----------------------------------------------------------- ! 
                ! Identifier showing whether explicit or implicit CIN is used !
                ! ----------------------------------------------------------- ! 

                ind_delcin(i) = 1._r8             

                ! --------------------------------------------------------- !
                ! Restore original output values of "iter_cin = 1" and exit !
                ! --------------------------------------------------------- !

                umf_out(i,0:mkx)         = umf_s(0:mkx)
                qvten_out(i,:mkx)        = qvten_s(:mkx)
                qlten_out(i,:mkx)        = qlten_s(:mkx)  
                qiten_out(i,:mkx)        = qiten_s(:mkx)
                sten_out(i,:mkx)         = sten_s(:mkx)
                uten_out(i,:mkx)         = uten_s(:mkx)  
                vten_out(i,:mkx)         = vten_s(:mkx)
                qrten_out(i,:mkx)        = qrten_s(:mkx)
                qsten_out(i,:mkx)        = qsten_s(:mkx)  
                precip_out(i)            = precip_s
                snow_out(i)              = snow_s
                evapc_out(i,:mkx)        = evapc_s(:mkx)
                cush_inout(i)            = cush_s
                cufrc_out(i,:mkx)        = cufrc_s(:mkx)  
                slflx_out(i,0:mkx)       = slflx_s(0:mkx)  
                qtflx_out(i,0:mkx)       = qtflx_s(0:mkx)
                qcu_out(i,:mkx)          = qcu_s(:mkx)    
                qlu_out(i,:mkx)          = qlu_s(:mkx)  
                qiu_out(i,:mkx)          = qiu_s(:mkx)  
                cbmf_out(i)              = cbmf_s
                qc_out(i,:mkx)           = qc_s(:mkx)  
                rliq_out(i)              = rliq_s
                cnt_out(i)               = cnt_s
                cnb_out(i)               = cnb_s
                DO m = 1, ncnst
                   trten_out(i,:mkx,m)   = trten_s(:mkx,m)
                ENDDO

                ! ------------------------------------------------------------------------------ ! 
                ! Below are diagnostic output variables for detailed analysis of cumulus scheme. !
                ! The order of vertical index is reversed for this internal diagnostic output.   !
                ! ------------------------------------------------------------------------------ !   

                fer_out(i,mkx:1:-1)      = fer_s(:mkx)  
                fdr_out(i,mkx:1:-1)      = fdr_s(:mkx)  
                cinh_out(i)              = cin_s
                cinlclh_out(i)           = cinlcl_s
                qtten_out(i,mkx:1:-1)    = qtten_s(:mkx)
                slten_out(i,mkx:1:-1)    = slten_s(:mkx)
                ufrc_out(i,mkx:0:-1)     = ufrc_s(0:mkx)
                uflx_out(i,mkx:0:-1)     = uflx_s(0:mkx)  
                vflx_out(i,mkx:0:-1)     = vflx_s(0:mkx)  

                ufrcinvbase_out(i)       = ufrcinvbase_s
                ufrclcl_out(i)           = ufrclcl_s 
                winvbase_out(i)          = winvbase_s
                wlcl_out(i)              = wlcl_s
                plcl_out(i)              = plcl_s
                pinv_out(i)              = pinv_s    
                plfc_out(i)              = plfc_s    
                pbup_out(i)              = pbup_s
                ppen_out(i)              = ppen_s    
                qtsrc_out(i)             = qtsrc_s
                thlsrc_out(i)            = thlsrc_s
                thvlsrc_out(i)           = thvlsrc_s
                emfkbup_out(i)           = emfkbup_s
                cbmflimit_out(i)         = cbmflimit_s
                tkeavg_out(i)            = tkeavg_s
                zinv_out(i)              = zinv_s
                rcwp_out(i)              = rcwp_s
                rlwp_out(i)              = rlwp_s
                riwp_out(i)              = riwp_s

                wu_out(i,mkx:0:-1)       = wu_s(0:mkx)
                qtu_out(i,mkx:0:-1)      = qtu_s(0:mkx)
                thlu_out(i,mkx:0:-1)     = thlu_s(0:mkx)
                thvu_out(i,mkx:0:-1)     = thvu_s(0:mkx)
                uu_out(i,mkx:0:-1)       = uu_s(0:mkx)
                vu_out(i,mkx:0:-1)       = vu_s(0:mkx)
                qtu_emf_out(i,mkx:0:-1)  = qtu_emf_s(0:mkx)
                thlu_emf_out(i,mkx:0:-1) = thlu_emf_s(0:mkx)
                uu_emf_out(i,mkx:0:-1)   = uu_emf_s(0:mkx)
                vu_emf_out(i,mkx:0:-1)   = vu_emf_s(0:mkx)
                uemf_out(i,mkx:0:-1)     = uemf_s(0:mkx)

                dwten_out(i,mkx:1:-1)    = dwten_s(:mkx)
                diten_out(i,mkx:1:-1)    = diten_s(:mkx)
                flxrain_out(i,mkx:0:-1)  = flxrain_s(0:mkx)
                flxsnow_out(i,mkx:0:-1)  = flxsnow_s(0:mkx)
                ntraprd_out(i,mkx:1:-1)  = ntraprd_s(:mkx)
                ntsnprd_out(i,mkx:1:-1)  = ntsnprd_s(:mkx)

                excessu_arr_out(i,mkx:1:-1)  = excessu_arr_s(:mkx)
                excess0_arr_out(i,mkx:1:-1)  = excess0_arr_s(:mkx)
                xc_arr_out(i,mkx:1:-1)       = xc_arr_s(:mkx)
                aquad_arr_out(i,mkx:1:-1)    = aquad_arr_s(:mkx)
                bquad_arr_out(i,mkx:1:-1)    = bquad_arr_s(:mkx)
                cquad_arr_out(i,mkx:1:-1)    = cquad_arr_s(:mkx)
                bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr_s(:mkx)
                bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr_s(:mkx)

                DO m = 1, ncnst
                   trflx_out(i,mkx:0:-1,m)   = trflx_s(0:mkx,m)  
                   tru_out(i,mkx:0:-1,m)     = tru_s(0:mkx,m)
                   tru_emf_out(i,mkx:0:-1,m) = tru_emf_s(0:mkx,m)
                ENDDO

                id_exit = .FALSE.
                go to 333

             ENDIF

          ENDIF

          ! ------------------------------------------------------------------ !
          ! Define a release level, 'prel' and release layer, 'krel'.          !
          ! 'prel' is the lowest level from which buoyancy sorting occurs, and !
          ! 'krel' is the layer index containing 'prel' in it, similar to  the !
          ! previous definitions of 'kinv', 'klcl', and 'klfc'.    In order to !
          ! ensure that only PBL scheme works within the PBL,  if LCL is below !
          ! PBL top height, then 'krel = kinv', while if LCL is above  PBL top !
          ! height, then 'krel = klcl'.   Note however that regardless of  the !
          ! definition of 'krel', cumulus convection induces fluxes within PBL !
          ! through 'fluxbelowinv'.  We can make cumulus convection start from !
          ! any level, even within the PBL by appropriately defining 'krel'  & !
          ! 'prel' here. Then it must be accompanied by appropriate definition !
          ! of source air properties, CIN, and re-setting of 'fluxbelowinv', & !
          ! many other stuffs.                                                 !
          ! Note that even when 'prel' is located above the PBL top height, we !
          ! still have cumulus convection between PBL top height and 'prel':   !
          ! we simply assume that no lateral mixing occurs in this range.      !
          ! ------------------------------------------------------------------ !

          IF( klcl .LT. kinv ) THEN
             krel    = kinv
             prel    = ps0(krel-1)
             thv0rel = thv0bot(krel) 
          ELSE
             krel    = klcl
             prel    = plcl 
             thv0rel = thv0lcl
          ENDIF

          ! --------------------------------------------------------------------------- !
          ! Calculate cumulus base mass flux ('cbmf'), fractional area ('ufrcinv'), and !
          ! and mean vertical velocity (winv) of cumulus updraft at PBL top interface.  !
          ! Also, calculate updraft fractional area (ufrclcl) and vertical velocity  at !
          ! the LCL (wlcl). When LCL is below PBLH, cinlcl = 0 and 'ufrclcl = ufrcinv', !
          ! and 'wlcl = winv.                                                           !
          ! Only updrafts strong enough to overcome CIN can rise over PBL top interface.! 
          ! Thus,  in order to calculate cumulus mass flux at PBL top interface, 'cbmf',!
          ! we need to know 'CIN' ( the strength of potential energy barrier ) and      !
          ! 'sigmaw' ( a standard deviation of updraft vertical velocity at the PBL top !
          ! interface, a measure of turbulentce strength in the PBL ).   Naturally, the !
          ! ratio of these two variables, 'mu' - normalized CIN by TKE- is key variable !
          ! controlling 'cbmf'.  If 'mu' becomes large, only small fraction of updrafts !
          ! with very strong TKE can rise over the PBL - both 'cbmf' and 'ufrc' becomes !
          ! small, but 'winv' becomes large ( this can be easily understood by PDF of w !
          ! at PBL top ).  If 'mu' becomes small, lots of updraft can rise over the PBL !
          ! top - both 'cbmf' and 'ufrc' becomes large, but 'winv' becomes small. Thus, !
          ! all of the key variables associated with cumulus convection  at the PBL top !
          ! - 'cbmf', 'ufrc', 'winv' where 'cbmf = rho*ufrc*winv' - are a unique functi !
          ! ons of 'mu', normalized CIN. Although these are uniquely determined by 'mu',! 
          ! we usually impose two comstraints on 'cbmf' and 'ufrc': (1) because we will !
          ! simply assume that subsidence warming and drying of 'kinv-1' layer in assoc !
          ! iation with 'cbmf' at PBL top interface is confined only in 'kinv-1' layer, !
          ! cbmf must not be larger than the mass within the 'kinv-1' layer. Otherwise, !
          ! instability will occur due to the breaking of stability con. If we consider !
          ! semi-Lagrangian vertical advection scheme and explicitly consider the exten !
          ! t of vertical movement of each layer in association with cumulus mass flux, !
          ! we don't need to impose this constraint. However,  using a  semi-Lagrangian !
          ! scheme is a future research subject. Note that this constraint should be ap !
          ! plied for all interfaces above PBL top as well as PBL top interface.   As a !
          ! result, this 'cbmf' constraint impose a 'lower' limit on mu - 'mumin0'. (2) !
          ! in order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to   be !
          ! valid, cumulus updraft fractional area should be much smaller than 1.    In !
          ! current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole vertical !
          ! layers where cumulus convection occurs. At the PBL top interface,  the same !
          ! constraint is made by imposing another lower 'lower' limit on mu, 'mumin1'. !
          ! After that, also limit 'ufrclcl' to be smaller than 'rmaxfrac' by 'mumin2'. !
          ! --------------------------------------------------------------------------- !

          ! --------------------------------------------------------------------------- !
          ! Calculate normalized CIN, 'mu' satisfying all the three constraints imposed !
          ! on 'cbmf'('mumin0'), 'ufrc' at the PBL top - 'ufrcinv' - ( by 'mumin1' from !
          ! a parameter sentence), and 'ufrc' at the LCL - 'ufrclcl' ( by 'mumin2').    !
          ! Note that 'cbmf' does not change between PBL top and LCL  because we assume !
          ! that buoyancy sorting does not occur when cumulus updraft is unsaturated.   !
          ! --------------------------------------------------------------------------- !

          IF( use_CINcin ) THEN       
             wcrit = SQRT( 2._r8 * cin * rbuoy )      
          ELSE
             wcrit = SQRT( 2._r8 * cinlcl * rbuoy )   
          ENDIF
          sigmaw = SQRT( rkfre * tkeavg + epsvarw )
          mu = wcrit/sigmaw/1.4142_r8                  
          IF( mu .GE. 3._r8 ) THEN
             ! write(iulog,*) 'mu >= 3'
             id_exit = .TRUE.
             go to 333
          ENDIF
          rho0inv = ps0(kinv-1)/(r*thv0top(kinv-1)*exns0(kinv-1))
          cbmf = (rho0inv*sigmaw/2.5066_r8)*EXP(-mu**2)
          ! 1. 'cbmf' constraint
          cbmflimit = 0.9_r8*dp0(kinv-1)/g/dt
          mumin0 = 0._r8
          IF( cbmf .GT. cbmflimit ) mumin0 = SQRT(-LOG(2.5066_r8*cbmflimit/rho0inv/sigmaw))
          ! 2. 'ufrcinv' constraint
          mu = MAX(MAX(mu,mumin0),mumin1)
          ! 3. 'ufrclcl' constraint      
          mulcl = SQRT(2._r8*cinlcl*rbuoy)/1.4142_r8/sigmaw
          mulclstar = SQRT(MAX(0._r8,2._r8*(EXP(-mu**2)/2.5066_r8)**2*(1._r8/erfc(mu)**2-0.25_r8/rmaxfrac**2)))
          IF( mulcl .GT. 1.e-8_r8 .AND. mulcl .GT. mulclstar ) THEN
             mumin2 = compute_mumin2(mulcl,rmaxfrac,mu)
             IF( mu .GT. mumin2 ) THEN
                WRITE(iulog,*) 'Critical error in mu calculation in UW_ShCu'
                CALL endrun
             ENDIF
             mu = MAX(mu,mumin2)
             IF( mu .EQ. mumin2 ) limit_ufrc(i) = 1._r8
          ENDIF
          IF( mu .EQ. mumin0 ) limit_cbmf(i) = 1._r8
          IF( mu .EQ. mumin1 ) limit_ufrc(i) = 1._r8

          ! ------------------------------------------------------------------- !    
          ! Calculate final ['cbmf','ufrcinv','winv'] at the PBL top interface. !
          ! Note that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
          ! 'ufrclcl' are smaller than ufrcmax with no instability.             !
          ! ------------------------------------------------------------------- !

          cbmf = (rho0inv*sigmaw/2.5066_r8)*EXP(-mu**2)                       
          winv = sigmaw*(2._r8/2.5066_r8)*EXP(-mu**2)/erfc(mu)
          ufrcinv = cbmf/winv/rho0inv

          ! ------------------------------------------------------------------- !
          ! Calculate ['ufrclcl','wlcl'] at the LCL. When LCL is below PBL top, !
          ! it automatically becomes 'ufrclcl = ufrcinv' & 'wlcl = winv', since !
          ! it was already set to 'cinlcl=0' if LCL is below PBL top interface. !
          ! Note 'cbmf' at the PBL top is the same as 'cbmf' at the LCL.  Note  !
          ! also that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
          ! 'ufrclcl' are smaller than ufrcmax and there is no instability.     !
          ! By construction, it must be 'wlcl > 0' but for assurance, I checked !
          ! this again in the below block. If 'ufrclcl < 0.1%', just exit.      !
          ! ------------------------------------------------------------------- !

          wtw = winv * winv - 2._r8 * cinlcl * rbuoy
          IF( wtw .LE. 0._r8 ) THEN
             ! write(iulog,*) 'wlcl < 0 at the LCL'
             exit_wtw(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          ENDIF
          wlcl = SQRT(wtw)
          ufrclcl = cbmf/wlcl/rho0inv
          wrel = wlcl
          IF( ufrclcl .LE. 0.0001_r8 ) THEN
             ! write(iulog,*) 'ufrclcl <= 0.0001' 
             exit_ufrc(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          ENDIF
          ufrc(krel-1) = ufrclcl

          ! ----------------------------------------------------------------------- !
          ! Below is just diagnostic output for detailed analysis of cumulus scheme !
          ! ----------------------------------------------------------------------- !

          ufrcinvbase        = ufrcinv
          winvbase           = winv
          umf(kinv-1:krel-1) = cbmf   
          wu(kinv-1:krel-1)  = winv   

          ! -------------------------------------------------------------------------- ! 
          ! Define updraft properties at the level where buoyancy sorting starts to be !
          ! happening, i.e., by definition, at 'prel' level within the release layer.  !
          ! Because no lateral entrainment occurs upto 'prel', conservative scalars of ! 
          ! cumulus updraft at release level is same as those of source air.  However, ! 
          ! horizontal momentums of source air are modified by horizontal PGF forcings ! 
          ! from PBL top interface to 'prel'.  For this case, we should add additional !
          ! horizontal momentum from PBL top interface to 'prel' as will be done below !
          ! to 'usrc' and 'vsrc'. Note that below cumulus updraft properties - umf, wu,!
          ! thlu, qtu, thvu, uu, vu - are defined all interfaces not at the layer mid- !
          ! point. From the index notation of cumulus scheme, wu(k) is the cumulus up- !
          ! draft vertical velocity at the top interface of k layer.                   !
          ! Diabatic horizontal momentum forcing should be treated as a kind of 'body' !
          ! forcing without actual mass exchange between convective updraft and        !
          ! environment, but still taking horizontal momentum from the environment to  !
          ! the convective updrafts. Thus, diabatic convective momentum transport      !
          ! vertically redistributes environmental horizontal momentum.                !
          ! -------------------------------------------------------------------------- !

          emf(krel-1)  = 0._r8
          umf(krel-1)  = cbmf
          wu(krel-1)   = wrel
          thlu(krel-1) = thlsrc
          qtu(krel-1)  = qtsrc
          CALL conden(prel,thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
          IF( id_check .EQ. 1 ) THEN
             exit_conden(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          ENDIF
          thvu(krel-1) = thj * ( 1._r8 + zvir*qvj - qlj - qij )       

          uplus = 0._r8
          vplus = 0._r8
          IF( krel .EQ. kinv ) THEN
             uplus = PGFc * ssu0(kinv) * ( prel - ps0(kinv-1) )
             vplus = PGFc * ssv0(kinv) * ( prel - ps0(kinv-1) )
          ELSE
             DO k = kinv, MAX(krel-1,kinv)
                uplus = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
                vplus = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
             END DO
             uplus = uplus + PGFc * ssu0(krel) * ( prel - ps0(krel-1) )
             vplus = vplus + PGFc * ssv0(krel) * ( prel - ps0(krel-1) )
          END IF
          uu(krel-1) = usrc + uplus
          vu(krel-1) = vsrc + vplus      

          DO m = 1, ncnst
             tru(krel-1,m)  = trsrc(m)
          ENDDO

          ! -------------------------------------------------------------------------- !
          ! Define environmental properties at the level where buoyancy sorting occurs !
          ! ('pe', normally, layer midpoint except in the 'krel' layer). In the 'krel' !
          ! layer where buoyancy sorting starts to occur, however, 'pe' is defined     !
          ! differently because LCL is regarded as lower interface for mixing purpose. !
          ! -------------------------------------------------------------------------- !

          pe      = 0.5_r8 * ( prel + ps0(krel) )
          dpe     = prel - ps0(krel)
          exne    = exnf(pe)
          thvebot = thv0rel
          thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
          qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
          ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
          ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )
          DO m = 1, ncnst
             tre(m) = tr0(krel,m)  + sstr0(krel,m) * ( pe - p0(krel) )
          ENDDO

          !-------------------------! 
          ! Buoyancy-Sorting Mixing !
          !-------------------------!------------------------------------------------ !
          !                                                                           !
          !  In order to complete buoyancy-sorting mixing at layer mid-point, and so  ! 
          !  calculate 'updraft mass flux, updraft w velocity, conservative scalars'  !
          !  at the upper interface of each layer, we need following 3 information.   ! 
          !                                                                           !
          !  1. Pressure where mixing occurs ('pe'), and temperature at 'pe' which is !
          !     necessary to calculate various thermodynamic coefficients at pe. This !
          !     temperature is obtained by undiluted cumulus properties lifted to pe. ! 
          !  2. Undiluted updraft properties at pe - conservative scalar and vertical !
          !     velocity -which are assumed to be the same as the properties at lower !
          !     interface only for calculation of fractional lateral entrainment  and !
          !     detrainment rate ( fer(k) and fdr(k) [Pa-1] ), respectively.    Final !
          !     values of cumulus conservative scalars and w at the top interface are !
          !     calculated afterward after obtaining fer(k) & fdr(k).                 !
          !  3. Environmental properties at pe.                                       !
          ! ------------------------------------------------------------------------- !

          ! ------------------------------------------------------------------------ ! 
          ! Define cumulus scale height.                                             !
          ! Cumulus scale height is defined as the maximum height cumulus can reach. !
          ! In case of premitive code, cumulus scale height ('cush')  at the current !
          ! time step was assumed to be the same as 'cush' of previous time step.    !
          ! However, I directly calculated cush at each time step using an iterative !
          ! method. Note that within the cumulus scheme, 'cush' information is  used !
          ! only at two places during buoyancy-sorting process:                      !
          ! (1) Even negatively buoyancy mixtures with strong vertical velocity      !
          !     enough to rise up to 'rle*scaleh' (rle = 0.1) from pe are entrained  !
          !     into cumulus updraft,                                                !  
          ! (2) The amount of mass that is involved in buoyancy-sorting mixing       !
          !      process at pe is rei(k) = rkm/scaleh/rho*g [Pa-1]                   !
          ! In terms of (1), I think critical stopping distance might be replaced by !
          ! layer thickness. In future, we will use rei(k) = (0.5*rkm/z0(k)/rho/g).  !
          ! In the premitive code,  'scaleh' was largely responsible for the jumping !
          ! variation of precipitation amount.                                       !
          ! ------------------------------------------------------------------------ !   

          scaleh = tscaleh
          IF( tscaleh .LT. 0.0_r8 ) scaleh = 1000._r8 

          ! Save time : Set iter_scaleh = 1. This will automatically use 'cush' from the previous time step
          !             at the first implicit iteration. At the second implicit iteration, it will use
          !             the updated 'cush' by the first implicit cin. So, this updating has an effect of
          !             doing one iteration for cush calculation, which is good. 
          !             So, only this setting of 'iter_scaleh = 1' is sufficient-enough to save computation time.
          ! OK

          DO iter_scaleh = 1, 3

             ! ---------------------------------------------------------------- !
             ! Initialization of 'kbup' and 'kpen'                              !
             ! ---------------------------------------------------------------- !
             ! 'kbup' is the top-most layer in which cloud buoyancy is positive !
             ! both at the top and bottom interface of the layer. 'kpen' is the !
             ! layer upto which cumulus panetrates ,i.e., cumulus w at the base !
             ! interface is positive, but becomes negative at the top interface.!
             ! Here, we initialize 'kbup' and 'kpen'. These initializations are !  
             ! not trivial but important, expecially   in calculating turbulent !
             ! fluxes without confliction among several physics as explained in !
             ! detail in the part of turbulent fluxes calculation later.   Note !
             ! that regardless of whether 'kbup' and 'kpen' are updated or  not !
             ! during updraft motion,  penetrative entrainments are dumped down !
             ! across the top interface of 'kbup' later.      More specifically,!
             ! penetrative entrainment heat and moisture fluxes are  calculated !
             ! from the top interface of 'kbup' layer  to the base interface of !
             ! 'kpen' layer. Because of this, initialization of 'kbup' & 'kpen' !
             ! influence the convection system when there are not updated.  The !  
             ! below initialization of 'kbup = krel' assures  that  penetrative !
             ! entrainment fluxes always occur at interfaces above the PBL  top !
             ! interfaces (i.e., only at interfaces k >=kinv ), which seems  to !
             ! be attractable considering that the most correct fluxes  at  the !
             ! PBL top interface can be ontained from the 'fluxbelowinv'  using !
             ! reconstructed PBL height.                                        ! 
             ! The 'kbup = krel'(after going through the whole buoyancy sorting !
             ! proces during updraft motion) implies that cumulus updraft  from !
             ! the PBL top interface can not reach to the LFC,so that 'kbup' is !
             ! not updated during upward. This means that cumulus updraft   did !
             ! not fully overcome the buoyancy barrier above just the PBL top.  !
             ! If 'kpen' is not updated either ( i.e., cumulus cannot rise over !
             ! the top interface of release layer),penetrative entrainment will !
             ! not happen at any interfaces.  If cumulus updraft can rise above !
             ! the release layer but cannot fully overcome the buoyancy barrier !
             ! just above PBL top interface, penetratve entrainment   occurs at !
             ! several above interfaces, including the top interface of release ! 
             ! layer. In the latter case, warming and drying tendencies will be !
             ! be initiated in 'krel' layer. Note current choice of 'kbup=krel' !
             ! is completely compatible with other flux physics without  double !
             ! or miss counting turbulent fluxes at any interface. However, the !
             ! alternative choice of 'kbup=krel-1' also has itw own advantage - !
             ! when cumulus updraft cannot overcome buoyancy barrier just above !
             ! PBL top, entrainment warming and drying are concentrated in  the !
             ! 'kinv-1' layer instead of 'kinv' layer for this case. This might !
             ! seems to be more dynamically reasonable, but I will choose the   !
             ! 'kbup = krel' choice since it is more compatible  with the other !
             ! parts of the code, expecially, when we chose ' use_emf=.false. ' !
             ! as explained in detail in turbulent flux calculation part.       !
             ! ---------------------------------------------------------------- ! 

             kbup    = krel
             kpen    = krel

             ! ------------------------------------------------------------ !
             ! Since 'wtw' is continuously updated during vertical motion,  !
             ! I need below initialization command within this 'iter_scaleh'!
             ! do loop. Similarily, I need initializations of environmental !
             ! properties at 'krel' layer as below.                         !
             ! ------------------------------------------------------------ !

             wtw     = wlcl * wlcl
             pe      = 0.5_r8 * ( prel + ps0(krel) )
             dpe     = prel - ps0(krel)
             exne    = exnf(pe)
             thvebot = thv0rel
             thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
             qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
             ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
             ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )
             DO m = 1, ncnst
                tre(m) = tr0(krel,m)  + sstr0(krel,m)  * ( pe - p0(krel) )
             ENDDO

             ! ----------------------------------------------------------------------- !
             ! Cumulus rises upward from 'prel' ( or base interface of  'krel' layer ) !
             ! until updraft vertical velocity becomes zero.                           !
             ! Buoyancy sorting is performed via two stages. (1) Using cumulus updraft !
             ! properties at the base interface of each layer,perform buoyancy sorting !
             ! at the layer mid-point, 'pe',  and update cumulus properties at the top !
             ! interface, and then  (2) by averaging updated cumulus properties at the !
             ! top interface and cumulus properties at the base interface,   calculate !
             ! cumulus updraft properties at pe that will be used  in buoyancy sorting !
             ! mixing - thlue, qtue and, wue.  Using this averaged properties, perform !
             ! buoyancy sorting again at pe, and re-calculate fer(k) and fdr(k). Using !
             ! this recalculated fer(k) and fdr(k),  finally calculate cumulus updraft !
             ! properties at the top interface - thlu, qtu, thvu, uu, vu. In the below,!
             ! 'iter_xc = 1' performs the first stage, while 'iter_xc= 2' performs the !
             ! second stage. We can increase the number of iterations, 'nter_xc'.as we !
             ! want, but a sample test indicated that about 3 - 5 iterations  produced !
             ! satisfactory converent solution. Finally, identify 'kbup' and 'kpen'.   !
             ! ----------------------------------------------------------------------- !

             DO k = krel, mkx - 1 ! Here, 'k' is a layer index.

                km1 = k - 1

                thlue = thlu(km1)
                qtue  = qtu(km1)    
                wue   = wu(km1)
                wtwb  = wtw  

                DO iter_xc = 1, niter_xc

                   wtw = wu(km1) * wu(km1)

                   ! ---------------------------------------------------------------- !
                   ! Calculate environmental and cumulus saturation 'excess' at 'pe'. !
                   ! Note that in order to calculate saturation excess, we should use ! 
                   ! liquid water temperature instead of temperature  as the argument !
                   ! of "qsat". But note normal argument of "qsat" is temperature.    !
                   ! ---------------------------------------------------------------- !

                   CALL conden(pe,thle,qte,thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   thv0j    = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                   rho0j    = pe / ( r * thv0j * exne )
                   qsat_arg = thle*exne     
                   status   = qsat(qsat_arg,pe,es(1),qs(1),gam(1),1)
                   excess0  = qte - qs(1)

                   CALL conden(pe,thlue,qtue,thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   ! ----------------------------------------------------------------- !
                   ! Detrain excessive condensate larger than 'criqc' from the cumulus ! 
                   ! updraft before performing buoyancy sorting. All I should to do is !
                   ! to update 'thlue' &  'que' here. Below modification is completely !
                   ! compatible with the other part of the code since 'thule' & 'qtue' !
                   ! are used only for buoyancy sorting. I found that as long as I use !
                   ! 'niter_xc >= 2',  detraining excessive condensate before buoyancy !
                   ! sorting has negligible influence on the buoyancy sorting results. !   
                   ! ----------------------------------------------------------------- !
                   IF( (qlj + qij) .GT. criqc ) THEN
                      exql  = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
                      exqi  = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
                      qtue  = qtue - exql - exqi
                      thlue = thlue + (xlv/cp/exne)*exql + (xls/cp/exne)*exqi 
                   ENDIF
                   CALL conden(pe,thlue,qtue,thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   thvj     = thj * ( 1._r8 + zvir * qvj - qlj - qij )
                   tj       = thj * exne ! This 'tj' is used for computing thermo. coeffs. below
                   qsat_arg = thlue*exne
                   status   = qsat(qsat_arg,pe,es(1),qs(1),gam(1),1)
                   excessu  = qtue - qs(1)

                   ! ------------------------------------------------------------------- !
                   ! Calculate critical mixing fraction, 'xc'. Mixture with mixing ratio !
                   ! smaller than 'xc' will be entrained into cumulus updraft.  Both the !
                   ! saturated updrafts with 'positive buoyancy' or 'negative buoyancy + ! 
                   ! strong vertical velocity enough to rise certain threshold distance' !
                   ! are kept into the updraft in the below program. If the core updraft !
                   ! is unsaturated, we can set 'xc = 0' and let the cumulus  convection !
                   ! still works or we may exit.                                         !
                   ! Current below code does not entrain unsaturated mixture. However it !
                   ! should be modified such that it also entrain unsaturated mixture.   !
                   ! ------------------------------------------------------------------- !

                   ! ----------------------------------------------------------------- !
                   ! cridis : Critical stopping distance for buoyancy sorting purpose. !
                   !          scaleh is only used here.                                !
                   ! ----------------------------------------------------------------- !

                   cridis = rle*scaleh                 ! Original code
                   ! cridis = 1._r8*(zs0(k) - zs0(k-1))  ! New code

                   ! ---------------- !
                   ! Buoyancy Sorting !
                   ! ---------------- !                   

                   ! ----------------------------------------------------------------- !
                   ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
                   ! ----------------------------------------------------------------- !

                   IF( ( excessu .LE. 0._r8 .AND. excess0 .LE. 0._r8 ) .OR. ( excessu .GE. 0._r8 .AND. excess0 .GE. 0._r8 ) ) THEN
                      xc = MIN(1._r8,MAX(0._r8,1._r8-2._r8*rbuoy*g*cridis/wue**2._r8*(1._r8-thvj/thv0j)))
                      ! Below 3 lines are diagnostic output not influencing
                      ! numerical calculations.
                      aquad = 0._r8
                      bquad = 0._r8
                      cquad = 0._r8
                   ELSE
                      ! -------------------------------------------------- !
                      ! Case 2 : When either cumulus or env. is saturated. !
                      ! -------------------------------------------------- !
                      xsat    = excessu / ( excessu - excess0 );
                      thlxsat = thlue + xsat * ( thle - thlue );
                      qtxsat  = qtue  + xsat * ( qte - qtue );
                      CALL conden(pe,thlxsat,qtxsat,thj,qvj,qlj,qij,qse,id_check,qsat)
                      IF( id_check .EQ. 1 ) THEN
                         exit_conden(i) = 1._r8
                         id_exit = .TRUE.
                         go to 333
                      END IF
                      thvxsat = thj * ( 1._r8 + zvir * qvj - qlj - qij )               
                      ! -------------------------------------------------- !
                      ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
                      ! -------------------------------------------------- ! 
                      DO kk = 1, 2 
                         IF( kk .EQ. 1 ) THEN
                            thv_x0 = thvj;
                            thv_x1 = ( 1._r8 - 1._r8/xsat ) * thvj + ( 1._r8/xsat ) * thvxsat;
                         ELSE
                            thv_x1 = thv0j;
                            thv_x0 = ( xsat / ( xsat - 1._r8 ) ) * thv0j + ( 1._r8/( 1._r8 - xsat ) ) * thvxsat;
                         ENDIF
                         aquad =  wue**2;
                         bquad =  2._r8*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv0j - 2._r8*wue**2;
                         cquad =  2._r8*rbuoy*g*cridis*(thv_x0 -  thv0j)/thv0j +       wue**2;
                         IF( kk .EQ. 1 ) THEN
                            IF( ( bquad**2-4._r8*aquad*cquad ) .GE. 0._r8 ) THEN
                               CALL roots(aquad,bquad,cquad,xs1,xs2,status)
                               x_cu = MIN(1._r8,MAX(0._r8,MIN(xsat,MIN(xs1,xs2))))
                            ELSE
                               x_cu = xsat;
                            ENDIF
                         ELSE 
                            IF( ( bquad**2-4._r8*aquad*cquad) .GE. 0._r8 ) THEN
                               CALL roots(aquad,bquad,cquad,xs1,xs2,status)
                               x_en = MIN(1._r8,MAX(0._r8,MAX(xsat,MIN(xs1,xs2))))
                            ELSE
                               x_en = 1._r8;
                            ENDIF
                         ENDIF
                      ENDDO
                      IF( x_cu .EQ. xsat ) THEN
                         xc = MAX(x_cu, x_en);
                      ELSE
                         xc = x_cu;
                      ENDIF
                   ENDIF

                   ! ------------------------------------------------------------------------ !
                   ! Compute fractional lateral entrainment & detrainment rate in each layers.!
                   ! The unit of rei(k), fer(k), and fdr(k) is [Pa-1].  Alternative choice of !
                   ! 'rei(k)' is also shown below, where coefficient 0.5 was from approximate !
                   ! tuning against the BOMEX case.                                           !
                   ! In order to prevent the onset of instability in association with cumulus !
                   ! induced subsidence advection, cumulus mass flux at the top interface  in !
                   ! any layer should be smaller than ( 90% of ) total mass within that layer.!
                   ! I imposed limits on 'rei(k)' as below,  in such that stability condition ! 
                   ! is always satisfied.                                                     !
                   ! Below limiter of 'rei(k)' becomes negative for some cases, causing error.!
                   ! So, for the time being, I came back to the original limiter.             !
                   ! ------------------------------------------------------------------------ !
                   ee2    = xc**2
                   ud2    = 1._r8 - 2._r8*xc + xc**2
                   ! rei(k) = ( rkm / scaleh / g / rho0j )        ! Default.
                   rei(k) = ( 0.5_r8 * rkm / z0(k) / g /rho0j ) ! Alternative.
                   IF( xc .GT. 0.5_r8 ) rei(k) = MIN(rei(k),0.9_r8*LOG(dp0(k)/g/dt/umf(km1) + 1._r8)/dpe/(2._r8*xc-1._r8))
                   fer(k) = rei(k) * ee2
                   fdr(k) = rei(k) * ud2

                   ! ------------------------------------------------------------------------------ !
                   ! Iteration Start due to 'maxufrc' constraint [ ****************************** ] ! 
                   ! ------------------------------------------------------------------------------ !

                   ! -------------------------------------------------------------------------- !
                   ! Calculate cumulus updraft mass flux and penetrative entrainment mass flux. !
                   ! Note that  non-zero penetrative entrainment mass flux will be asigned only !
                   ! to interfaces from the top interface of 'kbup' layer to the base interface !
                   ! of 'kpen' layer as will be shown later.                                    !
                   ! -------------------------------------------------------------------------- !

                   umf(k) = umf(km1) * EXP( dpe * ( fer(k) - fdr(k) ) )
                   emf(k) = 0._r8    

                   ! --------------------------------------------------------- !
                   ! Compute cumulus updraft properties at the top interface.  !
                   ! Also use Tayler expansion in order to treat limiting case !
                   ! --------------------------------------------------------- !

                   IF( fer(k)*dpe .LT. 1.e-4_r8 ) THEN
                      thlu(k) = thlu(km1) + ( thle + ssthl0(k) * dpe / 2._r8 - thlu(km1) ) * fer(k) * dpe
                      qtu(k)  =  qtu(km1) + ( qte  +  ssqt0(k) * dpe / 2._r8 -  qtu(km1) ) * fer(k) * dpe
                      uu(k)   =   uu(km1) + ( ue   +   ssu0(k) * dpe / 2._r8 -   uu(km1) ) * fer(k) * dpe - PGFc * ssu0(k) * dpe
                      vu(k)   =   vu(km1) + ( ve   +   ssv0(k) * dpe / 2._r8 -   vu(km1) ) * fer(k) * dpe - PGFc * ssv0(k) * dpe
                      DO m = 1, ncnst
                         tru(k,m)  =  tru(km1,m) + ( tre(m)  + sstr0(k,m) * dpe / 2._r8  -  tru(km1,m) ) * fer(k) * dpe
                      ENDDO
                   ELSE
                      thlu(k) = ( thle + ssthl0(k) / fer(k) - ssthl0(k) * dpe / 2._r8 ) -          &
                           ( thle + ssthl0(k) * dpe / 2._r8 - thlu(km1) + ssthl0(k) / fer(k) ) * EXP(-fer(k) * dpe)
                      qtu(k)  = ( qte  +  ssqt0(k) / fer(k) -  ssqt0(k) * dpe / 2._r8 ) -          &  
                           ( qte  +  ssqt0(k) * dpe / 2._r8 -  qtu(km1) +  ssqt0(k) / fer(k) ) * EXP(-fer(k) * dpe)
                      uu(k) =   ( ue + ( 1._r8 - PGFc ) * ssu0(k) / fer(k) - ssu0(k) * dpe / 2._r8 ) - &
                           ( ue +     ssu0(k) * dpe / 2._r8 -   uu(km1) + ( 1._r8 - PGFc ) * ssu0(k) / fer(k) ) * EXP(-fer(k) * dpe)
                      vu(k) =   ( ve + ( 1._r8 - PGFc ) * ssv0(k) / fer(k) - ssv0(k) * dpe / 2._r8 ) - &
                           ( ve +     ssv0(k) * dpe / 2._r8 -   vu(km1) + ( 1._r8 - PGFc ) * ssv0(k) / fer(k) ) * EXP(-fer(k) * dpe)
                      DO m = 1, ncnst
                         tru(k,m)  = ( tre(m)  + sstr0(k,m) / fer(k) - sstr0(k,m) * dpe / 2._r8 ) - &  
                              ( tre(m)  + sstr0(k,m) * dpe / 2._r8 - tru(km1,m) + sstr0(k,m) / fer(k) ) * EXP(-fer(k) * dpe)
                      ENDDO
                   END IF

                   !------------------------------------------------------------------- !
                   ! Expel some of cloud water and ice from cumulus  updraft at the top !
                   ! interface.  Note that this is not 'detrainment' term  but a 'sink' !
                   ! term of cumulus updraft qt ( or one part of 'source' term of  mean !
                   ! environmental qt ). At this stage, as the most simplest choice, if !
                   ! condensate amount within cumulus updraft is larger than a critical !
                   ! value, 'criqc', expels the surplus condensate from cumulus updraft !
                   ! to the environment. A certain fraction ( e.g., 'frc_sus' ) of this !
                   ! expelled condesnate will be in a form that can be suspended in the !
                   ! layer k where it was formed, while the other fraction, '1-frc_sus' ! 
                   ! will be in a form of precipitatble (e.g.,can potentially fall down !
                   ! across the base interface of layer k ). In turn we should describe !
                   ! subsequent falling of precipitable condensate ('1-frc_sus') across !
                   ! the base interface of the layer k, &  evaporation of precipitating !
                   ! water in the below layer k-1 and associated evaporative cooling of !
                   ! the later, k-1, and falling of 'non-evaporated precipitating water !
                   ! ( which was initially formed in layer k ) and a newly-formed preci !
                   ! pitable water in the layer, k-1', across the base interface of the !
                   ! lower layer k-1.  Cloud microphysics should correctly describe all !
                   ! of these process.  In a near future, I should significantly modify !
                   ! this cloud microphysics, including precipitation-induced downdraft !
                   ! also.                                                              !
                   ! ------------------------------------------------------------------ !

                   CALL conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   IF( (qlj + qij) .GT. criqc ) THEN
                      exql    = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
                      exqi    = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
                      ! ---------------------------------------------------------------- !
                      ! It is very important to re-update 'qtu' and 'thlu'  at the upper ! 
                      ! interface after expelling condensate from cumulus updraft at the !
                      ! top interface of the layer. As mentioned above, this is a 'sink' !
                      ! of cumulus qt (or equivalently, a 'source' of environmentasl qt),!
                      ! not a regular convective'detrainment'.                           !
                      ! ---------------------------------------------------------------- !
                      qtu(k)  = qtu(k) - exql - exqi
                      thlu(k) = thlu(k) + (xlv/cp/exns0(k))*exql + (xls/cp/exns0(k))*exqi 
                      ! ---------------------------------------------------------------- !
                      ! Expelled cloud condensate into the environment from the updraft. ! 
                      ! After all the calculation later, 'dwten' and 'diten' will have a !
                      ! unit of [ kg/kg/s ], because it is a tendency of qt. Restoration !
                      ! of 'dwten' and 'diten' to this correct unit through  multiplying !
                      ! 'umf(k)*g/dp0(k)' will be performed later after finally updating !
                      ! 'umf' using a 'rmaxfrac' constraint near the end of this updraft !
                      ! buoyancy sorting loop.                                           !
                      ! ---------------------------------------------------------------- !
                      dwten(k) = exql   
                      diten(k) = exqi
                   ELSE
                      dwten(k) = 0._r8
                      diten(k) = 0._r8
                   ENDIF
                   ! ----------------------------------------------------------------- ! 
                   ! Update 'thvu(k)' after detraining condensate from cumulus updraft.!
                   ! ----------------------------------------------------------------- ! 
                   CALL conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
                   IF( id_check .EQ. 1 ) THEN
                      exit_conden(i) = 1._r8
                      id_exit = .TRUE.
                      go to 333
                   END IF
                   thvu(k) = thj * ( 1._r8 + zvir * qvj - qlj - qij )

                   ! ----------------------------------------------------------- ! 
                   ! Calculate updraft vertical velocity at the upper interface. !
                   ! In order to calculate 'wtw' at the upper interface, we use  !
                   ! 'wtw' at the lower interface. Note  'wtw'  is continuously  ! 
                   ! updated as cumulus updraft rises.                           !
                   ! ----------------------------------------------------------- !

                   bogbot = rbuoy * ( thvu(km1) / thvebot  - 1._r8 ) ! Cloud buoyancy at base interface
                   bogtop = rbuoy * ( thvu(k) / thv0top(k) - 1._r8 ) ! Cloud buoyancy at top  interface

                   delbog = bogtop - bogbot
                   drage  = fer(k) * ( 1._r8 + rdrag )
                   expfac = EXP(-2._r8*drage*dpe)

                   wtwb = wtw
                   IF( drage*dpe .GT. 1.e-3_r8 ) THEN
                      wtw = wtw*expfac + (delbog + (1._r8-expfac)*(bogbot + delbog/(-2._r8*drage*dpe)))/(rho0j*drage)
                   ELSE
                      wtw = wtw + dpe * ( bogbot + bogtop ) / rho0j
                   ENDIF

                   ! Force the plume rise at least to klfc of the undiluted plume.
                   ! Because even the below is not complete, I decided not to include this.

                   ! if( k .le. klfc ) then
                   !     wtw = max( 1.e-2_r8, wtw )
                   ! endif 

                   ! -------------------------------------------------------------- !
                   ! Repeat 'iter_xc' iteration loop until 'iter_xc = niter_xc'.    !
                   ! Also treat the case even when wtw < 0 at the 'kpen' interface. !
                   ! -------------------------------------------------------------- !  

                   IF( wtw .GT. 0._r8 ) THEN   
                      thlue = 0.5_r8 * ( thlu(km1) + thlu(k) )
                      qtue  = 0.5_r8 * ( qtu(km1)  +  qtu(k) )         
                      wue   = 0.5_r8 *   SQRT( MAX( wtwb + wtw, 0._r8 ) )
                   ELSE
                      go to 111
                   ENDIF

                ENDDO ! End of 'iter_xc' loop  

111             CONTINUE

                ! --------------------------------------------------------------------------- ! 
                ! Add the contribution of self-detrainment  to vertical variations of cumulus !
                ! updraft mass flux. The reason why we are trying to include self-detrainment !
                ! is as follows.  In current scheme,  vertical variation of updraft mass flux !
                ! is not fully consistent with the vertical variation of updraft vertical w.  !
                ! For example, within a given layer, let's assume that  cumulus w is positive !
                ! at the base interface, while negative at the top interface. This means that !
                ! cumulus updraft cannot reach to the top interface of the layer. However,    !
                ! cumulus updraft mass flux at the top interface is not zero according to the !
                ! vertical tendency equation of cumulus mass flux.   Ideally, cumulus updraft ! 
                ! mass flux at the top interface should be zero for this case. In order to    !
                ! assures that cumulus updraft mass flux goes to zero when cumulus updraft    ! 
                ! vertical velocity goes to zero, we are imposing self-detrainment term as    !
                ! below by considering layer-mean cloud buoyancy and cumulus updraft vertical !
                ! velocity square at the top interface. Use of auto-detrainment term will  be !
                ! determined by setting 'use_self_detrain=.true.' in the parameter sentence.  !
                ! --------------------------------------------------------------------------- !

                IF( use_self_detrain ) THEN
                   autodet = MIN( 0.5_r8*g*(bogbot+bogtop)/(MAX(wtw,0._r8)+1.e-4_r8), 0._r8 ) 
                   umf(k)  = umf(k) * EXP( 0.637_r8*(dpe/rho0j/g) * autodet )   
                END IF
                IF( umf(k) .EQ. 0._r8 ) wtw = -1._r8

                ! -------------------------------------- !
                ! Below block is just a dignostic output !
                ! -------------------------------------- ! 

                excessu_arr(k) = excessu
                excess0_arr(k) = excess0
                xc_arr(k)      = xc
                aquad_arr(k)   = aquad
                bquad_arr(k)   = bquad
                cquad_arr(K)   = cquad
                bogbot_arr(k)  = bogbot
                bogtop_arr(k)  = bogtop

                ! ------------------------------------------------------------------- !
                ! 'kbup' is the upper most layer in which cloud buoyancy  is positive ! 
                ! both at the base and top interface.  'kpen' is the upper most layer !
                ! up to cumulus can reach. Usually, 'kpen' is located higher than the !
                ! 'kbup'. Note we initialized these by 'kbup = krel' & 'kpen = krel'. !
                ! As explained before, it is possible that only 'kpen' is updated,    !
                ! while 'kbup' keeps its initialization value. For this case, current !
                ! scheme will simply turns-off penetrative entrainment fluxes and use ! 
                ! normal buoyancy-sorting fluxes for 'kbup <= k <= kpen-1' interfaces,!
                ! in order to describe shallow continental cumulus convection.        !
                ! ------------------------------------------------------------------- !

                ! if( bogbot .gt. 0._r8 .and. bogtop .gt. 0._r8 ) then 
                ! if( bogtop .gt. 0._r8 ) then          
                IF( bogtop .GT. 0._r8 .AND. wtw .GT. 0._r8 ) THEN 
                   kbup = k
                END IF

                IF( wtw .LE. 0._r8 ) THEN
                   kpen = k
                   go to 45
                END IF

                wu(k) = SQRT(wtw)
                IF( wu(k) .GT. 100._r8 ) THEN
                   exit_wu(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                ENDIF

                ! ---------------------------------------------------------------------------- !
                ! Iteration end due to 'rmaxfrac' constraint [ ***************************** ] ! 
                ! ---------------------------------------------------------------------------- !

                ! ---------------------------------------------------------------------- !
                ! Calculate updraft fractional area at the upper interface and set upper ! 
                ! limit to 'ufrc' by 'rmaxfrac'. In order to keep the consistency  among !
                ! ['ufrc','umf','wu (or wtw)'], if ufrc is limited by 'rmaxfrac', either !
                ! 'umf' or 'wu' should be changed. Although both 'umf' and 'wu (wtw)' at !
                ! the current upper interface are used for updating 'umf' & 'wu'  at the !
                ! next upper interface, 'umf' is a passive variable not influencing  the !
                ! buoyancy sorting process in contrast to 'wtw'. This is a reason why we !
                ! adjusted 'umf' instead of 'wtw'. In turn we updated 'fdr' here instead !
                ! of 'fer',  which guarantees  that all previously updated thermodynamic !
                ! variables at the upper interface before applying 'rmaxfrac' constraint !
                ! are already internally consistent,  even though 'ufrc'  is  limited by !
                ! 'rmaxfrac'. Thus, we don't need to go through interation loop again.If !
                ! If we update 'fer' however, we should go through above iteration loop. !
                ! ---------------------------------------------------------------------- !

                rhos0j  = ps0(k) / ( r * 0.5_r8 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )
                ufrc(k) = umf(k) / ( rhos0j * wu(k) )
                IF( ufrc(k) .GT. rmaxfrac ) THEN
                   limit_ufrc(i) = 1._r8 
                   ufrc(k) = rmaxfrac
                   umf(k)  = rmaxfrac * rhos0j * wu(k)
                   fdr(k)  = fer(k) - LOG( umf(k) / umf(km1) ) / dpe
                ENDIF

                ! ------------------------------------------------------------ !
                ! Update environmental properties for at the mid-point of next !
                ! upper layer for use in buoyancy sorting.                     !
                ! ------------------------------------------------------------ ! 

                pe      = p0(k+1)
                dpe     = dp0(k+1)
                exne    = exn0(k+1)
                thvebot = thv0bot(k+1)
                thle    = thl0(k+1)
                qte     = qt0(k+1)
                ue      = u0(k+1)
                ve      = v0(k+1) 
                DO m = 1, ncnst
                   tre(m)  = tr0(k+1,m)
                ENDDO

             END DO   ! End of cumulus updraft loop from the 'krel' layer to 'kpen' layer.

             ! ------------------------------------------------------------------------------- !
             ! Up to this point, we finished all of buoyancy sorting processes from the 'krel' !
             ! layer to 'kpen' layer: at the top interface of individual layers, we calculated !
             ! updraft and penetrative mass fluxes [ umf(k) & emf(k) = 0 ], updraft fractional !
             ! area [ ufrc(k) ],  updraft vertical velocity [ wu(k) ],  updraft  thermodynamic !
             ! variables [thlu(k),qtu(k),uu(k),vu(k),thvu(k)]. In the layer,we also calculated !
             ! fractional entrainment-detrainment rate [ fer(k), fdr(k) ], and detrainment ten !
             ! dency of water and ice from cumulus updraft [ dwten(k), diten(k) ]. In addition,!
             ! we updated and identified 'krel' and 'kpen' layer index, if any.  In the 'kpen' !
             ! layer, we calculated everything mentioned above except the 'wu(k)' and 'ufrc(k)'!
             ! since a real value of updraft vertical velocity is not defined at the kpen  top !
             ! interface (note 'ufrc' at the top interface of layer is calculated from 'umf(k)'!
             ! and 'wu(k)'). As mentioned before, special treatment is required when 'kbup' is !
             ! not updated and so 'kbup = krel'.                                               !
             ! ------------------------------------------------------------------------------- !

             ! ------------------------------------------------------------------------------ !
             ! During the 'iter_scaleh' iteration loop, non-physical ( with non-zero values ) !
             ! values can remain in the variable arrays above (also 'including' in case of wu !
             ! and ufrc at the top interface) the 'kpen' layer. This can happen when the kpen !
             ! layer index identified from the 'iter_scaleh = 1' iteration loop is located at !
             ! above the kpen layer index identified from   'iter_scaleh = 3' iteration loop. !
             ! Thus, in the following calculations, we should only use the values in each     !
             ! variables only up to finally identified 'kpen' layer & 'kpen' interface except ! 
             ! 'wu' and 'ufrc' at the top interface of 'kpen' layer.    Note that in order to !
             ! prevent any problems due to these non-physical values, I re-initialized    the !
             ! values of [ umf(kpen:mkx), emf(kpen:mkx), dwten(kpen+1:mkx), diten(kpen+1:mkx),! 
             ! fer(kpen:mkx), fdr(kpen+1:mkx), ufrc(kpen:mkx) ] to be zero after 'iter_scaleh'!
             ! do loop.                                                                       !
             ! ------------------------------------------------------------------------------ !

45           CONTINUE

             ! ------------------------------------------------------------------------------ !
             ! Calculate 'ppen( < 0 )', updarft penetrative distance from the lower interface !
             ! of 'kpen' layer. Note that bogbot & bogtop at the 'kpen' layer either when fer !
             ! is zero or non-zero was already calculated above.                              !
             ! It seems that below qudarature solving formula is valid only when bogbot < 0.  !
             ! Below solving equation is clearly wrong ! I should revise this !               !
             ! ------------------------------------------------------------------------------ ! 

             IF( drage .EQ. 0._r8 ) THEN
                aquad =  ( bogtop - bogbot ) / ( ps0(kpen) - ps0(kpen-1) )
                bquad =  2._r8 * bogbot
                cquad = -wu(kpen-1)**2 * rho0j
                CALL roots(aquad,bquad,cquad,xc1,xc2,status)
                IF( status .EQ. 0 ) THEN
                   IF( xc1 .LE. 0._r8 .AND. xc2 .LE. 0._r8 ) THEN
                      ppen = MAX( xc1, xc2 )
                      ppen = MIN( 0._r8,MAX( -dp0(kpen), ppen ) )  
                   ELSEIF( xc1 .GT. 0._r8 .AND. xc2 .GT. 0._r8 ) THEN
                      ppen = -dp0(kpen)
                      WRITE(iulog,*) 'Warning : UW-Cumulus penetrates upto kpen interface'
                   ELSE
                      ppen = MIN( xc1, xc2 )
                      ppen = MIN( 0._r8,MAX( -dp0(kpen), ppen ) )  
                   ENDIF
                ELSE
                   ppen = -dp0(kpen)
                   WRITE(iulog,*) 'Warning : UW-Cumulus penetrates upto kpen interface'
                ENDIF
             ELSE 
                ppen = compute_ppen(wtwb,drage,bogbot,bogtop,rho0j,dp0(kpen))
             ENDIF
             IF( ppen .EQ. -dp0(kpen) .OR. ppen .EQ. 0._r8 ) limit_ppen(i) = 1._r8

             ! -------------------------------------------------------------------- !
             ! Re-calculate the amount of expelled condensate from cloud updraft    !
             ! at the cumulus top. This is necessary for refined calculations of    !
             ! bulk cloud microphysics at the cumulus top. Note that ppen < 0._r8   !
             ! In the below, I explicitly calculate 'thlu_top' & 'qtu_top' by       !
             ! using non-zero 'fer(kpen)'.                                          !    
             ! -------------------------------------------------------------------- !

             IF( fer(kpen)*(-ppen) .LT. 1.e-4_r8 ) THEN
                thlu_top = thlu(kpen-1) + ( thl0(kpen) + ssthl0(kpen) * (-ppen) / 2._r8 - thlu(kpen-1) ) * &
                fer(kpen) * (-ppen)
                qtu_top  =  qtu(kpen-1) + (  qt0(kpen) +  ssqt0(kpen) * (-ppen) / 2._r8  - qtu(kpen-1) ) * &
                fer(kpen) * (-ppen)
             ELSE
                thlu_top = ( thl0(kpen) + ssthl0(kpen) / fer(kpen) - ssthl0(kpen) * (-ppen) / 2._r8 ) - &
                     ( thl0(kpen) + ssthl0(kpen) * (-ppen) / 2._r8 - thlu(kpen-1) + ssthl0(kpen) / fer(kpen) ) &
                      * EXP(-fer(kpen) * (-ppen))
                qtu_top  = ( qt0(kpen)  +  ssqt0(kpen) / fer(kpen) -  ssqt0(kpen) * (-ppen) / 2._r8 ) - &  
                     ( qt0(kpen)  +  ssqt0(kpen) * (-ppen) / 2._r8 -  qtu(kpen-1) +  ssqt0(kpen) / fer(kpen) ) * &
                     EXP(-fer(kpen) * (-ppen))
             END IF

             CALL conden(ps0(kpen-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check,qsat)
             IF( id_check .EQ. 1 ) THEN
                exit_conden(i) = 1._r8
                id_exit = .TRUE.
                go to 333
             END IF
             exntop = ((ps0(kpen-1)+ppen)/p00)**rovcp
             IF( (qlj + qij) .GT. criqc ) THEN
                dwten(kpen) = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
                diten(kpen) = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
                qtu_top  = qtu_top - dwten(kpen) - diten(kpen)
                thlu_top = thlu_top + (xlv/cp/exntop)*dwten(kpen) + (xls/cp/exntop)*diten(kpen) 
             ELSE
                dwten(kpen) = 0._r8
                diten(kpen) = 0._r8
             ENDIF

             ! ----------------------------------------------------------------------- !
             ! Calculate cumulus scale height as the top height that cumulus can reach.!
             ! ----------------------------------------------------------------------- !

             rhos0j = ps0(kpen-1)/(r*0.5_r8*(thv0bot(kpen)+thv0top(kpen-1))*exns0(kpen-1))  
             cush   = zs0(kpen-1) - ppen/rhos0j/g
             scaleh = cush 

          END DO   ! End of 'iter_scaleh' loop.   

          ! -------------------------------------------------------------------- !   
          ! The 'forcedCu' is logical identifier saying whether cumulus updraft  !
          ! overcome the buoyancy barrier just above the PBL top. If it is true, !
          ! cumulus did not overcome the barrier -  this is a shallow convection !
          ! with negative cloud buoyancy, mimicking  shallow continental cumulus !
          ! convection. Depending on 'forcedCu' parameter, treatment of heat  &  !
          ! moisture fluxes at the entraining interfaces, 'kbup <= k < kpen - 1' !
          ! will be set up in a different ways, as will be shown later.          !
          ! -------------------------------------------------------------------- !

          IF( kbup .EQ. krel ) THEN 
             forcedCu = .TRUE.
             limit_shcu(i) = 1._r8
          ELSE
             forcedCu = .FALSE.
             limit_shcu(i) = 0._r8
          ENDIF

          ! ------------------------------------------------------------------ !
          ! Filtering of unerasonable cumulus adjustment here.  This is a very !
          ! important process which should be done cautiously. Various ways of !
          ! filtering are possible depending on cases mainly using the indices !
          ! of key layers - 'klcl','kinv','krel','klfc','kbup','kpen'. At this !
          ! stage, the followings are all possible : 'kinv >= 2', 'klcl >= 1', !
          ! 'krel >= kinv', 'kbup >= krel', 'kpen >= krel'. I must design this !
          ! filtering very cautiously, in such that none of  realistic cumulus !
          ! convection is arbitrarily turned-off. Potentially, I might turn-off! 
          ! cumulus convection if layer-mean 'ql > 0' in the 'kinv-1' layer,in !
          ! order to suppress cumulus convection growing, based at the Sc top. ! 
          ! This is one of potential future modifications. Note that ppen < 0. !
          ! ------------------------------------------------------------------ !

          cldhgt = ps0(kpen-1) + ppen
          IF( forcedCu ) THEN
             ! write(iulog,*) 'forcedCu - did not overcome initial buoyancy barrier'
             exit_cufilter(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          END IF
          ! Limit 'additional shallow cumulus' for DYCOMS simulation.
          ! if( cldhgt.ge.88000._r8 ) then
          !     id_exit = .true.
          !     go to 333
          ! end if

          ! ------------------------------------------------------------------------------ !
          ! Re-initializing some key variables above the 'kpen' layer in order to suppress !
          ! the influence of non-physical values above 'kpen', in association with the use !
          ! of 'iter_scaleh' loop. Note that umf, emf,  ufrc are defined at the interfaces !
          ! (0:mkx), while 'dwten','diten', 'fer', 'fdr' are defined at layer mid-points.  !
          ! Initialization of 'fer' and 'fdr' is for correct writing purpose of diagnostic !
          ! output. Note that we set umf(kpen)=emf(kpen)=ufrc(kpen)=0, in consistent  with !
          ! wtw < 0  at the top interface of 'kpen' layer. However, we still have non-zero !
          ! expelled cloud condensate in the 'kpen' layer.                                 !
          ! ------------------------------------------------------------------------------ !

          umf(kpen:mkx)     = 0._r8
          emf(kpen:mkx)     = 0._r8
          ufrc(kpen:mkx)    = 0._r8
          dwten(kpen+1:mkx) = 0._r8
          diten(kpen+1:mkx) = 0._r8
          fer(kpen+1:mkx)   = 0._r8
          fdr(kpen+1:mkx)   = 0._r8

          ! ------------------------------------------------------------------------ !
          ! Calculate downward penetrative entrainment mass flux, 'emf(k) < 0',  and !
          ! thermodynamic properties of penetratively entrained airs at   entraining !
          ! interfaces. emf(k) is defined from the top interface of the  layer  kbup !
          ! to the bottom interface of the layer 'kpen'. Note even when  kbup = krel,!
          ! i.e.,even when 'kbup' was not updated in the above buoyancy  sorting  do !
          ! loop (i.e., 'kbup' remains as the initialization value),   below do loop !
          ! of penetrative entrainment flux can be performed without  any conceptual !
          ! or logical problems, because we have already computed all  the variables !
          ! necessary for performing below penetrative entrainment block.            !
          ! In the below 'do' loop, 'k' is an interface index at which non-zero 'emf'! 
          ! (penetrative entrainment mass flux) is calculated. Since cumulus updraft !
          ! is negatively buoyant in the layers between the top interface of 'kbup'  !
          ! layer (interface index, kbup) and the top interface of 'kpen' layer, the !
          ! fractional lateral entrainment, fer(k) within these layers will be close !
          ! to zero - so it is likely that only strong lateral detrainment occurs in !
          ! thses layers. Under this situation,we can easily calculate the amount of !
          ! detrainment cumulus air into these negatively buoyanct layers by  simply !
          ! comparing cumulus updraft mass fluxes between the base and top interface !
          ! of each layer: emf(k) = emf(k-1)*exp(-fdr(k)*dp0(k))                     !
          !                       ~ emf(k-1)*(1-rei(k)*dp0(k))                       !
          !                emf(k-1)-emf(k) ~ emf(k-1)*rei(k)*dp0(k)                  !
          ! Current code assumes that about 'rpen~10' times of these detrained  mass !
          ! are penetratively re-entrained down into the 'k-1' interface. And all of !
          ! these detrained masses are finally dumped down into the top interface of !
          ! 'kbup' layer. Thus, the amount of penetratively entrained air across the !
          ! top interface of 'kbup' layer with 'rpen~10' becomes too large.          !
          ! Note that this penetrative entrainment part can be completely turned-off !
          ! and we can simply use normal buoyancy-sorting involved turbulent  fluxes !
          ! by modifying 'penetrative entrainment fluxes' part below.                !
          ! ------------------------------------------------------------------------ !

          ! -----------------------------------------------------------------------!
          ! Calculate entrainment mass flux and conservative scalars of entraining !
          ! free air at interfaces of 'kbup <= k < kpen - 1'                       !
          ! ---------------------------------------------------------------------- !

          DO k = 0, mkx
             thlu_emf(k) = thlu(k)
             qtu_emf(k)  = qtu(k)
             uu_emf(k)   = uu(k)
             vu_emf(k)   = vu(k)
             DO m = 1, ncnst
                tru_emf(k,m)  = tru(k,m)
             ENDDO
          END DO

          DO k = kpen - 1, kbup, -1  ! Here, 'k' is an interface index at which
             ! penetrative entrainment fluxes are calculated. 

             rhos0j = ps0(k) / ( r * 0.5_r8 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )

             IF( k .EQ. kpen - 1 ) THEN

                ! ------------------------------------------------------------------------ ! 
                ! Note that 'ppen' has already been calculated in the above 'iter_scaleh'  !
                ! loop assuming zero lateral entrainmentin the layer 'kpen'.               !
                ! ------------------------------------------------------------------------ !       

                ! -------------------------------------------------------------------- !
                ! Calculate returning mass flux, emf ( < 0 )                           !
                ! Current penetrative entrainment rate with 'rpen~10' is too large and !
                ! future refinement is necessary including the definition of 'thl','qt'! 
                ! of penetratively entrained air.  Penetratively entrained airs across !
                ! the 'kpen-1' interface is assumed to have the properties of the base !
                ! interface of 'kpen' layer. Note that 'emf ~ - umf/ufrc = - w * rho'. !
                ! Thus, below limit sets an upper limit of |emf| to be ~ 10cm/s, which !
                ! is very loose constraint. Here, I used more restricted constraint on !
                ! the limit of emf, assuming 'emf' cannot exceed a net mass within the !
                ! layer above the interface. Similar to the case of warming and drying !
                ! due to cumulus updraft induced compensating subsidence,  penetrative !
                ! entrainment induces compensating upwelling -     in order to prevent !  
                ! numerical instability in association with compensating upwelling, we !
                ! should similarily limit the amount of penetrative entrainment at the !
                ! interface by the amount of masses within the layer just above the    !
                ! penetratively entraining interface.                                  !
                ! -------------------------------------------------------------------- !

                IF( ( umf(k)*ppen*rei(kpen)*rpen ) .LT. -0.1_r8*rhos0j )         limit_emf(i) = 1._r8
                IF( ( umf(k)*ppen*rei(kpen)*rpen ) .LT. -0.9_r8*dp0(kpen)/g/dt ) limit_emf(i) = 1._r8             

                emf(k) = MAX( MAX( umf(k)*ppen*rei(kpen)*rpen, -0.1_r8*rhos0j), -0.9_r8*dp0(kpen)/g/dt)
                thlu_emf(k) = thl0(kpen) + ssthl0(kpen) * ( ps0(k) - p0(kpen) )
                qtu_emf(k)  = qt0(kpen)  + ssqt0(kpen)  * ( ps0(k) - p0(kpen) )
                uu_emf(k)   = u0(kpen)   + ssu0(kpen)   * ( ps0(k) - p0(kpen) )     
                vu_emf(k)   = v0(kpen)   + ssv0(kpen)   * ( ps0(k) - p0(kpen) )   
                DO m = 1, ncnst
                   tru_emf(k,m)  = tr0(kpen,m)  + sstr0(kpen,m)  * ( ps0(k) - p0(kpen) )
                ENDDO

             ELSE ! if(k.lt.kpen-1). 

                ! --------------------------------------------------------------------------- !
                ! Note we are coming down from the higher interfaces to the lower interfaces. !
                ! Also note that 'emf < 0'. So, below operation is a summing not subtracting. !
                ! In order to ensure numerical stability, I imposed a modified correct limit  ! 
                ! of '-0.9*dp0(k+1)/g/dt' on emf(k).                                          !
                ! --------------------------------------------------------------------------- !

                IF( use_cumpenent ) THEN  ! Original Cumulative Penetrative Entrainment

                   IF( ( emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen ) .LT. -0.1_r8*rhos0j )        limit_emf(i) = 1
                   IF( ( emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen ) .LT. -0.9_r8*dp0(k+1)/g/dt ) limit_emf(i) = 1         
                   emf(k) = MAX(MAX(emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen, -0.1_r8*rhos0j), -0.9_r8*dp0(k+1)/g/dt )    
                   IF( ABS(emf(k)) .GT. ABS(emf(k+1)) ) THEN
                      thlu_emf(k) = ( thlu_emf(k+1) * emf(k+1) + thl0(k+1) * ( emf(k) - emf(k+1) ) ) / emf(k)
                      qtu_emf(k)  = ( qtu_emf(k+1)  * emf(k+1) + qt0(k+1)  * ( emf(k) - emf(k+1) ) ) / emf(k)
                      uu_emf(k)   = ( uu_emf(k+1)   * emf(k+1) + u0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
                      vu_emf(k)   = ( vu_emf(k+1)   * emf(k+1) + v0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
                      DO m = 1, ncnst
                         tru_emf(k,m)  = ( tru_emf(k+1,m)  * emf(k+1) + tr0(k+1,m)  * ( emf(k) - emf(k+1) ) ) / emf(k)
                      ENDDO
                   ELSE   
                      thlu_emf(k) = thl0(k+1)
                      qtu_emf(k)  =  qt0(k+1)
                      uu_emf(k)   =   u0(k+1)
                      vu_emf(k)   =   v0(k+1)
                      DO m = 1, ncnst
                         tru_emf(k,m)  =  tr0(k+1,m)
                      ENDDO
                   ENDIF

                ELSE ! Alternative Non-Cumulative Penetrative Entrainment

                   IF( ( -umf(k)*dp0(k+1)*rei(k+1)*rpen ) .LT. -0.1_r8*rhos0j )        limit_emf(i) = 1
                   IF( ( -umf(k)*dp0(k+1)*rei(k+1)*rpen ) .LT. -0.9_r8*dp0(k+1)/g/dt ) limit_emf(i) = 1         
                   emf(k) = MAX(MAX(-umf(k)*dp0(k+1)*rei(k+1)*rpen, -0.1_r8*rhos0j), -0.9_r8*dp0(k+1)/g/dt )    
                   thlu_emf(k) = thl0(k+1)
                   qtu_emf(k)  =  qt0(k+1)
                   uu_emf(k)   =   u0(k+1)
                   vu_emf(k)   =   v0(k+1)
                   DO m = 1, ncnst
                      tru_emf(k,m)  =  tr0(k+1,m)
                   ENDDO

                ENDIF

             ENDIF

             ! ---------------------------------------------------------------------------- !
             ! In this GCM modeling framework,  all what we should do is to calculate  heat !
             ! and moisture fluxes at the given geometrically-fixed height interfaces -  we !
             ! don't need to worry about movement of material height surface in association !
             ! with compensating subsidence or unwelling, in contrast to the bulk modeling. !
             ! In this geometrically fixed height coordinate system, heat and moisture flux !
             ! at the geometrically fixed height handle everything - a movement of material !
             ! surface is implicitly treated automatically. Note that in terms of turbulent !
             ! heat and moisture fluxes at model interfaces, both the cumulus updraft  mass !
             ! flux and penetratively entraining mass flux play the same role -both of them ! 
             ! warms and dries the 'kbup' layer, cools and moistens the 'kpen' layer,   and !
             ! cools and moistens any intervening layers between 'kbup' and 'kpen' layers.  !
             ! It is important to note these identical roles on turbulent heat and moisture !
             ! fluxes of 'umf' and 'emf'.                                                   !
             ! When 'kbup' is a stratocumulus-topped PBL top interface,  increase of 'rpen' !
             ! is likely to strongly diffuse stratocumulus top interface,  resulting in the !
             ! reduction of cloud fraction. In this sense, the 'kbup' interface has a  very !
             ! important meaning and role : across the 'kbup' interface, strong penetrative !
             ! entrainment occurs, thus any sharp gradient properties across that interface !
             ! are easily diffused through strong mass exchange. Thus, an initialization of ! 
             ! 'kbup' (and also 'kpen') should be done very cautiously as mentioned before. ! 
             ! In order to prevent this stron diffusion for the shallow cumulus convection  !
             ! based at the Sc top, it seems to be good to initialize 'kbup = krel', rather !
             ! that 'kbup = krel-1'.                                                        !
             ! ---------------------------------------------------------------------------- !

          END DO

          !------------------------------------------------------------------ !
          !                                                                   ! 
          ! Compute turbulent heat, moisture, momentum flux at all interfaces !
          !                                                                   !
          !------------------------------------------------------------------ !
          ! It is very important to note that in calculating turbulent fluxes !
          ! below, we must not double count turbulent flux at any interefaces.!
          ! In the below, turbulent fluxes at the interfaces (interface index !
          ! k) are calculated by the following 4 blocks in consecutive order: !
          !                                                                   !
          ! (1) " 0 <= k <= kinv - 1 "  : PBL fluxes.                         !
          !     From 'fluxbelowinv' using reconstructed PBL height. Currently,!
          !     the reconstructed PBLs are independently calculated for  each !
          !     individual conservative scalar variables ( qt, thl, u, v ) in !
          !     each 'fluxbelowinv',  instead of being uniquely calculated by !
          !     using thvl. Turbulent flux at the surface is assumed to be 0. !
          ! (2) " kinv <= k <= krel - 1 " : Non-buoyancy sorting fluxes       !
          !     Assuming cumulus mass flux  and cumulus updraft thermodynamic !
          !     properties (except u, v which are modified by the PGFc during !
          !     upward motion) are conserved during a updraft motion from the !
          !     PBL top interface to the release level. If these layers don't !
          !     exist (e,g, when 'krel = kinv'), then  current routine do not !
          !     perform this routine automatically. So I don't need to modify !
          !     anything.                                                     ! 
          ! (3) " krel <= k <= kbup - 1 " : Buoyancy sorting fluxes           !
          !     From laterally entraining-detraining buoyancy sorting plumes. ! 
          ! (4) " kbup <= k < kpen-1 " : Penetrative entrainment fluxes       !
          !     From penetratively entraining plumes,                         !
          !                                                                   !
          ! In case of normal situation, turbulent interfaces  in each groups !
          ! are mutually independent of each other. Thus double flux counting !
          ! or ambiguous flux counting requiring the choice among the above 4 !
          ! groups do not occur normally. However, in case that cumulus plume !
          ! could not completely overcome the buoyancy barrier just above the !
          ! PBL top interface and so 'kbup = krel' (.forcedCu=.true.) ( here, !
          ! it can be either 'kpen = krel' as the initialization, or ' kpen > !
          ! krel' if cumulus updraft just penetrated over the top of  release !
          ! layer ). If this happens, we should be very careful in organizing !
          ! the sequence of the 4 calculation routines above -  note that the !
          ! routine located at the later has the higher priority.  Additional ! 
          ! feature I must consider is that when 'kbup = kinv - 1' (this is a !
          ! combined situation of 'kbup=krel-1' & 'krel = kinv' when I  chose !
          ! 'kbup=krel-1' instead of current choice of 'kbup=krel'), a strong !
          ! penetrative entrainment fluxes exists at the PBL top interface, & !
          ! all of these fluxes are concentrated (deposited) within the layer ! 
          ! just below PBL top interface (i.e., 'kinv-1' layer). On the other !
          ! hand, in case of 'fluxbelowinv', only the compensating subsidence !
          ! effect is concentrated in the 'kinv-1' layer and 'pure' turbulent !
          ! heat and moisture fluxes ( 'pure' means the fluxes not associated !
          ! with compensating subsidence) are linearly distributed throughout !
          ! the whole PBL. Thus different choice of the above flux groups can !
          ! produce very different results. Output variable should be written !
          ! consistently to the choice of computation sequences.              !
          ! When the case of 'kbup = krel(-1)' happens,another way to dealing !
          ! with this case is to simply ' exit ' the whole cumulus convection !
          ! calculation without performing any cumulus convection.     We can !
          ! choose this approach by specifying a condition in the  'Filtering !
          ! of unreasonable cumulus adjustment' just after 'iter_scaleh'. But !
          ! this seems not to be a good choice (although this choice was used !
          ! previous code ), since it might arbitrary damped-out  the shallow !
          ! cumulus convection over the continent land, where shallow cumulus ! 
          ! convection tends to be negatively buoyant.                        !
          ! ----------------------------------------------------------------- !  

          ! --------------------------------------------------- !
          ! 1. PBL fluxes :  0 <= k <= kinv - 1                 !
          !    All the information necessary to reconstruct PBL ! 
          !    height are passed to 'fluxbelowinv'.             !
          ! --------------------------------------------------- !

          xsrc  = qtsrc
          xmean = qt0(kinv)
          xtop  = qt0(kinv+1) + ssqt0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
          xbot  = qt0(kinv-1) + ssqt0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )        
          CALL fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
          qtflx(0:kinv-1) = xflx(0:kinv-1)

          xsrc  = thlsrc
          xmean = thl0(kinv)
          xtop  = thl0(kinv+1) + ssthl0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
          xbot  = thl0(kinv-1) + ssthl0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )        
          CALL fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
          slflx(0:kinv-1) = cp * exns0(0:kinv-1) * xflx(0:kinv-1)

          xsrc  = usrc
          xmean = u0(kinv)
          xtop  = u0(kinv+1) + ssu0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
          xbot  = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
          CALL fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
          uflx(0:kinv-1) = xflx(0:kinv-1)

          xsrc  = vsrc
          xmean = v0(kinv)
          xtop  = v0(kinv+1) + ssv0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
          xbot  = v0(kinv-1) + ssv0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
          CALL fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
          vflx(0:kinv-1) = xflx(0:kinv-1)

          DO m = 1, ncnst
             xsrc  = trsrc(m)
             xmean = tr0(kinv,m)
             xtop  = tr0(kinv+1,m) + sstr0(kinv+1,m) * ( ps0(kinv)   - p0(kinv+1) )
             xbot  = tr0(kinv-1,m) + sstr0(kinv-1,m) * ( ps0(kinv-1) - p0(kinv-1) )        
             CALL fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
             trflx(0:kinv-1,m) = xflx(0:kinv-1)
          ENDDO

          ! -------------------------------------------------------------- !
          ! 2. Non-buoyancy sorting fluxes : kinv <= k <= krel - 1         !
          !    Note that when 'krel = kinv', below block is never executed !
          !    as in a desirable, expected way ( but I must check  if this !
          !    is the case ). The non-buoyancy sorting fluxes are computed !
          !    only when 'krel > kinv'.                                    !
          ! -------------------------------------------------------------- !          

          uplus = 0._r8
          vplus = 0._r8
          DO k = kinv, krel - 1
             kp1 = k + 1
             qtflx(k) = cbmf * ( qtsrc  - (  qt0(kp1) +  ssqt0(kp1) * ( ps0(k) - p0(kp1) ) ) )          
             slflx(k) = cbmf * ( thlsrc - ( thl0(kp1) + ssthl0(kp1) * ( ps0(k) - p0(kp1) ) ) ) * cp * exns0(k)
             uplus    = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
             vplus    = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
             uflx(k)  = cbmf * ( usrc + uplus -  (  u0(kp1)  +   ssu0(kp1) * ( ps0(k) - p0(kp1) ) ) ) 
             vflx(k)  = cbmf * ( vsrc + vplus -  (  v0(kp1)  +   ssv0(kp1) * ( ps0(k) - p0(kp1) ) ) )
             DO m = 1, ncnst
                trflx(k,m) = cbmf * ( trsrc(m)  - (  tr0(kp1,m) +  sstr0(kp1,m) * ( ps0(k) - p0(kp1) ) ) )
             ENDDO
          END DO

          ! ------------------------------------------------------------------------ !
          ! 3. Buoyancy sorting fluxes : krel <= k <= kbup - 1                       !
          !    In case that 'kbup = krel - 1 ' ( or even in case 'kbup = krel' ),    ! 
          !    buoyancy sorting fluxes are not calculated, which is consistent,      !
          !    desirable feature.                                                    !  
          ! ------------------------------------------------------------------------ !

          DO k = krel, kbup - 1      
             kp1 = k + 1
             slflx(k) = cp * exns0(k) * umf(k) * ( thlu(k) - ( thl0(kp1) + ssthl0(kp1) * ( ps0(k) - p0(kp1) ) ) )
             qtflx(k) = umf(k) * ( qtu(k) - ( qt0(kp1) + ssqt0(kp1) * ( ps0(k) - p0(kp1) ) ) )
             uflx(k)  = umf(k) * ( uu(k) - ( u0(kp1) + ssu0(kp1) * ( ps0(k) - p0(kp1) ) ) )
             vflx(k)  = umf(k) * ( vu(k) - ( v0(kp1) + ssv0(kp1) * ( ps0(k) - p0(kp1) ) ) )
             DO m = 1, ncnst
                trflx(k,m) = umf(k) * ( tru(k,m) - ( tr0(kp1,m) + sstr0(kp1,m) * ( ps0(k) - p0(kp1) ) ) )
             ENDDO
          END DO

          ! ------------------------------------------------------------------------- !
          ! 4. Penetrative entrainment fluxes : kbup <= k <= kpen - 1                 !
          !    The only confliction that can happen is when 'kbup = kinv-1'. For this !
          !    case, turbulent flux at kinv-1 is calculated  both from 'fluxbelowinv' !
          !    and here as penetrative entrainment fluxes.  Since penetrative flux is !
          !    calculated later, flux at 'kinv - 1 ' will be that of penetrative flux.!
          !    However, turbulent flux calculated at 'kinv - 1' from penetrative entr.!
          !    is less attractable,  since more reasonable turbulent flux at 'kinv-1' !
          !    should be obtained from 'fluxbelowinv', by considering  re-constructed ! 
          !    inversion base height. This conflicting problem can be solved if we can!
          !    initialize 'kbup = krel', instead of kbup = krel - 1. This choice seems!
          !    to be more reasonable since it is not conflicted with 'fluxbelowinv' in!
          !    calculating fluxes at 'kinv - 1' ( for this case, flux at 'kinv-1' is  !
          !    always from 'fluxbelowinv' ), and flux at 'krel-1' is calculated from  !
          !    the non-buoyancy sorting flux without being competed with penetrative  !
          !    entrainment fluxes. Even when we use normal cumulus flux instead of    !
          !    penetrative entrainment fluxes at 'kbup <= k <= kpen-1' interfaces,    !
          !    the initialization of kbup=krel perfectly works without any conceptual !
          !    confliction. Thus it seems to be much better to choose 'kbup = krel'   !
          !    initialization of 'kbup', which is current choice.                     !
          !    Note that below formula uses conventional updraft cumulus fluxes for   !
          !    shallow cumulus which did not overcome the first buoyancy barrier above!
          !    PBL top while uses penetrative entrainment fluxes for the other cases  !
          !    'kbup <= k <= kpen-1' interfaces. Depending on cases, however, I can   !
          !    selelct different choice.                                              !
          ! ------------------------------------------------------------------------------------------------------------------ !
          !   if( forcedCu ) then                                                                                              !
          !       slflx(k) = cp * exns0(k) * umf(k) * ( thlu(k) - ( thl0(kp1) + ssthl0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
          !       qtflx(k) =                 umf(k) * (  qtu(k) - (  qt0(kp1) +  ssqt0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
          !       uflx(k)  =                 umf(k) * (   uu(k) - (   u0(kp1) +   ssu0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
          !       vflx(k)  =                 umf(k) * (   vu(k) - (   v0(kp1) +   ssv0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
          !       do m = 1, ncnst                                                                                              !
          !          trflx(k,m) = umf(k) * ( tru(k,m) - ( tr0(kp1,m) + sstr0(kp1,m) * ( ps0(k) - p0(kp1) ) ) )                 !
          !       enddo                                                                                                        !
          !   else                                                                                                             !
          !       slflx(k) = cp * exns0(k) * emf(k) * ( thlu_emf(k) - ( thl0(k) + ssthl0(k) * ( ps0(k) - p0(k) ) ) )           !
          !       qtflx(k) =                 emf(k) * (  qtu_emf(k) - (  qt0(k) +  ssqt0(k) * ( ps0(k) - p0(k) ) ) )           !
          !       uflx(k)  =                 emf(k) * (   uu_emf(k) - (   u0(k) +   ssu0(k) * ( ps0(k) - p0(k) ) ) )           !
          !       vflx(k)  =                 emf(k) * (   vu_emf(k) - (   v0(k) +   ssv0(k) * ( ps0(k) - p0(k) ) ) )           !
          !       do m = 1, ncnst                                                                                              !
          !          trflx(k,m) = emf(k) * ( tru_emf(k,m) - ( tr0(k,m) + sstr0(k,m) * ( ps0(k) - p0(k) ) ) )                   !
          !       enddo                                                                                                        !
          !   endif                                                                                                            !
          !                                                                                                                    !
          !   if( use_uppenent ) then ! Combined Updraft + Penetrative Entrainment Flux                                        !
          !       slflx(k) = cp * exns0(k) * umf(k) * ( thlu(k)     - ( thl0(kp1) + ssthl0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
          !                  cp * exns0(k) * emf(k) * ( thlu_emf(k) - (   thl0(k) +   ssthl0(k) * ( ps0(k) - p0(k) ) ) )       !
          !       qtflx(k) =                 umf(k) * (  qtu(k)     - (  qt0(kp1) +  ssqt0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
          !                                  emf(k) * (  qtu_emf(k) - (    qt0(k) +    ssqt0(k) * ( ps0(k) - p0(k) ) ) )       !                   
          !       uflx(k)  =                 umf(k) * (   uu(k)     - (   u0(kp1) +   ssu0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
          !                                  emf(k) * (   uu_emf(k) - (     u0(k) +     ssu0(k) * ( ps0(k) - p0(k) ) ) )       !                      
          !       vflx(k)  =                 umf(k) * (   vu(k)     - (   v0(kp1) +   ssv0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
          !                                  emf(k) * (   vu_emf(k) - (     v0(k) +     ssv0(k) * ( ps0(k) - p0(k) ) ) )       !                     
          !       do m = 1, ncnst                                                                                              !
          !          trflx(k,m) = umf(k) * ( tru(k,m) - ( tr0(kp1,m) + sstr0(kp1,m) * ( ps0(k) - p0(kp1) ) ) ) + &             ! 
          !                       emf(k) * ( tru_emf(k,m) - ( tr0(k,m) + sstr0(k,m) * ( ps0(k) - p0(k) ) ) )                   ! 
          !       enddo                                                                                                        !
          ! ------------------------------------------------------------------------------------------------------------------ !

          DO k = kbup, kpen - 1      
             kp1 = k + 1
             slflx(k) = cp * exns0(k) * emf(k) * ( thlu_emf(k) - ( thl0(k) + ssthl0(k) * ( ps0(k) - p0(k) ) ) )
             qtflx(k) =                 emf(k) * (  qtu_emf(k) - (  qt0(k) +  ssqt0(k) * ( ps0(k) - p0(k) ) ) ) 
             uflx(k)  =                 emf(k) * (   uu_emf(k) - (   u0(k) +   ssu0(k) * ( ps0(k) - p0(k) ) ) ) 
             vflx(k)  =                 emf(k) * (   vu_emf(k) - (   v0(k) +   ssv0(k) * ( ps0(k) - p0(k) ) ) )
             DO m = 1, ncnst
                trflx(k,m) = emf(k) * ( tru_emf(k,m) - ( tr0(k,m) + sstr0(k,m) * ( ps0(k) - p0(k) ) ) ) 
             ENDDO
          END DO

          ! ------------------------------------------- !
          ! Turn-off cumulus momentum flux as an option !
          ! ------------------------------------------- !

          IF( .NOT. use_momenflx ) THEN
             uflx(0:mkx) = 0._r8
             vflx(0:mkx) = 0._r8
          ENDIF

          ! -------------------------------------------------------- !
          ! Condensate tendency by compensating subsidence/upwelling !
          ! -------------------------------------------------------- !

          uemf(0:mkx)         = 0._r8
          DO k = 0, kinv - 2  ! Assume linear updraft mass flux within the PBL.
             uemf(k) = cbmf * ( ps0(0) - ps0(k) ) / ( ps0(0) - ps0(kinv-1) ) 
          END DO
          uemf(kinv-1:krel-1) = cbmf
          uemf(krel:kbup-1)   = umf(krel:kbup-1)
          uemf(kbup:kpen-1)   = emf(kbup:kpen-1) ! Only use penetrative entrainment flux consistently.

          comsub(1:mkx) = 0._r8
          DO k = 1, kpen
             comsub(k)  = 0.5_r8 * ( uemf(k) + uemf(k-1) ) 
          END DO

          DO k = 1, kpen
             IF( comsub(k) .GE. 0._r8 ) THEN
                IF( k .EQ. mkx ) THEN
                   thlten_sub = 0._r8
                   qtten_sub  = 0._r8
                   qlten_sub  = 0._r8
                   qiten_sub  = 0._r8
                   nlten_sub  = 0._r8
                   niten_sub  = 0._r8
                ELSE
                   thlten_sub = g * comsub(k) * ( thl0(k+1) - thl0(k) ) / ( p0(k) - p0(k+1) )
                   qtten_sub  = g * comsub(k) * (  qt0(k+1) -  qt0(k) ) / ( p0(k) - p0(k+1) )
                   qlten_sub  = g * comsub(k) * (  ql0(k+1) -  ql0(k) ) / ( p0(k) - p0(k+1) )
                   qiten_sub  = g * comsub(k) * (  qi0(k+1) -  qi0(k) ) / ( p0(k) - p0(k+1) )
                   nlten_sub  = g * comsub(k) * (  tr0(k+1,ixnumliq) -  tr0(k,ixnumliq) ) / ( p0(k) - p0(k+1) )
                   niten_sub  = g * comsub(k) * (  tr0(k+1,ixnumice) -  tr0(k,ixnumice) ) / ( p0(k) - p0(k+1) )
                ENDIF
             ELSE
                IF( k .EQ. 1 ) THEN
                   thlten_sub = 0._r8
                   qtten_sub  = 0._r8
                   qlten_sub  = 0._r8
                   qiten_sub  = 0._r8
                   nlten_sub  = 0._r8
                   niten_sub  = 0._r8
                ELSE
                   thlten_sub = g * comsub(k) * ( thl0(k) - thl0(k-1) ) / ( p0(k-1) - p0(k) )
                   qtten_sub  = g * comsub(k) * (  qt0(k) -  qt0(k-1) ) / ( p0(k-1) - p0(k) )
                   qlten_sub  = g * comsub(k) * (  ql0(k) -  ql0(k-1) ) / ( p0(k-1) - p0(k) )
                   qiten_sub  = g * comsub(k) * (  qi0(k) -  qi0(k-1) ) / ( p0(k-1) - p0(k) )
                   nlten_sub  = g * comsub(k) * (  tr0(k,ixnumliq) -  tr0(k-1,ixnumliq) ) / ( p0(k-1) - p0(k) )
                   niten_sub  = g * comsub(k) * (  tr0(k,ixnumice) -  tr0(k-1,ixnumice) ) / ( p0(k-1) - p0(k) )
                ENDIF
             ENDIF
             thl_prog = thl0(k) + thlten_sub * dt
             qt_prog  = MAX( qt0(k) + qtten_sub * dt, 1.e-12_r8 )
             CALL conden(p0(k),thl_prog,qt_prog,thj,qvj,qlj,qij,qse,id_check,qsat)
             IF( id_check .EQ. 1 ) THEN
                id_exit = .TRUE.
                go to 333
             ENDIF
             ! qlten_sink(k) = ( qlj - ql0(k) ) / dt
             ! qiten_sink(k) = ( qij - qi0(k) ) / dt
             qlten_sink(k) = MAX( qlten_sub, - ql0(k) / dt ) ! For consistency with prognostic macrophysics scheme
             qiten_sink(k) = MAX( qiten_sub, - qi0(k) / dt ) ! For consistency with prognostic macrophysics scheme
             nlten_sink(k) = MAX( nlten_sub, - tr0(k,ixnumliq) / dt ) 
             niten_sink(k) = MAX( niten_sub, - tr0(k,ixnumice) / dt )
          END DO

          ! --------------------------------------------- !
          !                                               !
          ! Calculate convective tendencies at each layer ! 
          !                                               !
          ! --------------------------------------------- !

          ! ----------------- !
          ! Momentum tendency !
          ! ----------------- !

          DO k = 1, kpen
             km1 = k - 1 
             uten(k) = ( uflx(km1) - uflx(k) ) * g / dp0(k)
             vten(k) = ( vflx(km1) - vflx(k) ) * g / dp0(k) 
             uf(k)   = u0(k) + uten(k) * dt
             vf(k)   = v0(k) + vten(k) * dt
             ! do m = 1, ncnst
             !    trten(k,m) = ( trflx(km1,m) - trflx(k,m) ) * g / dp0(k)
             !  ! Limit trten(k,m) such that negative value is not developed.
             !  ! This limitation does not conserve grid-mean tracers and future
             !  ! refinement is required for tracer-conserving treatment.
             !    trten(k,m) = max(trten(k,m),-tr0(k,m)/dt)              
             ! enddo
          END DO

          ! ----------------------------------------------------------------- !
          ! Tendencies of thermodynamic variables.                            ! 
          ! This part requires a careful treatment of bulk cloud microphysics.!
          ! Relocations of 'precipitable condensates' either into the surface ! 
          ! or into the tendency of 'krel' layer will be performed just after !
          ! finishing the below 'do-loop'.                                    !        
          ! ----------------------------------------------------------------- !

          rliq    = 0._r8
          rainflx = 0._r8
          snowflx = 0._r8

          DO k = 1, kpen

             km1 = k - 1

             ! ------------------------------------------------------------------------------ !
             ! Compute 'slten', 'qtten', 'qvten', 'qlten', 'qiten', and 'sten'                !
             !                                                                                !
             ! Key assumptions made in this 'cumulus scheme' are :                            !
             ! 1. Cumulus updraft expels condensate into the environment at the top interface !
             !    of each layer. Note that in addition to this expel process ('source' term), !
             !    cumulus updraft can modify layer mean condensate through normal detrainment !
             !    forcing or compensating subsidence.                                         !
             ! 2. Expelled water can be either 'sustaining' or 'precipitating' condensate. By !
             !    definition, 'suataining condensate' will remain in the layer where it was   !
             !    formed, while 'precipitating condensate' will fall across the base of the   !
             !    layer where it was formed.                                                  !
             ! 3. All precipitating condensates are assumed to fall into the release layer or !
             !    ground as soon as it was formed without being evaporated during the falling !
             !    process down to the desinated layer ( either release layer of surface ).    !
             ! ------------------------------------------------------------------------------ !

             ! ------------------------------------------------------------------------- !     
             ! 'dwten(k)','diten(k)' : Production rate of condensate  within the layer k !
             !      [ kg/kg/s ]        by the expels of condensate from cumulus updraft. !
             ! It is important to note that in terms of moisture tendency equation, this !
             ! is a 'source' term of enviromental 'qt'.  More importantly,  these source !
             ! are already counted in the turbulent heat and moisture fluxes we computed !
             ! until now, assuming all the expelled condensate remain in the layer where ! 
             ! it was formed. Thus, in calculation of 'qtten' and 'slten' below, we MUST !
             ! NOT add or subtract these terms explicitly in order not to double or miss !
             ! count, unless some expelled condensates fall down out of the layer.  Note !
             ! this falling-down process ( i.e., precipitation process ) and  associated !
             ! 'qtten' and 'slten' and production of surface precipitation flux  will be !
             ! treated later in 'zm_conv_evap' in 'convect_shallow_tend' subroutine.     ! 
             ! In below, we are converting expelled cloud condensate into correct unit.  !
             ! I found that below use of '0.5 * (umf(k-1) + umf(k))' causes conservation !
             ! errors at some columns in global simulation. So, I returned to originals. !
             ! This will cause no precipitation flux at 'kpen' layer since umf(kpen)=0.  !
             ! ------------------------------------------------------------------------- !

             dwten(k) = dwten(k) * 0.5_r8 * ( umf(k-1) + umf(k) ) * g / dp0(k) ! [ kg/kg/s ]
             diten(k) = diten(k) * 0.5_r8 * ( umf(k-1) + umf(k) ) * g / dp0(k) ! [ kg/kg/s ]  

             ! dwten(k) = dwten(k) * umf(k) * g / dp0(k) ! [ kg/kg/s ]
             ! diten(k) = diten(k) * umf(k) * g / dp0(k) ! [ kg/kg/s ]

             ! --------------------------------------------------------------------------- !
             ! 'qrten(k)','qsten(k)' : Production rate of rain and snow within the layer k !
             !     [ kg/kg/s ]         by cumulus expels of condensates to the environment.!         
             ! This will be falled-out of the layer where it was formed and will be dumped !
             ! dumped into the release layer assuming that there is no evaporative cooling !
             ! while precipitable condensate moves to the relaes level. This is reasonable ! 
             ! assumtion if cumulus is purely vertical and so the path along which precita !
             ! ble condensate falls is fully saturared. This 're-allocation' process of    !
             ! precipitable condensate into the release layer is fully described in this   !
             ! convection scheme. After that, the dumped water into the release layer will !
             ! falling down across the base of release layer ( or LCL, if  exact treatment ! 
             ! is required ) and will be allowed to be evaporated in layers below  release !
             ! layer, and finally non-zero surface precipitation flux will be calculated.  !
             ! This latter process will be separately treated 'zm_conv_evap' routine.      !
             ! --------------------------------------------------------------------------- !

             qrten(k) = frc_rasn * dwten(k)
             qsten(k) = frc_rasn * diten(k) 

             ! ----------------------------------------------------------------------- !         
             ! 'rainflx','snowflx' : Cumulative rain and snow flux integrated from the ! 
             !     [ kg/m2/s ]       release leyer to the 'kpen' layer. Note that even !
             ! though wtw(kpen) < 0 (and umf(kpen) = 0) at the top interface of 'kpen' !
             ! layer, 'dwten(kpen)' and diten(kpen)  were calculated after calculating !
             ! explicit cloud top height. Thus below calculation of precipitation flux !
             ! is correct. Note that  precipitating condensates are formed only in the !
             ! layers from 'krel' to 'kpen', including the two layers.                 !
             ! ----------------------------------------------------------------------- !

             rainflx = rainflx + qrten(k) * dp0(k) / g
             snowflx = snowflx + qsten(k) * dp0(k) / g

             ! ------------------------------------------------------------------------ !
             ! 'slten(k)','qtten(k)'                                                    !
             !  Note that 'slflx(k)' and 'qtflx(k)' we have calculated already included !
             !  all the contributions of (1) expels of condensate (dwten(k), diten(k)), !
             !  (2) mass detrainment ( delta * umf * ( qtu - qt ) ), & (3) compensating !
             !  subsidence ( M * dqt / dz ). Thus 'slflx(k)' and 'qtflx(k)' we computed ! 
             !  is a hybrid turbulent flux containing one part of 'source' term - expel !
             !  of condensate. In order to calculate 'slten' and 'qtten', we should add !
             !  additional 'source' term, if any. If the expelled condensate falls down !
             !  across the base of the layer, it will be another sink (negative source) !
             !  term.  Note also that we included frictional heating terms in the below !
             !  calculation of 'slten'.                                                 !
             ! ------------------------------------------------------------------------ !

             slten(k) = ( slflx(km1) - slflx(k) ) * g / dp0(k)
             IF( k .EQ. 1 ) THEN
                slten(k) = slten(k) - g / 4._r8 / dp0(k) * (                            &
                     uflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +     & 
                     vflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)))
             ELSEIF( k .GE. 2 .AND. k .LE. kpen-1 ) THEN
                slten(k) = slten(k) - g / 4._r8 / dp0(k) * (                            &
                     uflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +     &
                     uflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +   &
                     vflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)) +     &
                     vflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
             ELSEIF( k .EQ. kpen ) THEN
                slten(k) = slten(k) - g / 4._r8 / dp0(k) * (                            &
                     uflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +   &
                     vflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
             ENDIF
             qtten(k) = ( qtflx(km1) - qtflx(k) ) * g / dp0(k)

             ! ---------------------------------------------------------------------------- !
             ! Compute condensate tendency, including reserved condensate                   !
             ! We assume that eventual detachment and detrainment occurs in kbup layer  due !
             ! to downdraft buoyancy sorting. In the layer above the kbup, only penetrative !
             ! entrainment exists. Penetrative entrained air is assumed not to contain any  !
             ! condensate.                                                                  !
             ! ---------------------------------------------------------------------------- !

             ! Compute in-cumulus condensate at the layer mid-point.

             IF( k .LT. krel .OR. k .GT. kpen ) THEN
                qlu_mid = 0._r8
                qiu_mid = 0._r8
                qlj     = 0._r8
                qij     = 0._r8
             ELSEIF( k .EQ. krel ) THEN 
                CALL conden(prel,thlu(krel-1),qtu(krel-1),thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                ENDIF
                qlubelow = qlj       
                qiubelow = qij       
                CALL conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                qlu_mid = 0.5_r8 * ( qlubelow + qlj ) * ( prel - ps0(k) )/( ps0(k-1) - ps0(k) )
                qiu_mid = 0.5_r8 * ( qiubelow + qij ) * ( prel - ps0(k) )/( ps0(k-1) - ps0(k) )
             ELSEIF( k .EQ. kpen ) THEN 
                CALL conden(ps0(k-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                qlu_mid = 0.5_r8 * ( qlubelow + qlj ) * ( -ppen )        /( ps0(k-1) - ps0(k) )
                qiu_mid = 0.5_r8 * ( qiubelow + qij ) * ( -ppen )        /( ps0(k-1) - ps0(k) )
                qlu_top = qlj
                qiu_top = qij
             ELSE
                CALL conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                qlu_mid = 0.5_r8 * ( qlubelow + qlj )
                qiu_mid = 0.5_r8 * ( qiubelow + qij )
             ENDIF
             qlubelow = qlj       
             qiubelow = qij       

             ! 1. Sustained Precipitation

             qc_l(k) = ( 1._r8 - frc_rasn ) * dwten(k) ! [ kg/kg/s ]
             qc_i(k) = ( 1._r8 - frc_rasn ) * diten(k) ! [ kg/kg/s ]

             ! 2. Detrained Condensate

             IF( k .LE. kbup ) THEN 
                qc_l(k) = qc_l(k) + g * 0.5_r8 * ( umf(k-1) + umf(k) ) * fdr(k) * qlu_mid ! [ kg/kg/s ]
                qc_i(k) = qc_i(k) + g * 0.5_r8 * ( umf(k-1) + umf(k) ) * fdr(k) * qiu_mid ! [ kg/kg/s ]
                qc_lm   =         - g * 0.5_r8 * ( umf(k-1) + umf(k) ) * fdr(k) * ql0(k)  
                qc_im   =         - g * 0.5_r8 * ( umf(k-1) + umf(k) ) * fdr(k) * qi0(k)
                ! Below 'nc_lm', 'nc_im' should be used only when frc_rasn = 1.
                nc_lm   =         - g * 0.5_r8 * ( umf(k-1) + umf(k) ) * fdr(k) * tr0(k,ixnumliq)  
                nc_im   =         - g * 0.5_r8 * ( umf(k-1) + umf(k) ) * fdr(k) * tr0(k,ixnumice)
             ELSE
                qc_lm   = 0._r8
                qc_im   = 0._r8
                nc_lm   = 0._r8
                nc_im   = 0._r8
             ENDIF

             ! 3. Detached Updraft 

             IF( k .EQ. kbup ) THEN
                qc_l(k) = qc_l(k) + g * umf(k) * qlj     / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                qc_i(k) = qc_i(k) + g * umf(k) * qij     / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                qc_lm   = qc_lm   - g * umf(k) * ql0(k)  / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                qc_im   = qc_im   - g * umf(k) * qi0(k)  / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                nc_lm   = nc_lm   - g * umf(k) * tr0(k,ixnumliq)  / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                nc_im   = nc_im   - g * umf(k) * tr0(k,ixnumice)  / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
             ENDIF

             ! 4. Cumulative Penetrative entrainment detrained in the 'kbup' layer
             !    Explicitly compute the properties detrained penetrative entrained airs in k = kbup layer.

             IF( k .EQ. kbup ) THEN
                CALL conden(p0(k),thlu_emf(k),qtu_emf(k),thj,qvj,ql_emf_kbup,qi_emf_kbup,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   id_exit = .TRUE.
                   go to 333
                ENDIF
                IF( ql_emf_kbup .GT. 0._r8 ) THEN
                   nl_emf_kbup = tru_emf(k,ixnumliq)
                ELSE
                   nl_emf_kbup = 0._r8
                ENDIF
                IF( qi_emf_kbup .GT. 0._r8 ) THEN
                   ni_emf_kbup = tru_emf(k,ixnumice)
                ELSE
                   ni_emf_kbup = 0._r8
                ENDIF
                qc_lm   = qc_lm   - g * emf(k) * ( ql_emf_kbup - ql0(k) ) / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                qc_im   = qc_im   - g * emf(k) * ( qi_emf_kbup - qi0(k) ) / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                nc_lm   = nc_lm   - g * emf(k) * ( nl_emf_kbup - tr0(k,ixnumliq) ) / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
                nc_im   = nc_im   - g * emf(k) * ( ni_emf_kbup - tr0(k,ixnumice) ) / ( ps0(k-1) - ps0(k) ) ! [ kg/kg/s ]
             ENDIF

             qlten_det   = qc_l(k) + qc_lm
             qiten_det   = qc_i(k) + qc_im

             ! --------------------------------------------------------------------------------- !
             ! 'qlten(k)','qiten(k)','qvten(k)','sten(k)'                                        !
             ! Note that falling of precipitation will be treated later.                         !
             ! The prevension of negative 'qv,ql,qi' will be treated later in positive_moisture. !
             ! --------------------------------------------------------------------------------- ! 

             IF( use_expconten ) THEN
                IF( use_unicondet ) THEN
                   qc_l(k) = 0._r8
                   qc_i(k) = 0._r8 
                   qlten(k) = frc_rasn * dwten(k) + qlten_sink(k) + qlten_det
                   qiten(k) = frc_rasn * diten(k) + qiten_sink(k) + qiten_det
                ELSE 
                   qlten(k) = qc_l(k) + frc_rasn * dwten(k) + ( MAX( 0._r8, ql0(k) + ( qc_lm + qlten_sink(k) ) * &
                   dt ) - ql0(k) ) / dt
                   qiten(k) = qc_i(k) + frc_rasn * diten(k) + ( MAX( 0._r8, qi0(k) + ( qc_im + qiten_sink(k) ) * &
                   dt ) - qi0(k) ) / dt
                   trten(k,ixnumliq) = MAX( nc_lm + nlten_sink(k), - tr0(k,ixnumliq) / dt )
                   trten(k,ixnumice) = MAX( nc_im + niten_sink(k), - tr0(k,ixnumice) / dt )
                ENDIF
             ELSE
                IF( use_unicondet ) THEN
                   qc_l(k) = 0._r8
                   qc_i(k) = 0._r8 
                ENDIF
                qlten(k) = dwten(k) + ( qtten(k) - dwten(k) - diten(k) ) * ( ql0(k) / qt0(k) )
                qiten(k) = diten(k) + ( qtten(k) - dwten(k) - diten(k) ) * ( qi0(k) / qt0(k) )
             ENDIF

             qvten(k) = qtten(k) - qlten(k) - qiten(k)
             sten(k)  = slten(k) + xlv * qlten(k) + xls * qiten(k)

             ! -------------------------------------------------------------------------- !
             ! 'rliq' : Verticall-integrated 'suspended cloud condensate'                 !
             !  [m/s]   This is so called 'reserved liquid water'  in other subroutines   ! 
             ! of CAM3, since the contribution of this term should not be included into   !
             ! the tendency of each layer or surface flux (precip)  within this cumulus   !
             ! scheme. The adding of this term to the layer tendency will be done inthe   !
             ! 'stratiform_tend', just after performing sediment process there.           !
             ! The main problem of these rather going-back-and-forth and stupid-seeming   ! 
             ! approach is that the sediment process of suspendened condensate will not   !
             ! be treated at all in the 'stratiform_tend'.                                !
             ! Note that 'precip' [m/s] is vertically-integrated total 'rain+snow' formed !
             ! from the cumulus updraft. Important : in the below, 1000 is rhoh2o ( water !
             ! density ) [ kg/m^3 ] used for unit conversion from [ kg/m^2/s ] to [ m/s ] !
             ! for use in stratiform.F90.                                                 !
             ! -------------------------------------------------------------------------- ! 

             qc(k)  =  qc_l(k) +  qc_i(k)   
             rliq   =  rliq    + qc(k) * dp0(k) / g / 1000._r8    ! [ m/s ]

          END DO

          precip  =  rainflx + snowflx                       ! [ kg/m2/s ]
          snow    =  snowflx                                 ! [ kg/m2/s ] 

          ! ---------------------------------------------------------------- !
          ! Now treats the 'evaporation' and 'melting' of rain ( qrten ) and ! 
          ! snow ( qsten ) during falling process. Below algorithms are from !
          ! 'zm_conv_evap' but with some modification, which allows separate !
          ! treatment of 'rain' and 'snow' condensates. Note that I included !
          ! the evaporation dynamics into the convection scheme for complete !
          ! development of cumulus scheme especially in association with the ! 
          ! implicit CIN closure. In compatible with this internal treatment !
          ! of evaporation, I should modify 'convect_shallow',  in such that !
          ! 'zm_conv_evap' is not performed when I choose UW PBL-Cu schemes. !                                          
          ! ---------------------------------------------------------------- !

          evpint_rain    = 0._r8 
          evpint_snow    = 0._r8
          flxrain(0:mkx) = 0._r8
          flxsnow(0:mkx) = 0._r8
          ntraprd(:mkx)  = 0._r8
          ntsnprd(:mkx)  = 0._r8

          DO k = mkx, 1, -1  ! 'k' is a layer index : 'mkx'('1') is the top ('bottom') layer

             ! ----------------------------------------------------------------------------- !
             ! flxsntm [kg/m2/s] : Downward snow flux at the top of each layer after melting.! 
             ! snowmlt [kg/kg/s] : Snow melting tendency.                                    !
             ! Below allows melting of snow when it goes down into the warm layer below.     !
             ! ----------------------------------------------------------------------------- !

             IF( t0(k) .GT. 273.16_r8 ) THEN
                snowmlt = MAX( 0._r8, flxsnow(k) * g / dp0(k) ) 
             ELSE
                snowmlt = 0._r8
             ENDIF

             ! ----------------------------------------------------------------- !
             ! Evaporation rate of 'rain' and 'snow' in the layer k, [ kg/kg/s ] !
             ! where 'rain' and 'snow' are coming down from the upper layers.    !
             ! I used the same evaporative efficiency both for 'rain' and 'snow'.!
             ! Note that evaporation is not allowed in the layers 'k >= krel' by !
             ! assuming that inside of cumulus cloud, across which precipitation !
             ! is falling down, is fully saturated.                              !
             ! The asumptions in association with the 'evplimit_rain(snow)' are  !
             !   1. Do not allow evaporation to supersate the layer              !
             !   2. Do not evaporate more than the flux falling into the layer   !
             !   3. Total evaporation cannot exceed the input total surface flux !
             ! ----------------------------------------------------------------- !

             status = qsat(t0(k),p0(k),es(1),qs(1),gam(1), 1)          
             subsat = MAX( ( 1._r8 - qv0(k)/qs(1) ), 0._r8 )
             IF( noevap_krelkpen ) THEN
                IF( k .GE. krel ) subsat = 0._r8
             ENDIF

             evprain  = kevp * subsat * SQRT(flxrain(k)+snowmlt*dp0(k)/g) 
             evpsnow  = kevp * subsat * SQRT(MAX(flxsnow(k)-snowmlt*dp0(k)/g,0._r8))

             evplimit = MAX( 0._r8, ( qw0_in(i,k) - qv0(k) ) / dt ) 

             evplimit_rain = MIN( evplimit,      ( flxrain(k) + snowmlt * dp0(k) / g ) * g / dp0(k) )
             evplimit_rain = MIN( evplimit_rain, ( rainflx - evpint_rain ) * g / dp0(k) )
             evprain = MAX(0._r8,MIN( evplimit_rain, evprain ))

             evplimit_snow = MIN( evplimit,   MAX( flxsnow(k) - snowmlt * dp0(k) / g , 0._r8 ) * g / dp0(k) )
             evplimit_snow = MIN( evplimit_snow, ( snowflx - evpint_snow ) * g / dp0(k) )
             evpsnow = MAX(0._r8,MIN( evplimit_snow, evpsnow ))

             IF( ( evprain + evpsnow ) .GT. evplimit ) THEN
                tmp1 = evprain * evplimit / ( evprain + evpsnow )
                tmp2 = evpsnow * evplimit / ( evprain + evpsnow )
                evprain = tmp1
                evpsnow = tmp2
             ENDIF

             evapc(k) = evprain + evpsnow

             ! ------------------------------------------------------------- !
             ! Vertically-integrated evaporative fluxes of 'rain' and 'snow' !
             ! ------------------------------------------------------------- !

             evpint_rain = evpint_rain + evprain * dp0(k) / g
             evpint_snow = evpint_snow + evpsnow * dp0(k) / g

             ! -------------------------------------------------------------- !
             ! Net 'rain' and 'snow' production rate in the layer [ kg/kg/s ] !
             ! -------------------------------------------------------------- !         

             ntraprd(k) = qrten(k) - evprain + snowmlt
             ntsnprd(k) = qsten(k) - evpsnow - snowmlt

             ! -------------------------------------------------------------------------------- !
             ! Downward fluxes of 'rain' and 'snow' fluxes at the base of the layer [ kg/m2/s ] !
             ! Note that layer index increases with height.                                     !
             ! -------------------------------------------------------------------------------- !

             flxrain(k-1) = flxrain(k) + ntraprd(k) * dp0(k) / g
             flxsnow(k-1) = flxsnow(k) + ntsnprd(k) * dp0(k) / g
             flxrain(k-1) = MAX( flxrain(k-1), 0._r8 )
             IF( flxrain(k-1) .EQ. 0._r8 ) ntraprd(k) = -flxrain(k) * g / dp0(k)
             flxsnow(k-1) = MAX( flxsnow(k-1), 0._r8 )         
             IF( flxsnow(k-1) .EQ. 0._r8 ) ntsnprd(k) = -flxsnow(k) * g / dp0(k)

             ! ---------------------------------- !
             ! Calculate thermodynamic tendencies !
             ! --------------------------------------------------------------------------- !
             ! Note that equivalently, we can write tendency formula of 'sten' and 'slten' !
             ! by 'sten(k)  = sten(k) - xlv*evprain  - xls*evpsnow - (xls-xlv)*snowmlt' &  !
             !    'slten(k) = sten(k) - xlv*qlten(k) - xls*qiten(k)'.                      !
             ! The above formula is equivalent to the below formula. However below formula !
             ! is preferred since we have already imposed explicit constraint on 'ntraprd' !
             ! and 'ntsnprd' in case that flxrain(k-1) < 0 & flxsnow(k-1) < 0._r8          !
             ! Note : In future, I can elborate the limiting of 'qlten','qvten','qiten'    !
             !        such that that energy and moisture conservation error is completely  !
             !        suppressed.                                                          !
             ! Re-storation to the positive condensate will be performed later below       !
             ! --------------------------------------------------------------------------- !

             qlten(k) = qlten(k) - qrten(k)
             qiten(k) = qiten(k) - qsten(k)
             qvten(k) = qvten(k) + evprain  + evpsnow
             qtten(k) = qlten(k) + qiten(k) + qvten(k)
             IF( ( qv0(k) + qvten(k)*dt ) .LT. qmin(1) .OR. &
                  ( ql0(k) + qlten(k)*dt ) .LT. qmin(2) .OR. &
                  ( qi0(k) + qiten(k)*dt ) .LT. qmin(3) ) THEN
                limit_negcon(i) = 1._r8
             END IF
             sten(k)  = sten(k) - xlv*evprain  - xls*evpsnow - (xls-xlv)*snowmlt
             slten(k) = sten(k) - xlv*qlten(k) - xls*qiten(k)

             !  slten(k) = slten(k) + xlv * ntraprd(k) + xls * ntsnprd(k)         
             !  sten(k)  = slten(k) + xlv * qlten(k)   + xls * qiten(k)

          END DO

          ! ------------------------------------------------------------- !
          ! Calculate final surface flux of precipitation, rain, and snow !
          ! Convert unit to [m/s] for use in 'check_energy_chng'.         !  
          ! ------------------------------------------------------------- !

          precip  = ( flxrain(0) + flxsnow(0) ) / 1000._r8
          snow    =   flxsnow(0) / 1000._r8       

          ! --------------------------------------------------------------------------- !
          ! Until now, all the calculations are done completely in this shallow cumulus !
          ! scheme. If you want to use this cumulus scheme other than CAM3, then do not !
          ! perform below block. However, for compatible use with the other subroutines !
          ! in CAM3, I should subtract the effect of 'qc(k)' ('rliq') from the tendency !
          ! equation in each layer, since this effect will be separately added later in !
          ! in 'stratiform_tend' just after performing sediment process there. In order !
          ! to be consistent with 'stratiform_tend', just subtract qc(k)  from tendency !
          ! equation of each layer, but do not add it to the 'precip'. Apprently,  this !
          ! will violate energy and moisture conservations.    However, when performing !
          ! conservation check in 'tphysbc.F90' just after 'convect_shallow_tend',   we !
          ! will add 'qc(k)' ( rliq ) to the surface flux term just for the purpose  of !
          ! passing the energy-moisture conservation check. Explicit adding-back of 'qc'!
          ! to the individual layer tendency equation will be done in 'stratiform_tend' !
          ! after performing sediment process there. Simply speaking, in 'tphysbc' just !
          ! after 'convect_shallow_tend', we will dump 'rliq' into surface as a  'rain' !
          ! in order to satisfy energy and moisture conservation, and  in the following !
          ! 'stratiform_tend', we will restore it back to 'qlten(k)' ( 'ice' will go to !  
          ! 'water' there) from surface precipitation. This is a funny but conceptually !
          ! entertaining procedure. One concern I have for this complex process is that !
          ! output-writed stratiform precipitation amount will be underestimated due to !
          ! arbitrary subtracting of 'rliq' in stratiform_tend, where                   !
          ! ' prec_str = prec_sed + prec_pcw - rliq' and 'rliq' is not real but fake.   ! 
          ! However, as shown in 'srfxfer.F90', large scale precipitation amount (PRECL)!
          ! that is writed-output is corrected written since in 'srfxfer.F90',  PRECL = !
          ! 'prec_sed + prec_pcw', without including 'rliq'. So current code is correct.!
          ! Note also in 'srfxfer.F90', convective precipitation amount is 'PRECC =     ! 
          ! prec_zmc(i) + prec_cmf(i)' which is also correct.                           !
          ! --------------------------------------------------------------------------- !

          DO k = 1, kpen       
             qtten(k) = qtten(k) - qc(k)
             qlten(k) = qlten(k) - qc_l(k)
             qiten(k) = qiten(k) - qc_i(k)
             slten(k) = slten(k) + ( xlv * qc_l(k) + xls * qc_i(k) )
             ! ---------------------------------------------------------------------- !
             ! Since all reserved condensates will be treated  as liquid water in the !
             ! 'check_energy_chng' & 'stratiform_tend' without an explicit conversion !
             ! algorithm, I should consider explicitly the energy conversions between !
             ! 'ice' and 'liquid' - i.e., I should convert 'ice' to 'liquid'  and the !
             ! necessary energy for this conversion should be subtracted from 'sten'. ! 
             ! Without this conversion here, energy conservation error come out. Note !
             ! that there should be no change of 'qvten(k)'.                          !
             ! ---------------------------------------------------------------------- !
             sten(k)  = sten(k)  - ( xls - xlv ) * qc_i(k)
          END DO

          ! --------------------------------------------------------------- !
          ! Prevent the onset-of negative condensate at the next time step  !
          ! Potentially, this block can be moved just in front of the above !
          ! block.                                                          ! 
          ! --------------------------------------------------------------- !

          ! Modification : I should check whether this 'positive_moisture_single' routine is
          !                consistent with the one used in UW PBL and cloud macrophysics schemes.
          ! Modification : Below may overestimate resulting 'ql, qi' if we use the new 'qc_l', 'qc_i'
          !                in combination with the original computation of qlten, qiten. However,
          !                if we use new 'qlten,qiten', there is no problem.

          qv0_star(:mkx) = qv0(:mkx) + qvten(:mkx) * dt
          ql0_star(:mkx) = ql0(:mkx) + qlten(:mkx) * dt
          qi0_star(:mkx) = qi0(:mkx) + qiten(:mkx) * dt
          s0_star(:mkx)  =  s0(:mkx) +  sten(:mkx) * dt
          CALL positive_moisture_single( xlv, xls, mkx, dt, qmin(1), qmin(2), qmin(3), dp0, qv0_star,&
           ql0_star, qi0_star, s0_star, qvten, qlten, qiten, sten )
          qtten(:mkx)    = qvten(:mkx) + qlten(:mkx) + qiten(:mkx)
          slten(:mkx)    = sten(:mkx)  - xlv * qlten(:mkx) - xls * qiten(:mkx)

          ! --------------------- !
          ! Tendencies of tracers !
          ! --------------------- !

          DO m = 4, ncnst

             IF( m .NE. ixnumliq .AND. m .NE. ixnumice ) THEN

                trmin = qmin(m)
                !#ifdef MODAL_AERO
                !          do mm = 1, ntot_amode
                !             if( m .eq. numptr_amode(mm) ) then
                !                 trmin = 1.e-5_r8
                !                 goto 55
                !             endif              
                !          enddo
                !       55 continue
                !#endif 
                trflx_d(0:mkx) = 0._r8
                trflx_u(0:mkx) = 0._r8           
                DO k = 1, mkx-1
                   IF( cnst_get_type_byind(m) .EQ. 'wet' ) THEN
                      pdelx = dp0(k)
                   ELSE
                      pdelx = dpdry0(k)
                   ENDIF
                   km1 = k - 1
                   dum = ( tr0(k,m) - trmin ) *  pdelx / g / dt + trflx(km1,m) - trflx(k,m) + trflx_d(km1)
                   trflx_d(k) = MIN( 0._r8, dum )
                ENDDO
                DO k = mkx, 2, -1
                   IF( cnst_get_type_byind(m) .EQ. 'wet' ) THEN
                      pdelx = dp0(k)
                   ELSE
                      pdelx = dpdry0(k)
                   ENDIF
                   km1 = k - 1
                   dum = ( tr0(k,m) - trmin ) * pdelx / g / dt + trflx(km1,m) - trflx(k,m) + &
                        trflx_d(km1) - trflx_d(k) - trflx_u(k) 
                   trflx_u(km1) = MAX( 0._r8, -dum ) 
                ENDDO
                DO k = 1, mkx
                   IF( cnst_get_type_byind(m) .EQ. 'wet' ) THEN
                      pdelx = dp0(k)
                   ELSE
                      pdelx = dpdry0(k)
                   ENDIF
                   km1 = k - 1
                   ! Check : I should re-check whether '_u', '_d' are correctly ordered in 
                   !         the below tendency computation.
                   trten(k,m) = ( trflx(km1,m) - trflx(k,m) + & 
                        trflx_d(km1) - trflx_d(k) + &
                        trflx_u(km1) - trflx_u(k) ) * g / pdelx
                ENDDO

             ENDIF

          ENDDO

          ! ---------------------------------------------------------------- !
          ! Cumpute default diagnostic outputs                               !
          ! Note that since 'qtu(krel-1:kpen-1)' & 'thlu(krel-1:kpen-1)' has !
          ! been adjusted after detraining cloud condensate into environment ! 
          ! during cumulus updraft motion,  below calculations will  exactly !
          ! reproduce in-cloud properties as shown in the output analysis.   !
          ! ---------------------------------------------------------------- ! 

          CALL conden(prel,thlu(krel-1),qtu(krel-1),thj,qvj,qlj,qij,qse,id_check,qsat)
          IF( id_check .EQ. 1 ) THEN
             exit_conden(i) = 1._r8
             id_exit = .TRUE.
             go to 333
          END IF
          qcubelow = qlj + qij
          qlubelow = qlj       
          qiubelow = qij       
          rcwp     = 0._r8
          rlwp     = 0._r8
          riwp     = 0._r8

          ! --------------------------------------------------------------------- !
          ! In the below calculations, I explicitly considered cloud base ( LCL ) !
          ! and cloud top height ( ps0(kpen-1) + ppen )                           !
          ! ----------------------------------------------------------------------! 
          DO k = krel, kpen ! This is a layer index
             ! ------------------------------------------------------------------ ! 
             ! Calculate cumulus condensate at the upper interface of each layer. !
             ! Note 'ppen < 0' and at 'k=kpen' layer, I used 'thlu_top'&'qtu_top' !
             ! which explicitly considered zero or non-zero 'fer(kpen)'.          !
             ! ------------------------------------------------------------------ ! 
             IF( k .EQ. kpen ) THEN 
                CALL conden(ps0(k-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check,qsat)
             ELSE
                CALL conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
             ENDIF
             IF( id_check .EQ. 1 ) THEN
                exit_conden(i) = 1._r8
                id_exit = .TRUE.
                go to 333
             END IF
             ! ---------------------------------------------------------------- !
             ! Calculate in-cloud mean LWC ( qlu(k) ), IWC ( qiu(k) ),  & layer !
             ! mean cumulus fraction ( cufrc(k) ),  vertically-integrated layer !
             ! mean LWP and IWP. Expel some of in-cloud condensate at the upper !
             ! interface if it is largr than criqc. Note cumulus cloud fraction !
             ! is assumed to be twice of core updraft fractional area. Thus LWP !
             ! and IWP will be twice of actual value coming from our scheme.    !
             ! ---------------------------------------------------------------- !
             qcu(k)   = 0.5_r8 * ( qcubelow + qlj + qij )
             qlu(k)   = 0.5_r8 * ( qlubelow + qlj )
             qiu(k)   = 0.5_r8 * ( qiubelow + qij )
             cufrc(k) = ( ufrc(k-1) + ufrc(k) )
             IF( k .EQ. krel ) THEN
                cufrc(k) = ( ufrclcl + ufrc(k) )*( prel - ps0(k) )/( ps0(k-1) - ps0(k) )
             ELSE IF( k .EQ. kpen ) THEN
                cufrc(k) = ( ufrc(k-1) + 0._r8 )*( -ppen )        /( ps0(k-1) - ps0(k) )
                IF( (qlj + qij) .GT. criqc ) THEN           
                   qcu(k) = 0.5_r8 * ( qcubelow + criqc )
                   qlu(k) = 0.5_r8 * ( qlubelow + criqc * qlj / ( qlj + qij ) )
                   qiu(k) = 0.5_r8 * ( qiubelow + criqc * qij / ( qlj + qij ) )
                ENDIF
             ENDIF
             rcwp = rcwp + ( qlu(k) + qiu(k) ) * ( ps0(k-1) - ps0(k) ) / g * cufrc(k)
             rlwp = rlwp +   qlu(k)            * ( ps0(k-1) - ps0(k) ) / g * cufrc(k)
             riwp = riwp +   qiu(k)            * ( ps0(k-1) - ps0(k) ) / g * cufrc(k)
             qcubelow = qlj + qij
             qlubelow = qlj
             qiubelow = qij
          END DO
          ! ------------------------------------ !      
          ! Cloud top and base interface indices !
          ! ------------------------------------ !
          cnt = REAL( kpen, r8 )
          cnb = REAL( krel - 1, r8 )

          ! ------------------------------------------------------------------------- !
          ! End of formal calculation. Below blocks are for implicit CIN calculations ! 
          ! with re-initialization and save variables at iter_cin = 1._r8             !
          ! ------------------------------------------------------------------------- !

          ! --------------------------------------------------------------- !
          ! Adjust the original input profiles for implicit CIN calculation !
          ! --------------------------------------------------------------- !

          IF( iter .NE. iter_cin ) THEN 

             ! ------------------------------------------------------------------- !
             ! Save the output from "iter_cin = 1"                                 !
             ! These output will be writed-out if "iter_cin = 1" was not performed !
             ! for some reasons.                                                   !
             ! ------------------------------------------------------------------- !

             qv0_s(:mkx)           = qv0(:mkx) + qvten(:mkx) * dt
             ql0_s(:mkx)           = ql0(:mkx) + qlten(:mkx) * dt
             qi0_s(:mkx)           = qi0(:mkx) + qiten(:mkx) * dt
             s0_s(:mkx)            = s0(:mkx)  +  sten(:mkx) * dt 
             u0_s(:mkx)            = u0(:mkx)  +  uten(:mkx) * dt
             v0_s(:mkx)            = v0(:mkx)  +  vten(:mkx) * dt 
             qt0_s(:mkx)           = qv0_s(:mkx) + ql0_s(:mkx) + qi0_s(:mkx)
             t0_s(:mkx)            = t0(:mkx)  +  sten(:mkx) * dt / cp
             DO m = 1, ncnst
                tr0_s(:mkx,m)      = tr0(:mkx,m) + trten(:mkx,m) * dt
             ENDDO

             umf_s(0:mkx)          = umf(0:mkx)
             qvten_s(:mkx)         = qvten(:mkx)
             qlten_s(:mkx)         = qlten(:mkx)  
             qiten_s(:mkx)         = qiten(:mkx)
             sten_s(:mkx)          = sten(:mkx)
             uten_s(:mkx)          = uten(:mkx)  
             vten_s(:mkx)          = vten(:mkx)
             qrten_s(:mkx)         = qrten(:mkx)
             qsten_s(:mkx)         = qsten(:mkx)  
             precip_s              = precip
             snow_s                = snow
             evapc_s(:mkx)         = evapc(:mkx)
             cush_s                = cush
             cufrc_s(:mkx)         = cufrc(:mkx)  
             slflx_s(0:mkx)        = slflx(0:mkx)  
             qtflx_s(0:mkx)        = qtflx(0:mkx)  
             qcu_s(:mkx)           = qcu(:mkx)  
             qlu_s(:mkx)           = qlu(:mkx)  
             qiu_s(:mkx)           = qiu(:mkx)  
             fer_s(:mkx)           = fer(:mkx)  
             fdr_s(:mkx)           = fdr(:mkx)  
             cin_s                 = cin
             cinlcl_s              = cinlcl
             cbmf_s                = cbmf
             rliq_s                = rliq
             qc_s(:mkx)            = qc(:mkx)
             cnt_s                 = cnt
             cnb_s                 = cnb
             qtten_s(:mkx)         = qtten(:mkx)
             slten_s(:mkx)         = slten(:mkx)
             ufrc_s(0:mkx)         = ufrc(0:mkx) 

             uflx_s(0:mkx)         = uflx(0:mkx)  
             vflx_s(0:mkx)         = vflx(0:mkx)  

             ufrcinvbase_s         = ufrcinvbase
             ufrclcl_s             = ufrclcl 
             winvbase_s            = winvbase
             wlcl_s                = wlcl
             plcl_s                = plcl
             pinv_s                = ps0(kinv-1)
             plfc_s                = plfc        
             pbup_s                = ps0(kbup)
             ppen_s                = ps0(kpen-1) + ppen        
             qtsrc_s               = qtsrc
             thlsrc_s              = thlsrc
             thvlsrc_s             = thvlsrc
             emfkbup_s             = emf(kbup)
             cbmflimit_s           = cbmflimit
             tkeavg_s              = tkeavg
             zinv_s                = zs0(kinv-1)
             rcwp_s                = rcwp
             rlwp_s                = rlwp
             riwp_s                = riwp

             wu_s(0:mkx)           = wu(0:mkx)
             qtu_s(0:mkx)          = qtu(0:mkx)
             thlu_s(0:mkx)         = thlu(0:mkx)
             thvu_s(0:mkx)         = thvu(0:mkx)
             uu_s(0:mkx)           = uu(0:mkx)
             vu_s(0:mkx)           = vu(0:mkx)
             qtu_emf_s(0:mkx)      = qtu_emf(0:mkx)
             thlu_emf_s(0:mkx)     = thlu_emf(0:mkx)
             uu_emf_s(0:mkx)       = uu_emf(0:mkx)
             vu_emf_s(0:mkx)       = vu_emf(0:mkx)
             uemf_s(0:mkx)         = uemf(0:mkx)

             dwten_s(:mkx)         = dwten(:mkx)
             diten_s(:mkx)         = diten(:mkx)
             flxrain_s(0:mkx)      = flxrain(0:mkx)
             flxsnow_s(0:mkx)      = flxsnow(0:mkx)
             ntraprd_s(:mkx)       = ntraprd(:mkx)
             ntsnprd_s(:mkx)       = ntsnprd(:mkx)

             excessu_arr_s(:mkx)   = excessu_arr(:mkx)
             excess0_arr_s(:mkx)   = excess0_arr(:mkx)
             xc_arr_s(:mkx)        = xc_arr(:mkx)
             aquad_arr_s(:mkx)     = aquad_arr(:mkx)
             bquad_arr_s(:mkx)     = bquad_arr(:mkx)
             cquad_arr_s(:mkx)     = cquad_arr(:mkx)
             bogbot_arr_s(:mkx)    = bogbot_arr(:mkx)
             bogtop_arr_s(:mkx)    = bogtop_arr(:mkx)

             DO m = 1, ncnst
                trten_s(:mkx,m)    = trten(:mkx,m)
                trflx_s(0:mkx,m)   = trflx(0:mkx,m)
                tru_s(0:mkx,m)     = tru(0:mkx,m)
                tru_emf_s(0:mkx,m) = tru_emf(0:mkx,m)
             ENDDO

             ! ----------------------------------------------------------------------------- ! 
             ! Recalculate environmental variables for new cin calculation at "iter_cin = 2" ! 
             ! using the updated state variables. Perform only for variables necessary  for  !
             ! the new cin calculation.                                                      !
             ! ----------------------------------------------------------------------------- !

             qv0(:mkx)   = qv0_s(:mkx)
             ql0(:mkx)   = ql0_s(:mkx)
             qi0(:mkx)   = qi0_s(:mkx)
             s0(:mkx)    = s0_s(:mkx)
             t0(:mkx)    = t0_s(:mkx)

             qt0(:mkx)   = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
             thl0(:mkx)  = (t0(:mkx) - xlv*ql0(:mkx)/cp - xls*qi0(:mkx)/cp)/exn0(:mkx)
             thvl0(:mkx) = (1._r8 + zvir*qt0(:mkx))*thl0(:mkx)

             ssthl0      = slope(mkx,thl0,p0) ! Dimension of ssthl0(:mkx) is implicit
             ssqt0       = slope(mkx,qt0 ,p0)
             ssu0        = slope(mkx,u0  ,p0)
             ssv0        = slope(mkx,v0  ,p0)
             DO m = 1, ncnst
                sstr0(:mkx,m) = slope(mkx,tr0(:mkx,m),p0)
             ENDDO

             DO k = 1, mkx

                thl0bot = thl0(k) + ssthl0(k) * ( ps0(k-1) - p0(k) )
                qt0bot  = qt0(k)  + ssqt0(k)  * ( ps0(k-1) - p0(k) )
                CALL conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                thv0bot(k)  = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                thvl0bot(k) = thl0bot * ( 1._r8 + zvir*qt0bot )

                thl0top = thl0(k) + ssthl0(k) * ( ps0(k) - p0(k) )
                qt0top  =  qt0(k) + ssqt0(k)  * ( ps0(k) - p0(k) )
                CALL conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check,qsat)
                IF( id_check .EQ. 1 ) THEN
                   exit_conden(i) = 1._r8
                   id_exit = .TRUE.
                   go to 333
                END IF
                thv0top(k)  = thj * ( 1._r8 + zvir*qvj - qlj - qij )
                thvl0top(k) = thl0top * ( 1._r8 + zvir*qt0top )

             END DO

          ENDIF               ! End of 'if(iter .ne. iter_cin)' if sentence. 

       END DO                ! End of implicit CIN loop (cin_iter)      

       ! ----------------------- !
       ! Update Output Variables !
       ! ----------------------- !

       umf_out(i,0:mkx)             = umf(0:mkx)
       slflx_out(i,0:mkx)           = slflx(0:mkx)
       qtflx_out(i,0:mkx)           = qtflx(0:mkx)
       qvten_out(i,:mkx)            = qvten(:mkx)
       qlten_out(i,:mkx)            = qlten(:mkx)
       qiten_out(i,:mkx)            = qiten(:mkx)
       sten_out(i,:mkx)             = sten(:mkx)
       uten_out(i,:mkx)             = uten(:mkx)
       vten_out(i,:mkx)             = vten(:mkx)
       qrten_out(i,:mkx)            = qrten(:mkx)
       qsten_out(i,:mkx)            = qsten(:mkx)
       precip_out(i)                = precip
       snow_out(i)                  = snow
       evapc_out(i,:mkx)            = evapc(:mkx)
       cufrc_out(i,:mkx)            = cufrc(:mkx)
       qcu_out(i,:mkx)              = qcu(:mkx)
       qlu_out(i,:mkx)              = qlu(:mkx)
       qiu_out(i,:mkx)              = qiu(:mkx)
       cush_inout(i)                = cush
       cbmf_out(i)                  = cbmf
       rliq_out(i)                  = rliq
       qc_out(i,:mkx)               = qc(:mkx)
       cnt_out(i)                   = cnt
       cnb_out(i)                   = cnb

       DO m = 1, ncnst
          trten_out(i,:mkx,m)       = trten(:mkx,m)
       ENDDO

       ! ------------------------------------------------- !
       ! Below are specific diagnostic output for detailed !
       ! analysis of cumulus scheme                        !
       ! ------------------------------------------------- !

       fer_out(i,mkx:1:-1)          = fer(:mkx)  
       fdr_out(i,mkx:1:-1)          = fdr(:mkx)  
       cinh_out(i)                  = cin
       cinlclh_out(i)               = cinlcl
       qtten_out(i,mkx:1:-1)        = qtten(:mkx)
       slten_out(i,mkx:1:-1)        = slten(:mkx)
       ufrc_out(i,mkx:0:-1)         = ufrc(0:mkx)
       uflx_out(i,mkx:0:-1)         = uflx(0:mkx)  
       vflx_out(i,mkx:0:-1)         = vflx(0:mkx)  

       ufrcinvbase_out(i)           = ufrcinvbase
       ufrclcl_out(i)               = ufrclcl 
       winvbase_out(i)              = winvbase
       wlcl_out(i)                  = wlcl
       plcl_out(i)                  = plcl
       pinv_out(i)                  = ps0(kinv-1)
       plfc_out(i)                  = plfc    
       pbup_out(i)                  = ps0(kbup)        
       ppen_out(i)                  = ps0(kpen-1) + ppen            
       qtsrc_out(i)                 = qtsrc
       thlsrc_out(i)                = thlsrc
       thvlsrc_out(i)               = thvlsrc
       emfkbup_out(i)               = emf(kbup)
       cbmflimit_out(i)             = cbmflimit
       tkeavg_out(i)                = tkeavg
       zinv_out(i)                  = zs0(kinv-1)
       rcwp_out(i)                  = rcwp
       rlwp_out(i)                  = rlwp
       riwp_out(i)                  = riwp

       wu_out(i,mkx:0:-1)           = wu(0:mkx)
       qtu_out(i,mkx:0:-1)          = qtu(0:mkx)
       thlu_out(i,mkx:0:-1)         = thlu(0:mkx)
       thvu_out(i,mkx:0:-1)         = thvu(0:mkx)
       uu_out(i,mkx:0:-1)           = uu(0:mkx)
       vu_out(i,mkx:0:-1)           = vu(0:mkx)
       qtu_emf_out(i,mkx:0:-1)      = qtu_emf(0:mkx)
       thlu_emf_out(i,mkx:0:-1)     = thlu_emf(0:mkx)
       uu_emf_out(i,mkx:0:-1)       = uu_emf(0:mkx)
       vu_emf_out(i,mkx:0:-1)       = vu_emf(0:mkx)
       uemf_out(i,mkx:0:-1)         = uemf(0:mkx)

       dwten_out(i,mkx:1:-1)        = dwten(:mkx)
       diten_out(i,mkx:1:-1)        = diten(:mkx)
       flxrain_out(i,mkx:0:-1)      = flxrain(0:mkx)
       flxsnow_out(i,mkx:0:-1)      = flxsnow(0:mkx)
       ntraprd_out(i,mkx:1:-1)      = ntraprd(:mkx)
       ntsnprd_out(i,mkx:1:-1)      = ntsnprd(:mkx)

       excessu_arr_out(i,mkx:1:-1)  = excessu_arr(:mkx)
       excess0_arr_out(i,mkx:1:-1)  = excess0_arr(:mkx)
       xc_arr_out(i,mkx:1:-1)       = xc_arr(:mkx)
       aquad_arr_out(i,mkx:1:-1)    = aquad_arr(:mkx)
       bquad_arr_out(i,mkx:1:-1)    = bquad_arr(:mkx)
       cquad_arr_out(i,mkx:1:-1)    = cquad_arr(:mkx)
       bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr(:mkx)
       bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr(:mkx)

       DO m = 1, ncnst
          trflx_out(i,mkx:0:-1,m)   = trflx(0:mkx,m)  
          tru_out(i,mkx:0:-1,m)     = tru(0:mkx,m)
          tru_emf_out(i,mkx:0:-1,m) = tru_emf(0:mkx,m)
       ENDDO

333    IF(id_exit) THEN ! Exit without cumulus convection

          exit_UWCu(i) = 1._r8

          ! --------------------------------------------------------------------- !
          ! Initialize output variables when cumulus convection was not performed.!
          ! --------------------------------------------------------------------- !

          umf_out(i,0:mkx)             = 0._r8   
          slflx_out(i,0:mkx)           = 0._r8
          qtflx_out(i,0:mkx)           = 0._r8
          qvten_out(i,:mkx)            = 0._r8
          qlten_out(i,:mkx)            = 0._r8
          qiten_out(i,:mkx)            = 0._r8
          sten_out(i,:mkx)             = 0._r8
          uten_out(i,:mkx)             = 0._r8
          vten_out(i,:mkx)             = 0._r8
          qrten_out(i,:mkx)            = 0._r8
          qsten_out(i,:mkx)            = 0._r8
          precip_out(i)                = 0._r8
          snow_out(i)                  = 0._r8
          evapc_out(i,:mkx)            = 0._r8
          cufrc_out(i,:mkx)            = 0._r8
          qcu_out(i,:mkx)              = 0._r8
          qlu_out(i,:mkx)              = 0._r8
          qiu_out(i,:mkx)              = 0._r8
          cush_inout(i)                = -1._r8
          cbmf_out(i)                  = 0._r8   
          rliq_out(i)                  = 0._r8
          qc_out(i,:mkx)               = 0._r8
          cnt_out(i)                   = 1._r8
          cnb_out(i)                   = REAL(mkx, r8)

          fer_out(i,mkx:1:-1)          = 0._r8  
          fdr_out(i,mkx:1:-1)          = 0._r8  
          cinh_out(i)                  = -1._r8 
          cinlclh_out(i)               = -1._r8 
          qtten_out(i,mkx:1:-1)        = 0._r8
          slten_out(i,mkx:1:-1)        = 0._r8
          ufrc_out(i,mkx:0:-1)         = 0._r8
          uflx_out(i,mkx:0:-1)         = 0._r8  
          vflx_out(i,mkx:0:-1)         = 0._r8  

          ufrcinvbase_out(i)           = 0._r8 
          ufrclcl_out(i)               = 0._r8 
          winvbase_out(i)              = 0._r8    
          wlcl_out(i)                  = 0._r8    
          plcl_out(i)                  = 0._r8    
          pinv_out(i)                  = 0._r8     
          plfc_out(i)                  = 0._r8     
          pbup_out(i)                  = 0._r8    
          ppen_out(i)                  = 0._r8    
          qtsrc_out(i)                 = 0._r8    
          thlsrc_out(i)                = 0._r8    
          thvlsrc_out(i)               = 0._r8    
          emfkbup_out(i)               = 0._r8
          cbmflimit_out(i)             = 0._r8    
          tkeavg_out(i)                = 0._r8    
          zinv_out(i)                  = 0._r8    
          rcwp_out(i)                  = 0._r8    
          rlwp_out(i)                  = 0._r8    
          riwp_out(i)                  = 0._r8    

          wu_out(i,mkx:0:-1)           = 0._r8    
          qtu_out(i,mkx:0:-1)          = 0._r8        
          thlu_out(i,mkx:0:-1)         = 0._r8         
          thvu_out(i,mkx:0:-1)         = 0._r8         
          uu_out(i,mkx:0:-1)           = 0._r8        
          vu_out(i,mkx:0:-1)           = 0._r8        
          qtu_emf_out(i,mkx:0:-1)      = 0._r8         
          thlu_emf_out(i,mkx:0:-1)     = 0._r8         
          uu_emf_out(i,mkx:0:-1)       = 0._r8          
          vu_emf_out(i,mkx:0:-1)       = 0._r8    
          uemf_out(i,mkx:0:-1)         = 0._r8    

          dwten_out(i,mkx:1:-1)        = 0._r8    
          diten_out(i,mkx:1:-1)        = 0._r8    
          flxrain_out(i,mkx:0:-1)      = 0._r8     
          flxsnow_out(i,mkx:0:-1)      = 0._r8    
          ntraprd_out(i,mkx:1:-1)      = 0._r8    
          ntsnprd_out(i,mkx:1:-1)      = 0._r8    

          excessu_arr_out(i,mkx:1:-1)  = 0._r8    
          excess0_arr_out(i,mkx:1:-1)  = 0._r8    
          xc_arr_out(i,mkx:1:-1)       = 0._r8    
          aquad_arr_out(i,mkx:1:-1)    = 0._r8    
          bquad_arr_out(i,mkx:1:-1)    = 0._r8    
          cquad_arr_out(i,mkx:1:-1)    = 0._r8    
          bogbot_arr_out(i,mkx:1:-1)   = 0._r8    
          bogtop_arr_out(i,mkx:1:-1)   = 0._r8    

          DO m = 1, ncnst
             trten_out(i,:mkx,m)       = 0._r8
             trflx_out(i,mkx:0:-1,m)   = 0._r8  
             tru_out(i,mkx:0:-1,m)     = 0._r8
             tru_emf_out(i,mkx:0:-1,m) = 0._r8
          ENDDO

       END IF

    END DO                  ! end of big i loop for each column.

    ! ---------------------------------------- !
    ! Writing main diagnostic output variables !
    ! ---------------------------------------- !


    RETURN

  END SUBROUTINE compute_uwshcu

  ! ------------------------------ !
  !                                ! 
  ! Beginning of subroutine blocks !
  !                                !
  ! ------------------------------ !

  SUBROUTINE getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
    ! ----------------------------------------------------------- !
    ! Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and !
    ! 'cinlcl, plfc' if any. Assume 'thv' is linear in each layer !
    ! both for cumulus and environment. Note that this subroutine !
    ! only include positive CIN in calculation - if there are any !
    ! negative CIN, it is assumed to be zero.    This is slightly !
    ! different from 'single_cin' below, where both positive  and !
    ! negative CIN are included.                                  !
    ! ----------------------------------------------------------- !
    REAL(r8) pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin,frc

    IF( thvubot .GT. thv0bot .AND. thvutop .GT. thv0top ) THEN
       plfc = pbot
       RETURN
    ELSEIF( thvubot .LE. thv0bot .AND. thvutop .LE. thv0top ) THEN 
       cin  = cin - ( (thvubot/thv0bot - 1._r8) + (thvutop/thv0top - 1._r8)) * (pbot - ptop) /        &
            ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)) )
    ELSEIF( thvubot .GT. thv0bot .AND. thvutop .LE. thv0top ) THEN 
       frc  = ( thvutop/thv0top - 1._r8 ) / ( (thvutop/thv0top - 1._r8) - (thvubot/thv0bot - 1._r8) )
       cin  = cin - ( thvutop/thv0top - 1._r8 ) * ( (ptop + frc*(pbot - ptop)) - ptop ) /             &
            ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)) )
    ELSE            
       frc  = ( thvubot/thv0bot - 1._r8 ) / ( (thvubot/thv0bot - 1._r8) - (thvutop/thv0top - 1._r8) )
       plfc = pbot - frc * ( pbot - ptop )
       cin  = cin - ( thvubot/thv0bot - 1._r8)*(pbot - plfc)/                                         & 
            ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top * exnf(ptop)))
    ENDIF

    RETURN
  END SUBROUTINE getbuoy

  FUNCTION single_cin(pbot,thv0bot,ptop,thv0top,thvubot,thvutop)
    ! ------------------------------------------------------- !
    ! Function to calculate a single layer CIN by summing all ! 
    ! positive and negative CIN.                              !
    ! ------------------------------------------------------- ! 
    REAL(r8) :: single_cin
    REAL(r8)    pbot,thv0bot,ptop,thv0top,thvubot,thvutop 

    single_cin = ( (1._r8 - thvubot/thv0bot) + (1._r8 - thvutop/thv0top)) * ( pbot - ptop ) / &
         ( pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)) )
    RETURN
  END FUNCTION single_cin


  SUBROUTINE conden(p,thl,qt,th,qv,ql,qi,rvls,id_check,qsat)
    ! --------------------------------------------------------------------- !
    ! Calculate thermodynamic properties from a given set of ( p, thl, qt ) !
    ! --------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL(r8), INTENT(in)  :: p
    REAL(r8), INTENT(in)  :: thl
    REAL(r8), INTENT(in)  :: qt
    REAL(r8), INTENT(out) :: th
    REAL(r8), INTENT(out) :: qv
    REAL(r8), INTENT(out) :: ql
    REAL(r8), INTENT(out) :: qi
    REAL(r8), INTENT(out) :: rvls
    INTEGER , INTENT(out) :: id_check
    INTEGER , EXTERNAL    :: qsat
    REAL(r8)              :: tc,temps
    REAL(r8)              :: leff, nu, qc
    INTEGER               :: iteration
    REAL(r8)              :: es(1)              ! Saturation vapor pressure
    REAL(r8)              :: qs(1)              ! Saturation spec. humidity
    REAL(r8)              :: gam(1)             ! (L/cp)*dqs/dT
    INTEGER               :: status             ! Return status of qsat call

    tc   = thl*exnf(p)
    ! Modification : In order to be compatible with the dlf treatment in stratiform.F90,
    !                we may use ( 268.15, 238.15 ) with 30K ramping instead of 20 K,
    !                in computing ice fraction below. 
    !                Note that 'cldwat_fice' uses ( 243.15, 263.15 ) with 20K ramping for stratus.
    nu   = MAX(MIN((268._r8 - tc)/20._r8,1.0_r8),0.0_r8)  ! Fraction of ice in the condensate. 
    leff = (1._r8 - nu)*xlv + nu*xls                      ! This is an estimate that hopefully speeds convergence

    ! --------------------------------------------------------------------------- !
    ! Below "temps" and "rvls" are just initial guesses for iteration loop below. !
    ! Note that the output "temps" from the below iteration loop is "temperature" !
    ! NOT "liquid temperature".                                                   !
    ! --------------------------------------------------------------------------- !

    temps  = tc
    status = qsat(temps,p,es(1),qs(1),gam(1), 1)
    rvls   = qs(1)

    IF( qs(1) .GE. qt ) THEN  
       id_check = 0
       qv = qt
       qc = 0._r8
       ql = 0._r8
       qi = 0._r8
       th = tc/exnf(p)
    ELSE 
       DO iteration = 1, 10
          temps  = temps + ( (tc-temps)*cp/leff + qt - rvls )/( cp/leff + ep2*leff*rvls/r/temps/temps )
          status = qsat(temps,p,es(1),qs(1),gam(1),1)
          rvls   = qs(1)
       END DO
       qc = MAX(qt - qs(1),0._r8)
       qv = qt - qc
       ql = qc*(1._r8 - nu)
       qi = nu*qc
       th = temps/exnf(p)
       IF( ABS((temps-(leff/cp)*qc)-tc) .GE. 1._r8 ) THEN
          id_check = 1
       ELSE
          id_check = 0
       END IF
    END IF

    RETURN
  END SUBROUTINE conden

  SUBROUTINE roots(a,b,c,r1,r2,status)
    ! --------------------------------------------------------- !
    ! Subroutine to solve the second order polynomial equation. !
    ! I should check this subroutine later.                     !
    ! --------------------------------------------------------- !
    REAL(r8), INTENT(in)  :: a
    REAL(r8), INTENT(in)  :: b
    REAL(r8), INTENT(in)  :: c
    REAL(r8), INTENT(out) :: r1
    REAL(r8), INTENT(out) :: r2
    INTEGER , INTENT(out) :: status
    REAL(r8)              :: q

    status = 0

    IF( a .EQ. 0._r8 ) THEN                            ! Form b*x + c = 0
       IF( b .EQ. 0._r8 ) THEN                        ! Failure: c = 0
          status = 1
       ELSE                                           ! b*x + c = 0
          r1 = -c/b
       ENDIF
       r2 = r1
    ELSE
       IF( b .EQ. 0._r8 ) THEN                        ! Form a*x**2 + c = 0
          IF( a*c .GT. 0._r8 ) THEN                  ! Failure: x**2 = -c/a < 0
             status = 2  
          ELSE                                       ! x**2 = -c/a 
             r1 = SQRT(-c/a)
          ENDIF
          r2 = -r1
       ELSE                                            ! Form a*x**2 + b*x + c = 0
          IF( (b**2 - 4._r8*a*c) .LT. 0._r8 ) THEN   ! Failure, no real roots
             status = 3
          ELSE
             q  = -0.5_r8*(b + SIGN(1.0_r8,b)*SQRT(b**2 - 4._r8*a*c))
             r1 =  q/a
             r2 =  c/q
          ENDIF
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE roots

  FUNCTION slope(mkx,field,p0)
    ! ------------------------------------------------------------------ !
    ! Function performing profile reconstruction of conservative scalars !
    ! in each layer. This is identical to profile reconstruction used in !
    ! UW-PBL scheme but from bottom to top layer here.     At the lowest !
    ! layer near to surface, slope is defined using the two lowest layer !
    ! mid-point values. I checked this subroutine and it is correct.     !
    ! ------------------------------------------------------------------ !
    REAL(r8)             :: slope(mkx)
    INTEGER,  INTENT(in) :: mkx
    REAL(r8), INTENT(in) :: field(mkx)
    REAL(r8), INTENT(in) :: p0(mkx)

    REAL(r8)             :: below
    REAL(r8)             :: above
    INTEGER              :: k

    below = ( field(2) - field(1) ) / ( p0(2) - p0(1) )
    DO k = 2, mkx
       above = ( field(k) - field(k-1) ) / ( p0(k) - p0(k-1) )
       IF( above .GT. 0._r8 ) THEN
          slope(k-1) = MAX(0._r8,MIN(above,below))
       ELSE 
          slope(k-1) = MIN(0._r8,MAX(above,below))
       END IF
       below = above
    END DO
    slope(mkx) = slope(mkx-1)

    RETURN
  END FUNCTION slope

  FUNCTION qsinvert(qt,thl,psfc,qsat)
    ! ----------------------------------------------------------------- !
    ! Function calculating saturation pressure ps (or pLCL) from qt and !
    ! thl ( liquid potential temperature,  NOT liquid virtual potential ! 
    ! temperature) by inverting Bolton formula. I should check later if !
    ! current use of 'leff' instead of 'xlv' here is reasonable or not. !
    ! ----------------------------------------------------------------- !
    REAL(r8)          :: qsinvert    
    REAL(r8)             qt, thl, psfc
    REAL(r8)             ps, Pis, Ts, err, dlnqsdT, dTdPis
    REAL(r8)             dPisdps, dlnqsdps, derrdps, dps 
    REAL(r8)             Ti, rhi, TLCL, PiLCL, psmin, dpsmax
    INTEGER              i
    INTEGER, EXTERNAL :: qsat
    REAL(r8)          :: es(1)                     ! saturation vapor pressure
    REAL(r8)          :: qs(1)                     ! saturation spec. humidity
    REAL(r8)          :: gam(1)                    ! (L/cp)*dqs/dT
    INTEGER           :: status                    ! return status of qsat call
    REAL(r8)          :: leff, nu

    psmin  = 100._r8*100._r8 ! Default saturation pressure [Pa] if iteration does not converge
    dpsmax = 1._r8           ! Tolerance [Pa] for convergence of iteration

    ! ------------------------------------ !
    ! Calculate best initial guess of pLCL !
    ! ------------------------------------ !

    Ti       =  thl*(psfc/p00)**rovcp
    status   =  qsat(Ti,psfc,es(1),qs(1),gam(1),1)
    rhi      =  qt/qs(1)      
    IF( rhi .LE. 0.01_r8 ) THEN
       WRITE(iulog,*) 'Source air is too dry and pLCL is set to psmin in uwshcu.F90' 
       qsinvert = psmin
       RETURN
    END IF
    TLCL     =  55._r8 + 1._r8/(1._r8/(Ti-55._r8)-LOG(rhi)/2840._r8); ! Bolton's formula. MWR.1980.Eq.(22)
    PiLCL    =  TLCL/thl
    ps       =  p00*(PiLCL)**(1._r8/rovcp)

    DO i = 1, 10
       Pis      =  (ps/p00)**rovcp
       Ts       =  thl*Pis
       status   =  qsat(Ts,ps,es(1),qs(1),gam(1),1)
       err      =  qt - qs(1)
       nu       =  MAX(MIN((268._r8 - Ts)/20._r8,1.0_r8),0.0_r8)        
       leff     =  (1._r8 - nu)*xlv + nu*xls                   
       dlnqsdT  =  gam(1)*(cp/leff)/qs(1)
       dTdPis   =  thl
       dPisdps  =  rovcp*Pis/ps 
       dlnqsdps = -1._r8/(ps - (1._r8 - ep2)*es(1))
       derrdps  = -qs(1)*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
       dps      = -err/derrdps
       ps       =  ps + dps
       IF( ps .LT. 0._r8 ) THEN
          WRITE(iulog,*) 'pLCL iteration is negative and set to psmin in uwshcu.F90', qt, thl, psfc 
          qsinvert = psmin
          RETURN    
       END IF
       IF( ABS(dps) .LE. dpsmax ) THEN
          qsinvert = ps
          RETURN
       END IF
    END DO
    WRITE(iulog,*) 'pLCL does not converge and is set to psmin in uwshcu.F90', qt, thl, psfc 
    qsinvert = psmin
    RETURN
  END FUNCTION qsinvert

  REAL(r8) FUNCTION compute_alpha(del_CIN,ke)
    ! ------------------------------------------------ !
    ! Subroutine to compute proportionality factor for !
    ! implicit CIN calculation.                        !   
    ! ------------------------------------------------ !
    REAL(r8) :: del_CIN, ke
    REAL(r8) :: x0, x1

    INTEGER  :: iteration

    x0 = 0._r8
    DO iteration = 1, 10
       x1 = x0 - (EXP(-x0*ke*del_CIN) - x0)/(-ke*del_CIN*EXP(-x0*ke*del_CIN) - 1._r8)
       x0 = x1
    END DO
    compute_alpha = x0

    RETURN

  END FUNCTION compute_alpha

  REAL(r8) FUNCTION compute_mumin2(mulcl,rmaxfrac,mulow)
    ! --------------------------------------------------------- !
    ! Subroutine to compute critical 'mu' (normalized CIN) such ! 
    ! that updraft fraction at the LCL is equal to 'rmaxfrac'.  !
    ! --------------------------------------------------------- !  
    REAL(r8) :: mulcl, rmaxfrac, mulow
    REAL(r8) :: x0, x1, ex, ef, exf, f, fs
    INTEGER  :: iteration

    x0 = mulow
    DO iteration = 1, 10
       ex = EXP(-x0**2)
       ef = erfc(x0)
       ! if(x0.ge.3._r8) then
       !    compute_mumin2 = 3._r8 
       !    goto 20
       ! endif 
       exf = ex/ef
       f  = 0.5_r8*exf**2 - 0.5_r8*(ex/2._r8/rmaxfrac)**2 - (mulcl*2.5066_r8/2._r8)**2
       fs = (2._r8*exf**2)*(exf/SQRT(3.141592_r8)-x0) + (0.5_r8*x0*ex**2)/(rmaxfrac**2)
       x1 = x0 - f/fs     
       x0 = x1
    END DO
    compute_mumin2 = x0

20  RETURN

  END FUNCTION compute_mumin2

  REAL(r8) FUNCTION compute_ppen(wtwb,D,bogbot,bogtop,rho0j,dpen)
    ! ----------------------------------------------------------- !
    ! Subroutine to compute critical 'ppen[Pa]<0' ( pressure dis. !
    ! from 'ps0(kpen-1)' to the cumulus top where cumulus updraft !
    ! vertical velocity is exactly zero ) by considering exact    !
    ! non-zero fer(kpen).                                         !  
    ! ----------------------------------------------------------- !  
    REAL(r8) :: wtwb, D, bogbot, bogtop, rho0j, dpen
    REAL(r8) :: x0, x1, f, fs, SB, s00
    INTEGER  :: iteration

    ! Buoyancy slope
    SB = ( bogtop - bogbot ) / dpen
    ! Sign of slope, 'f' at x = 0
    ! If 's00>0', 'w' increases with height.
    s00 = bogbot / rho0j - D * wtwb

    IF( D*dpen .LT. 1.e-8 ) THEN
       IF( s00 .GE. 0._r8 ) THEN
          x0 = dpen       
       ELSE
          x0 = MAX(0._r8,MIN(dpen,-0.5_r8*wtwb/s00))
       ENDIF
    ELSE
       IF( s00 .GE. 0._r8 ) THEN
          x0 = dpen
       ELSE 
          x0 = 0._r8
       ENDIF
       DO iteration = 1, 5
          f  = EXP(-2._r8*D*x0)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
               (SB*x0+bogbot-SB/(2._r8*D))/(D*rho0j)
          fs = -2._r8*D*EXP(-2._r8*D*x0)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
               (SB)/(D*rho0j)
          x1 = x0 - f/fs     
          x0 = x1
       END DO

    ENDIF

    compute_ppen = -MAX(0._r8,MIN(dpen,x0))

  END FUNCTION compute_ppen

  SUBROUTINE fluxbelowinv(cbmf,ps0,mkx,kinv,dt,xsrc,xmean,xtopin,xbotin,xflx)   
    ! ------------------------------------------------------------------------- !
    ! Subroutine to calculate turbulent fluxes at and below 'kinv-1' interfaces.!
    ! Check in the main program such that input 'cbmf' should not be zero.      !  
    ! If the reconstructed inversion height does not go down below the 'kinv-1' !
    ! interface, then turbulent flux at 'kinv-1' interface  is simply a product !
    ! of 'cmbf' and 'qtsrc-xbot' where 'xbot' is the value at the top interface !
    ! of 'kinv-1' layer. This flux is linearly interpolated down to the surface !
    ! assuming turbulent fluxes at surface are zero. If reconstructed inversion !
    ! height goes down below the 'kinv-1' interface, subsidence warming &drying !
    ! measured by 'xtop-xbot', where  'xtop' is the value at the base interface !
    ! of 'kinv+1' layer, is added ONLY to the 'kinv-1' layer, using appropriate !
    ! mass weighting ( rpinv and rcbmf, or rr = rpinv / rcbmf ) between current !
    ! and next provisional time step. Also impose a limiter to enforce outliers !
    ! of thermodynamic variables in 'kinv' layer  to come back to normal values !
    ! at the next step.                                                         !
    ! ------------------------------------------------------------------------- !            
    INTEGER,  INTENT(in)                     :: mkx, kinv 
    REAL(r8), INTENT(in)                     :: cbmf, dt, xsrc, xmean, xtopin, xbotin
    REAL(r8), INTENT(in),  DIMENSION(0:mkx)  :: ps0
    REAL(r8), INTENT(out), DIMENSION(0:mkx)  :: xflx  
    INTEGER k
    REAL(r8) rcbmf, rpeff, dp, rr, pinv_eff, xtop, xbot, pinv, xtop_ori, xbot_ori

    xflx(0:mkx) = 0._r8
    dp = ps0(kinv-1) - ps0(kinv) 

    IF( ABS(xbotin-xtopin) .LE. 1.e-13_r8 ) THEN
       xbot = xbotin - 1.e-13_r8
       xtop = xtopin + 1.e-13_r8
    ELSE
       xbot = xbotin
       xtop = xtopin
    ENDIF
    ! -------------------------------------- !
    ! Compute reconstructed inversion height !
    ! -------------------------------------- !
    xtop_ori = xtop
    xbot_ori = xbot
    rcbmf = ( cbmf * g * dt ) / dp                  ! Can be larger than 1 : 'OK'      
    rpeff = ( xmean - xtop ) / ( xbot - xtop ) 
    rpeff = MIN( MAX(0._r8,rpeff), 1._r8 )          ! As of this, 0<= rpeff <= 1   
    IF( rpeff .EQ. 0._r8 .OR. rpeff .EQ. 1._r8 ) THEN
       xbot = xmean
       xtop = xmean
    ENDIF
    ! Below two commented-out lines are the old code replacing the above 'if' block.   
    ! if(rpeff.eq.1) xbot = xmean
    ! if(rpeff.eq.0) xtop = xmean    
    rr       = rpeff / rcbmf
    pinv     = ps0(kinv-1) - rpeff * dp             ! "pinv" before detraining mass
    pinv_eff = ps0(kinv-1) + ( rcbmf - rpeff ) * dp ! Effective "pinv" after detraining mass
    ! ----------------------------------------------------------------------- !
    ! Compute turbulent fluxes.                                               !
    ! Below two cases exactly converges at 'kinv-1' interface when rr = 1._r8 !
    ! ----------------------------------------------------------------------- !
    DO k = 0, kinv - 1
       xflx(k) = cbmf * ( xsrc - xbot ) * ( ps0(0) - ps0(k) ) / ( ps0(0) - pinv )
    END DO
    IF( rr .LE. 1._r8 ) THEN
       xflx(kinv-1) =  xflx(kinv-1) - ( 1._r8 - rr ) * cbmf * ( xtop_ori - xbot_ori )
    ENDIF

    RETURN
  END SUBROUTINE fluxbelowinv

  SUBROUTINE positive_moisture_single( xlv, xls, mkx, dt, qvmin, qlmin, qimin, dp, qv, ql, qi, s, qvten, qlten, qiten, sten )
    ! ------------------------------------------------------------------------------- !
    ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
    ! force them to be larger than minimum value by (1) condensating water vapor      !
    ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
    ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
    ! Update final state variables and tendencies associated with this correction.    !
    ! If any condensation happens, update (s,t) too.                                  !
    ! Note that (qv,ql,qi,s) are final state variables after applying corresponding   !
    ! input tendencies and corrective tendencies                                      !
    ! ------------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER,  INTENT(in)     :: mkx
    REAL(r8), INTENT(in)     :: xlv, xls
    REAL(r8), INTENT(in)     :: dt, qvmin, qlmin, qimin
    REAL(r8), INTENT(in)     :: dp(mkx)
    REAL(r8), INTENT(inout)  :: qv(mkx), ql(mkx), qi(mkx), s(mkx)
    REAL(r8), INTENT(inout)  :: qvten(mkx), qlten(mkx), qiten(mkx), sten(mkx)
    INTEGER   k
    REAL(r8)  dql, dqi, dqv, sum, aa, dum 

    DO k = mkx, 1, -1        ! From the top to the 1st (lowest) layer from the surface
       dql = MAX(0._r8,1._r8*qlmin-ql(k))
       dqi = MAX(0._r8,1._r8*qimin-qi(k))
       qlten(k) = qlten(k) +  dql/dt
       qiten(k) = qiten(k) +  dqi/dt
       qvten(k) = qvten(k) - (dql+dqi)/dt
       sten(k)  = sten(k)  + xlv * (dql/dt) + xls * (dqi/dt)
       ql(k)    = ql(k) +  dql
       qi(k)    = qi(k) +  dqi
       qv(k)    = qv(k) -  dql - dqi
       s(k)     = s(k)  +  xlv * dql + xls * dqi
       dqv      = MAX(0._r8,1._r8*qvmin-qv(k))
       qvten(k) = qvten(k) + dqv/dt
       qv(k)    = qv(k)   + dqv
       IF( k .NE. 1 ) THEN 
          qv(k-1)    = qv(k-1)    - dqv*dp(k)/dp(k-1)
          qvten(k-1) = qvten(k-1) - dqv*dp(k)/dp(k-1)/dt
       ENDIF
       qv(k) = MAX(qv(k),qvmin)
       ql(k) = MAX(ql(k),qlmin)
       qi(k) = MAX(qi(k),qimin)
    END DO
    ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture. 
    IF( dqv .GT. 1.e-20_r8 ) THEN
       sum = 0._r8
       DO k = 1, mkx
          IF( qv(k) .GT. 2._r8*qvmin ) sum = sum + qv(k)*dp(k)
       ENDDO
       aa = dqv*dp(1)/MAX(1.e-20_r8,sum)
       IF( aa .LT. 0.5_r8 ) THEN
          DO k = 1, mkx
             IF( qv(k) .GT. 2._r8*qvmin ) THEN
                dum      = aa*qv(k)
                qv(k)    = qv(k) - dum
                qvten(k) = qvten(k) - dum/dt
             ENDIF
          ENDDO
       ELSE 
          WRITE(iulog,*) 'Full positive_moisture is impossible in uwshcu'
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE positive_moisture_single

  SUBROUTINE findsp ( ncol,pver, q, t, p, tsp, qsp)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !     find the wet bulb temperature for a given t and q
    !     in a longitude height section
    !     wet bulb temp is the temperature and spec humidity that is 
    !     just saturated and has the same enthalpy
    !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
    !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
    !
    ! Method: 
    ! a Newton method is used
    ! first guess uses an algorithm provided by John Petch from the UKMO
    ! we exclude points where the physical situation is unrealistic
    ! e.g. where the temperature is outside the range of validity for the
    !      saturation vapor pressure, or where the water vapor pressure
    !      exceeds the ambient pressure, or the saturation specific humidity is 
    !      unrealistic
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------

    !use wv_saturation, only: estblf, hlatv, tmin, hlatf, rgasv, pcf, &
    !                         cp, epsqs, ttrice

    !
    !     input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: pver                  ! number of atmospheric columns

    REAL(r8), INTENT(in) :: q(ncol,pver)        ! water vapor (kg/kg)
    REAL(r8), INTENT(in) :: t(ncol,pver)        ! temperature (K)
    REAL(r8), INTENT(in) :: p(ncol,pver)        ! pressure    (Pa)
    !
    ! output arguments
    !
    REAL(r8), INTENT(out) :: tsp(ncol,pver)      ! saturation temp (K)
    REAL(r8), INTENT(out) :: qsp(ncol,pver)      ! saturation mixing ratio (kg/kg)
    !
    ! local variables
    !
    INTEGER i                 ! work variable
    INTEGER k                 ! work variable
    LOGICAL lflg              ! work variable
    INTEGER iter              ! work variable
    INTEGER l                 ! work variable
    LOGICAL :: error_found

    REAL(r8) omeps                ! 1 minus epsilon
    REAL(r8) trinv                ! work variable
    REAL(r8) es                   ! sat. vapor pressure
    REAL(r8) desdt                ! change in sat vap pressure wrt temperature
    !     real(r8) desdp                ! change in sat vap pressure wrt pressure
    REAL(r8) dqsdt                ! change in sat spec. hum. wrt temperature
    REAL(r8) dgdt                 ! work variable
    REAL(r8) g                    ! work variable
    REAL(r8) weight(ncol)        ! work variable
    REAL(r8) hlatsb               ! (sublimation)
    REAL(r8) hlatvp               ! (vaporization)
    REAL(r8) hltalt(ncol,pver)   ! lat. heat. of vap.
    REAL(r8) tterm                ! work var.
    REAL(r8) qs                   ! spec. hum. of water vapor
    REAL(r8) tc                   ! crit temp of transition to ice
    REAL(r8) tt0 

    ! work variables
    REAL(r8) t1, q1, dt, dq
    REAL(r8) dtm, dqm
    REAL(r8) qvd, a1, tmp
    REAL(r8) rair
    REAL(r8) r1b, c1, c2, c3
    REAL(r8) denom
    REAL(r8) dttol
    REAL(r8) dqtol
    INTEGER doit(ncol) 
    REAL(r8) enin(ncol), enout(ncol)
    REAL(r8) tlim(ncol)

    omeps = 1.0_r8 - epsqs
    trinv = 1.0_r8/ttrice
    a1 = 7.5_r8*LOG(10._r8)
    rair =  287.04_r8
    c3 = rair*a1/cp
    dtm = 0._r8    ! needed for iter=0 blowup with f90 -ei
    dqm = 0._r8    ! needed for iter=0 blowup with f90 -ei
    dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
    dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
    tt0 = 273.15_r8  ! Freezing temperature 
    !  tmin = 173.16 ! the coldest temperature we can deal with
    !
    ! max number of times to iterate the calculation
    iter = 50
    !
    DO k = 1,pver

       !
       ! first guess on the wet bulb temperature
       !
       DO i = 1,ncol

          !#ifdef DEBUG
          !         if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
          !            write(iulog,*) ' '
          !            write(iulog,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
          !         endif
          !#endif
          ! limit the temperature range to that relevant to the sat vap pres tables
          IF (WACCM_MOZART )THEN
             tlim(i) = MIN(MAX(t(i,k),173._r8),373._r8)
          ELSE
             tlim(i) = MIN(MAX(t(i,k),128._r8),373._r8)
          ENDIF
          es = estblf(tlim(i))
          denom = p(i,k) - omeps*es
          qs = epsqs*es/denom
          doit(i) = 0
          enout(i) = 1._r8
          ! make sure a meaningful calculation is possible
          IF (p(i,k) > 5._r8*es .AND. qs > 0._r8 .AND. qs < 0.5_r8) THEN
             !
             ! Saturation specific humidity
             !
             qs = MIN(epsqs*es/denom,1._r8)
             !
             ! "generalized" analytic expression for t derivative of es
             !  accurate to within 1 percent for 173.16 < t < 373.16
             !
             ! Weighting of hlat accounts for transition from water to ice
             ! polynomial expression approximates difference between es over
             ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
             ! -40): required for accurate estimate of es derivative in transition
             ! range from ice to water also accounting for change of hlatv with t
             ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
             !
             tc     = tlim(i) - tt0
             lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
             weight(i) = MIN(-tc*trinv,1.0_r8)
             hlatsb = hlatv + weight(i)*hlatf
             hlatvp = hlatv - 2369.0_r8*tc
             IF (tlim(i) < tt0) THEN
                hltalt(i,k) = hlatsb
             ELSE
                hltalt(i,k) = hlatvp
             END IF
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

             ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp =  q(i,k) - qs
             c1 = hltalt(i,k)*c3
             c2 = (tlim(i) + 36._r8)**2
             r1b    = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
             !#ifdef DEBUG
             !             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
             !                write(iulog,*) ' relative humidity ', q(i,k)/qs
             !                write(iulog,*) ' first guess ', tsp(i,k)
             !             endif
             !#endif
             es = estblf(tsp(i,k))
             qsp(i,k) = MIN(epsqs*es/(p(i,k) - omeps*es),1._r8)
          ELSE
             doit(i) = 1
             tsp(i,k) = tlim(i)
             qsp(i,k) = q(i,k)
             enin(i) = 1._r8
          ENDIF
       END DO   ! end do i
       !
       ! now iterate on first guess
       !
       DO l = 1, iter
          dtm = 0
          dqm = 0
          DO i = 1,ncol
             IF (doit(i) == 0) THEN
                es = estblf(tsp(i,k))
                !
                ! Saturation specific humidity
                !
                qs = MIN(epsqs*es/(p(i,k) - omeps*es),1._r8)
                !
                ! "generalized" analytic expression for t derivative of es
                ! accurate to within 1 percent for 173.16 < t < 373.16
                !
                ! Weighting of hlat accounts for transition from water to ice
                ! polynomial expression approximates difference between es over
                ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                ! -40): required for accurate estimate of es derivative in transition
                ! range from ice to water also accounting for change of hlatv with t
                ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                !
                tc     = tsp(i,k) - tt0
                lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
                weight(i) = MIN(-tc*trinv,1.0_r8)
                hlatsb = hlatv + weight(i)*hlatf
                hlatvp = hlatv - 2369.0_r8*tc
                IF (tsp(i,k) < tt0) THEN
                   hltalt(i,k) = hlatsb
                ELSE
                   hltalt(i,k) = hlatvp
                END IF
                IF (lflg) THEN
                   tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
                ELSE
                   tterm = 0.0_r8
                END IF
                desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
                !              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
                g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
                dgdt = -(cp + hltalt(i,k)*dqsdt)
                t1 = tsp(i,k) - g/dgdt
                dt = ABS(t1 - tsp(i,k))/t1
                tsp(i,k) = MAX(t1,tmin)
                es = estblf(tsp(i,k))
                q1 = MIN(epsqs*es/(p(i,k) - omeps*es),1._r8)
                dq = ABS(q1 - qsp(i,k))/MAX(q1,1.e-12_r8)
                qsp(i,k) = q1
                !#ifdef DEBUG
                !               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                !                  write(iulog,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
                !               endif
                !#endif
                dtm = MAX(dtm,dt)
                dqm = MAX(dqm,dq)
                ! if converged at this point, exclude it from more iterations
                IF (dt < dttol .AND. dq < dqtol) THEN
                   doit(i) = 2
                ENDIF
                enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
                ! bail out if we are too near the end of temp range
                IF ( WACCM_MOZART )THEN
                   IF (tsp(i,k) < 174.16_r8) THEN
                      doit(i) = 4
                   ENDIF

                ELSE
                   IF (tsp(i,k) < 130.16_r8) THEN
                      doit(i) = 4
                   ENDIF

                ENDIF
                !                  doit(i) = 4
                !               endif
             ELSE
             ENDIF
          END DO              ! do i = 1,ncol

          IF (dtm < dttol .AND. dqm < dqtol) THEN
             go to 10
          ENDIF

       END DO                 ! do l = 1,iter
10     CONTINUE

       error_found = .FALSE.
       IF (dtm > dttol .OR. dqm > dqtol) THEN
          DO i = 1,ncol
             IF (doit(i) == 0) error_found = .TRUE.
          END DO
          IF (error_found) THEN
             DO i = 1,ncol
                IF (doit(i) == 0) THEN
                   WRITE(iulog,*) ' findsp not converging at point i, k ', i, k
                   WRITE(iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                   WRITE(iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                   CALL endrun ('FINDSP')
                ENDIF
             END DO
          ENDIF
       ENDIF
       DO i = 1,ncol
          IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
             error_found = .TRUE.
          ENDIF
       END DO
       IF (error_found) THEN
          DO i = 1,ncol
             IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
                WRITE(iulog,*) ' the enthalpy is not conserved for point ', &
                     i, k, enin(i), enout(i)
                WRITE(iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                WRITE(iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                CALL endrun ('FINDSP')
             ENDIF
          END DO
       ENDIF

    END DO                    ! level loop (k=1,pver)

    RETURN
  END SUBROUTINE findsp

  ! ------------------------ !
  !                          ! 
  ! End of subroutine blocks !
  !                          !
  ! ------------------------ !
  !------------------------------------------------------------------------------------------

  FUNCTION DERFC(X)
    !--------------------------------------------------------------------
    !
    ! This subprogram computes approximate values for erfc(x).
    !   (see comments heading CALERF).
    !
    !   Author/date: W. J. Cody, January 8, 1985
    !
    !--------------------------------------------------------------------
    INTEGER, PARAMETER :: rk = r8 ! 8 byte real

    ! argument
    REAL(rk), INTENT(in) :: X

    ! return value
    REAL(rk) :: DERFC

    ! local variables
    INTEGER :: JINT = 1
    !------------------------------------------------------------------

    CALL CALERF_r8(X, DERFC, JINT)
  END FUNCTION DERFC

  !------------------------------------------------------------------------------------------

  FUNCTION ERFC_r4(X)
    !--------------------------------------------------------------------
    !
    ! This subprogram computes approximate values for erfc(x).
    !   (see comments heading CALERF).
    !
    !   Author/date: W. J. Cody, January 8, 1985
    !
    !--------------------------------------------------------------------
    INTEGER, PARAMETER :: rk = r4 ! 4 byte real

    ! argument
    REAL(rk), INTENT(in) :: X

    ! return value
    REAL(rk) :: ERFC_r4

    ! local variables
    INTEGER :: JINT = 1
    !------------------------------------------------------------------

    CALL CALERF_r4(X, ERFC_r4, JINT)
  END FUNCTION ERFC_r4
  !------------------------------------------------------------------
  !
  ! 6 December 2006 -- B. Eaton
  ! The following comments are from the original version of CALERF.
  ! The only changes in implementing this module are that the function
  ! names previously used for the single precision versions have been
  ! adopted for the new generic interfaces.  To support these interfaces
  ! there is now both a single precision version (calerf_r4) and a
  ! double precision version (calerf_r8) of CALERF below.  These versions
  ! are hardcoded to use IEEE arithmetic.
  !
  !------------------------------------------------------------------
  !
  ! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
  !   for a real argument  x.  It contains three FUNCTION type
  !   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
  !   and one SUBROUTINE type subprogram, CALERF.  The calling
  !   statements for the primary entries are:
  !
  !                   Y=ERF(X)     (or   Y=DERF(X)),
  !
  !                   Y=ERFC(X)    (or   Y=DERFC(X)),
  !   and
  !                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
  !
  !   The routine  CALERF  is intended for internal packet use only,
  !   all computations within the packet being concentrated in this
  !   routine.  The function subprograms invoke  CALERF  with the
  !   statement
  !
  !          CALL CALERF(ARG,RESULT,JINT)
  !
  !   where the parameter usage is as follows
  !
  !      Function                     Parameters for CALERF
  !       call              ARG                  Result          JINT
  !
  !     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
  !     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
  !     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
  !
  !   The main computation evaluates near-minimax approximations
  !   from "Rational Chebyshev approximations for the error function"
  !   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
  !   transportable program uses rational functions that theoretically
  !   approximate  erf(x)  and  erfc(x)  to at least 18 significant
  !   decimal digits.  The accuracy achieved depends on the arithmetic
  !   system, the compiler, the intrinsic functions, and proper
  !   selection of the machine-dependent constants.
  !
  !*******************************************************************
  !*******************************************************************
  !
  ! Explanation of machine-dependent constants
  !
  !   XMIN   = the smallest positive floating-point number.
  !   XINF   = the largest positive finite floating-point number.
  !   XNEG   = the largest negative argument acceptable to ERFCX;
  !            the negative of the solution to the equation
  !            2*exp(x*x) = XINF.
  !   XSMALL = argument below which erf(x) may be represented by
  !            2*x/sqrt(pi)  and above which  x*x  will not underflow.
  !            A conservative value is the largest machine number X
  !            such that   1.0 + X = 1.0   to machine precision.
  !   XBIG   = largest argument acceptable to ERFC;  solution to
  !            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
  !            W(x) = exp(-x*x)/[x*sqrt(pi)].
  !   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
  !            machine precision.  A conservative value is
  !            1/[2*sqrt(XSMALL)]
  !   XMAX   = largest acceptable argument to ERFCX; the minimum
  !            of XINF and 1/[sqrt(pi)*XMIN].
  !
  !   Approximate values for some important machines are:
  !
  !                          XMIN       XINF        XNEG     XSMALL
  !
  !  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
  !  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
  !  IEEE (IBM/XT,
  !    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
  !  IEEE (IBM/XT,
  !    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
  !  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
  !  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
  !  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
  !  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
  !
  !
  !                          XBIG       XHUGE       XMAX
  !
  !  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
  !  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
  !  IEEE (IBM/XT,
  !    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
  !  IEEE (IBM/XT,
  !    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
  !  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
  !  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
  !  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
  !  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
  !
  !*******************************************************************
  !*******************************************************************
  !
  ! Error returns
  !
  !  The program returns  ERFC = 0      for  ARG .GE. XBIG;
  !
  !                       ERFCX = XINF  for  ARG .LT. XNEG;
  !      and
  !                       ERFCX = 0     for  ARG .GE. XMAX.
  !
  !
  ! Intrinsic functions required are:
  !
  !     ABS, AINT, EXP
  !
  !
  !  Author: W. J. Cody
  !          Mathematics and Computer Science Division
  !          Argonne National Laboratory
  !          Argonne, IL 60439
  !
  !  Latest modification: March 19, 1990
  !
  !------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  SUBROUTINE CALERF_r8(ARG, RESULT, JINT)

    !------------------------------------------------------------------
    !  This version uses 8-byte reals
    !------------------------------------------------------------------
    INTEGER, PARAMETER :: rk = r8

    ! arguments
    REAL(rk), INTENT(in)  :: arg
    INTEGER,  INTENT(in)  :: jint
    REAL(rk), INTENT(out) :: RESULT

    ! local variables
    INTEGER :: I

    REAL(rk) :: X, Y, YSQ, XNUM, XDEN, DEL

    !------------------------------------------------------------------
    !  Mathematical constants
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: ZERO   = 0.0E0_rk
    REAL(rk), PARAMETER :: FOUR   = 4.0E0_rk
    REAL(rk), PARAMETER :: ONE    = 1.0E0_rk
    REAL(rk), PARAMETER :: HALF   = 0.5E0_rk
    REAL(rk), PARAMETER :: TWO    = 2.0E0_rk
    REAL(rk), PARAMETER :: SQRPI  = 5.6418958354775628695E-1_rk
    REAL(rk), PARAMETER :: THRESH = 0.46875E0_rk
    REAL(rk), PARAMETER :: SIXTEN = 16.0E0_rk

    !------------------------------------------------------------------
    !  Machine-dependent constants: IEEE single precision values
    !------------------------------------------------------------------
    !S      real, parameter :: XINF   =  3.40E+38
    !S      real, parameter :: XNEG   = -9.382E0
    !S      real, parameter :: XSMALL =  5.96E-8 
    !S      real, parameter :: XBIG   =  9.194E0
    !S      real, parameter :: XHUGE  =  2.90E3
    !S      real, parameter :: XMAX   =  4.79E37

    !------------------------------------------------------------------
    !  Machine-dependent constants: IEEE double precision values
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: XINF   =   1.79E308_r8
    REAL(rk), PARAMETER :: XNEG   = -26.628E0_r8
    REAL(rk), PARAMETER :: XSMALL =   1.11E-16_r8
    REAL(rk), PARAMETER :: XBIG   =  26.543E0_r8
    REAL(rk), PARAMETER :: XHUGE  =   6.71E7_r8
    REAL(rk), PARAMETER :: XMAX   =   2.53E307_r8

    !------------------------------------------------------------------
    !  Coefficients for approximation to  erf  in first interval
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: A(5) = (/ 3.16112374387056560E00_rk, 1.13864154151050156E02_rk, &
         3.77485237685302021E02_rk, 3.20937758913846947E03_rk, &
         1.85777706184603153E-1_rk /)
    REAL(rk), PARAMETER :: B(4) = (/ 2.36012909523441209E01_rk, 2.44024637934444173E02_rk, &
         1.28261652607737228E03_rk, 2.84423683343917062E03_rk /)

    !------------------------------------------------------------------
    !  Coefficients for approximation to  erfc  in second interval
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: C(9) = (/ 5.64188496988670089E-1_rk, 8.88314979438837594E00_rk, &
         6.61191906371416295E01_rk, 2.98635138197400131E02_rk, &
         8.81952221241769090E02_rk, 1.71204761263407058E03_rk, &
         2.05107837782607147E03_rk, 1.23033935479799725E03_rk, &
         2.15311535474403846E-8_rk /)
    REAL(rk), PARAMETER :: D(8) = (/ 1.57449261107098347E01_rk, 1.17693950891312499E02_rk, &
         5.37181101862009858E02_rk, 1.62138957456669019E03_rk, &
         3.29079923573345963E03_rk, 4.36261909014324716E03_rk, &
         3.43936767414372164E03_rk, 1.23033935480374942E03_rk /)

    !------------------------------------------------------------------
    !  Coefficients for approximation to  erfc  in third interval
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: P(6) = (/ 3.05326634961232344E-1_rk, 3.60344899949804439E-1_rk, &
         1.25781726111229246E-1_rk, 1.60837851487422766E-2_rk, &
         6.58749161529837803E-4_rk, 1.63153871373020978E-2_rk /)
    REAL(rk), PARAMETER :: Q(5) = (/ 2.56852019228982242E00_rk, 1.87295284992346047E00_rk, &
         5.27905102951428412E-1_rk, 6.05183413124413191E-2_rk, &
         2.33520497626869185E-3_rk /)

    !------------------------------------------------------------------
    X = ARG
    Y = ABS(X)
    IF (Y .LE. THRESH) THEN
       !------------------------------------------------------------------
       !  Evaluate  erf  for  |X| <= 0.46875
       !------------------------------------------------------------------
       YSQ = ZERO
       IF (Y .GT. XSMALL) YSQ = Y * Y
       XNUM = A(5)*YSQ
       XDEN = YSQ
       DO I = 1, 3
          XNUM = (XNUM + A(I)) * YSQ
          XDEN = (XDEN + B(I)) * YSQ
       END DO
       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
       IF (JINT .NE. 0) RESULT = ONE - RESULT
       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
       GO TO 80
    ELSE IF (Y .LE. FOUR) THEN
       !------------------------------------------------------------------
       !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
       !------------------------------------------------------------------
       XNUM = C(9)*Y
       XDEN = Y
       DO I = 1, 7
          XNUM = (XNUM + C(I)) * Y
          XDEN = (XDEN + D(I)) * Y
       END DO
       RESULT = (XNUM + C(8)) / (XDEN + D(8))
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
    ELSE
       !------------------------------------------------------------------
       !  Evaluate  erfc  for |X| > 4.0
       !------------------------------------------------------------------
       RESULT = ZERO
       IF (Y .GE. XBIG) THEN
          IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 30
          IF (Y .GE. XHUGE) THEN
             RESULT = SQRPI / Y
             GO TO 30
          END IF
       END IF
       YSQ = ONE / (Y * Y)
       XNUM = P(6)*YSQ
       XDEN = YSQ
       DO I = 1, 4
          XNUM = (XNUM + P(I)) * YSQ
          XDEN = (XDEN + Q(I)) * YSQ
       END DO
       RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
       RESULT = (SQRPI -  RESULT) / Y
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
    END IF
30  CONTINUE
    !------------------------------------------------------------------
    !  Fix up for negative argument, erf, etc.
    !------------------------------------------------------------------
    IF (JINT .EQ. 0) THEN
       RESULT = (HALF - RESULT) + HALF
       IF (X .LT. ZERO) RESULT = -RESULT
    ELSE IF (JINT .EQ. 1) THEN
       IF (X .LT. ZERO) RESULT = TWO - RESULT
    ELSE
       IF (X .LT. ZERO) THEN
          IF (X .LT. XNEG) THEN
             RESULT = XINF
          ELSE
             YSQ = AINT(X*SIXTEN)/SIXTEN
             DEL = (X-YSQ)*(X+YSQ)
             Y = EXP(YSQ*YSQ) * EXP(DEL)
             RESULT = (Y+Y) - RESULT
          END IF
       END IF
    END IF
80  CONTINUE
  END SUBROUTINE CALERF_r8

  !------------------------------------------------------------------------------------------

  SUBROUTINE CALERF_r4(ARG, RESULT, JINT)

    !------------------------------------------------------------------
    !  This version uses 4-byte reals
    !------------------------------------------------------------------
    INTEGER, PARAMETER :: rk = r4

    ! arguments
    REAL(rk), INTENT(in)  :: arg
    INTEGER,  INTENT(in)  :: jint
    REAL(rk), INTENT(out) :: RESULT

    ! local variables
    INTEGER :: I

    REAL(rk) :: X, Y, YSQ, XNUM, XDEN, DEL

    !------------------------------------------------------------------
    !  Mathematical constants
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: ZERO   = 0.0E0_rk
    REAL(rk), PARAMETER :: FOUR   = 4.0E0_rk
    REAL(rk), PARAMETER :: ONE    = 1.0E0_rk
    REAL(rk), PARAMETER :: HALF   = 0.5E0_rk
    REAL(rk), PARAMETER :: TWO    = 2.0E0_rk
    REAL(rk), PARAMETER :: SQRPI  = 5.6418958354775628695E-1_rk
    REAL(rk), PARAMETER :: THRESH = 0.46875E0_rk
    REAL(rk), PARAMETER :: SIXTEN = 16.0E0_rk

    !------------------------------------------------------------------
    !  Machine-dependent constants: IEEE single precision values
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: XINF   =  3.40E+38_r4
    REAL(rk), PARAMETER :: XNEG   = -9.382E0_r4
    REAL(rk), PARAMETER :: XSMALL =  5.96E-8_r4 
    REAL(rk), PARAMETER :: XBIG   =  9.194E0_r4
    REAL(rk), PARAMETER :: XHUGE  =  2.90E3_r4
    REAL(rk), PARAMETER :: XMAX   =  4.79E37_r4

    !------------------------------------------------------------------
    !  Coefficients for approximation to  erf  in first interval
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: A(5) = (/ 3.16112374387056560E00_rk, 1.13864154151050156E02_rk, &
         3.77485237685302021E02_rk, 3.20937758913846947E03_rk, &
         1.85777706184603153E-1_rk /)
    REAL(rk), PARAMETER :: B(4) = (/ 2.36012909523441209E01_rk, 2.44024637934444173E02_rk, &
         1.28261652607737228E03_rk, 2.84423683343917062E03_rk /)

    !------------------------------------------------------------------
    !  Coefficients for approximation to  erfc  in second interval
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: C(9) = (/ 5.64188496988670089E-1_rk, 8.88314979438837594E00_rk, &
         6.61191906371416295E01_rk, 2.98635138197400131E02_rk, &
         8.81952221241769090E02_rk, 1.71204761263407058E03_rk, &
         2.05107837782607147E03_rk, 1.23033935479799725E03_rk, &
         2.15311535474403846E-8_rk /)
    REAL(rk), PARAMETER :: D(8) = (/ 1.57449261107098347E01_rk, 1.17693950891312499E02_rk, &
         5.37181101862009858E02_rk, 1.62138957456669019E03_rk, &
         3.29079923573345963E03_rk, 4.36261909014324716E03_rk, &
         3.43936767414372164E03_rk, 1.23033935480374942E03_rk /)

    !------------------------------------------------------------------
    !  Coefficients for approximation to  erfc  in third interval
    !------------------------------------------------------------------
    REAL(rk), PARAMETER :: P(6) = (/ 3.05326634961232344E-1_rk, 3.60344899949804439E-1_rk, &
         1.25781726111229246E-1_rk, 1.60837851487422766E-2_rk, &
         6.58749161529837803E-4_rk, 1.63153871373020978E-2_rk /)
    REAL(rk), PARAMETER :: Q(5) = (/ 2.56852019228982242E00_rk, 1.87295284992346047E00_rk, &
         5.27905102951428412E-1_rk, 6.05183413124413191E-2_rk, &
         2.33520497626869185E-3_rk /)

    !------------------------------------------------------------------
    X = ARG
    Y = ABS(X)
    IF (Y .LE. THRESH) THEN
       !------------------------------------------------------------------
       !  Evaluate  erf  for  |X| <= 0.46875
       !------------------------------------------------------------------
       YSQ = ZERO
       IF (Y .GT. XSMALL) YSQ = Y * Y
       XNUM = A(5)*YSQ
       XDEN = YSQ
       DO I = 1, 3
          XNUM = (XNUM + A(I)) * YSQ
          XDEN = (XDEN + B(I)) * YSQ
       END DO
       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
       IF (JINT .NE. 0) RESULT = ONE - RESULT
       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
       GO TO 80
    ELSE IF (Y .LE. FOUR) THEN
       !------------------------------------------------------------------
       !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
       !------------------------------------------------------------------
       XNUM = C(9)*Y
       XDEN = Y
       DO I = 1, 7
          XNUM = (XNUM + C(I)) * Y
          XDEN = (XDEN + D(I)) * Y
       END DO
       RESULT = (XNUM + C(8)) / (XDEN + D(8))
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
    ELSE
       !------------------------------------------------------------------
       !  Evaluate  erfc  for |X| > 4.0
       !------------------------------------------------------------------
       RESULT = ZERO
       IF (Y .GE. XBIG) THEN
          IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 30
          IF (Y .GE. XHUGE) THEN
             RESULT = SQRPI / Y
             GO TO 30
          END IF
       END IF
       YSQ = ONE / (Y * Y)
       XNUM = P(6)*YSQ
       XDEN = YSQ
       DO I = 1, 4
          XNUM = (XNUM + P(I)) * YSQ
          XDEN = (XDEN + Q(I)) * YSQ
       END DO
       RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
       RESULT = (SQRPI -  RESULT) / Y
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
    END IF
30  CONTINUE
    !------------------------------------------------------------------
    !  Fix up for negative argument, erf, etc.
    !------------------------------------------------------------------
    IF (JINT .EQ. 0) THEN
       RESULT = (HALF - RESULT) + HALF
       IF (X .LT. ZERO) RESULT = -RESULT
    ELSE IF (JINT .EQ. 1) THEN
       IF (X .LT. ZERO) RESULT = TWO - RESULT
    ELSE
       IF (X .LT. ZERO) THEN
          IF (X .LT. XNEG) THEN
             RESULT = XINF
          ELSE
             YSQ = AINT(X*SIXTEN)/SIXTEN
             DEL = (X-YSQ)*(X+YSQ)
             Y = EXP(YSQ*YSQ) * EXP(DEL)
             RESULT = (Y+Y) - RESULT
          END IF
       END IF
    END IF
80  CONTINUE
  END SUBROUTINE CALERF_r4


  !
  ! Common block and statement functions for saturation vapor pressure
  ! look-up procedure, J. J. Hack, February 1990
  !
  ! $Id$
  !
  !module wv_saturation
  !  use shr_kind_mod, only: r8 => shr_kind_r8
  !  use abortutils,   only: endrun
  !  use cam_logfile,  only: iulog
  !
  !  implicit none
  !  private
  !  save
  !!
  !! Public interfaces
  !!
  !  public gestbl   ! Initialization subroutine
  !  public estblf   ! saturation pressure table lookup
  !  public aqsat    ! Returns saturation vapor pressure
  !  public aqsatd   ! Same as aqsat, but also returns a temperature derivitive
  !  public vqsatd   ! Vector version of aqsatd
  !  public fqsatd   ! Function version of vqsatd
  !  public qsat_water  ! saturation mixing ration with respect to liquid water
  !  public vqsat_water ! vector version of qsat_water
  !  public qsat_ice    ! saturation mixing ration with respect to ice
  !  public vqsat_ice   ! vector version of qsat_ice
  !  public vqsatd_water
  !  public aqsat_water
  !  public vqsatd2_water         ! Variant of vqsatd_water to print out dqsdT
  !  public vqsatd2_water_single  ! Single value version of vqsatd2_water
  !  public vqsatd2
  !  public vqsatd2_single
  !  public polysvp
  !
  ! Data used by cldwat
  !
  !  public hlatv, tmin, hlatf, rgasv, pcf, cp, epsqs, ttrice
  !
  ! Data
  !
  !  integer plenest  ! length of saturation vapor pressure table
  !  parameter (plenest=250)
  !
  ! Table of saturation vapor pressure values es from tmin degrees
  ! to tmax+1 degrees k in one degree increments.  ttrice defines the
  ! transition region where es is a combination of ice & water values
  !
  !  real(r8) estbl(plenest)      ! table values of saturation vapor pressure
  !  real(r8) tmin       ! min temperature (K) for table
  !  real(r8) tmax       ! max temperature (K) for table
  !  real(r8) ttrice     ! transition range from es over H2O to es over ice
  !  real(r8) pcf(6)     ! polynomial coeffs -> es transition water to ice
  !  real(r8) epsqs      ! Ratio of h2o to dry air molecular weights 
  !  real(r8) rgasv      ! Gas constant for water vapor
  !  real(r8) hlatf      ! Latent heat of vaporization
  !  real(r8) hlatv      ! Latent heat of fusion
  !  real(r8) cp         ! specific heat of dry air
  !  real(r8) tmelt      ! Melting point of water (K)
  !  logical icephs  ! false => saturation vapor press over water only

!!!CONTAINS

  REAL(r8) FUNCTION estblf( td )
    !
    ! Saturation vapor pressure table lookup
    !
    REAL(r8), INTENT(in) :: td         ! Temperature for saturation lookup
    !
    REAL(r8) :: e       ! intermediate variable for es look-up
    REAL(r8) :: ai
    INTEGER  :: i
    !
    e = MAX(MIN(td,tmax),tmin)   ! partial pressure
    i = INT(e-tmin)+1
    ai = AINT(e-tmin)
    estblf = (tmin+ai-e+1._r8)* &
         estbl(i)-(tmin+ai-e)* &
         estbl(i+1)
  END FUNCTION estblf

  SUBROUTINE gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   , &
       latvap  ,latice  ,rh2o    ,cpair   ,tmeltx   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Builds saturation vapor pressure table for later lookup procedure.
    ! 
    ! Method: 
    ! Uses Goff & Gratch (1946) relationships to generate the table
    ! according to a set of free parameters defined below.  Auxiliary
    ! routines are also included for making rapid estimates (well with 1%)
    ! of both es and d(es)/dt for the particular table configuration.
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use spmd_utils, only: masterproc
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: tmn           ! Minimum temperature entry in es lookup table
    REAL(r8), INTENT(in) :: tmx           ! Maximum temperature entry in es lookup table
    REAL(r8), INTENT(in) :: epsil         ! Ratio of h2o to dry air molecular weights
    REAL(r8), INTENT(in) :: trice         ! Transition range from es over range to es over ice
    REAL(r8), INTENT(in) :: latvap        ! Latent heat of vaporization
    REAL(r8), INTENT(in) :: latice        ! Latent heat of fusion
    REAL(r8), INTENT(in) :: rh2o          ! Gas constant for water vapor
    REAL(r8), INTENT(in) :: cpair         ! Specific heat of dry air
    REAL(r8), INTENT(in) :: tmeltx        ! Melting point of water (K)
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(r8) t             ! Temperature
    INTEGER n          ! Increment counter
    INTEGER lentbl     ! Calculated length of lookup table
    INTEGER itype      ! Ice phase: 0 -> no ice phase
    !            1 -> ice phase, no transition
    !           -x -> ice phase, x degree transition
    LOGICAL ip         ! Ice phase logical flag
    !
    !-----------------------------------------------------------------------
    !
    ! Set es table parameters
    !
    tmin   = tmn       ! Minimum temperature entry in table
    tmax   = tmx       ! Maximum temperature entry in table
    ttrice = trice     ! Trans. range from es over h2o to es over ice
    icephs = ip        ! Ice phase (true or false)
    !
    ! Set physical constants required for es calculation
    !
    epsqs  = epsil
    hlatv  = latvap
    hlatf  = latice
    rgasv  = rh2o
    cp     = cpair
    !tmelt  = tmeltx
    !
    lentbl = INT(tmax-tmin+2.000001_r8)
    IF (lentbl .GT. plenest) THEN
       WRITE(iulog,9000) tmax, tmin, plenest
       CALL endrun ('GESTBL')    ! Abnormal termination
    END IF
    !
    ! Begin building es table.
    ! Check whether ice phase requested.
    ! If so, set appropriate transition range for temperature
    !
    IF (icephs) THEN
       IF (ttrice /= 0.0_r8) THEN
          itype = -ttrice
       ELSE
          itype = 1
       END IF
    ELSE
       itype = 0
    END IF
    !
    t = tmin - 1.0_r8
    DO n=1,lentbl
       t = t + 1.0_r8
       CALL gffgch(t,estbl(n),itype)
    END DO
    !
    DO n=lentbl+1,plenest
       estbl(n) = -99999.0_r8
    END DO
    !
    ! Table complete -- Set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
    ! is valid in the range -40 < t < 0 (degrees C).
    !
    !                  --- Degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_r8
    pcf(2) = -5.47288442819e+00_r8
    pcf(3) = -3.67471858735e-01_r8
    pcf(4) = -8.95963532403e-03_r8
    pcf(5) = -7.78053686625e-05_r8
    !
    !                  --- Degree 6 approximation ---
    !
    !-----pcf(1) =  7.63285250063e-02
    !-----pcf(2) = -5.86048427932e+00
    !-----pcf(3) = -4.38660831780e-01
    !-----pcf(4) = -1.37898276415e-02
    !-----pcf(5) = -2.14444472424e-04
    !-----pcf(6) = -1.36639103771e-06
    !
    !if (masterproc) then
    !WRITE(iulog,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
    !end if

    RETURN
    !
9000 FORMAT('GESTBL: FATAL ERROR *********************************',/, &
         ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
         ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
         ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
    !
  END SUBROUTINE gestbl

  !
  !-----------------------------------------------------------------------

  SUBROUTINE gffgch(t       ,es      ,itype   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Computes saturation vapor pressure over water and/or over ice using
    ! Goff & Gratch (1946) relationships. 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! T (temperature), and itype are input parameters, while es (saturation
    ! vapor pressure) is an output parameter.  The input parameter itype
    ! serves two purposes: a value of zero indicates that saturation vapor
    ! pressures over water are to be returned (regardless of temperature),
    ! while a value of one indicates that saturation vapor pressures over
    ! ice should be returned when t is less than freezing degrees.  If itype
    ! is negative, its absolute value is interpreted to define a temperature
    ! transition region below freezing in which the returned
    ! saturation vapor pressure is a weighted average of the respective ice
    ! and water value.  That is, in the temperature range 0 => -itype
    ! degrees c, the saturation vapor pressures are assumed to be a weighted
    ! average of the vapor pressure over supercooled water and ice (all
    ! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use shr_kind_mod, only: r8 => shr_kind_r8
    !   use physconst, only: tmelt
    !   use abortutils, only: endrun

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: t          ! Temperature
    !
    ! Output arguments
    !
    INTEGER, INTENT(inout) :: itype   ! Flag for ice phase and associated transition

    REAL(r8), INTENT(out) :: es         ! Saturation vapor pressure
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(r8) e1         ! Intermediate scratch variable for es over water
    REAL(r8) e2         ! Intermediate scratch variable for es over water
    REAL(r8) eswtr      ! Saturation vapor pressure over water
    REAL(r8) f          ! Intermediate scratch variable for es over water
    REAL(r8) f1         ! Intermediate scratch variable for es over water
    REAL(r8) f2         ! Intermediate scratch variable for es over water
    REAL(r8) f3         ! Intermediate scratch variable for es over water
    REAL(r8) f4         ! Intermediate scratch variable for es over water
    REAL(r8) f5         ! Intermediate scratch variable for es over water
    REAL(r8) ps         ! Reference pressure (mb)
    REAL(r8) t0         ! Reference temperature (freezing point of water)
    REAL(r8) term1      ! Intermediate scratch variable for es over ice
    REAL(r8) term2      ! Intermediate scratch variable for es over ice
    REAL(r8) term3      ! Intermediate scratch variable for es over ice
    REAL(r8) tr         ! Transition range for es over water to es over ice
    REAL(r8) ts         ! Reference temperature (boiling point of water)
    REAL(r8) weight     ! Intermediate scratch variable for es transition
    INTEGER itypo   ! Intermediate scratch variable for holding itype
    !
    !-----------------------------------------------------------------------
    !
    ! Check on whether there is to be a transition region for es
    !
    IF (itype < 0) THEN
       tr    = ABS(float(itype))
       itypo = itype
       itype = 1
    ELSE
       tr    = 0.0_r8
       itypo = itype
    END IF
    IF (tr > 40.0_r8) THEN
       WRITE(6,900) tr
       CALL endrun ('GFFGCH')                ! Abnormal termination
    END IF
    !
    IF(t < (tmelt - tr) .AND. itype == 1) go to 10
    !
    ! Water
    !
    ps = 1013.246_r8
    ts = 373.16_r8
    e1 = 11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*LOG10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = LOG10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8
    eswtr = es
    !
    IF(t >= tmelt .OR. itype == 0) go to 20
    !
    ! Ice
    !
10  CONTINUE
    t0    = tmelt
    term1 = 2.01889049_r8/(t0/t)
    term2 = 3.56654_r8*LOG(t0/t)
    term3 = 20.947031_r8*(t0/t)
    es    = 575.185606e10_r8*EXP(-(term1 + term2 + term3))
    !
    IF (t < (tmelt - tr)) go to 20
    !
    ! Weighted transition between water and ice
    !
    weight = MIN((tmelt - t)/tr,1.0_r8)
    es = weight*es + (1.0_r8 - weight)*eswtr
    !
20  CONTINUE
    itype = itypo
    RETURN
    !
900 FORMAT('GFFGCH: FATAL ERROR ******************************',/, &
         'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
         ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
         ' 40.0 DEGREES C',/, ' TR = ',f7.2)
    !
  END SUBROUTINE gffgch

  SUBROUTINE aqsat(t       ,p       ,es      ,qs        ,ii      , &
       ILEN    ,kk      ,kstart  ,kend      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
    ! This routine is useful for evaluating only a selected region in the
    ! vertical.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ii             ! I dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: kk             ! K dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: ILEN           ! Length of vectors in I direction which
    INTEGER, INTENT(in) :: kstart         ! Starting location in K direction
    INTEGER, INTENT(in) :: kend           ! Ending location in K direction
    REAL(r8), INTENT(in) :: t(ii,kk)          ! Temperature
    REAL(r8), INTENT(in) :: p(ii,kk)          ! Pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(ii,kk)         ! Saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(ii,kk)         ! Saturation specific humidity
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(r8) omeps             ! 1 - 0.622
    INTEGER i, k           ! Indices
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative values
          ! that can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE aqsat

  !++xl
  SUBROUTINE aqsat_water(t       ,p       ,es      ,qs        ,ii      , &
       ILEN    ,kk      ,kstart  ,kend      )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
    ! This routine is useful for evaluating only a selected region in the
    ! vertical.
    !
    ! Method:
    ! <Describe the algorithm(s) used in the routine.>
    ! <Also include any applicable external references.>
    !
    ! Author: J. Hack
    !
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ii             ! I dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: kk             ! K dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: ILEN           ! Length of vectors in I direction which
    INTEGER, INTENT(in) :: kstart         ! Starting location in K direction
    INTEGER, INTENT(in) :: kend           ! Ending location in K direction
    REAL(r8), INTENT(in) :: t(ii,kk)          ! Temperature
    REAL(r8), INTENT(in) :: p(ii,kk)          ! Pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(ii,kk)         ! Saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(ii,kk)         ! Saturation specific humidity
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(r8) omeps             ! 1 - 0.622
    INTEGER i, k           ! Indices
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          !        es(i,k) = estblf(t(i,k))
          es(i,k) = polysvp(t(i,k),0)
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative values
          ! that can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE aqsat_water
  !--xl


  SUBROUTINE aqsatd(t       ,p       ,es      ,qs      ,gam     , &
       ii      ,ILEN    ,kk      ,kstart  ,kend    )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g).   
    ! 
    ! Method: 
    ! Differs from aqsat by also calculating and returning
    ! gamma (l/cp)*(d(qsat)/dT)
    ! Input arrays temperature and pressure (dimensioned ii,kk).
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ii            ! I dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: ILEN          ! Vector length in I direction
    INTEGER, INTENT(in) :: kk            ! K dimension of arrays t, p, es, qs
    INTEGER, INTENT(in) :: kstart        ! Starting location in K direction
    INTEGER, INTENT(in) :: kend          ! Ending location in K direction

    REAL(r8), INTENT(in) :: t(ii,kk)         ! Temperature
    REAL(r8), INTENT(in) :: p(ii,kk)         ! Pressure

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(ii,kk)        ! Saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(ii,kk)        ! Saturation specific humidity
    REAL(r8), INTENT(out) :: gam(ii,kk)       ! (l/cp)*(d(qs)/dt)
    !
    !---------------------------Local workspace-----------------------------
    !
    LOGICAL lflg          ! True if in temperature transition region
    INTEGER i             ! i index for vector calculations
    INTEGER k             ! k index
    REAL(r8) omeps            ! 1. - 0.622
    REAL(r8) trinv            ! Reciprocal of ttrice (transition range)
    REAL(r8) tc               ! Temperature (in degrees C)
    REAL(r8) weight           ! Weight for es transition from water to ice
    REAL(r8) hltalt           ! Appropriately modified hlat for T derivatives
    REAL(r8) hlatsb           ! hlat weighted in transition region
    REAL(r8) hlatvp           ! hlat modified for t changes above freezing
    REAL(r8) tterm            ! Account for d(es)/dT in transition region
    REAL(r8) desdt            ! d(es)/dT
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative qs
          ! values which can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    !
    DO k=kstart,kend
       DO i=1,ILEN
          !
          ! Weighting of hlat accounts for transition from water to ice
          ! polynomial expression approximates difference between es over
          ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
          ! -40): required for accurate estimate of es derivative in transition
          ! range from ice to water also accounting for change of hlatv with t
          ! above freezing where constant slope is given by -2369 j/(kg c) =cpv - cw
          !
          tc     = t(i,k) - tmelt
          lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
          weight = MIN(-tc*trinv,1.0_r8)
          hlatsb = hlatv + weight*hlatf
          hlatvp = hlatv - 2369.0_r8*tc
          IF (t(i,k) < tmelt) THEN
             hltalt = hlatsb
          ELSE
             hltalt = hlatvp
          END IF
          IF (lflg) THEN
             tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
          ELSE
             tterm = 0.0_r8
          END IF
          desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
          gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
          IF (qs(i,k) == 1.0_r8) gam(i,k) = 0.0_r8
       END DO
    END DO
    !
    go to 20
    !
    ! No icephs or water to ice transition
    !
10  DO k=kstart,kend
       DO i=1,ILEN
          !
          ! Account for change of hlatv with t above freezing where
          ! constant slope is given by -2369 j/(kg c) = cpv - cw
          !
          hlatvp = hlatv - 2369.0_r8*(t(i,k)-tmelt)
          IF (icephs) THEN
             hlatsb = hlatv + hlatf
          ELSE
             hlatsb = hlatv
          END IF
          IF (t(i,k) < tmelt) THEN
             hltalt = hlatsb
          ELSE
             hltalt = hlatvp
          END IF
          desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
          gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
          IF (qs(i,k) == 1.0_r8) gam(i,k) = 0.0_r8
       END DO
    END DO
    !
20  RETURN
  END SUBROUTINE aqsatd

  SUBROUTINE vqsatd(t       ,p       ,es      ,qs      ,gam      , &
       len     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg   ! true if in temperature transition region
    !
    INTEGER i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(r8) tc        ! temperature (in degrees C)
    REAL(r8) weight    ! weight for es transition from water to ice
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) tterm     ! account for d(es)/dT in transition region
    REAL(r8) desdt     ! d(es)/dT
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = estblf(t(i))
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    DO i=1,len
       !
       ! Weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition
       ! range from ice to water also accounting for change of hlatv with t
       ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
       weight = MIN(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       IF (lflg) THEN
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       ELSE
          tterm = 0.0_r8
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    RETURN
    !
    ! No icephs or water to ice transition
    !
10  DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       IF (icephs) THEN
          hlatsb = hlatv + hlatf
       ELSE
          hlatsb = hlatv
       END IF
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd

  !++xl
  SUBROUTINE vqsatd_water(t       ,p       ,es      ,qs      ,gam      , &
       len     )

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)

    !
    !--------------------------Local Variables------------------------------
    !
    !
    INTEGER i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) desdt     ! d(es)/dT
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = polysvp(t(i),0)
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! No icephs or water to ice transition
    !
    DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       hlatsb = hlatv
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd_water

  FUNCTION polysvp (T,TYPE)
    !  Compute saturation vapor pressure by using
    ! function from Goff and Gatch (1946)

    !  Polysvp returned in units of pa.
    !  T is input in units of K.
    !  type refers to saturation with respect to liquid (0) or ice (1)

    !REAL(r8) dum

    REAL(r8) T,polysvp

    INTEGER TYPE

    ! ice

    IF (TYPE.EQ.1) THEN

       ! Goff Gatch equation (good down to -100 C)

       polysvp = 10._r8**(-9.09718_r8*(273.16_r8/t-1._r8)-3.56654_r8* &
            LOG10(273.16_r8/t)+0.876793_r8*(1._r8-t/273.16_r8)+ &
            LOG10(6.1071_r8))*100._r8

    END IF

    ! Goff Gatch equation, uncertain below -70 C

    IF (TYPE.EQ.0) THEN
       polysvp = 10._r8**(-7.90298_r8*(373.16_r8/t-1._r8)+ &
            5.02808_r8*LOG10(373.16_r8/t)- &
            1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/373.16_r8))-1._r8)+ &
            8.1328e-3_r8*(10._r8**(-3.49149_r8*(373.16_r8/t-1._r8))-1._r8)+ &
            LOG10(1013.246_r8))*100._r8
    END IF


  END FUNCTION polysvp
  !--xl

  INTEGER FUNCTION fqsatd(t    ,p    ,es    ,qs   ,gam   , len     )
    !----------------------------------------------------------------------- 
    ! Purpose: 
    ! This is merely a function interface vqsatd.
    !------------------------------Arguments--------------------------------
    ! Input arguments
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure
    ! Output arguments
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    ! Call vqsatd
    CALL vqsatd(t       ,p       ,es      ,qs      ,gam  , len     )
    fqsatd = 1
    RETURN
  END FUNCTION fqsatd

  REAL(r8) FUNCTION qsat_water(t,p)
    !  saturation mixing ratio w/respect to liquid water
    REAL(r8) t ! temperature
    REAL(r8) p ! pressure (Pa)
    REAL(r8) es ! saturation vapor pressure (Pa)
    REAL(r8) ps, ts, e1, e2, f1, f2, f3, f4, f5, f
    !  real(r8) t0inv ! 1/273.
    !  data t0inv/0.003663/
    !  save t0inv
    !  es = 611.*exp(hlatv/rgasv*(t0inv-1./t))

    ps = 1013.246_r8
    ts = 373.16_r8
    e1 = 11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*LOG10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = LOG10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8

    qsat_water = epsqs*es/(p-(1.-epsqs)*es) ! saturation w/respect to liquid only
    IF(qsat_water < 0.) qsat_water = 1.

  END FUNCTION qsat_water

  SUBROUTINE vqsat_water(t,p,qsat_water,len)
    !  saturation mixing ratio w/respect to liquid water
    INTEGER, INTENT(in)  :: len
    REAL(r8) t(len) ! temperature
    REAL(r8) p(len) ! pressure (Pa)
    REAL(r8) qsat_water(len)
    REAL(r8) es ! saturation vapor pressure (Pa)
    REAL(r8), PARAMETER :: t0inv = 1._r8/273._r8
    REAL(r8) coef
    INTEGER :: i

    coef = hlatv/rgasv
    DO i=1,len
       es = 611._r8*EXP(coef*(t0inv-1./t(i)))
       qsat_water(i) = epsqs*es/(p(i)-(1.-epsqs)*es) ! saturation w/respect to liquid only
       IF(qsat_water(i) < 0.) qsat_water(i) = 1.
    ENDDO

    RETURN

  END SUBROUTINE vqsat_water

  REAL(r8) FUNCTION qsat_ice(t,p)
    !  saturation mixing ratio w/respect to ice
    REAL(r8) t ! temperature
    REAL(r8) p ! pressure (Pa)
    REAL(r8) es ! saturation vapor pressure (Pa)
    REAL(r8), PARAMETER :: t0inv = 1._r8/273._r8
    es = 611.*EXP((hlatv+hlatf)/rgasv*(t0inv-1./t))
    qsat_ice = epsqs*es/(p-(1.-epsqs)*es) ! saturation w/respect to liquid only
    IF(qsat_ice < 0.) qsat_ice = 1.

  END FUNCTION qsat_ice

  SUBROUTINE vqsat_ice(t,p,qsat_ice,len)
    !  saturation mixing ratio w/respect to liquid water
    INTEGER,INTENT(in) :: len
    REAL(r8) t(len) ! temperature
    REAL(r8) p(len) ! pressure (Pa)
    REAL(r8) qsat_ice(len)
    REAL(r8) es ! saturation vapor pressure (Pa)
    REAL(r8), PARAMETER :: t0inv = 1._r8/273._r8
    REAL(r8) coef
    INTEGER :: i

    coef = (hlatv+hlatf)/rgasv
    DO i=1,len
       es = 611.*EXP(coef*(t0inv-1./t(i)))
       qsat_ice(i) = epsqs*es/(p(i)-(1.-epsqs)*es) ! saturation w/respect to liquid only
       IF(qsat_ice(i) < 0.) qsat_ice(i) = 1.
    ENDDO

    RETURN

  END SUBROUTINE vqsat_ice

  ! Sungsu
  ! Below two subroutines (vqsatd2_water,vqsatd2_water_single) are by Sungsu
  ! Replace 'gam -> dqsdt'
  ! Sungsu

  SUBROUTINE vqsatd2_water(t       ,p       ,es      ,qs      ,dqsdt      , &
       len     )

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    ! real(r8), intent(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    ! Sungsu
    REAL(r8), INTENT(out) :: dqsdt(len)  ! (d(qs)/dt)
    ! End by Sungsu

    !
    !--------------------------Local Variables------------------------------
    !
    !
    INTEGER i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) desdt     ! d(es)/dT

    ! Sungsu
    REAL(r8) gam(len)  ! (l/cp)*(d(qs)/dt)
    ! End by Sungsu

    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = polysvp(t(i),0)
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! No icephs or water to ice transition
    !
    DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       hlatsb = hlatv
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
       ! Sungsu
       dqsdt(i) = (cp/hltalt)*gam(i)
       ! End by Sungsu
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd2_water

  SUBROUTINE vqsatd2_water_single(t       ,p       ,es      ,qs      ,dqsdt)

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !

    REAL(r8), INTENT(in) :: t       ! temperature
    REAL(r8), INTENT(in) :: p       ! pressure

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs   ! saturation specific humidity
    ! real(r8), intent(out) :: gam  ! (l/cp)*(d(qs)/dt)
    ! Sungsu
    REAL(r8), INTENT(out) :: dqsdt  ! (d(qs)/dt)
    ! End by Sungsu

    !
    !--------------------------Local Variables------------------------------
    !
    !
    !INTEGER i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) desdt     ! d(es)/dT

    ! Sungsu
    REAL(r8) gam  ! (l/cp)*(d(qs)/dt)
    ! End by Sungsu

    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    !  do i=1,len
    es = polysvp(t,0)
    !
    ! Saturation specific humidity
    !
    qs = epsqs*es/(p - omeps*es)
    !
    ! The following check is to avoid the generation of negative
    ! values that can occur in the upper stratosphere and mesosphere
    !
    qs = MIN(1.0_r8,qs)
    !
    IF (qs < 0.0_r8) THEN
       qs = 1.0_r8
       es = p
    END IF
    !  end do
    !
    ! No icephs or water to ice transition
    !
    !  do i=1,len
    !
    ! Account for change of hlatv with t above freezing where
    ! constant slope is given by -2369 j/(kg c) = cpv - cw
    !
    hlatvp = hlatv - 2369.0_r8*(t-tmelt)
    hlatsb = hlatv
    IF (t < tmelt) THEN
       hltalt = hlatsb
    ELSE
       hltalt = hlatvp
    END IF
    desdt  = hltalt*es/(rgasv*t*t)
    gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
    IF (qs == 1.0_r8) gam = 0.0_r8
    ! Sungsu
    dqsdt = (cp/hltalt)*gam
    ! End by Sungsu
    !  end do
    !
    RETURN
    !
  END SUBROUTINE vqsatd2_water_single


  SUBROUTINE vqsatd2(t       ,p       ,es      ,qs      ,dqsdt      , &
       len     )
    !----------------------------------------------------------------------- 
    ! Sungsu : This is directly copied from 'vqsatd' but 'dqsdt' is output
    !          instead of gam for use in Sungsu's equilibrium stratiform
    !          macrophysics scheme.
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(r8), INTENT(in) :: t(len)       ! temperature
    REAL(r8), INTENT(in) :: p(len)       ! pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    ! real(r8), intent(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    ! Sungsu
    REAL(r8), INTENT(out) :: dqsdt(len)  ! (d(qs)/dt)
    ! End by Sungsu 

    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg   ! true if in temperature transition region
    !
    INTEGER i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(r8) tc        ! temperature (in degrees C)
    REAL(r8) weight    ! weight for es transition from water to ice
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) tterm     ! account for d(es)/dT in transition region
    REAL(r8) desdt     ! d(es)/dT

    ! Sungsu
    REAL(r8) gam(len)  ! (l/cp)*(d(qs)/dt)
    ! End by Sungsu
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = estblf(t(i))
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    DO i=1,len
       !
       ! Weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition
       ! range from ice to water also accounting for change of hlatv with t
       ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
       weight = MIN(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       IF (lflg) THEN
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       ELSE
          tterm = 0.0_r8
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
       ! Sungsu
       dqsdt(i) = (cp/hltalt)*gam(i)
       ! End by Sungsu
    END DO
    RETURN
    !
    ! No icephs or water to ice transition
    !
10  DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       IF (icephs) THEN
          hlatsb = hlatv + hlatf
       ELSE
          hlatsb = hlatv
       END IF
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
       ! Sungsu
       dqsdt(i) = (cp/hltalt)*gam(i)
       ! End by Sungsu
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd2


  ! Below routine is by Sungsu

  SUBROUTINE vqsatd2_single(t       ,p       ,es      ,qs      ,dqsdt)
    !----------------------------------------------------------------------- 
    ! Sungsu : This is directly copied from 'vqsatd' but 'dqsdt' is output
    !          instead of gam for use in Sungsu's equilibrium stratiform
    !          macrophysics scheme.
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(r8), INTENT(in) :: t       ! temperature
    REAL(r8), INTENT(in) :: p       ! pressure
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: es     ! saturation vapor pressure
    REAL(r8), INTENT(out) :: qs     ! saturation specific humidity
    ! real(r8), intent(out) :: gam    ! (l/cp)*(d(qs)/dt)
    ! Sungsu
    REAL(r8), INTENT(out) :: dqsdt  ! (d(qs)/dt)
    ! End by Sungsu 

    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg   ! true if in temperature transition region
    !
    !  integer i      ! index for vector calculations
    !
    REAL(r8) omeps     ! 1. - 0.622
    REAL(r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(r8) tc        ! temperature (in degrees C)
    REAL(r8) weight    ! weight for es transition from water to ice
    REAL(r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(r8) hlatsb    ! hlat weighted in transition region
    REAL(r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(r8) tterm     ! account for d(es)/dT in transition region
    REAL(r8) desdt     ! d(es)/dT

    ! Sungsu
    REAL(r8) gam       ! (l/cp)*(d(qs)/dt)
    ! End by Sungsu
    !
    !-----------------------------------------------------------------------
    !
    omeps = 1.0_r8 - epsqs

    !  do i=1,len

    es = estblf(t)
    !
    ! Saturation specific humidity
    !
    qs = epsqs*es/(p - omeps*es)
    !
    ! The following check is to avoid the generation of negative
    ! values that can occur in the upper stratosphere and mesosphere
    !
    qs = MIN(1.0_r8,qs)
    !
    IF (qs < 0.0_r8) THEN
       qs = 1.0_r8
       es = p
    END IF

    !  end do
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice

    !  do i=1,len
    !
    ! Weighting of hlat accounts for transition from water to ice
    ! polynomial expression approximates difference between es over
    ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
    ! -40): required for accurate estimate of es derivative in transition
    ! range from ice to water also accounting for change of hlatv with t
    ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
    !
    tc     = t - tmelt
    lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
    weight = MIN(-tc*trinv,1.0_r8)
    hlatsb = hlatv + weight*hlatf
    hlatvp = hlatv - 2369.0_r8*tc
    IF (t < tmelt) THEN
       hltalt = hlatsb
    ELSE
       hltalt = hlatvp
    END IF
    IF (lflg) THEN
       tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
    ELSE
       tterm = 0.0_r8
    END IF
    desdt  = hltalt*es/(rgasv*t*t) + tterm*trinv
    gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
    IF (qs == 1.0_r8) gam = 0.0_r8
    ! Sungsu
    dqsdt = (cp/hltalt)*gam
    ! End by Sungsu
    !  end do
    RETURN
    !
    ! No icephs or water to ice transition
    !

10  CONTINUE

    !10 do i=1,len
    !
    ! Account for change of hlatv with t above freezing where
    ! constant slope is given by -2369 j/(kg c) = cpv - cw
    !
    hlatvp = hlatv - 2369.0_r8*(t-tmelt)
    IF (icephs) THEN
       hlatsb = hlatv + hlatf
    ELSE
       hlatsb = hlatv
    END IF
    IF (t < tmelt) THEN
       hltalt = hlatsb
    ELSE
       hltalt = hlatvp
    END IF
    desdt  = hltalt*es/(rgasv*t*t)
    gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
    IF (qs == 1.0_r8) gam = 0.0_r8
    ! Sungsu
    dqsdt = (cp/hltalt)*gam
    ! End by Sungsu

    !  end do
    !
    RETURN
    !
  END SUBROUTINE vqsatd2_single


  !end module wv_saturation 

  SUBROUTINE endrun (msg)
    !-------------------------------------------------------------------------------------------
    ! Purpose: Print an optional string and abort
    !-------------------------------------------------------------------------------------------
    !
    ! Arguments
    !
    CHARACTER(len=*), INTENT(in), OPTIONAL :: msg    ! string to be printed

    IF (PRESENT (msg)) THEN
       WRITE(6,*)'ENDRUN:', msg
    ELSE
       WRITE(6,*)'ENDRUN: called without a message string'
    END IF

    ! tcraig - abort causes a core dump, use stop instead
    !      call abort ()

    STOP

  END SUBROUTINE endrun

END MODULE uwshcu

