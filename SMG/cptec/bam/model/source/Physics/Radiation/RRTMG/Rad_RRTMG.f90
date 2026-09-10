
!  Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) (coarse mode) included
!  Fine aerosol mode is included
!  Modifications (15) are marked by
!  !tar begin.....!tar end

MODULE Rad_RRTMG

  USE radsw, ONLY : rad_rrtmg_sw, radsw_init 
  USE radlw, ONLY : rad_rrtmg_lw, radlw_init 

    IMPLICIT NONE
  SAVE

  PRIVATE
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers
  PUBLIC :: Init_Rad_RRTMG
  PUBLIC :: Run_Rad_RRTMG_SW
  PUBLIC :: Run_Rad_RRTMG_LW
  PUBLIC :: Finalize_Rad_RRTMG
  REAL(r8),PARAMETER :: pref=101300.0_r8  ! reference pressure at layer midpoints (Pa)
  INTEGER            :: rrtmg_levs=-1               ! number of pressure levels greate than 1.e-4_r8 mbar
  INTEGER ::    nbndsw=14
  INTEGER ::    nbndlw=16
  ! Molecular weights
  REAL(r8), PARAMETER :: mwdry       = 28.966_R8       ! molecular weight dry air ~ kg/kmole ~ g/mol
  REAL(r8), PARAMETER :: mwh2o       = 18.016_R8       ! molecular weight water vapor
  REAL(r8), PARAMETER :: mwco2       =  44.0_r8        ! molecular weight co2
  REAL(r8), PARAMETER :: mwn2o       =  44.0_r8        ! molecular weight n2o
  REAL(r8), PARAMETER :: mwch4       =  16.0_r8        ! molecular weight ch4
  REAL(r8), PARAMETER :: mwf11       = 136.0_r8        ! molecular weight cfc11
  REAL(r8), PARAMETER :: mwf12       = 120.0_r8        ! molecular weight cfc12
  REAL(r8), PARAMETER :: mwo3        =  48.0_r8        ! molecular weight O3
  REAL(r8), PARAMETER :: mwso2       =  64.0_r8        ! molecular weight SO2
  REAL(r8), PARAMETER :: mwso4       =  96.0_r8        ! molecular weight SO4
  REAL(r8), PARAMETER :: mwh2o2      =  34.0_r8        ! molecular weight H2O2
  REAL(r8), PARAMETER :: mwdms       =  62.0_r8        ! molecular weight DMS

  REAL(r8), PARAMETER :: amdw  = 1.607793_r8! Molecular weight of dry air / water vapor
  REAL(r8), PARAMETER :: amdc  = 0.658114_r8! Molecular weight of dry air / carbon dioxide
  REAL(r8), PARAMETER :: amdo  = 0.603428_r8! Molecular weight of dry air / ozone
  REAL(r8), PARAMETER :: amdm  = 1.805423_r8! Molecular weight of dry air / methane
  REAL(r8), PARAMETER :: amdn  = 0.658090_r8! Molecular weight of dry air / nitrous oxide
  REAL(r8), PARAMETER :: amdo2 = 0.905140_r8! Molecular weight of dry air / oxygen
  REAL(r8), PARAMETER :: amdc1 = 0.210852_r8! Molecular weight of dry air / CFC11
  REAL(r8), PARAMETER :: amdc2 = 0.239546_r8! Molecular weight of dry air / CFC12
  REAL(R8), PARAMETER :: SHR_CONST_STEBOL  = 5.67e-8_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
  REAL(KIND=r8), ALLOCATABLE :: a_h(:)
  REAL(KIND=r8), ALLOCATABLE :: b_h(:)
  INTEGER :: imca_sw
  INTEGER :: imca_lw
  CHARACTER(LEN=5) :: FeedBackOptics_cld='     '

CONTAINS

  !--------------------------------------------------------------------------------
  ! sets the number of model levels RRTMG operates
  !--------------------------------------------------------------------------------

  SUBROUTINE Init_Rad_RRTMG(a_hybr,b_hybr,kMax,nbndsw_in,nbndlw_in)
    IMPLICIT NONE
    INTEGER , INTENT(in) :: kMax
    REAL(KIND=r8)   , INTENT(in) :: a_hybr(kMax+1)
    REAL(KIND=r8)   , INTENT(in) :: b_hybr(kMax+1)
    INTEGER , INTENT(in) :: nbndsw_in
    INTEGER , INTENT(in) :: nbndlw_in
    REAL(r8) :: pref_mid(kMax)  ! reference pressure at layer midpoints (Pa)
    INTEGER  :: k
    nbndsw =nbndsw_in
    nbndlw =nbndlw_in
    ALLOCATE(a_h(kMax+1))
    ALLOCATE(b_h(kMax+1))
    a_h = a_hybr / 100._r8   !  in  mb !
    b_h = b_hybr
    imca_sw=0
    imca_lw=1
    !FeedBackOptics_cld='WRF'
    FeedBackOptics_cld='CAM5'

    DO k=1,kMax
!      pref_mid(k)=sl(kMax-k+1)*pref
!      SB  changed to hybrid (already top to bottom)
       pref_mid(k)=0.5_r8*(a_hybr(k)+a_hybr(k+1)+pref * &
                           (b_hybr(k)+b_hybr(k+1)))
    END DO
    !0.01 mb -- > 1.0 Pa
    rrtmg_levs = 0!count( pref_mid(:) > 1.e-2_r8 ) ! pascals (1e-4 mbar)
    DO k=kMax,1,-1
       IF (pref_mid(k) > 1.e-2_r8)rrtmg_levs=rrtmg_levs+1
       !PRINT*,'rrtmg_levs=',rrtmg_levs
    END DO
    IF(rrtmg_levs == 0)rrtmg_levs=kMax
    rrtmg_levs=kMax
    CALL radsw_init()
    CALL radlw_init(pref_mid,kMax)
  END SUBROUTINE Init_Rad_RRTMG


  !--------------------------------------------------------------------------------
  ! creates (alloacates) an rrtmg_state object
  !--------------------------------------------------------------------------------

  SUBROUTINE rrtmg_state_create(nCols,kMax,gps,gt,tg,pmidmb,pintmb,tlay  ,tlev  )
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: gps    (nCols)!mb
    REAL(KIND=r8), INTENT(IN   ) :: gt     (nCols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: tg     (nCols)

    REAL(KIND=r8), INTENT(OUT  ) :: pmidmb (nCols,rrtmg_levs)!mb
    REAL(KIND=r8), INTENT(OUT  ) :: pintmb (nCols,rrtmg_levs+1)!mb
    REAL(KIND=r8), INTENT(OUT  ) :: tlay   (nCols,rrtmg_levs)
    REAL(KIND=r8), INTENT(OUT  ) :: tlev   (nCols,rrtmg_levs+1)

    !
    !  LOCAL 
    !
    REAL(KIND=r8) :: tint   (nCols,kMax+1)    ! Model interface temperature
    REAL(KIND=r8) :: pint   (nCols,kMax+1)   ! Interface pressures  
    REAL(KIND=r8) :: lnpint (nCols,kMax+1)
    REAL(KIND=r8) :: pmid   (nCols,kMax)
    REAL(KIND=r8) :: lnpmid (nCols,kMax)

    REAL(r8) :: dy
    INTEGER  :: k
    INTEGER  :: kk
    INTEGER  :: i
    pmidmb=0.0_r8;pintmb=0.0_r8;tlay  =0.0_r8;tlev  =0.0_r8;tint  =0.0_r8
    pint  =0.0_r8;lnpint=0.0_r8;pmid  =0.0_r8;lnpmid=0.0_r8
    DO k=1,kMax+1,1
       DO i=1,nCols
!         pint    (i,k) = MAX(si((kMax+1)-k+1)*gps(i) ,1.0e-12_r8)
          pint    (i,k) = MAX(a_h(k)+b_h(k)*gps(i) ,1.0e-12_r8)
          lnpint  (i,k) =  LOG(pint  (i,k))
       END DO
    END DO
    DO k=1,kMax,1
       DO i=1,nCols
!         pmid    (i,k) = MAX(sl(kMax-k+1)*gps(i) ,1.0e-12_r8)
          pmid    (i,k) = MAX(0.5_r8*(a_h(k)+a_h(k+1)+(b_h(k)+b_h(k+1))* &
                                                       gps(i)),1.0e-12_r8)
          lnpmid  (i,k) = LOG(pmid(i,k))
       END DO
    END DO



    ! Calculate interface temperatures (following method
    ! used in radtpl for the longwave), using surface upward flux and
    ! stebol constant in mks units

    DO i = 1,nCols
       tint(i,1) = gt(i,1)
       tint(i,kMax+1) = tg(i)!sqrt(sqrt(lwup(i)/SHR_CONST_STEBOL))
       DO k = 2,kMax
          dy = (lnpint(i,k) - lnpmid(i,k)) / (lnpmid(i,k-1) - lnpmid(i,k))
          tint(i,k) = gt(i,k) - dy * (gt(i,k) - gt(i,k-1))
       END DO
    END DO

    DO k = 1, rrtmg_levs

       kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
       DO i = 1,nCols

          pmidmb(i,k) = pmid(i,kk) ! mbar
          pintmb(i,k) = pint(i,kk) ! mbar

          tlay(i,k)   = gt  (i,kk)
          tlev(i,k)   = tint(i,kk)
       END DO
    ENDDO
    DO i = 1,nCols

       ! bottom interface
       pintmb(i,rrtmg_levs+1) = pint(i,kMax+1)  ! mbar
       tlev  (i,rrtmg_levs+1) = tint(i,kMax+1)

       ! top layer thickness
       IF (rrtmg_levs==kMax+1) THEN
          pmidmb(i,1) = 0.5_r8 * pintmb(i,2) 
          pintmb(i,1) = 1.e-4_r8 ! mbar
       ENDIF
    ENDDO

  END SUBROUTINE rrtmg_state_create

  !--------------------------------------------------------------------------------
  ! updates the concentration fields
  !--------------------------------------------------------------------------------
  SUBROUTINE rrtmg_state_update(&
       nCols   , &! INTEGER , INTENT(IN   ) :: nCols
       kMax    , &! INTEGER , INTENT(IN   ) :: kMax
       sp_hum  , &! real(r8), INTENT(IN   ) :: sp_hum (nCols,kMax) ! specific  humidity   !( kg/kg).
       n2o     , &! real(r8), INTENT(IN   ) :: n2o    (nCols,kMax) ! nitrous oxide mass mixing ratio  !( mol/mol).
       ch4     , &! real(r8), INTENT(IN   ) :: ch4    (nCols,kMax) ! methane  mass mixing ratio !( mol/mol).
       o2      , &! real(r8), INTENT(IN   ) :: o2     (nCols,kMax) ! O2       mass mixing ratio !( mol/mol).
       cfc11   , &! real(r8), INTENT(IN   ) :: cfc11  (nCols,kMax) ! cfc11  mass mixing ratio !( mol/mol).
       cfc12   , &! real(r8), INTENT(IN   ) :: cfc12  (nCols,kMax) ! cfc12  mass mixing ratio !( mol/mol).
       o3      , &! real(r8), INTENT(IN   ) :: o3     (nCols,kMax) ! Ozone  mass mixing ratio !( kg/kg).
       co2     , &! real(r8), INTENT(IN   ) :: co2    (nCols,kMax) ! co2   mass mixing ratio !( mol/mol).
       
       h2ovmr  , &! real(r8), INTENT(OUT  ) :: h2ovmr  (nCols,rrtmg_levs) 
       o3vmr   , &! real(r8), INTENT(OUT  ) :: o3vmr   (nCols,rrtmg_levs) 
       co2vmr  , &! real(r8), INTENT(OUT  ) :: co2vmr  (nCols,rrtmg_levs) 
       ch4vmr  , &! real(r8), INTENT(OUT  ) :: ch4vmr  (nCols,rrtmg_levs) 
       o2vmr   , &! real(r8), INTENT(OUT  ) :: o2vmr   (nCols,rrtmg_levs) 
       n2ovmr  , &! real(r8), INTENT(OUT  ) :: n2ovmr  (nCols,rrtmg_levs) 
       cfc11vmr, &! real(r8), INTENT(OUT  ) :: cfc11vmr(nCols,rrtmg_levs) 
       cfc12vmr, &! real(r8), INTENT(OUT  ) :: cfc12vmr(nCols,rrtmg_levs) 
       cfc22vmr, &! real(r8), INTENT(OUT  ) :: cfc22vmr(nCols,rrtmg_levs) 
       ccl4vmr   )! real(r8), INTENT(OUT  ) :: ccl4vmr (nCols,rrtmg_levs) 

    IMPLICIT NONE

    INTEGER , INTENT(IN   ) :: nCols
    INTEGER , INTENT(IN   ) :: kMax
    REAL(r8), INTENT(IN   ) :: sp_hum  (nCols,kMax)        ! specific  humidity   !( kg/kg).
    REAL(r8), INTENT(IN   ) :: n2o     (nCols,kMax)        ! nitrous oxide mass mixing ratio  !( mol/mol).
    REAL(r8), INTENT(IN   ) :: ch4     (nCols,kMax)        ! methane  mass mixing ratio !( mol/mol).
    REAL(r8), INTENT(IN   ) :: o2      (nCols,kMax)        ! O2       mass mixing ratio !( mol/mol).
    REAL(r8), INTENT(IN   ) :: cfc11   (nCols,kMax)        ! cfc11  mass mixing ratio !( mol/mol).
    REAL(r8), INTENT(IN   ) :: cfc12   (nCols,kMax)        ! cfc12  mass mixing ratio !( mol/mol).
    REAL(r8), INTENT(IN   ) :: o3      (nCols,kMax)        ! Ozone  mass mixing ratio !( kg/kg).
    REAL(r8), INTENT(IN   ) :: co2     (nCols,kMax)        ! co2   mass mixing ratio !( mol/mol).

    REAL(r8), INTENT(OUT  ) :: h2ovmr  (nCols,rrtmg_levs)  ! Molecular weight of dry air / water vapor
    REAL(r8), INTENT(OUT  ) :: o3vmr   (nCols,rrtmg_levs)  ! Molecular weight of dry air / ozone
    REAL(r8), INTENT(OUT  ) :: co2vmr  (nCols,rrtmg_levs)  ! Molecular weight of dry air / carbon dioxide 
    REAL(r8), INTENT(OUT  ) :: ch4vmr  (nCols,rrtmg_levs)  ! Molecular weight of dry air / methane
    REAL(r8), INTENT(OUT  ) :: o2vmr   (nCols,rrtmg_levs)  ! Molecular weight of dry air / oxygen
    REAL(r8), INTENT(OUT  ) :: n2ovmr  (nCols,rrtmg_levs)  ! Molecular weight of dry air / nitrous oxide
    REAL(r8), INTENT(OUT  ) :: cfc11vmr(nCols,rrtmg_levs)  ! Molecular weight of dry air / CFC11
    REAL(r8), INTENT(OUT  ) :: cfc12vmr(nCols,rrtmg_levs)  ! Molecular weight of dry air / CFC12
    REAL(r8), INTENT(OUT  ) :: cfc22vmr(nCols,rrtmg_levs) 
    REAL(r8), INTENT(OUT  ) :: ccl4vmr (nCols,rrtmg_levs) 


    INTEGER  :: i, kk, k


    ! Get specific humidity         call rad_cnst_get_gas(icall,'H2O', pstate, pbuf, sp_hum)
    ! Get oxygen mass mixing ratio. call rad_cnst_get_gas(icall,'O2',  pstate, pbuf, o2)
    ! Get ozone mass mixing ratio.  call rad_cnst_get_gas(icall,'O3',  pstate, pbuf, o3)
    ! Get CO2 mass mixing ratio     call rad_cnst_get_gas(icall,'CO2', pstate, pbuf, co2)
    ! Get N2O mass mixing ratio     call rad_cnst_get_gas(icall,'N2O', pstate, pbuf, n2o)
    ! Get CH4 mass mixing ratio     call rad_cnst_get_gas(icall,'CH4', pstate, pbuf, ch4)
    ! Get CFC mass mixing ratios    call rad_cnst_get_gas(icall,'CFC11', pstate, pbuf, cfc11)
    !                               call rad_cnst_get_gas(icall,'CFC12', pstate, pbuf, cfc12)
    h2ovmr  =0.0_r8;o3vmr=0.0_r8;co2vmr  =0.0_r8;ch4vmr  =0.0_r8
    o2vmr=0.0_r8;n2ovmr  =0.0_r8;cfc11vmr=0.0_r8;cfc12vmr=0.0_r8
    cfc22vmr=0.0_r8;ccl4vmr =0.0_r8
    DO k = 1, rrtmg_levs

       kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
       DO i = 1,nCols
          ch4vmr  (i,k)   = ch4(i,kk) * amdm
          h2ovmr  (i,k)   = (sp_hum(i,kk) / (1._r8 - sp_hum(i,kk)))! * amdw
          o3vmr   (i,k)   = o3   (i,kk)/amdo
         ! REAL (kind=r8), PARAMETER :: co2vmr_def = 370.0e-6  !( mol/mol).
         ! co2 kg/kg      =  mol/mol
          co2vmr  (i,k)   = co2  (i,kk) * amdc !  REAL(r8), PARAMETER :: amdc  = 0.658114_r8! Molecular weight of dry air / carbon dioxide
          ch4vmr  (i,k)   = ch4  (i,kk) * amdm
          o2vmr   (i,k)   = o2   (i,kk) * amdo2
          n2ovmr  (i,k)   = n2o  (i,kk) * amdn
          cfc11vmr(i,k)   = cfc11(i,kk) * amdc1
          cfc12vmr(i,k)   = cfc12(i,kk) * amdc2
          cfc22vmr(i,k)   = 0._r8
          ccl4vmr (i,k)   = 0._r8
       END DO
    ENDDO

  END SUBROUTINE rrtmg_state_update

  SUBROUTINE Run_Rad_RRTMG_SW(&
                                ! Model Info and flags
       ncols , kmax  , nls   , noz   , &
       icld  , inalb , s0    , cosz  , &
       ratio ,&
                                ! Atmospheric fields
       pl20  , dpl   , tl    , ql    , &
       o3l   , co2l  , gps   , imask , tg    , &
                                ! SURFACE:  albedo
       alvdf , alndf , alvdr , alndr , &
                                ! SW Radiation fields 
       swinc ,                         &
       radvbc, radvdc, radnbc, radndc, &
       radvbl, radvdl, radnbl, radndl, &
       dswclr, dswtop, ssclr , ss    , &
       aslclr, asl   ,                 &
                                ! Cloud field
       cld   , clu   , fice  , &
       rei_in   , rel_in   , taud_in  , &
       cicewp_in,cliqwp_in , &
       E_cld_tau_in,    E_cld_tau_w_in  ,&
       E_cld_tau_w_g_in,E_cld_tau_w_f_in ,cldfprime_in  , &
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013)
               ifaeros,aod,asy,ssa,z_aer,topog, &      
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013)
                aodF,asyF,ssaF,z_aerF )
!tar end       

!       
    IMPLICIT NONE
    ! Input arguments
    !==========================================================================
    !
    ! ___________________________
    ! MODEL INFORMATION AND FLAGS
    !
    ! ncols.....Number of grid points on a gaussian latitude circle
    ! kmax......Number of sigma levels
    ! nls.......Number of layers in the stratosphere.
    ! noz.......
    ! icld......Select three types of cloud emissivity/optical depth
    !           =1 : old cloud emissivity (optical depth) setting
    !                ccu :  0.05 *dp
    !                css :  0.025*dp       for ice cloud t<253.0
    !                       0.05 *dp       for ice cloud t>253.0
    !           =2 : new cloud emissivity (optical depth) setting
    !                ccu : (0.16)*dp
    !                css :  0.0                              t<-82.5
    !                      (2.0e-6*(t-tcrit)**2)*dp    -82.5<t<-10.0
    !                      (6.949e-3*(t-273)+.08)*dp   -10.0<t<  0.0
    !                      (0.08)*dp                   -10.0<t<  0.0
    !           =3    : ccm3 based cloud emissivity
    ! inalb.....Select two types of surface albedo
    !           =1 : input two  types surface albedo (2 diffused)
    !                direct beam albedos are calculated by the subr.
    !           =2 : input four types surfc albedo (2 diff,2 direct)
    ! s0........Solar constant  at proper sun-earth distance
    ! cosz......Cosine of solar zenith angle
    !
    ! __________________
    ! ATMOSPHERIC FIELDS
    !
    !    flip arrays (k=1 means top of atmosphere)
    !
    ! pl20......Flip array of pressure at bottom of layers (mb)
    ! dpl.......Flip array of pressure difference bettween levels
    ! tl........Flip array of temperature (K)
    ! ql........Flip array of specific humidity (g/g)
    ! o3l.......Flip array of ozone mixing ratio (g/g)
    ! gps.......Surface pressure (mb)
    ! imask.....Sea/Land mask
    !
    ! ________________
    ! SURFACE:  albedo
    !
    ! alvdf.....Visible diffuse surface albedo
    ! alndf.....Near-ir diffuse surface albedo
    ! alvdr.....Visible beam surface albedo
    ! alndr.....Near-ir beam surface albedo
    !
    ! ________________________________________
    ! SW Radiation fields 
    !
    ! swinc.....Incident SW at top 

    ! radvbc....Down Sfc SW flux visible beam    (all-sky)
    ! radnbc....Down Sfc SW flux Near-IR beam    (all-sky)
    ! radvdc....Down Sfc SW flux visible diffuse (all-sky)
    ! radndc....Down Sfc SW flux Near-IR diffuse (all-sky)

    ! radvbl....Down Sfc SW flux visible beam    (clear)  
    ! radnbl....Down Sfc SW flux Near-IR beam    (clear)  
    ! radvdl....Down Sfc SW flux visible diffuse (clear)  
    ! radndl....Down Sfc SW flux Near-IR diffuse (clear)  

    ! dswclr....Net SW flux at TOA (clear)   = Abs by Atmos + Sfc
    ! dswtop....Net SW flux at TOA (all-sky) = Abs by Atmos + Sfc
    ! ssclr.....Net SW flux at SFC (clear)   = Abs by Sfc
    ! ss........Net SW flux at SFC (all-sky) = Abs by Sfc

    ! asl.......Heating rate due to shortwave         (K/s)
    ! aslclr....Heating rate due to shortwave (clear) (K/s)
    !
    ! ___________
    ! Cloud field
    !
    ! cld.......Supersaturation cloud fraction
    ! clu.......Convective cloud fraction     
    ! fice......Fraction of cloud water in the form of ice
    ! rei.......Ice particle Effective Radius (microns)
    ! rel.......Liquid particle Effective Radius (microns)
    ! taud......Shortwave cloud optical depth
    !
    !tar begin
    !
    ! climate aerosol parameters of coarse mode (Kinne, 2013)
    ! ifaeros= 0 (original aerosol) =2 (climate aerosol)
    ! aod....aerosol optical depth
    ! asy....asymmetry factor
    ! ssa....single scattering albedo
    ! z_aer.. vertical profile of AOD
    ! topog..topography field  
    !    
    !tar end
    !
    !tar begin 
    ! climate aerosol parameters of fine mode (Kinne, 2013)       
    ! aodF....aerosol optical depth
    ! asyF....asymmetry factor
    ! ssaF....single scattering albedo
    ! z_aerF.. vertical profile of AOD
    !tar end        
    !==========================================================================

    ! Model Info and flags
    INTEGER      ,    INTENT(in   ) :: ncols
    INTEGER      ,    INTENT(in   ) :: kmax
    INTEGER      ,    INTENT(in   ) :: nls
    LOGICAL      ,    INTENT(in   ) :: noz
    INTEGER      ,    INTENT(in   ) :: icld
    INTEGER      ,    INTENT(in   ) :: inalb
    REAL(KIND=r8),    INTENT(in   ) :: s0
    REAL(KIND=r8),    INTENT(in   ) :: cosz   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: ratio

    ! Atmospheric fields
    REAL(KIND=r8),    INTENT(in   ) :: pl20   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: dpl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tl     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: ql     (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: o3l    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: co2l   (ncols,kmax)  !mol/mol
    REAL(KIND=r8),    INTENT(in   ) :: gps    (ncols)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tg     (ncols)

    REAL(r8) :: n2o    (nCols,kMax) ! nitrous oxide mass mixing ratio
    REAL(r8) :: ch4    (nCols,kMax) ! methane  mass mixing ratio
    REAL(r8) :: o2     (nCols,kMax) ! O2    mass mixing ratio
    REAL(r8) :: cfc11  (nCols,kMax) ! cfc11  mass mixing ratio
    REAL(r8) :: cfc12  (nCols,kMax) ! cfc12  mass mixing ratio
    REAL(r8) :: co2    (nCols,kMax) ! co2   mass mixing ratio

    ! SURFACE:  albedo
    REAL(KIND=r8),    INTENT(in   ) :: alvdf  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: alndf  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: alvdr  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: alndr  (ncols)
    ! SW Radiation fields 
    REAL(KIND=r8),    INTENT(out) :: swinc  (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvbc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvdc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radnbc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radndc (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvbl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radvdl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radnbl (ncols)
    REAL(KIND=r8),    INTENT(out) :: radndl (ncols)
    REAL(KIND=r8),    INTENT(out) :: aslclr (ncols,kmax)
    REAL(KIND=r8),    INTENT(out) :: asl    (ncols,kmax)
    REAL(KIND=r8),    INTENT(out) :: ss     (ncols)
    REAL(KIND=r8),    INTENT(out) :: ssclr  (ncols)
    REAL(KIND=r8),    INTENT(out) :: dswtop (ncols)
    REAL(KIND=r8),    INTENT(out) :: dswclr (ncols)

    ! Cloud field and Microphysics
    REAL(KIND=r8),    INTENT(in   ) :: cld    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: clu    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: fice   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: rei_in    (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: rel_in    (ncols,kmax) 
    REAL(KIND=r8),    INTENT(inout) :: taud_in   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) ::  cicewp_in(ncols,kmax) 
    REAL(KIND=r8),    INTENT(in   ) ::  cliqwp_in(ncols,kmax) 
    REAL(KIND=r8), INTENT(in) ::  cldfprime_in (ncols,kmax)

    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_in    (nbndsw, nCols, kMax)      ! cloud optical depth
    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_w_in  (nbndsw, nCols, kMax)      ! cloud optical 
    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_w_g_in(nbndsw, nCols, kMax)      ! cloud optical 
    REAL(KIND=r8) ,   INTENT(in   ):: E_cld_tau_w_f_in(nbndsw, nCols, kMax)      ! cloud optical 

!tar begin
!    Climate aerosol parameters of coarse mode (Kinne, 2013)
    INTEGER,    INTENT(in   ) :: ifaeros
    REAL(KIND=r8),    INTENT(IN) :: aod(nCols,14) 
    REAL(KIND=r8),    INTENT(IN) :: asy(nCols,14)
    REAL(KIND=r8),    INTENT(IN) :: ssa(nCols,14)
    REAL(KIND=r8),    INTENT(IN) :: z_aer(nCols,40)
    REAL(KIND=r8),    INTENT(IN) :: topog (nCols)
!tar end
!
!tar begin
!    Climate aerosol parameters of fine mode (Kinne, 2013)
    REAL(KIND=r8),    INTENT(IN) :: aodF(nCols,14) 
    REAL(KIND=r8),    INTENT(IN) :: asyF(nCols,14)
    REAL(KIND=r8),    INTENT(IN) :: ssaF(nCols,14)
    REAL(KIND=r8),    INTENT(IN) :: z_aerF(nCols,40)
!tar end    
!
    ! local  variables

    !type(rrtmg_state_t), intent(in) :: r_state
    REAL(r8) :: r_state_h2ovmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_o3vmr   (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_co2vmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_ch4vmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_o2vmr   (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_n2ovmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_cfc11vmr(nCols,rrtmg_levs) 
    REAL(r8) :: r_state_cfc12vmr(nCols,rrtmg_levs) 
    REAL(r8) :: r_state_cfc22vmr(nCols,rrtmg_levs) 
    REAL(r8) :: r_state_ccl4vmr (nCols,rrtmg_levs) 

    REAL(r8) :: r_state_pmidmb  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_pintmb  (nCols,rrtmg_levs+1) 
    REAL(r8) :: r_state_tlay    (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_tlev    (nCols,rrtmg_levs+1) 

    REAL(r8) :: state_pmid      (nCols,rrtmg_levs)    ! Level pressure (Pascals)
    REAL(r8) :: cldfprime       (nCols,rrtmg_levs)    ! Fractional cloud cover
    REAL(KIND=r8) :: rei    (ncols,rrtmg_levs)
    REAL(KIND=r8) :: rel    (ncols,rrtmg_levs) 
    REAL(KIND=r8) :: cicewp (ncols,rrtmg_levs) 
    REAL(KIND=r8) :: cliqwp (ncols,rrtmg_levs) 
    REAL(KIND=r8) :: taud   (ncols,rrtmg_levs)
    REAL(KIND=r8) :: E_cld_tau    (nbndsw, nCols, rrtmg_levs)! cloud optical depth
    REAL(KIND=r8) :: E_cld_tau_w  (nbndsw, nCols, rrtmg_levs)! cloud optical 
    REAL(KIND=r8) :: E_cld_tau_w_g(nbndsw, nCols, rrtmg_levs)! cloud optical 
    REAL(KIND=r8) :: E_cld_tau_w_f(nbndsw, nCols, rrtmg_levs)! cloud optical 

    REAL(r8) :: E_aer_tau       (nCols, 0:rrtmg_levs, nbndsw) ! aerosol optical depth
    REAL(r8) :: E_aer_tau_w     (nCols, 0:rrtmg_levs, nbndsw) ! aerosol OD * ssa
    REAL(r8) :: E_aer_tau_w_g   (nCols, 0:rrtmg_levs, nbndsw) ! aerosol OD * ssa * asm
    REAL(r8) :: E_aer_tau_w_f   (nCols, 0:rrtmg_levs, nbndsw) ! aerosol OD * ssa * fwd

    REAL(r8) :: eccf                ! Eccentricity factor (1./earth-sun dist^2)
    REAL(r8) :: E_coszrs(nCols)     ! Cosine solar zenith angle
    REAL(r8) :: E_asdir (nCols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
    REAL(r8) :: E_aldir (nCols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
    REAL(r8) :: E_asdif (nCols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
    REAL(r8) :: E_aldif (nCols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad
    REAL(r8) :: sfac    (nbndsw)    ! factor to account for solar variability in each band 


    ! Output arguments

    REAL(r8) :: solin (nCols)      ! Incident solar flux                   ! swinc.....Incident SW at top 
    REAL(r8) :: qrs   (nCols,rrtmg_levs) ! Solar heating rate                    ! asl.......Heating rate due to shortwave         (K/s)
    REAL(r8) :: qrsc  (nCols,rrtmg_levs) ! Clearsky solar heating rate           ! aslclr....Heating rate due to shortwave (clear) (K/s)

    REAL(r8) :: fsntoa(nCols)      ! Net solar flux at TOA                 ! dswclr....Net SW flux at TOA (clear)   = Abs by Atmos + Sfc
    REAL(r8) :: fsutoa(nCols)      ! Upward solar flux at TOA              ! dswtop....Net SW flux at TOA (all-sky) = Abs by Atmos + Sfc
    REAL(r8) :: fsnrtoac(nCols)    ! Clear sky near-IR flux absorbed at toa

    REAL(r8) :: fsds(nCols)        ! Flux shortwave downwelling surface    ! ssclr.....Net SW flux at SFC (clear)   = Abs by Sfc
    REAL(r8) :: fsns  (nCols)      ! Surface absorbed solar flux           ! ss........Net SW flux at SFC (all-sky) = Abs by Sfc

    REAL(r8) :: fsnt  (nCols)      ! Total column absorbed solar flux

    REAL(r8) :: fsnsc(nCols)       ! Clear sky surface absorbed solar flux       ! radvbl....Down Sfc SW flux visible beam    (clear) 
    REAL(r8) :: fsdsc(nCols)       ! Clear sky surface downwelling solar flux    ! radnbl....Down Sfc SW flux Near-IR beam    (clear) 
    REAL(r8) :: fsntc(nCols)       ! Clear sky total column absorbed solar flx   ! radvdl....Down Sfc SW flux visible diffuse (clear) 
    REAL(r8) :: fsntoac(nCols)     ! Clear sky net solar flx at TOA              ! radndl....Down Sfc SW flux Near-IR diffuse (clear) 

    REAL(r8) :: sols(nCols)        ! Direct  solar rad on surface (< 0.7)      ! radvbc....Down Sfc SW flux visible beam    (all-sky)
    REAL(r8) :: solsd(nCols)       ! Diffuse solar rad on surface (< 0.7)      ! radvdc....Down Sfc SW flux visible diffuse (all-sky) 
    REAL(r8) :: soll(nCols)        ! Direct  solar rad on surface (>= 0.7)     ! radnbc....Down Sfc SW flux Near-IR beam    (all-sky)
    REAL(r8) :: solld(nCols)       ! Diffuse solar rad on surface (>= 0.7)     ! radndc....Down Sfc SW flux Near-IR diffuse (all-sky)

    REAL(r8) :: solscl(nCols)      ! Clear sky  Direct solar rad on surface (< 0.7)
    REAL(r8) :: sollcl(nCols)      ! Clear sky  Direct solar rad on surface (>= 0.7)
    REAL(r8) :: solsdcl(nCols)     ! Clear sky  Diffuse solar rad on surface (< 0.7)
    REAL(r8) :: solldcl(nCols)     ! Clear sky  Diffuse solar rad on surface (>= 0.7)


    REAL(r8) :: fsnirtoa(nCols)     ! Near-IR flux absorbed at toa
    REAL(r8) :: fsnrtoaq(nCols)     ! Net near-IR flux at toa >= 0.7 microns




    REAL(r8) :: fns (nCols,rrtmg_levs+1)  ! net flux at interfaces
    REAL(r8) :: fcns(nCols,rrtmg_levs+1)  ! net clear-sky flux at interfaces

    REAL(r8) :: su(nCols,rrtmg_levs+1,nbndsw)  ! shortwave spectral flux up
    REAL(r8) :: sd(nCols,rrtmg_levs+1,nbndsw)  ! shortwave spectral flux down
    logical  :: old_convert
    !  ---  parameter constants for gas volume mixing ratioes
    ! Gases other than N2, O2, Ar, and H2O are present in the atmosphere at extremely low 
    ! concentrations and are called trace gases. Despite their low concentrations, these 
    ! trace gases can be of critical importance for the greenhouse effect, the ozone layer, 
    ! smog, and other environmental issues. Mixing ratios of trace gases are commonly given
    ! in units of parts per million volume ( ppmv or simply ppm), parts per billion volume ( ppbv or ppb), 
    ! or parts per trillion volume ( pptv or ppt); 1 ppmv = 1x10-6 mol/mol, 1 ppbv = 1x10-9 mol/mol, 
    ! and 1 pptv = 1x10-12 mol/mol. For example, the present-day CO2 concentration is 365 ppmv (365x10-6 mol/mol).
    !
    !
    !            1 gas volume X
    !  1 ppm = --------------------= 1e-6 mol/mol
    !             1e6 air volumes
    !  
    !  CH4
    !   
    !  1 mole C= 12.01 g
    !  4 moles H= 4.04 g
    !  Therefore, 1 mole of CH4 = 12.01g + 4.04g = 16.05g/mol
    !
    !  mwdry       = 28.966_R8       ! molecular weight dry air ~ kg/kmole ~ g/mol

    REAL (kind=r8), PARAMETER :: co2vmr_def = 370.0e-6  !( mol/mol).
    REAL (kind=r8), PARAMETER :: n2ovmr_def = 0.31e-6   !( mol/mol).
    REAL (kind=r8), PARAMETER :: ch4vmr_def = 1.50e-6   !( mol/mol).
    REAL (kind=r8), PARAMETER :: o2vmr_def  = 0.209     !( mol/mol).
    REAL (kind=r8), PARAMETER :: covmr_def  = 1.50e-8   !( mol/mol).
    REAL (kind=r8), PARAMETER :: f11vmr_def = 3.520e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: f12vmr_def = 6.358e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: f22vmr_def = 1.500e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: cl4vmr_def = 1.397e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: f113vmr_def= 8.2000e-11!( mol/mol).  ! gfdl 1999 value

    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: kk
    swinc  =0.0_r8;radvbc =0.0_r8;radvdc =0.0_r8;radnbc =0.0_r8;radndc =0.0_r8;radvbl =0.0_r8;
    radvdl =0.0_r8;radnbl =0.0_r8;radndl =0.0_r8;aslclr =0.0_r8;asl    =0.0_r8;ss     =0.0_r8;
    ssclr  =0.0_r8;dswtop =0.0_r8;dswclr =0.0_r8;state_pmid=0.0_r8;cldfprime=0.0_r8;
    E_aer_tau=0.0_r8;E_aer_tau_w=0.0_r8;E_aer_tau_w_g=0.0_r8;
    E_aer_tau_w_f=0.0_r8;E_coszrs=0.0_r8;E_asdir =0.0_r8;E_aldir =0.0_r8;E_asdif =0.0_r8;E_aldif =0.0_r8;
    sfac    =0.0_r8;solin =0.0_r8;qrs   =0.0_r8;qrsc  =0.0_r8;
    fsntoa=0.0_r8;fsutoa=0.0_r8;fsnrtoac=0.0_r8;fsds=0.0_r8;fsns  =0.0_r8;fsnt =0.0_r8;
    fsnsc=0.0_r8;fsdsc=0.0_r8;fsntc=0.0_r8;fsntoac=0.0_r8;
    sols=0.0_r8;solsd=0.0_r8;soll=0.0_r8;solld=0.0_r8;
    solscl=0.0_r8;sollcl=0.0_r8;solsdcl=0.0_r8;solldcl=0.0_r8;
    fsnirtoa=0.0_r8;fsnrtoaq=0.0_r8;fns =0.0_r8;fcns=0.0_r8;su=0.0_r8;
    rei   =0.0_r8;rel   =0.0_r8;cicewp=0.0_r8;cliqwp=0.0_r8;taud  =0.0_r8;
    E_cld_tau    =0.0_r8;E_cld_tau_w  =0.0_r8;E_cld_tau_w_g=0.0_r8;E_cld_tau_w_f=0.0_r8
    ! alvdf.....Visible diffuse surface albedo
    ! alndf.....Near-ir diffuse surface albedo
    ! alvdr.....Visible beam surface albedo
    ! alndr.....Near-ir beam surface albedo
    DO i=1,ncols
       E_asdir (i)        = alvdr(i) ! 0.2-0.7 micro-meter srfc alb: direct rad visible beam 
       E_aldir (i)        = alndr(i) ! 0.7-5.0 micro-meter srfc alb: direct rad  Near-IR beam    
       E_asdif (i)        = alvdf(i) ! 0.2-0.7 micro-meter srfc alb: diffuse rad visible diffuse 
       E_aldif (i)        = alndf(i) ! 0.7-5.0 micro-meter srfc alb: diffuse rad Near-IR diffuse 
       E_coszrs(i)        = cosz (i) ! Cosine solar zenith angle
    END DO
    DO k=1,nbndsw
       sfac(k)=1.00_r8
    END DO
    eccf=ratio
    CALL rrtmg_state_create(&
         nCols                                    ,&
         kMax                                     ,&
         gps            (1:nCols)                 ,&!mb
         tl             (1:nCols,1:kMax)          ,&
         tg             (1:nCols)                 ,&
         r_state_pmidmb (1:nCols,1:rrtmg_levs)    ,&!mb
         r_state_pintmb (1:nCols,1:rrtmg_levs+1)  ,&!mb
         r_state_tlay   (1:nCols,1:rrtmg_levs)    ,&
         r_state_tlev   (1:nCols,1:rrtmg_levs+1)   )

    DO k=1,rrtmg_levs
       DO i=1,ncols
          state_pmid(i,k)=r_state_pmidmb(i,k)!mb 
       END DO
    END DO

    n2o   = n2ovmr_def   !( mol/mol).
    ch4   = ch4vmr_def
    o2    = o2vmr_def
    cfc11 = f11vmr_def
    cfc12 = f12vmr_def
    DO k=1,kMax
       DO i=1,ncols
          co2(i,k)=co2l(i,k) !mol/mol
         END DO
    END DO
   ! co2   = co2vmr_def

    CALL rrtmg_state_update(&
         nCols                                 , &! INTEGER , INTENT(IN   ) :: nCols
         kMax                                  , &! INTEGER , INTENT(IN   ) :: kMax
         ql              (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: sp_hum  (nCols,kMax) ! specific  humidity   !( kg/kg).
         n2o             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: n2o     (nCols,kMax) ! nitrous oxide mass mixing ratio  !( mol/mol).
         ch4             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: ch4     (nCols,kMax) ! methane  mass mixing ratio !( mol/mol).
         o2              (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: o2      (nCols,kMax) ! O2       mass mixing ratio !( mol/mol).
         cfc11           (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: cfc11   (nCols,kMax) ! cfc11  mass mixing ratio !( mol/mol).
         cfc12           (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: cfc12   (nCols,kMax) ! cfc12  mass mixing ratio !( mol/mol).
         o3l             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: o3      (nCols,kMax) ! Ozone  mass mixing ratio !( kg/kg).
         co2             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: co2     (nCols,kMax) ! co2   mass mixing ratio !( mol/mol).
         r_state_h2ovmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: h2ovmr  (nCols,rrtmg_levs) 
         r_state_o3vmr   (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: o3vmr   (nCols,rrtmg_levs) 
         r_state_co2vmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: co2vmr  (nCols,rrtmg_levs) 
         r_state_ch4vmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: ch4vmr  (nCols,rrtmg_levs) 
         r_state_o2vmr   (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: o2vmr   (nCols,rrtmg_levs) 
         r_state_n2ovmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: n2ovmr  (nCols,rrtmg_levs) 
         r_state_cfc11vmr(1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: cfc11vmr(nCols,rrtmg_levs) 
         r_state_cfc12vmr(1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: cfc12vmr(nCols,rrtmg_levs) 
         r_state_cfc22vmr(1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: cfc22vmr(nCols,rrtmg_levs) 
         r_state_ccl4vmr (1:nCols,1:rrtmg_levs)  )! real(r8), INTENT(OUT  ) :: ccl4vmr (nCols,rrtmg_levs) 

    DO k = 1, rrtmg_levs
       kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
       DO i = 1,nCols
          !cldfprime(i,k)=MAX(cld(i,kk),clu(i,kk)) 
          cldfprime    (i, k) = cldfprime_in(i, kk)
          rei   (i,k)=  rei_in(i,kk)    
          rel   (i,k)=  rel_in(i,kk)   
          cicewp(i,k)=  cicewp_in(i,kk)
          cliqwp(i,k)=  cliqwp_in(i,kk)
          taud  (i,k)=  taud_in(i,kk)
       END DO
    END DO
    DO k = 1, rrtmg_levs
       kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
       DO i = 1,nCols
          DO j=1,nbndsw
             E_cld_tau    (j, i, k) = E_cld_tau_in    (j, i, kk)! cloud optical depth
             E_cld_tau_w  (j, i, k) = E_cld_tau_w_in  (j, i, kk)! cloud optical 
             E_cld_tau_w_g(j, i, k) = E_cld_tau_w_g_in(j, i, kk)! cloud optical 
             E_cld_tau_w_f(j, i, k) = E_cld_tau_w_f_in(j, i, kk)! cloud optical 
          END DO
       END DO
    END DO

    DO j=1,nbndsw
       DO k = 0, rrtmg_levs
          kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
          DO i=1,nCols
                            !(nCols, 0:rrtmg_levs, nbndsw)
             E_aer_tau      (i,k,j) = 0.0_r8
             E_aer_tau_w    (i,k,j) = 0.0_r8
             E_aer_tau_w_g  (i,k,j) = 0.0_r8
             E_aer_tau_w_f  (i,k,j) = 0.0_r8
          END DO
       END DO
    END DO
    old_convert = .FALSE. 
    CALL rad_rrtmg_sw(&
         imca_sw                                     , &!
         nCols                                    , &!integer , intent(in) :: nCols
         rrtmg_levs                               , &!integer , intent(in) :: kMax
         rrtmg_levs+1                             , &!integer , intent(in) :: kMax+1
         nCols                                    , &!integer , intent(in) :: ncol              ! number of atmospheric columns
         rrtmg_levs                               , &!integer , intent(in) :: rrtmg_levs        ! number of levels rad is applied
         FeedBackOptics_cld                       , &
         r_state_h2ovmr   (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_h2ovmr  (nCols,rrtmg_levs) 
         r_state_o3vmr    (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_o3vmr   (nCols,rrtmg_levs) 
         r_state_co2vmr   (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_co2vmr  (nCols,rrtmg_levs) 
         r_state_ch4vmr   (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_ch4vmr  (nCols,rrtmg_levs) 
         r_state_o2vmr    (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_o2vmr   (nCols,rrtmg_levs) 
         r_state_n2ovmr   (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_n2ovmr  (nCols,rrtmg_levs) 
         r_state_cfc11vmr (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_cfc11vmr(nCols,rrtmg_levs) 
         r_state_cfc12vmr (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_cfc12vmr(nCols,rrtmg_levs) 
         r_state_cfc22vmr (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_cfc22vmr(nCols,rrtmg_levs) 
         r_state_ccl4vmr  (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_ccl4vmr (nCols,rrtmg_levs) 
         r_state_pmidmb   (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_pmidmb  (nCols,rrtmg_levs) 
         r_state_pintmb   (1:nCols,1:rrtmg_levs+1), &!real(r8), intent(in) :: r_state_pintmb  (nCols,rrtmg_levs+1)
         r_state_tlay     (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: r_state_tlay    (nCols,rrtmg_levs) 
         r_state_tlev     (1:nCols,1:rrtmg_levs+1), &!real(r8), intent(in) :: r_state_tlev    (nCols,rrtmg_levs+1)
         state_pmid       (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: E_pmid     (nCols,rrtmg_levs)   ! Level pressure (hPascals)
         cldfprime        (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(in) :: E_cld     (nCols,rrtmg_levs)   ! Fractional cloud cover
         E_aer_tau        (1:nCols, 0:rrtmg_levs,1: nbndsw), &!real(r8), intent(in) :: E_aer_tau       (nCols, 0:rrtmg_levs, nbndsw)! aerosol optical depth
         E_aer_tau_w      (1:nCols, 0:rrtmg_levs,1: nbndsw), &!real(r8), intent(in) :: E_aer_tau_w     (nCols, 0:rrtmg_levs, nbndsw)! aerosol OD * ssa
         E_aer_tau_w_g    (1:nCols, 0:rrtmg_levs,1: nbndsw), &!real(r8), intent(in) :: E_aer_tau_w_g   (nCols, 0:rrtmg_levs, nbndsw)! aerosol OD * ssa * asm
         E_aer_tau_w_f    (1:nCols, 0:rrtmg_levs,1: nbndsw), &!real(r8), intent(in) :: E_aer_tau_w_f   (nCols, 0:rrtmg_levs, nbndsw)! aerosol OD * ssa * fwd
         eccf                                     , &!real(r8), intent(in) :: eccf                         ! Eccentricity factor (1./earth-sun dist^2)
         E_coszrs         (1:nCols)               , &!real(r8), intent(in) :: E_coszrs   (nCols)! Cosine solar zenith angle
         solin            (1:nCols)               , &!real(r8), intent(out):: solin   (nCols)! Incident solar flux
         sfac             (1:nbndsw)              , &!real(r8), intent(in) :: sfac    (nbndsw)! factor to account for solar variability in each band 
         E_asdir          (1:nCols)               , &!real(r8), intent(in) :: E_asdir   (nCols)! 0.2-0.7 micro-meter srfc alb: direct rad
         E_asdif          (1:nCols)               , &!real(r8), intent(in) :: E_asdif   (nCols)! 0.2-0.7 micro-meter srfc alb: diffuse rad
         E_aldir          (1:nCols)               , &!real(r8), intent(in) :: E_aldir   (nCols)! 0.7-5.0 micro-meter srfc alb: direct rad
         E_aldif          (1:nCols)               , &!real(r8), intent(in) :: E_aldif   (nCols)! 0.7-5.0 micro-meter srfc alb: diffuse rad
         qrs               (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(out):: qrs   (nCols,rrtmg_levs) ! Solar heating rate
         qrsc              (1:nCols,1:rrtmg_levs)  , &!real(r8), intent(out):: qrsc  (nCols,rrtmg_levs) ! Clearsky solar heating rate
         fsnt              (1:nCols)               , &!real(r8), intent(out):: fsnt         (nCols)      ! Total column absorbed solar flux
         fsntc             (1:nCols)               , &!real(r8), intent(out):: fsntc         (nCols)      ! Clear sky total column absorbed solar flx
         fsntoa            (1:nCols)               , &!real(r8), intent(out):: fsntoa         (nCols)      ! Net solar flux at TOA
         fsutoa            (1:nCols)               , &!real(r8), intent(out):: fsutoa         (nCols)      ! Upward solar flux at TOA
         fsntoac           (1:nCols)               , &!real(r8), intent(out):: fsntoac         (nCols)      ! Clear sky net solar flx at TOA
         fsnirtoa          (1:nCols)               , &!real(r8), intent(out):: fsnirtoa        (nCols)      ! Near-IR flux absorbed at toa
         fsnrtoac          (1:nCols)               , &!real(r8), intent(out):: fsnrtoac        (nCols)      ! Clear sky near-IR flux absorbed at toa
         fsnrtoaq          (1:nCols)               , &!real(r8), intent(out):: fsnrtoaq        (nCols)      ! Net near-IR flux at toa >= 0.7 microns
         fsns              (1:nCols)               , &!real(r8), intent(out):: fsns         (nCols)      ! Surface absorbed solar flux
         fsnsc             (1:nCols)               , &!real(r8), intent(out):: fsnsc         (nCols)      ! Clear sky surface absorbed solar flux
         fsdsc             (1:nCols)               , &!real(r8), intent(out):: fsdsc         (nCols)      ! Clear sky surface downwelling solar flux
         fsds              (1:nCols)               , &!real(r8), intent(out):: fsds         (nCols)      ! Flux shortwave downwelling surface
         sols              (1:nCols)               , &!real(r8), intent(out):: sols         (nCols)      ! Direct solar rad on surface  (<  0.7)
         soll              (1:nCols)               , &!real(r8), intent(out):: soll         (nCols)      ! Direct solar rad on surface  (>= 0.7)
         solsd             (1:nCols)               , &!real(r8), intent(out):: solsd         (nCols)      ! Diffuse solar rad on surface (<  0.7)
         solld             (1:nCols)               , &!real(r8), intent(out):: solld         (nCols)      ! Diffuse solar rad on surface (>= 0.7)
         solscl            (1:nCols)               , &!real(r8), intent(out):: solscl         (nCols)      ! Clear sky Direct solar rad on surface  (<  0.7)
         sollcl            (1:nCols)               , &!real(r8), intent(out):: sollcl         (nCols)      ! Clear sky Direct solar rad on surface  (>= 0.7)
         solsdcl           (1:nCols)               , &!real(r8), intent(out):: solsdcl         (nCols)      ! Clear sky Diffuse solar rad on surface (<  0.7)
         solldcl           (1:nCols)               , &!real(r8), intent(out):: solldcl         (nCols)      ! Clear sky Diffuse solar rad on surface (>= 0.7)
         fns               (1:nCols,1:kMax+1)      , &!real(r8), intent(out):: fns             (nCols,kMax+1)! net flux at interfaces
         fcns              (1:nCols,1:kMax+1)      , &!real(r8), intent(out):: fcns            (nCols,kMax+1)! net clear-sky flux at interfaces
         su                (1:nCols,1:rrtmg_levs+1,1:nbndsw), &!real(r8), pointer :: su      (nCols,rrtmg_levs+1,nbndsw)! shortwave spectral flux up
         sd                (1:nCols,1:rrtmg_levs+1,1:nbndsw), &!real(r8), pointer :: sd      (nCols,rrtmg_levs+1,nbndsw)! shortwave spectral flux down
         rei               (1:ncols,1:rrtmg_levs)    , &
         rel               (1:ncols,1:rrtmg_levs)    , &
         cicewp            (1:ncols,1:rrtmg_levs)    , &
         cliqwp            (1:ncols,1:rrtmg_levs)    , &
         taud              (1:ncols,1:rrtmg_levs)    , &
         E_cld_tau         (1:nbndsw, 1:nCols, 1:rrtmg_levs)  , &!real(r8), intent(in) :: E_cld_tau       (nbndsw, nCols, rrtmg_levs)      ! cloud optical depth
         E_cld_tau_w       (1:nbndsw, 1:nCols, 1:rrtmg_levs)  , &!real(r8), intent(in) :: E_cld_tau_w     (nbndsw, nCols, rrtmg_levs)      ! cloud optical 
         E_cld_tau_w_g     (1:nbndsw, 1:nCols, 1:rrtmg_levs)  , &!real(r8), intent(in) :: E_cld_tau_w_g   (nbndsw, nCols, rrtmg_levs)      ! cloud optical 
         E_cld_tau_w_f     (1:nbndsw, 1:nCols, 1:rrtmg_levs)  , &!real(r8), intent(in) :: E_cld_tau_w_f   (nbndsw, nCols, rrtmg_levs)      ! cloud optical 
         old_convert                                          , &! logical , optional, intent(in)   :: old_convert
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013)
               ifaeros,aod,asy,ssa,z_aer,topog, &      
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013)
                aodF,asyF,ssaF,z_aerF )
!tar end
!

    DO i = 1,nCols
       swinc  (i)=solin           (i)      ! Incident solar flux                    ! swinc.....Incident SW at top 
       radvbc (i)=sols            (i)      ! Direct  solar rad on surface  (<  0.7) ! radvbc....Down Sfc SW flux visible beam    (all-sky)
       radvdc (i)=solsd           (i)      ! Diffuse solar rad on surface  (<  0.7) ! 
       radnbc (i)=soll            (i)      ! Direct  solar rad on surface  (>= 0.7) ! 
       radndc (i)=solld           (i)      ! Diffuse solar rad on surface  (>= 0.7) ! 
       radvbl (i)=solscl          (i)      ! Clear sky Direct solar rad on surface  ! (<  0.7) ! radvbl....Down Sfc SW flux visible beam    (clear)  
       radvdl (i)=solsdcl         (i)      ! Clear sky Diffuse solar rad on surface ! (<  0.7) ! radvdl....Down Sfc SW flux visible diffuse (clear)  
       radnbl (i)=sollcl          (i)      ! Clear sky Direct solar rad on surface  ! (>= 0.7) ! radnbl....Down Sfc SW flux Near-IR beam    (clear)  
       radndl (i)=solldcl         (i)      ! Clear sky Diffuse solar rad on surface ! (>= 0.7) ! radndl....Down Sfc SW flux Near-IR diffuse (clear)  
       ss     (i) = fsns          (i)      !           Surface absorbed solar flux  ! ss........Net SW flux at SFC (all-sky) = Abs by Sfc
       ssclr  (i)=fsnsc           (i)      ! Clear sky surface absorbed solar flux  ! ssclr.....Net SW flux at SFC (clear)   = Abs by Sfc
       dswtop (i)=fsntoa          (i)      !           Net solar flux at TOA        ! dswtop....Net SW flux at TOA (all-sky) = Abs by Atmos + Sfc
       dswclr (i)=fsntoac         (i)      ! Clear sky net solar flux at TOA        ! dswclr....Net SW flux at TOA (clear)   = Abs by Atmos + Sfc
    END DO

    DO k = 1, rrtmg_levs
       kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
       DO i = 1,nCols
          aslclr (i,kk) = qrsc (i,k)  ! asl.......Heating rate due to shortwave         (K/s)
          asl    (i,kk) = qrs  (i,k)  ! aslclr....Heating rate due to shortwave (clear) (K/s)
       END DO
    END DO

  END SUBROUTINE Run_Rad_RRTMG_SW


  !---------------------------------------------------------------------------
  SUBROUTINE Run_Rad_RRTMG_LW( &
                                ! Model Info and flags
       ncols         ,kmax      ,                             &
                                ! Atmospheric fields
       Pbot          ,Pmid      ,DP            ,Te          , &
       Qe            ,O3        ,co2l        ,Tg          , &
       gps           ,&
                                ! SURFACE
       imask         ,                                        &
                                ! LW Radiation fields 
       lw_toa_up_clr ,lw_toa_up ,lw_cool_clr    ,lw_cool    , &
       lw_sfc_net_clr,lw_sfc_net,lw_sfc_down_clr,lw_sfc_down, &
                                ! Cloud field
       cld           ,clu       ,fice           ,             &
       rei_in           ,rel_in       ,lmixr  ,tauc   ,cicewp_in  ,cliqwp_in ,cldfprime_in                       )

    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! INPUT VARIABLES
    !---------------------------------------------------------------------------

    ! Number of atmospheric columns to solve
    INTEGER, INTENT(in) :: ncols
    ! Number of atmospheric levels
    INTEGER, INTENT(in) :: kmax

    REAL(KIND=r8), INTENT(in) ::  gps(ncols)  !Surface Pressure (mb)
    REAL(KIND=r8), INTENT(in) ::  Tg (ncols) ! Ground surface temperature (K)
    REAL(KIND=r8), INTENT(in) ::  Pbot(ncols,kmax)  ! Pressure at bottom of layers (mb)
    REAL(KIND=r8), INTENT(in) ::  Pmid(ncols,kmax)  ! Pressure at Middle of Layer(mb)
    REAL(KIND=r8), INTENT(in) ::  DP(ncols,kmax) ! Pressure difference bettween levels (mb)
    REAL(KIND=r8), INTENT(in) ::  Te (ncols,kmax) ! Temperature at middle of Layer (K)
    REAL(KIND=r8), INTENT(in) ::  Qe (ncols,kmax) ! Specific Humidity at middle of layer (g/g)
    REAL(KIND=r8), INTENT(in) ::  O3 (ncols,kmax)! Ozone Mixing ratio at middle of layer (g/g)
    REAL(KIND=r8), INTENT(in) ::  cld (ncols,kmax)! Large scale cloud amount in layers
    REAL(KIND=r8), INTENT(in) ::  clu (ncols,kmax) ! Cumulus cloud amount in layers
    REAL(KIND=r8), INTENT(in) ::  fice (ncols,kmax)! fractional amount of cloud that is ice
    REAL(KIND=r8), INTENT(in) ::  rei_in (ncols,kmax)! Ice particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(in) ::  rel_in (ncols,kmax)! Liquid particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(in) ::  lmixr(ncols,kmax) ! ice/water mixing ratio
    REAL(KIND=r8), INTENT(in) ::  tauc(nbndlw,ncols,kmax)
    REAL(KIND=r8), INTENT(in) ::  cicewp_in (ncols,kmax)
    REAL(KIND=r8), INTENT(in) ::  cliqwp_in (ncols,kmax)
    REAL(KIND=r8), INTENT(in) ::  cldfprime_in (ncols,kmax)
    REAL(KIND=r8), INTENT(in) ::  co2l (ncols,kmax)! co2 concentration in ppmv(??)
    INTEGER(KIND=i8), INTENT(IN) :: imask (ncols)  ! vegetagion mask (0=ocean, -1=sea ice and 13=land ice)
    REAL(r8) :: n2o    (nCols,kMax) ! nitrous oxide mass mixing ratio
    REAL(r8) :: ch4    (nCols,kMax) ! methane  mass mixing ratio
    REAL(r8) :: o2     (nCols,kMax) ! O2    mass mixing ratio
    REAL(r8) :: cfc11  (nCols,kMax) ! cfc11  mass mixing ratio
    REAL(r8) :: cfc12  (nCols,kMax) ! cfc12  mass mixing ratio
    REAL(r8) :: co2    (nCols,kMax) ! co2   mass mixing ratio

    !---------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    !---------------------------------------------------------------------------


    REAL(KIND=r8), INTENT(inout) :: lw_toa_up_clr  (ncols)    ! Upward TOA longwave flux (clear)
    REAL(KIND=r8), INTENT(inout) :: lw_toa_up      (ncols)         ! Upward TOA longwave flux
    REAL(KIND=r8), INTENT(inout) :: lw_sfc_net_clr (ncols)    ! Net Surface longwave flux (clear)
    REAL(KIND=r8), INTENT(inout) :: lw_sfc_net     (ncols)        ! Net Surface longwave flux
    REAL(KIND=r8), INTENT(inout) :: lw_sfc_down_clr(ncols)   ! Downward Surface longwave flux (clear)
    REAL(KIND=r8), INTENT(inout) :: lw_sfc_down    (ncols)       ! Downward Surface longwave flux
    REAL(KIND=r8), INTENT(inout) :: lw_cool_clr    (ncols,kmax) ! Cooling rate (K/s, clear)
    REAL(KIND=r8), INTENT(inout) :: lw_cool        (ncols,kmax)     ! Cooling rate (K/s)


!    INTEGER, PARAMETER :: nbndlw=16

    !type(rrtmg_state_t), intent(in) :: r_state
    REAL(r8) :: r_state_o3vmr   (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_co2vmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_ch4vmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_o2vmr   (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_n2ovmr  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_cfc11vmr(nCols,rrtmg_levs) 
    REAL(r8) :: r_state_cfc12vmr(nCols,rrtmg_levs) 
    REAL(r8) :: r_state_cfc22vmr(nCols,rrtmg_levs) 
    REAL(r8) :: r_state_ccl4vmr (nCols,rrtmg_levs) 

    REAL(r8) :: r_state_pmidmb  (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_pintmb  (nCols,rrtmg_levs+1) 
    REAL(r8) :: r_state_tlay    (nCols,rrtmg_levs) 
    REAL(r8) :: r_state_tlev    (nCols,rrtmg_levs+1) 
    REAL(r8) :: r_state_h2ovmr  (nCols,rrtmg_levs) 

    REAL(r8) :: state_pmid(nCols,rrtmg_levs)    ! Level pressure (Pascals)


    ! Output arguments

    REAL(r8) :: tauc_lw      (nbndlw,nCols,rrtmg_levs)   ! Cloud longwave optical depth by band
    REAL(r8) :: c_cld_lw_abs (nCols,rrtmg_levs,nbndlw) ! cloud absorption optics depth   (LW)
    REAL(r8) :: qrl          (nCols,rrtmg_levs)     !, intent(out) Longwave heating rate
    REAL(r8) :: qrlc         (nCols,rrtmg_levs)     !, intent(out) Clearsky longwave heating rate
    REAL(r8) :: flns         (nCols)    !, intent(out)Surface cooling flux
    REAL(r8) :: flnsc        (nCols)    !, intent(out)Clear sky surface cooing
    REAL(r8) :: flnt         (nCols)    !, intent(out)Net outgoing flux
    REAL(r8) :: flntc        (nCols)    !, intent(out)Net clear sky outgoing flux
    REAL(r8) :: flwds        (nCols)         !, intent(out)cam_out%flwds,  Down longwave flux at surface
    REAL(r8) :: flut         (nCols)    !, intent(out)Upward flux at top of model
    REAL(r8) :: flutc        (nCols)    !, intent(out)Upward clear-sky flux at top of model
    REAL(r8) :: fnl          (nCols,rrtmg_levs+1)! , intent(out)    ! net flux at interfaces
    REAL(r8) :: fcnl         (nCols,rrtmg_levs+1)    ! , intent(out)clear sky net flux at interfaces
    REAL(r8) :: fldsc        (nCols)         ! , intent(out)Down longwave clear flux at surface

    REAL(r8) :: lu(nCols,rrtmg_levs+1,nbndlw)  ! shortwave spectral flux up
    REAL(r8) :: ld(nCols,rrtmg_levs+1,nbndlw)  ! shortwave spectral flux down
    REAL(KIND=r8) ::  rei     (nCols,rrtmg_levs)! Ice particle Effective Radius (microns)
    REAL(KIND=r8) ::  rel     (nCols,rrtmg_levs)! Liquid particle Effective Radius (microns)
    REAL(KIND=r8) ::  cicewp  (nCols,rrtmg_levs)
    REAL(KIND=r8) ::  cliqwp  (nCols,rrtmg_levs)
    REAL(KIND=r8) :: cldfprime (ncols,rrtmg_levs)


    REAL (kind=r8), PARAMETER :: co2vmr_def = 370.0e-6  !( mol/mol).
    REAL (kind=r8), PARAMETER :: n2ovmr_def = 0.31e-6   !( mol/mol).
    REAL (kind=r8), PARAMETER :: ch4vmr_def = 1.50e-6   !( mol/mol).
    REAL (kind=r8), PARAMETER :: o2vmr_def  = 0.209     !( mol/mol).
    REAL (kind=r8), PARAMETER :: covmr_def  = 1.50e-8   !( mol/mol).
    REAL (kind=r8), PARAMETER :: f11vmr_def = 3.520e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: f12vmr_def = 6.358e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: f22vmr_def = 1.500e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: cl4vmr_def = 1.397e-10 !( mol/mol).  ! aer 2003 value
    REAL (kind=r8), PARAMETER :: f113vmr_def= 8.2000e-11!( mol/mol).  ! gfdl 1999 value
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: kk
    lw_toa_up_clr  =0.0_r8
    lw_toa_up      =0.0_r8
    lw_sfc_net_clr =0.0_r8
    lw_sfc_net     =0.0_r8
    lw_sfc_down_clr=0.0_r8
    lw_sfc_down    =0.0_r8
    lw_cool_clr    =0.0_r8
    lw_cool        =0.0_r8
    state_pmid  =0.0_r8
    cldfprime    =0.0_r8
    tauc_lw        =0.0_r8
    c_cld_lw_abs   =0.0_r8
    qrl            =0.0_r8
    qrlc           =0.0_r8
    flns           =0.0_r8
    flnsc          =0.0_r8
    flnt           =0.0_r8
    flntc          =0.0_r8
    flwds          =0.0_r8
    flut           =0.0_r8
    flutc          =0.0_r8
    fnl            =0.0_r8
    fcnl           =0.0_r8
    fldsc          =0.0_r8
    rei     =0.0_r8
    rel     =0.0_r8
    cicewp  =0.0_r8
    cliqwp  =0.0_r8

    CALL rrtmg_state_create(&
         nCols                                    ,&
         kMax                                     ,&
         gps            (1:nCols)                 ,&!mb
         Te             (1:nCols,1:kMax)          ,&
         Tg             (1:nCols)                 ,&
         r_state_pmidmb (1:nCols,1:rrtmg_levs)    ,&!mb
         r_state_pintmb (1:nCols,1:rrtmg_levs+1)  ,&!mb
         r_state_tlay   (1:nCols,1:rrtmg_levs)    ,&
         r_state_tlev   (1:nCols,1:rrtmg_levs+1)   )



    n2o   = n2ovmr_def   !( mol/mol).
    ch4   = ch4vmr_def !( mol/mol).
    o2    = o2vmr_def !( mol/mol).
    cfc11 = f11vmr_def !( mol/mol).
    cfc12 = f12vmr_def !( mol/mol).
    DO k=1,kMax
       DO i=1,ncols
          co2(i,k)=co2l(i,k) !mol/mol
         END DO
    END DO
    !co2   = co2vmr_def !( mol/mol).
    CALL rrtmg_state_update(&
         nCols                                 , &! INTEGER , INTENT(IN   ) :: nCols
         kMax                                  , &! INTEGER , INTENT(IN   ) :: kMax
         qe              (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: sp_hum  (nCols,kMax) ! specific  humidity   !( kg/kg).
         n2o             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: n2o     (nCols,kMax) ! nitrous oxide mass mixing ratio  !( mol/mol).
         ch4             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: ch4     (nCols,kMax) ! methane  mass mixing ratio !( mol/mol).
         o2              (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: o2      (nCols,kMax) ! O2       mass mixing ratio !( mol/mol).
         cfc11           (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: cfc11   (nCols,kMax) ! cfc11  mass mixing ratio !( mol/mol).
         cfc12           (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: cfc12   (nCols,kMax) ! cfc12  mass mixing ratio !( mol/mol).
         O3              (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: o3      (nCols,kMax) ! Ozone  mass mixing ratio !( kg/kg).
         co2             (1:nCols,1:kMax)      , &! real(r8), INTENT(IN   ) :: co2     (nCols,kMax) ! co2   mass mixing ratio !( mol/mol).
         r_state_h2ovmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: h2ovmr  (nCols,rrtmg_levs) 
         r_state_o3vmr   (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: o3vmr   (nCols,rrtmg_levs) 
         r_state_co2vmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: co2vmr  (nCols,rrtmg_levs) 
         r_state_ch4vmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: ch4vmr  (nCols,rrtmg_levs) 
         r_state_o2vmr   (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: o2vmr   (nCols,rrtmg_levs) 
         r_state_n2ovmr  (1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: n2ovmr  (nCols,rrtmg_levs) 
         r_state_cfc11vmr(1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: cfc11vmr(nCols,rrtmg_levs) 
         r_state_cfc12vmr(1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: cfc12vmr(nCols,rrtmg_levs) 
         r_state_cfc22vmr(1:nCols,1:rrtmg_levs), &! real(r8), INTENT(OUT  ) :: cfc22vmr(nCols,rrtmg_levs) 
         r_state_ccl4vmr (1:nCols,1:rrtmg_levs)  )! real(r8), INTENT(OUT  ) :: ccl4vmr (nCols,rrtmg_levs) 

    DO j=1,nbndlw
       DO k = 1, rrtmg_levs
          kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
          DO i = 1,nCols
             c_cld_lw_abs    (i,k,j)    = 0.0_r8
          END DO
       END DO
    END DO

    DO j=1,nbndlw
       DO k = 1, rrtmg_levs
          kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
          DO i = 1,nCols
             tauc_lw  (j,i,k)=tauc(j,i,kk)
             !cldfprime(i,k)=MAX(cld(i,kk),clu(i,kk)) 
             cldfprime(i,k)=cldfprime_in(i,kk)    
             rei   (i,k)=  rei_in(i,kk)    
             rel   (i,k)=  rel_in(i,kk)   
             cicewp(i,k)=  cicewp_in(i,kk)
             cliqwp(i,k)=  cliqwp_in(i,kk)
          END DO
       END DO
    END DO
    DO k=1,rrtmg_levs
       DO i=1,ncols
          state_pmid(i,k)=r_state_pmidmb(i,k)!              mb
       END DO
    END DO

    CALL rad_rrtmg_lw( &
         imca_lw                                         , &!
         nCols                                           , &
         rrtmg_levs                                      , &
         r_state_pmidmb  (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_pmidmb  (nCols,kMax)
         r_state_pintmb  (1:nCols,1:rrtmg_levs+1)        , &! real(r8), intent(in ) :: r_state_pintmb  (nCols,kMax+1)
         r_state_tlay    (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_tlay    (nCols,kMax) 
         r_state_tlev    (1:nCols,1:rrtmg_levs+1)        , &! real(r8), intent(in ) :: r_state_tlev    (nCols,kMax+1)
         r_state_h2ovmr  (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_h2ovmr  (nCols,kMax)
         r_state_o3vmr   (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_o3vmr   (nCols,kMax)
         r_state_co2vmr  (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_co2vmr  (nCols,kMax)
         r_state_ch4vmr  (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_ch4vmr  (nCols,kMax)
         r_state_o2vmr   (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_o2vmr   (nCols,kMax)
         r_state_n2ovmr  (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_n2ovmr  (nCols,kMax)
         r_state_cfc11vmr(1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_cfc11vmr(nCols,kMax)
         r_state_cfc12vmr(1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_cfc12vmr(nCols,kMax)
         r_state_cfc22vmr(1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_cfc22vmr(nCols,kMax)
         r_state_ccl4vmr (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: r_state_ccl4vmr (nCols,kMax)
         state_pmid      (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: pmid            (nCols,kMax)          ! Level pressure (hPascals)
         c_cld_lw_abs    (1:nCols,1:rrtmg_levs,1:nbndlw) , &! real(r8), intent(in ) :: aer_lw_abs      (nCols,kMax,nbndlw)   ! aerosol absorption optics depth (LW)
         cldfprime       (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(in ) :: cld             (nCols,kMax)          ! Cloud cover
         cicewp          (1:nCols,1:rrtmg_levs)          , &!
         cliqwp          (1:nCols,1:rrtmg_levs)          , &!
         rei             (1:nCols,1:rrtmg_levs)          , &!
         rel             (1:nCols,1:rrtmg_levs)          , &!
         tauc_lw         (1:nbndlw,1:nCols,1:rrtmg_levs) , &! real(r8), intent(in ) :: tauc_lw         (nbndlw,nCols,kMax)   ! Cloud longwave optical depth by band
         qrl             (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(out) :: qrl             (nCols,kMax)          ! Longwave heating rate
         qrlc            (1:nCols,1:rrtmg_levs)          , &! real(r8), intent(out) :: qrlc            (nCols,kMax)          ! Clearsky longwave heating rate
         flns            (1:nCols)                       , &! real(r8), intent(out) :: flns            (nCols)               ! Surface cooling flux
         flnt            (1:nCols)                       , &! real(r8), intent(out) :: flnt            (nCols)               ! Net outgoing flux
         flnsc           (1:nCols)                       , &! real(r8), intent(out) :: flnsc           (nCols)               ! Clear sky surface cooing
         flntc           (1:nCols)                       , &! real(r8), intent(out) :: flntc           (nCols)! Net clear sky outgoing flux
         flwds           (1:nCols)                       , &! real(r8), intent(out) :: flwds           (nCols)! Down longwave flux at surface
         flut            (1:nCols)                       , &! real(r8), intent(out) :: flut            (nCols)! Upward flux at top of model
         flutc           (1:nCols)                       , &! real(r8), intent(out) :: flutc           (nCols)         ! Upward clear-sky flux at top of model
         fnl             (1:nCols,1:rrtmg_levs+1)        , &! real(r8), intent(out) :: fnl             (nCols,kMax+1)! net flux at interfaces
         fcnl            (1:nCols,1:rrtmg_levs+1)        , &! real(r8), intent(out) :: fcnl            (nCols,kMax+1)! clear sky net flux at interfaces
         fldsc           (1:nCols)                       , &! real(r8), intent(out) :: fldsc           (nCols)! Down longwave clear flux at surface
         lu              (1:nCols,1:rrtmg_levs+1,1:nbndlw), &! real(r8), intent(out) :: lu             (nCols,kMax+1,nbndlw)    ! longwave spectral flux up
         ld              (1:nCols,1:rrtmg_levs+1,1:nbndlw), &! real(r8), intent(out) :: ld             (nCols,kMax+1,nbndlw)    ! longwave spectral flux down
         FeedBackOptics_cld                               , &
         nCols                                            , &
         rrtmg_levs                                       , &
         rrtmg_levs+1      )

    DO i = 1,nCols
       lw_toa_up_clr  (i) = flutc(i)          ! Upward clear-sky flux at top of model  ! Upward TOA longwave flux (clear)
       lw_toa_up      (i) = flut (i)          ! Upward           flux at top of model  ! Upward TOA longwave flux
       lw_sfc_net_clr (i) = flntc(i)          ! Net clear sky outgoing flux            ! Net Surface longwave flux (clear)
       lw_sfc_net     (i) = flnt (i)          ! Net outgoing flux                      ! Net Surface longwave flux
       lw_sfc_down_clr(i) = fldsc(i)          ! Down longwave clear flux at surface    ! Downward Surface longwave flux (clear)
       lw_sfc_down    (i) = flwds(i)          ! Down longwave flux at surface          ! Downward Surface longwave flux
      ! WRITE(*,*)fldsc(i), flwds(i) 
    END DO

    DO k = 1, rrtmg_levs
       kk = MAX(k + ((kMax+1)-rrtmg_levs)-1,1)
       DO i = 1,nCols
          lw_cool_clr (i,kk)  = qrlc (i,k)     ! Clearsky longwave heating rate  ! Cooling rate (K/s, clear)
          lw_cool     (i,kk)  = qrl  (i,k)     ! Cooling rate (K/s)
       
       END DO
    END DO


  END SUBROUTINE Run_Rad_RRTMG_LW

  SUBROUTINE Finalize_Rad_RRTMG()
    IMPLICIT NONE

  END SUBROUTINE Finalize_Rad_RRTMG

END MODULE Rad_RRTMG
