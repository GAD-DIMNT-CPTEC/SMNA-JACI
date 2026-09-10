MODULE GwddSchemeCPTEC
  IMPLICIT NONE
  SAVE

   PRIVATE
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(31)! Kind for 128-bits Real Numbers

  REAL(r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
  REAL(R8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
  REAL(R8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
  REAL(R8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_R8       ! molecular weight dry air ~ kg/kmole
  REAL(r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  REAL(R8),PARAMETER :: SHR_CONST_CPWV   = 1.810e3_R8      ! specific heat of water vap ~ J/kg/K
  REAL(R8),PARAMETER :: SHR_CONST_G      = 9.80616_R8      ! acceleration of gravity ~ m/s^2
  REAL(R8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  REAL(R8),PARAMETER :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  REAL(r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg

  REAL (KIND=r8), PARAMETER   :: cp    =                  1004.6_r8! specific heat of air           (j/kg/k)
  REAL (KIND=r8), PARAMETER   :: gasr  =                  287.05_r8! gas constant of dry air        (j/kg/k)
  REAL(r8), PUBLIC, PARAMETER :: gravit = shr_const_g      ! gravitational acceleration

  REAL(r8), PUBLIC, PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: rair = shr_const_rdair    ! Gas constant for dry air (J/K/kg)
  REAL(r8), PUBLIC, PARAMETER :: cpwv  = shr_const_cpwv
  REAL(r8), PUBLIC, PARAMETER :: zvir = SHR_CONST_RWV/rair - 1          ! rh2o/rair - 1
  REAL(KIND=r8), PARAMETER :: Pi                 = 3.14159265358979323846_r8


  ! gw_vert ---- gw_satn
  !              |
  !              gw_wake
  !
  !*L------------------COMDECK C_A----------------------------------------
  ! History:
  ! Version  Date      Comment.
  !  5.0  07/05/99  Replace variable A by more meaningful name for
  !                 conversion to C-P 'C' dynamics grid. R. Rawlins
  !  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

  ! Mean radius of Earth in m.
  REAL(KIND=r8), PARAMETER  :: Earth_Radius = 6371229.0_r8
  REAL(KIND=r8), PARAMETER :: recip_a2  =1.0_r8/(earth_radius*earth_radius)     ! 1/(radius of earth)^2

  !*L------------------COMDECK C_G----------------------------------------
  ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

  REAL(KIND=r8), PARAMETER :: G = 9.80665_r8
  ! lambda=100000 m
  ! kwv =2*pi/lambda
  ! 
  ! kwv    = 6.28e-5_r8          ! 100 km wave length
  !GWD_FRC_max = 6.0
  !GWD_FRC_min = 2.0
  REAL(KIND=r8), PARAMETER :: Lambda_GWAVE_max = 2.0E+05_r8!highest wave length
  REAL(KIND=r8), PARAMETER :: Lambda_GWAVE_min = 2.5E+02_r8!lowet wave length
  REAL(KIND=r8)  :: kay     ! surface stress constant ( m**-1)

  REAL(KIND=r8) ,PARAMETER  :: frc=5.5_r8!            ! critical Froude number
  LOGICAL        ,PARAMETER :: l_fix_gwsatn=.TRUE.! if true then invoke minor bug fixes in gwsatn
  LOGICAL        ,PARAMETER :: l_gwd_40km=.TRUE.!.FALSE.  ! if true then don't apply GWD above 40km 
  LOGICAL        ,PARAMETER :: l_taus_scale=.FALSE.! if true then surface stress is made to
    !                                            ! depend on the low level Froude number
    ! Number of standard deviations above the mean orography of top
    !  of sub-grid mountains
  REAL(KIND=r8),PARAMETER :: NSIGMA = 3.50_r8
  REAL(KIND=r8) :: akappa

  !*----------------------------------------------------------------------


  !*----------------------------------------------------------------------
   PUBLIC :: Init_Gwave
   PUBLIC :: g_wave
CONTAINS
  SUBROUTINE Init_Gwave(         )
    REAL(KIND=r8) :: ku!is the hisghest wavenumber lying outside the range of trapped lee waves
    REAL(KIND=r8) :: kl!is the lowest unresolved wave number of model
    REAL(KIND=r8) :: pow1
    REAL(KIND=r8) :: pow2
    
    INTEGER :: k
    akappa=gasr/cp
    !
    pow1=1.0_r8/2.0_r8
    pow2=3.0_r8/2.0_r8
    ku=(2.0_r8*Pi)/Lambda_GWAVE_min   !is the hisghest wavenumber lying outside the range of trapped lee waves
    kl=(2.0_r8*Pi)/Lambda_GWAVE_max   !is the lowest unresolved wave number of model
    
    kay=(1.0_r8/3.0_r8)*(((ku**(pow2)) - (kl**(pow2)))/ &
                           ((ku**(pow1)) - (kl**(pow1))))   
   
    !WRITE(0,*)'Init_Gwave',kay
  END SUBROUTINE Init_Gwave


  SUBROUTINE g_wave( &
       kMax         , & !INTEGER, INTENT(IN   ) :: kMax
       nCols        , & !INTEGER, INTENT(IN   ) :: nCols ! number of points per row
         prsi ,prsl  ,phii ,phil    ,&
       gt           , & !REAL(KIND=r8), INTENT(IN        ) :: gt          (nCols,kMax)         !REAL(r8), INTENT(in) :: gt (nCols,kMax) 
       gq           , & !REAL(KIND=r8), INTENT(IN        ) :: gq          (nCols,kMax)         !REAL(r8), INTENT(in) :: gq (nCols,kMax) 
       gu           , & !REAL(KIND=r8), INTENT(IN        ) :: gu          (nCols,kMax)         !REAL(r8), INTENT(in) :: gu (nCols,kMax) 
       gv           , & !REAL(KIND=r8), INTENT(IN        ) :: gv          (nCols,kMax)         !REAL(r8), INTENT(in) :: gv (nCols,kMax) 
       topo         , & !REAL(KIND=r8), INTENT(IN   ) :: topo        (nCols)
       !colrad       , & !REAL(KIND=r8), INTENT(IN   ) :: colrad      (nCols)        
       sd_orog      , & !REAL(KIND=r8) ,INTENT(in   ) :: sd_orog     (nCols)! standard deviation of orography (m)
       orog_grad_xx , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_xx(nCols)! dh/dx squared gradient orography
       orog_grad_yy , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_yy(nCols)! dh/dy squared gradient orography
       orog_grad_xy , & !REAL(KIND=r8), INTENT(in   ) :: orog_grad_xy(nCols)! (dh/dx)(dh/dy) gradient orography
       imask        , & !INTEGER       , INTENT(in   ) :: imask       (nCols)! index for land points
       timestep     , & !REAL(KIND=r8) ,INTENT(in   ) :: timestep            ! timestep (s)
       du_dt            , & !REAL(KIND=r8) ,INTENT(out):: du_dt(nCols,kMax)   ! total GWD du/dt on land/theta
       dv_dt            , & !REAL(KIND=r8) ,INTENT(out):: dv_dt(nCols,kMax)   ! total GWD dv/dt on land/theta
       iret           ) !INTEGER        ,INTENT(OUT  ) :: iret               ! return code : iret=0 normal exit

    IMPLICIT NONE
    !
    ! Description:
    ! 1) Interpolate winds to theta points
    !    gather data for land points only
    ! 2) Call surface stress routine
    ! 3) Calculate stress profiles due to different components of the
    !    scheme. Calculate associated wind increments.
    ! 4) Interpolate acceleration to wind points and update winds
    ! 5) Gather diagnostics from land points and interpolate from p to
    !    u,v staggering on 'c' grid
    !
    ! current code owner: S.Webster
    !
    ! history:
    ! version   date     comment
    ! -------   ----     -------
    !  5.2   15/11/00   original deck.          Stuart Webster
    !  5.3   16/10/01   Remove code no longer required because of
    !                   simplifications to scheme. Stuart Webster
    !  5.3   27/07/01   Permit a LAM model to run with GWD scheme on.
    !                                                      S. Cusack
    !  5.4   28/08/02   Introduce numerical limiter for flow blocking
    !                   scheme.                            S. Webster
    !  6.2   17/05/06   Remove t_to_p reference. P.Selwood.
    !  6.2   21/02/06   Introduce surface stress diagnostics. S.Webster
    !  6.2   21/02/06   Pass l_taus_scale switch thru' to GWSURF4A
    !                   and l_fix_gwsatn thru to GWSATN4A. S. Webster
    !
    !
    !   language: fortran 77 + common extensions.
    !   this code is written to umdp3 v6 programming standards.
    !
    ! suitable for single column use, with calls to: uv_to_p removed
    !                                                p_to_uv removed
    ! suitable for rotated grids
    !
    ! global variables (*called comdecks etc...):
    !
    ! SUBROUTINE ARGUMENTS
    !
    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: nCols ! number of points per row
    REAL(KIND=r8), INTENT(in   ) :: prsi   (ncols,kMax+1)  !     prsi     - real, pressure at layer interfaces [Pa]
    REAL(KIND=r8), INTENT(in   ) :: prsl   (ncols,kMax)    !     prsl     - real, mean layer presure [Pa]
    REAL(KIND=r8), INTENT(in   ) :: phii   (nCols,kMax+1) !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8), INTENT(in   ) :: phil   (nCols,kMax)   !===>  PHIL(K) INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8), INTENT(IN   ) :: gt          (nCols,kMax)        !REAL(r8), INTENT(in) :: gt (nCols,kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: gq          (nCols,kMax)        !REAL(r8), INTENT(in) :: gq (nCols,kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: gu          (nCols,kMax)        !REAL(r8), INTENT(in) :: gu (nCols,kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: gv          (nCols,kMax)        !REAL(r8), INTENT(in) :: gv (nCols,kMax)   
    REAL(KIND=r8), INTENT(IN   ) :: topo        (nCols)
    !REAL(KIND=r8), INTENT(IN   ) :: colrad      (nCols)        
    INTEGER(KIND=i8), INTENT(in   ) :: imask       (nCols)! index for land points
    REAL(KIND=r8), INTENT(in   ) :: orog_grad_xx(nCols)! dh/dx squared gradient orography
    REAL(KIND=r8), INTENT(in   ) :: orog_grad_yy(nCols)! dh/dy squared gradient orography
    REAL(KIND=r8), INTENT(in   ) :: orog_grad_xy(nCols)! (dh/dx)(dh/dy) gradient orography
    REAL(KIND=r8) ,INTENT(in   ) :: sd_orog     (nCols)! standard deviation of orography (m)
    REAL(KIND=r8) ,INTENT(in   ) :: timestep            ! timestep (s)
    INTEGER        ,INTENT(OUT  ) :: iret               ! return code : iret=0 normal exit
    !  Start of full field diagnostic arrays. This space is allocated
    ! in GWD_CTL2 only if a diagnostic is called.
    !
    REAL(KIND=r8) ,INTENT(out):: du_dt(nCols,kMax)   ! total GWD du/dt on land/theta
    REAL(KIND=r8) ,INTENT(out):: dv_dt(nCols,kMax)   ! total GWD dv/dt on land/theta




    REAL(KIND=r8)    :: theta          (nCols,kMax)! primary theta field (K)
    REAL(KIND=r8)    :: rho            (1:nCols,kMax)! density*(radius of earth)^2


    REAL(KIND=r8)    :: r_rho_levels     (1:nCols,  kMax)! height of rho level above earth's centre (m)
    REAL(KIND=r8)    :: r_theta_levels   (1:nCols,0:kMax)! height of theta level above earth's centre (m)
    !
    !
    ! END of full field diagnostic arrays
    ! Below are the stash flags for calculating diagnostics:
    !
    INTEGER :: land_points         ! number of land points
    LOGICAL :: stress_ud_on          !u              stress
    LOGICAL :: stress_ud_p_on    !u on press               stress
    LOGICAL :: stress_vd_on          !v              stress
    LOGICAL :: stress_ud_satn_on !u satn stress
    LOGICAL :: stress_vd_satn_on !v satn stress
    LOGICAL :: stress_ud_wake_on !u wake stress
    LOGICAL :: stress_vd_wake_on !v wake stress
    LOGICAL :: du_dt_satn_on          !u accel (saturation)
    LOGICAL :: du_dt_satn_p_on   !u accel (saturation)
    LOGICAL :: dv_dt_satn_on          !v accel (saturation)
    LOGICAL :: du_dt_wake_on          !u accel blocked flow
    LOGICAL :: dv_dt_wake_on          !v accel blocked flow
    LOGICAL :: u_s_d_on          !u_s_d   diag switch
    LOGICAL :: v_s_d_on          !v_s_d   diag switch
    LOGICAL :: nsq_s_d_on          !nsq_s_d diag switch
    LOGICAL :: fr_d_on           !fr_d    switch
    LOGICAL :: bld_d_on          !bld_d   switch
    LOGICAL :: bldt_d_on          !bldt_d   switch
    LOGICAL :: num_lim_d_on          !num_lim_d switch
    LOGICAL :: num_fac_d_on          !num_fac_d switch
    LOGICAL :: tausx_d_on          !tausx_d switch
    LOGICAL :: tausy_d_on          !tausy_d switch
    LOGICAL :: taus_scale_d_on  !taus_scale_d switch


    !
    ! The integers below are set to size nCols if the corresponding
    ! diagnostic is called or to 1 if it is not. These are set in GWD_CTL2
    !
    INTEGER :: points_stress_ud                                    
    INTEGER :: points_stress_vd                                    
    INTEGER :: points_stress_ud_satn                              
    INTEGER :: points_stress_vd_satn                              
    INTEGER :: points_stress_ud_wake                              
    INTEGER :: points_stress_vd_wake                              
    INTEGER :: points_du_dt_satn                                    
    INTEGER :: points_dv_dt_satn                                    
    INTEGER :: points_du_dt_wake                                    
    INTEGER :: points_dv_dt_wake                                    
    INTEGER :: points_u_s_d                                            
    INTEGER :: points_v_s_d                                            
    INTEGER :: points_nsq_s_d                                     
    INTEGER :: points_fr_d                                            
    INTEGER :: points_bld_d                                            
    INTEGER :: points_bldt_d                                      
    INTEGER :: points_num_lim_d                                    
    INTEGER :: points_num_fac_d                                    
    INTEGER :: points_tausx_d                                     
    INTEGER :: points_tausy_d                                     
    INTEGER :: points_taus_scale_d
    !

    REAL(KIND=r8) :: stress_ud     (nCols    ,0:kMax) !u   total stress
    REAL(KIND=r8) :: stress_vd     (nCols    ,0:kMax) !v   total stress
    REAL(KIND=r8) :: stress_ud_satn(nCols    ,0:kMax)!u   satn  stress
    REAL(KIND=r8) :: stress_vd_satn(nCols    ,0:kMax)!v   satn  stress
    REAL(KIND=r8) :: stress_ud_wake(nCols    ,0:kMax)!u   wake  stress
    REAL(KIND=r8) :: stress_vd_wake(nCols    ,0:kMax)!v   wake  stress
    REAL(KIND=r8) :: du_dt_satn    (nCols    ,kMax)!u acceln (saturation)
    REAL(KIND=r8) :: dv_dt_satn    (nCols    ,kMax)!v acceln (saturation)
    REAL(KIND=r8) :: du_dt_land    (nCols,kMax)   ! total GWD du/dt on land/theta
    REAL(KIND=r8) :: dv_dt_land    (nCols,kMax)   ! total GWD dv/dt on land/theta
    REAL(KIND=r8) :: du_dt_wake    (nCols    ,kMax)!u acceln (blocked flow)
    REAL(KIND=r8) :: dv_dt_wake    (nCols    ,kMax)!v acceln (blocked flow)

    REAL(KIND=r8) :: u_s_d         (nCols  ) ! u_s  diag at theta pts
    REAL(KIND=r8) :: v_s_d         (nCols  ) ! v_s  diag at theta pts
    REAL(KIND=r8) :: nsq_s_d       (nCols  ) ! nsq_s diag at theta pts
    REAL(KIND=r8) :: fr_d          (nCols  )        ! Fr diag at theta pts
    REAL(KIND=r8) :: bld_d         (nCols  )        ! blocked layer depth at theta pts
    REAL(KIND=r8) :: bldt_d        (nCols  )        ! % of time blocked layer diagnosed
    REAL(KIND=r8) :: num_lim_d     (nCols  ) ! % of time numerical limiter invoked
    REAL(KIND=r8) :: num_fac_d     (nCols  ) ! % redn. of flow-blocking stress after numerical limiter invoked
    REAL(KIND=r8) :: tausx_d       (nCols)!x-component of surface stress
    REAL(KIND=r8) :: tausy_d       (nCols)!y-component of surface stress
    REAL(KIND=r8) :: taus_scale_d  (nCols)! Factor surface stress scaled by
    ! if Froude no. dependence is on.

    !      LOGICAL ,intent(in ):: at_extremity(4)    ! indicates if this processor is at north,
    !                                             ! south, east or west of the processor grid

    ! FLDTYPE definitions for the different field types recognised on the
    ! decomposition
    INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
    INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
    INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
    INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
    INTEGER,PARAMETER:: fld_type_comp_wave  = 4
    ! Compressed WAM Wave Field
    INTEGER,PARAMETER:: fld_type_full_wave  = 5
    ! Uncompressed WAM Wave Field
    INTEGER,PARAMETER:: fld_type_rim_wave   = 6
    ! Boundary data for WAM Wave Field
    INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
    INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
    ! FLDTYPE end
    ! DOMTYP contains different model domain types
    !
    ! Author : P.Burton
    ! History:
    ! Version  Date      Comment.
    ! 5.0      15/04/99  New comdeck
    ! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

    INTEGER,PARAMETER:: mt_global        = 1
    INTEGER,PARAMETER:: mt_lam           = 2
    INTEGER,PARAMETER:: mt_cyclic_lam    = 3
    INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
    INTEGER,PARAMETER:: mt_single_column = 5
    ! DOMTYP end

    !--------------------------------------------------------------------
    ! LOCAL DYNAMIC ARRAYS:
    !--------------------------------------------------------------------

    INTEGER :: k_top     (nCols) ! model level at mountain tops -
    !                                 ! exact definition given in gwsurf
    INTEGER :: k_top_max          ! max(k_top)

    ! parameters for mpp

    INTEGER,PARAMETER :: pnorth= 1 ! north processor address in the neighbor array
    INTEGER,PARAMETER :: peast = 2 ! east processor address in the neighbor array
    INTEGER,PARAMETER :: psouth= 3 ! south processor address in the neighbor array
    INTEGER,PARAMETER :: pwest  = 4! west processor address in the neighbor array
    INTEGER,PARAMETER :: nodomain = -1! value in neighbor array if the domain has
    !  no neighbor in this direction. otherwise
    !  the value will be the tid of the neighbor

    INTEGER         :: i,k,points! loop counters in routine

    ! Work arrays

    REAL(KIND=r8) :: work_u           (nCols   ,kMax)
    REAL(KIND=r8) :: work_v           (nCols   ,kMax)
    !                                  REAL :: work_halo      (1:nCols,0:kMax)
    !                                  REAL :: work_on_v_grid  (nCols,kMax)

    REAL(KIND=r8) :: up_land          (nCols,kMax)! interpolated u on theta grid
    REAL(KIND=r8) :: vp_land          (nCols,kMax)! interpolated V on theta grid
    REAL(KIND=r8) :: theta_land       (nCols,kMax)! land theta field on theta levels
    REAL(KIND=r8) :: r_rho_levels_land(nCols,kMax)!  land field of heights of
    !                                    !  rho levels above z=0

    REAL(KIND=r8) :: r_theta_levels_land(nCols,kMax)!  land field of heights of
    !                                    !  theta levels above z=0

    REAL(KIND=r8) :: rho_land        (nCols,kMax)  ! density at land points

    REAL(KIND=r8) :: s_x_lin_stress   (nCols)   ! 'surface'  x_lin_stress land pnts
    REAL(KIND=r8) :: s_y_lin_stress   (nCols)   ! 'surface'  y_lin_stress land pnts
    REAL(KIND=r8) :: s_x_wake_stress  (nCols)  ! 'surface' x_wake_stress land pts
    REAL(KIND=r8) :: s_y_wake_stress  (nCols)  ! 'surface' y_wake_stress land pts
    REAL(KIND=r8) :: s_x_orog         (nCols)  ! 'surface' x_orog on land points
    REAL(KIND=r8) :: s_y_orog         (nCols)         ! 'surface' y_orog on land points
    REAL(KIND=r8) :: lift             (nCols)             ! depth of blocked layer
    REAL(KIND=r8) :: fr               (nCols)               ! low level froude number
    REAL(KIND=r8) :: rho_s            (nCols)            ! low level density

    REAL(KIND=r8) :: orog_grad_xx_land(nCols)! dh/dx squared gradient orography
    REAL(KIND=r8) :: orog_grad_yy_land(nCols)! dh/dy squared gradient orography
    REAL(KIND=r8) :: orog_grad_xy_land(nCols)! (dh/dx)(dh/dy) gradient orography
    REAL(KIND=r8) :: sd_orog_land     (nCols)! standard deviation of orography (m)

    !
    ! Land points arrays below are for the total GWD stress and
    ! its 4 individual components. (x and y components for each)
    !
    REAL(KIND=r8) :: stress_ud_land          ( nCols , 0:kMax )   
    REAL(KIND=r8) :: stress_vd_land          ( nCols , 0:kMax )   
    REAL(KIND=r8) :: stress_ud_satn_land     ( nCols , 0:kMax )   
    REAL(KIND=r8) :: stress_vd_satn_land     ( nCols , 0:kMax )   
    REAL(KIND=r8) :: stress_ud_wake_land     ( nCols , 0:kMax )   
    REAL(KIND=r8) :: stress_vd_wake_land     ( nCols , 0:kMax )
    !
    !  Land point arrays below are for the 4 individual components of
    !  the GWD wind increment.  (x and y components for each)
    !
    REAL(KIND=r8) :: du_dt_satn_land ( nCols , kMax )   
    REAL(KIND=r8) :: dv_dt_satn_land ( nCols , kMax )   
    REAL(KIND=r8) :: du_dt_wake_land ( nCols , kMax )   
    REAL(KIND=r8) :: dv_dt_wake_land ( nCols , kMax )

    !
    !  Land point arrays below are for the 9 GWD 'surface' diagnostics.
    !
    REAL(KIND=r8) :: u_s_d_land     ( nCols        )
    REAL(KIND=r8) :: v_s_d_land     ( nCols        )
    REAL(KIND=r8) :: nsq_s_d_land   ( nCols        )
    REAL(KIND=r8) :: fr_d_land      ( nCols        )
    REAL(KIND=r8) :: bld_d_land     ( nCols        )
    REAL(KIND=r8) :: bldt_d_land    ( nCols        )
    REAL(KIND=r8) :: num_lim_d_land ( nCols      )
    REAL(KIND=r8) :: num_fac_d_land ( nCols      )
    REAL(KIND=r8) :: tausx_d_land   ( nCols        )
    REAL(KIND=r8) :: tausy_d_land   ( nCols        )
    REAL(KIND=r8) :: taus_scale_d_land ( nCols )
    LOGICAL :: l_drag(nCols)           ! whether point has a non-zero stress or not
    REAL(KIND=r8) :: pmid   (nCols,kMax)
    REAL(KIND=r8) :: state_t     (nCols,kMax)  
    REAL(KIND=r8) :: state_q     (nCols,kMax)  
    REAL(KIND=r8) :: state_u     (nCols,kMax)  
    REAL(KIND=r8) :: state_v     (nCols,kMax)  
    REAL(KIND=r8) :: state_pmid  (nCols,kMax)  
    REAL(KIND=r8) :: state_pint  (nCols,kMax+1)  
    REAL(KIND=r8) :: state_lnpint(nCols,kMax+1)  
    REAL(KIND=r8) :: state_pdel  (nCols,kMax)  
    REAL(KIND=r8) :: state_rpdel (nCols,kMax)
    REAL(KIND=r8) :: state_lnpmid(nCols,kMax)
    REAL(KIND=r8) :: state_zi   (1:nCols,1:kMax+1)   
    REAL(KIND=r8) :: state_zm   (1:nCols,1:kMax)
    REAL(KIND=r8) :: sigkiv     (1:nCols,1:kMax)

    k_top=0
    k_top_max=0 
    theta         =0.0_r8
    rho           =0.0_r8     
    r_rho_levels  =0.0_r8 
    r_theta_levels=0.0_r8
    stress_ud     (1:nCols    ,0:kMax)=0.0_r8
    stress_vd     (1:nCols    ,0:kMax)=0.0_r8
    stress_ud_satn(1:nCols    ,0:kMax)=0.0_r8
    stress_vd_satn(1:nCols    ,0:kMax)=0.0_r8
    stress_ud_wake(1:nCols    ,0:kMax)=0.0_r8
    stress_vd_wake(1:nCols    ,0:kMax)=0.0_r8
    du_dt_satn =0.0_r8
    dv_dt_satn =0.0_r8
    du_dt_land =0.0_r8
    dv_dt_land =0.0_r8
    du_dt_wake =0.0_r8
    dv_dt_wake =0.0_r8
    u_s_d =0.0_r8
    v_s_d =0.0_r8
    nsq_s_d =0.0_r8
    fr_d  =0.0_r8
    bld_d =0.0_r8
    bldt_d =0.0_r8
    num_lim_d =0.0_r8
    num_fac_d =0.0_r8
    tausx_d =0.0_r8
    tausy_d =0.0_r8
    taus_scale_d  =0.0_r8
    r_theta_levels_land=0.0_r8
    s_x_lin_stress   =0.0_r8
    s_y_lin_stress   =0.0_r8
    s_x_wake_stress  =0.0_r8
    s_y_wake_stress  =0.0_r8
    s_x_orog=0.0_r8
    s_y_orog=0.0_r8
    lift=0.0_r8
    fr=0.0_r8
    rho_s=0.0_r8
    orog_grad_xx_land=0.0_r8
    orog_grad_yy_land=0.0_r8
    orog_grad_xy_land=0.0_r8
    sd_orog_land=0.0_r8
    stress_ud_land     (1: nCols , 0:kMax )= 0.0_r8   
    stress_vd_land     (1: nCols , 0:kMax )= 0.0_r8    
    stress_ud_satn_land(1: nCols , 0:kMax )= 0.0_r8
    stress_vd_satn_land(1: nCols , 0:kMax )= 0.0_r8
    stress_ud_wake_land(1: nCols , 0:kMax )= 0.0_r8
    stress_vd_wake_land(1: nCols , 0:kMax )= 0.0_r8
    du_dt_satn_land=0.0_r8
    dv_dt_satn_land=0.0_r8
    du_dt_wake_land=0.0_r8
    dv_dt_wake_land=0.0_r8
    
    taus_scale_d_land=0.0_r8
    
    land_points=0
    DO i=1,ncols
       IF(imask(i).GE.1_i8) THEN
          land_points=land_points+1
       END IF
    END DO

    !
    ! END of full field diagnostic arrays
    ! Below are the stash flags for calculating diagnostics:
    !
    stress_ud_on=.TRUE.                !u stress
    points_stress_ud=land_points
    stress_ud     (1:land_points   ,0:kMax)=0.0_r8 !u total stress


    stress_ud_p_on=.TRUE.        !u on press             stress

    stress_vd_on=.TRUE.        !v            stress
    points_stress_vd=land_points
    stress_vd     (1:land_points   ,0:kMax)=0.0_r8 !v   total stress

    stress_ud_satn_on=.TRUE. !u satn stress
    points_stress_ud_satn=land_points
    stress_ud_satn(1:land_points    ,0:kMax)=0.0_r8!u   satn  stress

    stress_vd_satn_on=.TRUE. !v satn stress
    points_stress_vd_satn=land_points
    stress_vd_satn(1:land_points    ,0:kMax)=0.0_r8!v   satn  stress

    stress_ud_wake_on=.TRUE. !u wake stress
    points_stress_ud_wake=land_points
    stress_ud_wake(1:land_points    ,0:kMax)=0.0_r8!u   wake  stress

    stress_vd_wake_on=.TRUE. !v wake stress
    points_stress_vd_wake=land_points
    stress_vd_wake(1:land_points    ,0:kMax)=0.0_r8!v   wake  stress

    du_dt_satn_on=.TRUE.        !u accel (saturation)
    points_du_dt_satn=land_points
    du_dt_satn    (1:land_points    ,1:kMax)=0.0_r8!u acceln (saturation)

    du_dt_satn_p_on=.TRUE.        !u accel (saturation)

    dv_dt_satn_on=.TRUE.        !v accel (saturation)
    points_dv_dt_satn=land_points
    dv_dt_satn    (1:land_points    ,1:kMax)=0.0_r8!v acceln (saturation)

    du_dt_wake_on=.TRUE.        !u accel blocked flow
    points_du_dt_wake=land_points
    du_dt_wake    (1:land_points    ,1:kMax)=0.0_r8!u acceln (blocked flow)

    dv_dt_wake_on=.TRUE.        !v accel blocked flow
    points_dv_dt_wake=land_points
    dv_dt_wake    (1:land_points    ,1:kMax)=0.0_r8!v acceln (blocked flow)

    u_s_d_on=.TRUE.          !u_s_d   diag switch
    points_u_s_d=land_points
    u_s_d         (1:land_points  )=0.0_r8 ! u_s  diag at theta pts

    v_s_d_on=.TRUE.          !v_s_d   diag switch
    points_v_s_d=land_points
    v_s_d         (1:land_points  ) =0.0_r8! v_s  diag at theta pts

    nsq_s_d_on=.TRUE.        !nsq_s_d diag switch
    points_nsq_s_d=land_points
    nsq_s_d       (1:land_points  ) =0.0_r8! nsq_s diag at theta pts

    fr_d_on=.TRUE.                !fr_d         switch
    points_fr_d=land_points
    fr_d          (1:land_points  ) =0.0_r8   ! Fr diag at theta pts

    bld_d_on=.TRUE.          !bld_d   switch
    points_bld_d=land_points
    bld_d         (1:land_points  )=0.0_r8    ! blocked layer depth at theta pts

    bldt_d_on=.TRUE.         !bldt_d   switch
    points_bldt_d=land_points
    bldt_d        (1:land_points  )=0.0_r8    ! % of time blocked layer diagnosed

    num_lim_d_on=.TRUE.        !num_lim_d switch
    points_num_lim_d=land_points
    num_lim_d     (1:land_points  )=0.0_r8 ! % of time numerical limiter invoked

    num_fac_d_on=.TRUE.        !num_fac_d switch
    points_num_fac_d=land_points
    num_fac_d     (1:land_points  )=0.0_r8 ! % redn. of flow-blocking stress
    ! after numerical limiter invoked
    tausx_d_on=.TRUE.        !tausx_d switch
    points_tausx_d=land_points
    tausx_d     (1:land_points)=0.0_r8!x-component of surface stress

    tausy_d_on=.TRUE.        !tausy_d switch
    points_tausy_d=land_points
    tausy_d     (1:land_points)=0.0_r8!y-component of surface stress

    taus_scale_d_on=.TRUE.        !taus_scale_d switch
    points_taus_scale_d=land_points
    taus_scale_d(1:land_points)=0.0_r8! Factor surface stress scaled by
    ! if Froude no. dependence is on.

    !-----------------------------------------------------------------------------
    DO i=1,nCols
       !state_pint       (i,kMax+1) = gps(i)*si(1)
       state_pint       (i,kMax+1) = prsi(i,1)
    END DO
    DO k=kMax+1,1,-1
       DO i=1,nCols
          !state_pint    (i,k)      = MAX(si(kMax+2-k)*gps(i) ,0.0001_r8)
          state_pint    (i,k)      = MAX(prsi(i,kMax+2-k),0.0001_r8)
       END DO
    END DO


    DO k=1,kMax
       DO i=1,nCols
          state_t (i,kMax+1-k) =  gt (i,k)
          state_q (i,kMax+1-k) =  gq (i,k)
          state_u (i,kMax+1-k) =  gu (i,k)
          state_v (i,kMax+1-k) =  gv (i,k)
          !state_pmid(i,kMax+1-k) = sl(k)*gps (i)
          state_pmid(i,kMax+1-k) = prsl(i,k)
       END DO
    END DO

    DO k=1,kMax
       DO i=1,nCols    
          state_pdel    (i,k) = MAX(state_pint(i,k+1) - state_pint(i,k),0.000000005_r8)
          state_rpdel   (i,k) = 1.0_r8/MAX((state_pint(i,k+1) - state_pint(i,k)),0.00000000005_r8)
          state_lnpmid  (i,k) = LOG(state_pmid(i,k))        
       END DO
    END DO
    DO k=1,kMax+1
       DO i=1,nCols
           state_lnpint(i,k) =  LOG(state_pint  (i,k))
       END DO
    END DO


    ! Derive new temperature and geopotential fields

    CALL geopotential_t(                                 &
         state_lnpint(1:nCols,1:kMax+1)   , state_pint (1:nCols,1:kMax+1)  , &
         state_pmid  (1:nCols,1:kMax)     , state_pdel  (1:nCols,1:kMax)   , state_rpdel(1:nCols,1:kMax)   , &
         state_t     (1:nCols,1:kMax)     , state_q     (1:nCols,1:kMax)   , rair   , gravit , zvir   ,&
         state_zi    (1:nCols,1:kMax+1)   , state_zm    (1:nCols,1:kMax)   , nCols   ,nCols, kMax)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Parei aqui
    !
    DO i=1,nCols
       R_THETA_LEVELS(i,0) = EARTH_RADIUS + MAX(topo(i),0.0_r8)+ state_zi(i,kMax+1)
    END DO

    DO k = 1, kmax
       DO i = 1, ncols
          !             --   --  -(R/Cp)        --   --  -(R/Cp)     --   --   (R/Cp)  
          !            |  P    |               |       |            |  P0   |       
          !sigkiv(k) = |-------|            == |sl(k)  |         => |-------|       
          !            |  P0   |               |       |            |  P    |       
          !             --   --                 --   --              --   --       

          !Pbl_PotTep(i,k)=PK_sigkiv(i,k)*gt(i,k)*(1.0_r8+eps*gq(i,k))
          !Pbl_PotTep(i,k)=*gt(i,k)*(1.0_r8+eps*gq(i,k))
          sigkiv(i,k)=(1.0_r8/((prsl(i,k)/prsi(i,1))**akappa))
       END DO
    END DO

    DO k=1,kMax
       DO i=1,nCols
          work_u               (i,k) = gu    (i,k)
          work_v               (i,k) = gv    (i,k)
          R_THETA_LEVELS       (i,k) = (EARTH_RADIUS + MAX(topo(i),0.0_r8)) + state_zi    (i,kMax+1-k)
          r_rho_levels         (i,k) = (EARTH_RADIUS + MAX(topo(i),0.0_r8)) + state_zm    (i,kMax+1-k)
          THETA                (i,k) = sigkiv(i,k)*gt(i,k)
          !pmid                 (i,k) = gps(i)*sl(k)
          pmid                 (i,k) = prsl(i,k)

          !PRINT*,R_THETA_LEVELS       (i,k),state_zi    (i,kMax+1-k)

          !j/kg/kelvin
          !
          ! P = rho * R * T
          !
          !            P
          ! rho  = -------
          !          R * T
          !
          !           1             R * T
          !  rrho = ----- =        -------
          !          rho              P
          !
          !
          rho(i,k) = (pmid(i,k)/(gasr*gt(i,k)))*(earth_radius*earth_radius)! density*(radius of earth)^2
       END DO
    END DO

    !------------------------------------------------------------------
    !l    1.2  gather winds at land points
    !------------------------------------------------------------------

    DO k=1,kMax
       land_points=0
       DO i=1,ncols
          IF(imask(i).GE.1_i8) THEN
             land_points=land_points+1
             !------------------------------------------------------------------
             !l    1.2  gather winds at land points
             !------------------------------------------------------------------
             up_land(land_points,k) =work_u(i,k)
             vp_land(land_points,k) =work_v(i,k)
             !------------------------------------------------------------------
             !l    1.3  gather theta, rho and heights at land points
             !------------------------------------------------------------------

             r_rho_levels_land  (land_points,k)= r_rho_levels  (i,k) - r_theta_levels(i,0)
             r_theta_levels_land(land_points,k)= r_theta_levels(i,k) - r_theta_levels(i,0)
             rho_land           (land_points,k)= rho  (i,k)*recip_a2
             theta_land         (land_points,k)= theta(i,k)
          END IF
       END DO
    END DO

    land_points=0
    DO i=1,ncols
       IF(imask(i).GE.1_i8) THEN
          land_points=land_points+1
          orog_grad_xx_land(land_points)  = orog_grad_xx(i)
          orog_grad_yy_land(land_points)  = orog_grad_yy(i)
          orog_grad_xy_land(land_points)  = orog_grad_xy(i)
          sd_orog_land     (land_points)  = sd_orog     (i)        
       END IF
    END DO
    !------------------------------------------------------------------
    !l    2. calculate anisotropic 'surface' stress,CALL gw_surf
    !------------------------------------------------------------------
   ! PRINT*,'first gw_surf ' 
    IF (land_points  >   0) THEN
       points=land_points
       ! DEPENDS ON: gw_surf
       CALL gw_surf(&
            r_theta_levels_land(1:points,1:kMax), &!REAL,intent(in):: r_theta_levels  (points,kMax)! heights on theta levels
            rho_land           (1:points,1:kMax), &!REAL,intent(in):: rho                  (points,kMax)! density on rho levels
            theta_land         (1:points,1:kMax), &!REAL,intent(in):: theta           (points,kMax)! theta on theta levels
            up_land            (1:points,1:kMax), &!REAL,intent(in):: u                  (points,kMax)! u on rho levels
            vp_land            (1:points,1:kMax), &!REAL,intent(in):: v                  (points,kMax)! v on rho levels
            timestep                                   , &!REAL,intent(in):: timestep                   ! timestep (s)
            sd_orog_land       (1:points)       , &!REAL,intent(in):: sd_orog         (points)! standard deviation of the sub-grid orography
            orog_grad_xx_land  (1:points)       , &!REAL,intent(in):: sigma_xx        (points)! (dh/dx)^2 grid box average of the
            orog_grad_xy_land  (1:points)       , &!REAL,intent(in):: sigma_xy        (points)! (dh/dx)*(dh/dy)
            orog_grad_yy_land  (1:points)       , &!REAL,intent(in):: sigma_yy        (points)! (dh/dy)^2
            s_x_lin_stress     (1:points)       , &!REAL,intent(out):: s_x_lin_stress (points)! 'surface' linear stress in x-dirn
            s_y_lin_stress     (1:points)       , &!REAL,intent(out):: s_y_lin_stress (points)! 'surface' linear stress in y-dirn
            s_x_wake_stress    (1:points)       , &!REAL,intent(out):: s_x_wake_stress(points)! wake stress in x-dirn
            s_y_wake_stress    (1:points)       , &!REAL,intent(out):: s_y_wake_stress(points)! wake stress in y-dirn
            s_x_orog           (1:points)       , &!REAL,intent(out):: s_x_orog       (points)! 'surface' stress/orog in x-dirn
            s_y_orog           (1:points)       , &!REAL,intent(out):: s_y_orog       (points)! 'surface' stress/orog in y-dirn
            kMax                                , &!INTEGER,intent(in) :: kMax ! number of model levels
            land_points                         , &!INTEGER,intent(in) :: points ! number of gwd points on pe
            kay                                 , &!REAL ,intent(in):: kay                    !  gwd surface stress constant (m-1)
            rho_s              (1:points)       , &!REAL,intent(out):: rho_s(points)          ! density - av from z=0 to nsigma*sd_orog
            l_taus_scale                        , &!LOGICAL,intent(in):: l_taus_scale            ! true allows variation of surface stress on the low level Froude number
            k_top              (1:points)       , &!INTEGER ,intent(out):: k_top(points)              ! topmost model level including mountains.
            k_top_max                           , &!INTEGER ,intent(out):: k_top_max              ! max(k_top)
            lift               (1:points)       , &! REAL,intent(out):: lift(points)           ! blocked layer depth
            l_drag             (1:points)       , &!LOGICAL,intent(out):: l_drag(points)             ! whether a point has a non-zero surface stress or not
            fr                 (1:points)       , &!REAL,intent(out):: fr(points)                  ! low level froude number for linear
            frc                                 , &!REAL ,intent(in):: frc                    !  critical froude number below which hydraulic jumps are triggered
            u_s_d_land         (1:points)       , &!REAL,intent(out):: u_s_d(points_u_s_d)    ! 0-nsigma*sd_orog u_s diag
            u_s_d_on                            , &!LOGICAL,intent(in):: u_s_d_on
            points_u_s_d                        , &!INTEGER ,intent(in) ::points_u_s_d
            v_s_d_land         (1:points)       , &!REAL,intent(out):: v_s_d(points_v_s_d)    ! 0-nsigma_sd_orog v_s diag
            v_s_d_on                            , &!LOGICAL,intent(in):: v_s_d_on  
            points_v_s_d                        , &!INTEGER ,intent(in) ::points_v_s_d
            nsq_s_d_land       (1:points)       , &!REAL,intent(out):: nsq_s_d(points_nsq_s_d)! 0-nsigma*sd_orog nsq_s diagnostic
            nsq_s_d_on                          , &!LOGICAL,intent(in):: nsq_s_d_on          
            points_nsq_s_d                      , &!INTEGER ,intent(in) ::points_nsq_s_d            
            fr_d_land          (1:points)       , &!REAL,intent(out):: fr_d(points_fr_d)          ! Froude no. diagnostic
            fr_d_on                             , &!LOGICAL,intent(in):: fr_d_on         
            points_fr_d                                , &!INTEGER ,intent(in) ::points_fr_d
            bld_d_land         (1:points)       , &!REAL,intent(out):: bld_d(points_bld_d)    ! blocked layer depth diagnostic
            bld_d_on                            , &!LOGICAL,intent(in):: bld_d_on
            points_bld_d                        , &!INTEGER ,intent(in) ::points_bld_d
            bldt_d_land        (1:points)       , &!REAL,intent(out):: bldt_d(points_bldt_d)  ! %  of time blocked flow param. invoked
            bldt_d_on                           , &!LOGICAL,intent(in):: bldt_d_on
            points_bldt_d                       , &!INTEGER ,intent(in) ::points_bldt_d
            num_lim_d_land     (1:points)       , &!REAL,intent(out):: num_lim_d(points_num_lim_d)! % of time numerical limiter invoked
            num_lim_d_on                        , &!LOGICAL,intent(in):: num_lim_d_on        
            points_num_lim_d                    , &!INTEGER ,intent(in) ::points_num_lim_d
            num_fac_d_land     (1:points)       , &!REAL,intent(out):: num_fac_d(points_num_fac_d)! % reduction of flow blocking stressafter numerical limiter invoked
            num_fac_d_on                        , &!LOGICAL,intent(in):: num_fac_d_on            
            points_num_fac_d                    , &!INTEGER ,intent(in) ::points_num_fac_d
            tausx_d_land       (1:points)       , &!REAL,intent(out):: tausx_d(points_tausx_d) ! x-component of total surface stress
            tausx_d_on                          , &!LOGICAL,intent(in):: tausx_d_on      
            points_tausx_d                      , &!INTEGER ,intent(in) ::points_tausx_d
            tausy_d_land       (1:points)       , &!REAL,intent(out):: tausy_d(points_tausy_d) ! y-component of total surface stress
            tausy_d_on                          , &!LOGICAL,intent(in):: tausy_d_on             
            points_tausy_d                      , &!INTEGER ,intent(in) ::points_tausy_d
            taus_scale_d_land  (1:points)       , &!REAL,intent(out):: taus_scale_d(points_taus_scale_d)! scaling factor for surface stress when Fr dependence of surface stress invoked.
            taus_scale_d_on                     , &!LOGICAL,intent(in):: taus_scale_d_on
            points_taus_scale_d                   )!INTEGER ,intent(in) ::points_taus_scale_d

!PRINT*,'end gw_surf '  
       !------------------------------------------------------------------
       !l    3. calculate stress profile and accelerations,
       !l       CALL gw_vert
       !------------------------------------------------------------------
!PRINT*,'first gw_vert ' 
       ! DEPENDS ON: gw_vert
       CALL gw_vert(            &
            rho_land             (1:points,1:kMax), &!REAL,intent(in):: rho           (points,kMax)! density on rho levels
            r_rho_levels_land    (1:points,1:kMax), &!REAL,intent(in):: r_rho_levels  (points,kMax)! heights on rho levels
            r_theta_levels_land  (1:points,1:kMax), &!REAL,intent(in):: r_theta_levels(points,kMax)! heights on theta levels
            theta_land           (1:points,1:kMax), &!REAL,intent(in):: theta         (points,kMax)! theta field
            up_land              (1:points,1:kMax), &!REAL,intent(in):: u             (points,kMax)! u field
            vp_land              (1:points,1:kMax), &!REAL,intent(in):: v             (points,kMax)! v field
            kMax                                  , &!INTEGER ,intent(in) :: kMax                     ! number of model levels
            land_points                           , &!INTEGER ,intent(in) :: points                     ! number of points
            kay                                   , &!REAL,intent(in):: kay                        ! GWD hydrostatic constant
            sd_orog_land         (1:points)       , &!REAL,intent(in):: sd_orog        (points)  ! standard deviation of the sub-grid orography
            s_x_lin_stress       (1:points)       , &!REAL,intent(in):: s_x_lin_stress (points) ! 'surface' lin  x_stress
            s_y_lin_stress       (1:points)       , &!REAL,intent(in):: s_y_lin_stress (points)! 'surface' lin  y_stress
            s_x_wake_stress      (1:points)       , &!REAL,intent(in):: s_x_wake_stress(points)! 'surface' x_wake_stress
            s_y_wake_stress      (1:points)       , &!REAL,intent(in):: s_y_wake_stress(points)! 'surface' y_wake_stress
            s_x_orog             (1:points)       , &!REAL,intent(in):: s_x_orog   (points) ! 'surface' x_stress term
            s_y_orog             (1:points)       , &!REAL,intent(in):: s_y_orog   (points) ! 'surface' y_stress term
            du_dt_land           (1:points,1:kMax), &!REAL,intent(out):: du_dt     (points,kMax) ! total GWD du/dt
            dv_dt_land           (1:points,1:kMax), &!REAL,intent(out):: dv_dt     (points,kMax)  ! total GWD dv/dt
            k_top                (1:points)       , &!INTEGER,intent(in):: k_top   (points)         ! model level at mountain tops
            k_top_max                             , &!INTEGER,intent(in):: k_top_max! max(k_top)
            !lift                 (1:points)       , &!REAL   ,intent(in):: lift    (points) ! blocked layer depth
            l_drag               (1:points)       , &!LOGICAL,intent(in):: l_drag  (points)! true if a non-zero surface stress
            fr                   (1:points)       , &!REAL   ,intent(in):: fr      (points)! low level froude number
            rho_s                (1:points)       , &!REAL   ,intent(in):: rho_s   (points)! Low level density calculated in gwsurf
            l_fix_gwsatn                          , &!LOGICAL,intent(in):: l_fix_gwsatn ! switch to include minor bug fixes  
            l_gwd_40km                            , &!LOGICAL,intent(in):: l_gwd_40km     ! switch to turn off GWD above 40km
            ! diagnostics
            stress_ud_land       (1:points,0:kMax), &!REAL   ,intent(out) :: stress_ud         (points_stress_ud,0:kMax)! x total stress diag
            points_stress_ud                      , &!INTEGER,intent(in ) :: points_stress_ud
            stress_ud_on                          , &!LOGICAL,intent(in ) :: stress_ud_on 
            stress_ud_p_on                        , &!LOGICAL,intent(in ) :: stress_ud_p_on   
            stress_vd_land       (1:points,0:kMax), &!REAL   ,intent(out) :: stress_vd         (points_stress_vd,0:kMax)! y total stress diag
            points_stress_vd                      , &!INTEGER,intent(in ) :: points_stress_vd
            stress_vd_on                          , &!LOGICAL,intent(in ) :: stress_vd_on
            stress_ud_satn_land  (1:points,0:kMax), &!REAL   ,intent(out) :: stress_ud_satn   (points_stress_ud_satn,0:kMax)! x saturation stress diag
            points_stress_ud_satn                 , &!INTEGER,intent(in ) :: points_stress_ud_satn
            stress_ud_satn_on                     , &!LOGICAL,intent(in ) :: stress_ud_satn_on
            stress_vd_satn_land  (1:points,0:kMax), &!REAL   ,intent(out) :: stress_vd_satn   (points_stress_vd_satn,0:kMax)! y saturation stress diag
            points_stress_vd_satn                 , &!INTEGER,intent(in ) :: points_stress_vd_satn        
            stress_vd_satn_on                     , &!LOGICAL,intent(in ) :: stress_vd_satn_on
            stress_ud_wake_land  (1:points,0:kMax), &!REAL   ,intent(out) :: stress_ud_wake   (points_stress_ud_wake,0:kMax)! x blocked flow stress diag
            points_stress_ud_wake                 , &!INTEGER,intent(in ) :: points_stress_ud_wake
            stress_ud_wake_on                     , &!LOGICAL,intent(in ) :: stress_ud_wake_on
            stress_vd_wake_land  (1:points,0:kMax), &!REAL   ,intent(out) :: stress_vd_wake   (points_stress_vd_wake,0:kMax)! y blocked flow stress diag
            points_stress_vd_wake                 , &!INTEGER,intent(in ) :: points_stress_vd_wake    
            stress_vd_wake_on                     , &!LOGICAL,intent(in ) :: stress_vd_wake_on
            du_dt_satn_land      (1:points,1:kMax), &!REAL   ,intent(out) :: du_dt_satn       (points_du_dt_satn,kMax)! du/dt diagnostic (saturation)
            points_du_dt_satn                     , &!INTEGER,intent(in ) :: points_du_dt_satn    
            du_dt_satn_on                         , &!LOGICAL,intent(in ) :: du_dt_satn_on
            du_dt_satn_p_on                       , &!LOGICAL,intent(in ) :: du_dt_satn_p_on
            dv_dt_satn_land      (1:points,1:kMax), &!REAL   ,intent(out) :: dv_dt_satn       (points_dv_dt_satn,kMax)! dv/dt diagnostic (saturation)
            points_dv_dt_satn                     , &!INTEGER,intent(in ) :: points_dv_dt_satn
            dv_dt_satn_on                         , &!LOGICAL,intent(in ) :: dv_dt_satn_on
            du_dt_wake_land      (1:points,1:kMax), &!REAL   ,intent(out) :: du_dt_wake       (points_du_dt_wake,kMax)! du/dt diagnostic (blocked flow)
            points_du_dt_wake                     , &!INTEGER,intent(in ) :: points_du_dt_wake
            du_dt_wake_on                         , &!LOGICAL,intent(in ) :: du_dt_wake_on
            dv_dt_wake_land      (1:points,1:kMax), &!REAL   ,intent(out) :: dv_dt_wake       (points_dv_dt_wake,kMax)! dv/dt diagnostic (blocked flow)
            points_dv_dt_wake                     , &!INTEGER,intent(in ) :: points_dv_dt_wake
            dv_dt_wake_on                           )!LOGICAL,intent(in ) :: dv_dt_wake_on

!PRINT*,'end gw_vert ' 
       DO k=1,kMax
          land_points=0
          DO i=1,ncols
             IF(imask(i).GE.1_i8) THEN
                land_points=land_points+1
                !------------------------------------------------------------------
                !l    1.2  gather winds at land points
                !------------------------------------------------------------------
                du_dt(i,k) =du_dt_land(land_points,k)
                dv_dt(i,k) =dv_dt_land(land_points,k)
              !PRINT*,k,du_dt(i,k),dv_dt(i,k)

             END IF
          END DO
       END DO
    END IF ! on land_points > 0


    iret=0

  END SUBROUTINE g_wave
  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  ! subroutine gw_vert to calculate the vertical profile of gravity wave
  !            stress and associated wind increments
  !
  SUBROUTINE gw_vert(       &
       rho                   , &!REAL,intent(in):: rho(points,kMax)! density on rho levels
       r_rho_levels          , &!REAL,intent(in):: r_rho_levels(points,kMax)! heights on rho levels
       r_theta_levels        , &!REAL,intent(in):: r_theta_levels(points,kMax)! heights on theta levels
       theta                 , &!REAL,intent(in):: theta(points,kMax)  ! theta field
       u                     , &!REAL,intent(in):: u(points,kMax)! u field
       v                     , &!REAL,intent(in):: v(points,kMax)! v field
       kMax                  , &!INTEGER ,intent(in) :: kMax                ! number of model levels
       points                , &!INTEGER ,intent(in) :: points                ! number of points
       kay                   , &!REAL,intent(in):: kay                   ! GWD hydrostatic constant
       sd_orog                    , &!REAL,intent(in):: sd_orog(points)  ! standard deviation of the sub-grid orography
       s_x_lin_stress        , &!REAL,intent(in):: s_x_lin_stress(points) ! 'surface' lin  x_stress
       s_y_lin_stress        , &!REAL,intent(in):: s_y_lin_stress(points)! 'surface' lin  y_stress
       s_x_wake_stress       , &!REAL,intent(in):: s_x_wake_stress(points)! 'surface' x_wake_stress
       s_y_wake_stress       , &!REAL,intent(in):: s_y_wake_stress(points)! 'surface' y_wake_stress
       s_x_orog              , &!REAL,intent(in):: s_x_orog(points) ! 'surface' x_stress term
       s_y_orog              , &!REAL,intent(in):: s_y_orog(points) ! 'surface' y_stress term
       du_dt                 , &!REAL,intent(out):: du_dt(points,kMax) ! total GWD du/dt
       dv_dt                 , &!REAL,intent(out):: dv_dt(points,kMax)  ! total GWD dv/dt
       k_top                 , &!INTEGER ,intent(in) :: k_top(points)         ! model level at mountain tops
       k_top_max             , &!INTEGER ,intent(in) :: k_top_max             ! max(k_top)
       !lift                  , &!REAL,intent(in):: lift    (points) ! blocked layer depth
       l_drag                , &!LOGICAL,intent(in):: l_drag(points)        ! true if a non-zero surface stress
       fr                    , &!REAL,intent(in):: fr(points)            ! low level froude number
       rho_s                 , &!REAL,intent(in):: rho_s(points)         ! Low level density calculated in gwsurf
       l_fix_gwsatn          , &!LOGICAL,intent(in):: l_fix_gwsatn ! switch to include minor bug fixes  
       l_gwd_40km            , &!LOGICAL,intent(in):: l_gwd_40km     ! switch to turn off GWD above 40km
       ! diagnostics
         stress_ud             , &!REAL,intent(out):: stress_ud     (points_stress_ud,0:kMax)! x total stress diag
         points_stress_ud      , &!INTEGER,intent(in):: points_stress_ud
         stress_ud_on          , &!LOGICAL,intent(in) :: stress_ud_on 
         stress_ud_p_on        , &!LOGICAL,intent(in) :: stress_ud_p_on   
         stress_vd             , &!REAL,intent(out):: stress_vd     (points_stress_vd,0:kMax)! y total stress diag
         points_stress_vd      , &!INTEGER,intent(in):: points_stress_vd
         stress_vd_on          , &!LOGICAL,intent(in) :: stress_vd_on
         stress_ud_satn        , &!REAL,intent(out):: stress_ud_satn(points_stress_ud_satn,0:kMax)! x saturation stress diag
         points_stress_ud_satn , &!INTEGER,intent(in):: points_stress_ud_satn
         stress_ud_satn_on            , &!LOGICAL,intent(in) :: stress_ud_satn_on
         stress_vd_satn        , &!REAL,intent(out):: stress_vd_satn(points_stress_vd_satn,0:kMax)! y saturation stress diag
         points_stress_vd_satn , &!INTEGER,intent(in):: points_stress_vd_satn      
         stress_vd_satn_on            , &!LOGICAL,intent(in) :: stress_vd_satn_on
         stress_ud_wake        , &!REAL,intent(out):: stress_ud_wake(points_stress_ud_wake,0:kMax)! x blocked flow stress diag
         points_stress_ud_wake , &!INTEGER,intent(in):: points_stress_ud_wake
         stress_ud_wake_on            , &!LOGICAL,intent(in) :: stress_ud_wake_on
         stress_vd_wake        , &!REAL,intent(out):: stress_vd_wake(points_stress_vd_wake,0:kMax)! y blocked flow stress diag
         points_stress_vd_wake , &!INTEGER,intent(in):: points_stress_vd_wake    
         stress_vd_wake_on            , &!LOGICAL,intent(in) :: stress_vd_wake_on
         du_dt_satn            , &!REAL,intent(out):: du_dt_satn(points_du_dt_satn,kMax)! du/dt diagnostic (saturation)
         points_du_dt_satn     , &!INTEGER,intent(in):: points_du_dt_satn          
         du_dt_satn_on         , &!LOGICAL,intent(in) :: du_dt_satn_on
         du_dt_satn_p_on       , &!LOGICAL,intent(in) :: du_dt_satn_p_on
         dv_dt_satn            , &!REAL,intent(out):: dv_dt_satn(points_dv_dt_satn,kMax)! dv/dt diagnostic (saturation)
         points_dv_dt_satn     , &!INTEGER,intent(in):: points_dv_dt_satn
         dv_dt_satn_on         , &!LOGICAL,intent(in) :: dv_dt_satn_on
         du_dt_wake            , &!REAL,intent(out):: du_dt_wake(points_du_dt_wake,kMax)! du/dt diagnostic (blocked flow)
         points_du_dt_wake     , &!INTEGER,intent(in):: points_du_dt_wake
         du_dt_wake_on         , &!LOGICAL,intent(in) :: du_dt_wake_on
         dv_dt_wake            , &!REAL,intent(out):: dv_dt_wake(points_dv_dt_wake,kMax)! dv/dt diagnostic (blocked flow)
         points_dv_dt_wake     , &!INTEGER,intent(in):: points_dv_dt_wake
         dv_dt_wake_on           )!LOGICAL,intent(in) :: dv_dt_wake_on

    IMPLICIT NONE

    ! Description:
    !     calculates the gwd stress profiles and wind increments.
    !     1. calculate stress profile and wind increments for
    !        linear hydrostatic waves.
    !     2. calculate stress profile and wind increments for
    !        the blocked flow.
    !
    ! current code owner: S. Webster
    !
    ! history:
    ! version  date      comment
    !  5.2   15/11/00   original code.  Stuart Webster.
    !
    !  5.3   16/10/01   Remove lee waves and hydraulic jumps as part
    !                   of the simpler 4A scheme. Stuart Webster
    !
    !  6.2   21/02/06   Introduce minor bug fixes using l_fix_gwsatn.
    !                                             Stuart Webster
    !
    ! code description:
    ! language: fortran 77 + common extensions
    ! this code is written to umdp3 v6 programming standards.
    ! suitable for single column use,rotated grids

    ! global variables

    ! local constants
    ! none

    !
    ! SUBROUTINE ARGUMENTS:
    !

    INTEGER ,INTENT(in) :: kMax                ! number of model levels
    INTEGER ,INTENT(in) :: points                ! number of points
    INTEGER ,INTENT(in) :: k_top(points)         ! model level at mountain tops
    !                                               ! full definition in gwsurf
    INTEGER ,INTENT(in) :: k_top_max             ! max(k_top)

    !
    ! Integers below determine size of diagnostic arrays. They are
    ! set to points if called and to 1 if not. This is done in GWD_CTL2
    !
    INTEGER,INTENT(in):: points_stress_ud
    INTEGER,INTENT(in):: points_stress_vd              
    INTEGER,INTENT(in):: points_stress_ud_satn      
    INTEGER,INTENT(in):: points_stress_vd_satn      
    INTEGER,INTENT(in):: points_stress_ud_wake      
    INTEGER,INTENT(in):: points_stress_vd_wake      
    INTEGER,INTENT(in):: points_du_dt_satn              
    INTEGER,INTENT(in):: points_dv_dt_satn              
    INTEGER,INTENT(in):: points_du_dt_wake              
    INTEGER,INTENT(in):: points_dv_dt_wake

    REAL(KIND=r8),INTENT(in):: r_rho_levels(points,kMax)! heights on rho levels
    REAL(KIND=r8),INTENT(in):: r_theta_levels(points,kMax)! heights on theta levels
    REAL(KIND=r8),INTENT(in):: rho(points,kMax)! density on rho levels
    REAL(KIND=r8),INTENT(in):: theta(points,kMax)  ! theta field
    REAL(KIND=r8),INTENT(in):: u(points,kMax)! u field
    REAL(KIND=r8),INTENT(in):: v(points,kMax)! v field

    REAL(KIND=r8),INTENT(in):: s_x_lin_stress(points) ! 'surface' lin  x_stress
    REAL(KIND=r8),INTENT(in):: s_y_lin_stress(points)! 'surface' lin  y_stress

    REAL(KIND=r8),INTENT(in):: s_x_wake_stress(points)! 'surface' x_wake_stress
    REAL(KIND=r8),INTENT(in):: s_y_wake_stress(points)! 'surface' y_wake_stress

    REAL(KIND=r8),INTENT(in):: s_x_orog(points) ! 'surface' x_stress term
    REAL(KIND=r8),INTENT(in):: s_y_orog(points) ! 'surface' y_stress term
    !REAL(KIND=r8),INTENT(in):: lift    (points) ! blocked layer depth
    REAL(KIND=r8),INTENT(in):: sd_orog(points)  ! standard deviation of the sub-grid orography
    REAL(KIND=r8),INTENT(in):: fr(points)            ! low level froude number
    REAL(KIND=r8),INTENT(in):: kay                   ! GWD hydrostatic constant
    REAL(KIND=r8),INTENT(in):: rho_s(points)         ! Low level density calculated in gwsurf
    REAL(KIND=r8),INTENT(out):: du_dt(points,kMax) ! total GWD du/dt
    REAL(KIND=r8),INTENT(out):: dv_dt(points,kMax)  ! total GWD dv/dt

    !
    ! Diagnostics
    !
    REAL(KIND=r8),INTENT(out):: stress_ud     (points_stress_ud,0:kMax)! x total stress diag
    REAL(KIND=r8),INTENT(out):: stress_vd     (points_stress_vd,0:kMax)! y total stress diag

    REAL(KIND=r8),INTENT(out):: stress_ud_satn(points_stress_ud_satn,0:kMax)! x saturation stress diag
    REAL(KIND=r8),INTENT(out):: stress_vd_satn(points_stress_vd_satn,0:kMax)! y saturation stress diag
    REAL(KIND=r8),INTENT(out):: stress_ud_wake(points_stress_ud_wake,0:kMax)! x blocked flow stress diag
    REAL(KIND=r8),INTENT(out):: stress_vd_wake(points_stress_vd_wake,0:kMax)! y blocked flow stress diag

    REAL(KIND=r8),INTENT(out):: du_dt_satn(points_du_dt_satn,kMax)! du/dt diagnostic (saturation)
    REAL(KIND=r8),INTENT(out):: dv_dt_satn(points_dv_dt_satn,kMax)! dv/dt diagnostic (saturation)

    REAL(KIND=r8),INTENT(out):: du_dt_wake(points_du_dt_wake,kMax)! du/dt diagnostic (blocked flow)
    REAL(KIND=r8),INTENT(out):: dv_dt_wake(points_dv_dt_wake,kMax)! dv/dt diagnostic (blocked flow)

    LOGICAL,INTENT(in):: l_drag(points)        ! true if a non-zero surface stress

    LOGICAL,INTENT(in):: l_fix_gwsatn ! switch to include minor bug fixes  
    LOGICAL,INTENT(in):: l_gwd_40km     ! switch to turn off GWD above 40km

    !
    !  Diagnostic switches
    !
    LOGICAL,INTENT(in) :: stress_ud_on     
    LOGICAL,INTENT(in) :: stress_ud_p_on   
    LOGICAL,INTENT(in) :: stress_vd_on     
    LOGICAL,INTENT(in) :: stress_ud_satn_on
    LOGICAL,INTENT(in) :: stress_vd_satn_on
    LOGICAL,INTENT(in) :: stress_ud_wake_on
    LOGICAL,INTENT(in) :: stress_vd_wake_on
    LOGICAL,INTENT(in) :: du_dt_satn_on    
    LOGICAL,INTENT(in) :: du_dt_satn_p_on  
    LOGICAL,INTENT(in) :: dv_dt_satn_on    
    LOGICAL,INTENT(in) :: du_dt_wake_on    
    LOGICAL,INTENT(in) :: dv_dt_wake_on

    !--------------------------------------------------------------------
    ! LOCAL ARRAYS AND SCALARS
    !--------------------------------------------------------------------

    INTEGER i,k                   ! loop counter in routine

    !      REAL(KIND=r8)                                                              &
    !     & unit_x(points)                                                   &
    !                             ! x_compnt of unit stress vector
    !     &,unit_y(points)        ! y_compnt of unit stress vector

    ! function and subroutine calls:
    !      EXTERNAL gw_satn,gw_wake

    !-------------------------------------------------------------------
    !   1.0 start  preliminaries
    ! initialise increment and increment diagnostics
    !------------------------------------------------------------
    DO k=1,kMax

       DO i=1,points
          du_dt(i,k)=0.0_r8
          dv_dt(i,k)=0.0_r8
       END DO

       IF( du_dt_satn_on .OR. du_dt_satn_p_on ) THEN
          DO i=1,points
             du_dt_satn(i,k)=0.0_r8
          END DO
       END IF

       IF( dv_dt_satn_on ) THEN
          DO i=1,points
             dv_dt_satn(i,k)=0.0_r8
          END DO
       END IF

       IF( du_dt_wake_on ) THEN
          DO i=1,points
             du_dt_wake(i,k)=0.0_r8
          END DO
       END IF

       IF( dv_dt_wake_on ) THEN
          DO i=1,points
             dv_dt_wake(i,k)=0.0_r8
          END DO
       END IF

    END DO ! kMax

    !
    !  and now for stress diagnostics
    !
    DO k=0,kMax

       IF (stress_ud_on .OR. stress_ud_p_on ) THEN
          DO i=1,points
             stress_ud(i,k) = 0.0_r8
          END DO
       END IF

       IF (stress_vd_on ) THEN
          DO i=1,points
             stress_vd(i,k) = 0.0_r8
          END DO
       END IF

       IF (stress_ud_satn_on ) THEN
          DO i=1,points
             stress_ud_satn(i,k) = 0.0_r8
          END DO
       END IF

       IF (stress_vd_satn_on ) THEN
          DO i=1,points
             stress_vd_satn(i,k) = 0.0_r8
          END DO
       END IF

       IF (stress_ud_wake_on ) THEN
          DO i=1,points
             stress_ud_wake(i,k) = 0.0_r8
          END DO
       END IF

       IF (stress_vd_wake_on ) THEN
          DO i=1,points
             stress_vd_wake(i,k) = 0.0_r8
          END DO
       END IF

    END DO

    !---------------------------------------------------------------------
    !  2. launch linear hydrostatic waves from level k_top and
    !     calculate stress profile due to wave saturation effects.
    !---------------------------------------------------------------------
    ! DEPENDS ON: gw_satn
    CALL gw_satn(&
         rho                 (1:points,1:kMax) , & !REAL,intent(in):: rho(points,kMax)! density
         r_rho_levels        (1:points,1:kMax) , & !REAL,intent(in):: r_rho_levels(points,kMax)! heights above z=0 on rho levels
         r_theta_levels      (1:points,1:kMax) , & !REAL,intent(in):: r_theta_levels(points,kMax)! heights above z=0 on theta levels
         theta               (1:points,1:kMax) , & !REAL,intent(in):: theta(points,kMax) ! theta field
         u                   (1:points,1:kMax) , & !REAL,intent(in):: u(points,kMax) ! u field
         v                   (1:points,1:kMax) , & !REAL,intent(in):: v(points,kMax) ! v field
         s_x_lin_stress      (1:points)        , & !REAL,intent(in):: s_x_stress(points)! linear surface x_stress
         s_y_lin_stress      (1:points)        , & !REAL,intent(in):: s_y_stress(points)        ! linear surface y_stress
         kMax                                  , & !INTEGER,intent(in):: kMax                   ! number of model levels
         points                                , & !INTEGER,intent(in):: points                   ! number of points
         kay                                   , & !REAL,intent(in):: kay                  ! GWD Constant (m-1)
         sd_orog             (1:points)        , & !REAL,intent(in):: sd_orog(points)        ! standard deviation of orography
         s_x_orog            (1:points)        , & !REAL,intent(in):: s_x_orog(points)! 'surface' x_orog
         s_y_orog            (1:points)        , & !REAL,intent(in):: s_y_orog(points)! 'surface' y_orog
         du_dt               (1:points,1:kMax) , & !REAL,intent(out):: du_dt(points,kMax)! total GWD du/dt
         dv_dt               (1:points,1:kMax) , & !REAL,intent(out):: dv_dt(points,kMax) ! total GWD dv/dt
         k_top               (1:points)        , & !INTEGER,intent(in):: k_top(points)           ! model level at mountain tops
         k_top_max                             , & !INTEGER,intent(in):: k_top_max            ! max(k_top)
         fr                  (1:points)        , & !REAL,intent(in):: fr(points)                !in low level Froude number
         l_drag              (1:points)        , & !LOGICAL,intent(in):: l_drag(points)           ! whether a point has a
         rho_s               (1:points)        , & !REAL,intent(in):: rho_s(points)     ! surface density (calculated in gwsurf)
         l_fix_gwsatn                          , & !LOGICAL,intent(in):: l_fix_gwsatn           ! Switch at vn6.2 - if true then
         l_gwd_40km                            , & !LOGICAL,intent(in):: l_gwd_40km           ! Switch at vn6.6 - if true then
         ! diagnostics
         stress_ud           (1:points_stress_ud,0:kMax) , & !REAL,intent(out):: stress_ud(points_stress_ud,0:kMax)! u-stress diagnostic
         points_stress_ud                              , & !INTEGER,intent(in):: points_stress_ud
         stress_ud_on                                  , & !LOGICAL,intent(in):: stress_ud_on
         stress_ud_p_on                                , & !LOGICAL,intent(in):: stress_ud_p_on           
         stress_vd           (1:points_stress_vd,0:kMax) , & !REAL,intent(out):: stress_vd(points_stress_vd,0:kMax)! v-stress diagnostic
         points_stress_vd                              , & !INTEGER,intent(in):: points_stress_vd 
         stress_vd_on                                  , & !LOGICAL,intent(in):: stress_vd_on  
         stress_ud_satn       (1:points_stress_ud_satn,0:kMax), & !REAL,intent(out):: stress_ud_satn(points_stress_ud_satn,0:kMax)! x saturation stress diag
         points_stress_ud_satn                         , & !INTEGER,intent(in):: points_stress_ud_satn
         stress_ud_satn_on                             , & !LOGICAL,intent(in):: stress_ud_satn_on    
         stress_vd_satn       (1:points_stress_vd_satn,0:kMax), & !REAL,intent(out):: stress_vd_satn(points_stress_vd_satn,0:kMax)! y saturation stress diag
         points_stress_vd_satn                         , & !INTEGER,intent(in):: points_stress_vd_satn
         stress_vd_satn_on                             , & !LOGICAL,intent(in):: stress_vd_satn_on    
         du_dt_satn           (1:points_du_dt_satn,1:kMax) , & !REAL,intent(out):: du_dt_satn(points_du_dt_satn,kMax)! Saturation du/dt
         points_du_dt_satn                             , & !INTEGER,intent(in):: points_du_dt_satn        
         du_dt_satn_on                                 , & !LOGICAL,intent(in):: du_dt_satn_on   
         du_dt_satn_p_on                               , & !LOGICAL,intent(in):: du_dt_satn_p_on
         dv_dt_satn           (1:points_dv_dt_satn,1:kMax)  , & !REAL,intent(out):: dv_dt_satn(points_dv_dt_satn,kMax)! Saturation dv/dt
         points_dv_dt_satn                             , & !INTEGER,intent(in):: points_dv_dt_satn
         dv_dt_satn_on                                   ) !LOGICAL,intent(in):: dv_dt_satn_on

    !
    !------------------------------------------------------------------
    ! 3 apply uniform stress reduction between surface and k_top
    !     for points where blocked flow was diagnosed.
    !------------------------------------------------------------------

    ! DEPENDS ON: gw_wake
    CALL gw_wake( &
         s_x_wake_stress       (1:points), &!REAL,intent(in):: s_x_stress(points)    ! surface x_stress
         s_y_wake_stress       (1:points), &!REAL,intent(in):: s_y_stress(points)    ! surface y_stress
         kMax                             , &!INTEGER,intent(in):: kMax                   ! number of model levels
         rho                   (1:points,1:kMax), &!REAL,intent(in):: rho(points,kMax)  ! density on rho levels
         !                         r_rho_levels          , &!REAL,intent(in):: r_rho_levels(points,kMax)! heights above z=0 on rho levels
         r_theta_levels        (1:points,1:kMax), &!REAL,intent(in):: r_theta_levels(points,kMax)! heights above z=0 on theta levels
         points                                  , &!INTEGER,intent(in):: points                   ! number of GWD points
         k_top                 (1:points), &!INTEGER,intent(in):: k_top(points)           ! level of blocked layer top
         l_drag                (1:points), &!LOGICAL,intent(in  ) :: l_drag(points)      ! true if a non-zero surface stress
         du_dt                 (1:points,1:kMax), &!REAL,intent(inout):: du_dt(points,kMax)! total GWD du/dt
         dv_dt                 (1:points,1:kMax), &!REAL,intent(inout):: dv_dt(points,kMax)! total GWD dv/dt
         ! diagnostics
         stress_ud             (1:points_stress_ud,0:kMax), &!REAL  ,INTENT(OUT  ) :: stress_ud(points_stress_ud,0:kMax)! u stress  diagnostic 
         points_stress_ud      , &!INTEGER,intent(in   ) :: points_stress_ud
         stress_ud_on          , &!LOGICAL,intent(in   ) :: stress_ud_on
         stress_ud_p_on        , &!LOGICAL,intent(in   ) :: stress_ud_p_on
         stress_vd             (1:points_stress_vd,0:kMax), &!REAL  ,INTENT(OUT  ) :: stress_vd(points_stress_vd,0:kMax)! v stress  diagnostic
         points_stress_vd      , &!INTEGER,intent(in   ) :: points_stress_vd        
         stress_vd_on          , &!LOGICAL,intent(in   ) :: stress_vd_on  
         stress_ud_wake        (1:points_stress_ud_wake,0:kMax), &!REAL  ,INTENT(OUT  ) :: stress_ud_wake(points_stress_ud_wake,0:kMax)! u blocked flow stress diag
         points_stress_ud_wake , &!INTEGER,intent(in   ) :: points_stress_ud_wake
         stress_ud_wake_on     , &!LOGICAL,intent(in   ) :: stress_ud_wake_on
         stress_vd_wake        (1:points_stress_vd_wake,0:kMax), &!REAL  ,INTENT(OUT  ) :: stress_vd_wake(points_stress_vd_wake,0:kMax)! v blocked flow stress diag
         points_stress_vd_wake , &!INTEGER,intent(in   ) :: points_stress_vd_wake
         stress_vd_wake_on     , &!LOGICAL,intent(in   ) :: stress_vd_wake_on
         du_dt_wake            (1:points_du_dt_wake,1:kMax), &!REAL  ,INTENT(OUT  ) :: du_dt_wake(points_du_dt_wake,kMax)! u-acceln  diagnostic
         points_du_dt_wake     , &!INTEGER,intent(in   ) :: points_du_dt_wake
         du_dt_wake_on         , &!LOGICAL,intent(in   ) :: du_dt_wake_on 
         dv_dt_wake            (1:points_dv_dt_wake,1:kMax), &!REAL  ,INTENT(OUT  ) :: dv_dt_wake(points_dv_dt_wake,kMax)! v-acceln  diagnostic
         points_dv_dt_wake     , &!INTEGER,intent(in   ) :: points_dv_dt_wake
         dv_dt_wake_on           )!LOGICAL,intent(in   ) :: dv_dt_wake_on


    RETURN
  END SUBROUTINE gw_vert





  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  ! subroutine gw_wake to deposit blocked flow stress uniformly
  ! from ground up to sub-grid mountain tops
  !
  SUBROUTINE gw_wake(        &
       s_x_stress           , &!REAL,intent(in):: s_x_stress(points)    ! surface x_stress
       s_y_stress           , &!REAL,intent(in):: s_y_stress(points)    ! surface y_stress
       kMax               , &!INTEGER,intent(in):: kMax              ! number of model levels
       rho                  , &!REAL,intent(in):: rho(points,kMax)  ! density on rho levels
       !        r_rho_levels         , &!REAL,intent(in):: r_rho_levels(points,kMax)! heights above z=0 on rho levels
    r_theta_levels       , &!REAL,intent(in):: r_theta_levels(points,kMax)! heights above z=0 on theta levels
         points               , &!INTEGER,intent(in):: points              ! number of GWD points
         k_top                , &!INTEGER,intent(in):: k_top(points)       ! level of blocked layer top
         l_drag               , &!LOGICAL,intent(in  ) :: l_drag(points)      ! true if a non-zero surface stress
         du_dt                , &!REAL,intent(inout):: du_dt(points,kMax)! total GWD du/dt
         dv_dt                , &!REAL,intent(inout):: dv_dt(points,kMax)! total GWD dv/dt
         ! diagnostics
    stress_ud            , &!REAL  ,INTENT(OUT  ) :: stress_ud(points_stress_ud,0:kMax)! u stress  diagnostic 
         points_stress_ud     , &!INTEGER,intent(in   ) :: points_stress_ud
         stress_ud_on         , &!LOGICAL,intent(in   ) :: stress_ud_on
         stress_ud_p_on             , &!LOGICAL,intent(in   ) :: stress_ud_p_on
         stress_vd            , &!REAL  ,INTENT(OUT  ) :: stress_vd(points_stress_vd,0:kMax)! v stress  diagnostic
         points_stress_vd     , &!INTEGER,intent(in   ) :: points_stress_vd     
         stress_vd_on         , &!LOGICAL,intent(in   ) :: stress_vd_on  
         stress_ud_wake       , &!REAL  ,INTENT(OUT  ) :: stress_ud_wake(points_stress_ud_wake,0:kMax)! u blocked flow stress diag
         points_stress_ud_wake, &!INTEGER,intent(in   ) :: points_stress_ud_wake
         stress_ud_wake_on    , &!LOGICAL,intent(in   ) :: stress_ud_wake_on
         stress_vd_wake       , &!REAL  ,INTENT(OUT  ) :: stress_vd_wake(points_stress_vd_wake,0:kMax)! v blocked flow stress diag
         points_stress_vd_wake, &!INTEGER,intent(in   ) :: points_stress_vd_wake
         stress_vd_wake_on    , &!LOGICAL,intent(in   ) :: stress_vd_wake_on
         du_dt_wake           , &!REAL  ,INTENT(OUT  ) :: du_dt_wake(points_du_dt_wake,kMax)! u-acceln  diagnostic
         points_du_dt_wake    , &!INTEGER,intent(in   ) :: points_du_dt_wake
         du_dt_wake_on        , &!LOGICAL,intent(in   ) :: du_dt_wake_on 
         dv_dt_wake           , &!REAL  ,INTENT(OUT  ) :: dv_dt_wake(points_dv_dt_wake,kMax)! v-acceln  diagnostic
         points_dv_dt_wake    , &!INTEGER,intent(in   ) :: points_dv_dt_wake
         dv_dt_wake_on          )!LOGICAL,intent(in   ) :: dv_dt_wake_on


    IMPLICIT NONE

    ! Description: deposit blocked flow stress uniformly between
    !              the surface and the blocked layer top (k_top)
    !
    ! current code owner: S. Webster
    !
    ! history:
    ! version  date      comment
    !  5.2   15/11/00   Original code.      Stuart Webster
    !  5.3   16/10/01   Change k_bot to k_top. Stuart Webster
    !  5.4   28/08/02   Correct header comments. Stuart Webster
    !
    !
    ! code description:
    ! language: fortran 77 + common extensions
    ! this code is written to umdp3 v6 programming standards.
    ! suitable for single column use,rotated grids

    ! local constants
    ! none

    !
    ! SUBROUTINE ARGUMENTS
    !
    INTEGER,INTENT(in):: kMax              ! number of model levels
    INTEGER,INTENT(in):: points              ! number of GWD points
    !
    ! integers below set diagnostic array sizes to points if diagnostic
    ! is called or to 1 if not. These are set in GWD_CTL2
    !
    INTEGER,INTENT(in):: points_stress_ud
    INTEGER,INTENT(in):: points_stress_vd     
    INTEGER,INTENT(in):: points_stress_ud_wake
    INTEGER,INTENT(in):: points_stress_vd_wake
    INTEGER,INTENT(in):: points_du_dt_wake    
    INTEGER,INTENT(in):: points_dv_dt_wake
    INTEGER,INTENT(in):: k_top(points)       ! level of blocked layer top
    REAL(KIND=r8),INTENT(in):: s_x_stress(points)    ! surface x_stress
    REAL(KIND=r8),INTENT(in):: s_y_stress(points)  ! surface y_stress
    !      REAL(KIND=r8),intent(in):: r_rho_levels(points,kMax)! heights above z=0 on rho levels
    REAL(KIND=r8),INTENT(in):: r_theta_levels(points,kMax)! heights above z=0 on theta levels
    REAL(KIND=r8),INTENT(in):: rho(points,kMax)  ! density on rho levels
    REAL(KIND=r8),INTENT(inout):: du_dt(points,kMax)! total GWD du/dt
    REAL(KIND=r8),INTENT(inout):: dv_dt(points,kMax)! total GWD dv/dt

    !
    ! Diagnostic arrays
    !
    REAL(KIND=r8), INTENT(OUT  ) :: du_dt_wake(points_du_dt_wake,kMax)! u-acceln  diagnostic
    REAL(KIND=r8), INTENT(OUT  ) :: dv_dt_wake(points_dv_dt_wake,kMax)! v-acceln  diagnostic
    REAL(KIND=r8), INTENT(INOUT) :: stress_ud(points_stress_ud,0:kMax)! u stress  diagnostic
    REAL(KIND=r8), INTENT(INOUT) :: stress_vd(points_stress_vd,0:kMax)! v stress  diagnostic
    REAL(KIND=r8), INTENT(OUT  ) :: stress_ud_wake(points_stress_ud_wake,0:kMax)! u blocked flow stress diag
    REAL(KIND=r8), INTENT(OUT  ) :: stress_vd_wake(points_stress_vd_wake,0:kMax)! v blocked flow stress diag
    LOGICAL,INTENT(in  ) :: l_drag(points)      ! true if a non-zero surface stress

    !
    ! switches for diagnostics.
    !
    LOGICAL,INTENT(in) :: stress_ud_on  
    LOGICAL,INTENT(in) :: stress_ud_p_on
    LOGICAL,INTENT(in) :: stress_vd_on  
    LOGICAL,INTENT(in) :: stress_ud_wake_on
    LOGICAL,INTENT(in) :: stress_vd_wake_on
    LOGICAL,INTENT(in) :: du_dt_wake_on 
    LOGICAL,INTENT(in) :: dv_dt_wake_on


    !----------------------------------------------------------------
    ! LOCAL ARRAYS AND SCALARS
    !----------------------------------------------------------------
    REAL(KIND=r8) :: x_stress   (points,2)! x_stresses (layer boundaries)
    REAL(KIND=r8) :: y_stress   (points,2)! y_stresses (layer boundaries)
    REAL(KIND=r8) :: dz_x_stress(points)! x component of stress gradient
    REAL(KIND=r8) :: dz_y_stress(points)! y component of stress gradient
    REAL(KIND=r8) :: delta_z              ! difference in height across layer(s)
    INTEGER :: i,k                  ! loop counters
    INTEGER :: kk,kl,ku             ! level counters in routine

     x_stress   =0.0_r8
     y_stress   =0.0_r8
     dz_x_stress=0.0_r8
     dz_y_stress=0.0_r8
     delta_z    =0.0_r8
     delta_z    =0.0_r8
    ! function and subroutine calls
    ! none

    !-------------------------------------------------------------------
    !   1. start level  preliminaries
    !-------------------------------------------------------------------

    kl=1
    ku=2

    DO i=1,points
       IF( l_drag(i) ) THEN

          delta_z = r_theta_levels(i,k_top(i))

          dz_x_stress(i) = s_x_stress(i) / delta_z
          dz_y_stress(i) = s_y_stress(i) / delta_z

       END IF ! if l_drag

    END DO  ! points


    IF( stress_ud_on .OR. stress_ud_p_on .OR. stress_ud_wake_on) THEN

       DO i=1,points
          IF( l_drag(i) ) THEN
             x_stress(i,kl) = s_x_stress(i)
          END IF
       END DO

       IF( stress_ud_on .OR. stress_ud_p_on) THEN
          DO i=1,points
             IF( l_drag(i) ) THEN
                stress_ud(i,0) = stress_ud(i,0) + x_stress(i,kl)
             END IF
          END DO
       END IF

       IF( stress_ud_wake_on ) THEN
          DO i=1,points
             IF( l_drag(i) ) THEN
                stress_ud_wake(i,0) = x_stress(i,kl)
             END IF
          END DO
       END IF

    END IF  !  stress_ud_on .or. stress_ud_p_on .or. stress_ud_wake_on


    IF( stress_vd_on .OR. stress_vd_wake_on) THEN

       DO i=1,points
          IF( l_drag(i) ) THEN
             y_stress(i,kl) = s_y_stress(i)
          END IF
       END DO

       IF( stress_vd_on ) THEN
          DO i=1,points
             IF( l_drag(i) ) THEN
                stress_vd(i,0) = stress_vd(i,0) + y_stress(i,kl)
             END IF
          END DO
       END IF

       IF( stress_vd_wake_on ) THEN
          DO i=1,points
             IF( l_drag(i) ) THEN
                stress_vd_wake(i,0) = y_stress(i,kl)
             END IF
          END DO
       END IF

    END IF  !  stress_vd_on .or. stress_vd_wake_on

    !------------------------------------------------------------------
    !    2 loop levels
    !      k is the rho level counter
    !------------------------------------------------------------------

    DO k=1,kMax


       DO i=1,points
          IF( l_drag(i) .AND. k <= k_top(i) ) THEN
             !
             !  note that a constant stress drop across a level in height coords =>
             !  a non-uniform drag because of rho variation.
             !
             du_dt(i,k) =  du_dt(i,k) - dz_x_stress(i)/rho(i,k)
             dv_dt(i,k) =  dv_dt(i,k) - dz_y_stress(i)/rho(i,k)
             !PRINT*,k, du_dt(i,k),dv_dt(i,k),rho(i,k)
          END IF   ! if l_drag(i) .and. k<=k_top(i)

       END DO ! points

       ! diagnostics
       IF( du_dt_wake_on ) THEN
          DO i=1,points
             IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                du_dt_wake(i,k) = - dz_x_stress(i)/rho(i,k)
             END IF
          END DO
       END IF

       IF( dv_dt_wake_on ) THEN
          DO i=1,points
             IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                dv_dt_wake(i,k) = - dz_y_stress(i)/rho(i,k)
             END IF
          END DO
       END IF

       IF( stress_ud_on .OR. stress_ud_p_on .OR. stress_ud_wake_on) THEN
          DO i=1,points
             IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                IF ( k  ==  1 ) THEN
                   delta_z = r_theta_levels(i,k)
                ELSE
                   delta_z= r_theta_levels(i,k) - r_theta_levels(i,k-1)
                END IF
                x_stress(i,ku) = x_stress(i,kl)
                x_stress(i,ku) = x_stress(i,kl)-dz_x_stress(i)*delta_z
             END IF   ! l_drag & k < k_top
          END DO

          IF( stress_ud_on .OR. stress_ud_p_on) THEN
             DO i=1,points
                IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                   stress_ud(i,k) = stress_ud(i,k) + x_stress(i,ku)
                END IF   ! l_drag & k < k_top
             END DO
          END IF       ! stress_ud on

          IF( stress_ud_wake_on ) THEN
             DO i=1,points
                IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                   stress_ud_wake(i,k) = x_stress(i,ku)
                END IF   ! l_drag & k < k_top
             END DO
          END IF       ! stress_ud_wake on

       END IF      ! stress_ud_on .or. stress_ud_p_on .or. stress_wake_ud_on


       IF( stress_vd_on .OR. stress_vd_wake_on) THEN
          DO i=1,points
             IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                IF ( k  ==  1 ) THEN
                   delta_z = r_theta_levels(i,k)
                ELSE
                   delta_z= r_theta_levels(i,k) - r_theta_levels(i,k-1)
                END IF
                y_stress(i,ku) = y_stress(i,kl)
                y_stress(i,ku) = y_stress(i,kl)-dz_y_stress(i)*delta_z
             END IF   ! l_drag & k < k_top
          END DO

          IF( stress_vd_on ) THEN
             DO i=1,points
                IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                   stress_vd(i,k) = stress_vd(i,k) + y_stress(i,ku)
                END IF   ! l_drag & k < k_top
             END DO
          END IF       ! stress_vd on

          IF( stress_vd_wake_on ) THEN
             DO i=1,points
                IF( l_drag(i) .AND. k <= k_top(i) ) THEN
                   stress_vd_wake(i,k) = y_stress(i,ku)
                END IF   ! l_drag & k < k_top
             END DO
          END IF       ! stress_vd_wake on

       END IF         ! stress_vd_on .or. stress_wake_vd_on


       ! swap storage for lower and upper layers
       kk=kl
       kl=ku
       ku=kk

    END DO
    !   end loop kMax

    RETURN
  END SUBROUTINE gw_wake




  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  !  subroutine gw_surf to calculate surface stress vector for gwd.
  !             also, calculate blocked layer surface stress.
  !
  SUBROUTINE gw_surf(&
       r_theta_levels   , &!REAL ,intent(in):: r_theta_levels(points,kMax)! heights on theta levels
       rho              , &!REAL ,intent(in):: rho(points,kMax)     ! density on rho levels
       theta            , &!REAL ,intent(in):: theta(points,kMax)   ! theta on theta levels
       u                , &!REAL ,intent(in):: u(points,kMax)       ! u on rho levels
       v                , &!REAL ,intent(in):: v(points,kMax)       ! v on rho levels
       timestep         , &!REAL ,intent(in):: timestep                  ! timestep (s)
       sd_orog          , &!REAL,intent(in):: sd_orog(points)        ! standard deviation of the sub-grid orography
       sigma_xx         , &!REAL,intent(in):: sigma_xx(points)       ! (dh/dx)^2 grid box average of the
       sigma_xy         , &!REAL,intent(in):: sigma_xy(points)       ! (dh/dx)*(dh/dy)
       sigma_yy         , &!REAL,intent(in):: sigma_yy(points)       ! (dh/dy)^2
       s_x_lin_stress   , &!REAL,intent(out):: s_x_lin_stress(points)! 'surface' linear stress in x-dirn
       s_y_lin_stress   , &!REAL,intent(out):: s_y_lin_stress(points)! 'surface' linear stress in y-dirn
       s_x_wake_stress  , &!REAL,intent(out):: s_x_wake_stress(points)! wake stress in x-dirn
       s_y_wake_stress  , &!REAL,intent(out):: s_y_wake_stress(points)! wake stress in y-dirn
       s_x_orog         , &!REAL,intent(out):: s_x_orog(points)       ! 'surface' stress/orog in x-dirn
       s_y_orog         , &!REAL    ,intent(out):: s_y_orog(points)       ! 'surface' stress/orog in y-dirn
       kMax           , &!INTEGER ,intent(in) :: kMax                 ! number of model levels
       points           , &!INTEGER ,intent(in) :: points                 ! number of gwd points on pe
       kay              , &!REAL ,intent(in):: kay                    !  gwd surface stress constant (m-1)
       rho_s            , &!REAL,intent(out):: rho_s(points)          ! density - av from z=0 to nsigma*sd_orog
       l_taus_scale     , &!LOGICAL,intent(in):: l_taus_scale           ! true allows variation of surface stress on the low level Froude number
       k_top            , &!INTEGER ,intent(out):: k_top(points)          ! topmost model level including mountains.
       k_top_max        , &!INTEGER ,intent(out):: k_top_max              ! max(k_top)
       lift             , &! REAL,intent(out):: lift(points)           ! blocked layer depth
       l_drag           , &!LOGICAL,intent(out):: l_drag(points)         ! whether a point has a non-zero surface stress or not
       fr               , &!REAL,intent(out):: fr(points)             ! low level froude number for linear
       frc              , &!REAL ,intent(in):: frc                    !  critical froude number below which hydraulic jumps are triggered
       u_s_d            , &!REAL,intent(out):: u_s_d(points_u_s_d)    ! 0-nsigma*sd_orog u_s diag
       u_s_d_on         , &!LOGICAL,intent(in):: u_s_d_on
       points_u_s_d     , &!INTEGER ,intent(in) ::points_u_s_d
       v_s_d            , &!REAL,intent(out):: v_s_d(points_v_s_d)    ! 0-nsigma_sd_orog v_s diag
       v_s_d_on         , &!LOGICAL,intent(in):: v_s_d_on  
       points_v_s_d     , &!INTEGER ,intent(in) ::points_v_s_d
       nsq_s_d          , &!REAL,intent(out):: nsq_s_d(points_nsq_s_d)! 0-nsigma*sd_orog nsq_s diagnostic
       nsq_s_d_on       , &!LOGICAL,intent(in):: nsq_s_d_on          
       points_nsq_s_d   , &!INTEGER ,intent(in) ::points_nsq_s_d        
       fr_d             , &!REAL,intent(out):: fr_d(points_fr_d)      ! Froude no. diagnostic
       fr_d_on          , &!LOGICAL,intent(in):: fr_d_on     
       points_fr_d      , &!INTEGER ,intent(in) ::points_fr_d
       bld_d            , &!REAL,intent(out):: bld_d(points_bld_d)    ! blocked layer depth diagnostic
       bld_d_on         , &!LOGICAL,intent(in):: bld_d_on
       points_bld_d     , &!INTEGER ,intent(in) ::points_bld_d
       bldt_d           , &!REAL,intent(out):: bldt_d(points_bldt_d)  ! %  of time blocked flow param. invoked
       bldt_d_on        , &!LOGICAL,intent(in):: bldt_d_on
       points_bldt_d    , &!INTEGER ,intent(in) ::points_bldt_d
       num_lim_d        , &!REAL,intent(out):: num_lim_d(points_num_lim_d)! % of time numerical limiter invoked
       num_lim_d_on     , &!LOGICAL,intent(in):: num_lim_d_on       
       points_num_lim_d , &!INTEGER ,intent(in) ::points_num_lim_d
       num_fac_d        , &!REAL,intent(out):: num_fac_d(points_num_fac_d)! % reduction of flow blocking stressafter numerical limiter invoked
       num_fac_d_on     , &!LOGICAL,intent(in):: num_fac_d_on           
       points_num_fac_d , &!INTEGER ,intent(in) ::points_num_fac_d
       tausx_d          , &!REAL,intent(out):: tausx_d(points_tausx_d) ! x-component of total surface stress
       tausx_d_on       , &!LOGICAL,intent(in):: tausx_d_on      
       points_tausx_d   , &!INTEGER ,intent(in) ::points_tausx_d
       tausy_d          , &!REAL,intent(out):: tausy_d(points_tausy_d) ! y-component of total surface stress
       tausy_d_on        , &!LOGICAL,intent(in):: tausy_d_on             
       points_tausy_d   , &!INTEGER ,intent(in) ::points_tausy_d
       taus_scale_d     , &!REAL,intent(out):: taus_scale_d(points_taus_scale_d)! scaling factor for surface stress when Fr dependence of surface stress invoked.
       taus_scale_d_on        , &!LOGICAL,intent(in):: taus_scale_d_on
       points_taus_scale_d)!INTEGER ,intent(in) ::points_taus_scale_d

    IMPLICIT NONE

    ! Description:
    !              calculates the total surface stress using the
    !              linear hydrostatic GWD surface stress equation
    !              for a 2-d hill in the absence of rotation and
    !              friction. The blocked layer depth is then used
    !              to partition this stress into a linear hydrostatic
    !              GWD component and a blocked flow component.
    !              nsigma*sd_orog is taken to be the top of the sub-grid
    !              mountains and the low level U and N are calculated
    !              as averages over this layer. This U and N are used
    !              in the calculation of the surface stress and in the
    !              determination of the blocked layer and low level
    !              Froude number.
    !              From vn6.2 it is possible to allow the surface drag to
    !              be dependent on the low level Froude Number. This is
    !              motivated by the results of
    !              Wells et al.2005 (QJRMS, 131, 1321-1338).
    !
    ! Current code owner: S. Webster
    !
    ! history:
    ! version  date      comment
    !  5.2   15/11/00   original code. loosely based on gwsurf3a at vn5.1.
    !                   major overhaul though. see above for more details.
    !                                                Stuart Webster
    !
    !  5.3   16/10/01   Another major overhaul. Total stress now calculated
    !                   using the usual linear hydrostatic GWD surface
    !                   stress equation and this is then split into a
    !                   blocked flow stress and a linear hydrostatic
    !                   GWD stress.
    !                                                Stuart Webster
    !  5.4   28/08/02   Introduce numerical limiter and associated
    !                   diagnostics for flow blocking scheme.
    !                                                Stuart Webster
    !  6.2   15/08/05   Free format fixes. P.Selwood
    !  6.2   21/02/06   Introduce surface stress diagnostics. S. Webster
    !  6.2   21/02/06   Introduce dependence of the total surface drag on
    !                   the low level Froude number. Stuart Webster
    !
    ! code description:
    ! language: fortran 77 + common extensions
    ! this code is written to umdp3 v6 programming standards.
    ! suitable for single column use,rotated grids

    ! global variables
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    !      REAL(KIND=r8), Parameter :: G = 9.80665

    !*----------------------------------------------------------------------
    !*L------------------COMDECK C_PI---------------------------------------
    !LL
    !LL 4.0 19/09/95  New value for PI. Old value incorrect
    !LL               from 12th decimal place. D. Robinson
    !LL 5.1 7/03/00   Fixed/Free format P.Selwood
    !LL

    ! Pi
    REAL(KIND=r8), PARAMETER :: Pi                 = 3.14159265358979323846_r8

    ! Conversion factor degrees to radians
    REAL(KIND=r8), PARAMETER :: Pi_Over_180        = Pi/180.0_r8

    ! Conversion factor radians to degrees
    REAL(KIND=r8), PARAMETER :: Recip_Pi_Over_180  = 180.0_r8/Pi

    !*----------------------------------------------------------------------

    ! local constants
    !
    !  Description: This comdeck defines the constants for the 3B and 4A
    !               versions of the Gravity Wave Drag Code. These are
    !               tuneable parameters but are unlikely to be changed.
    !
    !  History:
    !  Version    Date     Comment
    !  -------    ----     -------
    !    3.4     18/10/94  Original Version    J.R. Mitchell
    !    4.3      7/03/97  Remove KAY_LEE (now set in RUNCNST) S.Webster
    !    4.5     03/08/98  Add GAMMA_SATN (Used in 06_3B). D. Robinson
    !    5.0     14/07/99  Remove redundant switches/variables: only keep
    !                      s-version 06_3B parameters. R Rawlins
    !    5.2     15/11/00  Set parameters for the 4A scheme.
    !                                                  Stuart Webster.
    !    5.3     16/10/01  Partition 3B and 4A scheme parameters.
    !                                                  Stuart Webster.
    !     5.3     21/09/01  Set parameters for the spectral (middle
    !                       atmosphere) gravity wave forcing scheme.
    !                       Warner and McIntyre, J. Atm. Sci., 2001,
    !                       Scaife et al., Geophys. Res. Lett., 2000
    !                       Scaife et al, J. Atm. Sci., 2002 give
    !                       further details.
    !                                                    Adam Scaife
    !    5.4     Alters c_gwave.h include file to change from launch level
    !            to launch eta to make less model level dependent.
    !                                                    Adam Scaife
    !    5.5     25/02/03  Remove 3B GWD parameter settings. Stuart Webster
    !
    !    6.2     16/06/05  Move CCL parameter to gw_ussp. Andrew Bushell
    !
    !
    ! Number of standard deviations above the mean orography of top
    !  of sub-grid mountains
    !REAL(KIND=r8),PARAMETER :: NSIGMA = 2.5_r8


    ! SPECTRAL GRAVITY WAVE SCHEME PARAMETERS:

    ! LMINL = 1/max wavelength at launch
    ! LSTAR = 1/characteristic wavelength
    ! ETALAUNCH = eta for model launch

    REAL(KIND=r8),PARAMETER:: LMINL = 1.0_r8/20000.0_r8
    REAL(KIND=r8),PARAMETER:: LSTAR = 1.0_r8/4300.0_r8
    REAL(KIND=r8),PARAMETER:: ETALAUNCH = 0.045_r8

    !
    ! SUBROUTINE ARGUMENTS
    !
    INTEGER ,INTENT(in) :: kMax                 ! number of model levels
    INTEGER ,INTENT(in) :: points                 ! number of gwd points on pe
    INTEGER ,INTENT(out):: k_top(points)          ! topmost model level including mountains.
    !                             ! defined to be the kth rho level in
    !                             ! which the nsigma*sd_orog height lies,
    !                             ! i.e.nsigma*sd_orog lies between
    !                             ! r_theta_levels(k-1) and
    !                             ! r_theta_levels(k).

    INTEGER ,INTENT(out):: k_top_max              ! max(k_top)

    !
    ! Dimensions of diagnostic arrays.
    ! dimension = points if diag called,
    !           = 1      if diag not called.
    ! These are set in GWD_CTL2
    !
    INTEGER ,INTENT(in) ::points_nsq_s_d        
    INTEGER ,INTENT(in) ::points_u_s_d        
    INTEGER ,INTENT(in) ::points_v_s_d        
    INTEGER ,INTENT(in) ::points_fr_d        
    INTEGER ,INTENT(in) ::points_bld_d        
    INTEGER ,INTENT(in) ::points_bldt_d        
    INTEGER ,INTENT(in) ::points_num_lim_d        
    INTEGER ,INTENT(in) ::points_num_fac_d        
    INTEGER ,INTENT(in) ::points_tausx_d        
    INTEGER ,INTENT(in) ::points_tausy_d        
    INTEGER ,INTENT(in) ::points_taus_scale_d

    REAL(KIND=r8) ,INTENT(in):: r_theta_levels(points,kMax)! heights on theta levels
    REAL(KIND=r8) ,INTENT(in):: rho(points,kMax)     ! density on rho levels
    REAL(KIND=r8) ,INTENT(in):: theta(points,kMax)   ! theta on theta levels
    REAL(KIND=r8) ,INTENT(in):: u(points,kMax)       ! u on rho levels
    REAL(KIND=r8) ,INTENT(in):: v(points,kMax)       ! v on rho levels
    REAL(KIND=r8) ,INTENT(in):: timestep               ! timestep (s)
    REAL(KIND=r8) ,INTENT(in):: kay                    !  gwd surface stress constant (m-1)
    REAL(KIND=r8) ,INTENT(in):: frc                    !  critical froude number below which
    !                                             !  hydraulic jumps are triggered

    REAL(KIND=r8),INTENT(in):: sd_orog(points)        ! standard deviation of the sub-grid
    !                                            ! orography
    REAL(KIND=r8),INTENT(in):: sigma_xx(points)       ! (dh/dx)^2 grid box average of the
    !                             ! the hi-res source dataset. here,
    !                             ! dh is the height change of the hi-res
    !                             ! data and dx is the distance between
    !                             ! adjacent data points.

    REAL(KIND=r8),INTENT(in):: sigma_xy(points)       ! (dh/dx)*(dh/dy)
    REAL(KIND=r8),INTENT(in):: sigma_yy(points)       ! (dh/dy)^2
    REAL(KIND=r8),INTENT(out):: s_x_lin_stress  (points) ! 'surface' linear stress in x-dirn
    REAL(KIND=r8),INTENT(out):: s_y_lin_stress  (points) ! 'surface' linear stress in y-dirn
    REAL(KIND=r8),INTENT(out):: s_x_wake_stress (points)! wake stress in x-dirn
    REAL(KIND=r8),INTENT(out):: s_y_wake_stress (points)! wake stress in y-dirn

    REAL(KIND=r8),INTENT(out):: s_x_orog        (points)       ! 'surface' stress/orog in x-dirn
    REAL(KIND=r8),INTENT(out):: s_y_orog        (points)       ! 'surface' stress/orog in y-dirn

    REAL(KIND=r8),INTENT(out):: rho_s           (points)          ! density - av from z=0 to nsigma*sd_orog
    REAL(KIND=r8),INTENT(out):: lift            (points)           ! blocked layer depth
    REAL(KIND=r8),INTENT(out):: fr              (points)             ! low level froude number for linear
    !                                             ! hydrostatic waves - so nsigma*sd_orog
    !                                             ! is replaced by nsigma*sd_orog-lift
    REAL(KIND=r8),INTENT(out):: u_s_d           (points_u_s_d)    ! 0-nsigma*sd_orog u_s diag
    REAL(KIND=r8),INTENT(out):: v_s_d           (points_v_s_d)    ! 0-nsigma_sd_orog v_s diag
    REAL(KIND=r8),INTENT(out):: nsq_s_d         (points_nsq_s_d)! 0-nsigma*sd_orog nsq_s diagnostic
    REAL(KIND=r8),INTENT(out):: fr_d            ( points_fr_d)      ! Froude no. diagnostic
    REAL(KIND=r8),INTENT(out):: bld_d           (points_bld_d)    ! blocked layer depth diagnostic
    REAL(KIND=r8),INTENT(out):: bldt_d          (points_bldt_d)  ! %  of time blocked flow param. invoked

    REAL(KIND=r8),INTENT(out):: num_lim_d       (points_num_lim_d)! % of time numerical limiter invoked
    REAL(KIND=r8),INTENT(out):: num_fac_d       (points_num_fac_d)! % reduction of flow blocking stress
    ! after numerical limiter invoked
    REAL(KIND=r8),INTENT(out):: tausx_d         (points_tausx_d) ! x-component of total surface stress
    REAL(KIND=r8),INTENT(out):: tausy_d         (points_tausy_d) ! y-component of total surface stress

    REAL(KIND=r8),INTENT(out):: taus_scale_d    (points_taus_scale_d)! scaling factor for surface stress when
    !                                                   ! Fr dependence of surface stress invoked.

    LOGICAL,INTENT(in):: l_taus_scale           ! true allows variation of surface stress
    !                             ! on the low level Froude number

    LOGICAL,INTENT(out):: l_drag(points)         ! whether a point has a non-zero surface
    !                             ! stress or not

    !
    ! Diagnostic switches
    !
    LOGICAL,INTENT(in):: u_s_d_on
    LOGICAL,INTENT(in):: v_s_d_on                          
    LOGICAL,INTENT(in):: nsq_s_d_on                        
    LOGICAL,INTENT(in):: fr_d_on                           
    LOGICAL,INTENT(in):: bld_d_on                          
    LOGICAL,INTENT(in):: bldt_d_on                         
    LOGICAL,INTENT(in):: num_lim_d_on                      
    LOGICAL,INTENT(in):: num_fac_d_on                      
    LOGICAL,INTENT(in):: tausx_d_on                        
    LOGICAL,INTENT(in):: tausy_d_on                        
    LOGICAL,INTENT(in):: taus_scale_d_on

    !------------------------------------------------------------------
    ! LOCAL ARRAYS AND SCALARS
    !------------------------------------------------------------------

    REAL(KIND=r8):: u_s(points)     ! u-winds - av from z=0 to nsigma*sd_orog
    REAL(KIND=r8):: v_s(points)     ! v-winds - av from z=0 to nsigma*sd_orog

    REAL(KIND=r8):: nsq_s(points)   ! n squared av from z=0 to nsigma*sd_orog
    REAL(KIND=r8):: fr_for_diag(points)    ! low level froude number

    REAL(KIND=r8):: num_lim(points)        ! 0 if limiter not invoked, 100 if it is.

    REAL(KIND=r8):: num_fac(points)   ! percentage reduction of the
    !                             ! flow-blocking stress after the limiter
    !                             ! was invoked

    REAL(KIND=r8):: taus_scale(points)     ! factor by which surf stress is scaled
    !                             ! when Froude no. dependency is invoked.

    INTEGER :: i,k                    ! loop counter in routine

    REAL(KIND=r8) :: speed        !     wind speed / wind speed in dirn stress
    REAL(KIND=r8) :: speedcalc    ! numerator of calcuation for speed
    REAL(KIND=r8) :: s_stress_sq  ! denominater of calculation for speed
    REAL(KIND=r8) :: s_wake_stress! amplitude of surface wake stress
    REAL(KIND=r8) :: s_wake_limit ! numerical limit for s_wake_stress
    REAL(KIND=r8) :: n                      ! brunt_vaisala frequency
    REAL(KIND=r8) :: r_frc                  ! reciprocal of the critical froude number
    REAL(KIND=r8) :: calc         ! calculation for surface stress magnitude
    REAL(KIND=r8) :: calc1                  ! as per calc

    REAL(KIND=r8) :: n_squared              ! n^2 on rho levels

    REAL(KIND=r8) :: dzb             ! relevant depth to vertical average
    !                             ! of current theta level
    REAL(KIND=r8) :: dzt                    ! height from surface to kth theta level
    !                             ! or nsigma*sd_orog

    LOGICAL l_cont(points)         ! level continue


    ! function and subroutine calls: none
    speed      = 0.0_r8 
    speedcalc   = 0.0_r8  
    s_stress_sq  = 0.0_r8 
    s_wake_stress = 0.0_r8
    s_wake_limit  = 0.0_r8
    n       = 0.0_r8
    r_frc       = 0.0_r8
    calc       = 0.0_r8
    calc1       = 0.0_r8
    n_squared     = 0.0_r8
    dzb             = 0.0_r8
    dzt       = 0.0_r8
    s_x_lin_stress   = 0.0_r8
    s_y_lin_stress   = 0.0_r8
    s_x_wake_stress  = 0.0_r8
    s_y_wake_stress  = 0.0_r8

    s_x_orog            = 0.0_r8
    s_y_orog            = 0.0_r8

    rho_s            = 0.0_r8
    lift            = 0.0_r8
    fr              = 0.0_r8
    u_s_d            = 0.0_r8
    v_s_d            = 0.0_r8
    nsq_s_d            = 0.0_r8
    fr_d            = 0.0_r8
    bld_d            = 0.0_r8
    bldt_d            = 0.0_r8
    num_lim_d            = 0.0_r8
    num_fac_d        = 0.0_r8
    tausx_d         = 0.0_r8
    tausy_d         = 0.0_r8
    taus_scale_d  = 0.0_r8
    !---------------------------------------------------------------------
    ! 1.0 initialisation
    !---------------------------------------------------------------------

    k_top_max = 1
    r_frc     = 1.0_r8/frc

    DO i=1,points
       nsq_s(i)       = 0.0_r8
    END DO

    DO i=1,points
       rho_s(i)       = 0.0_r8
    END DO

    DO i=1,points
       u_s(i)         = 0.0_r8
    END DO

    DO i=1,points
       v_s(i)         = 0.0_r8
    END DO

    DO i=1,points
       l_cont(i)     =.TRUE.
    END DO

    DO i=1,points
       s_x_lin_stress(i)  = 0.0_r8
    END DO

    DO i=1,points
       s_y_lin_stress(i)  = 0.0_r8
    END DO

    DO i=1,points
       s_x_wake_stress(i) = 0.0_r8
    END DO

    DO i=1,points
       s_y_wake_stress(i) = 0.0_r8
    END DO

    DO i=1,points
       l_drag(i) = .TRUE.
    END DO

    DO i=1,points
       num_lim(i) = 0.0_r8
    END DO

    DO i=1,points
       num_fac(i) = 0.0_r8
    END DO

    DO i=1,points
       taus_scale(i) = 1.0_r8
    END DO

    !------------------------------------------------------------------
    ! 2.
    !    Calculation of the average surface quantities and
    !    the nsigma*sd_orog point values. The surface
    !    is the average from 0 up to nsigma*sd_orog.
    !    All interpolations now done in terms of height.
    !    k_top is defined to be the rho level containing nsigma*sd_orog
    !    dzb=dzt when k=2 because include level 1 depth in z, but
    !    don't calculate u and n at rho level 1.
    !
    !---------------------------------------------------------------------
    DO i=1,points
       IF (sd_orog(i)  <=  0.0_r8 ) THEN
          l_drag(i) = .FALSE.
          l_cont(i) = .FALSE.
          k_top(i)  = 2
       END IF
    END DO

    DO k=2,kMax-1
       DO i=1,points
          IF ( l_cont(i) ) THEN

             IF ( r_theta_levels(i,k)  <   nsigma*sd_orog(i) ) THEN
                dzt = r_theta_levels(i,k)
                IF ( k  ==  2 ) THEN
                   dzb = dzt
                ELSE
                   dzb = r_theta_levels(i,k) -  r_theta_levels(i,k-1)
                END IF
             ELSE

                dzt = nsigma * sd_orog(i)
                IF ( k  ==  2 ) THEN
                   dzb = dzt
                ELSE
                   dzb = nsigma*sd_orog(i) -  r_theta_levels(i,k-1)
                END IF

                k_top(i)   =  k
                IF( k_top_max  <   k ) k_top_max=k
                l_cont(i) = .FALSE.

             END IF   ! r_theta_levels(i,k) >  nsigma*sd_orog(i)

             n_squared = 2.0_r8*g*(theta(i,k)-theta(i,k-1) )                 &
                  &                 /( (theta(i,k)+theta(i,k-1)) *                   &
                  &                    (r_theta_levels(i,k)-r_theta_levels(i,k-1)) )
             !
             ! Form z=0 to current level averages which
             ! when k=k_top become 0-nsigma*sd_orog averages
             !
             u_s(i)   = ( u_s(i)   * r_theta_levels(i,k-1) + u(i,k) * dzb    ) / dzt
             
             v_s(i)   = ( v_s(i)   * r_theta_levels(i,k-1) + v(i,k) * dzb    ) / dzt

             nsq_s(i) = ( nsq_s(i) * r_theta_levels(i,k-1) + n_squared* dzb  ) / dzt
             
             rho_s(i) = ( rho_s(i) * r_theta_levels(i,k-1) + rho(i,k) * dzb  ) / dzt

          END IF  !            l_cont

       END DO

    END DO

    !------------------------------------------------------------------
    ! 3.1 calculation of total surface stress and blocked
    !     layer depth, and the split of the total stress into a linear
    !     hydrostatic gravity wave stress and a blocked flow stress
    !-----------------------------------------------------------------
    DO i=1,points

       !
       ! First calculate the speed of the wind for the linear
       ! hydrostatic code
       !
       speed = u_s(i)*u_s(i) + v_s(i)*v_s(i)

       IF ( speed  <=  0.0_r8 ) THEN
          s_x_orog(i) = 0.0_r8
          s_y_orog(i) = 0.0_r8
       ELSE
          speed = SQRT(speed)
          s_x_orog(i)= (u_s(i)*sigma_xx(i) + v_s(i)*sigma_xy(i)) /speed
          s_y_orog(i)= (u_s(i)*sigma_xy(i) + v_s(i)*sigma_yy(i)) /speed
          s_stress_sq= s_x_orog(i)*s_x_orog(i) + s_y_orog(i)*s_y_orog(i)
          speedcalc  = u_s(i)*s_x_orog(i) + v_s(i)*s_y_orog(i)
          IF ( s_stress_sq  <=  0.0_r8 ) THEN
             speed    = 0.0_r8
          ELSE
             ! speed is the component of the wind perpendicular to the major
             ! axis of the orography.
             speed    = speedcalc / SQRT( s_stress_sq )
          END IF
       END IF

       IF ( nsq_s(i) >  0.0_r8 .AND. speed >  0.0_r8 .AND. l_drag(i) ) THEN

          n              = SQRT( nsq_s(i) )
          fr_for_diag(i) = speed / (n*nsigma*sd_orog(i))
          lift(i)        = nsigma * sd_orog(i) - r_frc * speed/n
          IF ( lift(i)  <   0.0_r8 ) THEN
             lift(i)      = 0.0_r8
          END IF

          !
          !  fr here is for the linear hydrostatic surface stress and
          !  then also for the critical stress calculation in gw_satn
          !
          fr(i) = speed / ( n*(nsigma*sd_orog(i)-lift(i)))
          !
          ! Calculate factor by which to scale surface stresses when dependence
          ! of surface stress on Froude number is invoked. Empirical expression
          ! is a simple fit to the experimental results plotted in Fig. 7 of
          ! Wells et al. (2005).
          !
          ! This function was modified at vn6.6 to keep ratio of drag:linear drag 
          ! equal to 1 for Froude numbers>2 
          !
          IF ( l_taus_scale ) THEN

             IF ( fr_for_diag(i) <= 1.0_r8 ) THEN
                taus_scale(i) = 1.6_r8*fr_for_diag(i)
             ELSE IF ( fr_for_diag(i) <= 2.0_r8 ) THEN
                taus_scale(i) = 2.2_r8 - 0.6_r8*fr_for_diag(i)
             END IF

          END IF

          calc  = kay * rho_s(i) * speed**3 * taus_scale(i) /           &
               &            (n*nsigma*nsigma*sd_orog(i)*sd_orog(i)*fr(i)*fr(i))
          s_x_lin_stress(i) = s_x_orog(i) * calc
          s_y_lin_stress(i) = s_y_orog(i) * calc

          !
          ! s_x_orog(i)*calc1 = linear hydrostatic prediction of the stress
          !
          calc1              = kay*rho_s(i)* speed * n * taus_scale(i)
          s_x_wake_stress(i) = s_x_orog(i) * calc1 - s_x_lin_stress(i)
          s_y_wake_stress(i) = s_y_orog(i) * calc1 - s_y_lin_stress(i)
          !
          ! Limit wake_stress so that numerical instability is not permitted
          !
          s_wake_stress = SQRT( s_x_wake_stress(i)*s_x_wake_stress(i)   &
               &               +s_y_wake_stress(i)*s_y_wake_stress(i))
          s_wake_limit  = 2*nsigma*sd_orog(i)*speed*rho_s(i)/timestep

          IF (s_wake_stress  >   s_wake_limit) THEN
             s_x_wake_stress(i) = s_x_wake_stress(i) * s_wake_limit      &
                  &                                              /s_wake_stress
             s_y_wake_stress(i) = s_y_wake_stress(i) * s_wake_limit      &
                  &                                              /s_wake_stress
             num_lim(i) = 100.0_r8
             num_fac(i) = 100.0_r8 * (s_wake_stress - s_wake_limit )         &
                  &                          / s_wake_stress
          END IF

       ELSE

          l_drag(i)      = .FALSE.
          lift(i)        = 0.0_r8
          fr_for_diag(i) = 0.0_r8
          fr(i)          = 0.0_r8

       END IF     ! speed or n or orog  >   0.0_r8

    END DO      ! i=points


    !-----------------------------------------------------------------
    ! 4 diagnostics
    !-----------------------------------------------------------------
    IF ( u_s_d_on ) THEN
       DO i=1,points
          u_s_d(i) = u_s(i)
       END DO
    END IF

    IF ( v_s_d_on ) THEN
       DO i=1,points
          v_s_d(i) = v_s(i)
       END DO
    END IF

    IF ( nsq_s_d_on ) THEN
       DO i=1,points
          IF ( nsq_s(i)  <   0.0_r8 ) nsq_s(i) = 0.0_r8
          nsq_s_d(i) = SQRT( nsq_s(i) )
       END DO
    END IF

    IF ( bld_d_on ) THEN
       DO i=1,points
          bld_d(i) = lift(i)
       END DO
    END IF

    IF ( bldt_d_on ) THEN
       DO i=1,points
          IF( s_x_wake_stress(i)  /=  0.0_r8 .OR.                          &
               &        s_y_wake_stress(i) /= 0.0_r8 )  THEN
             bldt_d(i) = 100.0_r8
          ELSE
             bldt_d(i) = 0.0_r8
          END IF
       END DO
    END IF

    IF ( fr_d_on ) THEN
       DO i=1,points
          fr_d(i) = fr_for_diag(i)
          !  limit Fr to sensible-ish values as previously had problems
          !  when creating time means which include very large values.
          IF (fr_d(i)  >   1000.0_r8) fr_d(i) = 1000.0_r8
       END DO
    END IF

    IF ( num_lim_d_on ) THEN
       DO i=1,points
          num_lim_d(i) = num_lim(i)
       END DO
    END IF

    IF ( num_fac_d_on ) THEN
       DO i=1,points
          num_fac_d(i) = num_fac(i)
       END DO
    END IF

    IF ( tausx_d_on ) THEN
       DO i=1,points
          tausx_d(i) = s_x_lin_stress(i) + s_x_wake_stress(i)
       END DO
    END IF

    IF ( tausy_d_on ) THEN
       DO i=1,points
          tausy_d(i) = s_y_lin_stress(i) + s_y_wake_stress(i)
       END DO
    END IF

    IF ( taus_scale_d_on ) THEN
       DO i=1,points
          taus_scale_d(i) = taus_scale(i)
       END DO
    END IF


    RETURN
  END SUBROUTINE gw_surf








  ! *****************************COPYRIGHT*******************************
  ! (C) Crown copyright Met Office. All rights reserved.
  ! For further details please refer to the file COPYRIGHT.txt
  ! which you should have received as part of this distribution.
  ! *****************************COPYRIGHT*******************************
  !
  ! SUBROUTINE GW_SATN: Calculates stress profile for linear hydrostatic
  !                     Waves
  !
  SUBROUTINE gw_satn(&
       rho                   , & !REAL,intent(in):: rho(points,kMax)! density
       r_rho_levels          , & !REAL,intent(in):: r_rho_levels(points,kMax)! heights above z=0 on rho levels
       r_theta_levels        , & !REAL,intent(in):: r_theta_levels(points,kMax)! heights above z=0 on theta levels
       theta                 , & !REAL,intent(in):: theta(points,kMax) ! theta field
       u                     , & !REAL,intent(in):: u(points,kMax) ! u field
       v                     , & !REAL,intent(in):: v(points,kMax) ! v field
       s_x_stress            , & !REAL,intent(in):: s_x_stress(points)! linear surface x_stress
       s_y_stress            , & !REAL,intent(in):: s_y_stress(points)   ! linear surface y_stress
       kMax                , & !INTEGER,intent(in):: kMax               ! number of model levels
       points                , & !INTEGER,intent(in):: points               ! number of points
       kay                   , & !REAL,intent(in):: kay                  ! GWD Constant (m-1)
       sd_orog               , & !REAL,intent(in):: sd_orog(points)      ! standard deviation of orography
       s_x_orog              , & !REAL,intent(in):: s_x_orog(points)! 'surface' x_orog
       s_y_orog              , & !REAL,intent(in):: s_y_orog(points)! 'surface' y_orog
       du_dt                 , & !REAL,intent(out):: du_dt(points,kMax)! total GWD du/dt
       dv_dt                 , & !REAL,intent(out):: dv_dt(points,kMax) ! total GWD dv/dt
       k_top                 , & !INTEGER,intent(in):: k_top(points)        ! model level at mountain tops
       k_top_max             , & !INTEGER,intent(in):: k_top_max            ! max(k_top)
       fr                    , & !REAL,intent(in):: fr(points)           !in low level Froude number
       l_drag                , & !LOGICAL,intent(in):: l_drag(points)       ! whether a point has a
       rho_s                 , & !REAL,intent(in):: rho_s(points)     ! surface density (calculated in gwsurf)
       l_fix_gwsatn          , & !LOGICAL,intent(in):: l_fix_gwsatn         ! Switch at vn6.2 - if true then
       l_gwd_40km            , & !LOGICAL,intent(in):: l_gwd_40km           ! Switch at vn6.6 - if true then
       ! diagnostics
    stress_ud             , & !REAL,intent(out):: stress_ud(points_stress_ud,0:kMax)! u-stress diagnostic
         points_stress_ud      , & !INTEGER,intent(in):: points_stress_ud
         stress_ud_on          , & !LOGICAL,intent(in):: stress_ud_on
         stress_ud_p_on        , & !LOGICAL,intent(in):: stress_ud_p_on       
         stress_vd             , & !REAL,intent(out):: stress_vd(points_stress_vd,0:kMax)! v-stress diagnostic
         points_stress_vd      , & !INTEGER,intent(in):: points_stress_vd        
         stress_vd_on          , & !LOGICAL,intent(in):: stress_vd_on  
         stress_ud_satn        , & !REAL,intent(out):: stress_ud_satn(points_stress_ud_satn,0:kMax)! x saturation stress diag
         points_stress_ud_satn , & !INTEGER,intent(in):: points_stress_ud_satn
         stress_ud_satn_on     , & !LOGICAL,intent(in):: stress_ud_satn_on    
         stress_vd_satn        , & !REAL,intent(out):: stress_vd_satn(points_stress_vd_satn,0:kMax)! y saturation stress diag
         points_stress_vd_satn , & !INTEGER,intent(in):: points_stress_vd_satn
         stress_vd_satn_on     , & !LOGICAL,intent(in):: stress_vd_satn_on    
         du_dt_satn            , & !REAL,intent(out):: du_dt_satn(points_du_dt_satn,kMax)! Saturation du/dt
         points_du_dt_satn     , & !INTEGER,intent(in):: points_du_dt_satn        
         du_dt_satn_on         , & !LOGICAL,intent(in):: du_dt_satn_on   
         du_dt_satn_p_on       , & !LOGICAL,intent(in):: du_dt_satn_p_on
         dv_dt_satn            , & !REAL,intent(out):: dv_dt_satn(points_dv_dt_satn,kMax)! Saturation dv/dt
         points_dv_dt_satn     , & !INTEGER,intent(in):: points_dv_dt_satn
         dv_dt_satn_on           ) !LOGICAL,intent(in):: dv_dt_satn_on

    IMPLICIT NONE

    ! Description:
    !
    !   Calculate stress profiles and hence drags due to linear
    !   hydrostatic gravity waves. in contrast to the 3a/3b schemes,
    !   linear hydrostatic waves may occur for any Froude number now.
    !   When the Froude number drops below Frc (as specified in c_gwave.h),
    !   the wave amplitude is reduced to be that of the air that is not
    !   blocked, rather than be the full depth of the sub-grid mountains.
    !   Thus, as the Froude number reduces to zero, so the linear
    !   hydrostatic wave amplitude reduces to zero.
    !
    !
    ! current code owner: S. Webster
    !
    ! history:
    ! version  date      comment
    !  5.2   15/11/00   Original code. Based on 5.1 GWSATN3B
    !                   code now includes the option to use the
    !                   turning of the wind algorithm. However,
    !                   by default, l_gwd_tow is set to false
    !                   in c_gwave.h. The saturation part of the code
    !                   is also now on a switch in the same header.
    !                   L_gwd_sat is true by default.
    !                                               Stuart Webster.
    ! 5.3    16/10/01   Change to nsigma*sd_orog for the mountain top
    !                   heights.                    Stuart Webster.
    ! 6.2    15/08/05   Free format fixes. P.Selwood
    ! 6.2    21/02/06   Two minor bug fixes - invoked using l_fix_gwsatn
    !                                               Stuart Webster
    !
    !
    ! Language: fortran 77 + common extensions
    ! This code is written to umdp3 v6 programming standards.
    ! Suitable for single column use,rotated grids

    ! Global variables
    !*L------------------COMDECK C_G----------------------------------------
    ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

    REAL(KIND=r8), PARAMETER :: G = 9.80665_r8

    !*----------------------------------------------------------------------

    ! Local constants
    !
    !  Description: This comdeck defines the constants for the 3B and 4A
    !               versions of the Gravity Wave Drag Code. These are
    !               tuneable parameters but are unlikely to be changed.
    !
    !  History:
    !  Version    Date     Comment
    !  -------    ----     -------
    !    3.4     18/10/94  Original Version    J.R. Mitchell
    !    4.3      7/03/97  Remove KAY_LEE (now set in RUNCNST) S.Webster
    !    4.5     03/08/98  Add GAMMA_SATN (Used in 06_3B). D. Robinson
    !    5.0     14/07/99  Remove redundant switches/variables: only keep
    !                      s-version 06_3B parameters. R Rawlins
    !    5.2     15/11/00  Set parameters for the 4A scheme.
    !                                                  Stuart Webster.
    !    5.3     16/10/01  Partition 3B and 4A scheme parameters.
    !                                                  Stuart Webster.
    !     5.3     21/09/01  Set parameters for the spectral (middle
    !                       atmosphere) gravity wave forcing scheme.
    !                       Warner and McIntyre, J. Atm. Sci., 2001,
    !                       Scaife et al., Geophys. Res. Lett., 2000
    !                       Scaife et al, J. Atm. Sci., 2002 give
    !                       further details.
    !                                                    Adam Scaife
    !    5.4     Alters c_gwave.h include file to change from launch level
    !            to launch eta to make less model level dependent.
    !                                                    Adam Scaife
    !    5.5     25/02/03  Remove 3B GWD parameter settings. Stuart Webster
    !
    !    6.2     16/06/05  Move CCL parameter to gw_ussp. Andrew Bushell
    !
    !
    ! Number of standard deviations above the mean orography of top
    !  of sub-grid mountains
    !REAL(KIND=r8),PARAMETER :: NSIGMA = 2.5_r8


    ! SPECTRAL GRAVITY WAVE SCHEME PARAMETERS:

    ! LMINL = 1/max wavelength at launch
    ! LSTAR = 1/characteristic wavelength
    ! ETALAUNCH = eta for model launch

    REAL(KIND=r8),PARAMETER:: LMINL = 1.0_r8/20000.0_r8
    REAL(KIND=r8),PARAMETER:: LSTAR = 1.0_r8/4300.0_r8
    REAL(KIND=r8),PARAMETER:: ETALAUNCH = 0.045_r8

    ! subroutine arguements;

    INTEGER,INTENT(in):: kMax               ! number of model levels
    INTEGER,INTENT(in):: points               ! number of points

    INTEGER,INTENT(in):: k_top(points)        ! model level at mountain tops
    INTEGER,INTENT(in):: k_top_max            ! max(k_top)

    !
    ! The integers below are set in GWD_CTL2 to points if the diagnostics
    ! are called and to 1 if not.
    !
    INTEGER,INTENT(in):: points_stress_ud
    INTEGER,INTENT(in):: points_stress_vd        
    INTEGER,INTENT(in):: points_stress_ud_satn
    INTEGER,INTENT(in):: points_stress_vd_satn
    INTEGER,INTENT(in):: points_du_dt_satn        
    INTEGER,INTENT(in):: points_dv_dt_satn

    REAL(KIND=r8),INTENT(in):: r_rho_levels(points,kMax)! heights above z=0 on rho levels

    REAL(KIND=r8),INTENT(in):: r_theta_levels(points,kMax)! heights above z=0 on theta levels

    REAL(KIND=r8),INTENT(in):: rho(points,kMax)! density
    REAL(KIND=r8),INTENT(in):: rho_s(points)     ! surface density (calculated in gwsurf)
    REAL(KIND=r8),INTENT(in):: theta(points,kMax) ! theta field
    REAL(KIND=r8),INTENT(in):: u(points,kMax) ! u field
    REAL(KIND=r8),INTENT(in):: v(points,kMax) ! v field

    REAL(KIND=r8),INTENT(in):: s_x_stress(points)! linear surface x_stress
    REAL(KIND=r8),INTENT(in):: s_y_stress(points)   ! linear surface y_stress

    REAL(KIND=r8),INTENT(in):: s_x_orog(points)! 'surface' x_orog
    REAL(KIND=r8),INTENT(in):: s_y_orog(points)! 'surface' y_orog

    REAL(KIND=r8),INTENT(in):: sd_orog(points)      ! standard deviation of orography
    REAL(KIND=r8),INTENT(in):: kay                  ! GWD Constant (m-1)

    REAL(KIND=r8),INTENT(in):: fr(points)           !in low level Froude number

    REAL(KIND=r8),INTENT(out):: du_dt          (points,kMax)! total GWD du/dt
    REAL(KIND=r8),INTENT(out):: dv_dt          (points,kMax) ! total GWD dv/dt

    ! Diagnostics

    REAL(KIND=r8),INTENT(out):: du_dt_satn     (points_du_dt_satn,kMax)        ! Saturation du/dt
    REAL(KIND=r8),INTENT(out):: dv_dt_satn     (points_dv_dt_satn,kMax)      ! Saturation dv/dt
    REAL(KIND=r8),INTENT(out):: stress_ud      (points_stress_ud,0:kMax)     ! u-stress diagnostic
    REAL(KIND=r8),INTENT(out):: stress_vd      (points_stress_vd,0:kMax)     ! v-stress diagnostic
    REAL(KIND=r8),INTENT(out):: stress_ud_satn (points_stress_ud_satn,0:kMax)! x saturation stress diag
    REAL(KIND=r8),INTENT(out):: stress_vd_satn (points_stress_vd_satn,0:kMax)! y saturation stress diag

    LOGICAL,INTENT(in):: l_drag(points)       ! whether a point has a
    !                           !non-zero surface stress

    LOGICAL,INTENT(in):: l_fix_gwsatn         ! Switch at vn6.2 - if true then
    !                           ! minor big fixes invoked

    LOGICAL,INTENT(in):: l_gwd_40km           ! Switch at vn6.6 - if true then
    !                           ! GWD not applied above z=40km

    !
    !  These logical switches determine whether a diagnostic is to be
    !  calculated.
    !
    LOGICAL,INTENT(in):: stress_ud_on
    LOGICAL,INTENT(in):: stress_ud_p_on       
    LOGICAL,INTENT(in):: stress_vd_on         
    LOGICAL,INTENT(in):: stress_ud_satn_on    
    LOGICAL,INTENT(in):: stress_vd_satn_on    
    LOGICAL,INTENT(in):: du_dt_satn_on        
    LOGICAL,INTENT(in):: du_dt_satn_p_on      
    LOGICAL,INTENT(in):: dv_dt_satn_on


    !-----------------------------------------------------------------
    ! LOCAL ARRAYS AND SCALARS
    !-----------------------------------------------------------------

    INTEGER :: i,k                  ! loop counters in routine

    INTEGER :: kk,kl,ku             ! level counters in routine

    REAL(KIND=r8) :: x_stress(points,2)! x_stresses (layer boundaries)
    REAL(KIND=r8) :: y_stress(points,2)   ! y_stresses (layer boundaries)

    REAL(KIND=r8) :: x_s_const(points) ! level independent constants for
    REAL(KIND=r8) :: y_s_const(points) ! calculation of critical stresses.

    REAL(KIND=r8) :: dzkrho(points)       ! depth of level k *rho
    !                           ! (for diagnostic calculations)

    REAL(KIND=r8) :: dzb      ! height difference across rho   layer
    REAL(KIND=r8) :: delta_z  ! height difference across theta layer
    REAL(KIND=r8) :: dzt                  ! height difference across 2 theta layers

    REAL(KIND=r8) :: dzu,dzl              ! height differences in half layers

    REAL(KIND=r8) :: ub                  ! u-wind at u layer boundary
    REAL(KIND=r8) :: vb                   ! v-wind at u layer boundary

    REAL(KIND=r8) :: n                   ! brunt_vaisala frequency
    REAL(KIND=r8) :: n_sq                 ! square of brunt_vaisala frequency

    REAL(KIND=r8) :: c_x_stress     ! critical x_stress
    REAL(KIND=r8) :: c_y_stress     ! critical y_stress
    REAL(KIND=r8) :: c_stress       ! critical stress amplitude
    REAL(KIND=r8) :: stress         ! current  stress amplitude
    REAL(KIND=r8) :: r_stress       ! inverse of current  stress amplitude
    REAL(KIND=r8) :: s_stress_sq    ! current  stress amplitude squared
    REAL(KIND=r8) :: stress_factor  ! 1. or 1.000001 when using more
    !                           ! numerically robust option

    REAL(KIND=r8) ::  speedcalc     ! dot product calculation for speed/stress

    REAL(KIND=r8) ::  maxgwdlev     ! set to 40km if l_gwd_40km is true
    !                           ! set to 1e10 if l_gwd_40km is false 
    !                           ! (so that if test is never true)  

    x_stress=0.0_r8
    y_stress=0.0_r8
    x_s_const=0.0_r8
    y_s_const=0.0_r8
    dzkrho=0.0_r8
    dzb     =0.0_r8
    delta_z =0.0_r8
    dzt     =0.0_r8
    dzu     =0.0_r8
    dzl       =0.0_r8
    ub     =0.0_r8
    vb     =0.0_r8

    n     =0.0_r8
    n_sq     =0.0_r8

    c_x_stress        =0.0_r8 
    c_y_stress       =0.0_r8  
    c_stress          =0.0_r8 
    stress           =0.0_r8  
    r_stress         =0.0_r8  
    s_stress_sq       =0.0_r8 
    stress_factor      =0.0_r8
    speedcalc     =0.0_r8
    maxgwdlev     =0.0_r8
    ! function and subroutine calls
    ! none
    du_dt=0.0_r8
    dv_dt=0.0_r8



    du_dt_satn         =0.0_r8
    dv_dt_satn         =0.0_r8
    stress_ud=0.0_r8
    stress_vd=0.0_r8
    stress_ud_satn     =0.0_r8
    stress_vd_satn     =0.0_r8
    !-------------------------------------------------------------------
    !   1. preliminaries
    !-------------------------------------------------------------------

    !
    ! Set critical stress_factor according to switch l_fix_gwsatn
    !
    IF ( L_fix_gwsatn ) THEN
       stress_factor = 1.000001_r8
    ELSE
       stress_factor = 1.0_r8
    END IF

    !
    ! Set maxgwdlev according to l_gwd_40km
    !
    IF ( L_gwd_40km ) THEN
       maxgwdlev = 40000.0_r8
    ELSE
       maxgwdlev = 1.e10_r8
    END IF


    kl=1
    ku=2

    DO i=1,points
       x_stress(i,kl) = s_x_stress(i)
       y_stress(i,kl) = s_y_stress(i)
    END DO

    !
    !  Calculate x_s_const.
    !  This includes all terms in the calculation of the critical
    !  stress that are independent of model level
    !
    DO i=1,points
       IF ( l_drag(i) .AND. kay  /=  0.0_r8 ) THEN
          s_stress_sq = SQRT( s_x_stress(i) * s_x_stress(i)             &
               &                 + s_y_stress(i) * s_y_stress(i) )
          !
          !  2 is from rho average in following loop over levels
          !  fr is passed in from gw_satn and scales the critical stress
          !  so that it is equal to the surface stress if rho, U and N
          !  were the same as in the surface stress calculation
          !
          r_stress = 1.0_r8 / ( 2.0_r8 *  s_stress_sq**3 * nsigma * sd_orog(i)  &
               &                      * nsigma * sd_orog(i) * fr(i) * fr(i)     )
          x_s_const(i) = s_x_orog(i) * r_stress
          y_s_const(i) = s_y_orog(i) * r_stress
       ELSE
          x_s_const(i) = 0.0_r8
          y_s_const(i) = 0.0_r8
       END IF
    END DO



    IF( stress_ud_on .OR. stress_ud_p_on ) THEN
       DO i=1,points
          stress_ud(i,0) = s_x_stress(i)
       END DO
    END IF

    IF( stress_vd_on ) THEN
       DO i=1,points
          stress_vd(i,0) = s_y_stress(i)
       END DO
    END IF

    IF( stress_ud_satn_on ) THEN
       DO i=1,points
          stress_ud_satn(i,0) = s_x_stress(i)
       END DO
    END IF

    IF( stress_vd_satn_on ) THEN
       DO i=1,points
          stress_vd_satn(i,0) = s_y_stress(i)
       END DO
    END IF

    !------------------------------------------------------------------
    !   2  Loop over levels calculating critical stresses at each theta
    !      level and reducing the actual stress to that level if it
    !      exceeds the critical stress. The change in stress
    !      across the model level is then used to determine the
    !      deceleration of the wind.
    !
    !      note: k is looping over rho levels
    !------------------------------------------------------------------

    DO k=2,kMax-1 ! min value of k_top is 2 so first check must
       !                     ! be at the top of that level

       DO i=1,points

          x_stress(i,ku) = x_stress(i,kl)
          y_stress(i,ku) = y_stress(i,kl)

          !
          ! First level for wave breaking test is the first
          ! one above the mountain tops.
          ! Top level for wavebreaking now either 40km or level below model lid. 
          !
          IF ( k>=k_top(i) .AND. r_theta_levels(i,k)<=maxgwdlev ) THEN

             IF  (  ( x_stress(i,kl)  /=  0.0_r8)                           &
                  &         .OR.( y_stress(i,kl)  /=  0.0_r8) )   THEN

                dzl        = r_theta_levels(i,k  ) -   r_rho_levels(i,k  )
                dzu        = r_rho_levels  (i,k+1) - r_theta_levels(i,k  )
                ! ub here actually is ub * dzb. Only divide speedcalc by dzb if
                ! needed in if test below
                ub         = dzu*u(i,k)       + dzl*u(i,k+1)
                vb         = dzu*v(i,k)       + dzl*v(i,k+1)
                speedcalc  = ub*s_x_stress(i) + vb*s_y_stress(i)
                !
                !   All wave stress deposited if dth/dz < 0 or wind more than
                !   pi/2 to direction of surface stress
                !
                IF ( theta(i,k+1)  <=  theta(i,k-1) .OR.                  &
                     &             speedcalc  <=  0.0_r8                 ) THEN

                   x_stress(i,ku) = 0.0_r8
                   y_stress(i,ku) = 0.0_r8

                ELSE

                   dzt        = r_theta_levels(i,k+1) -                    &
                        &                       r_theta_levels(i,k-1)
                   n_sq       = g * ( theta(i,k+1) - theta(i,k-1) )        &
                        &                         / ( theta(i,k) * dzt )
                   n          = SQRT( n_sq )

                   dzb        = r_rho_levels(i,k+1) - r_rho_levels(i,k)
                   speedcalc  = speedcalc / dzb

                   !
                   ! rho here is just the average of the full level values. Should
                   ! interpolate in z if this error wasn't so relatively small!
                   !

                   c_y_stress =  kay*(rho(i,k)+rho(i,k+1))*                &
                        &                        stress_factor * speedcalc**3/n
                   c_x_stress = x_s_const(i)   * c_y_stress
                   c_y_stress = y_s_const(i)   * c_y_stress

                   stress     = x_stress(i,kl) * x_stress(i,kl) +          &
                        &                       y_stress(i,kl) * y_stress(i,kl)
                   c_stress   = c_x_stress     * c_x_stress     +          &
                        &                       c_y_stress     * c_y_stress

                   IF ( stress  >   c_stress ) THEN
                      x_stress(i,ku) = x_stress(i,kl)                       &
                           &                            *SQRT(c_stress/stress)
                      y_stress(i,ku) = y_stress(i,kl)                       &
                           &                            *SQRT(c_stress/stress)
                   ELSE
                      x_stress(i,ku) = x_stress(i,kl)
                      y_stress(i,ku) = y_stress(i,kl)
                   END IF!      stress  <   c_stress

                END IF    !  n_sq > 0 and speedcalc > 0

                !------------------------------------------------------------------
                ! Calculate drag from vertical stress convergence
                ! Note that any stress drop at the first level above the mountains
                ! is applied uniformly from that level down to the ground.
                !------------------------------------------------------------------
                IF( k  ==  k_top(i) ) THEN
                   delta_z = r_theta_levels(i,k)
                   IF ( L_fix_gwsatn ) THEN
                      delta_z = delta_z * rho_s(i)/ rho(i,k)
                   END IF
                ELSE
                   delta_z = r_theta_levels(i,k) - r_theta_levels(i,k-1)
                END IF

                du_dt(i,k) = (x_stress(i,ku) - x_stress(i,kl))            &
                     &                    /(delta_z*rho(i,k))
                dv_dt(i,k) = (y_stress(i,ku) - y_stress(i,kl))            &
                     &                    /(delta_z*rho(i,k))
                !PRINT*, k,du_dt(i,k),dv_dt(i,k),delta_z
             END IF    ! stress x or y ne 0

          END IF    ! k ge k_top

       END DO

       ! diagnostics
       IF( stress_ud_on ) THEN
          DO i=1,points
             IF ( k  >=  k_top(i)   ) THEN
                stress_ud(i,k) = x_stress(i,ku)
             END IF
          END DO
       END IF

       IF( stress_vd_on ) THEN
          DO i=1,points
             IF ( k  >=  k_top(i)   ) THEN
                stress_vd(i,k) = y_stress(i,ku)
             END IF
          END DO
       END IF

       IF( stress_ud_satn_on ) THEN
          DO i=1,points
             IF ( k  >=  k_top(i)   ) THEN
                stress_ud_satn(i,k) = x_stress(i,ku)
             END IF
          END DO
       END IF

       IF( stress_vd_satn_on ) THEN
          DO i=1,points
             IF ( k  >=  k_top(i)   ) THEN
                stress_vd_satn(i,k) = y_stress(i,ku)
             END IF
          END DO
       END IF

       IF( du_dt_satn_on ) THEN
          DO i=1,points
             du_dt_satn(i,k) = du_dt(i,k)
          END DO
       END IF

       IF( dv_dt_satn_on ) THEN
          DO i=1,points
             dv_dt_satn(i,k) = dv_dt(i,k)
          END DO
       END IF

       ! swap storage for lower and upper layers
       kk=kl
       kl=ku
       ku=kk

    END DO
    !   end loop kMax

    !
    ! Top level diagnostics - assume now that no drag is applied at
    ! the top level
    !
    IF( stress_ud_on .OR. stress_ud_p_on ) THEN
       DO i=1,points
          stress_ud(i,kMax) = x_stress(i,kl)
       END DO
    END IF

    IF( stress_vd_on ) THEN
       DO i=1,points
          stress_vd(i,kMax) = y_stress(i,kl)
       END DO
    END IF

    IF( stress_ud_satn_on ) THEN
       DO i=1,points
          stress_ud_satn(i,kMax) = x_stress(i,kl)
       END DO
    END IF

    IF( stress_vd_satn_on ) THEN
       DO i=1,points
          stress_vd_satn(i,kMax) = y_stress(i,kl)
       END DO
    END IF

    !--------------------------------------------------------------------
    ! 3.1 set drags and stresses below k_top
    !     can only be Done now because need drag at k_top first
    !--------------------------------------------------------------------
    DO k=1,k_top_max-1
       DO i=1,points
          IF ( k  <   k_top(i) ) THEN
             du_dt(i,k) = du_dt(i,k_top(i))
             dv_dt(i,k) = dv_dt(i,k_top(i))
          END IF
       END DO
    END DO

    !
    ! saturation drag diagnostic
    !
    IF( du_dt_satn_on .OR. du_dt_satn_p_on ) THEN
       DO k=1,k_top_max-1
          DO i=1,points
             IF ( k  <   k_top(i) ) THEN
                du_dt_satn(i,k) = du_dt_satn(i,k_top(i))
             END IF
          END DO
       END DO
    END IF

    IF( dv_dt_satn_on ) THEN
       DO k=1,k_top_max-1
          DO i=1,points
             IF ( k  <   k_top(i) ) THEN
                dv_dt_satn(i,k) = dv_dt_satn(i,k_top(i))
             END IF
          END DO
       END DO
    END IF

    !
    !  stress diagnostics
    !
    IF( stress_ud_on .OR. stress_ud_p_on .OR. stress_ud_satn_on .OR.  &
         &    stress_vd_on .OR. stress_vd_satn_on) THEN

       DO k=1,k_top_max-1
          DO i=1,points
             IF ( k  <   k_top(i) ) THEN
                IF ( k  ==  1 ) THEN
                   dzkrho(i) =  r_theta_levels(i,1)* rho(i,k)
                ELSE
                   dzkrho(i) = (r_theta_levels(i,k)                      &
                        &                        -r_theta_levels(i,k-1))* rho(i,k)
                END IF
             END IF
          END DO

          IF( stress_ud_on ) THEN
             DO i=1,points
                IF ( k  <   k_top(i) ) THEN
                   stress_ud(i,k) = stress_ud(i,k-1)+du_dt(i,k)*dzkrho(i)
                END IF
             END DO
          END IF

          IF( stress_ud_satn_on ) THEN
             DO i=1,points
                IF ( k  <   k_top(i) ) THEN
                   stress_ud_satn(i,k)= stress_ud_satn(i,k-1)              &
                        &                              + du_dt(i,k) * dzkrho(i)
                END IF
             END DO
          END IF

          IF( stress_vd_on ) THEN
             DO i=1,points
                IF ( k  <   k_top(i) ) THEN
                   stress_vd(i,k)= stress_vd(i,k-1) + dv_dt(i,k) *dzkrho(i)
                END IF
             END DO
          END IF

          IF( stress_vd_satn_on ) THEN
             DO i=1,points
                IF ( k  <   k_top(i) ) THEN
                   stress_vd_satn(i,k)= stress_vd_satn(i,k-1)              &
                        &                              + dv_dt(i,k) * dzkrho(i)
                END IF
             END DO
          END IF

       END DO ! k=1,k_top_max-1

    END IF   ! stress_ud_on .or. stress_ud_p_on .or. stess_ud_satn_on .or.
    !              ! stress_vd_on .or. stess_vd_satn_on




    RETURN
  END SUBROUTINE gw_satn
  !===============================================================================
  SUBROUTINE geopotential_t(                                 &
       piln   ,  pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   ,nCols, kMax)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute the geopotential height (above the surface) at the midpoints and 
    ! interfaces using the input temperatures and pressures.
    !
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncol                  ! Number of longitudes
    INTEGER, INTENT(in) :: nCols
    INTEGER, INTENT(in) :: kMax
    REAL(r8), INTENT(in) :: piln (nCols,kMax+1)   ! Log interface pressures
    REAL(r8), INTENT(in) :: pint (nCols,kMax+1)   ! Interface pressures
    REAL(r8), INTENT(in) :: pmid (nCols,kMax)    ! Midpoint pressures
    REAL(r8), INTENT(in) :: pdel (nCols,kMax)    ! layer thickness
    REAL(r8), INTENT(in) :: rpdel(nCols,kMax)    ! inverse of layer thickness
    REAL(r8), INTENT(in) :: t    (nCols,kMax)    ! temperature
    REAL(r8), INTENT(in) :: q    (nCols,kMax)    ! specific humidity
    REAL(r8), INTENT(in) :: rair                 ! Gas constant for dry air
    REAL(r8), INTENT(in) :: gravit               ! Acceleration of gravity
    REAL(r8), INTENT(in) :: zvir                 ! rh2o/rair - 1

    ! Output arguments

    REAL(r8), INTENT(out) :: zi(nCols,kMax+1)     ! Height above surface at interfaces
    REAL(r8), INTENT(out) :: zm(nCols,kMax)      ! Geopotential height at mid level
    !
    !---------------------------Local variables-----------------------------
    !
    LOGICAL  :: fvdyn              ! finite volume dynamics
    INTEGER  :: i,k                ! Lon, level indices
    REAL(r8) :: hkk(nCols)         ! diagonal element of hydrostatic matrix
    REAL(r8) :: hkl(nCols)         ! off-diagonal element
    REAL(r8) :: rog                ! Rair / gravit
    REAL(r8) :: tv                 ! virtual temperature
    REAL(r8) :: tvfac              ! Tv/T
    !
    !-----------------------------------------------------------------------
    !
    rog = rair/gravit

    ! Set dynamics flag

    fvdyn = .FALSE.!dycore_is ('LR')

    ! The surface height is zero by definition.

    DO i = 1,ncol
       zi(i,kMax+1) = 0.0_r8
    END DO

    ! Compute zi, zm from bottom up. 
    ! Note, zi(i,k) is the interface above zm(i,k)

    DO k = kMax, 1, -1

       ! First set hydrostatic elements consistent with dynamics

       IF (fvdyn) THEN
          DO i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1.0_r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          END DO
       ELSE
          DO i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          END DO
       END IF

       ! Now compute tv, zm, zi

       DO i = 1,ncol
          tvfac   = 1.0_r8 + zvir * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv * hkl(i)
       END DO
    END DO

    RETURN
  END SUBROUTINE geopotential_t

END MODULE GwddSchemeCPTEC
