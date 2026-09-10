MODULE GwddSchemeUSSP
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

  ! Pi
  REAL(KIND=r8), PARAMETER :: Pi                 = 3.14159265358979323846_r8

  REAL(KIND=r8),PARAMETER:: LSTAR = 1.0_r8/4300.0_r8!1/characteristic wavelength
  REAL(KIND=r8),PARAMETER:: LMINL = 1.0_r8/20000.0_r8!1/max wavelength at launch
  REAL(KIND=r8),PARAMETER:: SigmaLAUNCH = 0.045_r8! SigmaLAUNCH = eta for model launch

  !
  !     Maximum number of iterations of Newton Raphson DO (While) loop
  INTEGER, PARAMETER :: MAXWHILE       = 9

  !
  !     Wavenumber at peak in spectrum
  REAL(KIND=r8), PARAMETER ::  MSTAR            = 2.0_r8 * PI * LSTAR
  !
  !     Reciprocal of mstar (m) and mstar^2 (m^2)
  REAL(KIND=r8), PARAMETER ::  RMSTAR           = 1.0_r8 / MSTAR
  REAL(KIND=r8), PARAMETER ::  RMSTARSQ         = RMSTAR * RMSTAR

  !
  !     Normalised minimum vertical wavenumber at launch
  REAL(KIND=r8), PARAMETER ::  MNLMIN           = LMINL / LSTAR

  !
  !     Equatorial planetary vorticity gradient parameter B_eq (/m /s )
  REAL(KIND=r8), PARAMETER :: BETA_EQ_RMSTAR    = 2.3E-11_r8 * RMSTAR
  !
  !  Common/Omega/Omega, two_omega
  REAL(KIND=r8), PARAMETER  :: Omega=7.27e-5_r8   ! rad/s
  REAL(KIND=r8), PARAMETER  :: two_omega=2*Omega
  !  Angular speed of Earth's rotation Omega to be initialised in SETCON

  ! Mean radius of Earth in m.
  REAL(KIND=r8), PARAMETER  :: Earth_Radius = 6371229.0_r8
  !*----------------------------------------------------------------------
  !*L------------------COMDECK C_G----------------------------------------
  ! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

  REAL(KIND=r8), PARAMETER :: G = 9.80665_r8
  REAL(KIND=r8), PARAMETER :: R_EARTH_RADIUS_SQ =1.0_r8 / (Earth_Radius*Earth_Radius)
  !
  ! ----------------------------------------------------------------------+-------
  !     Security parameters
  ! ----------------------------------------------------------------------+-------
  !
  !     Minimum allowed value of buoyancy frequency squared
  REAL(KIND=r8), PARAMETER ::  SQNMIN           = 1.0E-4_r8
  !
  ! ----------------------------------------------------------------------+-------
  !     Local parameters
  ! ----------------------------------------------------------------------+-------
  !
  !     Max number of directions, typically four.
  INTEGER, PARAMETER :: IDIR           = 4
  !
  !     Strength coefficient constant for Launch spectrum (CCL / A0)
  REAL(KIND=r8), PARAMETER :: CCL0 = 3.41910625e-9_r8
  !
  !
  !     Parameter beta in the launch spectrum total energy equation
  REAL(KIND=r8), PARAMETER :: BETA_E0           = 1.0227987125E-1_r8
  !
  !
  !     Azimuthal sector for launch spectrum integral Delta Phi / 2
  REAL(KIND=r8), PARAMETER :: DDPHIR2           = PI / IDIR
  !
  !     Parameter p in B_0(p) for launch spectrum intrinsic frequency
  !     NOTE: This parameter determines the intrinsic frequency spectrum
  !           shape and hence the integral form in 4.1, which is strictly
  !           valid only for p > 1. !!IF contemplating changes BE WARNED!!
  REAL(KIND=r8), PARAMETER :: PSAT              = 5.0_r8 / 3.0_r8
  !
  !     Psat - 1
  REAL(KIND=r8), PARAMETER :: PSATM1            = PSAT - 1.0_r8
  !
  !     2 - Psat
  REAL(KIND=r8), PARAMETER :: TWOMPSAT          = 2.0_r8 - PSAT
  !
  !     Power s of vertical wavenumber spectrum A_0(s,t) at low m
  REAL(KIND=r8), PARAMETER :: SS                = 1.0_r8
  !
  !     s + 1, s - 1
  REAL(KIND=r8), PARAMETER :: SSP1              = SS + 1.0_r8
  !
  !
  !     Power t=t_sat of vertical wavenumber spectrum at large m due to
  !     saturation by gravity wave breaking (and shape of chopping fn)
  REAL(KIND=r8), PARAMETER :: TT                = 3.0_r8
  !
  !     t - 1, t - 2, 1 / (t-2), (t-3) / (t-2), 2 - t
  REAL(KIND=r8), PARAMETER :: TTM1              = TT - 1.0_r8
  REAL(KIND=r8), PARAMETER :: TTM2              = TT - 2.0_r8
  REAL(KIND=r8), PARAMETER :: RTTM2              = 1.0_r8 / TTM2
  REAL(KIND=r8), PARAMETER :: TTRAT              = (TT - 3.0_r8) * RTTM2
  REAL(KIND=r8), PARAMETER :: TWOMTT              = 2.0_r8 - TT
  !
  !     s + t, 1 / (s+t)
  REAL(KIND=r8), PARAMETER :: SSPTT             = SS + TT
  REAL(KIND=r8), PARAMETER :: RSSPTT            = 1.0_r8 / SSPTT
  !
  !     Weight for (n+1)th guess at mNlX in iteration solution
  REAL(KIND=r8), PARAMETER :: MWEIGHT           = 0.8_r8

  !
  !     Minimum allowed non-zero value of Curvature Coefficient A
  REAL(KIND=r8), PARAMETER ::  ASECP            =  1.0E-20_r8
  REAL(KIND=r8), PARAMETER ::  ASECN            = -(ASECP)

  !
  ! ==Main Block==--------------------------------------------------------+-------
  ! --------------------------------
  !     Local Constants (a) Physical
  ! --------------------------------
  REAL(KIND=r8), PARAMETER :: COSPHI(1:4)=(/0.0_r8,-1.0_r8,  0.0_r8,1.0_r8/)
  REAL(KIND=r8), PARAMETER :: SINPHI(1:4)=(/1.0_r8, 0.0_r8 ,-1.0_r8,0.0_r8/)

  LOGICAL  , PARAMETER :: L_USSP_OPAQUE=.FALSE.  !IN  Switch for Opaque Upper boundary condition

  REAL(KIND=r8) :: akappa

  PUBLIC :: Init_GwUSSP
  PUBLIC :: gw_ussp
CONTAINS

  SUBROUTINE Init_GwUSSP(       )
    INTEGER :: k
    akappa=gasr/cp
    !
  END SUBROUTINE Init_GwUSSP

  SUBROUTINE gw_ussp( &
       kMax       , &
       nCols      , &
       prsi ,prsl  ,phii ,phil    ,&
       gt         , &
       gq         , &
       gu         , &
       gv         , &
       topo       , &
       colrad     , &
       G_X        , &
       G_Y          &
       )

    INTEGER, INTENT(IN   ) :: kMax
    INTEGER, INTENT(IN   ) :: nCols
    REAL(KIND=r8), INTENT(in   ) :: prsi   (ncols,kMax+1)  !     prsi     - real, pressure at layer interfaces [Pa]
    REAL(KIND=r8), INTENT(in   ) :: prsl   (ncols,kMax)    !     prsl     - real, mean layer presure [Pa]
    REAL(KIND=r8), INTENT(in   ) :: phii   (nCols,kMax+1) !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    REAL(KIND=r8), INTENT(in   ) :: phil   (nCols,kMax)   !===>  PHIL(K) INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(KIND=r8), INTENT(IN   ) :: gt    (nCols,kMax)       !REAL(r8), INTENT(in) :: gt (nCols,kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: gq    (nCols,kMax)       !REAL(r8), INTENT(in) :: gq (nCols,kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: gu    (nCols,kMax)       !REAL(r8), INTENT(in) :: gu (nCols,kMax)  
    REAL(KIND=r8), INTENT(IN   ) :: gv    (nCols,kMax)       !REAL(r8), INTENT(in) :: gv (nCols,kMax)   
    REAL(KIND=r8), INTENT(IN   ) :: topo  (nCols)
    REAL(KIND=r8), INTENT(IN   ) :: colrad(nCols)       
    REAL(KIND=r8), INTENT(OUT  ) :: G_X   (nCols,kMax) !Zonal component of wave-induced
    !force at longitude, level (m s^-2)
    !Note that in our notation G_X is
    !equivalent to DU_DT, the zonal
    !wind tendency (m s^-2)
    REAL(KIND=r8), INTENT(OUT  ) :: G_Y(nCols,kMax)!Meridional component of wave-induced
    !force at longitude, level (m s^-2)
    !Note that in our notation G_Y is
    !equivalent to DV_DT, the meridional
    !wind tendency (m s^-2)

    REAL(KIND=r8) :: SIN_THETA_LATITUDE(nCols)       ! P-GRID Latitudes
    REAL(KIND=r8) :: P_LAYER_BOUNDARIES(nCols,0:kMax)! Pa
    REAL(KIND=r8) :: THETA             (nCols,kMax)  ! IN Primary model array for theta primary theta field (K)
    REAL(KIND=r8) :: RHONT             (nCols,kMax)  ! Rho on Theta grid.
    REAL(KIND=r8) :: R_THETA_LEVELS    (nCols,0:kMax)! Distance of theta levels from Earth centre.
    REAL(KIND=r8) :: UONP(nCols,kMax)                ! U on Rho grid for local use  primary u field (ms**-1)
    REAL(KIND=r8) :: VONP(nCols,kMax)                ! V on Rho grid for local use  primary v field (ms**-1)

    REAL(KIND=r8) :: RHO_TH(nCols,kMax)! Rho on theta levels

    REAL(KIND=r8) :: NBV(nCols,kMax)! Buoyancy frequency
    !                             [Brunt Vaisala frequency] on half-levels
    !(rad s^-1)
    REAL(KIND=r8) :: UDOTK(nCols,kMax+1,IDIR)! Component of wind
    !                            in phi_jdir direction (m s^-1)
    REAL(KIND=r8) :: FPTOT(nCols,kMax,IDIR)! Pseudomomentum flux
    !                            integrated over
    !azimuthal sector (kg m^-1 s^-2)
    REAL(KIND=r8) :: G_G  (nCols,kMax,IDIR)!Wave-induced force per unit
    !mass due to azimuthal sectors (m s^-2)
    REAL(KIND=r8) :: RHOCL(nCols,IDIR)! [Rho . Cl]_klaunch
    REAL(KIND=r8) :: OMIN(nCols)! Either f_f or the equatorial minimum
    !                            frequency, whichever is less  (rad s^-1)
    REAL(KIND=r8) :: OMINRNBV    ! omega_min(launch) / N (k)
    REAL(KIND=r8) :: RHOCSK(nCols,kMax)! [Rho(z) . Csat(z)]_k
    REAL(KIND=r8) :: MGUESS_a(nCols,kMax+1,IDIR)! Starting value of vertical wavenumber for
    !                            crossing point search
    REAL(KIND=r8) :: DDU_a(nCols,kMax,IDIR)
    !                            Delta U=udotk(launch)-udotk(jlev)
    REAL(KIND=r8) :: ATTENUATION(nCols,kMax,IDIR)                        
    !                            Coefficient B in intersect point equation
    REAL(KIND=r8) :: CURVATURE(nCols,kMax,IDIR)                          
    !                            Term (A / B) in intersect point equation
    REAL(KIND=r8) :: ACOEFF(nCols,kMax,IDIR)                             
    !                            Coefficient A in intersect point equation
    REAL(KIND=r8) :: MINTERCEPT(nCols,kMax,IDIR)                         
    !                            Chop function B*[1 + (A/B)*mNlmin]^(t-2)
    REAL(KIND=r8) :: INTERCEPT1(nCols,kMax,IDIR)
    !                            Chop function B*[1 + (A/B)]^(t-2)
    REAL(KIND=r8) :: FPFAC(nCols,kMax,IDIR) ! Record of total flux
    REAL(KIND=r8) :: INDEXI(IDIR*nCols)! I location of chop type points
    !REAL(KIND=r8) :: INDEXJ(IDIR*nCols)! J location of chop type points
    REAL(KIND=r8) :: INDEXJD(IDIR*nCols)! JDIR location of chop type pnts
    REAL(KIND=r8) :: ATTE_C(IDIR*nCols)! Compressed attenuation array
    REAL(KIND=r8) :: CURV_C(IDIR*nCols)! Compressed curvature array
    REAL(KIND=r8) :: MNLX(IDIR*nCols,0:MAXWHILE)! Intersect mNlX estimates

    LOGICAL :: L_CHOP2             ! Indicates spectra with low m intersect

    INTEGER :: k
    INTEGER :: i
    !INTEGER :: j
    INTEGER :: JLEV
    INTEGER :: JDIR
    INTEGER :: ILAUNCH (nCols)  !Minimum level for launch level (all points)
    REAL(KIND=r8) ::  F_F!Inertial frequency at current latitude
    !(rad s^-1)
    REAL(KIND=r8) ::  MKILL! Wavenumber where Doppler transformed
    !                             spectrum is reduced to zero
    REAL(KIND=r8) ::  MNLY! High wavenumber intersection point
    REAL(KIND=r8) ::  TAIL_CHOP2B! Integral segment [(s+1) * mNLmin**(1-t)]
    REAL(KIND=r8) ::  HEAD_CHOP2A! Integral segment [(t-1) * mNLmin**(s+1)]
    REAL(KIND=r8) ::  A0_R_SP1TM1!A0/(s+1)(t-1) normalisation factor for the
    !                          ! launch spectrum vertical wavenumber
    REAL(KIND=r8) ::  A0_R_1MT   !-A0/(t-1) norm factor for launch spectrum m
    REAL(KIND=r8) :: CCS0        ! Constant component of saturation spectrum  
    REAL(KIND=r8) :: FPLUS ! Maximum range value of function f
    REAL(KIND=r8) :: GPLUS  ! Maximum range value of function g
    REAL(KIND=r8) :: FMINUS ! Minimum range value of function f
    REAL(KIND=r8) :: GMINUS ! Minimum range value of function g
    REAL(KIND=r8) :: FTERM  ! Intermediate  value of function f
    REAL(KIND=r8) :: GTERM  ! Intermediate  value of function g
    LOGICAL        :: L_FTHENG(IDIR*nCols) ! Indicate dir of spiral solution
    REAL(KIND=r8) :: WGTN(IDIR*nCols)! Weighting of n term in iter
    REAL(kind=r8) :: pmid(ncols,kmax)     ! midpoint pressures Pa

    INTEGER :: LAUNCHLEV (ncols)!Launch level at specific point
    INTEGER :: NCHOP2!Number of spectra with low m intersect
    INTEGER :: NNJD!Index values
    INTEGER :: NNI!Index values
    INTEGER :: JWHILE!Counter for while loop
    REAL(r8) :: state_t     (nCols,kMax)  
    REAL(r8) :: state_q     (nCols,kMax)  
    REAL(r8) :: state_u     (nCols,kMax)  
    REAL(r8) :: state_v     (nCols,kMax)  
    REAL(r8) :: state_pmid  (nCols,kMax)  
    REAL(r8) :: state_pint  (nCols,kMax+1)  
    REAL(r8) :: state_lnpint(nCols,kMax+1)  
    REAL(r8) :: state_pdel  (nCols,kMax)  
    REAL(r8) :: state_rpdel (nCols,kMax)
    REAL(r8) :: state_lnpmid(nCols,kMax)
    REAL(r8) :: state_zi        (1:nCols,1:kMax+1)   
    REAL(r8) :: state_zm        (1:nCols,1:kMax)
    REAL(KIND=r8)  :: sigkiv    (1:nCols,1:kMax)

    !-----------------------------------------------------------------------------
    DO i=1,nCols
       !state_pint       (i,kMax+1) = gps(i)*si(1)
       state_pint       (i,kMax+1) = prsi(i,1)
    END DO
    DO k=kMax+1,1,-1
       DO i=1,nCols
          !state_pint    (i,k)      = MAX(si(kMax+2-k)*gps(i) ,0.0001_r8)
          state_pint    (i,k) = MAX(prsi(i,kMax+2-k) ,0.0001_r8)
       END DO
    END DO


    DO k=1,kMax
       DO i=1,nCols
          state_t (i,kMax+1-k) =  gt (i,k)
          state_q (i,kMax+1-k) =  gq (i,k)
          state_u (i,kMax+1-k) =  gu (i,k)
          state_v (i,kMax+1-k) =  gv (i,k)
         ! state_pmid(i,kMax+1-k) = sl(k)*gps (i)
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

    !
    !..delsig     k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(2)
    !             k=3/2----si,ric,rf,km,kh,b,l -----------
    !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sl*** } delsig(1)
    !             k=1/2----si ----------------------------

    ! Derive new temperature and geopotential fields
!    PRINT*,'geopotential_t'
    CALL geopotential_t(                                 &
         state_lnpint(1:nCols,1:kMax+1)   , state_pint (1:nCols,1:kMax+1)  , &
         state_pmid  (1:nCols,1:kMax)     , state_pdel  (1:nCols,1:kMax)   , state_rpdel(1:nCols,1:kMax)   , &
         state_t     (1:nCols,1:kMax)     , state_q     (1:nCols,1:kMax)   , rair   , gravit , zvir   ,&
         state_zi    (1:nCols,1:kMax+1)   , state_zm    (1:nCols,1:kMax)   , nCols   ,nCols, kMax)


    FMINUS      = MNLMIN**SSPTT 
    TAIL_CHOP2B = SSP1 / (MNLMIN**TTM1)
    HEAD_CHOP2A = TTM1 * (MNLMIN**SSP1) 
    A0_R_SP1TM1 = 1.0_r8 / ( SS + TT - HEAD_CHOP2A )
    A0_R_1MT    = (-(SSP1) ) * A0_R_SP1TM1
    CCS0        = BETA_E0 * SIN(DDPHIR2) * PSATM1 / (PI * TWOMPSAT)

  
    !-----------------------------------------------------------------------------

    DO i=1,nCols
       SIN_THETA_LATITUDE(i)       = SIN(colrad(i))  ! P-GRID Latitudes
       P_LAYER_BOUNDARIES       (i,0) =  prsi(i,1)   !g ps(i)*si(1)
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
          P_LAYER_BOUNDARIES   (i,k) =prsi(i,k+1) ! MAX(si(k+1)*gps(i) ,0.0001_r8)
          R_THETA_LEVELS       (i,k) = R_THETA_LEVELS(i,0) + state_zi    (i,kMax+1-k)
          THETA                (i,k) = sigkiv(i,k)*gt(i,k)
          UONP                 (i,k) = gu(i,k)          !U on Rho grid for local use
          VONP                 (i,k) = gv(i,k)          !V on Rho grid for local use
          !pmid                 (i,k) = gps(i)*sl(k)
          pmid                 (i,k) = prsl(i,k)

          !j/kg/kelvin
          !
          ! P = rho * R * T
          !
          !            P
          ! rho  = -------
          !          R * T
          !
          !           1            R * T
          !  rrho = ----- =   -------
          !          rho             P
          !
          !
          RHONT(i,k) = pmid(i,k)/(gasr*gt(i,k))
       END DO
    END DO
! PRINT*,'Find model level for launch'
    !
    ! ----------------------------------------------------------------------+-------
    ! Find model level for launch
    ! ----------------------------------------------------------------------+-------
    !     Minimum value if variable height sources are used

    ILAUNCH=1
    DO k=1,kMax-1
       DO i=1,nCols
          ! SigmaLAUNCH = 0.045_r8! SigmaLAUNCH = eta for model launch

          !IF (sl(k) > SigmaLAUNCH .AND. sl(k+1) <= SigmaLAUNCH) THEN
          IF ((prsl(i,k)/prsi(i,1))> SigmaLAUNCH.AND. (prsl(i,k+1)/prsi(i,1)) <= SigmaLAUNCH) THEN
             !IF (SigmaLAUNCH-sl(k) > sl(k+1)-SigmaLAUNCH) THEN
             IF (SigmaLAUNCH-(prsl(i,k)/prsi(i,1)) > (prsl(i,k+1)/prsi(i,1))-SigmaLAUNCH) THEN
                ILAUNCH(i) = k
             ELSE
                ILAUNCH(i) = k + 1
             END IF
          END IF
      END DO
    END DO
    LAUNCHLEV = ILAUNCH

    Levels_do1a: DO K=1,kMax
       nCols_do1a: DO I=1,nCols
          ! ----------------------------------------------------------------------+-------
          !           Zero vertical divergence of pseudomomentum flux
          ! ----------------------------------------------------------------------+-------
          G_X(I,K) = 0.0_r8
          G_Y(I,K) = 0.0_r8
       END DO  nCols_do1a
    END DO  Levels_do1a
    !
    ! ----------------------------------------------------------------------+-------
    ! 1.0   Set variables that are to be defined on all model levels
    ! ----------------------------------------------------------------------+-------
    !
    Levels_do2: DO JLEV=2,(kMax - 1)
       ! ----------------------------------------------------------------------+-------
       ! 1.1   Density, buoyancy frequency and altitude for middle levels
       ! ----------------------------------------------------------------------+-------
       nCols_do2: DO I=1,nCols

          RHO_TH(I,JLEV) = RHONT(I,JLEV) * r_earth_radius_sq

          !           Buoyancy (Brunt-Vaisala) frequency calculation
          NBV(I,JLEV) = ( g*(THETA(I,JLEV+1)-THETA(I,JLEV-1))/  &
               &                          (THETA(I,JLEV) *                      &
               &       (R_THETA_LEVELS(I,JLEV+1)-R_THETA_LEVELS(I,JLEV-1))) )
          NBV(I,JLEV) = MAX(NBV(I,JLEV), SQNMIN)
          NBV(I,JLEV) = SQRT(NBV(I,JLEV))

       END DO  nCols_do2
    END DO  Levels_do2
! PRINT*,'Density, Buoyancy (Brunt-Vaisala) frequency at top and bottom'

    !
    ! ----------------------------------------------------------------------+-------
    !     Density, Buoyancy (Brunt-Vaisala) frequency at top and bottom
    ! ----------------------------------------------------------------------+-------
    nCols_do3: DO I=1,nCols
       !       Set rho and theta level value of Z_TH.
       RHO_TH(I,1)      = RHONT(I,1) * r_earth_radius_sq
       NBV(I,1)         = NBV(I,2)
       RHO_TH(I,kMax) = RHONT(I,kMax) * r_earth_radius_sq
       NBV(I,kMax)    = NBV(I,kMax-1)
    END DO  nCols_do3
    !
    Levels_do4: DO JLEV=kMax,2,-1
       ! ----------------------------------------------------------------------+-------
       ! 1.2   Set buoyancy frequency constant up to 1km altitude
       ! ----------------------------------------------------------------------+-------
       nCols_do4: DO I=1,nCols
          IF ( (R_THETA_LEVELS(I,JLEV) - EARTH_RADIUS) <  1.0E3_r8)    & !????????
               &      NBV(I,JLEV-1) = NBV(I,JLEV)
       END DO  nCols_do4
    END DO  Levels_do4
    !
    !
    ! ----------------------------------------------------------------------+-------
    ! 1.3  Compute component of wind U in each wave-propagation direction.
    !      U is the dot product of (u,v) with k_0 but n.b. UDOTK is half-way
    !      between rho levels.
    ! ----------------------------------------------------------------------+-------
    IDir_do1: DO JDIR=1,IDIR
       !
       Levels_do5: DO JLEV=1,kMax-1
          nCols_do5: DO I=1,nCols
             !           Assume theta levels are half way between rho levels.
             UDOTK(I,JLEV,JDIR) = &
                  &        0.5_r8*(UONP(I,JLEV) + UONP(I,JLEV+1))*COSPHI(JDIR) + &
                  &        0.5_r8*(VONP(I,JLEV) + VONP(I,JLEV+1))*SINPHI(JDIR)
          END DO  nCols_do5
       END DO  Levels_do5
       !
       ! ----------------------------------------------------------------------+-------
       !      Set wind component for top level, to be equal to that on the top
       !      Rho level, and total flux of horizontal pseudomomentum at bottom
       ! ----------------------------------------------------------------------+-------
       nCols_do5a: DO I=1,nCols
          UDOTK(I,kMax,JDIR) =                                    &
               &       UONP(I,kMax)*COSPHI(JDIR) + VONP(I,kMax)*SINPHI(JDIR)
          FPTOT(I,1,JDIR) = 0.0_r8
       END DO  nCols_do5a
       !
       ! ----------------------------------------------------------------------+-------
       ! 2.0  Initialize variables that need to be defined up to launch level
       ! ----------------------------------------------------------------------+-------
       !
       Levels_do6: DO JLEV=1,kMax
          nCols_do6: DO I=1,nCols
             ! ----------------------------------------------------------------------+-------
             !           Vertical divergence of total flux of horizontal pseudomom
             ! ----------------------------------------------------------------------+-------
             IF(JLEV >=1.and. JLEV <= ILAUNCH(i))THEN
                 G_G(I,JLEV,JDIR) = 0.0_r8
             END IF
          END DO  nCols_do6
       END DO  Levels_do6
       !
       ! ----------------------------------------------------------------------+-------
       ! 3.0  Initialize gravity wave spectrum variables for the launch level
       ! ----------------------------------------------------------------------+-------
       !
       nCols_do7: DO I=1,nCols
          ! ----------------------------------------------------------------------+-------
          !         Value at launch of (m*^2 * Total vertical flux of horizontal
          !         pseudomomentum ... UMDP 34: 1.14)
          ! ----------------------------------------------------------------------+-------
          RHOCL(I,JDIR) = RHO_TH(I,ILAUNCH(i)) * CCL0
       END DO  nCols_do7

    END DO  IDir_do1
    !
    ! ----------------------------------------------------------------------+-------
    ! 3.1 Compute minimum intrinsic frequency OMIN from inertial frequency
    !     squared F_F. See UMDP 34, eqn. (A.7).
    !-----------------------------------------------------------------------+-------
    nCols_do7a: DO I=1,nCols
       F_F = (TWO_OMEGA * SIN_THETA_LATITUDE(I))**2
       OMIN(I) = NBV(I,ILAUNCH(i)) * BETA_EQ_RMSTAR
       OMIN(I) = MAX( OMIN(I), F_F )
       OMIN(I) = SQRT(OMIN(I))
    END DO  nCols_do7a
    !
    ! ----------------------------------------------------------------------+-------
    ! 4.0 Calculations carried out at levels from Launch to Lid
    ! ----------------------------------------------------------------------+-------
    !
    !
    Levels_do8: DO JLEV=2,kMax
       nCols_do8: DO I=1,nCols
          IDir_do1a: DO JDIR=1,IDIR
             ! ----------------------------------------------------------------------+-------
             !           Total vertical flux of horizontal pseudomomentum (analytic
             !           integral under curve = rho_l * C_l / m*^2 ... UMDP 34: 1.14)
             ! ----------------------------------------------------------------------+-------
             IF (JLEV == LAUNCHLEV(i))  THEN
                FPTOT(I,JLEV,JDIR) = RMSTARSQ * RHOCL(I,JDIR)
             ELSE
                FPTOT(I,JLEV,JDIR) = 0.0_r8
             END IF
             FPTOT(I,JLEV,JDIR) = FPTOT(I,JLEV-1,JDIR) +           &
                  &                               FPTOT(I,JLEV,JDIR)
          END DO  IDir_do1a
       END DO  nCols_do8
    END DO  Levels_do8
    !
    ! ----------------------------------------------------------------------+-------
    ! 4.1 Compute [rho(z) . C(z)]_k for combined spectrum ...(UMDP 34: 1.16)
    ! ----------------------------------------------------------------------+-------
    !
    !
    !     IF (ABS(PSATM1) >= 0.1) THEN
    !     For current setting of parameter psat this test is always true
    Levels_do8a: DO JLEV=1,kMax
       nCols_do8a: DO I=1,nCols
         IF(JLEV >= ILAUNCH(i) .AND.JLEV <= kMax )THEN
          OMINRNBV = OMIN(I) / NBV(I,JLEV)
          RHOCSK(I,JLEV) = RHO_TH(I,JLEV) * CCS0 *              &
               &         (NBV(I,JLEV))**2 * (OMINRNBV**PSATM1) *                &
               &         (1.0_r8 - (OMINRNBV**TWOMPSAT)) / (1.0_r8 - (OMINRNBV**PSATM1))
          END IF
       END DO  nCols_do8a
    END DO  Levels_do8a
    !     ELSE
    !     Require a different functional form for normalisation factor B0
    !             BBS = 1.0 / ALOG(NBV(I,J,JLEV) / OMIN(I,J))
    !     END IF
    ! ----------------------------------------------------------------------+-------
    !     Loop over directions and levels and calculate horizontal component
    !     of the vertical flux of pseudomomentum for each azimuthal
    !     direction and for each altitude
    ! ----------------------------------------------------------------------+-------
    Levels_do92: DO JLEV=1,kMax
       L_CHOP2 = .FALSE.
       IDir_do2a: DO JDIR=1,IDIR
          nCols_do92: DO I=1,nCols
             IF( JLEV>=ILAUNCH(i)+1 .and.JLEV<= kMax)THEN
             ! ----------------------------------------------------------------------+-------
             !     Initialise MGUESS (start point for iterative searches if needed)
             ! ----------------------------------------------------------------------+-------
             MGUESS_a(I,JLEV,JDIR) = 0.0_r8
             Fptot_if1: IF (FPTOT(I,JLEV-1,JDIR) >  0.0_r8) THEN
                ! ----------------------------------------------------------------------+-------
                ! 4.2       Calculate variables that define the Chop Type Cases.
                ! ----------------------------------------------------------------------+-------
                DDU_a(I,JLEV,JDIR) = UDOTK(I,ILAUNCH(i),JDIR)            &
                     &                                - UDOTK(I,JLEV,JDIR)
                ! ----------------------------------------------------------------------+-------
                !             UMDP 34: 1.23 coefficient B
                ! ----------------------------------------------------------------------+-------
                ATTENUATION(I,JLEV,JDIR) =                              &
                     &           ( RHOCSK(I,JLEV) / RHOCL(I,JDIR) ) *               &
                     &           ( NBV(I,ILAUNCH(i)) / NBV(I,JLEV) )**TTM1
                ! ----------------------------------------------------------------------+-------
                !             UMDP 34: 1.22 coefficient A = (A/B) * B
                ! ----------------------------------------------------------------------+-------
                CURVATURE(I,JLEV,JDIR)  =  DDU_a(I,JLEV,JDIR) * MSTAR &
                     &                                     / NBV(I,ILAUNCH(i))
                !
                ACOEFF(I,JLEV,JDIR)     =  CURVATURE(I,JLEV,JDIR) *   &
                     &                                   ATTENUATION(I,JLEV,JDIR)
                !
                IF(TTM2 == 1.0_r8)THEN
                 MINTERCEPT(I,JLEV,JDIR) = ATTENUATION(I,JLEV,JDIR) * ( (MNLMIN * CURVATURE(I,JLEV,JDIR)) + 1.0_r8 )
                ELSE
                 MINTERCEPT(I,JLEV,JDIR) = ATTENUATION(I,JLEV,JDIR) * ( (MNLMIN * CURVATURE(I,JLEV,JDIR)) + 1.0_r8 )**TTM2
              END IF
              !
                Curv_if1: IF (CURVATURE(I,JLEV,JDIR) <  ASECN)  THEN
                   ! ----------------------------------------------------------------------+-------
                   !             Negative Doppler Shift : factor will hit zero (kill point)
                   ! ----------------------------------------------------------------------+-------
                   MKILL = 1.0_r8 / ABS(CURVATURE(I,JLEV,JDIR))
                   !
                   Mkill_if1: IF (MKILL <= MNLMIN )  THEN
                      ! ----------------------------------------------------------------------+-------
                      !               Chop Type IV : No flux propagates
                      ! ----------------------------------------------------------------------+-------
                      FPTOT(I,JLEV,JDIR) = 0.0_r8
                   ELSE
                      IF (MKILL >  1.0_r8)  THEN
                         IF(TTM2 == 1.0_r8)THEN
                            INTERCEPT1(I,JLEV,JDIR) = ATTENUATION(I,JLEV,JDIR)&
                              &              * ( 1.0_r8 + CURVATURE(I,JLEV,JDIR) )
                      ELSE
                            INTERCEPT1(I,JLEV,JDIR) = ATTENUATION(I,JLEV,JDIR)&
                              &              * ( 1.0_r8 + CURVATURE(I,JLEV,JDIR) )**TTM2
                      END IF    
                      ELSE
                         ! ----------------------------------------------------------------------+-------
                         !                 Doppler factor minimum (kill point) situated below
                         !                 mstar in the low-m part of the launch spectrum
                         ! ----------------------------------------------------------------------+-------
                         INTERCEPT1(I,JLEV,JDIR) = 0.0_r8
                      END IF
                      !
                      Lowend_if1: IF (INTERCEPT1(I,JLEV,JDIR) >= 1.0_r8) THEN
                         ! ----------------------------------------------------------------------+-------
                         !                 Chop Type I: Intersection in high wavenumber part only
                         ! ----------------------------------------------------------------------+-------
                         MNLY = ( ATTENUATION(I,JLEV,JDIR)**TTRAT -        &
                              &              ATTENUATION(I,JLEV,JDIR)) / ACOEFF(I,JLEV,JDIR)
                         !
                         FPTOT(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) *       &
                              &     (1.0_r8 - (A0_R_1MT * CURVATURE(I,JLEV,JDIR) * MNLY**TWOMTT))
                      ELSE
                         IF (MINTERCEPT(I,JLEV,JDIR) <= FMINUS)  THEN
                            ! ----------------------------------------------------------------------+-------
                            !                 Chop Type IIb: Low wavenumber intersect only below min
                            ! ----------------------------------------------------------------------+-------
                            FPTOT(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) *     &
                                 &            A0_R_SP1TM1 * TAIL_CHOP2B * MINTERCEPT(I,JLEV,JDIR) &
                                 &                 * ( (MNLMIN * CURVATURE(I,JLEV,JDIR)) + 1.0_r8 )
                         ELSE
                            ! ----------------------------------------------------------------------+-------
                            !                 Chop Type IIa: Low wavenumber intersect only
                            ! ----------------------------------------------------------------------+-------
                            L_CHOP2 = .TRUE.
                            MGUESS_a(I,JLEV,JDIR) = MIN(MKILL, 1.0_r8)
                            FPFAC(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) * A0_R_SP1TM1
                            FPTOT(I,JLEV,JDIR) = 0.0_r8
                         END IF
                      END IF  Lowend_if1
                      !
                   END IF  Mkill_if1
                   !
                ELSE IF (CURVATURE(I,JLEV,JDIR) >  ASECP)  THEN
                   ! ----------------------------------------------------------------------+-------
                   !             Positive Doppler Shift : non-zero factor (no kill point)
                   ! ----------------------------------------------------------------------+-------
                   IF(TTM2 == 1.0_r8)THEN
                      INTERCEPT1(I,JLEV,JDIR) = ATTENUATION(I,JLEV,JDIR)  &
                           &            * ( 1.0_r8 + CURVATURE(I,JLEV,JDIR) )
                   ELSE
                      INTERCEPT1(I,JLEV,JDIR) = ATTENUATION(I,JLEV,JDIR)  &
                           &            * ( 1.0_r8 + CURVATURE(I,JLEV,JDIR) )**TTM2                 
                 END IF
                 !
                   Chop3_if1: IF (INTERCEPT1(I,JLEV,JDIR) <  1.0_r8)  THEN
                      ! ----------------------------------------------------------------------+-------
                      !               Chop Type III: Intersection in both wavenumber parts
                      ! ----------------------------------------------------------------------+-------
                      FPFAC(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) * A0_R_SP1TM1
                      !
                      ! ----------------------------------------------------------------------+-------
                      !                 First find intersect in high wavenumber part
                      !                 UMDP 34: 1.25
                      ! ----------------------------------------------------------------------+-------
                      MNLY = ( ATTENUATION(I,JLEV,JDIR)**TTRAT -          &
                           &              ATTENUATION(I,JLEV,JDIR)) / ACOEFF(I,JLEV,JDIR)
                      !
                      FPTOT(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) *         &
                           &              A0_R_1MT * CURVATURE(I,JLEV,JDIR) * MNLY**TWOMTT
                      !
                      ! ----------------------------------------------------------------------+-------
                      !                 Then find intersect in low wavenumber part to reckon
                      !                 its flux contribution for addition when available
                      ! ----------------------------------------------------------------------+-------
                      IF (MINTERCEPT(I,JLEV,JDIR) <= FMINUS)  THEN
                         ! ----------------------------------------------------------------------+-------
                         !                 Chop Type IIIb: Low wavenumber intersect below min
                         ! ----------------------------------------------------------------------+-------
                         FPTOT(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) +       &
                              &                                   ( FPFAC(I,JLEV,JDIR) *       &
                              &                 TAIL_CHOP2B *  MINTERCEPT(I,JLEV,JDIR) *       &
                              &                ( (MNLMIN * CURVATURE(I,JLEV,JDIR)) + 1.0_r8 ) )
                      ELSE
                         ! ----------------------------------------------------------------------+-------
                         !                 Chop Type IIIa: Low wavenumber intersect
                         ! ----------------------------------------------------------------------+-------
                         L_CHOP2 = .TRUE.
                         MGUESS_a(I,JLEV,JDIR) = 1.0_r8
                      END IF
                      !
                      ! ----------------------------------------------------------------------+-------
                      !               ELSE Chop Type 0: No intersection (spectrum unaltered)
                      ! ----------------------------------------------------------------------+-------
                   END IF  Chop3_if1
                ELSE
                   ! ----------------------------------------------------------------------+-------
                   !             Negligible Doppler shift
                   ! ----------------------------------------------------------------------+-------
                   !               Strictly this is analytic solution mNLX.  UMDP 34: 1.27
                   MNLY = ATTENUATION(I,JLEV,JDIR)**RSSPTT
                   IF (MNLY <= MNLMIN)  THEN
                      ! ----------------------------------------------------------------------+-------
                      !               Chop Type IIb: Low wavenumber intersect only below min
                      ! ----------------------------------------------------------------------+-------
                      FPTOT(I,JLEV,JDIR) = FPTOT(I,JLEV,JDIR) *         &
                           &           A0_R_SP1TM1 * TAIL_CHOP2B * ATTENUATION(I,JLEV,JDIR)
                   ELSE
                      IF (MNLY <  1.0_r8)  FPTOT(I,JLEV,JDIR) =              &
                           ! ----------------------------------------------------------------------+-------
                           !               Chop Type IIc: Low wavenumber intersect only (analytic)
                           ! ----------------------------------------------------------------------+-------
                      &              FPTOT(I,JLEV,JDIR) * A0_R_SP1TM1 *                &
                           &               ( (SSPTT * (MNLY**SSP1)) - HEAD_CHOP2A )
                      ! ----------------------------------------------------------------------+-------
                      !                 ELSE Chop Type 0: No intersection (spectrum unaltered)
                      ! ----------------------------------------------------------------------+-------
                   END IF
                END IF  Curv_if1
                !
             END IF  Fptot_if1
             END IF
          END DO  nCols_do92
       END DO  IDir_do2a
       !
       Lchop2_if1: IF (L_CHOP2)  THEN
          ! ----------------------------------------------------------------------+-------
          !       Process low wavenumber contribution: evaluate intersect mNX
          ! ----------------------------------------------------------------------+-------
          NCHOP2 = 0
          !
          IDir_do2b: DO JDIR=1,IDIR
             nCols_do93: DO I=1,nCols
                IF (MGUESS_a(I,JLEV,JDIR) >  0.0_r8)  THEN
                   NCHOP2 = NCHOP2 + 1
                   !
                   INDEXJD(NCHOP2) = JDIR
                   !                  INDEXJ(NCHOP2)  = J
                   INDEXI(NCHOP2)  = I
                END IF
             END DO  nCols_do93
          END DO  IDir_do2b
          !         NVIEW(JLEV) = NCHOP2
          !         NSPIRA = 0
          !
          Nchop2_do1: DO I=1,NCHOP2
             ! ----------------------------------------------------------------------+-------
             !               Chop Type IIa : / Full solution required for mNlX
             !          or   Chop Type IIIa: ! ----------------------------------------------------------------------+-------
             NNJD = INT(INDEXJD(I))
             !            NNJ  = INDEXJ(I)
             NNI  = INT(INDEXI(I))
             !
             FPLUS  = MGUESS_a(NNI,JLEV,NNJD)**SSPTT
             GPLUS  = INTERCEPT1(NNI,JLEV,NNJD)
             !           FMINUS = MNLMIN**SSPTT    Defined as a constant
             GMINUS = MINTERCEPT(NNI,JLEV,NNJD)
             ATTE_C(I) = ATTENUATION(NNI,JLEV,NNJD)
             CURV_C(I) = CURVATURE(NNI,JLEV,NNJD)
             !
             FTERM = ( ((FMINUS / ATTE_C(I))**RTTM2) - 1.0_r8 ) / CURV_C(I)
             GTERM = GMINUS**RSSPTT
             L_FTHENG(I) = .FALSE.
             !
             Curv_if2: IF (CURVATURE(NNI,JLEV,NNJD) >  ASECP)  THEN
                ! ----------------------------------------------------------------------+-------
                !           Positive Doppler Shift
                ! ----------------------------------------------------------------------+-------
                WGTN(I) = 0.0_r8
             ELSE
                ! ----------------------------------------------------------------------+-------
                !           Negative Doppler Shift
                ! ----------------------------------------------------------------------+-------
                WGTN(I) = 1.0_r8 - MWEIGHT
                !
                IF (FPLUS <= GMINUS  .AND.  GPLUS >  FMINUS)  THEN
                   FTERM = (((FPLUS / ATTE_C(I))**RTTM2) - 1.0_r8)/ CURV_C(I)
                   GTERM = GPLUS**RSSPTT
                   L_FTHENG(I) = (GTERM  <   FTERM)
                   !
                ELSE IF (FPLUS >  GMINUS  .AND.  GPLUS <= FMINUS)  THEN
                   L_FTHENG(I) = (GTERM >= FTERM)
                   !
                ELSE IF (FPLUS <= GMINUS  .AND.  GPLUS <= FMINUS)  THEN
                   L_FTHENG(I) = .TRUE.
                   !
                   !             ELSE Use default settings
                END IF
             END IF  Curv_if2
             !
             IF (L_FTHENG(I))  THEN
                !             NSPIRA = NSPIRA + 1
                MNLX(I,0) = FTERM
             ELSE
                MNLX(I,0) = GTERM
             END IF
          END DO  Nchop2_do1
          !         NVIEW2(JLEV) = NSPIRA
          !
          Jwhile_do2: DO JWHILE=0,MAXWHILE-1
             Nchop2_do2: DO I=1,NCHOP2
                !
                IF (L_FTHENG(I))  THEN
                   ! ----------------------------------------------------------------------+-------
                   !           Obtain m_n+1 from g_n+1  = f_n (m_n)
                   ! ----------------------------------------------------------------------+-------
                   MNLX(I,JWHILE+1) = (                                    &
                        &           (((MNLX(I,JWHILE)**SSPTT) / ATTE_C(I))**RTTM2) - 1.0_r8 ) &
                        &           / CURV_C(I)
                ELSE
                   ! ----------------------------------------------------------------------+-------
                   !           Obtain m_n+1 from f_n+1  = g_n (m_n)
                   ! ----------------------------------------------------------------------+-------
                   IF(TTM2 == 1.0_r8)THEN
                       MNLX(I,JWHILE+1) = ( (ATTE_C(I) *                       &
                        &          ((1.0_r8 + (CURV_C(I) * MNLX(I,JWHILE)))))**RSSPTT )
                 ELSE
                       MNLX(I,JWHILE+1) = ( (ATTE_C(I) *                       &
                        &          ((1.0_r8 + (CURV_C(I) * MNLX(I,JWHILE)))**TTM2))**RSSPTT )
                 
                 END IF
                END IF
                !
                MNLX(I,JWHILE+1) = ((1.0_r8 - WGTN(I)) * MNLX(I,JWHILE+1)) + &
                     &                                  (WGTN(I)  * MNLX(I,JWHILE))
                !
             END DO  Nchop2_do2
          END DO  Jwhile_do2
          !
          !CDIR NODEP
          Nchop2_do3: DO I=1,NCHOP2
             NNJD = INT(INDEXJD(I))
             !            NNJ  = INDEXJ(I)
             NNI  = INT(INDEXI(I))
             !
             FPTOT(NNI,JLEV,NNJD) = FPTOT(NNI,JLEV,NNJD) +       &
                  &     (FPFAC(NNI,JLEV,NNJD) * ( ((MNLX(I,MAXWHILE)**SSP1) *    &
                  &      ( SSPTT + (SSP1 * MNLX(I,MAXWHILE) * CURV_C(I)) )) - HEAD_CHOP2A ))
             !
          END DO  Nchop2_do3
       END IF  Lchop2_if1
       !
       IDir_do2c: DO JDIR=1,IDIR
          nCols_do10: DO I=1,nCols
             !-----------------------------------------------------------------------+-------
             !         Now correct pseudomomentum flux in the evolved spectrum if the
             !         new value is non-physical (pseudomomentum flux cannot increase
             !         with altitude)
             !-----------------------------------------------------------------------+-------
             IF(JLEV==1)THEN
                FPTOT(I,JLEV,JDIR) =  MIN(FPTOT(I,JLEV,JDIR), FPTOT(I,JLEV,JDIR))
             ELSE
                FPTOT(I,JLEV,JDIR) =  MIN(FPTOT(I,JLEV,JDIR), FPTOT(I,JLEV-1,JDIR))
             END IF
             !
          END DO  nCols_do10
       END DO  IDir_do2c
       !
    END DO  Levels_do92

    !
    ! ----------------------------------------------------------------------+-------
    ! 4.5   If choosing Opaque Upper Boundary set fluxes to zero at top
    ! ----------------------------------------------------------------------+-------
    IF (L_USSP_OPAQUE) THEN
       IDir_do3: DO JDIR=1,IDIR
          nCols_do12: DO I=1,nCols
             FPTOT(I,kMax,JDIR) =  0.0_r8
          END DO  nCols_do12
       END DO  IDir_do3
    ENDIF
    !
    ! ----------------------------------------------------------------------+-------
    ! 5.0   Compute vertical divergence of pseudomomentum flux.
    ! ----------------------------------------------------------------------+-------
    IDir_do4: DO JDIR=1,IDIR
       Levels_do14: DO JLEV=1,kMax
          nCols_do14: DO I=1,nCols
             !           Pseudomomentum flux
             !FPTOT = kg m^-1 s^-2
             !G_G=    m s^-2)
             IF(JLEV>=ILAUNCH(i)+1.and. JLEV<=kMax)THEN
             G_G(I,JLEV,JDIR) =                                      &
                  &         g * (FPTOT(I,JLEV,JDIR) - FPTOT(I,JLEV-1,JDIR)) /    &
                  &          (P_LAYER_BOUNDARIES(I,JLEV) - P_LAYER_BOUNDARIES(I,JLEV-1))
             G_X(I,JLEV) = G_X(I,JLEV) + G_G(I,JLEV,JDIR) * COSPHI(JDIR)
             G_Y(I,JLEV) = G_Y(I,JLEV) + G_G(I,JLEV,JDIR) * SINPHI(JDIR)
             END IF
          END DO  nCols_do14
       END DO  Levels_do14
    END DO  IDir_do4

  END SUBROUTINE gw_ussp

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
END MODULE GwddSchemeUSSP
