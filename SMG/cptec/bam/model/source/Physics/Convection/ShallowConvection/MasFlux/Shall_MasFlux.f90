MODULE Shall_MasFlux
  USE Constants, ONLY :  i8
  
  USE Mod_GET_PRS, ONLY: GET_PHI

  IMPLICIT NONE
SAVE

  PRIVATE
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(13,60) ! the '60' maps to 64-bit real
  !  --- ...  Geophysics/Astronomy constants

  REAL(kind=r8),PARAMETER:: con_g      =9.80665e+0_r8     ! gravity           (m/s2)
  REAL(kind=r8),PARAMETER:: con_pi     =3.1415926535897931 ! pi
  REAL(kind=r8),PARAMETER:: con_rerth  =6.3712e+6      ! radius of earth   (m)
  !  --- ...  Thermodynamics constants


  REAL(kind=r8),PARAMETER:: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_rv     =4.6150e+2_r8      ! gas constant H2O    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_hvap   =2.5000e+6_r8      ! lat heat H2O cond   (J/kg)
  REAL(kind=r8),PARAMETER:: con_hfus   =3.3358e+5_r8      ! lat heat H2O fusion (J/kg)
  REAL(kind=r8),PARAMETER:: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cvap   =1.8460e+3_r8      ! spec heat H2O gas   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cliq   =4.1855e+3_r8      ! spec heat H2O liq   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_csol   =2.1060e+3_r8      ! spec heat H2O ice   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_ttp    =2.7316e+2_r8      ! temp at H2O 3pt     (K)
  REAL(kind=r8),PARAMETER:: con_psat   =6.1078e+2_r8      ! pres at H2O 3pt     (Pa)  
  real(kind=r8),PARAMETER:: con_t0c    =2.7315e+2_r8      ! temp at 0C          (K)

  !  Secondary constants

  REAL(kind=r8),PARAMETER:: con_rocp   =con_rd/con_cp
  REAL(kind=r8),PARAMETER:: con_fvirt  =con_rv/con_rd-1.0_r8 !(J/kg/K)/(J/kg/K)
  REAL(kind=r8),PARAMETER:: con_eps    =con_rd/con_rv
  REAL(kind=r8),PARAMETER:: con_epsm1  =con_rd/con_rv-1.0_r8
  REAL(kind=r8),PARAMETER :: PT01=0.01_r8 !1/100  cb
  REAL(kind=r8),PARAMETER :: rkap = con_rocp
  REAL(kind=r8),PARAMETER :: rk   = con_rocp

  !
  !   module funcphys
  INTEGER,PARAMETER:: nxpvs=7501

  REAL(r8) c1xpvs,c2xpvs,tbpvs(nxpvs)

  !   END module funcphys
   INTEGER  :: jcap


   PUBLIC :: Init_Shall_MasFlux
   PUBLIC :: Run_Shall_MasFlux
CONTAINS

  SUBROUTINE Init_Shall_MasFlux (kMax,trunc)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) ::  kMax
    INTEGER      , INTENT(IN   ) ::  trunc
    
    jcap=trunc
   
    CALL gfuncphys()

  END SUBROUTINE Init_Shall_MasFlux

  SUBROUTINE Run_Shall_MasFlux(nCols,kMax,dt,prsi_i ,prsl_i ,phii_i,phil_i ,tgrs,qgrs,ugrs,vgrs,omgb,qliq,qice,&
                               dudt,dvdt,dtdt,dqdt,dqldt,dqidt,kbot,ktop,kuo,noshal,mask,hpbl,sens,latheat)
    INTEGER      , INTENT(IN   ) :: nCols
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: dt

    REAL(KIND=r8), INTENT(in   ) :: prsi_i (1:nCols,1:kMax+1)  
    REAL(KIND=r8), INTENT(in   ) :: prsl_i (1:nCols,1:kMax  )  
    REAL(KIND=r8), INTENT(in   ) :: phii_i (1:nCols,1:kMax+1)  
    REAL(KIND=r8), INTENT(in   ) :: phil_i (1:nCols,1:kMax  )  
    REAL(KIND=r8), INTENT(INOUT) :: tgrs (nCols,kMax)!     tgrs     - real, layer mean temperature ( k )             ix,levs !
    REAL(KIND=r8), INTENT(INOUT) :: qgrs (nCols,kMax)!     qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
    REAL(KIND=r8), INTENT(IN   ) :: ugrs (nCols,kMax)!     ugrs,vgrs- real, u/v component of layer wind              ix,levs !
    REAL(KIND=r8), INTENT(IN   ) :: vgrs (nCols,kMax)
    REAL(KIND=r8), INTENT(IN   ) :: omgb (nCols,kMax) ! (Pa/s) !REAL(KIND=r8)   ,    INTENT(IN   ) :: omgb      (iMax,kMax) ! (Pa/s)
    REAL(KIND=r8), INTENT(INOUT) :: qliq (nCols,kMax)!     qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
    REAL(KIND=r8), INTENT(INOUT) :: qice (nCols,kMax)!     qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
    REAL(KIND=r8), INTENT(inout) :: dudt (nCols,kMax)
    REAL(KIND=r8), INTENT(inout) :: dvdt (nCols,kMax)
    INTEGER      , INTENT(INOUT) :: kbot (nCols)
    INTEGER      , INTENT(INOUT) :: ktop (nCols)
    INTEGER      , INTENT(IN   ) :: kuo(nCols)    !   kuo                     ! convection yes(1) or not(0) for shallow convection
    INTEGER      , INTENT(INOUT) :: noshal (nCols) 
    INTEGER(KIND=i8), INTENT(IN   ) :: mask    (1:nCols) 
    REAL(kind=r8), INTENT(IN   ) :: hpbl (nCols)  ! hpbl [m]
    REAL(KIND=r8), INTENT(in   ) :: sens    (1:nCols)
    REAL(KIND=r8), INTENT(in   ) :: latheat (1:nCols)
    REAL(KINd=r8), INTENT(OUT  ) :: dtdt (nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dqdt (nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dqldt(nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: dqidt(nCols,kMax)

    REAL(KIND=r8) :: ps   (nCols)        !cb
    REAL(kind=r8) :: heat (nCols)    !     hflx     - real, sensible heat flux  			im   !
    REAL(kind=r8) :: evap (nCols)    !     evap     - real, evaperation from latent heat flux		im   !
    REAL(kind=r8) :: rcs   (nCols)
    REAL(KIND=r8) :: prsi  (nCols,kMax+1) !     prsi     - real, pressure at layer interfaces             ix,levs+1
    REAL(KIND=r8) :: prsl  (nCols,kMax)   !     prsl     - real, mean layer presure                       ix,levs !
    REAL(KIND=r8) :: prsik (nCols,kMax+1)
    REAL(KIND=r8) :: prslk (nCols,kMax)
    REAL(KIND=r8) :: phii  (nCols,kMax+1)
    REAL(KIND=r8) :: phil  (nCols,kMax)
    REAL(KIND=r8) :: del   (nCols,kMax)
    REAL(KIND=r8) :: ql    (nCols,kMax,2)! 1-ice    ! 2-water
    REAL(kind=r8) :: dot   (nCols,kMax)!     vvel     - real, layer mean vertical velocity             ix,levs !
    REAL(kind=r8) :: rho   (nCols,kMax)
    REAL(kind=r8) :: rn    (nCols)
    REAL(kind=r8) :: q1    (nCols,kMax)
    REAL(kind=r8) :: t1    (nCols,kMax)
    REAL(kind=r8) :: u1    (nCols,kMax)
    REAL(kind=r8) :: v1    (nCols,kMax)
    REAL(kind=r8) :: slimsk(nCols)!           slmsk    - real, sea/land/ice mask (=0/1/2)                  im   !
    INTEGER       :: ncloud   !     ncld     - integer, number of cloud species                  1    !
    REAL(kind=r8) :: ud_mf (nCols,kMax)
    REAL(kind=r8) :: dt_mf (nCols,kMax)
    REAL(kind=r8) :: shalldudt (nCols,kMax)
    REAL(kind=r8) :: shalldvdt (nCols,kMax)
    REAL(kind=r8) ::  RKAPI 
    REAL(kind=r8) ::  RKAPP1
    REAL(KIND=r8) :: sik  (nCols,kMax+1)
    REAL(KIND=r8) :: sikp1(nCols,kMax+1)
    REAL(KIND=r8) :: DeltaP(nCols,kMax)
    REAL(KIND=r8) :: slk(nCols,kMax)
    REAL(kind=r8) :: pgrk(nCols)
    REAL(KIND=r8) :: delt,tem
    INTEGER       :: ntrac,i,k
    ntrac=0
    ncloud=1
    rcs=1.0_r8
    delt=2.0_r8*dt
    ! grell mask
    DO i=1,nCols      
       ps(i)=prsi_i(i,1)/1000.0_r8
       rcs  (i) = 1.0_r8       ! wind conversion
       heat (i) = sens       (i)   !     hflx     - real, sensible heat fluxim   !
       evap (i) = latheat    (i)   !     evap     - real, evaperation from latent heat fluxim   !
       IF(mask(i).GT.0_i8)THEN
          !mask2(i)=0 ! land
          slimsk(i) = 1.0_r8
       ELSE
          slimsk(i) = 0.0_r8
          !mask2(i)=1 ! seaice/water/ocean
       END IF
    END DO
    !tottracer=ntrac
    DO i=1, nCols
       prsi(i,kMax+1) =prsi_i(i,kMax+1)/1000.0_r8
       phii(i,kMax+1) =phii_i(i,kMax+1)*con_g
    END DO
    DO k=1,kMax
       DO i=1, nCols
           prsi (i,k)=prsi_i(i,k)/1000.0_r8
           prsl (i,k)=prsl_i(i,k)/1000.0_r8
           phii (i,k)=phii_i(i,k)*con_g
           phil (i,k)=phil_i(i,k)*con_g
           del(i,k) = (PRSI_i(i,k) - PRSI_i(i,k+1))/1000.0_r8  !cb
       END DO      
    END DO

    DO k=1,kMax
       DO i=1,nCols
         dtdt (i,k)= 0.0_r8
         dqdt (i,k)= 0.0_r8
         dqldt(i,k)= 0.0_r8
         dqidt(i,k)= 0.0_r8
         shalldudt(i,k)= 0.0_r8
         shalldvdt(i,k)= 0.0_r8
         ql(i,k,1) = qice (i,k)
         ql(i,k,2) = qliq (i,k)
         q1(i,k)   = qgrs (i,k) 
         t1(i,k)   = tgrs (i,k)
         u1(i,k)   = ugrs (i,k)
         v1(i,k)   = vgrs (i,k)
       END DO
    END DO
    DO k=1,kMax
      DO i=1,nCols
          DeltaP(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO

    RKAPI  = 1.0_r8 / RKAP
    RKAPP1 = 1.0_r8 + RKAP
    DO k=1,kMax+1
       DO i=1,nCols 
          sik(i,k)   = (prsi(i,k)/prsi(i,1)) ** rkap
          sikp1(i,k) = (prsi(i,k)/prsi(i,1)) ** rkapp1
       END DO
    END DO
    DO k=1,kMax
       DO i=1,nCols 
          tem        = rkapp1 * DeltaP(i,k)
          slk(i,k)   = (sikp1(i,k)-sikp1(i,k+1))/tem
       END DO
    END DO

    DO i=1,nCols
         pgrk(i)         = (prsi(i,1)*pt01) ** rk
         prsik(i,kMax+1) = sik(i,kMax+1) * pgrk(i)
    END DO
    DO k=1,kMax
      DO i=1,nCols
          prsik(i,k) = sik(i,k) * pgrk(i)
          prslk(i,k) = slk(i,k) * pgrk(i)
      END DO
    END DO
  
    DO k=1,kMax
       DO i=1,nCols
         !; p_at_units = "Pa" 
         !; t_at_units = "K" 
         !   RGAS = 287. ; J/(kg-K) => m2/(s2 K) 
         !; omega_at_units = "Pa/sec" 
         !   rho = p/(RGAS*t) ; density => kg/m3 
         rho(i,k) = (1000.0_r8*prsl(i,k))/(con_rd*tgrs (i,k))
         !   GRAV = 9.8 ; m/s2 

         !   w = -omega/(rho*GRAV) 

         dot(i,k) = -omgb (i,k)/(rho(i,k)*  con_g ) 
       END DO
    END DO
    !    prepare input, erase output
    !
    DO i=1,nCols
       noshal(i)=0
       kbot(i)=1
       ktop(i)=1
    END DO

    CALL shalcnv(                     &
       nCols                        , &! INTEGER      , INTENT(IN   ) :: im
       kMax                         , &! INTEGER      , INTENT(IN   ) :: km
       jcap                         , &! INTEGER      , INTENT(IN   ) :: jcap
       delt                         , &! REAL(kind=r8), INTENT(IN   ) :: delt
       del    (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(IN   ) :: del   (im,km)
       prsl   (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(IN   ) :: prsl  (im,km)
       ps     (1:nCols)             , &! REAL(kind=r8), INTENT(IN   ) :: ps    (im)
       phil   (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(IN   ) :: phil  (im,km)
       ql     (1:nCols,1:kMax,1:2)  , &! REAL(kind=r8), INTENT(INOUT) :: ql    (im,km,2)
       q1     (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(INOUT) :: q1    (im,km)
       t1     (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(INOUT) :: t1    (im,km)
       u1     (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(INOUT) :: u1    (im,km)
       v1     (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(INOUT) :: v1    (im,km)
       rcs    (1:nCols)             , &! REAL(kind=r8), INTENT(IN   ) :: rcs   (im)
       rn     (1:nCols)             , &! REAL(kind=r8), INTENT(OUT  ) :: rn    (im)
       kbot   (1:nCols)             , &! INTEGER      , INTENT(INOUT) :: kbot  (im)
       ktop   (1:nCols)             , &! INTEGER      , INTENT(INOUT) :: ktop  (im)
       kuo    (1:nCols)             , &! INTEGER      , INTENT(in) :: ktop  (im)
       noshal   (1:nCols)             , &! INTEGER      , INTENT(INOUT) :: noshal  (im) 
       slimsk (1:nCols)             , &! REAL(kind=r8), INTENT(IN   ) :: slimsk(im)
       dot    (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(IN   ) :: dot   (im,km)
       ncloud                       , &! INTEGER      , INTENT(IN   ) :: ncloud 
       hpbl   (1:nCols)             , &! REAL(kind=r8), INTENT(IN   ) :: hpbl  (im)
       heat   (1:nCols)             , &! REAL(kind=r8), INTENT(IN   ) :: heat  (im)
       evap   (1:nCols)             , &! REAL(kind=r8), INTENT(IN   ) :: evap  (im)
       ud_mf  (1:nCols,1:kMax)      , &! REAL(kind=r8), INTENT(OUT  ) :: ud_mf (im,km)
       dt_mf  (1:nCols,1:kMax)        )! REAL(kind=r8), INTENT(OUT  ) :: dt_mf (im,km)

    DO k=1,kMax
       DO i=1,nCols
          dtdt (i,k)= (t1(i,k) -tgrs (i,k))/delt
          tgrs (i,k) = t1(i,k)

          dqdt (i,k)= (q1(i,k) -qgrs (i,k))/delt
          qgrs (i,k) = q1(i,k)

          dqldt(i,k)= (ql(i,k,2) - qliq (i,k))/delt
          qliq (i,k) = ql(i,k,2) 

          dqidt(i,k)= (ql(i,k,1) -qice (i,k))/delt
          qice (i,k) = ql(i,k,1) 

          shalldudt (i,k) = (u1(i,k) - ugrs (i,k) )/delt
          shalldvdt (i,k) = (v1(i,k) - vgrs (i,k) )/delt
          dudt (i,k) = dudt (i,k) + shalldudt (i,k)
          dvdt (i,k) = dvdt (i,k) + shalldvdt (i,k)
       END DO
    END DO


  END SUBROUTINE Run_Shall_MasFlux

  SUBROUTINE shalcnv( &
       im        , &! INTEGER      , INTENT(IN   ) :: im
       km        , &! INTEGER      , INTENT(IN   ) :: km
       jcap      , &! INTEGER      , INTENT(IN   ) :: jcap
       delt      , &! REAL(kind=r8), INTENT(IN   ) :: delt
       del       , &! REAL(kind=r8), INTENT(IN   ) :: del   (im,km)
       prsl      , &! REAL(kind=r8), INTENT(IN   ) :: prsl  (im,km)
       ps        , &! REAL(kind=r8), INTENT(IN   ) :: ps    (im)
       phil      , &! REAL(kind=r8), INTENT(IN   ) :: phil  (im,km)
       ql        , &! REAL(kind=r8), INTENT(INOUT) :: ql    (im,km,2)
       q1        , &! REAL(kind=r8), INTENT(INOUT) :: q1    (im,km)
       t1        , &! REAL(kind=r8), INTENT(INOUT) :: t1    (im,km)
       u1        , &! REAL(kind=r8), INTENT(INOUT) :: u1    (im,km)
       v1        , &! REAL(kind=r8), INTENT(INOUT) :: v1    (im,km)
       rcs       , &! REAL(kind=r8), INTENT(IN   ) :: rcs   (im)
       rn        , &! REAL(kind=r8), INTENT(OUT  ) :: rn    (im)
       kbot      , &! INTEGER      , INTENT(INOUT) :: kbot  (im)
       ktop      , &! INTEGER      , INTENT(INOUT) :: ktop  (im)
       kuo       , &! INTEGER      , INTENT(INOUT) :: ktop  (im)
       noshal      , &! INTEGER      , INTENT(INOUT) :: noshal  (im) 
       slimsk    , &! REAL(kind=r8), INTENT(IN   ) :: slimsk(im)
       dot       , &! REAL(kind=r8), INTENT(IN   ) :: dot   (im,km)
       ncloud    , &! INTEGER      , INTENT(IN   ) :: ncloud 
       hpbl      , &! REAL(kind=r8), INTENT(IN   ) :: hpbl  (im)
       heat      , &! REAL(kind=r8), INTENT(IN   ) :: heat  (im)
       evap      , &! REAL(kind=r8), INTENT(IN   ) :: evap  (im)
       ud_mf     , &! REAL(kind=r8), INTENT(OUT  ) :: ud_mf (im,km)
       dt_mf       )! REAL(kind=r8), INTENT(OUT  ) :: dt_mf (im,km)
    !
    IMPLICIT NONE
    !
    INTEGER      , INTENT(IN   ) :: im
    INTEGER      , INTENT(IN   ) :: km
    INTEGER      , INTENT(IN   ) :: jcap
    REAL(kind=r8), INTENT(IN   ) :: delt
    REAL(kind=r8), INTENT(IN   ) :: del   (im,km)
    REAL(kind=r8), INTENT(IN   ) :: prsl  (im,km)
    REAL(kind=r8), INTENT(IN   ) :: ps    (im)
    REAL(kind=r8), INTENT(IN   ) :: phil  (im,km)
    REAL(kind=r8), INTENT(INOUT) :: ql    (im,km,2)
    REAL(kind=r8), INTENT(INOUT) :: q1    (im,km)
    REAL(kind=r8), INTENT(INOUT) :: t1    (im,km)
    REAL(kind=r8), INTENT(INOUT) :: u1    (im,km)
    REAL(kind=r8), INTENT(INOUT) :: v1    (im,km)
    REAL(kind=r8), INTENT(IN   ) :: rcs   (im)
    REAL(kind=r8), INTENT(OUT  ) :: rn    (im)
    INTEGER      , INTENT(INOUT) :: kbot(im)
    INTEGER      , INTENT(INOUT) :: ktop(im)
    INTEGER      , INTENT(in   ) :: kuo    (im)        !   kuo                     ! convection yes(1) or not(0) for shallow convection
    INTEGER      , INTENT(INOUT) :: noshal(im) 
    REAL(kind=r8), INTENT(IN   ) :: slimsk(im)
    REAL(kind=r8), INTENT(IN   ) :: dot   (im,km)!     vvel     - real, layer mean vertical velocity             im,levs !
    INTEGER      , INTENT(IN   ) :: ncloud !     ncld [1 - 0]    - integer, indicate cloud types                     1    !
    REAL(kind=r8), INTENT(IN   ) :: hpbl  (im)  ! hpbl [m]
    REAL(kind=r8), INTENT(IN   ) :: heat  (im)    !     hflx     - real, sensible heat flux                          im   !
    REAL(kind=r8), INTENT(IN   ) :: evap  (im)    !     evap     - real, evaperation from latent heat flux           im   !

    ! hchuang code change mass flux output
    REAL(kind=r8), INTENT(OUT  ) :: ud_mf (im,km)
    REAL(kind=r8), INTENT(OUT  ) :: dt_mf (im,km)
    !
    INTEGER       :: i
    INTEGER       :: indx
!    INTEGER       :: jmn
    INTEGER       :: k
    INTEGER       :: kk
!    INTEGER       :: latd
!    INTEGER       :: lond
    INTEGER       :: km1
    INTEGER       :: kpbl(im)
    !
!    REAL(kind=r8) :: alpha
!    REAL(kind=r8) :: alphal
!    REAL(kind=r8) :: alphas
    REAL(kind=r8) :: dellat
    REAL(kind=r8) :: desdt
!    REAL(kind=r8) :: deta
!    REAL(kind=r8) :: detad
!    REAL(kind=r8) :: dg
!    REAL(kind=r8) :: dh
!    REAL(kind=r8) :: dhh
!    REAL(kind=r8) :: dlnsig
    REAL(kind=r8) :: dp
    REAL(kind=r8) :: dq
    REAL(kind=r8) :: dqsdp
    REAL(kind=r8) :: dqsdt
    REAL(kind=r8) :: dt
    REAL(kind=r8) :: dt2
    !REAL(kind=r8) :: dtmax
    !REAL(kind=r8) :: dtmin
    REAL(kind=r8) :: dv1h
    REAL(kind=r8) :: dv1q
    REAL(kind=r8) :: dv2h
    REAL(kind=r8) :: dv2q
    REAL(kind=r8) :: dv1u
    REAL(kind=r8) :: dv1v
    REAL(kind=r8) :: dv2u
    REAL(kind=r8) :: dv2v
    REAL(kind=r8) :: dv3q
    REAL(kind=r8) :: dv3h
    REAL(kind=r8) :: dv3u
    REAL(kind=r8) :: dv3v
    REAL(kind=r8) :: clam
    REAL(kind=r8) :: dz
    REAL(kind=r8) :: e1
    REAL(kind=r8) :: es
    REAL(kind=r8) :: etah
    REAL(kind=r8) :: evef
    REAL(kind=r8) :: evfact
    REAL(kind=r8) :: evfactl
    REAL(kind=r8) :: factor
    !REAL(kind=r8) :: fjcap
    REAL(kind=r8) :: gamma
    REAL(kind=r8) :: pprime
    REAL(kind=r8) :: betaw
    REAL(kind=r8) :: qc
    REAL(kind=r8) :: qlk
    REAL(kind=r8) :: qrch
    REAL(kind=r8) :: qs
!    REAL(kind=r8) ::  rain
!    REAL(kind=r8) :: rfact
    REAL(kind=r8) :: shear
    REAL(kind=r8) :: tem1
!    REAL(kind=r8) :: tem2
    REAL(kind=r8) :: val
    REAL(kind=r8) :: val1
    REAL(kind=r8) :: val2
    REAL(kind=r8) :: w1
    REAL(kind=r8) :: w1l
    REAL(kind=r8) :: w1s
    REAL(kind=r8) :: w2
    REAL(kind=r8) :: w2l
    REAL(kind=r8) :: w2s
    REAL(kind=r8) :: w3
    REAL(kind=r8) :: w3l
    REAL(kind=r8) :: w3s
    REAL(kind=r8) :: w4
    REAL(kind=r8) :: w4l
    REAL(kind=r8) :: w4s
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: ptem
    REAL(kind=r8) :: ptem1
    REAL(kind=r8) :: pgcon
    !
    INTEGER       :: kb(im)
    INTEGER       :: kbcon(im)
    INTEGER       :: kbcon1(im)
    INTEGER       :: ktcon(im)
    INTEGER       :: kbm(im)
    INTEGER       :: kmax(im)
    !
    REAL(kind=r8) :: delhbar(im)
    REAL(kind=r8) :: delq(im)
    REAL(kind=r8) :: delq2(im)

    REAL(kind=r8) :: delqbar(im)
    REAL(kind=r8) :: delqev(im)
    REAL(kind=r8) :: deltbar(im)

    REAL(kind=r8) :: deltv(im)
    REAL(kind=r8) :: edt(im)

    REAL(kind=r8) :: wstar(im)
    REAL(kind=r8) :: sflx(im)

    REAL(kind=r8) :: pdot(im)
    REAL(kind=r8) :: po(im,km)

    REAL(kind=r8) :: qcond(im)
    REAL(kind=r8) :: qevap(im)
    REAL(kind=r8) :: hmax(im)

    REAL(kind=r8) :: rntot(im)
    REAL(kind=r8) :: vshear(im)

!    REAL(kind=r8) :: xk(im)
    REAL(kind=r8) :: xlamud(im)

    REAL(kind=r8) :: xmb(im)
    REAL(kind=r8) :: xmbmax(im)

    REAL(kind=r8) :: delubar(im)
    REAL(kind=r8) :: delvbar(im)

    REAL(kind=r8) :: cincr
    !1
    !   physical parameters

    REAL(kind=r8), PARAMETER :: c0=0.002_r8
!    REAL(kind=r8), PARAMETER :: cpoel=con_cp/con_hvap
!    REAL(kind=r8), PARAMETER :: delta=con_fvirt
    REAL(kind=r8), PARAMETER :: c1=3.e-4_r8
    REAL(kind=r8), PARAMETER :: el2orc=con_hvap*con_hvap/(con_rv*con_cp)
    REAL(kind=r8), PARAMETER :: elocp=con_hvap/con_cp
    REAL(kind=r8), PARAMETER :: h1=0.33333333_r8
    REAL(kind=r8), PARAMETER :: dthk=25.0_r8
    REAL(kind=r8), PARAMETER :: g =con_g
    REAL(kind=r8), PARAMETER :: fact1 =(con_cvap-con_cliq)/con_rv
    REAL(kind=r8), PARAMETER :: fact2 =con_hvap/con_rv-fact1*con_t0c
!    REAL(kind=r8), PARAMETER :: terr=0.0_r8
    REAL(kind=r8), PARAMETER :: cincrmax=180.0_r8
    REAL(kind=r8), PARAMETER :: cincrmin=120.0_r8




    !  local variables and arrays
    REAL(kind=r8) :: pfld(im,km)
    REAL(kind=r8) :: to(im,km)
    REAL(kind=r8) :: qo(im,km)
    REAL(kind=r8) :: uo(im,km)
    REAL(kind=r8) :: vo(im,km)
    REAL(kind=r8) :: qeso(im,km)
    !  cloud water
    !     real(kind=r8) qlko_ktcon(im), dellal(im,km), tvo(im,km),
    REAL(kind=r8) :: qlko_ktcon(im)
    REAL(kind=r8) :: dellal(im,km)
    REAL(kind=r8) :: dbyo(im,km)
    REAL(kind=r8) :: zo(im,km)
    REAL(kind=r8) :: xlamue(im,km)
    REAL(kind=r8) :: heo(im,km)
    REAL(kind=r8) :: heso(im,km)
    REAL(kind=r8) :: dellah(im,km)
    REAL(kind=r8) :: dellaq(im,km)

    REAL(kind=r8) :: dellau(im,km)
    REAL(kind=r8) :: dellav(im,km)
    REAL(kind=r8) :: hcko(im,km)
    REAL(kind=r8) :: ucko(im,km)
    REAL(kind=r8) :: vcko(im,km)
    REAL(kind=r8) :: qcko(im,km)

    REAL(kind=r8) :: eta(im,km)
    REAL(kind=r8) :: zi(im,km)
    REAL(kind=r8) :: pwo(im,km)

    REAL(kind=r8) :: tx1(im)
    !
    LOGICAL       :: totflg
    LOGICAL       :: cnvflg(im)
    LOGICAL       :: flg(im)
    !
    REAL(kind=r8), PARAMETER :: tf =233.16_r8
    REAL(kind=r8), PARAMETER :: tcr=263.16_r8
    REAL(kind=r8), PARAMETER :: tcrf=1.0_r8/(tcr-tf)
    !      parameter (tf=233.16_r8, tcr=263.16_r8, tcrf=1.0_r8/(tcr-tf))
    !
    !-----------------------------------------------------------------------
    !
    km1 = km - 1
    !
    !  compute surface buoyancy flux
    !
    DO i=1,im
       !  con_fvirt  =con_rv/con_rd-1.0_r8 !(J/kg/K)/(J/kg/K)

       sflx(i) = heat(i)+ con_fvirt*t1(i,1)*evap(i)
    ENDDO
    !
    !  initialize arrays
    !
    DO i=1,im
       noshal(i)=0
       IF(kuo(i) .EQ. 1) THEN
          noshal(i)=1             !   1 do not  shallow convection
       END IF
    END DO

    DO i=1,im
       cnvflg(i) = .TRUE.
       IF(kuo(i).EQ.1) cnvflg(i) = .FALSE.  !   deep convection yes(1) or not(0) for shallow convection
       IF(sflx(i).LE.0.0_r8) cnvflg(i) = .FALSE.
       IF(cnvflg(i)) THEN
          kbot(i)=km+1
          ktop(i)=0
       ENDIF
       rn(i)=0.0_r8
       kbcon(i)=km
       ktcon(i)=1
       kb(i)=km
       pdot(i) = 0.0_r8
       qlko_ktcon(i) = 0.0_r8
       edt(i)  = 0.0_r8
       vshear(i) = 0.0_r8
    ENDDO
    ! hchuang code change
    DO k = 1, km
       DO i = 1, im
          ud_mf(i,k) = 0.0_r8
          dt_mf(i,k) = 0.0_r8
       ENDDO
    ENDDO
    !!
    totflg = .TRUE.
    DO i=1,im
       totflg = totflg .AND. (.NOT. cnvflg(i))
    ENDDO
    IF(totflg) RETURN
    !!
    !
    dt2   = delt
    !val   =         1200.0_r8
    !dtmin = MAX(dt2, val )
    !val   =         3600.0_r8
    !dtmax = MAX(dt2, val )
    !  model tunable parameters are all here
    !alphal  = 0.5_r8
    !alphas  = 0.5_r8
    clam    = 0.3_r8
    betaw   = 0.03_r8
    !     evef    = 0.07_r8
    evfact  = 0.3_r8
    evfactl = 0.3_r8
    !
    pgcon   = 0.7_r8     ! Gregory et al. (1997, QJRMS)
    !pgcon   = 0.55_r8    ! Zhang & Wu (2003,JAS)
    !pgcon   = 0.9_r8    ! kubota (2015,JAS)
    !fjcap   = (float(jcap) / 126.0_r8) ** 2
    val     =           1.0_r8
    !fjcap   = MAX(fjcap,val)
    w1l     = -8.e-3_r8 
    w2l     = -4.e-2_r8
    w3l     = -5.e-3_r8 
    w4l     = -5.e-4_r8
    w1s     = -2.e-4_r8
    w2s     = -2.e-3_r8
    w3s     = -1.e-3_r8
    w4s     = -2.e-5_r8
    !
    !  define top layer for search of the downdraft originating layer
    !  and the maximum thetae for updraft
    !
    DO i=1,im
       kbm(i)   = km
       kmax(i)  = km
       tx1(i)   = 1.0_r8 / ps(i)
    ENDDO
    !     
    DO k = 1, km
       DO i=1,im
          IF (prsl(i,k)*tx1(i) .GT. 0.70_r8) kbm(i)   = k + 1
          IF (prsl(i,k)*tx1(i) .GT. 0.60_r8) kmax(i)  = k + 1
       ENDDO
    ENDDO
    DO i=1,im
       kbm(i)   = MIN(kbm(i),kmax(i))
    ENDDO
    !
    !  hydrostatic height assume zero terr and compute
    !  updraft entrainment rate as an inverse function of height
    !
    DO k = 1, km
       DO i=1,im
          zo(i,k) = phil(i,k) / g
       ENDDO
    ENDDO
    DO k = 1, km1
       DO i=1,im
          zi(i,k) = 0.5_r8*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
       ENDDO
    ENDDO
    DO i=1,im
       xlamue(i,km) = xlamue(i,km1)
    ENDDO
    !
    !  pbl height
    !
    DO i=1,im
       flg(i) = cnvflg(i)
       kpbl(i)= 1
    ENDDO
    DO k = 2, km1
       DO i=1,im
          IF (flg(i).AND.zo(i,k).LE.hpbl(i)) THEN
             kpbl(i) = k
          ELSE
             flg(i) = .FALSE.
          ENDIF
       ENDDO
    ENDDO
    DO i=1,im
       kpbl(i)= MIN(kpbl(i),kbm(i))
    ENDDO
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   convert surface pressure to mb from cb
    !
    DO k = 1, km
       DO i = 1, im
          IF (cnvflg(i) .AND. k .LE. kmax(i)) THEN
             pfld(i,k) = prsl(i,k) * 10.0_r8
             eta(i,k)  = 1.0_r8
             hcko(i,k) = 0.0_r8
             qcko(i,k) = 0.0_r8
             ucko(i,k) = 0.0_r8
             vcko(i,k) = 0.0_r8
             dbyo(i,k) = 0.0_r8
             pwo(i,k)  = 0.0_r8
             dellal(i,k) = 0.0_r8
             to(i,k)   = t1(i,k)
             qo(i,k)   = q1(i,k)
             uo(i,k)   = u1(i,k) * rcs(i)
             vo(i,k)   = v1(i,k) * rcs(i)
          ENDIF
       ENDDO
    ENDDO
    !
    !  column variables
    !  p is pressure of the layer (mb)
    !  t is temperature at t-dt (k)..tn
    !  q is mixing ratio at t-dt (kg/kg)..qn
    !  to is temperature at t+dt (k)... this is after advection and turbulan
    !  qo is mixing ratio at t+dt (kg/kg)..q1
    !
    DO k = 1, km
       DO i=1,im
          IF (cnvflg(i) .AND. k .LE. kmax(i)) THEN
             qeso(i,k) = 0.01_r8 * fpvs(to(i,k))      ! fpvs is in pa
             qeso(i,k) = con_eps * qeso(i,k) / (pfld(i,k) + con_epsm1*qeso(i,k))
             val1      =             1.e-8_r8
             qeso(i,k) = MAX(qeso(i,k), val1)
             val2      =           1.e-10_r8
             qo(i,k)   = MAX(qo(i,k), val2 )
             !           qo(i,k)   = min(qo(i,k),qeso(i,k))
             !           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          ENDIF
       ENDDO
    ENDDO
    !
    !  compute moist static energy
    !
    DO k = 1, km
       DO i=1,im
          IF (cnvflg(i) .AND. k .LE. kmax(i)) THEN
             !           tem       = g * zo(i,k) + con_cp * to(i,k)
             tem       = phil(i,k) + con_cp * to(i,k)
             heo(i,k)  = tem  + con_hvap * qo(i,k)
             heso(i,k) = tem  + con_hvap * qeso(i,k)
             !           heo(i,k)  = min(heo(i,k),heso(i,k))
          ENDIF
       ENDDO
    ENDDO
    !
    !  determine level with largest moist static energy within pbl
    !  this is the level where updraft starts
    !
    DO i=1,im
       IF (cnvflg(i)) THEN
          hmax(i) = heo(i,1)
          kb(i) = 1
       ENDIF
    ENDDO
    DO k = 2, km
       DO i=1,im
          IF (cnvflg(i).AND.k.LE.kpbl(i)) THEN
             IF(heo(i,k).GT.hmax(i)) THEN
                kb(i)   = k
                hmax(i) = heo(i,k)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    DO k = 1, km1
       DO i=1,im
          IF (cnvflg(i) .AND. k .LE. kmax(i)-1) THEN
             dz      = 0.5_r8 * (zo(i,k+1) - zo(i,k))
             dp      = 0.5_r8 * (pfld(i,k+1) - pfld(i,k))
             es      = 0.01_r8 * fpvs(to(i,k+1))      ! fpvs is in pa
             pprime  = pfld(i,k+1) + con_epsm1 * es
             qs      = con_eps * es / pprime
             dqsdp   = - qs / pprime
             desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
             dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
             gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
             dt      = (g * dz + con_hvap * dqsdp * dp) / (con_cp * (1.0_r8 + gamma))
             dq      = dqsdt * dt + dqsdp * dp
             to(i,k) = to(i,k+1) + dt
             qo(i,k) = qo(i,k+1) + dq
             po(i,k) = 0.5_r8 * (pfld(i,k) + pfld(i,k+1))
          ENDIF
       ENDDO
    ENDDO
    !
    DO k = 1, km1
       DO i=1,im
          IF (cnvflg(i) .AND. k .LE. kmax(i)-1) THEN
             qeso(i,k) = 0.01_r8 * fpvs(to(i,k))      ! fpvs is in pa
             qeso(i,k) = con_eps * qeso(i,k) / (po(i,k) + con_epsm1*qeso(i,k))
             val1      =             1.e-8_r8
             qeso(i,k) = MAX(qeso(i,k), val1)
             val2      =           1.e-10_r8
             qo(i,k)   = MAX(qo(i,k), val2 )
             !           qo(i,k)   = min(qo(i,k),qeso(i,k))
             heo(i,k)  = 0.5_r8 * g * (zo(i,k) + zo(i,k+1)) +   &
                  con_cp * to(i,k) + con_hvap * qo(i,k)
             heso(i,k) = 0.5_r8 * g * (zo(i,k) + zo(i,k+1)) +   &
                  con_cp * to(i,k) + con_hvap * qeso(i,k)
             uo(i,k)   = 0.5_r8 * (uo(i,k) + uo(i,k+1))
             vo(i,k)   = 0.5_r8 * (vo(i,k) + vo(i,k+1))
          ENDIF
       ENDDO
    ENDDO
    !
    !  look for the level of free convection as cloud base
    !
    DO i=1,im
       flg(i)   = cnvflg(i)
       IF(flg(i)) kbcon(i) = kmax(i)
    ENDDO
    DO k = 1, km1
       DO i=1,im
          IF (flg(i).AND.k.LT.kbm(i)) THEN
             IF(k.GT.kb(i).AND.heo(i,kb(i)).GT.heso(i,k)) THEN
                kbcon(i) = k
                flg(i)   = .FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    DO i=1,im
       IF(cnvflg(i)) THEN
          IF(kbcon(i).EQ.kmax(i)) cnvflg(i) = .FALSE.
       ENDIF
    ENDDO
    !!
    totflg = .TRUE.
    DO i=1,im
       totflg = totflg .AND. (.NOT. cnvflg(i))
    ENDDO
    IF(totflg) RETURN
    !!
    !
    !  determine critical convective inhibition
    !  as a function of vertical velocity at cloud base.
    !
    DO i=1,im
       IF(cnvflg(i)) THEN
          pdot(i)  = 10.0_r8* dot(i,kbcon(i))
       ENDIF
    ENDDO
    DO i=1,im
       IF(cnvflg(i)) THEN
          IF(slimsk(i).EQ.1.0_r8) THEN
             w1 = w1l
             w2 = w2l
             w3 = w3l
             w4 = w4l
          ELSE
             w1 = w1s
             w2 = w2s
             w3 = w3s
             w4 = w4s
          ENDIF
          IF(pdot(i).LE.w4) THEN
             ptem = (pdot(i) - w4) / (w3 - w4)
          ELSEIF(pdot(i).GE.-w4) THEN
             ptem = - (pdot(i) + w4) / (w4 - w3)
          ELSE
             ptem = 0.0_r8
          ENDIF
          val1    =             -1.0_r8
          ptem = MAX(ptem,val1)
          val2    =             1.0_r8
          ptem = MIN(ptem,val2)
          ptem = 1.0_r8 - ptem
          ptem1= 0.5_r8*(cincrmax-cincrmin)
          cincr = cincrmax - ptem * ptem1
          tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))
          IF(tem1.GT.cincr) THEN
             cnvflg(i) = .FALSE.
          ENDIF
       ENDIF
    ENDDO
    !!
    totflg = .TRUE.
    DO i=1,im
       totflg = totflg .AND. (.NOT. cnvflg(i))
    ENDDO
    IF(totflg) RETURN
    !!
    !
    !  assume the detrainment rate for the updrafts to be same as 
    !  the entrainment rate at cloud base
    !
    DO i = 1, im
       IF(cnvflg(i)) THEN
          xlamud(i) = xlamue(i,kbcon(i))
       ENDIF
    ENDDO
    !
    !  determine updraft mass flux for the subcloud layers
    !
    DO k = km1, 1, -1
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.LT.kbcon(i).AND.k.GE.kb(i)) THEN
                dz       = zi(i,k+1) - zi(i,k)
                ptem     = 0.5_r8*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
                eta(i,k) = eta(i,k+1) / (1.0_r8 + ptem * dz)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !  compute updraft cloud property
    !
    DO i = 1, im
       IF(cnvflg(i)) THEN
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
       ENDIF
    ENDDO
    !
    DO k = 2, km1
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.GT.kb(i).AND.k.LT.kmax(i)) THEN
                dz   = zi(i,k) - zi(i,k-1)
                tem  = 0.5_r8 * (xlamue(i,k)+xlamue(i,k-1)) * dz
                tem1 = 0.5_r8 * xlamud(i) * dz
                factor = 1.0_r8 + tem - tem1
                ptem = 0.5_r8 * tem + pgcon
                ptem1= 0.5_r8 * tem - pgcon
                hcko(i,k) = ((1.0_r8-tem1)*hcko(i,k-1)+tem*0.5_r8*   &
                     (heo(i,k)+heo(i,k-1)))/factor
                ucko(i,k) = ((1.0_r8-tem1)*ucko(i,k-1)+ptem*uo(i,k)   &
                     +ptem1*uo(i,k-1))/factor
                vcko(i,k) = ((1.0_r8-tem1)*vcko(i,k-1)+ptem*vo(i,k)   &
                     +ptem1*vo(i,k-1))/factor
                dbyo(i,k) = hcko(i,k) - heso(i,k)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !   taking account into convection inhibition due to existence of
    !    dry layers below cloud base
    !
    DO i=1,im
       flg(i) = cnvflg(i)
       kbcon1(i) = kmax(i)
    ENDDO
    DO k = 2, km1
       DO i=1,im
          IF (flg(i).AND.k.LT.kbm(i)) THEN
             IF(k.GE.kbcon(i).AND.dbyo(i,k).GT.0.0_r8) THEN
                kbcon1(i) = k
                flg(i)    = .FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i=1,im
       IF(cnvflg(i)) THEN
          IF(kbcon1(i).EQ.kmax(i)) cnvflg(i) = .FALSE.
       ENDIF
    ENDDO
    DO i=1,im
       IF(cnvflg(i)) THEN
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          IF(tem.GT.dthk) THEN
             cnvflg(i) = .FALSE.
          ENDIF
       ENDIF
    ENDDO
    !!
    totflg = .TRUE.
    DO i = 1, im
       totflg = totflg .AND. (.NOT. cnvflg(i))
    ENDDO
    IF(totflg) RETURN
    !!
    !
    !  determine convective cloud top as the level of zero buoyancy
    !   but limited to the level of 700 mb
    !
    DO i = 1, im
       flg(i) = cnvflg(i)
       IF(flg(i)) ktcon(i) = kbm(i)
    ENDDO
    DO k = 2, km1
       DO i=1,im
          IF (flg(i).AND.k .LT. kbm(i)) THEN
             IF(k.GT.kbcon1(i).AND.dbyo(i,k).LT.0.0_r8) THEN
                ktcon(i) = k
                flg(i)   = .FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !  turn off shallow convection if cloud top is less than pbl top
    !
    DO i=1,im
       IF(cnvflg(i)) THEN
          kk = kpbl(i)+1
          IF(ktcon(i).LE.kk) cnvflg(i) = .FALSE.
       ENDIF
    ENDDO
    !!
    totflg = .TRUE.
    DO i = 1, im
       totflg = totflg .AND. (.NOT. cnvflg(i))
    ENDDO
    IF(totflg) RETURN
    !!
    !
    !  compute updraft mass flux for the cloud layers
    !
    DO k = 2, km1
       DO i = 1, im
          IF(cnvflg(i)) THEN
             IF(k.GT.kbcon(i).AND.k.LE.ktcon(i)) THEN
                dz       = zi(i,k) - zi(i,k-1)
                ptem     = 0.5_r8*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
                eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !  specify upper limit of mass flux at cloud base
    !
    DO i = 1, im
       IF(cnvflg(i)) THEN
          !         xmbmax(i) = .1_r8
          !
          k = kbcon(i)
          dp = 1000.0_r8 * del(i,k)
          xmbmax(i) = dp / (g * dt2)
          !
          !         tem = dp / (g * dt2)
          !         xmbmax(i) = min(tem, xmbmax(i))
       ENDIF
    ENDDO
    !
    !  compute cloud moisture property and precipitation
    !
    DO k = 2, km
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.GT.kb(i).AND.k.LT.ktcon(i)) THEN
                dz    = zi(i,k) - zi(i,k-1)
                gamma = el2orc * qeso(i,k) / (to(i,k)**2)
                qrch = qeso(i,k)    &
                     + gamma * dbyo(i,k) / (con_hvap * (1.0_r8 + gamma))
                !j
                tem  = 0.5_r8 * (xlamue(i,k)+xlamue(i,k-1)) * dz
                tem1 = 0.5_r8 * xlamud(i) * dz
                factor = 1.0_r8 + tem - tem1
                qcko(i,k) = ((1.0_r8-tem1)*qcko(i,k-1)+tem*0.5_r8*   &
                     (qo(i,k)+qo(i,k-1)))/factor
                !j
                dq = eta(i,k) * (qcko(i,k) - qrch)
                !
                !             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
                !
                !  below lfc check if there is excess moisture to release latent heat
                !
                IF(dq.GT.0.0_r8) THEN
                   dp = 1000.0_r8 * del(i,k)
                   etah = 0.5_r8 * (eta(i,k) + eta(i,k-1))
                   qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                   qc = qlk + qrch
                   pwo(i,k) = etah * c0 * dz * qlk
                   dellal(i,k) = etah * c1 * dz * qlk * g / dp
                   qcko(i,k)= qc
                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !  this section is ready for cloud water
    !
    IF(ncloud.GT.0) THEN
       !
       !  compute liquid and vapor separation at cloud top
       !
       DO i = 1, im
          IF(cnvflg(i)) THEN
             k = ktcon(i)
             gamma = el2orc * qeso(i,k) / (to(i,k)**2)
             qrch = qeso(i,k) &
                  + gamma * dbyo(i,k) / (con_hvap * (1.0_r8 + gamma))
             dq = qcko(i,k-1) - qrch
             !
             !  check if there is excess moisture to release latent heat
             !
             IF(dq.GT.0.0_r8) THEN
                qlko_ktcon(i) = dq
                qcko(i,k-1) = qrch
             ENDIF
          ENDIF
       ENDDO
    ENDIF
    !
    !--- compute precipitation efficiency in terms of windshear
    !
    DO i = 1, im
       IF(cnvflg(i)) THEN
          vshear(i) = 0.0_r8
       ENDIF
    ENDDO
    DO k = 2, km
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.GT.kb(i).AND.k.LE.ktcon(i)) THEN
                shear= SQRT((uo(i,k)-uo(i,k-1)) ** 2   &
                     + (vo(i,k)-vo(i,k-1)) ** 2)
                vshear(i) = vshear(i) + shear
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, im
       IF(cnvflg(i)) THEN
          vshear(i) = 1.e3_r8 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591_r8-0.639_r8*vshear(i)   &
               +0.0953_r8*(vshear(i)**2)-0.00496_r8*(vshear(i)**3)
          edt(i)=1.0_r8-e1
          val =         0.9_r8
          edt(i) = MIN(edt(i),val)
          val =         0.0_r8
          edt(i) = MAX(edt(i),val)
       ENDIF
    ENDDO
    !
    !--- what would the change be, that a cloud with unit mass
    !--- will do to the environment?
    !
    DO k = 1, km
       DO i = 1, im
          IF(cnvflg(i) .AND. k .LE. kmax(i)) THEN
             dellah(i,k) = 0.0_r8
             dellaq(i,k) = 0.0_r8
             dellau(i,k) = 0.0_r8
             dellav(i,k) = 0.0_r8
          ENDIF
       ENDDO
    ENDDO
    !
    !--- changed due to subsidence and entrainment
    !
    DO k = 2, km1
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.GT.kb(i).AND.k.LT.ktcon(i)) THEN
                dp = 1000.0_r8 * del(i,k)
                dz = zi(i,k) - zi(i,k-1)
                !
                dv1h = heo(i,k)
                dv2h = 0.5_r8 * (heo(i,k) + heo(i,k-1))
                dv3h = heo(i,k-1)
                dv1q = qo(i,k)
                dv2q = 0.5_r8 * (qo(i,k) + qo(i,k-1))
                dv3q = qo(i,k-1)
                dv1u = uo(i,k)
                dv2u = 0.5_r8 * (uo(i,k) + uo(i,k-1))
                dv3u = uo(i,k-1)
                dv1v = vo(i,k)
                dv2v = 0.5_r8 * (vo(i,k) + vo(i,k-1))
                dv3v = vo(i,k-1)
                !
                tem  = 0.5_r8 * (xlamue(i,k)+xlamue(i,k-1))
                tem1 = xlamud(i)
                !j
                dellah(i,k) = dellah(i,k) +                   &
                     &     ( eta(i,k)*dv1h - eta(i,k-1)*dv3h               &
                     &    -  tem*eta(i,k-1)*dv2h*dz                        &
                     &    +  tem1*eta(i,k-1)*0.5_r8*(hcko(i,k)+hcko(i,k-1))*dz &
                     &         ) *g/dp
                !j
                dellaq(i,k) = dellaq(i,k) +                   &
                     &     ( eta(i,k)*dv1q - eta(i,k-1)*dv3q               &
                     &    -  tem*eta(i,k-1)*dv2q*dz                        &
                     &    +  tem1*eta(i,k-1)*0.5_r8*(qcko(i,k)+qcko(i,k-1))*dz &
                     &         ) *g/dp
                !j
                dellau(i,k) = dellau(i,k) +                   &
                     &     ( eta(i,k)*dv1u - eta(i,k-1)*dv3u               &
                     &    -  tem*eta(i,k-1)*dv2u*dz                        &
                     &    +  tem1*eta(i,k-1)*0.5_r8*(ucko(i,k)+ucko(i,k-1))*dz &
                     &    -  pgcon*eta(i,k-1)*(dv1u-dv3u)                  &
                     &         ) *g/dp
                !j
                dellav(i,k) = dellav(i,k) +                   &
                     &     ( eta(i,k)*dv1v - eta(i,k-1)*dv3v               &
                     &    -  tem*eta(i,k-1)*dv2v*dz                        &
                     &    +  tem1*eta(i,k-1)*0.5_r8*(vcko(i,k)+vcko(i,k-1))*dz &
                     &    -  pgcon*eta(i,k-1)*(dv1v-dv3v)                  &
                     &         ) *g/dp
                !j
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !------- cloud top
    !
    DO i = 1, im
       IF(cnvflg(i)) THEN
          indx = ktcon(i)
          dp = 1000.0_r8 * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *                     &
               &                     (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *                     &
               &                     (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *                     &
               &                     (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *                     &
               &                     (vcko(i,indx-1) - dv1v) * g / dp
          !
          !  cloud water
          !
          dellal(i,indx) = eta(i,indx-1) *                     &
               &                     qlko_ktcon(i) * g / dp
       ENDIF
    ENDDO
    !
    !  mass flux at cloud base for shallow convection
    !  (Grant, 2001)
    !
    DO i= 1, im
       IF(cnvflg(i)) THEN
          k = kbcon(i)
          !         ptem = g*sflx(i)*zi(i,k)/t1(i,1)
          ptem = g*sflx(i)*hpbl(i)/t1(i,1)
          wstar(i) = ptem**h1
          tem = po(i,k)*100.0_r8 / (con_rd*t1(i,k))
          xmb(i) = betaw*tem*wstar(i)
          xmb(i) = MIN(xmb(i),xmbmax(i))
       ENDIF
    ENDDO
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    DO k = 1, km
       DO i = 1, im
          IF (cnvflg(i) .AND. k .LE. kmax(i)) THEN
             qeso(i,k) = 0.01_r8 * fpvs(t1(i,k))      ! fpvs is in pa
             qeso(i,k) = con_eps * qeso(i,k) / (pfld(i,k) + con_epsm1*qeso(i,k))
             val     =             1.e-8_r8
             qeso(i,k) = MAX(qeso(i,k), val )
          ENDIF
       ENDDO
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    DO i = 1, im
       delhbar(i) = 0.0_r8
       delqbar(i) = 0.0_r8
       deltbar(i) = 0.0_r8
       delubar(i) = 0.0_r8
       delvbar(i) = 0.0_r8
       qcond(i) = 0.0_r8
    ENDDO
    DO k = 1, km
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.GT.kb(i).AND.k.LE.ktcon(i)) THEN
                dellat = (dellah(i,k) - con_hvap * dellaq(i,k)) / con_cp
                t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
                q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
                tem = 1.0_r8/rcs(i)
                u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
                v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
                dp = 1000.0_r8 * del(i,k)
                delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
                delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
                deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
                delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
                delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO k = 1, km
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.GT.kb(i).AND.k.LE.ktcon(i)) THEN
                qeso(i,k) = 0.01_r8 * fpvs(t1(i,k))      ! fpvs is in pa
                qeso(i,k) = con_eps * qeso(i,k)/(pfld(i,k) + con_epsm1*qeso(i,k))
                val     =             1.e-8_r8
                qeso(i,k) = MAX(qeso(i,k), val )
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    DO i = 1, im
       rntot(i) = 0.0_r8
       delqev(i) = 0.0_r8
       delq2(i) = 0.0_r8
       flg(i) = cnvflg(i)
    ENDDO
    DO k = km, 1, -1
       DO i = 1, im
          IF (cnvflg(i)) THEN
             IF(k.LT.ktcon(i).AND.k.GT.kb(i)) THEN
                rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * 0.001_r8 * dt2
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    ! evaporating rain
    !
    DO k = km, 1, -1
       DO i = 1, im
          IF (k .LE. kmax(i)) THEN
             deltv(i) = 0.0_r8
             delq(i) = 0.0_r8
             qevap(i) = 0.0_r8
             IF(cnvflg(i)) THEN
                IF(k.LT.ktcon(i).AND.k.GT.kb(i)) THEN
                   rn(i) = rn(i) + pwo(i,k) * xmb(i) * 0.001_r8 * dt2
                ENDIF
             ENDIF
             IF(flg(i).AND.k.LT.ktcon(i)) THEN
                evef = edt(i) * evfact
                IF(slimsk(i).EQ.1.0_r8) evef=edt(i) * evfactl
                !             if(slimsk(i).eq.1.0_r8) evef=.07_r8
                !             if(slimsk(i).ne.1.0_r8) evef = 0.0_r8
                qcond(i) = evef * (q1(i,k) - qeso(i,k))             &
                     / (1.0_r8 + el2orc * qeso(i,k) / t1(i,k)**2)
                dp = 1000.0_r8 * del(i,k)
                IF(rn(i).GT.0.0_r8.AND.qcond(i).LT.0.0_r8) THEN
                   qevap(i) = -qcond(i) * (1.0_r8-EXP(-0.32_r8*SQRT(dt2*rn(i))))
                   qevap(i) = MIN(qevap(i), rn(i)*1000.0_r8*con_g/dp)
                   delq2(i) = delqev(i) + 0.001_r8 * qevap(i) * dp / g
                ENDIF
                IF(rn(i).GT.0.0_r8.AND.qcond(i).LT.0.0_r8.AND.       &
                     delq2(i).GT.rntot(i)) THEN
                   qevap(i) = 1000.0_r8* g * (rntot(i) - delqev(i)) / dp
                   flg(i) = .FALSE.
                ENDIF
                IF(rn(i).GT.0.0_r8.AND.qevap(i).GT.0.0_r8) THEN
                   tem  = 0.001_r8 * dp / g
                   tem1 = qevap(i) * tem
                   IF(tem1.GT.rn(i)) THEN
                      qevap(i) = rn(i) / tem
                      rn(i) = 0.0_r8
                   ELSE
                      rn(i) = rn(i) - tem1
                   ENDIF
                   q1(i,k) = q1(i,k) + qevap(i)
                   t1(i,k) = t1(i,k) - elocp * qevap(i)
                   deltv(i) = - elocp*qevap(i)/dt2
                   delq(i) =  + qevap(i)/dt2
                   delqev(i) = delqev(i) + 0.001_r8*dp*qevap(i)/g
                ENDIF
                dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
                delqbar(i) = delqbar(i) + delq(i)*dp/g
                deltbar(i) = deltbar(i) + deltv(i)*dp/g
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !j
    !     do i = 1, im
    !     if(me.eq.31.and.cnvflg(i)) then
    !     if(cnvflg(i)) then
    !       print *, ' shallow delhbar, delqbar, deltbar = ',
    !    &             delhbar(i),con_hvap*delqbar(i),cp*deltbar(i)
    !       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
    !       print *, ' precip =', con_hvap*rn(i)*1000./dt2
    !       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
    !     endif
    !     enddo
    !j
    DO i = 1, im
       IF(cnvflg(i)) THEN
          IF(rn(i).LT.0.0_r8.OR..NOT.flg(i)) rn(i) = 0.0_r8
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          noshal(i) = 0
       ENDIF
    ENDDO
    !
    !  cloud water
    !
    IF (ncloud.GT.0) THEN
       !
       DO k = 1, km1
          DO i = 1, im
             IF (cnvflg(i)) THEN
                IF (k.GT.kb(i).AND.k.LE.ktcon(i)) THEN
                   tem  = dellal(i,k) * xmb(i) * dt2
                   tem1 = MAX(0.0_r8, MIN(1.0_r8, (tcr-t1(i,k))*tcrf))
                   IF (ql(i,k,2) .GT. -999.0_r8) THEN
                      ql(i,k,1) = ql(i,k,1) + tem * tem1              ! ice
                      ql(i,k,2) = ql(i,k,2) + tem *(1.0_r8-tem1)       ! water
                   ELSE
                      ql(i,k,1) = ql(i,k,1) + tem
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       !
    ENDIF
    !
    ! hchuang code change
    !
    DO k = 1, km
       DO i = 1, im
          IF(cnvflg(i)) THEN
             IF(k.GE.kb(i) .AND. k.LT.ktop(i)) THEN
                ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DO i = 1, im
       IF(cnvflg(i)) THEN
          k = ktop(i)-1
          dt_mf(i,k) = ud_mf(i,k)
       ENDIF
    ENDDO
    !!
    RETURN
  END  SUBROUTINE shalcnv
  
  
  !-------------------------------------------------------------------------------

      SUBROUTINE sig2press(nCols    ,& !
                           kMax     ,&!
                           pgr      ,&!
                           sl       ,&!
                           si       ,&!
                           prsi     ,&!
                           prsl     ,&!
                           prsik    ,&!
                           prslk      )
 
      IMPLICIT NONE
 
      INTEGER      , INTENT(IN   ) ::  nCols
      INTEGER      , INTENT(IN   ) ::  kMax
      REAL(kind=r8), INTENT(IN   ) ::  pgr(nCols)    !cb
      REAL(kind=r8), INTENT(IN   ) ::  sl(kMax)
      REAL(kind=r8), INTENT(IN   ) ::  si(kMax+1)
      REAL(kind=r8), INTENT(OUT  ) ::  prsi(nCols,kMax+1)
      REAL(kind=r8), INTENT(OUT  ) ::  prsl(nCols,kMax)
      REAL(kind=r8), INTENT(OUT  ) ::  prsik(nCols,kMax+1)
      REAL(kind=r8), INTENT(OUT  ) ::  prslk(nCols,kMax)
      
      REAL(kind=r8) :: pgrk(nCols)
      REAL(kind=r8) :: tem, rkapi, rkapp1
      REAL(kind=r8) :: slk(kMax) 
      REAL(kind=r8) :: sik(kMax+1)
      REAL(kind=r8) :: sikp1(kMax+1)
      INTEGER       :: i,k

      REAL(kind=r8), PARAMETER :: PT01=0.01
      REAL(kind=r8), PARAMETER:: con_rd     =2.8705e+2      ! gas constant air    (J/kg/K)
      REAL(kind=r8), PARAMETER:: con_cp     =1.0046e+3      ! spec heat air @p    (J/kg/K)
      REAL(kind=r8), PARAMETER:: con_rocp   =con_rd/con_cp
      REAL(kind=r8), PARAMETER:: rkap = con_rocp
      REAL(kind=r8), PARAMETER:: rk = con_rocp
!      REAL(kind=r8), PARAMETER :: cb2mb   = 10.0

      RKAPI  = 1.0 / RKAP
      RKAPP1 = 1.0 + RKAP
      DO k=1,kMax+1
        sik(k)   = si(k) ** rkap
        sikp1(k) = si(k) ** rkapp1
      END DO
      DO k=1,kMax
        tem      = rkapp1 * (si(k) - si(k+1))
        slk(k)   = (sikp1(k)-sikp1(k+1))/tem
        !sl(k)    = slk(k) ** rkapi
      END DO

      DO i=1,nCols
         prsi(i,kMax+1)  = si(kMax+1)*pgr(i)      ! prsi are now pressures
         pgrk(i)         = (pgr(i)*pt01) ** rk
         prsik(i,kMax+1) = sik(kMax+1) * pgrk(i)
      END DO
      DO k=1,kMax
        DO i=1,nCols
          prsi(i,k)  = si(k)*pgr(i)               ! prsi are now pressures
          prsl(i,k)  = sl(k)*pgr(i)
          prsik(i,k) = sik(k) * pgrk(i)
          prslk(i,k) = slk(k) * pgrk(i)
        END DO
      END DO
      RETURN
      END SUBROUTINE sig2press

  !-----------------------------------------------------------------------------------------  
  
    !-------------------------------------------------------------------------------
  SUBROUTINE gfuncphys()
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gfuncphys    Compute all physics function tables
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute all physics function tables.  Lookup tables are
    !   set up for computing saturation vapor pressure, dewpoint temperature,
    !   equivalent potential temperature, moist adiabatic temperature and humidity,
    !   pressure to the kappa, and lifting condensation level temperature.
    !
    ! Program History Log:
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gfuncphys
    !
    ! Subprograms called:
    !   gpvs        compute saturation vapor pressure table
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL gpvs()
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gfuncphys
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  SUBROUTINE gpvs()
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvs         Compute saturation vapor pressure table
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvs.
    !   Exact saturation vapor pressures are calculated in subprogram fpvsx.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gpvs
    !
    ! Subprograms called:
    !   (fpvsx)    inlinable function to compute saturation vapor pressure
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) :: xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvs-1)
    !   c1xpvs=1.-xmin/xinc
    c2xpvs=1.0_r8/xinc
    c1xpvs=1.0_r8-xmin*c2xpvs
    DO jx=1,nxpvs
       x=xmin+(jx-1)*xinc
       t=x
       tbpvs(jx)=fpvsx(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvs 
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvs(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvs         Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvs. See documentation for fpvsx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvs is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvs(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvs       Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) :: fpvs
    REAL(r8),INTENT(in):: t
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvs+c2xpvs*t,1.0_r8),REAL(nxpvs,r8))
    jx=INT(MIN(xj,nxpvs-1.0_r8))
    fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvs
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsx(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsx        Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The saturation vapor pressure over either liquid and ice is computed
    !   over liquid for temperatures above the triple point,
    !   over ice for temperatures 20 degress below the triple point,
    !   and a linear combination of the two for temperatures in between.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas, liquid and ice, and neglects the volume of the condensate.
    !   The model does account for the variation of the latent heat
    !   of condensation and sublimation with temperature.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   The reference for this computation is Emanuel(1994), pages 116-117.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvsx(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsx      Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) :: fpvsx
    REAL(r8),INTENT(in) :: t
    REAL(r8),PARAMETER :: tliq=con_ttp
    REAL(r8),PARAMETER :: tice=con_ttp-20.0
    REAL(r8),PARAMETER :: dldtl=con_cvap-con_cliq
    REAL(r8),PARAMETER :: heatl=con_hvap
    REAL(r8),PARAMETER :: xponal=-dldtl/con_rv
    REAL(r8),PARAMETER :: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
    REAL(r8),PARAMETER :: dldti=con_cvap-con_csol
    REAL(r8),PARAMETER :: heati=con_hvap+con_hfus
    REAL(r8),PARAMETER :: xponai=-dldti/con_rv
    REAL(r8),PARAMETER :: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
    REAL(r8) tr,w,pvl,pvi
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    IF(t.GE.tliq) THEN
       fpvsx=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
    ELSEIF(t.LT.tice) THEN
       fpvsx=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
    ELSE
       w=(t-tice)/(tliq-tice)
       pvl=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
       pvi=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
       fpvsx=w*pvl+(1.0_r8-w)*pvi
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION  fpvsx
  
  
  

END MODULE Shall_MasFlux


!PROGRAM Main
!  USE Shall_MasFlux
!END PROGRAM MAin
