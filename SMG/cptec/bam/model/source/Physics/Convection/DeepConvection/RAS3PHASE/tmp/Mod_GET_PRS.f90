MODULE Mod_GET_PRS
    IMPLICIT NONE
  SAVE

  PRIVATE
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(13,60) ! the '60' maps to 64-bit real
  REAL(kind=r8),PARAMETER:: con_cp     =1.0046e+3      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_rd     =2.8705e+2      ! gas constant air    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_rv     =4.6150e+2      ! gas constant H2O    (J/kg/K)

  REAL(kind=r8),PARAMETER:: con_fvirt  =con_rv/con_rd-1.
  REAL(kind=r8),PARAMETER:: con_rocp   =con_rd/con_cp


  REAL(kind=r8),PARAMETER   :: ri (0:20)=(/con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,&
       con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,con_rd,&
       con_rd/)
  REAL(kind=r8),PARAMETER   :: cpi(0:20)=(/con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,&
       con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,con_cp,&
       con_cp/)

  INTEGER :: thermodyn_id=1
  LOGICAL :: gen_coord_hybrid=.FALSE.
  PUBLIC  :: GET_PRS,GET_PHI2
CONTAINS
  !______________________________________________________________________________________________________________

  SUBROUTINE GET_PRS(&!
       ix     ,&!
       levs   ,&!
       ntrac  ,&!
       t      ,&!
       q      ,&!
       prsi   ,&!
       prki   ,&!
       prsl   ,&!
       prkl   ,&!
       phii   ,&!
       phil   ,&!
       del     )
    !
    !      USE tracer_const
    IMPLICIT NONE
    !
    INTEGER      , INTENT(IN   ) :: ix
    INTEGER      , INTENT(IN   ) :: levs
    INTEGER      , INTENT(IN   ) :: ntrac
    REAL(kind=r8), INTENT(IN   ) :: T(ix,levs)
    REAL(kind=r8), INTENT(IN   ) :: q(ix,levs,ntrac)
    REAL(kind=r8), INTENT(IN   ) :: prsi(ix,levs+1)
    REAL(kind=r8), INTENT(INOUT) :: prki(ix,levs+1)
    REAL(kind=r8), INTENT(INOUT) :: prsl(ix,levs)
    REAL(kind=r8), INTENT(INOUT) :: prkl(ix,levs)
    REAL(kind=r8), INTENT(OUT  ) :: phii(ix,levs+1)
    REAL(kind=r8), INTENT(OUT  ) :: phil(ix,levs)
    REAL(kind=r8), INTENT(OUT  ) :: del(ix,levs)

    REAL(kind=r8) :: xcp(ix,levs)
    REAL(kind=r8) :: xr(ix,levs)
    REAL(kind=r8) :: kappa(ix,levs)
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: dphib
    REAL(kind=r8) :: dphit
    REAL(kind=r8) :: dphi
    REAL (kind=r8), PARAMETER :: zero=0.0
    REAL (kind=r8), PARAMETER :: rkapi=1.0/con_rocp, rkapp1=1.0+con_rocp
    INTEGER       :: i, k
    !
    DO k=1,levs
       DO i=1,ix
          del(i,k) = PRSI(i,k) - PRSI(i,k+1)
       ENDDO
    ENDDO
    !
    IF( gen_coord_hybrid ) THEN                                       ! hmhj
       IF( thermodyn_id.EQ.3 ) THEN      ! Enthalpy case
          !
          ! hmhj : This is for generalized hybrid (Henry) with finite difference
          !        in the vertical and enthalpy as the prognostic (thermodynamic)
          !        variable.  However, the input "t" here is the temperature,
          !        not enthalpy (because this subroutine is called by gbphys where
          !        only temperature is available).
          !
          IF (prki(1,1) <= zero .OR. prkl(1,1) <= zero) THEN
             CALL GET_CPR(ix,levs,ntrac,q,xcp,xr)
             !
             DO k=1,levs
                DO i=1,ix
                   kappa(i,k) = xr(i,k)/xcp(i,k)
                   prsl(i,k)  = (PRSI(i,k) + PRSI(i,k+1))*0.5
                   prkl(i,k)  = (prsl(i,k)*0.01) ** kappa(i,k)
                ENDDO
             ENDDO
             DO k=2,levs
                DO i=1,ix
                   tem = 0.5 * (kappa(i,k) + kappa(i,k-1))
                   prki(i,k-1) = (prsi(i,k)*0.01) ** tem
                ENDDO
             ENDDO
             DO i=1,ix
                prki(i,1) = (prsi(i,1)*0.01) ** kappa(i,1)
             ENDDO
             k = levs + 1
             IF (prsi(1,k) .GT. 0.0) THEN
                DO i=1,ix
                   prki(i,k) = (prsi(i,k)*0.01) ** kappa(i,levs)
                ENDDO
             ENDIF
             !
             DO i=1,ix
                phii(i,1)   = 0.0           ! Ignoring topography height here
             ENDDO
             DO k=1,levs
                DO i=1,ix
                   TEM         = xr(i,k) * T(i,k)
                   DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM &
                        / (PRSI(i,k) + PRSI(i,k+1))
                   phil(i,k)   = phii(i,k) + DPHI
                   phii(i,k+1) = phil(i,k) + DPHI
                ENDDO
             ENDDO
          ENDIF
          IF (prsl(1,1) <= 0.0) THEN
             DO k=1,levs
                DO i=1,ix
                   prsl(i,k)  = (PRSI(i,k) + PRSI(i,k+1))*0.5
                ENDDO
             ENDDO
          ENDIF
          IF (phil(1,levs) <= 0.0) THEN ! If geopotential is not given, calculate
             DO i=1,ix
                phii(i,1)   = 0.0           ! Ignoring topography height here
             ENDDO
             CALL GET_R(ix,levs,ntrac,q,xr)
             DO k=1,levs
                DO i=1,ix
                   TEM         = xr(i,k) * T(i,k)
                   DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM &
                        / (PRSI(i,k) + PRSI(i,k+1))
                   phil(i,k)   = phii(i,k) + DPHI
                   phii(i,k+1) = phil(i,k) + DPHI
                ENDDO
             ENDDO
          ENDIF
       ELSE                                 ! gc Virtual Temp case
          IF (prki(1,1) <= zero .OR. prkl(1,1) <= zero) THEN
             DO k=1,levs
                DO i=1,ix
                   prsl(i,k) = (PRSI(i,k) + PRSI(i,k+1))*0.5
                   prkl(i,k) = (prsl(i,k)*0.01) ** con_rocp
                ENDDO
             ENDDO
             DO k=1,levs+1
                DO i=1,ix
                   prki(i,k) = (prsi(i,k)*0.01) ** con_rocp
                ENDDO
             ENDDO
             DO i=1,ix
                phii(i,1)   = 0.0           ! Ignoring topography height here
             ENDDO
             DO k=1,levs
                DO i=1,ix
                   TEM         = con_rd * T(i,k)*(1.0+con_fvirt*MAX(Q(i,k,1),zero))
                   DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM  &
                        / (PRSI(i,k) + PRSI(i,k+1))
                   phil(i,k)   = phii(i,k) + DPHI
                   phii(i,k+1) = phil(i,k) + DPHI
                ENDDO
             ENDDO
          ENDIF
          IF (prsl(1,1) <= 0.0) THEN
             DO k=1,levs
                DO i=1,ix
                   prsl(i,k)  = (PRSI(i,k) + PRSI(i,k+1))*0.5
                ENDDO
             ENDDO
          ENDIF
          IF (phil(1,levs) <= 0.0) THEN ! If geopotential is not given, calculate
             DO i=1,ix
                phii(i,1)   = 0.0         ! Ignoring topography height here
             ENDDO
             DO k=1,levs
                DO i=1,ix
                   TEM         = con_rd * T(i,k)*(1.0+con_fvirt*MAX(Q(i,k,1),zero))
                   DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM &
                        / (PRSI(i,k) + PRSI(i,k+1))
                   phil(i,k)   = phii(i,k) + DPHI
                   phii(i,k+1) = phil(i,k) + DPHI
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ELSE                                   ! Not gc Virtual Temp (Orig Joe)
       IF (prki(1,1) <= zero) THEN
          !                                      Pressure is in centibars!!!!
          DO i=1,ix
             prki(i,1) = (prsi(i,1)*0.01) ** con_rocp
          ENDDO
          DO k=1,levs
             DO i=1,ix
                prki(i,k+1) = (prsi(i,k+1)*0.01) ** con_rocp
                tem         = rkapp1 * del(i,k)
                prkl(i,k)   = (prki(i,k)*PRSI(i,k)-prki(i,k+1)*PRSI(i,k+1)) &
                     / tem
             ENDDO
          ENDDO

       ELSEIF (prkl(1,1) <= zero) THEN
          DO k=1,levs
             DO i=1,ix
                tem         = rkapp1 * del(i,k)
                prkl(i,k)   = (prki(i,k)*PRSI(i,k)-prki(i,k+1)*PRSI(i,k+1)) &
                     / tem
             ENDDO
          ENDDO
       ENDIF
       IF (prsl(1,1) <= 0.0) THEN
          DO k=1,levs
             DO i=1,ix
                PRSL(i,k)   = 100.0 * PRKL(i,k) ** rkapi
             ENDDO
          ENDDO
       ENDIF
       IF (phil(1,levs) <= 0.0) THEN ! If geopotential is not given, calculate
          DO i=1,ix
             phii(i,1)   = 0.0         ! Ignoring topography height here
          ENDDO
          DO k=1,levs
             DO i=1,ix
                TEM         = con_cp * T(i,k) * (1.0 + con_fvirt*MAX(Q(i,k,1),zero)) &
                     / PRKL(i,k)
                DPHIB       = (PRKI(i,k) - PRKL(i,k)) * TEM
                DPHIT       = (PRKL(i,k  ) - PRKI(i,k+1)) * TEM
                phil(i,k)   = phii(i,k) + DPHIB
                phii(i,k+1) = phil(i,k) + DPHIT
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    !
    RETURN
  END SUBROUTINE GET_PRS
  !______________________________________________________________________________________________________________

  SUBROUTINE GET_PHI(&
                    ix     , &!INTEGER	   , INTENT(IN   ) :: ix
                    levs   , &!INTEGER	   , INTENT(IN   ) :: levs
                    ntrac  , &!INTEGER	   , INTENT(IN   ) :: ntrac
                    t      , &!REAL(kind=r8), INTENT(IN   ) :: T(ix,levs)
                    q      , &!REAL(kind=r8), INTENT(IN   ) :: q(ix,levs,ntrac)
                    prsi   , &!REAL(kind=r8), INTENT(IN   ) :: prsi(ix,levs+1)
                    prki   , &!REAL(kind=r8), INTENT(IN   ) :: prki(ix,levs+1)
                    prkl   , &!REAL(kind=r8), INTENT(IN   ) :: prkl(ix,levs)
                    phii   , &!REAL(kind=r8), INTENT(IN   ) :: phii(ix,levs+1)
                    phil   , &!REAL(kind=r8), INTENT(IN   ) :: phil(ix,levs)
                    del      )
    ! 
    !      USE tracer_const
    IMPLICIT NONE
    !
    INTEGER      , INTENT(IN   ) :: ix
    INTEGER      , INTENT(IN   ) :: levs
    INTEGER      , INTENT(IN   ) :: ntrac
    REAL(kind=r8), INTENT(IN   ) :: T(ix,levs)
    REAL(kind=r8), INTENT(IN   ) :: q(ix,levs,ntrac)
    REAL(kind=r8), INTENT(IN   ) :: prsi(ix,levs+1)
    REAL(kind=r8), INTENT(IN   ) :: prki(ix,levs+1)
    REAL(kind=r8), INTENT(IN   ) :: prkl(ix,levs)
    REAL(kind=r8), INTENT(OUT  ) :: phii(ix,levs+1)
    REAL(kind=r8), INTENT(OUT  ) :: phil(ix,levs)
    REAL(kind=r8), INTENT(OUT  ) :: del (ix,levs)

    REAL(kind=r8) :: xr(ix,levs)
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: dphib
    REAL(kind=r8) :: dphit
    REAL(kind=r8) :: dphi
    REAL (kind=r8), PARAMETER :: zero=0.0
    INTEGER       :: i, k
    DO k=1,levs
       DO i=1,ix
          del(i,k) = PRSI(i,k) - PRSI(i,k+1)
       ENDDO
    ENDDO

    !
    DO i=1,ix
       phii(i,1)   = zero                     ! Ignoring topography height here
    ENDDO
    IF( gen_coord_hybrid ) THEN              ! hmhj
       IF( thermodyn_id.EQ.3 ) THEN           ! Enthalpy case
          CALL GET_R(ix,levs,ntrac,q,xr)
          DO k=1,levs
             DO i=1,ix
                TEM         = xr(i,k) * T(i,k)
                DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM &
                     /(PRSI(i,k) + PRSI(i,k+1))
                phil(i,k)   = phii(i,k) + DPHI
                phii(i,k+1) = phil(i,k) + DPHI
             ENDDO
          ENDDO
          !
       ELSE                                 ! gc Virtual Temp
          DO k=1,levs
             DO i=1,ix
                TEM         = con_rd * T(i,k) * (1.0 + con_fvirt*MAX(Q(i,k,1),zero))
                DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM &
                     /(PRSI(i,k) + PRSI(i,k+1))
                phil(i,k)   = phii(i,k) + DPHI
                phii(i,k+1) = phil(i,k) + DPHI
             ENDDO
          ENDDO
       ENDIF
    ELSE                                   ! Not gc Virt Temp (Orig Joe)
       DO k=1,levs
          DO i=1,ix
             TEM         = con_cp * T(i,k) * (1.0 + con_fvirt*MAX(Q(i,k,1),zero)) &
                  / PRKL(i,k)
             DPHIB       = (PRKI(i,k) - PRKL(i,k)) * TEM
             DPHIT       = (PRKL(i,k  ) - PRKI(i,k+1)) * TEM
             phil(i,k)   = phii(i,k) + DPHIB
             phii(i,k+1) = phil(i,k) + DPHIT
          ENDDO
       ENDDO
    ENDIF
    !
    RETURN
  END SUBROUTINE GET_PHI
  !______________________________________________________________________________________________________________

  SUBROUTINE GET_CPR(ix,levs,ntrac,q,xcp,xr)
    !
    !      USE tracer_const
    IMPLICIT NONE
    !
    REAL (kind=r8), PARAMETER :: zero=0.0
    INTEGER  ix, levs, ntrac
    REAL(kind=r8) q(ix,levs,ntrac)
    REAL(kind=r8) xcp(ix,levs),xr(ix,levs),sumq(ix,levs)
    INTEGER i, k, n
    !
    sumq = zero
    xr   = zero
    xcp  = zero
    DO n=1,ntrac
       IF( ri(n) > 0.0 ) THEN
          DO k=1,levs
             DO i=1,ix
                xr(i,k)   = xr(i,k)   + q(i,k,n) * ri(n)
                xcp(i,k)  = xcp(i,k)  + q(i,k,n) * cpi(n)
                sumq(i,k) = sumq(i,k) + q(i,k,n)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    DO k=1,levs
       DO i=1,ix
          xr(i,k)    = (1.-sumq(i,k))*ri(0)  + xr(i,k)
          xcp(i,k)   = (1.-sumq(i,k))*cpi(0) + xcp(i,k)
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE GET_CPR
  !______________________________________________________________________________________________________________

  SUBROUTINE GET_R(ix,levs,ntrac,q,xr)
    !
    !      USE tracer_const
    IMPLICIT NONE
    !
    REAL (kind=r8), PARAMETER :: zero=0.0
    INTEGER  ix, levs, ntrac
    REAL(kind=r8) q(ix,levs,ntrac)
    REAL(kind=r8) xr(ix,levs),sumq(ix,levs)
    INTEGER i, k, n
    !
    sumq = zero
    xr   = zero
    DO n=1,ntrac
       IF( ri(n) > 0.0 ) THEN
          DO k=1,levs
             DO i=1,ix
                xr(i,k)   = xr(i,k)   + q(i,k,n) * ri(n)
                sumq(i,k) = sumq(i,k) + q(i,k,n)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    DO k=1,levs
       DO i=1,ix
          xr(i,k)    = (1.-sumq(i,k))*ri(0)  + xr(i,k)
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE GET_R
  !______________________________________________________________________________________________________________
  SUBROUTINE GET_CP(ix,levs,ntrac,q,xcp)
    !
    !      USE tracer_const
    IMPLICIT NONE
    !
    REAL (kind=r8), PARAMETER :: zero=0.0
    INTEGER  ix, levs, ntrac
    REAL(kind=r8) q(ix,levs,ntrac)
    REAL(kind=r8) xcp(ix,levs),sumq(ix,levs)
    INTEGER i, k, n
    !
    sumq = zero
    xcp  = zero
    DO n=1,ntrac
       IF( cpi(n) > 0.0 ) THEN
          DO k=1,levs
             DO i=1,ix
                xcp(i,k)  = xcp(i,k)  + q(i,k,n) * cpi(n)
                sumq(i,k) = sumq(i,k) + q(i,k,n)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    DO k=1,levs
       DO i=1,ix
          xcp(i,k)   = (1.-sumq(i,k))*cpi(0) + xcp(i,k)
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE GET_CP
END MODULE Mod_GET_PRS
