MODULE Gwdd_ECMWF
  USE Utils, ONLY   : &
        IJtoIBJB     , &
        LinearIJtoIBJB

    IMPLICIT NONE
  SAVE

  PRIVATE
  !Kim, Y.-j. and Doyle, J. D. (2005), Extension of an orographic-drag parametrization 
  !                     scheme to incorporate orographic anisotropy and flow blocking. 
  !                 Q.J.R. Meteorol. Soc., 131: 1893-1921. doi: 10.1256/qj.04.160
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers
  REAL(kind=r8),PARAMETER :: con_g      =9.80665e+0_r8     ! gravity           (m/s2)
  REAL(kind=r8),PARAMETER :: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER :: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  REAL(kind=r8),PARAMETER :: con_rv     =4.6150e+2_r8      ! gas constant H2O    (J/kg/K)
  REAL(kind=r8),PARAMETER :: con_fvirt  =con_rv/con_rd-1.0_r8
  REAL(kind=r8),PARAMETER :: con_rocp   =con_rd/con_cp
  REAL(kind=r8),PARAMETER :: con_rerth  =6.3712e+6      ! radius of earth   (m)
  REAL(kind=r8),PARAMETER :: con_pi     =3.1415926535897931 ! pi

  REAL(kind=r8),PARAMETER :: CRITAC=0.10E-5_r8

  INTEGER, PARAMETER  :: nmtvr=14!     nmtvr    - integer, number of topographic variables such as  1    !
  !                         variance etc used in the GWD parameterization !

  REAL(kind=r8) :: cdmbgwd(2) ! Mtn Blking and GWD tuning factors
  !cdmbgwd  - real, multiplication factors for cdmb and gwd 2    !
  INTEGER :: IMX
  REAL(kind=r8), ALLOCATABLE :: HPRIME(:,:,:)
  INTEGER :: nlons         !     nlons(im)    - integer, number of total grid points in a latitude     !
  !                         circle through a point                   im   !
  !     lonf,
  INTEGER ::  latg  !latg- integer, number of lon/lat points                 1    !

  PUBLIC :: InitGwdd_ECMWF
  PUBLIC :: Run_Gwdd_ECMWF
CONTAINS

  SUBROUTINE InitGwdd_ECMWF(iMax,jMax,kMax,ibMax,jbMax,nfprt,&
                            nfvar,reducedGrid,fNameHPRIME,ibMaxPerJB)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: iMax
    INTEGER      , INTENT(IN   ) :: jMax
    INTEGER      , INTENT(IN   ) :: kMax
    INTEGER      , INTENT(IN   ) :: ibMax
    INTEGER      , INTENT(IN   ) :: jbMax
    INTEGER      , INTENT(IN   ) :: nfprt
    INTEGER      , INTENT(IN   ) :: nfvar
    LOGICAL      , INTENT(IN   ) :: reducedGrid
    CHARACTER(LEN=*), INTENT(IN   ) ::  fNameHPRIME
    INTEGER      , INTENT(IN   ) :: ibMaxPerJB(jbMax)
    REAL(KIND=r4) ::   buffer (iMax,jMax)
    REAL(KIND=r8) ::   var_in (iMax,jMax)
    REAL(KIND=r8) ::   var    (ibMax,jbMax)
    INTEGER :: k,i,j
    INTEGER :: LRecIn
    INTEGER :: ierr
    
    ALLOCATE(HPRIME(ibMax,jbMax,nmtvr));HPRIME=0.0_r8
    nlons=iMax
    latg =jMax
    
    buffer=0.0_r8
    INQUIRE (IOLENGTH=LRecIn) buffer 
    OPEN (UNIT=nfvar,FILE=TRIM(fNameHPRIME),FORM='UNFORMATTED', ACCESS='DIRECT', &
         ACTION='READ',RECL=LRecIn,STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(fNameHPRIME), ierr
       STOP "**(ERROR)**"
    END IF

    DO k=1,nmtvr
       READ(nfvar,rec=k) buffer
       var_in=REAL(buffer,KIND=r8)
       IF (reducedGrid) THEN
          CALL LinearIJtoIBJB(var_in,var)
       ELSE
          CALL IJtoIBJB(var_in,var)
       END IF
       DO j=1,jbMax 
          DO i=1,ibMaxPerJB(j)
             HPRIME(i,j,k) = var (i,j)
          END DO
       END DO
    END DO
    CLOSE(nfvar,STATUS='KEEP')
    IMX=iMax
    cdmbgwd(1)       = 1.0_r8       ! Mtn Blking and GWD tuning factors
    cdmbgwd(2)       = 1.0_r8       ! Mtn Blking and GWD tuning factors

    !  move water from vapor to liquid should the liquid amount be negative


  END SUBROUTINE InitGwdd_ECMWF

  SUBROUTINE Run_Gwdd_ECMWF(      &
       IM       , &
       KM       , &
       latco    , &
       dudt     , &
       dvdt     , &
       prsi ,prsl  ,&
       U1       , &
       V1       , &
       T1       , &
       Q1       , &
       pblh     , &
       cu_hr    , &
       cu_kbot  , &
       cu_ktop  , &
       cu_Kuo   , &
       colrad   , &
       DELTIM   )
    IMPLICIT NONE  
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: km
    INTEGER      , INTENT(IN   ) :: latco
    REAL(kind=r8), INTENT(IN   ) :: pblh(IM)                 ! Index for the PBL top layer!
    REAL(kind=r8), INTENT(IN   ) :: deltim

    REAL(kind=r8), INTENT(INOUT) :: dudt (IM,KM)
    REAL(kind=r8), INTENT(INOUT) :: dvdt (IM,KM)
    REAL(KIND=r8), INTENT(in   ) :: prsi   (im,km+1)  !     prsi     - real, pressure at layer interfaces [Pa]
    REAL(KIND=r8), INTENT(in   ) :: prsl   (im,km)    !     prsl     - real, mean layer presure [Pa]
    !REAL(KIND=r8), INTENT(in   ) :: phii   (im,km+1) !===>  PHIH(K+1)  INPUT GEOPOTENTIAL @ EDGES  IN MKS units (m)
    !REAL(KIND=r8), INTENT(in   ) :: phil   (im,km)   !===>  PHIL(K) INPUT GEOPOTENTIAL @ LAYERS IN MKS units (m)
    REAL(kind=r8), INTENT(IN   ) :: U1(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: V1(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: T1(IM,KM)  
    REAL(kind=r8), INTENT(IN   ) :: Q1(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: cu_hr    (im,km)
    INTEGER      , INTENT(IN   ) :: cu_kbot  (IM)
    INTEGER      , INTENT(IN   ) :: cu_ktop  (IM)
    INTEGER      , INTENT(IN   ) :: cu_Kuo   (IM)
    REAL(kind=r8), INTENT(IN   ) :: colrad   (IM)

    REAL(kind=r8) :: pgr  (IM)    ! pgr are  pressures(Pa)
    REAL(kind=r8) :: OC     (IM)
    REAL(kind=r8) :: OA4    (IM,4)
    REAL(kind=r8) :: CLX    (IM,4)
    REAL(kind=r8) :: ELVMAX2(IM)
    REAL(kind=r8) :: THETA  (IM)
    REAL(kind=r8) :: SIGMA  (IM)
    REAL(kind=r8) :: GAMMA  (IM)
    REAL(kind=r8) :: OroStdDev(IM)
    REAL(kind=r8) :: PHIL   (IM,KM)  !phil     - real, layer geopotential height ix,levs !
    REAL(kind=r8) :: PHII   (IM,KM+1)!phii     - real, interface geopotential height ix,levs+1
    !REAL(kind=r8) :: prsi   (IM,KM+1)!         - real, pressure at layer interfaces             ix,levs+1
    !REAL(kind=r8) :: prsl   (IM,KM)  !         - real, mean layer presure                       ix,levs !
    REAL(kind=r8) :: prslk  (IM,KM)   !     prslk    - real, Exner function at layer  		nCols,kMax !
    REAL(kind=r8) :: pgrk   (IM)
    REAL(kind=r8) :: DEL    (IM,KM) 
    REAL(kind=r8) :: DeltaP (IM,KM) 
    REAL(kind=r8) :: slk  (IM,KM)

    REAL(kind=r8) :: sikp1(iM,kM+1)

    REAL(kind=r8) :: DPHI
    REAL(kind=r8) :: TEM
    INTEGER       :: KPBL   (IM)                 ! Index for the PBL top layer!
    REAL(KIND=r8) :: ze     (IM,KM)
    REAL(KIND=r8) :: tv     (IM,KM)
    REAL(KIND=r8) :: delz   (IM,KM)
    REAL(KIND=r8) :: gwdcu  (IM,KM)
    REAL(KIND=r8) :: gwdcv  (IM,KM)
    REAL(KIND=r8) :: rbyg
    REAL (kind=r8), PARAMETER :: pt01=0.01_r8
    INTEGER :: i
    INTEGER :: k
    DO k=1,km
       DO i = 1,im
         DeltaP(i,k) = ((prsi(i,k)) - (prsi(i,k+1)))/prsi(i,1)
       END DO
    END DO

    DO k=1,km+1
       DO i = 1,im
          sikp1(i,k)   = ( prsi(i,k)/prsi(i,1)) ** (1.0_r8 + con_rocp)
       END DO
    END DO

    DO k=1,km
       DO i = 1,im
         slk(i,k)   = (sikp1(i,k)-sikp1(i,k+1))/((1.0_r8 + con_rocp) * (( prsi(i,k)/prsi(i,1)) - ( prsi(i,k+1)/prsi(i,1))))
       END DO
    END DO

    !
    !-- Terrain-specific inputs:
    !
    !          HPRIME(IM  )=HSTDV (I,J)   ! 1 -- standard deviation of orography
    !          OC    (IM  )=HCNVX (I,J)   ! 2 -- Normalized convexity
    !          OA4   (IM,1)=HASYW (I,J)   ! 3 -- orographic asymmetry in W-E plane
    !          OA4   (IM,2)=HASYS (I,J)   ! 4 -- orographic asymmetry in S-N plane
    !          OA4   (IM,3)=HASYSW(I,J)   ! 5 -- orographic asymmetry in SW-NE plane
    !          OA4   (IM,4)=HASYNW(I,J)   ! 6 -- orographic asymmetry in NW-SE plane
    !          CLX4  (IM,1)=HLENW (I,J)   ! 7 -- orographic length scale in W-E plane
    !          CLX4  (IM,2)=HLENS (I,J)   ! 8 -- orographic length scale in S-N plane
    !          CLX4  (IM,3)=HLENSW(I,J)   ! 9 -- orographic length scale in SW-NE plane
    !          CLX4  (IM,4)=HLENNW(I,J)   !10 -- orographic length scale in NW-SE plane
    !          THETA (IM  )=HANGL (I,J)   !11
    !          SIGMA (IM  )=HSLOP (I,J)   !12
    !          GAMMA (IM  )=HANIS (I,J)   !13
    !          ELVMAX2(IM  )=HZMAX (I,J)   !14

    DO i=1,im
       pgr(i) = prsi(i,1)
       !prsi(i,km+1)  = si(km+1)*pgr(i)           ! prsi are now pressures(Pa)
       pgrk(i)         = (pgr(i)*pt01) ** con_rocp
    ENDDO

    DO k=1,km
       DO i=1,im
          tv(i,k)=T1(i,k)*(1.0_r8+0.608_r8*Q1(i,k))
          !prsi(i,k)  = si(k)*pgr(i)               ! prsi are now pressures(Pa)
          !prsl(i,k)  = sl(k)*pgr(i)               ! vertical column of model pressure (Pa)
          prslk(i,k) = slk(i,k) * pgrk(i)
       ENDDO
    ENDDO

    DO i=1,im
       phii(i,1) = 0.0_r8           ! Ignoring topography height here
    END DO
    ! --- interpolate to max mtn height for index, iwklm(I) wk[gz]
    ! --- ELVMAX is limited to hncrit because to hi res topo30 orog.
    !             pkp1log =  phil(j,k+1) / con_g
    !             pklog   =  phil(j,k)   / con_g
    DO k=1,KM
       DO i=1,im
          TEM         = con_rd * T1(i,k)*(1.0_r8 + con_fvirt*MAX(Q1(i,k),0.0_r8))
          DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM/(PRSI(i,k) + PRSI(i,k+1))
          phil(i,k)   = phii(i,k) + DPHI
          phii(i,k+1) = phil(i,k) + DPHI
       ENDDO
    ENDDO
    DO k=1,KM
       DO i=1,im
          del(i,k) = PRSI(i,k) - PRSI(i,k+1)
       ENDDO
    ENDDO

    ze=0.0_r8
    DO i=1,im
      rbyg=con_rd/con_g*DeltaP(i,1)*0.5e0_r8
      delz(i,1)=MAX((rbyg * tv(i,1)),0.5_r8)      
      ze  (i,1)= delz(i,1)         !gt(i,1)*(con_rd/con_g)*(psur(i)-prsl(i,k)(i,1))/psur(i)
      ze  (i,1)= ze(i,1)
    END DO
    
    DO k=2,KM
       DO i=1,im
          delz(i,k)=0.5_r8*con_rd*(tv(i,k-1)+tv(i,k))* &
               LOG(prsl(i,k-1)/prsl(i,k))/con_g
          ze(i,k)=ze(i,k-1)+ delz(i,k)
       END DO
    END DO  
    KPBL=1
    DO k=1,KM
       DO i=1,im
          IF(ze(i,k) <= pblh(i))THEN
             KPBL(i)=k
          END IF 
       END DO
    END DO

    IF (nmtvr == 6) THEN

       DO i = 1, im
          OroStdDev(i) = hprime(i,latco,1)
          oc(i) = hprime(i,latco,2)
       ENDDO

       DO k = 1, 4
          DO i = 1, im
             oa4(i,k) = hprime(i,latco,k+2)
             clx(i,k) = 0.0_r8
          ENDDO
       ENDDO

    ELSEIF (nmtvr == 10) THEN

       DO i = 1, im
          OroStdDev(i) = hprime(i,latco,1)
          oc(i) = hprime(i,latco,2)
       ENDDO

       DO k = 1, 4
          DO i = 1, im
             oa4(i,k) = hprime(i,latco,k+2)
             clx(i,k) = hprime(i,latco,k+6)
          ENDDO
       ENDDO

    ELSEIF (nmtvr == 14) THEN

       DO i = 1, im
          OroStdDev(i) = hprime(i,latco,1)
          oc(i)        = hprime(i,latco,2)
       ENDDO

       DO k = 1, 4
          DO i = 1, im
             oa4(i,k) = hprime(i,latco,k+2)
             clx(i,k) = hprime(i,latco,k+6)
          ENDDO
       ENDDO

       DO i = 1, im
          theta(i)   = hprime(i,latco,11)
          gamma(i)   = hprime(i,latco,12)
          sigma(i)   = hprime(i,latco,13)
          ELVMAX2(i) = hprime(i,latco,14)
       ENDDO

    ELSE
       DO i = 1, im
          OroStdDev(i)  = 0.0_r8
          oc    (i)  = 0.0_r8
          theta (i)  = 0.0_r8
          gamma (i)  = 0.0_r8
          sigma (i)  = 0.0_r8
          ELVMAX2(i)  = 0.0_r8
       ENDDO
       DO k = 1, 4
          DO i = 1, im
             oa4   (i,k)  = 0.0_r8
             clx   (i,k)  = 0.0_r8
          ENDDO
       ENDDO

    ENDIF   ! end if_nmtvr

    CALL GWDPS( &
         IM                    , &
         KM                    , &
         dudt   (1:IM,1:KM)    , &
         dvdt   (1:IM,1:KM)    , &
         U1     (1:IM,1:KM)    , &
         V1     (1:IM,1:KM)    , &
         T1     (1:IM,1:KM)    , &
         Q1     (1:IM,1:KM)    , &
         KPBL   (1:IM)         , &
         PRSI   (1:IM,1:KM+1)  , &
         DEL    (1:IM,1:KM)    , &
         PRSL   (1:IM,1:KM)    , &
         PRSLK  (1:IM,1:KM)    , &
         PHII   (1:IM,1:KM+1)  , &
         PHIL   (1:IM,1:KM)    , &
         DELTIM                , &
         OroStdDev(1:IM)       , &
         OC     (1:IM)         , &
         OA4    (1:IM,1:4)     , &
         CLX    (1:IM,1:4)     , &
         THETA  (1:IM)         , &
         SIGMA  (1:IM)         , &
         GAMMA  (1:IM)         , &
         ELVMAX2(1:IM)           )


    CALL gbphys(&
       im         , &
       km         , &
       cu_hr      , &
       cu_kbot    , &
       cu_ktop    , &
       cu_Kuo     , &
       colrad     , &
       prsl       , &
       prsi       , &
       pgr        , &
       U1         , &
       V1         , &
       T1         , &
       Q1         , &
       gwdcu      , &
       gwdcv        )

!     dtf      - real, dynamics time step in seconds               1    !
!     dtp      - real, physics time step in seconds                1    !
    DO k=1,KM
       DO i=1,im
       dudt(i,k)= dudt(i,k) + gwdcu(i,k)
       dvdt(i,k)= dvdt(i,k) + gwdcv(i,k)
       
        IF(ABS(dudt(i,k)) > CRITAC )THEN
           dudt(i,k)= dudt(i,k)*CRITAC
        ELSE
           dudt(i,k)= dudt(i,k)
        END IF
        IF(ABS(dvdt(i,k)) > CRITAC )THEN
           dvdt(i,k)= dvdt(i,k)*CRITAC
        ELSE
           dvdt(i,k)= dvdt(i,k)
        END IF
       END DO
    END DO  

  END SUBROUTINE Run_Gwdd_ECMWF








  SUBROUTINE gbphys(im,km,cuhr,kbot,ktop,Kuo,colrad,prsl,prsi,pgr,ugrs,&
       vgrs,tgrs,qgrs,gwdcu,gwdcv)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: im
    INTEGER      , INTENT(IN   ) :: km
    REAL(KIND=r8), INTENT(IN   ) :: cuhr  (im,km)!           cuhr = temperature change due to deep convection
    REAL(KIND=r8), INTENT(IN   ) :: prsl (im,km)
    REAL(KIND=r8), INTENT(IN   ) :: prsi (im,km+1)
    INTEGER      , INTENT(IN   ) :: kbot  (im)
    INTEGER      , INTENT(IN   ) :: ktop  (im)
    INTEGER      , INTENT(IN   ) :: Kuo   (im)
    REAL(KIND=r8), INTENT(IN   ) :: colrad(im)
    REAL(KIND=r8), INTENT(IN   ) :: pgr   (im)
    REAL(KIND=r8), INTENT(IN   ) :: ugrs  (im,km)
    REAL(KIND=r8), INTENT(IN   ) :: vgrs  (im,km)
    REAL(KIND=r8), INTENT(IN   ) :: tgrs  (im,km)
    REAL(KIND=r8), INTENT(IN   ) :: qgrs  (im,km)
    REAL(kind=r8), INTENT(OUT  ) :: gwdcu (im,km)
    REAL(kind=r8), INTENT(OUT  ) :: gwdcv (im,km)


    REAL(kind=r8) :: qmax (im) !--- ...  calculate maximum convective heating rate            qmax [k/s]
    REAL(kind=r8) :: cumabs(im)
    REAL(kind=r8) :: coslat(im)
    REAL(kind=r8) :: dlength(im)
    REAL(kind=r8) :: cumchr(im,km)
    !REAL(kind=r8) :: prsi  (im,km+1)
    !REAL(kind=r8) :: prsl  (im,km)
    REAL(kind=r8) :: del   (im,km)
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    INTEGER       :: i,k,k1

    DO i = 1, im
       qmax(i)   = 0.0_r8
       cumabs(i) = 0.0_r8
    ENDDO
    DO k = 1, km
       DO i = 1, im
          !cuhr(i,k) = (gt0(i,k)-dtdt(i,k)) / dtp    ! moorthi
          cumchr(i,k)  = 0.0_r8
          gwdcu (i,k)  = 0.0_r8
          gwdcv (i,k)  = 0.0_r8
          IF (k >= kbot(i) .AND. k <= ktop(i)) THEN
             qmax(i)     = MAX(qmax(i),cuhr(i,k))
             cumabs(i)   = cuhr(i,k) + cumabs(i)
          ENDIF
       ENDDO
    ENDDO

    DO i = 1, im
       DO k = kbot(i), ktop(i)
          DO k1 = kbot(i), k
             cumchr(i,k) = cuhr(i,k1) + cumchr(i,k)
          ENDDO
          IF(cumabs(i) == 0.0_r8)THEN
             cumchr(i,k) = 0.0_r8
          ELSE
             cumchr(i,k) = cumchr(i,k) / cumabs(i)
          END IF 
       ENDDO
    ENDDO

    DO i = 1, im
       ! colrad.....colatitude  colrad=0 - 3.14 (0-180)from np to sp in radians
       !IF((((colrad(i)*180.0_r8)/3.1415926e0_r8)-90.0_r8)  > 0.0_r8 ) THEN
       coslat(i)   = COS(((colrad(i)))-(3.1415926e0_r8/2.0_r8))
       tem1        = con_rerth * (con_pi+con_pi)*coslat(i)/nlons
       tem2        = con_rerth *  con_pi/latg
       dlength(i)  = SQRT( tem1*tem1+tem2*tem2 )
    ENDDO

    !DO i=1,im
    !   prsi(i,km+1)  = MAX(si(km+1),10e-12_r8)*pgr(i)           ! prsi are now pressures(Pa)
    !ENDDO

    !DO k=1,km
    !   DO i=1,im
    !      prsi(i,k)  = si(k)*pgr(i)               ! prsi are now pressures(Pa)
    !      prsl(i,k)  = sl(k)*pgr(i)               ! vertical column of model pressure (Pa)
    !   ENDDO
    !ENDDO

    DO k=1,KM
       DO i=1,im
          del(i,k) = prsi(i,k) - prsi(i,k+1)
       ENDDO
    ENDDO

    CALL gwdc(im                , &
         km                     , &
         ugrs      (1:im,1:km)  , &
         vgrs      (1:im,1:km)  , &
         tgrs      (1:im,1:km)  , &
         qgrs      (1:im,1:km)  , &
         prsl      (1:im,1:km)  , &
         prsi      (1:im,1:km+1), &
         del       (1:im,1:km)  , &
         qmax      (1:im)       , &
         cumchr    (1:im,1:km)  ,  &
         ktop      (1:im)       , &
         kbot      (1:im)       , &
         Kuo       (1:im)       , &
         dlength   (1:im)       , &
         gwdcu     (1:im,1:km)  , &
         gwdcv     (1:im,1:km)    )

  END SUBROUTINE gbphys


  SUBROUTINE GWDPS( &
       IM       , &
       KM       , &
       dudt     , &
       dvdt     , &
       U1       , &
       V1       , &
       T1       , &
       Q1       , &
       KPBL     , &
       PRSI     , &
       DEL      , &
       PRSL     , &
       PRSLK    , &
       PHII     , &
       PHIL     , &
       DELTIM   , &
       HPRIME   , &
       OC       , &
       OA4      , &
       CLX4     , &
       THETA    , &
       SIGMA    , &
       GAMMA    , &
       ELVMAX2     )
    !
    !   ********************************************************************
    ! ----->  I M P L E M E N T A T I O N    V E R S I O N   <----------
    !
    !          --- Not in this code --  History of GWDP at NCEP----
    !              ----------------     -----------------------
    !  VERSION 3  MODIFIED FOR GRAVITY WAVES, LOCATION: .FR30(V3GWD)  *J*
    !---       3.1 INCLUDES VARIABLE SATURATION FLUX PROFILE CF ISIGST
    !---       3.G INCLUDES PS COMBINED W/ PH (GLAS AND GFDL)
    !-----         ALSO INCLUDED IS RI  SMOOTH OVER A THICK LOWER LAYER
    !-----         ALSO INCLUDED IS DECREASE IN DE-ACC AT TOP BY 1/2
    !-----     THE NMC GWD INCORPORATING BOTH GLAS(P&S) AND GFDL(MIGWD)
    !-----        MOUNTAIN INDUCED GRAVITY WAVE DRAG 
    !-----    CODE FROM .FR30(V3MONNX) FOR MONIN3
    !-----        THIS VERSION (06 MAR 1987)
    !-----        THIS VERSION (26 APR 1987)    3.G
    !-----        THIS VERSION (01 MAY 1987)    3.9
    !-----    CHANGE TO FORTRAN 77 (FEB 1989)     --- HANN-MING HENRY JUANG
    !-----    20070601 ELVMAX bug fix (*j*)
    !
    !   VERSION 4
    !                ----- This code -----
    !
    !-----   MODIFIED TO IMPLEMENT THE ENHANCED LOW TROPOSPHERIC GRAVITY
    !-----   WAVE DRAG DEVELOPED BY KIM AND ARAKAWA(JAS, 1995).
    !        Orographic Std Dev (hprime), Convexity (OC), Asymmetry (OA4)
    !        and Lx (CLX4) are input topographic statistics needed.
    !
    !-----   PROGRAMMED AND DEBUGGED BY HONG, ALPERT AND KIM --- JAN 1996.
    !-----   debugged again - moorthi and iredell --- may 1998.
    !-----
    !       Further Cleanup, optimization and modification
    !                                       - S. Moorthi May 98, March 99.
    !-----   modified for usgs orography data (ncep office note 424)
    !        and with several bugs fixed  - moorthi and hong --- july 1999.
    !
    !-----   Modified & implemented into NRL NOGAPS
    !                                       - Young-Joon Kim, July 2000
    !-----
    !   VERSION lm MB  (6): oz fix 8/2003
    !                ----- This code -----
    !
    !------   Changed to include the Lott and Miller Mtn Blocking
    !         with some modifications by (*j*)  4/02
    !        From a Principal Coordinate calculation using the
    !        Hi Res 8 minute orography, the Angle of the
    !        mtn with that to the East (x) axis is THETA, the slope
    !        parameter SIGMA. The anisotropy is in GAMMA - all  are input
    !        topographic statistics needed.  These are calculated off-line
    !        as a function of model resolution in the fortran code ml01rg2.f,
    !        with script mlb2.sh.   (*j*)
    !-----   gwdps_mb.f version (following lmi) elvmax < hncrit (*j*)
    !        MB3a expt to enhance elvmax mtn hgt see sigfac & hncrit
    !        gwdps_GWDFIX_v6.f FIXGWD GF6.0 20070608 sigfac=4.
    !-----
    !----------------------------------------------------------------------C
    !    USE
    !        ROUTINE IS CALLED FROM GBPHYS  (AFTER CALL TO MONNIN)
    !
    !    PURPOSE
    !        USING THE GWD PARAMETERIZATIONS OF PS-GLAS AND PH-
    !        GFDL TECHNIQUE.  THE TIME TENDENCIES OF U V
    !        ARE ALTERED TO INCLUDE THE EFFECT OF MOUNTAIN INDUCED
    !        GRAVITY WAVE DRAG FROM SUB-GRID SCALE OROGRAPHY INCLUDING
    !        CONVECTIVE BREAKING, SHEAR BREAKING AND THE PRESENCE OF
    !        CRITICAL LEVELS
    !
    !  INPUT
    !        A(IM,KM)  NON-LIN TENDENCY FOR V WIND COMPONENT
    !        B(IM,KM)  NON-LIN TENDENCY FOR U WIND COMPONENT
    !        U1(IM,KM) ZONAL WIND M/SEC  AT T0-DT
    !        V1(IM,KM) MERIDIONAL WIND M/SEC AT T0-DT
    !        T1(IM,KM) TEMPERATURE DEG K AT T0-DT
    !        Q1(IM,KM) SPECIFIC HUMIDITY AT T0-DT
    !
    !        DELTIM  TIME STEP    SECS
    !        SI(N)   P/PSFC AT BASE OF LAYER N
    !        SL(N)   P/PSFC AT MIDDLE OF LAYER N
    !        DEL(N)  POSITIVE INCREMENT OF P/PSFC ACROSS LAYER N
    !        KPBL(IM) is the index of the top layer of the PBL
    !        ipr &  for diagnostics
    !
    !  OUTPUT
    !        A, B    AS AUGMENTED BY TENDENCY DUE TO GWDPS
    !                OTHER INPUT VARIABLES UNMODIFIED.
    !   ********************************************************************
    IMPLICIT NONE

    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: km
    INTEGER      , INTENT(IN   ) :: KPBL(IM)                 ! Index for the PBL top layer!
    REAL(kind=r8), INTENT(IN   ) :: deltim

    REAL(kind=r8), INTENT(INOUT) :: dudt (IM,KM) 
    REAL(kind=r8), INTENT(INOUT) :: dvdt (IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: U1(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: V1(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: T1(IM,KM)  
    REAL(kind=r8), INTENT(IN   ) :: Q1(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: PRSI(IM,KM+1)
    REAL(kind=r8), INTENT(IN   ) :: DEL(IM,KM) 
    REAL(kind=r8), INTENT(IN   ) :: PRSL(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: PRSLK(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: PHIL(IM,KM)
    REAL(kind=r8), INTENT(IN   ) :: PHII(IM,KM+1)

    REAL(kind=r8), INTENT(IN   ) :: OC    (IM)
    REAL(kind=r8), INTENT(IN   ) :: OA4   (IM,4)
    REAL(kind=r8), INTENT(IN   ) :: CLX4  (IM,4)
    REAL(kind=r8), INTENT(IN   ) :: HPRIME(IM)

    REAL(kind=r8), INTENT(IN   ) :: ELVMAX2(IM)
    REAL(kind=r8), INTENT(IN   ) :: THETA(IM)
    REAL(kind=r8), INTENT(IN   ) :: SIGMA(IM)
    REAL(kind=r8), INTENT(IN   ) :: GAMMA(IM)

      !
    !     Some constants
    !
    REAL(kind=r8), PARAMETER :: pi=3.1415926535897931_r8
    REAL(kind=r8), PARAMETER :: rad_to_deg=180.0_r8/PI
    REAL(kind=r8), PARAMETER :: deg_to_rad=PI/180.0_r8
    
    REAL(kind=r8), PARAMETER :: DW2MIN=1.0_r8
    REAL(kind=r8), PARAMETER :: RIMIN=-100.0_r8
    REAL(kind=r8), PARAMETER :: RIC=0.25_r8
    REAL(kind=r8), PARAMETER :: BNV2MIN=1.0E-5_r8
    REAL(kind=r8), PARAMETER :: EFMAX=10.0_r8
    REAL(kind=r8), PARAMETER :: EFMIN=0.0_r8
    real(kind=r8), PARAMETER :: hpmax=200.0_r8 ! max standard deviation of orography
    !REAL(kind=r8), PARAMETER :: hpmax=1000.0_r8
    !REAL(kind=r8), PARAMETER :: hpmax=2400.0_r8
    REAL(kind=r8), PARAMETER :: hpmin=1.0_r8   ! min standard deviation of orography
    !
    REAL(kind=r8), PARAMETER :: FRC=1.0_r8
    REAL(kind=r8), PARAMETER :: CE=0.8_r8
    REAL(kind=r8), PARAMETER :: CEOFRC=CE/FRC
    REAL(kind=r8), PARAMETER :: frmax=100.0_r8
    REAL(kind=r8), PARAMETER :: CG=0.5_r8
    REAL(kind=r8), PARAMETER :: GMAX=1.0_r8
    REAL(kind=r8), PARAMETER :: VELEPS=1.0_r8
    REAL(kind=r8), PARAMETER :: FACTOP=0.5_r8
    REAL(kind=r8), PARAMETER :: RLOLEV=50000.0_r8
    !     real(kind=r8), PARAMETER :: parameter (RLOLEV=500.0_r8) 
    !     real(kind=r8), PARAMETER :: parameter (RLOLEV=0.5_r8)
    !
    ! --- for lm mtn blocking
    !     real(kind=r8), PARAMETER :: cdmb = 1.0_r8     ! non-dim sub grid mtn drag Amp (*j*)
    !REAL(kind=r8), PARAMETER :: hncrit=10000.0_r8   ! Max value in meters for ELVMAX (*j*)

    REAL(kind=r8), PARAMETER :: hncrit=8000.0_r8   ! Max value in meters for ELVMAX (*j*)
                                                   ! hncrit set to 8000m and sigfac added to
                                                   ! enhance elvmax mtn hgt
!    REAL(kind=r8), PARAMETER :: sigfac=4.0_r8     ! MB3a expt test for ELVMAX factor (*j*)
    REAL(kind=r8), PARAMETER :: sigfac=4.0_r8      ! MB3a expt test for ELVMAX factor (*j*)

    REAL(kind=r8), PARAMETER :: hminmt=50.0_r8     ! min mtn height (*j*)
    REAL(kind=r8), PARAMETER :: minwnd=1.1_r8      ! min wind component (*j*)

    !real(kind=r8), PARAMETER :: parameter (dpmin=00.0_r8)     ! Minimum thickness of the reference layer
    !real(kind=r8), PARAMETER :: parameter (dpmin=05.0_r8)     ! Minimum thickness of the reference layer
    !real(kind=r8), PARAMETER :: parameter (dpmin=20.0_r8)     ! Minimum thickness of the reference layer
    ! in centibars
    !REAL(kind=r8), PARAMETER :: dpmin=5000.0_r8   ! Minimum thickness of the reference layer in Pa
    REAL(kind=r8), PARAMETER :: dpmin=50.0_r8      ! Minimum thickness of the reference layer in Pa

    !
    INTEGER      , PARAMETER :: mdir=8
    REAL(kind=r8), PARAMETER :: FDIR=mdir/(PI+PI)

    INTEGER      , PARAMETER :: nwdir(1:mdir)=(/6,7,5,8,2,3,1,4/)
    !
    LOGICAL                  :: ICRILV(IM)
    !
    !----   MOUNTAIN INDUCED GRAVITY WAVE DRAG
    !
  ! for lm mtn blocking
    REAL(kind=r8) :: ELVMAX(IM)
    REAL(kind=r8) :: wk(IM)
    REAL(kind=r8) :: bnv2lm(IM,KM)
    REAL(kind=r8) :: PE(IM)
    REAL(kind=r8) :: EK(IM)
    REAL(kind=r8) :: ZBK(IM)
    REAL(kind=r8) :: UP(IM)
    REAL(kind=r8) :: DB(IM,KM)
    REAL(kind=r8) :: ANG(IM,KM)
    REAL(kind=r8) :: UDS(IM,KM)
    REAL(kind=r8) :: ZLEN
    REAL(kind=r8) :: DBTMP
    REAL(kind=r8) :: R
    REAL(kind=r8) :: PHIANG
    REAL(kind=r8) :: CDmb
    REAL(kind=r8) :: DBIM
    REAL(kind=r8) :: RDI
    REAL(kind=r8) :: DUSFC(IM)
    REAL(kind=r8) :: DVSFC(IM)
    REAL(kind=r8) :: TAUB(IM)
    REAL(kind=r8) :: XN(IM)
    REAL(kind=r8) :: YN(IM)
    REAL(kind=r8) :: UBAR(IM)
    REAL(kind=r8) :: VBAR(IM)
    REAL(kind=r8) :: ULOW(IM)
    REAL(kind=r8) :: OA(IM)
    REAL(kind=r8) :: CLX(IM)
    REAL(kind=r8) :: ROLL(IM)
    REAL(kind=r8) :: ULOI(IM)
    REAL(kind=r8) :: DTFAC(IM)
    REAL(kind=r8) :: XLINV(IM)
    REAL(kind=r8) :: DELKS(IM)
    REAL(kind=r8) :: DELKS1(IM)
    REAL(kind=r8) :: BNV2(IM,KM)
    REAL(kind=r8) :: TAUP(IM,KM+1)
    REAL(kind=r8) :: ri_n(IM,KM)
    REAL(kind=r8) :: TAUD(IM,KM)
    REAL(kind=r8) :: RO(IM,KM)
    REAL(kind=r8) :: VTK(IM,KM)
    REAL(kind=r8) :: VTJ(IM,KM)
    REAL(kind=r8) :: SCOR(IM)
    REAL(kind=r8) :: VELCO(IM,KM-1)
    REAL(kind=r8) :: bnv2bar(im)
    !
    INTEGER       :: kref(IM)
    INTEGER       :: kint(im)
    INTEGER       :: iwk(im)
    INTEGER       :: ipt(im)
    ! for lm mtn blocking
    INTEGER       :: kreflm(IM)
    INTEGER       :: iwklm(im)
    INTEGER       :: idxzb(im)
    INTEGER       :: ktrial
    INTEGER       :: klevm1
    !
    REAL(kind=r8) :: gor
    REAL(kind=r8) :: gocp
    REAL(kind=r8) :: fv
    REAL(kind=r8) :: gr2
    REAL(kind=r8) :: bnv
    REAL(kind=r8) :: fr
    REAL(kind=r8) :: brvf
    REAL(kind=r8) :: Gamma_Eff
    REAL(kind=r8) :: tem
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(kind=r8) :: temc
    REAL(kind=r8) :: temv
    REAL(kind=r8) :: wdir
    REAL(kind=r8) :: ti
    REAL(kind=r8) :: rdz
    REAL(kind=r8) :: dw2
    REAL(kind=r8) :: shr2
    REAL(kind=r8) :: bvf2   
    REAL(kind=r8) :: rdelks
    REAL(kind=r8) :: efact
    REAL(kind=r8) :: coefm
    REAL(kind=r8) :: gfobnv   
    REAL(kind=r8) :: scork
    REAL(kind=r8) :: rscor
    REAL(kind=r8) :: hd
    REAL(kind=r8) :: fro
    REAL(kind=r8) :: rim
    REAL(kind=r8) :: sira   
    REAL(kind=r8) :: dtaux
    REAL(kind=r8) :: dtauy
    REAL(kind=r8) :: pkp1log
    REAL(kind=r8) :: pklog
    INTEGER       :: kmm1
    INTEGER       :: kmm2
    INTEGER       :: lcap
    INTEGER       :: lcapp1
    INTEGER       :: kbps
    INTEGER       :: kbpsp1
    INTEGER       :: kbpsm1
    INTEGER       :: kmps
    INTEGER       :: idir
    INTEGER       :: nwd
    INTEGER       :: i
    INTEGER       :: j
    INTEGER       :: k
    INTEGER       :: klcap
    INTEGER       :: kp1
    INTEGER       :: kmpbl
    INTEGER       :: npt
    INTEGER       :: npr
    INTEGER       :: kmll
    INTEGER , PARAMETER :: ipr = 1
    ELVMAX=0.0_r8;    wk=0.0_r8;    bnv2lm=0.0_r8;    PE=0.0_r8
    EK=0.0_r8;    ZBK=0.0_r8;    UP=0.0_r8;    DB=0.0_r8
    ANG=0.0_r8;    UDS=0.0_r8;    ZLEN=0.0_r8;    DBTMP=0.0_r8
    R=0.0_r8;    PHIANG=0.0_r8;    CDmb=0.0_r8;    DBIM=0.0_r8
    RDI=0.0_r8;    DUSFC=0.0_r8;    DVSFC=0.0_r8;    TAUB=0.0_r8
    XN=0.0_r8;    YN=0.0_r8;    UBAR=0.0_r8;    VBAR=0.0_r8
    ULOW=0.0_r8;    OA=0.0_r8;    CLX=0.0_r8;    ROLL=0.0_r8
    ULOI=0.0_r8;    DTFAC=0.0_r8;    XLINV=0.0_r8;    DELKS=0.0_r8
    DELKS1=0.0_r8;    BNV2=0.0_r8;    TAUP=0.0_r8;    ri_n=0.0_r8
    TAUD=0.0_r8;    RO=0.0_r8;    VTK=0.0_r8;    VTJ=0.0_r8
    SCOR=0.0_r8;    VELCO=0.0_r8;    bnv2bar=0.0_r8;    gor=0.0_r8
    gocp=0.0_r8;fv=0.0_r8;gr2=0.0_r8;bnv=0.0_r8;    fr=0.0_r8;brvf=0.0_r8;Gamma_Eff=0.0_r8;tem=0.0_r8
    tem1=0.0_r8;tem2=0.0_r8;temc=0.0_r8;temv=0.0_r8;    wdir=0.0_r8;ti=0.0_r8;rdz=0.0_r8;dw2=0.0_r8
    shr2=0.0_r8;bvf2   =0.0_r8;rdelks=0.0_r8;    efact=0.0_r8;coefm=0.0_r8
    gfobnv   =0.0_r8;scork=0.0_r8;    rscor=0.0_r8;hd=0.0_r8
    fro=0.0_r8;rim=0.0_r8;    sira   =0.0_r8;dtaux=0.0_r8
    dtauy=0.0_r8;pkp1log=0.0_r8;pklog=0.0_r8
    kref=0;kint=0;iwk=0;ipt=0;kreflm=0;iwklm=0;idxzb=0
    ktrial=0;klevm1=0;kmm1=0;kmm2=0;    lcap=0;lcapp1=0;kbps=0;kbpsp1=0
    kbpsm1=0;kmps=0;idir=0;nwd=0;    i=0;j=0;k=0;klcap=0
    kp1=0;kmpbl=0;    npt=0;npr=0;    kmll=0
    
    ELVMAX=ELVMAX2

    !
    !     parameter (cdmb = 1.0_r8)     ! non-dim sub grid mtn drag Amp (*j*)
    ! non-dim sub grid mtn drag Amp (*j*)
    !     cdmb = 1.0_r8/float(IMX/192)
    !     cdmb = 192.0_r8/float(IMX)
    cdmb = 4.0_r8 * 192.0_r8/REAL(IMX)
    !cdmb = 4.0_r8 *  96.0_r8/REAL(IMX)
    !cdmb = 4.0_r8 * 48.0_r8/REAL(IMX)

    IF (cdmbgwd(1) >= 0.0_r8) cdmb = cdmb * cdmbgwd(1)
    !
    npr = 0
    DO I = 1, IM
       DUSFC(I) = 0.0_r8
       DVSFC(I) = 0.0_r8
    ENDDO
    !
    DO K = 1, KM
       DO I = 1, IM
          DB(I,K)  = 0.0_r8
          ANG(I,K) = 0.0_r8
          UDS(I,K) = 0.0_r8
       ENDDO
    ENDDO
    !
    RDI  = 1.0_r8 / con_rd
    GOR  = con_g/con_rd
    GR2  = con_g*GOR
    GOCP = con_g/con_cp
    FV   = con_rv/con_rd - 1
    !
    !     NCNT   = 0
    KMM1   = KM - 1
    KMM2   = KM - 2
    LCAP   = KM
    LCAPP1 = LCAP + 1
    !
    !
    IF ( NMTVR .EQ. 14) THEN 
       ! ----  for lm and gwd calculation points
       ipt = 0
       npt = 0
       DO I = 1,IM
          !                    50.0_r8                      1.0_r8
          IF ( (elvmax(i) .GT. HMINMT) .AND. (hprime(i) .GT. hpmin) )  THEN
             npt      = npt + 1
             ipt(npt) = i
             IF (ipr .EQ. i) npr = npt
          ENDIF
       ENDDO
       IF (npt .EQ. 0) RETURN     ! No gwd/mb calculation done!
       !
       !
       ! --- iwklm is the level above the height of the of the mountain.
       ! --- idxzb is the level of the dividing streamline.
       ! INITIALIZE DIVIDING STREAMLINE (DS) CONTROL VECTOR
       !
       DO i=1,npt
          iwklm(i)  = 2
          IDXZB(i)  = 0 
          kreflm(i) = 0
       ENDDO
       !
       ! start lm mtn blocking (mb) section
       !
       !..............................
       !..............................
       !
       !  (*j*)  11/03:  test upper limit on KMLL=km - 1
       !      then do not need hncrit -- test with large hncrit first.
       !       KMLL  = km / 2 ! maximum mtnlm height : # of vertical levels / 2
       KMLL = kmm1
       ! --- No mtn should be as high as KMLL (so we do not have to start at 
       ! --- the top of the model but could do calc for all levels).
       !
       !ELVMAX    -30000 <-> 30000
       !
       DO I = 1, npt
          j = ipt(i)
          ELVMAX(J) = MIN (ELVMAX(J) + sigfac * hprime(j), hncrit)
       ENDDO
       !
       DO K = 1,KMLL
          DO I = 1, npt
             j = ipt(i)
             ! --- interpolate to max mtn height for index, iwklm(I) wk[gz]
             ! --- ELVMAX is limited to hncrit because to hi res topo30 orog.
             pkp1log =  phil(j,k+1) / con_g
             pklog   =  phil(j,k)   / con_g
!!!-------     ELVMAX(J) = min (ELVMAX(J) + sigfac * hprime(j), hncrit)
             IF ( ( ELVMAX(j) .LE.  pkp1log ) .AND. ( ELVMAX(j) .GE.   pklog  ) ) THEN
                !     print *,' in gwdps_lm.f 1  =',k,ELVMAX(j),pklog,pkp1log,me
                ! ---        wk for diags but can be saved and reused.  
                wk(i)  = con_g * ELVMAX(j) / ( phil(j,k+1) - phil(j,k) )
                iwklm(I)  =  MAX(iwklm(I), k+1 ) 
                !     print *,' in gwdps_lm.f 2 npt=',npt,i,j,wk(i),iwklm(i),me
             ENDIF
             !
             ! ---        find at prsl levels large scale environment variables
             ! ---        these cover all possible mtn max heights
             VTJ(I,K)  = T1(J,K)  * (1.0_r8+FV*Q1(J,K))
             VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
             RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K) ! DENSITY Kg/M**3
          ENDDO
       ENDDO
       !
       ! testing for highest model level of mountain top
       !
       !         ihit = 2
       !         jhit = 0
       !        do i = 1, npt
       !        j=ipt(i)
       !          if ( iwklm(i) .gt. ihit ) then 
       !            ihit = iwklm(i)
       !            jhit = j
       !          endif
       !        enddo
       !     print *, ' mb: kdt,max(iwklm),jhit,phil,me=',
       !    &          kdt,ihit,jhit,phil(jhit,ihit),me

       klevm1 = KMLL - 1
       DO K = 1, klevm1  
          DO I = 1, npt
             j   = ipt(i)
             RDZ  = con_g   / ( phil(j,k+1) - phil(j,k) )
             ! ---                               Brunt-Vaisala Frequency
             BNV2LM(I,K) = (con_g+con_g) * RDZ * ( VTK(I,K+1)-VTK(I,K) ) &
                  / ( VTK(I,K+1)+VTK(I,K) )
             bnv2lm(i,k) = MAX( bnv2lm(i,k), bnv2min )
          ENDDO
       ENDDO
       !    print *,' in gwdps_lm.f 3 npt=',npt,j,RDZ,me
       !
       DO I = 1, npt
          J   = ipt(i)
          DELKS(I)  = 1.0_r8 / (PRSI(J,1) - PRSI(J,iwklm(i)))
          DELKS1(I) = 1.0_r8 / (PRSL(J,1) - PRSL(J,iwklm(i)))
          UBAR (I)  = 0.0_r8
          VBAR (I)  = 0.0_r8
          ROLL (I)  = 0.0_r8
          PE   (I)  = 0.0_r8
          EK   (I)  = 0.0_r8
          BNV2bar(I) = (PRSL(J,1)-PRSL(J,2)) * DELKS1(I) * BNV2LM(I,1)
       ENDDO

       ! --- find the dividing stream line height 
       ! --- starting from the level above the max mtn downward
       ! --- iwklm(i) is the k-index of mtn elvmax elevation
       DO Ktrial = KMLL, 1, -1
          DO I = 1, npt
             IF ( Ktrial .LT. iwklm(I) .AND. kreflm(I) .EQ. 0 ) THEN
                kreflm(I) = Ktrial
             ENDIF
          ENDDO
       ENDDO
       !
       ! --- in the layer kreflm(I) to 1 find PE (which needs N, ELVMAX)
       ! ---  make averages, guess dividing stream (DS) line layer.
       ! ---  This is not used in the first cut except for testing and
       ! --- is the vert ave of quantities from the surface to mtn top.
       !   
       DO I = 1, npt
          DO K = 1, Kreflm(I)
             J        = ipt(i)
             RDELKS     = DEL(J,K) * DELKS(I)
             UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K) ! trial Mean U below 
             VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K) ! trial Mean V below 
             ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K) ! trial Mean RO below 
             RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS1(I)
             BNV2bar(I) = BNV2bar(I) + BNV2lm(I,K) * RDELKS
             ! --- these vert ave are for diags, testing and GWD to follow (*j*).
          ENDDO
       ENDDO
       !
       ! --- integrate to get PE in the trial layer.
       ! --- Need the first layer where PE>EK - as soon as 
       ! --- IDXZB is not 0 we have a hit and Zb is found.
       !
       DO I = 1, npt
          J = ipt(i)
          DO K = iwklm(I), 1, -1
             PHIANG   =  ATAN2(V1(J,K),U1(J,K))*RAD_TO_DEG
             ANG(I,K) = ( THETA(J) - PHIANG )
             IF ( ANG(I,K) .GT.  90.0_r8 ) ANG(I,K) = ANG(I,K) - 180.0_r8
             IF ( ANG(I,K) .LT. -90.0_r8 ) ANG(I,K) = ANG(I,K) + 180.0_r8
             ANG(I,K) = ANG(I,K) * DEG_TO_RAD
             !
             UDS(I,K) = MAX(SQRT(U1(J,K)*U1(J,K) + V1(J,K)*V1(J,K)), minwnd)
             ! --- Test to see if we found Zb previously
             IF (IDXZB(I) .EQ. 0 ) THEN
                PE(I) = PE(I) + BNV2lm(I,K) *               &
                     ( con_g * ELVMAX(J) - phil(J,K) ) *          &
                     ( PHII(J,K+1) - PHII(J,K) ) / (con_g*con_g)
                ! --- KE
                ! --- Wind projected on the line perpendicular to mtn range, U(Zb(K)).
                ! --- kenetic energy is at the layer Zb
                ! --- THETA ranges from -+90deg |_ to the mtn "largest topo variations"
                UP(I)  =  UDS(I,K) * COS(ANG(I,K))
                EK(I)  = 0.5_r8 *  UP(I) * UP(I) 

                ! --- Dividing Stream lime  is found when PE =exceeds EK.
                IF ( PE(I) .GE.  EK(I) ) IDXZB(I) = K
                ! --- Then mtn blocked flow is between Zb=k(IDXZB(I)) and surface
                !
             ENDIF
          ENDDO
       ENDDO
       !
       !     print *,' in gwdps_lm.f 6  =',phiang,THETA(ipt(npt)),me
       !     print *,' in gwdps_lm.f 7  =',IDXZB(npt),PE(npt)
       !
       !
       DO I = 1, npt
          J    = ipt(i)
          ! --- Calc if N constant in layers (Zb guess) - a diagnostic only.
          ZBK(I) =  ELVMAX(J) - SQRT(UBAR(I)**2 + VBAR(I)**2)/BNV2bar(I)
       ENDDO
       !
       !
       ! --- The drag for mtn blocked flow
       ! 
       DO I = 1, npt
          J = ipt(i)
          ZLEN = 0.0_r8
          !      print *,' in gwdps_lm.f 9  =',i,j,IDXZB(i),me
          IF ( IDXZB(I) .GT. 0 ) THEN 
             DO K = IDXZB(I), 1, -1
                IF ( PHIL(J,IDXZB(I)) .GT.  PHIL(J,K) ) THEN
                   ZLEN = SQRT( ( PHIL(J,IDXZB(I)) - PHIL(J,K) ) / &
                        ( PHIL(J,K ) + con_g * hprime(J) ) )
                   ! --- lm eq 14:
                   R = (COS(ANG(I,K))**2 + GAMMA(J) * SIN(ANG(I,K))**2) / &
                        (gamma(J) * COS(ANG(I,K))**2 + SIN(ANG(I,K))**2)
                   ! --- (negitive of DB -- see sign at tendency)
                   DBTMP = 0.25_r8 *  CDmb *                                 &
                        MAX( 2.0_r8 - 1.0_r8 / R, 0.0_r8 ) * sigma(J) *            &
                        MAX(COS(ANG(I,K)), gamma(J)*SIN(ANG(I,K))) *   & 
                        ZLEN / hprime(J) 
                   DB(I,K) =  DBTMP * UDS(I,K)    
                   !
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       ! 
       !.............................
       !.............................
       ! end  mtn blocking section
       !
    ELSEIF ( NMTVR .NE. 14) THEN 
       ! ----  for mb not present and  gwd (nmtvr .ne .14) 
       ipt     = 0
       npt     = 0
       DO I = 1,IM
          IF ( hprime(i) .GT. hpmin )  THEN
             npt      = npt + 1
             ipt(npt) = i
             IF (ipr .EQ. i) npr = npt
          ENDIF
       ENDDO
       IF (npt .EQ. 0) RETURN     ! No gwd/mb calculation done!
       !
       DO i=1,npt
          IDXZB(i) = 0
       ENDDO
    ENDIF
    !
    !.............................
    !.............................
    !
    KMPBL  = km / 2 ! maximum pbl height : number of vertical levels / 2
    !
    !  Scale Gamma_Eff between IM=384*2 and 192*2 for T126/T170 and T62
    !
    IF (imx .GT. 0) THEN
       !
       ! Gamma_Eff is the effetive grid length, which was set to the length
       ! of grid boc in K96, but can be used pratically as a tunning coeffient
       !
       !Gamma_Eff = 1.0E-5_r8 * SQRT(FLOAT(IMX)/384.0_r8) !  this is inverse of Gamma_Eff!
       !Gamma_Eff = 1.0E-5_r8 * SQRT(FLOAT(IMX)/192.0_r8) !  this is inverse of Gamma_Eff!
       !Gamma_Eff = 0.5E-5_r8 * SQRT(FLOAT(IMX)/192.0_r8) !  this is inverse of Gamma_Eff!
       !Gamma_Eff = 1.0E-5_r8 * SQRT(FLOAT(IMX)/192.0_r8)/float(IMX/192)
       !Gamma_Eff = 1.0E-5_r8 / SQRT(FLOAT(IMX)/192.0_r8) !  this is inverse of Gamma_Eff!
       !Gamma_Eff = 0.5E-5_r8 / SQRT(REAL(IMX)/192.0_r8) !  this is inverse of Gamma_Eff!
       !Gamma_Eff = 0.5E-5_r8 / SQRT(REAL(IMX)/96.0_r8) !  this is inverse of Gamma_Eff!
       !
       ! 360    40000000
       !  x        y
       !
       ! x = 360/REAL(IMX)
       !
       ! y = (40000000 * (360/REAL(IMX))/360
       ! y = 40000000/REAL(IMX)
       !
       Gamma_Eff=1/(40000000.0_r8/REAL(IMX))
       !
       !Gamma_Eff = 2.0E-5_r8 * SQRT(FLOAT(IMX)/192.0_r8) !  this is inverse of Gamma_Eff!
       !Gamma_Eff = 2.5E-5_r8 * SQRT(FLOAT(IMX)/192.0_r8) !  this is inverse of Gamma_Eff!
    ENDIF
    IF (cdmbgwd(2) >= 0.0_r8) Gamma_Eff = Gamma_Eff * cdmbgwd(2)
    !
    DO K = 1,KM
       DO I =1,npt
          J         = ipt(i)
          VTJ(I,K)  = T1(J,K)  * (1.0_r8+FV*Q1(J,K))
          VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
          RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K) ! DENSITY TONS/M**3
          TAUP(I,K) = 0.0_r8
       ENDDO
    ENDDO
    DO K = 1,KMM1
       DO I =1,npt
          J         = ipt(i)
          TI        = 2.0_r8 / (T1(J,K)+T1(J,K+1))
          TEM       = TI  / (PRSL(J,K)-PRSL(J,K+1))
          RDZ       = con_g   / (phil(j,k+1) - phil(j,k))
          TEM1      = U1(J,K) - U1(J,K+1)
          TEM2      = V1(J,K) - V1(J,K+1)
          DW2       = TEM1*TEM1 + TEM2*TEM2
          SHR2      = MAX(DW2,DW2MIN) * RDZ * RDZ
          BVF2      = con_g*(GOCP+RDZ*(VTJ(I,K+1)-VTJ(I,K))) * TI
          ri_n(I,K) = MAX(BVF2/SHR2,RIMIN)   ! Richardson number
          !                                              Brunt-Vaisala Frequency
          !         TEM       = GR2 * (PRSL(J,K)+PRSL(J,K+1)) * TEM
          !         BNV2(I,K) = TEM * (VTK(I,K+1)-VTK(I,K))/(VTK(I,K+1)+VTK(I,K))
          BNV2(I,K) = (con_g+con_g) * RDZ * (VTK(I,K+1)-VTK(I,K)) &
               / (VTK(I,K+1)+VTK(I,K))
          bnv2(i,k) = MAX( bnv2(i,k), bnv2min )
       ENDDO
    ENDDO
    !      print *,' in gwdps_lm.f GWD:14  =',npt,kmm1,bnv2(npt,kmm1)
    !
    !     Apply 3 point smoothing on BNV2
    !
    !     do k=1,km
    !       do i=1,im
    !         vtk(i,k) = bnv2(i,k)
    !       enddo
    !     enddo
    !     do k=2,kmm1
    !       do i=1,im
    !         bnv2(i,k) = 0.25*(vtk(i,k-1)+vtk(i,k+1)) + 0.5*vtk(i,k)
    !       enddo
    !     enddo
    !
    !     Finding the first interface index above 50 hPa level
    !
    DO i=1,npt
       iwk(i) = 2
    ENDDO
    DO K=3,KMPBL
       DO I=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          IF (tem .LT. dpmin) iwk(i) = k
       ENDDO
    ENDDO
    !
    KBPS = 1
    KMPS = KM
    DO I=1,npt
       J         = ipt(i)
       kref(I)   = MAX(IWK(I), KPBL(J)+1 ) ! reference level 
       DELKS(I)  = 1.0_r8 / (PRSI(J,1) - PRSI(J,kref(I)))
       DELKS1(I) = 1.0_r8 / (PRSL(J,1) - PRSL(J,kref(I)))
       UBAR (I)  = 0.0_r8
       VBAR (I)  = 0.0_r8
       ROLL (I)  = 0.0_r8
       KBPS      = MAX(KBPS,  kref(I))
       KMPS      = MIN(KMPS,  kref(I))
       !
       BNV2bar(I) = (PRSL(J,1)-PRSL(J,2)) * DELKS1(I) * BNV2(I,1)
    ENDDO
    !      print *,' in gwdps_lm.f GWD:15  =',KBPS,KMPS
    KBPSP1 = KBPS + 1
    KBPSM1 = KBPS - 1
    DO K = 1,KBPS
       DO I = 1,npt
          IF (K .LT. kref(I)) THEN
             J          = ipt(i)
             RDELKS     = DEL(J,K) * DELKS(I)
             UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)   ! Mean U below kref
             VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)   ! Mean V below kref
             !
             ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)   ! Mean RO below kref
             RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS1(I)
             BNV2bar(I) = BNV2bar(I) + BNV2(I,K) * RDELKS
          ENDIF
       ENDDO
    ENDDO
    !      print *,' in gwdps_lm.f GWD:15B =',bnv2bar(npt)
    !
    !     FIGURE OUT LOW-LEVEL HORIZONTAL WIND DIRECTION AND FIND 'OA'
    !
    !             NWD  1   2   3   4   5   6   7   8
    !              WD  W   S  SW  NW   E   N  NE  SE
    !
    DO I = 1,npt
       J      = ipt(i)
       wdir   = ATAN2(UBAR(I),VBAR(I)) + pi
       idir   = MOD(NINT(fdir*wdir),mdir) + 1
       nwd    = nwdir(idir)
       OA(I)  = (1 - 2*INT( (NWD-1)/4 )) * OA4(J,MOD(NWD-1,4)+1)
       CLX(I) = CLX4(J,MOD(NWD-1,4)+1)
    ENDDO
    !
    !-----XN,YN            "LOW-LEVEL" WIND PROJECTIONS IN ZONAL
    !                                    & MERIDIONAL DIRECTIONS
    !-----ULOW             "LOW-LEVEL" WIND MAGNITUDE -        (= U)
    !-----BNV2             BNV2 = N**2
    !-----TAUB             BASE MOMENTUM FLUX
    !-----= -(RO * U**3/(N*XL)*GF(FR) FOR N**2 > 0
    !-----= 0.                        FOR N**2 < 0
    !-----FR               FROUDE    =   N*HPRIME / U
    !-----G                GMAX*FR**2/(FR**2+CG/OC)
    !
    !-----INITIALIZE SOME ARRAYS
    !
    DO I = 1,npt
       XN(I)     = 0.0_r8
       YN(I)     = 0.0_r8
       TAUB (I)  = 0.0_r8
       ULOW (I)  = 0.0_r8
       DTFAC(I)  = 1.0_r8
       ICRILV(I) = .FALSE. ! INITIALIZE CRITICAL LEVEL CONTROL VECTOR

       !
       !----COMPUTE THE "LOW LEVEL" WIND MAGNITUDE (M/S)
       !
       ULOW(I) = MAX(SQRT(UBAR(I)*UBAR(I) + VBAR(I)*VBAR(I)), 1.0_r8)
       ULOI(I) = 1.0_r8 / ULOW(I)
    ENDDO
    !
    DO  K = 1,KMM1
       DO  I = 1,npt
          J            = ipt(i)
          VELCO(I,K)   = 0.5_r8 * ((U1(J,K)+U1(J,K+1))*UBAR(I)   &
               +  (V1(J,K)+V1(J,K+1))*VBAR(I))
          VELCO(I,K)   = VELCO(I,K) * ULOI(I)
          !         IF ((VELCO(I,K).LT.VELEPS) .AND. (VELCO(I,K).GT.0.)) THEN
          !           VELCO(I,K) = VELEPS
          !         ENDIF
       ENDDO
    ENDDO
    !      
    !
    !   find the interface level of the projected wind where
    !   low levels & upper levels meet above pbl
    !
    DO i=1,npt
       kint(i) = km
    ENDDO
    DO k = 1,kmm1
       DO i = 1,npt
          IF (K .GT. kref(I)) THEN
             IF(velco(i,k) .LT. veleps .AND. kint(i) .EQ. km) THEN
                kint(i) = k+1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !  WARNING  KINT = KREF !!!!!!!!!
    DO i=1,npt
       kint(i) = kref(i)
    ENDDO
    !
    !Kim, Y.-j. and Doyle, J. D. (2005)
    !
    DO I = 1,npt
       J      = ipt(i)
       BNV    = SQRT( BNV2bar(I) )
       FR     = BNV     * ULOI(I) * MIN(HPRIME(J),hpmax)
       FR     = MIN(FR, FRMAX)
       XN(I)  = UBAR(I) * ULOI(I)
       YN(I)  = VBAR(I) * ULOI(I)
       !
       !     Compute the base level stress and store it in TAUB
       !     CALCULATE ENHANCEMENT FACTOR, NUMBER OF MOUNTAINS & ASPECT
       !     RATIO CONST. USE SIMPLIFIED RELATIONSHIP BETWEEN STANDARD
       !     DEVIATION & CRITICAL HGT
       !
       ! E is the 'anhancement factor' applied only at the reference level
       !   to representthe nonlinear enahncement of drag by
       !   resonant amplification ( Peltier and Clark, 1979 )
       !
       !      _    _ (Ce*Fro/Frc)
       !     |      |
       ! E = |OA + 2|
       !     |_    _|
       !
       EFACT    = (OA(I) + 2.0_r8) ** (CEOFRC*FR)
       !
       EFACT    = MIN( MAX(EFACT,EFMIN), EFMAX )
       !
       !
       ! Equation 7 [Kim, Y.-j. and Doyle, J. D. (2005)]
       !
       ! m if the 'number of mountains',  which estimate the bulk volume of sibgrid-scale
       ! orography associated with the nonlinearity of the flow Lx
       ! 
       !      _    _ (OA  + 1)
       !     |      |
       ! m = |1 + Lx|
       !     |_    _|
       !
       !
       COEFM    = (1.0_r8 + CLX(I)) ** (OA(I)+1.0_r8)
       !
       XLINV(I) = COEFM * Gamma_Eff
       !
       TEM      = FR * FR * OC(J)
       !
       !
       !                     
       !        Fr * Fr         
       ! G = ----------------
       !     Fr * Fr + Cg/OC 
       !
       !                     
       ! 1     Fr * Fr + Cg/OC  
       !--- = ----------------
       ! G       Fr * Fr      
       !
       !        Fr * Fr*OC + Cg  
       !       ---------------------
       ! 1           OC
       !--- = ----------------
       ! G       Fr * Fr      
       !
       ! 1       Fr * Fr*OC + Cg  
       !--- = ----------------
       ! G       Fr * Fr  * OC    
       !
       !        Fr * Fr * OC
       ! G = ------------------
       !       Fr * Fr*OC + Cg  
       !
       !
       ! G is an asymptotic function that provides a smooth transition between 2-D
       !   non-blocking and blocking cases as used by Pierrehumbert(1986) and include
       !   the influence of the vertical mountain aspect ratio through OC, empirically
       !   applying the original idea by Pierrehumbert (1986 [Ea 3.8])  
       !  
       GFOBNV   = GMAX  * TEM / ((TEM + CG)*BNV)  ! G/N0
       !
       TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I) * ULOW(I)  * GFOBNV  * EFACT         ! BASE FLUX Tau0
       !
       !tem      = min(HPRIME(I),hpmax)
       !TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * BNV * tem * tem
       !
       K        = MAX(1, kref(I)-1)
       TEM      = MAX(VELCO(I,K)*VELCO(I,K), 0.1_r8)
       SCOR(I)  = BNV2(I,K) / TEM  ! Scorer parameter below ref level
    ENDDO
    !                                                                       
    !----SET UP BOTTOM VALUES OF STRESS
    !
    DO K = 1, KBPS
       DO I = 1,npt
          IF (K .LE. kref(I)) TAUP(I,K) = TAUB(I)
       ENDDO
    ENDDO
    !
    !   Now compute vertical structure of the stress.
    !
    DO K = KMPS, KMM1                   ! Vertical Level K Loop!
       KP1 = K + 1
       DO I = 1, npt
          !
          !-----UNSTABLE LAYER IF RI < RIC
          !-----UNSTABLE LAYER IF UPPER AIR VEL COMP ALONG SURF VEL <=0 (CRIT LAY)
          !---- AT (U-C)=0. CRIT LAYER EXISTS AND BIT VECTOR SHOULD BE SET (.LE.)
          !
          IF (K .GE. kref(I)) THEN
             ICRILV(I) = ICRILV(I) .OR. ( ri_n(I,K) .LT. RIC)   &
                  .OR. (VELCO(I,K) .LE. 0.0_r8)
          ENDIF
       ENDDO
       !
       DO I = 1,npt
          IF (K .GE. kref(I))   THEN
             IF (.NOT.ICRILV(I) .AND. TAUP(I,K) .GT. 0.0_r8 ) THEN
                TEMV = 1.0_r8 / MAX(VELCO(I,K), 0.01_r8)
                !             IF (OA(I) .GT. 0. .AND.  PRSI(ipt(i),KP1).GT.RLOLEV) THEN
                IF (OA(I).GT.0.0_r8 .AND. kp1 .LT. kint(i)) THEN
                   SCORK   = BNV2(I,K) * TEMV * TEMV
                   RSCOR   = MIN(1.0_r8, SCORK / SCOR(I))
                   SCOR(I) = SCORK
                ELSE 
                   RSCOR   = 1.0_r8
                ENDIF
                !
                BRVF = SQRT(BNV2(I,K))        ! Brunt-Vaisala Frequency
                !             TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*VELCO(I,K)*0.5
                TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*0.5_r8   &
                     * MAX(VELCO(I,K),0.01_r8)
                HD   = SQRT(TAUP(I,K) / TEM1)
                FRO  = BRVF * HD * TEMV
                !
                !    RIM is the  MINIMUM-RICHARDSON NUMBER BY SHUTTS (1985)
                !
                TEM2   = SQRT(ri_n(I,K))
                TEM    = 1.0_r8 + TEM2 * FRO
                RIM    = ri_n(I,K) * (1.0_r8-FRO) / (TEM * TEM)
                !
                !    CHECK STABILITY TO EMPLOY THE 'SATURATION HYPOTHESIS'
                !    OF LINDZEN (1981) EXCEPT AT TROPOSPHERIC DOWNSTREAM REGIONS
                !
                !                                       ----------------------
                IF (RIM .LE. RIC .AND. &
                     !    &           (OA(I) .LE. 0. .OR.  PRSI(ipt(I),KP1).LE.RLOLEV )) THEN
                     (OA(I) .LE. 0.0_r8 .OR.  kp1 .GE. kint(i) )) THEN
                   TEMC = 2.0_r8 + 1.0_r8 / TEM2
                   HD   = VELCO(I,K) * (2.0_r8*SQRT(TEMC)-TEMC) / BRVF
                   TAUP(I,KP1) = TEM1 * HD * HD
                ELSE 
                   TAUP(I,KP1) = TAUP(I,K) * RSCOR
                ENDIF
                taup(i,kp1) = MIN(taup(i,kp1), taup(i,k))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !     DO I=1,IM
    !       taup(i,km+1) = taup(i,km)
    !     ENDDO
    !
    IF(LCAP .LE. KM) THEN
       DO KLCAP = LCAPP1, KM+1
          DO I = 1,npt
             SIRA          = PRSI(ipt(I),KLCAP) / PRSI(ipt(I),LCAP)
             TAUP(I,KLCAP) = SIRA * TAUP(I,LCAP)
          ENDDO
       ENDDO
    ENDIF
    !
    !     Calculate - (g/p*)*d(tau)/d(sigma) and Decel terms DTAUX, DTAUY
    !
    DO K = 1,KM
       DO I = 1,npt
          TAUD(I,K) = con_g * (TAUP(I,K+1) - TAUP(I,K)) / DEL(ipt(I),K)
       ENDDO
    ENDDO
    !
    !------LIMIT DE-ACCELERATION (MOMENTUM DEPOSITION ) AT TOP TO 1/2 VALUE
    !------THE IDEA IS SOME STUFF MUST GO OUT THE 'TOP'
    !
    DO KLCAP = LCAP, KM
       DO I = 1,npt
          TAUD(I,KLCAP) = TAUD(I,KLCAP) * FACTOP
       ENDDO
    ENDDO
    !
    !------IF THE GRAVITY WAVE DRAG WOULD FORCE A CRITICAL LINE IN THE
    !------LAYERS BELOW SIGMA=RLOLEV DURING THE NEXT DELTIM TIMESTEP,
    !------THEN ONLY APPLY DRAG UNTIL THAT CRITICAL LINE IS REACHED.
    !
    DO K = 1,KMM1
       DO I = 1,npt
          IF (K .GT. kref(I) .AND. PRSI(ipt(i),K) .GE. RLOLEV) THEN
             IF(TAUD(I,K).NE.0.0_r8) THEN
                TEM = DELTIM * TAUD(I,K)
                DTFAC(I) = MIN(DTFAC(I),ABS(VELCO(I,K)/TEM))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    DO K = 1,KM
       DO I = 1,npt
          J          = ipt(i)
          TAUD(I,K)  = TAUD(I,K) * DTFAC(I)
          
          IF(ABS(TAUD(I,K)) > CRITAC) TAUD(I,K)=sign(1.0_r8,TAUD(I,K))*CRITAC
          
          DTAUX      = TAUD(I,K) * XN(I)
          DTAUY      = TAUD(I,K) * YN(I)
          ! ---  lm mb (*j*)  changes overwrite GWD
          IF ( K .LT. IDXZB(I) .AND. IDXZB(I) .NE. 0 ) THEN
             DBIM = DB(I,K) / (1.0_r8 + DB(I,K)*DELTIM)
             dvdt(J,K)  =  dvdt(J,K) - DBIM * V1(J,K) 
             dudt(J,K)  =  dudt(J,K) - DBIM * U1(J,K) 
             !          if ( ABS(DBIM * U1(J,K)) .gt. .01 ) 
             !    & print *,' in gwdps_lmi.f KDT=',KDT,I,K,DB(I,K),
             !    &                      dbim,idxzb(I),U1(J,K),V1(J,K),me
             DUSFC(J)   = DUSFC(J) - DBIM * U1(J,K) * DEL(J,K)
             DVSFC(J)   = DVSFC(J) - DBIM * V1(J,K) * DEL(J,K)
          ELSE
             !
             dvdt(J,K)     = DTAUY     + dvdt(J,K)
             dudt(J,K)     = DTAUX     + dudt(J,K)
             DUSFC(J)      = DUSFC(J)  + DTAUX * DEL(J,K)
             DVSFC(J)      = DVSFC(J)  + DTAUY * DEL(J,K)
          ENDIF
       ENDDO
    ENDDO
    TEM    = -1.0_r8/con_g
    DO I = 1,npt
       J          = ipt(i)
       !       TEM    = (-1.E3/G)
       DUSFC(J) = TEM * DUSFC(J)
       DVSFC(J) = TEM * DVSFC(J)
    ENDDO
    !                                                                       
    !    MONITOR FOR EXCESSIVE GRAVITY WAVE DRAG TENDENCIES IF NCNT>0
    !
    !     IF(NCNT.GT.0) THEN
    !        IF(LAT.GE.38.AND.LAT.LE.42) THEN
    !CMIC$ GUARD 37
    !           DO 92 I = 1,IM
    !              IF(IKOUNT.GT.NCNT) GO TO 92
    !              IF(I.LT.319.OR.I.GT.320) GO TO 92
    !              DO 91 K = 1,KM
    !                 IF(ABS(TAUD(I,K)) .GT. CRITAC) THEN
    !                    IF(I.LE.IM) THEN
    !                       IKOUNT = IKOUNT+1
    !                       PRINT 123,I,LAT,KDT
    !                       PRINT 124,TAUB(I),BNV(I),ULOW(I),
    !    1                  GF(I),FR(I),ROLL(I),HPRIME(I),XN(I),YN(I)
    !                       PRINT 124,(TAUD(I,KK),KK = 1,KM)
    !                       PRINT 124,(TAUP(I,KK),KK = 1,KM+1)
    !                       PRINT 124,(ri_n(I,KK),KK = 1,KM)
    !                       DO 93 KK = 1,KMM1
    !                          VELKO(KK) =
    !    1                  0.5*((U1(I,KK)+U1(I,KK+1))*UBAR(I)+
    !    2                  (V1(I,KK)+V1(I,KK+1))*VBAR(I))*ULOI(I)
    !93                     CONTINUE
    !                       PRINT 124,(VELKO(KK),KK = 1,KMM1)
    !                       PRINT 124,(A    (I,KK),KK = 1,KM)
    !                       PRINT 124,(DTAUY(I,KK),KK = 1,KM)
    !                       PRINT 124,(B    (I,KK),KK = 1,KM)
    !                       PRINT 124,(DTAUX(I,KK),KK = 1,KM)
    !                       GO TO 92
    !                    ENDIF
    !                 ENDIF
    !91            CONTINUE
    !92         CONTINUE
    !CMIC$ END GUARD 37
    !123        FORMAT('  *** MIGWD PRINT *** I=',I3,' LAT=',I3,' KDT=',I3)
    !124        FORMAT(2X,  10E13.6)
    !        ENDIF
    !     ENDIF
    !
    !      print *,' in gwdps_lm.f 18  =',A(ipt(1),idxzb(1))
    !    &,                          B(ipt(1),idxzb(1)),me
    RETURN
  END SUBROUTINE GWDPS




  SUBROUTINE gwdc( im       , &
       km       , &
       u1       , &
       v1       , &
       t1       , &
       q1       , &
       pmid1    , &
       pint1    , &
       dpmid1   , &
       qmax     , &
       cumchr1  , &
       ktop     , &
       kbot     , &
       kuo      , &
       dlength  , &
       fu1      , &
       fv1      )

    !***********************************************************************
    !        ORIGINAL CODE FOR PARAMETERIZATION OF CONVECTIVELY FORCED
    !        GRAVITY WAVE DRAG FROM YONSEI UNIVERSITY, KOREA
    !        BASED ON THE THEORY GIVEN BY CHUN AND BAIK (JAS, 1998)
    !        MODIFIED FOR IMPLEMENTATION INTO THE GFS/CFS BY
    !        AKE JOHANSSON  --- AUG 2005
    !***********************************************************************

    IMPLICIT NONE

    !---------------------------- Arguments --------------------------------
    !
    !  Input variables
    !
    !  u        : midpoint zonal wind
    !  v        : midpoint meridional wind
    !  t        : midpoint temperatures
    !  pmid     : midpoint pressures
    !  pint     : interface pressures
    !  dpmid    : midpoint delta p ( pi(k)-pi(k-1) )
    !  qmax     : deep convective heating
    !  kcldtop  : Vertical level index for cloud top    ( mid level ) 
    !  kcldbot  : Vertical level index for cloud bottom ( mid level )
    !  kuo      : (0,1) dependent on whether convection occur or not
    !
    !  Output variables
    !
    !  fu1     : zonal wind tendency
    !  fv1     : meridional wind tendency
    !
    !-----------------------------------------------------------------------

    INTEGER :: im
    INTEGER :: km
    INTEGER :: ktop(im)
    INTEGER :: kbot(im)
    INTEGER :: kuo(im)
    INTEGER :: kcldtop(im)
    INTEGER :: kcldbot(im)

    REAL(kind=r8) :: dlength(im)
    REAL(kind=r8) :: qmax(im)
    REAL(kind=r8) :: cumchr1(im,km)
    REAL(kind=r8) :: cumchr(im,km)
    REAL(kind=r8) :: u1(im,km)
    REAL(kind=r8) :: v1(im,km)
    REAL(kind=r8) :: t1(im,km)
    REAL(kind=r8) :: q1(im,km)

    REAL(kind=r8) :: pmid1(im,km)
    REAL(kind=r8) :: dpmid1(im,km)
    REAL(kind=r8) :: pint1(im,km+1)
    REAL(kind=r8) :: fu1(im,km)
    REAL(kind=r8) :: fv1(im,km)
    REAL(kind=r8) :: u(im,km)
    REAL(kind=r8) :: v(im,km)
    REAL(kind=r8) :: t(im,km)
    REAL(kind=r8) :: spfh(im,km)

    REAL(kind=r8) :: pmid(im,km)
    REAL(kind=r8) :: dpmid(im,km)
    REAL(kind=r8) :: pint(im,km+1)


    !------------------------- Local workspace -----------------------------
    !
    !  i, k     : Loop index
    !  ii,kk    : Loop index
    !  cldbar   : Deep convective cloud coverage at the cloud top.
    !  ugwdc    : Zonal wind after GWDC paramterization
    !  vgwdc    : Meridional wind after GWDC parameterization
    !  plnmid   : Log(pmid) ( mid level )
    !  plnint   : Log(pint) ( interface level )
    !  dpint    : Delta pmid ( interface level )
    !  tauct    : Wave stress at the cloud top calculated using basic-wind
    !             parallel to the wind vector at the cloud top ( mid level )
    !  tauctx   : Wave stress at the cloud top projected in the east
    !  taucty   : Wave stress at the cloud top projected in the north
    !  qmax     : Maximum deep convective heating rate ( K s-1 ) in a  
    !             horizontal grid point calculated from cumulus para-
    !             meterization. ( mid level )
    !  wtgwc    : Wind tendency in direction to the wind vector at the cloud top level
    !             due to convectively generated gravity waves ( mid level )
    !  utgwc    : Zonal wind tendency due to convectively generated 
    !             gravity waves ( mid level )
    !  vtgwc    : Meridional wind tendency due to convectively generated
    !             gravity waves ( mid level )
    !  taugwci  : Profile of wave stress calculated using basic-wind
    !             parallel to the wind vector at the cloud top 
    !  taugwcxi : Profile of zonal component of gravity wave stress
    !  taugwcyi : Profile of meridional component of gravity wave stress 
    !
    !  taugwci, taugwcxi, and taugwcyi are defined at the interface level
    !
    !  bruni    : Brunt-Vaisala frequency ( interface level )
    !  brunm    : Brunt-Vaisala frequency ( mid level )
    !  rhoi     : Air density ( interface level )
    !  rhom     : Air density ( mid level )
    !  ti       : Temperature ( interface level )
    !  basicum  : Basic-wind profile. Basic-wind is parallel to the wind
    !             vector at the cloud top level. (mid level) 
    !  basicui  : Basic-wind profile. Basic-wind is parallel to the wind
    !             vector at the cloud top level. ( interface level )
    !  riloc    : Local Richardson number ( interface level )
    !  rimin    : Minimum Richardson number including both the basic-state
    !             and gravity wave effects ( interface level )
    !  gwdcloc  : Horizontal location where the GWDC scheme is activated.
    !  break    : Horizontal location where wave breaking is occurred.
    !  critic   : Horizontal location where critical level filtering is
    !             occurred.
    !  dogwdc   : Logical flag whether the GWDC parameterization is           
    !             calculated at a grid point or not.
    !  
    !  dogwdc is used in order to lessen CPU time for GWDC calculation.
    !
    !-----------------------------------------------------------------------

    INTEGER       :: i
    INTEGER       :: k
    INTEGER       :: k1
    INTEGER       :: kk
    INTEGER       :: kb

    REAL(kind=r8) :: cldbar(im)

    REAL(kind=r8) :: ugwdc(im,km)
    REAL(kind=r8) :: vgwdc(im,km)

    REAL(kind=r8) :: plnmid(im,km)
    REAL(kind=r8) :: plnint(im,km+1)
    REAL(kind=r8) :: dpint(im,km+1)

    REAL(kind=r8) :: tauct(im)
    REAL(kind=r8) :: tauctx(im)
    REAL(kind=r8) :: taucty(im)

    REAL(kind=r8) :: wtgwc(im,km)
    REAL(kind=r8) :: utgwc(im,km)
    REAL(kind=r8) :: vtgwc(im,km)

    REAL(kind=r8) :: taugwci(im,km+1)
    REAL(kind=r8) :: taugwcxi(im,km+1)
    REAL(kind=r8) :: taugwcyi(im,km+1)

    REAL(kind=r8) :: bruni(im,km+1)
    REAL(kind=r8) :: rhoi(im,km+1)
    REAL(kind=r8) :: ti(im,km+1)

    REAL(kind=r8) :: brunm(im,km)
    REAL(kind=r8) :: rhom(im,km)
    REAL(kind=r8) :: brunm1(im,km)
    REAL(kind=r8) :: rhom1(im,km)

    REAL(kind=r8) :: basicum(im,km)
    REAL(kind=r8) :: basicui(im,km+1)

    REAL(kind=r8) :: riloc(km+1)
    REAL(kind=r8) :: rimin(km+1)

    REAL(kind=r8) :: gwdcloc(im)
    REAL(kind=r8) :: break(im)
    REAL(kind=r8) :: critic(im)
    REAL(kind=r8) :: tem1
    REAL(kind=r8) :: tem2
    REAL(kind=r8) :: qtem

    LOGICAL       :: dogwdc(im)

    !-----------------------------------------------------------------------
    !
    !  ucltop    : Zonal wind at the cloud top ( mid level )
    !  vcltop    : Meridional wind at the cloud top ( mid level )
    !  windcltop : Wind speed at the cloud top ( mid level )
    !  shear     : Vertical shear of basic wind 
    !  cosphi    : Cosine of angle of wind vector at the cloud top
    !  sinphi    : Sine   of angle of wind vector at the cloud top
    !  c1        : Tunable parameter
    !  c2        : Tunable parameter
    !  dlength   : Grid spacing in the direction of basic wind at the cloud top
    !  nonlinct  : Nonlinear parameter at the cloud top
    !  nonlin    : Nonlinear parameter above the cloud top
    !  nonlins   : Saturation nonlinear parameter
    !  taus      : Saturation gravity wave drag
    !  n2        : Square of Brunt-Vaisala frequency
    !  dtdp      : dT/dp
    !  xstress   : Vertically integrated zonal momentum change due to GWDC
    !  ystress   : Vertically integrated meridional momentum change due to GWDC
    !  crit1     : Variable 1 for checking critical level
    !  crit2     : Variable 2 for checking critical level
    !  sum1      : Temporary variable
    !
    !-----------------------------------------------------------------------

    REAL(kind=r8) :: ucltop
    REAL(kind=r8) :: vcltop
    REAL(kind=r8) :: windcltop
    REAL(kind=r8) :: shear
    !REAL(kind=r8) :: kcldtopi
    REAL(kind=r8) :: cosphi
    REAL(kind=r8) :: sinphi
    REAL(kind=r8) :: angle
    REAL(kind=r8) :: nonlinct
    REAL(kind=r8) :: nonlin
    REAL(kind=r8) :: nonlins
    REAL(kind=r8) :: taus 

    !-----------------------------------------------------------------------
    REAL(kind=r8), PARAMETER :: c1=1.41_r8
    REAL(kind=r8), PARAMETER :: c2=-0.38_r8
    REAL(kind=r8), PARAMETER :: ricrit=0.25_r8
    REAL(kind=r8), PARAMETER :: n2min=1.e-32_r8
    REAL(kind=r8), PARAMETER :: zero=0.0_r8
    REAL(kind=r8), PARAMETER :: one=1.0_r8
    REAL(kind=r8), PARAMETER :: taumin=1.0e-20_r8
    REAL(kind=r8), PARAMETER :: tauctmax=-20.0_r8
    REAL(kind=r8), PARAMETER :: qmin=1.0e-10_r8
    REAL(kind=r8), PARAMETER :: shmin=1.0e-20_r8
    REAL(kind=r8), PARAMETER :: rimax=1.0e+20_r8
    REAL(kind=r8), PARAMETER :: rimaxm=0.99e+20_r8
    REAL(kind=r8), PARAMETER :: rimaxp=1.01e+20_r8
    REAL(kind=r8), PARAMETER :: rilarge=0.9e+20_r8
    REAL(kind=r8), PARAMETER :: riminx=-1.0e+20_r8
    REAL(kind=r8), PARAMETER :: riminm=-1.01e+20_r8
    REAL(kind=r8), PARAMETER :: riminp=-0.99e+20_r8
    REAL(kind=r8), PARAMETER :: rismall=-0.9e+20_r8

    REAL(kind=r8) :: n2
    REAL(kind=r8) :: dtdp
    REAL(kind=r8) :: xstress
    REAL(kind=r8) :: ystress
    REAL(kind=r8) :: crit1
    REAL(kind=r8) :: crit2
    REAL(kind=r8) :: pi
    REAL(kind=r8) :: p1
    REAL(kind=r8) :: p2
    !-----------------------------------------------------------------------
    !        Write out incoming variables
    !-----------------------------------------------------------------------

    !     fhourpr = zero
    !     if (lprnt) then
    !       if (fhour.ge.fhourpr) then
    !         print *,' '
    !         write(*,*) 'Inside GWDC raw input start print at fhour = ',
    !    &               fhour
    !         write(*,*) 'im  IM  KM  ',im,im,km
    !         write(*,*) 'KBOT KTOP QMAX DLENGTH KUO  ',
    !    +     kbot(ipr),ktop(ipr),qmax(ipr),dlength(ipr),kuo(ipr)
    !         write(*,*) 'g  cp  rd  ',g,cp,rd

    !-------- Pressure levels ----------
    !         write(*,9100)
    !         ilev=km+1
    !         write(*,9110) ilev,(10.*pint1(ipr,ilev))
    !         do ilev=km,1,-1
    !           write(*,9120) ilev,(10.*pmid1(ipr,ilev)),
    !    &                         (10.*dpmid1(ipr,ilev))
    !           write(*,9110) ilev,(10.*pint1(ipr,ilev))
    !         enddo

    !-------- U1 V1 T1 ----------
    !         write(*,9130)
    !         do ilev=km,1,-1
    !           write(*,9140) ilev,U1(ipr,ilev),V1(ipr,ilev),T1(ipr,ilev)
    !         enddo

    !         print *,' '
    !         print *,' Inside GWDC raw input end print'
    !       endif
    !     endif

    !9100 FORMAT(//,14x,'PRESSURE LEVELS',//,&
    !         ' ILEV',6x,'PINT1',7x,'PMID1',6x,'DPMID1',/)
    !9110 FORMAT(i4,2x,f10.3)
    !9120 FORMAT(i4,12x,2(2x,f10.3))
    !9130 FORMAT(//,' ILEV',7x,'U1',10x,'V1',10x,'T1',/)
    !9140 FORMAT(i4,3(2x,f10.3))

    !-----------------------------------------------------------------------
    !        Create local arrays with reversed vertical indices
    !-----------------------------------------------------------------------

    DO k=1,km
       k1 = km - k + 1
       DO i=1,im
          u(i,k)      = u1(i,k1)
          v(i,k)      = v1(i,k1)
          t(i,k)      = t1(i,k1)
          spfh(i,k)   = MAX(q1(i,k1),qmin)
          pmid(i,k)   = pmid1(i,k1)
          dpmid(i,k)  = dpmid1(i,k1)
          cumchr(i,k) = cumchr1(i,k1)
       ENDDO
    ENDDO

    DO k=1,km+1
       k1 = km - k + 2
       DO i=1,im
          pint(i,k) = pint1(i,k1)
       ENDDO
    ENDDO

    DO i = 1, im
       kcldtop(i) = km - ktop(i) + 1
       kcldbot(i) = km - kbot(i) + 1
    ENDDO

    !      if (lprnt) then
    !        if (fhour.ge.fhourpr) then
    !          write(*,9200)
    !          do i=1,im
    !            write(*,9201) kuo(i),kcldbot(i),kcldtop(i)
    !          enddo
    !        endif
    !      endif

    !9200 FORMAT(//,'  Inside GWDC local variables start print',//,&
    !         2x,'KUO',2x,'KCLDBOT',2x,'KCLDTOP',//)
    !9201 FORMAT(i4,2x,i5,4x,i5)

    !***********************************************************************
    !
    !  Begin GWDC
    !
    !***********************************************************************

    pi     = 2.0_r8*ASIN(1.0_r8)

    !-----------------------------------------------------------------------
    !
    !  Initialize local variables
    !
    !-----------------------------------------------------------------------
    !                              PRESSURE VARIABLES
    !
    !  Interface 1 ======== pint(1)           *********
    !  Mid-Level 1 --------          pmid(1)            dpmid(1)
    !            2 ======== pint(2)           dpint(2)
    !            2 --------          pmid(2)            dpmid(2)
    !            3 ======== pint(3)           dpint(3)
    !            3 --------          pmid(3)            dpmid(3)
    !            4 ======== pint(4)           dpint(4)
    !            4 --------          pmid(4)            dpmid(4)
    !              ........
    !           17 ======== pint(17)          dpint(17) 
    !           17 --------          pmid(17)           dpmid(17)
    !           18 ======== pint(18)          dpint(18)
    !           18 --------          pmid(18)           dpmid(18)
    !           19 ======== pint(19)          *********
    !
    !-----------------------------------------------------------------------

    DO k = 1, km+1
       DO i = 1, im
          plnint(i,k)   = LOG(MAX(pint(i,k),0.00000001_r8))
          taugwci(i,k)  = zero
          taugwcxi(i,k) = zero
          taugwcyi(i,k) = zero
          bruni(i,k)    = zero
          rhoi(i,k)     = zero
          ti(i,k)       = zero
          basicui(i,k)  = zero
          riloc(k)      = zero
          rimin(k)      = zero
       ENDDO
    ENDDO

    DO k = 1, km
       DO i = 1, im
          plnmid(i,k)  = LOG(pmid(i,k))
          wtgwc(i,k)   = zero
          utgwc(i,k)   = zero
          vtgwc(i,k)   = zero
          ugwdc(i,k)   = zero
          vgwdc(i,k)   = zero
          brunm(i,k)   = zero
          rhom(i,k)    = zero
          basicum(i,k) = zero
       ENDDO
    ENDDO

    DO k = 2, km
       DO i = 1, im
          dpint(i,k) = pmid(i,k) - pmid(i,k-1)
       ENDDO
    ENDDO

    DO i = 1, im
       dpint(i,1)    = zero
       dpint(i,km+1) = zero
       tauct(i)      = zero
       tauctx(i)     = zero
       taucty(i)     = zero
       gwdcloc(i)    = zero
       break(i)      = zero
       critic(i)     = zero
    ENDDO

    !-----------------------------------------------------------------------
    !                              THERMAL VARIABLES
    !
    !  Interface 1 ========       TI(1)           RHOI(1)            BRUNI(1)
    !            1 -------- T(1)         RHOM(1)           BRUNM(1)
    !            2 ========       TI(2)           RHOI(2)            BRUNI(2)
    !            2 -------- T(2)         RHOM(2)           BRUNM(2)
    !            3 ========       TI(3)           RHOI(3)            BRUNI(3)
    !            3 -------- T(3)         RHOM(3)           BRUNM(3)
    !            4 ========       TI(4)           RHOI(4)            BRUNI(4)
    !            4 -------- T(4)         RHOM(4)           BRUNM(4)
    !              ........
    !           17 ========
    !           17 -------- T(17)        RHOM(17)          BRUNM(17)
    !           18 ========       TI(18)          RHOI(18)           BRUNI(18)
    !           18 -------- T(18)        RHOM(18)          BRUNM(18)
    !           19 ========       TI(19)          RHOI(19)           BRUNI(19)
    !
    !-----------------------------------------------------------------------

    DO k = 1, km
       DO i = 1, im
          rhom(i,k) = pmid(i,k) / (con_rd*t(i,k)*(1.0_r8+con_fvirt*spfh(i,k)))
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !
    !  Top interface temperature is calculated assuming an isothermal 
    !  atmosphere above the top mid level.
    !
    !-----------------------------------------------------------------------

    DO i = 1, im
       ti(i,1)    = t(i,1)
       rhoi(i,1)  = pint(i,1)/(con_rd*ti(i,1))
       bruni(i,1) = SQRT ( con_g*con_g / (con_cp*ti(i,1)) )
    ENDDO

    !-----------------------------------------------------------------------
    !
    !  Calculate interface level temperature, density, and Brunt-Vaisala
    !  frequencies based on linear interpolation of Temp in ln(Pressure)
    !
    !-----------------------------------------------------------------------

    DO k = 2, km
       DO i = 1, im
          tem1 = (plnmid(i,k)-plnint(i,k)) / (plnmid(i,k)-plnmid(i,k-1))
          tem2 = one - tem1
          ti(i,k)    = t(i,k-1)    * tem1 + t(i,k)    * tem2
          qtem       = spfh(i,k-1) * tem1 + spfh(i,k) * tem2
          rhoi(i,k)  = pint(i,k) / ( con_rd * ti(i,k)*(1.0_r8+con_fvirt*qtem) )
          dtdp       = (t(i,k)-t(i,k-1)) / (pmid(i,k)-pmid(i,k-1))
          n2         = con_g*con_g/ti(i,k)*( 1.0_r8/con_cp - rhoi(i,k)*dtdp ) 
          bruni(i,k) = SQRT (MAX (n2min, n2))
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !
    !  Bottom interface temperature is calculated assuming an isothermal
    !  atmosphere below the bottom mid level
    !
    !-----------------------------------------------------------------------

    DO i = 1, im
       ti(i,km+1)    = t(i,km)
       rhoi(i,km+1)  = pint(i,km+1)/(con_rd*ti(i,km+1)*(1.0_r8+con_fvirt*spfh(i,km)))
       bruni(i,km+1) = SQRT ( con_g*con_g / (con_cp*ti(i,km+1)) )
    ENDDO

    !-----------------------------------------------------------------------
    !
    !  Determine the mid-level Brunt-Vaisala frequencies.
    !             based on interpolated interface Temperatures [ ti ]
    !
    !-----------------------------------------------------------------------

    DO k = 1, km
       DO i = 1, im
          dtdp       = (ti(i,k+1)-ti(i,k)) / (pint(i,k+1)-pint(i,k))
          n2         = con_g*con_g/t(i,k)*( 1.0_r8/con_cp - rhom(i,k)*dtdp ) 
          brunm(i,k) = SQRT (MAX (n2min, n2))
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !        PRINTOUT
    !-----------------------------------------------------------------------

    !     if (lprnt) then
    !       if (fhour.ge.fhourpr) then

    !-------- Pressure levels ----------
    !         write(*,9101)
    !         do ilev=1,km
    !           write(*,9111) ilev,(0.01*pint(ipr,ilev)),
    !    &                         (0.01*dpint(ipr,ilev)),plnint(ipr,ilev)
    !           write(*,9121) ilev,(0.01*pmid(ipr,ilev)),
    !    &                         (0.01*dpmid(ipr,ilev)),plnmid(ipr,ilev)
    !         enddo
    !         ilev=km+1
    !         write(*,9111) ilev,(0.01*pint(ipr,ilev)),
    !    &                       (0.01*dpint(ipr,ilev)),plnint(ipr,ilev)

    !                2
    !-------- U V T N  ----------
    !         write(*,9102)
    !         do ilev=1,km
    !           write(*,9112) ilev,ti(ipr,ilev),(100.*bruni(ipr,ilev))
    !           write(*,9122) ilev,u(ipr,ilev),v(ipr,ilev),
    !    +                    t(ipr,ilev),(100.*brunm(ipr,ilev))
    !         enddo
    !         ilev=km+1
    !         write(*,9112) ilev,ti(ipr,ilev),(100.*bruni(ipr,ilev))

    !       endif
    !     endif

    !9101 FORMAT(//,14x,'PRESSURE LEVELS',//,  &
    !         ' ILEV',4x,'PINT',4x,'PMID',4x,'DPINT',3x,'DPMID',5x,'LNP',/)
    !9111 FORMAT(i4,1x,f8.2,9x,f8.2,9x,f8.2)
    !9121 FORMAT(i4,9x,f8.2,9x,f8.2,1x,f8.2)
    !9102 FORMAT(//' ILEV',5x,'U',7x,'V',5x,'TI',7x,'T',&
    !         5x,'BRUNI',3x,'BRUNM',//)
    !9112 FORMAT(i4,16x,f8.2,8x,f8.3)
    !9122 FORMAT(i4,2f8.2,8x,f8.2,8x,f8.3)

    !-----------------------------------------------------------------------
    !
    !  Set switch for no convection present
    !
    !-----------------------------------------------------------------------

    DO i = 1, im
       dogwdc(i) =.TRUE.
       IF (kuo(i) == 0 .OR. qmax(i) <= zero) dogwdc(i) =.FALSE.
    ENDDO

    !***********************************************************************
    !
    !        Big loop over grid points                    ONLY done if KUO=1
    !
    !***********************************************************************

    DO i = 1, im

       IF ( dogwdc(i) ) THEN                 !  For fast GWDC calculation

          kk        = kcldtop(i)
          kb        = kcldbot(i)
          cldbar(i) = 0.1_r8

          !-----------------------------------------------------------------------
          !
          !  Determine cloud top wind component, direction, and speed.
          !  Here, ucltop, vcltop, and windcltop are wind components and 
          !  wind speed at mid-level cloud top index
          !
          !-----------------------------------------------------------------------

          ucltop    = u(i,kcldtop(i))
          vcltop    = v(i,kcldtop(i))
          windcltop = SQRT( ucltop*ucltop + vcltop*vcltop )
          cosphi    = ucltop/windcltop
          sinphi    = vcltop/windcltop
          angle     = ACOS(cosphi)*180.0_r8/pi

          !-----------------------------------------------------------------------
          !
          !  Calculate basic state wind projected in the direction of the cloud 
          !  top wind.
          !  Input u(i,k) and v(i,k) is defined at mid level
          !
          !-----------------------------------------------------------------------

          DO k=1,km
             basicum(i,k) = u(i,k)*cosphi + v(i,k)*sinphi
          ENDDO

          !-----------------------------------------------------------------------
          !
          !  Basic state wind at interface level is also calculated
          !  based on linear interpolation in ln(Pressure)
          !
          !  In the top and bottom boundaries, basic-state wind at interface level
          !  is assumed to be vertically uniform.
          !
          !-----------------------------------------------------------------------

          basicui(i,1)   = basicum(i,1)
          DO k=2,km
             tem1 = (plnmid(i,k)-plnint(i,k)) / (plnmid(i,k)-plnmid(i,k-1))
             tem2 = one - tem1
             basicui(i,k) = basicum(i,k)*tem2 + basicum(i,k-1)*tem2
          ENDDO
          basicui(i,km+1) = basicum(i,km)

          !-----------------------------------------------------------------------
          !
          !  Calculate local richardson number 
          !
          !  basicum   : U at mid level
          !  basicui   : UI at interface level
          !
          !  Interface 1 ========       UI(1)            rhoi(1)  bruni(1)  riloc(1)
          !  Mid-level 1 -------- U(1)
          !            2 ========       UI(2)  dpint(2)  rhoi(2)  bruni(2)  riloc(2)
          !            2 -------- U(2)
          !            3 ========       UI(3)  dpint(3)  rhoi(3)  bruni(3)  riloc(3)
          !            3 -------- U(3)
          !            4 ========       UI(4)  dpint(4)  rhoi(4)  bruni(4)  riloc(4)
          !            4 -------- U(4)
          !              ........
          !           17 ========       UI(17) dpint(17) rhoi(17) bruni(17) riloc(17)
          !           17 -------- U(17)
          !           18 ========       UI(18) dpint(18) rhoi(18) bruni(18) riloc(18)
          !           18 -------- U(18)
          !           19 ========       UI(19)           rhoi(19) bruni(19) riloc(19)
          !
          !-----------------------------------------------------------------------     

          DO k=2,km
             shear     =  (basicum(i,k) - basicum(i,k-1))/dpint(i,k) * ( rhoi(i,k)*con_g )
             IF ( ABS(shear) .LT. shmin ) THEN
                riloc(k) = rimax
             ELSE
                riloc(k)  = (bruni(i,k)/shear) ** 2 
                IF (riloc(k) .GE. rimax ) riloc(k) = rilarge
             END IF
          ENDDO

          riloc(1)    = riloc(2)
          riloc(km+1) = riloc(km)

          !     if (lprnt.and.(i.eq.ipr)) then
          !       if (fhour.ge.fhourpr) then
          !         write(*,9104) ucltop,vcltop,windcltop,angle,kk
          !         do ilev=1,km
          !           write(*,9114) ilev,basicui(ipr,ilev),dpint(ipr,ilev),
          !    +      rhoi(ipr,ilev),(100.*bruni(ipr,ilev)),riloc(ilev)
          !           write(*,9124) ilev,(basicum(ipr,ilev))
          !         enddo
          !         ilev=km+1
          !         write(*,9114) ilev,basicui(ipr,ilev),dpint(ipr,ilev),
          !    +      rhoi(ipr,ilev),(100.*bruni(ipr,ilev)),riloc(ilev)
          !       endif
          !     endif

          !9104      FORMAT(//,'WIND VECTOR AT CLOUDTOP = (',f6.2,' , ',f6.2,' ) = ', &
          !               f6.2,' IN DIRECTION ',f6.2,4x,'KK = ',i2,//,&
          !               ' ILEV',2x,'BASICUM',2x,'BASICUI',4x,'DPINT',6x,'RHOI',5x,&
          !               'BRUNI',6x,'RI',/)
          !9114      FORMAT(i4,10x,f8.2,4(2x,f8.2))
          !9124      FORMAT(i4,1x,f8.2)

          !-----------------------------------------------------------------------
          !
          !  Calculate gravity wave stress at the interface level cloud top
          !      
          !  kcldtopi  : The interface level cloud top index
          !  kcldtop   : The midlevel cloud top index
          !  kcldbot   : The midlevel cloud bottom index
          !
          !  A : Find deep convective heating rate maximum
          !
          !      If kcldtop(i) is less than kcldbot(i) in a horizontal grid point,
          !      it can be thought that there is deep convective cloud. However,
          !      deep convective heating between kcldbot and kcldtop is sometimes
          !      zero in spite of kcldtop less than kcldbot. In this case,
          !      maximum deep convective heating is assumed to be 1.e-30. 
          !
          !  B : kk is the vertical index for interface level cloud top
          !
          !  C : Total convective fractional cover (cldbar) is used as the
          !      convective cloud cover for GWDC calculation instead of   
          !      convective cloud cover in each layer (concld).
          !                       a1 = cldbar*dlength
          !      You can see the difference between cldbar(i) and concld(i)
          !      in (4.a.2) in Description of the NCAR Community Climate    
          !      Model (CCM3).
          !      In NCAR CCM3, cloud fractional cover in each layer in a deep
          !      cumulus convection is determined assuming total convective
          !      cloud cover is randomly overlapped in each layer in the 
          !      cumulus convection.
          !
          !  D : Wave stress at cloud top is calculated when the atmosphere
          !      is dynamically stable at the cloud top
          !
          !  E : Cloud top wave stress and nonlinear parameter are calculated 
          !      using density, temperature, and wind that are defined at mid
          !      level just below the interface level in which cloud top wave
          !      stress is defined.
          !      Nonlinct is defined at the interface level.
          !  
          !  F : If the atmosphere is dynamically unstable at the cloud top,
          !      GWDC calculation in current horizontal grid is skipped.  
          !
          !  G : If mean wind at the cloud top is less than zero, GWDC
          !      calculation in current horizontal grid is skipped.
          !
          !  H : Maximum cloud top stress, tauctmax =  -20 N m^(-2),
          !      in order to prevent numerical instability.
          !
          !-----------------------------------------------------------------------
          !D
          IF ( basicui(i,kcldtop(i)) > zero ) THEN 
             !E
             IF ( riloc(kcldtop(i)) > ricrit ) THEN
                nonlinct  = ( con_g*qmax(i)*cldbar(i)*dlength(i) )/ &
                     (bruni(i,kcldtop(i))*t(i,kcldtop(i))*(basicum(i,kcldtop(i))**2))
                tauct(i)  = - (rhom(i,kcldtop(i))*(basicum(i,kcldtop(i))**2)) &
                     /   (bruni(i,kcldtop(i))*dlength(i))                  &
                     * basicum(i,kcldtop(i))*c1*c2*c2*nonlinct*nonlinct
                tauctx(i) = tauct(i)*cosphi
                taucty(i) = tauct(i)*sinphi
             ELSE
                !F
                tauct(i)  = zero
                tauctx(i) = zero 
                taucty(i) = zero
                go to 1000
             END IF
          ELSE
             !G
             tauct(i)  = zero
             tauctx(i) = zero 
             taucty(i) = zero
             go to 1000

          END IF
          !H
          IF ( tauct(i) .LT. tauctmax ) THEN
             tauct(i)  = tauctmax
             tauctx(i) = tauctmax*cosphi
             taucty(i) = tauctmax*sinphi
          END IF

          !      if (lprnt.and.(i.eq.ipr)) then
          !        if (fhour.ge.fhourpr) then
          !           write(*,9210) tauctx(ipr),taucty(ipr),tauct(ipr),angle,kk
          !        endif
          !      endif

          !9210      FORMAT(/,5x,'STRESS VECTOR = ( ',f8.3,' , ',f8.3,' ) = ',f8.3,&
          !               ' IN DIRECTION ',f6.2,4x,'KK = ',i2,/)

          !-----------------------------------------------------------------------
          !
          !  At this point, mean wind at the cloud top is larger than zero and
          !  local RI at the cloud top is larger than ricrit (=0.25)
          !
          !  Calculate minimum of Richardson number including both basic-state
          !  condition and wave effects.
          !
          !          g*Q_0*alpha*dx                  RI_loc*(1 - mu*|c2|)
          !  mu  =  ----------------  RI_min =  -----------------------------
          !           c_p*N*T*U^2                (1 + mu*RI_loc^(0.5)*|c2|)^2
          !
          !  Minimum RI is calculated for the following two cases
          !
          !  (1)   RIloc < 1.e+20  
          !  (2)   Riloc = 1.e+20  ----> Vertically uniform basic-state wind
          !
          !  RIloc cannot be smaller than zero because N^2 becomes 1.E-32 in the
          !  case of N^2 < 0.. Thus the sign of RINUM is determined by 
          !  1 - nonlin*|c2|.
          !
          !-----------------------------------------------------------------------

          DO k=kcldtop(i),1,-1

             IF ( k .NE. 1 ) THEN
                crit1 = ucltop*(u(i,k)+u(i,k-1))*0.5_r8
                crit2 = vcltop*(v(i,k)+v(i,k-1))*0.5_r8
             ELSE
                crit1 = ucltop*u(i,1)
                crit2 = vcltop*v(i,1)
             END IF

             IF((basicui(i,k) > zero).AND.(crit1 > zero).AND. (crit2 > zero)) THEN
                nonlin   = ( con_g*qmax(i)*cldbar(i)*dlength(i) )/    &
                     ( bruni(i,k)*ti(i,k)*(basicui(i,k)**2) )
                IF ( riloc(k)  <  rimaxm ) THEN
                   rimin(k) = riloc(k)*( 1 - nonlin*ABS(c2) ) /          &
                        ( 1 + nonlin*SQRT(riloc(k))*ABS(c2) )**2
                ELSE IF((riloc(k) > rimaxm).AND. (riloc(k) < rimaxp))THEN
                   rimin(k) = ( 1 - nonlin*ABS(c2) ) /( (nonlin**2)*(c2**2) ) 
                END IF
                IF ( rimin(k) <= riminx ) THEN
                   rimin(k) = rismall
                END IF
             ELSE
                rimin(k) = riminx
             END IF
          END DO

          !-----------------------------------------------------------------------
          !
          !  If minimum RI at interface cloud top is less than or equal to 1/4,
          !  GWDC calculation for current horizontal grid is skipped 
          !
          !-----------------------------------------------------------------------

          !-----------------------------------------------------------------------
          !
          !  Calculate gravity wave stress profile using the wave saturation
          !  hypothesis of Lindzen (1981).   
          !
          !  Assuming kcldtop(i)=10 and kcldbot=16,
          !
          !                             TAUGWCI  RIloc  RImin   UTGWC
          !
          !  Interface 1 ========       - 0.001         -1.e20
          !            1 --------                               0.000
          !            2 ========       - 0.001         -1.e20
          !            2 --------                               0.000
          !            3 ========       - 0.001         -1.e20
          !            3 --------                               -.xxx 
          !            4 ========       - 0.001  2.600  2.000
          !            4 --------                               0.000
          !            5 ========       - 0.001  2.500  2.000
          !            5 --------                               0.000
          !            6 ========       - 0.001  1.500  0.110
          !            6 --------                               +.xxx 
          !            7 ========       - 0.005  2.000  3.000
          !            7 --------                               0.000
          !            8 ========       - 0.005  1.000  0.222
          !            8 --------                               +.xxx
          !            9 ========       - 0.010  1.000  2.000
          !            9 --------                               0.000
          ! kcldtopi  10 ========  $$$  - 0.010 
          ! kcldtop   10 --------  $$$                          yyyyy
          !           11 ========  $$$  0
          !           11 --------  $$$
          !           12 ========  $$$  0
          !           12 --------  $$$
          !           13 ========  $$$  0
          !           13 --------  $$$
          !           14 ========  $$$  0
          !           14 --------  $$$
          !           15 ========  $$$  0
          !           15 --------  $$$
          !           16 ========  $$$  0
          ! kcldbot   16 --------  $$$
          !           17 ========       0
          !           17 -------- 
          !           18 ========       0
          !           18 -------- 
          !           19 ========       0
          !
          !-----------------------------------------------------------------------
          !
          !   Even though the cloud top level obtained in deep convective para-
          !   meterization is defined in mid-level, the cloud top level for
          !   the GWDC calculation is assumed to be the interface level just 
          !   above the mid-level cloud top vertical level index.
          !
          !-----------------------------------------------------------------------

          taugwci(i,kcldtop(i)) = tauct(i)                          !  *1

          DO k=kcldtop(i)-1,2,-1
             IF ( ABS(taugwci(i,k+1)) > taumin ) THEN                ! TAUGWCI
                IF ( riloc(k) > ricrit ) THEN                         ! RIloc
                   IF ( rimin(k) > ricrit ) THEN                       ! RImin
                      taugwci(i,k) = taugwci(i,k+1)
                   ELSE IF ((rimin(k) > riminp) .AND. &
                        (rimin(k) <= ricrit)) THEN
                      nonlins = (1.0_r8/ABS(c2))*( 2.0_r8*SQRT(2.0_r8 + 1.0_r8/SQRT(riloc(k)) ) &
                           - ( 2.0_r8 + 1.0_r8/SQRT(riloc(k)) )    )
                      taus    =  - ( rhoi(i,k)*( basicui(i,k)**2 ) )/            &
                           ( bruni(i,k)*dlength(i) ) *                   &
                           basicui(i,k)*c1*c2*c2*nonlins*nonlins
                      taugwci(i,k) = taus
                   ELSE IF((rimin(k) > riminm) .AND.(rimin(k) < riminp)) THEN
                      taugwci(i,k) = zero 
                   END IF                                              ! RImin
                ELSE

!!!!!!!!!! In the dynamically unstable environment, there is no gravity 
!!!!!!!!!! wave stress

                   taugwci(i,k) = zero    
                END IF                                                ! RIloc
             ELSE
                taugwci(i,k) = zero
             END IF                                                  ! TAUGWCI

             IF ( (basicum(i,k+1)*basicum(i,k) ) .LT. 0.0_r8 ) THEN
                taugwci(i,k+1) = zero
                taugwci(i,k)   = zero
             ENDIF

             IF (ABS(taugwci(i,k)) .GT. ABS(taugwci(i,k+1))) THEN
                taugwci(i,k) = taugwci(i,k+1)
             END IF

          END DO

!!!!!! Upper boundary condition to permit upward propagation of gravity  
!!!!!! wave energy at the upper boundary 

          taugwci(i,1) = taugwci(i,2)

          !-----------------------------------------------------------------------
          !
          !  Calculate zonal and meridional wind tendency 
          !
          !-----------------------------------------------------------------------

          DO k=1,km+1
             taugwcxi(i,k) = taugwci(i,k)*cosphi
             taugwcyi(i,k) = taugwci(i,k)*sinphi
          END DO

!!!!!! Vertical differentiation
!!!!!!
          DO k=1,kcldtop(i)-1
             tem1 = con_g / dpmid(i,k)
             wtgwc(i,k) = tem1 * (taugwci(i,k+1)  - taugwci(i,k))
             utgwc(i,k) = tem1 * (taugwcxi(i,k+1) - taugwcxi(i,k))
             vtgwc(i,k) = tem1 * (taugwcyi(i,k+1) - taugwcyi(i,k))
          END DO

          DO k=kcldtop(i),km
             wtgwc(i,k) = zero
             utgwc(i,k) = zero
             vtgwc(i,k) = zero
          END DO

          !-----------------------------------------------------------------------
          !
          !  Calculate momentum flux = stress deposited above cloup top
          !  Apply equal amount with opposite sign within cloud
          !
          !-----------------------------------------------------------------------

          xstress = zero
          ystress = zero
          DO k=1,kcldtop(i)-1
             xstress = xstress + utgwc(i,k)*dpmid(i,k)/con_g 
             ystress = ystress + vtgwc(i,k)*dpmid(i,k)/con_g  
          END DO

          !-----------------------------------------------------------------------
          !        ALT 1      ONLY UPPERMOST LAYER
          !-----------------------------------------------------------------------

          !C     kk = kcldtop(i)
          !C     tem1 = g / dpmid(i,kk)
          !C     utgwc(i,kk) = - tem1 * xstress
          !C     vtgwc(i,kk) = - tem1 * ystress

          !-----------------------------------------------------------------------
          !        ALT 2      SIN(KT-KB)
          !-----------------------------------------------------------------------

          kk = kcldtop(i)
          kb = kcldbot(i)
          DO k=kk,kb
             p1=pi/2.0_r8*(pint(i,k)-pint(i,kk))/    &
                  (pint(i,kb+1)-pint(i,kk))
             p2=pi/2.0_r8*(pint(i,k+1)-pint(i,kk))/  &
                  (pint(i,kb+1)-pint(i,kk))
             utgwc(i,k) = - con_g*xstress*(SIN(p2)-SIN(p1))/dpmid(i,k)
             vtgwc(i,k) = - con_g*ystress*(SIN(p2)-SIN(p1))/dpmid(i,k)
          ENDDO

          !-----------------------------------------------------------------------
          !        ALT 3      FROM KT to KB  PROPORTIONAL TO CONV HEATING
          !-----------------------------------------------------------------------

          !     do k=kcldtop(i),kcldbot(i)
          !     p1=cumchr(i,k)
          !     p2=cumchr(i,k+1)
          !     utgwc(i,k) = - g*xstress*(p1-p2)/dpmid(i,k)
          !     enddo

          !-----------------------------------------------------------------------
          !
          !  The GWDC should accelerate the zonal and meridional wind in the   
          !  opposite direction of the previous zonal and meridional wind, 
          !  respectively
          !
          !-----------------------------------------------------------------------

          !     do k=1,kcldtop(i)-1

          !      if (utgwc(i,k)*u(i,k) .gt. 0.0) then

          !-------------------- x-component-------------------

          !       write(6,'(a)')   
          !    +  '(GWDC) WARNING: The GWDC should accelerate the zonal wind '
          !       write(6,'(a,a,i3,a,i3)')   
          !    +  'in the opposite direction of the previous zonal wind', 
          !    +  ' at I = ',i,' and J = ',lat
          !       write(6,'(4(1x,e17.10))') u(i,kk),v(i,kk),u(i,k),v(i,k)
          !       write(6,'(a,1x,e17.10))') 'Vcld . V =',
          !    +  u(i,kk)*u(i,k)+v(i,kk)*v(i,k)

          !       if(u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k).gt.0.0)then
          !       do k1=1,km
          !         write(6,'(i2,36x,2(1x,e17.10))')
          !    +             k1,taugwcxi(i,k1),taugwci(i,k1)
          !         write(6,'(i2,2(1x,e17.10))') k1,utgwc(i,k1),u(i,k1) 
          !       end do
          !       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwcxi(i,km+1)
          !       end if

          !-------------------- Along wind at cloud top -----

          !       do k1=1,km
          !         write(6,'(i2,36x,2(1x,e17.10))')
          !    +             k1,taugwci(i,k1)
          !         write(6,'(i2,2(1x,e17.10))') k1,wtgwc(i,k1),basicum(i,k1) 
          !       end do
          !       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwci(i,km+1)

          !      end if

          !      if (vtgwc(i,k)*v(i,k) .gt. 0.0) then
          !       write(6,'(a)')
          !    +  '(GWDC) WARNING: The GWDC should accelerate the meridional wind'
          !       write(6,'(a,a,i3,a,i3)')
          !    +  'in the opposite direction of the previous meridional wind',
          !    +  ' at I = ',i,' and J = ',lat
          !       write(6,'(4(1x,e17.10))') u(i,kcldtop(i)),v(i,kcldtop(i)),
          !    +                            u(i,k),v(i,k)
          !       write(6,'(a,1x,e17.10))') 'Vcld . V =',
          !    +                    u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k)
          !       if(u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k).gt.0.0)then
          !       do k1=1,km
          !         write(6,'(i2,36x,2(1x,e17.10))')
          !    +                        k1,taugwcyi(i,k1),taugwci(i,k1)
          !         write(6,'(i2,2(1x,e17.10))') k1,vtgwc(i,k1),v(i,k1) 
          !       end do
          !       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwcyi(i,km+1)
          !       end if
          !      end if

          !     enddo

1000      CONTINUE

       END IF   ! DO GWDC CALCULATION

    END DO   ! I-LOOP 

    !***********************************************************************

    !      if (lprnt) then
    !        if (fhour.ge.fhourpr) then
    !-------- UTGWC VTGWC ----------
    !          write(*,9220)
    !          do ilev=1,km
    !            write(*,9221) ilev,(86400.*utgwc(ipr,ilev)),
    !     +                         (86400.*vtgwc(ipr,ilev))
    !          enddo
    !        endif
    !      endif

    !9220 FORMAT(//,14x,'TENDENCY DUE TO GWDC',//,&
    !         ' ILEV',6x,'UTGWC',7x,'VTGWC',/)
    !9221 FORMAT(i4,2(2x,f10.3))

    !-----------------------------------------------------------------------
    !
    !  For GWDC performance analysis        
    !
    !-----------------------------------------------------------------------

    DO i = 1, im
       kk=kcldtop(i)

       IF ( dogwdc(i) .AND. (ABS(taugwci(i,kk)).GT.taumin) ) THEN

          gwdcloc(i) = one

          DO k = 1, kk-1
             IF ( ABS(taugwci(i,k)-taugwci(i,kk)).GT.taumin ) THEN
                break(i) = 1.0_r8
                go to 2000
             ENDIF
          ENDDO
2000      CONTINUE

          DO k = 1, kk-1

             IF ( ( ABS(taugwci(i,k)).LT.taumin ) .AND.  &
                  ( ABS(taugwci(i,k+1)).GT.taumin ) .AND.  &
                  ( basicum(i,k+1)*basicum(i,k) .LT. 0.0_r8 ) ) THEN
                critic(i) = 1.0_r8
                !         print *,i,k,' inside GWDC  taugwci(k) = ',taugwci(i,k)
                !         print *,i,k+1,' inside GWDC  taugwci(k+1) = ',taugwci(i,k+1)
                !         print *,i,k,' inside GWDC  basicum(k) = ',basicum(i,k)
                !         print *,i,k+1,' inside GWDC  basicum(k+1) = ',basicum(i,k+1)
                !         print *,i,' inside GWDC  critic = ',critic(i)
                GOTO 2010
             ENDIF
          ENDDO
2010      CONTINUE

       ENDIF

    ENDDO

    !-----------------------------------------------------------------------
    !        Convert back local GWDC Tendency arrays to GFS model vertical indices
    !        Outgoing (FU1,FV1)=(utgwc,vtgwc)
    !-----------------------------------------------------------------------

    DO k=1,km
       k1=km-k+1
       DO i=1,im
          fu1(i,k1)    = utgwc(i,k)
          fv1(i,k1)    = vtgwc(i,k)
          brunm1(i,k1) = brunm(i,k)
          rhom1(i,k1)  = rhom(i,k)
       ENDDO
    ENDDO

    !      if (lprnt) then
    !        if (fhour.ge.fhourpr) then
    !!-------- UTGWC VTGWC ----------
    !          write(*,9225)
    !          do ilev=km,1,-1
    !            write(*,9226) ilev,(86400.*fu1(ipr,ilev)),
    !     +                         (86400.*fv1(ipr,ilev))
    !          enddo
    !        endif
    !      endif

    !9225 FORMAT(//,14x,'TENDENCY DUE TO GWDC - TO GBPHYS',//,&
    !         ' ILEV',6x,'UTGWC',7x,'VTGWC',/)
    !9226 FORMAT(i4,2(2x,f10.3))

    RETURN
  END SUBROUTINE gwdc




END MODULE Gwdd_ECMWF
