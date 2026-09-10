!
!  $Author: pkubota $
!  $Date: 2007/02/15 17:03:10 $
!  $Revision: 1.5 $
!
MODULE Shall_Tied

USE Constants, ONLY :  &
       cp                 , &
       hl                 , &
       gasr               , &
       grav               , &
       rmwmd              , &
       rmwmdi             , &
       e0c                , &
       delq               , &
       p00                , &
       r8
USE Options, ONLY :       &
       rccmbl            , &
       sthick            , &
       sacum             , &
       acum0             , &
       tbase             , &
       ki                , &
       mlrg              , &
       is                , &
       doprec            , &
       cflric            , &
       ifilt             , &
       dt                , &
       kt                , &
       ktp               , &
       jdt               , &
       nfprt

IMPLICIT NONE
SAVE

  PRIVATE
PUBLIC :: InitShall_Tied
PUBLIC :: shalv2

  REAL(KIND=r8) :: aa(15)
  REAL(KIND=r8) :: ad(15)
  REAL(KIND=r8) :: ac(15)
  REAL(KIND=r8) :: actop
  REAL(KIND=r8) :: thetae(151,181)
  REAL(KIND=r8) :: tfmthe(431,241)
  REAL(KIND=r8) :: qfmthe(431,241)
  REAL(KIND=r8) :: ess
!  INTEGER :: kbase
!  INTEGER :: kcr
!  REAL(KIND=r8), ALLOCATABLE :: dels  (:)
!  REAL(KIND=r8), ALLOCATABLE :: gams  (:)
!  REAL(KIND=r8), ALLOCATABLE :: gammod(:)
!  REAL(KIND=r8), ALLOCATABLE :: delmod(:)
  REAL(KIND=r8) :: rlocp
  REAL(KIND=r8) :: rgrav
  REAL(KIND=r8) :: rlrv
  REAL(KIND=r8) :: const1
  REAL(KIND=r8) :: const2
  REAL(KIND=r8) :: xx1

  REAL(KIND=r8), PARAMETER :: xkapa=0.2857143_r8


CONTAINS

SUBROUTINE InitShall_Tied()
    !INTEGER, INTENT(IN) :: kmax
    !REAL(KIND=r8),    INTENT(IN) :: si (kmax+1)
    !REAL(KIND=r8),    INTENT(IN) :: del(kmax  )
    !REAL(KIND=r8),    INTENT(IN) :: sl (kmax  )
    !REAL(KIND=r8),    INTENT(IN) :: cl (kmax  )
    !CALL InitShalv2(si, del, sl, cl, kmax)

END SUBROUTINE InitShall_Tied


  ! shalv2 :subgrid scale shallow cumulus parameterization - enhanced
  !         vertical temperature and moisture diffusion returning adjusted
  !         temperature and specific humidity.


  
  SUBROUTINE shalv2( tin, qin,prsi ,prsl , deltim, &
       ktop, plcl, kuo, kmaxp, kctop1, kcbot1, noshal1, &
       newr, ncols, kmax,dtdt ,dqdt )
    !
    !**************************************************************************
  !
    !         use:
    !         sr to compute the parameterized effects of shallow
    !         convective clouds (vector version).  routine used in
    !         gwater  moist convection kuolcl.
    !         
    !         purpose:
    !         sub-grid-scale shallow convective cloud parameterization.
    !         this routine computes the effects of shallow convection
    !         based on tiedtke (1984), ecmwf workshop on convection in
    !         large-scale numerical models.
    !         tapered k profile in cloud  developed by caplan and long.
    !         
    !         srs called ... none
    !
    !         
    !         jpg modifications nov 6 85
    ! 
    ! input
    !     
    !     si    : p/ps at sigma level interfaces (array from surface up)
    !     del   : del p positive across sigma lyr  (array from surface up)
    !     sl    : p/ps at sigma layers               (array from surface up)
    !     cl    : 1 - sl                     (array from surface up)
    !     tin   : temperature       (longitude, height array, deg k)
    !     qin   : specific humidity (longitude, height array, gm/gm)
    !     ps    : surface pressure  (longitude array        , cb   )
    !     deltim: timestep (  sec  )
    !     ktop  : cloud tops (sigma layer index) passed from conkuo
    !     plcl  : pressure   (cb) at lcl   passed from conkuo
    !     kuo   : flag to indicate that deep convection was done
    !             kuo, ktop and plcl are longitude arrays
    !     newr  : flag for ccm3 based cloud radiation
    !     output
    !
    !     qin   : updated specific humidity    (gm/gm)
    !     tin   : updated temperature          (deg k)
    !
    !
    !     external variables
    !**************************************************************************
    ! ncols.......Number of grid points on a gaussian latitude circle
    ! kmax.......Number of sigma levels
    ! kmaxp......kmaxp=kmax+1
    ! deltim.....deltim=dt where dt time interval,usually =delt,but changes
    !            in nlnmi (dt=1.) and at dead start(delt/4,delt/2)
    ! si.........si(l)=1.0-ci(l).
    ! ci.........sigma value at each level.
    ! del........sigma spacing for each layer computed in routine "setsig".
    ! ps.........surface pressure
    ! sl.........sigma value at midpoint of
    !                                         each layer : (k=287/1005)
    !
    !                                                                     1
    !                                             +-                   + ---
    !                                             !     k+1         k+1!  k
    !                                             !si(l)   - si(l+1)   !
    !                                     sl(l) = !--------------------!
    !                                             !(k+1) (si(l)-si(l+1)!
    !                                             +-                  -+
    ! kcbot1....
    ! kctop1....
    ! noshal1...
    ! cp........Specific heat of air           (j/kg/k)
    ! hl........heat of evaporation of water     (j/kg)
    ! gasr......gas constant of dry air        (j/kg/k)
    ! grav......grav   gravity constant        (m/s**2)
    !**************************************************************************
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    INTEGER, INTENT(in   ) :: kmaxp
    REAL(KIND=r8),    INTENT(in   ) :: deltim
    REAL(KIND=r8),    INTENT(in   ) :: prsi   ( ncols , kmax +1)
    REAL(KIND=r8),    INTENT(in   ) :: prsl   ( ncols , kmax )
    REAL(KIND=r8),    INTENT(inout) :: plcl   ( ncols )
    INTEGER, INTENT(in   ) :: ktop   ( ncols )
    INTEGER, INTENT(in   ) :: kuo    ( ncols )
    REAL(KIND=r8),    INTENT(inout) :: tin    ( ncols , kmax )
    REAL(KIND=r8),    INTENT(inout) :: qin    ( ncols , kmax )
    LOGICAL, INTENT(in   ) :: newr
    INTEGER, INTENT(inout) :: kcbot1 (ncols)
    INTEGER, INTENT(inout  ) :: kctop1 (ncols)
    INTEGER, INTENT(inout  ) :: noshal1(ncols)
    REAL(KINd=r8)   ,    INTENT(OUT  ) :: dtdt (ncols,kMax)
    REAL(KIND=r8)   ,    INTENT(OUT  ) :: dqdt (ncols,kMax)

    REAL(KIND=r8)    :: dk(    ncols,kmax-1)
    REAL(KIND=r8)    :: terp  (ncols,kmax-1)
    REAL(KIND=r8)    :: f     (ncols,kmax  )
    REAL(KIND=r8)    :: ff    (ncols,kmax  )
    REAL(KIND=r8)    :: g     (ncols,kmax  )
    REAL(KIND=r8)    :: gg    (ncols,kmax  )
    REAL(KIND=r8)    :: tnew  (ncols,kmax  )
    REAL(KIND=r8)    :: qnew  (ncols,kmax  )
    REAL(KIND=r8)    :: a     (ncols,kmax  )
    REAL(KIND=r8)    :: b     (ncols,kmax  )
    REAL(KIND=r8)    :: c     (ncols,kmax  )
    REAL(KIND=r8)    :: DeltaP(ncols,kmax  )
    REAL(KIND=r8)    :: ud    (ncols,kmax-1)
    REAL(KIND=r8)    :: dels  (ncols,kmax-1)
    REAL(KIND=r8)    :: delmod(ncols,kmax-1)
    REAL(KIND=r8)    :: gams  (ncols,kmax-1)
    REAL(KIND=r8)    :: gammod(ncols,kmax-1)
    INTEGER          :: kbase (ncols)
    INTEGER          :: kcr   (ncols)
    REAL(KIND=r8)    :: ps     ( ncols )
    REAL(KIND=r8)    :: sl     (ncols,kmax  )
    REAL(KIND=r8)    :: si     (ncols,kmax+1)


    REAL(KIND=r8)    :: rec   (ncols)
    INTEGER :: icheck(ncols)
    INTEGER :: kctop (ncols)
    INTEGER :: kcbot (ncols)
    INTEGER :: noshal(ncols)
    LOGICAL :: searching(ncols)

    INTEGER ::  n
    INTEGER ::  i
    INTEGER ::  l
    INTEGER ::  k

    REAL(KIND=r8)          :: dt2,gr,ggrr,dryl
    INTEGER       :: kmaxm
    INTEGER       :: levhm1
    INTEGER       :: loleh1
    INTEGER       :: lolem1
    INTEGER       :: lonlev
    INTEGER       :: lonleh

    dt2 = deltim* 2.0e0_r8
    kmaxm = kmax  -1
    levhm1= kmax -1
    loleh1= ncols *levhm1
    lolem1 =  ncols * kmaxm
    lonlev =  ncols * kmax
    lonleh =  ncols * kmax

    !
    !     zero  arrays
    !
    DO k=1,kMax
      DO i=1,ncols
          dtdt (i,k)=0.0_r8
          dqdt (i,k)=0.0_r8
      END DO
    END DO

    DO i=1,ncols
       noshal(i)=0
       icheck(i)=0
       ps     ( i )=prsi(i,1)/1000.0_r8
    END DO
    gr=2.0e0_r8*grav/gasr
    ggrr=gr*gr
    dryl=grav/cp
    DO k=1,kMax
      DO i=1,ncols
          sl(i,k)=prsl(i,k)/prsi(i,1) 
          DeltaP(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO
    DO k=1,kMax+1
      DO i=1,ncols
          si(i,k)=prsi(i,k)/prsi(i,1) 
      END DO
    END DO

    !cl(k) = 1.0_r8 - sl(k)


    DO k=1,kMax-1
      DO i=1,ncols
          !dels(k)=ggrr*si(k+1)*si(k+1)/ ( DeltaP(i,k)*( (1-sl(k+1)) - (1-sl(k)) )    )
          !dels(k)=ggrr*si(k+1)*si(k+1)/ ( DeltaP(i,k)*( 1-sl(k+1) - 1+sl(k) )    )
          !dels(k)=ggrr*si(k+1)*si(k+1)/ ( DeltaP(i,k)*( -sl(k+1) +sl(k) )    )
          !dels(k)=ggrr*si(k+1)*si(k+1)/ ( DeltaP(i,k)*( sl(k) - sl(k+1) ))
          !              ggrr*si(k+1)*si(k+1)
          !dels(k) = ------------------------------------
          !           ( DeltaP(i,k)*( sl(k) - sl(k+1)))

          !              ggrr*prsi(i,k+1)/prsi(i,1)  *   prsi(i,k+1)/prsi(i,1)
          !dels(k) = ------------------------------------------------------------------
          !           ( DeltaP(i,k)*( prsl(i,k)/prsi(i,1) - prsl(i,k+1)/prsi(i,1)))


          !              ggrr*(prsi(i,k+1)  *   prsi(i,k+1))/(prsi(i,1)*prsi(i,1))
          !dels(k) = ----------------------------------------------------------------
          !           ( DeltaP(i,k)*( prsl(i,k) - prsl(i,k+1))/prsi(i,1))


          !              ggrr*(prsi(i,k+1)  *   prsi(i,k+1))/(prsi(i,1)*prsi(i,1))
          !dels(k) = --------------------------------------------------------------------------------
          !           ( (prsi(i,k) - prsi(i,k+1))/prsi(i,1))*( prsl(i,k) - prsl(i,k+1))/prsi(i,1))


          !              ggrr*(prsi(i,k+1)  *   prsi(i,k+1))/(prsi(i,1)*prsi(i,1))
          !dels(k) = -------------------------------------------------------------------------------
          !           ( (prsi(i,k) - prsi(i,k+1))*( prsl(i,k) - prsl(i,k+1)))/(prsi(i,1)*prsi(i,1))


          !              ggrr*(prsi(i,k+1)  *   prsi(i,k+1))
          !dels(i,k) = -------------------------------------------------------------------------------
          !           ( (prsi(i,k) - prsi(i,k+1))*( prsl(i,k) - prsl(i,k+1)))

                       
          dels(i,k) =  (ggrr*(prsi(i,k+1)*prsi(i,k+1)))/( (prsi(i,k) - prsi(i,k+1))*( prsl(i,k) - prsl(i,k+1)))


          !delmod(k)=gr*dryl*si(k+1)/DeltaP(i,k)
          !
          !delmod(k)=gr*dryl*(prsi(i,k+1)/prsi(i,1))/DeltaP(i,k)
          !
          !delmod(k)=gr*dryl*(prsi(i,k+1)/prsi(i,1))/((prsi(i,k) - prsi(i,k+1))/prsi(i,1))
       
          !           gr*dryl*prsi(i,k+1)/prsi(i,1)
          !delmod(k)= --------------------------------------------
          !          (prsi(i,k) - prsi(i,k+1))/prsi(i,1)

          !           gr*dryl*prsi(i,k+1)
          !delmod(k)= --------------------------------------------
          !          (prsi(i,k) - prsi(i,k+1))

                     
          delmod(i,k)= (gr*dryl*prsi(i,k+1))/(prsi(i,k) - prsi(i,k+1))


          gams(i,k)=DeltaP(i,k)*dels(i,k)/DeltaP(i,k+1)

          gammod(i,k)=DeltaP(i,k)*delmod(i,k)/DeltaP(i,k+1)

      END DO
    END DO
    kbase=1
    kcr  =1
    DO k = 1, kMax-1
       DO i=1,ncols


          !IF(sl(k).GT.0.7e0_r8) THEN
          IF(prsl(i,k).GT.70000.0_r8) THEN
             kbase(i)=k
             kcr  (i)=k
          END IF
       END DO
    END DO


    dk=0.0_r8
    DO i=1,ncols
       IF(kuo(i) .EQ. 1) THEN
          noshal(i)=1
       END IF
    END DO
    !
    !     get cloud base .. overwrite plcl with normalized value
    !
    DO i=1,ncols
       !PRINT*,'pkubota=',plcl(i),ps(i),plcl(i)/ps(i)
       plcl(i)=plcl(i)/ps(i)
    END DO

    !    DO  l=1, ncols
    !       DO n=2, kmax
    !          kcbot(l)=n-1
    !          IF ( plcl(l) .GE. sl(n) ) THEN
    !             EXIT
    !          END IF
    !       END DO
    !    END DO

    searching = .TRUE.
    kcbot = kmax-1
    DO n=2, kmax
       DO  l=1, ncols
          IF (searching(l) .AND. plcl(l) >=  sl(l,n)) THEN
             searching(l) = .FALSE.
             kcbot(l)=n-1
          END IF
       END DO
    END DO
    IF (newr) THEN
       kcbot1(1:ncols) = kcbot(1:ncols)
    END IF
    !
    !     set cloud tops
    !
    DO i=1,ncols
       kctop(i) = MIN(kcr(i),ktop(i))
       IF(newr)kctop1(i)=kctop(i)
    END DO
    !
    !     check for too high cloud base
    !     remem. that plcl has been divided by ps
    !
    DO i=1,ncols
       IF(plcl(i).LT. sl(i,kbase(i))) THEN
          noshal(i) = 1
       END IF
    END DO

    DO i=1,ncols
       IF(kcbot(i).GE.kctop(i)) THEN
          noshal(i) = 1
          kcbot (i) = kctop(i)
       END IF
    END DO
    !
    !     test for moist convective instability
    !
    DO k=1,kmax-1
       DO i=1,ncols
          a(i,k+1) =  tin(i,k+1) - tin(i,k)
          b(i,k+1) =  qin(i,k+1) - qin(i,k)
          c(i,k+1) =  tin(i,k+1) + tin(i,k)
          a(i,k  ) =  &
               (cp/(sl(i,k+1)-sl(i,k)))*a(i,k+1)-&
               (0.5e0_r8*gasr/si(i,k+1))*c(i,k+1)+&
               (hl/(sl(i,k+1)-sl(i,k)))*b(i,k+1)
       END DO
    END DO

    DO k = 1, kmax
       DO i=1, ncols
          IF (k >= kcbot(i) .AND. k <= kctop(i) .AND. a(i,k) > 0.0e0_r8) THEN
             icheck(i)=1
          END IF
       END DO
    END DO

    DO i=1,ncols
       IF(icheck(i).EQ.0) THEN
          noshal(i)=1
       END IF
    END DO
    !
    !     set mixing coefficient dk (m**2/s) profile
    !     dk(n) is value at top of layer n
    !     n.b.  dk(kctop) .ne. 0
    !
    !DO l=1, ncols
    !   IF (noshal(l) .NE. 1) THEN
    !      dk(l,kctop(l)) = 1.0e0_r8
    !      dk(l,kcbot(l)) = 1.5e0_r8
    !      kbetw= kctop(l)-(kcbot(l)+1)
    !      IF(kbetw .GE. 1) dk(l,kctop(l)-1)=3.0e0_r8
    !      IF(kbetw .GT. 1) THEN
    !         DO k=kcbot(l)+1, kctop(l)-2
    !            dk(l,k)=5.0e0_r8
    !         END DO
    !      END IF
    !   END IF
    !END DO
    DO l=1, ncols
       IF (noshal(l) /= 1) THEN
          dk(l,kctop(l)) = 1.0e0_r8
          dk(l,kcbot(l)) = 1.5e0_r8
       END IF
    END DO
    DO l=1, ncols
       IF (noshal(l) /= 1 .AND. kctop(l) >= kcbot(l)+2) THEN
          dk(l,kctop(l)-1)=3.0e0_r8
       END IF
    END DO
    DO k = 1, kmax
       DO l = 1, ncols
          IF (&
               noshal(l) /= 1         .AND. &
               kctop(l) >  kcbot(l)+2 .AND. &
               k >= kcbot(l)+1        .AND. &
               k <= kctop(l)-2               ) THEN
             dk(l,k)=5.0e0_r8
          END IF
       END DO
    END DO
    !
    !     compute adiabatic lapse rate terms for temperature eq.
    !
    DO k=1,kmaxm
       DO i=1,ncols
          terp(i,k) = 1.0e0_r8/(tin(i,k)+tin(i,k+1))
       END DO
    END DO
    !
    !     n.b. terp(i,n) is valid at top of layer n
    !
    DO k=1, kmaxm
       DO i=1,ncols
          ff(i,k) = (dt2*delmod(i,k))*terp(i,k)*dk(i,k)
          gg(i,k) = (dt2*gammod(i,k))/(dt2*delmod(i,k))*ff(i,k)
       END DO

    END DO
    !
    !     compute elements of tridiagonal matrix a,b,c
    !
    a=0.0_r8 ! call reset(a,lonlev)
    b=0.0_r8 ! call reset(b,lonlev)
    c=0.0_r8 ! call reset(c,lonlev)

    DO k=1,kmaxm
       DO i=1,ncols
          terp(i,k) = terp(i,k)*terp(i,k)*dk(i,k)
       END DO
    END DO

    DO n=1, kmaxm
       DO i=1,ncols
          c(i,n) = (-dels(i,n) *dt2)*terp(i,n)
       END DO
    END DO

    DO n=2, kmax
       DO i=1,ncols
          a(i,n) = (-gams(i,n-1)*dt2) *terp(i,n-1)
       END DO
    END DO

    DO k=1,kmax
       DO i=1,ncols
          b(i,k)=  1.0e0_r8-a(i,k)-c(i,k)
       END DO
    END DO
    !
    !     compute forcing terms f for temperature and
    !     g for water vapor
    !
    DO i=1,ncols
       f(i,   1)=tin(i,   1)+ff(i,    1)
       f(i,kmax)=tin(i,kmax)-gg(i,kmaxm)
    END DO

    DO n=2, kmaxm
       DO i=1,ncols
          f(i,n)=tin(i,n)+ff(i,n)-gg(i,n-1)
       END DO
    END DO

    DO k=1,kmax
       DO i=1,ncols
          g(i,k)=qin(i,k)
       END DO
    END DO
   !
    !     solution of tridiagonal problems
    !     forward part
    !
    DO i=1,ncols
       rec (i)=1.0e0_r8/b(i,1)
       ud(i,1)=c(i,1)*rec(i)
       ff(i,1)=f(i,1)*rec(i)
       gg(i,1)=g(i,1)*rec(i)
    END DO

    DO k=2, kmaxm
       DO i=1,ncols
          rec(i)=1.0e0_r8/( b(i,k)-a(i,k)*ud(i,k-1))
          ud(i,k)=c(i,k)*rec(i)
          ff(i,k)=rec(i)*(f(i,k)-a(i,k)*ff(i,k-1))
       END DO

       IF(k .LE. levhm1)THEN
          DO i=1,ncols
             gg(i,k)=rec(i)*(g(i,k)-a(i,k)*gg(i,k-1))
          END DO
       END IF

    END DO
    !
    !     now determine solutions
    !
    DO i=1,ncols
       rec(i)=1.0e0_r8/(b(i,kmax)-a(i,kmax)*ud(i,kmaxm))
       tnew(i,kmax)=rec(i)*(f(i,kmax)-a(i,kmax)*ff(i,kmaxm))
    END DO

    DO n= kmaxm ,1,-1
       DO i=1,ncols
          tnew(i,n)=ff(i,n)-tnew(i,n+1)*ud(i,n)
       END DO
    END DO
    !     water vapor solution
    !
    DO i=1,ncols
       rec(i)=1.0e0_r8/(b(i,kmax)-a(i,kmax)*ud(i,levhm1))
       qnew(i,kmax)=rec(i)*(g(i,kmax)-a(i,kmax)*gg(i,levhm1))
       IF(newr)noshal1(i)=noshal(i)
    END DO

    DO k=kmax-1,1,-1
       DO i=1,ncols
          qnew(i,k)=gg(i,k)-ud(i,k)*qnew(i,k+1)
       END DO
    END DO
    !
    !     store solutions
    !
    DO k=1,kmax
       DO i=1,ncols
          IF(noshal(i).EQ.0)THEN
             dtdt (i,k)=(tnew(i,k)-tin(i,k))/dt2
             tin(i,k)=tnew(i,k)
          END IF
       END DO
    END DO

    DO k=1,kmax
       DO i=1,ncols
          IF(noshal(i).EQ.0)THEN
             dqdt (i,k)=(qnew(i,k)-qin(i,k))/dt2
             qin(i,k)=qnew(i,k)
          END IF
       END DO
    END DO
    !
    !     end of computation of shallow convection
    !
!999 FORMAT(' PARAMETER SETTING IS WRONG IN SUBR.SHALV2')
  END SUBROUTINE shalv2

!SUBROUTINE InitShalv2(si, del, sl, cl, kmax)
!    INTEGER, INTENT(IN) :: kmax
!    REAL(KIND=r8),    INTENT(IN) :: si (kmax+1)
!    REAL(KIND=r8),    INTENT(IN) :: del(kmax  )
!    REAL(KIND=r8),    INTENT(IN) :: sl (kmax  )
!    REAL(KIND=r8),    INTENT(IN) :: cl (kmax  )!       cl(k) = 1.0_r8 - sl(k)
    
!   INTEGER :: n 
!   REAL(KIND=r8)    :: gr
!   REAL(KIND=r8)    :: ggrr
!   REAL(KIND=r8)    :: dryl
! 
!   ALLOCATE (dels  (kmax-1))
!   ALLOCATE (gams  (kmax-1))
!   ALLOCATE (gammod(kmax-1))
!   ALLOCATE (delmod(kmax-1))

!   gr=2.0e0_r8*grav/gasr
!   ggrr=gr*gr
!   dryl=grav/cp
!   DO n = 1, kmax-1
!      dels(n)=ggrr*si(n+1)*si(n+1)/ ( &
!           del(n)*( cl(n+1) - cl(n) )    )
!      delmod(n)=gr*dryl*si(n+1)/del(n)
!      gams(n)=del(n)*dels(n)/del(n+1)
!      gammod(n)=del(n)*delmod(n)/del(n+1)
!      IF(sl(n).GT.0.7e0_r8) THEN
!         kbase=n
!         kcr=n
!      END IF
!   END DO
!  END SUBROUTINE InitShalv2

END MODULE Shall_Tied
