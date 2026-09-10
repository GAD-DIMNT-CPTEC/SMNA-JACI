!
!  $Author: pkubota $
!  $Date: 2010/04/20 20:18:04 $
!  $Revision: 1.7 $
!
MODULE Cu_Kuolcl

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

USE Parallelism, ONLY: &
     MsgOne, &
     FatalError

IMPLICIT NONE
SAVE

  PRIVATE
  PUBLIC :: Init_Cu_Kuolcl
  PUBLIC :: RunCu_Kuolcl

  REAL(KIND=r8) :: aa(15)
  REAL(KIND=r8) :: ad(15)
  REAL(KIND=r8) :: ac(15)
  REAL(KIND=r8) :: actop

  !----------------------------------------------------------------------
  ! MSTAD2 - Tables and dimensions
  !----------------------------------------------------------------------

  ! hmjb - Eh melhor deixar estas matrizes com tamanho variavel,
  !   assim eh facil extender as tabelas caso seja necessario
  !   para evitar os erros do tipo mstad2 ier /= 0

  ! Temperature in Kelvin
  INTEGER      , PARAMETER :: the_ntemp = 171      ! Number of steps
  REAL(KIND=r8), PARAMETER :: the_dtemp = 1.000_r8 ! Step in temperature
  REAL(KIND=r8), PARAMETER :: the_temp0 = 160.0_r8 ! First temperature

  ! Pressure in ??
  INTEGER      , PARAMETER :: the_npres = 181      ! Number of steps
  REAL(KIND=r8), PARAMETER :: the_dpres = 0.005_r8 ! Step in pressure
  REAL(KIND=r8), PARAMETER :: the_dpinv = 1.0_r8/the_dpres
  REAL(KIND=r8), PARAMETER :: the_pres0 = 0.300_r8 ! First pressure

  REAL(KIND=r8) :: thetae(the_ntemp,the_npres)

  ! Temperature in Kelvin
  INTEGER      , PARAMETER :: fm_nthe = 431      ! Number of steps
  REAL(KIND=r8), PARAMETER :: fm_dthe = 1.000_r8 ! Step in theta_e
  REAL(KIND=r8), PARAMETER :: fm_the0 = 170.0_r8 ! First theta_e

  ! Pressure in ??
  INTEGER      , PARAMETER :: fm_npres = 241      ! Number of steps
  REAL(KIND=r8), PARAMETER :: fm_dpres = 0.005_r8 ! Step in pressure
  REAL(KIND=r8), PARAMETER :: fm_dpinv = 1.0_r8/fm_dpres
  REAL(KIND=r8), PARAMETER :: fm_pres0 = 0.000_r8 ! First pressure

  REAL(KIND=r8) :: tfmthe(fm_nthe,fm_npres)
  REAL(KIND=r8) :: qfmthe(fm_nthe,fm_npres)

  !----------------------------------------------------------------------

  REAL(KIND=r8) :: ess
  INTEGER :: kbase
  INTEGER :: kcr
  REAL(KIND=r8), ALLOCATABLE :: dels  (:)
  REAL(KIND=r8), ALLOCATABLE :: gams  (:)
  REAL(KIND=r8), ALLOCATABLE :: gammod(:)
  REAL(KIND=r8), ALLOCATABLE :: delmod(:)
  REAL(KIND=r8) :: rlocp
  REAL(KIND=r8) :: rgrav
  REAL(KIND=r8) :: rlrv
  REAL(KIND=r8) :: const1
  REAL(KIND=r8) :: const2
  REAL(KIND=r8) :: xx1
  REAL(KIND=r8), PARAMETER  :: xkapa=0.2857143_r8
  REAL(KIND=r8), PARAMETER  :: d1=0.6107042e0_r8
  REAL(KIND=r8), PARAMETER  :: d2=4.441157e-2_r8
  REAL(KIND=r8), PARAMETER  :: d3=1.432098e-3_r8
  REAL(KIND=r8), PARAMETER  :: d4=2.651396e-5_r8
  REAL(KIND=r8), PARAMETER  :: d5=3.009998e-7_r8
  REAL(KIND=r8), PARAMETER  :: d6=2.008880e-9_r8
  REAL(KIND=r8), PARAMETER  :: d7=6.192623e-12_r8


CONTAINS

  SUBROUTINE Init_Cu_Kuolcl()

    CALL InitMstad2()

  END SUBROUTINE Init_Cu_Kuolcl
 
 
  REAL(KIND=r8) FUNCTION es(t)
    REAL(KIND=r8), INTENT(IN) :: t
    REAL(KIND=r8) :: tx

    ! statement function

    tx = t - tbase
    IF (tx >= -50.0_r8) THEN
       es = d1 + tx*(d2 + tx*(d3 + tx*(d4 + tx*(d5 + tx*(d6 + d7*tx)))))
    ELSE
       es=0.00636e0_r8*EXP(25.6e0_r8*(tx+50.e0_r8)/(tbase-50.e0_r8))
    END IF
  END FUNCTION es




  SUBROUTINE RunCu_Kuolcl(dt   ,prsi ,prsl  , qn  , qn1 , &
                          tn1  ,dq   ,geshem, kuo , plcl, kktop, &
                          kkbot,ncols,kmax,dtdt,dqdt)
    !
    !==========================================================================
    !==========================================================================
    ! imx.......=ncols+1 or ncols+2   :this dimension instead of ncols
    !              is used in order to avoid bank conflict of memory
    !              access in fft computation and make it efficient. the
    !              choice of 1 or 2 depends on the number of banks and
    !              the declared type of grid variable (REAL(KIND=r8)*4,REAL(KIND=r8)*8)
    !              to be fourier transformed.
    !              cyber machine has the symptom.
    !              cray machine has no bank conflict, but the argument
    !              'imx' in subr. fft991 cannot be replaced by ncols
    ! kmax.......Number of sigma levels
    ! ncols.......Number of grid points on a gaussian latitude circle
    ! dt.........time interval,usually =delt,but changes
    !            in nlnmi (dt=1.) and at dead start(delt/4,delt/2)
    ! nkuo.......index used by routine "conkuo"
    !            to accumulate number of
    !            points for which trouble was
    !            encountered while computing
    !            a moist adiabat.
    ! ps.........sfc pres (cb)
    ! del........sigma spacing for each layer computed in routine "setsig".
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
    ! si........si(l)=1.0-ci(l).
    ! ci........sigma value at each level.
    ! qn........qn is q (specific humidit) at the n-1 time level
    ! qn1.......qn1 is q (specific humidit) at the n+1 time level
    ! tn1.......tn1 is tmp (temperature) at the n+1 time level
    ! dq.........specific humidit difference between levels
    ! geshem.....set aside convective precip in separate array for diagnostics
    ! msta.......
    ! kuo........flag to indicate that deep convection was done
    !            kuo, ktop and plcl are longitude arrays
    ! plcl.......pressure at the lcl
    ! kktop......ktop (g.t.e 1 and l.t.e. km) is highest lvl for which
    !            tin is colder than moist adiabat given by the
    !            (ktop=1 denotes tin(k=2) is already g.t.e. moist adb.)
    !            allowance is made for perhaps one level below ktop
    !            at which tin was warmer than tmst.
    ! kkbot......is the first regular level above the lcl
    ! cp.........Specific heat of air           (j/kg/k)
    ! hl.........heat of evaporation of water     (j/kg)
    ! gasr.......gas constant of dry air        (j/kg/k)
    ! g..........grav   gravity constant        (m/s**2)
    ! rmwmd......fracao molar entre a agua e o ar
    ! sthick.....upper limit for originating air for lcl.  replaces kthick.
    ! sacum......top level for integrated moisture convergence test. replaces
    !            kacum
    ! acum0......threshold moisture convergence such that integrated moisture
    !            convergence > - acum0 for convection to occur.
    ! tbase......constant tbase =  273.15e00
    ! ki.........lowest level from which parcels can be lifted to find lcl
    INTEGER      ,    INTENT(in   ) :: kmax
    INTEGER      ,    INTENT(in   ) :: ncols
    REAL(KIND=r8),    INTENT(in   ) :: dt
    REAL(KIND=r8),    INTENT(IN   ) :: prsi   (ncols,kMax+1)   !
    REAL(KIND=r8),    INTENT(IN   ) :: prsl   (ncols,kMax  )   !
    REAL(KIND=r8),    INTENT(in   ) :: qn     (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: qn1    (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: tn1    (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: dq     (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: geshem (ncols)
    INTEGER      ,    INTENT(inout) :: kuo    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: plcl   (ncols)
    INTEGER      ,    INTENT(inout) :: kktop  (ncols)
    INTEGER      ,    INTENT(inout) :: kkbot  (ncols)
    REAL(KINd=r8),    INTENT(out  ) :: dtdt   (nCols,kMax)
    REAL(KIND=r8),    INTENT(out  ) :: dqdt   (nCols,kMax)

    !
    !     these are for monitoring of gpv in gfidi.
    !
    REAL(KIND=r8)    :: press (ncols,kmax)
    REAL(KIND=r8)    :: tin   (ncols,kmax)
    REAL(KIND=r8)    :: qin   (ncols,kmax)
    REAL(KIND=r8)    :: tmst  (ncols,kmax)
    REAL(KIND=r8)    :: qmst  (ncols,kmax)
    REAL(KIND=r8)    :: dtkuo (ncols,kmax)
    REAL(KIND=r8)    :: dqkuo (ncols,kmax)
    REAL(KIND=r8)    :: esat  (ncols)
    REAL(KIND=r8)    :: dtvirt(ncols,kmax)
    REAL(KIND=r8)    :: deltaq(ncols,kmax)
    REAL(KIND=r8)    :: tpar  (ncols)
    REAL(KIND=r8)    :: espar (ncols)
    REAL(KIND=r8)    :: qspar (ncols)
    REAL(KIND=r8)    :: qpar  (ncols)
    REAL(KIND=r8)    :: qex1  (ncols)
    REAL(KIND=r8)    :: tlcl  (ncols)
    REAL(KIND=r8)    :: qexces(ncols)
    REAL(KIND=r8)    :: dqdp  (ncols)
    REAL(KIND=r8)    :: deltap(ncols)
    REAL(KIND=r8)    :: unstab(ncols)
    REAL(KIND=r8)    :: water (ncols)
    REAL(KIND=r8)    :: q1    (ncols)
    REAL(KIND=r8)    :: q2    (ncols)
    REAL(KIND=r8)    :: qsatsm(ncols)
    REAL(KIND=r8)    :: qsum  (ncols)
    REAL(KIND=r8)    :: qsatk (ncols)
    REAL(KIND=r8)    :: x     (ncols)
    REAL(KIND=r8)    :: ubar  (ncols)
    REAL(KIND=r8)    :: b     (ncols)
    REAL(KIND=r8)    :: qeff1 (ncols)
    REAL(KIND=r8)    :: qeff2 (ncols)
    REAL(KIND=r8)    :: pcpwat(ncols)
    REAL(KIND=r8)    :: hnew  (ncols)
    REAL(KIND=r8)    :: slcl  (ncols)
    REAL(KIND=r8)    :: localAcum(ncols)
    LOGICAL          :: llift(ncols)
    LOGICAL          :: lconv(ncols)
    INTEGER          :: ll   (ncols)
    REAL(KIND=r8)    :: sigtop
    INTEGER          :: ksgtop(ncols)
    REAL(KIND=r8)    :: cappa
    REAL(KIND=r8)    :: rdt
    REAL(KIND=r8)    :: cpovl 
    REAL(KIND=r8)    :: tempx(ncols)
    REAL(KIND=r8)    :: coodsl(ncols,kmax)
    REAL(KIND=r8)    :: DeltaPr(ncols,kmax)
    REAL(KIND=r8)    :: ps     (ncols)

    INTEGER :: kthick(ncols)
    INTEGER :: kacum(ncols)
    INTEGER :: i
    INTEGER :: k
    INTEGER :: kk

    !INTEGER, SAVE :: ifp
    !DATA ifp/1/
    cappa=gasr/cp
    DO k=1,kMax
      DO i=1,nCols
         dtdt  (i,k)=0.0_r8
         dqdt  (i,k)=0.0_r8
      END DO
    END DO

    IF(dt .EQ. 0.0e0_r8) RETURN
    !
    !     set default values
    !
    DO i=1,ncols
       kktop(i)=1
       ksgtop(i)=1
       kthick(i)=1
       kacum(i)=1
       kkbot(i)=1
       kuo  (i)=0
       ll(i)=0
       llift(i)=.FALSE.
       lconv(i)=.FALSE.
       ps(i)=prsi(i,1)/1000.0_r8! convert Pa to cb 
    END DO
    DO k=1,kMax
      DO i=1,nCols
          DeltaPr(i,k) = (prsi(i,k) - prsi(i,k+1))/prsi(i,1)
      END DO
    END DO

    !
    !     define kthick in terms of sigma values
    !
    ! thick=0.65e0_r8 ! sthick; upper limit for originating air for lcl.
    ! replaces kthick.

    DO k=1,kmax
       !
       !IF(prsi(i,k)/prsi(i,1) >=  65000.0_r8/100000.0_r8 .and. prsi(i,k+1)/prsi(i,1) <  65000.0_r8/100000.0_r8 )
       !IF(prsi(i,k) >=  65000.0_r8 .and. prsi(i,k+1) <  65000.0_r8 ) THEN
       !
       DO i=1,ncols
          IF(prsi(i,k) >=  65000.0_r8 .and. prsi(i,k+1) <  65000.0_r8 ) THEN
          !IF(si(k).GE.sthick.AND.si(k+1).LT.sthick) THEN
             kthick(i) = k
!             go to 12
          END IF
       END DO
    END DO
!12  CONTINUE

   DO k=1,kmax
       ! sacum=0.46e0_r8 ! sacum; top level for integrated moisture 
                                                          ! convergence test. replaces
                                                          ! kacum
       !
       !
       !
       !IF(prsi(i,k)/prsi(i,1) >=  46000.0_r8/100000.0_r8 .and. prsi(i,k+1)/prsi(i,1) <  46000.0_r8/100000.0_r8 )
       !IF(prsi(i,k) >=  46000.0_r8 .and. prsi(i,k+1) <  46000.0_r8 ) THEN
       DO i=1,ncols
          IF(prsi(i,k) >=  46000.0_r8 .and. prsi(i,k+1) <  46000.0_r8 ) THEN
       !   IF(si(k).GE.sacum.AND.si(k+1).LT.sacum) THEN
             kacum(i) = k
!            go to 14
          END IF
       END DO
    END DO
!14  CONTINUE

    sigtop=0.075_r8
    DO k=1,kmax
       DO i=1,ncols
          !IF(prsi(i,k)/prsi(i,1) >=  7500.0_r8/100000.0_r8 .and. prsi(i,k+1)/prsi(i,1) <  7500.0_r8/100000.0_r8 )
          IF(prsi(i,k) >=  7500.0_r8 .and. prsi(i,k+1) <  7500.0_r8 ) THEN
          !IF(si(k).GE.sigtop.AND.si(k+1).LT.sigtop) THEN
             ksgtop(i) = k
             !go to 16
          END IF
       END DO
    END DO
!16  CONTINUE
    !IF (ifp .EQ. 1) THEN
    !   ifp=0
    !   WRITE (*, '(3(a,i5))') ' kthick= ', kthick, &
    !        ' kacum= ', kacum, ' ksgtop= ', ksgtop
    !ENDIF
    rdt=1.0e0_r8/dt
    cpovl= cp/hl
    kk = kmax
    !
    !     qn is q at the n-1 time level
    !     del=del sigma (del p ovr psfc)
    !     ps=sfc pres (cb)
    !     calculate dq only up to levh - set to zero above levh
    !
    DO k=1,kmax
       DO i=1,ncols
          coodsl(i,k) =  prsl(i,k)/prsi(i,1)
          dq    (i,k) =  qn1(i,k) - qn(i,k)
       END DO
    END DO
    !
    !     net time rate moisture inflow ..lwst kacum.. lyrs
    !     acum is a leapfrog quantity, like dq
    !
    !     acum0=-2.0e-8_r8! acum0; threshold moisture convergence such that 
    !                              integrated moisture
    !                              convergence > - acum0 for convection to occur.

    DO i=1, ncols
       localAcum(i)= acum0
       ! ps.........sfc pres (cb)
       hnew(i) = 1.0e-2_r8*ps(i)
    END DO

    DO k=1,kMax
       DO i=1, ncols
          IF(k>=1 .AND. k <= kacum(i))THEN
             localAcum(i)=localAcum(i)+rdt*dq(i,k)*DeltaPr(i,k)
          END IF
       END DO
    END DO

    DO i=1,ncols
       IF(localAcum(i).LT.0.0_r8) kuo(i) =  2
    END DO

    DO k=1,kmax
       DO i=1,ncols
         ! ps.........sfc pres (cb)
         ! press(i,k)=sl(k)*ps(i)
          press(i,k)=prsl(i,k)/1000.0_r8 !convet Pa to cb
       END DO
    END DO
    !
    !     tin, qin are prelim tmp and q at n+1
    !     zero q if it is  negative
    !
    DO k=1,kmax
       DO i=1,ncols
          qin(i,k) = qn1(i,k)
          IF (qn1(i,k).LE.0.0e0_r8) qin(i,k) = 1.0e-12_r8
       END DO
    END DO

    DO k=1,kmax
       DO i=1,ncols
          tin(i,k) = tn1(i,k)
       END DO
    END DO
    !ki=1     ! ki  ; lowest level from which parcels can be lifted to find lcl
    DO i=1,ncols
       qex1(i) = 0.0e0_r8
       qpar(i) = qin(i,ki)
    END DO
    !
    !     lift parcel from k=ki until it becomes saturated
    !
    DO k =1,kMax
       DO i=1,ncols
          IF(k >= ki .AND. k <= kthick(i))THEN
             IF(.NOT.llift(i)) THEN
                tpar(i) = tin(i,ki)* &
                          EXP(cappa*LOG(press(i,k)/press(i,ki)))
                tempx(i) = tpar(i) - tbase
                IF (tempx(i) >= -50.0_r8) THEN
                   espar(i) = d1 + tempx(i)*(d2 + tempx(i)*(d3 +    tempx(i)*&
                       (d4 + tempx(i)*(d5 + tempx(i)*(d6 + d7*tempx(i))))))
                ELSE
                   espar(i)=0.00636e0_r8*EXP(25.6e0_r8*(tempx(i)+50.e0_r8)/(tbase-50.e0_r8))
                END IF
                qspar(i) =rmwmd*espar(i)/ &
                     (press(i,k)-(1.0e0_r8-rmwmd)*espar(i))
                qexces(i) = qpar(i) - qspar(i)
                !
                !     if parcel not saturated,  try next level
                !
                IF (qexces(i).LT.0.0e0_r8) THEN
                   qex1(i) = qexces(i)
                   !
                   !     saturated - go down and find p,t, sl at the lcl;
                   !     if sat exists in first layer (k=ki), use this as lcl
                   !
                ELSE IF (k.EQ.ki) THEN
                   plcl(i) = press(i,k)
                   tlcl(i) = tpar(i)
                   tlcl(i) = tlcl(i) + 1.0e0_r8
                   slcl(i) = plcl(i)/ps(i)
                   ll(i)   = k
                   kkbot(i) = ll(i)
                   llift(i) = .TRUE.
                ELSE
                   dqdp(i) = (qexces(i)-qex1(i))/ &
                        (press(i,k-1)-press(i,k))
                   deltap(i)= qexces(i)/dqdp(i)
                   plcl(i)=press(i,k) + deltap(i)
                   tlcl(i)=  tpar(i) * (1.0e0_r8+2.0e0_r8*cappa*deltap(i) &
                        /(press(i,k)+press(i,k-1)))
                   tlcl(i) = tlcl(i) + 1.0e0_r8
                   slcl(i) = plcl(i)/ps(i)
                   ll(i)   = k
                   kkbot(i) = ll(i)
                   llift(i) = .TRUE.
                   !
                   !     give parcel a one-degree goose
                   !
                ENDIF
                !
                !     lifting cond level found - get out of loop
                !
             END IF
          END IF
       END DO
    END DO
    !
    !     quit if parcel still unsat at k = kthick - - -
    !     in this case,set low value of plcl as signal to shalmano
    !
    DO i=1,ncols
       IF(.NOT.llift(i)) THEN
          plcl(i) = 1.0e0_r8
          kuo  (i) = 5
       ENDIF
    END DO

    CALL mstad2(hnew  ,coodsl,tin   ,tmst  ,qmst  ,kktop , &
         slcl  ,ll    ,qin   ,tlcl  ,llift ,ncols  ,kmax)
    !
    !     tmst and qmst,k=1...ktop contain cloud temps and specific humid.
    !     store values of plcl and ktop to pass to shalmano
    !

    lconv = llift

    !
    !     test 2...thickness of unstable region must exceed 300 mb or 30000 Pa
    !
    DO i=1,ncols
       IF(lconv(i)) THEN
           unstab(i) = prsl(i,ll(i)) - prsl(i,kktop(i)) !  /prsi(i,1)
          !unstab(i) = sl(ll(i)) - sl(kktop(i))
          IF (unstab(i).LT. 30000.0_r8 )kuo(i)=7
          !IF (unstab(i).LT.0.3e0_r8)kuo(i)=7
          !
          !     here  kuo requirements are met, first compute water=tons water
          !     subst accum in zero leap-frog time step per sq m in lyrs 1-k
          !     avl  water  (water=water*g/ps).gt.zero)...
          !
          water(i)=0.0e0_r8
       ENDIF
    END DO

    DO k=1,kmax
       DO i=1,ncols
          IF(lconv(i).AND.k.LE.kktop(i)) THEN
             water(i)=water(i)+  dq(i,k)*DeltaPr(i,k)
          ENDIF
       END DO
    END DO

    DO i=1,ncols
       IF(lconv(i)) THEN
          IF(water(i).LE.0.0e0_r8) kuo(i)=8
          IF(kuo(i).GT.0) lconv(i)=.FALSE.
       ENDIF
    END DO

    DO i=1,ncols
       IF(lconv(i)) THEN
          q1(i)=0.0e0_r8
          q2(i)=0.0e0_r8
          qsatsm(i)=0.0e0_r8
          qsum(i)  =0.0e0_r8
       ENDIF
    END DO
    !
    !     calculate four vertical averages - sat deficit of environment,
    !     virt temp excess of cloud, and  q, qsat for environment ---
    !
    DO k=1,kmax
       DO i=1,ncols
          IF(lconv(i).AND.k.GE.ll(i).AND.k.LE.kktop(i)) THEN
             tempx(i) = tin(i,k) - tbase
             IF (tempx(i) >= -50.0_r8) THEN
                esat(i) = d1 + tempx(i)*(d2 + tempx(i)*(d3 +    tempx(i)*&
                    (d4 + tempx(i)*(d5 + tempx(i)*(d6 + d7*tempx(i))))))
             ELSE
                esat(i)=0.00636e0_r8*EXP(25.6e0_r8*(tempx(i)+50.e0_r8)/(tbase-50.e0_r8))
             END IF
             qsatk(i)=rmwmd*esat(i)/(press(i,k)-(1.0e0_r8-rmwmd)*esat(i))
             x(i) = qsatk(i) - qin(i,k)
             deltaq(i,k) = x(i)
             q1(i) = q1(i) + x(i)*DeltaPr(i,k)
             qsum(i) = qsum(i) + qin(i,k)*DeltaPr(i,k)
             qsatsm(i) = qsatsm(i) + qsatk(i)*DeltaPr(i,k)
             x(i)=tmst(i,k)-tin(i,k)+0.61e0_r8*tin(i,k)*(qmst(i,k)-qin(i,k))
             dtvirt(i,k) = x(i)
             q2(i)=q2(i)+x(i)*DeltaPr(i,k)
          ENDIF
       END DO
    END DO

    DO i=1,ncols
       IF(lconv(i)) THEN
          q2(i)=q2(i)*cpovl
          IF (q1(i).LE.0.0e0_r8) q1(i) = 1.0e-9_r8
          IF (q2(i).LE.0.0e0_r8) q2(i) = 1.0e-9_r8
          ubar(i) = qsum(i)/qsatsm(i)
          IF (ubar(i).GE.1.0e0_r8) ubar(i) =0.999e0_r8
          b(i) = 1.0e0_r8-ubar(i)
          IF (b(i).GT.1.0e0_r8) b(i)=1.0e0_r8
          qeff1(i) = water(i) * b(i)/q1(i)
          qeff2(i) = water(i) *(1.0e0_r8-b(i))/q2(i)
          IF(qeff1(i).LT.0.002e0_r8) lconv(i)=.FALSE.
          qeff1(i)=MIN(qeff1(i),1.0e0_r8)
          qeff2(i)=MIN(qeff2(i),1.0e0_r8)
       ENDIF
    END DO
    !
    ! exclude convective clouds from top of model
    !
    DO i=1,ncols
       ! if (kktop(i) .eq. kmax) kktop(i)=kmax-1
       kktop(i)=MIN(kktop(i),ksgtop(i))
    END DO

    DO k=1,kmax
       DO i=1,ncols
          IF(lconv(i).AND.k.LE.kktop(i)) THEN
             IF(k.LT.ll(i)) THEN
                dtkuo(i,k) = 0.0e0_r8
                dqkuo(i,k) = 0.0e0_r8
             ELSE
                dqkuo(i,k) = qeff1(i) * deltaq(i,k)
                dtkuo(i,k) = qeff2(i) * dtvirt(i,k)
             ENDIF
             !
             !     start loop at 1 to record loss of dq in diagnostics for k lt. ll
             !
             tin(i,k)=tin(i,k)+dtkuo(i,k)
             qin(i,k)=dqkuo(i,k)+ qin(i,k)
          ENDIF
       END DO
    END DO
    !
    !     calculate convective precipitation from latent heating
    !
    DO i=1,ncols
       IF(lconv(i)) THEN
          pcpwat(i) = qeff2(i)*q2(i)
          !
          !     set aside convective precip in separate array for diagnostics
          !
          geshem(i) = geshem(i)+ps(i)*pcpwat(i)/(grav*2.0e0_r8)
       ENDIF
    END DO
    !
    !     remove all water vapor that was set aside in 'water' sum---
    !
    DO k=1,kmax
       DO i=1,ncols
          IF(lconv(i).AND.k.LE.kktop(i)) THEN
             qin(i,k) = qin(i,k) -dq(i,k)
          ENDIF
       END DO
    END DO
    !
    !     restore basic fields
    !
    DO k=1,kmax
       DO i=1,ncols
          IF(lconv(i).AND.k.LE.kktop(i)) THEN
             dtdt (i,k)=(tin(i,k)-tn1(i,k))/dt
             tn1(i,k)=tin(i,k)
             dqdt (i,k)=(qin(i,k)-qn1(i,k))/dt
             qn1(i,k)=qin(i,k)
          ENDIF
       END DO
    END DO

    DO i=1,ncols
       IF(lconv(i)) THEN
          kuo(i)=1
       ENDIF
    END DO
  END SUBROUTINE RunCu_Kuolcl
 
   
  ! mstad2 :this subroutine lifts the parcel from the lcl and computes the
  !         parcel temperature at each level; temperature is retrieved from
  !         lookup tables defined in the first call.
  !         This version-mstad2-lifts parcel from lifting condensation
  !         level, where pres is slcl and temp is tlcl
  !         level ll is the first regular level above the lcl
  !         this version contains double resolution tables
  !         and uses virtual temp in the buoyancy test
  !         new sat vapor pressure used
  




  SUBROUTINE mstad2(ps, coodsl, tin, tmst, qmst, ktop,  &
       slcl, ll, qin, tlcl, llift, ncols, kmax)
       
    !
    ! this routine accepts input data of
    !         ps                  = surface pressure dvded by 100 cbs
    !         coodsl(k=1,--,km)   = factor such that p at lvl k = ps*coodsl(k)
    !         tin(1 + (k-1)*it ) = input temps in a column
    !         km                 = number of levels 
    ! and returns                
    !         the = equiv pot temp calc from parcel p and t at its lcl
    !         ktop (g.t.e 1 and l.t.e. km) is highest lvl for which
    !              tin is colder than moist adiabat given by the
    !              (ktop=1 denotes tin(k=2) is already g.t.e. moist adb.)
    !                   allowance is made for perhaps one level below ktop
    !                  at which tin was warmer than tmst. 
    !         tmst(k) and qmst(k) for k=2,--,ktop, are temp and sat spec
    !         hums on moist adb the.
    !         pressure in layer one must lie between 50 and 110 cbs.
    !         temp in layer one must lie between 220 and 330 degrees
    !         the resulting the must lie between 220 and 500 degrees.
    !         ( the is tested for this possibility, an error return of
    !         ier=0 denoting okeh conditions, a return of ier=1 denoting
    !         violation of this range.)
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kmax......Number of sigma levels
    ! cp........Specific heat of air           (j/kg/k)
    ! hl........heat of evaporation of water     (j/kg)
    ! gasr......gas constant of dry air        (j/kg/k)
    ! tbase.....temperature of fusion of ice
    ! tmst......temperature on moist adb
    ! qmst......spec hums on moist adb
    ! qin.......spec hums in a column
    ! slcl......pressure of  lifting condensation level (in sigma)
    ! tlcl......temp of  lifting condensation level
    ! ll........is the first regular level above the lcl
    ! llift
    !==========================================================================
    IMPLICIT NONE

    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: coodsl(ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tin  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: tmst (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: qmst (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: qin  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: ps    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: slcl (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tlcl (ncols)
    INTEGER, INTENT(inout) :: ktop (ncols)
    INTEGER, INTENT(in   ) :: ll   (ncols)
    LOGICAL, INTENT(in   ) :: llift(ncols)
    !
    !  local arrays (temporary storage)
    !
    REAL(KIND=r8)          :: pp    (ncols)
    REAL(KIND=r8)          :: ti    (ncols)
    INTEGER       :: jt    (ncols)
    REAL(KIND=r8)          :: x     (ncols)
    REAL(KIND=r8)          :: pk    (ncols)
    INTEGER       :: kp    (ncols)
    REAL(KIND=r8)          :: y     (ncols)
    REAL(KIND=r8)          :: yy    (ncols)
    REAL(KIND=r8)          :: xx    (ncols)
    REAL(KIND=r8)          :: the   (ncols)
    REAL(KIND=r8)          :: tk    (ncols)
    INTEGER       :: kt    (ncols)
    REAL(KIND=r8)          :: pi    (ncols)
    REAL(KIND=r8)          :: qlcl  (ncols)
    INTEGER       :: ip    (ncols)
    INTEGER       :: lstb  (ncols)
    INTEGER       :: lcld  (ncols)
    REAL(KIND=r8)          :: tvdiff(ncols)
    INTEGER       :: ktop1 (ncols)


    INTEGER       :: i
    INTEGER       :: k
    CHARACTER(LEN=*), PARAMETER :: h="**(mstad2)**"
    CHARACTER(LEN=256) :: line

    !
    !     surface press is ps, pressure at other lvls is ps times coodsl(i,k)
    !     compute theta e  = the  , for layer one
    !

    ! Get lower index and interpolation weight for:
    !   Tlcl in the Theta_e(T,P) table.

    DO i=1,ncols
       IF(llift(i)) THEN
          ! get position in the matrix
          ti(i)=MAX((tlcl(i)-the_temp0),0.0_r8)/the_dtemp + 1
          ! round to smaller int => get lower index
          jt(i)=INT(ti(i))
          ! difference between int and real values is the 
          ! interpolation weight
          x(i)=ti(i)-jt(i)
          xx(i)=1.0e0_r8-x(i)

          ! check for temperature range
          IF(jt(i).LT.1) THEN
             WRITE(line,'(g12.4," < ",g12.4)') ti(i),the_temp0
             CALL FatalError(h//"Theta_e table: temp="//line//"=low limit")
          ENDIF
          IF(jt(i).GE.the_ntemp) THEN
             WRITE(line,'(g12.4," > ",g12.4)') ti(i),the_temp0+(the_ntemp-1)*the_dtemp
             CALL FatalError(h//"Theta_e table: temp"//line//"high limit")
          ENDIF
       ENDIF
    END DO

    ! Get lower index and interpolation weight for:
    !   Plcl in the Theta_e(T,P) table.

    DO i=1,ncols
       IF(llift(i)) THEN
          ! Transform sigma_lcl to pressure_lcl
          pp(i)=ps(i)*slcl(i)
          ! get position in the matrix
          ! This should be:
          !    pk(i)=(pp(i)-the_pres0)/the_dpres + 1
          ! But gives diferent numerical result than this: 
          pk(i)=pp(i)*the_dpinv - the_pres0*the_dpinv + 1
          ! round to smaller int => get lower index
          kp(i)=INT(pk(i))
          ! difference between int and real values is the 
          ! interpolation weight
          y(i)=pk(i)-kp(i)
          yy(i)=1.0e0_r8-y(i)

          ! check for pressure range
          IF(kp(i).LT.1) THEN
             WRITE(line,'(g12.4," < ",g12.4)') pp(i),the_pres0
             CALL FatalError(h//"Theta_e table: pres="//line//"=low limit")
          ENDIF
          IF(kp(i).GE.the_npres) THEN
             WRITE(line,'(g12.4," > ",g12.4)') pp(i),the_pres0+(the_npres-1)*the_dpres
             CALL FatalError(h//"Theta_e table: pres"//line//"high limit")
          ENDIF
       ENDIF
    END DO

    !
    ! INTERPOLATE to find Theta_e at LCL
    !

    DO i=1,ncols
       IF(llift(i)) THEN
          the(i)=xx(i)*(yy(i)*thetae(jt(i)  ,kp(i)  )+   &
                         y(i)*thetae(jt(i)  ,kp(i)+1)  ) &
                 +x(i)*(yy(i)*thetae(jt(i)+1,kp(i)  )+   &
                         y(i)*thetae(jt(i)+1,kp(i)+1)  )
       ENDIF
    END DO

    !
    ! Get t and q on mst adiabat (tmst and qmst) for layers 
    !   which are colder than the mst adiabat
    !

    ! Get lower index and interpolation weight for:
    !   Theta_e_lcl in the tfmthe(The,P) and qfmthe(The,P) tables

    DO i=1,ncols
       IF(llift(i)) THEN
          ! get position in the matrix
          tk(i)=(the(i)-fm_the0)/fm_dthe + 1
          ! round to smaller int => get lower index
          kt(i)=INT(tk(i))
          ! difference between int and real values is the 
          ! interpolation weight
          y(i)=tk(i)-kt(i)
          yy(i)=1.0e0_r8-y(i)

          ! check for theta_e range
          IF(kt(i).LT.1) CALL FatalError(h//&
               "tfmthe_e table: the < low limit")
          IF(kt(i).GE.fm_nthe) CALL FatalError(h//&
               "tfmthe_e table: the > high limit")
       ENDIF
    END DO

    ! Get lower index and interpolation weight for:
    !   Plcl in the tfmthe(The,P) and qfmthe(The,P) tables

    DO i=1,ncols
       IF(llift(i)) THEN
          ! get position in the matrix
          ! This should be:
          !    pi(i)=(pp(i)-fm_pres0)/fm_dpres + 1
          ! But gives diferent numerical result than this: 
          pi(i)=pp(i)*fm_dpinv - fm_pres0*fm_dpinv + 1
          ! round to smaller int => get lower index
          ip(i)=INT(pi(i))
          ! difference between int and real values is the 
          ! interpolation weight
          x(i)=pi(i)-ip(i)
          xx(i)=1.0_r8-x(i)

          ! check for pressure range
          IF(ip(i).LT.1) CALL FatalError(h//&
               "tfmthe_e table: pres < low limit")
          IF(ip(i).GE.fm_npres) CALL FatalError(h//&
               "tfmthe_e table: pres > high limit")
       ENDIF
    END DO

    !
    ! INTERPOLATE to get qlcl
    !
    DO i=1,ncols
       IF(llift(i)) THEN
          qlcl(i)=(1.0e0_r8-x(i))*(yy(i)*qfmthe(kt(i),ip(i)) &
               +y(i)*qfmthe(kt(i)+1,ip(i)) ) &
               +x(i)*(yy(i)*qfmthe(kt(i),ip(i)+1) &
               +y(i)*qfmthe(kt(i)+1,ip(i)+1) )
       ENDIF
    END DO

    !
    ! Initialize qmst and tmst to their LCL values
    !
    DO k = 1,kmax
       DO i=1,ncols
          tmst(i,k) = tlcl(i)
          qmst(i,k) = qlcl(i)
       END DO
    END DO

    !
    !     we will allow one stable layer (with tin > tmst) to
    !     interrupt a sequence of unstable layers.
    !
    DO i=1,ncols
       IF(llift(i)) THEN
          ktop(i)=ll(i)-1
          IF (ktop(i) <= 0) ktop(i)=1
          lstb(i)=0
          lcld(i)=0
       ENDIF
    END DO

    ! Start from the ground and go upward
    ! For each layer, test if 
    DO k=1,kmax
       DO i=1,ncols
          IF(llift(i).AND.lcld(i).EQ.0.AND.k.GE.ll(i)) THEN
             ! Get layer pressure
             pp(i)=ps(i)*coodsl(i,k)
             ! Get position in the matrix
             ! This should be:
             !    pi(i)=(pp(i)-fm_pres0)/fm_dpres + 1
             ! But gives diferent numerical result than this: 
             pi(i)=pp(i)*fm_dpinv - fm_pres0*fm_dpinv + 1
             ! round to smaller int => get lower index
             ip(i)=INT(pi(i))
             ! difference between int and real values is the 
             ! interpolation weight
             x(i)=pi(i)-ip(i)

             ! INTERPOLATE
             tmst(i,k)=(1.0e0_r8-x(i))*(yy(i)*tfmthe(kt(i),ip(i)) &
                  +y(i)*tfmthe(kt(i)+1,ip(i))) &
                  +x(i)*(yy(i)*tfmthe(kt(i),ip(i)+1) &
                  +y(i)*tfmthe(kt(i)+1,ip(i)+1))
             qmst(i,k)=(1.0e0_r8-x(i))*(yy(i)*qfmthe(kt(i),ip(i)) &
                  +y(i)*qfmthe(kt(i)+1,ip(i))) &
                  +x(i)*(yy(i)*qfmthe(kt(i),ip(i)+1) &
                  +y(i)*qfmthe(kt(i)+1,ip(i)+1))
             !
             !     buoyancy test with virtual temp correction
             !
             tvdiff(i) = tmst(i,k)-tin(i,k)  &
                  + 0.61e0_r8 *tin(i,k)*(qmst(i,k) - qin(i,k))
             IF(tvdiff(i).LE.0.0_r8) THEN
                IF(lstb(i).EQ.0) THEN
                   lstb(i) = 1
                   ktop1(i) = ktop(i)
                   ktop(i)  = ktop(i) + 1
                ELSE
                   lcld(i) = 1
                   IF((ktop(i)-ktop1(i)-1).LE.0) THEN
                      ktop(i) = ktop1(i)
                   ENDIF
                ENDIF
             ELSE
                ktop(i)=ktop(i)+1
             ENDIF

          ENDIF
       END DO
    END DO

  END SUBROUTINE mstad2

  SUBROUTINE InitMstad2()
    REAL(KIND=r8)          :: t
    REAL(KIND=r8)          :: esat
    REAL(KIND=r8)          :: kappa
    REAL(KIND=r8)          :: eps
    REAL(KIND=r8)          :: el
    REAL(KIND=r8)          :: p
    REAL(KIND=r8)          :: ratio
    REAL(KIND=r8)          :: power
    REAL(KIND=r8)          :: pd
    REAL(KIND=r8)          :: pdkap
    REAL(KIND=r8)          :: thee
    REAL(KIND=r8)          :: fun
    REAL(KIND=r8)          :: dfun
    REAL(KIND=r8)          :: chg
    INTEGER       :: i
    INTEGER       :: k
    INTEGER       :: l
    REAL(KIND=r8), PARAMETER :: rv=461.5e0_r8
    REAL(KIND=r8), PARAMETER :: cl=4187.e0_r8

    ! Convergence limit for Newton method
    REAL(KIND=r8) :: crit
    crit=2.0_r8**(-23.0_r8)

    ! set up tables of theta e of p and t, and of t of theta e and p

    kappa=gasr/cp
    eps=gasr/rv

    ! Construct table to get theta_e from temperature and pressure
    ! The table is thetae(i,k) where:
    !
    ! t= the_temp0, temp0+dtemp, ..., temp0+dtemp*(ntemp-1)
    ! p= the_pres0, pres0+dpres, ..., pres0+dpres*(npres-1)
    !
    ! These are parameters defined in module preamble

    !$OMP PARALLEL DO PRIVATE(t,el,esat,k,p,ratio,power,pd,pdkap)
    DO i=1,the_ntemp
       t=the_temp0+real(i-1,r8)*the_dtemp

       el=hl-(cl-cp)*(t-tbase)
       esat = es(t)

       p=the_pres0
       DO k=1,the_npres

          ratio=eps*esat/(1.e2_r8*p-esat)
          power=el*ratio/(cp*t)
          pd=p-1.0e-2_r8*esat
          pdkap=kappa*LOG(pd)
          pdkap=EXP(pdkap)
          thetae(i,k)=t* EXP (power)/(pdkap)

          p=p+the_dpres
       END DO
    END DO 
    !$OMP END PARALLEL DO

    ! Construct table to get temperature and saturation vapor
    ! pressure from thetae and p. The tables are tfmthe(i,k), 
    ! and qfmthe(i,k) where:
    ! 
    ! t_i= fm_temp0, temp0+dtemp, ..., temp0+dtemp*(ntemp-1)
    ! p_k= fm_pres0, pres0+dpres, ..., pres0+dpres*(npres-1)

    !$OMP PARALLEL DO PRIVATE(p,t,l,ess,ratio,el,power,pdkap,fun,dfun,chg,thee)
    DO i=1,fm_nthe

       thee=fm_the0+REAL(i-1,r8)*fm_dthe

       ! Set values for p=0.0
       !   tfmthe(i,1) is identically 0.0e0
       !   qfmthe(i,1) is identically 0.0e0

       tfmthe(i,1) = 0.0e0_r8
       qfmthe(i,1) = 0.0e0_r8

       ! Now set values for p > 0.0

       p=fm_dpres
       DO k=2,fm_npres

          ! first guess

          IF (p < 0.025_r8) THEN
             t=100.0_r8
          ELSE IF (p < 0.05_r8) THEN
             t=tbase
          ELSE
             t=300.0_r8
          END IF

          ! newton iteration method

          DO l=1,100
             !
             !     sat vapor pressure
             !
             ess = es(t)
             ratio=eps*ess/(1.e2_r8*p-ess)
             el=hl-(cl-cp)*(t-tbase)
             power=ratio*el/(cp*t)
             pdkap=p-1.0e-2_r8*ess
             pdkap=kappa*LOG(pdkap)
             pdkap=EXP(pdkap)
             fun=t* EXP (power)/pdkap
             dfun=(fun/t)*(1.0e0_r8+(ratio/(cp*t) )*((cp-cl)*t &
                  +(p/(p-1.0e-2_r8*ess))*(el*el/(rv*t))))
             chg=(thee-fun)/dfun
             t=t+chg
             IF( ABS(chg) .LT. crit) EXIT
          END DO
          tfmthe(i,k)=t
          !
          !     sat vapor pressure
          !
          ess = es(t)
          qfmthe(i,k)=eps*ess/(1.e2_r8*p-(1.0e0_r8-eps)*ess)
          p=p+fm_dpres
       END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE InitMstad2

END MODULE Cu_Kuolcl
