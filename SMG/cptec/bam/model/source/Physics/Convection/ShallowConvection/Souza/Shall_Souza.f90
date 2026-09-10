!
!  $Author: panetta $
!  $Date: 2007/07/20 13:58:43 $
!  $Revision: 1.7 $
!
MODULE Shall_Souza

  USE Constants, ONLY :  &
       cp                 , &
       hl                 , &
       gasr               , &
       grav               , &
       rmwmd              , &
       rmwmdi             , &
       e0c                , &
       delq               , &
       r8                 , &
       i8
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
  !!PUBLIC :: InitShall_Souza
  PUBLIC :: shallsouza

  REAL(KIND=r8) :: aa(15)
  REAL(KIND=r8) :: ad(15)
  REAL(KIND=r8) :: ac(15)
  REAL(KIND=r8) :: actop
  REAL(KIND=r8) :: thetae(151,181)
  REAL(KIND=r8) :: tfmthe(431,241)
  REAL(KIND=r8) :: qfmthe(431,241)
  REAL(KIND=r8) :: ess
  INTEGER :: kbase
  INTEGER :: kcr
  REAL(KIND=r8) , ALLOCATABLE :: dels  (:)
  REAL(KIND=r8) , ALLOCATABLE :: gams  (:)
  REAL(KIND=r8) , ALLOCATABLE :: gammod(:)
  REAL(KIND=r8) , ALLOCATABLE :: delmod(:)
  REAL(KIND=r8) :: rlocp
  REAL(KIND=r8) :: rgrav
  REAL(KIND=r8) :: rlrv
  REAL(KIND=r8)  :: const1
  REAL(KIND=r8)  :: const2
  REAL(KIND=r8)  :: xx1

  REAL(KIND=r8) , PARAMETER :: xkapa=0.2857143_r8


CONTAINS


  SUBROUTINE shallsouza(te,qe,ps,sl,tfz,qfz,dt,nCols,kMax,kuo,  &
       noshal, klcl,ktop,par6,par7)
    ! Shallow cumulus heating and moistening tendencies
    ! Enio Pereira de Souza 12/Jul/2001
    ! modified and adapted to GCM by Silvio Nilo Figueroa Jun-2004
    !
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8) ,    INTENT(INOUT) :: te    (nCols,kMax)  ! air temperature (K)
    REAL(KIND=r8) ,    INTENT(INOUT) :: qe    (nCols,kMax)  ! specific humidity (kg/kg)
    INTEGER, INTENT(INOUT) :: klcl  (nCols)
    INTEGER, INTENT(INOUT) :: ktop  (nCols)
    INTEGER, INTENT(INOUT) :: noshal(nCols)
    REAL(KIND=r8) ,    INTENT(IN   ) :: ps    (nCols)       ! surface pressure (cb)
    REAL(KIND=r8) ,    INTENT(IN   ) :: sl    (kMax)       ! sigma layers
    REAL(KIND=r8) ,    INTENT(IN   ) :: tfz   (nCols)       ! sensible heat flux (W/m^2)
    REAL(KIND=r8) ,    INTENT(IN   ) :: qfz   (nCols)       ! latent heat flux (W/m^2)
    REAL(KIND=r8) ,    INTENT(IN   ) :: dt                 ! time step (s)
    REAL(KIND=r8) ,    INTENT(IN   ) :: par6
    REAL(KIND=r8) ,    INTENT(IN   ) :: par7
    !
    INTEGER :: kzi  (nCols)
    INTEGER :: kuo  (nCols)
    LOGICAL :: llift(nCols)
    LOGICAL :: lconv(nCols)
    !
    REAL(KIND=r8)  :: press   (nCols,kMax)
    REAL(KIND=r8)  :: tv      (nCols,kMax)
    REAL(KIND=r8)  :: ze      (nCols,kMax)
    REAL(KIND=r8)  :: den     (nCols,kMax)
    REAL(KIND=r8)  :: delz    (nCols,kMax)
    REAL(KIND=r8)  :: pi      (nCols,kMax)
    REAL(KIND=r8)  :: th      (nCols,kMax)
    REAL(KIND=r8)  :: thv     (nCols,kMax)
    REAL(KIND=r8)  :: es      (nCols,kMax)
    REAL(KIND=r8)  :: qes     (nCols,kMax)
    REAL(KIND=r8)  :: se      (nCols,kMax)
    REAL(KIND=r8)  :: uhe     (nCols,kMax)
    REAL(KIND=r8)  :: uhes    (nCols,kMax)
    REAL(KIND=r8)  :: uhc     (nCols,kMax)
    REAL(KIND=r8)  :: dtdt    (nCols,kMax)
    REAL(KIND=r8)  :: dqdt    (nCols,kMax)
    REAL(KIND=r8)  :: gamma   (nCols,kMax)
    REAL(KIND=r8)  :: dlamb   (nCols,kMax)
    REAL(KIND=r8)  :: dldzby2 (nCols,kMax)
    REAL(KIND=r8)  :: sc      (nCols,kMax)
    REAL(KIND=r8)  :: sc0     (nCols,kMax)
    REAL(KIND=r8)  :: qc      (nCols,kMax)
    REAL(KIND=r8)  :: ql      (nCols,kMax)
    REAL(KIND=r8)  :: tvm     (nCols,kMax)
    REAL(KIND=r8)  :: scv     (nCols,kMax)
    REAL(KIND=r8)  :: sev     (nCols,kMax)
    REAL(KIND=r8)  :: sc0v    (nCols,kMax)
    REAL(KIND=r8)  :: scvm    (nCols,kMax)
    REAL(KIND=r8)  :: sevm    (nCols,kMax)
    REAL(KIND=r8)  :: sc0vm   (nCols,kMax)
    REAL(KIND=r8)  :: buoy1   (nCols)
    REAL(KIND=r8)  :: buoy2   (nCols)
    REAL(KIND=r8)  :: cape    (nCols)
    REAL(KIND=r8)  :: tcape   (nCols)
    REAL(KIND=r8)  :: emp     (nCols,kMax)
    REAL(KIND=r8)  :: efic    (nCols)
    REAL(KIND=r8)  :: fin     (nCols)
    REAL(KIND=r8)  :: tcold   (nCols)
    REAL(KIND=r8)  :: thot    (nCols)
    REAL(KIND=r8)  :: sigwshb (nCols)
    REAL(KIND=r8)  :: wc      (nCols,kMax)
    REAL(KIND=r8)  :: wssc    (nCols,kMax)
    REAL(KIND=r8)  :: wqsc    (nCols,kMax)
    REAL(KIND=r8)  :: dsdt    (nCols,kMax)
    REAL(KIND=r8)  :: tpar    (nCols)
    REAL(KIND=r8)  :: qex1    (nCols)
    REAL(KIND=r8)  :: qpar    (nCols)
    REAL(KIND=r8)  :: espar   (nCols)
    REAL(KIND=r8)  :: qspar   (nCols)
    REAL(KIND=r8)  :: qexces  (nCols)
    REAL(KIND=r8)  :: dqdp    (nCols)
    REAL(KIND=r8)  :: deltap  (nCols)
    REAL(KIND=r8)  :: plcl    (nCols)
    REAL(KIND=r8)  :: tlcl    (nCols)
    INTEGER    :: icount(nCols)

    !snf
    !REAL :: rho
    !REAL :: dp
    !--------
    !
    ! constants
    !
    !
    REAL(KIND=r8) , PARAMETER :: grav=9.806_r8
    REAL(KIND=r8) , PARAMETER :: vlat=2.5e6_r8
    REAL(KIND=r8) , PARAMETER :: cp=1004.5_r8
    REAL(KIND=r8) , PARAMETER :: rd=287.0_r8
    REAL(KIND=r8) , PARAMETER :: rcp=0.286_r8
    REAL(KIND=r8) , PARAMETER :: p00=100000.0_r8
    REAL(KIND=r8) , PARAMETER :: gammad=0.00976_r8
    REAL(KIND=r8) , PARAMETER :: es00=611.2_r8
    REAL(KIND=r8) , PARAMETER :: epslon=0.622_r8
    REAL(KIND=r8) , PARAMETER :: ummeps=0.378_r8
    REAL(KIND=r8) , PARAMETER :: ta0=273.15_r8
    REAL(KIND=r8) , PARAMETER :: co1=21709759.15_r8
    REAL(KIND=r8) , PARAMETER :: co2=29.65_r8
    REAL(KIND=r8) , PARAMETER :: co3=17.67_r8
    REAL(KIND=r8) , PARAMETER :: c0=0.0_r8
    REAL(KIND=r8) , PARAMETER :: dlamb0=1.0e-6_r8

    REAL(KIND=r8)             :: zref
    REAL(KIND=r8) , PARAMETER :: dthv=2.0_r8
    REAL(KIND=r8) , PARAMETER :: fifty=50.0_r8
    !
    INTEGER    :: i
    INTEGER    :: k
    INTEGER    :: ki

    ! 300 400  560 570 800  was tested
    zref=par6
    !
    DO i=1,nCols
       klcl    (i)=1
       ktop    (i)=1
       noshal  (i)=1
       tcape   (i)=0.0_r8
       sigwshb (i)=0.0_r8
       cape    (i)=0.0_r8
       efic    (i)=0.0_r8
       buoy2   (i)=0.0_r8
       buoy1   (i)=0.0_r8
    END DO

    DO k=1,kMax
       DO i=1,nCols
          dtdt (i,k)=0.0_r8
          dqdt (i,k)=0.0_r8
          dsdt (i,k)=0.0_r8
          ze   (i,k)=0.0_r8
          ql   (i,k)=0.0_r8
          press(i,k)=0.0_r8
          tv   (i,k)=0.0_r8
          th   (i,k)=0.0_r8
       END DO
    END DO
    !
    ! begining of a long loop in index " i "
    !
    !
    ! constructing height profile
    !
    DO k=1,kMax
       DO i=1,nCols
          press(i,k)=ps(i)*sl(k)
          press(i,k)=press(i,k)*1000.0_r8
          tv(i,k)=te(i,k)*(1.0_r8+0.608_r8*qe(i,k))
          th(i,k)=te(i,k)*EXP(rcp*LOG(p00/press(i,k)))
          pi(i,k)=cp*EXP(rcp*LOG(press(i,k)/p00))
          thv(i,k)=th(i,k)*(1.0_r8+0.608_r8*qe(i,k))
          es(i,k)=es00*EXP(co3*(te(i,k)-ta0)/(te(i,k)-co2))
          qes(i,k)=epslon*es(i,k)/(press(i,k)-ummeps*es(i,k))
          !snf added
          IF(qes(i,k) .LE. 1.0e-08_r8) qes(i,k)=1.0e-08_r8
          !-----
          den(i,k)=press(i,k)/(rd*tv(i,k))
       END DO
    END DO

    DO i=1,nCols
       ze(i,1)=29.25_r8*te(i,1)*(101.3_r8-ps(i))/ps(i)
       ze(i,1)=MAX(0.0_r8,ze(i,1))
    END DO

    DO k=2,kMax
       DO i=1,nCols
          delz(i,k)=0.5_r8*rd*(tv(i,k-1)+tv(i,k))* &
               LOG(press(i,k-1)/press(i,k))/grav
          ze(i,k)=ze(i,k-1)+delz(i,k)
       END DO
    END DO
    !!snf added
    !-----------------
    !        DO k=1,kMax
    !        DO i=1,nCols
    !            rho=press(i,k)/(287.0_r8*tv(i,k))
    !            dp=ps(i)*del(k)*1000.0_r8
    !            delz(i,k)=dp/(9.81_r8*rho)
    !        enddo
    !        enddo
    !
    !
    !        DO i=1,nCols
    !            ze(i,1)=delz(i,1)/2.0_r8
    !            ze(i,1)=MAX(0.0_r8,ze(i,1))
    !        dO k=2,kMax
    !            ze(i,k)=ze(i,k-1)+delz(i,k)
    !        END DO
    !        END DO
    !------------
    ki=1
    !
    DO i=1,nCols
       llift(i)=.FALSE.
       qex1(i) = 0.0_r8
       qpar(i) = qe(i,ki)
    END DO
    !
    !
    !     lift parcel from k=ki until it becomes saturated
    !
    DO k = ki, kMax
       DO i = 1, nCols
          IF(kuo(i) /= 1 .or. tfz(i) >= 0.0_r8)THEN
             IF(.NOT.llift(i)) THEN
                tpar(i) = te(i,ki)* &
                     EXP(rcp*LOG(press(i,k)/press(i,ki)))
                espar(i) = es2(tpar(i))
                qspar(i) = epslon*espar(i)*1000.0_r8/ &
                     (press(i,k)-(1-epslon)*espar(i)*1000.0_r8)
                qexces(i) = qpar(i) - qspar(i)
                !
                !     if parcel not saturated,  try next level
                !
                IF (qexces(i).LT.0.0_r8) THEN
                   qex1(i) = qexces(i)
                   !
                   !     saturated - go down and find p,t, sl at the lcl;
                   !     if sat. exists in first layer (k=ki), use this as lcl
                   !
                ELSE IF (k.EQ.ki) THEN
                   plcl(i)  = press(i,k)
                   tlcl(i)  = tpar(i)
                   tlcl(i)  = tlcl(i) + 1.0_r8
                   klcl(i)  = k
                   llift(i) = .TRUE.
                ELSE
                   dqdp(i)   = (qexces(i)-qex1(i))/ &
                        (press(i,k-1)-press(i,k))
                   deltap(i) = qexces(i)/dqdp(i)
                   plcl(i)   = press(i,k) + deltap(i)
                   tlcl(i)   = tpar(i) * (1.0_r8+2.0_r8*rcp*deltap(i) &
                        /(press(i,k)+press(i,k-1)))
                   tlcl(i)   = tlcl(i) + 1.0_r8
                   klcl(i)   = k
                   llift(i) = .TRUE.
                   !
                   !     give parcel a one-degree goose
                   !
                END IF
                !
                !     lifting cond level found - get out of loop
                !
             END IF
          END IF
       END DO
    END DO
    !
    ! testing the difference in height between lcl and zi, if zi < z(lcl)
    ! there will be no shallow cumulus
    !
    kzi=0
    DO k=3,kMax
       DO i=1,nCols
          IF((thv(i,k)-thv(i,k-1)) > dthv .AND. kzi(i) == 0 )THEN
             kzi(i)=k
          END IF
       END DO
    END DO
    DO i=1,nCols
       lconv(i) = kuo(i) == 1 .or. tfz(i) <= 0.0_r8 .or. .NOT.llift(i) .or. &
            kzi(i) < klcl(i) .or. klcl(i) <= 2
    END DO
    !
    !     quit if parcel still unsat at k = kthick - - -
    !     in this case,set low value of plcl as signal to shalmano
    !
    DO k=1,kMax
       DO i=1,nCols
          IF (.NOT.lconv(i)) THEN
             !
             ! static energy profiles
             !
             se   (i,k)=cp*te(i,k)+grav*ze(i,k)
             uhe  (i,k)=se(i,k)+vlat*qe(i,k)
             uhes (i,k)=se(i,k)+vlat*qes(i,k)
             gamma(i,k)=co1*press(i,k)*qes(i,k)*qes(i,k)/es(i,k)
             gamma(i,k)=gamma(i,k)/(te(i,k)*te(i,k))
          END IF
       END DO
    END DO

    DO k=2,MAXVAL(klcl)
       DO i=1,nCols
          IF (.NOT.lconv(i) .AND. k <= klcl(i) ) THEN
             !
             ! vertical profile of the entrainment rate and cloud moist
             ! static energy
             !
             uhc(i,1)=uhe(i,1)
             uhc(i,k)=uhe(i,1)
          END IF
       END DO
    END DO

    DO k=MINVAL(klcl),kMax
       DO i=1,nCols
          IF (.NOT.lconv(i).AND. k >= klcl(i)+1) THEN
             dlamb  (i,k)=EXP(LOG(dlamb0)+2.3_r8*ze(i,k)/zref)
             dlamb  (i,k)=MIN(0.1_r8,dlamb(i,k))
             dldzby2(i,k)=dlamb(i,k)*delz(i,k)/2.0_r8
             uhc    (i,k)=(uhc(i,k-1)-dldzby2(i,k)*(uhc(i,k-1)-uhe(i,k)- &
                  uhe(i,k-1)))/(1+dldzby2(i,k))
          END IF
       END DO
    END DO

    !
    ! calculating cloud variables qc, ql, sc
    !
    DO k=1,kMax
       DO i=1,nCols
          IF (.NOT.lconv(i)) THEN
             sc (i,k)=se (i,k)+(uhc(i,k)-uhes(i,k))/(1+gamma(i,k))
             sc0(i,k)=se (i,k)+(uhc(i,1)-uhes(i,k))/(1+gamma(i,k))
             qc (i,k)=qes(i,k)+gamma(i,k)*(uhc(i,k)-uhes(i,k))/ &
                  (vlat*(1+gamma(i,k)))
             ql (i,k)=0.0_r8
          END IF
       END DO
    END DO

    DO k=MINVAL(klcl),kMax
       DO i=1,nCols
          IF (.NOT.lconv(i) .AND. k >= klcl(i)+1) THEN
             ql(i,k)= ql(i,k-1)-(qc(i,k)-qc(i,k-1))-dlamb(i,k)* &
                  (qc(i,k)-qe(i,k))*delz(i,k)- &
                  (c0+dlamb(i,k))*ql(i,k-1)*delz(i,k)
             ql(i,k)=MAX(0.00000001_r8,ql(i,k))
          END IF
       END DO
    END DO

    !
    ! determining cloud top based on integrated buoyancy
    !
    DO k=1,kMax
       DO i=1,nCols
          IF (.NOT.lconv(i)) THEN
             scv (i,k)=sc (i,k)+cp*te(i,k)*(0.608_r8*qc(i,k)-ql(i,k))
             sev (i,k)=se (i,k)+0.608_r8*cp*te(i,k)*qe(i,k)
             sc0v(i,k)=sc0(i,k)+0.608_r8*cp*te(i,k)*qe(i,k)
          END IF
       END DO
    END DO

    DO k=2,kMax
       DO i=1,nCols
          IF (.NOT.lconv(i)) THEN
             scvm (i,k)=(scv(i,k)+scv(i,k-1))/2.0_r8
             sevm (i,k)=(sev(i,k)+sev(i,k-1))/2.0_r8
             sc0vm(i,k)=(sc0v(i,k)+sc0v(i,k-1))/2.0_r8
             tvm  (i,k)=(tv(i,k)+tv(i,k-1))/2.0_r8
          END IF
       END DO
    END DO

    DO i=1,nCols
       IF (.NOT.lconv(i)) THEN
          !
          ! determination of the integrated buoyancy between surface and lcl
          ! The calculation assumes that the surface flux is in W/m2. Therefore
          ! we divide it by density*cp in order to convert it to Km/s
          !
          buoy1(i)=0.0_r8
          buoy1(i)=tfz(i)*(1.0_r8+0.608_r8*qe(i,1))/(cp*den(i,1))
          buoy1(i)=grav*ze(i,klcl(i))*buoy1(i)/tv(i,1)
          buoy1(i)=EXP(0.29_r8+0.6667_r8*LOG(buoy1(i)))
          !
          ! checking wether the parcel is able to sustain positive
          ! buoyancy one level above lcl
          !
          cape(i)=0.0_r8
          buoy2(i)=0.0_r8
          buoy2(i)=buoy2(i)+gammad*(scvm(i,klcl(i)+1)-sevm(i,klcl(i)+1))* &
               delz(i,klcl(i)+1)/tvm(i,klcl(i)+1)
       END IF
    END DO

    DO i=1,nCols
       lconv(i) = lconv(i) .or. (buoy1(i)+buoy2(i)) <= 0.0_r8
    END DO



    DO k=MINVAL(klcl),kMax
       DO i=1,nCols
          IF (.NOT.lconv(i) .and. k >= klcl(i)+2) THEN
             !
             ! calculating cloud top and cape
             !
             buoy2(i)=buoy2(i)+gammad*(scvm(i,k)-sevm(i,k))* &
                  delz(i,k)/tvm(i,k)

             IF ((buoy1(i) + buoy2(i)) <= 0.0_r8 .AND. ktop(i) == 1) THEN
                ktop(i)=k-1
             END IF

          END IF
       END DO
    END DO
    DO k=MINVAL(klcl),MAXVAL(ktop)
       DO i=1,nCols
          IF (.NOT.lconv(i) .AND. k >= klcl(i) .AND. k <= ktop(i)) THEN
             emp (i,k)=(sc0vm(i,k)-sevm(i,k))
             emp (i,k)=MAX(0.0_r8,emp(i,k))
             cape(i  )=cape(i)+gammad*emp(i,k)*delz(i,k)/tvm(i,k)
          END IF
       END DO
    END DO
    !
    ! calculating the cloud base mass flux
    !
    ! se quer aumentar fluxo de massa usar tcold(i)=te(i,klcl)
    !
    DO i=1,nCols
       IF (.NOT.lconv(i)) THEN
          thot (i) = te(i,1)
          tcold(i) = te(i,klcl(i))
       END IF
    END DO

    !    icount=1
    !    DO i=1,nCols
    !      IF (.NOT.lconv(i)) THEN
    !        DO k=3,ktop(i)
    !          tcold (i)=tcold (i)+te(i,k)
    !          icount(i)=icount(i)+1
    !        END DO
    !      END IF
    !    END DO

    icount=1
    DO k=3,MAXVAL(ktop)
       DO i=1,nCols
          IF (.NOT.lconv(i) .AND.k <= ktop(i)) THEN
             tcold (i)=tcold (i)+te(i,k)
             icount(i)=icount(i)+1
          END IF
       END DO
    END DO

    DO i=1,nCols
       IF (.NOT.lconv(i)) THEN
          tcold(i)=tcold(i)/icount(i)
          efic(i)=(thot(i)-tcold(i))/thot(i)
       END IF
    END DO

    DO i=1,nCols
       lconv(i)=lconv(i).or.efic(i).LE.0.0_r8.OR.cape(i).LT.40.0_r8
    END DO

    DO i=1,nCols
       IF (.NOT.lconv(i)) THEN
          !
          fin(i)=tfz(i)+qfz(i)
          !
          !snf
          !       tcape parece esta forte......  trocar...
          !       tcape(i)=cape(i)
          !snf    tcape(i)=2*cape(i)
          !       tcape(i)=1.25*cape(i)
          !GOOD   tcape(i)=1.6 *cape(i)
          !
          tcape(i)=par7*cape(i)
          !
          ! noshall parameter define where shallow convection works
          !
          noshal(i)=0
          !
          ! fluxo na base da nuvem
          !
          sigwshb(i)=efic(i)*fin(i)/(den(i,klcl(i))*tcape(i))
       END IF
    END DO

    wssc=0.0_r8
    wqsc=0.0_r8

    DO k= MINVAL(klcl),MAXVAL(ktop)
       DO i=1,nCols
          IF (.NOT.lconv(i) .AND. k >= klcl(i) .AND. k <= ktop(i)) THEN
             wc(i,k)=sigwshb(i)*((ze(i,ktop(i))-ze(i,k))/ &
                  (ze(i,ktop(i))-ze(i,klcl(i))))
             wssc(i,k)=wc(i,k)*(sc(i,k)-vlat*ql(i,k)-se(i,k))
             wqsc(i,k)=wc(i,k)*(qc(i,k)+ql(i,k)-qe(i,k))
          END IF
       END DO
    END DO

    DO k=MINVAL(klcl),MAXVAL(ktop)
       DO i=1,nCols
          IF (.NOT.lconv(i) .AND. k >= klcl(i)+1 .AND. k <= ktop(i)-1) THEN
             dsdt(i,k)=-(wssc(i,k+1)-wssc(i,k-1))/(ze(i,k+1)-ze(i,k-1))
             dtdt(i,k)=  dsdt(i,k)/pi(i,k)
             dqdt(i,k)=-(wqsc(i,k+1)-wqsc(i,k-1))/(ze(i,k+1)-ze(i,k-1))
          END IF
       END DO
    END DO
    !
    ! updating temperature and moisture fields due to shallow convection
    !
    DO k=1,kMax
       DO i=1,nCols
          te(i,k)=te(i,k)+dtdt(i,k)*2.0_r8*dt
          qe(i,k)=qe(i,k)+dqdt(i,k)*2.0_r8*dt
       END DO
    END DO

  END SUBROUTINE shallsouza

  REAL(KIND=r8)  FUNCTION es2(t)
    REAL(KIND=r8) , INTENT(IN) :: t
    REAL(KIND=r8)  :: tx
    REAL(KIND=r8) , PARAMETER  :: d1=0.6107042e0_r8
    REAL(KIND=r8) , PARAMETER  :: d2=4.441157e-2_r8
    REAL(KIND=r8) , PARAMETER  :: d3=1.432098e-3_r8
    REAL(KIND=r8) , PARAMETER  :: d4=2.651396e-5_r8
    REAL(KIND=r8) , PARAMETER  :: d5=3.009998e-7_r8
    REAL(KIND=r8) , PARAMETER  :: d6=2.008880e-9_r8
    REAL(KIND=r8) , PARAMETER  :: d7=6.192623e-12_r8

    ! statement function

    tx = t - tbase
    IF (tx >= -50.0_r8) THEN
       es2 = d1 + tx*(d2 + tx*(d3 + tx*(d4 + tx*(d5 + tx*(d6 + d7*tx)))))
    ELSE
       es2=0.00636e0_r8*EXP(25.6e0_r8*(tx+50.e0_r8)/(tbase-50.e0_r8))
    END IF
  END FUNCTION es2

END MODULE Shall_Souza
