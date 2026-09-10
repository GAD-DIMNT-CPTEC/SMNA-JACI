!
!  $Author: pkubota $
!  $Date: 2006/11/13 12:56:50 $
!  $Revision: 1.2 $
!
MODULE LrgSclPrec


  USE Constants, ONLY :  &
       cp                 , &
       hl                 , &
       gasr               , &
       grav               , &
       rmwmd              , &
       rmwmdi             , &
       e0c                , &
       r8

  USE Options, ONLY :       &
       tbase              , &
       nfprt

    IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: InitLrgscl
  PUBLIC :: lrgscl 
  REAL(KIND=r8) :: rlocp
  REAL(KIND=r8) :: rgrav
  REAL(KIND=r8) :: rlrv
  REAL(KIND=r8) :: const1
  REAL(KIND=r8) :: const2
  REAL(KIND=r8) :: xx1
  REAL(KIND=r8)   , PARAMETER :: tcrit    =   273.15_r8

CONTAINS

  SUBROUTINE InitLrgscl
    REAL(KIND=r8)    :: cpol
    REAL(KIND=r8)    :: const
    rlocp  = hl/cp
    cpol   = cp/hl
    rgrav  = 1.0e0_r8/grav
    rlrv   = -hl/(rmwmdi*gasr)
    const  = 0.1e0_r8*e0c*EXP(-rlrv/tbase)
    const1 = const*rmwmd
    const2 = const*(rmwmd-1.0e0_r8)
    xx1    = rlrv/cpol
  END SUBROUTINE InitLrgscl





  ! lrgscl :calculates the precipitation resulting from large-scale
  !         processes as well as the adjusted temperature and specific
  !         humidity.



  SUBROUTINE lrgscl(geshem, tf, qs, qf, ps, del, sl, dt, &
       mlrg, latco, ncols, kMax)
    !
    !
    !***********************************************************************
    !
    ! lrgscl is called by the subroutine gwater.
    !
    ! lrgscl calls no subroutines.
    !
    !***********************************************************************
    !
    ! argument(dimensions)                       description
    !
    !        geshem(ncols)             input : accumulated precipitation (m)
    !                                         before adding large-scale
    !                                         precipitation for current
    !                                         time step (gaussian).
    !                                output : accumulated precipitation (m)
    !                                         after adding large-scale
    !                                         precipitation for current
    !                                         time step (gaussian).
    !        tf    (imx,kMax)         input : temperature prediction for the
    !                                         current time step before
    !                                         performing large-scale moist
    !                                         processes (gaussian).
    !                                output : temperature prediction for the
    !                                         current time step after
    !                                         performing large-scale moist
    !                                         processes (gaussian).
    !        qs    (imx,kMax)       output : saturation specific humidity
    !                                         (gaussian).
    !        qf    (imx,kMax)        input : specific humidity for the
    !                                         current time step before
    !                                         performing large-scale moist
    !                                         processes (gaussian).
    !                                output : specific humidity for the
    !                                         current time step after
    !                                         performing large-scale moist
    !                                         processes (gaussian).
    !        ps    (imx)              input : predicted surface pressure (cb
    !                                         (gaussian).
    !        del   (kMax)             input : sigma spacing for each layer
    !                                         computed in routine "setsig".
    !        sl    (kMax)             input : sigma values at center of each
    !                                         layer. computed in routine
    !                                         "setsig".
    !        prec  (imx,2)                  : falling precipitation for
    !                                         this time step (gaussian).
    !        super (imx)                    : quantity of precipitable water
    !                                         water in a layer (gaussian).
    !        dpovg (imx,kMax)              : a temporary storage array used
    !                                         in computing saturation
    !                                         specific humidity. also used
    !                                         in computing the precipitable
    !                                         from the specific humidity
    !                                         (gaussian).
    !        evap  (imx)                    : a temporary storage array for
    !                                         evaporation
    !        amtevp(imx,2)                  : the same as above
    !
    !***********************************************************************
    !
    !  ncols......Number of grid points on a gaussian latitude circle  
    !  kMax......Number of sigma levels  
    !  imx.......ncols+1 or ncols+2   :this dimension instead of ncols
    !              is used in order to avoid bank conflict of memory
    !              access in fft computation and make it efficient. the
    !              choice of 1 or 2 depends on the number of banks and
    !              the declared type of grid variable (REAL(KIND=r8)*4,REAL(KIND=r8)*8)
    !              to be fourier transformed.
    !              cyber machine has the symptom.
    !              cray machine has no bank conflict, but the argument
    !              'imx' in subr. fft991 cannot be replaced by ncols     
    !  dt........time interval,usually =delt,but changes
    !            in nlnmi (dt=1.) and at dead start(delt/4,delt/2)
    !  mlrg......mlrg=1 ;output of pre-adjusted & post adjusted temp. &
    !            s.h. in lrgscl
    !  latco.....latitude 
    !
    !***  qf=q(n+1),qs=sat. q at t(n+1)=tf. ps=surf. press(cb)
    !***********************************************************************  
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(inout) :: geshem(ncols)
    REAL(KIND=r8),    INTENT(inout) :: tf    (ncols,kMax) 
    REAL(KIND=r8),    INTENT(inout) :: qf    (ncols,kMax)
    REAL(KIND=r8),    INTENT(inout) :: qs    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: ps    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: del   (kMax)
    REAL(KIND=r8),    INTENT(in   ) :: sl    (kMax)
    REAL(KIND=r8),    INTENT(in   ) :: dt
    INTEGER, INTENT(in   ) :: mlrg
    INTEGER, INTENT(in   ) :: latco
    !
    REAL(KIND=r8)    :: prec  (ncols,2)
    REAL(KIND=r8)    :: dpovg (ncols,kMax)
    REAL(KIND=r8)    :: super (ncols)
    REAL(KIND=r8)    :: fact  (ncols)
    REAL(KIND=r8)    :: evap  (ncols)
    REAL(KIND=r8)    :: amtevp(ncols,2)
    REAL(KIND=r8)    :: q     (ncols,kMax)
    REAL(KIND=r8)    :: t     (ncols,kMax)
    REAL(KIND=r8)    :: dtsqrt
    REAL(KIND=r8)    :: rhsat
    REAL(KIND=r8)    :: qcond
    REAL(KIND=r8)    :: qevap
    !REAL(KIND=r8)    :: esft
    REAL(KIND=r8)    :: aa,cc,expcut,X
    INTEGER :: k
    INTEGER :: i
    INTEGER :: l
    !

    dtsqrt=SQRT(dt)
    X=0.0_r8
    aa=MINEXPONENT(X)
    cc=MAXEXPONENT(X)

    DO k=1,kMax
       DO i=1,ncols
          q (i,k)   = qf(i,k)
          t (i,k)   = tf(i,k)
          !esft=es5(tf(i,k))
          !qs(i,k) = 0.622_r8*esft/(1000.0_r8*(sl(k)*ps(i))-esft)
          !IF(qs(i,k) <= 1.0e-12_r8  )      qs(i,k)=1.0e-12_r8
          expcut=MIN(MAX(rlrv/tf(i,k),aa),cc)
          qs(i,k)   = EXP(expcut)
          IF(sl(k)*ps(i) + const2*qs(i,k) == 0.0_r8)THEN  
              qs(i,k)   = const1*qs(i,k)/MAX(sl(k)*ps(i) + const2*qs(i,k),1.0e-12_r8)
          ELSE
              qs(i,k)   = const1*qs(i,k)/MAX(sl(k)*ps(i) + const2*qs(i,k),1.0e-12_r8)
          END IF
          dpovg(i,k)= (del(k)*rgrav)*ps(i)
       END DO
    END DO

    prec=0.0_r8
    !
    !     pcpn process.....top lyr downward
    !
    rhsat=0.80e0_r8

    DO k=1, kMax
       l= kMax +1-k
       IF(l == 1)rhsat=0.90e0_r8
       DO i=1, ncols
          super(i)=qf(i,l)-qs(i,l)
          !
          !     compute wet-bulb adjustment to t and q, and augment
          !     precipitation falling through column.
          !
          qcond=MAX(0.0e0_r8,super(i))/(1.0e0_r8-xx1*qs(i,l)/(tf(i,l)*tf(i,l)))
          qcond=MAX(0.0e0_r8,qcond)
          tf(i,l)=tf(i,l)+rlocp*qcond
          qf(i,l)=qf(i,l)-qcond
          prec(i,1)=prec(i,1)+qcond*dpovg(i,l)
       END DO
       !
       !     finished with super-saturated point
       !
       evap=0.0_r8
       DO i=1,ncols
          IF(super(i) .LE. 0.0e0_r8.AND.prec(i,1).GT.0.0e0_r8) THEN
             evap(i)=rhsat*qs(i,l)-qf(i,l)
          END IF
          fact(i)=0.32e0_r8*dtsqrt*SQRT(prec(i,1))
       END DO

       amtevp =0.0_r8

       DO i=1,ncols
          IF(evap(i).GT.0.0e0_r8) THEN
             amtevp(i,2)=fact(i)*evap(i)
             prec  (i,2)=prec(i,1)/dpovg(i,l)
          END IF
          amtevp(i,1)=MIN(amtevp(i,2),prec(i,2))
          !
          !     monitor
          !
          qevap     = MAX(0.0e0_r8,amtevp(i,1))/ &
               (1.0e0_r8-xx1*qs(i,l)/(tf(i,l)*tf(i,l)))
          tf(i,l)   = tf(i,l)-rlocp*qevap
          qf(i,l)   = qf(i,l)+qevap
          prec(i,1) = prec(i,1)-qevap*dpovg(i,l)
          prec(i,1) = MAX(0.0e0_r8,prec(i,1))
       END DO

    END DO
    !
    !     pcpn reaches ground level....factor of .5 since pcpn is for
    !     two (leapfrog) time-steps....accum pcpn (geshem) is not leapfrogg
    IF(mlrg.EQ.1) THEN
       !
       !     monitor
       !
       DO i=1,ncols
          IF(prec(i,1).GT.0.0e0_r8) THEN
             WRITE(nfprt,999)i,latco,(t (i,k),k=1,kMax)
             WRITE(nfprt,888)        (tf(i,k),k=1,kMax)
             WRITE(nfprt,777)        (q (i,k),k=1,kMax)
             WRITE(nfprt,666)        (qf(i,k),k=1,kMax)
          END IF
       END DO

    END IF

    DO i=1,ncols
       geshem(i)=geshem(i)+prec(i,1)*0.5e0_r8
    END DO
666 FORMAT(' QLN '         ,8X,10E12.5)
777 FORMAT(' QLO '         ,8X,10E12.5)
888 FORMAT(' TLN '         ,8X,10E12.5)
999 FORMAT(' TLO ',I3,1X,I3,1X,10E12.5)
  END SUBROUTINE lrgscl

  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    REAL(KIND=r8), INTENT(IN) :: t
    REAL(KIND=r8)            :: ae  (2)
    REAL(KIND=r8)            :: be  (2)
    REAL(KIND=r8)            :: ht  (2)

    !
    ht(1)=hl/cp
    
    ht(2)=2.834e6_r8/cp
    
    be(1)=0.622_r8*ht(1)/0.286_r8
    
    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)
    
    be(2)=0.622_r8*ht(2)/0.286_r8
    
    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)


    IF (t <= tcrit) THEN
       es5 = EXP(ae(2)-be(2)/t)
    ELSE
       es5 = EXP(ae(1)-be(1)/t)
    END IF
  END FUNCTION es5

END MODULE LrgSclPrec
