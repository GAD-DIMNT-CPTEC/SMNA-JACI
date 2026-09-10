!
!  $Author: pkubota $
!  $Date: 2010/04/20 20:18:04 $
!  $Revision: 1.5 $
!
MODULE ModConRas

  !                 |Init_Cu_RelAraSch
  !                 |
  !   gwater2-------|
  !                 |
  !                 |RunCu_RelAraSch-------|ras--------|qsat
  !                 |             |           |
  !                 |             |           |cloud-------|acritn
  !                 |             |                        |
  !                 |             |                        |rncl
  !                 |             |rnevp------|qstar9
  !                 |
  !                 |shllcl-------|mstad2
  !                 |



 USE Constants, ONLY :  &
       cp                 , &
       hl                 , &
       gasr               , &
       grav               , &
       rmwmd              , &
       r8   
USE Options, ONLY :       &
       sthick            , &
       tbase             , &
       ki

USE Parallelism, ONLY: &
     MsgOne, &
     FatalError

  IMPLICIT NONE
SAVE

  PRIVATE

  PUBLIC :: Init_Cu_RelAraSch
  PUBLIC :: RunCu_RelAraSch
  PUBLIC :: shllcl
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

  REAL(KIND=r8), PARAMETER :: xkapa=0.2857143_r8
  REAL(KIND=r8), PARAMETER  :: d1=.6107042e0_r8
  REAL(KIND=r8), PARAMETER  :: d2=4.441157e-2_r8
  REAL(KIND=r8), PARAMETER  :: d3=1.432098e-3_r8
  REAL(KIND=r8), PARAMETER  :: d4=2.651396e-5_r8
  REAL(KIND=r8), PARAMETER  :: d5=3.009998e-7_r8
  REAL(KIND=r8), PARAMETER  :: d6=2.008880e-9_r8
  REAL(KIND=r8), PARAMETER  :: d7=6.192623e-12_r8

CONTAINS

  SUBROUTINE Init_Cu_RelAraSch()

    CALL InitAcritn()
    CALL InitMstad2()
    
  END SUBROUTINE Init_Cu_RelAraSch



  REAL(KIND=r8) FUNCTION es3(t)
    REAL(KIND=r8), INTENT(IN) :: t
    REAL(KIND=r8) :: tx

    ! statement function

    tx = t - tbase
    IF (tx >= -50.0_r8) THEN
       es3 = d1 + tx*(d2 + tx*(d3 + tx*(d4 + tx*(d5 + tx*(d6 + d7*tx)))))
    ELSE
       es3=0.00636e0_r8*EXP(25.6e0_r8*(tx+50.e0_r8)/(tbase-50.e0_r8))
    END IF
  END FUNCTION es3





  SUBROUTINE InitAcritn()
    REAL(KIND=r8),    PARAMETER :: actp=1.7_r8
    REAL(KIND=r8),    PARAMETER :: facm=1.0_r8
    REAL(KIND=r8)    :: au(15) 
    REAL(KIND=r8)    :: ph(15) 
    REAL(KIND=r8)    :: tem
    INTEGER :: l

    ph =(/150.0_r8, 200.0_r8, 250.0_r8, 300.0_r8, 350.0_r8, 400.0_r8, 450.0_r8, 500.0_r8, &
          550.0_r8, 600.0_r8, 650.0_r8, 700.0_r8, 750.0_r8, 800.0_r8, 850.0_r8/)

    aa =(/ 1.6851_r8, 1.1686_r8, 0.7663_r8, 0.5255_r8, 0.4100_r8, 0.3677_r8, &
           0.3151_r8, 0.2216_r8, 0.1521_r8, 0.1082_r8, 0.0750_r8, 0.0664_r8, &
           0.0553_r8, 0.0445_r8, 0.0633_r8/)

    ad = 0.0_r8
    ac = 0.0_r8

    actop   = actp*facm

    DO l=1,15
       aa(l) = aa(l)*facm
    END DO

    DO l=2,15
       tem   = ph(l) - ph(l-1)
       au(l) = aa(l-1) / tem
       ad(l) = aa(l)   / tem
       ac(l) = ph(l)*au(l) - ph(l-1)*ad(l)
       ad(l) = ad(l) - au(l)
    END DO
  END SUBROUTINE InitAcritn



  ! acritn :relaxed arakawa-schubert.



  SUBROUTINE acritn(len, pl, plb, acr)
    INTEGER, INTENT(IN ) :: len
    REAL(KIND=r8),    INTENT(IN ) :: pl(len) 
    REAL(KIND=r8),    INTENT(IN ) :: plb(len) 
    REAL(KIND=r8),    INTENT(INOUT) :: acr(len)

    INTEGER :: i
    INTEGER :: iwk

    DO i = 1, len
       iwk = INT((pl(i) * 0.02_r8 - 0.999999999_r8))
       IF (iwk > 1) THEN
          IF (iwk <= 15) THEN
             acr(i) = ac(iwk) + pl(i) * ad(iwk)
          ELSE
             acr(i) = aa(15)
          END IF
       ELSE   
          acr(i) = actop
       END IF
       acr(i) = acr(i) * (plb(i) - pl(i))
    END DO
  END SUBROUTINE acritn



  ! qsat   :relaxed arakawa-schubert.



  SUBROUTINE qsat(tt,p,q,dqdt,ldqdt)
    REAL(KIND=r8),    INTENT(IN )    :: tt
    REAL(KIND=r8),    INTENT(IN )    :: p
    REAL(KIND=r8),    INTENT(INOUT)    :: q
    REAL(KIND=r8),    INTENT(INOUT)    :: dqdt  
    LOGICAL, INTENT(IN )    :: ldqdt

    REAL(KIND=r8), PARAMETER :: airmw = 28.97_r8   
    REAL(KIND=r8), PARAMETER :: h2omw = 18.01_r8   
    REAL(KIND=r8), PARAMETER :: one   = 1.0_r8   
    REAL(KIND=r8), PARAMETER :: esfac = h2omw/airmw        
    REAL(KIND=r8), PARAMETER :: erfac = (one-esfac)/esfac 
    REAL(KIND=r8), PARAMETER :: b6  = 6.136820929e-11_r8*esfac
    REAL(KIND=r8), PARAMETER :: b5  = 2.034080948e-8_r8 *esfac
    REAL(KIND=r8), PARAMETER :: b4  = 3.031240396e-6_r8 *esfac
    REAL(KIND=r8), PARAMETER :: b3  = 2.650648471e-4_r8 *esfac
    REAL(KIND=r8), PARAMETER :: b2  = 1.428945805e-2_r8 *esfac
    REAL(KIND=r8), PARAMETER :: b1  = 4.436518521e-1_r8 *esfac
    REAL(KIND=r8), PARAMETER :: b0  = 6.107799961e+0_r8 *esfac
    REAL(KIND=r8), PARAMETER :: c1  = b1   
    REAL(KIND=r8), PARAMETER :: c2  = b2*2.0_r8
    REAL(KIND=r8), PARAMETER :: c3  = b3*3.0_r8
    REAL(KIND=r8), PARAMETER :: c4  = b4*4.0_r8
    REAL(KIND=r8), PARAMETER :: c5  = b5*5.0_r8
    REAL(KIND=r8), PARAMETER :: c6  = b6*6.0_r8
    REAL(KIND=r8), PARAMETER :: tmin=223.15_r8
    REAL(KIND=r8), PARAMETER :: tmax=323.15_r8 
    REAL(KIND=r8), PARAMETER :: tice=273.16_r8
    REAL(KIND=r8) :: t
    REAL(KIND=r8) :: d

    t =  MIN ( MAX (tt,tmin),tmax) - tice
    q = (t*(t*(t*(t*(t*(t*b6+b5)+b4)+b3)+b2)+b1)+b0)
    d = one / (p-erfac*q)
    q = q * d

    IF (ldqdt)  THEN
       dqdt = (t*(t*(t*(t*(t*c6+c5)+c4)+c3)+c2)+c1)
       dqdt = (dqdt + erfac*q) * d
    END IF

  END SUBROUTINE qsat



  ! rncl   :relaxed arakawa-schubert



  SUBROUTINE rncl(len, pl, rno, clf)
    INTEGER, INTENT(IN ) :: len
    REAL(KIND=r8),    INTENT(IN ) :: pl (len)  
    REAL(KIND=r8),    INTENT(INOUT) :: rno(len)
    REAL(KIND=r8),    INTENT(INOUT) :: clf(len)

    REAL(KIND=r8), PARAMETER :: p5=500.0_r8
    REAL(KIND=r8), PARAMETER :: p8=800.0_r8
    REAL(KIND=r8), PARAMETER :: pt8=0.8_r8
    REAL(KIND=r8), PARAMETER :: pfac=0.2_r8/(p8-p5)
    REAL(KIND=r8), PARAMETER :: cucld=0.5_r8
    INTEGER :: i

    DO i = 1, len
       rno(i) = 1.0_r8
       clf(i) = cucld
       IF (pl(i) >= p5 .AND. pl(i) <= p8) THEN
          rno(i) = (p8-pl(i))*pfac + pt8
       ELSE IF (pl(i) > p8) THEN
          rno(i) = pt8
       END IF
    END DO
  END SUBROUTINE rncl





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
    REAL(KIND=r8)          :: crit   
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
       esat = es3(t)

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

    ! construct table to get temperature from thetae and p
    ! the table is tfmthe(i,k), where
    ! thetae=220,230,--,500 = 210+10*i, and
    ! p=press/100 cb = 0,.1_r8,.2_r8,--,1.1_r8 =  0.1_r8*(k-1)
    ! for p=0, tfmthe(i,1) is identically 0.0e0_r8

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
             ess = es3(t)
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
          ess = es3(t)
          qfmthe(i,k)=eps*ess/(1.e2_r8*p-(1.0e0_r8-eps)*ess)
          p=p+fm_dpres
       END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE InitMstad2




  ! mstad2 :this subroutine lifts the parcel from the lcl and computes the
  !         parcel temperature at each level; temperature is retrieved from
  !         lookup tables defined in the first call.
  !         This version-mstad2-lifts parcel from lifting condensation
  !         level, where pres is slcl and temp is tlcl
  !         level ll is the first regular level above the lcl
  !         this version contains double resolution tables
  !         and uses virtual temp in the buoyancy test
  !         new sat vapor pressure used





  SUBROUTINE mstad2(ps, sig, tin, tmst, qmst, ktop,  ier, & 
       slcl, ll, qin, tlcl, llift, ncols, kmax)
       
    !
    ! this routine accepts input data of
    !         ps                  = surface pressure dvded by 100 cbs
    !         sig(k=1,--,km)   = factor such that p at lvl k = ps*sig(k)
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
    REAL(KIND=r8),    INTENT(in   ) :: sig(kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tin  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: tmst (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: qmst (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: qin  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: ps    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: slcl (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tlcl (ncols)
    INTEGER, INTENT(inout) :: ktop (ncols)
    INTEGER, INTENT(inout) :: ier  (ncols)
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
    !     surface press is ps, pressure at other lvls is ps times sig(k)
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
             ier(i) = 1
             WRITE(line,'(g12.4," < ",g12.4)') ti(i),the_temp0
             CALL FatalError(h//"Theta_e table: temp="//line//"=low limit")
          ENDIF
          IF(jt(i).GE.the_ntemp) THEN
             ier(i) = 1
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
             ier(i) = 1
             CALL FatalError(h//"Theta_e table: pres="//line//"=low limit")
          ENDIF
          IF(kp(i).GE.the_npres) THEN
             ier(i) = 1
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
          IF(kt(i).LT.1) THEN
             ier(i) = 1
             CALL FatalError(h//&
               "tfmthe_e table: the < low limit")
          END IF 
          IF(kt(i).GE.fm_nthe) THEN
             ier(i) = 1
             CALL FatalError(h//&
               "tfmthe_e table: the > high limit")
          END IF
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
          IF(ip(i).LT.1)THEN
             ier(i) = 1
            CALL FatalError(h//&
               "tfmthe_e table: pres < low limit")
          END IF
          IF(ip(i).GE.fm_npres)THEN
             ier(i) = 1
            CALL FatalError(h//&
               "tfmthe_e table: pres > high limit")
          END IF
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
             pp(i)=ps(i)*sig(k)
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







  SUBROUTINE qstar9 (t     ,p     ,n     ,qstt  ,layr  ,kmax  )
    !
    !
    ! function vqsat
    !     purpose
    !     vector computation of saturation mixing ratio
    !     usage
    !     qsatt(1;n) = vqsat(t(1;n),p(1;n),n;qsatt(1;n))
    !     description of parameters
    !     t        temperature vector (deg k)
    !     p        pressure vector (mb)
    !     vqsat    saturation mixing ratio vector
    !     kmax.....Number of sigma levels  
    !     n........Number of grid points on a gaussian latitude circle
    !     qstt.....saturation mixing ratio
    !     layr.....=2,3,4,5,...,kmax
    INTEGER, INTENT(in   ) :: kmax  
    INTEGER, INTENT(in   ) :: n
    REAL(KIND=r8)   , INTENT(in   ) :: t(n)
    REAL(KIND=r8)   , INTENT(in   ) :: p(n)
    REAL(KIND=r8)   , INTENT(inout) :: qstt(n)
    INTEGER, INTENT(in   ) :: layr

    REAL(KIND=r8)                   :: e1(n)
    REAL(KIND=r8)                   :: e2(n)
    REAL(KIND=r8)                   :: tq(n)
    INTEGER                :: i1(n)
    INTEGER                :: i2(n)
    REAL(KIND=r8)   , PARAMETER     :: zp622  = 0.622_r8  
    REAL(KIND=r8)   , PARAMETER     :: z1p0s1 = 1.00001_r8 
    REAL(KIND=r8)   , PARAMETER     :: z1p622 = 1.622_r8  
    REAL(KIND=r8)   , PARAMETER     :: z138p9 = 138.90001_r8 
    REAL(KIND=r8)   , PARAMETER     :: z198p9 = 198.99999_r8 
    REAL(KIND=r8)   , PARAMETER     :: z200   = 200.0_r8    
    REAL(KIND=r8)   , PARAMETER     :: z337p9 = 337.9_r8    
    REAL(KIND=r8)   , PARAMETER     :: est(139) = (/ &
         0.31195e-02_r8, 0.36135e-02_r8, 0.41800e-02_r8, &
         0.48227e-02_r8, 0.55571e-02_r8, 0.63934e-02_r8, 0.73433e-02_r8, &
         0.84286e-02_r8, 0.96407e-02_r8, 0.11014e-01_r8, 0.12582e-01_r8, &
         0.14353e-01_r8, 0.16341e-01_r8, 0.18574e-01_r8, 0.21095e-01_r8, &
         0.23926e-01_r8, 0.27096e-01_r8, 0.30652e-01_r8, 0.34629e-01_r8, &
         0.39073e-01_r8, 0.44028e-01_r8, 0.49546e-01_r8, 0.55691e-01_r8, &
         0.62508e-01_r8, 0.70077e-01_r8, 0.78700e-01_r8, 0.88128e-01_r8, &
         0.98477e-01_r8, 0.10983e+00_r8, 0.12233e+00_r8, 0.13608e+00_r8, &
         0.15121e+00_r8, 0.16784e+00_r8, 0.18615e+00_r8, 0.20627e+00_r8, &
         0.22837e+00_r8, 0.25263e+00_r8, 0.27923e+00_r8, 0.30838e+00_r8, &
         0.34030e+00_r8, 0.37520e+00_r8, 0.41334e+00_r8, 0.45497e+00_r8, &
         0.50037e+00_r8, 0.54984e+00_r8, 0.60369e+00_r8, 0.66225e+00_r8, &
         0.72589e+00_r8, 0.79497e+00_r8, 0.86991e+00_r8, 0.95113e+00_r8, &
         0.10391e+01_r8, 0.11343e+01_r8, 0.12372e+01_r8, 0.13484e+01_r8, &
         0.14684e+01_r8, 0.15979e+01_r8, 0.17375e+01_r8, 0.18879e+01_r8, &
         0.20499e+01_r8, 0.22241e+01_r8, 0.24113e+01_r8, 0.26126e+01_r8, &
         0.28286e+01_r8, 0.30604e+01_r8, 0.33091e+01_r8, 0.35755e+01_r8, &
         0.38608e+01_r8, 0.41663e+01_r8, 0.44930e+01_r8, 0.48423e+01_r8, &
         0.52155e+01_r8, 0.56140e+01_r8, 0.60394e+01_r8, 0.64930e+01_r8, &
         0.69767e+01_r8, 0.74919e+01_r8, 0.80406e+01_r8, 0.86246e+01_r8, &
         0.92457e+01_r8, 0.99061e+01_r8, 0.10608e+02_r8, 0.11353e+02_r8, &
         0.12144e+02_r8, 0.12983e+02_r8, 0.13873e+02_r8, 0.14816e+02_r8, &
         0.15815e+02_r8, 0.16872e+02_r8, 0.17992e+02_r8, 0.19176e+02_r8, &
         0.20428e+02_r8, 0.21750e+02_r8, 0.23148e+02_r8, 0.24623e+02_r8, &
         0.26180e+02_r8, 0.27822e+02_r8, 0.29553e+02_r8, 0.31378e+02_r8, &
         0.33300e+02_r8, 0.35324e+02_r8, 0.37454e+02_r8, 0.39696e+02_r8, &
         0.42053e+02_r8, 0.44531e+02_r8, 0.47134e+02_r8, 0.49869e+02_r8, &
         0.52741e+02_r8, 0.55754e+02_r8, 0.58916e+02_r8, 0.62232e+02_r8, &
         0.65708e+02_r8, 0.69351e+02_r8, 0.73168e+02_r8, 0.77164e+02_r8, &
         0.81348e+02_r8, 0.85725e+02_r8, 0.90305e+02_r8, 0.95094e+02_r8, &
         0.10010e+03_r8, 0.10533e+03_r8, 0.11080e+03_r8, 0.11650e+03_r8, &
         0.12246e+03_r8, 0.12868e+03_r8, 0.13517e+03_r8, 0.14193e+03_r8, &
         0.14899e+03_r8, 0.15634e+03_r8, 0.16400e+03_r8, 0.17199e+03_r8, &
         0.18030e+03_r8, 0.18895e+03_r8, 0.19796e+03_r8, 0.20733e+03_r8, &
         0.21708e+03_r8, 0.22722e+03_r8, 0.23776e+03_r8, 0.24871e+03_r8/)
    REAL(KIND=r8)                   :: qfac(kmax)
    INTEGER                :: i
    REAL(KIND=r8)                   :: a1622

    a1622   = 1.0e0_r8  / z1p622
    qfac    = 1.0e0_r8

    DO i = 1,n
       tq(i) = t(i) - z198p9
       IF(t(i).LT.z200   ) tq(i) = z1p0s1
       IF(t(i).GT.z337p9 ) tq(i) = z138p9
       i1(i) = INT(tq(i))
       i2(i) = i1(i) + 1
    END DO

    DO i = 1,n
       e1(i) =  est(  i1(i) )
       e2(i) =  est(  i2(i) )
    END DO

    DO i = 1,n
       qstt(i)   = tq(i) - i1(i)
       qstt(i)   = qstt(i) * ( e2(i)-e1(i) )
       qstt(i)   = qstt(i) +   e1(i)
       e1(i) = p(i) * a1622
       IF( e1(i).LT.qstt(i) ) qstt(i) = e1(i)
       qstt(i)   = zp622 *  qstt(i) / p(i)
       qstt(i) = qstt(i)*qfac(layr)
    END DO

  END SUBROUTINE qstar9



  ! rnevp  :relaxed arakawa-schubert



  SUBROUTINE rnevp(ft    ,fq    ,pain  ,fp    ,sig   ,dsig  ,clfric, &
       rcon  ,rlar  ,dtc3  ,ncols ,kmax  ,nlst)
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kmax
    !
    !     default parameter statements for dimensioning resolution
    !     
    INTEGER, INTENT(in   ) :: nlst
    REAL(KIND=r8),    INTENT(in   ) :: dtc3
    !     
    !     input model fields
    !
    REAL(KIND=r8),    INTENT(in   ) :: fp  (ncols) 
    REAL(KIND=r8),    INTENT(inout) :: fq  (ncols,kmax) 
    REAL(KIND=r8),    INTENT(inout) :: ft  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: dsig(kmax) 
    REAL(KIND=r8),    INTENT(in   ) :: sig (kmax)
    !     
    !     integers and logicals are listed first, then half precisions
    !
    REAL(KIND=r8),    INTENT(in   ) :: clfric
    REAL(KIND=r8),    INTENT(in   ) :: pain (ncols,kmax ) 
    !
    !     this is need to use the dao type re-evaporation
    ! 
    REAL(KIND=r8),    INTENT(inout) :: rcon(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rlar(ncols)   
    !     
    !     integers and logicals are listed first, then half precisions
    !
    LOGICAL :: detprc
    LOGICAL :: timeit
    REAL(KIND=r8)    :: shsat    (ncols,kmax )
    !
    !     this is need to use the dao type re-evaporation
    !
    REAL(KIND=r8)    :: temp2 (ncols,kmax)
    REAL(KIND=r8)    :: tl    (ncols,kmax)
    REAL(KIND=r8)    :: ql    (ncols,kmax)
    REAL(KIND=r8)    :: rain  (ncols,kmax)
    REAL(KIND=r8)    :: pl    (ncols,kmax)
    REAL(KIND=r8)    :: amass (ncols,kmax)
    REAL(KIND=r8)    :: clfrac(ncols,kmax)
    REAL(KIND=r8)    :: tempc (ncols,kmax)
    REAL(KIND=r8)    :: cvt   (ncols,kmax)
    REAL(KIND=r8)    :: cvq   (ncols,kmax)
    REAL(KIND=r8)    :: evp9  (ncols,kmax)
    REAL(KIND=r8)    :: plk   (ncols,kmax)
    REAL(KIND=r8)    :: art   (ncols)
    INTEGER :: im 
    INTEGER :: nlay  
    INTEGER :: nlayp1
    INTEGER :: nlaym1
    INTEGER :: nlaym2
    INTEGER :: ls 
    INTEGER :: ls1 
    REAL(KIND=r8)    :: grav  
    REAL(KIND=r8)    :: grav2 
    REAL(KIND=r8)    :: pi 
    REAL(KIND=r8)    :: pi2 
    REAL(KIND=r8)    :: gamfac
    REAL(KIND=r8)    :: cp 
    REAL(KIND=r8)    :: velta 
    REAL(KIND=r8)    :: ptop  
    REAL(KIND=r8)    :: hltm  
    REAL(KIND=r8)    :: pcon  
    REAL(KIND=r8)    :: critl 
    REAL(KIND=r8)    :: flim  
    REAL(KIND=r8)    :: dhtop 
    INTEGER :: l 
    INTEGER :: i 
    REAL(KIND=r8)    :: rphf  
    REAL(KIND=r8)    :: elocp 
    INTEGER :: n 
    REAL(KIND=r8)    :: relax 
    REAL(KIND=r8)    :: rnfrac
    REAL(KIND=r8)    :: rpow  
    REAL(KIND=r8)    :: exparg

    !
    !     to save storage space use the temp array from comp3 for local store
    !                                                                     
    !     temp level 1 has accumulated evaporation of convective rain
    !     temp level 2 has saturation mixing ratio                   
    !     temp level 3 has local temperature                         
    !     temp level 4 has local mixing ratio                        
    !     temp level 5 has change in mixing ratio (at each iteration)
    !     temp level 6 has dqdt, one of the terms for the iteration  
    !     temp level 7 is  a temporary for actual evaporation of rain
    !     
    !     
    !     flip arrays and transfer to local variables
    !     
    im=ncols
    nlay=kmax
    nlayp1=nlay+1
    nlaym1=nlay-1
    nlaym2=nlay-2
    ls=nlst
    ls1=nlst+1
    detprc=.TRUE.
    timeit=.TRUE.
    grav=9.81_r8
    grav2=0.01_r8*grav
    pi=4.0_r8*ATAN(1.0_r8)
    pi2=2.0_r8*pi
    gamfac=1.348e7_r8
    cp=1003.0_r8
    velta=0.608_r8
    ptop=0.01_r8
    hltm=2.56e6_r8
    pcon=2.e-3_r8
    critl=0.0_r8
    flim=765.0_r8
    dhtop=0.95_r8
    DO l = 1,kmax
       DO i = 1,ncols
          tl(i,l) = ft(i,kmax -l+1)
          ql(i,l) = fq(i,kmax-l+1)
          pl(i,l) = sig(kmax-l+1)*fp(i)*10.0_r8
          plk(i,l) = EXP(xkapa* LOG(pl(i,l)))
          amass(i,l) = 0.01_r8*grav/(10.0_r8*fp(i)*dsig(kmax-l+1))
          temp2(i,l)= pl(i,l) *0.001_r8
          temp2(i,l)=SQRT(temp2(i,l))
       END DO
    END DO
    DO l=1,kmax
       DO i=1,ncols
          rain(i,l)=pain(i,l)
       END DO
    END DO
    rphf=3600.0_r8/dtc3
    elocp=hltm/cp
    DO i = 1,im
       rcon(i) = 0.0e0_r8
       rlar(i) = 0.0e0_r8
    END DO
    DO l=1,nlay
       DO i = 1,im
          evp9 (i,l)  = 0.0e0_r8
          tempc(i,l) = 0.0e0_r8
          cvt  (i,l)   = 0.0e0_r8
          cvq  (i,l)   = 0.0e0_r8
       END DO
    END DO
    ! 
    !     do loop for moisture evaporation ability and convec evaporation
    !
    DO l = 2,nlay
       DO i = 1,im
          tempc(i,3) = tl(i,l)
          tempc(i,4) = ql(i,l)
       END DO
       DO n = 1,3
          CALL qstar9(tempc(1,3),pl(1,l),im   ,tempc(1,2),l    ,kmax  )
          IF(n .EQ. 1) THEN
             relax = 0.5_r8
          ELSE
             relax = 1.0_r8
          END IF
          DO i = 1,im
             tempc(i,5) = tempc(i,2) - tempc(i,4)
             tempc(i,6) = tempc(i,2)*gamfac/tempc(i,3)**2.0_r8
             tempc(i,5) = tempc(i,5)/(1.0e0_r8+tempc(i,6))
             tempc(i,4) = tempc(i,4) + tempc(i,5)*relax
             tempc(i,3) = tempc(i,3) - tempc(i,5)*elocp*relax
          END DO
       END DO
       DO i = 1,im
          shsat(i,l) = tempc(i,4)
          evp9(i,l) = (tempc(i,4) - ql(i,l))/amass(i,l)
          rcon(i) = rcon(i) + rain(i,l)
          clfrac(i,l) = (1.0_r8-(ql(i,l)/tempc(i,4)))*clfric
          tempc(i,8) = MAX(tempc(i,8),clfrac(i,l))
          art(i) = 0.0e0_r8
          IF(rcon(i).GT.0.0e0_r8 .AND. evp9(i,l).GT.0.0e0_r8) THEN
             rnfrac = tempc(i,8)
             rnfrac = MIN(rnfrac,1.0e0_r8)
             rpow = EXP(0.578_r8* LOG(rcon(i)*rphf*temp2(i,l)))
             exparg = -1.04e-4_r8*dtc3*rpow
             art(i) = 1.0_r8 - (EXP(exparg))
             tempc(i,7) = evp9(i,l)*art(i)*rnfrac
             IF(tempc(i,7) .GE. rcon(i)) tempc(i,7) = rcon(i)
             rcon(i)   = rcon(i) - tempc(i,7)
             evp9(i,l) = evp9(i,l) - tempc(i,7)
             cvq(i,l) = cvq(i,l) + tempc(i,7)*amass(i,l)
             cvt(i,l) = cvt(i,l) - tempc(i,7)*amass(i,l)*elocp
             !
             !     update as rain falls and reevaporates
             !
             tl(i,l) = tl(i,l) + cvt(i,l)
             ql(i,l) = ql(i,l) + cvq(i,l)
          ELSE
             tempc(i,7) = 0.0e0_r8
          ENDIF
       END DO
    END DO
    !
    !     flip arrays for return (only flip t and q)
    !
    DO l = 1,nlay
       DO i = 1,im
          ft(i,l) = tl(i,kmax -l+1)
          fq(i,l) = ql(i,kmax-l+1)
       END DO
    END DO
  END SUBROUTINE rnevp



  ! cloud  :relaxed arakawa-schubert.



  SUBROUTINE cloud(  &
       len   ,lenc  ,k     ,ic    ,rasalf,setras,frac  ,alhl  ,rkap  , &
       poi   ,qoi   ,uoi   ,voi   ,prs   ,prj   ,pcu   ,cln   ,tcu   , &
       qcu   ,alf   ,bet   ,gam   ,prh   ,pri   ,hol   ,eta   ,hst   , &
       qol   ,gmh   ,tx1   ,tx2   ,tx3   ,tx4   ,tx5   ,tx6   ,tx7   , &
       tx8   ,alm   ,wfn   ,akm   ,qs1   ,clf   ,uht   ,vht   ,wlq   , &
       ia    ,i1    ,i2    ,cmb2pa,rhmax )
    ! cmb2pa....Parameter cmb2pa = 100.0_r8
    ! rhmax.....Parameter rhmax  = 0.9999_r8    
    ! ncols......Number of grid points on a gaussian latitude circle    
    ! actp......Parameter actp   = 1.7_r8   
    ! facm......Parameter facm   = 1.0_r8      
    ! p5........Parameter p5     = 500.0_r8  
    ! p8........Parameter p8     = 800.0_r8  
    ! pt8.......Parameter pt8    = 0.8_r8 
    ! pt2.......Parameter pt2    = 0.2_r8 
    ! pfac......Parameter pfac   = pt2/(p8-p5)    
    ! cucld.....Parameter cucld  = 0.5_r8   
    ! len.......ncols
    ! lenc......ncols
    ! k.........kmax
    ! ic
    ! rasalf
    ! setras
    ! frac......frac=1.0_r8 
    ! alhl......alhl=2.52e6_r8
    ! rkap......rkap=r/cp=287.05e0_r8/1004.6e0_r8
    ! poi.......array pot. temp. 
    ! qoi.......array humidity.
    ! prs
    ! uoi
    ! voi
    ! prj
    ! tcu
    ! qcu
    ! cln
    ! alf
    ! bet
    ! gam
    ! prh
    ! pri
    ! akm
    ! wfn
    ! hol
    ! qol
    ! eta
    ! hst
    ! gmh
    ! alm
    ! wlq
    ! qs1
    ! tx1
    ! tx2
    ! tx3
    ! tx4
    ! tx5
    ! tx6
    ! tx7
    ! tx8
    ! uht
    ! vht
    ! clf
    ! pcu
    ! ia
    ! i1
    ! i2
    !
    !==========================================================================
    REAL(KIND=r8),    INTENT(in   ) :: cmb2pa
    REAL(KIND=r8),    INTENT(in   ) :: rhmax  

    INTEGER, INTENT(in   ) :: len
    INTEGER, INTENT(in   ) :: lenc
    INTEGER, INTENT(in   ) :: k
    INTEGER, INTENT(in   ) :: ic
    REAL(KIND=r8),    INTENT(in   ) :: rasalf
    LOGICAL, INTENT(in   ) :: setras
    REAL(KIND=r8),    INTENT(in   ) :: frac 
    REAL(KIND=r8),    INTENT(in   ) :: alhl
    REAL(KIND=r8),    INTENT(in   ) :: rkap

    REAL(KIND=r8),    INTENT(in   ) :: poi(len,k)
    REAL(KIND=r8),    INTENT(in   ) :: qoi(len,k) 
    REAL(KIND=r8),    INTENT(in   ) :: prs(len,k+1)
    REAL(KIND=r8),    INTENT(in   ) :: uoi(len,k)
    REAL(KIND=r8),    INTENT(in   ) :: voi(len,k)
    REAL(KIND=r8),    INTENT(in   ) :: prj(len,k+1)
    REAL(KIND=r8),    INTENT(inout  ) :: tcu(len,k) 
    REAL(KIND=r8),    INTENT(inout  ) :: qcu(len,k)
    REAL(KIND=r8),    INTENT(inout  ) :: cln(len)

    REAL(KIND=r8),    INTENT(inout) :: alf(len,k)
    REAL(KIND=r8),    INTENT(inout) :: bet(len,k) 
    REAL(KIND=r8),    INTENT(inout) :: gam(len,k)
    REAL(KIND=r8),    INTENT(inout) :: prh(len,k)
    REAL(KIND=r8),    INTENT(inout) :: pri(len,k)
    REAL(KIND=r8),    INTENT(inout) :: akm(lenc)  
    REAL(KIND=r8),    INTENT(inout) :: wfn(lenc)

    REAL(KIND=r8),    INTENT(inout) :: hol(lenc,k) 
    REAL(KIND=r8),    INTENT(inout) :: qol(lenc,k)  
    REAL(KIND=r8),    INTENT(inout) :: eta(lenc,k) 
    REAL(KIND=r8),    INTENT(inout) :: hst(lenc,k)
    REAL(KIND=r8),    INTENT(inout) :: gmh(lenc,k) 
    REAL(KIND=r8),    INTENT(inout) :: alm(lenc) 
    REAL(KIND=r8),    INTENT(inout) :: wlq(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: qs1(lenc)
    REAL(KIND=r8),    INTENT(inout) :: tx1(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: tx2(lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx3(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: tx4(lenc)
    REAL(KIND=r8),    INTENT(inout) :: tx5(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: tx6(lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx7(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: tx8(lenc)
    REAL(KIND=r8),    INTENT(inout) :: uht(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: vht(lenc) 
    REAL(KIND=r8),    INTENT(inout  ) :: clf(lenc)   
    REAL(KIND=r8),    INTENT(inout) :: pcu(lenc)

    INTEGER, INTENT(inout) :: ia(lenc)
    INTEGER, INTENT(inout) :: i1(lenc)
    INTEGER, INTENT(inout) :: i2(lenc)

    REAL(KIND=r8)    :: rkapp1
    REAL(KIND=r8)    :: albcp
    REAL(KIND=r8)    :: onebcp
    REAL(KIND=r8)    :: onebg
    REAL(KIND=r8)    :: cpbg
    REAL(KIND=r8)    :: twobal
    REAL(KIND=r8)    :: etamn
    INTEGER :: km1 
    INTEGER :: ic1 
    INTEGER :: l   
    INTEGER :: i   
    INTEGER :: j   
    REAL(KIND=r8)    :: tem1 
    INTEGER :: len1 
    INTEGER :: len2 
    INTEGER :: isav 
    INTEGER :: len11 
    INTEGER :: ii    
    REAL(KIND=r8)    :: tem   
    INTEGER :: lena  
    INTEGER :: lenb  
    INTEGER :: lena1 
    INTEGER :: ksl     

    rkapp1 = 1.0_r8  + rkap
    onebcp = 1.0_r8  / cp
    albcp  = alhl * onebcp
    onebg  = 1.0_r8  / grav
    cpbg   = cp   * onebg
    twobal = 2.0_r8 / alhl
    km1 = k  - 1
    ic1 = ic + 1
    !     
    !     settiing alf, bet, gam, prh, and pri : done only when setras=.t.
    !     
    IF (setras) THEN
       DO l=1,k
          DO i=1,lenc
             prh(i,l) = (prj(i,l+1)*prs(i,l+1) - prj(i,l)*prs(i,l)) &
                  / ((prs(i,l+1)-prs(i,l)) * rkapp1)
          END DO
       END DO

       DO l=1,k
          DO i=1,lenc
             tx5(i) = poi(i,l) * prh(i,l)
             tx1(i) = (prs(i,l) + prs(i,l+1)) * 0.5_r8
             tx3(i) = tx5(i)
             CALL qsat(tx3(i), tx1(i), tx2(i), tx4(i), .TRUE.)
             alf(i,l) = tx2(i) - tx4(i) * tx5(i)
             bet(i,l) = tx4(i) * prh(i,l)
             gam(i,l) = 1.0_r8 / ((1.0_r8 + tx4(i)*albcp) * prh(i,l))
             pri(i,l) = (cp/cmb2pa) / (prs(i,l+1) - prs(i,l))
          END DO
       END DO
    ENDIF

    DO l=1,k
       DO i=1,len
          tcu(i,l) = 0.0_r8
          qcu(i,l) = 0.0_r8
       END DO
    END DO

    DO i=1,lenc
       tx1(i)   = prj(i,k+1) * poi(i,k)
       qs1(i)   = alf(i,k) + bet(i,k)*poi(i,k)
       qol(i,k) = MIN(qs1(i)*rhmax,qoi(i,k))
       hol(i,k) = tx1(i)*cp + qol(i,k)*alhl
       eta(i,k) = 0.0e0_r8
       tx2(i)   = (prj(i,k+1) - prj(i,k)) * poi(i,k) * cp
    END DO

    IF (ic .LT. km1) THEN
       DO l=km1,ic1,-1
          DO i=1,lenc
             qs1(i)   = alf(i,l) + bet(i,l)*poi(i,l)
             qol(i,l) = MIN(qs1(i)*rhmax,qoi(i,l))
             tem1     = tx2(i) + prj(i,l+1) * poi(i,l) * cp 
             hol(i,l) = tem1 + qol(i,l )* alhl
             hst(i,l) = tem1 + qs1(i)   * alhl
             tx1(i)   = (prj(i,l+1) - prj(i,l)) * poi(i,l)
             eta(i,l) = eta(i,l+1) + tx1(i)*cpbg
             tx2(i)   = tx2(i)     + tx1(i)*cp
          END DO
       END DO
    ENDIF

    DO i=1,lenc
       hol(i,ic) = tx2(i)
       qs1(i)    = alf(i,ic) + bet(i,ic)*poi(i,ic)
       qol(i,ic) = MIN(qs1(i)*rhmax,qoi(i,ic)) 
       tem1      = tx2(i) + prj(i,ic1) * poi(i,ic) * cp 
       hol(i,ic) = tem1 + qol(i,ic) * alhl
       hst(i,ic) = tem1 + qs1(i)    * alhl
       tx3(i   ) = (prj(i,ic1) - prh(i,ic)) * poi(i,ic)
       eta(i,ic) = eta(i,ic1) + cpbg * tx3(i)
    END DO

    DO i=1,lenc
       tx2(i) = hol(i,k)  - hst(i,ic)
       tx1(i) = 0.0e0_r8
    END DO
    !     
    !     entrainment parameter
    !     
    DO l=ic,km1
       DO i=1,lenc
          tx1(i) = tx1(i) + (hst(i,ic) - hol(i,l)) *  &
               (eta(i,l) - eta(i,l+1))
       END DO
    END DO

    len1 = 0
    len2 = 0
    isav = 0

    DO i=1,lenc
       IF (tx1(i) .GT. 0.0e0_r8 .AND. tx2(i) .GT. 0.0e0_r8) THEN
          len1      = len1 + 1
          ia(len1)  = i
          alm(len1) = tx2(i) / tx1(i)
       ENDIF
    END DO

    len2 = len1

    IF (ic1 .LT. k) THEN
       DO i=1,lenc
          IF (tx2(i) .LE. 0.0_r8 .AND. (hol(i,k) .GT. hst(i,ic1))) THEN
             len2      = len2 + 1
             ia(len2)  = i
             alm(len2) = 0.0_r8
          ENDIF
       END DO
    ENDIF

    IF (len2 .EQ. 0) THEN 
       DO j=1,k
          DO i=1,lenc
             hst(i,j) = 0.0_r8
             qol(i,j) = 0.0_r8
          END DO
       END DO 
       DO i=1,lenc
          pcu(i) = 0.0_r8
       END DO
       RETURN
    ENDIF
    len11 = len1 + 1
    !     
    !     normalized massflux
    !     
    DO i=1,len2
       eta(i,k) = 1.0_r8
       ii       = ia(i)
       tx2(i)   = 0.5_r8 * (prs(ii,ic) + prs(ii,ic1))
       tx4(i)   = prs(ii,k)
    END DO

    DO i=len11,len2
       wfn(i)   = 0.0_r8
       ii       = ia(i)
       IF (hst(ii,ic1) .LT. hst(ii,ic)) THEN
          tx6(i) = (hst(ii,ic1)-hol(ii,k))/(hst(ii,ic1)-hst(ii,ic))
       ELSE
          tx6(i) = 0.0_r8
       ENDIF
       tx2(i) = 0.5_r8 * (prs(ii,ic1)+prs(ii,ic1+1)) * (1.0_r8-tx6(i)) &
            + tx2(i)      * tx6(i)
    END DO

    CALL acritn(len2, tx2, tx4, tx3)

    DO l=km1,ic,-1
       DO i=1,len2
          tx1(i) = eta(ia(i),l)
       END DO

       DO i=1,len2
          eta(i,l) = 1.0_r8 + alm(i) * tx1(i)
       END DO
    END DO
    !     
    !     cloud workfunction
    !     
    IF (len1 .GT. 0) THEN
       DO i=1,len1
          ii = ia(i)
          wfn(i) = - gam(ii,ic) * (prj(ii,ic1) - prh(ii,ic)) &
               *  hst(ii,ic) * eta(i,ic1)
       END DO
    ENDIF

    DO i=1,len2
       ii = ia(i)
       tx1(i) = hol(ii,k)
    END DO

    IF (ic1 .LE. km1) THEN
       DO l=km1,ic1,-1
          DO i=1,len2
             ii = ia(i)
             tem = tx1(i) + (eta(i,l) - eta(i,l+1)) * hol(ii,l)
             pcu(i) = prj(ii,l+1) - prh(ii,l)
             tem1   = eta(i,l+1) * pcu(i)
             tx1(i) = tx1(i)*pcu(i)
             pcu(i) = prh(ii,l) - prj(ii,l)
             tem1   = (tem1 + eta(i,l) * pcu(i)) * hst(ii,l)
             tx1(i) = tx1(i) + tem*pcu(i)
             wfn(i) = wfn(i) + (tx1(i) - tem1) * gam(ii,l)
             tx1(i) = tem
          END DO
       END DO
    ENDIF
    lena = 0
    IF (len1 .GT. 0) THEN
       !cdir nodep
       DO i=1,len1
          ii = ia(i)
          wfn(i) = wfn(i) + tx1(i) * gam(ii,ic) *  &
               (prj(ii,ic1)-prh(ii,ic)) - tx3(i)
          IF (wfn(i) .GT. 0.0_r8) THEN
             lena = lena + 1
             i1(lena) = ia(i)
             i2(lena) = i
             tx1(lena) = wfn(i)
             tx2(lena) = qs1(ia(i))
             tx6(lena) = 1.0_r8
          ENDIF
       END DO
    ENDIF
    lenb = lena
    DO i=len11,len2
       wfn(i) = wfn(i) - tx3(i)
       IF (wfn(i) .GT. 0.0_r8 .AND. tx6(i) .GT. 0.0_r8) THEN
          lenb = lenb + 1
          i1(lenb)  = ia(i)
          i2(lenb)  = i
          tx1(lenb) = wfn(i)
          tx2(lenb) = qs1(ia(i))
          tx4(lenb) = tx6(i)
       ENDIF
    END DO
    IF (lenb .LE. 0) THEN
       DO j=1,k
          DO i=1,lenc
             hst(i,j) = 0.0_r8
             qol(i,j) = 0.0_r8
          END DO
       END DO
       DO i=1,lenc
          pcu(i) = 0.0_r8
       END DO
       RETURN
    ENDIF
    DO i=1,lenb
       wfn(i) = tx1(i)
       qs1(i) = tx2(i)
    END DO
    DO l=ic,k
       DO i=1,lenb
          tx1(i) = eta(i2(i),l)
       END DO
       DO i=1,lenb
          eta(i,l) = tx1(i)
       END DO
    END DO
    lena1 = lena + 1
    DO i=1,lena
       ii = i1(i)
       tx8(i) = hst(ii,ic) - hol(ii,ic)
    END DO
    DO i=lena1,lenb
       ii = i1(i)
       tx6(i) = tx4(i)
       tem    = tx6(i) * (hol(ii,ic)-hol(ii,ic1)) + hol(ii,ic1)
       tx8(i) = hol(ii,k) - tem
       tem1   = tx6(i) * (qol(ii,ic)-qol(ii,ic1)) + qol(ii,ic1)
       tx5(i) = tem    - tem1 * alhl
       qs1(i) = tem1   + tx8(i)*(1.0e0_r8/alhl)
       tx3(i) = hol(ii,ic)
    END DO
    DO i=1,lenb
       ii = i1(i)
       wlq(i) = qol(ii,k) - qs1(i)     * eta(i,ic)
       uht(i) = uoi(ii,k) - uoi(ii,ic) * eta(i,ic)
       vht(i) = voi(ii,k) - voi(ii,ic) * eta(i,ic)
       tx7(i) = hol(ii,k)
    END DO
    DO l=km1,ic,-1
       DO i=1,lenb
          ii = i1(i)
          tem    = eta(i,l) - eta(i,l+1)
          wlq(i) = wlq(i) + tem * qol(ii,l)
          uht(i) = uht(i) + tem * uoi(ii,l)
          vht(i) = vht(i) + tem * voi(ii,l)
       END DO
    END DO
    !     
    !     calculate gs and part of akm (that requires eta)
    !     
    DO i=1,lenb
       ii = i1(i)
       tem        = (poi(ii,km1) - poi(ii,k)) /  &
            (prh(ii,k) - prh(ii,km1))
       hol(i,k)   = tem *(prj(ii,k)-prh(ii,km1))*prh(ii,k)*pri(ii,k)
       hol(i,km1) = tem *(prh(ii,k)-prj(ii,k))*prh(ii,km1)*pri(ii,km1)
       akm(i)     = 0.0e0_r8
       tx2(i)     = 0.5_r8 * (prs(ii,ic) + prs(ii,ic1))
    END DO

    IF (ic1 .LE. km1) THEN
       DO l=km1,ic1,-1
          DO i=1,lenb
             ii = i1(i)
             tem      = (poi(ii,l-1) - poi(ii,l)) * eta(i,l) &
                  / (prh(ii,l) - prh(ii,l-1))
             hol(i,l)   = tem * (prj(ii,l)-prh(ii,l-1)) * prh(ii,l) &
                  *  pri(ii,l)  + hol(i,l)
             hol(i,l-1) = tem * (prh(ii,l)-prj(ii,l)) * prh(ii,l-1) &
                  * pri(ii,l-1)
             akm(i)   = akm(i) - hol(i,l) &
                  * (eta(i,l)   * (prh(ii,l)-prj(ii,l)) + &
                  eta(i,l+1) * (prj(ii,l+1)-prh(ii,l))) / prh(ii,l)
          END DO
       END DO
    ENDIF

    CALL rncl(lenb, tx2, tx1, clf)

    DO i=1,lenb
       tx2(i) = (1.0e0_r8 - tx1(i)) * wlq(i)
       wlq(i) = tx1(i) * wlq(i)
       tx1(i) = hol(i,ic)
    END DO

    DO i=lena1, lenb
       ii = i1(i)
       tx1(i) = tx1(i) + (tx5(i)-tx3(i)+qol(ii,ic)*alhl)* &
            (pri(ii,ic)/cp)
    END DO

    DO i=1,lenb
       hol(i,ic) = tx1(i) - tx2(i) * albcp * pri(i1(i),ic)
    END DO

    IF (lena .GT. 0) THEN
       DO i=1,lena
          ii = i1(i)
          akm(i) = akm(i) - eta(i,ic1) * (prj(ii,ic1) - prh(ii,ic))  &
               * tx1(i) / prh(ii,ic)
       END DO
    ENDIF
    !     
    !     calculate gh
    !     
    DO i=1,lenb
       ii = i1(i)
       tx3(i)   =  qol(ii,km1) - qol(ii,k)
       gmh(i,k) = hol(i,k) + tx3(i) * pri(ii,k) * (albcp*0.5e0_r8)
       akm(i)   = akm(i) + gam(ii,km1)*(prj(ii,k)-prh(ii,km1))  &
            * gmh(i,k)
    END DO

    IF (ic1 .LE. km1) THEN
       DO l=km1,ic1,-1
          DO i=1,lenb
             ii = i1(i)
             tx2(i) = tx3(i)
             tx3(i) = (qol(ii,l-1) - qol(ii,l)) * eta(i,l)
             tx2(i) = tx2(i) + tx3(i)
             gmh(i,l) = hol(i,l) + tx2(i)   * pri(ii,l) * (albcp*0.5e0_r8)
          END DO
       END DO
    ENDIF

    DO i=lena1,lenb
       tx3(i) = tx3(i) + twobal &
            * (tx7(i) - tx8(i) - tx5(i) - qol(i1(i),ic)*alhl)
    END DO

    DO i=1,lenb
       gmh(i,ic) = tx1(i) + pri(i1(i),ic) * onebcp &
            * (tx3(i)*(alhl*0.5e0_r8) + eta(i,ic) * tx8(i))
    END DO
    !     
    !     calculate hc part of akm
    !     
    IF (ic1 .LE. km1) THEN
       DO i=1,lenb
          tx1(i) = gmh(i,k)
       END DO
       DO l=km1,ic1,-1
          DO i=1,lenb
             ii = i1(i)
             tx1(i) = tx1(i) + (eta(i,l) - eta(i,l+1)) * gmh(i,l)
             tx2(i) = gam(ii,l-1) * (prj(ii,l) - prh(ii,l-1))
          END DO
          IF (l .EQ. ic1) THEN
             DO i=lena1,lenb
                tx2(i) = 0.0e0_r8
             END DO
          ENDIF
          DO i=1,lenb
             ii = i1(i)
             akm(i) = akm(i) + tx1(i) *  &
                  (tx2(i) + gam(ii,l)*(prh(ii,l)-prj(ii,l)))
          END DO
       END DO
    ENDIF

    DO i=lena1,lenb
       ii = i1(i)
       tx2(i) = 0.5_r8 * (prs(ii,ic) + prs(ii,ic1)) &
            + 0.5_r8*(prs(ii,ic+2) - prs(ii,ic)) * (1.0e0_r8-tx6(i))
       tx1(i) = prs(ii,ic1)
       tx5(i) = 0.5_r8 * (prs(ii,ic1) + prs(ii,ic+2))
       IF ((tx2(i) .GE. tx1(i)) .AND. (tx2(i) .LT. tx5(i))) THEN
          tx6(i)     = 1.0e0_r8 - (tx2(i) - tx1(i)) / (tx5(i) - tx1(i))
          tem        = pri(ii,ic1) / pri(ii,ic)
          hol(i,ic1) = hol(i,ic1) + hol(i,ic) * tem
          hol(i,ic)  = 0.0e0_r8
          gmh(i,ic1) = gmh(i,ic1) + gmh(i,ic) * tem
          gmh(i,ic)  = 0.0e0_r8
       ELSEIF (tx2(i) .LT. tx1(i)) THEN
          tx6(i) = 1.0_r8
       ELSE
          tx6(i) = 0.0_r8
       ENDIF
    END DO

    DO i=1,lenc
       pcu(i) = 0.0_r8
    ENDDO

    DO i=1,lenb
       ii = i1(i)
       IF (akm(i) .LT. 0.0e0_r8 .AND. wlq(i) .GE. 0.0_r8) THEN
          wfn(i) = - tx6(i) * wfn(i) * rasalf / akm(i)
       ELSE
          wfn(i) = 0.0e0_r8
       ENDIF
       tem       = (prs(ii,k+1)-prs(ii,k))*(cmb2pa*frac)
       wfn(i)    = MIN(wfn(i), tem)
       !
       !     compute cloud amount
       !
       etamn=0.0_r8
       IF(km1.GT.ic) THEN
          DO ksl= km1,ic,-1
             etamn=etamn+eta(i,ksl)
          END DO
          etamn=etamn/real(km1-ic,kind=r8)
          tx1(i)=wfn(i)*864.0_r8*etamn
       ELSE
          tx1(i)=wfn(i)*864.0_r8
       ENDIF
       !
       !     precipitation
       !
       pcu(ii) =  wlq(i) * wfn(i) * onebg
       !     
       !     cumulus friction at the bottom layer
       !     
       tx4(i)   = wfn(i) * (1.0_r8/alhl)
       tx5(i)   = wfn(i) * onebcp
    END DO

    DO i=1,lenb
       ii = i1(i)
       cln(ii) = tx1(i)
    END DO
    !     
    !     theta and q change due to cloud type ic
    !     
    DO l=ic,k
       DO i=1,lenb
          ii = i1(i)
          tem       = (gmh(i,l) - hol(i,l)) * tx4(i)
          tem1      =  hol(i,l) * tx5(i)
          tcu(ii,l) = tem1 / prh(ii,l)
          qcu(ii,l) = tem
       END DO
    END DO

    DO l=1,k
       DO i=1,lenc
          hst(i,l) = 0.0_r8
          qol(i,l) = 0.0_r8
       END DO
    END DO
  END SUBROUTINE cloud



  ! ras     :relaxed arakawa-schubert.



  SUBROUTINE ras( &
       len   ,lenc  ,k     ,dt    ,ncrnd ,krmax ,frac  , &
       rasal , botop,alhl  ,rkap  ,poi   , qoi  ,uoi   ,voi   , &
       prs   ,prj   ,cup   ,cln   ,q1    ,q2    ,alf   ,bet   ,gam   , &
       prh   ,pri   ,hoi   ,eta   ,tcu   ,qcu   ,hst   ,qol   ,gmh   , &
       tx1   ,tx2   ,tx3   ,tx4   ,tx5   ,tx6   ,tx7   ,tx8   ,tx9   , &
       wfn   ,akm   ,qs1   ,clf   ,uht   ,vht   ,wlq   ,pcu   ,ia    , &
       i1    ,i2    ,kmax  ,kmaxm1,kmaxp1)

    !==========================================================================
    !
    ! ras     :relaxed arakawa-schubert.
    !==========================================================================
    ! kmax......Number of sigma levels       
    ! kmaxm1....Parameter kmaxm1 = kmax-1
    ! kmaxp1....Parameter kmaxp1 = kmax+1
    ! icm.......Parameter icm    = 100      
    ! daylen....Parameter daylen = 86400.0_r8
    ! cp........specific heat of air           (j/kg/k) 
    ! grav......gas constant of dry air        (j/kg/k)
    ! cmb2pa....Parameter cmb2pa = 100.0_r8
    ! rhmax.....Parameter rhmax  = 0.9999_r8 
    ! ncols......Number of grid points on a gaussian latitude circle   
    ! actp......Parameter actp   = 1.7_r8 
    ! facm......Parameter facm   = 1.0_r8   
    ! p5........Parameter p5     = 500.0_r8     
    ! p8........Parameter p8     = 800.0_r8         
    ! pt8.......Parameter pt8    = 0.8_r8  
    ! pt2.......Parameter pt2    = 0.2_r8
    ! pfac......Parameter pfac   = pt2/(p8-p5)    
    ! cucld.....Parameter cucld  = 0.5_r8  
    ! len.......ncols
    ! lenc......ncols
    ! k.........kmax
    ! rkap......rkap=r/cp=287.05e0_r8/1004.6e0_r8
    ! alhl......alhl=2.52e6_r8
    ! dt........time interval,usually =delt,but changes
    !           in nlnmi (dt=1.) and at dead start(delt/4,delt/2) 
    ! frac......frac=1.
    ! krmax.....=nls
    ! nls.......Number of layers in the stratosphere.    
    ! ncrnd.....ncrnd=0
    ! poi.......array pot. temp.  
    ! qoi.......array humidity.
    ! prs  
    ! prj  
    ! uoi  
    ! voi  
    ! rasal
    ! q1.......heating diagnostics 
    ! q2.......moitening diagnostics       
    ! cln  
    ! cup......rainfall diagnostic (mm/day)
    ! tcu  
    ! qcu  
    ! alf  
    ! bet  
    ! gam  
    ! gmh  
    ! eta  
    ! hoi  
    ! hst  
    ! qol  
    ! prh  
    ! pri  
    ! tx1  
    ! tx2  
    ! tx3  
    ! tx4  
    ! tx5  
    ! tx6  
    ! tx7  
    ! tx8  
    ! tx9  
    ! wfn  
    ! akm  
    ! qs1  
    ! wlq  
    ! pcu  
    ! uht  
    ! vht  
    ! clf  
    ! ia   
    ! i1   
    ! i2   
    ! botop.......botop=.true.
    !==========================================================================
    INTEGER, INTENT(in   ) :: kmax  
    INTEGER, INTENT(in   ) :: kmaxm1
    INTEGER, INTENT(in   ) :: kmaxp1

    INTEGER, INTENT(in   ) :: len
    INTEGER, INTENT(in   ) :: lenc
    INTEGER, INTENT(inout) :: k
    REAL(KIND=r8),    INTENT(in   ) :: rkap
    REAL(KIND=r8),    INTENT(in   ) :: alhl
    REAL(KIND=r8),    INTENT(in   ) :: dt
    REAL(KIND=r8),    INTENT(in   ) :: frac
    INTEGER, INTENT(in   ) :: krmax
    INTEGER, INTENT(in   ) :: ncrnd

    REAL(KIND=r8),    INTENT(inout) :: poi  (len,k)
    REAL(KIND=r8),    INTENT(inout) :: qoi  (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: prs  (len,k+1) 
    REAL(KIND=r8),    INTENT(inout) :: prj  (len,k+1)
    REAL(KIND=r8),    INTENT(in   ) :: uoi  (len,k)
    REAL(KIND=r8),    INTENT(in   ) :: voi  (len,k)
    REAL(KIND=r8),    INTENT(in   ) :: rasal(k-1)

    REAL(KIND=r8),    INTENT(inout) :: q1   (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: q2   (len,k) 
    REAL(KIND=r8),    INTENT(inout  ) :: cln  (len,k)
    REAL(KIND=r8),    INTENT(inout) :: cup  (len,k)
    REAL(KIND=r8),    INTENT(inout) :: tcu  (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: qcu  (len,k)

    REAL(KIND=r8),    INTENT(inout) :: alf  (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: bet  (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: gam  (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: gmh  (lenc,k)
    REAL(KIND=r8),    INTENT(inout) :: eta  (lenc,k) 
    REAL(KIND=r8),    INTENT(inout) :: hoi  (lenc,k) 
    REAL(KIND=r8),    INTENT(inout) :: hst  (lenc,k) 
    REAL(KIND=r8),    INTENT(inout) :: qol  (lenc,k)
    REAL(KIND=r8),    INTENT(inout) :: prh  (len,k) 
    REAL(KIND=r8),    INTENT(inout) :: pri  (len,k)

    REAL(KIND=r8),    INTENT(inout) :: tx1  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx2  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx3  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx4  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx5  (lenc)
    REAL(KIND=r8),    INTENT(inout) :: tx6  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx7  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx8  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: tx9  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: wfn  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: akm  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: qs1  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: wlq  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: pcu  (lenc)
    REAL(KIND=r8),    INTENT(inout) :: uht  (lenc) 
    REAL(KIND=r8),    INTENT(inout) :: vht  (lenc) 
    REAL(KIND=r8),    INTENT(inout  ) :: clf  (lenc)
    INTEGER, INTENT(inout) :: ia   (lenc)   
    INTEGER, INTENT(inout) :: i1   (lenc) 
    INTEGER, INTENT(inout) :: i2   (lenc)
    LOGICAL, INTENT(in   ) :: botop    

    REAL(KIND=r8)   , PARAMETER :: cmb2pa=100.0_r8
    REAL(KIND=r8)   , PARAMETER :: rhmax =0.9999_r8
    INTEGER, PARAMETER :: icm   =100 
    REAL(KIND=r8)   , PARAMETER :: daylen=86400.0_r8

    REAL(KIND=r8)                   :: rkapp1
    REAL(KIND=r8)                   :: onebcp
    REAL(KIND=r8)                   :: albcp
    INTEGER                :: km1
    REAL(KIND=r8)                   :: onbdt
    REAL(KIND=r8)                   :: fracs
    INTEGER                :: kcr
    INTEGER                :: kfx
    INTEGER                :: ncmx
    INTEGER                :: nc
    INTEGER                :: ib
    REAL(KIND=r8)                   :: rasalf
    INTEGER                :: ic(icm)
    LOGICAL                :: setras
    INTEGER                :: i
    INTEGER                :: l

    setras = .FALSE.
    rkapp1=1.0_r8 + rkap
    onebcp = 1.0_r8 / cp
    albcp = alhl * onebcp
    !
    !     dgd modify k value
    !
    k=kmaxm1
    DO i=1,lenc
       poi(i,kmaxm1)=poi(i,kmax)*(prs(i,kmaxp1)-prs(i,kmax))+ &
            poi(i,kmaxm1)*(prs(i,kmax)-prs(i,kmaxm1))
       poi(i,kmaxm1)=poi(i,kmaxm1)/(prs(i,kmaxp1)-prs(i,kmaxm1))
       qoi(i,kmaxm1)=qoi(i,kmax)*(prs(i,kmaxp1)-prs(i,kmax))+ &
            qoi(i,kmaxm1)*(prs(i,kmax)-prs(i,kmaxm1))
       qoi(i,kmaxm1)=qoi(i,kmaxm1)/(prs(i,kmaxp1)-prs(i,kmaxm1))
    END DO

    DO i=1,lenc
       prj(i,kmaxp1)=prj(i,kmaxp1)
       prj(i,kmax)=prj(i,kmaxp1)
       prs(i,kmaxp1)=prs(i,kmaxp1)
       prs(i,kmax)=prs(i,kmaxp1)
    END DO
    DO l=1,k
       DO i=1,lenc
          prh(i,l) = (prj(i,l+1)*prs(i,l+1) - prj(i,l)*prs(i,l)) &
               / ((prs(i,l+1)-prs(i,l)) * rkapp1)
       END DO
    END DO
    DO l=1,k
       DO i=1,lenc
          tx5(i) = poi(i,l) * prh(i,l)
          tx1(i) = (prs(i,l) + prs(i,l+1)) * 0.5_r8
          tx3(i) = tx5(i)
          CALL qsat(tx3(i), tx1(i), tx2(i), tx4(i), .TRUE.)
          alf(i,l) = tx2(i) - tx4(i) * tx5(i)
          bet(i,l) = tx4(i) * prh(i,l)
          gam(i,l) = 1.0_r8 / ((1.0_r8 + tx4(i)*albcp) * prh(i,l))
          pri(i,l) = (cp/cmb2pa) / (prs(i,l+1) - prs(i,l))
       END DO
    END DO
    !
    !     done modification
    !
    km1    = k  - 1
    onbdt  = 1.0_r8 / dt
    fracs  = frac  * onbdt
    !
    !     set number of clouds to adjust during this call to ras, ncmx,
    !     and the cloud calling sequence, ic.  this allows various
    !     combinations of randomly and sequentially called clouds.
    !
    kcr   = MIN(km1,krmax)
    kfx   = km1 - kcr
    ncmx  = kfx + ncrnd
    IF (kfx .GT. 0) THEN
       IF (botop) THEN
          DO nc=1,kfx
             ic(nc) = k - nc
          END DO
       ELSE
          DO nc=kfx,1,-1
             ic(nc) = k - nc
          END DO
       END IF
    END IF
    !
    !     this area commented until machine independent random number
    !     generator can be found.  parameters set in RunCu_RelAraSch to use 
    !     non random clouds
    !     
    !     IF (ncrnd .gt. 0) THEN
    !     cray rng setup
    !     call ranset(iseed)
    !     dec rng setup
    !     call srand(iseed)
    !     do 30 i=1,ncrnd
    !     cray rng
    !     irnd = (ranf()-0.0005_r8)*(kcr-krmin+1)
    !     dec rng
    !     irnd = (rand()-0.0005_r8)*(kcr-krmin+1)
    !     ic(kfx+i) = irnd + krmin
    !     30    continue
    !     END IF
    !     
    !     loop over clouds to be adjusted during this call
    !
    DO nc=1,ncmx
       ib = ic(nc)
       rasalf = rasal(ib) * onbdt
       CALL cloud( &
            len   ,lenc  ,k     ,ib    ,rasalf,setras,fracs ,alhl  ,rkap  , &
            poi   ,qoi   ,uoi   ,voi   ,prs   ,prj   ,pcu   ,cln(1,ib),tcu, &
            qcu   ,alf   ,bet   ,gam   ,prh   ,pri   ,hoi   ,eta   ,hst   , &
            qol   ,gmh   ,tx1   ,tx2   ,tx3   ,tx4   ,tx5   ,tx6   ,tx7   , &
            tx8   ,tx9   ,wfn   ,akm   ,qs1   ,clf   ,uht   ,vht   ,wlq   , &
            ia    ,i1    ,i2    ,cmb2pa,rhmax )
       DO l=ib,k
          DO i=1,lenc
             !
             !     update pot. temp. and humidity.
             !
             poi(i,l) = poi(i,l) + tcu(i,l) * dt
             qoi(i,l) = qoi(i,l) + qcu(i,l) * dt
             !
             !     heating and moitening diagnostics
             !
             q1(i,l)  = q1(i,l)  + tcu(i,l) * prh(i,l) * daylen
             q2(i,l)  = q2(i,l)  + qcu(i,l) * daylen
          END DO
       END DO
       !
       !     rainfall diagnostic (mm/day)
       !
       DO i=1,lenc
          cup(i,ib) = cup(i,ib) + pcu(i) * daylen
       END DO
    END DO
  END SUBROUTINE ras




  ! RunCu_RelAraSch :used to interface with the relaxed arakawa schubert code of 
  !         moorthi and suarez.



  SUBROUTINE RunCu_RelAraSch (dtwrk ,dqwrk ,sl    ,si    ,fpn   ,ktop  ,kbot  ,rrr   , &
       hrar  ,qrar  ,dt    ,ftn   ,fqn   ,del   ,kuo   ,cldm  , &
       cflric,kmaxp1,kmaxm1,ncols  ,kmax  ,nls   )
    !
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle     
    ! kmax......Number of sigma levels     
    ! nls.......Number of layers in the stratosphere.    
    ! nlst......Parameter nlst   = 01    
    ! kmaxm1....Parameter kmaxm1 = kmax-1
    ! kmaxp1....Parameter kmaxp1 = kmax+1
    ! icm.......Parameter icm    = 100 
    ! daylen....Parameter daylen = 86400.0 
    ! cmb2pa....Parameter cmb2pa = 100.0
    ! rhmax.....Parameter rhmax  = 0.9999 
    ! dt........time interval,usually =delt,but changes
    !           in nlnmi (dt=1.) and at dead start(delt/4,delt/2)   
    ! cflric....parameter used by relaxed arakawa-schubert
    ! fpn.......sfc pres (cb)
    ! sl........sigma value at midpoint of
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
    ! dtwrk
    ! dqwrk
    ! ktop......ktop (g.t.e 1 and l.t.e. km) is highest lvl for which
    !           tin is colder than moist adiabat given by the
    !           (ktop=1 denotes tin(k=2) is already g.t.e. moist adb.)
    !           allowance is made for perhaps one level below ktop
    !           at which tin was warmer than tmst. 
    ! kbot......is the first regular level above the lcl    
    ! ftn.......temperature field in the spectral grid
    ! fqn.......specific humidit field in the spectral grid  
    ! hrar......this array is needed for the heating from ras scheme
    ! qrar......this array is needed for the mostening from ras scheme 
    ! rrr  
    ! del ......sigma spacing for each layer computed in routine "setsig". 
    ! kuo.......flag to indicate that deep convection was done
    !           kuo, ktop and plcl are longitude arrays       
    ! cldm 
    !==========================================================================
    INTEGER, INTENT(in   ) :: kmaxp1
    INTEGER, INTENT(in   ) :: kmaxm1
    INTEGER, INTENT(in   ) :: ncols    
    INTEGER, INTENT(in   ) :: kmax    
    INTEGER, INTENT(in   ) :: nls     
    !
    !     default parameter statements for dimensioning resolution
    !             
    REAL(KIND=r8),    INTENT(in   ) :: dt
    REAL(KIND=r8),    INTENT(in   ) :: cflric
    REAL(KIND=r8),    INTENT(in   ) :: fpn  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sl   (kmax)
    REAL(KIND=r8),    INTENT(in   ) :: si   (kmaxp1)
    REAL(KIND=r8),    INTENT(in   ) :: dtwrk(ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: dqwrk(ncols,kmax)
    INTEGER, INTENT(inout) :: ktop (ncols)
    INTEGER, INTENT(inout) :: kbot (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ftn  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: fqn  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: hrar (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: qrar (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: rrr  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: del  (kmax)
    INTEGER, INTENT(inout) :: kuo  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cldm (ncols)
    !     
    !     set up variables
    !
    REAL(KIND=r8)                   :: prs  (ncols,kmaxp1)
    REAL(KIND=r8)                   :: prj  (ncols,kmaxp1)
    REAL(KIND=r8)                   :: rasal(kmaxm1)
    REAL(KIND=r8)                   :: poi  (ncols,kmax)
    REAL(KIND=r8)                   :: qoi  (ncols,kmax)
    REAL(KIND=r8)                   :: uoi  (ncols,kmax)
    REAL(KIND=r8)                   :: voi  (ncols,kmax)
    REAL(KIND=r8)                   :: rrr1 (ncols)  
    REAL(KIND=r8)                   :: cln  (ncols,kmax)
    REAL(KIND=r8)                   :: q1   (ncols,kmax)
    REAL(KIND=r8)                   :: q2   (ncols,kmax)
    REAL(KIND=r8)                   :: tcu  (ncols,kmax)
    REAL(KIND=r8)                   :: qcu  (ncols,kmax)
    REAL(KIND=r8)                   :: rnew (ncols)
    REAL(KIND=r8)                   :: rdum (ncols)
    REAL(KIND=r8)                   :: rvd  (ncols,kmax)
    !     
    !     work arrays needed for relaxed arakawa-schubert
    !        
    REAL(KIND=r8)                 :: alf(ncols,kmax) 
    REAL(KIND=r8)                 :: bet(ncols,kmax) 
    REAL(KIND=r8)                 :: gam(ncols,kmax)
    REAL(KIND=r8)                 :: gmh(ncols,kmax)
    REAL(KIND=r8)                 :: eta(ncols,kmax) 
    REAL(KIND=r8)                 :: hoi(ncols,kmax) 
    REAL(KIND=r8)                 :: hst(ncols,kmax)
    REAL(KIND=r8)                 :: qol(ncols,kmax)
    REAL(KIND=r8)                 :: prh(ncols,kmax) 
    REAL(KIND=r8)                 :: pri(ncols,kmax)

    REAL(KIND=r8)                 :: tx1(ncols) 
    REAL(KIND=r8)                 :: tx2(ncols)
    REAL(KIND=r8)                 :: tx3(ncols) 
    REAL(KIND=r8)                 :: tx4(ncols) 
    REAL(KIND=r8)                 :: tx5(ncols)
    REAL(KIND=r8)                 :: tx6(ncols) 
    REAL(KIND=r8)                 :: tx7(ncols) 
    REAL(KIND=r8)                 :: tx8(ncols) 
    REAL(KIND=r8)                 :: tx9(ncols)
    REAL(KIND=r8)                 :: wfn(ncols) 
    REAL(KIND=r8)                 :: akm(ncols) 
    REAL(KIND=r8)                 :: qs1(ncols) 
    REAL(KIND=r8)                 :: wlq(ncols) 
    REAL(KIND=r8)                 :: pcu(ncols)
    REAL(KIND=r8)                 :: uht(ncols) 
    REAL(KIND=r8)                 :: vht(ncols) 
    REAL(KIND=r8)                 :: clf(ncols)
    INTEGER              :: ia (ncols)   
    INTEGER              :: i1 (ncols) 
    INTEGER              :: i2 (ncols) 
    LOGICAL              :: botop
    INTEGER              :: len  
    INTEGER              :: lenc   
    INTEGER              :: k1  
    REAL(KIND=r8)                 :: twodt  
    INTEGER              :: krmin  
    INTEGER              :: krmax  
    INTEGER              :: ncrnd  
    REAL(KIND=r8)                 :: alhl   
    REAL(KIND=r8)                 :: rkap   
    REAL(KIND=r8)                 :: frac   
    INTEGER              :: k  
    REAL(KIND=r8)                 :: cnst   
    INTEGER              :: i  

    INTEGER, PARAMETER     ::  nlst   =01

    !
    !     define constants
    !
    botop=.TRUE.
    len=ncols
    lenc=ncols
    k1=kmax
    twodt=2.0_r8*dt
    krmin=nls
    krmax=krmin
    ncrnd=0
    alhl=2.52e6_r8
    rkap=gasr/cp
    frac=1.0_r8
    !     
    !     define variables 
    !     
    !     now need to set up arrays for call to arakawa -schubert and 
    !     flip over since the subroutine uses an inverted vertical
    !     coordinate compared to that in the rest of the model
    !     
    !     set up arrays
    !
    DO k=1,kmaxp1
       cnst=si(k)**rkap
       DO i=1,ncols
          prj(i,kmax+2-k)=cnst
          prs(i,kmax+2-k)=si(k)*fpn(i)*10.0_r8
       END DO
    END DO
    DO  k=1,kmax
       cnst=((1.0_r8/sl(k))**rkap)
       DO i=1,ncols
          poi(i,kmax+1-k)=dtwrk(i,k)*cnst
          qoi(i,kmax+1-k)=dqwrk(i,k)
          hrar(i,k)=0.0_r8
          qrar(i,k)=0.0_r8
          uoi(i,k)=0.0_r8
          voi(i,k)=0.0_r8
          cln(i,k)=0.0_r8
          q1(i,k)=0.0_r8
          q2(i,k)=0.0_r8
          tcu(i,k)=0.0_r8
          qcu(i,k)=0.0_r8
          alf(i,k)=0.0_r8
          bet(i,k)=0.0_r8
          gam(i,k)=0.0_r8
          gmh(i,k)=0.0_r8
          eta(i,k)=0.0_r8
          hoi(i,k)=0.0_r8
          hst(i,k)=0.0_r8
          qol(i,k)=0.0_r8 
          prh(i,k)=0.0_r8
          pri(i,k)=0.0_r8
          rvd(i,k)=0.0_r8
       END DO
    END DO
    DO k=1,kmaxm1
       rasal(k)=twodt/7200.0_r8
    END DO
    DO i=1,ncols
       rrr(i)=0.0_r8
       rrr1(i)=0.0_r8
       ktop(i)=3
       kbot(i)=3
       tx9(i)=0.0_r8
       !
       !     needed for cloud fraction based on mass flux
       !
       cldm(i)=0.0_r8
    END DO
    !     
    !     call relaxed arakawa-schubert
    !     
    CALL ras( &
         len   ,lenc  ,k1    ,twodt ,ncrnd ,krmax ,frac  , &
         rasal ,botop ,alhl  ,rkap  ,poi   ,qoi   ,uoi   ,voi   , &
         prs   ,prj   ,rvd   ,cln   ,q1  ,q2 ,alf   ,bet   ,gam   , &
         prh   ,pri   ,hoi   ,eta   ,tcu   ,qcu ,hst   ,qol   ,gmh   , &
         tx1   ,tx2   ,tx3   ,tx4   ,tx5   ,tx6 ,tx7   ,tx8   ,tx9   , &
         wfn   ,akm   ,qs1   ,clf   ,uht   ,vht ,wlq   ,pcu   ,ia    , &
         i1    ,i2    ,kmax  ,kmaxm1,kmaxp1)
    !     
    !     now need to assign and flip output arrays
    !     
    DO k=1,kmax
       DO i=1,ncols
          IF(cln(i,kmaxp1-k).GT.0.0_r8) THEN
             ktop(i)=k
             !
             !     needed for cloud fraction based on mass flux
             !
             cldm(i)=cldm(i)+cln(i,kmaxp1-k)
          ELSE
             ktop(i)=ktop(i)
             cldm(i)=cldm(i)
          END IF
       END DO
    END DO
    DO k=1,kmax
       DO i=1,ncols
          hrar(i,k)=q1(i,kmaxp1-k)/86400.0_r8
          qrar(i,k)=q2(i,kmaxp1-k)/86400.0_r8
          ftn(i,k)=dtwrk(i,k)+twodt*q1(i,kmaxp1-k)/86400.0_r8
          fqn(i,k)=dqwrk(i,k)+twodt*q2(i,kmaxp1-k)/86400.0_r8
          rvd(i,k)=twodt*rvd(i,k)/86400.0_r8
       END DO
    END DO
    DO i=1,ncols
       cldm(i)=0.035_r8*LOG(1.0_r8+cldm(i))
    END DO
    !
    !     now for subcloud layer
    !
    DO i=1,ncols
       hrar(i,1)=hrar(i,2)
       qrar(i,1)=qrar(i,2)
       ftn(i,1)=dtwrk(i,1)+twodt*hrar(i,2)
       fqn(i,1)=dqwrk(i,1)+twodt*qrar(i,2)
    END DO
    !
    !     now for re-evaporation of fallingrain
    !
    CALL rnevp(ftn   ,fqn   ,rvd   ,fpn   ,sl    ,del   ,cflric, &
         rnew  ,rdum  ,twodt ,ncols ,kmax  ,nlst)
    !
    !        now need to modify values for re-evaporation
    !
    DO i=1,ncols
       rrr(i)=0.5_r8*1.0_r8/1000.0_r8*rnew(i)
       IF(rrr(i).GT.0.0_r8) kuo(i)=1
    END DO

    DO k=1,kmax
       DO i=1,ncols
          hrar(i,k)=(ftn(i,k)-dtwrk(i,k))/twodt
          qrar(i,k)=(fqn(i,k)-dqwrk(i,k))/twodt
       END DO
    END DO
  END SUBROUTINE RunCu_RelAraSch




  SUBROUTINE shllcl(dt, ps, sl, si, qn1, tn1, kuo, &
       plcl, kktop, kkbot, ncols, kmax)
    !
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle  
    ! kmax......Number of sigma levels  
    ! msta
    ! dt........time interval,usually =delt,but changes
    !           in nlnmi (dt=1.) and at dead start(delt/4,delt/2)   
    ! ps........surface pressure      
    ! sl........sigma value at midpoint of
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
    ! qn1.......qn1 is q (specific humidit) at the n+1 time level 
    ! tn1.......tn1 is tmp (temperature) at the n+1 time level 
    ! kktop......ktop (g.t.e 1 and l.t.e. km) is highest lvl for which
    !            tin is colder than moist adiabat given by the
    !            (ktop=1 denotes tin(k=2) is already g.t.e. moist adb.)
    !            allowance is made for perhaps one level below ktop
    !            at which tin was warmer than tmst.
    ! kuo........flag to indicate that deep convection was done
    !            kuo, ktop and plcl are longitude arrays     
    ! plcl.......pressure at the lcl 
    ! kkbot......is the first regular level above the lcl   
    ! cp........Specific heat of air      (j/kg/k) 
    ! hl........heat of evaporation of water     (j/kg) 
    ! gasr......gas constant of dry air        (j/kg/k)  
    ! rmwmd.....fracao molar entre a agua e o ar       
    ! sthick....upper limit for originating air for lcl.  replaces kthick.      
    ! ki........lowest level from which parcels can be lifted to find lcl            
    !==========================================================================

    INTEGER, INTENT(IN   ) :: ncols
    INTEGER, INTENT(IN   ) :: kmax
    REAL(KIND=r8),    INTENT(IN   ) :: dt
    REAL(KIND=r8),    INTENT(IN   ) :: ps   (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: sl   (kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: si   (kmax+1)
    REAL(KIND=r8),    INTENT(IN   ) :: qn1  (ncols,kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: tn1  (ncols,kmax)
    INTEGER, INTENT(INOUT  ) :: kktop(ncols)
    INTEGER, INTENT(INOUT  ) :: kuo  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: plcl (ncols)
    INTEGER, INTENT(INOUT  ) :: kkbot(ncols)

    REAL(KIND=r8)    :: press (ncols,kmax)
    REAL(KIND=r8)    :: tin   (ncols,kmax)
    REAL(KIND=r8)    :: qin   (ncols,kmax)         
    REAL(KIND=r8)    :: tmst  (ncols,kmax)
    REAL(KIND=r8)    :: qmst  (ncols,kmax)
    REAL(KIND=r8)    :: tpar  (ncols)
    REAL(KIND=r8)    :: espar (ncols)
    REAL(KIND=r8)    :: qspar (ncols)
    REAL(KIND=r8)    :: qpar  (ncols)
    REAL(KIND=r8)    :: qex1  (ncols)
    REAL(KIND=r8)    :: tlcl  (ncols)
    REAL(KIND=r8)    :: qexces(ncols)
    REAL(KIND=r8)    :: dqdp  (ncols)
    REAL(KIND=r8)    :: deltap(ncols)
    REAL(KIND=r8)    :: hnew  (ncols)
    REAL(KIND=r8)    :: slcl  (ncols)
    INTEGER :: ier   (ncols)         
    LOGICAL :: llift (ncols)
    LOGICAL :: lconv (ncols)
    INTEGER :: ll    (ncols)       

    REAL(KIND=r8)    :: cappa
    REAL(KIND=r8)    :: rdt
    REAL(KIND=r8)    :: cpovl
    REAL(KIND=r8)    :: kk
    INTEGER :: i
    INTEGER :: k
    INTEGER :: kthick

    cappa=gasr/cp

    IF(dt.EQ.0.0e0_r8) RETURN
    !
    !     set default values
    !
    DO i=1,ncols
       kktop(i)=1
       kkbot(i)=1
       llift(i)=.FALSE.
       lconv(i)=.FALSE.
    END DO
    !     
    !     define kthick in terms of sigma values
    !     
    DO k=1,kmax
       IF(si(k).GE.sthick.AND.si(k+1).LT.sthick) THEN
          kthick = k
          EXIT
       ENDIF
    END DO
    rdt=1.0e0_r8/dt
    cpovl= cp/hl
    kk = kmax
    !     
    !     qn is q at the n-1 time level
    !     del=del sigma (del p ovr psfc)
    !     ps=sfc pres (cb)
    !     
    DO i=1, ncols
       hnew(i) = 1.0e-2_r8*ps(i)
    END DO

    DO  k=1,kmax
       DO  i=1,ncols
          press(i,k)=sl(k)*ps(i)
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

    DO i=1,ncols
       qex1(i) = 0.0e0_r8
       qpar(i) = qin(i,ki)
    END DO
    !     
    !     lift parcel from k=ki until it becomes saturated
    !
    DO k =ki,kthick
       DO i=1,ncols
          IF(.NOT.llift(i)) THEN
             tpar(i) = tin(i,ki)* &
                  EXP(cappa*LOG(press(i,k)/press(i,ki)))
             espar(i) = es3(tpar(i))
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
          ENDIF
       END DO
    END DO
    !     
    !     quit if parcel still unsat at k = kthick
    !     in this case,set low value of plcl as signal to shalmano
    !
    DO i=1,ncols
       IF(.NOT.llift(i)) THEN
          plcl(i) = 1.0e0_r8
          kuo  (i) = 5
       ENDIF
    END DO
    CALL mstad2(hnew, sl, tin, tmst, qmst,  kktop, ier, &
         slcl, ll, qin, tlcl, llift, ncols, kmax)
    !
    !     tmst and qmst,k=1...ktop contain cloud temps and specific humid.
    !     store values of plcl and ktop to pass to shalmano
    !
    DO i=1,ncols
       IF(llift(i)) THEN
          IF(ier(i).NE.0) THEN
             kuo(i) =6
             kktop(i) = 1
             plcl(i) = 10.e0_r8
             lconv(i)=.FALSE.
          ELSE
             lconv(i)=.TRUE.
          ENDIF
       ENDIF
    END DO
  END SUBROUTINE shllcl

END MODULE ModConRas
