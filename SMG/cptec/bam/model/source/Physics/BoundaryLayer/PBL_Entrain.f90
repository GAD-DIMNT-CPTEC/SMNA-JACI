
!  $Id: GEOS_Utilities.F90,v 1.16 2007/10/04 20:38:23 dasilva Exp $

!BOP

! !MODULE: GEOS_Utils -- A Module to containing computational utilities

MODULE GEOS_UtilsMod
  USE Constants, ONLY :     &
       cp,            &
       grav,          &
       gasr,          &
       r8

  ! !USES:

  !  use MAPL_ConstantsMod

  ! #include "qsatlqu.code"
  ! #include "qsatice.code"
  ! #include "esatlqu.code"
  ! #include "esatice.code"
  ! #include "trisolve.code"
  ! #include "trilu.code"

  IMPLICIT NONE
  SAVE
  PRIVATE


  ! !PUBLIC MEMBER FUNCTIONS:

  PUBLIC GEOS_QsatSet

  PUBLIC GEOS_QsatLQU
  PUBLIC GEOS_QsatICE
  PUBLIC GEOS_Qsat
  PUBLIC GEOS_DQsat

  PUBLIC GEOS_TRILU
  PUBLIC GEOS_TRISOLVE

  !EOP

  INTERFACE GEOS_QsatICE
     MODULE PROCEDURE QSATICE0
     MODULE PROCEDURE QSATICE1
     MODULE PROCEDURE QSATICE2
     MODULE PROCEDURE QSATICE3
  END INTERFACE

  INTERFACE GEOS_QsatLQU
     MODULE PROCEDURE QSATLQU0
     MODULE PROCEDURE QSATLQU1
     MODULE PROCEDURE QSATLQU2
     MODULE PROCEDURE QSATLQU3
  END INTERFACE

  INTERFACE GEOS_DQsat
     MODULE PROCEDURE DQSAT0
     MODULE PROCEDURE DQSAT1
     MODULE PROCEDURE DQSAT2
     MODULE PROCEDURE DQSAT3
  END INTERFACE

  INTERFACE GEOS_Qsat
     MODULE PROCEDURE QSAT0
     MODULE PROCEDURE QSAT1
     MODULE PROCEDURE QSAT2
     MODULE PROCEDURE QSAT3
  END INTERFACE

  INTERFACE GEOS_TRILU
     MODULE PROCEDURE GEOS_TRILU1
     MODULE PROCEDURE GEOS_TRILU2
     MODULE PROCEDURE GEOS_TRILU3
  END INTERFACE

  INTERFACE GEOS_TRISOLVE
     MODULE PROCEDURE GEOS_TRISOLVE1
     MODULE PROCEDURE GEOS_TRISOLVE2
     MODULE PROCEDURE GEOS_TRISOLVE3
  END INTERFACE
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_H2OMW  = 18.01_r8                  ! kg/Kmole
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_AIRMW  = 28.97_r8                  ! kg/Kmole
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_TICE   = 273.16_r8                 ! K


  REAL(KIND=r8),    PARAMETER :: ESFAC = MAPL_H2OMW/MAPL_AIRMW
  REAL(KIND=r8),    PARAMETER :: MAX_MIXING_RATIO = 1.0_r8  
  REAL(KIND=r8),    PARAMETER :: ZEROC   = MAPL_TICE

  REAL(KIND=r8),    PARAMETER :: TMINTBL    =  150.0_r8
  REAL(KIND=r8),    PARAMETER :: TMAXTBL    =  333.0_r8
  INTEGER, PARAMETER :: DEGSUBS    =  100
  REAL(KIND=r8),    PARAMETER :: ERFAC      = (DEGSUBS/ESFAC)
  REAL(KIND=r8),    PARAMETER :: DELTA_T    =  1.0_r8 / DEGSUBS
  INTEGER, PARAMETER :: TABLESIZE  =  (TMAXTBL-TMINTBL)*DEGSUBS + 1
  REAL(KIND=r8),    PARAMETER :: TMIX       = -20.0_r8

  LOGICAL      :: UTBL       = .TRUE.
  INTEGER      :: TYPE       =  1

  LOGICAL      :: FIRST      = .TRUE.

  REAL(KIND=r8)      :: ESTFRZ
  REAL(KIND=r8)      :: ESTLQU

  REAL(KIND=r8)      :: ESTBLE(TABLESIZE)
  REAL(KIND=r8)      :: ESTBLW(TABLESIZE)
  REAL(KIND=r8)      :: ESTBLX(TABLESIZE)

  REAL(KIND=r8),    PARAMETER :: TMINSTR = -95.0_r8
  REAL(KIND=r8),    PARAMETER :: TSTARR1 = -75.0_r8
  REAL(KIND=r8),    PARAMETER :: TSTARR2 = -65.0_r8
  REAL(KIND=r8),    PARAMETER :: TSTARR3 = -50.0_r8
  REAL(KIND=r8),    PARAMETER :: TSTARR4 = -40.0_r8
  REAL(KIND=r8),    PARAMETER :: TMAXSTR = +60.0_r8

  REAL(KIND=r8),  PARAMETER :: B6 = 6.136820929E-11_r8*100.0_r8
  REAL(KIND=r8),  PARAMETER :: B5 = 2.034080948E-8_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: B4 = 3.031240396E-6_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: B3 = 2.650648471E-4_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: B2 = 1.428945805E-2_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: B1 = 4.436518521E-1_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: B0 = 6.107799961E+0_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI6= 1.838826904E-10_r8*100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI5= 4.838803174E-8_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI4= 5.824720280E-6_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI3= 4.176223716E-4_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI2= 1.886013408E-2_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI1= 5.034698970E-1_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: BI0= 6.109177956E+0_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S16= 0.516000335E-11_r8*100.0_r8
  REAL(KIND=r8),  PARAMETER :: S15= 0.276961083E-8_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S14= 0.623439266E-6_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S13= 0.754129933E-4_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S12= 0.517609116E-2_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S11= 0.191372282E+0_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S10= 0.298152339E+1_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S26= 0.314296723E-10_r8*100.0_r8
  REAL(KIND=r8),  PARAMETER :: S25= 0.132243858E-7_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S24= 0.236279781E-5_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S23= 0.230325039E-3_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S22= 0.129690326E-1_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S21= 0.401390832E+0_r8 *100.0_r8
  REAL(KIND=r8),  PARAMETER :: S20= 0.535098336E+1_r8 *100.0_r8


  REAL(KIND=r8), PARAMETER  :: DI(0:3)=(/ 57518.5606E08_r8, 2.01889049_r8, 3.56654_r8, 20.947031_r8 /)
  REAL(KIND=r8), PARAMETER  :: CI(0:3)=(/ 9.550426_r8, -5723.265_r8, 3.53068_r8, -.00728332_r8 /)
  REAL(KIND=r8), PARAMETER  :: DL(6)=(/  -7.902980_r8, 5.02808_r8, -1.3816_r8, 11.344_r8, 8.1328_r8, -3.49149_r8 /)
  REAL(KIND=r8), PARAMETER  :: LOGPS = 3.005714898_r8  ! log10(1013.246)
  REAL(KIND=r8), PARAMETER  :: TS = 373.16_r8
  REAL(KIND=r8), PARAMETER  :: CL(0:9)=(/54.842763_r8, -6763.22_r8, -4.21000_r8, .000367_r8, &
       .0415_r8, 218.8_r8,  53.878000_r8, -1331.22_r8, -9.44523_r8, .014025_r8  /)


  REAL(KIND=r8) :: TMINLQU    =  ZEROC - 40.0_r8
  REAL(KIND=r8) :: TMINICE    =  ZEROC + TMINSTR




CONTAINS

  !BOPI

  ! !IROUTINE: GEOS_QsatLqu Computes saturation specific humidity over
  !            liquid water.
  ! !IROUTINE: GEOS_QsatIce Computes saturation specific humidity over
  !            frozen water.

  ! !INTERFACE:

  !    function GEOS_QsatLqu(TL,PL,DQ) result(QS)
  !    function GEOS_QsatIce(TL,PL,DQ) result(QS)
  !
  ! Overloads:
  !
  !      real,                               intent(IN)               :: TL
  !      logical,                  optional, intent(IN)               :: PL
  !      real,                     optional, intent(OUT)              :: DQ
  !      real                                                         :: QS
  !
  !      real,                               intent(IN)               :: TL(:)
  !      logical,                  optional, intent(IN)               :: PL(:)
  !      real,                     optional, intent(OUT)              :: DQ(:)
  !      real, dimension(size(TL,1))                                  :: QS
  !
  !      REAL(KIND=r8),                               intent(IN)               :: TL(:,:)
  !      logical,                  optional, intent(IN)               :: PL(:,:)
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: DQ(:,:)
  !      REAL(KIND=r8), dimension(size(TL,1),size(TL,2))                       :: QS
  !
  !      REAL(KIND=r8),                               intent(IN)               :: TL(:,:,:)
  !      logical,                  optional, intent(IN)               :: PL(:,:,:)
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: DQ(:,:,:)
  !      REAL(KIND=r8), dimension(size(TL,1),size(TL,2),size(TL,3))            :: QS
  !
  !

  ! !DESCRIPTION:  Uses various formulations of the saturation
  !                vapor pressure to compute the saturation specific 
  !    humidity and, optionally, its derivative with respect to temperature
  !    for temperature TL and pressure PL. If PL is not present
  !    it returns the saturation vapor pressure and, optionally, its derivative. 
  !
  !    All pressures are in Pascals and all temperatures in Kelvins.
  !
  !    The choice of saturation vapor pressure formulation is controlled by  GEOS_QsatSet.
  !    Three choices are currently supported: The CAM formulation,
  !    Murphy and Koop (2005, QJRMS), and the Staar formulation from NSIPP-1.
  !    The default is Starr. All three are valid up to 333K. Above the 
  !    freezing point, GEOS_QsatIce returns values at the freezing point.
  !    Murphy and Koop is valid down to 150K, for both liquid and ice.
  !    The other two are valid down to 178K for ice and 233K for super-cooled liquid. 
  !
  !    Another choice is whether to use the exact formulation
  !    or a table look-up. This can also be controlled with GEOS_QsatSet.
  !    The default is to do a table look-up. The tables are generated
  !    at 0.1C intervals, controlled by parameter DEGSUBS=10.
  ! 
  !    
  !EOPI


  FUNCTION QSATLQU0(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ
    REAL(KIND=r8)    :: QS

    REAL(KIND=r8)    :: TI,W
    REAL(KIND=r8)    :: DD
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    INTEGER :: IT
    REAL(KIND=r8)    :: TX 
    REAL(KIND=r8)    :: PX 
    REAL(KIND=r8)    :: EX 
    REAL(KIND=r8)    :: DX 
    !#define TX TL
    TX=TL
    !#define PX PL
    IF(PRESENT(PL)) PX=PL
    !#define EX QS   
    !if(present(PL)) PX=PL
    !#define DX DQ
    !IF(PRESENT(DQ)) DX=DQ
    !#include "qsatlqu.code"

    IF(UTBL) THEN

       IF(FIRST) THEN
          FIRST = .FALSE.
          CALL ESINIT
       END IF

       IF    (TX<=TMINLQU) THEN
          EX=ESTLQU
          IF(PRESENT(DQ)) DDQ = 0.0_r8
       ELSEIF(TX>=TMAXTBL  ) THEN
          EX=ESTBLW(TABLESIZE)
          IF(PRESENT(DQ)) DDQ = 0.0_r8
       ELSE
          TT  = (TX - TMINTBL)*DEGSUBS+1
          IT  = INT(TT)
          DDQ = ESTBLW(IT+1) - ESTBLW(IT)
          EX  = ((TT-IT)*DDQ + ESTBLW(IT))
       END IF

       IF(PRESENT(PL)) THEN
          IF(PX > EX) THEN
             DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
             EX = EX*DD
             IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
          ELSE
             EX  = MAX_MIXING_RATIO
             IF(PRESENT(DQ)) DX = 0.0_r8
          END IF
       ELSE
          IF(PRESENT(DQ)) DX = DDQ
       END IF

    ELSE  ! Exact Formulation

       IF    (TX<TMINLQU) THEN
          TI = TMINLQU
       ELSEIF(TX>TMAXTBL) THEN
          TI = TMAXTBL
       ELSE
          TI = TX
       END IF

       !#include "esatlqu.code"
       IF    (TYPE==1) THEN
          TT = TI-ZEROC       !  Starr polynomial fit
          EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
       ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
          TT = TS/TI
          EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
               DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
               DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
               LOGPS + 2.0_r8                                               )
       ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
          EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
               TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
       ENDIF

       !#include "esatlqu.code"
       IF(PRESENT(DQ)) THEN
          IF    (TX<TMINLQU) THEN
             DDQ = 0.0_r8
          ELSEIF(TX>TMAXTBL) THEN
             DDQ = 0.0_r8
          ELSE
             IF(PX>EX) THEN
                DD = EX
                TI = TX + DELTA_T
                !#include "esatlqu.code"
                IF    (TYPE==1) THEN
                   TT = TI-ZEROC       !  Starr polynomial fit
                   EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
                ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
                   TT = TS/TI
                   EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                        DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                        DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                        LOGPS + 2.0_r8                                               )
                ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                   EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                        TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
                ENDIF

                !#include "esatlqu.code"
                DDQ = EX-DD
                EX  = DD
             ENDIF
          END IF
       END IF

       IF(PRESENT(PL)) THEN
          IF(PX > EX) THEN
             DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
             EX = EX*DD
             IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
          ELSE
             EX = MAX_MIXING_RATIO
             IF(PRESENT(DQ)) DX = 0.0_r8
          END IF
       ELSE
          IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
       END IF

    ENDIF ! not table
    QS   = EX 
    IF(PRESENT(DQ)) DQ=DX

    !#undef  DX
    !#undef  TX
    !#undef  EX
    !#undef  PX
    RETURN
  END FUNCTION QSATLQU0

  FUNCTION QSATLQU1(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL(:)
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL(:)
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ(:)
    REAL(KIND=r8)    :: QS(SIZE(TL,1))
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: DX
    INTEGER :: I
    REAL(KIND=r8)    :: TI,W  
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    REAL(KIND=r8)    :: DD
    INTEGER :: IT
    DO I=1,SIZE(TL,1)
       !#define TX TL
       TX=TL(I)
       !#define PX PL
       IF(PRESENT(PL)) PX=PL(I)
       !#define EX QS   
       !if(present(PL)) PX=PL
       !#define DX DQ
       !IF(PRESENT(DQ)) DX=DQ(I)
       !#include "qsatlqu.code"


       !#define TX TL(I)
       !#define PX PL(I)
       !#define EX QS(I)
       !#define DX DQ(I)
       !#include "qsatlqu.code"

       IF(UTBL) THEN

          IF(FIRST) THEN
             FIRST = .FALSE.
             CALL ESINIT
          END IF

          IF    (TX<=TMINLQU) THEN
             EX=ESTLQU
             IF(PRESENT(DQ)) DDQ = 0.0_r8
          ELSEIF(TX>=TMAXTBL  ) THEN
             EX=ESTBLW(TABLESIZE)
             IF(PRESENT(DQ)) DDQ = 0.0_r8
          ELSE
             TT  = (TX - TMINTBL)*DEGSUBS+1
             IT  = INT(TT)
             DDQ = ESTBLW(IT+1) - ESTBLW(IT)
             EX  = ((TT-IT)*DDQ + ESTBLW(IT))
          END IF

          IF(PRESENT(PL)) THEN
             IF(PX > EX) THEN
                DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
                EX = EX*DD
                IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
             ELSE
                EX  = MAX_MIXING_RATIO
                IF(PRESENT(DQ)) DX = 0.0_r8
             END IF
          ELSE
             IF(PRESENT(DQ)) DX = DDQ
          END IF

       ELSE  ! Exact Formulation

          IF    (TX<TMINLQU) THEN
             TI = TMINLQU
          ELSEIF(TX>TMAXTBL) THEN
             TI = TMAXTBL
          ELSE
             TI = TX
          END IF


          !#include "esatlqu.code"
          IF    (TYPE==1) THEN
             TT = TI-ZEROC       !  Starr polynomial fit
             EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
          ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
             TT = TS/TI
             EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                  DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                  DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                  LOGPS + 2.0_r8                                               )
          ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
             EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                  TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
          ENDIF

          !#include "esatlqu.code"

          IF(PRESENT(DQ)) THEN
             IF    (TX<TMINLQU) THEN
                DDQ = 0.0_r8
             ELSEIF(TX>TMAXTBL) THEN
                DDQ = 0.0_r8
             ELSE
                IF(PX>EX) THEN
                   DD = EX
                   TI = TX + DELTA_T
                   !#include "esatlqu.code"
                   IF    (TYPE==1) THEN
                      TT = TI-ZEROC       !  Starr polynomial fit
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
                   ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
                      TT = TS/TI
                      EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                           DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                           DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                           LOGPS + 2.0_r8                                               )
                   ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                      EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                           TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
                   ENDIF

                   !#include "esatlqu.code"
                   DDQ = EX-DD
                   EX  = DD
                ENDIF
             END IF
          END IF

          IF(PRESENT(PL)) THEN
             IF(PX > EX) THEN
                DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
                EX = EX*DD
                IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
             ELSE
                EX = MAX_MIXING_RATIO
                IF(PRESENT(DQ)) DX = 0.0_r8
             END IF
          ELSE
             IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
          END IF

       ENDIF ! not table

       !#undef  DX
       !#undef  TX
       !#undef  PX
       !#undef  EX
       QS(I)   = EX
       IF(PRESENT(DQ)) DQ(I)=DX

    END DO
  END FUNCTION QSATLQU1

  FUNCTION QSATLQU2(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL(:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ(:,:)
    REAL(KIND=r8)    :: QS(SIZE(TL,1),SIZE(TL,2))
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: DX

    INTEGER :: I, J
    REAL(KIND=r8)    :: TI,W  
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    REAL(KIND=r8)    :: DD
    INTEGER :: IT
    DO J=1,SIZE(TL,2)
       DO I=1,SIZE(TL,1)
          !#define TX TL
          TX=TL(I,J)
          !#define PX PL
          IF(PRESENT(PL)) PX=PL(I,J)
          !#define EX QS   
          !if(present(PL)) PX=PL
          !#define DX DQ
          !IF(PRESENT(DQ)) DX=DQ(I,J)

          !#define TX TL(I,J)
          !#define PX PL(I,J)
          !#define EX QS(I,J)
          !#define DX DQ(I,J)
          !#include "qsatlqu.code"

          IF(UTBL) THEN

             IF(FIRST) THEN
                FIRST = .FALSE.
                CALL ESINIT
             END IF

             IF    (TX<=TMINLQU) THEN
                EX=ESTLQU
                IF(PRESENT(DQ)) DDQ = 0.0_r8
             ELSEIF(TX>=TMAXTBL  ) THEN
                EX=ESTBLW(TABLESIZE)
                IF(PRESENT(DQ)) DDQ = 0.0_r8
             ELSE
                TT  = (TX - TMINTBL)*DEGSUBS+1
                IT  = INT(TT)
                DDQ = ESTBLW(IT+1) - ESTBLW(IT)
                EX  = ((TT-IT)*DDQ + ESTBLW(IT))
             END IF

             IF(PRESENT(PL)) THEN
                IF(PX > EX) THEN
                   DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
                   EX = EX*DD
                   IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                ELSE
                   EX  = MAX_MIXING_RATIO
                   IF(PRESENT(DQ)) DX = 0.0_r8
                END IF
             ELSE
                IF(PRESENT(DQ)) DX = DDQ
             END IF

          ELSE  ! Exact Formulation

             IF    (TX<TMINLQU) THEN
                TI = TMINLQU
             ELSEIF(TX>TMAXTBL) THEN
                TI = TMAXTBL
             ELSE
                TI = TX
             END IF
             !#include "esatlqu.code"
             IF    (TYPE==1) THEN
                TT = TI-ZEROC       !  Starr polynomial fit
                EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
             ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
                TT = TS/TI
                EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                     DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                     DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                     LOGPS + 2.0_r8                                               )
             ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                     TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
             ENDIF

             !#include "esatlqu.code"

             IF(PRESENT(DQ)) THEN
                IF    (TX<TMINLQU) THEN
                   DDQ = 0.0_r8
                ELSEIF(TX>TMAXTBL) THEN
                   DDQ = 0.0_r8
                ELSE
                   IF(PX>EX) THEN
                      DD = EX
                      TI = TX + DELTA_T
                      !#include "esatlqu.code"
                      IF    (TYPE==1) THEN
                         TT = TI-ZEROC       !  Starr polynomial fit
                         EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
                      ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
                         TT = TS/TI
                         EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                              DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                              DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                              LOGPS + 2.0_r8                                               )
                      ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                         EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                              TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
                      ENDIF
                      !#include "esatlqu.code"
                      DDQ = EX-DD
                      EX  = DD
                   ENDIF
                END IF
             END IF

             IF(PRESENT(PL)) THEN
                IF(PX > EX) THEN
                   DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
                   EX = EX*DD
                   IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                ELSE
                   EX = MAX_MIXING_RATIO
                   IF(PRESENT(DQ)) DX = 0.0_r8
                END IF
             ELSE
                IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
             END IF

          ENDIF ! not table

          !#include "qsatlqu.code"
          !#undef  DX
          !#undef  TX
          !#undef  PX
          !#undef  EX
          QS (I,J)  = EX
          IF(PRESENT(DQ)) DQ(I,J)=DX

       END DO
    END DO
  END FUNCTION QSATLQU2

  FUNCTION QSATLQU3(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL(:,:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ(:,:,:)
    REAL(KIND=r8)    :: QS(SIZE(TL,1),SIZE(TL,2),SIZE(TL,3))
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: DX
    INTEGER :: I, J, K
    REAL(KIND=r8)    :: TI,W  
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    REAL(KIND=r8)    :: DD
    INTEGER :: IT
    DO K=1,SIZE(TL,3)
       DO J=1,SIZE(TL,2)
          DO I=1,SIZE(TL,1)
             !#define TX TL
             TX=TL(I,J,K)
             !#define PX PL
             IF(PRESENT(PL)) PX=PL(I,J,K)
             !#define EX QS   
             !if(present(PL)) PX=PL
             !#define DX DQ
             !IF(PRESENT(DQ)) DX=DQ(I,J,K)
             !#define TX TL(I,J,K)
             !#define PX PL(I,J,K)
             !#define EX QS(I,J,K)
             !#define DX DQ(I,J,K)
             !#include "qsatlqu.code"

             IF(UTBL) THEN

                IF(FIRST) THEN
                   FIRST = .FALSE.
                   CALL ESINIT
                END IF

                IF    (TX<=TMINLQU) THEN
                   EX=ESTLQU
                   IF(PRESENT(DQ)) DDQ = 0.0_r8
                ELSEIF(TX>=TMAXTBL  ) THEN
                   EX=ESTBLW(TABLESIZE)
                   IF(PRESENT(DQ)) DDQ = 0.0_r8
                ELSE
                   TT  = (TX - TMINTBL)*DEGSUBS+1
                   IT  = INT(TT)
                   DDQ = ESTBLW(IT+1) - ESTBLW(IT)
                   EX  = ((TT-IT)*DDQ + ESTBLW(IT))
                END IF

                IF(PRESENT(PL)) THEN
                   IF(PX > EX) THEN
                      DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
                      EX = EX*DD
                      IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                   ELSE
                      EX  = MAX_MIXING_RATIO
                      IF(PRESENT(DQ)) DX = 0.0_r8
                   END IF
                ELSE
                   IF(PRESENT(DQ)) DX = DDQ
                END IF

             ELSE  ! Exact Formulation

                IF    (TX<TMINLQU) THEN
                   TI = TMINLQU
                ELSEIF(TX>TMAXTBL) THEN
                   TI = TMAXTBL
                ELSE
                   TI = TX
                END IF

                !#include "esatlqu.code"
                IF    (TYPE==1) THEN
                   TT = TI-ZEROC       !  Starr polynomial fit
                   EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
                ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
                   TT = TS/TI
                   EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                        DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                        DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                        LOGPS + 2.0_r8                                               )
                ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                   EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                        TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
                ENDIF

                !#include "esatlqu.code"
                IF(PRESENT(DQ)) THEN
                   IF    (TX<TMINLQU) THEN
                      DDQ = 0.0_r8
                   ELSEIF(TX>TMAXTBL) THEN
                      DDQ = 0.0_r8
                   ELSE
                      IF(PX>EX) THEN
                         DD = EX
                         TI = TX + DELTA_T
                         !#include "esatlqu.code"
                         IF    (TYPE==1) THEN
                            TT = TI-ZEROC       !  Starr polynomial fit
                            EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)
                         ELSEIF(TYPE==2) THEN   !  Fit used in CAM.
                            TT = TS/TI
                            EX = 10.0_r8**(  DL(1)*(TT - 1.0_r8) + DL(2)*LOG10(TT)             + &
                                 DL(3)*(10.0_r8**(DL(4)*(1.0_r8 - (1.0_r8/TT))) - 1.0_r8)/10000000.0_r8 + &
                                 DL(5)*(10.0_r8**(DL(6)*(TT -   1.0_r8    )) - 1.0_r8)/1000.0_r8     + &
                                 LOGPS + 2.0_r8                                               )
                         ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                            EX = EXP(                     (CL(0) + CL(1)/TI + CL(2)*log(TI) + CL(3)*TI) + &
                                 TANH(CL(4)*(TI-CL(5))) * (CL(6) + CL(7)/TI + CL(8)*log(TI) + CL(9)*TI)   )
                         ENDIF

                         !#include "esatlqu.code"
                         DDQ = EX-DD
                         EX  = DD
                      ENDIF
                   END IF
                END IF

                IF(PRESENT(PL)) THEN
                   IF(PX > EX) THEN
                      DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
                      EX = EX*DD
                      IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                   ELSE
                      EX = MAX_MIXING_RATIO
                      IF(PRESENT(DQ)) DX = 0.0_r8
                   END IF
                ELSE
                   IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
                END IF

             ENDIF ! not table

             !#include "qsatlqu.code"
             !#undef  DX
             !#undef  TX
             !#undef  PX
             !#undef  EX
             QS(I,J,K)=EX 
             IF(PRESENT(DQ)) DQ(I,J,K)=DX
          END DO
       END DO
    END DO
  END FUNCTION QSATLQU3



  FUNCTION QSATICE0(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ
    REAL(KIND=r8)    :: QS
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: DX
    REAL(KIND=r8)    :: TI,W
    REAL(KIND=r8)    :: DD
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    INTEGER :: IT
    !#define TX TL
    TX=TL
    !#define PX PL
    IF(PRESENT(PL)) PX=PL
    !#define EX QS   
    !if(present(PL)) PX=PL
    !#define DX DQ

    !#define TX TL
    !#define PX PL
    !#define EX QS
    !#define DX DQ
    !#include "qsatice.code"

    IF(UTBL) THEN

       IF(FIRST) THEN
          FIRST = .FALSE.
          CALL ESINIT
       END IF

       IF    (TX<=TMINTBL) THEN
          EX=ESTBLE(1)
          IF(PRESENT(DQ)) DDQ = 0.0_r8
       ELSEIF(TX>=ZEROC  ) THEN
          EX=ESTFRZ
          IF(PRESENT(DQ)) DDQ = 0.0_r8
       ELSE
          TT  = (TX - TMINTBL)*DEGSUBS+1
          IT  = INT(TT)
          DDQ = ESTBLE(IT+1) - ESTBLE(IT)
          EX  =  ((TT-IT)*DDQ + ESTBLE(IT))
       END IF


       IF(PRESENT(PL)) THEN
          IF(PX > EX) THEN
             DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
             EX = EX*DD
             IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
          ELSE
             EX  = MAX_MIXING_RATIO
             IF(PRESENT(DQ)) DX = 0.0_r8
          END IF
       ELSE
          IF(PRESENT(DQ)) DX = DDQ
       END IF

    ELSE  ! Exact Formulation

       IF    (TX<TMINICE) THEN
          TI = TMINICE
       ELSEIF(TX>ZEROC  ) THEN
          TI = ZEROC
       ELSE
          TI = TX
       END IF

       !#include "esatice.code"

       IF( TYPE==1) THEN     ! Use Starr formulation
          TT = TI - ZEROC
          IF    (TT < TSTARR1                   ) THEN
             EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
          ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
             W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
             EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                  + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
          ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
             EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
          ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
             W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
             EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                  + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
          ELSE
             EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
          ENDIF
       ELSEIF(TYPE==2) THEN !  Fit used in CAM.
          TT = ZEROC/TI
          EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
       ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
          EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
       ENDIF

       !#include "esatice.code"
       IF(PRESENT(DQ)) THEN
          IF    (TX<TMINICE) THEN
             DDQ = 0.0_r8
          ELSEIF(TX>ZEROC  ) THEN
             DDQ = 0.0_r8
          ELSE
             IF(PX>EX) THEN
                DD = EX
                TI = TX + DELTA_T
                !#include "esatice.code"

                IF( TYPE==1) THEN     ! Use Starr formulation
                   TT = TI - ZEROC
                   IF    (TT < TSTARR1                   ) THEN
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
                   ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                      W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                      EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                           + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                   ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                   ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                      W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                      EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                           + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                   ELSE
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                   ENDIF
                ELSEIF(TYPE==2) THEN !  Fit used in CAM.
                   TT = ZEROC/TI
                   EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
                ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                   EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
                ENDIF

                !#include "esatice.code"
                DDQ = EX-DD
                EX  = DD
             ENDIF
          END IF
       END IF

       IF(PRESENT(PL)) THEN
          IF(PX > EX) THEN
             DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
             EX = EX*DD
             IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
          ELSE
             EX = MAX_MIXING_RATIO
             IF(PRESENT(DQ)) DX = 0.0_r8
          END IF
       ELSE
          IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
       END IF


    END IF




    !#include "qsatice.code"
    !#undef  DX
    !#undef  TX
    !#undef  EX
    !#undef  PX
    QS=EX 
    IF(PRESENT(DQ)) DQ=DX
    RETURN
  END FUNCTION QSATICE0

  FUNCTION QSATICE1(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL(:)
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL(:)
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ(:)
    REAL(KIND=r8)    :: QS(SIZE(TL,1))
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: DX
    INTEGER :: I
    REAL(KIND=r8)    :: TI,W  
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    REAL(KIND=r8)    :: DD
    INTEGER :: IT
    DO I=1,SIZE(TL,1)
       !#define TX TL
       TX=TL(I)
       !#define PX PL
       IF(PRESENT(PL)) PX=PL(I)
       !#define EX QS   
       !if(present(PL)) PX=PL
       !#define DX DQ
       !IF(PRESENT(DQ)) DX=DQ(I)

       !#define TX TL(I)
       !#define PX PL(I)
       !#define EX QS(I)
       !#define DX DQ(I)
       !#include "qsatice.code"

       IF(UTBL) THEN

          IF(FIRST) THEN
             FIRST = .FALSE.
             CALL ESINIT
          END IF

          IF    (TX<=TMINTBL) THEN
             EX=ESTBLE(1)
             IF(PRESENT(DQ)) DDQ = 0.0_r8
          ELSEIF(TX>=ZEROC  ) THEN
             EX=ESTFRZ
             IF(PRESENT(DQ)) DDQ = 0.0_r8
          ELSE
             TT  = (TX - TMINTBL)*DEGSUBS+1
             IT  = INT(TT)
             DDQ = ESTBLE(IT+1) - ESTBLE(IT)
             EX  =  ((TT-IT)*DDQ + ESTBLE(IT))
          END IF


          IF(PRESENT(PL)) THEN
             IF(PX > EX) THEN
                DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
                EX = EX*DD
                IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
             ELSE
                EX  = MAX_MIXING_RATIO
                IF(PRESENT(DQ)) DX = 0.0_r8
             END IF
          ELSE
             IF(PRESENT(DQ)) DX = DDQ
          END IF

       ELSE  ! Exact Formulation

          IF    (TX<TMINICE) THEN
             TI = TMINICE
          ELSEIF(TX>ZEROC  ) THEN
             TI = ZEROC
          ELSE
             TI = TX
          END IF

          !#include "esatice.code"

          IF( TYPE==1) THEN     ! Use Starr formulation
             TT = TI - ZEROC
             IF    (TT < TSTARR1                   ) THEN
                EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
             ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                     + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
             ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
             ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                     + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
             ELSE
                EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
             ENDIF
          ELSEIF(TYPE==2) THEN !  Fit used in CAM.
             TT = ZEROC/TI
             EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
          ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
             EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
          ENDIF

          !#include "esatice.code"
          IF(PRESENT(DQ)) THEN
             IF    (TX<TMINICE) THEN
                DDQ = 0.0_r8
             ELSEIF(TX>ZEROC  ) THEN
                DDQ = 0.0_r8
             ELSE
                IF(PX>EX) THEN
                   DD = EX
                   TI = TX + DELTA_T
                   !#include "esatice.code"

                   IF( TYPE==1) THEN     ! Use Starr formulation
                      TT = TI - ZEROC
                      IF    (TT < TSTARR1                   ) THEN
                         EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
                      ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                         W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                         EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                              + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                      ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                         EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                      ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                         W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                         EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                              + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                      ELSE
                         EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                      ENDIF
                   ELSEIF(TYPE==2) THEN !  Fit used in CAM.
                      TT = ZEROC/TI
                      EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
                   ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                      EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
                   ENDIF

                   !#include "esatice.code"
                   DDQ = EX-DD
                   EX  = DD
                ENDIF
             END IF
          END IF

          IF(PRESENT(PL)) THEN
             IF(PX > EX) THEN
                DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
                EX = EX*DD
                IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
             ELSE
                EX = MAX_MIXING_RATIO
                IF(PRESENT(DQ)) DX = 0.0_r8
             END IF
          ELSE
             IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
          END IF


       END IF




       !#include "qsatice.code"
       !#undef  DX
       !#undef  TX
       !#undef  PX
       !#undef  EX
       QS(I)=EX 
       IF(PRESENT(DQ)) DQ(I)=DX
    END DO
  END FUNCTION QSATICE1

  FUNCTION QSATICE2(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL(:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ(:,:)
    REAL(KIND=r8)    :: QS(SIZE(TL,1),SIZE(TL,2))
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: DX

    INTEGER :: I, J
    REAL(KIND=r8)    :: TI,W  
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    REAL(KIND=r8)    :: DD
    INTEGER :: IT
    DO J=1,SIZE(TL,2)
       DO I=1,SIZE(TL,1)
          !#define TX TL
          TX=TL(I,J)
          !#define PX PL
          IF(PRESENT(PL)) PX=PL(I,J)
          !#define EX QS   
          !if(present(PL)) PX=PL
          !#define DX DQ
          !IF(PRESENT(DQ)) DX=DQ(I,J)


          !#define TX TL(I,J)
          !#define PX PL(I,J)
          !#define EX QS(I,J)
          !#define DX DQ(I,J)
          !#include "qsatice.code"

          IF(UTBL) THEN

             IF(FIRST) THEN
                FIRST = .FALSE.
                CALL ESINIT
             END IF

             IF    (TX<=TMINTBL) THEN
                EX=ESTBLE(1)
                IF(PRESENT(DQ)) DDQ = 0.0_r8
             ELSEIF(TX>=ZEROC  ) THEN
                EX=ESTFRZ
                IF(PRESENT(DQ)) DDQ = 0.0_r8
             ELSE
                TT  = (TX - TMINTBL)*DEGSUBS+1
                IT  = INT(TT)
                DDQ = ESTBLE(IT+1) - ESTBLE(IT)
                EX  =  ((TT-IT)*DDQ + ESTBLE(IT))
             END IF


             IF(PRESENT(PL)) THEN
                IF(PX > EX) THEN
                   DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
                   EX = EX*DD
                   IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                ELSE
                   EX  = MAX_MIXING_RATIO
                   IF(PRESENT(DQ)) DX = 0.0_r8
                END IF
             ELSE
                IF(PRESENT(DQ)) DX = DDQ
             END IF

          ELSE  ! Exact Formulation

             IF    (TX<TMINICE) THEN
                TI = TMINICE
             ELSEIF(TX>ZEROC  ) THEN
                TI = ZEROC
             ELSE
                TI = TX
             END IF


             !#include "esatice.code"

             IF( TYPE==1) THEN     ! Use Starr formulation
                TT = TI - ZEROC
                IF    (TT < TSTARR1                   ) THEN
                   EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
                ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                   W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                   EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                        + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                   EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                   W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                   EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                        + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                ELSE
                   EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                ENDIF
             ELSEIF(TYPE==2) THEN !  Fit used in CAM.
                TT = ZEROC/TI
                EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
             ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
             ENDIF

             !#include "esatice.code"

             IF(PRESENT(DQ)) THEN
                IF    (TX<TMINICE) THEN
                   DDQ = 0.0_r8
                ELSEIF(TX>ZEROC  ) THEN
                   DDQ = 0.0_r8
                ELSE
                   IF(PX>EX) THEN
                      DD = EX
                      TI = TX + DELTA_T

                      !#include "esatice.code"

                      IF( TYPE==1) THEN     ! Use Starr formulation
                         TT = TI - ZEROC
                         IF    (TT < TSTARR1                   ) THEN
                            EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
                         ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                            W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                            EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                                 + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                         ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                            EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                         ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                            W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                            EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                                 + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                         ELSE
                            EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                         ENDIF
                      ELSEIF(TYPE==2) THEN !  Fit used in CAM.
                         TT = ZEROC/TI
                         EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
                      ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                         EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
                      ENDIF

                      !#include "esatice.code"
                      DDQ = EX-DD
                      EX  = DD
                   ENDIF
                END IF
             END IF

             IF(PRESENT(PL)) THEN
                IF(PX > EX) THEN
                   DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
                   EX = EX*DD
                   IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                ELSE
                   EX = MAX_MIXING_RATIO
                   IF(PRESENT(DQ)) DX = 0.0_r8
                END IF
             ELSE
                IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
             END IF


          END IF

          !#include "qsatice.code"
          !#undef  DX
          !#undef  TX
          !#undef  PX
          !#undef  EX 
          QS(I,J)= EX   
          IF(PRESENT(DQ)) DQ(I,J)=DX
       END DO
    END DO
  END FUNCTION QSATICE2

  FUNCTION QSATICE3(TL,PL,DQ) RESULT(QS)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(IN) :: PL(:,:,:)
    REAL(KIND=r8), OPTIONAL,    INTENT(OUT):: DQ(:,:,:)
    REAL(KIND=r8)    :: QS(SIZE(TL,1),SIZE(TL,2),SIZE(TL,3))
    REAL(KIND=r8)    :: TX
    REAL(KIND=r8)    :: PX
    REAL(KIND=r8)    :: EX
    REAL(KIND=r8)    :: DX
    INTEGER :: I, J, K
    REAL(KIND=r8)    :: TI,W  
    REAL(KIND=r8)    :: TT
    REAL(KIND=r8)    :: DDQ
    REAL(KIND=r8)    :: DD
    INTEGER :: IT
    DO K=1,SIZE(TL,3)
       DO J=1,SIZE(TL,2)
          DO I=1,SIZE(TL,1)

             !#define TX TL
             TX=TL(I,J,K)
             !#define PX PL
             IF(PRESENT(PL)) PX=PL(I,J,K)
             !#define EX QS   
             !if(present(PL)) PX=PL
             !#define DX DQ
             !IF(PRESENT(DQ)) DX=DQ(I,J,K)

             !#define TX TL(I,J,K)
             !#define PX PL(I,J,K)
             !#define EX QS(I,J,K)
             !#define DX DQ(I,J,K)
             !#include "qsatice.code"

             IF(UTBL) THEN

                IF(FIRST) THEN
                   FIRST = .FALSE.
                   CALL ESINIT
                END IF

                IF    (TX<=TMINTBL) THEN
                   EX=ESTBLE(1)
                   IF(PRESENT(DQ)) DDQ = 0.0_r8
                ELSEIF(TX>=ZEROC  ) THEN
                   EX=ESTFRZ
                   IF(PRESENT(DQ)) DDQ = 0.0_r8
                ELSE
                   TT  = (TX - TMINTBL)*DEGSUBS+1
                   IT  = INT(TT)
                   DDQ = ESTBLE(IT+1) - ESTBLE(IT)
                   EX  =  ((TT-IT)*DDQ + ESTBLE(IT))
                END IF


                IF(PRESENT(PL)) THEN
                   IF(PX > EX) THEN
                      DD = (ESFAC/(PX - (1.0_r8-ESFAC)*EX))
                      EX = EX*DD
                      IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                   ELSE
                      EX  = MAX_MIXING_RATIO
                      IF(PRESENT(DQ)) DX = 0.0_r8
                   END IF
                ELSE
                   IF(PRESENT(DQ)) DX = DDQ
                END IF

             ELSE  ! Exact Formulation

                IF    (TX<TMINICE) THEN
                   TI = TMINICE
                ELSEIF(TX>ZEROC  ) THEN
                   TI = ZEROC
                ELSE
                   TI = TX
                END IF

                !#include "esatice.code"

                IF( TYPE==1) THEN     ! Use Starr formulation
                   TT = TI - ZEROC
                   IF    (TT < TSTARR1                   ) THEN
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
                   ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                      W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                      EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                           + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                   ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                   ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                      W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                      EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                           + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                   ELSE
                      EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                   ENDIF
                ELSEIF(TYPE==2) THEN !  Fit used in CAM.
                   TT = ZEROC/TI
                   EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
                ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                   EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
                ENDIF

                !#include "esatice.code"

                IF(PRESENT(DQ)) THEN
                   IF    (TX<TMINICE) THEN
                      DDQ = 0.0_r8
                   ELSEIF(TX>ZEROC  ) THEN
                      DDQ = 0.0_r8
                   ELSE
                      IF(PX>EX) THEN
                         DD = EX
                         TI = TX + DELTA_T
                         !#include "esatice.code"

                         IF( TYPE==1) THEN     ! Use Starr formulation
                            TT = TI - ZEROC
                            IF    (TT < TSTARR1                   ) THEN
                               EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
                            ELSEIF(TT >= TSTARR1 .AND. TT < TSTARR2) THEN
                               W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
                               EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
                                    + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                            ELSEIF(TT >= TSTARR2 .AND. TT < TSTARR3) THEN
                               EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
                            ELSEIF(TT >= TSTARR3 .AND. TT < TSTARR4) THEN
                               W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
                               EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
                                    + (1.0_r8-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                            ELSE
                               EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
                            ENDIF
                         ELSEIF(TYPE==2) THEN !  Fit used in CAM.
                            TT = ZEROC/TI
                            EX = DI(0) * EXP( -(DI(1)/TT + DI(2)*LOG(TT) + DI(3)*TT) )
                         ELSEIF(TYPE==3) THEN   !  Murphy and Koop (2005, QJRMS)
                            EX = EXP( CI(0)+ CI(1)/TI + CI(2)*log(TI) + CI(3)*TI )
                         ENDIF

                         !#include "esatice.code"
                         DDQ = EX-DD
                         EX  = DD
                      ENDIF
                   END IF
                END IF

                IF(PRESENT(PL)) THEN
                   IF(PX > EX) THEN
                      DD = ESFAC/(PX - (1.0_r8-ESFAC)*EX)
                      EX = EX*DD
                      IF(PRESENT(DQ)) DX = DDQ*ERFAC*PX*DD*DD
                   ELSE
                      EX = MAX_MIXING_RATIO
                      IF(PRESENT(DQ)) DX = 0.0_r8
                   END IF
                ELSE
                   IF(PRESENT(DQ)) DX = DDQ*(1.0_r8/DELTA_T)
                END IF


             END IF




             !#include "qsatice.code"
             !#undef  DX
             !#undef  TX
             !#undef  PX
             !#undef  EX
             QS(I,J,K)=EX 
             IF(PRESENT(DQ)) DQ(I,J,K)=DX
          END DO
       END DO
    END DO
  END FUNCTION QSATICE3



  !==============================================
  !==============================================

  !  Traditional Qsat and Dqsat (these are deprecated)

  !==============================================
  !==============================================

  !BOPI

  ! !IROUTINE: GEOS_Qsat -- Computes satuation specific humidity.

  ! !INTERFACE:

  !    function GEOS_Qsat(TL,PL,RAMP,PASCALS,DQSAT) result(QSAT)
  !
  ! Overloads:
  !
  !      REAL(KIND=r8),                      intent(IN)                        :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: DQSAT
  !      REAL(KIND=r8)                                                         :: QSAT
  !
  !      REAL(KIND=r8), dimension(:),        intent(IN)                        :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: DQSAT(:)
  !      REAL(KIND=r8), dimension(size(PL,1))                                  :: QSAT
  !
  !      REAL(KIND=r8), dimension(:,:),      intent(IN)                        :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: DQSAT(:,:)
  !      REAL(KIND=r8), dimension(size(PL,1),size(PL,2))                       :: QSAT
  !
  !      REAL(KIND=r8), dimension(:,:,:),    intent(IN)                        :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: DQSAT(:,:,:)
  !      REAL(KIND=r8), dimension(size(PL,1),size(PL,2),size(PL,3))            :: QSAT
  !

  ! !DESCRIPTION:  Uses various formulations of the saturation
  !                vapor pressure to compute the saturation specific 
  !    humidity for temperature TL and pressure PL.
  !
  !    For temperatures <= TMIX (-20C)
  !    the calculation is done over ice; for temperatures >= ZEROC (0C) the calculation
  !    is done over liquid water; and in between these values,
  !    it interpolates linearly between the two.
  !
  !    The optional RAMP is the width of this
  !    ice/water ramp (i.e., TMIX = ZEROC-RAMP); its default is 20.
  !
  !    If PASCALS is true, PL is
  !    assumed to be in Pa; if false or not present, it is assumed to be in mb.
  !
  !    The choice of saturation vapor pressure formulation is a compile-time
  !    option. Three choices are currently supported: The CAM formulation,
  !    Murphy and Koop (2005, QJRMS), and Staars formulation from NSIPP-1.
  !
  !    Another compile time choice is whether to use the exact formulation
  !    or a table look-up.
  !    If UTBL is true, tabled values of the saturation vapor pressures
  !    are used. These tables are automatically generated at a 0.1K resolution
  !    for whatever vapor pressure formulation is being used.
  ! 
  !    
  !EOPI


  FUNCTION QSAT0(TL,PL,RAMP,PASCALS,DQSAT) RESULT(QSAT)
    REAL(KIND=r8),   INTENT(IN) :: TL, PL
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: DQSAT
    REAL(KIND=r8)    :: QSAT

    REAL(KIND=r8)    :: URAMP, DD, QQ, TI, DQ, PP
    INTEGER :: IT

    IF(PRESENT(RAMP)) THEN
       URAMP = -ABS(RAMP)
    ELSE
       URAMP = TMIX
    END IF

    IF(PRESENT(PASCALS)) THEN
       IF(PASCALS) THEN
          PP = PL
       ELSE
          PP = PL*100.0_r8
       END IF
    ELSE
       PP = PL*100.0_r8
    END IF

    IF(URAMP==TMIX .OR. URAMP==0.0_r8 .AND. UTBL) THEN

       IF(FIRST) THEN
          FIRST = .FALSE.
          CALL ESINIT
       END IF

       IF    (TL<=TMINTBL) THEN
          TI = TMINTBL
       ELSEIF(TL>=TMAXTBL-.001_r8) THEN
          TI = TMAXTBL-.001_r8
       ELSE
          TI = TL
       END IF

       TI = (TI - TMINTBL)*DEGSUBS+1
       IT = INT(TI)

       IF(URAMP==TMIX) THEN
          DQ    = ESTBLX(IT+1) - ESTBLX(IT)
          QSAT  = (TI-IT)*DQ + ESTBLX(IT)
       ELSE
          DQ    = ESTBLE(IT+1) - ESTBLE(IT)
          QSAT  = (TI-IT)*DQ + ESTBLE(IT)
       ENDIF

       IF(PRESENT(DQSAT)) DQSAT = DQ*DEGSUBS

       IF(PP <= QSAT) THEN
          QSAT = MAX_MIXING_RATIO
          IF(PRESENT(DQSAT)) DQSAT = 0.0_r8
       ELSE
          DD = 1.0_r8/(PP - (1.0_r8-ESFAC)*QSAT)
          QSAT = ESFAC*QSAT*DD
          IF(PRESENT(DQSAT)) DQSAT = ESFAC*DQSAT*PP*(DD*DD)
       END IF

    ELSE

       TI = TL - ZEROC

       IF    (TI <= URAMP) THEN
          QSAT  =  QSATICE0(TL,PP,DQ=DQSAT)
       ELSEIF(TI >= 0.0_r8  ) THEN
          QSAT  =  QSATLQU0(TL,PP,DQ=DQSAT)
       ELSE
          QSAT  =  QSATICE0(TL,PP,DQ=DQSAT)
          QQ    =  QSATLQU0(TL,PP,DQ=DQ   )
          TI    =  TI/URAMP
          QSAT  =  TI*(QSAT - QQ) +  QQ
          IF(PRESENT(DQSAT)) DQSAT = TI*(DQSAT-DQ) + DQ
       END IF

    END IF

  END FUNCTION QSAT0

  FUNCTION QSAT1(TL,PL,RAMP,PASCALS,DQSAT) RESULT(QSAT)
    REAL(KIND=r8),              INTENT(IN) :: TL(:), PL(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: DQSAT(:)
    REAL(KIND=r8) :: QSAT(SIZE(TL,1))
    INTEGER :: I
    DO I=1,SIZE(TL,1)
       IF (PRESENT(DQSAT)) THEN
          QSAT(I) = QSAT0(TL(I),PL(I),RAMP,PASCALS,DQSAT(I))
       ELSE
          QSAT(I) = QSAT0(TL(I),PL(I),RAMP,PASCALS)
       END IF
    END DO
  END FUNCTION QSAT1

  FUNCTION QSAT2(TL,PL,RAMP,PASCALS,DQSAT) RESULT(QSAT)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:), PL(:,:)
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: DQSAT(:,:)
    REAL(KIND=r8) :: QSAT(SIZE(TL,1),SIZE(TL,2))
    INTEGER :: I, J
    DO J=1,SIZE(TL,2)
       DO I=1,SIZE(TL,1)
          IF (PRESENT(DQSAT)) THEN
             QSAT(I,J) = QSAT0(TL(I,J),PL(I,J),RAMP,PASCALS,DQSAT(I,J))
          ELSE
             QSAT(I,J) = QSAT0(TL(I,J),PL(I,J),RAMP,PASCALS)
          END IF
       END DO
    END DO
  END FUNCTION QSAT2

  FUNCTION QSAT3(TL,PL,RAMP,PASCALS,DQSAT) RESULT(QSAT)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:,:), PL(:,:,:)
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: DQSAT(:,:,:)
    REAL(KIND=r8) :: QSAT(SIZE(TL,1),SIZE(TL,2),SIZE(TL,3))
    INTEGER :: I, J, K
    DO K=1,SIZE(TL,3)
       DO J=1,SIZE(TL,2)
          DO I=1,SIZE(TL,1)
             IF (PRESENT(DQSAT)) THEN
                QSAT(I,J,K) = QSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS,DQSAT(I,J,K))
             ELSE
                QSAT(I,J,K) = QSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS)
             END IF
          END DO
       END DO
    END DO
  END FUNCTION QSAT3

  !=======================================================================================

  !BOPI

  ! !IROUTINE: GEOS_DQsat -- Computes derivative satuation specific humidity wrt temperature.

  ! !INTERFACE:

  !    function GEOS_DQsat(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
  !
  ! Overloads:
  !
  !      REAL(KIND=r8),                               intent(IN)               :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8),                     optional, intent(OUT)              :: QSAT
  !      REAL(KIND=r8)                                                         :: DQSAT
  !
  !      REAL(KIND=r8), dimension(:),                 intent(IN)               :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8), dimension(:),       optional, intent(OUT)              :: QSAT
  !      REAL(KIND=r8), dimension(size(PL,1))                                  :: DQSAT
  !
  !      REAL(KIND=r8), dimension(:,:),               intent(IN)               :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8), dimension(:,:),     optional, intent(OUT)              :: QSAT
  !      REAL(KIND=r8), dimension(size(PL,1),size(PL,2))                       :: DQSAT
  !
  !      REAL(KIND=r8), dimension(:,:,:),             intent(IN)               :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8), dimension(:,:,:),   optional, intent(OUT)              :: QSAT
  !      REAL(KIND=r8), dimension(size(PL,1),size(PL,2),size(PL,3))            :: DQSAT
  !
  !      REAL(KIND=r8), dimension(:,:,:,:),           intent(IN)               :: TL, PL
  !      logical,                  optional, intent(IN)               :: PASCALS
  !      REAL(KIND=r8),                     optional, intent(IN)               :: RAMP
  !      REAL(KIND=r8), dimension(:,:,:,:), optional, intent(OUT)              :: QSAT
  !      REAL(KIND=r8), dimension(size(PL,1),size(PL,2),size(PL,3),size(PL,4)) :: DQSAT

  ! !DESCRIPTION:  Differentiates the approximations used
  !                by GEOS_Qsat with respect to temperature,
  !    using the same scheme to handle ice. Arguments are as in 
  !    GEOS_Qsat, with the addition of QSAT, which is the saturation specific
  !    humidity. This is for economy, in case both qsat and dqsat are 
  !    required.
  !                

  !EOPI


  FUNCTION DQSAT0(TL,PL,RAMP,PASCALS,QSAT) RESULT(DQSAT)
    REAL(KIND=r8),   INTENT(IN) :: TL, PL
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: QSAT
    REAL(KIND=r8)    :: DQSAT
    REAL(KIND=r8)    :: URAMP, TT, WW, DD, DQQ, QQ, TI, DQI, QI, PP
    INTEGER :: IT

    IF(PRESENT(RAMP)) THEN
       URAMP = -ABS(RAMP)
    ELSE
       URAMP = TMIX
    END IF

    IF(PRESENT(PASCALS)) THEN
       IF(PASCALS) THEN
          PP = PL
       ELSE
          PP = PL*100.0_r8
       END IF
    ELSE
       PP = PL*100.0_r8
    END IF

    IF(URAMP==TMIX .OR. URAMP==0.0_r8 .AND. UTBL) THEN

       IF(FIRST) THEN
          FIRST = .FALSE.
          CALL ESINIT
       END IF

       IF    (TL<=TMINTBL) THEN
          TI = TMINTBL
       ELSEIF(TL>=TMAXTBL-.001_r8) THEN
          TI = TMAXTBL-.001_r8
       ELSE
          TI = TL
       END IF

       TT = (TI - TMINTBL)*DEGSUBS+1
       IT = INT(TT)

       IF(URAMP==TMIX) THEN
          DQQ =  ESTBLX(IT+1) - ESTBLX(IT)
          QQ  =  (TT-IT)*DQQ + ESTBLX(IT)
       ELSE
          DQQ =  ESTBLE(IT+1) - ESTBLE(IT)
          QQ  =  (TT-IT)*DQQ + ESTBLE(IT)
       ENDIF

       IF(PP <= QQ) THEN
          IF(PRESENT(QSAT)) QSAT = MAX_MIXING_RATIO
          DQSAT = 0.0_r8
       ELSE
          DD = 1.0_r8/(PP - (1.0_r8-ESFAC)*QQ)
          IF(PRESENT(QSAT)) QSAT = ESFAC*QQ*DD
          DQSAT = (ESFAC*DEGSUBS)*DQQ*PP*(DD*DD)
       END IF

    ELSE

       TI = TL - ZEROC

       IF    (TI <= URAMP) THEN
          QQ  = QSATICE0(TL,PP,DQ=DQSAT)
          IF(PRESENT(QSAT)) QSAT  = QQ
       ELSEIF(TI >= 0.0_r8  ) THEN
          QQ  = QSATLQU0(TL,PP,DQ=DQSAT)
          IF(PRESENT(QSAT)) QSAT  = QQ
       ELSE
          QQ  = QSATLQU0(TL,PP,DQ=DQQ)
          QI  = QSATICE0(TL,PP,DQ=DQI)
          TI  = TI/URAMP
          DQSAT = TI*(DQI - DQQ) + DQQ
          IF(PRESENT(QSAT)) QSAT  = TI*(QI - QQ) +  QQ
       END IF

    END IF

  END FUNCTION DQSAT0

  FUNCTION DQSAT1(TL,PL,RAMP,PASCALS,QSAT) RESULT(DQSAT)
    REAL(KIND=r8),              INTENT(IN) :: TL(:), PL(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: QSAT(:)
    REAL(KIND=r8) :: DQSAT(SIZE(TL,1))
    INTEGER :: I
    DO I=1,SIZE(TL,1)
       IF (PRESENT(QSAT)) THEN
          DQSAT(I) = DQSAT0(TL(I),PL(I),RAMP,PASCALS,QSAT(I))
       ELSE
          DQSAT(I) = DQSAT0(TL(I),PL(I),RAMP,PASCALS)
       ENDIF
    END DO
  END FUNCTION DQSAT1

  FUNCTION DQSAT2(TL,PL,RAMP,PASCALS,QSAT) RESULT(DQSAT)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:), PL(:,:)
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: QSAT(:,:)
    REAL(KIND=r8) :: DQSAT(SIZE(TL,1),SIZE(TL,2))
    INTEGER :: I, J
    DO J=1,SIZE(TL,2)
       DO I=1,SIZE(TL,1)
          IF (PRESENT(QSAT)) THEN
             DQSAT(I,J) = DQSAT0(TL(I,J),PL(I,J),RAMP,PASCALS,QSAT(I,J))
          ELSE
             DQSAT(I,J) = DQSAT0(TL(I,J),PL(I,J),RAMP,PASCALS)
          END IF
       END DO
    END DO
  END FUNCTION DQSAT2

  FUNCTION DQSAT3(TL,PL,RAMP,PASCALS,QSAT) RESULT(DQSAT)
    REAL(KIND=r8),              INTENT(IN) :: TL(:,:,:), PL(:,:,:)
    LOGICAL, OPTIONAL, INTENT(IN) :: PASCALS
    REAL(KIND=r8),    OPTIONAL, INTENT(IN) :: RAMP
    REAL(KIND=r8),    OPTIONAL, INTENT(OUT):: QSAT(:,:,:)
    REAL(KIND=r8) :: DQSAT(SIZE(TL,1),SIZE(TL,2),SIZE(TL,3))
    INTEGER :: I, J, K
    DO K=1,SIZE(TL,3)
       DO J=1,SIZE(TL,2)
          DO I=1,SIZE(TL,1)
             IF (PRESENT(QSAT)) THEN
                DQSAT(I,J,K) = DQSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS,QSAT(I,J,K))
             ELSE
                DQSAT(I,J,K) = DQSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS)
             END IF
          END DO
       END DO
    END DO
  END FUNCTION DQSAT3

  !==============================================

  !BOPI

  ! !IROUTINE: GEOS_QsatSet -- Sets behavior of GEOS_QsatLqu an GEOS_QsatIce

  ! !INTERFACE:

  SUBROUTINE GEOS_QsatSet(USETABLE,FORMULATION)
    LOGICAL, OPTIONAL, INTENT(IN) :: USETABLE
    INTEGER, OPTIONAL, INTENT(IN) :: FORMULATION

    ! !DESCRIPTION: GEOS_QsatSet can be used to modify 
    !  the behavior of GEOS_QsatLqu an GEOS_QsatIce 
    !  from its default setting.

    !  If {\tt \bf USETABLE} is true, tabled values of the saturation vapor pressures are used.
    !  These tables are automatically generated at a 0.1K resolution for whatever
    !  vapor pressure formulation is being used. The default is to use the table.

    !  {\tt \bf FORMULATION} sets the saturation vapor pressure function.
    !  Three formulations of saturation vapor pressure are supported: 
    !  the Starr code that was in NSIPP-1 (FORMULATION==1), the formulation in  CAM 
    !  (FORMULATION==2), and Murphy and Koop (2005, QJRMS) (FORMULATION==3).
    !  The default is FORMULATION=1.

    !  If appropriate, GEOS_QsatSet also initializes the tables. If GEOS_QsatSet is
    !  not called and tables are required, they will be initialized the first time
    !  a Qsat function is called.

    !EOPI

    IF(PRESENT(UseTable   )) UTBL = UseTable
    IF(PRESENT(Formulation)) TYPE = MAX(MIN(Formulation,3),1)

    IF(TYPE==3)  THEN ! Murphy and Koop (2005, QJRMS)
       TMINICE    =  MAX(TMINTBL,110.)
       TMINLQU    =  MAX(TMINTBL,123.)
    ELSE
       TMINLQU    =  ZEROC - 40.0
       TMINICE    =  ZEROC + TMINSTR
    ENDIF

    IF(UTBL) CALL ESINIT

    RETURN
  END SUBROUTINE GEOS_QsatSet

  !=======================================================================================

  SUBROUTINE ESINIT

    ! Saturation vapor pressure table initialization. This is invoked if UTBL is true 
    ! on the first call to any qsat routine or whenever GEOS_QsatSet is called 
    ! N.B.--Tables are in Pa

    INTEGER :: I
    REAL(KIND=r8)    :: T
    LOGICAL :: UT

    UT = UTBL
    UTBL=.FALSE.

    DO I=1,TABLESIZE

       T = (I-1)*DELTA_T + TMINTBL

       ESTBLW(I) = QSATLQU0(T)

       IF(T>ZEROC) THEN
          ESTBLE(I) = ESTBLW(I)
       ELSE
          ESTBLE(I) = QSATICE0(T)
       END IF

       T = T-ZEROC

       IF(T>=TMIX .AND. T<0.0_r8) THEN
          ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
       ELSE
          ESTBLX(I) = ESTBLE(I)
       END IF

    END DO

    ESTFRZ = QSATLQU0(ZEROC  )
    ESTLQU = QSATLQU0(TMINLQU)

    UTBL = UT

  END SUBROUTINE ESINIT












  !*************************************************************************
  !*************************************************************************

  !  Tridiagonal solvers

  !*************************************************************************
  !*************************************************************************

  !BOP

  ! !IROUTINE:  VTRISOLVE -- Solves for tridiagonal system that has been decomposed by VTRILU


  ! !INTERFACE:

  !  subroutine GEOS_TRISOLVE ( A,B,C,Y,YG )

  ! !ARGUMENTS:

  !    REAL(KIND=r8), dimension([:,[:,]] :),  intent(IN   ) ::  A, B, C
  !    REAL(KIND=r8), dimension([:,[:,]] :),  intent(INOUT) ::  Y
  !    REAL(KIND=r8), dimension([:,[:,]] :),  intent(IN   ) ::  YG

  ! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
  !   $LU x = f$. This is done by first solving $L g = f$ for $g$, and 
  !   then solving $U x = g$ for $x$. The solutions are:
  ! $$
  ! \begin{array}{rcl}
  ! g_1 & = & f_1, \\
  ! g_k & = & \makebox[2 in][l]{$f_k - g_{k-1} \hat{a}_{k}$,}  k=2, K, \\
  ! \end{array}
  ! $$
  ! and  
  ! $$
  ! \begin{array}{rcl}
  ! x_K & = & g_K /\hat{b}_K, \\
  ! x_k & = & \makebox[2 in][l]{($g_k - c_k g_{k+1}) / \hat{b}_{k}$,}  k=K-1, 1 \\
  ! \end{array}
  ! $$
  !  
  !  On input A contains the $\hat{a}_k$, the lower diagonal of $L$,
  !   B contains the $1/\hat{b}_k$, inverse of the  main diagonal of $U$,
  !   C contains the $c_k$, the upper diagonal of $U$. The forcing, $f_k$ is
  !   
  !   It returns the
  !   solution in the r.h.s input vector, Y. A has the multiplier from the
  !   decomposition, B the 
  !   matrix (U), and C the upper diagonal of the original matrix and of U.
  !   YG is the LM+1 (Ground) value of Y.

  !EOP



  !BOP

  ! !IROUTINE:  VTRILU --  Does LU decomposition of tridiagonal matrix.

  ! !INTERFACE:

  !  subroutine GEOS_TRILU  ( A,B,C )

  ! !ARGUMENTS:

  !    REAL(KIND=r8), dimension ([:,[:,]] :), intent(IN   ) ::  C
  !    REAL(KIND=r8), dimension ([:,[:,]] :), intent(INOUT) ::  A, B

  ! !DESCRIPTION: {\tt VTRILU} performs an $LU$ decomposition on
  ! a tridiagonal matrix $M=LU$.
  !
  ! $$
  ! M = \left( \begin{array}{ccccccc}
  !      b_1 & c_1 & & & & & \\
  !      a_2 & b_2 & c_2 & & & &  \\
  !      &  \cdot& \cdot & \cdot & & &  \\
  !      & & \cdot& \cdot & \cdot & &  \\
  !      &&  & \cdot& \cdot & \cdot &  \\
  !      &&&& a_{K-1} & b_{K-1} & c_{K-1}   \\
  !      &&&&& a_{K} & b_{K}
  !    \end{array} \right)
  ! $$
  !
  ! $$
  ! \begin{array}{lr}
  ! L = \left( \begin{array}{ccccccc}
  !      1 &&&&&& \\
  !      \hat{a}_2 & 1 & &&&&  \\
  !      &  \cdot& \cdot &  & & &  \\
  !      & & \cdot& \cdot &  &&  \\
  !      &&  & \cdot& \cdot &  &  \\
  !      &&&& \hat{a}_{K-1} & 1 &   \\
  !      &&&&& \hat{a}_{K} & 1
  !    \end{array} \right)
  ! &
  ! U = \left( \begin{array}{ccccccc}
  !      \hat{b}_1 & c_1 &&&&& \\
  !       & \hat{b}_2 & c_2 &&&&  \\
  !      &  & \cdot & \cdot & & &  \\
  !      & & & \cdot & \cdot &&  \\
  !      &&  & & \cdot & \cdot &  \\
  !      &&&&  & \hat{b}_{K-1} & c_{K-1}   \\
  !      &&&&&  & \hat{b}_{K}
  !    \end{array} \right)
  ! \end{array}
  ! $$
  !
  ! On input, A, B, and C contain, $a_k$, $b_k$, and $c_k$
  ! the lower, main, and upper diagonals of the matrix, respectively.
  ! On output, B contains $1/\hat{b}_k$, the inverse of the main diagonal of $U$,
  ! and A contains $\hat{a}_k$,
  ! the lower diagonal of $L$. C contains the upper diagonal of the original matrix and of $U$.
  !
  ! The new diagonals $\hat{a}_k$ and $\hat{b}_k$ are:
  ! $$
  ! \begin{array}{rcl}
  ! \hat{b}_1 & = & b_1, \\
  ! \hat{a}_k & = & \makebox[2 in][l]{$a_k / \hat{b}_{k-1}$,}  k=2, K, \\
  ! \hat{b}_k & = & \makebox[2 in][l]{$b_k - c_{k-1} \hat{a}_k$,} k=2, K. 
  ! \end{array}
  ! $$
  !EOP



  !#define DIMS
  SUBROUTINE GEOS_TRILU1 ( A,B,C )
    !#include "trilu.code"
    REAL(KIND=r8), DIMENSION ( :), INTENT(IN   ) ::  C
    REAL(KIND=r8), DIMENSION ( :), INTENT(INOUT) ::  A, B

    INTEGER :: LM, L

    LM = SIZE(A,SIZE(SHAPE(A)))

    B( 1) = 1.0_r8 /  B( 1)

    DO L = 2,LM
       A( L) = A( L) * B( L-1)
       B( L) = 1.0_r8 / ( B( L) - C( L-1) * A( L) )
    ENDDO

    RETURN

  END SUBROUTINE GEOS_TRILU1

  SUBROUTINE GEOS_TRISOLVE1 ( A,B,C,Y,YG )
    !#include "trisolve.code"
    REAL(KIND=r8), DIMENSION( :),  INTENT(IN   ) ::  A, B, C
    REAL(KIND=r8), DIMENSION( :),  INTENT(INOUT) ::  Y
    REAL(KIND=r8), DIMENSION( :),  INTENT(IN   ) ::  YG

    INTEGER :: LM, L

    LM = SIZE(A,SIZE(SHAPE(A)))

    ! Sweep down, modifying rhs with multiplier A

    Y( 1) = Y( 1) - A( 1)*YG( 1)

    DO L = 2,LM
       Y( L) = Y( L) - Y( L-1) * A( L)
    ENDDO

    ! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    Y( LM)   = (Y( LM) - C( LM) * YG(  2))*B( LM)

    DO L = LM-1,1,-1
       Y( L) = (Y( L ) - C( L ) * Y( L+1))*B( L )
    ENDDO

    RETURN

  END SUBROUTINE GEOS_TRISOLVE1
  !#undef DIMS

  !#define DIMS :,
  SUBROUTINE GEOS_TRILU2 ( A,B,C )
    !#include "trilu.code"
    REAL(KIND=r8), DIMENSION (:, :), INTENT(IN   ) ::  C
    REAL(KIND=r8), DIMENSION (:, :), INTENT(INOUT) ::  A, B

    INTEGER :: LM, L

    LM = SIZE(A,SIZE(SHAPE(A)))

    B(:, 1) = 1.0_r8 /  B(:, 1)

    DO L = 2,LM
       A(:, L) = A(:, L) * B(:, L-1)
       B(:, L) = 1.0_r8 / ( B(:, L) - C(:, L-1) * A(:, L) )
    ENDDO

    RETURN

  END SUBROUTINE GEOS_TRILU2

  SUBROUTINE GEOS_TRISOLVE2 ( A,B,C,Y,YG )
    !#include "trisolve.code"
    REAL(KIND=r8), DIMENSION(:, :),  INTENT(IN   ) ::  A, B, C
    REAL(KIND=r8), DIMENSION(:, :),  INTENT(INOUT) ::  Y
    REAL(KIND=r8), DIMENSION(:, :),  INTENT(IN   ) ::  YG

    INTEGER :: LM, L

    LM = SIZE(A,SIZE(SHAPE(A)))

    ! Sweep down, modifying rhs with multiplier A

    Y(:, 1) = Y(:, 1) - A(:, 1)*YG(:, 1)

    DO L = 2,LM
       Y(:, L) = Y(:, L) - Y(:, L-1) * A(:, L)
    ENDDO

    ! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    Y(:, LM)   = (Y(:, LM) - C(:, LM) * YG(:,  2))*B(:, LM)

    DO L = LM-1,1,-1
       Y(:, L) = (Y(:, L ) - C(:, L ) * Y(:, L+1))*B(:, L )
    ENDDO

    RETURN

  END SUBROUTINE GEOS_TRISOLVE2
  !#undef DIMS

  !#define DIMS :,:,
  SUBROUTINE GEOS_TRILU3 ( A,B,C )
    !#include "trilu.code"
    REAL(KIND=r8), DIMENSION (:,:, :), INTENT(IN   ) ::  C
    REAL(KIND=r8), DIMENSION (:,:, :), INTENT(INOUT) ::  A, B

    INTEGER :: LM, L

    LM = SIZE(A,SIZE(SHAPE(A)))

    B(:,:, 1) = 1.0_r8 /  B(:,:, 1)

    DO L = 2,LM
       A(:,:, L) = A(:,:, L) * B(:,:, L-1)
       B(:,:, L) = 1.0_r8 / ( B(:,:, L) - C(:,:, L-1) * A(:,:, L) )
    ENDDO

    RETURN

  END SUBROUTINE GEOS_TRILU3

  SUBROUTINE GEOS_TRISOLVE3 ( A,B,C,Y,YG )
    !#include "trisolve.code"
    REAL(KIND=r8), DIMENSION(:,:, :),  INTENT(IN   ) ::  A, B, C
    REAL(KIND=r8), DIMENSION(:,:, :),  INTENT(INOUT) ::  Y
    REAL(KIND=r8), DIMENSION(:,:, :),  INTENT(IN   ) ::  YG

    INTEGER :: LM, L

    LM = SIZE(A,SIZE(SHAPE(A)))

    ! Sweep down, modifying rhs with multiplier A

    Y(:,:, 1) = Y(:,:, 1) - A(:,:, 1)*YG(:,:, 1)

    DO L = 2,LM
       Y(:,:, L) = Y(:,:, L) - Y(:,:, L-1) * A(:,:, L)
    ENDDO

    ! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    Y(:,:, LM)   = (Y(:,:, LM) - C(:,:, LM) * YG(:,:,  2))*B(:,:, LM)

    DO L = LM-1,1,-1
       Y(:,:, L) = (Y(:,:, L ) - C(:,:, L ) * Y(:,:, L+1))*B(:,:, L )
    ENDDO

    RETURN

  END SUBROUTINE GEOS_TRISOLVE3
  !#undef DIMS



END MODULE GEOS_UtilsMod



MODULE PBL_Entrain
  USE Constants, ONLY :     &
       cp,            &
       grav,          &
       gasr,          &
       r8

  USE GEOS_UtilsMod
  IMPLICIT NONE
  PRIVATE

  !-----------------------------------------------------------------------
  !
  !      set default values to some instance independent parameters       
  !

  REAL(KIND=r8) :: akmax       =  1.e4_r8 ! maximum value for a diffusion coefficient 
  ! (m2/s)

  REAL(KIND=r8) :: zcldtopmax  =  3.e3_r8 ! maximum altitude for cloud top of 
  ! radiatively driven convection (m)    
  REAL(KIND=r8)   , PARAMETER, PUBLIC :: MAPL_H2OMW2  = 18.01_r8                  ! kg/Kmole
  REAL(KIND=r8)   , PARAMETER, PUBLIC :: MAPL_RUNIV  = 8314.3_r8                 ! J/(Kmole K)
  REAL(KIND=r8)   , PARAMETER, PUBLIC :: MAPL_AIRMW2  = 28.97_r8                  ! kg/Kmole
  REAL(KIND=r8)   , PARAMETER, PUBLIC :: MAPL_RGAS   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
  REAL(KIND=r8)   , PARAMETER, PUBLIC :: MAPL_ALHS   = 2.8368E6_r8               ! J/kg
  REAL(KIND=r8)   , PARAMETER, PUBLIC :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0_r8   ! --
  REAL(KIND=r8), parameter, public :: MAPL_KAPPA  = 2.0_r8/7.0_r8                ! --


  REAL(KIND=r8)   , PARAMETER  :: missing_value = 0.0_r8  ! MAPL_UNDEF
  REAL(KIND=r8)   , PARAMETER  :: tfreeze  =  273.16_r8                 ! K
  REAL(KIND=r8)   , PARAMETER  :: rdgas    =  MAPL_RGAS
  REAL(KIND=r8)   , PARAMETER  :: hls      =  MAPL_ALHS
  REAL(KIND=r8)   , PARAMETER  :: d608     =  MAPL_VIREPS

  REAL(KIND=r8)   , PARAMETER  ::   vonkarm  =  0.4_r8!MAPL_KARMAN
  REAL(KIND=r8)   , PARAMETER  ::   cp_air=  1004.4_r8!MAPL_KARMAN
  REAL(KIND=r8)   , PARAMETER  ::   hlv=2.6e6_r8
  REAL(KIND=r8)   , PARAMETER  ::   ramp     =  20.0_r8
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_GRAV   = 9.80_r8                   ! m^2/s
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_KARMAN = 0.40_r8                   ! --
  REAL(KIND=r8), parameter, public :: MAPL_P00    = 100000.0_r8               ! Pa
  REAL(KIND=r8), parameter, public :: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K)

  INTEGER, PARAMETER  :: ESMF_MAXSTR=255
  
  
  PUBLIC :: PBLNASA
CONTAINS


 SUBROUTINE PBLNASA(&
                 IM         , &!INTEGER      , INTENT(IN   ) :: IM
                 LM         , &!INTEGER      , INTENT(IN   ) :: LM
                 pcnst      , &!INTEGER      , INTENT(IN   ) :: pcnst
                 DT         , &!REAL(KIND=r8), INTENT(IN   ) :: DT
                 kvm        , &!REAL(KIND=r8), INTENT(INOUT) :: kvm        (IM,LM + 1) 
                 kvh        , &!REAL(KIND=r8), INTENT(INOUT) :: kvh        (IM,LM + 1 )
                 psomc      , &!REAL(KIND=r8), INTENT(IN   ) :: psomc      (IM,LM)        !(a
                 zm         , &!REAL(KIND=r8), INTENT(IN   ) :: zm        (IM,LM)             ! u wind input 
                 um1        , &!REAL(KIND=r8), INTENT(IN   ) :: um1       (IM,LM)             ! u wind input
                 vm1        , &!REAL(KIND=r8), INTENT(IN   ) :: vm1       (IM,LM)             ! v wind input
                 tm1        , &!REAL(KIND=r8), INTENT(IN   ) :: tm1       (IM,LM)             ! temperature input
                 qm1        , &!REAL(KIND=r8), INTENT(IN   ) :: qm1       (IM,LM,pcnst)  ! moisture and trace constituent input
                 pmidm1     , &!REAL(KIND=r8), INTENT(IN   ) :: pmidm1    (IM,LM)          ! midpoint pressures
                 pintm1     , &!REAL(KIND=r8), INTENT(IN   ) :: pintm1    (IM,LM + 1)    ! interface pressures
                 LwCoolRateC, &!REAL(KIND=r8), INTENT(IN   ) :: LwCoolRateC     (IM,LM) !clearsky_air_temperature_tendency_lw K s-1
                 LwCoolRate , &!REAL(KIND=r8), INTENT(IN   ) :: LwCoolRate      (IM,LM)!             air_temperature_tendency_due_to_longwave K s-1
                 cldtot     , &!REAL(KIND=r8), INTENT(IN   ) :: cldtot    (IM,LM)        
                 qliq       , &!REAL(KIND=r8), INTENT(IN   ) :: qliq      (IM,LM)        
                 bstar      , &!REAL(KIND=r8), INTENT(IN   ) :: bstar     (IM)   !surface_bouyancy_scale m s-2
                 ustar      , &!REAL(KIND=r8), INTENT(IN   ) :: USTAR      (IM) !surface_velocity_scale m s-1
                 FRLAND      ) !REAL(KIND=r8), INTENT(IN   ) :: FRLAND     (IM)
   IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: LM
    INTEGER      , INTENT(IN   ) :: pcnst
    REAL(KIND=r8), INTENT(IN   ) :: DT
    REAL(KIND=r8), INTENT(IN   ) :: bstar     (IM)   !surface_bouyancy_scale m s-2
    REAL(KIND=r8), INTENT(IN   ) :: USTAR      (IM) !surface_velocity_scale m s-1
    REAL(KIND=r8), INTENT(IN   ) :: FRLAND     (IM)
    REAL(KIND=r8), INTENT(IN   ) :: psomc      (IM,LM)        !(a
    REAL(KIND=r8), INTENT(INOUT) :: kvm        (IM,LM + 1) 
    REAL(KIND=r8), INTENT(INOUT) :: kvh        (IM,LM + 1 )
    REAL(KIND=r8), INTENT(IN   ) :: LwCoolRate      (IM,LM)!             air_temperature_tendency_due_to_longwave K s-1
    REAL(KIND=r8), INTENT(IN   ) :: LwCoolRateC     (IM,LM) !clearsky_air_temperature_tendency_lw K s-1
    REAL(KIND=r8), INTENT(IN   ) :: cldtot    (IM,LM)        
    REAL(KIND=r8), INTENT(IN   ) :: qliq      (IM,LM)        
    REAL(KIND=r8), INTENT(IN   ) :: zm        (IM,LM)             ! u wind input 
    REAL(KIND=r8), INTENT(IN   ) :: um1       (IM,LM)             ! u wind input
    REAL(KIND=r8), INTENT(IN   ) :: vm1       (IM,LM)             ! v wind input
    REAL(KIND=r8), INTENT(IN   ) :: tm1       (IM,LM)             ! temperature input
    REAL(KIND=r8), INTENT(IN   ) :: qm1       (IM,LM,pcnst)  ! moisture and trace constituent input
    REAL(KIND=r8), INTENT(IN   ) :: pmidm1    (IM,LM)          ! midpoint pressures
    REAL(KIND=r8), INTENT(IN   ) :: pintm1    (IM,LM + 1)    ! interface pressures
  
  REAL(KIND=r8) :: RADLW      (IM,LM)!        air_temperature_tendency_due_to_longwave K s-1
  REAL(KIND=r8) :: RADLWC     (IM,LM) !clearsky_air_temperature_tendency_lw K s-1
  REAL(KIND=r8) :: KH (IM,0:LM) !, INTENT(OUT)Heat diffusivity at base of each layer  (m+2 s-1).
  REAL(KIND=r8) :: KM (IM,0:LM) !, INTENT(OUT)Momentum diffusivity at base of each layer (m+2 s-1).
!  REAL(KIND=r8) :: DU (IM,0:LM) !, INTENT(OUT)Magnitude of wind shear (s-1).
  REAL(KIND=r8) :: Z  (IM,LM) !, INTENT(IN )Height of layer center above the surface (m).
  REAL(KIND=r8) :: ZLE(IM,0:LM) !, INTENT(IN )Height of layer base above the surface (m).
  REAL(KIND=r8) :: THV(IM,LM) !, INTENT(IN )Virtual potential temperature at layer center (K).
  REAL(KIND=r8) :: T  (IM,LM) !, INTENT(IN )Eastward velocity at layer center (m s-1).
  REAL(KIND=r8) :: U  (IM,LM) !, INTENT(IN )Eastward velocity at layer center (m s-1).
  REAL(KIND=r8) :: V  (IM,LM) !, INTENT(IN )Northward velocity at layer center (m s-1).
  REAL(KIND=r8) :: Q  (IM,LM) !, INTENT(IN )specific_humidity' ().
  REAL(KIND=r8) :: TH (IM,LM) !, INTENT(IN )potential_temperature (K).
  REAL(KIND=r8) :: PLE(IM,0:LM) !, INTENT(IN )air_pressure (Pa).
  REAL(KIND=r8) :: LOUIS        !, INTENT(IN )Louis scheme parameters (usually 5).
  REAL(KIND=r8) :: MINSHEAR     !, INTENT(IN )Min shear allowed in Ri calculation (s-1).
  REAL(KIND=r8) :: MINTHICK     !, INTENT(IN )Min layer thickness (m).
  REAL(KIND=r8) :: LAMBDAM       !, INTENT(IN )Blackadar(1962) length scale parameter for momentum (m).
  REAL(KIND=r8) :: LAMBDAM2      !, INTENT(IN )Second Blackadar parameter for momentum (m).
  REAL(KIND=r8) :: LAMBDAH    !, INTENT(IN )Blackadar(1962) length scale parameter for heat (m).
  REAL(KIND=r8) :: LAMBDAH2   !, INTENT(IN )Second Blackadar parameter for heat (m).
  REAL(KIND=r8) :: ZKMENV      !, INTENT(IN )Transition height for Blackadar param for momentum (m)
  REAL(KIND=r8) :: ZKHENV      !, INTENT(IN )Transition height for Blackadar param for heat        (m)
  REAL(KIND=r8) :: KHMMAX      !, INTENT(IN )Maximum allowe diffusivity (m+2 s-1).
  REAL(KIND=r8) :: PRANDTLSFC 
  REAL(KIND=r8) :: PRANDTLRAD 
  REAL(KIND=r8) :: BETA_SURF 
  REAL(KIND=r8) :: BETA_RAD  
  REAL(KIND=r8) :: TPFAC_SURF
  REAL(KIND=r8) :: ENTRATE_SURF 
  REAL(KIND=r8) :: PCEFF_SURF   
  REAL(KIND=r8) :: KHRADFAC    
  REAL(KIND=r8) :: KHSFCFAC    
  REAL(KIND=r8) :: EKM                  (IM,LM)
  REAL(KIND=r8) :: EKH                  (IM,LM)
  REAL(KIND=r8) :: KHSFC           (IM,LM)
  REAL(KIND=r8) :: KHRAD           (IM,LM)
  REAL(KIND=r8) :: ZCLD            (IM)    
  REAL(KIND=r8) :: ZRADML          (IM) 
  REAL(KIND=r8) :: ZRADBS          (IM)
  REAL(KIND=r8) :: ZSML            (IM)
  REAL(KIND=r8) :: ZCLDTOP          (IM)
  REAL(KIND=r8) :: WESFC           (IM)
  REAL(KIND=r8) :: WERAD           (IM)
  REAL(KIND=r8) :: DBUOY           (IM)
  REAL(KIND=r8) :: VSCSFC          (IM)
  REAL(KIND=r8) :: VSCRAD          (IM)
  REAL(KIND=r8) :: KERAD           (IM)
  REAL(KIND=r8) :: VSCBRV          (IM)
  REAL(KIND=r8) :: WEBRV           (IM)
  REAL(KIND=r8) :: DSIEMS          (IM)
  REAL(KIND=r8) :: CHIS            (IM)
  REAL(KIND=r8) :: DELSINV          (IM)
  REAL(KIND=r8) :: SMIXT           (IM)
  REAL(KIND=r8) :: CLDRF           (IM)
  REAL(KIND=r8) :: RADRCODE          (IM)
  !REAL(KIND=r8) :: ALH_DIAG        (IM,0:LM)  ! Blackadar Length Scale diagnostic (m) [Optional] 
  REAL(KIND=r8) :: LWCRT(IM,1:LM)  !cloudy_LW_radiation_tendency_used_by_Lock_scheme K s-1

  !   TH = TH + DTS
  !  QH = QH + DQS

  !SH sensible_heat_flux W m-2
  !
  !DSH derivative_of_sensible_heat_wrt_dry_static_energy W m-2 K-1
  !
  !EVAP evaporation kg m-2 s-1
  !
  !DEVAP derivative_of_evaporation_wrt_QS kg m-2 s-1
  !
  !TA     surface_air_temperature (K)
  !
  !TA = TH - (SH/MAPL_CP + DSH  *DTS)/CT
  !
  !QH -> QHAT  effective_surface_specific_humidity
  !
  !QA = QH - (EVAP       + DEVAP*DQS)/CQ
  !
  !UU =>speed surface_wind_speedm s-1
  !
  !CQ surface_exchange_coefficient_for_moisture'kg m-2 s-1'
  !
  !CM surface_exchange_coefficient_for_momentum'kg m-2 s-1'
  !
  !DZ surface_layer_height'm',      
  !
  !BSTAR = (MAPL_GRAV/(RHOS*sqrt(CM*max(UU,1.e-30)/RHOS))) * (CT*(TH-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ*(QH-QA))
  REAL(KIND=r8), dimension(IM,LM)           ::  TV, DMI, PLO, QS, DQS, QL, QI, QA
  REAL(KIND=r8), dimension(IM,0:LM)         :: PKE
  REAL(KIND=r8), dimension(IM,1:LM-1)       :: TVE, RDZ
  integer                             :: LOCK_ON

  INTEGER :: RC,L,I,K
  RC=0
  LOUIS       =5.0_r8            !, INTENT(IN )Louis scheme parameters (usually 5).
  MINSHEAR     =0.0030_r8  !, INTENT(IN )Min shear allowed in Ri calculation (s-1).
  MINTHICK    =0.1_r8           !, INTENT(IN )Min layer thickness (m).
  LAMBDAM        =160.0_r8    !, INTENT(IN )Blackadar(1962) length scale parameter for momentum (m).
  LAMBDAM2        =160.0_r8    !, INTENT(IN )Second Blackadar parameter for momentum (m).
  LAMBDAH    =160.0_r8   !, INTENT(IN )Blackadar(1962) length scale parameter for heat (m).
  LAMBDAH2   =160.0_r8   !, INTENT(IN )Second Blackadar parameter for heat (m).
  ZKMENV     =3000.0_r8    !, INTENT(IN )Transition height for Blackadar param for momentum (m)
  ZKHENV     =3000.0_r8    !, INTENT(IN )Transition height for Blackadar param for heat     (m)
  KHMMAX     =100.0_r8     !, INTENT(IN )Maximum allowe diffusivity (m+2 s-1).
  PRANDTLSFC =1.0_r8
  PRANDTLRAD =0.75_r8
  BETA_SURF    =0.23_r8
  BETA_RAD     =0.23_r8
  TPFAC_SURF   =10.0_r8
  ENTRATE_SURF =1.e-3_r8
  PCEFF_SURF   =0.05_r8
  KHRADFAC    =0.85_r8
  KHSFCFAC    =0.85_r8

  ! LOUIS        =5.0_r8
  ! LAMBDAM      =160.0_r8
  ! LAMBDAM2     =160.0_r8
  ! LAMBDAH      =160.0_r8
  ! LAMBDAH2     =160.0_r8
  ! ZKMENV       =3000.0_r8
  ! ZKHENV       =3000.0_r8
  ! MINTHICK     =0.1_r8
  ! MINSHEAR     =0.0030_r8
  ! C_B          =2.50E-6_r8
  ! LAMBDA_B     =1500._r8
  ! AKHMMAX      =100._r8
   LOCK_ON      =1       
  ! PRANDTLSFC   =1.0_r8
  ! PRANDTLRAD   =0.75_r8
  ! BETA_RAD          =0.23_r8
  ! BETA_SURF         =0.23_r8
  ! KHRADFAC          =0.85_r8
  ! KHSFCFAC          =0.85_r8
  ! TPFAC_SURF        =10.0_r8
  ! ENTRATE_SURF =1.e-3_r8
  ! PCEFF_SURF   =0.05_r8
  DO L=1,LM
     DO I=1,IM          
        T     (i,L  ) = tm1(i,L)
        Q     (i,L  ) = qm1(I,L,1)!, INTENT(IN )specific_humidity' ().
        U     (i,L  ) = um1(I,L)        ! u wind input
        V     (i,L  ) = vm1(I,L)        ! v wind input
        RADLW (I,L  ) = LwCoolRate (i,LM+1-L ) 
        RADLWC(I,L  ) = LwCoolRateC(i,LM+1-L ) 
        IF(pcnst>1)THEN
           IF(pcnst == 2) QL    (I,L  ) = qm1(I,L,2)
           IF(pcnst == 3) QI    (I,L  ) = qm1(I,L,3)
        ELSE
           QL    (I,L  ) = 0.0_r8 
           QI    (I,L  ) = 0.0_r8
        END IF
        QA    (I,L  ) = cldtot     (I,LM+1-L)
     END DO
  END DO
  DO L=0,LM
     DO I=1,IM          
        KH (I,L) =   kvh  (I,L + 1)
        KM (I,L) =   kvm  (I,L + 1)
        PLE(I,L) =  pintm1(i,L + 1) !, INTENT(IN )air_pressure (Pa).
     END DO
  END DO
  DO L=1,LM
     DO I=1,IM          
             TH (i,L) =T (i,L) * psomc (i,L)!(MAPL_P00/PLE (i,L))**MAPL_KAPPA
     END DO
  END DO
!  DO  L=1,LM
!     DO I=1,IM          
!       WRITE(*,'(I5,4F12.5)')L,kvh  (I,L + 1), kvm  (I,L + 1)
!     END DO
!  END DO

! Compute the edge heights using Arakawa-Suarez hydrostatic equation
!---------------------------------------------------------------------------
  PKE(:,LM) = 1.0_r8 
  DO L=0,LM-1
     DO I=1,IM
              PKE(i,l) = psomc (i,L+1)
     END DO
  END DO
 
  ZLE(:,LM) = 0.00005_r8
  do L=LM,1,-1
     DO I=1,IM
             ZLE(i,L-1) = zm(i,L)
     END DO
  end do
  
! Layer height, pressure, and virtual temperatures
!-------------------------------------------------
! First add up clouds to get total fraction, ice, and liquid
   !  QL  = QLCN + QLLS
   !  QI  = QICN + QILS
   !  QA  = CLCN + CLLS
   !  QI(:,:,1:LM) =0.0
  do L=LM,1,-1
     DO I=1,IM
          Z  (I,L)   = 0.5_r8*(ZLE(I,L-1) + ZLE(I,L))
          PLO(I,L)   = 0.5_r8*(PLE(I,L-1) + PLE(I,L))
          TV (I,L)   = T (I,L) *( 1.0_r8 + MAPL_VIREPS *Q(I,L) - QL(I,L) - QI(I,L) )
          THV(I,L)   = TV(I,L)*(TH(I,L)/T(I,L))
      END DO
  END DO

  do L=LM-1,1,-1
      DO I=1,IM
          TVE(i,L) = (TV (i,L) + TV (i,L+1))*0.5_r8
      END DO
  END DO

! Miscellaneous factors
!----------------------

  do L=LM-1,1,-1
    DO I=1,IM
          RDZ(I,L) = PLE(I,L) / ( MAPL_RGAS * TVE(I,L) )
          RDZ(I,L) = RDZ(I,L) / (Z(I, L)-Z(I,L+1))
     END DO
  END DO

      do L=LM,1,-1
        DO I=1,IM
            DMI(I,L ) = (MAPL_GRAV*DT)/(PLE(i,L)-PLE(i,L-1))
         END DO
      END DO

! Need DQSAT for Lock scheme
      DO k=1,LM
         DO I=1,IM
            DQS(I,K)      = GEOS_DQSAT(T(I,K), PLO(I,K), qsat=QS(I,K), PASCALS=.true. ) 
         END DO
      END DO
!===> Running 1-2-1 smooth of bottom 5 levels of Virtual Pot. Temp.
      DO I=1,IM
         THV(I,LM) = THV(I,LM-1)*0.25_r8 + THV(I,LM)*0.75_r8
      END DO

      do L = LM-1,LM-5,-1
         DO I=1,IM
            THV(I,L) = THV(I,L-1)*0.25_r8 + THV(I,L)*0.50_r8 + THV(I,L+1)*0.25_r8 
         END DO
      END DO

!   ...then add Lock.
!--------------------
      do L = 1,LM
         DO I=1,IM
             LWCRT(i,l) = RADLW(i,l) - RADLWC(i,l)
         END DO
      END DO
   CALL entrain(  &
       IM                           ,& ! INTENT(in)     :: IM 
       LM                          ,& ! INTENT(in)     :: LM
       RADLW       (1:IM,1:LM)  ,& ! INTENT(in)     :: RADLW         (:,:,:) 
       USTAR       (1:IM     )  ,& ! INTENT(in)     :: USTAR         (:,:)   
       BSTAR       (1:IM     )  ,& ! INTENT(in)     :: BSTAR         (:,:)   
       FRLAND      (1:IM     )  ,& ! INTENT(in)     :: FRLAND         (:,:)   
       T           (1:IM,1:LM)  ,& ! INTENT(in)        :: T                 (:,:,:) 
       Q           (1:IM,1:LM)  ,& ! INTENT(in)        :: Q                 (:,:,:) 
       QL          (1:IM,1:LM)  ,& ! INTENT(in)        :: QL                 (:,:,:) 
       QI          (1:IM,1:LM)  ,& ! INTENT(in)        :: QI                 (:,:,:) 
       QA          (1:IM,1:LM)  ,& ! INTENT(in)        :: QA                 (:,:,:) 
       U           (1:IM,1:LM)  ,& ! INTENT(in)     :: U              (:,:,:) 
       V           (1:IM,1:LM)  ,& ! INTENT(in)        :: V                (:,:,:) 
       Z           (1:IM,1:LM)  ,& ! INTENT(in)        :: Z                (:,:,:) 
       PLO         (1:IM,1:LM)  ,& ! INTENT(in)        :: PLO                (:,:,:) 
       ZLE           (1:IM,0:LM)  ,& ! INTENT(in)        :: ZLE                (:,:,:) 
       PLE           (1:IM,0:LM)  ,& ! INTENT(in)        :: PLE                (:,:,:) 
!       QS           (1:IM,1:LM)  ,& ! INTENT(in)        :: QS                (:,:,:)     
       DQS           (1:IM,1:LM)  ,& ! INTENT(in)        :: DQS                (:,:,:) 
       KM           (1:IM,0:LM)  ,& ! INTENT(inout)  :: KM                (:,:,:) 
       KH           (1:IM,0:LM)  ,& ! INTENT(inout)  :: KH                (:,:,:) 
       EKM           (1:IM,1:LM)  ,& ! INTENT(out)       :: EKM                (:,:,:) 
       EKH           (1:IM,1:LM)  ,& ! INTENT(out)       :: EKH                (:,:,:) 
       KHSFC           (1:IM,1:LM)  ,& ! INTENT(out)       :: KHSFC         (:,:,:) 
       KHRAD           (1:IM,1:LM)  ,& ! INTENT(out)       :: KHRAD         (:,:,:) 
       ZCLD           (1:IM)        ,& ! INTENT(out)       :: ZCLD            (:,:)      
       ZRADML           (1:IM)        ,& ! INTENT(out)       :: ZRADML          (:,:)   
       ZRADBS           (1:IM)        ,& ! INTENT(out)       :: ZRADBS          (:,:)
       ZSML           (1:IM)        ,& ! INTENT(out)       :: ZSML            (:,:)  
       ZCLDTOP     (1:IM)        ,& ! INTENT(out)       :: ZCLDTOP          (:,:) 
       WESFC           (1:IM)        ,& ! INTENT(out)       :: WESFC           (:,:)
       WERAD           (1:IM)        ,& ! INTENT(out)       :: WERAD           (:,:)
       DBUOY           (1:IM)        ,& ! INTENT(out)       :: DBUOY           (:,:)
       VSCSFC           (1:IM)        ,& ! INTENT(out)       :: VSCSFC          (:,:)
       VSCRAD           (1:IM)        ,& ! INTENT(out)       :: VSCRAD          (:,:)
       KERAD           (1:IM)        ,& ! INTENT(out)       :: KERAD           (:,:)
       VSCBRV           (1:IM)        ,& ! INTENT(out)       :: VSCBRV          (:,:)
       WEBRV           (1:IM)        ,& ! INTENT(out)       :: WEBRV           (:,:)
       DSIEMS           (1:IM)        ,& ! INTENT(out)       :: DSIEMS          (:,:)
       CHIS           (1:IM)        ,& ! INTENT(out)       :: CHIS            (:,:)
       DELSINV     (1:IM)        ,& ! INTENT(out)       :: DELSINV          (:,:)
       SMIXT           (1:IM)        ,& ! INTENT(out)       :: SMIXT           (:,:)
       CLDRF           (1:IM)        ,& ! INTENT(out)       :: CLDRF           (:,:)
       RADRCODE    (1:IM)        ,& ! INTENT(out)       :: RADRCODE          (:,:)
       PRANDTLSFC                    ,& ! INTENT(in)         :: PRANDTLSFC  
       PRANDTLRAD                    ,& ! INTENT(in)         :: PRANDTLRAD  
       BETA_SURF                     ,& ! INTENT(in)         :: BETA_SURF        
       BETA_RAD                      ,& ! INTENT(in)         :: BETA_RAD        
       TPFAC_SURF                    ,& ! INTENT(in)         :: TPFAC_SURF  
       ENTRATE_SURF                  ,& ! INTENT(in)         :: ENTRATE_SURF
       PCEFF_SURF                    ,& ! INTENT(in)         :: PCEFF_SURF  
       KHRADFAC                      ,& ! INTENT(in)         :: KHRADFAC        
       KHSFCFAC                       ) ! INTENT(in)         :: KHSFCFAC        

  DO L=0,LM
     DO I=1,IM          
          kvm(I,L + 1)  = KM (I,L)  
          kvh(I,L + 1)  = KH (I,L)  
     END DO
  END DO

 END SUBROUTINE PBLNASA
 

  !
  SUBROUTINE entrain(  &
       IM             ,& !           INTENT(in)          :: ie                                          
       LM             ,&                                                                          
       tdtlw_in       ,& !         INTENT(in)     :: tdtlw_in  (:,:,:)         (1:IM,1:JM,1:LM)  
       u_star         ,& !         INTENT(in)     :: u_star         (:,:)         (1:IM,1:JM     )  
       b_star         ,& !         INTENT(in)     :: b_star         (:,:)         (1:IM,1:JM     )  
       frland         ,& !         INTENT(in)     :: frland         (:,:)         (1:IM,1:JM     )  
       t              ,& !         INTENT(in)     :: t         (:,:,:)                (1:IM,1:JM,1:LM)  
       qv             ,& !         INTENT(in)     :: qv         (:,:,:)                (1:IM,1:JM,1:LM)  
       ql             ,& !         INTENT(in)     :: ql         (:,:,:)                (1:IM,1:JM,1:LM)  
       qi             ,& !         INTENT(in)     :: qi         (:,:,:)                (1:IM,1:JM,1:LM)  
       qa             ,& !         INTENT(in)     :: qa         (:,:,:)                (1:IM,1:JM,1:LM)  
       u              ,& !           INTENT(in)          :: u   (:,:,:)                (1:IM,1:JM,1:LM)  
       v              ,& !           INTENT(in)          :: v   (:,:,:)                (1:IM,1:JM,1:LM)  
       zfull          ,& !         INTENT(in)     :: zfull         (:,:,:)       (1:IM,1:JM,1:LM)  
       pfull          ,& !         INTENT(in)     :: pfull         (:,:,:)       (1:IM,1:JM,1:LM)  
       zhalf          ,& !         INTENT(in)          :: zhalf         (:,:,:)       (1:IM,1:JM,0:LM)  
       phalf          ,& !         INTENT(in)     :: phalf         (:,:,:)       (1:IM,1:JM,0:LM)  
!       qs             ,& !         INTENT(in)          :: qs          (:,:,:)       (1:IM,1:JM,1:LM)  
       dqs            ,& !         INTENT(in)          :: dqs          (:,:,:)       (1:IM,1:JM,1:LM)  
       diff_m         ,& !         INTENT(inout)  :: diff_m   (:,:,:)                (1:IM,1:JM,0:LM)  
       diff_t         ,& !         INTENT(inout)  :: diff_t   (:,:,:)                (1:IM,1:JM,0:LM)  
       k_m_entr       ,& !         INTENT(out)    :: k_m_entr  (:,:,:)         (1:IM,1:JM,1:LM)  
       k_t_entr       ,& !         INTENT(out)    :: k_t_entr  (:,:,:)         (1:IM,1:JM,1:LM)  
       k_sfc_diag     ,& !         INTENT(out)    :: k_rad_diag(:,:,:)         (1:IM,1:JM,1:LM)  
       k_rad_diag     ,& !         INTENT(out)    :: k_sfc_diag(:,:,:)         (1:IM,1:JM,1:LM)  
       zcloud         ,& !         INTENT(out)    :: zcloud   (:,:)                  (1:IM,1:JM)         
       zradml         ,& !         INTENT(out)    :: zradml   (:,:)                  (1:IM,1:JM)         
       zradbase       ,& !         INTENT(out)    :: zradbase  (:,:)               (1:IM,1:JM)         
       zsml           ,& !         INTENT(out)    :: zsml     (:,:)                 (1:IM,1:JM)         
       zcldtop_diag   ,& !         INTENT(out)    :: zcldtop_diag    (:,:)     (1:IM,1:JM)         
       wentr_sfc_diag ,& !         INTENT(out)    :: wentr_sfc_diag  (:,:)     (1:IM,1:JM)         
       wentr_rad_diag ,& !         INTENT(out)    :: wentr_rad_diag  (:,:)     (1:IM,1:JM)         
       del_buoy_diag  ,& !           INTENT(out)    :: del_buoy_diag   (:,:)     (1:IM,1:JM)         
       vsfc_diag      ,& !         INTENT(out)    :: vsfc_diag         (:,:)               (1:IM,1:JM)         
       vrad_diag      ,& !           INTENT(out)    :: vrad_diag        (:,:)               (1:IM,1:JM)         
       kentrad_diag   ,& !           INTENT(out)    :: kentrad_diag         (:,:) (1:IM,1:JM)         
       vbrv_diag      ,& !           INTENT(out)    :: vbrv_diag        (:,:)               (1:IM,1:JM)         
       wentr_brv_diag ,& !           INTENT(out)    :: wentr_brv_diag  (:,:)     (1:IM,1:JM)         
       dsiems_diag    ,& !           INTENT(out)    :: dsiems_diag         (:,:) (1:IM,1:JM)         
       chis_diag      ,& !           INTENT(out)    :: chis_diag        (:,:)               (1:IM,1:JM)         
       delsinv_diag   ,& !         INTENT(out)    :: delsinv_diag    (:,:)     (1:IM,1:JM)         
       slmixture_diag ,& !         INTENT(out)    :: slmixture_diag  (:,:)     (1:IM,1:JM)         
       cldradf_diag   ,& !         INTENT(out)    :: cldradf_diag    (:,:)     (1:IM,1:JM)         
       radrcode       ,& !         INTENT(out)    :: radrcode         (:,:)               (1:IM,1:JM)         
       prandtlsfc     ,& !         INTENT(in)          :: prandtlsfc                                  
       prandtlrad     ,& !         INTENT(in)          :: prandtlrad                                  
       beta_surf      ,& !         INTENT(in)          :: beta_surf                                         
       beta_rad       ,& !         INTENT(in)          :: beta_rad                                         
       tpfac_sfc      ,& !         INTENT(in)          :: tpfac_sfc                                         
       entrate_sfc    ,& !         INTENT(in)          :: entrate_sfc                                 
       pceff_sfc      ,& !         INTENT(in)          :: pceff_sfc                                         
       khradfac       ,& !         INTENT(in)          :: khradfac                                         
       khsfcfac       ,& !         INTENT(in)          :: khsfcfac                                         
       kbot            ) ! INTENT(in ) , OPTIONAL :: kbot(:,:)

    !-----------------------------------------------------------------------
    !
    !      variables
    !
    !      -----
    !      input
    !      -----
    !
    !      is,ie,js,je  i,j indices marking the slab of model working on                             
    !      time      variable needed for netcdf diagnostics !!!*** REMOVED JTB 12/18/2003
    !
    !      convect   is surface based moist convection occurring in this
    !                grid box?
    !      u_star    friction velocity (m/s)
    !      b_star    buoyancy scale (m/s**2)
    !
    !      three dimensional fields on model full levels, REAL(KIND=r8)s dimensioned
    !      (:,:,pressure), third index running from top of atmosphere to 
    !      bottom
    !          
    !      t         temperature (K)
    !      qv        water vapor specific humidity  (kg vapor/kg air)
    !      ql        liquid water specific humidity (kg cond/kg air)
    !      qi        ice water specific humidity    (kg cond/kg air)
    !      qa        cloud fraction 
    !      qs        saturation specific humidity (kg/kg)
    !      dqs       derivative of qs w/ respect to temp 
    !      zfull     height of full levels (m)
    !      pfull     pressure (Pa)
    !      u         zonal wind (m/s)
    !      v         meridional wind (m/s)
    !
    !      the following two fields are on the model half levels, with
    !      size(zhalf,3) = size(t,3) +1, zhalf(:,:,size(zhalf,3)) 
    !      must be height of surface (if you are not using eta-model)
    !
    !      zhalf     height at half levels (m)
    !      phalf     pressure at half levels (Pa)
    !
    !      ------------
    !      input/output
    !      ------------
    !
    !      the following variables are defined at half levels and are
    !      dimensions 1:nlev
    !
    !      diff_t   input and output heat diffusivity (m2/sec)
    !      diff_m   input and output momentum diffusivity (m2/sec)
    !
    !      The diffusivity coefficient output from the routine includes
    !      the modifications to use the internally calculated diffusivity
    !      coefficients.
    !
    !      ------
    !      output
    !      ------
    !
    !      The following variables are defined at half levels and are
    !      dimensions 1:nlev.
    !
    !      k_t_entr  heat diffusivity coefficient (m**2/s)
    !      k_m_entr  momentum diffusivity coefficient (m**2/s)
    !      zsml      height of surface driven mixed layer (m)
    !
    !      --------------
    !      optional input
    !      --------------
    !
    !      kbot      integer indicating the lowest true layer of atmosphere
    !                this is used only for eta coordinate model
    !
    !      --------
    !      internal
    !      --------
    !
    !
    !      General variables
    !      -----------------
    !
    !      slv         virtual static energy (J/kg)      
    !      density     air density (kg/m3)
    !      hleff       effective latent heat of vaporization/sublimation 
    !                  (J/kg)
    !
    !
    !
    !      Variables related to surface driven convective layers
    !      -----------------------------------------------------
    !
    !      vsurf       surface driven buoyancy velocity scale (m/s)
    !      vshear      surface driven shear velocity scale (m/s)
    !      wentr_pbl   surface driven entrainment rate (m/s)
    !      convpbl     1 is surface driven convective layer present
    !                  0 otherwise
    !      pblfq       1 if the half level is part of a surface driven
    !                  layer, 0 otherwise
    !      k_m_troen   momentum diffusion coefficient (m2/s)
    !      k_t_troen   heat diffusion coefficient (m2/s)
    !
    !
    !      Variables related to cloud top driven radiatively driven layers
    !      ---------------------------------------------------------------
    !
    !      zradbase    height of base of radiatively driven mixed layer (m)
    !      zradtop     height of top of radiatively driven mixed layer (m)
    !      zradml      depth of radiatively driven mixed layer (m)
    !      vrad        radiatively driven velocity scale (m/s)
    !      radf        longwave jump at cloud top (W/m2) -- the radiative 
    !                  forcing for cloud top driven mixing.
    !      wentr_rad   cloud top driven entrainment (m/s)
    !      svpcp       cloud top value of liquid water virtual static energy 
    !                  divided by cp (K)
    !      radpbl      1 if cloud top radiatively driven layer is present
    !                  0 otherwise
    !      radfq       1 if the half level is part of a radiatively driven
    !                  layer, 0 otherwise
    !      k_rad       radiatively driven diffusion coefficient (m2/s)
    !
    !      Diagnostic variables
    !      --------------------
    !
    !      fqinv       1 if an inversion occurs at altitudes less 
    !                    than 3000 m, 0 otherwise
    !      zinv        altitude of inversion base (m)
    !      invstr      strength of inversion in slv/cp (K)
    !
    !-----------------------------------------------------------------------

    INTEGER,         INTENT(in)     :: IM
    INTEGER,         INTENT(in)     :: LM
    REAL(KIND=r8),            INTENT(in)     :: tdtlw_in  (1:IM,1:LM)   
    REAL(KIND=r8),            INTENT(in)     :: u_star         (1:IM     )  
    REAL(KIND=r8),            INTENT(in)     :: b_star         (1:IM     )  
    REAL(KIND=r8),            INTENT(in)     :: frland         (1:IM     )  
    REAL(KIND=r8),            INTENT(in)     :: t         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: qv         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: ql         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: qi         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: qa         (1:IM,1:LM)
!    REAL(KIND=r8),            INTENT(in)     :: qs         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: dqs         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: u         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: v         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: zfull         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: pfull         (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(in)     :: zhalf         (1:IM,0:LM) 
    REAL(KIND=r8),            INTENT(in)     :: phalf         (1:IM,0:LM) 
    REAL(KIND=r8),            INTENT(inout)  :: diff_m         (1:IM,0:LM)
    REAL(KIND=r8),            INTENT(inout)  :: diff_t         (1:IM,0:LM)
    REAL(KIND=r8),            INTENT(out)    :: k_m_entr       (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(out)    :: k_t_entr       (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(out)    :: k_rad_diag     (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(out)    :: k_sfc_diag     (1:IM,1:LM)
    REAL(KIND=r8),            INTENT(out)    :: zsml           (1:IM)
    REAL(KIND=r8),            INTENT(out)    :: zradml         (1:IM)
    REAL(KIND=r8),            INTENT(out)    :: zcloud         (1:IM) 
    REAL(KIND=r8),            INTENT(out)    :: zradbase       (1:IM) 
    INTEGER        , INTENT(in ) , OPTIONAL :: kbot(1:IM)
    REAL(KIND=r8),            INTENT(in)     :: prandtlsfc
    REAL(KIND=r8),            INTENT(in)     :: prandtlrad
    REAL(KIND=r8),            INTENT(in)     :: beta_surf
    REAL(KIND=r8),            INTENT(in)     :: beta_rad
    REAL(KIND=r8),            INTENT(in)     :: khradfac
    REAL(KIND=r8),            INTENT(in)     :: tpfac_sfc
    REAL(KIND=r8),            INTENT(in)     :: entrate_sfc
    REAL(KIND=r8),            INTENT(in)     :: pceff_sfc
    REAL(KIND=r8),            INTENT(in)     :: khsfcfac


    REAL(KIND=r8),            INTENT(out) :: wentr_rad_diag  (1:IM)
    REAL(KIND=r8),            INTENT(out) :: wentr_sfc_diag  (1:IM)
    REAL(KIND=r8),            INTENT(out) :: del_buoy_diag   (1:IM)
    REAL(KIND=r8),            INTENT(out) :: vrad_diag       (1:IM)
    REAL(KIND=r8),            INTENT(out) :: kentrad_diag    (1:IM)
    REAL(KIND=r8),            INTENT(out) :: vbrv_diag       (1:IM)
    REAL(KIND=r8),            INTENT(out) :: wentr_brv_diag  (1:IM)
    REAL(KIND=r8),            INTENT(out) :: dsiems_diag     (1:IM)
    REAL(KIND=r8),            INTENT(out) :: chis_diag       (1:IM)
    REAL(KIND=r8),            INTENT(out) :: zcldtop_diag    (1:IM)
    REAL(KIND=r8),            INTENT(out) :: delsinv_diag    (1:IM)
    REAL(KIND=r8),            INTENT(out) :: slmixture_diag  (1:IM)
    REAL(KIND=r8),            INTENT(out) :: cldradf_diag    (1:IM)
    REAL(KIND=r8),            INTENT(out) :: vsfc_diag       (1:IM)
    REAL(KIND=r8),            INTENT(out) :: radrcode        (1:IM)



    INTEGER  :: i
    INTEGER  :: j
    INTEGER  :: k
    INTEGER  :: ibot
    INTEGER  :: itmp
    INTEGER  :: nlev
    INTEGER  :: nlat
    INTEGER  :: nlon
    INTEGER  :: ipbl (1:IM)
    INTEGER  :: kmax (1:IM)
    INTEGER  :: kcldtop(1:IM)
    INTEGER  :: kcldbot(1:IM)
    INTEGER  :: kcldtop2(1:IM)
    !LOGICAL  :: used
    LOGICAL  :: keeplook
    LOGICAL  :: do_jump_exit(1:IM)
    REAL(KIND=r8)     :: maxradf(1:IM)
    !REAL(KIND=r8)     :: tmpradf
    REAL(KIND=r8)     :: stab
    REAL(KIND=r8)     :: maxqdot
    !REAL(KIND=r8)     :: tmpqdot
    REAL(KIND=r8)     :: wentrmax
    REAL(KIND=r8)     :: maxinv
    REAL(KIND=r8)     :: qlcrit
    REAL(KIND=r8)     :: ql00
    REAL(KIND=r8)     :: qlm1
    REAL(KIND=r8)     :: Abuoy
    REAL(KIND=r8)     :: Ashear
    REAL(KIND=r8)     :: wentr_tmp(1:IM)
    REAL(KIND=r8)     :: hlf
    REAL(KIND=r8)     :: k_entr_tmp(1:IM)
    REAL(KIND=r8)     :: tmpjump
    REAL(KIND=r8)     :: critjump(1:IM)
    REAL(KIND=r8)     :: radperturb(1:IM)
    REAL(KIND=r8)     :: buoypert(1:IM)
    REAL(KIND=r8)     :: tmp1(1:IM)
    REAL(KIND=r8)     :: tmp2(1:IM)
    REAL(KIND=r8)     :: slmixture(1:IM)
    REAL(KIND=r8)     :: vsurf3     (1:IM)
    REAL(KIND=r8)     :: vshear3    (1:IM)
    REAL(KIND=r8)     :: vrad3 (1:IM)
    REAL(KIND=r8)     :: vbr3  (1:IM)
    REAL(KIND=r8)     :: dsiems(1:IM)
    REAL(KIND=r8)     :: dslvcptmp
    REAL(KIND=r8)     :: ztmp
    REAL(KIND=r8)     :: zradtop   (1:IM)
    REAL(KIND=r8)     :: vrad      (1:IM)
    REAL(KIND=r8)     :: radf      (1:IM)
    REAL(KIND=r8)     :: svpcp     (1:IM)
    REAL(KIND=r8)     :: chis      (1:IM)
    REAL(KIND=r8)     :: vbrv      (1:IM)
    REAL(KIND=r8)     :: vsurf     (1:IM)
    REAL(KIND=r8)     :: vshear    (1:IM)
    REAL(KIND=r8)     :: wentr_rad (1:IM)
    REAL(KIND=r8)     :: wentr_pbl (1:IM)
    REAL(KIND=r8)     :: wentr_brv (1:IM)
    REAL(KIND=r8)     :: convpbl   (1:IM)
    REAL(KIND=r8)     :: radpbl    (1:IM)
    REAL(KIND=r8)     :: zinv      (1:IM)
    REAL(KIND=r8)     :: fqinv     (1:IM)
    REAL(KIND=r8)     :: invstr    (1:IM)
    REAL(KIND=r8)     :: slv       (1:IM,1:LM)
    REAL(KIND=r8)     :: density   (1:IM,1:LM)
    REAL(KIND=r8)     :: qc        (1:IM,1:LM)
    REAL(KIND=r8)     :: hleff     (1:IM,1:LM)
    REAL(KIND=r8)     :: radfq     (1:IM,1:LM)
    REAL(KIND=r8)     :: pblfq     (1:IM,1:LM)
    !REAL(KIND=r8)    :: rtmp      (1:IM,1:LM+1)
    REAL(KIND=r8)     :: k_m_troen (1:IM,1:LM)
    REAL(KIND=r8)     :: k_t_troen (1:IM,1:LM)
    REAL(KIND=r8)     :: k_rad     (1:IM,1:LM)
    LOGICAL           :: TEST      (1:IM)
    LOGICAL           :: CONDIC    (1:IM)

    !-----------------------------------------------------------------------
    !
    !      initialize variables
    k_m_entr    = 0.0_r8
    k_t_entr    = 0.0_r8
    k_rad_diag  = 0.0_r8
    k_sfc_diag  = 0.0_r8
    zsml        = 0.0_r8 
    zradml      = 0.0_r8 
    zcloud      = 0.0_r8 
    zradbase    = 0.0_r8 
    
    wentr_rad_diag   = 0.0_r8
    wentr_sfc_diag   = 0.0_r8
    del_buoy_diag    = 0.0_r8
    vrad_diag        = 0.0_r8
    kentrad_diag     = 0.0_r8
    vbrv_diag        = 0.0_r8
    wentr_brv_diag   = 0.0_r8
    dsiems_diag      = 0.0_r8
    chis_diag        = 0.0_r8
    zcldtop_diag     = 0.0_r8
    delsinv_diag     = 0.0_r8
    slmixture_diag   = 0.0_r8
    cldradf_diag     = 0.0_r8
    vsfc_diag        = 0.0_r8
    radrcode         = 0.0_r8

    convpbl    = 0.0_r8
    wentr_pbl  = 0.0_r8
    vsurf      = 0.0_r8
    vshear     = 0.0_r8
    pblfq      = 0.0_r8
    k_t_troen  = 0.0_r8
    k_m_troen  = 0.0_r8
    radpbl     = 0.0_r8
    svpcp      = 0.0_r8
    zradbase   = 0.0_r8
    zradtop    = 0.0_r8
    zradml     = 0.0_r8
    zsml       = 0.0_r8
    zcloud     = 0.0_r8
    vrad       = 0.0_r8
    radf       = 0.0_r8
    radfq      = 0.0_r8
    wentr_rad  = 0.0_r8
    k_rad      = 0.0_r8
    k_t_entr   = 0.0_r8
    k_m_entr   = 0.0_r8
    fqinv      = 0.0_r8
    zinv       = 0.0_r8
    invstr     = 0.0_r8
    zsml       = 0.0_r8   ! note that this must be zero as this is 
    ! indicates stable surface layer and this
    ! value is output for use in gravity
    ! wave drag scheme
    qlcrit     = 1.0e-6_r8
    Abuoy      = 0.23_r8
    Ashear     = 25.0_r8
    wentrmax   = 0.05_r8
    TEST=.TRUE.
    CONDIC=.TRUE.


    !--------------------------------------------------------------------------
    ! Initialize optional outputs

    !IF (ASSOCIATED(wentr_sfc_diag)) wentr_sfc_diag = missing_value
    !IF (ASSOCIATED(wentr_rad_diag)) wentr_rad_diag = missing_value
    !IF (ASSOCIATED(del_buoy_diag))  del_buoy_diag  = missing_value
    !IF (ASSOCIATED(vrad_diag))      vrad_diag      = missing_value
    !IF (ASSOCIATED(vsfc_diag))      vsfc_diag      = missing_value
    !IF (ASSOCIATED(kentrad_diag))   kentrad_diag   = missing_value  
    !IF (ASSOCIATED(chis_diag))      chis_diag      = missing_value
    !IF (ASSOCIATED(vbrv_diag))      vbrv_diag      = missing_value
    !IF (ASSOCIATED(dsiems_diag))    dsiems_diag    = missing_value
    !IF (ASSOCIATED(wentr_brv_diag)) wentr_brv_diag = missing_value
    !IF (ASSOCIATED(zcldtop_diag))   zcldtop_diag   = missing_value
    !IF (ASSOCIATED(g_diag)) slmixture_diag = missing_value
    !IF (ASSOCIATED(delsinv_diag))   delsinv_diag   = missing_value
    !IF (ASSOCIATED(cldradf_diag))   cldradf_diag   = missing_value
    !IF (ASSOCIATED(radrcode))       radrcode       = missing_value
    wentr_sfc_diag = missing_value
    wentr_rad_diag = missing_value
    del_buoy_diag  = missing_value
    vrad_diag           = missing_value
    vsfc_diag           = missing_value
    kentrad_diag   = missing_value
    chis_diag           = missing_value
    vbrv_diag           = missing_value
    dsiems_diag    = missing_value
    wentr_brv_diag = missing_value
    zcldtop_diag   = missing_value
    slmixture_diag = missing_value
    delsinv_diag   = missing_value
    cldradf_diag   = missing_value
    radrcode           = missing_value
    !-----------------------------------------------------------------------
    !
    !      Sizes
    !

    nlev = LM
    nlon = IM   

    !-----------------------------------------------------------------------
    !
    !      set up specific humidities and static energies  
    !      compute airdensity
    !
    DO k=1,LM
       DO i=1,IM
          IF( T(i,k) <= TFREEZE-ramp) THEN
             HLEFF(i,k) = HLS
          END IF
       END DO
   END DO
    
   DO k=1,LM
      DO i=1,IM
            IF ((T(i,k) > TFREEZE-ramp) .AND. (T(i,k) < TFREEZE) ) THEN
                HLEFF(i,k) =  ( (T(i,k) - TFREEZE+ramp)*HLV + (TFREEZE -T(i,k) )*HLS   ) / ramp
            END IF
      END DO
   END DO     

   DO k=1,LM
      DO i=1,IM
            IF( T(i,k) >= TFREEZE  ) THEN
                HLEFF(i,k) = HLV
            END IF
      END DO
   END DO         

   !--------------------------------------------------------------------------
   !      Compute total cloud condensate - qc. These are grid box mean values.
   DO k=1,LM
      DO i=1,IM
            qc (i,k) = (ql(i,k) + qi(i,k))
            slv(i,k) = cp_air * t(i,k) *(1+d608*qv(i,k) - qc(i,k)) + grav*zfull(i,k) - hleff(i,k)*qc(i,k)
      END DO
   END DO         


    !!       slv     = cp_air*t + grav*zfull - hleff*qc
    !!       slv     = slv*(1+d608*(qv+qc))
   DO k=1,LM
      DO i=1,IM
             density (i,k)= pfull(i,k)/(RDGAS*(t(i,k) *(1.0_r8+d608*qv(i,k)-qc(i,k))))               
      END DO
   END DO         


    !--------------------------
    ! 
    !      big loop over points
    !
    ibot = nlev

    DO i=1,IM

       !---------------
       ! reset indices

       ipbl   (i) = -1
       kcldtop(i) = -1
    END DO
    

    !------------------------------------------------------
    ! Find depth of surface driven mixing by raising a 
    ! parcel from the surface with some excess buoyancy
    ! to its level of neutral buoyancy.  Note the use
    ! slv as the density variable permits one to goes
    ! through phase changes to find parcel top

    CALL mpbl_depth(&
             IM                        , &
         ibot                        , &
             tpfac_sfc                , &
             entrate_sfc                , &
             pceff_sfc                , & 
             t           (1:IM,1:ibot), &
             qv           (1:IM,1:ibot), &
             u           (1:IM,1:ibot), &
             v           (1:IM,1:ibot), &
             zfull     (1:IM,1:ibot), &
             pfull     (1:IM,1:ibot), &
             b_star    (1:IM)        , &
             u_star    (1:IM)        , &
             ipbl           (1:IM)        , &
         zsml           (1:IM)          )

       DO i=1,IM
          !-----------------------------------------------------------
          !
          ! SURFACE DRIVEN CONVECTIVE LAYERS
          !
          ! Note this part is done only if b_star > 0., that is,
          ! upward surface buoyancy flux


          IF (b_star(i) .GT. 0.0_r8) THEN


             !------------------------------------------------------
             ! Define velocity scales vsurf and vshear
             !           
             ! vsurf   =  (u_star*b_star*zsml)**(1/3)
             ! vshear  =  (Ashear**(1/3))*u_star

             vsurf3 (i)  = u_star(i)*b_star(i)*zsml(i)
             vshear3(i)  = Ashear*u_star(i)*u_star(i)*u_star(i)

             vsurf  (i) = vsurf3(i)**(1.0_r8/3.0_r8)
             !IF (ASSOCIATED(vsfc_diag))      vsfc_diag(i)      = vsurf(i)
             vsfc_diag(i)      = vsurf(i)
             !------------------------------------------------------
             ! Following Lock et al. 2000, limit height of surface
             ! well mixed layer if interior stable interface is
             ! found.  An interior stable interface is diagnosed if
             ! the slope between 2 full levels is greater than critjump

             critjump(i) = 2.0_r8
           END IF
      END DO
      
           

      TEST=.TRUE.
      DO k = ibot, MINVAL(ipbl)+1, -1
         DO i=1,IM
            !-----------------------------------------------------------
            !
            ! SURFACE DRIVEN CONVECTIVE LAYERS
            !
            ! Note this part is done only if b_star > 0., that is,
            ! upward surface buoyancy flux
            IF (b_star(i) .GT. 0.0_r8) THEN
               IF (TEST(i) .AND. ipbl(i) .LT. ibot) THEN 

                  tmpjump =(slv(i,k-1)-slv(i,k))/cp_air 
             
                  IF (tmpjump .GT. critjump(i)) THEN
                     ipbl(i) = k
                     zsml(i) = zhalf(i,ipbl(i))
                     TEST(i) = .FALSE.
                  END IF
               END IF
            END IF                
        END DO
      END DO
      

      DO i=1,IM

          !-----------------------------------------------------------
          !
          ! SURFACE DRIVEN CONVECTIVE LAYERS
          !
          ! Note this part is done only if b_star > 0., that is,
          ! upward surface buoyancy flux


          IF (b_star(i) .GT. 0.0_r8) THEN

             !-------------------------------------
             ! compute entrainment rate
             !

             tmp1(i) = grav*MAX(0.1_r8,(slv(i,ipbl(i)-1)-slv(i,ipbl(i)))/   &
                         cp_air)/(slv(i,ipbl(i))/cp_air)
             tmp2(i) = ((vsurf3(i)+vshear3(i))**(2.0_r8/3.0_r8)) / zsml(i)

             wentr_tmp(i)= MIN( wentrmax,  MAX(0.0_r8, (beta_surf *        &
                             (vsurf3(i) + vshear3(i))/zsml(i))/         &
                             (tmp1(i)+tmp2(i)) ) )

             !----------------------------------------
             ! fudgey adjustment of entrainment to reduce it
             ! for shallow boundary layers, and increase for 
             ! deep ones
             IF ( zsml(i) .LT. 1600.0_r8 ) THEN 
                wentr_tmp(i) = wentr_tmp(i) * ( zsml(i) / 800.0_r8 )
             ELSE
                wentr_tmp(i) = 2.0_r8*wentr_tmp(i)
             ENDIF
             !-----------------------------------------


             k_entr_tmp(i) = wentr_tmp(i)*(zfull(i,ipbl(i)-1) - zfull(i,ipbl(i)))  
             k_entr_tmp(i) = MIN ( k_entr_tmp(i), akmax )

             pblfq(i,ipbl(i):ibot) = 1.0_r8
             convpbl(i)         = 1.0_r8
             wentr_pbl(i)       = wentr_tmp(i)
             k_t_troen(i,ipbl(i))  = k_entr_tmp(i)
             k_m_troen(i,ipbl(i))  = k_entr_tmp(i)
             k_t_entr (i,ipbl(i))  = k_t_entr(i,ipbl(i)) + k_entr_tmp(i)
             k_m_entr (i,ipbl(i))  = k_m_entr(i,ipbl(i)) + k_entr_tmp(i)

             !IF (ASSOCIATED(wentr_sfc_diag)) wentr_sfc_diag(i) = wentr_tmp(i)
             wentr_sfc_diag(i) = wentr_tmp(i)
             !------------------------------------------------------
             ! compute diffusion coefficients in the interior of
             ! the PBL
          END IF
       END DO             


       CALL diffusivity_pbl2(&
            ibot                    , &
            IM                      , &
            LM                      , &
            b_star     (1:IM)     , &
            ipbl       (1:IM)     , &
            zsml       (1:IM)     , &
            khsfcfac                    , &
            k_entr_tmp (1:IM)     , & 
            vsurf      (1:IM)     , &
            frland     (1:IM)     , & 
            zhalf      (1:IM,0:LM), &
            k_m_troen  (1:IM,1:LM), &
            k_t_troen  (1:IM,1:LM)   )
            
            
      DO i=1,IM
          !-----------------------------------------------------------
          !
          ! SURFACE DRIVEN CONVECTIVE LAYERS
          !
          ! Note this part is done only if b_star > 0., that is,
          ! upward surface buoyancy flux
          IF (b_star(i) .GT. 0.0_r8) THEN
             IF (ipbl(i) .LT. ibot) THEN


                k_t_entr(i,(ipbl(i)+1):ibot) =                    & 
                     k_t_entr(i,(ipbl(i)+1):ibot) +               &
                     k_t_troen(i,(ipbl(i)+1):ibot)

                k_m_entr(i,(ipbl(i)+1):ibot) =                    & 
                     k_m_entr(i,(ipbl(i)+1):ibot) +               &
                     k_m_troen(i,(ipbl(i)+1):ibot)*prandtlsfc
             END IF
          END IF
       END DO

       !-----------------------------------------------------------
       !
       ! NEGATIVELY BUOYANT PLUMES DRIVEN BY 
       ! LW RADIATIVE COOLING AND/OR BUOYANCY REVERSAL
       !
       ! This part is done only if a level kcldtop can be 
       ! found with: 
       !
       !    qc(kcldtop)>=qlcrit.and.qc(kcldtop-1)<qlcrit
       !
       ! below zcldtopmax

       kmax = ibot+1
       TEST =.TRUE.

       DO k = 1, ibot
          DO i=1,IM
             IF(TEST(i) .AND. zhalf(i,k) < zcldtopmax) THEN
                kmax(i) = k
                TEST(i)=.FALSE.
             END IF
          END DO
       END DO

          !-----------------------------------------------------------
          ! Find cloud top and bottom using GRID BOX MEAN or IN-CLOUD 
          ! value of qc.  Decision occurs where qc is calculated
          DO i=1,IM
             kcldtop(i)  = ibot+1
          END DO
          
          TEST =.TRUE.
          DO k = ibot,MINVAL(kmax),-1
             DO i=1,IM
                qlm1 = qc(i,k-1)  ! qc one level UP
                ql00 = qc(i,k)
                stab = slv(i,k-1) - slv(i,k) 
                IF (TEST(i) .AND. ( ql00  .GE. qlcrit ) .AND. ( qlm1 .LT. qlcrit) .AND. (stab .GT. 0.0_r8) ) THEN
                   kcldtop(i)  = k   
                   TEST(i)=.FALSE.
                END IF
             END DO
          END DO


          DO i=1,IM
             IF (kcldtop(i) .GE. ibot+1) THEN 
             !IF (ASSOCIATED(radrcode)) radrcode(i)=1.
                 radrcode(i)=1.0_r8
                !go to 55
                CONDIC(i)=.FALSE.
             ENDIF
             IF(CONDIC(i)) THEN 
                kcldtop2(i)=MIN( kcldtop(i)+1,nlev)
                ! Look one level further down in case first guess is a thin diffusive veil
                IF( (qc(i,kcldtop(i)) .LT. 10.0_r8*qlcrit ) .AND. (qc(i,kcldtop2(i)) .GE. 10.0_r8*qc(i,kcldtop(i)) ) ) THEN
                   kcldtop(i)=kcldtop2(i)
                ENDIF
             END IF
          END DO

          kcldbot  = ibot+1
                 TEST  =.TRUE.
          DO k = ibot,MINVAL(kcldtop(:)),-1
             DO i=1,IM
                IF(CONDIC(i)) THEN 
                   qlm1 = qc(i,k-1)  ! qc one level UP
                   ql00 = qc(i,k)
                   IF (TEST(i).AND. ( ql00  .LT. qlcrit ) .AND. ( qlm1 .GE. qlcrit) ) THEN
                      kcldbot(i)  = k   
                      TEST(i) =.fALSE.
                   END IF
                END IF
             END DO
          END DO

          DO i=1,IM
             IF (CONDIC(i).AND.kcldtop(i) .EQ. kcldbot(i)) THEN 
                !IF (ASSOCIATED(radrcode)) radrcode(i)=2. 
                radrcode(i)=2.0_r8
                CONDIC(i) =.FALSE.
             END IF
          END DO

          DO i=1,IM
             IF(CONDIC(i)) THEN 

                ! With diffusion of ql, qi "cloud top" found via these quantities may be above radiation max
                kcldtop2(i)=MIN( kcldtop(i)+2,nlev)
                maxradf(i) = MAXVAL( -1.0_r8*tdtlw_in(i,kcldtop2(i):kcldtop2(i)) )

                maxradf(i) = maxradf(i)*cp_air*     &
                             ( (phalf(i,kcldtop(i)+1)-phalf(i,kcldtop(i))) / grav )

                maxradf(i) = MAX( maxradf(i) , 0.0_r8 ) ! do not consider cloud tops that are heating

                !-----------------------------------------------------------
                ! Calculate optimal mixing fraction - chis - for buoyancy 
                ! reversal.  Use effective heat of evap/subl *tion.  Ignore 
                ! diffs across cldtop
                hlf = hleff(i,kcldtop(i))

                tmp1(i) = ( slv(i,kcldtop(i)-1)  -  hlf*qc(i,kcldtop(i)-1) ) - &
                         ( slv(i,kcldtop(i))    -  hlf*qc(i,kcldtop(i))   )
                tmp1(i) = dqs(i,kcldtop(i))*tmp1(i)/cp_air

                tmp2(i) = ( qv(i,kcldtop(i)-1)   +  qc(i,kcldtop(i)-1) ) - &
                           ( qv(i,kcldtop(i))     +  qc(i,kcldtop(i))   )  

                chis(i) = -qc(i,kcldtop(i))*( 1 + hlf * dqs(i,kcldtop(i)) / cp_air )

                IF ( ( tmp2(i) - tmp1(i) ) >= 0.0_r8 ) THEN
                   chis(i) = 0.0_r8
                ELSE
                   chis(i) = chis(i) / ( tmp2(i) - tmp1(i) ) 
                ENDIF

                IF ( chis(i) .GT. 1.0_r8 ) chis(i)=1.0_r8

                slmixture(i) = ( 1.0_r8-chis(i) )* ( slv(i,kcldtop(i))    -  hlf*qc(i,kcldtop(i))   )   &
                               +       chis(i)  * ( slv(i,kcldtop(i)-1)  -  hlf*qc(i,kcldtop(i)-1) )


                !-----------------------------------------------------------
                ! compute temperature of parcel at cloud top, svpcp.
                svpcp(i) = slmixture(i) /cp_air

                buoypert(i)   = ( slmixture(i) - slv(i,kcldtop(i)) )/cp_air
   
                !-----------------------------------------------------------
                ! calculate my best guess at the LCs' D parameter attributed 
                ! to Siems et al.
                stab       = slv(i,kcldtop(i)-1) - slv(i,kcldtop(i))
                IF (stab .EQ. 0.0_r8) THEN 
                  dsiems(i)     =  ( slv(i,kcldtop(i)) - slmixture(i) ) ! / 1.0_r8  ! arbitrary, needs to be re-thought 
                ELSE
                  dsiems(i)     =  ( slv(i,kcldtop(i)) - slmixture(i) ) / stab
                END IF
                dsiems (i)    =  MIN( dsiems(i), 10.0_r8 )
                dsiems (i)    =  MAX( dsiems(i),  0.0_r8 )
                radf(i)  = maxradf(i)
                zradtop(i) = zhalf(i,kcldtop(i))

                !-----------------------------------------------------------
                ! find depth of radiatively driven convection 

                !-----------------------------------------------------------
                ! Expose radperturb and other funny business outside of radml_depth
                radperturb (i)  = MIN( maxradf(i)/100.0_r8 , 0.3_r8 ) ! dim. argument based on 100m deep cloud over 1000s
                do_jump_exit(i) = .TRUE.
                critjump(i)     = 0.3_r8
                svpcp(i)        = svpcp(i) - radperturb(i)
             END IF
          END DO

          CALL radml_depth( &
                              IM                            , &
                              LM                            , &
                              ibot                          , &
                              kcldtop     (1:IM)            , &
                              svpcp       (1:IM)            , &
                              zradtop     (1:IM)            , &
                              critjump    (1:IM)            , &
                              do_jump_exit(1:IM)            , &
                              slv         (1:IM,1:LM)/cp_air, &
                              zfull       (1:IM,1:LM)       , &
                              zhalf       (1:IM,0:LM)       , &
                              zradbase    (1:IM)            , &
                              zradml      (1:IM)            , &
                              CONDIC      (1:IM)              )       
                    
          DO i=1,IM
             IF(CONDIC(i))THEN          
                IF (kcldtop(i) >= ibot) THEN 
                   zradbase(i) = 0.0_r8
                   zradml(i)   = zradtop(i)  
                END IF

                zcloud(i) = zhalf(i,kcldtop(i)) - zhalf(i,kcldbot(i))

                IF (zradml(i) .LE. 0.0_r8 ) THEN 
                   !IF (ASSOCIATED(radrcode)) radrcode(i)=3.
                   radrcode(i)=3.0_r8
                   CONDIC(i) =.FALSE.
                   !go to 55   ! break out here if zradml<=0.0
                END IF
                IF(CONDIC(i))THEN          

                   !-----------------------------------------------------------
                   ! compute radiation driven scale
                   !
                   ! Vrad**3 = g*zradml*radf/density/slv

                    vrad3(i) = grav*zradml(i)*maxradf(i)/density(i,kcldtop(i))/ &
                          slv(i,kcldtop(i))   


                   !-----------------------------------------------------------
                   ! compute entrainment rate
                   !

                   !-----------------------------------------------------------
                   ! tmp1 here should be the buoyancy jump at cloud top
                   ! SAK has it w/ resp to parcel property - svpcp. Im not 
                   ! sure about that.
                   tmp1(i) = grav*MAX(0.1_r8,((slv(i,kcldtop(i)-1)/cp_air)-          &
                             svpcp(i)))/(slv(i,kcldtop(i))/cp_air)

                   !-----------------------------------------------------------
                   ! Straightforward buoyancy jump across cloud top
                    tmp1(i) = grav*MAX( 0.1_r8, ( slv(i,kcldtop(i)-1)-slv(i,kcldtop(i)) )/cp_air )          &
                         / ( slv(i,kcldtop(i)) /cp_air )

                   !-----------------------------------------------------------
                   ! compute buoy rev driven scale
                    vbr3(i)  = ( MAX( tmp1(i)*zcloud(i), 0.0_r8)**3 )
                    vbr3(i)  = Abuoy*(chis(i)**2)*MAX(dsiems(i),0.0_r8)*SQRT( vbr3(i) ) 

                   !----------------------------------------
                   ! adjust velocity scales to prevent jumps 
                   ! near zradtop=zcldtopmax
                   IF ( zradtop(i) .GT. zcldtopmax-500.0_r8 ) THEN 
                      vrad3(i) = vrad3(i)*(zcldtopmax - zradtop(i))/500.0_r8  
                      vbr3(i)  = vbr3(i) *(zcldtopmax - zradtop(i))/500.0_r8  
                   ENDIF

                   vrad3(i)=MAX( vrad3(i), 0.0_r8 ) ! these REAL(KIND=r8)ly should not be needed
                   vbr3(i) =MAX( vbr3 (i),  0.0_r8 )
                    !-----------------------------------------



                   vrad(i) = vrad3(i) ** (1.0_r8/3._r8)    
                   vbrv(i) = vbr3(i)**(1.0_r8/3.0_r8)

                   tmp2(i) = (  vrad(i)**2 + vbrv(i)**2  ) / MAX(zradml(i),1.0_r8)
                   wentr_rad(i) = MIN(wentrmax,beta_rad*(vrad3(i)+vbr3(i))/zradml(i)/  &
                                (tmp1(i)+tmp2(i)))

                   wentr_brv(i) =  beta_rad*vbr3(i)/zradml(i)/(tmp1(i)+tmp2(i))


                   !----------------------------------------
                   ! fudgey adjustment of entrainment to reduce it
                   ! for shallow boundary layers, and increase for 
                   ! deep ones
                   IF ( zradtop(i) .LT. 2400.0_r8 ) THEN 
                      wentr_rad(i) = wentr_rad(i) * ( zradtop(i) / 800.0_r8 )
                   ELSE
                      wentr_rad(i) = 3.0_r8*wentr_rad(i)
                   END IF
                   !-----------------------------------------

                   k_entr_tmp(i) = MIN ( akmax, wentr_rad(i)* (zfull(i,kcldtop(i)-1)-zfull(i,kcldtop(i))) )

                   radfq(i,kcldtop(i))     = 1.0_r8
                   radpbl(i)            = 1.0_r8
                   k_rad(i,kcldtop(i))     = k_entr_tmp(i)
                   k_t_entr (i,kcldtop(i)) = k_t_entr(i,kcldtop(i)) + k_entr_tmp(i)
                   k_m_entr (i,kcldtop(i)) = k_m_entr(i,kcldtop(i)) + k_entr_tmp(i)


             !IF (ASSOCIATED(del_buoy_diag))  del_buoy_diag(i)  = tmp1
             !IF (ASSOCIATED(vrad_diag))      vrad_diag(i)      = vrad(i)
             !IF (ASSOCIATED(kentrad_diag))   kentrad_diag(i,j)   = k_entr_tmp(i,j)!

             !IF (ASSOCIATED(chis_diag))      chis_diag(i,j)      = chis(i,j)
             !IF (ASSOCIATED(vbrv_diag))      vbrv_diag(i,j)      = vbrv(i,j)
             !IF (ASSOCIATED(dsiems_diag))    dsiems_diag(i,j)    = dsiems
             !IF (ASSOCIATED(wentr_brv_diag)) wentr_brv_diag(i,j) = wentr_brv(i,j)
             !IF (ASSOCIATED(wentr_rad_diag)) wentr_rad_diag(i,j) = wentr_rad(i,j)

             !IF (ASSOCIATED(zcldtop_diag))   zcldtop_diag(i,j)   = zhalf(i,j,kcldtop) 
             !IF (ASSOCIATED(slmixture_diag)) slmixture_diag(i,j) = buoypert
             !IF (ASSOCIATED(delsinv_diag))   delsinv_diag(i,j)   = ( slv(i,j,kcldtop-1) - slv(i,j,kcldtop) )/cp_air
             !IF (ASSOCIATED(cldradf_diag))   cldradf_diag(i,j)   = radf(i,j)

                   del_buoy_diag(i)  = tmp1(i)
                     vrad_diag(i)      = vrad(i)
                   kentrad_diag(i)   = k_entr_tmp(i)

                   chis_diag(i)      = chis(i)
                   vbrv_diag(i)      = vbrv(i)
                   dsiems_diag(i)    = dsiems(i)
                   wentr_brv_diag(i) = wentr_brv(i)
                   wentr_rad_diag(i) = wentr_rad(i)

                   zcldtop_diag(i)   = zhalf(i,kcldtop(i)) 
                   slmixture_diag(i) = buoypert(i)
                   delsinv_diag(i)   = ( slv(i,kcldtop(i)-1) - slv(i,kcldtop(i)) )/cp_air
                   cldradf_diag(i)   = radf(i)

                   !-----------------------------------------------------------
                   ! handle case of radiatively driven top being the same top
                   ! as surface driven top

                   IF (ipbl(i) .EQ. kcldtop(i) .AND. ipbl(i) .GT. 0) THEN

                      tmp2(i) = ((vbr3(i)+vrad3(i)+vsurf3(i)+vshear3(i))**(2.0_r8/3.0_r8)) / zradml(i)
   
                      wentr_rad(i) = MIN( wentrmax,  MAX(0.0_r8,              &
                        ((beta_surf *(vsurf3(i) + vshear3(i))+beta_rad*(vrad3(i)+vbr3(i)) )/ &
                         zradml(i))/(tmp1(i)+tmp2(i)) ) )

                      wentr_pbl(i)       = wentr_rad(i)

                      k_entr_tmp(i) = MIN ( akmax, wentr_rad(i)*(zfull(i,kcldtop(i)-1)-zfull(i,kcldtop(i))) )

                      pblfq(i,ipbl(i))        = 1.0_r8
                      radfq(i,kcldtop(i))     = 1.0_r8
                      radpbl(i)            = 1.0_r8
                      k_rad(i,kcldtop(i))     = k_entr_tmp(i)
                      k_t_troen(i,ipbl(i))    = k_entr_tmp(i)
                      k_m_troen(i,ipbl(i))    = k_entr_tmp(i)
                      k_t_entr (i,kcldtop(i)) = k_entr_tmp(i)
                      k_m_entr (i,kcldtop(i)) = k_entr_tmp(i)
                   END IF
                END IF
             END IF
          END DO 
             !-----------------------------------------------------------
             ! if there are any interior layers to calculate diffusivity
 
             DO k = MINVAL(kcldtop(:))+1,ibot
                DO i=1,IM
                   IF(CONDIC(i))THEN          
                      IF ( kcldtop(i) .LT. ibot ) THEN   
 
                         ztmp = MAX(0.0_r8,(zhalf(i,k)-zradbase(i))/zradml(i) )

                         IF (ztmp.GT.0.0_r8) THEN

                            radfq(i,k)     = 1.0_r8
                            k_entr_tmp(i) = khradfac*vonkarm*( vrad(i)+vbrv(i) )*ztmp* zradml(i)*ztmp*((1.0_r8-ztmp)**0.5_r8)
                            k_entr_tmp(i)       = MIN ( k_entr_tmp(i), akmax )
                            k_rad    (i,k) = k_entr_tmp(i)
                            k_t_entr (i,k) = k_t_entr(i,k) + k_entr_tmp(i)
                            k_m_entr (i,k) = k_m_entr(i,k) + k_entr_tmp(i)*prandtlrad
                         END IF
                      END IF
                   END IF
                END DO
             END DO   

          !-----------------------------------------------------------
          ! handle special case of zradbase < zsml
          !
          ! in this case there should be no entrainment from the 
          ! surface.
          DO i=1,IM
             IF (CONDIC(i) .and. zradbase(i) .LT. zsml(i) .AND. convpbl(i) .EQ. 1.0_r8 .AND. ipbl(i) .GT. kcldtop(i)) THEN
                wentr_pbl(i)           = 0.0_r8
                pblfq    (i,ipbl(i)) = 0.0_r8
                k_t_entr (i,ipbl(i)) = k_t_entr (i,ipbl(i)) - k_t_troen(i,ipbl(i))
                k_m_entr (i,ipbl(i)) = k_m_entr (i,ipbl(i)) - k_m_troen(i,ipbl(i))          
                k_t_troen(i,ipbl(i)) = 0.0_r8
                k_m_troen(i,ipbl(i)) = 0.0_r8 
             END IF
          END DO



          k_sfc_diag = k_t_troen    
          k_rad_diag = k_rad    

          !-----------------------------------------------------------
          !
          ! Modify diffusivity coefficients using MAX( A , B )        
 
          DO k = 2, ibot     
              DO i=1,IM
                 k_t_entr(i,k)=MIN(300.0_r8,MAX( 1.00_r8,k_t_entr(i,k))) !(m/sec**2)
                 k_m_entr(i,k)=MIN(300.0_r8,MAX( 0.10_r8,k_m_entr(i,k))) !(m/sec**2)
                 diff_t(i,k) = MAX( k_t_entr(i,k) ,  diff_t(i,k) )
                 diff_m(i,k) = MAX( k_m_entr(i,k) ,  diff_m(i,k) )
                 k_t_entr(i,k) = MAX( k_t_entr(i,k) ,  0.0_r8 )
                 k_m_entr(i,k) = MAX( k_m_entr(i,k) ,  0.0_r8 )
             END DO
         END DO
    !END DO



    !-----------------------------------------------------------------------
    ! 
    !      subroutine end
    !

  END SUBROUTINE entrain


  !======================================================================= 
  !
  !  Subroutine to calculate pbl depth
  !

  !======================================================================= 
  !
  !  Subroutine to calculate pbl depth
  !
  SUBROUTINE mpbl_depth(&
       im          , &! IM                        , &
       ibot        , &! ibot                        , &
       tpfac       , &! tpfac_sfc                , &
       entrate     , &! entrate_sfc                , &
       pceff       , &! pceff_sfc                , & 
       t           , &! t         (1:IM,1:ibot), &
       q           , &! qv         (1:IM,1:ibot), &
       u           , &! u         (1:IM,1:ibot), &
       v           , &! v         (1:IM,1:ibot), &
       z           , &! zfull         (1:IM,1:ibot), &
       p           , &! pfull         (1:IM,1:ibot), &
       b_star      , &! b_star         (1:IM)       , &
       u_star      , &! u_star         (1:IM)       , &
       ipbl        , &! ipbl         (1:IM)       , &
       ztop          )! zsml         (1:IM)         )

    !
    !  -----
    !  INPUT
    !  -----
    !
    !  t             temperature (K)
    !  q             specific humidity (g/g)
    !  u             zonal wind (m/s)
    !  v             meridional wind (m/s)
    !  b_star        buoyancy scale (m s-2)
    !  u_star        surface velocity scale (m/s)
    !       
    !  ------
    !  OUTPUT
    !  ------
    !
    !  ipbl          half level containing pbl height
    !  ztop          pbl height (m)
    INTEGER, INTENT(IN   ) :: im          ! IM                      , &
    INTEGER, INTENT(IN   ) :: ibot        ! ibot                     , &
    REAL(KIND=r8),    INTENT(in   ) :: t     (1:IM,1:ibot)
    REAL(KIND=r8),    INTENT(in   ) :: z     (1:IM,1:ibot)
    REAL(KIND=r8),    INTENT(in   ) :: q     (1:IM,1:ibot)
    REAL(KIND=r8),    INTENT(in   ) :: p     (1:IM,1:ibot)
    REAL(KIND=r8),    INTENT(in   ) :: u     (1:IM,1:ibot)
    REAL(KIND=r8),    INTENT(in   ) :: v     (1:IM,1:ibot)
    REAL(KIND=r8),    INTENT(in   ) :: b_star(1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: u_star(1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: tpfac
    REAL(KIND=r8),    INTENT(in   ) :: entrate 
    REAL(KIND=r8),    INTENT(in   ) :: pceff
    INTEGER      ,    INTENT(out  ) :: ipbl  (1:IM)
    REAL(KIND=r8),    INTENT(out  ) :: ztop  (1:IM)


    REAL(KIND=r8)     :: tep  (1:IM)
    REAL(KIND=r8)     :: z1   (1:IM)
    REAL(KIND=r8)     :: z2   (1:IM)
    REAL(KIND=r8)     :: t1   (1:IM)
    REAL(KIND=r8)     :: t2   (1:IM)
    REAL(KIND=r8)     :: qp   (1:IM)
    REAL(KIND=r8)     :: pp   (1:IM)
    REAL(KIND=r8)     :: qsp  (1:IM) 
    REAL(KIND=r8)     :: dqp  (1:IM) 
    REAL(KIND=r8)     :: dqsp (1:IM)
    REAL(KIND=r8)     :: u1   (1:IM)
    REAL(KIND=r8)     :: v1   (1:IM)
    REAL(KIND=r8)     :: du   (1:IM)
    REAL(KIND=r8)     :: ws
    REAL(KIND=r8)     :: k_t_ref
    REAL(KIND=r8)     :: entfr     (1:IM)
    REAL(KIND=r8)     :: tpfac_x 
    REAL(KIND=r8)     :: entrate_x (1:IM)
    REAL(KIND=r8)     :: vscale    (1:IM)
    INTEGER  :: k,nlev,i
    LOGICAL  :: TEST      (1:IM)

    REAL(KIND=r8)     :: qst( 1:IM )

    ipbl =0;    ztop =0.0_r8;    tep  =0.0_r8;
    z1   =0.0_r8;    z2   =0.0_r8;    t1   =0.0_r8;    t2   =0.0_r8;
    qp   =0.0_r8;    pp   =0.0_r8;    qsp  =0.0_r8;    dqp  =0.0_r8;
    dqsp =0.0_r8;    u1   =0.0_r8;    v1   =0.0_r8;    du   =0.0_r8;
    ws=0.0_r8;    k_t_ref=0.0_r8;    entfr=0.0_r8;    tpfac_x =0.0_r8;
    entrate_x=0.0_r8;    vscale  =0.0_r8;  
    qst=0.0_r8;

          !-----------------------------------------------------------
          !
          ! SURFACE DRIVEN CONVECTIVE LAYERS
          !
          ! Note this part is done only if b_star > 0., that is,
          ! upward surface buoyancy flux
    nlev = ibot
    TEST  =.TRUE.
    DO i=1,IM
       IF (b_star(i) .GT. 0.0_r8) THEN

  
          ztop(i) = 0._r8

          !calculate surface parcel properties

         tep(i)  = t(i,nlev)
         qp (i)  = q(i,nlev)
         z1 (i)  = z(i,nlev)

         !--------------------------------------------
         ! wind dependence of plume character. 
         ! 
         !    actual_entrainment_rate_at_z  ~ entrate * [ 1.0 +  |U(z)-U(0)| / vscale ]
         ! 
         ! entrate:  tunable param from rc file
         ! vscale:   tunable param hardwired here.

         vscale(i)    = 5.0_r8 ! m s-1


        !---------------------------------------------
        ! tpfac scales up bstar by inv. ratio of
        ! heat-bubble area to stagnant area

        tep(i)  = tep(i) * (1.0_r8+ tpfac * b_star(i)/grav)

        !search for level where this is exceeded              

        ztop(i) = z1(i)
        t1  (i) = t(i,nlev)
        v1  (i) = v(i,nlev)
        u1  (i) = u(i,nlev)
       END IF
    END DO
        
    DO k = nlev-1 , 2, -1
       DO i=1,IM
          IF (b_star(i) .GT. 0.0_r8) THEN

             z2(i) = z(i,k)
             t2(i) = t(i,k)
             pp(i) = p(i,k)

             du(i) = SQRT ( ( u(i,k) - u1(i) )**2 + ( v(i,k) - v1(i) )**2 )

             entrate_x(i) = entrate * ( 1.0_r8 + du(i) / vscale(i) )

             entfr(i) = MIN( entrate_x(i) *(z2(i)-z1(i)), 0.99_r8 )

             qp(i)  = qp(i)  + entfr(i)*(q(i,k)-qp(i))

             ! dry adiabatic ascent through one layer.
             ! Static energy conserved. 
             tep(i) = tep(i) - grav*( z2(i)-z1(i) )/cp_air

             ! Environmental air entrained
             tep(i) = tep(i) + entfr(i)*(t(i,k)-tep(i))

             dqsp(i) = geos_dqsat(tep(i) , pp(i) , qsat=qsp(i),  pascals=.TRUE. )

             dqp(i) = MAX( qp(i) - qsp(i), 0.0_r8 )/(1.0_r8+(hlv/cp_air)*dqsp(i) )
             qp(i)  = qp(i) - dqp(i)
             tep(i) = tep(i)  + pceff * hlv * dqp(i)/cp_air  ! "Precipitation efficiency" basically means fraction
             ! of condensation heating that gets applied to parcel


             ! If parcel temperature (tep) colder than env (t2)
             ! OR if entrainment too big, declare this the PBL top
             IF ( TEST(i) .AND. (t2(i) .GE. tep(i)) .OR. ( entfr(i) .GE. 0.9899_r8 ) ) THEN
                ztop(i) = 0.5_r8*(z2(i)+z1(i))
                ipbl(i) = k+1
                TEST(i) = .FALSE.
             END IF

             z1(i) = z2(i)
             t1(i) = t2(i)
             
          END IF

        END DO
     END DO

    RETURN

  END SUBROUTINE mpbl_depth





  SUBROUTINE radml_depth( &   
       IM           , &             !IM                    , &
       LM          , &             !LM                    , &
       ibot           , &             !ibot                   , &
       kcldtop     , &             !kcldtop         (1:IM,j)     , &
       svp         , &             !svpcp         (1:IM,j)     , &
       zt          , &             !zradtop         (1:IM,j)     , &
       critjump    , &             !critjump         (1:IM,j), &
       do_jump_exit, &             !do_jump_exit(1:IM,j), &
       t           , &             !slv         (1:IM,j,kcldtop(1:IM,j):ibot)/cp_air,                    &
       zf          , &             !zfull         (1:IM,j,kcldtop(1:IM,j):ibot),                      &
       zh          , &             !zhalf         (1:IM,j,kcldtop(1:IM,j):ibot), &
       zb          , &             !zradbase         (1:IM,j),          &
       zml         , &              !zradml         (1:IM,j))
       CONDIC        ) 
        
    !=======================================================================

    !======================================================================= 
    !
    !  Subroutine to calculate bottom and depth of radiatively driven mixed
    !  layer
    !
    !

    !
    !  -----
    !  INPUT
    !  -----
    !
    !  svp    cloud top liquid water virtual static energy divided by cp (K)
    !  zt     top of radiatively driven layer (m)
    !  t      liquid water virtual static energy divided by cp (K)
    !  zf     full level height above ground (m)
    !  zh     half level height above ground (m)
    !       
    !  ------
    !  OUTPUT
    !  ------
    !
    !  zb      base height of radiatively driven mixed layer (m)
    !  zml     depth of radiatively driven mixed layer (m)

    INTEGER, INTENT(IN) :: IM  
    INTEGER, INTENT(IN) :: LM  
    INTEGER, INTENT(IN) :: ibot        
    INTEGER, INTENT(IN) :: kcldtop    (1:IM) 
    REAL(KIND=r8)   ,   INTENT(in ) :: svp     (1:IM) 
    REAL(KIND=r8)   ,   INTENT(in ) :: zt      (1:IM) 
    REAL(KIND=r8)   ,   INTENT(in ) :: critjump(1:IM) 
    REAL(KIND=r8)   ,   INTENT(in ) :: t       (1:IM,1:LM) 
    REAL(KIND=r8)   ,   INTENT(in ) :: zf      (1:IM,1:LM) 
    REAL(KIND=r8)   ,   INTENT(in ) :: zh      (1:IM,0:LM) 
    REAL(KIND=r8)   ,   INTENT(out) :: zb      (1:IM) 
    REAL(KIND=r8)   ,   INTENT(out) :: zml     (1:IM) 
    LOGICAL         ,   INTENT(in ) :: do_jump_exit(1:IM) 
    LOGICAL         ,   INTENT(in ) :: CONDIC(1:IM) 

    REAL(KIND=r8)    :: svpar  (1:IM) 
    REAL(KIND=r8)    :: h1     (1:IM) 
    REAL(KIND=r8)    :: h2     (1:IM) 
    REAL(KIND=r8)    :: t1     (1:IM) 
    REAL(KIND=r8)    :: t2     (1:IM) 
    REAL(KIND=r8)    :: entrate
    REAL(KIND=r8)    :: entfr
    INTEGER :: k,i
    !INTEGER :: nlev
   zb =0.0_r8;   zml=0.0_r8;   svpar  =0.0_r8;   h1  =0.0_r8; h2  =0.0_r8;   t1  =0.0_r8;   t2  =0.0_r8;   entrate=0.0_r8
   entfr=0.0_r8
    DO i=1,IM
       IF (CONDIC(i) .and. kcldtop(i) .LT. ibot) THEN 

          !initialize zml
          zml(i) = 0.0_r8

          !compute # of levels
          !nlev = 1:IM

          !calculate cloud top parcel properties
          svpar(i)  = svp (i)
          h1   (i)  = zf  (i,1)
          t1   (i)  = t   (i,1)
      END IF
    END DO
    entrate = 0.2_r8/200.0_r8
    ! If above is false keep looking
    DO k =  MINVAL(kcldtop),ibot
       DO i=1,IM
          IF (CONDIC(i) .and. kcldtop(i) .LT. ibot) THEN 

             !search for level where parcel is warmer than env             
             ! first cut out if parcel is already warmer than
             ! cloudtop. 

             IF (t1(i).LT.svpar(i)) THEN
                 zb(i)  = h1(i)
                 zml(i) = 0.00_r8
                 RETURN
             END IF

             h2(i) = zf(i,k)
             t2(i) = t(i,k)

             IF (t2(i).LT.svpar(i)) THEN
                IF ( ABS(t1(i) - t2(i) ) .GT. 0.2_r8 ) THEN 
                   zb(i) = h2(i) + (h1(i) - h2(i))*(svpar(i) - t2(i))/(t1(i) - t2(i))
                   zb(i) = MAX( zb(i) , 0.0_r8 )  ! final protection against interp problems
                ELSE
                   zb(i) = h2(i)
                END IF
                zml(i) = zt(i) - zb(i)
                RETURN
             END IF

             IF (do_jump_exit(i) .AND. (t1(i)-t2(i)) .GT. critjump(i) .AND. k .GT. 2) THEN
                zb(i) = zh(i,k)
                zml(i) = zt(i) - zb(i)
                RETURN
             END IF

             entfr = MIN( entrate*(h1(i)-h2(i)), 1.0_r8 )
             svpar(i) = svpar(i) + entfr*(t2(i)-svpar(i))

             h1(i) = h2(i)
             t1(i) = t2(i)
          END IF
       END DO
    ENDDO
    DO i=1,IM
       IF (CONDIC(i) .and.kcldtop(i) .LT. ibot) THEN 

          zb(i) = 0.0_r8
          zml(i) = zt(i)
       END IF
    END DO

    RETURN
  END SUBROUTINE radml_depth



  SUBROUTINE diffusivity_pbl2(&
       ibot          , &          !ibot                    , &
       IM          , &          !IM                          , &
       LM          , &          !LM                          , &
       b_star          , &          !b_star     (1:IM,j)          , &
       ipbl          , &          !ipbl       (1:IM,j)          , &
       h          , &          !zsml       (1:IM,j)          , &
       kfac       , &          !khsfcfac                  , &
       k_ent      , &          !k_entr_tmp (1:IM,j)          , & 
       vsurf      , &          !vsurf      (1:IM,j)          , &
       frland     , &          !frland     (1:IM,j)          , & 
       zm         , &          !zhalf      (1:IM,j,0:LM), &
       k_m        , &          !k_m_troen  (1:IM,j,1:LM), &
       k_t           )          !k_t_troen  (1:IM,j,1:LM)   )

    !=======================================================================
    !========================================================================  
    !       Subroutine to return the vertical K-profile of diffusion 
    !       coefficients for the surface driven convective mixed layer.
    !       This code returns to form used Lock et al.. Does not match
    !       to surface layer.  
    !    
    !   call diffusivity_pbl2(h, k_ent, vsurf, zm, k_m, k_t)
    !                
    !      h:      Depth of surface driven mixed layer (m) 
    !      k_ent:  PBL top entrainment diffusivity (m+2 s-1)
    !      vsurf:  PBL top entrainment velocity scale (m s-1)
    !      zm:     Half level heights relative to the ground (m),        DIM[1:lm+1]
    !      k_m:    Momentum diffusion coefficient (m+2 s-1),             DIM[1:lm] (edges)
    !      k_t:    Heat and tracer diffusion coefficient (m+2 s-1)       DIM[1:lm] (edges)

    INTEGER, INTENT(in   ) :: ibot
    INTEGER, INTENT(in   ) :: IM
    INTEGER, INTENT(in   ) :: LM
    REAL(KIND=r8),    INTENT(in   ) :: b_star (1:IM)
    INTEGER, INTENT(in   ) :: ipbl   (1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: h      (1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: k_ent  (1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: vsurf  (1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: frland (1:IM)
    REAL(KIND=r8),    INTENT(in   ) :: kfac 
    REAL(KIND=r8),    INTENT(in   ) :: zm     (1:IM,0:LM)
    REAL(KIND=r8),    INTENT(out  ) :: k_m    (1:IM,1:LM)
    REAL(KIND=r8),    INTENT(out  ) :: k_t    (1:IM,1:LM)
    REAL(KIND=r8)     :: EE(1:IM)
    REAL(KIND=r8)     :: hin(1:IM)
    REAL(KIND=r8)     :: kfacx (1:IM)
    INTEGER  :: k,i
    INTEGER  :: kk
    k_m =0.0_r8; k_t=0.0_r8; EE=0.0_r8 ;hin=0.0_r8;kfacx=0.0_r8
    !lm = SIZE(zm,1)-1
    DO i=1,im

       !-----------------------------------------------------------
       !
       ! SURFACE DRIVEN CONVECTIVE LAYERS
       !
       ! Note this part is done only if b_star > 0., that is,
       ! upward surface buoyancy flux
       IF (b_star(i) .GT. 0.) THEN
          IF (ipbl(i) .LT. ibot) THEN

             ! Kluge. Raise KHs over land 
             !---------------------------
             IF ( frland(i) < 0.5_r8 ) THEN 
                kfacx(i) = kfac 
             ELSE
                kfacx(i) = kfac*2.0_r8
             END IF

             hin(i)    = 0.0_r8 ! 200._r8  ! "Organization" scale for plume (m).

          END IF 
       END IF
    END DO   


    DO k=1,LM
       DO i=1,im

       !-----------------------------------------------------------
       !
       ! SURFACE DRIVEN CONVECTIVE LAYERS
       !
       ! Note this part is done only if b_star > 0., that is,
       ! upward surface buoyancy flux
          IF (b_star(i) .GT. 0.0_r8) THEN
             IF (ipbl(i) .LT. ibot) THEN
                k_m(i,k) = 0.0_r8
                k_t(i,k) = 0.0_r8
             END IF
          END IF
       END DO                
    END DO

    !! factor = (zm(k)/hinner)* (1.0 -(zm(k)-hinner)/(h-hinner))**2
    DO i=1,im

       !-----------------------------------------------------------
       !
       ! SURFACE DRIVEN CONVECTIVE LAYERS
       !
       ! Note this part is done only if b_star > 0., that is,
       ! upward surface buoyancy flux
          IF (b_star(i) .GT. 0.0_r8) THEN
             IF (ipbl(i) .LT. ibot) THEN
                 IF ( vsurf(i)*h(i) .GT. 0.0_r8 ) THEN
                     EE (i) = 1.0_r8 - SQRT( k_ent(i) / ( kfacx(i) * vonkarm * vsurf(i) * h(i) ) )
                     EE (i) = MAX( EE(i) , 0.7_r8 )  ! If EE is too small, then punt, as LCs
                 END IF
             END IF
          END IF
    END DO


    DO k=1,lm
       DO i=1,im
          !-----------------------------------------------------------
          !
          ! SURFACE DRIVEN CONVECTIVE LAYERS
          !
          ! Note this part is done only if b_star > 0., that is,
          ! upward surface buoyancy flux
          IF (b_star(i) .GT. 0.0_r8) THEN
             IF (ipbl(i) .LT. ibot) THEN
                IF ( vsurf(i)*h(i) .GT. 0.0_r8 ) THEN
                   IF( ( zm(i,k) .LE. h(i) ) .AND.  ( zm(i,k) .GT. hin(i) )  ) THEN
                      !WHERE ( ( zm(1:lm) .LE. h ) .AND.  ( zm(1:lm) .GT. hin(i) )  )
                      k_t(i,k) = kfacx(i) * vonkarm * vsurf(i) * ( zm(i,k)-hin(i) ) * &
                      ( 1.0_r8 - EE(i)*( (zm(i,k)-hin(i))/(h(i)-hin(i)) ))**2
                      k_m(i,k) = k_t(i,k)
                     !END WHERE
                   END IF
                END IF
             END IF
          END IF   
       END DO
    END DO

    RETURN
  END SUBROUTINE diffusivity_pbl2
END MODULE PBL_Entrain
