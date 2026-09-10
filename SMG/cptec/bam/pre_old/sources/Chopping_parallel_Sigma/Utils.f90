!
!  $Author: pkubota $
!  $Date: 2011/04/07 16:00:31 $
!  $Revision: 1.18 $
!
MODULE Utils

  ! CreateAssocLegFunc
  ! DestroyAssocLegFunc
  !
  ! CreateGaussQuad    ------------------| CreateLegPol
  !
  ! CreateGridValues
  !
  ! DestroyGaussQuad   ------------------| DestroyLegPol
  !
  ! iminv
  !
  ! Rg                 ------------------| Balanc
  !                                  |
  !                                  | Orthes
  !                                  |
  !                                  | Ortran
  !                                  |
  !                                  | Hqr2    ------|Hqr3
  !                                  |
  !                                  | Balbak
  !                                  |
  !                                  | Znorma
  !
  ! Tql2
  ! Tred2
  ! tmstmp2
  ! InitTimeStamp
  !
  ! TimeStamp  --------------------------| caldat
  !
  ! IBJBtoIJ_R (Interface)
  ! IJtoIBJB_R (Interface)
  ! IBJBtoIJ_I (Interface)
  ! IJtoIBJB_I (Interface)
  !
  ! SplineIJtoIBJB_R2D (Interface) ------| CyclicCubicSpline
  !
  ! SplineIBJBtoIJ_R2D (Interface) ------| CyclicCubicSpline
  !
  ! LinearIJtoIBJB_R2D (Interface) ------| CyclicLinear
  !
  ! LinearIBJBtoIJ_R2D (Interface) ------| CyclicLinear
  !
  ! NearestIJtoIBJB_I2D (Interface)------| CyclicNearest_i
  !
  ! NearestIJtoIBJB_R2D (Interface)------| CyclicNearest_r
  !
  ! NearestIBJBtoIJ_I2D (Interface)------| CyclicNearest_i
  !
  ! NearestIBJBtoIJ_R2D (Interface)------| CyclicNearest_r
  !
  ! FreqBoxIJtoIBJB_I2D (Interface)------| CyclicFreqBox_i
  !
  ! FreqBoxIJtoIBJB_R2D (Interface)------| CyclicFreqBox_r
  !
  ! SeaMaskIJtoIBJB_R2D (Interface)------| CyclicSeaMask_r
  !
  ! SeaMaskIBJBtoIJ_R2D (Interface)------| CyclicSeaMask_r
  !
  ! AveBoxIJtoIBJB_R2D (Interface) ------| CyclicAveBox_r
  !
  ! AveBoxIBJBtoIJ_R2D (Interface) ------| CyclicAveBox_r
  !
  ! vfirec     --------------------------| vfinit
  !
  !
  !-------------------------------------------------------------------
  !  ASSOCIATED LEGENDRE FUNCTIONS
  !  Module computes and stores Associated Legendre
  !  Functions and Epslon.
  !
  !  Module exports three routines:
  !     CreateAssocLegFunc  initializes module and compute functions;
  !     DestroyAssocLegFunc destroys module;
  !
  !  Module usage:
  !  CreateAssocLegFunc should be invoked once, prior to any other
  !  module routine. It computes and hides the function values.
  !  DestroyAssocLegFunc destroys all internal info.
  !
  !  Module use values exported by Sizes and procedures from Auxiliary



  USE Sizes, ONLY:  &
       mnMax,       &
       mnMap,       &
       mMax,        &
       nMax,        &
       nExtMax,     &
       mnExtMax,    &
       mymnExtMax,  &
       iMax,        &
       jMax,        &
       kMax,        &
       jMaxHalf,    &
       jbMax,       &
       ibMax,       &
       ibPerIJ,     &
       jbPerIJ,     &
       iPerIJB,     &
       jPerIJB,     &
       ibMaxPerJB,  &
       iMaxPerJ,    &
       lm2m,        &
       mExtMap,     &
       nExtMap,     &
       mymMax,      &
       mymnExtMap,  &
       mymExtMap,   &
       mynExtMap,   &
       myfirstlat,  &
       myfirstlon,  &
       mylastlon,   &
       mylastlat,   &
       mnExtMap

  USE Parallelism, ONLY: &
       myId, &
       MsgOne, &
       FatalError


  USE InputParameters, ONLY : r8, r16, nfprt, pai, twomg, emrad, &
                              ImaxOut, JmaxOut, KmaxInp, KmaxInpp, KmaxOut, KmaxOutp, &
                              Rd, Cp, GEps, Gama, Grav, Rv

  USE InputArrays, ONLY : DelSigmaOut, SigLayerOut, SigInterOut, &
                          DelSInp, SigLInp, SigIInp, &
                          gTopoInp, gPsfcInp, gTvirInp, gPresInp, &
                          gTopoOut, gPsfcOut, gLnPsOut

  IMPLICIT NONE

  PRIVATE

  !
  !  LEGANDRE POLINOMIAL AND ITS ROOTS
  !
  !  Module exports four routines:
  !     CreateLegPol  initializes module;
  !     DestroyLegPol destroys module;
  !     LegPol        computes polinomial
  !     LegPolRoots   computes roots of even degree Legandre Pol
  !
  !  Module does not export (or require) any data value.
  !


  PUBLIC :: NewSigma
  PUBLIC :: SigmaInp
  PUBLIC :: NewPs
  PUBLIC :: CreateLegPol
  PUBLIC :: DestroyLegPol
  PUBLIC :: LegPol
  PUBLIC :: LegPolRootsandWeights
  PUBLIC :: CreateGridValues
  PUBLIC :: CreateGaussQuad
  PUBLIC :: DestroyGaussQuad
  PUBLIC :: GaussColat
  PUBLIC :: SinGaussColat
  PUBLIC :: CosGaussColat
  PUBLIC :: AuxGaussColat
  PUBLIC :: GaussPoints
  PUBLIC :: GaussWeights
  PUBLIC :: colrad
  PUBLIC :: colrad2D
  PUBLIC :: cos2lat
  PUBLIC :: ercossin
  PUBLIC :: fcor
  PUBLIC :: cosiv
  PUBLIC :: allpolynomials
  PUBLIC :: coslatj
  PUBLIC :: sinlatj
  PUBLIC :: coslat
  PUBLIC :: sinlat
  PUBLIC :: coslon
  PUBLIC :: sinlon
  PUBLIC :: longit
  PUBLIC :: rcl
  PUBLIC :: rcs2
  PUBLIC :: lonrad
  PUBLIC :: lati
  PUBLIC :: long
  PUBLIC :: cosz
  PUBLIC :: cos2d
  PUBLIC :: NoBankConflict
  PUBLIC :: CreateAssocLegFunc
  PUBLIC :: DestroyAssocLegFunc
  PUBLIC :: Reset_Epslon_To_Local
  PUBLIC :: Epslon
  PUBLIC :: LegFuncS2F
  PUBLIC :: IBJBtoIJ
  PUBLIC :: IJtoIBJB
  PUBLIC :: SplineIJtoIBJB
  PUBLIC :: SplineIBJBtoIJ
  PUBLIC :: LinearIJtoIBJB
  PUBLIC :: LinearIBJBtoIJ
  PUBLIC :: NearestIBJBtoIJ
  PUBLIC :: NearestIJtoIBJB
  PUBLIC :: FreqBoxIJtoIBJB
  PUBLIC :: SeaMaskIBJBtoIJ
  PUBLIC :: SeaMaskIJtoIBJB
  PUBLIC :: AveBoxIBJBtoIJ
  PUBLIC :: AveBoxIJtoIBJB
  PUBLIC :: CyclicNearest_r
  PUBLIC :: CyclicLinear
  PUBLIC :: CyclicLinear_ABS

  PUBLIC :: vfirec

  INTERFACE IBJBtoIJ
     MODULE PROCEDURE IBJBtoIJ_R, IBJBtoIJ_I
  END INTERFACE

  INTERFACE IJtoIBJB
     MODULE PROCEDURE &
          IJtoIBJB_R, IJtoIBJB_I, &
          IJtoIBJB3_R, IJtoIBJB3_I
  END INTERFACE

  INTERFACE SplineIBJBtoIJ
     MODULE PROCEDURE SplineIBJBtoIJ_R2D
  END INTERFACE

  INTERFACE SplineIJtoIBJB
     MODULE PROCEDURE SplineIJtoIBJB_R2D
  END INTERFACE

  INTERFACE LinearIBJBtoIJ
     MODULE PROCEDURE LinearIBJBtoIJ_R2D
  END INTERFACE

  INTERFACE LinearIJtoIBJB
     MODULE PROCEDURE LinearIJtoIBJB_R2D
  END INTERFACE

  INTERFACE NearestIBJBtoIJ
     MODULE PROCEDURE NearestIBJBtoIJ_I2D, NearestIBJBtoIJ_R2D
  END INTERFACE

  INTERFACE NearestIJtoIBJB
     MODULE PROCEDURE &
          NearestIJtoIBJB_I2D, NearestIJtoIBJB_R2D, &
          NearestIJtoIBJB_I3D, NearestIJtoIBJB_R3D
  END INTERFACE

  INTERFACE SeaMaskIBJBtoIJ
     MODULE PROCEDURE  SeaMaskIBJBtoIJ_R2D
  END INTERFACE

  INTERFACE SeaMaskIJtoIBJB
     MODULE PROCEDURE SeaMaskIJtoIBJB_R2D
  END INTERFACE

  INTERFACE FreqBoxIJtoIBJB
     MODULE PROCEDURE FreqBoxIJtoIBJB_I2D, FreqBoxIJtoIBJB_R2D
  END INTERFACE

  INTERFACE AveBoxIBJBtoIJ
     MODULE PROCEDURE AveBoxIBJBtoIJ_R2D
  END INTERFACE

  INTERFACE AveBoxIJtoIBJB
     MODULE PROCEDURE AveBoxIJtoIBJB_R2D
  END INTERFACE

  INTERFACE NoBankConflict
     MODULE PROCEDURE NoBankConflictS, NoBankConflictV
  END INTERFACE

  !  Module usage:
  !     CreateGaussQuad  should be invoked once, before any other routine
  !                      of this module, to set up maximum degree
  !                      of base functions (say, n);
  !                      it computes and hides n Gaussian Points and Weights
  !                      over interval [-1:1];
  !                      it also creates and uses Mod LegPol
  !     DestroyGaussQuad destrois hidden data structure and leaves module ready
  !                      for re-start, if desired, with another maximum degree.

  REAL(KIND=r8), ALLOCATABLE :: GaussColat(:)
  REAL(KIND=r8), ALLOCATABLE :: SinGaussColat(:)
  REAL(KIND=r8), ALLOCATABLE :: CosGaussColat(:)
  REAL(KIND=r8), ALLOCATABLE :: AuxGaussColat(:)
  REAL(KIND=r8), ALLOCATABLE :: GaussPoints(:)
  REAL(KIND=r8), ALLOCATABLE :: GaussWeights(:)
  REAL(KIND=r8), ALLOCATABLE :: auxpol(:,:)
  REAL(KIND=r8), ALLOCATABLE :: colrad(:)
  REAL(KIND=r8), ALLOCATABLE :: rcs2(:)
  REAL(KIND=r8), ALLOCATABLE :: colrad2D(:,:)
  REAL(KIND=r8), ALLOCATABLE :: cos2lat(:,:)
  REAL(KIND=r8), ALLOCATABLE :: ercossin(:,:)
  REAL(KIND=r8), ALLOCATABLE :: fcor(:,:)
  REAL(KIND=r8), ALLOCATABLE :: cosiv(:,:)
  REAL(KIND=r8), ALLOCATABLE :: coslatj(:)
  REAL(KIND=r8), ALLOCATABLE :: sinlatj(:)
  REAL(KIND=r8), ALLOCATABLE :: coslat(:,:)
  REAL(KIND=r8), ALLOCATABLE :: sinlat(:,:)
  REAL(KIND=r8), ALLOCATABLE :: coslon(:,:)
  REAL(KIND=r8), ALLOCATABLE :: sinlon(:,:)
  REAL(KIND=r8), ALLOCATABLE :: longit(:,:)
  REAL(KIND=r8), ALLOCATABLE :: rcl(:,:)
  REAL(KIND=r8), ALLOCATABLE :: lonrad(:,:)
  REAL(KIND=r8), ALLOCATABLE :: lati(:)
  REAL(KIND=r8), ALLOCATABLE :: long(:)
  REAL(KIND=r8), ALLOCATABLE :: cosz(:)
  REAL(KIND=r8), ALLOCATABLE :: cos2d(:,:)

  !  Module Hided Data:
  !     maxDegree is the degree of the base functions (n)
  !     created specifies if module was created or not



  LOGICAL           :: created=.FALSE.
  LOGICAL           :: reducedgrid=.FALSE.
  LOGICAL           :: allpolynomials
  INTEGER           :: maxDegree=-1



  REAL(KIND=r8),    ALLOCATABLE :: Epslon(:)
  REAL(KIND=r8),    ALLOCATABLE :: LegFuncS2F(:,:)
  REAL(KIND=r8),    ALLOCATABLE :: Square(:)
  REAL(KIND=r8),    ALLOCATABLE :: Den(:)


  !  Module Hidden data


  INTEGER                        :: nAuxPoly
  REAL(KIND=r16), ALLOCATABLE, DIMENSION(:) :: AuxPoly1, AuxPoly2
  LOGICAL, PARAMETER   :: dumpLocal=.FALSE.

  ! Index mappings to/from diagonal from/to column for
  ! 'extended' spectral representations

  INTEGER,  ALLOCATABLE :: ExtDiagPerCol(:)   ! diag=DiagPerCol(col )
  INTEGER,  ALLOCATABLE :: ExtColPerDiag(:)   ! col =ColPerDiag(diag)

  ! date are always in the form yyyymmddhh ( year, month, day, hour). hour in
  ! is in 0-24 form

  INTEGER, PRIVATE :: JulianDayInitIntegration

CONTAINS



  !CreateAssocLegFunc  should be invoked once, before any other routine,
  !                    to compute and store Associated Legendre Functions
  !                    at Gaussian Points for all
  !                    Legendre Orders and Degrees (defined at Sizes)

SUBROUTINE NewSigma

  ! Given Delta Sigma Values Computes
  ! Sigma Interface and Layer Values

  IMPLICIT NONE

  INTEGER :: k

  REAL (KIND=r8) :: RdoCp, CpoRd, RdoCp1

  ! Compute New Sigma Interface Values
  SigInterOut(1)=1.0_r8
  DO k=1,KmaxOut
     SigInterOut(k+1)=SigInterOut(k)-DelSigmaOut(k)
  END DO
  SigInterOut(KmaxOutp)=0.0_r8

  ! Compute New Sigma Layer Values
  RdoCp=Rd/Cp
  CpoRd=Cp/Rd
  RdoCp1=RdoCp+1.0_r8
  DO k=1,KmaxOut
     SigLayerOut(k)=((SigInterOut(k)**RdoCp1-SigInterOut(k+1)**RdoCp1)/ &
                    (RdoCp1*(SigInterOut(k)-SigInterOut(k+1))))**CpoRd
  END DO

END SUBROUTINE NewSigma


SUBROUTINE SigmaInp

  ! Given Delta Sigma Values Computes
  ! Sigma Interface and Layer Values

  IMPLICIT NONE

  REAL (KIND=r8) :: RdoCp, CpoRd, RdoCp1

  INTEGER :: k

  ! Compute New Sigma Interface Values
  SigIInp(1)=1.0_r8
  DO k=1,KmaxInp
     SigIInp(k+1)=SigIInp(k)-DelSInp(k)
  END DO
  SigIInp(KmaxInpp)=0.0_r8

  ! Compute New Sigma Layer Values
  RdoCp=Rd/Cp
  CpoRd=Cp/Rd
  RdoCp1=RdoCp+1.0_r8
  DO k=1,KmaxInp
     SigLInp(k)=((SigIInp(k)**RdoCp1-SigIInp(k+1)**RdoCp1)/ &
                (RdoCp1*(SigIInp(k)-SigIInp(k+1))))**CpoRd
  END DO

END SUBROUTINE SigmaInp


SUBROUTINE NewPs

  ! Computes a New Surface Pressure given a New Orography.
  ! The New Pressure is computed assuming a Hydrostatic Balance
  ! and a Constant Temperature Lapse Rate.  Below Ground, the
  ! Lapse Rate is assumed to be -6.5 K/km.

  IMPLICIT NONE

  INTEGER :: ls, k, i, j

  REAL (KIND=r8) :: goRd, GammaLR

  REAL (KIND=r8), DIMENSION (Ibmax,Jbmax) :: Zu

  goRd=Grav/Rd

  ! Compute Surface Pressure Below the Original Ground
  ls=0
  k=1
  GammaLR=Gama
  DO j=1,Jbmax  
     DO i=1,Ibmaxperjb(j)
        Zu(i,j)=gTopoInp(i,j)-gTvirInp(i,k,j)/GammaLR* &
                ((gPsfcInp(i,j)/gPresInp(i,k,j))**(-GammaLR/goRd)-1.0_r8)
        IF (gTopoOut(i,j) <= Zu(i,j)) THEN
           IF (ABS(GammaLR) > GEps) THEN
              gPsfcOut(i,j)=gPresInp(i,k,j)*(1.0_r8+GammaLR/gTvirInp(i,k,j)* &
                            (gTopoOut(i,j)-Zu(i,j)))**(-goRd/GammaLR)
           ELSE
              gPsfcOut(i,j)=gPresInp(i,k,j)*EXP(-goRd/gTvirInp(i,k,j)* &
                            (gTopoOut(i,j)-Zu(i,j)))
           END IF
        ELSE
           gPsfcOut(i,j)=0.0_r8
           ls=ls+1
        END IF
     END DO
  END DO

  ! Compute Surface Pressure Above the Original Ground
  DO k=2,KmaxInp
     IF (ls > 0) THEN
        DO j=1,Jbmax
           DO i=1,Ibmaxperjb(j)
              IF (gPsfcOut(i,j) == 0.0_r8) THEN
                 GammaLR=-goRd*LOG(gTvirInp(i,k-1,j)/gTvirInp(i,k,j))/ &
                         LOG(gPresInp(i,k-1,j)/gPresInp(i,k,j))
                 IF (ABS(GammaLR) > GEps) THEN
                    Zu(i,j)=Zu(i,j)-gTvirInp(i,k,j)/GammaLR* &
                            ((gPresInp(i,k-1,j)/gPresInp(i,k,j))** &
                            (-GammaLR/goRd)-1.0_r8)
                 ELSE
                    Zu(i,j)=Zu(i,j)+gTvirInp(i,k,j)/ &
                            goRd*LOG(gPresInp(i,k-1,j)/gPresInp(i,k,j))
                 END IF
                 IF (gTopoOut(i,j) <= Zu(i,j)) THEN
                    IF (ABS(GammaLR) > GEps) THEN
                       gPsfcOut(i,j)=gPresInp(i,k,j)*(1.0_r8+GammaLR/gTvirInp(i,k,j)* &
                                     (gTopoOut(i,j)-Zu(i,j)))**(-goRd/GammaLR)
                    ELSE
                       gPsfcOut(i,j)=gPresInp(i,k,j)* &
                                     EXP(-goRd/gTvirInp(i,k,j)*(gTopoOut(i,j)-Zu(i,j)))
                    END IF
                    ls=ls-1
                 END IF
              END IF
           END DO
        END DO
     END IF
  END DO

  ! Compute Surface Pressure Over the Top
  IF (ls > 0) THEN
     k=KmaxInp
     GammaLR=0.0_r8
     DO j=1,Jbmax
        DO i=1,Ibmaxperjb(j)
           IF (gPsfcOut(i,j) == 0.0_r8) THEN
              gPsfcOut(i,j)=gPresInp(i,k,j)* &
                            EXP(-goRd/gTvirInp(i,k,j)*(gTopoOut(i,j)-Zu(i,j)))
           END IF
        END DO
     END DO
  END IF

  ! ln(ps) in cBar (from mBar)     
  DO j=1,Jbmax        
     DO i=1,Ibmaxperjb(j)
      gLnPsOut(i,j)=LOG(gPsfcOut(i,j)/10.0_r8)
    END DO
  END DO  
END SUBROUTINE NewPs




  SUBROUTINE CreateAssocLegFunc(allpolynomials)
    LOGICAL, INTENT(IN) :: allpolynomials
    INTEGER:: m, n, mn, j, mp, np, jp, mglobalp
    CHARACTER(LEN=*), PARAMETER :: h="**(CreateAssocLegFunc)**"

    IF (allpolynomials) THEN
      ALLOCATE (Epslon(mnExtMax))
      ALLOCATE (LegFuncS2F(jMaxHalf, mnExtMax))
     ELSE
      ALLOCATE (Epslon(mymnExtMax))
      ALLOCATE (LegFuncS2F(jMaxHalf, mymnExtMax))
      ALLOCATE (auxpol(jMaxHalf, mMax))
    ENDIF

    ALLOCATE (Square(nExtMax))
    ALLOCATE (Den(nExtMax))
    DO n = 1, nExtMax
       Square(n) = REAL((n-1)*(n-1),r8)
       Den(n)    = 1.0_r8/(4.0_r8*Square(n) - 1.0_r8)
    END DO
    IF (allpolynomials) THEN
       DO mn = 1, mnExtMax
          m = mExtMap(mn)
          n = nExtMap(mn)
          Epslon(mn) = SQRT((Square(n)-Square(m))*Den(n))
       END DO
     ELSE
       DO mn = 1, mymnExtMax
          m = lm2m(mymExtMap(mn))
          n = mynExtMap(mn)
          Epslon(mn) = SQRT((Square(n)-Square(m))*Den(n))
       END DO
    END IF


    IF (allpolynomials) THEN
       LegFuncS2F(1:jMaxHalf, mnExtMap(1,1)) = SQRT(0.5_r8)
       DO m = 2, mMax
          DO j = 1, jMaxHalf
             LegFuncS2F(j, mnExtMap(m,m)) =   &
                  SQRT(1.0_r8 + 0.5_r8/REAL(m-1,r8)) * &
                  SinGaussColat(j)          * &
                  LegFuncS2F(j, mnExtMap(m-1,m-1))
          END DO
       END DO
       !$OMP PARALLEL DO PRIVATE(mp,jp)
       DO mp = 1, mMax
          DO jp = 1, jMaxHalf
             LegFuncS2F(jp, mnExtMap(mp,mp+1)) =   &
                  SQRT(1.0_r8 + 2.0_r8*REAL(mp,r8))     * &
                  GaussPoints(jp)              * &
                  LegFuncS2F(jp, mnExtMap(mp,mp))
          END DO
       END DO
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(mp,jp,np)
       DO mp = 1, mMax
          DO np = mp+2, nExtMax
             DO jp = 1, jMaxHalf
                LegFuncS2F(jp, mnExtMap(mp,np)) =          &
                     ( GaussPoints(jp)                 * &
                       LegFuncS2F(jp, mnExtMap(mp,np-1)) - &
                       Epslon(mnExtMap(mp,np-1))        * &
                       LegFuncS2F(jp, mnExtMap(mp,np-2))   &
                      ) / Epslon(mnExtMap(mp,np))
             END DO
          END DO
       END DO
       !$OMP END PARALLEL DO
       CALL Reset_Epslon_To_Local ()
     ELSE
       auxpol(1:jMaxHalf, 1) = SQRT(0.5_r8)
       DO m = 2, mMax
          DO j = 1, jMaxHalf
             auxpol(j, m) =   &
                  SQRT(1.0_r8 + 0.5_r8/REAL(m-1,r8)) * &
                  SinGaussColat(j)          * &
                  auxpol(j, m-1)
          END DO
       END DO
       !$OMP PARALLEL DO PRIVATE(mp,jp,mglobalp)
       DO mp = 1, mymMax
          mglobalp = lm2m(mp)
          DO jp = 1, jMaxHalf
             LegFuncS2F(jp, mymnExtMap(mp,mglobalp)) = &
                  auxpol(jp, mglobalp)
          END DO
       END DO
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(mp,jp,mglobalp)
       DO mp = 1, mymMax
          mglobalp = lm2m(mp)
          DO jp = 1, jMaxHalf
             LegFuncS2F(jp, mymnExtMap(mp,mglobalp+1)) =   &
                  SQRT(1.0_r8 + 2.0_r8*REAL(mglobalp,r8))     * &
                  GaussPoints(jp)              * &
                  LegFuncS2F(jp, mymnExtMap(mp,mglobalp))
          END DO
       END DO
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(mp,jp,mglobalp,np)
       DO mp = 1, mymMax
          mglobalp = lm2m(mp)
          DO np = mglobalp+2, nExtMax
             DO jp = 1, jMaxHalf
                LegFuncS2F(jp, mymnExtMap(mp,np)) =          &
                     ( GaussPoints(jp)                 * &
                       LegFuncS2F(jp, mymnExtMap(mp,np-1)) - &
                       Epslon(mymnExtMap(mp,np-1))        * &
                       LegFuncS2F(jp, mymnExtMap(mp,np-2))   &
                      ) / Epslon(mymnExtMap(mp,np))
             END DO
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE CreateAssocLegFunc


  SUBROUTINE Reset_Epslon_To_Local ()
  INTEGER:: m, n, mn

!
!   Convert Epslon to local mpi structure
!
    DO mn = 1, mymnExtMax
       m = lm2m(mymExtMap(mn))
       n = mynExtMap(mn)
       Epslon(mn) = SQRT((Square(n)-Square(m))*Den(n))
    END DO

  END SUBROUTINE Reset_Epslon_To_Local


  !DestroyAssocLegFunc  Deallocates all stored values



  SUBROUTINE DestroyAssocLegFunc()
    CHARACTER(LEN=*), PARAMETER :: h="**(DestroyAssocLegFunc)**"
    DEALLOCATE (Epslon, LegFuncS2F)
    DEALLOCATE (Square, Den)
  END SUBROUTINE DestroyAssocLegFunc



  !  AUXILIARY PROCEDURES

  !  Module exports two routines:
  !     NoBankConflict  given input integer (size of an array at
  !                     any dimension) returns the next integer
  !                     that should dimension the array to avoid
  !                     memory bank conflicts. The vector version
  !                     returns a vector of integers, given a vector
  !                     of input integers.

  !  Module does not require any other module.
  !  Module does not export any value.



  FUNCTION NoBankConflictS(s) RESULT(p)
    INTEGER, INTENT(IN) :: s
    INTEGER             :: p
    IF ((MOD(s,2)==0) .AND. (s/=0)) THEN
       p = s + 1
    ELSE
       p = s
    END IF
  END FUNCTION NoBankConflictS
  FUNCTION NoBankConflictV(s) RESULT(p)
    INTEGER, INTENT(IN) :: s(:)
    INTEGER             :: p(SIZE(s))
    WHERE ((MOD(s,2)==0) .AND. (s/=0))
       p = s + 1
    ELSEWHERE
       p = s
    END WHERE
  END FUNCTION NoBankConflictV

  !  GAUSSIAN POINTS AND WEIGHTS FOR QUADRATURE
  !  OVER LEGENDRE POLINOMIALS BASE FUNCTIONS

  !  Module exports three routines:
  !     CreateGaussQuad   initializes module
  !     DestroyGaussQuad  destroys module

  !  Module export two arrays:
  !     GaussPoints and GaussWeights
  !
  !  Module uses Module LegPol (Legandre Polinomials);
  !  transparently to the user, that does not have
  !  to create and/or destroy LegPol






  !CreateGaussQuad  computes and hides 'degreeGiven' Gaussian Points
  !                 and Weights over interval [-1:1];
  !                 creates and uses Mod LegPol.



  SUBROUTINE CreateGaussQuad (degreeGiven)
    INTEGER, INTENT(IN) :: degreeGiven
    REAL(KIND=r16), ALLOCATABLE :: FVals(:)
    CHARACTER(LEN=*), PARAMETER :: h="**(CreateGaussQuad)**"
    CHARACTER(LEN=256) :: line
    CHARACTER(LEN=10) :: c1
    INTEGER :: j

    !check invocation sequence and input data

    IF (degreeGiven <=0) THEN
       WRITE(c1,"(i10)") degreeGiven
       WRITE(nfprt,"(a,' invoked with degree ',a)") h, TRIM(ADJUSTL(c1))
       STOP
    ELSE
       maxDegree = degreeGiven
    END IF

    !allocate areas

    ALLOCATE (GaussColat(maxDegree/2))
    ALLOCATE (SinGaussColat(maxDegree/2))
    ALLOCATE (CosGaussColat(maxDegree/2))
    ALLOCATE (AuxGaussColat(maxDegree/2))
    ALLOCATE (FVals(maxDegree/2))
    ALLOCATE (GaussPoints(maxDegree))
    ALLOCATE (GaussWeights(maxDegree))
    ALLOCATE (colrad(maxDegree))
    ALLOCATE (rcs2(maxDegree/2))

    !create ModLegPol

    CALL CreateLegPol (maxDegree)

    !Gaussian Points are the roots of legandre polinomial of degree maxDegree

    CALL LegPolRootsandWeights(maxDegree)
    GaussWeights(maxDegree/2+1:maxDegree) = GaussWeights(maxDegree/2:1:-1)
    GaussColat = ACOS(CosGaussColat)
    SinGaussColat = SIN(GaussColat)
    AuxGaussColat = 1.0_r16/(SinGaussColat*SinGaussColat)
    GaussPoints(1:maxDegree/2) = CosGaussColat
    GaussPoints(maxDegree/2+1:maxDegree) = -GaussPoints(maxDegree/2:1:-1)

    DO j=1,maxDegree/2
       colrad(j)=    GaussColat(j)
       colrad(maxDegree+1-j)=pai-GaussColat(j)
       rcs2(j)=AuxGaussColat(j)
    END DO

  END SUBROUTINE CreateGaussQuad


  SUBROUTINE CreateGridValues

    INTEGER :: i, j, ib, jb, jhalf
    REAL (KIND=r8)  :: rd
    ALLOCATE (colrad2D(ibMax,jbMax))
    ALLOCATE (rcl     (ibMax,jbMax))
    rcl=0.0_r8
    ALLOCATE (coslatj (jMax))
    ALLOCATE (sinlatj (jMax))
    ALLOCATE (coslat  (ibMax,jbMax))
    ALLOCATE (sinlat  (ibMax,jbMax))
    ALLOCATE (coslon  (ibMax,jbMax))
    ALLOCATE (sinlon  (ibMax,jbMax))
    ALLOCATE (cos2lat (ibMax,jbMax))
    ALLOCATE (ercossin(ibMax,jbMax))
    ALLOCATE (fcor    (ibMax,jbMax))
    ALLOCATE (cosiv   (ibMax,jbMax))
    ALLOCATE (longit  (ibMax,jbMax))
    ALLOCATE (lonrad  (ibMax,jbMax))
    ALLOCATE (lati    (jMax))
    ALLOCATE (long    (iMax))
    ALLOCATE (cosz    (jMax))
    ALLOCATE (cos2d   (ibMax,jbMax))

    rd=45.0_r8/ATAN(1.0_r8)
    DO j=1,jMax
       lati(j)=90.0_r8-colrad(j)*rd
    END DO
    DO i=1,imax
       long(i) = (i-1)*360.0_r8/REAL(iMax,r8)
    ENDDO

    !$OMP PARALLEL DO PRIVATE(jhalf)
    DO j=1,jMax/2
       jhalf = jMax-j+1
       sinlatj (j)=COS(colrad(j))
       coslatj (j)=SIN(colrad(j))
       sinlatj (jhalf)= - sinlatj(j)
       coslatj (jhalf)= coslatj(j)
    ENDDO
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO PRIVATE(ib,i,j,jhalf)
    DO jb = 1, jbMax
       DO ib = 1, ibMaxPerJB(jb)
          j = jPerIJB(ib,jb)
          i = iPerIJB(ib,jb)
          jhalf = MIN(j, jMax-j+1)
          colrad2D(ib,jb)=colrad(j)
          sinlat  (ib,jb)=sinlatj(j)
          coslat  (ib,jb)=coslatj(j)
          rcl     (ib,jb)=rcs2(jhalf)
          longit  (ib,jb)=(i-1)*pai*2.0_r8/iMaxPerJ(j)
          lonrad  (ib,jb)=(i-1)*360.0_r8/REAL(iMaxPerJ(j),r8)
          sinlon  (ib,jb)=SIN(longit(ib,jb))
          coslon  (ib,jb)=COS(longit(ib,jb))
          cos2lat (ib,jb)=1.0_r8/rcl(ib,jb)
!**(JP)** mudei
!          ercossin(ib,jb)=sinlat(ib,jb)*rcl(ib,jb)/er
!          fcor    (ib,jb)=twomg*sinlat(ib,jb)
!          cosiv   (ib,jb)=1./coslat(ib,jb)
! para a forma antiga, so para bater binario no euleriano;
! impacto no SemiLagrangeano eh desconhecido
          ercossin(ib,jb)=COS(colrad(j))*rcl(ib,jb)/emrad
          fcor    (ib,jb)=twomg*COS(colrad(j))
          cosiv   (ib,jb)=SQRT(rcl(ib,jb))
!**(JP)** fim de alteracao
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE CreateGridValues

  !DestroyGaussQuad  destroy module internal data;
  !                  destroys module LegPol;
  !                  get ready for new module usage



  SUBROUTINE DestroyGaussQuad
    CHARACTER(LEN=*), PARAMETER :: h="**(DestroyGaussQuad)**"

    !check invocation sequence

    maxDegree = -1

    !destroy ModLegPol

    CALL DestroyLegPol()

    !deallocate module areas
    !
    DEALLOCATE (GaussColat)
    DEALLOCATE (SinGaussColat)
    DEALLOCATE (AuxGaussColat)
    DEALLOCATE (GaussPoints)
    DEALLOCATE (GaussWeights)
  END SUBROUTINE DestroyGaussQuad


  !  Module usage:
  !     CreateLegPol should be invoked once, before any other routine;
  !     LegPol can be invoked after CreateLegPol, as much as required;
  !     LegPolRoots can be invoked after CreateLegPol, as much as required;
  !     DestroyLegPol should be invoked at the end of the computation.
  !     CreateLegPol estabilishes the maximum degree for which LegPol
  !     and LegPolRoots should be invoked. If this maximum
  !     degree has to be changed, module should be destroied and
  !     created again.
  !
  !     LegPol and LegPolRoots deal with colatitudes;
  !     LegPol take colatitudes as abcissas and
  !     LegPolRoots produces colatitudes as roots.
  !     Colatitudes are expressed in radians, with 0 at the
  !     North Pole, pi/2 at the Equator and pi at the South Pole.
  !
  !
  !  Module Hided Data:
  !     AuxPoly1, AuxPoly2 are constants used to evaluate LegPol
  !     nAuxPoly is the maximum degree of the Polinomial to be computed
  !     created specifies if module was created or not
  !

  !
  !  CreateLegPol does module initialization;
  !             should be executed once, prior to LegPol;
  !             input argument is the maximum degree of LegPol invocations
  !
  SUBROUTINE CreateLegPol (maxDegree)
    INTEGER, INTENT(IN) :: maxDegree
    INTEGER :: i
    !
    !  Check Correction
    !
    IF ( maxDegree <= 0) THEN
       WRITE (0,"('**(CreateLegPol)** maxDegree <= 0; maxDegree =')") maxDegree
       STOP
    ELSE IF (created) THEN
       WRITE (0,"('**(CreateLegPol)** invoked twice without destruction')")
       STOP
    END IF
    !
    !  Allocate and Compute Constants
    !
    nAuxPoly = maxDegree
    ALLOCATE (AuxPoly1(nAuxPoly))
    ALLOCATE (AuxPoly2(nAuxPoly))
    DO i = 1, nAuxPoly
       AuxPoly1(i) = REAL(2*i-1,r16)/REAL(i,r16)
       AuxPoly2(i) = REAL(1-i,r16)/REAL(i,r16)
    END DO
    created = .TRUE.
  END SUBROUTINE CreateLegPol
  !
  ! DestroyLegPol destroys module initialization;
  !               required if maximum degree of LegPol has to be changed;
  !               in this case, invoke DestroyLegPol and CreateLegPol with
  !               new degree
  !
  SUBROUTINE DestroyLegPol
    IF (.NOT. created) THEN
       WRITE (0,"('**(DestroyLegPol)** invoked without initialization')")
       STOP
    END IF
    DEALLOCATE (AuxPoly1, AuxPoly2)
    created = .FALSE.
  END SUBROUTINE DestroyLegPol
  !
  ! LegPol computes Legandre Polinomial of degree 'degree' at a
  !        vector of colatitudes 'Col'
  !        degree should be <= maximum degree estabilished by CreateLegPol
  !
  FUNCTION LegPol (degree, Col)
    INTEGER,                      INTENT(IN ) :: degree
    REAL(KIND=r16),    DIMENSION(:),        INTENT(IN ) :: Col
    REAL(KIND=r16),    DIMENSION(SIZE(Col))             :: LegPol
    !
    !  Auxiliary Variables
    !
    REAL(KIND=r16), DIMENSION(SIZE(Col)) :: P0     ! Polinomial of degree i - 2
    REAL(KIND=r16), DIMENSION(SIZE(Col)) :: P1     ! Polinomial of degree i - 1
    REAL(KIND=r16), DIMENSION(SIZE(Col)) :: X      ! Cosine of colatitude
    INTEGER :: iDegree                   ! loop index
    INTEGER :: left                      ! loop iterations before unrolling
    !
    !  Check Correctness
    !
    IF (.NOT. created) THEN
       WRITE (0,"('**(LegPol)** invoked without initialization')")
       STOP
    END IF
    !
    !  Case degree >=2 and degree <= maximum degree
    !
    IF ((degree >= 2) .AND. (degree <= nAuxPoly)) THEN
       !
       !  Initialization
       !
       left = MOD(degree-1,6)
       X  = COS(Col)
       P0 = 1.0_r16
       P1 = X
       !
       !  Apply recurrence relation
       !     Loop upper bound should be 'degree';
       !     it is not due to unrolling
       !
       DO iDegree = 2, left+1
          LegPol = X*P1*AuxPoly1(iDegree) + P0*AuxPoly2(iDegree)
          P0 = P1; P1 = LegPol
       END DO
       !
       !  Unroll recurrence relation, to speed up computation
       !
       DO iDegree = left+2, degree, 6
          LegPol = X*P1    *AuxPoly1(iDegree  ) + P0    *AuxPoly2(iDegree  )
          P0     = X*LegPol*AuxPoly1(iDegree+1) + P1    *AuxPoly2(iDegree+1)
          P1     = X*P0    *AuxPoly1(iDegree+2) + LegPol*AuxPoly2(iDegree+2)
          LegPol = X*P1    *AuxPoly1(iDegree+3) + P0    *AuxPoly2(iDegree+3)
          P0     = X*LegPol*AuxPoly1(iDegree+4) + P1    *AuxPoly2(iDegree+4)
          P1     = X*P0    *AuxPoly1(iDegree+5) + LegPol*AuxPoly2(iDegree+5)
       END DO
       LegPol = P1
    !
    !  Case degree == 0
    !
    ELSE IF (degree == 0) THEN
       LegPol = 1.0_r16
    !
    !  Case degree == 1
    !
    ELSE IF (degree == 1) THEN
       LegPol = COS(Col)
    !
    !  Case degree <= 0 or degree > maximum degree
    !
    ELSE
       WRITE(nfprt,"('**(LegPol)** invoked with degree ',i6,&
            &' out of bounds')") degree
       STOP
    END IF
  END FUNCTION LegPol
  !
  !  LegPolRootsandweights  computes the roots of the Legandre Polinomial
  !              of even degree 'degree' and respective Gaussian Weights
  !              Roots are expressed as colatitudes at interval
  !              (0, pi/2)
  !
  SUBROUTINE LegPolRootsandWeights (degree)
    INTEGER, INTENT(IN) :: degree
    REAL(KIND=r16),                 DIMENSION(degree/2) :: Col
    REAL(KIND=r16),    ALLOCATABLE, DIMENSION(:)      :: XSearch, FSearch
    REAL(KIND=r16), DIMENSION(degree/2) :: Pol    ! Polinomial of degree i
    REAL(KIND=r16), DIMENSION(degree/2) :: P0     ! Polinomial of degree i-2
    REAL(KIND=r16), DIMENSION(degree/2) :: P1     ! Polinomial of degree i-1
    REAL(KIND=r16), DIMENSION(degree/2) :: X
    REAL(KIND=r16), DIMENSION(degree/2) :: XC
    LOGICAL, ALLOCATABLE, DIMENSION(:)      :: Mask
    INTEGER :: i                               ! loop index
    INTEGER :: nPoints                         ! to start bissection
    INTEGER, PARAMETER :: multSearchStart=4 ! * factor to start bissection
    REAL(KIND=r16)    :: step
    REAL(KIND=r16)    :: pi, scale
    INTEGER, PARAMETER  :: itmax=10     ! maximum number of newton iterations
    INTEGER, PARAMETER  :: nDigitsOut=2 ! precision digits of gaussian points:
                                        ! machine epsilon - nDigitsOut
    REAL(KIND=r16)    :: rootPrecision   ! relative error in gaussian points

    INTEGER :: halfDegree
    INTEGER :: iDegree, it               ! loop index
    INTEGER :: left                      ! loop iterations before unrolling
    !
    !  Check Correctness
    !
    IF (.NOT. created) THEN
       WRITE (0,"('**(LegPolRootsandWeights)** invoked without initialization')")
       STOP
    ELSE IF ( (degree <= 0) .OR. (degree > nAuxPoly) ) THEN
       WRITE (0,"('**(LegPolRootsandWeights)** invoked with degree ',i6,&
            &' out of bounds')") degree
       STOP
    ELSE IF ( MOD(degree,2) .NE. 0) THEN
       WRITE (0,"('**(LegPolRootsandWeights)** invoked with odd degree ',i6)") degree
       STOP
    END IF
    !
    !  Initialize Constants
    !
    pi = 4.0_r16 * ATAN(1.0_r16)
    rootPrecision = EPSILON(1.0_r16)*10.0_r16**(nDigitsOut)
    !
    !  LegPolRoots uses root simmetry with respect to pi/2.
    !  It finds all roots in the interval [0,pi/2]
    !  Remaining roots are simmetric
    !
    halfDegree = degree/2
    !
    !  bissection method to find roots:
    !  get equally spaced points in interval [0,pi/2]
    !  to find intervals containing roots
    !
    nPoints = multSearchStart*halfDegree
    step = pi/(2.0_r16*REAL(nPoints,r16))
    ALLOCATE (XSearch(nPoints))
    ALLOCATE (FSearch(nPoints))
    ALLOCATE (Mask   (nPoints-1))
    DO i = 1, nPoints
       XSearch(i) = step*REAL(i-1,r16)
    END DO
    FSearch = LegPol (degree, XSearch)
    !
    !  select intervals containing roots
    !
    Mask = FSearch(1:nPoints-1)*FSearch(2:nPoints) < 0.0
    !
    !  are there enough intervals?
    !
    IF (COUNT(Mask) .NE. halfDegree) THEN
       WRITE(nfprt,"('**(LegPolRoots)** ',i6,' bracketing intervals to find '&
            &,i6,' roots')") COUNT(Mask), halfDegree
       STOP
    END IF
    !
    !  extract intervals containing roots
    !
    Col = 0.5_r16 * (PACK(XSearch(1:nPoints-1), Mask) + &
                  PACK(XSearch(2:nPoints  ), Mask) )
    DEALLOCATE (XSearch)
    DEALLOCATE (FSearch)
    DEALLOCATE (Mask)
    !
    !    loop while there is a root to be found
    !
    it = 1
    X  = COS(Col)
    scale = 2.0_r16/REAL(degree*degree,r16)
    left = MOD(degree-1,6)
!   WRITE (UNIT=*, FMT='(/,I6,2I5,1P2E16.8)') degree, r16, r16, rootPrecision,ATAN(1.0_r16)/REAL(degree,r16)
    DO
       IF (it.gt.itmax) THEN
         WRITE (0,"('**(LegPolRootsandWeights)** failed to converge  ',i6,&
                   &' itmax')") itmax
         EXIT
       END IF
       !
       !   initialization
       !
       it = it + 1
       P0 = 1.0_r16
       P1 = X
       !
       !   Apply recurrence relation
       !     Loop upper bound should be 'degree';
       !     it is not due to unrolling
       !
       DO iDegree = 2, left+1
          Pol = X*P1*AuxPoly1(iDegree) + P0*AuxPoly2(iDegree)
          P0 = P1; P1 = Pol
       END DO
       !
       !   Unroll recurrence relation, to speed up computation
       !
       DO iDegree = left+2, degree, 6
          Pol = X*P1    *AuxPoly1(iDegree  ) + P0    *AuxPoly2(iDegree  )
          P0     = X*Pol*AuxPoly1(iDegree+1) + P1    *AuxPoly2(iDegree+1)
          P1     = X*P0    *AuxPoly1(iDegree+2) + Pol*AuxPoly2(iDegree+2)
          Pol = X*P1    *AuxPoly1(iDegree+3) + P0    *AuxPoly2(iDegree+3)
          P0     = X*Pol*AuxPoly1(iDegree+4) + P1    *AuxPoly2(iDegree+4)
          P1     = X*P0    *AuxPoly1(iDegree+5) + Pol*AuxPoly2(iDegree+5)
       END DO
       XC = P1 * (1.0_r16-X*X) / (REAL(degree,r16) * (P0-X*P1))
       X = X - XC
       IF (MAXVAL(ABS(XC/X)).LT.rootPrecision) THEN
          GaussWeights(1:halfDegree) = REAL(scale * (1.0_r16-X*X) / (P0*P0), KIND=r8)
          GaussWeights(halfDegree+1:Degree)=GaussWeights(halfDegree:1:-1)
          EXIT
       END IF
    END DO
    !
    !
    CosGaussColat = REAL(X,KIND=r8)
    !
  END SUBROUTINE LegPolRootsandWeights

 !
 ! maps (ib,jb) into (i,j)
 !
 SUBROUTINE IBJBtoIJ_R(var_in,var_out)
    REAL(KIND=r8)   , INTENT(IN   ) :: var_in (:,:)
    REAL(KIND=r8)   , INTENT(OUT  ) :: var_out(:,:)
    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: ib
    INTEGER                :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1,jbmax
       DO ib = 1,ibmaxPerJB(jb)
          i = iPerIJB(ib,jb)
          j = jPerIJB(ib,jb)-myfirstlat+1
          var_out(i,j)=var_in(ib,jb)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IBJBtoIJ_R
 !
 ! maps (i,j) into (ib,jb)
 !
 SUBROUTINE IJtoIBJB_R(var_in,var_out)
    REAL(KIND=r8)   , INTENT(IN  ) :: var_in (:,:)
    REAL(KIND=r8)   , INTENT(OUT ) :: var_out(:,:)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1,ibMaxPerJB(jb)
         i = iPerIJB(ib,jb)
         j = jPerIJB(ib,jb)
         var_out(ib,jb)=var_in(i,j)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB_R

 ! 3D version by hmjb
 SUBROUTINE IJtoIBJB3_R(var_in,var_out)
    REAL(KIND=r8)   , INTENT(IN  ) :: var_in (iMax,kMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: var_out(ibMax,kMax,jbMax)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1,ibMaxPerJB(jb)
         i = iPerIJB(ib,jb)
         j = jPerIJB(ib,jb)
         var_out(ib,:,jb)=var_in(i,:,j)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB3_R

 SUBROUTINE IBJBtoIJ_I(var_in,var_out)
    INTEGER, INTENT(IN   ) :: var_in (:,:)
    INTEGER, INTENT(OUT  ) :: var_out(:,:)
    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: ib
    INTEGER                :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1,jbmax
       DO ib = 1,ibmaxPerJB(jb)
          i = iPerIJB(ib,jb)
          j = jPerIJB(ib,jb)-myfirstlat+1
          var_out(i,j)=var_in(ib,jb)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IBJBtoIJ_I
 !
 ! maps (i,j) into (ib,jb)
 !
 SUBROUTINE IJtoIBJB_I(var_in,var_out)
    INTEGER, INTENT(IN  ) :: var_in (:,:)
    INTEGER, INTENT(OUT ) :: var_out(:,:)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1, ibMaxPerJB(jb)
         i = iPerIJB(ib,jb)
         j = jPerIJB(ib,jb)
         var_out(ib,jb)=var_in(i,j)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB_I
 ! 3D version by hmjb
 SUBROUTINE IJtoIBJB3_I(var_in,var_out)
    INTEGER, INTENT(IN  ) :: var_in (iMax,kMax,jMax)
    INTEGER, INTENT(OUT ) :: var_out(ibMax,kMax,jbMax)
    INTEGER               :: i
    INTEGER               :: j
    INTEGER               :: ib
    INTEGER               :: jb
   !$OMP PARALLEL DO PRIVATE(i,j,ib)
    DO jb = 1, jbMax
      DO ib = 1, ibMaxPerJB(jb)
         i = iPerIJB(ib,jb)
         j = jPerIJB(ib,jb)
         var_out(ib,:,jb)=var_in(i,:,j)
       END DO
    END DO
   !$OMP END PARALLEL DO
 END SUBROUTINE IJtoIBJB3_I



  SUBROUTINE SplineIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8) :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast 
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**SplineIJtoIBJB**"

    PRINT *, h

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,ib,jb,i,fout)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicCubicSpline(iMax, iMaxPerJ(j), &
               FieldIn(1,j), FOut, ifirst, ilast)
          DO i=ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE SplineIJtoIBJB_R2D




  SUBROUTINE SplineIBJBtoIJ_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (ibMax,jbMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(iMax,jMax)

    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jl
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**SplineIBJBtoIJ**"

    PRINT *, h
    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
       DO j = myfirstlat,mylastlat
          jl = j-myfirstlat+1
          iFirst = ibPerIJ(1,j)
          jb = jbPerIJ(1,j)
!         CALL CyclicCubicSpline(iMaxPerJ(j), iMax, &
!              FieldIn(iFirst,jb), FieldOut(1,jl))
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)-myfirstlat+1
             FieldOut(i,j)=FieldIn(ib,jb)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE SplineIBJBtoIJ_R2D





  ! CyclicCubicSpline:
  ! Cubic Spline interpolation on cyclic, equally spaced data over [0:2pi].
  ! Input data:
  !   DimIn: How many input points: abcissae are supposed to be x(i),
  !          i=1,...,DimIn, with x(i) = 2*pi*(i-1)/DimIn
  !   DimOut: How many output points: abcissae are supposed to be x(i),
  !          i=1,...,DimOut, with x(i) = 2*pi*(i-1)/DimOut
  !   FieldIn: function values at the DimIn abcissae
  ! Output data:
  !   FieldOut: function values at the DimOut abcissae
  ! Requirements:
  !   DimIn >= 4. Program Stops if DimIn < 4.


  SUBROUTINE CyclicCubicSpline (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    INTEGER, INTENT(IN ) :: DimIn
    INTEGER, INTENT(IN ) :: DimOut
    INTEGER, INTENT(IN ) :: ifirst
    INTEGER, INTENT(IN ) :: ilast
    REAL(KIND=r8),    INTENT(IN ) :: FieldIn(DimIn)
    REAL(KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)


    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8) :: alfa(DimIn+2)
    REAL(KIND=r8) :: beta(DimIn+2)
    REAL(KIND=r8) :: gama(DimIn+2)
    REAL(KIND=r8) :: delta(DimIn+2)
    REAL(KIND=r8) :: der(DimIn+2)
    REAL(KIND=r8) :: c1(DimIn+2)
    REAL(KIND=r8) :: c2(DimIn+2)
    REAL(KIND=r8) :: ratio
    REAL(KIND=r8) :: dxm
    REAL(KIND=r8) :: dxm2
    REAL(KIND=r8) :: dx
    REAL(KIND=r8) :: dx2
    REAL(KIND=r8) :: pi
    REAL(KIND=r8) :: hIn
    REAL(KIND=r8) :: hIn2
    REAL(KIND=r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst,ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 4

    IF (DimIn < 4) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 2) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi   = 4.0_r8*ATAN(1.0_r8)
    hIn  = (2.0_r8*pi)/REAL(DimIn, r8)
    hOut = (2.0_r8*pi)/REAL(DimOut,r8)
    hIn2 = hIn*hIn

    ! tridiagonal system initialization:
    ! indices of alfa, beta, gama and delta are equation numbers,
    ! which are derivatives on a diagonal system

    alfa(2) = 0.0_r8
    beta(2) = 6.0_r8
    gama(2) = 0.0_r8
    delta(2) = (6.0_r8/hIn2)*&
            (FieldIn(1)-2.0_r8*FieldIn(2)+FieldIn(3))
    DO iIn = 3, DimIn-1
       alfa(iIn) = 1.0_r8
       beta(iIn) = 4.0_r8
       gama(iIn) = 1.0_r8
       delta(iIn) = (6.0_r8/hIn2)*&
            (FieldIn(iIn-1)-2.0_r8*FieldIn(iIn)+FieldIn(iIn+1))
    END DO
    alfa(DimIn) = 1.0_r8
    beta(DimIn) = 4.0_r8
    gama(DimIn) = 1.0_r8
    delta(DimIn) = (6.0_r8/hIn2)*&
         (FieldIn(DimIn-1)-2.0_r8*FieldIn(DimIn)+FieldIn(1))
    alfa(DimIn+1) = 0.0_r8
    beta(DimIn+1) = 6.0_r8
    gama(DimIn+1) = 0.0_r8
    delta(DimIn+1) = (6.0_r8/hIn2)*&
         (FieldIn(DimIn)-2.0_r8*FieldIn(1)+FieldIn(2))

    ! backward elimination

    DO iIn = 3, DimIn
       beta(iIn) = beta(iIn) - (gama(iIn-1)*alfa(iIn))/beta(iIn-1)
       delta(iIn) = delta(iIn) - (delta(iIn-1)*alfa(iIn))/beta(iIn-1)
    END DO

    ! forward substitution

    der(DimIn+1) = delta(DimIn+1)/beta(DimIn+1)
    DO iIn = DimIn, 2, -1
       der(iIn) = (delta(iIn) - gama(iIn)*der(iIn+1))/beta(iIn)
    END DO

    der(1) = 2.0_r8*der(2) - der(3)
    der(DimIn+2) = 2.0_r8*der(DimIn+1) - der(DimIn)

    ! interpolation

    DO iIn = 1, DimIn+2
       c1(iIn) = der(iIn)/(6.0_r8*hIn)
    END DO
    DO iIn = 1, DimIn
       c2(iIn) = FieldIn(iIn)/hIn - der(iIn)*(hIn/6.0_r8)
    END DO
    c2(DimIn+1) = FieldIn(1)/hIn - der(DimIn+1)*(hIn/6.0_r8)
    c2(DimIn+2) = FieldIn(2)/hIn - der(DimIn+2)*(hIn/6.0_r8)
    DO iOut = ifirst,ilast
       iIn = INT(REAL(iOut-1,r8)*ratio) + 2
       dxm = REAL(iOut-1,r8)*hOut - REAL(iIn-2,r8)*hIn
       dxm2 = dxm*dxm
       dx = REAL(iIn-1,r8)*hIn - REAL(iOut-1,r8)*hOut
       dx2 = dx*dx
       FieldOut(iOut) = &
            (dxm2*c1(iIn)+c2(iIn))*dxm + &
            (dx2*c1(iIn-1)+c2(iIn-1))*dx
    END DO
  END SUBROUTINE CyclicCubicSpline






  SUBROUTINE LinearIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8) :: Fout(imax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**LinearIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,i,ib,fout,jb)
       DO j = myfirstlat,mylastlat
          ifirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicLinear(iMax, iMaxPerJ(j), &
               FieldIn(1,j), FOut, ifirst,ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE LinearIJtoIBJB_R2D




  SUBROUTINE LinearIBJBtoIJ_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (ibMax,jbMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(iMax,jMax)

    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jl
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**LinearIBJBtoIJ**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
       DO j = myfirstlat,mylastlat
          jl = j-myfirstlat+1
          iFirst = ibPerIJ(1,j)
          jb = jbPerIJ(1,j)
!         CALL CyclicLinear(iMaxPerJ(j), iMax, &
!              FieldIn(iFirst,jb), FieldOut(1,jl))
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)-myfirstlat+1
             FieldOut(i,j)=FieldIn(ib,jb)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE LinearIBJBtoIJ_R2D





  ! CyclicLinear:
  ! Linear interpolation on cyclic, equally spaced data over [0:2pi].
  ! Input data:
  !   DimIn: How many input points: abcissae are supposed to be x(i),
  !          i=1,...,DimIn, with x(i) = 2*pi*(i-1)/DimIn
  !   DimOut: How many output points: abcissae are supposed to be x(i),
  !          i=1,...,DimOut, with x(i) = 2*pi*(i-1)/DimOut
  !   FieldIn: function values at the DimIn abcissae
  ! Output data:
  !   FieldOut: function values at the DimOut abcissae
  ! Requirements:
  !   DimIn >= 1. Program Stops if DimIn < 1.


  SUBROUTINE CyclicLinear (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    INTEGER, INTENT(IN ) :: DimIn
    INTEGER, INTENT(IN ) :: DimOut
    INTEGER, INTENT(IN ) :: ifirst
    INTEGER, INTENT(IN ) :: ilast
    REAL(KIND=r8),    INTENT(IN ) :: FieldIn(DimIn)
    REAL(KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)


    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8) :: c(DimIn+1)
    REAL(KIND=r8) :: ratio
    REAL(KIND=r8) :: dxm
    REAL(KIND=r8) :: dx
    REAL(KIND=r8) :: pi
    REAL(KIND=r8) :: hIn
    REAL(KIND=r8) :: hInInv
    REAL(KIND=r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst,ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi = 4.0_r8*ATAN(1.0_r8)
    hIn  = (2.0_r8*pi)/REAL(DimIn, r8)
    hInInv = 1.0_r8/hIn
    hOut = (2.0_r8*pi)/REAL(DimOut,r8)

    ! interpolation

    DO iIn = 1, DimIn
       c(iIn) = FieldIn(iIn)*hInInv
    END DO
    c(DimIn+1) = FieldIn(1)*hInInv

    DO iOut = ifirst,ilast
       iIn = INT(REAL(iOut-1,r8)*ratio) + 2
       dxm = REAL(iOut-1,r8)*hOut - REAL(iIn-2,r8)*hIn
       dx  = REAL(iIn-1,r8)*hIn - REAL(iOut-1,r8)*hOut
       FieldOut(iOut) = dxm*c(iIn)+dx*c(iIn-1)
    END DO
  END SUBROUTINE CyclicLinear

  SUBROUTINE CyclicLinear_ABS (DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    INTEGER, INTENT(IN ) :: DimIn
    INTEGER, INTENT(IN ) :: DimOut
    INTEGER, INTENT(IN ) :: ifirst
    INTEGER, INTENT(IN ) :: ilast
    REAL(KIND=r8),    INTENT(IN ) :: FieldIn(DimIn)
    REAL(KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)


    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8) :: c(DimIn+1)
    REAL(KIND=r8) :: ratio
    REAL(KIND=r8) :: dxm
    REAL(KIND=r8) :: dx
    REAL(KIND=r8) :: pi
    REAL(KIND=r8) :: hIn
    REAL(KIND=r8) :: hInInv
    REAL(KIND=r8) :: hOut

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst,ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear_ABS)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi = 4.0_r8*ATAN(1.0_r8)
    hIn  = (2.0_r8*pi)/REAL(DimIn, r8)
    hInInv = 1.0_r8/hIn
    hOut = (2.0_r8*pi)/REAL(DimOut,r8)

    ! interpolation

    DO iIn = 1, DimIn
       c(iIn) = ABS(FieldIn(iIn))*hInInv
    END DO
    c(DimIn+1) = ABS(FieldIn(1))*hInInv
    DO iOut = ifirst,ilast
       iIn = INT(REAL(iOut-1,r8)*ratio) + 2
       dxm = REAL(iOut-1,r8)*hOut - REAL(iIn-2,r8)*hIn
       dx  = REAL(iIn-1,r8)*hIn - REAL(iOut-1,r8)*hOut
       FieldOut(iOut) = dxm*c(iIn)+dx*c(iIn-1)
    END DO
  END SUBROUTINE CyclicLinear_ABS

  SUBROUTINE NearestIJtoIBJB_I2D(FieldIn,FieldOut)
    INTEGER, INTENT(IN  ) :: FieldIn (iMax,jMax)
    INTEGER, INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    INTEGER            :: Fout(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicNearest_i(iMax, iMaxPerJ(j), &
               FieldIn(1,j), FOut, ifirst, ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE NearestIJtoIBJB_I2D



  SUBROUTINE NearestIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: Fout(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast 
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicNearest_r(iMax, iMaxPerJ(j), &
               FieldIn(1,j), FOut, ifirst, ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE NearestIJtoIBJB_R2D

  ! 3D version by hmjb
  SUBROUTINE NearestIJtoIBJB_I3D(FieldIn,FieldOut)
    INTEGER, INTENT(IN  ) :: FieldIn (iMax,kMax,jMax)
    INTEGER, INTENT(OUT ) :: FieldOut(ibMax,kMax,jbMax)

    INTEGER            :: FOut(iMax)
    INTEGER            :: i,k
    INTEGER            :: ifirst
    INTEGER            :: ilast 
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb,k)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          DO k=1,kMax 
             CALL CyclicNearest_i(iMax, iMaxPerJ(j), &
                  FieldIn(1,k,j), FOut, ifirst, ilast)
             DO i = ifirst,ilast
                ib = ibperij(i,j)
                jb = jbperij(i,j)
                FieldOut(ib,k,jb) = Fout(i)
             END DO
          ENDDO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,:,jb)=FieldIn(i,:,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE NearestIJtoIBJB_I3D


  ! 3D version by hmjb
  SUBROUTINE NearestIJtoIBJB_R3D(FieldIn,FieldOut)
    REAL(KIND=r8), INTENT(IN  ) :: FieldIn (iMax,kMax,jMax)
    REAL(KIND=r8), INTENT(OUT ) :: FieldOut(ibMax,kMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i,k
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb,k)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          DO k=1,kMax 
             CALL CyclicNearest_r(iMax, iMaxPerJ(j), &
                  FieldIn(1,k,j), FOut, ifirst, ilast)
             DO i = ifirst,ilast
                ib = ibperij(i,j)
                jb = jbperij(i,j)
                FieldOut(ib,k,jb) = Fout(i)
             END DO
          ENDDO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,:,jb)=FieldIn(i,:,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE NearestIJtoIBJB_R3D


  SUBROUTINE NearestIBJBtoIJ_I2D(FieldIn,FieldOut)
    INTEGER, INTENT(IN  ) :: FieldIn (ibMax,jbMax)
    INTEGER, INTENT(OUT ) :: FieldOut(iMax,jMax)

    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jl
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIBJBtoIJ**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
       DO j = myfirstlat,mylastlat
          jl = j-myfirstlat+1
          iFirst = ibPerIJ(1,j)
          jb = jbPerIJ(1,j)
!         CALL CyclicNearest_i(iMaxPerJ(j), iMax, &
!              FieldIn(iFirst,jb), FieldOut(1,jl))
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)-myfirstlat+1
             FieldOut(i,j)=FieldIn(ib,jb)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE NearestIBJBtoIJ_I2D


  SUBROUTINE NearestIBJBtoIJ_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(iMax,jMax)

    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jl
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**NearestIBJBtoIJ**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
       DO j = myfirstlat,mylastlat
          jl = j-myfirstlat+1
          iFirst = ibPerIJ(1,j)
          jb = jbPerIJ(1,j)
!         CALL CyclicNearest_r(iMaxPerJ(j), iMax, &
!              FieldIn(iFirst,jb), FieldOut(1,jl))
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)-myfirstlat+1
             FieldOut(i,j)=FieldIn(ib,jb)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE NearestIBJBtoIJ_R2D

 SUBROUTINE CyclicNearest_i(DimIn, DimOut, FieldIn, FieldOut,ifirst,ilast)

    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    INTEGER,    INTENT(IN ) :: FieldIn (DimIn)
    INTEGER,    INTENT(OUT) :: FieldOut(DimOut)

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8) :: ratio
    REAL(KIND=r8) :: pi
    REAL(KIND=r8) :: hIn
    REAL(KIND=r8) :: hInInv
    REAL(KIND=r8) :: hOut
    REAL(KIND=r8) :: difalfa
    REAL(KIND=r8) :: alfaIn
    REAL(KIND=r8) :: alfaOut
    INTEGER :: mplon(DimOut)

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst,ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicNearest)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)
    hIn    = (2.0_r8*pi)/REAL(DimIn, r8)
    hOut   = (2.0_r8*pi)/REAL(DimOut,r8)
    hInInv = 1.0_r8/hIn

    ! interpolation


    DO iOut = ifirst,ilast
       alfaOut=(iOut-1)*hOut
       difalfa=1000E+12
       alfaIn =0.0
       DO iIn = 1, DimIn
          difalfa=min(ABS(alfaIn-alfaOut),difalfa)
          alfaIn=alfaIn+hIn
       END DO
       alfaIn =0.0
       DO iIn = 1, DimIn
          IF (ABS(alfaIn-alfaOut) == difalfa ) THEN
             mplon(iOut) = iIn
             FieldOut(iOut) = FieldIn(mplon(iOut))
          END IF
          alfaIn=alfaIn+hIn
       END DO
    END DO
    RETURN
 END SUBROUTINE CyclicNearest_i



SUBROUTINE CyclicNearest_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)

    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    REAL   (KIND=r8),    INTENT(IN ) :: FieldIn (DimIn)
    REAL   (KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL   (KIND=r8) :: ratio
    REAL   (KIND=r8) :: pi
    REAL   (KIND=r8) :: hIn
    REAL   (KIND=r8) :: hOut
    REAL   (KIND=r8) :: difalfa
    REAL   (KIND=r8) :: difinout
    REAL   (KIND=r8) :: alfaIn
    REAL   (KIND=r8) :: alfaOut
    INTEGER :: lonout

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicNearest)** ERROR: Few input data points"
    END IF

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst, ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)
    hIn    = (2.0_r8*pi)/REAL(DimIn, r8)
    hOut   = (2.0_r8*pi)/REAL(DimOut,r8)

    ! interpolation

    DO iOut = ifirst,ilast
       alfaOut=(iOut-1)*hOut
       alfaIn =0.0
       difalfa = min (alfaOut,2.0_r8*pi-alfaOut)
       lonout = 1
       DO iIn = 2, DimIn
          alfaIN = (iIn-1)*hIn
          difinout = ABS(alfaIn-alfaOut)
          IF (difinout.lt.difalfa) THEN
             difalfa = difinout
             lonout = iIn
          ENDIF
       END DO
       FieldOut(iOut) = FieldIn(lonout)
    END DO
    RETURN
 END SUBROUTINE CyclicNearest_r

 SUBROUTINE CyclicNearest_r2(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)

    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    REAL   (KIND=r8),    INTENT(IN ) :: FieldIn (DimIn)
    REAL   (KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL   (KIND=r8) :: ratio
    REAL   (KIND=r8) :: pi
    REAL   (KIND=r8) :: hIn
    REAL   (KIND=r8) :: hInInv
    REAL   (KIND=r8) :: hOut
    REAL   (KIND=r8) :: difalfa
    REAL   (KIND=r8) :: alfaIn
    REAL   (KIND=r8) :: alfaOut
    INTEGER :: mplon(DimOut)

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst, ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicNearest)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)
    hIn    = (2.0_r8*pi)/REAL(DimIn, r8)
    hOut   = (2.0_r8*pi)/REAL(DimOut,r8)
    hInInv = 1.0_r8/hIn

    ! interpolation

    DO iOut = ifirst,ilast
       alfaOut=(iOut-1)*hOut
       difalfa=1000E+12
       alfaIn =0.0
       DO iIn = 1, DimIn
          difalfa=min(ABS(alfaIn-alfaOut),difalfa)
          alfaIn=alfaIn+hIn
       END DO
       alfaIn =0.0
       DO iIn = 1, DimIn
          IF (ABS(alfaIn-alfaOut) == difalfa ) THEN
             mplon(iOut) = iIn
             FieldOut(iOut) = FieldIn(mplon(iOut))
          END IF
          alfaIn=alfaIn+hIn
       END DO
    END DO
    RETURN
 END SUBROUTINE CyclicNearest_r2
 
 
  SUBROUTINE FreqBoxIJtoIBJB_I2D(FieldIn,FieldOut)
    INTEGER, INTENT(IN  ) :: FieldIn (iMax,jMax)
    INTEGER, INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    INTEGER   :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**FreqBoxIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicFreqBox_i(iMax, iMaxPerJ(j), &
          FieldIn(1:,j), FOut, ifirst, ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE FreqBoxIJtoIBJB_I2D

  SUBROUTINE FreqBoxIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**FreqBoxIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicFreqBox_r(iMax, iMaxPerJ(j), &
          FieldIn(1:,j), FOut, ifirst, ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE FreqBoxIJtoIBJB_R2D






  SUBROUTINE CyclicFreqBox_i(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    IMPLICIT NONE
    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    INTEGER,    INTENT(IN ) :: FieldIn (DimIn)
    INTEGER,    INTENT(OUT) :: FieldOut(DimOut)
    INTEGER,    PARAMETER   :: ncat   = 13 !number of catagories found

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8)    :: ratio
    REAL(KIND=r8)    :: pi
    REAL(KIND=r8)    :: hIn
    REAL(KIND=r8)    :: hOut
    INTEGER :: ici
    INTEGER :: ico
    INTEGER :: ioi
    INTEGER :: ioo
    INTEGER :: i
    INTEGER :: k
    LOGICAL :: flgin  (5)
    LOGICAL :: flgout (5)
    REAL   (KIND=r8) :: dwork  (2*(DimIn+DimOut+2))
    REAL   (KIND=r8) :: wtlon  (   DimIn+DimOut+2 )
    INTEGER :: mplon  (DimIn+DimOut+2,2)
    REAL   (KIND=r8) :: work   (ncat,DimOut)
    REAL   (KIND=r8) :: wrk2   (     DimOut)
    INTEGER :: undef  =0
    REAL   (KIND=r8) :: dof
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: lond
    REAL   (KIND=r8) :: wln
    INTEGER :: lni
    INTEGER :: lno
    INTEGER :: nc
    INTEGER :: mm
    INTEGER :: n
    INTEGER :: nn
    INTEGER :: nd
    INTEGER :: mdist  (7)
    INTEGER :: ndist  (7)
    REAL   (KIND=r8) :: fm
    INTEGER :: nx
    INTEGER :: kl
    REAL   (KIND=r8) :: b      (5)
    REAL   (KIND=r8) :: fr
    REAL   (KIND=r8) :: cmx
    INTEGER :: kmx
    REAL   (KIND=r8) :: fq
    REAL   (KIND=r8) :: fmk
    REAL   (KIND=r8) :: frk
    INTEGER :: iq
    INTEGER :: nq
    INTEGER :: ns
    INTEGER :: nxk
    INTEGER :: klass  (ncat)
    DATA KLASS/6*1,2,2,3,2,3,4,5/

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst, ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicFreqBox)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicFreqBox)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1)=.true.
    flgin (2)=.true.
    flgin (3)=.false.
    flgin (4)=.true.
    flgin (5)=.true.
    flgout(1)=.true.
    flgout(2)=.true.
    flgout(3)=.false.
    flgout(4)=.true.
    flgout(5)=.true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi=DimIn+DimOut+2
    hIn   =(2.0_r8*pi)/REAL(DimIn, r8)
    IF (flgin(5) .OR. flgin(4)) THEN
       ici=0
       dof=0.5_r8
    ELSE
       ici=1
       dof=0.0_r8
    END IF
    DO i=1,DimIn
       dwork(i+ioi)= (dof+DBLE(i-1))*hIn
    END DO
    !
    !     output grid longitudes
    !
    ioo=2*DimIn+DimOut+3
    hOut  =(2.0_r8*pi)/REAL(DimOut,r8)

    IF (flgout(5) .OR. flgout(4)) THEN
       ico=0
       dof=0.5_r8
    ELSE
       ico=1
       dof=0.0_r8
    END IF
    DO i=1,DimOut
       dwork(i+ioo)= (dof+DBLE(i-1))*hOut
    END DO
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1=1
    i2=1
    i3=1
    DO
       IF (dwork(i1+ioi) == dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1) THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF(i1-ici > DimIn/2)mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i2=i2+1
          i3=i3+1
       ELSE IF (dwork(i1+ioi) < dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i3=i3+1
       ELSE
          dwork(i3)=dwork(i2+ioo)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i2=i2+1
          i3=i3+1
       END IF
       IF ((i1 > DimIn) .OR. (i2 > DimOut)) EXIT
    END DO

    IF (i1 > DimIn) i1=1
    IF (i2 > DimOut) i2=1
    DO
       IF (i2 /= 1) THEN
          dwork(i3)=dwork(i2+ioo)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=1
          IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+1
             IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
          END IF
          mplon(i3-1,2)=i2-ico
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+i2-ico
             IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
          END IF
          i2=i2+1
          IF (i2 > DimOut)i2=1
          i3=i3+1
       END IF
       IF (i1 /= 1)THEN
          dwork(i3)=dwork(i1+ioi)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=i1-ici
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+i1-ici
             IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
          END IF
          mplon(i3-1,2)=1
          IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+1
             IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
          END IF
          i1=i1+1
          IF (i1 > DimIn)i1=1
          i3=i3+1
       END IF
       IF (i1 == 1 .AND. i2 == 1) EXIT
    END DO
    wtlon(i3-1)=2.0_r8*pi+dwork(1)-dwork(i3-1)
    mplon(i3-1,1)=1
    IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
    IF (.NOT.flgin(2)) THEN
       mplon(i3-1,1)=DimIn/2+1
       IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
    END IF
    mplon(i3-1,2)=1
    IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
    IF (.NOT.flgout(2)) THEN
       mplon(i3-1,2)=DimOut/2+1
       IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
    END IF
    lond=i3-1

    ! interpolation

    mdist(1:7)        =0
    ndist(1:7)        =0
    FieldOut(ifirst:ilast)   =0.0_r8
    wrk2  (1:DimOut)     =0.0_r8
    work(1:ncat,1:DimOut)=0.0_r8

   DO i=1,lond

      wln=wtlon(i)
      lni=mplon(i,1)
      lno=mplon(i,2)

      IF (FieldIn(lni) == undef)    CYCLE

      nc =FieldIn(lni)

      IF (nc > ncat.or.lno > DimOut) THEN
        WRITE(nfprt,*)nc,lno,i,lni
        STOP 'ERROR IN nc,lno,i and lni point'
      END IF

      IF (nc.lt.1.or.lno.lt.1) THEN
        WRITE(nfprt,*)nc,lno,i,lni
        STOP 'ERROR IN nc,lno,i and lni point'
      END IF
      work(nc ,lno)=work(nc,lno)+wln
      wrk2(lno)=wrk2(lno)+wln
   END DO

   fq=1.0_r8
   nd=0
   ns=0

   DO i=ifirst,ilast

      FieldOut(i)=undef

      IF (wrk2(i) == 0.0_r8) CYCLE

      fm =0.0_r8
      nx =undef
      mm =0
      nn=1
      b(1)=0.0_r8
      b(2)=0.0_r8
      b(3)=0.0_r8
      b(4)=0.0_r8
      b(5)=0.0_r8

        DO n=1,ncat
          fr=work(n,i)/wrk2(i)

     IF (fm < fr) THEN
              fm=fr
              nx=n
            END IF

          kl   =klass(n)
          b(kl)=b(kl)+fr

   IF (fr > 0.5_r8) nn=0
          IF (work(n,i).ne.0.0_r8) mm=mm+1
        END DO

      cmx=0.0_r8
      kmx=0

        DO k=1,5
          IF (b(k) > cmx) THEN
            cmx=b(k)
            kmx=k
          END IF
        END DO

      IF (klass(nx) == kmx) THEN
        FieldOut(i)=nx
        nd=nd+1
        IF (fm.ne.0.0_r8.and.fm < fq) THEN
          fq=fm
          iq=i
          nq=nx
        END IF

      ELSE

        fmk=0.0_r8

        DO n=1,ncat
          IF (klass(n).ne.kmx) CYCLE
            frk=work(n,i)/wrk2(i)
            IF (fmk.lt.frk) THEN
              fmk=frk
              nxk=n
            END IF
        END DO

        FieldOut(i)=nxk
        ns=ns+1
        IF (fmk.ne.0.0_r8.and.fm.lt.fq) THEN
          fq=fmk
          iq=i
          nq=nxk
        END IF
      END IF

      IF (mm.gt.7.and.mm.gt.0)mm=7

      mdist(mm)=mdist(mm)+1
      ndist(mm)=ndist(mm)+nn
   END DO
   RETURN
 END SUBROUTINE CyclicFreqBox_i


 SUBROUTINE CyclicFreqBox_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)
    IMPLICIT NONE
    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    REAL   (KIND=r8),    INTENT(IN ) :: FieldIn (DimIn)
    REAL   (KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)
    INTEGER,    PARAMETER   :: ncat   = 13 !number of catagories found

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8)    :: ratio
    REAL(KIND=r8)    :: pi
    REAL(KIND=r8)    :: hIn
    REAL(KIND=r8)    :: hOut
    INTEGER :: ici
    INTEGER :: ico
    INTEGER :: ioi
    INTEGER :: ioo
    INTEGER :: i
    INTEGER :: k
    LOGICAL :: flgin  (5)
    LOGICAL :: flgout (5)
    REAL   (KIND=r8) :: dwork  (2*(DimIn+DimOut+2))
    REAL   (KIND=r8) :: wtlon  (   DimIn+DimOut+2 )
    INTEGER :: mplon  (DimIn+DimOut+2,2)
    REAL   (KIND=r8) :: work   (ncat,DimOut)
    REAL   (KIND=r8) :: wrk2   (     DimOut)
    REAL   (KIND=r8) :: undef  =0.0
    REAL   (KIND=r8) :: dof
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: lond
    REAL   (KIND=r8) :: wln
    INTEGER :: lni
    INTEGER :: lno
    INTEGER :: nc
    INTEGER :: mm
    INTEGER :: n
    INTEGER :: nn
    INTEGER :: nd
    INTEGER :: mdist  (7)
    INTEGER :: ndist  (7)
    REAL   (KIND=r8) :: fm
    INTEGER :: nx
    INTEGER :: kl
    REAL   (KIND=r8) :: b      (5)
    REAL   (KIND=r8) :: fr
    REAL   (KIND=r8) :: cmx
    INTEGER :: kmx
    REAL   (KIND=r8) :: fq
    REAL   (KIND=r8) :: fmk
    REAL   (KIND=r8) :: frk
    INTEGER :: iq
    INTEGER :: nq
    INTEGER :: ns
    INTEGER :: nxk
    INTEGER :: klass  (ncat)
    DATA KLASS/6*1,2,2,3,2,3,4,5/

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst, ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1)=.true.
    flgin (2)=.true.
    flgin (3)=.false.
    flgin (4)=.true.
    flgin (5)=.true.
    flgout(1)=.true.
    flgout(2)=.true.
    flgout(3)=.false.
    flgout(4)=.true.
    flgout(5)=.true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi=DimIn+DimOut+2
    hIn   =(2.0_r8*pi)/REAL(DimIn, r8)
    IF (flgin(5) .OR. flgin(4)) THEN
       ici=0
       dof=0.5_r8
    ELSE
       ici=1
       dof=0.0_r8
    END IF
    DO i=1,DimIn
       dwork(i+ioi)= (dof+DBLE(i-1))*hIn
    END DO
    !
    !     output grid longitudes
    !
    ioo=2*DimIn+DimOut+3
    hOut  =(2.0_r8*pi)/REAL(DimOut,r8)

    IF (flgout(5) .OR. flgout(4)) THEN
       ico=0
       dof=0.5_r8
    ELSE
       ico=1
       dof=0.0_r8
    END IF
    DO i=1,DimOut
       dwork(i+ioo)= (dof+DBLE(i-1))*hOut
    END DO
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1=1
    i2=1
    i3=1
    DO
       IF (dwork(i1+ioi) == dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1) THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF(i1-ici > DimIn/2)mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i2=i2+1
          i3=i3+1
       ELSE IF (dwork(i1+ioi) < dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i3=i3+1
       ELSE
          dwork(i3)=dwork(i2+ioo)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i2=i2+1
          i3=i3+1
       END IF
       IF ((i1 > DimIn) .OR. (i2 > DimOut)) EXIT
    END DO

    IF (i1 > DimIn) i1=1
    IF (i2 > DimOut) i2=1
    DO
       IF (i2 /= 1) THEN
          dwork(i3)=dwork(i2+ioo)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=1
          IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+1
             IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
          END IF
          mplon(i3-1,2)=i2-ico
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+i2-ico
             IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
          END IF
          i2=i2+1
          IF (i2 > DimOut)i2=1
          i3=i3+1
       END IF
       IF (i1 /= 1)THEN
          dwork(i3)=dwork(i1+ioi)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=i1-ici
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+i1-ici
             IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
          END IF
          mplon(i3-1,2)=1
          IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+1
             IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
          END IF
          i1=i1+1
          IF (i1 > DimIn)i1=1
          i3=i3+1
       END IF
       IF (i1 == 1 .AND. i2 == 1) EXIT
    END DO
    wtlon(i3-1)=2.0_r8*pi+dwork(1)-dwork(i3-1)
    mplon(i3-1,1)=1
    IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
    IF (.NOT.flgin(2)) THEN
       mplon(i3-1,1)=DimIn/2+1
       IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
    END IF
    mplon(i3-1,2)=1
    IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
    IF (.NOT.flgout(2)) THEN
       mplon(i3-1,2)=DimOut/2+1
       IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
    END IF
    lond=i3-1

    ! interpolation

    mdist(1:7)           =0
    ndist(1:7)           =0
    FieldOut(ifirst:ilast)   =0.0_r8
    wrk2  (1:DimOut)     =0.0_r8
    work(1:ncat,1:DimOut)=0.0_r8

   DO i=1,lond

      wln=wtlon(i)
      lni=mplon(i,1)
      lno=mplon(i,2)

      IF (FieldIn(lni) == undef)    CYCLE

      nc =INT(FieldIn(lni))

      IF (nc > ncat.or.lno > DimOut) THEN
        WRITE(nfprt,*)nc,lno,i,lni
        STOP 'ERROR IN nc,lno,i and lni point'
      END IF

      IF (nc.lt.1.or.lno.lt.1) THEN
        WRITE(nfprt,*)nc,lno,i,lni
        STOP 'ERROR IN nc,lno,i and lni point'
      END IF
      work(nc ,lno)=work(nc,lno)+wln
      wrk2(lno)=wrk2(lno)+wln
   END DO

   fq=1.0_r8
   nd=0
   ns=0

   DO i=ifirst,ilast

      FieldOut(i)=undef

      IF (wrk2(i) == 0.0_r8) CYCLE

      fm =0.0_r8
      nx =undef
      mm =0
      nn=1
      b(1)=0.0_r8
      b(2)=0.0_r8
      b(3)=0.0_r8
      b(4)=0.0_r8
      b(5)=0.0_r8

        DO n=1,ncat
          fr=work(n,i)/wrk2(i)

     IF (fm < fr) THEN
              fm=fr
              nx=n
            END IF

          kl   =klass(n)
          b(kl)=b(kl)+fr

   IF (fr > 0.5_r8) nn=0
          IF (work(n,i).ne.0.0_r8) mm=mm+1
        END DO

      cmx=0.0_r8
      kmx=0

        DO k=1,5
          IF (b(k) > cmx) THEN
            cmx=b(k)
            kmx=k
          END IF
        END DO

      IF (klass(nx) == kmx) THEN
        FieldOut(i)=nx
        nd=nd+1
        IF (fm.ne.0.0_r8.and.fm < fq) THEN
          fq=fm
          iq=i
          nq=nx
        END IF

      ELSE

        fmk=0.0_r8

        DO n=1,ncat
          IF (klass(n).ne.kmx) CYCLE
            frk=work(n,i)/wrk2(i)
            IF (fmk.lt.frk) THEN
              fmk=frk
              nxk=n
            END IF
        END DO

        FieldOut(i)=nxk
        ns=ns+1
        IF (fmk.ne.0.0_r8.and.fm.lt.fq) THEN
          fq=fmk
          iq=i
          nq=nxk
        END IF
      END IF

      IF (mm.gt.7.and.mm.gt.0)mm=7

      mdist(mm)=mdist(mm)+1
      ndist(mm)=ndist(mm)+nn
   END DO
   RETURN
 END SUBROUTINE CyclicFreqBox_r




  SUBROUTINE SeaMaskIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**SeaMaskIJtoIBJB**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicSeaMask_r(iMax, iMaxPerJ(j), &
          FieldIn(1,j), FOut, ifirst, ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE SeaMaskIJtoIBJB_R2D



  SUBROUTINE SeaMaskIBJBtoIJ_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(iMax,jMax)

    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jl
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**SeaMaskIBJBtoIJ**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
       DO j = myfirstlat,mylastlat
          jl = j-myfirstlat+1
          iFirst = ibPerIJ(1,j)
          jb = jbPerIJ(1,j)
!         CALL CyclicSeaMask_r(iMaxPerJ(j), iMax, &
!              FieldIn(iFirst,jb), FieldOut(1,jl))
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)-myfirstlat+1
             FieldOut(i,j)=FieldIn(ib,jb)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE SeaMaskIBJBtoIJ_R2D


  SUBROUTINE CyclicSeaMask_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)

    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    REAL   (KIND=r8),    INTENT(IN ) :: FieldIn (DimIn)
    REAL   (KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8)    :: ratio
    REAL(KIND=r8)    :: pi
    REAL(KIND=r8)    :: hIn
    REAL(KIND=r8)    :: hOut
    INTEGER :: ici
    INTEGER :: ico
    INTEGER :: ioi
    INTEGER :: ioo
    INTEGER :: i
    LOGICAL :: flgin  (5)
    LOGICAL :: flgout (5)
    REAL   (KIND=r8) :: dwork  (2*(DimIn+DimOut+2))
    REAL   (KIND=r8) :: wtlon  (   DimIn+DimOut+2 )
    INTEGER :: mplon  (DimIn+DimOut+2,2)
    REAL   (KIND=r8) :: work   (DimOut)
    REAL   (KIND=r8) :: undef  =290.0_r8
    REAL   (KIND=r8) :: dof
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: lond
    REAL   (KIND=r8) :: wln
    INTEGER :: lni
    INTEGER :: lno

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst, ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1)=.true.
    flgin (2)=.true.
    flgin (3)=.false.
    flgin (4)=.true.
    flgin (5)=.true.
    flgout(1)=.true.
    flgout(2)=.true.
    flgout(3)=.false.
    flgout(4)=.true.
    flgout(5)=.true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi=DimIn+DimOut+2
    hIn   =(2.0_r8*pi)/REAL(DimIn, r8)
    IF (flgin(5) .OR. flgin(4)) THEN
       ici=0
       dof=0.5_r8
    ELSE
       ici=1
       dof=0.0_r8
    END IF
    DO i=1,DimIn
       dwork(i+ioi)= (dof+DBLE(i-1))*hIn
    END DO
    !
    !     output grid longitudes
    !
    ioo=2*DimIn+DimOut+3
    hOut  =(2.0_r8*pi)/REAL(DimOut,r8)

    IF (flgout(5) .OR. flgout(4)) THEN
       ico=0
       dof=0.5_r8
    ELSE
       ico=1
       dof=0.0_r8
    END IF
    DO i=1,DimOut
       dwork(i+ioo)= (dof+(i-1))*hOut
    END DO
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1=1
    i2=1
    i3=1
    DO
       IF (dwork(i1+ioi) == dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1) THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF(i1-ici > DimIn/2)mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i2=i2+1
          i3=i3+1
       ELSE IF (dwork(i1+ioi) < dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i3=i3+1
       ELSE
          dwork(i3)=dwork(i2+ioo)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i2=i2+1
          i3=i3+1
       END IF
       IF ((i1 > DimIn) .OR. (i2 > DimOut)) EXIT
    END DO

    IF (i1 > DimIn) i1=1
    IF (i2 > DimOut) i2=1
    DO
       IF (i2 /= 1) THEN
          dwork(i3)=dwork(i2+ioo)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=1
          IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+1
             IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
          END IF
          mplon(i3-1,2)=i2-ico
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+i2-ico
             IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
          END IF
          i2=i2+1
          IF (i2 > DimOut)i2=1
          i3=i3+1
       END IF
       IF (i1 /= 1)THEN
          dwork(i3)=dwork(i1+ioi)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=i1-ici
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+i1-ici
             IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
          END IF
          mplon(i3-1,2)=1
          IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+1
             IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
          END IF
          i1=i1+1
          IF (i1 > DimIn)i1=1
          i3=i3+1
       END IF
       IF (i1 == 1 .AND. i2 == 1) EXIT
    END DO
    wtlon(i3-1)=2.0_r8*pi+dwork(1)-dwork(i3-1)
    mplon(i3-1,1)=1
    IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
    IF (.NOT.flgin(2)) THEN
       mplon(i3-1,1)=DimIn/2+1
       IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
    END IF
    mplon(i3-1,2)=1
    IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
    IF (.NOT.flgout(2)) THEN
       mplon(i3-1,2)=DimOut/2+1
       IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
    END IF
    lond=i3-1

    ! interpolation

    FieldOut(1:DimOut)   =0.0_r8
    work(1:DimOut)=0.0_r8

    DO i=1,lond
      lni=mplon(i,1)
      IF (FieldIn(lni) < 0.0_r8) THEN ! OBS 0.0 Kelvin valor minino de
          ! temperatura usada na interpolacao
         wln = wtlon(i)
         lno = mplon(i,2)
         FieldOut(lno)= FieldOut(lno)+FieldIn(lni)*wln
         work(lno) = work(lno)+wln
      END IF
    END DO

    DO i=ifirst,ilast
      IF (work(i) == 0.0_r8)THEN
         FieldOut(i)=undef
      ELSE
         FieldOut(i)=FieldOut(i)/work(i)
      END IF
    END DO

 END SUBROUTINE CyclicSeaMask_r



 SUBROUTINE AveBoxIJtoIBJB_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (iMax,jMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    REAL(KIND=r8)      :: FOut(iMax)
    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ilast
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**AveBoxIJtoIBJB_R2D**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,ilast,fout,i,ib,jb)
       DO j = myfirstlat,mylastlat
          iFirst = myfirstlon(j)
          ilast  = mylastlon(j)
          CALL CyclicAveBox_r(iMax, iMaxPerJ(j), &
               FieldIn(1,j), FOut, ifirst, ilast)
          DO i = ifirst,ilast
             ib = ibperij(i,j)
             jb = jbperij(i,j)
             FieldOut(ib,jb) = Fout(i)
          END DO
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)
             FieldOut(ib,jb)=FieldIn(i,j)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
 END SUBROUTINE AveBoxIJtoIBJB_R2D



 SUBROUTINE AveBoxIBJBtoIJ_R2D(FieldIn,FieldOut)
    REAL(KIND=r8)   , INTENT(IN  ) :: FieldIn (ibMax,jbMax)
    REAL(KIND=r8)   , INTENT(OUT ) :: FieldOut(iMax,jMax)

    INTEGER            :: i
    INTEGER            :: iFirst
    INTEGER            :: ib
    INTEGER            :: j
    INTEGER            :: jl
    INTEGER            :: jb

    CHARACTER(LEN=*), PARAMETER :: h="**AveBoxIBJBtoIJ_R2D**"

    IF (reducedGrid) THEN
       !$OMP PARALLEL DO PRIVATE(iFirst,jb,jl)
       DO j = myfirstlat,mylastlat
          jl = j-myfirstlat+1
          iFirst = ibPerIJ(1,j)
          jb = jbPerIJ(1,j)
!         CALL CyclicAveBox_r(iMaxPerJ(j), iMax, &
!              FieldIn(iFirst,jb), FieldOut(1,jl))
       END DO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO PRIVATE(ib,i,j)
       DO jb = 1, jbMax
          DO ib = 1, ibMaxPerJB(jb)
             i = iPerIJB(ib,jb)
             j = jPerIJB(ib,jb)-myfirstlat+1
             FieldOut(i,j)=FieldIn(ib,jb)
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF
 END SUBROUTINE AveBoxIBJBtoIJ_R2D


 SUBROUTINE CyclicAveBox_r(DimIn, DimOut, FieldIn, FieldOut, ifirst, ilast)

    INTEGER,    INTENT(IN ) :: DimIn
    INTEGER,    INTENT(IN ) :: DimOut
    INTEGER,    INTENT(IN ) :: ifirst
    INTEGER,    INTENT(IN ) :: ilast
    REAL   (KIND=r8),    INTENT(IN ) :: FieldIn (DimIn)
    REAL   (KIND=r8),    INTENT(OUT) :: FieldOut(DimOut)

    INTEGER :: iIn
    INTEGER :: iOut
    INTEGER :: iRatio
    REAL(KIND=r8)    :: ratio
    REAL(KIND=r8)    :: pi
    REAL(KIND=r8)    :: hIn
    REAL(KIND=r8)    :: hOut
    INTEGER :: ici
    INTEGER :: ico
    INTEGER :: ioi
    INTEGER :: ioo
    INTEGER :: i
    LOGICAL :: flgin  (5)
    LOGICAL :: flgout (5)
    REAL   (KIND=r8) :: dwork  (2*(DimIn+DimOut+2))
    REAL   (KIND=r8) :: wtlon  (   DimIn+DimOut+2 )
    INTEGER :: mplon  (DimIn+DimOut+2,2)
    REAL   (KIND=r8) :: work   (DimOut)
    REAL   (KIND=r8) :: undef  =-999.0_r8
    REAL   (KIND=r8) :: dof
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: lond
    REAL   (KIND=r8) :: wln
    INTEGER :: lni
    INTEGER :: lno

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst, ilast
          iIn =iRatio*(iOut-1) + 1
          FieldOut(iOut) = FieldIn(iIn)
       END DO
       RETURN
    END IF

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear)** ERROR: Few input data points"
    END IF

    ! protection: output data should fit into input data + 2 intervals

    ratio = REAL(DimIn,r8)/REAL(DimOut,r8)
    iIn = INT(REAL(DimOut-1,r8)*ratio) + 2
    IF (iIn > DimIn + 1) THEN
       STOP "**(CyclicCubicSpline)** ERROR: Output data out of input interval"
    END IF

    ! initialization

    pi     = 4.0_r8*ATAN(1.0_r8)

    !
    !     flags: (in or out)
    !     1     start at north pole (true) start at south pole (false)
    !     2     start at prime meridian (true) start at i.d.l. (false)
    !     3     latitudes are at center of box (true)
    !           latitudes are at edge (false) north edge if 1=true
    !                                south edge if 1=false
    !     4     longitudes are at center of box (true)
    !           longitudes are at western edge of box (false)
    !     5     gaussian (true) regular (false)
    !
    flgin (1)=.true.
    flgin (2)=.true.
    flgin (3)=.false.
    flgin (4)=.true.
    flgin (5)=.true.
    flgout(1)=.true.
    flgout(2)=.true.
    flgout(3)=.false.
    flgout(4)=.true.
    flgout(5)=.true.

    !
    !     latitudes done, now do longitudes
    !
    !     input grid longitudes
    !
    ioi=DimIn+DimOut+2
    hIn   =(2.0_r8*pi)/REAL(DimIn, r8)
    IF (flgin(5) .OR. flgin(4)) THEN
       ici=0
       dof=0.5_r8
    ELSE
       ici=1
       dof=0.0_r8
    END IF
    DO i=1,DimIn
       dwork(i+ioi)= (dof+DBLE(i-1))*hIn
    END DO
    !
    !     output grid longitudes
    !
    ioo=2*DimIn+DimOut+3
    hOut  =(2.0_r8*pi)/REAL(DimOut,r8)

    IF (flgout(5) .OR. flgout(4)) THEN
       ico=0
       dof=0.5_r8
    ELSE
       ico=1
       dof=0.0_r8
    END IF
    DO i=1,DimOut
       dwork(i+ioo)= (dof+(i-1))*hOut
    END DO
    !
    !     produce single ordered set of longitudes for both grids
    !     determine longitude weighting and index mapping
    !
    i1=1
    i2=1
    i3=1
    DO
       IF (dwork(i1+ioi) == dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1) THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF(i1-ici > DimIn/2)mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i2=i2+1
          i3=i3+1
       ELSE IF (dwork(i1+ioi) < dwork(i2+ioo)) THEN
          dwork(i3)=dwork(i1+ioi)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i1=i1+1
          i3=i3+1
       ELSE
          dwork(i3)=dwork(i2+ioo)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=dwork(i3)-dwork(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.flgin(2)) THEN
                mplon(i3-1,1)=DimIn/2+i1-ici
                IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.flgout(2)) THEN
                mplon(i3-1,2)=DimOut/2+i2-ico
                IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
             END IF
          END IF
          i2=i2+1
          i3=i3+1
       END IF
       IF ((i1 > DimIn) .OR. (i2 > DimOut)) EXIT
    END DO

    IF (i1 > DimIn) i1=1
    IF (i2 > DimOut) i2=1
    DO
       IF (i2 /= 1) THEN
          dwork(i3)=dwork(i2+ioo)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=1
          IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+1
             IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
          END IF
          mplon(i3-1,2)=i2-ico
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+i2-ico
             IF (i2-ico > DimOut/2) mplon(i3-1,2)=i2-ico-DimOut/2
          END IF
          i2=i2+1
          IF (i2 > DimOut)i2=1
          i3=i3+1
       END IF
       IF (i1 /= 1)THEN
          dwork(i3)=dwork(i1+ioi)
          wtlon(i3-1)=dwork(i3)-dwork(i3-1)
          mplon(i3-1,1)=i1-ici
          IF (.NOT.flgin(2)) THEN
             mplon(i3-1,1)=DimIn/2+i1-ici
             IF (i1-ici > DimIn/2) mplon(i3-1,1)=i1-ici-DimIn/2
          END IF
          mplon(i3-1,2)=1
          IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
          IF (.NOT.flgout(2)) THEN
             mplon(i3-1,2)=DimOut/2+1
             IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
          END IF
          i1=i1+1
          IF (i1 > DimIn)i1=1
          i3=i3+1
       END IF
       IF (i1 == 1 .AND. i2 == 1) EXIT
    END DO
    wtlon(i3-1)=2.0_r8*pi+dwork(1)-dwork(i3-1)
    mplon(i3-1,1)=1
    IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn
    IF (.NOT.flgin(2)) THEN
       mplon(i3-1,1)=DimIn/2+1
       IF (.NOT.(flgin(4) .OR. flgin(5))) mplon(i3-1,1)=DimIn/2
    END IF
    mplon(i3-1,2)=1
    IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut
    IF (.NOT.flgout(2)) THEN
       mplon(i3-1,2)=DimOut/2+1
       IF (.NOT.(flgout(4) .OR. flgout(5))) mplon(i3-1,2)=DimOut/2
    END IF
    lond=i3-1

    ! interpolation

    FieldOut(1:DimOut)   =0.0_r8
    work(1:DimOut)=0.0_r8

    DO i=1,lond
      lni=mplon(i,1)
      IF (FieldIn(lni) /= undef) THEN
         wln = wtlon(i)
         lno = mplon(i,2)
         FieldOut(lno)= FieldOut(lno)+FieldIn(lni)*wln
         work(lno) = work(lno)+wln
      END IF
    END DO

    DO i=ifirst,ilast
      IF (work(i) == 0.0_r8)THEN
         FieldOut(i)=undef
      ELSE
         FieldOut(i)=FieldOut(i)/work(i)
      END IF
    END DO

 END SUBROUTINE CyclicAveBox_r
!
!------------------------------- VFORMAT ----------------------------------
!
 SUBROUTINE vfirec(iunit,a,n,type)

  INTEGER, INTENT(IN)  :: iunit  !#TO deve ser kind default
  INTEGER, INTENT(IN)  :: n
  REAL(KIND=r8), INTENT(OUT)    :: a(n)
  CHARACTER(len=* ), INTENT(IN) :: type
  !
  ! local
  !
  CHARACTER(len=1 ) :: vc(0:63)
  CHARACTER(len=80) :: line
  CHARACTER(len=1 ) :: cs
  INTEGER           :: ich0
  INTEGER           :: ich9
  INTEGER           :: ichcz
  INTEGER           :: ichca
  INTEGER           :: ichla
  INTEGER           :: ichlz
  INTEGER           :: i
  INTEGER           :: nvalline
  INTEGER           :: nchs
  INTEGER           :: ic
  INTEGER           :: ii
  INTEGER           :: isval
  INTEGER           :: iii
  INTEGER           :: ics
  INTEGER           :: nn
  INTEGER           :: nbits
  INTEGER           :: nc
  REAL(KIND=r8)              :: bias
  REAL(KIND=r8)              :: fact
  REAL(KIND=r8)              :: facti
  REAL(KIND=r8)              :: scfct

  vc='0'
  IF (vc(0).ne.'0') CALL vfinit(vc)

  ich0 =ichar('0')
  ich9 =ichar('9')
  ichcz=ichar('Z')
  ichlz=ichar('z')
  ichca=ichar('A')
  ichla=ichar('a')

  READ (iunit,'(2i8,2e20.10)')nn,nbits,bias,fact

  IF (nn.ne.n) THEN
    PRINT*,' Word count mismatch on vfirec record '
    PRINT*,' Words on record - ',nn
    PRINT*,' Words expected  - ',n
    STOP 'vfirec'
  END IF

  nvalline=(78*6)/nbits
  nchs=nbits/6

  DO i=1,n,nvalline
    READ(iunit,'(a78)') line
    ic=0
    DO ii=i,i+nvalline-1
      isval=0
      IF(ii.gt.n) EXIT
      DO iii=1,nchs
         ic=ic+1
         cs=line(ic:ic)
         ics=ichar(cs)
         IF (ics.le.ich9) THEN
            nc=ics-ich0
         ELSE IF (ics.le.ichcz) THEN
            nc=ics-ichca+10
         ELSE
            nc=ics-ichla+36
         END IF
         isval=ior(ishft(nc,6*(nchs-iii)),isval)
      END DO ! loop iii
        a(ii)=isval
    END DO ! loop ii

  END DO ! loop i

  facti=1.0_r8/fact

  IF (type.eq.'LIN') THEN
    DO i=1,n

      a(i)=a(i)*facti-bias

      !print*,'VFM=',i,a(i)
    END DO
  ELSE IF (type.eq.'LOG') THEN
    scfct=2.0_r8**(nbits-1)
    DO i=1,n
        a(i)=sign(1.0_r8,a(i)-scfct)  &
           *(10.0_r8**(abs(20.0_r8*(a(i)/scfct-1.0_r8))-10.0_r8))
    END DO
  END IF
 END SUBROUTINE vfirec
!--------------------------------------------------------
 SUBROUTINE vfinit(vc)
   CHARACTER(len=1), INTENT(OUT  ) :: vc   (*)
   CHARACTER(len=1)                :: vcscr(0:63)
   INTEGER                         :: n

   DATA vcscr/'0','1','2','3','4','5','6','7','8','9'   &
              ,'A','B','C','D','E','F','G','H','I','J'  &
              ,'K','L','M','N','O','P','Q','R','S','T'  &
              ,'U','V','W','X','Y','Z','a','b','c','d'  &
              ,'e','f','g','h','i','j','k','l','m','n'  &
              ,'o','p','q','r','s','t','u','v','w','x'  &
              ,'y','z','{','|'/

  DO n=0,63
      vc(n)=vcscr(n)
  END DO
 END SUBROUTINE vfinit

END MODULE Utils
