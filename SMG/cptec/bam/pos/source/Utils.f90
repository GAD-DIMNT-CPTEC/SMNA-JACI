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

  USE Constants, ONLY : Ndv, CvLHEv, r8, r16, nfprt, pai, twomg, emrad

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
  PUBLIC :: glat
  PUBLIC :: long
  PUBLIC :: cosz
  PUBLIC :: cos2d
  PUBLIC :: NoBankConflict
  PUBLIC :: CreateAssocLegFunc
  PUBLIC :: DestroyAssocLegFunc
  PUBLIC :: Reset_Epslon_To_Local
  PUBLIC :: Epslon
  PUBLIC :: LegFuncS2F
  PUBLIC :: scase
  PUBLIC :: getpoint

  INTERFACE NoBankConflict
     MODULE PROCEDURE NoBankConflictS, NoBankConflictV
  END INTERFACE

  INTERFACE scase
     MODULE PROCEDURE scase2D, scase3D
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
  REAL(KIND=r8), ALLOCATABLE :: glat(:)
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
    ALLOCATE (glat    (jMax))
    ALLOCATE (long    (iMax))
    ALLOCATE (cosz    (jMax))
    ALLOCATE (cos2d   (ibMax,jbMax))

    rd=45.0_r8/ATAN(1.0_r8)
    DO j=1,jMax/2
       lati(j)=90.0_r8-colrad(j)*rd
       lati(jMax-j+1)=-lati(j)
    END DO
    glat=lati
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

  SUBROUTINE scase2D (kflo, nfe, iclcd, gausin)

    IMPLICIT NONE

    INTEGER, PARAMETER :: Ndi=150
    INTEGER, PARAMETER :: Ndp=Ndi+Ndv

    INTEGER, INTENT(IN) :: kflo
    INTEGER, INTENT(IN) :: nfe(Ndp)
    INTEGER, INTENT(IN) :: iclcd(Ndv)

    REAL (KIND=r8), TARGET, INTENT(INOUT) :: gausin(:,:)

    IF (iclcd(nfe(kflo)) == 1) THEN
       gausin=ABS(gausin)
    ELSE IF (iclcd(nfe(kflo)) == 2) THEN
       gausin=gausin/CvLHEv
    ELSE
    END IF

  END SUBROUTINE scase2D

  SUBROUTINE scase3D (kflo, nfe, iclcd, gausin)

    IMPLICIT NONE

    INTEGER, PARAMETER :: Ndi=150
    INTEGER, PARAMETER :: Ndp=Ndi+Ndv

    INTEGER, INTENT(IN) :: kflo
    INTEGER, INTENT(IN) :: nfe(Ndp)
    INTEGER, INTENT(IN) :: iclcd(Ndv)

    REAL (KIND=r8), TARGET, INTENT(INOUT) :: gausin(:,:,:)

    IF (iclcd(nfe(kflo)) == 1) THEN
       gausin=ABS(gausin)
    ELSE IF (iclcd(nfe(kflo)) == 2) THEN
       gausin=gausin/CvLHEv
    ELSE
    END IF

  END SUBROUTINE scase3D

 INTEGER FUNCTION getpoint(coord,locate)
  IMPLICIT NONE
  REAL(KIND=r8), INTENT(IN   ) :: coord(:)
  REAL(KIND=r8), INTENT(IN   ) :: locate
  REAL(KIND=r8)                :: GeomDist(SIZE(coord)) 
  INTEGER             :: npoint(1) 
  GeomDist =ABS(coord-locate)  
  npoint   =MINLOC(GeomDist)
  getpoint =npoint(1)
  !PRINT*,'locate=',locate,'getpoint=',getpoint, 'coord',coord
 END FUNCTION getpoint

END MODULE Utils
