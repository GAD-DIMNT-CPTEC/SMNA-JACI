!
!  $Author: pkubota $
!  $Date: 2011/04/07 16:00:31 $
!  $Revision: 1.18 $
!
MODULE Utils

  ! CreateAssocLegFunc
  ! DestroyAssocLegFunc
  !
  ! DumpAssocLegFunc   ------------------| DumpMatrix(interface)
  !
  ! CreateGaussQuad    ------------------| CreateLegPol
  !
  ! CreateGridValues
  !
  ! DestroyGaussQuad   ------------------| DestroyLegPol
  !
  ! DumpGaussQuad      ------------------| DumpMatrix(interface)
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
  !     DumpAssocLegFunc    Dumps module info
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

  USE Constants, ONLY : &
       i4,         &
       i8,         &
       r8,         &
       pai,        &
       twomg,      &
       er,         &
       r16

  USE Options, ONLY: &
      reducedGrid, nfprt, nscalars, fNameGauss,nAeros, &
      calndr,&
      julday,&
      caldat

  USE IOLowLevel, ONLY : &
      ReadGauss,         &
      WriteGauss

  USE Parallelism, ONLY: &
       myId, &
       MsgOne, &
       FatalError

   IMPLICIT NONE
  SAVE       


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
  PUBLIC :: DumpGaussQuad
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
  PUBLIC :: cel_area
  PUBLIC :: total_mass
  PUBLIC :: total_flux
  PUBLIC :: massconsrv
  PUBLIC :: allpolynomials
  PUBLIC :: fconsrv
  PUBLIC :: fconsrv_flux
  PUBLIC :: totmas
  PUBLIC :: totflux
  PUBLIC :: nTtimes 
  PUBLIC :: aTfluxco2
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
  PUBLIC :: vmax
  PUBLIC :: vaux
  PUBLIC :: vmaxVert
  PUBLIC :: NoBankConflict
  PUBLIC :: DumpMatrix
  PUBLIC :: CreateAssocLegFunc
  PUBLIC :: DestroyAssocLegFunc
  PUBLIC :: DumpAssocLegFunc
  PUBLIC :: Reset_Epslon_To_Local
  PUBLIC :: Epslon
  PUBLIC :: LegFuncS2F
  PUBLIC :: Iminv
  PUBLIC :: Rg
  PUBLIC :: Tql2
  PUBLIC :: Tred2
  PUBLIC :: tmstmp2
  PUBLIC :: InitTimeStamp
  PUBLIC :: TimeStamp
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
  PUBLIC :: dcol
  PUBLIC :: scol
  PUBLIC :: vfirec

  INTERFACE Tql2
     MODULE PROCEDURE Tql2_I, Tql2_I8
  END INTERFACE

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
  INTERFACE DumpMatrix
     MODULE PROCEDURE &
          DumpMatrixReal1D, DumpMatrixReal2D, DumpMatrixReal3D, &
          DumpMatrixInteger1D, DumpMatrixInteger2D, DumpMatrixInteger3D
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
  REAL(KIND=r8), ALLOCATABLE :: cel_area(:)
  REAL(KIND=r8), ALLOCATABLE :: total_mass(:)
  REAL(KIND=r8), ALLOCATABLE :: total_flux(:)
  REAL(KIND=r8), ALLOCATABLE :: massconsrv(:)
  REAL(KIND=r8), ALLOCATABLE :: fconsrv(:,:)
  REAL(KIND=r8), ALLOCATABLE :: fconsrv_flux(:,:)
  REAL(KIND=r8), ALLOCATABLE :: totmas(:)
  REAL(KIND=r8), ALLOCATABLE :: totflux(:)
  INTEGER      , ALLOCATABLE :: nTtimes(:) 
  REAL(KIND=r8), ALLOCATABLE :: aTfluxco2(:)
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
  REAL(KIND=r8), ALLOCATABLE :: vmax(:,:)
  REAL(KIND=r8), ALLOCATABLE :: vaux(:)
  REAL(KIND=r8), ALLOCATABLE :: vmaxVert(:)
  REAL(KIND=r8), ALLOCATABLE :: dcol(:)
  REAL(KIND=r8), ALLOCATABLE :: scol(:)
  !  Module Hided Data:
  !     maxDegree is the degree of the base functions (n)
  !     created specifies if module was created or not



  LOGICAL           :: created=.FALSE.
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



  !DumpAssocLegFunc  Dumps all stored values



  SUBROUTINE DumpAssocLegFunc()
    CHARACTER(LEN=*), PARAMETER :: h="**(DumpAssocLegFunc)**"
    WRITE(nfprt,"(a,' dumping stored values')") h
    CALL DumpMatrix('Epslon', Epslon)
    CALL DumpMatrix('LegFuncS2F',LegFuncS2F)
  END SUBROUTINE DumpAssocLegFunc

  !  AUXILIARY PROCEDURES

  !  Module exports two routines:
  !     NoBankConflict  given input integer (size of an array at
  !                     any dimension) returns the next integer
  !                     that should dimension the array to avoid
  !                     memory bank conflicts. The vector version
  !                     returns a vector of integers, given a vector
  !                     of input integers.
  !     DumpMatrix      Dumps input matrix, for ranks 1 to 3, of type
  !                     integer or real. Output is limited to 10 columns.

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
  SUBROUTINE DumpMatrixReal1D(name, m)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(KIND=r8), INTENT(IN) :: m(:)
    INTEGER :: n1, i1, i1h
    CHARACTER(LEN=10) :: c1
    n1=SIZE(m,1)
    WRITE(c1,"(i10)") n1
    WRITE(nfprt,"(' Dump Matrix ',a,'[',a,']')") &
         TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1))
    i1h = MIN(10,n1)
    WRITE(c1,"(i10)") i1h
    WRITE(nfprt,"(a,'[1:',a,']=',1P,10e9.1)") &
            TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)),&
            (m(i1), i1=1,i1h)
  END SUBROUTINE DumpMatrixReal1D
  SUBROUTINE DumpMatrixReal2D(name, m)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(KIND=r8), INTENT(IN) :: m(:,:)
    INTEGER :: n1, n2, i1, i2, i2h
    CHARACTER(LEN=10) :: c1, c2
    n1=SIZE(m,1); n2=SIZE(m,2)
    WRITE(c1,"(i10)") n1
    WRITE(c2,"(i10)") n2
    WRITE(nfprt,"(' Dump Matrix ',a,'[',a,',',a,']')") &
         TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
         TRIM(ADJUSTL(c2))
    DO i1 = 1, n1
       WRITE(c1,"(i10)") i1
       i2h = MIN(10,n2)
       WRITE(c2,"(i10)") i2h
       WRITE(nfprt,"(a,'[',a,',1:',a,']=',1P,10e9.1)") &
            TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
            TRIM(ADJUSTL(c2)), &
            (m(i1,i2), i2=1,i2h)
    END DO
  END SUBROUTINE DumpMatrixReal2D
  SUBROUTINE DumpMatrixReal3D(name, m)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(KIND=r8), INTENT(IN) :: m(:,:,:)
    INTEGER :: n1, n2, n3, i1, i2, i3, i2h
    CHARACTER(LEN=10) :: c1, c2, c3
    n1=SIZE(m,1); n2=SIZE(m,2); n3=SIZE(m,3)
    WRITE(c1,"(i10)") n1
    WRITE(c2,"(i10)") n2
    WRITE(c3,"(i10)") n3
    WRITE(nfprt,"(' Dump Matrix ',a,'[',a,',',a,',',a,']')") &
         TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
         TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3))
    DO i3 = 1, n3
       WRITE(c3,"(i10)") i3
       PRINT *, ''
       DO i1 = 1, n1
          WRITE(c1,"(i10)") i1
          i2h = MIN(10,n2)
          WRITE(c2,"(i10)") i2h
          WRITE(nfprt,"(a,'[',a,',1:',a,',',a,']=',1P,10e9.1)") &
               TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
               TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3)), &
               (m(i1,i2,i3), i2=1,i2h)
       END DO
    END DO
  END SUBROUTINE DumpMatrixReal3D
  SUBROUTINE DumpMatrixInteger1D(name, m)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: m(:)
    INTEGER :: n1, i1, i1h
    CHARACTER(LEN=10) :: c1
    n1=SIZE(m,1)
    WRITE(c1,"(i10)") n1
    WRITE(nfprt,"(' Dump Matrix ',a,'[',a,']')") &
         TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1))
    i1h = MIN(10,n1)
    WRITE(c1,"(i10)") i1h
    WRITE(nfprt,"(a,'[1:',a,']=',10i8)") &
            TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)),&
            (m(i1), i1=1,i1h)
  END SUBROUTINE DumpMatrixInteger1D
  SUBROUTINE DumpMatrixInteger2D(name, m)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: m(:,:)
    INTEGER :: n1, n2, i1, i2, i2h
    CHARACTER(LEN=10) :: c1, c2
    n1=SIZE(m,1); n2=SIZE(m,2)
    WRITE(c1,"(i10)") n1
    WRITE(c2,"(i10)") n2
    WRITE(nfprt,"(' Dump Matrix ',a,'[',a,',',a,']')") &
         TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
         TRIM(ADJUSTL(c2))
    DO i1 = 1, n1
       WRITE(c1,"(i10)") i1
       i2h = MIN(10,n2)
       WRITE(c2,"(i10)") i2h
       WRITE(nfprt,"(a,'[',a,',1:',a,']=',10i8)") &
            TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
            TRIM(ADJUSTL(c2)), &
            (m(i1,i2), i2=1,i2h)
    END DO
  END SUBROUTINE DumpMatrixInteger2D
  SUBROUTINE DumpMatrixInteger3D(name, m)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: m(:,:,:)
    INTEGER :: n1, n2, n3, i1, i2, i3, i2h
    CHARACTER(LEN=10) :: c1, c2, c3
    n1=SIZE(m,1); n2=SIZE(m,2); n3=SIZE(m,3)
    WRITE(c1,"(i10)") n1
    WRITE(c2,"(i10)") n2
    WRITE(c3,"(i10)") n3
    WRITE(nfprt,"(' Dump Matrix ',a,'[',a,',',a,',',a,']')") &
         TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
         TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3))
    DO i3 = 1, n3
       WRITE(c3,"(i10)") i3
       PRINT *, ''
       DO i1 = 1, n1
          WRITE(c1,"(i10)") i1
          i2h = MIN(10,n2)
          WRITE(c2,"(i10)") i2h
          WRITE(nfprt,"(a,'[',a,',1:',a,',',a,']=',10i8)") &
               TRIM(ADJUSTL(name)), TRIM(ADJUSTL(c1)), &
               TRIM(ADJUSTL(c2)), TRIM(ADJUSTL(c3)), &
               (m(i1,i2,i3), i2=1,i2h)
       END DO
    END DO
  END SUBROUTINE DumpMatrixInteger3D


  !  GAUSSIAN POINTS AND WEIGHTS FOR QUADRATURE
  !  OVER LEGENDRE POLINOMIALS BASE FUNCTIONS

  !  Module exports three routines:
  !     CreateGaussQuad   initializes module
  !     DestroyGaussQuad  destroys module
  !     DumpGaussQuade    dumps Gaussian Quadrature

  !  Module export two arrays:
  !     GaussPoints and GaussWeights
  !
  !  Module uses Module LegPol (Legandre Polinomials);
  !  transparently to the user, that does not have
  !  to create and/or destroy LegPol






  !CreateGaussQuad  computes and hides 'degreeGiven' Gaussian Points
  !                 and Weights over interval [-1:1];
  !                 creates and uses Mod LegPol.



  SUBROUTINE CreateGaussQuad (degreeGiven, gaussGiven)
    INTEGER, INTENT(IN) :: degreeGiven
    LOGICAL, INTENT(IN) :: gaussGiven
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

    IF (gaussGiven) THEN
       Call ReadGauss(CosGaussColat,GaussWeights,fNamegauss,maxDegree)
    ELSE
       CALL LegPolRootsandWeights(maxDegree)
       GaussWeights(maxDegree/2+1:maxDegree) = GaussWeights(maxDegree/2:1:-1)
       IF (myId == 0) THEN
          Call WriteGauss(CosGaussColat,GaussWeights,fNamegauss,maxDegree)
       END IF
    END IF
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

    WRITE (line, FMT='(3F16.12)') &
         MINVAL(ABS(colrad(1:(maxDegree) -1)-colrad(2:(maxDegree)))*180.0_r16/pai), &
         MAXVAL(ABS(colrad(1:(maxDegree) -1)-colrad(2:(maxDegree)))*180.0_r16/pai), &
         180.0_r16/REAL(maxDegree,r16)
    CALL MsgOne(h,line)

    WRITE (line, FMT='(3F16.3)') &
         112.0_r16*MINVAL(ABS(colrad(1:(maxDegree)-1)-colrad(2:(maxDegree)))*180.0_r16/pai), &
         112.0_r16*MAXVAL(ABS(colrad(1:(maxDegree)-1)-colrad(2:(maxDegree)))*180.0_r16/pai), &
         112.0_r16*(180.0_r16/REAL(maxDegree,r16))
    CALL MsgOne(h,line)

    WRITE (line, FMT='(F24.18,1PE24.12)') SUM(GaussWeights(maxDegree/2+1:maxDegree)),&
         1.0_r16-SUM(GaussWeights(maxDegree/2+1:maxDegree))
    CALL MsgOne(h,line)


  END SUBROUTINE CreateGaussQuad


  SUBROUTINE CreateGridValues

    INTEGER :: i, j, ib, jb, jhalf
    REAL(KIND=r8)    :: sinjm(0:jmax/2)
    REAL(KIND=r8)    :: colb(jMax) 
    REAL(KIND=r8)    :: glat(jMax) 
    ALLOCATE (colrad2D(ibMax,jbMax))
    ALLOCATE (rcl     (ibMax,jbMax))
    rcl=0.0_r8
    ALLOCATE (cel_area(jMax))
    ALLOCATE (total_mass(0:nscalars))
    total_mass=0.0_r8
    ALLOCATE (total_flux(0:nAeros))
    total_flux=0.0_r8
    ALLOCATE (massconsrv(jmax))
    ALLOCATE (fconsrv(nscalars,jmax))
    ALLOCATE (fconsrv_flux(nAeros,jmax))
    ALLOCATE (totmas(nscalars));totmas=0.0_r8
    ALLOCATE (totflux(nAeros)) ;totflux=0.0_r8
    ALLOCATE (nTtimes(nAeros)) ;nTtimes=0.0_r8
    ALLOCATE (aTfluxco2(nAeros)) ;aTfluxco2=0.0_r8
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
    ALLOCATE (vmax    (kMax,jbMax))
    ALLOCATE (vaux(kMax))
    ALLOCATE (vmaxVert(kMax))
    ALLOCATE (dcol(jMax)) 
    ALLOCATE (scol(jMax))
    DO j=1,jMax
       lati(j)=colrad(j)
    END DO
    DO i=1,imax
       long(i) = (i-1)*360.0_r8/REAL(iMax,r8)
    ENDDO
    sinjm(0) = 1.0_r8
    DO j=1,jMax/2
       sinjm(j) = cos((lati(j)+lati(j+1))/2.0_r8)
    END DO
    colb(1) = 0.0_r8
    DO j = 2, jMax
       colb(j) = 0.5_r8*(colrad(j)+colrad(j-1))
    END DO
    DO j = 1, jMax-1
       dcol(j) = colb(j+1)-colb(j)
    END DO
    ! theta = 90.0_r8-(180.0_r8/pai)*colrad(i) ! colatitude -> latitude
    ! the 180 degrees are divided into 37 bands with 5deg each
    ! except for the first and last, which have 2.5 deg
    ! The centers of the bands are located at:
    !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
    dcol(jMax) = pai-colb(jMax)
    DO j = 1, jMax
       glat(j) = 90.0_r8 - (180.0_r8/pai)*colrad(j) 
       scol(j) = SIN(colrad(j))
    END DO


    !$OMP PARALLEL DO PRIVATE(jhalf)
    DO j=1,jMax/2
       jhalf = jMax-j+1
       sinlatj (j)=COS(colrad(j))
       coslatj (j)=SIN(colrad(j))
       sinlatj (jhalf)= - sinlatj(j)
       coslatj (jhalf)= coslatj(j)
       cel_area(j) = 0.5_r8 / iMaxPerJ(j) * (sinjm(j-1)-sinjm(j))
       cel_area(jhalf) = cel_area(j)
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
          ercossin(ib,jb)=COS(colrad(j))*rcl(ib,jb)/er
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


  !DumpGaussPoints  dumps gaussian points and weights




  SUBROUTINE DumpGaussQuad()
    CHARACTER(LEN=*), PARAMETER :: h="**(DumpGaussQuad)**"
    CHARACTER(LEN=10) :: c1
    WRITE(c1,"(i10)") maxDegree
    WRITE(nfprt,"(a,' created with maxDegree=',a)") h, TRIM(ADJUSTL(c1))
    CALL DumpMatrix('GaussPoints',GaussPoints)
    CALL DumpMatrix('GaussWeights',GaussWeights)
  END SUBROUTINE DumpGaussQuad
  !
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
       AuxPoly2(i) = REAL(1-i  ,r16)/REAL(i,r16)
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
    WRITE (UNIT=*, FMT='(/,I6,2I5,1P2E16.8)') degree, r16, r16, rootPrecision,ATAN(1.0_r16)/REAL(degree,r16)
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

  !        LINEAR ALGEBRA PROCEDURES
  !
  !        Module exports several routines for matrix inversion
  !        and eigenvalue and eigenvector computations
  !

  !  Module does not require any other module.
  !  Module does not export any value.

SUBROUTINE Iminv (a,n,d,l,m)
  REAL(KIND=r8), INTENT(INOUT) :: a(*)
  REAL(KIND=r8), INTENT(OUT) :: d
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: l(n), m(n)
  !
  !        Iminv computes the inverse of matrix a through a gauss-jordan
  !        algorithm. The output matrix overwrites the input.
  !
  !        ...............................................................
  !
  !        search for largest element
  !
  REAL(KIND=r8) :: biga, hold
  INTEGER :: nk, k, j, iz, i, ij, ki, ji, jk, kj, jr, jq, jp, kk, ik
!cdir novector
  d=1.0_r8
  nk=-n
  DO k=1,n
     nk=nk+n
     l(k)=k
     m(k)=k
     kk=nk+k
     biga=a(kk)
     DO j=k,n
        iz=n*(j-1)
        DO i=k,n
           ij=iz+i
           IF( ABS (biga) <  ABS (a(ij))) THEN
              biga=a(ij)
              l(k)=i
              m(k)=j
           END IF
        END DO
     END DO
     !
     !        interchange rows
     !
     j=l(k)
     IF (j > k) THEN
        ki=k-n
!cdir nodep
        DO i=1,n
           ki=ki+n
           hold=-a(ki)
           ji=ki-k+j
           a(ki)=a(ji)
           a(ji) =hold
        END DO
     END IF
     !
     !        interchange columns
     !
     i=m(k)
     IF (i > k) THEN
        jp=n*(i-1)
!cdir nodep
        DO j=1,n
           jk=nk+j
           ji=jp+j
           hold=-a(jk)
           a(jk)=a(ji)
           a(ji) =hold
        END DO
     END IF
     !
     !        divide column by minus pivot (value of pivot element is
     !        contained in biga)
     !
     IF(biga .EQ. 0.0_r8) THEN
        d=0.0_r8
        RETURN
     END IF
     DO i=1,n
        IF(i .NE. k) THEN
           ik=nk+i
           a(ik)=a(ik)/(-biga)
        END IF
     END DO
     !
     !        reduce matrix
     !
     DO i=1,n
        ik=nk+i
        ij=i-n
!cdir nodep
        DO j=1,n
           ij=ij+n
           IF (i .EQ. k) CYCLE
           IF (j .EQ. k) CYCLE
           kj=ij-i+k
           a(ij)=a(ik)*a(kj)+a(ij)
        END DO
     END DO
     !
     !        divide row by pivot
     !
     kj=k-n
     DO j=1,n
        kj=kj+n
        IF (j .EQ. k) CYCLE
        a(kj)=a(kj)/biga
     END DO
     !
     !        product of pivots
     !
     d=d*biga
     !
     !        replace pivot by reciprocal
     !
     a(kk)=1.0_r8/biga
  END DO
  !
  !        final row and column interchange
  !
  k=n
  DO
     k=(k-1)
     IF (k .LE. 0) RETURN
     i=l(k)
     IF (i > k) THEN
        jq=n*(k-1)
        jr=n*(i-1)
!cdir nodep
        DO j=1,n
           jk=jq+j
           hold=a(jk)
           ji=jr+j
           a(jk)=-a(ji)
           a(ji) =hold
        END DO
     END IF
     j=m(k)
     IF (j .GT. k) THEN
        ki=k-n
!cdir nodep
        DO i=1,n
           ki=ki+n
           hold=a(ki)
           ji=ki-k+j
           a(ki)=-a(ji)
           a(ji) =hold
        END DO
     END IF
  END DO
END SUBROUTINE iminv
SUBROUTINE Balanc(nm,n,a,low,igh,scal)
  !
  !   ** balanc balances a real general matrix, and isolates
  !   **        eigenvalues whenever possible.
  !
  INTEGER, INTENT(IN)  :: nm, n
  REAL(KIND=r8), INTENT(INOUT)  :: a(nm,*)
  REAL(KIND=r8), INTENT(OUT)    :: scal(*)
  INTEGER, INTENT(OUT) :: low, igh
  !
  LOGICAL :: noconv
  INTEGER :: i,j,k,l,m,jj,iexc
  REAL(KIND=r8)    :: c,f,g,r,s,b2,radi
  !
  !   ** radi  is a machine dependent parameter specifying
  !            the base of the machine floating pont representation.
  !
  radi = 2.0_r8
  b2 = radi * radi
  k=1
  l=n
      GOTO 100
  !
  !   ** In-line procedure for row and column exchange.
  !
   20 scal (M)=REAL(J,r8)
      IF (J .EQ. M) GOTO 50
  !
      DO 30 I=1,L
      F=A(I,J)
      A(I,J)=A(I,M)
      A(I,M)=F
   30 CONTINUE
  !
      DO 40 I=K,N
      F=A(J,I)
      A(J,I)=A(M,I)
      A(M,I)=F
   40 CONTINUE
  !
   50 GOTO (80,130) IEXC
  !
  !   ** Search for rows isolating an eigenvalue and push them down.
  !
   80 IF (L .EQ. 1) GOTO 280
      L=L-1
  !
  !   ** For J = L  step -1 until 1 DO -- .
  !
  100 DO 120 JJ=1,L
      J=L+1-JJ
  !
      DO 110 I=1,L
      IF (I .EQ. J) GOTO 110
      IF (A(J,I) .NE. 0.) GOTO 120
  110 CONTINUE
  !
      M=L
      IEXC=1
      GOTO 20
  120 CONTINUE
  !
      GOTO 140
  !
  !   ** Search for columns isolating an eigenvalue and push them left.
  !
  130 K=K+1
  !
  140 DO 170 J=K,L
  !
      DO 150 I=K,L
      IF (I .EQ. J) GOTO 150
      IF (A(I,J) .NE. 0.) GOTO 170
  150 CONTINUE
  !
      M=K
      IEXC=2
      GOTO 20
  170 CONTINUE
  !
  !   ** Now balance the submatrix in rows K to L.
  !
      DO 180 I=K,L
      scal (I)=1.0_r8
  180 CONTINUE
  !
  !   ** Interative loop for norm reduction.
  !
  190 NOCONV=.FALSE.
  !
      DO 270 I=K,L
      C=0.0_r8
      R=0.0_r8
  !
      DO 200 J=K,L
      IF (J .EQ. I) GOTO 200
      C=C+ABS(A(J,I))
      R=R+ABS(A(I,J))
  200 CONTINUE
  !
  !   ** Guard against zero C or R due to underflow.
  !
      IF (C.EQ.0.0_r8 .OR. R.EQ.0.0_r8) GOTO 270
      G=R/radi
      F=1.0_r8
      S=C+R
  210 IF (C .GE. G) GOTO 220
      F=F*radi
      C=C*B2
      GOTO 210
  220 G=R*radi
  230 IF (C .LT.G) GOTO 240
      F=F/radi
      C=C/B2
      GOTO 230
  !
  !   ** Now balance
  !
  240 IF ((C+R)/F .GE. 0.95_r8*S) GOTO 270
      G=1./F
      scal (I)=scal (I)*F
      NOCONV=.TRUE.
  !
      DO 250 J=K,N
      A(I,J)=A(I,J)*G
  250 CONTINUE
  !
      DO 260 J=1,L
      A(J,I)=A(J,I)*F
  260 CONTINUE
  !
  270 CONTINUE
  !
      IF(NOCONV) GOTO 190
  !
  280 LOW=K
      IGH=L
  !

END SUBROUTINE Balanc
SUBROUTINE Balbak(nm,n,low,igh,scal,m,z)
  !
  !   ** BALBAK forms the eigenvectors of a real general matrix
  !             from the eigenvectors of that matrix
  !             transformed by BALANC.
  !
  INTEGER, INTENT(IN) :: nm,n,low,igh,m
  REAL(KIND=r8), INTENT(INOUT) :: z(nm,n)
  REAL(KIND=r8), INTENT(IN) :: scal(n)
  !
  INTEGER :: I,J,K,II
  REAL(KIND=r8) :: S
  !
  IF (M .EQ. 0) GOTO 200
  IF (IGH .EQ. LOW) GOTO 120
  !
  DO 110 I=LOW,IGH
    S=scal (I)
  !
  !   ** Left hand eigenvectors are back transformed
  !      if the foregoing statment is replaced by S = 1.0 / SCALE(I) .
  !
  DO 100 J=1,M
    Z(I,J)=Z(I,J)*S
  100 CONTINUE
  !
  110 CONTINUE
  !
  !   ** For I=LOW-1 step -1 until 1,
  !          IGH+1 step 1 until N DO -- .
  !
  120 DO 140 II=1,N
    I=II
    IF (I.GE.LOW .AND. I.LE.IGH) GOTO 140
    IF (I .LT. LOW) I=LOW-II
    K=INT(scal(I))
    IF (K .EQ. I) GOTO 140
  !
    DO 130 J=1,M
        S=Z(I,J)
        Z(I,J)=Z(K,J)
        Z(K,J)=S
  130 CONTINUE
  !
  140 CONTINUE
  !
  200 RETURN
END SUBROUTINE Balbak
SUBROUTINE Hqr2(nm,n,low,igh,h,wr,wi,z,ierr,matz,machep,tol,*)
  !
  !    ** HQR2 computes the eigenvalues and/or eigenvectors
  !            of a real upper Hessemberg matrix using the
  !            QR method.
  !
  REAL(KIND=r8), INTENT(IN) :: machep, tol
  INTEGER, INTENT(IN) :: nm, n, low, igh, matz
  INTEGER, INTENT(OUT) :: ierr
  REAL(KIND=r8), INTENT(INOUT) :: h(nm,n),z(nm,n)
  REAL(KIND=r8), INTENT(OUT) :: wr(n),wi(n)
  !
  INTEGER :: i,j,k,l,m,en,ll,mm,na,its,mp2,enm2
  !
  REAL(KIND=r8) :: P,Q,R,S,T,X,W,Y,ZZ,NORM
  !
  LOGICAL :: NOTLAS
  !
  !    ** MACHEP is a machine dependent parameter specifying the
  !              relative precision of the floating point arithmetic.
  !              It must be recomputed and replaced for the specific
  !              machine in use.
  !
      IERR=0
      NORM=0.0_r8
      K = 1
  !
  !    ** Store roots isolated by BALANC and compute matrix norm.
  !
      DO 50 I=1,N
  !
      DO 40 J=K,N
      IF (ABS(H(I,J)) .LT. TOL) H(I,J)=TOL
      NORM=NORM+ABS(H(I,J))
   40 CONTINUE
  !
      K=I
      IF (I.GE.LOW .AND. I.LE.IGH) GOTO 50
      WR(I)=H(I,I)
      WI(I)=0.0_r8
   50 CONTINUE
  !
      EN=IGH
      T=0.0_r8
  !
  !    ** Search for next eigenvalues.
  !
   60 IF (EN .LT. LOW) GOTO 340
      ITS=0
      NA=EN-1
      ENM2=NA-1
  !
  !    ** Look for single small sub-diagonal element.
  !       For L=EN step -1 until LOW DO -- .
  !
   70 DO 80 LL=LOW,EN
      L=EN+LOW-LL
      IF (L .EQ. LOW) GOTO 100
      S=ABS(H(L-1,L-1))+ABS(H(L,L))
      IF (S .EQ. 0.0_r8) S=NORM
      IF (ABS(H(L,L-1)) .LE. MACHEP*S) GOTO 100
   80 CONTINUE
  !
  !    ** Form shift.
  !
  100 X=H(EN,EN)
      IF (L .EQ. EN) GOTO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF (L .EQ. NA) GOTO 280
      IF (ITS .EQ. 30) GOTO 1000
      IF (ITS.NE.10 .AND. ITS.NE.20) GOTO 130
  !
  !    Form exceptional shift.
  !
      T=T+X
  !
      DO 120 I=LOW,EN
      H(I,I)=H(I,I)-X
  120 CONTINUE
  !
      S=ABS(H(EN,NA))+ABS(H(NA,ENM2))
      X=0.75_r8*S
      Y=X
      W=-0.4375_r8*S*S
  130 ITS=ITS+1
  !
  !    ** Look for two consecutive small sub-diagonal elements.
  !       For M=EN-2 step -1 until L DO -- .
  !
      DO 140 MM=L,ENM2
      M=ENM2+L-MM
      ZZ=H(M,M)
      R=X-ZZ
      S=Y-ZZ
      P=(R*S-W)/H(M+1,M)+H(M,M+1)
      Q=H(M+1,M+1)-ZZ-R-S
      R=H(M+2,M+1)
      S=ABS(P)+ABS(Q)+ABS(R)
      P=P/S
      Q=Q/S
      R=R/S
      IF (M .EQ. L) GOTO 150
      IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. MACHEP*ABS(P)* &
         (ABS(H(M-1,M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GOTO 150
  140 CONTINUE
  !
  150 MP2=M+2
  !
      DO 160 I=MP2,EN
      H(I,I-2)=0.0_r8
      IF (I .EQ. MP2) GOTO 160
      H(I,I-3)=0.0_r8
  160 CONTINUE
  !
  !    ** Double QR step involving rows L to END and columns M to EN
  !
      DO 260 K=M,NA
      NOTLAS=K .NE. NA
      IF (K .EQ. M) GOTO 170
      P=H(K,K-1)
      Q=H(K+1,K-1)
      R=0.
      IF(NOTLAS) R=H(K+2,K-1)
      X=ABS(P)+ABS(Q)+ABS(R)
      IF (X .EQ. 0.0_r8) GOTO 260
      P=P/X
      Q=Q/X
      R=R/X
  170 S=SIGN(SQRT(P*P+Q*Q+R*R),P)
      IF (K .EQ. M) GOTO 180
      H(K,K-1)=-S*X
      GOTO 190
  180 IF (L .NE. M) H(K,K-1)=-H(K,K-1)
  190 P=P+S
      X=P/S
      Y=Q/S
      ZZ=R/S
      Q=Q/P
      R=R/P
  !
  !    ** Row modification
  !
      DO 210 J=K,N
      P=H(K,J)+Q*H(K+1,J)
      IF (.NOT.NOTLAS) GOTO 200
      P=P+R*H(K+2,J)
      H(K+2,J)=H(K+2,J)-P*ZZ
  200 H(K+1,J)=H(K+1,J)-P*Y
      H(K,J)=H(K,J)-P*X
  210 CONTINUE
  !
      J=MIN(EN,K+3)
  !
  !    ** Column modification
  !
      DO 230 I=1,J
      P=X*H(I,K)+Y*H(I,K+1)
      IF (.NOT.NOTLAS) GOTO 220
      P=P+ZZ*H(I,K+2)
      H(I,K+2)=H(I,K+2)-P*R
  220 H(I,K+1)=H(I,K+1)-P*Q
      H(I,K)=H(I,K)-P
  230 CONTINUE
  !
      IF (MATZ .EQ. 0) GOTO 260
  !
  !    ** Accumulate transformations
  !
      DO 250 I=LOW,IGH
      P=X*Z(I,K)+Y*Z(I,K+1)
      IF (.NOT.NOTLAS) GOTO 240
      P=P+ZZ*Z(I,K+2)
      Z(I,K+2)=Z(I,K+2)-P*R
  240 IF (ABS(P) .LT. TOL) P=TOL
      Z(I,K+1)=Z(I,K+1)-P*Q
      Z(I,K)=Z(I,K)-P
  250 CONTINUE
  !
  260 CONTINUE
  !
      GOTO 70
  !
  !    ** One root found.
  !
  270 H(EN,EN)=X+T
      WR(EN)=H(EN,EN)
      WI(EN)=0.0_r8
      EN=NA
      GOTO 60
  !
  !    ** Two roots found.
  !
  280 P=(Y-X)*0.5_r8
      Q=P*P+W
      ZZ=SQRT(ABS(Q))
      H(EN,EN)=X+T
      X=H(EN,EN)
      H(NA,NA)=Y+T
      IF (Q .LT. 0.0_r8) GOTO 320
  !
  !    ** Real pair.
  !
      ZZ=P+SIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF (ZZ .NE. 0.0_r8) WR(EN)=X-W/ZZ
      WI(NA)=0.0_r8
      WI(EN)=0.0_r8
  !
      IF (MATZ .EQ. 0) GOTO 330
  !
      X=H(EN,NA)
      S=ABS(X)+ABS(ZZ)
      P=X/S
      Q=ZZ/S
      R=SQRT(P*P+Q*Q)
      P=P/R
      Q=Q/R
  !
  !    ** Row  modification.
  !
      DO 290 J=NA,N
      ZZ = H(NA,J)
      H(NA,J)=Q*ZZ+P*H(EN,J)
      H(EN,J)=Q*H(EN,J)-P*ZZ
  290 CONTINUE
  !
  !    ** Column modification.
  !
      DO 300 I=1,EN
      ZZ=H(I,NA)
      H(I,NA)=Q*ZZ+P*H(I,EN)
      H(I,EN)=Q*H(I,EN)-P*ZZ
  300 CONTINUE
  !
  !    ** Accumulate transformations.
  !
      DO 310 I=LOW,IGH
      ZZ=Z(I,NA)
      Z(I,NA)=Q*ZZ+P*Z(I,EN)
      Z(I,EN)=Q*Z(I,EN)-P*ZZ
  310 CONTINUE
  !
      GOTO 330
  !
  !    ** Complex pair
  !
  320 WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
  !
  330 EN=ENM2
  !
      GOTO 60
  !
  340 IF (MATZ .EQ. 0) RETURN
  !
  !    ** All roots found.
  !       Backsubstitute to find vectors of upper triangular form.
  !
      IF (NORM .EQ. 0.0_r8) GOTO 1001
  !
      CALL HQR3(NM,N,LOW,IGH,H,WR,WI,Z,MACHEP,NORM)
  !
      GOTO 1001
  !
  !    ** Set error - no convergence to an
  !       eigenvalue after 30 iterations
  !
 1000 IERR = EN
 1001 RETURN
END SUBROUTINE Hqr2
SUBROUTINE Hqr3(nm,n,low,igh,h,wr,wi,z,machep,norm)
  !
  !    ** HQR3 backsubstitutes to find
  !            vectors of upper triangular form.
  !
  INTEGER, INTENT(IN) ::  nm,n,low,igh
  REAL(KIND=r8), INTENT(INOUT) ::  h(nm,n)
  REAL(KIND=r8), INTENT(IN) ::  wr(n),wi(n),machep,norm
  REAL(KIND=r8), INTENT(OUT) :: z(nm,n)

  !
  INTEGER ::  i,j,k,m,en,ii,jj,na,nn,enm2
  REAL(KIND=r8) :: P,Q,R,S,T,X,W,Y,AR,AI,BR,BI,RA,SA,ZZ
  !
  !    ** For EN=N step -1 until 1 DO -- .
  !
      DO 800 NN=1,N
      EN=N+1-NN
      P=WR(EN)
      Q=WI(EN)
      NA=EN-1
      IF (Q) 710,600,800
  !
  !    ** Real vector.
  !
  600 M=EN
      H(EN,EN)=1.0_r8
      IF (NA .EQ. 0) GOTO 800
  !
  !    ** For I=EN-1 step -1 until 1 DO -- .
  !
      DO 700 II=1,NA
      I=EN-II
      W=H(I,I)-P
      R=H(I,EN)
      IF (M .GT. NA) GOTO 620
  !
      DO 610 J=M,NA
      R=R+H(I,J)*H(J,EN)
  610 CONTINUE
  !
  620 IF (WI(I) .GE. 0.0_r8) GOTO 630
      ZZ=W
      S=R
      GOTO 700
  630 M=I
      IF (WI(I) .NE. 0.0_r8) GOTO 640
      T=W
      IF (W .EQ. 0.0_r8) T=MACHEP*NORM
      H(I,EN)=-R/T
      GOTO 700
  !
  !    ** Solve real equations.
  !
  640 X=H(I,I+1)
      Y=H(I+1,I)
      Q=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
      T=(X*S-ZZ*R)/Q
      H(I,EN)=T
      IF (ABS(X) .LE. ABS(ZZ)) GOTO 650
      H(I+1,EN)=(-R-W*T)/X
      GOTO 700
  650 H(I+1,EN)=(-S-Y*T)/ZZ
  !
  700 CONTINUE
  !
  !    ** End real vector.
  !
      GOTO 800
  !
  !    ** Complex vector.
  !
  710 M=NA
  !
  !    ** Last vector component chosen imaginary so that
  !       eigenvector matrix is triangular.
  !
      IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GOTO 720
      H(NA,NA)=Q/H(EN,NA)
      H(NA,EN)=-(H(EN,EN)-P)/H(EN,NA)
      GOTO 730
  720 H(NA,NA)=ZD(0.0_r8,-H(NA,EN),H(NA,NA)-P,Q)
      H(NA,EN)=ZD(-H(NA,EN),0.0_r8,H(NA,NA)-P,Q)
  730 H(EN,NA)=0.0_r8
      H(EN,EN)=1.0_r8
      ENM2=NA-1
      IF (ENM2 .EQ. 0) GOTO 800
  !
  !    ** For I=EN-2 step -1 until 1 DO -- .
  !
      DO 790 II=1,ENM2
      I=NA-II
      W=H(I,I)-P
      RA=0.0_r8
      SA=H(I,EN)
  !
      DO 760 J=M,NA
      RA=RA+H(I,J)*H(J,NA)
      SA=SA+H(I,J)*H(J,EN)
  760 CONTINUE
  !
      IF (WI(I) .GE. 0.0_r8) GOTO 770
      ZZ=W
      R=RA
      S=SA
      GOTO 790
  770 M=I
      IF (WI(I) .NE. 0.0_r8) GOTO 780
      H(I,NA)=ZD(-RA,-SA,W,Q)
      H(I,EN)=ZD(-SA,RA,W,Q)
      GOTO 790
  !
  !    ** Solve complex equations.
  !
  780 X=H(I,I+1)
      Y=H(I+1,I)
      AR=X*R-ZZ*RA+Q*SA
      AI=X*S-ZZ*SA-Q*RA
      BR=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
      BI=(WR(I)-P)*2.0_r8*Q
      IF (BR.EQ.0.0_r8 .AND. BI.EQ.0.0_r8) BR=MACHEP*NORM*&
         (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
      H(I,NA)=ZD(AR,AI,BR,BI)
      H(I,EN)=ZD(AI,-AR,BR,BI)
      IF (ABS(X) .LE. (ABS(ZZ)+ABS(Q))) GOTO 785
      H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X
      H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X
      GOTO 790
  785 H(I+1,NA)=ZD(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q)
      H(I+1,EN)=ZD(-S-Y*H(I,EN),R+Y*H(I,NA),ZZ,Q)
  790 CONTINUE
  !
  !    ** End complex vector.
  !
  800 CONTINUE
  !
  !    ** End back substitution.
  !       vectors of isolated roots.
  !
      DO 840 I=1,N
      IF (I.GE.LOW .AND. I.LE.IGH) GOTO 840
      DO 820 J=I,N
      Z(I,J)=H(I,J)
  820 CONTINUE
  !
  840 CONTINUE
  !
  !    ** Multiply by transformations matrix to give
  !       vectors of original full matrix.
  !       For J=N step -1 until LOW DO -- .
  !
      DO 880 JJ=LOW,N
      J=N+LOW-JJ
      M=MIN(J,IGH)
  !
      DO 880 I=LOW,IGH
      ZZ=0.0_r8
  !
      DO 860 K=LOW,M
      ZZ=ZZ+Z(I,K)*H(K,J)
  860 CONTINUE
  !
      Z(I,J)=ZZ
  880 CONTINUE
  !
      CONTAINS
         FUNCTION ZD(A1,A2,A3,A4)
            REAL(KIND=r8), INTENT(in) :: A1,A2,A3,A4
            REAL(KIND=r8) :: ZD

            ZD=(A1*A3+A2*A4)/(A3*A3+A4*A4)
         END FUNCTION ZD
END SUBROUTINE Hqr3
SUBROUTINE ortran(nm,n,low,igh,a,ort,z,tolh)
  !
  !   ** ORTRAN accumulates the orthogonal similarity
  !             tranformations used in the reduction of a real
  !             general matrix to upper Hessemberg form
  !             by ORTHES.
  !
  INTEGER, INTENT(IN) ::  nm,n,low,igh
  REAL(KIND=r8), INTENT(IN) ::  a(nm,n), tolh
  REAL(KIND=r8), INTENT(OUT) ::  z(nm,n),ort(n)
  !
  INTEGER ::  i,j,kl,mm,mp,mp1
  REAL(KIND=r8) :: g
  !
  !   ** Initialize Z to identity matrix
  !
      DO 80 I=1,N
      DO 60 J=1, N
      Z(I,J)=0.0_r8
   60 CONTINUE
      Z(I,I)=1.0_r8
   80 CONTINUE
  !
      KL=IGH-LOW-1
      IF (KL .LT. 1) GOTO 200
  !
  !   ** For MP = IGH-1 step -1 until LOW+1 DO --.
  !
      DO 140 MM=1,KL
      MP=IGH-MM
      IF (A(MP,MP-1) .EQ. 0.0_r8) GOTO 140
      MP1=MP+1
  !
      DO 100 I=MP1,IGH
      ORT(I)=A(I,MP-1)
  100 CONTINUE
  !
      DO 130 J=MP,IGH
      G=0.0_r8
  !
      DO 110 I=MP,IGH
      G=G+ORT(I)*Z(I,J)
  110 CONTINUE
  !
  !   ** Divisor below is negative of H formed in ORTHES
  !      double division avoids possible underflow.
  !
      G=(G/ORT(MP))/A(MP,MP-1)
  !
      DO 120 I=MP,IGH
      Z(I,J)=Z(I,J)+G*ORT(I)
      IF (ABS(Z(I,J)) .LT. TOLH) Z(I,J)=TOLH
  120 CONTINUE
  !
  130 CONTINUE
  !
  140 CONTINUE
  !
  200 RETURN
END SUBROUTINE Ortran

SUBROUTINE Orthes(nm,n,low,igh,a,ort,tolh)
  !
  !   ** ORTHES reduces a real general matrix to upper
  !             Hessemberg form using orthogonal similarity.
  !
  INTEGER, INTENT(IN) :: nm, n, low, igh
  REAL(KIND=r8), INTENT(INOUT) :: a(nm,n)
  REAL(KIND=r8), INTENT(OUT) :: ort(igh)
  REAL(KIND=r8), INTENT(IN) :: tolh
  !
  INTEGER :: I,J,M,II,JJ,LA,MP,KP1
  REAL(KIND=r8) :: F,G,H,SCAL
  !
      LA=IGH-1
      KP1=LOW+1
      IF (LA .LT. KP1) GOTO 200
  !
      DO 180 M=KP1,LA
      H=0.0_r8
      ORT(M)=0.0_r8
      SCAL=0.0_r8
  !
  !   ** Scale column.
  !
      DO 90 I=M,IGH
      SCAL=SCAL+ABS(A(I,M-1))
   90 CONTINUE
  !
      IF (SCAL .EQ. 0.0_r8) GOTO 180
      MP=M+IGH
  !
  !   ** For I = IGH step -1 until M DO -- .
  !
      DO 100 II=M,IGH
      I=MP-II
      ORT(I)=A(I,M-1)/SCAL
      H=H+ORT(I)*ORT(I)
  100 CONTINUE
  !
      G=-SIGN(SQRT(H),ORT(M))
      H=H-ORT(M)*G
      ORT(M)=ORT(M)-G
  !
  !   ** Form (I-(U*UT)/H) * A .
  !
      DO 130 J=M,N
      F=0.0_r8
  !
  !   ** For I = IGH step -1 until M DO -- .
  !
      DO 110 II=M,IGH
      I=MP-II
      F=F+ORT(I)*A(I,J)
  110 CONTINUE
  !
      F=F/H
  !
      DO 120 I=M,IGH
      A(I,J)=A(I,J)-F*ORT(I)
      IF (ABS(A(I,J)) .LT. TOLH) A(I,J)=TOLH
  120 CONTINUE
  !
  130 CONTINUE
  !
  !   ** Form (I-(U*UT)/H) * A * (I-(U*UT)/H) .
  !
      DO 160 I=1,IGH
      F=0.0_r8
  !
  !   ** For J = IGH step -1 until M DO -- .
  !
      DO 140 JJ=M,IGH
      J=MP-JJ
      F=F+ORT(J)*A(I,J)
  140 CONTINUE
  !
      F=F/H
  !
      DO 150 J=M,IGH
      A(I,J)=A(I,J)-F*ORT(J)
      IF (ABS(A(I,J)) .LT. TOLH) A(I,J)=TOLH
  150 CONTINUE
  !
  160 CONTINUE
  !
      ORT(M)=SCAL*ORT(M)
      A(M,M-1)=SCAL*G
      IF (ABS(A(M,M-1)) .LT. TOLH) A(M,M-1)=TOLH
  180 CONTINUE
  !
  200 RETURN
END SUBROUTINE Orthes


SUBROUTINE Rg(nm,n,a,wr,wi,matz,z,ierr,eps,scal,ort)
  !
  !   ** RG calculates the eigenvalues and/or eigenvectors
  !         of a real general matrix.
  !
  !   ** Arguments:
  !
  !   ** NM - row dimension of matrix A at the calling routine: Input
  !           Integer variable.
  !
  !   ** N - current dimension of matrix A: Input
  !           Integer variable; must be .LE. NM.
  !
  !   ** A - real matrix (destroyed): Input
  !          Real array with dimensions A(NM,N).
  !
  !   ** WR - real part of the eigenvalues: Output
  !           Real vector with dimensions WR(N).
  !
  !   ** WI - imaginary part of the eigenvalues: Output
  !           Real vector with dimensions WI(N).
  !
  !           OBS: There is nor ordenation for the eigenvalues, except
  !                that for the conjuate complex pairs are put together
  !                and the pair with real positive imaginary part comes
  !                in first place.
  !
  !   ** MATZ - integer variable: Input
  !             Meaning:
  !   ** MATZ = 0 - only eigenvalues non-filtered
  !   ** MATZ = 1 - eigenvalues and eigenvectors normalized and filtered
  !   ** MATZ = 2 - eigenvalues and eigenvectors non-norm. and non-filt.
  !   ** MATZ = 3 - eigenvalues and eigenvectors normalized, non-filt.
  !   **            and without zeroes for .LE. TOLx (see ZNORMA).
  !
  !   ** Z - Eigenvectors: real and imaginary parts, so that:
  !          a) for a real J-th eigenvalue WR(J).NE.0.AND.WI(J).EQ.0,
  !             the J-th eigenvector is (Z(I,J),I=1,N);
  !          b) for a imaginary J-th eigenvalue with WI(J).NE.0,
  !             the (J+1)-th eigenvalue is its conjugate complex;
  !             the J-th eigenvector has real part (Z(I,J),I=1,N) and
  !             imaginary part (Z(I,J+1),I=1,N), and the (J+1)-th
  !             eigenvector has real part (Z(I,J),I=1,N) and
  !             imaginary part (-Z(I,J+1),I=1,N).
  !          Real array with dimensions Z(NM,N): Output
  !
  !   ** IERR - is a integer variable: Output, indicating:
  !             - if N .GT. NM, then the routine RG stop calculations
  !               and returns with IERR=10*N;
  !             - if 30 iteractions is exceeded for the J-th eigenvalue
  !               computation, then the routine RG stop calculations
  !               and returns with IERR=J and the J+1, J+2, ..., N
  !               eigenvalues are computed, but none eigenvector is
  !               computed;
  !             - for a normal termination IERR is set zero.
  !
  !   ** EPS - is a machine dependent parameter specifying the
  !            relative precision of the floating point arithmetic.
  !            It must be recomputed for the specific machine in use.
  !            It is 2**(-20) for 32 bitsand 2**(-50) for 64 bits.
  !            Real variable: Input.
  !
  !   ** TOLH - tolerance value to filter the Hessemberg matrix.
  !             Real variable: Local.
  !
  !   ** TOLW - tolerance value to filter the eigenvalues.
  !             Real variable: Local.
  !
  !   ** TOLZ - tolerance value to filter the eigenvectors.
  !             Real variable: Local.
  !
  !   ** SCALE - working real vector with dimensions SCALE(N).
  !
  !   ** ORT - working real vector with dimensions ORT(N).
  !
  INTEGER, INTENT(IN) ::  nm,n,matz
  REAL(KIND=r8), INTENT(INOUT) :: a(nm,n)
  REAL(KIND=r8), INTENT(OUT) :: z(nm,n),wr(n),wi(n),scal(n),ort(n)
  REAL(KIND=r8), INTENT(IN) :: eps
  INTEGER, INTENT(OUT) :: ierr
  !
  INTEGER :: low,igh
  REAL(KIND=r8) :: TOLH,TOLW,TOLZ
  !
      TOLH=EPS
      TOLW=EPS
      TOLZ=EPS
  !
      IF (N .GT. NM) THEN
      IERR=10*N
      RETURN
      ENDIF
  !
  !   Performing the balance of the input real general matrix
  !   (in place).
  !
      CALL BALANC(NM,N,A,LOW,IGH,SCAL)
  !
  !   Performing the redution of the balanced matrix (in place) to
  !   the Hessemberg superior form. It is used similarity orthogonal
  !   transformations.
  !
      CALL ORTHES(NM,N,LOW,IGH,A,ORT,TOLH)
  !
      IF (MATZ .NE. 0) THEN
  !
  !   Saving the transformations above for eigenvector computations.
  !
      CALL ORTRAN(NM,N,LOW,IGH,A,ORT,Z,TOLH)
      ENDIF
  !
  !   Computing the eigenvalues/eigenvectors of the Hessemberg matrix
  !   using the QR method.
  !
      CALL HQR2(NM,N,LOW,IGH,A,WR,WI,Z,IERR,MATZ,EPS,TOLH,*10)
  !
      IF (IERR .EQ. 0) THEN
  !
  !   Back-transforming the eigenvectors of the Hessembeg matrix to
  !   the eigenvectors of the original input matrix.
  !
      CALL BALBAK(NM,N,LOW,IGH,SCAL,N,Z)
  !
  !   Normalizing and filtering the eigenvectors (See MATZ above and
  !   comments inside ZNORMA routine).
  !
      CALL ZNORMA(NM,N,WR,WI,Z,MATZ,A,TOLW,TOLZ)
      ENDIF
  !
   10 IF (IERR .EQ. 0) RETURN
      WRITE(nfprt,20) IERR
   20 FORMAT(/,1X,'**** The',I4,1X,'-th Eigenvalue Did Not Converge ',&
             '****',//)
  !

END SUBROUTINE Rg
SUBROUTINE Tql2_I(nm,n,d,e,z,eps,ierr)
  !
  !   Abstract: Computes the Eigenvalues and Eigenvectors of a real
  !             Symmetric Tridiagonal Matrix Using the QL Method.
  !
  INTEGER, INTENT(IN) ::  nm
  INTEGER(KIND=i4), INTENT(IN) ::  n
  INTEGER, INTENT(OUT) ::  ierr
  REAL(KIND=r8), INTENT(INOUT) ::  d(n),e(n),z(nm,n)
  REAL(KIND=r8), INTENT(IN) :: eps
  !
  REAL(KIND=r8) :: B,C,F,G,H,P,R,S
  INTEGER :: itm=50
  INTEGER ::  i,j,k,l,m,ii,l1,mml
  !
  !   EPS is a Machine Dependent Parameter Specifying the
  !       Relative Precision of the Floating Point Arithmetic.
  !
  !   EPS = 2.0 ** -50 - 64 BITS
  !   EPS = 2.0 ** -20 - 32 BITS
  !
      IERR=0
  !
      DO 100 I=2,N
      E(I-1)=E(I)
  100 CONTINUE
  !
      F=0.0_r8
      B=0.0_r8
      E(N)=0.0_r8
  !
      DO 240 L=1,N
      J=0
      H=EPS*(ABS(D(L))+ABS(E(L)))
      IF (B .LT. H) B=H
  !
  !   Look for Small Sub-diagonal Element.
  !
      DO 110 M=L,N
      IF (ABS(E(M)) .LE. B) GOTO 120
  !
  !   E(N) is Always Zero, so there is No Exit
  !        Through the Bottom of the Loop.
  !
  110 CONTINUE
  !
  120 IF (M .EQ. L) GOTO 220
  130 IF (J .EQ. ITM) THEN
  !
  !   No Convergence to an Eigenvalue after 50 Iterations.
  !
      IERR=L
      WRITE(nfprt,400) L
  400 FORMAT(/,' *** The',I4,'-th Eigenvalue Did Not Converge ***',/)
      RETURN
      ENDIF
  !
      J=J+1
  !
  !   Form Shift.
  !
      L1=L+1
      G=D(L)
      P=(D(L1)-G)/(2.0_r8*E(L))
      R=SQRT(P*P+1.0_r8)
      D(L)=E(L)/(P+SIGN(R,P))
      H=G-D(L)
  !
      DO 140 I=L1,N
      D(I)=D(I)-H
  140 CONTINUE
  !
      F=F+H
  !
  !   QL Transformation.
  !
      P=D(M)
      C=1.0_r8
      S=0.0_r8
      MML=M-L
  !
  !   For I=M-1 Step -1 Until L DO -- .
  !
      DO 200 II=1,MML
      I=M-II
      G=C*E(I)
      H=C*P
  !
      IF (ABS(P) .LT. ABS(E(I))) THEN
      C=P/E(I)
      R=SQRT(C*C+1.0_r8)
      E(I+1)=S*E(I)*R
      S=1.0_r8/R
      C=C*S
      ELSE
      C=E(I)/P
      R=SQRT(C*C+1.0_r8)
      E(I+1)=S*P*R
      S=C/R
      C=1.0_r8/R
      ENDIF
  !
      P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
  !
  !   Form Vector.
  !
      DO 180 K=1,N
      H=Z(K,I+1)
      Z(K,I+1)=S*Z(K,I)+C*H
      Z(K,I)=C*Z(K,I)-S*H
  180 CONTINUE
  !
  200 CONTINUE
  !
      E(L)=S*P
      D(L)=C*P
      IF (ABS(E(L)) .GT. B) GOTO 130
  220 D(L)=D(L)+F
  240 CONTINUE
  !
  !   Order Eigenvalues and Eigenvectors.
  !
      DO 300 II=2,N
      I=II-1
      K=I
      P=D(I)
  !
      DO 260 J=II,N
      IF (D(J) .LT. P) THEN
      K=J
      P=D(J)
      ENDIF
  260 CONTINUE
  !
      IF (K .NE. I) THEN
      D(K)=D(I)
      D(I)=P
  !
      DO 280 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
      Z(J,K)=P
  280 CONTINUE
      ENDIF
  !
  300 CONTINUE
  !

END SUBROUTINE Tql2_I
SUBROUTINE Tql2_I8(nm,n,d,e,z,eps,ierr)
  !
  !   Abstract: Computes the Eigenvalues and Eigenvectors of a real
  !             Symmetric Tridiagonal Matrix Using the QL Method.
  !
  INTEGER, INTENT(IN) ::  nm
  INTEGER(KIND=i8), INTENT(IN) ::  n
  INTEGER, INTENT(OUT) ::  ierr
  REAL(KIND=r8), INTENT(INOUT) ::  d(n),e(n),z(nm,n)
  REAL(KIND=r8), INTENT(IN) :: eps
  !
  REAL(KIND=r8) :: B,C,F,G,H,P,R,S
  INTEGER :: itm=50
  INTEGER ::  i,j,k,l,m,ii,l1,mml
  !
  !   EPS is a Machine Dependent Parameter Specifying the
  !       Relative Precision of the Floating Point Arithmetic.
  !
  !   EPS = 2.0 ** -50 - 64 BITS
  !   EPS = 2.0 ** -20 - 32 BITS
  !
      IERR=0
  !
      DO 100 I=2,N
      E(I-1)=E(I)
  100 CONTINUE
  !
      F=0.0_r8
      B=0.0_r8
      E(N)=0.0_r8
  !
      DO 240 L=1,N
      J=0
      H=EPS*(ABS(D(L))+ABS(E(L)))
      IF (B .LT. H) B=H
  !
  !   Look for Small Sub-diagonal Element.
  !
      DO 110 M=L,N
      IF (ABS(E(M)) .LE. B) GOTO 120
  !
  !   E(N) is Always Zero, so there is No Exit
  !        Through the Bottom of the Loop.
  !
  110 CONTINUE
  !
  120 IF (M .EQ. L) GOTO 220
  130 IF (J .EQ. ITM) THEN
  !
  !   No Convergence to an Eigenvalue after 50 Iterations.
  !
      IERR=L
      WRITE(nfprt,400) L
  400 FORMAT(/,' *** The',I4,'-th Eigenvalue Did Not Converge ***',/)
      RETURN
      ENDIF
  !
      J=J+1
  !
  !   Form Shift.
  !
      L1=L+1
      G=D(L)
      P=(D(L1)-G)/(2.0_r8*E(L))
      R=SQRT(P*P+1.0_r8)
      D(L)=E(L)/(P+SIGN(R,P))
      H=G-D(L)
  !
      DO 140 I=L1,N
      D(I)=D(I)-H
  140 CONTINUE
  !
      F=F+H
  !
  !   QL Transformation.
  !
      P=D(M)
      C=1.0_r8
      S=0.0_r8
      MML=M-L
  !
  !   For I=M-1 Step -1 Until L DO -- .
  !
      DO 200 II=1,MML
      I=M-II
      G=C*E(I)
      H=C*P
  !
      IF (ABS(P) .LT. ABS(E(I))) THEN
      C=P/E(I)
      R=SQRT(C*C+1.0_r8)
      E(I+1)=S*E(I)*R
      S=1.0_r8/R
      C=C*S
      ELSE
      C=E(I)/P
      R=SQRT(C*C+1.0_r8)
      E(I+1)=S*P*R
      S=C/R
      C=1.0_r8/R
      ENDIF
  !
      P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
  !
  !   Form Vector.
  !
      DO 180 K=1,N
      H=Z(K,I+1)
      Z(K,I+1)=S*Z(K,I)+C*H
      Z(K,I)=C*Z(K,I)-S*H
  180 CONTINUE
  !
  200 CONTINUE
  !
      E(L)=S*P
      D(L)=C*P
      IF (ABS(E(L)) .GT. B) GOTO 130
  220 D(L)=D(L)+F
  240 CONTINUE
  !
  !   Order Eigenvalues and Eigenvectors.
  !
      DO 300 II=2,N
      I=II-1
      K=I
      P=D(I)
  !
      DO 260 J=II,N
      IF (D(J) .LT. P) THEN
      K=J
      P=D(J)
      ENDIF
  260 CONTINUE
  !
      IF (K .NE. I) THEN
      D(K)=D(I)
      D(I)=P
  !
      DO 280 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
      Z(J,K)=P
  280 CONTINUE
      ENDIF
  !
  300 CONTINUE
  !

END SUBROUTINE Tql2_I8
SUBROUTINE Tred2(nm,n,a,d,e,z)
  !
  INTEGER, INTENT(IN) :: nm
  INTEGER(KIND=i8), INTENT(IN) :: n
  REAL(KIND=r8), INTENT(IN) :: a(nm,n)
  REAL(KIND=r8), INTENT(OUT) :: d(n),e(n),z(nm,n)
  !
  REAL(KIND=r8) :: F,G,H,HH,SCAL
  INTEGER ::  I,J,K,L,II,JP1
  !
      DO 100 I=1,N
      DO 100 J=1,I
      Z(I,J)=A(I,J)
  100 CONTINUE
  !
      IF (N .EQ. 1) GOTO 320
  !   for I=N step -1 until 2 do
      DO 300 II=2,N
      I=N+2-II
      L=I-1
      H=0.0_r8
      SCAL=0.0_r8
      IF (L .LT. 2) GOTO 130
  !   scale row (algol tol then not needed)
      DO 120 K=1,L
      SCAL=SCAL+ABS(Z(I,K))
  120 CONTINUE
  !
      IF (SCAL .NE. 0.0_r8 ) GOTO 140
  130 E(I)=Z(I,L)
      GOTO 290
  !
  140 DO 150 K=1,L
      Z(I,K)=Z(I,K)/SCAL
      H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
  !
      F=Z(I,L)
      G=-SIGN(SQRT(H),F)
      E(I)=SCAL*G
      H=H-F*G
      Z(I,L)=F-G
      F=0.0_r8
  !
      DO 240 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.0_r8
  !   form element of A*U
      DO 180 K=1,J
      G=G+Z(J,K)*Z(I,K)
  180 CONTINUE
  !
      JP1=J+1
      IF (L .LT. JP1) GO TO 220
  !
      DO 200 K=JP1,L
      G=G+Z(K,J)*Z(I,K)
  200 CONTINUE
  !   form element of P
  220 E(J)=G/H
      F=F+E(J)*Z(I,J)
  240 CONTINUE
  !
      HH=F/(H+H)
  !   form reduced A
      DO 260 J=1,L
      F=Z(I,J)
      G=E(J)-HH*F
      E(J)=G
  !
      DO 260 K=1,J
      Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
  !
  290 D(I)=H
  300 CONTINUE
  !
  320 D(1)=0.0_r8
      E(1)=0.0_r8
  !   accumulation of transformation matrices
      DO 500 I=1,N
      L=I-1
      IF (D(I) .EQ. 0.0_r8 ) GOTO 380
  !
      DO 360 J=1,L
      G=0.0_r8
  !
      DO 340 K=1,L
      G=G+Z(I,K)*Z(K,J)
  340 CONTINUE
  !
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  !
  380 D(I)=Z(I,I)
      Z(I,I)=1.0_r8
      IF (L .LT. 1) GOTO 500
  !
      DO 400 J=1,L
      Z(I,J)=0.0_r8
      Z(J,I)=0.0_r8
  400 CONTINUE
  !
  500 CONTINUE
  !

END SUBROUTINE Tred2
SUBROUTINE Znorma(nm,n,wr,wi,z,matz,h,tolw,tolz)
  !
  !   ** ZNORMA normalizes and filters the eigenvectors and filters the
  !             eigenvalues.
  !
  !      It sets ZZ = a + b * i, corresponding to the maximum absolute
  !      value of the eigenvector, to:
  !
  !      a) 1.0 + 0.0   * i - if B .EQ. 0.0
  !      b) 1.0 + (b/a) * i - if  ABS(a) .GE. ABS(b)
  !      c) 1.0 + (a/b) * i - if  ABS(a) .LT. ABS(b)
  !
  INTEGER, INTENT(IN) :: nm,n,matz
  REAL(KIND=r8), INTENT(INOUT) :: z(nm,n),wr(n),wi(n)
  REAL(KIND=r8), INTENT(OUT) :: h(nm,n)
  REAL(KIND=r8), INTENT(IN) :: tolw, tolz
  !
  REAL(KIND=r8) :: ZZ,DIV
  INTEGER:: I,J,IC,J1
  !
      IF (MATZ .EQ. 2) RETURN
  !
      DO 900 J=1,N
      DO 900 I=1,N
      H(I,J)=Z(I,J)
  900 CONTINUE
  !
      DO 910 J=1,N
  !
      IF (WI(J) .EQ. 0.0_r8 ) THEN
  !
      ZZ=0.0_r8
      DO 940 I=1,N
      ZZ=MAX(ZZ,ABS(H(I,J)))
      IF (ZZ .EQ. ABS(H(I,J))) IC=I
  940 CONTINUE
  !
      DO 950 I=1,N
      Z(I,J)=H(I,J)/H(IC,J)
  950 CONTINUE
  !
      ELSEIF (WI(J) .GT. 0.0_r8 ) THEN
  !
      ZZ=0.0
      J1=J+1
      DO 960 I=1,N
      DIV=H(I,J)*H(I,J)+H(I,J1)*H(I,J1)
      ZZ=MAX(ZZ,DIV)
      IF (ZZ .EQ. DIV) IC=I
  960 CONTINUE
      IF (ABS(H(IC,J)) .LT. ABS(H(IC,J1))) THEN
      DIV=1.0_r8/H(IC,J1)
      ELSE
      DIV=1.0_r8/H(IC,J)
      ENDIF
      IF (DIV .NE. 0.0 ) THEN
      DO 970 I=1,N
      Z(I,J)=H(I,J)*DIV
      Z(I,J1)=H(I,J1)*DIV
  970 CONTINUE
      ENDIF
  !
      ENDIF
  !
  910 CONTINUE
  !
      IF (MATZ .EQ. 3) RETURN
  !
      DIV=0.0_r8
      DO 980 J=1,N
      ZZ=SQRT(WR(J)*WR(J)+WI(J)*WI(J))
      DIV=MAX(DIV,ZZ)
  980 CONTINUE
      IF (DIV .LE. 0.0_r8 ) DIV=1.0_r8
  !
      DO 990 J=1,N
      IF (ABS(WR(J)/DIV) .LT. TOLW) WR(J)=0.0_r8
      IF (ABS(WI(J)/DIV) .LT. TOLW) WI(J)=0.0_r8
      DO 990 I=1,N
      IF (ABS(Z(I,J)) .LT. TOLZ) Z(I,J)=0.0_r8
  990 CONTINUE
  !
END SUBROUTINE Znorma



  SUBROUTINE tmstmp2(id, ifday, tod, ihr, iday, mon, iyr)
    !
    !
    !==========================================================================
    !    id(4).......date of current data
    !                id(1)....hour(00/12)
    !                id(2)....month
    !                id(3)....day of month
    !                id(4)....year
    !    ifday.......model forecast day
    !    tod.........todx=tod+swint*f3600, model forecast time of
    !                day in seconds
    !                swint....sw subr. call interval in hours
    !                swint has to be less than or equal to trint
    !                              and mod(trint,swint)=0
    !                f3600=3.6e3
    !    ihr.........hour(00/12)
    !    iday........day of month
    !    mon.........month
    !    iyr.........year
    !    yrl.........length of year in days
    !    monl(12)....length of each month in days
    !==========================================================================
    !

    INTEGER, INTENT(in ) :: id(4)
    INTEGER, INTENT(in ) :: ifday
    REAL(KIND=r8),    INTENT(in ) :: tod
    INTEGER, INTENT(out) :: ihr
    INTEGER, INTENT(out) :: iday
    INTEGER, INTENT(out) :: mon
    INTEGER, INTENT(out) :: iyr
    INTEGER          :: ioptn     
    INTEGER          :: day
    INTEGER          :: month
    INTEGER          :: year
    INTEGER          :: LenYearbyDay
    INTEGER          :: kday
    INTEGER          :: idaymn
    REAL(KIND=r8)    :: ctim
    REAL(KIND=r8)    :: hrmodl
    INTEGER          :: monl(12)

    REAL(KIND=r8), PARAMETER :: yrl =   365.2500
    REAL(KIND=r8), PARAMETER ::  ep = .015625
    DATA MONL/31,28,31,30,31,30,&
         31,31,30,31,30,31/


    ctim=tod+id(1)*3600.0_r8

    IF (ctim >= 86400.e0_r8) THEN
       kday=1
       ctim=ctim-86400.e0_r8
    ELSE
       kday=0
    END IF
    !
    !     adjust time to reduce round off error in divsion
    !
    iday = id(3) + ifday + kday
    hrmodl = (ctim+ep)/3600.0_r8
    ihr = INT(hrmodl,KIND=i4)
    mon = id(2)
    iyr = id(4)
    DO
       ioptn=1    
       day  =31
       month  =12
       year=iyr
       CALL calndr (ioptn,  day, month, year, LenYearbyDay)
       idaymn = monl(mon)
       IF (LenYearbyDay == 366 .AND. mon == 2) &
            idaymn=29
       IF (iday <= idaymn) RETURN
       iday = iday - idaymn
       mon = mon + 1
       IF (mon < 13) CYCLE
       mon = 1
       iyr = iyr + 1
    END DO
  END SUBROUTINE tmstmp2

  SUBROUTINE InitTimeStamp(dateinit_s,idate)
    CHARACTER(len=10), INTENT(out) ::  dateinit_s
    INTEGER,           INTENT(in ) :: idate(4)

    !local variables

    INTEGER :: hhi, mmi, ddi, yyyyi

    yyyyi = idate(4)
    mmi   = idate(2)
    ddi   = idate(3)
    hhi   = idate(1)

    dateinit_s='          '
    WRITE(dateinit_s,'(i4.4,3i2.2)') yyyyi, mmi, ddi, hhi

    ! computes the julian day of this calendar date

    juliandayinitintegration = julday(mmi, ddi, yyyyi)
  END SUBROUTINE InitTimeStamp






  SUBROUTINE TimeStamp(datenow_s, idatec, jdt, dt)
    CHARACTER(len=10), INTENT(out  ) :: datenow_s
    INTEGER,           INTENT(inout) :: idatec(4)
    INTEGER,           INTENT(in   )  :: jdt
    REAL(KIND=r8),              INTENT(in   )  :: dt

    INTEGER  :: hhc, mmc, ddc, yyyyc, juliandaynow

    juliandaynow = juliandayinitintegration + (INT(dt)*jdt)/(24*3600)
    CALL caldat(juliandaynow, mmc, ddc, yyyyc)
    hhc = MOD(INT(dt)*jdt/3600,24) +  idatec (1)
    datenow_s='          '
    WRITE(datenow_s,'(i4.4, 3i2.2)' ) yyyyc, mmc, ddc, hhc
    idatec = (/hhc,mmc,ddc,yyyyc /)
  END SUBROUTINE TimeStamp
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
    INTEGER(KIND=i8), INTENT(IN  ) :: var_in (:,:)
    INTEGER(KIND=i8), INTENT(OUT ) :: var_out(:,:)
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
    INTEGER(KIND=i8), INTENT(IN  ) :: var_in (iMax,kMax,jMax)
    INTEGER(KIND=i8), INTENT(OUT ) :: var_out(ibMax,kMax,jbMax)
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
    REAL(KIND=r8) :: c(0:DimIn+2)
    REAL(KIND=r8) :: ratio
    REAL(KIND=r8) :: dxm
    REAL(KIND=r8) :: dx
    REAL(KIND=r8) :: pi
    REAL(KIND=r8) :: hIn
    REAL(KIND=r8) :: hInInv
    REAL(KIND=r8) :: hOut

    FieldOut=0.0_r8;c=0.0_r8;ratio=0.0_r8;dxm=0.0_r8;dx=0.0_r8;pi=0.0_r8;hIn=0.0_r8;hInInv=0.0_r8;hOut=0.0_r8;iIn=0;iOut=0;iRatio=0    

    ! protection: input data size should be at least 1
    IF (DimOut < 1) RETURN
    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear)** ERROR: Few input data points"
    END IF

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst,ilast
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

    pi = 4.0_r8*ATAN(1.0_r8)
    hIn  = (2.0_r8*pi)/REAL(DimIn, r8)
    hOut = (2.0_r8*pi)/REAL(DimOut,r8)
    hInInv = 1.0_r8/hIn

    ! interpolation

    c(0) = FieldIn(DimIn)*hInInv 
    DO iIn = 1, DimIn
       c(iIn) = FieldIn(iIn)*hInInv
    END DO
    c(DimIn+1) = FieldIn(1)*hInInv
    c(DimIn+2) = FieldIn(2)*hInInv

    DO iOut = ifirst,ilast
       IF(iOut <1) cycle
       iIn = INT(REAL(iOut-1,r8)*ratio) + 2
       IF(iIn> DimIn+2 )CYCLE
       IF(iIn>= 1 )THEN
       dxm = REAL(iOut-1,r8)*hOut - REAL( iIn-2,r8)*hIn
       dx  = REAL( iIn-1,r8)*hIn  - REAL(iOut-1,r8)*hOut
       FieldOut(iOut) = dxm*c(iIn)+dx*c(iIn-1)
       END IF
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
    REAL(KIND=r8) :: c(0:DimIn+2)
    REAL(KIND=r8) :: ratio
    REAL(KIND=r8) :: dxm
    REAL(KIND=r8) :: dx
    REAL(KIND=r8) :: pi
    REAL(KIND=r8) :: hIn
    REAL(KIND=r8) :: hInInv
    REAL(KIND=r8) :: hOut

    FieldOut=0.0_r8;c=0.0_r8;ratio=0.0_r8;dxm=0.0_r8;dx=0.0_r8;pi=0.0_r8;hIn=0.0_r8;hInInv=0.0_r8;hOut=0.0_r8;iIn=0;iOut=0;iRatio=0

    ! protection: input data size should be at least 1

    IF (DimIn < 1) THEN
       STOP "**(CyclicLinear_ABS)** ERROR: Few input data points"
    END IF

    ! case every output data abcissae is some input data abcissae

    iRatio = DimIn/DimOut
    IF (iRatio*DimOut == DimIn) THEN
       DO iOut = ifirst,ilast
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

    pi = 4.0_r8*ATAN(1.0_r8)
    hIn  = (2.0_r8*pi)/REAL(DimIn, r8)
    hInInv = 1.0_r8/hIn
    hOut = (2.0_r8*pi)/REAL(DimOut,r8)

    ! interpolation

    c(0) = ABS(FieldIn(DimIn))*hInInv
    DO iIn = 1, DimIn
       c(iIn) = ABS(FieldIn(iIn))*hInInv
    END DO
    c(DimIn+1) = ABS(FieldIn(1))*hInInv
    c(DimIn+2) = ABS(FieldIn(2))*hInInv

    DO iOut = ifirst,ilast
       IF(iOut <1) cycle 
       iIn = INT(REAL(iOut-1,r8)*ratio) + 2
       IF(iIn> DimIn+2 .or. iIn < 1 )CYCLE
       IF(iIn>= 1 )THEN
       dxm = REAL(iOut-1,r8)*hOut - REAL(iIn-2,r8)*hIn
       dx  = REAL(iIn-1,r8)*hIn - REAL(iOut-1,r8)*hOut
       FieldOut(iOut) = dxm*c(iIn)+dx*c(iIn-1)
       END IF
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
    INTEGER(KIND=i8), INTENT(IN  ) :: FieldIn (iMax,jMax)
    INTEGER(KIND=i8), INTENT(OUT ) :: FieldOut(ibMax,jbMax)

    INTEGER(KIND=i8)   :: FOut(iMax)
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
          FieldIn(1:iMax,j), FOut(1:iMax), ifirst, ilast)
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
    INTEGER(KIND=i8),    INTENT(IN ) :: FieldIn (DimIn)
    INTEGER(KIND=i8),    INTENT(OUT) :: FieldOut(DimOut)
    INTEGER,    PARAMETER   :: ncat   = 15 !number of catagories found

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
    REAL   (KIND=r8) :: wtlon  (2*(DimIn+DimOut+2))
    INTEGER          :: mplon  (2*(DimIn+DimOut+2),2)
    REAL   (KIND=r8) :: work   (ncat,DimOut)
    REAL   (KIND=r8) :: wrk2   (     DimOut)
    INTEGER :: undef  
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
    REAL   (KIND=r8) :: b      (0:5)
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
    INTEGER, PARAMETER :: klass  (0:ncat)=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
    undef  =0;mdist=0;ndist=0;b=0.0_r8;ici=0;ico=0;ioi=0;ioo=0
    dwork=0.0_r8;wtlon=0.0_r8;mplon=0;work=0.0_r8;wrk2=0.0_r8;FieldOut=0
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
       ! 0 - > 2pi
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
       IF ((i1 > DimIn) .OR. (i2 > DimOut) .OR. (i3 > (2*(DimIn+DimOut+2)))) EXIT
    END DO

    IF (i1 > DimIn) i1=1
    IF (i2 > DimOut) i2=1
    IF (i3 > (2*(DimIn+DimOut+2)))i3=1

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
    FieldOut(ifirst:ilast)   =0 
    wrk2  (1:DimOut)     =0.0_r8
    work(1:ncat,1:DimOut)=0.0_r8

   DO i=1,lond

      wln=wtlon(i)
      lni=mplon(i,1)
      lno=mplon(i,2)
      IF (lni == 0)    CYCLE
      IF (FieldIn(lni) == undef)    CYCLE

      nc = INT(FieldIn(lni),KIND=i4)

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
      b(0)=0.0_r8
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
    INTEGER,    PARAMETER   :: ncat   = 100 !number of catagories found

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
    REAL   (KIND=r8) :: undef 
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
    REAL   (KIND=r8) :: b      (0:5)
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
    INTEGER, PARAMETER :: klass  (0:ncat)=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
                                            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
                                            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
                                            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
                                            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
    undef  =0.0_r8;mdist=0;ndist=0;b=0.0_r8;ici=0;ico=0;ioi=0;ioo=0
    dwork=0.0_r8;wtlon=0.0_r8;mplon=0;work=0.0_r8;wrk2=0.0_r8;FieldOut=0.0_r8
 
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
      nx =INT(undef)
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
    REAL   (KIND=r8) :: undef
    REAL   (KIND=r8) :: dof
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: lond
    REAL   (KIND=r8) :: wln
    INTEGER :: lni
    INTEGER :: lno
    undef  =290.0_r8
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
    REAL   (KIND=r8) :: wtlon  (2*(DimIn+DimOut+2))
    INTEGER          :: mplon  (2*(DimIn+DimOut+2),2)
    REAL   (KIND=r8) :: work   (DimOut)
    REAL   (KIND=r8) :: undef 
    REAL   (KIND=r8) :: dof
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: i3
    INTEGER :: lond
    REAL   (KIND=r8) :: wln
    INTEGER :: lni
    INTEGER :: lno
    undef  =-999.0_r8
    mplon=0;wtlon=0.0;ico=0;ioi=0;ioo=0
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
       IF ((i1 > DimIn) .OR. (i2 > DimOut) .or. (i3 > (2*(DimIn+DimOut+2))) ) EXIT
    END DO

    IF (i1 > DimIn) i1=1
    IF (i2 > DimOut) i2=1
    IF (i3 > (2*(DimIn+DimOut+2))) i3=1
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
      IF(lni <1)cycle
      IF (FieldIn(lni) /= undef) THEN
         wln = wtlon(i)
         lno = mplon(i,2)
         IF(lno < 1) cycle
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
