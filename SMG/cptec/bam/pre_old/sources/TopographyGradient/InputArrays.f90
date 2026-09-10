MODULE InputArrays

  USE InputParameters, ONLY: r4, r8, Imax, Jmax, Kmax, Mnwv2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GetArrays, ClsArrays

  INTEGER (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                     InitialDate, CurrentDate

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  SigmaInteface, SigmaLayer, qTopoInp

  REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  GrADSOut

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  qTopo, qTopoS, CosLatInv

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  Topo, DTopoDx, DTopoDy

CONTAINS


SUBROUTINE GetArrays

  IMPLICIT NONE

  ALLOCATE (InitialDate(4), CurrentDate(4))
  ALLOCATE (GrADSOut(Imax,Jmax))
  ALLOCATE (SigmaInteface(Kmax+1), SigmaLayer(Kmax), qTopoInp(Mnwv2))
  ALLOCATE (qTopo(Mnwv2), qTopoS(Mnwv2), CosLatInv(Jmax))
  ALLOCATE (Topo(Imax,Jmax), DTopoDx(Imax,Jmax), DTopoDy(Imax,Jmax))

END SUBROUTINE GetArrays


SUBROUTINE ClsArrays

  IMPLICIT NONE

  DEALLOCATE (InitialDate, CurrentDate)
  DEALLOCATE (GrADSOut)
  DEALLOCATE (SigmaInteface, SigmaLayer, qTopoInp)
  DEALLOCATE (qTopo, qTopoS, CosLatInv)
  DEALLOCATE (Topo, DTopoDx, DTopoDy)

END SUBROUTINE ClsArrays


END MODULE InputArrays
