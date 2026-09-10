!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE InputArrays

  USE InputParameters, ONLY: i4, r4, r8, KmaxInp, KmaxInpp, NTracers, &
                             ImaxOut, JmaxOut, KmaxOut, KmaxOutp, &
                             Mnwv2Inp, Mnwv2Out, Mnwv3Inp, Mnwv3Out

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GetArrays, ClsArrays, GetSpHuTracers, ClsSpHuTracers

  INTEGER (KIND=i4), DIMENSION (4), PUBLIC :: DateInitial, DateCurrent

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  DelSInp, SigIInp, SigLInp, SigIOut, SigLOut

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: qWorkInp, qWorkOut

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  DelSigmaInp, SigInterInp, SigLayerInp, &
                  DelSigmaOut, SigInterOut, SigLayerOut

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  qTopoInp, qLnPsInp, qTopoOut, qLnPsOut, qTopoRec

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  qDivgInp, qVortInp, qTvirInp, &
                  qDivgOut, qVortOut, qTvirOut, &
                  qUvelInp, qVvelInp

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: &
                  qSpHuInp, qSpHuOut

  REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: gWorkOut

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  gTopoInp, gTopoOut, gTopoDel, &
                  gLnPsInp, gPsfcInp, gLnPsOut, gPsfcOut

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: &
                  gUvelInp, gVvelInp, gTvirInp, &
                  gDivgInp, gVortInp, gPresInp, &
                  gUvelOut, gVvelOut, gTvirOut, &
                  gPresOut

  REAL (KIND=r8), DIMENSION (:,:,:,:), ALLOCATABLE, PUBLIC :: &
                  gSpHuInp, gSpHuOut

  INTEGER :: Mnwv2Rec

CONTAINS


SUBROUTINE GetArrays

  IMPLICIT NONE

  Mnwv2Rec=MAX(Mnwv2Inp,Mnwv2Out)
  ALLOCATE (qTopoRec (Mnwv2Rec))
  ALLOCATE (qWorkInp (Mnwv2Inp))
  ALLOCATE (qWorkOut (Mnwv2Out))
  ALLOCATE (SigIInp (KmaxInpp), SigLInp (KmaxInp), DelSInp (KmaxInp))
  ALLOCATE (SigIOut (KmaxOutp), SigLOut (KmaxOut))
  ALLOCATE (DelSigmaInp (KmaxInp), SigInterInp (KmaxInpp), SigLayerInp (KmaxInp))
  ALLOCATE (DelSigmaOut (KmaxOut), SigInterOut (KmaxOutp), SigLayerOut (KmaxOut))
  ALLOCATE (qTopoInp (Mnwv2Inp), qLnPsInp (Mnwv2Inp))
  ALLOCATE (qTopoOut (Mnwv2Out), qLnPsOut (Mnwv2Out))
  ALLOCATE (qDivgInp (Mnwv2Inp,KmaxInp), &
            qVortInp (Mnwv2Inp,KmaxInp))
  ALLOCATE (qDivgOut (Mnwv2Out,KmaxOut), &
            qVortOut (Mnwv2Out,KmaxOut))
  ALLOCATE (qUvelInp (Mnwv3Inp,KmaxInp), &
            qVvelInp (Mnwv3Inp,KmaxInp))
  ALLOCATE (qTvirInp (Mnwv2Inp,KmaxInp))
  ALLOCATE (qTvirOut (Mnwv2Out,KmaxOut))
  ALLOCATE (gWorkOut (ImaxOut,JmaxOut))
  ALLOCATE (gTopoInp (ImaxOut,JmaxOut), &
            gTopoOut (ImaxOut,JmaxOut), &
            gTopoDel (ImaxOut,JmaxOut), &
            gLnPsInp (ImaxOut,JmaxOut), &
            gLnPsOut (ImaxOut,JmaxOut), &
            gPsfcInp (ImaxOut,JmaxOut), &
            gPsfcOut (ImaxOut,JmaxOut))
  ALLOCATE (gUvelInp (ImaxOut,JmaxOut,KmaxInp), &
            gVvelInp (ImaxOut,JmaxOut,KmaxInp))
  ALLOCATE (gDivgInp (ImaxOut,JmaxOut,KmaxInp), &
            gVortInp (ImaxOut,JmaxOut,KmaxInp))
  ALLOCATE (gTvirInp (ImaxOut,JmaxOut,KmaxInp), &
            gPresInp (ImaxOut,JmaxOut,KmaxInp))
  ALLOCATE (gUvelOut (ImaxOut,JmaxOut,KmaxOut), &
            gVvelOut (ImaxOut,JmaxOut,KmaxOut))
  ALLOCATE (gTvirOut (ImaxOut,JmaxOut,KmaxOut), &
            gPresOut (ImaxOut,JmaxOut,KmaxOut))

END SUBROUTINE GetArrays


SUBROUTINE ClsArrays

  IMPLICIT NONE

  DEALLOCATE (qTopoRec)
  DEALLOCATE (qWorkInp)
  DEALLOCATE (qWorkOut)
  DEALLOCATE (DelSigmaInp, SigInterInp, SigLayerInp)
  DEALLOCATE (DelSigmaOut, SigInterOut, SigLayerOut)
  DEALLOCATE (qTopoInp, qLnPsInp)
  DEALLOCATE (qTopoOut, qLnPsOut)
  DEALLOCATE (qDivgInp, qVortInp)
  DEALLOCATE (qDivgOut, qVortOut)
  DEALLOCATE (qUvelInp, qVvelInp)
  DEALLOCATE (qTvirInp)
  DEALLOCATE (qTvirOut)
  DEALLOCATE (gWorkOut)
  DEALLOCATE (gTopoInp, gTopoOut, gTopoDel)
  DEALLOCATE (gLnPsInp, gLnPsOut, gPsfcInp, gPsfcOut)
  DEALLOCATE (gUvelInp, gVvelInp)
  DEALLOCATE (gDivgInp, gVortInp)
  DEALLOCATE (gTvirInp)
  DEALLOCATE (gPresInp)
  DEALLOCATE (gUvelOut, gVvelOut)
  DEALLOCATE (gTvirOut)
  DEALLOCATE (gPresOut)

END SUBROUTINE ClsArrays


SUBROUTINE GetSpHuTracers

  IMPLICIT NONE

  ALLOCATE (qSpHuInp (Mnwv2Inp,KmaxInp,NTracers))
  ALLOCATE (qSpHuOut (Mnwv2Out,KmaxOut,NTracers))
  ALLOCATE (gSpHuInp (ImaxOut,JmaxOut,KmaxInp,NTracers))
  ALLOCATE (gSpHuOut (ImaxOut,JmaxOut,KmaxOut,NTracers))

END SUBROUTINE GetSpHuTracers


SUBROUTINE ClsSpHuTracers

  IMPLICIT NONE

  DEALLOCATE (qSpHuInp, qSpHuOut)
  DEALLOCATE (gSpHuInp, gSpHuOut)

END SUBROUTINE ClsSpHuTracers


END MODULE InputArrays
