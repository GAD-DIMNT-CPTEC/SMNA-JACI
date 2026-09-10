!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE InputArrays

  USE InputParameters, ONLY: i4, r4, r8, KmaxInp, KmaxInpp, NTracers, &
                             ImaxOut, JmaxOut, KmaxOut, KmaxOutp
  USE Sizes,           ONLY: ibMax, jbMax,kMax, kMaxloc_out, kmaxloc_in, &
                             jMax, jMaxHalf, jMinPerM, iMaxPerJ, &
                             mMax, nExtMax, mnMax, mnExtMax, mnExtMap, &
                             mymnMax, mymnExtMax, mnMax_out

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GetArrays, ClsArrays, GetSpHuTracers, ClsSpHuTracers

  INTEGER (KIND=i4), DIMENSION (4), PUBLIC :: DateInitial, DateCurrent

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  DelSInp, SigIInp, SigLInp, SigIOut, SigLOut

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: qWorkInp, qWorkprOut, qtorto
  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: qWorkOut, qWorkOut1

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  DelSigmaInp, SigInterInp, SigLayerInp, &
                  DelSigmaOut, SigInterOut, SigLayerOut

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  qTopoInp, qLnPsInp, qTopoOut, qLnPsOut, qTopoRec

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  qDivgInp, qVortInp, qTvirInp, &
                  qDivgOut, qVortOut, qTvirOut, &
                  qUvelInp, qVvelInp, qUvelOut, qVvelOut

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: &
                  qSpHuInp, qSpHuOut

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: gWorkOut
  REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: gWorkprout

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

  INTEGER :: Mnwv2Rec, Mnwv2Out, Mnwv2Inp, Mnwv3Inp, Mnwv3Out

CONTAINS


SUBROUTINE GetArrays

  IMPLICIT NONE

  Mnwv2Inp = 2*mymnMax
  Mnwv2Out = Mnwv2Inp
  Mnwv3Inp = 2*mymnExtMax
  Mnwv3Out = Mnwv3Inp
  Mnwv2Rec=Mnwv2Inp
  ALLOCATE (qTopoRec (Mnwv2Rec))
  ALLOCATE (qWorkInp (2*mnmax))
  ALLOCATE (qWorkOut (2*mnmax_Out,KmaxOut))
  ALLOCATE (qWorkOut1(2*mnmax_Out,KmaxOut))
  ALLOCATE (qWorkprOut (2*mnmax_Out))
  ALLOCATE (qtorto     (2*mnmax_Out))
  ALLOCATE (SigIInp (KmaxInpp), SigLInp (KmaxInp), DelSInp (KmaxInp))
  ALLOCATE (SigIOut (KmaxOutp), SigLOut (KmaxOut))
  ALLOCATE (DelSigmaInp (KmaxInp), SigInterInp (KmaxInpp), SigLayerInp (KmaxInp))
  ALLOCATE (DelSigmaOut (KmaxOut), SigInterOut (KmaxOutp), SigLayerOut (KmaxOut))
  ALLOCATE (qTopoInp (Mnwv2Inp), qLnPsInp (Mnwv2Inp))
  ALLOCATE (qTopoOut (Mnwv2Out), qLnPsOut (Mnwv2Out))
  ALLOCATE (qDivgInp (Mnwv2Inp,Kmaxloc_In), &
            qVortInp (Mnwv2Inp,Kmaxloc_In))
  ALLOCATE (qDivgOut (Mnwv2Out,Kmaxloc_Out), &
            qVortOut (Mnwv2Out,Kmaxloc_Out))
  ALLOCATE (qUvelInp (Mnwv3Out,Kmaxloc_In), &
            qVvelInp (Mnwv3Out,Kmaxloc_In))
  ALLOCATE (qUvelOut (Mnwv3Out,Kmaxloc_Out), &
            qVvelOut (Mnwv3Out,Kmaxloc_Out))
  ALLOCATE (qTvirInp (Mnwv2Inp,Kmaxloc_In))
  ALLOCATE (qTvirOut (Mnwv2Out,Kmaxloc_Out))
  ALLOCATE (gWorkOut (ImaxOut,JmaxOut,max(KmaxOut,KmaxInp)))
  ALLOCATE (gWorkprOut (ImaxOut,JmaxOut))
  ALLOCATE (gTopoInp (Ibmax,jbmax),     &
            gTopoOut (Ibmax,jbmax),     &
            gTopoDel (Ibmax,jbmax),     &
            gLnPsInp (Ibmax,jbmax),     &
            gLnPsOut (Ibmax,jbmax),     &
            gPsfcInp (Ibmax,jbmax),     &
            gPsfcOut (Ibmax,jbmax))    
  ALLOCATE (gUvelInp (Ibmax, KmaxInp, Jbmax) , &
            gVvelInp (Ibmax, KmaxInp, Jbmax))
  ALLOCATE (gDivgInp (Ibmax, KmaxInp, Jbmax) , &
            gVortInp (Ibmax, KmaxInp, Jbmax))
  ALLOCATE (gTvirInp (Ibmax, KmaxInp, Jbmax) , &
            gPresInp (Ibmax, KmaxInp, Jbmax))
  ALLOCATE (gUvelOut (Ibmax, KmaxOut, Jbmax) , &
            gVvelOut (Ibmax, KmaxOut, Jbmax))
  ALLOCATE (gTvirOut (Ibmax, KmaxOut, Jbmax) , &
            gPresOut (Ibmax, KmaxOut, Jbmax))

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
  DEALLOCATE (qUvelOut, qVvelOut)
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

  ALLOCATE (qSpHuInp (Mnwv2Inp,KmaxInp,NTracers+2))
  ALLOCATE (qSpHuOut (Mnwv2Out,KmaxOut,NTracers+2))
  ALLOCATE (gSpHuInp (Ibmax,KmaxInp,Jbmax,NTracers+2))
  ALLOCATE (gSpHuOut (Ibmax,KmaxOut,Jbmax,NTracers+2))

END SUBROUTINE GetSpHuTracers


SUBROUTINE ClsSpHuTracers

  IMPLICIT NONE

  DEALLOCATE (qSpHuInp, qSpHuOut)
  DEALLOCATE (gSpHuInp, gSpHuOut)

END SUBROUTINE ClsSpHuTracers


END MODULE InputArrays
