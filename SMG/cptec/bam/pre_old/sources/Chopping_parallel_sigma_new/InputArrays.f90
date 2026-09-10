!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE InputArrays

  USE InputParameters, ONLY: i4, r4, r8, KmaxInp, KmaxInpp, NTracers, &
                             ImaxOut, JmaxOut, KmaxOut, KmaxOutp,ImaxInp,jmaxInp
  USE Sizes,           ONLY: ibMax, jbMax,kMax, kMaxloc_out, kmaxloc_in, &
                             iMax,jMax, jMaxHalf, jMinPerM, iMaxPerJ, &
                             mMax, nExtMax, mnMax, mnExtMax, mnExtMap, &
                             mymnMax, mymnExtMax, mnMax_out,mnMax_out

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GetArrays, ClsArrays, GetSpHuTracers, ClsSpHuTracers

  INTEGER (KIND=i4), DIMENSION (4), PUBLIC :: DateInitial, DateCurrent

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  DelSInp, SigIInp, SigLInp, SigIOut, SigLOut

  REAL (KIND=r4), DIMENSION (:), ALLOCATABLE, PUBLIC :: qWorkInp, qWorkprOut, qtorto
  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: qWorkOut, qWorkOut1,qWorkInOut, qWorkInOut1

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  DelSigmaInp, SigInterInp, SigLayerInp, &
                  DelSigmaOut, SigInterOut, SigLayerOut

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
                  qTopoInp, qLnPsInp, qTopoOut,qTopoOutSpec, qLnPsOut

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  qDivgInp, qVortInp, qTvirInp, gWorkprInp,&
                  qDivgOut, qVortOut, qTvirOut,&
                  qUvelInp, qVvelInp, qUvelOut,qVvelOut

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: &
                  qSpHuInp, qSpHuOut

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: gWorkOut
  REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: gWorkprout

  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE, PUBLIC :: &
                  gTopoInp, gTopoOut,gTopoOutGaus, gTopoOutGaus8,gTopoDel, &
                  gLnPsInp, gPsfcInp, gLnPsOut, gPsfcOut, gpresaux

  REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE, PUBLIC :: &
                  gUvelInp, gVvelInp, gTvirInp, &
                  gDivgInp, gVortInp, gPresInp, gPresInpp, &
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
!  ALLOCATE (qTopoRec    (Mnwv2Rec))
  ALLOCATE (qWorkInp    (2*mnmax))
!  ALLOCATE (qWorkInpAux (2*mnMax_out))
  ALLOCATE (qWorkOut    (2*mnmax_Out,KmaxOut))
  ALLOCATE (qWorkOut1   (2*mnmax_Out,KmaxOut))
  ALLOCATE (qWorkInOut  (2*mnmax_Out,KmaxInp))
  ALLOCATE (qWorkInOut1 (2*mnmax_Out,KmaxInp))
  ALLOCATE (qWorkprOut  (2*mnmax_Out))
  ALLOCATE (qtorto      (2*mnmax_Out))
  ALLOCATE (SigIInp     (KmaxInpp), SigLInp (KmaxInpp), DelSInp (KmaxInp))
  ALLOCATE (SigIOut     (KmaxOutp), SigLOut (KmaxOutp))
  ALLOCATE (DelSigmaInp (KmaxInpp), SigInterInp (KmaxInpp), SigLayerInp (KmaxInpp))
  ALLOCATE (DelSigmaOut (KmaxOut) , SigInterOut (KmaxOutp), SigLayerOut (KmaxOutp))
  ALLOCATE (qTopoInp    (Mnwv2Inp), qLnPsInp (Mnwv2Inp))
  ALLOCATE (qTopoOut    (Mnwv2Out),qTopoOutSpec(Mnwv2Out), qLnPsOut (Mnwv2Out))
  ALLOCATE (qDivgInp (Mnwv2Inp,Kmaxloc_In), &
            qVortInp (Mnwv2Inp,Kmaxloc_In))
  ALLOCATE (qDivgOut (Mnwv2Out,Kmaxloc_Out),&
            qVortOut (Mnwv2Out,Kmaxloc_Out))
  ALLOCATE (qUvelInp (Mnwv3Out,Kmaxloc_In), &
            qVvelInp (Mnwv3Out,Kmaxloc_In))
  ALLOCATE (qUvelOut (Mnwv3Out,Kmaxloc_Out), &
            qVvelOut (Mnwv3Out,Kmaxloc_Out))
  ALLOCATE (qTvirInp (Mnwv2Inp,Kmaxloc_In))
  ALLOCATE (qTvirOut (Mnwv2Out,Kmaxloc_Out))
  ALLOCATE (gWorkOut (ImaxOut,JmaxOut,max(KmaxOut,KmaxInp)))
  ALLOCATE (gWorkprOut (ImaxOut,JmaxOut), &
            gWorkprInp (ImaxInp,jmaxInp))
  ALLOCATE (gTopoInp (Ibmax,jbmax),     &
            gTopoOut (Ibmax,jbmax),     &
            gTopoOutGaus(Ibmax,jbmax),     &
            gTopoOutGaus8(ImaxOut,JmaxOut),     &
            gTopoDel (Ibmax,jbmax),     &
            gLnPsInp (Ibmax,jbmax),     &
            gLnPsOut (Ibmax,jbmax),     &
            gPsfcInp (Ibmax,jbmax),     &
            gPsfcOut (Ibmax,jbmax))    
  ALLOCATE (gpresaux (Ibmax, KmaxOutp))
  ALLOCATE (gUvelInp (Ibmax, KmaxInp, Jbmax) , &
            gVvelInp (Ibmax, KmaxInp, Jbmax))
  ALLOCATE (gDivgInp (Ibmax, KmaxInp, Jbmax) , &
            gVortInp (Ibmax, KmaxInp, Jbmax))
  ALLOCATE (gTvirInp (Ibmax, KmaxInp, Jbmax) , &
            gPresInp (Ibmax, KmaxInp, Jbmax) , &
            gPresInpp(Ibmax, KmaxInpp,Jbmax))
  ALLOCATE (gUvelOut (Ibmax, KmaxOut, Jbmax) , &
            gVvelOut (Ibmax, KmaxOut, Jbmax))
  ALLOCATE (gTvirOut (Ibmax, KmaxOut, Jbmax) , &
            gPresOut (Ibmax, KmaxOut, Jbmax))

 gPsfcOut=0.0_r8

!  DateInitial=0;DateCurrent=0

!  DelSInp=0.0_r4; SigIInp=0.0_r4;SigLInp=0.0_r4;SigIOut=0.0_r4; SigLOut=0.0_r4

!  qWorkInp=0.0_r4; qWorkprOut=0.0_r4; qtorto=0.0_r4;
!  qWorkOut=0.0_r8; qWorkOut1=0.0_r8;qWorkInOut=0.0_r8; qWorkInOut1=0.0_r8

!  DelSigmaInp=0.0_r8; SigInterInp=0.0_r8; SigLayerInp=0.0_r8
!  DelSigmaOut=0.0_r8; SigInterOut=0.0_r8; SigLayerOut=0.0_r8

!  qTopoInp=0.0_r8; qLnPsInp=0.0_r8; qTopoOut=0.0_r8;qTopoOutSpec=0.0_r8
!  qLnPsOut=0.0_r8

!  qDivgInp=0.0_r8; qVortInp=0.0_r8; qTvirInp=0.0_r8; gWorkprInp=0.0_r8
!  qDivgOut=0.0_r8; qVortOut=0.0_r8
!  qTvirOut=0.0_r8
!  qUvelInp=0.0_r8; qVvelInp=0.0_r8; qUvelOut=0.0_r8
!  qVvelOut=0.0_r8


!  gWorkOut=0.0_r8
!  gWorkprout=0.0_r4

!  gTopoInp=0.0_r8; gTopoOut=0.0_r8;gTopoOutGaus=0.0_r8
!  gTopoOutGaus8=0.0_r8;gTopoDel=0.0_r8
!  gLnPsInp=0.0_r8; gPsfcInp=0.0_r8; gLnPsOut=0.0_r8; gPsfcOut=0.0_r8; gpresaux=0.0_r8

!  gUvelInp=0.0_r8; gVvelInp=0.0_r8; gTvirInp=0.0_r8
!  gDivgInp=0.0_r8; gVortInp=0.0_r8; gPresInp=0.0_r8; gPresInpp=0.0_r8
!  gUvelOut=0.0_r8; gVvelOut=0.0_r8; gTvirOut=0.0_r8
!  gPresOut=0.0_r8;



END SUBROUTINE GetArrays


SUBROUTINE ClsArrays

  IMPLICIT NONE

!  DEALLOCATE (qTopoRec )
  DEALLOCATE (qWorkInp )
!  DEALLOCATE (qWorkInpAux )
  DEALLOCATE (qWorkOut    )
  DEALLOCATE (qWorkOut1   )
  DEALLOCATE (qWorkInOut  )
  DEALLOCATE (qWorkInOut1 )
  DEALLOCATE (qWorkprOut  )
  DEALLOCATE (qtorto      )
  DEALLOCATE (SigIInp , SigLInp , DelSInp)
  DEALLOCATE (SigIOut , SigLOut )
  DEALLOCATE (DelSigmaInp , SigInterInp , SigLayerInp )
  DEALLOCATE (DelSigmaOut , SigInterOut , SigLayerOut )
  DEALLOCATE (qTopoInp , qLnPsInp )
  DEALLOCATE (qTopoOut ,qTopoOutSpec, qLnPsOut )
  DEALLOCATE (qDivgInp , &
              qVortInp )
  DEALLOCATE (qDivgOut , &
              qVortOut  )
  DEALLOCATE (qUvelInp , &
              qVvelInp )
  DEALLOCATE (qUvelOut , &
              qVvelOut )
  DEALLOCATE (qTvirInp )
  DEALLOCATE (qTvirOut )
  DEALLOCATE (gWorkOut )
  DEALLOCATE (gWorkprOut , &
              gWorkprInp )
  DEALLOCATE (gTopoInp ,     &
            gTopoOut ,     &
            gTopoOutGaus,     &
            gTopoOutGaus8,     &
            gTopoDel ,     &
            gLnPsInp ,     &
            gLnPsOut ,     &
            gPsfcInp ,     &
            gPsfcOut )    
  DEALLOCATE (gpresaux )
  DEALLOCATE (gUvelInp  , &
              gVvelInp )
  DEALLOCATE (gDivgInp  , &
              gVortInp )
  DEALLOCATE (gTvirInp  , &
              gPresInp  , &
              gPresInpp)
  DEALLOCATE (gUvelOut  , &
              gVvelOut )
  DEALLOCATE (gTvirOut  , &
              gPresOut )

END SUBROUTINE ClsArrays


SUBROUTINE GetSpHuTracers

  IMPLICIT NONE
  IF (.NOT.ALLOCATED(qSpHuInp))THEN
      ALLOCATE (qSpHuInp (Mnwv2Inp,KmaxInp,NTracers+2))
  ELSE
      DEALLOCATE (qSpHuInp)
      ALLOCATE (qSpHuInp (Mnwv2Inp,KmaxInp,NTracers+2))
  END IF

  IF (.NOT.ALLOCATED(qSpHuOut))THEN
      ALLOCATE (qSpHuOut (Mnwv2Out,KmaxOut,NTracers+2))
  ELSE
      DEALLOCATE (qSpHuOut)
      ALLOCATE (qSpHuOut (Mnwv2Out,KmaxOut,NTracers+2))
  END IF

  IF (.NOT.ALLOCATED(gSpHuInp))THEN
      ALLOCATE (gSpHuInp (Ibmax,KmaxInp,Jbmax,NTracers+2))
  ELSE
      DEALLOCATE (gSpHuInp)
      ALLOCATE (gSpHuInp (Ibmax,KmaxInp,Jbmax,NTracers+2))
  END IF

  IF (.NOT.ALLOCATED(gSpHuOut))THEN
      ALLOCATE (gSpHuOut (Ibmax,KmaxOut,Jbmax,NTracers+2))
  ELSE
      DEALLOCATE (gSpHuOut)
      ALLOCATE (gSpHuOut (Ibmax,KmaxOut,Jbmax,NTracers+2))
  END IF
  
END SUBROUTINE GetSpHuTracers


SUBROUTINE ClsSpHuTracers

  IMPLICIT NONE
  IF (ALLOCATED(qSpHuInp))THEN
      DEALLOCATE (qSpHuInp)
  END IF

  IF (ALLOCATED(qSpHuOut))THEN
      DEALLOCATE (qSpHuOut)
  END IF

  IF (ALLOCATED(gSpHuInp))THEN
      DEALLOCATE (gSpHuInp)
  END IF

  IF (ALLOCATED(gSpHuOut))THEN
      DEALLOCATE (gSpHuOut)
  END IF
  

END SUBROUTINE ClsSpHuTracers


END MODULE InputArrays
