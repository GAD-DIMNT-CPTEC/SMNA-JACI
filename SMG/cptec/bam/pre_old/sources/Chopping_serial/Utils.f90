!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE Utils

  USE InputParameters, ONLY : r8, &
                              ImaxOut, JmaxOut, KmaxInp, KmaxInpp, KmaxOut, KmaxOutp, &
                              Rd, Cp, GEps, Gamma, Grav, Rv

  USE InputArrays, ONLY : DelSigmaOut, SigLayerOut, SigInterOut, &
                          DelSInp, SigLInp, SigIInp, &
                          gTopoInp, gPsfcInp, gTvirInp, gPresInp, &
                          gTopoOut, gPsfcOut, gLnPsOut

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: NewSigma, SigmaInp, NewPs


CONTAINS


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

  REAL (KIND=r8), DIMENSION (ImaxOut,JmaxOut) :: Zu

  goRd=Grav/Rd

  ! Compute Surface Pressure Below the Original Ground
  ls=0
  k=1
  GammaLR=Gamma
  DO j=1,JmaxOut
     DO i=1,ImaxOut
        Zu(i,j)=gTopoInp(i,j)-gTvirInp(i,j,k)/GammaLR* &
                ((gPsfcInp(i,j)/gPresInp(i,j,k))**(-GammaLR/goRd)-1.0_r8)
        IF (gTopoOut(i,j) <= Zu(i,j)) THEN
           IF (ABS(GammaLR) > GEps) THEN
              gPsfcOut(i,j)=gPresInp(i,j,k)*(1.0_r8+GammaLR/gTvirInp(i,j,k)* &
                            (gTopoOut(i,j)-Zu(i,j)))**(-goRd/GammaLR)
           ELSE
              gPsfcOut(i,j)=gPresInp(i,j,k)*EXP(-goRd/gTvirInp(i,j,k)* &
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
        DO j=1,JmaxOut
           DO i=1,ImaxOut
              IF (gPsfcOut(i,j) == 0.0_r8) THEN
                 GammaLR=-goRd*LOG(gTvirInp(i,j,k-1)/gTvirInp(i,j,k))/ &
                         LOG(gPresInp(i,j,k-1)/gPresInp(i,j,k))
                 IF (ABS(GammaLR) > GEps) THEN
                    Zu(i,j)=Zu(i,j)-gTvirInp(i,j,k)/GammaLR* &
                            ((gPresInp(i,j,k-1)/gPresInp(i,j,k))** &
                            (-GammaLR/goRd)-1.0_r8)
                 ELSE
                    Zu(i,j)=Zu(i,j)+gTvirInp(i,j,k)/ &
                            goRd*LOG(gPresInp(i,j,k-1)/gPresInp(i,j,k))
                 END IF
                 IF (gTopoOut(i,j) <= Zu(i,j)) THEN
                    IF (ABS(GammaLR) > GEps) THEN
                       gPsfcOut(i,j)=gPresInp(i,j,k)*(1.0_r8+GammaLR/gTvirInp(i,j,k)* &
                                     (gTopoOut(i,j)-Zu(i,j)))**(-goRd/GammaLR)
                    ELSE
                       gPsfcOut(i,j)=gPresInp(i,j,k)* &
                                     EXP(-goRd/gTvirInp(i,j,k)*(gTopoOut(i,j)-Zu(i,j)))
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
     DO j=1,JmaxOut
        DO i=1,ImaxOut
           IF (gPsfcOut(i,j) == 0.0_r8) THEN
              gPsfcOut(i,j)=gPresInp(i,j,k)* &
                            EXP(-goRd/gTvirInp(i,j,k)*(gTopoOut(i,j)-Zu(i,j)))
           END IF
        END DO
     END DO
  END IF

  ! ln(ps) in cBar (from mBar)
  gLnPsOut=LOG(gPsfcOut/10.0_r8)

END SUBROUTINE NewPs


END MODULE Utils
