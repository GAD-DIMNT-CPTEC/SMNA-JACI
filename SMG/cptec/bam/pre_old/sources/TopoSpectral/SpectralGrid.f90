!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE SpectralGrid

  USE InputParameters, ONLY : r8, Imax, Imx, Jmax, Mnwv2

  USE FastFourierTransform, ONLY : InvFFT, DirFFT

  USE LegendreTransform, ONLY : Spec2Four, Four2Spec

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Grid2SpecCoef, SpecCoef2Grid


CONTAINS


  SUBROUTINE SpecCoef2Grid (Ldim, qf, gf)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Ldim

    REAL (KIND=r8), INTENT(IN OUT) :: qf(Mnwv2,Ldim)

    REAL (KIND=r8), INTENT(OUT) :: gf(Imax,Jmax,Ldim)

    REAL (KIND=r8) :: Four(Imx,Jmax,Ldim)

    CALL Spec2Four (qf, Four)
    CALL InvFFT (Four, gf)

  END SUBROUTINE SpecCoef2Grid


  SUBROUTINE Grid2SpecCoef (Ldim, gf, qf)

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  Ldim

    REAL (KIND=r8), INTENT(IN) :: gf(Imax,Jmax,Ldim)

    REAL (KIND=r8), INTENT(OUT) :: qf(Mnwv2,Ldim)

    REAL (KIND=r8) :: Four(Imx,Jmax,Ldim)

    CALL DirFFT (gf, Four)
    CALL Four2Spec (Four, qf)

  END SUBROUTINE Grid2SpecCoef


END MODULE SpectralGrid
