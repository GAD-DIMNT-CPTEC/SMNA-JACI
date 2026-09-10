!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE SpectralGrid

  USE InputParameters, ONLY : r8, Imax, Imx, Jmax, Mnwv2

  USE FastFourierTransform, ONLY : InvFFT, DirFFT

  USE LegendreTransform, ONLY : transs, Spec2Four, Four2Spec

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Grid2SpecCoef, SpecCoef2Grid


CONTAINS


  SUBROUTINE SpecCoef2Grid (Ldim, qf, f)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Ldim

    REAL (KIND=r8), INTENT(IN OUT) :: qf(Mnwv2,Ldim)

    REAL (KIND=r8), INTENT(OUT) :: f(Imax,Jmax,Ldim)

    REAL (KIND=r8) :: Four(Imx,Jmax,Ldim)

    INTEGER :: isign=-1

    CALL transs (Ldim, isign, qf)
    CALL Spec2Four (qf, Four)
    CALL InvFFT (Four, f)

  END SUBROUTINE SpecCoef2Grid


  SUBROUTINE Grid2SpecCoef (Ldim, f, qf)

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  Ldim

    REAL (KIND=r8), INTENT(IN) :: f(Imax,Jmax,Ldim)

    REAL (KIND=r8), INTENT(OUT) :: qf(Mnwv2,Ldim)

    REAL (KIND=r8) :: Four(Imx,Jmax,Ldim)

    CALL DirFFT (f, Four)
    CALL Four2Spec (Four, qf)

  END SUBROUTINE Grid2SpecCoef


END MODULE SpectralGrid
