MODULE Spectral2Grid

  USE InputParameters, ONLY : r8, Imax, Imx, Jmax, Mend1, Mend2, Mnwv2, JmaxHf

  USE FastFourierTransform, ONLY : InvFFT

  USE LegendreTransform, ONLY : la0, Spec2Four

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Transpose, SpecCoef2Grid, SpecCoef2GridD


CONTAINS


  SUBROUTINE Transpose (qf)

    IMPLICIT NONE

    ! qf(Mnwv2) input: spectral representation of a
    !                  global field coluMnwise storage
    !          output: spectral representation of a
    !                  global field diagonalwise storage

    REAL (KIND=r8), INTENT(INOUT) :: qf(Mnwv2)

    REAL (KIND=r8) :: qw(Mnwv2)

    INTEGER :: k
    INTEGER :: l
    INTEGER :: lx
    INTEGER :: mn
    INTEGER :: mm
    INTEGER :: nlast
    INTEGER :: nn

    l=0
    DO mm=1,Mend1
       nlast=Mend2-mm
       DO nn=1,nlast
          l=l+1
          lx=la0(mm,nn)
          qw(2*lx-1)=qf(2*l-1)
          qw(2*lx)=qf(2*l)
       END DO
    END DO
    DO mn=1,Mnwv2
       qf(mn)=qw(mn)
    END DO

  END SUBROUTINE Transpose


  SUBROUTINE SpecCoef2Grid (qf, gf)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: qf(Mnwv2)

    REAL (KIND=r8), INTENT(OUT) :: gf(Imax,Jmax)

    REAL (KIND=r8) :: Four(Imx,Jmax)

    CALL Spec2Four (qf, Four, .FALSE.)
    CALL InvFFT (Four, gf)

  END SUBROUTINE SpecCoef2Grid


  SUBROUTINE SpecCoef2GridD (qf, gf)

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: qf(Mnwv2)

    REAL (KIND=r8), INTENT(OUT) :: gf(Imax,Jmax)

    REAL (KIND=r8) :: Four(Imx,Jmax)

    CALL Spec2Four (qf, Four, .TRUE.)
    CALL InvFFT (Four, gf)
    gf(:,JmaxHf+1:Jmax)=-gf(:,JmaxHf+1:Jmax)

  END SUBROUTINE SpecCoef2GridD


END MODULE Spectral2Grid
