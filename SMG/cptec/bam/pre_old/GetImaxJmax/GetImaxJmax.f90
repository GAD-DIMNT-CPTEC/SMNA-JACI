PROGRAM GetIJ

   IMPLICIT NONE

   INTEGER, PARAMETER :: Kr=8

   INTEGER :: Mend, Imax, Jmax

   REAL (KIND=Kr) :: Dl, Dx, Dkm=112.0_Kr

   LOGICAL :: LinearGrid

   CHARACTER (LEN=1) :: CLG

   CHARACTER (LEN=5) :: CMend

   CALL GETARG (1,CMend)
   CALL GETARG (2,CLG)

   READ (CMend, '(I5)') Mend
   IF (CLG == 'L') THEN
      LinearGrid=.TRUE.
      WRITE (UNIT=*, FMT='(/,A)') &
            ' Linear Triangular Truncation : '
   ELSE
      LinearGrid=.FALSE.
      WRITE (UNIT=*, FMT='(/,A)') &
            ' Quadratic Triangular Truncation : '
   END IF

   CALL GetImaxJmax ()

   Dl=360.0_Kr/REAL(Imax,Kr)
   Dx=Dl*Dkm

   WRITE (UNIT=*, FMT='(/,3(A,I5,/))') &
         ' Mend : ', Mend, ' Imax : ', Imax, ' Jmax : ', Jmax
   WRITE (UNIT=*, FMT='(A,F13.9,A)')   ' Dl: ', Dl, ' Degrees'
   WRITE (UNIT=*, FMT='(A,F13.2,A,/)') ' Dx: ', Dx, ' km'


CONTAINS


SUBROUTINE GetImaxJmax ()

  IMPLICIT NONE

  INTEGER, PARAMETER :: r8=8

  INTEGER :: Nx, Nm, N2m, N3m, N5m, n2, n3, n5, j, n, Check, Jfft

  INTEGER, SAVE :: Lfft=40000

  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: Ifft

  N2m=CEILING(LOG(REAL(Lfft,r8))/LOG(2.0_r8))
  N3m=CEILING(LOG(REAL(Lfft,r8))/LOG(3.0_r8))
  N5m=CEILING(LOG(REAL(Lfft,r8))/LOG(5.0_r8))
  Nx=N2m*(N3m+1)*(N5m+1)

  ALLOCATE (Ifft (Nx))
  Ifft=0

  n=0
  DO n2=1,N2m
     Jfft=(2**n2)
     IF (Jfft > Lfft) EXIT
     DO n3=0,N3m
        Jfft=(2**n2)*(3**n3)
        IF (Jfft > Lfft) EXIT
        DO n5=0,N5m
           Jfft=(2**n2)*(3**n3)*(5**n5)
           IF (Jfft > Lfft) EXIT
           n=n+1
           Ifft(n)=Jfft
        END DO
     END DO
  END DO
  Nm=n

  n=0
  DO 
     Check=0
     n=n+1
     DO j=1,Nm-1
        IF (Ifft(j) > Ifft(j+1)) THEN
           Jfft=Ifft(j)
           Ifft(j)=Ifft(j+1)
           Ifft(j+1)=Jfft
           Check=1
        END IF
     END DO
     IF (Check == 0) EXIT
  END DO

  IF (LinearGrid) THEN
     Jfft=2
  ELSE
     Jfft=3
  END IF
  Imax=Jfft*Mend+1
  DO n=1,Nm
     IF (Ifft(n) >= Imax) THEN
        Imax=Ifft(n)
        EXIT
     END IF
  END DO
  Jmax=Imax/2
  IF (MOD(Jmax, 2) /= 0) Jmax=Jmax+1

  DEALLOCATE (Ifft)

END SUBROUTINE GetImaxJmax


END PROGRAM GetIJ
