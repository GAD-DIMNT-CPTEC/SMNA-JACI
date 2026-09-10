!
!  $Author: pkubota $
!  $Date: 2007/09/18 13:32:24 $
!  $Revision: 1.3 $
!
MODULE PrblSize

  USE Constants, ONLY : r8, nfinp, nfprt,plevs
  USE Parallelism, ONLY : myid

  IMPLICIT NONE

  PRIVATE

  INTEGER, PUBLIC :: nGaus
  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: Pmand
  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: aLnPmd

  PUBLIC :: CreatePrblSize
  INTEGER, PUBLIC :: trunc_in, nlon, nlat, vert_in, vert_out

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: PmandInput


CONTAINS


  SUBROUTINE CreatePrblSize ()

    USE Constants, ONLY : Trunc, Lev, nMand

    IMPLICIT NONE

    READ (Trunc(3:6), FMT='(I4)') trunc_in

    READ (Lev(2:4), FMT='(I3)') vert_in

    CALL GetImaxJmax (trunc_in, nlon, nlat)

    nGaus= nlon*nlat
    CALL CreatePressureLevels ()

  END SUBROUTINE CreatePrblSize


  SUBROUTINE CreatePressureLevels ()

    USE Constants, ONLY : nMand

    IMPLICIT NONE

    INTEGER :: k, lmax

    CHARACTER (LEN=9) :: FmtMand='(   F8.2)'

    IF (nMand <= 0) THEN
      nMand=1
      ALLOCATE (PmandInput(nMand))
      PmandInput(nMand)=0.0_r8
    ELSE
      ALLOCATE (PmandInput(nMand))
      IF (myid.eq.0) THEN
         IF (nMand < 10) THEN
           WRITE (FmtMand(2:4), FMT='(I3.3)') nMand
         ELSE
           WRITE (FmtMand(2:4), FMT='(I3.3)') 10
         END IF
      END IF
      DO k=1,nMand
         PmandInput(k)=plevs(k)
      END DO
    END IF

    IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,/)') ' PmandInput:'
    IF (myid.eq.0) WRITE (UNIT=nfprt, FMT=FmtMand) PmandInput
    IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(A)') ' '

    IF (PmandInput(nMand) > 0.0_r8) THEN
      Lmax=nMand
    ELSE
      Lmax=18
    END IF
    IF (myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,2(A,I5),A,F8.2,/)') &
                     ' Lmax = ', Lmax, '   nMand = ', nMand, &
                     '   PmandInput(nMand) = ', PmandInput(nMand)

    ALLOCATE (Pmand(Lmax), aLnPmd(Lmax))

    IF (PmandInput(nMand) == 0.0_r8) THEN
       Pmand = (/ 1000.0_r8, 925.0_r8, 850.0_r8, 775.0_r8, 700.0_r8, 500.0_r8, &
                   400.0_r8, 300.0_r8, 250.0_r8, 200.0_r8, 150.0_r8, 100.0_r8, &
                    70.0_r8,  50.0_r8,  30.0_r8,  20.0_r8,  10.0_r8,   3.0_r8 /)
    ELSE
       Pmand=PmandInput
    END IF

    DO k=1,Lmax
      aLnPmd(k)=LOG(Pmand(k))
    ENDDO
    vert_out = lmax

  END SUBROUTINE CreatePressureLevels


  SUBROUTINE GetImaxJmax (Mend, Imax, Jmax)

    USE Constants, ONLY : Linear

    IMPLICIT NONE

    INTEGER, INTENT (IN) :: Mend
    INTEGER, INTENT (OUT) :: Imax, Jmax

    INTEGER :: Nx, Nm, N2m, N3m, N5m, &
                         n2, n3, n5, j, n, Check, Jfft

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

    IF (Linear) THEN
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
    IF (MOD(Jmax,2) /= 0) Jmax=Jmax+1

    DEALLOCATE (Ifft)

  END SUBROUTINE GetImaxJmax


END MODULE PrblSize
