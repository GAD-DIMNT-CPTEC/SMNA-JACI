MODULE MiscMod

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GetLongitudes
   PUBLIC :: GetImaxJmax
   PUBLIC :: GetGaussianLatitudes

   !
   !precisao dos dados do modelo BAM
   integer, public, parameter :: I4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
   integer, public, parameter :: I8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
   integer, public, parameter :: R4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
   integer, public, parameter :: R8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
   integer, public, parameter :: strlen = 1024

   !Logical Units 
   integer, parameter :: stderr = 0 ! Error Unit
   integer, parameter :: stdinp = 5 ! Input Unit
   integer, parameter :: stdout = 6 ! Output Unit


   ! Constants
   real(kind=r8), public, parameter :: rd     = 45_r8/ATAN(1.0_r8) ! convert to radian
   real(kind=r8), public, parameter :: EMRad   = 6.37E6_r8         ! Earth Mean Radius (m)
   real(kind=r8), public, parameter :: EMRad1  = 1.0_r8/EMRad      ! 1/EMRad (1/m)
   real(kind=r8), public, parameter :: EMRad12 = EMRad1*EMRad1     ! EMRad1**2 (1/m2)
   real(kind=r8), public, parameter :: EMRad2  = EMRad*EMRad       ! EMRad**2 (m2)

   integer :: iMax
   integer :: jMax
!   integer :: kMax
!   integer :: Mend
   integer :: Mend1
   integer :: Mend2
   integer :: Mend3
   integer :: MnWv0
   integer :: MnWv1
   integer :: MnWv2
   integer :: MnWv3
   integer :: jMaxHf
   integer :: iMx


CONTAINS


SUBROUTINE InitParameters ( mend )

  IMPLICIT NONE

  integer, intent(in) :: mend

  
  CALL GetImaxJmax (Mend, Imax, Jmax)

  Mend1 = Mend+1
  Mend2 = Mend+2
  Mend3 = Mend+3
  Mnwv2 = Mend1*Mend2
  Mnwv0 = Mnwv2/2
  Mnwv3 = Mnwv2+2*Mend1
  Mnwv1 = Mnwv3/2

  Imx    = Imax+2
  Jmaxhf = Jmax/2

END SUBROUTINE InitParameters


SUBROUTINE GetImaxJmax (Mend, Imax, Jmax)

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

  Imax=3*Mend+1
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

SUBROUTINE GetGaussianLatitudes (Jmax, Lat)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: Jmax ! Number of Gaussian Latitudes

   REAL (KIND=r8), DIMENSION (Jmax), INTENT(OUT) :: &
                   Lat ! Gaussian Latitudes In Degree

   INTEGER :: j

   REAL (KIND=r8) :: eps, rd, dCoLatRadz, CoLatRad, dCoLatRad, p2, p1

   eps=1.0e-12_r8
   rd=45.0_r8/ATAN(1.0_r8)
   dCoLatRadz=((180.0_r8/REAL(Jmax,r8))/rd)/10.0_r8
   CoLatRad=0.0_r8
   DO j=1,Jmax/2
      dCoLatRad=dCoLatRadz
      DO WHILE (dCoLatRad > eps)
         CALL LegendrePolynomial (Jmax,CoLatRad,p2)
         DO
            p1=p2
            CoLatRad=CoLatRad+dCoLatRad
            CALL LegendrePolynomial (Jmax,CoLatRad,p2)
            IF (SIGN(1.0_r8,p1) /= SIGN(1.0_r8,p2)) EXIT
         END DO
         CoLatRad=CoLatRad-dCoLatRad
         dCoLatRad=dCoLatRad*0.25_r8
      END DO
      Lat(j)=90.0_r8-CoLatRad*rd
      Lat(Jmax-j+1)=-Lat(j)
      CoLatRad=CoLatRad+dCoLatRadz
   END DO

END SUBROUTINE GetGaussianLatitudes


SUBROUTINE LegendrePolynomial (N, CoLatRad, LegPol)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N ! Order of the Ordinary Legendre Function

   REAL (KIND=r8), INTENT(IN) :: CoLatRad ! Colatitude (In Radians)

   REAL (KIND=r8), INTENT(OUT) :: LegPol ! Value of The Ordinary Legendre Function

   INTEGER :: i

   REAL (KIND=r8) :: x, y1, y2, g, y3

   x=COS(CoLatRad)
   y1=1.0_r8
   y2=X
   DO i=2,N
      g=x*y2
      y3=g-y1+g-(g-y1)/REAL(i,r8)
      y1=y2
      y2=y3
   END DO
   LegPol=y3

END SUBROUTINE LegendrePolynomial

SUBROUTINE GetLongitudes (Imax, Lon0, Lon)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: Imax ! Number of Longitudes

   REAL (KIND=r8), INTENT(In) :: Lon0 ! First Longitude In Degree

   REAL (KIND=r8), DIMENSION (Imax), INTENT(OUT) :: &
                   Lon ! Longitudes In Degree

   INTEGER :: i

   REAL (KIND=r8) :: dx

   dx=360.0_r8/REAL(Imax,r8)
   DO i=1,Imax
      Lon(i)=Lon0+REAL(i-1,r8)*dx
   END DO

END SUBROUTINE GetLongitudes


END MODULE MiscMod
