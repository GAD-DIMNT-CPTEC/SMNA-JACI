!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE LinearInterpolation

   ! Logics Assumes that Input and Output Data 
   ! First Point is Near North Pole and Greenwhich

   USE InputParameters, ONLY: r4, r8, Idim, Jdim, Imax, Jmax, &
                              Lon0, Lat0, Undef

   IMPLICIT NONE

   INTEGER, DIMENSION (:), ALLOCATABLE :: &
            LowerLon, UpperLon, LowerLat, UpperLat

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: &
         LonIn, LatIn, LonOut, LatOut

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: &
         LeftLowerWgt, LeftUpperWgt, &
         RightLowerWgt, RightUpperWgt

   PUBLIC :: InitLinearInterpolation, DoLinearInterpolation


CONTAINS


SUBROUTINE InitLinearInterpolation ()

   IMPLICIT NONE

   ALLOCATE (LonIn(Idim), LatIn(Jdim))
   ALLOCATE (LonOut(Imax), LatOut(Jmax))
   ALLOCATE (LowerLon(Imax), LowerLat(Jmax))
   ALLOCATE (UpperLon(Imax), UpperLat(Jmax))
   ALLOCATE (LeftLowerWgt(Imax,Jmax))
   ALLOCATE (LeftUpperWgt(Imax,Jmax))
   ALLOCATE (RightLowerWgt(Imax,Jmax))
   ALLOCATE (RightUpperWgt(Imax,Jmax))

   CALL GetLongitudes (Idim, Lon0, LonIn)
   CALL GetLongitudes (Imax, 0.0_r8, LonOut)
   CALL GetRegularLatitudes (Jdim, Lat0, LatIn)
   CALL GetGaussianLatitudes (Jmax, LatOut)

   CALL HorizontalInterpolationWeights ()

END SUBROUTINE InitLinearInterpolation


SUBROUTINE DoLinearInterpolation (VarIn, VarOut)

   IMPLICIT NONE

   REAL (KIND=r8), INTENT (IN), DIMENSION (Idim,Jdim) :: VarIn

   REAL (KIND=r8), INTENT (OUT), DIMENSION (Imax,Jmax) :: VarOut


   CALL HorizontalInterpolation (VarIn, VarOut)

END SUBROUTINE DoLinearInterpolation


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


SUBROUTINE GetRegularLatitudes (Jmax, Lat0, Lat)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: Jmax ! Number of Regular Latitudes

   REAL (KIND=r8), INTENT(In) :: Lat0 ! First Latitude In Degree

   REAL (KIND=r8), DIMENSION (Jmax), INTENT(OUT) :: &
                   Lat ! Regular Latitudes In Degree

   INTEGER ::j 

   REAL (KIND=r8) :: dy

   dy=2.0_r8*Lat0/REAL(Jmax-1,r8)
   DO j=1,Jmax
      Lat(j)=Lat0-REAL(j-1,r8)*dy
   END DO

END SUBROUTINE GetRegularLatitudes


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


SUBROUTINE HorizontalInterpolationWeights ()

   IMPLICIT NONE

   INTEGER :: i, j, iLon, jLat

   REAL (KIND=r8) :: dLon, LowerLonIn, UpperLonIn, dbx, LatOutA, &
                     dlat, LowerLatIn, UpperLatIn, dby

   INTEGER, DIMENSION (1) :: iLonm, jLatm

   REAL (KIND=r8), DIMENSION (Imax) :: dx, ddx

   REAL (KIND=r8), DIMENSION (Jmax) :: dy, ddy

   DO i=1,Imax
      iLonm=MINLOC(ABS(LonIn(1:Idim)-LonOut(i)))
      iLon=iLonm(1)
      dLon=LonOut(i)-LonIn(iLon)
      IF (dLon < 0.0_r8) THEN
         IF (iLon /= 1) THEN
            LowerLon(i)=iLon-1
            LowerLonIn=LonIn(iLon-1)
         ELSE
            LowerLon(i)=Idim
            LowerLonIn=LonIn(Idim)-360.0_r8
         END IF
         UpperLon(i)=iLon
         UpperLonIn=LonIn(iLon)
      ELSE IF (dLon == 0.0_r8) THEN
         LowerLon(i)=iLon
         LowerLonIn=LonIn(iLon)
         UpperLon(i)=iLon
         UpperLonIn=LonIn(iLon)
      ELSE
         LowerLon(i)=iLon
         LowerLonIn=LonIn(iLon)
         IF (iLon /= Idim) THEN
            UpperLon(i)=iLon+1
            UpperLonIn=LonIn(iLon+1)
         ELSE
            UpperLon(i)=1
            UpperLonIn=LonIn(1)+360.0_r8
         END IF
      END IF
      dx(i)=LonOut(i)-LowerLonIn
      ddx(i)=UpperLonIn-LowerLonIn
   END DO

   DO j=1,Jmax
      jLatm=MINLOC(ABS(LatIn(1:Jdim)-LatOut(j)))
      jLat=jLatm(1)
      LatOutA=LatOut(j)
      IF (LatOutA > LatIn(1)) LatOutA=LatIn(1)
      IF (LatOutA < LatIn(Jdim)) LatOutA=LatIn(Jdim)
      dlat=LatOutA-LatIn(jLat)
      IF (dlat > 0.0_r8) THEN
         IF (jLat /= 1) THEN
            LowerLat(j)=jLat-1
            LowerLatIn=LatIn(jLat-1)
         ELSE
            LowerLat(j)=1
            LowerLatIn=LatIn(1)
         END IF
         UpperLat(j)=jLat
         UpperLatIn=LatIn(jLat)
      ELSE IF (dlat == 0.0_r8) THEN
         LowerLat(j)=jLat
         LowerLatIn=LatIn(jLat)
         UpperLat(j)=jLat
         UpperLatIn=LatIn(jLat)
      ELSE
         LowerLat(j)=jLat
         LowerLatIn=LatIn(jLat)
         IF (jLat /= Jdim) THEN
            UpperLat(j)=jLat+1
            UpperLatIn=LatIn(jLat+1)
         ELSE
            UpperLat(j)=Jdim
            UpperLatIn=LatIn(Jdim)
         END IF
      END IF
      dy(j)=LatOutA-LowerLatIn
      ddy(j)=UpperLatIn-LowerLatIn
   END DO

   DO j=1,Jmax
      DO i=1,Imax
         IF (ddx(i) == 0.0_r8 .AND. ddy(j) == 0.0_r8) THEN
            LeftLowerWgt(i,j)=1.0_r8
            LeftUpperWgt(i,j)=0.0_r8
            RightLowerWgt(i,j)=0.0_r8
            RightUpperWgt(i,j)=0.0_r8
         ELSE IF (ddx(i) == 0.0_r8) THEN
            LeftUpperWgt(i,j)=dy(j)/ddy(j)
            LeftLowerWgt(i,j)=1.0_r8-LeftUpperWgt(i,j)
            RightLowerWgt(i,j)=0.0_r8
            RightUpperWgt(i,j)=0.0_r8
         ELSE IF (ddy(j) == 0.0_r8) THEN
            RightLowerWgt(i,j)=dx(i)/ddx(i)
            LeftLowerWgt(i,j)=1.0_r8-RightLowerWgt(i,j)
            LeftUpperWgt(i,j)=0.0_r8
            RightUpperWgt(i,j)=0.0_r8
         ELSE
            dbx=dx(i)/ddx(i)
            dby=dy(j)/ddy(j)
            RightUpperWgt(i,j)=dbx*dby
            LeftLowerWgt(i,j)=1.0_r8-dbx-dby+RightUpperWgt(i,j)
            LeftUpperWgt(i,j)=dby-RightUpperWgt(i,j)
            RightLowerWgt(i,j)=dbx-RightUpperWgt(i,j)
         END IF
      END DO
   END DO

END SUBROUTINE HorizontalInterpolationWeights


SUBROUTINE HorizontalInterpolation (VarIn, VarOut)

!  i - Longitude Index
!  j - Latitude Index

!  Input Grid Box That Contains The Output Value (i,j):
!  LowerLon(i,j) - Lower Input Longitude Index
!  UpperLon(i,j) - Upper Input Longitude Index
!  LowerLat(i,j) - Lower Input Latitude  Index
!  UpperLat(i,j) - Upper Input Latitude  Index

!  Pre-Calculated Weights for Linear Horizontal Interpolation:
!  LeftLowerWgt  - for Left-Lower  Corner of Box
!  LeftUpperWgt  - for Left-Upper  Corner of Box
!  RightLowerWgt - for Right-Lower Corner of Box
!  RightUpperWgt - for Right-Upper Corner of Box

   IMPLICIT NONE 
   REAL (KIND=r8), INTENT (IN), DIMENSION (Idim,Jdim) :: VarIn

   REAL (KIND=r8), INTENT (OUT), DIMENSION (Imax,Jmax) :: VarOut

   INTEGER :: j, i, il, iu, jl, ju

   DO j=1,Jmax
      DO i=1,Imax
         il=LowerLon(i)
         jl=LowerLat(j)
         iu=UpperLon(i)
         ju=UpperLat(j)
         IF (VarIn(il,jl) == Undef .AND. VarIn(il,ju) == Undef .AND. &
             VarIn(iu,jl) == Undef .AND. VarIn(iu,ju) == Undef) THEN
            Varout(i,j)=Undef
         ELSE
         Varout(i,j)=0.0_r8
            IF (VarIn(il,jl) /= Undef) &
               VarOut(i,j)=VarOut(i,j)+LeftLowerWgt(i,j)*VarIn(il,jl)
            IF (VarIn(il,ju) /= Undef) &
               VarOut(i,j)=VarOut(i,j)+LeftUpperWgt(i,j)*VarIn(il,ju)
            IF (VarIn(iu,jl) /= Undef) &
               VarOut(i,j)=VarOut(i,j)+RightLowerWgt(i,j)*VarIn(iu,jl)
            IF (VarIn(iu,ju) /= Undef) &
               VarOut(i,j)=VarOut(i,j)+RightUpperWgt(i,j)*VarIn(iu,ju)
         END IF
      END DO
   END DO

END SUBROUTINE HorizontalInterpolation


END MODULE LinearInterpolation
