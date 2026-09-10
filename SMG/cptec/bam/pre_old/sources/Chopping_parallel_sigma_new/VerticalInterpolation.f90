!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE VerticalInterpolation

  USE InputParameters, ONLY: r8, nfprt, Gama, Grav, Rd, Lc, Rv,Cp

  IMPLICIT NONE
 
  PRIVATE

  PUBLIC :: VertSigmaInter

CONTAINS


SUBROUTINE VertSigmaInter (ibdim,Im, Km1, Km2, Nt, &
                           p1, u1, v1, t1, q1, &
                           p2, u2, v2, t2, q2)

  ! From NCEP Early 2003
  ! Vertically Interpolate Upper-Air Fields:

  ! Wind, Temperature, Humidity and other Tracers are Interpolated.
  ! The Interpolation is Cubic Lagrangian in Log Pressure
  ! with a Monotonic Constraint in the Center of the Domain.
  ! In the Outer Intervals it is Linear in Log Pressure.
  ! Outside the Domain, Fields are Generally Held Constant,
  ! Except for Temperature and Humidity Below the Input DOmain,
  ! Where the Temperature Lapse Rate is Held Fixed at -6.5 K/km
  ! and the Relative Humidity is Held Constant.

  ! Input Argument List:
  !   Im           First Dimension
  !   Km1          Number of Input Levels
  !   Km2          Number of Output Levels
  !   Nt           Number of Tracers
  !   p1           Input Pressures
  !                (Ordered from Bottom to Top of Atmosphere)
  !   u1           Input Zonal Wind
  !   v1           Input Meridional Wind
  !   t1           Input Temperature (K)
  !   q1           Input Tracers (Specific Humidity First)
  !   p2           Output Pressures

  ! Output Argument List:
  !   u2            Output Zonal Wind
  !   v2            Output Meridional Wind
  !   t2            Output Temperature (K)
  !   q2            Output Tracers (Specific Humidity First)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: ibdim, Im, Km1, Km2, Nt

  REAL (KIND=r8), DIMENSION (Ibdim,Km1), INTENT (IN) :: p1, u1, v1, t1

  REAL (KIND=r8), DIMENSION (Ibdim,Km1,Nt), INTENT (IN) :: q1

  REAL (KIND=r8), DIMENSION (Ibdim,Km2), INTENT (IN) :: p2

  REAL (KIND=r8), DIMENSION (Ibdim,Km2), INTENT (OUT) :: u2, v2, t2

  REAL (KIND=r8), DIMENSION (Ibdim,Km2,Nt), INTENT (OUT) :: q2

  INTEGER :: k, i, n

  REAL (KIND=r8) :: dltdz, dlpvdrt, dz,RdByCp

  REAL (KIND=r8), DIMENSION (Ibdim,Km1) :: z1

  REAL (KIND=r8), DIMENSION (Ibdim,Km2) :: z2

  REAL (KIND=r8), DIMENSION (Ibdim,Km1,3+Nt) :: c1

  REAL (KIND=r8), DIMENSION (Ibdim,Km2,3+Nt) :: c2
  
  RdByCp=Rd/Cp

  dltdz=Gama*Rd/Grav
  dlpvdrt=-Lc/Rv

  ! Compute Log Pressure Interpolating Coordinate and
  ! Copy Input Wind, Temperature, Humidity and other Tracers

  DO k=1,Km1
    DO i=1,Im
      z1(i,k)=-LOG(p1(i,k))
      c1(i,k,1)=u1(i,k)
      c1(i,k,2)=v1(i,k)
      c1(i,k,3)=t1(i,k)
      c1(i,k,4)=q1(i,k,1)
    ENDDO
  ENDDO
  DO n=2,Nt
    DO k=1,Km1
      DO i=1,Im
        c1(i,k,3+n)=q1(i,k,n)
      ENDDO
    ENDDO
  ENDDO
  DO k=1,Km2
    DO i=1,Im
      z2(i,k)=-LOG(p2(i,k))
    ENDDO
  ENDDO

  ! Perform Lagrangian One-Dimensional Interpolation that is
  ! 4th-Order in Interior, 
  ! 2nd-Order in Outside Intervals and
  ! 1st-Order for Extrapolation

  CALL terp3 (ibdim,Im, Km1, Km2, 3+Nt, z1, c1, z2, c2)

  ! Copy Output Wind, Temperature, Specific Humidity and other Tracers
  ! Except Below the Input Domain, Let Temperature Increase with a
  ! Fixed Lapse Rate and Let the Relative Humidity Remain Constant

  DO k=1,Km2
    DO i=1,Im
      u2(i,k)=c2(i,k,1)
      v2(i,k)=c2(i,k,2)
      dz=z2(i,k)-z1(i,1)
      IF (dz >= 0.0_r8) THEN
        t2(i,k)  =c2(i,k,3)
        q2(i,k,1)=c2(i,k,4)
      ELSE
      !  dltdz=Gama*Rd/Grav
      !  Gama=-6.5E-3_r8
        IF(Km2 == Km1)THEN
           IF(k > 1.and.k < Km2)THEN
             ! t2(i,k)  =t1(i,1  )*EXP(0.5_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
             ! q2(i,k,1)=q1(i,1,1)*EXP(0.5_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
              t2(i,k)  =0.25_r8*c2(i,k-1,3)+0.5_r8*c2(i,k,3)+0.25_r8*c2(i,k+1,3)
              q2(i,k,1)=0.25_r8*c2(i,k-1,4)+0.5_r8*c2(i,k,4)+0.25_r8*c2(i,k+1,4)
           ELSE IF(k == 1)THEN
              t2(i,k)  =t1(i,1)!*EXP(dltdz*dz)
              q2(i,k,1)=q1(i,1,1)!*EXP(dlpvdrt*(1.0_r8/t2(i,k)-1.0_r8/t1(i,1))-dz)
           ELSE IF(k == Km2)THEN
              t2(i,k)  =t1(i,k)  !*EXP(dltdz*dz)
              q2(i,k,1)=q1(i,k,1)!*EXP(dlpvdrt*(1.0_r8/t2(i,k)-1.0_r8/t1(i,1))-dz)
           END IF
        ELSE
           IF(k > 1.and.k < Km2)THEN
             ! t2(i,k)  =t1(i,1  )*EXP(0.5_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
             ! q2(i,k,1)=q1(i,1,1)*EXP(0.5_r8*RdByCp*ABS(log(p1(i,k))-log(p2(i,k))))  
              t2(i,k)  =0.25_r8*c2(i,k-1,3)+0.5_r8*c2(i,k,3)+0.25_r8*c2(i,k+1,3)
              q2(i,k,1)=0.25_r8*c2(i,k-1,4)+0.5_r8*c2(i,k,4)+0.25_r8*c2(i,k+1,4)
           ELSE IF(k == 1)THEN
              t2(i,k)  =t1(i,1)!*EXP(dltdz*dz)
              q2(i,k,1)=q1(i,1,1)!*EXP(dlpvdrt*(1.0_r8/t2(i,k)-1.0_r8/t1(i,1))-dz)
           ELSE
              t2(i,k)  =t1(i,1)*EXP(dltdz*dz)
              q2(i,k,1)=q1(i,1,1)*EXP(dlpvdrt*(1.0_r8/t2(i,k)-1.0_r8/t1(i,1))-dz)
           END IF
        END IF
      ENDIF
    ENDDO
  ENDDO
!  DO n=2,nT
!    DO k=1,Km2
!      DO i=1,Im
!        q2(i,k,n)=c2(i,k,3+n)
!      ENDDO
!    ENDDO
!  ENDDO

  ! Copy Output Tracers
 
  DO n=2,nT
    DO k=1,Km2
      DO i=1,Im
        dz=z2(i,k)-z1(i,1)
        IF (dz >= 0.0_r8) THEN
           q2(i,k,n)=c2(i,k,3+n)
        ELSE
           IF(Km2 == Km1)THEN
              IF(k > 1.and.k < Km2)THEN
                 q2(i,k,n)=0.25_r8*c2(i,k-1,3+n)+0.5_r8*c2(i,k,3+n)+0.25_r8*c2(i,k+1,3+n)
              ELSE IF(k == 1)THEN
                 q2(i,k,n)=c1(i,1,3+n)
              ELSE IF(k == Km2)THEN
                 q2(i,k,n)=c1(i,k,3+n)
              END IF
           ELSE
              IF(k > 1.and.k < Km2)THEN
                 q2(i,k,n)=0.25_r8*c2(i,k-1,3+n)+0.5_r8*c2(i,k,3+n)+0.25_r8*c2(i,k+1,3+n)
              ELSE IF(k == 1)THEN
                 q2(i,k,n)=c1(i,1,3+n)
              ELSE
                 q2(i,k,n)=(c2(i,k,3+n)+c1(i,Km1,3+n))/2.0_r8
              END IF
           END IF 
        END IF
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE VertSigmaInter


SUBROUTINE terp3 (ibdim,Im, Km1, Km2, Nm, z1, q1, z2, q2)

  ! From NCEP Early 2003
  ! Cubically Interpolate Field(s) in One Dimension Along the Column(s).
  ! The Interpolation is Cubic Lagrangian with a Monotonic Constraint
  ! in the Center of the Domain. In the Outer Intervals it is Linear.
  ! Outside the Domain, Fields Are Held Constant.

  ! Input Argument List:
  !   Im           Number of Columns
  !   Km1          Number of Input Points in Each Column
  !   Km2          Number of Output Points in Each Column
  !   Nm           Number of Fields per Column
  !   z1           Input Coordinate Values in which to Interpolate
  !                (z1 Must Be Strictly Monotonic in Either Direction)
  !   q1           Input Fields to Interpolate
  !   z2           Output Coordinate Values to which to Interpolate
  !                (z2 Need Not Be Monotonic)
    
  ! Output Argument List:
  !   q2            Output Interpolated Fields

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: ibdim,Im, Km1, Km2, Nm

  REAL (KIND=r8), DIMENSION (ibdim,Km1), INTENT (IN) :: z1

  REAL (KIND=r8), DIMENSION (ibdim,Km2), INTENT (IN) :: z2

  REAL (KIND=r8), DIMENSION (ibdim,Km1,Nm), INTENT (IN) :: q1

  REAL (KIND=r8), DIMENSION (ibdim,Km2,Nm), INTENT (OUT) :: q2

  INTEGER :: k1, k2, i, n

  REAL (KIND=r8) :: z1a, z1b, z1c, z1d, z2s, q1a, q1b, q1c, q1d, q2s

  INTEGER, DIMENSION (Ibdim,Km2) :: k1s

  REAL (KIND=r8), DIMENSION (Ibdim) :: ffa, ffb, ffc, ffd

  ! Find the Surrounding Input Interval for Each Output Point

  CALL rsearch (ibdim, Im, Km1, Km2, z1, z2, k1s)

  ! Generally Interpolate Cubically with Monotonic Constraint
  ! From Two Nearest Input Points on Either Side of the Output Point,
  ! But Within the Two Edge Intervals Interpolate Linearly.
  ! Keep the Output Fields Constant Outside the Input Domain.

  DO k2=1,Km2
    DO i=1,Im
      k1=k1s(i,k2)
      IF (k1 == 1 .OR. k1 == Km1-1) THEN
        z2s=z2(i,k2)
        z1a=z1(i,k1)
        z1b=z1(i,k1+1)
        ffa(i)=(z2s-z1b)/(z1a-z1b)
        ffb(i)=(z2s-z1a)/(z1b-z1a)
      ELSEIF (k1 > 1 .AND. k1 < Km1-1) THEN
        z2s=z2(i,k2)
        z1a=z1(i,k1-1)
        z1b=z1(i,k1)
        z1c=z1(i,k1+1)
        z1d=z1(i,k1+2)
        ffa(i)=(z2s-z1b)/(z1a-z1b)*(z2s-z1c)/(z1a-z1c)*(z2s-z1d)/(z1a-z1d)
        ffb(i)=(z2s-z1a)/(z1b-z1a)*(z2s-z1c)/(z1b-z1c)*(z2s-z1d)/(z1b-z1d)
        ffc(i)=(z2s-z1a)/(z1c-z1a)*(z2s-z1b)/(z1c-z1b)*(z2s-z1d)/(z1c-z1d)
        ffd(i)=(z2s-z1a)/(z1d-z1a)*(z2s-z1b)/(z1d-z1b)*(z2s-z1c)/(z1d-z1c)
      ENDIF
    ENDDO

  ! Interpolate

    DO n=1,Nm
      DO i=1,Im
        k1=k1s(i,k2)
        IF (k1 == 0) THEN
          q2s=q1(i,1,n)
        ELSEIF (k1 == Km1) THEN
          q2s=q1(i,Km1,n)
        ELSEIF (k1 == 1 .OR. k1 == Km1-1) THEN
          q1a=q1(i,k1,n)
          q1b=q1(i,k1+1,n)
          q2s=ffa(i)*q1a+ffb(i)*q1b
        ELSE
          q1a=q1(i,k1-1,n)
          q1b=q1(i,k1,n)
          q1c=q1(i,k1+1,n)
          q1d=q1(i,k1+2,n)
          q2s=MIN(MAX(ffa(i)*q1a+ffb(i)*q1b+ffc(i)*q1c+ffd(i)*q1d, &
              MIN(q1b,q1c)),MAX(q1b,q1c))
        ENDIF
        q2(i,k2,n)=q2s
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE terp3


SUBROUTINE rsearch (ibdim, Im, Km1, Km2, z1, z2, l2)

  ! From NCEP Early 2003
  ! Search for a Surrounding Real Interval

  ! Searches Monotonic Sequences of Real Numbers for Intervals
  ! that Surround a Given Search SeT of Real Numbers
  ! The Sequences Must Be Monotonically Ascending
  ! The Input Sequences and Sets and the Output Locations 
  ! May Be Arbitrarily Dimensioned

  ! Input Argument List:
  !   Im           Number of Sequences to Search
  !   Km1          Number of Points in Each Sequence
  !   Km2          Number of Points to Search For
  !                in Each Respective Sequence
  !   z1           Sequence Values to Search
  !                (z1 Must Be Monotonically Ascending)
  !   z2           Set of Values to Search For
  !                (z2 Need Not Be Monotonic)

  ! Output Argument List:
  !   l2            Interval Locations Having Values from 0 to Km1
  !                 (z2 Will Be Between z1(l2) and z1(l2+1))

  ! Remarks:

  ! Returned Values of 0 or Km1 Indicate That The Given Search Value
  ! is Outside the Range of the Sequence.

  ! If a Search Value is Identical to One of the Sequence Values
  ! Then the Location Returned Points to the Identical value.
  ! If the Sequence Is Not Strictly Monotonic and a Search Value Is
  ! Identical to More Than One of the Sequence Values, Then the
  ! Location Returned May Point to Any of the Identical Values.

  ! To Be Exact, For Each i From 1 To Im and For Each k From 1 To Km2,
  ! z=z2(i,k) is the Search Value and l=l2(i,k) is the Location Returned.
  ! If l=0, Then z is Less Than the Start Point z1(i,1).
  ! If l=Km1, Then z is Greater Than or Equal to the End Point
  ! z1(i,Km1).
  ! Otherwise z is Between the Values z1(i,l) And z1(i,l+1) and 
  ! May Equal the Former.

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibdim, Im, Km1, Km2

  REAL (KIND=r8), DIMENSION (Ibdim,Km1), INTENT(IN) :: z1

  REAL (KIND=r8), DIMENSION (Ibdim,Km2), INTENT(IN) :: z2

  INTEGER, DIMENSION (Ibdim,Km2), INTENT(OUT) :: l2

  INTEGER :: i, k2

  INTEGER, DIMENSION (Km2) :: indx, rc

  ! Find the Surrounding Input Interval for Each Output Point

  DO i=1,Im

    IF (z1(i,1) <= z1(i,Km1)) THEN

    ! Input Coordinate is Monotonically Ascending

      CALL bsrch (Km1, Km2, z1(i,:), z2(i,:), indx, rc)

      DO k2=1,Km2
        l2(i,k2)=indx(k2)-rc(k2)
      ENDDO

    ELSE

    ! Input Coordinate is Monotonically Descending

    WRITE (UNIT=nfprt, FMT='(/,A)') ' Warnning: '
    WRITE (UNIT=nfprt, FMT='(A)')   ' Input Coordinate is Monotonically Descending'
    WRITE (UNIT=nfprt, FMT='(A)')   ' The Implemented Binary Search Does Not Allowed That'
    WRITE (UNIT=nfprt, FMT='(A,/)') ' Stopping Computation at SUBROUTINE rsearch '
    STOP

    ENDIF

  ENDDO

END SUBROUTINE rsearch


SUBROUTINE bsrch (n, m, x, y, indx, rc)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: m, n

  REAL (KIND=r8), DIMENSION (n), INTENT (IN) :: x
  REAL (KIND=r8), DIMENSION (m), INTENT (IN) :: y

  INTEGER, DIMENSION (m), INTENT (OUT) :: indx, rc

  INTEGER :: i, j

  out: DO j=1,m

    IF (y(j) < x(1)) THEN
      indx(j)=1
      rc(j)=1
      CYCLE out
    END IF
    IF (y(j) > x(n)) THEN
      indx(j)=n+1
      rc(j)=1
      CYCLE out
    END IF
    DO i=1,n
      IF (y(j) == x(i)) THEN
        indx(j)=i
        rc(j)=0
        CYCLE out
      END IF
    END DO
    DO i=1,n-1
      IF (y(j) > x(i) .AND. y(j) < x(i+1)) THEN
        indx(j)=i+1
        rc(j)=1
        CYCLE out
      END IF
    END DO

  END DO out

END SUBROUTINE bsrch


END MODULE VerticalInterpolation
