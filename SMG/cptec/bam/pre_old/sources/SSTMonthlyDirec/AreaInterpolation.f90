!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
MODULE AreaInterpolation

  USE InputParameters, ONLY : r4, r8, IdimIn=>Idim, JdimIn=>Jdim, &
                              IdimOut=>Imax, JdimOut=>Jmax, Undef, &
                              PolarMean, FlagInput, Flagoutput, MaskInput

  ! IdimIn  number of longitude points for the input grid
  ! JdimIn  number of latitude points for the input grid
  ! IdimOut number of longitude points for the output grid
  ! JdimOut number of latitude points for the output grid
  ! Undef  Undefined value which if found in input array causes
  !        that location to be ignored in interpolation.  Used
  !        as the output value for output points with no defined
  !        and/or unmasked data
  ! PolarMean flag to performe average at poles
  ! FlagInput:  input  grid flags 
  ! FlagOutput: output grid flags 
  ! flags: (input or output)
  !   1   start at north pole (true) start at south pole (false)
  !   2   start at prime meridian (true) start at i.d.l. (false)
  !   3   latitudes are at center of box (true)
  !       latitudes are at edge (false) north edge if 1=true
  !                                     south edge if 1=false
  !   4   longitudes are at center of box (true)
  !       longitudes are at western edge of box (false)
  !   5   gaussian (true) regular (false)
  ! MaskInput  input grid mask to confine interpolation of input
  !        data to certain areas (1=interpolate, 0=don't)

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitAreaInterpolation, DoAreaInterpolation

  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE, PUBLIC :: gLats

  ! lons   longitude dimension for weights
  ! lats   latitude dimension for weights
  ! lond   total number of longitude weights
  ! latd   total number of latitude weigths
  ! lwrk   dimension for working area
  INTEGER :: lons, lats, lond, latd, lwrk

  ! mplon  longitude index mapping from input (,1) to output (,2)
  ! mplat  latitude index mapping from input (,1) to output (,2)
  INTEGER, DIMENSION (:,:), ALLOCATABLE :: mplon, mplat

  ! wtlon  area weights in the longitudinal direction
  ! wtlat  area weights in the latitudianl direction
  REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: wtlon, wtlat


CONTAINS


  SUBROUTINE InitAreaInterpolation ()

    IMPLICIT NONE

    lons=IdimIn+IdimOut+2
    lats=JdimIn+JdimOut+2
    lwrk=MAX(2*lons,2*lats)

    IF (ALLOCATED(mplon)) DEALLOCATE(mplon)
    IF (ALLOCATED(mplat)) DEALLOCATE(mplat)
    ALLOCATE (mplon(lons,2), mplat(lats,2))
    IF (ALLOCATED(wtlon)) DEALLOCATE(wtlon)
    IF (ALLOCATED(wtlat)) DEALLOCATE(wtlat)
    ALLOCATE (wtlon(lons), wtlat(lats))

    ! Computing Weights For Horizontal Interpolation

    CALL GetAreaInterpolationWeights ()

  END SUBROUTINE InitAreaInterpolation


  SUBROUTINE DoAreaInterpolation (FieldInput, FieldOutput)

    ! Interpolation Subroutine.
    ! Should only be called after Subroutine GetAreaInterpolationWeights 
    ! has been called.
    ! After interpolation user is responsible for output masking
    ! and pole interpolation for type 2 (pole centered) grids.
    ! Pole interpolation can be done with PolarMean=.true.
    ! 
    ! subroutine arguments:
    ! 
    ! FieldInput input field to be interpolated
    ! FieldOutput output field resulting from interpolation

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: FieldInput(IdimIn,JdimIn)

    REAL (KIND=r8), INTENT(OUT) :: FieldOutput(IdimOut,JdimOut)

    INTEGER :: j, i, lti, lto, lni, lno

    REAL (KIND=r8) :: wlt, wln

    REAL (KIND=r8) :: work(IdimOut,JdimOut)

    FieldOutput=0.0_r8
    work=0.0_r8

    DO j=1,latd
       wlt=wtlat(j)
       lti=mplat(j,1)
       lto=mplat(j,2)
       DO i=1,lond
          lni=mplon(i,1)
          IF (MaskInput(lni,lti) /= 0) THEN
             IF (FieldInput(lni,lti) /= Undef) THEN
                wln=wtlon(i)
                lno=mplon(i,2)
                FieldOutput(lno,lto)=FieldOutput(lno,lto)+FieldInput(lni,lti)*wlt*wln
                work(lno,lto)=work(lno,lto)+wlt*wln
             END IF
          END IF
       END DO
    END DO

    DO j=1,JdimOut
       DO i=1,IdimOut
          IF (work(i,j) == 0.0_r8)THEN
             FieldOutput(i,j)=Undef
          ELSE
             FieldOutput(i,j)=FieldOutput(i,j)/work(i,j)
          END IF
       END DO
    END DO

    IF (PolarMean) THEN
      FieldOutput(:,1)=SUM(FieldOutput(:,1))/REAL(IdimOut,r8)
      FieldOutput(:,JdimOut)=SUM(FieldOutput(:,JdimOut))/REAL(IdimOut,r8)
    END IF

  END SUBROUTINE DoAreaInterpolation


  SUBROUTINE GetAreaInterpolationWeights ()

    ! Interpolation Weight Calculation
    ! This Subroutine should be called once to determine the area
    ! weights and index mapping between a pair of grids on a sphere.
    ! The weights and map indices are used by subroutine
    ! DoAreaInterpolation to perfom the actual interpolation.

    IMPLICIT NONE

    REAL (KIND=r8) :: WorkWeights(lwrk)

    INTEGER :: lath, joi, joo, j, j1, j2, j3, &
                         ioi, ici, i, ioo, ico, i1, i2, i3

    REAL (KIND=r4) :: eps
    REAL (KIND=r8) :: dpi, rad, drltm, drltp, dlat, dof, delrdi, delrdo

    dpi=4.0_r8*ATAN(1.0_r8)
    rad=180.0_r8/dpi
    eps=EPSILON(1.0_r4)

    ! input grid latitudes

    joi=JdimIn+JdimOut+2
    IF (FlagInput(5)) THEN

       ! gaussian grid case

       lath=JdimIn/2
       CALL GaussianLatitudes (lath, WorkWeights)
       DO j=2,JdimIn
          IF (j <= lath) THEN
             drltm=-dpi/2.0_r8+WorkWeights(j-1)
             drltp=-dpi/2.0_r8+WorkWeights(j)
          ELSE IF(j > lath+1)THEN
             drltm=dpi/2.0_r8-WorkWeights(JdimIn-j+2)
             drltp=dpi/2.0_r8-WorkWeights(JdimIn-j+1)
          ELSE
             drltm=0.0_r8
             drltp=0.0_r8
          END IF
          WorkWeights(j+joi)=SIN((drltm+drltp)/2.0_r8)
       END DO
       WorkWeights(1+joi)=-1.0_r8
       WorkWeights(lath+1+joi)=0.0_r8
       WorkWeights(JdimIn+1+joi)=1.0_r8
    ELSE

       ! regular grid case

       IF (FlagInput(3)) THEN
          dlat=dpi/REAL(JdimIn-1,r8)
          dof=-(dpi+dlat)/2.0_r8
       ELSE
          dlat=dpi/REAL(JdimIn,r8)
          dof=-dpi/2.0_r8
       END IF
       DO j=2,JdimIn
          WorkWeights(joi+j)=SIN(dof+dlat*REAL(j-1,r8))
       END DO
       WorkWeights(1+joi)=-1.0_r8
       WorkWeights(JdimIn+1+joi)=1.0_r8
    END IF

    ! output grid latitudes

    joo=2*JdimIn+JdimOut+3

    IF (FlagOutput(5)) THEN

       ! gaussian grid case

       lath=JdimOut/2
       CALL GaussianLatitudes (lath, WorkWeights)
       IF (ALLOCATED(gLats)) DEALLOCATE(gLats)
       ALLOCATE (gLats(JdimOut))
       DO j=1,lath
         gLats(j)=90.0_r8-rad*WorkWeights(j)
         gLats(JdimOut-j+1)=-gLats(j)
       END DO
       DO j=2,JdimOut
          IF (j <= lath) THEN
             drltm=-dpi/2.0_r8+WorkWeights(j-1)
             drltp=-dpi/2.0_r8+WorkWeights(j)
          ELSE IF (j > lath+1) THEN
             drltm=dpi/2.0_r8-WorkWeights(JdimOut-j+2)
             drltp=dpi/2.0_r8-WorkWeights(JdimOut-j+1)
          ELSE
             drltm=0.0_r8
             drltp=0.0_r8
          END IF
          WorkWeights(j+joo)=SIN((drltm+drltp)/2.0_r8)
       END DO
       WorkWeights(1+joo)=-1.0_r8
       WorkWeights(lath+1+joo)=0.0_r8
       WorkWeights(JdimOut+1+joo)=1.0_r8
    ELSE

       ! regular grid case

       IF (FlagOutput(3)) THEN
          dlat=dpi/REAL(JdimOut-1,r8)
          dof=-(dpi+dlat)/2.0_r8
       ELSE
          dlat=dpi/REAL(JdimOut,r8)
          dof=-dpi/2.0_r8
       END IF
       IF (ALLOCATED(gLats)) DEALLOCATE(gLats)
       ALLOCATE (gLats(JdimOut))
       DO j=1,JdimOut
         gLats(j)=90.0_r8-rad*dlat*REAL(j-1,r8)
       END DO
       DO j=2,JdimOut
          WorkWeights(joo+j)=SIN(dof+dlat*REAL(j-1,r8))
       END DO
       WorkWeights(1+joo)=-1.0_r8
       WorkWeights(JdimOut+1+joo)=1.0_r8
    END IF

    ! produce single ordered set of sin(lat) for both grids
    ! determine latitude weighting and index mapping

    j1=1
    j2=1
    j3=1
    DO
       IF (ABS(WorkWeights(j1+joi)-WorkWeights(j2+joo)) < eps) THEN
          WorkWeights(j3)=WorkWeights(j1+joi)
          IF (j3 /= 1) THEN
             wtlat(j3-1)=WorkWeights(j3)-WorkWeights(j3-1)
             mplat(j3-1,1)=j1-1
             IF (FlagInput(1)) mplat(j3-1,1)=JdimIn+2-j1
             mplat(j3-1,2)=j2-1
             IF (FlagOutput(1)) mplat(j3-1,2)=JdimOut+2-j2
          END IF
          j1=j1+1
          j2=j2+1
          j3=j3+1
       ELSE IF (WorkWeights(j1+joi) < WorkWeights(j2+joo)) THEN
          WorkWeights(j3)=WorkWeights(j1+joi)
          IF (j3 /= 1) THEN
             wtlat(j3-1)=WorkWeights(j3)-WorkWeights(j3-1)
             mplat(j3-1,1)=j1-1
             IF (FlagInput(1)) mplat(j3-1,1)=JdimIn+2-j1
             mplat(j3-1,2)=j2-1
             IF (FlagOutput(1)) mplat(j3-1,2)=JdimOut+2-j2
          END IF
          j1=j1+1
          j3=j3+1
       ELSE
          WorkWeights(j3)=WorkWeights(j2+joo)
          IF (j3 /= 1)THEN
             wtlat(j3-1)=WorkWeights(j3)-WorkWeights(j3-1)
             mplat(j3-1,1)=j1-1
             IF (FlagInput(1)) mplat(j3-1,1)=JdimIn+2-j1
             mplat(j3-1,2)=j2-1
             IF (FlagOutput(1)) mplat(j3-1,2)=JdimOut+2-j2
          END IF
          j2=j2+1
          j3=j3+1
       END IF
       IF (.NOT.(j1 <= JdimIn+1 .AND. j2 <= JdimOut+1)) EXIT
    END DO
    latd=j3-2

    ! latitudes done, now do longitudes

    ! input grid longitudes

    ioi=IdimIn+IdimOut+2
    delrdi=(2.0_r8*dpi)/REAL(IdimIn,r8)
    IF (FlagInput(5) .OR. FlagInput(4)) THEN
       ici=0
       dof=0.5_r8
    ELSE
       ici=1
       dof=0.0_r8
    END IF
    DO i=1,IdimIn
       WorkWeights(i+ioi)= (dof+REAL(i-1,r8))*delrdi
    END DO

    ! output grid longitudes

    ioo=2*IdimIn+IdimOut+3
    delrdo=(2.0_r8*dpi)/REAL(IdimOut,r8)
    IF (FlagOutput(5) .OR. FlagOutput(4)) THEN
       ico=0
       dof=0.5_r8
    ELSE
       ico=1
       dof=0.0_r8
    END IF
    DO i=1,IdimOut
       WorkWeights(i+ioo)= (dof+REAL(i-1,r8))*delrdo
    END DO

    ! produce single ordered set of longitudes for both grids
    ! determine longitude weighting and index mapping

    i1=1
    i2=1
    i3=1
    DO
       IF (ABS(WorkWeights(i1+ioi)-WorkWeights(i2+ioo)) < eps) THEN
          WorkWeights(i3)=WorkWeights(i1+ioi)
          IF (i3 /= 1) THEN
             wtlon(i3-1)=WorkWeights(i3)-WorkWeights(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.FlagInput(2)) THEN
                mplon(i3-1,1)=IdimIn/2+i1-ici
                IF(i1-ici > IdimIn/2)mplon(i3-1,1)=i1-ici-IdimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.FlagOutput(2)) THEN
                mplon(i3-1,2)=IdimOut/2+i2-ico
                IF (i2-ico > IdimOut/2) mplon(i3-1,2)=i2-ico-IdimOut/2
             END IF
          END IF
          i1=i1+1
          i2=i2+1
          i3=i3+1
       ELSE IF (WorkWeights(i1+ioi) < WorkWeights(i2+ioo)) THEN
          WorkWeights(i3)=WorkWeights(i1+ioi)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=WorkWeights(i3)-WorkWeights(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.FlagInput(2)) THEN
                mplon(i3-1,1)=IdimIn/2+i1-ici
                IF (i1-ici > IdimIn/2) mplon(i3-1,1)=i1-ici-IdimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.FlagOutput(2)) THEN
                mplon(i3-1,2)=IdimOut/2+i2-ico
                IF (i2-ico > IdimOut/2) mplon(i3-1,2)=i2-ico-IdimOut/2
             END IF
          END IF
          i1=i1+1
          i3=i3+1
       ELSE
          WorkWeights(i3)=WorkWeights(i2+ioo)
          IF (i3 /= 1)THEN
             wtlon(i3-1)=WorkWeights(i3)-WorkWeights(i3-1)
             mplon(i3-1,1)=i1-ici
             IF (.NOT.FlagInput(2)) THEN
                mplon(i3-1,1)=IdimIn/2+i1-ici
                IF (i1-ici > IdimIn/2) mplon(i3-1,1)=i1-ici-IdimIn/2
             END IF
             mplon(i3-1,2)=i2-ico
             IF (.NOT.FlagOutput(2)) THEN
                mplon(i3-1,2)=IdimOut/2+i2-ico
                IF (i2-ico > IdimOut/2) mplon(i3-1,2)=i2-ico-IdimOut/2
             END IF
          END IF
          i2=i2+1
          i3=i3+1
       END IF
       IF (.NOT.(i1 <= IdimIn .AND. i2 <= IdimOut)) EXIT
    END DO
    IF (i1 > IdimIn) i1=1
    IF (i2 > IdimOut) i2=1
    DO
       IF (i2 /= 1) THEN
          WorkWeights(i3)=WorkWeights(i2+ioo)
          wtlon(i3-1)=WorkWeights(i3)-WorkWeights(i3-1)
          mplon(i3-1,1)=1
          IF (.NOT.(FlagInput(4) .OR. FlagInput(5))) mplon(i3-1,1)=IdimIn
          IF (.NOT.FlagInput(2)) THEN
             mplon(i3-1,1)=IdimIn/2+1
             IF (.NOT.(FlagInput(4) .OR. FlagInput(5))) mplon(i3-1,1)=IdimIn/2
          END IF
          mplon(i3-1,2)=i2-ico
          IF (.NOT.FlagOutput(2)) THEN
             mplon(i3-1,2)=IdimOut/2+i2-ico
             IF (i2-ico > IdimOut/2) mplon(i3-1,2)=i2-ico-IdimOut/2
          END IF
          i2=i2+1
          IF (i2 > IdimOut)i2=1
          i3=i3+1
       END IF
       IF (i1 /= 1)THEN
          WorkWeights(i3)=WorkWeights(i1+ioi)
          wtlon(i3-1)=WorkWeights(i3)-WorkWeights(i3-1)
          mplon(i3-1,1)=i1-ici
          IF (.NOT.FlagInput(2)) THEN
             mplon(i3-1,1)=IdimIn/2+i1-ici
             IF (i1-ici > IdimIn/2) mplon(i3-1,1)=i1-ici-IdimIn/2
          END IF
          mplon(i3-1,2)=1
          IF (.NOT.(FlagOutput(4) .OR. FlagOutput(5))) mplon(i3-1,2)=IdimOut
          IF (.NOT.FlagOutput(2)) THEN
             mplon(i3-1,2)=IdimOut/2+1
             IF (.NOT.(FlagOutput(4) .OR. FlagOutput(5))) mplon(i3-1,2)=IdimOut/2
          END IF
          i1=i1+1
          IF (i1 > IdimIn)i1=1
          i3=i3+1
       END IF
       IF (.NOT.(i1 /=1 .OR. i2 /=1)) EXIT
    END DO
    wtlon(i3-1)=2.0_r8*dpi+WorkWeights(1)-WorkWeights(i3-1)
    mplon(i3-1,1)=1
    IF (.NOT.(FlagInput(4) .OR. FlagInput(5))) mplon(i3-1,1)=IdimIn
    IF (.NOT.FlagInput(2)) THEN
       mplon(i3-1,1)=IdimIn/2+1
       IF (.NOT.(FlagInput(4) .OR. FlagInput(5))) mplon(i3-1,1)=IdimIn/2
    END IF
    mplon(i3-1,2)=1
    IF (.NOT.(FlagOutput(4) .OR. FlagOutput(5))) mplon(i3-1,2)=IdimOut
    IF (.NOT.FlagOutput(2)) THEN
       mplon(i3-1,2)=IdimOut/2+1
       IF (.NOT.(FlagOutput(4) .OR. FlagOutput(5))) mplon(i3-1,2)=IdimOut/2
    END IF
    lond=i3-1

  END SUBROUTINE GetAreaInterpolationWeights


  SUBROUTINE GaussianLatitudes (Lath, CoLatitude)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Lath

    REAL (KIND=r8), INTENT(OUT) :: CoLatitude(Lath)

    INTEGER :: Lats, j

    REAL (KIND=r8) :: eps, dGcolIn, Gcol, dGcol, p1, p2

    eps=EPSILON(1.0_r8)*100.0_r8
    Lats=2*Lath
    dGcolIn=ATAN(1.0_r8)/REAL(Lats,r8)
    Gcol=0.0_r8
    DO j=1,Lath
       dGcol=dGcolIn
       DO
          CALL LegendrePolynomial (Lats, Gcol, p2)
          DO
             p1=p2
             Gcol=Gcol+dGcol
             CALL LegendrePolynomial (Lats, Gcol, p2)
             IF (SIGN(1.0_r8,p1) /= SIGN(1.0_r8,p2)) EXIT
          END DO
          IF (dGcol <= eps) EXIT
          Gcol=Gcol-dGcol
          dGcol=dGcol*0.25_r8
       END DO
       CoLatitude(j)=Gcol
    END DO

  END SUBROUTINE GaussianLatitudes


  SUBROUTINE LegendrePolynomial (N, Colatitude, Pln)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N

    REAL (KIND=r8), INTENT(IN) :: Colatitude

    REAL (KIND=r8), INTENT(OUT) :: Pln

    INTEGER :: i

    REAL (KIND=r8) :: x, y1, y2, y3, g

    x=COS(Colatitude)
    y1=1.0_r8
    y2=x
    DO i=2,N
       g=x*y2
       y3=g-y1+g-(g-y1)/REAL(i,r8)
       y1=y2
       y2=y3
    END DO
    Pln=y3

  END SUBROUTINE LegendrePolynomial


END MODULE AreaInterpolation
