!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM SnowClima

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nflsm, nfclm, nfout, nfctl, &
                               Imax, Jmax, MonthBefore, MonthAfter, FactorA, FactorB, &
                               DirPreOut, DirModelIn, NameLSM, NameAlb, &
                               VarNameS, nLats, mskfmt, TimeGrADS, &
                               Date, Exts, DirMain, VarName, &
                               IcePoints, GrADS, &
                               InitInputParameters
 
   IMPLICIT NONE

   INTEGER :: i, j, LRecIn, LRecOut,ios

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: gLats

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SnowModel

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: AlbedoBefore

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: AlbedoAfter

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Albedo

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SnowOut

   CALL InitInputParameters ()

   WRITE (nLats(3:7), '(I5.5)') Jmax

   WRITE (mskfmt(2:7), '(I6)') Imax

   ALLOCATE (LandSeaMask(Imax,Jmax))
   ALLOCATE (SnowModel(Imax,Jmax))
   ALLOCATE (AlbedoBefore(Imax,Jmax), AlbedoAfter(Imax,Jmax))
   ALLOCATE (SnowOut(Imax,Jmax), Albedo(Imax,Jmax))

   ! Land Sea Mask : Input
   OPEN (UNIT=nflsm, FILE=TRIM(DirMain)//DirPreOut//NameLSM//nLats, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//NameLSM//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nflsm, FMT=mskfmt) LandSeaMask
   CLOSE (UNIT=nflsm)

   ! Climatological Albedo at Model Grid : Input
   INQUIRE (IOLENGTH=LRecIn) AlbedoBefore
   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//DirPreOut//NameAlb//nLats, &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//NameAlb//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfclm, REC=MonthBefore) AlbedoBefore
   READ  (UNIT=nfclm, REC=MonthAfter) AlbedoAfter
   CLOSE (UNIT=nfclm)

   ! Computes Climatological Snow Based on Climatological Albedo
   DO j=1,Jmax
     DO i=1,Imax
       SnowOut(i,j)=0.0_r8
       ! Linear Interpolation of Albedo in Time
       Albedo(i,j)=FactorA*REAL(AlbedoBefore(i,j),r8)+ &
                   FactorB*REAL(AlbedoAfter(i,j),r8)
       IF (LandSeaMask(i,j) /= 0) THEN
         IF (.NOT.IcePoints) THEN
           IF (Albedo(i,j) >= 0.40_r8 .AND. Albedo(i,j) < 0.49_r8) THEN
             SnowOut(i,j)=5.0_r8
           END IF
           IF (Albedo(i,j) >= 0.49_r8 .AND. Albedo(i,j) < 0.69_r8) THEN
             SnowOut(i,j)=10.0_r8
           END IF
           IF (Albedo(i,j) >= 0.69_r8 .AND. Albedo(i,j) < 0.75_r8) THEN
             SnowOut(i,j)=20.0_r8
           END IF
         END IF
         IF (Albedo(i,j) >= 0.75_r8) THEN
           SnowOut(i,j)=3000.0_r8
         END IF
       END IF
     END DO
   END DO

   ! IEEE-32 Bits Output
   SnowModel=REAL(SnowOut,r4)

   ! Climatological Snow at Model Grid and
   ! Time of Initial Condititon : Model Input
   INQUIRE (IOLENGTH=LRecOut) SnowModel
   OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirModelIn//VarNameS//Date//Exts//nLats, &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//VarNameS//Date//Exts//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfout, REC=1) SnowModel
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      ! Climatological Snow at Model Grid and
      ! Time of Initial Condititon : GrADS
      ! Getting Gaussian Latitudes
      CALL GaussianLatitudes ()
      ! Write GrADS Control File
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//VarName//Date//Exts//nLats//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarName//Date//Exts//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
            TRIM(DirMain)//DirModelIn//VarNameS//Date//Exts//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Climatological Snow on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      WRITE (UNIT=nfctl, FMT='(8F10.5)') gLats(Jmax:1:-1)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR '//TimeGrADS//' 6HR'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'SNOW 0 99 Climatological Snow Depth [kg/m2]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** SnowClima ENDS NORMALLY ***"

CONTAINS


SUBROUTINE GaussianLatitudes ()

   IMPLICIT NONE

   INTEGER :: JmaxH, j

   REAL (KIND=r8) :: eps, rad, dGcolIn, Gcol, dGcol, p1, p2

   ALLOCATE (gLats(Jmax))

   eps=EPSILON(1.0_r8)*100.0_r8
   JmaxH=Jmax/2
   dGcolIn=ATAN(1.0_r8)/REAL(Jmax,r8)
   rad=45.0_r8/ATAN(1.0_r8)
   Gcol=0.0_r8
   DO j=1,JmaxH
      dGcol=dGcolIn
      DO
         CALL LegendrePolynomial (Jmax, Gcol, p2)
         DO
            p1=p2
            Gcol=Gcol+dGcol
            CALL LegendrePolynomial (Jmax, Gcol, p2)
            IF (SIGN(1.0_r8,p1) /= SIGN(1.0_r8,p2)) EXIT
         END DO
         IF (dGcol <= eps) EXIT
         Gcol=Gcol-dGcol
         dGcol=dGcol*0.25_r8
      END DO
      gLats(j)=90.0_r8-rad*Gcol
      gLats(Jmax-j+1)=-gLats(j)
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


END PROGRAM SnowClima
