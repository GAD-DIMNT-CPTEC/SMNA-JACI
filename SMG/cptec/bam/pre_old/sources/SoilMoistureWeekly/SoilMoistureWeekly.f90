!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
! 
PROGRAM SoilMoistureWeekly

   USE InputParameters, ONLY: r4, r8, &
                              nferr, nfclm, nfout, nfctl, &
                              Imax, Jmax, Idim, Jdim,kdim, &
                              nLats, DirPreOut, DirModelIn, &
                              VarName, VarNameOut, DirMain, GrADS, Linear, &
                              InitInputParameters,Date

   USE LinearInterpolation, ONLY: gLatsL=>LatOut, &
       InitLinearInterpolation, DoLinearInterpolation

   USE AreaInterpolation, ONLY: gLatsA=>gLats, &
       InitAreaInterpolation, DoAreaInterpolation

   IMPLICIT NONE
   ! Linear (.TRUE.) or Area Weighted (.FALSE.) Interpolation
   ! Interpolate Regular To Gaussian
   ! Regular Input Data is Assumed to be Oriented with
   ! the North Pole and Greenwich as the First Point
   ! Gaussian Output Data is Interpolated to be Oriented with
   ! the North Pole and Greenwich as the First Point
   ! Input for the AGCM is Assumed IEEE 32 Bits Big Endian

   INTEGER :: LRecIn, LRecOut, ios,Month,k,i,j
   CHARACTER (LEN=12) :: TimeGrADS='  Z         '
   CHARACTER (LEN=3), DIMENSION(12) :: MonthChar = &
             (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC' /)

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureWeeklyIn
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureWeeklyInput
   REAL (KIND=r4), DIMENSION (:,:,:), ALLOCATABLE :: SoilMoistureWeeklyOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureWeeklyOutput

   CALL InitInputParameters ()
   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (SoilMoistureWeeklyIn    (Idim,Jdim))
   ALLOCATE (SoilMoistureWeeklyInput (Idim,Jdim))
   ALLOCATE (SoilMoistureWeeklyOut   (Imax,Jmax,kdim))
   ALLOCATE (SoilMoistureWeeklyOutput(Imax,Jmax))

   ! Read In Input SoilMoistureWeekly

   INQUIRE (IOLENGTH=LRecIn) SoilMoistureWeeklyIn (1:Idim,1:Jdim)
   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.'//TRIM(Date(1:8)), &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', &
         STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.'//TRIM(Date(1:8)), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO k=1,kdim
      READ  (UNIT=nfclm, REC=k) SoilMoistureWeeklyIn(1:Idim,1:Jdim)
      SoilMoistureWeeklyInput(1:Idim,1:Jdim)=REAL(SoilMoistureWeeklyIn(1:Idim,1:Jdim),r8)

      ! Interpolate Input Regular Grid Deep Soil Temperature To Gaussian Grid on Output

      IF (Linear) THEN
         CALL DoLinearInterpolation (SoilMoistureWeeklyInput(:,:), SoilMoistureWeeklyOutput(:,:))
      ELSE
         CALL DoAreaInterpolation (SoilMoistureWeeklyInput(:,:), SoilMoistureWeeklyOutput(:,:))
      END IF
      SoilMoistureWeeklyOut(1:Imax,1:Jmax,k)=REAL(SoilMoistureWeeklyOutput(1:Imax,1:Jmax),r4)
      DO j=1,jMax
         DO i=1,iMax
            IF (SoilMoistureWeeklyOut(i,j,k) <=0.0_r4)THEN
                SoilMoistureWeeklyOut(i,j,k) = 0.9_r4
            END IF
         END DO
      END DO
   END DO   
   
   CLOSE (UNIT=nfclm)

   ! Generates ZoRL Input File for the AGCM Model

   INQUIRE (IOLENGTH=LRecOut) SoilMoistureWeeklyOut(1:Imax,1:Jmax,1)
   OPEN (FILE=TRIM(DirMain)//DirModelIn//TRIM(VarNameOut)//'.'//TRIM(Date(1:8))//nLats, &
         UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//TRIM(VarNameOut)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO k=1,kdim
      WRITE (UNIT=nfout, REC=k) SoilMoistureWeeklyOut(1:Imax,1:Jmax,k)
   END DO
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
   TimeGrADS(1:2)=Date(9:10)
   TimeGrADS(4:5)=Date(7:8)
   TimeGrADS(9:12)=Date(1:4)
   READ (Date(5:6), FMT='(I2)') Month
   TimeGrADS(6:8)=MonthChar(Month)

      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//'.'//TRIM(Date(1:8))//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//'.'//TRIM(Date(1:8))//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirModelIn//TRIM(VarNameOut)//'.'//TRIM(Date(1:8))//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'UNDEF -999.0'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE SoilMoistureWeekly on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 8 LEVELS 1 2 3 4 5 6 7 8'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR '//TimeGrADS//' 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'SoilMoisture 8 99 SoilMoistureWeekly [kg/m3]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** SoilMoistureWeekly ENDS NORMALLY ***"

END PROGRAM SoilMoistureWeekly
