!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM SNOWWeekly

   USE InputParameters, ONLY: r4, r8, &
                              nferr, nfclm, nfout, nfctl, &
                              Imax, Jmax, Idim, Jdim, &
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

   INTEGER :: LRecIn, LRecOut, ios,Month
   CHARACTER (LEN=12) :: TimeGrADS='  Z         '
   CHARACTER (LEN=3), DIMENSION(12) :: MonthChar = &
             (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC' /)

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SNOWWeeklyIn
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SNOWWeeklyInput
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SNOWWeeklyOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SNOWWeeklyOutput

   CALL InitInputParameters ()
   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (SNOWWeeklyIn(Idim,Jdim))
   ALLOCATE (SNOWWeeklyInput(Idim,Jdim))
   ALLOCATE (SNOWWeeklyOut(Imax,Jmax))
   ALLOCATE (SNOWWeeklyOutput(Imax,Jmax))

   ! Read In Input SNOWWeekly

   INQUIRE (IOLENGTH=LRecIn) SNOWWeeklyIn
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
   READ  (UNIT=nfclm, REC=1) SNOWWeeklyIn
   SNOWWeeklyInput=REAL(SNOWWeeklyIn,r8)
   CLOSE (UNIT=nfclm)

   ! Interpolate Input Regular Grid Deep Soil Temperature To Gaussian Grid on Output

   IF (Linear) THEN
      CALL DoLinearInterpolation (SNOWWeeklyInput(:,:), SNOWWeeklyOutput(:,:))
   ELSE
      CALL DoAreaInterpolation (SNOWWeeklyInput(:,:), SNOWWeeklyOutput(:,:))
   END IF
   SNOWWeeklyOut(:,:)=REAL(SNOWWeeklyOutput(:,:),r4)
   
   ! Generates ZoRL Input File for the AGCM Model

   INQUIRE (IOLENGTH=LRecOut) SNOWWeeklyOut
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
   WRITE (UNIT=nfout, REC=1) SNOWWeeklyOut
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE SNOWWeekly on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR '//TimeGrADS//' 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'snow 0 99 SNOWWeekly [kg/m3]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** SNOWWeekly ENDS NORMALLY ***"

END PROGRAM SNOWWeekly
