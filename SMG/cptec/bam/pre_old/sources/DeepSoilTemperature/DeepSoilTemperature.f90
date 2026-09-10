!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM DeepSoilTemperature

   USE InputParameters, ONLY: r4, r8, Undef, &
                              nferr, nfclm, nfout, nfctl, &
                              Imax, Jmax, Idim, Jdim, &
                              nLats, DirPreOut, DirModelIn, &
                              VarName, VarNameOut, DirMain, GrADS, Linear, &
                              InitInputParameters

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

   INTEGER :: LRecIn, LRecOut, ios

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: DeepSoilTempIn
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: DeepSoilTempInput
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: DeepSoilTempOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: DeepSoilTempOutput

   CALL InitInputParameters ()
   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (DeepSoilTempIn(Idim,Jdim))
   ALLOCATE (DeepSoilTempInput(Idim,Jdim))
   ALLOCATE (DeepSoilTempOut(Imax,Jmax))
   ALLOCATE (DeepSoilTempOutput(Imax,Jmax))

   ! Read In Input Deep Soil Temperature

   INQUIRE (IOLENGTH=LRecIn) DeepSoilTempIn
   OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
         UNIT=nfclm, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfclm, REC=1) DeepSoilTempIn
   DeepSoilTempInput=REAL(DeepSoilTempIn,r8)
   CLOSE (UNIT=nfclm)

   ! Interpolate Input Regular Grid Deep Soil Temperature To Gaussian Grid on Output

   IF (Linear) THEN
      CALL DoLinearInterpolation (DeepSoilTempInput, DeepSoilTempOutput)
   ELSE
      CALL DoAreaInterpolation (DeepSoilTempInput, DeepSoilTempOutput)
   END IF
   DeepSoilTempOut=REAL(DeepSoilTempOutput,r4)

   ! Generates Tg1, Tg2, Tg3 Input File for the AGCM Model
   ! At this Climatological Input It Is Assumed that Tg1=Tg2=Tg3

   INQUIRE (IOLENGTH=LRecOut) DeepSoilTempOut
   OPEN (FILE=TRIM(DirMain)//DirModelIn//VarNameOut//nLats, &
         UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarNameOut//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfout, Rec=1) DeepSoilTempOut
   WRITE (UNIT=nfout, Rec=2) DeepSoilTempOut
   WRITE (UNIT=nfout, Rec=3) DeepSoilTempOut
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//DirPreOut//VarNameOut//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarNameOut//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirModelIn//VarNameOut//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Soil Temperature on a Gaussian Grid'
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 3'
      WRITE (UNIT=nfctl, FMT='(A)') 'DST1 0 99 Deep Soil Temperature [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'DST2 0 99 Deep Soil Temperature [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'DST3 0 99 Deep Soil Temperature [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** DeepSoilTemperature ENDS NORMALLY ***"

END PROGRAM DeepSoilTemperature
