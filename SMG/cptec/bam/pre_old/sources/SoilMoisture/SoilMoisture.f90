!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM SoilMoisture

   USE InputParameters, ONLY: r4, r8, Undef, SeaValue, &
                              nferr, nflsm, nfclm, nfout, nfctl, &
                              Imax, Jmax, Idim, Jdim, &
                              nLats, mskfmt, NameLSMask, DirPreOut, DirModelIn, &
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

   INTEGER :: Month, LRecIn, LRecOut, ios

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureIn
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureInput
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: SoilMoistureOutput

   CALL InitInputParameters ()
   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (SoilMoistureInput(Idim,Jdim))
   ALLOCATE (SoilMoistureOutput(Imax,Jmax))
   ALLOCATE (SoilMoistureIn(Idim,Jdim))
   ALLOCATE (SoilMoistureOut(Imax,Jmax))
   ALLOCATE (LandSeaMask(Imax,Jmax))

   ! Land Sea Mask : Input
   OPEN (FILE=TRIM(DirMain)//DirPreOut//NameLSMask//nLats, &
         UNIT=nflsm, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//NameLSMask//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nflsm, FMT=mskfmt) LandSeaMask
   CLOSE (UNIT=nflsm)

   INQUIRE (IOLENGTH=LRecIn) SoilMoistureIn
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

   INQUIRE (IOLENGTH=LRecOut) SoilMoistureOut
   OPEN (FILE=TRIM(DirMain)//DirModelIn//VarNameOut//nLats, &
         UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//VarNameOut//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   ! Read In Input SoilMoisture

   DO Month=1,12

      READ (UNIT=nfclm, REC=Month) SoilMoistureIn
      SoilMoistureInput=REAL(SoilMoistureIn,r8)
  
      ! Interpolate Input Regular SoilMoisture To Gaussian Grid Output
  
      IF (Linear) THEN
         CALL DoLinearInterpolation (SoilMoistureInput, SoilMoistureOutput)
      ELSE
         CALL DoAreaInterpolation (SoilMoistureInput, SoilMoistureOutput)
      END IF

      ! Adjust Value Over Sea (150)
      WHERE (LandSeaMask == 0) SoilMoistureOutput=SeaValue
      SoilMoistureOut=REAL(SoilMoistureOutput,r4)
  
      ! Write Out Adjusted Interpolated SoilMoisture
  
      WRITE (UNIT=nfout, REC=Month) SoilMoistureOut
  
   END DO

   CLOSE (UNIT=nfclm)
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Soil Moisture on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF  1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  1'
      WRITE (UNIT=nfctl, FMT='(A)') 'SOMO  0 99 SoilMoisture [cm]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** SoilMoisture ENDS NORMALLY ***"

END PROGRAM SoilMoisture
