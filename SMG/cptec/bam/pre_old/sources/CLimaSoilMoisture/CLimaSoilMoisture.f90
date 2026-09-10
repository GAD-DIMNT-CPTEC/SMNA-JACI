!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM SoilMoisture

   USE InputParameters, ONLY: r4, r8, &
                              nferr, nfclm, nfout, nfctl, &
                              Imax, Jmax, Idim,Layer, Jdim, &
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

   INTEGER :: LRecIn, LRecOut, ios,k,j

   REAL (KIND=r4), DIMENSION (:,:,:,:), ALLOCATABLE :: SoilMoistureIn
   REAL (KIND=r8), DIMENSION (:,:,:,:), ALLOCATABLE :: SoilMoistureInput
   REAL (KIND=r4), DIMENSION (:,:,:,:), ALLOCATABLE :: SoilMoistureOut
   REAL (KIND=r8), DIMENSION (:,:,:,:), ALLOCATABLE :: SoilMoistureOutput

   CALL InitInputParameters ()
   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (SoilMoistureIn(Idim,Jdim,Layer,12))
   ALLOCATE (SoilMoistureInput(Idim,Jdim,Layer,12))
   ALLOCATE (SoilMoistureOut(Imax,Jmax,Layer,12))
   ALLOCATE (SoilMoistureOutput(Imax,Jmax,Layer,12))

   ! Read In Input SoilMoisture

   INQUIRE (IOLENGTH=LRecIn) SoilMoistureIn
   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', &
         STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfclm, REC=1) SoilMoistureIn
   SoilMoistureInput=REAL(SoilMoistureIn,r8)
   CLOSE (UNIT=nfclm)

   ! Interpolate Input Regular Grid Deep Soil SoilMoisture To Gaussian Grid on Output

   DO k=1,12
      DO j=1,Layer
         IF (Linear) THEN
            CALL DoLinearInterpolation (SoilMoistureInput(:,:,j,k), SoilMoistureOutput(:,:,j,k))
         ELSE
            CALL DoAreaInterpolation (SoilMoistureInput(:,:,j,k), SoilMoistureOutput(:,:,j,k))
         END IF
         SoilMoistureOut(:,:,j,k)=REAL(SoilMoistureOutput(:,:,j,k),r4)
      END DO
   END DO
   ! Generates ZoRL Input File for the AGCM Model

   INQUIRE (IOLENGTH=LRecOut) SoilMoistureOut
   OPEN (FILE=TRIM(DirMain)//DirModelIn//TRIM(VarNameOut)//nLats, &
         UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//TRIM(VarNameOut)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfout, REC=1) SoilMoistureOut
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirModelIn//TRIM(VarNameOut)//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'UNDEF -999.0'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE SoilMoisture on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 3 LEVELS 0 1 2'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'soilms 3 99 SoilMoisture [C]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** SoilMoisture ENDS NORMALLY ***"

END PROGRAM SoilMoisture
