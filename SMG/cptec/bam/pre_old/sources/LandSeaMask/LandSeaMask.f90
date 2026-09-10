!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM LandSeaMask

   USE InputParameters, ONLY: r4, r8, Undef, &
                              nferr, nfclm, nfoua, nfoub, nfctl, &
                              Imax, Jmax, Idim, Jdim, &
                              nLats, mskfmt, DirPreOut, VarName, &
                              VarNameIn, VarNameOut, DirMain, GrADS, &
                              InitInputParameters

   USE AreaInterpolation, ONLY: gLats, InitAreaInterpolation, &
                                        DoAreaInterpolation

   IMPLICIT NONE

   ! Horizontal Area Interpolator
   ! Interpolate Regular To Gaussian
   ! Regular Input Data is Assumed to be Oriented with
   ! the North Pole and Greenwich as the First Point
   ! Gaussian Output Data is Interpolated to be Oriented with
   ! the North Pole and Greenwich as the First Point
   ! Input for the AGCM is Assumed IEEE 32 Bits Big Endian

   INTEGER :: m, j, LRecIn, LRecOut, ios

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMaskOut

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: WaterIn
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: WaterOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: WaterInput
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: WaterOutput

   CALL InitInputParameters ()

   CALL InitAreaInterpolation ()

   ALLOCATE (WaterInput(Idim,Jdim))
   ALLOCATE (WaterOutput(Imax,Jmax))
   ALLOCATE (WaterIn(Idim,Jdim))
   ALLOCATE (WaterOut(Imax,Jmax))
   ALLOCATE (LandSeaMaskOut(Imax,Jmax))

   ! Read In Input Water

   INQUIRE (IOLENGTH=LRecIn) WaterIn(:,1)
   OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameIn)//'.dat', &
         UNIT=nfclm, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarNameIn)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO j=1,Jdim
      READ (UNIT=nfclm, REC=j) WaterIn(:,j)
   END DO
   CLOSE (UNIT=nfclm)
   WaterInput=REAL(WaterIn,r8)
  
   ! Interpolate Input Water To Output Water
  
   CALL DoAreaInterpolation (WaterInput, WaterOutput)
   WaterOut=REAL(WaterOutput,r4)

   ! Generate Land Sea Mask (= 1 over Land and = 0 over Sea)
   WHERE (WaterOutput < 50.0_r8)
        LandSeaMaskOut=1
   ELSEWHERE
        LandSeaMaskOut=0
   ENDWHERE

   OPEN (FILE=TRIM(DirMain)//DirPreOut//VarName//nLats, &
         UNIT=nfoua, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfoua, FMT=mskfmt) LandSeaMaskOut
   CLOSE (UNIT=nfoua)
  
   IF (GrADS) THEN
      ! Write Out Adjusted Interpolated Water
      INQUIRE (IOLENGTH=LRecOut) WaterOut
      OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//nLats//'.dat', &
            UNIT=nfoub, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//nLats//'.dat', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfoub, REC=1) WaterOut
      WRITE (UNIT=nfoub, REC=2) REAL(LandSeaMaskOut,r4)
      CLOSE (UNIT=nfoub)

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
             TRIM(DirMain)//DirPreOut//TRIM(VarNameOut)//nLats//'.dat'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Land Sea Mask on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0/REAL(Imax,8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      WRITE (UNIT=nfctl, FMT='(8F10.5)') gLats(Jmax:1:-1)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 2'
      WRITE (UNIT=nfctl, FMT='(A)') 'WPER  0 99 Percentage of Water [%]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask [1-Land 0-Sea]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** LandSeaMask ENDS NORMALLY ***"

END PROGRAM LandSeaMask
