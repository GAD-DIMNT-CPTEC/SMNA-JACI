!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM DeepSoilTemperatureClima

   ! First Point of Input and Output Data is at North Pole and Greenwhich

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER :: Idim, Jdim, LRec, ios

   LOGICAL :: GrADS

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=24) :: VarName='DeepSoilTemperatureClima'

   CHARACTER (LEN=11) :: FileBCs='tgdeep.form'

   CHARACTER (LEN=528) :: DirBCs

   CHARACTER (LEN=528) :: DirMain

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: DeepSoilTemp

   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfclm=10   ! To Read Formatted Climatological DeepSoilTemp Data
   INTEGER :: nfout=20   ! To Write Unformatted Climatological DeepSoilTemp Data
   INTEGER :: nfctl=30   ! To Write Output Data Description

   NAMELIST /InputDim/ Idim, Jdim, GrADS, DirBCs, DirMain

   Idim=138
   Jdim=116
   GrADS=.TRUE.
   DirBCs='./ '
   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//VarName//'.nml', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//VarName//'.nml', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfinp, NML=InputDim)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')  ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A,I6)') '    Idim = ', Idim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '    Jdim = ', Jdim
   WRITE (UNIT=nfprt, FMT='(A,L6)') '   GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')    '  DirBCs = '//TRIM(DirBCs)
   WRITE (UNIT=nfprt, FMT='(A)')    ' DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'

   ALLOCATE (DeepSoilTemp(Idim,Jdim))

   OPEN (UNIT=nfclm, FILE=TRIM(DirBCs)//FileBCs, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirBCs)//FileBCs, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfclm, FMT='(5E15.8)') DeepSoilTemp
   CLOSE (UNIT=nfclm)

   INQUIRE (IOLENGTH=LRec) DeepSoilTemp
   OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//VarName//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarName//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfout, REC=1) DeepSoilTemp
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//VarName//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarName//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//VarName//'.dat'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'UNDEF -999.0'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE CLimatological Deep Soil Temperature'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Idim,' LINEAR ',0.0_r4,360.0_r4/REAL(Idim,r4)
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'YDEF ',Jdim,' LINEAR ',-90.0_r4,180.0_r4/REAL(Jdim-1,r4)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'DSTP 0 99 Deep Soil Temperature [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** DeepSoilTemperatureClima ENDS NORMALLY ***"


END PROGRAM DeepSoilTemperatureClima
