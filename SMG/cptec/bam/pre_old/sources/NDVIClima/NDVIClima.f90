!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM NDVIClima

   ! First Point of Input and Output Data is at North Pole and Greenwhich
   ! Over Sea Value is 0.001 cm

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER :: Idim, Jdim, Idim1, LRec, ios

   LOGICAL :: GrADS

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=20) :: VarName='NDVIClima'

   CHARACTER (LEN=41) :: FileBCs='ndviclm.form'

   CHARACTER (LEN=528) :: DirBCs

   CHARACTER (LEN=528) :: DirMain

   REAL (KIND=r4), DIMENSION (:,:,:), ALLOCATABLE :: NDVIC

   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfclm=10   ! To Read Formatted Climatological NDVIC Data
   INTEGER :: nfout=20   ! To Write Unformatted Climatological NDVIC Data
   INTEGER :: nfctl=30   ! To Write Output Data Description

   NAMELIST /InputDim/ Idim, Jdim, GrADS, DirBCs, DirMain

   Idim=360
   Jdim=180
   GrADS=.TRUE.
   DirBCs='./ '
   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//TRIM(VarName)//'.nml', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//TRIM(VarName)//'.nml', &
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

   Idim1=Idim!-1
   ALLOCATE (NDVIC(Idim,Jdim,12))
   INQUIRE (IOLENGTH=LRec) NDVIC(:,:,:)
   OPEN(nfclm,FILE=TRIM(DirBCs)//TRIM(FileBCs),FORM='unformatted',&
       ACCESS='DIRECT',recl=LRec,ACTION='READ',STATUS='OLD',&
       IOSTAT=ios)

   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirBCs)//TRIM(FileBCs), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfclm,rec=1) NDVIC
   CLOSE (UNIT=nfclm)

   INQUIRE (IOLENGTH=LRec) NDVIC(:,:,:)
   OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfout, REC=1) NDVIC(1:Idim1,1:Jdim,1:12)
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'UNDEF -999.0'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE CLimatological NDVI'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Idim1,' LINEAR ',0.0_r4,360.0_r4/REAL(Idim1,r4)
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'YDEF ',Jdim,' LINEAR ',-90.0_r4,180.0_r4/REAL(Jdim,r4)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'ndvi 0 99 NDVI [%]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** NDVIClima ENDS NORMALLY ***"

END PROGRAM NDVIClima
