!
!  $Author: bonatti $
!  $Date: 2007/09/18 18:07:15 $
!  $Revision: 1.2 $
!
PROGRAM SoilMoistureClima

   ! First Point of Initial Data is at North Pole and I. D. Line
   ! First Point of Output  Data is at North Pole and Greenwhich

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER :: Idim, Jdim, Month, LRec, ios

   LOGICAL :: GrADS

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=17) :: VarName='SoilMoistureClima'

   CHARACTER (LEN=11) :: FileBCs='soilms.form'

   CHARACTER (LEN=528) :: DirBCs

   CHARACTER (LEN=528) :: DirMain

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SoilMoisture

   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfclm=10   ! To Read Formatted Climatological Soil Moisture Data
   INTEGER :: nfout=20   ! To Write Unformatted Climatological Soil Moisture Data
   INTEGER :: nfctl=30   ! To Write Output Data Description

   NAMELIST /InputDim/ Idim, Jdim, GrADS, DirBCs, DirMain

   Idim=360
   Jdim=181
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

   ALLOCATE (SoilMoisture(Idim,Jdim))

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

   INQUIRE (IOLENGTH=LRec) SoilMoisture
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

   DO Month=1,12
      READ (UNIT=nfclm, FMT='(5E15.8)') SoilMoisture
      CALL FlipMatrix (SoilMoisture, Idim, Jdim)
      WRITE (UNIT=nfout, REC=Month) SoilMoisture
   END DO

   CLOSE (UNIT=nfclm)
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE CLimatological Soil Moisture'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Idim,' LINEAR ',0.0_r4,360.0_r4/REAL(Idim,r4)
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'YDEF ',Jdim,' LINEAR ',-90.0_r4,180.0_r4/REAL(Jdim-1,r4)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF  1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  1'
      WRITE (UNIT=nfctl, FMT='(A)') 'SOMO  0 99 SoilMoisture [cm]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** SoilMoistureClima ENDS NORMALLY ***"

CONTAINS


SUBROUTINE FlipMatrix (h, Idim, Jdim)

   ! Flips a Matrix Over I.D.L. and Greenwitch

   ! Input:
   ! h(Idim,Jdim) - Matrix to be Flipped
   !   Idim       - Column Dimension of h(Idim,Jdim)
   !   Jdim       - Row Dimension of h(Idim,Jdim)

   ! Output:
   ! h(Idim,Jdim) - Flipped Matrix


   REAL (KIND=r4), INTENT(INOUT) :: h(Idim,Jdim)

   INTEGER, INTENT(IN) :: Idim

   INTEGER, INTENT(IN) :: Jdim

   INTEGER :: Idimd, Idimd1

   REAL (KIND=r4) :: wk(Idim,Jdim)

   Idimd=Idim/2
   Idimd1=Idimd+1

   wk=h
   h(1:Idimd,:)=wk(Idimd1:Idim,:)
   h(Idimd1:Idim,:)=wk(1:Idimd,:)

END SUBROUTINE FlipMatrix


END PROGRAM SoilMoistureClima
