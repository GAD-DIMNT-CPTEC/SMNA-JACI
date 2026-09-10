!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM PorceClayMaskSiB2Clima

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER :: Idim, Jdim, LRecOut, LRecGad, ios,i,j,k

   LOGICAL :: GrADS

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=59) :: VarName='PorceClayMaskSiB2Clima'

   CHARACTER (LEN=50) :: VarNameG='PorceClayMaskSiB2ClimaG'

   CHARACTER (LEN=51) :: FileBCs='claymsk.form'

   CHARACTER (LEN=12) :: DirBCs='pre/databcs/'

   CHARACTER (LEN=52) :: NameNML='PorceClayMaskSiB2Clima.nml'

   CHARACTER (LEN=528) :: DirMain
   INTEGER :: Layer
   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfclm=40   ! To Read Formatted Climatological PorcClay Mask
   INTEGER :: nfvgm=50   ! To Write Unformatted Climatological PorcClay Mask
   INTEGER :: nfout=60   ! To Write GrADS Climatological PorcClay Mask
   INTEGER :: nfctl=70   ! To Write GrADS Control File

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: PorcClayMask

   REAL (KIND=r4), DIMENSION (:,:,:), ALLOCATABLE :: PorcClayMask2
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: ClayMaskGad
   NAMELIST /InputDim/ Idim, Jdim,Layer, GrADS, DirMain

   Idim=720
   Jdim=360
   Layer=6
   GrADS=.TRUE.
   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//TRIM(NameNML), &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//TRIM(NameNML), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfinp, NML=InputDim)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')  ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Idim = ', Idim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Jdim = ', Jdim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '     Layer= ', Layer
   WRITE (UNIT=nfprt, FMT='(A,L6)') '    GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')    '  DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'

   ALLOCATE (PorcClayMask(Idim,Jdim), PorcClayMask2(Idim,Jdim,Layer),ClayMaskGad(Idim,Jdim))
   INQUIRE (IOLENGTH=LRecOut) PorcClayMask2
   OPEN (UNIT=nfclm,FILE=TRIM(DirMain)//DirBCs//TRIM(FileBCs),&
        form='UNFORMATTED',ACCESS='DIRECT',recl=LRecOut,ACTION='READ',&
        STATUS='OLD',IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirBCs//TRIM(FileBCs), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfclm, rec=1) PorcClayMask2
   WHERE( PorcClayMask2 >100  )
      PorcClayMask2=0.0
   END WHERE
   DO k=1,Layer
      DO j=1,Jdim
         DO i=1,Idim
           ClayMaskGad(i,j)= ClayMaskGad(i,j) + PorcClayMask2(i,j,k)
         END DO
      END DO
   END DO
   PorcClayMask=INT(ClayMaskGad/REAL(Layer,kind=r4))
   CLOSE (UNIT=nfclm)
   CALL FlipMatrix ()

   INQUIRE (IOLENGTH=LRecOut) PorcClayMask
   OPEN (UNIT=nfvgm, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfvgm, REC=1) PorcClayMask
   CLOSE (UNIT=nfvgm)

   IF (GrADS) THEN
      INQUIRE (IOLENGTH=LRecGad) ClayMaskGad
      OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//'.dat', &
            FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecGad, ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//'.dat', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      ClayMaskGad=REAL(PorcClayMask,r4)
      WRITE (UNIT=nfout, REC=1) ClayMaskGad
      CLOSE (UNIT=nfout)

      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//'.dat'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'UNDEF -999.0'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE SiB2 PorcClay Mask'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Idim,' LINEAR ',0.0_r4,360.0_r4/REAL(Idim,r4)
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'YDEF ',Jdim,' LINEAR ',-89.5_r4,179.0_r4/REAL(Jdim-1,r4)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 0 '
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'VEGM 1 99 PorcClay Mask [No Dim]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)

   END IF
   DEALLOCATE (PorcClayMask, PorcClayMask2,ClayMaskGad)

PRINT *, "*** PorcClayMaskSSiB ENDS NORMALLY ***"

CONTAINS


SUBROUTINE FlipMatrix ()

   ! Flips Over The Rows of a Matrix, After Flips Over
   ! I.D.L. and Greenwitch

   ! Input:
   ! PorceClayMaskSiB2Clima(Idim,Jdim) - Matrix to be Flipped

   ! Output:
   ! PorceClayMaskSiB2Clima(Idim,Jdim) - Flipped Matrix

   INTEGER :: Idimd, Idimd1

   INTEGER :: wk(Idim,Jdim)

   Idimd=Idim/2
   Idimd1=Idimd+1

   wk=PorcClayMask
   PorcClayMask(1:Idimd,:)=wk(Idimd1:Idim,:)
   PorcClayMask(Idimd1:Idim,:)=wk(1:Idimd,:)

   wk=PorcClayMask
   !PorcClayMask(:,1:Jdim,:)=wk(:,Jdim:1:-1,:)

END SUBROUTINE FlipMatrix


END PROGRAM PorceClayMaskSiB2Clima
