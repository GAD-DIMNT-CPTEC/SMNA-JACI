!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM TopoSpectral

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nftpi, nftpo, nfout, nfctl, &
                               Imax, Jmax, Mnwv2, &
                               Trunc, nLats, VarNameT, &
                               DirMain, DirPreOut,DirModelIn, &
                               GrADS, &
                               InitInputParameters

   USE FastFourierTransform, ONLY : CreateFFT

   USE LegendreTransform, ONLY : CreateGaussRep, CreateSpectralRep, CreateLegTrans, gLats

   USE SpectralGrid, ONLY : Grid2SpecCoef, SpecCoef2Grid

   IMPLICIT NONE

   INTEGER :: LRec, ios

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: TopogIn, TopogOut

   REAL (KIND=r4), DIMENSION (:), ALLOCATABLE :: CoefTopOut

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: CoefTop

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Topog

   CALL InitInputParameters ()
   CALL CreateSpectralRep ()
   CALL CreateGaussRep ()
   CALL CreateFFT ()
   CALL CreateLegTrans ()

   ALLOCATE (Topog(Imax,Jmax), TopogIn(Imax,Jmax))
   ALLOCATE (CoefTop(Mnwv2), CoefTopOut(Mnwv2))
   IF (GrADS) ALLOCATE (TopogOut(Imax,Jmax))

   ! Read in Gaussian Grid Topography
   INQUIRE (IOLENGTH=LRec) TopogIn
   OPEN (FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//nLats, &
         UNIT=nftpi, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nftpi, REC=1) TopogIn
   CLOSE (UNIT=nftpi)

   Topog=REAL(TopogIn,r8)
   ! SpectraL Coefficient of Topography
   CALL Grid2SpecCoef (1, Topog, CoefTop)

   ! Write Out Topography Spectral Coefficient
   OPEN (FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//Trunc, &
         UNIT=nftpo, FORM='UNFORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//Trunc, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   CoefTopOut=REAL(CoefTop,r4)
   WRITE (UNIT=nftpo) CoefTopOut
   CLOSE (UNIT=nftpo)

   ! Write In Gaussian Grid Topography
   OPEN (FILE=TRIM(DirMain)//TRIM(DirModelIn)//TRIM(VarNameT)//'Rec'//nLats, &
            UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
       WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//'Rec'//nLats, &
             ' returned IOStat = ', ios
       STOP  ' ** (Error) **'
   END IF
   WRITE (UNIT=nfout, REC=1) TopogIn
   CLOSE(nfout,STATUS='KEEP')
   IF (GrADS) THEN
      ! SpectraL Coefficient of Topography
      CALL SpecCoef2Grid (1, CoefTop, Topog)
      TopogOut=REAL(Topog,r4)

      ! Write Out Recomposed Gaussian Grid Topography
      OPEN (FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//'Rec'//nLats, &
            UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//'Rec'//nLats, &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfout, REC=1) TopogIn
      WRITE (UNIT=nfout, REC=2) TopogOut
      CLOSE (UNIT=nfout)

      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//'Rec'//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//'Rec'//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
            TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//'Rec'//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Topography on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      WRITE (UNIT=nfctl, FMT='(8F10.5)') gLats(Jmax:1:-1)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 2'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPI 0 99 Interpolated Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPR 0 99 Recomposed Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)

   END IF
PRINT *, "*** TopoSpectral ENDS NORMALLY ***"

END PROGRAM TopoSpectral
