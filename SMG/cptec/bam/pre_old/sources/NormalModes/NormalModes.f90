!
!  $Author: bonatti $
!  $Date: 2009/03/20 18:00:00 $
!  $Revision: 1.1.1.1 $
!
PROGRAM NormalModes

   USE InputParameters, ONLY: r8, nferr, nfprt, nfmod, nfsig, &
                              Kmax, go, HnCut, &
                              FileNormalModes, FileDeltaSigma, &
                              DirDat, DirSig, &
                              InitInputParameters

   USE VerticalModes, ONLY: GetVerticalModes

   USE HorizontalModes, ONLY: GetHorizontalModes

   IMPLICIT NONE

   ! Compute Vertical and Triangular Truncation Horizontal Modes

   INTEGER :: Mods, ios

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: gh, DeltaSigma

   CALL InitInputParameters ()

   ! File To Write Selected Normal Modes (Hn > HnCut ; Period > PerCut)
   OPEN (UNIT=nfmod, FILE=TRIM(DirDat)//FileNormalModes, &
         FORM='UNFORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
            TRIM(DirDat)//FileNormalModes, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   ! Vertical Modes
   ALLOCATE (gh(Kmax), DeltaSigma(Kmax))
   !OPEN  (UNIT=nfsig, FILE=TRIM(DirSig)//FileDeltaSigma, STATUS='OLD')
   ! File To Read Delta Sigma for Kmax
   OPEN (UNIT=nfsig, FILE=TRIM(DirSig)//FileDeltaSigma, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
            TRIM(DirDat)//FileDeltaSigma, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfsig, FMT='(5F9.6)') DeltaSigma
   CLOSE (UNIT=nfsig)
   WRITE (UNIT=nfprt, FMT='(A)') ' DeltaSigma :'
   WRITE (UNIT=nfprt, FMT='(5F9.6)') DeltaSigma
   WRITE (UNIT=nfprt, FMT='(A)') ' '
   CALL GetVerticalModes (gh, DeltaSigma)

   ! Get Number of Modes With Equivalent Depth > HnCut
   CALL SetNumberOfModes ()

   ! Horizontal Modes
   CALL GetHorizontalModes (Mods, gh)

   CLOSE (UNIT=nfmod)


CONTAINS


SUBROUTINE SetNumberOfModes ()

   IMPLICIT NONE

   INTEGER :: n

   ! Get Number of Modes With Equivalent Depth > HnCut
   WRITE (UNIT=nfprt, FMT='(A)') ' '
   Mods=0
   DO n=1,Kmax
      IF (gh(n)/go > HnCut) THEN
         Mods=Mods+1
         WRITE (UNIT=nfprt, FMT='(A,I5,2(A,F10.2))') &
               ' n = ', n, ' HnCut = ', HnCut, ' Hn = ', gh(n)/go
      END IF
   END DO
   IF (Mods <= 0 .OR. Mods > Kmax) THEN
      WRITE (UNIT=nferr, FMT='(A)') &
            ' ERROR: The Equivalent Heights of Normal Modes is Wrong'
      STOP ' (Mods <= 0 .OR. Mods > Kmax)'
   END IF
   WRITE (UNIT=nfprt, FMT='(A,I5,2(A,F10.2))') &
         ' n = ', Mods+1, ' HnCut = ', HnCut, ' Hn = ', gh(Mods+1)/go
   WRITE (UNIT=nfprt, FMT='(/,A,I5,/)') ' Mods = ',Mods

END SUBROUTINE SetNumberOfModes


END PROGRAM NormalModes
