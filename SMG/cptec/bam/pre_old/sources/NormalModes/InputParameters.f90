!
!  $Author: bonatti $
!  $Date: 2009/03/20 18:00:00 $
!  $Revision: 1.1.1.1 $
!
MODULE InputParameters

   IMPLICIT NONE

   PRIVATE

   ! Selecting Kinds
   ! Kind for 64-bits Integer Numbers
   INTEGER, PARAMETER, PUBLIC :: i8=SELECTED_INT_KIND(14)
   ! Kind for 64-bits Real Numbers
   INTEGER, PARAMETER, PUBLIC :: r8=SELECTED_REAL_KIND(15)

   INTEGER, PUBLIC :: Mend, Kmax, Mend1, KmaxM, KmaxP, &
                      NxSym, NxAsy, MakeVec

   REAL (KIND=r8), PUBLIC :: Eps, HnCut, PerCut, twoOmega, SqRt2

   CHARACTER (LEN=256), PUBLIC :: DirMain, DirDat , DirSig

   LOGICAL, PUBLIC :: PrtOut

   REAL (KIND=r8), PUBLIC :: &
         go=9.80665_r8,   & ! Gravity Acceleration (m/s2)
         Rd=287.05_r8,    & ! Dry Air Gas Constant (m2/s2/K)
         Cp=1004.6_r8,    & ! Dry Air Heat Capacity at Constant Pressure (m2/s2/K)
         RE=6370000.0_r8, & ! Mean Earth Radius (m)
         Ps=1013.0_r8       ! Pessure Value of Reference (hPa)

   CHARACTER (LEN=13), PUBLIC :: FileNormalModes='NMI.T    L   '
   CHARACTER (LEN=15), PUBLIC :: FileDeltaSigma='DeltaSigma.L   '

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfinp=5    ! Standard Read In
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nfsig=10   ! To Read Delta Sigma Data
   INTEGER, PUBLIC :: nfmod=20   ! To Write Selected Normal Modes

   PUBLIC :: InitInputParameters


CONTAINS


SUBROUTINE InitInputParameters ()

   IMPLICIT NONE

   REAL (KIND=r8) :: Pi, PerCutDays

   INTEGER ::  Modd, Lend, Kend, ios

   CHARACTER (LEN=15) :: NameNML='NormalModes.nml'

   NAMELIST /NorModNML/ Mend,Kmax, HnCut, PerCutDays, &
                        PrtOut, DirMain

   Mend=62
   Kmax=28
   HnCut=1000.0_r8
   PerCutDays=2.0_r8
   PrtOut=.FALSE.
   DirMain='./ '

   OPEN  (UNIT=nfinp, FILE='./'//NameNML, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
             ' ** (Error) ** Open file ', &
              './'//NameNML, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   READ  (UNIT=nfinp, NML=NorModNML)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')    ' &NorModNML'
   WRITE (UNIT=nfprt, FMT='(A,I4)')   '        Mend = ', Mend
   WRITE (UNIT=nfprt, FMT='(A,I4)')   '        Kmax = ', Kmax
   WRITE (UNIT=nfprt, FMT='(A,F8.2)') '       HnCut = ', HnCut
   WRITE (UNIT=nfprt, FMT='(A,F6.2)') '  PerCutDays = ', PerCutDays
   WRITE (UNIT=nfprt, FMT='(A,L6)')   '      PrtOut = ', PrtOut
   WRITE (UNIT=nfprt, FMT='(2A)')     '     DirMain = ', TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')    ' /'

   Eps=1.0_r8
   !Eps=4.0_r8*EPSILON(Eps)
   Eps=EPSILON(Eps)
   Pi=4.0_r8*ATAN(1.0_r8)

   Mend1=Mend+1
   KmaxM=Kmax-1
   KmaxP=Kmax+1

   Modd=Mend-(Mend/2)*2
   Lend=(Mend+Modd)/2
   Kend=Lend+1-Modd
   NxSym=Lend+2*Kend
   NxAsy=2*Lend+Kend

   MakeVec=1
   PerCut=PerCutDays*86400.0_r8/(2.0_r8*Pi)

   twoOmega=4.0_r8*Pi/86400.0_r8

   SqRt2=1.0_r8/SQRT(2.0_r8)

   DirDat=TRIM(DirMain)//'model/datain/ '
   DirSig=TRIM(DirMain)//'pre/datain/ '

   WRITE (FileNormalModes(6:9),   FMT='(I4.4)') Mend
   WRITE (FileNormalModes(11:13), FMT='(I3.3)') Kmax
   WRITE (FileDeltaSigma(13:15),  FMT='(I3.3)') Kmax

   WRITE (UNIT=nfprt, FMT='(A25,1PG12.5,/)') &
         ' Normal Mode Precision = ', Eps
   WRITE (UNIT=nfprt, FMT='(2A)')   '          DirDat = ', TRIM(DirDat)
   WRITE (UNIT=nfprt, FMT='(2A)')   '          DirSig = ', TRIM(DirSig)
   WRITE (UNIT=nfprt, FMT='(2A)')   '  FileDeltaSigma = ', FileDeltaSigma
   WRITE (UNIT=nfprt, FMT='(2A,/)') ' FileNormalModes = ', FileNormalModes

END SUBROUTINE InitInputParameters


END MODULE InputParameters
