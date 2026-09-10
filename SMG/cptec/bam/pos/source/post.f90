!
!  $Author: pkubota $
!  $Date: 2007/10/10 20:28:03 $
!  $Revision: 1.3 $
!
PROGRAM GLobalModelPostProcessing

  ! Reads Spectral Forecast Coefficients of Topography,
  ! Log of Sfc Pressure, Temperature, Divergence, Vorticity
  ! and Humidity in Sigma Layers. Converts These Values
  ! to Selected Mandatory Pressure Levels.

  USE Parallelism, ONLY:   &
       CreateParallelism,  &
       DestroyParallelism, &
       MsgOne,             &
       FatalError,         &
       unitDump,           &
       myId_four,          &
       myId

  USE Constants, ONLY : InitParameters, nfprt, nFFrs, nFBeg, nFEnd,datalib,rfd
  USE Conversion, ONLY : CreateConversion
  USE Sizes, ONLY : kmax, lmax, mmax
  USE Init,  ONLY : Initall
  USE InputArrays,  ONLY : GetArrays
  USE PrblSize, ONLY : CreatePrblSize,trunc_in,nlon,nlat,vert_in,vert_out
  USE GaussSigma, ONLY : CreateHybrCoor
  USE RegInterp, ONLY : InitAreaInterpolation
  USE FileAccess, ONLY : InitFiles
  USE PostLoop, ONLY : postgl,InitPostLoop
  USE tables, ONLY: tables_readed,Init_tables
  USE PhysicalFunctions, ONLY: InitPhysicalFunctions
  USE Watches

  IMPLICIT NONE


  INTEGER :: nFile
  
!!MARCELO 1
  LOGICAL, PARAMETER :: instrument=.TRUE.
  INTEGER :: nThreads=0
  INTEGER :: iThread
  TYPE(Watch), ALLOCATABLE :: wt(:)
!!MARCELO 1


  INCLUDE 'mpif.h'

  ! engage MPI

  CALL CreateParallelism()

  CALL InitParameters ()
  CALL InitPhysicalFunctions()
  CALL InitFiles ()
  CALL InitPostLoop()
  CALL CreatePrblSize ()
  CALL InitAll(trunc_in, nlon, nlat, vert_in, vert_out)
  CALL GetArrays
  CALL CreateConversion ()
  CALL CreateHybrCoor   ()
  CALL InitAreaInterpolation ()
  !If grib tables not read proceed
  
  IF(.not. tables_readed) CALL Init_tables(datalib,rfd)

  IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,I3.3,A,I2.2,/)') &
        ' Post-Processing Resolution: T', Mmax, 'L',Kmax

  ! Do Post-Processing for Files nFFrs to nFEnd
  DO nFile=nFFrs,nFEnd
     IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(A,I5)') ' nFile = ',nFile
     CALL postgl (nFFrs, nFBeg, nFEnd, nFile)
     IF(myid.eq.0) WRITE (UNIT=nfprt, FMT='(/,A,/)') ' Advanced Time Loop'
  END DO

!!MARCELO 3
  IF (instrument) THEN
     iThread = 0
     CALL DumpWatch(wt(iThread), unitDump,'TempoPosProcessamento')
     CALL DestroyWatch(wt(iThread))
  END IF
!!MARCELO 3

  CALL DestroyParallelism("*Post-Processing ends normally*")

END PROGRAM GLobalModelPostProcessing
