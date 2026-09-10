MODULE Parallelism

  ! Deals with parallelism; exports
  !   MPI procedures for
  !       initialization 
  !       normal finalization 
  !       error finalization
  !   MPI process data:
  !       how many processes in MPI computation
  !       which processes is this
  ! Creates one file per process for printouts,
  ! unscrambling outputs. Provides procedures for:
  !       writing messages to output and to the file
  !       writing messages only to the file


  IMPLICIT NONE

  PRIVATE

  ! public data

  INTEGER,          PUBLIC :: maxNodes      ! # MPI processes in the computation
  INTEGER,          PUBLIC :: myId          ! MPI process rank
  INTEGER,          PUBLIC :: maxNodes_four ! # MPI processes in fourier group  
  INTEGER,          PUBLIC :: myId_four     ! MPI process rank in fourier group
  INTEGER,          PUBLIC :: mygroup_four  ! fourier group
  INTEGER,          PUBLIC :: unitDump      ! this process dumping file unit
  INTEGER,          PUBLIC :: COMM_FOUR     ! Communicator of Fourier Group

  ! public procedures

  PUBLIC :: CreateParallelism
  PUBLIC :: CreateFourierGroup
  PUBLIC :: FatalError
  PUBLIC :: DestroyParallelism
  PUBLIC :: MsgOut
  PUBLIC :: MsgDump
  PUBLIC :: MsgOne

  ! private data

  CHARACTER(LEN=4) :: cNProc      ! maxNodes in characters, for printing
  CHARACTER(LEN=4) :: cThisProc   ! myId in characters, for printing
  INTEGER, PARAMETER :: stdout=6

  INCLUDE 'mpif.h'

CONTAINS



  !*** Create Parallelism ***



  SUBROUTINE CreateParallelism()
    INTEGER :: ierror
    INTEGER :: CompName, numThreads
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: Pname
    CHARACTER(LEN=*), PARAMETER :: h="**(CreateParallelism)**"
    CHARACTER(LEN=4)  :: caux
    CHARACTER(LEN=14) :: nameDump

    !$    INTEGER,   EXTERNAL   ::  OMP_GET_MAX_THREADS


    ! inicialize MPI

    CALL MPI_INIT(ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE(caux,"(i4)") ierror
       CALL FatalError(h//" MPI_INIT returns "//caux)
    END IF

    ! number of MPI processes

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, maxNodes , ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE(caux,"(i4)") ierror
       CALL FatalError(h//" MPI_COMM_SIZE returns "//caux)
    END IF

    ! number of this process on MPI computation

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myId, ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE(caux,"(i4)") ierror
       CALL FatalError(h//" MPI_COMM_RANK returns "//caux)
    END IF

    ! process name

    CALL MPI_GET_PROCESSOR_NAME(Pname, CompName, ierror)
    IF (ierror /= MPI_SUCCESS) THEN
       WRITE(caux,"(i4)") ierror
       CALL FatalError(h//" MPI_GET_PROCESSOR_NAME retorna "//caux)
    END IF

    ! OMP parallelism

    numThreads = 1
    !$   numThreads = OMP_GET_MAX_THREADS()

    ! Number of processes for file name

    WRITE(cNProc ,"(i4.4)") maxNodes*numThreads
    WRITE(cThisProc,"(i4.4)") myId

    ! select unit for dumping

    unitDump = 30

    ! generate dump file name and open file

    nameDump="Dump."//cThisProc//"."//cNProc
    OPEN(unitDump, FILE=nameDump, STATUS='replace')

    ! Number of processes for file name

    WRITE(cNProc ,"(i4)") maxNodes
    WRITE(cThisProc,"(i4)") myId
    WRITE (caux,"(i4)") numThreads

    ! tell the world I'm alive

    CALL MsgOne(h," Process "//TRIM(ADJUSTL(cThisProc))//" ("//&
         &Pname(1:CompName)//") among "//TRIM(ADJUSTL(cNProc))//&
         &" processes with "//TRIM(ADJUSTL(caux))//&
         &" threads is alive")
    CALL MsgDump(h," Process "//TRIM(ADJUSTL(cThisProc))//" ("//&
         &Pname(1:CompName)//") among "//TRIM(ADJUSTL(cNProc))//&
         &" processes with "//TRIM(ADJUSTL(caux))//&
         &" threads is alive")
  END SUBROUTINE CreateParallelism


  SUBROUTINE CreateFourierGroup(mygroup_four,myid_four)
 
    INTEGER, INTENT(IN) :: mygroup_four 
    INTEGER, INTENT(IN) :: myid_four
    
    INTEGER :: ierror

    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,mygroup_four,myid_four,COMM_FOUR,ierror)

  END SUBROUTINE CreateFourierGroup



  SUBROUTINE Msg(unit, h, message)
    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=*), INTENT(IN) :: h, message
    
    INTEGER :: thisThread
    INTEGER, PARAMETER :: maxLineLen=128!stdout record size
    CHARACTER(LEN=LEN(h)+LEN(message)) :: fullMsg!full message
    INTEGER :: fullMsgLen,lineCount,lines,first,last

!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
    thisThread = 0
!$  thisThread = OMP_GET_THREAD_NUM()
    IF (thisThread == 0) then
       fullMsg=TRIM(h)//TRIM(message)
       fullMsgLen=len_Trim(fullMsg)
       lineCount=FullMsgLen/maxLineLen
       IF(lineCount*maxLineLen < fullMsgLen)THEN
          lineCount=lineCount +1
       END IF
       DO lines=1,lineCount
          first=(lines-1)*maxLineLen+1
          last =Min(lines*MaxLineLen,FullMsgLen)
          WRITE(unit,"(a)") fullMsg(first:last)
!$        CALL FLUSH(unit)       
       END DO
    END IF
  END SUBROUTINE Msg



  !*** Dump messsage only at dump file ***



  SUBROUTINE MsgDump(h, message)
    CHARACTER(LEN=*), INTENT(IN) :: h, message
    CALL Msg(unitDump, h, message)
  END SUBROUTINE MsgDump



  !*** Dump message at stdout and dump file ***



  SUBROUTINE MsgOut(h, message)
    CHARACTER(LEN=*), INTENT(IN) :: h, message
    CALL Msg(stdout, h, message)
  END SUBROUTINE MsgOut



  !*** Dump message at stdout and dump file ***



  SUBROUTINE MsgOne(h, message)
    CHARACTER(LEN=*), INTENT(IN) :: h, message
    IF (myId == 0) THEN
       CALL Msg(stdout, h, message)
    END IF
  END SUBROUTINE MsgOne



  !*** Dump error message everywhere and destroy parallelism ***



  SUBROUTINE FatalError(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER :: ierror=-1
    INTEGER :: ierr
    CHARACTER(LEN=11) :: h="**(ERROR)**"
    CALL MsgOut(h, message)
    CALL MsgDump(h, message)
    CALL MPI_ABORT(MPI_COMM_WORLD, ierror, ierr)
    STOP
  END SUBROUTINE FatalError



  !*** Destroy parallelism ***



  SUBROUTINE DestroyParallelism(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER ierror
    CHARACTER(LEN=24) :: h="**(DestroyParallelism)**"
    CALL MsgOne(h, " Process "//TRIM(ADJUSTL(cThisProc))//" "//message)
    CALL MsgDump(h, " Process "//TRIM(ADJUSTL(cThisProc))//" "//message)
    CLOSE(unitDump)
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    CALL MPI_FINALIZE(ierror)
  END SUBROUTINE DestroyParallelism
END MODULE Parallelism
