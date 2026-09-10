MODULE Watches

  ! Execution time instrumentation
  ! Measures execution time of selected parts of the program,
  ! associating one watch to each part. The associated watch
  ! (an integer) should be turned on as execution enters the
  ! selected part and turned off as execution leaves the part.
  ! To turn watch "tp" on, call ChangeWatch(tp). There is no
  ! procedure to turn a watch off; instead, turn on another
  ! watch (invoking ChangeWatch). If needed, create a default
  ! watch to be turned on, for parts of the code where measurement
  ! is not important.
  ! 
  ! Usage: 
  !   CreateWatches(nW, firstW)
  !     creates nW watches; starts watch firstW
  !   NameWatch(tp, name)
  !     gives name "name" to watch tp
  !   ChangeWatch (tp)
  !     turn off current watch; turn on watch tp (0<=tp<=nW)
  !   DumpWatches (unit)
  !     dump all watches on fortran io unit "unit"
  !   DestroyWatches()
  !     destroy all watches
  !
  ! Execution order:
  !   create watches and  name each watch at the beginning of main program;
  !   change watch at the border of

  USE Parallelism, ONLY:&
       FatalError, myid, maxnodes


   IMPLICIT NONE
  SAVE       

  PRIVATE
  PUBLIC :: CreateWatches
  PUBLIC :: NameWatch
  PUBLIC :: ChangeWatch
  PUBLIC :: DumpWatches
  PUBLIC :: DestroyWatches


  INCLUDE 'mpif.h'
  INTEGER                                      :: nWatches 
  INTEGER                                      :: currWatch, ierr
  REAL             , ALLOCATABLE, DIMENSION(:) :: TimeCpu
  REAL             , ALLOCATABLE, DIMENSION(:) :: TimeWall
  REAL             , ALLOCATABLE, DIMENSION(:) :: Tmin, tmax, tmean
  CHARACTER(LEN=16), ALLOCATABLE, DIMENSION(:) :: Name
  REAL                                         :: timeUnit
CONTAINS



  ! CriateWatches:
  !   Criates "nW" watches (nW > 0)
  !   Starts watch "firstW"  (firstW <= nW)



  SUBROUTINE CreateWatches(nW, firstW)
    INTEGER, INTENT(IN) :: nW
    INTEGER, INTENT(IN) :: firstW

    INTEGER             :: ierr
    INTEGER             :: cnt
    INTEGER             :: cntRate
    REAL                :: tc, tw
    CHARACTER(LEN=*), PARAMETER :: h='**(CreateWatches)**'
    CHARACTER(LEN=10)   :: c0, c1
    ierr=0
    cnt=0
    cntRate=0
    tc=0.0; tw=0.0
    ! check input arguments

    IF (nW <=0) THEN
       WRITE(c0,"(i10)") nW
       CALL FatalError(h//" number of watches ("//TRIM(ADJUSTL(c0))//") should be > 0")
       STOP 1
    ELSE IF (firstW > nW) THEN
       WRITE(c0,"(i10)") nW
       WRITE(c1,"(i10)") firstW
       CALL FatalError(h//" first watch ("//TRIM(ADJUSTL(c1))//&
            &") > number of watches ("//TRIM(ADJUSTL(c0))//")")
       STOP 1
    END IF

    ! initialize module data

    nWatches=nW
    ALLOCATE (TimeCpu (0:nWatches), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" allocate TimeCpu returned error "//TRIM(ADJUSTL(c0)))
    END IF
    ALLOCATE (Tmin (0:nWatches), STAT=ierr)
    ALLOCATE (Tmax (0:nWatches), STAT=ierr)
    ALLOCATE (Tmean(0:nWatches), STAT=ierr)
    TimeCpu = 0.0
    Tmin = 0.0
    Tmax = 0.0
    Tmean= 0.0

    ALLOCATE (TimeWall(0:nWatches), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" allocate TimeWall returned error "//TRIM(ADJUSTL(c0)))
    END IF
    TimeWall= 0.0

    ALLOCATE (Name    (0:nWatches), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" allocate Name returned error "//TRIM(ADJUSTL(c0)))
    END IF
    Name = '          '
    Name(0) = "Total "

    ! get current cpu and wall time

    CALL Cpu_Time(tc)
    CALL System_Clock(cnt, cntRate)
    IF (cntRate == 0) THEN
       CALL FatalError(h//" no clock at this instalation!!!")
    END IF
    timeUnit = 1.0/REAL(cntRate)
    tw = REAL(cnt)*timeUnit

    ! turn on total watch

    TimeCpu(0) = - tc
    TimeWall(0) = - tw

    ! turn on currWatch

    currWatch=firstW
    TimeCpu (currWatch)  = - tc
    TimeWall(currWatch)  = - tw
  END SUBROUTINE CreateWatches



  ! Destroy all watches



  SUBROUTINE DestroyWatches()
    CHARACTER(LEN=*), PARAMETER :: h='**(DestroyWatches)**'
    INTEGER                     :: ierr
    CHARACTER(LEN=10)           :: c0

    DEALLOCATE (TimeCpu, STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" deallocate TimeCpu returned error "//TRIM(ADJUSTL(c0)))
    END IF

    DEALLOCATE (TimeWall, STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" deallocate TimeWall returned error "//TRIM(ADJUSTL(c0)))
    END IF

    DEALLOCATE (Name, STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" deallocate Name returned error "//TRIM(ADJUSTL(c0)))
    END IF
  END SUBROUTINE DestroyWatches



  ! Dump Watches into unit 'unit' or into stdout if "unit" not present



  SUBROUTINE DumpWatches(unit)
    INTEGER, INTENT(IN), OPTIONAL :: unit
    !$ INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

    CHARACTER(LEN=*), PARAMETER :: h='**(DumpWatches)**'
    INTEGER tp, cnt, numthreads
    REAL :: tc, tw

    ! get current cpu and wall time

    CALL Cpu_Time(tc)
    CALL System_Clock(cnt)
    tw = REAL(cnt)*timeUnit

    ! turn off watches 0 and currWatch during dump

    TimeCpu (        0)  = TimeCpu (        0) + tc
    TimeCpu (currWatch)  = TimeCpu (currWatch) + tc
    TimeWall(        0)  = TimeWall(        0) + tw
    TimeWall(currWatch)  = TimeWall(currWatch) + tw

    ! dump

     CALL MPI_REDUCE(TimeCpu,tmax,nWatches+1,MPI_REAL, MPI_MAX, 0,&
                                MPI_COMM_WORLD,ierr)
     CALL MPI_REDUCE(TimeCpu,tmin,nWatches+1,MPI_REAL, MPI_MIN, 0,&
                                MPI_COMM_WORLD,ierr)
     CALL MPI_REDUCE(TimeCpu,tmean,nWatches+1,MPI_REAL, MPI_SUM, 0,&
                                MPI_COMM_WORLD,ierr)
     tmean = tmean / maxnodes
    numthreads = 1
    !$ numthreads = OMP_GET_MAX_THREADS()
    IF (PRESENT(unit)) THEN
       WRITE(unit,"(a,'     TimeCpu (mean per thread)  TimeWall')") h
       DO tp = 1, nWatches
          WRITE (unit,"(2x,a,3f12.3)") &
               Name(tp), TimeCpu(tp), TimeCpu(tp)/numthreads, TimeWall(tp)
       END DO
       WRITE (unit,"(2x,a,3f12.3)") &
            Name( 0), TimeCpu( 0), TimeCpu( 0)/numthreads, TimeWall( 0)
!      CALL FLUSH(unit)
    ELSE
       WRITE(*,"(a,'     TimeCpu (mean per thread)  TimeWall')") h
       DO tp = 1, nWatches
          WRITE (*,"(2x,a,3f12.3)") &
               Name(tp), TimeCpu(tp), TimeCpu(tp)/numthreads, TimeWall(tp)
       END DO
       WRITE (*,"(2x,a,3f12.3)") &
            Name( 0), TimeCpu( 0), TimeCpu( 0)/numthreads, TimeWall( 0)
    END IF
    IF (myid.eq.0) THEN
      IF (PRESENT(unit)) THEN
         WRITE(unit,"(a,'   Max time     Min Time      Average     Imbalance  %')") h
         DO tp = 1, nWatches
            WRITE (unit,"(2x,a,4f12.3)") Name(tp), Tmax(tp),&
                 Tmin(tp),tmean(tp),(Tmax(tp)-tmean(tp))/MAX(tmean(tp),1.0e-20)*100 
         END DO
         CALL FLUSH(unit)
      ELSE
         WRITE(*,"(a,'   Max time     Min Time     Average     Imbalance % ') ") h
         DO tp = 1, nWatches
            WRITE (*,"(2x,a,4f12.3)") Name(tp), Tmax(tp),&
                 Tmin(tp),tmean(tp),(Tmax(tp)-tmean(tp))/MAX(tmean(tp),1.0e-20)*100 
         END DO
      END IF
    END IF


    ! restore watches 0 and currWatch

    TimeCpu (        0)  = TimeCpu (        0) - tc
    TimeCpu (currWatch)  = TimeCpu (currWatch) - tc
    TimeWall(        0)  = TimeWall(        0) - tw
    TimeWall(currWatch)  = TimeWall(currWatch) - tw
  END SUBROUTINE DumpWatches



  ! Make watch "tp" the current watch



  SUBROUTINE ChangeWatch(tp)
    INTEGER, INTENT(IN) :: tp

    INTEGER :: cnt
    REAL :: tc, tw
    CHARACTER(LEN=*), PARAMETER :: h='**(ChangeWatch)**'
    CHARACTER(LEN=10) :: c0, c1

    ! avoid turning on current watch;
    ! check if input argument in bounds

    IF (tp == currWatch) THEN
       RETURN
    ELSE IF (tp <= 0 .OR. tp > nWatches) THEN
       WRITE(c0,"(i10)") tp
       WRITE(c1,"(i10)") nWatches
       CALL FatalError(h//" input watch ("//TRIM(ADJUSTL(c0))//&
            &") should be > 0 and < number of watches ("//&
            &TRIM(ADJUSTL(c1)))
    END IF

    ! get current cpu and wall time

    CALL Cpu_Time(tc)
    CALL System_Clock(cnt)
    tw = REAL(cnt)*timeUnit

    ! turn off currWatch

    TimeCpu (currWatch)  = TimeCpu (currWatch) + tc
    TimeWall(currWatch)  = TimeWall(currWatch) + tw

    ! change currWatch and turn it on

    currWatch=tp
    TimeCpu (currWatch)  = TimeCpu (currWatch) - tc
    TimeWall(currWatch)  = TimeWall(currWatch) - tw
  END SUBROUTINE ChangeWatch


  ! Set up name for watch


  SUBROUTINE NameWatch (tp, nameW)
    INTEGER, INTENT(IN) :: tp 
    CHARACTER(LEN=*), INTENT(IN) :: nameW 

    CHARACTER(LEN=*), PARAMETER :: h='**(NameWatch)**'
    CHARACTER(LEN=10) :: c0, c1

    ! check if input argument in bounds

    IF (tp <= 0 .OR. tp > nWatches) THEN
       WRITE(c0,"(i10)") tp
       WRITE(c1,"(i10)") nWatches
       CALL FatalError(h//" input watch ("//TRIM(ADJUSTL(c0))//&
            &") should be > 0 and < number of watches ("//&
            &TRIM(ADJUSTL(c1)))
    END IF
    Name(tp) = nameW
  END SUBROUTINE NameWatch
END MODULE Watches
