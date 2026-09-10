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
  !   CreateWatch(nW, firstW)
  !     creates nW watches; starts watch firstW
  !   NameWatch(tp, name)
  !     gives name "name" to watch tp
  !   ChangeWatch (tp)
  !     turn off current watch; turn on watch tp (0<=tp<=nW)
  !   DumpWatch (unit)
  !     dump all watches on fortran io unit "unit"
  !   DestroyWatch()
  !     destroy all watches
  !
  ! Execution order:
  !   create watches and  name each watch at the beginning of main program;
  !   change watch at the border of


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Watch
  PUBLIC :: CreateWatch
  PUBLIC :: NameWatch
  PUBLIC :: ChangeWatch
  PUBLIC :: DumpWatch
  PUBLIC :: DestroyWatch


  type Watch
     INTEGER                    :: nWatches 
     INTEGER                    :: currWatch
     REAL             , POINTER :: TimeCpu(:)
     REAL             , POINTER :: TimeWall(:)
     CHARACTER(LEN=16), POINTER :: Name(:)
     REAL                       :: timeUnit
  end type Watch
CONTAINS



  ! CriateWatch:
  !   Criates "nW" watches (nW > 0)
  !   Starts watch "firstW"  (firstW <= nW)



  FUNCTION CreateWatch(nW, firstW) result(wt)
    INTEGER, INTENT(IN) :: nW
    INTEGER, INTENT(IN) :: firstW
    TYPE(Watch)         :: wt

    INTEGER             :: ierr
    INTEGER             :: cnt
    INTEGER             :: cntRate
    REAL                :: tc, tw
    CHARACTER(LEN=*), PARAMETER :: h='**(CreateWatch)**'
    CHARACTER(LEN=10)   :: c0, c1

    ! check input arguments

    IF (nW <=0) THEN
       WRITE(c0,"(i10)") nW
       CALL FatalError(h//" number of watches ("//TRIM(ADJUSTL(c0))//") should be > 0")
       STOP 1
    ELSE IF (firstW > nW) THEN
       WRITE(c0,"(i10)") nW
       WRITE(c1,"(i10)") firstW
       WRITE(0,'(A)')
       CALL FatalError(h//" first watch ("//TRIM(ADJUSTL(c1))//&
            &") > number of watches ("//TRIM(ADJUSTL(c0))//")")
       STOP 1
    END IF

    ! initialize module data

    wt%nWatches=nW
    ALLOCATE (wt%TimeCpu (0:nW), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" allocate TimeCpu returned error "//TRIM(ADJUSTL(c0)))
    END IF
    wt%TimeCpu = 0.0

    ALLOCATE (wt%TimeWall(0:nW), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" allocate TimeWall returned error "//TRIM(ADJUSTL(c0)))
    END IF
    wt%TimeWall= 0.0

    ALLOCATE (wt%Name    (0:nW), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" allocate Name returned error "//TRIM(ADJUSTL(c0)))
    END IF
    wt%Name = '          '
    wt%Name(0) = "Total "

    ! get current cpu and wall time

    CALL Cpu_Time(tc)
    CALL System_Clock(cnt, cntRate)
    IF (cntRate == 0) THEN
       CALL FatalError(h//" no clock at this instalation!!!")
    END IF
    wt%timeUnit = 1.0/REAL(cntRate)
    tw = REAL(cnt)*wt%timeUnit

    ! turn on total watch

    wt%TimeCpu(0) = - tc
    wt%TimeWall(0) = - tw

    ! turn on currWatch

    wt%currWatch=firstW
    wt%TimeCpu (wt%currWatch)  = - tc
    wt%TimeWall(wt%currWatch)  = - tw
  END FUNCTION CreateWatch



  ! Destroy all watches



  SUBROUTINE DestroyWatch(wt)
    TYPE(Watch), INTENT(INOUT) :: wt

    CHARACTER(LEN=*), PARAMETER :: h='**(DestroyWatch)**'
    INTEGER                     :: ierr
    CHARACTER(LEN=10)           :: c0

    DEALLOCATE (wt%TimeCpu, STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" deallocate TimeCpu returned error "//TRIM(ADJUSTL(c0)))
    END IF

    DEALLOCATE (wt%TimeWall, STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" deallocate TimeWall returned error "//TRIM(ADJUSTL(c0)))
    END IF

    DEALLOCATE (wt%Name, STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(c0,"(i10)") ierr
       CALL FatalError(h//" deallocate Name returned error "//TRIM(ADJUSTL(c0)))
    END IF
  END SUBROUTINE DestroyWatch



  ! Dump Watches into unit 'unit' or into stdout if "unit" not present



  SUBROUTINE DumpWatch(wt, unit,NAME)
    TYPE(Watch), INTENT(INOUT) :: wt
    INTEGER, INTENT(IN), OPTIONAL :: unit
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: NAME
    CHARACTER(LEN=*), PARAMETER :: h='**(DumpWatch)**'
    INTEGER tp, cnt
    REAL :: tc, tw

    ! get current cpu and wall time

    CALL Cpu_Time(tc)
    CALL System_Clock(cnt)
    tw = REAL(cnt)*wt%timeUnit

    ! turn off watches 0 and currWatch during dump

    wt%TimeCpu (           0)  = wt%TimeCpu (           0) + tc
    wt%TimeCpu (wt%currWatch)  = wt%TimeCpu (wt%currWatch) + tc
    wt%TimeWall(           0)  = wt%TimeWall(           0) + tw
    wt%TimeWall(wt%currWatch)  = wt%TimeWall(wt%currWatch) + tw

    ! dump

    IF (PRESENT(unit)) THEN
       OPEN(unit,FILE=TRIM(NAME),STATUS='UNKNOWN')
       WRITE(unit,"(a,'       TimeCpu   TimeWall')") h
       DO tp = 1, wt%nWatches
          WRITE (unit,"(2x,a,2f12.3)") &
               wt%Name(tp), wt%TimeCpu(tp), wt%TimeWall(tp)
       END DO
       WRITE (unit,"(2x,a,2f12.3)") &
            wt%Name( 0), wt%TimeCpu( 0), wt%TimeWall( 0)
    ELSE
       WRITE(*,"(a,'          TimeCpu   TimeWall')") h
       DO tp = 1, wt%nWatches
          WRITE (*,"(2x,a,2f12.3)") &
               wt%Name(tp), wt%TimeCpu(tp), wt%TimeWall(tp)
       END DO
       WRITE (*,"(2x,a,2f12.3)") &
            wt%Name( 0), wt%TimeCpu( 0), wt%TimeWall( 0)
    END IF

    ! restore watches 0 and currWatch

    wt%TimeCpu (           0)  = wt%TimeCpu (           0) - tc
    wt%TimeCpu (wt%currWatch)  = wt%TimeCpu (wt%currWatch) - tc
    wt%TimeWall(           0)  = wt%TimeWall(           0) - tw
    wt%TimeWall(wt%currWatch)  = wt%TimeWall(wt%currWatch) - tw
  END SUBROUTINE DumpWatch



  ! Make watch "tp" the current watch



  SUBROUTINE ChangeWatch(wt, tp)
    TYPE(Watch), INTENT(INOUT) :: wt
    INTEGER, INTENT(IN) :: tp

    INTEGER :: cnt
    REAL :: tc, tw
    CHARACTER(LEN=*), PARAMETER :: h='**(ChangeWatch)**'
    CHARACTER(LEN=10) :: c0, c1

    ! avoid turning on current watch;
    ! check if input argument in bounds

    IF (tp == wt%currWatch) THEN
       RETURN
    ELSE IF (tp <= 0 .OR. tp > wt%nWatches) THEN
       WRITE(c0,"(i10)") tp
       WRITE(c1,"(i10)") wt%nWatches
       CALL FatalError(h//" input watch ("//TRIM(ADJUSTL(c0))//&
            &") should be > 0 and < number of watches ("//&
            &TRIM(ADJUSTL(c1))//")")
    END IF

    ! get current cpu and wall time

    CALL Cpu_Time(tc)
    CALL System_Clock(cnt)
    tw = REAL(cnt)*wt%timeUnit

    ! turn off currWatch

    wt%TimeCpu (wt%currWatch)  = wt%TimeCpu (wt%currWatch) + tc
    wt%TimeWall(wt%currWatch)  = wt%TimeWall(wt%currWatch) + tw

    ! change currWatch and turn it on

    wt%currWatch=tp
    wt%TimeCpu (wt%currWatch)  = wt%TimeCpu (wt%currWatch) - tc
    wt%TimeWall(wt%currWatch)  = wt%TimeWall(wt%currWatch) - tw
  END SUBROUTINE ChangeWatch


  ! Set up name for watch


  SUBROUTINE NameWatch (wt, tp, nameW)
    TYPE(Watch), INTENT(INOUT) :: wt
    INTEGER, INTENT(IN) :: tp 
    CHARACTER(LEN=*), INTENT(IN) :: nameW 

    CHARACTER(LEN=*), PARAMETER :: h='**(NameWatch)**'
    CHARACTER(LEN=10) :: c0, c1

    ! check if input argument in bounds

    IF (tp <= 0 .OR. tp > wt%nWatches) THEN
       WRITE(c0,"(i10)") tp
       WRITE(c1,"(i10)") wt%nWatches
       CALL FatalError(h//" input watch ("//TRIM(ADJUSTL(c0))//&
            &") should be > 0 and < number of watches ("//&
            &TRIM(ADJUSTL(c1))//")")
    END IF
    wt%Name(tp) = nameW
  END SUBROUTINE NameWatch


 SUBROUTINE FatalError(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER :: ierror=-1
    INTEGER :: ierr
    CHARACTER(LEN=10) :: h="**(ERROR)**"
    WRITE(0,'(A)')message
    STOP
  END SUBROUTINE FatalError

END MODULE Watches
