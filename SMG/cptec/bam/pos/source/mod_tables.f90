!
! Author: Luiz Flavio (PAD)
! 
!
MODULE TABLES

  !Vars
  INTEGER :: size_tb(3)
  TYPE tb1
    CHARACTER(LEN=4)  :: name
    CHARACTER(LEN=40) :: title
    CHARACTER(LEN=16) :: unit
    CHARACTER(LEN=5)  :: tipo
    CHARACTER(LEN=7)  :: level
    CHARACTER(LEN=2)  :: coment1
    CHARACTER(LEN=2)  :: coment2
    INTEGER           :: id
    INTEGER           :: dec_scal_fact
    INTEGER           :: precision

  END TYPE tb1
  
  TYPE(tb1),ALLOCATABLE,DIMENSION(:) :: table1
  
  !Levels
  TYPE tb2
    CHARACTER(LEN=7)  :: level_type
    CHARACTER(LEN=26) :: level_descr
    CHARACTER(LEN=8)  :: unit
    CHARACTER(LEN=6)  :: vert
    CHARACTER(LEN=4)  :: positive
    INTEGER           :: default
    INTEGER           :: id
    INTEGER           :: p1
    INTEGER           :: p2
  END TYPE tb2
  
  TYPE(tb2),ALLOCATABLE,DIMENSION(:) :: table2
  
  !Center
  TYPE tb3
    CHARACTER(LEN=10) :: center
    INTEGER           :: id
    INTEGER           :: grib_center
    INTEGER           :: sub_center
  END TYPE tb3
  
  TYPE(tb3),ALLOCATABLE,DIMENSION(:) :: table3
  
  LOGICAL :: tables_readed=.FALSE.
  PUBLIC :: Init_tables  
 CONTAINS

  SUBROUTINE Init_tables(datalib,rfd)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: datalib
    CHARACTER(LEN=*) :: rfd
    INTEGER :: ierr,i,k,is
    CHARACTER(LEN=200) :: lixo
    CHARACTER :: ic
    is=MAX(1, LEN_TRIM(datalib))

    DO k=1,3    
    
      WRITE(ic,FMT='(I1.1)') k

      OPEN(UNIT=63,FILE=TRIM(datalib)//'/'//'tab'//ic//'_'//TRIM(rfd(SCAN (rfd,'d')+2:))//'.dat',&
      ACCESS='sequential',STATUS='old')
    
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
    
      i=0;ierr=0
    
      DO WHILE (ierr==0)
        i=i+1
        READ(63,FMT='(A200)',IOSTAT=ierr) lixo
      END DO
      size_tb(k)=i-1 
    
      CLOSE(UNIT=63)
    END DO
    ALLOCATE(table1(size_tb(1)),table2(size_tb(2)),table3(size_tb(3)))
    CALL ReadTables(datalib,rfd)
  END SUBROUTINE Init_Tables
 
  SUBROUTINE ReadTables(datalib,rfd)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: datalib
    CHARACTER(LEN=*) :: rfd
    INTEGER :: k,i,is
    CHARACTER(LEN=200) :: lixo
    CHARACTER :: ic

    is=MAX(1, LEN_TRIM(datalib))

    DO k=1,3    
    
      WRITE(ic,FMT='(I1.1)') k
!
! bug 
!
      OPEN(UNIT=63,FILE=TRIM(datalib)//'/'//'tab'//ic//'_'//TRIM(rfd(SCAN (rfd,'d')+2:))//'.dat',&
      ACCESS='sequential',STATUS='old')
    
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
      READ (63,FMT='(A200)') lixo
 
      DO i=1,size_tb(k)
!        PRINT *,'i=',i
        IF(k==1) THEN
!           WRITE(*,*)'tab1'
!PRES | 001 | PRESSURE                                   | HPA              | float | sfc    | -999 | 16 |   |   |
!A4  3X I3 3X A40                                       5X A16             3X   A5 3X A7    2X   I4 3XI2 3X   A2  A2
           READ(63,FMT='(A4,3X,I3.3,3X,A40,5X,A16,3X,A5,3X,A7,2X,I4.4,3X,I2.2,3X,A2,3X,A2)')  &
                        table1(i)%name,         &
                        table1(i)%id,           &
                        table1(i)%title,        &
                        table1(i)%unit,         &
                        table1(i)%tipo,         &
                        table1(i)%level,        &
                        table1(i)%dec_scal_fact,&
                        table1(i)%precision,    &
                        table1(i)%coment1,      &
                        table1(i)%coment2
!PK           WRITE(*,FMT='(A4,3X,I3,3X,A40,5X,A16,3X,A5,3X,A7,2X,I4,3X,I2,3X,A2,3X,A2)')  &
!PK                        table1(i)%name,         &
!PK                        table1(i)%id,           &
!PK                        table1(i)%title,        &
!PK                        table1(i)%unit,         &
!PK                        table1(i)%tipo,         &
!PK                        table1(i)%level,        &
!PK                        table1(i)%dec_scal_fact,&
!PK                        table1(i)%precision,    &
!PK                        table1(i)%coment1,      &
!PK                         table1(i)%coment2
        ELSE IF(k==2) THEN
!sfc100m  | 100 meter above earth surf| m        | single |   up |  105 | 0 |  0 | 100
!A7         A26                         A8         A6         A4    I4    I2   I2  I3    
           READ(63,FMT='(A7,4X,A26,2X,A8,3X,A6,3X,A4,3X,I4.4,2X,I2.2,2X,I3.3,3X,I3.3)') &
                        table2(i)%level_type,   &
                        table2(i)%level_descr,  &
                        table2(i)%unit,         &
                        table2(i)%vert,         &
                        table2(i)%positive,     &
                        table2(i)%default,      &
                        table2(i)%id,           &
                        table2(i)%p1,           &
                        table2(i)%p2  
!           WRITE(*,*)'tab2'
!PK           WRITE(*,FMT='(A8,3X,A26,2X,A8,3X,A6,3X,A4,3X,I4.4,2X,I2.2,2X,I3.3,3X,I3.3)') &
!PK                        table2(i)%level_type,   &
!PK                        table2(i)%level_descr,  &
!PK                        table2(i)%unit,         &
!PK                        table2(i)%vert,         &
!PK                        table2(i)%positive,     &
!PK                        table2(i)%default,      &
!PK                        table2(i)%id,           &
!PK                        table2(i)%p1,           &
!PK                        table2(i)%p2  

        ELSE
           READ(63,FMT='(A10,2X,I3.3,3X,I3.3,3X,I1.1)') &
                        table3(i)%center,       &
                        table3(i)%id,           &
                        table3(i)%grib_center,  &
                        table3(i)%sub_center
        END IF

      END DO
    END DO  
    tables_readed=.TRUE.

  END SUBROUTINE ReadTables
  
END MODULE TABLES
