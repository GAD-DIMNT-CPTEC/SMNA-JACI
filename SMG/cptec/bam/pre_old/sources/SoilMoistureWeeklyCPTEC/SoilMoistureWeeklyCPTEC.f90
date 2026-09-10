!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM SoilMoistureWeeklyCPTEC

   ! First Point of Initial Data is near South Pole and near Greenwhich
   ! First Point of Output  Data is near North Pole and near Greenwhich

   IMPLICIT NONE

   INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers
   INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers

   INTEGER :: Idim, Jdim, kdim,Month, LRec, ios,k,i,j

   LOGICAL :: GrADS

   CHARACTER (LEN=11) :: DirPreIn='pre/datain/'

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=18) :: VarName='SoilMoistureWeekly'

!   CHARACTER (LEN=17) :: FileSoilMoisture='gdas1.T  Z.snogrd'
!                                         'GL_SM.GPNR.2014032912.gra'
   CHARACTER (LEN=25) :: FileSoilMoisture='GL_SM.GPNR.          .vfm'

   CHARACTER (LEN=12) :: TimeGrADS='  Z         '

   CHARACTER (LEN=27) :: NameNML='SoilMoistureWeeklyCPTEC.nml'

   CHARACTER (LEN=10) :: Date

   CHARACTER (LEN=528) :: DirMain

   CHARACTER (LEN=3), DIMENSION(12) :: MonthChar = &
             (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC' /)

   REAL (KIND=r4), DIMENSION (:,:,:), ALLOCATABLE :: SoilMoisture 
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SoilMoisture_out 

   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfsno=10   ! To Read Unformatted Weekly SoilMoisture Data
   INTEGER :: nfout=20   ! To Write Unformatted Weekly SoilMoisture Data
   INTEGER :: nfctl=30   ! To Write Output Data Description

   NAMELIST /InputDim/ Idim, Jdim,kdim, GrADS, Date, DirMain

   Idim=1440
   Jdim=720
   kdim=8
   GrADS=.TRUE.
   Date='          '
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
   WRITE (UNIT=nfprt, FMT='(A,I6)') '    Idim = ', Idim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '    Jdim = ', Jdim
   WRITE (UNIT=nfprt, FMT='(A,I6)') '    Kdim = ', Kdim
   WRITE (UNIT=nfprt, FMT='(A,L6)') '   GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')    '    Date = '//Date
   WRITE (UNIT=nfprt, FMT='(A)')    ' DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'

   ALLOCATE (SoilMoisture(Idim,Jdim,kdim))
   ALLOCATE (SoilMoisture_out(Idim,Jdim))

   FileSoilMoisture(12:21)=Date(1:10)! hour
   TimeGrADS(1:2)       =Date(9:10)! hour
   TimeGrADS(4:5)       =Date(7:8) ! day
   TimeGrADS(9:12)      =Date(1:4) ! year
   READ (Date(5:6), FMT='(I2)') Month
   TimeGrADS(6:8)=MonthChar(Month)


   OPEN (UNIT=nfsno, FILE=TRIM(DirMain)//TRIM(DirPreIn)//TRIM(FileSoilMoisture), &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
       ' ** (Error) ** Open file ', &
         TRIM(DirMain)//TRIM(DirPreIn)//TRIM(FileSoilMoisture), &
       ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
  

    !DO k=1,kdim
       !
       !print*,'================== GLSM for k====================',k
       !
       CALL vfirec(nfsno,SoilMoisture(1,1,1),Idim*jdim*kdim,'LIN')


    !END DO

   CLOSE (UNIT=nfsno)
 

   INQUIRE (IOLENGTH=LRec)SoilMoisture (1:Idim,1:Jdim,1)
   OPEN (UNIT=nfout, FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.'//Date(1:8), &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.'//Date(1:8), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO k=1,kdim
       DO j=1,Jdim
          DO i=1,Idim
             SoilMoisture_out(i,j)=MAX(SoilMoisture(i,Jdim+1-j,k),0.0_r4)
          END DO
       END DO
       SoilMoisture_out=CSHIFT (SoilMoisture_out, SHIFT=SIZE(SoilMoisture_out,1)/2,DIM = 1) 
      WRITE (UNIT=nfout, REC=k)SoilMoisture_out
   END DO
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.'//Date(1:8)//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.'//Date(1:8)//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.'//Date(1:8)
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'UNDEF -999.0'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE CPTEC Weekly SoilMoisture'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Idim,' LINEAR ',0.5_r4,360.0_r4/REAL(Idim,r4)
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'YDEF ',Jdim,' LINEAR ',-89.5_r4,180.0_r4/REAL(Jdim-1,r4)
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 8 LEVELS 1 2 3 4 5 6 7 8'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR '//TimeGrADS//' 6HR'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 1'
      WRITE (UNIT=nfctl, FMT='(A)') 'SoilMoisture 8 99 CPTEC Weekly SoilMoisture [kg/m2]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** SoilMoistureWeeklyCPTEC ENDS NORMALLY ***"

CONTAINS

!
!------------------------------- VFORMAT ----------------------------------
!
 SUBROUTINE vfirec(iunit,aa,n,type)

  INTEGER, INTENT(IN)  :: iunit  !#TO deve ser kind default
  INTEGER, INTENT(IN)  :: n
  REAL(KIND=r4), INTENT(OUT)    :: aa(n)
  CHARACTER(len=* ), INTENT(IN) :: type

  !
  ! local
  !
  CHARACTER(len=1 ) :: vc(0:63)
  CHARACTER(len=80) :: line
  CHARACTER(len=1 ) :: cs
  INTEGER           :: ich0
  INTEGER           :: ich9
  INTEGER           :: ichcz
  INTEGER           :: ichca
  INTEGER           :: ichla
  INTEGER           :: ichlz
  INTEGER           :: i
  INTEGER           :: nvalline
  INTEGER           :: nchs
  INTEGER           :: ic
  INTEGER           :: ii
  INTEGER           :: isval
  INTEGER           :: iii
  INTEGER           :: ics
  INTEGER           :: nn
  INTEGER           :: nbits
  INTEGER           :: nc
  REAL(KIND=r8)     :: bias
  REAL(KIND=r8)     :: fact
  REAL(KIND=r8)     :: facti
  REAL(KIND=r8)     :: scfct
  REAL(KIND=r8)     :: a(n)

  vc='0'
  IF (vc(0).ne.'0') CALL vfinit(vc)

  ich0 =ichar('0')
  ich9 =ichar('9')
  ichcz=ichar('Z')
  ichlz=ichar('z')
  ichca=ichar('A')
  ichla=ichar('a')

  READ (iunit,'(2i8,2e20.10)')nn,nbits,bias,fact

  IF (nn.ne.n) THEN
    PRINT*,' Word count mismatch on vfirec record '
    PRINT*,' Words on record - ',nn
    PRINT*,' Words expected  - ',n
    STOP 'vfirec'
  END IF

  nvalline=(78*6)/nbits
  nchs=nbits/6

  DO i=1,n,nvalline
    READ(iunit,'(a78)') line
    ic=0
    DO ii=i,i+nvalline-1
      isval=0
      IF(ii.gt.n) EXIT
      DO iii=1,nchs
         ic=ic+1
         cs=line(ic:ic)
         ics=ichar(cs)
         IF (ics.le.ich9) THEN
            nc=ics-ich0
         ELSE IF (ics.le.ichcz) THEN
            nc=ics-ichca+10
         ELSE
            nc=ics-ichla+36
         END IF
         isval=ior(ishft(nc,6*(nchs-iii)),isval)
      END DO ! loop iii
        a(ii)=isval
    END DO ! loop ii

  END DO ! loop i

  facti=1.0_r8/fact

  IF (type.eq.'LIN') THEN
    DO i=1,n

      a(i)=a(i)*facti-bias

      !print*,'VFM=',i,a(i)
    END DO
  ELSE IF (type.eq.'LOG') THEN
    scfct=2.0_r8**(nbits-1)
    DO i=1,n
        a(i)=sign(1.0_r8,a(i)-scfct)  &
           *(10.0_r8**(abs(20.0_r8*(a(i)/scfct-1.0_r8))-10.0_r8))
    END DO
  END IF
  aa=REAL(a,kind=r4)
 END SUBROUTINE vfirec

!--------------------------------------------------------
 SUBROUTINE vfinit(vc)
   CHARACTER(len=1), INTENT(OUT  ) :: vc   (*)
   CHARACTER(len=1)                :: vcscr(0:63)
   INTEGER                         :: n

   DATA vcscr/'0','1','2','3','4','5','6','7','8','9'   &
              ,'A','B','C','D','E','F','G','H','I','J'  &
              ,'K','L','M','N','O','P','Q','R','S','T'  &
              ,'U','V','W','X','Y','Z','a','b','c','d'  &
              ,'e','f','g','h','i','j','k','l','m','n'  &
              ,'o','p','q','r','s','t','u','v','w','x'  &
              ,'y','z','{','|'/

  DO n=0,63
      vc(n)=vcscr(n)
  END DO
 END SUBROUTINE vfinit


END PROGRAM SoilMoistureWeeklyCPTEC
