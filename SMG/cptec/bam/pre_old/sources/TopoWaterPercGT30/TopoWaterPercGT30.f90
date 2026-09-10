!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM TopoWaterPercGT30

   ! First Point of Initial Data is at North Pole and I. D. Line
   ! First Point of Output  Data is at North Pole and Greenwhich

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER, PARAMETER :: NumBox=33, Imax=43200, Jmax=21600

   REAL (KIND=r4), PARAMETER :: TpMinSea=0.0_r4, TpMinLand=1.0_r4, &
                                WpMin=0.0_r4, WpMax=100_r4, &
                                Undef=-9999.0_r4

   INTEGER :: n, io, im, jo, jm, j, m, mx, mr, ms, ma, mb, &
              ix, jx, i1, i2, j1, j2, nr, LRecOut, LRecGad, ios

   REAL (KIND=r4) :: dxy, lono, lato, long, latg

   LOGICAL :: GrADS

   CHARACTER (LEN=3) :: gd

   CHARACTER (LEN=8) :: VarNameT='TopoGT30'

   CHARACTER (LEN=9) :: VarNameW='WaterGT30'

   CHARACTER (LEN=13) :: VarNameG='TopoWaterGT30'

   CHARACTER (LEN=17) :: VarName='TopoWaterPercGT30'

   CHARACTER (LEN=20) :: DirPreBCs='pre/databcs/GTOPO30/'

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=528) :: DirMain

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: Tp, Wp

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: Topog, Water
   
   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfclm=10   ! To Read GT30 Topography Data
   INTEGER :: nftop=20   ! To Write GT30 Topography Data
   INTEGER :: nfwat=30   ! To Write Water Percentage Data
   INTEGER :: nftpw=40   ! To Write GrADS GT30 Topography and Water Percentage Data
   INTEGER :: nfctl=50   ! To Write Output Data Description

   CHARACTER (LEN=7), DIMENSION(NumBox) :: TopoName = (/&
             'W180N90', 'W140N90', 'W100N90', 'W060N90', 'W020N90', &
             'E020N90', 'E060N90', 'E100N90', 'E140N90', 'W180N40', &
             'W140N40', 'W100N40', 'W060N40', 'W020N40', 'E020N40', &
             'E060N40', 'E100N40', 'E140N40', 'W180S10', 'W140S10', &
             'W100S10', 'W060S10', 'W020S10', 'E020S10', 'E060S10', &
             'E100S10', 'E140S10', 'W180S60', 'W120S60', 'W060S60', &
             'W000S60', 'E060S60', 'E120S60' /)

   INTEGER, DIMENSION(NumBox) :: ImBox = (/ &
            4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, &
            4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, &
            4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, &
            7200, 7200, 7200, 7200, 7200, 7200 /)

   INTEGER, DIMENSION(NumBox) :: JmBox = (/ &
            6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, &
            6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, &
            6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, &
            3600, 3600, 3600, 3600, 3600, 3600 /)

   NAMELIST /InputDim/ GrADS, DirMain

   GrADS=.TRUE.
   DirMain='./ '

   OPEN (UNIT=nfinp, FILE='./'//VarName//'.nml', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              './'//VarName//'.nml', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfinp, NML=InputDim)
   CLOSE (UNIT=nfinp)

   WRITE (UNIT=nfprt, FMT='(/,A)')  ' &InputDim'
   WRITE (UNIT=nfprt, FMT='(A,L6)') '    GrADS = ', GrADS
   WRITE (UNIT=nfprt, FMT='(A)')    '  DirMain = '//TRIM(DirMain)
   WRITE (UNIT=nfprt, FMT='(A,/)')  ' /'
   ALLOCATE(Topog(Imax,Jmax))
   ALLOCATE(Water(Imax,Jmax))
   DO n=1,NumBox
      ALLOCATE (Tp(ImBox(n),JmBox(n)), Wp(ImBox(n),JmBox(n)))
      OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//DirPreBCs//TopoName(n)//'.dat', &
            FORM='UNFORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='READ', STATUS='OLD', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreBCs//TopoName(n)//'.dat', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      READ  (UNIT=nfclm) Tp
      CLOSE (UNIT=nfclm)
      WHERE (Tp == Undef)
         ! Over Sea: Wp=100% and Tp=0m
         Wp=WpMax
         Tp=TpMinSea
      ELSEWHERE
         ! Over Land: Wp=0% and Tp=MAX(1m,Tp)
         Wp=WpMin
         Tp=MAX(TpMinLand,Tp)
      ENDWHERE
      SELECT CASE (n)
         CASE (1:9)
            io=1+(n-1)*ImBox(1)
            im=n*ImBox(1)
            jo=1
            jm=JmBox(1)
         CASE (10:18)
            io=1+(n-10)*ImBox(10)
            im=(n-9)*ImBox(10)
            jo=1+JmBox(1)
            jm=JmBox(1)+JmBox(10)
         CASE (19:27)
            io=1+(n-19)*ImBox(19)
            im=(n-18)*ImBox(19)
            jo=1+JmBox(1)+JmBox(10)
            jm=JmBox(1)+JmBox(10)+JmBox(19)
         CASE (28:NumBox)
            io=1+(n-28)*ImBox(28)
            im=(n-27)*ImBox(28)
            jo=1+JmBox(1)+JmBox(10)+JmBox(19)
            jm=JmBox(1)+JmBox(10)+JmBox(19)+JmBox(28)
      END SELECT
      Topog(io:im,jo:jm)=Tp(1:ImBox(n),1:JmBox(n))
      Water(io:im,jo:jm)=Wp(1:ImBox(n),1:JmBox(n))
      WRITE (UNIT=nfprt, FMT='(7I8)') n, io, im, jo, jm, ImBox(n), JmBox(n)
      WRITE (UNIT=nfprt, FMT='(1P2G12.5)') MINVAL(Tp), MAXVAL(Tp)
      WRITE (UNIT=nfprt, FMT='(1P2G12.5)') MINVAL(Wp), MAXVAL(Wp)
      DEALLOCATE (Tp, Wp)
   END DO

   CALL FlipMatrix (Topog,Imax,Jmax)
   CALL FlipMatrix (Water,Imax,Jmax)

   INQUIRE (IOLENGTH=LRecOut) Topog(:,1)
   OPEN (UNIT=nftop, FILE=TRIM(DirMain)//DirPreOut//VarNameT//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarNameT//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO j=1,Jmax
      WRITE (UNIT=nftop, REC=j) Topog(:,j)
   END DO
   CLOSE (UNIT=nftop)

   OPEN (UNIT=nfwat, FILE=TRIM(DirMain)//DirPreOut//VarNameW//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarNameW//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF 
   DO j=1,Jmax
      WRITE (UNIT=nfwat, REC=j) Water(:,j)
   END DO
   CLOSE (UNIT=nfwat)

   IF (GrADS) THEN
      mx=128
      mr=16
      ms=mx/mr
      ix=Imax/mr
      jx=Jmax/ms
      dxy=360.0_r4/REAL(Imax,r4)
      lono=0.5_r4*dxy
      lato=90.0_r4-0.5_r4*dxy
      m=0
      DO ma=1,mr
         i1=1+(ma-1)*ix
         i2=ma*ix
         long=lono+REAL(i1-1,r4)*dxy
         DO mb=1,ms
            m=m+1
            j1=1+(mb-1)*jx
            j2=mb*jx
            latg=lato-REAL(j2-1,r4)*dxy
            WRITE (gd, FMT='(I3.3)') m
            INQUIRE (IOLENGTH=LRecGad) Topog(i1:i2,1)
            OPEN (UNIT=nftpw, FILE=TRIM(DirMain)//DirPreOut//VarNameG//gd//'.dat', &
                  FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecGad, &
                  ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
            IF (ios /= 0) THEN
               WRITE (UNIT=nferr, FMT='(3A,I4)') &
                     ' ** (Error) ** Open file ', &
                       TRIM(DirMain)//DirPreOut//VarNameG//gd//'.dat', &
                     ' returned IOStat = ', ios
               STOP  ' ** (Error) **'
            END IF
            nr=0
            DO j=j1,j2
               nr=nr+1
               WRITE (UNIT=nftpw, REC=nr) Topog(i1:i2,j)
            END DO
            DO j=j1,j2
               nr=nr+1
               WRITE (UNIT=nftpw, REC=nr) Water(i1:i2,j)
            END DO
            CLOSE (UNIT=nftpw)
            WRITE (UNIT=nfprt, FMT='(A,4I6,2F14.7)') &
                  ' '//gd//' ', i1, i2, j1, j2, long, latg
            dxy=360.0_r4/REAL(Imax,r4)
            OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//VarNameG//gd//'.ctl', &
                  FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
                  STATUS='REPLACE', IOSTAT=ios)
            IF (ios /= 0) THEN
               WRITE (UNIT=nferr, FMT='(3A,I4)') &
                     ' ** (Error) ** Open file ', &
                       TRIM(DirMain)//DirPreOut//VarNameG//gd//'.ctl', &
                     ' returned IOStat = ', ios
               STOP  ' ** (Error) **'
            END IF
            WRITE (UNIT=nfctl, FMT='(A)') 'DSET ^'// &
                  TRIM(DirMain)//DirPreOut//VarNameG//gd//'.dat'
            WRITE (UNIT=nfctl, FMT='(A)') '*'
            WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
            WRITE (UNIT=nfctl, FMT='(A)') '*'
            WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
            WRITE (UNIT=nfctl, FMT='(A)') '*'
            WRITE (UNIT=nfctl, FMT='(A)') 'TITLE GT30 Topography and Water Percentage'
            WRITE (UNIT=nfctl, FMT='(A)') '*'
            WRITE (UNIT=nfctl, FMT='(A,I5,A,2F15.10)') 'XDEF ', ix, ' LINEAR ', long, dxy
            WRITE (UNIT=nfctl, FMT='(A,I5,A,2F15.10)') 'YDEF ', jx, ' LINEAR ', latg, dxy
            WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
            WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
            WRITE (UNIT=nfctl, FMT='(A)') '*'
            WRITE (UNIT=nfctl, FMT='(A)') 'VARS 2'
            WRITE (UNIT=nfctl, FMT='(A)') 'TOPO 0 99 Topography [m]'
            WRITE (UNIT=nfctl, FMT='(A)') 'WPER 0 99 Percentage of Water [%]'
            WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
            CLOSE (UNIT=nfctl)
         END DO
      END DO
   END IF
   DEALLOCATE(Topog)
   DEALLOCATE(Water)

PRINT *, "*** TopoWaterPercGT30 ENDS NORMALLY ***"

CONTAINS

SUBROUTINE FlipMatrix (h, Imax, Jmax)
 
   ! Flips a Matrix Over I.D.L. and Greenwitch

   ! Input:
   ! h(Imax,Jmax) - Matrix to be Flipped
   !   Imax       - Column Dimension of h(Imax,Jmax)
   !   Jmax       - Row Dimension of h(Imax,Jmax)

   ! Output:
   ! h(Imax,Jmax) - Flipped Matrix
   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER, INTENT(IN) :: Imax

   INTEGER, INTENT(IN) :: Jmax

   REAL (KIND=r4), INTENT(INOUT) :: h(Imax,Jmax)


   INTEGER :: Imaxd, Imaxd1, j

   REAL (KIND=r4) :: wk(Imax/2)

   Imaxd=Imax/2
   Imaxd1=Imaxd+1

   DO j=1,Jmax
      wk=h(Imaxd1:Imax,j)
      h(Imaxd1:Imax,j)=h(1:Imaxd,j)
      h(1:Imaxd,j)=wk
   END DO

END SUBROUTINE FlipMatrix
END PROGRAM TopoWaterPercGT30
