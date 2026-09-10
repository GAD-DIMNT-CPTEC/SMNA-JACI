!
!  $Author: bonatti $
!  $Date: 2007/09/18 18:07:15 $
!  $Revision: 1.2 $
!
PROGRAM TopoWaterPercNavy

   ! First Point of Initial Data is at South Pole and Greenwitch
   ! First Point of Output  Data is at North Pole and Greenwitch

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
            r4 = SELECTED_REAL_KIND(6) ! Kind for 32-bits Real Numbers

   INTEGER, PARAMETER :: ImBox=30, JmBox=30, NLonBox=72, NLatBox=36

   INTEGER, PARAMETER :: Imax=NLonBox*ImBox, Jmax=NLatBox*JmBox

   REAL (KIND=r4) :: Undef=-9999.0_r4

   INTEGER :: nj, ni, io, im, jo, jm, j, LRec, ios

   REAL (KIND=r4) :: dxy

   LOGICAL :: GrADS

   CHARACTER (LEN=8) :: VarNameT='TopoNavy'

   CHARACTER (LEN=9) :: VarNameW='WaterNavy'

   CHARACTER (LEN=11) :: VarNameBCs='navytm.form'

   CHARACTER (LEN=13) :: VarNameG='TopoWaterNavy'

   CHARACTER (LEN=17) :: VarName='TopoWaterPercNavy'

   CHARACTER (LEN=12) :: DirPreBCs='pre/databcs/'

   CHARACTER (LEN=12) :: DirPreOut='pre/dataout/'

   CHARACTER (LEN=528) :: DirMain

   REAL (KIND=r4), DIMENSION (ImBox,JmBox) :: Tp, Wp

   REAL (KIND=r4), DIMENSION (Imax,Jmax) :: Topog, Water

   INTEGER :: nferr=0    ! Standard Error Print Out
   INTEGER :: nfinp=5    ! Standard Read In
   INTEGER :: nfprt=6    ! Standard Print Out
   INTEGER :: nfclm=10   ! To Read Navy Topography Data
   INTEGER :: nftop=20   ! To Write Navy Topography Data
   INTEGER :: nfwat=30   ! To Write Water Percentage Data
   INTEGER :: nftpw=40   ! To Write GrADS Navy Topography and Water Percentage Data
   INTEGER :: nfctl=50   ! To Write Output Data Description

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

   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//DirPreBCs//VarNameBCs, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreBCs//VarNameBCs, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO nj=1,NLatBox
      jo=1+JmBox*(nj-1)
      jm=jo+JmBox-1
      DO ni=1,NLonBox
         io=1+ImBox*(ni-1)
         im=io+ImBox-1
         CALL GetTopographyBox (ImBox, JmBox, Tp, Wp)
         Topog(io:im,jo:jm)=Tp(1:ImBox,1:JmBox)
         Water(io:im,jo:jm)=Wp(1:ImBox,1:JmBox)
      END DO
   END DO
   CLOSE (UNIT=nfclm)

   CALL FlipMatrix (Topog,Imax,Jmax)
   CALL FlipMatrix (Water,Imax,Jmax)

   INQUIRE (IOLENGTH=LRec) topog(:,1)
   OPEN (UNIT=nftop, FILE=TRIM(DirMain)//DirPreOut//VarNameT//'.dat', &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, &
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
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, &
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
      OPEN (UNIT=nftpw, FILE=TRIM(DirMain)//DirPreOut//VarNameG//'.dat', &
            FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarNameG//'.dat', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      DO j=1,Jmax
         WRITE (UNIT=nftpw, REC=j) Topog(:,j)
      END DO
      DO j=1,Jmax
         WRITE (UNIT=nftpw, REC=j+Jmax) Water(:,j)
      END DO
      CLOSE (UNIT=nftpw)

      dxy=360.0_r4/REAL(Imax,r4)
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//VarNameG//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarNameG//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET ^'// &
            TRIM(DirMain)//DirPreOut//VarNameG//'.dat'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Navy Topography and Water Percentage'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F12.7,F10.7)') 'XDEF ', Imax, ' LINEAR ', &
                                                      0.5_r4*dxy, dxy
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F12.7,F10.7)') 'YDEF ', Jmax, ' LINEAR ', &
                                                      -90.0_r4+0.5_r4*dxy, dxy
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 2'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO 0 99 Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'WPER 0 99 Percentage of Water [%]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** TopoWaterPercNavy ENDS NORMALLY ***"

CONTAINS


SUBROUTINE GetTopographyBox (ImBox, JmBox, Tp, Wp)
 
   ! Reads Modified Navy Data
   ! Input Data Contains Only Terrain Height (Index 1)
   ! and Percentage of Water (Index 2)

   ! Input:
   ! ImBox - Number of Longitude Points in a Box
   ! JmBox - Number of Latitude Points in a Box

   ! Output:
   ! Tp(ImBox,JmBox) - Terrain Height
   ! Wp(ImBox,JmBox) - Percentage of Water

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: ImBox
   INTEGER, INTENT(IN) :: JmBox

   REAL (KIND=r4), INTENT(OUT) :: Tp(ImBox,JmBox)
   REAL (KIND=r4), INTENT(OUT) :: Wp(ImBox,JmBox)

   INTEGER :: lon
   INTEGER :: lat
   INTEGER :: i, j

   REAL (KIND=r4) :: zfct

   INTEGER :: itop(2,ImBox,JmBox)

   LOGICAL, PARAMETER :: IWant=.TRUE., IWant2=.TRUE.

   zfct=100.0_r4*(1200.0_r4/3937.0_r4)

   ! Lat - Northern Most Colatitude of Grid Box (5 Dg X 5 Dg)
   ! Lon - Western Most Longitude
   READ (UNIT=nfclm, FMT='(2I3)') lat, lon
   READ (UNIT=nfclm, FMT='(30I3)') itop

   lat=90-lat

   DO j=1,JmBox
     DO i=1,ImBox

       ! Corrections (By J.P.Bonatti, 18 Nov 1999)

       ! itop(1,i,j) - Is the Normalized Topography:
       !               Must Be Between 100 and 320
       ! itop(2,i,j) - Is The Percentage of Water:
       !               Must Be Between 0 and 100

       ! Correction Were Done Analysing The Surrounding Values

       IF (IWant) THEN

         ! Wrong Value
         IF (itop(1,i,j) == 330) THEN
           itop(1,i,j)=300
         END IF

         ! Wrong Value
         IF (itop(1,i,j) == 340) THEN
           itop(1,i,j)=300
         END IF

         ! Wrong Value
         IF (itop(1,i,j) == 357) THEN
           itop(1,i,j)=257
         END IF

         ! Wrong Value
         IF (itop(1,i,j) == 433) THEN
           itop(1,i,j)=113
         END IF

         ! May Be Missing Value
         IF (itop(1,i,j) == 511) THEN
           itop(1,i,j)=100
         END IF

         ! May Be Missing Value
         IF (itop(2,i,j) == 127) THEN
           IF (itop(1,i,j) == 100) THEN
             itop(2,i,j)=100
           ELSE
             itop(2,i,j)=1
           END IF
         END IF

         IF (IWant2) THEN

           ! To Avoid Negative Topography
           IF (itop(1,i,j) < 100) THEN
             itop(1,i,j)=100
           END IF

           ! To Avoid Non-Desirable Values of Surface
           ! Temperature Over Water Inside Continents
           IF (itop(1,i,j) > 100 .AND. itop(2,i,j) >= 50) THEN
             itop(2,i,j)=49
           END IF

         END IF

       END IF

       Tp(i,j)=zfct*REAL(itop(1,i,j)-100,r4)
       Wp(i,j)=REAL(itop(2,i,j),r4)

     END DO
   END DO

END SUBROUTINE GetTopographyBox


SUBROUTINE FlipMatrix (h, Imax, Jmax)
 
   ! Flips Over The Southern and Northern Hemispheres

   ! Input:
   ! h(Imax,Jmax) - Matrix to be Flipped
   !   Imax       - Column Dimension of h(Imax,Jmax)
   !   Jmax       - Row Dimension of h(Imax,Jmax)

   ! Output:
   ! h(Imax,Jmax) - Flipped Matrix

   IMPLICIT NONE

   REAL (KIND=r4), INTENT(INOUT) :: h(Imax,Jmax)
   INTEGER, INTENT(IN) :: Imax
   INTEGER, INTENT(IN) :: Jmax

   INTEGER :: Jmaxd, Jmaxd1

   REAL (KIND=r4) :: wk(Imax,Jmax)

   Jmaxd=Jmax/2
   Jmaxd1=Jmaxd+1

   wk=h
   h(:,1:Jmaxd)=wk(:,Jmax:Jmaxd1:-1)
   h(:,Jmaxd1:Jmax)=wk(:,Jmaxd:1:-1)

END SUBROUTINE FlipMatrix


END PROGRAM TopoWaterPercNavy
