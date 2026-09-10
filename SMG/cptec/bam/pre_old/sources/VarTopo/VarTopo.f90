!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM VarTopo

   USE InputParameters, ONLY: r4, r8, Undef, VarMin, VarMax, &
                              nferr, nfclm, nfoua, nfoub, nfouc, nfctl, &
                              Imax, Jmax, Idim, Jdim, &
                              nLats, VarNameT, VarNameV, DirPreOut, DirModelIn, &
                              DirMain, VarName, VarNameG, GrADS, Linear, &
                              InitInputParameters

   USE LinearInterpolation, ONLY: gLatsL=>LatOut, &
       InitLinearInterpolation, DoLinearInterpolation

   USE AreaInterpolation, ONLY: gLatsA=>gLats, &
       InitAreaInterpolation, DoAreaInterpolation

   IMPLICIT NONE

   ! Horizontal Area Interpolator
   ! Interpolate Regular To Gaussian
   ! Regular Input Data is Assumed to be Oriented with
   ! the North Pole and Greenwich as the First Point

   ! Set Undefined Value for Input Data at Locations which
   ! are not to be Included in Interpolation

   INTEGER :: m, j, LRecIn, LRecOut, ios

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: TopoIn
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: FieldOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: TopoInput
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: TopoOutput
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: VarTopInput
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: VarTopOutput

   CALL InitInputParameters ()

   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (TopoIn(Idim,Jdim))
   ALLOCATE (FieldOut(Imax,Jmax))
   ALLOCATE (TopoInput(Idim,Jdim))
   ALLOCATE (TopoOutput(Imax,Jmax))
   ALLOCATE (VarTopInput(Idim,Jdim))
   ALLOCATE (VarTopOutput(Imax,Jmax))

   ! Read In Input Topo
   INQUIRE (IOLENGTH=LRecIn) TopoIn(:,1)
   OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
         UNIT=nfclm, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO j=1,Jdim
      READ (UNIT=nfclm, REC=j) TopoIn(:,j)
   END DO
   CLOSE (UNIT=nfclm)
   TopoInput=REAL(TopoIn,r8)
   VarTopInput=TopoInput*TopoInput
  
   ! Interpolate Input Regular Grid Albedo To Gaussian Grid on Output
   IF (Linear) THEN
      CALL DoLinearInterpolation (TopoInput, TopoOutput)
      CALL DoLinearInterpolation (VarTopInput, VarTopOutput)
   ELSE
      CALL DoAreaInterpolation (TopoInput, TopoOutput)
      CALL DoAreaInterpolation (VarTopInput, VarTopOutput)
   END IF
   VarTopOutput=MAX(VarTopOutput-TopoOutput*TopoOutput,VarMin)
   !VarTopOutput=MIN(VarTopOutput,VarMax)

   INQUIRE (IOLENGTH=LRecOut) FieldOut

   ! Write Out Adjusted Interpolated Topography Variance (m2)
   OPEN (FILE=TRIM(DirMain)//DirModelIn//VarNameV//nLats, &
         UNIT=nfoua, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//VarNameV//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   FieldOut=REAL(VarTopOutput,r4)
   WRITE (UNIT=nfoua, REC=1) FieldOut
   CLOSE (UNIT=nfoua)

   ! Write Out Adjusted Interpolated Topography (m)
   OPEN (FILE=TRIM(DirMain)//DirPreOut//VarNameT//nLats, &
         UNIT=nfoub, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarNameT//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   FieldOut=REAL(TopoOutput,r4)
   WRITE (UNIT=nfoub, REC=1) FieldOut
   CLOSE (UNIT=nfoub)
  
   IF (GrADS) THEN

      ! Write Out Adjusted Interpolated Topography and Variance
      OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//nLats, &
            UNIT=nfouc, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//nLats, &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      FieldOut=REAL(TopoOutput,r4)
      WRITE (UNIT=nfouc, REC=1) FieldOut
      FieldOut=REAL(VarTopOutput,r4)
      WRITE (UNIT=nfouc, REC=2) FieldOut
      CLOSE (UNIT=nfouc)

      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
            TRIM(DirMain)//DirPreOut//TRIM(VarNameG)//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Topography and Variance on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS 2'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO 0 99 Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'VART 0 99 Variance of Topography [m2]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)

   END IF
PRINT *, "*** VarTopo ENDS NORMALLY ***"

END PROGRAM VarTopo
