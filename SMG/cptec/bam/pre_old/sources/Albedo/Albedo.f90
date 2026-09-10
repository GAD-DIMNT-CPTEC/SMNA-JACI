!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM Albedo

   USE InputParameters, ONLY: r4, r8, Undef, &
                              nferr, nfclm, nfout, nfctl, &
                              Imax, Jmax, Idim, Jdim, &
                              nLats, VarNameOut, DirPreOut, &
                              VarName, DirMain, GrADS, Linear, &
                              InitInputParameters

   USE LinearInterpolation, ONLY: gLatsL=>LatOut, &
       InitLinearInterpolation, DoLinearInterpolation

   USE AreaInterpolation, ONLY: gLatsA=>gLats, &
       InitAreaInterpolation, DoAreaInterpolation

   IMPLICIT NONE

   ! Linear (.TRUE.) or Area Weighted (.FALSE.) Interpolation
   ! Interpolate Regular To Gaussian
   ! Regular Input Data is Assumed to be Oriented with
   ! the North Pole and Greenwich as the First Point
   ! Gaussian Output Data is Interpolated to be Oriented with
   ! the North Pole and Greenwich as the First Point
   ! Input for the AGCM is Assumed IEEE 32 Bits Big Endian

   INTEGER :: Month, LRecIn, LRecOut, ios

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: AlbedoIn
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: AlbedoOut
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: AlbedoInput
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: AlbedoOutput

   CALL InitInputParameters ()

   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (AlbedoInput(Idim,Jdim))
   ALLOCATE (AlbedoOutput(Imax,Jmax))
   ALLOCATE (AlbedoIn(Idim,Jdim))
   ALLOCATE (AlbedoOut(Imax,Jmax))

   INQUIRE (IOLENGTH=LRecIn) AlbedoIn
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

   INQUIRE (IOLENGTH=LRecOut) AlbedoOut
   OPEN (FILE=TRIM(DirMain)//DirPreOut//VarNameOut//nLats, &
         UNIT=nfout, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarNameOut//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   ! Read In Input Albedo

   DO Month=1,12

      READ (UNIT=nfclm, REC=Month) AlbedoIn
      AlbedoInput=REAL(Albedoin,r8)
  
      ! Interpolate Input Regular Grid Albedo To Gaussian Grid on Output
  
      IF (Linear) THEN
         CALL DoLinearInterpolation (AlbedoInput, AlbedoOutput)
      ELSE
         CALL DoAreaInterpolation (AlbedoInput, AlbedoOutput)
      END IF
      AlbedoOut=REAL(AlbedoOutput,r4)
  
      ! Write Out Adjusted Interpolated Albedo
  
      WRITE (UNIT=nfout, REC=Month) AlbedoOut
  
   END DO

   CLOSE (UNIT=nfclm)
   CLOSE (UNIT=nfout)

   IF (GrADS) THEN
      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//DirPreOut//VarNameOut//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarNameOut//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//VarNameOut//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Albedo on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
            'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF  1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  1'
      WRITE (UNIT=nfctl, FMT='(A)') 'ALBE  0 99 Albedo [No Dim]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'
      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** Albedo ENDS NORMALLY ***"

END PROGRAM Albedo
