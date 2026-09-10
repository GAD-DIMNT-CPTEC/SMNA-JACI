PROGRAM TopographyGradient

  USE InputParameters, ONLY: r4, r8, rad, dLon, EMRad1, &
                             nferr, nfprt, nfinp, nfout, nfctl, &
                             Imax, Jmax, Kmax,Idim,Jdim, Mend1, Mend2, &
                             Undef, ForecastDay, TimeOfDay, Months, &
                             TGrADS, FileInp, FileOut,DirMain,DirModelIn,DirPreOut, &
                             InitParameters,nfclm,nfoub,Linear,VarNameT,VarName,nLats,GrADS
   USE LinearInterpolation, ONLY: gLatsL=>LatOut, &
       InitLinearInterpolation, DoLinearInterpolation

   USE AreaInterpolation, ONLY: gLatsA=>gLats, &
       InitAreaInterpolation, DoAreaInterpolation

  USE InputArrays, ONLY: qTopoInp, GrADSOut, qtopo, qTopoS, Topo, &
                         DTopoDx, DTopoDy, CosLatInv, &
                         InitialDate, CurrentDate, SigmaInteface, SigmaLayer, &
                         GetArrays, ClsArrays

  USE FastFourierTransform, ONLY : CreateFFT

  USE LegendreTransform, ONLY : CreateGaussRep, CreateSpectralRep, &
                                CreateLegTrans, gLats, snnp1

  USE Spectral2Grid, ONLY : Transpose, SpecCoef2Grid, SpecCoef2GridD
 
  IMPLICIT NONE

  INTEGER :: i, j,ll, mm, nn, LRec,LRecIn,LRecOut,ios
  INTEGER, PARAMETER :: nVar=14

  REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: HPRIMEIn
  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: HPRIMEInput
  REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: HPRIMEOutput
  REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: FieldOut

  CALL InitParameters ()

  IF (Linear) THEN
     CALL InitLinearInterpolation ()
  ELSE
     CALL InitAreaInterpolation ()
  END IF

  ALLOCATE (HPRIMEIn(Idim,Jdim))
  ALLOCATE (HPRIMEInput(Idim,Jdim))
  ALLOCATE (HPRIMEOutput(Imax,Jmax))
  ALLOCATE (FieldOut(Imax,Jmax))

  CALL GetArrays ()

  CALL CreateSpectralRep ()
  CALL CreateGaussRep ()
  CALL CreateFFT ()
  CALL CreateLegTrans ()
  
   ! Read In Input Topo DATA
   INQUIRE (IOLENGTH=LRecIn) HPRIMEIn(:,:)
   OPEN (FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.dat', &
         UNIT=nfclm, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarName)//'.dat', &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   ! WRITE In Input Topo DATA
   INQUIRE (IOLENGTH=LRecOut) FieldOut(:,:)
   OPEN (FILE=TRIM(DirMain)//TRIM(DirModelIn)//TRIM(VarNameT)//nLats, &
         UNIT=nfoub, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, &
         ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirModelIn)//TRIM(VarNameT)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   DO j=1,nVar
      READ (UNIT=nfclm, REC=j) HPRIMEIn(:,:)
      HPRIMEInput=REAL(HPRIMEIn,r8)
      ! Interpolate Input Regular Grid Albedo To Gaussian Grid on Output
      IF (Linear) THEN
         CALL DoLinearInterpolation (HPRIMEInput(:,:), HPRIMEOutput(:,:))
      ELSE
         CALL DoAreaInterpolation   (HPRIMEInput(:,:), HPRIMEOutput(:,:))
      END IF
      FieldOut=REAL(HPRIMEOutput,r4)
      WRITE (UNIT=nfoub, REC=j) FieldOut
   END DO

   CLOSE (UNIT=nfclm)
   CLOSE (UNIT=nfoub)
   
   IF (GrADS) THEN
      ! Write GrADS Control File
      OPEN (FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//nLats//'.ctl', &
            UNIT=nfctl, FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//TRIM(DirPreOut)//TRIM(VarNameT)//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '//TRIM(DirMain)//TRIM(DirModelIn)//TRIM(VarNameT)//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Topography and Variance on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)')'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF 1 LEVELS 1000'
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      WRITE (UNIT=nfctl, FMT='(A,i5)')'VARS ',nVar
      WRITE (UNIT=nfctl, FMT='(A)')'HSTDV  0 99 standard deviation of orography'
      WRITE (UNIT=nfctl, FMT='(A)')'HCNVX  0 99 Normalized convexity'
      WRITE (UNIT=nfctl, FMT='(A)')'HASYW  0 99 orographic asymmetry in W-E plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HASYS  0 99 orographic asymmetry in S-N plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HASYSW 0 99 orographic asymmetry in SW-NE plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HASYNW 0 99 orographic asymmetry in NW-SE plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HLENW  0 99 orographic length scale in W-E plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HLENS  0 99 orographic length scale in S-N plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HLENSW 0 99 orographic length scale in SW-NE plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HLENNW 0 99 orographic length scale in NW-SE plane'
      WRITE (UNIT=nfctl, FMT='(A)')'HANGL  0 99 angle of the mountain range w/r/t east'
      WRITE (UNIT=nfctl, FMT='(A)')'HSLOP  0 99 slope of orography'
      WRITE (UNIT=nfctl, FMT='(A)')'HANIS  0 99 anisotropy/aspect ratio'
      WRITE (UNIT=nfctl, FMT='(A)')'HZMAX  0 99 max height above mean orography'
      WRITE (UNIT=nfctl, FMT='(A)')'ENDVARS'
      CLOSE (UNIT=nfctl)

   END IF

  
   ! Write Out Adjusted Interpolated Topography data 

  CosLatInV=1.0_r8/COS(rad*gLats)
  WRITE(*,*)TRIM(DirMain)//TRIm(DirModelIn)//TRIM(FileInp)
  WRITE(*,*)TRIM(DirMain)//TRIM(DirModelIn)//TRIM(FileOut)
  OPEN(UNIT=nfinp, FILE=TRIM(DirMain)//TRIm(DirModelIn)//TRIM(FileInp), &
                   FORM='UNFORMATTED', STATUS='OLD')

  INQUIRE (IOLENGTH=LRec) GrADSOut
  OPEN(UNIT=nfout, FILE=TRIM(DirMain)//TRIM(DirModelIn)//TRIM(FileOut), &
                   FORM='UNFORMATTED', STATUS='REPLACE', &
                   ACCESS='DIRECT', RECL=LRec)

  ! Read Topography Coeficients From Model Spectral Initial Condition
  
  READ (UNIT=nfinp) ForecastDay, TimeOfDay, InitialDate, CurrentDate, &
                    SigmaInteface, SigmaLayer
  write (*,*) ForecastDay, TimeOfDay
  write (*,*) InitialDate
  write (*,*) CurrentDate
  READ (UNIT=nfinp) qTopoInp
  qTopo=REAL(qTopoInp,r8)
  CALL Transpose (qTopo)
  CALL SpecCoef2Grid (qTopo, Topo)
  DO i=1,Imax
     GrADSOut(i,Jmax:1:-1)=REAL(Topo(i,1:Jmax),r4)
  END DO
  write (*,*) minval(gradsout)
  write (*,*) maxval(gradsout)
  WRITE (UNIT=nfout, REC=1) GrADSOut

  ! Zonal Gradient off Topography
  ll=0
  DO mm=1,Mend1
     DO nn=1,Mend2-mm
        ll=ll+1
        qTopoS(2*ll-1)=-REAL(nn-1,r8)*qTopo(2*ll)
        qTopoS(2*ll  )=+REAL(nn-1,r8)*qTopo(2*ll-1)
     END DO
  END DO
  CALL SpecCoef2Grid (qTopoS, DTopoDx)
  DO i=1,Imax
     GrADSOut(i,Jmax:1:-1)=REAL(DTopoDx(i,1:Jmax)*CosLatInV(1:Jmax)*EMRad1,r4)
  END DO
  WRITE (UNIT=nfout, REC=2) GrADSOut

 ! Meridional Gradient off Topography
   CALL SpecCoef2GridD (qTopo, DTopoDy)
   DO i=1,Imax
      GrADSOut(i,Jmax:1:-1)=REAL(DTopoDy(i,1:Jmax)*CosLatInV(1:Jmax)*EMRad1,r4)
   END DO
   WRITE (UNIT=nfout, REC=3) GrADSOut

  CLOSE(UNIT=nfinp)
  CLOSE(UNIT=nfout)

  TGrADS='  Z         '
  WRITE (TGrADS(1:2), FMT='(I2.2)') InitialDate(1)
  WRITE (TGrADS(4:5), FMT='(I2.2)') InitialDate(3)
  TGrADS(6:8)=Months(InitialDate(2))
  WRITE (TGrADS(9:12), FMT='(I4.4)') InitialDate(4)

  OPEN(UNIT=nfctl, FILE=TRIM(DirMain)//TRIM(DirPreOut)//TRIM(FileOut)//'.ctl', &
                   FORM='FORMATTED', STATUS='REPLACE')
  WRITE (UNIT=nfctl, FMT='(A)') 'dset '//TRIM(DirMain)//TRIM(DirModelIn)//TRIM(FileOut)
  WRITE (UNIT=nfctl, FMT='(A)') 'options big_endian'
  WRITE (UNIT=nfctl, FMT='(A,1P,G15.7)') 'undef ', Undef
  WRITE (UNIT=nfctl, FMT='(A)') 'title Topography and its Gradient'
  WRITE (UNIT=nfctl, FMT='(A,I5,A,F10.5)') 'xdef ', Imax, ' linear    0.0 ', dLon
  WRITE (UNIT=nfctl, FMT='(A,I5,A,5F10.5)') 'ydef ',Jmax, ' levels ', gLats(Jmax:Jmax-4:-1)
  WRITE (UNIT=nfctl, FMT='(18X,5F10.5)') gLats(Jmax-5:1:-1)
  WRITE (UNIT=nfctl, FMT='(A)') 'zdef 1 levels 1000'
  WRITE (UNIT=nfctl, FMT='(3A)') 'tdef 1 linear ', TGrADS, ' 06hr'
  WRITE (UNIT=nfctl, FMT='(A)') 'vars 3'
  WRITE (UNIT=nfctl, FMT='(A)') 'topo 1 99 Topography (m)'
  WRITE (UNIT=nfctl, FMT='(A)') 'dtpx 1 99 Zonal Gradient of Topography (m/m)'
  WRITE (UNIT=nfctl, FMT='(A)') 'dtpy 1 99 Meridional Gradient of Topography (m/m)'
  WRITE (UNIT=nfctl, FMT='(A)') 'endvars'
  CLOSE (UNIT=nfctl)

  CALL ClsArrays ()

END PROGRAM TopographyGradient
