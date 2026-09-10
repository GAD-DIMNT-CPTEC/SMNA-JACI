!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM OCMClima

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nficn, nflsm, &
                               nfclm,nfocm, nfsto, nfout, nfctl, &
                               Mend, Kmax, Imax, Jmax, Idim, Jdim, Mnwv2, &
                               Year, Month, Day, Hour, nMon,nlevl,&
                               To, OCMSeaIce, Lat0, &
                               OCMOpenWater, OCMSeaIceThreshold, LapseRate, &
                               Trunc, nLats, mskfmt, VarName, NameLSM, &
                               DirPreOut, DirModelIn, DirClmOCM, FileClmOCM, &
                               GrADS, Linear,FileClmOCMData, &
                               MonthChar, DateICn, Preffix, Suffix, DirMain, &
                               InitInputParameters

   USE FastFourierTransform, ONLY : CreateFFT

   USE LegendreTransform, ONLY : CreateGaussRep, CreateSpectralRep, CreateLegTrans

   USE SpectralGrid, ONLY : SpecCoef2Grid

   USE LinearInterpolation, ONLY: gLatsL=>LatOut, &
       InitLinearInterpolation, DoLinearInterpolation

   USE AreaInterpolation, ONLY: gLatsA=>gLats, &
       InitAreaInterpolation, DoAreaInterpolation

   IMPLICIT NONE

   ! Reads the 1x1 OCM Global Monthly OI Climatology From NCEP,
   ! Interpolates it Using Area Weigth or Bi-Linear Into a Gaussian Grid

   INTEGER :: j, i, m, nr, LRec, ios,it,k,irec,im,ierr

   INTEGER :: ForecastDay

   REAL (KIND=r4) :: TimeOfDay

   REAL (KIND=r8) :: RGOCMMax, RGOCMMin, GGOCMMax, GGOCMMin, MGOCMMax, MGOCMMin

   INTEGER :: ICnDate(4), CurrentDate(4), Header(8)

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask, SeaIceMask

   REAL (KIND=r4), DIMENSION (:), ALLOCATABLE :: CoefTopIn

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: WrOut
   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: aux
   INTEGER       , DIMENSION (:,:), ALLOCATABLE :: auxI

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: CoefTop

   REAL (KIND=r8)   , ALLOCATABLE :: otemp_in  (:,:,:,:)
   REAL (KIND=r8)   , ALLOCATABLE :: salt_in   (:,:,:,:)
   REAL (KIND=r8)   , ALLOCATABLE :: bathy_in   (:,:)
   REAL (KIND=r8)   , ALLOCATABLE :: waterqual_in  (:,:)

   REAL (KIND=r8)   , ALLOCATABLE :: otemp_out (:,:,:,:)
   REAL (KIND=r8)   , ALLOCATABLE :: salt_out  (:,:,:,:)
   REAL (KIND=r8)   , ALLOCATABLE :: bathy_out  (:,:)
   REAL (KIND=r8)   , ALLOCATABLE :: waterqual_out  (:,:)

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Topog, OCMIn, OCMGaus, &
                   SeaIceFlagIn, SeaIceFlagOut

   CHARACTER (LEN=6), DIMENSION (:), ALLOCATABLE :: OCMLabel

   CALL InitInputParameters ()
   CALL CreateSpectralRep ()
   CALL CreateGaussRep ()
   CALL CreateFFT ()
   CALL CreateLegTrans ()

   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (LandSeaMask(Imax,Jmax), SeaIceMask(Imax,Jmax))
   ALLOCATE (CoefTopIn(Mnwv2))
   ALLOCATE (CoefTop(Mnwv2), Topog(Imax,Jmax))
   ALLOCATE (OCMIn(Idim,Jdim), OCMGaus(Imax,Jmax), WrOut(Imax,Jmax))
   ALLOCATE (SeaIceFlagIn(Idim,Jdim), SeaIceFlagOut(Imax,Jmax))
   ALLOCATE (OCMLabel(Jdim))
   ALLOCATE (aux (Idim,Jdim))
   ALLOCATE (otemp_in        (Idim,Jdim,nlevl,nMon))
   ALLOCATE (salt_in         (Idim,Jdim,nlevl,nMon))
   ALLOCATE (bathy_in        (Idim,Jdim))
   ALLOCATE (waterqual_in    (Idim,Jdim))

   ALLOCATE(otemp_out       (Imax,Jmax,nlevl,nMon))
   ALLOCATE(salt_out        (Imax,Jmax,nlevl,nMon))
   ALLOCATE(bathy_out       (Imax,Jmax))
   ALLOCATE(waterqual_out   (Imax,Jmax))
   ALLOCATE(auxI  (Imax,Jmax))
   INQUIRE(IOLENGTH=LRec)aux(1:Idim,1:Jdim)

   OPEN(nfocm,FILE=TRIM(DirMain)//TRIM(DirClmOCM)//TRIM(FileClmOCMData), &
         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LRec,STATUS='OLD',&
         ACTION='READ', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//TRIM(DirClmOCM)//TRIM(FileClmOCMData), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF


   irec=0
   DO it=1, nMon
      DO k=1,nlevl
         irec=irec+1
         READ(nfocm,rec=irec)aux
         otemp_in(1:Idim,1:Jdim,k,it)=aux(1:Idim,1:Jdim)
      END DO
      DO k=1,nlevl
         irec=irec+1
         READ(nfocm,rec=irec)aux!
         salt_in(1:Idim,1:Jdim,k,it)=aux(1:Idim,1:Jdim)
      END DO
      irec=irec+1
      READ(nfocm,rec=irec)aux
       bathy_in=aux
      irec=irec+1
      READ(nfocm,rec=irec)aux
       waterqual_in=aux
   END DO
   irec=0
   ! Read in SpectraL Coefficient of Topography from ICn
   ! to Ensure that Topography is the Same as Used by Model
   OPEN (FILE=TRIM(DirMain)//DirModelIn//Preffix//DateICn//Suffix//Trunc, &
         UNIT=nficn, FORM='UNFORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//Preffix//DateICn//Suffix//Trunc, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nficn) ForecastDay, TimeOfDay, ICnDate, CurrentDate
   READ  (UNIT=nficn) CoefTopIn
   CLOSE (UNIT=nficn)
   CoefTop=REAL(CoefTopIn,r8)
   CALL SpecCoef2Grid (1, CoefTop, Topog)

   ! Read in Land-Sea Mask Data Set
   OPEN (UNIT=nflsm, FILE=TRIM(DirMain)//DirPreOut//NameLSM//nLats, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//NameLSM//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nflsm, FMT=mskfmt) LandSeaMask
   CLOSE (UNIT=nflsm)

   ! Open File for Land-Sea Mask and OCM Data to Global Model Input
   INQUIRE (IOLENGTH=LRec) WrOut
   OPEN (UNIT=nfsto, FILE=TRIM(DirMain)//DirModelIn//VarName//DateICn(1:8)//nLats, &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//VarName//DateICn(1:8)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   ! Write out Land-Sea Mask to OCM Data Set
   ! The LSMask will be Transfered by Model to Post-Processing

   WrOut=REAL(1-2*LandSeaMask,r4)
   WRITE (UNIT=nfsto, REC=1) WrOut

   IF (GrADS) THEN
      OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats, &
            FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats, &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
   END IF

   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//TRIM(DirClmOCM)//FileClmOCM, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)

   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirClmOCM)//FileClmOCM, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   ! Loop Through Months
   irec=0
   DO m=1,12
      PRINT*,'month=',m
      READ (UNIT=nfclm, FMT='(8I5)') Header
      WRITE (UNIT=nfprt, FMT='(/,1X,9I5,/)') m, Header
      READ (UNIT=nfclm, FMT='(16F5.2)') OCMIn
      IF (MAXVAL(OCMIn) < 100.0_r8) OCMIn=OCMIn+To

      ! Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
      ! Over Open Water Set Input OCM = MIN of OCMOpenWater
      ! Over Non Ice Points Before Interpolation
      IF (OCMSeaIce < 100.0_r8) OCMSeaIce=OCMSeaIce+To
      DO j=1,Jdim
         DO i=1,Idim
            SeaIceFlagIn(i,j)=0.0_r8
            IF (OCMIn(i,j) < OCMSeaIce) THEN
               SeaIceFlagIn(i,j)=1.0_r8
            ELSE
               OCMIn(i,j)=MAX(OCMIn(i,j),OCMOpenWater)
            END IF
         END DO
      END DO
      ! Min And Max Values of Input OCM
      RGOCMMax=MAXVAL(OCMIn)
      RGOCMMin=MINVAL(OCMIn)

      ! Interpolate Flag from 1x1 Grid to Gaussian Grid, Fill SeaIceMask=1
      ! Over Interpolated Points With 50% or More Sea Ice, =0 Otherwise
      PRINT*,'SeaIceFlagOut='

      IF (Linear) THEN
         CALL DoLinearInterpolation (SeaIceFlagIn, SeaIceFlagOut)
      ELSE
         CALL DoAreaInterpolation (SeaIceFlagIn, SeaIceFlagOut)
      END IF
      SeaIceMask=INT(SeaIceFlagOut+0.5_r8)
      WHERE (LandSeaMask == 1) SeaIceMask=0
      PRINT*,'OCMIn='

      ! Interpolate OCM from 1x1 Grid to Gaussian Grid
      IF (Linear) THEN
         CALL DoLinearInterpolation (OCMIn, OCMGaus)
      ELSE
         CALL DoAreaInterpolation (OCMIn, OCMGaus)
      END IF
       PRINT*,'otemp_out='


       ! Interpolate OCM from 1x1 Grid to Gaussian Grid
      DO k=1,nlevl
         IF (Linear) THEN
            CALL DoLinearInterpolation (otemp_in(:,:,k,m), otemp_out(:,:,k,m))
         ELSE
            CALL DoAreaInterpolation   (otemp_in(:,:,k,m), otemp_out(:,:,k,m))
         END IF
      END DO
      PRINT*,'salt_in='

      ! Interpolate OCM from 1x1 Grid to Gaussian Grid
      DO k=1,nlevl
         IF (Linear) THEN
            CALL DoLinearInterpolation (salt_in(:,:,k,m), salt_out(:,:,k,m))
         ELSE
            CALL DoAreaInterpolation   (salt_in(:,:,k,m), salt_out(:,:,k,m))
         END IF
      END DO
      PRINT*,'bathy_in='

      ! Interpolate OCM from 1x1 Grid to Gaussian Grid
      IF (Linear) THEN
         CALL DoLinearInterpolation (bathy_in, bathy_out)
      ELSE
         CALL DoAreaInterpolation (bathy_in, bathy_out)
      END IF
      
      ! Interpolate OCM from 1x1 Grid to Gaussian Grid
      PRINT*,'waterqual_in='

      IF (Linear) THEN
         CALL DoLinearInterpolation (waterqual_in, waterqual_out)
      ELSE
         CALL DoAreaInterpolation (waterqual_in, waterqual_out)
      END IF

      ! Min and Max Values of Gaussian Grid
      GGOCMMax=MAXVAL(OCMGaus)
      GGOCMMin=MINVAL(OCMGaus)
      PRINT*,'Min and Max Values of Gaussian Grid='

      DO j=1,Jmax
         DO i=1,Imax
            IF (LandSeaMask(i,j) == 1) THEN
               ! Set OCM = Undef Over Land
               OCMGaus(i,j)=Undef
               waterqual_out(i,j) = 0
               IF(Topog(i,j) <=0.0)THEN
                  bathy_out(i,j)=1
               ELSE
                  bathy_out(i,j)=Topog(i,j)
               END IF
            ELSE IF (SeaIceMask(i,j) == 1) THEN
               ! Set OCM Sea Ice Threshold Minus 1 Over Sea Ice
               bathy_out(i,j)=bathy_out(i,j)
               OCMGaus(i,j)=OCMSeaIceThreshold-1.0_r8
               waterqual_out(i,j) = INT(waterqual_out(i,j))
            ELSE
               bathy_out(i,j)=bathy_out(i,j)
               waterqual_out(i,j) = INT(waterqual_out(i,j))
               ! Correct OCM for Topography, Do Not Create or
               ! Destroy Sea Ice Via Topography Correction
               OCMGaus(i,j)=OCMGaus(i,j)-Topog(i,j)*LapseRate
               IF (OCMGaus(i,j) < OCMSeaIceThreshold) &
                  OCMGaus(i,j)=OCMSeaIceThreshold+0.2_r8
            END IF
         END DO
      END DO
      PRINT*,'Min and Max Values of Corrected Gaussian Grid OCM Excluding Land Points'

      ! Min and Max Values of Corrected Gaussian Grid OCM Excluding Land Points
      MGOCMMax=MAXVAL(OCMGaus,MASK=OCMGaus/=Undef)
      MGOCMMin=MINVAL(OCMGaus,MASK=OCMGaus/=Undef)

      ! Write out Gaussian Grid Weekly OCM
      WrOut=REAL(OCMGaus,r4)
      WRITE (UNIT=nfsto, REC=m+1) WrOut

      WRITE (UNIT=nfprt, FMT='(/,3(A,I2.2),A,I4)') &
            ' Hour = ', Hour, ' Day = ', Day, &
            ' Month = ', Month, ' Year = ', Year

      WRITE (UNIT=nfprt, FMT='(/,A,3(A,2F8.2,/))') &
            ' Mean Weekly OCM Interpolation :', &
            ' Regular  Grid OCM: Min, Max = ', RGOCMMin, RGOCMMax, &
            ' Gaussian Grid OCM: Min, Max = ', GGOCMMin, GGOCMMax, &
            ' Masked G Grid OCM: Min, Max = ', MGOCMMin, MGOCMMax

      IF (GrADS) THEN
         !         nr=1+4*(m-1)
         WrOut=REAL(Topog,r4)
         irec=irec+1
         WRITE (UNIT=nfout, REC=irec) WrOut
         WrOut=REAL(1-2*LandSeaMask,r4)
         irec=irec+1
         WRITE (UNIT=nfout, REC=irec) WrOut
         WrOut=REAL(SeaIceMask,r4)
         irec=irec+1
         WRITE (UNIT=nfout, REC=irec) WrOut
         WrOut=REAL(OCMGaus,r4)
         irec=irec+1
         WRITE (UNIT=nfout, REC=irec) WrOut
         WrOut=REAL(bathy_out,r4)
         irec=irec+1
         WRITE (UNIT=nfout, REC=irec) WrOut
         WrOut=REAL(waterqual_out,r4)
         irec=irec+1
         WRITE (UNIT=nfout, REC=irec) WrOut
         DO k=1,nlevl
            WrOut=REAL(otemp_out(:,:,nlevl-k+1,m),r4)+To
            irec=irec+1
            WRITE(nfout,rec=irec) WrOut
         END DO
         DO k=1,nlevl
            WrOut=REAL(salt_out(:,:,nlevl-k+1,m),r4)
            irec=irec+1
            WRITE(nfout,rec=irec) WrOut
         END DO
      END IF

   ! End Loop Through Months
   END DO




   CLOSE (UNIT=nfclm)
   CLOSE (UNIT=nfsto)
   CLOSE (UNIT=nfout)
   ! Open File for Land-Sea Mask and OCM Data to Global Model Input
   !INQUIRE (IOLENGTH=LRec) WrOut
   !OPEN (UNIT=nfsto, FILE=TRIM(DirMain)//DirModelIn//VarName//DateICn(1:8)//nLats, &
   !      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec, ACTION='WRITE', &
   !      STATUS='REPLACE', IOSTAT=ios)
   !IF (ios /= 0) THEN
   !   WRITE (UNIT=nferr, FMT='(3A,I4)') &
   !         ' ** (Error) ** Open file ', &
   !           TRIM(DirMain)//DirModelIn//VarName//DateICn(1:8)//nLats, &
   !         ' returned IOStat = ', ios
   !   STOP  ' ** (Error) **'
   !END IF

    !
    !   read in ocean bathymetry (<0 over ocean; >=0 over land) 
    !   and determines the bottom of the model at each point, 
    !   lbottom indicates the first inactive layer
    INQUIRE (IOLENGTH=LRec) WrOut(1:Imax,1:Jmax)
    OPEN(92, file=TRIM(DirMain)//DirModelIn//'/'//'ocean_depth'//nLats, form='unformatted', &
         access='direct',recl=LRec, status='REPLACE', action='write', IOSTAT=ios)
    IF (ios /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(DirMain)//DirModelIn//'/'//'ocean_depth',ios
       STOP "**(ERROR)**"
    END IF
    irec=0
    irec=irec+1
    WrOut=REAL(bathy_out,r4)
    PRINT*,'aaa',MAXVAL(WrOut),MINVAL(WrOut)
    WRITE(92,rec=1)WrOut
    CLOSE (92)
    !   global annual average optical water type from the map of 
    !   siminot and le treut (1986,jgr).
    !   water types  -   numerical value in file:
    !    land               0
    !    i                  1
    !    ii                 2
    !    iii                3
    !    ia                 4
    !    ib                 5
    INQUIRE (IOLENGTH=LRec) WrOut(1:Imax,1:Jmax)
    OPEN(52, file=TRIM(DirMain)//DirModelIn//'/'//'water_type'//nLats, form='unformatted', &
         access='direct',recl=LRec, status='REPLACE', action='write', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(DirMain)//DirModelIn//'/'//'water_type', ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    irec=irec+1
    auxI=INT(waterqual_out)
    PRINT*,'aaa',MAXVAL(auxI),MINVAL(auxI)
    WRITE(52,rec=irec)auxI
    CLOSE (52)


!   observed monthly ltm mean temp and salinity 
!   will be used in relaxation              
    INQUIRE (IOLENGTH=LRec) WrOut(1:Imax,1:Jmax)
    OPEN(62, file=TRIM(DirMain)//DirModelIn//'/'//'temp_ltm_month'//nLats, form='unformatted', &
         access='direct',recl=LRec, status='REPLACE', action='write', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(DirMain)//DirModelIn//'/'//'temp_ltm_month.dat', ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO im=1,12
      DO k=1,nlevl
          irec=irec+1
          WrOut=REAL(otemp_out(:,:,nlevl-k+1,im),r4)+To
          WRITE(62,rec=irec)WrOut
       END DO
    END DO
    CLOSE(62)

!   observed monthly ltm mean temp and salinity 
!   will be used in relaxation              
    INQUIRE (IOLENGTH=LRec) WrOut(1:Imax,1:Jmax)
    OPEN(72, file=TRIM(DirMain)//DirModelIn//'/'//'salt_ltm_month'//nLats, form='unformatted', &
         access='direct',recl=LRec, status='REPLACE', action='write', IOSTAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(UNIT=nfprt,FMT="('**(ERROR)** Open file ',a,' returned iostat=',i4)") &
            TRIM(DirMain)//DirModelIn//'/'//'salt_ltm_month.dat', ierr
       STOP "**(ERROR)**"
    END IF
    irec=0
    DO im=1,12
      DO k=1,nlevl
          irec=irec+1
          WrOut=REAL(salt_out(:,:,nlevl-k+1,im),r4)
          WRITE(72,rec=irec)WrOut
       END DO
    END DO
    CLOSE(72)

   IF (GrADS) THEN
      ! Write GrADS Control File
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Monthly Climatological OI OCM on a Gaussian Grid'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',Imax,' LINEAR ',0.0_r8,360.0_r8/REAL(Imax,r8)
      WRITE (UNIT=nfctl, FMT='(A,I5,A)') 'YDEF ',Jmax,' LEVELS '
      IF (Linear) THEN
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsL(Jmax:1:-1)
      ELSE
         WRITE (UNIT=nfctl, FMT='(8F10.5)') gLatsA(Jmax:1:-1)
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'ZDEF    19 LEVELS  1000 900 800 700 600 500'
      WRITE (UNIT=nfctl, FMT='(A)') '                    400 300 250 200 150 125'
      WRITE (UNIT=nfctl, FMT='(A)') '                    100  75  50  30  20  10 0'

      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2007 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  8'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      WRITE (UNIT=nfctl, FMT='(A)') 'OCMC  0 99 Climatological OCM Topography Corrected [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'BATM  0 99 bathymetric [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'WAQL  0 99 waterqual   [0-5]'
      WRITE (UNIT=nfctl, FMT='(A)') 'TEMP 19 99 otemp       [C]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SALT 19 99 salinity    [g/kg]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'

      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** OCMClima ENDS NORMALLY ***"

END PROGRAM OCMClima
