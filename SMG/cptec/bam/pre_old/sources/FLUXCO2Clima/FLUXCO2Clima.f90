!
!  $Author: tomita $
!  $Date: 2007/08/01 20:09:58 $
!  $Revision: 1.1.1.1 $
!
PROGRAM SSTClima

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nficn, nflsm, &
                               nfclm,nfclm2, nfsto, nfout, nfctl, &
                               Mend, Kmax, Imax, Jmax, Idim, Jdim, Mnwv2, &
                               Year, Month, Day, Hour, &
                               To, CO2SeaIce,SSTSeaIce, Lat0, &
                               CO2OpenWater, CO2SeaIceThreshold, LapseRate, &
                               SSTOpenWater, SSTSeaIceThreshold, &
                               Trunc, nLats, mskfmt, VarName, NameLSM, &
                               DirPreOut, DirModelIn, DirClmFluxCO2, FileClmFluxCO2, &
			       DirClmSST,FileClmSST, &
                               GrADS, Linear, &
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

   ! Reads the 1x1 CO2 Global Monthly OI Climatology From NCEP,
   ! Interpolates it Using Area Weigth or Bi-Linear Into a Gaussian Grid

   INTEGER :: j, i, m, nr, LRec, ios

   INTEGER :: ForecastDay

   REAL (KIND=r4) :: TimeOfDay

   REAL (KIND=r8) :: RGCO2Max, RGCO2Min, GGCO2Max, GGCO2Min, MGCO2Max, MGCO2Min
   REAL (KIND=r8) :: RGSSTMax, RGSSTMin
   INTEGER :: ICnDate(4), CurrentDate(4),Header(5)

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask, SeaIceMask

   REAL (KIND=r4), DIMENSION (:), ALLOCATABLE :: CoefTopIn

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: WrOut,WrIn

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: CoefTop

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Topog, CO2In,SSTIn, CO2Gaus,SSTGaus,&
                   SeaIceFlagIn, SeaIceFlagOut

   CHARACTER (LEN=6), DIMENSION (:), ALLOCATABLE :: CO2Label
   INTEGER, PARAMETER :: lmon(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
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
   ALLOCATE (CO2In(Idim,Jdim),SSTIn(Idim,Jdim),WrIn(Idim,Jdim), CO2Gaus(Imax,Jmax),SSTGaus(Imax,Jmax), WrOut(Imax,Jmax))
   ALLOCATE (SeaIceFlagIn(Idim,Jdim), SeaIceFlagOut(Imax,Jmax))
   ALLOCATE (CO2Label(Jdim))

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

   ! Open File for Land-Sea Mask and CO2 Data to Global Model Input
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
   ! Write out Land-Sea Mask to CO2 Data Set
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

   OPEN (UNIT=nfclm2, FILE=TRIM(DirMain)//TRIM(DirClmSST)//FileClmSST, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirClmSST)//FileClmSST, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF


   INQUIRE(IOLENGTH=lrec)WrIn
   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//TRIM(DirClmFluxCO2)//FileClmFluxCO2, &
         FORM='UNFORMATTED', ACCESS='DIRECT',RECL=lrec, &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirClmFluxCO2)//FileClmFluxCO2, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF

   ! Loop Through Months
   DO m=1,12
      READ (UNIT=nfclm2, FMT='(8I5)') Header
      WRITE (UNIT=nfprt, FMT='(/,1X,9I5,/)') m, Header
      READ (UNIT=nfclm2, FMT='(16F5.2)') SSTIn

      READ (nfclm,rec=m) WrIn
      CO2In=WrIn
      IF (MAXVAL(SSTIn) < 100.0_r8) SSTIn=SSTIn+To

      ! Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
      ! Over Open Water Set Input CO2 = MIN of CO2OpenWater
      ! Over Non Ice Points Before Interpolation
      IF (SSTSeaIce < 100.0_r8) SSTSeaIce=SSTSeaIce+To
      DO j=1,Jdim
         DO i=1,Idim
            SeaIceFlagIn(i,j)=0.0_r8
            IF (SSTIn(i,j) < SSTSeaIce) THEN
               SeaIceFlagIn(i,j)=1.0_r8
            ELSE
               SSTIn(i,j)=MAX(SSTIn(i,j),SSTOpenWater)
            END IF
         END DO
      END DO
      ! Min And Max Values of Input CO2
      RGSSTMax=MAXVAL(SSTIn)
      RGSSTMin=MINVAL(SSTIn)

      ! Interpolate Flag from 1x1 Grid to Gaussian Grid, Fill SeaIceMask=1
      ! Over Interpolated Points With 50% or More Sea Ice, =0 Otherwise
      IF (Linear) THEN
         CALL DoLinearInterpolation (SeaIceFlagIn, SeaIceFlagOut)
      ELSE
         CALL DoAreaInterpolation (SeaIceFlagIn, SeaIceFlagOut)
      END IF
 
      SeaIceMask=INT(SeaIceFlagOut+0.5_r8)
      WHERE (LandSeaMask == 1) SeaIceMask=0
 
      ! Interpolate SST from 1x1 Grid to Gaussian Grid
      IF (Linear) THEN
         CALL DoLinearInterpolation (SSTIn, SSTGaus)
      ELSE
         CALL DoAreaInterpolation (SSTIn, SSTGaus)
      END IF

      ! Interpolate CO2 from 1x1 Grid to Gaussian Grid
      IF (Linear) THEN
         CALL DoLinearInterpolation (CO2In, CO2Gaus)
      ELSE
         CALL DoAreaInterpolation (CO2In, CO2Gaus)
      END IF
      ! Min and Max Values of Gaussian Grid
      GGCO2Max=MAXVAL(CO2Gaus)
      GGCO2Min=MINVAL(CO2Gaus)

      DO j=1,Jmax
         DO i=1,Imax
            IF (LandSeaMask(i,j) == 1) THEN
               ! Set SST = Undef Over Land
               SSTGaus(i,j)=Undef
            ELSE IF (SeaIceMask(i,j) == 1) THEN
               ! Set SST Sea Ice Threshold Minus 1 Over Sea Ice
               SSTGaus(i,j)=SSTSeaIceThreshold-1.0_r8
            ELSE
               ! Correct SST for Topography, Do Not Create or
               ! Destroy Sea Ice Via Topography Correction
               SSTGaus(i,j)=SSTGaus(i,j)-Topog(i,j)*LapseRate
               IF (SSTGaus(i,j) < SSTSeaIceThreshold) &
                  SSTGaus(i,j)=SSTSeaIceThreshold+0.2_r8
            END IF
         END DO
      END DO

      DO j=1,Jmax
         DO i=1,Imax
            IF (LandSeaMask(i,j) == 1) THEN
               ! Set CO2 = Undef Over Land
               CO2Gaus(i,j)=Undef
            ELSE IF (SeaIceMask(i,j) == 1) THEN
               ! Set CO2 Sea Ice Threshold Minus 1 Over Sea Ice
               CO2Gaus(i,j)=CO2SeaIceThreshold
            ELSE
               ! Correct CO2 for Topography, Do Not Create or
               ! Destroy Sea Ice Via Topography Correction
               CO2Gaus(i,j)=(CO2Gaus(i,j)+0.03333334_r8 +0.00136_r8)/(86400.0_r8*REAL(lmon(m),kind=r8))!-Topog(i,j)*LapseRate
               !IF (CO2Gaus(i,j) < CO2SeaIceThreshold) &
               !   CO2Gaus(i,j)=CO2SeaIceThreshold+0.2_r8
            END IF
         END DO
      END DO
      ! Min and Max Values of Corrected Gaussian Grid CO2 Excluding Land Points
      MGCO2Max=MAXVAL(CO2Gaus,MASK=CO2Gaus/=Undef)
      MGCO2Min=MINVAL(CO2Gaus,MASK=CO2Gaus/=Undef)

      ! Write out Gaussian Grid Weekly CO2
      WrOut=REAL(CO2Gaus,r4)
      WRITE (UNIT=nfsto, REC=m+1) WrOut

      WRITE (UNIT=nfprt, FMT='(/,3(A,I2.2),A,I4)') &
            ' Hour = ', Hour, ' Day = ', Day, &
            ' Month = ', Month, ' Year = ', Year

      WRITE (UNIT=nfprt, FMT='(/,A,3(A,2F8.2,/))') &
            ' Mean Weekly CO2 Interpolation :', &
            ' Regular  Grid SST: Min, Max = ', RGSSTMin, RGSSTMax, &
            ' Gaussian Grid CO2: Min, Max = ', GGCO2Min, GGCO2Max, &
            ' Masked G Grid CO2: Min, Max = ', MGCO2Min, MGCO2Max

      IF (GrADS) THEN
         nr=1+5*(m-1)
         WrOut=REAL(Topog,r4)
         WRITE (UNIT=nfout, REC=nr) WrOut
         WrOut=REAL(1-2*LandSeaMask,r4)
         WRITE (UNIT=nfout, REC=nr+1) WrOut
         WrOut=REAL(SeaIceMask,r4)
         WRITE (UNIT=nfout, REC=nr+2) WrOut
         WrOut=REAL(SSTGaus,r4)
         WRITE (UNIT=nfout, REC=nr+3) WrOut
         WrOut=REAL(CO2Gaus,r4)
         WRITE (UNIT=nfout, REC=nr+4) WrOut
      END IF

   ! End Loop Through Months
   END DO
   CLOSE (UNIT=nfclm)
   CLOSE (UNIT=nfsto)
   CLOSE (UNIT=nfout)

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
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Monthly Climatological OI CO2 on a Gaussian Grid'
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF 12 LINEAR JAN2007 1MO'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  5'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SSTC  0 99 Climatological SST Topography Corrected [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'CO2C  0 99 Climatological CO2 Topography Corrected [kg/m2/s]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'

      CLOSE (UNIT=nfctl)
   END IF
PRINT *, "*** SSTClima ENDS NORMALLY ***"

END PROGRAM SSTClima
