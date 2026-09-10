PROGRAM SSTWeekly

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nficn, nflsm, nfsti, &
                               nfclm, nfsto, nfout, nfctl, &
                               Mend, Kmax, Imax, Jmax, Idim, Jdim, Mnwv2, &
                               Year, Month, Day, Hour, MonthLength, &
                               GrADSTime, Trunc, nLats, mskfmt, VarName, NameLSM, &
                               DirPreOut, DirModelIn, DirClmSST, FileClmSST, &
                               To, SSTSeaIce, LatClimSouth, LatClimNorth, Lat0, &
                               SSTOpenWater, SSTSeaIceThreshold, LapseRate, &
                               ClimWindow, GrADS, Linear, &
                               DateICn, Preffix, Suffix, DirMain, &
                               InitInputParameters

   USE FastFourierTransform, ONLY : CreateFFT

   USE LegendreTransform, ONLY : CreateGaussRep, CreateSpectralRep, CreateLegTrans

   USE SpectralGrid, ONLY : SpecCoef2Grid

   USE LinearInterpolation, ONLY: gLatsL=>LatOut, &
       InitLinearInterpolation, DoLinearInterpolation

   USE AreaInterpolation, ONLY: gLatsA=>gLats, &
       InitAreaInterpolation, DoAreaInterpolation

   IMPLICIT NONE

   ! Reads the Mean Weekly 1x1 SST Global From NCEP,
   ! Interpolates it Using Area Weigth Into a Gaussian Grid

   INTEGER :: j, i, m, js, jn, js1, jn1, ja, jb, LRecIn, LRecOut, ios

   INTEGER :: ForecastDay

   REAL (KIND=r4) :: TimeOfDay

   REAL (KIND=r8) :: RGSSTMax, RGSSTMin, GGSSTMax, GGSSTMin, MGSSTMax, MGSSTMin

   INTEGER :: ICnDate(4), CurrentDate(4)

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask, SeaIceMask

   REAL (KIND=r4), DIMENSION (:), ALLOCATABLE :: CoefTopIn

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SSTWklIn, WrOut

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: CoefTop

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Topog, SSTClim, &
                   SSTIn, SSTGaus, SeaIceFlagIn, SeaIceFlagOut

   CHARACTER (LEN=6), DIMENSION (:), ALLOCATABLE :: SSTLabel

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
   ALLOCATE (CoefTopIn(Mnwv2), SSTWklIn(Idim,Jdim))
   ALLOCATE (CoefTop(Mnwv2), Topog(Imax,Jmax), SSTClim(Idim,Jdim))
   ALLOCATE (SSTIn(Idim,Jdim), SSTGaus(Imax,Jmax), WrOut(Imax,Jmax))
   ALLOCATE (SeaIceFlagIn(Idim,Jdim), SeaIceFlagOut(Imax,Jmax))
   ALLOCATE (SSTLabel(Jdim))

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

   ! Read Mean Weekly 1 deg x 1 deg SST
   INQUIRE (IOLENGTH=LRecIn) SSTWklIn
   OPEN (UNIT=nfsti, FILE=TRIM(DirMain)//DirPreOut//VarName//'.'//DateICn(1:8), &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecIn, ACTION='READ', &
         STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirPreOut//VarName//'.'//DateICn(1:8), &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   READ  (UNIT=nfsti, REC=1) SSTWklIn
   CLOSE (UNIT=nfsti)
   IF (MAXVAL(SSTWklIn) < 100.0_r8) SSTWklIn=SSTWklIn+To

   ! Get SSTClim Climatological and Index for High Latitude
   ! Substitution of SST Actual by Climatology
   IF (ClimWindow) THEN
      CALL SSTClimatological ()
      IF (MAXVAL(SSTClim) < 100.0_r8) SSTClim=SSTClim+To
      CALL SSTClimaWindow ()
      jn1=jn-1
      js1=js+1
      ja=jn
      jb=js
   ELSE
      jn=0
      js=Jdim+1
      jn1=0
      js1=Jdim+1
      ja=1
      jb=Jdim
   END IF

   DO j=1,Jdim
      IF (j >= jn .AND. j <= js) THEN
         SSTLabel(j)='Weekly'
         DO i=1,Idim
            SSTIn(i,j)=REAL(SSTWklIn(i,j),r8)
         END DO
      ELSE
         SSTLabel(j)='Climat'
         DO i=1,Idim
            SSTIn(i,j)=SSTClim(i,j)
         END DO
      END IF
   END DO
   IF (jn1 >= 1) WRITE (UNIT=nfprt, FMT='(6(I4,1X,A))') (j,SSTLabel(j),j=1,jn1)
   WRITE (UNIT=nfprt, FMT='(6(I4,1X,A))') (j,SSTLabel(j),j=ja,jb)
   IF (js1 <= Jdim) WRITE (UNIT=nfprt, FMT='(6(I4,1X,A))') (j,SSTLabel(j),j=js1,Jdim)

   ! Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
   ! Over Open Water Set Input SST = MIN of SSTOpenWater
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
   ! Min And Max Values of Input SST
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
   ! Min and Max Values of Gaussian Grid
   GGSSTMax=MAXVAL(SSTGaus)
   GGSSTMin=MINVAL(SSTGaus)

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
   ! Min and Max Values of Corrected Gaussian Grid SST Excluding Land Points
   MGSSTMax=MAXVAL(SSTGaus,MASK=SSTGaus/=Undef)
   MGSSTMin=MINVAL(SSTGaus,MASK=SSTGaus/=Undef)

   ! Write out Land-Sea Mask and SST Data to Global Model Input
   INQUIRE (IOLENGTH=LRecOut) WrOut
   OPEN (UNIT=nfsto, FILE=TRIM(DirMain)//DirModelIn//VarName//DateICn(1:8)//nLats, &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//VarName//DateICn(1:8)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   ! Write out Land-Sea Mask to SST Data Set
   ! The LSMask will be Transfered by Model to Post-Processing
   WrOut=REAL(1-2*LandSeaMask,r4)
   WRITE (UNIT=nfsto, REC=1) WrOut
   ! Write out Gaussian Grid Weekly SST
   WrOut=REAL(SSTGaus,r4)
   WRITE (UNIT=nfsto, REC=2) WrOut
   CLOSE (UNIT=nfsto)

   WRITE (UNIT=nfprt, FMT='(/,3(A,I2.2),A,I4)') &
         ' Hour = ', Hour, ' Day = ', Day, &
         ' Month = ', Month, ' Year = ', Year

   WRITE (UNIT=nfprt, FMT='(/,A,3(A,2F8.2,/))') &
         ' Mean Weekly SST Interpolation :', &
         ' Regular  Grid SST: Min, Max = ', RGSSTMin, RGSSTMax, &
         ' Gaussian Grid SST: Min, Max = ', GGSSTMin, GGSSTMax, &
         ' Masked G Grid SST: Min, Max = ', MGSSTMin, MGSSTMax

   IF (GrADS) THEN
      
      OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats, &
            FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
            STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//VarName//DateICn(1:8)//nLats, &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WrOut=REAL(Topog,r4)
      WRITE (UNIT=nfout, REC=1) WrOut
      WrOut=REAL(1-2*LandSeaMask,r4)
      WRITE (UNIT=nfout, REC=2) WrOut
      WrOut=REAL(SeaIceMask,r4)
      WRITE (UNIT=nfout, REC=3) WrOut
      WrOut=REAL(SSTGaus,r4)
      WRITE (UNIT=nfout, REC=4) WrOut
      CLOSE (UNIT=nfout)
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Weekly SST on a Gaussian Grid'
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
      WRITE (UNIT=nfctl, FMT='(A)') 'TDEF  1 LINEAR '//GrADSTime//' 6HR'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  4'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SSTW  0 99 Weekly SST Topography Corrected [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'

      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** SSTWeekly ENDS NORMALLY ***"

CONTAINS


SUBROUTINE SSTClimatological ()

   IMPLICIT NONE

   ! 1950-1979 1 Degree x 1 Degree SST 
   ! Global NCEP OI Monthly Climatology
   ! Grid Orientation (SSTR):
   ! (1,1) = (0.5_r8W,89.5_r8N)
   ! (Idim,Jdim) = (0.5_r8E,89.5_r8S)

   INTEGER :: m, MonthBefore, MonthAfter

   REAL (KIND=r8) :: DayHour, DayCorrection, FactorBefore, FactorAfter

   INTEGER :: Header(8)

   REAL (KIND=r8) :: SSTBefore(Idim,Jdim), SSTAfter(Idim,Jdim)

   DayHour=REAL(Day,r8)+REAL(Hour,r8)/24.0_r8
   MonthBefore=Month-1
   IF (DayHour > (1.0_r8+REAL(MonthLength(Month),r8)/2.0_r8)) &
       MonthBefore=Month
   MonthAfter=MonthBefore+1
   IF (MonthBefore < 1) MonthBefore=12
   IF (MonthAfter > 12) MonthAfter=1
   DayCorrection=REAL(MonthLength(MonthBefore),r8)/2.0_r8-1.0_r8
   IF (MonthBefore == Month) DayCorrection=-DayCorrection-2.0_r8
   FactorAfter=2.0_r8*(DayHour+DayCorrection)/ &
               REAL(MonthLength(MonthBefore)+MonthLength(MonthAfter),r8)
   FactorBefore=1.0_r8-FactorAfter

   WRITE (UNIT=nfprt, FMT='(/,A)') ' From SSTClimatological:'
   WRITE (UNIT=nfprt, FMT='(/,A,I4,3(A,I2.2))') &
         ' Year = ', Year, ' Month = ', Month, &
         ' Day = ', Day, ' Hour = ', Hour
   WRITE (UNIT=nfprt, FMT='(/,2(A,I2))') &
         ' MonthBefore = ', MonthBefore, ' MonthAfter = ', MonthAfter
   WRITE (UNIT=nfprt, FMT='(/,2(A,F9.6),/)') &
         ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter

   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//TRIM(DirClmSST)//FileClmSST, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirClmSST)//FileClmSST, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO m=1,12
      READ (UNIT=nfclm, FMT='(8I5)') Header
      WRITE (UNIT=nfprt, FMT='(/,1X,9I5,/)') m, Header
      READ (UNIT=nfclm, FMT='(16F5.2)') SSTClim
      IF (m == MonthBefore) THEN
         SSTBefore=SSTClim
      END IF
      IF (m == MonthAfter) THEN
         SSTAfter=SSTClim
      END IF
   END DO
   CLOSE (UNIT=nfclm)

   ! Linear Interpolation in Time for Year, Month, Day and Hour 
   ! of the Initial Condition
   SSTClim=FactorBefore*SSTBefore+FactorAfter*SSTAfter

END SUBROUTINE SSTClimatological


SUBROUTINE SSTClimaWindow ()

   IMPLICIT NONE

   INTEGER :: j

   REAL (KIND=r8) :: Lat, dLat

   ! Get Indices to Use CLimatological SST Out of LatClimSouth to LatClimNorth
   js=0
   jn=0
   dLat=2.0_r8*Lat0/REAL(Jdim-1,r8)
   DO j=1,Jdim
      Lat=Lat0-REAL(j-1,r8)*dLat
      IF (Lat > LatClimSouth) js=j
      IF (Lat > LatClimNorth) jn=j
   END DO
   js=js+1

   WRITE (UNIT=nfprt, FMT='(/,A,/)')' From SSTClimaWindow:'
   WRITE (UNIT=nfprt, FMT='(A,I3,A,F7.3)') &
         ' js = ', js, ' LatClimSouth=', Lat0-REAL(js-1,r8)*dLat
   WRITE (UNIT=nfprt, FMT='(A,I3,A,F7.3,/)') &
         ' jn = ', jn, ' LatClimNorth=', Lat0-REAL(jn-1,r8)*dLat

END SUBROUTINE SSTClimaWindow


END PROGRAM SSTWeekly
