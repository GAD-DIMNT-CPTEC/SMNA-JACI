PROGRAM CO2MonthlyDirec

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nficn, nflsm, nfsti, &
                               nfclm, nfsto, nfout, nfctl,nfclm2, &
                               Mend, Kmax, Imax, Jmax, Idim, Jdim, Mnwv2, &
                               Year, Month, Day, Hour, MonthLength, &
                               GrADSTime, Trunc, nLats, mskfmt, VarName, NameLSM, &
                               DirPreOut, DirModelIn, DirClmCO2, DirClmSST,FileClmCO2,DirObsCO2,FileClmSST, &
                               To, CO2SeaIce, LatClimSouth, LatClimNorth, Lat0, &
                               CO2OpenWater, CO2SeaIceThreshold, LapseRate,SSTOpenWater, SSTSeaIce,&
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

   ! Reads the Mean Weekly 1x1 CO2 Global From NCEP,
   ! Interpolates it Using Area Weigth Into a Gaussian Grid

   INTEGER :: j, i, m, js, jn, js1, jn1, ja, jb, LRecIn, LRecOut, ios

   INTEGER :: ForecastDay

   REAL (KIND=r4) :: TimeOfDay

   REAL (KIND=r8) :: RGCO2Max, RGCO2Min, GGCO2Max, GGCO2Min, MGCO2Max, MGCO2Min

   INTEGER :: ICnDate(4), CurrentDate(4),Header(5)

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask, SeaIceMask

   REAL (KIND=r4), DIMENSION (:), ALLOCATABLE :: CoefTopIn

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: CO2WklIn_Sea,CO2WklIn_land, WrOut

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: CoefTop
   REAL(KIND=r4), ALLOCATABLE :: LABSGRID(:)
   REAL (KIND=r8), DIMENSION (:,:,:), ALLOCATABLE :: SSTIn
   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Topog, CO2Clim,WrIn, &
                   CO2In_Sea,CO2In_Land,CO2Gaus, CO2Gaus_Sea,CO2Gaus_LAnd, SeaIceFlagIn, SeaIceFlagOut
   INTEGER, PARAMETER :: lmon(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

   CHARACTER (LEN=6), DIMENSION (:), ALLOCATABLE :: CO2Label
   INTEGER, PARAMETER :: NSX=55*55
   INTEGER :: LABS(8)
   INTEGER :: LWORD,IDS
   INTEGER :: LRECL,LREC
   INTEGER :: LENGHT,LRECM,LRECN
   INTEGER :: IREC, IREC2,IREC3
   CHARACTER(LEN= 16) :: NMCO2
   INTEGER :: NCO2
   INTEGER :: NS,opn
   CHARACTER (LEN= 10) :: NDCO2(NSX)
   CHARACTER (LEN=256) :: DRCO2
   NAMELIST /FNCO2NML/ NCO2,NDCO2,DRCO2
   DATA  NMCO2 /'oico2.          '/
   DATA  LWORD /1/
   

   CALL InitInputParameters ()
   CALL CreateSpectralRep ()
   CALL CreateGaussRep ()
   CALL CreateFFT ()
   CALL CreateLegTrans ()

   OPEN(11,FILE=TRIM(DirMain)//TRIM(DirObsCO2)//'/co2mtd.nml',STATUS='OLD')
   READ(11,FNCO2NML)

   IDS=INDEX(DRCO2//' ',' ')-1   
   IREC=0
   IREC2=0
   IF (LWORD .GT. 0) THEN
   LRECL=IMAX*JMAX*LWORD
   LENGHT=4*LRECL*(NCO2+1)
   ELSE
   LRECL=IMAX*JMAX/ABS(LWORD)
   LENGHT=4*LRECL*(NCO2+1)*ABS(LWORD)
   ENDIF
   LRECM=4+NCO2*10
   LRECN=IMAX*JMAX*4
   
   IF (LRECM .GT. LRECN)STOP ' ERROR: HEADER EXCEED RESERVED RECORD SPACE'

   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (LandSeaMask(Imax,Jmax), SeaIceMask(Imax,Jmax))
   ALLOCATE (CoefTopIn(Mnwv2), CO2WklIn_sea(Idim,Jdim),CO2WklIn_Land(Idim,Jdim))
   ALLOCATE (CoefTop(Mnwv2), Topog(Imax,Jmax), CO2Clim(Idim,Jdim))
   ALLOCATE (CO2In_Sea(Idim,Jdim),CO2In_Land(Idim,Jdim), CO2Gaus_Sea(Imax,Jmax))
   ALLOCATE (CO2Gaus_Land(Imax,Jmax),CO2Gaus(Imax,Jmax), WrOut(Imax,Jmax))
   ALLOCATE (SeaIceFlagIn(Idim,Jdim), SeaIceFlagOut(Imax,Jmax))
   ALLOCATE (CO2Label(Jdim))
   ALLOCATE (LABSGRID(Idim*Jdim),SSTIn(Idim,Jdim,12),WrIn(Idim,Jdim))
  
   ! Write out Land-Sea Mask and CO2 Data to Global Model Input
    
   INQUIRE (IOLENGTH=LRecOut) WrOut
   OPEN (UNIT=nfsto, FILE=TRIM(DirMain)//DirModelIn//TRIM(VarName)//DateICn(1:8)//nLats, &
         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
         STATUS='REPLACE', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirMain)//DirModelIn//TRIM(VarName)//DateICn(1:8)//nLats, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
    
   IREC=IREC+1
   WRITE(*,'(I5,5X,A)')IREC,'NDCO2'
   WRITE(nfsto,REC=IREC)NCO2,(NDCO2(NS),NS=1,NCO2)
 
   IF (GrADS) THEN
      INQUIRE (IOLENGTH=LRecOut) WrOut
      OPEN (UNIT=nfout, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//DateICn(1:8)//nLats//'.bin', &
        FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRecOut, ACTION='WRITE', &
        STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE (UNIT=nferr, FMT='(3A,I4)') &
           ' ** (Error) ** Open file ', &
             TRIM(DirMain)//DirPreOut//TRIM(VarName)//DateICn(1:8)//nLats, &
           ' returned IOStat = ', ios
        STOP  ' ** (Error) **'
      END IF
   END IF
   !
   ! Read in SpectraL Coefficient of Topography from ICn
   ! to Ensure that Topography is the Same as Used by Model
   !
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
 
   ! Read in sst Mask Data Set


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
!
!    WRITE OUT LAND-SEA MASK TO UNIT nfsto CO2 DATA SET
!    THIS RECORD WILL BE TRANSFERED BY MODEL TO POST-P
!
   WRITE(*,'(/,6I8,/)')IMAX,JMAX,LRECL,LWORD,NCO2,LENGHT
   IREC=IREC+1
   WRITE(*,'(I5,5X,A)')IREC,'LSMK'
   ! Write out Land-Sea Mask to CO2 Data Set
   ! The LSMask will be Transfered by Model to Post-Processing
   WrOut=REAL(1-2*LandSeaMask,r4)
   !WRITE (UNIT=nfsto, REC=1) WrOut
   WRITE(nfsto,REC=IREC)WrOut

   DO m=1,12
      READ (UNIT=nfclm2, FMT='(8I5)') Header
      WRITE (UNIT=nfprt, FMT='(/,1X,9I5,/)') m, Header
      READ (UNIT=nfclm2, FMT='(16F5.2)') WrIn
      SSTIn(:,:,m)=WrIn(:,:)
      IF (MAXVAL(SSTIn(:,:,m)) < 100.0_r8) SSTIn(:,:,m)=SSTIn(:,:,m)+To

   END DO
!**************************************************
!
!    LOOP OVER CO2 FILES
!
   DO NS=1,NCO2
!
!     INPUT:  UNIT 50 - weekly CO2's
!
      NMCO2(7:16)=NDCO2(NS)
!      OPEN(75,FILE=DRCO2(1:IDS)//NMCO2,STATUS='UNKNOWN')
      WRITE(*,*)DRCO2(1:IDS)//NMCO2
      INQUIRE(IOLENGTH=lrec)CO2WklIn_sea
      OPEN(75,FILE=DRCO2(1:IDS)//NMCO2,&
          ACTION='READ',FORM='UNFORMATTED',&
          ACCESS='DIRECT',STATUS='OLD',RECL=lrec,IOSTAT=opn)
!
!     READ 1 DEG X 1 DEG CO2 - DEGREE CELSIUS
!

       IF(opn /= 0)THEN    
          PRINT*,'ERROR AT OPEN OF THE FILE ',DRCO2(1:IDS)//NMCO2,'STATUS=',opn
       ELSE
         irec3=1
         READ(75,rec=irec3)LABSGRID
         irec3=irec3+1
         READ(75,rec=irec3)CO2WklIn_Land
         irec3=irec3+1
         READ(75,rec=irec3)CO2WklIn_Sea
       END IF
       CLOSE(75,STATUS='KEEP')
       LABS(1) = INT(LABSGRID(1))
       LABS(2) = INT(LABSGRID(2))
       LABS(3) = INT(LABSGRID(3))
       LABS(4) = INT(LABSGRID(4))
       LABS(5) = INT(LABSGRID(5))
       LABS(6) = INT(LABSGRID(6))
       LABS(7) = INT(LABSGRID(7))
       LABS(8) = INT(LABSGRID(8))
!      READ (75,'(7I5,I10)')LABS
!      WRITE(* ,'(/,7I5,I10,/)')LABS
!      READ(75,'(20F4.2)')CO2WklIn
!      CLOSE(75)
!
!     Get CO2Clim Climatological and Index for High Latitude
!     Substitution of CO2 Actual by Climatology
!
      IF (ClimWindow) THEN
         CALL CO2Climatological ()
         !IF (MAXVAL(CO2Clim) < -100.0_r8) CO2Clim=CO2Clim+To
         CALL CO2ClimaWindow ()
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
            CO2Label(j)='Observ'
            DO i=1,Idim
               CO2In_Land(i,j)=REAL(CO2WklIn_LAnd(i,j),r8)
               CO2In_Sea(i,j)=REAL(CO2WklIn_Sea(i,j),r8)
            END DO
         ELSE
            CO2Label(j)='Climat'
            DO i=1,Idim
               CO2In_Land(i,j)=0.0_r8!CO2Clim(i,j)
               CO2In_Sea(i,j)=0.0_r8!CO2Clim(i,j)
            END DO
         END IF
      END DO

      IF (jn1 >= 1) WRITE (UNIT=nfprt, FMT='(6(I4,1X,A))') (j,CO2Label(j),j=1,jn1)
      WRITE (UNIT=nfprt, FMT='(6(I4,1X,A))') (j,CO2Label(j),j=ja,jb)
      IF (js1 <= Jdim) WRITE (UNIT=nfprt, FMT='(6(I4,1X,A))') (j,CO2Label(j),j=js1,Jdim)
!
!     Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
!     Over Open Water Set Input CO2 = MIN of CO2OpenWater
!     Over Non Ice Points Before Interpolation
!     
      SeaIceFlagIn=0.0_r8
      IF (SSTSeaIce < 100.0_r8) SSTSeaIce=SSTSeaIce+To
      DO j=1,Jdim
         DO i=1,Idim
            SeaIceFlagIn(i,j)=0.0_r8
            IF (SSTIn(i,j,LABS(2)) < SSTSeaIce) THEN
!PRINT*,SSTIn(i,j,LABS(2)) , SSTSeaIce
               SeaIceFlagIn(i,j)=1.0_r8
            ELSE
               SSTIn(i,j,LABS(2))=MAX(SSTIn(i,j,LABS(2)),SSTOpenWater)
            END IF
         END DO
      END DO
!      
!     Min And Max Values of Input CO2
!
      RGCO2Max=MAXVAL(CO2In_Sea)
      RGCO2Min=MINVAL(CO2In_Sea)
!
! Interpolate Flag from 1x1 Grid to Gaussian Grid, Fill SeaIceMask=1
! Over Interpolated Points With 50% or More Sea Ice, =0 Otherwise
!
      SeaIceFlagOut=0.0_r8
      IF (Linear) THEN
         CALL DoLinearInterpolation (SeaIceFlagIn, SeaIceFlagOut)
      ELSE
         CALL DoAreaInterpolation (SeaIceFlagIn, SeaIceFlagOut)
      END IF
      SeaIceMask=INT(SeaIceFlagOut+0.5_r8)
      WHERE (LandSeaMask == 1) SeaIceMask=0

      ! Interpolate CO2 from 1x1 Grid to Gaussian Grid
      CO2Gaus_Land=0.0_r8
      IF (Linear) THEN
         CALL DoLinearInterpolation (CO2In_Land, CO2Gaus_Land)
      ELSE
         CALL DoAreaInterpolation (CO2In_Land, CO2Gaus_Land)
      END IF

      ! Interpolate CO2 from 1x1 Grid to Gaussian Grid
      CO2Gaus_Sea=0.0_r8
      IF (Linear) THEN
         CALL DoLinearInterpolation (CO2In_Sea, CO2Gaus_Sea)
      ELSE
         CALL DoAreaInterpolation (CO2In_Sea, CO2Gaus_Sea)
      END IF

      ! Min and Max Values of Gaussian Grid
       GGCO2Max=MAXVAL(CO2Gaus_Sea)
       GGCO2Min=MINVAL(CO2Gaus_Sea)

       DO j=1,Jmax
          DO i=1,Imax
             IF (LandSeaMask(i,j) == 1) THEN
               ! Set CO2 = Undef Over Land
               !CO2Gaus(i,j)=Undef
               CO2Gaus(i,j)=(CO2Gaus_Land(i,j))/(86400.0_r8*REAL(lmon(LABS(2)),kind=r8))!-Topog(i,j)*LapseRate       
             ELSE IF (SeaIceMask(i,j) == 1) THEN
               ! Set CO2 Sea Ice Threshold Minus 1 Over Sea Ice
               CO2Gaus(i,j)=CO2SeaIceThreshold !-1.0_r8
             ELSE
               ! Correct CO2 for Topography, Do Not Create or
               ! Destroy Sea Ice Via Topography Correction
               !CO2Gaus_Sea(i,j)=CO2Gaus_Sea(i,j)!-Topog(i,j)*LapseRate
               CO2Gaus(i,j)=(CO2Gaus_Sea(i,j)+0.03333334_r8 +0.00136_r8)/(86400.0_r8*REAL(lmon(LABS(2)),kind=r8))!-Topog(i,j)*LapseRate
               !IF (CO2Gaus_Sea(i,j) < CO2SeaIceThreshold) &
               !   CO2Gaus_Sea(i,j)=CO2SeaIceThreshold+0.2_r8
            END IF
         END DO
      END DO


      ! Min and Max Values of Corrected Gaussian Grid CO2 Excluding Land Points
      MGCO2Max=MAXVAL(CO2Gaus,MASK=CO2Gaus/=Undef)
      MGCO2Min=MINVAL(CO2Gaus,MASK=CO2Gaus/=Undef)
      !
      !     WRITE OUT GAUSSIAN GRID CO2
      !
      !
      ! Write out Gaussian Grid Monthly CO2
      !
      WrOut=REAL(CO2Gaus,r4)
      IREC=IREC+1
      WRITE(*,'(I5,5X,A)')IREC,NDCO2(NS)
      WRITE(nfsto,REC=IREC)WrOut

      WRITE (UNIT=nfprt, FMT='(/,3(A,I2.2),A,I4)') &
            ' Hour = ', Hour, ' Day = ', LABS(3), &
            ' Month = ', LABS(2), ' Year = ', LABS(1)

      WRITE (UNIT=nfprt, FMT='(/,A,3(A,2F8.2,/))') &
         ' Mean Weekly CO2 Interpolation :', &
         ' Regular  Grid CO2: Min, Max = ', RGCO2Min, RGCO2Max, &
         ' Gaussian Grid CO2: Min, Max = ', GGCO2Min, GGCO2Max, &
         ' Masked G Grid CO2: Min, Max = ', MGCO2Min, MGCO2Max

      IF (GrADS) THEN
         
         WrOut=REAL(Topog,r4)
         IREC2=IREC2+1
         WRITE (UNIT=nfout, REC=IREC2) WrOut
         
         WrOut=REAL(1-2*LandSeaMask,r4)
         IREC2=IREC2+1
         WRITE (UNIT=nfout, REC=IREC2) WrOut

         WrOut=REAL(SeaIceMask,r4)
         IREC2=IREC2+1
         WRITE (UNIT=nfout, REC=IREC2) WrOut
     
         WrOut=REAL(CO2Gaus,r4)
         IREC2=IREC2+1
         WRITE (UNIT=nfout, REC=IREC2) WrOut
      END IF
   END DO
 
   CLOSE (UNIT=nfout)
   CLOSE(nfsto)

   IF (GrADS) THEN
      
      ! Write GrADS Control File
      OPEN (UNIT=nfctl, FILE=TRIM(DirMain)//DirPreOut//TRIM(VarName)//DateICn(1:8)//nLats//'.ctl', &
            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
            ACTION='WRITE', STATUS='REPLACE', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (UNIT=nferr, FMT='(3A,I4)') &
               ' ** (Error) ** Open file ', &
                 TRIM(DirMain)//DirPreOut//TRIM(VarName)//DateICn(1:8)//nLats//'.ctl', &
               ' returned IOStat = ', ios
         STOP  ' ** (Error) **'
      END IF
      WRITE (UNIT=nfctl, FMT='(A)') 'DSET '// &
             TRIM(DirMain)//DirPreOut//TRIM(VarName)//DateICn(1:8)//nLats//'.bin'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'OPTIONS YREV BIG_ENDIAN'
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A,1PG12.5)') 'UNDEF ', Undef
      WRITE (UNIT=nfctl, FMT='(A)') '*'
      WRITE (UNIT=nfctl, FMT='(A)') 'TITLE Weekly CO2 on a Gaussian Grid'
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
      WRITE (UNIT=nfctl, FMT='(A6,I6,A)') 'TDEF  ',NCO2,' LINEAR '//GrADSTime//' 1Mo'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  4'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      WRITE (UNIT=nfctl, FMT='(A)') 'CO2W  0 99 Weekly CO2 Topography Corrected [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'

      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** CO2MonthlyDirec ENDS NORMALLY ***"

CONTAINS


SUBROUTINE CO2Climatological ()

   IMPLICIT NONE

   ! 1950-1979 1 Degree x 1 Degree CO2 
   ! Global NCEP OI Monthly Climatology
   ! Grid Orientation (CO2R):
   ! (1,1) = (0.5_r8W,89.5_r8N)
   ! (Idim,Jdim) = (0.5_r8E,89.5_r8S)

   INTEGER :: m, MonthBefore, MonthAfter

   REAL (KIND=r8) :: DayHour, DayCorrection, FactorBefore, FactorAfter

   INTEGER :: Header(8)

   REAL (KIND=r8) :: CO2Before(Idim,Jdim), CO2After(Idim,Jdim)
   DayHour=0.0_r8
   !DayHour=REAL(Day,r8)+REAL(Hour,r8)/24.0_r8
   MonthBefore=LABS(2)
   !IF (DayHour > (1.0_r8+REAL(MonthLength(Month),r8)/2.0_r8)) &
   !    MonthBefore=Month
   MonthAfter=MonthBefore+1
   IF (MonthBefore < 1) MonthBefore=12
   IF (MonthAfter > 12) MonthAfter=1
   DayCorrection=REAL(MonthLength(MonthBefore),r8)/2.0_r8-1.0_r8
   !IF (MonthBefore == Month) DayCorrection=-DayCorrection-2.0_r8
   DayCorrection=-DayCorrection-2.0_r8
   FactorAfter=2.0_r8*(DayHour+DayCorrection)/ &
               REAL(MonthLength(MonthBefore)+MonthLength(MonthAfter),r8)
   FactorBefore=1.0_r8-FactorAfter

   WRITE (UNIT=nfprt, FMT='(/,A)') ' From CO2Climatological:'
   WRITE (UNIT=nfprt, FMT='(/,A,I4,3(A,I2.2))') &
         ' Year = ', LABS(2), ' Month = ', LABS(2), &
         ' Day = ',  LABS(3), ' Hour = ', Hour
   WRITE (UNIT=nfprt, FMT='(/,2(A,I2))') &
         ' MonthBefore = ', MonthBefore, ' MonthAfter = ', MonthAfter
   WRITE (UNIT=nfprt, FMT='(/,2(A,F9.6),/)') &
         ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter

   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//TRIM(DirClmCO2)//FileClmCO2, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirClmCO2)//FileClmCO2, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO m=1,12
      READ (UNIT=nfclm, FMT='(8I5)') Header
      WRITE (UNIT=nfprt, FMT='(/,1X,9I5,/)') m, Header
      READ (UNIT=nfclm, FMT='(16F5.2)') CO2Clim
      IF (m == MonthBefore) THEN
         CO2Before=CO2Clim
      END IF
      IF (m == MonthAfter) THEN
         CO2After=CO2Clim
      END IF
   END DO
   CLOSE (UNIT=nfclm)

   ! Linear Interpolation in Time for Year, Month, Day and Hour 
   ! of the Initial Condition
   CO2Clim=FactorBefore*CO2Before+FactorAfter*CO2After

END SUBROUTINE CO2Climatological

SUBROUTINE CO2Climatological_O ()

   IMPLICIT NONE

   ! 1950-1979 1 Degree x 1 Degree CO2 
   ! Global NCEP OI Monthly Climatology
   ! Grid Orientation (CO2R):
   ! (1,1) = (0.5_r8W,89.5_r8N)
   ! (Idim,Jdim) = (0.5_r8E,89.5_r8S)

   INTEGER :: m, MonthBefore, MonthAfter

   REAL (KIND=r8) :: DayHour, DayCorrection, FactorBefore, FactorAfter

   INTEGER :: Header(8)

   REAL (KIND=r8) :: CO2Before(Idim,Jdim), CO2After(Idim,Jdim)

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

   WRITE (UNIT=nfprt, FMT='(/,A)') ' From CO2Climatological:'
   WRITE (UNIT=nfprt, FMT='(/,A,I4,3(A,I2.2))') &
         ' Year = ', Year, ' Month = ', Month, &
         ' Day = ', Day, ' Hour = ', Hour
   WRITE (UNIT=nfprt, FMT='(/,2(A,I2))') &
         ' MonthBefore = ', MonthBefore, ' MonthAfter = ', MonthAfter
   WRITE (UNIT=nfprt, FMT='(/,2(A,F9.6),/)') &
         ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter

   OPEN (UNIT=nfclm, FILE=TRIM(DirMain)//TRIM(DirClmCO2)//FileClmCO2, &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', &
         ACTION='READ', STATUS='OLD', IOSTAT=ios)
   IF (ios /= 0) THEN
      WRITE (UNIT=nferr, FMT='(3A,I4)') &
            ' ** (Error) ** Open file ', &
              TRIM(DirClmCO2)//FileClmCO2, &
            ' returned IOStat = ', ios
      STOP  ' ** (Error) **'
   END IF
   DO m=1,12
      READ (UNIT=nfclm, FMT='(8I5)') Header
      WRITE (UNIT=nfprt, FMT='(/,1X,9I5,/)') m, Header
      READ (UNIT=nfclm, FMT='(16F5.2)') CO2Clim
      IF (m == MonthBefore) THEN
         CO2Before=CO2Clim
      END IF
      IF (m == MonthAfter) THEN
         CO2After=CO2Clim
      END IF
   END DO
   CLOSE (UNIT=nfclm)

   ! Linear Interpolation in Time for Year, Month, Day and Hour 
   ! of the Initial Condition
   CO2Clim=FactorBefore*CO2Before+FactorAfter*CO2After

END SUBROUTINE CO2Climatological_O



SUBROUTINE CO2ClimaWindow ()

   IMPLICIT NONE

   INTEGER :: j

   REAL (KIND=r8) :: Lat, dLat

   ! Get Indices to Use CLimatological CO2 Out of LatClimSouth to LatClimNorth
   js=0
   jn=0
   dLat=2.0_r8*Lat0/REAL(Jdim-1,r8)
   DO j=1,Jdim
      Lat=Lat0-REAL(j-1,r8)*dLat
      IF (Lat > LatClimSouth) js=j
      IF (Lat > LatClimNorth) jn=j
   END DO
   js=js+1

   WRITE (UNIT=nfprt, FMT='(/,A,/)')' From CO2ClimaWindow:'
   WRITE (UNIT=nfprt, FMT='(A,I3,A,F7.3)') &
         ' js = ', js, ' LatClimSouth=', Lat0-REAL(js-1,r8)*dLat
   WRITE (UNIT=nfprt, FMT='(A,I3,A,F7.3,/)') &
         ' jn = ', jn, ' LatClimNorth=', Lat0-REAL(jn-1,r8)*dLat

END SUBROUTINE CO2ClimaWindow


END PROGRAM CO2MonthlyDirec
