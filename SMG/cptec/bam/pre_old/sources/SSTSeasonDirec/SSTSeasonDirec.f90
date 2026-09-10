PROGRAM SSTMonthlyDirec

   USE InputParameters, ONLY : r4, r8, Undef, &
                               nferr, nfprt, nficn, nflsm, nfsti, &
                               nfclm, nfsto, nfout, nfctl, &
                               Mend, Kmax, Imax, Jmax, Idim, Jdim, Mnwv2, &
                               Year, Month, Day, Hour, MonthLength, FirstYear,LastYear,nMonth,MonthChar,&
                               YearInitialCondition ,MonthInitialCondition,FirstYearClim,LastYearClim,&
                               GrADSTime, Trunc, nLats, mskfmt, VarName, NameLSM, &
                               DirPreOut, DirModelIn, DirClmSST, FileClmSST,DirObsSST, &
                               To, SSTSeaIce, LatClimSouth, LatClimNorth, Lat0, &
                               SSTOpenWater, SSTSeaIceThreshold, LapseRate, &
                               ClimWindow, GrADS, Linear, &
                               DateICn, DatePerc,Preffix, Suffix, DirMain, &
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

   INTEGER, DIMENSION (:,:), ALLOCATABLE :: LandSeaMask, SeaIceMask, IntSSTWklIn

   REAL (KIND=r4), DIMENSION (:), ALLOCATABLE :: CoefTopIn

   REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: SSTWklIn, WrOut

   REAL(KIND=4),ALLOCATABLE :: sstpers (:,:,:,:)
   REAL(KIND=4),ALLOCATABLE :: sstanoma(:,:,:,:)
   REAL(KIND=4),ALLOCATABLE :: sststime(:,:,:,:)
   REAL(KIND=4),ALLOCATABLE :: sstclima(:,:,:)
   REAL(KIND=4),ALLOCATABLE :: sstcount(:,:,:)

   REAL (KIND=r8), DIMENSION (:), ALLOCATABLE :: CoefTop

   REAL (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: Topog, SSTClim, &
                   SSTIn, SSTGaus, SeaIceFlagIn, SeaIceFlagOut

   CHARACTER (LEN=6), DIMENSION (:), ALLOCATABLE :: SSTLabel
   INTEGER, PARAMETER :: NSX=550*550
   INTEGER :: LABS(8)
   INTEGER :: LWORD,IDS
   INTEGER :: LRECL
   INTEGER :: LENGHT,LRECM,LRECN
   INTEGER :: IREC, IREC2,nREC
   CHARACTER(LEN= 16) :: NMSST
   INTEGER :: NSST
   INTEGER :: NS
   INTEGER :: nYearsFirstSerie
   INTEGER :: nMonthFirstSerie 
   INTEGER :: nYearsLastSerie
   INTEGER :: nMonthLastSerie 
   
   INTEGER :: FirstMonthSerie 
   INTEGER :: LastMonthSerie 

   INTEGER :: iyear
   INTEGER :: imes 
   INTEGER :: nyear
   INTEGER :: nmes 

   LOGICAL :: lexist,test
   CHARACTER (LEN= 10) :: NDSST(NSX)
   CHARACTER (LEN=256) :: DRSST
   NAMELIST /FNSSTNML/ NSST,NDSST,DRSST
   DATA  NMSST /'oisst.          '/
   DATA  LWORD /1/
   test=.true.

   CALL InitInputParameters ()
   CALL CreateSpectralRep ()
   CALL CreateGaussRep ()
   CALL CreateFFT ()
   CALL CreateLegTrans ()

   OPEN(11,FILE=TRIM(DirMain)//TRIM(DirObsSST)//'/sstmtd.nml',STATUS='OLD')
   READ(11,FNSSTNML)

   IDS=INDEX(DRSST//' ',' ')-1   
   IREC=0
   IREC2=0
   IF (LWORD .GT. 0) THEN
   LRECL=IMAX*JMAX*LWORD
   LENGHT=4*LRECL*(NSST+1)
   ELSE
   LRECL=IMAX*JMAX/ABS(LWORD)
   LENGHT=4*LRECL*(NSST+1)*ABS(LWORD)
   ENDIF
   LRECM=4+NSST*10
   LRECN=IMAX*JMAX*4
   
   IF (LRECM .GT. LRECN)STOP ' ERROR: HEADER EXCEED RESERVED RECORD SPACE'

   IF (Linear) THEN
      CALL InitLinearInterpolation ()
   ELSE
      CALL InitAreaInterpolation ()
   END IF

   ALLOCATE (LandSeaMask(Imax,Jmax), SeaIceMask(Imax,Jmax))
   ALLOCATE (CoefTopIn(Mnwv2), SSTWklIn(Idim,Jdim), IntSSTWklIn(Idim,Jdim))
   ALLOCATE (CoefTop(Mnwv2), Topog(Imax,Jmax), SSTClim(Idim,Jdim))
   ALLOCATE (SSTIn(Idim,Jdim), SSTGaus(Imax,Jmax), WrOut(Imax,Jmax))
   ALLOCATE (SeaIceFlagIn(Idim,Jdim), SeaIceFlagOut(Imax,Jmax))
   ALLOCATE (SSTLabel(Jdim))
   ALLOCATE(sstpers (1:Idim,1:Jdim,FirstYear:LastYear+1,1:nMonth));sstpers=0.0
   ALLOCATE(sstanoma(1:Idim,1:Jdim,FirstYear:LastYear+1,1:nMonth));sstanoma=0.0
   ALLOCATE(sststime(1:Idim,1:Jdim,FirstYear:LastYear+1,1:nMonth));sststime=0.0
   ALLOCATE(sstclima(1:Idim,1:Jdim,1:nMonth));sstclima=0.0
   ALLOCATE(sstcount(1:Idim,1:Jdim,1:nMonth));sstcount=0.0

  
   ! Write out Land-Sea Mask and SST Data to Global Model Input
    
   INQUIRE (IOLENGTH=LRecOut) WrOut
   OPEN (UNIT=nfsto, FILE=TRIM(DirMain)//DirModelIn//TRIM(VarName)//DateICn(1:8)//'-'//DatePerc(1:8)//nLats, &
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
   WRITE(*,'(I5,5X,A)')IREC,'NDSST'
   WRITE(nfsto,REC=IREC)NSST,(NDSST(NS),NS=1,NSST)
 
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
!    WRITE OUT LAND-SEA MASK TO UNIT nfsto SST DATA SET
!    THIS RECORD WILL BE TRANSFERED BY MODEL TO POST-P
!
   WRITE(*,'(/,6I8,/)')IMAX,JMAX,LRECL,LWORD,NSST,LENGHT
   IREC=IREC+1
   WRITE(*,'(I5,5X,A)')IREC,'LSMK'
   ! Write out Land-Sea Mask to SST Data Set
   ! The LSMask will be Transfered by Model to Post-Processing
   WrOut=REAL(1-2*LandSeaMask,r4)
   !WRITE (UNIT=nfsto, REC=1) WrOut
   WRITE(nfsto,REC=IREC)WrOut
!**************************************************
!
!    LOOP OVER SST FILES
!
   DO NS=1,NSST
!
!     INPUT:  UNIT 50 - weekly sst's
!
      NMSST(7:16)=NDSST(NS)
      INQUIRE (FILE=DRSST(1:IDS)//NMSST,exist=lexist)
      IF(lexist) THEN
         OPEN(75,FILE=DRSST(1:IDS)//NMSST,STATUS='UNKNOWN')
         !
         !     READ 1 DEG X 1 DEG SST - DEGREE CELSIUS
         !
         READ (75,'(7I5,I10)')LABS
         !WRITE(* ,'(/,7I5,I10,/)')LABS
         READ(75,'(20I6)')IntSSTWklIn
         SSTWklIn=IntSSTWklIn/100.0
         iyear=LABS(1)
         imes =LABS(2)
         nYearsLastSerie= iyear
         nMonthLastSerie= imes 
         IF(test)THEN
            nYearsFirstSerie =iyear
            nMonthFirstSerie =imes
            WRITE (GrADSTime(6:8), FMT='(A3)') MonthChar(imes)
            WRITE (GrADSTime(9:12), FMT='(I4.4)') iyear
            test=.FALSE.
         END IF
         DO j=1,Jdim
            DO i=1,Idim
               sststime(i,j,iyear,imes) = SSTWklIn(i,j)
            END DO
         END DO
         WRITE(*,*)DRSST(1:IDS)//NMSST,MAXVAL(sststime(:,:,iyear,imes)),MINVAL(sststime(:,:,iyear,imes))
         CLOSE(75,STATUS='KEEP')
      ELSE
         PRINT*,'* The  ', DRSST(1:IDS)//NMSST, ' file not exist*'
         STOP
      END IF
   END DO

   IF(FirstYearClim >= nYearsFirstSerie .and. LastYearClim <= nYearsLastSerie)THEN
      WRITE(*,*)'OK CLIMATOLOGY PERIOD ',FirstYearClim,' >= ', nYearsFirstSerie,' .and. ',LastYearClim,' <=', nYearsLastSerie
   ELSE
      WRITE(*,*)'ERROR CLIMATOLOGY PERIOD ',FirstYearClim,' >= ', nYearsFirstSerie,' .and. ',LastYearClim,' <=', nYearsLastSerie
      STOP
   END IF


   DO iyear=FirstYearClim,LastYearClim
      DO imes=1,nMonth
         DO j=1,Jdim
            DO i=1,Idim
               sstclima(i,j,imes) = sstclima(i,j,imes) + sststime(i,j,iyear,imes)
               sstcount(i,j,imes) = sstcount(i,j,imes) + 1.0
            END DO
         END DO
      END DO
   END DO

   DO imes=1,nMonth
      DO j=1,Jdim
         DO i=1,Idim
            IF(sstcount(i,j,imes) /= 0.0 )THEN
               sstclima(i,j,imes) = sstclima(i,j,imes) / sstcount(i,j,imes)
            END IF
         END DO
      END DO
   END DO
   
   DO iyear=nYearsFirstSerie,nYearsLastSerie
      FirstMonthSerie=1 
      LastMonthSerie =12
      IF(iyear == nYearsFirstSerie)THEN
         FirstMonthSerie=nMonthFirstSerie
         LastMonthSerie =12
      END IF
      IF(iyear == nYearsLastSerie )THEN
         FirstMonthSerie=1
         LastMonthSerie=nMonthLastSerie
      END IF
      DO imes=FirstMonthSerie,LastMonthSerie
         DO j=1,Jdim
            DO i=1,Idim 
               IF(sstcount(i,j,imes) /= 0.0  )THEN
                  sstanoma(i,j,iyear,imes) = sststime(i,j,iyear,imes) - sstclima(i,j,imes) 
               ELSE
                  sstanoma(i,j,iyear,imes) = 0.0
               END IF
            END DO
         END DO
      END DO
   END DO
   IF(YearInitialCondition >= nYearsFirstSerie .and. YearInitialCondition<=nYearsLastSerie)THEN
      WRITE(*,*)'OK PERSI SST PERIOD ',YearInitialCondition,' >= ', nYearsFirstSerie,' .and. ',YearInitialCondition,' <=', nYearsLastSerie

   ELSE
      WRITE(*,*)'ERROR PERSI SST PERIOD ',YearInitialCondition,' >= ', nYearsFirstSerie,' .and. ',YearInitialCondition,' <=', nYearsLastSerie
      STOP
   END IF
   IF(FirstYear == nYearsFirstSerie .and. LastYear == nYearsLastSerie)THEN
      WRITE(*,*)'OK PERSI SST PERIOD ',FirstYear,' == ', nYearsFirstSerie,' .and. ',LastYear,' ==', nYearsLastSerie

   ELSE
      WRITE(*,*)'ERROR PERSI SST PERIOD ',FirstYear,' == ', nYearsFirstSerie,' .and. ',LastYear,' ==', nYearsLastSerie
      STOP
   END IF

   DO iyear=FirstYear,LastYear+1
      FirstMonthSerie=1 
      LastMonthSerie =12
      IF(iyear == nYearsFirstSerie)THEN
         FirstMonthSerie=nMonthFirstSerie
      END IF

      DO imes=FirstMonthSerie,LastMonthSerie

         IF(iyear>=YearInitialCondition ) THEN
            IF(iyear>YearInitialCondition )THEN
                     PRINT*,'CLIM1+ANOM',iyear,YearInitialCondition,imes,MonthInitialCondition
               DO j=1,Jdim
                  DO i=1,Idim
                     IF(sstcount(i,j,imes) /= 0.0 )THEN
                        sstpers(i,j,iyear,imes) = sstclima(i,j,imes) + sstanoma(i,j,YearInitialCondition,MonthInitialCondition)   !sst pesistida climatologia+anomaliajan
                     END IF
                     IF(sstpers(i,j,iyear,imes) < -1.79)THEN
                        sstpers(i,j,iyear,imes) = -1.79        !sst pesistida ice
                     END IF
                  END DO
               END DO
            ELSEIF(iyear==YearInitialCondition .and. imes == MonthInitialCondition)THEN
                    PRINT*,'OBS2-OBS1',iyear,YearInitialCondition,imes,MonthInitialCondition
               DO j=1,Jdim
                  DO i=1,Idim
                     IF(sstcount(i,j,imes) /= 0.0 )THEN
                        sstpers(i,j,iyear,imes) = sststime(i,j,iyear,imes)   !sst observada
                     END IF
                     IF(sstpers(i,j,iyear,imes) < -1.79)THEN
                        sstpers(i,j,iyear,imes) = -1.79        !sst pesistida ice
                     END IF
                  END DO
               END DO
            ELSEIF(iyear==YearInitialCondition .and. imes < MonthInitialCondition)THEN
                    PRINT*,'OBS2-OBS2',iyear,YearInitialCondition,imes,MonthInitialCondition
               DO j=1,Jdim
                  DO i=1,Idim
                     IF(sstcount(i,j,imes) /= 0.0 )THEN
                        sstpers(i,j,iyear,imes) = sststime(i,j,iyear,imes)   !sst observada
                     END IF
                     IF(sstpers(i,j,iyear,imes) < -1.79)THEN
                        sstpers(i,j,iyear,imes) = -1.79        !sst pesistida ice
                     END IF
                  END DO
               END DO
            ELSEIF(iyear==YearInitialCondition .and. imes > MonthInitialCondition)THEN 
                     PRINT*,'CLIM3+ANOM',iyear,YearInitialCondition,imes,MonthInitialCondition
               DO j=1,Jdim
                  DO i=1,Idim
                     IF(sstcount(i,j,imes) /= 0.0 )THEN
                        sstpers(i,j,iyear,imes) = sstclima(i,j,imes) + sstanoma(i,j,YearInitialCondition,MonthInitialCondition)   !sst pesistida climatologia+anomaliajan
                     END IF
                     IF(sstpers(i,j,iyear,imes) < -1.79)THEN
                        sstpers(i,j,iyear,imes) = -1.79        !sst pesistida ice
                     END IF
                  END DO
               END DO
            ELSE
                    PRINT*,'OBS2-OBS2',iyear,YearInitialCondition,imes,MonthInitialCondition
               DO j=1,Jdim
                  DO i=1,Idim
                     IF(sstcount(i,j,imes) /= 0.0 )THEN
                        sstpers(i,j,iyear,imes) = sststime(i,j,iyear,imes)   !sst observada
                     END IF
                     IF(sstpers(i,j,iyear,imes) < -1.79)THEN
                        sstpers(i,j,iyear,imes) = -1.79        !sst pesistida ice
                     END IF
                  END DO
               END DO
            END IF
         ELSE
                  PRINT*,'OBS',iyear,YearInitialCondition,imes,MonthInitialCondition
            DO j=1,Jdim
               DO i=1,Idim
                  IF(sstcount(i,j,imes) /= 0.0 )THEN
                     sstpers(i,j,iyear,imes) = sststime(i,j,iyear,imes)   !sst observada
                  END IF
                  IF(sstpers(i,j,iyear,imes) < -1.79)THEN
                     sstpers(i,j,iyear,imes) = -1.79        !sst pesistida ice
                  END IF
               END DO
            END DO

         END IF
      END DO
   END DO

   DO imes=1,nMonth
      DO j=1,Jdim
         DO i=1,Idim
            IF(sstclima(i,j,imes) < -1.79)THEN
               sstclima(i,j,imes) = -1.79        !sst pesistida climatologia+anomaliajan
            END IF
         END DO
      END DO
   END DO

!**************************************************
!
!    LOOP OVER SST FILES
!
   nREC=0
   NS=0
   DO iyear=FirstYear,LastYear+1
      FirstMonthSerie=1 
      LastMonthSerie =12
      IF(iyear == nYearsFirstSerie)THEN
         FirstMonthSerie=nMonthFirstSerie
      END IF

      DO imes=FirstMonthSerie,LastMonthSerie
         nREC=nREC+1
         NS=NS+1
         IF(NS <=NSST )THEN
!
!     INPUT:  UNIT 50 - weekly sst's
!
            NMSST(7:16)=NDSST(NS)
            OPEN(75,FILE=DRSST(1:IDS)//NMSST,STATUS='UNKNOWN')
            WRITE(*,*)DRSST(1:IDS)//NMSST
!
!     READ 1 DEG X 1 DEG SST - DEGREE CELSIUS
!
            READ (75,'(7I5,I10)')LABS
            WRITE(* ,'(/,7I5,I10,/)')LABS
            !READ(75,'(20I6)')IntSSTWklIn
            !SSTWklIn=IntSSTWklIn/100.0
            nyear=LABS(1)
            nmes =LABS(2)
            CLOSE(75)
         ELSE
            LABS(1)=iyear
            LABS(2)=imes
            nyear=iyear
            nmes =imes 
         END IF 
         DO j=1,Jdim
            DO i=1,Idim
               SSTWklIn(i,j)=sstpers(i,j,nyear,nmes)
            END DO
         END DO

!
!     Get SSTClim Climatological and Index for High Latitude
!     Substitution of SST Actual by Climatology
!
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
               SSTLabel(j)='Observ'
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
!
!     Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
!     Over Open Water Set Input SST = MIN of SSTOpenWater
!     Over Non Ice Points Before Interpolation
!     
      !PRINT*,SSTIn
      !STOP
         IF (SSTSeaIce < 100.0_r8) SSTSeaIce=SSTSeaIce+To
         IF (MAXVAL(SSTIn) < 100.0_r8) SSTIn=SSTIn+To
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
!      
!     Min And Max Values of Input SST
!
         RGSSTMax=MAXVAL(SSTIn)
         RGSSTMin=MINVAL(SSTIn)
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

      ! Interpolate SST from 1x1 Grid to Gaussian Grid
         SSTGaus=0.0_r8
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
         !
         !     WRITE OUT GAUSSIAN GRID SST
         !
         !
            ! Write out Gaussian Grid Monthly SST
      !
         WrOut=REAL(SSTGaus,r4)
         IREC=IREC+1
         WRITE(*,'(I5,5X,A)')IREC,NDSST(NS)
         WRITE(nfsto,REC=IREC)WrOut

         WRITE (UNIT=nfprt, FMT='(/,3(A,I2.2),A,I4)') &
               ' Hour = ', Hour, ' Day = ', LABS(3), &
               ' Month = ', LABS(2), ' Year = ', LABS(1)

         WRITE (UNIT=nfprt, FMT='(/,A,3(A,2F8.2,/))') &
            ' Mean Weekly SST Interpolation :', &
            ' Regular  Grid SST: Min, Max = ', RGSSTMin, RGSSTMax, &
            ' Gaussian Grid SST: Min, Max = ', GGSSTMin, GGSSTMax, &
            ' Masked G Grid SST: Min, Max = ', MGSSTMin, MGSSTMax

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

            WrOut=REAL(SSTGaus,r4)
            IREC2=IREC2+1
            WRITE (UNIT=nfout, REC=IREC2) WrOut
         END IF
      END DO
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
      WRITE (UNIT=nfctl, FMT='(A6,I6,A)') 'TDEF  ',nREC,' LINEAR '//GrADSTime//' 1Mo'
      WRITE (UNIT=nfctl, FMT='(A)') 'VARS  4'
      WRITE (UNIT=nfctl, FMT='(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      WRITE (UNIT=nfctl, FMT='(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      WRITE (UNIT=nfctl, FMT='(A)') 'SSTW  0 99 Weekly SST Topography Corrected [K]'
      WRITE (UNIT=nfctl, FMT='(A)') 'ENDVARS'

      CLOSE (UNIT=nfctl)
   END IF

PRINT *, "*** SSTMonthlyDirec ENDS NORMALLY ***"

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

   WRITE (UNIT=nfprt, FMT='(/,A)') ' From SSTClimatological:'
   WRITE (UNIT=nfprt, FMT='(/,A,I4,3(A,I2.2))') &
         ' Year = ', LABS(2), ' Month = ', LABS(2), &
         ' Day = ',  LABS(3), ' Hour = ', Hour
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

SUBROUTINE SSTClimatological_O ()

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

END SUBROUTINE SSTClimatological_O



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


END PROGRAM SSTMonthlyDirec
