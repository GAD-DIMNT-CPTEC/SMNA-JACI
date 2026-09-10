! Eduardo Khamis - 18-06-2021
module Mod_CO2MonthlyDirec

   use Mod_FastFourierTransform, only : createFFT, destroyFFT
   use Mod_LegendreTransform, only : createGaussRep, createSpectralRep, createLegTrans, &
                                     destroyLegendreObjects
   use Mod_SpectralGrid, only : transp, specCoef2Grid
   use Mod_LinearInterpolation, only: gLatsL=>latOut, &
       initLinearInterpolation, doLinearInterpolation
   use Mod_AreaInterpolation, only: gLatsA=>gLats, &
       initAreaInterpolation, doAreaInterpolation
   use Mod_Get_xMax_yMax, only : getxMaxyMax
   use Mod_Messages, only : msgWarningOut
   use Mod_FileManager
   use Mod_Namelist, only : varCommonNameListData

   implicit none

   public :: getNameCO2MonthlyDirec, initCO2MonthlyDirec, generateCO2MonthlyDirec, shouldRunCO2MonthlyDirec

   private

   include 'files.h'
   include 'pre.h'
   include 'precision.h'
   include 'messages.h'

   ! input parameters

   type CO2MonthlyDirecNameListData
    integer            :: zMax                      ! Number of Layers of the Initial Condition for the Global Model 
    integer            :: xDim                      ! Number of Longitudes For Climatological CO2 data
    integer            :: yDim                      ! Number of Latitudes For Climatological CO2 data
    real(kind = p_r8)  :: co2SeaIce = -0.000_p_r8   ! sem uso? co2SeaIceThreshold no lugar?
    real(kind = p_r8)  :: sstSeaIce = -1.749_p_r8   ! SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
    real(kind = p_r8)  :: latClimSouth =-50.0_p_r8  ! Southern Latitude For Climatological CO2 Data
    real(kind = p_r8)  :: latClimNorth = 60.0_p_r8  ! Northern Latitude For Climatological CO2 Data
    logical            :: climWindow = .false.      ! Flag to Climatological CO2 Data Window
    logical            :: linear = .true.           ! Flag for Linear (T) or Area Weighted (F) Interpolation
    logical            :: linearGrid = .false.      ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
    character(len = 7) :: preffix = 'GANLNMC'       ! Preffix of the Initial Condition for the Global Model
    character(len = 6) :: suffix = 'S.unf.'         ! Suffix of the Initial Condition for the Global Model 
    character(len = maxPathLength) :: dirObsCO2namelist = './' 
   end type CO2MonthlyDirecNameListData

   type(varCommonNameListData)       :: varCommon
   type(CO2MonthlyDirecNameListData) :: var
   namelist /CO2MonthlyDirecNameList/   var

   integer :: xmx, yMaxHf, mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3
   integer :: year, month, day, hour

   real (kind = p_r8) :: lon0, lat0, p_to, co2OpenWater, sstOpenWater,&
                             co2SeaIceThreshold, lapseRate

   logical :: flagInput(5), flagOutput(5)

   integer, dimension (:,:), allocatable :: maskInput

   character (len = 12)  :: gradsTime='  Z         '
   character (len = 10)  :: trunc='T     L   '
   character (len = 7)   :: nLats='.G     '
   character (len = 10)  :: mskfmt = '(      I1)'
   character (len = 255) :: varName='CO2MonthlyDirec'
   character (len = 16)  :: nameLSM='ModelLandSeaMask'
   character (len = 11)  :: fileClmCO2='co2aoi.form'
   character (len = 11)  :: fileClmSST='ersst.form'
   character (len = 23)  :: nameNML='CO2MonthlyDirec.nml'
   character (len = 3), dimension (12) :: monthChar = &
             (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

   integer, dimension (12) :: monthLength = &
            (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

   integer :: nfprt=6    ! Standard print Out
   integer :: nficn      ! To read Topography from Initial Condition
   integer :: nflsm      ! To read Formatted Land Sea Mask
   !integer :: nfsti=30   ! To read Unformatted 1x1 CO2
   integer :: nfclm      ! To read Formatted Climatological CO2  ! duas vezes para arquivos diferentes sem fechar?
   integer :: nfsto      ! To write Unformatted Gaussian Grid CO2
   integer :: nfout      ! To write GrADS Topography, Land Sea, Se Ice and Gauss CO2
   integer :: nfctl      ! To write GrADS Control file
   integer :: nfclm2     ! To read Formatted Climatological SST
   integer :: nfnml      ! To read namelist dirObsCO2namelist

   ! local variables

   ! Reads the Mean Weekly 1x1 CO2 Global From NCEP,
   ! Interpolates it Using Area Weigth Into a Gaussian Grid

   integer :: j, i, m, js, jn, js1, jn1, ja, jb, lRecIn, lRecOut, ios
   integer :: forecastDay
   real (kind = p_r4) :: timeOfDay
   real (kind = p_r8) :: rgCO2Max, rgCO2Min, ggCO2Max, ggCO2Min, mgCO2Max, mgCO2Min
   integer :: iCnDatee(4), currentDate(4), headerAux(5)
   integer, dimension (:,:), allocatable :: landSeaMask, seaIceMask

   real (kind = p_r4), dimension (:), allocatable :: coefTopIn
   real (kind = p_r4), dimension (:,:), allocatable :: co2WklIn_Sea, co2WklIn_land, wrOut
   real (kind = p_r8), dimension (:), allocatable :: coefTop
   real (kind = p_r4), allocatable :: labsGrid(:)
   real (kind = p_r8), dimension (:,:,:), allocatable :: sstIn
   real (kind = p_r8), dimension (:,:), allocatable :: topog, co2Clim,wrIn, &
                   co2In_Sea,co2In_Land,co2Gaus, co2Gaus_Sea,co2Gaus_Land, seaIceFlagIn, seaIceFlagOut
   integer, parameter :: lmon(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

   character (len = 6), dimension (:), allocatable :: co2Label
   integer, parameter :: NSX=55*55
   integer :: LABS(8)
   integer :: LWORD,IDS
   integer :: LRECL,LREC
   integer :: LENGHT,LRECM,LRECN
   integer :: IREC, IREC2,IREC3
   character (len =  16) :: NMCO2
   integer :: NCO2
   integer :: NS,opn
   character (len =  10) :: NDCO2(NSX)
   character (len = 256) :: DRCO2
   namelist /FNCO2NML/ NCO2,NDCO2,DRCO2
   data  NMCO2 /'oico2.          '/
   data  LWORD /1/
   
   character (len = maxPathLength) :: inFileName, outFileName

   character(len = *), parameter :: headerMsg = 'CO2 Monthly Direc          | '


   !namelist /InputDim/ mEnd, zMax, xDim, yDim, &
   !                    co2SeaIce, latClimSouth, latClimNorth, &
   !                    climWindow, linear, linearGrid,grads, &
   !                    dateICn, preffix, suffix, dirMain


contains


  ! Eduardo Khamis - 16-05-2021
  function getNameCO2MonthlyDirec() result(returnModuleName)
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "CO2MonthlyDirec"
  end function getNameCO2MonthlyDirec


  ! Eduardo Khamis - 16-05-2021
  function shouldRunCO2MonthlyDirec() result(shouldRun)
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunCO2MonthlyDirec


  ! Eduardo Khamis - 16-05-2021
  function getOutFileName() result(co2MonthlyDirecOutFilename)
    implicit none
    character(len = maxPathLength) :: co2MonthlyDirecOutFilename

    co2MonthlyDirecOutFilename = trim(varCommon%dirModelIn) // trim(varName) // varCommon%date(1:8) // nLats
  end function getOutFileName


  ! Eduardo Khamis - 16-05-2021
  subroutine initCO2MonthlyDirec (nameListFileUnit, varCommon_)
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
    
    read(unit = nameListFileUnit, nml = CO2MonthlyDirecNameList)
    varCommon = varCommon_

    !call printNameList()

    call getxMaxyMax (varCommon%mEnd, varCommon%xMax, varCommon%yMax)
 
    mEnd1=varCommon%mEnd+1
    mEnd2=varCommon%mEnd+2
    mnwv2=mEnd1*mEnd2
    mnwv0=mnwv2/2
    mnwv3=mnwv2+2*mEnd1
    mnwv1=mnwv3/2
 
    xmx=varCommon%xMax+2
    yMaxHf=varCommon%yMax/2
 
    p_to=273.15_p_r8
    co2OpenWater=-1.7_p_r8+p_to
    co2SeaIceThreshold=0.0_p_r8
    sstOpenWater=-1.7_p_r8+p_to
    var%sstSeaIce=-1.749_p_r8   ! SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
 
    lapseRate=0.0065_p_r8
 
    read (varCommon%date, fmt = '(I4,3I2)') year, month, day, hour
    if (mod(year,4) == 0) monthLength(2)=29
 
    ! For Linear Interpolation
    lon0=0.5_p_r8  ! Start Near Greenwhich
    lat0=89.5_p_r8 ! Start Near North Pole
 
    ! For Area Weighted Interpolation
    allocate (maskInput(var%xDim,var%yDim))
    maskInput=1
    flagInput(1)=.true.   ! Start at North Pole
    flagInput(2)=.true.   ! Start at Prime Meridian
    flagInput(3)=.false.  ! Latitudes Are at North Edge of Box
    flagInput(4)=.false.  ! Longitudes Are at Western Edge of Box
    flagInput(5)=.false.  ! Regular Grid
    flagOutput(1)=.true.  ! Start at North Pole
    flagOutput(2)=.true.  ! Start at Prime Meridian
    flagOutput(3)=.false. ! Latitudes Are at North Edge of Box
    flagOutput(4)=.true.  ! Longitudes Are at Center of Box
    flagOutput(5)=.true.  ! Gaussian Grid
    if(var%linearGrid)then
       write (trunc(2: 6), fmt = '(A1,I4.4)')'L',varCommon%mEnd
       write (trunc(8:10), fmt = '(I3.3)') var%zMax
    else
       write (trunc(2: 6), fmt = '(A1,I4.4)')'Q',varCommon%mEnd
       write (trunc(8:10), fmt = '(I3.3)') var%zMax
    end if
 
    write (nLats(3:7), '(I5.5)') varCommon%yMax
 
    write (mskfmt(2:7), '(I6)') varCommon%xMax
 
    write (gradsTime(1:2), fmt = '(I2.2)') hour
    write (gradsTime(4:5), fmt = '(I2.2)') day
    write (gradsTime(6:8), fmt = '(A3)') monthChar(month)
    write (gradsTime(9:12), fmt = '(I4.4)') year

  end subroutine initCO2MonthlyDirec


  ! Eduardo Khamis - 16-05-2021
  subroutine printNameList()
    implicit none

    write (unit = nfprt, fmt = '(/,A)')    ' &CO2MonthlyDirecNameList'
    write (unit = nfprt, fmt = '(A,I6)')   '          Mend = ', varCommon%mEnd
    write (unit = nfprt, fmt = '(A,I6)')   '          Kmax = ', var%zMax
    write (unit = nfprt, fmt = '(A,I6)')   '          Idim = ', var%xDim
    write (unit = nfprt, fmt = '(A,I6)')   '          Jdim = ', var%yDim
    write (unit = nfprt, fmt = '(A,F6.3)') '     CO2SeaIce = ', var%co2SeaIce
    write (unit = nfprt, fmt = '(A,F6.1)') '  LatClimSouth = ', var%latClimSouth
    write (unit = nfprt, fmt = '(A,F6.1)') '  LatClimNorth = ', var%latClimNorth
    write (unit = nfprt, fmt = '(A,L6)')   '    ClimWindow = ', var%climWindow
    write (unit = nfprt, fmt = '(A,L6)')   '        Linear = ', var%linear
    write (unit = nfprt, fmt = '(A,L6)')   '    LinearGrid = ', var%linearGrid
    write (unit = nfprt, fmt = '(A,L6)')   '         GrADS = ', varCommon%grads
    write (unit = nfprt, fmt = '(A)')      '       DateICn = '//varCommon%date
    write (unit = nfprt, fmt = '(A)')      '       Preffix = '//var%preffix
    write (unit = nfprt, fmt = '(A,/)')    ' /'
  end subroutine printNameList
        
  
  ! Eduardo Khamis - 16-05-2021        
  function generateCO2MonthlyDirec() result(isExecOk)
     implicit none
     logical :: isExecOk

     isExecOk = .false.

     call createSpectralRep (mEnd1, mEnd2, mnwv1)
     call createGaussRep (varCommon%yMax, yMaxHf)
     call createFFT (varCommon%xMax)
     call createLegTrans (mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)
  
     inFileName = trim(var%dirObsCO2namelist)//'/co2mtd.nml'
     nfnml = openFile(trim(inFileName), 'formatted', 'sequential', -1, 'read', 'old')
     if (nfnml < 0) return
     read(nfnml,FNCO2NML)
  
     IDS=index(DRCO2//' ',' ')-1   
     IREC=0
     IREC2=0
     if (LWORD .GT. 0) then
     LRECL=varCommon%xMax*varCommon%yMax*LWORD
     LENGHT=4*LRECL*(NCO2+1)
     else
     LRECL=varCommon%xMax*varCommon%yMax/abs(LWORD)
     LENGHT=4*LRECL*(NCO2+1)*abs(LWORD)
     endif
     LRECM=4+NCO2*10
     LRECN=varCommon%xMax*varCommon%yMax*4
     
     if (LRECM .GT. LRECN)stop ' ERROR: HEADER EXCEED RESERVED RECORD SPACE'
  
     if (var%linear) then
        call initLinearInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, lat0, lon0)
     else
        call initAreaInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
     end if

     call allocateData()  
    
     ! Write out Land-Sea Mask and CO2 Data to Global Model Input
      
     inquire (iolength = lRecOut) wrOut
     outFileName = trim(varCommon%dirModelIn)//trim(varName)//varCommon%date(1:8)//nLats
     nfsto = openFile(trim(outFileName), 'unformatted', 'direct', lRecOut, 'write', 'replace')
     if (nfsto < 0) return     
      
     IREC=IREC+1
     write(*,'(I5,5X,A)')IREC,'NDCO2'
     write(nfsto,REC=IREC)NCO2,(NDCO2(NS),NS=1,NCO2)
   
     if (varCommon%grads) then
        inquire (iolength = lRecOut) wrOut
        outFileName = trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats//'.bin'
        nfout = openFile(trim(outFileName), 'unformatted', 'direct', lRecOut, 'write', 'replace')
        if (nfout < 0) return        
     end if
     !
     ! Read in SpectraL Coefficient of Topography from ICn
     ! to Ensure that Topography is the Same as Used by Model
     !
     inFileName = trim(varCommon%dirModelIn)//trim(var%preffix)//varCommon%date//var%suffix//trunc
     nficn = openFile(trim(inFileName), 'unformatted', 'sequential', -1, 'read', 'old')
     if (nficn < 0) return     
     read  (unit = nficn) forecastDay, timeOfDay, iCnDatee, currentDate
     read  (unit = nficn) coefTopIn
     close (unit = nficn)
     coefTop=real(coefTopIn,p_r8)
     call transp(mnwv2, mEnd1, mEnd2, coefTop)
     call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)
     
   
     ! Read in sst Mask Data Set
  
     inFileName = trim(varCommon%dirClmSST)//fileClmSST
     nfclm2 = openFile(trim(inFileName), 'formatted', 'sequential', -1, 'read', 'old')
     if (nfclm2 < 0) return     
  
     ! Read in Land-Sea Mask Data Set
  
     inFileName = trim(varCommon%dirPreOut)//trim(nameLSM)//nLats
     nflsm = openFile(trim(inFileName), 'unformatted', 'sequential', -1, 'read', 'old') 
     if (nflsm < 0) return     
     read  (unit = nflsm) landSeaMask
     close (unit = nflsm)
  !
  !    WRITE OUT LAND-SEA MASK TO UNIT nfsto CO2 DATA SET
  !    THIS RECORD WILL BE TRANSFERED BY MODEL TO POST-P
  !
     write(*,'(/,6I8,/)')varCommon%xMax,varCommon%yMax,LRECL,LWORD,NCO2,LENGHT
     IREC=IREC+1
     write(*,'(I5,5X,A)')IREC,'LSMK'
     ! Write out Land-Sea Mask to CO2 Data Set
     ! The LSMask will be Transfered by Model to Post-Processing
     wrOut=real(1-2*landSeaMask,p_r4)
     !WRITE (UNIT=nfsto, REC=1) WrOut
     write (nfsto, REC=IREC) wrOut
  
     do m=1,12
        read (unit = nfclm2, fmt = '(8I5)') headerAux
        write (unit = nfprt, fmt = '(/,1X,9I5,/)') m, headerAux
        read (unit = nfclm2, fmt = '(16F5.2)') wrIn
        sstIn(:,:,m)=wrIn(:,:)
        if (maxval(sstIn(:,:,m)) < 100.0_p_r8) sstIn(:,:,m)=sstIn(:,:,m)+p_to
  
     end do
  !**************************************************
  !
  !    LOOP OVER CO2 FILES
  !
     do NS=1,NCO2
  !
  !     INPUT:  UNIT 50 - weekly CO2's
  !
        NMCO2(7:16)=NDCO2(NS)
  !      OPEN(75,FILE=DRCO2(1:IDS)//NMCO2,STATUS='UNKNOWN')
        write(*,*)DRCO2(1:IDS)//NMCO2
        inquire(iolength = lrec)CO2WklIn_sea
        open(75,FILE=DRco2(1:IDS)//NMCO2,&
            ACTION='READ',FORM='UNFORMATTED',&
            ACCESS='DIRECT',STATUS='OLD',RECL=lrec,IOSTAT=opn)
  !
  !     READ 1 DEG X 1 DEG CO2 - DEGREE CELSIUS
  !
  
         if(opn /= 0)then    
            print*,'ERROR AT open OF THE file ',DRCO2(1:IDS)//NMCO2,'STATUS=',opn
         else
           irec3=1
           read(75,rec=irec3)labsGrid
           irec3=irec3+1
           read(75,rec=irec3)co2WklIn_Land
           irec3=irec3+1
           read(75,rec=irec3)co2WklIn_Sea
         end if
         close(75,STATUS='KEEP')
         LABS(1) = int(labsGrid(1))
         LABS(2) = int(labsGrid(2))
         LABS(3) = int(labsGrid(3))
         LABS(4) = int(labsGrid(4))
         LABS(5) = int(labsGrid(5))
         LABS(6) = int(labsGrid(6))
         LABS(7) = int(labsGrid(7))
         LABS(8) = int(labsGrid(8))
  !      READ (75,'(7I5,I10)')LABS
  !      WRITE(* ,'(/,7I5,I10,/)')LABS
  !      READ(75,'(20F4.2)')CO2WklIn
  !      CLOSE(75)
  !
  !     Get CO2Clim Climatological and Index for High Latitude
  !     Substitution of CO2 Actual by Climatology
  !
        if (var%climWindow) then
           call co2Climatological ()
           !IF (MAXVAL(CO2Clim) < -100.0_r8) CO2Clim=CO2Clim+To
           call co2ClimaWindow ()
           jn1=jn-1
           js1=js+1
           ja=jn
           jb=js
        else
           jn=0
           js=var%yDim+1
           jn1=0
           js1=var%yDim+1
           ja=1
           jb=var%yDim
        end if
  
        do j=1,var%yDim
           if (j >= jn .and. j <= js) then
              co2Label(j)='Observ'
              do i=1,var%xDim
                 co2In_Land(i,j)=real(co2WklIn_Land(i,j),p_r8)
                 co2In_Sea(i,j)=real(co2WklIn_Sea(i,j),p_r8)
              end do
           else
              co2Label(j)='Climat'
              do i=1,var%xDim
                 co2In_Land(i,j)=0.0_p_r8!co2Clim(i,j)
                 co2In_Sea(i,j)=0.0_p_r8!co2Clim(i,j)
              end do
           end if
        end do
  
        if (jn1 >= 1) write (unit = nfprt, fmt = '(6(I4,1X,A))') (j,co2Label(j),j=1,jn1)
        write (unit = nfprt, fmt = '(6(I4,1X,A))') (j,co2Label(j),j=ja,jb)
        if (js1 <= var%yDim) write (unit = nfprt, fmt = '(6(I4,1X,A))') (j,co2Label(j),j=js1,var%yDim)
  !
  !     Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
  !     Over Open Water Set Input CO2 = MIN of CO2OpenWater
  !     Over Non Ice Points Before Interpolation
  !     
        seaIceFlagIn=0.0_p_r8
        if (var%sstSeaIce < 100.0_p_r8) var%sstSeaIce=var%sstSeaIce+p_to
        do j=1,var%yDim
           do i=1,var%xDim
              seaIceFlagIn(i,j)=0.0_p_r8
              if (sstIn(i,j,LABS(2)) < var%sstSeaIce) then
  !PRINT*,SSTIn(i,j,LABS(2)) , SSTSeaIce
                 seaIceFlagIn(i,j)=1.0_p_r8
              else
                 sstIn(i,j,LABS(2))=max(sstIn(i,j,LABS(2)),sstOpenWater)
              end if
           end do
        end do
  !      
  !     Min And Max Values of Input CO2
  !
        rgCO2Max=maxval(co2In_Sea)
        rgCO2Min=minval(co2In_Sea)
  !
  ! Interpolate Flag from 1x1 Grid to Gaussian Grid, Fill SeaIceMask=1
  ! Over Interpolated Points With 50% or More Sea Ice, =0 Otherwise
  !
        seaIceFlagOut=0.0_p_r8
        if (var%linear) then
           call doLinearInterpolation (seaIceFlagIn, seaIceFlagOut)
        else
           call doAreaInterpolation (seaIceFlagIn, seaIceFlagOut)
        end if
        seaIceMask=int(seaIceFlagOut+0.5_p_r8)
        where (landSeaMask == 1) seaIceMask=0
  
        ! Interpolate CO2 from 1x1 Grid to Gaussian Grid
        co2Gaus_Land=0.0_p_r8
        if (var%linear) then
           call doLinearInterpolation (co2In_Land, co2Gaus_Land)
        else
           call doAreaInterpolation (co2In_Land, co2Gaus_Land)
        end if
  
        ! Interpolate CO2 from 1x1 Grid to Gaussian Grid
        co2Gaus_Sea=0.0_p_r8
        if (var%linear) then
           call doLinearInterpolation (co2In_Sea, co2Gaus_Sea)
        else
           call doAreaInterpolation (co2In_Sea, co2Gaus_Sea)
        end if
  
        ! Min and Max Values of Gaussian Grid
         ggCO2Max=maxval(co2Gaus_Sea)
         ggCO2Min=minval(co2Gaus_Sea)
  
         do j=1,varCommon%yMax
            do i=1,varCommon%xMax
               if (landSeaMask(i,j) == 1) then
                 ! Set CO2 = Undef Over Land
                 !CO2Gaus(i,j)=Undef
                 co2Gaus(i,j)=(co2Gaus_Land(i,j))/(86400.0_p_r8*real(lmon(LABS(2)),kind=p_r8))!-Topog(i,j)*LapseRate       
               else if (seaIceMask(i,j) == 1) then
                 ! Set CO2 Sea Ice Threshold Minus 1 Over Sea Ice
                 co2Gaus(i,j)=co2SeaIceThreshold !-1.0_r8
               else
                 ! Correct CO2 for Topography, Do Not Create or
                 ! Destroy Sea Ice Via Topography Correction
                 !CO2Gaus_Sea(i,j)=CO2Gaus_Sea(i,j)!-Topog(i,j)*LapseRate
                 co2Gaus(i,j)=(co2Gaus_Sea(i,j)+0.03333334_p_r8 +0.00136_p_r8)/(86400.0_p_r8*real(lmon(LABS(2)),kind=p_r8))!-Topog(i,j)*LapseRate
                 !IF (CO2Gaus_Sea(i,j) < CO2SeaIceThreshold) &
                 !   CO2Gaus_Sea(i,j)=CO2SeaIceThreshold+0.2_r8
              end if
           end do
        end do
  
  
        ! Min and Max Values of Corrected Gaussian Grid CO2 Excluding Land Points
        mgCO2Max=maxval(co2Gaus,MASK=co2Gaus/=p_undef)
        mgCO2Min=minval(co2Gaus,MASK=co2Gaus/=p_undef)
        !
        !     WRITE OUT GAUSSIAN GRID CO2
        !
        !
        ! Write out Gaussian Grid Monthly CO2
        !
        wrOut=real(co2Gaus,p_r4)
        IREC=IREC+1
        write (*,'(I5,5X,A)')IREC,NDCO2(NS)
        write (nfsto, REC=IREC) wrOut
  
        write (unit = nfprt, fmt = '(/,3(A,I2.2),A,I4)') &
              ' Hour = ', hour, ' Day = ', LABS(3), &
              ' Month = ', LABS(2), ' Year = ', LABS(1)
  
        write (unit = nfprt, fmt = '(/,A,3(A,2F8.2,/))') &
           ' Mean Weekly CO2 Interpolation :', &
           ' Regular  Grid CO2: min, max = ', rgCO2Min, rgCO2Max, &
           ' Gaussian Grid CO2: min, max = ', ggCO2Min, ggCO2Max, &
           ' Masked G Grid CO2: min, max = ', mgCO2Min, mgCO2Max
  
        if (varCommon%grads) then
           
           wrOut=real(topog,p_r4)
           IREC2=IREC2+1
           write (unit = nfout, REC=IREC2) wrOut
           
           wrOut=real(1-2*landSeaMask,p_r4)
           IREC2=IREC2+1
           write (unit = nfout, REC=IREC2) wrOut
  
           wrOut=real(seaIceMask,p_r4)
           IREC2=IREC2+1
           write (unit = nfout, REC=IREC2) wrOut
       
           wrOut=real(co2Gaus,p_r4)
           IREC2=IREC2+1
           write (unit = nfout, REC=IREC2) wrOut
        end if
     end do
   
     close (unit = nfout)
     close(nfsto)
  
     if (varCommon%grads) then
        
        ! Write GrADS Control File
        outFileName = trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats//'.ctl'
        nfctl = openFile(trim(outFileName), 'formatted', 'sequential', -1, 'write', 'replace')
        if (nfctl < 0) return        
        write (unit = nfctl, fmt = '(A)') 'DSET '// &
               trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats//'.bin'
        write (unit = nfctl, fmt = '(A)') '*'
        write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
        write (unit = nfctl, fmt = '(A)') '*'
        write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
        write (unit = nfctl, fmt = '(A)') '*'
        write (unit = nfctl, fmt = '(A)') 'TITLE Weekly CO2 on a Gaussian Grid'
        write (unit = nfctl, fmt = '(A)') '*'
        write (unit = nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
                            'XDEF ',varCommon%xMax,' LINEAR ',0.0_p_r8,360.0_p_r8/real(varCommon%xMax,p_r8)
        write (unit = nfctl, fmt = '(A,I5,A)') 'YDEF ',varCommon%yMax,' LEVELS '
        if (var%linear) then
           write (unit = nfctl, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
        else
           write (unit = nfctl, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
        end if
        write (unit = nfctl, fmt = '(A)') 'ZDEF  1 LEVELS 1000'
        write (unit = nfctl, fmt = '(A6,I6,A)') 'TDEF  ',NCO2,' LINEAR '//gradsTime//' 1Mo'
        write (unit = nfctl, fmt = '(A)') 'VARS  4'
        write (unit = nfctl, fmt = '(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
        write (unit = nfctl, fmt = '(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
        write (unit = nfctl, fmt = '(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
        write (unit = nfctl, fmt = '(A)') 'CO2W  0 99 Weekly CO2 Topography Corrected [K]'
        write (unit = nfctl, fmt = '(A)') 'ENDVARS'
  
        close (unit = nfctl)
     end if

     call deallocateData()  
     call destroyFFT()
     call destroyLegendreObjects()
    
     isExecOk = .true.

  end function generateCO2MonthlyDirec
  

  ! Eduardo Khamis - 16-05-2021
  subroutine allocateData()
    implicit none

    allocate (landSeaMask(varCommon%xMax,varCommon%yMax), seaIceMask(varCommon%xMax,varCommon%yMax))
    allocate (coefTopIn(mnwv2), co2WklIn_sea(var%xDim,var%yDim), co2WklIn_Land(var%xDim,var%yDim))
    allocate (coefTop(mnwv2), topog(varCommon%xMax,varCommon%yMax), co2Clim(var%xDim,var%yDim))
    allocate (co2In_Sea(var%xDim,var%yDim), co2In_Land(var%xDim,var%yDim), co2Gaus_Sea(varCommon%xMax,varCommon%yMax))
    allocate (co2Gaus_Land(varCommon%xMax,varCommon%yMax), co2Gaus(varCommon%xMax,varCommon%yMax), wrOut(varCommon%xMax,varCommon%yMax))
    allocate (seaIceFlagIn(var%xDim,var%yDim), seaIceFlagOut(varCommon%xMax,varCommon%yMax))
    allocate (co2Label(var%yDim))
    allocate (labsGrid(var%xDim*var%yDim), sstIn(var%xDim,var%yDim,12), wrIn(var%xDim,var%yDim)) 

  end subroutine allocateData


  ! Eduardo Khamis - 16-05-2021
  subroutine deallocateData()
    implicit none

    deallocate (landSeaMask, seaIceMask)
    deallocate (coefTopIn, co2WklIn_sea, co2WklIn_Land)
    deallocate (coefTop, topog, co2Clim)
    deallocate (co2In_Sea, co2In_Land, co2Gaus_Sea)
    deallocate (co2Gaus_Land, co2Gaus, wrOut)
    deallocate (seaIceFlagIn, seaIceFlagOut)
    deallocate (co2Label)
    deallocate (labsGrid, sstIn, wrIn) 

  end subroutine deallocateData
  

  ! Eduardo Khamis - 16-05-2021
  subroutine co2Climatological ()
  
     implicit none
  
     ! 1950-1979 1 Degree x 1 Degree CO2 
     ! Global NCEP OI Monthly Climatology
     ! Grid Orientation (CO2R):
     ! (1,1) = (0.5_r8W,89.5_r8N)
     ! (Idim,Jdim) = (0.5_r8E,89.5_r8S)
  
     integer :: m, monthBefore, monthAfter
  
     real (kind = p_r8) :: dayHour, dayCorrection, FactorBefore, FactorAfter
  
     integer :: headerAux(8)
  
     real (kind = p_r8) :: co2Before(var%xDim,var%yDim), co2After(var%xDim,var%yDim)
     dayHour=0.0_p_r8
     !DayHour=REAL(Day,r8)+REAL(Hour,r8)/24.0_r8
     monthBefore=LABS(2)
     !IF (DayHour > (1.0_r8+REAL(monthLength(month),r8)/2.0_r8)) &
     !    monthBefore=month
     monthAfter=monthBefore+1
     if (monthBefore < 1) monthBefore=12
     if (monthAfter > 12) monthAfter=1
     dayCorrection=real(monthLength(monthBefore),p_r8)/2.0_p_r8-1.0_p_r8
     !IF (monthBefore == month) DayCorrection=-DayCorrection-2.0_r8
     dayCorrection=-dayCorrection-2.0_p_r8
     FactorAfter=2.0_p_r8*(dayHour+dayCorrection)/ &
                 real(monthLength(monthBefore)+monthLength(monthAfter),p_r8)
     FactorBefore=1.0_p_r8-FactorAfter
  
     write (unit = nfprt, fmt = '(/,A)') ' From CO2Climatological:'
     write (unit = nfprt, fmt = '(/,A,I4,3(A,I2.2))') &
           ' Year = ', LABS(2), ' Month = ', LABS(2), &
           ' Day = ',  LABS(3), ' Hour = ', hour
     write (unit = nfprt, fmt = '(/,2(A,I2))') &
           ' MonthBefore = ', monthBefore, ' MonthAfter = ', monthAfter
     write (unit = nfprt, fmt = '(/,2(A,F9.6),/)') &
           ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter
  
     inFileName = trim(varCommon%dirClmCO2)//fileClmCO2
     nfclm = openFile(trim(inFileName), 'formatted', 'sequential', -1, 'read', 'old')
     if (nfclm < 0) return     
     do m=1,12
        read (unit = nfclm, fmt = '(8I5)') headerAux
        write (unit = nfprt, fmt = '(/,1X,9I5,/)') m, headerAux
        read (unit = nfclm, fmt = '(16F5.2)') co2Clim
        if (m == monthBefore) then
           co2Before=co2Clim
        end if
        if (m == monthAfter) then
           co2After=co2Clim
        end if
     end do
     close (unit = nfclm)
  
     ! Linear Interpolation in Time for Year, Month, Day and Hour 
     ! of the Initial Condition
     co2Clim=FactorBefore*co2Before+FactorAfter*co2After
  
  end subroutine co2Climatological
 

  ! Eduardo Khamis - 16-05-2021
  subroutine co2Climatological_O ()
  
     implicit none
  
     ! 1950-1979 1 Degree x 1 Degree CO2 
     ! Global NCEP OI Monthly Climatology
     ! Grid Orientation (CO2R):
     ! (1,1) = (0.5_r8W,89.5_r8N)
     ! (Idim,Jdim) = (0.5_r8E,89.5_r8S)
  
     integer :: m, monthBefore, monthAfter
  
     real (kind = p_r8) :: dayHour, dayCorrection, FactorBefore, FactorAfter
  
     integer :: headerAux(8)
  
     real (kind = p_r8) :: co2Before(var%xDim,var%yDim), co2After(var%xDim,var%yDim)
  
     dayHour=real(day,p_r8)+real(hour,p_r8)/24.0_p_r8
     monthBefore=month-1
     if (dayHour > (1.0_p_r8+real(monthLength(month),p_r8)/2.0_p_r8)) &
         monthBefore=month
     monthAfter=monthBefore+1
     if (monthBefore < 1) monthBefore=12
     if (monthAfter > 12) monthAfter=1
     dayCorrection=real(monthLength(monthBefore),p_r8)/2.0_p_r8-1.0_p_r8
     if (monthBefore == month) dayCorrection=-dayCorrection-2.0_p_r8
     FactorAfter=2.0_p_r8*(dayHour+dayCorrection)/ &
                 real(monthLength(monthBefore)+monthLength(monthAfter),p_r8)
     FactorBefore=1.0_p_r8-FactorAfter
  
     write (unit = nfprt, fmt = '(/,A)') ' From CO2Climatological:'
     write (unit = nfprt, fmt = '(/,A,I4,3(A,I2.2))') &
           ' Year = ', year, ' Month = ', month, &
           ' Day = ', day, ' Hour = ', hour
     write (unit = nfprt, fmt = '(/,2(A,I2))') &
           ' MonthBefore = ', monthBefore, ' MonthAfter = ', monthAfter
     write (unit = nfprt, fmt = '(/,2(A,F9.6),/)') &
           ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter
  
     inFileName = trim(varCommon%dirClmCO2)//fileClmCO2
     nfclm = openFile(trim(inFileName), 'formatted', 'sequential', -1, 'read', 'old')
     if (nfclm < 0) return     
     do m=1,12
        read (unit = nfclm, fmt = '(8I5)') headerAux
        write (unit = nfprt, fmt = '(/,1X,9I5,/)') m, headerAux
        read (unit = nfclm, fmt = '(16F5.2)') co2Clim
        if (m == monthBefore) then
           co2Before=co2Clim
        end if
        if (m == monthAfter) then
           co2After=co2Clim
        end if
     end do
     close (unit = nfclm)
  
     ! Linear Interpolation in Time for Year, Month, Day and Hour 
     ! of the Initial Condition
     co2Clim=FactorBefore*co2Before+FactorAfter*co2After
  
  end subroutine co2Climatological_O
  

  ! Eduardo Khamis - 16-05-2021
  subroutine co2ClimaWindow ()
  
     implicit none
  
     integer :: j
  
     real (kind = p_r8) :: lat, dLat
  
     ! Get Indices to Use CLimatological CO2 Out of LatClimSouth to LatClimNorth
     js=0
     jn=0
     dLat=2.0_p_r8*lat0/real(var%yDim-1,p_r8)
     do j=1,var%yDim
        lat=lat0-real(j-1,p_r8)*dLat
        if (lat > var%latClimSouth) js=j
        if (lat > var%latClimNorth) jn=j
     end do
     js=js+1
  
     write (unit = nfprt, fmt = '(/,A,/)')' From CO2ClimaWindow:'
     write (unit = nfprt, fmt = '(A,I3,A,F7.3)') &
           ' js = ', js, ' LatClimSouth=', lat0-real(js-1,p_r8)*dLat
     write (unit = nfprt, fmt = '(A,I3,A,F7.3,/)') &
           ' jn = ', jn, ' LatClimNorth=', lat0-real(jn-1,p_r8)*dLat
  
  end subroutine co2ClimaWindow


end module Mod_CO2MonthlyDirec
