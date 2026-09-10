!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SSTMonthlyDirec </br></br>
!#
!# **Brief**: Module responsible for reading reads the NCEP Global Monthly 
!# Climatic Sea Surface Temperature (TSM) datasets, interpolates them using the
!# area weight in a Gaussian grid of the model, and writes the output file </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/sstaoi.form </br>
!# &bull; pre/datasst/oiv2monthly/sstmtd.nml </br>
!# &bull; pre/dataout/ModelLandSeaMask.GZZZZZ (Ex.: pre/dataout/ModelLandSeaMask.G00450) </br>
!# &bull; model/datain/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ (Ex.: model/datain/GANLNMC2015043000S.unf.TQ0299L064)
!# </br></br>
!# 
!# **Files out:**
!#
!# &bull; model/datain/SSTMonthlyDirecYYYYMMDD.GZZZZZ </br>
!# &bull; pre/dataout/SSTMonthlyDirecYYYYMMDD.GZZZZZ.bin
!# </br></br>
!# 
!# When in the names of these files there is the prefix SMT instead of NMC, is
!# because they suffered dampening (smoothing). </br></br>
!# 
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti      - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita               - version: 1.1.1 </li>
!#  <li>17-09-2019 - Eduardo Khamis       - version: 2.0.1 </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> </br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> </br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_SSTMonthlyDirec

  use Mod_FastFourierTransform, only : createFFT, destroyFFT
  use Mod_LegendreTransform, only : createGaussRep, createSpectralRep, createLegTrans, &
                                    destroyLegendreObjects
  use Mod_SpectralGrid, only : transp, specCoef2Grid
  use Mod_LinearInterpolation, only: gLatsL=>LatOut, &
      initLinearInterpolation, doLinearInterpolation
  use Mod_AreaInterpolation, only: gLatsA=>gLats, &
      initAreaInterpolation, doAreaInterpolation
  use Mod_Get_xMax_yMax, only : getxMaxyMax
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgOut, msgInLineFormatOut, msgNewLine

  implicit none

  public :: generateSSTMonthlyDirec
  public :: getNameSSTMonthlyDirec
  public :: initSSTMonthlyDirec
  public :: shouldRunSSTMonthlyDirec

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  ! input parameters
  integer :: xmx, yMaxHf, mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3
  integer :: year, month, day, hour

  real(kind = p_r8) :: p_lon0, p_lat0, p_to, sstOpenWater, &
                            sstSeaIceThreshold, lapseRate

  logical :: polarMean
  logical :: flagInput(5), flagOutput(5)

  integer, dimension (:,:), allocatable :: maskInput

  character(len = 12) :: gradsTime='  Z         '
                                   ! TQ0062L028
  character(len = 10) :: trunc='T     L   '
  character(len = 7) :: nLats='.G     '
  character(len = 10) :: mskfmt = '(      I1)'
  character(len = 255) :: varName='SSTMonthlyDirec'
  character(len = 16) :: nameLSM='ModelLandSeaMask'
  character(len = 11) :: fileClmSST='sstaoi.form'

  character(len = 3), dimension (12) :: monthChar = &
            (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

  integer, dimension (12) :: monthLength = &
           (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  integer :: nfprt=6 ! Standard print Out
  integer :: nficn   ! To read Topography from Initial Condition
  integer :: nflsm   ! To read Formatted Land Sea Mask
  integer :: nfsti=30   ! To read Unformatted 1x1 SST
  integer :: nfclm   ! To read Formatted Climatological SST
  integer :: nfsto   ! To write Unformatted Gaussian Grid SST
  integer :: nfout   ! To write GrADS Topography, Land Sea, Se Ice and Gauss SST
  integer :: nfctl   ! To write GrADS Control file

  ! local variables

  ! Reads the Mean Weekly 1x1 SST Global From NCEP,
  ! Interpolates it Using Area Weigth Into a Gaussian Grid

  integer :: j, i, m, js, jn, js1, jn1, ja, jb, lRecIn, lRecOut, ios
  integer :: forecastDay
  integer :: iCnDate(4), currentDate(4)
  integer, dimension (:,:), allocatable :: landSeaMask, seaIceMask, intSSTWklIn

  real(kind = p_r4) :: timeOfDay
  real(kind = p_r8) :: rgSSTMax, rgSSTMin, ggSSTMax, ggSSTMin, mgSSTMax, mgSSTMin
  real(kind = p_r4), dimension (:), allocatable :: coefTopIn
  real(kind = p_r4), dimension (:,:), allocatable :: sstWklIn, wrOut
  real(kind = p_r8), dimension (:), allocatable :: coefTop
  real(kind = p_r8), dimension (:,:), allocatable :: topog, sstClim, &
                  sstIn, sstGaus, seaIceFlagIn, seaIceFlagOut

  character(len = 6), dimension (:), allocatable :: sstLabel
  integer, parameter :: nsx=550*550
  integer :: lAbs(8)
  integer :: lWord,ids
  integer :: lRecL
  integer :: length,lRecM,lRecN
  integer :: iRec, iRec2
  character(len = 16) :: nmSST
  integer :: nSST
  integer :: ns
  character(len =  10) :: ndSST(nsx)
  character(len = 256) :: drSST
  namelist /fnSSTNML/ nSST,ndSST,drSST
  data  nmSST /'oisst.          '/
  data  lWord /1/

  type SSTMonthlyDirecNameListData
    integer            :: zmax          
    !# Number of Layers 
    integer            :: xDim         
    !# Number of Longitudes
    integer            :: yDim          
    !# Number of Latitudes
    character(len = maxPathLength) :: dirSSTnamelist = './' 
    !# Observed SST namelist file "sstmtd.nml"
    real(kind = p_r8)  :: sstSeaIce = -1.749   
    !# SST Value in Celsius Degree Over Sea Ice (-1.749 NCEP, -1.799 CAC)
    real(kind = p_r8)  :: latClimSouth =-50.0  ! Southern Latitude For Climatological SST Data
    real(kind = p_r8)  :: latClimNorth = 60.0  ! Northern Latitude For Climatological SST Data
    logical            :: climWindow = .false. ! Flag to Climatological SST Data Window
    logical            :: linear = .true.      ! Flag for Linear (T) or Area Weighted (F) Interpolation
    logical            :: linearGrid = .false. ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
    character(len = 7) :: preffix = 'GANLNMC'  ! Preffix of the Initial Condition for the Global Model
    character(len = 6) :: suffix = 'S.unf.'    ! Suffix of the Initial Condition for the Global Model
  end type SSTMonthlyDirecNameListData

  type(varCommonNameListData)       :: varCommon
  type(SSTMonthlyDirecNameListData) :: var
  namelist /SSTMonthlyDirecNameList/   var

  character(len = *), parameter :: headerMsg = 'SST Monthly Direc          | '


  contains


  function getNameSSTMonthlyDirec() result(returnModuleName)
    !# Returns SSTMonthlyDirec Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SSTMonthlyDirec Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "SSTMonthlyDirec"
  end function getNameSSTMonthlyDirec


  function shouldRunSSTMonthlyDirec() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunSSTMonthlyDirec


  function getOutFileName() result(sstMonthlyDirecOutFilename)
    !# Gets SSTMonthlyDirec Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SSTMonthlyDirec Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: sstMonthlyDirecOutFilename

    sstMonthlyDirecOutFilename = trim(varCommon%dirModelIn) // trim(varName) // varCommon%date(1:8) // nLats
  end function getOutFileName


  subroutine initSSTMonthlyDirec(nameListFileUnit, varCommon_)
    !# Initialization of SSTMonthlyDirec module
    !# ---
    !# @info
    !# **Brief:** Initialization of SSTMonthlyDirec module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo 
    implicit none
 
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
   
    read(unit = nameListFileUnit, nml = SSTMonthlyDirecNameList)
    varCommon = varCommon_ 

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
    sstOpenWater=-1.7_p_r8+p_to
    sstSeaIceThreshold=271.2_p_r8
 
    lapseRate=0.0065_p_r8
 
    read (varCommon%date, fmt = '(I4,3I2)') year, month, day, hour
    if (mod(year,4) == 0) monthLength(2)=29
 
    ! For Linear Interpolation
    p_lon0=0.5_p_r8  ! Start Near Greenwhich
    p_lat0=89.5_p_r8 ! Start Near North Pole
 
    ! For Area Weighted Interpolation
    allocate (maskInput(var%xDim,var%yDim))
    maskInput=1
    polarMean=.false.
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
 
  end subroutine initSSTMonthlyDirec
  
  
  function generateSSTMonthlyDirec() result(isExecOk)
    !# Generates SSTMonthlyDirec
    !# ---
    !# @info
    !# **Brief:** Generates SSTMonthlyDirec. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none

    logical :: isExecOk
    integer :: inFileUnit

    isExecOk = .false. 

    !call printNameList() 

    call createSpectralRep (mEnd1, mEnd2, mnwv1)
    call createGaussRep (varCommon%yMax, yMaxHf)
    call createFFT (varCommon%xMax)
    call createLegTrans (mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)
 
    inFileUnit = openFile(trim(var%dirSSTnamelist)//'/sstmtd.nml', 'formatted', 'sequential', -1, 'read', 'old')
    if(inFileUnit < 0) return
    read(inFileUnit,fnSSTNML)
    close(inFileUnit)

 
    ids=index(drSST//' ',' ')-1   
    iRec=0
    iRec2=0
    if (lWord .GT. 0) then
    lRecL=varCommon%xMax*varCommon%yMax*lWord
    length=4*lRecL*(nSST+1)
    else
    lRecL=varCommon%xMax*varCommon%yMax/abs(lWord)
    length=4*lRecL*(nSST+1)*abs(lWord)
    endif
    lRecM=4+nSST*10
    lRecN=varCommon%xMax*varCommon%yMax*4
    
    if (lRecM .GT. lRecN)stop ' ERROR: HEADER EXCEED RESERVED RECORD SPACE'
 
    if (var%linear) then
       call initLinearInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, p_lat0, p_lon0)
    else
       call initAreaInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
    end if
 
    allocate (landSeaMask(varCommon%xMax,varCommon%yMax), seaIceMask(varCommon%xMax,varCommon%yMax))
    allocate (coefTopIn(mnwv2), sstWklIn(var%xDim,var%yDim), intSSTWklIn(var%xDim,var%yDim))
    allocate (coefTop(mnwv2), topog(varCommon%xMax,varCommon%yMax), sstClim(var%xDim,var%yDim))
    allocate (sstIn(var%xDim,var%yDim), sstGaus(varCommon%xMax,varCommon%yMax), wrOut(varCommon%xMax,varCommon%yMax))
    allocate (seaIceFlagIn(var%xDim,var%yDim), seaIceFlagOut(varCommon%xMax,varCommon%yMax))
    allocate (sstLabel(var%yDim))
 
   
    ! Write out Land-Sea Mask and SST Data to Global Model Input
     
    inquire (iolength=lRecOut) wrOut
    nfsto = openFile(trim(varCommon%dirModelIn)//trim(varName)//varCommon%date(1:8)//nLats, 'unformatted', 'direct', lRecOut, 'write', 'replace') 
    if(nfsto < 0) return
    iRec=iRec+1
    write(*,'(I5,5X,A)')iRec,'NDSST'
    write(nfsto,REC=iRec)nSST,(ndSST(ns),ns=1,nSST)
  
    if (varCommon%grads) then
       inquire (iolength=lRecOut) wrOut
       nfout = openFile(trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats//'.bin', 'unformatted', 'direct', lRecOut, 'write', 'replace')
       if(nfout < 0) return
    end if
    !
    ! Read in SpectraL Coefficient of Topography from ICn
    ! to Ensure that Topography is the Same as Used by Model
    !
    nficn = openFile(trim(varCommon%dirModelIn)//var%preffix//varCommon%date//var%suffix//trunc, 'unformatted', 'sequential', -1, 'read', 'old')
    if(nficn < 0) return
    read(unit = nficn) forecastDay, timeOfDay, iCnDate, currentDate
    read(unit = nficn) coefTopIn
    close(unit = nficn)
    coefTop=real(coefTopIn,p_r8)
    call transp(mnwv2, mEnd1, mEnd2, coefTop)
    call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)
 
    ! Read in Land-Sea Mask Data Set
 
    nflsm = openFile(trim(varCommon%dirPreOut)//nameLSM//nLats, 'unformatted', 'sequential', -1, 'read', 'old') !!! lRecIn ?
    if(nflsm < 0) return
    read(unit = nflsm) landSeaMask
    close(unit = nflsm)
 !
 !    WRITE OUT LAND-SEA MASK TO UNIT nfsto SST DATA SET
 !    THIS RECORD WILL BE TRANSFERED BY MODEL TO POST-P
 !
    write(*,'(/,6I8,/)')varCommon%xMax,varCommon%yMax,lRecL,lWord,nSST,length
    iRec=iRec+1
    write(*,'(I5,5X,A)')iRec,'LSMK'
    ! Write out Land-Sea Mask to SST Data Set
    ! The LSMask will be Transfered by Model to Post-Processing
    wrOut=real(1-2*landSeaMask,p_r4)
    !WRITE(unit = nfsto, REC=1) wrOut
    write(nfsto,REC=iRec)wrOut
 
 !**************************************************
 !
 !    LOOP OVER SST FILES
 !
    do ns=1,nSST
 !
 !     INPUT:  UNIT 50 - weekly sst's
 !
       nmSST(7:16)=ndSST(ns)
       inFileUnit = openFile(drSST(1:ids)//nmSST, 'formatted', 'sequential', -1, 'read', 'unknown')
       if(inFileUnit < 0) return
       write(*,*)drSST(1:ids)//nmSST
 !
 !     READ 1 DEG X 1 DEG SST - DEGREE CELSIUS
 !
       read (inFileUnit,'(7I5,I10)') lAbs
       write(* ,'(/,7I5,I10,/)') lAbs
       read (inFileUnit,'(20I6)') intSSTWklIn
       sstWklIn = real(intSSTWklIn) / 100.0_p_r4
       close(inFileUnit)

       do j = 1, var%yDim
          do i = 1, var%xDim
             if (sstWklIn(i,j) < -1.79_p_r4) then
                sstWklIn(i,j) = -1.79_p_r4        !sea ice temperature
             endif
          enddo
       enddo
       
 !
 !     Get SSTClim Climatological and Index for High Latitude
 !     Substitution of SST Actual by Climatology
 !
       if (var%climWindow) then
          call sstClimatological ()
          if (maxval(sstClim) < 100.0_p_r8) sstClim=sstClim+p_to
          call sstClimaWindow ()
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
             sstLabel(j)='Observ'
             do i=1,var%xDim
                sstIn(i,j)=real(sstWklIn(i,j),p_r8)
             end do
          else
             sstLabel(j)='Climat'
             do i=1,var%xDim
                sstIn(i,j)=sstClim(i,j)
             end do
          end if
       end do
       if (jn1 >= 1) write(unit = nfprt, fmt = '(6(I4,1X,A))') (j,sstLabel(j),j=1,jn1)
       write(unit = nfprt, fmt = '(6(I4,1X,A))') (j,sstLabel(j),j=ja,jb)
       if (js1 <= var%yDim) write(unit = nfprt, fmt = '(6(I4,1X,A))') (j,sstLabel(j),j=js1,var%yDim)
 !
 !     Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
 !     Over Open Water Set Input SST = MIN of SSTOpenWater
 !     Over Non Ice Points Before Interpolation
 !     
       !PRINT*,SSTIn
       !STOP
       if (var%sstSeaIce < 100.0_p_r8) var%sstSeaIce=var%sstSeaIce+p_to
       if (maxval(sstIn) < 100.0_p_r8) sstIn=sstIn+p_to
       do j=1,var%yDim
          do i=1,var%xDim
             seaIceFlagIn(i,j)=0.0_p_r8
             if (sstIn(i,j) < var%sstSeaIce) then
                seaIceFlagIn(i,j)=1.0_p_r8
             else
                sstIn(i,j)=MAX(sstIn(i,j),sstOpenWater)
             end if
          end do
       end do
 !      
 !     Min And Max Values of Input SST
 !
       rgSSTMax=maxval(sstIn)
       rgSSTMin=minval(sstIn)
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
 
       ! Interpolate SST from 1x1 Grid to Gaussian Grid
       sstGaus=0.0_p_r8
       if (var%linear) then
          call doLinearInterpolation (sstIn, sstGaus)
       else
          call doAreaInterpolation (sstIn, sstGaus)
       end if
 
       ! Min and Max Values of Gaussian Grid
        ggSSTMax=maxval(sstGaus)
        ggSSTMin=minval(sstGaus)
 
        do j=1,varCommon%yMax
           do i=1,varCommon%xMax
              if (landSeaMask(i,j) == 1) then
                ! Set SST = Undef Over Land
                sstGaus(i,j)=p_undef
              else if (seaIceMask(i,j) == 1) then
                ! Set SST Sea Ice Threshold Minus 1 Over Sea Ice
                sstGaus(i,j)=sstSeaIceThreshold-1.0_p_r8
              else
                ! Correct SST for Topography, Do Not Create or
                ! Destroy Sea Ice Via Topography Correction
                sstGaus(i,j)=sstGaus(i,j)-topog(i,j)*lapseRate
                if (sstGaus(i,j) < sstSeaIceThreshold) &
                   sstGaus(i,j)=sstSeaIceThreshold+0.2_p_r8
             end if
          end do
       end do
 
 
       ! Min and Max Values of Corrected Gaussian Grid SST Excluding Land Points
       mgSSTMax=maxval(sstGaus,mask=sstGaus/=p_undef)
       mgSSTMin=minval(sstGaus,mask=sstGaus/=p_undef)
       !
       !     WRITE OUT GAUSSIAN GRID SST
       !
       !
       ! Write out Gaussian Grid Monthly SST
       !
       wrOut=real(sstGaus,p_r4)
       iRec=iRec+1
       write(*,'(I5,5X,A)')iRec,ndSST(ns)
       write(nfsto,REC=iRec)wrOut
 
       write(unit = nfprt, fmt = '(/,3(A,I2.2),A,I4)') &
             ' Hour = ', hour, ' Day = ', lAbs(3), &
             ' Month = ', lAbs(2), ' Year = ', lAbs(1)
 
       write(unit = nfprt, fmt = '(/,A,3(A,2F8.2,/))') &
          ' Mean Weekly SST Interpolation :', &
          ' Regular  Grid SST: min, max = ', rgSSTMin, rgSSTMax, &
          ' Gaussian Grid SST: min, max = ', ggSSTMin, ggSSTMax, &
          ' Masked G Grid SST: min, max = ', mgSSTMin, mgSSTMax
 
       if (varCommon%grads) then
          
          wrOut=real(topog,p_r4)
          iRec2=iRec2+1
          write(unit = nfout, REC=iRec2) wrOut
          
          wrOut=real(1-2*landSeaMask,p_r4)
          iRec2=iRec2+1
          write(unit = nfout, REC=iRec2) wrOut
 
          wrOut=real(seaIceMask,p_r4)
          iRec2=iRec2+1
          write(unit = nfout, REC=iRec2) wrOut

          wrOut=real(sstGaus,p_r4)
          iRec2=iRec2+1
          write(unit = nfout, REC=iRec2) wrOut
       end if
    end do
  
    close(unit = nfout)
    close(nfsto)
 
    if (varCommon%grads) then
       
       ! Write GrADS Control File
       nfctl = openFile(trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
       if(nfctl < 0) return 
       write(unit = nfctl, fmt = '(A)') 'DSET ^'// &
              trim(varName)//varCommon%date(1:8)//nLats//'.bin'
       write(unit = nfctl, fmt = '(A)') '*'
       write(unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
       write(unit = nfctl, fmt = '(A)') '*'
       write(unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
       write(unit = nfctl, fmt = '(A)') '*'
       write(unit = nfctl, fmt = '(A)') 'TITLE Weekly SST on a Gaussian Grid'
       write(unit = nfctl, fmt = '(A)') '*'
       write(unit = nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
                           'XDEF ',varCommon%xMax,' LINEAR ',0.0_p_r8,360.0_p_r8/real(varCommon%xMax,p_r8)
       write(unit = nfctl, fmt = '(A,I5,A)') 'YDEF ',varCommon%yMax,' LEVELS '
       if (var%linear) then
          write(unit = nfctl, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
       else
          write(unit = nfctl, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
       end if
       write(unit = nfctl, fmt = '(A)') 'ZDEF  1 LEVELS 1000'
       write(unit = nfctl, fmt = '(A6,I6,A)') 'TDEF  ',nSST,' LINEAR '//gradsTime//' 1Mo'
       write(unit = nfctl, fmt = '(A)') 'VARS  4'
       write(unit = nfctl, fmt = '(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
       write(unit = nfctl, fmt = '(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
       write(unit = nfctl, fmt = '(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
       write(unit = nfctl, fmt = '(A)') 'SSTW  0 99 Weekly SST Topography Corrected [K]'
       write(unit = nfctl, fmt = '(A)') 'ENDVARS'
 
       close(unit = nfctl)
    end if
 
    deallocate (landSeaMask, seaIceMask)
    deallocate (coefTopIn, sstWklIn, intSSTWklIn)
    deallocate (coefTop, topog, sstClim)
    deallocate (sstIn, sstGaus, wrOut)
    deallocate (seaIceFlagIn, seaIceFlagOut)
    deallocate (sstLabel)

    call destroyFFT()
    call destroyLegendreObjects()

    isExecOk = .true. 
 
  end function generateSSTMonthlyDirec

  
  subroutine printNameList()
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none

    write(unit = nfprt, fmt = '(/,A)')    ' &InputDim'
    write(unit = nfprt, fmt = '(A,I6)')   '          Mend = ', varCommon%mEnd
    write(unit = nfprt, fmt = '(A,I6)')   '          Kmax = ', var%zMax
    write(unit = nfprt, fmt = '(A,I6)')   '          Idim = ', var%xDim
    write(unit = nfprt, fmt = '(A,I6)')   '          Jdim = ', var%yDim
    write(unit = nfprt, fmt = '(A,F6.3)') '     SSTSeaIce = ', var%sstSeaIce
    write(unit = nfprt, fmt = '(A,F6.1)') '  LatClimSouth = ', var%latClimSouth
    write(unit = nfprt, fmt = '(A,F6.1)') '  LatClimNorth = ', var%latClimNorth
    write(unit = nfprt, fmt = '(A,L6)')   '    ClimWindow = ', var%climWindow
    write(unit = nfprt, fmt = '(A,L6)')   '        Linear = ', var%linear
    write(unit = nfprt, fmt = '(A,L6)')   '    LinearGrid = ', var%linearGrid
    write(unit = nfprt, fmt = '(A,L6)')   '         GrADS = ', varCommon%grads
    write(unit = nfprt, fmt = '(A)')      '       DateICn = '//varCommon%date
    write(unit = nfprt, fmt = '(A)')      '       Preffix = '//var%preffix
    write(unit = nfprt, fmt = '(A)')      '     dirSSTnamelist = '//trim(var%dirSSTnamelist)
    write(unit = nfprt, fmt = '(A,/)')    ' /'

    call msgInLineFormatOut(' ' // headerMsg, '(A)')
    call msgInLineFormatOut('   Imax = ', '(A)')
    call msgInLineFormatOut(varCommon%xMax, '(I6)')
    call msgInLineFormatOut('   Jmax = ', '(A)')
    call msgInLineFormatOut(varCommon%yMax, '(I6)')
    call msgNewLine();

  end subroutine printNameList
  
 
  subroutine SSTClimatological ()
    !# SSTClimatological
    !# ---
    !# @info
    !# **Brief:** SSTClimatological. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none
 
    ! 1950-1979 1 Degree x 1 Degree SST 
    ! Global NCEP OI Monthly Climatology
    ! Grid Orientation (SSTR):
    ! (1,1) = (0.5_p_r8W,89.5_p_r8N)
    ! (Idim,Jdim) = (0.5_p_r8E,89.5_p_r8S)
 
    integer :: m, monthBefore, monthAfter
 
    real(kind = p_r8) :: dayHour, dayCorrection, FactorBefore, FactorAfter
 
    integer :: Header(8)
 
    real(kind = p_r8) :: sstBefore(var%xDim,var%yDim), sstAfter(var%xDim,var%yDim)
    dayHour=0.0_p_r8
    !dayHour=real(day,p_r8)+real(Hour,p_r8)/24.0_p_r8
    monthBefore=lAbs(2)
    !IF (dayHour > (1.0_p_r8+real(monthLength(month),p_r8)/2.0_p_r8)) &
    !    monthBefore=month
    monthAfter=monthBefore+1
    if (monthBefore < 1) monthBefore=12
    if (monthAfter > 12) monthAfter=1
    dayCorrection=real(monthLength(monthBefore),p_r8)/2.0_p_r8-1.0_p_r8
    !IF (monthBefore == month) dayCorrection=-dayCorrection-2.0_p_r8
    dayCorrection=-dayCorrection-2.0_p_r8
    FactorAfter=2.0_p_r8*(dayHour+dayCorrection)/ &
                real(monthLength(monthBefore)+monthLength(monthAfter),p_r8)
    FactorBefore=1.0_p_r8-FactorAfter
 
    write(unit = nfprt, fmt = '(/,A)') ' From SSTClimatological:'
    write(unit = nfprt, fmt = '(/,A,I4,3(A,I2.2))') &
          ' Year = ', lAbs(2), ' Month = ', lAbs(2), &
          ' Day = ',  lAbs(3), ' Hour = ', hour
    write(unit = nfprt, fmt = '(/,2(A,I2))') &
          ' MonthBefore = ', monthBefore, ' MonthAfter = ', monthAfter
    write(unit = nfprt, fmt = '(/,2(A,F9.6),/)') &
          ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter
 
    nfclm = openFile(trim(varCommon%dirClmSST)//fileClmSST, 'formatted', 'sequential', -1, 'read', 'old')
    if(nfclm < 0) return          
    do m=1,12
       read(unit = nfclm, fmt = '(8I5)') Header
       write(unit = nfprt, fmt = '(/,1X,9I5,/)') m, Header
       read(unit = nfclm, fmt = '(16F5.2)') sstClim
       if (m == monthBefore) then
          sstBefore=sstClim
       end if
       if (m == monthAfter) then
          sstAfter=sstClim
       end if
    end do
    close(unit = nfclm)
 
    ! Linear Interpolation in Time for Year, Month, Day and Hour 
    ! of the Initial Condition
    sstClim=FactorBefore*sstBefore+FactorAfter*sstAfter
 
  end subroutine SSTClimatological
 
 
  subroutine SSTClimatological_O ()
    !# SSTClimatological_O
    !# ---
    !# @info
    !# **Brief:** SSTClimatological_O. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none
 
    ! 1950-1979 1 Degree x 1 Degree SST 
    ! Global NCEP OI Monthly Climatology
    ! Grid Orientation (SSTR):
    ! (1,1) = (0.5_p_r8W,89.5_p_r8N)
    ! (Idim,Jdim) = (0.5_p_r8E,89.5_p_r8S)
 
    integer :: m, monthBefore, monthAfter
 
    real(kind = p_r8) :: dayHour, dayCorrection, FactorBefore, FactorAfter
 
    integer :: Header(8)
 
    real(kind = p_r8) :: sstBefore(var%xDim,var%yDim), sstAfter(var%xDim,var%yDim)
 
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
 
    write(unit = nfprt, fmt = '(/,A)') ' From SSTClimatological:'
    write(unit = nfprt, fmt = '(/,A,I4,3(A,I2.2))') &
          ' Year = ', year, ' Month = ', month, &
          ' Day = ', day, ' Hour = ', hour
    write(unit = nfprt, fmt = '(/,2(A,I2))') &
          ' MonthBefore = ', monthBefore, ' MonthAfter = ', monthAfter
    write(unit = nfprt, fmt = '(/,2(A,F9.6),/)') &
          ' FactorBefore = ', FactorBefore, ' FactorAfter = ', FactorAfter

    nfclm = openFile(trim(varCommon%dirClmSST)//fileClmSST, 'formatted', 'sequential', -1, 'read', 'old')
    if(nfclm < 0) return      
    do m=1,12
       read(unit = nfclm, fmt = '(8I5)') Header
       write(unit = nfprt, fmt = '(/,1X,9I5,/)') m, Header
       read(unit = nfclm, fmt = '(16F5.2)') sstClim
       if (m == monthBefore) then
          sstBefore=sstClim
       end if
       if (m == monthAfter) then
          sstAfter=sstClim
       end if
    end do
    close(unit = nfclm)
 
    ! Linear Interpolation in Time for Year, Month, Day and Hour 
    ! of the Initial Condition
    sstClim=FactorBefore*sstBefore+FactorAfter*sstAfter
 
  end subroutine SSTClimatological_O
  
  
  subroutine SSTClimaWindow ()
    !# Gets Indices to Use CLimatological SST Out of LatClimSouth to LatClimNorth
    !# ---
    !# @info
    !# **Brief:** Gets Indices to Use CLimatological SST Out of LatClimSouth to
    !# LatClimNorth. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo
    implicit none
 
    integer :: j
 
    real(kind = p_r8) :: lat, dLat
 
    js=0
    jn=0
    dLat=2.0_p_r8*p_lat0/real(var%yDim-1,p_r8)
    do j=1,var%yDim
       lat=p_lat0-real(j-1,p_r8)*dLat
       if (lat > var%latClimSouth) js=j
       if (lat > var%latClimNorth) jn=j
    end do
    js=js+1
 
    write(unit = nfprt, fmt = '(/,A,/)')' From SSTClimaWindow:'
    write(unit = nfprt, fmt = '(A,I3,A,F7.3)') &
          ' js = ', js, ' LatClimSouth=', p_lat0-real(js-1,p_r8)*dLat
    write(unit = nfprt, fmt = '(A,I3,A,F7.3,/)') &
          ' jn = ', jn, ' LatClimNorth=', p_lat0-real(jn-1,p_r8)*dLat
 
  end subroutine SSTClimaWindow
  

end module Mod_SSTMonthlyDirec
