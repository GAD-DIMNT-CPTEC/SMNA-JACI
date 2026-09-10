! Eduardo Khamis - 22-03-2021
module Mod_OCMClima

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

  public :: getNameOCMClima, initOCMClima, generateOCMClima, shouldRunOCMClima

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  ! input parameters

  type OCMClimaNameListData
    integer :: zMax
    integer :: xDim
    integer :: yDim
    real (kind = p_r8) :: ocmSeaIce = -1.749
    logical :: linear = .false.
    logical :: linearGrid = .false.
    character (len = 128) :: nameLSM = 'ModelLandSeaMask'
    character (len = 7) :: preffix = 'GANLSMT'
    character (len = 6) :: suffix = 'S.unf.'
  end type OCMClimaNameListData

  type(varCommonNameListData) :: varCommon
  type(OCMClimaNameListData)  :: var
  namelist /OCMClimaNameList/    var

  integer :: xmx, yMaxHf, &
             mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3

  integer :: year, month, day, hour

  real (kind = p_r8) :: lon0, lat0, p_to, ocmOpenWater, &
                            ocmSeaIceThreshold, lapseRate

  logical :: flagInput(5), flagOutput(5)

  character (len = 10)  :: trunc='T     L   '
  character (len = 7)   :: nLats='.G     '
  character (len = 10)  :: mskfmt = '(      I1)'
  character (len = 8)   :: varName='OCMClima'
  character (len = 11)  :: fileClmOCM='ersst.form'
  character (len = 20)  :: fileClmOCMData='ocm_data.bin'

  integer, parameter :: nMon  = 12
  integer, parameter :: nlevl = 19

  character (len = 3), dimension (12) :: monthChar = &
            (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

  integer :: nfprt = 6 ! Standard print Out
  integer :: nficn     ! To read Topography from Initial Condition
  integer :: nflsm     ! To read Formatted Land Sea Mask
  integer :: nfclm     ! To read Formatted Climatological OCM
  integer :: nfocm     ! To read Formatted Climatological OCM
  integer :: nfsto     ! To write Unformatted Gaussian Grid OCM
  integer :: nfout     ! To write GrADS Topography, Land Sea, Se Ice and Gauss OCM
  integer :: nfoce
  integer :: nfwat
  integer :: nftem
  integer :: nfsal

  ! internal variables

  ! Reads the 1x1 OCM Global Monthly OI Climatology From NCEP,
  ! Interpolates it Using Area Weigth or Bi-Linear Into a Gaussian Grid

  integer :: j, i, m, nr, lRec, ios, it, k, irec, im, ierr

  integer :: forecastDay

  real (kind = p_r4) :: timeOfDay

  real (kind = p_r8) :: rgOCMMax, rgOCMMin, ggOCMMax, ggOCMMin, mgOCMMax, mgOCMMin

  integer :: iCnDate(4), currentDate(4), headerAux(8)

  integer, dimension (:,:), allocatable :: landSeaMask, seaIceMask

  real (kind = p_r4), dimension (:),   allocatable :: coefTopIn

  real (kind = p_r4), dimension (:,:), allocatable :: wrOut
  real (kind = p_r4), dimension (:,:), allocatable :: aux
  integer           , dimension (:,:), allocatable :: auxI

  real (kind = p_r8), dimension   (:), allocatable :: coefTop

  real (kind = p_r8)   , allocatable :: otemp_in      (:,:,:,:)
  real (kind = p_r8)   , allocatable :: salt_in       (:,:,:,:)
  real (kind = p_r8)   , allocatable :: bathy_in      (:,:)
  real (kind = p_r8)   , allocatable :: waterqual_in  (:,:)

  real (kind = p_r8)   , allocatable :: otemp_out     (:,:,:,:)
  real (kind = p_r8)   , allocatable :: salt_out      (:,:,:,:)
  real (kind = p_r8)   , allocatable :: bathy_out     (:,:)
  real (kind = p_r8)   , allocatable :: waterqual_out (:,:)

  real (kind = p_r8), dimension (:,:), allocatable :: topog, ocmIn, ocmGaus, &
                  seaIceFlagIn, seaIceFlagOut

  character (len = 6), dimension (:), allocatable :: ocmLabel

  character (len = maxPathLength) :: inFileName, outFileName 

  character (len = *), parameter :: header = 'OCM Clima          | '



contains



  ! Eduardo Khamis - 22-03-2021        
  function getNameOCMClima() result(returnModuleName)
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "OCMClima"
  end function getNameOCMClima


  ! Eduardo Khamis - 22-03-2021
  function shouldRunOCMClima() result(shouldRun)
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunOCMClima


  ! Eduardo Khamis - 22-03-2021
  function getOutFileName() result(ocmClimaOutFilename)
    implicit none
    character(len = maxPathLength) :: ocmClimaOutFilename

    ocmClimaOutFilename = trim(varCommon%dirModelIn) // trim(varName) // varCommon%date(1:8) // nLats
  end function getOutFileName


  ! Eduardo Khamis - 22-03-2021
  subroutine initOCMClima (nameListFileUnit, varCommon_)
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
    
    read(unit = nameListFileUnit, nml = OCMClimaNameList)
    varCommon = varCommon_

    !call printNameList()
    call getxMaxyMax (varCommon%mEnd, varCommon%xMax, varCommon%yMax)
    write (nLats(3:7), '(I5.5)') varCommon%yMax

    mEnd1 = varCommon%mEnd + 1
    mEnd2 = varCommon%mEnd + 2
    mnwv2 = mEnd1 * mEnd2
    mnwv0 = mnwv2 / 2
    mnwv3 = mnwv2 + 2 * mEnd1
    mnwv1 = mnwv3 / 2

    xmx = varCommon%xMax + 2
    yMaxHf = varCommon%yMax / 2

    p_to = 273.15_p_r8
    ocmOpenWater = -1.7_p_r8 + p_to
    ocmSeaIceThreshold = 271.2_p_r8

    lapseRate = 0.0065_p_r8

    read (varCommon%date, fmt = '(I4,3I2)') year, month, day, hour

    ! For Linear Interpolation
    lon0 = 0.5_p_r8  ! Start Near Greenwhich
    lat0 = 89.5_p_r8 ! Start Near North Pole

    ! For Area Weighted Interpolation
    flagInput(1) = .true.   ! Start at North Pole
    flagInput(2) = .true.   ! Start at Prime Meridian
    flagInput(3) = .false.  ! Latitudes Are at North Edge of Box
    flagInput(4) = .false.  ! Longitudes Are at Western Edge of Box
    flagInput(5) = .false.  ! Regular Grid
    flagOutput(1) = .true.  ! Start at North Pole
    flagOutput(2) = .true.  ! Start at Prime Meridian
    flagOutput(3) = .false. ! Latitudes Are at North Edge of Box
    flagOutput(4) = .true.  ! Longitudes Are at Center of Box
    flagOutput(5) = .true.  ! Gaussian Grid

    if (var%linearGrid) then
      trunc(2:2) = 'L'
    else
      trunc(2:2) = 'Q'
    end if
    write (trunc(3:6),  fmt = '(I4.4)') varCommon%mEnd
    write (trunc(8:10), fmt = '(I3.3)') var%zMax

    write (mskfmt(2:7),  '(I6)') varCommon%xMax

    write (nLats(3:7), '(I5.5)') varCommon%yMax

  end subroutine initOCMClima


  ! Eduardo Khamis - 22-03-2021
  subroutine printNameList()
    implicit none
    
    write (unit = nfprt, fmt = '(/,A)')    ' &OCMClimaNameList'
    write (unit = nfprt, fmt = '(A,I6)')   '       mEnd = ', varCommon%mEnd
    write (unit = nfprt, fmt = '(A,I6)')   '       zMmax = ', var%zMax
    write (unit = nfprt, fmt = '(A,I6)')   '       xDim = ', var%xDim
    write (unit = nfprt, fmt = '(A,I6)')   '       yDim = ', var%yDim
    write (unit = nfprt, fmt = '(A,F6.3)') '  ocmSeaIce = ', var%ocmSeaIce
    write (unit = nfprt, fmt = '(A,L6)')   '     linear = ', var%linear
    write (unit = nfprt, fmt = '(A,L6)')   ' linearGrid = ', var%linearGrid
    write (unit = nfprt, fmt = '(A,L6)')   '      grads = ', varCommon%grads
    write (unit = nfprt, fmt = '(A)')      '    dateICn = ' // varCommon%date
    write (unit = nfprt, fmt = '(A)')      '    preffix = ' // var%preffix
    write (unit = nfprt, fmt = '(A)')      '    suffix = ' // var%suffix
    write (unit = nfprt, fmt = '(A)')      '    dirPreOut = ' // trim(varCommon%dirPreOut)
    write (unit = nfprt, fmt = '(A)')      '    dirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = nfprt, fmt = '(A)')      '    dirClmSST = ' // trim(varCommon%dirClmSST)
    write (unit = nfprt, fmt = '(A,/)')    ' /'
  end subroutine printNameList


  function readLandSeaMask() result(isReadOk)
    implicit none
    logical :: isReadOk

    isReadOk = .false.

    inFileName = trim(varCommon%dirPreOut)//trim(var%nameLSM)//nLats
    nflsm = openFile(trim(inFileName), 'unformatted', 'sequential', -1, 'read', 'old')
    if (nflsm < 0) return
    read(unit=nflsm) landSeaMask
    close(unit=nflsm)

    isReadOk = .true.

  end function readLandSeaMask


! Eduardo Khamis - 22-03-2021
   function generateOCMClima() result(isExecOk)
      implicit none
      logical :: isExecOk

      isExecOk = .false.

      write (unit = nfprt, fmt = '(/,A,I6)') '   Imax = ', varCommon%xMax
      write (unit = nfprt, fmt = '(A,I6,/)') '   Jmax = ', varCommon%yMax
   
      call createSpectralRep (mEnd1, mEnd2, mnwv1)
      call createGaussRep (varCommon%yMax, yMaxHf)
      call createFFT (varCommon%xMax)
      call createLegTrans (mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)
   
      if (var%linear) then
         call initLinearInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, lat0, lon0)
      else
         call initAreaInterpolation (var%xDim, var%yDim, varCommon%xMax, varCommon%yMax, flagInput, flagOutput)
      end if
   
      call allocateData()

      inquire (iolength = lRec) aux(1:var%xDim, 1:var%yDim)
      inFileName = trim(varCommon%dirClmOCM)//trim(fileClmOCMData)
      nfocm = openFile(trim(inFileName), 'unformatted', 'direct', lRec, 'read', 'old')
      if (nfocm < 0) return   
   
   
      irec=0
      do it=1, nMon
         do k=1,nlevl
            irec=irec+1
            read(nfocm,rec=irec)aux
            otemp_in(1:var%xDim,1:var%yDim,k,it)=aux(1:var%xDim,1:var%yDim)
         end do
         do k=1,nlevl
            irec=irec+1
            read(nfocm,rec=irec)aux!
            salt_in(1:var%xDim,1:var%yDim,k,it)=aux(1:var%xDim,1:var%yDim)
         end do
         irec=irec+1
         read(nfocm,rec=irec)aux
          bathy_in=aux
         irec=irec+1
         read(nfocm,rec=irec)aux
          waterqual_in=aux
      end do
      irec=0
      ! Read in SpectraL Coefficient of Topography from ICn
      ! to Ensure that Topography is the Same as Used by Model
      inFileName = trim(varCommon%dirModelIn) // var%preffix // varCommon%date // var%suffix // trunc
      nficn = openFile(trim(inFileName), 'unformatted', 'sequential', -1, 'read', 'old') 
      if(nficn < 0) return
      read  (unit = nficn) forecastDay, timeOfDay, iCnDate, currentDate
      read  (unit = nficn) coefTopIn
      close (unit = nficn)     
      coefTop=real(coefTopIn,p_r8)
      call transp(mnwv2, mEnd1, mEnd2, coefTop)
      call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)
   
      ! Read in Land-Sea Mask Data Set
      if (.not. readlandSeaMask()) then
        call msgWarningOut(header, "Error reading LandSeaMask file")
        return
      end if
   
      ! Open File for Land-Sea Mask and OCM Data to Global Model Input
      inquire (iolength = lRec) wrOut
      outFileName = trim(varCommon%dirModelIn)//trim(varName)//varCommon%date(1:8)//nLats 
      nfsto = openFile(trim(outFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
      if (nfsto < 0) return
   
      ! Write out Land-Sea Mask to OCM Data Set
      ! The LSMask will be Transfered by Model to Post-Processing
   
      wrOut=real(1-2*landSeaMask,p_r4)
      write (unit = nfsto, rec = 1) wrOut
   
      if (varCommon%grads) then
         outFileName = trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats
         nfout = openFile(trim(outFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
         if (nfout < 0) return
      end if
   
      inFileName = trim(varCommon%dirClmOCM)//fileClmOCM
      nfclm = openFile(trim(inFileName), 'formatted', 'sequential', -1, 'read', 'old')
      if (nfclm < 0) return
   
      ! Loop Through Months
      irec=0
      do m=1,12
         print*,'month=',m
         read (unit = nfclm, fmt = '(8I5)') headerAux
         write (unit = nfprt, fmt = '(/,1X,9I5,/)') m, headerAux
         read (unit = nfclm, fmt = '(16F5.2)') ocmIn
         if (maxval(ocmIn) < 100.0_p_r8) ocmIn=ocmIn+p_to
   
         ! Convert Sea Ice into SeaIceFlagIn = 1 and Over Ice, = 0
         ! Over Open Water Set Input OCM = MIN of OCMOpenWater
         ! Over Non Ice Points Before Interpolation
         if (var%ocmSeaIce < 100.0_p_r8) var%ocmSeaIce=var%ocmSeaIce+p_to
         do j=1,var%yDim
            do i=1,var%xDim
               seaIceFlagIn(i,j)=0.0_p_r8
               if (ocmIn(i,j) < var%ocmSeaIce) then
                  seaIceFlagIn(i,j)=1.0_p_r8
               else
                  ocmIn(i,j)=max(ocmIn(i,j),ocmOpenWater)
               end if
            end do
         end do
         ! Min And Max Values of Input OCM
         rgOCMMax=maxval(ocmIn)
         rgOCMMin=minval(ocmIn)
   
         ! Interpolate Flag from 1x1 Grid to Gaussian Grid, Fill SeaIceMask=1
         ! Over Interpolated Points With 50% or More Sea Ice, =0 Otherwise
         print*,'SeaIceFlagOut='
   
         if (var%linear) then
            call doLinearInterpolation (seaIceFlagIn, seaIceFlagOut)
         else
            call doAreaInterpolation (seaIceFlagIn, seaIceFlagOut)
         end if
         seaIceMask=int(seaIceFlagOut+0.5_p_r8)
         where (landSeaMask == 1) seaIceMask=0
         print*,'OCMIn='
   
         ! Interpolate OCM from 1x1 Grid to Gaussian Grid
         if (var%linear) then
            call doLinearInterpolation (ocmIn, ocmGaus)
         else
            call doAreaInterpolation (ocmIn, ocmGaus)
         end if
          print*,'otemp_out='
   
   
          ! Interpolate OCM from 1x1 Grid to Gaussian Grid
         do k=1,nlevl
            if (var%linear) then
               call doLinearInterpolation (otemp_in(:,:,k,m), otemp_out(:,:,k,m))
            else
               call doAreaInterpolation   (otemp_in(:,:,k,m), otemp_out(:,:,k,m))
            end if
         end do
         print*,'salt_in='
   
         ! Interpolate OCM from 1x1 Grid to Gaussian Grid
         do k=1,nlevl
            if (var%linear) then
               call doLinearInterpolation (salt_in(:,:,k,m), salt_out(:,:,k,m))
            else
               call doAreaInterpolation   (salt_in(:,:,k,m), salt_out(:,:,k,m))
            end if
         end do
         print*,'bathy_in='
   
         ! Interpolate OCM from 1x1 Grid to Gaussian Grid
         if (var%linear) then
            call doLinearInterpolation (bathy_in, bathy_out)
         else
            call doAreaInterpolation (bathy_in, bathy_out)
         end if
         
         ! Interpolate OCM from 1x1 Grid to Gaussian Grid
         print*,'waterqual_in='
   
         if (var%linear) then
            call doLinearInterpolation (waterqual_in, waterqual_out)
         else
            call doAreaInterpolation (waterqual_in, waterqual_out)
         end if
   
         ! Min and Max Values of Gaussian Grid
         ggOCMMax=maxval(ocmGaus)
         ggOCMMin=minval(ocmGaus)
         print*,'Min and Max Values of Gaussian Grid='
   
         do j=1,varCommon%yMax
            do i=1,varCommon%xMax
               if (landSeaMask(i,j) == 1) then
                  ! Set OCM = Undef Over Land
                  ocmGaus(i,j)=p_undef
                  waterqual_out(i,j) = 0
                  if(topog(i,j) <=0.0)then
                     bathy_out(i,j)=1
                  else
                     bathy_out(i,j)=topog(i,j)
                  end if
               else if (seaIceMask(i,j) == 1) then
                  ! Set OCM Sea Ice Threshold Minus 1 Over Sea Ice
                  bathy_out(i,j)=bathy_out(i,j)
                  ocmGaus(i,j)=ocmSeaIceThreshold-1.0_p_r8
                  waterqual_out(i,j) = int(waterqual_out(i,j))
               else
                  bathy_out(i,j)=bathy_out(i,j)
                  waterqual_out(i,j) = int(waterqual_out(i,j))
                  ! Correct OCM for Topography, Do Not Create or
                  ! Destroy Sea Ice Via Topography Correction
                  ocmGaus(i,j)=ocmGaus(i,j)-topog(i,j)*lapseRate
                  if (ocmGaus(i,j) < ocmSeaIceThreshold) &
                     ocmGaus(i,j)=ocmSeaIceThreshold+0.2_p_r8
               end if
            end do
         end do
         print*,'Min and Max Values of Corrected Gaussian Grid OCM Excluding Land Points'
   
         ! Min and Max Values of Corrected Gaussian Grid OCM Excluding Land Points
         mgOCMMax=maxval(ocmGaus,MASK=ocmGaus/=p_undef)
         mgOCMMin=minval(ocmGaus,MASK=ocmGaus/=p_undef)
   
         ! Write out Gaussian Grid Weekly OCM
         wrOut=real(ocmGaus,p_r4)
         write (unit = nfsto, rec = m+1) wrOut
   
         write (unit = nfprt, fmt = '(/,3(A,I2.2),A,I4)') &
               ' Hour = ', hour, ' Day = ', day, &
               ' Month = ', month, ' Year = ', year
   
         write (unit = nfprt, fmt = '(/,A,3(A,2F8.2,/))') &
               ' Mean Weekly OCM Interpolation :', &
               ' Regular  Grid OCM: Min, Max = ', rgOCMMin, rgOCMMax, &
               ' Gaussian Grid OCM: Min, Max = ', ggOCMMin, ggOCMMax, &
               ' Masked G Grid OCM: Min, Max = ', mgOCMMin, mgOCMMax
   
         if (varCommon%grads) then
            !         nr=1+4*(m-1)
            wrOut=real(topog,p_r4)
            irec=irec+1
            write (unit = nfout, rec = irec) wrOut
            wrOut=real(1-2*landSeaMask,p_r4)
            irec=irec+1
            write (unit = nfout, rec = irec) wrOut
            wrOut=real(seaIceMask,p_r4)
            irec=irec+1
            write (unit = nfout, rec = irec) wrOut
            wrOut=real(ocmGaus,p_r4)
            irec=irec+1
            write (unit = nfout, rec = irec) wrOut
            wrOut=real(bathy_out,p_r4)
            irec=irec+1
            write (unit = nfout, rec = irec) wrOut
            wrOut=real(waterqual_out,p_r4)
            irec=irec+1
            write (unit = nfout, rec = irec) wrOut
            do k=1,nlevl
               wrOut=real(otemp_out(:,:,nlevl-k+1,m),p_r4)+p_to
               irec=irec+1
               write(nfout,rec = irec) wrOut
            end do
            do k=1,nlevl
               wrOut=real(salt_out(:,:,nlevl-k+1,m),p_r4)
               irec=irec+1
               write(nfout,rec = irec) wrOut
            end do
         end if
   
      ! End Loop Through Months
      end do
   
   
   
   
      close (unit = nfclm)
      close (unit = nfsto)
      close (unit = nfout)
   
      !
      !   read in ocean bathymetry (<0 over ocean; >=0 over land) 
      !   and determines the bottom of the model at each point, 
      !   lbottom indicates the first inactive layer
      inquire (iolength = lRec) wrOut(1:varCommon%xMax, 1:varCommon%yMax)
      outFileName = trim(varCommon%dirModelIn)//'/'//'ocean_depth'//nLats 
      nfoce = openFile(trim(outFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
      if (nfoce < 0) return
      irec = 0
      irec = irec + 1
      wrOut = real(bathy_out, p_r4)
      print *, 'aaa', maxval(wrOut), minval(wrOut)
      write(nfoce, rec = 1) wrOut
      close (nfoce)
      !   global annual average optical water type from the map of 
      !   siminot and le treut (1986,jgr).
      !   water types  -   numerical value in file:
      !    land               0
      !    i                  1
      !    ii                 2
      !    iii                3
      !    ia                 4
      !    ib                 5
      inquire (iolength = lRec) wrOut(1:varCommon%xMax, 1:varCommon%yMax)
      outFileName = trim(varCommon%dirModelIn)//'/'//'water_type'//nLats
      nfwat = openFile(trim(outFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
      if (nfwat < 0) return
      irec = 0
      irec = irec + 1
      auxI = int(waterqual_out)
      print *, 'aaa', maxval(auxI), minval(auxI)
      write (nfwat,rec = irec) auxI
      close (nfwat)
   
   
   !  observed monthly ltm mean temp and salinity 
   !  will be used in relaxation              
      inquire (iolength = lRec) wrOut(1:varCommon%xMax, 1:varCommon%yMax)
      outFileName = trim(varCommon%dirModelIn)//'/'//'temp_ltm_month'//nLats
      nftem = openFile(trim(outFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
      if (nftem < 0) return
      irec = 0
      do im = 1, 12
        do k = 1, nlevl
            irec = irec + 1
            wrOut = real(otemp_out(:, :, nlevl - k + 1, im), p_r4) + p_to
            write(nftem, rec = irec) wrOut
         end do
      end do
      close (nftem)
   
   !  observed monthly ltm mean temp and salinity 
   !  will be used in relaxation              
      inquire (iolength = lRec) wrOut(1:varCommon%xMax, 1:varCommon%yMax)
      outFileName = trim(varCommon%dirModelIn)//'/'//'salt_ltm_month'//nLats
      nfsal = openFile(trim(outFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
      if (nfsal < 0) return
      irec = 0
      do im =  1, 12
        do k = 1, nlevl
            irec = irec + 1
            wrOut = real(salt_out(:, :, nlevl - k + 1, im), p_r4)
            write(nfsal, rec = irec) wrOut
         end do
      end do
      close (nfsal)
   
     if (varCommon%grads .and. .not. generateGrads()) then
       call msgWarningOut(header, "Error while generating grads files") 
       return
     end if

     call deallocateData()
     call destroyFFT()
     call destroyLegendreObjects()

     isExecOk = .true.

   end function generateOCMClima


   !  Eduardo Khamis - 23-03-2021
   function generateGrads() result(isGradsOK)
      implicit none
      logical :: isGradsOk
      integer :: gradsFileUnit

      isGradsOk = .false.

      ! Write GrADS Control File
      outFileName = trim(varCommon%dirPreOut)//trim(varName)//varCommon%date(1:8)//nLats//'.ctl'
      gradsFileUnit = openFile(trim(outFileName), 'formatted', 'sequential', -1, 'write', 'replace')
      if (gradsFileUnit < 0) return

      write (unit = gradsFileUnit, fmt = '(A)') 'DSET ^'//trim(varName)//varCommon%date(1:8)//nLats
      write (unit = gradsFileUnit, fmt = '(A)') '*'
      write (unit = gradsFileUnit, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
      write (unit = gradsFileUnit, fmt = '(A)') '*'
      write (unit = gradsFileUnit, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
      write (unit = gradsFileUnit, fmt = '(A)') '*'
      write (unit = gradsFileUnit, fmt = '(A)') 'TITLE Monthly Climatological OI OCM on a Gaussian Grid'
      write (unit = gradsFileUnit, fmt = '(A)') '*'
      write (unit = gradsFileUnit, fmt = '(A,I5,A,F8.3,F15.10)') &
                          'XDEF ',varCommon%xMax,' LINEAR ',0.0_p_r8,360.0_p_r8/real(varCommon%xMax,p_r8)
      write (unit = gradsFileUnit, fmt = '(A,I5,A)') 'YDEF ',varCommon%yMax,' LEVELS '
      if (var%linear) then
         write (unit = gradsFileUnit, fmt = '(8F10.5)') gLatsL(varCommon%yMax:1:-1)
      else
         write (unit = gradsFileUnit, fmt = '(8F10.5)') gLatsA(varCommon%yMax:1:-1)
      end if
      write (unit = gradsFileUnit, fmt = '(A)') 'ZDEF    19 LEVELS  1000 900 800 700 600 500'
      write (unit = gradsFileUnit, fmt = '(A)') '                    400 300 250 200 150 125'
      write (unit = gradsFileUnit, fmt = '(A)') '                    100  75  50  30  20  10 0'
   
      write (unit = gradsFileUnit, fmt = '(A)') 'TDEF 12 LINEAR JAN2007 1MO'
      write (unit = gradsFileUnit, fmt = '(A)') 'VARS  8'
      write (unit = gradsFileUnit, fmt = '(A)') 'TOPO  0 99 Model Recomposed Topography [m]'
      write (unit = gradsFileUnit, fmt = '(A)') 'LSMK  0 99 Land Sea Mask    [1:Sea -1:Land]'
      write (unit = gradsFileUnit, fmt = '(A)') 'SIMK  0 99 Sea Ice Mask     [1:SeaIce 0:NoIce]'
      write (unit = gradsFileUnit, fmt = '(A)') 'OCMC  0 99 Climatological OCM Topography Corrected [K]'
      write (unit = gradsFileUnit, fmt = '(A)') 'BATM  0 99 bathymetric [m]'
      write (unit = gradsFileUnit, fmt = '(A)') 'WAQL  0 99 waterqual   [0-5]'
      write (unit = gradsFileUnit, fmt = '(A)') 'TEMP 19 99 otemp       [C]'
      write (unit = gradsFileUnit, fmt = '(A)') 'SALT 19 99 salinity    [g/kg]'
      write (unit = gradsFileUnit, fmt = '(A)') 'ENDVARS'
   
      close (unit = gradsFileUnit)

      isGradsOk = .true.

   end function generateGrads


   !  Eduardo Khamis - 23-03-2021
   subroutine allocateData()
      implicit none

      allocate (landSeaMask(varCommon%xMax,varCommon%yMax), seaIceMask(varCommon%xMax,varCommon%yMax))
      allocate (coefTopIn(mnwv2))
      allocate (coefTop(mnwv2), topog(varCommon%xMax,varCommon%yMax))
      allocate (ocmIn(var%xDim,var%yDim), ocmGaus(varCommon%xMax,varCommon%yMax), wrOut(varCommon%xMax,varCommon%yMax))
      allocate (seaIceFlagIn(var%xDim,var%yDim), seaIceFlagOut(varCommon%xMax,varCommon%yMax))
      allocate (ocmLabel(var%yDim))
      allocate (aux (var%xDim,var%yDim))
      allocate (otemp_in        (var%xDim,var%yDim,nlevl,nMon))
      allocate (salt_in         (var%xDim,var%yDim,nlevl,nMon))
      allocate (bathy_in        (var%xDim,var%yDim))
      allocate (waterqual_in    (var%xDim,var%yDim))
   
      allocate(otemp_out       (varCommon%xMax,varCommon%yMax,nlevl,nMon))
      allocate(salt_out        (varCommon%xMax,varCommon%yMax,nlevl,nMon))
      allocate(bathy_out       (varCommon%xMax,varCommon%yMax))
      allocate(waterqual_out   (varCommon%xMax,varCommon%yMax))
      allocate(auxI  (varCommon%xMax,varCommon%yMax))

   end subroutine allocateData


   !  Eduardo Khamis - 23-03-2021
   subroutine deallocateData()
      implicit none

      deallocate (landSeaMask, seaIceMask     )
      deallocate (coefTopIn                   )
      deallocate (coefTop, topog              )
      deallocate (ocmIn, ocmGaus, wrOut       )
      deallocate (seaIceFlagIn, seaIceFlagOut )
      deallocate (ocmLabel       )
      deallocate (aux            )
      deallocate (otemp_in       )
      deallocate (salt_in        )
      deallocate (bathy_in       )
      deallocate (waterqual_in   )
   
      deallocate (otemp_out      )
      deallocate (salt_out       )
      deallocate (bathy_out      )
      deallocate (waterqual_out  )
      deallocate (auxI           )

   end subroutine deallocateData
   

end module Mod_OCMClima
