!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_AreaInterpolation </br></br>
!#
!# **Brief**: Module responsible for determining the area weights and index 
!# mapping between a pair of grids on a sphere. The weights and map indices are
!# used to perfom the actual interpolation. After interpolation user is 
!# responsible for output masking and pole interpolation for type 2 
!# (pole centered) grids. </br></br>
!#
!# **Files in:**
!#
!# &bull;  ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ?
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.2 <br><br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti  - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita           - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel Lamosa    - version: 2.0.0 - Module creation </li>
!#  <li>26-04-2019 - Denis Eiras      - version: 2.0.1 </li>
!#  <li>03-06-2019 - Eduardo Khamis   - version: 2.0.2 </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time </li>
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


module Mod_AreaInterpolation

  implicit none
  private

  public :: initAreaInterpolation, doAreaInterpolation, gLats

  !parameters
  include 'precision.h'

  !input variables
  logical :: polarMean    
  !# Flag to performe average at poles 
  integer :: xDim        
  !# Number of longitude points for the input grid  
  integer :: yDim       
  !# Number of latitude points for the input grid  
  integer :: xMax       
  !# Number of longitude points for the output grid    
  integer :: yMax          
  !# number of latitude points for the output grid

  ! flags: (input or output)
  !   1   start at north pole (true) start at south pole (false)
  !   2   start at prime meridian (true) start at i.d.l. (false)
  !   3   latitudes are at center of box (true)
  !       latitudes are at edge (false) north edge if 1=true
  !                                     south edge if 1=false
  !   4   longitudes are at center of box (true)
  !       longitudes are at western edge of box (false)
  !   5   gaussian (true) regular (false)
  
  logical :: flagInput(5)  
  !# Input  grid flags
  logical :: flagOutput(5) 
  !# Output grid flags
   
  !# input grid mask to confine interpolation of input data to certain areas (1=interpolate, 0=don't)
  integer, dimension(:, :), allocatable :: maskInput   
  
  integer :: lons 
  !# longitude dimension for weights
  integer :: lats 
  !# latitude dimension for weights
  integer :: lond 
  !# total number of longitude weights
  integer :: latd 
  !# total number of latitude weigths
  integer :: lwrk 
  !# dimension for working area
  real(kind = p_r8), dimension(:), allocatable :: gLats 
  !# ?

  integer, dimension(:, :), allocatable :: mplon 
  !# longitude index mapping from input (,1) to output (,2)
  integer, dimension(:, :), allocatable :: mplat 
  !# latitude index mapping from input (,1) to output (,2)
  real(kind = p_r8), dimension(:), allocatable :: wtlon 
  !# area weights in the longitudinal direction
  real(kind = p_r8), dimension(:), allocatable :: wtlat 
  !# area weights in the latitudianl direction

contains


  subroutine initAreaInterpolation(xDim1, yDim1, xMax1, yMax1, flagInput_, flagOutput_)
    !# Initializes Area Interpolation
    !# ---
    !# @info
    !# **Brief:** Initializes Area Interpolation </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    
    integer, intent(in) :: xDim1 
    !# xDim
    integer, intent(in) :: yDim1 
    !# yDim
    integer, intent(in) :: xMax1 
    !# xMax
    integer, intent(in) :: yMax1 
    !# yMax
    logical, intent(in) :: flagInput_(5)  
    !# Input  grid flags
    logical, intent(in) :: flagOutput_(5) 
    !# Output grid flags
    
    integer :: idx

    xDim = xDim1
    yDim = yDim1
    xMax = xMax1
    yMax = yMax1

    !Load parameters
    if(allocated(maskInput)) deallocate(maskInput)
    allocate(maskInput(xDim, yDim))
    maskInput = 1
    ! always equals 1, for all pre modules until now ...

    polarMean = .false.
    ! always false for all pre modules until now ...

    ! flags: (input or output)
    !   1   start at north pole (true) start at south pole (false)
    !   2   start at prime meridian (true) start at i.d.l. (false)
    !   3   latitudes are at center of box (true)
    !       latitudes are at edge (false) north edge if 1=true
    !                                     south edge if 1=false
    !   4   longitudes are at center of box (true)
    !       longitudes are at western edge of box (false)
    !   5   gaussian (true) regular (false)

    do idx = 1, 5
      flagInput(idx) = flagInput_(idx)
      flagOutput(idx) = flagOutput_(idx)
    enddo

    lons = xDim + xMax + 2
    lats = yDim + yMax + 2
    lwrk = max(2 * lons, 2 * lats)

    if(allocated(mplon)) deallocate(mplon)
    if(allocated(mplat)) deallocate(mplat)
    allocate(mplon(lons, 2))
    allocate(mplat(lats, 2))
    if (allocated(wtlon)) deallocate(wtlon)
    if (allocated(wtlat)) deallocate(wtlat)
    allocate(wtlon(lons))
    allocate(wtlat(lats))

    call getAreaInterpolationWeights()
    ! Computing Weights For Horizontal Interpolation
    
  end subroutine initAreaInterpolation


  subroutine doAreaInterpolation (fieldInput, fieldOutput)
    !# Do Area Interpolation
    !# ---
    !# @info
    !# **Brief:** Should only be called after Subroutine GetAreaInterpolationWeights
    !# has been called. </br>
    !# After interpolation user is responsible for output masking and pole 
    !# interpolation for type 2 (pole centered) grids. Pole interpolation can be
    !# done with polarMean=.true. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo

    real(kind = p_r8), dimension(:, :), intent(in) :: fieldInput  
    !# input field (xDim, yDim) to be rank interpolated
    real(kind = p_r8), dimension(:, :), intent(out) :: fieldOutput 
    !# output field (xMax, yMax) resulting from rank interpolation
    
    ! auxiliary variables
    integer :: j     
    integer :: i     
    integer :: lti   
    integer :: lni   
    integer :: lto   
    integer :: lno   
    real(kind = p_r8) :: wlt 
    real(kind = p_r8) :: wln 
    real(kind = p_r8) :: work(xMax, yMax) 

    fieldOutput = 0.0_p_r8
    work = 0.0_p_r8

    do j = 1, latd
      wlt = wtlat(j)
      lti = mplat(j, 1)
      lto = mplat(j, 2)
      do i = 1, lond
        lni = mplon(i, 1)
        if(maskInput(lni, lti) /= 0) then
          if(fieldInput(lni, lti) /= p_undef) then
            wln = wtlon(i)
            lno = mplon(i, 2)
            fieldOutput(lno, lto) = fieldOutput(lno, lto) + fieldInput(lni, lti) * wlt * wln
            work(lno, lto) = work(lno, lto) + wlt * wln
          end if
        end if
      end do
    end do

    do j = 1, yMax
      do i = 1, xMax
        if(work(i, j) == 0.0_p_r8) then
          fieldOutput(i, j) = p_undef
        else
          fieldOutput(i, j) = fieldOutput(i, j) / work(i, j)
        end if
      end do
    end do

    if(polarMean) then
      fieldOutput(:, 1) = sum(fieldOutput(:, 1)) / real(xMax, p_r8)
      fieldOutput(:, yMax) = sum(fieldOutput(:, yMax)) / real(xMax, p_r8)
    end if

  end subroutine doAreaInterpolation


  subroutine getAreaInterpolationWeights()
    !# Gets the Interpolation Weight Calculation
    !# ---
    !# @info
    !# **Brief:** This Subroutine should be called once to determine the area
    !# weights and index mapping between a pair of grids on a sphere. </br>
    !# The weights and map indices are used by subroutine doAreaInterpolation
    !# to perfom the actual interpolation. </br>
    !# **Authors**: </br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo

    ! auxiliary variables
    integer :: latH 
    integer :: joi  
    integer :: joo  
    integer :: j    
    !# loop iterator
    integer :: j1   
    integer :: j2   
    integer :: j3   
    integer :: ioi  
    integer :: ici  
    integer :: i    
    !# loop iterator
    integer :: ioo  
    integer :: ico  
    integer :: i1   
    integer :: i2   
    integer :: i3   

    ! parameters
    real(kind = p_r4), parameter :: p_eps = epsilon(1.0_p_r4)         
    real(kind = p_r8), parameter :: p_dpi = 4.0_p_r8 * atan(1.0_p_r8) 
    real(kind = p_r8) :: rad = 180.0_p_r8 / p_dpi        

    ! local variables
    real(kind = p_r8) :: workWeights(lwrk) 
    real(kind = p_r8) :: drltm             
    real(kind = p_r8) :: drltp             
    real(kind = p_r8) :: dlat              
    real(kind = p_r8) :: dof               
    real(kind = p_r8) :: delrdi            
    real(kind = p_r8) :: delrdo            


    ! input grid latitudes
    joi = yDim + yMax + 2
    if(flagInput(5)) then

      ! gaussian grid case
      latH = yDim / 2
      call gaussianLatitudes(latH, workWeights)
      do j = 2, yDim
        if(j <= latH) then
          drltm = -p_dpi / 2.0_p_r8 + workWeights(j - 1)
          drltp = -p_dpi / 2.0_p_r8 + workWeights(j)
        else if(j > latH + 1) then
          drltm = p_dpi / 2.0_p_r8 - workWeights(yDim - j + 2)
          drltp = p_dpi / 2.0_p_r8 - workWeights(yDim - j + 1)
        else
          drltm = 0.0_p_r8
          drltp = 0.0_p_r8
        end if
        workWeights(j + joi) = sin((drltm + drltp) / 2.0_p_r8)
      end do
      workWeights(1 + joi) = -1.0_p_r8
      workWeights(latH + 1 + joi) = 0.0_p_r8
      workWeights(yDim + 1 + joi) = 1.0_p_r8
    else

      ! regular grid case
      if(flagInput(3)) then
        dlat = p_dpi / real(yDim - 1, p_r8)
        dof = -(p_dpi + dlat) / 2.0_p_r8
      else
        dlat = p_dpi / real(yDim, p_r8)
        dof = -p_dpi / 2.0_p_r8
      end if

      do j = 2, yDim
        workWeights(joi + j) = sin(dof + dlat * real(j - 1, p_r8))
      end do
      workWeights(1 + joi) = -1.0_p_r8
      workWeights(yDim + 1 + joi) = 1.0_p_r8
    end if

    ! output grid latitudes
    joo = 2 * yDim + yMax + 3
    if(flagOutput(5)) then

      ! gaussian grid case
      latH = yMax / 2
      call gaussianLatitudes(latH, workWeights)
      if(allocated(gLats)) deallocate(gLats)
      allocate (gLats(yMax))
      do j = 1, latH
        gLats(j) = 90.0_p_r8 - rad * workWeights(j)
        gLats(yMax - j + 1) = -gLats(j)
      end do
      do j = 2, yMax
        if (j <= latH) then
          drltm = -p_dpi / 2.0_p_r8 + workWeights(j - 1)
          drltp = -p_dpi / 2.0_p_r8 + workWeights(j)
        else if (j > latH + 1) then
          drltm = p_dpi / 2.0_p_r8 - workWeights(yMax - j + 2)
          drltp = p_dpi / 2.0_p_r8 - workWeights(yMax - j + 1)
        else
          drltm = 0.0_p_r8
          drltp = 0.0_p_r8
        end if
        workWeights(j + joo) = sin((drltm + drltp) / 2.0_p_r8)
      end do
      workWeights(1 + joo) = -1.0_p_r8
      workWeights(latH + 1 + joo) = 0.0_p_r8
      workWeights(yMax + 1 + joo) = 1.0_p_r8
    else

      ! regular grid case
      if(flagOutput(3)) then
        dlat = p_dpi / real(yMax - 1, p_r8)
        dof = -(p_dpi + dlat) / 2.0_p_r8
      else
        dlat = p_dpi / real(yMax, p_r8)
        dof = -p_dpi / 2.0_p_r8
      end if
      if(allocated(gLats)) deallocate(gLats)
      allocate(gLats(yMax))
      do j = 1, yMax
        gLats(j) = 90.0_p_r8 - rad * dlat * real(j - 1, p_r8)
      end do
      do j = 2, yMax
        workWeights(joo + j) = sin(dof + dlat * real(j - 1, p_r8))
      end do
      workWeights(1 + joo) = -1.0_p_r8
      workWeights(yMax + 1 + joo) = 1.0_p_r8
    end if

    ! produce single ordered set of sin(lat) for both grids
    ! determine latitude weighting and index mapping
    j1 = 1
    j2 = 1
    j3 = 1
    do
      if(abs(workWeights(j1 + joi) - workWeights(j2 + joo)) < p_eps) then
        workWeights(j3) = workWeights(j1 + joi)
        if(j3 /= 1) then
          wtlat(j3 - 1) = workWeights(j3) - workWeights(j3 - 1)
          mplat(j3 - 1, 1) = j1 - 1
          if(flagInput(1)) mplat(j3 - 1, 1) = yDim + 2 - j1
          mplat(j3 - 1, 2) = j2 - 1
          if(flagOutput(1)) mplat(j3 - 1, 2) = yMax + 2 - j2
        end if
        j1 = j1 + 1
        j2 = j2 + 1
        j3 = j3 + 1
      else if(workWeights(j1 + joi) < workWeights(j2 + joo)) then
        workWeights(j3) = workWeights(j1 + joi)
        if(j3 /= 1) then
          wtlat(j3 - 1) = workWeights(j3) - workWeights(j3 - 1)
          mplat(j3 - 1, 1) = j1 - 1
          if(flagInput(1)) mplat(j3 - 1, 1) = yDim + 2 - j1
          mplat(j3 - 1, 2) = j2 - 1
          if (flagOutput(1)) mplat(j3 - 1, 2) = yMax + 2 - j2
        end if
        j1 = j1 + 1
        j3 = j3 + 1
      else
        workWeights(j3) = workWeights(j2 + joo)
        if(j3 /= 1) then
          wtlat(j3 - 1) = workWeights(j3) - workWeights(j3 - 1)
          mplat(j3 - 1, 1) = j1 - 1
          if(flagInput(1)) mplat(j3 - 1, 1) = yDim + 2 - j1
          mplat(j3 - 1, 2) = j2 - 1
          if(flagOutput(1)) mplat(j3 - 1, 2) = yMax + 2 - j2
        end if
        j2 = j2 + 1
        j3 = j3 + 1
      end if
      if(.not.(j1 <= yDim + 1 .and. j2 <= yMax + 1)) exit
    end do
    latd = j3 - 2

    ! latitudes done, now do longitudes
    ! input grid longitudes
    ioi = xDim + xMax + 2
    delrdi = (2.0_p_r8 * p_dpi) / real(xDim, p_r8)
    if(flagInput(5) .or. flagInput(4)) then
      ici = 0
      dof = 0.5_p_r8
    else
      ici = 1
      dof = 0.0_p_r8
    end if
    do i = 1, xDim
      workWeights(i + ioi) = (dof + real(i - 1, p_r8)) * delrdi
    end do

    ! output grid longitudes
    ioo = 2 * xDim + xMax + 3
    delrdo = (2.0_p_r8 * p_dpi) / real(xMax, p_r8)
    if(flagOutput(5) .or. flagOutput(4)) then
      ico = 0
      dof = 0.5_p_r8
    else
      ico = 1
      dof = 0.0_p_r8
    end if

    do i = 1, xMax
      workWeights(i + ioo) = (dof + real(i - 1, p_r8)) * delrdo
    end do

    ! produce single ordered set of longitudes for both grids
    ! determine longitude weighting and index mapping
    i1 = 1
    i2 = 1
    i3 = 1
    do
      if(abs(workWeights(i1 + ioi) - workWeights(i2 + ioo)) < p_eps) then
        workWeights(i3) = workWeights(i1 + ioi)
        if (i3 /= 1) then
          wtlon(i3 - 1) = workWeights(i3) - workWeights(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if(.not.flagInput(2)) then
            mplon(i3 - 1, 1) = xDim / 2 + i1 - ici
            if(i1 - ici > xDim / 2) mplon(i3 - 1, 1) = i1 - ici - xDim / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if(.not.flagOutput(2)) then
            mplon(i3 - 1, 2) = xMax / 2 + i2 - ico
            if(i2 - ico > xMax / 2) mplon(i3 - 1, 2) = i2 - ico - xMax / 2
          end if
        end if
        i1 = i1 + 1
        i2 = i2 + 1
        i3 = i3 + 1
      else if(workWeights(i1 + ioi) < workWeights(i2 + ioo)) then
        workWeights(i3) = workWeights(i1 + ioi)
        if(i3 /= 1) then
          wtlon(i3 - 1) = workWeights(i3) - workWeights(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if(.not.flagInput(2)) then
            mplon(i3 - 1, 1) = xDim / 2 + i1 - ici
            if(i1 - ici > xDim / 2) mplon(i3 - 1, 1) = i1 - ici - xDim / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if (.not.flagOutput(2)) then
            mplon(i3 - 1, 2) = xMax / 2 + i2 - ico
            if(i2 - ico > xMax / 2) mplon(i3 - 1, 2) = i2 - ico - xMax / 2
          end if
        end if
        i1 = i1 + 1
        i3 = i3 + 1
      else
        workWeights(i3) = workWeights(i2 + ioo)
        if (i3 /= 1)then
          wtlon(i3 - 1) = workWeights(i3) - workWeights(i3 - 1)
          mplon(i3 - 1, 1) = i1 - ici
          if(.not.flagInput(2)) then
            mplon(i3 - 1, 1) = xDim / 2 + i1 - ici
            if(i1 - ici > xDim / 2) mplon(i3 - 1, 1) = i1 - ici - xDim / 2
          end if
          mplon(i3 - 1, 2) = i2 - ico
          if(.not.flagOutput(2)) then
            mplon(i3 - 1, 2) = xMax / 2 + i2 - ico
            if(i2 - ico > xMax / 2) mplon(i3 - 1, 2) = i2 - ico - xMax / 2
          end if
        end if
        i2 = i2 + 1
        i3 = i3 + 1
      end if
      if(.not.(i1 <= xDim .and. i2 <= xMax)) exit
    end do
    if(i1 > xDim) i1 = 1
    if(i2 > xMax) i2 = 1
    do
      if(i2 /= 1) then
        workWeights(i3) = workWeights(i2 + ioo)
        wtlon(i3 - 1) = workWeights(i3) - workWeights(i3 - 1)
        mplon(i3 - 1, 1) = 1
        if(.not.(flagInput(4) .or. flagInput(5))) mplon(i3 - 1, 1) = xDim
        if(.not.flagInput(2)) then
          mplon(i3 - 1, 1) = xDim / 2 + 1
          if(.not.(flagInput(4) .or. flagInput(5))) mplon(i3 - 1, 1) = xDim / 2
        end if
        mplon(i3 - 1, 2) = i2 - ico
        if(.not.flagOutput(2)) then
          mplon(i3 - 1, 2) = xMax / 2 + i2 - ico
          if(i2 - ico > xMax / 2) mplon(i3 - 1, 2) = i2 - ico - xMax / 2
        end if
        i2 = i2 + 1
        if(i2 > xMax) i2 = 1
        i3 = i3 + 1
      end if
      if(i1 /= 1) then
        workWeights(i3) = workWeights(i1 + ioi)
        wtlon(i3 - 1) = workWeights(i3) - workWeights(i3 - 1)
        mplon(i3 - 1, 1) = i1 - ici
        if(.not.flagInput(2)) then
          mplon(i3 - 1, 1) = xDim / 2 + i1 - ici
          if(i1 - ici > xDim / 2) mplon(i3 - 1, 1) = i1 - ici - xDim / 2
        end if
        mplon(i3 - 1, 2) = 1
        if(.not.(flagOutput(4) .or. flagOutput(5))) mplon(i3 - 1, 2) = xMax
        if(.not.flagOutput(2)) then
          mplon(i3 - 1, 2) = xMax / 2 + 1
          if(.not.(flagOutput(4) .or. flagOutput(5))) mplon(i3 - 1, 2) = xMax / 2
        end if
        i1 = i1 + 1
        if(i1 > xDim) i1 = 1
        i3 = i3 + 1
      end if
      if(.not.(i1 /= 1 .or. i2 /= 1)) exit
    end do
    wtlon(i3 - 1) = 2.0_p_r8 * p_dpi + workWeights(1) - workWeights(i3 - 1)
    mplon(i3 - 1, 1) = 1
    if(.not.(flagInput(4) .or. flagInput(5))) mplon(i3 - 1, 1) = xDim
    if(.not.flagInput(2)) then
      mplon(i3 - 1, 1) = xDim / 2 + 1
      if(.not.(flagInput(4) .or. flagInput(5))) mplon(i3 - 1, 1) = xDim / 2
    end if
    mplon(i3 - 1, 2) = 1
    if(.not.(flagOutput(4) .or. flagOutput(5))) mplon(i3 - 1, 2) = xMax
    if(.not.flagOutput(2)) then
      mplon(i3 - 1, 2) = xMax / 2 + 1
      if(.not.(flagOutput(4) .or. flagOutput(5))) mplon(i3 - 1, 2) = xMax / 2
    end if
    lond = i3 - 1
  end subroutine getAreaInterpolationWeights
  

  subroutine gaussianLatitudes(latH, coLatitude)
    !# Gaussian Latitudes
    !# ---
    !# @info
    !# **Brief:** Gaussian Latitudes </br>
    !# **Authors**:</br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 <br>
    !# @endinfo
    
    integer, intent(in) :: latH              
    real(kind = p_r8), intent(out) :: coLatitude(latH)  

    real(kind = p_r8), parameter :: p_eps = epsilon(1.0_p_r8) * 100.0_p_r8 
    integer :: lats     
    integer :: j        
    !# loop iterator
    real(kind = p_r8) :: dGcolIn  
    real(kind = p_r8) :: gCol     
    real(kind = p_r8) :: dGcol    
    real(kind = p_r8) :: p1       
    real(kind = p_r8) :: p2       

    lats = 2 * latH
    dGcolIn = atan(1.0_p_r8) / real(lats, p_r8)
    gCol = 0.0_p_r8
    do j = 1, latH
      dGcol = dGcolIn
      do
        call legendrePolynomial(lats, gCol, p2)
        do
          p1 = p2
          gCol = gCol + dGcol
          call LegendrePolynomial(lats, gCol, p2)
          if(sign(1.0_p_r8, p1) /= sign(1.0_p_r8, p2)) exit
        end do
        if(dGcol <= p_eps) exit
        gCol = gCol - dGcol
        dGcol = dGcol * 0.25_p_r8
      end do
      coLatitude(j) = gCol
    end do

  end subroutine GaussianLatitudes

  subroutine legendrePolynomial(n, coLatitude, pln)
    !# Legendre Polynomial
    !# ---
    !# @info
    !# **Brief:** Legendre Polynomial </br>
    !# **Authors**:</br>
    !# &bull; Daniel Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo

    ! parameters
    integer, intent(in) :: n        
    real(kind = p_r8), intent(in) :: coLatitude 
    real(kind = p_r8), intent(out) :: pln        

    ! auxiliary variables
    integer :: i 
    !# loop iterator 
    real(kind = p_r8) :: x  
    real(kind = p_r8) :: y1 
    real(kind = p_r8) :: y2 
    real(kind = p_r8) :: y3 
    real(kind = p_r8) :: g  

    x = cos(coLatitude)
    y1 = 1.0_p_r8
    y2 = x
    do i = 2, n
      g = x * y2
      y3 = g - y1 + g - (g - y1) / real(i, p_r8)
      y1 = y2
      y2 = y3
    end do
    pln = y3
  end subroutine legendrePolynomial


end module Mod_AreaInterpolation
