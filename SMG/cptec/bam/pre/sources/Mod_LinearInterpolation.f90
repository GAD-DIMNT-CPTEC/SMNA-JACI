!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_LinearInterpolation </br></br>
!#
!# **Brief**: Linear Interpolator. Logics Assumes that Input and Output Data First
!# Point is Near North Pole and Greenwhich. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti   - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita            - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa  - version: 2.0.0 </li>
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

module Mod_LinearInterpolation

  implicit none
  private

  public :: latOut, initLinearInterpolation, doLinearInterpolation

  include 'precision.h'

  !parameters
  real(kind = p_r8), parameter :: p_Lat0 = 90.0_p_r8 
  !# Start at North Pole

  !input variables
  integer :: xDim          
  !# Number of longitude points for the input grid
  integer :: yDim          
  !# Number of latitude points for the input grid
  integer :: xMax          
  !# Number of longitude points for the output grid
  integer :: yMax          
  !# number of latitude points for the output grid

  integer, dimension(:), allocatable :: lowerLon 
  !# Lower Input Longitude Index
  integer, dimension(:), allocatable :: upperLon 
  !# Upper Input Longitude Index
  integer, dimension(:), allocatable :: lowerLat 
  !# Lower Input Latitude  Index
  integer, dimension(:), allocatable :: upperLat 
  !# Upper Input Latitude  Index
  real(kind = p_r8), dimension(:), allocatable :: lonIn  
  real(kind = p_r8), dimension(:), allocatable :: latIn  
  real(kind = p_r8), dimension(:), allocatable :: lonOut 
  real(kind = p_r8), dimension(:), allocatable :: latOut 

  real(kind = p_r8), dimension(:, :), allocatable :: leftLowerWgt 
  !# for Left-Lower  Corner of Box
  real(kind = p_r8), dimension(:, :), allocatable :: leftUpperWgt 
  !# for Left-Upper  Corner of Box
  real(kind = p_r8), dimension(:, :), allocatable :: rightLowerWgt 
  !# for Right-Lower Corner of Box
  real(kind = p_r8), dimension(:, :), allocatable :: rightUpperWgt 
  !# for Right-Upper Corner of Box


contains


  subroutine initLinearInterpolation(xDimIn, yDimIn, xMaxIn, yMaxIn, lat0, lon0)
    !# Initializes Linear Interpolation
    !# ---
    !# @info
    !# **Brief:** Allocates data and initializes the linear interpolation. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: xDimIn
    integer, intent(in) :: yDimIn
    integer, intent(in) :: xMaxIn
    integer, intent(in) :: yMaxIn
    real(kind = p_r8), intent(in) :: lat0
    real(kind = p_r8), intent(in) :: lon0

    xDim = xDimIn
    yDim = yDimIn
    xMax = xMaxIn
    yMax = yMaxIn

    call deallocateData()
    call allocateData()
    call getLongitudes(xDim, lon0, lonIn)
    call getLongitudes(xMax, 0.0_p_r8, lonOut)
    call getRegularLatitudes(yDim, lat0, latIn)
    call getGaussianLatitudes(yMax, latOut)
    call horizontalInterpolationWeights()

  end subroutine initLinearInterpolation


  subroutine allocateData()
    !# Allocates data
    !# ---
    !# @info
    !# **Brief:** Allocates data. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    allocate(lonIn(xDim))
    allocate(latIn(yDim))
    allocate(lonOut(xMax))
    allocate(latOut(yMax))
    allocate(lowerLon(xMax))
    allocate(lowerLat(yMax))
    allocate(upperLon(xMax))
    allocate(upperLat(yMax))
    allocate(leftLowerWgt(xMax, yMax))
    allocate(leftUpperWgt(xMax, yMax))
    allocate(rightLowerWgt(xMax, yMax))
    allocate(rightUpperWgt(xMax, yMax))
  end subroutine allocateData


  subroutine deallocateData()
    !# Deallocates data
    !# ---
    !# @info
    !# **Brief:** Deallocates data. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 </br>
    !# @endinfo
    implicit none
    if(allocated(lonIn)) then
      deallocate(lonIn)
      deallocate(latIn)
      deallocate(lonOut)
      deallocate(latOut)
      deallocate(lowerLon)
      deallocate(lowerLat)
      deallocate(upperLon)
      deallocate(upperLat)
      deallocate(leftLowerWgt)
      deallocate(leftUpperWgt)
      deallocate(rightLowerWgt)
      deallocate(rightUpperWgt)
    end if
  end subroutine deallocateData


  subroutine doLinearInterpolation(varIn, varOut)
    !# Does Linear Horizontal Interpolation
    !# ---
    !# @info
    !# **Brief:** Does Linear Horizontal Interpolation. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in), dimension(xDim, yDim) :: varIn  
    real(kind = p_r8), intent(out), dimension(xMax, yMax) :: varOut 

    call horizontalInterpolation(varIn, varOut)
  end subroutine doLinearInterpolation


  subroutine getLongitudes(xMax, lon0, lon)
    !# Gets Longitudes
    !# ---
    !# @info
    !# **Brief:** Gets Longitudes given the number of longitudes and the first
    !# longitude. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: xMax 
    !# Number of Longitudes
    real(kind = p_r8), intent(in) :: lon0 
    !# First Longitude In Degree
    real(kind = p_r8), intent(out), dimension(xMax) :: lon 
    !# Longitudes In Degree

    integer :: i
    real(kind = p_r8) :: dx

    dx = 360.0_p_r8 / real(xMax, p_r8)
    do i = 1, xMax
      lon(i) = lon0 + real(i - 1, p_r8) * dx
    end do

  end subroutine getLongitudes


  subroutine getRegularLatitudes(yMax, lat0, lat)
    !# Gets Regular Latitudes
    !# ---
    !# @info
    !# **Brief:** Gets Regular Latitudes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: yMax 
    !# Number of Regular Latitudes
    real(kind = p_r8), intent(in) :: lat0 
    !# First Latitude In Degree
    real(kind = p_r8), intent(out), dimension(yMax) :: lat 
    !# Regular Latitudes In Degree

    integer :: j
    real(kind = p_r8) :: dy

    dy = 2.0_p_r8 * lat0 / real(yMax - 1, p_r8)
    do j = 1, yMax
      lat(j) = lat0 - real(j - 1, p_r8) * dy
    end do

  end subroutine getRegularLatitudes


  subroutine getGaussianLatitudes(yMax, lat)
    !# Gets Gaussian Latitudes
    !# ---
    !# @info
    !# **Brief:** Gets Gaussian Latitudes given Number of Gaussian Latitudes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: yMax 
    !# Number of Gaussian Latitudes
    real(kind = p_r8), intent(out), dimension(yMax) :: lat 
    !# Gaussian Latitudes In Degree
    
    integer :: j
    real(kind = p_r8) :: eps
    real(kind = p_r8) :: rd
    real(kind = p_r8) :: dCoLatRadz
    real(kind = p_r8) :: coLatRad
    real(kind = p_r8) :: dCoLatRad
    real(kind = p_r8) :: p2
    real(kind = p_r8) :: p1

    eps = 1.0e-12_p_r8
    rd = 45.0_p_r8 / atan(1.0_p_r8)
    dCoLatRadz = ((180.0_p_r8 / real(yMax, p_r8)) / rd) / 10.0_p_r8
    coLatRad = 0.0_p_r8
    do j = 1, yMax / 2
      dCoLatRad = dCoLatRadz
      do while(dCoLatRad > eps)
        call legendrePolynomial(yMax, coLatRad, p2)
        do
          p1 = p2
          coLatRad = coLatRad + dCoLatRad
          call legendrePolynomial(yMax, coLatRad, p2)
          if(sign(1.0_p_r8, p1) /= sign(1.0_p_r8, p2)) exit
        end do
        coLatRad = coLatRad - dCoLatRad
        dCoLatRad = dCoLatRad * 0.25_p_r8
      end do
      lat(j) = 90.0_p_r8 - coLatRad * rd
      lat(yMax - j + 1) = -Lat(j)
      coLatRad = coLatRad + dCoLatRadz
    end do
  end subroutine getGaussianLatitudes


  subroutine legendrePolynomial(n, coLatRad, legPol)
    !# Legendre Polynomial
    !# ---
    !# @info
    !# **Brief:** Legendre Polynomial. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: n 
    !# Order of the Ordinary Legendre Function
    real(kind = p_r8), intent(in) :: coLatRad 
    !# Colatitude (In Radians)
    real(kind = p_r8), intent(out) :: legPol 
    !# Value of The Ordinary Legendre Function

    integer :: i
    real(kind = p_r8) :: x
    real(kind = p_r8) :: y1
    real(kind = p_r8) :: y2
    real(kind = p_r8) :: g
    real(kind = p_r8) :: y3

    x = cos(coLatRad)
    y1 = 1.0_p_r8
    y2 = x
    do i = 2, n
      g = x * y2
      y3 = g - y1 + g - (g - y1) / real(i, p_r8)
      y1 = y2
      y2 = y3
    end do
    legPol = y3

  end subroutine legendrePolynomial


  subroutine horizontalInterpolationWeights()
    !# Horizontal Interpolation Weights
    !# ---
    !# @info
    !# **Brief:** Calculates the Horizontal Interpolation Weights. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer :: i
    integer :: j
    integer :: iLon
    integer :: jLat

    real(kind = p_r8) :: dLon
    real(kind = p_r8) :: lowerLonIn
    real(kind = p_r8) :: upperLonIn
    real(kind = p_r8) :: dbx
    real(kind = p_r8) :: latOutA
    real(kind = p_r8) :: dLat
    real(kind = p_r8) :: lowerLatIn
    real(kind = p_r8) :: upperLatIn
    real(kind = p_r8) :: dby

    integer, dimension(1) :: iLonm
    integer, dimension(1) :: jLatm

    real(kind = p_r8), dimension(xMax) :: dx
    real(kind = p_r8), dimension(xMax) :: ddx
    real(kind = p_r8), dimension(yMax) :: dy
    real(kind = p_r8), dimension(yMax) :: ddy

    do i = 1, xMax
      iLonm = minloc(abs(lonIn(1:xDim) - lonOut(i)))
      iLon = iLonm(1)
      dLon = lonOut(i) - lonIn(iLon)
      if(dLon < 0.0_p_r8) then
        if(iLon /= 1) then
          lowerLon(i) = iLon - 1
          lowerLonIn = lonIn(iLon - 1)
        else
          lowerLon(i) = xDim
          lowerLonIn = lonIn(xDim) - 360.0_p_r8
        end if
        upperLon(i) = iLon
        upperLonIn = lonIn(iLon)
      else if(dLon == 0.0_p_r8) then
        lowerLon(i) = iLon
        lowerLonIn = lonIn(iLon)
        upperLon(i) = iLon
        upperLonIn = lonIn(iLon)
      else
        lowerLon(i) = iLon
        lowerLonIn = lonIn(iLon)
        if(iLon /= xDim) then
          upperLon(i) = iLon + 1
          upperLonIn = lonIn(iLon + 1)
        else
          upperLon(i) = 1
          upperLonIn = lonIn(1) + 360.0_p_r8
        end if
      end if
      dx(i) = lonOut(i) - lowerLonIn
      ddx(i) = upperLonIn - lowerLonIn
    end do

    do j = 1, yMax
      jLatm = minloc(abs(latIn(1:yDim) - latOut(j)))
      jLat = jLatm(1)
      latOutA = latOut(j)
      if(latOutA > latIn(1))    latOutA = latIn(1)
      if(latOutA < latIn(yDim)) latOutA = latIn(yDim)
      dLat = latOutA - latIn(jLat)
      if(dLat > 0.0_p_r8) then
        if(jLat /= 1) then
          lowerLat(j) = jLat - 1
          lowerLatIn = latIn(jLat - 1)
        else
          lowerLat(j) = 1
          lowerLatIn = latIn(1)
        end if
        upperLat(j) = jLat
        upperLatIn = latIn(jLat)
      else if(dLat == 0.0_p_r8) then
        lowerLat(j) = jLat
        lowerLatIn = latIn(jLat)
        upperLat(j) = jLat
        upperLatIn = latIn(jLat)
      else
        lowerLat(j) = jLat
        lowerLatIn = latIn(jLat)
        if(jLat /= yDim) then
          upperLat(j) = jLat + 1
          upperLatIn = latIn(jLat + 1)
        else
          upperLat(j) = yDim
          upperLatIn = latIn(yDim)
        end if
      end if
      dy(j) = latOutA - lowerLatIn
      ddy(j) = upperLatIn - lowerLatIn
    end do

    do j = 1, yMax
      do i = 1, xMax
        if(ddx(i) == 0.0_p_r8 .and. ddy(j) == 0.0_p_r8) then
          leftLowerWgt(i, j) = 1.0_p_r8
          leftUpperWgt(i, j) = 0.0_p_r8
          rightLowerWgt(i, j) = 0.0_p_r8
          rightUpperWgt(i, j) = 0.0_p_r8
        else if(ddx(i) == 0.0_p_r8) then
          leftUpperWgt(i, j) = dy(j) / ddy(j)
          leftLowerWgt(i, j) = 1.0_p_r8 - leftUpperWgt(i, j)
          rightLowerWgt(i, j) = 0.0_p_r8
          rightUpperWgt(i, j) = 0.0_p_r8
        else if(ddy(j) == 0.0_p_r8) then
          rightLowerWgt(i, j) = dx(i) / ddx(i)
          leftLowerWgt(i, j) = 1.0_p_r8 - rightLowerWgt(i, j)
          leftUpperWgt(i, j) = 0.0_p_r8
          rightUpperWgt(i, j) = 0.0_p_r8
        else
          dbx = dx(i) / ddx(i)
          dby = dy(j) / ddy(j)
          rightUpperWgt(i, j) = dbx * dby
          leftLowerWgt(i, j) = 1.0_p_r8 - dbx - dby + rightUpperWgt(i, j)
          leftUpperWgt(i, j) = dby - rightUpperWgt(i, j)
          rightLowerWgt(i, j) = dbx - rightUpperWgt(i, j)
        end if
      end do
    end do

  end subroutine horizontalInterpolationWeights


  subroutine horizontalInterpolation(varIn, varOut)
    !# Horizontal Interpolation
    !# ---
    !# @info
    !# **Brief:** Calculates the Horizontal Interpolation Weights. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    
    real(kind = p_r8), intent(in), dimension(xDim, yDim) :: varIn  
    real(kind = p_r8), intent(out), dimension(xMax, yMax) :: varOut
    !# subroutine return

    integer :: j    ! j - Latitude Index
    integer :: i    ! i - Longitude Index
    integer :: il
    integer :: iu
    integer :: jl
    integer :: ju
    
    !# Input Grid Box That Contains The Output Value (i,j):\n
    do j = 1, yMax
      do i = 1, xMax
        il = lowerLon(i)     ! lowerLon(i,j) - Lower Input Longitude Index\n
        jl = lowerLat(j)     ! lowerLat(i,j) - Lower Input Latitude  Index\n
        iu = upperLon(i)     ! upperLon(i,j) - Upper Input Longitude Index
        ju = upperLat(j)     ! upperLat(i,j) - Upper Input Latitude  Index
        if(varIn(il, jl) == p_undef .and. varIn(il, ju) == p_undef .and. &
          varIn(iu, jl) == p_undef .and. varIn(iu, ju) == p_undef) then
          varOut(i, j) = p_undef
        else
          varOut(i, j) = 0.0_p_r8
          ! Pre-Calculated Weights for Linear Horizontal Interpolation:
          ! leftLowerWgt  - for Left-Lower  Corner of Box
          ! leftUpperWgt  - for Left-Upper  Corner of Box
          ! rightLowerWgt - for Right-Lower Corner of Box
          ! rightUpperWgt - for Right-Upper Corner of Box
          if(varIn(il, jl) /= p_undef) varOut(i, j) = varOut(i, j) + leftLowerWgt(i, j) * &
            varIn(il, jl)
          if(varIn(il, ju) /= p_undef) varOut(i, j) = varOut(i, j) + leftUpperWgt(i, j) * &
            varIn(il, ju)
          if(varIn(iu, jl) /= p_undef) varOut(i, j) = varOut(i, j) + rightLowerWgt(i, j) * &
            varIn(iu, jl)
          if(varIn(iu, ju) /= p_undef) varOut(i, j) = varOut(i, j) + rightUpperWgt(i, j) * &
            varIn(iu, ju)
        end if
      end do
    end do

  end subroutine horizontalInterpolation

end module Mod_LinearInterpolation
