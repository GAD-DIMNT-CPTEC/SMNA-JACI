!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_LegendreTransform </br></br>
!#
!# **Brief**: LEGENDRE TRANSFORMS </br></br>
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
!# **Version**: 2.0.1 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!# <table>
!#  <tr><td style="width:85px"><li>13-11-2004 </li></td> <td style="width:130px">- Jose P. Bonatti </td>  <td>- version: 1.0.0 </td></tr>
!#  <tr><td><li>01-08-2007 </li></td>                    <td>- Tomita </td>                               <td>- version: 1.1.1 </td></tr>
!#  <tr><td><li>01-04-2018 </li></td>                    <td>- Daniel M. Lamosa</br>- Barbara Yamada </td><td>- version: 2.0.0 </td></tr>
!#  <tr><td><li>08-04-2019 </li></td>                    <td>- Eduardo Khamis </td>                       <td>- version: 2.0.1 </td></tr>
!# </table>
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

module Mod_LegendreTransform

  implicit none
  private

  public :: createLegTrans, createSpectralRep, createGaussRep, &
    transs, spec2Four, four2Spec, splitTrans, gLats, la0, destroyLegendreObjects, &
     emRad, emRad1, emRad12, emRad2

  !parameters
  integer, parameter :: p_r8 = selected_real_kind(15)        
  !# kind for 64bits
  
  !input variables
  real(kind = p_r8), parameter :: emRad   = 6.37E6_p_r8      
  !# Earth Mean Radius (m)
  real(kind = p_r8), parameter :: emRad1  = 1.0_p_r8/emRad   
  !# 1/emRad (1/m)
  real(kind = p_r8), parameter :: emRad12 = emRad1**2        
  !# emRad1**2 (1/m2)
  real(kind = p_r8), parameter :: emRad2  = emRad**2         
  !# emRad**2 (m2)

  integer, parameter :: p_nferr = 0

  integer, allocatable, dimension(:) :: lenDiag         
  !
  integer, allocatable, dimension(:) :: lenDiagExt      
  !
  integer, allocatable, dimension(:) :: lastPrevDiag    
  !
  integer, allocatable, dimension(:) :: lastPrevDiagExt 
  !
  integer, allocatable, dimension(:, :) :: la0             
  !
  integer, allocatable, dimension(:, :) :: la1             
  !

  real(kind = p_r8), allocatable, dimension(:) :: eps             
  !
  real(kind = p_r8), allocatable, dimension(:) :: colRad          
  !
  real(kind = p_r8), allocatable, dimension(:) :: rCs2            
  !
  real(kind = p_r8), allocatable, dimension(:) :: wgt             
  !
  real(kind = p_r8), allocatable, dimension(:) :: gLats           
  !
  real(kind = p_r8), allocatable, dimension(:, :) :: legS2F          
  !
  real(kind = p_r8), allocatable, dimension(:, :) :: legExtS2F       
  !
  real(kind = p_r8), allocatable, dimension(:, :) :: legDerS2F       
  !
  real(kind = p_r8), allocatable, dimension(:, :) :: legF2S          
  !
  real(kind = p_r8), allocatable, dimension(:, :) :: legDerNS        
  !
  real(kind = p_r8), allocatable, dimension(:, :) :: legDerEW        
  !

  logical :: created = .false.

  interface splitTrans
    module procedure splitTrans2D, splitTrans3D
  end interface

  interface spec2Four
    module procedure spec2Four1D, spec2Four2D
  end interface

  interface four2Spec
    module procedure four2Spec1D, four2Spec2D
  end interface


contains


  subroutine createSpectralRep(mEnd1In, mEnd2In, mnwv1In)
    !# Creates Spectral Representation
    !# ---
    !# @info
    !# **Brief:** Creates Spectral Representation. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: mEnd2In
    integer, intent(in) :: mnwv1In

    integer :: l          
    integer :: mm         
    integer :: nn         
    real(kind = p_r8) :: am 
    real(kind = p_r8) :: an 

    allocate(la0(mEnd1In, mEnd1In))
    allocate(la1(mEnd1In, mEnd2In))
    allocate(eps(mnwv1In))

    l = 0
    do nn = 1, mEnd1In
      do mm = 1, mEnd2In - nn
        l = l + 1
        la0(mm, nn) = l
      end do
    end do
    l = 0
    do mm = 1, mEnd1In
      l = l + 1
      la1(mm, 1) = l
    end do
    do nn = 2, mEnd2In
      do mm = 1, mEnd1In + 2 - nn
        l = l + 1
        la1(mm, nn) = l
      end do
    end do

    do l = 1, mEnd1In
      eps(l) = 0.0_p_r8
    end do
    l = mEnd1In
    do nn = 2, mEnd2In
      do mm = 1, mEnd1In + 2 - nn
        l = l + 1
        am = mm - 1
        an = mm + nn - 2
        eps(l) = sqrt((an * an - am * am) / (4.0_p_r8 * an * an - 1.0_p_r8))
      end do
    end do

  end subroutine createSpectralRep


  subroutine createGaussRep(yMaxIn, yMaxHfIn)
    !# Creates Gaussian Representation
    !# ---
    !# @info
    !# **Brief:** Creates Gaussian Representation. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    integer, intent(in) :: yMaxIn
    integer, intent(in) :: yMaxHfIn

    allocate(colRad(yMaxHfIn))
    allocate(rCs2(yMaxHfIn))
    allocate(wgt(yMaxHfIn))
    allocate(gLats(yMaxIn))

    call gaussianLatitudes(yMaxIn, yMaxHfIn)

  end subroutine createGaussRep


  subroutine gaussianLatitudes(yMaxIn, yMaxHfIn)
    !# Calculates Gaussian Latitudes and Gaussian Weights
    !# ---
    !# @info
    !# **Brief:** Calculates Gaussian Latitudes and Gaussian Weights for Use in
    !# Grid-Spectral and Spectral-Grid Transforms. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    
    !    note: pgi failed to compile epsil, scal and dgColIn as parameter
    !    real(kind=p_r8), parameter :: epsil   = epsilon(1.0_p_r8) * 100.0_p_r8    
    !    real(kind=p_r8), parameter :: scal    = 2.0_p_r8 / (real(yMax, p_r8)&
    !                                          * real(yMax, p_r8))                 ! scale
    !    real(kind=p_r8), parameter :: dgColIn = atan(1.0_p_r8) / real(yMax, p_r8) 
    integer, intent(in) :: yMaxIn
    integer, intent(in) :: yMaxHfIn

    real(kind = p_r8) :: epsil, scal, dgColIn

    integer :: j       
    !# Loop iterator
    real(kind = p_r8) :: gCol    
    real(kind = p_r8) :: dgCol 
    real(kind = p_r8) :: p2      
    real(kind = p_r8) :: p1      
    real(kind = p_r8) :: rad     

    epsil = epsilon(1.0_p_r8) * 100.0_p_r8
    rad = 45.0_p_r8 / atan(1.0_p_r8)
    scal = 2.0_p_r8 / (real(yMaxIn, p_r8) * real(yMaxIn, p_r8))
    dgColIn = atan(1.0_p_r8) / real(yMaxIn, p_r8)

    gCol = 0.0_p_r8
    do j = 1, yMaxHfIn
      dgCol = dgColIn
      do
        call legendrePolynomial(yMaxIn, gCol, p2)
        do
          p1 = p2
          gCol = gCol + dgCol
          call legendrePolynomial(yMaxIn, gCol, p2)
          if(sign(1.0_p_r8, p1) /= sign(1.0_p_r8, p2)) exit
        end do
        if(dgCol <= epsil) exit
        gCol = gCol - dgCol
        dgCol = dgCol * 0.25_p_r8
      end do
      colRad(j) = gCol
      gLats(j) = 90.0_p_r8 - rad * gCol
      gLats(yMaxIn - j + 1) = -gLats(j)
      call legendrePolynomial(yMaxIn - 1, gCol, p1)
      wgt(j) = scal * (1.0_p_r8 - cos(gCol) * cos(gCol)) / (p1 * p1)
      rCs2(j) = 1.0_p_r8 / (sin(gCol) * sin(gCol))
    end do

  end subroutine gaussianLatitudes


  subroutine legendrePolynomial(n, colatitude, pln)
    !# Calculates the Value of the Ordinary Legendre Function
    !# ---
    !# @info
    !# **Brief:** Calculates the Value of the Ordinary Legendre Function of Given
    !# Order at a Specified Colatitude. Used to Determine Gaussian Latitudes. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: n          
    real(kind = p_r8), intent(in) :: colatitude 
    real(kind = p_r8), intent(out) :: pln        

    integer :: i  ! Loop iterator
    real(kind = p_r8) :: x  
    real(kind = p_r8) :: y1 
    real(kind = p_r8) :: y2 
    real(kind = p_r8) :: y3 
    real(kind = p_r8) :: g  

    x = cos(colatitude)
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


  subroutine createLegTrans(mnwv0In, mnwv1In, mnwv2In, mnwv3In, mEnd1In, mEnd2In, yMaxHfIn)
    !# Creates Legendre Transforms
    !# ---
    !# @info
    !# **Brief:** Creates Legendre Transforms. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    integer, intent(in) :: mnwv0In
    integer, intent(in) :: mnwv1In
    integer, intent(in) :: mnwv2In
    integer, intent(in) :: mnwv3In
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: mEnd2In
    integer, intent(in) :: yMaxHfIn

    character(len = *), parameter :: h = "**(createLegTrans)**" 

    integer :: diag 

    if(created) then
      write(unit = p_nferr, fmt = '(2a)') h, ' already created'
      stop
    else
      created = .true.
    end if

    ! Associated Legendre Functions
    allocate(legS2F   (mnwv2In, yMaxHfIn))
    allocate(legDerS2F(mnwv2In, yMaxHfIn))
    allocate(legExtS2F(mnwv3In, yMaxHfIn))
    allocate(legF2S   (mnwv2In, yMaxHfIn))
    allocate(legDerNS (mnwv2In, yMaxHfIn))
    allocate(legDerEW (mnwv2In, yMaxHfIn))
    call legPols(mnwv0In, mnwv1In, mnwv2In, mEnd1In, mEnd2In, yMaxHfIn)

    ! diagonal length
    allocate(lenDiag(mEnd1In))
    allocate(lenDiagExt(mEnd2In))
    do diag = 1, mEnd1In
      lenDiag(diag) = 2 * (mEnd1In + 1 - diag)
    end do
    lenDiagExt(1) = 2 * mEnd1In
    do diag = 2, mEnd2In
      lenDiagExt(diag) = 2 * (mEnd1In + 2 - diag)
    end do

    ! last element previous diagonal
    allocate(lastPrevDiag(mEnd1In))
    allocate(lastPrevDiagExt(mEnd2In))
    do diag = 1, mEnd1In
      lastPrevDiag(diag) = (diag - 1) * (2 * mEnd1In + 2 - diag)
    end do
    lastPrevDiagExt(1) = 0
    do diag = 2, mEnd2In
      lastPrevDiagExt(diag) = (diag - 1) * (2 * mEnd1In + 4 - diag) - 2
    end do

  end subroutine createLegTrans


  subroutine legPols(mnwv0In, mnwv1In, mnwv2In, mEnd1In, mEnd2In, yMaxHfIn)
    !# Legendre Polynomials
    !# ---
    !# @info
    !# **Brief:** Legendre Polynomials. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    integer, intent(in) :: mnwv0In
    integer, intent(in) :: mnwv1In
    integer, intent(in) :: mnwv2In
    integer, intent(in) :: mEnd1In
    integer, intent(in) :: mEnd2In
    integer, intent(in) :: yMaxHfIn

    integer :: j  
    integer :: l  
    integer :: nn 
    integer :: mm 
    integer :: mn 
    integer :: lx 

    real(kind = p_r8) :: pln(mnwv1In)
    real(kind = p_r8) :: dpln(mnwv0In)
    real(kind = p_r8) :: der(mnwv0In)
    real(kind = p_r8) :: plnwcs(mnwv0In)

    do j = 1, yMaxHfIn
      call pln2(mnwv1In, mEnd1In, mEnd2In, yMaxHfIn, pln, colrad, j, eps, la1)
      l = 0
      do nn = 1, mEnd1In
        do mm = 1, mEnd2In - nn
          l = l + 1
          lx = la1(mm, nn)
          legS2F(2 * l - 1, j) = pln(lx)
          legS2F(2 * l, j) = pln(lx)
        end do
      end do
      do mn = 1, mnwv2In
        legF2S(mn, j) = legS2F(mn, j) * wgt(j)
      end do
      do mn = 1, mnwv1In
        legExtS2F(2 * mn - 1, j) = pln(mn)
        legExtS2F(2 * mn, j) = pln(mn)
      end do
      call plnder(mnwv0In, mnwv1In, mEnd1In, mEnd2In, pln, dpln, der, plnwcs, rcs2(j), wgt(j), eps, la1)
      do mn = 1, mnwv0In
        legDerS2F(2 * mn - 1, j) = dpln(mn)
        legDerS2F(2 * mn, j) = dpln(mn)
        legDerNS(2 * mn - 1, j) = der(mn)
        legDerNS(2 * mn, j) = der(mn)
        legDerEW(2 * mn - 1, j) = plnwcs(mn)
        legDerEW(2 * mn, j) = plnwcs(mn)
      end do
    end do

  end subroutine legPols


  subroutine pln2 (mnwv1, mEnd1, mEnd2, yMaxHf, sln, colrad, lat, eps, la1)
    !# Calculates the associated legendre functions
    !# ---
    !# @info
    !# **Brief:** Calculates the associated legendre functions at one specified
    !# latitude. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: mnwv1
    integer, intent(in) :: mEnd1
    integer, intent(in) :: mEnd2
    integer, intent(in) :: yMaxHf

    real(kind = p_r8), intent(out) :: sln(mnwv1)
    real(kind = p_r8), intent(in) :: colrad(yMaxHf)
    real(kind = p_r8), intent(in) :: eps(mnwv1)

    integer, intent(in) :: lat
    integer, intent(in) :: la1(mEnd1, mEnd2)

    integer :: mm, nn, lx, ly, lz

    real(kind = p_r8) :: colr, sinlat, coslat, prod

    logical, save :: first = .true.

    real(kind = p_r8), allocatable, dimension(:), save :: x, y
    real(kind = p_r8), save :: rthf

    if (first) then
      allocate (x(mEnd1))
      allocate (y(mEnd1))
      first = .false.
      do mm = 1, mEnd1
        x(mm) = sqrt(2.0_p_r8 * mm + 1.0_p_r8)
        y(mm) = sqrt(1.0_p_r8 + 0.5_p_r8 / real(mm, p_r8))
      end do
      rthf = sqrt(0.5_p_r8)
    endif
    colr = colrad(lat)
    sinlat = cos(colr)
    coslat = sin(colr)
    prod = 1.0_p_r8
    do mm = 1, mEnd1
      sln(mm) = rthf * prod
      !     line below should only be used where exponent range is limted
      !     if(prod < flim) prod=0.0_p_r8
      prod = prod * coslat * y(mm)
    end do

    do mm = 1, mEnd1
      sln(mm + mEnd1) = x(mm) * sinlat * sln(mm)
    end do
    do nn = 3, mEnd2
      do mm = 1, mEnd1 + 2 - nn
        lx = la1(mm, nn)
        ly = la1(mm, nn - 1)
        lz = la1(mm, nn - 2)
        sln(lx) = (sinlat * sln(ly) - eps(ly) * sln(lz)) / eps(lx)
      end do
    end do

  end subroutine pln2


  subroutine plnder (mnwv0, mnwv1, mEnd1, mEnd2, pln, dpln, der, plnwcs, rcs2l, wgtl, eps, la1)
    !# Calculates zonal and meridional pseudo-derivatives
    !# ---
    !# @info
    !# **Brief:** Calculates zonal and meridional pseudo-derivatives as well as
    !# laplacians of the associated legendre functions. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: mnwv0
    integer, intent(in) :: mnwv1
    integer, intent(in) :: mEnd1
    integer, intent(in) :: mEnd2

    real(kind = p_r8), intent(inout) :: pln(mnwv1)
    real(kind = p_r8), intent(out) :: dpln(mnwv0)
    real(kind = p_r8), intent(out) :: der(mnwv0)
    real(kind = p_r8), intent(out) :: plnwcs(mnwv0)
    real(kind = p_r8), intent(in) :: rcs2l
    real(kind = p_r8), intent(in) :: wgtl
    real(kind = p_r8), intent(in) :: eps(mnwv1)

    integer, intent(in) :: la1(mEnd1, mEnd2)

    integer :: n, l, nn, mm, mn, lm, l0, lp

    real(kind = p_r8) :: raa, wcsa
    real(kind = p_r8) :: x(mnwv1)

    logical, save :: first = .true.

    real(kind = p_r8), allocatable, save :: an(:)
    ! 
    !     compute pln derivatives
    ! 
    if (first) then
      allocate (an(mEnd2))
      do n = 1, mEnd2
        an(n) = real(n - 1, p_r8)
      end do
      first = .false.
    end if
    raa = wgtl * emRad12
    wcsa = rcs2l * wgtl * emRad1
    l = 0
    do mm = 1, mEnd1
      l = l + 1
      x(l) = an(mm)
    end do
    do nn = 2, mEnd2
      do mm = 1, mEnd1 + 2 - nn
        l = l + 1
        x(l) = an(mm + nn - 1)
      end do
    end do
    l = mEnd1
    do nn = 2, mEnd1
      do mm = 1, mEnd2 - nn
        l = l + 1
        lm = la1(mm, nn - 1)
        l0 = la1(mm, nn)
        lp = la1(mm, nn + 1)
        der(l) = x(lp) * eps(l0) * pln(lm) - x(l0) * eps(lp) * pln(lp)
      end do
    end do
    do mm = 1, mEnd1
      der(mm) = -x(mm) * eps(mm + mEnd1) * pln(mm + mEnd1)
    end do
    do mn = 1, mnwv0
      dpln(mn) = der(mn)
      der(mn) = wcsa * der(mn)
    end do
    l = 0
    do nn = 1, mEnd1
      do mm = 1, mEnd2 - nn
        l = l + 1
        l0 = la1(mm, nn)
        plnwcs(l) = an(mm) * pln(l0)
      end do
    end do
    do mn = 1, mnwv0
      plnwcs(mn) = wcsa * plnwcs(mn)
    end do
    do nn = 1, mEnd1
      do mm = 1, mEnd2 - nn
        l0 = la1(mm, nn)
        lp = la1(mm, nn + 1)
        pln(l0) = x(l0) * x(lp) * raa * pln(l0)
      end do
    end do

  end subroutine plnder


  subroutine transs (mnwv2, mEnd1, mEnd2, lDim, signal, a)
    !# Transposes scalar arrays of spectral coefficients
    !# ---
    !# @info
    !# **Brief:** After input, transposes scalar arrays of spectral coefficients
    !# by swapping the order of the subscripts representing the degree and order
    !# of the associated legendre functions. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none
    !  argument(dimensions)                description
    ! 
    !     lDim                 input: number of layers. 
    !     a(mnwv2,lDim)        input: spectral representation of a global field
    !                                 at "n" levels. 
    !                                 signal=+1 diagonalwise storage 
    !                                 signal=-1 coluMnwise   storage 
    !                         output: spectral representation of a global field
    !                                 at "n" levels. 
    !                                 signal=+1 coluMnwise   storage 
    !                                 signal=-1 diagonalwise storage 
    integer, intent(in) :: mnwv2
    integer, intent(in) :: mEnd1
    integer, intent(in) :: mEnd2
    integer, intent(in) :: lDim
    integer, intent(in) :: signal
    real(kind = p_r8), intent(inout) :: a(mnwv2, lDim)

    real(kind = p_r8) :: w(mnwv2)

    integer :: k
    integer :: l
    integer :: lx
    integer :: mn
    integer :: mm
    integer :: nlast
    integer :: nn

    if (signal == 1) then
      do k = 1, lDim
        l = 0
        do mm = 1, mEnd1
          nlast = mEnd2 - mm
          do nn = 1, nlast
            l = l + 1
            lx = la0(mm, nn)
            w(2 * l - 1) = a(2 * lx - 1, k)
            w(2 * l) = a(2 * lx, k)
          end do
        end do
        do mn = 1, mnwv2
          a(mn, k) = w(mn)
        end do
      end do
    else
      do k = 1, lDim
        l = 0
        do mm = 1, mEnd1
          nlast = mEnd2 - mm
          do nn = 1, nlast
            l = l + 1
            lx = la0(mm, nn)
            w(2 * lx - 1) = a(2 * l - 1, k)
            w(2 * lx) = a(2 * l, k)
          end do
        end do
        do mn = 1, mnwv2
          a(mn, k) = w(mn)
        end do
      end do
    end if

  end subroutine transs


  subroutine sumSpec (nMax, mMax, mnwv, xmx, yMax, yMaxHf, zMax, &
    spec, leg, four, dLength, lastPrev)
    !# Fourier representation from spectral representation
    !# ---
    !# @info
    !# **Brief:** Fourier representation from spectral representation. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: nMax           
    !# assoc leg func degrees (=trunc+1 or +2)
    integer, intent(in) :: mMax           
    !# assoc leg func waves (=trunc+1)
    integer, intent(in) :: mnwv           
    !# spectral coef, real+imag 
    !  integer, intent(in) :: xMax           !# longitudes * 2 (real+imag four coef)
    integer, intent(in) :: xmx            
    !# longitudes + 2 (real+imag four coef)
    integer, intent(in) :: yMax           
    !# latitudes (full sphere)
    integer, intent(in) :: yMaxHf         
    !# latitudes (hemisphere)
    integer, intent(in) :: zMax           
    !# verticals
    integer, intent(in) :: dLength(nMax)  
    !# diagonal length (real+imag)
    integer, intent(in) :: lastPrev(nMax) 
    !# last element previous diagonal (real+imag)

    real(kind = p_r8), intent(in) :: spec(mnwv, zMax)      
    !# spectral field
    real(kind = p_r8), intent(in) :: leg (mnwv, yMaxHf)    
    !# associated legendre function
!    real(kind = p_r8), intent(out) :: four(xMax, yMax, zMax) !# full fourier field
    real(kind = p_r8), intent(out) :: four(xmx, yMax, zMax) 
    !# full fourier field

    integer :: j, k, ele, diag

    real(kind = p_r8) :: oddDiag(2 * mMax, yMaxHf, zMax)
    real(kind = p_r8) :: evenDiag(2 * mMax, yMaxHf, zMax)

    ! initialize diagonals

    oddDiag = 0.0_p_r8
    evenDiag = 0.0_p_r8

    ! sum odd diagonals (n+m even)

    do k = 1, zMax
      do j = 1, yMaxHf
        do diag = 1, nMax, 2
          do ele = 1, dLength(diag)
            oddDiag(ele, j, k) = oddDiag(ele, j, k) + &
              leg(ele + lastPrev(diag), j) * spec(ele + lastPrev(diag), k)
          end do
        end do
      end do
    end do

    ! sum even diagonals (n+m odd)

    do k = 1, zMax
      do j = 1, yMaxHf
        do diag = 2, nMax, 2
          do ele = 1, dLength(diag)
            evenDiag(ele, j, k) = evenDiag(ele, j, k) + &
              leg(ele + lastPrev(diag), j) * spec(ele + lastPrev(diag), k)
          end do
        end do
      end do
    end do

    ! use even-odd properties

    !$cdir nodep
    do k = 1, zMax
      do j = 1, yMaxHf
        do ele = 1, 2 * mMax
          four(ele, j, k) = oddDiag(ele, j, k) + evenDiag(ele, j, k)
          four(ele, yMax + 1 - j, k) = oddDiag(ele, j, k) - evenDiag(ele, j, k)
        end do
      end do
    end do
!    four(2 * mMax + 1:xMax, :, :) = 0.0_p_r8
    four(2 * mMax + 1:xmx, :, :) = 0.0_p_r8

  end subroutine sumSpec


  subroutine sumFour (nMax, mMax, mnwv, xmx, yMax, yMaxHf, zMax, &
    spec, leg, four, dLength, lastPrev)
    !# Spectral representation from fourier representation
    !# ---
    !# @info
    !# **Brief:** Spectral representation from fourier representation. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: nMax           
    !# assoc leg func degrees (=trunc+1 or +2)
    integer, intent(in) :: mMax           
    !# assoc leg func waves (=trunc+1)
    integer, intent(in) :: mnwv           
    !# spectral coef, real+imag 
!    integer, intent(in) :: xMax           !# longitudes * 2 (real+imag four coef)
    integer, intent(in) :: xmx            
    !# longitudes + 2 (real+imag four coef)
    integer, intent(in) :: yMax           
    !# latitudes (full sphere)
    integer, intent(in) :: yMaxHf         
    !# latitudes (hemisphere)
    integer, intent(in) :: zMax           
    !# verticals
    integer, intent(in) :: dLength(nMax)  
    !# diagonal length (real+imag)
    integer, intent(in) :: lastPrev(nMax) 
    !# last element previous diagonal (real+imag)

    real(kind = p_r8), intent(out) :: spec(mnwv, zMax)      
    !# spectral field
    real(kind = p_r8), intent(in) :: leg (mnwv, yMaxHf)    
    !# associated legendre function
!    real(kind = p_r8), intent(in) :: four(xMax, yMax, zMax) ! full fourier field
    real(kind = p_r8), intent(in) :: four(xmx, yMax, zMax) 
    !# full fourier field

    integer :: j, jj, k, ele, diag

    real(kind = p_r8), dimension(2 * mMax, yMaxHf, zMax) :: fourEven, fourOdd

    ! initialize result

    spec = 0.0_p_r8

    ! use even-odd properties

    do k = 1, zMax
      do j = 1, yMaxHf
        jj = yMax - j + 1
        do ele = 1, 2 * mMax
          fourEven(ele, j, k) = four(ele, j, k) + four(ele, jj, k)
          fourOdd (ele, j, k) = four(ele, j, k) - four(ele, jj, k)
        end do
      end do
    end do

    ! sum odd diagonals (n+m even)

    do k = 1, zMax
      do j = 1, yMaxHf
        do diag = 1, nMax, 2
          do ele = 1, dLength(diag)
            spec(ele + lastPrev(diag), k) = spec(ele + lastPrev(diag), k) + &
              fourEven(ele, j, k) * leg(ele + lastPrev(diag), j)
          end do
        end do
      end do
    end do

    ! sum even diagonals (n+m odd)

    do k = 1, zMax
      do j = 1, yMaxHf
        do diag = 2, nMax, 2
          do ele = 1, dLength(diag)
            spec(ele + lastPrev(diag), k) = spec(ele + lastPrev(diag), k) + &
              fourOdd(ele, j, k) * leg(ele + lastPrev(diag), j)
          end do
        end do
      end do
    end do

  end subroutine sumFour


  subroutine spec2Four2D (mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf, spec, four, der)
    !# Spectral representation to fourier representation 2D
    !# ---
    !# @info
    !# **Brief:** Spectral representation to fourier 2D. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: mnwv2
    integer, intent(in) :: mnwv3
    integer, intent(in) :: mEnd1
    integer, intent(in) :: mEnd2
    integer, intent(in) :: yMaxHf
    real(kind = p_r8), intent(in) :: spec(:, :)
    real(kind = p_r8), intent(out) :: four(:, :, :)

    logical, intent(in), optional :: der

    integer :: s1, s2, f1, f2, f3

    logical :: extended, derivate

    character(len = *), parameter :: h = "**(spec2Four2D)**"

    if (.not. created) then
      write(unit = p_nferr, fmt = '(2A)') h, &
        ' Module not created; invoke InitLegTrans prior to this call'
      stop
    end if

    s1 = size(spec, 1); s2 = size(spec, 2)
    f1 = size(four, 1); f2 = size(four, 2); f3 = size(four, 3)

    extended = .false.
    if (s1 == mnwv2) then
      extended = .false.
    else if (s1 == mnwv3) then
      extended = .true.
    else
      write(unit = p_nferr, fmt = '(2A,I10)') h, &
        ' wrong first dim of spec: ', s1
    end if

    if (s2 /= f3) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' vertical layers of spec and four dissagre :', s2, f3
      stop
    end if

    if (f1 < 2 * mEnd1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' first dimension of four too small: ', f1, 2 * mEnd1
      stop
    end if

    if (f2 /= 2 * yMaxHf) then
      write(unit = p_nferr, fmt = '(2A,I10,A,I10)') h, &
        ' second dimension of four is ', f2, '; should be ', 2 * yMaxHf
      stop
    end if

    if (present(der)) then
      derivate = der
    else
      derivate = .false.
    end if

    if (derivate .and. extended) then
      write(unit = p_nferr, fmt = '(2A)') h, &
        ' derivative cannot be applied to extended gaussian field'
      stop
    end if

    if (extended) then
      call sumSpec (mEnd2, mEnd1, mnwv3, f1, f2, yMaxHf, f3, &
        spec, legExtS2F, four, lenDiagExt, lastPrevDiagExt)
    else if (derivate) then
      call sumSpec (mEnd1, mEnd1, mnwv2, f1, f2, yMaxHf, f3, &
        spec, legDerS2F, four, lenDiag, lastPrevDiag)
    else
      call sumSpec (mEnd1, mEnd1, mnwv2, f1, f2, yMaxHf, f3, &
        spec, legS2F, four, lenDiag, lastPrevDiag)
    end if

  end subroutine spec2Four2D

  
  subroutine spec2Four1D (mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf, spec, four, der)
    !# Spectral representation to fourier representation 1D
    !# ---
    !# @info
    !# **Brief:** Spectral representation to fourier representation 1D. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: mnwv2
    integer, intent(in) :: mnwv3
    integer, intent(in) :: mEnd1
    integer, intent(in) :: mEnd2
    integer, intent(in) :: yMaxHf
    real(kind = p_r8), intent(in) :: spec(:)
    real(kind = p_r8), intent(out) :: four(:, :)

    logical, intent(in), optional :: der

    integer :: s1, f1, f2, f3

    logical :: extended, derivate

    character(len = *), parameter :: h = "**(spec2Four1D)**"

    if (.not. created) then
      write(unit = p_nferr, fmt = '(2A)') h, &
        ' Module not created; invoke InitLegTrans prior to this call'
      stop
    end if

    s1 = size(spec, 1)
    f1 = size(four, 1); f2 = size(four, 2)

    if (s1 == mnwv2) then
      extended = .false.
    else if (s1 == mnwv3) then
      extended = .true.
    else
      ! TODO ... check "else" extended = ?
      extended = .false.
      write(unit = p_nferr, fmt = '(2A,I10)') h, &
        ' wrong first dim of spec : ', s1
    end if

    if (f1 < 2 * mEnd1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' first dimension of four too small: ', f1, 2 * mEnd1
      stop
    end if

    if (f2 /= 2 * yMaxHf) then
      write(unit = p_nferr, fmt = '(2A,I10,A,I10)') h, &
        ' second dimension of four is ', f2, '; should be ', 2 * yMaxHf
      stop
    end if

    if (present(der)) then
      derivate = der
    else
      derivate = .false.
    end if

    if (derivate .and. extended) then
      write(unit = p_nferr, fmt = '(2A)') h, &
        ' derivative cannot be applied to extended gaussian field'
      stop
    end if

    f3 = 1
    if (extended) then
      call sumSpec (mEnd2, mEnd1, mnwv3, f1, f2, yMaxHf, f3, &
        spec, legExtS2F, four, lenDiagExt, lastPrevDiagExt)
    else if (derivate) then
      call sumSpec (mEnd1, mEnd1, mnwv2, f1, f2, yMaxHf, f3, &
        spec, legDerS2F, four, lenDiag, lastPrevDiag)
    else
      call sumSpec (mEnd1, mEnd1, mnwv2, f1, f2, yMaxHf, f3, &
        spec, legS2F, four, lenDiag, lastPrevDiag)
    end if

  end subroutine spec2Four1D


  subroutine four2Spec2D (mnwv2, mEnd1, yMaxHf, four, spec)
    !# Fourier representation to spectral representation 2D
    !# ---
    !# @info
    !# **Brief:** Fourier representation to spectral representation 2D. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: mnwv2
    integer, intent(in) :: mEnd1
    integer, intent(in) :: yMaxHf
    real(kind = p_r8), intent(out) :: spec(:, :)
    real(kind = p_r8), intent(in) :: four(:, :, :)

    integer :: s1, s2, f1, f2, f3

    character(len = *), parameter :: h = "**(four2Spec2D)**"

    if (.not. created) then
      write(unit = p_nferr, fmt = '(2A)') h, &
        ' Module not created; invoke InitLegTrans prior to this call'
      stop
    end if

    s1 = size(spec, 1); s2 = size(spec, 2)
    f1 = size(four, 1); f2 = size(four, 2); f3 = size(four, 3)

    if (s1 /= mnwv2) then
      write(unit = p_nferr, fmt = '(2A,I10)') h, &
        ' wrong first dim of spec: ', s1
    end if

    if (s2 /= f3) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' vertical layers of spec and four dissagre: ', s2, f3
      stop
    end if

    if (f1 < 2 * mEnd1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' first dimension of four too small: ', f1, 2 * mEnd1
      stop
    end if

    if (f2 /= 2 * yMaxHf) then
      write(unit = p_nferr, fmt = '(2A,I10,A,I10)') h, &
        ' second dimension of four is ', f2, '; should be ', 2 * yMaxHf
      stop
    end if

    call sumFour (mEnd1, mEnd1, mnwv2, f1, f2, yMaxHf, f3, &
      spec, legF2S, four, lenDiag, lastPrevDiag)

  end subroutine four2Spec2D


  subroutine four2Spec1D (mnwv2, mEnd1, yMaxHf, four, spec)
    !# Fourier representation to spectral representation 1D
    !# ---
    !# @info
    !# **Brief:** Fourier representation to spectral representation 1D. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    integer, intent(in) :: mnwv2
    integer, intent(in) :: mEnd1
    integer, intent(in) :: yMaxHf
    real(kind = p_r8), intent(out) :: spec(:)
    real(kind = p_r8), intent(in) :: four(:, :)

    integer :: s1, f1, f2, f3

    character(len = *), parameter :: h = "**(four2Spec1D)**"

    if (.not. created) then
      write(unit = p_nferr, fmt = '(2A)') h, &
        ' Module not created; invoke InitLegTrans prior to this call'
      stop
    end if

    s1 = size(spec, 1)
    f1 = size(four, 1); f2 = size(four, 2)

    if (s1 /= mnwv2) then
      write(unit = p_nferr, fmt = '(2A,I10)') h, &
        ' wrong first dim of spec: ', s1
    end if

    if (f1 < 2 * mEnd1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' first dimension of four too small: ', f1, 2 * mEnd1
      stop
    end if

    if (f2 /= 2 * yMaxHf) then
      write(unit = p_nferr, fmt = '(2A,I10,A,I10)') h, &
        ' second dimension of four is ', f2, '; should be ', 2 * yMaxHf
      stop
    end if

    f3 = 1
    call sumFour (mEnd1, mEnd1, mnwv2, f1, f2, yMaxHf, f3, &
      spec, legF2S, four, lenDiag, lastPrevDiag)

  end subroutine four2Spec1D


  subroutine splitTrans3D (full, north, south)
    !# Split Transforms 3D
    !# ---
    !# @info
    !# **Brief:** Split Transforms 3D. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    real(kind = p_r8), intent(in) :: full (:, :, :)
    real(kind = p_r8), intent(out) :: north(:, :, :)
    real(kind = p_r8), intent(out) :: south(:, :, :)

    integer :: if1, in1, is1
    integer :: if2, in2, is2
    integer :: if3, in3, is3
    integer :: i, j, k

    character(len = *), parameter :: h = "**(splitTrans3D)**"

    if1 = size(full, 1); in1 = size(north, 1); is1 = size(south, 1)
    if2 = size(full, 2); in2 = size(north, 2); is2 = size(south, 2)
    if3 = size(full, 3); in3 = size(north, 3); is3 = size(south, 3)

    if (in1 /= is1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' dim 1 of north and south dissagree: ', in1, is1
      stop
    end if
    if (in2 /= is2) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' dim 2 of north and south dissagree: ', in2, is2
      stop
    end if
    if (in3 /= is3) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' dim 3 of north and south dissagree: ', in3, is3
      stop
    end if

    if (in1 < if1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' first dimension of north too small: ', in1, if1
      stop
    end if
    if (if2 /= 2 * in3) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' second dimension of full /= 2*third dimension of north: ', if2, 2 * in3
      stop
    end if
    if (if3 /= in2) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' second dimension of north and third dimension of full dissagree: ', in2, if3
      stop
    end if

    north = 0.0_p_r8
    south = 0.0_p_r8
    do k = 1, in2
      do j = 1, in3
        do i = 1, if1
          north(i, k, j) = full(i, j, k)
          south(i, k, j) = full(i, if2 - j + 1, k)
        end do
      end do
    end do

  end subroutine splitTrans3D


  subroutine splitTrans2D (full, north, south)
    !# Split Transforms 2D
    !# ---
    !# @info
    !# **Brief:** Split Transforms 2D. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    implicit none

    real(kind = p_r8), intent(in) :: full (:, :)
    real(kind = p_r8), intent(out) :: north(:, :)
    real(kind = p_r8), intent(out) :: south(:, :)

    character(len = *), parameter :: h = "**(splitTrans2D)**"

    integer :: if1, in1, is1
    integer :: if2, in2, is2
    integer :: i, j

    if1 = size(full, 1); in1 = size(north, 1); is1 = size(south, 1)
    if2 = size(full, 2); in2 = size(north, 2); is2 = size(south, 2)

    if (in1 /= is1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' dim 1 of north and south dissagree: ', in1, is1
      stop
    end if
    if (in2 /= is2) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' dim 2 of north and south dissagree: ', in2, is2
      stop
    end if

    if (in1 < if1) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' first dimension of north too small: ', in1, if1
      stop
    end if
    if (if2 /= 2 * in2) then
      write(unit = p_nferr, fmt = '(2A,2I10)') h, &
        ' second dimension of full /= 2*second dimension of north: ', if2, 2 * in2
      stop
    end if

    north = 0.0_p_r8
    south = 0.0_p_r8
    do j = 1, in2
      do i = 1, if1
        north(i, j) = full(i, j)
        south(i, j) = full(i, if2 - j + 1)
      end do
    end do

  end subroutine splitTrans2D
  

  subroutine destroyLegendreObjects()
    !# Destroys Legendre Objects
    !# ---
    !# @info
    !# **Brief:** Destroys Legendre Objects. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: abr/2019 </br>
    !# @endinfo
    deallocate(lenDiag)        
    deallocate(lenDiagExt)     
    deallocate(lastPrevDiag)   
    deallocate(lastPrevDiagExt)
    deallocate(la0)            
    deallocate(la1)            
  
    deallocate(eps)            
    deallocate(colRad)         
    deallocate(rCs2)           
    deallocate(wgt)            
    deallocate(gLats)          
    deallocate(legS2F)         
    deallocate(legExtS2F)      
    deallocate(legDerS2F)      
    deallocate(legF2S)         
    deallocate(legDerNS)       
    deallocate(legDerEW)       

    created = .false.

  end subroutine destroyLegendreObjects

end module Mod_LegendreTransform
