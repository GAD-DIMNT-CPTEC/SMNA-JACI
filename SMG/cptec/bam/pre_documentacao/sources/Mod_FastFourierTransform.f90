!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_FastFourierTransform </br></br>
!#
!# **Brief**: Computes direct and inverse FFT(Fast Fourier Transform) transforms
!# of sequences of n real numbers for restricted values of n. Input data for the
!# direct transform and output data from the inverse transform is supposed to
!# represent a periodic function with period n+1. Where n should be in the form
!# n=2*m e m=2**i * 3**j * 5**k, with i>=1 and j, k>=0. </br>
!# 
!# The direct FFT of a sequence of n real numbers fIn(1:n) is a sequence of n+1
!# real numbers fOut(1:n+1) such that
!#   <table align="center">
!#   <tr><td style="width:300px">fOut(2*k+1) = 1/m * SUM(fIn(j+1)*COS(PI*j*k/m)),</td> <td>k=0,..., m   </td></tr>
!#   <tr><td>fOut(2*k+2) = 1/m * SUM(fIn(j+1)*SIN(PI*j*k/m)),</td> <td>k=0,..., m-1 </td></tr>
!#   </table>
!# where both summations are taken for j=0,..., n-1. </br>
!#  
!# The inverse FFT of a sequence of n+1 real numbers fOut(1:n+1) is a sequence 
!# of n real numbers fIn(1:n) such that
!#   <div align="center">fIn(j+1) = 0.5_r8*(fOut(1)*COS(0.0r) + fOut(n+1)*COS(PI*j))
!#   + SUM(fOut(2*k+1)*COS(PI*j*k/m)) + SUM(fOut(2*k+2)*SIN(PI*j*k/m)) </div>
!# where both summations are taken for k=1,..., m-1. </br></br>
!# 
!# USAGE</br>
!# <ul type="disc">
!# <li>Procedure createFFT(n) creates hidden data structures to compute sets of FFTs
!# of length n. Should be invoked prior to the first FFT of that size. </li>
!# 
!# <li>Procedures dirFFT and invFFT compute multiple direct and inverse FFTs of size
!# n. Successive invocations of these procedures use the same hidden data
!# structures created by a single call to createFFT. </li>
!# 
!# <li>destroyFFT destroys internal data structures. Usefull if one has to change the
!# size of FFTs. In that case, createFFT should be invoked with the new size,
!# after the FFT of the previous size was destroied. </li>
!# 
!# <li>This configuration don't allow simultaneous computing of FFT of multiple sizes.
!# A data type is required for that case. </li>
!# </ul> </br>
!# 
!# ALGORITHM </br>
!# 
!# The fundamental FFT algorithm uses 2, 3, 4 and 5 as bases and operates over m
!# complex numbers, expressed as pairs of real numbers. The algorithm requires 
!# m = 2**i * 3**j * 5**k, with i >= 1 and at j, k >= 0 </br>
!# 
!# The direct and the inverse procedures use the fundamental FFT algorithm over
!# complexes builded from even and odd functions of the input data set. Properties
!# of even and odd functions are used to unscramble the desired results (in the
!# direct case) or to scramble input (in the inverse case) from/to the fundamental
!# algorithm. </br>
!# 
!# Periodicity is required to compute even/odd functions. </br></br>
!# 
!# REFERENCES </br>
!#  <ul type="disc">
!# <li>Rader, C. M., "Discrete Fourier Transforms When the Number of Data Samples
!# Is Prime", Proceedings of the IEEE, June 1968, pp 1107-1008. </li>
!# 
!# <li>Singleton, R. C., "Algol Procedures for the Fast Fourier Transform", CACM,
!# November 1969, pp 773-775. </li>
!# 
!# <li>Brigham, E. O., "The Fast Fourier Transform", Prentice-Hall. </li>
!# </ul></br>
!# 
!# CODE HISTORY </br>
!#  
!# Converted from ancient Fortran 77 1-D code by CPTEC in 2001. The conversion
!# procedure transposes input data sets to achieve higher efficiency on memory
!# accesses. Performance increases with the number of transforms computed
!# simultaneously. </br></br>
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
!# <table>
!#  <tr><td style="width:80px"><li>13-11-2004</li></td>  <td style="width:130px">- Jose P. Bonatti </td>  <td>- version: 1.0.0 </td></tr>
!#  <tr><td><li>01-08-2007</li></td>                     <td>- Tomita </td>                               <td>- version: 1.1.1 </td></tr>
!#  <tr><td><li>01-04-2018</li></td>                     <td>- Daniel M. Lamosa</br>- Bárbara Yamada </td>  <td>- version: 2.0.0 </td></tr></li>
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

module Mod_FastFourierTransform

  implicit none
  private

  public :: createFFT, destroyFFT, dirFFT, invFFT

  interface invFFT
    module procedure invD1, invD2, invD3
  end interface

  interface dirFFT
    module procedure dirD1, dirD2, dirD3
  end interface

  ! parameters
  integer, parameter :: p_r8 = selected_real_kind(15)              
  !# kind for 64bits
  integer, parameter :: p_nferr = 0                                
  !# Standard Error Print Out
  integer, parameter :: nBase = 4                                  
  integer, parameter :: base(nBase) = (/ 4, 2, 3, 5 /)             
  integer, parameter :: permutation(nBase) = (/ 2, 3, 1, 4 /)      
  real(kind = p_r8), parameter :: radi = atan(1.0_p_r8) / 45.0_p_r8 
  real(kind = p_r8), parameter :: sin60 = sin(60.0_p_r8 * radi)     
  real(kind = p_r8), parameter :: sin36 = sin(36.0_p_r8 * radi)     
  real(kind = p_r8), parameter :: sin72 = sin(72.0_p_r8 * radi)     
  real(kind = p_r8), parameter :: cos36 = cos(36.0_p_r8 * radi)     
  real(kind = p_r8), parameter :: cos72 = cos(72.0_p_r8 * radi)     

  ! Internal variables
  logical :: created = .false. !
  integer :: nGiven            !

  integer, allocatable :: factors(:) !
  real(kind = p_r8), allocatable :: trigs(:)   !


contains


  subroutine createFFT(xMax)
    !# Creates Fast Fourier Transform
    !# ---
    !# @info
    !# **Brief:** Initializes internal values for Fast Fourier Transforms (FFTs)
    !# of size nIn. Creates hidden data structures to compute sets of FFTs of
    !# length n. Should be invoked prior to the first FFT of that size. 
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: xMax ! 
    character(len = 15), parameter :: h = "**(createFFT)**"  
    !# header

    if(created) then
      write(unit = p_nferr, fmt = '(2a)') h, ' invoked without destroying previous FFT'
      stop
    end if
    created = .true.
    nGiven = xMax
    call factorize(nGiven)
    call trigFactors(nGiven)
  end subroutine createFFT


  subroutine factorize(nIn)
    !# Factorizes
    !# ---
    !# @info
    !# **Brief:** Factorizes nIn/2 in powers of 4, 3, 2, 5, if possible.
    !# Otherwise, stops with error message. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: nIn  !

    character(len = 15), parameter :: h = "**(factorize)**" ! header
    character(len = 15) :: charInt                        ! Character representation of integer
    integer :: powers(nBase)  !
    integer :: nOut           !
    integer :: sumPowers      !
    integer :: ifac           !
    integer :: i              !
    integer :: j              !

    call nextPossibleSize(nIn, nOut, powers)

    if(nIn /= nOut) then
      write(charInt, fmt = '(i15)') nIn
      write(unit = p_nferr, fmt = '(4a)') h, ' FFT size = ', trim(adjustl(charInt)), ' not factorizable '
      write(charInt, fmt = '(i15)') nOut
      write(unit = p_nferr, fmt = '(3A)') h, ' Next factorizable FFT size is ', trim(adjustl(charInt))
      stop
    end if

    sumPowers = sum(powers)
    allocate(factors(sumPowers + 1))
    factors(1) = sumPowers
    ifac = 1
    do i = 1, nBase
      j = permutation(i)
      factors(ifac + 1:ifac + powers(j)) = base(j)
      ifac = ifac + powers(j)
    end do
  end subroutine factorize


  subroutine nextPossibleSize(nIn, nOut, powers)
    !# Next possible FFT size
    !# ---
    !# @info
    !# **Brief:** Next possible Fast Fourier Transform size and its factors. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: nIn           
    integer, intent(out) :: nOut          
    !# positive even
    integer, intent(out) :: powers(nBase) 

    character(len = 22), parameter :: h = "**(nextPossibleSize)**" 
    !# header
    character(len = 15) :: charNIn                               
    !# Character representation of nIn
    integer :: i    
    !# Loop iterator
    integer :: left 
    !# portion of nOut/2 yet to be factorized

    integer, parameter :: limit = huge(nIn) - 1 
    !# Maximum representable integer

    if(nIn <= 0) then
      write(charNIn, fmt = '(i15)') nIn
      write(unit = p_nferr, fmt = '(3a)') h, ' Meaningless FFT size = ', trim(adjustl(charNIn))
      stop
    else if(mod(nIn, 2) == 0) then
      nOut = nIn
    else
      nOut = nIn + 1
    end if

    ! Loop over evens, starting from nOut, looking for next factorizable even/2
    do
      left = nOut / 2
      powers = 0

      ! factorize nOut/2
      do i = 1, nBase
        do
          if(mod(left, base(i)) == 0) then
            powers(i) = powers(i) + 1
            left = left / base(i)
          else
            exit
          end if
        end do
      end do

      if(left == 1) then
        exit
      else if(nOut < limit) then
        nOut = nOut + 2
      else
        write(charNIn, fmt = '(i15)') nIn
        write(unit = p_nferr, fmt = '(4a)') h, ' Next factorizable FFT size > ', &
          trim(adjustl(charNIn)), ' is not representable in this machine'
        stop
      end if
    end do

  end subroutine nextPossibleSize


  subroutine trigFactors (nIn)
    !# TrigFactors
    !# ---
    !# @info
    !# **Brief:** Sin and Cos required to compute Fast Fourier Transform (FFT) of size nIn. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: nIn 

    integer :: nn             
    integer :: nh             
    integer :: i              
    !# Loop iterator
    real(kind = p_r8) :: pi     
    real(kind = p_r8) :: del    
    real(kind = p_r8) :: angle  
    nn = nIn / 2
    nh = (nn + 1) / 2
    allocate(trigs(2 * (nn + nh)))

    pi = 2.0_p_r8 * asin(1.0_p_r8)
    del = (2.0_p_r8 * pi) / real(nn, p_r8)

    do i = 1, 2 * nn, 2
      angle = 0.5_p_r8 * real(i - 1, p_r8) * del
      trigs(i) = cos(angle)
      trigs(i + 1) = sin(angle)
    end do

    del = 0.5_p_r8 * del
    do i = 1, 2 * nh, 2
      angle = 0.5_p_r8 * real(i - 1, p_r8) * del
      trigs(2 * nn + i) = cos(angle)
      trigs(2 * nn + i + 1) = sin(angle)
    end do
  end subroutine trigFactors


  subroutine destroyFFT()
    !# Destroys Fast Fourier Transform
    !# ---
    !# @info
    !# **Brief:** Destroys internal data structures. Usefull if one has to change
    !# the size of Fast Fourier Transforms (FFTs). In that case, createFFT should
    !# be invoked with the new size, after the FFT of the previous size was 
    !# destroied. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    character(len = 16), parameter :: h = "**(destroyFFT)**" 
    !# header
    if(.not. created) then
      write(unit = p_nferr, fmt = '(2a)') h, ' there is no FFT to destroy'
      stop
    end if
    created = .false.
    deallocate(factors)
    deallocate(trigs)
  end subroutine destroyFFT


  subroutine splitFour(fin, a, b, ldin, n, nh, lot)
    !# Split Four
    !# ---
    !# @info
    !# **Brief:** Split Four. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: ldin          
    integer, intent(in) :: n             
    integer, intent(in) :: nh            
    integer, intent(in) :: lot           
    real(kind = p_r8), intent(in) :: fin(ldin, lot) 
    real(kind = p_r8), intent(out) :: a(lot, nh)    
    real(kind = p_r8), intent(out) :: b(lot, nh)    

    integer :: i           
    !# Loop iterator
    integer :: j           
    !# Loop iterator
    real(kind = p_r8) :: c  
    real(kind = p_r8) :: s  

    do i = 1, lot
      a(i, 1) = fin(1, i) + fin(n + 1, i)
      b(i, 1) = fin(1, i) - fin(n + 1, i)
    end do

    do j = 2, (nh + 1) / 2
      c = trigs(n + 2 * j - 1)
      s = trigs(n + 2 * j)
      do i = 1, lot
        a(i, j) = (fin(2 * j - 1, i) + fin(n + 3 - 2 * j, i)) &
          - (s * (fin(2 * j - 1, i) - fin(n + 3 - 2 * j, i)) &
            + c * (fin(2 * j, i) + fin(n + 4 - 2 * j, i)))

        a(i, nh + 2 - j) = (fin(2 * j - 1, i) + fin(n + 3 - 2 * j, i)) &
          + (s * (fin(2 * j - 1, i) - fin(n + 3 - 2 * j, i)) &
            + c * (fin(2 * j, i) + fin(n + 4 - 2 * j, i)))

        b(i, j) = (c * (fin(2 * j - 1, i) - fin(n + 3 - 2 * j, i))  &
          - s * (fin(2 * j, i) + fin(n + 4 - 2 * j, i))) &
          + (fin(2 * j, i) - fin(n + 4 - 2 * j, i))

        b(i, nh + 2 - j) = (c * (fin(2 * j - 1, i) - fin(n + 3 - 2 * j, i))  &
          - s * (fin(2 * j, i) + fin(n + 4 - 2 * j, i))) &
          - (fin(2 * j, i) - fin(n + 4 - 2 * j, i))
      end do
    end do

    if((nh>=2) .and. (mod(nh, 2)==0)) then
      do i = 1, lot
        a(i, nh / 2 + 1) = 2.0_p_r8 * fin(nh + 1, i)
        b(i, nh / 2 + 1) = -2.0_p_r8 * fin(nh + 2, i)
      end do
    end if
  end subroutine splitFour


  subroutine joinFour(a, b, fout, ldout, n, nh, lot)
    !# Join Four
    !# ---
    !# @info
    !# **Brief:** Join Four. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: ldout           
    integer, intent(in) :: n               
    integer, intent(in) :: nh              
    integer, intent(in) :: lot             
    real(kind = p_r8), intent(in) :: a(lot, nh)       
    real(kind = p_r8), intent(in) :: b(lot, nh)       
    real(kind = p_r8), intent(out) :: fout(ldout, lot) 

    integer :: i             
    !# Loop iterator
    integer :: j             
    !# Loop iterator
    real(kind = p_r8) :: scalR 
    !# scale
    real(kind = p_r8) :: scalH 
    !# half of scale
    real(kind = p_r8) :: c     
    real(kind = p_r8) :: s     

    scalR = 1.0_p_r8 / real(n, p_r8)
    scalH = 0.5_p_r8 * scalR

    do i = 1, lot
      fout(1, i) = scalR * (a(i, 1) + b(i, 1))
      fout(n + 1, i) = scalR * (a(i, 1) - b(i, 1))
      fout(2, i) = 0.0_p_r8
    end do

    do j = 2, (nh + 1) / 2
      c = trigs(n + 2 * j - 1)
      s = trigs(n + 2 * j)
      do i = 1, lot
        fout(2 * j - 1, i) = scalH * ((a(i, j) + a(i, nh + 2 - j))  &
          + (c * (b(i, j) + b(i, nh + 2 - j))  &
            + s * (a(i, j) - a(i, nh + 2 - j))))

        fout(n + 3 - 2 * j, i) = scalH * ((a(i, j) + a(i, nh + 2 - j))  &
          - (c * (b(i, j) + b(i, nh + 2 - j))  &
            + s * (a(i, j) - a(i, nh + 2 - j))))

        fout(2 * j, i) = scalH * ((c * (a(i, j) - a(i, nh + 2 - j))  &
          - s * (b(i, j) + b(i, nh + 2 - j))) &
          + (b(i, nh + 2 - j) - b(i, j)))

        fout(n + 4 - 2 * j, i) = scalH * ((c * (a(i, j) - a(i, nh + 2 - j)) &
          - s * (b(i, j) + b(i, nh + 2 - j)))&
          - (b(i, nh + 2 - j) - b(i, j)))
      end do
    end do

    if((nh>=2) .and. (mod(nh, 2)==0)) then
      do i = 1, lot
        fout(nh + 1, i) = scalR * a(i, nh / 2 + 1)
        fout(nh + 2, i) = -scalR * b(i, nh / 2 + 1)
      end do
    end if

    fout(n + 2:ldout, :) = 0.0_p_r8

  end subroutine joinFour


  subroutine splitGaus(fin, a, b, ldin, nh, lot)
    !# Split Gaus
    !# ---
    !# @info
    !# **Brief:** Split Gaus. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: ldin          
    integer, intent(in) :: nh            
    integer, intent(in) :: lot           
    real(kind = p_r8), intent(in) :: fin(ldin, lot) 
    real(kind = p_r8), intent(out) :: a(lot, nh)    
    real(kind = p_r8), intent(out) :: b(lot, nh)    

    integer :: i 
    !# Loop iterator
    integer :: j 
    !# Loop iterator

    do j = 1, nh
      do i = 1, lot
        a(i, j) = fin(2 * j - 1, i)
        b(i, j) = fin(2 * j, i)
      end do
    end do
  end subroutine splitGaus


  subroutine joinGaus (a, b, fout, ldout, nh, lot)
    !# Join Gaussian
    !# ---
    !# @info
    !# **Brief:** Join Gaussian. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: ldout           
    integer, intent(in) :: nh              
    integer, intent(in) :: lot             
    real(kind = p_r8), intent(out) :: fout(ldout, lot) 
    real(kind = p_r8), intent(in) :: a(lot, nh)       
    real(kind = p_r8), intent(in) :: b(lot, nh)       

    integer :: i 
    !# Loop iterator
    integer :: j 
    !# Loop iterator

    do j = 1, nh
      do i = 1, lot
        fout(2 * j - 1, i) = a(i, j)
        fout(2 * j, i) = b(i, j)
      end do
    end do

    fout(2 * nh + 1:ldout, :) = 0.0_p_r8

  end subroutine joinGaus


  subroutine onePass(a, b, c, d, lot, nh, ifac, la)
    !# One Pass
    !# ---
    !# @info
    !# **Brief:** One Pass. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: lot        
    integer, intent(in) :: nh         
    !# PROD(factor(1:K))
    integer, intent(in) :: ifac       
    !# factor(k)
    integer, intent(in) :: la         
    !# PROD(factor(1:k-1))
    real(kind = p_r8), intent(in) :: a(lot, nh)  
    real(kind = p_r8), intent(in) :: b(lot, nh)  
    real(kind = p_r8), intent(out) :: c(lot, nh) 
    real(kind = p_r8), intent(out) :: d(lot, nh) 

    integer :: m     
    integer :: jump  
    integer :: i     
    !# Loop iterator
    integer :: j     
    !# Loop iterator
    integer :: k     
    !# Loop iterator
    integer :: ia    
    integer :: ja    
    integer :: ib    
    integer :: jb    
    integer :: kb    
    integer :: ic    
    integer :: jc    
    integer :: kc    
    integer :: id    
    integer :: jd    
    integer :: kd    
    integer :: ie    
    integer :: je    
    integer :: ke    

    real(kind = p_r8) :: c1      
    real(kind = p_r8) :: s1      
    real(kind = p_r8) :: c2      
    real(kind = p_r8) :: s2      
    real(kind = p_r8) :: c3      
    real(kind = p_r8) :: s3      
    real(kind = p_r8) :: c4      
    real(kind = p_r8) :: s4      
    real(kind = p_r8) :: wka     
    real(kind = p_r8) :: wkb     
    real(kind = p_r8) :: wksina  
    real(kind = p_r8) :: wksinb  
    real(kind = p_r8) :: wkaacp  
    real(kind = p_r8) :: wkbacp  
    real(kind = p_r8) :: wkaacm  
    real(kind = p_r8) :: wkbacm  

    m = nh / ifac
    jump = (ifac - 1) * la

    ia = 0
    ib = m
    ic = 2 * m
    id = 3 * m
    ie = 4 * m

    ja = 0
    jb = la
    jc = 2 * la
    jd = 3 * la
    je = 4 * la

    if(ifac == 2) then
      do j = 1, la
        do i = 1, lot
          c(i, j + ja) = a(i, j + ia) + a(i, j + ib)
          c(i, j + jb) = a(i, j + ia) - a(i, j + ib)
          d(i, j + ja) = b(i, j + ia) + b(i, j + ib)
          d(i, j + jb) = b(i, j + ia) - b(i, j + ib)
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        c1 = trigs(kb + 1)
        s1 = trigs(kb + 2)
        ja = ja + jump
        jb = jb + jump
        do j = k + 1, k + la
          do i = 1, lot
            wka = a(i, j + ia) - a(i, j + ib)
            c(i, j + ja) = a(i, j + ia) + a(i, j + ib)
            wkb = b(i, j + ia) - b(i, j + ib)
            d(i, j + ja) = b(i, j + ia) + b(i, j + ib)
            c(i, j + jb) = c1 * wka - s1 * wkb
            d(i, j + jb) = s1 * wka + c1 * wkb
          end do
        end do
      end do
    elseif(ifac == 3) then
      do j = 1, la
        do i = 1, lot
          wka = a(i, j + ib) + a(i, j + ic)
          wksina = sin60 * (a(i, j + ib) - a(i, j + ic))
          wkb = b(i, j + ib) + b(i, j + ic)
          wksinb = sin60 * (b(i, j + ib) - b(i, j + ic))
          c(i, j + ja) = a(i, j + ia) + wka
          c(i, j + jb) = (a(i, j + ia) - 0.5_p_r8 * wka) - wksinb
          c(i, j + jc) = (a(i, j + ia) - 0.5_p_r8 * wka) + wksinb
          d(i, j + ja) = b(i, j + ia) + wkb
          d(i, j + jb) = (b(i, j + ia) - 0.5_p_r8 * wkb) + wksina
          d(i, j + jc) = (b(i, j + ia) - 0.5_p_r8 * wkb) - wksina
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        kc = kb + kb
        c1 = trigs(kb + 1)
        s1 = trigs(kb + 2)
        c2 = trigs(kc + 1)
        s2 = trigs(kc + 2)
        ja = ja + jump
        jb = jb + jump
        jc = jc + jump
        do j = k + 1, k + la
          do i = 1, lot
            wka = a(i, j + ib) + a(i, j + ic)
            wksina = sin60 * (a(i, j + ib) - a(i, j + ic))
            wkb = b(i, j + ib) + b(i, j + ic)
            wksinb = sin60 * (b(i, j + ib) - b(i, j + ic))
            c(i, j + ja) = a(i, j + ia) + wka
            d(i, j + ja) = b(i, j + ia) + wkb
            c(i, j + jb) = c1 * ((a(i, j + ia) - 0.5_p_r8 * wka) - wksinb) &
              - s1 * ((b(i, j + ia) - 0.5_p_r8 * wkb) + wksina)
            d(i, j + jb) = s1 * ((a(i, j + ia) - 0.5_p_r8 * wka) - wksinb) &
              + c1 * ((b(i, j + ia) - 0.5_p_r8 * wkb) + wksina)
            c(i, j + jc) = c2 * ((a(i, j + ia) - 0.5_p_r8 * wka) + wksinb) &
              - s2 * ((b(i, j + ia) - 0.5_p_r8 * wkb) - wksina)
            d(i, j + jc) = s2 * ((a(i, j + ia) - 0.5_p_r8 * wka) + wksinb) &
              + c2 * ((b(i, j + ia) - 0.5_p_r8 * wkb) - wksina)
          end do
        end do
      end do
    elseif(ifac == 4) then
      do j = 1, la
        do i = 1, lot
          wkaacp = a(i, j + ia) + a(i, j + ic)
          wkaacm = a(i, j + ia) - a(i, j + ic)
          wkbacp = b(i, j + ia) + b(i, j + ic)
          wkbacm = b(i, j + ia) - b(i, j + ic)
          c(i, j + ja) = wkaacp + (a(i, j + ib) + a(i, j + id))
          c(i, j + jc) = wkaacp - (a(i, j + ib) + a(i, j + id))
          d(i, j + jb) = wkbacm + (a(i, j + ib) - a(i, j + id))
          d(i, j + jd) = wkbacm - (a(i, j + ib) - a(i, j + id))
          d(i, j + ja) = wkbacp + (b(i, j + ib) + b(i, j + id))
          d(i, j + jc) = wkbacp - (b(i, j + ib) + b(i, j + id))
          c(i, j + jb) = wkaacm - (b(i, j + ib) - b(i, j + id))
          c(i, j + jd) = wkaacm + (b(i, j + ib) - b(i, j + id))
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        kc = kb + kb
        kd = kc + kb
        c1 = trigs(kb + 1)
        s1 = trigs(kb + 2)
        c2 = trigs(kc + 1)
        s2 = trigs(kc + 2)
        c3 = trigs(kd + 1)
        s3 = trigs(kd + 2)
        ja = ja + jump
        jb = jb + jump
        jc = jc + jump
        jd = jd + jump
        do j = k + 1, k + la
          do i = 1, lot
            wkaacp = a(i, j + ia) + a(i, j + ic)
            wkbacp = b(i, j + ia) + b(i, j + ic)
            wkaacm = a(i, j + ia) - a(i, j + ic)
            wkbacm = b(i, j + ia) - b(i, j + ic)
            c(i, j + ja) = wkaacp + (a(i, j + ib) + a(i, j + id))
            d(i, j + ja) = wkbacp + (b(i, j + ib) + b(i, j + id))
            c(i, j + jc) = c2 * (wkaacp - (a(i, j + ib) + a(i, j + id))) &
              - s2 * (wkbacp - (b(i, j + ib) + b(i, j + id)))
            d(i, j + jc) = s2 * (wkaacp - (a(i, j + ib) + a(i, j + id))) &
              + c2 * (wkbacp - (b(i, j + ib) + b(i, j + id)))
            c(i, j + jb) = c1 * (wkaacm - (b(i, j + ib) - b(i, j + id))) &
              - s1 * (wkbacm + (a(i, j + ib) - a(i, j + id)))
            d(i, j + jb) = s1 * (wkaacm - (b(i, j + ib) - b(i, j + id))) &
              + c1 * (wkbacm + (a(i, j + ib) - a(i, j + id)))
            c(i, j + jd) = c3 * (wkaacm + (b(i, j + ib) - b(i, j + id))) &
              - s3 * (wkbacm - (a(i, j + ib) - a(i, j + id)))
            d(i, j + jd) = s3 * (wkaacm + (b(i, j + ib) - b(i, j + id))) &
              + c3 * (wkbacm - (a(i, j + ib) - a(i, j + id)))
          end do
        end do
      end do
    elseif(ifac == 5) then
      do j = 1, la
        do i = 1, lot
          c(i, j + ja) = a(i, j + ia) + (a(i, j + ib) + a(i, j + ie)) + (a(i, j + ic) + a(i, j + id))
          d(i, j + ja) = b(i, j + ia) + (b(i, j + ib) + b(i, j + ie)) + (b(i, j + ic) + b(i, j + id))
          c(i, j + jb) = (a(i, j + ia)                       &
            + cos72 * (a(i, j + ib) + a(i, j + ie)) &
            - cos36 * (a(i, j + ic) + a(i, j + id)))&
            - (sin72 * (b(i, j + ib) - b(i, j + ie)) &
              + sin36 * (b(i, j + ic) - b(i, j + id)))
          c(i, j + je) = (a(i, j + ia)                       &
            + cos72 * (a(i, j + ib) + a(i, j + ie)) &
            - cos36 * (a(i, j + ic) + a(i, j + id)))&
            + (sin72 * (b(i, j + ib) - b(i, j + ie)) &
              + sin36 * (b(i, j + ic) - b(i, j + id)))
          d(i, j + jb) = (b(i, j + ia)                       &
            + cos72 * (b(i, j + ib) + b(i, j + ie)) &
            - cos36 * (b(i, j + ic) + b(i, j + id)))&
            + (sin72 * (a(i, j + ib) - a(i, j + ie)) &
              + sin36 * (a(i, j + ic) - a(i, j + id)))
          d(i, j + je) = (b(i, j + ia)                       &
            + cos72 * (b(i, j + ib) + b(i, j + ie)) &
            - cos36 * (b(i, j + ic) + b(i, j + id)))&
            - (sin72 * (a(i, j + ib) - a(i, j + ie)) &
              + sin36 * (a(i, j + ic) - a(i, j + id)))
          c(i, j + jc) = (a(i, j + ia)                       &
            - cos36 * (a(i, j + ib) + a(i, j + ie)) &
            + cos72 * (a(i, j + ic) + a(i, j + id)))&
            - (sin36 * (b(i, j + ib) - b(i, j + ie)) &
              - sin72 * (b(i, j + ic) - b(i, j + id)))
          c(i, j + jd) = (a(i, j + ia)                       &
            - cos36 * (a(i, j + ib) + a(i, j + ie)) &
            + cos72 * (a(i, j + ic) + a(i, j + id)))&
            + (sin36 * (b(i, j + ib) - b(i, j + ie)) &
              - sin72 * (b(i, j + ic) - b(i, j + id)))
          d(i, j + jc) = (b(i, j + ia)                       &
            - cos36 * (b(i, j + ib) + b(i, j + ie)) &
            + cos72 * (b(i, j + ic) + b(i, j + id)))&
            + (sin36 * (a(i, j + ib) - a(i, j + ie)) &
              - sin72 * (a(i, j + ic) - a(i, j + id)))
          d(i, j + jd) = (b(i, j + ia)                       &
            - cos36 * (b(i, j + ib) + b(i, j + ie)) &
            + cos72 * (b(i, j + ic) + b(i, j + id)))&
            - (sin36 * (a(i, j + ib) - a(i, j + ie)) &
              - sin72 * (a(i, j + ic) - a(i, j + id)))
        end do
      end do
      do k = la, m - 1, la
        kb = k + k
        kc = kb + kb
        kd = kc + kb
        ke = kd + kb
        c1 = trigs(kb + 1)
        s1 = trigs(kb + 2)
        c2 = trigs(kc + 1)
        s2 = trigs(kc + 2)
        c3 = trigs(kd + 1)
        s3 = trigs(kd + 2)
        c4 = trigs(ke + 1)
        s4 = trigs(ke + 2)
        ja = ja + jump
        jb = jb + jump
        jc = jc + jump
        jd = jd + jump
        je = je + jump
        do j = k + 1, k + la
          do i = 1, lot
            c(i, j + ja) = a(i, j + ia)               &
              + (a(i, j + ib) + a(i, j + ie)) &
              + (a(i, j + ic) + a(i, j + id))
            d(i, j + ja) = b(i, j + ia)               &
              + (b(i, j + ib) + b(i, j + ie)) &
              + (b(i, j + ic) + b(i, j + id))
            c(i, j + jb) = c1 * ((a(i, j + ia)                   &
              + cos72 * (a(i, j + ib) + a(i, j + ie))   &
              - cos36 * (a(i, j + ic) + a(i, j + id)))  &
              - (sin72 * (b(i, j + ib) - b(i, j + ie))   &
                + sin36 * (b(i, j + ic) - b(i, j + id)))) &
              - s1 * ((b(i, j + ia)                   &
                + cos72 * (b(i, j + ib) + b(i, j + ie))   &
                - cos36 * (b(i, j + ic) + b(i, j + id)))  &
                + (sin72 * (a(i, j + ib) - a(i, j + ie))   &
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + jb) = s1 * ((a(i, j + ia)                   &
              + cos72 * (a(i, j + ib) + a(i, j + ie))   &
              - cos36 * (a(i, j + ic) + a(i, j + id)))  &
              - (sin72 * (b(i, j + ib) - b(i, j + ie))   &
                + sin36 * (b(i, j + ic) - b(i, j + id)))) &
              + c1 * ((b(i, j + ia)                   &
                + cos72 * (b(i, j + ib) + b(i, j + ie))   &
                - cos36 * (b(i, j + ic) + b(i, j + id)))  &
                + (sin72 * (a(i, j + ib) - a(i, j + ie))   &
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            c(i, j + je) = c4 * ((a(i, j + ia)                   &
              + cos72 * (a(i, j + ib) + a(i, j + ie))   &
              - cos36 * (a(i, j + ic) + a(i, j + id)))  &
              + (sin72 * (b(i, j + ib) - b(i, j + ie))   &
                + sin36 * (b(i, j + ic) - b(i, j + id)))) &
              - s4 * ((b(i, j + ia)                   &
                + cos72 * (b(i, j + ib) + b(i, j + ie))   &
                - cos36 * (b(i, j + ic) + b(i, j + id)))  &
                - (sin72 * (a(i, j + ib) - a(i, j + ie))   &
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + je) = s4 * ((a(i, j + ia)                   &
              + cos72 * (a(i, j + ib) + a(i, j + ie))   &
              - cos36 * (a(i, j + ic) + a(i, j + id)))  &
              + (sin72 * (b(i, j + ib) - b(i, j + ie))   &
                + sin36 * (b(i, j + ic) - b(i, j + id)))) &
              + c4 * ((b(i, j + ia)                   &
                + cos72 * (b(i, j + ib) + b(i, j + ie))   &
                - cos36 * (b(i, j + ic) + b(i, j + id)))  &
                - (sin72 * (a(i, j + ib) - a(i, j + ie))   &
                  + sin36 * (a(i, j + ic) - a(i, j + id))))
            c(i, j + jc) = c2 * ((a(i, j + ia)                   &
              - cos36 * (a(i, j + ib) + a(i, j + ie))   &
              + cos72 * (a(i, j + ic) + a(i, j + id)))  &
              - (sin36 * (b(i, j + ib) - b(i, j + ie))   &
                - sin72 * (b(i, j + ic) - b(i, j + id)))) &
              - s2 * ((b(i, j + ia)                   &
                - cos36 * (b(i, j + ib) + b(i, j + ie))   &
                + cos72 * (b(i, j + ic) + b(i, j + id)))  &
                + (sin36 * (a(i, j + ib) - a(i, j + ie))   &
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + jc) = s2 * ((a(i, j + ia)                   &
              - cos36 * (a(i, j + ib) + a(i, j + ie))   &
              + cos72 * (a(i, j + ic) + a(i, j + id)))  &
              - (sin36 * (b(i, j + ib) - b(i, j + ie))   &
                - sin72 * (b(i, j + ic) - b(i, j + id)))) &
              + c2 * ((b(i, j + ia)                   &
                - cos36 * (b(i, j + ib) + b(i, j + ie))   &
                + cos72 * (b(i, j + ic) + b(i, j + id)))  &
                + (sin36 * (a(i, j + ib) - a(i, j + ie))   &
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
            c(i, j + jd) = c3 * ((a(i, j + ia)                   &
              - cos36 * (a(i, j + ib) + a(i, j + ie))   &
              + cos72 * (a(i, j + ic) + a(i, j + id)))  &
              + (sin36 * (b(i, j + ib) - b(i, j + ie))   &
                - sin72 * (b(i, j + ic) - b(i, j + id)))) &
              - s3 * ((b(i, j + ia)                   &
                - cos36 * (b(i, j + ib) + b(i, j + ie))   &
                + cos72 * (b(i, j + ic) + b(i, j + id)))  &
                - (sin36 * (a(i, j + ib) - a(i, j + ie))   &
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
            d(i, j + jd) = s3 * ((a(i, j + ia)                   &
              - cos36 * (a(i, j + ib) + a(i, j + ie))   &
              + cos72 * (a(i, j + ic) + a(i, j + id)))  &
              + (sin36 * (b(i, j + ib) - b(i, j + ie))   &
                - sin72 * (b(i, j + ic) - b(i, j + id)))) &
              + c3 * ((b(i, j + ia)                   &
                - cos36 * (b(i, j + ib) + b(i, j + ie))   &
                + cos72 * (b(i, j + ic) + b(i, j + id)))  &
                - (sin36 * (a(i, j + ib) - a(i, j + ie))   &
                  - sin72 * (a(i, j + ic) - a(i, j + id))))
          end do
        end do
      end do
    endif

  end subroutine onePass


  subroutine dir(fin, fout, ldin, ldout, n, lot)
    !# Computes multiple direct FFT
    !# ---
    !# @info
    !# **Brief:** Procedures dirFFT compute multiple direct Fast Fourier Transforms
    !# (FFTs) of size n. Successive invocations of these procedures use the same
    !# hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: ldin  !
    integer, intent(in) :: ldout !
    integer, intent(in) :: n     !
    integer, intent(in) :: lot   ! 
    real(kind = p_r8), intent(in) :: fin(ldin, lot)   !
    real(kind = p_r8), intent(out) :: fout(ldout, lot) !

    character(len = 12), parameter :: h = "**(dir)**" !
    character(len = 15) :: charInt1      !
    character(len = 15) :: charInt2      !

    logical :: ab2cd !
    integer :: nh    !
    integer :: nfax  !
    integer :: la    !
    integer :: k     !

    real(kind = p_r8), dimension(lot, n / 2) :: a !
    real(kind = p_r8), dimension(lot, n / 2) :: b !
    real(kind = p_r8), dimension(lot, n / 2) :: c !
    real(kind = p_r8), dimension(lot, n / 2) :: d !

    if(.not. created) then
      write(unit = p_nferr, fmt = '(2a)') h, ' FFT was not created'
      stop
    else if(n /= nGiven) then
      write(charInt1, fmt = '(i15)') n
      write(charInt2, fmt = '(i15)') nGiven
      write(unit = p_nferr, fmt = '(4a)') h, &
        ' FFT invoked with size ', trim(adjustl(charInt1)), &
        ' but created with size ', trim(adjustl(charInt2))
      stop
    else if(ldout < n + 1) then
      write(charInt1, fmt = '(i15)') ldout
      write(charInt2, fmt = '(i15)') n + 1
      write(unit = p_nferr, fmt = '(4a)') h, &
        ' Output field has first dimension ', trim(adjustl(charInt1)), &
        '; should be at least ', trim(adjustl(charInt2))
      stop
    else if(ldin < n) then
      write(charInt1, fmt = '(i15)') ldin
      write(charInt2, fmt = '(i15)') n
      write(unit = p_nferr, fmt = '(4a)') h, &
        ' Input field has first dimension ', trim(adjustl(charInt1)), &
        '; should be at least ', trim(adjustl(charInt2))
      stop
    end if

    nfax = factors(1)
    nh = n / 2

    call splitGaus(fin, a, b, ldin, nh, lot)

    la = 1
    ab2cd = .true.
    do k = 1, nfax
      if(ab2cd) then
        call onePass(a, b, c, d, lot, nh, factors(k + 1), la)
        ab2cd = .false.
      else
        call onePass(c, d, a, b, lot, nh, factors(k + 1), la)
        ab2cd = .true.
      end if
      la = la * factors(k + 1)
    end do

    if(ab2cd) then
      call joinFour(a, b, fout, ldout, n, nh, lot)
    else
      call joinFour(c, d, fout, ldout, n, nh, lot)
    end if

  end subroutine dir


  subroutine dirD1(fin, fout)
    !# Computes multiple direct FFT
    !# ---
    !# @info
    !# **Brief:** Procedures dirFFT compute multiple direct Fast Fourier Transforms
    !# (FFTs) of size n. Successive invocations of these procedures use the same
    !# hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in) :: fin(:)  !
    real(kind = p_r8), intent(out) :: fout(:) !

    integer :: din  !
    integer :: dout !

    din = size(fin)
    dout = size(fout)
    call dir(fin, fout, din, dout, nGiven, 1)

  end subroutine dirD1


  subroutine dirD2(fin, fout)
    !# Computes multiple direct FFT
    !# ---
    !# @info
    !# **Brief:** Procedures dirFFT compute multiple direct Fast Fourier Transforms
    !# (FFTs) of size n. Successive invocations of these procedures use the same
    !# hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in) :: fin(:, :)  !
    real(kind = p_r8), intent(out) :: fout(:, :) !

    character(len = 11), parameter :: h = "**(dirD2)**" !

    integer :: d2in  !
    integer :: d2out !
    integer :: din   !
    integer :: dout  !

    d2in = size(fin, 2)
    d2out = size(fout, 2)
    if(d2in /= d2out) then
      write(unit = p_nferr, fmt = '(2a,2i10)') h, &
        ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
      stop
    end if

    din = size(fin, 1)
    dout = size(fout, 1)
    call dir(fin, fout, din, dout, nGiven, d2in)

  end subroutine dirD2


  subroutine dirD3(fin, fout)
    !# Computes multiple direct FFT
    !# ---
    !# @info
    !# **Brief:** Procedures dirFFT compute multiple direct Fast Fourier Transforms
    !# (FFTs) of size n. Successive invocations of these procedures use the same
    !# hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in) :: fin(:, :, :)  !
    real(kind = p_r8), intent(out) :: fout(:, :, :) !

    character(len = 11), parameter :: h = "**(dirD3)**" !

    integer :: d2in  !
    integer :: d2out !
    integer :: d3in  !
    integer :: d3out !
    integer :: din   !
    integer :: dout  !
    integer :: dio   !

    d2in = size(fin, 2)
    d2out = size(fout, 2)
    if(d2in /= d2out) then
      write(unit = p_nferr, fmt = '(2a,2i10)') h, &
        ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
      stop
    end if

    d3in = size(fin, 3)
    d3out = size(fout, 3)
    if(d3in /= d3out) then
      write(unit = p_nferr, fmt = '(2a,2i10)') h, &
        ' Error: dim 3 of fin, fout, differ: ', d3in, d3out
      stop
    end if

    din = size(fin, 1)
    dout = size(fout, 1)
    dio = d2in * d3in
    call dir(fin, fout, din, dout, nGiven, dio)

  end subroutine dirD3


  subroutine inv(fin, fout, ldin, ldout, n, lot)
    !# Computes multiple inverse FFT
    !# ---
    !# @info
    !# **Brief:** Procedures invFFT compute multiple inverse Fast Fourier
    !# Transforms (FFTs) of size n. Successive invocations of these procedures
    !# use the same hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    integer, intent(in) :: ldin                     !
    integer, intent(in) :: ldout                    !
    integer, intent(in) :: n                        !
    integer, intent(in) :: lot                      !
    real(kind = p_r8), intent(in) :: fin (ldin, lot)  !
    real(kind = p_r8), intent(out) :: fout(ldout, lot) !

    character(len = 12), parameter :: h = "**(inv)**" !
    character(len = 15) :: charInt1                   !
    character(len = 15) :: charInt2                   !

    logical :: ab2cd !
    integer :: nh    !
    integer :: nfax  !
    integer :: la    !
    integer :: k     !

    real(kind = p_r8), dimension(lot, n / 2) :: a !
    real(kind = p_r8), dimension(lot, n / 2) :: b !
    real(kind = p_r8), dimension(lot, n / 2) :: c !
    real(kind = p_r8), dimension(lot, n / 2) :: d !

    if(.not. created) then
      write(unit = p_nferr, fmt = '(4a)') h, ' FFT was not created'
      stop
    else if(n /= nGiven) then
      write(charInt1, fmt = '(i15)') n
      write(charInt2, fmt = '(i15)') nGiven
      write(unit = p_nferr, fmt = '(4a)') h, &
        ' FFT invoked with size ', trim(adjustl(charInt1)), &
        ' but created with size ', trim(adjustl(charInt2))
      stop
    else if(ldin < n + 1) then
      write(charInt1, fmt = '(i15)') ldin
      write(charInt2, fmt = '(i15)') n + 1
      write(unit = p_nferr, fmt = '(4a)') h, &
        ' Input field has first dimension ', trim(adjustl(charInt1)), &
        '; should be at least ', trim(adjustl(charInt2))
      stop
    else if(ldout < n) then
      write(charInt1, fmt = '(i15)') ldout
      write(charInt2, fmt = '(i15)') n
      write(unit = p_nferr, fmt = '(4a)') h, &
        ' Output field has first dimension ', trim(adjustl(charInt1)), &
        '; should be at least ', trim(adjustl(charInt2))
      stop
    end if

    nfax = factors(1)
    nh = n / 2

    call splitFour(fin, a, b, ldin, n, nh, lot)

    la = 1
    ab2cd = .true.
    do k = 1, nfax
      if(ab2cd) then
        call onePass(a, b, c, d, lot, nh, factors(k + 1), la)
        ab2cd = .false.
      else
        call onePass(c, d, a, b, lot, nh, factors(k + 1), la)
        ab2cd = .true.
      end if
      la = la * factors(k + 1)
    end do

    if(ab2cd) then
      call joinGaus(a, b, fout, ldout, nh, lot)
    else
      call joinGaus(c, d, fout, ldout, nh, lot)
    end if

  end subroutine inv


  subroutine invD1(fin, fout)
    !# Computes multiple inverse FFT
    !# ---
    !# @info
    !# **Brief:** Procedures invFFT compute multiple inverse Fast Fourier
    !# Transforms (FFTs) of size n. Successive invocations of these procedures
    !# use the same hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in) :: fin(:)  !
    real(kind = p_r8), intent(out) :: fout(:) !

    integer :: din  !
    integer :: dout !

    din = size(fin)
    dout = size(fout)
    call inv(fin, fout, din, dout, nGiven, 1)

  end subroutine invD1


  subroutine invD2(fin, fout)
    !# Computes multiple inverse FFT
    !# ---
    !# @info
    !# **Brief:** Procedures invFFT compute multiple inverse Fast Fourier
    !# Transforms (FFTs) of size n. Successive invocations of these procedures
    !# use the same hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in) :: fin(:, :)  !
    real(kind = p_r8), intent(out) :: fout(:, :) !

    character(len = 11), parameter :: h = "**(invD2)**" !

    integer :: d2in  !
    integer :: d2out !
    integer :: din   !
    integer :: dout  !

    d2in = size(fin, 2)
    d2out = size(fout, 2)
    if(d2in /= d2out) then
      write(unit = p_nferr, fmt = '(2a,2i10)') h, &
        ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
      stop
    end if

    din = size(fin, 1)
    dout = size(fout, 1)
    call inv(fin, fout, din, dout, nGiven, d2in)

  end subroutine invD2


  subroutine invD3(fin, fout)
    !# invD3
    !# ---
    !# @info
    !# **Brief:** Procedures invFFT compute multiple inverse Fast Fourier
    !# Transforms (FFTs) of size n. Successive invocations of these procedures
    !# use the same hidden data structures created by a single call to createFFT. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Barbara A. G. P. Yamada </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    real(kind = p_r8), intent(in) :: fin(:, :, :)  !
    real(kind = p_r8), intent(out) :: fout(:, :, :) !

    character(len = 11), parameter :: h = "**(invD3)**" !

    integer :: d2in  !
    integer :: d2out !
    integer :: d3in  !
    integer :: d3out !
    integer :: din   !
    integer :: dout  !
    integer :: dio   !

    d2in = size(fin, 2)
    d2out = size(fout, 2)
    if(d2in /= d2out) then
      write(unit = p_nferr, fmt = '(2a,2i10)') h, &
        ' Error: dim 2 of fin, fout, differ: ', d2in, d2out
      stop
    end if

    d3in = size(fin, 3)
    d3out = size(fout, 3)
    if(d3in /= d3out) then
      write(unit = p_nferr, fmt = '(2a,2i10)') h, &
        ' Error: dim 3 of fin, fout, differ: ', d3in, d3out
      stop
    end if

    din = size(fin, 1)
    dout = size(fout, 1)
    dio = d2in * d3in
    call inv(fin, fout, din, dout, nGiven, dio)

  end subroutine invD3

end module Mod_FastFourierTransform
