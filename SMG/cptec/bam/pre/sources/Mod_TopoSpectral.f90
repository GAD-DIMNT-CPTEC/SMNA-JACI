!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_TopoSpectral </br></br>
!#
!# **Brief**: Module responsible for generating the spectral topography and
!# recomposed topography to the Gaussian grid of the model </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/dataout/Topography.GZZZZZ  (Ex.: pre/dataout/Topography.G00450) </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/Topography.TQXXXX (Ex.: pre/dataout/Topography.TQ0299) </br>
!# &bull; pre/dataout/TopographyRec.GZZZZZ.dat or (pre/dataout/TopographyRec.GZZZZZ) </br>
!# &bull; pre/dataout/TopographyRec.GZZZZZ.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita          - version: 1.1.1 </li>
!#  <li>03-06-2019 - Eduardo Khamis  - version: 2.0.0 </li>
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

module Mod_TopoSpectral

  use Mod_FastFourierTransform, only : createFFT, destroyFFT
  use Mod_LegendreTransform, only : createGaussRep, createSpectralRep, createLegTrans, gLats, destroyLegendreObjects
  use Mod_SpectralGrid, only : transp, grid2SpecCoef, specCoef2Grid
  use Mod_Get_xMax_yMax, only : getxMaxyMax
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData

  implicit none

  public :: getNameTopoSpectral, initTopoSpectral, generateTopoSpectral, shouldRunTopoSpectral

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !  input parameters

  type TopoSpectralNameListData
    logical :: linearGrid = .false.
    character (len = maxPathLength) :: varNameT = 'Topography' 
    !# prefix for input and output Gaussian Grid Topography Filename
  end type TopoSpectralNameListData

  type(varCommonNameListData)    :: varCommon
  type(TopoSpectralNameListData) :: var
  namelist /TopoSpectralNameList/   var

  ! internal variables

  integer :: xmx, yMaxHf, &
    mEnd1, mEnd2, mnwv0, mnwv1, mnwv2, mnwv3
  character (len = 7) :: trunc = '.T     '
  character (len = 7) :: nLats = '.G     '

  integer :: nftpi = 10   ! To read Reagular Grid Topography
  integer :: nftpo = 20   ! To write Unformatted Topography Sprectral Coefficients
  integer :: nfout = 30   ! To write GrADS Topography data on a Gaussian Grid
  integer :: nfctl = 40   ! To write Output data Description

  integer :: lRec

  real (kind = p_r4), dimension (:, :), allocatable :: topogIn, topogOut
  real (kind = p_r4), dimension (:), allocatable :: coefTopOut
  real (kind = p_r8), dimension (:), allocatable :: coefTop
  real (kind = p_r8), dimension (:, :), allocatable :: topog


contains


  function getNameTopoSpectral() result(returnModuleName)
    !# Returns TopoSpectral Module Name
    !# ---
    !# @info
    !# **Brief:** Returns TopoSpectral Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "TopoSpectral"
  end function getNameTopoSpectral


  subroutine initTopoSpectral(nameListFileUnit, varCommon_)
    !# Initialization of TopoSpectral module
    !# ---
    !# @info
    !# **Brief:** Initialization of TopoSpectral module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
  
    !      integer :: ios
    !
    !      mEnd=62            ! Spectral Resolution Horizontal Truncation
    !      linearGrid=.false. ! Flag to Set Linear (T) or Quadratic Gaussian Grid (T)
    !      grads=.true.       ! Flag to Get Recomposed Topography
    !      dirMain='./ '      ! Main Datain/Dataout Directory
    !
    !      open (unit=nfinp, FILE='./'//nameNML, &
    !            FORM='FORMATTED', ACCESS='SEQUENTIAL', &
    !            ACTION='READ', status='OLD', IOSTAT=ios)
    !      if (ios /= 0) then
    !         write (unit=nferr, fmt='(3A,I4)') &
    !               ' ** (Error) ** open file ', &
    !                 './'//nameNML, &
    !               ' returned iostat = ', ios
    !         stop  ' ** (Error) **'
    !      end if
    !      read  (unit=nfinp, NML=topoSpectralNameList)
    !      close (unit=nfinp)

    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_
        
    read(unit = nameListFileUnit, nml = TopoSpectralNameList)
    varCommon = varCommon_
 
    call getxMaxyMax (varCommon%mEnd, varCommon%xMax, varCommon%yMax)

    if (var%linearGrid) then
      trunc(3:3) = 'L'
    else
      trunc(3:3) = 'Q'
    end if
    write (trunc(4:7), fmt = '(I4.4)') varCommon%mEnd
    write (nLats(3:7), '(I5.5)') varCommon%yMax

    mEnd1 = varCommon%mEnd + 1
    mEnd2 = varCommon%mEnd + 2
    mnwv2 = mEnd1 * mEnd2
    mnwv0 = mnwv2 / 2
    mnwv3 = mnwv2 + 2 * mEnd1
    mnwv1 = mnwv3 / 2

    xmx = varCommon%xMax + 2
    yMaxHf = varCommon%yMax / 2

  end subroutine initTopoSpectral


  function shouldRunTopoSpectral() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getGaussianGridTopoOutFileName()) .or. .not. fileExists(getSpectralCoeficientTopoOutFileName())
  end function shouldRunTopoSpectral


  function getSpectralCoeficientTopoOutFileName() result(topoSpectralOutFilename)
    !# Gets Spectral Coefficient of Topography Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets Spectral Coefficient of Topography Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: topoSpectralOutFilename

    topoSpectralOutFilename = trim(varCommon%dirPreOut) // trim(var%varNameT) // trunc
  end function getSpectralCoeficientTopoOutFileName


  function getGaussianGridTopoOutFileName() result(gaussianGridTopoOutFilename)
    !# Gets Gaussian Grid Topography Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets Gaussian Grid Topography Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: gaussianGridTopoOutFilename
    
    gaussianGridTopoOutFilename = trim(varCommon%dirModelIn) // trim(var%varNameT) // 'Rec' // nLats
  end function getGaussianGridTopoOutFileName


  subroutine printNameList()
    !# Prints NameList
    !# ---
    !# @info
    !# **Brief:** Prints NameList. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    write (unit = p_nfprt, fmt = '(/,A)')    ' &TopoSpectralNameList'
    write (unit = p_nfprt, fmt = '(A,I6)')   '        mEnd = ', varCommon%mEnd
    write (unit = p_nfprt, fmt = '(A,L6)')   '  linearGrid = ', var%linearGrid
    write (unit = p_nfprt, fmt = '(A,L6)')   '       grads = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')      '     dirPreOut  = ' // trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A)')      '     dirModelIn = ' // trim(varCommon%dirModelIn)
    write (unit = p_nfprt, fmt = '(A,/)')    ' /'

  end subroutine printNameList


  function generateTopoSpectral() result(isExecOk)
    !# Generates TopoSpectral
    !# ---
    !# @info
    !# **Brief:** Generates TopoSpectral. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk

    isExecOk = .false.

    call createSpectralRep (mEnd1, mEnd2, mnwv1)
    call createGaussRep (varCommon%yMax, yMaxHf)
    call createFFT (varCommon%xMax)
    call createLegTrans (mnwv0, mnwv1, mnwv2, mnwv3, mEnd1, mEnd2, yMaxHf)

    allocate (topog(varCommon%xMax, varCommon%yMax), topogIn(varCommon%xMax, varCommon%yMax))
    allocate (coefTop(mnwv2), coefTopOut(mnwv2))
    if (varCommon%grads) allocate (topogOut(varCommon%xMax, varCommon%yMax))

    ! Read in Gaussian Grid Topography
    inquire (iolength = lRec) topogIn
    nftpi = openFile(trim(varCommon%dirPreOut) // trim(var%varNameT) // nLats, 'unformatted', 'direct', lRec, 'read', 'old')
    if(nftpi < 0) return
    read  (unit = nftpi, rec = 1) topogIn
    close (unit = nftpi)

    topog = real(topogIn, p_r8)
    ! SpectraL Coefficient of Topography
    call grid2SpecCoef (mnwv2, mEnd1, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, topog, coefTop)

    ! Write Out Topography Spectral Coefficient
    nftpo = openFile(getSpectralCoeficientTopoOutFileName(), 'unformatted', 'sequential', -1, 'write', 'replace')
    if(nftpo < 0) return
    coefTopOut = real(coefTop, p_r4)
    write (unit = nftpo) coefTopOut
    close (unit = nftpo)

    ! Write In Gaussian Grid Topography
    nfout = openFile(getGaussianGridTopoOutFileName(), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(nfout < 0) return
    write (unit = nfout, rec = 1) topogIn
    close(nfout, status = 'keep')
    if (varCommon%grads) then
      ! SpectraL Coefficient of Topography
      !         call transp(mnwv2, mEnd1, mEnd2, coefTop) ! need this?
      call specCoef2Grid (mnwv2, mnwv3, mEnd1, mEnd2, varCommon%xMax, varCommon%yMax, xmx, yMaxHf, coefTop, topog)
      topogOut = real(topog, p_r4)

      ! Write Out Recomposed Gaussian Grid Topography
      nfout = openFile(trim(varCommon%dirPreOut) // trim(var%varNameT) // 'Rec' // nLats, &
        'unformatted', 'direct', lRec, 'write', 'replace')
      if(nfout < 0) return
      write (unit = nfout, rec = 1) topogIn
      write (unit = nfout, rec = 2) topogOut
      close (unit = nfout)

      ! Write GrADS Control File
      nfctl = openFile(trim(varCommon%dirPreOut) // trim(var%varNameT) // 'Rec' // nLats // '.ctl', &
      'formatted', 'sequential', -1, 'write', 'replace')
      if (nfctl < 0) return

      write (unit = nfctl, fmt = '(A)') 'DSET ^' // &
        trim(var%varNameT) // 'Rec' // nLats
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', p_undef
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A)') 'TITLE Topography on a Gaussian Grid'
      write (unit = nfctl, fmt = '(A)') '*'
      write (unit = nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
        'XDEF ', varCommon%xMax, ' LINEAR ', 0.0_p_r8, 360.0_p_r8 / real(varCommon%xMax, p_r8)
      write (unit = nfctl, fmt = '(A,I5,A)') 'YDEF ', varCommon%yMax, ' LEVELS '
      write (unit = nfctl, fmt = '(8F10.5)') gLats(varCommon%yMax:1:-1)
      write (unit = nfctl, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
      write (unit = nfctl, fmt = '(A)') 'TDEF 1 LINEAR JAN2005 1MO'
      write (unit = nfctl, fmt = '(A)') 'VARS 2'
      write (unit = nfctl, fmt = '(A)') 'TOPI 0 99 Interpolated Topography [m]'
      write (unit = nfctl, fmt = '(A)') 'TOPR 0 99 Recomposed Topography [m]'
      write (unit = nfctl, fmt = '(A)') 'ENDVARS'
      close (unit = nfctl)

    end if

    call destroyFFT()
    call destroyLegendreObjects()
    
    isExecOk = .true.
  end function generateTopoSpectral

end module Mod_TopoSpectral
