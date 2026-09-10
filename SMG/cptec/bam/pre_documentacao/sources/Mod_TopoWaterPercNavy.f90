!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_TopoWaterPercNavy </br></br>
!#
!# **Brief**: Module responsible for extracting the data of topograph and
!# porcent of water from originary data US-Navy 10’ (file navytm.form) and
!# generate the files from read in the GrADSS </br></br>
!# 
!# The topografy and water/continent porcentage informations are derived from
!# two data sets of land elevation obtained from U.S. Navy – Fleet Numerical 
!# Oceanography Center (FNOC) and do U.S. Geological Survey’s – Earth Resources
!# Observation and Science (EROS) Data Center. The Navy-FNOC data set began to 
!# be created in the mid 60's. After, these data were reprocessed by NCAR and 
!# NOAA-NGDC. For each area of 10x10 minutes of degrees, the US-Navy 10' data
!# set includes mean land height, maximum and minimum heights, number and 
!# orientation of significant ridges, and percentage of water and urban areas. 
!# In the pre-processing, topography data are normalized to the values between
!# 100 to 320 and water percentage values between 0 to 100. Corrections are 
!# also applied to avoid negative topography and undesirable surface temperature
!# values on areas covered with water on the continent. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/databcs/navytm.form </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/TopoNavy.dat </br>
!# &bull; pre/dataout/WaterNavy.dat </br>
!# &bull; pre/dataout/TopoWaterNavy.dat </br>
!# &bull; pre/dataout/TopoWaterNavy.ctl
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 1.2.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>18-09-2007 - Jose P. Bonatti - version: 1.0.0 </li>
!#  <li>04-06-2019 - Eduardo Khamis  - version: 1.2.0 </li>
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

module Mod_TopoWaterPercNavy
   ! First Point of Initial Data is at South Pole and Greenwitch
   ! First Point of Output  Data is at North Pole and Greenwitch
   use Mod_FileManager
   use Mod_Namelist, only : varCommonNameListData

   implicit none

   public :: generateTopoWaterPercNavy
   public :: getNameTopoWaterPercNavy
   public :: initTopoWaterPercNavy 
   public :: shouldRunTopoWaterPercNavy

   private
   include 'files.h'
   include 'pre.h'
   include 'precision.h'
   include 'messages.h'

   integer, parameter :: imBox=30, jmBox=30, nLonBox=72, nLatBox=36
   integer, parameter :: xMax=nLonBox*imBox, yMax=nLatBox*jmBox

   real (kind = p_r4) :: undef=-9999.0_p_r4

   integer :: nj, ni, io, im, jo, jm, j, lRec, ios

   real (kind = p_r4) :: dxy

   real (kind = p_r4), dimension (imBox,jmBox) :: tp, wp
   real (kind = p_r4), dimension (xMax,yMax) :: topog, water

   
   integer :: nfclm   
   !# To read Navy Topography data
   integer :: nftop   
   !# To write Navy Topography data
   integer :: nfwat   
   !# To write Water Percentage data
   integer :: nftpw   
   !# To write GrADS Navy Topography and Water Percentage data
   integer :: nfctl   
   !# To write Output data Description

   type TopoWaterPercNavyNameListData
     character (len = maxPathLength) :: varNameT   = 'TopoNavy'      
     !# file name prefix for writing Navy Topography data
     character (len = maxPathLength) :: varNameW   = 'WaterNavy'     
     !# file name prefix for writing Water Percentage data
     character (len = maxPathLength) :: varNameBCs = 'navytm.form'   
     !# file name prefix for reading Navy Topography data
     character (len = maxPathLength) :: varNameG   = 'TopoWaterNavy' 
     !# file name prefix for writing GrADS Navy Topography and Water Percentage data
   end type TopoWaterPercNavyNameListData

   type(varCommonNameListData)         :: varCommon
   type(TopoWaterPercNavyNameListData) :: var
   namelist /topoWaterPercNavyNameList/   var


contains


  function getNameTopoWaterPercNavy() result(returnModuleName)
    !# Returns TopoWaterPercNavy Module Name
    !# ---
    !# @info
    !# **Brief:** Returns TopoWaterPercNavy Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "TopoWaterPercNavy"
  end function getNameTopoWaterPercNavy
 
 
  subroutine initTopoWaterPercNavy(nameListFileUnit, varCommon_)
    !# Initalization of TopoWaterPercNavy module
    !# ---
    !# @info
    !# **Brief:** Initialization of TopoWaterPercNavy module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = TopoWaterPercNavyNameList)
    varCommon = varCommon_
  end subroutine initTopoWaterPercNavy
 
 
  function shouldRunTopoWaterPercNavy() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    logical :: shouldRun

    shouldRun = .not. ( fileExists(getTopoWaterPercNavyFileName()) .or. fileExists(getWaterNavyFileName()) )
  end function shouldRunTopoWaterPercNavy
 
 
  function getTopoWaterPercNavyFileName() result(topoWaterPercNavyOutFilename)
    !# Gets TopoWaterPercNavy Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets TopoWaterPercNavy Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    character(len = maxPathLength) :: topoWaterPercNavyOutFilename

    topoWaterPercNavyOutFilename = trim(varCommon%dirPreOut) // trim(var%varNameT) // '.dat'
  end function getTopoWaterPercNavyFileName


  function getWaterNavyFileName() result(waterNavyFileName)
    !# Gets Water Navy out Filename
    !# ---
    !# @info
    !# **Brief:** Gets Water Navy out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    character(len = maxPathLength) :: waterNavyFileName

    waterNavyFileName = trim(varCommon%dirPreOut) // trim(var%varNameW) // '.dat'
  end function getWaterNavyFileName   


  function generateTopoWaterPercNavy() result(isExecOk)
    !# Generates TopoWaterPercNavy
    !# ---
    !# @info
    !# **Brief:** Generates TopoWaterPercNavy. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer :: ios
    logical :: isExecOk

    isExecOk = .false.
   
    write (unit = p_nfprt, fmt = '(/,A)')  ' &InputDim'
    write (unit = p_nfprt, fmt = '(A,L6)') '    GrADS = ', varCommon%grads
    write (unit = p_nfprt, fmt = '(A)')    '  DirPreBCs = '//trim(varCommon%dirBCs)
    write (unit = p_nfprt, fmt = '(A)')    '  DirPreOut = '//trim(varCommon%dirPreOut)
    write (unit = p_nfprt, fmt = '(A,/)')  ' /'
   

    nfclm = openFile(trim(varCommon%dirBCs)//trim(var%varNameBCs), 'formatted', 'sequential', -1, 'read', 'old')
    if (nfclm < 0) return
    do nj=1,nLatBox
       jo=1+jmBox*(nj-1)
       jm=jo+jmBox-1
       do ni=1,nLonBox
          io=1+imBox*(ni-1)
          im=io+imBox-1
          call getTopographyBox (imBox, jmBox, tp, wp)
          topog(io:im,jo:jm)=tp(1:imBox,1:jmBox)
          water(io:im,jo:jm)=wp(1:imBox,1:jmBox)
       end do
    end do
    close (unit = nfclm)
   
    call flipMatrix (xMax,yMax,topog)
    call flipMatrix (xMax,yMax,water)
   

    inquire (iolength=lRec) topog(:,1)
    nftop = openFile(getTopoWaterPercNavyFileName(), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(nftop < 0) return
    do j=1,yMax
       write (unit = nftop, rec = j) topog(:,j)
    end do
    close (unit = nftop)
   
    nfwat = openFile(getWaterNavyFileName(), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(nfwat < 0) return
    do j=1,yMax
       write (unit = nfwat, rec = j) water(:,j)
    end do
    close (unit = nfwat)
   
    if (varCommon%grads) then
       call generateGrads()
    end if
   
    isExecOk = .true.
  end function generateTopoWaterPercNavy


  subroutine generateGrads()
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer :: ios
    character(len=maxPathLength) gradsBaseName
    character(len=maxPathLength) gradsPathBaseName

    gradsBaseName = var%varNameG
    gradsPathBaseName = trim(varCommon%dirPreOut)//trim(gradsBaseName)


    nftpw = openFile(trim(gradsPathBaseName)//'.dat', 'unformatted', 'direct', lRec, 'write', 'replace')
    if(nftpw < 0) return
    do j=1,yMax
      write (unit = nftpw, rec = j) topog(:,j)
    end do
    do j=1,yMax
      write (unit = nftpw, rec = j+yMax) water(:,j)
    end do
    close (unit = nftpw)

    dxy=360.0_p_r4/real(xMax,p_r4)
    nfctl = openFile(trim(gradsPathBaseName)//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
    if(nfctl < 0) return

    write (unit = nfctl, fmt = '(A)') 'DSET ^'// trim(gradsBaseName)//'.dat'
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A,1PG12.5)') 'UNDEF ', undef
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A)') 'TITLE Navy Topography and Water Percentage'
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A,I5,A,F12.7,F10.7)') 'XDEF ', xMax, ' LINEAR ', &
                                                    0.5_p_r4*dxy, dxy
    write (unit = nfctl, fmt = '(A,I5,A,F12.7,F10.7)') 'YDEF ', yMax, ' LINEAR ', &
                                                    -90.0_p_r4+0.5_p_r4*dxy, dxy
    write (unit = nfctl, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
    write (unit = nfctl, fmt = '(A)') 'TDEF 1 LINEAR JAN2005 1MO'
    write (unit = nfctl, fmt = '(A)') '*'
    write (unit = nfctl, fmt = '(A)') 'VARS 2'
    write (unit = nfctl, fmt = '(A)') 'TOPO 0 99 Topography [m]'
    write (unit = nfctl, fmt = '(A)') 'WPER 0 99 Percentage of Water [%]'
    write (unit = nfctl, fmt = '(A)') 'ENDVARS'
    close (unit = nfctl)
  end subroutine generateGrads

  
  subroutine getTopographyBox (imBox, jmBox, tp, wp)
    !# Gets TopographyBox
    !# ---
    !# @info
    !# **Brief:** Gets TopographyBox. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    
    ! Reads Modified Navy Data
    ! Input Data Contains Only Terrain Height (Index 1)
    ! and Percentage of Water (Index 2)
   
    ! Input:
    ! ImBox - Number of Longitude Points in a Box
    ! JmBox - Number of Latitude Points in a Box
   
    ! Output:
    ! Tp(ImBox,JmBox) - Terrain Height
    ! Wp(ImBox,JmBox) - Percentage of Water
   
    implicit none
   
    integer, intent(in) :: imBox
    integer, intent(in) :: jmBox
   
    real (kind = p_r4), intent(OUT) :: tp(imBox,jmBox)
    real (kind = p_r4), intent(OUT) :: wp(imBox,jmBox)
   
    integer :: lon
    integer :: lat
    integer :: i, j
   
    real (kind = p_r4) :: zfct
  
    integer :: itop(2,imBox,jmBox)
   
    logical, parameter :: IWant=.true., IWant2=.true.
   
    zfct=100.0_p_r4*(1200.0_p_r4/3937.0_p_r4)
   
    ! Lat - Northern Most Colatitude of Grid Box (5 Dg X 5 Dg)
    ! Lon - Western Most Longitude
    read (unit = nfclm, fmt = '(2I3)') lat, lon
    read (unit = nfclm, fmt = '(30I3)') itop
   
    lat=90-lat
   
    do j=1,jmBox
      do i=1,imBox
   
        ! Corrections (By J.P.Bonatti, 18 Nov 1999)
   
        ! itop(1,i,j) - Is the Normalized Topography:
        !               Must Be Between 100 and 320
        ! itop(2,i,j) - Is The Percentage of Water:
        !               Must Be Between 0 and 100
   
        ! Correction Were Done Analysing The Surrounding Values
   
        if (IWant) then
   
          ! Wrong Value
          if (itop(1,i,j) == 330) then
            itop(1,i,j)=300
          end if
   
          ! Wrong Value
          if (itop(1,i,j) == 340) then
            itop(1,i,j)=300
          end if
   
          ! Wrong Value
          if (itop(1,i,j) == 357) then
            itop(1,i,j)=257
          end if
   
          ! Wrong Value
          if (itop(1,i,j) == 433) then
            itop(1,i,j)=113
          end if
   
          ! May Be Missing Value
          if (itop(1,i,j) == 511) then
            itop(1,i,j)=100
          end if
   
          ! May Be Missing Value
          if (itop(2,i,j) == 127) then
            if (itop(1,i,j) == 100) then
              itop(2,i,j)=100
            else
              itop(2,i,j)=1
              end if
          end if
   
          if (IWant2) then
   
            ! To Avoid Negative Topography
            if (itop(1,i,j) < 100) then
              itop(1,i,j)=100
            end if
   
            ! To Avoid Non-Desirable Values of Surface
            ! Temperature Over Water Inside Continents
            if (itop(1,i,j) > 100 .and. itop(2,i,j) >= 50) then
              itop(2,i,j)=49
            end if
   
          end if
   
        end if
   
        tp(i,j)=zfct*real(itop(1,i,j)-100,p_r4)
        wp(i,j)=real(itop(2,i,j),p_r4)
   
      end do
    end do
   
  end subroutine getTopographyBox
   

  subroutine flipMatrix (xMax, yMax, h)
    !# Flips Matrix
    !# ---
    !# @info
    !# **Brief:** Flips Matrix. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    
    ! Flips Over The Southern and Northern Hemispheres
   
    ! Input:
    !   Imax       - Column Dimension of h(Imax,Jmax)
    !   Jmax       - Row Dimension of h(Imax,Jmax)
    ! h(Imax,Jmax) - Matrix to be Flipped
   
    ! Output:
    ! h(Imax,Jmax) - Flipped Matrix
   
    implicit none
   
    integer, intent(in) :: xMax
    integer, intent(in) :: yMax
    real (kind = p_r4), intent(inout) :: h(xMax,yMax)
   
    integer :: yMaxd, yMaxd1
   
    real (kind = p_r4) :: wk(xMax,yMax)
   
    yMaxd=yMax/2
    yMaxd1=yMaxd+1
   
    wk=h
    h(:,1:yMaxd)=wk(:,yMax:yMaxd1:-1)
    h(:,yMaxd1:yMax)=wk(:,yMaxd:1:-1)
   
  end subroutine flipMatrix
   

end module Mod_TopoWaterPercNavy
