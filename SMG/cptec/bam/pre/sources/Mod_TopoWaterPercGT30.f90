!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_TopoWaterPercGT30 </br></br>
!#
!# **Brief**: Module responsible for extracting the data of topograph and porcent
!# of water from originary data Gtopo30, which are subdivided in 33 files
!# correspondents to regions (or latitude zones) of the globe (files of the type
!# ElonNlat.dat, WlonNlat.dat, ElonSlat.dat and WlonSlat.dat, on what lon and lat
!# are the longitude and the latitude of the data) </br></br>
!# 
!# GTOPO30 is a land elevation data set with global cover (90°S-90°N, 180ºW-180°E)
!# and horizontal grid of 30 seconds of degree (0.008333333 degrees), resulting in
!# a Digital Elevation Model (DEM) with dimensions between of 21600 rows e 43200
!# columns. The data represent the elevation in meters, above mean sea level and
!# the values range from 407 to 8752 meters. Ocean areas are masked as "no data"
!# and encoded with a value -9999 (undef). Smalls islands in the oceans with areas
!# less than 1 km2 are not represented by DEM.
!# The ocean areas, coded as undef, assume values of 100% for water percentage
!# and 0 for topografy. On the continent, water percentage is defined as 0% and
!# the topografy varies from a minimum value of one meter (1m).
!# For simplification, the original data has been converted to the intermediate
!# resolution of 0.16º. The data set derived for this intermediate resolution is
!# composed of the variables average elevation above sea level and water/continent
!# percentage. </br></br>
!# 
!#
!# **Files in:**
!#
!# &bull; pre/databcs/gtopo30/W180N90.dat </br>
!# &bull; pre/databcs/gtopo30/W140N90.dat </br>
!# &bull; pre/databcs/gtopo30/W100N90.dat </br>
!# &bull; pre/databcs/gtopo30/W060N90.dat </br>
!# &bull; pre/databcs/gtopo30/W020N90.dat </br>
!# &bull; pre/databcs/gtopo30/E020N90.dat </br>
!# &bull; pre/databcs/gtopo30/E060N90.dat </br>
!# &bull; pre/databcs/gtopo30/E100N90.dat </br>
!# &bull; pre/databcs/gtopo30/E140N90.dat </br>
!# &bull; pre/databcs/gtopo30/W180N40.dat </br>
!# &bull; pre/databcs/gtopo30/W140N40.dat </br>
!# &bull; pre/databcs/gtopo30/W100N40.dat </br>
!# &bull; pre/databcs/gtopo30/W060N40.dat </br>
!# &bull; pre/databcs/gtopo30/W020N40.dat </br>
!# &bull; pre/databcs/gtopo30/E020N40.dat </br> 
!# &bull; pre/databcs/gtopo30/E060N40.dat </br>
!# &bull; pre/databcs/gtopo30/E100N40.dat </br>
!# &bull; pre/databcs/gtopo30/E140N40.dat </br>
!# &bull; pre/databcs/gtopo30/W180S10.dat </br>
!# &bull; pre/databcs/gtopo30/W140S10.dat </br>
!# &bull; pre/databcs/gtopo30/W100S10.dat </br>
!# &bull; pre/databcs/gtopo30/W060S10.dat </br>
!# &bull; pre/databcs/gtopo30/W020S10.dat </br>
!# &bull; pre/databcs/gtopo30/E020S10.dat </br>
!# &bull; pre/databcs/gtopo30/E060S10.dat </br>
!# &bull; pre/databcs/gtopo30/E100S10.dat </br>
!# &bull; pre/databcs/gtopo30/E140S10.dat </br>
!# &bull; pre/databcs/gtopo30/W180S60.dat </br>
!# &bull; pre/databcs/gtopo30/W120S60.dat </br>
!# &bull; pre/databcs/gtopo30/W060S60.dat </br>
!# &bull; pre/databcs/gtopo30/W000S60.dat </br>
!# &bull; pre/databcs/gtopo30/E060S60.dat </br>
!# &bull; pre/databcs/gtopo30/E120S60.dat
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; pre/dataout/TopoGT30.dat </br>
!# &bull; pre/dataout/WaterGT30.dat </br>
!# &bull; pre/dataout/TopoWaterGT30"gd".dat </br>
!# &bull; pre/dataout/TopoWaterGT30"gd".ctl
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
!#  <li>04-06-2019 - Eduardo Khamis  - version: 2.0.0 </li>
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

module Mod_TopoWaterPercGT30
   ! First Point of Initial Data is at North Pole and I. D. Line
   ! First Point of Output  Data is at North Pole and Greenwhich
   use Mod_FileManager
   use Mod_Namelist, only : varCommonNameListData

   implicit none

   public :: generateTopoWaterPercGT30
   public :: getNameTopoWaterPercGT30
   public :: initTopoWaterPercGT30
   public :: shouldRunTopoWaterPercGT30

   private
   include 'files.h'
   include 'pre.h'
   include 'precision.h'
   include 'messages.h'

   integer, parameter :: numBox=33, xMax=43200, yMax=21600

   real (kind = p_r4), parameter :: tpMinSea=0.0_p_r4, tpMinLand=1.0_p_r4, &
                                wpMin=0.0_p_r4, wpMax=100_p_r4, &
                                undef=-9999.0_p_r4

   integer :: n, io, im, jo, jm, j, m, mx, mr, ms, ma, mb, &
              ix, jx, i1, i2, j1, j2, nr, lRecOut

   real (kind = p_r4) :: dxy, lono, lato, long, latg

   character (len = 3) :: gd

   real (kind = p_r4), dimension (:,:), allocatable :: tp, wp
   real (kind = p_r4), dimension (:,:), allocatable :: topog, water
   
   integer :: nferr=0    ! Standard Error print Out
   integer :: nfclm=10   ! To read GT30 Topography data
   integer :: nftop=20   ! To write GT30 Topography data
   integer :: nfwat=30   ! To write Water Percentage data
   integer :: nftpw=40   ! To write GrADS GT30 Topography and Water Percentage data
   integer :: nfctl=50   ! To write Output data Description

   character (len = 7), dimension(numBox) :: topoName = (/&
             'W180N90', 'W140N90', 'W100N90', 'W060N90', 'W020N90', &
             'E020N90', 'E060N90', 'E100N90', 'E140N90', 'W180N40', &
             'W140N40', 'W100N40', 'W060N40', 'W020N40', 'E020N40', &
             'E060N40', 'E100N40', 'E140N40', 'W180S10', 'W140S10', &
             'W100S10', 'W060S10', 'W020S10', 'E020S10', 'E060S10', &
             'E100S10', 'E140S10', 'W180S60', 'W120S60', 'W060S60', &
             'W000S60', 'E060S60', 'E120S60' /)

   integer, dimension(numBox) :: imBox = (/ &
            4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, &
            4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, &
            4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, 4800, &
            7200, 7200, 7200, 7200, 7200, 7200 /)

   integer, dimension(numBox) :: jmBox = (/ &
            6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, &
            6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, &
            6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, &
            3600, 3600, 3600, 3600, 3600, 3600 /)

   type TopoWaterPercGT30NameListData
     character (len = maxPathLength) :: varNameT = 'TopoGT30'      
     !# TopoGT30 prefix out filename 
     character (len = maxPathLength) :: varNameW = 'WaterGT30'     
     !# WaterGT30 prefix out filename 
     character (len = maxPathLength) :: varNameG = 'TopoWaterGT30' 
     !# grADS prefix out filename  
     character (len = maxPathLength) :: dirBCs = './'              
     !# local Boundary Conditions directory
   end type TopoWaterPercGT30NameListData

   type(varCommonNameListData)         :: varCommon
   type(TopoWaterPercGT30NameListData) :: var
   namelist /TopoWaterPercGT30NameList/   var


contains


  function getNameTopoWaterPercGT30() result(returnModuleName)
    !# Returns TopoWaterPercGT30 Module Name
    !# ---
    !# @info
    !# **Brief:** Returns TopoWaterPercGT30 Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "TopoWaterPercGT30"
  end function getNameTopoWaterPercGT30
 
 
  subroutine initTopoWaterPercGT30(nameListFileUnit, varCommon_)
    !# Initialization of TopoWaterPercGT30 module
    !# ---
    !# @info
    !# **Brief:** Initialization of TopoWaterPercGT30 module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = TopoWaterPercGT30NameList)
    varCommon = varCommon_
  end subroutine initTopoWaterPercGT30
 
 
  function shouldRunTopoWaterPercGT30() result(shouldRun)
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

    shouldRun = .not. ( fileExists(getTopoGT30OutFileName()) .or. fileExists(getWaterGT30OutFileName()) )
  end function shouldRunTopoWaterPercGT30
 
 
  function getTopoGT30OutFileName() result(topoGT30OutFilename)
    !# Gets TopoGT30 Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets TopoGT30 Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    character(len = maxPathLength) :: topoGT30OutFilename

    topoGT30OutFilename = trim(varCommon%dirPreOut) // trim(var%varNameT) // '.dat'
  end function getTopoGT30OutFileName
 

  function getWaterGT30OutFileName() result(waterGT30OutFilename)
    !# Gets WaterGT30 Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets WaterGT30 Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    character(len = maxPathLength) :: waterGT30OutFilename

    waterGT30OutFilename = trim(varCommon%dirPreOut) // trim(var%varNameW) // '.dat'
  end function getWaterGT30OutFileName



  function generateTopoWaterPercGT30() result(isExecOk)
    !# Gets WaterGT30 Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets WaterGT30 Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    implicit none
    integer                         :: ios
    logical :: isExecOk

    isExecOk = .false.
   
    write (unit=p_nfprt, fmt='(/,A)')  ' &InputDim'
    write (unit=p_nfprt, fmt='(A,L6)') '      GrADS = ', varCommon%grads
    write (unit=p_nfprt, fmt='(A)')    '  DirPreBCs = '//trim(var%dirBCs)
    write (unit=p_nfprt, fmt='(A)')    '  DirPreOut = '//trim(varCommon%dirPreOut)
    write (unit=p_nfprt, fmt='(A,/)')  ' /'

    allocate(topog(xMax,yMax))
    allocate(water(xMax,yMax))

    do n=1,numBox
      allocate (tp(imBox(n),jmBox(n)), wp(imBox(n),jmBox(n)))

      nfclm = openFile(trim(var%dirBCs)//topoName(n)//'.dat', 'unformatted', 'sequential', -1, 'read', 'old')
      if(nfclm < 0) return
      read  (unit=nfclm) tp
      close (unit=nfclm)
      where (tp == undef)
        ! Over Sea: Wp=100% and Tp=0m
        wp=wpMax
        tp=tpMinSea
      elsewhere
        ! Over Land: Wp=0% and Tp=MAX(1m,Tp)
        wp=wpMin
        tp=max(tpMinLand,tp)
      endwhere
      select case (n)
        case (1:9)
          io=1+(n-1)*imBox(1)
          im=n*imBox(1)
          jo=1
          jm=jmBox(1)
        case (10:18)
          io=1+(n-10)*imBox(10)
          im=(n-9)*imBox(10)
          jo=1+jmBox(1)
          jm=jmBox(1)+jmBox(10)
        case (19:27)
          io=1+(n-19)*imBox(19)
          im=(n-18)*imBox(19)
          jo=1+jmBox(1)+jmBox(10)
          jm=jmBox(1)+jmBox(10)+jmBox(19)
        case (28:numBox)
          io=1+(n-28)*imBox(28)
          im=(n-27)*imBox(28)
          jo=1+jmBox(1)+jmBox(10)+jmBox(19)
          jm=jmBox(1)+jmBox(10)+jmBox(19)+jmBox(28)
      end select
      topog(io:im,jo:jm)=tp(1:imBox(n),1:jmBox(n))
      water(io:im,jo:jm)=wp(1:imBox(n),1:jmBox(n))
      write (unit=p_nfprt, fmt='(7I8)') n, io, im, jo, jm, imBox(n), jmBox(n)
      write (unit=p_nfprt, fmt='(1P2G12.5)') minval(tp), maxval(tp)
      write (unit=p_nfprt, fmt='(1P2G12.5)') minval(wp), maxval(wp)
      deallocate (tp, wp)
    end do
   
    call flipMatrix (topog,xMax,yMax)
    call flipMatrix (water,xMax,yMax)
   
    inquire (iolength=lRecOut) topog(:,1)
    nftop = openFile(getTopoGT30OutFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
    if(nftop < 0) return
      do j=1,yMax
        write (unit=nftop, rec=j) topog(:,j)
      end do
      close (unit=nftop)
   
      nfwat = openFile(getWaterGT30OutFileName(), 'unformatted', 'direct', lRecOut, 'write', 'replace')
      if(nfwat < 0) return
      do j=1,yMax
        write (unit=nfwat, rec=j) water(:,j)
      end do
      close (unit=nfwat)
   
      if (varCommon%grads) then
        call generateGrads() 
    end if

    deallocate(topog)
    deallocate(water)
   
    isExecOk = .true.
  end function generateTopoWaterPercGT30
  

  subroutine generateGrads()
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Jose P. Bonatti </br>
    !# **Date**: nov/2004 </br>
    !# @endinfo
    implicit none
    integer :: ios
    integer :: lRecGrads
    character(len=maxPathLength) gradsBaseName
    character(len=maxPathLength) gradsPathBaseName
  
    mx=128
    mr=16
    ms=mx/mr
    ix=xMax/mr
    jx=yMax/ms
    dxy=360.0_p_r4/real(xMax,p_r4)
    lono=0.5_p_r4*dxy
    lato=90.0_p_r4-0.5_p_r4*dxy
    m=0
    do ma=1,mr
      i1=1+(ma-1)*ix
      i2=ma*ix
      long=lono+real(i1-1,p_r4)*dxy
      do mb=1,ms
        m=m+1
        j1=1+(mb-1)*jx
        j2=mb*jx
        latg=lato-real(j2-1,p_r4)*dxy
        write (gd, fmt='(I3.3)') m

        gradsBaseName = trim(var%varNameG)//gd
        gradsPathBaseName = trim(varCommon%dirPreOut)//trim(gradsBaseName)
        inquire (iolength=lRecGrads) topog(i1:i2,1)
        nftpw = openFile(trim(gradsPathBaseName)//'.dat', 'unformatted', 'direct', lRecGrads, 'write', 'replace')
        if(nftpw < 0) return
        nr=0
        do j=j1,j2
          nr=nr+1
          write (unit=nftpw, rec=nr) topog(i1:i2,j)
        end do
        do j=j1,j2
          nr=nr+1
          write (unit=nftpw, rec=nr) water(i1:i2,j)
        end do
        close (unit=nftpw)
        write (unit=p_nfprt, fmt='(A,4I6,2F14.7)') &
          ' '//gd//' ', i1, i2, j1, j2, long, latg
          dxy=360.0_p_r4/real(xMax,p_r4)

        nfctl = openFile(trim(gradsPathBaseName)//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
        if(nfctl < 0) return
          write (unit=nfctl, fmt='(A)') 'DSET ^'// trim(gradsBaseName)//'.dat'
          write (unit=nfctl, fmt='(A)') '*'
          write (unit=nfctl, fmt='(A)') 'OPTIONS YREV BIG_ENDIAN'
          write (unit=nfctl, fmt='(A)') '*'
          write (unit=nfctl, fmt='(A,1PG12.5)') 'UNDEF ', undef
          write (unit=nfctl, fmt='(A)') '*'
          write (unit=nfctl, fmt='(A)') 'TITLE GT30 Topography and Water Percentage'
          write (unit=nfctl, fmt='(A)') '*'
          write (unit=nfctl, fmt='(A,I5,A,2F15.10)') 'XDEF ', ix, ' LINEAR ', long, dxy
          write (unit=nfctl, fmt='(A,I5,A,2F15.10)') 'YDEF ', jx, ' LINEAR ', latg, dxy
          write (unit=nfctl, fmt='(A)') 'ZDEF 1 LEVELS 1000'
          write (unit=nfctl, fmt='(A)') 'TDEF 1 LINEAR JAN2005 1MO'
          write (unit=nfctl, fmt='(A)') '*'
          write (unit=nfctl, fmt='(A)') 'VARS 2'
          write (unit=nfctl, fmt='(A)') 'TOPO 0 99 Topography [m]'
          write (unit=nfctl, fmt='(A)') 'WPER 0 99 Percentage of Water [%]'
          write (unit=nfctl, fmt='(A)') 'ENDVARS'
          close (unit=nfctl)
      end do
    end do
  end subroutine generateGrads


  subroutine flipMatrix (h, xMax, yMax)  
    !# Flips Matrix
    !# ---
    !# @info
    !# **Brief:** Flips Matrix. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jun/2019 </br>
    !# @endinfo
    ! Flips a Matrix Over I.D.L. and Greenwitch
  
    ! Input:
    ! h(Imax,Jmax) - Matrix to be Flipped
    !   Imax       - Column Dimension of h(Imax,Jmax)
    !   Jmax       - Row Dimension of h(Imax,Jmax)
   
    ! Output:
    ! h(Imax,Jmax) - Flipped Matrix
    implicit none

    integer, parameter :: p_r4 = selected_real_kind(6) ! kind for 32-bits real Numbers
    integer, intent(in) :: xMax
    integer, intent(in) :: yMax

    real (kind = p_r4), intent(inout) :: h(xMax,yMax)
   
    integer :: xMaxd, xMaxd1, j
   
    real (kind = p_r4) :: wk(xMax/2)
   
    xMaxd=xMax/2
    xMaxd1=xMaxd+1
   
    do j=1,yMax
      wk=h(xMaxd1:xMax,j)
      h(xMaxd1:xMax,j)=h(1:xMaxd,j)
      h(1:xMaxd,j)=wk
    end do
   
  end subroutine flipMatrix

   
end module Mod_TopoWaterPercGT30
