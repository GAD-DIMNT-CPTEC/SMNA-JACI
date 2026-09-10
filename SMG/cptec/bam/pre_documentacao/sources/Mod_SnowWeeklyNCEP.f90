!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_SnowWeeklyNCEP </br></br>
!#
!# **Brief**: Module responsible for generating the initial condition of water
!# storage</br></br>
!# 
!# This subroutine uses the accumulated snow data provided by National Centers
!# for Environmental Prediction (NCEP) as a second way to generate the initial 
!# condition of water storage due to precipitation interception. This data is a
!# estimate obtained from satellite information that has a frequency of one day.
!# Nowadays this file is not being available by the centre.</br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/datain/gdas1.ThhZ.snogrd.YYYYMMDDHH (Ex.: pre/datain/gdas1.T00Z.snogrd.2004032600)
!# </br></br>
!# 
!# **intermediate files:*
!# &bull; pre/dataout/SNOWWeekly.YYYYMMDD (Ex.: pre/datain/SNOWWeekly.2004032600)</br>
!# &bull; pre/dataout/SNOWWeekly.YYYYMMDD.ctl
!# </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/SNOWWeeklyYYYYMMDD.GZZZZZ (Ex.: model/datain/SNOWWeekly20150430.G00450) </br>
!# &bull; pre/dataout/SNOWWeeklyYYYYMMDD.GZZZZZ.ctl 
!# </br></br>
!# 
!# The SNOWWeeklyYYYYMMDD.GZZZZZ file that is in the model/datain directory is 
!# the file that the model will use in its execution, has the specific format 
!# (unformatted, direct access, 64-byte integer) and the SNOWWeeklyYYYYMMDD.GZZZZZ
!# file that is in the pre/dataout directory is specific to viewing in GrADS 
!# software.</br></br>
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
!#  <li>05-09-2019 - Eduardo Khamis    - version: 2.0.0 </li>
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

module Mod_SnowWeeklyNCEP
  ! First Point of Initial Data is near South Pole and near Greenwhich
  ! First Point of Output  Data is near North Pole and near Greenwhich
  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgOut

  implicit none

  public :: generateSnowWeeklyNCEP
  public :: getNameSnowWeeklyNCEP
  public :: initSnowWeeklyNCEP
  public :: shouldRunSnowWeeklyNCEP

  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  integer :: year, month, lRec

  character(len = 3), dimension(12) :: monthChar = &
            (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
               'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC' /)

  real (kind = p_r4), dimension (:,:), allocatable :: snow 

  integer :: nferr=0    
  !# Standard Error print Out
  integer :: nfinp=5    
  !# Standard read in
  integer :: nfprt=6    
  !# Standard print Out
  integer :: p_nfsno=10   
  !# To read Unformatted Weekly SNOW data
  integer :: p_nfout=20   
  !# To write Unformatted Weekly SNOW data
  integer :: p_nfctl=30   
  !# To write Output data Description

  type SnowWeeklyNCEPNameListData
    integer :: xDim
    integer :: yDim
    character(len = maxPathLength) :: varName='SNOWWeekly'         
    !# output prefix file name
    character(len = maxPathLength) :: fileSnow='gdas1.T  Z.snogrd' 
    !# input prefix file name
    character(len = maxPathLength) :: timeGrADS='  Z         '     
    !# grADS time mask for ctl file
    character(len = 10) :: dateSnowWeeklyNCEP
  end type SnowWeeklyNCEPNameListData

  type(varCommonNameListData)      :: varCommon
  type(SnowWeeklyNCEPNameListData) :: var
  namelist /SnowWeeklyNCEPNameList/   var

  character(len = *), parameter :: headerMsg = 'Snow Weekly NCEP          | '


  contains


  function getNameSnowWeeklyNCEP() result(returnModuleName)
    !# Returns SnowWeeklyNCEP Module Name
    !# ---
    !# @info
    !# **Brief:** Returns SnowWeeklyNCEP Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo  
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "SnowWeeklyNCEP"
  end function getNameSnowWeeklyNCEP


  subroutine initSnowWeeklyNCEP(nameListFileUnit, varCommon_)
    !# Initialization of SnowWeeklyNCEP module
    !# ---
    !# @info
    !# **Brief:** Initialization of SnowWeeklyNCEP module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo  
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = SnowWeeklyNCEPNameList)
    varCommon = varCommon_
  end subroutine initSnowWeeklyNCEP


  function shouldRunSnowWeeklyNCEP() result(shouldRun)
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
  end function shouldRunSnowWeeklyNCEP


  function getOutFileName() result(snowWeeklyNCEPOutFilename)
    !# Gets SnowWeeklyNCEP Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets SnowWeeklyNCEP Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo    
    implicit none
    character(len = maxPathLength) :: snowWeeklyNCEPOutFilename

    snowWeeklyNCEPOutFilename = trim(varCommon%dirPreOut) // trim(var%varName) // trim(var%dateSnowWeeklyNCEP(1:8))
  end function getOutFileName


  function generateSnowWeeklyNCEP() result(isExecOk)
    !# Generates SnowWeeklyNCEP
    !# ---
    !# @info
    !# **Brief:** Generates Snow Weekly NCEP output. This subroutine is the main
    !# method for use this module. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo  
    implicit none
    logical :: isExecOk
    character (len = maxPathLength) :: fileInFileName , fileOutFileName

    isExecOk = .false.

    !call printNameList()

    read (var%dateSnowWeeklyNCEP(1:4), fmt = '(I4)') year
    if (year >= 2015) then
      write(unit = nfprt, fmt = '(A)') "*** WARNING: SnowWeeklyNCEP will not run for year >= 2015 ***"
      return
    endif
 
    allocate(snow(var%xDim,var%yDim))
 
    var%fileSnow(8:9)  = var%dateSnowWeeklyNCEP(9:10)
    var%timeGrADS(1:2) = var%dateSnowWeeklyNCEP(9:10)
    var%timeGrADS(4:5) = var%dateSnowWeeklyNCEP(7:8)
    var%timeGrADS(9:12)= var%dateSnowWeeklyNCEP(1:4)
    read (var%dateSnowWeeklyNCEP(5:6), fmt = '(I2)') month
    var%timeGrADS(6:8) = monthChar(month)

    fileInFileName = trim(varCommon%dirPreIn)//trim(var%fileSnow)//'.'//var%dateSnowWeeklyNCEP
    call msgOut(headerMsg, "Opening as input: " // trim(fileInFileName))
    p_nfsno = openFile(trim(fileInFileName), 'unformatted', 'sequential', -1, 'read', 'old')
    if(p_nfsno < 0) return
    read (unit = p_nfsno) snow 
    close(unit = p_nfsno)
 
    inquire (iolength=lRec) snow
    fileOutFileName = trim(varCommon%dirPreOut)//trim(var%varName)//'.'//var%dateSnowWeeklyNCEP(1:8)
    call msgOut(headerMsg, "Opening as output: " // trim(fileOutFileName))
    p_nfout = openFile(trim(fileOutFileName), 'unformatted', 'direct', lRec, 'write', 'replace')
    if(p_nfout < 0) return
    write(unit = p_nfout, rec = 1) snow 
    close(unit = p_nfout)
 
    if (varCommon%grads) then
      if( .not. generateGrads()) return
    end if

    isExecOk = .true.
  end function generateSnowWeeklyNCEP


  function generateGrads() result(isExecOk)
    !# Generates Grads
    !# ---
    !# @info
    !# **Brief:** Generates Grads. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: sep/2019 </br>
    !# @endinfo 
     implicit none
     logical :: isExecOk

     isExecOk = .false.

     p_nfctl = openFile(trim(varCommon%dirPreOut)//trim(var%varName)//'.'//var%dateSnowWeeklyNCEP(1:8)//'.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
     if(p_nfctl < 0) return

     write(unit = p_nfctl, fmt = '(A)') 'DSET ^'// &
            trim(varCommon%dirPreOut)//trim(var%varName)//'.'//var%dateSnowWeeklyNCEP(1:8)
     write(unit = p_nfctl, fmt = '(A)') '*'
     write(unit = p_nfctl, fmt = '(A)') 'OPTIONS YREV BIG_ENDIAN'
     write(unit = p_nfctl, fmt = '(A)') '*'
     write(unit = p_nfctl, fmt = '(A)') 'UNDEF -999.0'
     write(unit = p_nfctl, fmt = '(A)') '*'
     write(unit = p_nfctl, fmt = '(A)') 'TITLE NCEP Weekly SNOW'
     write(unit = p_nfctl, fmt = '(A)') '*'
     write(unit = p_nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
           'XDEF ',var%xDim,' LINEAR ',0.5_p_r4,360.0_p_r4/real(var%xDim,p_r4)
     write(unit = p_nfctl, fmt = '(A,I5,A,F8.3,F15.10)') &
           'YDEF ',var%yDim,' LINEAR ',-89.5_p_r4,180.0_p_r4/real(var%yDim-1,p_r4)
     write(unit = p_nfctl, fmt = '(A)') 'ZDEF 1 LEVELS 1000'
     write(unit = p_nfctl, fmt = '(A)') 'TDEF 1 LINEAR '//trim(var%timeGrADS)//' 6HR'
     write(unit = p_nfctl, fmt = '(A)') '*'
     write(unit = p_nfctl, fmt = '(A)') 'VARS 1'
     write(unit = p_nfctl, fmt = '(A)') 'SNOW 0 99 NCEP Weekly SNOW [kg/m2]'
     write(unit = p_nfctl, fmt = '(A)') 'ENDVARS'
     close(unit = p_nfctl)
     isExecOk = .true.
  end function generateGrads


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

    write(unit = nfprt, fmt = '(/,A)')  ' &InputDim'
    write(unit = nfprt, fmt = '(A,I6)') '        Idim = ', var%xDim
    write(unit = nfprt, fmt = '(A,I6)') '        Jdim = ', var%yDim
    write(unit = nfprt, fmt = '(A,L6)') '       GrADS = ', varCommon%grads
    write(unit = nfprt, fmt = '(A)')    '    fileSnow = '//trim(var%fileSnow)
    write(unit = nfprt, fmt = '(A)')    '   timeGrADS = '//trim(var%timeGrADS)
    write(unit = nfprt, fmt = '(A)')    '        Date = '//var%dateSnowWeeklyNCEP
    write(unit = nfprt, fmt = '(A,/)')  ' /'
  end subroutine printNameList


end module Mod_SnowWeeklyNCEP
