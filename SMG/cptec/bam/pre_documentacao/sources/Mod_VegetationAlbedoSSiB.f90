!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_VegetationAlbedoSSiB </br></br>
!#
!# **Brief**: Module responsible for the extraction of vegetation and soil 
!# parameters for SSiB </br></br>
!# 
!# The vegetation and soil parameters are obtained from the morphological and 
!# physiological data of the SSiB vegetation (sibalb.form and sibveg.form files).
!# </br></br>
!# 
!# **Files in:**
!#
!# &bull; pre/databcs/sibalb.form </br>
!# &bull; pre/databcs/sibveg.form </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/VegetationSSiB </br>
!# &bull; model/datain/AlbedoSSiB
!# </br></br>
!#
!# **Author**: Jose P. Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Jose P. Bonatti  - version: 1.0.0 </li>
!#  <li>01-08-2007 - Tomita           - version: 1.1.1 </li>
!#  <li>01-04-2018 - Daniel M. Lamosa - version: 2.0.0 </li>
!#  <li>04-02-2020 - Eduardo Khamis   - version: 2.1.0 </li>
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

module Mod_VegetationAlbedoSSiB

  use Mod_FileManager
  use Mod_Namelist, only : varCommonNameListData
  use Mod_Messages, only : msgwarningOut

  implicit none
  
  public :: generateVegetationAlbedoSSiB
  public :: getNameVegetationAlbedoSSiB
  public :: initVegetationAlbedoSSiB
  public :: shouldRunVegetationAlbedoSSiB
  
  private
  include 'files.h'
  include 'pre.h'
  include 'precision.h'
  include 'messages.h'

  !input variables ---------------------------------------------------------------------------
  type VegetationAlbedoSSiBNameListData
    character(len=maxPathLength) :: albedoBCs='sibalb.form'           
    !# SSiB albedo input file name
    character(len=maxPathLength) :: albedoOut='AlbedoSSiB'            
    !# prefix output file Albedo
    character(len=maxPathLength) :: vegetationBCs='sibveg.form'       
    !# vegetation bcs input file name
    character(len=maxPathLength) :: vegetationOut='VegetationSSiB'    
    !# prefix output file Vegetation
  end type VegetationAlbedoSSiBNameListData

  type(varCommonNameListData)            :: varCommon
  type(VegetationAlbedoSSiBNameListData) :: var
  namelist /VegetationAlbedoSSiBNameList/   var

  character(len=*), parameter :: header = 'Vegetation Albedo SSiB          | '


  contains


  function getNameVegetationAlbedoSSiB() result(returnModuleName)
    !# Returns VegetationAlbedoSSiB Module Name
    !# ---
    !# @info
    !# **Brief:** Returns VegetationAlbedoSSiB Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName

    returnModuleName = "VegetationAlbedoSSiB"
  end function getNameVegetationAlbedoSSiB


  subroutine initVegetationAlbedoSSiB(nameListFileUnit, varCommon_)
    !# Initialization of VegetationAlbedoSSiB module
    !# ---
    !# @info
    !# **Brief:** Initialization of VegetationAlbedoSSiB module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    integer, intent(in) :: nameListFileUnit
    type(varCommonNameListData), intent(in) :: varCommon_

    read(unit = nameListFileUnit, nml = VegetationAlbedoSSiBNameList)
    varCommon = varCommon_
  end subroutine initVegetationAlbedoSSiB


  function shouldRunVegetationAlbedoSSiB() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it does
    !# not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getVegetationOutFileName()) .or. .not. fileExists(getAlbedoOutFileName())
  end function shouldRunVegetationAlbedoSSiB


  function getVegetationOutFileName() result(vegetationSSiBOutFilename)
    !# Gets VegetationSSiB Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets VegetationSSiB Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: vegetationSSiBOutFilename

    vegetationSSiBOutFilename = trim(varCommon%dirModelIn) // trim(var%vegetationOut)
  end function getVegetationOutFileName


  function getAlbedoOutFileName() result(albedoSSiBOutFilename)
    !# Gets AlbedoSSiB Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets AlbedoSSiB Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Eduardo Khamis </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    character(len = maxPathLength) :: albedoSSiBOutFilename

    albedoSSiBOutFilename = trim(varCommon%dirModelIn) // trim(var%albedoOut)
  end function getAlbedoOutFileName


  function generateVegetationAlbedoSSiB() result(isExecOk)
    !# Generates Vegetation Albedo SSiB output file
    !# ---
    !# @info
    !# **Brief:** Generates Vegetation Albedo SSiB output file. This subroutine is the
    !# main method for use this module. Only file name of namelist is needed to use it. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# &bull; Eduardo Khamis - changing subroutine to function </br>
    !# **Date**: jan/2020 </br>
    !# @endinfo
    implicit none
    logical :: isExecOk

    isExecOk = .false.
    
    if (.not. vegMorphoPhysioSSiBData()) then
      call msgWarningOut(header, "Error reading sibveg.form or writing VegetationSSiB file")
      return
    end if
    if (.not. albedoSSiBData()) then
      call msgWarningOut(header, "Error reading sibalg.form or writing AlbedoSSiB file")
      return
    end if


    isExecOk = .true.    
  end function generateVegetationAlbedoSSiB


  function vegMorphoPhysioSSiBData () result(isReturnOk)
    !# Generates Vegetation Morphological and Physiological SSiB output file
    !# ---
    !# @info
    !# **Brief:** Generates Vegetation Morphological and Physiological SSiB output file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReturnOk

    integer, parameter :: p_ityp = 13 
    integer, parameter :: p_imon = 12 
    integer, parameter :: p_icg  = 2  
    integer, parameter :: p_iwv  = 3  
    integer, parameter :: p_ild  = 2  
    integer, parameter :: p_idp  = 3  
    integer, parameter :: p_ibd  = 2  

    integer :: i    
    integer :: ntyp 
    integer :: iv   
    integer :: im   
    integer :: k    
    integer :: mmm  
    integer :: vgiFileUnit
    integer :: vgoFileUnit

    ! Vegetation and Soil Parameters
    ! real(kind=p_r4) tran(p_ityp, p_icg, p_iwv, p_ild)
    ! real(kind=p_r4) ref(p_ityp, p_icg, p_iwv, p_ild)
    real(kind=p_r4) rstpar(p_ityp, p_icg, p_iwv)
    ! real(kind=p_r4) soref(p_ityp, p_iwv)
    real(kind=p_r4) chil(p_ityp, p_icg)
    real(kind=p_r4) topt(p_ityp, p_icg)
    real(kind=p_r4) tll(p_ityp, p_icg)
    real(kind=p_r4) tu(p_ityp, p_icg)
    real(kind=p_r4) defac(p_ityp, p_icg)
    real(kind=p_r4) ph1(p_ityp, p_icg)
    real(kind=p_r4) ph2(p_ityp, p_icg)
    real(kind=p_r4) rootd(p_ityp, p_icg)
    ! real(kind=p_r4) rlmax(p_ityp, p_icg)
    ! real(kind=p_r4) rootca(p_ityp, p_icg)
    ! real(kind=p_r4) rplant(p_ityp, p_icg)
    ! real(kind=p_r4) rdres(p_ityp, p_icg)
    real(kind=p_r4) bee(p_ityp)
    real(kind=p_r4) phsat(p_ityp)
    real(kind=p_r4) satco(p_ityp)
    real(kind=p_r4) poros(p_ityp)
    ! real(kind=p_r4) slope(p_ityp)
    real(kind=p_r4) zdepth(p_ityp, p_idp)
    real(kind=p_r4) green(p_ityp, p_imon, p_icg)
    real(kind=p_r4) xcover(p_ityp, p_imon, p_icg)
    real(kind=p_r4) zlt(p_ityp, p_imon, p_icg)
    real(kind=p_r4) x0x(p_ityp, p_imon)
    real(kind=p_r4) xd(p_ityp, p_imon)
    real(kind=p_r4) z2(p_ityp, p_imon)
    real(kind=p_r4) z1(p_ityp, p_imon)
    real(kind=p_r4) xdc(p_ityp, p_imon)
    real(kind=p_r4) xbc(p_ityp, p_imon)
    ! real(kind=p_r4) rootl(p_ityp, p_imon, p_icg)

   rstpar=0.0_p_r4   ;   chil=0.0_p_r4     
   topt=0.0_p_r4     ;   tll=0.0_p_r4      
   tu=0.0_p_r4       ;   defac=0.0_p_r4    
   ph1=0.0_p_r4      ;   ph2 =0.0_p_r4     
   rootd  =0.0_p_r4  ;   bee =0.0_p_r4     
   phsat=0.0_p_r4    ;   satco =0.0_p_r4   
   poros =0.0_p_r4   ;   zdepth =0.0_p_r4  
   green =0.0_p_r4   ;   xcover=0.0_p_r4   
   zlt=0.0_p_r4      ;   x0x =0.0_p_r4     
   xd =0.0_p_r4      ;   z2  =0.0_p_r4     
   z1=0.0_p_r4       ;   xdc =0.0_p_r4     
   xbc =0.0_p_r4     


    isReturnOk = .true.

    ! Open and read VegetationBCs --------------------------------------------------
    vgiFileUnit = openFile(trim(varCommon%dirBCs)//trim(var%vegetationBCs), 'formatted', 'sequential', -1, 'read', 'old')
    if (vgiFileUnit < 0) return

    read (unit=vgiFileUnit, fmt='(a)')
    do
      read(unit=vgiFileUnit, fmt='(i3)') ntyp
      read(unit=vgiFileUnit, fmt=*) ((rstpar(ntyp,iv,im), im=1,3), iv=1,2), &
                                (chil(ntyp,iv),   iv=1,2),              &
                                (topt(ntyp,iv),   iv=1,2),              &
                                (tll(ntyp,iv),    iv=1,2),              &
                                (tu(ntyp,iv),     iv=1,2),              &
                                (defac(ntyp,iv),  iv=1,2),              &
                                (ph1(ntyp,iv),    iv=1,2),              &
                                (ph2(ntyp,iv),    iv=1,2),              &
                                (rootd (ntyp,iv), iv=1,2),              &
                                bee(ntyp),                              &
                                phsat(ntyp),                            &
                                satco(ntyp),                            &
                                poros(ntyp),                            &
                                (zdepth(ntyp,k), k=1,3)
      do i=1, 12
        read (unit=vgiFileUnit, fmt=*) mmm, &
                                   (zlt(ntyp,i,iv), iv=1,2),    &
                                   (green(ntyp,i,iv), iv=1,2),  &
                                   z2(ntyp,i),                  &
                                   z1(ntyp,i),                  &
                                   (xcover(ntyp,i,iv), iv=1,2), &
                                   x0x(ntyp,i),                 &
                                   xd(ntyp,i),                  &
                                   xbc(ntyp,i),                 &
                                   xdc(ntyp,i)
      end do
      if(ntyp == p_ityp) exit
    end do
    close(unit=vgiFileUnit)
    ! ------------------------------------------------------------------------------
    
    ! Open and write Vegetation ----------------------------------------------------
    vgoFileUnit = openFile(trim(varCommon%dirModelIn)//trim(var%vegetationOut), 'unformatted', 'sequential', -1, 'write', 'replace')
    if(vgoFileUnit < 0) return

    write(unit=vgoFileUnit) rstpar, chil, topt, tll, tu, defac, ph1, ph2, &
                        rootd, bee, phsat, satco, poros, zdepth
    write(unit=vgoFileUnit) green, xcover, zlt, x0x, xd, z2, z1, xdc, xbc
    close(unit=vgoFileUnit)
    
    ! Print information in screen. To be commented by Daniel M. Lamosa
    !do ntyp=1, p_ityp
    !  write(unit=p_nfprt, fmt='(I3)') ntyp
    !  write(unit=p_nfprt, fmt='(6f12.4,/,2f8.2,/,6f8.1,/,2f10.4,4f8.2,/,&
    !                          &2f6.1,/,2f10.4,/,e12.5,f10.4,/,3e12.5)') &
    !                      ((rstpar(ntyp, iv, im), im=1,3), iv=1,2),     &
    !                      (chil(ntyp, iv), iv=1,2),                     &
    !                      (topt(ntyp, iv), iv=1,2),                     &
    !                      (tll(ntyp, iv), iv=1,2),                      &
    !                      (tu(ntyp, iv), iv=1,2),                       &
    !                      (defac(ntyp, iv), iv=1,2),                    &
    !                      (ph1(ntyp, iv), iv=1,2),                      &
    !                      (ph2(ntyp, iv), iv=1,2),                      &
    !                      (rootd(ntyp, iv), iv=1,2),                    &
    !                      bee(ntyp),                                    &
    !                      phsat(ntyp),                                  &
    !                      satco(ntyp),                                  &
    !                      poros(ntyp),                                  &
    !                      (zdepth(ntyp,k), k=1,3)
    !  do i=1, 12
    !    write(unit=p_nfprt, fmt='(i3,/,2f10.5,2f10.5,2f8.4,/,&
    !                             &2f10.5,2f10.5,/,2f10.2)')  &
    !                        i,                               &
    !                        (zlt(ntyp, i, iv), iv=1,2),      &
    !                        (green(ntyp,i,iv), iv=1,2),      &
    !                        z2(ntyp, i),                     &
    !                        z1(ntyp, i),                     &
    !                        (xcover(ntyp, i, iv), iv=1,2),   &
    !                        x0x(ntyp, i),                    &
    !                        xd(ntyp, i),                     &
    !                        xbc(ntyp, i),                    &
    !                        xdc(ntyp, i)
    !  end do
    !end do
    
    ! ------------------------------------------------------------------------------

    isReturnOk = .true.
  end function vegMorphoPhysioSSiBData


  function albedoSSiBData() result(isReturnOk)
    !# Writes Albedo SSiB output file
    !# ---
    !# @info
    !# **Brief:** Writes Albedo SSiB output file. </br>
    !# **Authors**: </br>
    !# &bull; Daniel M. Lamosa </br>
    !# **Date**: mar/2018 </br>
    !# @endinfo
    implicit none
    logical :: isReturnOk

    integer, parameter :: p_ityp = 13 
    integer, parameter :: p_imon = 12 
    integer, parameter :: p_icg  = 2  
    integer, parameter :: p_njj  = 6  
    integer, parameter :: p_ild  = 2  
    integer, parameter :: p_ibd  = 2  
    integer, parameter :: p_nj   = 9  
    integer, parameter :: p_nk   = 3  

    integer :: itp 
    integer :: l   
    integer :: j   
    integer :: i   
    integer :: k   
    integer :: ml  
    integer :: mon 
    integer :: aliFileUnit
    integer :: aloFileUnit
   
    ! Albedo Parameters
    real(kind=p_r4) :: cedfu(p_ityp, p_imon, p_nj)
    real(kind=p_r4) :: cedir(p_ityp, p_imon, p_nj, 3)
    real(kind=p_r4) :: cedfu1(2, p_ityp, p_imon, p_njj, 3)
    real(kind=p_r4) :: cedir1(2, p_ityp, p_imon, p_njj, p_nk, 3)
    real(kind=p_r4) :: cedfu2(2, p_ityp, p_imon, p_njj, 3)
    real(kind=p_r4) :: cedir2(2, p_ityp, p_imon, p_njj, p_nk, 3)
    real(kind=p_r4) :: cledfu(p_ityp, p_imon, p_nj)
    real(kind=p_r4) :: cledir(p_ityp, p_imon, p_nj, 3)
    real(kind=p_r4) :: cether(p_ityp, p_imon, 2)
    real(kind=p_r4) :: xmiu(p_imon, p_nk)
    real(kind=p_r4) :: xmiw(p_imon, p_nk)

    isReturnOk = .false.

    cedfu = 0.0_p_r4
    cedir = 0.0_p_r4
    cedfu1 = 0.0_p_r4
    cedir1 = 0.0_p_r4
    cedfu2 = 0.0_p_r4
    cedir2 = 0.0_p_r4
    cledfu = 0.0_p_r4
    cledir = 0.0_p_r4
    xmiu = 0.0_p_r4
    cether = 0.0_p_r4
    xmiw = 0.0_p_r4

    ! open and read sib albedo (sibalb.form) --------------------------------------------
    aliFileUnit = openFile(trim(varCommon%dirBCs)//trim(var%albedoBCs), 'formatted', 'sequential', -1, 'read', 'old')
    if (aliFileUnit < 0) return

    do itp=1, 12
      do l=1, 2
        read(unit=aliFileUnit, fmt='(6e12.5)') (cether(itp, mon, l), mon=1,12)
      end do
      do j=1, p_nj
        read(unit=aliFileUnit, fmt='(6e12.5)') (cledfu(itp, mon, j), mon=1,12)
        do l=1, 3
          read(unit=aliFileUnit, fmt='(6e12.5)') (cledir(itp,mon,j,l), mon=1,12)
        end do
      end do
      do i=1, p_nk
        read(unit=aliFileUnit, fmt='(6e12.5)') (xmiu(mon, i), mon=1,12)
      end do
      do j=1, p_nj
        read(unit=aliFileUnit, fmt='(6e12.5)') (cedfu(itp, mon, j), mon=1,12)
        do l=1, 3
          read(unit=aliFileUnit, fmt='(6e12.5)') (cedir(itp, mon, j, l), mon=1,12)
        end do
      end do
      do ml=1, 2
        do j=1, p_njj
          do l=1, 3
            read(unit=aliFileUnit, fmt='(6e12.5)') (cedfu1(ml, itp, mon, j, l), mon=1,12)
            read(unit=aliFileUnit, fmt='(6e12.5)') (cedfu2(ml, itp, mon, j, l), mon=1,12)
          end do
        end do
        do j=1, p_njj
          do k=1, p_nk
            do l=1, 3
              read(unit=aliFileUnit, fmt='(6e12.5)') (cedir1(ml, itp, mon, j, k, l), mon=1,12)
              read(unit=aliFileUnit, fmt='(6e12.5)') (cedir2(ml, itp, mon, j, k, l), mon=1,12)
            end do
          end do
        end do
      end do
    end do
    do l=1, 2
      read(unit=aliFileUnit, fmt='(7e12.5)') (cether(13, mon, l), mon=1,7)
    end do
    do j=1, p_nj
      read(unit=aliFileUnit, fmt='(7e12.5)') (cledfu(13, mon, j), mon=1,7)
      do l=1,3
        read(unit=aliFileUnit, fmt='(7e12.5)') (cledir(13, mon, j, l), mon=1,7)
      end do
    end do
    do i=1, p_nk
      read(unit=aliFileUnit, fmt='(7e12.5)') (xmiw(mon, i), mon=1,7)
    end do
    
    ! Print data in screen. Commented by Daniel M. Lamosa
    !do i=1, p_nk
    !  write(unit=p_nfprt, fmt='(7e12.5)') (xmiw(mon, i), mon=1,7)
    !end do
    
    do j=1, p_nj
      read(unit=aliFileUnit, fmt='(7e12.5)') (cedfu(13, mon, j), mon=1,7)
      do l=1, 3
         read(unit=aliFileUnit, fmt='(7e12.5)') (cedir(13, mon, j, l), mon=1,7)
      end do
    end do
    do ml=1, 2
      do j=1, p_njj
        do l=1, 3
          read(unit=aliFileUnit, fmt='(7e12.5)') (cedfu1(ml, 13, mon, j, l), mon=1,7)
          read(unit=aliFileUnit, fmt='(7e12.5)') (cedfu2(ml, 13, mon, j, l), mon=1,7)
        end do
      end do
      do j=1, p_njj
        do k=1, p_nk
          do l=1, 3
            read(unit=aliFileUnit, fmt='(7e12.5)') (cedir1(ml, 13, mon, j, k, l), mon=1,7)
            read(unit=aliFileUnit, fmt='(7e12.5)') (cedir2(ml, 13, mon, j, k, l), mon=1,7)
          end do
        end do
      end do
    end do
    close(unit=aliFileUnit)
    ! ----------------------------------------------------------------------------------
    
    ! Open and write Albedo output file
    aloFileUnit = openFile(trim(varCommon%dirModelIn)//trim(var%albedoOut), 'unformatted', 'sequential', -1, 'write', 'replace')
    if (aloFileUnit < 0) return

    write(unit=aloFileUnit) cedfu, cedir, cedfu1, cedir1, cedfu2, cedir2, &
                        cledfu, cledir, xmiu, cether, xmiw
    close(unit=aloFileUnit)


    isReturnOk = .true.
  end function albedoSSiBData


end module Mod_VegetationAlbedoSSiB
