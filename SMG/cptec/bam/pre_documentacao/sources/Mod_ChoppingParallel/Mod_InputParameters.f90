!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_InputParameters </br></br>
!#
!# **Brief**: Module containing most of Chopping parameters. </br></br>
!#
!# **Files in:**
!#
!# &bull; ? </br></br>
!#
!# **Files out:**
!#
!# &bull; ? </br></br>
!# 
!# **Author**: Paulo Bonatti </br>
!#
!# **Version**: 2.1.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2008 - Paulo Bonatti  - version: 1.3.0 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 - some adaptations for modularizing Chopping </li>
!#  <li>09-10-2019 - Eduardo Khamis - version: 2.1.0 - changing for operational Chopping </li>
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


module Mod_InputParameters

  use nemsio_module_mpi
  use nemsio_gfs
  USE netcdf

  use Mod_Parallelism_Group_Chopping, only : &
    myId &
    , mpiCommGroup &
    , fatalError

  use Mod_FileManager, only: &
    openFile &
    , getFileUnit

  implicit none

  private
  include 'precision.h'
  include 'messages.h'
  include 'mpif.h'
  include 'files.h'

  type ChoppingNameListData
    integer :: mEndInp, kMaxInp, mEndOut, kMaxOut, mEndMin, mEndCut, iter, nProc_vert, ibdim_size, tamBlock
    real (kind = p_r8) :: smthPerCut
    logical :: getOzone, getTracers, grADS, grADSOnly, gdasOnly, smoothTopo, rmGANL, linearGrid, givenfouriergroups
    character (len = 10) :: dateLabel
    character (len = 2) :: utc
    character (len = 16) :: nCepName
    character (len = 16) :: StrFormat
    character (len = 500) :: dataGDAS, gdasInp
    character(len = maxPathLength) :: dirInp, dirOut, dirTop, dirSig, dirGrd, dGDInp
    character(len = maxPathLength) :: dirPreOut = './'                       
    !# output data directory
  end type ChoppingNameListData

  public :: ChoppingNameListData
  public :: InitParameters
  public :: fillDataOutFileName
  public :: check

  real (kind = p_r8), parameter, public :: pai = 3.14159265358979_p_r8
  real (kind = p_r8), parameter, public :: twomg = 1.458492e-4_p_r8

  integer, public :: ImaxInp, JmaxInp, &
    ImaxOut, JmaxOut

  integer, public :: Mnwv2Inp, Mnwv3Inp, &
    Mend1Out, Mend2Out, Mend3Out, &
    Mnwv2Out, Mnwv3Out, Mnwv0Out, Mnwv1Out, &
    ImxOut, JmaxhfOut, KmaxInpp, KmaxOutp, &
    NTracers, Kdim, ICaseRec, ICaseDec, &
    MFactorFourier, MTrigsFourier, Iter

  integer (kind = p_i8), public :: IDVCInp
  integer (kind = p_i8), public :: IDSLInp
  integer (kind = p_i8), public :: ForecastDay

  real (kind = p_r4), public :: TimeOfDay

  real (kind = p_r8), public :: cTv, SmthPerCut

  character (len = 500), public :: DataCPT, DataInp, DataOut, DataOup, DataTop, DataTopG, DataSig, DataSigInp
  character (len = 500), public :: OzonInp, TracInp, OzonOut, TracOut

  real (kind = p_r8), public :: EMRad   
  !# Earth Mean Radius (m)
  real (kind = p_r8), public :: EMRad1  
  !# 1/EMRad (1/m)
  real (kind = p_r8), public :: EMRad12 
  !# EMRad1**2 (1/m2)
  real (kind = p_r8), public :: EMRad2  
  !# EMRad**2 (m2)
  real (kind = p_r8), public :: Grav    
  !# Gravity Acceleration (m2/s)
  real (kind = p_r8), public :: Rd      
  !# Dry Air Gas Constant (m2/s2/K)
  real (kind = p_r8), public :: Rv      
  !# Water Vapour Air Gas Constant (m2/s2/K)
  real (kind = p_r8), public :: Cp      
  !# Dry Air Heat Capacity (m2/s2/K)
  real (kind = p_r8), public :: Lc      
  !# Latent Heat of Condensation (m2/s2)
  real (kind = p_r8), public :: Gama    
  !# Mean Atmospheric Lapse Rate (K/m)
  real (kind = p_r8), public :: GEps    
  !# precision For Constante Lapse Rate (No dim)
  real (kind = p_r8), public, parameter :: p0mb = 1000.0_p_r8       
  !# surface pressure (mbar)
  real (kind = p_r8), public, parameter :: p0Pa = 100000.0_p_r8       
  !# surface pressure (Pa)

  real (kind = p_r8), public :: RoCp    
  !# Rd over Cp
  real (kind = p_r8), public :: RoCp1   
  !# RoCp + 1
  real (kind = p_r8), public :: RoCpr   
  !# 1 / RoCp

  character (len = 10), public :: TruncInp, TruncOut
  character (len = 3), dimension (12), public :: MonChar

  character(len = *), parameter :: headerMsg = 'Chopping            | '

contains


  function InitParameters (itertrc, dtoChp)  result(isExecOk)
    !# Initializes parameters
    !# ---
    !# @info
    !# **Brief:** Initializes parameters. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none
    integer, intent(in) :: itertrc
    type(ChoppingNameListData), intent(inout) :: dtoChp
    logical :: isExecOk

    !    Namelist parameters
    !    integer :: dtoChp%mEndInp, dtoChp%kMaxInp, , intent(in) dtoChp%mEndOut, dtoChp%kMaxOut, mEndMin, mEndCut, iter, nProc_vert, ibdim_size, tamBlock
    !    real (kind = p_r8) , intent(in) :: smthPerCut
    !    logical   :: getOzone, getTracers, , intent(in) grADS, grADSOnly, dtoChp%gdasOnly, smoothTopo, rmGANL, linearGrid, givenfouriergroups
    !    character (len = 10) , intent(in) :: dtoChp%dateLabel
    !    character (len = 16) , intent(in) :: dtoChp%ncepName,dtoChp%StrFormat
    !    character (len = 3) , intent(in) :: prefix
    !    character (len = 500), , intent(in) public :: dtoChp%DataGDAS
    !    character(len = maxPathLength) , intent(in) :: dirPreIn = './'                        
    ! input data directory
    !    character(len = maxPathLength) , intent(in) :: dirPreOut = './'                       
    ! output data directory
    !    character (len = 2) , intent(in) :: dtoChp%utc
    !    character (len = 500), , intent(in) public :: dtoChp%DataGDAS

    integer :: im, jm, itrc, nficr

    character (len = 10) :: TrGrdOut
    character (len = 6) :: TrGrdOutGaus

    logical :: ExistGANL

    isExecOk = .false.   
    MonChar = (/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
      'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /)

    ! TODO - check if is desired to have default values
    !    dtoChp%mEndOut = 213     ! Spectral Horizontal Resolution of Output data
    !    dtoChp%kMaxOut = 42      ! Number of Layers of Output data
    !    MendMin = 127     ! Minimum Spectral Resolution For Doing Topography Smoothing
    !    MendCut = 0       ! Spectral Resolution Cut Off for Topography Smoothing
    !    Iter = 10         ! Number of Iteractions in Topography Smoothing
    !    nproc_vert = 1         ! Number of processors to be used in the vertical
    ! (if givenfouriergroups set to TRUE)
    !    ibdim_size = 192       ! size of basic data block (ibmax)
    !    tamBlock = 512     ! number of fft's allocated in each block
    !    SmthPerCut = 0.12_p_r8 ! Percentage for Topography Smoothing
    !    GetOzone = .false.   ! Flag to Produce Ozone Files
    !    GetTracers = .false. ! Flag to Produce Tracers Files
    !    GrADS = .true.       ! Flag to Produce GrADS Files
    !    dtoChp%DataGDAS = "Grid"
    !    GrADSOnly = .false.  ! Flag to only Produce GrADS Files (do not Produce Inputs for Model)
    !    dtoChp%gdasOnly = .false.   ! Flag to only Produce Input CPTEC Analysis file
    !    SmoothTopo = .true.  ! Flag to Performe Topography Smoothing
    !    RmGANL = .false.     ! Flag to Remove GANL file if Desired
    !    LinearGrid = .false. ! Flag to Set Linear (T) or Quadratic Gaussian Grid (F)
    !    givenfouriergroups = .false.! false if processor division should be automatic
    !    dtoChp%dateLabel = 'yyyymmddhh' ! Date: yyyymmddhh or dtoChp%dateLabel='        hh'
    !       If Year (yyyy), Month (mm) and Day (dd) Are Unknown
    !    dtoChp%utc = 'hh'               ! dtoChp%utc Hour: hh, Must Be Given if dtoChp%dateLabel='          ', else dtoChp%utc=' '
    !    dtoChp%ncepName = 'gdas1 '      ! NCEP Analysis Preffix for Input file Name
    !    dtoChp%StrFormat = 'old '      ! Get NCEP Analysis and dumping to 'old' format or 'new' format
    !    DirMain = './ '          ! Main User data Directory
    !    DirHome = './ '          ! Home User Sources Directory
    !    prefix = 'NMC'
    ! dtoChp%mEndInp : Spectral Horizontal Resolution of Input Data
    ! dtoChp%kMaxInp : Number of Layers of Input Data

        if (myid.eq.0) then
          write (unit = p_nfprt, fmt = '(A)')      ' '
          write (unit = p_nfprt, fmt = '(A)')      ' &ChopNML'
          write (unit = p_nfprt, fmt = '(A,I5)')   '  dtoChp%mEndInp    = ', dtoChp%mEndInp
          write (unit = p_nfprt, fmt = '(A,I5)')   '  dtoChp%mEndOut    = ', dtoChp%mEndOut
          write (unit = p_nfprt, fmt = '(A,I5)')   '  dtoChp%kMaxOut    = ', dtoChp%kMaxOut
          write (unit = p_nfprt, fmt = '(A,I5)')   '  MendMin    = ', dtoChp%MendMin
          write (unit = p_nfprt, fmt = '(A,I5)')   '  MendCut    = ', dtoChp%MendCut
          write (unit = p_nfprt, fmt = '(A,I5)')   '  Iter       = ', dtoChp%Iter
          write (unit = p_nfprt, fmt = '(A,I5)')   '  nproc_vert = ', dtoChp%nproc_vert
          write (unit = p_nfprt, fmt = '(A,I5)')   '  ibdim_size = ', dtoChp%ibdim_size
          write (unit = p_nfprt, fmt = '(A,I5)')   '  tamBlock   = ', dtoChp%tamBlock
          write (unit = p_nfprt, fmt = '(A,F7.3)') '  SmthPerCut = ', dtoChp%SmthPerCut
          write (unit = p_nfprt, fmt = '(A,L6)')   '  GetOzone   = ', dtoChp%GetOzone
          write (unit = p_nfprt, fmt = '(A,L6)')   '  GetTracers = ', dtoChp%GetTracers
          write (unit = p_nfprt, fmt = '(A,L6)')   '  GrADS      = ', dtoChp%GrADS
          write (unit = p_nfprt, fmt = '(A,A )')   '  dtoChp%DataGDAS   = ', dtoChp%DataGDAS
          write (unit = p_nfprt, fmt = '(A,L6)')   '  GrADSOnly  = ', dtoChp%GrADSOnly
          write (unit = p_nfprt, fmt = '(A,L6)')   '  dtoChp%gdasOnly   = ', dtoChp%gdasOnly
          write (unit = p_nfprt, fmt = '(A,L6)')   '  SmoothTopo = ', dtoChp%SmoothTopo
          write (unit = p_nfprt, fmt = '(A,L6)')   '  LinearGrid = ', dtoChp%LinearGrid
          write (unit = p_nfprt, fmt = '(A,L6)')   '  givenfouriergroups = ', dtoChp%givenfouriergroups
          write (unit = p_nfprt, fmt = '(2A)')     '  dtoChp%dateLabel  = ', dtoChp%dateLabel
          write (unit = p_nfprt, fmt = '(2A)')     '  dtoChp%utc        = ', dtoChp%utc
          write (unit = p_nfprt, fmt = '(2A)')     '  dtoChp%ncepName   = ', trim(dtoChp%ncepName)
          write (unit = p_nfprt, fmt = '(2A)')     '  dtoChp%StrFormat   = ', trim(dtoChp%StrFormat)
     !     write (unit = p_nfprt, fmt = '(2A)')     '  DirMain    = ', trim(DirMain)
     !     write (unit = p_nfprt, fmt = '(2A)')     '  DirHome    = ', trim(DirHome)
     !     write (unit = p_nfprt, fmt = '(2A)')     '  prefix     = ', trim(dtoChp%prefix)
          write (unit = p_nfprt, fmt = '(A)')      ' /'
        endif
    !
    !    DGDInp = trim(DirMain) // 'pre/datain/ '
    !    dtoChp%dirInp = trim(DirMain) // 'model/datain/ '
    !    dtoChp%dirOut = trim(DirMain) // 'model/datain/ '
    !    dtoChp%dirTop = trim(DirMain) // 'pre/dataTop/ '
    !    !  dtoChp%dirTop=trim(DirMain)//'pre/dataout/ '
    !    ! dtoChp%dirSig=trim(DirHome)//'sources/Chopping/ '
    !    dtoChp%dirSig = trim(DirMain) // 'pre/datain/ '
    !    dtoChp%dirGrd = trim(DirMain) // 'pre/dataout/ '

    if(trim(dtoChp%DataGDAS) == 'Spec')then
      if (dtoChp%dateLabel == '          ') then
        dtoChp%gDASInp = trim(dtoChp%ncepName) // '.T' // dtoChp%utc // 'Z.SAnl'
        ! Input NCEP Analysis file Name Without dtoChp%dateLabel
      else
        dtoChp%utc = dtoChp%dateLabel(9:10)
        dtoChp%gDASInp = trim(dtoChp%ncepName) // '.T' // dtoChp%utc // 'Z.SAnl.' // dtoChp%dateLabel 
        ! Input NCEP Analysis file Name
      endif
    else if(trim(dtoChp%DataGDAS) == 'Grid')then
      if (dtoChp%dateLabel == '          ') then
        dtoChp%gDASInp = trim(dtoChp%ncepName) // '.T' // dtoChp%utc // 'Z.atmanl.nemsio.' 
        ! Input NCEP Analysis file Name
      else
        dtoChp%utc = dtoChp%dateLabel(9:10)
        dtoChp%gDASInp = trim(dtoChp%ncepName) // '.T' // dtoChp%utc // 'Z.atmanl.nemsio.' // dtoChp%dateLabel 
        ! Input NCEP Analysis file Name
      endif
    else if(trim(dtoChp%DataGDAS) == 'Netcdf')then
      if (dtoChp%dateLabel == '          ') then
        dtoChp%gDASInp=trim(dtoChp%ncepName)//'.T'//dtoChp%utc//'Z.atmanl.netcdf.' ! Input NCEP Analysis File Name
      else
        dtoChp%utc=dtoChp%dateLabel(9:10)
        dtoChp%gDASInp=TRIM(dtoChp%ncepName)//'.T'//dtoChp%utc//'Z.atmanl.netcdf.'//dtoChp%dateLabel ! Input NCEP Analysis File Name
      end if
    endif

    if(.not.dtoChp%gdasOnly)then
      if(.not. GetGDASDateLabelRes (dtoChp)) return
    endif

    if (myid.eq.0) then
      write (unit = p_nfprt, fmt = '(/,A,I5)')   '  mEndInp    = ', dtoChp%mEndInp
      write (unit = p_nfprt, fmt = '(A,I5,/)')   '  kMaxInp    = ', dtoChp%kMaxInp
    endif

    if(trim(dtoChp%DataGDAS) == 'Spec') then
      call GetImaxJmax (dtoChp%mEndInp, ImaxInp, JmaxInp, dtoChp%linearGrid)
      call GetImaxJmax (dtoChp%mEndOut, ImaxOut, JmaxOut, dtoChp%linearGrid)

    else if (trim(dtoChp%DataGDAS) == 'Grid' .or. trim(dtoChp%DataGDAS) == 'Netcdf') then
      if(itertrc ==1) then
        do itrc = 1, 40000
          call GetImaxJmax (itrc, im, jm, dtoChp%linearGrid)
          if((im > ImaxInp) .and. (jm > JmaxInp))exit
        enddo
        ImaxOut = ImaxInp
        JmaxOut = JmaxInp
        dtoChp%mEndInp = itrc
        dtoChp%mEndOut = itrc
      else
        if(.not. dtoChp%gdasOnly) then
          do itrc = 1, 40000
            call GetImaxJmax (itrc, im, jm, dtoChp%linearGrid)
            if((im > ImaxInp) .and. (jm > JmaxInp))exit
          enddo
          dtoChp%mEndInp = itrc
          call GetImaxJmax (dtoChp%mEndOut, ImaxOut, JmaxOut, dtoChp%linearGrid)
        else
          call GetImaxJmax (dtoChp%mEndOut, ImaxOut, JmaxOut, dtoChp%linearGrid)
        endif
      endif
    endif

    if (dtoChp%linearGrid) then
      TruncInp = 'TL    L   '
      TruncOut = 'TL    L   '
    else
      TruncInp = 'TQ    L   '
      TruncOut = 'TQ    L   '
    endif

    write (TruncInp(3:6), fmt = '(I4.4)') dtoChp%mEndInp
    write (TruncInp(8:10), fmt = '(I3.3)') dtoChp%kMaxInp
    write (TruncOut(3:6), fmt = '(I4.4)') dtoChp%mEndOut
    write (TruncOut(8:10), fmt = '(I3.3)') dtoChp%kMaxOut

    TrGrdOut = 'G     L   '
    TrGrdOutGaus = 'G     '
    write (TrGrdOutGaus(2:6), fmt = '(I5.5)') JmaxOut
    write (TrGrdOut(2:6), fmt = '(I5.5)') JmaxOut
    write (TrGrdOut(8:10), fmt = '(I3.3)') dtoChp%kMaxOut

    call fillDataOutFileName(dtoChp)                                
    ! Output Topo-Smoothed CPTEC Analysis file Name
    if (dtoChp%smoothTopo) then
      DataInp = 'GANLSMT' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
      ! Input Topo-Smoothed CPTEC Analysis file Name
      DataOup = 'GANLSMT' // dtoChp%dateLabel // 'S.unf.' // TruncOut 
      ! Output Topo-Smoothed CPTEC Analysis file Name
      OzonInp = 'OZONSMT' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
      ! Input Ozone file Name
      TracInp = 'TRACSMT' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
      ! Input Tracers file Name
      OzonOut = 'OZONSMT' // dtoChp%dateLabel // 'S.grd.' // TrGrdOut 
      ! Grid Ouput Ozone file Name
      TracOut = 'TRACSMT' // dtoChp%dateLabel // 'S.grd.' // TrGrdOut 
      ! Grid Ouput Tracers file Name
    else
      DataInp = 'GANLNMC' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
      ! Input No Topo-Smoothed CPTEC Analysis file Name
      DataOup = 'GANLNMC' // dtoChp%dateLabel // 'S.unf.' // TruncOut 
      ! Output Topo-Smoothed CPTEC Analysis file Name
      OzonInp = 'OZONNMC' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
      ! Input Ozone file Name
      TracInp = 'TRACNMC' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
      ! Input Tracers file Name
      OzonOut = 'OZONNMC' // dtoChp%dateLabel // 'S.grd.' // TrGrdOut 
      ! Grid Ouput Ozone file Name
      TracOut = 'TRACNMC' // dtoChp%dateLabel // 'S.grd.' // TrGrdOut 
      ! Grid Ouput Tracers file Name
    endif
    DataCPT = 'GANLCPT' // dtoChp%dateLabel // 'S.unf.' // TruncInp 
    ! Input CPTEC No Topo-Smoothed Analysis file Name

    DataTop = 'Topography2.' // TruncOut(1:6)         
    ! Input Topography Spec file Name
    DataTopG = 'Topography2.' // trim(TrGrdOutGaus)   
    ! Input Topography Grid file Name
    DataSigInp = 'HybridLevels.' // TruncInp(7:10)     
    ! Delta Sigma file Input
    ! DataSig = 'HybridLevels.' // TruncOut(7:10)        
    ! HybridLevels file Input
    DataSig = 'HybridLevels.' // TruncOut(7:10)        
    ! Delta Sigma file Output

    inquire (file = trim(dtoChp%dirInp) // trim(DataInp), exist = ExistGANL)
    if (ExistGANL .and. dtoChp%rmGAnl) then
      if (myid.eq.0) then
        nficr = getFileUnit()
        open    (unit = nficr, file = trim(dtoChp%dirInp) // trim(DataInp))
        close   (unit = nficr, status = 'DELETE')
        inquire (file = trim(dtoChp%dirInp) // trim(DataInp), exist = ExistGANL)
        write   (unit = p_nfprt, fmt = '(/,A)') ' file Removed if false: '
        write   (unit = p_nfprt, fmt = '(A,L6)')  trim(dtoChp%dirInp) // trim(DataInp), ExistGANL
      else
        ExistGANL = .false.
      endif
    endif

    if (myid.eq.0) then
      write (unit = p_nfprt, fmt = '(A)')      ' '
      write (unit = p_nfprt, fmt = '(A,I5)')   '  ImaxInp    = ', ImaxInp
      write (unit = p_nfprt, fmt = '(A,I5)')   '  JmaxInp    = ', JmaxInp
      write (unit = p_nfprt, fmt = '(A,I5)')   '  ImaxOut    = ', ImaxOut
      write (unit = p_nfprt, fmt = '(A,I5)')   '  JmaxOut    = ', JmaxOut

      write (unit = p_nfprt, fmt = '(A)')      ' '
      write (unit = p_nfprt, fmt = '(2A)')     '  DataCPT    = ', trim(DataCPT)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirInp     = ', trim(dtoChp%dirInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  DataInp    = ', trim(DataInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirOut     = ', trim(dtoChp%dirOut)
      write (unit = p_nfprt, fmt = '(2A)')     '  DataOut    = ', trim(DataOut)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirTop     = ', trim(dtoChp%dirTop)
      write (unit = p_nfprt, fmt = '(2A)')     '  DataTop    = ', trim(DataTop)
      write (unit = p_nfprt, fmt = '(2A)')     '  DataTopG   = ', trim(DataTopG)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirSig     = ', trim(dtoChp%dirSig)
      write (unit = p_nfprt, fmt = '(2A)')     '  DataSig    = ', trim(DataSig)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirGrd     = ', trim(dtoChp%dirGrd)
      write (unit = p_nfprt, fmt = '(2A)')     '  DGDInp     = ', trim(dtoChp%dgdInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  gDASInp    = ', trim(dtoChp%gDASInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  OzonInp    = ', trim(OzonInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  TracInp    = ', trim(TracInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  OzonOut    = ', trim(OzonOut)
      write (unit = p_nfprt, fmt = '(2A)')     '  TracOut    = ', trim(TracOut)

      write (unit = p_nfprt, fmt = '(/,A)')    '  Numerical precision (kind): '
      write (unit = p_nfprt, fmt = '(A,I5)')   '          p_i8 = ', p_i8
      write (unit = p_nfprt, fmt = '(A,I5)')   '          p_r4 = ', p_r4
      write (unit = p_nfprt, fmt = '(A,I5)')   '          p_r8 = ', p_r8
      write (unit = p_nfprt, fmt = '(A)')      ' '
    endif

    KmaxInpp = dtoChp%kMaxInp + 1
    KmaxOutp = dtoChp%kMaxOut + 1

    NTracers = 1
    Kdim = 1
    ICaseRec = -1
    ICaseDec = 1

    if (dtoChp%mEndCut <= 0 .or. dtoChp%mEndCut > dtoChp%mEndOut) dtoChp%mEndCut = dtoChp%mEndOut

    IDVCInp = 0_p_i8
    IDSLInp = 0_p_i8

    ForecastDay = 0_p_i8
    TimeOfDay = 0.0_p_r4
    EMRad = 6.37E6_p_r8
    EMRad1 = 1.0_p_r8 / EMRad
    EMRad12 = EMRad1 * EMRad1
    EMRad2 = EMRad * EMRad

    Grav = 9.80665_p_r8
    Rd = 287.05_p_r8
    Rv = 461.50_p_r8
    Cp = 1004.6_p_r8
    Lc = 2.5E6_p_r8
    Gama = -6.5E-3_p_r8
    GEps = 1.E-9_p_r8
    cTv = Rv / Rd - 1.0_p_r8

    RoCp = Rd / Cp
    RoCp1 = RoCp + 1.0_p_r8
    RoCpr = 1.0_p_r8 / RoCp
    isExecOk = .true.
  end function InitParameters


  subroutine fillDataOutFileName(dtoChp)
    !# Fills Data Out File Name
    !# ---
    !# @info
    !# **Brief:** Fills Data Out File Name. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    type(ChoppingNameListData), intent(in) :: dtoChp

    if (dtoChp%linearGrid) then
      TruncOut = 'TL    L   '
    else
      TruncOut = 'TQ    L   '
    endif

    write (TruncOut(3:6), fmt = '(I4.4)') dtoChp%mEndOut
    write (TruncOut(8:10), fmt = '(I3.3)') dtoChp%kMaxOut

    if (dtoChp%smoothTopo) then
      DataOut = 'GANLSMT' // dtoChp%dateLabel // 'S.unf.' // TruncOut
    else
      DataOut = 'GANLNMC' // dtoChp%dateLabel // 'S.unf.' // TruncOut
    endif

  end subroutine fillDataOutFileName



  subroutine GetImaxJmax (Mend, Imax, Jmax, linG)
    !# Gets Imax and Jmax
    !# ---
    !# @info
    !# **Brief:** Gets Imax and Jmax. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none

    integer, intent (in) :: Mend
    integer, intent (OUT) :: Imax, Jmax
    logical, intent(in), optional :: linG

    logical :: lineargrid

    integer :: Nx, Nm, N2m, N3m, N5m, &
      n2, n3, n5, j, n, Check, Jfft

    integer, save :: Lfft = 40000

    integer, dimension (:), allocatable, save :: Ifft

    if (present(linG)) then
      linearGrid = linG
    else
      linearGrid = .false.
    end if

    N2m = ceiling(log(real(Lfft, p_r8)) / log(2.0_p_r8))
    N3m = ceiling(log(real(Lfft, p_r8)) / log(3.0_p_r8))
    N5m = ceiling(log(real(Lfft, p_r8)) / log(5.0_p_r8))
    Nx = N2m * (N3m + 1) * (N5m + 1)

    allocate (Ifft (Nx))
    Ifft = 0

    n = 0
    do n2 = 1, N2m
      Jfft = (2**n2)
      if (Jfft > Lfft) exit
      do n3 = 0, N3m
        Jfft = (2**n2) * (3**n3)
        if (Jfft > Lfft) exit
        do n5 = 0, N5m
          Jfft = (2**n2) * (3**n3) * (5**n5)
          if (Jfft > Lfft) exit
          n = n + 1
          Ifft(n) = Jfft
        enddo
      enddo
    enddo
    Nm = n

    n = 0
    do
      Check = 0
      n = n + 1
      do j = 1, Nm - 1
        if (Ifft(j) > Ifft(j + 1)) then
          Jfft = Ifft(j)
          Ifft(j) = Ifft(j + 1)
          Ifft(j + 1) = Jfft
          Check = 1
        endif
      enddo
      if (Check == 0) exit
    enddo

    if (LinearGrid) then
      Jfft = 2
    else
      Jfft = 3
    endif
    Imax = Jfft * Mend + 1
    do n = 1, Nm
      if (Ifft(n) >= Imax) then
        Imax = Ifft(n)
        exit
      endif
    enddo
    Jmax = Imax / 2
    if (mod(Jmax, 2) /= 0) Jmax = Jmax + 1

    deallocate (Ifft)

  end subroutine GetImaxJmax


  function GetGDASDateLabelRes (dtoChp)  result(isExecOk)
    use netcdf

    !# Gets GDAS Date Label Res
    !# ---
    !# @info
    !# **Brief:** Gets GDAS Date Label Res. </br>
    !# **Authors**: </br>
    !# &bull; Paulo Bonatti </br>
    !# **Date**: nov/2008 <br>
    !# @endinfo
    implicit none
    type(ChoppingNameListData), intent(inout) :: dtoChp
    logical :: isExecOk

    integer (kind = p_i8), dimension (4) :: Date

    real (kind = p_r4) :: TimODay

    character (len = 1), dimension (32) :: Descriptor
    real (kind = p_r4), dimension (2 * 100 + 1) :: SiSl
    real (kind = p_r4), dimension (44) :: Extra
    character (len=528)   ::FILE_NAME

    !GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-
    integer :: iret
    type(nemsio_gfile) :: gfile
    integer :: ios
    integer             :: ncid,ivar,zdimid
    integer             :: ndims_in
    integer             :: nvars_in
    integer             :: ngatts_in ,itime
    integer             :: unlimdimid_in
    integer, parameter :: nFile  =   1
    integer, parameter :: nMaxFile  =   1
    integer, parameter :: nvar  =   22
    INTEGER, PARAMETER :: MAX_ATT_LEN      = 800
    CHARACTER(LEN=2)   :: cDate1
    CHARACTER(LEN=2)   :: cDate2
    CHARACTER(LEN=4)   :: cDate4
    CHARACTER(LEN=2)   :: cDate3
    CHARACTER(LEN=2)   :: citime
    integer, parameter :: natrb =   8
    INTEGER            :: NX, NY,  NZ  ,NZ2  , NT
    CHARACTER(LEN=50)  :: xname, yname,zname,zname2, tname
    ! We recommend that each variable carry a "units" attribute.
    TYPE :: attribute_time
       character (len =  9) :: long_name             = "long_name"
       character (len =  9) :: units                 = "units" ;
       character (len = 14) :: cartesian_axis        = "cartesian_axis"
       character (len = 13) :: calendar_type         = "calendar_type"
       character (len =  9) :: calendar              = "calendar" ;

       character (len = 99) :: var_long_name           != "Time" ;
       character (len = 99) :: var_units               != "days since 1948-01-01 00:00:00" ;
       character (len = 99) :: var_cartesian_axis      != "T" ;
       character (len = 99) :: var_calendar_type       != "JULIAN"
       character (len = 99) :: var_calendar            != "JULIAN" ;
    END TYPE attribute_time
    TYPE(attribute_time) :: Coordtime

    CHARACTER (len = 26), PARAMETER :: ATR_NAME(1:natrb)=RESHAPE((/&
        "hydrostatic     ",&
        "ncnsto          ",&
        "ak              ",&
        "bk              ",&
        "source          ",&
        "grid            ",&
        "im              ",&
        "jm              "/),(/natrb/))
    TYPE :: GBL_attribute
      INTEGER                       :: IntegerGlobal     !"= "UNCLASSIFIED" ;
      CHARACTER (len = MAX_ATT_LEN) :: StringGlobal(natrb)         !"= "UNCLASSIFIED" ;
      REAL                          :: FloatGlobal_ak (128)      !"= "UNCLASSIFIED" ;
      REAL                          :: FloatGlobal_bk (128)      !"= "UNCLASSIFIED" ;
      INTEGER                       :: IntegerGlobal_im     !"= "UNCLASSIFIED" ;
      INTEGER                       :: IntegerGlobal_jm     !"= "UNCLASSIFIED" ;
    END TYPE GBL_attribute
    CHARACTER (len = 20), PARAMETER :: VAR_NAME(1:nvar,1:nMaxFile)=RESHAPE((/&
        "grid_xt       ","lon           ","grid_yt       ","lat           ","pfull         ","phalf         ","time          ",&
        "clwmr         ","delz          ","dpres         ","dzdt          ","grle          ","hgtsfc        ","icmr          ",&
        "o3mr          ","pressfc       ","rwmr          ","snmr          ","spfh          ","tmp           ","ugrd          ",&
        "vgrd          "/),(/nvar,nMaxFile/))
    TYPE(GBL_attribute), ALLOCATABLE :: Agbl(:)
    integer                , ALLOCATABLE :: var_varid    (:,:)

    type(nemsio_head) :: gfshead
    integer :: nfnmc  

    !GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-

    isExecOk = .false.
    if(myid.eq.0) write (unit = p_nfprt, fmt = '(/,A)') ' Getting Date Label from GDAS NCEP File'
    !GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-

    if(trim(dtoChp%DataGDAS) == 'Spec')then
      nfnmc = openFile(trim(trim(dtoChp%dgdInp) // trim(dtoChp%gDASInp)), 'unformatted', 'sequential', -1, 'read', 'old')
      if(nfnmc < 0) return
      ! Descriptor dtoChp%dateLabel (See NMC Office Note 85) (Not Used at CPTEC)
      read (unit = nfnmc) Descriptor

      ! TimODay : Time of Day in Seconds
      ! Date(1) : dtoChp%utc Hour
      ! Date(2) : Month
      ! Date(3) : Day
      ! Date(4) : Year
      ! Sigma Interfaces and Layers - SiSl(1:201)
      ! Extra Information           - Extra(1:44)
      !  o ID Sigma Structure       - Extra(18)
      !    = 1 Phillips
      !    = 2 Mean
      !  o ID Vertical Coordinate   - Extra(19)
      !    = 1 Sigma (0 for old files)
      !    = 2 Sigma-p
      read  (unit = nfnmc) TimODay, Date, SiSl, Extra
      close (unit = nfnmc)

      ! dtoChp%dateLabel='yyyymmddhh'
      write (dtoChp%dateLabel(1:4), fmt = '(I4.4)') Date(4)
      write (dtoChp%dateLabel(5:6), fmt = '(I2.2)') Date(2)
      write (dtoChp%dateLabel(7:8), fmt = '(I2.2)') Date(3)
      write (dtoChp%dateLabel(9:10), fmt = '(I2.2)') Date(1)

      ! dtoChp%utc
      write (dtoChp%utc(1:2), fmt = '(I2.2)') Date(1)

      ! Resolution
      dtoChp%mEndInp = int(Extra(1))
      dtoChp%kMaxInp = int(Extra(2))

      IDSLInp = int(Extra(18), p_i8)
      IDVCInp = int(Extra(19), p_i8)

    else if(trim(dtoChp%DataGDAS) == 'Grid') then
      ! 
      !Inicializa a lib nemsio
      call nemsio_init(iret = iret)
      ! Abertura do arquivo nemsio

      call nemsio_open(gfile, trim(dtoChp%dgdInp) // trim(dtoChp%gDASInp), 'read', mpiCommGroup, iret = iret)

      !  open (read) nemsio grid file headers

      call nemsio_getfilehead(gfile, &
        idate = gfshead%idate, nfhour = gfshead%nfhour, nfminute = gfshead%nfminute, &
        nfsecondn = gfshead%nfsecondn, nfsecondd = gfshead%nfsecondd, &
        version = gfshead%version, nrec = gfshead%nrec, dimx = gfshead%dimx, &
        dimy = gfshead%dimy, dimz = gfshead%dimz, jcap = gfshead%jcap, &
        ntrac = gfshead%ntrac, ncldt = gfshead%ncldt, nsoil = gfshead%nsoil, &
        idsl = gfshead%idsl, idvc = gfshead%idvc, idvm = gfshead%idvm, &
        idrt = gfshead%idrt, extrameta = gfshead%extrameta, &
        nmetavari = gfshead%nmetavari, nmetavarr = gfshead%nmetavarr, &
        nmetavarl = gfshead%nmetavarl, nmetavarr8 = gfshead%nmetavarr8, &
        nmetaaryi = gfshead%nmetaaryi, nmetaaryr = gfshead%nmetaaryr, &
        iret = ios)

      if(myid.eq.0) then
        print*, 'idate=', gfshead%idate
        print*, 'nfhour=', gfshead%nfhour
        print*, 'nfminute=', gfshead%nfminute
        print*, 'nfsecondn=', gfshead%nfsecondn
        print*, 'nfsecondd=', gfshead%nfsecondd
        print*, 'version=', gfshead%version
        print*, 'nrec=', gfshead%nrec
        print*, 'dimx=', gfshead%dimx
        print*, 'dimy=', gfshead%dimy
        print*, 'dimz=', gfshead%dimz
        print*, 'jcap=', gfshead%jcap
        print*, 'ntrac=', gfshead%ntrac
        print*, 'ncldt=', gfshead%ncldt
        print*, 'nsoil=', gfshead%nsoil
        print*, 'idsl=', gfshead%idsl
        print*, 'idvc=', gfshead%idvc
        print*, 'idvm=', gfshead%idvm
        print*, 'idrt=', gfshead%idrt
        print*, 'extrameta=', gfshead%extrameta
        print*, 'nmetavari=', gfshead%nmetavari
        print*, 'nmetavarr=', gfshead%nmetavarr
        print*, 'nmetavarl=', gfshead%nmetavarl
        print*, 'nmetavarr8=', gfshead%nmetavarr8
        print*, 'nmetaaryi=', gfshead%nmetaaryi
        print*, 'nmetaaryr=', gfshead%nmetaaryr
        print*, 'iret=', ios

      endif

      TimODay = gfshead%nfsecondn
      Date   (1) = gfshead%idate(4)
      Date   (3) = gfshead%idate(3)
      Date   (2) = gfshead%idate(2)
      Date   (4) = gfshead%idate(1)

      ! SiSl
      ! Extra

      ! TimODay : Time of Day in Seconds
      ! Date(1) : dtoChp%utc Hour
      ! Date(2) : Month
      ! Date(3) : Day
      ! Date(4) : Year
      ! Sigma Interfaces and Layers - SiSl(1:201)
      ! Extra Information           - Extra(1:44)
      !  o ID Sigma Structure       - Extra(18)
      !    = 1 Phillips
      !    = 2 Mean
      !  o ID Vertical Coordinate   - Extra(19)
      !    = 1 Sigma (0 for old files)
      !    = 2 Sigma-p


      !Fecha o arquivo nemsio

      call nemsio_close(gfile, iret = iret)

      !Finaliza

      call nemsio_finalize()

      ! dtoChp%dateLabel='yyyymmddhh'
      write (dtoChp%dateLabel(1:4), fmt = '(I4.4)') Date(4)
      write (dtoChp%dateLabel(5:6), fmt = '(I2.2)') Date(2)
      write (dtoChp%dateLabel(7:8), fmt = '(I2.2)') Date(3)
      write (dtoChp%dateLabel(9:10), fmt = '(I2.2)') Date(1)

      ! dtoChp%utc
      write (dtoChp%utc(1:2), fmt = '(I2.2)') Date(1)

      ! Resolution
      dtoChp%mEndInp = gfshead%jcap
      dtoChp%kMaxInp = gfshead%dimz
      ImaxInp = gfshead%dimx
      JmaxInp = gfshead%dimy
      IDSLInp = gfshead%idsl
      IDVCInp = gfshead%idvc

      ! Resolution
    ELSE IF(TRIM(dtoChp%DataGDAS) == 'Netcdf') THEN
      !
      If (.not. Allocated( Agbl      ) )   ALLOCATE(Agbl      (1:nMaxFile))
      If (.not. Allocated( var_varid ) )   ALLOCATE(var_varid (1:nvar,1:nMaxFile))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Open the input file.

      FILE_NAME=TRIM(dtoChp%dgdInp)//TRIM(dtoChp%gDASInp)

      PRINT*,TRIM(FILE_NAME)
      if(.not. check( nf90_open  (TRIM(FILE_NAME)    , nf90_nowrite, ncid) )) return

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !
      ! There are a number of inquiry functions in netCDF which can be
      ! used to learn about an unknown netCDF file. NF90_INQ tells how many
      ! netCDF variables, dimensions, and global attributes are in the
      ! file; also the dimension id of the unlimited dimension, if there
      ! is one.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(.not. check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in))) return

      !Inquire about the dimensions
      !:-------:-------:-------:-------:-------:-------:-------:
      if(.not. check(nf90_inquire_dimension(ncid,1,xname ,NX))) return
      if(.not. check(nf90_inquire_dimension(ncid,2,yname ,NY))) return
      if(.not. check(nf90_inquire_dimension(ncid,3,zname ,NZ))) return
      if(.not. check(nf90_inquire_dimension(ncid,4,zname2,NZ2))) return
      if(.not. check(nf90_inquire_dimension(ncid,5,tname ,NT))) return
      !'
      !
      ! Get the varids of the pressure and temperature netCDF variables.
      !
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(1)    , Agbl(nFile)%StringGlobal(1)))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(2)    , Agbl(nFile)%IntegerGlobal))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(3)    , Agbl(nFile)%FloatGlobal_ak))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(4)    , Agbl(nFile)%FloatGlobal_bk))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(5)    , Agbl(nFile)%StringGlobal(3)))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(6)    , Agbl(nFile)%StringGlobal(4)))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(7)    , Agbl(nFile)%IntegerGlobal_im))) return
      if(.not. check( nf90_get_att(ncid, nf90_global,  ATR_NAME(8)    , Agbl(nFile)%IntegerGlobal_jm))) return

      ! Get the varids of the latitude and longitude coordinate variables.

      DO ivar=1,nvar
         if(.not. check( nf90_inq_varid(ncid, VAR_NAME(ivar,nFile) , var_varid(ivar,nFile)) )) return
      END DO


  !
      if(.not. check( nf90_get_att(ncid    ,var_varid(7,nFile),  Coordtime%long_name            , Coordtime%var_long_name         ) )) return
      if(.not. check( nf90_get_att(ncid    ,var_varid(7,nFile),  Coordtime%units                , Coordtime%var_units             ) )) return
      if(.not. check( nf90_get_att(ncid    ,var_varid(7,nFile),  Coordtime%cartesian_axis       , Coordtime%var_cartesian_axis    ) )) return
      if(.not. check( nf90_get_att(ncid    ,var_varid(7,nFile),  Coordtime%calendar_type        , Coordtime%var_calendar_type     ) )) return
      if(.not. check( nf90_get_att(ncid    ,var_varid(7,nFile),  Coordtime%calendar             , Coordtime%var_calendar          ) )) return

      !hours since 2021-03-22 12:00:00
      cDate4 =   Coordtime%var_units(13:16)
      cDate2 =   Coordtime%var_units(18:19)
      cDate3 =   Coordtime%var_units(21:22)
      cDate1 =   Coordtime%var_units(24:25)
      citime =   Coordtime%var_units(30:31)
      READ(cDate4,'(I4.4)')Date   (4)
      READ(cDate2,'(I2.2)')Date   (2)
      READ(cDate3,'(I2.2)')Date   (3)
      READ(cDate1,'(I2.2)')Date   (1)
      READ(citime,'(I2.2)')itime
      TimODay         =itime


      !   SiSl
      !   Extra

      ! TimODay : Time of Day in Seconds
      ! Date(1) : UTC Hour
      ! Date(2) : Month
      ! Date(3) : Day
      ! Date(4) : Year
      ! Sigma Interfaces and Layers - SiSl(1:201)
      ! Extra Information	    - Extra(1:44)
      !  o ID Sigma Structure	    - Extra(18)
      !    = 1 Phillips
      !    = 2 Mean
      !  o ID Vertical Coordinate   - Extra(19)
      !    = 1 Sigma (0 for old files)
      !    = 2 Sigma-p

      if(.not. check( nf90_inq_dimid(ncid, "pfull", zdimid))) return



      ! DateLabel='yyyymmddhh'
      WRITE (dtoChp%dateLabel(1: 4), FMT='(I4.4)') Date(4)
      WRITE (dtoChp%dateLabel(5: 6), FMT='(I2.2)') Date(2)
      WRITE (dtoChp%dateLabel(7: 8), FMT='(I2.2)') Date(3)
      WRITE (dtoChp%dateLabel(9:10), FMT='(I2.2)') Date(1)

      ! UTC
      WRITE (dtoChp%utc(1:2), FMT='(I2.2)') Date(1)

      ! Resolution
      dtoChp%mEndInp=gfshead%jcap

      dtoChp%kMaxInp=NZ
      ImaxInp=NX
      JmaxInp=NY
    !  o ID Sigma Structure       - Extra(18)
    !    = 1 Phillips
    !    = 2 Mean
    !  o ID Vertical Coordinate   - Extra(19)
    !    = 1 Sigma (0 for old files)
    !    = 2 Sigma-p
    !IDSLInp=INT(Extra(18),i4)
    !IDVCInp=INT(Extra(19),i4)

      IDSLInp=0

      IDVCInp=2

      if(.not. check( nf90_close(ncid) )) return

      !Finaliza
      If ( Allocated( Agbl      ) )   DEALLOCATE(Agbl    )
      If ( Allocated( var_varid ) )   DEALLOCATE(var_varid )

    endif
    isExecOk = .true.
  end function GetGDASDateLabelRes


function check(status) result(isExecOk)
    implicit none
    logical :: isExecOk

    integer, INTENT ( in) :: status
    isExecOk = .false.
    if(status /= nf90_noerr) then
      call fatalError(headerMsg, "Error checking Netcdf")
      return
    end if
    isExecOk = .true.

  end function check


end module Mod_InputParameters
