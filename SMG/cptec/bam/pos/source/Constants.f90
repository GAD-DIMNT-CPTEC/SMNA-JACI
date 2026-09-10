!
!  $Author: pkubota $
!  $Date: 2007/10/26 17:07:02 $
!  $Revision: 1.7 $
!
MODULE Constants

   IMPLICIT NONE

   PRIVATE

                               ! Selecting Kinds:
                               ! r4 : Kind for 32-bits Real Numbers
                               ! i4 : Kind for 32-bits Integer Numbers
                               ! r8 : Kind for 64-bits Real Numbers
   INTEGER, PARAMETER, PUBLIC :: r4 = SELECTED_REAL_KIND(6)
   INTEGER, PARAMETER, PUBLIC :: i4 = SELECTED_INT_KIND(9)
   INTEGER, PARAMETER, PUBLIC :: r8 = SELECTED_REAL_KIND(15)
   INTEGER, PARAMETER, PUBLIC :: r16 = SELECTED_REAL_KIND(15)
   REAL (KIND=r8), PARAMETER, PUBLIC :: pai=3.14159265358979_r8
   REAL (KIND=r8), PARAMETER, PUBLIC :: twomg=1.458492e-4_r8

                     ! nffrs : value to indicated if model use or not initialization
                     !         or to indicate if it is a cold or warm start run:
                     !         nffrs=-1 - for runs with normal mode initialization
                     !         nffrs=0  - for runs without normal mode initialization
                     !         nffrs=1  - for warm start runs
                     ! nfbeg : number of the first forecasted file to be post-processed
                     ! nfend : number of forecasted files to be post-processed
   INTEGER, PUBLIC :: nFFrs=-1 
   INTEGER, PUBLIC :: nFBeg=1
   INTEGER, PUBLIC :: nFEnd=1

                     ! nmand : number of pressure levels listed below in format 10f8.2
                     !         nmand=-1 means the use of default 18-levels
   INTEGER, PUBLIC :: nmand=-1

                     ! RegIntIn : flag to interpolate outputs on regular grid (.TRUE.)
                     !            .FALSE. to get outputs on Gaussian grid
   INTEGER, PUBLIC :: ibdim_size = 12  ! size of basic data block (ibmax)
   INTEGER, PUBLIC :: tamBlock = 512 ! number of fft's allocated in each block
   INTEGER, PUBLIC :: nproc_vert = 1 !Number of processors to be used in the 
                                ! vertical (if givenfouriergroups set to TRUE)
   LOGICAL, PUBLIC :: givenfouriergroups =.TRUE.  ! False if processor division
                                                  ! should be automatic
   LOGICAL, PUBLIC :: RegIntIn=.FALSE.
                     !
                     ! ENS : flag for ensemble products
                     !
   LOGICAL, PUBLIC :: ENS=.FALSE.

           ! Flag to Gaussian grid type
   LOGICAL, PUBLIC :: Linear=.FALSE.

                     ! Flag to binary or grib format output
   LOGICAL, PUBLIC :: Binary=.FALSE.

   REAL  (KIND=r8), PUBLIC :: res=-1 ! horizontal resolution

   CHARACTER (LEN=6 ), PUBLIC :: trunc =' ' ! horizontal truncation = TQxxxx
   CHARACTER (LEN=4 ), PUBLIC :: lev   =' ' ! vertical layers = Lxx
   CHARACTER (LEN=10), PUBLIC :: labeli=' ' ! first run forecast label (yyyymmddhh)
   CHARACTER (LEN=10), PUBLIC :: labelf=' ' ! final run forecast label (yyyymmddhh)
   CHARACTER (LEN=3 ), PUBLIC :: prefx =' ' ! preffix for input and output files
   CHARACTER (LEN=3 ), PUBLIC :: prefy ='POS' ! preffix for output files
   CHARACTER (LEN=1 ), PUBLIC :: req   =' ' ! flag to select requested field file:
                                            ! p for rfd.pnt (wheather list)
                                            ! s for rfd.sfc (surface list)
                                            ! c for rfd.clm (climatological list)
                                            ! e for rfd.eta (etafiles list)
                                            ! nothing for rfd (default list)
   INTEGER            , PUBLIC :: kpds13 = 11  !  interval in hours to output diagnostics  
   CHARACTER (LEN=256), PUBLIC :: datain=' '  ! input directory
   CHARACTER (LEN=256), PUBLIC :: datalib=' ' ! lib directory
   CHARACTER (LEN=256), PUBLIC :: dataout=' ' ! output directory

   REAL (KIND=r8), PUBLIC :: Undef  ! Missing or Undefined Value
   REAL (KIND=r8), PUBLIC :: CvMbPa ! Converts MiliBar to Pascal (Pa/Mb)
   REAL (KIND=r8), PUBLIC :: CvLHEv ! Converts Latent Heat on Evaporation

   REAL (KIND=r8), PUBLIC :: EMRad   ! Earth Mean Radio (m)
   REAL (KIND=r8), PUBLIC :: EMRad1  ! 1 by Earth Mean Radio (1/m)
   REAL (KIND=r8), PUBLIC :: EMRad12 ! 1 by Square of Earth Mean Radio (1/m2)
   REAL (KIND=r8), PUBLIC :: EMRad2  ! Square of Earth Mean Radio (m2)

   REAL (KIND=r8), PUBLIC :: Po     ! Reference Surface Pressure Value (Mb)
   REAL (KIND=r8), PUBLIC :: P5     ! Reference Pressure Value (Mb)
   REAL (KIND=r8), PUBLIC :: Pt     ! Reference Top Pressure Value (Mb)

   REAL (KIND=r8), PUBLIC :: a    ! Constant Used in Tetens Formulae for Es
   REAL (KIND=r8), PUBLIC :: b    ! Constant Used in Tetens Formulae for Es
   REAL (KIND=r8), PUBLIC :: To   ! Constant Used in Tetens Formulae for Es (K)
   REAL (KIND=r8), PUBLIC :: Eo   ! Constant Used in Tetens Formulae for Es (Mb)

   REAL (KIND=r8), PUBLIC :: Tref  ! Temperature Reference Value (K)
   REAL (KIND=r8), PUBLIC :: Zref  ! Geopotential Height Reference Value (m)
   REAL (KIND=r8), PUBLIC :: TVVTa ! Atmosphere Lapse Rate a
   REAL (KIND=r8), PUBLIC :: TVVTb ! Atmosphere Lapse Rate b

   REAL (KIND=r8), PUBLIC :: PRHcut ! Pressure Level Cut for Specific Humidity (Mb)
   REAL (KIND=r8), PUBLIC :: RHmin  ! Minimum Value for Relative Humidity
   REAL (KIND=r8), PUBLIC :: RHmax  ! Maximum Value for Relative Humidity
   REAL (KIND=r8), PUBLIC :: SHmin  ! Minimum Value for Specific Humidity (kg/kg)

   REAL (KIND=r8), PUBLIC :: Grav   ! Earth Gravity Aceleration (m/s2)
   REAL (KIND=r8), PUBLIC :: Rd     ! Dry Air Gas Constant (m2/s2/K)
   REAL (KIND=r8), PUBLIC :: Rv     ! Water Vapor Gas Constant (m2/s2/K)
   REAL (KIND=r8), PUBLIC :: Cp     ! Dry Air Gas Specific Heat Capcity at Pressure Constant (m2/s2/K)
   REAL (KIND=r8), PUBLIC :: Cv     ! Dry Air Gas Specific Heat Capcity at Volumn Constant (m2/s2/K)

   REAL (KIND=r8), PUBLIC :: Eps  ! Rd by Rv
   REAL (KIND=r8), PUBLIC :: Eps1 ! 1-Eps
   REAL (KIND=r8), PUBLIC :: CTv  ! Constant Used to Convert Tv into T, or vice-versa.

   REAL (KIND=r8), PUBLIC :: RdByCp   ! Rd by Cp
   REAL (KIND=r8), PUBLIC :: RdByCp1  ! Rd by Cp + 1
   REAL (KIND=r8), PUBLIC :: CpByRd   ! Cp by Rd
   REAL (KIND=r8), PUBLIC :: RdByGrav ! Rd by Grav
   REAL (KIND=r8), PUBLIC :: GravByRd ! -Grav by Rd

   INTEGER, PARAMETER, PUBLIC :: ndv=58 !pkubota
   INTEGER, PARAMETER, PUBLIC :: mdf=3

   INTEGER, PUBLIC :: kdv(ndv)
   INTEGER, PUBLIC :: nvv(ndv)
   INTEGER, PUBLIC :: mkdv(ndv)
   INTEGER, PUBLIC :: lif(ndv)
   INTEGER, PUBLIC :: iclcd(ndv)
   INTEGER, PUBLIC :: nudv(ndv)
   INTEGER, PUBLIC :: nureq(ndv,mdf)
   CHARACTER (LEN=7  ),PUBLIC :: rfd   =' ' ! post processed file name of required fields

   CHARACTER (LEN=40), PUBLIC :: chrdv(ndv)
   CHARACTER (LEN=40), PUBLIC :: chreq(ndv,mdf)
   LOGICAL, PUBLIC :: postclim
   LOGICAL, PUBLIC :: ExtrapoAdiabatica=.FALSE.
   ! Files Units

   INTEGER, PUBLIC :: nferr=0    ! Standard Error Print Out
   INTEGER, PUBLIC :: nfprt=6    ! Standard Print Out
   INTEGER, PUBLIC :: nfinp=7    ! Read In Namelist
   INTEGER, PUBLIC :: nffct=10   ! To Read Forecasting Output Data
   INTEGER, PUBLIC :: nfpos=11   ! To Write Post-Processed Data
   INTEGER, PUBLIC :: nfctl=12   ! To Write Post-Processed Description
   iNTEGER, PUBLIC :: nfmdf=13   ! To Read Model Output Filenames
   INTEGER, PUBLIC :: nfppf=14   ! To Write Post-Processing Output Filenames
   INTEGER, PUBLIC :: nfdir=31   ! To Read Forecast Output File Structure
   INTEGER, PUBLIC :: nfrfd=32   ! To Read Requested Post-processed Fields

   INTEGER, PUBLIC :: newlat0
   INTEGER, PUBLIC :: newlat1
   INTEGER, PUBLIC :: newlon0
   INTEGER, PUBLIC :: newlon1
   REAL (KIND=r8), PUBLIC :: RecLat(2)=(/-90.0_r8, 90.0_r8/) !
   REAL (KIND=r8), PUBLIC :: RecLon(2)=(/  0.0_r8,360.0_r8/) ! 
   LOGICAL, PUBLIC :: RunRecort=.FALSE.
   REAL (KIND=r8), PUBLIC :: plevs(160) ! 

  PUBLIC :: InitParameters

CONTAINS

SUBROUTINE InitParameters ()

   IMPLICIT NONE

   INTEGER :: ierr
   NAMELIST /PosInput/ nffrs, nfbeg, nfend, nmand, RegIntIn, Linear, &
                       trunc, lev, labeli, labelf,kpds13, prefx,prefy, req,  &
                       datain,datalib,dataout,&
                       Binary,postclim,res,ENS,ExtrapoAdiabatica,RecLat,RecLon,RunRecort,&
                       givenfouriergroups,nproc_vert,ibdim_size,tamBlock

   NAMELIST /PressureLevel/plevs 
   OPEN(unit=nfinp, file="POSTIN-GRIB", action="read", status="old", iostat=ierr)

   READ (nfinp, NML=PosInput, iostat=ierr)
   print *, "trunc=", trunc

   plevs=0.0_r8
   READ (nfinp, NML=PressureLevel, iostat=ierr)

   CALL InitConstants ()

   CALL InitDFT ()

END SUBROUTINE InitParameters


SUBROUTINE InitConstants ()

   IMPLICIT NONE

   REAL (KIND=r8) :: MWWater   ! Molecular Weight of Water (kg/kmol)
   REAL (KIND=r8) :: MWN2      ! Molecular Weight of Nitrogen (kg/kmol)
   REAL (KIND=r8) :: PN2       ! Atmosphere Nitrogen Percetage
   REAL (KIND=r8) :: MWO2      ! Molecular Weight of Oxigen (kg/kmol)
   REAL (KIND=r8) :: PO2       ! Atmosphere Oxigen Percetage
   REAL (KIND=r8) :: MWAr      ! Molecular Weight of Argon (kg/kmol)
   REAL (KIND=r8) :: PAr       ! Atmosphere Argon Percetage
   REAL (KIND=r8) :: MWCO2     ! Molecular Weight of Carbon Gas (kg/kmol)
   REAL (KIND=r8) :: PCO2      ! Atmosphere Carbon Gas Percetage
   REAL (KIND=r8) :: MWDAir    ! Mean Molecular Weight of Dry Air (kg/kmol)
   REAL (KIND=r8) :: Avogrado  ! Avogrado's Constant (10**-20/kmol)
   REAL (KIND=r8) :: Boltzmann ! Boltzmann's Constant (J/K)
   REAL (KIND=r8) :: REstar    ! Universal Gas Constant (10**23 J/K/kmol)

   Undef=-2.56E+33_r8

   CvMbPa=100.0_r8
   CvLHEv=28.9_r8

   EMRad=6.37E6_r8
   EMRad1=1.0_r8/EMRad
   EMRad12=EMRad1*EMRad1
   EMRad2=EMRad*EMRad

   Po=1000.0_r8
   P5=500.0_r8
   Pt=0.1_r8

   a=17.2693882_r8
   b=35.86_r8
   To=273.16_r8
   Eo=6.1078_r8

   Tref=290.66_r8
   Zref=75.0_r8
   TVVTa=0.0065_r8
   TVVTb=0.0050_r8

   PRHcut=2.0_r8
   RHmin=0.00000001_r8
   RHmax=0.99999999_r8
   SHmin=1.0e-7_r8

   MWWater=18.016_r8
   MWN2=28.0134_r8
   PN2=0.7809_r8
   MWO2=31.9988_r8
   PO2=0.2095_r8
   MWAr=39.9480_r8
   PAr=0.0093_r8
   MWCO2=44.0103_r8
   PCO2=0.0003_r8
   MWDAir=PN2*MWN2+PO2*MWO2+PAr*MWAr+PCO2*MWCO2

   Avogrado=6022.52_r8
   Boltzmann=1.38054_r8
   REstar=Avogrado*Boltzmann

   Grav=9.80665_r8
   Rd=REstar/MWDAir
   Rv=REstar/MWWater
   Cp=3.5_r8*Rd
   Cv=Cp-Rd

   Eps=Rd/Rv
   Eps1=1.0_r8-Eps
   CTv=Eps1/Eps

   RdByCp=Rd/Cp
   RdByCp1=RdByCp+1.0_r8
   CpByRd=Cp/Rd
   RdByGrav=Rd/Grav
   GravByRd=-Grav/Rd

   IF(Linear)THEN
       IF(trunc(2:2) /= 'L' ) THEN
          WRITE(nfprt,*)'error in trunc=',TRIM(trunc)
          STOP
       END IF
   ELSE
       IF(trunc(2:2) /= 'Q' ) THEN
          WRITE(nfprt,*)'error in trunc=',TRIM(trunc)
          STOP
       END IF
   END IF 


END SUBROUTINE InitConstants

SUBROUTINE InitDFT ()

   IMPLICIT NONE

   chrdv(1:30)= (/ &
        'SURFACE PRESSURE                        ', &
        'MASK                                    ', &
        'SURFACE ZONAL WIND (U)                  ', &
        'ZONAL WIND (U)                          ', &
        'SURFACE MERIDIONAL WIND (V)             ', &
        'MERIDIONAL WIND (V)                     ', &
        'OMEGA                                   ', &
        'DIVERGENCE                              ', &
        'VORTICITY                               ', &
        'STREAM FUNCTION                         ', &
        'ZONAL WIND PSI                          ', &
        'MERIDIONAL WIND PSI                     ', &
        'VELOCITY POTENTIAL                      ', &
        'ZONAL WIND CHI                          ', &
        'MERIDIONAL WIND CHI                     ', &
        'GEOPOTENTIAL HEIGHT                     ', &
        'SEA LEVEL PRESSURE                      ', &
        'SURFACE ABSOLUTE TEMPERATURE            ', &
        'ABSOLUTE TEMPERATURE                    ', &
        'SURFACE RELATIVE HUMIDITY               ', &
        'RELATIVE HUMIDITY                       ', &
        'SPECIFIC HUMIDITY                       ', &
        'INST. PRECIP. WATER                     ', &
        'POTENTIAL TEMPERATURE                   ', &
        'SURFACE TEMPERATURE                     ', &
        'TIME MEAN SURFACE PRESSURE              ', &
        'TIME MEAN MASK                          ', &
        'TIME MEAN SURFACE ZONAL WIND (U)        ', &
        'TIME MEAN ZONAL WIND (U)                ', &
        'TIME MEAN SURFACE MERIDIONAL WIND (V)   ' /)
   chrdv(31:ndv)= (/ &
        'TIME MEAN MERIDIONAL WIND (V)           ', &
        'TIME MEAN DERIVED OMEGA                 ', &
        'TIME MEAN DIVERGENCE                    ', &
        'TIME MEAN VORTICITY                     ', &
        'TIME MEAN STREAM FUNCTION               ', &
        'TIME MEAN ZONAL WIND PSI                ', &
        'TIME MEAN MERIDIONAL WIND PSI           ', &
        'TIME MEAN VELOCITY POTENTIAL            ', &
        'TIME MEAN ZONAL WIND CHI                ', &
        'TIME MEAN MERIDIONAL WIND CHI           ', &
        'TIME MEAN GEOPOTENTIAL HEIGHT           ', &
        'TIME MEAN SEA LEVEL PRESSURE            ', &
        'TIME MEAN SURFACE ABSOLUTE TEMPERATURE  ', &
        'TIME MEAN ABSOLUTE TEMPERATURE          ', &
        'TIME MEAN SURFACE RELATIVE HUMIDITY     ', &
        'TIME MEAN RELATIVE HUMIDITY             ', &
        'TIME MEAN SPECIFIC HUMIDITY             ', &
        'TIME MEAN PRECIP. WATER                 ', &
        'TIME MEAN POTENTIAL TEMPERATURE         ', &
        'TIME MEAN SURFACE TEMPERATURE           ', &
        'EVAPORATION                             ', &
        'VERT. INT. HORIZ. MOIS. DIFFUSION       ', &
        'GEOMETRIC MEAN SURFACE PRESSURE         ', &
        'VERT. INT. HORIZ. MOIS. DIFFUSION       ', &
        'TIME MEAN DEEP SOIL TEMPERATURE         ', &
        'GROUND/SURFACE COVER TEMPERATURE        ', &
        'CANOPY TEMPERATURE                      ', &
        '2 METRE DEWPOINT TEMPERATURE            '/) !36 !PK---

   kdv(1:ndv)= (/ &
          1,   1,   3,   3,   3,   3,   3,   3,   3,   2, &
          2,   2,   2,   2,   2,   2,   2,   3,   3,   3, &
          3,   3,   2,   3,   1,   1,   1,   3,   3,   3, &
          3,   3,   3,   3,   2,   2,   2,   2,   2,   2, &
          2,   2,   3,   3,   3,   3,   3,   2,   3,   1, &
          1,   2,   1,   2,   1,   1,   1,   1 /)!PK

   nvv(1:ndv)= (/ &
          1,   2,   1,   2,   1,   2,   2,   2,   2,   2, &
          2,   2,   2,   2,   2,   2,   1,   1,   2,   1, &
          2,   2,   1,   2,   1,   1,   2,   1,   2,   1, &
          2,   2,   2,   2,   2,   2,   2,   2,   2,   2, &
          2,   1,   1,   2,   1,   2,   2,   1,   2,   1, &
          1,   1,   1,   1,   1,   1,   1,   1 /)!PK

   iclcd(1:ndv)= (/ &
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
          0,   0,   0,   0,   1,   0,   0,   0,   0,   0, &
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
          0,   0,   0,   0,   0,   0,   0,   0,   0,   1, &
          2,   3,   4,   3,   1,   1,   1,   1 /)!PK

   nudv(1:ndv)= (/ &
        132,   0,  60,  60,  60,  60, 151,  50,  50,  90, &
         60,  60,  90,  60,  60,  10, 134,  40,  40,   0, &
          0,   0, 110,  40,  40, 132,   0,  60,  60,  60, &
         60, 151,  50,  50,  90,  60,  60,  90,  60,  60, &
         10, 134,  40,  40,   0,   0,   0, 110,  40,  40, &
        121, 121, 132, 121,  40,  40,  40,  40 /)!PK

   chreq( 1: 1,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq( 2: 2,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq( 3: 3,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq( 4: 4,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq( 5: 5,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq( 6: 6,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq( 7: 7,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq( 8: 8,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq( 9: 9,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        'VORTICITY                                 ' /),(/1,mdf/))
   chreq(10:10,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'VORTICITY                                 ', &
        '                                          ' /),(/1,mdf/))
   chreq(11:11,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'VORTICITY                                 ', &
        '                                          ' /),(/1,mdf/))
   chreq(12:12,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'VORTICITY                                 ', &
        '                                          ' /),(/1,mdf/))
   chreq(13:13,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        '                                          ' /),(/1,mdf/))
   chreq(14:14,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        '                                          ' /),(/1,mdf/))
   chreq(15:15,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'DIVERGENCE                                ', &
        '                                          ' /),(/1,mdf/))
   chreq(16:16,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'VIRTUAL TEMPERATURE                       ', &
        '                                          ' /),(/1,mdf/))
   chreq(17:17,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'VIRTUAL TEMPERATURE                       ', &
        '                                          ' /),(/1,mdf/))
   chreq(18:18,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        'VIRTUAL TEMPERATURE                       ' /),(/1,mdf/))
   chreq(19:19,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        'VIRTUAL TEMPERATURE                       ' /),(/1,mdf/))
   chreq(20:20,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        'VIRTUAL TEMPERATURE                       ' /),(/1,mdf/))
   chreq(21:21,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        'VIRTUAL TEMPERATURE                       ' /),(/1,mdf/))
   chreq(22:22,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        'VIRTUAL TEMPERATURE                       ' /),(/1,mdf/))
   chreq(23:23,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        '                                          ' /),(/1,mdf/))
   chreq(24:24,1:mdf)= RESHAPE((/ &
        'LN SURFACE PRESSURE                       ', &
        'SPECIFIC HUMIDITY                         ', &
        'VIRTUAL TEMPERATURE                       ' /),(/1,mdf/))
   chreq(25:25,1:mdf)= RESHAPE((/ &
        'SURFACE TEMPERATURE                       ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(26:26,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(27:27,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(28:28,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(29:29,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(30:30,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(31:31,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(32:32,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(33:33,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(34:34,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        'TIME MEAN VORTICITY                       ' /),(/1,mdf/))
   chreq(35:35,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN VORTICITY                       ', &
        '                                          ' /),(/1,mdf/))
   chreq(36:36,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN VORTICITY                       ', &
        '                                          ' /),(/1,mdf/))
   chreq(37:37,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN VORTICITY                       ', &
        '                                          ' /),(/1,mdf/))
   chreq(38:38,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        '                                          ' /),(/1,mdf/))
   chreq(39:39,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        '                                          ' /),(/1,mdf/))
   chreq(40:40,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN DIVERGENCE                      ', &
        '                                          ' /),(/1,mdf/))
   chreq(41:41,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ', &
        '                                          ' /),(/1,mdf/))
   chreq(42:42,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ', &
        '                                          ' /),(/1,mdf/))
   chreq(43:43,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ' /),(/1,mdf/))
   chreq(44:44,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ' /),(/1,mdf/))
   chreq(45:45,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ' /),(/1,mdf/))
   chreq(46:46,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ' /),(/1,mdf/))
   chreq(47:47,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ' /),(/1,mdf/))
   chreq(48:48,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        '                                          ' /),(/1,mdf/))
   chreq(49:49,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'TIME MEAN SPECIFIC HUMIDITY               ', &
        'TIME MEAN VIRTUAL TEMPERATURE             ' /),(/1,mdf/))
   chreq(50:50,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE TEMPERATURE             ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(51:51,1:mdf)= RESHAPE((/ &
        'LATENT HEAT FLUX FROM SURFACE             ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(52:52,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'HORIZONTAL MOISTURE DIFFUSION             ', &
        '                                          ' /),(/1,mdf/))
   chreq(53:53,1:mdf)= RESHAPE((/ &
        'TIME MEAN LOG SURFACE PRESSURE            ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(54:54,1:mdf)= RESHAPE((/ &
        'TIME MEAN SURFACE PRESSURE                ', &
        'HORIZONTAL MOISTURE DIFFUSION             ', &
        '                                          ' /),(/1,mdf/))
   chreq(55:55,1:mdf)= RESHAPE((/ &
        'TIME MEAN DEEP SOIL TEMPERATURE           ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(56:56,1:mdf)= RESHAPE((/ &
        'GROUND/SURFACE COVER TEMPERATURE          ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))
   chreq(57:57,1:mdf)= RESHAPE((/ &
        'CANOPY TEMPERATURE                        ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))

   chreq(58:58,1:mdf)= RESHAPE((/ &
        'TEMPERATURE AT 2-M FROM SURFACE           ', &
        '                                          ', &
        '                                          ' /),(/1,mdf/))

   nureq(1:ndv,1:mdf)= RESHAPE((/ &
         142, 142, 142, 142, 142, 142, 142, 142, 142, 142, &
         142, 142, 142, 142, 142, 142, 142, 142, 142, 142, &
         142, 142, 142, 142,  40, 132, 132, 132, 132, 132, &
         132, 132, 132, 132, 132, 132, 132, 132, 132, 132, &
         132, 132, 132, 132, 132, 132, 132, 132, 132,  40, &
         170, 132, 142, 132,  40,  40,  40,  40,                &
           0,   0,  50,  50,  50,  50,  50,  50,  50,  50, &
          50,  50,  50,  50,  50,  40,  40,   0,   0,   0, &
           0,   0,   0,   0,   0,   0,   0,  50,  50,  50, &
          50,  50,  50,  50,  50,  50,  50,  50,  50,  50, &
          40,  40,   0,   0,   0,   0,   0,   0,   0,   0, &
           0,  51,   0,  52,   0,   0,   0,   0,                &
          0,   0,   50,  50,  50,  50,  50,  50,  50,   0, &
          0,   0,    0,   0,   0,   0,   0,  40,  40,  40, &
         40,  40,    0,  40,   0,   0,   0,  50,  50,  50, &
         50,  50,   50,  50,   0,   0,   0,   0,   0,   0, &
          0,   0,   40,  40,  40,  40,  40,   0,  40,   0, &
          0,   0,    0,   0,   0,   0,   0,   0            /), &
          (/ndv,mdf/))

END SUBROUTINE InitDFT


END MODULE Constants
