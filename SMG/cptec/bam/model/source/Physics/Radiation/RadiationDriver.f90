!===============================================================================
! SVN $Id: shr_orb_mod.F90 6752 2007-10-04 21:02:15Z jwolfe $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/branch_tags/
!cesm1_0_rel_tags/cesm1_0_rel01_share3_100616/shr/shr_orb_mod.F90 $
!===============================================================================
! Modified by Tarassova 2015
!  Climate aerosol (Kinne, 2013) coarse mode included
!  Fine aerosol mode 2000 is included
!  Modifications (21) are marked by 
!  !tar begin   and  !tar end
!
MODULE shr_orb_mod

    IMPLICIT NONE
  SAVE

  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------

  INTEGER,PARAMETER :: SHR_KIND_R8 = SELECTED_REAL_KIND(15) ! 8 byte real
  INTEGER,PARAMETER :: SHR_KIND_R4 = SELECTED_REAL_KIND( 6) ! 4 byte real
  INTEGER,PARAMETER :: SHR_KIND_RN = KIND(1.0)              ! native real
  INTEGER,PARAMETER :: SHR_KIND_I8 = SELECTED_INT_KIND (14) ! 8 byte integer
  INTEGER,PARAMETER :: SHR_KIND_I4 = SELECTED_INT_KIND ( 9) ! 4 byte integer
  INTEGER,PARAMETER :: SHR_KIND_IN = KIND(1)                ! native integer
  INTEGER,PARAMETER :: SHR_KIND_CS = 80                     ! short char
  INTEGER,PARAMETER :: SHR_KIND_CL = 256                    ! long char
  INTEGER,PARAMETER :: SHR_KIND_CX = 512                    ! extra-long char
  INTEGER(SHR_KIND_IN),PARAMETER,PRIVATE :: R8 = SHR_KIND_R8 ! rename for local readability only
  !----------------------------------------------------------------------------
  ! physical constants (all data public)
  !----------------------------------------------------------------------------

  REAL(R8),PARAMETER :: SHR_CONST_PI      = 3.14159265358979323846_R8  ! pi

  !----------------------------------------------------------------------------
  ! PUBLIC: Interfaces and global data
  !----------------------------------------------------------------------------
  PUBLIC :: shr_orb_cosz
  PUBLIC :: shr_orb_params
  PUBLIC :: shr_orb_decl
  PUBLIC :: shr_orb_print

  REAL   (SHR_KIND_R8),PUBLIC,PARAMETER :: SHR_ORB_UNDEF_REAL = 1.e36_SHR_KIND_R8 ! undefined real 
  INTEGER(SHR_KIND_IN),PUBLIC,PARAMETER :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int
  ! low-level shared variables for logging, these may not be parameters
  !----------------------------------------------------------------------------


  INTEGER(SHR_KIND_IN) :: s_loglev = 0
  INTEGER(SHR_KIND_IN) :: s_logunit  = 6


  !----------------------------------------------------------------------------
  ! PRIVATE: by default everything else is private to this module
  !----------------------------------------------------------------------------

  REAL   (SHR_KIND_R8),PARAMETER :: pi                 = SHR_CONST_PI
  REAL   (SHR_KIND_R8),PARAMETER :: SHR_ORB_ECCEN_MIN  =   0.0_SHR_KIND_R8 ! min value for eccen
  REAL   (SHR_KIND_R8),PARAMETER :: SHR_ORB_ECCEN_MAX  =   0.1_SHR_KIND_R8 ! max value for eccen
  REAL   (SHR_KIND_R8),PARAMETER :: SHR_ORB_OBLIQ_MIN  = -90.0_SHR_KIND_R8 ! min value for obliq
  REAL   (SHR_KIND_R8),PARAMETER :: SHR_ORB_OBLIQ_MAX  = +90.0_SHR_KIND_R8 ! max value for obliq
  REAL   (SHR_KIND_R8),PARAMETER :: SHR_ORB_MVELP_MIN  =   0.0_SHR_KIND_R8 ! min value for mvelp
  REAL   (SHR_KIND_R8),PARAMETER :: SHR_ORB_MVELP_MAX  = 360.0_SHR_KIND_R8 ! max value for mvelp


  LOGICAL, PARAMETER                :: log_print=.FALSE. ! Flags print of status/error

  REAL   (SHR_KIND_R8) :: eccen     ! orbital eccentricity
  REAL   (SHR_KIND_R8) :: obliq      ! obliquity in degrees
  REAL   (SHR_KIND_R8) :: mvelp      ! moving vernal equinox long
  REAL   (SHR_KIND_R8) :: obliqr     ! Earths obliquity in rad
  REAL   (SHR_KIND_R8) :: lambm0     ! Mean long of perihelion at
  REAL   (SHR_KIND_R8) :: mvelpp     ! moving vernal equinox long
  REAL   (SHR_KIND_R8) :: eccf       ! Earth-sun distance factor
  REAL   (SHR_KIND_R8) :: declin     ! Solar declination (radians)

  !===============================================================================
CONTAINS
  !===============================================================================

  REAL(SHR_KIND_R8) FUNCTION shr_orb_cosz(jday,lat,lon)

    !----------------------------------------------------------------------------
    !
    ! FUNCTION to return the cosine of the solar zenith angle.
    ! Assumes 365.0 days/year.
    !
    !--------------- Code History -----------------------------------------------
    !
    ! Original Author: Brian Kauffman
    ! Date:            Jan/98
    ! History:         adapted from statement FUNCTION in share/orb_cosz.h
    !
    !----------------------------------------------------------------------------

    REAL   (SHR_KIND_R8),INTENT(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    REAL   (SHR_KIND_R8),INTENT(in) :: lat    ! Centered latitude (radians)
    REAL   (SHR_KIND_R8),INTENT(in) :: lon    ! Centered longitude (radians)

    !----------------------------------------------------------------------------

    shr_orb_cosz = SIN(lat)*SIN(declin) - &
         &              COS(lat)*COS(declin)*COS(jday*2.0_SHR_KIND_R8*pi + lon)

  END FUNCTION shr_orb_cosz

  !===============================================================================

  SUBROUTINE shr_orb_params( &
      iyear_AD )
    !-------------------------------------------------------------------------------
    !
    ! Calculate earths orbital parameters using Dave Threshers formula which 
    ! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term 
    ! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy 
    ! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
    !
    !------------------------------Code history-------------------------------------
    !
    ! Original Author: Erik Kluzek
    ! Date:            Oct/97
    !
    !-------------------------------------------------------------------------------

    !----------------------------- Arguments ------------------------------------
    INTEGER(SHR_KIND_IN),INTENT(in)    :: iyear_AD  ! Year to calculate orbit for

    !------------------------------ Parameters ----------------------------------
    INTEGER(SHR_KIND_IN),PARAMETER :: poblen =47 ! # of elements in series wrt obliquity
    INTEGER(SHR_KIND_IN),PARAMETER :: pecclen=19 ! # of elements in series wrt eccentricity
    INTEGER(SHR_KIND_IN),PARAMETER :: pmvelen=78 ! # of elements in series wrt vernal equinox
    REAL   (SHR_KIND_R8),PARAMETER :: psecdeg = 1.0_SHR_KIND_R8/3600.0_SHR_KIND_R8 ! arc sec to deg conversion

    REAL   (SHR_KIND_R8) :: degrad = pi/180._SHR_KIND_R8   ! degree to radian conversion factor
    REAL   (SHR_KIND_R8) :: yb4_1950AD         ! number of years before 1950 AD

    ! Cosine series data for computation of obliquity: amplitude (arc seconds),
    ! rate (arc seconds/year), phase (degrees).

    REAL   (SHR_KIND_R8), PARAMETER :: obamp(poblen) =  & ! amplitudes for obliquity cos series
         &      (/   -2462.2214466_SHR_KIND_R8, -857.3232075_SHR_KIND_R8, -629.3231835_SHR_KIND_R8,   &
         &            -414.2804924_SHR_KIND_R8, -311.7632587_SHR_KIND_R8,  308.9408604_SHR_KIND_R8,   &
         &            -162.5533601_SHR_KIND_R8, -116.1077911_SHR_KIND_R8,  101.1189923_SHR_KIND_R8,   &
         &             -67.6856209_SHR_KIND_R8,   24.9079067_SHR_KIND_R8,   22.5811241_SHR_KIND_R8,   &
         &             -21.1648355_SHR_KIND_R8,  -15.6549876_SHR_KIND_R8,   15.3936813_SHR_KIND_R8,   &
         &              14.6660938_SHR_KIND_R8,  -11.7273029_SHR_KIND_R8,   10.2742696_SHR_KIND_R8,   &
         &               6.4914588_SHR_KIND_R8,    5.8539148_SHR_KIND_R8,   -5.4872205_SHR_KIND_R8,   &
         &              -5.4290191_SHR_KIND_R8,    5.1609570_SHR_KIND_R8,    5.0786314_SHR_KIND_R8,   &
         &              -4.0735782_SHR_KIND_R8,    3.7227167_SHR_KIND_R8,    3.3971932_SHR_KIND_R8,   &
         &              -2.8347004_SHR_KIND_R8,   -2.6550721_SHR_KIND_R8,   -2.5717867_SHR_KIND_R8,   &
         &              -2.4712188_SHR_KIND_R8,    2.4625410_SHR_KIND_R8,    2.2464112_SHR_KIND_R8,   &
         &              -2.0755511_SHR_KIND_R8,   -1.9713669_SHR_KIND_R8,   -1.8813061_SHR_KIND_R8,   &
         &              -1.8468785_SHR_KIND_R8,    1.8186742_SHR_KIND_R8,    1.7601888_SHR_KIND_R8,   &
         &              -1.5428851_SHR_KIND_R8,    1.4738838_SHR_KIND_R8,   -1.4593669_SHR_KIND_R8,   &
         &               1.4192259_SHR_KIND_R8,   -1.1818980_SHR_KIND_R8,    1.1756474_SHR_KIND_R8,   &
         &              -1.1316126_SHR_KIND_R8,    1.0896928_SHR_KIND_R8/)

    REAL   (SHR_KIND_R8), PARAMETER :: obrate(poblen) = & ! rates for obliquity cosine series
         &        (/  31.609974_SHR_KIND_R8, 32.620504_SHR_KIND_R8, 24.172203_SHR_KIND_R8,   &
         &            31.983787_SHR_KIND_R8, 44.828336_SHR_KIND_R8, 30.973257_SHR_KIND_R8,   &
         &            43.668246_SHR_KIND_R8, 32.246691_SHR_KIND_R8, 30.599444_SHR_KIND_R8,   &
         &            42.681324_SHR_KIND_R8, 43.836462_SHR_KIND_R8, 47.439436_SHR_KIND_R8,   &
         &            63.219948_SHR_KIND_R8, 64.230478_SHR_KIND_R8,  1.010530_SHR_KIND_R8,   &
         &             7.437771_SHR_KIND_R8, 55.782177_SHR_KIND_R8,  0.373813_SHR_KIND_R8,   &
         &            13.218362_SHR_KIND_R8, 62.583231_SHR_KIND_R8, 63.593761_SHR_KIND_R8,   &
         &            76.438310_SHR_KIND_R8, 45.815258_SHR_KIND_R8,  8.448301_SHR_KIND_R8,   &
         &            56.792707_SHR_KIND_R8, 49.747842_SHR_KIND_R8, 12.058272_SHR_KIND_R8,   &
         &            75.278220_SHR_KIND_R8, 65.241008_SHR_KIND_R8, 64.604291_SHR_KIND_R8,   &
         &             1.647247_SHR_KIND_R8,  7.811584_SHR_KIND_R8, 12.207832_SHR_KIND_R8,   &
         &            63.856665_SHR_KIND_R8, 56.155990_SHR_KIND_R8, 77.448840_SHR_KIND_R8,   &
         &             6.801054_SHR_KIND_R8, 62.209418_SHR_KIND_R8, 20.656133_SHR_KIND_R8,   &
         &            48.344406_SHR_KIND_R8, 55.145460_SHR_KIND_R8, 69.000539_SHR_KIND_R8,   &
         &            11.071350_SHR_KIND_R8, 74.291298_SHR_KIND_R8, 11.047742_SHR_KIND_R8,   &
         &             0.636717_SHR_KIND_R8, 12.844549_SHR_KIND_R8/)

    REAL   (SHR_KIND_R8), PARAMETER :: obphas(poblen) = & ! phases for obliquity cosine series
         &      (/    251.9025_SHR_KIND_R8, 280.8325_SHR_KIND_R8, 128.3057_SHR_KIND_R8,   &
         &            292.7252_SHR_KIND_R8,  15.3747_SHR_KIND_R8, 263.7951_SHR_KIND_R8,   &
         &            308.4258_SHR_KIND_R8, 240.0099_SHR_KIND_R8, 222.9725_SHR_KIND_R8,   &
         &            268.7809_SHR_KIND_R8, 316.7998_SHR_KIND_R8, 319.6024_SHR_KIND_R8,   &
         &            143.8050_SHR_KIND_R8, 172.7351_SHR_KIND_R8,  28.9300_SHR_KIND_R8,   &
         &            123.5968_SHR_KIND_R8,  20.2082_SHR_KIND_R8,  40.8226_SHR_KIND_R8,   &
         &            123.4722_SHR_KIND_R8, 155.6977_SHR_KIND_R8, 184.6277_SHR_KIND_R8,   &
         &            267.2772_SHR_KIND_R8,  55.0196_SHR_KIND_R8, 152.5268_SHR_KIND_R8,   &
         &             49.1382_SHR_KIND_R8, 204.6609_SHR_KIND_R8,  56.5233_SHR_KIND_R8,   &
         &            200.3284_SHR_KIND_R8, 201.6651_SHR_KIND_R8, 213.5577_SHR_KIND_R8,   &
         &             17.0374_SHR_KIND_R8, 164.4194_SHR_KIND_R8,  94.5422_SHR_KIND_R8,   &
         &            131.9124_SHR_KIND_R8,  61.0309_SHR_KIND_R8, 296.2073_SHR_KIND_R8,   &
         &            135.4894_SHR_KIND_R8, 114.8750_SHR_KIND_R8, 247.0691_SHR_KIND_R8,   &
         &            256.6114_SHR_KIND_R8,  32.1008_SHR_KIND_R8, 143.6804_SHR_KIND_R8,   &
         &             16.8784_SHR_KIND_R8, 160.6835_SHR_KIND_R8,  27.5932_SHR_KIND_R8,   &
         &            348.1074_SHR_KIND_R8,  82.6496_SHR_KIND_R8/)

    ! Cosine/sine series data for computation of eccentricity and fixed vernal 
    ! equinox longitude of perihelion (fvelp): amplitude, 
    ! rate (arc seconds/year), phase (degrees).

    REAL   (SHR_KIND_R8), PARAMETER :: ecamp (pecclen) = & ! ampl for eccen/fvelp cos/sin series
         &      (/   0.01860798_SHR_KIND_R8,  0.01627522_SHR_KIND_R8, -0.01300660_SHR_KIND_R8,   &
         &           0.00988829_SHR_KIND_R8, -0.00336700_SHR_KIND_R8,  0.00333077_SHR_KIND_R8,   &
         &          -0.00235400_SHR_KIND_R8,  0.00140015_SHR_KIND_R8,  0.00100700_SHR_KIND_R8,   &
         &           0.00085700_SHR_KIND_R8,  0.00064990_SHR_KIND_R8,  0.00059900_SHR_KIND_R8,   &
         &           0.00037800_SHR_KIND_R8, -0.00033700_SHR_KIND_R8,  0.00027600_SHR_KIND_R8,   &
         &           0.00018200_SHR_KIND_R8, -0.00017400_SHR_KIND_R8, -0.00012400_SHR_KIND_R8,   &
         &           0.00001250_SHR_KIND_R8/)

    REAL   (SHR_KIND_R8), PARAMETER :: ecrate(pecclen) = & ! rates for eccen/fvelp cos/sin series
         &      (/    4.2072050_SHR_KIND_R8,  7.3460910_SHR_KIND_R8, 17.8572630_SHR_KIND_R8,  &
         &           17.2205460_SHR_KIND_R8, 16.8467330_SHR_KIND_R8,  5.1990790_SHR_KIND_R8,  &
         &           18.2310760_SHR_KIND_R8, 26.2167580_SHR_KIND_R8,  6.3591690_SHR_KIND_R8,  &
         &           16.2100160_SHR_KIND_R8,  3.0651810_SHR_KIND_R8, 16.5838290_SHR_KIND_R8,  &
         &           18.4939800_SHR_KIND_R8,  6.1909530_SHR_KIND_R8, 18.8677930_SHR_KIND_R8,  &
         &           17.4255670_SHR_KIND_R8,  6.1860010_SHR_KIND_R8, 18.4174410_SHR_KIND_R8,  &
         &            0.6678630_SHR_KIND_R8/)

    REAL   (SHR_KIND_R8), PARAMETER :: ecphas(pecclen) = & ! phases for eccen/fvelp cos/sin series
         &      (/    28.620089_SHR_KIND_R8, 193.788772_SHR_KIND_R8, 308.307024_SHR_KIND_R8,  &
         &           320.199637_SHR_KIND_R8, 279.376984_SHR_KIND_R8,  87.195000_SHR_KIND_R8,  &
         &           349.129677_SHR_KIND_R8, 128.443387_SHR_KIND_R8, 154.143880_SHR_KIND_R8,  &
         &           291.269597_SHR_KIND_R8, 114.860583_SHR_KIND_R8, 332.092251_SHR_KIND_R8,  &
         &           296.414411_SHR_KIND_R8, 145.769910_SHR_KIND_R8, 337.237063_SHR_KIND_R8,  &
         &           152.092288_SHR_KIND_R8, 126.839891_SHR_KIND_R8, 210.667199_SHR_KIND_R8,  &
         &            72.108838_SHR_KIND_R8/)

    ! Sine series data for computation of moving vernal equinox longitude of 
    ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).      

    REAL   (SHR_KIND_R8), PARAMETER :: mvamp (pmvelen) = & ! amplitudes for mvelp sine series 
         &      (/   7391.0225890_SHR_KIND_R8, 2555.1526947_SHR_KIND_R8, 2022.7629188_SHR_KIND_R8,  &
         &          -1973.6517951_SHR_KIND_R8, 1240.2321818_SHR_KIND_R8,  953.8679112_SHR_KIND_R8,  &
         &           -931.7537108_SHR_KIND_R8,  872.3795383_SHR_KIND_R8,  606.3544732_SHR_KIND_R8,  &
         &           -496.0274038_SHR_KIND_R8,  456.9608039_SHR_KIND_R8,  346.9462320_SHR_KIND_R8,  &
         &           -305.8412902_SHR_KIND_R8,  249.6173246_SHR_KIND_R8, -199.1027200_SHR_KIND_R8,  &
         &            191.0560889_SHR_KIND_R8, -175.2936572_SHR_KIND_R8,  165.9068833_SHR_KIND_R8,  &
         &            161.1285917_SHR_KIND_R8,  139.7878093_SHR_KIND_R8, -133.5228399_SHR_KIND_R8,  &
         &            117.0673811_SHR_KIND_R8,  104.6907281_SHR_KIND_R8,   95.3227476_SHR_KIND_R8,  &
         &             86.7824524_SHR_KIND_R8,   86.0857729_SHR_KIND_R8,   70.5893698_SHR_KIND_R8,  &
         &            -69.9719343_SHR_KIND_R8,  -62.5817473_SHR_KIND_R8,   61.5450059_SHR_KIND_R8,  &
         &            -57.9364011_SHR_KIND_R8,   57.1899832_SHR_KIND_R8,  -57.0236109_SHR_KIND_R8,  &
         &            -54.2119253_SHR_KIND_R8,   53.2834147_SHR_KIND_R8,   52.1223575_SHR_KIND_R8,  &
         &            -49.0059908_SHR_KIND_R8,  -48.3118757_SHR_KIND_R8,  -45.4191685_SHR_KIND_R8,  &
         &            -42.2357920_SHR_KIND_R8,  -34.7971099_SHR_KIND_R8,   34.4623613_SHR_KIND_R8,  &
         &            -33.8356643_SHR_KIND_R8,   33.6689362_SHR_KIND_R8,  -31.2521586_SHR_KIND_R8,  &
         &            -30.8798701_SHR_KIND_R8,   28.4640769_SHR_KIND_R8,  -27.1960802_SHR_KIND_R8,  &
         &             27.0860736_SHR_KIND_R8,  -26.3437456_SHR_KIND_R8,   24.7253740_SHR_KIND_R8,  &
         &             24.6732126_SHR_KIND_R8,   24.4272733_SHR_KIND_R8,   24.0127327_SHR_KIND_R8,  &
         &             21.7150294_SHR_KIND_R8,  -21.5375347_SHR_KIND_R8,   18.1148363_SHR_KIND_R8,  &
         &            -16.9603104_SHR_KIND_R8,  -16.1765215_SHR_KIND_R8,   15.5567653_SHR_KIND_R8,  &
         &             15.4846529_SHR_KIND_R8,   15.2150632_SHR_KIND_R8,   14.5047426_SHR_KIND_R8,  &
         &            -14.3873316_SHR_KIND_R8,   13.1351419_SHR_KIND_R8,   12.8776311_SHR_KIND_R8,  &
         &             11.9867234_SHR_KIND_R8,   11.9385578_SHR_KIND_R8,   11.7030822_SHR_KIND_R8,  &
         &             11.6018181_SHR_KIND_R8,  -11.2617293_SHR_KIND_R8,  -10.4664199_SHR_KIND_R8,  &
         &             10.4333970_SHR_KIND_R8,  -10.2377466_SHR_KIND_R8,   10.1934446_SHR_KIND_R8,  &
         &            -10.1280191_SHR_KIND_R8,   10.0289441_SHR_KIND_R8,  -10.0034259_SHR_KIND_R8/)

    REAL   (SHR_KIND_R8), PARAMETER :: mvrate(pmvelen) = & ! rates for mvelp sine series 
         &      (/    31.609974_SHR_KIND_R8, 32.620504_SHR_KIND_R8, 24.172203_SHR_KIND_R8,   &
         &             0.636717_SHR_KIND_R8, 31.983787_SHR_KIND_R8,  3.138886_SHR_KIND_R8,   &
         &            30.973257_SHR_KIND_R8, 44.828336_SHR_KIND_R8,  0.991874_SHR_KIND_R8,   &
         &             0.373813_SHR_KIND_R8, 43.668246_SHR_KIND_R8, 32.246691_SHR_KIND_R8,   &
         &            30.599444_SHR_KIND_R8,  2.147012_SHR_KIND_R8, 10.511172_SHR_KIND_R8,   &
         &            42.681324_SHR_KIND_R8, 13.650058_SHR_KIND_R8,  0.986922_SHR_KIND_R8,   &
         &             9.874455_SHR_KIND_R8, 13.013341_SHR_KIND_R8,  0.262904_SHR_KIND_R8,   &
         &             0.004952_SHR_KIND_R8,  1.142024_SHR_KIND_R8, 63.219948_SHR_KIND_R8,   &
         &             0.205021_SHR_KIND_R8,  2.151964_SHR_KIND_R8, 64.230478_SHR_KIND_R8,   &
         &            43.836462_SHR_KIND_R8, 47.439436_SHR_KIND_R8,  1.384343_SHR_KIND_R8,   &
         &             7.437771_SHR_KIND_R8, 18.829299_SHR_KIND_R8,  9.500642_SHR_KIND_R8,   &
         &             0.431696_SHR_KIND_R8,  1.160090_SHR_KIND_R8, 55.782177_SHR_KIND_R8,   &
         &            12.639528_SHR_KIND_R8,  1.155138_SHR_KIND_R8,  0.168216_SHR_KIND_R8,   &
         &             1.647247_SHR_KIND_R8, 10.884985_SHR_KIND_R8,  5.610937_SHR_KIND_R8,   &
         &            12.658184_SHR_KIND_R8,  1.010530_SHR_KIND_R8,  1.983748_SHR_KIND_R8,   &
         &            14.023871_SHR_KIND_R8,  0.560178_SHR_KIND_R8,  1.273434_SHR_KIND_R8,   &
         &            12.021467_SHR_KIND_R8, 62.583231_SHR_KIND_R8, 63.593761_SHR_KIND_R8,   &
         &            76.438310_SHR_KIND_R8,  4.280910_SHR_KIND_R8, 13.218362_SHR_KIND_R8,   &
         &            17.818769_SHR_KIND_R8,  8.359495_SHR_KIND_R8, 56.792707_SHR_KIND_R8,   &
         &            8.448301_SHR_KIND_R8,  1.978796_SHR_KIND_R8,  8.863925_SHR_KIND_R8,   &
         &             0.186365_SHR_KIND_R8,  8.996212_SHR_KIND_R8,  6.771027_SHR_KIND_R8,   &
         &            45.815258_SHR_KIND_R8, 12.002811_SHR_KIND_R8, 75.278220_SHR_KIND_R8,   &
         &            65.241008_SHR_KIND_R8, 18.870667_SHR_KIND_R8, 22.009553_SHR_KIND_R8,   &
         &            64.604291_SHR_KIND_R8, 11.498094_SHR_KIND_R8,  0.578834_SHR_KIND_R8,   &
         &             9.237738_SHR_KIND_R8, 49.747842_SHR_KIND_R8,  2.147012_SHR_KIND_R8,   &
         &             1.196895_SHR_KIND_R8,  2.133898_SHR_KIND_R8,  0.173168_SHR_KIND_R8/)

    REAL   (SHR_KIND_R8), PARAMETER :: mvphas(pmvelen) = & ! phases for mvelp sine series
         &      (/    251.9025_SHR_KIND_R8, 280.8325_SHR_KIND_R8, 128.3057_SHR_KIND_R8,   &
         &            348.1074_SHR_KIND_R8, 292.7252_SHR_KIND_R8, 165.1686_SHR_KIND_R8,   &
         &            263.7951_SHR_KIND_R8,  15.3747_SHR_KIND_R8,  58.5749_SHR_KIND_R8,   &
         &             40.8226_SHR_KIND_R8, 308.4258_SHR_KIND_R8, 240.0099_SHR_KIND_R8,   &
         &            222.9725_SHR_KIND_R8, 106.5937_SHR_KIND_R8, 114.5182_SHR_KIND_R8,   &
         &            268.7809_SHR_KIND_R8, 279.6869_SHR_KIND_R8,  39.6448_SHR_KIND_R8,   &
         &            126.4108_SHR_KIND_R8, 291.5795_SHR_KIND_R8, 307.2848_SHR_KIND_R8,   &
         &             18.9300_SHR_KIND_R8, 273.7596_SHR_KIND_R8, 143.8050_SHR_KIND_R8,   &
         &            191.8927_SHR_KIND_R8, 125.5237_SHR_KIND_R8, 172.7351_SHR_KIND_R8,   &
         &            316.7998_SHR_KIND_R8, 319.6024_SHR_KIND_R8,  69.7526_SHR_KIND_R8,   &
         &            123.5968_SHR_KIND_R8, 217.6432_SHR_KIND_R8,  85.5882_SHR_KIND_R8,   &
         &            156.2147_SHR_KIND_R8,  66.9489_SHR_KIND_R8,  20.2082_SHR_KIND_R8,   &
         &            250.7568_SHR_KIND_R8,  48.0188_SHR_KIND_R8,   8.3739_SHR_KIND_R8,   &
         &             17.0374_SHR_KIND_R8, 155.3409_SHR_KIND_R8,  94.1709_SHR_KIND_R8,   &
         &            221.1120_SHR_KIND_R8,  28.9300_SHR_KIND_R8, 117.1498_SHR_KIND_R8,   &
         &            320.5095_SHR_KIND_R8, 262.3602_SHR_KIND_R8, 336.2148_SHR_KIND_R8,   &
         &            233.0046_SHR_KIND_R8, 155.6977_SHR_KIND_R8, 184.6277_SHR_KIND_R8,   &
         &            267.2772_SHR_KIND_R8,  78.9281_SHR_KIND_R8, 123.4722_SHR_KIND_R8,   &
         &            188.7132_SHR_KIND_R8, 180.1364_SHR_KIND_R8,  49.1382_SHR_KIND_R8,   &
         &            152.5268_SHR_KIND_R8,  98.2198_SHR_KIND_R8,  97.4808_SHR_KIND_R8,   &
         &            221.5376_SHR_KIND_R8, 168.2438_SHR_KIND_R8, 161.1199_SHR_KIND_R8,   &
         &             55.0196_SHR_KIND_R8, 262.6495_SHR_KIND_R8, 200.3284_SHR_KIND_R8,   &
         &            201.6651_SHR_KIND_R8, 294.6547_SHR_KIND_R8,  99.8233_SHR_KIND_R8,   &
         &            213.5577_SHR_KIND_R8, 154.1631_SHR_KIND_R8, 232.7153_SHR_KIND_R8,   &
         &            138.3034_SHR_KIND_R8, 204.6609_SHR_KIND_R8, 106.5938_SHR_KIND_R8,   &
         &            250.4676_SHR_KIND_R8, 332.3345_SHR_KIND_R8,  27.3039_SHR_KIND_R8/)

    !---------------------------Local variables----------------------------------
    INTEGER(SHR_KIND_IN) :: i       ! Index for series summations
    REAL   (SHR_KIND_R8) :: obsum   ! Obliquity series summation
    REAL   (SHR_KIND_R8) :: cossum  ! Cos series summation for eccentricity/fvelp
    REAL   (SHR_KIND_R8) :: sinsum  ! Sin series summation for eccentricity/fvelp
    REAL   (SHR_KIND_R8) :: fvelp   ! Fixed vernal equinox long of perihelion
    REAL   (SHR_KIND_R8) :: mvsum   ! mvelp series summation
    REAL   (SHR_KIND_R8) :: beta    ! Intermediate argument for lambm0
    REAL   (SHR_KIND_R8) :: years   ! Years to time of interest ( pos <=> future)
    REAL   (SHR_KIND_R8) :: eccen2  ! eccentricity squared
    REAL   (SHR_KIND_R8) :: eccen3  ! eccentricity cubed

    !-------------------------- Formats -----------------------------------------
    CHARACTER(*),PARAMETER :: svnID  = "SVN " // &
         "$Id: shr_orb_mod.F90 6752 2007-10-04 21:02:15Z jwolfe $"
    CHARACTER(*),PARAMETER :: svnURL = "SVN <unknown URL>" 
    !  character(*),parameter :: svnURL = "SVN " // &
    !  "$URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/branch_tags/
    !cesm1_0_rel_tags/cesm1_0_rel01_share3_100616/shr/shr_orb_mod.F90 $"
    CHARACTER(len=*),PARAMETER :: F00 = "('(shr_orb_params) ',4a)"
    CHARACTER(len=*),PARAMETER :: F01 = "('(shr_orb_params) ',a,i9)"
    CHARACTER(len=*),PARAMETER :: F02 = "('(shr_orb_params) ',a,f6.3)"
    CHARACTER(len=*),PARAMETER :: F03 = "('(shr_orb_params) ',a,es14.6)"

    !----------------------------------------------------------------------------
    ! radinp and algorithms below will need a degree to radian conversion factor

    IF ( log_print .AND. s_loglev > 0 ) THEN
       WRITE(s_logunit,F00) 'Calculate characteristics of the orbit:'
       WRITE(s_logunit,F00) svnID
       WRITE(s_logunit,F00) svnURL
    END IF

    ! Check for flag to use input orbit parameters

    IF ( iyear_AD == SHR_ORB_UNDEF_INT ) THEN

       ! Check input obliq, eccen, and mvelp to ensure reasonable

       IF( obliq == SHR_ORB_UNDEF_REAL )THEN
          WRITE(s_logunit,F00) 'Have to specify orbital parameters:'
          WRITE(s_logunit,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
          WRITE(s_logunit,F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
          WRITE(s_logunit,F00) 'obliq, eccen, mvelp specify the orbit directly:'
          WRITE(s_logunit,F00) 'The AMIP II settings (for a 1995 orbit) are: '
          WRITE(s_logunit,F00) ' obliq =  23.4441'
          WRITE(s_logunit,F00) ' eccen =   0.016715'
          WRITE(s_logunit,F00) ' mvelp = 102.7'
          CALL shr_sys_abort()
       ELSE IF ( log_print ) THEN
          WRITE(s_logunit,F00) 'Use input orbital parameters: '
       END IF
       IF( (obliq < SHR_ORB_OBLIQ_MIN).OR.(obliq > SHR_ORB_OBLIQ_MAX) ) THEN
          WRITE(s_logunit,F03) 'Input obliquity unreasonable: ', obliq
          CALL shr_sys_abort()
       END IF
       IF( (eccen < SHR_ORB_ECCEN_MIN).OR.(eccen > SHR_ORB_ECCEN_MAX) ) THEN
          WRITE(s_logunit,F03) 'Input eccentricity unreasonable: ', eccen
          CALL shr_sys_abort()
       END IF
       IF( (mvelp < SHR_ORB_MVELP_MIN).OR.(mvelp > SHR_ORB_MVELP_MAX) ) THEN
          WRITE(s_logunit,F03) 'Input mvelp unreasonable: ' , mvelp
          CALL shr_sys_abort()
       END IF
       eccen2 = eccen*eccen
       eccen3 = eccen2*eccen

    ELSE  ! Otherwise calculate based on years before present

       IF ( log_print .AND. s_loglev > 0) THEN
          WRITE(s_logunit,F01) 'Calculate orbit for year: ' , iyear_AD
       END IF
       yb4_1950AD = 1950.0_SHR_KIND_R8 - REAL(iyear_AD,SHR_KIND_R8)
       IF ( ABS(yb4_1950AD) .GT. 1000000.0_SHR_KIND_R8 )THEN
          WRITE(s_logunit,F00) 'orbit only valid for years+-1000000'
          WRITE(s_logunit,F00) 'Relative to 1950 AD'
          WRITE(s_logunit,F03) '# of years before 1950: ',yb4_1950AD
          WRITE(s_logunit,F01) 'Year to simulate was  : ',iyear_AD
          CALL shr_sys_abort()
       END IF

       ! The following calculates the earths obliquity, orbital eccentricity
       ! (and various powers of it) and vernal equinox mean longitude of
       ! perihelion for years in the past (future = negative of years past),
       ! using constants (see parameter section) given in the program of:
       !
       ! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
       ! of Daily Insolation.  Contribution 18, Institute of Astronomy and
       ! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
       !
       ! and formulas given in the paper (where less precise constants are also
       ! given):
       !
       ! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
       ! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
       !
       ! The algorithm is valid only to 1,000,000 years past or hence.
       ! For a solution valid to 5-10 million years past see the above author.
       ! Algorithm below is better for years closer to present than is the
       ! 5-10 million year solution.
       !
       ! Years to time of interest must be negative of years before present
       ! (1950) in formulas that follow. 

       years = - yb4_1950AD

       ! In the summations below, cosine or sine arguments, which end up in
       ! degrees, must be converted to radians via multiplication by degrad.
       !
       ! Summation of cosine series for obliquity (epsilon in Berger 1978) in
       ! degrees. Convert the amplitudes and rates, which are in arc secs, into
       ! degrees via multiplication by psecdeg (arc seconds to degrees conversion
       ! factor).  For obliq, first term is Berger 1978 epsilon star; second
       ! term is series summation in degrees.

       obsum = 0.0_SHR_KIND_R8
       DO i = 1, poblen
          obsum = obsum + obamp(i)*psecdeg*COS((obrate(i)*psecdeg*years + &
               &       obphas(i))*degrad)
       END DO
       obliq = 23.320556_SHR_KIND_R8 + obsum

       ! Summation of cosine and sine series for computation of eccentricity 
       ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of 
       ! perihelion (fvelp; pi in Berger 1978), which is used for computation 
       ! of moving vernal equinox longitude of perihelion.  Convert the rates, 
       ! which are in arc seconds, into degrees via multiplication by psecdeg.

       cossum = 0.0_SHR_KIND_R8
       DO i = 1, pecclen
          cossum = cossum+ecamp(i)*COS((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
       END DO

       sinsum = 0.0_SHR_KIND_R8
       DO i = 1, pecclen
          sinsum = sinsum+ecamp(i)*SIN((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
       END DO

       ! Use summations to calculate eccentricity

       eccen2 = cossum*cossum + sinsum*sinsum
       eccen  = SQRT(eccen2)
       eccen3 = eccen2*eccen

       ! A series of cases for fvelp, which is in radians.

       IF (ABS(cossum) .LE. 1.0E-8_SHR_KIND_R8) THEN
          IF (sinsum .EQ. 0.0_SHR_KIND_R8) THEN
             fvelp = 0.0_SHR_KIND_R8
          ELSE IF (sinsum .LT. 0.0_SHR_KIND_R8) THEN
             fvelp = 1.5_SHR_KIND_R8*pi
          ELSE IF (sinsum .GT. 0.0_SHR_KIND_R8) THEN
             fvelp = .5_SHR_KIND_R8*pi
          ENDIF
       ELSE IF (cossum .LT. 0.0_SHR_KIND_R8) THEN
          fvelp = ATAN(sinsum/cossum) + pi
       ELSE IF (cossum .GT. 0.0_SHR_KIND_R8) THEN
          IF (sinsum .LT. 0.0_SHR_KIND_R8) THEN
             fvelp = ATAN(sinsum/cossum) + 2.0_SHR_KIND_R8*pi
          ELSE
             fvelp = ATAN(sinsum/cossum)
          ENDIF
       ENDIF

       ! Summation of sin series for computation of moving vernal equinox long
       ! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
       ! first term is fvelp in degrees; second term is Berger 1978 psi bar 
       ! times years and in degrees; third term is Berger 1978 zeta; fourth 
       ! term is series summation in degrees.  Convert the amplitudes and rates,
       ! which are in arc seconds, into degrees via multiplication by psecdeg.  
       ! Series summation plus second and third terms constitute Berger 1978
       ! psi, which is the general precession.

       mvsum = 0.0_SHR_KIND_R8
       DO i = 1, pmvelen
          mvsum = mvsum + mvamp(i)*psecdeg*SIN((mvrate(i)*psecdeg*years + &
               &       mvphas(i))*degrad)
       END DO
       mvelp = fvelp/degrad + 50.439273_SHR_KIND_R8*psecdeg*years + 3.392506_SHR_KIND_R8 + mvsum

       ! Cases to make sure mvelp is between 0 and 360.

       DO WHILE (mvelp .LT. 0.0_SHR_KIND_R8)
          mvelp = mvelp + 360.0_SHR_KIND_R8
       END DO
       DO WHILE (mvelp .GE. 360.0_SHR_KIND_R8)
          mvelp = mvelp - 360.0_SHR_KIND_R8
       END DO

    END IF  ! end of test on whether to calculate or use input orbital params

    ! Orbit needs the obliquity in radians

    obliqr = obliq*degrad

    ! 180 degrees must be added to mvelp since observations are made from the
    ! earth and the sun is considered (wrongly for the algorithm) to go around
    ! the earth. For a more graphic explanation see Appendix B in:
    !
    ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
    ! Periods.  J. of Geophysical Research 98:10,341-10,362.
    !
    ! Additionally, orbit will need this value in radians. So mvelp becomes
    ! mvelpp (mvelp plus pi)

    mvelpp = (mvelp + 180._SHR_KIND_R8)*degrad

    ! Set up an argument used several times in lambm0 calculation ahead.

    beta = SQRT(1._SHR_KIND_R8 - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in 
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 = 2._SHR_KIND_R8*((.5_SHR_KIND_R8*eccen + .125_SHR_KIND_R8*eccen3)*(1._SHR_KIND_R8 + beta)*SIN(mvelpp)  &
         &      - .250_SHR_KIND_R8*eccen2*(.5_SHR_KIND_R8    + beta)*SIN(2._SHR_KIND_R8*mvelpp)            &
         &      + .125_SHR_KIND_R8*eccen3*(1._SHR_KIND_R8/3._SHR_KIND_R8 + beta)*SIN(3._SHR_KIND_R8*mvelpp))

    IF ( log_print ) THEN
       WRITE(s_logunit,F03) '------ Computed Orbital Parameters ------'
       WRITE(s_logunit,F03) 'Eccentricity      = ',eccen
       WRITE(s_logunit,F03) 'Obliquity (deg)   = ',obliq
       WRITE(s_logunit,F03) 'Obliquity (rad)   = ',obliqr
       WRITE(s_logunit,F03) 'Long of perh(deg) = ',mvelp
       WRITE(s_logunit,F03) 'Long of perh(rad) = ',mvelpp
       WRITE(s_logunit,F03) 'Long at v.e.(rad) = ',lambm0
       WRITE(s_logunit,F03) '-----------------------------------------'
    END IF
     CALL shr_orb_print( iyear_AD )

  END SUBROUTINE shr_orb_params

  !===============================================================================

  SUBROUTINE shr_orb_decl(calday,delta,ratio)

    !-------------------------------------------------------------------------------
    !
    ! Compute earth/orbit parameters using formula suggested by
    ! Duane Thresher.
    !
    !---------------------------Code history----------------------------------------
    !
    ! Original version:  Erik Kluzek
    ! Date:              Oct/1997
    !
    !-------------------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    REAL   (SHR_KIND_R8),INTENT(in )  :: calday ! Calendar day, including fraction
    REAL   (SHR_KIND_R8),INTENT(out)  :: delta
    REAL   (SHR_KIND_R8),INTENT(out)  :: ratio
    !---------------------------Local variables-----------------------------
    REAL   (SHR_KIND_R8),PARAMETER :: dayspy = 365.0_SHR_KIND_R8  ! days per year
    REAL   (SHR_KIND_R8),PARAMETER :: ve     = 80.5_SHR_KIND_R8   ! Calday of vernal equinox
    ! assumes Jan 1 = calday 1

    REAL   (SHR_KIND_R8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
    REAL   (SHR_KIND_R8) ::   lmm    ! Intermediate argument involving lambm
    REAL   (SHR_KIND_R8) ::   lamb   ! Lambda, the earths long of perihelion
    REAL   (SHR_KIND_R8) ::   invrho ! Inverse normalized sun/earth distance
    REAL   (SHR_KIND_R8) ::   sinl   ! Sine of lmm

    ! Compute eccentricity factor and solar declination using
    ! day value where a round day (such as 213.0) refers to 0z at
    ! Greenwich longitude.
    !
    ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    ! 35:2362-2367.
    !
    ! To get the earths true longitude (position in orbit; lambda in Berger 
    ! 1978) which is necessary to find the eccentricity factor and declination,
    ! must first calculate the mean longitude (lambda m in Berger 1978) at
    ! the present day.  This is done by adding to lambm0 (the mean longitude
    ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
    ! an increment (declin lambda m in Berger 1978) that is the number of
    ! days past or before (a negative increment) the vernal equinox divided by
    ! the days in a model year times the 2*pi radians in a complete orbit.

    lambm = lambm0 + (calday - ve)*2._SHR_KIND_R8*pi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = SIN(lmm)
    lamb  = lambm  + eccen*(2._SHR_KIND_R8*sinl + eccen*(1.25_SHR_KIND_R8*SIN(2._SHR_KIND_R8*lmm)  &
         &     + eccen*((13.0_SHR_KIND_R8/12.0_SHR_KIND_R8)*SIN(3._SHR_KIND_R8*lmm) - 0.25_SHR_KIND_R8*sinl)))

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (declin)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
    ! rho will be used), and thus the eccentricity factor (eccf), can be 
    ! calculated from formulas given in Berger 1978.

    invrho = (1._SHR_KIND_R8 + eccen*COS(lamb - mvelpp)) / (1._SHR_KIND_R8 - eccen*eccen)

    ! Set solar declination and eccentricity factor

    declin  = ASIN(SIN(obliqr)*SIN(lamb))
    eccf   = invrho*invrho

    delta=declin
    !ratio   =  1.0/invrho
    ratio   = invrho
    RETURN

  END SUBROUTINE shr_orb_decl

  !===============================================================================

  SUBROUTINE shr_orb_print( &
    iyear_AD   )
    !-------------------------------------------------------------------------------
    !
    ! Print out the information on the Earths input orbital characteristics
    !
    !---------------------------Code history----------------------------------------
    !
    ! Original version:  Erik Kluzek
    ! Date:              Oct/1997
    !
    !-------------------------------------------------------------------------------

    !---------------------------Arguments----------------------------------------
    INTEGER(SHR_KIND_IN),INTENT(in) :: iyear_AD ! requested Year (AD)
    !REAL   (SHR_KIND_R8),INTENT(in) :: eccen    ! eccentricity (unitless) 
    ! (typically 0 to 0.1)
    !REAL   (SHR_KIND_R8),INTENT(in) :: obliq    ! obliquity (-90 to +90 degrees) 
    ! typically 22-26
    !REAL   (SHR_KIND_R8),INTENT(in) :: mvelp    ! moving vernal equinox at perhel
    ! (0 to 360 degrees)
    !-------------------------- Formats -----------------------------------------
    CHARACTER(len=*),PARAMETER :: F00 = "('(shr_orb_print) ',4a)"
    CHARACTER(len=*),PARAMETER :: F01 = "('(shr_orb_print) ',a,i9.4)"
    CHARACTER(len=*),PARAMETER :: F02 = "('(shr_orb_print) ',a,f6.3)"
    CHARACTER(len=*),PARAMETER :: F03 = "('(shr_orb_print) ',a,es14.6)"
    !----------------------------------------------------------------------------

    IF (s_loglev > 0) THEN
       IF ( iyear_AD .NE. SHR_ORB_UNDEF_INT ) THEN
          IF ( iyear_AD > 0 ) THEN
             WRITE(s_logunit,F01) 'Orbital parameters calculated for year: AD ',iyear_AD
          ELSE
             WRITE(s_logunit,F01) 'Orbital parameters calculated for year: BC ',iyear_AD
          END IF
       ELSE IF ( obliq /= SHR_ORB_UNDEF_REAL ) THEN
          WRITE(s_logunit,F03) 'Orbital parameters: '
          WRITE(s_logunit,F03) 'Obliquity (degree):              ', obliq
          WRITE(s_logunit,F03) 'Eccentricity (unitless):         ', eccen
          WRITE(s_logunit,F03) 'Long. of moving Perhelion (deg): ', mvelp
       ELSE
          WRITE(s_logunit,F03) 'Orbit parameters not set!'
       END IF
    ENDIF

  END SUBROUTINE shr_orb_print
  !===============================================================================
  SUBROUTINE shr_sys_abort()

    IMPLICIT NONE

!    CHARACTER(*)        ,OPTIONAL :: string  ! error message string
!    INTEGER(SHR_KIND_IN),OPTIONAL :: rc      ! error code
    STOP "shr_sys_abort"
  END SUBROUTINE shr_sys_abort
END MODULE shr_orb_mod


MODULE ModRadiationDriver

  ! InitRadiationDriver ---| InitRadiation
  !                        !
  !                        ! InitRadtim
  !                        !
  !                        ! InitGetoz
  !
  ! RadiationDriver ---| rqvirt
  !                    !
  !                    | spmrad ------| radtim
  !                    !              !
  !                    !              ! getoz
  !                    !              !
  !                    !              ! cldgen, cldgn2 or cldgn3
  !                    !              !
  !                    !              ! swrad, cliradsw or ukmet
  !                    !              !
  !                    !              ! lwrad, cliradlw or ukmet
  !                    !
  !                    ! RadDiagStor
  !                    !
  !                    ! RadGridHistStorage

!tar begin

  USE Parallelism, ONLY : myid 

!tar end

  USE Constants, ONLY :     &
       grav,cp,stefan, gm2dob, &
       hl, rmwmdi, gasr, e0c, rmwmd, delq, qmin, &
       r8, i8, pai, solcon

  USE Utils, ONLY: &
       tmstmp2
       
  USE Options, ONLY : &
       ilcon          ,&
       iccon          ,&
       ilwrad         ,& 
       iswrad         ,& 
       crdcld, tbase, nfprt, nferr, ifozone, nfctrl,iyear_AD,schemes,fNameCldOptSW ,&
       fNameCldOptLW ,jull,& 
!tar begin
       ifaeros
!tar end        

  USE Rad_COLA, ONLY:        &
       InitRadCOLA         , &
       swrad, lwrad

  USE Rad_Clirad, ONLY: cliradsw,InitCliradSW
  USE Rad_Cliradlw, ONLY: cliradlw,InitCliRadLW

  USE Rad_CliradTarasova, ONLY: CliradTarasova_sw,InitCliradTarasova_sw
  USE   Rad_CliradlwTarasova, ONLY:  CliradTarasova_lw,InitCliradTarasova_lw
  
  USE Rad_UKMO, ONLY: ukmo_swintf, ukmo_lwintf, InitRadUKMO

  USE Rad_RRTMG, ONLY:  Init_Rad_RRTMG,Run_Rad_RRTMG_SW,Run_Rad_RRTMG_LW,Finalize_Rad_RRTMG

  USE CloudOpticalProperty, ONLY: Init_Optical_Properties,RunCloudOpticalProperty2,Cloud_Micro_WRF,gethml

  USE Diagnostics, ONLY: &
       dodia       , &
       updia       , &
       StartStorDiag,&
       nDiag_cloudc, & ! cloud cover
       nDiag_lwdbot, & ! longwave downward at bottom
       nDiag_lwubot, & ! longwave upward at bottom
       nDiag_lwutop, & ! longwave upward at top
       nDiag_swdtop, & ! shortwave downward at top
       nDiag_swdbot, & ! shortwave downward at ground
       nDiag_swubot, & ! shortwave upward at bottom
       nDiag_swutop, & ! shortwave upward at top
       nDiag_swabea, & ! shortwave absorbed by the earth/atmosphere
       nDiag_swabgr, & ! shortwave absorbed by the ground
       nDiag_lwnetb, & ! net longwave at bottom
       nDiag_lwheat, & ! longwave heating
       nDiag_swheat, & ! shortwave heating
       nDiag_lwdbtc, & ! longwave downward at bottom (clear)
       nDiag_lwutpc, & ! longwave upward at top (clear)
       nDiag_swdbtc, & ! shortwave downward at ground (clear)
       nDiag_swubtc, & ! shortwave upward at bottom (clear)
       nDiag_swutpc, & ! shortwave upward at top (clear)
       nDiag_swaeac, & ! shortwave absorbed by the earth/atmosphere (clear)
       nDiag_swabgc, & ! shortwave absorbed by the ground (clear)
       nDiag_lwnbtc, & ! net longwave at bottom (clear)
       nDiag_vdtclc, & ! vertical dist total cloud cover
       nDiag_invcld, & ! inversion cloud
       nDiag_ssatcl, & ! supersaturation cloud
       nDiag_cnvcld, & ! convective cloud
       nDiag_shcvcl, & ! shallow convective cloud
       nDiag_clliwp, & ! cloud liquid water path
       nDiag_lwcemi, & ! longwave cloud emissivity
       nDiag_sclopd, & ! shortwave cloud optical depth
       nDiag_lwhtcl, & ! longwave heating (clear)
       nDiag_swhtcl, & ! shortwave heating (clear)
       nDiag_ozonmr, & ! ozone mass mixing ratio (g/g)
       nDiag_viozoc, & ! Vertically Integrated Ozone Content (Dobson units)
       nDiag_iceper, & ! Ice particle Effective Radius (microns)
       nDiag_liqper, & ! Liquid particle Effective Radius (microns)
       nDiag_co2aer    ! mix CO2 Concentration kg/kg

  USE GridHistory, ONLY:       &
       IsGridHistoryOn, StoreGridHistory, StoreMaskedGridHistory, dogrh, &
       nGHis_casrrs, nGHis_mofres, nGHis_vdtclc, nGHis_lwheat, nGHis_swheat, &
       nGHis_lwutop, nGHis_lwdbot, nGHis_coszen, nGHis_swdtop, nGHis_swdbvb, &
       nGHis_swdbvd, nGHis_swdbnb, nGHis_swdbnd, nGHis_vibalb, nGHis_vidalb, &
       nGHis_nibalb, nGHis_nidalb, nGHis_hcseai, nGHis_hsseai, nGHis_cloudc, &
       nGHis_dragcf, nGHis_nrdcan, nGHis_nrdgsc, nGHis_cascrs, nGHis_casgrs, &
       nGHis_canres, nGHis_gcovrs, nGHis_bssfrs, nGHis_ecairs, nGHis_tcairs, &
       nGHis_shfcan, nGHis_shfgnd, nGHis_tracan, nGHis_tragcv, nGHis_inlocp, &
       nGHis_inlogc, nGHis_bsevap, nGHis_canhea, nGHis_gcheat, nGHis_runoff, &
       nGHis_vdheat, nGHis_vduzon, nGHis_vdvmer, nGHis_vdmois, nGHis_swutop, &
       nGHis_lwubot, nGHis_ustres, nGHis_vstres, nGHis_sheatf, nGHis_lheatf, &
       nGHis_swdgrd, nGHis_swugrd, nGHis_cldlow, nGHis_cldmed, nGHis_cldHig
 USE shr_orb_mod ,Only : shr_orb_params,shr_orb_decl

 USE wv_saturation,Only : estblf,findsp

 USE PhysicalFunctions,Only : fpvs2es5

  USE Parallelism  , ONLY:       &
       MsgOne, myid  

 USE Watches, ONLY:  &
       ChangeWatch

   IMPLICIT NONE
  SAVE

 PRIVATE

 PUBLIC :: InitRadiationDriver
 PUBLIC :: RadiationDriver
 PUBLIC :: DestroyRadiationDriver

 ! General radiation routines
 PUBLIC :: COSZMED
 PUBLIC :: radtim
 PUBLIC :: getoz

 ! Ozone data
 INTEGER, PARAMETER :: nlm_getoz=18
 REAL(KIND=r8), ALLOCATABLE :: ozone(:,:,:)
 LOGICAL          :: first_getoz
 INTEGER          :: mon_getoz
 REAL(KIND=r8)    :: year_getoz
 LOGICAL          :: inter_getoz
 REAL(KIND=r8)    :: ozsig(18)    !ozsig(nlm_getoz)

 INTEGER :: monday(12)

 ! Usefull constants
 REAL(KIND=r8), PARAMETER    :: pai12 = pai/12.e0_r8
 REAL(KIND=r8), PARAMETER    :: pai2i = 1.0e0_r8/(2.0e0_r8*pai)
 REAL(KIND=r8), PARAMETER    :: fim24 = 24.0e0_r8/360.0e0_r8
 REAL(KIND=r8), PARAMETER    :: tmelt =273.16_r8
 REAL(KIND=r8), PARAMETER    :: pptop = 0.005_r8  ! Model-top presure
 INTEGER      , PARAMETER    :: nbndsw=14
 INTEGER      , PARAMETER    :: nbndlw=16

CONTAINS

 SUBROUTINE InitRadiationDriver(monl,yrl,kmax,a_hybr,b_hybr,dt,nls)
   IMPLICIT NONE
   INTEGER      , INTENT(in   ) :: monl(12)
   REAL(KIND=r8), INTENT(IN   ) :: yrl
   INTEGER      , INTENT(IN   ) :: kmax,nls
   REAL(KIND=r8), INTENT(in   ) :: a_hybr(kmax+1)
   REAL(KIND=r8), INTENT(in   ) :: b_hybr(kmax+1)
   REAL(KIND=r8), INTENT(in   ) :: dt
  
   IF (TRIM(iswrad).eq.'LCH'.OR.TRIM(ilwrad).eq.'HRS') THEN
      CALL InitRadCOLA()
   ENDIF
   IF (TRIM(iswrad).eq.'CRD'.OR.TRIM(ilwrad).eq.'CRD') THEN
      !CALL InitRadClirad
      CALL InitCliradSW()
      CALL InitCliRadLW()
   ENDIF
   IF (TRIM(iswrad).eq.'CRDTF'.OR.TRIM(ilwrad).eq.'CRDTF') THEN
      CALL InitCliradTarasova_sw()
      CALL InitCliradTarasova_lw()
   ENDIF
   IF (TRIM(iswrad).eq.'UKM'.OR.TRIM(ilwrad).eq.'UKM') THEN
      CALL InitRadUKMO(kmax,a_hybr,b_hybr,nls)
   ENDIF
   IF (TRIM(iswrad).eq.'RRTMG'.OR.TRIM(ilwrad).eq.'RRTMG') THEN
      CALL Init_Rad_RRTMG(a_hybr,b_hybr,kMax,nbndsw,nbndlw)
   ENDIF

   CALL shr_orb_params(iyear_AD)

   CALL InitRadtim(monl)

!  CALL InitGetoz(yrl,kmax,sl)
   CALL InitGetoz(yrl,kmax)
   
   CALL Init_Optical_Properties(dt,kmax,a_hybr,b_hybr,nbndsw,nbndlw,fNameCldOptSW ,fNameCldOptLW )

 END SUBROUTINE InitRadiationDriver

  SUBROUTINE InitRadtim(monl)
    INTEGER, INTENT(in ) :: monl(12)
    INTEGER              :: m
    monday(1)=0
    DO m=2,12
       monday(m)=monday(m-1)+monl(m-1)
    END DO
  END SUBROUTINE InitRadtim

! SUBROUTINE InitGetoz(yrl,kmax,sl)
  SUBROUTINE InitGetoz(yrl,kmax)
    IMPLICIT NONE
    REAL(KIND=r8),    INTENT(IN   ) :: yrl
    INTEGER, INTENT(IN   ) :: kmax
!   REAL(KIND=r8),    INTENT(in   ) :: sl (kmax)
    INTEGER                :: l, ll

    INTEGER, PARAMETER :: nl=37
    INTEGER, PARAMETER :: ns=4

    ALLOCATE(ozone(nlm_getoz,nl,ns))

    !
    !     four season climatological ozone data in nmc sigma layers
    !
    !     for seasonal variation
    !     season=1 - winter          season=2 - spring
    !     season=3 - summer          season=4 - fall
    !     unit of ozone mixing ratio is in ( 10**-4 g/g ).  the data is
    !     in 18 sigma layers from top to bottom.  for every layer, there
    !     are 37 latitudes at 5 degree interval from north pole to south
    !     pole.
    !     mrf86 18 layers
    !
    !
    !     1. winter
    !
    !     wint1(18,6)
    !
    ozone(1:18, 1:6, 1) = RESHAPE( (/ &
         .068467e0_r8,.052815e0_r8,.035175e0_r8,.022334e0_r8,.013676e0_r8,.007363e0_r8, &
         .003633e0_r8,.001582e0_r8,.001111e0_r8,.000713e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .069523e0_r8,.052249e0_r8,.034255e0_r8,.021379e0_r8,.012306e0_r8,.006727e0_r8, &
         .003415e0_r8,.001578e0_r8,.001072e0_r8,.000681e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070579e0_r8,.051684e0_r8,.033335e0_r8,.020423e0_r8,.010935e0_r8,.006091e0_r8, &
         .003197e0_r8,.001573e0_r8,.001034e0_r8,.000650e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074885e0_r8,.049987e0_r8,.030140e0_r8,.017894e0_r8,.009881e0_r8,.005543e0_r8, &
         .002907e0_r8,.001379e0_r8,.000961e0_r8,.000644e0_r8,.000512e0_r8,.000463e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .079190e0_r8,.048290e0_r8,.026945e0_r8,.015366e0_r8,.008826e0_r8,.004995e0_r8, &
         .002616e0_r8,.001184e0_r8,.000887e0_r8,.000637e0_r8,.000508e0_r8,.000486e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082443e0_r8,.047591e0_r8,.025358e0_r8,.014294e0_r8,.008233e0_r8,.004664e0_r8, &
         .002430e0_r8,.001068e0_r8,.000851e0_r8,.000644e0_r8,.000508e0_r8,.000474e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8/), &
         (/18,6/))
    !
    !     wint2(18,6)
    !
    ozone(1:18, 7:12, 1) = RESHAPE( (/ &
         .085695e0_r8,.046892e0_r8,.023772e0_r8,.013223e0_r8,.007640e0_r8,.004333e0_r8, &
         .002244e0_r8,.000951e0_r8,.000815e0_r8,.000650e0_r8,.000508e0_r8,.000463e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089618e0_r8,.042869e0_r8,.019963e0_r8,.010502e0_r8,.005966e0_r8,.003525e0_r8, &
         .001936e0_r8,.000906e0_r8,.000769e0_r8,.000625e0_r8,.000508e0_r8,.000452e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .093540e0_r8,.038846e0_r8,.016155e0_r8,.007781e0_r8,.004292e0_r8,.002716e0_r8, &
         .001628e0_r8,.000862e0_r8,.000724e0_r8,.000600e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .097097e0_r8,.034916e0_r8,.012983e0_r8,.006240e0_r8,.003666e0_r8,.002259e0_r8, &
         .001336e0_r8,.000730e0_r8,.000629e0_r8,.000549e0_r8,.000499e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .100654e0_r8,.030986e0_r8,.009812e0_r8,.004698e0_r8,.003041e0_r8,.001803e0_r8, &
         .001044e0_r8,.000599e0_r8,.000533e0_r8,.000499e0_r8,.000491e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .101724e0_r8,.026500e0_r8,.007228e0_r8,.003391e0_r8,.002058e0_r8,.001285e0_r8, &
         .000811e0_r8,.000531e0_r8,.000478e0_r8,.000449e0_r8,.000440e0_r8,.000421e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint3(18,6)
    !
    ozone(1:18, 13:18, 1) = RESHAPE( (/ &
         .102794e0_r8,.022015e0_r8,.004645e0_r8,.002084e0_r8,.001076e0_r8,.000767e0_r8, &
         .000577e0_r8,.000463e0_r8,.000423e0_r8,.000399e0_r8,.000389e0_r8,.000401e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .103456e0_r8,.018235e0_r8,.003195e0_r8,.001379e0_r8,.000771e0_r8,.000585e0_r8, &
         .000474e0_r8,.000411e0_r8,.000380e0_r8,.000362e0_r8,.000343e0_r8,.000348e0_r8, &
         .000346e0_r8,.000328e0_r8,.000317e0_r8,.000305e0_r8,.000302e0_r8,.000302e0_r8, &
         .104118e0_r8,.014455e0_r8,.001745e0_r8,.000674e0_r8,.000467e0_r8,.000403e0_r8, &
         .000370e0_r8,.000359e0_r8,.000337e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104106e0_r8,.012997e0_r8,.001479e0_r8,.000639e0_r8,.000468e0_r8,.000422e0_r8, &
         .000392e0_r8,.000372e0_r8,.000342e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104093e0_r8,.011539e0_r8,.001213e0_r8,.000604e0_r8,.000468e0_r8,.000442e0_r8, &
         .000414e0_r8,.000385e0_r8,.000347e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104087e0_r8,.010726e0_r8,.000971e0_r8,.000538e0_r8,.000440e0_r8,.000434e0_r8, &
         .000418e0_r8,.000397e0_r8,.000375e0_r8,.000343e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint4(18,6)
    !
    ozone(1:18, 19:24, 1) = RESHAPE( (/ &
         .102665e0_r8,.010977e0_r8,.001237e0_r8,.000590e0_r8,.000498e0_r8,.000479e0_r8, &
         .000458e0_r8,.000436e0_r8,.000421e0_r8,.000387e0_r8,.000326e0_r8,.000298e0_r8, &
         .000246e0_r8,.000227e0_r8,.000211e0_r8,.000200e0_r8,.000194e0_r8,.000186e0_r8, &
         .100892e0_r8,.012873e0_r8,.001886e0_r8,.000785e0_r8,.000643e0_r8,.000568e0_r8, &
         .000519e0_r8,.000487e0_r8,.000471e0_r8,.000437e0_r8,.000368e0_r8,.000305e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100534e0_r8,.013704e0_r8,.002028e0_r8,.000861e0_r8,.000701e0_r8,.000604e0_r8, &
         .000546e0_r8,.000513e0_r8,.000504e0_r8,.000462e0_r8,.000381e0_r8,.000307e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100218e0_r8,.015035e0_r8,.002537e0_r8,.001037e0_r8,.000790e0_r8,.000726e0_r8, &
         .000673e0_r8,.000628e0_r8,.000579e0_r8,.000512e0_r8,.000440e0_r8,.000374e0_r8, &
         .000307e0_r8,.000253e0_r8,.000227e0_r8,.000208e0_r8,.000194e0_r8,.000186e0_r8, &
         .099903e0_r8,.016365e0_r8,.003045e0_r8,.001214e0_r8,.000879e0_r8,.000848e0_r8, &
         .000801e0_r8,.000744e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000441e0_r8, &
         .000410e0_r8,.000358e0_r8,.000342e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .099547e0_r8,.017725e0_r8,.003693e0_r8,.001578e0_r8,.001125e0_r8,.000985e0_r8, &
         .000879e0_r8,.000795e0_r8,.000712e0_r8,.000643e0_r8,.000584e0_r8,.000521e0_r8, &
         .000482e0_r8,.000384e0_r8,.000351e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint5(18,6)
    !
    ozone(1:18, 25:30, 1) = RESHAPE( (/ &
         .099191e0_r8,.019085e0_r8,.004340e0_r8,.001943e0_r8,.001371e0_r8,.001122e0_r8, &
         .000957e0_r8,.000847e0_r8,.000770e0_r8,.000724e0_r8,.000669e0_r8,.000601e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .098107e0_r8,.020617e0_r8,.004758e0_r8,.002137e0_r8,.001516e0_r8,.001211e0_r8, &
         .000999e0_r8,.000848e0_r8,.000778e0_r8,.000730e0_r8,.000677e0_r8,.000603e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .097023e0_r8,.022148e0_r8,.005177e0_r8,.002332e0_r8,.001660e0_r8,.001300e0_r8, &
         .001041e0_r8,.000849e0_r8,.000786e0_r8,.000737e0_r8,.000686e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .093464e0_r8,.026177e0_r8,.008525e0_r8,.003892e0_r8,.002452e0_r8,.001609e0_r8, &
         .001116e0_r8,.000851e0_r8,.000809e0_r8,.000762e0_r8,.000690e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .089906e0_r8,.030206e0_r8,.011873e0_r8,.005453e0_r8,.003244e0_r8,.001918e0_r8, &
         .001192e0_r8,.000852e0_r8,.000832e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .080939e0_r8,.032414e0_r8,.014163e0_r8,.007241e0_r8,.004328e0_r8,.002522e0_r8, &
         .001481e0_r8,.000934e0_r8,.000861e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint6(18,6)
    !
    ozone(1:18, 31:36, 1) = RESHAPE( (/ &
         .071972e0_r8,.034622e0_r8,.016453e0_r8,.009029e0_r8,.005413e0_r8,.003127e0_r8, &
         .001770e0_r8,.001015e0_r8,.000890e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .069820e0_r8,.035028e0_r8,.016929e0_r8,.009389e0_r8,.005645e0_r8,.003260e0_r8, &
         .001843e0_r8,.001055e0_r8,.000905e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .067669e0_r8,.035434e0_r8,.017406e0_r8,.009749e0_r8,.005876e0_r8,.003393e0_r8, &
         .001916e0_r8,.001094e0_r8,.000920e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .065518e0_r8,.035975e0_r8,.017854e0_r8,.010100e0_r8,.006534e0_r8,.003985e0_r8, &
         .002321e0_r8,.001240e0_r8,.000966e0_r8,.000774e0_r8,.000640e0_r8,.000548e0_r8, &
         .000479e0_r8,.000384e0_r8,.000346e0_r8,.000316e0_r8,.000302e0_r8,.000302e0_r8, &
         .063367e0_r8,.036516e0_r8,.018302e0_r8,.010452e0_r8,.007192e0_r8,.004577e0_r8, &
         .002727e0_r8,.001387e0_r8,.001012e0_r8,.000762e0_r8,.000585e0_r8,.000490e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .061216e0_r8,.037359e0_r8,.019151e0_r8,.010633e0_r8,.006845e0_r8,.004382e0_r8, &
         .002691e0_r8,.001511e0_r8,.001061e0_r8,.000749e0_r8,.000568e0_r8,.000465e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint7(18)
    !
    ozone(1:18, 37, 1) = (/ &
         .059066e0_r8,.038201e0_r8,.019999e0_r8,.010813e0_r8,.006498e0_r8,.004188e0_r8, &
         .002656e0_r8,.001636e0_r8,.001110e0_r8,.000737e0_r8,.000551e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     2. spring
    !
    ozone(1:18, 1:6, 2) = RESHAPE( (/ &
         .074229e0_r8,.050084e0_r8,.030930e0_r8,.018676e0_r8,.011965e0_r8,.008165e0_r8, &
         .005428e0_r8,.003399e0_r8,.002098e0_r8,.001138e0_r8,.000780e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .074927e0_r8,.049459e0_r8,.029215e0_r8,.018025e0_r8,.011754e0_r8,.007786e0_r8, &
         .004972e0_r8,.002926e0_r8,.001817e0_r8,.001025e0_r8,.000758e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .075625e0_r8,.048835e0_r8,.027500e0_r8,.017375e0_r8,.011544e0_r8,.007407e0_r8, &
         .004516e0_r8,.002453e0_r8,.001536e0_r8,.000912e0_r8,.000737e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .077409e0_r8,.048159e0_r8,.026661e0_r8,.016596e0_r8,.010962e0_r8,.006972e0_r8, &
         .004160e0_r8,.002132e0_r8,.001391e0_r8,.000868e0_r8,.000686e0_r8,.000601e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .079194e0_r8,.047483e0_r8,.025822e0_r8,.015818e0_r8,.010380e0_r8,.006537e0_r8, &
         .003804e0_r8,.001811e0_r8,.001245e0_r8,.000825e0_r8,.000635e0_r8,.000570e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .084591e0_r8,.046553e0_r8,.025037e0_r8,.015156e0_r8,.009841e0_r8,.006124e0_r8, &
         .003534e0_r8,.001693e0_r8,.001170e0_r8,.000793e0_r8,.000631e0_r8,.000537e0_r8, &
         .000551e0_r8,.000509e0_r8,.000486e0_r8,.000516e0_r8,.000548e0_r8,.000446e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 2) = RESHAPE( (/ &
         .089988e0_r8,.045622e0_r8,.024253e0_r8,.014495e0_r8,.009303e0_r8,.005711e0_r8, &
         .003264e0_r8,.001574e0_r8,.001096e0_r8,.000762e0_r8,.000627e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .092863e0_r8,.042419e0_r8,.020704e0_r8,.012034e0_r8,.007417e0_r8,.004504e0_r8, &
         .002590e0_r8,.001334e0_r8,.000977e0_r8,.000731e0_r8,.000622e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .095737e0_r8,.039215e0_r8,.017155e0_r8,.009572e0_r8,.005532e0_r8,.003296e0_r8, &
         .001916e0_r8,.001094e0_r8,.000858e0_r8,.000699e0_r8,.000618e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .097501e0_r8,.035382e0_r8,.014856e0_r8,.008207e0_r8,.004619e0_r8,.002720e0_r8, &
         .001610e0_r8,.001012e0_r8,.000829e0_r8,.000687e0_r8,.000610e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .099264e0_r8,.031548e0_r8,.012557e0_r8,.006841e0_r8,.003705e0_r8,.002144e0_r8, &
         .001304e0_r8,.000930e0_r8,.000799e0_r8,.000675e0_r8,.000601e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .101718e0_r8,.026523e0_r8,.008473e0_r8,.004382e0_r8,.002392e0_r8,.001505e0_r8, &
         .001036e0_r8,.000836e0_r8,.000727e0_r8,.000618e0_r8,.000550e0_r8,.000494e0_r8, &
         .000501e0_r8,.000479e0_r8,.000473e0_r8,.000509e0_r8,.000541e0_r8,.000445e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 2) = RESHAPE( (/ &
         .104172e0_r8,.021499e0_r8,.004389e0_r8,.001922e0_r8,.001078e0_r8,.000865e0_r8, &
         .000767e0_r8,.000743e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104145e0_r8,.018082e0_r8,.003274e0_r8,.001493e0_r8,.000919e0_r8,.000762e0_r8, &
         .000678e0_r8,.000641e0_r8,.000584e0_r8,.000531e0_r8,.000495e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104118e0_r8,.014665e0_r8,.002159e0_r8,.001063e0_r8,.000759e0_r8,.000659e0_r8, &
         .000589e0_r8,.000539e0_r8,.000514e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .107719e0_r8,.013052e0_r8,.001822e0_r8,.000953e0_r8,.000701e0_r8,.000604e0_r8, &
         .000551e0_r8,.000525e0_r8,.000509e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .111320e0_r8,.011439e0_r8,.001485e0_r8,.000843e0_r8,.000642e0_r8,.000549e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .112375e0_r8,.011255e0_r8,.001357e0_r8,.000744e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 2) = RESHAPE( (/ &
         .109850e0_r8,.010424e0_r8,.001079e0_r8,.000567e0_r8,.000498e0_r8,.000479e0_r8, &
         .000463e0_r8,.000448e0_r8,.000418e0_r8,.000399e0_r8,.000389e0_r8,.000367e0_r8, &
         .000351e0_r8,.000328e0_r8,.000320e0_r8,.000337e0_r8,.000355e0_r8,.000304e0_r8, &
         .107002e0_r8,.009961e0_r8,.001025e0_r8,.000533e0_r8,.000497e0_r8,.000460e0_r8, &
         .000422e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107735e0_r8,.010146e0_r8,.001120e0_r8,.000576e0_r8,.000526e0_r8,.000477e0_r8, &
         .000430e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107021e0_r8,.012233e0_r8,.001533e0_r8,.000643e0_r8,.000556e0_r8,.000505e0_r8, &
         .000471e0_r8,.000448e0_r8,.000403e0_r8,.000362e0_r8,.000355e0_r8,.000296e0_r8, &
         .000251e0_r8,.000207e0_r8,.000180e0_r8,.000161e0_r8,.000152e0_r8,.000140e0_r8, &
         .106308e0_r8,.014320e0_r8,.001946e0_r8,.000709e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000473e0_r8,.000425e0_r8,.000423e0_r8,.000342e0_r8, &
         .000301e0_r8,.000257e0_r8,.000232e0_r8,.000212e0_r8,.000205e0_r8,.000209e0_r8, &
         .100592e0_r8,.015718e0_r8,.002411e0_r8,.001007e0_r8,.000802e0_r8,.000642e0_r8, &
         .000559e0_r8,.000526e0_r8,.000501e0_r8,.000474e0_r8,.000470e0_r8,.000439e0_r8, &
         .000430e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 2) = RESHAPE( (/ &
         .094877e0_r8,.017116e0_r8,.002877e0_r8,.001305e0_r8,.001018e0_r8,.000751e0_r8, &
         .000606e0_r8,.000539e0_r8,.000529e0_r8,.000524e0_r8,.000516e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .094163e0_r8,.020198e0_r8,.004594e0_r8,.001772e0_r8,.001077e0_r8,.000806e0_r8, &
         .000649e0_r8,.000565e0_r8,.000547e0_r8,.000537e0_r8,.000521e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .093449e0_r8,.023279e0_r8,.006312e0_r8,.002240e0_r8,.001135e0_r8,.000862e0_r8, &
         .000692e0_r8,.000591e0_r8,.000564e0_r8,.000549e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089886e0_r8,.026029e0_r8,.008558e0_r8,.003312e0_r8,.001655e0_r8,.001124e0_r8, &
         .000807e0_r8,.000631e0_r8,.000602e0_r8,.000568e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .086323e0_r8,.028778e0_r8,.010805e0_r8,.004383e0_r8,.002175e0_r8,.001386e0_r8, &
         .000923e0_r8,.000671e0_r8,.000640e0_r8,.000587e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082715e0_r8,.031096e0_r8,.013350e0_r8,.006131e0_r8,.003205e0_r8,.002043e0_r8, &
         .001304e0_r8,.000842e0_r8,.000734e0_r8,.000631e0_r8,.000555e0_r8,.000494e0_r8, &
         .000480e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 2) = RESHAPE( (/ &
         .079108e0_r8,.033415e0_r8,.015895e0_r8,.007878e0_r8,.004234e0_r8,.002700e0_r8, &
         .001686e0_r8,.001014e0_r8,.000829e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074807e0_r8,.034651e0_r8,.017056e0_r8,.008574e0_r8,.004769e0_r8,.002986e0_r8, &
         .001827e0_r8,.001079e0_r8,.000853e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070506e0_r8,.035887e0_r8,.018218e0_r8,.009270e0_r8,.005304e0_r8,.003271e0_r8, &
         .001969e0_r8,.001145e0_r8,.000878e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.037799e0_r8,.019680e0_r8,.009612e0_r8,.005481e0_r8,.003476e0_r8, &
         .002093e0_r8,.001123e0_r8,.000837e0_r8,.000631e0_r8,.000546e0_r8,.000447e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .064832e0_r8,.039712e0_r8,.021142e0_r8,.009954e0_r8,.005658e0_r8,.003681e0_r8, &
         .002218e0_r8,.001100e0_r8,.000796e0_r8,.000587e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063734e0_r8,.039842e0_r8,.022004e0_r8,.010859e0_r8,.005712e0_r8,.003589e0_r8, &
         .002155e0_r8,.001174e0_r8,.000856e0_r8,.000612e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 2) = (/ &
         .062636e0_r8,.039972e0_r8,.022867e0_r8,.011765e0_r8,.005766e0_r8,.003498e0_r8, &
         .002092e0_r8,.001248e0_r8,.000917e0_r8,.000637e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     3. summer
    !
    ozone(1:18, 1:6, 3) = RESHAPE( (/ &
         .059066e0_r8,.038201e0_r8,.019999e0_r8,.010813e0_r8,.006498e0_r8,.004188e0_r8, &
         .002656e0_r8,.001636e0_r8,.001110e0_r8,.000737e0_r8,.000551e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .061216e0_r8,.037359e0_r8,.019151e0_r8,.010633e0_r8,.006845e0_r8,.004382e0_r8, &
         .002691e0_r8,.001511e0_r8,.001061e0_r8,.000749e0_r8,.000568e0_r8,.000465e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063367e0_r8,.036516e0_r8,.018302e0_r8,.010452e0_r8,.007192e0_r8,.004577e0_r8, &
         .002727e0_r8,.001387e0_r8,.001012e0_r8,.000762e0_r8,.000585e0_r8,.000490e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .065518e0_r8,.035975e0_r8,.017854e0_r8,.010100e0_r8,.006534e0_r8,.003985e0_r8, &
         .002321e0_r8,.001240e0_r8,.000966e0_r8,.000774e0_r8,.000640e0_r8,.000548e0_r8, &
         .000479e0_r8,.000384e0_r8,.000346e0_r8,.000316e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.035434e0_r8,.017406e0_r8,.009749e0_r8,.005876e0_r8,.003393e0_r8, &
         .001916e0_r8,.001094e0_r8,.000920e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .069820e0_r8,.035028e0_r8,.016929e0_r8,.009389e0_r8,.005645e0_r8,.003260e0_r8, &
         .001843e0_r8,.001055e0_r8,.000905e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 3) = RESHAPE( (/ &
         .071972e0_r8,.034622e0_r8,.016453e0_r8,.009029e0_r8,.005413e0_r8,.003127e0_r8, &
         .001770e0_r8,.001015e0_r8,.000890e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .080939e0_r8,.032414e0_r8,.014163e0_r8,.007241e0_r8,.004328e0_r8,.002522e0_r8, &
         .001481e0_r8,.000934e0_r8,.000861e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .089906e0_r8,.030206e0_r8,.011873e0_r8,.005453e0_r8,.003244e0_r8,.001918e0_r8, &
         .001192e0_r8,.000852e0_r8,.000832e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .093464e0_r8,.026177e0_r8,.008525e0_r8,.003892e0_r8,.002452e0_r8,.001609e0_r8, &
         .001116e0_r8,.000851e0_r8,.000809e0_r8,.000762e0_r8,.000690e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .097023e0_r8,.022148e0_r8,.005177e0_r8,.002332e0_r8,.001660e0_r8,.001300e0_r8, &
         .001041e0_r8,.000849e0_r8,.000786e0_r8,.000737e0_r8,.000686e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .098107e0_r8,.020617e0_r8,.004758e0_r8,.002137e0_r8,.001516e0_r8,.001211e0_r8, &
         .000999e0_r8,.000848e0_r8,.000778e0_r8,.000730e0_r8,.000677e0_r8,.000603e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 3) = RESHAPE( (/ &
         .099191e0_r8,.019085e0_r8,.004340e0_r8,.001943e0_r8,.001371e0_r8,.001122e0_r8, &
         .000957e0_r8,.000847e0_r8,.000770e0_r8,.000724e0_r8,.000669e0_r8,.000601e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .099547e0_r8,.017725e0_r8,.003693e0_r8,.001578e0_r8,.001125e0_r8,.000985e0_r8, &
         .000879e0_r8,.000795e0_r8,.000712e0_r8,.000643e0_r8,.000584e0_r8,.000521e0_r8, &
         .000482e0_r8,.000384e0_r8,.000351e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .099903e0_r8,.016365e0_r8,.003045e0_r8,.001214e0_r8,.000879e0_r8,.000848e0_r8, &
         .000801e0_r8,.000744e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000441e0_r8, &
         .000410e0_r8,.000358e0_r8,.000342e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .100218e0_r8,.015035e0_r8,.002537e0_r8,.001037e0_r8,.000790e0_r8,.000726e0_r8, &
         .000673e0_r8,.000628e0_r8,.000579e0_r8,.000512e0_r8,.000440e0_r8,.000374e0_r8, &
         .000307e0_r8,.000253e0_r8,.000227e0_r8,.000208e0_r8,.000194e0_r8,.000186e0_r8, &
         .100534e0_r8,.013704e0_r8,.002028e0_r8,.000861e0_r8,.000701e0_r8,.000604e0_r8, &
         .000546e0_r8,.000513e0_r8,.000504e0_r8,.000462e0_r8,.000381e0_r8,.000307e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100892e0_r8,.012873e0_r8,.001886e0_r8,.000785e0_r8,.000643e0_r8,.000568e0_r8, &
         .000519e0_r8,.000487e0_r8,.000471e0_r8,.000437e0_r8,.000368e0_r8,.000305e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 3) = RESHAPE( (/ &
         .102665e0_r8,.010977e0_r8,.001237e0_r8,.000590e0_r8,.000498e0_r8,.000479e0_r8, &
         .000458e0_r8,.000436e0_r8,.000421e0_r8,.000387e0_r8,.000326e0_r8,.000298e0_r8, &
         .000246e0_r8,.000227e0_r8,.000211e0_r8,.000200e0_r8,.000194e0_r8,.000186e0_r8, &
         .104087e0_r8,.010726e0_r8,.000971e0_r8,.000538e0_r8,.000440e0_r8,.000434e0_r8, &
         .000418e0_r8,.000397e0_r8,.000375e0_r8,.000343e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104093e0_r8,.011539e0_r8,.001213e0_r8,.000604e0_r8,.000468e0_r8,.000442e0_r8, &
         .000414e0_r8,.000385e0_r8,.000347e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104106e0_r8,.012997e0_r8,.001479e0_r8,.000639e0_r8,.000468e0_r8,.000422e0_r8, &
         .000392e0_r8,.000372e0_r8,.000342e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104118e0_r8,.014455e0_r8,.001745e0_r8,.000674e0_r8,.000467e0_r8,.000403e0_r8, &
         .000370e0_r8,.000359e0_r8,.000337e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .103456e0_r8,.018235e0_r8,.003195e0_r8,.001379e0_r8,.000771e0_r8,.000585e0_r8, &
         .000474e0_r8,.000411e0_r8,.000380e0_r8,.000362e0_r8,.000343e0_r8,.000348e0_r8, &
         .000346e0_r8,.000328e0_r8,.000317e0_r8,.000305e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 3) = RESHAPE( (/ &
         .102794e0_r8,.022015e0_r8,.004645e0_r8,.002084e0_r8,.001076e0_r8,.000767e0_r8, &
         .000577e0_r8,.000463e0_r8,.000423e0_r8,.000399e0_r8,.000389e0_r8,.000401e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .101724e0_r8,.026500e0_r8,.007228e0_r8,.003391e0_r8,.002058e0_r8,.001285e0_r8, &
         .000811e0_r8,.000531e0_r8,.000478e0_r8,.000449e0_r8,.000440e0_r8,.000421e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .100654e0_r8,.030986e0_r8,.009812e0_r8,.004698e0_r8,.003041e0_r8,.001803e0_r8, &
         .001044e0_r8,.000599e0_r8,.000533e0_r8,.000499e0_r8,.000491e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .097097e0_r8,.034916e0_r8,.012983e0_r8,.006240e0_r8,.003666e0_r8,.002259e0_r8, &
         .001336e0_r8,.000730e0_r8,.000629e0_r8,.000549e0_r8,.000499e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .093540e0_r8,.038846e0_r8,.016155e0_r8,.007781e0_r8,.004292e0_r8,.002716e0_r8, &
         .001628e0_r8,.000862e0_r8,.000724e0_r8,.000600e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .089618e0_r8,.042869e0_r8,.019963e0_r8,.010502e0_r8,.005966e0_r8,.003525e0_r8, &
         .001936e0_r8,.000906e0_r8,.000769e0_r8,.000625e0_r8,.000508e0_r8,.000452e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 3) = RESHAPE( (/ &
         .085695e0_r8,.046892e0_r8,.023772e0_r8,.013223e0_r8,.007640e0_r8,.004333e0_r8, &
         .002244e0_r8,.000951e0_r8,.000815e0_r8,.000650e0_r8,.000508e0_r8,.000463e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082443e0_r8,.047591e0_r8,.025358e0_r8,.014294e0_r8,.008233e0_r8,.004664e0_r8, &
         .002430e0_r8,.001068e0_r8,.000851e0_r8,.000644e0_r8,.000508e0_r8,.000474e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .079190e0_r8,.048290e0_r8,.026945e0_r8,.015366e0_r8,.008826e0_r8,.004995e0_r8, &
         .002616e0_r8,.001184e0_r8,.000887e0_r8,.000637e0_r8,.000508e0_r8,.000486e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .074885e0_r8,.049987e0_r8,.030140e0_r8,.017894e0_r8,.009881e0_r8,.005543e0_r8, &
         .002907e0_r8,.001379e0_r8,.000961e0_r8,.000644e0_r8,.000512e0_r8,.000463e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .070579e0_r8,.051684e0_r8,.033335e0_r8,.020423e0_r8,.010935e0_r8,.006091e0_r8, &
         .003197e0_r8,.001573e0_r8,.001034e0_r8,.000650e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .069523e0_r8,.052249e0_r8,.034255e0_r8,.021379e0_r8,.012306e0_r8,.006727e0_r8, &
         .003415e0_r8,.001578e0_r8,.001072e0_r8,.000681e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 3) = (/ &
         .068467e0_r8,.052815e0_r8,.035175e0_r8,.022334e0_r8,.013676e0_r8,.007363e0_r8, &
         .003633e0_r8,.001582e0_r8,.001111e0_r8,.000713e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     4. fall
    !
    ozone(1:18, 1:6, 4) = RESHAPE( (/ &
         .062636e0_r8,.039972e0_r8,.022867e0_r8,.011765e0_r8,.005766e0_r8,.003498e0_r8, &
         .002092e0_r8,.001248e0_r8,.000917e0_r8,.000637e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063734e0_r8,.039842e0_r8,.022004e0_r8,.010859e0_r8,.005712e0_r8,.003589e0_r8, &
         .002155e0_r8,.001174e0_r8,.000856e0_r8,.000612e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .064832e0_r8,.039712e0_r8,.021142e0_r8,.009954e0_r8,.005658e0_r8,.003681e0_r8, &
         .002218e0_r8,.001100e0_r8,.000796e0_r8,.000587e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.037799e0_r8,.019680e0_r8,.009612e0_r8,.005481e0_r8,.003476e0_r8, &
         .002093e0_r8,.001123e0_r8,.000837e0_r8,.000631e0_r8,.000546e0_r8,.000447e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070506e0_r8,.035887e0_r8,.018218e0_r8,.009270e0_r8,.005304e0_r8,.003271e0_r8, &
         .001969e0_r8,.001145e0_r8,.000878e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074807e0_r8,.034651e0_r8,.017056e0_r8,.008574e0_r8,.004769e0_r8,.002986e0_r8, &
         .001827e0_r8,.001079e0_r8,.000853e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 4) = RESHAPE( (/ &
         .079108e0_r8,.033415e0_r8,.015895e0_r8,.007878e0_r8,.004234e0_r8,.002700e0_r8, &
         .001686e0_r8,.001014e0_r8,.000829e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .082715e0_r8,.031096e0_r8,.013350e0_r8,.006131e0_r8,.003205e0_r8,.002043e0_r8, &
         .001304e0_r8,.000842e0_r8,.000734e0_r8,.000631e0_r8,.000555e0_r8,.000494e0_r8, &
         .000480e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .086323e0_r8,.028778e0_r8,.010805e0_r8,.004383e0_r8,.002175e0_r8,.001386e0_r8, &
         .000923e0_r8,.000671e0_r8,.000640e0_r8,.000587e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089886e0_r8,.026029e0_r8,.008558e0_r8,.003312e0_r8,.001655e0_r8,.001124e0_r8, &
         .000807e0_r8,.000631e0_r8,.000602e0_r8,.000568e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .093449e0_r8,.023279e0_r8,.006312e0_r8,.002240e0_r8,.001135e0_r8,.000862e0_r8, &
         .000692e0_r8,.000591e0_r8,.000564e0_r8,.000549e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .094163e0_r8,.020198e0_r8,.004594e0_r8,.001772e0_r8,.001077e0_r8,.000806e0_r8, &
         .000649e0_r8,.000565e0_r8,.000547e0_r8,.000537e0_r8,.000521e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 4) = RESHAPE( (/ &
         .094877e0_r8,.017116e0_r8,.002877e0_r8,.001305e0_r8,.001018e0_r8,.000751e0_r8, &
         .000606e0_r8,.000539e0_r8,.000529e0_r8,.000524e0_r8,.000516e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .100592e0_r8,.015718e0_r8,.002411e0_r8,.001007e0_r8,.000802e0_r8,.000642e0_r8, &
         .000559e0_r8,.000526e0_r8,.000501e0_r8,.000474e0_r8,.000470e0_r8,.000439e0_r8, &
         .000430e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .106308e0_r8,.014320e0_r8,.001946e0_r8,.000709e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000473e0_r8,.000425e0_r8,.000423e0_r8,.000342e0_r8, &
         .000301e0_r8,.000257e0_r8,.000232e0_r8,.000212e0_r8,.000205e0_r8,.000209e0_r8, &
         .107021e0_r8,.012233e0_r8,.001533e0_r8,.000643e0_r8,.000556e0_r8,.000505e0_r8, &
         .000471e0_r8,.000448e0_r8,.000403e0_r8,.000362e0_r8,.000355e0_r8,.000296e0_r8, &
         .000251e0_r8,.000207e0_r8,.000180e0_r8,.000161e0_r8,.000152e0_r8,.000140e0_r8, &
         .107735e0_r8,.010146e0_r8,.001120e0_r8,.000576e0_r8,.000526e0_r8,.000477e0_r8, &
         .000430e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107002e0_r8,.009961e0_r8,.001025e0_r8,.000533e0_r8,.000497e0_r8,.000460e0_r8, &
         .000422e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 4) = RESHAPE( (/ &
         .109850e0_r8,.010424e0_r8,.001079e0_r8,.000567e0_r8,.000498e0_r8,.000479e0_r8, &
         .000463e0_r8,.000448e0_r8,.000418e0_r8,.000399e0_r8,.000389e0_r8,.000367e0_r8, &
         .000351e0_r8,.000328e0_r8,.000320e0_r8,.000337e0_r8,.000355e0_r8,.000304e0_r8, &
         .112375e0_r8,.011255e0_r8,.001357e0_r8,.000744e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .111320e0_r8,.011439e0_r8,.001485e0_r8,.000843e0_r8,.000642e0_r8,.000549e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .107719e0_r8,.013052e0_r8,.001822e0_r8,.000953e0_r8,.000701e0_r8,.000604e0_r8, &
         .000551e0_r8,.000525e0_r8,.000509e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104118e0_r8,.014665e0_r8,.002159e0_r8,.001063e0_r8,.000759e0_r8,.000659e0_r8, &
         .000589e0_r8,.000539e0_r8,.000514e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104145e0_r8,.018082e0_r8,.003274e0_r8,.001493e0_r8,.000919e0_r8,.000762e0_r8, &
         .000678e0_r8,.000641e0_r8,.000584e0_r8,.000531e0_r8,.000495e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 4) = RESHAPE( (/ &
         .104172e0_r8,.021499e0_r8,.004389e0_r8,.001922e0_r8,.001078e0_r8,.000865e0_r8, &
         .000767e0_r8,.000743e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .101718e0_r8,.026523e0_r8,.008473e0_r8,.004382e0_r8,.002392e0_r8,.001505e0_r8, &
         .001036e0_r8,.000836e0_r8,.000727e0_r8,.000618e0_r8,.000550e0_r8,.000494e0_r8, &
         .000501e0_r8,.000479e0_r8,.000473e0_r8,.000509e0_r8,.000541e0_r8,.000445e0_r8, &
         .099264e0_r8,.031548e0_r8,.012557e0_r8,.006841e0_r8,.003705e0_r8,.002144e0_r8, &
         .001304e0_r8,.000930e0_r8,.000799e0_r8,.000675e0_r8,.000601e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .097501e0_r8,.035382e0_r8,.014856e0_r8,.008207e0_r8,.004619e0_r8,.002720e0_r8, &
         .001610e0_r8,.001012e0_r8,.000829e0_r8,.000687e0_r8,.000610e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .095737e0_r8,.039215e0_r8,.017155e0_r8,.009572e0_r8,.005532e0_r8,.003296e0_r8, &
         .001916e0_r8,.001094e0_r8,.000858e0_r8,.000699e0_r8,.000618e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .092863e0_r8,.042419e0_r8,.020704e0_r8,.012034e0_r8,.007417e0_r8,.004504e0_r8, &
         .002590e0_r8,.001334e0_r8,.000977e0_r8,.000731e0_r8,.000622e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 4) = RESHAPE( (/ &
         .089988e0_r8,.045622e0_r8,.024253e0_r8,.014495e0_r8,.009303e0_r8,.005711e0_r8, &
         .003264e0_r8,.001574e0_r8,.001096e0_r8,.000762e0_r8,.000627e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .084591e0_r8,.046553e0_r8,.025037e0_r8,.015156e0_r8,.009841e0_r8,.006124e0_r8, &
         .003534e0_r8,.001693e0_r8,.001170e0_r8,.000793e0_r8,.000631e0_r8,.000537e0_r8, &
         .000551e0_r8,.000509e0_r8,.000486e0_r8,.000516e0_r8,.000548e0_r8,.000446e0_r8, &
         .079194e0_r8,.047483e0_r8,.025822e0_r8,.015818e0_r8,.010380e0_r8,.006537e0_r8, &
         .003804e0_r8,.001811e0_r8,.001245e0_r8,.000825e0_r8,.000635e0_r8,.000570e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .077409e0_r8,.048159e0_r8,.026661e0_r8,.016596e0_r8,.010962e0_r8,.006972e0_r8, &
         .004160e0_r8,.002132e0_r8,.001391e0_r8,.000868e0_r8,.000686e0_r8,.000601e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .075625e0_r8,.048835e0_r8,.027500e0_r8,.017375e0_r8,.011544e0_r8,.007407e0_r8, &
         .004516e0_r8,.002453e0_r8,.001536e0_r8,.000912e0_r8,.000737e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .074927e0_r8,.049459e0_r8,.029215e0_r8,.018025e0_r8,.011754e0_r8,.007786e0_r8, &
         .004972e0_r8,.002926e0_r8,.001817e0_r8,.001025e0_r8,.000758e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 4) = (/ &
         .074229e0_r8,.050084e0_r8,.030930e0_r8,.018676e0_r8,.011965e0_r8,.008165e0_r8, &
         .005428e0_r8,.003399e0_r8,.002098e0_r8,.001138e0_r8,.000780e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8/)

    ozsig(:) = (/ &
         .020747_r8,.073986_r8,.124402_r8,.174576_r8,.224668_r8,.274735_r8, &
         .324767_r8,.374806_r8,.424818_r8,.497450_r8,.593540_r8,.688125_r8, &
         .777224_r8,.856317_r8,.920400_r8,.960480_r8,.981488_r8,.995004_r8/)


    first_getoz=.TRUE.

    IF(first_getoz)THEN
       mon_getoz=INT(yrl/12.0_r8)
       year_getoz=yrl
!      IF(nlm_getoz.NE.kmax)THEN
!         inter_getoz=.TRUE.
!      ELSE
!         inter_getoz=.FALSE.
!         DO l=1,nlm_getoz
!            ll=nlm_getoz-l+1
!            IF(ABS(ozsig(l)-sl(ll)).GT.0.0001_r8)inter_getoz=.TRUE.
!         END DO
!      ENDIF
       inter_getoz=.TRUE.
       first_getoz=.FALSE.
    ENDIF
  END SUBROUTINE InitGetoz



 SUBROUTINE RadiationDriver(       &
      ! Run Flags
      first  , ifday, lcnvl , lthncl, nfin0 , nfin1 , nfcnv0,  &
      intcosz, kt   , mxrdcc,                                  &
      ! Time info
      yrl   , idatec , idate , tod   , jdt   , delt  ,         &
      trint , swint  ,                                         &
      ! Model Geometry
      colrad, lonrad, zenith, cos2d ,                          &
      ! Model information
      latco , ncols , kmax  , nls   , nlcs  , imask ,          &
      ! Atmospheric fields
      prsi ,prsl    ,phii   ,phil    ,&
      gps   , gtt   , gqq   , tsurf , omg   , tsea  ,          &
      QCF   ,QCL    , QCR   ,                                  &
      ! CONVECTION: convective clouds
      convts, convcs, convbs, convc , convt , convb ,          &
      ! SURFACE:  albedo
      AlbVisDiff, AlbNirDiff , AlbVisBeam , AlbNirBeam ,       &
      ! SW Radiation fields at last integer hour
      rSwToaDown,                                              &
      rVisDiff  , rNirDiff   , rVisBeam   , rNirBeam ,         &
      rVisDiffC , rNirDiffC  , rVisBeamC  , rNirBeamC,         &
      rSwSfcNet , rSwSfcNetC , SwSfcUp ,&
      ! SW Radiation fields at next integer hour
      ySwToaDown,                                              &
      yVisDiff  , yNirDiff   , yVisBeam   , yNirBeam ,         &
      yVisDiffC , yNirDiffC  , yVisBeamC  , yNirBeamC,         &
      ySwHeatRate, ySwHeatRateC,                               &
      ySwSfcNet , ySwSfcNetC , &
      ! Radiation field (Interpolated) at time = tod
      xVisDiff  , xNirDiff   , xVisBeam   , xNirBeam ,         &
      ! LW Radiation fields at last integer hour
      LwCoolRate, LwSfcDown  , LwSfcNet   , LwToaUp  ,         &
      LwCoolRateC, LwSfcDownC, LwSfcNetC  , LwToaUpC ,         &
      ! SSIB: Total radiation absorbed at ground
      slrad ,                                                  &
      ! SSIB INIT: Solar radiation with cos2
      ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff,  &
      ! Cloud field
      cldsav, CldCovTot,                                       &
      CldCovInv, CldCovSat, CldCovCon, CldCovSha,              &
      ! Microphysics
      CldLiqWatPath  , emisd , taud  , EFFCS    ,EFFIS  ,      &
      ! Chemistry
      o3mix  ,co2m ,dump,CLDF,&
!tar begin
!climate aerosol optical parameters of coarse mode
      aod,asy,ssa,z_aer,topog, &
!tar end 
!
!tar begin
!climate aerosol optical parameters of coarse mode
      aodF,asyF,ssaF,z_aerF)
!tar end 
!    
   IMPLICIT NONE
    !==========================================================================
    !
    ! _________
    ! RUN FLAGS
    !
    ! first.......control logical variable .true. or .false.
    ! ifday.......model forecast day
    ! lcnvl.......the lowest layer index where non-convective clouds can
    !             occur (ben says this should be 2 or more)
    !             constant lcnvl = 2
    ! lthncl......Minimum depth in mb of non-zero low level cloud
    !             consta lthncl=80
    ! nfin0.......input  file at time level t-dt
    ! nfin1.......input  file at time level t
    ! nfcnv0......initial information on convective clouds for int. radiation
    ! mxrdcc......use maximum random converage for radiative conv. clouds
    !             constant logical mxrdcc = .true.
    ! kt..........hour of present  time step
    !
    ! _________
    ! TIME INFO
    !
    ! yrl.........length of year in days
    ! idatec(4)...output : idatec(1)= current hour of day
    !                      idatec(2)= current day of month.
    !                      idatec(3)= current month of year.
    !                      idatec(4)= current year.
    ! idate(4)....output : idate(1) = initial hour of day
    !                      idate(2) = day of month.(???)
    !                      idate(3) = month of year.(???)
    !                      idate(4) = year.
    ! tod.........model forecast time of day in seconds
    ! jdt.........time step in getdia
    ! delt........time interval in sec (fixed throuh the integration)
    ! swint.......sw subr. call interval in hours
    !             swint has to be less than or equal to trint
    !                              and mod(trint,swint)=0
    ! trint.......ir subr. call interval in hours
    !
    ! ______________
    ! MODEL GEOMETRY
    !
    ! colrad.....colatitude  colrad=0-3.14 from np to sp in radians
    ! lonrad.....longitude in radians
    ! zenith.....cosine of solar zenith angle
    ! cos2
    ! sig........sigma coordinate at middle of layer
    ! sigml......sigma coordinate at bottom of layer
    ! delsig      k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(2)
    !             k=3/2----sigml,ric,rf,km,kh,b,l -----------
    !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(1)
    !             k=1/2----sigml ----------------------------
    !
    ! _________________
    ! MODEL INFORMATION
    !
    ! latco.......latitude
    ! ncols.....Number of grid points on a gaussian latitude circle
    ! kmax......Number of sigma levels
    ! nls..... .Number of layers in the stratosphere.
    ! nlcs......nlcs =   30
    ! imask......mascara continental
    !tar begin
    ! topog......topography field
    !tar end    
    !
    ! __________________
    ! ATMOSPHERIC FIELDS
    !
    ! gps.......Surface pressure in mb
    ! gtt.........gtt =  gtmp(imx,kmax) input  : temperature.
    ! gqq.........gqq = gq(imx,kmax)     input : specific humidity.
    ! tsurf
    ! omg.........omg   =  vertical velocity  (cb/sec)
    ! tsea
    !
    ! _____________________________
    ! CONVECTION: convective clouds
    !
    ! convc.......ncols convective cloud cover in 3 hr. avrage
    ! convt.......ncols convective cloud top  (sigma layer)
    ! convb.......ncols convective cloud base (sigma layer)
    ! convts
    ! convcs
    ! convbs
    !
    ! ________________
    ! SURFACE:  albedo
    !
    ! AlbVisBeam.......visible beam surface albedo
    ! AlbVisDiff.......visible diffuse surface albedo
    ! AlbNirBeam.......near-ir beam surface albedo
    ! AlbNirDiff.......near-ir diffuse surface albedo
    !
    ! ________________________________________
    ! SW Radiation fields at last integer hour
    !
    ! rSwToaDown......Incident SW at top (W/m^2)                
    ! rVisBeam....... Down Sfc SW flux visible beam    (all-sky)
    ! rVisDiff....... Down Sfc SW flux visible diffuse (all-sky)
    ! rNirBeam....... Down Sfc SW flux Near-IR beam    (all-sky)
    ! rNirDiff....... Down Sfc SW flux Near-IR diffuse (all-sky)
    ! rVisBeamC...... Down Sfc SW flux visible beam    (clear)  
    ! rVisDiffC...... Down Sfc SW flux visible diffuse (clear)  
    ! rNirBeamC...... Down Sfc SW flux Near-IR beam    (clear)  
    ! rNirDiffC...... Down Sfc SW flux Near-IR diffuse (clear)  
    !
    ! ________________________________________
    ! SW Radiation fields at next integer hour
    !
    ! ySwToaDown.....Incident SW at top 
    ! yVisBeam.......Down Sfc SW flux visible beam    (all-sky)
    ! yVisDiff.......Down Sfc SW flux visible diffuse (all-sky)
    ! yNirBeam.......Down Sfc SW flux Near-IR beam    (all-sky)
    ! yNirDiff.......Down Sfc SW flux Near-IR diffuse (all-sky)
    ! yVisBeamC......Down Sfc SW flux visible beam    (clear)  
    ! yVisDiffC......Down Sfc SW flux visible diffuse (clear)  
    ! yNirBeamC......Down Sfc SW flux Near-IR beam    (clear)  
    ! yNirDiffC......Down Sfc SW flux Near-IR diffuse (clear)  
    ! ySwHeatRate....Heating rate due to shortwave         (K/s)
    ! ySwHeatRateC...Heating rate due to shortwave (clear) (K/s)
    !
    ! ____________________________________________
    ! Radiation field (Interpolated) at time = tod
    !
    ! xVisBeam.......Down Sfc SW flux visible beam    (all-sky)
    ! xVisDiff.......Down Sfc SW flux visible diffuse (all-sky)
    ! xNirBeam.......Down Sfc SW flux Near-IR beam    (all-sky)
    ! xNirDiff.......Down Sfc SW flux Near-IR diffuse (all-sky)
    !
    ! ________________________________________
    ! LW Radiation fields at next integer 3-hour
    !
    ! LwCoolRate.....Cooling rate due to longwave  (all-sky) (K/s)
    ! LwSfcDown......Down Sfc LW flux              (all-sky) (W/m2)
    ! LwSfcNet.......Net Sfc LW                    (all-sky) (W/m2)
    ! LwToaUp........Longwave upward at top        (all-sky) (W/m2) 
    ! LwCoolRateC....Cooling rate due to longwave  (clear) (K/s)
    ! LwSfcDownC.....Down Sfc LW flux              (clear) (W/m2)
    ! LwSfcNetC......Net Sfc LW                    (clear) (W/m2)
    ! LwToaUpC.......Longwave upward at top        (clear) (W/m2)
    !
    ! ________________________________________
    ! SSIB: Total radiation absorbed at ground
    !
    ! slrad..........Total radiation absorbed at ground
    !
    ! ____________________________________
    ! SSIB INIT: Solar radiation with cos2
    !
    ! ssib_VisBeam.......Down Sfc SW flux visible beam    (all-sky)
    ! ssib_VisDiff.......Down Sfc SW flux visible diffuse (all-sky)
    ! ssib_NirBeam.......Down Sfc SW flux Near-IR beam    (all-sky)
    ! ssib_NirDiff.......Down Sfc SW flux Near-IR diffuse (all-sky)
    !
    ! ___________
    ! Cloud field
    !
    ! cldsav......Cloud cover
    ! CldCovTot......Total cloud cover (at each layer)
    ! CldCovInv......Inversion clouds                 
    ! CldCovSat......Saturation clouds                
    ! CldCovCon......Convection clouds                
    ! CldCovSha......Shallow convective clouds        
    !
    ! ____________
    ! Microphysics (if cld=3 or 4)
    !
    ! CldLiqWatPath...Cloud liquid water path (parametrized)
    ! emisd...........Emissivity 
    ! taud............Shortwave cloud optical depth
    !
    ! _________
    ! Chemistry
    !
    ! co2m......co2m is wgne standard value in mol/mol "co2m = /345.0e-6/
    ! o3mix
    !
    !tar begin
    !
    ! aod... aerosol optical depth of coarse mode
    ! asy... asymmetry factor of coarse mode
    ! ssa... single scattering albedo of coarse mode
    ! z_aer...aod vertical distribution
    !
    !tar end
    !
    !tar begin
    !
    ! aodF... aerosol optical depth of fine mode
    ! asyF... asymmetry factor of fine mode
    ! ssaF... single scattering albedo of fine mode
    ! z_aerF...aod vertical distribution of fine mode
    !
    !tar end    
    !==========================================================================

    ! Run Flags
    LOGICAL      ,    INTENT(in   ) :: first
    INTEGER      ,    INTENT(in   ) :: ifday
    INTEGER      ,    INTENT(in   ) :: lcnvl
    INTEGER      ,    INTENT(in   ) :: lthncl
    INTEGER      ,    INTENT(in   ) :: nfin0
    INTEGER      ,    INTENT(in   ) :: nfin1
    INTEGER      ,    INTENT(in   ) :: nfcnv0
    LOGICAL      ,    INTENT(IN   ) :: intcosz
    INTEGER      ,    INTENT(IN   ) :: kt
    LOGICAL      ,    INTENT(in   ) :: mxrdcc

    ! Time info
    REAL(KIND=r8),    INTENT(in   ) :: yrl
    INTEGER      ,    INTENT(in   ) :: idatec(4)
    INTEGER      ,    INTENT(in   ) :: idate(4)
    REAL(KIND=r8),    INTENT(in   ) :: tod
    INTEGER      ,    INTENT(IN   ) :: jdt
    REAL(KIND=r8),    INTENT(in   ) :: delt
    REAL(KIND=r8),    INTENT(in   ) :: trint
    REAL(KIND=r8),    INTENT(in   ) :: swint

    ! Model Geometry
    REAL(KIND=r8),    INTENT(in   ) :: colrad(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: lonrad(ncols)
    REAL(KIND=r8),    INTENT(inout) :: zenith  (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: cos2d (ncols)

    ! Model information
    INTEGER         , INTENT(in   ) :: latco
    INTEGER         , INTENT(IN   ) :: ncols
    INTEGER         , INTENT(IN   ) :: kmax
    INTEGER         , INTENT(IN   ) :: nls
    INTEGER         , INTENT(IN   ) :: nlcs
    INTEGER(KIND=i8), INTENT(IN   ) :: imask (ncols)
!tar begin
    REAL(KIND=r8),    INTENT(IN   ) :: topog (ncols)
!tar end    

    ! Atmospheric fields
    REAL(KIND=r8),    INTENT(in   ) :: prsi  (ncols,kMax+1)
    REAL(KIND=r8),    INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: phii  (nCols,kMax+1)
    REAL(KIND=r8),    INTENT(in   ) :: phil  (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gps   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gtt   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gqq   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tsurf (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: omg   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: QCF   (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: QCL   (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: QCR   (ncols,kMax)
    ! CONVECTION: convective clouds
    REAL(KIND=r8),    INTENT(in   ) :: convts(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convcs(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convbs(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convt (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convb (ncols)

    ! SURFACE:  albedo
    REAL(KIND=r8),    INTENT(in   ) :: AlbVisDiff (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: AlbNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: AlbVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: AlbNirBeam (ncols)

    ! SW Radiation fields at last integer hour
    REAL(KIND=r8),    INTENT(inout) :: rSwToaDown(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rSwSfcNet   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rSwSfcNetC  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: SwSfcUp (ncols)
    ! SW Radiation fields at next integer hour
    REAL(KIND=r8),    INTENT(inout) :: ySwToaDown(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: ySwHeatRate   (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: ySwHeatRateC  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: ySwSfcNet   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ySwSfcNetC  (ncols)

    ! Radiation field (Interpolated) at time = tod
    REAL(KIND=r8),    INTENT(inout) :: xVisDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: xNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: xVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: xNirBeam (ncols)

    ! LW Radiation fields at last integer hour
    REAL(KIND=r8),    INTENT(inout) :: LwCoolRate (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcDown  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcNet   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwToaUp    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwCoolRateC(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcDownC (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcNetC  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwToaUpC   (ncols)

    ! SSIB: Total radiation absorbed at ground
    REAL(KIND=r8),    INTENT(  OUT)  :: slrad(ncols)

    ! SSIB INIT: Solar radiation with cos2
    REAL(KIND=r8),    INTENT(  OUT)  :: ssib_VisBeam (ncols)
    REAL(KIND=r8),    INTENT(  OUT)  :: ssib_VisDiff (ncols)
    REAL(KIND=r8),    INTENT(  OUT)  :: ssib_NirBeam (ncols)
    REAL(KIND=r8),    INTENT(  OUT)  :: ssib_NirDiff (ncols)

    ! Cloud field
    REAL(KIND=r8),    INTENT(inout) :: cldsav(ncols)
    REAL(KIND=r8),    INTENT(inout) :: CldCovTot(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovInv(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovSat(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovCon(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovSha(ncols,kmax)

    ! Microphysics
    REAL(KIND=r8),    INTENT(inout) :: CldLiqWatPath  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: emisd (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: taud  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: EFFCS (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: EFFIS (ncols,kmax)
    REAL(KIND=r8) :: rei   (ncols,kmax)
    REAL(KIND=r8) :: rel   (ncols,kmax)

    ! Chemistry
    REAL(KIND=r8),    INTENT(INOUT) :: o3mix(ncols,kMax)   
    REAL(KIND=r8),    INTENT(INOUT) :: co2m(ncols,kMax)   !mol/mol
    REAL(KIND=r8),    INTENT(INOUT) :: dump(ncols,kMax) 
    REAL(KIND=r8),    INTENT(INOUT) :: CLDF(ncols,kMax)
!tar begin  
! climate aerosol optical parameters of coarse mode
!
    REAL(KIND=r8),    INTENT(inout) :: aod(ncols,14)
    REAL(KIND=r8),    INTENT(inout) :: asy(ncols,14)    
    REAL(KIND=r8),    INTENT(inout) :: ssa(ncols,14)    
    REAL(KIND=r8),    INTENT(inout) :: z_aer(ncols,40)    
!        
!tar end 
!
!tar begin  
! climate aerosol optical parameters of fine mode
!
    REAL(KIND=r8),    INTENT(inout) :: aodF(ncols,14)
    REAL(KIND=r8),    INTENT(inout) :: asyF(ncols,14)    
    REAL(KIND=r8),    INTENT(inout) :: ssaF(ncols,14)    
    REAL(KIND=r8),    INTENT(inout) :: z_aerF(ncols,40)    
!        
!tar end 
    !==========================================================================
    ! LOCAL VARIABLES
    !==========================================================================

    ! Radiation field (Interpolated) at time = tod
    REAL(KIND=r8) :: xSwToaDown(ncols)
    REAL(KIND=r8) :: xVisDiffC(ncols)
    REAL(KIND=r8) :: xNirDiffC(ncols)
    REAL(KIND=r8) :: xVisBeamC(ncols)
    REAL(KIND=r8) :: xNirBeamC(ncols)

    ! Interpolation weights
    REAL(KIND=r8)    :: fstp
    REAL(KIND=r8)    :: fstp1

    ! Some diagnostic of radiation
    REAL(KIND=r8) :: SwToaUp    (ncols) ! Upward Shortwave Flux at TOA        
    REAL(KIND=r8) :: SwToaUpC   (ncols) ! Upward Shortwave Flux at TOA (clear)
!    REAL(KIND=r8) :: SwSfcDown  (ncols) ! Downward Shortwave Flux at Surface
!    REAL(KIND=r8) :: SwSfcDownC (ncols) ! Downward Shortwave Flux at Surface (clear)

    REAL(KIND=r8) :: xSwSfcNet   (ncols) ! Shortwave Flux Absorbed at Surface
    REAL(KIND=r8) :: xSwSfcNetC  (ncols) ! Shortwave Flux Absorbed at Surface (clear)

    REAL(KIND=r8) :: rSwToaNet   (ncols) ! Shortwave Flux Absorbed at Surface
    REAL(KIND=r8) :: rSwToaNetC  (ncols) ! Shortwave Flux Absorbed at Surface (clear)

    REAL(KIND=r8) :: ySwToaNet   (ncols) ! Shortwave Flux Absorbed at Surface
    REAL(KIND=r8) :: ySwToaNetC  (ncols) ! Shortwave Flux Absorbed at Surface (clear)
    REAL(KIND=r8) :: ySwAtmAbs    (ncols) 
    REAL(KIND=r8) :: ySwAtmAbsC   (ncols) 
    REAL(KIND=r8) :: DeltaP (ncols,kmax)
    REAL(KIND=r8) :: gtmp   (ncols,kmax)
    REAL(KIND=r8) :: q      (ncols,kmax)
    REAL(KIND=r8) :: relhum (ncols,kmax)
    REAL(KIND=r8) :: CldLow (ncols) 
    REAL(KIND=r8) :: CldMed (ncols) 
    REAL(KIND=r8) :: CldHgh (ncols) 
    INTEGER :: modstp
    INTEGER :: intstp    
    LOGICAL :: ghl_local
    REAL(KIND=r8)    :: rfac 
    INTEGER :: i   
    INTEGER :: k   
    INTEGER :: icld
    CHARACTER(LEN=*), PARAMETER     :: h='**(RadiationDriver)**'

    ! inalb.......inalb Input two types of surface albedo
    !                >>> inalb= 1 : input two  types surfc albedo (2 diffused)
    !                             direct beam albedos are calculated by the subr.
    !                >>> inalb= 2 : input four types surfc albedo (2 diff,2 direct)

    INTEGER       :: inalb
 

    ghl_local = IsGridHistoryOn()

    DO i = 1,ncols
       CldLow (i)= 0.0_r8
       CldMed (i)= 0.0_r8
       CldHgh (i)= 0.0_r8
    END DO

    DO k=1,kmax
       DO i = 1,ncols
         DeltaP(i,k) = ((prsi(i,k)) - (prsi(i,k+1)))/prsi(i,1)
       END DO
    END DO

    !
    !     virtual temperature correction for temperature used in radiation
    !     and setting minimum specific humidity
    !
!    CALL rqvirt(ncols , kmax, sig, gps, gtt, gqq, gtmp, q, relhum )
    CALL rqvirt(ncols , kmax,prsl, gtt, gqq, gtmp, q, relhum )
 
    ! >>> inalb= 1    : input two  types surfc albedo (2 diffused)
    !                   direct beam albedos are calculated by the subr.
    ! >>> inalb= 2    : input four types surfc albedo (2 diff,2 direct)
    inalb=2
    IF(TRIM(iswrad).NE.'NON'.OR.TRIM(ilwrad).NE.'NON') THEN    

       icld=1

!       IF(TRIM(iccon).EQ.'ARA') icld=INT(crdcld) 

       IF(TRIM(iswrad).EQ.'CRD'.OR.TRIM(ilwrad).EQ.'CRD') icld=INT(crdcld)

       IF(TRIM(iswrad).EQ.'CRDTF'.OR.TRIM(ilwrad).EQ.'CRDTF') icld=INT(crdcld)

       IF(TRIM(iswrad).EQ.'UKM'.OR.TRIM(ilwrad).EQ.'UKM')icld=INT(crdcld)

       IF(TRIM(iswrad).EQ.'RRTMG'.OR.TRIM(ilwrad).EQ.'RRTMG')icld=INT(crdcld)
       !CALL MsgOne(h,'run spmrad')

       CALL spmrad (&
            ! Run Flags
            first  , ifday , lcnvl , lthncl, nfin0 , nfin1 , nfcnv0, &
            intcosz, mxrdcc, inalb , icld  ,                         &
            ! Time info
            yrl    , idatec, idate , tod   , delt  ,trint  , swint , &
            ! Model Geometry
            colrad , lonrad, zenith, cos2d , &
            ! Model information
!tar begin
            !ncols , kmax  , nls   , nlcs  , imask ,         &
            latco, ncols , kmax  , nls   , nlcs  , imask ,         &
!tar end
            ! Atmospheric fields
            prsi   ,prsl   ,phii   ,phil    ,&
            gtmp  , q     , tsurf , relhum, omg   , tsea  , &
            QCF    ,QCL    , QCR   ,                                 &
            ! CONVECTION: convective clouds
            convts , convcs, convbs, convc , convt , convb ,         &
            ! SURFACE:  albedo
            AlbVisDiff , AlbNirDiff , AlbVisBeam , AlbNirBeam ,      &
            ! SW Radiation fields at last integer hour
            rSwToaDown,                                              &
            rVisDiff   , rNirDiff   , rVisBeam   , rNirBeam   ,      &
            rVisDiffC  , rNirDiffC  , rVisBeamC  , rNirBeamC  ,      &
            rSwSfcNet  , rSwSfcNetC , rSwToaNet  , rSwToaNetC ,      &
            ! SW Radiation fields at next integer hour
            ySwToaDown,                                              &
            yVisDiff   , yNirDiff   , yVisBeam   , yNirBeam   ,      &
            yVisDiffC  , yNirDiffC  , yVisBeamC  , yNirBeamC  ,      &
            ySwHeatRate, ySwHeatRateC,                               &
            ySwSfcNet  , ySwSfcNetC , ySwToaNet  , ySwToaNetC ,      &
            ! SSIB INIT: Solar radiation with cos2
            ssib_VisDiff, ssib_NirDiff, ssib_VisBeam, ssib_NirBeam,  &
            ! LW Radiation fields at last integer hour
            LwCoolRate , LwSfcDown  , LwSfcNet   , LwToaUp    ,      &
            LwCoolRateC, LwSfcDownC , LwSfcNetC  , LwToaUpC   ,      &
            ! Cloud field
            CldLow   , CldMed   , CldHgh   , cldsav   , CldCovTot,   &
            CldCovInv, CldCovSat, CldCovCon, CldCovSha,              &
            ! Microphysics
            CldLiqWatPath , emisd , taud  , rei       ,rel    ,      &
            EFFCS ,EFFIS  ,        &
            ! Chemistry
            o3mix ,co2m,dump,CLDF,&
!tar begin
! climate aerosol optical parameters of coarse mode
            aod,asy,ssa,z_aer,topog, &
!tar end
!
!tar begin
! climate aerosol optical parameters of fine mode
            aodF,asyF,ssaF,z_aerF)
!tar end
!
       !CALL MsgOne(h,'end spmrad')

       !
       !     this is for radiation interpolation
       !
       !PK intstp=INT(3600.0_r8*swint/delt+0.01_r8)
       intstp=INT(swint/delt+0.01_r8)
       modstp=MOD(jdt-1,intstp)
       fstp  =REAL(modstp,r8)/REAL(intstp,r8)

       IF(jdt.LE.2.AND.kt.EQ.0) fstp=0.0_r8

       fstp1 =1.0_r8-fstp

       DO i=1,ncols
          xVisDiff (i)=fstp1*rVisDiff (i)+fstp*yVisDiff (i)
          xVisBeam (i)=fstp1*rVisBeam (i)+fstp*yVisBeam (i)
          xNirDiff (i)=fstp1*rNirDiff (i)+fstp*yNirDiff (i)
          xNirBeam (i)=fstp1*rNirBeam (i)+fstp*yNirBeam (i)
          xVisDiffC(i)=fstp1*rVisDiffC(i)+fstp*yVisDiffC(i)
          xVisBeamC(i)=fstp1*rVisBeamC(i)+fstp*yVisBeamC(i)
          xNirDiffC(i)=fstp1*rNirDiffC(i)+fstp*yNirDiffC(i)
          xNirBeamC(i)=fstp1*rNirBeamC(i)+fstp*yNirBeamC(i)
          xSwToaDown(i)=fstp1*rSwToaDown(i)+fstp*ySwToaDown(i)
          xSwSfcNet (i)=fstp1*rSwSfcNet (i)+fstp*ySwSfcNet (i)
          xSwSfcNetC(i)=fstp1*rSwSfcNetC(i)+fstp*ySwSfcNetC(i)
       END DO

       IF(modstp.EQ.intstp-1) THEN
          DO i=1,ncols
             rVisDiff (i)=yVisDiff (i)
             rVisBeam (i)=yVisBeam (i)
             rNirDiff (i)=yNirDiff (i)
             rNirBeam (i)=yNirBeam (i)
             rVisDiffC(i)=yVisDiffC(i)
             rVisBeamC(i)=yVisBeamC(i)
             rNirDiffC(i)=yNirDiffC(i)
             rNirBeamC(i)=yNirBeamC(i)
             rSwToaDown(i)=ySwToaDown(i)
             rSwSfcNet (i)=ySwSfcNet (i)
             rSwSfcNetC(i)=ySwSfcNetC(i)
          END DO
       END IF

       ! SSIB: Total radiation absorbed at ground 
       DO i=1,ncols
           SwSfcUp(i)=AlbVisDiff(i)*xVisDiff(i)+AlbVisBeam(i)*xVisBeam(i)+ &
                      AlbNirDiff(i)*xNirDiff(i)+AlbNirBeam(i)*xNirBeam(i)

          slrad(i) = LwSfcDown(i) & !+ SwSfcNet(i) 
               +(1.0_r8-AlbVisDiff(i))*xVisDiff(i) &
               +(1.0_r8-AlbVisBeam(i))*xVisBeam(i) &
               +(1.0_r8-AlbNirDiff(i))*xNirDiff(i) &
               +(1.0_r8-AlbNirBeam(i))*xNirBeam(i)

          slrad(i) = -14.3353e-04_r8*slrad(i)
       END DO
       
       ! essas variaveis sao locais, ou seja, sao calculadas aqui apenas para diagnostico!
       ySwAtmAbs=0.0_r8
       ySwAtmAbsC=0.0_r8
       DO k=1,kmax
          !rfac =cp*100.0_r8*delsig(k)/grav
          DO i=1,ncols
             rfac =cp*100.0_r8*DeltaP(i,k)/grav
             ySwAtmAbs (i)=ySwAtmAbs (i)+rfac*gps(i)*ySwHeatRate (i,k)
             ySwAtmAbsC(i)=ySwAtmAbsC(i)+rfac*gps(i)*ySwHeatRateC(i,k)
          END DO
       END DO

!     IF(dodia(nDiag_swutop)) CALL updia(ySwAtmAbs ,nDiag_swutop,latco)
!     IF(dodia(nDiag_swutpc)) CALL updia(ySwAtmAbsC,nDiag_swutpc,latco)

! isso parece errado pois os campos "x" foram interpolados para o timestep atual
! ... mas heating rate, usado acima para calcular a absorcao na atmosfera,
! permaneceu constante!!!
       DO i=1,ncols
          SwToaUp(i)=xSwToaDown(i) &
!               -SwSfcNet(i) &
               -(1.0_r8-AlbVisDiff(i))*xVisDiff(i) &
               -(1.0_r8-AlbVisBeam(i))*xVisBeam(i) &
               -(1.0_r8-AlbNirDiff(i))*xNirDiff(i) &
               -(1.0_r8-AlbNirBeam(i))*xNirBeam(i) &
               -ySwAtmAbs(i)
          SwToaUpC(i)=xSwToaDown(i) &
!               -SwSfcNetC(i) &
               -(1.0_r8-AlbVisDiff(i))*xVisDiffC(i) &
               -(1.0_r8-AlbVisBeam(i))*xVisBeamC(i) &
               -(1.0_r8-AlbNirDiff(i))*xNirDiffC(i) &
               -(1.0_r8-AlbNirBeam(i))*xNirBeamC(i) &
               -ySwAtmAbsC(i)
       END DO
       !
       ! Storage Diagnostics of Radiation  
       !
       IF( StartStorDiag ) THEN          
          CALL RadDiagStor(ncols,kmax,latco,&
               CldCovTot,CldCovInv,CldCovSat,CldCovCon,CldCovSha,&
               CldLiqWatPath,taud, rei       ,rel   ,emisd, &
               LwCoolRate,ySwHeatRate,ySwHeatRateC,LwCoolRateC,&
               cldsav,o3mix,co2m,&
               xSwToaDown,SwToaUp,SwToaUpC,LwSfcDown,LwToaUp,LwSfcDownC,LwToaUpC,LwSfcNet,LwSfcNetC,xVisDiff,& 
               xVisBeam,xNirDiff,xNirBeam,xVisDiffC,xVisBeamC,xNirDiffC,xNirBeamC,&
               AlbVisDiff,AlbVisBeam,AlbNirDiff,AlbNirBeam, &
               DeltaP, gps,dump)
!, SwSfcDown, SwSfcDownC  )

       END IF
       !
       !Storage GridHistory of Radiation 
       !    
       IF (ghl_local) THEN
          CALL RadGridHistStorage (ncols,kmax,latco,CldCovTot,LwCoolRate,ySwHeatRate,LwToaUp,&
                                LwSfcDown,zenith,xSwToaDown,xVisBeam,xVisDiff,xNirBeam,xNirDiff,& 
                                AlbVisBeam,AlbVisDiff,AlbNirBeam,AlbNirDiff,cldsav,SwToaUp,&
                                CldLow   , CldMed   , CldHgh    )
       END IF
      
   END IF
 END SUBROUTINE RadiationDriver 
 !
 ! Radiation GridHistory Storage 
 !    
 SUBROUTINE RadGridHistStorage (ncols,kmax,latco,CldCovTot,LwCoolRate,ySwHeatRate,LwToaUp,&
                                LwSfcDown,zenith,xSwToaDown,xVisBeam,xVisDiff,xNirBeam,xNirDiff,& 
                                AlbVisBeam,AlbVisDiff,AlbNirBeam,AlbNirDiff,cldsav,SwToaUp,&
                                CldLow   , CldMed   , CldHgh    )
   IMPLICIT NONE
  INTEGER, INTENT(IN   ) :: ncols
  INTEGER, INTENT(IN   ) :: kmax
  INTEGER, INTENT(IN   ) :: latco
  REAL(KIND=r8),    INTENT(in   ) :: CldCovTot(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: LwCoolRate   (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: ySwHeatRate   (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: LwToaUp(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: LwSfcDown(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: zenith(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xSwToaDown(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xVisBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xVisDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xNirBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xNirDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbVisBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbVisDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbNirBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbNirDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: cldsav(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: SwToaUp(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: CldLow (ncols) 
  REAL(KIND=r8),    INTENT(in   ) :: CldMed (ncols) 
  REAL(KIND=r8),    INTENT(in   ) :: CldHgh (ncols) 

  REAL(KIND=r8) :: bfr1(ncols)
  INTEGER       :: i
    
  IF (dogrh(nGHis_vdtclc,latco)) CALL StoreGridHistory(CldCovTot,nGHis_vdtclc,latco)
  IF (dogrh(nGHis_lwheat,latco)) CALL StoreGridHistory(   LwCoolRate,nGHis_lwheat,latco)
  IF (dogrh(nGHis_swheat,latco)) CALL StoreGridHistory(   ySwHeatRate,nGHis_swheat,latco)
  IF (dogrh(nGHis_lwutop,latco)) CALL StoreGridHistory(LwToaUp,nGHis_lwutop,latco)
  IF (dogrh(nGHis_lwdbot,latco)) CALL StoreGridHistory(LwSfcDown,nGHis_lwdbot,latco)
  IF (dogrh(nGHis_coszen,latco)) CALL StoreGridHistory(zenith,nGHis_coszen,latco)
  IF (dogrh(nGHis_swdtop,latco)) CALL StoreGridHistory(xSwToaDown,nGHis_swdtop,latco)
  IF (dogrh(nGHis_swdbvb,latco)) CALL StoreGridHistory( xVisBeam,nGHis_swdbvb,latco)
  IF (dogrh(nGHis_swdbvd,latco)) CALL StoreGridHistory( xVisDiff,nGHis_swdbvd,latco)
  IF (dogrh(nGHis_swdbnb,latco)) CALL StoreGridHistory( xNirBeam,nGHis_swdbnb,latco)
  IF (dogrh(nGHis_swdbnd,latco)) CALL StoreGridHistory( xNirDiff,nGHis_swdbnd,latco)
  IF (dogrh(nGHis_vibalb,latco)) CALL StoreGridHistory( AlbVisBeam,nGHis_vibalb,latco)
  IF (dogrh(nGHis_vidalb,latco)) CALL StoreGridHistory( AlbVisDiff,nGHis_vidalb,latco)
  IF (dogrh(nGHis_nibalb,latco)) CALL StoreGridHistory( AlbNirBeam,nGHis_nibalb,latco)
  IF (dogrh(nGHis_nidalb,latco)) CALL StoreGridHistory( AlbNirDiff,nGHis_nidalb,latco)
  IF (dogrh(nGHis_cloudc,latco)) CALL StoreGridHistory(cldsav,nGHis_cloudc,latco)
  IF (dogrh(nGHis_swutop,latco)) CALL StoreGridHistory(SwToaUp,nGHis_swutop,latco)
  IF (dogrh(nGHis_cldlow,latco)) CALL StoreGridHistory(CldLow,nGHis_cldlow,latco)
  IF (dogrh(nGHis_cldmed,latco)) CALL StoreGridHistory(CldMed,nGHis_cldmed,latco)
  IF (dogrh(nGHis_cldHig,latco)) CALL StoreGridHistory(CldHgh,nGHis_cldHig,latco)

     ! SW Sfc Down
  IF(dogrh(nGHis_swdgrd,latco))THEN
        DO i=1,ncols
           bfr1(i)=xVisDiff(i)+xVisBeam(i)+xNirDiff(i)+xNirBeam(i)
        END DO
        CALL StoreGridHistory(bfr1,nGHis_swdgrd,latco)
  END IF
     ! SW Sfc Up
  IF(dogrh(nGHis_swugrd,latco))THEN
        DO i=1,ncols
           bfr1(i)=AlbVisDiff(i)*xVisDiff(i)+AlbVisBeam(i)*xVisBeam(i)+ &
                   AlbNirDiff(i)*xNirDiff(i)+AlbNirBeam(i)*xNirBeam(i)
        END DO
        CALL StoreGridHistory(bfr1,nGHis_swugrd,latco)
  END IF
   
 END SUBROUTINE RadGridHistStorage
 !
 ! Radiation Diagnostics Storage 
 !
 SUBROUTINE RadDiagStor (ncols,kmax,latco,&
      CldCovTot,CldCovInv,CldCovSat,CldCovCon,CldCovSha,&
      CldLiqWatPath,taud, rei       ,rel   ,emisd, &
      LwCoolRate,ySwHeatRate,ySwHeatRateC,LwCoolRateC,&
      cldsav,o3mix,co2m,&
      xSwToaDown,SwToaUp,SwToaUpC,LwSfcDown,LwToaUp,LwSfcDownC,LwToaUpC,LwSfcNet,LwSfcNetC,xVisDiff,& 
      xVisBeam,xNirDiff,xNirBeam,xVisDiffC,xVisBeamC,xNirDiffC,xNirBeamC,&
      AlbVisDiff,AlbVisBeam,AlbNirDiff,AlbNirBeam, &
      DeltaP, gps,dump)
!, SwSfcDown, SwSfcDownC )

  IMPLICIT NONE
  INTEGER, INTENT(IN   ) :: ncols
  INTEGER, INTENT(IN   ) :: kmax
  INTEGER, INTENT(IN   ) :: latco
  REAL(KIND=r8),    INTENT(in   ) :: CldCovTot(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: CldCovInv(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: CldCovSat(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: CldCovCon(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: CldCovSha(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: CldLiqWatPath  (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: taud  (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: rei   (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: rel   (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: emisd (ncols,kmax)
  ! Longwave
  REAL(KIND=r8),    INTENT(in   ) :: LwCoolRate (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: LwCoolRateC(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: LwSfcDown (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: LwToaUp   (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: LwSfcDownC(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: LwSfcNet  (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: LwSfcNetC (ncols)

  ! Shortwave
  REAL(KIND=r8),    INTENT(in   ) :: ySwHeatRate (ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: ySwHeatRateC(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: xSwToaDown  (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: SwToaUp (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: SwToaUpC(ncols)

  REAL(KIND=r8),    INTENT(in   ) :: xVisDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xVisBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xNirDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xNirBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xVisDiffC(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xVisBeamC(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xNirDiffC(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: xNirBeamC(ncols)

  REAL(KIND=r8),    INTENT(in   ) :: cldsav(ncols)
  REAL(KIND=r8),    INTENT(in   ) :: o3mix (ncols,kMax)   
  REAL(KIND=r8),    INTENT(in   ) :: co2m(nCols,kMax)
  REAL(KIND=r8),    INTENT(in   ) :: LwToaUpC(ncols)

  REAL(KIND=r8),    INTENT(IN   ) :: DeltaP(ncols,kmax)
  REAL(KIND=r8),    INTENT(in   ) :: gps   (ncols)
  REAL(KIND=r8),    INTENT(INOUT ) :: dump(ncols,kmax)
!  REAL(KIND=r8),    INTENT(in   ) :: SwSfcUp    (ncols) ! Upward Shortwave Flux at Surface
!  REAL(KIND=r8),    INTENT(in   ) :: SwSfcUpC   (ncols) ! Upward Shortwave Flux at Surface (clear)
!  REAL(KIND=r8),    INTENT(in   ) :: SwSfcNet   (ncols) ! Shortwave Flux Absorbed at Surface
!  REAL(KIND=r8),    INTENT(in   ) :: SwSfcNetC  (ncols) ! Shortwave Flux Absorbed at Surface (clear)
!  REAL(KIND=r8),    INTENT(in   ) :: SwSfcDown  (ncols) ! Downward Shortwave Flux at Surface
!  REAL(KIND=r8),    INTENT(in   ) :: SwSfcDownC (ncols) ! Downward Shortwave Flux at Surface (clear)

  ! Surface Properties
  REAL(KIND=r8),    INTENT(in   ) :: AlbVisDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbVisBeam (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbNirDiff (ncols)
  REAL(KIND=r8),    INTENT(in   ) :: AlbNirBeam (ncols)
  
  REAL(KIND=r8) :: bfr1(ncols), rfac2
  INTEGER       :: i, k
     ! ====  CO2 
     
     IF(dodia(nDiag_co2aer)) CALL updia(co2m,nDiag_co2aer,latco)
     ! ==== OZONE

     IF(dodia(nDiag_ozonmr)) CALL updia(o3mix,nDiag_ozonmr,latco) !hmjb

     ! total integrated ozone (in dobson units)
     !IF(dodia(nDiag_viozoc)) THEN
        bfr1 = 0.0_r8
        DO k=1,kmax
          !rfac2=gm2dob*100.0_r8*delsig(k)/grav
           DO i=1,ncols
              rfac2=gm2dob*100.0_r8*DeltaP(i,k)/grav
              bfr1(i)=  bfr1(i)+rfac2*gps(i)*o3mix(i,k)
           END DO
        END DO
        !DO i=1,ncols
        !    dump(i,28)=bfr1(i)
        !END DO

     IF(dodia(nDiag_viozoc)) THEN
        CALL updia(bfr1,nDiag_viozoc,latco) !hmjb
     ENDIF

     ! ==== CLOUDS

     IF(dodia(nDiag_iceper)) CALL updia(rei*1.0e-6_r8,nDiag_iceper,latco)!microns --> meters

     IF(dodia(nDiag_liqper)) CALL updia(rel*1.0e-6_r8,nDiag_liqper,latco)!microns --> meters

     IF(dodia(nDiag_vdtclc)) CALL updia(CldCovTot,nDiag_vdtclc,latco)

     IF(dodia(nDiag_invcld)) CALL updia(CldCovInv,nDiag_invcld,latco)

     IF(dodia(nDiag_ssatcl)) CALL updia(CldCovSat,nDiag_ssatcl,latco)

     IF(dodia(nDiag_cnvcld)) CALL updia(CldCovCon,nDiag_cnvcld,latco)

     IF(dodia(nDiag_shcvcl)) CALL updia(CldCovSha,nDiag_shcvcl,latco)

     IF(dodia(nDiag_clliwp)) CALL updia(CldLiqWatPath,nDiag_clliwp,latco)

     IF(dodia(nDiag_sclopd)) CALL updia(taud,nDiag_sclopd,latco)

     IF(dodia(nDiag_lwcemi)) CALL updia(emisd,nDiag_lwcemi,latco)

     IF(dodia(nDiag_cloudc)) CALL updia(cldsav,nDiag_cloudc,latco)

     ! ==== LONG WAVE COOLING RATES

     IF(dodia(nDiag_lwheat)) CALL updia(LwCoolRate,nDiag_lwheat,latco)

     IF(dodia(nDiag_lwhtcl)) CALL updia(LwCoolRateC,nDiag_lwhtcl,latco) !hmjb

     ! ==== LONG WAVE ABSORPTION

     IF(dodia(nDiag_lwnetb)) CALL updia(LwSfcNet ,nDiag_lwnetb,latco)

     IF(dodia(nDiag_lwnbtc)) CALL updia(LwSfcNetC,nDiag_lwnbtc,latco)

     ! ==== LONG WAVE FLUXES

     IF(dodia(nDiag_lwdbot)) CALL updia(LwSfcDown,nDiag_lwdbot,latco)
  
     IF(dodia(nDiag_lwdbtc)) CALL updia(LwSfcDownC,nDiag_lwdbtc,latco)

     IF(dodia(nDiag_lwutop)) CALL updia(LwToaUp,nDiag_lwutop,latco)

     IF(dodia(nDiag_lwutpc)) CALL updia(LwToaUpC,nDiag_lwutpc,latco)

     ! ==== SW HEATING RATES

     IF(dodia(nDiag_swheat)) CALL updia(ySwHeatRate,nDiag_swheat,latco)

     IF(dodia(nDiag_swhtcl)) CALL updia(ySwHeatRateC,nDiag_swhtcl,latco) !hmjb

     ! ==== SW FLUXES

     ! SW Sfc Down
     IF(dodia(nDiag_swdbot))THEN
        DO i=1,ncols
           bfr1(i)=xVisDiff(i)+xVisBeam(i)+xNirDiff(i)+xNirBeam(i)
        END DO
        CALL updia(bfr1,nDiag_swdbot,latco)
     END IF

     ! SW Sfc Down (clear)
     IF(dodia(nDiag_swdbtc))THEN
        DO i=1,ncols
           bfr1(i)=xVisDiffC(i)+xVisBeamC(i)+xNirDiffC(i)+xNirBeamC(i)
        END DO
        CALL updia(bfr1,nDiag_swdbtc,latco)
     END IF

     ! SW Sfc Up
     IF(dodia(nDiag_swubot))THEN
        DO i=1,ncols
           bfr1(i)=AlbVisDiff(i)*xVisDiff(i)+AlbVisBeam(i)*xVisBeam(i)+ &
                   AlbNirDiff(i)*xNirDiff(i)+AlbNirBeam(i)*xNirBeam(i)
        END DO
        CALL updia(bfr1,nDiag_swubot,latco)
     END IF

     ! SW Sfc Up (clear)
     IF(dodia(nDiag_swubtc))THEN
        DO i=1,ncols
           bfr1(i)=AlbVisDiff(i)*xVisDiffC(i)+AlbVisBeam(i)*xVisBeamC(i)+ &
                   AlbNirDiff(i)*xNirDiffC(i)+AlbNirBeam(i)*xNirBeamC(i)
        END DO
        CALL updia(bfr1,nDiag_swubtc,latco)
     END IF

     IF(dodia(nDiag_swutop)) CALL updia(SwToaUp,nDiag_swutop,latco)

     IF(dodia(nDiag_swutpc)) CALL updia(SwToaUpC,nDiag_swutpc,latco)

     IF(dodia(nDiag_swdtop)) CALL updia(xSwToaDown,nDiag_swdtop,latco)

     ! ==== SW ABSORPTION

     ! shortwave absorbed by the earth and the atmosphere
     IF(dodia(nDiag_swabea))THEN
        DO i=1,ncols
           bfr1(i)=xSwToaDown(i)-SwToaUp(i)
        END DO
        CALL updia(bfr1,nDiag_swabea,latco)
     END IF

     ! shortwave absorbed by the earth and the atmosphere (clear)
     IF(dodia(nDiag_swaeac))THEN
        DO i=1,ncols
           bfr1(i)=xSwToaDown(i)-SwToaUpC(i)
        END DO
        CALL updia(bfr1,nDiag_swaeac,latco)
     END IF

     IF(dodia(nDiag_swabgr))THEN
        DO i=1,ncols
           bfr1(i)=(1.0_r8-AlbVisDiff(i))*xVisDiff(i) &
                  +(1.0_r8-AlbVisBeam(i))*xVisBeam(i) &
                  +(1.0_r8-AlbNirDiff(i))*xNirDiff(i) &
                  +(1.0_r8-AlbNirBeam(i))*xNirBeam(i)
        END DO
        CALL updia(bfr1,nDiag_swabgr,latco)
     END IF
!     IF(dodia(nDiag_swabgr)) CALL updia(SwSfcNet,nDiag_swabgr,latco)

     IF(dodia(nDiag_swabgc))THEN
        DO i=1,ncols
           bfr1(i)=(1.0_r8-AlbVisDiff(i))*xVisDiffC(i) &
                       +(1.0_r8-AlbVisBeam(i))*xVisBeamC(i) &
                       +(1.0_r8-AlbNirDiff(i))*xNirDiffC(i) &
                       +(1.0_r8-AlbNirBeam(i))*xNirBeamC(i)
        END DO
        CALL updia(bfr1,nDiag_swabgc,latco)
     END IF

 END SUBROUTINE RadDiagStor 

 SUBROUTINE DestroyRadiationDriver()
   IMPLICIT NONE
   
    
 END SUBROUTINE DestroyRadiationDriver




 SUBROUTINE rqvirt(ncols, kmax, prsl,tin, qin, tout, qout, relhum )
    !
    !==========================================================================
    ! rqvirt: Converts the current model virtual temperature to
    !         thermodynamic temperature for radiation calculation.
    !         Extrapolates moisture up into model dry layers using 
    !         exponential decrease up to a minimum value.
    !         This code requires moisture defined in all layers
    !
    !==========================================================================
    !
    !     input  - sfc pres
    !            - virtual temperature
    !            - moisture
    !
    !     output - thermodynamic temperature
    !            - moiture
    !            - relative humidity
    !
    !==========================================================================
    !
    !  -- Grid Info --
    !  ncols......Number of atmospheric columns
    !  kmax.......Number of vertical layers
    !  sig........Sigma coordinate at middle of layer
    !
    !  -- Atmospheric fields --
    !  gps........Surface pressure in mb
    !  tin........Virtual temperature in K
    !  qin........Moisture in g/g
    !
    !  tout.......Thermodynamic temperature in K
    !  qout.......Extrapolated moisture into model dry layers
    !  relhum.....Relative humidity relhum bound to 0-1
    !
    !  -- Physical Constants --     
    !  delq.......constant delq = 0.608e0
    !  qmin.......constant qmin = 1.0e-12
    !  hl.........heat of evaporation of water     (j/kg)
    !  gasr.......gas constant of dry air        (j/kg/k)
    !  rmwmd......fracao molar da agua e do ar seco
    !  rmwmdi.....fracao molar
    !  e0c........
    !  tbase......constant tbase =  273.15e00
    !==========================================================================

    ! Input Variables
    INTEGER,          INTENT(IN   ) :: ncols
    INTEGER,          INTENT(IN   ) :: kmax
!    REAL(KIND=r8),    INTENT(IN   ) :: sig   (kmax)

!    REAL(KIND=r8),    INTENT(IN   ) :: gps   (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: tin   (ncols,kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: qin   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: prsl  (ncols,kMax)

    ! Output Variables
    REAL(KIND=r8),    INTENT(INOUT) :: tout  (ncols,kmax)
    REAL(KIND=r8),    INTENT(INOUT) :: qout  (ncols,kmax)
    REAL(KIND=r8),    INTENT(INOUT) :: relhum(ncols,kmax)

    ! Local Variables
    REAL(KIND=r8) :: qsat(ncols,kmax)
    REAL(KIND=r8) :: prs (ncols,kmax)
    REAL(KIND=r8) :: tsp (nCols,kMax)      ! saturation temp (K)
    REAL(KIND=r8) :: qsp (nCols,kMax)      ! saturation mixing ratio (kg/kg)

    !REAL(KIND=r8) :: rrlrv
    REAL(KIND=r8) :: const
    REAL(KIND=r8) :: fac ! Saturation vapor pressure [mb]
    INTEGER :: k1
    INTEGER :: k
    INTEGER :: i


    !
    ! Get layer pressure 
    !
    !DO k1 = 1, MOD(kmax,4)
    !   prs(:,k1) = gps(:) * sig(k1)
    !END DO
    !DO k = k1, kmax, 4
    !   prs(:,k)   = gps(:) * sig(k)
    !   prs(:,k+1) = gps(:) * sig(k+1)
    !   prs(:,k+2) = gps(:) * sig(k+2)
    !   prs(:,k+3) = gps(:) * sig(k+3)
    !END DO
    DO k = 1, kmax
       DO i = 1, ncols
          prs(i,k) = prsl  (i,k)! Pa
       END DO
    END DO

    !
    ! Set minimum specific humidity
    ! 
    qout=qin
    WHERE(qout<qmin)
       qout=qmin
    END WHERE

    !
    ! Convert model virtual temp to thermodynmic temp
    !
    tout = tin / (1.0_r8+delq*qout)
     
    !call findsp (nCols,kMax, qout, tout, prsl, tsp, qsp)
    !
    ! Find saturated moisture
    !
    !rrlrv = -hl/(rmwmdi*gasr)
    !const = e0c*EXP(-rrlrv/tbase)
    DO k = 1, kmax
       DO i = 1, ncols
          !IF(tout(i,k)>270.0_r8)THEN
          !   fac = const*EXP(rrlrv/tout(i,k)) !mb
          !   qsat(i,k)=rmwmd* fac/MAX((prsl(i,k)/100.0_r8 - (1.0_r8-rmwmd)*fac),1.0e-12_r8)
          !ELSE 
             fac = fpvs2es5(tout(i,k))/100.0_r8 !Pa ->mb
             qsat(i,k)=rmwmd* fac/MAX(((prsl(i,k)/100.0_r8) - (1.0_r8-rmwmd)*fac),qmin)
          !END IF
          !IF(tout(i,k)>270.0_r8)THEN
          !   fac = const*EXP(rrlrv/tout(i,k))
          !    qsat(i,k)=rmwmd* fac/MAX((prsl(i,k)/100.0_r8 - (1.0_r8-rmwmd)*fac),1.0e-12_r8)
          !ELSE 
          !    qsat(i,k)=qsp(i,k)
          !END IF
       END DO
    END DO

    !
    ! Get relative humidity >=0 and <=1
    !
    DO k = 1, kmax
       DO i = 1, ncols
          relhum(i,k) = qout(i,k) / qsat(i,k)
          !IF (relhum(i,k) <= 0.0_r8) THEN
          !   qout(i,k) = qmin
          !   relhum(i,k) = qmin / qsat(i,k)
          !END IF
          IF (relhum(i,k) > 1.0_r8) THEN
             qout(i,k) = qsat(i,k)
             relhum(i,k) = 1.0_r8
          END IF
       END DO
    END DO
  END SUBROUTINE rqvirt

  SUBROUTINE COSZMED(idatec,tod,yrl,colrad,lonrad,cos2,ncols)
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: ncols
    INTEGER, INTENT(IN   )  :: idatec(4)
    REAL(KIND=r8)   , INTENT(IN   )  :: tod
    REAL(KIND=r8)   , INTENT(IN   )  :: yrl
    REAL(KIND=r8)   , INTENT(IN   )  :: colrad
    REAL(KIND=r8)   , INTENT(in   )  :: lonrad (ncols)
    REAL(KIND=r8)   , INTENT(INOUT  )  :: cos2

    REAL(KIND=r8)                    :: sindel
    REAL(KIND=r8)                    :: cosdel
    REAL(KIND=r8)                    :: cosmax
    REAL(KIND=r8)                    :: ctime
    REAL(KIND=r8)                    :: frh
    REAL(KIND=r8)                    :: btime
    REAL(KIND=r8)                    :: atime
    REAL(KIND=r8)                    :: zenith1
    REAL(KIND=r8)                    :: zenith2
    REAL(KIND=r8)                    :: zenith  (ncols)
    REAL(KIND=r8)                    :: f3600 =3.6e3_r8
    INTEGER                 :: ncount
    REAL(KIND=r8)                    :: alon=0.0_r8
    REAL(KIND=r8)              :: sdelt
    REAL(KIND=r8)              :: ratio
    REAL(KIND=r8)              :: etime
    REAL(KIND=r8)              :: xday
    INTEGER                 :: i
    CALL radtim(idatec,sdelt ,ratio ,etime ,tod   ,xday  ,yrl)

    sindel = SIN(sdelt)
    cosdel = COS(sdelt)
    ctime  = alon/15.0e0_r8
    cos2   = 0.0e0_r8
    cosmax = 0.0e0_r8
    btime  = 0.0e0_r8
    atime  = 0.0e0_r8
    zenith1  = sindel*COS(colrad)
    zenith2  = cosdel*SIN(colrad)
    frh=( MOD(tod+0.03125_r8,f3600)-0.03125_r8)/f3600
    ncount =0
    DO i=1,ncols
       btime       = fim24*lonrad(i)+ctime
       atime       = etime+pai12*(12.0_r8-idatec(1)-frh-btime)
       zenith(i)   = zenith1 + zenith2*COS(atime)
       IF(zenith(i).GT.0.0e0_r8) THEN
          ncount   =ncount+1
          cosmax     =cosmax+zenith(i)
       END IF
    END DO
    IF(ncount.EQ.0) ncount=1
    cos2=cosmax/REAL(ncount,r8)
  END SUBROUTINE COSZMED

  ! radtim :calculates the astronomical parameters: solar inclination,
  !         correction factor to local time, factor relating to the distance
  !         between earth and sun, and the julian day.
  SUBROUTINE radtim (id    ,delta ,ratio ,etime ,tod   ,xday  ,yrl)
    !
    !==========================================================================
    !
    !==========================================================================
    !  id(1)....hour(00/12)
    !  id(2)....month
    !  id(3)....day of month
    !  id(4)....year
    !  delta....solar inclination
    !  ratio....factor relating to the distance between the earth and the sun
    !  etime....correction factor to local time
    !  tod......model forecast time of day in seconds
    !  xday.....is julian day - 1 with fraction of day
    !  pai......constant pi=3.1415926
    !  yrl......length of year in days
    !  monl.....length of each month in days
    !==========================================================================
    INTEGER, INTENT(in ) :: id(4)
    REAL(KIND=r8),    INTENT(out) :: delta
    REAL(KIND=r8),    INTENT(out) :: ratio
    REAL(KIND=r8),    INTENT(out) :: etime
    REAL(KIND=r8),    INTENT(in ) :: tod
    REAL(KIND=r8),    INTENT(out) :: xday
    REAL(KIND=r8),    INTENT(in ) :: yrl

    REAL(KIND=r8),    PARAMETER :: day0=-1.0_r8
    REAL(KIND=r8),    PARAMETER :: f3600=3.6e3_r8
    REAL(KIND=r8)          :: psi
    INTEGER                :: yi,mi,di,hi,LenYearbyDay,nday2y
    REAL(KIND=r8)          :: fdayjul
    !
    !     id is now assumed to be the current date and hour
    !
    yi=id(4)
    mi=id(2)
    di=id(3)
    hi=id(1)
    CALL jull(yi,mi,di,hi,tod,fdayjul,nday2y,LenYearbyDay)
    xday=fdayjul
    IF (xday > day0) THEN
       psi=2.0e0_r8*pai*xday/yrl      
       IF(iyear_AD /= 0)THEN      
           CALL shr_orb_decl(xday,delta,ratio)
       ELSE
          delta=0.006918e0_r8-0.399912e0_r8*COS(   psi)+0.070257e0_r8*SIN(psi) &
               -0.006758e0_r8*COS(2.0e0_r8*psi)+0.000907e0_r8*SIN(2.0e0_r8*psi) &
               -0.002697e0_r8*COS(3.0e0_r8*psi)+0.001480e0_r8*SIN(3.0e0_r8*psi)
          ratio=1.000110e0_r8+0.034221e0_r8*COS(   psi)+0.001280e0_r8*SIN(psi) &
               +0.000719e0_r8*COS(2.0e0_r8*psi)+0.000077e0_r8*SIN(2.0e0_r8*psi)
       END IF  
       etime=0.000075e0_r8+0.001868e0_r8*COS(   psi)-0.032077e0_r8*SIN(psi) &
            -0.014615e0_r8*COS(2.0e0_r8*psi)-0.040849e0_r8*SIN(2.0e0_r8*psi)
    ELSE
       WRITE(nfprt,20)id,tod,xday
       WRITE(nferr,20)id,tod,xday
       STOP 2020
    END IF

20  FORMAT(' BAD DATE IN RADTIM.  ID=',4I5,' TOD=',G16.8,' XDAY=', &
         G16.8)
  END SUBROUTINE radtim
  
  !-----------------------------------------------------------------------
  ! Subroutine: SPMRAD
  !
  ! MAIN ROUTINE FOR RADIATION COMPUTATIONS
  ! CALLS SUBROUTINES:
  !
  !       RADTIM:   COMPUTES ASTRONOMICAL PARAMETERS
  !       GETOZ:    INTERPOLATES OZONE AMOUNT FROM CLIMATOLOGICAL VALUES
  !       CLDGEN:   COMPUTES CLOUD AMOUNTS
  !       SWRAD:    DOES SHORTWAVE RADIATION CALCULATIONS
  !       LWRAD:    DOES LONGWAVE RADIATION CALCULATIONS
  !
  !   COMPUTES COSINES OF SOLAR ZENITH ANGLE
  !   COMPUTES TOTAL CLOUD AMOUNT OF SUPERSATURATION CLOUDS
  !   PREPARES INPUT FLIPP ARRAYS FOR SWRAD AND LWRAD
  !
  !   OUTPUT:
  !
  !      the cooling rate due to long wave radiation, heating rate due to
  !      short wave radiation, downward longwave radiation at the bottom,
  !      and the following relating to downward surface fluxes: visible
  !      beam cloudy skies, visible diffuse cloudy skies, near-infrared
  !      beam cloudy skies, and near-infrared diffuse cloudy skies.
  !
  !
  !-----------------------------------------------------------------------

  SUBROUTINE spmrad ( &
       ! Run Flags
       first  , ifday , lcnvl , lthncl, nfin0 , nfin1 , nfcnv0, &
       intcosz, mxrdcc, inalb , icld  ,                         &
       ! Time info
       yrl    , idatec, idate , tod   , delt  , trint , swint , &
       ! Model Geometry
       colrad , lonrad, cosz  , cos2d , &
       ! Model information
!tar begin       
       !ncols , kmax  , nls   , nlcs  , imask ,         & 
       latco, ncols , kmax  , nls   , nlcs  , imask ,         & 
!tar end              
       ! Atmospheric fields
       prsi   ,prsl   ,phii   ,phil    ,&
       gtmp  , gwv   , gtg   , grh   , omg   , tsea  , &
       QCF    ,QCL    , QCR  ,                                  &
       ! CONVECTION: convective clouds
       convts , convcs, convbs, convc , convt , convb ,         &
       ! SURFACE:  albedo
       AlbVisDiff , AlbNirDiff , AlbVisBeam , AlbNirBeam ,      &
       ! SW Radiation fields at last integer hour
       rSwToaDown,                                              &
       rVisDiff , rNirDiff  , rVisBeam , rNirBeam ,             &
       rVisDiffC, rNirDiffC , rVisBeamC, rNirBeamC,             &
       rSwSfcNet, rSwSfcNetC, rSwToaNet  , rSwToaNetC ,         &
       ! SW Radiation fields at next integer hour
       ySwToaDown,                                              &
       yVisDiff   , yNirDiff   , yVisBeam , yNirBeam ,          &
       yVisDiffC  , yNirDiffC  , yVisBeamC, yNirBeamC,          &
       ySwHeatRate, ySwHeatRatec,                               &
       ySwSfcNet  , ySwSfcNetC , ySwToaNet  , ySwToaNetC ,      &
       ! SSIB INIT: Solar radiation with cos2
       ssib_VisDiff, ssib_NirDiff, ssib_VisBeam, ssib_NirBeam,  &
       ! LW Radiation fields at last integer hour
       LwCoolRate , LwSfcDown  , LwSfcNet , LwToaUp  ,          &
       LwCoolRatec, LwSfcDownC , LwSfcNetC, LwToaUpC ,          &
       ! Cloud field
       CldLow    , CldMed   , CldHgh   , CloudCover, CldCovTot , &
       CldCovInv , CldCovSat, CldCovCon, CldCovSha,             &
       ! Microphysics
       CldLiqWatPath  , emisd , taud  ,  rei      ,rel      , &
       EFFCS ,EFFIS  ,       &
       ! Chemistry
       o3mix ,co2m,dump,CLDF,&
!tar begin 
!climate aerosol optical parameters of coarse aerosol mode and topography field
       aod,asy,ssa,z_aer,topog, &
!tar end
!
!tar begin 
!climate aerosol optical parameters of coarse aerosol mode and topography field
       aodF,asyF,ssaF,z_aerF)
!tar end
!
    IMPLICIT NONE
    !==========================================================================
    !
    ! _________
    ! RUN FLAGS
    !
    ! first....control logical variable .true. or .false.
    ! ifday....model forecast day
    ! lcnvl....the lowest layer index where non-convective clouds can
    !          occur (ben says this should be 2 or more)
    !          constant lcnvl = 2
    ! lthncl...Minimum depth in mb of non-zero low level cloud
    !          constant lthncl=80
    ! nfin0....input  file at time level t-dt
    ! nfin1....input  file at time level t
    ! nfcnv0...initial information on convective clouds for int. radiation
    ! intcosz
    ! mxrdcc...use maximum random converage for radiative conv. clouds
    !          constant logical mxrdcc = .true.
    ! inalb....Select two types of surface albedo
    !          =1 : input two  types surface albedo (2 diffused)
    !               direct beam albedos are calculated by the subr.
    !          =2 : input four types surfc albedo (2 diff,2 direct)
    ! icld.....Select three types of cloud emissivity/optical depth
    !          =1 : old cloud emissivity (optical depth) setting
    !               ccu :  0.05 *dp
    !               css :  0.025*dp       for ice cloud t<253.0
    !                      0.05 *dp       for ice cloud t>253.0
    !          =2 : new cloud emissivity (optical depth) setting
    !               ccu : (0.16)*dp
    !               css :  0.0                              t<-82.5
    !                     (2.0e-6*(t-tcrit)**2)*dp    -82.5<t<-10.0
    !                     (6.949e-3*(t-273)+.08)*dp   -10.0<t<  0.0
    !                     (0.08)*dp                   -10.0<t<  0.0
    !          =3    : ccm3 based cloud emissivity
    !
    ! _________
    ! TIME INFO
    !
    ! yrl.........length of year in days
    !     idatec.......date of current data
    !     idatec(1)....hour(00/12)
    !     idatec(2)....month
    !     idatec(3)....day of month
    !     idatec(4)....year
    !   idate(1)=initial hour of day
    !   idate(2)=day of month.
    !   idate(3)=month of year.
    !   idate(4)=year.
    ! tod.........model forecast time of day in seconds
    ! delt........time interval in sec (fixed throuh the integration)
    ! trint.......ir subr. call interval in hours
    ! swint.......sw subr. call interval in hours
    !             swint has to be less than or equal to trint
    !                              and mod(trint,swint)=0
    !
    ! ______________
    ! MODEL GEOMETRY
    !
    ! colrad.....colatitude  colrad=0-3.14 from np to sp in radians
    ! lonrad.....longitude in radians
    ! cosz.......cosine of solar zenith angle
    ! cos2d
    ! sigmid.....sigma coordinate at middle of layer
    ! sigbot.....sigma coordinate at bottom of layer
    !
    ! _________________
    ! MODEL INFORMATION
    !
    ! latco.....latitude
    ! ncols.....Number of grid points on a gaussian latitude circle
    ! kmax......Number of sigma levels
    ! nls..... .Number of layers in the stratosphere.
    ! nlcs......nlcs =   30
    ! imask.....mascara continental
    !tar begin
    ! topog....topography field
    !tar end    
    !
    ! __________________
    ! ATMOSPHERIC FIELDS
    !
    ! gps......surface pressure in mb
    ! gtmp.....temperature in kelvin
    ! gwv......specific humidity in g/g
    ! gtg......ground surface temperature in kelvin
    ! grh......grh   =  relative humidity  (0-1)
    ! omg......omg   =  vertical velocity  (cb/sec)
    ! tsea.....effective surface radiative temperature ( tgeff )
    !
    ! _____________________________
    ! CONVECTION: convective clouds
    !
    ! convts
    ! convcs
    ! convbs
    ! convc....ncols convective cloud cover in 3 hr. avrage
    ! convt....ncols convective cloud top  (sigma layer)
    ! convb....ncols convective cloud base (sigma layer)
    !
    ! ________________
    ! SURFACE:  albedo
    !
    ! AlbVisDiff.....visible diffuse surface albedo
    ! AlbNirDiff.....near-ir diffuse surface albedo
    ! AlbVisBeam.....visible beam surface albedo
    ! AlbNirBeam.....near-ir beam surface albedo
    !
    ! ________________________________________
    ! SW Radiation fields at last integer hour
    !
    ! rSwToaDown......Incident SW at top (W/m^2)                
    ! rVisBeam....... Down Sfc SW flux visible beam    (all-sky)
    ! rVisDiff....... Down Sfc SW flux visible diffuse (all-sky)
    ! rNirBeam....... Down Sfc SW flux Near-IR beam    (all-sky)
    ! rNirDiff....... Down Sfc SW flux Near-IR diffuse (all-sky)
    ! rVisBeamC...... Down Sfc SW flux visible beam    (clear)  
    ! rVisDiffC...... Down Sfc SW flux visible diffuse (clear)  
    ! rNirBeamC...... Down Sfc SW flux Near-IR beam    (clear)  
    ! rNirDiffC...... Down Sfc SW flux Near-IR diffuse (clear)  
    !
    ! ________________________________________
    ! SW Radiation fields at next integer hour
    !
    ! ySwToaDown.....Incident SW at top 
    ! yVisBeam.......Down Sfc SW flux visible beam    (all-sky)
    ! yVisDiff.......Down Sfc SW flux visible diffuse (all-sky)
    ! yNirBeam.......Down Sfc SW flux Near-IR beam    (all-sky)
    ! yNirDiff.......Down Sfc SW flux Near-IR diffuse (all-sky)
    ! yVisBeamC......Down Sfc SW flux visible beam    (clear)  
    ! yVisDiffC......Down Sfc SW flux visible diffuse (clear)  
    ! yNirBeamC......Down Sfc SW flux Near-IR beam    (clear)  
    ! yNirDiffC......Down Sfc SW flux Near-IR diffuse (clear)  
    ! ySwHeatRate....Heating rate due to shortwave         (K/s)
    ! ySwHeatRateC...Heating rate due to shortwave (clear) (K/s)
    !
    ! ____________________________________________
    ! SSIB INIT: solar radiation with cos2d
    !
    ! ssib_VisBeam.......Down Sfc SW flux visible beam    (all-sky)
    ! ssib_VisDiff.......Down Sfc SW flux visible diffuse (all-sky)
    ! ssib_NirBeam.......Down Sfc SW flux Near-IR beam    (all-sky)
    ! ssib_NirDiff.......Down Sfc SW flux Near-IR diffuse (all-sky)
    !
    ! ________________________________________
    ! LW Radiation fields at next integer 3-hour
    !
    ! LwCoolRate.....Cooling rate due to longwave  (all-sky) (K/s)
    ! LwSfcDown......Down Sfc LW flux              (all-sky) (W/m2)
    ! LwSfcNet.......Net Sfc LW                    (all-sky) (W/m2)
    ! LwToaUp........Longwave upward at top        (all-sky) (W/m2) 
    ! LwCoolRateC....Cooling rate due to longwave  (clear) (K/s)
    ! LwSfcDownC.....Down Sfc LW flux              (clear) (W/m2)
    ! LwSfcNetC......Net Sfc LW                    (clear) (W/m2)
    ! LwToaUpC.......Longwave upward at top        (clear) (W/m2)
    !
    ! ___________
    ! Cloud field
    !
    ! CloudCover.....Total cloud cover
    ! CldCovTot......Total cloud cover (at each layer)
    ! CldCovInv......Inversion clouds                 
    ! CldCovSat......Saturation clouds                
    ! CldCovCon......Convection clouds                
    ! CldCovSha......Shallow convective clouds        
    !
    ! ____________
    ! Microphysics (if cld=3 or 4)
    !
    ! CldLiqWatPath...Cloud liquid water path (parametrized)
    ! emisd...........Emissivity 
    ! taud............Shortwave cloud optical depth
    ! rei.............Ice particle Effective Radius (microns)
    ! rel.............Liquid particle Effective Radius (microns)
    !
    ! _________
    ! Chemistry
    !
    ! co2m......co2m is wgne standard value in mol/mol "co2m = /345.0e-6/
    ! o3mix
    !
    !tar begin
    !climate aerosol optical parameters of coarse aerosol mode
    !    aod....aerosol optical depth
    !    asy....asymmetry factor
    !    ssa....single scattering albedo
    !    z_aer..aod vertical profile
    !tar end
    !
    !tar begin
    !climate aerosol optical parameters of fine aerosol mode
    !    aodF....aerosol optical depth
    !    asyF....asymmetry factor
    !    ssaF....single scattering albedo
    !    z_aerF..aod vertical profile
    !tar end    
    !
    !     visible = 0.0 to 0.7 mu  near-ir = 0.7 to 4.0 mu
    !
    !     o3.......layer ozone mixing ratio in g/g
    !              interpolated from climatological value
    !     cld......supersaturation cloud fraction
    !     clu......convective cloud fraction
    !               above two cloud fraction should be given either as
    !               climatological value or using diagnostic relation
    !               between cloud fraction and relative humidity etc.
    !               cloud amount in ir/sw subr. is max(cld,clu)
    !==========================================================================

    ! Run Flags
    LOGICAL      ,    INTENT(in   ) :: first
    INTEGER      ,    INTENT(in   ) :: ifday
    INTEGER      ,    INTENT(in   ) :: lcnvl
    INTEGER      ,    INTENT(in   ) :: lthncl
    INTEGER      ,    INTENT(in   ) :: nfin0
    INTEGER      ,    INTENT(in   ) :: nfin1
    INTEGER      ,    INTENT(in   ) :: nfcnv0
    LOGICAL      ,    INTENT(IN   ) :: intcosz
    LOGICAL      ,    INTENT(in   ) :: mxrdcc
    INTEGER      ,    INTENT(in   ) :: inalb
    INTEGER      ,    INTENT(in   ) :: icld

    ! Time info
    REAL(KIND=r8),    INTENT(in   ) :: yrl
    INTEGER      ,    INTENT(in   ) :: idatec(4)
    INTEGER      ,    INTENT(in   ) :: idate(4)
    REAL(KIND=r8),    INTENT(in   ) :: tod
    REAL(KIND=r8),    INTENT(in   ) :: delt
    REAL(KIND=r8),    INTENT(in   ) :: trint
    REAL(KIND=r8),    INTENT(in   ) :: swint

    ! Model Geometry
    REAL(KIND=r8),    INTENT(in   ) :: colrad(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: lonrad(ncols)
!tar begin
    REAL(KIND=r8) :: latrad(ncols)    ! latitude in degrees (local variables)
!tar end
    REAL(KIND=r8),    INTENT(inout) :: cosz  (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: cos2d (ncols)

    ! Model information
    !INTEGER         , INTENT(in   ) :: latco
!tar begin
    INTEGER         , INTENT(in   ) :: latco
!tar end
    INTEGER         , INTENT(IN   ) :: ncols
    INTEGER         , INTENT(IN   ) :: kmax
    INTEGER         , INTENT(IN   ) :: nls
    INTEGER         , INTENT(IN   ) :: nlcs
    INTEGER(KIND=i8), INTENT(IN   ) :: imask (ncols)
!tar begin
!topography field
    REAL(KIND=r8),    INTENT(IN   ) :: topog (ncols)
!tar end    

    ! Atmospheric fields
    REAL(KIND=r8),    INTENT(in   ) :: prsi  (ncols,kMax+1)
    REAL(KIND=r8),    INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: phii  (nCols,kMax+1)
    REAL(KIND=r8),    INTENT(in   ) :: phil  (nCols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gtmp  (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gwv   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: gtg   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: grh   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: omg   (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: QCF   (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: QCL   (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: QCR   (ncols,kMax)
    ! CONVECTION: convective clouds
    REAL(KIND=r8),    INTENT(in   ) :: convts(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convcs(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convbs(ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convc (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convt (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: convb (ncols)

    ! SURFACE:  albedo
    REAL(KIND=r8),    INTENT(in   ) :: AlbVisDiff (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: AlbNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: AlbVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: AlbNirBeam (ncols)

    ! SW Radiation fields at last integer hour
    REAL(KIND=r8),    INTENT(inout) :: rSwToaDown(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rVisBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rNirBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rSwSfcNet   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rSwSfcNetC  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rSwToaNet(ncols)
    REAL(KIND=r8),    INTENT(inout) :: rSwToaNetC(ncols)

    ! SW Radiation fields at next integer hour
    REAL(KIND=r8),    INTENT(inout) :: ySwToaDown(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirDiff (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirBeam (ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirDiffC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yVisBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: yNirBeamC(ncols)
    REAL(KIND=r8),    INTENT(inout) :: ySwHeatRate   (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: ySwHeatRateC  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: ySwSfcNet   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ySwSfcNetC  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ySwToaNet(ncols)
    REAL(KIND=r8),    INTENT(inout) :: ySwToaNetC(ncols)

    ! SSIB INIT: Solar radiation with cos2
    REAL(KIND=r8),    INTENT(out) :: ssib_VisDiff (ncols)
    REAL(KIND=r8),    INTENT(out) :: ssib_NirDiff (ncols)
    REAL(KIND=r8),    INTENT(out) :: ssib_VisBeam (ncols)
    REAL(KIND=r8),    INTENT(out) :: ssib_NirBeam (ncols)

    ! LW Radiation fields at last integer hour
    REAL(KIND=r8),    INTENT(inout) :: LwCoolRate (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcDown  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcNet   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwToaUp    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwCoolRateC(ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcDownC (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwSfcNetC  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: LwToaUpC   (ncols)

    ! Cloud field
    REAL(KIND=r8),    INTENT(inout) :: CldLow (ncols) 
    REAL(KIND=r8),    INTENT(inout) :: CldMed (ncols) 
    REAL(KIND=r8),    INTENT(inout) :: CldHgh (ncols) 

    REAL(KIND=r8),    INTENT(inout) :: CloudCover(ncols)
    REAL(KIND=r8),    INTENT(inout) :: CldCovTot (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovInv (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovSat (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovCon (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: CldCovSha (ncols,kmax)

    ! Microphysics
    REAL(KIND=r8),    INTENT(inout) :: CldLiqWatPath  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: emisd (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: taud  (ncols,kmax)
    REAL(KIND=r8),    INTENT(inout) :: rei   (ncols,kmax)! Ice particle Effective Radius (microns)
    REAL(KIND=r8),    INTENT(inout) :: rel   (ncols,kmax)! Liquid particle Effective Radius (microns)
    REAL(KIND=r8),    INTENT(in   ) :: EFFCS (ncols,kmax)
    REAL(KIND=r8),    INTENT(in   ) :: EFFIS (ncols,kmax)
    ! Chemistry
    REAL(KIND=r8),    INTENT(INOUT) :: o3mix(ncols,kMax)   
    REAL(KIND=r8),    INTENT(in   ) :: co2m(ncols,kMax)   !mol/mol
    REAL(KIND=r8),    INTENT(INOUT) :: dump(ncols,kMax)
    REAL(KIND=r8),    INTENT(INOUT) :: CLDF(ncols,kMax)
!tar begin
!climate aerosol optical parameters of coarse mode for RRTMG   
    REAL(KIND=r8),    INTENT(INOUT) :: aod(ncols,14) 
    REAL(KIND=r8),    INTENT(INOUT) :: asy(ncols,14)
    REAL(KIND=r8),    INTENT(INOUT) :: ssa(ncols,14)
    REAL(KIND=r8),    INTENT(INOUT) :: z_aer(ncols,40)     
!tar end
!
!tar begin
!climate aerosol optical parameters of fine mode for RRTMG    
    REAL(KIND=r8),    INTENT(INOUT) :: aodF(ncols,14) 
    REAL(KIND=r8),    INTENT(INOUT) :: asyF(ncols,14)
    REAL(KIND=r8),    INTENT(INOUT) :: ssaF(ncols,14)
    REAL(KIND=r8),    INTENT(INOUT) :: z_aerF(ncols,40)     
!tar end
!
! local aerosol parameters for CRDTF code

!tar begin
!climate aerosol optical parameters of coarse mode
    REAL(KIND=r8) :: aod8(ncols,8)
    REAL(KIND=r8) :: asy8(ncols,8)
    REAL(KIND=r8) :: ssa8(ncols,8)
!tar end
!
!tar begin
!climate aerosol optical parameters of fine mode
    REAL(KIND=r8) :: aodF8(ncols,8)
    REAL(KIND=r8) :: asyF8(ncols,8)
    REAL(KIND=r8) :: ssaF8(ncols,8)
!tar end
          
    !==========================================================================
    ! LOCAL VARIABLES
    !==========================================================================
    REAL(KIND=r8) :: c_cld_tau    (1:nbndsw,1:ncols,1:kmax) ! cloud extinction optical depth
    REAL(KIND=r8) :: c_cld_tau_w  (1:nbndsw,1:ncols,1:kmax) ! cloud single scattering albedo * tau
    REAL(KIND=r8) :: c_cld_tau_w_g(1:nbndsw,1:ncols,1:kmax) ! cloud assymetry parameter * w * tau
    REAL(KIND=r8) :: c_cld_tau_w_f(1:nbndsw,1:ncols,1:kmax) ! cloud forward scattered fraction * w * tau
    REAL(KIND=r8) :: c_cld_lw_abs (1:nbndlw,1:ncols,1:kmax) ! cloud absorption optics depth (LW)
    REAL(KIND=r8) :: cldfprime    (1:ncols,1:kmax)         ! combined cloud fraction (snow plus regular)
    REAL(KIND=r8) :: psurf   (ncols)

    REAL(KIND=r8) :: cos2  (ncols)
    REAL(KIND=r8) :: cos1  (ncols)

    ! FLIPPED atmospheric fields
    REAL(KIND=r8) :: FlipTe   (ncols,kmax)  ! Temperature (K)
    REAL(KIND=r8) :: FlipQe   (ncols,kmax)  ! Specific Humidity (g/g)
    REAL(KIND=r8) :: FlipPmid (ncols,kmax)  ! Pressure at middle of layer (mb)
    REAL(KIND=r8) :: FlipPbot (ncols,kmax)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8) :: FlipPInterface (ncols,kmax+1)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8) :: FlipDP   (ncols,kmax)  ! Layer thickness (mb) 
    REAL(KIND=r8) :: FlipO3   (ncols,kmax)
    REAL(KIND=r8) :: Flip_c_cld_tau    (1:nbndsw,1:ncols,1:kmax) ! cloud extinction optical depth
    REAL(KIND=r8) :: Flip_c_cld_tau_w  (1:nbndsw,1:ncols,1:kmax) ! cloud single scattering albedo * tau
    REAL(KIND=r8) :: Flip_c_cld_tau_w_g(1:nbndsw,1:ncols,1:kmax) ! cloud assymetry parameter * w * tau
    REAL(KIND=r8) :: Flip_c_cld_tau_w_f(1:nbndsw,1:ncols,1:kmax) ! cloud forward scattered fraction * w * tau
    REAL(KIND=r8) :: Flip_c_cld_lw_abs (1:nbndlw,1:ncols,1:kmax) ! cloud absorption optics depth (LW)

    !
    !     arrays temporarily for flipping sw and lw i/o
    !
    REAL(KIND=r8) :: cld    (ncols,kmax)
    REAL(KIND=r8) :: clu    (ncols,kmax)
    REAL(KIND=r8) :: asclr   (ncols,kmax)
    REAL(KIND=r8) :: asl     (ncols,kmax)
    INTEGER :: idatex(4)
    REAL(KIND=r8) :: clwp    (ncols,kmax)   ! Cloud Liquid Water Path
    REAL(KIND=r8) :: fice    (ncols,kmax)   ! fractional amount of cloud that is ice
    REAL(KIND=r8) :: lmixr   (ncols,kmax)   ! ice/water mixing ratio
    REAL(KIND=r8) :: cicewp (ncols,kmax)
    REAL(KIND=r8) :: cliqwp (ncols,kmax)
    REAL(KIND=r8) :: FlipTauCloudTot (ncols,kmax)
    REAL(KIND=r8) :: FlipTauCloudTotSW (ncols,kmax)
    REAL(KIND=r8) :: FlipCO2m    (ncols,kmax)   ! CO2 concentration mol/mol
    REAL(KIND=r8) :: FlipClwp    (ncols,kmax)   ! Cloud Liquid Water Path
    REAL(KIND=r8) :: FlipRei     (ncols,kmax)   ! Ice particle Effective Radius (microns)
    REAL(KIND=r8) :: FlipRel     (ncols,kmax)   ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8) :: FlipFice    (ncols,kmax)   ! fractional amount of cloud that is ice
    REAL(KIND=r8) :: FlipLMixR   (ncols,kmax)   ! ice/water mixing ratio
    REAL(KIND=r8) :: FlipTaud   (ncols,kmax) 
    REAL(KIND=r8) :: Flipcicewp (ncols,kmax) 
    REAL(KIND=r8) :: Flipcliqwp (ncols,kmax) 
    REAL(KIND=r8) :: Flip_cldfprime(ncols,kmax) 
    CHARACTER(len=4), PARAMETER :: covlp='RAND'
    LOGICAL         , PARAMETER :: noz  = .FALSE.

    REAL(KIND=r8)    :: date
    REAL(KIND=r8)    :: s0
    REAL(KIND=r8)    :: datex
    REAL(KIND=r8)    :: s0x
    REAL(KIND=r8)    :: todx
    REAL(KIND=r8)    :: ratio         ! factor relating to the distance between the earth and the sun
    CHARACTER(LEN=*), PARAMETER     :: h='**(RadiationDriver_SPMRAD)**'


    INTEGER :: ifdayx
    INTEGER :: jhr
    INTEGER :: jmon
    INTEGER :: jday
    INTEGER :: jyr
    INTEGER :: i
    INTEGER :: k
    INTEGER :: kk
    INTEGER :: kflip
    !REAL(KIND=r8)    :: pptop

    cos2     =0.0_r8;cos1     =0.0_r8;FlipTe    =0.0_r8;FlipQe   =0.0_r8
    FlipPmid =0.0_r8;FlipPbot =0.0_r8;FlipDP    =0.0_r8;FlipO3   =0.0_r8
    cld      =0.0_r8;clu      =0.0_r8;asclr     =0.0_r8;asl      =0.0_r8
    idatex   =0     ;clwp     =0.0_r8;rei       =0.0_r8;rel      =0.0_r8
    fice     =0.0_r8;lmixr    =0.0_r8; FlipClwp =0.0_r8;FlipRei  =0.0_r8
    FlipRel  =0.0_r8; FlipFice=0.0_r8; FlipLMixR=0.0_r8;rei      =0.0_r8
    rel      =0.0_r8;cicewp=0.0_r8;cliqwp=0.0_r8;Flipcicewp  =0.0_r8
    Flipcliqwp  =0.0_r8;cldfprime=0.0_r8;Flip_cldfprime=0.0_r8;FlipPInterface=0.0_r8
    FlipTauCloudTot=0.0_r8
    FlipTauCloudTotSW=0.0_r8
    CldLow=0.0_r8
    CldMed=0.0_r8
    CldHgh=0.0_r8
    DO i = 1,ncols
       psurf(i) =prsi  (i,1)/100.0_r8
    END DO 

    DO k=1,kmax
       DO i =1, ncols
          IF(prsl(i,k)/100.0_r8 > 850.0_r8)THEN
             CldLow (i) = MAX(CldLow(i), CldCovTot (i,k))
          ELSE IF(prsl(i,k)/100.0_r8 <= 850.0_r8 .and. prsl(i,k)/100.0_r8 >= 400.0_r8 )THEN
             CldMed (i) = MAX(CldMed(i), CldCovTot (i,k))
          ELSE
             CldHgh (i) = MAX(CldHgh(i), CldCovTot (i,k))
          END IF
       END DO
    ENDDO

!PK    IF(ABS( MOD((tod-delt)/3.6e3_r8+0.03125e0_r8,swint)).GT.0.0625e0_r8 &
!PK         .AND..NOT.first)RETURN
    IF(ABS( MOD((tod-delt)+0.03125e0_r8,swint)).GT.0.0625e0_r8 &
         .AND..NOT.first)RETURN


    !-----------------------------------------------------------------------
    ! Computation of astronomical parameters 
    !-----------------------------------------------------------------------

    ! At current full hour
    call astropar(ncols,colrad,lonrad,idatec,tod,yrl,date,s0,cosz,ratio)

    ! Advance 1hour
    ifdayx=ifday
!PK todx=tod+swint*3.6e3_r8
    todx=tod+swint!*3.6e3_r8

    IF(todx.GE.86400.0_r8)THEN
       todx= MOD(todx,86400.0_r8)
       ifdayx=ifdayx+1
    END IF

    CALL tmstmp2(idate ,ifdayx,todx  ,jhr   ,jday  ,jmon  ,jyr   )

    idatex(1)=jhr
    idatex(2)=jmon
    idatex(3)=jday
    idatex(4)=jyr

    ! At next full hour
    call astropar(ncols,colrad,lonrad,idatex,todx,yrl,datex,s0x,cos1,ratio)


    IF(intcosz)THEN
       cos2=cos2d
    ELSE
       cos2=cosz
    END IF

    !-----------------------------------------------------------------------
    ! Ozone
    !-----------------------------------------------------------------------

    !IF (ifozone.eq.0) THEN
    !   CALL getoz (ncols,kmax,prsl,prsi,colrad,date  ,FlipO3 )
    !   DO k=1,kmax
    !      kflip=kmax+1-k
    !      DO i = 1,ncols
    !         o3mix(i,k)=FlipO3(i,kflip)
    !      END DO
    !   END DO
    !ELSE
       ! If ifozone!=0, then we have to use the value read from file
       ! ... but we still have to flip it for the radiation routines
       DO k=1,kmax
          kflip=kmax+1-k
          DO i = 1,ncols
             FlipO3(i,kflip)=o3mix(i,k)
          END DO
       END DO
    !ENDIF

    !-----------------------------------------------------------------------
    ! Flip variables
    !-----------------------------------------------------------------------

    DO k=1,kmax+1
       kflip=kmax+2-k
       DO i = 1,ncols
          FlipPInterface(i,k)=prsi  (i,kflip)/100.0_r8!gps(i)*sigbot(kflip)
       END DO
    END DO

    DO k=1,kmax
       kflip=kmax+1-k
       DO i = 1,ncols
          FlipTe(i,k)=gtmp(i,kflip)
          FlipQe(i,k)=gwv (i,kflip)+ QCF (i,kflip)/(1.0_r8+QCF (i,kflip))+ QCL(i,kflip)/(1.0_r8+QCL(i,kflip))&
                                   + QCR (i,kflip)/(1.0_r8+QCR (i,kflip))
          FlipPbot(i,k)=prsi  (i,kflip)/100.0_r8!gps(i)*sigbot(kflip)
          FlipPmid(i,k)=prsl  (i,kflip)/100.0_r8!gps(i)*sigmid(kflip)
          FlipCO2m(i,k)=co2m(i,kflip)
       END DO
    END DO

    !pptop=0.005e0_r8
    DO i = 1,ncols
       FlipDP(i,1)=FlipPbot(i,1)-pptop
    END DO
    DO k=2,kmax
       DO i = 1,ncols
          FlipDP(i,k)=FlipPbot(i,k)-FlipPbot(i,k-1)
       END DO
    END DO

       !CALL MsgOne(h,'run Cloud_Field')

    !-----------------------------------------------------------------------
    ! Cloud Field and Microphysics (some day this will come from convection)
    !-----------------------------------------------------------------------

    CALL Cloud_Field(&
       ! Model info and flags
       ncols     , kmax     , nls      , icld     , covlp    , mxrdcc   , &
       lcnvl     , lthncl   ,colrad    , imask,                           &
       ! Atmospheric Fields
       prsi   ,prsl   ,phii   ,phil    ,&
       grh      , omg      , gtmp     ,            &
       ! Convective Clouds
       convc     , convt    , convb    , convts   , convcs   , convbs   , &
       CLDF      ,&
       ! Flipped Clouds for Radiation
       cld       , clu      ,                                             &
       ! Clouds for Diagnostics
       CloudCover, CldCovTot, CldCovInv, CldCovSat, CldCovCon, CldCovSha, &
       CldLow    ,CldMed    ,CldHgh  )
       !CALL MsgOne(h,'end Cloud_Field')

       !CALL MsgOne(h,'run icld')

    IF(icld.EQ.3.or.icld.EQ.4)THEN
       CALL Cloud_Micro_CCM3(&
            ! Model info
            ncols, kmax ,imask   , &
            ! Atmospheric Fields
            prsi   ,prsl   ,phii   ,phil    ,&
            gtmp , gwv   , tsea  , FlipPbot ,        &
            ! Cloud properties
            clwp , lmixr, fice  , rei   , rel   , taud       )

       CldLiqWatPath=clwp
    ENDIF

    IF(icld.EQ.5)THEN

       CALL Cloud_Micro_CAM5(&
            ! Model info
            ncols, kmax , imask    , &
            ! Atmospheric Fields
            prsi   ,prsl   ,phii   ,phil    ,&
            gtmp , gwv   , tsea  , FlipPbot ,        &
            QCF  ,QCL   ,                                   &
            ! Cloud properties
            clwp , lmixr, fice  , rei   , rel      , taud       )

       CldLiqWatPath=clwp
    ENDIF

    IF(icld.EQ.6)THEN
       CALL RunCloudOpticalProperty2 (       &
     !  iswrad    ,&
       ncols    ,&
       kmax     ,&
       ILCON    ,&
       imask    ,&
       colrad   ,&
       prsi   ,prsl   ,phii   ,phil    ,&
       gtmp     ,&
       gwv      ,&
       FlipPbot ,&
       QCF      ,&
       QCL      ,&
       QCR      ,&
       gtg      ,&
       tsea     ,&
       cld      ,&
       clu      ,&  
       EFFCS    ,& ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)! DROPLET EFFECTIVE RADIUS   (MICRON)
       EFFIS    ,& ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)! CLOUD ICE EFFECTIVE RADIUS (MICRON)
       ! Cloud properties
       CldLow  ,& 
       CldMed  ,& 
       CldHgh  ,& 
       clwp , lmixr, fice  , rei   , rel      , taud       ,cicewp,cliqwp,&
       c_cld_tau     ,c_cld_tau_w   ,c_cld_tau_w_g ,c_cld_tau_w_f ,c_cld_lw_abs,cldfprime  )
       !-----------------------------------------------------------------------
       ! Flip microphysics variables
       !-----------------------------------------------------------------------
       FlipTauCloudTotSW=0.0_r8
       DO kk=1,nbndsw
          DO k=1,kmax
             kflip=kmax+1-k
             DO i = 1,ncols
                Flip_cldfprime  (i,k)       =cldfprime     (   i,kflip)
                Flip_c_cld_tau    (kk,i,k)  =c_cld_tau     (kk,i,kflip)
                Flip_c_cld_tau_w  (kk,i,k)  =c_cld_tau_w   (kk,i,kflip)
                Flip_c_cld_tau_w_g(kk,i,k)  =c_cld_tau_w_g (kk,i,kflip)
                Flip_c_cld_tau_w_f(kk,i,k)  =c_cld_tau_w_f (kk,i,kflip)

                FlipTauCloudTotSW(i,k)=FlipTauCloudTotSW(i,k)+Flip_c_cld_tau(kk,i,k)
             END DO
          END DO
       END DO
       FlipTauCloudTot=0.0_r8
       DO kk=1,nbndlw
          DO k=1,kmax
             kflip=kmax+1-k
             DO i = 1,ncols
                Flip_c_cld_lw_abs (kk,i,k)  =c_cld_lw_abs  (kk,i,kflip)
                FlipTauCloudTot(i,k)  = FlipTauCloudTot(i,k) + Flip_c_cld_lw_abs (kk,i,k)
             END DO
          END DO
       END DO
       CldLiqWatPath=clwp
    ENDIF

    IF(icld.EQ.7)THEN

      CALL Cloud_Micro_WRF(&
       ! Model info
       ncols, kmax , imask   , ILCON,schemes,&
       !iswrad,ncols, kmax , imask   , ILCON,schemes,&
       ! Atmospheric Fields
       prsi   ,prsl   ,phii   ,phil    ,&
       gtmp   , gwv    , tsea  , FlipPbot ,        &
       QCF  ,QCL   ,  cld ,clu ,                                 &
       ! Cloud properties
       clwp , lmixr, fice  , rei   , rel   , taud ,cicewp,cliqwp,&
       c_cld_tau    ,      &
       c_cld_tau_w  ,      &
       c_cld_tau_w_g,      &
       c_cld_tau_w_f,      &
       c_cld_lw_abs ,      &
       cldfprime          )
       !-----------------------------------------------------------------------
       ! Flip microphysics variables
       !-----------------------------------------------------------------------
       FlipTauCloudTot=0.0_r8
       DO kk=1,nbndsw
          DO k=1,kmax
             kflip=kmax+1-k
             DO i = 1,ncols
                Flip_cldfprime  (i,k)       =cldfprime     (   i,kflip)
                Flip_c_cld_tau    (kk,i,k)  =c_cld_tau     (kk,i,kflip)
                Flip_c_cld_tau_w  (kk,i,k)  =c_cld_tau_w   (kk,i,kflip)
                Flip_c_cld_tau_w_g(kk,i,k)  =c_cld_tau_w_g (kk,i,kflip)
                Flip_c_cld_tau_w_f(kk,i,k)  =c_cld_tau_w_f (kk,i,kflip)

                FlipTauCloudTotSW(i,k)=FlipTauCloudTotSW(i,k)+Flip_c_cld_tau(kk,i,k)
             END DO
          END DO
       END DO

       FlipTauCloudTot=0.0_r8
       DO kk=1,nbndlw
          DO k=1,kmax
             kflip=kmax+1-k
             DO i = 1,ncols
                Flip_c_cld_lw_abs (kk,i,k)  =c_cld_lw_abs  (kk,i,kflip)
                FlipTauCloudTot(i,k)  = FlipTauCloudTot(i,k) + Flip_c_cld_lw_abs (kk,i,k)
             END DO
          END DO
       END DO
       CldLiqWatPath=clwp
    ENDIF

       !CALL MsgOne(h,'end icld')

    !-----------------------------------------------------------------------
    ! Flip microphysics variables
    !-----------------------------------------------------------------------

    DO k=1,kmax
       kflip=kmax+1-k
       DO i = 1,ncols
          FlipClwp  (i,k)=Clwp (i,kflip)
          FlipRei   (i,k)=Rei  (i,kflip)
          FlipRel   (i,k)=Rel  (i,kflip)
          Flipcicewp(i,k)=cicewp (i,kflip)
          Flipcliqwp(i,k)=cliqwp (i,kflip)
          FlipFice  (i,k)=Fice (i,kflip)
          FlipTaud  (i,k)=Taud (i,kflip)
          FlipLMixR (i,k)=LMixR(i,kflip)
       END DO
    END DO
    !WRITE(*,*) 'FlipTaud',MAXVAL(FlipTaud),MINVAL(FlipTaud)
    !WRITE(*,*) 'FlipRei',MAXVAL(FlipRei),MINVAL(FlipRei)
    !WRITE(*,*) 'FlipRel',MAXVAL(FlipRel),MINVAL(FlipRel)
    !WRITE(*,*) 'FlipClwp',MAXVAL(FlipClwp),MINVAL(FlipClwp)
!tar begin 
!latitudes in degrees:
    DO i=1,ncols
      latrad(i)=90.0_r8-(180.0_r8*colrad(i)/3.14159265_r8)
    END DO
!tar end     
!!
!tar begin
! Interpolation of aerosol optical parameters from RRTMG(14) to CRDTF (8)
  aod8=0.0_r8; asy8=0.0_r8; ssa8=0.0_r8
  aodF8=0.0_r8; asyF8=0.0_r8; ssaF8=0.0_r8
!
   CALL RRTM_CRDTF(ncols,aod,asy,ssa,aod8,asy8,ssa8)
   CALL RRTM_CRDTF(ncols,aodF,asyF,ssaF,aodF8,asyF8,ssaF8)
!
!tar end 
!tar print in spmrad before call of radiation codes
!!       IF(myid.EQ.24) then
!!     IF(idatec(1).EQ.6.AND.idatec(3).EQ.1) then
!    
!!     OPEN(unit=75,file='/scratchin/grupos/mcga/home/t.tarassova/OUTPUT_T/Tar.txt', &
!!     ACCESS='APPEND', STATUS='OLD')
!
!!     WRITE(75,*) 'in spmrad before Call of Rad Codes'
!!     WRITE(75,*) 'myid=', myid 
!!     WRITE(75,*) 'hour=', idatec(1), 'day=', idatec(3)    
!!     WRITE(75,*) 'ncols=', ncols 
!!     WRITE(75,*) 'ifaeros=', ifaeros 
!!     WRITE(75,*) 'kmax=', kmax, 'nls=', nls                       
!
!!     WRITE(75,*) 'idatec(4)=h,m,d,y', idatec(1),idatec(2),idatec(3),idatec(4)
!!     WRITE(75,*) 'colrad(1:ncols)=',(colrad(i),i=1,ncols)
!!     WRITE(75,*) 'latrad(1:ncols)=',(latrad(i),i=1,ncols)     
!!     WRITE(75,*) 'lonrad(1:ncols)=',(lonrad(i),i=1,ncols)
!!     WRITE(75,*) 'aod(1:ncols,3(0.512))=', (aod(i,3),i=1,ncols)     
!!     WRITE(75,*) 'asy(1:ncols,3(0.512))=', (asy(i,3),i=1,ncols)    
!!     WRITE(75,*) 'ssa(1:ncols,3(0.512))=', (ssa(i,3),i=1,ncols)

!     WRITE(75,*) 'aod(1:ncols,1(0.252))=', (aod(i,1),i=1,ncols)     
!     WRITE(75,*) 'asy(1:ncols,1(0.252))=', (asy(i,1),i=1,ncols)    
!     WRITE(75,*) 'ssa(1:ncols,1(0.252))=', (ssa(i,1),i=1,ncols) 
     
!     WRITE(75,*) 'aod(1:ncols,8(6.135))=', (aod(i,8),i=1,ncols)     
!     WRITE(75,*) 'asy(1:ncols,8(6.135))=', (asy(i,8),i=1,ncols)    
!     WRITE(75,*) 'ssa(1:ncols,8(6.135))=', (ssa(i,8),i=1,ncols)          
     
!!     WRITE(75,*) 'z_aer(1,1:40)=', (z_aer(1,i),i=1,40)
!!     WRITE(75,*) 'Pressure(1,1:kmax)=', FlipPbot(1,1:kmax)
!!     WRITE(75,*) 'Temp,K (1,1:kmax)=', FlipTe(1,1:kmax)
!!     WRITE(75,*) 'topog(1)=', topog(1) 
     
!!     WRITE(75,*) 'z_aer(70,1:40)=', (z_aer(70,i),i=1,40)
!!     WRITE(75,*) 'Pressure(70,1:kmax)=', FlipPbot(70,1:kmax)
!!     WRITE(75,*) 'Temp,K (70,1:kmax)=', FlipTe(70,1:kmax)
!!     WRITE(75,*) 'topog(70)=', topog(70)     
         
!!     WRITE(75,*) 'z_aer(ncols,1:40)=', (z_aer(ncols,i),i=1,40)
!!     WRITE(75,*) 'Pressure(ncols,1:kmax)=', FlipPbot(ncols,1:kmax)
!!     WRITE(75,*) 'Temp,K (ncols,1:kmax)=', FlipTe(ncols,1:kmax)
!!     WRITE(75,*) 'topog(ncols)=', topog(ncols)      
     
!!    WRITE(75,*) 'topog(ncols)=', topog(1:ncols) 
                     
! 
!!    CLOSE(75) 
!!       ENDIF 
!!            ENDIF      
!tar print
       !CALL MsgOne(h,'run shortwave')
    
    !-----------------------------------------------------------------------
    !     shortwave radiation - X-call = special call for SSIB init
    !-----------------------------------------------------------------------

    IF(first.AND.ifday.EQ.0.AND.tod.EQ.0.0_r8) THEN
       IF (TRIM(iswrad).eq.'LCH') THEN
          !old sw radiation routines
          CALL swrad(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz          , &
               icld        , inalb       , s0          , cos2         , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe       , &
               FlipO3      ,                                            &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam   , &
               ! SW Radiation fields 
               rSwToaDown  ,                                            &
               ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff , &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC    , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet    , &
               asclr       , asl         ,                              &
               ! Cloud field
               cld         , clu         , FlipTaud                       )
       ELSE IF (TRIM(iswrad).eq.'CRD') THEN
          !Clirad-sw radiation routines
          CALL cliradsw(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0          , cos2        , &
               schemes     , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    , psurf       , imask       , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff, &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime )
       ELSE IF (TRIM(iswrad).eq.'CRDTF') THEN
          !Clirad-sw radiation routines
          CALL CliradTarasova_sw(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0          , cos2        , &
               schemes     , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      ,FlipCO2m     , psurf       , imask       , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff, &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime,&
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013) and topography
               ifaeros,aod8,asy8,ssa8,z_aer,topog, &
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013) and topography
               aodF8,asyF8,ssaF8,z_aerF  )
!tar end
!
       ELSE IF (TRIM(iswrad).eq.'RRTMG') THEN
          !RRTMG-sw radiation routines
          CALL Run_Rad_RRTMG_SW(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0          , cos2        , &
               ratio       , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      ,FlipCO2m     , psurf       , imask       , gtg         , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff, &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime, & 
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013) and topography
               ifaeros,aod,asy,ssa,z_aer,topog, &      
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013) 
               aodF,asyF,ssaF,z_aerF  )
!tar end              
 
       ELSE IF (TRIM(iswrad).eq.'UKM') THEN
          !ESRAD-sw radiation routines
          CALL ukmo_swintf(&
               ! Model Info and flags
               ncols       , kmax        ,                             &
               ! Solar field
               cos2        , s0          ,                             &
               ! Atmospheric fields
               FlipPbot    , FlipPmid    , FlipDP      , FlipTe      , &
               FlipQe      , FlipO3      , FlipCO2m    , gtg         , &
               ! SURFACE:  albedo
               imask       ,                                           &
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff, &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    ,               &
               FlipRei     , FlipRel     , FlipLMixR                     )
       ELSE IF (TRIM(iswrad).eq.'NON') THEN
          WRITE(nfprt,*) 'WARN:: Skipping shortwave computation at idate=',idate
       ELSE
          WRITE(nfprt,*) 'ERROR:: WRONG OPTION iswrad=',TRIM(iswrad)
          STOP
       ENDIF

    END IF

    !-----------------------------------------------------------------------
    !     shortwave radiation - R-call = call at last integer hour
    !-----------------------------------------------------------------------

    IF(first.AND.(nfcnv0.EQ.0.OR.(nfcnv0.NE.0.AND.nfin0.EQ.nfin1)))THEN
       IF (TRIM(iswrad).eq.'LCH') THEN
          !old sw radiation routines
          CALL swrad(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz          , &
               icld        , inalb       , s0          , cosz         , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe       , &
               FlipO3      ,                                            &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam   , &
               ! SW Radiation fields 
               rSwToaDown  ,                                            &
               rVisBeam    , rVisDiff    , rNirBeam    , rNirDiff     , &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC    , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet    , &
               asclr       , asl         ,                              &
               ! Cloud field
               cld         , clu         , FlipTaud                       )
       ELSE IF (TRIM(iswrad).eq.'CRD') THEN
          !Clirad-sw radiation routines
          CALL cliradsw(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0          , cosz        , &
               schemes     , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    ,psurf         , imask       , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               rVisBeam    , rVisDiff    , rNirBeam    , rNirDiff    , &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime )
       ELSE IF (TRIM(iswrad).eq.'CRDTF') THEN
          !Clirad-sw radiation routines
          CALL CliradTarasova_sw(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0          , cosz        , &
               schemes     , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    ,psurf         , imask       , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               rVisBeam    , rVisDiff    , rNirBeam    , rNirDiff    , &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime,& 
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013)
               ifaeros,aod8,asy8,ssa8,z_aer,topog, &
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013)
               aodF8,asyF8,ssaF8,z_aerF)
!tar end

       ELSE IF (TRIM(iswrad).eq.'RRTMG') THEN
          !RRTMG-sw radiation routines
          CALL Run_Rad_RRTMG_SW(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0          , cosz        , &
               ratio       , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    ,psurf         , imask       , gtg         , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               rVisBeam    , rVisDiff    , rNirBeam    , rNirDiff    , &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  ,Flipcliqwp   , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime, & 
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013) and topography
               ifaeros,aod,asy,ssa,z_aer,topog, &      
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013) 
               aodF,asyF,ssaF,z_aerF  )
!tar end   
       ELSE IF (TRIM(iswrad).eq.'UKM') THEN
          !ESRAD-sw radiation routines
          CALL ukmo_swintf(&
               ! Model Info and flags
               ncols       , kmax        ,                             &
               ! Solar field
               cosz        , s0          ,                             &
               ! Atmospheric fields
               FlipPbot    , FlipPmid    , FlipDP      , FlipTe      , &
               FlipQe      , FlipO3      , FlipCO2m    , gtg         , &
               ! SURFACE:  albedo
               imask       ,                                           &
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               rSwToaDown  ,                                           &
               rVisBeam    , rVisDiff    , rNirBeam    , rNirDiff    , &
               rVisBeamC   , rVisDiffC   , rNirBeamC   , rNirDiffC   , &
               rSwToaNetC  , rSwToaNet   , rSwSfcNetC  , rSwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    ,               &
               FlipRei     , FlipRel     , FlipLMixR                     )
       ELSE IF (TRIM(iswrad).eq.'NON') THEN
          IF(myid==0)WRITE(nfprt,*) 'WARN:: Skipping shortwave computation at idate=',idate
       ELSE
          IF(myid==0)WRITE(nfprt,*) 'ERROR:: WRONG OPTION iswrad=',TRIM(iswrad)
          STOP
       ENDIF

    END IF

    !-----------------------------------------------------------------------
    ! Flip microphysics variables
    !-----------------------------------------------------------------------

    DO k=1,kmax
       kflip=kmax+1-k
       DO i = 1,ncols
          Taud (i,kflip)=FlipTaud (i,k)
       END DO
    END DO

    !-----------------------------------------------------------------------
    !     shortwave radiation - Y-call = call at next integer hour
    !-----------------------------------------------------------------------

    IF (TRIM(iswrad).eq.'LCH') THEN
       !old sw radiation routines
          CALL swrad(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz          , &
               icld        , inalb       , s0x         , cos1         , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe       , &
               FlipO3      ,                                            &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam   , &
               ! SW Radiation fields 
               ySwToaDown  ,                                            &
               yVisBeam    , yVisDiff    , yNirBeam    , yNirDiff     , &
               yVisBeamC   , yVisDiffC   , yNirBeamC   , yNirDiffC    , &
               ySwToaNetC  , ySwToaNet   , ySwSfcNetC  , ySwSfcNet    , &
               asclr       , asl         ,                              &
               ! Cloud field
               cld         , clu         , FlipTaud                       )
    ELSE IF (TRIM(iswrad).eq.'CRD') THEN
       !Clirad-sw radiation routines
          CALL cliradsw(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0x         , cos1        , &
               schemes     , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    ,psurf         , imask       , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               ySwToaDown  ,                                           &
               yVisBeam    , yVisDiff    , yNirBeam    , yNirDiff    , &
               yVisBeamC   , yVisDiffC   , yNirBeamC   , yNirDiffC   , &
               ySwToaNetC  , ySwToaNet   , ySwSfcNetC  , ySwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud   , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime )
    ELSE IF (TRIM(iswrad).eq.'CRDTF') THEN
       !Clirad-sw radiation routines
          CALL CliradTarasova_sw(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0x         , cos1        , &
               schemes     , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    , psurf         , imask       , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               ySwToaDown  ,                                           &
               yVisBeam    , yVisDiff    , yNirBeam    , yNirDiff    , &
               yVisBeamC   , yVisDiffC   , yNirBeamC   , yNirDiffC   , &
               ySwToaNetC  , ySwToaNet   , ySwSfcNetC  , ySwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  , Flipcliqwp  , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime,&
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013)
               ifaeros,aod8,asy8,ssa8,z_aer,topog, &
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013)
               aodF8,asyF8,ssaF8,z_aerF)      
!tar end

       ELSE IF (TRIM(iswrad).eq.'RRTMG') THEN
          !RRTMG-sw radiation routines
          CALL Run_Rad_RRTMG_SW(&
               ! Model Info and flags
               ncols       , kmax        , nls         , noz         , &
               icld        , inalb       , s0x         , cos1        , &
               ratio       , &
               ! Atmospheric fields
               FlipPbot    , FlipDP      , FlipTe      , FlipQe      , &
               FlipO3      , FlipCO2m    ,psurf         , imask       , gtg         , &
               ! SURFACE:  albedo
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               ySwToaDown  ,                                           &
               yVisBeam    , yVisDiff    , yNirBeam    , yNirDiff    , &
               yVisBeamC   , yVisDiffC   , yNirBeamC   , yNirDiffC   , &
               ySwToaNetC  , ySwToaNet   , ySwSfcNetC  , ySwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    , &
               FlipRei     , FlipRel     , FlipTaud    , &
               Flipcicewp  ,Flipcliqwp   , &
               Flip_c_cld_tau    ,Flip_c_cld_tau_w     , &
               Flip_c_cld_tau_w_g,Flip_c_cld_tau_w_f   ,Flip_cldfprime, &  
!tar begin
! Climate aerosol parameters of coarse mode (Kinne, 2013) and topography
               ifaeros,aod,asy,ssa,z_aer,topog, &      
!tar end
!
!tar begin
! Climate aerosol parameters of fine mode (Kinne, 2013) 
               aodF,asyF,ssaF,z_aerF  )
!tar end

    ELSE IF (TRIM(iswrad).eq.'UKM') THEN
          !ESRAD-sw radiation routines
          CALL ukmo_swintf(&
               ! Model Info and flags
               ncols       , kmax        ,                             &
               ! Solar field
               cos1        , s0x         ,                             &
               ! Atmospheric fields
               FlipPbot    , FlipPmid    , FlipDP      , FlipTe      , &
               FlipQe      , FlipO3      , FlipCO2m    , gtg         , &
               ! SURFACE:  albedo
               imask       ,                                           &
               AlbVisDiff  , AlbNirDiff  , AlbVisBeam  , AlbNirBeam  , &
               ! SW Radiation fields 
               ySwToaDown  ,                                           &
               yVisBeam    , yVisDiff    , yNirBeam    , yNirDiff    , &
               yVisBeamC   , yVisDiffC   , yNirBeamC   , yNirDiffC   , &
               ySwToaNetC  , ySwToaNet   , ySwSfcNetC  , ySwSfcNet   , &
               asclr       , asl         ,                             &
               ! Cloud field
               cld         , clu         , FlipFice    ,               &
               FlipRei     , FlipRel     , FlipLMixR                     )
    ELSE IF (TRIM(iswrad).eq.'NON') THEN
       WRITE(nfprt,*) 'WARN:: Skipping shortwave computation at idate=',idate
    ELSE
       WRITE(nfprt,*) 'ERROR:: WRONG OPTION iswrad=',TRIM(iswrad)
       STOP
    ENDIF

    DO k=1,kmax
       DO i=1,ncols
          ySwHeatRate(i,k)=asl(i,kmax+1-k)    ! heating rate
          ySwHeatRatec(i,k)=asclr(i,kmax+1-k) ! heating rate clear
       END DO
    END DO
       !CALL MsgOne(h,'end shortwave')

       !CALL MsgOne(h,'begin longwave')

    !-----------------------------------------------------------------------
    !     longwave radiation
    !-----------------------------------------------------------------------
!PK    IF(ABS( MOD((tod-delt)/3.6e3_r8+0.03125e0_r8,trint)).LE.0.0625e0_r8 &
!PK         .OR.first)THEN
    IF(ABS( MOD((tod-delt)+0.03125e0_r8,trint)).LE.0.0625e0_r8 &
         .OR.first)THEN
       IF (TRIM(ilwrad).eq.'HRS') THEN
          ! New shortwave schemes need cld and clu separately, but old
          ! longwave scheme need them maximum overlaped. If you look
          ! above, you'll see that for CRD and UKM the new cldgn3()
          ! outputs both clu and cld... Therefore, before calling old
          ! longwave, we restore the state of the cloud field.
          if (mxrdcc) then
             cld = MAX(cld,clu)
             clu = 0.0_r8
          else
             clu = 0.0_r8
          endif
          
          CALL lwrad( &
               ! Model Info and flags
               ncols     , kmax     , nls     , nlcs    , noz      , icld    , &
               ! Atmospheric fields
               FlipPbot  , FlipPmid , FlipTe  , FlipQe  , FlipO3   , gtg     , &
               FlipCO2m,                                                &
               ! LW Radiation fields 
               LwToaUpC  , LwToaUp  , asclr   , asl     , LwSfcNetC, LwSfcNet, &
               LwSfcDownC, LwSfcDown,                                          &
               ! Cloud field and Microphysics
               cld      , clu     , Flipclwp  , FlipFice, FlipRei  , emisd       )
       ELSE IF (TRIM(ilwrad).eq.'CRD') THEN 
          CALL cliradlw ( &
               ! Model Info and flags
               ncols     , kmax      ,icld, &
              ! SURFACE:  albedo
               imask       ,                                           &
               ! Atmospheric fields
               FlipPmid  , FlipPInterface    , FlipTe    , FlipQe    , FlipO3    , &
               FlipCO2m  , gtg       , &   
               ! LW Radiation fields 
               LwToaUpC  , LwToaUp   , asclr     , asl       , LwSfcNetC , LwSfcNet  , &
               LwSfcDownC, LwSfcDown , &
               ! Cloud field and Microphysics
               cld       , clu        , FlipFice  ,             &
               FlipRei   , FlipRel  , FlipLMixR,Flip_c_cld_lw_abs,Flipcicewp  ,Flipcliqwp ,Flip_cldfprime  )

        !  WRITE(nfprt,*) 'ERROR:: CLIRAD_LW NOT YET IMPLEMENTED!'
          
       ELSE IF (TRIM(ilwrad).eq.'CRDTF') THEN 
          CALL CliradTarasova_lw ( &
               ! Model Info and flags
               ncols     , kmax      ,icld, &
               ! Atmospheric fields
               FlipPmid  , FlipPInterface    , FlipTe    , & 
               FlipQe    , FlipO3            , FlipCO2m  , gtg       , & 
               ! LW Radiation fields 
               LwToaUpC  , LwToaUp   , asclr     , asl       , &
               LwSfcNetC , LwSfcNet  , LwSfcDownC, LwSfcDown , &
               ! Cloud field and Microphysics
               cld       , clu      , FlipFice  ,             &
               FlipRei   , FlipRel  , FlipLMixR,Flip_c_cld_lw_abs,Flipcicewp  ,Flipcliqwp ,Flip_cldfprime  )

          !WRITE(nfprt,*) 'ERROR:: CLIRAD_LW NOT YET IMPLEMENTED!'
          

       ELSE IF (TRIM(ilwrad).eq.'RRTMG') THEN
          CALL Run_Rad_RRTMG_LW( &
               ! Model Info and flags
               ncols     , kmax     ,                         &
               ! Atmospheric fields
               FlipPbot  , FlipPmid , FlipDP    , FlipTe   ,  &
               FlipQe    , FlipO3   , FlipCO2m    , gtg      ,  &
               psurf           ,&
               ! SURFACE
               imask     ,                                    &
               ! LW Radiation fields 
               LwToaUpC  , LwToaUp  , asclr     , asl      ,  &
               LwSfcNetC , LwSfcNet , LwSfcDownC, LwSfcDown,  &
               ! Cloud field and Microphysics
               cld       , clu      , FlipFice  ,             &
               FlipRei   , FlipRel  , FlipLMixR,Flip_c_cld_lw_abs,Flipcicewp  ,Flipcliqwp ,Flip_cldfprime  )
       ELSE IF (TRIM(ilwrad).eq.'UKM') THEN
          CALL ukmo_lwintf( &
               ! Model Info and flags
               ncols     , kmax     ,                         &
               ! Atmospheric fields
               FlipPbot  , FlipPmid , FlipDP    , FlipTe   ,  &
               FlipQe    , FlipO3   , FlipCO2m  , gtg      ,  &
               ! SURFACE
               imask     ,                                    &
               ! LW Radiation fields 
               LwToaUpC  , LwToaUp  , asclr     , asl      ,  &
               LwSfcNetC , LwSfcNet , LwSfcDownC, LwSfcDown,  &
               ! Cloud field and Microphysics
               cld       , clu      , FlipFice  ,             &
               FlipRei   , FlipRel  , FlipLMixR                 )

       ENDIF

       IF (TRIM(ilwrad).ne.'NON') THEN
          DO k=1,kmax
             DO i = 1,ncols
                LwCoolRate (i,k)=asl  (i,kmax+1-k)    ! cooling rate
                LwCoolRatec(i,k)=asclr(i,kmax+1-k) ! cooling rate clear
             END DO
          END DO
       ENDIF

    ENDIF
       !CALL MsgOne(h,'end longwave')

  END SUBROUTINE spmrad

  SUBROUTINE cldgn2 ( &
       covlp ,prsl,prsi,grh   ,omg   ,gtmp  ,css   ,ccu   , &
       cdin  ,cstc  ,ccon  ,cson  ,mxrdcc,lcnvl ,convc ,convt ,convb , &
       convts,convcs,convbs,ncols ,kmax  ,nls                  )

    ! parameters and input variables:                                     c
    !        covlp = 'maxi'         maximum overlap of convective cloud   c
    !                               or thick low cloud used in ir subr.   c
    !              = 'rand'         random  overlap of convective cloud   c
    !                               or thick low cloud used in ir subr.   c
    !        date  =  julian day of model forecast date                   c
    !        grh   =  relative humidity  (0-1)                            c
    !        omg   =  vertical velocity  (cb/sec)                         c
    !        prsl  =  pressure   (Pa)                                     c
    !        sigmid   =  sigma coordinate at middle of layer              c
    !        gtmp  =  layer temperature (k)                               c
    ! output variables:                                                   c
    !        css   =  ncols*kmax supersatuation cloud cover fraction      c
    !        ccu   =  ncols*kmax convective cloud cover fraction          c
    !
    ! values from subr-gwater:                                            c
    !       convc  =  ncols convective cloud cover in 3 hr. avrage        c
    !       convt  =  ncols convective cloud top  (sigma layer)           c
    !       convb  =  ncols convective cloud base (sigma layer)           c
    !==========================================================================
    !
    !
    !   ncols.....Number of grid points on a gaussian latitude circle
    !   kmax......Number of grid points at vertical
    !   nls..... .Number of layers in the stratosphere.
    !   cdin
    !   cstc......cstc=clow change necessary in order to properly mark inv
    !             cloud height
    !   ccon......convc  =  ncols convective cloud cover in 3 hr. avrage
    !   cson
    !   cp........Specific heat of air (j/kg/k)
    !   gasr......Constant of dry air      (j/kg/k)
    !   mxrdcc....use maximum random converage for radiative conv. clouds
    !             constant logical mxrdcc = .true.
    !   lcnvl.....the lowest layer index where non-convective clouds can
    !             occur (ben says this should be 2 or more)
    !             constant lcnvl = 2
    !   convts
    !   convcs
    !   convbs
    !
    !
    INTEGER,          INTENT(in ) :: ncols
    INTEGER,          INTENT(in ) :: kmax
    INTEGER,          INTENT(in ) :: nls
    CHARACTER(len=4), INTENT(in ) :: covlp
    REAL(KIND=r8),             INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),             INTENT(in   ) :: prsi  (ncols,kMax+1)
    REAL(KIND=r8),             INTENT(in   ) :: grh (ncols,kmax)
    REAL(KIND=r8),             INTENT(in   ) :: omg (ncols,kmax)
    REAL(KIND=r8),             INTENT(in   ) :: gtmp(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: css (ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: ccu (ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cdin(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cstc(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: ccon(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cson(ncols,kmax)

    LOGICAL,          INTENT(in ) :: mxrdcc
    INTEGER,          INTENT(in ) :: lcnvl
    REAL(KIND=r8),             INTENT(in ) :: convc (ncols)
    REAL(KIND=r8),             INTENT(in ) :: convt (ncols)
    REAL(KIND=r8),             INTENT(in ) :: convb (ncols)

    REAL(KIND=r8),             INTENT(in ) :: convts(ncols)
    REAL(KIND=r8),             INTENT(in ) :: convcs(ncols)
    REAL(KIND=r8),             INTENT(in ) :: convbs(ncols)

    REAL(KIND=r8)                          :: pll   (ncols,kmax)
    REAL(KIND=r8)                          :: pii   (ncols,kmax+1)
    REAL(KIND=r8)                          :: conv  (ncols)
    REAL(KIND=r8)                          :: clow  (ncols)
    REAL(KIND=r8)                          :: cmid  (ncols)
    REAL(KIND=r8)                          :: chigh (ncols)
    INTEGER                       :: ktop  (ncols)
    INTEGER                       :: kbot  (ncols)
    INTEGER                       :: klow  (ncols)
    INTEGER                       :: kmid  (ncols)
    INTEGER                       :: khigh (ncols)
    REAL(KIND=r8)                          :: pt    (ncols)
    REAL(KIND=r8)                          :: cinv  (ncols)
    REAL(KIND=r8)                          :: delomg(ncols)
    REAL(KIND=r8)                          :: csat  (ncols,kmax)
    INTEGER                       :: ktops (ncols)
    INTEGER                       :: kbots (ncols)
    REAL(KIND=r8)                          :: convs (ncols)
    REAL(KIND=r8)                          :: convh (ncols)
    INTEGER                       :: invb  (ncols)
    REAL(KIND=r8)                          :: dthdpm(ncols)
    INTEGER                       :: i
    INTEGER                       :: k
    INTEGER                       :: kl
    INTEGER                       :: lon
    REAL(KIND=r8)                          :: zrth
    REAL(KIND=r8)                          :: zrths
    REAL(KIND=r8)                          :: arcp
    REAL(KIND=r8)                          :: dthdp
    !
    !     dthdpm-----min (d(theta)/d(p))
    !     invb ------the base of inversion layer
    !     dthdpc-----criterion of (d(theta)/d(p))
    !
    REAL(KIND=r8), PARAMETER               :: dthdpc = -0.4e-1_r8

    REAL(KIND=r8), PARAMETER               :: f6p67= 6.67e0_r8
    REAL(KIND=r8), PARAMETER               :: f400p= 4.0e2_r8
    REAL(KIND=r8), PARAMETER               :: f700p= 8.0e2_r8!7.5e2_r8!7.0e2_r8
    REAL(KIND=r8), PARAMETER               :: f0p6 = 0.6e0_r8
    REAL(KIND=r8), PARAMETER               :: f0p2 = 0.2e0_r8
    REAL(KIND=r8), PARAMETER               :: f0p8 = 0.8e0_r8
    REAL(KIND=r8), PARAMETER               :: f5m5 = 5.0e-5_r8
    REAL(KIND=r8), PARAMETER               :: f1e4 = 1.0e+4_r8


    DO k = 1, kmax
       DO i = 1, ncols
          css (i,k) = 0.0e0_r8
          ccu (i,k) = 0.0e0_r8
          csat(i,k) = 0.0e0_r8
          cdin(i,k) = 0.0e0_r8
          cstc(i,k) = 0.0e0_r8
          ccon(i,k) = 0.0e0_r8
          cson(i,k) = 0.0e0_r8
       END DO
    END DO
    !
    !     the clouds generation scheme is based on j. slingo              c
    !     (1984 ecmwf workshop).  the scheme generates 4 type of clouds   c
    !     of convective, high, middle and low clouds.                     c
    !
    DO kl = 1, kmax
       DO lon = 1,ncols
          pll(lon,kl) = prsl(lon,kl)/100.0_r8!gps(lon) * sigmid(kl)
       END DO
    END DO
    DO kl = 1, kmax+1
       DO lon = 1,ncols
          pii(lon,kl) = prsi(lon,kl)/100.0_r8!gps(lon) * sigmid(kl)
       END DO
    END DO
    !
    !     initialization
    !
    DO lon = 1,ncols
       conv(lon)  = convc(lon)
       convh(lon) = convc(lon)
       clow(lon)  = 0.0_r8
       cinv(lon)  = 0.0_r8
       cmid(lon)  = 0.0_r8
       chigh(lon) = 0.0_r8
       ktop(lon)  = INT(convt(lon)+0.5_r8)
       IF (ktop(lon).LT.1.OR.ktop(lon).GT.kmax) ktop(lon)=1
       kbot(lon)  = INT(convb(lon)+0.5_r8)
       IF (kbot(lon).LT.1.OR.kbot(lon).GT.kmax) kbot(lon)=1
       klow(lon)  = 1
       kmid(lon)  = 1
       khigh(lon) = 1
       convs(lon) = convcs(lon)
       ktops(lon) = INT(convts(lon)+0.5_r8)
       IF (ktops(lon).LT.1.OR.ktops(lon).GT.kmax) ktops(lon)=1
       kbots(lon) = INT(convbs(lon)+0.5_r8)
       IF (kbots(lon).LT.1.OR.kbots(lon).GT.kmax) kbots(lon)=1
       !
       !     1. define convective cloud ***   done in subr-gwater
       !     cloud top and base are defined by kuo scheme: convt, convb
       !     cloud amount is calculated from precipitation rate : convc
       !     single layer clouds conputations start here, from bottom up
       !
       !     define high clouds due to strong convection
       !
       pt(lon) = pll(lon,ktop(lon)) !gps(lon) * sigmid(ktop(lon))
       IF ((pt(lon) <= f400p).AND. (conv(lon) > 0.0_r8)) THEN
          chigh(lon) = 2.0_r8*convh(lon)
          chigh(lon) = MIN(chigh(lon),1.0_r8)
          khigh(lon) = ktop(lon) + 1
       END IF
       IF (ktop(lon)-kbot(lon) >= 1) THEN
          zrth=1.0_r8/REAL(ktop(lon)-kbot(lon),r8)
          conv(lon)=1.0_r8-(1.0_r8-conv(lon))**zrth
       END IF
       IF (ktops(lon)-kbots(lon) >= 1) THEN
          zrths=1.0_r8/REAL(ktops(lon)-kbots(lon),r8)
          convs(lon)=1.0_r8-(1.0_r8-convs(lon))**zrths
       END IF
    END DO
    !
    !     compute low stratus associated with inversions, based on ecwmf's
    !     scheme, with lower criterion of d(theta)/d(p)
    !
    arcp = gasr / cp

    DO lon = 1, ncols
       invb(lon) = MIN(kmax,kmax-nls)
       dthdpm(lon) = 0.0_r8
    END DO

    DO kl = 2, kmax
       DO lon = 1, ncols
          IF (pll(lon,kl) > f700p) THEN
             dthdp = (gtmp(lon,kl-1)*(1000.0_r8/pll(lon,kl-1))**arcp &
                     -gtmp(lon,kl  )*(1000.0_r8/pll(lon,kl  ))**arcp)/ &
                  (pll(lon,kl-1) - pll(lon,kl))
             IF (dthdp < 0.0_r8) THEN
                invb(lon) = MIN(kl-1,invb(lon))
                IF(dthdp.LT.dthdpc) THEN
                   IF(dthdp.LT.dthdpm(lon)) THEN
                      dthdpm(lon)=dthdp
                      klow(lon)=kl
                   END IF
                END IF
             END IF
          END IF
       END DO
    END DO
    !
    !     klow change above necessary to mark inversion cloud height
    !
    DO lon = 1, ncols
       IF (dthdpm(lon) < dthdpc .AND. grh(lon,invb(lon)) > f0p6) THEN
          cinv(lon) = - f6p67 * (dthdpm(lon)-dthdpc)
          cinv(lon) = MAX(cinv(lon),0.0_r8)
          cinv(lon) = MIN(cinv(lon),1.0_r8)
          IF (grh(lon,invb(lon)) < f0p8) &
               cinv(lon) = cinv(lon)* &
               (1.0_r8-(f0p8 - grh(lon,invb(lon)))/f0p2)
       END IF
       clow(lon)=cinv(lon)
       IF (conv(lon) <= 0.0_r8) THEN
          cdin(lon,klow(lon))=cinv(lon)
       END IF
    END DO
    !
    !     clow change necessary in order to properly mark inv cloud height
    !
    !     main loop for cloud amount determination
    !
    DO kl = lcnvl+1, MIN(kmax,kmax-nls)
       DO lon = 1, ncols
          !
          !     general define cloud due to saturation
          !
          csat(lon,kl)=(grh(lon,kl) - 0.9_r8) / .1
          csat(lon,kl)=(MAX(csat(lon,kl), 0.0_r8)) ** 2
          IF (pll(lon,kl) < 700.0_r8) THEN
             IF (omg(lon,kl) >= f5m5) THEN
                csat(lon,kl)=0.0_r8
             ELSE IF (omg(lon,kl) >= -f5m5) THEN
                delomg(lon) = (omg(lon,kl)+f5m5)*f1e4
                csat(lon,kl)=csat(lon,kl)*  &
                     (1.0_r8 - delomg(lon) * delomg(lon))
             END IF
          END IF
          csat(lon,kl)=MIN(csat(lon,kl), 1.0_r8)
          cstc(lon,kl)=csat(lon,kl)
       END DO
    END DO
    DO lon = 1, ncols
       DO kl = 1, kmax
          IF (cdin(lon,kl) > 0.0_r8) THEN
             css(lon,kmax+1-kl) = cdin(lon,kl)
          ELSE
             css(lon,kmax+1-kl) = csat(lon,kl)
          END IF
       END DO
    END DO
    !
    !     convective cloud is maximum overlaping in ir subr.
    !
    IF (covlp == 'MAXI') THEN
       DO lon = 1, ncols
          DO kl = kbot(lon), ktop(lon)
             ccu(lon,kmax-kl+1) = conv(lon)
          END DO
       END DO
       !
       !     convective cloud is random overlaping in ir subr.
       !
    ELSE IF (covlp == 'RAND') THEN

       DO lon = 1, ncols
          DO kl = kbots(lon), ktops(lon)
             IF ((conv(lon) <= 0.0_r8) .AND. (cdin(lon,kl) >= 0.0_r8)) THEN
                ccon(lon,kl)=convs(lon)
                cson(lon,kl)=convs(lon)
             END IF
          END DO
       END DO

       DO lon = 1, ncols
          DO kl = 1, kmax
             IF ((kl >= kbot(lon)) .AND. (kl <= ktop(lon))) THEN
                ccon(lon,kl)=ccon(lon,kl)+conv(lon)
             ELSE
                ccon(lon,kl)=ccon(lon,kl)
             END IF
             ccon(lon,kl)=MIN(ccon(lon,kl),1.0_r8)
             css(lon,kmax-kl+1) = MIN(1.0_r8,css(lon,kmax-kl+1))
          END DO
       END DO

       DO lon = 1, ncols
          ccon(lon,khigh(lon))=ccon(lon,khigh(lon))+chigh(lon)
       END DO

       DO kl = 1, kmax
          DO lon = 1, ncols
             IF(mxrdcc)THEN
                css(lon,kmax-kl+1) =  &
                     (1.0_r8-ccon(lon,kl))*css(lon,kmax-kl+1) + &
                     ccon(lon,kl)
                css(lon,kmax-kl+1) = MIN(1.0_r8,css(lon,kmax-kl+1))
             END IF
          END DO
       END DO

    END IF
  END SUBROUTINE cldgn2






  SUBROUTINE cldgen (covlp ,prsl,prsi,grh   ,omg   ,gtmp  ,css   , &
       ccu   ,cdin  ,cstc  ,ccon  ,cson  ,mxrdcc,lcnvl ,lthncl, &
       convc ,convt ,convb ,ncols ,kmax  ,nls            )

    !==========================================================================
    ! cldgen :perform clouds generation scheme based on j. slingo
    !         (1984 ecmwf workshop); the scheme generates 4 type of clouds
    !         of convective, high, middle and low clouds.
    !==========================================================================
    ! parameters and input variables:
    !        covlp = 'maxi'         maximum overlap of convective cloud
    !                               or thick low cloud used in ir subr.
    !              = 'rand'         random  overlap of convective cloud
    !                               or thick low cloud used in ir subr.
    !        date  =  julian day of model forecast date
    !        grh   =  relative humidity  (0-1)
    !        omg   =  vertical velocity  (cb/sec)
    !        prsl  =  pressure   (Pa)
    !        sigmid   =  sigma coordinate at middle of layer
    !        gtmp  =  layer temperature (k)
    !     output variables:
    !        css   =  ncols*kmax supersatuation cloud cover fraction
    !        ccu   =  ncols*kmax convective cloud cover fraction
    !---------------------------------------------------------------------
    ! values from subr-gwater
    !       convc  =  ncols*kmax convective cloud cover in 3 hr. avrage
    !       convt  =  ncols*kmax convective cloud top  (sigma layer)
    !       convb  =  ncols*kmax convective cloud base (sigma layer)
    !       prcp1,prcp2,prcp3,prcpt,toplv,botlv: are used in subr "gwater"
    !==========================================================================
    !
    !   ncols......Number of grid points on a gaussian latitude circle
    !   kmax......Number of grid points at vertical
    !   nls..... .Number of layers in the stratosphere.
    !   cdin
    !   cstc......cstc=clow change necessary in order to properly mark inv
    !             cloud height
    !   ccon......convc  =  ncols*kmax convective cloud cover in 3 hr. avrage
    !   cson
    !   cp........Specific heat of air (j/kg/k)
    !   gasr......Constant of dry air      (j/kg/k)
    !   mxrdcc....use maximum random converage for radiative conv. clouds
    !             constant logical mxrdcc = .true.
    !   lcnvl.....the lowest layer index where non-convective clouds can
    !             occur (ben says this should be 2 or more)
    !             constant lcnvl = 2
    !   lthncl....Minimum depth in mb of non-zero low level cloud
    !             consta lthncl=80
    !==========================================================================
    INTEGER,          INTENT(IN   ) :: ncols
    INTEGER,          INTENT(in ) :: kmax
    INTEGER,          INTENT(in ) :: nls

    CHARACTER(len=4), INTENT(in   ) :: covlp
    REAL(KIND=r8),             INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),             INTENT(in   ) :: prsi  (ncols,kMax+1)
!    REAL(KIND=r8),             INTENT(in   ) :: gps (ncols)
!    REAL(KIND=r8),             INTENT(in   ) :: sigmid (kmax)
    REAL(KIND=r8),             INTENT(in   ) :: grh (ncols,kmax)
    REAL(KIND=r8),             INTENT(in   ) :: omg (ncols,kmax)
    REAL(KIND=r8),             INTENT(in   ) :: gtmp(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: css (ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: ccu (ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cdin(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cstc(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: ccon(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cson(ncols,kmax)

    LOGICAL      ,             INTENT(in   ) :: mxrdcc
    INTEGER      ,             INTENT(in   ) :: lcnvl
    INTEGER      ,             INTENT(in   ) :: lthncl
    REAL(KIND=r8),             INTENT(in   ) :: convc(ncols)
    REAL(KIND=r8),             INTENT(in   ) :: convt(ncols)
    REAL(KIND=r8),             INTENT(in   ) :: convb(ncols)
    !
    !     dthdpm-----min (d(theta)/d(p))
    !     invb ------the base of inversion layer
    !     dthdpc-----criterion of (d(theta)/d(p))
    !
    REAL(KIND=r8),    PARAMETER :: dthdpc = -0.4e-1_r8
    REAL(KIND=r8)    :: pll   (ncols,kmax)
    REAL(KIND=r8)    :: pii   (ncols,kmax+1)
    REAL(KIND=r8)    :: conv  (ncols)
    REAL(KIND=r8)    :: clow  (ncols)
    REAL(KIND=r8)    :: cmid  (ncols)
    REAL(KIND=r8)    :: chigh (ncols)
    INTEGER          :: ktop  (ncols)
    INTEGER          :: kbot  (ncols)
    INTEGER          :: klow  (ncols)
    INTEGER          :: kmid  (ncols)
    INTEGER          :: khigh (ncols)
    REAL(KIND=r8)    :: pt    (ncols)
    REAL(KIND=r8)    :: cx    (ncols)
    REAL(KIND=r8)    :: cinv  (ncols)
    REAL(KIND=r8)    :: delomg(ncols)
    REAL(KIND=r8)    :: csat  (ncols,kmax)
    REAL(KIND=r8)    :: dthdpm(ncols)
    INTEGER          :: invb  (ncols)

    REAL(KIND=r8), PARAMETER :: f700p= 8.0e2_r8!7.5e2_r8!7.0e2_r8
    REAL(KIND=r8), PARAMETER :: f400p= 4.0e2_r8
    REAL(KIND=r8), PARAMETER :: f6p67= 6.67e0_r8
    REAL(KIND=r8), PARAMETER :: f0p9 = 0.9e0_r8
    REAL(KIND=r8), PARAMETER :: f0p8 = 0.8e0_r8
    REAL(KIND=r8), PARAMETER :: f0p4 = 0.4e0_r8
    REAL(KIND=r8), PARAMETER :: f0p3 = 0.3e0_r8
    REAL(KIND=r8), PARAMETER :: f0p2 = 0.2e0_r8
    REAL(KIND=r8), PARAMETER :: f0p6 = 0.6e0_r8
    REAL(KIND=r8), PARAMETER :: f1e4 = 1.0e+4_r8
    REAL(KIND=r8), PARAMETER :: f5m5 = 5.0e-5_r8

    INTEGER :: i
    INTEGER :: k
    INTEGER :: kl
    INTEGER :: lon
    
    REAL(KIND=r8)    :: arcp
    REAL(KIND=r8)    :: dthdp
    REAL(KIND=r8)    :: thklow
    pll    = 0.0_r8
    conv   = 0.0_r8
    clow   = 0.0_r8
    cmid   = 0.0_r8
    chigh  = 0.0_r8
    ktop   = 0
    kbot   = 0
    klow   = 0
    kmid   = 0
    khigh  = 0
    pt     = 0.0_r8
    cx     = 0.0_r8
    cinv   = 0.0_r8
    delomg = 0.0_r8
    csat   = 0.0_r8
    dthdpm = 0.0_r8
    invb   = 0
    arcp = 0.0_r8
    dthdp = 0.0_r8
    thklow = 0.0_r8
    DO k = 1, kmax
       DO i = 1, ncols
          css(i,k) = 0.0_r8
          ccu(i,k) = 0.0_r8
          csat(i,k) = 0.0_r8
          cdin(i,k) = 0.0_r8
          cstc(i,k) = 0.0_r8
          ccon(i,k) = 0.0_r8
          cson(i,k) = 0.0_r8
       END DO
    END DO
    !
    !     the clouds generation scheme is based on j. slingo
    !     (1984 ecmwf workshop).  the scheme generates 4 type of clouds
    !     of convective, high, middle and low clouds.
    !
    DO kl = 1, kmax
       DO lon = 1,ncols
          pll(lon,kl) = prsl(lon,kl)/100.0_r8 !gps(lon) * sigmid(kl)
       END DO
    END DO
    DO kl = 1, kmax+1
       DO lon = 1,ncols
          pii(lon,kl) = prsi(lon,kl)/100.0_r8!gps(lon) * sigmid(kl)
       END DO
    END DO
    !
    !     initialization
    !
    DO  lon = 1,ncols
       conv(lon)  = convc(lon)
       clow(lon)  = 0.0_r8
       cinv(lon)  = 0.0_r8
       cmid(lon)  = 0.0_r8
       chigh(lon) = 0.0_r8
       ktop(lon)  = INT(convt(lon)+0.5_r8)
       IF (ktop(lon).LT.1.OR.ktop(lon).GT.kmax) ktop(lon)=1
       kbot(lon)  = INT(convb(lon)+0.5_r8)
       IF (kbot(lon).LT.1.OR.kbot(lon).GT.kmax) kbot(lon)=1
       klow(lon)  = 1
       kmid(lon)  = 1
       khigh(lon) = 1
       !
       !     1. define convective cloud
       !     done in subr-gwater
       !     cloud top and base are defined by kuo scheme: convt, convb
       !     cloud amount is calculated from precipitation rate : convc
       !     *** single layer clouds conputations start here, from bottom up
       !
       !     define high clouds due to strong convection
       !
       pt(lon) =  pll(lon,ktop(lon)) !gps(lon) * sigmid(ktop(lon))
       IF ((pt(lon).LE.f400p).AND.&
            (conv(lon).GE.f0p4)) THEN
          chigh(lon) = 2.0_r8 * (conv(lon) - f0p3)
          chigh(lon) = MIN(chigh(lon),1.0_r8)
          khigh(lon) = ktop(lon) + 1
          ccon(lon,khigh(lon))=chigh(lon)
       END IF
    END DO
    !
    !     compute low stratus associated with inversions, based on ecwmf's
    !     scheme, with lower criterion of d(theta)/d(p)
    !
    arcp = gasr / cp

    DO lon=1,ncols
       invb(lon) =  MIN(kmax, kmax-nls)
       dthdpm(lon) = 0.0_r8
    END DO

    DO kl=2,kmax
       DO lon=1,ncols
          IF (pll(lon,kl) .GT. f700p) THEN
             dthdp = (gtmp(lon,kl-1)*(1000.0_r8/pll(lon,kl-1))**arcp &
                     -gtmp(lon,kl  )*(1000.0_r8/pll(lon,kl  ))**arcp)/&
                  (pll(lon,kl-1) - pll(lon,kl))
             IF(dthdp.LT.0.0_r8) THEN
                invb(lon) = MIN(kl-1,invb(lon))
                IF(dthdp.LT.dthdpc) THEN
                   IF(dthdp.LT.dthdpm(lon)) THEN
                      dthdpm(lon)=dthdp
                      klow(lon)=kl
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       END DO
    END DO
    !
    !     klow change above necessary to mark inversion cloud height
    !
    DO lon=1,ncols
       IF(dthdpm(lon).LT.dthdpc.AND.grh(lon,invb(lon)).GT.f0p6)THEN
          cinv(lon) = - f6p67 * (dthdpm(lon)-dthdpc)
          cinv(lon) = MAX(cinv(lon),0.0_r8)
          cinv(lon) = MIN(cinv(lon),1.0_r8)
          IF (grh(lon,invb(lon)) .LT. f0p8)&
               cinv(lon) = cinv(lon)*&
               (1.0_r8-(f0p8 - grh(lon,invb(lon)))/f0p2)
       ENDIF
       clow(lon)=cinv(lon)
       cdin(lon,klow(lon))=cinv(lon)
    END DO
    !
    !     clow change necessary in order to properly mark inv cloud height
    !
    !     main loop for cloud amount determination
    !
    DO kl = lcnvl+1, MIN(kmax, kmax-nls)
       DO lon = 1,ncols
          !
          !     general define cloud due to saturation
          !
          IF (pll(lon,kl) .GT. f400p) THEN
             cx(lon) = (grh(lon,kl) - f0p8) / f0p2
          ELSE
             cx(lon) = (grh(lon,kl) - f0p9) / 0.1_r8
          END IF
          cx(lon) = (MAX(cx(lon), 0.0_r8)) ** 2
          cx(lon) =  MIN(cx(lon), 1.0_r8)
          !
          !     start vertical process from bottom to top
          !
          IF (pll(lon,kl) .GT. f700p) THEN
             !
             !     2. define low cloud ***
             !     low cloud is defined one layer thick ranging from layer 3 to 700mb
             !     there are two type possible generating mechanisms. due boundary
             !     t inversion type and associated with vertical motion.
             !
             !     define low super satuated clouds but adjusted by vertical motion
             !
             IF (omg(lon,kl) .GE. f5m5) THEN
                cx(lon) = 0.0_r8
             ELSE IF (omg(lon,kl) .GE. -f5m5) THEN
                delomg(lon) = (omg(lon,kl)+f5m5)*f1e4
                cx(lon) = cx(lon) * (1.0_r8 - delomg(lon) * delomg(lon))
             END IF
             IF (cx(lon) .GT. clow(lon)) THEN
                klow(lon) = kl
                clow(lon) = cx(lon)
             END IF
          ELSE IF (pll(lon,kl) .GT. f400p) THEN
             !
             !     3. define middle cloud ***
             !     middle cloud is defined one layer thick between 700 and 400 mb.
             !
             !     define middle clouds only in supersaturate type
             !
             IF (cx(lon) .GT. cmid(lon)) THEN
                kmid(lon) = kl
                cmid(lon) = cx(lon)
             END IF
          ELSE
             !
             !     4. define high cloud
             !     high cloud is defined only one layer thick from 400 mb and up.
             !
             !     define high clouds due to satuation
             !
             IF (cx(lon) .GT. chigh(lon)) THEN
                khigh(lon) = kl
                chigh(lon) = cx(lon)
             END IF
             !
             !     end of vertical computation
             !
          END IF
       END DO
    END DO
    DO lon = 1,ncols
       css(lon,kmax-khigh(lon)+1) = chigh(lon)
       cstc(lon,khigh(lon))=chigh(lon)
       css(lon,kmax-kmid(lon)+1) = cmid(lon)
       cstc(lon,kmid(lon))=cmid(lon)
       !
       !     for very thin low cloud adding its thickness
       !     pressure thickness of low cloud layer
       !
       IF(klow(lon).GE.lcnvl) THEN
          thklow=0.5_r8*(pll(lon,klow(lon)-1)-pll(lon,klow(lon)+1))
          IF(thklow.LE.REAL(lthncl,r8)) THEN
             css(lon,kmax-klow(lon)) = clow(lon)
             cstc(lon,klow(lon)+1)=clow(lon)
          ENDIF
          css(lon,kmax-klow(lon)+1) = clow(lon)
          cstc(lon,klow(lon))=clow(lon)
       ENDIF
    END DO
    !
    !     convective cloud is maximum overlaping in ir subr.
    !
    IF (covlp .EQ. 'MAXI') THEN
       DO lon = 1,ncols
          DO kl = kbot(lon), ktop(lon)
             ccu(lon,kmax-kl+1) = conv(lon)
          END DO
       END DO
       !
       !     convective cloud is random overlaping in ir subr.
       !
    ELSE IF (covlp .EQ. 'RAND') THEN
       DO lon = 1,ncols
          DO kl = kbot(lon), ktop(lon)
             css(lon,kmax-kl+1) = MIN(1.0_r8,css(lon,kmax-kl+1))
             IF(mxrdcc)THEN
                css(lon,kmax-kl+1) = MAX(conv(lon),css(lon,kmax-kl+1))
                ccon(lon,kl)=conv(lon)
             END IF
          END DO
       END DO
    END IF
  END SUBROUTINE cldgen


  !hmjb This is the same as cldgen(), but with some bugs corrected.
  !The main change is that now both convective and layer clouds
  !are outputed separately. A new subroutine was necessary
  !to keep the original results with kuo. However,  I think
  !kuo should be used with this correction as well... but for
  !that futher studies about this impact must be done.
  !   This version is intended for, and was tested with, the new
  !shortwave radiation codes: clirad and ukmet
  !
  SUBROUTINE cldgn3 (prsl   ,prsi,grh   ,omg   ,gtmp  ,css   , ccu   , &
       cdin  ,cstc  ,ccon  ,cson  ,lcnvl ,lthncl, convc ,convt ,convb ,&
       ncols ,kmax  ,nls      )

    !==========================================================================
    ! cldgen :perform clouds generation scheme based on j. slingo
    !         (1984 ecmwf workshop); the scheme generates 4 type of clouds
    !         of convective, high, middle and low clouds.
    !==========================================================================
    ! parameters and input variables:
    !        covlp = 'maxi'         maximum overlap of convective cloud
    !                               or thick low cloud used in ir subr.
    !              = 'rand'         random  overlap of convective cloud
    !                               or thick low cloud used in ir subr.
    !        date  =  julian day of model forecast date
    !        grh   =  relative humidity  (0-1)
    !        omg   =  vertical velocity  (cb/sec)
    !        prsl  =  pressure   (Pa)
    !        sigmid   =  sigma coordinate at middle of layer
    !        gtmp  =  layer temperature (k)
    !     output variables:
    !        css   =  ncols*kmax supersatuation cloud cover fraction
    !        ccu   =  ncols*kmax convective cloud cover fraction
    !---------------------------------------------------------------------
    ! values from subr-gwater
    !       convc  =  ncols convective cloud cover in 3 hr. avrage
    !       convt  =  ncols convective cloud top  (sigma layer)
    !       convb  =  ncols convective cloud base (sigma layer)
    !==========================================================================
    !
    !   ncols......Number of grid points on a gaussian latitude circle
    !   kmax......Number of grid points at vertical
    !   nls..... .Number of layers in the stratosphere.
    !   cdin
    !   cstc......cstc=clow change necessary in order to properly mark inv
    !             cloud height
    !   ccon......convc  =  ncols convective cloud cover in 3 hr. avrage
    !   cson
    !   cp........Specific heat of air (j/kg/k)
    !   gasr......Constant of dry air      (j/kg/k)
    !   mxrdcc....use maximum random converage for radiative conv. clouds
    !             constant logical mxrdcc = .true.
    !   lcnvl.....the lowest layer index where non-convective clouds can
    !             occur (ben says this should be 2 or more)
    !             constant lcnvl = 2
    !   lthncl....Minimum depth in mb of non-zero low level cloud
    !             consta lthncl=80
    !==========================================================================
    INTEGER,          INTENT(IN   ) :: ncols
    INTEGER,          INTENT(in ) :: kmax
    INTEGER,          INTENT(in ) :: nls
    REAL(KIND=r8),             INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),             INTENT(in   ) :: prsi  (ncols,kMax+1)
!    REAL(KIND=r8),             INTENT(in   ) :: gps (ncols)
!    REAL(KIND=r8),             INTENT(in   ) :: sigmid (kmax)
    REAL(KIND=r8),             INTENT(in   ) :: grh (ncols,kmax)
    REAL(KIND=r8),             INTENT(in   ) :: omg (ncols,kmax)
    REAL(KIND=r8),             INTENT(in   ) :: gtmp(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: css (ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: ccu (ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cdin(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cstc(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: ccon(ncols,kmax)
    REAL(KIND=r8),             INTENT(inout) :: cson(ncols,kmax)

    INTEGER      ,             INTENT(in ) :: lcnvl
    INTEGER      ,             INTENT(in ) :: lthncl
    REAL(KIND=r8),             INTENT(in ) :: convc(ncols)
    REAL(KIND=r8),             INTENT(in ) :: convt(ncols)
    REAL(KIND=r8),             INTENT(in ) :: convb(ncols)

    !
    !     dthdpm-----min (d(theta)/d(p))
    !     invb ------the base of inversion layer
    !     dthdpc-----criterion of (d(theta)/d(p))
    !
    REAL(KIND=r8),    PARAMETER :: dthdpc = -0.4e-1_r8
    REAL(KIND=r8)    :: pll   (ncols,kmax)
    REAL(KIND=r8)    :: pii   (ncols,kmax+1)
    REAL(KIND=r8)    :: conv  (ncols)
    REAL(KIND=r8)    :: clow  (ncols)
    REAL(KIND=r8)    :: cmid  (ncols)
    REAL(KIND=r8)    :: chigh (ncols)
    INTEGER :: ktop  (ncols)
    INTEGER :: kbot  (ncols)
    INTEGER :: klow  (ncols)
    INTEGER :: kmid  (ncols)
    INTEGER :: khigh (ncols)
    REAL(KIND=r8)    :: pt    (ncols)
    REAL(KIND=r8)    :: cx    (ncols)
    REAL(KIND=r8)    :: cinv  (ncols)
    REAL(KIND=r8)    :: delomg(ncols)
    REAL(KIND=r8)    :: csat  (ncols,kmax)
    REAL(KIND=r8)    :: dthdpm(ncols)
    INTEGER :: invb  (ncols)

    REAL(KIND=r8), PARAMETER :: f700p= 8.0e2_r8!7.5e2_r8! 7.0e2_r8
    REAL(KIND=r8), PARAMETER :: f400p= 4.0e2_r8
    REAL(KIND=r8), PARAMETER :: f6p67= 6.67e0_r8
    REAL(KIND=r8), PARAMETER :: f0p9 = 0.9e0_r8
    REAL(KIND=r8), PARAMETER :: f0p8 = 0.8e0_r8
    REAL(KIND=r8), PARAMETER :: f0p4 = 0.4e0_r8
    REAL(KIND=r8), PARAMETER :: f0p3 = 0.3e0_r8
    REAL(KIND=r8), PARAMETER :: f0p2 = 0.2e0_r8
    REAL(KIND=r8), PARAMETER :: f0p6 = 0.6e0_r8
    REAL(KIND=r8), PARAMETER :: f1e4 = 1.0e+4_r8
    REAL(KIND=r8), PARAMETER :: f5m5 = 5.0e-5_r8

    INTEGER :: i
    INTEGER :: k
    INTEGER :: kl
    INTEGER :: lon
    REAL(KIND=r8)    :: arcp
    REAL(KIND=r8)    :: dthdp
    REAL(KIND=r8)    :: thklow


    DO k = 1, kmax
       DO i = 1, ncols
          css(i,k) = 0.0_r8
          ccu(i,k) = 0.0_r8
          csat(i,k) = 0.0_r8
          cdin(i,k) = 0.0_r8
          cstc(i,k) = 0.0_r8
          ccon(i,k) = 0.0_r8
          cson(i,k) = 0.0_r8
       END DO
    END DO
    !
    !     the clouds generation scheme is based on j. slingo
    !     (1984 ecmwf workshop).  the scheme generates 4 type of clouds
    !     of convective, high, middle and low clouds.
    !
    DO kl = 1, kmax
       DO lon = 1,ncols
          pll(lon,kl) = prsl(lon,kl)/100.0_r8!gps(lon) * sigmid(kl)
       END DO
    END DO
    DO kl = 1, kmax+1
       DO lon = 1,ncols
          pii(lon,kl) = prsi(lon,kl)/100.0_r8!gps(lon) * sigmid(kl)
       END DO
    END DO
    !
    !     initialization
    !
    DO  lon = 1,ncols
       conv(lon)  = convc(lon)
       clow(lon)  = 0.0_r8
       cinv(lon)  = 0.0_r8
       cmid(lon)  = 0.0_r8
       chigh(lon) = 0.0_r8
       ktop(lon)  = INT( convt(lon)+0.5_r8 )
       IF (ktop(lon).LT.1.OR.ktop(lon).GT.kmax) ktop(lon)=1
       kbot(lon)  = INT(convb(lon)+0.5_r8)
       IF (kbot(lon).LT.1.OR.kbot(lon).GT.kmax) kbot(lon)=1
       klow(lon)  = 1
       kmid(lon)  = 1
       khigh(lon) = 1
       !
       !     1. define convective cloud
       !     done in subr-gwater
       !     cloud top and base are defined by kuo scheme: convt, convb
       !     cloud amount is calculated from precipitation rate : convc
       !     *** single layer clouds conputations start here, from bottom up
       !
       !     define high clouds due to strong convection
       !
       pt(lon) =pll(lon,ktop(lon)) ! gps(lon) * sigmid(ktop(lon))
       IF ((pt(lon).LE.f400p).AND.&
            (conv(lon).GE.f0p4)) THEN
          chigh(lon) = 2.0_r8 * (conv(lon) - f0p3)
          chigh(lon) = MIN(chigh(lon),1.0_r8)
          khigh(lon) = MIN(ktop(lon) + 1,kmax)
          ccon(lon,khigh(lon))=chigh(lon)
       END IF
    END DO
    !
    !     compute low stratus associated with inversions, based on ecwmf's
    !     scheme, with lower criterion of d(theta)/d(p)
    !
    arcp = gasr / cp

    DO lon=1,ncols
       invb(lon) =  MIN(kmax, kmax-nls)
       dthdpm(lon) = 0.0_r8
    END DO

    DO kl=2,kmax
       DO lon=1,ncols
          IF (pll(lon,kl) .GT. f700p) THEN
             dthdp = (gtmp(lon,kl-1)*(1000.0_r8/pll(lon,kl-1))**arcp &
                     -gtmp(lon,kl  )*(1000.0_r8/pll(lon,kl  ))**arcp)/&
                  (pll(lon,kl-1) - pll(lon,kl))
             IF(dthdp.LT.0.0_r8) THEN
                invb(lon) = MIN(kl-1,invb(lon))
                IF(dthdp.LT.dthdpc) THEN
                   IF(dthdp.LT.dthdpm(lon)) THEN
                      dthdpm(lon)=dthdp
                      klow(lon)=kl
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       END DO
    END DO
    !
    !     klow change above necessary to mark inversion cloud height
    !
    DO lon=1,ncols
       IF(dthdpm(lon).LT.dthdpc.AND.grh(lon,invb(lon)).GT.f0p6)THEN
          cinv(lon) = - f6p67 * (dthdpm(lon)-dthdpc)
          cinv(lon) = MAX(cinv(lon),0.0_r8)
          cinv(lon) = MIN(cinv(lon),1.0_r8)
          IF (grh(lon,invb(lon)) .LT. f0p8)&
               cinv(lon) = cinv(lon)*&
               (1.0_r8-(f0p8 - grh(lon,invb(lon)))/f0p2)
       ENDIF
       clow(lon)=cinv(lon)
       cdin(lon,klow(lon))=cinv(lon)
    END DO
    !
    !     clow change necessary in order to properly mark inv cloud height
    !
    !     main loop for cloud amount determination
    !
    DO kl = lcnvl+1, MIN(kmax, kmax-nls)
       DO lon = 1,ncols
          !
          !     general define cloud due to saturation
          !
          IF (pll(lon,kl) .GT. f400p) THEN
             cx(lon) = (grh(lon,kl) - f0p8) / f0p2
          ELSE
             cx(lon) = (grh(lon,kl) - f0p9) / 0.1_r8
          END IF
          cx(lon) = (MAX(cx(lon), 0.0_r8)) ** 2
          cx(lon) =  MIN(cx(lon), 1.0_r8)
          !
          !     start vertical process from bottom to top
          !
          IF (pll(lon,kl) .GT. f700p) THEN
             !
             !     2. define low cloud ***
             !     low cloud is defined one layer thick ranging from layer 3 to 700mb
             !     there are two type possible generating mechanisms. due boundary
             !     t inversion type and associated with vertical motion.
             !
             !     define low super satuated clouds but adjusted by vertical motion
             !
             IF (omg(lon,kl) .GE. f5m5) THEN
                cx(lon) = 0.0_r8
             ELSE IF (omg(lon,kl) .GE. -f5m5) THEN
                delomg(lon) = (omg(lon,kl)+f5m5)*f1e4
                cx(lon) = cx(lon) * (1.0_r8 - delomg(lon) * delomg(lon))
             END IF
             IF (cx(lon) .GT. clow(lon)) THEN
                klow(lon) = kl
                clow(lon) = cx(lon)
             END IF
          ELSE IF (pll(lon,kl) .GT. f400p) THEN
             !
             !     3. define middle cloud ***
             !     middle cloud is defined one layer thick between 700 and 400 mb.
             !
             !     define middle clouds only in supersaturate type
             !
             IF (cx(lon) .GT. cmid(lon)) THEN
                kmid(lon) = kl
                cmid(lon) = cx(lon)
             END IF
          ELSE
             !
             !     4. define high cloud
             !     high cloud is defined only one layer thick from 400 mb and up.
             !
             !     define high clouds due to satuation
             !
             IF (cx(lon) .GT. chigh(lon)) THEN
                khigh(lon) = kl
                chigh(lon) = cx(lon)
             END IF
             !
             !     end of vertical computation
             !
          END IF
       END DO
    END DO
    DO lon = 1,ncols
       css (lon,kmax-khigh(lon)+1) = chigh(lon)
       cstc(lon,khigh(lon))=chigh(lon)
       css (lon,kmax-kmid(lon)+1) = cmid(lon)
       cstc(lon,kmid(lon))=cmid(lon)
       !
       !     for very thin low cloud adding its thickness
       !     pressure thickness of low cloud layer
       !
       IF(klow(lon).GE.lcnvl) THEN
          thklow=0.5_r8*(pll(lon,klow(lon)-1)-pll(lon,klow(lon)+1))
          IF(thklow.LE.REAL(lthncl,r8)) THEN
             css(lon,kmax-klow(lon)) = clow(lon)
             cstc(lon,klow(lon)+1)=clow(lon)
          ENDIF
          css(lon,kmax-klow(lon)+1) = clow(lon)
          cstc(lon,klow(lon))=clow(lon)
       ENDIF
    END DO
    !hmjb - Correction of many bugs
    !          The way it was before, convective went out as 0.0 and
    !       strat had the maximum of convc and strat. Now this is fixed
    !
    ! Bugs.
    ! line 04 Should always be executed. shortwave subroutine need convective clouds
    ! line 13 Should always be executed as it is a numerical check only
    ! line 15 This is done inside shortwave subroutines, so we don't need to do it here.
    !         However, the values modified there will not go inside longwave subroutine,
    !         and they should... must fix this in spmrad
    ! line 16 ccon is only used for output... So it should always receive conv
    !
    !OLD CODE-------------------------------------------------------------
    !    !
    !    !     convective cloud is maximum overlaping in ir subr.
    !    !
    !01    IF (covlp .EQ. 'MAXI') THEN
    !02       DO lon = 1,ncols
    !03          DO kl = kbot(lon), ktop(lon)
    !04             ccu(lon,kmax-kl+1) = conv(lon)
    !05          END DO
    !06       END DO
    !07       !
    !08       !     convective cloud is random overlaping in ir subr.
    !09       !
    !10    ELSE IF (covlp .EQ. 'RAND') THEN
    !11       DO lon = 1,ncols
    !12          DO kl = kbot(lon), ktop(lon)
    !13             css(lon,kmax-kl+1) = MIN(1.0_r8,css(lon,kmax-kl+1))
    !14             IF(mxrdcc)THEN
    !15                css(lon,kmax-kl+1) = MAX(conv(lon),css(lon,kmax-kl+1))
    !16                ccon(lon,kl)=conv(lon)
    !17             END IF
    !18          END DO
    !19       END DO
    !20    END IF
    !NEW CODE-------------------------------------------------------------
    DO lon = 1,ncols
       DO kl = kbot(lon), ktop(lon)
          ccu(lon,kmax-kl+1) = conv(lon)
          css(lon,kmax-kl+1) = MIN(1.0_r8,css(lon,kmax-kl+1))
          ccon(lon,kl)=conv(lon)
       END DO
    END DO

    RETURN
  END SUBROUTINE cldgn3


  ! getoz  :interpolates climatological ozone data into a given latitude
  !         and model julian date.
  SUBROUTINE getoz (ncols,kmax,prsl,prsi,colrad,date,o3l)
    !
    ! input parameters and variables:
    !     ncols  =  number of atmospheric columns
    !     kmax   =  number of atmospheric layers
    !     sigmid =  sigma coordinate at middle of layer
    !     colrad =  colatitude of each column (0-3.14 from np to sp in radians)
    !     date   =  model julian date
    !
    ! tabulated data
    !     ozone  =  climatological ozone mixing ratio in 18 sigma layers 
    !               and in 5 degree latitude interval
    !
    ! output variables:
    !     o3l   =  18 layers ozone mixing ratio in given lat and date
    !
    !==========================================================================
    ! :: kmax.....Number of grid points at vertical
    ! :: sigmid.......sigma coordinate at middle of layer
    ! :: pai......constant pi=3.1415926
    ! :: yrl......length of year in days
    !==========================================================================
    
    ! Input variables
    INTEGER      ,    INTENT(IN   ) :: ncols
    INTEGER      ,    INTENT(IN   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: prsi  (ncols,kMax+1)
    !REAL(KIND=r8),    INTENT(IN   ) :: sigmid (kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: colrad(ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: date

    ! Output variable
    REAL(KIND=r8),    INTENT(INOUT) :: o3l(ncols,kmax)

    ! Local variables
    REAL(KIND=r8) :: a1   (nlm_getoz)
    REAL(KIND=r8) :: a2   (nlm_getoz)
    REAL(KIND=r8) :: a3   (nlm_getoz)
    REAL(KIND=r8) :: a4   (nlm_getoz)
    REAL(KIND=r8) :: b1   (nlm_getoz)
    REAL(KIND=r8) :: b2   (nlm_getoz)
    REAL(KIND=r8) :: b3   (nlm_getoz)
    REAL(KIND=r8) :: b4   (nlm_getoz)
    REAL(KIND=r8) :: do3a (nlm_getoz)
    REAL(KIND=r8) :: do3b (nlm_getoz)
    REAL(KIND=r8) :: ozo3l(ncols,nlm_getoz)

    REAL(KIND=r8), PARAMETER :: rlag = 14.8125e0_r8

    INTEGER :: l
    INTEGER :: la   (ncols)
    INTEGER :: ll   (ncols)
!    INTEGER :: kmx
    INTEGER :: imon
    INTEGER :: isea
    INTEGER :: nsur

    INTEGER :: k
    INTEGER :: i
    INTEGER :: kk    (ncols,kmax)
    REAL(KIND=r8)    :: theta
    REAL(KIND=r8)    :: flat
    REAL(KIND=r8)    :: rang
    REAL(KIND=r8)    :: rsin1(ncols)
    REAL(KIND=r8)    :: rcos1(ncols)
    REAL(KIND=r8)    :: rcos2(ncols)
    REAL(KIND=r8)    :: rate (ncols)
    REAL(KIND=r8)    :: aa
    REAL(KIND=r8)    :: bb
    LOGICAL :: notfound(ncols,kMax)
    
    nsur=1
    
    !kmx=nlm_getoz
    !
    !     find closest place in the data according to input slat.
    !
    IF(date.GT.year_getoz) date=date-year_getoz

    imon=INT(date/mon_getoz) + 1

    IF(imon.LT.1)imon=1

    isea=imon/3 + 1

    IF(isea.EQ.5) isea=1
    IF(isea.GT.5) THEN
       WRITE(nfprt,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       WRITE(nferr,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       STOP 9954
    END IF
    DO i=1,ncols
       theta = 90.0_r8-(180.0_r8/pai)*colrad(i) ! colatitude -> latitude
       ! the 180 degrees are divided into 37 bands with 5deg each
       ! except for the first and last, which have 2.5 deg
       ! The centers of the bands are located at:
       !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
       flat  = 0.2_r8*theta ! indexing the latitudes: goes from -18. to +18.
       ! find the latitude index before and after each latitude
       la(i)    = INT(19.501e0_r8-flat )!
       ll(i)    = INT(19.001e0_r8-flat)

       !
       !     find sin and cos coefficients for time interpolation.
       !
       rang=2.0e0_r8*pai*(date-rlag)/year_getoz
       rsin1(i)=SIN(rang)
       rcos1(i)=COS(rang)
       rcos2(i)=COS(2.0e0_r8*rang)
       rate(i)=REAL(19-ll(i),r8)-flat
       !
       !     ozone interpolation in latitude and time
       !
    END DO
    DO k=1,nlm_getoz
       DO i=1,ncols
          a1(k) =2.5e-1_r8*(ozone(k,la(i),1)+ozone(k,la(i),2)+ &
               ozone(k,la(i),3)+ozone(k,la(i),4))
          a2(k) =0.5e0_r8*(ozone(k,la(i),2)-ozone(k,la(i),4))
          a3(k) =0.5e0_r8*(ozone(k,la(i),1)-ozone(k,la(i),3))
          a4(k) =2.5e-1_r8*(ozone(k,la(i),1)+ozone(k,la(i),3)- &
               ozone(k,la(i),2)-ozone(k,la(i),4))
          b1(k) =2.5e-1_r8*(ozone(k,ll(i),1)+ozone(k,ll(i),2)+ &
               ozone(k,ll(i),3)+ozone(k,ll(i),4))
          b2(k) =0.5e0_r8*(ozone(k,ll(i),2)-ozone(k,ll(i),4))
          b3(k) =0.5e0_r8*(ozone(k,ll(i),1)-ozone(k,ll(i),3))
          b4(k) =2.5e-1_r8*(ozone(k,ll(i),1)+ozone(k,ll(i),3)- &
               ozone(k,ll(i),2)-ozone(k,ll(i),4))
          do3a(k)=a1(k)+rsin1(i)*a2(k)+rcos1(i)*a3(k)+rcos2(i)*a4(k)
          do3b(k)=b1(k)+rsin1(i)*b2(k)+rcos1(i)*b3(k)+rcos2(i)*b4(k)
          ozo3l(i,k)=do3a(k)+rate(i)*(do3b(k)-do3a(k))
          ozo3l(i,k)=1.0e-04_r8*ozo3l(i,k)
       END DO
    END DO
    !ozsig(:) = (/ &
    !     .020747_r8,.073986_r8,.124402_r8,.174576_r8,.224668_r8,.274735_r8, &
    !     .324767_r8,.374806_r8,.424818_r8,.497450_r8,.593540_r8,.688125_r8, &
    !     .777224_r8,.856317_r8,.920400_r8,.960480_r8,.981488_r8,.995004_r8/)

    IF(inter_getoz)THEN
       DO l=1,kmax
          DO i=1,nCols
            !notfound(i,l) = sigmid(l) > ozsig(1)
             notfound(i,l) = prsl(i,l)/prsl(i,nsur) > ozsig(1)
             IF (notfound(i,l)) THEN
                kk(i,l)=nlm_getoz
             ELSE
                kk(i,l)=2
             END IF
          END DO
       END DO
       DO l=1,kmax
          DO i=1,nCols
             IF (notfound(i,l)) THEN
                DO k=2,nlm_getoz
                  !IF(sigmid(l).GT.ozsig(k-1).AND.sigmid(l).LE.ozsig(k))THEN
                   IF( prsl(i,l)/prsl(i,nsur) > ozsig(k-1).AND. prsl(i,l)/prsl(i,nsur) <= ozsig(k))THEN

                      kk(i,l)=k
                      EXIT
                   END IF
                END DO
             END IF
          END DO
       END DO
       !DO l = 1, kmax
       !   DO i= 1, ncols
       !      kk(i,l) = kk(i,l)
       !   END DO
       !END DO
    END IF
    IF(inter_getoz)THEN
       DO l=1,kmax
          DO i=1,ncols
             aa=(ozo3l(i,kk(i,l))-ozo3l(i,kk(i,l)-1))/(ozsig(kk(i,l))-ozsig(kk(i,l)-1))
             bb= ozo3l(i,kk(i,l)-1)-aa*ozsig(kk(i,l)-1)
             o3l(i,kmax+1-l)=bb+aa*prsl(i,l)/prsl(i,nsur)
          END DO
       END DO
    END IF
    IF(.NOT.inter_getoz)THEN
       DO l=1,nlm_getoz
          DO i=1,ncols
             o3l(i,l)=ozo3l(i,l)
          END DO
       END DO
    ENDIF
  END SUBROUTINE getoz

  SUBROUTINE Cloud_Field(&
       ! Model info and flags
       ncols     , kmax     , nls      , icld     , covlp    , mxrdcc   , &
       lcnvl     , lthncl   , colrad   , imask    ,                       &
       ! Atmospheric Fields
       prsi   ,prsl   ,phii   ,phil    ,&
       Qrel     , omg      , Te       ,            &
       ! Convective Clouds
       convc     , convt    , convb    , convts   , convcs   , convbs   , &
       CLDF      ,&       
       ! Flipped Clouds for Radiation
       cld       , clu      ,                                             &
       ! Clouds for Diagnostics
       CloudCover, CldCovTot, CldCovInv, CldCovSat, CldCovCon, CldCovSha, &
       CldLow    ,CldMed    ,CldHgh  )

    IMPLICIT NONE

    !--------------------------------------------------------------------------------- 
    ! Input
    !--------------------------------------------------------------------------------- 

    ! Model info and flags
    INTEGER, INTENT(IN) :: ncols
    INTEGER, INTENT(IN) :: kmax
    INTEGER, INTENT(IN) :: nls
    INTEGER, INTENT(IN) :: icld
    CHARACTER(len=4), INTENT(IN) :: covlp
    LOGICAL, INTENT(IN) :: mxrdcc                     
    INTEGER, INTENT(IN) :: lcnvl                      
    INTEGER, INTENT(IN) :: lthncl                     
    REAL(KIND=r8),    INTENT(in   ) :: colrad(ncols)      
    INTEGER(KIND=i8), INTENT(IN   ) :: imask (ncols)

    ! Atmospheric Fields
    REAL(KIND=r8),    INTENT(in   ) :: prsi  (ncols,kMax+1)
    REAL(KIND=r8),    INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: phii  (nCols,kMax+1)
    REAL(KIND=r8),    INTENT(in   ) :: phil  (nCols,kMax)
!    REAL(KIND=r8), INTENT(IN) :: sig  (kmax)
!    REAL(KIND=r8), INTENT(IN) :: Ps  (ncols)  
    REAL(KIND=r8), INTENT(IN) :: Qrel  (ncols,kmax)  
    REAL(KIND=r8), INTENT(IN) :: omg  (ncols,kmax) 
    REAL(KIND=r8), INTENT(IN) :: Te (ncols,kmax)

    ! Convective Clouds
    REAL(KIND=r8), INTENT(IN) :: convc (ncols)        
    REAL(KIND=r8), INTENT(IN) :: convt (ncols)        
    REAL(KIND=r8), INTENT(IN) :: convb (ncols)        
    REAL(KIND=r8), INTENT(IN) :: convts(ncols)        
    REAL(KIND=r8), INTENT(IN) :: convcs(ncols)        
    REAL(KIND=r8), INTENT(IN) :: convbs(ncols)        
    REAL(KIND=r8), INTENT(IN) :: CLDF  (nCols,kMax)
    !--------------------------------------------------------------------------------- 
    ! Output
    !--------------------------------------------------------------------------------- 

    ! Flipped Clouds for Radiation
    REAL(KIND=r8), INTENT(INOUT) :: cld (ncols,kmax) ! = max(Inv,Sat,Sha)
    REAL(KIND=r8), INTENT(INOUT) :: clu (ncols,kmax) ! = Con

    ! Clouds Diagnostics
    REAL(KIND=r8), INTENT(INOUT) :: CloudCover(ncols)
    REAL(KIND=r8), INTENT(INOUT) :: CldCovTot (ncols,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: CldCovInv (ncols,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: CldCovSat (ncols,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: CldCovCon (ncols,kmax)
    REAL(KIND=r8), INTENT(INOUT) :: CldCovSha (ncols,kmax)
    REAL(KIND=r8), INTENT(inout) :: CldLow (ncols) 
    REAL(KIND=r8), INTENT(inout) :: CldMed (ncols) 
    REAL(KIND=r8), INTENT(inout) :: CldHgh (ncols) 



    !--------------------------------------------------------------------------------- 
    ! Local Variables
    !--------------------------------------------------------------------------------- 
    !  ---  pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
    REAL(kind=r8), PARAMETER :: con_pi     =3.1415926535897931_r8 ! pi
    INTEGER, PARAMETER ::iflip= 1!      iflip           : control flag for in/out vertical indexing      !
    !   iflip           : control flag for in/out vertical indexing         !
    !                     =0: index from toa to surface                     !
    !                     =1: index from surface to toa                     !
    !   iovr            : control flag for cloud overlap                    !
    !                     =0 random overlapping clouds                      !
    !                     =1 max/ran overlapping clouds                     !
    INTEGER, PARAMETER :: iovr =0 !      iovrsw/iovrlw   : control flag for cloud overlap (sw/lw rad)     !
    !                        =0 random overlapping clouds                   !
    !                        =1 max/ran overlapping clouds                  !
    REAL (kind=r8) :: clds   (ncols,5)  !real   , intent(out) :: clds(:,:)fraction of clouds for low, mid, hi, tot, bl   
    INTEGER  :: mtop   (ncols,3)  !integer, intent(out) :: mtop(:,:)vertical indices for low, mid, hi cloud tops        
    INTEGER  :: mbot   (ncols,3)  !integer, intent(out) :: mbot(:,:)vertical indices for low, mid, hi cloud bases

    REAL (kind=r8), PARAMETER :: ptopc(1:4,1:2)=  RESHAPE( (/ 1050.0_r8, 650.0_r8, 400.0_r8, 0.0_r8,  &
       1050.0_r8, 750.0_r8, 500.0_r8, 0.0_r8/),(/4,2/) )
    REAL (kind=r8) :: tem1
    !-- Aux variables
    REAL (kind=r8) :: ptop1(ncols,4),xlat(ncols)
    REAL (kind=r8) :: plyr (ncols,kmax),cldtot (ncols,kmax),cldcnv(ncols,kmax)
    INTEGER :: i,k,kflip,id

    !--------------------------------------------------------------------------------- 
    !--------------------------------------------------------------------------------- 

    cld    = 0.0_r8
    clu    = 0.0_r8
    CldCovTot = 0.0_r8
    CldCovInv = 0.0_r8
    CldCovSat = 0.0_r8
    CldCovCon = 0.0_r8
    CldCovSha = 0.0_r8
    CloudCover = 0.0_r8

    IF (TRIM(iswrad).eq.'LCH'.OR.TRIM(iswrad).eq.'NON') THEN
       IF(icld.EQ.1) &  
            CALL cldgen (covlp ,prsl,prsi,Qrel   ,omg   ,Te  ,cld  , &
            clu  ,CldCovInv  ,CldCovSat  ,CldCovCon  ,CldCovSha  ,mxrdcc,lcnvl ,lthncl, &
            convc ,convt ,convb ,ncols ,kmax  ,nls     )
       
       IF(icld.EQ.3) &
            CALL cldgn2 ( covlp ,prsl,prsi,Qrel   ,omg   ,Te  ,cld  , clu  ,&
            CldCovInv  ,CldCovSat  ,CldCovCon  ,CldCovSha  ,mxrdcc,lcnvl ,convc ,convt ,convb , &
            convts,convcs,convbs,ncols ,kmax  ,nls            )
    ELSE
       CALL cldgn3 (prsl   ,prsi,Qrel   ,omg   ,Te  ,cld  , &
            clu  ,CldCovInv  ,CldCovSat  ,CldCovCon  ,CldCovSha  ,lcnvl ,lthncl, &
            convc ,convt ,convb ,ncols ,kmax  ,nls     )
    ENDIF
    !SUBROUTINE cldgn3 (prsl   ,prsi,grh   ,omg   ,gtmp  ,css   , ccu   , &
    !      cld -> css   =  ncols*kmax supersatuation cloud cover fraction
    !      clu -> ccu   =  ncols*kmax convective cloud cover fraction
    IF(TRIM(ILCON).EQ.'LSC' .OR. TRIM(ILCON).EQ.'YES' ) THEN
      ! CLDF=0.0_r8
    ELSE
    IF((TRIM(ISWRAD) == 'RRTMG'.and. TRIM(ILWRAD) == 'RRTMG')) THEN
       DO k=1,kmax
          kflip=kmax+1-k
          DO i=1,ncols
             IF(CldCovSat(i,k) > 0.01_r8)THEN
                IF (imask(i) .lt. 1_i8) THEN
                   !ocean  
                 !  cld(i,kflip) = MAX(CLDF  (i,k) , cld(i,kflip))          !supersatuation cloud cover fraction
                   cld(i,kflip) = (CLDF  (i,k) + cld(i,kflip))/2.0_r8      !supersatuation cloud cover fraction
                ELSE
                   !land
                   cld(i,kflip) = (CLDF  (i,k) + cld(i,kflip))/2.0_r8      !supersatuation cloud cover fraction
                END IF
             ELSE
                cld(i,kflip) = 0.0_r8
             END IF  
                IF (imask(i) .lt. 1_i8) THEN
                   !ocean  
                  ! clu(i,kflip)    = MAX(CLDF  (i,k) , clu(i,kflip))        !convective cloud cover fraction
                   clu(i,kflip)    = (CLDF  (i,k) + clu(i,kflip))/2.0_r8        !convective cloud cover fraction

                ELSE
                   !land
                   clu(i,kflip)    = (CLDF  (i,k) + clu(i,kflip))/2.0_r8        !convective cloud cover fraction
                END IF
          END DO
       END DO
    END IF
    IF((TRIM(ISWRAD) == 'CRD'  .and. TRIM(ILWRAD) == 'CRD') ) THEN
       DO k=1,kmax
          kflip=kmax+1-k
          DO i=1,ncols
             IF(CldCovSat(i,k) > 0.01_r8)THEN
                IF (imask(i) .lt. 1_i8) THEN
                   !ocean  
!                   cld(i,kflip) = MAX(CLDF  (i,k) , cld(i,kflip))          !supersatuation cloud cover fraction
                   cld(i,kflip) = (CLDF  (i,k) + cld(i,kflip))/2.0_r8      !supersatuation cloud cover fraction
                ELSE
                   !land
                   cld(i,kflip) = (CLDF  (i,k) + cld(i,kflip))/2.0_r8      !supersatuation cloud cover fraction
                END IF
             ELSE
                cld(i,kflip) = 0.0_r8
             END IF  
                IF (imask(i) .lt. 1_i8) THEN
                   !ocean  
                   !clu(i,kflip)    = MAX(CLDF  (i,k) , clu(i,kflip))        !convective cloud cover fraction
                   clu(i,kflip)    = (CLDF  (i,k) + clu(i,kflip))/2.0_r8        !convective cloud cover fraction

                ELSE
                   !land
                   clu(i,kflip)    = (CLDF  (i,k) + clu(i,kflip))/2.0_r8        !convective cloud cover fraction
                END IF
          END DO
       END DO
    END IF
    IF((TRIM(ISWRAD) == 'CRDTF'  .and. TRIM(ILWRAD) == 'CRDTF')) THEN
       DO k=1,kmax
          kflip=kmax+1-k
          DO i=1,ncols
             IF(CldCovSat(i,k) > 0.01_r8)THEN
                IF (imask(i) .lt. 1_i8) THEN
                   !ocean  
!                   cld(i,kflip) = MAX(CLDF  (i,k) , cld(i,kflip))          !supersatuation cloud cover fraction
                   cld(i,kflip) = (CLDF  (i,k) + cld(i,kflip))/2.0_r8      !supersatuation cloud cover fraction

                ELSE
                   !land
                   cld(i,kflip) = (CLDF  (i,k) + cld(i,kflip))/2.0_r8      !supersatuation cloud cover fraction
                END IF
             ELSE
                cld(i,kflip) = 0.0_r8
             END IF  
                IF (imask(i) .lt. 1_i8) THEN
                   !ocean  
!                   clu(i,kflip)    = MAX(CLDF  (i,k) , clu(i,kflip))        !convective cloud cover fraction
                   clu(i,kflip)    = (CLDF  (i,k) + clu(i,kflip))/2.0_r8        !convective cloud cover fraction

                ELSE
                   !land
                   clu(i,kflip)    = (CLDF  (i,k) + clu(i,kflip))/2.0_r8        !convective cloud cover fraction
                END IF
          END DO
       END DO
    END IF
    END IF
    !     
    !     compute cloud cover diagnostic
    !     
    DO i=1,ncols
       CloudCover(i) = 1.0e0_r8
    END DO

    if (TRIM(iswrad).eq.'LCH'.OR.TRIM(iswrad).eq.'NON') then
       DO k=1,kmax
          kflip=kmax+1-k
          DO i=1,ncols
             CldCovTot(i,k)=cld(i,kflip) !->CldCovTot (total cloud cover)
             CloudCover(i) = CloudCover(i)*(1.0e0_r8-cld(i,k))
          END DO
       END DO
    else
       DO k=1,kmax
          kflip=kmax+1-k
          DO i=1,ncols
             CldCovTot(i,k)=max(cld(i,kflip),clu(i,kflip)) !->CldCovTot (total cloud cover)
             CloudCover(i) = CloudCover(i)*(1.0e0_r8-CldCovTot(i,k))
          END DO
       END DO
    endif

    DO i=1,ncols
       CloudCover(i) = 1.0e0_r8-CloudCover(i)
    END DO


    DO i = 1, ncols
       xlat(i)=(((colrad(i)))-(3.1415926e0_r8/2.0_r8))
    ENDDO

    !  ---  find top pressure for each cloud domain for given latitude
    !       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
    !  ---  i=1,2 are low-lat (<45 degree) and pole regions)

    DO id = 1, 4
       tem1 = ptopc(id,2) - ptopc(id,1)

       DO i =1, ncols
          ptop1(i,id) = ptopc(id,1) +                                   &
               tem1 * MAX( 0.0_r8, 4.0_r8*ABS(xlat(i))/con_pi-1.0_r8 )
       ENDDO
    ENDDO
    DO k=1,kmax
       kflip=kmax+1-k
       DO i=1,ncols
          plyr(i,k) = 10.0_r8 * prsl(i,k)/100.0_r8 !Pa -- > mb
          cldtot(i,k) =CldCovTot(i,k)
          cldcnv(i,k) =clu(i,kflip)
       ENDDO
    ENDDO

!          pll(lon,kl) = prsl(lon,kl)/100.0_r8!gps(lon) * sigmid(kl)

    !  ---  compute low, mid, high, total, and boundary layer cloud fractions
    !       and clouds top/bottom layer indices for low, mid, and high clouds.
    !       The three cloud domain boundaries are defined by ptopc.  The cloud
    !       overlapping method is defined by control flag 'iovr', which is
    !  ---  also used by the lw and sw radiation programs.

    CALL gethml ( &
                                !  ---  inputs:
         plyr   (1:ncols,1:kmax) , &!real   , intent(in) :: plyr(:,:) model layer mean pressure in mb (100Pa)
         ptop1  (1:ncols,1:4)    , &!real   , intent(in) :: ptop1(:,:)pressure limits of cloud domain interfaces  in mb (100Pa)
         cldtot (1:ncols,1:kmax) , &!real   , intent(in) :: cldtot(:,:)total or straiform cloud profile in fraction
         cldcnv (1:ncols,1:kmax) , &!real   , intent(in) :: cldcnv(:,:)convective cloud (for diagnostic scheme only)
         ncols                   , &!integer, intent(in) :: IX      horizontal dimention               
         kmax                    , &!integer, intent(in) :: NLAY    vertical layer dimensions            
         iflip                   , &!integer, intent(in) :: iflip  control flag for in/out vertical indexing
         iovr                    , &!integer, intent(in) :: iovr  control flag for cloud overlap          
                                !  ---  outputs:
         clds   (1:ncols,1:5)    , &!real   , intent(out) :: clds(:,:)fraction of clouds for low, mid, hi, tot, bl   
         mtop   (1:ncols,1:3)    , &!integer, intent(out) :: mtop(:,:)vertical indices for low, mid, hi cloud tops        
         mbot   (1:ncols,1:3)      )!integer, intent(out) :: mbot(:,:)vertical indices for low, mid, hi cloud bases

       DO i =1, ncols
          CldLow(i)=clds(i,1)
          CldMed(i)=clds(i,2)
          CldHgh(i)=clds(i,3)
       END DO

    DO k=1,kmax
       DO i =1, ncols
          IF(prsl(i,k)/100.0_r8 > 850.0_r8)THEN
             CldLow (i) = MAX(CldLow(i), cldtot (i,k)) 
          ELSE IF(prsl(i,k)/100.0_r8 <= 850.0_r8 .and. prsl(i,k)/100.0_r8 >= 400.0_r8 )THEN
             CldMed (i) = MAX(CldMed(i), cldtot (i,k)) 
          ELSE
             CldHgh (i) = MAX(CldHgh(i), cldtot (i,k)) 
          END IF
       END DO
    ENDDO

  END SUBROUTINE Cloud_Field

  SUBROUTINE Cloud_Micro_CAM5(&
       ! Model info
       ncols, kmax , imask   , &
       ! Atmospheric Fields
       prsi   ,prsl   ,phii   ,phil    ,&
       Te   , Qe    , tsea  , FlipPbot ,        &
       QCF  ,QCL   ,                                   &
       ! Cloud properties
       clwp , lmixr, fice  , rei   , rel   , taud       )
    IMPLICIT NONE

    ! As in the CCM2, cloud optical properties in the CCM3 are accounted for using
    ! the Slingo (1989) parameterization for liquid water droplet clouds. This
    ! scheme relates the extinction optical depth, the single-scattering albedo,
    ! and the asymmetry parameter to the cloud liquid water path and cloud drop
    ! effective radius. The latter two microphysical cloud properties were
    ! statically specified in the CCM2. In particular, in-cloud liquid water paths
    ! were evaluated from a prescribed, meridionally and height varying, but
    ! time independent, cloud liquid water density profile, rho_l(z), which
    ! was analytically determined on the basis of a meridionally specified
    ! liquid water scale height (e.g. see Kiehl et al., 1994; Kiehl, 1991).
    ! The cloud drop effective radius was simplly specified to be 10microns
    ! for all clouds. The CCM3 continues to diagnose cloud optical properties,
    ! but relaxes the rigid CCM2 framework. CCM3 employs the same exponentially
    ! decaying vertical profile for in-cloud water concentration
    !
    !             rho_l=rho_l^0*exp(-z/h_l)               eq 4.a.11
    !
    ! , where rho_l^0=0.21g/m3. Instead of specifying a zonally symmetric meridional
    ! dependence for the cloud water scale heigh, h_l, it is locally diagnosed
    ! as a function of the vertically integrated water vapor (precipitable water) 
    !
    !          h_l=700 ln [1+\frac{1}{g} \int_pT^ps q dp]  eq 4.a.12
    !
    ! hmjb> It is not explained, but the units of h_l must be meters, the same 
    ! hmjb> of the height, z.
    !
    ! The cloud water path (CWP) is determined by integrating the liquid
    ! water concentration using
    !
    !                 cwp = int rho_l dz     eq. 4.a.13
    ! 
    ! Which can be analytically evaluated for an arbitrary layer k as
    !
    !  rho_l^0 h_l [exp(-z_bot(k)/h_l) - exp(-z_top(k)/h_l)]   eq. 4.a.14
    !
    ! Where z_bot and z_top are the heights of the k'th layer interfaces.
    !
    ! hmjb> It is not explained, but the units of clwp must be g/m2
    ! hmjb> since it is the integral of rho_l*dz (eq.4.a.13)
    !
    ! CCM3 Documentation, pg 50
    ! Observational studies have shown a distinct difference between
    ! maritime and continental effective cloud drop size, r_e, for warm
    ! clouds. For this reason, the CCM3 differentiates between the cloud
    ! drop effective radius for clouds diagnosed over maritime and
    ! continental regimes (Kiehl, 1994). Over the ocean, the cloud drop
    ! effective radius for liquid water clouds, r_el, is specified to be
    ! 10microns, as in the CCM3. Over land masses r_el is determinedusing
    !
    ! r_el = 5 microns             T > -10oC
    !      = 5-5(t+10)/20 microns  -30oC <= T <= -10oC     eq. 4.a.14.1
    !      = r_ei                  T < -30oC
    !
    ! An ice particle effective radius, r_ei, is also diagnosed by CCM3,
    ! which at the moment amounts to a specification of ice radius as a
    ! function of normalized pressure
    !
    ! r_ei = 10 microns                                 p/ps > p_I^high   
    !      = r_ei^max - (r_ei^max - r_ei^min)           p/ps <= p_I^high       eq. 4.a.15.1
    !            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    ! where r_ei^max=30microns, r_ei^min=10microns, p_I^high=0.4 and p_I^low=0.0
    !
    ! hmjb>> I think there is a typo in the equation, otherwise the 
    ! hmjb>> expression for r_ei is not a continuous funcion of p/ps.
    ! hmjb>> For p/ps=p_I^high, r_ei should be r_ei^min and not r_ei^max.
    ! hmjb>> The correct equation is:
    ! hmjb>> r_ei = 10 microns                                 p/ps > p_I^high   
    ! hmjb>>      = R_EI^MIN - (r_ei^max - r_ei^min)           p/ps <= p_I^high
    ! hmjb>>            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    !--------------------------------------------------------------------------------- 
    ! Input/Output Variables
    !--------------------------------------------------------------------------------- 

    ! Model info
    INTEGER         , INTENT(IN) :: ncols  
    INTEGER         , INTENT(IN) :: kmax   
    !REAL(KIND=r8)   , INTENT(IN) :: sigbot(kmax+1)  ! Sigma cordinate at bottom of layer
    !REAL(KIND=r8)   , INTENT(IN) :: delsig(kmax)  ! Layer thickness (sigma)
    INTEGER(KIND=i8), INTENT(IN) :: imask (ncols) ! Ocean/Land mask

    ! Atmospheric Fields
    REAL(KIND=r8), INTENT(in   ) :: prsi  (ncols,kMax+1)
    REAL(KIND=r8), INTENT(in   ) :: prsl  (ncols,kMax)
    REAL(KIND=r8), INTENT(in   ) :: phii  (nCols,kMax+1)
    REAL(KIND=r8), INTENT(in   ) :: phil  (nCols,kMax)
    !REAL(KIND=r8), INTENT(IN) :: Ps  (ncols)      ! Surface pressure (mb)
    REAL(KIND=r8), INTENT(IN) :: Te  (ncols,kmax) ! Temperature (K)
    REAL(KIND=r8), INTENT(IN) :: Qe  (ncols,kmax) ! Specific Humidity (g/g)
    REAL(KIND=r8), INTENT(IN) :: tsea(ncols)
    REAL(KIND=r8), INTENT(IN) :: FlipPbot (ncols,kmax)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8), INTENT(IN) :: QCF(ncols,kMax)
    REAL(KIND=r8), INTENT(IN) :: QCL(ncols,kMax)

    ! Cloud properties
    REAL(KIND=r8), INTENT(OUT) :: clwp (ncols,kmax) ! Cloud Liquid Water Path
    REAL(KIND=r8), INTENT(OUT) :: lmixr(ncols,kmax) ! Ice/Water mixing ratio
    REAL(KIND=r8), INTENT(OUT) :: fice (ncols,kmax) ! Fractional amount of cloud that is ice
    REAL(KIND=r8), INTENT(OUT) :: rei  (ncols,kmax) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: rel  (ncols,kmax) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: taud (ncols,kmax) ! Shortwave cloud optical depth
    REAL(KIND=r8) :: emis(ncols,kMax)       ! cloud emissivity (fraction)

    !--------------------------------------------------------------------------------- 
    ! Parameters
    !--------------------------------------------------------------------------------- 

    REAL(KIND=r8), PARAMETER :: abarl=2.261e-2_r8
    REAL(KIND=r8), PARAMETER :: bbarl=1.4365_r8
    REAL(KIND=r8), PARAMETER :: abari=3.448e-3_r8
    REAL(KIND=r8), PARAMETER :: bbari=2.431_r8
    REAL(KIND=r8), PARAMETER :: kabsl=0.090361_r8               ! longwave liquid absorption coeff (m**2/g)

    REAL(KIND=r8), PARAMETER :: clwc0   = 0.21_r8 ! Reference liquid water concentration (g/m3)        
    REAL(KIND=r8), PARAMETER :: reimin  = 10.0_r8 ! Minimum of Ice particle efective radius (microns)  
    REAL(KIND=r8), PARAMETER :: reirnge = 20.0_r8 ! Range of Ice particle efective radius (microns)   
    REAL(KIND=r8), PARAMETER :: sigrnge = 0.4_r8  ! Normalized pressure range                         
    REAL(KIND=r8), PARAMETER :: sigmax  = 0.4_r8  ! Normalized pressure maximum                       

    !REAL(KIND=r8), PARAMETER :: pptop = 0.005_r8       ! Model-top presure                                 

    !--------------------------------------------------------------------------------- 
    ! Local Variables
    !--------------------------------------------------------------------------------- 
    REAL(KIND=r8) :: Psurf  (ncols) ! Height at middle of layer (m)
    REAL(KIND=r8) :: DeltaP  (ncols,kmax) ! Height at middle of layer (m)

    REAL(KIND=r8) :: landm       (1:ncols)   ! Land fraction ramped
    REAL(KIND=r8) :: icefrac     (1:ncols)   ! Ice fraction
    REAL(KIND=r8) :: snowh       (1:ncols)   ! Snow depth over land, water equivalent (m)
    REAL(KIND=r8) :: ocnfrac     (1:nCols)   ! Ocean fraction
    REAL(KIND=r8) :: landfrac    (1:ncols)   ! Land fraction

    REAL(KIND=r8) :: hl     (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: rhl    (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: pw     (ncols)        ! precipitable water (kg/m2)
    REAL(KIND=r8) :: Zibot  (ncols,kmax+1) ! Height at middle of layer (m)
    REAL(KIND=r8) :: emziohl(ncols,kmax+1) ! exponential of Minus zi Over hl (no dim)

!    REAL(KIND=r8) :: tauxcl(ncols,kmax)    ! extinction optical depth of liquid phase
!    REAL(KIND=r8) :: tauxci(ncols,kmax)    ! extinction optical depth of ice phase

    !-- Aux variables

    INTEGER :: i,k
!    REAL(KIND=r8) :: weight
    REAL(KIND=r8) :: kabs                   ! longwave absorption coeff (m**2/g)
    REAL(KIND=r8) :: kabsi                  ! ice absorption coefficient
    

    !--------------------------------------------------------------------------------- 
    !--------------------------------------------------------------------------------- 

    clwp=0.0_r8
    lmixr=0.0_r8
    fice=0.0_r8
    rei=0.0_r8
    rel=0.0_r8
    pw=0.0_r8

    landm=0.0_r8
    snowh=0.0_r8
    icefrac=0.0_r8
    landfrac=0.0_r8
    ocnfrac=0.0_r8

    DO k=1,kmax
       DO i = 1,ncols
         DeltaP(i,k) = ((prsi(i,k)) - (prsi(i,k+1)))/prsi(i,1)
       END DO
    END DO

    DO i=1,nCols
       Psurf(i) =prsi  (i,1)/100.0_r8
    END DO
 
    DO i=1,nCols
       IF(schemes == 1 .and. imask(i) == 13_i8) snowh(i)=5.0_r8
       IF(schemes == 2 .and. imask(i) == 13_i8) snowh(i)=5.0_r8
       IF(schemes == 3 .and. imask(i) == 15_i8) snowh(i)=5.0_r8    
       IF(imask(i) >   0_i8)THEN
          ! land
          icefrac  (i)=0.0_r8
          landfrac (i)=1.0_r8
          ocnfrac  (i)=0.0_r8
       ELSE
          ! water/ocean
          landfrac  (i) =0.0_r8
          ocnfrac   (i) =1.0_r8
          IF(ocnfrac(i).GT.0.01_r8.AND.ABS(tsea(i)).LT.260.0_r8) THEN
             icefrac(i) = 1.0_r8
             ocnfrac(i) = 1.0_r8
          ENDIF
       END IF
    END DO

    ! Heights corresponding to sigma at middle of layer: sig(k)
    ! Assuming isothermal atmosphere within each layer
    DO i=1,nCols
       Zibot(i,1) = 0.0_r8
       DO k=2,kMax
          Zibot(i,k) = Zibot(i,k-1) + (gasr/grav)*Te(i,k-1)* &
    !               LOG(sigbot(k-1)/sigbot(k))
               LOG(FlipPbot(i,kmax+2-k)/FlipPbot(i,kmax+1-k))
       END DO
    END DO
    
    DO i=1,ncols
       Zibot(i,kmax+1)=Zibot(i,kmax)+(gasr/grav)*Te(i,kmax)* &
            LOG(FlipPbot(i,1)/pptop)
    END DO

    ! precitable water, pw = sum_k { delsig(k) . Qe(k) } . Ps . 100 / g
    !                   pw = sum_k { Dp(k) . Qe(k) } / g
    !
    ! 100 is to change from mbar to pascal
    ! Dp(k) is the difference of pressure (N/m2) between bottom and top of layer
    ! Qe(k) is specific humidity in (g/g)
    ! gravity is m/s2 => so pw is in Kg/m2
    DO k=1,kmax
       DO i = 1,ncols
          pw(i) = pw(i) + DeltaP(i,k)*Qe(i,k)
       END DO
    END DO
    DO i = 1,ncols
       pw(i)=100.0_r8*pw(i)*Psurf(i)/grav
    END DO
    !
    ! diagnose liquid water scale height from precipitable water
    DO i=1,ncols
       hl(i)  = 700.0_r8*LOG(MAX(pw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,kmax+1
       DO i=1,ncols
    !          emziohl(i,k) = EXP(-zibot(i,k)/hl(i))
          emziohl(i,k) = EXP(-zibot(i,k)*rhl(i))
       END DO
    END DO
    !    DO i=1,ncols
    !       emziohl(i,kmax+1) = 0.0_r8
    !    END DO

    ! The units are g/m2.
    DO k=1,kmax
       DO i=1,ncols
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k) - emziohl(i,k+1))
       END DO
    END DO

    ! If we want to calculate the 'droplets/cristals' mixing ratio, we need
    ! to find the amount of dry air in each layer. 
    !
    !             dry_air_path = int rho_air dz  
    !
    ! This can be simply done using the hydrostatic equation:
    !
    !              dp/dz   = -rho grav
    !              dp/grav = -rho dz
    !
    !
    ! The units are g/m2. The factor 1e5 accounts for the change
    !  mbar to Pa and kg/m2 to g/m2. 
    DO k=1,kmax
       DO i=1,ncols
          lmixr(i,k)=clwp(i,k)*grav*1.0e-5_r8/DeltaP(i,k)/Psurf(i)
       END DO
    END DO

    !
    ! Cloud water and ice particle sizes, saved in physics buffer for radiation
    ! Author: Byron Boville  Sept 06, 2002, assembled from existing subroutines
    !
    CALL cldefr(&
         ncols                             , & !INTEGER , INTENT(in ) :: pcols                ! number of atmospheric columns
         kmax                              , & !INTEGER , INTENT(in ) :: kMax                 ! number of vertical levels
         Te         (1:ncols,1:kMax)       , & !REAL(r8), INTENT(in ) :: t       (pcols,kMax) ! Temperature
         rel        (1:ncols,1:kMax)       , & !REAL(r8), INTENT(out) :: rel     (pcols,kMax) ! Liquid effective drop size (microns)
         rei        (1:ncols,1:kMax)       , & !REAL(r8), INTENT(out) :: rei     (pcols,kMax) ! Ice effective drop size (microns)
         landm      (1:ncols)              , & !REAL(r8), INTENT(in ) :: landm   (pcols)      !
         icefrac    (1:ncols)              , & !REAL(r8), INTENT(in ) :: icefrac (pcols)      ! Ice fraction
         snowh      (1:ncols)                ) !REAL(r8), INTENT(in ) :: snowh   (pcols)      ! Snow depth over land, water equivalent (m)


    ! define fractional amount of cloud that is ice
    ! if warmer than -10 degrees c then water phase
    ! docs CCM3, eq 4.a.16.1     
    ! allcld_liq = state%q(:,:,ixcldliq)
    ! allcld_ice = state%q(:,:,ixcldice)
    IF  (TRIM(ILCON) == 'YES' .OR. TRIM(ILCON) == 'LSC' )THEN
       DO k=1,kmax
          DO i=1,ncols
             fice(i,k)=MAX(MIN((263.16_r8-Te(i,k))*0.05_r8,1.0_r8),0.0_r8)
             !fice(i,k) = allcld_ice(i,k) /max(1.e-10_r8,(allcld_ice(i,k) + allcld_liq(i,k)))  
          END DO
       END DO    
    ELSE IF ( TRIM(ILCON) == 'MIC'.or. TRIM(ILCON).EQ.'HWRF' .or. TRIM(ILCON).EQ.'HGFS'.or.&
              TRIM(ILCON).EQ.'UKMO' .or. TRIM(ILCON).EQ.'MORR' .or.TRIM(ILCON).EQ.'HUMO') THEN
       DO k=1,kmax
          DO i=1,ncols
             fice(i,k) = MIN(MAX( QCF(i,k) /max(1.e-10_r8,(QCF(i,k) + QCL(i,k))),0.000_r8),1.0_r8)
             !fice(i,k) = allcld_ice(i,k) /max(1.e-10_r8,(allcld_ice(i,k) + allcld_liq(i,k)))
          END DO
       END DO
    END IF

    ! Compute optical depth from liquid water
    DO k=1,kMax
       DO i=1,nCols
          !note that optical properties for ice valid only
          !in range of 13 > rei > 130 micron (Ebert and Curry 92)
          !if ( microp_scheme .eq. 'MG' ) then
          kabsi = 0.005_r8 + 1.0_r8/min(max(13._r8,rei(i,k)),130._r8)
          !else if ( microp_scheme .eq. 'RK' ) then
          !   kabsi = 0.005_r8 + 1._r8/rei(i,k)
          !END IF
          !     (m**2/g)
          kabs = kabsl*(1.0_r8-fice(i,k)) + kabsi*fice(i,k) 
          ! cloud emissivity (fraction)
          emis(i,k) = 1.0_r8 - exp(-1.66_r8*kabs*clwp(i,k))
          ! cloud optical depth
          taud(i,k) = kabs*clwp(i,k)! g/m2
       END DO
    END DO


  END SUBROUTINE Cloud_Micro_CAM5

  !===============================================================================
  SUBROUTINE cldefr( &
       ncols    , &!INTEGER, INTENT(in) :: ncol                  ! number of atmospheric columns
       kMax    , &!INTEGER, INTENT(in) :: kMax                  ! number of vertical levels
       t       , &!REAL(r8), INTENT(in) :: t       (pcols,kMax)        ! Temperature
       rel     , &!REAL(r8), INTENT(out) :: rel(pcols,kMax)      ! Liquid effective drop size (microns)
       rei     , &!REAL(r8), INTENT(out) :: rei(pcols,kMax)      ! Ice effective drop size (microns)
       landm   , &!REAL(r8), INTENT(in) :: landm   (pcols)
       icefrac , &!REAL(r8), INTENT(in) :: icefrac (pcols)       ! Ice fraction
       snowh     )!REAL(r8), INTENT(in) :: snowh   (pcols)         ! Snow depth over land, water equivalent (m)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud water and ice particle size 
    ! 
    ! Method: 
    ! use empirical formulas to construct effective radii
    ! 
    ! Author: J.T. Kiehl, B. A. Boville, P. Rasch
    ! 
    !-----------------------------------------------------------------------

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: ncols                  ! number of atmospheric columns
    INTEGER, INTENT(in) :: kMax                  ! number of vertical levels

    REAL(r8), INTENT(in) :: icefrac (ncols)       ! Ice fraction
    REAL(r8), INTENT(in) :: t       (ncols,kMax)  ! Temperature
    REAL(r8), INTENT(in) :: landm   (ncols)
    REAL(r8), INTENT(in) :: snowh   (ncols)       ! Snow depth over land, water equivalent (m)
    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: rel(ncols,kMax)      ! Liquid effective drop size (microns)
    REAL(r8), INTENT(out) :: rei(ncols,kMax)      ! Ice effective drop size (microns)
    !

    !++pjr
    ! following Kiehl
    CALL reltab(ncols,kMax,  t,  landm, icefrac, rel, snowh)

    ! following Kristjansson and Mitchell
    CALL reitab(ncols,kMax, t, rei)
    !--pjr
    !
    !
    RETURN
  END SUBROUTINE cldefr

  !===============================================================================
  SUBROUTINE reltab(nCols,kMax, t, landm, icefrac, rel, snowh)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud water size
    ! 
    ! Method: 
    ! analytic formula following the formulation originally developed by J. T. Kiehl
    ! 
    ! Author: Phil Rasch
    ! 
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER  , INTENT(in) :: nCols
    INTEGER, INTENT(in) :: kMax                  ! number of vertical levels

    REAL(r8), INTENT(in) :: icefrac(nCols)       ! Ice fraction
    REAL(r8), INTENT(in) :: snowh(nCols)         ! Snow depth over land, water equivalent (m)
    REAL(r8), INTENT(in) :: landm(nCols)         ! Land fraction ramping to zero over ocean
    REAL(r8), INTENT(in) :: t(nCols,kMax)        ! Temperature

    !
    ! Output arguments
    !
    REAL(r8), INTENT(out) :: rel(nCols,kMax)      ! Liquid effective drop size (microns)
    !
    !---------------------------Local workspace-----------------------------
    !
    INTEGER i,k               ! Lon, lev indices
    REAL(r8) rliqland         ! liquid drop size if over land
    REAL(r8) rliqocean        ! liquid drop size if over ocean
    REAL(r8) rliqice          ! liquid drop size if over sea ice
    !
    !-----------------------------------------------------------------------
    !
    rliqice   = 14.0_r8
    !rliqocean = 14.0_r8
    rliqocean = 10.0_r8
    rliqland  = 8.0_r8
    DO k=1,kMax
       DO i=1,nCols
          ! jrm Reworked effective radius algorithm
          ! Start with temperature-dependent value appropriate for continental air
          ! Note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * MIN(1.0_r8,MAX(0.0_r8,(tmelt-t(i,k))*0.05_r8))
          ! Modify for snow depth over land
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * MIN(1.0_r8,MAX(0.0_r8,snowh(i)*10.0_r8))
          ! Ramp between polluted value over land to clean value over ocean.
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * MIN(1.0_r8,MAX(0.0_r8,1.0_r8-landm(i)))
          ! Ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k) = rel(i,k) + (rliqice-rel(i,k)) * MIN(1.0_r8,MAX(0.0_r8,icefrac(i)))
          ! end jrm
       END DO
    END DO
  END SUBROUTINE reltab



  !===============================================================================
  SUBROUTINE reitab(nCols,kMax, t, re)
    !

    INTEGER  , INTENT(in) :: nCols
    INTEGER  , INTENT(in) :: kMax
    REAL(r8), INTENT(in ) :: t(nCols,kMax)
    REAL(r8), INTENT(out) :: re(nCols,kMax)
    REAL(r8) :: corr
    INTEGER :: i
    INTEGER :: k
    INTEGER :: index
    !
    !       Tabulated values of re(T) in the temperature interval
    !       180 K -- 274 K; hexagonal columns assumed:
    !
    REAL(KIND=r8), PARAMETER :: retab(95)=(/                                                 &
         5.92779_r8, 6.26422_r8, 6.61973_r8, 6.99539_r8, 7.39234_r8,        &
         7.81177_r8, 8.25496_r8, 8.72323_r8, 9.21800_r8, 9.74075_r8, 10.2930_r8,        &
         10.8765_r8, 11.4929_r8, 12.1440_r8, 12.8317_r8, 13.5581_r8, 14.2319_r8,         &
         15.0351_r8, 15.8799_r8, 16.7674_r8, 17.6986_r8, 18.6744_r8, 19.6955_r8,        &
         20.7623_r8, 21.8757_r8, 23.0364_r8, 24.2452_r8, 25.5034_r8, 26.8125_r8,        &
         27.7895_r8, 28.6450_r8, 29.4167_r8, 30.1088_r8, 30.7306_r8, 31.2943_r8,         &
         31.8151_r8, 32.3077_r8, 32.7870_r8, 33.2657_r8, 33.7540_r8, 34.2601_r8,         &
         34.7892_r8, 35.3442_r8, 35.9255_r8, 36.5316_r8, 37.1602_r8, 37.8078_r8,        &
         38.4720_r8, 39.1508_r8, 39.8442_r8, 40.5552_r8, 41.2912_r8, 42.0635_r8,        &
         42.8876_r8, 43.7863_r8, 44.7853_r8, 45.9170_r8, 47.2165_r8, 48.7221_r8,        &
         50.4710_r8, 52.4980_r8, 54.8315_r8, 57.4898_r8, 60.4785_r8, 63.7898_r8,        &
         65.5604_r8, 71.2885_r8, 75.4113_r8, 79.7368_r8, 84.2351_r8, 88.8833_r8,        &
         93.6658_r8, 98.5739_r8, 103.603_r8, 108.752_r8, 114.025_r8, 119.424_r8,         &
         124.954_r8, 130.630_r8, 136.457_r8, 142.446_r8, 148.608_r8, 154.956_r8,        &
         161.503_r8, 168.262_r8, 175.248_r8, 182.473_r8, 189.952_r8, 197.699_r8,        &
         205.728_r8, 214.055_r8, 222.694_r8, 231.661_r8, 240.971_r8, 250.639_r8/)        
    !
    !
    DO k=1,kMax
       DO i=1,nCols
          index = INT(t(i,k)-179.0_r8)
          index = MIN(MAX(index,1),94)
          corr = t(i,k) - INT(t(i,k))
          re(i,k) = retab(index)*(1.0_r8-corr)                &
               +retab(index+1)*corr
          !           re(i,k) = amax1(amin1(re(i,k),30.0_r8),10.0_r8)
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE reitab

  SUBROUTINE Cloud_Micro_CCM3(&
       ! Model info
       ncols, kmax ,  imask   , &
       ! Atmospheric Fields
       prsi   ,prsl   ,phii   ,phil    ,&
       Te   , Qe    , tsea  , FlipPbot ,        &
       ! Cloud properties
       clwp , lmixr, fice  , rei   , rel   , taud       )
    IMPLICIT NONE

    ! As in the CCM2, cloud optical properties in the CCM3 are accounted for using
    ! the Slingo (1989) parameterization for liquid water droplet clouds. This
    ! scheme relates the extinction optical depth, the single-scattering albedo,
    ! and the asymmetry parameter to the cloud liquid water path and cloud drop
    ! effective radius. The latter two microphysical cloud properties were
    ! statically specified in the CCM2. In particular, in-cloud liquid water paths
    ! were evaluated from a prescribed, meridionally and height varying, but
    ! time independent, cloud liquid water density profile, rho_l(z), which
    ! was analytically determined on the basis of a meridionally specified
    ! liquid water scale height (e.g. see Kiehl et al., 1994; Kiehl, 1991).
    ! The cloud drop effective radius was simplly specified to be 10microns
    ! for all clouds. The CCM3 continues to diagnose cloud optical properties,
    ! but relaxes the rigid CCM2 framework. CCM3 employs the same exponentially
    ! decaying vertical profile for in-cloud water concentration
    !
    !             rho_l=rho_l^0*exp(-z/h_l)               eq 4.a.11
    !
    ! , where rho_l^0=0.21g/m3. Instead of specifying a zonally symmetric meridional
    ! dependence for the cloud water scale heigh, h_l, it is locally diagnosed
    ! as a function of the vertically integrated water vapor (precipitable water) 
    !
    !          h_l=700 ln [1+\frac{1}{g} \int_pT^ps q dp]  eq 4.a.12
    !
    ! hmjb> It is not explained, but the units of h_l must be meters, the same 
    ! hmjb> of the height, z.
    !
    ! The cloud water path (CWP) is determined by integrating the liquid
    ! water concentration using
    !
    !                 cwp = int rho_l dz     eq. 4.a.13
    ! 
    ! Which can be analytically evaluated for an arbitrary layer k as
    !
    !  rho_l^0 h_l [exp(-z_bot(k)/h_l) - exp(-z_top(k)/h_l)]   eq. 4.a.14
    !
    ! Where z_bot and z_top are the heights of the k'th layer interfaces.
    !
    ! hmjb> It is not explained, but the units of clwp must be g/m2
    ! hmjb> since it is the integral of rho_l*dz (eq.4.a.13)
    !
    ! CCM3 Documentation, pg 50
    ! Observational studies have shown a distinct difference between
    ! maritime and continental effective cloud drop size, r_e, for warm
    ! clouds. For this reason, the CCM3 differentiates between the cloud
    ! drop effective radius for clouds diagnosed over maritime and
    ! continental regimes (Kiehl, 1994). Over the ocean, the cloud drop
    ! effective radius for liquid water clouds, r_el, is specified to be
    ! 10microns, as in the CCM3. Over land masses r_el is determinedusing
    !
    ! r_el = 5 microns             T > -10oC
    !      = 5-5(t+10)/20 microns  -30oC <= T <= -10oC     eq. 4.a.14.1
    !      = r_ei                  T < -30oC
    !
    ! An ice particle effective radius, r_ei, is also diagnosed by CCM3,
    ! which at the moment amounts to a specification of ice radius as a
    ! function of normalized pressure
    !
    ! r_ei = 10 microns                                 p/ps > p_I^high   
    !      = r_ei^max - (r_ei^max - r_ei^min)           p/ps <= p_I^high       eq. 4.a.15.1
    !            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    ! where r_ei^max=30microns, r_ei^min=10microns, p_I^high=0.4 and p_I^low=0.0
    !
    ! hmjb>> I think there is a typo in the equation, otherwise the 
    ! hmjb>> expression for r_ei is not a continuous funcion of p/ps.
    ! hmjb>> For p/ps=p_I^high, r_ei should be r_ei^min and not r_ei^max.
    ! hmjb>> The correct equation is:
    ! hmjb>> r_ei = 10 microns                                 p/ps > p_I^high   
    ! hmjb>>      = R_EI^MIN - (r_ei^max - r_ei^min)           p/ps <= p_I^high
    ! hmjb>>            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    !--------------------------------------------------------------------------------- 
    ! Input/Output Variables
    !--------------------------------------------------------------------------------- 

    ! Model info
    INTEGER         , INTENT(IN) :: ncols  
    INTEGER         , INTENT(IN) :: kmax   
    !REAL(KIND=r8)   , INTENT(IN) :: sigmid(kmax)  ! Sigma cordinate at middle of layer
    !REAL(KIND=r8)   , INTENT(IN) :: sigbot(kmax+1)  ! Sigma cordinate at bottom of layer
    !REAL(KIND=r8)   , INTENT(IN) :: delsig(kmax)  ! Layer thickness (sigma)
    INTEGER(KIND=i8), INTENT(IN) :: imask (ncols) ! Ocean/Land mask

    ! Atmospheric Fields
    REAL(KIND=r8), INTENT(in   ) :: prsi  (ncols,kMax+1)  !  pressure (PA)
    REAL(KIND=r8), INTENT(in   ) :: prsl  (ncols,kMax)  !  pressure (PA)
    REAL(KIND=r8), INTENT(in   ) :: phii  (nCols,kMax+1)  !  Height (m)
    REAL(KIND=r8), INTENT(in   ) :: phil  (nCols,kMax)  !  Height (m)
    !REAL(KIND=r8), INTENT(IN) :: Ps  (ncols)      ! Surface pressure (mb)
    REAL(KIND=r8), INTENT(IN) :: Te  (ncols,kmax) ! Temperature (K)
    REAL(KIND=r8), INTENT(IN) :: Qe  (ncols,kmax) ! Specific Humidity (g/g)
    REAL(KIND=r8), INTENT(IN) :: tsea(ncols)
    REAL(KIND=r8), INTENT(IN) :: FlipPbot (ncols,kmax)  ! Pressure at bottom of layer (mb)

    ! Cloud properties
    REAL(KIND=r8), INTENT(OUT) :: clwp (ncols,kmax) ! Cloud Liquid Water Path
    REAL(KIND=r8), INTENT(OUT) :: lmixr(ncols,kmax) ! Ice/Water mixing ratio
    REAL(KIND=r8), INTENT(OUT) :: fice (ncols,kmax) ! Fractional amount of cloud that is ice
    REAL(KIND=r8), INTENT(OUT) :: rei  (ncols,kmax) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: rel  (ncols,kmax) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: taud (ncols,kmax) ! Shortwave cloud optical depth

    !--------------------------------------------------------------------------------- 
    ! Parameters
    !--------------------------------------------------------------------------------- 

    REAL(KIND=r8), PARAMETER :: abarl=2.261e-2_r8
    REAL(KIND=r8), PARAMETER :: bbarl=1.4365_r8
    REAL(KIND=r8), PARAMETER :: abari=3.448e-3_r8
    REAL(KIND=r8), PARAMETER :: bbari=2.431_r8

    REAL(KIND=r8), PARAMETER :: clwc0   = 0.21_r8 ! Reference liquid water concentration (g/m3)        
    REAL(KIND=r8), PARAMETER :: reimin  = 10.0_r8 ! Minimum of Ice particle efective radius (microns)  
    REAL(KIND=r8), PARAMETER :: reirnge = 20.0_r8 ! Range of Ice particle efective radius (microns)   
    REAL(KIND=r8), PARAMETER :: sigrnge = 0.4_r8  ! Normalized pressure range                         
    REAL(KIND=r8), PARAMETER :: sigmax  = 0.4_r8  ! Normalized pressure maximum                       

    !REAL(KIND=r8), PARAMETER :: pptop = 0.005_r8  ! Model-top presure                                 

    !--------------------------------------------------------------------------------- 
    ! Local Variables
    !--------------------------------------------------------------------------------- 
    REAL(KIND=r8) :: Psurf     (ncols)        !  Surface pressure (mb)
    REAL(KIND=r8) :: DeltaP    (ncols,kmax)        !  Surface pressure (mb)

    REAL(KIND=r8) :: hl     (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: rhl    (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: pw     (ncols)        ! precipitable water (kg/m2)
    REAL(KIND=r8) :: Zibot  (ncols,kmax+1) ! Height at middle of layer (m)
    REAL(KIND=r8) :: emziohl(ncols,kmax+1) ! exponential of Minus zi Over hl (no dim)

    REAL(KIND=r8) :: tauxcl(ncols,kmax)    ! extinction optical depth of liquid phase
    REAL(KIND=r8) :: tauxci(ncols,kmax)    ! extinction optical depth of ice phase

    !-- Aux variables

    INTEGER :: i,k
    REAL(KIND=r8) :: weight

    !--------------------------------------------------------------------------------- 
    !--------------------------------------------------------------------------------- 

    clwp=0.0_r8
    lmixr=0.0_r8
    fice=0.0_r8
    rei=0.0_r8
    rel=0.0_r8
    pw=0.0_r8
    DO k=1,kmax
       DO i = 1,ncols
         DeltaP(i,k) = ((prsi(i,k)) - (prsi(i,k+1)))/prsi(i,1)
       END DO
    END DO

    DO i=1,nCols
       Psurf(i) =prsi  (i,1)/100.0_r8
    END DO
    ! Heights corresponding to sigma at middle of layer: sig(k)
    ! Assuming isothermal atmosphere within each layer
    DO i=1,nCols
       Zibot(i,1) = 0.0_r8
       DO k=2,kMax
          Zibot(i,k) = Zibot(i,k-1) + (gasr/grav)*Te(i,k-1)* &
!               LOG(sigbot(k-1)/sigbot(k))
               LOG(FlipPbot(i,kmax+2-k)/FlipPbot(i,kmax+1-k))
       END DO
    END DO
    DO i=1,ncols
       Zibot(i,kmax+1)=Zibot(i,kmax)+(gasr/grav)*Te(i,kmax)* &
            LOG(FlipPbot(i,1)/pptop)
    END DO

! precitable water, pw = sum_k { delsig(k) . Qe(k) } . Ps . 100 / g
!                   pw = sum_k { Dp(k) . Qe(k) } / g
!
! 100 is to change from mbar to pascal
! Dp(k) is the difference of pressure (N/m2) between bottom and top of layer
! Qe(k) is specific humidity in (g/g)
! gravity is m/s2 => so pw is in Kg/m2
    DO k=1,kmax
       DO i = 1,ncols
          !pw(i) = pw(i) + delsig(k)*Qe(i,k)
          pw(i) = pw(i) + DeltaP(i,k)*Qe(i,k)
       END DO
    END DO
    DO i = 1,ncols
       pw(i)=100.0_r8*pw(i)*Psurf(i)/grav
    END DO
    !
    ! diagnose liquid water scale height from precipitable water
    DO i=1,ncols
       hl(i)  = 700.0_r8*LOG(MAX(pw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,kmax+1
       DO i=1,ncols
!          emziohl(i,k) = EXP(-zibot(i,k)/hl(i))
          emziohl(i,k) = EXP(-zibot(i,k)*rhl(i))
       END DO
    END DO
!    DO i=1,ncols
!       emziohl(i,kmax+1) = 0.0_r8
!    END DO

    ! The units are g/m2.
    DO k=1,kmax
       DO i=1,ncols
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k) - emziohl(i,k+1))
       END DO
    END DO

! If we want to calculate the 'droplets/cristals' mixing ratio, we need
! to find the amount of dry air in each layer. 
!
!             dry_air_path = int rho_air dz  
!
! This can be simply done using the hydrostatic equation:
!
!              dp/dz = -rho grav
!            dp/grav = -rho dz
!
!
! The units are g/m2. The factor 1e5 accounts for the change
!  mbar to Pa and kg/m2 to g/m2. 
    DO k=1,kmax
       DO i=1,ncols
          !lmixr(i,k)=clwp(i,k)*grav*1.0e-5_r8/delsig(k)/Psurf(i)
          lmixr(i,k)=clwp(i,k)*grav*1.0e-5_r8/DeltaP(i,k)/Psurf(i)

       END DO
    END DO

    ! determine Ice particle Effective Radius (rei)
    ! as function of normalized pressure 
    ! docs CCM3, eq 4.a.15.1
    DO k=1,kmax
      ! weight   = MIN((sigmid(k)-sigmax)/sigrnge,0.0_r8)
       DO i=1,ncols
          weight   = MIN(((prsl(i,k)/prsl(i,1))-sigmax)/sigrnge,0.0_r8)
          rei(i,k) = reimin - reirnge*weight
       END DO
    END DO

    ! define fractional amount of cloud that is ice
    ! if warmer than -10 degrees c then water phase
    ! docs CCM3, eq 4.a.16.1     
    DO k=1,kmax
       DO i=1,ncols
          fice(i,k)=MAX(MIN((263.16_r8-Te(i,k))*0.05_r8,1.0_r8),0.0_r8)
       END DO
    END DO    

    ! determine Liquid particle Effective Radius (rel) 
    ! as function of normalized pressure
    ! docs CCM3, eq 4.a.15.1
    DO k=1,kmax
       DO i=1,ncols
          IF (imask(i) .lt. 1_i8) THEN
             !ocean  
             rel(i,k) = 10.0_r8
          ELSE
             !land
             rel(i,k) = 5.0_r8+5.0_r8*fice(i,k)
             IF(fice(i,k) == 1.0_r8) rel(i,k) = rei(i,k)
          END IF
       END DO
    END DO    

    ! Compute optical depth from liquid water
    DO k=1,kmax
       DO i=1,ncols
          ! ccm3 manual, page 53, eqs 4.b.3 and 4.b.7
          tauxcl(i,k) = clwp(i,k)*(abarl + bbarl/rel(i,k))*(1.0_r8-fice(i,k))
          tauxci(i,k) = clwp(i,k)*(abari + bbari/rei(i,k))*fice(i,k)
          IF (tsea(i) > 0.0_r8) THEN
             taud(i,k)=0.70_r8*(tauxcl(i,k)+tauxci(i,k))
          ELSE
             taud(i,k)=1.00_r8*(tauxcl(i,k)+tauxci(i,k))
          END IF
       END DO
    END DO


  END SUBROUTINE Cloud_Micro_CCM3

  SUBROUTINE ASTROPAR(ncols,colrad,lonrad,id,tod,yrl,date,solar,cosine,ratio)
    IMPLICIT NONE

    ! Input Model Info
    INTEGER      , INTENT(IN ) :: ncols
    REAL(KIND=r8), INTENT(IN ) :: colrad(ncols)
    REAL(KIND=r8), INTENT(IN ) :: lonrad(ncols)

    ! Input time Info 
    INTEGER      , INTENT(IN ) :: id(4)
    REAL(KIND=r8), INTENT(IN ) :: tod
    REAL(KIND=r8), INTENT(IN ) :: yrl

    ! Ouput Variables
    REAL(KIND=r8), INTENT(OUT) :: date          ! julian day
    REAL(KIND=r8), INTENT(OUT) :: solar         ! solar constant
    REAL(KIND=r8), INTENT(OUT) :: cosine(ncols) ! cosine of solar zenith angle
    REAL(KIND=r8), INTENT(OUT) :: ratio         ! factor relating to the distance between the earth and the sun

    ! Local Variables
    REAL(KIND=r8)    :: delta
    REAL(KIND=r8)    :: etime
    REAL(KIND=r8)    :: frh
    REAL(KIND=r8)    :: atime
    REAL(KIND=r8)    :: sindel
    REAL(KIND=r8)    :: cosdel
    REAL(KIND=r8)    :: coslat (ncols)
    REAL(KIND=r8)    :: sinlat (ncols)
    INTEGER :: i

    !-----------------------------------------------------------------------
    !.. delta ;solar inclination
    !.. etime ;correction factor to local time
    !.. ratio ;factor relating to the distance between the earth and the sun
    !.. date  ;julian day
    !
    CALL radtim (id, delta, ratio, etime, tod, date, yrl)
    solar=solcon*ratio

    sindel=SIN(delta)
    cosdel=COS(delta)
    DO i=1,ncols
       coslat(i)  = cosdel*SIN(colrad(i))
       sinlat(i)  = sindel*COS(colrad(i))
    END DO
    frh =( MOD(tod +0.03125_r8,3.6e3_r8)-0.03125_r8)/3.6e3_r8
    DO i=1,ncols
       atime =etime +pai12*(12.0_r8-id(1)-lonrad(i)*fim24-frh )
       cosine(i)=sinlat(i)  + coslat(i)  * COS(atime)
    END DO

  END SUBROUTINE ASTROPAR
  
  
  SUBROUTINE RRTM_CRDTF(ncols,aod,asy,ssa,aod_clrd,asy_clrd,ssa_clrd)
! 
!
! written by Tarasova on November 2015
! re-calculates aerosol optical parameters at the Clirad SW 2007 wavelengths from
! that given at the RRTM SW wavelengths inside the Global Model in RadiationDriver, Subroutine spmrad:
!
!    3.462, 2.789, 2.325, 2.046, 1.784, 1.462, 1.271, 1.010, 0.702, 0.533, 0.393, 0.304, 0.232, 8.021
!
!    Inverted: 0.232, 0.304, 0.393, 0.533, 0.702, 1.010, 1.271, 1.462, 1.784, 2.046, 2.325, 2.789, 3.462, 8.021
!
!  data prepared for RRTM in Global model (RadiationDriver, Subroutine spmrad):
!
!  aod(ncols,14)   "aerosol optical depth"
!  asy(ncols,14)   "asymmetry factor"
!  ssa(ncols,14)   "single scattering albedo" 
!  
!  
!  recalculated at the new wavelengths  CLIRAD_SW2007 (mcm):  (RadiationDriver, Subroutine spmrad)
!
!   0.252, 0.313, 0.512, 0.772, 0.960, 5.61, 1.745, 6.135
!
!  aod_clrd(ncols,8)
!  asy_clrd(ncols,8)
!  ssa_clrd(ncols,8)


 IMPLICIT NONE

!INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers 
!
  INTEGER  ::  i, j, jinv, ncols
!
  REAL(KIND=r8), INTENT(in) :: aod(ncols,14)
  REAL(KIND=r8), INTENT(in) :: asy(ncols,14)
  REAL(KIND=r8), INTENT(in) :: ssa(ncols,14)
!
  REAL(KIND=r8), INTENT(out) :: aod_clrd(ncols,8)
  REAL(KIND=r8), INTENT(out) :: asy_clrd(ncols,8)
  REAL(KIND=r8), INTENT(out) :: ssa_clrd(ncols,8)
!
!
  REAL(KIND=r8) :: rrtmsw(ncols,14)
!
  REAL(KIND=r8) :: clrdsw(ncols,8)
!

!__________________________________________________
    aod_clrd=0.0_r8
    clrdsw=0.0_r8
    rrtmsw=0.0_r8
    DO j=1,13
       DO i=1,ncols
          jinv=14-j
          rrtmsw(i,j)=aod(i,jinv)
       END DO
    END DO
    DO i=1,ncols
       rrtmsw(i,14)=aod(i,14)
    END DO

    CALL interp_sw2007(ncols,clrdsw, rrtmsw)

    DO j=1,8
       DO i=1,ncols
          aod_clrd(i,j)=clrdsw(i,j)
       END DO
    END DO
!__________________________________________________
    ssa_clrd=0.0_r8
    clrdsw=0.0_r8
    rrtmsw=0.0_r8
    DO j=1,13
       DO i=1,ncols
          jinv=14-j
          rrtmsw(i,j)=ssa(i,jinv)
       END DO
    END DO
    DO i=1,ncols
       rrtmsw(i,14)=ssa(i,14)
    END DO

    CALL interp_sw2007(ncols,clrdsw,rrtmsw)

    DO j=1,8
       DO i=1,ncols
          ssa_clrd(i,j) =clrdsw(i,j)
       END DO
    END DO
!__________________________________________________
    asy_clrd=0.0_r8
    clrdsw=0.0_r8
    rrtmsw=0.0_r8
    DO j=1,13
       DO i=1,ncols
          jinv=14-j
          rrtmsw(i,j)=asy(i,jinv)
       END DO
    END DO
    DO i=1,ncols
       rrtmsw(i,14)=asy(i,14)
    END DO

    CALL interp_sw2007(ncols,clrdsw, rrtmsw)

    DO j=1,8
       DO i=1,ncols
          asy_clrd(i,j) =clrdsw(i,j)
       END DO
    END DO

!__________________________________________________
!
  END SUBROUTINE RRTM_CRDTF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE interp_sw2007(ncols,x,y)

! Tarasova: calculates aerosol optical parameters at the Clirad sw 2007 wavelengths:

    IMPLICIT NONE
!
!INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers   
    INTEGER      , INTENT(in) :: ncols
    REAL(KIND=r8), INTENT(in), DIMENSION(ncols,14) :: y

    REAL(KIND=r8), INTENT(out), DIMENSION(ncols,8) :: x

    REAL(KIND=r8), PARAMETER, DIMENSION(14) :: ly= (/0.232_r8, 0.304_r8, 0.393_r8, 0.533_r8, 0.702_r8, &
                      1.010_r8, 1.271_r8, 1.462_r8, 1.784_r8, 2.046_r8, 2.325_r8, 2.789_r8, 3.462_r8, 8.021_r8/)

    REAL(KIND=r8), PARAMETER, DIMENSION(8) :: lx= (/0.252_r8, 0.313_r8, 0.512_r8, 0.772_r8, 0.960_r8, & 
                      5.61_r8, 1.745_r8, 6.135_r8/)

    INTEGER  :: i,j,m

    INTEGER, PARAMETER, DIMENSION(8)  :: jk= (/1,2,3,5,5,13,8,13/)

    DO j=1,8
       DO i=1,ncols
          m=jk(j)
          x(i,j)=y(i,m)+(y(i,m)-y(i,m+1))/(ly(m)-ly(m+1))*(lx(j)-ly(m))
       END DO
    END DO

  END SUBROUTINE interp_sw2007


END MODULE ModRadiationDriver
